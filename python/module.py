#!/usr/bin/env python

from icecube import icetray, phys_services
#from icecube.ShowerLLH import GridSearch
from icecube.ShowerLLH import ShowerLLHFitParams
from icecube import dataclasses as dc
from SupportFunctions import b2e

import numpy as np
import time


class GridLLH(icetray.I3ConditionalModule):

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('DirectionFit', 'I3Particle to get direction from',
                          'ShowerPlane')
        self.AddParameter('RecoPulses', 'Pulse series for reconstruction',
                          'IceTopHLCVEMPulses')
        self.AddParameter('LLHTables', 'Input LLH tables', 0)
        self.AddParameter('binDict', 'Input bins', 0)
        self.AddParameter('grids', 'Input grids for LLH method', 0)
        self.AddParameter('runMode', 'Run in flexible(recommended) or static',
                          'flexible')
        self.AddParameter('saturated', 'Turn on saturation of DOMs', False)

    def Configure(self):
        self.DirectionFit = self.GetParameter('DirectionFit')
        self.recoPulses = self.GetParameter('RecoPulses')
        self.LLHTables = self.GetParameter('LLHTables')
        self.bins = self.GetParameter('binDict')
        self.grids = self.GetParameter('grids')
        self.runMode = self.GetParameter('runMode')
        self.saturated = self.GetParameter('saturated')

        if self.runMode == 'flexible':
            from SupportFunctions import getLLH_C as getLLH
            for comp in self.LLHTables.keys():
                self.LLHTables[comp] = self.LLHTables[comp].flatten()
        if self.runMode == 'static':
            from SupportFunctions import getLLH_fast as getLLH
        self.getLLH = getLLH
        self.SAT_LG = 90000.

        pass

    def Geometry(self, frame):
        self.geometry = frame['I3Geometry']
        self.PushFrame(frame)

    def Calibration(self, frame):
        self.calibration = frame['I3Calibration']
        self.vemcal_map = self.calibration.vem_cal
        self.PushFrame(frame)

    def DetectorStatus(self, frame):
        self.status = frame['I3DetectorStatus']
        self.tankCheck = True
        index = 0
        self.om2index = {}
        self.tankpositions = []
        self.snowheight = []
        # Build om2index, tank positions, and snowheight lists
        for number, station in self.geometry.stationgeo:
            for tank in station:
                check = False
                for om in tank.omkey_list:
                    if om not in self.status.dom_status.keys():
                        continue
                    if self.status.dom_status[om].pmt_hv < 5e-7:   # Dom is on
                        continue
                    self.om2index[om] = index
                    check = True
                # If any of the doms pass, save necessary data
                if check:
                    self.tankpositions.append(tank.position)
                    self.snowheight.append(tank.snowheight)
                    index += 1

        ntanks = len(self.tankpositions)
        if ntanks < 5:  # or some arbitrary low number
            print('Warning: not enough tanks found in run')
            self.tankCheck = False

        # Convert tank positions and snow heights to arrays
        tankx = [self.tankpositions[i].x for i in range(ntanks)]
        tanky = [self.tankpositions[i].y for i in range(ntanks)]
        tankz = [self.tankpositions[i].z for i in range(ntanks)]
        self.tankxyz = np.array([tankx, tanky, tankz]).T
        self.snowheight = np.asarray(self.snowheight)

        self.PushFrame(frame)

    def Physics(self, frame):

        # Checks necessary for the reconstruction to run
        if not self.tankCheck or not frame.Has(self.DirectionFit):
            self.PushFrame(frame)
            return
        c1 = frame[self.DirectionFit].fit_status.name == 'OK'
        c2 = frame.Has(self.recoPulses)
        if not c1 or not c2:
            self.PushFrame(frame)
            return

        # Load LLH tables
        hists = self.LLHTables
        # binNames = self.bins.keys()
        Ebins = self.bins['E']
        Cbins = self.bins['C']
        cmax = 1.01 * Cbins[-2]

        # Get list of good charges
        VEMcharges = np.zeros(self.snowheight.shape)
        vemPulses = frame[self.recoPulses]
        if vemPulses.__class__ == dc.I3RecoPulseSeriesMapMask:
            vemPulses = vemPulses.apply(frame)
        for om, pulses in vemPulses:
            if om not in self.om2index.keys():
                continue
            charge = pulses[0].charge
            VEMcharges[self.om2index[om]] = charge

            # Deal with saturated DOMs
            if self.saturated:
                if charge != charge:
                    continue
                vemCalib = self.vemcal_map[om]
                pe_per_vem = vemCalib.pe_per_vem / vemCalib.corr_factor
                gain = self.status.dom_status[om].dom_gain_type
                lg_sat = self.SAT_LG / pe_per_vem
                if (charge > lg_sat) and (gain == dc.I3DOMStatus.Low):
                    charge = cmax

        # Shower properties
        goodcharges = np.logical_not(np.isnan(VEMcharges))
        Dir = frame[self.DirectionFit].dir
        theta = Dir.zenith
        phi = Dir.azimuth
        Tank_xyz = self.tankxyz[goodcharges]
        unbinned_vals = {}
        unbinned_vals['E'] = [0] * goodcharges.sum()
        unbinned_vals['Z'] = [theta] * goodcharges.sum()
        unbinned_vals['S'] = self.snowheight[goodcharges]
        unbinned_vals['D'] = [0] * goodcharges.sum()
        unbinned_vals['C'] = VEMcharges[goodcharges]

        # Bin quantities in advance
        binnedVals = [np.digitize(unbinned_vals[i], self.bins[i]) - 1
                      for i in self.bins if self.bins[i] is not None]
        # Skip this if running the static fast version
        if self.runMode == 'flexible':
            binnedVals = np.transpose(binnedVals).flatten().astype(np.int32)

        # Grid search function to run
        def gridLLH(grid, hist):
            RecoMaxLLHs, argmax = self.getLLH(grid, hist, self.bins,
                                              binnedVals, theta, phi, Tank_xyz)
            return RecoMaxLLHs, argmax

        # Perform grid search for each composition
        depthMax = len(self.grids.keys())
        for comp in hists.keys():

            LLHs = [0 for i in range(depthMax)]
            depth = 0
            temp_grid = self.grids[0]
            while depth < depthMax:
                LLHs[depth], Ebin = gridLLH(np.asarray(temp_grid), hists[comp])
                xy_max = temp_grid[LLHs[depth].argmax()]
                depth += 1
                if depth < depthMax:
                    temp_grid = self.grids[depth][xy_max]

            # Setup I3Particle with position and direction
            ML_x, ML_y = xy_max
            ML_z = 1947
            track = dc.I3Particle()
            track.pos.x, track.pos.y, track.pos.z = (ML_x, ML_y, ML_z)
            track.dir.set_theta_phi(theta, phi)
            track.shape = track.InfiniteTrack
            track.fit_status = dc.I3Particle.OK

            # Given most likely position, find most likely energy bin
            new_grid = np.asarray([xy_max])
            maxLLH, Ebin = gridLLH(new_grid, hists[comp])
            maxLLH, Ebin = maxLLH[0], Ebin[0]
            track.energy = b2e(Ebins, Ebin)

            # Params outdated but potentially useful for showing llh-space
            #coarse_llhs = icetray.vector_double()
            #med_llhs = icetray.vector_double()
            #fine_llhs = icetray.vector_double()
            # coarse_llhs.extend(LLHs[0])
            # med_llhs.extend(LLHs[1])
            # fine_llhs.extend(LLHs[2])

            params = ShowerLLHFitParams()
            params.maxLLH = float(maxLLH)
            #params.coarse_llhs = coarse_llhs
            #params.med_llhs = med_llhs
            #params.fine_llhs = fine_llhs

            ## WRITE TO FRAME ##
            frame['ShowerLLH_' + comp] = track
            #frame['maxLLH_'+comp] = dc.I3Double(maxLLH)
            frame['ShowerLLHParams_' + comp] = params

        self.PushFrame(frame)
        return

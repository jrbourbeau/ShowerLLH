#!/usr/bin/env python

#=========================================================================
# File Name     : LLHTable.py
# Description   :
# Creation Date : 06-07-2016
# Created By    : James Bourbeau
#=========================================================================

import os
import re
import numpy as np

from llhbins import LLHBins
from icecube import icetray, phys_services
from icecube import dataclasses as dc

import support_functions.simFunctions as simFunctions


class LLHTable(object):

    def __init__(self, bintype=None, table=None):
        self.bins = LLHBins(bintype=bintype)
        self.count_table = table

    def __str__(self):
        return 'bintype = {}, count_table = {}'.format(self.bins.bintype, self.count_table)

    def make_LLH_tables(self, filelist, outfile):
        ''' Converts the count tables provided in filelist to log-probability tables '''
        if isinstance(filelist, basestring):
            filelist = [filelist]
        cumulative, norm = {}, {}
        # Sum all CountTables
        for i, file in enumerate(filelist):
            # print('Loading '+file)
            counts_table = np.load(file)
            counts_table = counts_table.item()
            # Get composition from simulation
            sim = re.split('_|\.', os.path.basename(file))[1]
            comp = simFunctions.sim2comp(sim, full=True)
            # Add to cumulative table
            if comp not in cumulative.keys():
                cumulative[comp] = np.zeros(
                    counts_table['counts'].shape, dtype=float)
            cumulative[comp] += counts_table['counts']
        # Normalize and log tables (log10(Hist/Hist.sum()))
        for comp in cumulative.keys():
            cumulative[comp] += .1    # baseline so zeros don't give errors
            norm[comp] = np.zeros(cumulative[comp].shape, dtype=float)
            for idx in np.ndindex(cumulative[comp].shape[:-1]):
                norm[comp][idx] = np.log10(
                    cumulative[comp][idx] / sum(cumulative[comp][idx]))

        # Write LLH table to file
        if self.bins is None:
            raise SystemExit(
                'LLH bins are of type NoneType. Please specify them...')
        else:
            np.save(outfile, {'bintype': self.bins.bintype,
                              'llhtables': norm})
            print('{} saved.'.format(outfile))

    def save_tables(self, outfile='~/LLHTable.npy'):
        ''' Write LLH table to file '''
        if self.bins is None or self.llhtable is None:
            raise SystemExit(
                'Either the LLH bins or LLH table are of type NoneType. Please specify them...')
        else:
            np.save(outfile, {'bins': self.bins, 'llhtables': self.llhtable})
            print('{} saved.'.format(outfile))


class FillHist(icetray.I3Module):
    ''' Icetray module to fill counts histograms from simulation '''

    def __init__(self, context):
        icetray.I3Module.__init__(self, context)
        self.AddOutBox('OutBox')
        self.AddParameter('binDict', 'Input LLH bins', None)
        self.AddParameter('bintype', 'Input LLH bin type', None)
        self.AddParameter('recoPulses',
                          'Input RecoPulseSeries to use for reconstruction',
                          'CleanedHLCTankPulses')
        self.AddParameter('outFile', 'Input name for output file', 'test')

    def Configure(self):
        self.bins = self.GetParameter('binDict')
        self.bintype = self.GetParameter('bintype')
        self.recoPulses = self.GetParameter('recoPulses')
        self.outFile = self.GetParameter('outFile')
        pass

    def Geometry(self, frame):
        # Load I3Geometry to be used later
        self.geometry = frame['I3Geometry']
        self.PushFrame(frame)

    def Calibration(self, frame):
        self.calibration = frame['I3Calibration']
        self.PushFrame(frame)

    def DetectorStatus(self, frame):
        # Load I3DetectorStatus to gather positions and snow heights
        # for active tanks
        self.om2index = {}
        self.tankpositions = []
        self.snowheight = []
        self.status = frame['I3DetectorStatus']
        self.dom = self.status.dom_status
        index = 0
        for number, station in self.geometry.stationgeo:
            for tank in station:
                omList = tank.omkey_list
                if not any([om in self.dom.keys() for om in omList]):
                    continue
                for om in omList:
                    self.om2index[om] = index
                self.tankpositions.append(tank.position)
                self.snowheight.append(tank.snowheight)
                index += 1

        if len(self.tankpositions) < 5:  # or some arbitrary low number
            raise ValueError("Not enough tanks found!")

        self.hist = np.zeros(
            [len(self.bins[i]) - 1 for i in self.bins if self.bins[i] is not None])
        self.snowheight = np.array(self.snowheight)
        self.PushFrame(frame)

    def Physics(self, frame):

        # Get shower information
        try:  # the first frame doesn't have a particle?
            # MC = frame['MCPrimary']
            MC = frame['I3MCTree'].primaries[0]
            VEMpulses = frame[self.recoPulses]
            if VEMpulses.__class__ == dc.I3RecoPulseSeriesMapMask:
                VEMpulses = VEMpulses.apply(frame)
        except KeyError:
            self.PushFrame(frame)
            return
        # Closest approach distances
        ClosestDist = phys_services.I3Calculator.closest_approach_distance
        dists = np.array([ClosestDist(MC, xyz) for xyz in self.tankpositions])

        # Saturation values
        sat_lg_pe = 90000.
        vemcal_map = self.calibration.vem_cal

        # Calculate charges in tanks
        VEMcharges = np.zeros(self.snowheight.shape)
        saturated = np.zeros(self.snowheight.shape, dtype='bool')
        for om, pulses in VEMpulses:

            # Skip OMKey(39,62)
            if om not in self.om2index.keys():
                continue

            # Only record first pulse
            charge = pulses[0].charge
            # CleanedHLCTankPulses has one charge per tank
            idx = self.om2index[om]
            VEMcharges[idx] = charge

            # Check saturation
            vemCalib = vemcal_map[om]
            pe_per_vem = vemCalib.pe_per_vem / vemCalib.corr_factor
            lg_sat = sat_lg_pe / pe_per_vem
            gain = self.dom[om].dom_gain_type
            # Note: CleanedHLCTankPulses only returns one pulse per tank. If
            # the high-gain DOM is not saturated, then neither is the low-gain
            # DOM. If the high-gain DOM is saturated, then the low-gain DOM
            # pulse is returned. So we only really need to worry about
            # low-gain DOM saturation.
            if (charge > lg_sat) and (gain == dc.I3DOMStatus.Low):
                saturated[idx] = True

        # Pick out the tanks with valid charges
        goodcharges = np.logical_not(np.isnan(VEMcharges))
        # Set saturated DOMs to fixed value
        cbins = self.bins['C']
        # Make sure they'll fall in last bin
        VEMcharges[goodcharges * saturated] = 1.01 * cbins[-2]

        unbinned_vals = {}
        unbinned_vals['E'] = [np.log10(MC.energy)] * goodcharges.sum()
        unbinned_vals['Z'] = [MC.dir.zenith] * goodcharges.sum()
        unbinned_vals['S'] = self.snowheight[goodcharges]
        unbinned_vals['D'] = dists[goodcharges]
        unbinned_vals['C'] = VEMcharges[goodcharges]

        # Bin shower
        all_unbinned = [unbinned_vals[k]
                        for k in self.bins if self.bins[k] is not None]
        all_edges = [self.bins[i]
                     for i in self.bins if self.bins[i] is not None]
        event_hist = np.histogramdd(all_unbinned, all_edges)[0]
        self.hist += event_hist
        self.PushFrame(frame)
        return

    def Finish(self):
        # Record count tables and bin values in outFile
        if self.outFile:
            d = {}
            d['bintype'] = self.bintype
            d['counts'] = self.hist
            np.save(self.outFile, d)
        return


def merge_counts_tables(filelist, outfile):
    merged_counts_table = {}
    for i, file in enumerate(filelist):
        counts_table = np.load(file)
        counts_table = counts_table.item()
        if i == 0:
            merged_counts_table['bintype'] = counts_table['bintype']
            merged_counts_table['counts'] = np.zeros(
                counts_table['counts'].shape)
        merged_counts_table['counts'] += counts_table['counts']

    np.save(outfile, merged_counts_table)

#!/usr/bin/env python

#=========================================================================
# File Name     : LLHGrid.py
# Description   :
# Creation Date : 06-07-2016
# Created By    : James Bourbeau
#=========================================================================

import numpy as np

from llhbins import LLHBins
from icecube import icetray, phys_services, dataio
from icecube import dataclasses as dc

import support_functions.simFunctions as simFunctions
# from support_functions.llhtools import getVertex


class LLHGrid(object):

    def __init__(self, config=None):
        self.config = config
        self.GCD_file = simFunctions.getGCD(config)
        self.tank_xy = self.get_tank_coords()
        self.search_grid = self.get_search_grids()

    def __str__(self):
        return 'config = {}, GCD_file = {}'.format(self.config, self.GCD_file)

    def get_tank_coords(self):
        tank_xy = {}
        # Extract Geometry and DetectorStatus
        print('Loading {}...'.format(self.GCD_file))
        i3file = dataio.I3File(self.GCD_file)
        finished, geo_check, status_check = False, False, False
        while not finished:
            frame = i3file.pop_frame()
            if frame == None:
                print 'Did not find both Geometry and DetectorStatus in file'
                raise
            if frame.Has('I3Geometry') and not geo_check:
                geo = frame.Get('I3Geometry')
                geo_check = True
            if frame.Has('I3DetectorStatus') and not status_check:
                stat = frame.Get('I3DetectorStatus')
                status_check = True
            if geo_check and status_check:
                finished = True

        # add tanks that are active to tankpositions
        tank_xy = []
        for number, station in geo.stationgeo:
            for tank in station:
                om1, om2 = tank.omkey_list
                dst = stat.dom_status
                if om1 not in dst.keys() or om2 not in dst.keys():
                    continue
                elif dst[om1].pmt_hv > 5e-7 and dst[om2].pmt_hv > 5e-7:
                    tank_xy.append((tank.position.x, tank.position.y))

        return tank_xy

    def get_hex(self, l, nside, theta, xy0):
        ''' Create hex with spacing l, nside spots per side, tilt theta, and center xy0 '''

        # Create a baseline
        x0, y0 = xy0
        xx = np.array([x0 + l * np.cos(theta) *
                       n for n in range(-nside + 1, nside)])
        yy = np.array([y0 + l * np.sin(theta) *
                       n for n in range(-nside + 1, nside)])

        # Expand into a parallelogram
        dx = l * np.cos(theta + np.pi / 3)
        dy = l * np.sin(theta + np.pi / 3)
        xx_bins = np.array([xx + i * dx for i in range(-nside + 1, nside)])
        yy_bins = np.array([yy + i * dy for i in range(-nside + 1, nside)])

        # Refine shape to hexagon
        nmid = 2 * nside - 1
        k = np.array([[i + j for j in range(nmid)] for i in range(nmid)])
        spots = ((k >= nside - 1) * (k <= 2 * nmid - nside - 1))

        return zip(xx_bins[spots], yy_bins[spots])

    def get_search_grids(self):
        ''' Create grid for ShowerLLH iterative grid search '''

        steps = [125, 20, 5]
        nhex = [6 for i in range(len(steps) - 1)]
        boundaryDist = 2.33 * steps[0]
        # steps = kwargs['steps']
        # nhex = kwargs['nhex']

        l = steps[0]    # grid spacing
        Z = 1947        # uniform plane for depth

        tank_x, tank_y = np.transpose(self.tank_xy)

        # pick out the positions of the tanks on the bottom edge and fit
        bot = (tank_y < -390)
        x_bot, y_bot = tank_x[bot], tank_y[bot]
        m, b = np.polyfit(x_bot, y_bot, 1)
        theta = np.arctan(m)

        # Create oversized hexagon centered at (0,0)
        big_nside = int((len(x_bot) / 2 * 125) / l) * 2
        grid_xy = self.get_hex(l, big_nside, theta, [0, 0])

        # Reduce to shape of detector
        grid = {}
        grid[0] = []
        for x, y in grid_xy:
            dmin = np.sqrt((tank_x - x)**2 + (tank_y - y)**2).min()
            if dmin <= boundaryDist:
                grid[0].append((x, y))

        # Recursive function for building nested grids
        def fillGrid(coord, depth, grid):
            hexCoords = self.get_hex(steps[depth], nhex[depth - 1], theta, coord)
            try:
                grid[depth][coord] = hexCoords
            except KeyError:
                grid[depth] = {coord: hexCoords}
            if depth != len(steps) - 1:
                for newCoord in hexCoords:
                    fillGrid(newCoord, depth + 1, grid)

        # Build nested grids
        if len(steps) > 1:
            for coord in grid[0]:
                fillGrid(coord, 1, grid)

        return grid

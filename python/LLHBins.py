#!/usr/bin/env python

#=========================================================================
# File Name     : LLHBins.py
# Description   : Provide a class for ShowerLLH bins
# Creation Date : 06-07-2016
# Created By    : James Bourbeau
#=========================================================================

import collections
import numpy as np

class LLHBins(object):

    def __init__(self, bintype=None):
        self.bintype = bintype
        self.bins = self.getbins(bintype)

    def __str__(self):
        return '( {} | bintype={})'.format(self.bins, self.bintype)

    def getbins(self, bintype, order='EZSDC'):
        """ Create ordered dictionary for likelihood bins """

        if bintype is None:
            return None

        bins = collections.OrderedDict()
        deg_to_rad = np.pi / 180.
        # Energy bins in Log10(Energy/GeV)
        bins['E'] = np.arange(4, 9.501, 0.05)
        # Zenith bins in radians (equal solid-angle bins, cutoff at 40 deg_to_rads)
        bins['Z'] = np.linspace(1, np.cos(40*deg_to_rad), 4)
        bins['Z'] = np.append(np.arccos(bins['Z']), np.pi/2)
        if bintype in ['logdist','nozenith']:
            bins['Z'] = None
        # Snow bins in meters
        bins['S'] = np.array([-1, .001, .5, .85])
        bins['S'] = np.append(bins['S'], np.inf)
        # Distance bins in meters
        bins['D'] = np.append(np.arange(0,600,10), np.arange(600,1051,50))
        if bintype in ['logdist']:
            bins['D'] = 10**np.linspace(np.log10(5), np.log10(1050), 149)
            bins['D'] = np.append(0, bins['D'])
        bins['D'] = np.append(bins['D'], np.inf)
        # Charge bins in VEM
        bins['C'] = 10**np.linspace(-3, 4.5, 46)[10:38]
        bins['C'] = np.append(0, bins['C'])
        bins['C'] = np.append(bins['C'], np.inf)

        return bins

if __name__ == "__main__":
    b = LLHBins(bintype='logdist')
    print('b = {}'.format(b))

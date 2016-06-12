#!/usr/bin/env python

#=========================================================================
# File Name     : LLHTable.py
# Description   :
# Creation Date : 06-07-2016
# Created By    : James Bourbeau
#=========================================================================

import numpy as np
from LLHBins import LLHBins

class LLHTable(object):

    def __init__(self, bintype=None, table=None):
        self.bins = LLHBins(bintype=bintype)
        self.count_table = table
        self.llhtable = None

    def __str__(self):
        return 'LLHBins = {}, count_table = {}'.format(self.llhbins, self.count_table)

    def add_count_table(self, table):
        self.count_table = table

    def counts_to_log_prob(self):
        self.llhtable = []

    def save_tables(self, outfile='~/LLHTable.npy'):
        ''' Write LLH table to file '''
        if self.bins is None or self.llhtable is None:
            raise SystemExit('Either the LLH bins or LLH table are of type NoneType. Please specify them...')
        else:
            np.save(outfile, {'bins': self.bins, 'llhtables': self.llhtable})
            print('{} saved.'.format(outfile))


if __name__ == "__main__":
    a = LLHTable(bintype='standard')
    print('a = {}'.format(a))

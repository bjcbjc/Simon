# -*- coding: utf-8 -*-
"""
Created on Thu Jul 25 14:47:08 2013

@author: bjchen
"""

from numpy import *
from baseStrLib import BaseCountString

class SampleRnaAlleleCount:

    
    def __init__(self, filename):
        """
        data: #record array
        alleleCount;
        """
    
        self.data = loadtxt(filename, delimiter='\t', \
            dtype=dtype([('chrm', 'object'), ('pos', 'int64'), ('ref', 'object'), \
            ('numRead', 'f'), ('baseCountString', 'object')]), converters={0: lambda s: s.strip('chr')})
        self.alleleCount = zeros( (self.data.size, \
            BaseCountString.baseLabels.shape[0], BaseCountString.baseLabels.shape[1] ))
        self.ntlabels = BaseCountString.baseLabels
        
        for vi in range(self.alleleCount.shape[0]):
            self.alleleCount[vi, :, :] = BaseCountString.baseCountStringToMatrix( \
                BaseCountString.baseCountStringRemoveIndels( self.data[vi]['baseCountString']) )
                
    def collapseStrand(self):
        if self.ntlabels.ndim == 2:
            self.ntlabels = self.ntlabels[0]
            self.alleleCount = sum(self.alleleCount, 1)
        
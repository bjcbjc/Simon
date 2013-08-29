# -*- coding: utf-8 -*-
"""
Created on Fri Jul 26 17:14:35 2013

@author: bjchen
"""

#import vcf
import numpy
import re
import traceback

class VCFobj:

    def flattenstruct(self, variantCalls):
        #return dict
        calls = {}
        for vcall in variantCalls:
            calls[vcall.sample] = {}
            for fd, val in vcall.data._asdict().items():
                if val != None:
                    calls[vcall.sample][fd] = val
        return calls

    def parse(self, filename):
        f = open(filename, 'r')
        line= f.readline().strip()
        data = []
        flagfields = []
        strfields = []
        while line:
            if line[0] == '#':
                if line[:6] == '#CHROM':
                    header = line.strip('#').split('\t')
                    fltidx = header.index('FILTER')
                    infoidx = header.index('INFO')
                    fmtidx = header.index('FORMAT')
                    samples = header[fmtidx+1:]
                elif 'INFO=<ID' in line:
                    if 'Type=Flag' in line: flagfields.append(line.split('<ID=')[1].split(',')[0])
                    if 'Type=String' in line or 'Number=.' in line or 'Number=A' in line: 
                        strfields.append(line.split('<ID=')[1].split(',')[0])
                line = f.readline().strip()
                continue            
            line = line.split('\t')
            linedata = line[:fltidx+1]
            info = line[infoidx]
            for s in strfields:
                info = re.sub('%s=(?P<sname>[\w\-\,.]+)'%s, '%s="\g<sname>"'%s, info)
            for s in flagfields:
                info = info.replace(s, '%s=True'%s)
            info = info.replace(';',',').replace('1000genomes.AF', 'kgenomes_AF')
            try:
                linedata.append( eval('dict(%s)'%info) )
            except:
                print info
                traceback.print_exc()
                exit()

            fmt = line[fmtidx].split(':')            
            calls = map(lambda(a): a.split(':'), line[fmtidx+1:])
            calledidx = filter(lambda(a): len(calls[a]) == len(fmt), range(len(samples)))
            linedata.append(fmt)
            linedata.append(calledidx)
            linedata.append( map(lambda(a): calls[a], calledidx) )
            
            data.append(tuple(linedata))
            line = f.readline().strip()
        f.close()
        header = header[:fltidx+1]
        header.extend(['INFO','FORMAT','validSampleIdx', 'callStat'])
        return header, samples, data
    
    def __init__(self, filename):
        
#        f = open(filename, 'r')
#        vobj = vcf.Reader(f)
#        copyfields = filter(lambda(a): a[0] != '_', vars(vobj).keys())
#        for key in copyfields:
#            setattr(self, key, getattr(vobj, key))
#        
        dtype = [('CHROM', 'O'), ('POS', int), ('ID','O'), ('REF','O'), ('ALT', 'O'), \
            ('QUAL', 'O'), ('FILTER','O'), ('INFO','O'), ('FORMAT','O'), \
            ('validSampleIdx', 'O'), ('callStat', 'O')  ]
#        self.records = []
#        for r in vobj:
#            self.records.append( \
#                (r.CHROM, r.POS, r.REF, r.ALT, r.ID, r.QUAL, r.FILTER, r.INFO, self.flattenstruct(r.samples)) )
#        self.records = numpy.array( self.records, dtype=dtype)
#        del vobj  
#        f.close()
        self.header, self.samples, self.records = self.parse(filename)
        self.records = numpy.array(self.records, dtype=dtype)
        
        
        
        
#!/bin/python

from sys import argv, exit
import re
import time

def parseHeader(header):
    header = filter(lambda(l): '<ID=' in l, header)
    token = map(lambda(l): re.search('##(\w+)=<ID=(\w+),Number=([\w.+]),Type=(\w+),Description="([\S\s]+)">', l), header)
    token = filter(lambda(l): l!=None, token)
    token = map(lambda(l): l.groups(), token)
    info = filter(lambda(l): l[0]=='INFO', token)
    format = filter(lambda(l): l[0] == 'FORMAT', token)
    return info, format
        


fn = argv[1]
if len(argv) < 3:
    outfnhead = fn
else:
    outfnhead = argv[2]

t = open(fn).readlines()
header = filter(lambda(l):l[0]=='#', t)
text = [ l.strip('\n').split('\t') for l in filter(lambda(l): l[0]!='#', t) ]
infoheader, formatheader = parseHeader(header)

tic = time.time()
columnheader = header[-1][1:].strip('\n').split('\t')
infoindex = columnheader.index('INFO')
infolines = [ l[infoindex] for l in text ]
#infodata = map(lambda(l): dict(re.findall('(\w+)=([^;]+)',l)), infolines )
prog = re.compile('(\w+)=([^;]+)')
infodata = [ dict( prog.findall(l) ) for l in infolines ]
formatdata = [l[infoindex+1:] for l in text]
print 'parse info, %s'%(time.time()-tic)

numheader = '%s\t%s'%(columnheader[1], columnheader[5])
numlines = [ '\t'.join([l[1], l[5] if l[5]!='.' else 'NaN']) for l in text ]
textheader = '\t'.join(map(lambda(i): columnheader[i], [0, 2, 3, 4, 6]))
textlines = [ '\t'.join([l[0], l[2], l[3], l[4], l[6]]) for l in text ]

del text

nvariant = len(numlines)
tic = time.time()
#info
for dummy1, name, num, datatype, dummy2 in infoheader:
    if num.isdigit():
        num = int(num)
    else:
        num = -1
    if datatype == 'Flag':
        miskey = '0'
        addcol = [ '1' if name in l else miskey for l in infolines ] #map(lambda(l): '1' if name in l else miskey, infolines)
        numeric = True
        if num == 0: 
            num = 1
    elif datatype == 'String':
        miskey = ''
        addcol = [ l[name] if name in l.keys() else miskey for l in infodata ] #map(lambda(l): l[name] if name in l.keys() else miskey, infodata)
        numeric = False
    else:
        if num != -1:
            miskey = '\t'.join( ['NaN']*num )
            addcol = [ l[name].replace(',','\t') if name in l.keys() else miskey for l in infodata ] #map(lambda(l): l[name].replace(',','\t') if name in l.keys() else miskey, infodata)
            numeric = True
        else:
            miskey = ''
            addcol = [ l[name] if name in l.keys() else miskey for l in infodata ] #map(lambda(l): l[name] if name in l.keys() else miskey, infodata)
            numeric = False
    if numeric:
        if not all( l==miskey for l in addcol):
            numheader = numheader + '\t' + '\t'.join( [name]*num)
            #numlines = [ '%s\t%s'%l for l in zip(numlines,addcol)]
            for i in xrange(nvariant):
                numlines[i] = numlines[i] + '\t' + addcol[i]
    else:
        if not all( l==miskey for l in addcol):
            textheader = textheader + '\t' + name
            #textlines = [ '%s\t%s'%l for l in zip(textlines,addcol)]
            for i in xrange(nvariant):
                textlines[i] = textlines[i] + '\t' + addcol[i]

print 'add info columns, %s'%(time.time()-tic)

del infolines
del infodata

#format
samplename = columnheader[ columnheader.index('FORMAT')+1: ]
samplename = [ l.replace('-','_') for l in samplename ]
nsample = len(samplename)
if nsample != len(formatdata[0])-1:
    print 'nsample inconst, exit'
    exit()

tic = time.time()
#formatdata = map(lambda(l): dict( zip(l[0].split(':'), zip(*map(lambda(x): x.split(':'), l[1:])))), formatdata)
formatdata = [ dict( zip(l[0].split(':'), zip(*[ x.split(':') for x in l[1:]]))) for l in formatdata ]
print 'parse format, %s'%(time.time()-tic)

tic = time.time()
for dummy1, name, num, datatype, dummy2 in formatheader:
    if num.isdigit():
        num = int(num) * nsample
    else:
        num = -1    

    if datatype == 'String':
        miskey = ''
        addcol = ['\t'.join(l[name]) if name in l.keys() else miskey for l in formatdata ] #map(lambda(l): '\t'.join(l[name]) if name in l.keys() else miskey, formatdata)
        numeric = False
    else:
        if num != -1:
            miskey = '\t'.join( ['NaN']*num )
            addcol = [ '\t'.join(l[name]).replace(',','\t') if name in l.keys() else miskey for l in formatdata ] #map(lambda(l): '\t'.join(l[name]).replace(',','\t') if name in l.keys() else miskey, formatdata)
            numeric = True
        else:
            newnum = set([ x.count(',') for l in formatdata if name in l for x in l[name] if x!= '.'  ]) #count from only non empty
            if len(newnum) == 1:
                numPerSamp = newnum.pop() + 1
                num = numPerSamp * nsample
                miskey = '\t'.join( ['NaN']*num )
                miskeyPerSamp = '\t'.join( ['NaN']*numPerSamp )
                addcol = [ '\t'.join([ x if x!='.' else miskeyPerSamp for x in l[name]]).replace(',','\t') if name in l else miskey for l in formatdata]
#                addcol = [ ('\t'.join(l[name]).replace(',','\t') if l['name']!='.' else miskey) if name in l.keys() else miskey for l in formatdata ] 
                numeric = True
            else:
                miskey = ''
                addcol = [ '\t'.join(l[name]) if name in l.keys() else miskey for l in formatdata ] #map(lambda(l): l[name] if name in l.keys() else miskey, formatdata)
                numeric = False

    samphead  = [ b for a in samplename for b in [a]*(num/nsample) ]

#    if name=='BQ':
#        print numeric, num, miskey
#        print formatdata[0][name]
#        print addcol[0]
#        exit()

    if numeric:
        if not all( l==miskey for l in addcol):
            numheader = numheader + '\t' + '\t'.join( [ a+'_'+b for a,b in zip(samphead, [name]*num) ] )
            #numlines = [ '%s\t%s'%l for l in zip(numlines,addcol)]
            for i in xrange(nvariant):
                numlines[i] = numlines[i] + '\t' + addcol[i]
    else:
        if not all( l==miskey for l in addcol):
            if num == -1: 
                num = nsample
                samphead  = [ b for a in samplename for b in [a]*(num/nsample) ]
            textheader = textheader + '\t' + '\t'.join( [ a+'_'+b for a,b in zip(samphead, [name]*num) ] )
#            if name=='BQ':
#                print num
#                print samphead
#                print textheader
#                exit()
            #textlines = [ '%s\t%s'%l for l in zip(textlines,addcol)]
            for i in xrange(nvariant):
                textlines[i] = textlines[i] + '\t' + addcol[i]

print 'add format columns, %s'%(time.time()-tic)

f = open(outfnhead + '.mtx', 'w')
f.write(numheader + '\n')
for line in numlines:
    f.write(line + '\n')
f.close()

f = open(outfnhead + '.cell', 'w')
f.write(textheader + '\n')
for line in textlines:
    f.write(line + '\n')
f.close()


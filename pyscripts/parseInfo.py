

from sys import argv
import os

fn = argv[1]
outputdir = argv[2]
if len(argv) > 3: 
    maxalt = int(argv[3]) - 1
else:
    maxalt = 0

if outputdir[-1] != '/': outputdir = outputdir + '/'

t = open(fn).readlines()
headers = filter(lambda(l): l[0] == '#', t)

fmtnames = []
infonum = []

for line in headers:
    if 'INFO=' in line:
        infonames.append( line.split('ID=')[1].split(',')[0] )
        expectnumdata = line.split('Number=')[1].split(',')[0]
        datatype = line.split('Type=')[1].split(',')[0]

        infoexpand.num(False)
        if datatype in ['Integer', 'Float', 'Flag']:
            if expectnumdata == '1':
                infoexpand[-1] = True

headers = headers[-1].strip().split('\t')
t = map(lambda(l): l.strip().split('\t'), filter(lambda(l): l[0]!='#', t))

infoidx = headers.index('INFO')
infohash = {}
#### unfinished
infocolumn = map(lambda(l): l[infoidx].split(':'), t)
numattr = map(len, infocolumn)
t = map(lambda(l): l[infoidx+1:], t)

nsample = len(t[0])
rowidx = range(len(t))
flatdata = []
for sampi in range(nsample):
    infodata = map(lambda(l): l[sampi].split(':'), t)
    numdata = map(len, infodata)
    validrowidx = filter(lambda(i): numdata[i] == numattr[i], rowidx)
    # varidx, sampidx, fieldname, fielddata
    map(lambda(i): flatdata.extend( map(lambda(a): [i+1, sampi+1, a[0], a[1]], zip(infocolumn[i], infodata[i])) ), validrowidx)

print 'len(flatdata)= %d'%len(flatdata)
#write out tab file for each attribute: [varidx, sampidx, data], #loc x #sample
for fdidx in range(len(infonames)):
    field = infonames[fdidx]
    subtable = filter(lambda(l): l[2] == field, flatdata)
    if len(subtable) > 0:
        f = open('%s%s.%s'%(outputdir, os.path.basename(argv[1]), field.replace('.','_')), 'w')
        if infoexpand[fdidx]:
            for line in subtable:
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t')))
        elif infoexpandByAllele[fdidx]:
            for line in subtable:
                nalt = line[3].count(',')
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t') + '\tNaN'*(maxalt-nalt)) )
        elif infoexpandByAllelePlusRef[fdidx]:
            for line in subtable:
                nalt = line[3].count(',')
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t') + '\tNaN'*(maxalt-nalt+1)) )
        else:
            for line in subtable:
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3]))
        f.close()

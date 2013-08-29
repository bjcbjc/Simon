

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
fmtexpand = []
fmtexpandByAllele = []
fmtexpandByAllelePlusRef = []
for line in headers:
    if 'FORMAT=' in line:
        fmtnames.append( line.split('ID=')[1].split(',')[0] )
        expectnumdata = line.split('Number=')[1].split(',')[0]
        datatype = line.split('Type=')[1].split(',')[0]

        fmtexpand.append(False)
        fmtexpandByAllele.append(False)
        fmtexpandByAllelePlusRef.append(False)

        if datatype in ['Integer', 'Float', 'Flag']:
            if expectnumdata.isdigit():
                fmtexpand[-1] = True
            elif  expectnumdata == 'A' and maxalt != 0:
                fmtexpandByAllele[-1] = True
            elif 'allelic depth' in line.lower() and maxalt != 0:
                fmtexpandByAllelePlusRef[-1] = True


headers = headers[-1].strip().split('\t')
t = map(lambda(l): l.strip().split('\t'), filter(lambda(l): l[0]!='#', t))

fmtidx = headers.index('FORMAT')
fmtcolumn = map(lambda(l): l[fmtidx].split(':'), t)
numattr = map(len, fmtcolumn)
t = map(lambda(l): l[fmtidx+1:], t)

nsample = len(t[0])
rowidx = range(len(t))
flatdata = []
for sampi in range(nsample):
    fmtdata = map(lambda(l): l[sampi].split(':'), t)
    numdata = map(len, fmtdata)
    validrowidx = filter(lambda(i): numdata[i] == numattr[i], rowidx)
    # varidx, sampidx, fieldname, fielddata
    map(lambda(i): flatdata.extend( map(lambda(a): [i+1, sampi+1, a[0], a[1]], zip(fmtcolumn[i], fmtdata[i])) ), validrowidx)

print 'len(flatdata)= %d'%len(flatdata)
#write out tab file for each attribute: [varidx, sampidx, data], #loc x #sample
for fdidx in range(len(fmtnames)):
    field = fmtnames[fdidx]
    subtable = filter(lambda(l): l[2] == field, flatdata)
    if len(subtable) > 0:
        f = open('%s%s.%s'%(outputdir, os.path.basename(argv[1]), field.replace('.','_')), 'w')
        if fmtexpand[fdidx]:
            for line in subtable:
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t')))
        elif fmtexpandByAllele[fdidx]:
            for line in subtable:
                nalt = line[3].count(',')
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t') + '\tNaN'*(maxalt-nalt)) )
        elif fmtexpandByAllelePlusRef[fdidx]:
            for line in subtable:
                nalt = line[3].count(',')
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3].replace(',','\t') + '\tNaN'*(maxalt-nalt+1)) )
        else:
            for line in subtable:
                f.write('%s\t%s\t%s\n'%(line[0], line[1], line[3]))
        f.close()

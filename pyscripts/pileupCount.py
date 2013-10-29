#!/usr/bin/python

import argparse
import subprocess
import traceback
import re
from collections import Counter, OrderedDict
from randstr import randstr
from baseStrLib import BaseString
import time

def chrmReplaceFunc(matchobj):
    if matchobj.group('chrm') == 'MT': return 'chrM'
    elif matchobj.group('chrm') == 'chrM': return 'MT'
    elif 'chr' not in matchobj.group('chrm'): return 'chr'+matchobj.group('chrm')
    else: return matchobj.group('chrm')[3:]
    

    
def createTempLocFile(args, validChrms):
    #if vcf file, create a temporary file
    f = open(args.f)
    chrmformat = f.readline().strip('\n>').split()[0]
    f.close()

    if chrmformat not in validChrms:
        print 'inconsistent chrm format, reference: %s, bam: %s'%(chrmformat, validChrms[0])
        exit(1)
        
    f = open(args.loc)
    line = f.readline()
    f.close()
    removelater = []
    if 'fileformat=VCF' in line:
        locfile = args.tmpdir + randstr() + '.loc'
        #create temporary location file
        subprocess.call('grep -v ^# %s | cut -f1,2 > %s'%(args.loc, locfile), shell=True)
        removelater.append( locfile )
    elif len(line.split('\t')) > 2:
        locfile = args.tmpdir + randstr() + '.loc'
        #create temporary location file
        subprocess.call('cut -f1,2 %s > %s'%(args.loc, locfile), shell=True)
        removelater.append( locfile )
    else:
        locfile = args.loc

    #check location chrm format and filter unavailable contigs
    f = open(locfile)
    loc = f.readlines()
    f.close()
    locChrms = list(set(map(lambda(l): l.split()[0], loc)))    
    
    bamWithChr, locWithChr = False, False
    reprog = re.compile('^chr')
    if any(map(lambda(l): reprog.match(l), validChrms)):
        bamWithChr = True
        replaceProg = re.compile('(?P<chrm>^[0-9]{1,2}|^[XYMT]{1,2})')
    if any(map(lambda(l): reprog.match(l), locChrms)):
        locWithChr = True
        replaceProg = re.compile('(?P<chrm>^chr')
    if bamWithChr != locWithChr:
        loc = [ re.sub(replaceProg, chrmReplaceFunc, line) for line in loc]

    loc = filter(lambda(line): line.split()[0] in validChrms, loc)
    #write another temporary file
    locfile = args.tmpdir + randstr() + '.loc'
    print 'tmpfile: %d lines, '%len(loc), locfile
    f = open(locfile, 'w')
    for line in loc: f.write(line)
    f.close()
    removelater.append( locfile )
    
    return locfile, removelater

def pileupToCount(samout, out, minread, stranded=True, tableformat=False):
    #pileup special chars in read_bases: ^: start of read, followed by an additional char for mapping qual
    # $: end of read
    # > or < suggests non-coding region
    excludePattern = re.compile('(\^.)|(\$)')
    if stranded:
        baseLabels = ['>', 'A', 'T', 'C', 'G', 'N', '<', 'a', 't', 'c', 'g', 'n']
    else:
        baseLabels = ['>', 'A', 'T', 'C', 'G', 'N']

    if tableformat:
        out.write('\t'.join(['chr', 'pos', 'ref', '#read']) + '\t' + '\t'.join(baseLabels) + '\n')
    line = samout.readline()
    count = 1
    outbuffer = ''
    while line:
        line = line.split() #chr, pos, ref, #read, read_bases, read_qual

        if int(line[3]) >= minread:
            if not stranded: 
                line[4] = line[4].upper().replace(',','.').replace('<','>')
            
            line[4] = line[4].replace('.', line[2].upper()).replace(',', line[2].lower())
            indels, line[4] = BaseString.matchIndels(line[4])

            #count
            baseCount = Counter( re.sub(excludePattern, '', line[4] ))   #dict with counts and letters
            baseCount.update(indels)

            if tableformat:
                countstring = '\t'.join( [ '%s'%baseCount[b] for b in baseLabels ] )
                #out.write('\t'.join(line[:4]) + '\t' + countstring + '\n')
                outbuffer = outbuffer + '\t'.join(line[:4]) + '\t' + countstring + '\n'
            else:
                baseCount = OrderedDict( sorted(baseCount.items(), key=lambda t:t[0]) )
                #decide what to write to the output: chr, pos, ref, total_read, base_count_string
                baseCountString = ''.join([ k+'%s'%v for k, v in baseCount.iteritems()])
                #out.write('\t'.join(line[:4]) + '\t' + baseCountString + '\n')            
                outbuffer = outbuffer + '\t'.join(line[:4]) + '\t' + baseCountString + '\n'
        count = count + 1
        if count%100 == 0:
            out.write(outbuffer)
            outbuffer = ''
        line = samout.readline()
    if outbuffer != '':
        out.write(outbuffer)
    return

def bamChrms(args):
    sampipe = subprocess.Popen(args.samtools + ' view -H ' + args.bam, shell=True, stdout=subprocess.PIPE)
    samheader = sampipe.stdout.read().split('\n')
    sampipe.stdout.close()
    samheader = filter(lambda(l): l[:3]=='@SQ', samheader)
    chrms = map(lambda(l): l.split()[1].replace('SN:',''), samheader)
    return chrms

def closefiles(files):
    for f in files:
        if f != None: f.close()
    return

argp = argparse.ArgumentParser(prog='pileupCount.py')
argp.add_argument('-d', default=10000, metavar='INT', type=int, help='INT for samtools mpileup, read maximally INT reads per bam [10,000]' )
argp.add_argument('-q', default=0, metavar='INT',  type=int, help='INT for samtools mpileup, min mapping quality for an alignment to be used [0]' )
argp.add_argument('-Q', default=13, metavar='INT', type=int, help='INT for samtools mpileup, mim base quality for a base to be considered [13]' )
argp.add_argument('-f', type=str, default='/nethome/bjchen/Projects/Simon/h37_hg19_chrOnly.fa', metavar='file', help='faidx-indexed reference fasta' )
argp.add_argument('-loc', type=str, metavar='file', default='', help='file, file can be VCF or bed, [chr, pos] as columns')
argp.add_argument('-reg', type=str, metavar='file', default='',help='file, file contains regions for allele counts; mutual exclusive with -loc; one output file will be generated by each region')
argp.add_argument('-bam', type=str, metavar='file', required=True, help='file, bam file for pileup')
argp.add_argument('-r',  metavar='INT', required=True, type=int, default=10, help='min reads required for each base')
argp.add_argument('-o', type=str, default='test.out', metavar='file', help='file name for the output, [<bam>.out]')
argp.add_argument('-tableformat', action='store_true', help='make table output instead; this option is turned off by default')
argp.add_argument('-stranded', action='store_true', help='a switch to specify if it is strand-specific data, [false]' )
argp.add_argument('-samtools', type=str, default='/data/NYGC/Software/samtools/samtools-0.1.19/samtools', help='samtools path, [/data/NYGC/Software/samtools/samtools-0.1.19]')
argp.add_argument('-tmpdir', type=str, default='./', help='path for temporary files, [./]')
#argp.add_argument('-refbasecheck', nargs='?', action='store_true', help='a switch to turn on double checking of reference bases in vcf and in genome; only valid if location file is in vcf format, [false]')
argp.add_argument('-samtools_mpileup_option', nargs='?', metavar='value', help='options that are passed to samtools mpileup')

args = argp.parse_known_args()
passon = ''
if isinstance(args, tuple):
    passon = ' '.join(args[1])
    args = args[0]
if args.tmpdir[-1] != '/': args.tmpdir = args.tmpdir + '/'


argtb = vars(args)
sampara = ['d', 'q', 'Q', 'f']

if args.loc != '' and args.reg != '':
    print 'Both -loc and -reg specified; ignore -loc.'
    args.loc = ''
elif args.loc == '' and args.reg == '':
    print 'Must specify -loc or -reg'
    exit()

validChrms = bamChrms(args)
if args.reg != '':
    regions = [l.strip() for l in open(args.reg).readlines()]
    if any([ l.split(':')[0] not in validChrms for l in regions ]):
        print 'Some chromosomes specified in regions are not found in bam'
        exit()
    cmd = args.samtools + ' mpileup -r %s '
else:
    locfile, tmpfiles = createTempLocFile(args, validChrms)
    cmd = args.samtools + ' mpileup -l %s '%(locfile)

for switch in sampara:
    cmd = cmd + '-%s %s '%(switch, argtb[switch])

cmd = cmd + ' ' + passon + ' ' + args.bam    
print 'Run command: ' + cmd

samout, samerr, out = None, None, None
try:
    #run samtools and pipe the output for further process (save space)
    if args.loc != '':
        sampipe = subprocess.Popen(cmd, shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
        samout = sampipe.stdout #, sampipe.stderr

        out = open(args.o, 'w')
        pileupToCount(samout, out, args.r, stranded=args.stranded, tableformat=args.tableformat)
        
        if len(tmpfiles) > 0:
            #remove temp file
            subprocess.call( 'rm -f ' + ' '.join(tmpfiles), shell=True )
        closefiles([samout, samerr, out])
    else: #regions; call mpileup with one region at a time; this is faster than -l
        for subr in regions:
            sampipe = subprocess.Popen(cmd%(subr), shell=True, stdout=subprocess.PIPE)#, stderr=subprocess.PIPE)
            samout = sampipe.stdout #, sampipe.stderr

            subr = subr.replace(',','')
            outfntag = '.'.join( re.split('[:-]', subr) )
            out = open(args.o + '.' + outfntag, 'w')
            pileupToCount(samout, out, args.r, stranded=args.stranded, tableformat=args.tableformat)
            closefiles([samout, samerr, out])
except:
    if len(tmpfiles) > 0:
        #remove temp file
        subprocess.call(  'rm -f ' + ' '.join(tmpfiles), shell=True )
    closefiles([samout, samerr, out])
    traceback.print_exc()


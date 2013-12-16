from sys import argv
from os import system, popen
from string import Template
import os.path
import copy
import cmdGenerator
import jobFactory
import configRobot
import itertools
import re
from randstr import randstr


matlabcmd = '/usr/local/MATLAB/R2013b/bin/matlab -nodisplay -r "addpath(genpath(\'/nethome/bjchen/BJLib/Matlabox\')); addpath(genpath(\'%s/functions\')); %s; quit"\n'

#read all available parameters for programs
def getAvailableParas():
    t = open('/nethome/bjchen/BJLib/pylib/qsub/AllAvailableParas.txt').readlines()
    t = map(lambda(l):l.strip(), t)
    t = map(lambda(l):l.split('\t'), t)
    tb = {}
    for l in t:
        tb[ l[0] ] = l[1:]
    return tb


def filterVCFByBed(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    intersect =  'bedops -e -1  <(cat {vcfpath}{vcf} | /data/NYGC/Software/python/bin/python2.7 /data/NYGC/Software/bedops/bin/vcf2bed -) <( gawk \'$4>={mincov} {{print $0}}\' {bedpath}{bed}  ) > {output}.cov{mincov}.bed'

    #intersect = '/data/NYGC/Software/bedtools/bedtools-2.17.0/bin/bedtools intersect -a {vcfpath}{vcf} -b <(gawk \'{{OFS="\\t"; if($4 >={mincov}){{ print $0}}}}\' {bedpath}{bed}) -u -sorted > {output}.filtered.0.vcf'
    union = 'sort -u {output1}.cov{mincov}.bed {output2}.cov{mincov}.bed > {output}.cov{mincov}.bed'
    intersect2 = 'bedops -e -1 <( cat {output1}.cov{mincov}.bed) <( cat {output2}.cov{mincov}.bed) > {output}.intersect.cov{mincov}.bed'

    cmdset = configRobot.makeParasList(cmdset, ['vcf','bed', 'mincov'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    vcfs, beds = configRobot.popParas(cmdset, ['vcf', 'bed'])
    mincovs = cmdset.pop('mincov')

    mem = cmdset['mem']
    sgeopt = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for vcf, bed in zip(vcfs, beds): #one job per group
        bed = [l.strip() for l in bed.split(',')]
        cmdset['vcf'] = vcf

        jobprefix = prefix + '_' + vcf
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, vcf, randstr()]) + '/'
        outputbase = sampletmpoutpath + vcf
        CMDs = []
        if 'toShell' in cmdset.keys():
            CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )
        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
        for mincov in mincovs:
            cmdset['mincov'] = mincov
            for i in range(1,len(bed)+1):
                cmdset['bed'] = cmdset['bedname'].format(bed=bed[i-1])
                cmdset['output%d'%i] = outputbase + '_' + cmdset['bed'] 
                CMDs.append( intersect.format(output=cmdset['output%d'%i], **cmdset ) )
            if len(bed) > 1:
                CMDs.append( union.format( output=outputbase, **cmdset ) )
                CMDs.append( intersect2.format( output=outputbase, **cmdset) )
        CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}*.bed {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, outputpath=outputpath) )) 
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager

def snpeff(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    snpeffcmd = '{java} -Xmx{mem} -jar {snpeffpath}snpEff.jar eff -v -i vcf -o vcf -stats {output}.snpEff_summary.html -noLog -c {snpeffpath}/snpEff.config {annotation} {input}.variants.filtered.vcf > {output}.variants.snpEff.vcf'

    cmdset = configRobot.makeParasList(cmdset, ['vcf','sample'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    vcfs, bams = configRobot.popParas(cmdset, ['vcf', 'bam'])

    mem = cmdset['mem']
    vcfname = cmdset.pop('vcfname')
    bamname = cmdset.pop('bamname')
    outputname = cmdset.pop('output')
    if 'num_core' not in cmdset: cmdset['num_core'] = 1
    if int(cmdset['num_core']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['num_core']]
    else:
        sgeopt = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for vcf, bam in zip(vcfs, bams): #one job per group
        cmdset['vcf'] = vcf
        cmdset['bam'] = bam
        cmdset['inputvcf'] = inputpath + vcfname.format( **cmdset )
        cmdset['inputbam'] = inputpath + bamname.format( **cmdset )
        jobprefix = prefix + outputname.format(**cmdset)
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, outputname.format(**cmdset), randstr()]) + '/'
        cmdset['output'] = sampletmpoutpath + outputname.format( **cmdset )
        CMDs = []
        if 'toShell' in cmdset.keys():
            CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )
        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
        for subcmd in [vcfconvert, rmmismatch, rmrepeat, rmintron, rmhomopoly, blat, rmedit]:
            CMDs.append( subcmd.format( **cmdset ) )
        CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}*txt {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, outputpath=outputpath) )) 
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager


def snpirfilter_virmidall(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    catvcf1 = '''egrep -v ^# {vcfpath}{vcf}.germ.virmid.vcf | gawk '{{OFS="\\t"; $7=$7"_GERM"; print $0}}' > {tmpdata}'''
    catvcf2 = '''egrep -v ^# {vcfpath}{vcf}.loh.virmid.vcf | gawk '{{OFS="\\t"; $7=$7"_LOH"; print $0}}' >> {tmpdata}'''
    catvcf3 = 'egrep -v ^# {vcfpath}{vcf}.som.virmid.vcf  >> {tmpdata}'
    headvcf = 'egrep ^# {vcfpath}{vcf}.som.virmid.vcf > {vcfpath}{vcf}.all.virmid.vcf'
    mergevcf = 'sort -k1,1 -k2,2n {tmpdata} >> {vcfpath}{vcf}.all.virmid.vcf'

    vcf_somatic_caller_convert = 'bash {snpirpath}/convertVCF_somatic_caller.sh {inputvcf} {output}.txt {callername}'
    rmmismatch = '{perl} {snpirpath}/filter_mismatch_first6bp.pl -infile {output}.txt -outfile {output}.rmhex.txt -bamfile {inputbam}'
    rmrepeat = '''awk '{{OFS="\\t";$2=$2-1"\\t"$2;print $0}}' {output}.rmhex.txt | {bedtoolspath}/intersectBed -a stdin -b {repeatmask} -v | cut -f1,3-7 > {output}.rmhex.rmsk.txt'''
    rmintron = '{perl} {snpirpath}/filter_intron_near_splicejuncts.pl -infile {output}.rmhex.rmsk.txt -outfile {output}.rmhex.rmsk.rmintron.txt -genefile {geneannotation}'
    rmhomopoly = '{perl} {snpirpath}/filter_homopolymer_nucleotides.pl -infile {output}.rmhex.rmsk.rmintron.txt -outfile {output}.rmhex.rmsk.rmintron.rmhom.txt -refgenome {reference} '
    blat = '{perl} {snpirpath}/BLAT_candidates.pl -infile {output}.rmhex.rmsk.rmintron.rmhom.txt -outfile {output}.rmhex.rmsk.rmintron.rmhom.rmblat.txt -bamfile {inputbam} -refgenome {reference} '
    rmedit = '''awk '{{OFS="\\t";$2=$2-1"\\t"$2;print $0}}' {output}.rmhex.rmsk.rmintron.rmhom.rmblat.txt | {bedtoolspath}/intersectBed -a stdin -b {rnaedit} -v > {output}.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed'''

    filters = [catvcf1, catvcf2, catvcf3, headvcf, mergevcf, vcf_somatic_caller_convert, rmmismatch, rmrepeat, rmintron, rmhomopoly, blat, rmedit]

    if 'removeedit' not in cmdset:
        cmdset['removeedit'] = True
    else:
        if cmdset['removeedit'] == 'True': cmdset['removeedit'] = True
        else: cmdset['removeedit'] = False
    if not cmdset['removeedit']:
        filters = filters[:-1]

    cmdset = configRobot.makeParasList(cmdset, ['vcf','sample', 'vcfname'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    vcfpath = cmdGenerator.checkPath(cmdset.pop('vcfpath'))
    cmdset['vcfpath'] = vcfpath
    bampath = cmdGenerator.checkPath(cmdset.pop('bampath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))#, create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    vcfs, bams = configRobot.popParas(cmdset, ['vcf', 'bam'])

    mem = cmdset['mem']
    vcfnames = cmdset.pop('vcfname')
    bamname = cmdset.pop('bamname')
    outputname = cmdset.pop('output')
    if 'num_core' not in cmdset: cmdset['num_core'] = 1
    if int(cmdset['num_core']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['num_core']]
    else:
        sgeopt = []

    nvcf = len(vcfnames)
    if 'callername' in cmdset:
        callernames = cmdset.pop('callername')
        if type(callernames) == type('caller'): callernames = [callernames]*nvcf
    else:
        callernames = ['']*nvcf

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for vcf, bam in zip(vcfs, bams): 
        for vcfname, callername in zip(vcfnames, callernames): #one job per group
            CMDs = []
            cmdset['callername'] = callername
            cmdset['vcf'] = vcf
            cmdset['bam'] = bam
            cmdset['vcfname'] = vcfname.format( **cmdset)
            cmdset['inputvcf'] = vcfpath.format( **cmdset) + vcfname.format( **cmdset )
            cmdset['inputbam'] = bampath.format( **cmdset)  + bamname.format( **cmdset )
            jobprefix = prefix + outputname.format(**cmdset)
            sampletmpoutpath = tmpoutpath + '_'.join([prefix, outputname.format(**cmdset), randstr()]) + '/'
            cmdset['tmpdata'] = sampletmpoutpath + 'tmpdata'
            cmdset['output'] = sampletmpoutpath + outputname.format( **cmdset )
            if 'toShell' in cmdset.keys():
                CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )
            CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
            for subcmd in filters:
                CMDs.append( subcmd.format( **cmdset ) )
            CMDs.append( cmdGenerator.formatCmd( 'rm -f {tmpdata}'.format(**cmdset)) )
            CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}* {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, outputpath=outputpath.format(**cmdset)) )) 
            CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath.format(**cmdset))) )

            jobmanager.createJob(jobprefix, CMDs, outpath = outputpath.format(**cmdset), outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager



def snpirfilter(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    outputnametb = {'vcfconvert': '{output}.txt', 'rmmismatch':'{output}.rmhex.txt', 'rmrepeat':'{output}.rmhex.rmsk.txt', 'rmintron':'{output}.rmhex.rmsk.rmintron.txt', 'rmhomopoly':'{output}.rmhex.rmsk.rmintron.rmhom.txt', 'blat':'{output}.rmhex.rmsk.rmintron.rmhom.rmblat.txt', 'rmedit':'{output}.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed'}
    vcfconvert = 'bash {snpirpath}/convertVCF.sh {inputvcf} {output}.txt {varqual} '
    vcf_somatic_caller_convert = 'bash {snpirpath}/convertVCF_somatic_caller.sh {inputvcf} {output}.txt {callername}'
    rmmismatch = '{perl} {snpirpath}/filter_mismatch_first6bp.pl -infile {vcfconvert} -outfile {output}.rmhex.txt -bamfile {inputbam}'
    rmrepeat = '''awk '{{OFS="\\t";$2=$2-1"\\t"$2;print $0}}' {output}.rmhex.txt | {bedtoolspath}/intersectBed -a stdin -b {repeatmask} -v | cut -f1,3-7 > {output}.rmhex.rmsk.txt'''
    rmintron = '{perl} {snpirpath}/filter_intron_near_splicejuncts.pl -infile {output}.rmhex.rmsk.txt -outfile {output}.rmhex.rmsk.rmintron.txt -genefile {geneannotation}'
    rmhomopoly = '{perl} {snpirpath}/filter_homopolymer_nucleotides.pl -infile {output}.rmhex.rmsk.rmintron.txt -outfile {output}.rmhex.rmsk.rmintron.rmhom.txt -refgenome {reference} '
    blat = '{perl} {snpirpath}/BLAT_candidates.pl -infile {output}.rmhex.rmsk.rmintron.rmhom.txt -outfile {output}.rmhex.rmsk.rmintron.rmhom.rmblat.txt -bamfile {inputbam} -refgenome {reference} '
    rmedit = '''awk '{{OFS="\\t";$2=$2-1"\\t"$2;print $0}}' {output}.rmhex.rmsk.rmintron.rmhom.rmblat.txt | {bedtoolspath}/intersectBed -a stdin -b {rnaedit} -v > {output}.rmhex.rmsk.rmintron.rmhom.rmblat.rmedit.bed'''

    filters = [vcfconvert, rmmismatch, rmrepeat, rmintron, rmhomopoly] #, blat, rmedit]
    if 'somatic_caller' not in cmdset:
        cmdset['somatic_caller'] = False
    else:
        if cmdset['somatic_caller'] == 'True': cmdset['somatic_caller'] = True
        else: cmdset['somatic_caller'] = False
    if 'removeedit' not in cmdset:
        cmdset['removeedit'] = True
    else:
        if cmdset['removeedit'] == 'True': cmdset['removeedit'] = True
        else: cmdset['removeedit'] = False
    if 'runBLAT' not in cmdset:
        cmdset['runBLAT'] = True
    else:
        if cmdset['runBLAT'] == 'True': cmdset['runBLAT'] = True
        else: cmdset['runBLAT'] = False
    if cmdset['somatic_caller']:
        filters[0] = vcf_somatic_caller_convert

    if cmdset['runBLAT']: filters.append(blat)
    if cmdset['removeedit']: filters.append(rmedit)

    filternameset = ['vcfconvert', 'rmmismatch', 'rmrepat', 'rmintron', 'rmhomopoly', 'blat', 'rmedit']
    if 'reStartAt' in cmdset:
        reStartIdx = filternameset.index(cmdset['reStartAt'])
        filters =  filters[reStartIdx:] 

    cmdset = configRobot.makeParasList(cmdset, ['vcf','sample', 'vcfname'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    vcfpath = cmdGenerator.checkPath(cmdset.pop('vcfpath'))
    bampath = cmdGenerator.checkPath(cmdset.pop('bampath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))#, create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    vcfs, bams = configRobot.popParas(cmdset, ['vcf', 'bam'])

    mem = cmdset['mem']
    vcfnames = cmdset.pop('vcfname')
    bamname = cmdset.pop('bamname')
    outputname = cmdset.pop('output')
    if 'num_core' not in cmdset: cmdset['num_core'] = 1
    if int(cmdset['num_core']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['num_core']]
    else:
        sgeopt = []

    nvcf = len(vcfnames)
    if 'callername' in cmdset:
        callernames = cmdset.pop('callername')
        if type(callernames) == type('caller'): callernames = [callernames]*nvcf
    else:
        callernames = ['']*nvcf

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for vcf, bam in zip(vcfs, bams): 
        for vcfnamepattern, callername in zip(vcfnames, callernames): #one job per group
            if '*' in vcfnamepattern:
                f = popen('ls %s'%(vcfpath.format(**cmdset)+vcfnamepattern.format(vcf=vcf)))
                vcfnamelist = [l.replace(vcfpath.format(**cmdset),'') for l in f.read().split()]
                f.close()
            else:
                vcfnamelist = [vcfnamepattern]
            for vcfname in vcfnamelist:
                vcfname = vcfname.replace('.rmhex.rmsk.rmintron.rmhom.txt','')
                CMDs = []
                cmdset['callername'] = callername
                cmdset['vcf'] = vcf
                cmdset['bam'] = bam
                cmdset['vcfname'] = vcfname.format( **cmdset)
                cmdset['inputvcf'] = vcfpath.format( **cmdset) + vcfname.format( **cmdset )
                cmdset['inputbam'] = bampath.format( **cmdset)  + bamname.format( **cmdset )
                jobprefix = prefix + outputname.format(**cmdset)
                sampletmpoutpath = tmpoutpath + '_'.join([prefix, outputname.format(**cmdset), randstr()]) + '/'
                cmdset['output'] = sampletmpoutpath + outputname.format( **cmdset )
                if 'toShell' in cmdset.keys():
                    CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )
                CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
                if 'reStartAt' in cmdset:
                    if reStartIdx > 0:
                        CMDs.append( 'cp %s%s %s'%(outputpath.format(**cmdset), outputnametb[filternameset[reStartIdx-1]].format(output=outputname.format(**cmdset)), sampletmpoutpath) )
                for subcmd in filters:
                    CMDs.append( subcmd.format( **cmdset ) )
                CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}* {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, outputpath=outputpath.format(**cmdset)) )) 
                CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath.format(**cmdset))) )

                jobmanager.createJob(jobprefix, CMDs, outpath = outputpath.format(**cmdset), outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager



def GATK_genotyper(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    gatk = '{java} -Djava.io.tmpdir={javatmpdir} -Xmx{mem} -jar {gatkjar} -R {reference} '
    genotyper = '-T UnifiedGenotyper -U ALLOW_N_CIGAR_READS {input} -o {output}.vcf --dbsnp {DBSNP} -nt {num_core} -nct {num_thread_core} -glm {-glm} -stand_call_conf {-stand_call_conf} -stand_emit_conf {-stand_emit_conf} '

    cmdset = configRobot.makeParasList(cmdset, ['group'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    groups, groupnames = configRobot.popParas(cmdset, ['group', 'groupname'])

    mem = cmdset['mem']
    bamname = cmdset.pop('input')

    if 'num_core' not in cmdset: cmdset['num_core'] = 1
    if int(cmdset['num_core']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['num_core']]
    else:
        sgeopt = []

    if 'passon' not in cmdset:
        cmdset['passon'] = ''

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for sample, name in zip(groups, groupnames): #one job per group
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, name, randstr()]) + '/'
        jobprefix = prefix + name
        cmdset['input'] = ' '.join( ['-I ' + inputpath + bamname.format(sample=x) for x in re.split('[,\s]+', sample)] )
        cmdset['output'] = sampletmpoutpath + name
        cmdset['javatmpdir'] = sampletmpoutpath
        CMDs = []
        if 'toShell' in cmdset.keys():
            CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )

        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
        CMDs.append( cmdGenerator.formatCmd( (gatk+genotyper).format(**cmdset) + cmdset['passon']) )
        CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}*vcf* {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, outputpath=outputpath) )) 
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager


def recalAlign(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    ftread = 'samtools view -bh -q {-q} -F {-F} {input} -o {output}.flt.bam '
    bamindex = 'samtools index {output}.flt.bam '
    gatk = '{java} -Djava.io.tmpdir={javatmpdir} -Xmx{mem} -jar {gatkjar} -R {reference} '
    indeltarget = '-T RealignerTargetCreator -U ALLOW_N_CIGAR_READS -I {output}.flt.bam -o {output}.indel.intervals -known {KGMILLS} -known {KGINDEL} -nt {num_thread} '
    indelrealign = '-T IndelRealigner -U ALLOW_N_CIGAR_READS -I {output}.flt.bam -known {KGMILLS} -known {KGINDEL} --targetIntervals {output}.indel.intervals -o {output}.flt.realigned.bam ' #-compress 0 '
    baserecal = '-T BaseRecalibrator -U ALLOW_N_CIGAR_READS -I {output}.flt.realigned.bam -knownSites {DBSNP} -knownSites {KGMILLS} -knownSites {KGINDEL} -o {output}.recal.table -nct {num_thread} '
    baserecalbam = '-T PrintReads -U ALLOW_N_CIGAR_READS -I {output}.flt.realigned.bam --BQSR {output}.recal.table -o {output}.flt.realigned.recal.bam -nct {num_thread}'

    cmdset = configRobot.makeParasList(cmdset, ['sample'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    samples = configRobot.popParas(cmdset, ['sample'])
    mem = cmdset['mem']
    input = cmdset.pop('input')

    if 'num_thread' not in cmdset: cmdset['num_thread'] = '1'
    if int(cmdset['num_thread']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['num_thread']]
    else:
        sgeopt = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for sample in samples: #one job per sample
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, sample, randstr()]) + '/'
        jobprefix = prefix + sample
        cmdset['input'] = input.format(sample=inputpath+sample)
        cmdset['output'] = sampletmpoutpath + input.format(sample=sample).replace('.bam','')
        cmdset['javatmpdir'] = sampletmpoutpath

        CMDs = []
        if 'toShell' in cmdset.keys():
            CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )

        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )

        CMDs.append( ftread.format(**cmdset) )
        CMDs.append( bamindex.format(**cmdset) )
        CMDs.append( (gatk+indeltarget).format(**cmdset) )
        CMDs.append( (gatk+indelrealign).format(**cmdset) )
        CMDs.append( (gatk+baserecal).format(**cmdset) )
        CMDs.append( (gatk+baserecalbam).format(**cmdset) )

        CMDs.append( cmdGenerator.formatCmd( 'mv {sampletmpoutpath}{sample}*realigned*  {sampletmpoutpath}{sample}*.intervals {sampletmpoutpath}{sample}*.table {outputpath}'.format(sampletmpoutpath=sampletmpoutpath, sample=sample, outputpath=outputpath) )) 
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager


def bwaalign(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    fastqstream = '<(zcat {fastq} | sed -e "s/ \([12]\):\([YN]\)/:\\1 \\1:\\2/")' #retain pair info
    bwaaln = ' {bwaprog} aln -Y -t {-t} {reference} {fastq} > {outputfile}.sai'
    bwasamse = ' {bwaprog} samse -n 4 '
    bwasamse_input = ' {reference} {outputfile}.sai {fastq} '
    coorconvert = '{java} -Djava.io.tmpdir={javatmpdir} -Xmx{mem} -cp /data/NYGC/Resources/SNPiR/ convertCoordinates - | samtools view -Sb - > {outputfile}.cord.bam'
    picard = '{java} -Djava.io.tmpdir={javatmpdir} -Xmx{mem} -jar {picardpath}{picardprog} TMP_DIR={tmpdir} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true '
    picard_addRG = ' RGID={fqid} RGLB={sample} RGPL=Illumina RGPU={barcode} RGSM={sample} RGCN=NYGC RGDS={fqinfo} SORT_ORDER=coordinate INPUT={outputfile}.cord.bam OUTPUT={outputfile}.sorted.bam '
    picard_merge = ' OUTPUT={sample}.sorted.bam ASSUME_SORTED=true USE_THREADING=true '
    picard_mdup = ' INPUT={sample}.sorted.bam OUTPUT={sample}.longheader.bam METRICS_FILE={sample}.sorted.mdup.metrics ASSUME_SORTED=true '
    reheader = 'samtools reheader <(samtools view -H {sample}.longheader.bam | egrep -v ''chr\w+-'') {sample}.longheader.bam > {sample}.sorted.mdup.bam'
    samindex = 'samtools index {sample}.sorted.mdup.bam'


    cmdset = configRobot.makeParasList(cmdset, ['sample'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    samples = configRobot.popParas(cmdset, ['sample'])
    mem = cmdset['mem']

    if '-t' not in cmdset: cmdset['-t'] = 1
    if int(cmdset['-t']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['-t']]
    else:
        sgeopt = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for sample in samples: #one job per sample
        sampleinputpath = inputpath.format(sample=sample)
        f = popen('ls %s*.fastq.gz'%sampleinputpath)
        readfns = f.read().split()
        f.close()
        fqid = 0
        jobprefix = prefix + sample
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, sample, randstr()]) + '/'
        cmdset['javatmpdir'] = sampletmpoutpath
        CMDs = []
        if 'toShell' in cmdset.keys():
            CMDs.append( cmdGenerator.formatCmd( cmdset['toShell'] ) )

        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )

        allbams = []
        for fastq in readfns:
            tokens = fastq.split('/')[-1].split('_')
            barcode = tokens[1]
            fqinfo = '_'.join(tokens[2:]).replace('.fastq.gz','')
            fqid += 1
            outputfile = sampletmpoutpath + fastq.split('/')[-1].replace('.fastq.gz', '')
            CMDs.append( bwaaln.format(fastq=fastqstream.format(fastq=fastq), outputfile=outputfile, **cmdset) )
            CMDs.append( (bwasamse + bwasamse_input + ' | ' + coorconvert).format( \
                        outputfile=outputfile, fastq=fastqstream.format(fastq=fastq), **cmdset) )
#            CMDs.append( bwasamse.format(**cmdset) + \
#                         bwasamse_input.format(outputfile=outputfile, fastq=fastqstream.format(fastq=fastq), **cmdset) + \
#                         ' | ' + \
#                         coorconvert.format(outputfile=outputfile, **cmdset) )
            CMDs.append( picard.format(picardprog='AddOrReplaceReadGroups.jar', tmpdir=sampletmpoutpath, **cmdset) + \
                       picard_addRG.format(fqid=fqid, sample=sample, barcode=barcode, fqinfo=fqinfo, outputfile=outputfile) )
            allbams.append(outputfile +'.sorted.bam')

        CMDs.append( picard.format(picardprog='MergeSamFiles.jar', tmpdir=sampletmpoutpath, **cmdset) + \
                   picard_merge.format(sample=sampletmpoutpath+sample) + ' INPUT=' + ' INPUT='.join(allbams) )
        CMDs.append( picard.format(picardprog='MarkDuplicates.jar', tmpdir=sampletmpoutpath, **cmdset) + \
                     picard_mdup.format(sample=sampletmpoutpath+sample) )
        CMDs.append( reheader.format(sample=sampletmpoutpath+sample) )
        CMDs.append( samindex.format(sample=sampletmpoutpath+sample) )

        CMDs.append( cmdGenerator.formatCmd( 'mv %s%s.sorted.mdup* %s'%(sampletmpoutpath, sample, outputpath) )) 
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager
                    

def varPositionInRead(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, prefix, bam, vcf = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix', 'bam', 'vcf'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath) 
    bampath = cmdGenerator.checkPath(cmdset.pop('bampath'))
    vcfpath = cmdGenerator.checkPath(cmdset.pop('vcfpath'))
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    pycmd = '/data/NYGC/Software/python/Python-2.7.3/python ' + cmdGenerator.checkPath(cmdset.pop('programpath')) + cmdset.pop('pyprog')
    
    if type(bam) != type([]) and type(bam) != type(()):
        bam = [bam]
    if type(vcf) != type([]) and type(vcf) != type(()):
        vcf = [vcf]

    for sampi in range(len(bam)):
        paraset = copy.deepcopy(cmdset)        
        bamfn = bam[sampi]
        vcffn = vcf[sampi]
        jobtag = vcffn.replace('_Aligned.out.WithReadGroup.sorted','').replace('.bam','') + '_' + bamfn.replace('_Aligned.out.WithReadGroup.sorted','')
        jobfnprefix = prefix + '_' + jobtag

        outputfile = jobtag + '.varreadpos.gz'

        CMDs = []
        sampletmpoutpath = tmpoutpath + randstr() + '/'
        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
    
        paraset['-o'] =  sampletmpoutpath + outputfile
        paraset['-bam'] = bampath + bamfn
        paraset['-vcf'] = vcfpath + vcffn

        CMDs.append( cmdGenerator.formatCmd( pycmd, paraset ) )
        CMDs.append( cmdGenerator.formatCmd( 'mv %s %s'%(paraset['-o'], outputpath) ) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobfnprefix, CMDs, outpath = outputpath, outfn = jobfnprefix, trackcmd=False)
    return jobmanager


def bowtie(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    bowtie_setting = ' --local --fast-local '

    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    prog = cmdset.pop('prog')

    if int(cmdset['-p']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['-p']]
    else:
        sgeopt = []


    outputbam = cmdset.pop('bamfn')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    f = popen('ls %s*R1*.fastq.gz'%inputpath)
    readfns = f.read().split()
    f.close()

    sampletmpoutpath = tmpoutpath + prefix + '_' + randstr() + '/'
    jobprefix = prefix 
    
    CMDs = []
    if 'toShell' in cmdset.keys():
        CMDs.append( cmdGenerator.formatCmd( cmdset.pop('toShell') ) )

    CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )

    callcmd = prog + bowtie_setting + '-1 %s -2 %s'%(','.join(readfns), ','.join(readfns).replace('R1','R2'))
    callcmd = cmdGenerator.formatCmd(callcmd, cmdset) + ' | samtools view -Sb - > %s'%(sampletmpoutpath + outputbam)
    CMDs.append(callcmd)

    CMDs.append( cmdGenerator.formatCmd( 'mv %s* %s'%(sampletmpoutpath, outputpath) )) 
    CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

    jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeopt=sgeopt)
    return jobmanager


def mismatchCountByRead(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath) 
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    pycmd = '/data/NYGC/Software/python/Python-2.7.3/python ' + cmdGenerator.checkPath(cmdset.pop('programpath')) + cmdset.pop('pyprog')
    
    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    #expand if needed
    if any( ['*' in l for l in samples] ):
        fstr = samples
        samples = []
        for fkey in fstr:
            f = popen('ls %s%s'%(inputpath, fkey))
            samples.extend([ l.replace(inputpath, '') for l in  f.read().split()])
            f.close()

    for sampi in range(len(samples)):
        sample = samples[sampi]
        paraset = copy.deepcopy(cmdset)        
        jobfnprefix = prefix + '_' + sample
        if bam == '=sample':
            inputfile = sample
            outputfile = sample + '.mismatch.stat'
        else: 
            print 'not implemented'
            exit()

        CMDs = []
        sampletmpoutpath = tmpoutpath + randstr() + '/'
        CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
    
        paraset['-o'] =  sampletmpoutpath + outputfile
        paraset['-bam'] = inputpath + inputfile
        
        CMDs.append( cmdGenerator.formatCmd( pycmd, paraset ) )
        CMDs.append( cmdGenerator.formatCmd( 'mv %s %s'%(paraset['-o'], outputpath) ) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobfnprefix, CMDs, outpath = outputpath, outfn = jobfnprefix, trackcmd=False)
    return jobmanager


def replaceMQ(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    cmd, mem, time, prefix, addext = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix', 'addext'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    bampattern, replaceval = configRobot.popParas(cmdset, ['bampattern', 'replaceval'])
    if 'sampara' in cmdset.keys():
        sampara = cmdset.pop('sampara')
        if sampara[0] != '"': sampara = '"' + sampara + '"'
    else:
        sampara = ''
    f = popen('ls %s%s'%(inputpath, bampattern))
    bam = map(lambda(l): l.split('/')[-1], f.read().split())
    f.close()

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    for bamfile in bam:
        jobprefix = prefix + bamfile
        CMDs = []
        jobtmppath = tmpoutpath + randstr() + '/'

        inputfile = inputpath + bamfile
        outputfile = outputpath + bamfile.replace('.bam', '.%s.bam'%(addext))
        tmpoutfile = jobtmppath + bamfile.replace('.bam', '.%s.bam'%(addext))

        if 'toShell' in cmdset.keys():
            CMDs.append( cmdset['toShell'] ) 

        CMDs.append( cmdGenerator.checkPathOnNode(jobtmppath) )
        CMDs.append( cmdGenerator.formatCmd( 'bash /nethome/bjchen/BJLib/shlib/replaceFieldBam.sh', inputfile, '5', replaceval, sampara, '>', tmpoutfile) )
        CMDs.append( cmdGenerator.formatCmd('mv %s %s'%(tmpoutfile, outputfile)) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath)) )

        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False)
    return jobmanager


def batchreadvcf(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    vcffnkey = cmdset.pop('vcffnkey')
    if ' ' in vcffnkey:
        tmp = vcffnkey.split(' ')
        vcffnkey = "{'" + "','".join(tmp) + "'}"
    else:
        vcffnkey = "'" + vcffnkey + "'"

    vcfpath = cmdGenerator.checkPath(cmdset.pop('vcfpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    matlabworkpath = cmdGenerator.checkPath(cmdset.pop('matlabworkpath'))
    call, njob, outfnhead = configRobot.popParas(cmdset, ['call', 'njob', 'outfnhead'])
    njob = int(njob)
    if 'runOnServer' in cmdset.keys():
        runOnServer = configRobot.popParas(cmdset, ['runOnServer'])
    else:
        runOnServer = ''

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for jobidx in range(1, njob+1):
        jobprefix = prefix + '%02d'%jobidx
        CMDs = []
        functioncall = call + "(%d, %d, %s, '%s', '%s', '%s');"%(jobidx, njob, vcffnkey, vcfpath,outputpath,outfnhead)

        CMDs.append( cmdGenerator.formatCmd( matlabcmd%(matlabworkpath, functioncall) ) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, outputpath) ) )
        jobmanager.createJob(jobprefix, CMDs, outpath = outputpath, outfn = jobprefix, trackcmd=False, sgeJob=False, runOnServer=runOnServer)
    return jobmanager

def varCall_samtools(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['sample'])
    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    caller = cmdset.pop('caller')
    bam = cmdset.pop('bam')
    sample = cmdset.pop('sample') #list of samples

    if 'iterparas' in cmdset.keys():
        iterparas = cmdset.pop('iterparas')
    else:
        iterparas = ["''"]

    if 'specifyname' in cmdset.keys():
        cmdset = configRobot.makeParasList(cmdset, ['specifyname'])
        specifyname = configRobot.popParas(cmdset, ['specifyname'])
        if len(specifyname) != len(sample):
            print '#specifyname should be the same as #sample lists'
            exit()
    else:
        specifyname = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    callcmd = programpath + caller
    if 'mpileup' not in callcmd: callcmd = callcmd + ' mpileup '

    for sampidx in range(len(sample)):
        if bam == '=sample':
            sampleinputpath = inputpath
            sampleoutputpath = outputpath
        else:
            sampleinputpath = inputpath + sample + '/'
            sampleoutputpath = outputpath + sample + '/'
            print "not implemented"

        sampletmpoutpath = tmpoutpath + prefix + '/'

        #if multiple samples, add input path for each sample
        if sample[sampidx].count('.bam') > 1:
            bamfiles = ' '.join( [sampleinputpath + x for x in re.split('[,\s]+', sample[sampidx])] )
        else:
            bamfiles = sampleinputpath + sample[sampidx]

        for addpara in iterparas:
            paraset = copy.deepcopy(cmdset)
            if len(specifyname) == 0:
                jobprefix = prefix + '_' + re.sub('[,\s]+', '_', sample[sampidx])
                if len(jobprefix) > 60:
                    print 'jobprefix is too long, use specifyname instead'
                    exit()
            else:
                jobprefix = prefix + specifyname[sampidx]

            if addpara != "''": 
                paraset[addpara] = ''
                jobprefix = jobprefix + '_' + addpara.replace(' ','')

            CMDs = []
            if 'toShell' in paraset.keys():
                CMDs.append( cmdGenerator.formatCmd( paraset.pop('toShell') ) )

            CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )

            #eg. samtools mpileup -DSugf ref.fa aln.bam | bcftools view -vcg - > var.vcf

            CMDs.append( cmdGenerator.formatCmd( programpath+caller, paraset, bamfiles, ' | bcftools view -vcg - > %s%s.vcf'%(sampletmpoutpath, jobprefix)) )

            CMDs.append( cmdGenerator.formatCmd( 'mv %s.vcf %s'%(sampletmpoutpath+jobprefix, sampleoutputpath) )) 
            CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, sampleoutputpath)) )

            jobmanager.createJob(jobprefix, CMDs, outpath = sampleoutputpath, outfn = jobprefix, trackcmd=False)
    return jobmanager



def varCall(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    callertag = cmdset.pop('callertag').lower()
    if callertag == 'samtools':
        jobmanager = varCall_samtools(cmdset, runmode)
        return jobmanager

    #caller sample and output format string
    #flag: options before template or not
    callertb = {'mutect': [Template(' --input_file:normal $inputpath$normalsample --input_file:tumor $inputpath$tumorsample --out $outputpath$outbase.txt --coverage_file $outputpath$outbase.coverage.wig.txt --vcf $outputpath$outbase.vcf '), False], 'varscan': [Template('samtools mpileup -B -f $reference -q $minq $inputpath$normalsample > $outputpath$normalsample.pileup \n\nsamtools mpileup -B -f $reference -q $minq $inputpath$tumorsample > $outputpath$tumorsample.pileup\n\n$callcmd somatic $outputpath$normalsample.pileup $outputpath$tumorsample.pileup --output-snp $outputpath$outbase.snv --output-indel $outputpath$outbase.indel --output-vcf '), False], 'sniper': [Template(' $inputpath$tumorsample $inputpath$normalsample $outputpath$outbase.vcf '), True], 'virmid': [Template(' -D $inputpath$tumorsample -N $inputpath$normalsample -w $outputpath -o $outbase '), False], 'strelka':[Template(' --tumor $inputpath$tumorsample --normal $inputpath$normalsample --ref $reference --config $strelka_config --output-dir $outputpath \n\nmake -C $outputpath'), False] }

    cmdset = configRobot.makeParasList(cmdset, ['normalsample', 'tumorsample'])
    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    caller = cmdset.pop('caller')
    bam = cmdset.pop('bam')
    normalsample = cmdset.pop('normalsample')
    tumorsample = cmdset.pop('tumorsample')

    if 'iterparas' in cmdset.keys():
        iterparas = cmdset.pop('iterparas')
    else:
        iterparas = ["''"]

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '.jar' in caller:
        callcmd = '$java -version \n\n$java -Xmx%dg '%(int(mem.replace('G',''))-1) 
        if '-Djava.io.tmpdir' in cmdset.keys():
            callcmd = callcmd + ' -Djava.io.tmpdir%s '%(configRobot.popParas(cmdset, ['-Djava.io.tmpdir']))
        callcmd = callcmd + ' -jar ' + programpath + caller
    else:
        callcmd = programpath + caller

    template = callertb[callertag][0]
    optionsbeforetemp = callertb[callertag][1]

    for sampidx in range(len(normalsample)):
        if bam == '=sample':
            sampleinputpath = inputpath
            sampleoutputpath = outputpath
        else:
            sampleinputpath = inputpath + sample + '/'
            sampleoutputpath = outputpath + sample + '/'

        sampletmpoutpath = tmpoutpath + prefix + '_' + randstr() + '/'

        for addpara in iterparas:
            paraset = copy.deepcopy(cmdset)
            tempset = {}
            tempset['inputpath'] = sampleinputpath
            tempset['outputpath'] = sampletmpoutpath
            tempset['normalsample'] = normalsample[sampidx]
            tempset['tumorsample'] = tumorsample[sampidx]
            if callertag == 'varscan':
                tempset['reference'] = paraset.pop('-f')
                tempset['minq'] = paraset.pop('-q')
                tempset['callcmd'] = callcmd
                callcmd = ''
            elif callertag == 'strelka':
                tempset['reference'] = paraset.pop('--ref')
                tempset['strelka_config'] = paraset.pop('--config')
            jobprefix = prefix + '_' + normalsample[sampidx] + '_' + tumorsample[sampidx]

            if addpara != "''": 
                paraset[addpara] = ''
                jobprefix = jobprefix + '_' + addpara.replace(' ','')

            tempset['outbase'] = jobprefix + '.' + callertag

            CMDs = []
            if 'toShell' in paraset.keys():
                CMDs.append( cmdGenerator.formatCmd( paraset.pop('toShell') ) )

            if callertag != 'strelka':
                CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )
            else:
                CMDs.append( cmdGenerator.formatCmd('echo "TMPJOBDIR=%s"'%(sampletmpoutpath)))

            if callertag == 'virmid':
                CMDs.append( cmdGenerator.formatCmd('cp %s %s'%(sampleinputpath+tempset['normalsample'], sampletmpoutpath)) ) 
                CMDs.append( cmdGenerator.formatCmd('cp %s %s'%(sampleinputpath+tempset['tumorsample'], sampletmpoutpath)) )
                tempset['inputpath'] = sampletmpoutpath

            if optionsbeforetemp:
                varcallcmds = cmdGenerator.formatCmd( callcmd, paraset, template.substitute(tempset) ) 
            else:
                varcallcmds = cmdGenerator.formatCmd( callcmd, template.substitute(tempset), paraset ) 
            for vcmd in varcallcmds.split('\n\n'):
                CMDs.append( vcmd )

            if callertag == 'strelka':
                CMDs.append( cmdGenerator.formatCmd( 'mv %sresults/all.somatic.snvs.vcf  %s'%(tempset['outputpath'], sampleoutputpath+tempset['outbase']+'.snv.vcf') )) 
                CMDs.append( cmdGenerator.formatCmd( 'mv %sresults/all.somatic.indels.vcf  %s'%(tempset['outputpath'], sampleoutputpath+tempset['outbase']+'.indel.vcf') )) 
            elif callertag == 'virmid':
                CMDs.append( cmdGenerator.formatCmd( 'mv %s*.report %s'%(tempset['outputpath'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd( 'mv %s*.germ %s'%(tempset['outputpath'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd( 'mv %s*.gm %s'%(tempset['outputpath'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd( 'mv %s*.som %s'%(tempset['outputpath'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd( 'mv %s*.loh %s'%(tempset['outputpath'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd('rm -f %s'%(sampletmpoutpath+tempset['normalsample'])) ) 
                CMDs.append( cmdGenerator.formatCmd('rm -f %s'%(sampletmpoutpath+tempset['tumorsample'])) )
            else:
                CMDs.append( cmdGenerator.formatCmd( 'mv %s* %s'%(tempset['outputpath']+tempset['outbase'], sampleoutputpath) )) 
                CMDs.append( cmdGenerator.formatCmd( 'rm -f %s*.pileup'%(tempset['outputpath']) ) )
            CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobprefix, jobmanager.ext, sampleoutputpath)) )

            jobmanager.createJob(jobprefix, CMDs, outpath = sampleoutputpath, outfn = jobprefix, trackcmd=False)
            if callcmd == '' and 'callcmd' in tempset.keys():
                callcmd = tempset['callcmd']
    return jobmanager




def pileupCount(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['-r','-d','-q','stranded'])
    cmd, mem, time, samples, prefix, bam, loc = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam', 'locfile'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath) 
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    locfilepath = cmdGenerator.checkPath(cmdset.pop('locfilepath'))
    cmdset['-tmpdir'] = cmdGenerator.checkPath(cmdset.pop('-tmpdir'))

    opt_r, opt_d, opt_q, opt_stranded = configRobot.popParas(cmdset, ['-r', '-d', '-q', 'stranded'])

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    setuppathcmd = cmdGenerator.formatCmd('source ./pathsetup')
    pycmd = '/data/NYGC/Software/python/Python-2.7.3/python ' + cmdGenerator.checkPath(cmdset.pop('programpath')) + cmdset.pop('pyprog')
    
    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]

    iterlist = list(itertools.product( opt_r, opt_d, opt_q, opt_stranded ))

    for sampi in range(len(samples)):
        sample = samples[sampi]
        paraset = copy.deepcopy(cmdset)        
        paraset['-loc'] = locfilepath + loc[sampi]
        
        for r, d, q, stranded in iterlist:
            sampleRunStr = sample + '.r%s.d%dk.q%s.s%s'%(r, int(d)/1000, q, stranded[0])
            jobfnprefix = prefix + '_' + sampleRunStr
            if bam == '=sample':
                inputfile = sample
                outputfile = sample + sampleRunStr + '.count'
            else: 
                print 'not implemented'
                exit()
            
            paraset['-r'] = r
            paraset['-d'] = d
            paraset['-q'] = q
            if stranded == 'True': paraset['-stranded'] = ' '
            paraset['-o'] = paraset['-tmpdir'] + outputfile
            paraset['-bam'] = inputpath + inputfile
            countcmd = cmdGenerator.formatCmd( pycmd, paraset )
            mvfilecmd = cmdGenerator.formatCmd( 'mv %s %s'%(paraset['-o'], paraset['-o'].replace(paraset['-tmpdir'], outputpath) ))
            mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath))

            jobmanager.createJob(jobfnprefix, [setuppathcmd, countcmd, mvfilecmd, mvscriptcmd], outpath = outputpath, outfn = jobfnprefix)
    return jobmanager



def markDup(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, shellsetup = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix','toShell'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'), create=createpath)
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))
    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-1)
    mdupjar = 'MarkDuplicates.jar'

    for sample in samples:
        jobname = prefix + '_' + sample
        paraset = copy.deepcopy(cmdset)
        if bam == '=sample':
            sampleinputpath = inputpath
            sampletmpoutpath = tmpoutpath
            sampleoutputpath = outputpath
            bamfile = sample
        else:
            sampleinputpath = inputpath + sample
            sampletmpoutpath = cmdGenerator.checkPath(tmpoutpath + sample, create=createpath)
            sampleoutputpath = cmdGenerator.checkPath(outputpath + sample, create=createpath)
            bamfile = bam

        paraset['INPUT'] = '=%s/%s'%(sampleinputpath, bamfile)
        paraset['OUTPUT'] = '=%s/%s'%(sampletmpoutpath, bamfile.replace('.bam', '.mdup.bam'))
        paraset['METRICS_FILE'] = '=%s/%s'%(sampletmpoutpath, bamfile.replace('.bam', '.mdup.metrics.txt'))
        paraset = configRobot.validParas(paraset, availParas[mdupjar])
        CMDs = []
        
        CMDs.append( cmdGenerator.formatCmd( shellsetup ) )
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+mdupjar, paraset) )
        CMDs.append( cmdGenerator.formatCmd( 'mv %s* %s'%(paraset['OUTPUT'].replace('=','').replace('.bam',''), sampleoutputpath )))
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, sampleoutputpath)) )
        jobmanager.createJob(jobname, CMDs, outpath = sampleinputpath, outfn = jobname)
    return jobmanager





###main
def main(cmd = '', module = '', config = 'projconfig.txt', runmode = 'test'):
    global availParas
    cmdsets = configRobot.readConfig(config)
    setnames = cmdsets.keys()
    availParas = getAvailableParas()

    for setname in setnames:
        if setname != module and cmdsets[setname]['cmd'] != cmd:
            continue

        jobmanager = eval( '%s(cmdsets[setname], runmode)'%cmdsets[setname]['cmd'] )

        if runmode == 'run' and jobmanager != '':
            submitted, skipped = jobmanager.submitJob()
            jobmanager.removeJobFn(status='skipped')
            del jobmanager
            jobmanager = ''


if __name__ == '__main__':
    main( **dict( map(lambda(l): l.split('='), argv[1:]) ) )

    #usage:
    #python rnaseq-pipeline.py <option>=<val>
    #options:
    #  cmd = cmd defined in each set in configure file
    #  module = setname defined in configure file
    #  runmode = {'run', 'test'}; 'test' will only generate script files; 'run' will generate script files and submit the jobs
    #
    #  cmd and module are used to specify which command and/or which module(set) to run
    #






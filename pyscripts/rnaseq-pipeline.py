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


def recalAlign(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    ftread = 'samtools view -bh -q {-q} -F {-F} {input} -o {output}.flt.bam '
    bamindex = 'samtools index {output}.flt.bam '
    gatk = '{java} -Xmx{mem} -jar {gatkjar} -R {reference} '
    indeltarget = '-T RealignerTargetCreator -I {output}.flt.bam -o {output}.indel.intervals -known {KGMILLS} -known {KGINDEL} '
    indelrealign = '-T IndelRealigner -I {output}.flt.bam -known {KGMILLS} -known {KGINDEL} --targetIntervals {output}.indel.intervals -o {output}.flt.realigned.bam -compress 0 '
    baserecal = '-T BaseRecalibrator -I {output}.flt.realigned.bam -knownSites {DBSNP} -knownSites {KGMILLS} -knownSites {KGINDEL} -o {output}.recal.table -nct {num_thread} '
    baserecalbam = '-T PrintReads -I {output}.flt.realigned.bam --BQSR {output}.recal.table -o {output}.flt.realigned.recal.bam -nct {num_thread}'

    cmdset = configRobot.makeParasList(cmdset, ['sample'])
    prefix, time  = configRobot.popParas(cmdset, ['prefix', 'time'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    tmpoutpath = cmdGenerator.checkPath(cmdset.pop('tmpoutpath'))
    samples = configRobot.popParas(cmdset, ['sample'])
    mem = cmdset['mem']
    input = cmdset.pop('input')

    if '-t' not in cmdset: cmdset['-t'] = 1
    if int(cmdset['-t']) > 1: #multiple threads per job
        sgeopt = ['-pe make ' + cmdset['-t']]
    else:
        sgeopt = []

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    for sample in samples: #one job per sample
        sampletmpoutpath = tmpoutpath + '_'.join([prefix, sample, randstr()]) + '/'
        jobprefix = prefix + sample
        cmdset['input'] = input.format(sample=inputpath+sample)
        cmdset['output'] = sampletmpoutpath + input.format(sample=sample).replace('.bam','')

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

        CMDs.append( cmdGenerator.formatCmd( 'mv %s%s.realigned* %s'%(sampletmpoutpath, sample, outputpath) )) 
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
    coorconvert = '{java} -Xmx{mem} -cp /data/NYGC/Resources/SNPiR/ convertCoordinates - | samtools view -Sb - > {outputfile}.cord.bam'
    picard = '{java} -Xmx{mem} -jar {picardpath}{picardprog} TMP_DIR={tmpdir} VALIDATION_STRINGENCY=SILENT CREATE_INDEX=true '
    picard_addRG = ' RGID={fqid} RGLB={sample} RGPL=Illumina RGPU={barcode} RGSM={sample} RGCN=NYGC RGDS={fqinfo} SORT_ORDER=coordinate INPUT={outputfile}.cord.bam OUTPUT={outputfile}.sorted.bam '
    picard_merge = ' OUTPUT={sample}.sorted.bam ASSUME_SORTED=true USE_THREADING=true '
    picard_mdup = ' INPUT={sample}.sorted.bam OUTPUT={sample}.sorted.mdup.bam METRICS_FILE={sample}.sorted.mdup.metrics ASSUME_SORTED=true '


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
    callertb = {'mutect': [Template(' --input_file:normal $inputpath$normalsample --input_file:tumor $inputpath$tumorsample --out $outputpath$outbase.txt --coverage_file $outputpath$outbase.coverage.wig.txt --vcf $outputpath$outbase.vcf '), False], 'varscan': [Template('samtools mpileup -B -f $reference -q $minq $inputpath$normalsample > $outputpath$normalsample.pileup \n\nsamtools mpileup -B -f $reference -q $minq $inputpath$tumorsample > $outputpath$tumorsample.pileup\n\n$callcmd somatic $outputpath$normalsample.pileup $outputpath$tumorsample.pileup --output-snp $outputpath$outbase.snv --output-indel $outputpath$outbase.indel --output-vcf '), False], 'sniper': [Template(' $inputpath$tumorsample $inputpath$normalsample $outputpath$outbase.vcf '), True] }

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
            jobprefix = prefix + '_' + normalsample[sampidx] + '_' + tumorsample[sampidx]

            if addpara != "''": 
                paraset[addpara] = ''
                jobprefix = jobprefix + '_' + addpara.replace(' ','')

            tempset['outbase'] = jobprefix + '.' + callertag

            CMDs = []
            if 'toShell' in paraset.keys():
                CMDs.append( cmdGenerator.formatCmd( paraset.pop('toShell') ) )

            CMDs.append( cmdGenerator.checkPathOnNode(sampletmpoutpath) )

            if optionsbeforetemp:
                CMDs.append( cmdGenerator.formatCmd( callcmd, paraset, template.substitute(tempset) ) )
            else:
                CMDs.append( cmdGenerator.formatCmd( callcmd, template.substitute(tempset), paraset ) )

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


def DESeqPair(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmdset = configRobot.makeParasList(cmdset, ['meta', 'group1', 'group2'])
    cmd, mem, time, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'prefix'])
    group1, group2 = configRobot.popParas(cmdset, ['group1', 'group2'])
    template = open(cmdset.pop('template')).read()
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    cmdset['inputpath'] = inputpath
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    cmdset['outputpath'] = outputpath
    cmdset['prefix'] = prefix

    meta = configRobot.popParas(cmdset, ['meta'])
    if meta[0] == "''" and len(meta) == 1:
        cmdset['meta'] = 'c()'
    else:
        cmdset['meta'] = 'c(\'' + '\', \''.join(meta) + '\')'

    if cmdset['countfnprefix'] == "''": cmdset['countfnprefix'] = ''
    if cmdset['countfnsuffix'] == "''": cmdset['countfnsuffix'] = ''

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')
    for i in range(len(group1)):
        paraset = copy.deepcopy(cmdset)        
        paraset['gr1'] = group1[i]
        paraset['gr2'] = group2[i]
        jobfnprefix = prefix + '_' + group1[i] + '_' + group2[i]

        f = open('./%s.R'%(jobfnprefix), 'w')
        f.write(template%paraset)
        f.close()

        deseqcmd = cmdGenerator.formatCmd('Rscript', './%s.R'%(jobfnprefix))
        mvRscriptcmd = cmdGenerator.formatCmd('mv ./%s.R %s'%(jobfnprefix, outputpath))
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath))

        jobmanager.createJob(jobfnprefix, [setuppathcmd, deseqcmd, mvRscriptcmd, mvscriptcmd], outpath = outputpath, outfn = jobfnprefix)
    return jobmanager
    
    

def HTSeqCount(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath) 
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    gtf = configRobot.popParas(cmdset, ['GTF'])

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]

    for sample in samples:
        paraset = copy.deepcopy(cmdset)        
        jobfnprefix = prefix + '_' + sample
        if bam == '=sample':
            inputfile = sample
            outputfile = sample + '.count'
        else: 
            inputfile = sample + '/' + bam
            outputfile = sample + '/' + bam + '.count'
            cmdGenerator.checkPath(outputpath + sample, create=createpath)

        samcmd = 'samtools view -h %s | '%(inputpath+inputfile)
        htseq = 'python -m HTSeq.scripts.count -q '
        countcmd = cmdGenerator.formatCmd(samcmd, htseq, paraset, '-', gtf, ' > %s'%(outputpath + outputfile))
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, outputpath))

        jobmanager.createJob(jobfnprefix, [setuppathcmd, countcmd, mvscriptcmd], outpath = outputpath, outfn = jobfnprefix)
    return jobmanager
        

#mapping
#pipeline component: run tophat to map data
def tophat(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix'])
    paired = configRobot.popParas(cmdset, ['paired'])
    readext, outpath, genome = configRobot.popParas(cmdset, ['readext', 'outputpath', 'genome'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools\necho $BOWTIE2_INDEXES')

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]

    cmdGenerator.checkPath(outpath + '%s/'%prefix, create=createpath)

    for sample in samples:
        if paired == 'paired' or paired == 'yes':
            reads = map(lambda(i): inputpath + sample + '_%d'%i + readext, [1,2])
        elif paired == 'single' or paired == 'no':
            reads = inputpath + sample + readext
        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = outpath + '%s/'%prefix + sample
        jobfnprefix = prefix + '_' + sample
        
        tophatcmd = cmdGenerator.formatCmd(cmd, paraset, genome, reads)
        mvscriptcmd = cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, paraset['-o']))

        if int(paraset['-p']) > 1: #multiple threads per job
            sgeopt = ['-pe smp ' + paraset['-p']]
        else:
            sgeopt = []

        #need to create the output directory first, otherwise SGE complains cannot put the stdout in its path
        cmdGenerator.checkPath(paraset['-o'], create=createpath)
        jobmanager.createJob(jobfnprefix, [setuppathcmd, tophatcmd, mvscriptcmd], outpath = paraset['-o'], outfn = jobfnprefix, sgeopt=sgeopt)
    return jobmanager


def cufflinks(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    outputpath = cmdGenerator.checkPath(outputpath + '%s/'%prefix, create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    for sample in samples:
        jobname = prefix + '_' + sample
        CMD = []
        CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = outputpath + sample
        paraset = configRobot.validParas(paraset, availParas['cufflinks'])
        cmdGenerator.checkPath(paraset['-o'], create=createpath)
        CMD.append( cmdGenerator.formatCmd(cmd, paraset, inputpath+'%s/'%sample+bam) )
        
        CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['-o'])) )
        sgeopt = []
        if '-p' in paraset.keys():
            if int(paraset['-p']) > 1: #multi threads
                sgeopt = ['-pe smp ' + paraset['-p']]
        elif '--num-threads' in paraset.keys():
            if int(paraset['--num-threads']) > 1:
                sgeopt = ['-pe smp ' + paraset['-p']]
        jobmanager.createJob(jobname, CMD, outpath=paraset['-o'], outfn=jobname, sgeopt=sgeopt)
    return jobmanager


def cuffcompare(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, gtf)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset['-o'] = outputpath + paraset['-o']
    paraset = configRobot.validParas(paraset, availParas['cuffcompare'])
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )

    jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager



def cuffmerge(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = '"'
    for sample in samples:
        sampletext = sampletext + '%s%s/%s\\n'%(inputpath, sample, gtf)
    sampletext = sampletext + '"'

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
    paraset = copy.deepcopy(cmdset)
    paraset['-o'] = outputpath + paraset['-o']

    CMD.append( cmdGenerator.formatCmd('echo', sampletext, '>', paraset['-o'] + '.samples') )
    
    paraset = configRobot.validParas(paraset, availParas['cuffmerge'])
    cmdGenerator.checkPath(paraset['-o'], create=createpath)
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, paraset['-o'] + '.samples') )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['-o'])) )
    CMD.append( cmdGenerator.formatCmd('rm -f', paraset['-o'] + '.samples') )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['-o'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager

def cuffdiff(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, bam)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset = configRobot.validParas(paraset, availParas['cuffdiff'])
    cmdGenerator.checkPath(paraset['--output-dir'], create=createpath)
    CMD.append( cmdGenerator.formatCmd(cmd, paraset, gtf, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['--output-dir'])) )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['--output-dir'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager


def cuffdiff_v1(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, gtf, bam = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'sample', 'prefix', 'gtf', 'bam'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))

    if type(samples) != type([]) and type(samples) != type(()):
        samples = [samples]
    
    sampletext = ''
    for sample in samples:
        sampletext = sampletext + '%s%s/%s '%(inputpath, sample, bam)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    jobname = prefix
    CMD = []
    CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

    paraset = copy.deepcopy(cmdset)
    paraset = configRobot.validParas(paraset, availParas['cuffdiff'])
    cmdGenerator.checkPath(paraset['--output-dir'], create=createpath)
    CMD.append( cmdGenerator.formatCmd('/ifs/home/c2b2/dp_lab/bc2252/SeqTool/cufflinks_1/cuffdiff', paraset, gtf, sampletext) )
    CMD.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, paraset['--output-dir'])) )
    sgeopt = []
    if '-p' in paraset.keys():
        if int(paraset['-p']) > 1: #multi threads
            sgeopt = ['-pe smp ' + paraset['-p']]
    elif '--num-threads' in paraset.keys():
        if int(paraset['--num-threads']) > 1:
            sgeopt = ['-pe smp ' + paraset['-p']]
    jobmanager.createJob(jobname, CMD, outpath=paraset['--output-dir'], outfn=jobname, sgeopt=sgeopt, trackcmd=False)
    return jobmanager



def countmismatches(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    direct = ['forward', 'reverse']
    for sample in samples:
        for d in direct:
            jobname = prefix+'_'+sample+'_'+d
            CMD = []
            CMD.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
            CMD.append( cmdGenerator.formatCmd('python','countmismatches.py',inputpath+sample+'/'+bam, outputpath+sample+'.%s.countmis'%d, d) )
            CMD.append( cmdGenerator.formatCmd('mv', '%s%s'%(jobname, jobmanager.ext), outputpath) )
            jobmanager.createJob(jobname, CMD, outpath=outputpath, outfn=jobname, trackcmd=False)
    return jobmanager
    

def filetersingleton(cmdset, runmode='test'): 
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time, samples, prefix, bam = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix', 'bam'])
    if 'TMP_DIR' in cmdset.keys():
        TMP_DIR = cmdset.pop('TMP_DIR')
    else:
        TMP_DIR = ''

    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setuppathcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')


    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-1)
    samview = 'samtools view -b -h -F 8'
    reorder = 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'
    RG = 'AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT RGLB=dUTP RGPL=illumina RGPU=1'
    mdup = 'MarkDuplicates.jar'
    

    for sample in samples:
        jobfnprefix = prefix + '_' + sample
        paraset = copy.deepcopy(cmdset)
        if TMP_DIR != '': paraset['TMP_DIR'] = '%s'%TMP_DIR
        tmpbam = []

        CMDs = []
        CMDs.append(setuppathcmd)

        paraset['INPUT'] = '=%s/%s'%(inputpath+sample, bam)
        paraset['OUTPUT'] = '=%s/%s.reorder.bam'%(inputpath+sample, bam.replace('.bam',''))
        tmpbam.append(paraset['OUTPUT'].strip('='))
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+reorder, paraset) )
        
        paraset['INPUT'] = '=%s/%s.reorder.bam'%(inputpath+sample, bam.replace('.bam',''))
        paraset['OUTPUT'] = '=%s/%s.reorder.addRG.bam'%(inputpath+sample, bam.replace('.bam',''))
        paraset['RGSM'] = '=%s'%sample        
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+RG, paraset) )

        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = '%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam','.addRG'))
        CMDs.append( cmdGenerator.formatCmd(samview, paraset, inputpath+sample+'/'+bam.replace('.bam','.addRG.bam')) )

        
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam','.addRG'))
        paraset['OUTPUT'] = '%s/%s.mdup.bam'%(paraset['INPUT'].replace('.bam', ''))
        CMDs.append( cmdGenerator.formatCmd(javacmd, programpath+mdup, paraset) )

        paraset = copy.deepcopy(cmdset)
        CMDs.append( cmdGenerator.formatCmd('samtools index', bam.replace('.bam', '.reorder.addRG.filter.mdup.bam')) )


        CMDs.append( cmdGenerator.formatCmd('rm -f', tmpbam) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobfnprefix, jobmanager.ext, inputpath+sample)) )

        jobmanager.createJob(jobfnprefix, CMDs, outpath = inputpath+sample, outfn = jobfnprefix)
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


def preGATK(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    picardpath = cmdGenerator.checkPath(cmdset.pop('picardpath'))
    gatkpath = cmdGenerator.checkPath(cmdset.pop('gatkpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)
        

    samview = 'samtools view -b -h -F 264'
    reorder = picardpath + 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'
    RG = picardpath + 'AddOrReplaceReadGroups.jar VALIDATION_STRINGENCY=LENIENT RGLB=dUTP RGPL=illumina RGPU=1'
    mdupjar = picardpath + 'MarkDuplicates.jar'
    GATK = gatkpath + 'GenomeAnalysisTK.jar '
    createTg = '-T RealignerTargetCreator '
    realign = '-T IndelRealigner '

    idxcmd = 'samtools index'
    clearup = 'rm -f '


    for sample in samples:
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        jobname = prefix + '_' + sample

        #filter
        paraset = copy.deepcopy(cmdset)
        paraset['-o'] = '%s/%s.filter.bam'%(inputpath+sample, bam.replace('.bam',''))
        lastoutput = paraset['-o']
        del paraset['-R']
        del paraset['-filterMBQ']
        #paraset = configRobot.validParas(paraset, availParas['samtools'])
        CMDs.append( cmdGenerator.formatCmd(samview, paraset, inputpath+sample+'/'+bam ) )

        #reorder by chrm
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.reorder.bam'%(lastoutput.replace('.bam',''))
        paraset['REFERENCE'] = '=%s'%paraset['-R']
        paraset = configRobot.validParas(paraset, availParas['ReorderSam.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, reorder, paraset) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        lastoutput = paraset['OUTPUT'].strip('=')
        
        #add RG
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.addRG.bam'%(lastoutput.replace('.bam',''))
        paraset['RGSM'] = '=%s'%sample        
        paraset = configRobot.validParas(paraset, availParas['AddOrReplaceReadGroups.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, RG, paraset) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        lastoutput = paraset['OUTPUT'].strip('=')

        #mark duplicates
        paraset = copy.deepcopy(cmdset)
        paraset['INPUT'] = '=%s'%lastoutput
        paraset['OUTPUT'] = '=%s.mdup.bam'%(lastoutput.replace('.bam', ''))
        paraset['METRICS_FILE'] = '=%s/%s'%(inputpath + sample, prefix + '_mdupmetrics.txt')
        paraset = configRobot.validParas(paraset, availParas['MarkDuplicates.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, mdupjar, paraset) )
        CMDs.append( cmdGenerator.formatCmd(idxcmd, paraset['OUTPUT'].strip('=')) )
        lastoutput = paraset['OUTPUT'].strip('=')

        #create intervals
        paraset = copy.deepcopy(cmdset)
        paraset['-I'] = lastoutput
        paraset['-o'] = lastoutput.replace('.bam', '.intervals')
        CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+createTg, paraset) )

        #realign
        paraset['-targetIntervals'] = paraset['-o']
        paraset['-o'] = lastoutput.replace('.bam', '.realign.bam')
        CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+realign, paraset) )

        #clear up
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput) )
        CMDs.append( cmdGenerator.formatCmd(clearup, lastoutput.replace('.bam', '.intervals')) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, inputpath+sample)) )
        jobmanager.createJob(jobname, CMDs, outpath = inputpath+sample, outfn = jobname)
    
    return jobmanager


def picardReorderSam(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)        
    reorder = programpath + 'ReorderSam.jar VALIDATION_STRINGENCY=LENIENT'


    for sample in samples:
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )

        jobname = prefix + '_' + sample

        #reorder by chrm
        paraset = copy.deepcopy(cmdset)
        if bam == '=sample':
            inputfile = sample
        else:
            inputfile = sample + '/' + bam
        paraset['INPUT'] = '=%s'%(inputpath + inputfile)
        paraset['OUTPUT'] = '=%s.reorder.bam'%(outputpath + inputfile.replace('.bam',''))
        paraset = configRobot.validParas(paraset, availParas['ReorderSam.jar'])
        CMDs.append( cmdGenerator.formatCmd(javacmd, reorder, paraset) )
        
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
        jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
    
    return jobmanager


def GATK_genotyper(cmdset, runmode='test'):
    if runmode == 'test':
        createpath = False
    else:
        createpath = True
    
    cmd, mem, time, samples, prefix = configRobot.popParas(cmdset,['cmd', 'mem', 'time', 'sample', 'prefix'])
    
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    gatkpath = cmdGenerator.checkPath(cmdset.pop('gatkpath'))

    bam = configRobot.popParas(cmdset, 'bam')
    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)
        
    GATK = gatkpath + 'GenomeAnalysisTK.jar '
    genotyper = '-T UnifiedGenotyper '


    CMDs = []
    CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
    
    jobname = prefix

    paraset = copy.deepcopy(cmdset)
    paraset['-I'] = inputpath + samples[0] + '/' + bam
    for si in range(1,len(samples)):
        paraset['-I'] = paraset['-I'] + ' -I ' + inputpath + samples[si] + '/' + bam

    CMDs.append( cmdGenerator.formatCmd(javacmd, GATK+genotyper, paraset) )

    CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
    jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
    
    return jobmanager


def picardQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    metrics = {'CollectRnaSeqMetrics.jar': 'RnaSeq', 'CollectMultipleMetrics.jar': '', 'EstimateLibraryComplexity.jar': 'Lib', 'CollectGcBiasMetrics.jar': 'GC'}
    metrickeys = ['CollectRnaSeqMetrics.jar', 'CollectMultipleMetrics.jar', 'EstimateLibraryComplexity.jar', 'CollectGcBiasMetrics.jar']

    cmd, mem, time, bam, prefix = configRobot.popParas(cmdset, ['cmd', 'mem', 'time', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    if 'runmetric' in cmdset.keys():
        cmdset = configRobot.makeParasList(cmdset, ['runmetric'])
        runmetric = configRobot.popParas(cmdset, ['runmetric'])
        metrickeys = list(set(metrickeys).intersection(runmetric))

    samples = cmdset.pop('sample')

    pathsetup = 'source ./pathsetup'
    javacmd = 'java -Xmx%dg -jar'%(int(mem.replace('G',''))-2)

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    for sample in samples:
        jobname = prefix + '_' + sample.split('/')[-1]
        allcmds = []
            
        paraset = copy.deepcopy(cmdset)
        if bam == '=sample':
            paraset['INPUT'] = '=%s'%(inputpath+sample)
        else:
            paraset['INPUT'] = '=%s/%s'%(inputpath+sample, bam)
        if '/' in sample: sample = sample.split('/')[-1]

        paraset['TMP_DIR'] = paraset['TMP_DIR'] + prefix + '_' + sample + '/'
        cmdGenerator.checkPath(paraset['TMP_DIR'].strip('='), create=createpath)

        for metric in metrickeys:
            if 'MultipleMetrics' in metric:
                paraset['OUTPUT'] = '=%s'%(outputpath + sample + metrics[metric])
            else:
                paraset['OUTPUT'] = '=%s.txt'%(outputpath + sample + '.' + metrics[metric])
            paraset['CHART_OUTPUT'] = '%s'%(paraset['OUTPUT'].replace('.txt', '.pdf'))
            paraset['SUMMARY_OUTPUT'] = '%s'%(paraset['OUTPUT'].replace('.txt', '.summary.txt'))

            #filter out parameters that are not supported
            metricparaset = configRobot.validParas(paraset, availParas[metric])
            allcmds.append(cmdGenerator.formatCmd(pathsetup))
            allcmds.append(cmdGenerator.formatCmd(javacmd, programpath + metric, metricparaset))

        allcmds.append(cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)))
        allcmds.append(cmdGenerator.formatCmd('rm -Rf', paraset['TMP_DIR'].strip('=')))
        jobmanager.createJob(jobname, allcmds, outpath = outputpath, outfn = jobname)
    return jobmanager

    
    
def RNASeQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time = configRobot.popParas(cmdset, ['cmd', 'mem', 'time'])
    samples, bam, prefix = configRobot.popParas(cmdset, ['sample', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))
    if '-Djava.io.tmpdir' in cmdset.keys():
        javacmd = 'java ' + '-Djava.io.tmpdir' + cmdGenerator.checkPath(cmdset.pop('-Djava.io.tmpdir'))
    else:
        javacmd = 'java'
    javacmd = javacmd + ' -Xmx%dg -jar'%(int(mem.replace('G',''))-2)


    #need to generate -s, -o
    #generate a temperory file
    samplestr = '"Sample ID\\tBam File\\tNotes\\n'
    for sample in samples:
        samplestr = samplestr + '%s\\t%s\\t%s\\n'%(sample, inputpath + sample + '/' + bam, sample)
    samplestr = samplestr + '"'
    samplefile = '%s.samples'%prefix

    cmdset['-s'] = samplefile
    cmdset['-o'] = inputpath + 'RNA-SeQC/'

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))

    setupcmd = cmdGenerator.formatCmd('source ~/libraries/setup_seqtools')
    createsamplefile = cmdGenerator.formatCmd('echo', samplestr, '>', samplefile)
    removesamplefile = cmdGenerator.formatCmd('rm -f', samplefile)
    mvscriptcmd = cmdGenerator.formatCmd('mv %s%s %s'%(prefix, jobmanager.ext, cmdset['-o']))
    qccmd = cmdGenerator.formatCmd(javacmd, programpath+cmd, cmdset)

    cmdGenerator.checkPath(cmdset['-o'], create=createpath)
    jobmanager.createJob(prefix, [setupcmd, createsamplefile, qccmd, removesamplefile, mvscriptcmd], outfn = prefix, outpath = cmdset['-o'], trackcmd=False)
    return jobmanager


def RSeQC(cmdset, runmode='test'):
    global availParas
    if runmode == 'test':
        createpath = False
    else:
        createpath = True

    cmd, mem, time = configRobot.popParas(cmdset, ['cmd', 'mem', 'time'])
    samples, bam, prefix = configRobot.popParas(cmdset, ['sample', 'bam', 'prefix'])
    inputpath = cmdGenerator.checkPath(cmdset.pop('inputpath'))
    outputpath = cmdGenerator.checkPath(cmdset.pop('outputpath'), create=createpath)
    programpath = cmdGenerator.checkPath(cmdset.pop('programpath'))

    jobmanager = jobFactory.jobManager(mem=mem, time=time, overwrite=cmdset.pop('overwrite'))
    programs = ['inner_distance.py', 'junction_annotation.py', 'junction_saturation.py', 'read_GC.py', 'read_duplication.py']

    for sample in samples:
        jobname = prefix + '_' + sample
        CMDs = []
        CMDs.append( cmdGenerator.formatCmd('source ~/libraries/setup_seqtools') )
        for prog in programs:
            paraset = copy.deepcopy(cmdset)
            paraset['-i'] = inputpath + sample + '/' + bam
            paraset['-o'] = outputpath + sample + '.%s'%(prog.replace('.py', ''))
            paraset = configRobot.validParas(paraset, availParas[prog])
            if '-o' not in paraset.keys():
                paraset['>'] = outputpath + sample + '.%s'%(prog.replace('.py', ''))                            

            CMDs.append( cmdGenerator.formatCmd('python', programpath+prog, paraset) )
        CMDs.append( cmdGenerator.formatCmd('mv ./%s%s %s'%(jobname, jobmanager.ext, outputpath)) )
        jobmanager.createJob(jobname, CMDs, outpath = outputpath, outfn = jobname)
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

        jobmanager = ''


        if 'RNA-SeQC' in cmdsets[setname]['cmd']:
            jobmanager = RNASeQC(cmdsets[setname], runmode)
        else:
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






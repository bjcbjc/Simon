
%read RNA calls from bwa

%from somartic callers
samplepair = {'132540-10N_132540-1T', '132540-1T'; ...
    '138381-4N_138381-2T', '138381-2T'};
allsitefile = {'RNA/132540.allsites', 'RNA/138381.allsites'};
fnsuffix = { ...
    'mutect', 'bwa', '.mutect.vcf'; ...
    'mutect', 'star', '.mdup.MQ100.mutect.vcf'; ...
    'varscan', 'bwa', '.varscan.snv.vcf'; ...
    'varscan', 'star', '.mdup.varscan.snv.vcf'; ...
    'varscan', 'bwa', '.varscan.indel.vcf'; ...
    'varscan', 'star', '.mdup.varscan.indel.vcf'; ...
    'strelka', 'bwa', '.strelka.snv.vcf'; ...
    'strelka', 'star', '.mdup.strelka.snv.vcf'; ...
    'strelka', 'bwa', '.strelka.indel.vcf'; ...
    'strelka', 'star', '.mdup.strelka.indel.vcf'; ...
    'virmid', 'bwa', '.all.virmid.vcf'; ...
    'virmid', 'star', '.mdup.all.virmid.vcf'; ...
    'gatk', 'bwa', '.vcf'; ...
    'gatk', 'star', '.vcf'};

steps = {'rmhex', 'rmsk', 'rmintron', 'rmhom', 'rmblat', 'rmedit'};

res = cell(1, size(samplepair, 1));
for sampidx = 1:size(samplepair, 1)
    tic;
    f = fopen(allsitefile{sampidx}, 'r');
    t = textscan(f, '%s  %f');
    fclose(f);
    
    res{sampidx}.filterName = {'null', 'som', 'hex', 'rpt', 'intron', 'hom', 'blat', 'edit'};
    res{sampidx}.locidx = gloc2index( numericchrm(t{1}), t{2});
    N = length(res{sampidx}.locidx);
    
    fprintf('read all callers...\n');
    for vcfidx = 1:size(fnsuffix, 1)
        datadir = ['RNA/' fnsuffix{vcfidx,1}, '_', fnsuffix{vcfidx,2}, '/'];
        fdname = [ fnsuffix{vcfidx,1}, '_', fnsuffix{vcfidx,2}];
        if ~isempty(strfind( fnsuffix{vcfidx,3}, 'indel'))
            fdname = [fdname, '_indel'];
        end
        
        res{sampidx}.(fdname).validAfterFilter = false(length(res{sampidx}.locidx), 2+length(steps));
        res{sampidx}.(fdname).numReadRef = zeros(size(res{sampidx}.(fdname).validAfterFilter));
        res{sampidx}.(fdname).numReadAlt = zeros(size(res{sampidx}.(fdname).validAfterFilter));
                       
        if ~strcmp(fnsuffix{vcfidx,1}, 'gatk')
            vcffn = [datadir samplepair{sampidx, 1}, fnsuffix{vcfidx,3}];
            [locidx, filter, error] = readTumorSomaticCalls(...
                vcffn, fnsuffix{vcfidx,1});
            filter = strcmp(filter, 'PASS');
            if ~isempty(error)
                fprintf( '%s\n', error);
                break
            end
        else
            vcffn = [datadir samplepair{sampidx,2}, fnsuffix{vcfidx,3}];
            [vcfdata, status] = VCFFUNC.extract( ...
                vcffn, ...
                'returncol', [1,2,6], 'parseformatstr', '%s %f %f');                        
            if status == 0
                fdname = [ fnsuffix{vcfidx,1}, '_', fnsuffix{vcfidx,2}];
                locidx = gloc2index(numericchrm( vcfdata{1}), vcfdata{2});
                filter = vcfdata{3} >= 20;
            else
                fprintf('error extract gatk, %s, %s\n', samplepair{sampidx,2}, fnsuffix{vcfidx,2});
            end
        end
    
        res{sampidx}.(fdname).validAfterFilter( ismember(res{sampidx}.locidx, locidx), 1) = true;
        res{sampidx}.(fdname).validAfterFilter( ismember(res{sampidx}.locidx, locidx(filter)), 2) = true;
                        
        for stepidx = 1:length(steps)
            if strcmpi(steps{stepidx}, 'rmedit')
                fn = [vcffn, '.', strjoin(steps(1:stepidx), '.'), '.bed'];
            else
                fn = [vcffn, '.', strjoin(steps(1:stepidx), '.'), '.txt'];            
            end
            if strcmp(fnsuffix{vcfidx,1}, 'gatk')
                fn = strrep(fn, '.vcf', '');
            end
            filterres = readSNPiRTxtFile(fn);
            [~, idx] = ismember( filterres.locidx, res{sampidx}.locidx);
            res{sampidx}.(fdname).validAfterFilter(idx(idx~=0), stepidx+2) = true;
            res{sampidx}.(fdname).numReadRef(idx(idx~=0),stepidx+2) = filterres.numReadRef(idx~=0);
            res{sampidx}.(fdname).numReadAlt(idx(idx~=0),stepidx+2) = filterres.numReadAlt(idx~=0);            
        end
        res{sampidx}.(fdname).validAfterFilter = sparse(res{sampidx}.(fdname).validAfterFilter);
        res{sampidx}.(fdname).numReadRef = sparse(res{sampidx}.(fdname).numReadRef);
        res{sampidx}.(fdname).numReadAlt = sparse(res{sampidx}.(fdname).numReadAlt);
    end
    toc;
end

for sampidx = 1:size(samplepair,1)
    rnavar = res{sampidx};
    save(sprintf('data/rnavar.%s.withfilter.mat',samplepair{sampidx,1}), 'rnavar');
end

    
    






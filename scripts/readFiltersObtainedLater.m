
%get filters obtained later
samplepair = {'132540-1T--132540-10N', '138381-2T--138381-4N'};
matfnsample = {'132540-10N_132540-1T', '138381-4N_138381-2T'};
fnsuffix = { ...
    'strelka', '.indel.strelka.v0.4.7.vcf'; ...    
    'strelka', '.snv.strelka.v0.4.7.vcf'; ...
    'varscan', '.indel.varscan.v2.3.6.vcf'; ...
    'varscan', '.snv.varscan.somatic.v2.3.6.vcf'; ...
    'mutect', '.snv.mutect.somatic.v1.1.4.vcf'; ...
    'virmid', '.snv.virmid.v1.0.2.vcf'};

startstep = 5;
steps = {'rmhex', 'rmsk', 'rmintron', 'rmhom', 'rmblat', 'rmedit'};
aligner = {'bwa'}; %, 'star'}; % use RNA's coverage to filter calls
filterfndir = 'DNA/vcf/snpirfiltered/TMPSPLIT/SUBSPLIT/TMP/';

for sampidx = 2%1:length(samplepair)
    tic;
    newdnavar = loadStructData(sprintf('data/dnavar.%s.withfilter.mat', matfnsample{sampidx}));
    %read SNPiR filter results
    for vcfidx = [1, 2, 6]
        for aidx = 1:length(aligner)
            fdname = [ fnsuffix{vcfidx,1}, '_', aligner{aidx}];
            if ~isempty(strfind( fnsuffix{vcfidx,2}, 'indel'))
                fdname = [fdname, '_indel'];
            end

            newdnavar.(fdname).validAfterFilter = full(newdnavar.(fdname).validAfterFilter);
            newdnavar.(fdname).numReadRef = full(newdnavar.(fdname).numReadRef);
            newdnavar.(fdname).numReadAlt = full(newdnavar.(fdname).numReadAlt);


            for stepidx = startstep:length(steps)
                if strcmpi(steps{stepidx}, 'rmedit')
                    fn = [filterfndir, samplepair{sampidx}, fnsuffix{vcfidx,2}, '.', strjoin(steps(1:stepidx), '.'), '.bed'];
                else
                    fn = [filterfndir, samplepair{sampidx}, fnsuffix{vcfidx,2}, '.', strjoin(steps(1:stepidx), '.'), '.txt'];
                end
                if ~exist(fn, 'file')
                    fprintf('%s: skip after %s\n', fdname, steps{stepidx});
                    break;
                end
                filterres = readSNPiRTxtFile(fn);
                [~, idx] = ismember( filterres.locidx, newdnavar.locidx);
                newdnavar.(fdname).validAfterFilter(idx(idx~=0), stepidx+2) = true;
                newdnavar.(fdname).numReadRef(idx(idx~=0),stepidx+2) = filterres.numReadRef(idx~=0);
                newdnavar.(fdname).numReadAlt(idx(idx~=0),stepidx+2) = filterres.numReadAlt(idx~=0);
            end
            newdnavar.(fdname).validAfterFilter = sparse(newdnavar.(fdname).validAfterFilter);
            newdnavar.(fdname).numReadRef = sparse(newdnavar.(fdname).numReadRef);
            newdnavar.(fdname).numReadAlt = sparse(newdnavar.(fdname).numReadAlt);
        end
    end
    save(sprintf('data/newdnavar.%s.withfilter.mat', matfnsample{sampidx}), 'newdnavar');
    clear newdnavar
    toc;
end

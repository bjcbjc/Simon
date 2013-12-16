
%read DNA calls, tumor only

samplepair = {'132540-1T--132540-10N', '138381-2T--138381-4N'};
fnsuffix = { ...
    'strelka', '.indel.strelka.v0.4.7.vcf'; ...    
    'strelka', '.snv.strelka.v0.4.7.vcf'; ...
    'varscan', '.indel.varscan.v2.3.6.vcf'; ...
    'varscan', '.snv.varscan.somatic.v2.3.6.vcf'; ...
    'mutect', '.snv.mutect.somatic.v1.1.4.vcf'; ...
    'virmid', '.snv.virmid.v1.0.2.vcf'};

steps = {'rmhex', 'rmsk', 'rmintron', 'rmhom', 'rmblat', 'rmedit'};
aligner = {'bwa', 'star'}; % use RNA's coverage to filter calls
datadir = 'DNA/vcf/originals/'; 
covlimit = [1, 6, 11];

% res = cell(1, length(samplepair));
for sampidx = 1:length(samplepair)
    tic;
%     fprintf('read all callers...\n');
%     for vcfidx = 1:size(fnsuffix, 1)              
%         [locidx, filter, error] = readTumorSomaticCalls(...
%             [datadir samplepair{sampidx}, fnsuffix{vcfidx,2}], fnsuffix{vcfidx,1});
%         filter = strcmp(filter, 'PASS');
%         if ~isempty(error)
%             fprintf( '%s\n', error);
%             break
%         else
%             for aidx = 1:length(aligner) % duplicate, because later will be filtered by RNA's mapped reads
%                 fdname = [fnsuffix{vcfidx,1}, '_', aligner{aidx}];
%                 if ~isempty(strfind( fnsuffix{vcfidx,2}, 'indel'))
%                     fdname = [fdname, '_indel'];
%                 end
%                 res{sampidx}.(fdname).locidx = locidx;
%                 res{sampidx}.(fdname).filter = filter;
%             end
%         end
%     end
%     
%     if ~isempty(error) break; end
%     
%     fprintf('combines...\n');
%     fds = fieldnames(res{sampidx});
%     res{sampidx}.locidx = [];
%     for fdidx = 1:length(fds)
%         res{sampidx}.locidx = union(res{sampidx}.locidx, res{sampidx}.(fds{fdidx}).locidx);
%     end

%     res{sampidx}.covlimit = covlimit;
%     %read filtered by RNA coverage
%     for vcfidx = 1:size(fnsuffix, 1)
%         for aidx = 1:length(aligner)
%             beddir = ['DNA/vcf/filterByRNACov_' aligner{aidx}, '/'];
%             fdname = [ fnsuffix{vcfidx,1}, '_', aligner{aidx}];
%             if ~isempty(strfind( fnsuffix{vcfidx,2}, 'indel'))
%                 fdname = [fdname, '_indel'];
%             end
%             
%             res{sampidx}.(fdname).covfilter = false(length(res{sampidx}.locidx), length(res{sampidx}.covlimit));                        
%             for covidx = 1:length(res{sampidx}.covlimit)                
%                 fn = [beddir, samplepair{sampidx}, fnsuffix{vcfidx,2}, sprintf('.sorted.intersect.cov%d.bed',res{sampidx}.covlimit(covidx))];
%                 f = fopen(fn, 'r');
%                 data = textscan(f, '%s %*s %f %*[^\n]');
%                 fclose(f);
%                 res{sampidx}.(fdname).covfilter( ...
%                     ismember(res{sampidx}.locidx, gloc2index(numericchrm(data{1}), data{2})), covidx) = true;
%             end            
%         end
%     end
    
%     res{sampidx}.filterName = {'null', 'som', 'hex', 'rpt', 'intron', 'hom'};%, 'blat', 'edit'};    
    %read SNPiR filter results
    for vcfidx = 4 %1:size(fnsuffix, 1)
        for aidx = 1:length(aligner)            
            fdname = [ fnsuffix{vcfidx,1}, '_', aligner{aidx}];
            if ~isempty(strfind( fnsuffix{vcfidx,2}, 'indel'))
                fdname = [fdname, '_indel'];
            end
            
            res{sampidx}.(fdname).validAfterFilter = full(res{sampidx}.(fdname).validAfterFilter);
            res{sampidx}.(fdname).numReadRef = full(res{sampidx}.(fdname).numReadRef);
            res{sampidx}.(fdname).numReadAlt = full(res{sampidx}.(fdname).numReadAlt);
%             res{sampidx}.(fdname).validAfterFilter = false(length(res{sampidx}.locidx), 2+length(steps));
%             res{sampidx}.(fdname).validAfterFilter( ...
%                 ismember(res{sampidx}.locidx, res{sampidx}.(fdname).locidx), 1) = true;
%             res{sampidx}.(fdname).validAfterFilter( ...
%                 ismember(res{sampidx}.locidx, res{sampidx}.(fdname).locidx(res{sampidx}.(fdname).filter)), 2) = true;
%             res{sampidx}.(fdname) = rmfield(res{sampidx}.(fdname), {'locidx', 'filter'});             
%             res{sampidx}.(fdname).numReadRef = zeros(size(res{sampidx}.(fdname).validAfterFilter));
%             res{sampidx}.(fdname).numReadAlt = zeros(size(res{sampidx}.(fdname).validAfterFilter));
            
            for stepidx = 1:length(steps)
                if strcmpi(steps{stepidx}, 'rmedit')
                    fn = ['DNA/vcf/snpirfiltered/', samplepair{sampidx}, fnsuffix{vcfidx,2}, '.', strjoin(steps(1:stepidx), '.'), '.bed'];
                else
                    fn = ['DNA/vcf/snpirfiltered/', samplepair{sampidx}, fnsuffix{vcfidx,2}, '.', strjoin(steps(1:stepidx), '.'), '.txt'];
                end
                if ~exist(fn, 'file')
                    fprintf('%s: skip after %s\n', fdname, steps{stepidx});
                    break;
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
    end
    toc;
end


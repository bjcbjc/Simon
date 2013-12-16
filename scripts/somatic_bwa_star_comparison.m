
do.readRnaVar = true;
do.readDnaVar = true;
% do.readSNPiR = false;
% do.readSNPiRmat = true;
do.countByFilter = [0, true]; %if 1x2, the second used to indicate save figures
do.compareAligner = [0, 1];
do.compareSample = [0, 1];
do.compareSNPiR = [1, 1];
% do.NToverlapBar = [0, 1];
% do.PMoverlap = [0,0];
% do.DRoverlapBar = [0, 1];
% do.DRSomaticoverlap = [0, 1];
%do.DRoverlapPercent = [1, 1];
para.figdir = 'figures/rnavar/somatic_snpir/';
para.dnaMinCov = 11; %coverage is based on RNA data
para.rnaMinCov = 11;
%para.covFilterAligner = 'bwa';
% para.somaticAsDRdiff = false;
% para.figdir = [para.figdir, sprintf('DnaCov%d_RnaCov%d/',para.dnaMinCov, para.rnaMinCov)];

sample = {'Sample_132540-10N', 'Sample_132540-1T', ...
    'Sample_137064-1T', 'Sample_139508-1T', 'Sample_138381-4N', 'Sample_138381-2T'};
samplePair = {'Sample_132540-10N', 'Sample_132540-1T'; ...
    'Sample_132540-10N', 'Sample_137064-1T'; ...
    'Sample_132540-10N', 'Sample_139508-1T'; ...
    'Sample_138381-4N', 'Sample_138381-2T'};

if do.readDnaVar && ~exist('DnaVar', 'var')
    DnaVar.sample = {'132540', '138381'};
    DnaVar.res = cell(size(DnaVar.sample));
    for i = 1:length(DnaVar.sample)
        tmp = loadStructData(['data/dnavar.' DnaVar.sample{i} '.mat']);        
        DnaVar.res{i} = tmp.tumor;
        DnaVar.res{i}.covlimit = tmp.covlimit;        
    end
    clear tmp
end


if do.readRnaVar && ~exist('RnaVar', 'var')
    RnaVar.sample = {'132540', '138381'};
    RnaVar.res = cell(1, length(RnaVar.sample));
    for i = 1:length(RnaVar.sample)
        RnaVar.res{i} = loadStructData(['data/rnavar.' RnaVar.sample{i} '.mat']);
    end
end

% figure('visible', 'on');
% set(gcf, 'position', [0,0, 3000, 800]);
if do.countByFilter(1)
    for i = 1:length(RnaVar.sample)
        fds = fieldnames(RnaVar.res{i});
        fds(ismember(fds, {'locidx','filterName'})) = [];
        for j = 1:length(fds)
            if mod(j,2) == 1
                subplot(2, 6, ceil(j/2));                
            else
                subplot(2, 6, 6+j/2);
            end            
            nread = RnaVar.res{i}.(fds{j}).numReadRef(:,end) + RnaVar.res{i}.(fds{j}).numReadAlt(:,end);
            count = full(sum(RnaVar.res{i}.(fds{j}).validAfterFilter, 1));
            count(1, end+1:end+2) = [sum(RnaVar.res{i}.(fds{j}).validAfterFilter(:,end) & nread >= 6), ...
                sum(RnaVar.res{i}.(fds{j}).validAfterFilter(:,end) & nread >= 11) ];
            count(2,1:end-2) = sum( bsxfun(@and, RnaVar.res{i}.(fds{j}).validAfterFilter(:,2), ...
                RnaVar.res{i}.(fds{j}).validAfterFilter), 1);
            tmp = RnaVar.res{i}.(fds{j}).validAfterFilter(:,end) & ...
                RnaVar.res{i}.(fds{j}).validAfterFilter(:,2);
            count(2,end-1:end) = [ sum(tmp & nread>=6), sum(tmp & nread >= 11) ];
            count(1,:) = count(1,:);
            count(:, 2) = [];
            filtername = [RnaVar.res{i}.filterName, {'6+x', '11+x'}];
            filtername(2) = [];
            nfilter = length(filtername);
            [ax,h1,h2] = plotyy(1:nfilter, count(2,:), 1:nfilter, count(2,:)./count(1,:));
            set(h1, 'linestyle', '-', 'marker', 'x');
            set(h2, 'linestyle', '-', 'marker', 'x');
            set(get(ax(1),'ylabel'), 'string', '#somatic');
            set(get(ax(2),'ylabel'), 'string', '%somatic');
            set(ax(2),'xtick',[]);
            th = title(strrep(fds{j},'_', ' '));
            set(gca, 'xtick', 1:nfilter, 'xticklabel', filtername);
            xlabel('added filter');
            xlim(ax(1),[0.8, nfilter+0.2]);
            xlim(ax(2),[0.8, nfilter+0.2]);
            if j == 5
                thpos = get(th, 'position');
                thext = get(th, 'extent');
                text(thpos(1), thpos(2)+thext(4)*1.3, RnaVar.sample{i},'fontsize',12);
            end
        end        
        if do.countByFilter(2)            
            set(gcf, 'paperpositionmode', 'auto');
            saveas(gcf, [para.figdir, RnaVar.sample{i}, '.somatic_ratio_with_filter.png'], 'png');
        end
    end
end

if do.compareAligner(1)
    for i = 1:length(RnaVar.sample)
        fds = fieldnames(RnaVar.res{i});
        fds(ismember(fds, {'locidx','filterName'})) = [];
        fds(cellfun(@isempty, strfind(fds, 'bwa'))) = [];
        for j = 1:length(fds)            
            if mod(j,2) == 1
                subplot(2, 6, ceil(j/2));                
            else
                subplot(2, 6, 6+j/2);
            end  
            bwa = fds{j};
            star = strrep(fds{j}, '_bwa', '_star');
            filtername = [RnaVar.res{i}.filterName, {'6+x', '11+x'}];
            filtername(2) = [];
            nfilter = length(RnaVar.res{i}.filterName);            
            nread_bwa = RnaVar.res{i}.(bwa).numReadRef(:,end) + RnaVar.res{i}.(bwa).numReadAlt(:,end);
            nread_star = RnaVar.res{i}.(star).numReadRef(:,end) + RnaVar.res{i}.(star).numReadAlt(:,end);
            filteridx = setdiff(1:nfilter,2);
            covlimit = [6, 11];
            for snvtype = 1:2    
                bardata = zeros(nfilter+1, 3);
                if snvtype == 1
                    typefilter_A = true(size(RnaVar.res{i}.locidx));
                    typefilter_B = typefilter_A;
                else
                    typefilter_A = RnaVar.res{i}.(bwa).validAfterFilter(:,2);
                    typefilter_B = RnaVar.res{i}.(star).validAfterFilter(:,2);
                end
                for ftidx = 1:nfilter+1
                    if ftidx < nfilter
                        A = RnaVar.res{i}.locidx( typefilter_A & ...
                            RnaVar.res{i}.(bwa).validAfterFilter(:,filteridx(ftidx)));
                        B = RnaVar.res{i}.locidx( typefilter_B & ...
                            RnaVar.res{i}.(star).validAfterFilter(:,filteridx(ftidx)));
                    else
                        A = RnaVar.res{i}.locidx( typefilter_A & ...
                            RnaVar.res{i}.(bwa).validAfterFilter(:,filteridx(end)) ...
                            & nread_bwa >= covlimit(ftidx-nfilter+1) );
                        B = RnaVar.res{i}.locidx( typefilter_B & ...
                            RnaVar.res{i}.(star).validAfterFilter(:,filteridx(end)) ...
                            & nread_star >= covlimit(ftidx-nfilter+1) );
                    end
                    bardata(ftidx,:) = [sum(~ismember(A,B)), sum(ismember(A,B)), sum(~ismember(B,A))];
                end
                
                subplot(2,6,(snvtype-1)*6+j);
                bar(bardata, 'stacked');            
                if j == 1 && snvtype == 1
                    legend({'bwa', 'overlap', 'star'});
                end
                if snvtype == 1
                    ylabel('# var');
                else
                    ylabel('#somatic');
                end
                th = title(strrep(fds{j},'_', ' '));
                set(gca, 'xtick', 1:nfilter+1, 'xticklabel', filtername);                
                xlim([0.5 nfilter+1.5]);
                ylim([0, sum(bardata(1,:))])
                if j == 3 && snvtype == 1
                    thpos = get(th, 'position');
                    thext = get(th, 'extent');
                    text(thpos(1), thpos(2)+thext(4)*1.3, RnaVar.sample{i},'fontsize',12);
                end
            end
        end        
        if do.compareAligner(2)            
            set(gcf, 'paperpositionmode', 'auto');
            saveas(gcf, [para.figdir, RnaVar.sample{i}, '.aligner_overlap.png'], 'png');
        end
    end
end

% if do.DRoverlapPercent(1)                    
%     caller = {'mutect', 'varscan', 'strelka', 'virmid'};    
%     aligner = {'bwa', 'star'};
%     for i = 1:length(RnaVar.sample)
%         covidx = find(DnaVar.res{i}.covlimit==para.dnaMinCov)+1;
%         for aidx = 1:length(aligner)
%             for ci = 1:length(caller)
%                 tmp.(caller{ci}) = DnaVar.res{dnaSampIdx(i)}.locidx(...
%                     DnaVar.res{dnaSampIdx(i)}.(para.covFilterAligner).(caller{ci})(:,covidx) );
%             end
%         
%         tmp.all = intersect(tmp.mutect, intersect(tmp.sniper, tmp.varscan));
%         tmp.union = union(tmp.mutect, union(tmp.sniper, tmp.varscan));
%         tmp.any2 = union(intersect(tmp.mutect, tmp.sniper), ...
%             union(intersect(tmp.mutect, tmp.varscan), intersect(tmp.sniper, tmp.varscan)));
%         
%         nfilter = length(SNPiRres.res{snpirSampIdx(i)}.filterName)+1;
%         fds = fieldnames(tmp);
%         bardata = zeros(3,length(fds),nfilter);
%         for fti = 1:nfilter
%             if fti == 1
%                 ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,1) + ...
%                     SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,1);
%                 slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx(ntotalread >= para.rnaMinCov);
%             else
%                 ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,fti-1) + ...
%                     SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,fti-1);
%                 slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx( ...
%                     SNPiRres.res{snpirSampIdx(i)}.validAfterFilter(:,fti-1) & ...
%                     ntotalread >= para.rnaMinCov);
%             end
%             for fdi = 1:length(fds)                
%                 dlocidx = tmp.(fds{fdi});
%                 bardata(:,fdi,fti) = [sum(~ismember(dlocidx,slocidx)), ...
%                     sum(ismember(dlocidx, slocidx)), sum(~ismember(slocidx, dlocidx))];
%             end
%         end
%         for fdi = 1:length(fds)            
%             subplot(2,3,fdi);
%             pdata = squeeze(bardata(:,fdi,:))'; %nfilter x 3
%             pdata = [pdata(:,2)./sum(pdata(:,1:2),2), pdata(:,2)./sum(pdata(:,2:3),2)];
%             [ax,h1,h2] = plotyy(1:nfilter, pdata(:,1), 1:nfilter, pdata(:,2));
%             set(h1, 'linestyle', '-', 'marker', 'x');
%             set(h2, 'linestyle', '-', 'marker', 'x');
%             set(get(ax(1),'ylabel'), 'string', 'in DNA');
%             set(get(ax(2),'ylabel'), 'string', 'in RNA-snpir');
%             set(ax(2),'xtick',[]);
%             th = title(fds{fdi});
%             xticklabel = strrep(SNPiRres.res{snpirSampIdx(i)}.filterName, 'rm', '');
%             xticklabel = strrep(xticklabel, 'sk', 'rpt');
%             set(gca, 'xtick', 1:nfilter, 'xticklabel', [{'none'}, xticklabel]);
%             xlabel('added filter');
% 
%             if fdi == 2
%                 thpos = get(th, 'position');
%                 thext = get(th, 'extent');
%                 text(thpos(1), thpos(2)+thext(4)*1.3, strrep(overlapsample{i}, 'Sample_', ''),'fontsize',12);
%             end
%         end
%         if do.DRoverlapPercent(2)
%             saveas(fh, [para.figdir, strrep(overlapsample{i}, 'Sample_', ''), '.DNA_RNA_ratio.png'], 'png');
%         end
%     end
% end

if do.compareSample(1)
    aligner = 'star';
    caller = {'mutect', 'varscan', 'strelka'};%, 'virmid'};
    filteridx = setdiff(1:length(RnaVar.res{1}.filterName), 2);
    nfilter = length(filteridx);
    filterName = [RnaVar.res{1}.filterName(filteridx), {'6+x', '11+x'}];
    for cidx = 1:length(caller)
        for ftidx = 1:length(filterName)
            bardata = zeros(2, 2);
            if ftidx < nfilter
                useidx = filteridx(ftidx);
            else
                useidx = filteridx(end);
            end
            for sampi = 1:2
                bardata(sampi, 2) = sum( ...
                    RnaVar.res{sampi}.([caller{cidx} '_' aligner]).validAfterFilter(:, useidx) ...
                    & RnaVar.res{sampi}.([caller{cidx} '_' aligner]).validAfterFilter(:, 2));
                bardata(sampi, 1) = sum( ...
                    RnaVar.res{sampi}.([caller{cidx} '_' aligner]).validAfterFilter(:, useidx)) ...
                    - bardata(sampi,2);
            end
            subplot(3,3,ftidx);
            barh(bardata, 'stacked', 'barwidth', 0.5);
            if ftidx == 1
                legend({'non-som', 'somatic'}, 'location', 'SE');
            end
            set(gca, 'ytick', 1:2, 'yticklabel', RnaVar.sample);
            ylim([0.5 2.5]);
            if cidx ~= 3
                xlim([ min(bardata(:,1))*0.9, max(sum(bardata,2))]);
            end
            xlabel('# var');
            th = title(filterName{ftidx});
            if ftidx == 2
                thpos = get(th, 'position');
                thext = get(th, 'extent');
                text(thpos(1), thpos(2)+thext(4)*1.3, caller{cidx},'fontsize',12);
            end
        end
        if do.compareSample(2)
            saveas(gcf, [para.figdir, caller{cidx}, '_', aligner, '_sample_varcount.png'], 'png');
        end
    end
end

if do.compareSNPiR(1)
    [overlapsample, dnaSampIdx, snpirSampIdx] = intersect(DnaVar.sample, SNPiRres.sample);
    rmi = ~cellfun(@isempty, strfind(overlapsample, 'N'));
    overlapsample(rmi) = [];
    dnaSampIdx(rmi) = [];
    snpirSampIdx(rmi) = [];
    nsample = length(overlapsample);    
    caller = {'mutect', 'varscan', 'strelka', 'virmid'};
    for i = 1:nsample        
        covidx = find(DnaVar.res{dnaSampIdx(i)}.covlimit==para.dnaMinCov)+1;        
        for ci = 1:3 %get somatic only
            tmp.(caller{ci}) = DnaVar.res{dnaSampIdx(i)}.locidx(...
                DnaVar.res{dnaSampIdx(i)}.(para.covFilterAligner).(caller{ci})(:,covidx) & ...
                DnaVar.res{dnaSampIdx(i)}.([caller{ci} '_somatic']));
        end
        
        tmp.all = intersect(tmp.mutect, intersect(tmp.sniper, tmp.varscan));
        tmp.union = union(tmp.mutect, union(tmp.sniper, tmp.varscan));
        tmp.any2 = union(intersect(tmp.mutect, tmp.sniper), ...
            union(intersect(tmp.mutect, tmp.varscan), intersect(tmp.sniper, tmp.varscan)));
        
        %find normal sample in SNPiR
        normsampidx = find(strcmp(samplePair(strcmp(samplePair(:,2), overlapsample{i}),1), SNPiRres.sample));
        
        nfilter = length(SNPiRres.res{snpirSampIdx(i)}.filterName)+1;
        fds = fieldnames(tmp);
        bardata = zeros(3,length(fds),nfilter);
        for fti = 1:nfilter
            for fdi = 1:length(fds)
                if fti == 1
                    ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,1) + ...
                        SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,1);
                    slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx(ntotalread >= para.rnaMinCov);
                    if para.somaticAsDRdiff 
                        slocidx = setdiff(slocidx, ...
                            SNPiRres.res{normsampidx}.locidx);
                    end
                else
                    ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,fti-1) + ...
                        SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,fti-1);
                    slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx( ...
                        SNPiRres.res{snpirSampIdx(i)}.validAfterFilter(:,fti-1) & ...
                        ntotalread >= para.rnaMinCov);
                    if para.somaticAsDRdiff 
                        slocidx = setdiff(slocidx, ...
                            SNPiRres.res{normsampidx}.locidx( ...
                            SNPiRres.res{normsampidx}.validAfterFilter(:,fti-1)));
                    end
                end
                dlocidx = tmp.(fds{fdi});
                bardata(:,fdi,fti) = [sum(~ismember(dlocidx,slocidx)), ...
                    sum(ismember(dlocidx, slocidx)), sum(~ismember(slocidx, dlocidx))];
            end
        end
        clf
        for fdi = 1:length(fds)            
            subplot(2,3,fdi);
            pdata = squeeze(bardata(:,fdi,:))'; %nfilter x 3
            pdata = [pdata(:,2)./sum(pdata(:,1:2),2), pdata(:,2)./sum(pdata(:,2:3),2)];
            [ax,h1,h2] = plotyy(1:nfilter, pdata(:,1), 1:nfilter, pdata(:,2));
            set(h1, 'linestyle', '-', 'marker', 'x');
            set(h2, 'linestyle', '-', 'marker', 'x');
            set(get(ax(1),'ylabel'), 'string', 'in DNA');
            set(get(ax(2),'ylabel'), 'string', 'in RNA-snpir');
            set(ax(2),'xtick',[]);
%             set(ax(1), 'ylim', [0, max(pdata(:,1))+0.1]);
%             set(ax(2), 'ylim', [0, max(pdata(:,2))+0.1]);
%             bar(squeeze(bardata(:,fdi,:))', 'stacked')
            th = title(fds{fdi});
%             if fdi == 4
%                 legend({'DNA', 'overlap', 'RNA'});
%             end
            xticklabel = strrep(SNPiRres.res{snpirSampIdx(i)}.filterName, 'rm', '');
            xticklabel = strrep(xticklabel, 'sk', 'rpt');
            set(gca, 'xtick', 1:nfilter, 'xticklabel', [{'none'}, xticklabel]);
            xlabel('added filter');
%            ylabel('variant count');
%             xlim([0.5, nfilter+0.5])
%             ylim([0 sum(bardata(:,fdi,1))])            
            if fdi == 2
                thpos = get(th, 'position');
                thext = get(th, 'extent');
                text(thpos(1), thpos(2)+thext(4)*1.3, strrep(overlapsample{i}, 'Sample_', ''),'fontsize',12);
            end
        end
%         suptitle(strrep(overlapsample{i}, 'Sample_', ''));
        if do.DRSomaticoverlap(2)
            saveas(fh, [para.figdir, strrep(overlapsample{i}, 'Sample_', ''), '.DNA_RNA_somatic_overlap.png'], 'png');
        end
    end
    if do.compareSNPiR(2)
    end
end


do.readDnaVar = true;
do.readSNPiR = false;
do.readSNPiRmat = true;
do.countByFilter = [0, true]; %if 1x2, the second used to indicate save figures
do.NToverlapBar = [0, 1];
do.PMoverlap = [0,0];
do.DRoverlapBar = [0, 1];
do.DRSomaticoverlap = [0, 1];
do.DRoverlapPercent = [1, 1];
para.figdir = 'figures/rnavar/snpir/';
para.dnaMinCov = 11; %coverage is based on RNA data
para.rnaMinCov = 11;
para.covFilterAligner = 'bwa';
para.somaticAsDRdiff = false;
para.figdir = [para.figdir, sprintf('DnaCov%d_RnaCov%d/',para.dnaMinCov, para.rnaMinCov)];

sample = {'Sample_132540-10N', 'Sample_132540-1T', ...
    'Sample_137064-1T', 'Sample_139508-1T', 'Sample_138381-4N', 'Sample_138381-2T'};
samplePair = {'Sample_132540-10N', 'Sample_132540-1T'; ...
    'Sample_132540-10N', 'Sample_137064-1T'; ...
    'Sample_132540-10N', 'Sample_139508-1T'; ...
    'Sample_138381-4N', 'Sample_138381-2T'};

if do.readDnaVar
    fns = listfilename('data/dnavar*.mat');
    DnaVar.res = cell(length(fns)*2, 1);
    DnaVar.sample = cell(1, length(DnaVar.res));
    for i = 1:length(fns)
        tmp = loadStructData(['data/' fns{i}]);
        DnaVar.sample{i*2-1} = ['Sample_', tmp.sample{1}];
        DnaVar.sample{i*2} = ['Sample_', tmp.sample{2}];
        DnaVar.res{i*2-1} = tmp.normal;
        DnaVar.res{i*2-1}.covlimit = tmp.covlimit;
        DnaVar.res{i*2} = tmp.tumor;
        DnaVar.res{i*2}.covlimit = tmp.covlimit;
    end
    clear tmp
end

if do.readSNPiR
    SNPiRres.sample = sample;
    SNPiRres.res = cell(length(sample),1);
    for i = 1:length(sample)
        snpirRes = readSNPiRResult('RNA/snpir/vcf/', sample{i});
        save(sprintf('data/snpirres.%s.mat', sample{i}), 'snpirRes');
        SNPiRres.res{i} = snpirRes;
    end
end

if do.readSNPiRmat && ~exist('SNPiRres', 'var')
    SNPiRres.sample = sample;
    SNPiRres.res = cell(length(sample),1);
    for i = 1:length(sample)
        SNPiRres.res{i} = loadStructData(sprintf('data/snpirres.%s.mat', sample{i}));
    end
end

if do.countByFilter(1)
    if ~exist('fh', 'var')
        fh = figure('position', [0, 0, 1500, 600]);
        set(fh, 'paperpositionmode', 'auto');
    end
    nsample = length(sample);    
    for i = 1:nsample
        subplot(2,3, i);
        nfilter = length(SNPiRres.res{i}.filterName)+1;
        ntotalread = SNPiRres.res{i}.numReadRef(:,end) + SNPiRres.res{i}.numReadAlt(:,end);
        bar(1:nfilter+2, [length(SNPiRres.res{i}.locidx), full(sum(SNPiRres.res{i}.validAfterFilter,1)), ...
            full(sum(SNPiRres.res{i}.validAfterFilter(:,end) & ntotalread > 5)), ...
            full(sum(SNPiRres.res{i}.validAfterFilter(:,end) & ntotalread > 10))] );
        xticklabel = strrep(SNPiRres.res{i}.filterName, 'rm', '');
        xticklabel = strrep(xticklabel, 'sk', 'rpt');
        set(gca, 'xtick', 1:nfilter+2, 'xticklabel', [{'none'}, xticklabel, {'6+x', '11+x'}]);
        xlabel('added filter');
        ylabel('variant count');
        xlim([0.5, nfilter+2.5])
        ylim([0 length(SNPiRres.res{i}.locidx)])
        title(strrep(sample{i}, 'Sample_', ''));
    end
    if do.countByFilter(2)
        saveas(fh, [para.figdir 'filter_counts.png'], 'png');
    end
end

if do.NToverlapBar(1)
    if ~exist('fh', 'var')
        fh = figure('position', [0, 0, 800, 600]);
        set(fh, 'paperpositionmode', 'auto');        
    else        
        if ishghandle(fh)
            oldfigloc = get(fh, 'position');
            set(fh, 'position',  + [oldfigloc(1:2),800,600]);
        else
            fh = figure('position', [0, 0, 800, 600]);
            set(fh, 'paperpositionmode', 'auto');
        end
    end
    npair = size(samplePair, 1);
    bardata = zeros(nfilter, 3, npair);
    for i = 1:npair
        subplot(2,2,i);
        nfilter = length(SNPiRres.res{i}.filterName)+1;        
        [~, sampidx] = ismember(samplePair(i,:), SNPiRres.sample);        
        bardata(1, 1, i) = sum(~ismember(SNPiRres.res{sampidx(1)}.locidx, SNPiRres.res{sampidx(2)}.locidx));
        bardata(1, 2, i) = sum(ismember(SNPiRres.res{sampidx(1)}.locidx, SNPiRres.res{sampidx(2)}.locidx));
        bardata(1, 3, i) = sum(~ismember(SNPiRres.res{sampidx(2)}.locidx, SNPiRres.res{sampidx(1)}.locidx));
        for j = 2:nfilter
            a = SNPiRres.res{sampidx(1)}.locidx( SNPiRres.res{sampidx(1)}.validAfterFilter(:, j-1));
            b = SNPiRres.res{sampidx(2)}.locidx( SNPiRres.res{sampidx(2)}.validAfterFilter(:, j-1));
            bardata(j, :, i) = [ sum(~ismember(a,b)), sum(ismember(a,b)), sum(~ismember(b,a)) ];
        end
        bar(1:nfilter, bardata(:,:,i), 'stacked' );
        legend({'normal', 'overlap', 'tumor'});
        xticklabel = strrep(SNPiRres.res{sampidx(1)}.filterName, 'rm', '');
        xticklabel = strrep(xticklabel, 'sk', 'rpt');
        set(gca, 'xtick', 1:nfilter, 'xticklabel', [{'none'}, xticklabel]);
        xlabel('added filter');
        ylabel('variant count');
        xlim([0.5, nfilter+0.5])
        ylim([0 sum(bardata(1,:,i))])
        title(strjoin(strrep(sample(sampidx), 'Sample_',''), ' vs '));
    end
    if do.NToverlapBar(2)
        saveas(fh, [para.figdir 'normal_tumor_overlap_counts.png'], 'png');
    end
end

if do.PMoverlap(1)
end

if do.DRoverlapBar(1)    
    if ~exist('DnaVar', 'var')
        fns = listfilename('data/dnavar*.mat');
        DnaVar.res = cell(length(fns)*2, 1);
        DnaVar.sample = cell(1, length(DnaVar.res));
        for i = 1:length(fns)
            tmp = loadStructData(['data/' fns{i}]);
            DnaVar.sample{i*2-1} = ['Sample_', tmp.sample{1}];            
            DnaVar.sample{i*2} = ['Sample_', tmp.sample{2}];
            DnaVar.res{i*2-1} = tmp.normal;
            DnaVar.res{i*2-1}.covlimit = tmp.covlimit;
            DnaVar.res{i*2} = tmp.tumor;
            DnaVar.res{i*2}.covlimit = tmp.covlimit;
        end
        clear tmp
    end
    
    if ~exist('fh', 'var')
        fh = figure('position', [0, 0, 1200, 600]);
        set(fh, 'paperpositionmode', 'auto');        
    else        
        if ishghandle(fh)
            oldfigloc = get(fh, 'position');
            set(fh, 'position',  + [oldfigloc(1:2),1200,600]);
        else
            fh = figure('position', [0, 0, 1200, 600]);
            set(fh, 'paperpositionmode', 'auto');
        end
    end
    [overlapsample, dnaSampIdx, snpirSampIdx] = intersect(DnaVar.sample, SNPiRres.sample);
    nsample = length(overlapsample);    
    caller = {'mutect', 'sniper', 'varscan'};    
    for i = 1:nsample        
        covidx = find(DnaVar.res{dnaSampIdx(i)}.covlimit==para.dnaMinCov)+1;
        for ci = 1:3
            tmp.(caller{ci}) = DnaVar.res{dnaSampIdx(i)}.locidx(...
                DnaVar.res{dnaSampIdx(i)}.(para.covFilterAligner).(caller{ci})(:,covidx) );
        end        
        
        tmp.all = intersect(tmp.mutect, intersect(tmp.sniper, tmp.varscan));
        tmp.union = union(tmp.mutect, union(tmp.sniper, tmp.varscan));
        tmp.any2 = union(intersect(tmp.mutect, tmp.sniper), ...
            union(intersect(tmp.mutect, tmp.varscan), intersect(tmp.sniper, tmp.varscan)));
        
        nfilter = length(SNPiRres.res{snpirSampIdx(i)}.filterName)+1;
        fds = fieldnames(tmp);
        bardata = zeros(3,length(fds),nfilter);
        for fti = 1:nfilter
            if fti == 1
                ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,1) + ...
                    SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,1);
                slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx(ntotalread >= para.rnaMinCov);
            else
                ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,fti-1) + ...
                    SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,fti-1);
                slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx( ...
                    SNPiRres.res{snpirSampIdx(i)}.validAfterFilter(:,fti-1) & ...
                    ntotalread >= para.rnaMinCov);
            end
            for fdi = 1:length(fds)                
                dlocidx = tmp.(fds{fdi});
                bardata(:,fdi,fti) = [sum(~ismember(dlocidx,slocidx)), ...
                    sum(ismember(dlocidx, slocidx)), sum(~ismember(slocidx, dlocidx))];
            end
        end
        for fdi = 1:length(fds)
            subplot(2,3,fdi);
            bar(squeeze(bardata(:,fdi,:))', 'stacked')
            title(fds{fdi});
            if fdi == 4
                legend({'DNA', 'overlap', 'RNA'});
            end
            xticklabel = strrep(SNPiRres.res{snpirSampIdx(i)}.filterName, 'rm', '');
            xticklabel = strrep(xticklabel, 'sk', 'rpt');
            set(gca, 'xtick', 1:nfilter, 'xticklabel', [{'none'}, xticklabel]);
            xlabel('added filter');
            ylabel('variant count');
            xlim([0.5, nfilter+0.5])
            ylim([0 sum(bardata(:,fdi,1))])            
        end
        suptitle(strrep(overlapsample{i}, 'Sample_', ''));
        if do.DRoverlapBar(2)
            saveas(fh, [para.figdir, strrep(overlapsample{i}, 'Sample_', ''), '.DNA_RNA_overlap.png'], 'png');
        end
    end
end

if do.DRSomaticoverlap(1)    
    if ~exist('DnaVar', 'var')
        fns = listfilename('data/dnavar*.mat');
        DnaVar.res = cell(length(fns)*2, 1);
        DnaVar.sample = cell(1, length(DnaVar.res));
        for i = 1:length(fns)
            tmp = loadStructData(['data/' fns{i}]);
            DnaVar.sample{i*2-1} = ['Sample_', tmp.sample{1}];            
            DnaVar.sample{i*2} = ['Sample_', tmp.sample{2}];
            DnaVar.res{i*2-1} = tmp.normal;
            DnaVar.res{i*2-1}.covlimit = tmp.covlimit;
            DnaVar.res{i*2} = tmp.tumor;
            DnaVar.res{i*2}.covlimit = tmp.covlimit;
        end
        clear tmp
    end
    
    if ~exist('fh', 'var')
        fh = figure('position', [0, 0, 1500, 600]);
        set(fh, 'paperpositionmode', 'auto');        
    else        
        if ishghandle(fh)
            oldfigloc = get(fh, 'position');
            set(fh, 'position',  + [oldfigloc(1:2),1500,600]);
            set(fh, 'paperpositionmode', 'auto');
        else
            fh = figure('position', [0, 0, 1500, 600]);
            set(fh, 'paperpositionmode', 'auto');
        end
    end
    
    [overlapsample, dnaSampIdx, snpirSampIdx] = intersect(DnaVar.sample, SNPiRres.sample);
    rmi = ~cellfun(@isempty, strfind(overlapsample, 'N'));
    overlapsample(rmi) = [];
    dnaSampIdx(rmi) = [];
    snpirSampIdx(rmi) = [];
    nsample = length(overlapsample);    
    caller = {'mutect', 'sniper', 'varscan'};
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
end

if do.DRoverlapPercent(1)    
    if ~exist('DnaVar', 'var')
        fns = listfilename('data/dnavar*.mat');
        DnaVar.res = cell(length(fns)*2, 1);
        DnaVar.sample = cell(1, length(DnaVar.res));
        for i = 1:length(fns)
            tmp = loadStructData(['data/' fns{i}]);
            DnaVar.sample{i*2-1} = ['Sample_', tmp.sample{1}];            
            DnaVar.sample{i*2} = ['Sample_', tmp.sample{2}];
            DnaVar.res{i*2-1} = tmp.normal;
            DnaVar.res{i*2-1}.covlimit = tmp.covlimit;
            DnaVar.res{i*2} = tmp.tumor;
            DnaVar.res{i*2}.covlimit = tmp.covlimit;
        end
        clear tmp
    end
    
    if ~exist('fh', 'var')
        fh = figure('position', [0, 0, 1500, 600]);
        set(fh, 'paperpositionmode', 'auto');        
    else        
        if ishghandle(fh)
            oldfigloc = get(fh, 'position');
            set(fh, 'position',  + [oldfigloc(1:2),1500,600]);
        else
            fh = figure('position', [0, 0, 1500, 600]);
            set(fh, 'paperpositionmode', 'auto');
        end
    end
    [overlapsample, dnaSampIdx, snpirSampIdx] = intersect(DnaVar.sample, SNPiRres.sample);
    nsample = length(overlapsample);    
    caller = {'mutect', 'sniper', 'varscan'};    
    for i = 1:nsample        
        covidx = find(DnaVar.res{dnaSampIdx(i)}.covlimit==para.dnaMinCov)+1;
        for ci = 1:3
            tmp.(caller{ci}) = DnaVar.res{dnaSampIdx(i)}.locidx(...
                DnaVar.res{dnaSampIdx(i)}.(para.covFilterAligner).(caller{ci})(:,covidx) );
        end        
        
        tmp.all = intersect(tmp.mutect, intersect(tmp.sniper, tmp.varscan));
        tmp.union = union(tmp.mutect, union(tmp.sniper, tmp.varscan));
        tmp.any2 = union(intersect(tmp.mutect, tmp.sniper), ...
            union(intersect(tmp.mutect, tmp.varscan), intersect(tmp.sniper, tmp.varscan)));
        
        nfilter = length(SNPiRres.res{snpirSampIdx(i)}.filterName)+1;
        fds = fieldnames(tmp);
        bardata = zeros(3,length(fds),nfilter);
        for fti = 1:nfilter
            if fti == 1
                ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,1) + ...
                    SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,1);
                slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx(ntotalread >= para.rnaMinCov);
            else
                ntotalread = SNPiRres.res{snpirSampIdx(i)}.numReadRef(:,fti-1) + ...
                    SNPiRres.res{snpirSampIdx(i)}.numReadAlt(:,fti-1);
                slocidx = SNPiRres.res{snpirSampIdx(i)}.locidx( ...
                    SNPiRres.res{snpirSampIdx(i)}.validAfterFilter(:,fti-1) & ...
                    ntotalread >= para.rnaMinCov);
            end
            for fdi = 1:length(fds)                
                dlocidx = tmp.(fds{fdi});
                bardata(:,fdi,fti) = [sum(~ismember(dlocidx,slocidx)), ...
                    sum(ismember(dlocidx, slocidx)), sum(~ismember(slocidx, dlocidx))];
            end
        end
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
            th = title(fds{fdi});
            xticklabel = strrep(SNPiRres.res{snpirSampIdx(i)}.filterName, 'rm', '');
            xticklabel = strrep(xticklabel, 'sk', 'rpt');
            set(gca, 'xtick', 1:nfilter, 'xticklabel', [{'none'}, xticklabel]);
            xlabel('added filter');

            if fdi == 2
                thpos = get(th, 'position');
                thext = get(th, 'extent');
                text(thpos(1), thpos(2)+thext(4)*1.3, strrep(overlapsample{i}, 'Sample_', ''),'fontsize',12);
            end
        end
        if do.DRoverlapPercent(2)
            saveas(fh, [para.figdir, strrep(overlapsample{i}, 'Sample_', ''), '.DNA_RNA_ratio.png'], 'png');
        end
    end
end


if exist('do', 'var')
    olddo = do;
end

do.sample = '138381';

do.plotTestRnaExp = {1, 1, '%s_TestRnaExp.%s.png'};
do.plotRnaExp = {0, 1, '%s_RnaExp.%s.png'};
do.plotRnaAlleleCount = {0, 1, '%s_RnaAlleleCount.%s.png'};
do.plotAF = {0, 1, '%s_AF.%s.png'};
do.DnaSom3Caller_Rna2Aligner_Venn = {0, 1, '%s_Dna_Rna_2Alinger.%s.png'}; %
do.Rna_Aligner_Venn = {0, 1, '%s_Rna_Alinger.%s.png'}; %in DNA/RNA; 
do.DnaSom3_Rna4Caller_Ratio = {0, 1, '%s_DnaSom_Rna4Caller.%s.png'}; %in DNA/RNA; bar
do.DnaSom3_RnaSom1_RatioCurve = {0, 1, '%s_DnaSom_RnaSomCaller.%s.png'}; %in DNA/RNA; 
do.DnaSom3Caller_Rna1Caller_Venn = {0, 1, '%s_DnaSom3Caller_vs_RnaPerCaller.%s.venn.png'}; % Venn, subplot by caller (x2, no filter and all filters)
do.DnaSom3Caller_RnaGatk_RatioCurve = {0, 1, '%s_DnaSom_RnaGatk.%s.png'}; % compare Dna somatic calls with RNA GATK, subplot by filter
do.RnaCaller_Venn = {0, 1, '%s_Rna.%s.venn.png'}; %compare calls in RNA between callers, subplot by filter

para.sample = {'138381', '138381-4N', '138381-2T';
    '132540', '132540-10N', '132540-1T' };
para.aligner = '_bwa';
para.caller = {'mutect', 'strelka', 'virmid', 'varscan', 'gatk'};
para.covlimit = [1, 6, 11];
para.RnaUseSomatic = true;
para.DnaUseSomatic = true;
para.globalCov = 1;
para.DnaFilterByRnaCov = para.globalCov;
para.DnaFilterByDnaCov = para.globalCov;
para.RnaGatkFilterByNormCov = para.globalCov;
para.RnaFilterByRnaCov = para.globalCov;
para.filterExonUtr = false;
para.applyDnaFilter = true;
para.DnaTrueSet = {'union', 'any2', 'all'};
para.datatype = {'DNA', 'RNA'};
para.figdir = 'figures/rnavar/DnaRnaComparison/';
para.barcolor = [0.3, 0.3, 0.8; 0.8, 0.4, 0.4; 0.3, 0.7, 0.3]; %only in DNA, overlap, only in RNA
para.figsize = [0,0, 3000, 1200];
fig = figure('position', para.figsize,'visible','off');
set(fig, 'paperpositionmode', 'auto');

if exist('olddo', 'var')
    if ~strcmpi(olddo.sample, do.sample)
        sampi = find(strcmp(para.sample(:,1), do.sample));
        if ~exist('dnavar', 'var');
            load(sprintf('data/dnavar.%s_%s.withfilter.mat', para.sample{sampi,2}, para.sample{sampi,3}));
        end
        if ~exist('rnavar', 'var');
            load(sprintf('data/rnavar.%s_%s.withfilter.mat', para.sample{sampi,2}, para.sample{sampi,3}));
        end
        if ~exist('normalRnaGatk', 'var');
                normalRnaGatk = loadStructData(sprintf('data/snpirres.Sample_%s.mat',para.sample{sampi,2}));
        end
    end
    clear olddo
else
    sampi = find(strcmp(para.sample(:,1), do.sample));    
    load(sprintf('data/dnavar.%s_%s.withfilter.mat', para.sample{sampi,2}, para.sample{sampi,3}));
    load(sprintf('data/rnavar.%s_%s.withfilter.mat', para.sample{sampi,2}, para.sample{sampi,3}));
    normalRnaGatk = loadStructData(sprintf('data/snpirres.Sample_%s.mat',para.sample{sampi,2}));    
end

if para.applyDnaFilter
    fltfiletag = 'fltDR';
else
    fltfiletag = 'fltR';
end
if para.filterExonUtr
    fltfiletag = [fltfiletag 'ExonUtr'];
else
    fltfiletag = [fltfiletag 'AllReg'];
end
if para.RnaUseSomatic
    fltfiletag = [fltfiletag 'RnaSom'];
else
    fltfiletag = [fltfiletag 'RnaAll'];
end
if strcmp(para.aligner, '_bwa')
    fltfiletag = [fltfiletag 'Bwa'];
else
    fltfiletag = [fltfiletag 'Star'];
end

if ~exist('exonUtrLocIdx', 'var')
    load GENCODE.exonutrlocidx.v16.mat
end

exonUtrDnaFilter = ismember(dnavar.locidx, exonUtrLocIdx);
exonUtrRnaFilter = ismember(rnavar.locidx, exonUtrLocIdx);
normalRnaFilter = ~ismember(rnavar.locidx, normalRnaGatk.locidx);
normalDnaFilter = ~ismember(dnavar.locidx, normalRnaGatk.locidx);


if do.plotTestRnaExp{1}    
    
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    %fltIdx = nfilter;        
    rnaAlleleCount = loadStructData('data/rnaAlleleCount.138381.mat');
    rnaAlleleCount = rnaAlleleCount.fromfltbam;
    rnaAlleleCount.DP = sum(rnaAlleleCount.count(:, ~ismember(rnaAlleleCount.ntbase, {'>'})), 2);
    
    rnaExp = loadStructData('data/varMappedGene.138381.dna.rna.union.mat');
    
    for loopidx = 1:length(loopCov) 
        for truesetIdx = 1:3 % one file
            % DNA
            for fltIdx = 2:nfilter
                trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                    para.DnaTrueSet{truesetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                D = dnavar.locidx(trueflt);                
                for rnaSetIdx = truesetIdx %1:4
                    %RNA
                    if rnaSetIdx < 4
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                            para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                    else
                        
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
                            'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);
                    end
                    R = rnavar.locidx( rnaflt );
                    
                    U = union(D, R);
                    inD = ismember(U,D);
                    inR = ismember(U,R);
                    [~, uidx] = ismember(U, rnaAlleleCount.locidx);
                    groupdata = cell(length(U), 1);
                    groupdata(inD & inR) = {'overlap'};
                    groupdata(inD & ~inR) = {'only in DNA'};
                    groupdata(~inD & inR) = {'only in RNA'};
                    
                    DP = zeros(length(uidx), 1);
                    altcount = zeros(length(uidx), 1);
                    rnaexp = zeros(length(uidx), 1);                    
                    DP(uidx~=0) = rnaAlleleCount.DP(uidx(uidx~=0));
                    altcount(uidx~=0) = rnaAlleleCount.altcount(uidx(uidx~=0));
                    [~, expidx] = ismember(U, rnaExp.locidx);
                    rnaexp(expidx~=0) = log2( rnaExp.geneNormalizedCountMax( ...
                        expidx(expidx~=0), strcmp(rnaExp.sample, '138381-2T'))+1);                    
                    subplot(2,3,fltIdx-1);
                    dataidx = ~strcmp(groupdata, 'only in DNA');
                    pval = kruskalwallis_figurecontrol(rnaexp(dataidx), groupdata(dataidx), fig);
                    title(sprintf('%s, pval=%0.3g', ...
                        dnavar.filterName{validFilterIdx(fltIdx)}, pval), 'fontsize', 14);                    

                    ylabel( 'log2 gene normalized count', 'fontsize', 14);
                end
            end
            if do.plotTestRnaExp{2}
                fltfiletag = sprintf('DnaSom%s_Cov%d%s',para.DnaTrueSet{truesetIdx},loopCov(loopidx),para.aligner);
                saveas(fig, sprintf(['%s/' do.plotTestRnaExp{3}],para.figdir, ...
                    do.sample,sprintf('%s',fltfiletag)), 'png');
            end
        end
    end
end


if do.plotRnaExp{1}    
    groupColor = {'only in DNA', para.barcolor(1,:); ...
        'overlap', para.barcolor(2,:); ...        
        'only in RNA', para.barcolor(3,:) };
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    %fltIdx = nfilter;        
    rnaAlleleCount = loadStructData('data/rnaAlleleCount.138381.mat');
    rnaAlleleCount = rnaAlleleCount.fromfltbam;
    rnaAlleleCount.DP = sum(rnaAlleleCount.count(:, ~ismember(rnaAlleleCount.ntbase, {'>'})), 2);
    
    rnaExp = loadStructData('data/varMappedGene.138381.dna.rna.union.mat');
    
    for loopidx = 1:length(loopCov) 
        for truesetIdx = 1:3 % one file
            % DNA
            for fltIdx = 2:nfilter
                trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                    para.DnaTrueSet{truesetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                D = dnavar.locidx(trueflt);                
                for rnaSetIdx = truesetIdx %1:4
                    %RNA
                    if rnaSetIdx < 4
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                            para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                    else
                        
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
                            'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);
                    end
                    R = rnavar.locidx( rnaflt );
                    
                    U = union(D, R);
                    inD = ismember(U,D);
                    inR = ismember(U,R);
                    [~, uidx] = ismember(U, rnaAlleleCount.locidx);
                    groupdata = cell(length(U), 1);
                    groupdata(inD & inR) = {'overlap'};
                    groupdata(inD & ~inR) = {'only in DNA'};
                    groupdata(~inD & inR) = {'only in RNA'};
                    
                    DP = zeros(length(uidx), 1);
                    altcount = zeros(length(uidx), 1);
                    rnaexp = zeros(length(uidx), 1);                    
                    DP(uidx~=0) = rnaAlleleCount.DP(uidx(uidx~=0));
                    altcount(uidx~=0) = rnaAlleleCount.altcount(uidx(uidx~=0));
                    [~, expidx] = ismember(U, rnaExp.locidx);
                    rnaexp(expidx~=0) = log2( rnaExp.geneNormalizedCountMax( ...
                        expidx(expidx~=0), strcmp(rnaExp.sample, '138381-2T'))+1);                    
                    subplot(2,3,fltIdx-1);
                    %h = gscatter(DP, altcount, groupdata);
                    h = gscatter(rnaexp, altcount, groupdata);
                    legend('location', 'NW');
                    for hidx = 1:length(h)
                        set(h(hidx), 'color', groupColor{strcmp(get(h(hidx),'DisplayName'), groupColor(:,1)),2});
                    end
                    title(dnavar.filterName{validFilterIdx(fltIdx)}, 'fontsize', 14);                    

                    xlabel( 'log2 gene normalized count', 'fontsize', 14);
                    ylabel( '# ALT read', 'fontsize', 14);
                    
                end
            end
            if do.plotRnaExp{2}
                fltfiletag = sprintf('DnaSom%s_Cov%d%s',para.DnaTrueSet{truesetIdx},loopCov(loopidx),para.aligner);
                saveas(fig, sprintf(['%s/' do.plotRnaExp{3}],para.figdir, ...
                    do.sample,sprintf('%s',fltfiletag)), 'png');
            end
        end
    end
end


if do.plotRnaAlleleCount{1}    
    groupColor = {'only in DNA', para.barcolor(1,:); ...
        'overlap', para.barcolor(2,:); ...        
        'only in RNA', para.barcolor(3,:) };
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    %fltIdx = nfilter;        
    rnaAlleleCount = loadStructData('data/rnaAlleleCount.138381.mat');
    rnaAlleleCount = rnaAlleleCount.fromfltbam;
    rnaAlleleCount.DP = sum(rnaAlleleCount.count(:, ~ismember(rnaAlleleCount.ntbase, {'>'})), 2);
    
    for loopidx = 1:length(loopCov) 
        for truesetIdx = 1:3 % one file
            % DNA
            for fltIdx = 2:nfilter
                trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                    para.DnaTrueSet{truesetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                D = dnavar.locidx(trueflt);                
                for rnaSetIdx = truesetIdx %1:4
                    %RNA
                    if rnaSetIdx < 4
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                            para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                    else
                        
                        rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
                            'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);
                    end
                    R = rnavar.locidx( rnaflt );
                    
                    U = union(D, R);
                    inD = ismember(U,D);
                    inR = ismember(U,R);
                    [~, uidx] = ismember(U, rnaAlleleCount.locidx);
                    groupdata = cell(length(U), 1);
                    groupdata(inD & inR) = {'overlap'};
                    groupdata(inD & ~inR) = {'only in DNA'};
                    groupdata(~inD & inR) = {'only in RNA'};
                    
                    DP = zeros(length(uidx), 1);
                    altcount = zeros(length(uidx), 1);
                    DP(uidx~=0) = rnaAlleleCount.DP(uidx(uidx~=0));
                    altcount(uidx~=0) = rnaAlleleCount.altcount(uidx(uidx~=0));
                    %subplot(4,4, (rnaSetIdx-1)*4+plotidx)
                    subplot(2,3,fltIdx-1);
                    h = gscatter(DP, altcount, groupdata);
                    legend('location', 'NW');
                    for hidx = 1:length(h)
                        set(h(hidx), 'color', groupColor{strcmp(get(h(hidx),'DisplayName'), groupColor(:,1)),2});
                    end
                    title(dnavar.filterName{validFilterIdx(fltIdx)}, 'fontsize', 14);                    
%                     if rnaSetIdx < 4
%                         title(sprintf('RNA %s', para.DnaTrueSet{rnaSetIdx}), 'fontsize', 14);
%                     else
%                         title('RNA any3', 'fontsize', 14);
%                     end
                    xlabel( '# total read', 'fontsize', 14);
                    ylabel( '# ALT read', 'fontsize', 14);
                    
                end
            end
            if do.plotRnaAlleleCount{2}
                fltfiletag = sprintf('DnaSom%s_Cov%d%s',para.DnaTrueSet{truesetIdx},loopCov(loopidx),para.aligner);
                saveas(fig, sprintf(['%s/' do.plotRnaAlleleCount{3}],para.figdir, ...
                    do.sample,sprintf('%s',fltfiletag)), 'png');
            end
        end
    end
end

if do.plotAF{1}        
    groupColor = {'truely missed in RNA', para.barcolor(1,:); ...
        'called in RNA', para.barcolor(2,:); ...
        'no ALT RNA reads', [0.6, 0.6, 0.6]; ...
        'truely missed in DNA', para.barcolor(3,:); ...
        'called in DNA', para.barcolor(2,:); ...
        'no ALT DNA reads', [0.6, 0.6, 0.6]};
    ylabelName = {'RNA AF', '# RNA ALT reads', 'DNA AF', '# DNA ALT reads'};
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    fltIdx = nfilter;        
    DnaCallRnaCount = loadStructData('data/rnacount.som.138381-2T.mat');    
    RnaCallDnaCount = loadStructData('data/dnacount.rnacall.138381-2T.mat');
    DnaCallRnaCount.AF = DnaCallRnaCount.altcount ./ (DnaCallRnaCount.refcount + DnaCallRnaCount.altcount);
    RnaCallDnaCount.AF = RnaCallDnaCount.altcount ./ (RnaCallDnaCount.refcount + RnaCallDnaCount.altcount);
    for loopidx = 2%1:length(loopCov) 
        for truesetIdx = 1:3 % one file
            % DNA
            trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                para.DnaTrueSet{truesetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);            
            D = dnavar.locidx(trueflt);      
            
            [~, didx] = ismember(D, DnaCallRnaCount.locidx);
            if any(didx == 0)
                error('unknown rna calls\n');
            end
            for rnaSetIdx = 1:4
                %RNA
                if rnaSetIdx < 4
                    rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                        para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);
                else
                
                    rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
                        'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);                
                end
                R = rnavar.locidx( rnaflt );
                [~, ridx] = ismember(R, RnaCallDnaCount.locidx);
                if any(ridx == 0)
                    error('unkown rna calls\n');
                end
                
                plotdata = {DnaCallRnaCount.AF(didx), DnaCallRnaCount.altcount(didx), ...
                    RnaCallDnaCount.AF(ridx), RnaCallDnaCount.altcount(ridx)};
                groupdata = cell(1,2);
                groupdata{1} = cell(length(D), 1);
                groupdata{1}(:) = {'truely missed in RNA'};
                groupdata{1}(ismember(D,R)) = {'called in RNA'};
                groupdata{1}(DnaCallRnaCount.AF(didx)==0) = {'no ALT RNA reads'};
                groupdata{2} = cell(length(R), 1);
                groupdata{2}(:) = {'truely missed in DNA'};
                groupdata{2}(ismember(R,D)) = {'called in DNA'};
                groupdata{2}(RnaCallDnaCount.AF(ridx)==0) = {'no ALT DNA reads'};
                for plotidx = 1:length(plotdata)                
                    subplot(4,4, (rnaSetIdx-1)*4+plotidx)
                    [~, vidx] = sort(plotdata{plotidx});
                    h = gscatter(1:length(plotdata{plotidx}), plotdata{plotidx}(vidx), groupdata{ceil(plotidx/2)}(vidx));
                    legend('location', 'NW');
                    for hidx = 1:length(h)
                        set(h(hidx), 'color', groupColor{strcmp(get(h(hidx),'DisplayName'), groupColor(:,1)),2});
                    end
                    if rnaSetIdx < 4
                        title(sprintf('RNA %s', para.DnaTrueSet{rnaSetIdx}), 'fontsize', 14);
                    else
                        title('RNA any3', 'fontsize', 14);
                    end
                    ylabel( ylabelName{plotidx}, 'fontsize', 14);
                end                
            end
            if do.plotAF{2}
                fltfiletag = sprintf('DnaSom%s_Cov%d%s',para.DnaTrueSet{truesetIdx},loopCov(loopidx),para.aligner);
                saveas(fig, sprintf(['%s/' do.plotAF{3}],para.figdir, ...
                    do.sample,sprintf('%s',fltfiletag)), 'png');
            end
        end
    end
end

if do.DnaSom3_Rna4Caller_Ratio{1}        
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    fltIdx = nfilter;    
    rnaSetName = [strcat(para.DnaTrueSet, ' som'), {'union+GATK','any3'}];
    for loopidx = 1:length(loopCov) %one row        
        for truesetIdx = 1:3 %subplot col
            countdata = zeros(length(rnaSetName), 3); %RNA somatic (union, any2, all), somatic-union + GATK  
            % DNA
            trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                para.DnaTrueSet{truesetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);            
            for rnaSetIdx = 1:3
                %RNA
                rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                    para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), validFilterIdx(fltIdx), true);            
                overlap = sum(ismember(dnavar.locidx(trueflt), rnavar.locidx(rnaflt)));
                countdata(rnaSetIdx, :) = [sum(trueflt)-overlap, overlap, sum(rnaflt)-overlap];
                if strcmp(para.DnaTrueSet{rnaSetIdx}, 'union')
                     %GATK
                     gatk = VarFilter.CovAndNumReadFilter(rnavar, ['gatk', para.aligner], loopCov(loopidx), validFilterIdx(fltIdx));
                     gatk = gatk & normalRnaFilter;
                     if validFilterIdx(fltIdx) == 1
                         gatk = gatk & VarFilter.CumValidAfterFilter(rnavar, ['gatk', para.aligner], validFilterIdx(fltIdx), false) ;
                     else
                         gatk = gatk & VarFilter.CumValidAfterFilter(rnavar, ['gatk', para.aligner], validFilterIdx(fltIdx), true) ;
                     end
                    rnaflt = rnaflt & gatk;
                    overlap = sum(ismember(dnavar.locidx(trueflt), rnavar.locidx(rnaflt)));
                    countdata(4, :) = [sum(trueflt)-overlap, overlap, sum(rnaflt)-overlap];
                end
            end
            rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
                'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);
            overlap = sum(ismember(dnavar.locidx(trueflt), rnavar.locidx(rnaflt)));
            countdata(5, :) = [sum(trueflt)-overlap, overlap, sum(rnaflt)-overlap];
            
%             rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller([1:3 5]), para.aligner), ...
%                 'any3', loopCov(loopidx), validFilterIdx(fltIdx), true);
%             rnaextra = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), '_star'), ...
%                 'all', loopCov(loopidx)+5, validFilterIdx(fltIdx), false) & ...
%                 VarFilter.TrueSet(rnavar, strcat('gatk_star'), ...
%                 'all', loopCov(loopidx), validFilterIdx(fltIdx), true) & ...
%                 VarFilter.TrueSet(rnavar, strcat('gatk_bwa'), ...
%                 'all', loopCov(loopidx), validFilterIdx(fltIdx), true) & ...
%                 normalRnaFilter & exonUtrRnaFilter;
%             rnaflt = rnaflt | rnaextra;
%             overlap = sum(ismember(dnavar.locidx(trueflt), rnavar.locidx(rnaflt)));
%             countdata(6, :) = [sum(trueflt)-overlap, overlap, sum(rnaflt)-overlap];
%             
%             rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
%                 'any2', loopCov(loopidx), validFilterIdx(fltIdx), true);
%             rnaextra = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
%                 'union', loopCov(loopidx)+5, validFilterIdx(fltIdx), true) & ...
%                 VarFilter.TrueSet(rnavar, strcat(para.caller(5), para.aligner), ...
%                 'all', loopCov(loopidx)+5, validFilterIdx(fltIdx), true) & normalRnaFilter & exonUtrRnaFilter;
%             rnaflt = rnaflt | rnaextra;
%             overlap = sum(ismember(dnavar.locidx(trueflt), rnavar.locidx(rnaflt)));
%             countdata(7, :) = [sum(trueflt)-overlap, overlap, sum(rnaflt)-overlap];
            
            subplot(2, 3, (loopidx-1)*3+truesetIdx);
            h = bar(1:size(countdata,1), countdata, 'stacked');
            set(h, 'barwidth', 0.6, 'edgecolor','none');
            if loopidx == 2 && truesetIdx == 1
                legend({'DNA', 'overlap', 'RNA'},'fontsize', 16); %,'location','SE');
            end
            for hidx = 1:length(h)
                set(h(hidx), 'facecolor', para.barcolor(hidx,:));
            end
            xlim([0.5 size(countdata,1)+0.5]);
            set(gca, 'xtick', 1:size(countdata,1), 'xticklabel', rnaSetName, 'fontsize', 16);
            ylabel('# variant', 'fontsize', 14);                        
            title(sprintf('DNA %s, Cov %d', para.DnaTrueSet{truesetIdx}, loopCov(loopidx)), 'fontsize', 16);
            fprintf('DNA %s, Cov %d\n', para.DnaTrueSet{truesetIdx}, loopCov(loopidx));
            disp(bsxfun(@rdivide, countdata(:,2), bsxfun(@plus, countdata(:,2), countdata(:, [1 3]))));
        end
    end
    if  do.DnaSom3_Rna4Caller_Ratio{2}        
        saveas(fig, sprintf(['%s/' do.DnaSom3_Rna4Caller_Ratio{3}],para.figdir, ...
            do.sample,sprintf('%srate',fltfiletag)), 'png');        
    end
        
end

if do.DnaSom3_RnaSom1_RatioCurve{1}
    set(fig, 'position', [0,0, 2400, 1600]);
    fig(2) = figure('position', [0,0, 2400, 1600],'visible','off');
    set(fig(2), 'paperpositionmode', 'auto');
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    for loopidx = 1:length(loopCov) %one file        
        for cidx = 1:4 %subplot row            
            rnafd = [para.caller{cidx} para.aligner];
            for truesetIdx = 1:4 %subplot col
                ratiodata = zeros(nfilter, 4, 2); %dim2: all, som, allexonutr, somexonutr; dim3:DNA, RNA
                %somatic or all only apply to RNA; DNA always uses somatic
                %calls
                for fltIdx = 1:nfilter
                    % DNA
                    if truesetIdx < 4
                        trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                            para.DnaTrueSet{truesetIdx}, loopCov(loopidx), max(2,validFilterIdx(fltIdx)), para.DnaUseSomatic);
                    else %compare to the somatic caller
                        trueflt = VarFilter.TrueSet(dnavar, rnafd, ...
                            'all', loopCov(loopidx), max(2,validFilterIdx(fltIdx)), para.DnaUseSomatic);
                    end
                                        
                    %RNA
                    rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx));
                    
                    for useSom = 0:1
                        D = trueflt;
                        R = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), useSom) ;
                        overlap = sum(ismember(dnavar.locidx(D), rnavar.locidx(R)));
                        ratiodata(fltIdx, useSom+1, :) = overlap ./ [sum(D), sum(R)];
                        
                        D = trueflt & exonUtrDnaFilter ;
                        R = rnaflt & exonUtrRnaFilter & ...
                            VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), useSom) ;
                        overlap = sum(ismember(dnavar.locidx(D), rnavar.locidx(R)));
                        ratiodata(fltIdx, useSom+3, :) = overlap ./ [sum(D), sum(R)];
                    end                        
                end
                for datatype = 1:2 %DNA, RNA
                    set(0,'CurrentFigure',fig(datatype));
                    subplot(4, 4, (cidx-1)*4+truesetIdx);
                    plotInDnaRnaRatio(ratiodata(:,:,datatype), [], ...
                        rnavar.filterName(validFilterIdx), ...
                        'ylabel', sprintf('%% in %s calls', para.datatype{datatype}), ... 
                        'linename', {'all', 'som', 'ExonUtr', 'som ExonUtr'});
%                     plotInDnaRnaRatio(ratiodata(:,:,1), ratiodata(:,:,2), ...
%                         rnavar.filterName(validFilterIdx));%, 'linename', {'all', 'som', 'ExonUtr', 'som ExonUtr'});
                    if truesetIdx < 4
                        title(sprintf('DNA %s, RNA %s', para.DnaTrueSet{truesetIdx}, para.caller{cidx}), 'fontsize', 14);
                    else
                        title(sprintf('DNA %s, RNA %s', para.caller{cidx}, para.caller{cidx}), 'fontsize', 14);
                    end
                end
            end
        end
        if  do.DnaSom3_RnaSom1_RatioCurve{2}
            for datatype = 1:2
                saveas(fig(datatype), sprintf(['%s/' do.DnaSom3_RnaSom1_RatioCurve{3}],para.figdir, ...
                    do.sample,sprintf('%sCov%d.%srate',fltfiletag,loopCov(loopidx),para.datatype{datatype})), 'png');
            end
        end
    end
    set(fig, 'position', para.figsize);
    close(fig(2));
    fig(2) = [];
end

if do.DnaSom3Caller_Rna1Caller_Venn{1}  % Venn, subplot by caller (filter per file)
    if strcmp(para.aligner, '_bwa')
        prefix = 'Bwa';
    else
        prefix = 'Star';
    end
    usesom_exonutr = [false, false; true, false; false, true; true, true];    
    exontag = {'AllReg', 'ExonUtr'};
    somtag = {'RnaAll', 'RnaSom'};
    setdata = cell(1,4);        
    loopCov = [1, 6];    
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2); % 2 is somatic
    nfilter = length(validFilterIdx);    
    for loopidx = 1:length(loopCov) %one file    
        for somexonidx = 1:4
            para.RnaUseSomatic = usesom_exonutr(somexonidx, 1);
            para.filterExonUtr = usesom_exonutr(somexonidx, 2);            
            for fltIdx = 1:nfilter
                fltfiletag = [prefix, exontag{para.filterExonUtr+1}, ...
                    somtag{para.RnaUseSomatic+1}, sprintf('Cov%d',loopCov(loopidx)), ...
                    rnavar.filterName{validFilterIdx(fltIdx)}];
                for i = 1:3
                    dnafd = [para.caller{i}, para.aligner];
                    flt = VarFilter.CovAndNumReadFilter(dnavar, dnafd, loopCov(loopidx), max(2,validFilterIdx(fltIdx)));
                    flt = flt & VarFilter.CumValidAfterFilter(dnavar, dnafd, max(2,validFilterIdx(fltIdx)), para.DnaUseSomatic);
                    if para.filterExonUtr
                        flt = flt & exonUtrDnaFilter;
                    end
                    setdata{i} = dnavar.locidx( flt );
                end

                for i = 1:length(para.caller)
                    rnafd = [para.caller{i} para.aligner];
                    if strcmp(para.caller{i}, 'gatk')
                        rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx));
                        rnaflt = rnaflt & normalRnaFilter ;
                    else
                        rnaflt = VarFilter.CovAndNumReadFilter(rnavar,rnafd, loopCov(loopidx), validFilterIdx(fltIdx));
                    end
                    if para.RnaUseSomatic
                        rnaflt = rnaflt & rnavar.(rnafd).validAfterFilter(:,2);
                    end
                    if para.filterExonUtr
                        rnaflt = rnaflt & exonUtrRnaFilter;
                    end

                    if strcmp(para.caller{i}, 'gatk') && validFilterIdx(fltIdx) == 1
                            rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar,rnafd,nfilter,false);                        
                    else
                        rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar,rnafd,nfilter,para.RnaUseSomatic);
                    end
                    
                    setdata{4} = rnavar.locidx( rnaflt );
                    subplot(2,3,i);
                    fixvenn(setdata, 'label', [para.caller(1:3), {'rna'}],'numfontsize',14,'labelfontsize',14);
                    title([para.caller{i} '+' rnavar.filterName{validFilterIdx(fltIdx)}], 'fontsize',14,'fontweight','bold');
                end
                if do.DnaSom3Caller_Rna1Caller_Venn{2}                    
                    saveas(fig, sprintf(['%s/' do.DnaSom3Caller_Rna1Caller_Venn{3}],para.figdir, do.sample, fltfiletag), 'png');
                end
            end            
        end        
    end
end

if do.DnaSom3Caller_RnaGatk_RatioCurve{1}
    set(fig, 'position', [0,0, 1800, 300]);
    fig(2) = figure('position', [0,0, 1800, 300],'visible','off');
    set(fig(2), 'paperpositionmode', 'auto');
    loopCov = [1, 6];    
    validFilterIdx = 1:length(dnavar.filterName); % 2 is somatic
    nfilter = length(validFilterIdx);    
    rnafd = ['gatk' para.aligner];
    for loopidx = 1:length(loopCov) %one file                        
        for truesetIdx = 1:3 %subplot col
            ratiodata = zeros(nfilter, 4, 2); %dim2: all, som, allexonutr, somexonutr; dim3:DNA, RNA
            ratiodata(:, [2 4], :) = NaN;            
            for fltIdx = 1:nfilter
                % DNA                
                trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                    para.DnaTrueSet{truesetIdx}, loopCov(loopidx), max(2,validFilterIdx(fltIdx)), para.DnaUseSomatic);                                
                %RNA
                rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx));
                rnaflt = rnaflt & normalRnaFilter;
                
                D = trueflt;
                if validFilterIdx(fltIdx) == 1
                    R = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), false) ;
                else
                    R = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), true) ;
                end
                overlap = sum(ismember(dnavar.locidx(D), rnavar.locidx(R)));
                ratiodata(fltIdx, 1, :) = overlap ./ [sum(D), sum(R)];
                
                D = trueflt & exonUtrDnaFilter ;
                if validFilterIdx(fltIdx) == 1
                    R = rnaflt & exonUtrRnaFilter & ...
                        VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), false) ;
                else
                    R = rnaflt & exonUtrRnaFilter & ...
                        VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), true) ;
                end
                overlap = sum(ismember(dnavar.locidx(D), rnavar.locidx(R)));
                ratiodata(fltIdx, 3, :) = overlap ./ [sum(D), sum(R)];                
            end
            for datatype = 1:2 %DNA, RNA
                set(0,'CurrentFigure',fig(datatype));
                subplot(1, 3, truesetIdx);
                plotInDnaRnaRatio(ratiodata(:,:,datatype), [], ...
                    strrep(rnavar.filterName(validFilterIdx), 'som', 'qual'), ...
                    'ylabel', sprintf('%% in %s calls', para.datatype{datatype}));
                    %'linename', {'all', 'som', 'ExonUtr', 'som ExonUtr'});                
                title(sprintf('DNA %s, RNA GATK', para.DnaTrueSet{truesetIdx}), 'fontsize', 14);                
            end
        end        
        if  do.DnaSom3Caller_RnaGatk_RatioCurve{2}
            for datatype = 1:2
                saveas(fig(datatype), sprintf(['%s/' do.DnaSom3Caller_RnaGatk_RatioCurve{3}],para.figdir, ...
                    do.sample,sprintf('%sCov%d.%srate',fltfiletag,loopCov(loopidx),para.datatype{datatype})), 'png');
            end
        end
    end
    set(fig, 'position', para.figsize);
    close(fig(2));
    fig(2) = [];
end


if do.RnaCaller_Venn{1} %compare calls in RNA between callers, subplot by filter    
    set(fig, 'position', [0,0, 2400, 1600]);
    if strcmp(para.aligner, '_bwa')
        prefix = 'Bwa';
    else
        prefix = 'Star';
    end
    usesom_exonutr = [false, false; true, false; false, true; true, true];    
    exontag = {'AllReg', 'ExonUtr'};
    somtag = {'RnaAll', 'RnaSom'};
    loopCov = [1 6];
    validFilterIdx = setdiff(1:length(dnavar.filterName), 2);
    nfilter = length(validFilterIdx);
    for loopidx = 1:length(loopCov)
        for somexonidx = 1:4
            para.RnaUseSomatic = usesom_exonutr(somexonidx, 1);
            para.filterExonUtr = usesom_exonutr(somexonidx, 2);
            fltfiletag = [prefix, exontag{para.filterExonUtr+1}, somtag{para.RnaUseSomatic+1}, sprintf('Cov%d',loopCov(loopidx))];
            setdata = cell(1,4);
            for fltIdx = 1:nfilter
                for cidx = 1:3
                    rnafd = [para.caller{cidx} para.aligner];
                    rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx)); 
                    rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), para.RnaUseSomatic) ;                    
                    if para.filterExonUtr
                        rnaflt = rnaflt & exonUtrRnaFilter;
                    end
                    setdata{cidx} = rnavar.locidx(rnaflt);
                end
                
                %GATK
                rnafd = [para.caller{end} para.aligner];
                rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx));                 
                rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), true);
                if para.RnaUseSomatic
                    rnaflt = rnaflt & normalRnaFilter;
                end                
                if para.filterExonUtr
                    rnaflt = rnaflt & exonUtrRnaFilter;
                end
                setdata{4} = rnavar.locidx( rnaflt );
                subplot(2,4, fltIdx);
                fixvenn(setdata, 'label', [para.caller(1:3), {'gatk'}],'numfontsize',14,'labelfontsize',14);                
                title(rnavar.filterName{validFilterIdx(fltIdx)}, 'fontsize',14,'fontweight','bold');                
            end
            if do.RnaCaller_Venn{2}
                saveas(fig, sprintf(['%s/' do.RnaCaller_Venn{3}],para.figdir, do.sample,fltfiletag), 'png');
            end
        end
    end
    set(fig, 'position', para.figsize);
end

if do.Rna_Aligner_Venn{1} %compare calls in RNA between callers, subplot by filter    
    set(fig, 'position', [0,0, 2400, 1600]);    
    usesom_exonutr = [false, false; true, false; false, true; true, true];    
    exontag = {'AllReg', 'ExonUtr'};
    somtag = {'RnaAll', 'RnaSom'};
    loopCov = [1 6];    
    nfilter = length(rnavar.filterName);
    callerIdx = [1:3 5];
    aligners = {'_bwa', '_star'};
    for loopidx = 1:length(loopCov)
        for somexonidx = 2%1:4
            para.RnaUseSomatic = usesom_exonutr(somexonidx, 1);
            para.filterExonUtr = usesom_exonutr(somexonidx, 2);
            fltfiletag = [exontag{para.filterExonUtr+1}, somtag{para.RnaUseSomatic+1}, sprintf('Cov%d',loopCov(loopidx))];
            setdata = cell(1,2); 
            set4data = cell(1,4);
            for cidx = 1:length(callerIdx)
                for aidx = 1:2
                    rnafd = [para.caller{callerIdx(cidx)} aligners{aidx}];
                    rnaflt = VarFilter.TrueSet(rnavar, rnafd, ...
                        'all', loopCov(loopidx), nfilter, para.RnaUseSomatic);
                    if para.filterExonUtr
                        rnaflt = rnaflt & exonUtrRnaFilter;
                    end
                    if strcmp(para.caller{callerIdx(cidx)}, 'gatk')
                        rnaflt = rnaflt & normalRnaFilter;
                        if aidx == 1
                            set4data{1} = rnavar.locidx( rnaflt);
                        else
                            set4data{4} = rnavar.locidx( rnaflt);
                        end
                    end                    
                    setdata{aidx} = rnavar.locidx( rnaflt );
                end
                subplot(2,4,cidx);
                fixvenn(setdata, 'label', strrep(aligners, '_', ''),'numfontsize',14,'labelfontsize',14);
                title(para.caller{callerIdx(cidx)}, 'fontsize', 14);
            end
            
            for truesetIdx = 1:3
                for aidx = 1:2
                    rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), aligners{aidx}), ...
                        para.DnaTrueSet{truesetIdx}, loopCov(loopidx), nfilter, para.RnaUseSomatic);
                    if para.filterExonUtr
                        rnaflt = rnaflt & exonUtrRnaFilter;
                    end
                    setdata{aidx} = rnavar.locidx(rnaflt);                    
                    if strcmp(para.DnaTrueSet{truesetIdx}, 'union')
                        set4data{aidx+1} = rnavar.locidx( rnaflt );
                    end
                end
                subplot(2,4, length(callerIdx)+truesetIdx);
                fixvenn(setdata, 'label', strrep(aligners, '_', ''),'numfontsize',14,'labelfontsize',14);
                title(para.DnaTrueSet{truesetIdx}, 'fontsize', 14);
            end
            
            subplot(2,4,8);
            fixvenn(set4data, 'label', {'gatk bwa', 'union bwa', 'union star', 'gatk star'},'numfontsize',14,'labelfontsize',14);            
            if do.Rna_Aligner_Venn{2}
                saveas(fig, sprintf(['%s/' do.Rna_Aligner_Venn{3}],para.figdir, do.sample,fltfiletag), 'png');
            end
        end
    end
    set(fig, 'position', para.figsize);
end


if do.DnaSom3Caller_Rna2Aligner_Venn{1}  % Venn
    usesom_exonutr = [false, false; true, false; false, true; true, true];    
    exontag = {'AllReg', 'ExonUtr'};
    somtag = {'RnaAll', 'RnaSom'};
    setdata = cell(1,4);        
    loopCov = [1, 6];        
    nfilter = length(rnavar.filterName);
    for loopidx = 1:length(loopCov) %one file    
        for somexonidx = 2%1:4
            para.RnaUseSomatic = usesom_exonutr(somexonidx, 1);
            para.filterExonUtr = usesom_exonutr(somexonidx, 2);                        
            fltfiletag = [exontag{para.filterExonUtr+1}, ...
                somtag{para.RnaUseSomatic+1}, sprintf('Cov%d',loopCov(loopidx)), ...
                rnavar.filterName{nfilter}];
            setdata = cell(1,4);
            setlabel = {'any2', 'union', 'gatk', ''};
            % DNA
            trueflt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), ...
                setlabel{1}, loopCov(loopidx), nfilter, true);  
            if para.filterExonUtr
                truelft = trueflt & exonUtrDnaFilter;
            end
            setdata{1} = dnavar.locidx( trueflt );
                        
            rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), para.aligner), ...
                setlabel{2}, loopCov(loopidx), nfilter, para.RnaUseSomatic );
            if para.filterExonUtr
                rnaflt = rnaflt & exonUtrRnaFilter;
            end
            setdata{2} = rnavar.locidx( rnaflt );
            setlabel{2} = [setlabel{2} strrep(para.aligner, '_', ' ')];
            
            rnaflt = VarFilter.TrueSet(rnavar, ['gatk', para.aligner], ...
                'all', loopCov(loopidx), nfilter, para.RnaUseSomatic ) & normalRnaFilter;
            if para.filterExonUtr
                rnaflt = rnaflt & exonUtrRnaFilter;
            end
            setdata{3} = rnavar.locidx( rnaflt );
            setlabel{3} = [setlabel{3} strrep(para.aligner, '_', ' ')];
            
            for rnaSetIdx = 1:3            
                %RNA
                rnaflt = VarFilter.TrueSet(rnavar, strcat(para.caller(1:3), '_star'), ...
                    para.DnaTrueSet{rnaSetIdx}, loopCov(loopidx), nfilter, para.RnaUseSomatic);            
                if para.filterExonUtr
                    rnaflt = rnaflt & exonUtrRnaFilter;
                end
                setdata{4} = rnavar.locidx( rnaflt );
                setlabel{4} = [para.DnaTrueSet{rnaSetIdx} ' star'];
                subplot(2,2,rnaSetIdx);
                fixvenn(setdata, 'label', setlabel,'numfontsize',14,'labelfontsize',14);                
            end
            
            rnaflt = VarFilter.TrueSet(rnavar, 'gatk_star', ...
                'all', loopCov(loopidx), nfilter, para.RnaUseSomatic ) & normalRnaFilter;
            if para.filterExonUtr
                rnaflt = rnaflt & exonUtrRnaFilter;
            end
            setdata{4} = rnavar.locidx( rnaflt );
            setlabel{4} = 'gatk_star';
            subplot(2,2,4);
            fixvenn(setdata, 'label', setlabel, 'numfontsize', 14, 'labelfontsize', 14);
                
            if do.DnaSom3Caller_Rna2Aligner_Venn{2}
                saveas(fig, sprintf(['%s/' do.DnaSom3Caller_Rna2Aligner_Venn{3}],para.figdir, do.sample, fltfiletag), 'png');
            end            
        end        
    end
end


close(fig);


% if do.DnaSom_RnaGatk{1} % Venn, compare Dna somatic calls with RNA GATK, subplot by filter
%     fig(2) = figure('position', get(fig(1), 'position'),'visible','off');
%     set(fig(2), 'paperpositionmode', 'auto');
%     clf
%     if strcmp(para.aligner, '_bwa')
%         prefix = 'Bwa';
%     else
%         prefix = 'Star';
%     end
%     setdata = cell(1,4);            
%     loopCov = [1, 6];
%     use_exonutr = [false; true];
%     exontag = {'AllReg', 'ExonUtr'};    
%     rnafd = ['gatk' para.aligner];
%     validFilterIdx = 1:length(rnavar.filterName);
%     nfilter = length(validFilterIdx);
%     for loopidx = 1:length(loopCov)
%         for somexonidx = 1:size(use_exonutr, 1)  %per file
%             para.filterExonUtr = use_exonutr(somexonidx);
%             if ~para.applyDnaFilter        
%                 for di = 1:3
%                     dnafd = [para.caller{di}, para.aligner];
%                     flt = VarFilter.CovAndNumReadFilter(dnavar, dnafd, loopCov(loopidx), 2);                    
%                     flt = flt & VarFilter.CumValidAfterFilter(dnavar,dnafd,2,para.DnaUseSomatic);                    
%                     setdata{di} = dnavar.locidx( flt );
%                 end
%             end
%             for fltIdx = 1:nfilter
%                 if para.applyDnaFilter
%                     for di = 1:3
%                         dnafd = [para.caller{di}, para.aligner];
%                         flt = VarFilter.CovAndNumReadFilter(dnavar, dnafd, loopCov(loopidx), max(2,validFilterIdx(fltIdx)));                             
%                         flt = flt & VarFilter.CumValidAfterFilter(dnavar, dnafd, max(2,validFilterIdx(fltIdx)), para.DnaUseSomatic);                        
%                         if para.filterExonUtr
%                             flt = flt & exonUtrDnaFilter;
%                         end
%                         setdata{di} = dnavar.locidx( flt );
%                     end
%                 end
%                 
%                 rnaflt = VarFilter.CovAndNumReadFilter(rnavar, rnafd, loopCov(loopidx), validFilterIdx(fltIdx));
%                 rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, validFilterIdx(fltIdx), true);
%                 rnaflt = rnaflt & normalRnaFilter;
%                 if para.filterExonUtr
%                     rnaflt = rnaflt & exonUtrRnaFilter;
%                 end
%                 
%                 setdata{4} = rnavar.locidx( rnaflt );
%                 subplot(2,4,fltIdx);
%                 fixvenn(setdata, 'label', [para.caller(1:3), {'gatk'}],'numfontsize',14,'labelfontsize',14);                
%                 if strcmp(validFilterIdx(fltIdx), 'som')
%                     title('qual', 'fontsize',14,'fontweight','bold');
%                 else
%                     title(rnavar.filterName{validFilterIdx(fltIdx)}, 'fontsize',14,'fontweight','bold');
%                 end                
%             end
%             if do.DnaSom_RnaGatk{2}
%                 fltfiletag = [prefix, exontag{para.filterExonUtr+1}, sprintf('Cov%d',loopCov(loopidx))];
%                 saveas(fig, sprintf(['%s/' do.DnaSom_RnaGatk{3}],para.figdir, do.sample, fltfiletag), 'png');
%             end
% 
%             clf            
%             for i = 1:length(para.DnaTrueSet)
%                 if ~para.applyDnaFilter
%                     flt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), para.DnaTrueSet{i});
%                     for di = 1:3
%                         dnafd = [para.caller{di}, para.aligner];
%                         flt = VarFilter.CovAndNumReadFilter(dnavar, dnafd, loopCov(loopidx), 2);  
%                         flt = flt & VarFilter.CumValidAfterFilter(dnavar, dnafd, 2, para.DnaUseSomatic);                        
%                         if para.filterExonUtr
%                             flt = flt & exonUtrDnaFilter;
%                         end
%                     end
%                     trueset = dnavar.locidx(flt);
%                 end
%                 ratio = zeros(nfilter, 2);
%                 
%                 for j = 1:nfilter
%                     rnaflt = VarFilter.CovFilter(rnavar, rnafd, para.RnaGatkFilterByNormCov);
%                     rnaflt = rnaflt & VarFilter.CumValidAfterFilter(rnavar, rnafd, j, para.RnaUseSomatic);
%                     rnaflt = rnaflt & VarFilter.NumReadFilter(rnavar, rnafd, para.RnaFilterByRnaCov, j);
%                     rnaflt = rnaflt & normalRnaFilter;
%                     
%                     if para.RnaUseSomatic
%                         rnaflt = rnaflt & exonUtrRnaFilter;
%                     end
%                     setdata{4} = rnavar.locidx( rnaflt );
%                     if para.applyDnaFilter
%                         flt = VarFilter.TrueSet(dnavar, strcat(para.caller(1:3), para.aligner), para.DnaTrueSet{i});
%                         for di = 1:3
%                             dnafd = [para.caller{di}, para.aligner];
%                             flt = flt & VarFilter.CovFilter(dnavar, dnafd, para.DnaFilterByRnaCov);
%                             flt = flt & VarFilter.CumValidAfterFilter(dnavar, dnafd, max(j,2), para.DnaUseSomatic);
%                             flt = flt & VarFilter.NumReadFilter(dnavar, dnafd, para.DnaFilterByDnaCov, j);
%                             if para.filterExonUtr
%                                 flt = flt & exonUtrDnaFilter;
%                             end
%                         end
%                         trueset = dnavar.locidx(flt);
%                     end
%                     ratio(j, 1) = length(intersect(setdata{4}, trueset)) ./ length(trueset);
%                     ratio(j, 2) = length(intersect(setdata{4}, trueset)) ./ length(setdata{4});
%                 end
%                 subplot(2,2,i);
%                 [ax,h1,h2] = plotyy(1:nfilter+1, ratio(:,1), 1:nfilter+1, ratio(:,2));
%                 set(h1, 'linestyle', '-', 'marker', 'x');
%                 set(h2, 'linestyle', '-', 'marker', 'x');
%                 set(get(ax(1),'ylabel'), 'string', 'in DNA');
%                 set(get(ax(2),'ylabel'), 'string', 'in RNA-GATK');
%                 set(ax(2),'xtick',[]);
%                 set(ax(1), 'xlim', [0.5, nfilter+1.5]);
%                 set(ax(2), 'xlim', [0.5, nfilter+1.5]);
%                 
%                 th = title(para.DnaTrueSet{i});
%                 xticklabel = strrep(rnavar.filterName(1:nfilter), 'rm', '');
%                 xticklabel = strrep(xticklabel, 'sk', 'rpt');
%                 set(gca, 'xtick', 1:nfilter+1, 'xticklabel', [xticklabel, {'ExonUtr'}]);
%             end
%             if do.DnaSom_RnaGatk{2}
%                 saveas(fig, sprintf(['%s/' do.DnaSom_RnaGatk{4}],para.figdir, do.sample,fltfiletag), 'png');
%             end
%         end
%     end
%     close(fig(2));
%     fig(2) = [];
% end

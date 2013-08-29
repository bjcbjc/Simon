
patientid = '138381'; %'132540'; %
if ~exist('patientData', 'var')
    load(sprintf('data/patient.%s.mat', patientid));
end
if ~exist('GENCODE', 'var')
    load data/GENCODE.mat
end
if ~exist('DESeqData', 'var')
    load data/DESeq.NormalTumor.mat;    
end
if ~exist('CENTROMERE', 'var')
    load CENTROMERE.hg19.mat
end

if ~isfield(DESeqData, 'loc')
    DESeqData.loc = GENOMEFUNC.idQuery(GENCODE, DESeqData.gid, {'chrm', 'start', 'end'});
end

samplelabel = {'normal', 'tumor'};
samplecolor = {'b','r'};
minread = 11;

datatype = {'RNA', 'DNA'};

callerColidx = find(strcmp(patientData.variantAttr_cell_name, 'set'));

caller = {'mutect', 'sniper', 'varscan', 'mutect-sniper', 'mutect-varscan', ...
    'sniper-varscan', 'mutect-sniper-varscan'};
callerColor = [0.8, 0.1, 0.1; 0.1, 0.6, 0.1; 0.1, 0.1, 0.8; 0.8, 0.6, 0.2; 0.8, 0.3, 0.6; 0.2, 0.8, 0.7; 0,0,0];

%columns: normal, tumor
totalReadDNA = full(sum(patientData.ALLCountDNA{1}(:, 2:7),2)); %exclude those in introns
totalReadDNA(:,2) = sum(patientData.ALLCountDNA{2}(:, 2:7),2);
totalReadRNA = full(sum(patientData.ALLCountRNA{1}(:, 2:7),2)); %exclude those in introns
totalReadRNA(:,2) = sum(patientData.ALLCountRNA{2}(:, 2:7),2);

fprintf('total: %d, expressed: %d (N:%d, T:%d)\n', ...
    size(totalReadRNA,1), sum(any(totalReadRNA,2)), sum(totalReadRNA(:,1)~=0), sum(totalReadRNA(:,2)~=0));

distToCen = NaN(size(patientData.loc,1),1);
for i = 1:24
    varidx = patientData.loc(:,1)==i;
    cloc = CENTROMERE.loc(CENTROMERE.loc(:,1)==i,2:3);
    cloc = cloc(:)';
    distToCen(varidx) = min(abs(bsxfun(@minus, cloc, patientData.loc(varidx,2))),[],2);    
    clear varidx cloc
end

figure(1);
set(1, 'position', [200, 40, 1000, 650]);
set(1, 'paperpositionmode', 'auto');

%% pie chart: % of variants have reads in normal and/or tumor
clf
expressed = totalReadRNA > 0;
strtoken = cell(size(expressed));
strtoken(:) = {'-'};
strtoken(expressed(:,1),1) = {'Normal'};
strtoken(expressed(:,2),2) = {'Tumor'};
strtoken = strcat(strtoken(:,1), '+', strtoken(:,2));
strtoken = strrep(strtoken, '-+-', 'No reads');
strtoken = regexprep(strtoken, '^\-\+|\+\-$', '');
[count, pattern] = eleCounts(strtoken);
h = pieplus(count, 'inside', 1);
set(h(1:2:end), 'edgecolor', 'none');
set(h(2:2:end), 'fontsize', 14, 'fontweight', 'bold');
legend(pattern, 'fontsize', 14, 'location', 'EO');
title(sprintf('Patient %s\n%% variants have reads', patientid), 'fontsize', 14);
colormap([0.7 0.7 0.7; 0.3 0.3 0.8; 0.8 0.3 0.8; 0.8 0.3 0.3]);
saveas(1, sprintf('figures/RNA has-read pie, %s.eps', patientid), 'epsc2');

clear h strtoken count pattern

%% distribution of read-depth for variants
clf

colormap jet
subplot(2, 1, 1);
hist(log10(totalReadRNA+1), 30);
xlabel('log10 num reads', 'fontsize', 14);
ylabel('num variants', 'fontsize', 14);
legend(samplelabel, 'fontsize', 14);

subplot(2, 1, 2);
readlimit = [1, 6, 11, 101];
percent = [];
for i = 1:length(readlimit)
    percent = [percent; sum(totalReadRNA >= readlimit(i),1) ./ size(totalReadRNA,1)];
end
bar(percent*100, 'grouped');
set(gca, 'xtick', 1:length(readlimit), 'xticklabel', strcat('>=', cellfun(@num2str, num2cell(readlimit), 'unif', 0)));
set(gca, 'fontsize', 14);
ylabel('% variants', 'fontsize', 14);
legend(samplelabel, 'fontsize', 14);
saveas(1, sprintf('figures/distribution of num of RNA reads, %s.eps', patientid), 'epsc2');

clear readlimit percent

%% plot raw count, normal vs tumor

for dataidx = 1:length(datatype)
    subplot(2,2, (dataidx-1)*2+1 );
    plot(patientData.(['refCount' datatype{dataidx}])(:,1), ...
        patientData.(['refCount' datatype{dataidx}])(:,2), 'x', 'markersize',12);
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('REF read-count in %s, all variants', datatype{dataidx}), 'fontsize',14);
    axis square
    ylim([0 10000])
    xlim(ylim);
    hold on; plot(xlim, xlim, 'k--'); hold off;

    subplot(2,2, dataidx*2);
    plot(patientData.(['altCount' datatype{dataidx}])(:,1,1), ...
        patientData.(['altCount' datatype{dataidx}])(:,2,1), 'x', 'markersize', 12);
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('major ALT read-count in %s, all variants', datatype{dataidx}), 'fontsize', 14);
    axis square
    ylim([0 8000])
    xlim(ylim);
    hold on; plot(xlim, xlim, 'k--'); hold off;
end
saveas(1, sprintf('figures/raw count, all variants, normal vs tumor, %s.eps', patientid), 'epsc2');

%% choose subset of variants
% expressed = any(totalReadRNA>0,2); %expressed in either tumor or normal
expressed = all(totalReadRNA>=minread,2); %expressed in both tumor and normal

multAlt = ~cellfun(@isempty, patientData.alt(:,2));
gatkindel = strcmp(patientData.variantAttr_cell(:, ...
    strcmp(patientData.variantAttr_cell_name, 'set')), 'gatkIndel');
chrMT = patientData.loc(:,1) == 25;

% exclude multiple calls of ALT, indels, and chrMT for now
exclude = multAlt | gatkindel | chrMT;

subsetidx = expressed & ~exclude;

%%
DEgene = DESeqData.padj < 0.01;
DEloc = GENOMEFUNC.isWithinRegion(patientData.loc, DESeqData.loc(DEgene,:)); 
nonDEloc = GENOMEFUNC.isWithinRegion(patientData.loc, DESeqData.loc(DESeqData.padj >= 0.5,:));
[~, deidx] = GENOMEFUNC.isWithinRegion(patientData.loc, DESeqData.loc);
DEpvalAdj = NaN(size(patientData.loc,1),1);
DEpvalAdj(deidx~=0) = DESeqData.padj(deidx(deidx~=0));

DElabel = cell(length(DEloc),1);
DElabel(:) = {'NA'};
DElabel(DEpvalAdj<0.01) = {'padj < 0.01'};
DElabel(DEpvalAdj>=0.5) = {'padj >= 0.5'};
DElabel(DEpvalAdj>=0.01 & DEpvalAdj<0.5) = {'in-between'};
%there are some variants located in overlap of DE and non DE genes

%% calculate DNA, RNA allele-frequency 

%normalized by total reads at the bp
AF.DNA.ref = patientData.refCountDNA ./ totalReadDNA;
AF.DNA.alt = bsxfun(@rdivide, patientData.altCountDNA, totalReadDNA);
AF.RNA.ref = patientData.refCountRNA ./ totalReadRNA;
AF.RNA.alt = bsxfun(@rdivide, patientData.altCountRNA, totalReadRNA);

%normalized by ref+alt)
refAltSum = patientData.refCountDNA + nansum(patientData.altCountDNA, 3); 
alleleRatio.DNA.ref = patientData.refCountDNA ./ refAltSum; %loc x sample
alleleRatio.DNA.alt = bsxfun(@rdivide, patientData.altCountDNA, refAltSum);

refAltSum = patientData.refCountRNA + nansum(patientData.altCountRNA, 3); 
alleleRatio.RNA.ref = patientData.refCountRNA ./ refAltSum; %loc x sample
alleleRatio.RNA.alt = bsxfun(@rdivide, patientData.altCountRNA, refAltSum);
clear dnaAD refAltSum

RAratio.DNA = patientData.refCountDNA ./ (patientData.altCountDNA(:,:,1) + 1); %avoid divided by 0
RAratio.RNA = patientData.refCountRNA ./ (patientData.altCountRNA(:,:,1) + 1); %avoid divided by 0

%% plot raw count, normal vs tumor, color by DEseq-result; subset variant only

for dataidx = 1:length(datatype)
    subplot(2,3, (dataidx-1)*3+1 );
    h = gscatter(log10(patientData.(['refCount' datatype{dataidx}])(subsetidx,1)+1), ...
        log10(patientData.(['refCount' datatype{dataidx}])(subsetidx,2)+1), ...
        DElabel(subsetidx));
    set(h, 'markersize', 8, 'marker', 'o');
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 REF read-count in %s', datatype{dataidx}), 'fontsize',14);
    axis square    
    ylim(xlim);    
    if dataidx == 1
        h = legend('orientation', 'horizontal');
        set(h, 'position', [0.5 0.5 0.05 0.03], 'fontsize', 14); 
    else
        legend('off');
    end
    hold on; plot(xlim, xlim, 'k--'); hold off;    
    
    subplot(2,3, (dataidx-1)*3+2 );
    h = gscatter(log10(patientData.(['altCount' datatype{dataidx}])(subsetidx,1,1)+1), ...
        log10(patientData.(['altCount' datatype{dataidx}])(subsetidx,2,1)+1), ...
        DElabel(subsetidx));
    set(h, 'markersize',8, 'marker', 'o');
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 major ALT read-count in %s', datatype{dataidx}), 'fontsize', 14);
    axis square    
    ylim(xlim);
    hold on; plot(xlim, xlim, 'k--'); hold off;
    legend('off');
    
    subplot(2,3, dataidx*3 );
    h = gscatter(log10(eval(['totalRead' datatype{dataidx} '(subsetidx,1)'])+1), ...
        log10(eval(['totalRead' datatype{dataidx} '(subsetidx,2)'])+1), ...
        DElabel(subsetidx));
    set(h, 'markersize', 8, 'marker', 'o');
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 total read-count in %s', datatype{dataidx}), 'fontsize', 14);
    axis square
    ylim(xlim);
    hold on; plot(xlim, xlim, 'k--'); hold off;
    legend('off');
end
saveas(1, sprintf('figures/raw count, normal vs tumor, with DE-Seq calls %s.eps', patientid), 'epsc2');

%% plot raw count, normal vs tumor, color by DEseq-result, chrm, distance to centromere; only on subset of variants

for dataidx = 1
    subplot(1,3, 1 );
    h = gscatter(log10(eval(['totalRead' datatype{dataidx} '(subsetidx,1)'])+1), ...
        log10(eval(['totalRead' datatype{dataidx} '(subsetidx,2)'])+1), ...
        DElabel(subsetidx));
    set(h, 'markersize', 8, 'marker', 'o');
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 read-count in %s', datatype{dataidx}), 'fontsize',14);
    axis square    
    h = legend;
    set(h, 'position', [0.2, 0.15, 0.1,0.1]);
    ylim(xlim);    
    hold on; plot(xlim, xlim, 'k--'); hold off;
    
    
    subplot(1,3, 2);
    cmap = genColorMap('rkg',22);
    cmap = [0.2, 0.2, 0.9; cmap(1:18,:); 0.9,0.7,0.2; cmap(19:end,:)];
    scatter(log10(eval(['totalRead' datatype{dataidx} '(subsetidx,1)'])+1), ...
        log10(eval(['totalRead' datatype{dataidx} '(subsetidx,2)'])+1), ...
        30, patientData.loc(subsetidx,1));
    colormap(cmap);    
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 read-count in %s', datatype{dataidx}), 'fontsize', 14);
    axis square
    xlim([1, 4.5]);
    ylim(xlim);
    hold on; plot(xlim, xlim, 'k--'); hold off;
    freezeColors;
    h = colorbar('location','southoutside');    
    p = get(h, 'position'); p(2) = p(2) - 0.1;
    set(h, 'position', p);
    cbfreeze(h);        
    text(2, p(2) - 1, 'chromosome');
    
    subplot(1,3, 3);    
    scatter(log10(eval(['totalRead' datatype{dataidx} '(subsetidx,1)'])+1), ...
        log10(eval(['totalRead' datatype{dataidx} '(subsetidx,2)'])+1), ...
        30, log10(distToCen(subsetidx)+1));
    colormap(jet(20));    
    xlabel('normal', 'fontsize', 14);
    ylabel('tumor', 'fontsize', 14);
    title(sprintf('log10 read-count in %s', datatype{dataidx}), 'fontsize', 14);
    axis square
    xlim([1, 4.5]);
    ylim(xlim);
    hold on; plot(xlim, xlim, 'k--'); hold off;
    freezeColors;
    h = colorbar('location','southoutside');   
    p = get(h, 'position'); p(2) = p(2) - 0.1;
    set(h, 'position', p);
    cbfreeze(h);    
    text(2, p(2) - 1, 'distance to centromere');
end
saveas(1, sprintf('figures/RNA raw count, normal vs tumor, association with chrm1,19 and centromere %s.eps', patientid), 'epsc2');



%% plot DNA vs RNA ratio
clf
for sampi = 1:length(patientData.sample)
    subplot(2,2,sampi);
    plot(alleleRatio.DNA.ref(subsetidx,sampi), alleleRatio.RNA.ref(subsetidx,sampi), 'x', 'color', samplecolor{sampi});
    xlabel('DNA', 'fontsize', 14); 
    ylabel('RNA', 'fontsize', 14);
    title(sprintf('ref/(ref+alt), %s', samplelabel{sampi}), 'fontsize', 14);
    
    subplot(2,2,sampi+2);
    plot(alleleRatio.DNA.alt(subsetidx,sampi,1), alleleRatio.RNA.alt(subsetidx,sampi,1), 'x', 'color', samplecolor{sampi});
    xlabel('DNA', 'fontsize', 14);
    ylabel('RNA', 'fontsize', 14);
    title(sprintf('alt/(ref+alt), %s', samplelabel{sampi}), 'fontsize', 14);
end
saveas(1, sprintf('figures/DNA vs RNA ratio, %s.eps', patientid), 'epsc2');

%% plot DNA vs RNA allele-frequency

coloredByNormalALTRnaAF = true;
if coloredByNormalALTRnaAF
    oldcolormap = colormap;
    colormap(redgreencmap(50));
end

clf
pvalColor = [0.1, 0.1, 0.8; 0.7, 0.7, 0.7; 0.1, 0.6, 0.1; 0.8, 0.1, 0.1];
for sampi = 1:length(patientData.sample)
    %colormap(redbluecmap(20))
%     colormap(genColorMap([0.8, 0, 0.12; 0.02, 0.44, 0.7], 5));
    subplot(2,2,sampi);
    %plot(AF.DNA.ref(subsetidx,sampi), AF.RNA.ref(subsetidx,sampi), 'x', 'color', samplecolor{sampi});
    if coloredByNormalALTRnaAF
        scatter(AF.DNA.ref(subsetidx,sampi), AF.RNA.ref(subsetidx,sampi), 30, AF.RNA.alt(subsetidx,1,1));
    else
        h = gscatter(AF.DNA.ref(subsetidx,sampi), AF.RNA.ref(subsetidx,sampi), DElabel(subsetidx));    
        set(h, 'markersize', 6, 'marker', 'o','linewidth',1.2);
        if sampi == 1
            hleg = legend('location', 'NO', 'orientation', 'horizontal');
            set(hleg, 'position', [0.37, 0.96, 0.3, 0.04], 'fontsize', 14);
        else
            legend('off');
        end
    end
    xlabel('DNA', 'fontsize', 14); 
    ylabel('RNA', 'fontsize', 14);
    title(sprintf('REF freq, %s', samplelabel{sampi}), 'fontsize', 14);
    
    
    subplot(2,2,sampi+2);
    if coloredByNormalALTRnaAF
        scatter(AF.DNA.alt(subsetidx,sampi,1), AF.RNA.alt(subsetidx,sampi,1), 30, AF.RNA.alt(subsetidx,1,1));
    else
        h = gscatter(AF.DNA.alt(subsetidx,sampi,1), AF.RNA.alt(subsetidx,sampi,1), DElabel(subsetidx));
        set(h, 'markersize', 6, 'marker', 'o','linewidth',1.2);        
        legend('off');
    end
    xlabel('DNA', 'fontsize', 14);
    ylabel('RNA', 'fontsize', 14);
    title(sprintf('major ALT freq, %s', samplelabel{sampi}), 'fontsize', 14);
    
%     if sampi==2
%         colorbar('location','EastOutside');
%     end
    
end
if coloredByNormalALTRnaAF
    saveas(1, sprintf('figures/DNA vs RNA allele freq colored by normal RNA-ALT freq, %s.eps', patientid), 'epsc2');
else
    saveas(1, sprintf('figures/DNA vs RNA allele freq, %s.eps', patientid), 'epsc2');
end

if coloredByNormalALTRnaAF
    colormap(oldcolormap);
end
clear coloredByNormalALTRnaAF oldcolormap pvalColor

%% plot REF vs ALT for each sample, raw count

clf
for sampi = 1:length(patientData.sample)
    subplot(2,2,sampi);
    plot(patientData.refCountRNA(subsetidx,sampi), ...
        patientData.altCountRNA(subsetidx,sampi), 'x', 'markersize', 14, 'color', samplecolor{sampi});
    xlabel('REF read count', 'fontsize', 14);
    ylabel('ALT read count', 'fontsize', 14);
    title(sprintf('%s, %s', patientid, samplelabel{sampi}), 'fontsize', 14);
    hold on; plot(xlim, xlim, 'k--'); hold off
    
    subplot(2,2,sampi+2);
    plot(log10(patientData.refCountRNA(subsetidx,sampi)), ...
        log10(patientData.altCountRNA(subsetidx,sampi)), 'x', 'markersize', 14, 'color', samplecolor{sampi});
    xlabel('log10 REF read count', 'fontsize', 14);
    ylabel('log10 ALT read count', 'fontsize', 14);
    title(sprintf('%s, %s', patientid, samplelabel{sampi}), 'fontsize', 14);
    hold on; plot(xlim, xlim, 'k--'); hold off
end
% saveas(1, sprintf('figures/rna raw counts, ref vs alt, %s.eps', patientid), 'epsc2');

%% plot REF / ALT for each sample, freq

clf
for sampi = 1:length(patientData.sample)        
    subplot(1,2,sampi);
    plot(AF.RNA.ref(subsetidx,sampi), ...
        AF.RNA.alt(subsetidx,sampi,1), 'x', 'markersize', 14, 'color', samplecolor{sampi});
    xlabel('REF freq', 'fontsize', 14);
    ylabel('major ALT freq', 'fontsize', 14);
    title(sprintf('%s, %s', patientid, samplelabel{sampi}), 'fontsize', 14);
    axis square    
end
saveas(1, sprintf('figures/RNA allele freq, REF vs ALT, %s.eps', patientid), 'epsc2');

%color by caller groups
clf
for sampi = 1:length(patientData.sample)        
    subplot(1,2,sampi);
    
    h = gscatter(AF.RNA.ref(subsetidx,sampi), ...
        AF.RNA.alt(subsetidx,sampi,1), patientData.variantAttr_cell(subsetidx, callerColidx)); 
    [~, calleridx] = ismember(get(h, 'displayname'), caller);
    for hi = 1:length(h)
        set(h(hi), 'color', callerColor(calleridx(hi),:));
    end
    set(h, 'markersize', 12, 'marker', 'o', 'linewidth', 2);
    set(gca, 'fontsize', 14);
    xlabel('REF freq', 'fontsize', 14);
    ylabel('major ALT freq', 'fontsize', 14);
    title(sprintf('%s, %s', patientid, samplelabel{sampi}), 'fontsize', 14);
    axis square
    legend('location', 'NO');
end
saveas(1, sprintf('figures/RNA allele freq, REF vs ALT, colored by callers, %s.eps', patientid), 'epsc2');


%% REF/ALT


clf

for dataidx = 1:length(datatype)
    subplot(2,2, (dataidx-1)*2+1);
    h = gscatter(RAratio.(datatype{dataidx})(subsetidx, 1), ...
        RAratio.(datatype{dataidx})(subsetidx, 2), ...
        patientData.variantAttr_cell(subsetidx, callerColidx));
    set(h, 'marker', 'o', 'markersize', 12, 'linewidth',1);
    [~, calleridx] = ismember(get(h, 'displayname'), caller);
    for hi = 1:length(h)
        set(h(hi), 'color', callerColor(calleridx(hi),:));
    end
    set(gca, 'fontsize', 14);
    xlabel(sprintf('R/A in normal %s',datatype{dataidx}), 'fontsize', 14);
    ylabel(sprintf('R/A in tumor %s',datatype{dataidx}), 'fontsize', 14);
    lim = max(xlim, ylim);
    lim(1) = 0;
    xlim(lim);
    ylim(lim);
    axis square
    if dataidx == 1
        hleg = legend('location', 'NO', 'orientation', 'horizontal');
        set(hleg, 'position', [0.014, 0.95, 0.97, 0.049], 'fontsize', 14);
    else
        legend('off');
    end

    subplot(2,2,dataidx*2);    
    h = gscatter(RAratio.(datatype{dataidx})(subsetidx, 1), ...
        RAratio.(datatype{dataidx})(subsetidx, 2), ...
        patientData.variantAttr_cell(subsetidx, callerColidx));
    [~, calleridx] = ismember(get(h, 'displayname'), caller);
    for hi = 1:length(h)
        set(h(hi), 'color', callerColor(calleridx(hi),:));
    end
    set(h, 'marker', 'o', 'markersize', 12, 'linewidth', 1.2);
    set(gca, 'fontsize', 14);
    xlabel(sprintf('R/A in normal %s',datatype{dataidx}), 'fontsize', 14);
    ylabel(sprintf('R/A in tumor %s',datatype{dataidx}), 'fontsize', 14);
    xlim( [0, 100])
    ylim(xlim);
    axis square
    legend('off');
end
saveas(1, sprintf('figures/REF-ALT ratio, normal vs tumor, %s.eps', patientid), 'epsc2');

%% plot distribution of R/A by each caller group

% clf
% set(1, 'position', [200, 40, 1500, 500], 'visible', 'off', 'paperpositionmode','auto');
% % orient landscape
% for sampi = 1:length(patientData.sample)    
%     callergroups = unique(patientData.variantAttr_cell(subsetidx, callerColidx));    
%     for i = 1:length(callergroups)
%         vi = strcmp(patientData.variantAttr_cell(:, callerColidx), callergroups{i}) & subsetidx;
%         subplot(2, length(callergroups), i+(sampi-1)*length(callergroups))
%         hist(RAratio(vi, sampi), 10);
%         title(sprintf('%s, %s', callergroups{i}, samplelabel{sampi}), 'fontsize', 14);
%         xlabel('R/A ratio', 'fontsize', 12);
%         ylabel('# variants', 'fontsize', 12);
%         lim = xlim;
%         xlim([0 lim(2)]);
%     end    
% end
% % saveas(1, sprintf('figures/distribution of ALT-REF ratio by callers, %s.eps', patientid), 'epsc2');
% set(1, 'visible', 'on');
%% dbsnp vs bad A/R

% clf
% dbsnp = ~cellfun(@isempty, regexp(patientData.variantAttr_cell(:,1),'^rs'));
% subplot(2,2,1);
% hist(RAratio(subsetidx & dbsnp, 1), 50);
% ylim([0 1000]);
% xlabel('ALT/REF', 'fontsize', 14);
% title('DBSNP, normal', 'fontsize', 14);
% subplot(2,2,3);
% hist(RAratio(subsetidx & ~dbsnp, 1), 50);
% ylim([0 1000]);
% xlabel('ALT/REF', 'fontsize', 14);
% title('non-DBSNP, normal', 'fontsize', 14);
% subplot(2,2,2);
% hist(RAratio(subsetidx & dbsnp, 2), 50);
% ylim([0 1000]);
% xlabel('ALT/REF', 'fontsize', 14);
% title('DBSNP, tumor', 'fontsize', 14);
% subplot(2,2,4);
% hist(RAratio(subsetidx & ~dbsnp, 2), 50);
% ylim([0 1000]);
% xlabel('ALT/REF', 'fontsize', 14);
% title('non-DBSNP, tumor', 'fontsize', 14);

% saveas(1, sprintf('figures/hist ALT-REF ratio, separated by dbsnp, %s.eps', patientid), 'epsc2');

%% histogram of major ALT AF for each caller

clf
barcolors = [0.6, 0.6, 1; 1, 0.6, 0.6; 0.3, 0.3, 0.9; 0.9, 0.3, 0.3];
legendlabel = {'DNA-normal', 'DNA-tumor', 'RNA-normal', 'RNA-tumor'};
for i = 1:length(caller)
    idx = subsetidx & strcmp(patientData.variantAttr_cell(:, callerColidx), caller{i});
    subplot(3,3,i);
    if sum(idx) == 1
        tmp = [AF.DNA.alt(idx,:,1) AF.RNA.alt(idx,:,1)];
        for si = 1:4
            hold on;
            plot([tmp(si), tmp(si)], [0, 1], 'color', barcolors(si,:), 'linewidth', 4);            
            hold off;
        end
        ylim([0 2]);
        clear tmp
    else        
        [h, x] = hist([AF.DNA.alt(idx,:,1) AF.RNA.alt(idx,:,1)], 0:0.2:1);
        bar(x, h, 'hist');
        ylim([0, max(h(:))+1]);
        xlim([-0.1 1.1]);
        set(gca, 'box', 'off');
    end
    xlabel('ALT allele freq', 'fontsize', 12);
    ylabel('number of variants', 'fontsize', 12);
    title(caller{i}, 'fontsize', 14);
    if i == length(caller)
        h = legend(legendlabel);
        set(h, 'position',  [0.4273    0.1475    0.1235    0.1337]); 
    end
    if sum(idx) > 1
        h = findobj(gca, 'Type', 'patch');
        h = h(4:-1:1);
        for si = 1:4             
            set(h(si), 'edgecolor', 'none', 'facecolor', barcolors(si,:));
        end
    end
end
saveas(1, sprintf('figures/hist major ALT freq by callers, %s.eps',patientid), 'epsc2');



%% distribution of multiple alt

% altRatioInCall = bsxfun(@rdivide, patientData.altCount, sum(patientData.altCount,3));
% altRatioAll = bsxfun(@rdivide, patientData.altCount, notRefCount);
% 
% for sampi = 1:2
%     subplot(2,2,sampi)
%     vi = ~any(squeeze(isnan(altRatioInCall(:,sampi,:))),2) & multAlt(:,sampi);
%     bar(squeeze(altRatioInCall(vi,sampi,:)), 'stacked');
%     xlim([0 sum(vi)+1]);
%     title(sprintf('%s, ALT-ratio among calls', samplelabel{sampi}));
%     
%     subplot(2,2,sampi+2)
%     vi = ~any(squeeze(isnan(altRatioAll(:,sampi,:))),2) & multAlt(:,sampi);
%     bar(squeeze(altRatioAll(vi,sampi,:)), 'stacked');
%     xlim([0 sum(vi)+1]);
%     title(sprintf('%s, ALT-ratio among all non-ref', samplelabel{sampi}));
% end


%% which ones are indel-expressed



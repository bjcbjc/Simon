
if ~exist('RnaGeneReadCount', 'var')
    load data/RnaGeneReadCount.mat
end
if ~exist('ENSEMBL', 'var')
    load ENSEMBL.mat
end
if ~exist('DESeq_all', 'var')
    load data/DESeq.NormalTumor.all.mat
end


%%
t = regexp(RnaGeneReadCount.ensembl, '(\w+).(\d+)', 'tokens');
t = cellfun(@(x) x{1}, t, 'unif', 0);
t = vertcat(t{:});
if length(unique(t(:,1))) ~= size(t,1)
    error('ensembl id not unique\n');
end
RnaGeneReadCount.ensembl = t(:,1);
[~, i] = ismember(t(:,1), ENSEMBL.gene);
fprintf('%d ENSEMBL IDs not found\n', sum(i==0));
RnaGeneReadCount.ensembl(i==0) = [];
RnaGeneReadCount.count(i==0,:) = [];
RnaGeneReadCount.normcount(i==0,:) = [];

readcountfilter = 10;
minsamp = 10;
vi = sum(RnaGeneReadCount.count>readcountfilter,2)>minsamp-1;
RnaGeneReadCount.ensembl = RnaGeneReadCount.ensembl(vi);
RnaGeneReadCount.count = RnaGeneReadCount.count(vi,:);
RnaGeneReadCount.normcount = RnaGeneReadCount.normcount(vi,:);



gidEz = 5566; %PRKACA entrez
gidEsbl = unique( ENSEMBL.gene( ENSEMBL.entrez == gidEz) );
gidEsbl = gidEsbl{1};
%%
if ~exist('metaTB','var')
    load data/metaTB.mat
end
pairidx = find(ismember(metaTB.normalTumorPair(:,2), RnaGeneReadCount.sample) &ismember(metaTB.normalTumorPair(:,3),RnaGeneReadCount.sample));
[~, sampidx] = ismember( metaTB.normalTumorPair(pairidx,2:3), RnaGeneReadCount.sample);
logratio = log2((RnaGeneReadCount.normcount(:,sampidx(:,2))+1) ./ (RnaGeneReadCount.normcount(:,sampidx(:,1))+1));

%%

%calculate correlation
data = log(RnaGeneReadCount.normcount+1);
valididx = std(data,0,2) >= 0.2;
data = data(valididx,:);

Res.sample = RnaGeneReadCount.sample;
Res.ensembl = RnaGeneReadCount.ensembl(valididx);
Res.gene = GENOMEFUNC.idQuery(ENSEMBL, Res.ensembl, {'entrez'}, 'gene');
gidx = find(~cellfun(@isempty, strfind(Res.ensembl, gidEsbl)));
[Res.r, Res.pvalr] = corr(data',  data(gidx,:)');
[Res.rho, Res.pvalrho] = corr(data',  data(gidx,:)', 'type', 'spearman');


% valididx = std(logratio,0,2) >= 0.2;
logratio = logratio(valididx, :);
logRes.sample = RnaGeneReadCount.sample;
logRes.ensembl = RnaGeneReadCount.ensembl(valididx);
logRes.gene = GENOMEFUNC.idQuery(ENSEMBL, logRes.ensembl, {'entrez'}, 'gene');
gidx = find(~cellfun(@isempty, strfind(logRes.ensembl, gidEsbl)));
[logRes.r, logRes.pvalr] = corr(logratio',  logratio(gidx,:)');
[logRes.rho, logRes.pvalrho] = corr(logratio',  logratio(gidx,:)', 'type', 'spearman');


% save data/Prkaca.analysis.mat Res logRes
%%
% take best corr for each gene
structname = {'Res', 'logRes'};
for nameidx = 1:2
    S = eval(structname{nameidx});
    [S.ugene, ~, uidx] = unique(S.gene);
    ngene = length(S.ugene);
    rmi = isnan(S.ugene);
    S.ur = NaN(ngene, 1);
    S.upvalr = NaN(ngene, 1);
    S.urho = NaN(ngene, 1);
    S.upvalrho = NaN(ngene, 1);
    S.selectidx = NaN(ngene,1);
    for i = 1:ngene
        j = find(uidx == i);
        [S.upvalrho(i), mi] = min(S.pvalrho(j));
        S.urho(i) = S.rho( j(mi) );
        S.ur(i) = S.r( j(mi) );
        S.upvalr(i) = S.pvalr( j(mi) );
        S.selectidx(i) = j(mi);
    end
    S.ugene(rmi) = [];
    S.urho(rmi) = [];
    S.upvalrho(rmi) = [];
    S.ur(rmi) = [];
    S.upvalr(rmi) = [];
    S.selectidx(rmi) = [];
    eval(sprintf('%s=S;',structname{nameidx}));
    clear S
end

%%
%GSEA first
if ~exist('MSigDB', 'var')
    load MSigDB.mat
end
[~, vi, vj] = intersect(logRes.ugene, double(MSigDB.c2_cp.gene));
tic; 
[pvalue, qvalue, score,cutoff,num_annotated,num_non_annotated] = ...
    GSEA(logRes.ur(vi), MSigDB.c2_cp.mtx(vj,:), @abs);
toc;
c2cp.p = pvalue;
c2cp.q = qvalue;
c2cp.score = score;
c2cp.cutoff = cutoff;
c2cp.num_annotated = num_annotated;
c2cp.num_non_annotated = num_non_annotated;

%%
[~, vi, vj] = intersect(logRes.ugene, double(MSigDB.c6_all.gene));
tic; 
[pvalue, qvalue, score,cutoff,num_annotated,num_non_annotated] = ...
    GSEA(logRes.ur(vi), MSigDB.c6_all.mtx(vj,:), @abs);
toc;

%%
if ~exist('NCBI', 'var')
    load NCBI
end
exclude = all(RnaGeneReadCount.count<5,2) | isinf(DESeq_all.log2FoldChange) | isnan(DESeq_all.gene);
%get unique genes by padj from DESeq
valididx = find(~exclude);
[ugene, ~, uidx] = unique(DESeq_all.gene(~exclude));
foldchange = NaN(length(ugene),1);
for i = 1:length(foldchange)
    subsetidx = valididx(uidx==i);
    [~, mi] = min(DESeq_all.padj( subsetidx ));
    foldchange(i) = DESeq_all.log2FoldChange( subsetidx(mi) );
end
[~, si] = sort(foldchange);
table = GENOMEFUNC.entrez2Name(NCBI, ugene(si));
table(:,2) = numarray2strarray(foldchange(si));
% tabwrite('DESeqlog2foldchange.filtered5.bestpadj.symbol.rnk', table);

table(:,1) = numarray2strarray(ugene(si));
tabwrite('DESeqlog2foldchange.filtered5.bestpadj.entrez.rnk', table);

%%
exclude = all(RnaGeneReadCount.count<10,2) | isnan(RnaGeneReadCount.gene);
%get unique genes by padj from DESeq
valididx = find(~exclude);
[ugene, ~, uidx] = unique(RnaGeneReadCount.gene(~exclude));
foldchange = NaN(length(ugene),length(RnaGeneReadCount.logratio_sample));
for i = 1:length(foldchange)
    subsetidx = valididx(uidx==i);
    [~, mi] = min(DESeq_all.padj( subsetidx ));
    foldchange(i,:) = RnaGeneReadCount.logratio( subsetidx(mi),: );
end

dlmwriteplus('tables/logfoldchange.entrezgene.txt', foldchange, numarray2strarray(ugene), ['gene'; RnaGeneReadCount.logratio_sample]);

%%
if ~exist('MSigDB', 'var')
    load MSigDB
end
if ~exist('RnaGeneReadCount', 'var')
    load data/RnaGeneReadCount.mat    
end
if ~isfield(RnaGeneReadCount, 'logratio')
    if ~exist('metaTB','var')
        load data/metaTB.mat
    end
    pairidx = find(ismember(metaTB.normalTumorPair(:,2), RnaGeneReadCount.sample) &ismember(metaTB.normalTumorPair(:,3),RnaGeneReadCount.sample));
    [~, sampidx] = ismember( metaTB.normalTumorPair(pairidx,2:3), RnaGeneReadCount.sample);
    RnaGeneReadCount.logratio = log2((RnaGeneReadCount.normcount(:,sampidx(:,2))+1) ./ (RnaGeneReadCount.normcount(:,sampidx(:,1))+1));
    RnaGeneReadCount.logratio_sample = metaTB.normalTumorPair(pairidx,3);
    RnaGeneReadCount.logratio_sampledesc = strcat(metaTB.normalTumorPair(pairidx,1), '-', metaTB.normalTumorPair(pairidx,4));
    RnaGeneReadCount.logratio_sampledesc = strrep(RnaGeneReadCount.logratio_sampledesc, 'patient', 'P');
    RnaGeneReadCount.logratio_sampledesc = strrep(RnaGeneReadCount.logratio_sampledesc, 'primary', 'p');
    RnaGeneReadCount.logratio_sampledesc = strrep(RnaGeneReadCount.logratio_sampledesc, 'meta', 'm');
end

myc_annt = {'KIM_MYC_AMPLIFICATION_TARGETS_UP', ...
    'SCHUHMACHER_MYC_TARGETS_UP', ...
    'COLLER_MYC_TARGETS_UP', ...
    'DANG_REGULATED_BY_MYC_UP', ...
    'DANG_MYC_TARGETS_UP', ...
    'ACOSTA_PROLIFERATION_INDEPENDENT_MYC_TARGETS_UP', ...    
    'PID_MYC_ACTIVPATHWAY'};

keys = {{'pka'}, {'camp_','c-amp'}, {'_myc_'}, {'BENPORATH_MYC_TARGETS_WITH_EBOX'}};
labels = {'PKA', 'cAMP', 'MYC-target-up', 'MYC-target-up with E-box'};
excludecat = find(~cellfun(@isempty, strfind(MSigDB.ALL.categoryName, 'c4')));
excludecat = ismember(MSigDB.ALL.category, excludecat);
excludecat = excludecat | (sum(MSigDB.ALL.mtx)>500)';

catidx = false(length(MSigDB.ALL.name), length(keys));

for i = 1:length(keys)
    if i == 3
        catidx(:,i) = catidx(:,i) | ismember(MSigDB.ALL.name, myc_annt);
    else
        for j = 1:length(keys{i})
            catidx(:,i) = catidx(:,i) | ~cellfun(@isempty, strfind(lower(MSigDB.ALL.name), lower(keys{i}{j})));
            
        end         
    end
    catidx(:,i) = catidx(:,i) & ~excludecat;
end

    


%%   
set(gcf, 'paperpositionmode', 'auto');
validgene = all(RnaGeneReadCount.count(:, ~cellfun(@isempty, strfind(RnaGeneReadCount.sample,'T')) )>30,2);
for i = 1:length(keys)
    if i == 4
        gene = double(MSigDB.ALL.gene( any(MSigDB.ALL.mtx(:, catidx(:,i)),2) ...
            & any(MSigDB.ALL.mtx(:,catidx(:,3)),2) ));
    else
        gene = double(MSigDB.ALL.gene(any(MSigDB.ALL.mtx(:, catidx(:,i)),2)));
    end
    expidx = find(ismember(RnaGeneReadCount.gene, gene) & validgene);
    [~, si] = sort(sum(RnaGeneReadCount.logratio(expidx,:),2), 'descend');
    expidx = expidx(si);
    
    nanimagesc(RnaGeneReadCount.logratio(expidx, :), redbluecmap, 'zeroc', true);
    set(gca, 'xtick', 1:length(RnaGeneReadCount.logratio_sample), 'xticklabel', RnaGeneReadCount.logratio_sampledesc);
    if i == 1 || i == 4
        gridforimage('color', [0.7, 0.7, 0.7]);
        set(gca, 'ytick', 1:length(expidx), 'yticklabel', GENOMEFUNC.entrez2Name(NCBI, RnaGeneReadCount.gene(expidx)),'fontsize',12);
    end    
    colorbar
    title(labels{i}, 'fontsize', 14, 'fontweight', 'bold');
    saveas(gcf, sprintf('figures/PRKACA/%s_log2foldchange.png', labels{i}), 'png');
%     waitforbuttonpress;
end

%%
gidx = find(MSigDB.ALL.gene==5566);
catidx = find(MSigDB.ALL.mtx(gidx,:) & ~ismember(MSigDB.ALL.category, [6,7])' & sum(MSigDB.ALL.mtx)<200);
set(gcf, 'paperPositionMode', 'auto');
validgene = all(RnaGeneReadCount.count(:, ~cellfun(@isempty, strfind(RnaGeneReadCount.sample,'T')) )>30,2);
for i = 1:length(catidx)
    gene = double(MSigDB.ALL.gene(any(MSigDB.ALL.mtx(:, catidx(:,i)),2)));
    expidx = find(ismember(RnaGeneReadCount.gene, gene) & validgene);
    [~, si] = sort(sum(RnaGeneReadCount.logratio(expidx,:),2), 'descend');
    expidx = expidx(si);
    
    nanimagesc(RnaGeneReadCount.logratio(expidx, :), redbluecmap, 'zeroc', true);
    set(gca, 'xtick', 1:length(RnaGeneReadCount.logratio_sample), 'xticklabel', RnaGeneReadCount.logratio_sampledesc);
    if length(expidx) < 51
        gridforimage('color', [0.7, 0.7, 0.7]);
        set(gca, 'ytick', 1:length(expidx), 'yticklabel', GENOMEFUNC.entrez2Name(NCBI, RnaGeneReadCount.gene(expidx)),'fontsize',12);
    end    
    colorbar
    title(strrep(MSigDB.ALL.name{catidx(i)}, '_', ' '), 'fontsize', 14, 'fontweight', 'bold');
    saveas(gcf, sprintf('figures/PRKACA/%s.png', MSigDB.ALL.name{catidx(i)}), 'png');
end


%%
group = cell(length(RnaGeneReadCount.gene),1);
group(:) = {'rest'};
group(ismember(RnaGeneReadCount.gene, double(MSigDB.ALL.gene(any(MSigDB.ALL.mtx(:,catidx(:,4)),2))))) = {'MYC-EBOX'};
for i = 1:size(RnaGeneReadCount.logratio,2)
    subplot(4,5,i);
    boxplot(RnaGeneReadCount.logratio(:,i), group)
end
%%
kmres = parseText('tables/kmcluster86.txt','nrowname',0,'ncolname',1,'numeric',true);
kmres.gene = kmres.text(:,1);
kmres.name = kmres.colname(2:end);
kmres.mtx = kmres.text(:,2:end);
kmres = rmfield(kmres, {'colname', 'text'});

exclude = all(RnaGeneReadCount.count<10,2) | isnan(RnaGeneReadCount.gene);
%get unique genes by padj from DESeq
valididx = find(~exclude);
[ugene, ~, uidx] = unique(RnaGeneReadCount.gene(~exclude));
foldchange = NaN(length(ugene),length(RnaGeneReadCount.logratio_sample));
for i = 1:length(foldchange)
    subsetidx = valididx(uidx==i);
    [~, mi] = min(DESeq_all.padj( subsetidx ));
    foldchange(i,:) = RnaGeneReadCount.logratio( subsetidx(mi),: );
end

genename = GENOMEFUNC.entrez2Name(NCBI, ugene);
%%
for ki = 1:10%size(kmres.mtx,2)
    kmname = regexprep(kmres.name{ki}, '[-\s]', '');
    fn = sprintf(['figures/km86_foldchange/' kmname '.png']);
    [~, dataidx] = ismember(kmres.gene(kmres.mtx(:,ki)==1), ugene);
    writeHeatmapImage(foldchange(dataidx,:), fn, ...
        'rowlabel', genename(dataidx), ...
        'collabel', RnaGeneReadCount.logratio_sampledesc, ...
        'clim', [-4, 4], 'colormap', redbluecmap(20), 'title', kmname );
end

%%
clusmean = NaN(size(kmres.mtx,2), size(foldchange,2));
for ki = 1:size(kmres.mtx,2)
    [~, dataidx] = ismember(kmres.gene(kmres.mtx(:,ki)==1), ugene);
    clusmean(ki,:) = trimmean(foldchange(dataidx,:), 10);
end
gidx = find(ugene==5566);
[c, p] = corr(foldchange(gidx, :)', clusmean');

readcountfn = 'data/RnaGeneReadCount.mat';
% if ~exist(readcountfn, 'file')
    t = parseText('RNA/DESeq/Count_matrix.txt', 'nrowname', 1, 'ncolname',1, 'numeric', true);
    RnaGeneReadCount.sample = strrep(t.colname, 'Sample_', '');
    RnaGeneReadCount.count = t.text;
    RnaGeneReadCount.ensembl = t.rowname;
    
    tn = parseText('RNA/DESeq/Normalized_Count_matrix.txt', 'nrowname', 1, 'ncolname',1, 'numeric', true);
    tn.colname = strrep(tn.colname, 'Sample_', '');
    if nnz(~strcmp(t.rowname, tn.rowname))
        [~, ri] = ismember(RnaGeneReadCount.ensembl, tn.rowname);
    else
        ri = 1:length(RnaGeneReadCount.ensembl);
    end
    if nnz(strcmp(t.colname, tn.colname))
        [~, ci] = ismember(RnaGeneReadCount.sample, tn.colname);
    else
        ci = 1:length(RnaGeneReadCount.sample);
    end
    RnaGeneReadCount.normcount = tn.text(ri,ci);
    clear t    
    save(readcountfn, 'RnaGeneReadCount');
% else
%     load(readcountfn);
% end
%%
if ~exist('patienttb', 'var')
    patienttb = loadStructData('data/metaTB.mat');
    patienttb = patienttb.normalTumorPair;
end

if ~exist('RnaVarLoc', 'var')
    load data/RnaVarLoc.mat
end
%%
% read counts vs number of variants
callerparaset = fieldnames(RnaVarLoc);
callerparaset(~cellfun(@isempty, strfind(callerparaset, 'mutect'))) = [];
rmi = [];
for i = 1:length(callerparaset)
    if any(cellfun(@isempty, RnaVarLoc.(callerparaset{i}).lockey(:,2)))
        rmi = [rmi; i];
    end
end
callerparaset(rmi) = [];
npair = size(patienttb, 1);
nvariant.count = zeros(npair, 0);
nvariant.setname = {};
nfilteredvariant.count = zeros(npair, 0);
nfilteredvariant.setname = {};
for i = 1:length(callerparaset)
    nvariant.count(:, end+1) = cellfun(@(x) size(x,1), RnaVarLoc.(callerparaset{i}).lockey(:,2));
    nvariant.setname{end+1} = callerparaset{i};
    if ~isempty(strfind(callerparaset{i}, 'sniper'))     
        nfilteredvariant.count(:, end+1) = cellfun(@(x) sum(x(:,2)==2), RnaVarLoc.(callerparaset{i}).SS(:,2));
        nfilteredvariant.setname{end+1} = callerparaset{i};
    end
    if ~isempty(strfind(callerparaset{i}, 'varscan'))
        nfilteredvariant.count(:, end+1) = cellfun(@(x) sum(x==2), RnaVarLoc.(callerparaset{i}).SS(:,2));
        nfilteredvariant.setname{end+1} = callerparaset{i};
    end
end
%%
[~, sampidx] = ismember(patienttb(:,3), RnaGeneReadCount.sample);
[~, normsampidx] = ismember(patienttb(:,2), RnaGeneReadCount.sample);
% countstat = [mean( (RnaGeneReadCount.count(:,sampidx)+1)./ (RnaGeneReadCount.count(:,normsampidx)+1) ); ...
%     mean( (log10(RnaGeneReadCount.count(:,sampidx)+1)+1) ./ (log10(RnaGeneReadCount.count(:,normsampidx)+1)+1) ) ];

countstat = [mean(RnaGeneReadCount.count(:, sampidx)); ...    
    mean(log10(RnaGeneReadCount.count(:, sampidx)+1)) ];
normcountstat = [mean(RnaGeneReadCount.count(:, normsampidx)); ...    
    mean(log10(RnaGeneReadCount.count(:, normsampidx)+1)) ];
countstat = countstat ./ normcountstat;
clf
Ylabel = {'#var', '#filter var'};
Xlabel = {'mean T/N',  'mean log10(T/N)'};
for colidx = 1:2
    [~, si] = sort(countstat(colidx,:));
    for rowidx = 1:2    
        subplot(2,2, (rowidx-1)*2+colidx)        
        if rowidx == 1
            plot(countstat(colidx,si), nvariant.count(si,:), 'x-', 'linewidth', 1.2);
            legend(strrep(nvariant.setname, '_', ' '), 'location', 'NO', 'orientation','horizontal','fontsize',12);
        else
            plot(countstat(colidx,si), nfilteredvariant.count(si,:), 'x-', 'linewidth', 1.2);
            legend(strrep(nfilteredvariant.setname, '_', ' '), 'location', 'NO', 'orientation','horizontal','fontsize',12);
        end
        xlabel(Xlabel{colidx}, 'fontsize', 12);
        ylabel(Ylabel{rowidx}, 'fontsize', 12);
        xlim([min(countstat(colidx,si)) - 0.01, max(countstat(colidx,si)) + 0.01]);        
    end
end
set(gcf, 'position', [1400, 1200, 1000, 8000]);
% plot2svg('figures/rnavar/readcount vs variant, no mdup.svg', gcf);

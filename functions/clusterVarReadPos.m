%function clusterVarReadPos(tumorsample, Q, figdir)

Q = 0;

fn = '132540-10N-132540-1T.mutect-sniper-varscan.union.vcf_132540-1-T.dedup.realigned.recal.bam.varreadpos';
fn = sprintf('DNA/bam/varpos/pairmappedQ%d/%s', Q, fn);
tumorsample = regexp(fn, 'vcf_([\w\-]+).', 'tokens', 'once');
tumorsample = tumorsample{1};

varpos = parseText(fn,'nrowname',1, 'ncolname',0, 'numeric', true);

varsupportpcnt = bsxfun(@rdivide, varpos.text(:, 3:end), sum(varpos.text(:,3:end),2));
hassupport = ~any(isnan(varsupportpcnt),2);

kidx = NaN(size(hassupport));
warning off; 
kidx(hassupport,:) = kmeans(varsupportpcnt(hassupport,:), 60, 'emptyaction', 'drop'); 
warning on;

vcfData = loadStructData(sprintf('data/vcf-%s.mat', strtok(tumorsample, '-')));

vcf_varloc = [numericchrm(vcfData.variantAttr_cell(:,1)) vcfData.variantAttr_mtx(:,1)];
vcf_varindex = gloc2index(vcf_varloc);

varreadpos_loc = [numericchrm(varpos.rowname), varpos.text(:,1)];
varreadpos_index = gloc2index(varreadpos_loc);
[~, vcfidx] = ismember(varreadpos_index, vcf_varindex);

ssc = NaN(length(varreadpos_index), 1);
ssc(vcfidx~=0, :) = vcfData.formatData.Sniper_SSC(vcfidx(vcfidx~=0), 1);

% h = figure('position', [0 0 1500, 1200]);
% set(h, 'paperpositionmode', 'auto');
% scvi = find(~isnan(ssc));
% i = intersect(scvi, find(hassupport));
% boxplot(ssc(i), kidx(i))



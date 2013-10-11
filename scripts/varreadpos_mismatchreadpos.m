
fns = listfilename('RNA/bam/varpos/pairmappedQ0/*.varreadpos');
mismatchdatadir = 'RNA/bam/mismatchcount/first_batch/';
tumorname = regexp(fns, 'vcf_([\w\-]+).mdup', 'tokens');
tumorname = cellfun(@(x) x{1}, tumorname, 'unif', 0);
tumorname = vertcat(tumorname{:});

figure(1);
set(1, 'paperpositionmode', 'auto');
readcat = {'first-fwd-M', 'first-rev-M', 'second-fwd-M', 'second-rev-M'};
for ti = 1:length(tumorname)
    varpos = parseText(['RNA/bam/varpos/pairmappedQ0/' fns{ti}],'nrowname',1, 'ncolname',0, 'numeric', true);    
    [mismatch, nread] = plotAlignmentMismatchProfile( ...
        [mismatchdatadir tumorname{ti} '_Aligned.out.WithReadGroup.sorted.bam.mismatch.stat'], [], false);
    readlength = (size(varpos.text,2)-2)/4;
    varsupport = reshape(nansum( bsxfun(@rdivide, varpos.text(:, 3:end), sum(varpos.text(:,3:end),2)) ), readlength, 4);
%     varsupport = reshape(sum( varpos.text(:, 3:end) ), readlength, 4);
    [~, i] = ismember(readcat, mismatch.rowname);
    
    for j = 1:4
        subplot(2,2,j)
        [~, h1, h2] = plotyy(1:readlength, mismatch.text(i(j), 2:end), ...
            1:readlength, varsupport(:,j));
        set(h1, 'marker', 'o');
        set(h2, 'marker', 'o');
        title(readcat{j});
        legend({'mismatch', 'var support'});
    end
    saveas(1, sprintf('figures/mismatchcount/rna_star/first_batch/mismatch and var support, %s, %s.png',strtok(fns{ti}, '_'),tumorname{ti}), 'png');
end
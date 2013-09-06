
fns = listfilename('data/rnavar*var*indel*.mat');
for i = 1:length(fns)
    v_indel = loadStructData(['data/' fns{i}]);
    v_snv = loadStructData(['data/' strrep(fns{i}, 'indel', 'snv')]);
    if any(~strcmp(v_indel.sample, v_snv.sample))
        fprintf('%s, not the same samples\n', fns{i});
        continue;
    end
    lockey1 = v_indel.lockey();
    lockey2 = v_snv.lockey();
    if any(ismember(lockey1, lockey2))
        fprintf('%s, positions overlap between indels and snvs\n', fns{i});
        continue;
    end
    if any(~strcmp(v_indel.attrName_cell, v_snv.attrName_cell))
        fprintf('%s, cell name not the same\n', fns{i});
        continue;
    end
    if any(~strcmp(v_indel.attrName_mtx, v_snv.attrName_mtx))
        fprintf('%s, mtx name not the same\n', fns{i});
        continue;
    end
    fmt1 = fieldnames(v_indel.formatData);
    fmt2 = fieldnames(v_snv.formatData);
    if length(fmt1) ~= length(intersect(fmt1, fmt2))
        fprintf('%s, format data not the same\n', fns{i});
    end
    
    v_snv.variantAttr_cell = [v_snv.variantAttr_cell; v_indel.variantAttr_cell];
    v_snv.variantAttr_mtx = [v_snv.variantAttr_mtx; v_indel.variantAttr_mtx];
    for fmtidx = 1:length(fmt1)
        v_snv.formatData.(fmt1{fmtidx}) = [v_snv.formatData.(fmt1{fmtidx}); v_indel.formatData.(fmt1{fmtidx})];
    end
    v_snv.filename = {v_snv.filename; v_indel.filename};
    vcfData = v_snv;
    save(sprintf('data/%s', strrep(fns{i}, '.indel.', '.')), 'vcfData');
    clear vcfData v_indel v_snv lockey*
    fprintf('%d files done\n', i);
end


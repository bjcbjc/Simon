
vi = regexp(vcfdata.variantAttr_cell(:,1), '^([0-9]{1,2}|[XYMT]{1,2})');
vi = ~cellfun(@isempty, vi);

loc = vcfdata.variantAttr_cell(vi,1);
loc(strcmp(loc, 'X')) = {'23'};
loc(strcmp(loc, 'Y')) = {'24'};
loc(strcmp(loc, 'MT')) = {'25'};
loc = str2double(loc);
loc(:,2) = vcfdata.variantAttr_mtx(vi,1);

vcfdata.variantAttr_cell(~vi,:) = [];
vcfdata.variantAttr_mtx(~vi,:) = [];
vcfdata.validCalls(~vi,:) = [];


allele2vcfidx = zeros(size(aCountData.loc,1),1);
for chrm = 1:25
    asubidx = find(aCountData.loc(:,1)==chrm);
    vsubidx = find(loc(:,1)==chrm);
    [~, posidx] = ismember(aCountData.loc(asubidx,2), loc(vsubidx,2));
    allele2vcfidx( asubidx( posidx~=0) ) = vsubidx(posidx(posidx~=0));
end




vcfref = cellfun(@(x) x(1), vcfdata.variantAttr_cell(:, strcmp(vcfdata.attrName_cell(:,1), 'REF')));
vcfalt = cellfun(@(x) x(1), vcfdata.variantAttr_cell(:, strcmp(vcfdata.attrName_cell(:,1), 'ALT')));

tmp = aCountData.ref;
tmp(:,2) = cellstr(vcfref(allele2vcfidx));
tmp(:,3) = cellstr(vcfalt(allele2vcfidx));
tmp = upper(tmp);
vcfref_incorrect = find(~strcmp(tmp(:,1), tmp(:,2)));


inconsttable = [num2cell(aCountData.loc(vcfref_incorrect,:)), ...
    aCountData.ref(vcfref_incorrect), ...
    vcfdata.variantAttr_cell(allele2vcfidx(vcfref_incorrect), ismember(vcfdata.attrName_cell(:,1), {'REF','ALT'}))];



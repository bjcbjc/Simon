

fns = listfilename('RNA/138381.dna.rna.union.*.rnacount');
fntoken = strcat('loc.normal.', {'rnacount', 'MQ11.rnacount', 'excldup.rnacount', 'excldup.MQ11.rnacount'});
fdname = {'MQ00Dup1','MQ11Dup1', 'MQ00Dup0', 'MQ11Dup0'};


for i = 1:length(fdname)
    fn = fns(~cellfun(@isempty, strfind(fns, fntoken{i})));
    if length(fn) > 1
        error('should only map to one file name');
    end
    
    fn = fn{1};
    data = AlleleCountData.readTableFormatOutput(['RNA/' fn]);
    data = AlleleCountData.collapseStrand(data);
    [data.altcount, data.alt] = AlleleCountData.maxNonRef(data, false);
    data.indelTotal = AlleleCountData.indelCount(data, false, 'all');
    rnaAlleleCount.([fdname{i}, '_normal']) = data;
    
    fn = strrep(fn, '.normal', '');
    data = AlleleCountData.readTableFormatOutput(['RNA/' fn]);
    data = AlleleCountData.collapseStrand(data);
    [data.altcount, data.alt] = AlleleCountData.maxNonRef(data, false);
    data.indelTotal = AlleleCountData.indelCount(data, false, 'all');
    rnaAlleleCount.(fdname{i}) = data;
end

fns = {'138381.dna.rna.union.gatk.exon.bwa.loc.MQ20.Dup0.fromfltbam.normal.rnacount', ...
    '138381.dna.rna.union.gatk.exon.bwa.loc.MQ20.Dup0.fromfltbam.rnacount'};
for i = 1:2
    data = AlleleCountData.readTableFormatOutput(['RNA/' fns{i}]);
    data = AlleleCountData.collapseStrand(data);
    [data.altcount, data.alt] = AlleleCountData.maxNonRef(data, false);
    data.indelTotal = AlleleCountData.indelCount(data, false, 'all');
    if ~isempty(strfind(fns{i}, 'normal'))
        rnaAlleleCount.fromfltbam_normal = data;
    else
        rnaAlleleCount.fromfltbam = data;
    end
end

save data/rnaAlleleCount.138381.mat rnaAlleleCount

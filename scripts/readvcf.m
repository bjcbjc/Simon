
folder = 'DNA/vcf/';
fns = dir([folder '*.vcf']);

for i = 1:length(fns)
    pid = regexp(fns(i).name, '^(\w+)-', 'tokens', 'once');
    vcfData = VCF([folder fns(i).name], true, '/formatdatatmp/'); 
    samplegroup = cell(1,2);
    groupname = cell(1,2);
    for grpi = 1:2
        samplegroup{grpi} = vcfData.sample( (grpi-1)*4+1 : grpi*4 );
        gn = regexp(vcfData.sample{(grpi-1)*4+1}, '\w+\-(\w+)\-([TN])', 'tokens');
        groupname{grpi} = cellarray2str(gn{1}, '');
    end
    vcfData.consolidateFormat(samplegroup, groupname);
    save(sprintf('data/vcf-%s.new.mat',pid{1}),  'vcfData');
end

%%

folder = 'DNA/vcf/new/';
fns = dir([folder '*snpEff.kegg.oncSig.goBP.goCC.goMF.cm.cgn.mir.tft.vcf']);

for i = 1:length(fns)
    pid = regexp(fns(i).name, '^(\w+)-', 'tokens', 'once');
    vcfData = VCF([folder fns(i).name], true, '/formatdatatmp/'); 
    samplegroup = cell(1,2);
    groupname = cell(1,2);
    for grpi = 1:2
        samplegroup{grpi} = vcfData.sample( (grpi-1)*4+1 : grpi*4 );
        gn = regexp(vcfData.sample{(grpi-1)*4+1}, '\w+\-(\w+)\-([TN])', 'tokens');
        groupname{grpi} = cellarray2str(gn{1}, '');
    end
    vcfData.consolidateFormat(samplegroup, groupname);
    save(sprintf('data/vcf-%s.annt.mat',pid{1}),  'vcfData');
end
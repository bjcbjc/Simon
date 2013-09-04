
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

%%
f = dir('RNA/varcall/*.vcf');
fns = cell(size(f));
for i = 1:length(fns)
    fns{i} = f(i).name;
end
newfns = fns;
newfns = strrep(newfns, '_Aligned.out.WithReadGroup.sorted', '');
newfns = strrep(newfns, '_Sample_', '_');
newfns = strrep(newfns, 'mutect_', '');
newfns = strrep(newfns, 'varscan_', '');
newfns = strrep(newfns, 'sniper_', '');
newfns = strrep(newfns, '_-J-s0.01.sniper', '.sniper.-J-s0.01');
newfns = strrep(newfns, '.mdup.bam_', '_');
newfns = strrep(newfns, '.bam_', '_');


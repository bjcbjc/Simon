function batchreadvcf(jobidx, njob, vcfdir, outdir, outfnhead)

    f = dir([vcfdir '*.vcf']);
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
    newfns = strrep(newfns, '.vcf', '');
    
    filesize = [f.bytes];
    [~, si] = sort(filesize);    
    fns = fns(si);
    newfns = newfns(si);
        
    runidx = jobDivider(jobidx, njob, length(newfns));
    
    fprintf('processing %d files\n', length(runidx));
    for i = 1:length(runidx)
        if exist([outdir outfnhead '.' newfns{runidx(i)} '.mat'], 'file') == 2
            fprintf('skip %d(th) file, exists\n',i);
            continue;
        end        
        tic; vcfData = VCF([vcfdir fns{runidx(i)}],[],[],true, 5000); toc;
        save(sprintf('%s/%s.%s.mat', outdir, outfnhead, newfns{runidx(i)}), 'vcfData');
        clear vcfData
        fprintf('processed %d files: %s\n', i, fns{runidx(i)});
    end
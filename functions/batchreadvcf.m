function batchreadvcf(jobidx, njob, vcffnkey, vcfdir, outdir, outfnhead)

    if ~iscell(vcffnkey)
        vcffnkey = {vcffnkey};
    end
    fns = listfilename([vcfdir vcffnkey{1}]);
    for i = 2:length(vcffnkey)
        fns = [fns; listfilename([vcfdir vcffnkey{i}])];
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
    newfns = strrep(newfns, '.bam.', '.');
    newfns = strrep(newfns, '.q50_', '_');
    newfns = strrep(newfns, '.vcf', '');
    newfns = strrep(newfns, '.MQrep_', '_');
    newfns = strrep(newfns, '.MQ254_', '_');
    newfns = strrep(newfns, '.mdup_', '_');
    
    if length(unique(newfns)) ~=  length(newfns)
        fprintf('invalid new file names\n');
        return
    end
%     filesize = [f.bytes];
%     [~, si] = sort(filesize);    
%     fns = fns(si);
%     newfns = newfns(si);
        
    runidx = jobDivider(jobidx, njob, length(newfns));
    
    fprintf('processing %d files\n', length(runidx));
    for i = 1:length(runidx)
        if exist([outdir outfnhead '.' newfns{runidx(i)} '.mat'], 'file') == 2
            fprintf('skip %d(th) file, exists\n',i);
            continue;
        end        
        %tic; vcfData = VCF([vcfdir fns{runidx(i)}],[],[],true, 5000); toc;
        tic; vcfData = VCF([vcfdir fns{runidx(i)}],[],[],false, [], true); toc;
        if size(vcfData.variantAttr_cell,1) ~= size(vcfData.variantAttr_mtx,1)
            fprintf('skip invalid vcf, %s\n', newfns{runidx(i)});
        else
            save(sprintf('%s/%s.%s.mat', outdir, outfnhead, newfns{runidx(i)}), 'vcfData');
        end
        clear vcfData
        fprintf('processed %d files: %s\n', i, fns{runidx(i)});
    end
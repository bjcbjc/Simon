

patient = {'132540', '138381'; ...
    '10N', '4N'; ...
    '1T', '2T'};

chrmNum = {'X', '23'; 'Y', '24'; 'MT', '25'; 'M', '25'};

para.r = 1;
para.dRNA = 20;
para.dDNA = 200;
para.q = 50;


acountFolder = 'data/';

for patientidx = 1:size(patient, 2)
    tic;
    %load vcf, keep chrm only
    vcfData = loadStructData(sprintf('%svcf-%s.mat',acountFolder, patient{1, patientidx}));
    totalcalls = size(vcfData.variantAttr_cell,1);
    vi = regexp(vcfData.variantAttr_cell(:,1), '^([0-9]{1,2}|[XYMT]{1,2})');
    vi = cellfun(@isempty, vi);
    vcfData.variantAttr_cell(vi,:) = [];
    vcfData.variantAttr_mtx(vi,:) = [];    
    fds = fieldnames(vcfData.formatData);
    for fdidx = 1:length(fds)        
        vcfData.formatData.(fds{fdidx})(vi, :, :) = [];        
    end
    chrmidx = strcmp(vcfData.attrName_cell, 'CHROM');
    for chridx = 1:size(chrmNum,1)
        vcfData.variantAttr_cell(:, chrmidx) = strrep( ...
            vcfData.variantAttr_cell(:, chrmidx), chrmNum{chridx,1}, chrmNum{chridx,2});
    end
    vcflockey = strcat(vcfData.variantAttr_cell(:, chrmidx), '-', ...
        num2cellstr( vcfData.variantAttr_mtx(:, strcmp(vcfData.attrName_mtx, 'POS')) ) );
    
    
    %loop over samples to obtain allele read counts
    nsample = size(patient,1) - 1;
    aDataArray = cell(nsample, 2); %RNA, DNA
    lockey = cell(nsample,2);
    locidx = cell(nsample,2);
    datatype = {'RNA', 'DNA'};
    for sampidx = 1:nsample
        for dataidx = 1:length(datatype)
            matfn = sprintf('%sacount.%s.%s-%s.r%d.d%dk.q%d.mat', acountFolder, ...
                datatype{dataidx}, patient{1, patientidx}, ...
                patient{sampidx+1, patientidx}, para.r, para.(['d' datatype{dataidx}]), para.q);
            aDataArray{sampidx,dataidx} = loadStructData(matfn);
            aDataArray{sampidx,dataidx}.acount = sum(aDataArray{sampidx,dataidx}.acount, 3);
            aDataArray{sampidx,dataidx}.ref = upper(aDataArray{sampidx,dataidx}.ref);
            
            tmp = num2cellstr(aDataArray{sampidx,dataidx}.loc);
            lockey{sampidx,dataidx} = strcat(tmp(:,1), '-', tmp(:,2));            
        end
    end
    
    uloc = union(lockey{1,1}, lockey{2,1});
    
    fds = {'loc', 'numreads', 'ref', 'indels', 'acount'};
    for sampidx = 1:nsample %remove DNA acount entries
        rmi = ~ismember(lockey{sampidx,2}, uloc);
        for fdidx = 1:length(fds)
            aDataArray{sampidx,2}.(fds{fdidx})(rmi,:,:) = [];
        end
        lockey{sampidx,2}(rmi) = [];
    end
    
    nvariant = length(uloc);
    [~, toVcfIdx] = ismember(uloc, vcflockey);
    if any(toVcfIdx == 0)
        error('Some locations are not in VCF\n');
    end
    
    patientData.patient = aDataArray{1}.patient;
    patientData.sample = patient(2:end, patientidx);
    patientData.minRead = aDataArray{1}.minRead;
    patientData.maxDepth = aDataArray{1}.maxDepth;
    patientData.minMapQ = aDataArray{1}.minMapQ;
    patientData.loc = zeros(nvariant, 2);
    patientData.numReadRNA = zeros(nvariant, nsample);    
    patientData.numReadDNA = zeros(nvariant, nsample);
    patientData.ref = cell(nvariant, 1); %take from allele-count data
        
    %get multiple alt    
    multAlts = regexp( ...
        vcfData.variantAttr_cell(toVcfIdx, strcmp(vcfData.attrName_cell, 'ALT')), ...
        ',([ATCGNatcgn])', 'tokens');
    nalt = cellfun(@length, multAlts)+1;    
    patientData.alt = cell(nvariant, max(nalt) );    
    patientData.alt(:,1) = ... %take from vcf; but ALT can have multiple alleles
        upper( cellstr(cellfun(@(x) x(1), ...
        vcfData.variantAttr_cell(toVcfIdx, strcmp(vcfData.attrName_cell, 'ALT')))) );        
    rowidx = find(nalt>1);
    for ri = 1:length(rowidx)        
        patientData.alt(rowidx(ri), 2:nalt(rowidx(ri))) = horzcat(multAlts{rowidx(ri)}{:});
    end    
    for ci = 2:size(patientData.alt,2)
        patientData.alt(nalt==ci-1, ci:end) = {''};
    end
    patientData.alt = upper(patientData.alt);
    
    patientData.refCountRNA = zeros(nvariant, nsample);
    patientData.altCountRNA = zeros(nvariant, nsample, max(nalt));
    patientData.refCountDNA = zeros(nvariant, nsample);
    patientData.altCountDNA = zeros(nvariant, nsample, max(nalt));
    patientData.ntlabels = aDataArray{sampidx}.ntlabels(1,:);
    patientData.ALLCountDNA = cell(nsample, 1); %because want to store sparse matrix
    patientData.ALLCountRNA = cell(nsample, 2); %because want to store sparse matrix
    patientData.indelsRNA = cell(nvariant, nsample);
    patientData.indelsDNA = cell(nvariant, nsample);
    
    nnt = length(patientData.ntlabels);
    for sampidx = 1:nsample
        for dataidx = 1:length(datatype)
            [~, locidx{sampidx,dataidx}] = ismember(lockey{sampidx,dataidx}, uloc);
            
            %loc and ref will overwrite between samples; but they should be the
            %same between samples anway
            if dataidx == 1
                patientData.loc(locidx{sampidx,dataidx},:) = aDataArray{sampidx,dataidx}.loc;
                patientData.ref(locidx{sampidx,dataidx}) = aDataArray{sampidx,dataidx}.ref;
            end
            patientData.(['numRead' datatype{dataidx}])(locidx{sampidx,dataidx}, sampidx) ...
                = aDataArray{sampidx,dataidx}.numreads;
            patientData.(['indels' datatype{dataidx}])(locidx{sampidx,dataidx}, sampidx) ...
                = aDataArray{sampidx,dataidx}.indels;

            nvarSamp = length(aDataArray{sampidx,dataidx}.ref);
            [~, colidx] = ismember(aDataArray{sampidx,dataidx}.ref, patientData.ntlabels);
            patientData.(['refCount' datatype{dataidx}])(locidx{sampidx,dataidx}, sampidx) = ...
                aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], (1:nvarSamp)', colidx)  );

            for altidx = 1:size(patientData.alt,2)        
                [~, colidx] = ismember(patientData.alt(locidx{sampidx,dataidx},altidx), patientData.ntlabels);
                patientData.(['altCount' datatype{dataidx}])(locidx{sampidx,dataidx}(colidx~=0), sampidx, altidx) = ...
                    aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], find(colidx~=0), colidx(colidx~=0)) );
            end

            [rowidx, colidx] = find(aDataArray{sampidx,dataidx}.acount);
            patientData.(['ALLCount' datatype{dataidx}]){sampidx} = sparse( locidx{sampidx,dataidx}(rowidx), colidx, ...
                aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt],rowidx,colidx) ), ...
                nvariant ,nnt);
        end
    end
    
    %what to copy from vcf?
    [~, vcfFormatSampleOrder] = ismember(patient(2:end, patientidx), vcfData.sample);
    keep = ~ismember(vcfData.attrName_cell, {'CHROM', 'REF', 'ALT', 'FILTER'});
    patientData.variantAttr_cell_name = vcfData.attrName_cell(keep);
    patientData.variantAttr_cell = vcfData.variantAttr_cell( toVcfIdx, keep);
    keep = ~ismember(vcfData.attrName_mtx, {'POS'});
    patientData.variantAttr_mtx_name = vcfData.attrName_mtx(keep);
    patientData.variantAttr_mtx = vcfData.variantAttr_mtx( toVcfIdx, keep);
    fds = fieldnames(vcfData.formatData);
    for fdidx = 1:length(fds)
        patientData.formatData.(fds{fdidx}) = vcfData.formatData.(fds{fdidx})( toVcfIdx, vcfFormatSampleOrder, :);
    end
    patientData.vcfSample = vcfData.sample(vcfFormatSampleOrder);
    patientData.vcfDesc = vcfData.desc;
    
    %sort alt by count
    rowidx = nalt > 1;
    [~, ntreorder] = sort(squeeze(sum(patientData.altCountRNA(rowidx,:,:), 2)), 2, 'descend');
    rowidx = find(rowidx);
    fds = fds(ismember(fds, vcfData.allelicFormatField));
    [~, fdsDescIdx] = ismember(fds, vcfData.desc(:,1));
    altOnly = strcmp(vcfData.desc(fdsDescIdx, 2), 'A');
    altDataIdx = 2:max(nalt)+1;
    infoIdx = ismember(vcfData.desc(:,1), vcfData.attrName_cell) & ...
        strcmp(vcfData.desc(:,2), 'A') & ~strcmpi(vcfData.desc(:,1), '1000genomes.AF');
    [~, infoIdx] = ismember(vcfData.desc(infoIdx,1), patientData.variantAttr_cell_name);
    infoIdx(infoIdx==0) = [];
    for ri = 1:length(rowidx)
        patientData.altCountRNA(rowidx(ri),:,:) = ...
            patientData.altCountRNA(rowidx(ri),:,ntreorder(ri,:));
        patientData.altCountDNA(rowidx(ri),:,:) = ...
            patientData.altCountDNA(rowidx(ri),:,ntreorder(ri,:));
        patientData.alt(rowidx(ri),:) = patientData.alt(rowidx(ri), ntreorder(ri,:));        
        %for any field related to alt alleles, reorder the data
        for fdidx = 1:length(fds)
            if altOnly(fdidx)
                patientData.formatData.(fds{fdidx})(rowidx(ri),:,:) = ...
                    patientData.formatData.(fds{fdidx})(rowidx(ri), :, ntreorder(ri,:));
            else
                patientData.formatData.(fds{fdidx})(rowidx(ri),:,altDataIdx) = ...
                    patientData.formatData.(fds{fdidx})(rowidx(ri), :, altDataIdx(ntreorder(ri,:)));
            end
        end
        for infoii = 1:length(infoIdx)
            patientData.variantAttr_cell{rowidx(ri), infoIdx(infoii)} = ...
                patientData.variantAttr_cell{rowidx(ri), infoIdx(infoii)}( ntreorder( ri, 1:nalt(rowidx(ri)) ) );
        end
    end
        

    
    save(sprintf('data/patient.%s.q50.mat',patientData.patient), 'patientData');
    if patientidx < size(patient, 2)
        clear patientData
    end
    toc;
end
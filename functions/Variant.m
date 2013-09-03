classdef Variant < dataset
    
%     properties (Constant=true)
%         annotation = {'EFF', 'GO', 'CancerGeneNeighbor', 'CancerModule', ...
%             'OncogenicSignature', 'KEGG', 'Target'};
%     end
    
    properties        
        sample;
        minRead;
        maxDepth;
        minMapQ;
        
        ntlabel;
        vcfSample;
        vcfDesc;
                
        varSnpEff;
        varAnnt;
        
        ALLCountDNA;
        ALLCountRNA;
        
        strkeyDB;
        
        %preset of variants 
        presetidx;
        minReadDNA = 11;
        minReadRNA = 11;
        maxReadDNA = 1000;
        
        VarType;
    end
    
    properties (Access = private)
        %logical index for current selected variants        
        subsetdesc;
        subsetidx;
    end
    
    methods
        function obj = Variant(patientSampleInfo, DESeqData, CENTROMERE, varargin)                        
            para.r = 1;
            para.dRNA = 20;
            para.dDNA = 200;
            para.q = 0;
            para.aCountDir = 'data/';
            
            para = assignpara(para, varargin{:});
                        
            if ischar(CENTROMERE) && ~isempty(CENTROMERE)
                CENTROMERE = loadStructData(CENTROMERE);
            elseif ~isstruct(CENTROMERE)
                CENTROMERE = [];
            end
            if ischar(DESeqData) && ~isempty(DESeqData)
                DESeqData = loadStructData(DESeqData);
            elseif ~isstruct(DESeqData)
                DESeqData = [];
            end
            
            if ~isfield(DESeqData, 'loc')
                GENCODE = loadStructData('data/GENCODE.mat');
                DESeqData.loc = GENOMEFUNC.idQuery(GENCODE, DESeqData.gid, {'chrm', 'start', 'end'});
                clear GENCODE
            end
            
            %load vcf, keep chrm only
            [vcfData, vcflockey] = Variant.loadVcf(patientSampleInfo, para);
                        
            nsample = size(patientSampleInfo,1) - 1;
            %allele read counts
            [aDataArray, lockey, uloc] = Variant.loadAcountData(patientSampleInfo, para, vcflockey);
            
            nvariant = length(uloc);
            [~, toVcfIdx] = ismember(uloc, vcflockey);
            if any(toVcfIdx == 0)
                error('Some locations are not in VCF\n');
            end
            

            altColIdx = strcmp(vcfData.attrName_cell, 'ALT');
            %get multiple alt
            multAlts = regexp( ...
                vcfData.variantAttr_cell(toVcfIdx, altColIdx), ...
                ',([ATCGNatcgn])', 'tokens');
            nalt = cellfun(@length, multAlts)+1;
            maxnalt = max(nalt);
            [~, vcfFormatSampleOrder] = ismember(patientSampleInfo(2:end), vcfData.sample);
            keepcell = ~ismember(vcfData.attrName_cell, {'CHROM', 'REF', 'ALT', 'FILTER'});
            keepmtx = ~ismember(vcfData.attrName_mtx, {'POS'});
            
            fds = fieldnames(vcfData.formatData);
            existfdidx = find(ismember(fds, [vcfData.attrName_cell; vcfData.attrName_mtx]));
            nfds = length(fds);
            for fdidx = 1:length(existfdidx)
                newfd = ['fmt_' fds{existfdidx(fdidx)}];
                vcfData.formatData.(newfd) = vcfData.formatData.(fds{existfdidx(fdidx)});
                vcfData.formatData = rmfield(vcfData.formatData, fds{existfdidx(fdidx)});
                fds{existfdidx(fdidx)} = newfd;                
            end
            for fdidx = 1:nfds
                formatData.(fds{fdidx}) = vcfData.formatData.(fds{fdidx})( toVcfIdx, vcfFormatSampleOrder, :);                
            end
            
            formatDataCell = cell(1, nfds*2);
            formatDataCell(1:2:end) = struct2cell(formatData);
            formatDataCell(2:2:end) = fieldnames(formatData);
            formatDataForDS = arrayfun(@(fdidx) formatDataCell(fdidx:fdidx+1), 1:2:nfds*2, 'unif', 0);
            
            %allocate for speed
            obj@dataset({zeros(nvariant,2), 'loc'}, {cell(nvariant, 2), 'ref', 'alt_original'}, ...
                {cell(nvariant, maxnalt), 'alt'}, {zeros(nvariant, nsample, maxnalt), 'altCountRNA'}, ...
                {zeros(nvariant, nsample, maxnalt), 'altCountDNA'}, ...
                {zeros(nvariant, nsample), 'refCountRNA'}, {zeros(nvariant, nsample), 'refCountDNA'}, ...
                {repmat({''}, nvariant, nsample), 'indelRNA'}, {repmat({''}, nvariant, nsample), 'indelDNA'}, ...
                {vcfData.variantAttr_cell( toVcfIdx, keepcell), vcfData.attrName_cell{keepcell}}, ...
                {vcfData.variantAttr_mtx( toVcfIdx, keepmtx), vcfData.attrName_mtx{keepmtx}}, ...
                {zeros(nvariant, nsample), 'totalReadDNA'}, {zeros(nvariant, nsample), 'totalReadRNA'}, ...
                {zeros(nvariant, nsample), 'refDnaAF'}, {zeros(nvariant, nsample), 'refRnaAF'}, ...
                {zeros(nvariant, nsample, maxnalt), 'altDnaAF'}, {zeros(nvariant, nsample,maxnalt), 'altRnaAF'}, ...
                {zeros(nvariant, nsample), 'RAratioDNA'}, {zeros(nvariant, nsample), 'RAratioRNA'}, ...
                {NaN(nvariant, 2), 'distToCen', 'DEpvalAdj'}, ...
                formatDataForDS{:});
                        
            
            obj.sample = patientSampleInfo(2:end);
            obj.minRead = aDataArray{1}.minRead;
            obj.maxDepth = aDataArray{1}.maxDepth;
            obj.minMapQ = aDataArray{1}.minMapQ;
            obj.ntlabel = aDataArray{1}.ntlabels(1,:);
            obj.ALLCountDNA = cell(nsample, 1); %because want to store sparse matrix
            obj.ALLCountRNA = cell(nsample, 2); %because want to store sparse matrix                        
            obj.vcfSample = vcfData.sample(vcfFormatSampleOrder);
            obj.vcfDesc = vcfData.desc;
                        
            %fill in data
            obj = subsasgn(obj, substruct('.', 'alt', '()', {':', 1}),  ... %take from vcf; but ALT can have multiple alleles
                upper( cellstr(cellfun(@(x) x(1), ...
                vcfData.variantAttr_cell(toVcfIdx, altColIdx))) ) );
            rowidx = find(nalt>1);
            for ri = 1:length(rowidx)
                obj = subsasgn(obj, substruct('.','alt','()',{rowidx(ri), 2:nalt(rowidx(ri))}), horzcat(multAlts{rowidx(ri)}{:}) );
            end
            for ci = 2:maxnalt
                obj = subsasgn(obj, substruct('.', 'alt', '()', {nalt==ci-1, ci:maxnalt}), {''});
            end
            obj = subsasgn(obj, substruct('.', 'alt'), upper(subsref(obj, substruct('.', 'alt'))));
            obj = subsasgn(obj, substruct('.', 'alt_original'), vcfData.variantAttr_cell(toVcfIdx, altColIdx));
                        
            
            locidx = cell(nsample,2);
            nnt = length(obj.ntlabel);
            datatype = {'RNA', 'DNA'};
            for sampidx = 1:nsample
                for dataidx = 1:length(datatype)
                    [~, locidx{sampidx,dataidx}] = ismember(lockey{sampidx,dataidx}, uloc);
                    
                    %loc and ref will overwrite between samples; but they should be the
                    %same between samples anway
                    if dataidx == 1
                        obj = subsasgn(obj, substruct('.','loc','()', {locidx{sampidx,dataidx},':'}), aDataArray{sampidx,dataidx}.loc);
                        obj = subsasgn(obj, substruct('.', 'ref', '()', {locidx{sampidx,dataidx}}), aDataArray{sampidx,dataidx}.ref);
                    end                    
                    obj = subsasgn(obj, substruct('.', ['indels' datatype{dataidx}], '()', {locidx{sampidx,dataidx}, sampidx}),  ...
                        aDataArray{sampidx,dataidx}.indels );
                    
                    nvarSamp = length(aDataArray{sampidx,dataidx}.ref);
                    [~, colidx] = ismember(aDataArray{sampidx,dataidx}.ref, obj.ntlabel);
                    obj = subsasgn(obj, substruct('.', ['refCount' datatype{dataidx}], '()', {locidx{sampidx,dataidx}, sampidx}), ...
                        aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], (1:nvarSamp)', colidx)  ) );
                    
                    for altidx = 1:maxnalt
                        [~, colidx] = ismember(subsref(obj, substruct('.', 'alt', '()', {locidx{sampidx,dataidx},altidx})), obj.ntlabel);
                        obj = subsasgn(obj, substruct('.', ['altCount' datatype{dataidx}], '()', {locidx{sampidx,dataidx}(colidx~=0), sampidx, altidx}), ...
                            aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], find(colidx~=0), colidx(colidx~=0)) ) );
                    end
                    
                    [rowidx, colidx] = find(aDataArray{sampidx,dataidx}.acount);
                    obj.(['ALLCount' datatype{dataidx}]){sampidx} = sparse( locidx{sampidx,dataidx}(rowidx), colidx, ...
                        aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt],rowidx,colidx) ), ...
                        nvariant ,nnt);
                end
            end
            
            for i = 1:nsample
                obj = subsasgn(obj, substruct('.', 'totalReadDNA', '()', {':',i}), full(sum(obj.ALLCountDNA{i}(:, 2:7),2)) ); %exclude those in introns
                obj = subsasgn(obj, substruct('.', 'totalReadRNA', '()', {':',i}), full(sum(obj.ALLCountRNA{i}(:, 2:7),2)) ); %exclude those in introns
            end
            % calculate DNA, RNA allele-frequency
            obj = subsasgn(obj, substruct('.', 'refDnaAF','()',{':',':'}), subsref(obj, substruct('.', 'refCountDNA')) ./ subsref(obj, substruct('.', 'totalReadDNA')));
            obj = subsasgn(obj, substruct('.', 'altDnaAF','()',{':',':',':'}), bsxfun(@rdivide, subsref(obj, substruct('.', 'altCountDNA')), subsref(obj, substruct('.', 'totalReadDNA'))));
            obj = subsasgn(obj, substruct('.', 'refRnaAF','()',{':',':'}), subsref(obj, substruct('.', 'refCountRNA')) ./ subsref(obj, substruct('.', 'totalReadRNA')));
            obj = subsasgn(obj, substruct('.', 'altRnaAF','()',{':',':',':'}), bsxfun(@rdivide, subsref(obj, substruct('.', 'altCountRNA')), subsref(obj, substruct('.', 'totalReadRNA'))));
            obj = subsasgn(obj, substruct('.', 'RAratioDNA','()',{':',':'}), subsref(obj, substruct('.', 'refCountDNA')) ./ (subsref(obj, substruct('.', 'altCountDNA',  '()', {':',':',1})) + 1)); %avoid divided by 0
            obj = subsasgn(obj, substruct('.', 'RAratioRNA','()',{':',':'}), subsref(obj, substruct('.', 'refCountRNA')) ./ (subsref(obj, substruct('.', 'altCountRNA',  '()', {':',':',1})) + 1)); %avoid divided by 0
            
            if ~isempty(CENTROMERE)
                %calculate distance to centromere                
                for i = 1:24
                    varidx = obj.selfref('loc', '()', {':',1}) == i;
                    cloc = CENTROMERE.loc(CENTROMERE.loc(:,1)==i,2:3);
                    cloc = cloc(:)';
                    obj = subsasgn(obj, substruct('.', 'distToCen', '()', {varidx}), ...
                        min(abs(bsxfun(@minus, cloc, subsref(obj, substruct('.', 'loc', '()', {varidx,2})))),[],2));
                    clear varidx cloc
                end
            end
            
            if ~isempty(DESeqData)
                %get DESeq data
                [~, deidx] = GENOMEFUNC.isWithinRegion(obj.selfref('loc'), DESeqData.loc);                
                obj = subsasgn(obj, substruct('.', 'DEpvalAdj', '()', {deidx~=0}), DESeqData.padj(deidx(deidx~=0)));
            end
            
            %sort alt by count
            rowidx = nalt > 1;
            [~, ntreorder] = sort(squeeze(sum(subsref(obj, substruct('.', 'altCountRNA', '()', {rowidx,':',':'})), 2)), 2, 'descend');
            rowidx = find(rowidx);
            fds = fds(ismember(fds, vcfData.allelicFormatField));
            [~, fdsDescIdx] = ismember(fds, vcfData.desc(:,1));
            altOnly = strcmp(vcfData.desc(fdsDescIdx, 2), 'A');
            altDataIdx = 2:max(nalt)+1;
            infoIdx = ismember(vcfData.desc(:,1), vcfData.attrName_cell) & ...
                strcmp(vcfData.desc(:,2), 'A') & ~strcmpi(vcfData.desc(:,1), '1000genomes.AF');
            [~, infoIdx] = ismember(vcfData.desc(infoIdx,1), obj.Properties.VarNames);
            infoIdx(infoIdx==0) = [];
            for ri = 1:length(rowidx)
                obj = subsasgn(obj, substruct('.', 'altCountRNA', '()', {rowidx(ri),':',':'}), ...
                    obj.selfref('altCountRNA', '()', {rowidx(ri),':',ntreorder(ri,:)}));
                obj = subsasgn(obj, substruct('.', 'altCountDNA', '()', {rowidx(ri),':',':'}), ...
                    obj.selfref('altCountDNA', '()', {rowidx(ri),':',ntreorder(ri,:)}));
                obj = subsasgn(obj, substruct('.', 'alt', '()', {rowidx(ri),':'}), ...
                    obj.selfref('alt', '()', {rowidx(ri), ntreorder(ri,:)})) ;
                %for any field related to alt alleles, reorder the data
                for fdidx = 1:length(fds)
                    if altOnly(fdidx)
                        obj = subsasgn(obj, substruct('.', fds{fdidx}, '()', {rowidx(ri),':',':'}), ...
                            subsref(obj, substruct('.', fds{fdidx}, '()', {rowidx(ri), ':', ntreorder(ri,:)})) );
                    else
                        obj = subsasgn(obj, substruct('.', fds{fdidx}, '()', {rowidx(ri),':',altDataIdx}),  ...
                            subsref(obj, substruct('.', fds{fdidx}, '()', {rowidx(ri), ':', altDataIdx(ntreorder(ri,:))})) );
                    end
                end
                for infoii = 1:length(infoIdx)
                    obj = subsasgn(obj, substruct('.', obj.Properties.VarNames{infoIdx(infoii)}, '{}', {rowidx(ri)}), ...
                        subsref(obj, substruct('.', obj.Properties.VarNames{infoIdx(infoii)}, '{}', {rowidx(ri)}, '()', { ntreorder( ri, 1:nalt(rowidx(ri)) ) })) );
                end
            end
            
            %if snpeff is available
            if ~isempty(vcfData.snpEff)
                subSnpEff = vcfData.snpEff( toVcfIdx );
                neff = cellfun(@(x) size(x,1), subSnpEff);
                flatVarIdx = arrayfun(@(x) ones(neff(x),1)*x, (1:length(subSnpEff))', 'unif', 0);
                obj.varSnpEff = mat2dataset([vertcat(flatVarIdx{:}), strkey( vertcat(subSnpEff{:}), 'add' )], ...
                    'VarNames', ['varidx'; vcfData.snpEffFormat(:)]);
                
                obj.strkeyDB = strkey('dbname');
                %remove entries without annotations
                rmi = all(double(obj.varSnpEff(:, 2:end))==0,2);
                obj.varSnpEff(rmi,:) = [];
            end
            
            obj.VarType = cell(size(obj.Properties.VarNames));
            obj.VarType(:) = {'stat'};
            obj.VarType(ismember(obj.Properties.VarNames, [vcfData.attrName_cell; vcfData.attrName_mtx])) = {'info'};
            obj.VarType(ismember(obj.Properties.VarNames, fds)) = {'format'};        
            
            strfds = intersect(obj.vcfDesc(strcmp(obj.vcfDesc(:,3), 'String'),1), ...
                obj.Properties.VarNames);
            for i = 1:length(strfds)
                obj = subsasgn(obj, substruct('.', strfds{i}, '()', {cellfun(@isempty, obj.selfref(strfds{i}))}), {''});
            end
            
            
        end
       
        function obj = loadPresetIdx(obj)
            nvar = size(obj.selfref('loc'),1);
            %currently selected variants
            obj.subsetidx = true(nvar,1);
            
            % predefined sets of variants
            obj.presetidx.lowDNACov = all(obj.selfref('totalReadDNA') < obj.minReadDNA, 2);
            obj.presetidx.highDNACov = any(obj.selfref('totalReadDNA') > obj.maxReadDNA, 2);
            obj.presetidx.expressedEither = any(obj.selfref('totalReadRNA') >= obj.minReadRNA, 2); %expressed in either tumor or normal
            obj.presetidx.expressedBoth = all(obj.selfref('totalReadRNA') >= obj.minReadRNA, 2); %expressed in both tumor and normal
            
            obj.presetidx.multAlt = ~cellfun(@isempty, obj.selfref('alt','()', {':',2}));
            obj.presetidx.gatkindel = strcmp(obj.selfref('sets'), 'gatkIndel');
            obj.presetidx.chrMT = obj.selfref('loc', '()', {':',1}) == 25;
            obj.presetidx.DE01 = obj.selfref('DEpvalAdj') < 0.01;
            
            GT = strcat(obj.selfref('GT', '()', {':',1}), '-', obj.selfref('GT','()',{':',2}));
            GTset = strcat(GT, '-', obj.selfref('sets'));
            
            score = obj.selfref('Varscan_SSC');
            score(cellfun(@isempty, score)) = {'NaN'};
            score = str2double_fast(score);
            score(:,2) = obj.selfref('Sniper_SSC', '()', {':',2});
            
            obj.presetidx.distToCen500 = obj.selfref('distToCen') <= 500000;
            obj.presetidx.chr1_19_cen = (obj.selfref('loc','()',{':',1}) == 1 | ...
                obj.selfref('loc', '()', {':',1})==19) & ...
                obj.selfref('distToCen') < 10^6;
            obj.presetidx.homAltNormal = strcmp(GT, '1/1-0/1') | strcmp(GT, '1/1-0/0');
            obj.presetidx.varscanHet2Het = strcmp(GTset, '0/1-0/1-varscan');
            obj.presetidx.sniper30 = score(:,2) <= 30;
            obj.presetidx.chrY = obj.selfref('loc','()', {':',1}) == 24;            
        end
        
        function obj = selectVariant(obj, selection)
            % choose subset of variants
            % selection string example: 'expressedBoth & ~(multAlt | gatkindel |
            % chrMT)', each represents a logical index in obj.presetidx
            % selection can otherwise be a logical index to be directly
            % assigned to subsetidx
            
            if ischar(selection)
                obj.subsetdesc = selection;
                obj.subsetidx = eval(regexprep(selection, '([\w_]+)', 'obj.presetidx.$1'));
            elseif islogical(selection)
                if length(selection) == size(obj, 1)
                    obj.subsetdec = '';
                    obj.subsetidx = selection(:);
                else
                    fprintf('cannot set subsetidx: inconsistent dimension\n');
                end
            else
                fprintf('unknown selection input\n');
            end
            fprintf('select %d variants\n', sum(obj.subsetidx));
        end
                
        function availSet = availablePreset(obj)
            availSet = fieldnames(obj.presetidx);
        end
        
        function [idx, desc] = currentSet(obj)
            if all(obj.subsetidx)
                obj.subsetdesc = 'all';
            end
            idx = obj.subsetidx;
            desc = obj.subsetdesc;
        end
        
        function key = lockey(obj)
            key = strcat( arrayfun(@num2str, obj.selfref('loc', '()', {':',1}), 'unif', 0), '-', ...
                arrayfun(@num2str, obj.selfref('loc', '()', {':',2}), 'unif', 0) );
        end
        
        function writeText(obj, filename)            
            data = subsref(obj, substruct('()', {1,':'}));
            numeric = datasetfun(@isnumeric, data)';
            numeric(strcmp(obj.Properties.VarNames, 'loc')) = false;
            dsize = datasetfun(@(x) [size(x,2), size(x,3)], data, 'unif',0);
            dsize = cell2mat(dsize');
            [ur, ~, uidx] = unique(dsize, 'rows');
%             skip = find(all(bsxfun(@eq, ur, [2,2]),2) | all(bsxfun(@eq, ur, [2,4]),2));
            
            varnames = obj.Properties.VarNames(:)'; %row
            loc = double(obj.selfref('loc'));
            %numeric data, non-allelic
            singlenum = uidx == find( all(ur==1, 2)) & numeric;
            doublenum = uidx == find( all(bsxfun(@eq, ur, [2, 1]), 2)) & numeric;
            numdata = [loc, double(subsref(obj, substruct('()', {':', singlenum}))), ...
                double(datasetfun(@(x) x(:,1), subsref(obj, substruct('()', {':', doublenum})),'datasetoutput',true)), ...
                double(datasetfun(@(x) x(:,2), subsref(obj, substruct('()', {':', doublenum})),'datasetoutput',true)) ];
            numheader = [{'chrm', 'pos'}, varnames(singlenum), ...
                strcat(varnames(doublenum), '-', obj.sample{1}), ...
                strcat(varnames(doublenum), '-', obj.sample{2})];
            fprintf('print out numeric data: %s ...\n',strrep(filename, '.txt', '.num.txt'));
            dlmwriteplus(strrep(filename, '.txt', '.num.txt'), numdata, {}, numheader, sprintf('\t'), 'precision', '%.3f');
            
            %text data, non-allelic
            strings = datasetfun(@iscellstr, data)';
            singlestr = uidx == find( all(ur==1,2)) & strings;
            singlestr( ismember(varnames, { 'RPA','RU','VT','Varscan_SS','Varscan_SSC'} )) = false;
            doublestr = uidx == find( all(bsxfun(@eq, ur, [2,1]), 2)) & strings;
            doublestr( ismember(varnames, { 'indelRNA','indelDNA','AMQ','BQ','FREQ','FT','Varscan_DP4','IGT'} )) = false;
            strdata = [ arrayfun(@num2str, loc, 'unif', 0), ...
                cellstr( subsref(obj, substruct('()', {':', singlestr})) ), ...
                cellstr( datasetfun(@(x) x(:,1), subsref(obj, substruct('()', {':', doublestr})),'datasetoutput',true) ), ...
                cellstr( datasetfun(@(x) x(:,2), subsref(obj, substruct('()', {':', doublestr})),'datasetoutput',true) ) ];
            strheader = [{'chrm', 'pos'}, varnames(singlestr), ...
                strcat(varnames(doublestr), '-', obj.sample{1}), ...
                strcat(varnames(doublestr), '-', obj.sample{2})];
            fprintf('print out string data: %s ... \n', strrep(filename, '.txt', '.str.txt'));
            tabwrite(strrep(filename, '.txt', '.str.txt'), [strheader; strdata]);
            
            %allelic
            altidx = uidx == find(all(bsxfun(@eq, ur, [3,1]), 2));
            altdataidx = uidx == find(all(bsxfun(@eq, ur, [2,3]),2));
            data = [loc, ...
                double(datasetfun(@(x) x(:,1,1),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true)), ...
                double(datasetfun(@(x) x(:,1,2),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true)), ...
                double(datasetfun(@(x) x(:,1,3),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true)), ...
                double(datasetfun(@(x) x(:,2,1),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true)), ...
                double(datasetfun(@(x) x(:,2,2),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true)), ...
                double(datasetfun(@(x) x(:,2,3),subsref(obj, substruct('()', {':',altdataidx&numeric})), 'datasetoutput', true))];
            header = [{'alt1', 'alt2', 'alt3', 'chrm', 'loc'}, ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{1}, '-alt1'), ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{1}, '-alt2'), ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{1}, '-alt3'), ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{2}, '-alt1'), ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{2}, '-alt2'), ...
                strcat(varnames(altdataidx&numeric), '-', obj.sample{2}, '-alt3') ];
            
            fprintf('print out allelic data: %s ...\n', strrep(filename, '.txt', '.allelic.txt'));
            dlmwriteplus(strrep(filename, '.txt', '.allelic.txt'), data, ...
                cellstr(subsref(obj, substruct('()', {':',altidx}))), header, sprintf('\t'), 'precision', '%.3f');
            %tabwrite(strrep(filename, '.txt', '.allelic.txt'), [header; data]);
        end
        
        function [table, header] = getCurSelVarAnnotation(obj, fields, varidx)
            if nargin < 3, varidx = obj.subsetidx; end
            if nargin < 2, fields = []; end            
            
            varnames = obj.Properties.VarNames;
            varnames = varnames(:);
            if isempty(fields)                           
                varheader = {'loc'; 'totalReadDNA'; 'totalReadRNA'; ...
                    'altDnaAF'; 'altRnaAF'; 'DEpvalAdj'};
                varheader = [varheader; varnames(strcmpi(obj.VarType, 'info'))];
            else
                varheader = fields;
                varheader(~ismember(varheader, varnames)) = [];
            end
            
            addcol = 0;
            for i = 1:length(varheader)
                addcol = addcol + size(obj.selfref(varheader{i}),2);
            end
            
            if islogical(varidx)
                table = cell(sum(varidx), addcol);
                header = cell(1, size(table,2));
            else
                table = cell(length(varidx), addcol);
                header = cell(1, size(table,2));
            end
            
            coli = 1;
            for i = 1:length(varheader)
                addcol = size(obj.selfref(varheader{i}), 2);
                if isnumeric(obj.selfref(varheader{i}))
                    table(:, coli:coli+addcol-1) = arrayfun(@num2str, obj.selfref(varheader{i}, '()', {varidx,':',1}), 'unif', 0);
                else
                    table(:, coli:coli+addcol-1) = obj.selfref(varheader{i}, '()', {varidx,':',1});
                end
                header(coli:coli+addcol-1) = varheader(i);
                coli = coli + addcol;
            end
            
            table(cellfun(@isempty, table)) = {''};
            idx = cellfun(@isnumeric, table); 
            table(idx) = cellfun(@mat2str, table(idx), 'unif', 0);
            table = strrep(table, '''', '');
        end
        
        function [table, header] = getCurSelVarSnpEff(obj, fields, varidx, includeloc)
            if nargin < 4, includeloc = true; end
            if nargin < 3, varidx = obj.subsetidx; end
            if nargin < 2, fields = []; end 
            if isempty(fields)
                fields = setdiff(obj.varSnpEff.Properties.VarNames, {'varidx'});
            end
                       
            if islogical(varidx)
                varidx = find(varidx);
            end
            
            header = fields;
            rowidx = ismember(double(obj.varSnpEff.varidx), varidx);
            table = strkey( double(obj.varSnpEff(rowidx, fields)) );            
            
            if includeloc
                header = [{'chrm', 'loc'}, header];
                table = [ arrayfun(@num2str, ...
                    obj.selfref('loc', '()', {obj.varSnpEff.varidx(rowidx), ':'}), 'unif', 0), table];
            end
        end
        
        function [snpEffStat, table] = getCurSelVarSnpEffStat(obj, varidx, toFile)
            if nargin < 2, varidx = obj.subsetidx; end
            if nargin < 3, toFile = []; end
            
            if isempty(varidx), varidx = obj.subsetidx; end
            if islogical(varidx)
                rowidx = varidx(obj.varSnpEff.varidx);
            else
                rowidx = ismember(obj.varSnpEff.varidx, varidx);
            end
            fds = {'Effect','Effect_Impact','Functional_Class','Codon_Change',...
                'Amino_Acid_change','Transcript_BioType','Gene_Coding'};
            snpEffStat = cell(length(fds), 2);
            snpEffStat(:,1) = fds;
            for i = 1:length(fds)
                [count, tokens] = eleCounts(obj.varSnpEff.(fds{i})(rowidx));
                count(tokens==0) = [];
                tokens(tokens==0) = [];
                snpEffStat{i,2} = [strkey(tokens), num2cell(count)];
            end
            snpEffStat(cellfun(@isempty, snpEffStat(:,2)), :) = [];
            if ~isempty(toFile)
                ntokens = cellfun(@(x) size(x,1), snpEffStat(:,2));
                table = repmat({''}, max(ntokens)+1, 3*size(snpEffStat,1) - 1);
                table(1,1:3:end) = snpEffStat(:,1);
                for i = 1:size(snpEffStat,1)
                    table(2:ntokens(i)+1, (i-1)*3+1) = snpEffStat{i,2}(:,1);
                    table(2:ntokens(i)+1, (i-1)*3+2) = cellfun(@num2str, snpEffStat{i,2}(:,2), 'unif', 0);
                end
                tabwrite(toFile, table);
            else
                table = {};
            end
        end
        
        function data = selfref(obj, varname, varargin)
            data = subsref(obj, substruct('.', varname, varargin{:}));
        end
        
        function [varargout] = subsref(obj, S)
            if strcmp(S(1).type, '.')
                if isprop(obj, S(1).subs)
                    if length(S) == 1
                        varargout{1} = obj.(S(1).subs);
                    else
                        [varargout{1:nargout}] = subsref(obj.(S(1).subs), S(2:end));
                    end
                elseif ismethod(obj, S(1).subs)
                    if length(S) == 1
                        [varargout{1:nargout}] = eval(sprintf('%s(obj);', S(1).subs));
                    else
                        [varargout{1:nargout}] = eval(sprintf('%s(obj, S(2).subs{:});', S(1).subs));
                    end
                else
                    [varargout{1:nargout}] = subsref@dataset(obj, S);
                end
            else
                [varargout{1:nargout}] = subsref@dataset(obj, S);
            end
        end
        
        function obj = subsasgn(obj, S, B)
            if strcmp(S(1).type, '.') 
                if ismethod(obj, S(1).subs)
                    error('Cannot assign to method %s', S(1).subs);
                elseif ~isprop(obj, S(1).subs)
                    obj = subsasgn@dataset(obj, S, B);
                else
                    if length(S) == 1
                        obj.(S(1).subs) = B;
                    else
                        obj = subsasgn(obj.(S(1).subs), S(2:end));
                    end                    
                end
            else
                obj = subsasgn@dataset(obj, S, B);
            end                    
        end
        
    end
    methods (Static=true)
        function obj = loadobj(obj)
            obj = obj.loadPresetIdx();
        end
        
        function obj = saveobj(obj)
            obj.subsetidx = [];
            obj.subsetdesc = 'all';
            obj.presetidx = [];
            saveobj@super(obj);
        end
        
        function [vcfData, vcflockey] = loadVcf(patientSampleInfo, para)
            chrmNum = {'X', '23'; 'Y', '24'; 'MT', '25'; 'M', '25'};
            vcfData = loadStructData(sprintf('%svcf-%s.annt.mat',para.aCountDir, patientSampleInfo{1}));
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
            vcfData.attrName_cell = strrep(vcfData.attrName_cell, 'set', 'sets');
            vcfData.desc(:,1) = strrep(vcfData.desc(:,1), 'set', 'sets');
        end
        
        function [aDataArray, lockey, uloc] = loadAcountData(patientSampleInfo, para, vcflockey)
            if nargin < 3, vcflockey = []; end
            nsample = size(patientSampleInfo,1) - 1;
            aDataArray = cell(nsample, 2); %RNA, DNA
            lockey = cell(nsample,2);
            
            datatype = {'RNA', 'DNA'};
            for sampidx = 1:nsample
                for dataidx = 1:length(datatype)
                    matfn = sprintf('%sacount.%s.%s-%s.r%d.d%dk.q%d.mat', para.aCountDir, ...
                        datatype{dataidx}, patientSampleInfo{1}, ...
                        patientSampleInfo{sampidx+1}, para.r, para.(['d' datatype{dataidx}]), para.q);
                    aDataArray{sampidx,dataidx} = loadStructData(matfn);
                    aDataArray{sampidx,dataidx}.acount = sum(aDataArray{sampidx,dataidx}.acount, 3);
                    aDataArray{sampidx,dataidx}.ref = upper(aDataArray{sampidx,dataidx}.ref);
                    
                    tmp = num2cellstr(aDataArray{sampidx,dataidx}.loc);
                    lockey{sampidx,dataidx} = strcat(tmp(:,1), '-', tmp(:,2));
                end
            end
            
            if ~isempty(vcflockey)
                uloc = intersect(union(lockey{1,1}, lockey{2,1}), vcflockey);
            else
                uloc = union(lockey{1,1}, lockey{2,1});
            end
            
            fds = {'loc', 'numreads', 'ref', 'indels', 'acount'};
            for sampidx = 1:nsample %remove DNA acount entries, or RNA that are not in VCF
                rmi = ~ismember(lockey{sampidx,1}, uloc);
                for fdidx = 1:length(fds)
                    aDataArray{sampidx,1}.(fds{fdidx})(rmi,:,:) = [];
                end
                lockey{sampidx,1}(rmi) = [];
                
                rmi = ~ismember(lockey{sampidx,2}, uloc);
                for fdidx = 1:length(fds)
                    aDataArray{sampidx,2}.(fds{fdidx})(rmi,:,:) = [];
                end
                lockey{sampidx,2}(rmi) = [];
            end
        end
    end
    
    
end
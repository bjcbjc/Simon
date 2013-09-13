classdef VCF < dynamicprops
    
    
    properties
        filename;        
        sample;        
        desc;        
        maxnalt;
                                
        parseFormat = false;
        
        attrName_cell = {'CHROM'; 'ID'; 'REF'; 'ALT'; 'FILTER'};
        attrName_mtx = {'POS';'QUAL'};
        variantAttr_cell = {};
        variantAttr_mtx = [];
        formatData = {};
        allelicFormatField = {};
        pyparser = 'python /nethome/bjchen/Projects/Simon/pyscripts/parseFormat.py ';
        pyformatfiledir = '/formatdata/'
        snpEffFormat = '';
        snpEff = {};
                
    end
    
    methods
        function obj = VCF(fn, parseFormat, pyformatfiledir, goodVCF, blocksize, quick)
            if nargin >= 2
                obj.parseFormat = parseFormat;
            end
            if nargin == 3
                obj.pyformatfiledir = pyformatfiledir;
            end           
            if nargin < 4
                goodVCF = false;
            end
            if nargin < 5
                blocksize = 10000;
            end
            if nargin < 6
                quick = false;
            end
            
            obj.filename = fn;
            
            if quick
                obj = obj.readReshapedVCF();     
                obj.cleanup();
                return;
            elseif goodVCF
                if blocksize == 0
                    obj = obj.readGoodVCF0();
                else
                    obj = obj.readGoodVCF(blocksize);
                end
                obj.cleanup();
                return;
            end
            
            header = obj.obtainHeader();
            
            fprintf('read vcf file...\n');
            tic;
            txt = parseText(fn, 'skip', length(header), 'ncolname', 0, 'nrowname', 0);
            toc;
            txt = txt.text;
            
            obj.parseHeader(header);            
            
            fprintf('parse variant info...\n');            
            obj.parseVariantInfo(txt, header{end});            

            %get max number of ALT 
            colidx = strcmp(obj.attrName_cell, 'ALT');
            obj.maxnalt = max( cellfun(@(x) length(strfind(x,',')), obj.variantAttr_cell(:, colidx)) ) + 1;
                        
            if obj.parseFormat
                tokenidx = strfind(obj.filename, '/');
                if isempty(tokenidx)
                    basedir = '.';
                else
                    basedir = obj.filename(1:tokenidx(end));
                end
                if ~exist([basedir obj.pyformatfiledir], 'file')
                    system(sprintf('mkdir %s', [basedir obj.pyformatfiledir]));
                end
                fprintf('generate easy FORMAT data...\n');
                tic;
                system(sprintf('%s %s %s %d', obj.pyparser, obj.filename, [basedir obj.pyformatfiledir], obj.maxnalt));
                toc;
                
                fprintf('parse FORMAT data...\n');
                tic;
                obj.parseFormatData([basedir obj.pyformatfiledir], false);
                toc;
            end
            obj.cleanup();
        end
        
        function [header, l] = obtainHeader(obj)
            header = {};
            f = fopen(obj.filename, 'r');
            l = fgetl(f);
            while l ~= -1
                if strcmp(l(1), '#')
                    header = [header; l];
                else
                    break;
                end
                l = fgetl(f);
            end
            fclose(f);
        end
        
        function parseHeader(obj, header)
 
            if nargin < 2, header = obj.obtainHeader(); end
            if isempty(header), header = obj.obtainHeader(); end
            
            %{name, number, type, desc}
            infotable = obj.parseMeta(header, {'INFO'});
            obj.desc = infotable; %name, number, type, desc
            if ~isempty(infotable)
                obj.desc(:,end+1) = {'info'};
                numeric = (ismember(infotable(:,3), {'Float', 'Integer'}) ...
                    & strcmp(infotable(:,2), '1')) | strcmp(infotable(:,3), 'Flag');
                
                obj.attrName_mtx = [obj.attrName_mtx; infotable(numeric, 1)];
                obj.attrName_cell = [obj.attrName_cell; infotable(~numeric, 1)];
            end
            formattable = obj.parseMeta(header, {'FORMAT'});
            if ~isempty(formattable)
                formattable(:,end+1) = {'format'};
                obj.desc = [obj.desc; formattable];
            end
        end

        function parseVariantInfo(obj, text, header)
            numvar = size(text, 1);
                        
            h = regexp(header(2:end), '\t', 'split');
            obj.variantAttr_cell = cell(numvar, size(obj.attrName_cell,1));
            obj.variantAttr_mtx = NaN(numvar, size(obj.attrName_mtx,1));
            [~, idx] = ismember(obj.attrName_mtx, obj.desc(:,1));
            obj.variantAttr_mtx(:, strcmp(obj.desc(idx(idx~=0),3), 'Flag')) = 0;
            
            [~, colidx_cell] = ismember(obj.attrName_cell, h);
            obj.variantAttr_cell(:, colidx_cell~=0) = text(:, colidx_cell(colidx_cell~=0));
            
            [~, colidx_mtx] = ismember(obj.attrName_mtx, h);
            obj.variantAttr_mtx(:, colidx_mtx~=0) = str2double_fast(text(:, colidx_mtx(colidx_mtx~=0)));
                        
            idx = strcmp(h, 'INFO');
            %t = regexp(s, ';([\w\.]+)=*([\w\.\,\(\)\|\-]*)|^([\w\.]+)=*([\w\.\,\(\)\|\-]*)', 'tokens');
            %eff = regexp(a, '(\w+)\((\w+)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\|(\w*)\)','tokens');
            %effformat = regexp(v.desc{strcmp(v.desc(:,1), 'EFF'),4}, '(\w+)[\(\s]+(\w+)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*) \[','tokens');
            fprintf('cell info data; ');
            tic;
            for fdidx = find(colidx_cell==0)'
                if strcmpi(obj.attrName_cell{fdidx}, 'EFF')
                    val = regexp(text(:,idx), sprintf(';%s=(\\S+);|;%s=(\\S+)$|^%s=(\\S+)', ...
                        obj.attrName_cell{fdidx}, obj.attrName_cell{fdidx},obj.attrName_cell{fdidx}), 'tokens', 'once');
                else                    
                    val = regexp(text(:,idx), sprintf(';%s=([\\w.\\-\\,\\(\\)\\|\\/\\*]+)|^%s=([\\w.\\-\\,\\(\\)\\|\\/\\*]+)', ...
                        obj.attrName_cell{fdidx}, obj.attrName_cell{fdidx}), 'tokens', 'once');
                end
               vi = ~cellfun(@isempty, val);               
               obj.variantAttr_cell(vi, fdidx) = ...
                   cellfun(@(x) x{1}, val(vi), 'unif', 0); %flat nested cell arrays
               nameidx = strcmp(obj.desc(:,1), obj.attrName_cell{fdidx});
               if strcmp(obj.desc(nameidx, 2), 'A') && ...
                       ismember(obj.desc(nameidx, 3), {'Integer','Float'})
                   vi = find(vi);
                   singleval = cellfun(@isempty, strfind(obj.variantAttr_cell(vi, fdidx), ','));
                   obj.variantAttr_cell(vi(singleval), fdidx) = ...
                       num2cell(str2double_fast(obj.variantAttr_cell(vi(singleval), fdidx)));
                   obj.variantAttr_cell(vi(~singleval), fdidx) = ...
                       cellfun(@str2num, obj.variantAttr_cell(vi(~singleval), fdidx), 'unif', 0);
               end
            end
            toc;
            
            fprintf('mtx info data; ');
            tic;
            for fdidx = find(colidx_mtx==0)'                
                if strcmp( obj.desc(strcmp(obj.desc(:,1), obj.attrName_mtx{fdidx}), 3), 'Flag')
                    vi = ~cellfun(@isempty, strfind(text(:,idx), obj.attrName_mtx{fdidx}) );
                    obj.variantAttr_mtx(vi, fdidx) = 1;
                else                
                    val = regexp(text(:,idx), sprintf(';%s=([\\w\\.\\-\\,]+);|;%s=([\\w\\.\\-\\,]+)$|^%s=([\\w\\.\\-\\,]+)', ...
                        obj.attrName_mtx{fdidx},obj.attrName_mtx{fdidx},obj.attrName_mtx{fdidx}), 'tokens', 'once');
                    vi = ~cellfun(@isempty, val);
                    val(vi) = cellfun(@(x) x{1}, val(vi), 'unif', 0); %flat nested cell arrays
                    obj.variantAttr_mtx(vi, fdidx) = str2double_fast( val(vi));
                end
            end    
            toc;
            
            %sample data            
            idx = find(strcmp(h, 'FORMAT'));            
            obj.sample = h(idx+1:end)';                                    
            
            %if 'EFF' is in the data, parse it
            [~, effcolidx] = ismember('EFF', obj.desc(:,1));
            if effcolidx ~= 0
                fprintf('EFF data; ');
                tic;
                % get format of eff
                obj.snpEffFormat = regexp(obj.desc{effcolidx,4}, '(\w+)[\(\s]+(\w+)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*)[\|\s]+(\w*) \[','tokens','once');
                ncol = length(obj.snpEffFormat);                
                attridx = strcmpi(obj.attrName_cell, 'EFF');
                obj.snpEff = regexp(obj.variantAttr_cell(:,attridx), ',','split');
                obj.snpEff = cellfun(@(x) regexp(x, '[\(\)\|]', 'split'), obj.snpEff, 'unif', 0);
                for i = 1:length(obj.snpEff)
                    if isempty(obj.snpEff{i}{1}{1})
                        obj.snpEff{i} = repmat({''}, 1, ncol);
                    else
                        obj.snpEff{i} = cellfun(@(x) x(1:ncol), obj.snpEff{i}, 'unif', 0);
                        obj.snpEff{i} = vertcat(obj.snpEff{i}{:});
                    end
                end                
                
%                 obj.snpEff = regexp(obj.variantAttr_cell(:,attridx), ...
%                     '(\w+)\((\w+)\|(\w*)\|([\w\/]*)\|([\w\*]*)\|(\w*)\|([\w\.\-]*)\|(\w*)\|(\w*)\|([\w\.\-]*)\|(\w*)\|(\w*)\)','tokens');
%                 obj.snpEff = cellfun(@(x) vertcat(x{:}), obj.snpEff, 'unif', 0);
                
                obj.attrName_cell(attridx) = [];
                obj.variantAttr_cell(:, attridx) = [];                
                toc;
            end

            strinfotype = any(cellfun(@ischar, obj.variantAttr_cell),1);
            obj.variantAttr_cell(bsxfun(@and, cellfun(@isempty, obj.variantAttr_cell), ...
                strinfotype) ) = {''};
                
        end
        
        function table = parseMeta(obj, header, metatype)
            if nargin < 2, header = obj.obtainHeader(); end            
            if nargin < 3, metatype = {'INFO', 'FORMAT'}; end
            if isempty(header), header = obj.obtainHeader(); end 
            
            template = '=<ID=([\w\.]+),Number=([\w\.]+),Type=(\w+),Description="(.*)">';
            table = {};
            for i = 1:length(metatype)                
                subtable = regexp(header, [metatype{i} template],'tokens', 'once');
                subtable = subtable(~cellfun(@isempty, subtable));
                subtable = vertcat(subtable{:});                
                if isempty(subtable), continue; end
                subtable(:,1) = strrep(subtable(:,1), '.', '_');                
                table = [table; subtable];
            end
        end
        
        function parseFormatData(obj, formatdatadir, makesparse)
            if nargin < 3, makesparse = false; end
            headers = obj.obtainHeader();
            
            foldertokenidx = strfind(obj.filename, '/');
            if isempty(foldertokenidx)
                basefilename = obj.filename;
            else
                basefilename = obj.filename(foldertokenidx(end)+1:end);
            end
            formatFieldStruct = obj.getFormatStructure(headers);
            
            nvar = size(obj.variantAttr_mtx, 1);
            nsamp = length(obj.sample);
            %read numeric as matrix, sparse 
            % !!! but the sparse means missing values and real values can
            % be 0; switch NaN <-> 0 in storage
            for fdidx = 1:length(formatFieldStruct.numField)
                fn = sprintf('%s/%s.%s', formatdatadir, basefilename, formatFieldStruct.numField{fdidx});
                if ~exist(fn, 'file'), continue; end
                if formatFieldStruct.numFieldDim(fdidx) == 1
                    fmtdata = parseText(fn, 'nrowname', 0, 'ncolname', 0, ...
                        'ncol', 3, 'numeric', true, 'treatasempty', {'.'});
                    if makesparse
                        val = fmtdata.text(:,3);
                        nans = isnan(val);
                        val(val==0) = NaN;
                        val(nans) = 0;
                        obj.formatData.(formatFieldStruct.numField{fdidx}) = ...
                            sparse( fmtdata.text(:,1), fmtdata.text(:,2), ...
                            val, nvar, nsamp );
                    else
                        obj.formatData.(formatFieldStruct.numField{fdidx}) = NaN(nvar, nsamp);
                        obj.formatData.(formatFieldStruct.numField{fdidx})( ...
                            sub2ind([nvar, nsamp], fmtdata.text(:,1), fmtdata.text(:,2)) ) = fmtdata.text(:,3);
                    end
                elseif isnan(formatFieldStruct.numFieldDim(fdidx))
                    fmtdata = parseText(fn, 'nrowname', 0, 'ncolname', 0, ...
                        'ncol', 3, 'numericcol', [1, 2]);
                    obj.formatData.(formatFieldStruct.numField{fdidx}) = cell(nvar, nsamp);
                    obj.formatData.(formatFieldStruct.numField{fdidx})( ...
                        sub2ind([nvar, nsamp], fmtdata.numtext(:,1), fmtdata.numtext(:,2)) ) = ...
                        fmtdata.text(:,1);
                else %multicolumns
                    fmtdata = parseText(fn, 'nrowname', 0, 'ncolname', 0, ...
                        'numeric', true, 'treatasempty', {'.'});
                    if makesparse
                        val = fmtdata.text(:,3:end);
                        nans = isnan(val);
                        val(val==0) = NaN;
                        val(nans) = 0;
                        ndim = size(val,2);
                        obj.formatData.(formatFieldStruct.numField{fdidx}) = cell(ndim, 1);
                        for dimidx = 1:ndim
                            obj.formatData.(formatFieldStruct.numField{fdidx}){dimidx} = ...
                            sparse( fmtdata.text(:,1), fmtdata.text(:,2), ...
                            val(:, dimidx), nvar, nsamp );
                        end
                    else
                        [ndata, ndim] = size(fmtdata.text);
                        ndim = ndim - 2;
                        obj.formatData.(formatFieldStruct.numField{fdidx}) = NaN(nvar, nsamp, ndim);
                        for dimidx = 1:ndim
                            obj.formatData.(formatFieldStruct.numField{fdidx})( ...
                                sub2ind([nvar, nsamp, ndim], fmtdata.text(:,1), ...
                                fmtdata.text(:,2), repmat(dimidx,ndata,1)) ) = fmtdata.text(:, dimidx+2);
                        end
                    end
                end
            end
            
            %read string as cell arrays
            for fdidx = 1:length(formatFieldStruct.strField)
                fn = sprintf('%s/%s.%s', formatdatadir, basefilename, formatFieldStruct.strField{fdidx});
                if ~exist(fn, 'file'), continue; end
                fmtdata = parseText(fn, 'nrowname', 0, 'ncolname', 0, ...
                    'ncol', 3, 'numericcol', [1, 2]);
                obj.formatData.(formatFieldStruct.strField{fdidx}) = cell(nvar, nsamp);
                obj.formatData.(formatFieldStruct.strField{fdidx})( ...
                    sub2ind([nvar, nsamp], fmtdata.numtext(:,1), fmtdata.numtext(:,2)) ) = ...
                    fmtdata.text(:,1);                
                obj.formatData.(formatFieldStruct.strField{fdidx})( ...
                    cellfun(@isempty, obj.formatData.(formatFieldStruct.strField{fdidx})) ) = {''};
            end
        end
        
        function formatFieldStruct = getFormatStructure(obj, header)
            if nargin < 2, header = obj.obtainHeader(); end
                
            formatFieldTable = obj.parseMeta(header, {'FORMAT'});            
            
            rowidx = ismember(formatFieldTable(:,3), {'String','Character'});
            formatFieldStruct.strField = formatFieldTable(rowidx, 1);
            
            rowidx = ismember(formatFieldTable(:,3), {'Float','Integer'});
            if any(rowidx)
                formatFieldStruct.numField = formatFieldTable(rowidx, 1);
                formatFieldStruct.numFieldDim = str2double_fast(formatFieldTable(rowidx, 2));
                formatFieldStruct.numFieldDim( strcmp(formatFieldTable(rowidx,2), 'A') ) = obj.maxnalt;
                formatFieldStruct.numFieldDim( ~cellfun(@isempty, strfind(lower(formatFieldTable(rowidx,4)), 'allelic depth') ) ) = obj.maxnalt + 1;
            end
            
            rowidx = ismember(formatFieldTable(:,3), {'Flag'});
            if any(rowidx)
                formatFieldStruct.boolField = formatFieldTable(rowidx, 1);
                formatFieldStruct.boolFieldDim = str2double_fast(formatFieldTable(rowidx, 2));
                formatFieldStruct.boolFieldDim( strcmp(formatFieldTable(rowidx,2), 'A') ) = obj.maxnalt;
                formatFieldStruct.boolFieldDim( ~cellfun(@isempty, strfind(lower(formatFieldTable(rowidx,4)), 'allelic depth') ) ) = obj.maxnalt + 1;
            end
            
            allelic = ~cellfun(@isempty, strfind(lower(formatFieldTable(:,4)), 'allelic depth') ) | ...
                strcmp(formatFieldTable(:,2), 'A');                        
            obj.allelicFormatField = formatFieldTable(allelic,1);
        end
                
        function consolidateFormat(obj, samplegroup, groupname)
            %consolidate format data by samplegroup
            %samplegroup: cell, #group x 1, each group contains cell array
            %of 'sample id' to refer to obj.sample
            %
            %if a variant has data in more than one samples in a sample
            %group, the data from the latter sample will overwrite the
            %previous ones
            %
            ngroup = length(samplegroup);
            fds = fieldnames(obj.formatData);
            
            if length(groupname) ~= ngroup
                fprintf('error: #groupname ~= #groups\n');
                return
            end
            
            removeidx = [];
            for groupidx = 1:ngroup
                [~, sampidx] = ismember(samplegroup{groupidx}, obj.sample);
                sampidx(sampidx==0) = [];
                if length(sampidx) < 2, continue; end    
                
                for fdidx = 1:length(fds)
                    if iscell(obj.formatData.(fds{fdidx}))
                        for si = 2:length(sampidx) %loop over samples in this group
                            idx = ~cellfun(@isempty, obj.formatData.(fds{fdidx})(:, sampidx(si)));
                            %might overwrite data if they exists for previous
                            %samples in the same group
                            obj.formatData.(fds{fdidx})(idx, sampidx(1)) = ...
                                obj.formatData.(fds{fdidx})(idx, sampidx(si));
                        end
                    else
                        for si = 2:length(sampidx) %loop over samples in this group
                            idx = any(squeeze(~isnan(obj.formatData.(fds{fdidx})(:, sampidx(si), :))),2);
                            obj.formatData.(fds{fdidx})(idx, sampidx(1), :) = ...
                                obj.formatData.(fds{fdidx})(idx, sampidx(si), :);
                        end
                    end                    
                end 
                removeidx = [removeidx; reshape(sampidx(2:end), length(sampidx(2:end)), 1)];                    
            end
            %remove data
            if ~isempty(removeidx)
                for fdidx = 1:length(fds)
                    obj.formatData.(fds{fdidx})(:, removeidx, :) = [];
                    obj.formatData.(fds{fdidx}) = squeeze(obj.formatData.(fds{fdidx}));
                end
                obj.sample = groupname;
            end
        end
        
        function cleanup(obj)            
            rmfield = all(cellfun(@isempty, obj.variantAttr_cell), 1);
            obj.variantAttr_cell(:, rmfield) = [];
            obj.attrName_cell(rmfield) = [];            
            strcol = find(arrayfun(@(i) iscellstr(obj.variantAttr_cell(:,i)), 1:size(obj.variantAttr_cell,2)));
            if ~isempty(strcol)
                rmfield = all( strcmp('.', obj.variantAttr_cell(:, strcol)), 1);
                obj.variantAttr_cell(:, strcol(rmfield)) = [];
                obj.attrName_cell(strcol(rmfield), :) = [];
            end
            rmfield =  all(isnan(obj.variantAttr_mtx),1);
            obj.variantAttr_mtx(:, rmfield) = [];
            obj.attrName_mtx(rmfield) = [];
            
%             if ~isempty(obj.formatData)
%                 rmfield = squeeze(all(all(cellfun(@isempty, obj.formatData),1),2));
%                 obj.formatData(:, :, rmfield) = [];
%                 obj.formatfield(rmfield, :) = [];
%                 strcol = find(arrayfun(@(i) iscellstr(obj.formatData(:,:,i)), 1:size(obj.formatData,3)));
%                 if ~isempty(strcol)
%                     rmfield = squeeze(all(all( strcmp('.', obj.formatData(:,:,strcol)), 1), 2));
%                     obj.formatData(:,:,strcol(rmfield)) = [];
%                     obj.formatfield(strcol(rmfield)) = [];
%                 end
%             end
        end
        
        function obj = readGoodVCF(obj, blocksize)
            if nargin < 2, blocksize = 10000; end
            %FORMAT string is the same for all the records
            [header, firstvarline] = obj.obtainHeader();
            obj.parseHeader(header);
            fields = textscan(header{end}(2:end), '%s', 'delimiter', sprintf('\t'));
            fields = fields{1};
            nsample = length(fields) - 9;
            obj.sample = fields(10:end);
            [~, nvariant] = system(sprintf('grep -c -v ^# %s', obj.filename));
            nvariant = str2double(nvariant);
            
            fields = textscan(firstvarline, '%s', 'delimiter', sprintf('\t'));
            fields = fields{1};
            varformatstr = '%s %f %s %s %s %f %s %s'; %chrom, pos, id, ref, alt, qual, filter, info            
            formatfields = textscan(fields{9}, '%s', 'delimiter', ':');            
            formatfields = formatfields{1};
            nformatfields = length(formatfields);
            [~, formatidx] = ismember(formatfields, obj.desc(:,1));            
            for i = 1:nformatfields
                varformatstr = strcat(varformatstr, '%s');
            end
            sampleformatstr = '';
            for i = 1:nformatfields
                if ismember(lower(obj.desc(formatidx(i), 3)), {'float', 'integer'}) && ...
                    isnumericstring(obj.desc{formatidx(i),2})
                    ndataval = str2double(obj.desc{formatidx(i), 2});
                    if ndataval == 1
                        sampleformatstr = strcat(sampleformatstr, '%f');
                    else
                        sampleformatstr = strcat(sampleformatstr, '%s');
                    end
                    obj.formatData.(formatfields{i}) = NaN(nvariant, nsample, ndataval);
                else                    
                    sampleformatstr = strcat(sampleformatstr, '%s');
                    if ismember(formatfields{i}, {'GT', 'IGT'})
                        obj.formatData.(formatfields{i}) = NaN(nvariant, nsample);
                    else
                        obj.formatData.(formatfields{i}) = repmat({''}, nvariant, nsample);
                    end
                end
            end
            varformatstr = strcat(varformatstr, repmat(sampleformatstr, 1, nsample));
                           
            nmtx = length(obj.attrName_mtx);
            ncell = length(obj.attrName_cell);
            obj.variantAttr_mtx = NaN(nvariant, nmtx);
            obj.variantAttr_cell = repmat({''}, nvariant, ncell);            
            
            colheader = textscan(header{end}(2:end), '%s', 'delimiter', sprintf('\t'));
            colheader = colheader{1};
            colheader = colheader(:)';
            [~, colidx_mtx] = ismember(colheader, obj.attrName_mtx);
            [~, colidx_cell] = ismember(colheader, obj.attrName_cell);
            
            reptb = [strcat(obj.attrName_mtx, '='), strcat(arrayfun(@num2str, (1:nmtx)', 'unif',0),'=')];
            reptb = [reptb; [strcat(obj.attrName_cell, '='), strcat(arrayfun(@num2str, (1:ncell)'+nmtx,'unif',0), '=')]];
            
            lastreadvar = 0;
                        
            nheader = length(header);
            f = fopen(obj.filename);
%             tic;
            while ~feof(f)               
                d = textscan(f, varformatstr, blocksize, 'delimiter',sprintf('\t:'), 'headerlines', nheader, 'treatasempty', '.');
                nheader = 0; %the following iteration has no header lines
                nreadvar = length(d{1});
                
                fillidx = lastreadvar+1:lastreadvar+nreadvar;
                obj.variantAttr_mtx(fillidx, colidx_mtx(colidx_mtx~=0)) = cell2mat( ...
                    arrayfun(@(idx) d{idx}, find(colidx_mtx), 'unif', 0) );
                
                for i = 1:length(colidx_cell)
                    if colidx_cell(i) ~= 0
                        obj.variantAttr_cell(fillidx, colidx_cell(i)) = d{i};
                    end
                end
                
                %format data
                for i = 1:nformatfields
                    descidx = strcmp(formatfields{i}, obj.desc(:,1)) & strcmp(obj.desc(:,5), 'format');
                    tmp = d( 8+nformatfields+i:nformatfields:end );
                    tmp = horzcat(tmp{:});                    
                    if ismember(formatfields{i}, {'GT', 'IGT'})
                        obj.formatData.(formatfields{i})(fillidx, :) = strkey(tmp, 'add', 'lookup','/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat','savelater');
                    elseif ismember(obj.desc{descidx, 3}, {'Integer','Float'}) && ...
                            isnumericstring(obj.desc{descidx, 2})
                        ndataval = str2double(obj.desc{descidx,2});
                        if ndataval > 1
                            for sampi = 1:nsample
                                tmp(:, sampi) = strrep(tmp(:, sampi), '.', strjoin(repmat({'NaN'},1,ndataval), ','));
                                obj.formatData.(formatfields{i})(fillidx, sampi, :) = ...
                                    str2double_fast(tmp(:,sampi), strcat(strjoin(repmat({'%g'},1,ndataval),','), '#'));
                            end
                        else
                            obj.formatData.(formatfields{i})(fillidx, :) = tmp;
                        end
                    else
                        obj.formatData.(formatfields{i})(fillidx, :) = tmp;
                    end
                    clear tmp
                end
                
                %INFO
                d = d{8};
                
                if all(strcmp(d, '.')) %no info data
                    lastreadvar = lastreadvar+nreadvar;
                    continue;
                end
                
                %flag-info                
                flags = intersect(obj.desc(strcmp(obj.desc(:,3), 'Flag'), 1), obj.attrName_mtx);
                for i = 1:length(flags)
                    colidx = strcmp(obj.attrName_mtx, flags{i});
                    obj.variantAttr_mtx(fillidx, colidx) = double(...
                        ~cellfun(@isempty, strfind(d, flags{i})) );
                end
                                                
                for i = 1:size(reptb, 1)
                    d = strrep(d, reptb{i,1}, reptb{i,2});
                end
                                
                %infos = regexp(d, '([\w\.]+)=([\w\+\-\.\,]+)', 'tokens');
                infos = regexp(d, '(\d+)=([\w\+\-\.\,]+)', 'tokens');
                infos = cellfun(@(x) vertcat(x{:}), infos, 'unif', 0);
                nval = cellfun(@(x) size(x,1), infos);
                infos = vertcat(infos{:}); %flat the entire info data
                
                if ~isempty(infos)
                    varidx = cell2mat(arrayfun(@(i) repmat(i, nval(i-fillidx(1)+1), 1), fillidx', 'unif', 0));
                    infoidx = str2double_fast(infos(:,1));
                    infos = infos(:,2);

                    %num-info
                    curcolidx = infoidx <= nmtx;
                    obj.variantAttr_mtx( sub2ind([nvariant, nmtx], varidx(curcolidx), infoidx(curcolidx)) ) = ...
                        str2double_fast(infos( curcolidx ) );

                    %cell-info
                    curcolidx = ~curcolidx;
                    obj.variantAttr_cell( sub2ind([nvariant, ncell], varidx(curcolidx), infoidx(curcolidx)-nmtx) ) = infos( curcolidx ) ;                                                
                end
                lastreadvar = lastreadvar+nreadvar;
                clear d infos varidx infoidx
            end
%             toc;
            fclose(f);
            
            strkey('savenow');
            colidx = strcmp(obj.attrName_cell, 'ALT');
            obj.maxnalt = max( cellfun(@(x) length(strfind(x,',')), obj.variantAttr_cell(:, colidx)) ) + 1;
            obj.parseFormat = true;
        end
        
        function obj = readReshapedVCF(obj)
            header = obj.obtainHeader();
            obj.parseHeader(header);
            fields = textscan(header{end}(2:end), '%s', 'delimiter', sprintf('\t'));
            fields = fields{1};
            nsample = length(fields) - 9;
            obj.sample = fields(10:end);
            
            datatype = {'mtx', 'cell'};
            
            for dataidx = 1:2
                f = fopen([obj.filename '.' datatype{dataidx}], 'r');
                l = textscan(fgetl(f), '%s', 'delimiter', '\t');
                fclose(f);
                l = l{1};

                vartoken = regexp(l, '(\w+)_(\w+)', 'tokens');
                fmtidx = ~cellfun(@isempty, vartoken);
                vartoken = cellfun(@(x) x{1}, vartoken(fmtidx), 'unif', 0);
                vartoken = vertcat(vartoken{:}); %sample, varname
                if strcmp(datatype{dataidx}, 'mtx')
                    tb = readtable([obj.filename '.' datatype{dataidx}], 'filetype', 'text', 'delimiter', '\t', ...
                        'format', repmat('%f', 1, length(l)), 'treatasempty', '.');
                else
                    tb = readtable([obj.filename '.' datatype{dataidx}], 'filetype', 'text', 'delimiter', '\t', ...
                        'format', repmat('%s', 1, length(l)));
                end
                obj.(['attrName_' datatype{dataidx}]) = l(~fmtidx);
                obj.(['variantAttr_' datatype{dataidx}]) = tb{:, ~fmtidx};
                [fmtfds, ~, fmtuidx] = unique(vartoken(:,2));
                nvariant = size(obj.(['variantAttr_' datatype{dataidx}]),1);
                for i = 1:length(fmtfds)
                    colidx = find(fmtuidx == i) + sum(~fmtidx);
                    ndim = length(colidx)/nsample;
                    obj.formatData.(fmtfds{i}) = permute(reshape(tb{:, colidx}, nvariant, ndim, nsample), [1 3 2]);
                    colidx = find(fmtuidx == i);
                    if any(strcmp(obj.sample, vartoken(colidx([1 ndim+1]), 1)))
                        [~, sampi] = ismember(obj.sample, vartoken(colidx([1 ndim+1]), 1));
                        obj.formatData.(fmtfds{i}) = obj.formatData.(fmtfds{i})(:, sampi,:);
                    end                    
                    
%                     if ismember(fmtfds{i}, {'GT', 'IGT'})
%                         obj.formatData.(fmtfds{i}) = strkey( ...
%                             obj.formatData.(fmtfds{i}), 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat', 'savelater');
%                     end
                end
                clear tb                
            end            
%             strkey('savenow');
            colidx = strcmp(obj.attrName_cell, 'ALT');
            obj.maxnalt = max( cellfun(@(x) length(strfind(x,',')), obj.variantAttr_cell(:, colidx)) ) + 1;
            obj.parseFormat = true;
        end
        
        function keys = lockey(obj)
            keys = textscan(int2str(obj.variantAttr_mtx(:,1)'), '%s');
            keys = keys{1};
            keys = strcat(obj.variantAttr_cell(:,1), '-', keys);
%             keys = strcat(obj.variantAttr_cell(:,1), '-', ...
%                 arrayfun(@num2str, obj.variantAttr_mtx(:,1), 'unif', 0));
        end
    end
    
end


%         function hasvalue = parseText(obj, text, header)
%             numvar = size(text, 1);
%             %preallocate memory for headers
%             allocstr = '''CHROM'',[], ''POS'', [], ''ID'', [], ''REF'', [], ''ALT'', [], ''QUAL'', [], ''FILTER'', []';
%                         
%             fds = fieldnames(obj.desc);
%             for i = 1:length(fds)                
%                 allocstr = strcat(allocstr, sprintf(',''%s'',[]',fds{i}));
%             end          
%             %have to do two steps otherwise matlab complains!
%             tmp(numvar,1) = eval(sprintf('struct(%s);',allocstr));
%             obj.variants = tmp;
%             
%             fds = fieldnames(obj.variants);
%             hasvalue = false(length(fds), 1);
%             
%             h = textscan(strrep(header,'#',''), '%s');
%             h = h{1}; 
%             obj.sample = h{end};
%             
%             for i = 1:numvar
%                 %initialize all flag field = false
%                 for fi = 1:length(obj.flagfield)
%                     obj.variants(i).(obj.flagfield{fi}) = false;
%                 end
%                 %parse
%                 %s = textscan(text{i}, '%s');
%                 %s = s{1};
%                 s = text(i,:);
%                 for j = 1:length(s)
%                     if ismember(h{j}, fds)
%                         hasvalue(ismember(fds,h{j})) = true;
%                         if ismember(h{j}, obj.strfield)
%                             obj.variants(i).(h{j}) = s{j};
%                         else
%                             obj.variants(i).(h{j}) = str2double_fast(s{j});
%                         end
%                     elseif strcmpi(h{j}, 'INFO')
%                         t = textscan(s{j}, '%s', 'delimiter', ';');
%                         t = t{1};
%                         for ti = 1:length(t)
%                             data = textscan(t{ti}, '%s', 'delimiter', '=');
%                             data = data{1};
%                             if ~ismember(data{1}, fds)
%                                 continue
%                             end
%                             hasvalue(ismember(fds, data{1})) = true;
%                             if ismember(data{1}, obj.flagfield)
%                                 obj.variants(i).(data{1}) = true;
%                             elseif ismember(data{1}, obj.strfield)
%                                 obj.variants(i).(data{1}) = data{2};
%                             else
%                                 obj.variants(i).(data{1}) = eval(sprintf('[%s]'';',data{2}));
%                             end
%                         end
%                     elseif strcmpi(h{j}, 'FORMAT')
%                         infofd = textscan(s{j}, '%s', 'delimiter', ':');
%                         infofd = infofd{1};  
%                         for ii = 1:length(infofd)
%                             infofd{ii} = strcat(infofd{ii}, '_s');
%                         end
%                     elseif strcmpi(h{j}, obj.sample)
%                         infodata = textscan(s{j}, '%s', 'delimiter', ':');
%                         infodata = infodata{1};
%                         if length(infofd) ~= length(infodata)
%                             continue;
% %                             error('info length mismatch, record %d', i);
%                         end
%                         ifnumfd = ~ismember(infofd, obj.strfield);
%                         infodata(strcmp(infodata, '.')) = {''};
%                         infodata(ifnumfd) = strrep(infodata(ifnumfd), '%', '/100');
%                         infodata(ifnumfd) = strrep(infodata(ifnumfd), ',', ' ');
%                         infodata(ifnumfd) = strrep(infodata(ifnumfd), '.', '');
%                         for ii = 1:length(infofd)
%                             hasvalue(ismember(fds, infofd{ii})) = true;
%                             if ismember(infofd{ii}, obj.strfield)
%                                 obj.variants(i).(infofd{ii}) = infodata{ii};
%                             else
%                                 obj.variants(i).(infofd{ii}) = eval(sprintf('[%s]'';',infodata{ii}));
%                             end
%                         end
%                     end
%                 end
%             end            
%         end


%         function parseFormatData_old(obj, header)
%             h = textscan(strrep(header,'#',''), '%s');
%             h = h{1}; 
%             
%             formatinfo = regexp(headers, 'FORMAT=<ID=(\w+),Number=([\w\.]+),Type=(\w+),', 'tokens', 'once');
%             formatinfo = formatinfo(~cellfun(@isempty, formatinfo));
%             formatinfo = vertcat(formatinfo{:});
%             
%             idx = find(strcmp(h, 'FORMAT'));
%             nsample = length(h(idx+1:end));
%             
%             obj.formatData = cell(numvar, nsample, size(obj.formatfield,1));            
%             [~, orderidx] = cellfun(@(x) ismember(x, obj.formatfield(:,2)), obj.formatRawData(:,1), 'unif', 0);            
%             obj.formatRawData(:,1) = [];            
%             vi(:,1) = [];
%             obj.validCalls = sparse(vi);
%             for i = 1:numvar
%                 tmpdata = obj.formatRawData(i, vi(i,:));
%                 tmpdata = vertcat(tmpdata{:});
%                 obj.formatData(i, vi(i,:), orderidx{i}) = tmpdata;
%             end            
%             obj.formatData(cellfun(@isempty, obj.formatData)) = {''};
%             obj.formatData(:, :, obj.formatNum==1) = str2double_fast(obj.formatData(:, :, obj.formatNum==1));            
% %             obj.formatData(:, :, obj.formatNum > 1) = cellfun(@str2num, ...
% %                 obj.formatData(:, :, obj.formatNum > 1), 'unif', 0);                        
%             obj.formatRawData = {};
%         end
% 
%  function cleanup(obj, hasvalue)
%             fds = fieldnames(obj.variants);
%             for i = 1:length(fds)
%                 if ~hasvalue(i)
%                     obj.variants = rmfield(obj.variants, fds{i});
%                     obj.desc = rmfield(obj.desc, fds{i});
%                 end
%             end
%         end
%         
%         

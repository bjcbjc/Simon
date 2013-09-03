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
        function obj = VCF(fn, parseFormat, pyformatfiledir, goodVCF)
            if nargin >= 2
                obj.parseFormat = parseFormat;
            end
            if nargin == 3
                obj.pyformatfiledir = pyformatfiledir;
            end           
            if nargin < 4
                goodVCF = false;
            end
            
            obj.filename = fn;
            
            if goodVCF
                obj = obj.readGoodVCF();
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
            
            numeric = ismember(infotable(:,3), {'Float', 'Integer', 'Flag'}) ...
                & strcmp(infotable(:,2), '1');
                
            obj.attrName_mtx = [obj.attrName_mtx; infotable(numeric, 1)];
            obj.attrName_cell = [obj.attrName_cell; infotable(~numeric, 1)];
            
            formattable = obj.parseMeta(header, {'FORMAT'});
            obj.desc = [obj.desc; formattable];
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
        
        function obj = readGoodVCF(obj)
            %FORMAT string is the same for all the records
            [header, varExample] = obj.obtainHeader();
            obj.parseHeader(header);
            fields = textscan(varExample, '%s', 'delimiter', sprintf('\t'));
            fields = fields{1};
            nsample = length(fields) - 9;
            obj.sample = fields(10:end);
            
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
                    strcmp(obj.desc(formatidx(i), 2), '1')
                    sampleformatstr = strcat(sampleformatstr, '%f');
                else
                    sampleformatstr = strcat(sampleformatstr, '%s');
                end
            end
            varformatstr = strcat(varformatstr, repmat(sampleformatstr, 1, nsample));
            fprintf('read file\n');
            tic;
            f = fopen(obj.filename);
            d = textscan(f, varformatstr, 'delimiter',sprintf('\t:'), 'headerlines', length(header), 'treatasempty', '.');
            fclose(f);
            toc;
            
            nvariant = length(d{1});
            nmtx = length(obj.attrName_mtx);
            ncell = length(obj.attrName_cell);
            obj.variantAttr_mtx = NaN(nvariant, nmtx);
            obj.variantAttr_cell = repmat({''}, nvariant, ncell);
            fprintf('fill in basic info\n');
            tic;
            colheader = textscan(header{end}(2:end), '%s', 'delimiter', sprintf('\t'));
            colheader = colheader{1};
            colheader = colheader(:)';
            [~, colidx] = ismember(colheader, obj.attrName_mtx);
            obj.variantAttr_mtx(:, colidx(colidx~=0)) = cell2mat( ...
                arrayfun(@(idx) d{idx}, find(colidx), 'unif', 0) );
            [~, colidx] = ismember(colheader, obj.attrName_cell);
            for i = 1:length(colidx)
                if colidx(i) ~= 0
                    obj.variantAttr_cell(:, colidx(i)) = d{i};
                end
            end
            toc; 
            
            colidx = strcmp(obj.attrName_cell, 'ALT');
            obj.maxnalt = max( cellfun(@(x) length(strfind(x,',')), obj.variantAttr_cell(:, colidx)) ) + 1;
            
            fprintf('parse format data\n');
            tic;
            %format data
            for i = 1:nformatfields
                obj.formatData.(formatfields{i}) = ...
                    d( 8+nformatfields+i:nformatfields:end );
                obj.formatData.(formatfields{i}) = horzcat(obj.formatData.(formatfields{i}){:});
            end
            toc;
            
            
            %INFO
            d = d{8};
            
            %flag-info
            fprintf('fill in flag info\n');
            tic;
            flags = intersect(obj.desc(strcmp(obj.desc(:,3), 'Flag'), 1), obj.attrName_mtx);
            for i = 1:length(flags)
                colidx = strcmp(obj.attrName_mtx, flags{i});
                obj.variantAttr_mtx(:, colidx) = double(...
                    ~cellfun(@isempty, strfind(d, flags{i})) );
            end
            toc;
            
            fprintf('index info fields\n');
            tic;
            reptb = [strcat(obj.attrName_mtx, '='), strcat(arrayfun(@num2str, (1:nmtx)', 'unif',0),'=')];
            reptb = [reptb; [strcat(obj.attrName_cell, '='), strcat(arrayfun(@num2str, (1:ncell)'+nmtx,'unif',0), '=')]];
            for i = 1:size(reptb, 1)
                d = strrep(d, reptb{i,1}, reptb{i,2});
            end
            toc;
            
            fprintf('parse info data\n');
            tic;
            infos = regexp(d, '([\w\.]+)=([\w\+\-\.\,]+)', 'tokens');
            infos = cellfun(@(x) vertcat(x{:}), infos, 'unif', 0);
            nval = cellfun(@(x) size(x,1), infos);
            infos = vertcat(infos{:}); %flat the entire info data
            toc;
            
            fprintf('fill info data overhead\n');
            tic;            
            varidx = cell2mat(arrayfun(@(i) repmat(i, nval(i), 1), (1:nvariant)', 'unif', 0));
            infoidx = str2double_fast(infos(:,1));
            infos = infos(:,2); 
            toc;
            
            fprintf('fill in numeric info data\n');
            tic;
            %num-info
            curidx = infoidx <= nmtx;
            obj.variantAttr_mtx( sub2ind([nvariant, nmtx], varidx(curidx), infoidx(curidx)) ) = ...
                str2double_fast(infos( curidx ) );
            toc;
            
            fprintf('fill in non-numeric info data\n');
            tic;
            %cell-info
            curidx = ~curidx;
            obj.variantAttr_cell( sub2ind([nvariant, ncell], varidx(curidx), infoidx(curidx)-nmtx) ) = infos( curidx ) ;
            toc;
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

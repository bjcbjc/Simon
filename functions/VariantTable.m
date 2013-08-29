classdef VariantTable < handle
    
    properties
        filename;
        data;
        header;
        samplebit;
        callercount;
        SGD;
        RELATED;
    end

    
    
    methods 
        function obj = VariantTable(fn)
            t = parseText(fn,'nrowname',0,'ncolname',1, 'ncol',25);
            obj.filename = fn;
            obj.data = t.text;
            obj.header = t.colname;
            obj.samplebit = false(size(obj.data,1),5,2);
            [~,i] = ismember({'IP-5','SFP-5-A','SFP-5-B','FP-5-A','FP-5-B'}, obj.header);
            obj.samplebit(:,:,1) = logical(cellfun(@str2double, obj.data(:,i)));
            [~,i] = ismember({'IP-6','SFP-6-A','SFP-6-B','FP-6-A','FP-6-B'}, obj.header);
            obj.samplebit(:,:,2) = logical(cellfun(@str2double, obj.data(:,i)));
            [~,i] = ismember({'#GATK','#SAM','#FB'}, obj.header);
            obj.callercount = cellfun(@str2double, obj.data(:,i));
            load('/users/bjc/Desktop/Projects/DATA/Yeast/SGD/SGD.mat');
            obj.SGD = SGD;
            load RELATED;
            obj.RELATED = RELATED;
        end
        function d = getData(obj, columns)
            if ischar(columns)
                columns = {columns};
            end
            [~,i] = ismember(columns, obj.header);
            d = obj.data(:,i);
        end        
        function filterSingleCall(obj, drug)
            %for caffeine and rapamycin
            linecount = sum(logical(squeeze(sum(sum(obj.samplebit,1), 2))));
            withinlinecount = squeeze(sum(obj.samplebit,2)); % nrows x 2 lines
            valid = false(size(obj.data,1),1);
            valid(sum(obj.callercount > 0,2) >= 2) = true;
            valid(linecount == 2) = true;
            valid(any(withinlinecount >= 2, 2)) = true;
            load('/users/bjc/Desktop/Projects/DATA/Yeast/SGD/INT.mat');
            
            %check if within list interactions
            orfcol = strcmp(obj.header, 'ORF');            
            vi = ~strcmp(INT.ORF(:,1), INT.ORF(:,2));            
            vi = vi & ismember(INT.ORF(:,1), obj.data(:,orfcol))  ...
                & ismember(INT.ORF(:,2), obj.data(:,orfcol));
            fds = fieldnames(INT);
            for fi = 1:length(fds)
                INT.(fds{fi}) = INT.(fds{fi})(vi,:);
            end
            valid(ismember(obj.data(:, orfcol), INT.ORF(:,1)) | ...
                ismember(obj.data(:, orfcol), INT.ORF(:,2))) = true;
            
            %check if have drug annotation
            drug = obj.drugtranslate(drug);
            valid(ismember(obj.data(:,orfcol), obj.RELATED.([drug 'orf']))) = true;            
            if ~strcmp(drug, 'rapa')
                valid(ismember(obj.data(:,orfcol), obj.RELATED.dnaorf)) = true;
                valid(ismember(obj.data(:,orfcol), obj.RELATED.chemodrugorf)) = true;
                valid(ismember(obj.data(:,orfcol), obj.RELATED.uvxrayorf)) = true;
            end
            
            %apply filter
            obj.data(~valid, :) = [];
            obj.samplebit(~valid, :, :) = [];
            obj.callercount(~valid, :) = [];            
        end
        function addLineCounts(obj, colname)            
            linecount = sum(logical(squeeze(sum(obj.samplebit,2))),2);
            obj.header{end+1} = colname;
            obj.data(:,end+1) = cellfun(@num2str, num2cell(linecount), 'uniformoutput',false);
        end
        function addGeneDesc(obj)            
            orfcol = strcmp(obj.header, 'ORF');
            [~,i] = ismember(obj.data(:,orfcol), obj.SGD.orf);
            obj.header{end+1} = 'GENE DESC';
            obj.data(:,end+1) = {''};
            obj.data(i~=0,end) = obj.SGD.desc(i(i~=0));
        end
        function addInteraction(obj, withinlist, colname, drug) 
            drug = obj.drugtranslate(drug);
            load('/users/bjc/Desktop/Projects/DATA/Yeast/SGD/INT.mat');
            orfcol = strcmp(obj.header, 'ORF');
            %avoid self interaction
            vi = ~strcmp(INT.ORF(:,1), INT.ORF(:,2));
            if withinlist
                vi = vi & ismember(INT.ORF(:,1), obj.data(:,orfcol))  ...
                     & ismember(INT.ORF(:,2), obj.data(:,orfcol));
            else %limit to related genes
                drugfd = [drug 'orf'];
                vi = vi & (ismember(INT.ORF(:,1), obj.data(:,orfcol)) ...
                    | ismember(INT.ORF(:,2), obj.data(:,orfcol)) );
                vi = vi & (ismember(INT.ORF(:,1), obj.RELATED.(drugfd)) ...
                    | ismember(INT.ORF(:,2), obj.RELATED.(drugfd)) );
            end
            addColumn = cell(size(obj.data,1),1);
            edgeDensity = cell(size(obj.data,1),1);
            edgeDensity(:) = {''};
            fds = fieldnames(INT);
            for fi = 1:length(fds)
                INT.(fds{fi}) = INT.(fds{fi})(vi,:);
            end
            for i = 1:size(obj.data,1)
                uorf = union(INT.ORF( strcmp(INT.ORF(:,1),obj.data{i,orfcol}), 2), ...
                    INT.ORF( strcmp(INT.ORF(:,2),obj.data{i,orfcol}), 1) );
                itype = false(length(uorf), 2);
                
                ri = strcmp(obj.data(i,orfcol), INT.ORF(:,1));
                [~,ui] = ismember(INT.ORF(ri,2), uorf);
                itype(ui,1) = itype(ui,1) | strcmp(INT.type(ri),'genetic interactions');
                itype(ui,2) = itype(ui,2) | strcmp(INT.type(ri),'physical interactions');
                
                ri = strcmp(obj.data(i,orfcol), INT.ORF(:,2));
                [~,ui] = ismember(INT.ORF(ri,1), uorf);
                itype(ui,1) = itype(ui,1) | strcmp(INT.type(ri),'genetic interactions');
                itype(ui,2) = itype(ui,2) | strcmp(INT.type(ri),'physical interactions');
                gene = translategene(obj.SGD, uorf);
                [~, sgi] = sort(sum(itype, 2), 'descend');
                gene = gene(sgi);
                itype = itype(sgi, :);
                fieldstr = '';
                ngene = length(gene);
                for ii = 1:ngene
                    if all(itype(ii,:))
                        tailnote = '';
                    elseif itype(ii,1)
                        tailnote = '(g)';
                    else
                        tailnote = '(p)';
                    end
                    if ii < length(gene)
                        fieldstr = sprintf('%s%s%s; ',fieldstr, gene{ii}, tailnote);
                    else
                        fieldstr = sprintf('%s%s%s',fieldstr, gene{ii}, tailnote);
                    end
                end
                addColumn{i} = fieldstr;
                if withinlist
                    tmpmtx = zeros(ngene,ngene);
                    vi = ismember(INT.ORF(:,1), uorf) & ismember(INT.ORF(:,2),uorf);
                    [~,gi] = ismember(INT.ORF(vi,1),uorf);
                    [~,gj] = ismember(INT.ORF(vi,2),uorf);
                    tmpmtx(sub2ind([ngene, ngene], min(gi,gj), max(gi,gj))) = 1;
                    if ngene > 0                        
                        edgeDensity{i} = sprintf('%0.2g', nnz(tmpmtx) / (0.5*ngene*(ngene-1)));
                    end
                end
            end
            obj.header{end+1} = colname;
            obj.data(:, end+1) = addColumn;
            if withinlist
                obj.header{end+1} = 'INT DENSITY';
                obj.data(:, end+1) = edgeDensity;
            end
        end
        function addDrugAnnotation(obj, colname, drug)
            drug = obj.drugtranslate(drug);
            addColumn = cell(size(obj.data,1),1);
            addColumn(:) = {'No'};
            orfcol = strcmp(obj.header, 'ORF');
            addColumn(ismember(obj.data(:,orfcol), obj.RELATED.([drug 'orf']))) = {'Yes'};
            obj.header{end+1} = colname;
            obj.data(:, end+1) = addColumn;
        end
        function addGO(obj, colname, mGOsize, drug)
            drug = obj.drugtranslate(drug);
            load('/users/bjc/Desktop/Projects/DATA/Yeast/GO/GO.mat');
            %load('/users/bjc/Desktop/Projects/DATA/Yeast/GO/GO_forenrichment.mat');
            s = sum(GO.mtx);
            vci = s > mGOsize & s < 700;
            GO.mtx = GO.mtx(:, vci);
            GO.cond = GO.cond(vci);
            GO.spaceidx = GO.spaceidx(vci);
            GO.desc = GO.desc(vci);
            orfcol = strcmp(obj.header, 'ORF');
            [~, gi] = ismember(obj.data(:, orfcol), GO.orf);
            i = unique(gi(gi~=0));
            ci = find(any(GO.mtx(i,:)));
            s = sum(GO.mtx(i, ci));
            ci = ci(s>=2);%only consider annotation with >=2 of the listed genes
            s = s(s>=2);
            sa = sum(GO.mtx(:,ci));
            N = length(GO.orf);
            M = length(unique(obj.data(:, orfcol)));
            p = NaN(1,length(ci));
            for cii = 1:length(ci)
                p(cii) = sum(hygepdf(s(cii):M, N, M,sa(cii)));
            end
            if strcmpi(drug, 'nqo')
                [~,pc] = FDR(p, 0.05);
                if sum(p<=pc) > 0
                    keep = s > 50 | ( s > 2 & p <= pc);
                else
                    keep = s > 50 | p <= 0.01;
                end
            else
                keep = s > 3 | p < 0.01;                
            end
            fprintf('%d annotations obtained\n', full(sum(keep)));
            ci = ci( keep);
            p = p(keep);
            s = s(keep);
            addColumn = cell(size(obj.data,1),1);
            addColumn(:) = {''};
            for i = 1:length(addColumn)
                if gi(i) ~= 0
                    cii = find(GO.mtx(gi(i), ci)==1);
                    note = '';
                    for j = 1:length(cii)
                        note = sprintf('%s%s(%0.2g); ',note, GO.cond{ci(cii(j))},p(cii(j)));
                    end
                    addColumn{i} = note;
                end                
            end
            obj.header{end+1} = colname;
            obj.data(:, end+1) = addColumn;
        end        
        function drug = drugtranslate(obj, drug)
            if strcmpi(drug, '4NQO') || strcmpi(drug, 'NQO')
                drug = 'nqo';
            elseif strcmpi(drug, 'caffeine') || strcmpi(drug, 'caf')
                drug = 'caf';
            elseif strcmpi(drug, 'rapamycin') || strcmpi(drug, 'rapa')
                drug = 'rapa';
            elseif strcmpi(drug, 'dna')
                drug = 'dna';
            elseif strcmpi(drug, 'chemo')
                drug = 'chemodrug';
            elseif strcmpi(drug, 'uvxray')
                drug = 'uvxray';
            else
                error('unknown drug');
            end  
        end
        function outputTable(obj, fn, columnorder)
            [r c] = size(obj.data);
            if nargin < 3,
                columnorder = 1:c;
            else
                if any(columnorder > c)
                    error('column index out of range\n');
                end
            end
            
            f = fopen(fn, 'w');
            
            fprintf(f, '%s', obj.header{columnorder(1)});
            for ci = 2:length(columnorder)
                fprintf(f, '\t%s', obj.header{columnorder(ci)});
            end
            fprintf(f, '\n');
            for ri = 1:r
                fprintf(f, '%s', obj.data{ri,1});
                for ci = 2:length(columnorder)
                    fprintf(f, '\t%s', obj.data{ri,columnorder(ci)});
                end
                fprintf(f, '\n');
            end
            fclose(f);
        end
    end
    
end
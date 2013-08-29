classdef Patient < dynamicprops
    properties
        patient;
        sample;
        minRead;
        maxDepth;
        minMapQ;
        loc;        
        ref;
        alt;
        alt_original;
        refCountDNA;
        altCountDNA;
        ntlabels;
        ALLCountDNA;
        ALLCountRNA;
        indelsDNA;
        indelsRNA;
        refCountRNA;
        altCountRNA;
        variantAttr_cell_name;
        variantAttr_cell;
        variantAttr_mtx_name;
        variantAttr_mtx;
        formatData;
        vcfSample;
        vcfDesc;
        annt;
            
        totalReadDNA;
        totalReadRNA; 
        refDnaAF;
        altDnaAF;
        refRnaAF;
        altRnaAF;
        RAratioDNA;
        RAratioRNA;
        distToCen;
        
        %DESeq data
        DEpvalAdj;
        DElabel;
        
        %preset of variants 
        presetidx;
        minReadDNA = 11;
        minReadRNA = 11;
        maxReadDNA = 1000;
        
        callerColidx;        
        
        cdict;
        
%         caller = {'mutect', 'sniper', 'varscan', 'mutect-sniper', 'mutect-varscan', ...
%             'sniper-varscan', 'mutect-sniper-varscan'};
%         callerColor = [0.8, 0.1, 0.1; 0.1, 0.6, 0.1; 0.1, 0.1, 0.8; 0.8, 0.6, 0.2; 0.8, 0.3, 0.6; 0.2, 0.8, 0.7; 0,0,0];
    end
    
    properties (Access = private)
        %logical index for current selected variants        
        subsetdesc;
        subsetidx;
    end
    
    methods
        function obj = Patient(patientSampleInfo, DESeqData, CENTROMERE, varargin)
            obj.readData(patientSampleInfo, varargin{:});            
            
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

            obj.addStat(DESeqData, CENTROMERE);
            
            ckey = {'mutect', 'sniper', 'varscan', ...
                'mutect-sniper', 'mutect-varscan', 'sniper-varscan', 'mutect-sniper-varscan', ...
                'DE', 'non DE'};
            cmap = [0.8, 0.1, 0.1; 0.1, 0.6, 0.1; 0.1, 0.1, 0.8; ...
                0.8, 0.6, 0.2; 0.8, 0.3, 0.6; 0.2, 0.8, 0.7; 0.8,0.8,0.8; ...
                1, 0.1, 0.1; 0.1, 0.1, 1];
            obj.cdict = ColorDict(cmap, ckey);
        end
       
        function selectVariant(obj, selection)
            % choose subset of variants
            % selection string example: 'expressedBoth & ~(multAlt | gatkindel |
            % chrMT)', each represents a logical index in obj.presetidx
            % selection can otherwise be a logical index to be directly
            % assigned to subsetidx
            
            if ischar(selection)
                obj.subsetdesc = selection;
                obj.subsetidx = eval(regexprep(selection, '([\w_]+)', 'obj.presetidx.$1'));
            elseif islogical(selection)
                if length(selection) == length(obj.ref)
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
                
        function availSet = availabelPreset(obj)
            availSet = fieldnames(obj.presetidx);
        end
        
        function [idx, desc] = currentSet(obj)
            if all(obj.subsetidx)
                obj.subsetdesc = 'all';
            end
            idx = obj.subsetidx;
            desc = obj.subsetdesc;
        end
                                
        function fig = plotRAratio(obj, varargin)
            para.fig = '';
            para.save = 'REF-ALT ratio, normal vs tumor';
                        
            para = assignpara(para, varargin{:});
            fig = obj.setupFigure(para.fig);
            
            datatype = {'DNA', 'RNA'};
            for dataidx = 1:length(datatype)
                subplot(2,2, (dataidx-1)*2+1);
                h = gscatter(obj.(['RAratio' datatype{dataidx}])(obj.subsetidx, 1), ...
                    obj.(['RAratio' datatype{dataidx}])(obj.subsetidx, 2), ...
                    obj.variantAttr_cell(obj.subsetidx, obj.callerColidx));
                set(h, 'marker', 'o', 'markersize', 10, 'linewidth',1);
%                 [~, calleridx] = ismember(get(h, 'displayname'), obj.caller);
                for hi = 1:length(h)
                    set(h(hi), 'color', obj.cdict.color(get(h(hi), 'displayname')));
%                     set(h(hi), 'color', obj.callerColor(calleridx(hi),:));
                end
                set(gca, 'fontsize', 14);
                xlabel(sprintf('R/A in normal %s',datatype{dataidx}), 'fontsize', 14);
                ylabel(sprintf('R/A in tumor %s',datatype{dataidx}), 'fontsize', 14);
                lim = max(xlim, ylim);
                lim(1) = 0;
                xlim(lim);
                ylim(lim);
                axis square
                if dataidx == 1
                    hleg = legend('location', 'NO', 'orientation', 'horizontal');
                    set(hleg, 'position', [0.014, 0.95, 0.97, 0.049], 'fontsize', 14);
                else
                    legend('off');
                end
                
                subplot(2,2,dataidx*2);
                h = gscatter(obj.(['RAratio' datatype{dataidx}])(obj.subsetidx, 1), ...
                    obj.(['RAratio' datatype{dataidx}])(obj.subsetidx, 2), ...
                    obj.variantAttr_cell(obj.subsetidx, obj.callerColidx));
%                 [~, calleridx] = ismember(get(h, 'displayname'), obj.caller);
                for hi = 1:length(h)
                    set(h(hi), 'color', obj.cdict.color(get(h(hi), 'displayname')));
%                     set(h(hi), 'color', obj.callerColor(calleridx(hi),:));
                end
                set(h, 'marker', 'o', 'markersize', 12, 'linewidth', 1.2);
                set(gca, 'fontsize', 14);
                xlabel(sprintf('R/A in normal %s',datatype{dataidx}), 'fontsize', 14);
                ylabel(sprintf('R/A in tumor %s',datatype{dataidx}), 'fontsize', 14);
                xlim( [0, 100])
                ylim(xlim);
                axis square
                legend('off');
            end
            
            if para.save
                Patient.saveFigure('fig', fig, 'savedir', para.savedir, 'save', para.save, 'sampleid', obj.patient);                
            end
        end

        function [hscatter, fig] = plotScatterALTAF(obj, varargin)
            para = obj.figDefaultPara('scatter');            
            para.save = 'major ALT freq';
            
            para = assignpara(para, varargin{:});
            if isempty(para.colorByTitle) && ~isempty(para.colorBy)
                if ischar(para.colorBy)
                    para.colorByTile = para.colorBy;
                else
                    para.colorByTitle = 'data';
                end
            elseif isempty(para.colorByTitle) && isempty(para.colorBy)
                para.saveAppendColorBy = false;
            end
            parastr = paraPairArray(para);
            
            fig = obj.setupFigure(para.fig);
            figure(fig);
            
            
            samplelabel = {'normal', 'tumor'};            
            xfd = {'altDnaAF', 'altDnaAF', 'altDnaAF', 'altRnaAF'};
            yfd = {'altRnaAF', 'altRnaAF', 'altDnaAF', 'altRnaAF'};
            xsample = [1, 2, 1, 1];
            ysample = [1, 2, 2, 2];
            
            curcolormap = colormap();            
            
            for plotidx = 1:length(xfd)
                subplot(2,2, plotidx);
                xlim([-0.01 1.01]);
                ylim(xlim);
                [hscatter, plottype] = obj.plotScatter(obj.(xfd{plotidx})(obj.subsetidx,xsample(plotidx),1), ...
                    obj.(yfd{plotidx})(obj.subsetidx,ysample(plotidx),1), ...
                    'passpara', parastr{:});

                if xsample(plotidx) ~= ysample(plotidx)
                    xlabel(samplelabel{xsample(plotidx)}, 'fontsize', 14);
                    ylabel(samplelabel{ysample(plotidx)}, 'fontsize', 14);
                    title(upper(regexprep(xfd{plotidx}, 'alt(\w{3})AF', '$1')), 'fontsize', 14);
                else
                    xlabel(upper(regexprep(xfd{plotidx}, 'alt(\w{3})AF', '$1')), 'fontsize', 14);
                    ylabel(upper(regexprep(yfd{plotidx}, 'alt(\w{3})AF', '$1')), 'fontsize', 14);
                    title(samplelabel{xsample(plotidx)}, 'fontsize', 14);
                end
                if plotidx == 3 && strcmpi(plottype, 'gpatchscatter')
                    hleg = legend('location', 'NO', 'orientation', 'horizontal');
                    set(hleg, 'position', [0.37, 0, 0.3, 0.04], 'fontsize', 14);
                end
            end
                        
            supertitle = sprintf('major ALT allele frequency');
            if ~isempty(para.colorBy)
                supertitle = sprintf('%s, colored by %s',supertitle, para.colorByTitle);
            end
            suptitle(supertitle);
            if ~isempty(para.save)
                para.sampleid = obj.patient;
                if ~isempty(strfind(plottype, 'patch'))
                    para.format = 'svg';                    
                else
                    para.format = 'eps';                    
                end
                Patient.saveFigure(para);                
                colormap(curcolormap);
            end
        end
        
        function [hscatter, fig] = plotScatterDnaRnaAF(obj, varargin)
            para = obj.figDefaultPara('scatter');                     
            para.save = 'DNA vs RNA allele freq';
            
            para = assignpara(para, varargin{:});
            if isempty(para.colorByTitle) && ~isempty(para.colorBy) 
                if ischar(para.colorBy)
                    para.colorByTile = para.colorBy;
                else
                    para.colorByTitle = 'data';
                end
            elseif isempty(para.colorByTitle) && isempty(para.colorBy)
                para.saveAppendColorBy = false;
            end
            parastr = paraPairArray(para);
            
            fig = obj.setupFigure(para.fig);
            figure(fig);
            
            allele = {'ref', 'alt'};
            
            curcolormap = colormap();            
            
            for aidx = 1:length(allele)
                for sampi = 1:length(obj.sample)
                    xfd = [allele{aidx} 'DnaAF'];
                    yfd = [allele{aidx}, 'RnaAF'];
                    subplot(2,2, (aidx-1)*2+sampi);
                    xlim([-0.01,1.01]);
                    ylim(xlim);
                    [hscatter, plottype] = obj.plotScatter(obj.(xfd)(obj.subsetidx,sampi,1), ...
                        obj.(yfd)(obj.subsetidx,sampi,1), ...
                        'passpara', parastr{:});

                    if aidx == 2 && sampi == 1 && strcmpi(plottype, 'gpatchscatter')
                        h = legend('orientation', 'horizontal');
                        set(h, 'position', [0.37, 0, 0.3, 0.04], 'fontsize', 14);
                    end
                    xlabel('DNA', 'fontsize', 14);
                    ylabel('RNA', 'fontsize', 14);
                    title(sprintf('%s freq, %s', upper(allele{aidx}), obj.sample{sampi}), 'fontsize', 14);
                end
            end
                        
            supertitle = sprintf('DNA (X-axis) vs. RNA (Y-axis) allele frequency');
            if ~isempty(para.colorBy)
                supertitle = sprintf('%s, colored by %s',supertitle, para.colorByTitle);
            end
            suptitle(supertitle);
            if ~isempty(para.save)
                para.sampleid = obj.patient;
                if ~isempty(strfind(plottype, 'patch'))
                    para.format = 'svg';                    
                else
                    para.format = 'eps';                    
                end
                Patient.saveFigure(para);                                
                colormap(curcolormap);
            end
        end
        
        function fig = plotScatterNormalTumorDnaAF(obj, varargin)
            para.fig = [];            
            para.save = 'DNA AF normal vs tumor with genotype';            
            para.plotGT = {'0-0/1', '0/0-0/1', '0/0-1/1', '0/1-0/1', ...
                '1/1-0/0', '1/1-0/1'};
            
            para = assignpara(para, varargin{:});
            fig = obj.setupFigure(para.fig);
            
            GT = strcat(obj.formatData.GT(:,1), '-', obj.formatData.GT(:,2));
            for i = 1:min(6, length(para.plotGT))
                idx = obj.subsetidx & strcmp(GT, para.plotGT{i});
                subplot(2, 3, i);
                gscatter(obj.altDnaAF(idx, 1, 1), obj.altDnaAF(idx, 2, 1), ...
                    obj.variantAttr_cell(idx, obj.callerColidx));
                xlabel(sprintf('ALT AF in %s', obj.sample{1}), 'fontsize', 14);
                ylabel(sprintf('ALT AF in %s', obj.sample{2}), 'fontsize', 14);                
                axis square
                xlim([0, 1]);
                ylim(xlim);
                title(para.plotGT{i}, 'fontsize', 14);
                if ismember(i, 2:4)
                    h = legend('location', 'southeast');
                    if i == 2
                        set(h,'fontsize', 6);
                    end
                else
                    legend('location', 'northwest');
                end
                hold on; plot(xlim, xlim, 'k-', 'linewidth', 1.5); hold off
            end
            if ~isempty(para.save)
                Patient.saveFigure('savedir', para.savedir, 'save', para.save, 'sampleid', obj.patient, 'fig', fig);                
            end
        end
        
        function [hscatter, fig] = plotScatterNormalTumorCount(obj, varargin)        
            para = obj.figDefaultPara('scatter');                                    
            para.save = 'raw count, normal vs tumor';                                    
            
            para = assignpara(para, varargin{:});
            
            if isempty(para.colorByTitle) && ~isempty(para.colorBy) 
                if ischar(para.colorBy)
                    para.colorByTile = para.colorBy;
                else
                    para.colorByTitle = 'data';
                end
            elseif isempty(para.colorByTitle) && isempty(para.colorBy)
                para.saveAppendColorBy = false;
            end
            parastr = paraPairArray(para);
            
            fig = obj.setupFigure(para.fig);            
            figure(fig);
            curcolormap = colormap();
            
            datatype = {'DNA', 'RNA'};
            datafield = {'refCount', 'altCount', 'totalRead'};            
            
            for dataidx = 1:length(datatype)
                for fdidx = 1:length(datafield)
                    subplot(2,3, (dataidx-1)*3+fdidx );                    
                    fd = strcat(datafield{fdidx}, datatype{dataidx});
                    [hscatter, plottype] = obj.plotScatter(log10(obj.(fd)(obj.subsetidx,1,1)+1), ...
                        log10(obj.(fd)(obj.subsetidx,2,1)+1), ...
                        'passpara', parastr{:});

                    xlabel('normal', 'fontsize', 14);
                    ylabel('tumor', 'fontsize', 14);
                    title(sprintf('log10 %s in %s', datafield{fdidx}, datatype{dataidx}), 'fontsize',14);
                    %ylim(xlim);
                    axis square

                    if dataidx == 2 && fdidx == 1 && strcmpi(plottype, 'gpatchscatter') 
                        h = legend('orientation', 'horizontal');
                        set(h, 'position', [0.37, 0, 0.3, 0.04], 'fontsize', 14);
                    end
                    hold on; plot(xlim, xlim, 'k-', 'linewidth', 1.5); hold off;                
                end
            end
            
            supertitle = sprintf('Normal (X-axis) vs. Tumor (Y-axis) read counts');
            if ~isempty(para.colorBy)
                supertitle = sprintf('%s, colored by %s',supertitle, para.colorByTitle);
            end
            suptitle(supertitle);
            if ~isempty(para.save)
                para.sampleid = obj.patient;
                if ~isempty(strfind(plottype, 'patch'))
                    para.format = 'svg';                    
                else
                    para.format = 'eps';                    
                end
                Patient.saveFigure(para);                
                colormap(curcolormap);
            end
        end        
        
        function fig = plotHistMajorALTAF(obj, varargin)
            para.fig = [];            
            para.save = 'hist major ALT freq by callers';
            para.barcolors = [0.6, 0.6, 1; 1, 0.6, 0.6; 0.3, 0.3, 0.9; 0.9, 0.3, 0.3];
            para.legendlabel = {'DNA-normal', 'DNA-tumor', 'RNA-normal', 'RNA-tumor'};
            
            para = assignpara(para, varargin{:});
            fig = obj.setupFigure(para.fig);
            figure(fig);
                        
            for i = 1:length(obj.caller)
                idx = obj.subsetidx & strcmp(obj.variantAttr_cell(:, obj.callerColidx), obj.caller{i});
                subplot(3,3,i);
                if sum(idx) == 1
                    tmp = [obj.altDnaAF(idx,:,1), obj.altRnaAF(idx,:,1)];
                    for si = 1:4
                        hold on;
                        plot([tmp(si), tmp(si)], [0, 1], 'color', para.barcolors(si,:), 'linewidth', 4);
                        hold off;
                    end
                    ylim([0 2]);
                    clear tmp
                else
                    [h, x] = hist([obj.altDnaAF(idx,:,1), obj.altRnaAF(idx,:,1)], 0:0.2:1);
                    bar(x, h, 'hist');
                    ylim([0, max(h(:))+1]);
                    xlim([-0.1 1.1]);
                    set(gca, 'box', 'off');
                end
                xlabel('ALT allele freq', 'fontsize', 12);
                ylabel('number of variants', 'fontsize', 12);
                title(obj.caller{i}, 'fontsize', 14);
                if i == length(obj.caller)
                    h = legend(para.legendlabel);
                    set(h, 'position',  [0.4273    0.1475    0.1235    0.1337]);
                end
                if sum(idx) > 1
                    h = findobj(gca, 'Type', 'patch');
                    h = h(4:-1:1);
                    for si = 1:4
                        set(h(si), 'edgecolor', 'none', 'facecolor', para.barcolors(si,:));
                    end
                end
            end
            if ~isempty(para.save)
                Patient.saveFigure('fig', fig, 'savedir', para.savedir, 'save', para.save, ...
                    'sampleid', obj.patient, 'format', 'eps');                
            end
        end
                
                
        function [h, plottype] = plotScatter(obj, x, y, varargin)
            para = obj.figDefaultPara('scatter');            
            
            para = assignpara(para, varargin{:});
                                    
            if isempty(para.colorBy)                                
                h = para.func(x, y);          
                plottype = func2str(para.func);
                return
            elseif ischar(para.colorBy)
                if isprop(obj, para.colorBy)
                    para.colorBy = obj.(para.colorBy)(obj.currentSet);
                else
                    fprintf('unknown properties %s for colors\n', para.colorBy);
                end
            end
            
            if ~isempty(para.valueOrder)
                if ~isempty(para.colorBy) && isscalar(para.colorBy)
                    [para.colorBy, si] = sort(para.colorBy, para.valueOrder);
                    x = x(si);
                    y = y(si);
                end
            end
            
            if iscellstr(para.colorBy)
                if iscellstr(para.valueOrder)
                    groupcolor = obj.cdict.color(para.valueOrder);
                else
                    groupcolor = obj.cdict.color(unique(para.colorBy));
                end
                h = gPatchScatter(x, y, para.colorBy, 'facecolor', groupcolor, ...
                    'renderOrder', para.valueOrder);
%                 h = patchScatter(x, y, 'facecolor', obj.cdict.color(para.colorBy));
                plottype = 'gpatchscatter';
            else
                if ~isempty(para.clim)
                    para.colorBy(para.colorBy<para.clim(1)) = para.clim(1);
                    para.colorBy(para.colorBy>para.clim(1)) = para.clim(2);
                elseif ~isempty(para.climPrctile)
                    dlim = prctile(para.colorBy, para.climPrctile);
                    para.colorBy(para.colorBy<dlim(1)) = dlim(1);
                    para.colorBy(para.colorBy>dlim(2)) = dlim(2);
                end
                h = scatter(x, y, para.size, para.colorBy);
                colormap(para.colormap);                
                freezeColors;
                cbfreeze(colorbar);
                plottype = 'scatter';
            end
        end                
        
        
        function [table, header] = getCurSelVarAnnotation(obj, fields, varidx)
            if nargin < 3, varidx = obj.subsetidx; end
            if nargin < 2, fields = []; end
            if isempty(obj.annt)
                fprintf('No annotation is available\n');
                table = {};
                header = {};
                return
            end
            fds = fieldnames(obj.annt);
            if ~isempty(fields)
                [~, fdidx] = ismember(fields, fds);
                fields(fdidx==0) = [];
                fdidx(fdidx==0) = [];
            end
            if isempty(fields)
                fdidx = 1:length(fds);
            end
            
            varheader = {'loc'; 'totalReadDNA'; 'totalReadRNA'; ...
                'altDnaAF'; 'altRnaAF'; 'DEpvalAdj'};
            
            addcol = 0;
            for i = 1:length(varheader)
                addcol = addcol + size(obj.(varheader{i}),2);
            end
            
            if islogical(varidx)
                table = cell(sum(varidx), length(fdidx)+addcol+length(obj.variantAttr_cell_name)+length(obj.variantAttr_mtx_name));
                header = cell(1, size(table,2));
            else
                table = cell(length(varidx), length(fdidx)+addcol+length(obj.variantAttr_cell_name)+length(obj.variantAttr_mtx_name));
                header = cell(1, size(table,2));
            end
            coli = 1;
            for i = 1:length(varheader)
                addcol = size(obj.(varheader{i}), 2);
                table(:, coli:coli+addcol-1) = arrayfun(@num2str, obj.(varheader{i})(varidx,:,1), 'unif', 0);
                header(coli:coli+addcol-1) = varheader(i);
                coli = coli + addcol;
            end
            table(:, coli:coli+length(obj.variantAttr_mtx_name)-1) = ...
                arrayfun(@num2str, obj.variantAttr_mtx(varidx, :), 'unif', 0);
            header(coli:coli+length(obj.variantAttr_mtx_name)-1) = obj.variantAttr_mtx_name;
            coli = coli+length(obj.variantAttr_mtx_name);
            
            table(:, coli:coli+length(obj.variantAttr_cell_name)-1) = obj.variantAttr_cell(varidx,:);
            header(coli:coli+length(obj.variantAttr_cell_name)-1) = obj.variantAttr_cell_name;
            coli = coli+length(obj.variantAttr_cell_name);
            
            for i = 1:length(fdidx)                
                table(:, coli) = obj.annt.(fds{fdidx(i)})(varidx);
                header{coli} = fds{fdidx(i)};
                coli = coli + 1;
            end
            table(cellfun(@isempty, table)) = {'NaN'};
            idx = cellfun(@isscalar, table); 
            table(idx) = cellfun(@mat2str, table(idx), 'unif', 0);
        end
        
        function readAnnotationText(obj, filename)
            
            t = parseText(filename, 'ncolname', 1, 'nrowname', 0);
            chrm = t.text(:, strcmp(t.colname, 'CHROM'));
            chrm = strrep(chrm, 'X', '23');
            chrm = strrep(chrm, 'Y', '24');
            chrm = strrep(strrep(chrm, 'MT', '25'), 'M', '25');
            lockey = strcat(chrm, '-', t.text(:, strcmp(t.colname, 'POS')));
            plockey = strcat(arrayfun(@num2str, obj.loc(:,1), 'unif', 0), '-', ...
                arrayfun(@num2str, obj.loc(:,2), 'unif', 0) );
            [~, idx] = ismember(plockey, lockey);
            keep = false(length(t.colname), 1);
            tmp = lower(t.colname);
            keywords = lower({'Cosmic', 'Onco', 'SNPEFF', 'GO', 'KEGG', 'MiTarget', 'TF_Target'});            
            for i  = 1:length(keywords)
                keep = keep | ...
                    ~cellfun(@isempty, strfind(tmp, keywords{i})) ;
            end            
            keep = find(keep);
            t.colname = strrep(t.colname, '.', '_');
            
            nvar = size(obj.loc, 1);
            for i = 1:length(keep)
                obj.annt.(t.colname{keep(i)}) = cell(nvar, 1);                
                obj.annt.(t.colname{keep(i)})(idx~=0) = t.text(idx(idx~=0), keep(i));
                obj.annt.(t.colname{keep(i)})(idx==0) = {'NA'};
            end
        end        
        
        function addStat(obj, DESeqData, CENTROMERE)
            [nvar, nsamp] = size(obj.refCountDNA);
            obj.totalReadDNA = zeros(nvar, nsamp);
            obj.totalReadRNA = zeros(nvar, nsamp);
            for i = 1:nsamp
                obj.totalReadDNA(:,i) = full(sum(obj.ALLCountDNA{i}(:, 2:7),2)); %exclude those in introns
                obj.totalReadRNA(:,i) = full(sum(obj.ALLCountRNA{i}(:, 2:7),2)); %exclude those in introns
            end
            % calculate DNA, RNA allele-frequency
            obj.refDnaAF = obj.refCountDNA ./ obj.totalReadDNA;
            obj.altDnaAF = bsxfun(@rdivide, obj.altCountDNA, obj.totalReadDNA);
            obj.refRnaAF = obj.refCountRNA ./ obj.totalReadRNA;
            obj.altRnaAF = bsxfun(@rdivide, obj.altCountRNA, obj.totalReadRNA);            
            obj.RAratioDNA = obj.refCountDNA ./ (obj.altCountDNA(:,:,1) + 1); %avoid divided by 0
            obj.RAratioRNA = obj.refCountRNA ./ (obj.altCountRNA(:,:,1) + 1); %avoid divided by 0
            
            if ~isempty(CENTROMERE)
                %calculate distance to centromere
                obj.distToCen = NaN(nvar,1);
                for i = 1:24
                    varidx = obj.loc(:,1)==i;
                    cloc = CENTROMERE.loc(CENTROMERE.loc(:,1)==i,2:3);
                    cloc = cloc(:)';
                    obj.distToCen(varidx) = min(abs(bsxfun(@minus, cloc, obj.loc(varidx,2))),[],2);
                    clear varidx cloc
                end
            end
            
            if ~isempty(DESeqData)
                %get DESeq data
                [~, deidx] = GENOMEFUNC.isWithinRegion(obj.loc, DESeqData.loc);
                obj.DEpvalAdj = NaN(nvar,1);
                obj.DEpvalAdj(deidx~=0) = DESeqData.padj(deidx(deidx~=0));

                obj.DElabel = cell(nvar,1);
                obj.DElabel(:) = {'NA'};
                obj.DElabel(obj.DEpvalAdj<0.01) = {'padj < 0.01'};
                obj.DElabel(obj.DEpvalAdj>=0.5) = {'padj >= 0.5'};
                obj.DElabel(obj.DEpvalAdj>=0.01 & obj.DEpvalAdj<0.5) = {'in-between'};
            end
            
            %currently selected variants
            obj.subsetidx = true(nvar,1);
            
            % predefined sets of variants
            obj.presetidx.lowDNACov = all(obj.totalReadDNA < obj.minReadDNA, 2);
            obj.presetidx.highDNACov = any(obj.totalReadDNA > obj.maxReadDNA, 2);
            obj.presetidx.expressedEither = any(obj.totalReadRNA >= obj.minReadRNA, 2); %expressed in either tumor or normal
            obj.presetidx.expressedBoth = all(obj.totalReadRNA >= obj.minReadRNA, 2); %expressed in both tumor and normal
            
            obj.presetidx.multAlt = ~cellfun(@isempty, obj.alt(:,2));
            obj.presetidx.gatkindel = strcmp(obj.variantAttr_cell(:, ...
                strcmp(obj.variantAttr_cell_name, 'set')), 'gatkIndel');
            obj.presetidx.chrMT = obj.loc(:,1) == 25;
            obj.presetidx.DE01 = obj.DEpvalAdj < 0.01;
        end
        
        function readData(obj, patientSampleInfo, varargin)
            %patientSampleInfo = {'id', 'sample', 'sample'};
            %eg {'1381381', '4N', '2T'}
            para.r = 1;
            para.dRNA = 20;
            para.dDNA = 200;
            para.q = 0;
            para.aCountDir = 'data/';
            
            para = assignpara(para, varargin{:});
            chrmNum = {'X', '23'; 'Y', '24'; 'MT', '25'; 'M', '25'};
            
            
            %load vcf, keep chrm only
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
            
            
            %loop over samples to obtain allele read counts
            nsample = size(patientSampleInfo,1) - 1;
            aDataArray = cell(nsample, 2); %RNA, DNA
            lockey = cell(nsample,2);
            locidx = cell(nsample,2);
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
            
            obj.patient = aDataArray{1}.patient;
            obj.sample = patientSampleInfo(2:end);
            obj.minRead = aDataArray{1}.minRead;
            obj.maxDepth = aDataArray{1}.maxDepth;
            obj.minMapQ = aDataArray{1}.minMapQ;
            obj.loc = zeros(nvariant, 2);            
            obj.ref = cell(nvariant, 1); %take from allele-count data
            obj.alt_original = cell(nvariant, 1); %string
            
            altColIdx = strcmp(vcfData.attrName_cell, 'ALT');
            %get multiple alt
            multAlts = regexp( ...
                vcfData.variantAttr_cell(toVcfIdx, altColIdx), ...
                ',([ATCGNatcgn])', 'tokens');
            nalt = cellfun(@length, multAlts)+1;
            obj.alt = cell(nvariant, max(nalt) );
            obj.alt(:,1) = ... %take from vcf; but ALT can have multiple alleles
                upper( cellstr(cellfun(@(x) x(1), ...
                vcfData.variantAttr_cell(toVcfIdx, altColIdx))) );
            rowidx = find(nalt>1);
            for ri = 1:length(rowidx)
                obj.alt(rowidx(ri), 2:nalt(rowidx(ri))) = horzcat(multAlts{rowidx(ri)}{:});
            end
            for ci = 2:size(obj.alt,2)
                obj.alt(nalt==ci-1, ci:end) = {''};
            end
            obj.alt = upper(obj.alt);
            obj.alt_original = vcfData.variantAttr_cell(toVcfIdx, altColIdx);
            
            obj.refCountRNA = zeros(nvariant, nsample);
            obj.altCountRNA = zeros(nvariant, nsample, max(nalt));
            obj.refCountDNA = zeros(nvariant, nsample);
            obj.altCountDNA = zeros(nvariant, nsample, max(nalt));
            obj.ntlabels = aDataArray{sampidx}.ntlabels(1,:);
            obj.ALLCountDNA = cell(nsample, 1); %because want to store sparse matrix
            obj.ALLCountRNA = cell(nsample, 2); %because want to store sparse matrix
            obj.indelsRNA = cell(nvariant, nsample);
            obj.indelsDNA = cell(nvariant, nsample);
            
            nnt = length(obj.ntlabels);
            for sampidx = 1:nsample
                for dataidx = 1:length(datatype)
                    [~, locidx{sampidx,dataidx}] = ismember(lockey{sampidx,dataidx}, uloc);
                    
                    %loc and ref will overwrite between samples; but they should be the
                    %same between samples anway
                    if dataidx == 1
                        obj.loc(locidx{sampidx,dataidx},:) = aDataArray{sampidx,dataidx}.loc;
                        obj.ref(locidx{sampidx,dataidx}) = aDataArray{sampidx,dataidx}.ref;
                    end                    
                    obj.(['indels' datatype{dataidx}])(locidx{sampidx,dataidx}, sampidx) ...
                        = aDataArray{sampidx,dataidx}.indels;
                    
                    nvarSamp = length(aDataArray{sampidx,dataidx}.ref);
                    [~, colidx] = ismember(aDataArray{sampidx,dataidx}.ref, obj.ntlabels);
                    obj.(['refCount' datatype{dataidx}])(locidx{sampidx,dataidx}, sampidx) = ...
                        aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], (1:nvarSamp)', colidx)  );
                    
                    for altidx = 1:size(obj.alt,2)
                        [~, colidx] = ismember(obj.alt(locidx{sampidx,dataidx},altidx), obj.ntlabels);
                        obj.(['altCount' datatype{dataidx}])(locidx{sampidx,dataidx}(colidx~=0), sampidx, altidx) = ...
                            aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt], find(colidx~=0), colidx(colidx~=0)) );
                    end
                    
                    [rowidx, colidx] = find(aDataArray{sampidx,dataidx}.acount);
                    obj.(['ALLCount' datatype{dataidx}]){sampidx} = sparse( locidx{sampidx,dataidx}(rowidx), colidx, ...
                        aDataArray{sampidx,dataidx}.acount( sub2ind([nvarSamp, nnt],rowidx,colidx) ), ...
                        nvariant ,nnt);
                end
            end
            
            %what to copy from vcf?
            [~, vcfFormatSampleOrder] = ismember(patientSampleInfo(2:end), vcfData.sample);
            keep = ~ismember(vcfData.attrName_cell, {'CHROM', 'REF', 'ALT', 'FILTER'});
            obj.variantAttr_cell_name = vcfData.attrName_cell(keep);
            obj.variantAttr_cell = vcfData.variantAttr_cell( toVcfIdx, keep);
            keep = ~ismember(vcfData.attrName_mtx, {'POS'});
            obj.variantAttr_mtx_name = vcfData.attrName_mtx(keep);
            obj.variantAttr_mtx = vcfData.variantAttr_mtx( toVcfIdx, keep);
            fds = fieldnames(vcfData.formatData);
            for fdidx = 1:length(fds)
                obj.formatData.(fds{fdidx}) = vcfData.formatData.(fds{fdidx})( toVcfIdx, vcfFormatSampleOrder, :);
            end
            obj.vcfSample = vcfData.sample(vcfFormatSampleOrder);
            obj.vcfDesc = vcfData.desc;
            
            %sort alt by count
            rowidx = nalt > 1;
            [~, ntreorder] = sort(squeeze(sum(obj.altCountRNA(rowidx,:,:), 2)), 2, 'descend');
            rowidx = find(rowidx);
            fds = fds(ismember(fds, vcfData.allelicFormatField));
            [~, fdsDescIdx] = ismember(fds, vcfData.desc(:,1));
            altOnly = strcmp(vcfData.desc(fdsDescIdx, 2), 'A');
            altDataIdx = 2:max(nalt)+1;
            infoIdx = ismember(vcfData.desc(:,1), vcfData.attrName_cell) & ...
                strcmp(vcfData.desc(:,2), 'A') & ~strcmpi(vcfData.desc(:,1), '1000genomes.AF');
            [~, infoIdx] = ismember(vcfData.desc(infoIdx,1), obj.variantAttr_cell_name);
            infoIdx(infoIdx==0) = [];
            for ri = 1:length(rowidx)
                obj.altCountRNA(rowidx(ri),:,:) = ...
                    obj.altCountRNA(rowidx(ri),:,ntreorder(ri,:));
                obj.altCountDNA(rowidx(ri),:,:) = ...
                    obj.altCountDNA(rowidx(ri),:,ntreorder(ri,:));
                obj.alt(rowidx(ri),:) = obj.alt(rowidx(ri), ntreorder(ri,:));
                %for any field related to alt alleles, reorder the data
                for fdidx = 1:length(fds)
                    if altOnly(fdidx)
                        obj.formatData.(fds{fdidx})(rowidx(ri),:,:) = ...
                            obj.formatData.(fds{fdidx})(rowidx(ri), :, ntreorder(ri,:));
                    else
                        obj.formatData.(fds{fdidx})(rowidx(ri),:,altDataIdx) = ...
                            obj.formatData.(fds{fdidx})(rowidx(ri), :, altDataIdx(ntreorder(ri,:)));
                    end
                end
                for infoii = 1:length(infoIdx)
                    obj.variantAttr_cell{rowidx(ri), infoIdx(infoii)} = ...
                        obj.variantAttr_cell{rowidx(ri), infoIdx(infoii)}( ntreorder( ri, 1:nalt(rowidx(ri)) ) );
                end
            end
            
            obj.callerColidx = find(strcmp(obj.variantAttr_cell_name, 'set'));
            
            %if snpeff is available
            if ~isempty(vcfData.snpEff)
                snpEff = vcfData.snpEff( toVcfIdx );
                neff = cellfun(@(x) size(x,1), vcfData.snpEff);
                flatVarIdx = arrayfun(@(x) ones(neff(x),1)*x, (1:length(snpEff))');
                obj.annt.SnpEffName = ['varidx'; vcfData.snpEffFormat(:)];
                obj.annt.SnpEffMtx = [vertcat(flatVarIdx{:}), strkey( vertcat(snpEff{:}), 'add' )];
                obj.annt.SnpEff_metadb = strkey('dbname');
                %remove entries without annotations
                rmi = all(obj.annt.SnpEffMtx(:, 2:end)==0,2);
                obj.annt.SnpEffMtx(rmi,:) = [];
            end
        end
    end
    
    methods (Static=true)
        function fig = setupFigure(fig, position)
            if nargin < 2, fig = []; end
            if nargin < 3, position = [0, 0, 1000, 650]; end
            if ~isempty(fig)
                clf(fig);                
            else
                fig = figure('position', position);
                set(fig, 'paperpositionmode', 'auto');
            end
        end
        
        function para = figDefaultPara(type)            
            if strcmpi(type, 'scatter')
                para.func = @patchScatter; %only valid when colorBy is empty
                para.size = 40;
                para.colormap = jet(50);
                para.climPrctile = [5, 95]; %percentile for colormap if using scatter
                para.clim = []; %hard setting for colormap range; overwrite climPrctile
                para.fig = '';
                para.savedir = 'figures/';
                para.save = '';
                para.colorBy = '';
                para.colorByTitle = '';
                para.saveAppendColorBy = true;
                para.valueOrder = ''; %{'ascend', 'descend', {'labels'} (for group-scatter)
            else
                para = [];
            end
        end                
        
        function saveFigure(varargin)
            % saveFigure(handle, filename, format)
            % saveFigure(figpara)
            % saveFigure('para', value)
            %            
            para.fig = '';
            para.savedir = 'figures/';
            para.save = '';
            para.colorBy = '';
            para.colorByTitle = '';
            para.saveAppendColorBy = true;
            para.format = 'png';
            para.sampleid = '';
            
            filename = '';
            if ishandle(varargin{1})
                para.fig = varargin{1};
                filename = varargin{2};
                if length(varargin) > 2
                    para.format = varargin{3};
                end
            elseif length(varargin)==1 && isstruct(varargin{1})
                para = varargin{1};
            else
                para = assignpara(para, varargin{:});
            end
            
            if isempty(filename)
                filename = [para.savedir '/' para.save];
                if para.saveAppendColorBy
                    filename = [filename ', colored by ' para.colorByTitle];
                end
                if ~isempty(para.sampleid)
                    filename = [filename ', ' para.sampleid];
                end
            end
            
            if strcmpi(para.format, 'svg')
                plot2svg([filename '.svg'], fig);
            elseif strcmpi(para.format, 'eps')
                saveas(fig, [filename '.eps'], 'psc2');
            else
                saveas(fig, [filename, '.' para.format], para.format);
            end
            
        end
    end
end
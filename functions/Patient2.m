classdef Patient2 < handle
    properties
        patient;
        variant;
        
        cdict;
        
    end
        
    methods
        function obj = Patient2(patientSampleInfo, DESeqData, CENTROMERE, varargin)            
            obj.variant = Variant(patientSampleInfo, DESeqData, CENTROMERE, varargin{:});
            
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
            
            obj.variant = obj.variant.selectVariant(selection);            
        end
                
        function availSet = availableVarPreset(obj)
            availSet = obj.variant.availablePreset();
        end
        
        function [idx, desc] = currentVarSet(obj)
            [idx, desc] = obj.variant.currentSet();
        end
                                
        function fig = plotRAratio(obj, varargin)
            para.fig = '';
            para.save = 'REF-ALT ratio, normal vs tumor';
                        
            para = assignpara(para, varargin{:});
            fig = obj.setupFigure(para.fig);
            
            datatype = {'DNA', 'RNA'};
            for dataidx = 1:length(datatype)
                subplot(2,2, (dataidx-1)*2+1);
                h = gscatter(obj.variant.(['RAratio' datatype{dataidx}])(obj.currentVarSet, 1), ...
                    obj.variant.(['RAratio' datatype{dataidx}])(obj.currentVarSet, 2), ...
                    obj.variant.sets(obj.currentVarSet));
                set(h, 'marker', 'o', 'markersize', 10, 'linewidth',1);
                for hi = 1:length(h)
                    set(h(hi), 'color', obj.cdict.color(get(h(hi), 'displayname')));
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
                h = gscatter(obj.variant.(['RAratio' datatype{dataidx}])(obj.currentVarSet, 1), ...
                    obj.variant.(['RAratio' datatype{dataidx}])(obj.currentVarSet, 2), ...
                    obj.variant.sets(obj.currentVarSet));
                for hi = 1:length(h)
                    set(h(hi), 'color', obj.cdict.color(get(h(hi), 'displayname')));
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
                [hscatter, plottype] = obj.plotScatter(obj.variant.(xfd{plotidx})(obj.currentVarSet,xsample(plotidx),1), ...
                    obj.variant.(yfd{plotidx})(obj.currentVarSet,ysample(plotidx),1), ...
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
                    [hscatter, plottype] = obj.plotScatter(obj.variant.(xfd)(obj.currentVarSet,sampi,1), ...
                        obj.variant.(yfd)(obj.currentVarSet,sampi,1), ...
                        'passpara', parastr{:});

                    if aidx == 2 && sampi == 1 && strcmpi(plottype, 'gpatchscatter')
                        h = legend('orientation', 'horizontal');
                        set(h, 'position', [0.37, 0, 0.3, 0.04], 'fontsize', 14);
                    end
                    xlabel('DNA', 'fontsize', 14);
                    ylabel('RNA', 'fontsize', 14);
                    title(sprintf('%s freq, %s', upper(allele{aidx}), obj.variant.sample{sampi}), 'fontsize', 14);
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
            
            GT = strcat(obj.variant.GT(:,1), '-', obj.variant.GT(:,2));
            for i = 1:min(6, length(para.plotGT))
                idx = obj.currentVarSet & strcmp(GT, para.plotGT{i});
                subplot(2, 3, i);
                gscatter(obj.variant.altDnaAF(idx, 1, 1), obj.variant.altDnaAF(idx, 2, 1), ...
                    obj.variant.sets(idx));
                xlabel(sprintf('ALT AF in %s', obj.variant.sample{1}), 'fontsize', 14);
                ylabel(sprintf('ALT AF in %s', obj.variant.sample{2}), 'fontsize', 14);                
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
                    [hscatter, plottype] = obj.plotScatter(log10(obj.variant.(fd)(obj.currentVarSet,1,1)+1), ...
                        log10(obj.variant.(fd)(obj.currentVarSet,2,1)+1), ...
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
                        
            callers = unique(obj.variant.sets);
            for i = 1:length(callers)
                idx = obj.currentVarSet & strcmp(obj.variant.sets, callers{i});
                subplot(3,3,i);
                if sum(idx) == 1
                    tmp = [obj.variant.altDnaAF(idx,:,1), obj.variant.altRnaAF(idx,:,1)];
                    for si = 1:4
                        hold on;
                        plot([tmp(si), tmp(si)], [0, 1], 'color', para.barcolors(si,:), 'linewidth', 4);
                        hold off;
                    end
                    ylim([0 2]);
                    clear tmp
                else
                    [h, x] = hist([obj.variant.altDnaAF(idx,:,1), obj.variant.altRnaAF(idx,:,1)], 0:0.2:1);
                    bar(x, h, 'hist');
                    ylim([0, max(h(:))+1]);
                    xlim([-0.1 1.1]);
                    set(gca, 'box', 'off');
                end
                xlabel('ALT allele freq', 'fontsize', 12);
                ylabel('number of variants', 'fontsize', 12);
                title(callers{i}, 'fontsize', 14);
                if i == length(callers)
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
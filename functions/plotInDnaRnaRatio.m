function [ax, linehandle] = plotInDnaRnaRatio(data1, data2, xticklabel, varargin)
    para.color = [0, 0, 1; 0, 0.5, 0; 1, 0 0; 0.75, 0, 0.75];
    para.linewidth = 2;
    para.linestyle = {'-', '--'};
    para.ylabel = {'% in DNA calls', '% in RNA calls'};
    para.fontsize = 14;
    para.linename = {};
    
    para = assignpara(para, varargin{:});
    
    if ~iscell(para.ylabel) para.ylabel = {para.ylabel}; end
        
    nx = size(data1, 1);
    if nx ~= length(xticklabel)
        error('number of xticklabel does not match number of rows in data');
    end
    if ~isempty(data2)
        linehandle = cell(1,2);
        [ax, linehandle{1}, linehandle{2}] = plotyy(1:nx, data1, 1:nx, data2);
        set(ax(2),'xtick',[]);    
    else
        linehandle{1} = plot(1:nx, data1);
        ax = gca;
    end
    
    for i = 1:length(linehandle)
        set(linehandle{i}, 'linestyle', para.linestyle{i}, ...
            'linewidth', para.linewidth, 'marker', 'x');        
        set(get(ax(i),'ylabel'), 'string', para.ylabel{i}, 'fontsize', para.fontsize);
        set(ax(i), 'xlim', [0.5, nx+0.5]);   
        if i == 1
            M = max(data1(:)); m = min(data1(:));
        else
            M = max(data2(:)); m = min(data2(:));
        end
        margin = (M-m) / 10;
        set(ax(i), 'ylim', [m-margin, M+margin]);
        for j = 1:size(data1,2)
            set(linehandle{i}(j), 'color', para.color(j,:));
        end
    end        
    set(ax(1), 'xtick', 1:nx, 'xticklabel', xticklabel, 'fontsize', para.fontsize);    
    if ~isempty(para.linename)
        legend(para.linename, 'fontsize', para.fontsize, 'location', 'best');        
    end
end
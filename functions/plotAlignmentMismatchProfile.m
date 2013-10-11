function [t, nread] = plotAlignmentMismatchProfile( fns, outdir, doplot )
    if nargin < 3, doplot = true; end
    if ischar(fns)
        fns = {fns};
    end
    
    nread = 0;
    for fnidx = 1:length(fns)
        f = fopen(fns{fnidx}, 'r');
        nread = nread + fscanf(f, '%d', 1);
        fclose(f);
                
        %data format: first/second-fwd/rev-(r)[MIDS]
        % "r" is for read count of the number of MIDS; it's basically the
        %distribution of the number of MIDS for reads
        % Without "r", it's simply positional read counts of MIDS
        subt = parseText(fns{fnidx}, 'skip', 1, 'nrowname', 1, 'ncolname', 0, 'numeric', true);
        if fnidx == 1
           t = subt; 
           [~, si] = sort(t.rowname);
           t.rowname = t.rowname(si);
           t.text = t.text(si,:);
        else
            [~, si] = ismember(t.rowname, subt.rowname);
            t.text = t.text + subt.text(si,:);
        end        
    end
    if ~doplot
        return
    end
    n = size(t.text,2);
    %first column data: in positional, it's position 0, so all counts are
    %zeros; in read-distribution, it's 0 mismatches (or whichever the
    %category is), so most of the read counts should be here
    
    posIdx = ~cellfun(@isempty, regexp(t.rowname, '\w+\-\w+\-[DIMS]'));    
    forward = ~cellfun(@isempty, strfind(t.rowname, 'fwd'));
    first = ~cellfun(@isempty, strfind(t.rowname, 'first'));
    
    readlength = size(t.text,2) - 1;
    t.text = t.text ./ (nread/2) * 100;  %percentage
    
    plotcat = {'M', 'I', 'D', 'S'};
    catlabels = {'mismatch', 'insertion', 'deletion', 'soft-clip'};
    catidx = false(length(t.rowname), length(plotcat));
    for i = 1:length(plotcat)
        catidx(:,i) = ~cellfun(@isempty, strfind(t.rowname, plotcat{i}));
    end
    
    
    set(gcf, 'position', [0,0, 850, 900], 'paperpositionmode', 'auto', 'visible', 'off');    
    for i = 1:2 % column1: positional; column2: read-distribution        
        for j = 1:length(plotcat)
            subplot(4,2, (j-1)*2+i);            
            if i == 1
                idx = posIdx & catidx(:,j);                 
                graphdata = [t.text(idx & first & forward, 2:end) + t.text(idx & first & ~forward,n:-1:2); 
                    t.text(idx & ~first & forward, 2:end) + t.text(idx & ~first & ~forward,n:-1:2) ]'; 
            else
                idx = ~posIdx & catidx(:,j);
                graphdata = [sum(t.text(idx & first,2:end), 1); sum(t.text(idx & ~first,2:end), 1)]';
            end
            
            plot(graphdata, '.-', 'markersize', 20, 'linewidth', 1.5);
            m = min(graphdata(:));
            M = max(graphdata(:));
            mtick = floor(m/5)*5;
            ytickstep = (M-mtick) / 4;
            if ytickstep < 0.5
                pw = ceil(abs(log10(ytickstep)));
                ytickstep = round(10^pw*ytickstep)/10^pw;
            elseif ytickstep < 1
                ytickstep = 1;
            elseif ytickstep < 5
                ytickstep = round(ytickstep);
            elseif ytickstep > 5
                ytickstep = round(ytickstep/5)*5;
            end
            if readlength < 51
                xtickstep = 10;
            elseif readlength < 100
                xtickstep = 20;
            else
                xtickstep = 50;
            end
            set(gca, 'xtick', 10:xtickstep:readlength, 'ytick', mtick:ytickstep:M, 'fontsize', 12);
            xlim([0.9, readlength+0.1]);
            ylim([mtick-ytickstep/10, M+ytickstep/10]);
%             fprintf('%g, %g, %g, %g\n',m, mtick, M, tickstep);
            if i == 1
                xlabel('position', 'fontsize', 12, 'fontweight', 'bold');
                ylabel('% of reads', 'fontsize', 12, 'fontweight', 'bold');
            else
                xlabel(['number of ' catlabels{j}], 'fontsize', 12, 'fontweight', 'bold');
                ylabel('% of reads', 'fontsize', 12, 'fontweight', 'bold');
            end
        end
    end

    lastdir = strfind(fns{1}, '/');
    if ~isempty(lastdir)
        basefn = fns{1}(lastdir(end)+1:end);
    else
        basefn = fns{1};
    end
    
    h = suptitle(basefn);
    set(h, 'interpreter', 'none');
    
    h = legend({'read1', 'read2'}, 'fontsize', 12, 'orientation', 'horizontal');
    hpos = get(h, 'position');
    set(h, 'position', [0.4, 0.9, hpos(3), hpos(4)]);
    
    if outdir(end) ~= '/'
        outdir = [outdir '/'];
    end    
    plot2svg([outdir basefn '.mismatchStat.svg'], gcf);    
    
    numreads = t.text .* (nread/2) / 100;
    set(gcf, 'position', [0,0, 900, 700], 'paperpositionmode', 'auto', 'visible', 'off');    
    for typeidx = 1:length(plotcat)
        [~, i] = ismember(strcat({'first-fwd-', 'first-rev-', 'second-fwd-','second-rev-'}, plotcat{typeidx}), t.rowname);
        clf
        for j = 1:4
            subplot(2,2,j);
            plot(numreads(i(j),2:end),'.-', 'linewidth', 1.2);            
            title(t.rowname{i(j)}, 'fontsize', 12)
            xlabel('position', 'fontsize', 12);
            ylabel('# reads', 'fontsize', 12);                        
            
%             ax1 = gca;
%             ytick = get(ax1, 'ytick');            
%             yticklabel = ytick ./ (nread/2) * 100;
%             pw = ceil(abs(log10(yticklabel(1))));
%             yticklabel = round(10^pw.*yticklabel)/10^pw;
%             ax2 = axes('position', get(ax1, 'position'), ...
%                 'xaxislocation', 'bottom', 'yaxislocation', 'right', ...
%                 'color', 'none', 'xcolor', 'k', 'ycolor', 'k');
%             ylim(ax2, ylim(ax1));
%             set(ax2, 'ytick', ytick, 'yticklabel', yticklabel, 'xtick', []);
%             ylabel(ax2, '% of reads', 'fontsize', 12);
        end
        saveas(gcf, [outdir plotcat{typeidx} '_pos.' basefn '.png' ], 'png');
    end
    close(gcf);
end
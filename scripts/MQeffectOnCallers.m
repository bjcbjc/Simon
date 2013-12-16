

patientpair  = {'132540-10N_132540-1T', '132540-10N_137064-1T', '138381-4N_138381-2T'};
quallabel = {'', '.AS'};
caller = {'virmid'}; %{'varscan', 'sniper', 'strelka'};
callertable = {'varscan', '%s.mdup%s.varscan.snv.vcf', [1,2,8], {}, '%s %f %s'; ...
    'sniper', '%s.mdup%s.sniper.vcf', [1, 2], {'SS'}, '%s %f %f %f'; ...
    'strelka', '%s.mdup%s.strelka.snv.vcf', [1,2,7], {}, '%s %f %s'; ...
    'virmid', '%s.mdup%s.som.virmid.vcf', [1,2,7], {}, '%s %f %s'};
vcfdata = cell(1,length(quallabel));
for cidx = 1:length(caller)
    tic;
    for pidx = 1:length(patientpair)
        for qidx = 1:length(quallabel)
            tbidx = find(strcmp(callertable(:,1), caller{cidx}));
            if ~strcmp(caller{cidx}, 'virmid')
                [vcfdata{qidx}, status] = VCFFUNC.extract( ...
                    [sprintf('RNA/%s_star/', caller{cidx}), ...
                    sprintf(callertable{tbidx,2}, patientpair{pidx}, quallabel{qidx})], ...
                    'returncol', callertable{tbidx,3}, 'format', callertable{tbidx,4}, 'parseformatstr', callertable{tbidx,5});
            else
                [~, tumorname] = strtok(patientpair{pidx}, '_');
                [vcfdata{qidx}, status] = VCFFUNC.extract( ...
                    [sprintf('RNA/%s_star/', caller{cidx}), ...
                    sprintf(callertable{tbidx,2}, tumorname(2:end), quallabel{qidx})], ...
                    'returncol', callertable{tbidx,3}, 'format', callertable{tbidx,4}, 'parseformatstr', callertable{tbidx,5});
            end
            if isempty(tbidx)
                fprintf('unknown caller\n');
                continue
            end
            if status ~= 0
                fprintf('error extracting vcf\n');
                fprintf('%s\n', vcfdata{qidx});
                continue
            end            
            vcfdata{qidx}{1} = gloc2index( numericchrm(vcfdata{qidx}{1}), vcfdata{qidx}{2} );
            if strcmp(caller{cidx}, 'varscan')
                vcfdata{qidx}{3} = ~cellfun(@isempty, strfind(vcfdata{qidx}{3}, 'SS=2;'));
            elseif strcmp(caller{cidx}, 'sniper')
                vcfdata{qidx}(3) = [];
                vcfdata{qidx}{3} = vcfdata{qidx}{3} == 2;
            elseif strcmp(caller{cidx}, 'strelka') || strcmp(caller{cidx}, 'virmid')
                vcfdata{qidx}{3} = strcmp(vcfdata{qidx}{3}, 'PASS');                            
            end
        end
        subplot(2,3,pidx)
        fixvenn({vcfdata{1}{1}, vcfdata{2}{1}}, 'label', {'no MQ', 'STAR Q'}, 'labelfontsize',12);
        title(strrep(patientpair{pidx}, '_','-'), 'fontsize',12);
        
        subplot(2,3,pidx+3)
        fixvenn({vcfdata{1}{1}(vcfdata{1}{3}), vcfdata{2}{1}(vcfdata{2}{3})}, 'label', {'no MQ', 'STAR Q'}, 'labelfontsize',12);
        title([strrep(patientpair{pidx}, '_','-'), ' , somatic'], 'fontsize',12);        
    end
    suptitle(caller{cidx});
    set(gcf, 'paperpositionmode', 'auto');
    saveas(gcf, sprintf('figures/rnavar/MQeffect/%s.png', caller{cidx}), 'png');
    toc;
end

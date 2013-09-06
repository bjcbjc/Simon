

cleanmutect = false;
cleanvarscan = true;

if cleanmutect
    % clean up mutect vcf mat
    %
    % FA basically is AD(alt)/(AD(alt)+AD(ref)), so remove
    %
    fns = listfilename('data/rnavar*mutect.mat');
    for i = 1:length(fns)
        changed = false;
        vcfData = loadStructData(['data/' fns{i}]);
        if vcfData.maxnalt ~= 1
            fprintf('skip %s, more than one ALT at sites\n', fns{i});
            continue;
        end
        if isfield(vcfData.formatData, 'FA')
            vcfData.formatData = rmfield(vcfData.formatData, 'FA');
            changed = true;
        end
        nvar = size(vcfData.variantAttr_cell,1);
        nsamp = length(vcfData.sample);
        if iscell(vcfData.formatData.AD)
            tmp = NaN(nvar, nsamp, 2);
            for sampi = 1:nsamp
                tmp(:,sampi, :) = str2double_fast( ...
                    vcfData.formatData.AD(:,sampi), '%g,%g#');
            end
            vcfData.formatData.AD = tmp;
            changed = true;
        end
        if iscell(vcfData.formatData.BQ)
            vcfData.formatData.BQ = str2double_fast( ...
                strrep(vcfData.formatData.BQ(:,2), '.', 'NaN') );
            changed = true;
        end
        if all(vcfData.formatData.GT(:)==0)
            varformatstr = '%*s %f %*s %*s %*s %*f %*s %*s'; %chrom, pos, id, ref, alt, qual, filter, info            
            [~, nheader] = system(sprintf('head -1000 %s | grep -c ^# ', vcfData.filename));
            nheader = str2double(nheader);
            f = fopen(vcfData.filename);
            l = fgetl(f);
            while l(1) == '#'
                l = fgetl(f);
            end
            fields = textscan(l, '%s');
            fields = fields{1};
            formatfields = textscan(fields{9}, '%s', 'delimiter', ':');
            formatfields = formatfields{1};
            gt = find(strcmp(formatfields, 'GT'));
            for gtii = 1:length(formatfields)
                varformatstr = strcat(varformatstr, '%*s');
            end
            for sampi = 1:nsamp
                for gtii = 1:length(formatfields)
                    if gtii == gt
                        varformatstr = strcat(varformatstr, '%s');
                    else
                        varformatstr = strcat(varformatstr, '%*s');
                    end
                end
            end        
            frewind(f);
            d = textscan(f, varformatstr, 'delimiter',sprintf('\t:'), 'headerlines', nheader, 'treatasempty', '.');
            fclose(f);        
            if any(vcfData.variantAttr_mtx(:,1)~=d{1})
                fprintf('variant inconst, %s\n', fns{i});
            else
                vcfData.formatData.GT = strkey([d{2}, d{3}], 'add');
            end
            changed = true;
        end
        if changed
            save(['data/' fns{i}], 'vcfData');
        end
        clear vcfData
    end
end


if cleanvarscan
    %FREQ is just AD/(AD+RD), remove
    %info DP is just sum of format DPs
    fns = listfilename('data/rnavar*varscan.mat');
    for i = 1:length(fns)
        changed = false;
        vcfData = loadStructData(['data/' fns{i}]);
        if isfield(vcfData.formatData, 'FREQ')
            vcfData.formatData = rmfield(vcfData.formatData, 'FREQ');
            changed = true;
        end
        if iscell(vcfData.formatData.GT)
            vcfData.formatData.GT = strkey(vcfData, 'add');
            changed = true;
        end
        if iscell(vcfData.formatData.DP4)
            [nvar, nsamp] = size(vcfData.formatData.DP4);
            tmp = NaN(nvar, nsamp, 4);
            for sampi = 1:nsamp
                tmp(:,smpi,:) = str2double_fast( ...
                    vcfData.formatData.DP4(:,sampi), '%g,%g,%g,%g#');
            end
            vcfData.formatData.DP4 = tmp;
            changed = true;
        end
        if all(isnan(vcfData.formatData.GQ(:)))
            vcfData.formatData = rmfield(vcfData.formatData, 'GQ');
            changed = true;
        end
        [~, rmi] = ismember('DP', vcfData.attrName_mtx);
        if rmi ~= 0
            vcfData.attrName_mtx(rmi) = [];
            vcfData.variantAttr_mtx(:, rmi) = [];
            changed = true;
        end
        if changed
            save(['data/' fns{i}], 'vcfData');
        end
    end
end
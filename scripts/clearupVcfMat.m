

cleanmutect = false;
cleanvarscan = false;
cleansniper = true;

if cleanmutect
    % clean up mutect vcf mat
    %
    % FA basically is AD(alt)/(AD(alt)+AD(ref)), so remove
    %
    fns = listfilename('data/rnavar*q50*mutect*.mat');
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
        if iscell(vcfData.formatData.GT)
            vcfData.formatData.GT = strkey(vcfData.formatData.GT, 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');            
            changed = true;
        elseif all(vcfData.formatData.GT(:)==0)
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
                vcfData.formatData.GT = strkey([d{2}, d{3}], 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
            end
            changed = true;
        end
        if changed
            fprintf('update %s\n',fns{i});
            save(['data/' fns{i}], 'vcfData');
        end
        clear vcfData
    end
end


if cleanvarscan
    %FREQ is just AD/(AD+RD), remove
    %info DP is just sum of format DPs
    fns = listfilename('data/rnavar*130*varscan.mat');
    for i = 1:length(fns)
        changed = false;
        vcfData = loadStructData(['data/' fns{i}]);
        if isfield(vcfData.formatData, 'FREQ')
            vcfData.formatData = rmfield(vcfData.formatData, 'FREQ');
            changed = true;
        end
        if iscell(vcfData.formatData.GT)
            vcfData.formatData.GT = strkey(vcfData.formatData.GT, 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
            changed = true;
        end
        if iscell(vcfData.formatData.DP4)
            [nvar, nsamp] = size(vcfData.formatData.DP4);
            tmp = NaN(nvar, nsamp, 4);
            for sampi = 1:nsamp
                tmp(:,sampi,:) = str2double_fast( ...
                    vcfData.formatData.DP4(:,sampi), '%g,%g,%g,%g#');
            end
            vcfData.formatData.DP4 = tmp;
            changed = true;
        end
        if isfield(vcfData.formatData, 'GQ')
            if all(isnan(vcfData.formatData.GQ(:)))
                vcfData.formatData = rmfield(vcfData.formatData, 'GQ');
                changed = true;
            end
        end
        [~, rmi] = ismember('DP', vcfData.attrName_mtx);
        if rmi ~= 0
            vcfData.attrName_mtx(rmi) = [];
            vcfData.variantAttr_mtx(:, rmi) = [];
            changed = true;
        end
        if vcfData.maxnalt == 1
            %SScode = [0, 1, 2, 3, 5];
            %SSname = {'Reference', 'Germline', 'Somatic', 'LOH', 'Unknown'};
            [~, ssidx] = ismember('SS', vcfData.attrName_cell);
            if ssidx ~= 0
                tmp = str2double_fast(vcfData.variantAttr_cell(:, ssidx));
%                 tmp(tmp<5) = tmp(tmp<5) + 1;
%                 tmp = strkey( SSname(tmp), 'add');
                vcfData.variantAttr_mtx(:, end+1) = tmp;
                vcfData.attrName_mtx{end+1} = 'SS';                                
                vcfData.variantAttr_cell(:, ssidx) = [];
                vcfData.attrName_cell(ssidx) = [];
                clear tmp
                changed = true;
            end
            [~, sscidx] = ismember('SSC', vcfData.attrName_cell);
            if sscidx ~= 0
                tmp = str2double_fast(vcfData.variantAttr_cell(:, sscidx));
                vcfData.variantAttr_mtx(:, end+1) = tmp;
                vcfData.attrName_mtx{end+1} = 'SSC';
                vcfData.variantAttr_cell(:, sscidx) = [];
                vcfData.attrName_cell(sscidx) = [];
                clear tmp
                changed = true;
            end
        end
        if changed
            fprintf('update %s\n',fns{i});
            save(['data/' fns{i}], 'vcfData');
        end
    end
end


if cleansniper
    %FREQ is just AD/(AD+RD), remove
    %info DP is just sum of format DPs
    fns = listfilename('data/rnavar*130*sniper*.mat');
    for i = 1:length(fns)
        changed = false;
        vcfData = loadStructData(['data/' fns{i}]);        
        
        if iscell(vcfData.formatData.GT)
            fillGT = false;
        else
            fillGT = all(vcfData.formatData.GT(:) == 0);
        end
        if iscell(vcfData.formatData.IGT)
            fillIGT = false;
        else
            fillIGT = all(vcfData.formatData.IGT(:) == 0);
        end
        if fillGT || fillIGT            
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
            igt = find(strcmp(formatfields, 'IGT'));
            readfield = {};
            for gtii = 1:length(formatfields)
                varformatstr = strcat(varformatstr, '%*s');
            end
            for sampi = 1:nsamp
                for gtii = 1:length(formatfields)
                    if (gtii == gt && fillGT) || (gtii == igt && fillIGT)
                        varformatstr = strcat(varformatstr, '%s');
                        if gtii == gt
                            readfield{end+1} = 'GT';
                        else
                            readfield{end+1} = 'IGT';
                        end
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
                if fillGT
                    idx = find(strcmp(readfield, 'GT')) + 1;
                    vcfData.formatData.GT = strkey(horzcat(d{idx}), 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
                end
                if fillIGT
                    idx = find(strcmp(readfield, 'IGT')) + 1;
                    vcfData.formatData.IGT = strkey(horzcat(d{idx}), 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
                end
            end
            changed = true;
        end
        
        if iscell(vcfData.formatData.GT)
            vcfData.formatData.GT = strkey(vcfData.formatData.GT, 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
            changed = true;        
        end
        if iscell(vcfData.formatData.IGT)
            vcfData.formatData.IGT = strkey(vcfData.formatData.IGT, 'add', 'lookup', '/nethome/bjchen/Projects/Simon/data/DICTIONARY.mat');
            changed = true;
        end
        if iscell(vcfData.formatData.DP4)
            [nvar, nsamp] = size(vcfData.formatData.DP4);
            tmp = NaN(nvar, nsamp, 4);
            for sampi = 1:nsamp
                tmp(:,sampi,:) = str2double_fast( ...
                    vcfData.formatData.DP4(:,sampi), '%g,%g,%g,%g#');
            end
            vcfData.formatData.DP4 = tmp;
            changed = true;
        end
        if iscell(vcfData.formatData.BCOUNT)
            [nvar, nsamp] = size(vcfData.formatData.BCOUNT);
            tmp = NaN(nvar, nsamp, 4);
            for sampi = 1:nsamp
                tmp(:,sampi,:) = str2double_fast( ...
                    vcfData.formatData.BCOUNT(:,sampi), '%g,%g,%g,%g#');
            end
            vcfData.formatData.BCOUNT = tmp;
            changed = true;
        end
        
        if changed
            fprintf('update %s\n',fns{i});
            save(['data/' fns{i}], 'vcfData');
        end
    end
end
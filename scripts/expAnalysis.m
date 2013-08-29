
%% use new patient class
patient = {'132540', '138381'; ...
    '10N', '4N'; ...
    '1T', '2T'};

para.r = 1;
para.dRNA = 20;
para.dDNA = 200;
para.q = 50;
para.aCountDir = 'data/';
paraArray = paraPairArray(para);
if ~exist('DESeqData', 'var')
    DESeqData = loadStructData('data/DESeq.NormalTumor.mat');
end
if ~exist('CENTROMERE', 'var')
    CENTROMERE = loadStructData('CENTROMERE.hg19.mat');
end
tic; P132540 = Patient2(patient(:,1), DESeqData, CENTROMERE, paraArray{:}); toc;
tic; P138381 = Patient2(patient(:,2), DESeqData, CENTROMERE, paraArray{:}); toc;

save data/patient.132540.obj2.q50.mat P132540
save data/patient.138381.obj2.q50.mat P138381
%%
pid = {'132540', '138381'};
anntfile = {'132540-10-N-132540-1-T.mutect-sniper-varscan-gatkIndel.union.snpEff.va.vcf.txt', ...
    '138381-4-N-138381-2-T.mutect-sniper-varscan-gatkIndel.union.snpEff.va.vcf.txt'};
anntdir = '/nethome/bjchen/Projects/Simon/DNA/vcf/annotated/';

for i = 1:length(pid)
    curpatient = eval(['P' pid{i}]);
    GT = strcat(curpatient.formatData.GT(:,1), '-', curpatient.formatData.GT(:,2));
    GTset = strcat(GT, '-', curpatient.variantAttr_cell(:, curpatient.callerColidx));
    
    score = curpatient.variantAttr_cell(:,strcmp(curpatient.variantAttr_cell_name, 'Varscan_SSC'));    
    score(cellfun(@isempty, score)) = {'NaN'};
    score = str2double(score);        
    score(:,2) = curpatient.formatData.Sniper_SSC(:,2);
    
    curpatient.presetidx.distToCen500 = curpatient.distToCen <= 500000;
    curpatient.presetidx.chr1_19_cen = (curpatient.loc(:,1) == 1 | curpatient.loc(:,1)==19) & curpatient.distToCen < 10^6;
    curpatient.presetidx.homAltNormal = strcmp(GT, '1/1-0/1') | strcmp(GT, '1/1-0/0');
    curpatient.presetidx.varscanHet2Het = strcmp(GTset, '0/1-0/1-varscan');
    curpatient.presetidx.sniper30 = score(:,2) <= 30;
    curpatient.presetidx.chrY = curpatient.loc(:,1) == 24;
    curpatient.readAnnotationText([anntdir anntfile{i}]);
end
%%

pid = {'132540', '138381'};

filters = { 'expressedBoth', '~multAlt', '~gatkindel', '~chrMT', ...
    '~highDNACov', '~lowDNACov', ...
    '~distToCen500', '~chr1_19_cen', ...
    '~homAltNormal', '~varscanHet2Het', '~sniper30'};
filtersDPonly = { 'expressedBoth', '~multAlt', '~gatkindel', '~chrMT'};
filters(2:end) = strcat('&', filters(2:end));
filtersDPonly(2:end) = strcat('&', filtersDPonly(2:end));

doplot = true;
if ~exist('fig', 'var') && doplot
    fig = P132540.setupFigure();
end
if exist('fig', 'var')
    if ~ishandle(fig)
        fig = P132540.setupFigure();
    end
end
%%

usefilter = filtersDPonly;
figdir = 'figures/filterRnaRead10/';

for i = 1:length(pid)
    curpatient = eval(['P' pid{i}]);
    caller = regexp(curpatient.variantAttr_cell(:,curpatient.callerColidx), '^(\w)+-*','tokens', 'once');
    caller = cellfun(@(x) x{:}, caller, 'unif', 0);
    DElabel = repmat({'non DE'}, size(curpatient.loc,1), 1);
    DElabel(curpatient.DEpvalAdj < 0.01) = {'DE'};
    
    fprintf('========== current patient: %s ============\n', pid{i});
    if strcmp(pid{i}, '138381')
        curpatient.selectVariant( strcat(usefilter{:}, '& ~chrY') );
    else
        curpatient.selectVariant( strcat(usefilter{:}) );
    end
    [~, desc] = curpatient.currentSet();
    fprintf('filter: %s\n', desc);

    if doplot        
        curpatient.plotScatterNormalTumorCount('fig',fig, 'savedir', figdir);        
        curpatient.plotScatterNormalTumorCount('fig', fig, ...
            'colorBy',DElabel(curpatient.currentSet), ...
            'colorByTitle', 'DESeq','valueOrder',{'non DE','DE'}, 'savedir', figdir);
        
        curpatient.plotScatterDnaRnaAF('fig',fig, 'savedir', figdir);
        curpatient.plotScatterDnaRnaAF('fig',fig,'colorBy',caller(curpatient.currentSet), ...
            'colorByTitle', 'caller', 'savedir', figdir);
        curpatient.plotScatterDnaRnaAF('fig',fig, ...
            'colorBy',log10(curpatient.totalReadRNA(curpatient.currentSet,2)+1), ...
            'colorByTitle', 'log10(read in tumor)', ...
            'colormap', genColorMap([1 0.2 0.2; 0.6,0.6,0.6; 0.2,0.2,1], 10), ...
            'climPrctile', [10 90], 'savedir', figdir);
        curpatient.plotScatterDnaRnaAF('fig', fig, ...
            'colorBy',DElabel(curpatient.currentSet), ...
            'colorByTitle', 'DESeq','valueOrder',{'non DE','DE'}, 'savedir', figdir);
        
        curpatient.plotScatterALTAF('fig',fig,'colorBy',caller(curpatient.currentSet), ...
            'colorByTitle', 'caller', 'savedir', figdir);
        curpatient.plotScatterALTAF('fig', fig, ...
            'colorBy',log10(curpatient.totalReadRNA(curpatient.currentSet,2)+1), ...
            'colorByTitle', 'log10(read in tumor)', ...
            'colormap', genColorMap([1 0.2 0.2; 0.6,0.6,0.6; 0.2,0.2,1], 10), ...
            'climPrctile', [10 90], 'savedir', figdir);
        curpatient.plotScatterALTAF('fig', fig, ...
            'colorBy',DElabel(curpatient.currentSet), ...
            'colorByTitle', 'DESeq','valueOrder',{'non DE','DE'}, 'savedir', figdir);
                
%         curpatient.plotRAratio('fig', fig);
%         curpatient.plotHistMajorALTAF('fig',fig);
    end

end
%%
if ~exist('ENSEMBL', 'var')
    load ENSEMBL
end
if ~exist('NCBI', 'var')
    load NCBI
end


altRnaAFDiff = curpatient.altRnaAF(:,2,1) - curpatient.altRnaAF(:,1,1);
zNormalTumor = zscore(altRnaAFDiff(curpatient.currentSet));
idx = find(curpatient.currentSet);
idx = idx(abs(zNormalTumor) > 2);
[tableNormalTumor, header] = curpatient.getCurSelVarAnnotation([], idx);

entrez = GENOMEFUNC.idQuery(ENSEMBL, tableNormalTumor(:,strcmpi(header,'SNPEFF_TRANSCRIPT_ID')), {'entrez'}, 'transcript');
tableNormalTumor(:,end+1) = {''};
header{end+1} = 'Desc';
tableNormalTumor(~isnan(entrez), end) = GENOMEFUNC.idQuery(NCBI, entrez(~isnan(entrez)), 'desc');
keepcol = ~all(ismember(tableNormalTumor, {'NA','NaN','NONE'}),1);
tableNormalTumor = [header(keepcol); tableNormalTumor(:, keepcol)];
tabwrite(sprintf('tables/P%s, preliminary, abs z more than 2 from diff in altRnaAF.txt',curpatient.patient), tableNormalTumor);


plot(curpatient.altRnaAF(curpatient.currentSet,1,1), curpatient.altRnaAF(curpatient.currentSet,2,1), '.', 'markersize',25,'color',[0.3,0.3,0.3])
hold on; plot(curpatient.altRnaAF(idx,1,1), curpatient.altRnaAF(idx,2,1), '.', 'markersize',25,'color',[0.8,0.2,0.2]); hold off
xlabel('normal', 'fontsize', 18)
ylabel('tumor', 'fontsize', 18)
title('ALT freq in RNA', 'fontsize', 18)


altDnaRnaDiff = curpatient.altRnaAF(:,2,1) - curpatient.altDnaAF(:,2,1);
zDnaRna = zscore(altDnaRnaDiff(curpatient.currentSet));
idx = find(curpatient.currentSet);
idx = idx(abs(zDnaRna) > 2);
[tableDnaRna, header] = curpatient.getCurSelVarAnnotation([], idx);

entrez = GENOMEFUNC.idQuery(ENSEMBL, tableDnaRna(:,strcmpi(header,'SNPEFF_TRANSCRIPT_ID')), {'entrez'}, 'transcript');
tableDnaRna(:,end+1) = {''};
header{end+1} = 'Desc';
tableDnaRna(~isnan(entrez), end) = GENOMEFUNC.idQuery(NCBI, entrez(~isnan(entrez)), 'desc');
keepcol = ~all(ismember(tableDnaRna, {'NA','NaN','NONE'}),1);
tableDnaRna = [header(keepcol); tableDnaRna(:, keepcol)];
tabwrite(sprintf('tables/P%s, preliminary, abs z more than 2 from diff in altDnaRnaAF.txt', curpatient.patient), tableDnaRna);
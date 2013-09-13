
figdir = 'figures/rnavar/';
savefig = true;
set(gcf, 'position',  [1678, 1155, 875, 703]);


patienttb = {'patient5', '132540-10N', '132540-1T', 'primary'; ...
    'patient5', '132540-10N', '137064-1T', 'meta1'; ...
    'patient5', '132540-10N', '139508-1T', 'meta2'; ...
    'patient10', '138381-4N', '138381-2T', 'meta1'; ...
    'patient4', '130342-7N', '130342-1T', 'primary'; ...
    'patient4', '130342-7N', '130342-13T', 'meta1'; ...
    'patient4', '130342-7N', '130342-18T', 'meta2'; ...
    'patient4', '130342-7N', '130342-21T', 'meta3'; ...
    'patient4', '130342-7N', '135098-7T', 'meta4'} ;

npair = size(patienttb, 1);

callerset = {'mutect', 'mutect'; 'varscan', 'varscan'; ...
    'sniper', 'sniper'; 'mutect', 'mutect_MQ'; 'sniper', 'sniper_-J-s0.01'; ...
    'mutect', 'mutect_q50'};

datadir = 'data/';

%%
% index the location of all variant calls from RNA
for calleridx = 1:size(callerset,1)
    if isempty(strfind(callerset{calleridx,2}, '_'))
        keyfd = [callerset{calleridx,2} '_def'];
    else
        keyfd = strrep(strrep(callerset{calleridx,2}, '-',''), '.','');
    end
    if ~isfield(RnaVarLoc, keyfd)
        RnaVarLoc.(keyfd).lockey = cell(npair, 2); %mdup, and no-mdup
        RnaVarLoc.(keyfd).filter = cell(npair, 2);
        RnaVarLoc.(keyfd).GT = cell(npair, 2);
        if ~isempty(strfind(keyfd, 'varscan')) || ~isempty(strfind(keyfd, 'sniper'))
            RnaVarLoc.(keyfd).SS = cell(npair, 2);
        end
    else
        curnpair = size(RnaVarLoc.(keyfd).lockey, 1);
        RnaVarLoc.(keyfd).lockey(curnpair+1:npair, :) = cell(npair-curnpair, 2);
        RnaVarLoc.(keyfd).filter(curnpair+1:npair, :) = cell(npair-curnpair, 2);
        RnaVarLoc.(keyfd).GT(curnpair+1:npair, :) = cell(npair-curnpair, 2);
        RnaVarLoc.(keyfd).SS(curnpair+1:npair, :) = cell(npair-curnpair, 2);
    end
    
    for pairidx = 4:npair %1:npair
        if ~isempty(strfind(callerset{calleridx,2}, 'mutect_'))
            tmp = textscan(callerset{calleridx,2}, '%s', 'delimiter', '_');
            callersetname = [ tmp{1}{2} '.' tmp{1}{1}];
        else
            callersetname = strrep(callerset{calleridx,2},'_','.');
        end
        datafilename = [datadir ...
            sprintf('rnavar.%s_%s.mdup.%s.mat',patienttb{pairidx,2}, ...
            patienttb{pairidx,3}, callersetname)];
        if exist(datafilename, 'file')
            v = loadStructData(datafilename);
            RnaVarLoc.(keyfd).lockey{pairidx,1} = strkey(v.lockey(), ...
                'add', 'saveto', '~/Projects/Simon/data/DICTIONARY.mat', 'savelater');
            [~, colidx] = ismember('FILTER', v.attrName_cell);
            if colidx ~= 0
                RnaVarLoc.(keyfd).filter{pairidx,1} = strkey(strrep(v.variantAttr_cell(:, colidx),'.',''), ...
                    'add', 'saveto', '~/Projects/Simon/data/DICTIONARY.mat', 'savelater');
            end
            RnaVarLoc.(keyfd).GT{pairidx, 1} = v.formatData.GT;
            if isfield(RnaVarLoc.(keyfd), 'SS')
                if isfield(v.formatData, 'SS')
                    RnaVarLoc.(keyfd).SS{pairidx,1} = v.formatData.SS;
                else
                    RnaVarLoc.(keyfd).SS{pairidx,1} = v.variantAttr_mtx(:, strcmp(v.attrName_mtx, 'SS'));
                end
            end
        end
        
        datafilename = [datadir ...
            sprintf('rnavar.%s_%s.%s.mat',patienttb{pairidx,2}, ...
            patienttb{pairidx,3}, callersetname)];
        if exist(datafilename, 'file')
            v = loadStructData(datafilename);
            RnaVarLoc.(keyfd).lockey{pairidx,2} = strkey(v.lockey(), ...
                'add', 'saveto', '~/Projects/Simon/data/DICTIONARY.mat', 'savelater');
            [~, colidx] = ismember('FILTER', v.attrName_cell);
            if colidx ~= 0
                RnaVarLoc.(keyfd).filter{pairidx,2} = strkey(strrep(v.variantAttr_cell(:, colidx),'.',''), ...
                    'add', 'saveto', '~/Projects/Simon/data/DICTIONARY.mat', 'savelater');
            end
            RnaVarLoc.(keyfd).GT{pairidx, 2} = v.formatData.GT;
            if isfield(RnaVarLoc.(keyfd), 'SS')
                if isfield(v.formatData, 'SS')
                    RnaVarLoc.(keyfd).SS{pairidx,2} = v.formatData.SS;
                else
                    RnaVarLoc.(keyfd).SS{pairidx,2} = v.variantAttr_mtx(:, strcmp(v.attrName_mtx, 'SS'));
                end
            end
        end
    end
end
% strkey('savenow');
% save data/RnaVarLoc.mat RnaVarLoc
%%

pids = {'132540', '138381'};
DnaVarLoc.pid = pids;
DnaVarLoc.lockey = cell(2,1);
DnaVarLoc.call = cell(2,1);
DnaVarLoc.caller = {'mutect', 'sniper', 'varscan', 'gatkIndel'};
DnaVarLoc.sniper_SSC = cell(2,1);
DnaVarLoc.varscan_SSC = cell(2,1);
for pidx = 1:length(pids)
    pdata = loadStructData( [datadir sprintf('patient.%s.obj2.q50.mat', pids{pidx})]);
    DnaVarLoc.lockey{pidx} = strkey(pdata.variant.lockey(false), ...
        'add', 'saveto', '~/Projects/Simon/data/DICTIONARY.mat', 'savelater');
    DnaVarLoc.call{pidx} = false(length(DnaVarLoc.lockey{pidx}), length(DnaVarLoc.caller));
    for cidx = 1:length(DnaVarLoc.caller)
        DnaVarLoc.call{pidx}(:, cidx) = ~cellfun(@isempty, strfind(pdata.variant.sets, DnaVarLoc.caller{cidx}));
    end
    DnaVarLoc.sniper_SSC{pidx} = pdata.variant.Sniper_SSC(:,1);
    DnaVarLoc.varscan_SSC{pidx} = pdata.variant.Varscan_SSC;
end
strkey('savenow');
% save data/DnaVarLoc.mat DnaVarLoc
%%
% history of snv, patient 5's primary and metastatic samples, one figure each caller,
% 1x2 (mdup) or 2x2 (mdup x para)
ucaller = unique(callerset(:,1));
callerparaset = fieldnames(RnaVarLoc);
pairidx = find(strcmp(patienttb(:,1), 'patient5'));
label = patienttb(pairidx,4);
mduplabel = {'mdup', 'no-mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');

for uci = 1:length(ucaller)
    fdidx = find(~cellfun(@isempty, strfind(callerparaset, ucaller{uci})));        
    clf;
    for paraidx = 1:length(fdidx)
        for mdupidx = 1:2        
            subplot(2,2,(paraidx-1)*2+mdupidx);
            
            tmp = RnaVarLoc.(callerparaset{fdidx(paraidx)}).lockey(pairidx, mdupidx);
            for tmpidx = 1:length(tmp)
                if strcmp(ucaller{uci}, 'varscan')
                    tmp{tmpidx} = tmp{tmpidx}( ...
                        RnaVarLoc.(callerparaset{fdidx(paraidx)}).filter{pairidx(tmpidx),mdupidx}==passkey & ...
                        RnaVarLoc.(callerparaset{fdidx(paraidx)}).SS{pairidx(tmpidx),mdupidx} ==2 );
                elseif strcmp(ucaller{uci}, 'mutect')
                    tmp{tmpidx} = tmp{tmpidx}(RnaVarLoc.(callerparaset{fdidx(paraidx)}).filter{pairidx(tmpidx),mdupidx}==passkey);
                else
                    tmp{tmpidx} = tmp{tmpidx}( ...
                        RnaVarLoc.(callerparaset{fdidx(paraidx)}).SS{pairidx(tmpidx),mdupidx}(:,2) ==2 );
                end
            end
            vi = ~cellfun(@isempty, tmp);
            if any(vi)
                fixvenn(tmp(vi), ...
                    'label', label(vi), 'numfontsize', 10, 'labelfontsize',12);
            else
                cla
            end
            
            title([strrep(callerparaset{fdidx(paraidx)},'_', ' '), ...
                ' ', mduplabel{mdupidx}], 'fontsize', 12);
        end
    end
    if savefig
        filename = [figdir sprintf('venn patient5, somatic in primary-meta, %s', ucaller{uci})];
        plot2svg([filename '.svg'], gcf);
    end
end

%%
% effect of joint in sniper, two figure (mdup or no mdup), 2x2 (3 patient5 + 1
% patient 10)

callerparaset = fieldnames(RnaVarLoc);
fdidx = find(~cellfun(@isempty, strfind(callerparaset, 'sniper')));
setlabel = callerparaset(fdidx);
setlabel = strrep(setlabel, 'sniper_def', 'sniper');
setlabel = strrep(setlabel, 'sniper_Js001', 'sniper joint');

for mdupidx = 1:2
    clf;
    for pairidx = 1:size(patienttb, 1)        
        subplot(2,2,pairidx);
        tmp = cell(1,2);
        for paraidx = 1:length(fdidx)            
            tmp{paraidx} = RnaVarLoc.(callerparaset{fdidx(paraidx)}).lockey{pairidx,mdupidx};
            tmp{paraidx} = tmp{paraidx}(RnaVarLoc.(callerparaset{fdidx(paraidx)}).SS{pairidx,mdupidx}(:,2) == 2);
        end                                
        if length(tmp) > 3
            fixvenn(tmp, ...
                'label', setlabel, 'numfontsize', 7, 'labelfontsize',12);
            title([patienttb{pairidx,1}, ' ', patienttb{pairidx,4}], 'fontsize', 12);
        else
            fixvenn(tmp, ...
                'label', setlabel, 'numfontsize', 10, 'labelfontsize',12);
            title([patienttb{pairidx,1}, ' ', patienttb{pairidx,4}], 'fontsize', 12);
        end
    end
    if savefig
        filename = [figdir sprintf('venn effect of para, sniper')];
        plot2svg([filename '.svg'], gcf);
    end
end

%%
% overlap between samples: P5 x 3 + P10, one figure each caller, 1x2 (mdup) or 2x2 (mdup x para)
ucaller = unique(callerset(:,1));
callerparaset = fieldnames(RnaVarLoc);
mduplabel = {'mdup', 'no-mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');
label = strcat(patienttb(:,1), '-', patienttb(:,4));
label = strrep(label, 'patient', 'p');
for uci = 1:length(ucaller)
    fdidx = find(~cellfun(@isempty, strfind(callerparaset, ucaller{uci})));    
    if strcmp(ucaller{uci}, 'mutect')
        fdidx(strcmp(callerparaset(fdidx), 'mutect_def')) = [];            
    end
    clf;
    for paraidx = 1:length(fdidx)
        for mdupidx = 1:2        
            subplot(2,2,(paraidx-1)*2+mdupidx);
            
            tmp = RnaVarLoc.(callerparaset{fdidx(paraidx)}).lockey(:, mdupidx);
            for tmpidx = 1:length(tmp)
                if strcmp(ucaller{uci}, 'varscan')
                    tmp{tmpidx} = tmp{tmpidx}( ...
                        RnaVarLoc.(callerparaset{fdidx(paraidx)}).filter{tmpidx,mdupidx}==passkey & ...
                        RnaVarLoc.(callerparaset{fdidx(paraidx)}).SS{tmpidx,mdupidx} ==2 );
                elseif strcmp(ucaller{uci}, 'mutect')
                    tmp{tmpidx} = tmp{tmpidx}(RnaVarLoc.(callerparaset{fdidx(paraidx)}).filter{tmpidx,mdupidx}==passkey);
                else
                    tmp{tmpidx} = tmp{tmpidx}(RnaVarLoc.(callerparaset{fdidx(paraidx)}).SS{tmpidx,mdupidx}(:,2) ==2);
                end
                %                     tmp{tmpidx} = tmp{tmpidx}(RnaVarLoc.(callerparaset{fdidx(paraidx)}).filter{tmpidx,mdupidx}~=removekey);
            end
            
            fixvenn(tmp, ...
                'label', label, 'numfontsize', 12, 'labelfontsize',14);            
            title([strrep(callerparaset{fdidx(paraidx)},'_', ' '), ...
                ' ', mduplabel{mdupidx}], 'fontsize', 14);
        end
    end
    if savefig
        filename = [figdir sprintf('venn somatic overlap over all tumor samples, %s', ucaller{uci})];
        plot2svg([filename '.svg'], gcf);
    end
end
%%
% overlap between callers, one figure for each para
ucaller = {'mutect_MQ', 'sniper_def', 'varscan_def'};
mduplabel = {'mdup', 'no mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');
label = {'mutect', 'sniper', 'varscan'};

sslabel = {'somatic', 'LOH'};
ssStatus = 3;
for mdupidx = 1:2
    clf        
    for pairidx = 1:size(patienttb,1)
    tmp = cell(1,3);
        
        for i = 1:length(ucaller)
            if ~isempty(strfind(ucaller{i}, 'mutect'))
                rnavi = RnaVarLoc.(ucaller{i}).filter{pairidx, mdupidx}==passkey;                
            elseif ~isempty(strfind(ucaller{i}, 'varscan'))
                rnavi = RnaVarLoc.(ucaller{i}).filter{pairidx, mdupidx}==passkey;
                rnavi = rnavi & RnaVarLoc.(ucaller{i}).SS{pairidx, mdupidx} == ssStatus;                
            else
                rnavi = RnaVarLoc.(ucaller{i}).SS{pairidx, mdupidx}(:,2) == ssStatus ;                
            end
            tmp{i} = RnaVarLoc.(ucaller{i}).lockey{pairidx, mdupidx}(rnavi);
        end
        
        subplot(2,2,pairidx);
        fixvenn(tmp, ...
            'label', label, 'numfontsize', 12, 'labelfontsize',14);
        title([patienttb{pairidx,1}, ' ', patienttb{pairidx,4}, ' ', mduplabel{mdupidx}], 'fontsize', 14);
    end
    if savefig
        filename = [figdir sprintf('venn callers, %s, %s', mduplabel{mdupidx}, sslabel{ssStatus-1})];
        plot2svg([filename '.svg'], gcf);
    end
end


%%
% effect of MQ, one figure for mdup or no-mdup, mutect only, 2x2 (tumor samples)

mduplabel = {'mdup', 'no mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');
label = {'no MQ', 'MQ-pass', 'MQ-reject'};
for mdupidx = 1:2
    for pairidx = 1:size(patienttb,1)
        subplot(2,2,pairidx)
        pass = RnaVarLoc.mutect_MQ.filter{pairidx, mdupidx} == passkey;
        tmp = { RnaVarLoc.mutect_def.lockey{pairidx, mdupidx}, ...
            RnaVarLoc.mutect_MQ.lockey{pairidx, mdupidx}(pass), ...
            RnaVarLoc.mutect_MQ.lockey{pairidx, mdupidx}(~pass) };
        fixvenn(tmp, 'label', label, 'numfontsize', 12, 'labelfontsize', 14);
        title([patienttb{pairidx,1}, ' ', mduplabel{mdupidx}], 'fontsize', 14);
    end
    if savefig
        filename = [figdir sprintf('venn effect of MQ, %s', mduplabel{mdupidx})];
        plot2svg([filename '.svg'], gcf);
    end
end

%%
% effect of q50: uniquely mapped reads 

mdupidx = 2;
label = {'mutect', 'mutect MQ', 'mutect q50'};
mduplabel = {'mdup', 'no mdup'};
for pairidx = 1:size(patienttb, 1)
    subplot(2,2,pairidx)
    %pass = RnaVarLoc.mutect_MQ.filter{pairidx, mdupidx} == passkey;
    tmp = { RnaVarLoc.mutect_def.lockey{pairidx, mdupidx}, ...
        RnaVarLoc.mutect_MQ.lockey{pairidx, mdupidx}, ...
        RnaVarLoc.mutect_q50.lockey{pairidx, mdupidx}};
    fixvenn(tmp, 'label', label, 'numfontsize', 12, 'labelfontsize', 14);
    title([patienttb{pairidx,1}, ' ', mduplabel{mdupidx}], 'fontsize', 14);
end
    
if savefig
    filename = [figdir sprintf('venn effect of q50, %s', mduplabel{mdupidx})];
    plot2svg([filename '.svg'], gcf);
end

%%
% mutect rejected vs other callers, one figure each para

ucaller = {'mutect_MQ', 'sniper_def', 'varscan_def'};
mduplabel = {'mdup', 'no mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');
label = {'mutect-pass', 'sniper', 'varscan', 'mutect-reject'};

sslabel = {'somatic', 'LOH'};
ssStatus = 3;
for mdupidx = 1:2
    clf        
    for pairidx = 1:size(patienttb,1)
    tmp = cell(1,4);
        
        for i = 1:length(ucaller)
            if ~isempty(strfind(ucaller{i}, 'mutect'))
                rnavi = RnaVarLoc.(ucaller{i}).filter{pairidx, mdupidx}==passkey; 
                tmp{i} = RnaVarLoc.(ucaller{i}).lockey{pairidx, mdupidx}(rnavi);
                tmp{4} = RnaVarLoc.(ucaller{i}).lockey{pairidx, mdupidx}(~rnavi);
            elseif ~isempty(strfind(ucaller{i}, 'varscan'))
                rnavi = RnaVarLoc.(ucaller{i}).filter{pairidx, mdupidx}==passkey;
                rnavi = rnavi & RnaVarLoc.(ucaller{i}).SS{pairidx, mdupidx} == ssStatus;                
                tmp{i} = RnaVarLoc.(ucaller{i}).lockey{pairidx, mdupidx}(rnavi);
            else
                rnavi = RnaVarLoc.(ucaller{i}).SS{pairidx, mdupidx}(:,2) == ssStatus ;                
                tmp{i} = RnaVarLoc.(ucaller{i}).lockey{pairidx, mdupidx}(rnavi);
            end
        end
        
        subplot(2,2,pairidx);
        fixvenn(tmp, ...
            'label', label, 'numfontsize', 12, 'labelfontsize',14);
        title([patienttb{pairidx,1}, ' ', patienttb{pairidx,4}, ' ', mduplabel{mdupidx}, ' ', sslabel{ssStatus-1}], 'fontsize', 14);
    end
    if savefig
        filename = [figdir sprintf('venn callers with mutect-split, %s, %s', mduplabel{mdupidx}, sslabel{ssStatus-1})];
        plot2svg([filename '.svg'], gcf);
    end
end

%% 
% DNA vs RNA overlap of each caller, one figure each tumor-sample/mdup, 2x2 (3
% callers + sniper-J)

ucaller = fieldnames(RnaVarLoc);
ucaller(strcmp(ucaller, 'mutect_def')) = [];
mduplabel = {'mdup', 'no mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');

sslabel = {'somatic', 'LOH'};
ssStatus = 3;
for pidx = 1:length(DnaVarLoc.pid)
    dnavi = cellfun(@isempty, strfind(strkey(DnaVarLoc.lockey{pidx},'lookup', '~/Projects/Simon/data/DICTIONARY.mat'), 'chrM'));
    for mdupidx = 1:2
        clf
        rnadataidx = find(~cellfun(@isempty, strfind(patienttb(:,3), DnaVarLoc.pid{pidx})));
        
        for i = 1:length(ucaller)            
            rnavi = cellfun(@isempty, strfind(strkey(RnaVarLoc.(ucaller{i}).lockey{rnadataidx, mdupidx},'lookup', '~/Projects/Simon/data/DICTIONARY.mat'), 'chrM'));
            if ~isempty(strfind(ucaller{i}, 'mutect'))
                rnavi = rnavi & RnaVarLoc.(ucaller{i}).filter{rnadataidx, mdupidx}==passkey;
                caller = 'mutect';
            elseif ~isempty(strfind(ucaller{i}, 'varscan'))
                rnavi = rnavi & RnaVarLoc.(ucaller{i}).filter{rnadataidx, mdupidx}==passkey;
                rnavi = rnavi & RnaVarLoc.(ucaller{i}).SS{rnadataidx, mdupidx} ==ssStatus;
                caller = 'varscan';
            else
                rnavi = RnaVarLoc.(ucaller{i}).SS{rnadataidx, mdupidx}(:,2) ==ssStatus;
                caller = 'sniper';
            end
            subplot(2,2,i);
            
            fixvenn({DnaVarLoc.lockey{pidx}( dnavi & DnaVarLoc.call{pidx}(:, strcmp(DnaVarLoc.caller, caller))), ...
                RnaVarLoc.(ucaller{i}).lockey{rnadataidx, mdupidx}( rnavi ) }, ...
                'label', {'DNA', 'RNA'}, 'numfontsize', 12, 'labelfontsize',14);            
            title([patienttb{rnadataidx,1}, ' ', strrep(ucaller{i}, '_', ' '), ' ', mduplabel{mdupidx}], 'fontsize', 14);
        end
        if savefig
            filename = [figdir sprintf('venn DNA and RNA, %s, %s, %s', DnaVarLoc.pid{pidx}, mduplabel{mdupidx}, sslabel{ssStatus-1})];
            plot2svg([filename '.svg'], gcf);
        end
    end
end


%%
% DNA vs RNA overlap, across callers, one figure each tumor-sample, 3x3
% (DNA 3 ind callers x RNA 3 ind callers)
ucaller = unique(callerset(:,1));
mduplabel = {'mdup', 'no mdup'};
for pidx = 1:length(DnaVarLoc.pid)
    dnavi = cellfun(@isempty, strfind(strkey(DnaVarLoc.lockey{pidx},'lookup', '~/Projects/Simon/data/DICTIONARY.mat'), 'chrM'));
    for mdupidx = 1:2
        clf
        rnadataidx = find(~cellfun(@isempty, strfind(patienttb(:,3), DnaVarLoc.pid{pidx})));
        for i = 1:length(ucaller)            
            for j = 1:length(ucaller)
                if strcmp(ucaller{j}, 'mutect')
                    rnafd = 'mutect_MQ';
                else
                    rnafd = [ucaller{j} '_def'];
                end
                rnavi = cellfun(@isempty, strfind(strkey(RnaVarLoc.(rnafd).lockey{rnadataidx, mdupidx},'lookup', '~/Projects/Simon/data/DICTIONARY.mat'), 'chrM'));
                if strcmp(ucaller{j}, 'mutect')
                    rnavi = rnavi & RnaVarLoc.(rnafd).filter{rnadataidx, mdupidx}==passkey;                    
                elseif strcmp(ucaller{j}, 'varscan')
                    rnavi = rnavi & RnaVarLoc.(rnafd).filter{rnadataidx, mdupidx}==passkey;
                    rnavi = rnavi & RnaVarLoc.(rnafd).SS{rnadataidx, mdupidx} ==2;   
                else
                    rnavi = RnaVarLoc.(rnafd).SS{rnadataidx, mdupidx}(:,2) ==2;   
                end
                subplot(3,3,(i-1)*3+j);
                
                fixvenn({DnaVarLoc.lockey{pidx}( dnavi & DnaVarLoc.call{pidx}(:, strcmp(DnaVarLoc.caller, ucaller{i}))), ...
                    RnaVarLoc.(rnafd).lockey{rnadataidx, mdupidx}( rnavi ) }, ...
                    'label', {'DNA', 'RNA'}, 'numfontsize', 12, 'labelfontsize',14);
                if i == 1
                    title([patienttb{rnadataidx,1}, ' ', strrep(ucaller{j}, '_', ' '), ' ', mduplabel{mdupidx}], 'fontsize', 14);
                end
            end
        end
        if savefig
            filename = [figdir sprintf('venn DNA and RNA, across callers, %s, %s', DnaVarLoc.pid{pidx}, mduplabel{mdupidx})];
            plot2svg([filename '.svg'], gcf);
        end
    end
end

%%
% DNA vs RNA, DNA-1/2/3-callers overlap RNA, 4 sets: DNA + 3 callers,  2x2
% (DNA-union, DNA-any2, DNA-intersection)

%%
% DNA-intersect vs RNA intersect, 2x2 (tumor samples), perhaps RNA-editing

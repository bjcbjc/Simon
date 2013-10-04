
plotfig.D2sR2s = false;
plotfig.D3cR1c = false;


if ~exist('RnaVarLoc', 'var')
    load data/RnaVarLoc.mat
end
if ~exist('DnaVarLoc', 'var')
    load data/DnaVarLoc.mat
end


fnlist = {'132540-1-T--132540-10-N.indel.varscan.v2.3.6.vcf.loc', ...
    '138381-2-T--138381-4-N.indel.varscan.v2.3.6.vcf.loc'};
for fni = 1:length(fnlist)
    if ismember(fnlist{fni}, DnaVarLoc.fns)
        i = find(strcmp(DnaVarLoc.fns, fnlist{fni}));
        j = find(strcmp(DnaVarLoc.fns, strrep(fnlist{fni}, 'indel', 'snv')));
        if isempty(intersect(DnaVarLoc.locindex{i}, DnaVarLoc.locindex{j})) 
            DnaVarLoc.locindex{j} = [DnaVarLoc.locindex{j}; DnaVarLoc.locindex{i}];
            DnaVarLoc.loc{j} = [DnaVarLoc.loc{j}; DnaVarLoc.loc{i}];
            DnaVarLoc.fns(i) = [];
            DnaVarLoc.loc(i) = [];
            DnaVarLoc.locindex(i) = [];
        end
    end
end

if ~isfield(DnaVarLoc, 'pid')
    DnaVarLoc.pid = repmat({'132540'}, length(DnaVarLoc.fns),1);
    DnaVarLoc.pid(~cellfun(@isempty, strfind(DnaVarLoc.fns, '138381'))) = {'138381'};
end
if ~isfield(DnaVarLoc, 'caller')
    DnaVarLoc.caller = cell(length(DnaVarLoc.fns), 1);
    caller = {'varscan', 'somaticsniper', 'joint.somaticsniper', 'mutect'};
    for i = 1:length(caller)
        DnaVarLoc.caller(~cellfun(@isempty, strfind(DnaVarLoc.fns, caller{i}))) = caller(i);
    end
end
if ~isfield(DnaVarLoc, 'somatic')
    if ~exist('DnaVarLoc_somatic', 'var')
        load data/DnaVarLoc.somatic.mat
    end
    DnaVarLoc.somatic = cell(length(DnaVarLoc.fns),1);
    for i = 1:length(DnaVarLoc.fns)
        pidx = find(strcmp(DnaVarLoc_somatic.pid, DnaVarLoc.pid{i}));
        if strcmp(DnaVarLoc.caller{i}, 'joint.somaticsniper') || ...
            strcmp(DnaVarLoc.caller{i}, 'somaticsniper')
            cidx = find(strcmp(DnaVarLoc_somatic.caller, 'sniper'));
        else
            cidx = find(strcmp(DnaVarLoc_somatic.caller, DnaVarLoc.caller{i}));
        end
        
        DnaVarLoc.somatic{i} = sparse(ismember(DnaVarLoc.locindex{i}, ...
            DnaVarLoc_somatic.locindex{pidx}(DnaVarLoc_somatic.call{pidx}(:,cidx)) ));
    end
    clear DnaVarLoc_somatic
end

if ~isfield(DnaVarLoc, 'exonutr')
    if ~exist('GENCODE_exon_utr', 'var')
        load GENCODE_exon_utr.mat
    end    
    DnaVarLoc.exonutr = cell(length(DnaVarLoc.fns),1);
    for i = 1:length(DnaVarLoc.fns)
        DnaVarLoc.exonutr{i} = false(length(DnaVarLoc.locindex{i}),1);
        for chridx = 1:25
            nblk = ceil(max(DnaVarLoc.loc{i}( DnaVarLoc.loc{i}(:,1)==chridx,2)) /1e7);
            if isempty(nblk), continue; end
            for blkidx = 1:nblk
                lmin = (blkidx-1)*1e7+1;
                lmax = blkidx*1e7;
                vidx = DnaVarLoc.loc{i}(:,1) == chridx & ...
                    DnaVarLoc.loc{i}(:,2) >= lmin & DnaVarLoc.loc{i}(:,2) <= lmax;
                antidx = GENCODE_exon_utr.chrm == chridx & ...
                    GENCODE_exon_utr.start >= lmin & GENCODE_exon_utr.end <= lmax;
                DnaVarLoc.exonutr{i}(vidx) = any(bsxfun(@ge, DnaVarLoc.loc{i}(vidx,2), GENCODE_exon_utr.start(antidx)') ...
                    & bsxfun(@le, DnaVarLoc.loc{i}(vidx,2), GENCODE_exon_utr.end(antidx)'), 2);
            end
        end
    end
    clear GENCODE_exon_utr
end


figdir = 'figures/rnavar/overlapdna/';
savefig = true;
set(gcf, 'position',  [1678, 1155, 875, 703], 'visible', 'off', 'paperpositionmode', 'auto');
if ~exist('patienttb', 'var')
    patienttb = loadStructData('data/metaTB.mat');
    patienttb = patienttb.normalTumorPair;
end
npair = size(patienttb, 1);
upid = unique(DnaVarLoc.pid);

mduplabel = {'mdup', 'no mdup'};
passkey = strkey('PASS','lookup', '~/Projects/Simon/data/DICTIONARY.mat');
sslabel = {'somatic', 'LOH'};
ssStatus = 2;
chrMStartIndex = gloc2index([25, 1]);

if plotfig.D2sR2s
    ucaller = unique(DnaVarLoc.caller);
    for mdupidx = 1:2
        for pidx = 1:length(upid)
            clf
            for i = 1:length(ucaller)
                Didx = find(strcmp(DnaVarLoc.pid, upid{pidx}) & strcmp(DnaVarLoc.caller, ucaller{i}));
                dnavi = DnaVarLoc.locindex{Didx} >= chrMStartIndex;
                
                if strcmp(upid{pidx}, '132540')
                    Rpidx =find(strcmp(patienttb(:,3), '132540-1T'));
                else
                    Rpidx = find(strcmp(patienttb(:,3), '138381-2T'));
                end
                
                
                if strcmp(ucaller{i}, 'mutect')
                    Rfd = 'mutect_MQ100';
                    rnavi = RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx} >= chrMStartIndex;
                    rnasomatic = RnaVarLoc.(Rfd).filter{Rpidx, mdupidx}==passkey;
                elseif strcmp(ucaller{i}, 'joint.somaticsniper')
                    Rfd = 'sniper_Js001';
                    rnavi = RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx} >= chrMStartIndex;
                    rnasomatic = RnaVarLoc.(Rfd).SS{Rpidx, mdupidx}(:,2) ==ssStatus;
                elseif strcmp(ucaller{i}, 'somaticsniper')
                    Rfd = 'sniper_def';
                    rnavi = RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx} >= chrMStartIndex;
                    rnasomatic = RnaVarLoc.(Rfd).SS{Rpidx, mdupidx}(:,2) ==ssStatus;
                elseif strcmp(ucaller{i}, 'varscan')
                    Rfd = 'varscan_def';
                    rnavi = RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx} >= chrMStartIndex;
                    rnasomatic = RnaVarLoc.(Rfd).filter{Rpidx, mdupidx}==passkey & ...
                        RnaVarLoc.(Rfd).SS{Rpidx, mdupidx} ==ssStatus;
                else
                    Rfd = '';
                end
                
                subplot(2,2,i);
                
                fixvenn({DnaVarLoc.locindex{Didx}(~DnaVarLoc.somatic{Didx} & dnavi ), ...
                    DnaVarLoc.locindex{Didx}(DnaVarLoc.somatic{Didx} & dnavi ), ...
                    RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx}( rnasomatic & rnavi ), ...
                    RnaVarLoc.(Rfd).locindex{Rpidx, mdupidx}( ~rnasomatic & rnavi )}, ...
                    'label', {'DNA-rest', 'DNA-somatic', 'RNA-somatic', 'RNA-rest'}, 'numfontsize', 12, 'labelfontsize',14);
                title([patienttb{Rpidx,1}, ' ', ucaller{i}, ' ', mduplabel{mdupidx}], 'fontsize', 14);
            end
            if savefig
                filename = [figdir sprintf('venn DNA and RNA, %s, %s', DnaVarLoc.pid{Didx}, mduplabel{mdupidx})];
                plot2svg([filename '.svg'], gcf);
            end
        end
    end
end

if plotfig.D3cR1c
    mdupidx = 1;
    ucaller = {'mutect_MQ100', 'varscan_def', 'sniper_def', 'sniper_Js001'};
    for pidx = 1:length(upid)

        Didx = find(strcmp(DnaVarLoc.pid, upid{pidx}));
        [~, si] = ismember({'mutect', 'varscan','somaticsniper'}, DnaVarLoc.caller(Didx));
        Didx = Didx(si);
        tmp = cell(1,4);
        for i = 1:3
            dnavi = DnaVarLoc.locindex{Didx(i)} >= chrMStartIndex;            
            tmp{i} = DnaVarLoc.locindex{Didx(i)}(dnavi);
        end
            
        for i = 1:length(ucaller)
            if strcmp(upid{pidx}, '132540')
                Rpidx =find(strcmp(patienttb(:,3), '132540-1T'));
            else
                Rpidx = find(strcmp(patienttb(:,3), '138381-2T'));
            end
            
            rnavi = RnaVarLoc.(ucaller{i}).locindex{Rpidx, mdupidx} >= chrMStartIndex;
            
            subplot(2,2,i);
            tmp{4} = RnaVarLoc.(ucaller{i}).locindex{Rpidx, mdupidx}(rnavi);
            fixvenn(tmp, ...
                'label', {'DNA-mutect', 'DNA-varscan', 'DNA-sniper','RNA'}, 'numfontsize', 9, 'labelfontsize',12);
                title([patienttb{Rpidx,1}, ' ', ucaller{i}, ' ', mduplabel{mdupidx}], 'fontsize', 12);
        end
        if savefig
            filename = [figdir sprintf('venn DN-3-callers and RNA, %s, %s', DnaVarLoc.pid{Didx}, mduplabel{mdupidx})];
            plot2svg([filename '.svg'], gcf);
        end
    end
end
close(gcf);
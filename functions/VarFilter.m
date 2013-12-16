classdef VarFilter < handle
    
    methods (Static)
        function filter = TrueSet(varstruct, callerfd, method, mincov, snpirFilterIdx, UseSomaticFilter)                    
            if ~iscell(method)
                method = {method};
            end
            if ~iscell(callerfd)
                callerfd = {callerfd};
            end
            nloc = length(varstruct.locidx);
            filter = false(nloc, length(method));
            mtx = false(nloc, length(callerfd));
            for i = 1:length(callerfd)
                mtx(:,i) = VarFilter.CovAndNumReadFilter(varstruct, callerfd{i}, mincov, snpirFilterIdx);
                mtx(:,i) = mtx(:,i) & VarFilter.CumValidAfterFilter(varstruct, callerfd{i}, snpirFilterIdx, UseSomaticFilter);
            end
            %method: {'union', 'any2', 'all'};
            for i = 1:length(method)
                if strcmpi(method{i}, 'union')
                    filter(:,i) = any(mtx, 2);
                elseif strcmpi(method{i}, 'any2')
                    filter(:,i) = sum(mtx~=0, 2) > 1;
                elseif strcmpi(method{i}, 'any3')
                    filter(:,i) = sum(mtx~=0, 2) > 2;
                elseif strcmpi(method{i}, 'all')
                    filter(:,i) = all(mtx, 2);                
                end                    
            end
            filter = sparse(filter);
        end
        function filter = NumReadFilter(varstruct, callerfd, mincov, snpirFilterIdx)
            if iscell(callerfd)
                error('callerfd cannot be cell');
            end
            nSnpirFilter = size(varstruct.(callerfd).validAfterFilter,2);
            if nargin < 4
                snpirFilterIdx = 1:nSnpirFilter;
            end            
            filter = false(length(varstruct.locidx), length(snpirFilterIdx));
            cov = varstruct.(callerfd).numReadRef + varstruct.(callerfd).numReadAlt;
            for i = 1:length(snpirFilterIdx)
                idx = min(snpirFilterIdx(i), nSnpirFilter);
                if all(cov(:,idx)==0)
                    idx = 3;
                end
                filter(:,i) = cov(:,idx) >= mincov;
            end
            filter = sparse(filter);
        end
        function filter = CumValidAfterFilter(varstruct, callerfd, filteridx, UseSomaticFilter)
            if nargin < 4, UseSomaticFilter = true; end
            if length(filteridx) == 1
                if UseSomaticFilter
                    filter = all(varstruct.(callerfd).validAfterFilter(:, 1:max(2,filteridx)),2);
                else
                    filter = all(varstruct.(callerfd).validAfterFilter(:, setdiff(1:filteridx,2)),2);
                end
                return
            end
            filter = false(length(varstruct.locidx), length(filteridx));
            nSnpirFilter = size(varstruct.(callerfd).validAfterFilter,2);
            for i = 1:length(filteridx)
                idx = min(filteridx(i), nSnpirFilter);
                if ~UseSomaticFilter
                    filter(:, i) = all(varstruct.(callerfd).validAfter(:, setdiff(1:idx,2)),2);
                else
                    filter(:, i) = all(varstruct.(callerfd).validAfter(:, 1:max(2,idx)),2);
                end
            end
            filter = sparse(filter);
        end
        function filter = CovFilter(varstruct, callerfd, mincov)
            if mincov == 0
                filter = true(length(varstruct.locidx), 1);
                return
            end
            if ~ismember(mincov, varstruct.covlimit)
                error('mincov is not in varstruct.covlimit this filter is to filter variants by coverage of other sample');
            end
            filter = varstruct.(callerfd).covfilter(:, varstruct.covlimit==mincov);
        end
        function filter = CovAndNumReadFilter(varstruct, callerfd, mincov, snpirFilterIdx)
            % wrapper: call CovFilter and NumReadFilter
            filter = VarFilter.NumReadFilter(varstruct, callerfd, mincov, snpirFilterIdx);
            if isfield(  varstruct.(callerfd), 'covfilter')
                filter = filter & VarFilter.CovFilter(varstruct, callerfd, mincov);
            end
        end
    end    
    
end
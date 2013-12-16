function res = readSNPiRResult(fndir, fnbase)

    steps = {'','rmhex', 'rmsk', 'rmintron', 'rmhom', 'rmblat', 'rmedit'};
    
    f = fopen([fndir,fnbase, '.variant.vcf'], 'r');
    data = textscan(f, '%s %f %*[^\n]', 'commentstyle', '#');
    fclose(f);
    data{1} = numericchrm(data{1});
    res.locidx = gloc2index([data{1}, data{2}]);    
    nloc = length(res.locidx); 
    
    res.fndir = fndir;
    res.fnbase = fnbase;
    res.filterName = steps; 
    res.filterName{1} = 'qual';
    
    nfilter = length(res.filterName);
    res.validAfterFilter = false( nloc, nfilter );
    res.numReadRef = zeros( nloc, nfilter);
    res.numReadAlt = zeros( nloc, nfilter);
    for i = 1:length(steps)
        if strcmpi(steps{i}, 'rmedit')
            f = fopen([fndir, fnbase, strjoin(steps(1:i), '.'), '.bed'], 'r');
            data = textscan(f, '%s %*d %f %f,%f %*[^\n]');            
        else
            f = fopen([fndir, fnbase, strjoin(steps(1:i), '.'), '.txt'], 'r');
            data = textscan(f, '%s %f %f,%f %*[^\n]');            
        end
        [~, idx] = ismember( gloc2index( [numericchrm(data{1}), data{2}] ), res.locidx);
        res.validAfterFilter(idx, i) = true;
        res.numReadRef(idx,i) = data{3};
        res.numReadAlt(idx,i) = data{4};
        fclose(f);
    end
    res.validAfterFilter = sparse(res.validAfterFilter);
    res.numReadRef = sparse(res.numReadRef-res.numReadAlt);
    res.numReadAlt = sparse(res.numReadAlt);
end
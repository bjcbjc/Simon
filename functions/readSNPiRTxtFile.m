function res = readSNPiRTxtFile(fn)
    
    f = fopen(fn, 'r');
    if strcmpi(fn(end-3:end), '.bed')        
        data = textscan(f, '%s %*d %f %f,%f %*[^\n]');
    else        
        data = textscan(f, '%s %f %f,%f %*[^\n]');
    end
    fclose(f);
    
    res.locidx = gloc2index( [numericchrm(data{1}), data{2}] );    
    res.numReadRef = data{3};
    res.numReadAlt = data{4};
    
    rmi = isnan(res.locidx);
    if any(rmi)
        res.locidx(rmi) = [];
        res.numReadRef(rmi) = [];
        res.numReadAlt(rmi) = [];
    end    
    res.numReadRef = res.numReadRef-res.numReadAlt;
end


fns = dir('data/aCount*.mat');
for fnidx = 1:length(fns)
    aCountData = load(['data/' fns(fnidx).name]);
    aCountData = aCountData.aCountData;
    %filled in
    
    filename = regexp(aCountData.filename, '\S+/(\S)+$', 'tokens', 'once');
    filename = filename{1};
    
    tokens = regexp(filename, 'Sample_([\w\-]+)_\S+.r(\d+).d(\w+).q(\w+).s(\w)', 'tokens', 'once');
    aCountData.sample = tokens{1};
    tmp = textscan(tokens{1}, '%s', 'delimiter', '-');
    aCountData.patient = tmp{1}{1};
    aCountData.minRead = str2double(tokens{2});
    aCountData.maxDepth = str2double(strrep(tokens{3}, 'k', '000'));
    aCountData.minMapQ = str2double(tokens{4});
    if strcmpi(tokens{5}, 't')
        aCountData.stranded = true;
    else
        aCountData.stranded = false;
    end
    matfilename = sprintf('acount.%s.r%s.d%s.q%s.mat', ...
        aCountData.sample, tokens{2}, tokens{3}, tokens{4} );
    save(['data/' matfilename], 'aCountData');
    clear aCountData
end
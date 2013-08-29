
folder = 'DnaCount/';
fns = dir([folder '*d200k*.count']);
for i = range %1:length(fns)
    tic;
    aCountData = AlleleCountData([folder fns(i).name]);
    matfilename = aCountData.addparainfo();
    matfilename = strrep(matfilename, 'acount', 'acount.DNA');
    save([folder matfilename], 'aCountData');
    toc;
    clear aCountData matfilename
end
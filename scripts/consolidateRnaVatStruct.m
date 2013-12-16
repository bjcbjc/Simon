%load data/dnavar.138381-4N_138381-2T.withfilter.mat
fn = 'dnavar.132540-10N_132540-1T.withfilter.mat';
load(sprintf('data/%s',fn));
fds = fieldnames(dnavar);
fds(cellfun(@isempty, strfind(fds, 'bwa')) & cellfun(@isempty, strfind(fds, 'star')) ) = [];

DnaNumRead = dnavar;
for i = 1:length(fds)
    DnaNumRead.(fds{i}) = rmfield(DnaNumRead.(fds{i}), {'covfilter', 'validAfterFilter'});
end
fds(~cellfun(@isempty, strfind(fds,'star'))) = [];
DnaNumRead = rmfield(DnaNumRead, fds);
DnaNumRead = rmfield(DnaNumRead, 'covlimit');
save(sprintf('data/DnaNumRead.%s',strrep(fn,'dnavar.', '')), 'DnaNumRead');
clear DnaNumRead

fds = fieldnames(dnavar);
fds(cellfun(@isempty, strfind(fds, 'bwa')) & cellfun(@isempty, strfind(fds, 'star')) ) = [];
NumRead = full(dnavar.(fds{1}).numReadAlt + dnavar.(fds{1}).numReadRef);
for i = 2:length(fds)
    tmp = dnavar.(fds{i}).numReadAlt + dnavar.(fds{i}).numReadRef;
    idx = (isnan(NumRead) | NumRead == 0 ) & ~isnan(tmp) & tmp ~= 0;
    if nnz(idx)>0
        NumRead(idx) = tmp(idx);
    end
end
NumRead = sparse(NumRead);
for i = 1:length(fds)
    dnavar.(fds{i}) = rmfield(dnavar.(fds{i}), {'numReadAlt', 'numReadRef'});
end
dnavar.numRead = NumRead;

fds = unique(strrep(fds, '_star', '_bwa'));
covfilter = dnavar.(fds{1}).covfilter;
dnavar.(fds{1}) = rmfield(dnavar.(fds{1}), 'covfilter');
for i = 2:length(fds)
    covfilter = covfilter | dnavar.(fds{i}).covfilter;
    dnavar.(fds{i}) = rmfield(dnavar.(fds{i}), 'covfilter');
end
dnavar.bwa_covfilter = covfilter;
clear covfilter

fds = strrep(fds, '_bwa', '_star');
covfilter = dnavar.(fds{1}).covfilter;
dnavar.(fds{1}) = rmfield(dnavar.(fds{1}), 'covfilter');
for i = 2:length(fds)
    covfilter = covfilter | dnavar.(fds{i}).covfilter;
    dnavar.(fds{i}) = rmfield(dnavar.(fds{i}), 'covfilter');
end
dnavar.star_covfilter = covfilter;

for i = 1:length(fds)
    dnavar.(strrep(fds{i}, '_star', '')) = dnavar.(fds{i}).validAfterFilter;
end

fds = fieldnames(dnavar);
fds(cellfun(@isempty, strfind(fds, 'bwa')) & ...
    cellfun(@isempty, strfind(fds, 'star')) | ...
    ~cellfun(@isempty, strfind(fds, 'covfilter'))) = [];
dnavar = rmfield(dnavar, fds);

save(sprintf('data/%s', fn), 'dnavar');
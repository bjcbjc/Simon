function igvsnapshotscript(loc, figdir, scriptdir, scrname, append)

chrtag = {'X', 'Y', 'M'};
if ~append
    f = fopen([scriptdir, scrname], 'w');
else
    f = fopen([scriptdir, scrname], 'a');
end
fprintf(f, 'snapshotDirectory %s\n',figdir);
for i = 1:size(loc,1)
    if loc(i,1) > 22
        fprintf(f, 'goto chr%s:%d\n',chrtag{loc(i,1)-22}, loc(i,2));
    else
        fprintf(f, 'goto chr%d:%d\n',loc(i,1), loc(i,2));
    end
    fprintf(f, 'snapshot\n');    
end
fclose(f);


f = fopen('gencode.v16.annotation.gene.gtf','r');
t = textscan(f, '%s','delimiter','\n');
t = t{1};
fclose(f);

data = regexp(t, '^chr(\S+)\t(\S+)\tgene\t(\d+)\t(\d+)\t\S+\t([\+\-])\t\S+\tgene_id "(\S+)".*gene_type "(\S+)".*gene_status "(\S+)".*gene_name "(\S+)"', 'tokens', 'once');
data = vertcat(data{:});

GENCODE.id = data(:,6);
GENCODE.name = data(:,9);

chrmreplace = {'X', '23'; 'Y', '24'; 'M', '25'};
GENCODE.chrmstr = data(:,1);
for i = 1:size(chrmreplace, 1);
    GENCODE.chrmstr = strrep(GENCODE.chrmstr, chrmreplace{i,1}, chrmreplace{i,2});
end
GENCODE.chrm = str2double(GENCODE.chrmstr);
GENCODE.start = str2double(data(:,3));
GENCODE.end = str2double(data(:,4));
GENCODE.strand = data(:,5);

GENCODE.source = data(:,2);
GENCODE.type = data(:,7);
GENCODE.status = data(:,8);

save GENCODE GENCODE
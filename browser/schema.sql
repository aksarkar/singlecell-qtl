create table if not exists gene_info (gene string primary key, chr string, start int, end int, name string, strand string, source string);
create table if not exists log_mean (gene string, ind string, value real, primary key (gene, ind));
create table if not exists log_disp (gene string, ind string, value real, primary key (gene, ind));
create table if not exists umi (gene string, sample string, value int, primary key(gene, sample));
create table if not exists cell_to_ind (sample string primary key, ind string);
create table if not exists bulk (gene string, sample string, value real, primary key(gene, sample));

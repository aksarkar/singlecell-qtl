import numpy as np
import pandas as pd
import sqlite3

with sqlite3.connect('browser.db') as conn:
  conn.executescript("""
  create table if not exists gene_info (gene string primary key, chr string, start int, end int, name string, strand string, source string);
  create table if not exists cell_to_ind (sample string primary key, ind string);
  create table if not exists log_mean (gene string, ind string, value real, primary key (gene, ind));
  create table if not exists log_disp (gene string, ind string, value real, primary key (gene, ind));
  create table if not exists genotype (gene string, ind string, value real, primary key (gene, ind));
  create table if not exists umi (gene string, sample string, value int, primary key(gene, sample));
  create table if not exists bulk (gene string, sample string, value real, primary key(gene, sample));
  """)

gene_info = (pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-genes.txt.gz')
             .set_index('gene')
             .query('source == "H. sapiens"')
             .query('chr != "hsX"')
             .query('chr != "hsY"')
             .query('chr != "hsMT"'))
with sqlite3.connect('browser.db') as conn:
  gene_info.to_sql(name='gene_info', con=conn, if_exists='replace')

mean = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean.txt.gz', index_col='gene', sep=' ')
with sqlite3.connect('browser.db') as conn:
  (mean
   .reset_index()
   .melt(id_vars='gene', var_name='ind')
   .to_sql(name='log_mean', con=conn, index=False, if_exists='replace'))

disp = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-dispersion.txt.gz', index_col='gene', sep=' ')
with sqlite3.connect('browser.db') as conn:
  (disp
   .reset_index()
   .melt(id_vars='gene', var_name='ind')
   .to_sql(name='log_disp', con=conn, index=False, if_exists='replace'))

genotypes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/bulk-qtl-genotypes.txt.gz', index_col='gene', sep=' ')
with sqlite3.connect('browser.db') as conn:
  (genotypes
   .reset_index()
   .melt(id_vars='gene', var_name='ind').to_sql(name='genotype', con=conn, index=False, if_exists='replace'))

bulk = (pd.read_table('/project2/gilad/singlecell-qtl/bulk/counts_RNAseq_iPSC.txt', sep=' ', index_col='gene')
        .rename(index=lambda x: x.split('.')[0], columns=lambda x: 'NA{}'.format(x))
        .align(mean, axis=None, join='inner')[0])
bulk = np.log((bulk + 1) / bulk.sum(axis=0))
with sqlite3.connect('browser.db') as conn:
  (bulk
   .reset_index()
   .melt(id_vars='gene', var_name='ind')
   .to_sql(name='bulk', con=conn, index=False, if_exists='replace'))

annotations = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-annotation.txt')
keep_samples = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/quality-single-cells.txt', index_col=0, header=None)
annotations = annotations.loc[keep_samples.values.ravel()]
annotations['sample'] = annotations.apply(lambda x: '.'.join([x['chip_id'], str(x['experiment']), x['well']]), axis=1)
annotations['size'] = np.zeros(annotations.shape[0])
with sqlite3.connect('browser.db') as conn:
  for chunk in pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-counts.txt.gz', index_col=0, chunksize=100):
    chunk = chunk.align(mean, axis=None, join='inner')[0]
    annotations['size'] += chunk.sum(axis=0)
    chunk.reset_index().melt(id_vars='gene', var_name='ind').to_sql(name='umi', con=conn, index=False, if_exists='append')
  annotations[['sample', 'chip_id', 'size']].to_sql(name='annotation', con=conn, index=False, if_exists='replace')

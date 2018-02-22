"""Single cell eQTL browser"""

import numpy as np
import pandas as pd
import bokeh

def plot_umi(source):
  hist, edges = np.hist(source['umi'])
  p = bokeh.plotting.figure()
  p.quad(top=hist, bottom=0, left=edges[:-1], right=edges[1:])
  return p

# umi = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-counts.txt.gz', index_col=0)
# annotations = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-annotation.txt')
# keep_samples = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/quality-single-cells.txt', index_col=0, header=None)
# keep_genes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/genes-pass-filter.txt', index_col=0, header=None)
# umi = umi.loc[keep_genes.values.ravel(),keep_samples.values.ravel()]
# annotations = annotations.loc[keep_samples.values.ravel()]

# onehot = np.zeros((annotations.shape[0], len(categories)), dtype=np.float32)
# onehot[np.arange(onehot.shape[0]), annotations[key].apply(lambda x: categories.index(x))] = 1

# qtls = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean-qtls.txt.gz', compression='gzip', sep=' ')
# means = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean.txt.gz', compression='gzip', sep=' ')
# disps = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-dispersion.txt.gz', compression='gzip', sep=' ')

# genotypes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/bulk-qtl-genotypes.txt.gz', compression='gzip', sep=' ')

# source = bokeh.models.ColumnDataSource({'umi': umi.iloc[0]})
# bokeh.plotting.curdoc().add_root(bokeh.layouts.column(plot_umi(source)))

p = bokeh.plotting.figure()
p.line(x=np.arange(5), y=np.arange(5))
bokeh.plotting.curdoc().add_root(bokeh.layouts.column(p))

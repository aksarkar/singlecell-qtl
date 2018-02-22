"""Single cell eQTL browser

"""
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import numpy as np
import pandas as pd

gene_info = (pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/scqtl-genes.txt.gz')
             .set_index('gene')
             .query('source == "H. sapiens"')
             .query('chr != "hsX"')
             .query('chr != "hsY"')
             .query('chr != "hsMT"'))
means = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/zi2-mean.txt.gz', index_col='gene', sep=' ')
genotypes = pd.read_table('/home/aksarkar/projects/singlecell-qtl/data/bulk-qtl-genotypes.txt.gz', index_col='gene', sep=' ')
gene = genotypes.iloc[0].name

g, m = genotypes.loc[gene].align(means.loc[gene], join='inner')
sc_mean_by_geno = bokeh.plotting.figure(title=gene_info.loc[gene]['name'], width=600, height=400, tools=[])
sc_mean_by_geno.scatter(x=g.values, y=m.values, size=8)
sc_mean_by_geno.xaxis.axis_label = 'Centered genotype'
sc_mean_by_geno.yaxis.axis_label = 'Estimated single cell mean expression'

bulk = (pd.read_table('/project2/gilad/singlecell-qtl/bulk/counts_RNAseq_iPSC.txt', sep=' ', index_col='gene')
        .rename(index=lambda x: x.split('.')[0], columns=lambda x: 'NA{}'.format(x))
        .align(means, axis=None)[0])
bulk /= bulk.sum(axis=0)
bulk_mean_by_geno = bokeh.plotting.figure(title=gene_info.loc[gene]['name'], width=600, height=400, tools=[])
bulk_mean_by_geno.scatter(x=genotypes.loc[gene].values, y=bulk.loc[gene].values, size=8)
bulk_mean_by_geno.xaxis.axis_label = 'Centered genotype'
bulk_mean_by_geno.yaxis.axis_label = 'Bulk CPM'

layout = bokeh.layouts.layout([[sc_mean_by_geno, bulk_mean_by_geno]])

doc = bokeh.io.curdoc()
doc.title = 'scQTL browser'
doc.add_root(layout)

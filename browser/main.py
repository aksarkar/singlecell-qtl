"""Single cell eQTL browser

"""
import bokeh.io
import bokeh.layouts
import bokeh.models
import bokeh.plotting
import os.path
import numpy as np
import pandas as pd
import sqlite3

def update_gene(attr, old, new):
  selected = gene_data.selected['1d']['indices']
  if not selected:
    return
  with sqlite3.connect(os.path.join(os.path.dirname(__file__), 'browser.db')) as conn:
    gene = next(conn.execute('select gene from qtls where qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    ind_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      sql="""select genotype.ind, genotype.value as genotype, log_mean.value as mean, bulk.value as bulk 
          from log_mean, bulk, genotype 
          where log_mean.gene == bulk.gene and log_mean.ind == bulk.ind 
          and genotype.gene == bulk.gene and genotype.ind == bulk.ind and
          log_mean.gene == ?""",
      params=(gene,),
      con=conn))

def update_umi(attr, old, new):
  selected = ind_data.selected['1d']['indices']
  with sqlite3.connect(os.path.join(os.path.dirname(__file__), 'browser.db')) as conn:
    gene = next(conn.execute('select gene from qtls where qtls.gene == ?;', (gene_data.data['gene'][selected[0]],)))[0]
    if not selected:
      umi = pd.read_sql(
        """select umi.value, log_mean.value, log_disp.value, annotation.size 
        from annotation, log_mean, log_disp, umi 
        where annotation.chip_id == log_mean.ind and log_mean.ind == log_disp.ind 
        and umi.gene == ? and umi.sample == annotation.sample""",
        con=conn,
        params=(gene,))
    else:
      print(gene, ind_data.data['ind'][selected[0]])
      umi = pd.read_sql(
        """select umi.value, log_mean.value, log_disp.value, annotation.size 
        from annotation, log_mean, log_disp, umi 
        where annotation.chip_id == ? and annotation.chip_id == log_mean.ind 
        and log_mean.ind == log_disp.ind and umi.gene == ?
        and umi.sample == annotation.sample""",
        con=conn,
        params=(ind_data.data['ind'][selected[0]], gene,))
    counts, edges = np.histogram(umi['value'].values)
    umi_data.data = bokeh.models.ColumnDataSource.from_df(pd.DataFrame({'left': edges[:-1], 'right': edges[1:], 'count': counts}))

def init():
  with sqlite3.connect(os.path.join(os.path.dirname(__file__), 'browser.db')) as conn:
    gene_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      sql="""select gene_info.gene as gene, gene_info.name, qtls.id, qtls.p_bulk, qtls.beta_bulk, qtls.p_sc, qtls.beta_sc
      from gene_info, qtls 
      where gene_info.gene == qtls.gene
      order by p_bulk;""",
      con=conn))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['ind', 'genotype', 'mean', 'bulk']))

gene_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['gene', 'name', 'id', 'p_bulk', 'beta_bulk', 'p_sc', 'beta_sc']))
gene_data.on_change('selected', update_gene)

umi_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['left', 'right', 'count']))
dist_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['x', 'y']))

# These need to be module scope because bokeh.server looks there
qtls = bokeh.models.widgets.DataTable(
    source=gene_data,
    columns=[bokeh.models.widgets.TableColumn(field=x, title=x) for x in ['name', 'id', 'p_bulk', 'beta_bulk', 'p_sc', 'beta_sc']],
    width=1200,
    height=400)

sc_mean_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
sc_mean_by_geno.scatter(source=ind_data, x='genotype', y='mean', size=8)
sc_mean_by_geno.xaxis.axis_label = 'Centered genotype'
sc_mean_by_geno.yaxis.axis_label = 'Estimated single cell log mean expression'

bulk_mean_by_geno = bokeh.plotting.figure(width=400, height=400, tools=['tap'])
bulk_mean_by_geno.scatter(source=ind_data, x='genotype', y='bulk', size=8)
bulk_mean_by_geno.xaxis.axis_label = 'Centered genotype'
bulk_mean_by_geno.yaxis.axis_label = 'Bulk log CPM'

umi = bokeh.plotting.figure(width=400, height=400, tools=[])
umi.quad(source=umi_data, bottom=0, top='count', left='left', right='right')
umi.line(source=dist_data, x='x', y='y')
umi.xaxis.axis_label = 'Observed UMI'
umi.yaxis.axis_label = 'Number of cells'

layout = bokeh.layouts.layout([[qtls], [bulk_mean_by_geno, sc_mean_by_geno, umi]], sizing_mode='fixed')

doc = bokeh.io.curdoc()
doc.title = 'scQTL browser'
doc.add_root(layout)

init()

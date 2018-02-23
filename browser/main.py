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

def init(ind_data):
  with sqlite3.connect(os.path.join(os.path.dirname(__file__), 'browser.db')) as conn:
    gene = next(conn.execute('select gene from qtls order by p_beta limit 1;'))[0]
    ind_data.data = bokeh.models.ColumnDataSource.from_df(pd.read_sql(
      sql="""select genotype.value as genotype, log_mean.value as mean, bulk.value as bulk 
         from log_mean, bulk, genotype 
         where log_mean.gene == bulk.gene and log_mean.ind == bulk.ind 
         and genotype.gene == bulk.gene and genotype.ind == bulk.ind and
         log_mean.gene == ?""",
      params=(gene,),
      con=conn))

# These need to be separate because they have different dimension
ind_data = bokeh.models.ColumnDataSource(pd.DataFrame(columns=['genotype', 'mean', 'bulk']))
init(ind_data)

sc_mean_by_geno = bokeh.plotting.figure(width=600, height=400, tools=['tap'])
sc_mean_by_geno.scatter(source=ind_data, x='genotype', y='mean', size=8)
sc_mean_by_geno.xaxis.axis_label = 'Centered genotype'
sc_mean_by_geno.yaxis.axis_label = 'Estimated single cell log mean expression'

bulk_mean_by_geno = bokeh.plotting.figure(width=600, height=400, tools=['tap'])
bulk_mean_by_geno.scatter(source=ind_data, x='genotype', y='bulk', size=8)
bulk_mean_by_geno.xaxis.axis_label = 'Centered genotype'
bulk_mean_by_geno.yaxis.axis_label = 'Bulk log CPM'

layout = bokeh.layouts.layout([[sc_mean_by_geno, bulk_mean_by_geno]])

doc = bokeh.io.curdoc()
doc.title = 'scQTL browser'
doc.add_root(layout)


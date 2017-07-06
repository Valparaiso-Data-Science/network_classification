import pandas as pd
from bokeh.charts import Donut, show, output_file
from bokeh.palettes import d3
from bokeh.models import LabelSet

df = pd.read_csv('~/Downloads/network_classification/src/data/clean_data_with_new_chem.csv', index_col='Unnamed: 0')

d = Donut(df, label = 'Collection', color = d3['Category20'][16], plot_width = 600, plot_height=600)

output_file("donut.html")
show(d)


# Trying to create a heat map of correlations
#df_new = pd.DataFrame.copy(df)

#del df_new['Graph']
#del df_new['Collection']

#hm = figure()

#hm.rect(df_new, df_new, width=1, height=1)

#show(row(d, hm))
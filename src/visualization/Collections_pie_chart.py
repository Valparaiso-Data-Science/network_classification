import pandas as pd
from bokeh.charts import Donut, show, output_file
from bokeh.palettes import d3

df = pd.read_csv('~/Downloads/network_classification/src/data/clean_data_with_new_chem.csv', index_col='Unnamed: 0')

d = Donut(df, label = 'Collection', color = d3['Category20'][16], plot_width = 600, plot_height=600,
          hover_text='Collection')

output_file("donut.html")

show(d)


import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
from bokeh.charts import Bar, output_file, show


# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)
df_old = pd.read_csv('~/Downloads/network_classification/src/data/data_without_dimacs_or_bhoslib.csv', index_col=0)

# Plot pie chart with Collection column of df
#df['Collection'].value_counts().plot.pie(fontsize = 10.5, autopct='%1.0f%%')
#plt.ylabel('')
#
## Show plot
#plt.show()
my_colors = 'rgbkymc'
df_old['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)
plt.title('Counts of Collections Before Downsampling', fontsize=25)
#plt.show()


##plot bar chart with collection column of df
df['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)
plt.title('Counts of Collections After Downsampling', fontsize=25)

# Show plot
#plt.show()

# Plot in bokeh


dict_old = ColumnDataSource(dict(x=df_old['Collection'].value_counts(), y = df_old['Collection'].unique(),))
p=figure(title = 'Collections of graphs ', plot_width=1000)
p.vbar(bottom=0, top=df['Collection'].value_counts(), x=df['Collection'].unique(), width=100)
p.vbar(bottom=0, top=df_old['Collection'].value_counts(), x=df_old['Collection'].unique(), width=100)


import pandas as pd
import matplotlib.pyplot as plt
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
from bokeh.charts import Bar, output_file, show


# Read file
df = pd.read_csv('./../data/data_minmaxscale.csv', index_col=0)
df_old = pd.read_csv('./../data/data_without_dimacs_or_bhoslib.csv', index_col=0)

# Plot pie chart with Collection column of df
#df['Collection'].value_counts().plot.pie(fontsize = 10.5, autopct='%1.0f%%')
#plt.ylabel('')
#
## Show plot
#plt.show()
my_colors = 'rgbkymc'
p1 = plt.figure()
myfig=df_old['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)

#plt.title('Counts of Collections Before Downsampling', fontsize=25)
myfig.set_axisbelow(True)
myfig.yaxis.grid(color='gray', linestyle='dashed')


##plot bar chart with collection column of df
p2 = plt.figure()
myfig2 = df['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)
#plt.title('Counts of Collections After Downsampling', fontsize=25)
myfig2.set_axisbelow(True)
myfig2.yaxis.grid(color='gray', linestyle='dashed')

show(p1,p2)
#myfig2.show()

# Show plot
#plt.show()

# Plot in bokeh


#dict_old = ColumnDataSource(dict(x=df_old['Collection'].value_counts(), y = df_old['Collection'].unique(),))

#p=figure(title = 'Collections of graphs ', plot_width=1000)
#p.vbar(bottom=0, top=df['Collection'].value_counts(), x=df['Collection'].unique(), width=100)
#p.vbar(bottom=0, top=df_old['Collection'].value_counts(), x=df_old['Collection'].unique(), width=100)

#amounts = [119,114,62,48,36,22,19,16,13,11,9,9,8,6,2]
#names = ['Cheminformatics', 'Facebook Networks', 'Retweet Networks', 'Social Networks', 'Brain Networks',
#        'Web Graphs', 'Interaction Networks', 'Collaboration Networks', 'Recommendation Networks',
#         'Scientific Computing', 'Technological Networks', 'Biological Networks', 'Infrastructure Networks',
#         'Ecology Networks', 'Massive Network Data']

#source = ColumnDataSource(dict(x=amounts, y = names))

#p = Bar(df, values=amounts, label='Collection',group='Collection', title="Total MPG by CYL")

#show(p)
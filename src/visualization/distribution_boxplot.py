import pandas as pd
from bokeh.charts import BoxPlot, output_file, show
from bokeh.layouts import  column, row
from bokeh.models.widgets import Panel, Tabs

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_greater_than_20.csv', index_col=0)

collection = list(net['Collection'])
graph = list(net['Graph'])


#***************
# GENERATING PLOTS: BOKEH
#***************

pNodesC = BoxPlot(net, values='Nodes', label='Collection', color='Collection', title='Nodes', legend=None, plot_width=400, plot_height=400)
tNodesC = Panel(child=pNodesC, title='Nodes')
#edges by label:
pEdgesC = BoxPlot(net, values='Edges', label='Collection', color='Collection', title='Edges', legend=None, plot_width=400, plot_height=400)
tEdgesC = Panel(child=pEdgesC, title='Edges')
#density by label
pDensityC = BoxPlot(net, values='Density', label='Collection', color='Collection', title='Density', legend=None, plot_width=400, plot_height=400)
tDensityC = Panel(child=pDensityC, title='Density')
#max degree by label
pMaxDegC = BoxPlot(net, values='Maximum degree', label='Collection', color='Collection', title='Maximum Degree', legend=None, plot_width=400, plot_height=400)
tMaxDegC = Panel(child=pMaxDegC, title='Maximum Degree')
#min degree by label
pMinDegC = BoxPlot(net, values='Minimum degree', label='Collection', color='Collection', title='Minimum Degree', legend=None, plot_width=400, plot_height=400)
tMinDegC = Panel(child=pMinDegC, title='Minimum Degree')
#avg degree by label
pAvgDegC = BoxPlot(net, values='Average degree', label='Collection', color='Collection', title='Average Degree', legend=None, plot_width=400, plot_height=400)
tAvgDegC = Panel(child=pAvgDegC, title='Average Degree')
#assortativity by label
pAssC = BoxPlot(net, values='Assortativity', label='Collection', color='Collection', title='Assortativity', legend=None, plot_width=400, plot_height=400)
tAssC = Panel(child=pAssC, title='Assortativity')
#total triangles by label
pTotTriC = BoxPlot(net, values='Total triangles', label='Collection', color='Collection', title='Total Triangles', legend=None, plot_width=400, plot_height=400)
tTotTriC = Panel(child=pTotTriC, title='Total Triangles')
#average triangles by label
pAvgTriC = BoxPlot(net, values='Average triangles', label='Collection', color='Collection', title='Average Triangles', legend=None, plot_width=400, plot_height=400)
tAvgTriC = Panel(child=pAvgTriC, title='Average Triangles')
#maximum triangles by label
pMaxTriC = BoxPlot(net, values='Maximum triangles', label='Collection', color='Collection', title='Maximum Triangles', legend=None, plot_width=400, plot_height=400)
tMaxTriC = Panel(child=pMaxTriC, title='Maximum Triangles')
#average clustering coefficient by label
pAvgCCC = BoxPlot(net, values='Avg. clustering coef.', label='Collection', color='Collection', title='Average Clustering Coef.', legend=None, plot_width=400, plot_height=400)
tAvgCCC = Panel(child=pAvgCCC, title='Average Clustering Coef.')
#fraction of closed triangles by label
pFracCTC = BoxPlot(net, values='Frac. closed triangles', label='Collection', color='Collection', title='Fraction of Closed Triangles', legend=None, plot_width=400, plot_height=400)
tFracCTC = Panel(child=pFracCTC, title='Frac. Closed Triangles')
#maximum k-core by label
pMaxKC = BoxPlot(net, values='Maximum k-core', label='Collection', color='Collection', title='Maximum k-core', legend=None, plot_width=400, plot_height=400)
tMaxKC = Panel(child=pMaxKC, title='Max. k-Core')
#maximum clique by label
pMaxCliqueC = BoxPlot(net, values='Max. clique (lb)', label='Collection', color='Collection', title='Maximum Clique (lb)', legend=None, plot_width=400, plot_height=400)
tMaxCliqueC = Panel(child=pMaxCliqueC, title='Max. Clique (lb)')



#layout = Tabs(tabs=[tNodesC, tEdgesC, tDensityC, tMaxDegC, tMinDegC, tAvgDegC, tAssC, tTotTriC, tAvgTriC, tMaxTriC, tAvgCCC, tFracCTC, tMaxKC, tMaxCliqueC])
#changing font sizes of six graphs
pEdgesC.title.text_font_size = '20pt'
pAvgDegC.title.text_font_size = '20pt'
pAssC.title.text_font_size = '20pt'
pAvgCCC.title.text_font_size = '20pt'
pFracCTC.title.text_font_size = '20pt'
pMaxCliqueC.title.text_font_size = '20pt'

layout = column(row(pEdgesC, pAvgDegC, pAssC), row(pAvgCCC, pFracCTC, pMaxCliqueC))
output_file('dist_boxplot.html')
show(layout)
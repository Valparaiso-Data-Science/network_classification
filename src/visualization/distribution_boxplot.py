import pandas as pd
from bokeh.charts import BoxPlot, output_file, show
from bokeh.layouts import  column, row
from bokeh.models.widgets import Panel, Tabs
from bokeh.models.widgets.tables import DataTable

#net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_greater_than_20.csv', index_col=0)
net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv', index_col=0)
#net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)
net_tsne = pd.read_csv('~/PycharmProjects/network_classification/src/data/tsne_label_data.csv', index_col=0)

collection = list(net['Collection'])
graph = list(net['Graph'])

net['Label'] = list(net_tsne['Label'])
#***************
# GENERATING PLOTS: BOKEH
#***************

pNodesC = BoxPlot(net, values='Nodes', label='Label', color='Label', title='Nodes', legend=None, plot_width=400, plot_height=400)#,outliers=False)
tNodesC = Panel(child=pNodesC, title='Nodes')
#edges by label:
pEdgesC = BoxPlot(net, values='Edges', label='Label', color='Label', title='Edges', legend=None, plot_width=400, plot_height=400)#,outliers=False)
tEdgesC = Panel(child=pEdgesC, title='Edges')
#density by label
pDensityC = BoxPlot(net, values='Density', label='Label', color='Label', title='Density', legend=None, plot_width=400, plot_height=400)#,outliers=False)
tDensityC = Panel(child=pDensityC, title='Density')
#max degree by label
pMaxDegC = BoxPlot(net, values='Maximum degree', label='Label', color='Label', title='Maximum Degree', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tMaxDegC = Panel(child=pMaxDegC, title='Maximum Degree')
#min degree by label
pMinDegC = BoxPlot(net, values='Minimum degree', label='Label', color='Label', title='Minimum Degree', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tMinDegC = Panel(child=pMinDegC, title='Minimum Degree')
#avg degree by label
pAvgDegC = BoxPlot(net, values='Average degree', label='Label', color='Label', title='Average Degree', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tAvgDegC = Panel(child=pAvgDegC, title='Average Degree')
#assortativity by label
pAssC = BoxPlot(net, values='Assortativity', label='Label', color='Label', title='Assortativity', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tAssC = Panel(child=pAssC, title='Assortativity')
#total triangles by label
pTotTriC = BoxPlot(net, values='Total triangles', label='Label', color='Label', title='Total Triangles', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tTotTriC = Panel(child=pTotTriC, title='Total Triangles')
#average triangles by label
pAvgTriC = BoxPlot(net, values='Average triangles', label='Label', color='Label', title='Average Triangles', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tAvgTriC = Panel(child=pAvgTriC, title='Average Triangles')
#maximum triangles by label
pMaxTriC = BoxPlot(net, values='Maximum triangles', label='Label', color='Label', title='Maximum Triangles', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tMaxTriC = Panel(child=pMaxTriC, title='Maximum Triangles')
#average clustering coefficient by label
pAvgCCC = BoxPlot(net, values='Avg. clustering coef.', label='Label', color='Label', title='Average Clustering Coef.', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tAvgCCC = Panel(child=pAvgCCC, title='Average Clustering Coef.')
#fraction of closed triangles by label
pFracCTC = BoxPlot(net, values='Frac. closed triangles', label='Label', color='Label', title='Fraction of Closed Triangles', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tFracCTC = Panel(child=pFracCTC, title='Frac. Closed Triangles')
#maximum k-core by label
pMaxKC = BoxPlot(net, values='Maximum k-core', label='Label', color='Label', title='Maximum k-core', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tMaxKC = Panel(child=pMaxKC, title='Max. k-Core')
#maximum clique by label
pMaxCliqueC = BoxPlot(net, values='Max. clique (lb)', label='Label', color='Label', title='Maximum Clique (lb)', legend=None, plot_width=400, plot_height=400)#, outliers=False)
tMaxCliqueC = Panel(child=pMaxCliqueC, title='Max. Clique (lb)')

#**************
#TABLE TIME
#**************

#layout = Tabs(tabs=[tNodesC, tEdgesC, tDensityC, tMaxDegC, tMinDegC, tAvgDegC, tAssC, tTotTriC, tAvgTriC, tMaxTriC, tAvgCCC, tFracCTC, tMaxKC, tMaxCliqueC])
#changing font sizes of six graphs
pEdgesC.title.text_font_size = '20pt'
pAvgDegC.title.text_font_size = '20pt'
pAssC.title.text_font_size = '20pt'
pAvgCCC.title.text_font_size = '20pt'
pFracCTC.title.text_font_size = '18pt'
pMaxCliqueC.title.text_font_size = '20pt'

layout = column(row(pEdgesC, pAvgDegC, pAssC), row(pAvgCCC, pFracCTC, pMaxCliqueC))
output_file('boxplot_tsne_labels.html')
show(layout)
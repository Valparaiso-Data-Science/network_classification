import pandas as pd
from sklearn.cluster import KMeans
import matplotlib.pyplot as plt
from bokeh.charts import BoxPlot, output_file, show
from bokeh.models.widgets import Panel, Tabs

net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)

collection = list(net['Collection'])
graph = list(net['Graph'])

del net['Graph']
del net['Collection']
net_array = net.values
columns = net.columns
kmeans = KMeans(n_clusters=8, random_state=42)
labels = kmeans.fit_predict(net_array)

#***************
# TEST KMEANS
#***************
#print(labels)
net_ct = pd.DataFrame({'Labels':labels, 'Collection':collection})
ct = pd.crosstab( net_ct['Collection'], net_ct['Labels'])
print(ct)

#************
# SEPARATING LABELS
#************

net['Label'] = labels

c0 = net[net['Label']==0]
c1 = net[net['Label']==1]
c2 = net[net['Label']==2]
c3 = net[net['Label']==3]
c4 = net[net['Label']==4]
c5 = net[net['Label']==5]
c6 = net[net['Label']==6]
c7 = net[net['Label']==7]

#***************
# GENERATING PLOTS: BOKEH
#***************

#nodes by label:
pNodes = BoxPlot(net, values='Nodes', label='Label', title='Nodes by Label', color='Label')
tNodes = Panel(child=pNodes, title='Nodes')
#edges by label:
pEdges = BoxPlot(net, values='Edges', label='Label', title='Edges by Label', color='Label')
tEdges = Panel(child=pEdges, title='Edges')
#density by label
pDensity = BoxPlot(net, values='Density', label='Label', title='Density by Label', color='Label')
tDensity = Panel(child=pDensity, title='Density')
#max degree by label
pMaxDeg = BoxPlot(net, values='Maximum degree', label='Label', title='Maximum Degree by Label', color='Label')
tMaxDeg = Panel(child=pMaxDeg, title='Maximum Degree')
#min degree by label
pMinDeg = BoxPlot(net, values='Minimum degree', label='Label', title='Minimum Degree by Label', color='Label')
tMinDeg = Panel(child=pMinDeg, title='Minimum Degree')
#avg degree by label
pAvgDeg = BoxPlot(net, values='Average degree', label='Label', title='Average Degree by Label', color='Label')
tAvgDeg = Panel(child=pAvgDeg, title='Average Degree')
#assortativity by label
pAss = BoxPlot(net, values='Assortativity', label='Label', title='Assortativity by Label', color='Label')
tAss = Panel(child=pAss, title='Assortativity')
#total triangles by label
pTotTri = BoxPlot(net, values='Total triangles', label='Label', title='Total Triangles by Label', color='Label')
tTotTri = Panel(child=pTotTri, title='Total Triangles')
#average triangles by label
pAvgTri = BoxPlot(net, values='Average triangles', label='Label', title='Average Triangles by Label', color='Label')
tAvgTri = Panel(child=pAvgTri, title='Average Triangles')
#maximum triangles by label
pMaxTri = BoxPlot(net, values='Maximum triangles', label='Label', title='Maximum Triangles by Label', color='Label')
tMaxTri = Panel(child=pMaxTri, title='Maximum Triangles')
#average clustering coefficient by label
pAvgCC = BoxPlot(net, values='Avg. clustering coef.', label='Label', title='Average Clustering Coefficient by Label', color='Label')
tAvgCC = Panel(child=pAvgCC, title='Average Clustering Coef.')
#fraction of closed triangles by label
pFracCT = BoxPlot(net, values='Frac. closed triangles', label='Label', title='Fraction of Closed Triangles by Label', color='Label')
tFracCT = Panel(child=pFracCT, title='Frac. Closed Triangles')
#maximum k-core by label
pMaxK = BoxPlot(net, values='Maximum k-core', label='Label', title='Maximum k-core by Label', color='Label')
tMaxK = Panel(child=pMaxK, title='Max. k-Core')
#maximum clique by label
pMaxClique = BoxPlot(net, values='Max. clique (lb)', label='Label', title='Lower Bound of Maximum Clique by Label', color='Label')
tMaxClique = Panel(child=pMaxClique, title='Max. Clique (lb)')

layout = Tabs(tabs=[tNodes, tEdges, tDensity, tMaxDeg, tMinDeg, tAvgDeg, tAss, tTotTri, tAvgTri, tMaxTri, tAvgCC, tFracCT, tMaxK, tMaxClique])

output_file('boxplot.html')
show(layout)
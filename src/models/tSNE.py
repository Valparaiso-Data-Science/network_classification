#this program makes a tSNE visualization of data
#edited to make use of pca data
#currently using MinMaxScaler and learning rate of 1000.

import pandas as pd
from sklearn.manifold import TSNE
import numpy as np
from sklearn.preprocessing import MinMaxScaler, Imputer
from sklearn.pipeline import make_pipeline
import matplotlib.pyplot as plt
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import scipy


net = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)

#data for use with pca file:
#pca = pd.read_csv('~/PycharmProjects/network_classification/src/data/pca_for_tsne.csv', index_col=0)

#i need a list of the collections as numbers-useful in coloring tsne
graph_names = list(net['Graph'])
collection_names = list(net['Collection'])
net['Collection'] = net['Collection'].astype('category')
net['Collection'] = net['Collection'].cat.codes

collection_num = net['Collection']

del net['Collection']
del net['Graph']

#pca_array = pca.values
net_array = net.values


tsne = TSNE(metric=scipy.spatial.distance.canberra)
#random_state=987 for snake


tsne_features = tsne.fit_transform(net_array)
xs = tsne_features[:,0]
ys = tsne_features[:,1]

#plt.scatter(xs, ys, c=collection_num)
#plt.show()

#df = pd.DataFrame({"x" : xs, "y" : ys, "Category Number" : net['Collection'], "Category Name":collection_names, "Graph Name": graph_names})
#df.to_csv('~/PycharmProjects/network_classification/src/data/pca_tsne_data.csv')

data = {'x':xs, 'y':ys, 'Collection':collection_names, 'Graph':graph_names}
df = pd.DataFrame(data)

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Collection}")]

# Creating the scatter plot of the x and y coordinates
p=figure(title = 't-SNE ', plot_width=1000)

# Color the plot by collection
category = df['Collection']
all_categories = category.unique().tolist()
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Collection'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][16][i], size = 8, legend = graph)

# Creating scatter plot of centroids
#p.square_cross(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('canberra_tsne_plot.html')
show(p)
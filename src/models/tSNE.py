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


#net = pd.read_csv('~/PycharmProjects/network_classification/src/data/clean_data_with_new_chem.csv')

#data for use with pca file:
pca = pd.read_csv('~/PycharmProjects/network_classification/src/data/pca_for_tsne.csv')
#i need a list of the collections as numbers-useful in coloring tsne
graph_names = list(pca['Graph'])
collection_names = list(pca['Collection'])
pca['Collection'] = pca['Collection'].astype('category')
pca['Collection'] = pca['Collection'].cat.codes

del pca['Graph']

#graph_categories = []
#all_categories = net['Collection'].values

#all_cat_names = list(net['Collection'].values)


#creates a list of the existing categories, no repeats
#for category in all_categories:
#    if not category in graph_categories:
#        graph_categories.append( category )
#re-assigns each collection name in all_categories to an integer, corresponding with its position in graph_categories
#for i in range(len( all_categories )):
#    all_categories[i] = graph_categories.index( all_categories[i] )

#all_graph_names = list(net['Graph'])

#del net['Graph']
#del net['Collection']


#net_array = net.values
pca_array = pca.values

#commented out scaler for pca stuff
#min_max_scale = MinMaxScaler()
#imp = Imputer(missing_values='NaN', strategy='most_frequent', axis=0)

tsne = TSNE()

#pipeline = make_pipeline(min_max_scale, tsne)
#tsne_features = pipeline.fit_transform(net_array)

tsne_features = tsne.fit_transform(pca_array)
xs = tsne_features[:,0]
ys = tsne_features[:,1]

plt.scatter(xs, ys, c=pca['Collection'])
plt.show()

#df = pd.DataFrame({"x" : xs, "y" : ys, "Category Number" : pca['Collection'], "Category Name":collection_names, "Graph Name": graph_names})
#df.to_csv('~/PycharmProjects/network_classification/src/data/pca_tsne_data.csv')

data = {'x':xs, 'y':ys, 'Collection':collection_names, 'Graph':graph_names}
df = pd.DataFrame(data)

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Collection}")]

# Creating the scatter plot of the x and y coordinates
p=figure(title = 'PCA t-SNE ', plot_width=1000)

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
output_file('pca_tsne_plot.html')
show(p)
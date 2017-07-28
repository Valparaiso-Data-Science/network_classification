# This program creates plots using data ran on tsne and
# plots it with bokeh using kmeans to get the clusters

import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import matplotlib.pyplot as plt



# Read file
tsne_data = pd.read_csv('~/Downloads/network_classification/src/data/new_tsne_data.csv', index_col=0)

tsne_data = tsne_data[tsne_data['Category Name'] != 'Temporal Reachability']

# Make a copy of the data
tsne_new = pd.DataFrame.copy(tsne_data)


# Delete categorical columns
del tsne_new['Category Name']
del tsne_new['Graph Name']
del tsne_new['Category Number']


# Create array of data (only numerical columns)
tsne_array = tsne_new.values



#**************************
# Getting inertia graph
#**************************
def inertia_plot(data_array):
    ks = range(1, 20)
    inertias = []

    for k in ks:
        # Create a KMeans instance with k clusters: model
        model = KMeans(n_clusters=k)

        # Fit model to samples
        # Change to the data we need to check
        model.fit(data_array)

        # Append the inertia to the list of inertias
        inertias.append(model.inertia_)

    # Plot ks vs inertias
    plt.plot(ks, inertias, '-o')
    plt.title('Inertia vs. k')
    plt.xlabel('number of clusters, k')
    plt.ylabel('inertia')
    plt.xticks(ks)
    plt.show()



#**************************
# Creating plot of kmeans
# using tsne data in bokeh
#**************************
# Create KMeans with 8 clusters and fit data to model
kmeans = KMeans(n_clusters = 8)
kmeans.fit_transform(tsne_array)
labels = kmeans.predict(tsne_array)
centroids = kmeans.cluster_centers_

#**************************
#write out labels for use in boxplots/table
#**************************

#tsne_data['Label'] = labels
#tsne_data.to_csv('~/PycharmProjects/network_classification/src/data/tsne_label_data.csv')

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x = centroids[:,0]
centroids_y = centroids[:,1]

# Create cross tabulation of tsne data and print
df1 = pd.DataFrame({'labels':labels, 'Collection':tsne_data['Category Name']})
print("Crosstab for t-SNE data:\n")
ct = pd.crosstab(df1['Collection'], df1['labels'])
print(ct)

# Get category and graph names for new dataframe
category = tsne_data['Category Name']
names = tsne_data['Graph Name']

# Get all categories without repetitions
all_categories = category.unique().tolist()

# Assign the columns of tsne_array
xs = tsne_array[:,0]
ys = tsne_array[:, 1]
data = {'x': xs, 'y': ys, 'Category Name' : category, 'Graph': names}

# Create new pandas dataframe
df=pd.DataFrame(data)

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

# Creating the figure for the scatter plot
p=figure(title = 't-Distributed Stochastic Neighbor Embedding', plot_width=1000)
p.title.text_font_size = '25pt'

# Create scatter points and color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Category Name'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][17][i], size = 8, legend = graph)

# Creating scatter points of centroids
p.square(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"
#p.legend.label_text_font_size = "16pt"
#p.legend.background_fill_alpha = 0

# Save file and show plot
output_file('kmeans_centroids_plot.html')
show(p)


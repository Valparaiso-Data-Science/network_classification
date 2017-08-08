# This program runs k-means on scaled data and
# retrieves centroids. It then runs tsne on this new
# data and creates a plot of the clusters


import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
from sklearn.manifold import TSNE
import scipy


# Read file
raw = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

# Make a copy of the file
raw_new = pd.DataFrame.copy(raw)

# Delete categorical columns
del raw_new['Graph']
del raw_new['Collection']

# Create array of data frame values
raw_array = raw_new.values


# Create KMeans with 8 clusters and fit data to model
def kmeans_function(array):
    kmeans = KMeans(n_clusters = 8)
    kmeans.fit_transform(array)
    labels = kmeans.predict(array)
    centroids = kmeans.cluster_centers_
    #I DONT KNOW HOW RETURNS WORK IN PYTHON. CHECK THIS FOR ALL FUNCTIONS WITH RETURNS
    return labels, centroids

kmeans_function(raw_array)

# Create cross tabulation of raw data and print
def crosstab(labels, file):
    df = pd.DataFrame({'Labels':labels, 'Collection':file['Collection']})
    print("\nCrosstab for raw data:\n")
    ct = pd.crosstab(df['Collection'], df['Labels'])
    print(ct)

crosstab(labels, raw)

#df2 = pd.DataFrame({'labels':labels, 'Collection':file['Collection']})
#print("\nCrosstab for raw data:\n")
#ct = pd.crosstab(df2['Collection'], df2['labels'])
#print(ct)
#

# Get names of the numerical columns
column = raw_new.keys()

# Create new data frame with the centroids data from the kmeans in raw data
new_cent = pd.DataFrame(centroids, columns = column)
new_cent['Collection'] = 'Centroid'
new_cent['Graph'] = [0,1,2,3,4,5,6,7]

# Rearranging the order of the columns
new_cent = new_cent[['Graph', 'Collection', 'Nodes', 'Edges', 'Density', 'Maximum degree',
       'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles',
       'Average triangles', 'Maximum triangles', 'Avg. clustering coef.',
       'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)']]

# Concatenate raw data frame with centroids data frame
raw_with_cent = pd.concat([raw, new_cent])

# Create a copy of this data frame
raw_with_cent_copy = pd.DataFrame.copy(raw_with_cent)

# Delete categorical columns
del raw_with_cent_copy['Graph']
del raw_with_cent_copy['Collection']


# Get all the category names and graph names
category_new = raw_with_cent['Collection']
names_new = raw_with_cent['Graph']
all_new_categories = category_new.unique().tolist()


# Run tsne on the new data frame
def tsne_function(file):
    tsne = TSNE()
    tsne_features = tsne.fit_transform(file.values)
    xs_new = tsne_features[:, 0]
    ys_new = tsne_features[:, 1]
    data_new = {'x': xs_new, 'y': ys_new, 'Category Name': category_new, 'Graph': names_new}
    df_new = pd.DataFrame(data_new)

    return df_new

#tsne = TSNE()
#tsne_features = tsne.fit_transform(raw_with_cent_copy.values)
#xs_new = tsne_features[:, 0]
#ys_new = tsne_features[:, 1]
#
#data_new = {'x': xs_new, 'y': ys_new, 'Category Name' : category_new, 'Graph': names_new}
#
## Create new pandas dataframe
#df_new = pd.DataFrame(data_new)
#
#df_new.to_csv('~/Downloads/network_classification/src/models/tsne_14d_data.csv')

# Create new values that will use for label of centroids
values = np.array([-1, -1, -1, -1, -1, -1, -1, -1])
new_vals = np.append(np.array(labels), values)
new_vals = pd.DataFrame(new_vals)

# Add new column of labels to data frame
df_new['Label'] = new_vals


# Create hover tool for plot
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

# Creating the figure for the scatter plot
p=figure(title = 't-SNE in all dimensions', plot_width=1000)

all_new_labels = df_new['Label'].unique().tolist()

# Create scatter points and color the plot by collection
j = 0 # To use in loop
for i, graph in enumerate(all_new_categories):
    source = ColumnDataSource(df_new[df_new['Category Name'] == graph])
    if graph != 'Centroid':
        p.circle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
    # Plot the centroids
    else:
        p.square(x = 'x', y = 'y', source = source, color = 'black', size = 12, legend = graph)



# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('tsne_14d_plot.html')
show(p)

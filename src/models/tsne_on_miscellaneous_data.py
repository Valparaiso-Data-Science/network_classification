
import pandas as pd
from sklearn.manifold import TSNE
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline
import scipy

misc = pd.read_csv('~/Downloads/network_classification/src/data/miscellaneous_networks.csv', index_col=0)
del misc['Chromatic Number']
misc = misc.dropna()
scaled = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

df = pd.concat([scaled, misc])

df_new = pd.DataFrame.copy(df)

del df_new['Graph']
del df_new['Collection']
del df_new['Collection Hypothesis']

df_array = df_new.values


# Get names of the columns
column = df_new.keys()


# Get category and graph names for new dataframe
category = df['Collection']
names = df['Graph']

# Get all categories without repetitions
all_categories = category.unique().tolist()

# Run tsne on the new data frame
tsne = TSNE(metric=scipy.spatial.distance.canberra)
scaler = MinMaxScaler()
pipeline = make_pipeline(scaler, tsne)

tsne_features = pipeline.fit_transform(df_array)
xs = tsne_features[:, 0]
ys = tsne_features[:, 1]

data = {'x': xs, 'y': ys, 'Category Name' : category, 'Graph': names}

# Create new pandas dataframe
df_data=pd.DataFrame(data)

df_data_copy = pd.DataFrame.copy(df_data)

del df_data_copy['Category Name']
del df_data_copy['Graph']

#**************************
# Creating plot of kmeans
# using tsne data in bokeh
#**************************
# Create KMeans with 8 clusters and fit data to model
kmeans = KMeans(n_clusters = 8)
kmeans.fit_transform(df_data_copy)
labels = kmeans.predict(df_data_copy)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x = centroids[:,0]
centroids_y = centroids[:,1]

# Create cross tabulation of tsne data and print
df1 = pd.DataFrame({'labels':labels, 'Collection':df['Collection']})
print("Crosstab for t-SNE data:\n")
ct = pd.crosstab(df1['Collection'], df1['labels'])
print(ct)

# Add labels to data frame
df_data['Label'] = labels
all_labels = df_data['Label'].unique().tolist()

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}"), ("Cluster", "@Label")]

# Creating the figure for the scatter plot
p = figure(title = 't-SNE ', plot_width=1000)

# Create scatter points and color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df_data[df_data['Category Name'] == graph])
    p.circle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)

# Creating scatter points of centroids
p.square(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('kmeans_centroids_plot.html')
show(p)
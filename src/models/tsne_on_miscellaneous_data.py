
import pandas as pd
from sklearn.manifold import TSNE
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
from sklearn.cluster import KMeans
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline

misc = pd.read_csv('~/Downloads/network_classification/src/data/miscellaneous_networks.csv', index_col=0)
misc = misc.dropna()
del misc['Chromatic Number']
scaled = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

df = pd.concat([scaled, misc])

df_new = pd.DataFrame.copy(df)

del df_new['Collection']
del df_new['Graph']



df_array = df_new.values


# Get names of the columns
column = df_new.keys()

# Create new data frame with the centroids data from the kmeans in raw data
#new_cent = pd.DataFrame(centroids, columns = column)
#new_cent['Collection'] = 'Centroid'
#new_cent['Graph'] = [0,1,2,3,4,5,6,7]

# Rearranging the order of the columns
#new_cent = new_cent[['Graph', 'Collection', 'Nodes', 'Edges', 'Density', 'Maximum degree',
#       'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles',
#       'Average triangles', 'Maximum triangles', 'Avg. clustering coef.',
#       'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)']]

# Concatenate raw data frame with centroids data frame
#raw_with_cent = pd.concat([raw, new_cent])

# Create a copy of this data frame
#raw_with_cent_copy = pd.DataFrame.copy(raw_with_cent)

# Delete categorical columns
#del raw_with_cent_copy['Graph']
#del raw_with_cent_copy['Collection']

# Delete Nodes and Edges (Will this make a change?)
#del raw_with_cent_copy['Nodes']
#del raw_with_cent_copy['Edges']

# Get all the category names and graph names
#category_new = raw_with_cent['Collection']
#names_new = raw_with_cent['Graph']
#all_new_categories = category_new.unique().tolist()


# Get category and graph names for new dataframe
category = df['Collection']
names = df['Graph']

# Get all categories without repetitions
all_categories = category.unique().tolist()

# Run tsne on the new data frame
tsne = TSNE(random_state=42)
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
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

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
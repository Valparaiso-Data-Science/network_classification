import pandas as pd
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import matplotlib.pyplot as plt

# Read file
tsne_data = pd.read_csv('~/Downloads/network_classification/src/data/new_tsne_data.csv', index_col=0)
scaled = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

# Make a copy of the data
tsne_new = pd.DataFrame.copy(tsne_data)
scaled_new = pd.DataFrame.copy(scaled)

# Delete categorical columns
del tsne_new['Category Name']
del tsne_new['Graph Name']
del tsne_new['Category Number']
del scaled_new['Graph']
del scaled_new['Collection']
del scaled_new['Nodes']
del scaled_new['Edges']

# Create array of data (only numerical columns)
tsne_array = tsne_new.values
scaled_array = scaled_new.values

#**************************
# Getting inertia graph
#**************************
ks = range(1, 20)
inertias = []

for k in ks:
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k)

    # Fit model to samples
    model.fit(tsne_array)

    # Append the inertia to the list of inertias
    inertias.append(model.inertia_)

# Plot ks vs inertias
plt.plot(ks, inertias, '-o')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(ks)
#plt.show()

#**************************
# Creating plot of kmeans
# using tsne data in bokeh
#**************************
# Create KMeans with 8 clusters and fit data to model
kmeans = KMeans(n_clusters = 8)
kmeans.fit_transform(tsne_array)
labels = kmeans.predict(tsne_array)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x = centroids[:,0]
centroids_y = centroids[:,1]

# Create cross tabulation and print
df1 = pd.DataFrame({'labels':labels, 'Collection':tsne_data['Category Name']})
print("Crosstab for t-SNE data:\n")
ct = pd.crosstab(df1['Collection'], df1['labels'])
print(ct)

# Create second KMeans with scaled data
kmeans.fit_transform(scaled_array)
labels = kmeans.predict(scaled_array)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x_scaled = centroids[:,0]
centroids_y_scaled = centroids[:,1]

df2 = pd.DataFrame({'labels':labels, 'Collection':scaled['Collection']})
print("\nCrosstab for scaled data:\n")
ct = pd.crosstab(df2['Collection'], df2['labels'])
print(ct)

# Assign the columns of tsne_array
xs = tsne_array[:,0]
ys = tsne_array[:, 1]

# Get category and graph names for new dataframe
category = tsne_data['Category Name']
names = tsne_data['Graph Name']

# Get all categories without repetitions
all_categories = category.unique().tolist()

data = {'x': xs, 'y': ys, 'Category Name' : category, 'Graph': names}

# Create new pandas dataframe
df=pd.DataFrame(data)

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

# Creating the figure for the scatter plot
p=figure(title = 't-SNE ', plot_width=1000)

# Create scatter points and color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Category Name'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][16][i], size = 8, legend = graph)

# Creating scatter points of centroids
p.square(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('kmeans_centroids_plot.html')
show(p)




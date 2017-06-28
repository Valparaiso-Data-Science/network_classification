import pandas as pd
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import matplotlib.pyplot as plt

# Read file
net = pd.read_csv('~/PycharmProjects/network_classification/src/data//new_tsne_data.csv')

# Take out Remporal Reachability and Dynamic Networks
net = net[net['Category Name'] != 'Temporal Reachability']
net = net[net['Category Name'] != 'Dynamic Networks']

# Make a copy of the data
net_new = pd.DataFrame.copy(net)

# Delete categorical columns
del net_new['Category Name']
del net_new['Graph Name']
del net_new['Category Number']
del net_new['Unnamed: 0']

# Create array of new data (only numerical columns)
net_array = net_new.values

#**************************
# Getting inertia graph
#**************************
ks = range(1, 20)
inertias = []

for k in ks:
    # Create a KMeans instance with k clusters: model
    model = KMeans(n_clusters=k)

    # Fit model to samples
    model.fit(net_array)

    # Append the inertia to the list of inertias
    inertias.append(model.inertia_)

# Plot ks vs inertias
plt.plot(ks, inertias, '-o')
plt.xlabel('number of clusters, k')
plt.ylabel('inertia')
plt.xticks(ks)
plt.show()

#**************************
# Creating plot of kmeans
# using tsne data with bokeh
#**************************
# Create KMeans with 10 clusters and fit data to model
kmeans = KMeans(n_clusters = 8)
kmeans.fit_transform(net_array)
labels = kmeans.predict(net_array)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x = centroids[:,0]
centroids_y = centroids[:,1]

# Create cross tabulation and print
df1 = pd.DataFrame({'labels':labels, 'Collection':net['Category Name']})
ct = pd.crosstab(df1['Collection'], df1['labels'])
print(ct)

# Assign the columns of net_array
xs = net_array[:,0]
ys = net_array[:, 1]

# Get category and graph names for new dataframe
category = net['Category Name']
names = net['Graph Name']

# Get all categories without repetitions
all_categories = category.unique().tolist()

data = {'x': xs, 'y': ys, 'Category Name' : category, 'Graph': names}

# Create new pandas dataframe
df=pd.DataFrame(data)

# Create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

# Creating the scatter plot of the x and y coordinates
p=figure(title = 't-SNE ', plot_width=1000)

# Color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Category Name'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][16][i], size = 8, legend = graph)

# Creating scatter plot of centroids
p.square_cross(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('kmeans_centroids_plot.html')
show(p)




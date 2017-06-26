import pandas as pd
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import matplotlib.pyplot as plt

#read file
net = pd.read_csv('/Users/adrianaortiz/Downloads/new_tsne_data.csv')

#take out cheminformatics
#net = net[net['Category Name']!='Cheminformatics']
net = net[net['Category Name'] != 'Temporal Reachability']
net = net[net['Category Name'] != 'Dynamic Networks']

#make a copy of the data
net_new = pd.DataFrame.copy(net)

#delete categorical columns
del net_new['Category Name']
del net_new['Graph Name']
del net_new['Category Number']
del net_new['Unnamed: 0']

#create array of new data (only numerical columns)
net_array = net_new.values

#**************************
#getting inertia graph
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
# creating plot of kmeans
# using tsne data with bokeh
#**************************
#create KMeans with 10 clusters and fit data to model
kmeans = KMeans(n_clusters = 8)
kmeans.fit_transform(net_array)
labels = kmeans.predict(net_array)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x = centroids[:,0]
centroids_y = centroids[:,1]

# create cross tabulation and print
df1 = pd.DataFrame({'labels':labels, 'Collection':net['Category Name']})
ct = pd.crosstab(df1['Collection'], df1['labels'])
print(ct)

# assign the columns of net_array
xs = net_array[:,0]
ys = net_array[:, 1]
category = net['Category Name']
names = net['Graph Name']

# get all categories without repetitions
all_categories = category.unique().tolist()

data = {'x': xs, 'y': ys, 'Category Name' : category, 'Graph': names}

# create new dataframe
df=pd.DataFrame(data)

# create hover tool
hover = HoverTool()
hover.tooltips = [("Graph", "@Graph"),("Category", "@{Category Name}")]

# creating the scatter plot of the x and y coordinates
p=figure(title = 't-SNE ', plot_width=1000)

# color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Category Name'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][16][i], size=8, legend = graph)

# creating scatter plot of centroids
p.add_tools(hover)
p.square_cross(centroids_x, centroids_y, color='black', size=12, legend = 'Centroid')
p.legend.click_policy="hide"

# save file and show plot
output_file('kmeans_centroids_plot.html')
show(p)




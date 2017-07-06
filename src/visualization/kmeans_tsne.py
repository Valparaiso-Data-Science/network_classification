import pandas as pd
import numpy as np
from sklearn.cluster import KMeans
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE
from sklearn.preprocessing import MinMaxScaler
from sklearn.pipeline import make_pipeline

# Read file

tsne_data = pd.read_csv('~/Downloads/network_classification/src/data/new_tsne_data.csv', index_col=0)
raw = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

# Make a copy of the data
tsne_new = pd.DataFrame.copy(tsne_data)
raw_new = pd.DataFrame.copy(raw)

# Delete categorical columns
del tsne_new['Category Name']
del tsne_new['Graph Name']
del tsne_new['Category Number']
del raw_new['Graph']
del raw_new['Collection']
#del raw_new['Nodes']
#del raw_new['Edges']

# Create array of data (only numerical columns)
tsne_array = tsne_new.values
raw_array = raw_new.values


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
# Uncomment to show inertia graph
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
p=figure(title = 't-SNE ', plot_width=1000)

# Create scatter points and color the plot by collection
for i, graph in enumerate(all_categories):
    source = ColumnDataSource(df[df['Category Name'] == graph])
    p.circle(x='x', y='y', source = source, color = d3['Category20'][16][i], size = 8, legend = graph)

# Creating scatter points of centroids
p.square(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')
#p.square(centroids_x_raw, centroids_y_raw, color = 'blue', size = 12, legend = 'Centroid for raw data')

# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
#output_file('kmeans_centroids_plot.html')
#show(p)



#*******************************************
# Running the same thing in a different file
#*******************************************

# Create second KMeans with raw data
kmeans.fit_transform(raw_array)
labels = kmeans.predict(raw_array)
centroids = kmeans.cluster_centers_

# Assign the columns of centroids: centroids_x, centroids_y
centroids_x_raw = centroids[:,0]
centroids_y_raw = centroids[:,1]

# Create cross tabulation of raw data and print
df2 = pd.DataFrame({'labels':labels, 'Collection':raw['Collection']})
print("\nCrosstab for raw data:\n")
ct = pd.crosstab(df2['Collection'], df2['labels'])
print(ct)

# Get names of the columns
column = raw_new.keys()

# Create new data frame with the centroids data from the kmeans in raw data
new_cent = pd.DataFrame(centroids, columns = column)
new_cent['Collection'] = 'Centroid'
new_cent['Graph'] = [1,2,3,4,5,6,7,8]

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

# Delete Nodes and Edges (Will this make a change?)
#del raw_with_cent_copy['Nodes']
#del raw_with_cent_copy['Edges']

# Get all the category names and graph names
category_new = raw_with_cent['Collection']
names_new = raw_with_cent['Graph']
all_new_categories = category_new.unique().tolist()

# Run tsne on the new data frame
tsne = TSNE(random_state=42)
tsne_features = tsne.fit_transform(raw_with_cent_copy.values)
xs_new = tsne_features[:, 0]
ys_new = tsne_features[:, 1]

data_new = {'x': xs_new, 'y': ys_new, 'Category Name' : category_new, 'Graph': names_new}

# Create new pandas dataframe
df_new = pd.DataFrame(data_new)

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
p=figure(title = 'Scaled Data ', plot_width=1000)

j=0
# Create scatter points and color the plot by collection
for i, graph in enumerate(all_new_categories):
    source = ColumnDataSource(df_new[df_new['Category Name'] == graph])
    if graph != 'Centroid':
        if df_new.loc[j, 'Label'] == 0:
            p.circle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 1:
            p.triangle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 2:
            p.diamond(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 3:
            p.asterisk(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 4:
            p.cross(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 5:
            p.x(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        elif df_new.loc[j, 'Label'] == 6:
            p.inverted_triangle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        else:
            p.square(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)
        j += 1
    else:
        p.square(x = 'x', y = 'y', source = source, color = 'black', size = 12, legend = graph)


# Add tools and interactive legend
p.add_tools(hover)
p.legend.location = "top_left"
p.legend.click_policy="hide"

# Save file and show plot
output_file('raw_data_kmeans_centroids_plot.html')
show(p)

import pandas as pd
import numpy as np
from sklearn.preprocessing import Normalizer, Imputer
from sklearn.pipeline import make_pipeline
from sklearn.cluster import KMeans
from sklearn.manifold import TSNE
from bokeh.plotting import figure, output_file, show
from bokeh.models import HoverTool, ColumnDataSource
from bokeh.palettes import d3
import scipy

# Read file
tsne_data = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)


# Make a copy of the data
tsne_new = pd.DataFrame.copy(tsne_data)


# Delete categorical columns
del tsne_new['Collection']
del tsne_new['Graph']
#del tsne_new['Category Number']


# Create array of data (only numerical columns)
tsne_array = tsne_new.values

# Get category and graph names for new dataframe
category = tsne_data['Collection']
names = tsne_data['Graph']

# Get all categories without repetitions
all_categories = category.unique().tolist()

lrates = [100, 200, 300, 400, 500, 600, 700, 800, 900, 1000]
for l in lrates:
    tsne = TSNE(learning_rate=l, metric=scipy.spatial.distance.canberra)



    tsne_features = tsne.fit_transform(tsne_array)
    xs = tsne_features[:,0]
    ys = tsne_features[:,1]



    # **************************
    # Creating plot of kmeans
    # using tsne data in bokeh
    # **************************
    # Create KMeans with 8 clusters and fit data to model
    kmeans = KMeans(n_clusters=7)
    kmeans.fit_transform(tsne_features)
    labels = kmeans.predict(tsne_features)
    centroids = kmeans.cluster_centers_

    # **************************
    # write out labels for use in boxplots/table
    # **************************

    # tsne_data['Label'] = labels
    # tsne_data.to_csv('~/PycharmProjects/network_classification/src/data/tsne_label_data.csv')

    # Assign the columns of centroids: centroids_x, centroids_y
    centroids_x = centroids[:, 0]
    centroids_y = centroids[:, 1]

    df_labels = pd.read_csv('~/Downloads/network_classification/src/data/tsne_label_data.csv')

    # Create cross tabulation of tsne data and print
    df1 = pd.DataFrame({'labels': df_labels['Label'], 'Collection': df_labels['Category Name']})
    print("Crosstab for t-SNE data:\n")
    ct = pd.crosstab(df1['Collection'], df1['labels'])
    print(ct)

    data = {'x': xs, 'y': ys, 'Collection': category, 'Graph': names}
    df = pd.DataFrame(data)

    df.to_csv('~/Downloads/network_classification/src/data/tsne_canberra_coord_700_2.csv')

    # Create hover tool
    hover = HoverTool()
    hover.tooltips = [("Graph", "@Graph"), ("Category", "@{Collection}")]

    # Creating the scatter plot of the x and y coordinates
    m = figure(title = 't-SNE using Canberra Distance', plot_width=1000)

    # Color the plot by collection
    category = df['Collection']
    all_categories = category.unique().tolist()
    for i, graph in enumerate(all_categories):
        source = ColumnDataSource(df[df['Collection'] == graph])
        m.circle(x='x', y='y', source=source, color=d3['Category20'][16][i], size=8, legend=graph)

    # Creating scatter plot of centroids
    m.square_cross(centroids_x, centroids_y, color ='black', size = 12, legend = 'Centroid')

    # Add tools and interactive legend
    m.add_tools(hover)
    m.legend.location = "top_left"
    m.legend.click_policy = "hide"
    # Save file and show plot
    output_file(str(l) + '.html')
    show(m)

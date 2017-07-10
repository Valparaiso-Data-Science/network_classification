import pandas as pd

# Read file
df = pd.read_csv('~/Downloads/network_classification/data/raw/network-meta-statistics2.csv')

# Change column names
df.columns = (['Graph', 'Collection', 'Nodes', 'Edges', 'Density', 'Maximum degree',
       'Minimum degree', 'Average degree', 'Assortativity', 'Total triangles',
       'Average triangles', 'Maximum triangles', 'Avg. clustering coef.',
        'Frac. closed triangles', 'Maximum k-core', 'Max. clique (lb)', 'Chromatic Number'])

# Get only the miscellaneous data
misc_data = df[df['Collection'] == 'Miscellaneous Networks']

# Move data to a csv file
misc_data.to_csv('~/Downloads/network_classification/src/data/miscellaneous_networks.csv')
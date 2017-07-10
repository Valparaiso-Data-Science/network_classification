import pandas as pd
import matplotlib.pyplot as plt


# Read file
df = pd.read_csv('~/Downloads/network_classification/src/data/data_minmaxscale.csv', index_col=0)

# Plot pie chart with Collection column of df
df['Collection'].value_counts().plot.pie(fontsize = 10.5, autopct='%1.0f%%')
plt.ylabel('')

# Show plot
plt.show()
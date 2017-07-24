import pandas as pd
import matplotlib.pyplot as plt

# Read file
df = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_minmaxscale.csv', index_col=0)
df_old = pd.read_csv('~/PycharmProjects/network_classification/src/data/data_without_dimacs_or_bhoslib.csv', index_col=0)

# Plot pie chart with Collection column of df
#df['Collection'].value_counts().plot.pie(fontsize = 10.5, autopct='%1.0f%%')
#plt.ylabel('')
#
## Show plot
#plt.show()
my_colors = 'rgbkymc'
df_old['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)
plt.title('Counts of Collections Before Downsampling', fontsize=25)
plt.show()
##plot bar chart with collection column of df
df['Collection'].value_counts().plot.bar(color=my_colors, fontsize=15)
plt.title('Counts of Collections After Downsampling', fontsize=25)

# Show plot
plt.show()
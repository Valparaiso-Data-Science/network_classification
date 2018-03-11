# https://docs.python.org/2/library/csv.html
import csv

with open('miscNetworks.csv') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        fileName = row['Metadata']
        file = open(fileName+".bib","w")
        file.write("description network= "+row['Description Network']+"\n"+"collection= "+row['Collection']+ "\n"+ "tags="+ row['Tags']+"\n"+"edge weight= "+row['Edge Weights']+"\n"+ "Description= "+row['Description']+"\n"+"Citation= \""+row['Citation/Link']+'\"')
        file.close()
 

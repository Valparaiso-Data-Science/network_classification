# This code reads the miscellaneous.csv file 
# Converts each graph in the csv into a bibTex file

# https://docs.python.org/2/library/csv.html
import csv
import requests 

def doi2bib(doi):
  """
  Return a bibTeX string of metadata for a given DOI.
  """

  url = "http://dx.doi.org/" + doi

  headers = {"accept": "application/x-bibtex"}
  r = requests.get(url, headers = headers)

  return r.text
#print(doi2bib("10.1016/j.jctb.2006.04.002"))

with open('miscNetworks.csv') as csvfile:
	reader = csv.DictReader(csvfile)
	for row in reader:
		fileName = row['Metadata']
		file = open(fileName+".bib","w")
		file.write("description network= "+row['Description Network']+"\n"+"collection= "+row['Collection']+ "\n"+ "tags="+ row['Tags']+"\n"+"edge weight= "+row['Edge Weights']+"\n"+ "Description= "+row['Description']+"\n"+"Citation= \""+row['Citation/Link']+'\"'+"\n"+"\n")
		file.write(doi2bib(row['DOI']))
		file.close()
		 

import requests

def doi2bib(doi):
  """
  Return a bibTeX string of metadata for a given DOI.
  """

  url = "http://dx.doi.org/" + doi

  headers = {"accept": "application/x-bibtex"}
  r = requests.get(url, headers = headers)

  return r.text
print(doi2bib("10.1016/j.jctb.2006.04.002"))
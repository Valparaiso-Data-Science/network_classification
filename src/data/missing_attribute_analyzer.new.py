# Searches a CSV file and finds missing attributes
#
import csv
import os
import json
import fnmatch
import networkx as nx
import numpy as np

# data = np.genfromtxt('/home/cmorris5/Downloads/network-meta-statistics2.csv.xlsx', skip_header=0, delimiter=',', filling_values="0", dtype=None)
# np.savetxt('teststats.csv', data, delimiter=',', usecols = np.arrange(0,1000), fmt="%s")



infile = open('C:\\Users\Owner\\Documents\\VERUM\\Network stuff\\data_without_dimacs_or_bhoslib.csv', mode = 'rb')
reader = csv.reader(infile)
mydict = dict((rows[0], rows[1:]) for rows in reader)
print('mydict: ', mydict)

outfile = open('C:\\Users\\Owner\\Documents\\VERUM\\Network stuff\\testresults.csv', 'wb')
writer = csv.writer(outfile)
header = mydict['Graph']
completeRow = True  # used to see if a row is a complete row.  True by default
incompleteGraphs = []  # used to store incomplete graph names and attribute numbers
completeGraphs = []  # used to store complete graphs

# searches each key for missing data, if not missing data, it writes the complete row to the csv file.  If missing data, it add graph name and attributes missing to incomplete graphs
for keys in mydict:
    array = mydict[keys]  # may not be needed
    completeRow = True  # used to reset completeRow condition true for each iteration
    for i in range(len(array)):
        if array[i] is '':
            # print "The graph " + keys+ " is missing an the attribute " + header[i]
            completeRow = False
            incompleteGraphs.append(keys)
            incompleteGraphs.append(i)
    if completeRow is True:
        completeGraphs.append(keys)
        for j in array:
            completeGraphs.append(j)
        writer.writerow(completeGraphs)
        completeGraphs = []



print('incomplete Graphs: ', incompleteGraphs)
print(len(incompleteGraphs))
# queueCal(incompleteGraphs)

# searches a list (incompleteGraphs) for a string(graph name) and the numbers that follow(attributes to calculate), then queues the calculation
def queueCal(g):
    graph_name = None
    calculation = None
    nowComplete = False
    for i in g:
        if isinstance(i, str):
            graph_name = i
        elif isinstance(i, int):
            calculation = i
            # print "I am going to run calculation " + header[i] + " on " + str(graph_name)
            # x = doCalculation(graph_name,1)
            # print x
            global completeGraphs
            mydict[graph_name][i] = 777
            nowComplete = True
            for j in mydict[graph_name]:
                if j is "":
                    nowComplete = False
        if nowComplete is True:
            completeGraphs.append(graph_name)
            for k in mydict[graph_name]:
                completeGraphs.append(k)
            writer.writerow(completeGraphs)
            completeGraphs = []
            nowComplete = False


#def doCalculation(g, c):
 #   for root, dirs, files in os.walk("/export/home/cmorris5/Downloads/"):
  #      for names in root:
   #         if names.split(".")[0] is str(g):
    #            return "Yay, I found " + g
     #       else:
      #          return files


#queueCal(incompleteGraphs)
#xyz = doCalculation("soc-orkut", 1354)
#print(xyz)


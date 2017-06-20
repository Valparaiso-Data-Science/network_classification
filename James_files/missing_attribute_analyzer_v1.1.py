# Searches a CSV file and finds missing attributes
# Written by Charles Morris
# dictMax, dictMin, dictAverage taken from code by Casey Primozic
import csv
import os
import networkx as nx
import threading
import queue
import time

num_worker_threads = 2  # Adjust to desired number of threads
out_path = 'C:/Users/Owner/Documents/VERUM/Network stuff/testresults.csv'  # designate your own output path
# network_path = "/export/home/cmorris5/Downloads" #designate the location of the network path to walk (or can I add a way to do it from the current location by default?)

infile = open('C:/Users/Owner/Documents/VERUM/Network stuff/data_without_dimacs_or_bhoslib.csv',
              mode='rb')  # designate the location of the CSV file to read
reader = csv.reader(infile)
mydict = dict((rows[0], rows[1:]) for rows in reader)
infile.close()

queueAtt = []  # A list that contains a list of graphs and number corresponding to a calculation in doCalculation
incompleteGraphs = []  # used to store incomplete graph names and attribute numbers
completeGraphs = []  # used to store complete graphs

outfile = open(out_path, 'wb')
writer = csv.writer(outfile)
header = mydict['Graph']
completeRow = True
lock = threading.Lock()

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

outfile.close()


# searches a list (incompleteGraphs) for a string(graph name) and the numbers that follow(attributes to calculate)
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
            queueAtt.append([graph_name, i])
            global completeGraphs
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


# Takes a list from queueAtt and runs the correct calculation on the correct graph
def doCalculation(g):
    for root, dirs, files in os.walk(network_path):
        for names in files:
            if names.split(".")[0] == g[0]:
                graph = nx.read_adjlist(root + "/" + g[0] + ".txt")
                print("Calculating " + header[g[1]] + " for " + g[0])
                if g[1] is 0:
                    mydict[g[0]][g[1]] = "none"
                if g[1] is 1:
                    mydict[g[0]][g[1]] = nx.number_of_nodes(graph)
                elif g[1] is 2:
                    mydict[g[0]][g[1]] = nx.number_of_edges(graph)
                elif g[1] is 3:
                    mydict[g[0]][g[1]] = nx.density(graph)
                elif g[1] is 4:
                    newDict = dict(nx.degree(graph))
                    mydict[g[0]][g[1]] = dictMax(newDict)
                elif g[1] is 5:
                    newDict = dict(nx.degree(graph))
                    mydict[g[0]][g[1]] = dictMin(newDict)
                elif g[1] is 6:
                    newDict = dict(nx.degree(graph))
                    mydict[g[0]][g[1]] = dictAverage(newDict)
                elif g[1] is 7:
                    mydict[g[0]][g[1]] = nx.degree_assortativity_coefficient(graph)
                elif g[1] is 8:
                    newDict = dict(nx.triangles(graph))
                    mydict[g[0]][g[1]] = dictTotal(newDict)
                elif g[1] is 9:
                    newDict = dict(nx.triangles(graph))
                    mydict[g[0]][g[1]] = dictAverage(newDict)
                elif g[1] is 10:
                    newDict = dict(nx.triangles(graph))
                    mydict[g[0]][g[1]] = dictMax(newDict)
                elif g[1] is 11:
                    mydict[g[0]][g[1]] = nx.average_clustering(graph)
                elif g[1] is 12:  # Frac Closed Triangles (dont know formula)
                    mydict[g[0]][g[1]] = 'na'
                elif g[1] is 13:
                    mydict[g[0]][g[1]] = 'na'  # nx.k_core(graph)
                elif g[
                    1] is 14:  # max clique (max_clique doesnt exsist? and graph_clique_number returned different than in CSV file)
                    mydict[g[0]][g[1]] = 'na'
                elif g[1] is 15:  # Cant find Chromatic number in nx, will have to calculate upper bound?
                    mydict[g[0]][g[1]] = 'na'
                else:
                    mydict[g[0]][g[1]] = 'missing'


# Returns true or false after checking to see if a queued graph has been completed
def isComplete(graph):
    for keys in mydict[graph]:
        if keys is '':
            return False
    return True


# Dedicated to writing only completed queued graphs
def newWriter(graph):
    newCompleteGraphs = []
    with open(out_path, 'a') as outfile:
        writer = csv.writer(outfile)
        newCompleteGraphs.append(graph)
        for keys in mydict[graph]:
            newCompleteGraphs.append(keys)
        writer.writerow(newCompleteGraphs)


# Analysis Functions
def dictAverage(inDict):
    sum = 0
    for key in inDict:
        sum += inDict[key]
    if sum != 0:
        return sum / float(len(inDict))
    else:
        return 0


def dictMax(inDict):
    maxNum = 0
    if len(inDict) == 0:
        return 0
    for key in inDict:
        if maxNum < inDict[key]:
            maxNum = inDict[key]
    return maxNum


def dictMin(inDict):
    minNum = "null"
    if len(inDict) == 0:
        return 0
    for key in inDict:
        if minNum == "null":
            minNum = inDict[key]
        if minNum > inDict[key]:
            minNum = inDict[key]
    return minNum


def dictTotal(inDict):
    totalNum = 0
    if len(inDict) == 0:
        return 0
    for key in inDict:
        totalNum += inDict[key]
    return totalNum


queueCal(incompleteGraphs)  # generates list of graph and attributes needed to run


# Multithreading part of code
def worker():
    while True:
        item = q.get()
        if item is None:
            break
        doCalculation(item)
        if isComplete(item[0]) is True:
            lock.acquire()
            try:
                newWriter(item[0])
            finally:
                lock.release()
        q.task_done()


q = queue.Queue(0)
threads = []
for m in range(num_worker_threads):
    t = threading.Thread(target=worker)
    t.start()
    threads.append(t)
    time.sleep(.1)

for item in queueAtt:
    q.put(item)

q.join()

for i in range(num_worker_threads):
    q.put(None)
for t in threads:
    t.join()

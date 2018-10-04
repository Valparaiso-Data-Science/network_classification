GraphPACK: A Parallel Library of Graph Algorithms
===============================
In short, a parallel framework and library of graph algorithms for computing statistics and graph parameters in large-scale graphs. 

##Network Repository
This library is currently used to compute graph statistics for network repository (http://networkrepository.com) including global, distribution-based statistics, and other micro-level (local) statistics of vertices and edges. It is mainly used for extremely large graphs that are too slow to compute on the web/real-time using php/javascript graph codes. 


##Parallel Maximum Clique Enumeration

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem mce --threads 2

Another example:

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem mce -k 0 -t 2
	
	

##Parallel K-clique Enumeration

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem mce -k 30 --threads 2

Similarly:

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem kclique_enum -k 30 --threads 2
	
	



##Parallel Maximum Clique

This is the problem that graphpack defaults too!
Note that if -k is set, then k-clique problem is solved instead. 
Setting `-k 0` is fine, and indicates that you want to find the maximum clique.

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --threads 2

Similarly:

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem mc --threads 2





##Parallel K-clique

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem mc -k 30 --threads 2

Similarly:

	./graphpack --graph ../data/networks/socfb-Duke14.mtx --problem kclique -k 30 --threads 2
	
	
	
	

./graphpack --graph ../data/networks/soc-flickr.mtx -t 2 --problem mc --go kcore --lr hybrid --search_method repair --lo coloring --lp tcores --density_cutoff 0.95



Graph Coloring Algorithms
==============================================

## Coloring Algorithms
There are two main coloring algorithms:
+	vertex-centric: 	use a vertex-indexed array to assign colors directly to vertices, neighbors of vertices are explicitly searched
	-	well-suited for large sparse graphs, since O(|E|) if graph is connected
	
+	color-centric: 	vertices are explicitly assigned to color classes, color classes are explicitly searched
	- 	best for dense graphs, or coloring neighborhoods


Color-class centric coloring is specified as follows:

	./graphpack --graph graphs/soc-flickr.mtx --problem greedy_neigh_coloring --color_method color-centric

For vertex-centric coloring, use:

	./graphpack --graph graphs/soc-flickr.mtx --problem greedy_neigh_coloring --color_method vertex-centric

		--color_method color-centric
		--color_method vertex-centric
		
		
## Coloring Variants

There are two main coloring variants:

+	repair/recolor: 	if a vertex v is assigned a new color (never before used), an attempt is made to find a vertex in a previous color-class that can be swapped with a vertex in a color class formed after that vertex was assigned. If this situation arises, then the vertex v is swapped with the other vertex, and a new color class is avoided. Otherwise, v remains in the new color class.

+	simple/normal:	this procedure does not use the recolor algorithm, and simply assigns vertices via a sequential greedy coloring -- based on an ordering


For the greedy coloring method that uses repair/recolor algorithm, use:

	./graphpack --graph graphs/soc-flickr.mtx --problem greedy_neigh_coloring --color_method color-centric --color_type repair
	
Otherwise, for the normal greedy coloring without recolor, use:

	./graphpack --graph graphs/soc-flickr.mtx --problem greedy_neigh_coloring --color_method color-centric --color_type normal


## Vertex/edge ordering strategies

All the ordering method defined for maximum clique may be used to order the vertices/edges.
Most all ordering strategies are at most O(|E|) time.
The only exception is the ordering methods that use triangles. If the triangle counts are not precomputed, then it costs O(|E|^{3/2}) to compute triangle counts, and O(|E|) to order the vertices/edges using those counts.



## Neighborhood Coloring (Distance-Two Coloring)





Parameterized Multi-threaded Maximum Clique Framework: 
A Flexible Approach for Large Sparse and Dense Graphs
==============================================

##multi-threaded heuristic
+	ordering of vertices (order in which workers/processors search the vertices)
+	local neighborhood ordering: last vertex in ordering is greedily chosen to be explored

##pruning

##ordering
vertex ordering
edge/neighbor ordering (parallel)

##representation
+	global representation: optimizer decides automatically based on density and number of vertices
+	local representation: user may select csc only (sparse graphs), hybrid (csc + adj struct)
+	search representation: usually based on csc, but can be automatically adapted

##Search methods
+	repair methods
+	dynamic bounds: 


##Dynamic Local and Search Strategies
+	dynamic threshold: between 0 and 1, if 0, the bound is never used, if 1 then the bound is always used.
+	density threshold: between 0 and 1, if 0, the tightest local bound is used, if 1 then the bound is almost never used unless clique.


Parallel Maximum Clique (PMC) Library
=====================================

In short, a parameterized high performance library for computing maximum cliques in large sparse graphs.

Finding maximum cliques, k-cliques, and temporal strong components are in general NP-hard.
Yet, these can be computed fast in most social and information networks.
The PMC library is designed to be fast for solving these problems.
Algorithms in the PMC library are easily adaptable for use with a variety of orderings, heuristic strategies, and bounds.

* **Maximum clique:** 		 Given a simple undirected graph G and a number k, output the clique of largest size.
* **K-clique:** 	   		 In k-clique, the problem is to find a clique of size k if one exists.
* **Largest temporal-scc:** Given a temporal graph G, a temporal strong component is a set of vertices where all temporal paths exist between the vertices in that set. The Largest TSCC problem is to find the largest among all the temporal strong components.



Features
--------
0.  General framework for parallel maximum clique algorithms
1.	Optimized to be fast for large sparse graphs 
	+ 	Algorithms tested on networks of 1.8 billion edges
2.	Set of fast heuristics shown to give accurate approximations
3.	Algorithms for computing Temporal Strongly Connected Components (TSCC) of large dynamic networks
4.	Parameterized for computing k-cliques as fast as possible
5.  Includes a variety of tight linear time bounds for the maximum clique problem
6.  Ordering of vertices for each algorithm can be selected at runtime 
7.  Dynamically reduces the graph representation periodically as vertices are pruned or searched
	+   Lowers memory-requirements for massive graphs, increases speed, and has caching benefits


Synopsis
---------

### Setup
First, you'll need to compile the parallel maximum clique library.  

		$ cd path/to/graphpack/
		$ make

Afterwards, the following should work:  

		# compute maximum clique using the full algorithm `-a 0`
		./graphpack -f data/socfb-Texas84.mtx -a 0


*PMC* has been tested on Ubuntu linux (10.10 tested) and Mac OSX (Lion tested) with gcc-mp-4.7 and gcc-mp-4.5.4

Please let me know if you run into any issues.  

  
   
### Input file format
+ Matrix Market Coordinate Format (symmetric)  
For details see: <http://math.nist.gov/MatrixMarket/formats.html#MMformat>  

		%%MatrixMarket matrix coordinate pattern symmetric  
		4 4 6  
		2 1  
		3 1  
		3 2  
		4 1  
		4 2  
		4 3 


+ Edge list: These codes are designed to be as flexible as possible and accept many variations of edge lists. Note these codes may be slightly slower than the mtx reader. This is due to allowing flexible edge list formats. Hence, this reader must perform many checks to figure out the exact format of the params file, and performs any necessary preprocessing work that may be required.


	* Delimiters: The graphpack reader accepts comma, space, or tab delimited edge lists.
			
	* Comments: Comments are allowed and should be denoted with the first character of a newline as # or %

	* Weights: If an edge list contains weights on the 3rd column, they are simply ignored. A user may specify to read the weights by setting the wt parameter or by noting the graph is in fact a temporal graph.

	* Multigraph: When an edge list contains multiple edges, we simply remove the duplicate edges.

	* The edge list may also contain gaps in the vertex ids (non sequential vertex ids) and start from any positive integer. Self-loops are removed.

	* The edge list is assumed to be undirected. However, if a directed graph is given, it is simply treated as undirected.







Overview
---------

The parallel maximum clique algorithms use tight bounds that are fast to compute.
A few of those are listed below.

* K-cores
* Degree
* Neighborhood cores
* Greedy coloring

All bounds are dynamically updated.  

Examples of the three main maximum clique algorithms are given below.
Each essentially builds on the other.

	# uses the four basic k-core pruning steps
	./graphpack -f ../graphpack/data/output/socfb-Stanford3.mtx -a 2

	# k-core pruning and greedy coloring
	./graphpack -f ../graphpack/data/output/socfb-Stanford3.mtx -a 1
	
	# neighborhood core pruning (and ordering for greedy coloring)
	./graphpack -f ../graphpack/data/output/socfb-Stanford3.mtx -a 0





### Dynamic graph reduction
	
The reduction wait parameter `-r` below is set to be 1 second (default = 4 seconds).
 
	./graphpack -f data/sanr200-0-9.mtx -a 0 -t 2 -r 1

In some cases, it may make sense to turn off the explicit graph reduction. 
This is done by setting the reduction wait time '-r' to be very large.

	# Set the reduction wait parameter
	./graphpack -f data/socfb-Stanford3.mtx -a 0 -t 2 -r 999






### Orderings

The PMC algorithms are easily adapted to use various ordering strategies. 
To prescribe a vertex ordering, use the -o option with one of the following:
+ `deg`
+ `kcore`
+ `dual_deg`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;orders vertices by the sum of degrees from neighbors
+ `dual_kcore`&nbsp;&nbsp;orders vertices by the sum of core numbers from neighbors
+ `kcore_deg`&nbsp;&nbsp;&nbsp; vertices are ordered by the weight k(v)d(v)
+ `rand`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; randomized ordering of vertices



##### Direction of ordering

Vertices are searched by default in increasing order, to search vertices in decreasing order, use the `d` option:

	./graphpack -f data/p-hat700-2.mtx -a 0 -d




### Heuristic
The fast heuristic may also be customized to use various greedy selection strategies.
This is done by using `-h` with one of the following: 

+ `deg`
+ `kcore`
+ `kcore_deg`&nbsp;&nbsp;&nbsp; select vertex that maximizes k(v)d(v)
+ `rand`&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp; randomly select vertices


#### Terminate after applying the heuristic
Approximate the maximum clique using _ONLY_ the heuristic by not setting the exact algorithm via the `-a [num]` option.
For example:  

	./graphpack -f data/sanr200-0-9.mtx -h deg
	
#### Turning the heuristic off

	# heuristic is turned off by setting `-h 0`.
	./graphpack -f data/tscc_enron-only.mtx -h 0 -a 0



### K-clique

The parallel maximum clique algorithms have also been parameterized to find cliques of size k.
This routine is useful for many tasks in network analysis such as mining graphs and community detection.

	# Computes a clique of size 50 from the Stanford facebook network
	./graphpack -f data/socfb-Stanford3.mtx -a 0 -k 50


using `-o rand` to find potentially different cliques of a certain size

	# Computes a clique of size 36 from sanr200-0-9
	./graphpack -f data/sanr200-0-9.mtx -a 0 -k 36 -o rand


Documentation
--------------
To generate the document you must have doxygen and graphviz installed. On Mac OSX these can be stalled using *homebrew* with the following commands:

	# install doxygen and graphviz using homebrew on Mac OSX
	brew install doxygen
	brew install graphviz

Afterwards, the documentation is generated by simply typing `make docs` in the root directory of graphpack.
This creates the `./html` directory with the documentation.



Note that in the unlikely case that *Doxyfile* does not exist, simply type `doxygen -g` to generate a default config file.


Graph query optimization
------------------------

The graphpack uses a graph query optimizer uses fast graph statistics to perform the query as fast as possible.
In one simple example, suppose we are given a query that accesses relatively dense parts of the graph and is focused on checking neighbors or the existence of an edge, then our query optimizer would represent this subgraph as an adjacency matrix since it allows for O(1) checks on the existence of an edge, whereas if we had used compressed-sparse-columns then it would require log(n) using a binary search.
The query optimizer automatically selects the fastest plan for a given query. 
For queries on graphs, the fastest plan is highly dependent on the following:

* the properties of the graph itself, 
* the query, 
* the graph operations required for that query, 
* and the data structures used to represent the graph in memory.


Terms and conditions
--------------------
Please feel free to use these codes. We only ask that you cite:  

	TODO 

_These codes are research prototypes and may not work for you. No promises. But do email if you run into problems._


Copyright 2011-2017, Ryan A. Rossi and Nesreen K. Ahmed, All rights reserved.
# GraphPACK
# GraphPACK

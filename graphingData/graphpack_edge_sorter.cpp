/**
 ======================================================================================
 Name        : Graph PACKage (graphpack) - A Parallel Graph Algorithms Library
 Author      : Ryan A. Rossi (http://ryanrossi.com)
 	 	   Nesreen K. Ahmed (http://nesreenahmed.com)
 Description : A general high-performance parallel framework of graph algorithms for
 	 	   computing graph statistics of large graphs. The library is designed to
 	 	   be fast and flexible for large sparse graphs and dense graphs.

 Copyright (C) 2012-2017, Ryan A. Rossi and Nesreen Ahmed, All rights reserved.

 See http://www.GraphPACK.com for more information.
 ======================================================================================
 */

#include <cstdlib>
#include <cstring>
#include <vector>
#include "graphpack_graph.h"

using namespace std;
using namespace graphpack;


void graphpack_graph::degree_bucket_sort_parallel() {
    degree_bucket_sort(false);
}

/*
 * Each worker sorts the edges of a vertex, and
 * afterwards, grabs the next available vertex to sort
 *
 * Note: sort neighbors by degree (largest to smallest)
 */
void graphpack_graph::degree_bucket_sort_parallel(bool desc) {

//	int v, u, n, md, md_end, start, d, num;

	vector<int> tmp_edges(edges.size(),0);
//	tmp_edges.reserve(edges.size());
	int v = 0;
	double sec = get_time();
	#pragma omp parallel for schedule(dynamic) \
		shared(tmp_edges)
//	firstprivate(colors,ind,vs,es) \
//	private(u, P, C, C_max, mc_cur)
	for (v = 0; v < num_vertices(); ++v) {

		int u, n, md, md_end, start, d, num;

		n = vertices[v+1] - vertices[v] + 1;
		vector<int> vert(n);
		vector<int> pos(n);
		vector<int> deg(n);

		md = 0;
		for(u=1; u<n; u++) {
			deg[u] = degree[edges[vertices[v] + (u-1)]];
			if (deg[u] > md)
				md = deg[u];
		}

		md_end = md+1;
		vector < int > bin(md_end,0);

		for (u=1; u < n; u++)  bin[deg[u]]++;

		start = 1;
		for (d=0; d < md_end; d++) {
			num = bin[d];
			bin[d] = start;
			start = start + num;
		}

		for (u=1; u<n; u++) {
			pos[u] = bin[deg[u]];
			vert[pos[u]] = edges[vertices[v] + (u-1)];
			bin[deg[u]]++;
		}

//		int insert_pos = vertices[v];
		if (desc) {
			for (u=n-1; u>0; --u) {
				tmp_edges[vertices[v] + (u-1)] = vert[u];
			}

			cout << print_break();
			cout << "v = " << v <<endl;
			for (int ii=0; ii < num_vertices(); ii++)
			    for (u=vertices[ii]; u < vertices[ii+1]; u++) {
			        cout << "\t(" << ii << "," << tmp_edges[u] << ")" <<endl;
			    }

			// largest to smallest
//			tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.rbegin(),vert.rend()-1);
		}
		else {
			for (u=1; u<n; u++) {
				tmp_edges[vertices[v] + (u-1)] = vert[u];
			}

			//from smallest degree to largest
//			tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.begin()+1,vert.end());
		}
	}

	cout << "[mcpack: sorting neighbors]  |E| = " << edges.size();
	cout << ", |E_sorted| = " << tmp_edges.size() <<endl;
	cout << "sorting neighbors took " << get_time() - sec <<endl;
	edges = tmp_edges;
}


void graphpack_graph::degree_bucket_sort() {
    degree_bucket_sort(false);
}

// sort neighbors by degree (largest to smallest)
void graphpack_graph::degree_bucket_sort(bool desc) {

	int v, u, n, md, md_end, start, d, num;

	vector<int> tmp_edges;
	tmp_edges.reserve(edges.size());

//#pragma omp parallel for schedule(dynamic) \
//	shared(tmp_edges) \
//	firstprivate(colors,ind,vs,es) \
//	private(u, P, C, C_max, mc_cur)
	for (v = 0; v < num_vertices(); v++) {

		n = vertices[v+1] - vertices[v] + 1;
		vector<int> vert(n);
		vector<int> pos(n);
		vector<int> deg(n);

		md = 0;
		for(u=1; u<n; u++) {
			deg[u] = degree[edges[vertices[v] + (u-1)]];
			if (deg[u] > md)
				md = deg[u];
		}

		md_end = md+1;
		vector < int > bin(md_end,0);

		for (u=1; u < n; u++)  bin[deg[u]]++;

		start = 1;
		for (d=0; d < md_end; d++) {
			num = bin[d];
			bin[d] = start;
			start = start + num;
		}

		for (u=1; u<n; u++) {
			pos[u] = bin[deg[u]];
			vert[pos[u]] = edges[vertices[v] + (u-1)];
			bin[deg[u]]++;
		}

		if (desc) {
			// largest to smallest
			tmp_edges.insert(tmp_edges.end(),vert.rbegin(),vert.rend()-1);
		}
		else {
			//from smallest degree to largest
			tmp_edges.insert(tmp_edges.end(),vert.begin()+1,vert.end());
		}
	}

	cout << "[mcpack: sorting neighbors]  |E| = " << edges.size();
	cout << ", |E_sorted| = " << tmp_edges.size() <<endl;
	edges = tmp_edges;
}

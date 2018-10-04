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

#ifndef GRAPHPACK_MCE_UTILS_H_
#define GRAPHPACK_MCE_UTILS_H_

#include "graphpack_vertex.h"
#include <algorithm>
#include <set>

using namespace std;

namespace graphpack {
    inline
    static void add_clique(vector<int>& C, set< vector<int> >& cliques) {
        vector<int> clq(C.size());
        for (int j = 0; j < C.size(); j++) {
            clq[j] = C[j];
        }
        std::sort(clq.begin(), clq.end());
        cliques.insert(clq);
        clq.clear();
    }

    inline
    static void add_all_cliques(set< vector<int> >& cliques_local, set< vector<int> >& cliques) {

        //        cliques.insert(cliques_local);

        set< vector<int> >::iterator it;
        for( it = cliques_local.begin(); it != cliques_local.end(); it++) {
            const vector<int>& clq = (*it);
            cliques.insert(clq);
        }
    }

    static void print_all_maxcliques(set< vector<int> >& cliques) {
        set< vector<int> >::iterator it;
        cout << "size = " << cliques.size() << endl;
        for( it = cliques.begin(); it != cliques.end(); it++) {
            const vector<int>& clq = (*it);
            for (int j = 0; j < clq.size(); j++)
                cout << clq[j] << " ";
            printf("\n\n");
        }
    }

    static void check_maxcliques(graphpack_graph& G, set< vector<int> >& cliques, int & w) {

        set< vector<int> >::iterator it;
        for (it = cliques.begin(); it != cliques.end(); it++) {
            const vector<int>& clq = (*it);
            if (clq.size() != w) {
                printf("****************************************** \n");
                printf("ERROR: clique is NOT a maximum clique \n");
                printf("****************************************** \n");
            }
        }
    }

    static void add_maxclique_vertices(graphpack_graph &H, set< vector<int> >& C, int & w, vector<Vertex> &P) {
        printf("num vertices = %d \n", H.num_vertices());
        vector<int> ind(H.num_vertices());
        set< vector<int> >::iterator it;
        for( it = C.begin(); it != C.end(); it++) {
            const vector<int>& clq = (*it);
            for (int j = 0; j < clq.size(); j++) {
                int v = clq[j];
                if (ind[v] == 0) { // v has not yet been added
                    ind[v] = 1; // mark v
                    P.push_back(Vertex(v,0));
                }
            }
        }

        printf("number of unique vertices in the maximum cliques = %lu, omega = %d \n",P.size(), w);
    }

    inline static void print_mc_info(vector<int> &C_max, double &sec, bool verbose = true) {
        if (verbose) {
        cout << "*** [graphpack: thread " << omp_get_thread_num() + 1;
        cout << "]   current max clique = " << C_max.size();
        cout << ",  time = " << get_time() - sec << " sec" <<endl;
        }
    };

}
#endif

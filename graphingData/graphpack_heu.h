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

#ifndef GRAPHPACK_HEU_H_
#define GRAPHPACK_HEU_H_

#include "graphpack_headers.h"
#include "graphpack_graph.h"
#include "graphpack_utils.h"
#include "graphpack_params.h"
#include "graphpack_vertex.h"
#include "graphpack_framework.h"
#include "graphpack_neigh_coloring.h"
#include "graphpack_mce_utils.h"
#include <algorithm>

namespace graphpack {

    class graphpack_heu {
        public:
            vector<int>* E;
            vector<long long>* V;
            vector<int>* K;
            vector<int>* order;
            vector<int>* degree;
            vector<long long>* triangles;
            double sec;
            int ub;
            int block_size;
            int local_ordering;
            int global_ordering;
            bool global_small_to_large;
            bool local_small_to_large;
            bool verbose;

            graphpack_heu(graphpack_graph& G, params& params) {
                K = G.get_kcores();
                order = G.get_kcore_ordering();
                if (G.triangles.size() > 0) {
                    triangles = &G.triangles;
                }
                ub = params.ub;
//                strat = params.heu_strat;
                initialize();

                verbose = params.verbose;
//                global_ordering = set_ordering_enum(params.heu_global_ordering);
//                local_ordering = set_ordering_enum(params.heu_local_ordering);

                global_ordering = set_enum(params.heu_global_ordering);
                local_ordering = set_enum(params.heu_local_ordering);

                global_small_to_large = params.heu_global_s2l;
                local_small_to_large = params.heu_local_s2l;

                srand (time(NULL));

                if (local_ordering == DEGREE && G.degree.size() == 0) {
                    G.vertex_degrees();
                }
                else if (local_ordering == KCORE && G.kcore.size() == 0) {
                    G.compute_cores();
                }
                else if (local_ordering == TRIANGLES && G.t.size() == 0) {
                    printf("[ heuristic clique finder ]  precomputing triangles for local ordering\n");
                    G.compute_vertex_triangles(0,true);
                }
                else if (local_ordering == TRIANGLE_CORES && G.tri_core.size() == 0) {
                    printf("[ heuristic clique finder ]  precomputing triangle cores for local ordering\n");
                    G.compute_triangle_cores();
                }

                block_size = params.block_size;

            }

            graphpack_heu(graphpack_graph& G, int tmp_ub) {
                K = G.get_kcores();
                order = G.get_kcore_ordering();
                if (G.triangles.size() > 0) {
                    triangles = &G.triangles;
                }
                ub = tmp_ub;
//                strat = "kcore";
                initialize();
//                block_size = params.block_size;

            }

            void initialize() {
                sec = get_time();
                srand (time(NULL));
                block_size = 1;
            };

            int strategy(vector<int>& P);
            void set_strategy(string s) { local_ordering = set_ordering_enum(s); }
            int compute_heuristic(int v, graphpack_graph & G);

            static bool desc_heur(Vertex v,  Vertex u) {
                return (v.get_bound() > u.get_bound());
            }

            static bool incr_heur(Vertex v,  Vertex u) {
                return (v.get_bound() < u.get_bound());
            }

            int search(graphpack_graph& graph, vector<int>& C_max);
            int search(graphpack_graph& graph, vector<int>& C_max, params & p);
            int search_cores(graphpack_graph& graph, vector<int>& C_max, int lb);
            int search_bounds(graphpack_graph& graph, vector<int>& C_max);
            int search_bounds_ordering(graphpack_graph& G, vector<int>& C_max, params & p);
            int search_for_large_clique(graphpack_graph& G, vector<int>& C_max, params & p);

            int clique_compression(graphpack_graph& G,
                    vector<int>& C_max,
                    vector<int> & vertex_order_by_color,
                    set< vector<int> >& cliques,
                    params & p);

            int clique_compression_distance_two(graphpack_graph& G,
                    vector<int>& C_max,
                    vector<int> & vertex_order_by_color,
                    set< vector<int> >& cliques,
                    params & p);

            void branch(vector<Vertex>& P, int sz,
                    int& mc, vector<int>& C, vector<int>& ind);

            void clique_search(
                    vector<Vertex> &P,
                    vector<int>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    int& mc);

            void print_info(vector<int> C_max);

            int search_sparse(graphpack_graph& G, vector<int>& sol);
            void branch_sparse(vector<Vertex> &P,
                    vector<int>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc);

            int search_dense(graphpack_graph& G, vector<int>& sol);
            void branch_dense(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<int>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc,
                    bool** &adj);
    };
};
#endif

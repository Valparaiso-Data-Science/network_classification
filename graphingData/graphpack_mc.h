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

#ifndef GRAPHPACK_MC_H_
#define GRAPHPACK_MC_H_

#include <cstddef>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include "graphpack_headers.h"
#include "graphpack_utils.h"
#include "graphpack_graph.h"
#include "graphpack_params.h"
#include "graphpack_vertex.h"
#include "graphpack_neigh_cores.h"
#include "graphpack_neigh_coloring.h"
#include "graphpack_mce_utils.h"
#include "graphpack_framework.h"
#include <algorithm>

using namespace std;

namespace graphpack {

    class graphpack_mc {
        public:

            class OpCount {
                int i1, i2; //global vars
                public:
                    OpCount() : i1(0), i2(0) {} //constructor
                    void set_i1(const int ii)   { i1 = ii;      }
                    int get_i1() const          { return i1;    }
                    void set_i2(const int ii)   { i2 = ii;      }
                    int get_i2() const          { return i2;    }
                    void inc_i1()               { i1++;         }
            };

            vector<int>* edges;
            vector<long long>* vertices;
            vector<int>* bound;
            vector<int>* order;
            vector<int>* degree;

            set< vector<int> > cliques;

            int param_ub;
            int k_clique_size;
            int ub;
            int lb;
            int steps;
            int num_threads;
            int block_size;


            double time_limit;
            double sec;
            double wait_time;
            double threshold;


            bool not_reached_ub;
            bool time_expired_msg;
            bool is_enum;

            int edge_ordering;

            int global_ordering;
            int local_ordering;
            int search_ordering;

            int global_pruning;
            int local_pruning;
            int search_pruning;

            int global_representation;
            int local_representation;
            int search_representation;

            int search_method;
            int coloring_method;

            bool global_small_to_large;
            bool local_small_to_large;
            bool search_small_to_large;

            bool use_neigh_cores;



            void initialize() {
                steps = 0;
                edge_ordering = 0;
                not_reached_ub = true;
                time_expired_msg = true;

                block_size = 1;

                threshold = 0.025;
                is_enum = false;

                global_small_to_large = true;
                local_small_to_large = false;
                search_small_to_large = false;

            }

            void setup_bounds(params& params) {
                lb = params.lb;
                ub = params.ub;
                k_clique_size = params.k;
                if (k_clique_size == 0)  {
                    k_clique_size = ub;
                }
                time_limit = params.time_limit;
                wait_time = params.remove_time;
                sec = get_time();

                num_threads = params.threads;
                block_size = params.block_size;
                is_enum = params.is_enum;
            }

            graphpack_mc(graphpack_graph& G, params& params) {
                bound = G.get_kcores();
                order = G.get_kcore_ordering();
                setup_bounds(params);
                initialize();

                global_ordering = set_enum(params.global_ordering);
                local_ordering = set_enum(params.local_ordering);
                search_ordering = set_enum(params.search_ordering);

                cout << "global_ordering: " << global_ordering << ",  local_ordering: " << local_ordering << endl;

                global_small_to_large = params.global_small_to_large;
                local_small_to_large = params.local_small_to_large;
                search_small_to_large = params.search_small_to_large;

                global_pruning = set_enum(params.global_pruning);
                local_pruning = set_enum(params.local_pruning);
                search_pruning = set_enum(params.search_pruning);

                global_representation = set_enum_representation(params.global_representation);
                local_representation = set_enum_representation(params.local_representation);
                search_representation = set_enum_representation(params.search_representation);

                search_method = set_enum_search_method(params.search_method);
                coloring_method = set_enum_coloring_method(params.coloring_method);

                threshold = params.dynamic_cutoff;
                use_neigh_cores = params.use_neigh_cores;

                if (G.verbose) {
                    printf("search pruning = %d, ordering = %d, representation = %d, method = %d \n",
                            search_pruning, search_ordering, search_representation, search_method);
                }

            }

            ~graphpack_mc() {};


            /**
             * Parallel Maximum Clique Framework
             */
            int maxclique(graphpack_graph& G, vector<int>& sol, params & p);

            int maxclique_edge_neighborhood(graphpack_graph& G, vector<int>& sol, params & p);

            int greedy_coloring_neighborhoods(graphpack_graph& G, vector<int>& sol, params & p);
            int greedy_coloring_entire_graph(graphpack_graph& G, vector<int>& sol, params & p);

            //
            // Only CSC or Global Adj
            //
            void branch(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc);


            void branch_basic(
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    int* &pruned,
                    int& mc);

            void branch_static_order(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);

            void branch_repair(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);

            void branch_repair_only(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);



            void branch_dynamic(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps);


            //
            // Hybrid Global CSC + Adj
            //
            void branch_dense(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc,
                    bool** &adj);

            void branch_basic_dense(
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    int* &pruned,
                    int& mc,
                    bool** &adj);

            void branch_dense_static_order(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj);


            void branch_dense_repair(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj);

            void branch_dense_repair_only(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj);

            void branch_dense_dynamic(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps);



            //
            // Explicitly Induced Neighborhoods
            //
            void branch_neighborhood_repair(
                    graphpack_graph & N,
                    vector<Vertex> & P,
                    vector<Vertex> & ColOrd,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc);

            void branch_neighborhood_repair_only(
                    graphpack_graph & N,
                    vector<Vertex> & P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc);

            void branch_neighborhood_static_order(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<Vertex> &ColOrd,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc);


            void branch_neighborhood_cores(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);


            void branch_dynamic_neighborhood_cores(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps);

            void branch_dynamic_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps);


            void branch_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);


            void branch_basic_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max);


            /**
             * Neighborhood coloring induced
             */
            void neigh_coloring_induced(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<int>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int& mc,
                    bool** &adj);





            //
            // MAXIMUM CLIQUE ENUMERATION and K-CLIQUE ENUMERATION
            //

            int search_enum(graphpack_graph& G, set< vector<int> >& cliques, params &p);

            void branch_neighborhood_enum(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<short> &ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void branch_enum(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc,
                    vector< set< vector<int> > > &cliques_local);

            void branch_dense_enum(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);











            //
            // Only CSC or Global Adj
            //
            void clique_enum(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);


            void clique_enum_basic(
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_static_order(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_repair(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_repair_only(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);



            void clique_enum_dynamic(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps,
                    vector< set< vector<int> > > &cliques_local);


            //
            // Hybrid Global CSC + Adj
            //
            void clique_enum_dense(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& mc,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_basic_dense(
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    int* &pruned,
                    int& mc,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_dense_static_order(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);


            void clique_enum_dense_repair(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<Vertex> & ColOrd,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_dense_repair_only(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<short>& ind,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_dense_dynamic(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    bool** &adj,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps,
                    vector< set< vector<int> > > &cliques_local);



            //
            // Explicitly Induced Neighborhoods
            //
            void clique_enum_neighborhood_repair(
                    graphpack_graph & N,
                    vector<Vertex> & P,
                    vector<Vertex> & ColOrd,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_neighborhood_repair_only(
                    graphpack_graph & N,
                    vector<Vertex> & P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

            void clique_enum_neighborhood_static_order(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<Vertex> &ColOrd,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);



            void clique_enum_dynamic_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector<short>& ind,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< vector<OpCount> > &S,
                    vector<int> &levels,
                    vector<int> &steps,
                    vector< set< vector<int> > > &cliques_local);


            void clique_enum_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);



            void clique_enum_basic_neighborhood(
                    graphpack_graph & N,
                    vector<Vertex> &P,
                    vector<int>& C,
                    vector<int>& C_max,
                    vector< vector<int> >& colors,
                    int* &pruned,
                    int& max,
                    vector< set< vector<int> > > &cliques_local);

    };
};

#endif

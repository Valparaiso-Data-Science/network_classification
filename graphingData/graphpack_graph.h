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
#ifndef GRAPHPACK_GRAPH_H_
#define GRAPHPACK_GRAPH_H_


//#include <stdlib.h>     /* atoi */
#include <float.h>
#include <cstddef>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include <limits>
#include "math.h"
#include "graphpack_headers.h"
#include "graphpack_utils.h"
#include "graphpack_vertex.h"
//#include "graphpack_framework.h"


namespace graphpack {

    class graphpack_graph {
        private:

        public:
            // helper functions
            int read_mtx(const string& filename); // general mtx file and error-correcting
            int read_mtx_file(const string& filename);
            int read_edges(const string& filename);
            int read_metis(const string& filename);

            vector<int> edges;
            vector<long long> vertices;
            vector<double> wt;


            vector<int> degree;

            vector<int> bound;

//            vector<long long> tris_order; // triangle core ordering

            int min_degree;
            int max_degree;
            int avg_degree;
            bool is_gstats;
            string fn;
            bool** adj;
            bool is_dense;
            bool is_weighted;

            bool verbose;

            bool block_size;

            // user needs to specify if the actual vertex ids are important,
            // otherwise we discard them
            bool use_actual_labels;



            /**
             * @brief: for induced graph constructor, maps the vertex_id in H to its original id in G.
             * for v \in H (induced graph), vertex_lookup[v] returns previous vertex id (if any)
             * This is useful if the graph has been induced from another,
             * then can be used to get back to original vertex ids
             */

            vector<int> vertex_lookup;



            // constructor
//            graphpack_graph() {};
            graphpack_graph(const string& filename);
            graphpack_graph(bool graph_stats, const string& filename);
            graphpack_graph(const string& filename, bool make_adj);
            graphpack_graph(vector<long long> vs, vector<int> es) {
                edges = es;
                vertices = vs;
                vertex_degrees();
            }

//            graphpack_graph(graphpack_graph &G, vector<Vertex>& C);
            graphpack_graph(graphpack_graph &G, vector<int>& C);

            // ind is a map from v --> v', used in neighborhood branch for P
            graphpack_graph(vector<Vertex>& U, vector<int> & ind, vector<long long>& vs, vector<int>& es);

            graphpack_graph(vector<Vertex>& U, vector<long long>& vs, vector<int>& es);
            graphpack_graph(vector<int>& U, vector<long long>& vs, vector<int>& es);

            graphpack_graph(vector<int>& U, graphpack_graph& G);
            graphpack_graph(vector<Vertex>& U, graphpack_graph& G);
            // destructor
            ~graphpack_graph();

            graphpack_graph(int nverts, int nedges,
                    std::pair<int,int>* edges);

            graphpack_graph(
                    vector<Vertex>& U,
                    vector<long long>& vs,
                    vector<int>& es,
                    bool is_U_relabeled); //for P

            graphpack_graph(
                    vector<Vertex>& U,
                    vector<long long>& vs,
                    vector<int>& es,
                    int & mc
                    );

            graphpack_graph(
                    vector<Vertex>& U, vector<long long>& vs, vector<int>& es, int & mc, bool is_tcore);

            graphpack_graph(
                    vector<Vertex>& U, vector<long long>& vs, vector<int>& es, int & mc, bool is_tcore, bool is_enum);

            graphpack_graph(int nverts, int* heads, int* tails);


            /**
             * @brief   Functions for reading "almost" ANY graph given to us.
             *          Meant to be flexible, and assumes the user paramss a
             *          reasonable graph.
             *
             * @param filename
             */
            void read_graph(const string& filename);
            void read_edge_list(const string& filename);
            void detect_delim(string & line, string & delim);
            bool detect_weighted_graph(string & line, string & delim);
            void get_token(int & v, string & line, string & delim, size_t & pos, size_t & prev);
            void get_token(double & weight, string & line, string & delim, size_t & pos, size_t & prev, bool & is_weighted_graph);
            void remove_multiple_edges();

//
//
//            //! @brief given a ptr to a neighbor of a vertex, return the neighbor
//            inline int get_neighbor(long long j) {  return edges[j]; }
//
//            //! @brief given a ptr to a neighbor (of a vertex), return the weight of the neighbor
//            inline int get_weight(long long j) {  return wt[j]; }
////            inline int get_neighbor_weight(long long j) {  return edges[j]; }
//
//


            void clear() {
                vertices.clear();
                edges.clear();
                kcore.clear();
                e_v.clear();
                e_u.clear();
                vertex_lookup.clear();
                kcore_order.clear();
                degree.clear();
                t.clear();
                kappa.clear();
                tri_core.clear();
                tris_ordering.clear();
            }



            //! @brief given a ptr to a neighbor of a vertex, return the original vertex identifier given as params
            inline int get_params_vertex_id(long long j) {  return vertex_lookup[edges[j]]; }
//            inline int get_neighbor_vertex_id(long long j) {  return vertex_lookup[edges[j]]; }

            //! @brief given a vertex id, return the original vertex identifier given as params
            inline int get_params_vertex_id(int v) {  return vertex_lookup[v]; }
//            inline int get_vertex_id(int v) {  return vertex_lookup[v]; }





            void create_adj();
            void reduce_graph(int* &pruned);
            void reduce_graph(
                    vector<long long>& vs,
                    vector<int>& es,
                    int* &pruned,
                    int id,
                    int& mc);

            void reduce_graph(
                    vector<long long>& vs,
                    vector<int>& es,
                    vector<int> & pruned_vertices,
                    graphpack_graph& G,
                    int id,
                    int& mc);

            int num_vertices() { return vertices.size() - 1; }
            int num_edges() { return edges.size()/2; }
            vector <long long>* get_vertices(){ return &vertices; }
            vector<int>* get_edges(){ return &edges; }
            vector<int>* get_degree(){ return &degree; }
//            vector<long long>* get_triangles_array(){ return &tri_core; }
            vector<int> get_edges_array() { return edges; }
            vector<long long> get_vertices_array() { return vertices; };
            // reverse lookup: stores vertex id for edge, used in triangle cores
            vector<long long> e_v, e_u, eid;

            int vertex_degree(int v) { return vertices[v] - vertices[v+1]; }
            long long first_neigh(int v) { return vertices[v]; }
            long long last_neigh(int v) { return vertices[v+1]-1; }

            void sum_vertex_degrees();
            void vertex_degrees();
            void update_degrees();
            void update_degrees(bool flag);
            void update_degrees(int* &pruned, int& mc);
            void update_degrees(vector<long long> & vs, vector<int> & es, int & mc, int* & pruned);
            double density() { return (double)num_edges() / (num_vertices() * (num_vertices() - 1.0) / 2.0); }
            double density(int n, long long m) { return (double) m / (n * (n - 1.0) / 2.0); }
            int get_max_degree() { return max_degree; }
            int get_min_degree() { return min_degree; }
            double get_avg_degree() { return avg_degree; }
            bool is_dense_graph() { return is_dense; }



            void optimize_graph_ops(int adj_limit);

            void initialize();
            string get_file_extension(const string& filename);

            // set/get block size for parallel algorithms
            int get_block_size() { return block_size; }
            void set_block_size(int sz) { block_size = sz; }

            void basic_stats();
            void basic_stats_line();

            void init_edge_data_structs_no_adj();

            void bound_stats(int alg, int lb, graphpack_graph& G);

            // vertex sorter
            template<typename T>
            void compute_ordering(vector<T>& bound, vector<T>& order);

//            template<typename T>
//            void vertex_bucket_sort(vector<T>& bound, vector<T>& order);

            void vertex_bucket_sort(vector<int>& bound, vector<int>& order);

            template<typename VertexType>
            void vertex_bucket_sort(vector<VertexType>& P,  bool smallest_to_largest = true);

            void compute_ordering(string degree, vector<int>& order);


            // edge sorters
            void degree_bucket_sort();
            void degree_bucket_sort(bool desc);
            void degree_bucket_sort_parallel(bool desc);
            void degree_bucket_sort_parallel();

            int max_core;
            vector<int> kcore;
            vector<int> kcore_order;
            vector<int>* get_kcores() { return &kcore; }
            vector<int>* get_kcore_ordering() { return &kcore_order; }
            int get_max_core() { return max_core; }
            void update_kcores(int* &pruned);

            void compute_cores();
            void compute_cores(vector<int> &pruned);
            void induced_cores_ordering(
                    vector<long long>& V,
                    vector<int>& E,
                    int* &pruned);

            void induced_cores_ordering(
                    vector<long long>& V,
                    vector<int>& E,
                    vector<int> & pruned);




            /*----------------
             * triangle cores
             *----------------*/
            vector<long long> tri_core;
            // saves the triangle core ordering for later use
            vector <long long> tris_ordering;
            long long max_tri_core, max_t_edge;
            vector<long long>* get_triangle_cores(){ return &tri_core; }
            long long get_triangle_core_bound(){ return max_tri_core+2; }
//            long long get_max_edge_triangles() { return max_t_edge; }
            void compute_triangle_cores();
            void triangle_cores();
            void count_triangles();
            long long triangle_max_core();
            void print_tcore_info(vector<long long> bin, vector<long long> tris,
                    vector<long long> pos, int mt_end);

            void get_eid_in_triangle_core_order(long long e, long long &eid) { eid = tri_core[tris_ordering[e]];  }

            // takes an edge id, passes v and u back by reference
            void get_vertices_from_eid(const long long eid, long long &v, long long &u) {
                v = e_v[eid];
                u = e_u[eid];
            }

            double triangle_core_time;
            double fast_triangle_core_time;

            double edge_triangles;
            double fast_edge_triangles;

            void triangle_cores_approx(double & triangle_core_perf, double & triangle_perf);
            void triangle_cores_approx_adj(double & triangle_core_perf, double & triangle_perf);

            // parallel codes for triangle cores
            void parallel_triangle_cores();
            void parallel_count_triangles();
            void compute_parallel_triangle_cores();
            void test_total_triangles();
            void compute_max_triangle_core(int lb);

            // reverse lookup
            void map_edges_parallel();
            void map_edges_serial();

            void parallel_edge_kcores();

            /*
             * approximate triangle cores
             * fast for large sparse graphs
             */
            void triangle_cores_approx();
            void triangle_cores_approx(string strategy);
            void parallel_triangle_counting();
            void triangle_core_numbers();
            void triangle_core_numbers_optimized();

            // uses adj, fast check if edge exists, for neighborhood graphs or small graphs
            void triangle_cores_approx_adj();
            void triangle_cores_approx_adj(string strategy);
            void parallel_triangle_counting_adj();
            void triangle_core_numbers_adj();

            // serial triangle counting, but uses adj, mainly for neighborhood triangle cores
            void triangle_counting_adj();
            void triangle_counting_index(int & mc, vector<int> & pruned_edges);

            void triangle_cores_serial_adj();


            void triangle_core_numbers(int & lb);

            void tcore_pruner(graphpack_graph & H, int & lb, vector<int>& U_new);

            graphpack_graph k_triangle_cores(int & k);
            int k_triangle_cores_all(int& k);

            // for graphpack_repair_tcores_mctriangle
            void triangle_cores_serial_adj(int & mc);
//            void triangle_counting_adj(int & mc, vector<int> & pruned_edges);
            void triangle_core_numbers(int & mc, vector<int> & pruned_edges);

            void map_eid_to_vid(bool is_parallel = true);


            /**
             * Triangle Core Decomposition with Edge-CSC Representation
             * fixed the offset
             * see graphpack_graph_ecsc.cpp
             */
            // main routine, calls all other methods
            void triangle_cores_in_parallel(bool is_parallel = true);
            // maps edgeid to vertex pair, setting up edge-csc
            void map_edgeid_to_vertex_pair(bool parallel = true);
            // edge-based triangle counting
            void triangle_counting_in_parallel(bool is_parallel = true);
            void triangle_counting_in_parallel_adj(bool is_parallel = true);
            // main triangle core routine
            void triangle_core_numbers_parallel();



            // uses max to prune edges with insufficient triangle counts
            void k_triangle_counting(int & max, vector<int> & pruned_edges, bool is_parallel = true);


            // generalizing maxclique framework
            void triangle_counting_adj(int & mc, vector<int> & pruned_edges, bool is_parallel = true);
            void triangle_counting(int & mc, vector<int> & pruned_edges, bool is_parallel = true);



            void bin_triangle_core_numbers(); // for distribution
            void triangle_core_numbers_parallel(bool validate);
//            void validate_triangle_cores();
            void validate_triangle_cores(vector<long long> bin, vector<long long> tris, vector<long long> pos, int mt_end);


            /**
             * Densest Subgraph
             */
            void densest_subgraph(string ordering, int & subgraph_idx, double & max_density);
            void triangle_core_numbers_ordering();


            /*-------------------
             * vertex triangles
             *-------------------*/
            double mean_t, tri_ub;
//            uint64_t total_t;
            long long total_t;
            long long max_t;
            vector<long long>   t;
            // see compute_weighted_triangles()
            vector<double> weighted_triangles;
            vector<double> weighted_wedges;
            double total_weighted_triangles;
            double total_weighted_wedges;

            vector<double>  t_ub;
            vector<long long>   triangles;
            void compute_triangles(int block_size = 64);
            void compute_triangles_simple(int block_size = 64);
            // total triangles computed in serial, but counts are done in parallel
            void compute_weighted_triangles(int block_size = 512);
            // computes total triangles in parallel
            void compute_weighted_triangles_reduce(int block_size = 512);

            void triangle_stats();
            void triangle_bound();
            void total_triangles();

//            uint64_t
            long long num_triangles() { return total_t; }
            double get_avg_triangles() { return mean_t; }
            long long get_max_triangles() { return max_t; }
            double get_triangle_bound() { return tri_ub; }

            // weighted triangle defs
            double get_num_weighted_triangles() { return total_weighted_triangles; }
            double get_num_weighted_wedges() { return total_weighted_wedges; }
            double get_avg_weighted_triangles() { return total_weighted_triangles / (double)num_vertices(); }
            double get_avg_weighted_wedges() { return total_weighted_wedges / (double)num_vertices(); }


            long long get_max_edge_triangles() { return max_t_edge; }
            double get_avg_edge_triangles() { return double(max_t_edge) / total_t; }

            // clustering coefficients
            vector<double>  kappa;
            double mean_cc, global_cc;
            double get_cc_avg() { return mean_cc; }
            double get_cc_global() { return global_cc; }



            /*-------------------
             * assortativity
             *-------------------*/
            double r;
            void compute_assort();
            void compute_assort_parallel(int block_size = 64);
            double get_assortativity() { return r; }


            /**
             * Coloring number
             */
            int coloring_number(string vertex_ordering, bool decr_order);
            int color_vertices(vector< Vertex > & V);

            // sequential greedy coloring framework
            vector<int> color;
            vector<int> color_sizes;
            int k_coloring;

            // color_ptr[color] returns ptr to first colored vertex in colored_vertices
            vector<int> color_ptr;
            vector<int> colored_vertices;

            int max_independent_set_size;
            int min_independent_set_size;
            double avg_independent_set_size;


            bool conflicts(vector<vector<int> > & color, int k, int u);
            bool conflicts(vector<vector<int> > & color, int k, vector<short> & ind);




            void triangle_core_pruning(vector<Vertex> & P, int & lb);
            void triangle_core_pruning(vector<Vertex> & P, int & lb, string neigh_ordering);
//            void triangle_core_pruning(vector<Vertex> & P, int & mc);
            void triangle_core_pruning_only_adj(int & mc); /** O(|E|/2) */
            void triangle_core_pruning_only_adj_update_degree(vector<Vertex> & P, int & mc);





            // see clique utils
            void prune_vertex(int* &pruned, int i);

            int initial_pruning(graphpack_graph& G, int* &pruned, int lb);
            int initial_pruning(graphpack_graph& G, int* &pruned, int lb, bool** &adj);

            int initial_pruning(graphpack_graph& G, int* &pruned, int lb, bool** &adj, string pruning_strat, bool is_enum = false);
            int initial_pruning(graphpack_graph& G, int* &pruned, int lb, string pruning_strat, bool is_enum = false);

            int initial_pruning_enum(graphpack_graph& G, int* &pruned, int lb);
            int initial_pruning_enum(graphpack_graph& G, int* &pruned, int lb, bool** &adj);

            // new for tcore, inducing graph..
            int initial_pruning_induced(graphpack_graph& G, int* &pruned, int lb, string pruning_strat);

            void reduce_and_sort_edges(graphpack_graph& G, int* &pruned, int lb_idx);

            void order_vertices(vector<Vertex> &V, graphpack_graph &G,
                    int &lb_idx, int &lb, string vertex_ordering, bool decr_order);

            void order_vertices(vector<Vertex> &V, graphpack_graph &G,
                    string vertex_ordering, bool decr_order);

            void ordering_strategies(int & u, int & val, string & vertex_ordering);


            void order_vertices_tcore(vector<Vertex> &V, graphpack_graph &G,
                    int &lb_idx, int &lb, string vertex_ordering, bool decr_order,
                    vector<int> & pruned_vertices);

            void order_vertices(vector<Vertex> &V, string vertex_ordering, bool decr_order);


            void print_info(vector<int> &C_max, double &sec);
            void print_break();


            bool time_left(vector<int> &C_max, double sec,
                    double time_limit, bool &time_expired_msg);


            void graph_stats(graphpack_graph& G, int& mc, int id, double &sec);


            void graph_stats(graphpack_graph& G, vector<long long> & vs, vector<int> & es,
                    int& mc, int id, double &sec);


            void reduce_graph(
                    vector<long long>& vs,
                    vector<int>& es,
                    int* &pruned,
                    graphpack_graph& G,
                    int id,
                    int& mc);


            void reduce_graph(vector<int> & pruned_vertices, graphpack_graph &G);




            graphpack_graph preprocess_triangle_core(graphpack_graph & G, int & lb);
            void greedy_coloring_perf();

            template<typename T>
            void write_file(std::vector<T> &data, string &filename);

            template<typename T>
            void write_vector(std::vector<T> &data, string suffix);

            /**
             * Some debugging functions
             */
            void print_edges();
            void print_vertices_array();
            void print_base_stats();
            void print_degrees();
            void print_weighted_graph();

            void write_edge_list();

            void init_edge_data_structs();

            void init_edge_data_structs_pruned(int & mc, bool is_enum = false);

            void init_edge_data_structs_pruned(int & mc, vector<Vertex> & U, bool is_enum = false);

            void init_data_structs(int & mc, bool is_enum = false);

            bool clique_test(graphpack_graph& G, vector<int> C);

            void basic_stats(string prefix);

            void test_triangle_cores();
            int test_k_triangle_cores(int & k);

            void write_edge_list(string dir);
            void write_only_edge_list(string dir);

            void write_line(const char *filepath, const char *data);
            string get_filename_from_path(const string& s);

            // mcbound_triangles
            void compute_vertex_triangles(int max = 0, bool is_parallel = true);
            // this function also marks vertex as pruned in the pruned array
            void compute_vertex_triangles(int* &pruned, int max, bool is_parallel, bool is_enum = false);
            // for local neighborhoods
            void compute_vertex_triangles(vector<Vertex> &P, int* &pruned, int max, bool is_parallel);

            void compute_edge_triangles(int max = 0, bool is_parallel = true);

            template<typename T>
            void expansion(std::vector<T> &sets, int num_sets, std::vector<T> &property);

            template<typename T>
            void expansion(std::vector<Vertex> S, std::vector<T> &property);

    };

}
#endif

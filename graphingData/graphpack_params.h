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

#ifndef GRAPHPACK_PARAMS_H_
#define GRAPHPACK_PARAMS_H_

#include "graphpack_headers.h"
#include "graphpack_utils.h"

using namespace std;

class params {
    private:
        bool global_large_to_small;
        bool local_large_to_small;
        bool search_large_to_small;

    public:


        // instance variables
        int algorithm;
        int threads;
        int experiment;
        int lb;
        int ub;
        int param_ub;
        int adj_limit;
        int k;
        double time_limit;
        double remove_time;
        bool graph_stats;
        bool verbose;
        bool help;
        bool MCE;
        bool is_enum;

        string heu_strat;
        int heu_method;
        string format;
        string graph;
        string output;
        string edge_sorter;
        string vertex_search_order;
        string pruning_strat;

        string output_filename;

        string problem; // used for "compression"
        string type;
        string clique_ordering; // used in the compression problem, determines how the cliques are ordered when found

        string graph_reader;


        // adapt algorithms/representation via graph properties
        bool optimizer;

        // --heuristic, --heu = off, disables the fast heuristic
        string heuristic_type;
        string heu_global_ordering;
        string heu_local_ordering;

        bool heu_global_s2l;
        bool heu_local_s2l;

        /**
         * any bounds, different computation costs, depend on application/graph
         * todo: use more than one
         */
        string global_pruning;
        string local_pruning;
        string search_pruning;

        /**
         * any bounds, different computation costs, depend on application/graph
         */
        string global_ordering;
        string local_ordering;
        string search_ordering;

        // edge ordering (any bounds again)
        string global_edge_ordering; // helps cache misses
        bool global_edge_small_to_large;

        // direction of ordering: is smallest to largest?
        bool decreasing_order;
//        bool smallest_to_largest;
        bool global_small_to_large;
        bool local_small_to_large;
        bool search_small_to_large;



        /**
         * adj = adjacency matrix structure
         * csc = compressed sparse column
         * hybrid = both
         */
        string global_representation;
        string local_representation;
        string search_representation;

        string search_method;

        double dynamic_cutoff;
        double density_cutoff;
//        double global_density_cutoff;
//        double local_density_cutoff;

        /**
         *
         */
         bool is_explicit_reduce; // global_explicit_reduce, local_explicit_reduce
         bool is_edges_sorted;

//         bool is_shared; // information shared between threads



        /**
         * @brief Parallelization parameters
         * Type of scheduling: dynamic, static,
         * Block size to use: 1, 64, ...
         *
         * We use a block size of 1 for the majority of maxclique algorithms
         * For computing triangles we use a block size of 64
         */
        string schedule; // dynamic, static, ...
        int block_size;

        bool use_neigh_cores;

        string coloring_method;     // --coloring_method    [color-centric, vertex-centric]
        string coloring_type;       // --coloring_type      [repair, normal]
//        string color_sync;          // no need to synchronize, all workers must share max, since this information allows more pruning, and thus speeding up the algorithm!


        bool approx_tcore;
        bool batch_coloring;


        // fast parallel graph statistics
        struct network_stats {
                bool triangles;
                bool cc;
                bool global_cc;
                bool degree;
                bool density;
                bool assort;
                bool kcore;
                bool triangle_core;
                bool greedy_coloring;
                bool all;
        } network_stats;


        struct graph_type {
                bool weighted;      // is G weighted? 3rd column?
                bool temporal;      // is edge weights time? 3rd col?
                bool directed;      // is directed or undirected?
                bool self_loops;    // allow self-loops
                bool multigraph;    // allow multiple edges
        } graph_type;


        void initialize() {
            // default values
            algorithm = 0;
            threads = omp_get_max_threads();
            experiment = 0;
            // ub = lb = param_ub = k = 0
            ub = 0;
            lb = 0;
            param_ub = 0;
            k = 0; // k-clique problem
            verbose = false;
            graph_stats = false;
            help = false;
            MCE = false;
            is_enum = false;
            decreasing_order = false;
            heu_strat = "kcore";
            heu_method = 0;
            vertex_search_order = "deg";
            pruning_strat = "kcore";
            format = "mtx";
            graph = "graph/test/sample.mtx";
            output = "";
            edge_sorter = "";
//            coloring_method = "repair";
            coloring_method = "color-centric-pruning";
            coloring_type = "repair";

            output_filename = "output.txt";

            graph_reader = "general"; // general = limited assumptions on the graph, slowest, but user avoids having to set flags for defining how the graph is stored

            problem = "mc"; // compress
            type = "exact";
            clique_ordering = "deg"; // kcore, degree, wedges, ...

            approx_tcore = false; // use "-b approx"  (option)

            // parallelization parameters
            schedule = "dynamic"; // dynamic, static, ...
            block_size = 1; // 1 is best for maxclique, a larger values may be slightly faster, but large values also ignore order


            optimizer = true; // true=on, false=off

            // --heuristic, --heu = off, disables the fast heuristic
            heuristic_type = "pruned";  //"pruned", "inplace";//"fast";
            heu_global_ordering = "degree";
            heu_local_ordering = "deg_vol";
            heu_global_s2l = false;
            heu_local_s2l = true;

            /**
             * any bounds, different computation costs, depend on application/graph
             * also, how fast can they be updated and shared? degree is O(1) to update, whereas k-core is much slower
             *
             * todo: use more than one
             */
            global_pruning = "kcore"; // triangles, triangle cores
            local_pruning = "kcore";
            search_pruning = "coloring";

            /**
             * any bounds, different computation costs, depend on application/graph
             */
            global_ordering = "degree";
            local_ordering = "coloring";
            search_ordering = "coloring";

            global_edge_ordering = "degree"; // helps cache misses
            global_edge_small_to_large = false;

            // direction of ordering: is smallest to largest?
//            decreasing_order = false;
//            bool smallest_to_largest;

            global_small_to_large = true;
            local_small_to_large = false;
            search_small_to_large = false;

            global_large_to_small = false;
            local_large_to_small = false;
            search_large_to_small = false;

//            global_small_to_large = false;
//            local_small_to_large = true;
//            search_small_to_large = true;

            /**
             * adj = adjacency matrix structure
             * csc = compressed sparse column
             * hybrid = both
             */
            global_representation = "csc"; // csc useful for extremely large graphs, hybrid = relatively small and/or dense graphs
            local_representation = "hybrid";
            search_representation = "hybrid";

            search_method = "repair";

            dynamic_cutoff = 0.025;
            density_cutoff = 0.80;

            /**
             *
             */

//             is_explicit_reduce = true; // global_explicit_reduce, local_explicit_reduce
//             is_edges_sorted = true;

             // number of vertices before optimizer creates adj
             adj_limit = 10000;
             // terminate after time limit is reached
             time_limit = 60 * 60;           // max time to search
             // a single thread explicitly reduces the graph each 4 seconds
             remove_time = 4.0;              // time to wait before reducing graph
//             max_depth = 12

             use_neigh_cores = false;
             batch_coloring = false;

        }

        void usage(char *argv0) {
            const char *params =

                    "Usage: %s -a alg -f graphfile -t threads -o ordering -h heu_strat -u upper_bound -l lower_bound -r reduce_wait_time -w time_limit \n"
                    "\t-p, --problem                   : Clique variant to solve  \n"
                    "\t \n"
//                    "\t    Clique variants:\n"
//                    "\t        -p mc,  -p maxclique          Maximum Clique Problem:      Find a single maximum clique \n"
//                    "\t        -p mce, -p enum               Maximum Clique Enumeration:  Enumerate all maximum cliques \n"
//                    "\t        -p kclique      -k [num]      K-clique Problem:            Find clique of size k (both -p and -k options must be set) \n"
//                    "\t        -p kclique_enum -k [num]      K-clique Enumeration:        Enumerate all cliques of size k (both -p and -k must be set)\n"

                    "\t Clique variants               \n"
                    "\t========================================================================================================================\n"
                    "\t    Maximum Clique Problem:     -p mc,  -p maxclique       Find a single maximum clique \n"
                    "\t    Maximum Clique Enumeration: -p mce, -p enum            Enumerate all maximum cliques \n"
                    "\t    K-clique Problem:           -p kclique      -k [num]   Find clique of size k (both -p and -k options must be set) \n"
                    "\t    K-clique Enumeration:       -p kclique_enum -k [num]   Enumerate all cliques of size k (both -p and -k must be set)\n"
//                    "\t                                -p kclique -k [num] -enum   \n"
                    "\t\n"
                    "\tNote for the maximum clique problem and maximum clique enumeration, the -k option must NOT be set (or set to -k 0)\n"
                    "\t\n"
                    "\t-a, --algorithm                 : Algorithm for solving MAX-CLIQUE: 0 = full, 1 = no neighborhood cores, 2 = only basic k-core pruning steps  \n"
                    "\t-f, --graph                     : Input GRAPH file for computing the Maximum Clique (matrix market format or simple edge list). \n"
                    "\t-t, --threads                   : Number of THREADS for the algorithm to use (default = max). \n"
                    "\t-b, --block                     : Size of blocks assigned to workers (cores).  (default = 1). \n"
                    "\n"
                    "\n"
                    "\tParameterized Framework for Parallel Maximum Clique Algorithms\n"
                    "\t==============================================================\n"
                    "\n"
                    "\tPRUNING\n"
                    "\t--gp, --global_pruning             : Pruning of G, shared among cores/workers \n"
                    "\t--lp, --local_pruning              : Pruning local non-shared neighborhood pruning\n"
                    "\t--sp, --search_pruning             : Pruning used in the search\n"
                    "\n"
                    "\tORDERING\n"
                    "\t--go, --global_ordering             : Ordering of G, shared among cores/workers \n"
                    "\t--lo, --local_ordering              : Ordering local non-shared neighborhood pruning\n"
                    "\t--so, --search_ordering             : Ordering used in the search\n"
                    "\n"
                    "\tDIRECTION OF ORDERING\n"
                    "\t--gs2l, --global_small_to_large     : Direction of the ordering (default = false, smallest to largest value)  \n"
                    "\t--ls2l, --local_small_to_large      : Direction of the ordering (default = true, smallest to largest value)  \n"
                    "\t--ss2l, --search_small_to_large     : Direction of the ordering (default = true, smallest to largest value)  \n"
                    "\n"
                    "\tREPRESENTATION\n"
                    "\t--gr, --global_representation       : Representation of G, shared among cores/workers \n"
                    "\t--lr, --local_representation        : Representation local non-shared neighborhood pruning\n"
                    "\t--sr, --search_representation       : Representation used in the search\n"
                    "\n"
                    "\t\t'adj'    - adjacency matrix       : dense n by n matrix, O(|V|^2) storage cost\n"
                    "\t\t'csc'    - compressed sparse col  : large sparse graphs, O(|V|+|E|) storage cost\n"
                    "\t\t'hybrid' -  csc + adj             : small sparse and dense graphs, O(|V|^2 + |V| + |E|) storage cost\n"
                    "\n"
                    "\tExample: ./mcpack --gr csc --lr hybrid\n"
                    "\n"

                    "\tDYNAMIC COLORING METHOD FOR SEARCH\n"
                    "\t--dynamic_coloring                  : Coloring method to use in the dynamic search. options: normal or repair \n"

                    "\n"
                    "\tLOWER BOUND\n"
                    "\t--heuristic, --heu                  : Heuristic method (default = fast) \n"
                    "                                      Other choices: off - disables the heuristic \n"
                    "\t--hgo, --heu_go                     : Order of vertices in G \n"
                    "\t--hlo, --heu_lo                     : Local vertex ordering, determines selection of vertices for heuristic \n"
                    "\n"
                    "\tOrder direction\n"
                    "\t--heu_gs2l                          : smallest to largest ordering\n"
                    "\t--heu_ls2l\n"
                    "\n"
                    "\n"
                    "\tSEARCH METHOD\n"
                    "\t--search_method, --sm                  : Search method (default = repair) \n"
                    "\t    repair \n"
                    "\t    static_order \n"
                    "\t    cores \n"

                    "\n\n"
                    "\t COLORING OPTIONS\n"
                    "\t --coloring_method \t[color-centric-pruning (default), color-centric, vertex-centric-pruning, vertex-centric] \n"
                    "\t --coloring_type   \t[repair/recolor (default), off/simple]\n"
                    "\n"

                    "\n"
                    "\tOTHER PARAMTERS\n"
                    "\t--edge_ordering                     : Ordering for the edges in G (default = degree) \n"
                    "\t--density, --density_cutoff         : 0-1 threshold for applying tighter bounds (default = 0.80) \n"
                    "\t--reduce, --reduce_time             : Time to wait before explicitly reducing the graph (default = degree) \n"
                    "\t--time_limit                        : Terminate if runtime exceeds the time limit (default = 1 hour) \n"
                    "\t--adj_limit                         : Number of vertices for which adj representation is used (default = 10000) \n"
                    "\t--dynamic_cutoff                    : Dynamic cutoff of when to apply tighter bounds in the search\n"
                    "\t--ncores, --neigh_cores             : Set flag to use neigh cores for local pruning (default = 0 (off)) \n"


                    "\n"
//                    "\tOTHER PROBLEMS\n"
//                    "\t-k, --k_clique                      : Size of clique to return \n"
//                    "\t--mce, --enumerate                  : Maximum clique enumeration problem \n"


//                    global_edge_ordering = "degree"; // helps cache misses
//
//                    // direction of ordering: is smallest to largest?
//        //            decreasing_order = false;
//        //            bool smallest_to_largest;
//                    global_sm_to_lg = true;
//                    local_sm_to_lg = false;
//                    search_sm_to_lg = true;
//                    density_cutoff = 0.80;
//
//                    /**
//                     *
//                     */
//
//        //             is_explicit_reduce = true; // global_explicit_reduce, local_explicit_reduce
//        //             is_edges_sorted = true;
//                     k = 0; // k-clique problem
//                     // number of vertices before optimizer creates adj
//                     adj_limit = 10000;
//                     // terminate after time limit is reached
//                     time_limit = 60 * 60;           // max time to search
//                     // a single thread explicitly reduces the graph each 4 seconds
//                     remove_time = 4.0;              // time to wait before reducing graph
//        //             max_depth = 12
                    "\n"
                    "\n"
                    "\tUpper Bounds for ordering and pruning strategies\n"
                    "\t==============================================================\n"
                    "\t'degree',   'deg'                    : O(|V|)\n"
                    "\t'kcore',    'kcores'                 : O(|V|)\n"
                    "\t'coloring', 'color'                  : O(|E|)\n"
//                    "\t'wedges',   'wedge'                  : O(|V|)\n"
                    "\t'triangles','tri'                    : O(|E^{3/2}|)\n"
                    "\t'tcore',    'triangle_cores'         : O(|E^{3/2}|)\n"
                    "\n"
                    "\t Other methods for ordering include: \n"
                    "\t'triangle_vol', 'tri_vol'            : O(|E^{3/2}|)\n"
                    "\t'degree_vol',   'deg_vol'            : O(|V|)\n"
                    "\t'kcore_vol',    'kcore_vol'          : O(|E|)\n"
                    "\t'kcore_degree', 'kcore_deg'          : O(|E|)\n"
                    "\t'rand', 'random'                     : O(|V|)\n"
                    "\t'var',  'variance'                   : O(|V|)\n"
                    "\n"
                    "\tExample: ./mcpack --global_ordering deg --global_pruning kcore\n"
                    "\n"
                    "\n"
                    "\t==============================================================\n"
                    //"Usage: %s -a alg -f graphfile -t threads -o ordering -h heu_strat -u upper_bound -l lower_bound -r reduce_wait_time -w time_limit \n"
                    "\t-a algorithm                 : Algorithm for solving MAX-CLIQUE: 0 = full, 1 = no neighborhood cores, 2 = only basic k-core pruning steps  \n"
                    "\t-f graph file                : Input GRAPH file for computing the Maximum Clique (matrix market format or simple edge list). \n"
//                    "\t-o vertex search ordering    : Order in which vertices are searched (default = deg, [kcore, dual_deg, dual_kcore, kcore_deg, rand]) \n"
//                    "\t-d decreasing order          : Search vertices in DECREASING order. Note if '-d' is not set, then vertices are searched in increasing order by default. \n"
//                    "\t-e neigh/edge ordering       : Ordering of neighbors/edges (default = deg, [kcore, dual_deg, dual_kcore, kcore_deg, rand]) \n"
//                    "\t-h heuristic strategy        : Strategy for HEURISTIC method (default = kcore, [deg, dual_deg, dual_kcore, rand, 0 = skip heuristic]) \n"
//                    "\t-m heuristic method          : Heuristic METHOD (default = 0, [1 = fast pmcx-based heuristic, uses ordering from -o and -d ]) \n"
//                    "\t-u upper_bound               : UPPER-BOUND on clique size (default = K-cores).\n"
//                    "\t-l lower_bound               : LOWER-BOUND on clique size (default = Estimate using the Fast Heuristic). \n"
                    "\t-t threads                   : Number of THREADS for the algorithm to use (default = 1). \n"
                    "\t-r reduce_wait               : Number of SECONDS to wait before inducing the graph based on the unpruned vertices (default = 4 seconds). \n"
                    "\t-w time_limit                : Execution TIME LIMIT spent searching for max clique (default = 7 days) \n"
                    "\t-k clique size               : Solve K-CLIQUE problem: find clique of size k if it exists. Parameterized to be fast. \n"
                    "\t-s stats                     : Compute BOUNDS and other fast graph stats \n"
                    /*
                     * computing clique permutation/ordering for applications like graph compression
                     */
                    "\t--problem                    : Select the problem (default=mc, [compress]) \n"
                    "\t--clique_ordering            : Select how the cliques should be ordered for computing the clique ordering (default=deg, kcore, wedges)\n"// [applications = graph compression \n"
                    "\t--type                       : Type of compression: heuristic or exact maximum clique method \n"// [applications = graph compression \n"
                    "\t--output                     : Name of output file \n"
                    "\t-v verbose                   : Output additional details to the screen. \n"
                    "\t-? options                   : Print out this help menu. \n"
                    "\n\n"
                    "\tParallelized graph measures and statistics\n"
                    "\t------------------------------------------\n"
                    "\t--cc clustering coeff        : Compute local clustering coefficient \n"
                    "\t--global_cc                  : Compute global clustering coefficient \n"
                    "\t--triangles                  : Compute the number of triangles in the graph \n"
                    "\t--degree                     : Compute degrees of the graph \n"
                    "\t--kcore                      : Compute kcore decomposition \n"
                    "\t--triangle_core              : Compute triangle core decomposition \n"
                    "\t--greey_coloring             : Compute greedy coloring \n"
                    "\t--all                        : Compute ALL statistics \n"
                    "\n"
                    "Each of these can be outputed to a file by provided --stat_output <filename>\n"

                    "\n\n"
                    "\tGraph Readers: \n"
                    "\t   use '--reader [name]' to select the graph reader, the options for name are given below.\n"
                    "\t------------------------------------------\n"
                    "\tgeneral \t general graph reader makes limited assumptions on params file: \n"
                    "\t\t\tremoves self-loops, duplicate edges, and reads weights if stored (3rd column)\n"
                    "\tmtx (symmetric, sorted edges)\n"
//                    "Bounds and Orderings: \n"
//                    "triangle_core, kcore, degree, \n"
                    "";

            fprintf(stderr, params, argv0);
            exit(-1);
        }

        params() { initialize(); }
        params (int argc, char **argv) {
            if (argc <= 1) usage(argv[0]);
            initialize();

            map<string, int*> int_flags;
            int_flags["-t"] = int_flags["--threads"] = int_flags["--cores"] = &threads;
            int_flags["-l"] = int_flags["--lb"] = int_flags["--lowerbound"] = &lb;
            int_flags["-k"] = int_flags["--k_clique"] = int_flags["--kclique"] = int_flags["-u"] = int_flags["--upperbound"] = int_flags["--ub"] = &k;  // = &param_ub;
            int_flags["-m"] = int_flags["--heu_method"] = &heu_method; // 0 or 1
            int_flags["-b"] = int_flags["--block_size"] = &block_size;
            int_flags["-a"] = int_flags["--alg"] = int_flags["--algorithm"] = &algorithm;
            int_flags["--adj_limit"] = int_flags["--adj_size"] = &adj_limit;

            map<string, bool*> bool_flags;
            bool_flags["-?"] = bool_flags["--help"] = &help;
            bool_flags["-v"] = bool_flags["--verbose"] = &verbose;
            bool_flags["-s"] = bool_flags["--graph_stats"] = bool_flags["--bounds"] = &graph_stats;
            bool_flags["-d"] = bool_flags["--order_direction"] = bool_flags["--direction"] = &decreasing_order;
            bool_flags["-b"] = bool_flags["--approx_tcore"] = &approx_tcore;
            bool_flags["--enum"] = bool_flags["-enum"] = &is_enum;

            bool_flags["--gs2l"] = bool_flags["--global_small_to_large"] = &global_small_to_large;
            bool_flags["--ls2l"] = bool_flags["--local_small_to_large"] = &local_small_to_large;
            bool_flags["--ss2l"] = bool_flags["--search_small_to_large"] = &search_small_to_large;

            bool_flags["--gl2s"] = bool_flags["--global_large_to_small"] = &global_large_to_small;
            bool_flags["--ll2s"] = bool_flags["--local_large_to_small"] = &local_large_to_small;
            bool_flags["--sl2s"] = bool_flags["--search_large_to_small"] = &search_large_to_small;

            bool_flags["--heu_gs2l"] = bool_flags["--heu_global_small_to_large"] = &heu_global_s2l;
            bool_flags["--heu_ls2l"] = bool_flags["--heu_local_small_to_large"] = &heu_local_s2l;
            bool_flags["--ncores"] = bool_flags["--neigh_cores"] = &use_neigh_cores;
            bool_flags["--batch_coloring"] = bool_flags["--batch_colors"] = &batch_coloring;

            map<string, double*> double_flags;
            double_flags["-w"] = double_flags["--time_limit"] = double_flags["--timelimit"] = double_flags["--limit"] = &time_limit;
            double_flags["-r"] = double_flags["--delta"] = double_flags["--remove_time"] = &remove_time;
            double_flags["--threshold"] = double_flags["--thresh"] = double_flags["--density_cutoff"] = double_flags["--density"] = &density_cutoff;
            double_flags["--dynamic_cutoff"] = double_flags["--d_cutoff"] = &dynamic_cutoff;

            map<string, string*> string_flags;
            string_flags["-f"] = string_flags["--graph"] = string_flags["--file"] = &graph;
            string_flags["-h"] = string_flags["--heu_strat"] = string_flags["--heu_order"] = &heu_strat;
            string_flags["-e"] = string_flags["--edge_order"] = string_flags["--edge_ordering"] = &edge_sorter;
            string_flags["-o"] = string_flags["--vertex_ordering"] = string_flags["--vertex_order"] = string_flags["--ordering"] = &vertex_search_order;
            string_flags["--vertex_pruning"] = string_flags["--pruning"] = &pruning_strat;
            string_flags["-p"] = string_flags["--problem"] = string_flags["--prob"] = &problem;
            string_flags["--type"] = &type;
            string_flags["--clique_ordering"] = &clique_ordering; // used in the compression problem, determines how the cliques are ordered when found
            string_flags["--output"] = string_flags["--out"] = string_flags["--output_file"] = &output_filename;
            string_flags["--reader"] = string_flags["--graph_reader"] = string_flags["--graph_format"] = string_flags["--format"] = &graph_reader;

            string_flags["--gp"] = string_flags["--global_pruning"] = &global_pruning;
            string_flags["--lp"] = string_flags["--local_pruning"] = &local_pruning;
            string_flags["--sp"] = string_flags["--search_pruning"] = &search_pruning;

            string_flags["--go"] = string_flags["--global_ordering"] = &global_ordering;
            string_flags["--lo"] = string_flags["--local_ordering"] = &local_ordering;
            string_flags["--so"] = string_flags["--search_ordering"] = &search_ordering;

            string_flags["--gr"] = string_flags["--global_representation"] = &global_representation;
            string_flags["--lr"] = string_flags["--local_representation"] = &local_representation;
            string_flags["--sr"] = string_flags["--search_representation"] = &search_representation;

            string_flags["--heuristic"] = string_flags["--heu"] = &heuristic_type;
            string_flags["--heu_global_ordering"] = string_flags["--hgo"] = string_flags["--heu_go"] = &heu_global_ordering;
            string_flags["--heu_local_ordering"] = string_flags["--hlo"] = string_flags["--heu_lo"] = &heu_local_ordering;
            string_flags["--search_method"] = string_flags["--sm"] = &search_method;
            string_flags["--dynamic_coloring"] = string_flags["--coloring"] =
                    string_flags["--color_method"] = string_flags["--coloring_method"] =
                            string_flags["--cm"] = string_flags["--dc"] = &coloring_method;
            string_flags["--coloring_type"] = string_flags["--color_type"] = &coloring_type;
            string_flags["--edge_go"] = string_flags["--global_edge_ordering"] = string_flags["--edge_ordering"] = string_flags["--eo"] = &global_edge_ordering;


            validate(argc >= 2, "Error: graph must be supplied.");

            for (int i = 1; i < argc; ++i) {
                string flag(argv[i]);
                if (bool_flags.find(flag) != bool_flags.end()) {
                    *bool_flags[flag] = true;
                }
                if(i + 1 < argc) {
                    if (int_flags.find(flag) != int_flags.end()) {
                        *int_flags[flag] = atoi(argv[++i]);
                    }
                    else if (double_flags.find(flag) != double_flags.end()) {
                        *double_flags[flag] = atof(argv[++i]);
                    }
                    else if (string_flags.find(flag) != string_flags.end()) {
                        *string_flags[flag] = string(argv[++i]);
                    }
                    else {
                        printf("Invalid argument %s, ignoring.\n",flag.c_str());
                    }
                }
            }

            if (global_large_to_small) global_small_to_large = false;
            if (local_large_to_small) local_small_to_large = false;
            if (search_large_to_small) search_small_to_large = false;


            if (k != 0) {
                printf("solving k-clique ");
                if (problem == "mce" || problem == "enumeration")   printf("enumeration problem \n");
                else printf("problem \n");
                param_ub = ub = k-1;
            }
            if (algorithm > 9) MCE = true;

//            if (k > 0 && is_enum) {         // k-clique enumeration: if -k AND -enum/--enum flag set
//                problem = "kclique_enum";
//            }
//            else if (k == 0 && is_enum) {   // maxclique enumeration: if --enum flag set, but NOT -k
//                problem = "mce";
//            }
//            else if (k > 0 && !is_enum) { // kclique problem: if k>0, but --enum flag not set
//                problem = "kclique";
//            }
//            else { // maxclique problem
//                problem = "mc";
//            }

            for (int i=1; i < argc; i++) {
                // start of network stats/measures
                if(strcmp(argv[i],"--triangles") == 0)
                    network_stats.triangles = 1;
                else if(strcmp(argv[i],"--cc") == 0)
                    network_stats.cc = 1;
                else if(strcmp(argv[i],"--global_cc") == 0)
                    network_stats.global_cc = 1;
                else if(strcmp(argv[i],"--degree") == 0)
                    network_stats.degree = 1;
                else if(strcmp(argv[i],"--kcore") == 0)
                    network_stats.kcore = 1;
                else if(strcmp(argv[i],"--triangle_core") == 0)
                    network_stats.triangle_core = 1;
                else if(strcmp(argv[i],"--greedy_coloring") == 0)
                    network_stats.greedy_coloring = 1;
                else if(strcmp(argv[i],"--all") == 0) {
                    network_stats.all = 1;
                    network_stats.triangle_core = 1;
                    network_stats.kcore = 1;
                    network_stats.degree = 1;
                    network_stats.global_cc = 1;
                    network_stats.cc = 1;
                    network_stats.triangles = 1;
                    graph_stats = true;
                }
            }

            // both off, use default alg
            if (heu_strat == "0" && algorithm == -1) algorithm = 0;
            if (threads <= 0) threads = 1;

            if (!fexists(graph.c_str())) {
                usage(argv[0]);
                exit(-1);
            }

            FILE* fin = fopen(graph.c_str(), "r+t");
            if (fin == NULL) {
                usage(argv[0]);
                exit(-1);
            }
            fclose(fin);

            cout << "\n\nFile Name ------------------------ " << graph.c_str() << endl;
            if (!fexists(graph.c_str()) ) {
                cout << "File not found!" << endl;
                return;
            }
            omp_set_num_threads(threads);
            printf("[mcpack parallel]  %d worker nodes\n",threads);



            cout <<endl;
            cout << "======================================================================================" <<endl;
            cout << "         Parameterized Parallel Maximum Clique Framework" <<endl;
            cout << "======================================================================================" <<endl;
            cout << "problem = " << problem <<endl;
            print_line(80);
            cout << "heuristic type = " << heuristic_type <<endl;

            if (heu_global_s2l) cout << "heuristic global ordering = " << heu_global_ordering << ", smallest to largest"<<endl;
            else cout << "heuristic global ordering = " << heu_global_ordering << ", largest to smallest"<<endl;

            if (heu_local_s2l) cout << "heuristic local ordering = " << heu_local_ordering << ", smallest to largest" <<endl;// = " << heu_local_s2l <<endl;
            else cout << "heuristic local ordering = " << heu_local_ordering << ", largest to smallest" <<endl; // = " << heu_local_s2l <<endl;

            print_line(80);
            if (global_small_to_large) cout << "global ordering = " << global_ordering << ", smallest to largest" <<endl;
            else cout << "global ordering = " << global_ordering << ", largest to smallest" <<endl;

            cout << "global pruning = " << global_pruning <<endl;
            cout << "global representation = " << global_representation <<endl;

            if (global_edge_small_to_large) cout << "global edge ordering = " << global_edge_ordering << ", smallest to largest = " << global_edge_small_to_large <<endl;
            else cout << "global edge ordering = " << global_edge_ordering << ", largest to smallest" <<endl;

            print_line(80);
            if (local_small_to_large) cout << "local ordering = " << local_ordering << ", smallest to largest " <<endl;
            else cout << "local ordering = " << local_ordering << ", largest to smallest " <<endl;
            cout << "local pruning = " << local_pruning <<endl;
            cout << "local representation = " << local_representation <<endl;

            print_line(80);
            cout << "search method = " << search_method <<endl;
            if (search_small_to_large) cout << "search ordering = " << search_ordering << ", smallest to largest" <<endl;
            else cout << "search ordering = " << search_ordering << ", largest to smallest " <<endl;
            cout << "search pruning = " << search_pruning <<endl;
            cout << "search representation = " << search_representation <<endl;
            print_line(80);
            cout << "coloring_method (dynamic) = " << coloring_method <<endl;
            cout << "coloring_type = " << coloring_type <<endl;
            print_line(80);
            cout << "block size = " << block_size <<endl;
            cout << "workers (processing units, cores) = " << threads <<endl;
            print_line(80);
            cout << "dynamic density cutoff = " << density_cutoff <<endl;
            cout << "dynamic search threshold = " << dynamic_cutoff <<endl;
            cout << "neigh cores (local pruning/ordering) = " << use_neigh_cores <<endl;
            cout << "time limit = " << time_limit <<endl;
            cout << "reduce time = " << remove_time <<endl;
            cout << "adj limit = " << adj_limit <<endl;
            print_line(80);
            cout << "lb = " << lb <<endl;
            cout << "ub = " << ub <<endl;
            cout << "param_ub = " << param_ub <<endl;
            cout << "k = " << k <<endl;
            print_line(80);
            cout << "graph reader: " << graph_reader <<endl;
            print_line(80);
            cout <<endl;
        }
};
#endif

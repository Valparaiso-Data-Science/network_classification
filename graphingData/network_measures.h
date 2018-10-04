#ifndef NETWORK_MEASURES_H_
#define NETWORK_MEASURES_H_

using namespace std;

namespace graphpack {

    static void compute_network_stats(params & p, graphpack_graph & G) {

        string vertex_path = "data/vertex/";
        string vertex_dist_path = "data/dist/"; // distribution/binning

        string delim = ",";
        ostringstream stats_str;
        ostringstream col_names;
        string filename;

        stats_str << G.num_vertices() << delim << G.num_edges() << delim;
        col_names << "|V|" << delim << "|E|" << delim;

        stats_str << G.density() << delim;
        col_names << "density" << delim;

        G.vertex_degrees();

        stats_str << G.get_max_degree() << delim;
        stats_str << G.get_min_degree() << delim;
        stats_str << G.get_avg_degree() << delim;
        col_names << "d_max" << delim;
        col_names << "d_min" << delim;
        col_names << "d_avg" << delim;

        filename = vertex_path + extract_filename(p.graph) + string(".degree");
        write_results(G.degree,filename);


        vector<int> bin_deg(G.max_degree+1,0);
        bin_values(G.degree, bin_deg);
        filename = vertex_dist_path + extract_filename(p.graph) + string(".degree");
        write_results(bin_deg,filename,true);
        bin_deg.clear();

//        G.compute_assort_parallel(); // possible problems with assort > 1.0
        G.compute_assort();
        stats_str << G.get_assortativity() << delim;
        col_names << "r" << delim;


        G.compute_triangles();
        stats_str << G.num_triangles() << delim;
        stats_str << G.get_avg_triangles() << delim;
        stats_str << G.get_max_triangles() << delim;
        stats_str << G.get_cc_avg() << delim;
        stats_str << G.get_cc_global() << delim;
        col_names << "|T|" << delim;
        col_names << "T_avg" << delim;
        col_names << "T_max" << delim;
        col_names << "cc_avg" << delim;
        col_names << "cc_global" << delim;

        filename = vertex_path + extract_filename(p.graph) + string(".triangle");
        write_results(G.t,filename);

        vector<int> bin_tri(G.get_max_triangles()+1,0);
        bin_values(G.t, bin_tri);
        filename = vertex_dist_path + extract_filename(p.graph) + string(".triangle");
        write_results(bin_tri,filename,true);
        bin_tri.clear();



        G.compute_cores();
        stats_str << G.get_max_core()+1 << delim;
        col_names << "K" << delim;

        filename = vertex_path + extract_filename(p.graph) + string(".kcore");
        write_results(G.kcore,filename);

        vector<int> bin_kcore(G.get_max_core()+2,0);
        bin_values(G.kcore, bin_kcore);
        filename = vertex_dist_path + extract_filename(p.graph) + string(".kcore");
        write_results(bin_kcore,filename,true);
        bin_kcore.clear();



//        G.compute_parallel_triangle_cores();
//        stats_str << G.get_triangle_core_bound() << delim;
//        col_names << "T_core" << delim;
//
//        filename = vertex_path + extract_filename(p.graph) + string(".tcore");
//        write_results(G.tri_core,filename);
//
//        vector<int> bin_tcore(G.max_tri_core+1,0);
//        bin_values(G.tri_core, bin_tcore);
//        filename = vertex_dist_path + extract_filename(p.graph) + string(".tcore");
//        write_results(bin_tcore,filename,true);
//        bin_tcore.clear();


        // neighborhood density, triangles, kcore


        p.lb = G.get_max_core()+1; //ensure p.lb is set to the max core lowerbound
        vector<int> C;
//        graphpack_heu maxclique(G,p);
//        p.lb = maxclique.search(G, C, p);

        graphpack_heu maxclique(G,p);
        p.lb = maxclique.search(G, C, p);

//        graphpack_heu maxclique(G,p);
//        p.lb = maxclique.search_for_large_clique(G, C, p);

        stats_str << p.lb << delim;
        col_names << "omega_approx" << delim;

        p.lb = 0;

//        graphpack_mc finder(G,p);
//        int num_colors = finder.greedy_coloring_entire_graph(G,C,p);
//        stats_str << num_colors << delim;
//        col_names << "approx_colors" << delim;


//        filename = vertex_path + extract_filename(p.graph) + string(".coloring");
//        write_results(G.color,filename);
//
//
//        vector<int> bin_colors(num_colors+1,0);
//        bin_values(G.color, bin_colors);
//        filename = vertex_dist_path + extract_filename(p.graph) + string(".coloring");
//        write_results(bin_colors,filename,true);
//        bin_colors.clear();


//
//        int num_colors_dist_two = 0;
//        compute_greedy_distance_two_coloring(G, G.vertices, G.edges, num_colors_dist_two,p);
//        stats_str << num_colors_dist_two << delim;
//        col_names << "approx_dist2_colors";// << delim;



        // TODO:
        // For KCORE/TRICORE: num edges/verts in largest k-core, density, kcore with largest vertices
        // Denest subgraph
        // For coloring/dist2-coloring: largest indep set size, avg indep set size


        vector<string> stats_vec;
        stats_vec.push_back(col_names.str());
        stats_vec.push_back(stats_str.str());

        // p.graph contains fullpath to save global stats
        write_results(stats_vec,p.graph);
    }

    static void compute_graph_parameters(params & params, graphpack_graph & G) {

        double sec = get_time();

        /**
         * Pattern for processing output is simply:
         * "prefix measure name" followed by ":" and the value of the measure
         */
        string prefix = "network measures> ";
//        string stats = "";
//        string labels = "";
        ostringstream convert;
        ostringstream lbl;

//        stats += string(params.graph) + string(" \t");



        cout << prefix << "|V|: " << G.num_vertices() <<endl;
        cout << prefix << "|E|: " << G.num_edges() <<endl;


//        stats  += G.num_vertices() + "\t" + G.num_edges() + "\t";
//        labels += "|V| \t |E| \t";

        convert << G.num_vertices() << "\t" << G.num_edges() << "\t";
        lbl << "|V| \t |E| \t";

        /**
         * Density of the graph
         */
//        if (params.network_stats.density == 1 || params.network_stats.all == 1) {
            cout << prefix << "p: " << G.density() <<endl;

            convert << G.density() << "\t";
            lbl << "p \t";
//        }



        /**
         * vertex degrees and avg degree
         */
//        if (params.network_stats.degree == 1 || params.network_stats.all == 1) {
            G.vertex_degrees();
            cout << prefix << "d_max: " << G.get_max_degree() <<endl;
            cout << prefix << "d_avg: " << G.get_avg_degree() <<endl;

            convert << G.get_max_degree() << "\t" << G.get_avg_degree() << "\t";
            lbl << "d_max \t d_avg \t";
//        }

        /**
         * Assortative
         */
//        if (params.network_stats.assort == 1 || params.network_stats.all == 1) {
//            G.compute_assort();
            G.compute_assort_parallel();
//            stats += G.get_assortativity() + "\t";
//            labels += "r \t";

            cout << prefix << "r: " << G.get_assortativity() <<endl;
            convert << G.get_assortativity() << "\t";
            lbl << "r \t";
//        }

        /**
         * Triangles, local clustering coefficient,
         * global clustering coefficient, ...
         */
//        if (params.network_stats.triangles == 1 ||
//                params.network_stats.cc == 1 ||
//                params.network_stats.global_cc == 1 || params.network_stats.all == 1) {
            G.compute_triangles();
            cout << prefix << "|T|: " << G.num_triangles() <<endl;
            cout << prefix << "T_avg: " << G.get_avg_triangles() <<endl;
            cout << prefix << "T_max: " << G.get_max_triangles() <<endl;
            cout << prefix << "cc_avg: " << G.get_cc_avg() <<endl;
            cout << prefix << "cc_global: " << G.get_cc_global() <<endl;
//            stats += G.num_triangles() + "\t" + G.get_avg_triangles() + "\t" + G.get_max_triangles() + "\t" + G.get_cc_avg() + "\t" + G.get_cc_global();
//            labels += "|T| \t T_avg \t T_max \t cc_avg \t cc_global \t";

            convert << G.num_triangles() << "\t" << G.get_avg_triangles();
            convert << "\t" << G.get_max_triangles() << "\t" << G.get_cc_avg() << "\t" << G.get_cc_global() << "\t";
            lbl << "|T| \t T_avg \t T_max \t cc_avg \t cc_global \t";
//        }




        /**
         * K-core decomposition
         */
//        if (params.network_stats.kcore == 1 || params.network_stats.all == 1) {
            G.compute_cores();
            cout << prefix << "K: " << G.get_max_core()+1 <<endl;
//            stats += G.get_max_core()+1 + "\t";
//            labels += "K \t";

            convert << G.get_max_core()+1 << "\t";
            lbl << "K \t";
//        }

        /**
         * Triangle core decomposition
         */
//        if (params.network_stats.triangle_core == 1 || params.network_stats.all == 1) {
//            G.compute_triangle_cores();
            G.compute_parallel_triangle_cores();
            cout << prefix << "T_core: " << G.get_triangle_core_bound() <<endl;
//            stats += G.get_triangle_core_bound() + "\t";
//            labels += "T_core \t";

            convert << G.get_triangle_core_bound() << "\t";
            lbl << "T_core \t";
//        }


//
        //! lower-bound of max clique
        vector<int> C;
        if (params.heu_method == 0 || params.network_stats.all == 1) {
            graphpack_heu maxclique(G,params);
            params.lb = maxclique.search(G, C, params);
        }
        cout << "Heuristic found clique of size " << params.lb;
        cout << " in " << get_time() - sec << " seconds" <<endl;
        cout << "[mcpack: heuristic]  ";
        print_max_clique(C);

        convert << params.lb << "\t";
        lbl << "omega_approx \t";
//




        /**
         * Perform a greedy coloring of the graph
         * Coloring number approximation
         */
//        if (params.network_stats.greedy_coloring == 1) {
//            G.coloring_number("degree","desc");
//        }


        /**
         * Performance in terms of runtime seconds
         */
        cout << prefix << "time: " << get_time() - sec <<endl;

        vector<string> stats_vec, labels_vec;
        stats_vec.reserve(10);
        labels_vec.reserve(10);

        stats_vec.push_back(convert.str());
        labels_vec.push_back(lbl.str());

        string fn = string("stats-") + string(params.graph);
        string fn_labels = string("labels-stats-") + string(params.graph);
        cout << "writing graph stats to " << fn <<endl;
        write_results(stats_vec,fn);
        write_results(labels_vec,fn_labels);

    }

    static void compute_network_measures(params & params, graphpack_graph & G) {

        double sec = get_time();

        /**
         * Pattern for processing output is simply:
         * "prefix measure name" followed by ":" and the value of the measure
         */
        string prefix = "network measures> ";
//        string stats = "";
//        string labels = "";
//        ostringstream Convert;

        cout << prefix << "|V|: " << G.num_vertices() <<endl;
        cout << prefix << "|E|: " << G.num_edges() <<endl;

//        stats += string(params.threads) + "\t" + string(G.num_vertices()) + "\t" + string(G.num_edges()) + "\t";
//        labels += "threads \t |V| \t |E| \t";

        /**
         * Density of the graph
         */
        if (params.network_stats.density == 1) {
            cout << prefix << "p: " << G.density() <<endl;
//            stats += G.density() + "\t";
//            labels += "p \t";
        }



        /**
         * vertex degrees and avg degree
         */
        if (params.network_stats.degree == 1) {
            G.vertex_degrees();
            cout << prefix << "d_max: " << G.get_max_degree() <<endl;
            cout << prefix << "d_avg: " << G.get_avg_degree() <<endl;
//            stats += G.get_max_degree() + "\t" + G.get_avg_degree() + "\t";
//            labels += "d_max \t d_avg \t";
        }

        /**
         * Assortative
         */
        if (params.network_stats.assort == 1) {
            G.compute_assort();
//            stats += G.get_assortativity() + "\t";
//            labels += "r \t";
        }

        /**
         * Triangles, local clustering coefficient,
         * global clustering coefficient, ...
         */
        if (params.network_stats.triangles == 1 ||
                params.network_stats.cc == 1 ||
                params.network_stats.global_cc == 1) {
            G.compute_triangles();
            cout << prefix << "|T|: " << G.num_triangles() <<endl;
            cout << prefix << "T_avg: " << G.get_avg_triangles() <<endl;
            cout << prefix << "T_max: " << G.get_max_triangles() <<endl;
            cout << prefix << "cc_avg: " << G.get_cc_avg() <<endl;
            cout << prefix << "cc_global: " << G.get_cc_global() <<endl;
//            stats += G.num_triangles() + "\t" + G.get_avg_triangles() + "\t" + G.get_max_triangles() + "\t" + G.get_cc_avg() + "\t" + G.get_cc_global();
//            labels += "|T| \t T_avg \t T_max \t cc_avg \t cc_global \t";
        }




        /**
         * K-core decomposition
         */
        if (params.network_stats.kcore == 1) {
            G.compute_cores();
            cout << prefix << "K: " << G.get_max_core()+1 <<endl;
//            stats += G.get_max_core()+1 + "\t";
//            labels += "K \t";
        }

        /**
         * Triangle core decomposition
         */
        if (params.network_stats.triangle_core == 1) {
            G.compute_triangle_cores();
            cout << prefix << "T_core: " << G.get_triangle_core_bound() <<endl;
//            stats += G.get_triangle_core_bound() + "\t";
//            labels += "T_core \t";
        }




        /**
         * Perform a greedy coloring of the graph
         * Coloring number approximation
         */
//        if (params.network_stats.greedy_coloring == 1) {
//            G.coloring_number("degree","desc");
//        }


        /**
         * Performance in terms of runtime seconds
         */
        cout << prefix << "time: " << get_time() - sec <<endl;

    }
}
#endif /* NETWORK_MEASURES_H_ */

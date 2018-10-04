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

#include "graphpack.h"
#include "network_measures.h" // for NR

using namespace std;
using namespace graphpack;

/**
 * NR STATS: USE THIS FUNCTION FOR COMPUTING STATS OF LARGER GRAPHS (WHICH ARE UNABLE TO BE COMPUTED VIA WEB)
 * NOTE: ENSURE ~/Desktop/NetworkRepository/graphpack_with_nr/ is USED
 *
 * STAT FILES SAVED TO DATA DIR {GLOBAL, DIST, VERTEX, EDGE}, e.g., "data/global/" should contain ".txt" files
 *
 * EXAMPLE:
 * ./graphpack -p stats -f graphs/ca-netscience.mtx
 *
 * @param p
 */
static void batch_network_stats(params & p) {
	cout << "\nBATCH NETWORK STATS\n" <<endl;
	vector<string> dirs;

	//    string path = "/Volumes/Macintosh HD/Users/ryarossi/graphs_proc/";
	string path = "/Users/rrossi/Desktop/NetworkRepository/data/graphs/";
//	string path = "../"; // apollo -> for bn

	//    dirs.push_back(path + "");

//	dirs.push_back(path + "socfb/");
//	dirs.push_back(path + "bn/");

//
	dirs.push_back(path + "misc/");
//	dirs.push_back(path + "eco/");
	dirs.push_back(path + "web/");
	dirs.push_back(path + "dimacs10/");
//
////
////	//    string path = "../../Desktop/BigGraphAnalytics/data/graphs/";
////	//    dirs.push_back(path + "bio/");
////	//    dirs.push_back(path + "ca/");
	    dirs.push_back(path + "dynamic_graphs/");
//	    dirs.push_back(path + "ia/");
//	    dirs.push_back(path + "inf/");
////	//    dirs.push_back(path + "misc/");
//	    dirs.push_back(path + "rec/");
	    dirs.push_back(path + "retweet_graphs/");
	    dirs.push_back(path + "rt/");
	    dirs.push_back(path + "sc/");
	    dirs.push_back(path + "tech/"); // tech-p2p not finished... needs to be fixed! reprocessed.
	    dirs.push_back(path + "soc/");
////	//    dirs.push_back(path + "socfb/");
////	//    dirs.push_back(path + "misc/");
//	    dirs.push_back(path + "tscc/");
//	    dirs.push_back(path + "web/");
////	    dirs.push_back(path + "rand/");
////	    dirs.push_back(path + "massive/");
////	    dirs.push_back(path + "new_networks/");
//	    dirs.push_back(path + "dimacs/");
//	    dirs.push_back(path + "bhoslib/");
//	//    dirs.push_back(path + "dimacs10/");


	string global_path = "data/global/";
	string vertex_path = "data/vertex/";
	string vertex_dist_path = "data/dist/"; /// distribution/binning
	string edge_path = "data/edge/";

	for (int d = 0; d < dirs.size(); ++d) {
//	for (int d=dirs.size()-1; d>=0; --d) {
		string dir = dirs[d];
		vector<string> files;

		/// returns the files in that directory
		getdir (dir,files);

		for (int i = 0; i < files.size(); ++i) {
			cout << files[i] << endl;
			if (files[i] == "." || files[i] == ".." || files[i] == ".DS_Store" || files[i] == ".DS_Store.dimacs" || files[i] == "keller6-FAIL.mtx") { continue; }
			if (files[i]=="discogs_style.edges" || files[i]=="discogs_lgenre.edges" || files[i]=="discogs_genre.edges") { continue; }
			//! read graph
			string fn = string(dir) + string(files[i]);
			cout << "                 Input graph: " << get_filesize(fn.c_str()) << " bytes" <<endl;

			int size_ub = 590000000; // 590mb
			if (get_filesize(fn.c_str()) > size_ub) continue; // ignore file if larger than SIZE_UB

			graphpack_graph G(p.graph_stats,fn);
			G.verbose = p.verbose;

			if (G.vertices.size() == 0 || G.edges.size() == 0) continue;
			G.optimize_graph_ops(p.adj_limit);

			string filename = global_path + extract_filename(string(files[i])) + string(".txt");
			p.graph = filename;

			/// NOTE: Final global stats path - speedup stats computations by avoiding graphs that already have stats
//			string nr_path = "~/Desktop/NetworkRepository/data/network_info/global/";
			string nr_path = "./data/global/";
			string global_stats_fn = nr_path + extract_filename(string(files[i])) + string(".txt");
			if (fexists(global_stats_fn.c_str())==false) { /// graph does not have global stats, compute them!!
				/// compute global parameters
				compute_network_stats(p,G);
			}
		}
	}

}


int main(int argc, char *argv[]) {

    //! parse command args
    params p(argc, argv);

    //! read graph
    graphpack_graph G(p.graph_stats,p.graph);
    G.verbose = p.verbose;

    // Automatically choose the optimal graph representation
    G.optimize_graph_ops(p.adj_limit);

    if (p.graph_stats) { G.bound_stats(p.algorithm, p.lb, G); }

    //! ensure wait time is greater than the time to recompute the graph data structures
    if (G.num_edges() > 1000000000 && p.remove_time < 120)  p.remove_time = 120;
    else if (G.num_edges() > 250000000 && p.remove_time < 10) p.remove_time = 10;
    cout << "explicit reduce is set to " << p.remove_time << " seconds" <<endl;

    //! upper-bound of max clique
    double sec = tic();
    G.compute_cores();
    p.ub = G.get_max_core() + 1;
    cout << "K: " << p.ub <<endl;
    cout << "K_time: " << get_time() - sec <<endl;


    if (p.problem == "stats") {
	    cout << "calling batch_network_stats(p)" <<endl;
	    batch_network_stats(p);
	    print_line(100); cout << "FINISHED COMPUTING STATS... terminating..." <<endl; print_line(100);
	    return 0;
    }

    //! lower-bound of max clique
    vector<int> C; int heu_clique_size = 0;
    if (p.heuristic_type != "off") {
        if (p.heuristic_type == "inplace" || p.heuristic_type == "pruned") {
            graphpack_heu maxclique(G,p);
            p.lb = maxclique.search_for_large_clique(G, C, p);
        }
        else if (p.heu_method == 0) {
            graphpack_heu maxclique(G,p);
            p.lb = maxclique.search(G, C, p);
            heu_clique_size = p.lb;
        }
        cout << "Heuristic found clique of size " << p.lb;
        cout << " in " << get_time() - sec << " seconds" <<endl;
        cout << "[graphpack: heuristic]  ";
        print_max_clique(C);
        if (!G.clique_test(G,C))  { print_line(80); cout << "ERROR: C is NOT a clique." <<endl; print_line(80); }
    }

    if (p.lb == 0) { if (p.k >= p.ub) { p.lb = p.ub; } else { p.lb = p.k-1; } }

    /**
     * Maximum Clique Problem
     * K-clique Problem
     */
    if (p.problem == "maxclique" || p.problem == "mc" ||
            ((p.problem == "kclique" || p.problem == "k-clique") && p.k >= 0)) {

        //! check solution found by heuristic
        if (p.lb == p.ub && p.k == 0)
            cout << "Heuristic found optimal solution." << endl;
        else {
            // k-clique problem (k > 0) and lb larger than k
            if (p.k > 0 && p.lb > p.k) {    p.lb = p.k - 1; }
            if (p.k == 0) { print_line(80); printf("MAXIMUM CLIQUE PROBLEM \n");  print_line(80); }
            else          { print_line(80); printf("K-CLIQUE PROBLEM \n");        print_line(80); }

            //! find a single maximum clique
            graphpack_mc finder(G,p);
            p.lb = finder.maxclique(G,C,p);

            toc(sec);
            cout << "Time taken: " << sec << " SEC" << endl;
            cout << "Size (omega): " << p.lb << endl;
            print_max_clique(C);

            // if heu found optimal, then we already verified the clique (note: maxclique may also explicitly remove verts)
            if (p.lb != heu_clique_size && !G.clique_test(G,C)) { print_line(80); cout << "ERROR: C is NOT a clique." <<endl; print_line(80); }
        }
    }
    /**
     * Maximum Clique Enumeration
     * K-clique Enumeration
     */
    else if (p.problem == "mce" || p.problem == "enumeration" || p.problem == "enum" ||
            p.problem == "kclique_enum" || p.problem == "k-clique_enum") {
        set< vector<int> > cliques;
        p.is_enum = true;
        if (p.k == 0) { // MAX CLIQUE ENUMERATION
            print_line(80); printf("MAXIMUM CLIQUE ENUMERATION PROBLEM \n"); print_line(80);

            //! now enumerate those cliques
            graphpack_mc enumerator(G,p);
            p.lb = enumerator.search_enum(G,cliques,p);
        }
        else { // K-CLIQUE ENUMERATION
            print_line(80); printf("K-CLIQUE ENUMERATION PROBLEM \n"); print_line(80);

            // ensure lb is no larger than k, otherwise cliques of size k will be missed.
            if (p.lb > p.k) p.lb = p.k;
            graphpack_mc enumerator(G,p);
            p.lb = enumerator.search_enum(G,cliques,p);
        }
        toc(sec);
        cout << "Time MCE: " << sec << " SEC" << endl;
        cout << "Size MCE (omega): " << p.lb << endl;
        cout << "Number of maximum cliques (mu): " << cliques.size() <<endl;

        check_maxcliques(G,cliques,p.lb);
        if (p.verbose) { print_n_maxcliques(cliques, cliques.size()); }
        else print_n_maxcliques(cliques, 1);
    }
    return 0;
}

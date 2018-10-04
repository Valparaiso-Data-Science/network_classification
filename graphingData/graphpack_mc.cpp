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

#include "graphpack_mc.h"

using namespace std;
using namespace graphpack;

/**
 * maximum clique problem and k-clique problem
 */
int graphpack_mc::maxclique(graphpack_graph& G, vector<int>& sol, params & p) {

    vertices = G.get_vertices();
    edges = G.get_edges();
    degree = G.get_degree();

//    if (p.problem == "mce") lb = lb - 1;
//    if (p.problem == "kclique") { lb=lb-1; cout << "[maxclique]  k-clique problem, using lower bound " <<  lb <<endl; }
    if (p.k>0) { lb=lb-1; cout << "[maxclique]  k-clique problem, using lower bound " <<  lb <<endl; }

    int* pruned = new int[G.num_vertices()];
    memset(pruned, 0, G.num_vertices() * sizeof(int));
    int max = lb;
    int u = 0;

    // initial pruning
    int lb_idx = 0;
    switch (global_pruning) {
        case KCORE: {
            printf("pruning method: k-core \n");
            lb_idx = pruning(G, max, p.global_pruning, pruned, false, is_enum);
            break;
        }
        case TRIANGLES: {
            // selects between parallel and sequential or csc only vs. hybrid of csc/adj
            printf("pruning method: k-core and triangles \n");

            // first prune using k-cores
            lb_idx = pruning(G, max, "kcore", pruned, false, is_enum);

            // then prune via triangles
            G.compute_vertex_triangles(pruned, max, true);
            break;
        }
        case TRIANGLES_ONLY: {
            printf("pruning method: triangles only \n");
            // selects between parallel and sequential (or csc only vs. hybrid of csc/adj)
            G.compute_vertex_triangles(pruned, max, true, is_enum);
            break;
        }
        case TRIANGLE_CORES: {

             vector<int> W;
             for (int i = G.num_vertices()-1; i >= 0; i--) {
                 if (G.kcore[G.kcore_order[i]] >= max)  W.push_back(G.kcore_order[i]);
             }
             graphpack_graph H(W,G.vertices,G.edges);
             H.basic_stats();

             H.map_edgeid_to_vertex_pair(true);

             vector<int> pruned_edges(G.num_edges()+1,0);
             H.k_triangle_counting(max,pruned_edges, true);
             printf("[triangle cores]  argmax tr(u) = %lld,  sum_u tr(u) = %llu \n", H.max_t_edge, H.total_t);
             H.triangle_core_numbers_parallel();

             /**
              * this function explicitly removes edges, but only marks vertices as pruned
              */
             triangle_core_remove_edges(H,pruned,pruned_edges,max);
             H.compute_cores();
             printf("updated max core = %d \n", H.max_core);
             G = H;
             break;
        }
    } // end of global_pruning


    cout << "[graphpack: initial pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
    G.reduce_graph(pruned);
    G.update_degrees();
    cout << "[graphpack: initial pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;

    if (p.global_edge_ordering == "kcore") {
        printf("ordering edges by k-cores \n");
        edge_bucket_sort_parallel(G, G.kcore, p);
    }
    else {
        printf("ordering edges by degree (default) \n");
        edge_bucket_sort_parallel(G, G.degree, p);
    }

    bool** adj;
    if (G.is_dense_graph()) {
        adj = G.adj;
        cout << "[mcpack] the graph G is dense! " <<endl;
    }
    else    cout << "[mcpack] the graph G is sparse! " <<endl;


    // set to worst case bound of cores/coloring
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);

    vector<int> C, C_max;
    C.reserve(G.get_max_degree()+1);
    C_max.reserve(G.get_max_degree()+1);

    // init the neigh coloring array
    vector< vector<int> > colors(G.get_max_degree()+1);
    for (int i = 0; i < G.get_max_degree()+1; i++)  colors[i].reserve(G.get_max_degree()+1);

    // order verts for our search routine
    vector<Vertex> V;
    V.reserve(G.num_vertices());
    compute_ordering(V,G,global_ordering,p.global_small_to_large,lb_idx,lb);
    cout << "|V| = " << V.size() << ", f(0) = " << V[0].get_bound() << ", f(|V|) = " << V[V.size()-1].get_bound() <<endl;

    vector<short> ind(G.num_vertices(),0);
    vector<int> es = G.get_edges_array();
    vector<long long> vs = G.get_vertices_array();

    vector<double> induce_time(num_threads,get_time());
    for (int t = 0; t < num_threads; ++t)  induce_time[t] = induce_time[t] + t; // t/4;


    vector< vector<OpCount> > S; // setup dynamic search data structures
    S.resize(num_threads+1);
    if (search_method == DYNAMIC || search_method == DYNAMIC_KCORE) {
        printf("setting up dynamic search data structures \n");
        for (int t = 0; t < num_threads; ++t) S[t].resize(G.get_max_core() + 1);
        cout << "setting up dynamic step counting structs" <<endl;
        for (int t = 0; t < num_threads; ++t)
            for (int i=0; i < G.get_max_core() + 1; ++i) {
                S[t][i].set_i1(0);
                S[t][i].set_i2(0);
            }
    }
    cout << "num threads = " << num_threads <<endl;
    vector<int> levels(num_threads+1,0);
    vector<int> steps(num_threads+1,0);

    #pragma omp parallel for schedule(dynamic,block_size) shared(pruned, G, T, V, max, C_max, induce_time) \
        firstprivate(colors,ind,vs,es) private(u, P, C)
    for (int i = 0; i < (V.size()) - (max-1); ++i) {
        if (not_reached_ub) {
            if (G.time_left(C_max,sec,time_limit,time_expired_msg)) {

                u = V[i].get_id();
                if ((*bound)[u] > max) {
                    P.push_back(V[i]);
                    for (long long j = vs[u]; j < vs[u + 1]; ++j)
                        if (!pruned[es[j]])
                            if ((*bound)[es[j]] > max) {
                                P.push_back(Vertex(es[j], (*degree)[es[j]]));
                            }

                    switch(local_representation) {
                        case HYBRID: {

                            if (P.size() > max) { // degree bound
                                // neighborhood core ordering and pruning
                                if (use_neigh_cores) {
                                    neigh_cores_bound(vs,es,P,ind,max);

                                    if (! (P.size() > max && P[0].get_bound() >= max) ) {
                                        printf("BREAKING DUE TO _LOCAL_ NEIGHBORHOOD CORES \n");
                                        break;
                                    }
                                }

                                if (P.size() > max) {

                                    graphpack_graph N(P,vs,es,max,p.local_pruning == "tcores"); // true if tcores, otherwise false

                                    if (P.size() > max) { // if (P.size() <= max) continue;

                                        if (N.density() > p.density_cutoff) { // for flickr and other dense graphs
                                            if (G.verbose) printf("neighborhood density = %lg \n", N.density());

                                            switch(local_pruning) {
                                                case TRIANGLES: {
                                                    N.compute_vertex_triangles(P,pruned,max,false);
                                                    break;
                                                }
                                                case TRIANGLE_CORES: {
                                                    N.triangle_cores_serial_adj(max);

                                                    if (N.get_triangle_core_bound() > max) {

                                                        // only prune if neighborhood triangle-core bound passes
                                                        N.triangle_core_pruning(P,max);
                                                        if (G.verbose) printf("p = %lg, T_v = %lld, |P| = %lu \n", N.density(), N.get_triangle_core_bound(), P.size());
                                                        // only sort if P has enough vertices after pruning via triangle-cores
                                                        if (P.size() <= max) { P = T; }
                                                    } // bound
                                                    else { // set P to be empty, so the check below fails
                                                        P = T;
                                                    }
                                                    break;
                                                }
                                            }
                                        }

                                        if (P.size() > max) {
                                            // note: vertices in P are mapped using old vertex identifiers!
                                            greedy_neighborhood_coloring(N.vertices, N.edges, P, C, C_max, colors, max, N.adj);

                                            if (P.back().get_bound() > max) {

                                                // todo: fix problem
                                                if (local_ordering != COLORING) {
                                                    printf("local ordering = %d \n", local_ordering);
                                                    compute_local_ordering(P,N.vertices,N.edges,local_ordering,p.local_small_to_large,0,max);
                                                }

                                                switch(search_method) {
                                                    case REPAIR: { // repair + static ordering
                                                        vector<Vertex> ColOrd = P; // copy P
                                                        branch_neighborhood_repair(N, P, ColOrd, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                    case REPAIR_ONLY: {
                                                        branch_neighborhood_repair_only(N, P, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                    case STATIC_ORDER: { // static ordering only
                                                        vector<Vertex> ColOrd = P;
                                                        branch_neighborhood_static_order(N, P, ColOrd, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                    case KCORE_SEARCH: {
                                                        // kcore pruning + greedy coloring
                                                        branch_neighborhood_cores(N, P, C, C_max, ind, colors, pruned, max);
                                                        break;
                                                    }
                                                    case DYNAMIC_KCORE: {
                                                        // kcore pruning + greedy coloring
                                                        branch_dynamic_neighborhood_cores(N, P, C, C_max, ind, colors, pruned, max, S, levels, steps);
                                                        break;
                                                    }
                                                    case DYNAMIC: {
                                                        // kcore pruning + greedy coloring
                                                        branch_dynamic_neighborhood(N, P, C, C_max, ind, colors, pruned, max, S, levels, steps);
                                                        break;
                                                    }
                                                    case DYNAMIC_COLORING: {
                                                        // only greedy coloring
                                                        branch_neighborhood(N, P, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                    case BASIC: {
                                                        branch_basic_neighborhood(N, P, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                    default: {
                                                        branch_basic_neighborhood(N, P, C, C_max, colors, pruned, max);
                                                        break;
                                                    }
                                                } // select search method
                                            }
                                        }
                                    }
                                }
                            }
                            break;
                        } //end of HYBRID case
                        case CSC: {
                            if (P.size() > max) { // degree bound
                                // neighborhood core ordering and pruning
                                if (use_neigh_cores) {
                                    neigh_cores_bound(vs,es,P,ind,max);
                                    if (! (P.size() > max && P[0].get_bound() >= max) ) { break; }
                                }

                                if (P.size() > max) {

                                    if (G.is_dense_graph())
                                        neigh_coloring_dense(vs,es,P,ind,C,C_max,colors,max, adj);
                                    else neigh_coloring_bound(vs,es,P,ind,C,C_max,colors,pruned,max);

                                    if (P.back().get_bound() > max) {

                                        switch(local_pruning) {
                                            case TRIANGLES: {
                                                G.compute_vertex_triangles(P,pruned,max,false);
                                                break;
                                            }
                                        }

                                        if (local_ordering != COLORING)
                                            compute_local_ordering(P,G,vs,es,local_ordering,p.local_small_to_large,0,max);


                                        if (G.is_dense_graph()) {
                                            switch(search_method) {
                                                case REPAIR: { // repair + static ordering
                                                    vector<Vertex> ColOrd = P; // copy P
                                                    branch_dense_repair(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max, adj);
                                                    break;
                                                }
                                                case REPAIR_ONLY: {
                                                    // dynamic coloring and repair with dynamic ordering
                                                    branch_dense_repair_only(vs,es, P, ind, C, C_max, colors, pruned, max, adj);
                                                    break;
                                                }
                                                case STATIC_ORDER: { // static ordering only
                                                    vector<Vertex> ColOrd = P;
                                                    branch_dense_static_order(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max, adj);
                                                    break;
                                                }
                                                case DYNAMIC_KCORE: {
                                                    // kcore pruning + greedy coloring
                                                    break;
                                                }
                                                case DYNAMIC: {
                                                    // kcore pruning + greedy coloring
                                                    branch_dense_dynamic(vs, es, P, C, C_max, ind, colors, pruned, max, adj, S, levels, steps);
                                                    break;
                                                }
                                                case DYNAMIC_COLORING: {
                                                    branch_dense(vs,es,P, ind, C, C_max, colors, pruned, max, adj);
                                                    break;
                                                }
                                                default: {
                                                    // only greedy coloring
                                                    branch_basic_dense(P, ind, C, C_max, pruned, max, adj);
                                                    break;
                                                }
                                            } // select search method
                                        }
                                        else {
                                            switch(search_method) {
                                                case REPAIR: { // repair + static ordering
                                                    vector<Vertex> ColOrd = P; // copy P
                                                    branch_repair(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max);
                                                    break;
                                                }
                                                case REPAIR_ONLY: {
                                                    branch_repair_only(vs,es, P, ind, C, C_max, colors, pruned, max);
                                                    break;
                                                }
                                                case STATIC_ORDER: { // static ordering only
                                                    vector<Vertex> ColOrd = P; // copy P
                                                    branch_static_order(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max);
                                                    break;
                                                }
                                                case DYNAMIC_KCORE: {
                                                    // kcore pruning + greedy coloring
                                                    break;
                                                }
                                                case DYNAMIC: {
                                                    // kcore pruning + greedy coloring
                                                    branch_dynamic(vs,es, P, ind, C, C_max, colors, pruned, max, S, levels, steps);
                                                    break;
                                                }
                                                case DYNAMIC_COLORING: {
                                                    // only greedy coloring
                                                    branch(vs,es,P, ind, C, C_max, colors, pruned, max);
                                                    break;
                                                }
                                                case BASIC: {
                                                    // no greedy coloring
                                                    branch_basic(P, ind, C, C_max, pruned, max);
                                                    break;
                                                }
                                                default: {
                                                    // no greedy coloring
                                                    branch_basic(P, ind, C, C_max, pruned, max);
                                                    break;
                                                }
                                            } // select search method
                                        }
                                    }
                                }

                            }
                            break;
                        } //end of CSC case
                    }
                    P = T;
                } // first bound
                pruned[u] = 1;
                (*degree)[u] = 0;
                (*bound)[u] = 0;

                if (G.is_dense_graph()) {
                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                        adj[u][es[j]] = false;
                        adj[es[j]][u] = false;
                        (*degree)[es[j]]--;
                    }
                }
                else {
                    for (long long j = vs[u]; j < vs[u + 1]; j++) { (*degree)[es[j]]--; }
                }

                // dynamically reduce graph in a thread-safe manner
                if ((get_time() - induce_time[omp_get_thread_num()]) > wait_time) {
                    G.reduce_graph( vs, es, pruned, G, i+lb_idx, max);
                    G.graph_stats(G, vs, es, max, i+lb_idx, sec);
                    induce_time[omp_get_thread_num()] = get_time();
                }
            }
        }
    }
    if (pruned) delete[] pruned;

    if (search_method == DYNAMIC_KCORE || search_method == DYNAMIC) {
        int total_steps = 0;
        for (int t = 0; t < num_threads; ++t) { total_steps += steps[t]; };
        printf("total steps = %d \n", total_steps);

        if (G.verbose) {
            for (int i = 0; i < max; ++i) {
                for (int t = 0; t < num_threads; ++t) {
                    cout << "tid = " << t+1 << ", depth = " << i << ", i1 = " << S[t][i].get_i1() << ", i2 = " << S[t][i].get_i2();
                    cout << ", level_ratio = " << (double)S[t][i].get_i2()/(double)S[t][i].get_i1() << ", ratio = " << (double)S[t][levels[t]].get_i1()/steps[t] <<endl;
                }
                printf("\n");
            }
        }
        printf("total steps = %d \n", total_steps);
    }
    print_line(80);
    if (p.k > 0 && max > p.k) max = p.k; // ensure clique of k returned
    printf("[mc] max = %d \n", max);

    // check solution from heuristic
    if (max == sol.size()) {
        cout << "solution from heu is optimal, heu = " << sol.size() << ", max = " << max <<endl;   G.print_break();
        return sol.size();
    }
    sol.swap(C_max);
    G.print_break();
    return sol.size();
}

//
// Neighborhood is Explicitly Induced.
// Hybrid using both CSC and Adjacency Matrix
//
void graphpack_mc::branch_neighborhood(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    branch_neighborhood(N, R, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) { // maximal clique, check if largest so far
                    // obtain lock
#pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_basic_neighborhood(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.size() > max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    branch_basic_neighborhood(N, R, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_neighborhood_static_order(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<Vertex> & ColOrd,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (N.adj[v][w])
                            if (!pruned[N.vertex_lookup[w]])
                                newColOrd.push_back(ColOrd[k]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_static(N.vertices, N.edges, R, newColOrd,
                            C, C_max, colors, max, N.adj);
                    branch_neighborhood_static_order(N, R, newColOrd, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_neighborhood_repair(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<Vertex> & ColOrd,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (N.adj[v][w])
                            if (!pruned[N.vertex_lookup[w]])
                                newColOrd.push_back(ColOrd[k]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair(N.vertices, N.edges, R, newColOrd,
                            C, C_max, colors, max, N.adj);
                    branch_neighborhood_repair(N, R, newColOrd, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
#pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_neighborhood_repair_only(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());


                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair_only(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    branch_neighborhood_repair_only(N, R, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
#pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

            }
            else return;
        }
    }
}

void graphpack_mc::branch_dynamic_neighborhood(
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
        vector<int> &steps) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {

                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {

                        switch (search_pruning) {
                            case KCORE: {
                                neigh_cores_bound(N.vertices,N.edges,R,ind,max);
                                break;
                            }
                            case TRIANGLES: {
                                // creates t array in N, does not sort!
                                N.compute_vertex_triangles(R,pruned,max,false);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }

                        if (search_ordering != COLORING && search_ordering != KCORE){ // && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,N,search_ordering,search_small_to_large,0,max);
                        }
                    }


                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    }
                    else {
                        greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    }

                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    branch_dynamic_neighborhood(N, R, C, C_max, ind, colors, pruned, max, S, levels, steps);
                    levels[tid]--;
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_neighborhood_cores(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<int>& C,
        vector<int>& C_max,
        vector<short>& ind,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    neigh_cores_bound(N.vertices,N.edges,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                    greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);
                    branch_neighborhood_cores(N, R, C, C_max, ind, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
#pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_dynamic_neighborhood_cores(
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
        vector<int> &steps) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {
                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {
                        neigh_cores_bound(N.vertices,N.edges,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                    }
                    greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj);

                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    branch_dynamic_neighborhood_cores(N, R, C, C_max, ind, colors, pruned, max, S, levels, steps);
                    levels[tid]--;
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}



//
// CSC ONLY -- Neighborhood is not explicitly induced.
//
void graphpack_mc::branch(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());
                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1) bound check
                    neigh_coloring_bound(vs, es, R, ind, C, C_max, colors, pruned, max);
                    // search reordered R
                    branch(vs, es, R, ind, C, C_max, colors, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }

                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}


void graphpack_mc::branch_basic(
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size()) {
            // terminating condition
            if (C.size() + P.size() > max) {
                int v = P.back().get_id();  // remove v from P
                C.push_back(v);             // add v to C

                vector<Vertex> R;
                R.reserve(P.size());

                // mark neighbors of v
                for (long long j = (*vertices)[v]; j < (*vertices)[v + 1]; j++)   ind[(*edges)[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);

                for (long long j = (*vertices)[v]; j < (*vertices)[v + 1]; j++)  ind[(*edges)[j]] = 0;

                if (R.size()) {
                    branch_basic(R, ind, C, C_max, pruned, max);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }

                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_static_order(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<Vertex> & ColOrd,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);


                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (ind[w])
                            if (!pruned[w])
                                if ((*bound)[w] > max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1) bound check
                    greedy_neighborhood_coloring_static(vs,es,R,newColOrd,C,C_max,colors,max, ind); // csc only -- sparse
//                    neigh_coloring_bound(vs, es, R, ind, C, C_max, colors, pruned, max);
                    // search reordered R
                    branch_static_order(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max);
                }


                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_repair(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<Vertex> & ColOrd,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);


                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (ind[w])
                            if (!pruned[w])
                                if ((*bound)[w] > max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1) bound check
                    greedy_neighborhood_coloring_repair(vs,es,R,newColOrd,C,C_max,colors,max, ind); // csc only -- sparse
//                    neigh_coloring_bound(vs, es, R, ind, C, C_max, colors, pruned, max);
                    // search reordered R
                    branch_repair(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max);
                }


                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_repair_only(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);


                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1) bound check
                    greedy_neighborhood_coloring_repair_only(vs,es,R,C,C_max,colors,max, ind); // csc only -- sparse
                    branch_repair_only(vs, es, R, ind, C, C_max, colors, pruned, max);
                }


                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_dynamic(
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
        vector<int> &steps) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (ind[P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 0;

                if (R.size() > 0) {

                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {

                        switch (search_pruning) {
                            case KCORE: {
                                // prunes verts from R and sorts automatically
                                neigh_cores_bound(vs,es,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                                break;
                            }
                            case TRIANGLES: {
                                // saves triangles in R, does not sort!
                                compute_vertex_triangles_dense(R,pruned,max,false, vs, es);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }
                        if (search_ordering != COLORING && search_ordering != KCORE && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,vs,es,search_ordering,search_small_to_large,0,max);
                        }
                    }

                    // color graph induced by R and sort for O(1) bound check
                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, ind);
                    }
                    else {
                        neigh_coloring_bound(vs, es, R, ind, C, C_max, colors, pruned, max);
                    }



                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    branch_dynamic(vs, es, R, ind, C, C_max, colors, pruned, max, S, levels, steps);
                    levels[tid]--;
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}





//
// CSC/ADJ ONLY -- Neighborhood is not explicitly induced.
//
void graphpack_mc::branch_dense_static_order(
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
        bool** &adj) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);
                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned, since threads dynamically update it
                    if (adj[v][P[k].get_id()])
                        if ((*bound)[P[k].get_id()] > max)
                            R.push_back(P[k]);

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (adj[v][w])
                            if (!pruned[w])
                                newColOrd.push_back(ColOrd[k]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    neigh_coloring_dense(vs, es, R, ind, C, C_max, colors, max, adj);
                    greedy_neighborhood_coloring_static(vs, es, R, newColOrd, C, C_max, colors, max, adj);

                    branch_dense_static_order(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, adj);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_basic_dense(
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        int* &pruned,
        int& max,
        bool** &adj) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.size() > max) {
                int v = P.back().get_id();   C.push_back(v);
                vector<Vertex> R;    R.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned
                    if (adj[v][P[k].get_id()])
                        if ((*bound)[P[k].get_id()] > max)
                            R.push_back(P[k]);

                if (R.size() > 0) {
                    branch_basic_dense(R, ind, C, C_max, pruned, max, adj);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }

                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_dense(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        bool** &adj) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);
                vector<Vertex> R;
                R.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned, since threads dynamically update it
                    if (adj[v][P[k].get_id()])
                        if ((*bound)[P[k].get_id()] > max)
                            R.push_back(P[k]);

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    neigh_coloring_dense(vs, es, R, ind, C, C_max, colors, max, adj);
                    branch_dense(vs, es, R, ind, C, C_max, colors, pruned, max, adj);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }

                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}

void graphpack_mc::branch_dense_repair(
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
        bool** &adj) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();

                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (adj[v][u] && !pruned[u])
                        if ((*bound)[u] > max)
                            R.push_back(P[k]);
                }

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (adj[v][w])
                            if (!pruned[w])
                                if ((*bound)[w] > max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair(vs, es, R, newColOrd, C, C_max, colors, max, adj);
                    branch_dense_repair(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, adj);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_dense_repair_only(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        bool** &adj) {

    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {

                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (adj[v][u] && !pruned[u])
                        if ((*bound)[u] > max)
                            R.push_back(P[k]);
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, adj);
                    branch_dense_repair_only(vs, es, R, ind, C, C_max, colors, pruned, max, adj);
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();
            }
            else return;
        }
    }
}

void graphpack_mc::branch_dense_dynamic(
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
        vector<int> &steps) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


    // stop early if ub is reached
    if (not_reached_ub) {
        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() > max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] > max)
                                R.push_back(P[k]);
                }

                if (R.size() > 0) {

                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {

                        switch (search_pruning) {
                            case KCORE: {
                                neigh_cores_bound(vs,es,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                                break;
                            }
                            case TRIANGLES: {
                                // creates t array in N, does not sort!
                                compute_vertex_triangles_dense(R,pruned,max,false, vs, es, adj);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }

                        if (search_ordering != COLORING && search_ordering != KCORE && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,vs,es,search_ordering,search_small_to_large,0,max);
                        }
                    }

                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, adj);
                    }
                    else {
                        greedy_neighborhood_coloring(vs, es, R, C, C_max, colors, max, adj);
                    }

                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    branch_dense_dynamic(vs, es, R, C, C_max, ind, colors, pruned, max, adj, S, levels, steps);
                    levels[tid]--;
                }
                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
    }
}




/**
 * MAXIMUM CLIQUE ENUMERATION and K-CLIQUE ENUMERATION
 */
int graphpack_mc::search_enum(graphpack_graph& G, set< vector<int> >& sol, params &p) {

    vertices = G.get_vertices();
    edges = G.get_edges();
    degree = G.get_degree();

    lb = lb - 1;

    int* pruned = new int[G.num_vertices()];
    memset(pruned, 0, G.num_vertices() * sizeof(int));
    int max = lb;
    int u = 0;

    printf("enumerating cliques k = %d, lb = %d, ub = %d \n", k_clique_size, lb, ub);

    // initial pruning
    int lb_idx = 0;
    switch (global_pruning) {
        case KCORE: {
            printf("pruning method: k-core \n");
            lb_idx = pruning(G, max, p.global_pruning, pruned, false, is_enum);
            break;
        }
        case TRIANGLES: {
            // selects between parallel and sequential or csc only vs. hybrid of csc/adj
            printf("pruning method: k-core and triangles \n");

            // first prune using k-cores
            lb_idx = pruning(G, max, "kcore", pruned, false, is_enum);

            // then prune via triangles
            G.compute_vertex_triangles(pruned, max, true);
            break;
        }
        case TRIANGLES_ONLY: {
            printf("pruning method: triangles only \n");
            // selects between parallel and sequential
            // or csc only vs. hybrid of csc/adj
            G.compute_vertex_triangles(pruned, max, true, is_enum);
            break;
        }
        case TRIANGLE_CORES: {
            vector<int> W;
            for (int i = G.num_vertices()-1; i >= 0; i--) {
                if (G.kcore[G.kcore_order[i]] >= max)
                    W.push_back(G.kcore_order[i]);
            }
            graphpack_graph H(W,G.vertices,G.edges);
            H.basic_stats();

            H.map_edgeid_to_vertex_pair(true);

            vector<int> pruned_edges(G.num_edges()+1,0);
            H.k_triangle_counting(max,pruned_edges, true);
            printf("[triangle cores]  argmax tr(u) = %lld,  sum_u tr(u) = %llu \n", H.max_t_edge, H.total_t);

            // compute the triangle cores
            H.triangle_core_numbers_parallel();

            /**
             * this function explicitly removes edges,
             * but only marks vertices as pruned
             */
            triangle_core_remove_edges(H,pruned,pruned_edges,max);
            H.compute_cores();
            printf("updated max core = %d \n", H.max_core);
            G = H;
            break;
        }
    } // end of global_pruning switch


    cout << "[maxpack: initial k-core pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
    G.reduce_graph(pruned);
    G.update_degrees();
    cout << "[graphpack: initial k-core pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;



    if (p.global_edge_ordering == "kcore") {
        printf("ordering edges by k-cores \n");
        edge_bucket_sort_parallel(G, G.kcore, p);
    }
    else {
        printf("ordering edges by degree (default) \n");
        edge_bucket_sort_parallel(G, G.degree, p);
    }

    // copy G.adj to local adj if possible
    bool** adj;
    if (G.is_dense_graph()) {
        adj = G.adj;
        cout << "[mcpack] the graph G is dense! " <<endl;
    }
    else  cout << "[mcpack] the graph G is sparse! " <<endl;

    // set to worst case bound of cores/coloring
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);

    vector<int> C, C_max;
    C.reserve(G.get_max_degree()+1);
    C_max.reserve(G.get_max_degree()+1);

    // init the neigh coloring array
    vector< vector<int> > colors(G.get_max_degree()+1);
    for (int i = 0; i < G.get_max_degree()+1; i++)  colors[i].reserve(G.get_max_degree()+1);


    // order verts for our search routine
    vector<Vertex> V;
    V.reserve(G.num_vertices());
    compute_ordering(V,G,global_ordering,p.global_small_to_large,lb_idx,lb);
    cout << "|V| = " << V.size() << ", f(0) = " << V[0].get_bound() << ", f(|V|) = " << V[V.size()-1].get_bound() <<endl;

    vector<short> ind(G.num_vertices(),0);
    vector<int> es = G.get_edges_array();
    vector<long long> vs = G.get_vertices_array();

    vector<double> induce_time(num_threads,get_time());
    for (int t = 0; t < num_threads; ++t)  induce_time[t] = induce_time[t] + t; // t/4;


    // setup dynamic search data structures
    vector< vector<OpCount> > S;
    S.resize(num_threads+1);
    if (search_method == DYNAMIC || search_method == DYNAMIC_KCORE) {
        printf("setting up dynamic search data structures \n");
        for (int t = 0; t < num_threads; ++t) S[t].resize(G.get_max_core() + 1);
        cout << "setting up dynamic step counting structs" <<endl;
        for (int t = 0; t < num_threads; ++t) {
            for (int i=0; i < G.get_max_core() + 1; ++i) {
                S[t][i].set_i1(0);
                S[t][i].set_i2(0);
            }
        }
    }
    cout << "num threads = " << num_threads <<endl;

    vector<int> levels(num_threads+1,0);
    vector<int> steps(num_threads+1,0);


    vector< set< vector<int> > > cliques_local(num_threads);

    #pragma omp parallel for schedule(dynamic,block_size) shared(pruned, G, T, V, max, C_max, induce_time) \
        firstprivate(colors,ind,vs,es) private(u, P, C)
    for (int i = 0; i < (V.size()); ++i) {
        if (G.time_left(C_max,sec,time_limit,time_expired_msg)) {

            u = V[i].get_id();
            if ((*bound)[u] >= max) {
                P.push_back(V[i]);
                for (long long j = vs[u]; j < vs[u + 1]; ++j)
                    if (!pruned[es[j]])
                        if ((*bound)[es[j]] >= max)
                            P.push_back(Vertex(es[j], (*degree)[es[j]]));


                switch(local_representation) {
                    case HYBRID: {
                        if (P.size() >= max) { // degree bound

                            // neighborhood core ordering and pruning
                            if (use_neigh_cores) {
                                neigh_cores_enum(vs,es,P,ind,max);
                                if (!(P.size() >= max && P[0].get_bound()+1 >= max)) {
                                    P = T;
                                    break;
                                }
                            }

                            if (P.size() >= max) {
                                graphpack_graph N(P,vs,es,max,p.local_pruning == "tcores",is_enum); // true if tcores, otherwise false

                                if (P.size() >= max) { // if (P.size() <= max) continue;

                                    if (N.density() > p.density_cutoff) { // for flickr and other dense graphs
                                        if (G.verbose) printf("neighborhood density = %lg \n", N.density());

                                        switch(local_pruning) {
                                            case TRIANGLES: {
                                                N.compute_vertex_triangles(P,pruned,max,false);
                                                break;
                                            }
                                            case TRIANGLE_CORES: {
                                                N.triangle_cores_serial_adj(max);

                                                if (N.get_triangle_core_bound() >= max) {

                                                    // only prune if neighborhood triangle-core bound passes
                                                    N.triangle_core_pruning(P,max);
                                                    if (G.verbose) printf("p = %lg, T_v = %lld, |P| = %lu \n", N.density(), N.get_triangle_core_bound(), P.size());
                                                    // only sort if P has enough vertices after pruning via triangle-cores
                                                    if (P.size() <= max) { P = T; }
                                                } // bound
                                                else { // set P to be empty, so the check below fails
                                                    P = T;
                                                }
                                                break;
                                            }
                                        }
                                    }




                                    if (P.size() >= max) {
                                        if (G.is_dense_graph())
                                            neigh_coloring_dense_enum(N.vertices,N.edges,P,ind,C,C_max,colors,max, N.adj);
                                        else
                                            neigh_coloring_enum(N.vertices,N.edges,P,ind,C,C_max,colors,pruned,max);


                                        if (P.back().get_bound() >= max) {

                                            if (local_ordering != COLORING) {
                                                compute_local_ordering(P,N,local_ordering,p.local_small_to_large,0,max);
                                            }

                                            switch(search_method) {
                                                case REPAIR: { // repair + static ordering
                                                    //                                                    branch_neighborhood_enum(N,P,ind,C,C_max,colors,pruned,max,cliques_local);
                                                    vector<Vertex> ColOrd = P; // copy P
                                                    clique_enum_neighborhood_repair(N, P, ColOrd, C, C_max, colors, pruned, max, cliques_local);
                                                    break;
                                                }
                                                case REPAIR_ONLY: {
                                                    clique_enum_neighborhood_repair_only(N, P, C, C_max, colors, pruned, max, cliques_local);
                                                    break;
                                                }
                                                case STATIC_ORDER: { // static ordering only
                                                    vector<Vertex> ColOrd = P;
                                                    clique_enum_neighborhood_static_order(N, P, ColOrd, C, C_max, colors, pruned, max, cliques_local);
                                                    break;
                                                }
                                                case DYNAMIC: { // static ordering only
                                                    clique_enum_dynamic_neighborhood(N, P, C, C_max, ind, colors, pruned, max, S, levels, steps, cliques_local);
                                                    break;
                                                }
                                                case DYNAMIC_COLORING: {
                                                    // only greedy coloring
                                                    clique_enum_neighborhood(N, P, C, C_max, colors, pruned, max, cliques_local);
                                                    break;
                                                }
                                                case BASIC: {
                                                    clique_enum_basic_neighborhood(N,P,C,C_max,colors,pruned,max,cliques_local);
                                                    break;
                                                }
                                                default: {
                                                    clique_enum_basic_neighborhood(N,P,C,C_max,colors,pruned,max,cliques_local);
                                                    break;
                                                }
                                            } // select search method
                                        }
                                    }
                                }
                            }
                        }
                        break;
                    } //end of HYBRID case
                    case CSC: {
                        if (P.size() >= max) { // degree bound
                            // neighborhood core ordering and pruning
                            if (use_neigh_cores) {
                                neigh_cores_enum(vs,es,P,ind,max);
                                if (! (P.size() >= max && P[0].get_bound()+1 >= max) ) {
                                    P = T;
                                    break;
                                }
                            }

                            if (P.size() >= max) {

                                if (G.is_dense_graph()) neigh_coloring_dense_enum(vs,es,P,ind,C,C_max,colors,max, adj);
                                else                    neigh_coloring_enum(vs,es,P,ind,C,C_max,colors,pruned,max);


                                if (P.back().get_bound() >= max) {

                                    if (p.local_ordering != "coloring")
                                        compute_local_ordering(P,G,vs,es,local_ordering,p.local_small_to_large,0,max);

                                    if (G.is_dense_graph()) {
                                        switch(search_method) {
                                            case REPAIR: { // repair + static ordering
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_dense_repair(vs,es, P, ColOrd,ind, C, C_max, colors, pruned, max, adj, cliques_local);
                                                break;
                                            }
                                            case REPAIR_ONLY: {
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_dense_repair_only(vs,es, P,ind, C, C_max, colors, pruned, max, adj, cliques_local);
                                                break;
                                            }
                                            case STATIC_ORDER: { // static ordering only
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_dense_static_order(vs,es, P, ColOrd,ind, C, C_max, colors, pruned, max, adj, cliques_local);
                                                break;
                                            }
                                            case DYNAMIC: {
                                                clique_enum_dense_dynamic(vs,es,P,C,C_max,ind,colors,pruned,max,adj,S,levels,steps,cliques_local);
                                                break;
                                            }
                                            case DYNAMIC_COLORING: {
                                                clique_enum_dense(vs,es,P,ind,C,C_max,colors,pruned,max,adj,cliques_local);
                                                break;
                                            }
                                            case BASIC: {
                                                // only greedy coloring
                                                clique_enum_basic_dense(P,ind,C,C_max,pruned,max,adj,cliques_local);
                                                break;
                                            }
                                            default: {
                                                // only greedy coloring
                                                clique_enum_basic_dense(P,ind,C,C_max,pruned,max,adj,cliques_local);
                                                break;
                                            }
                                        } // switch
                                    } // if dense
                                    else {
                                        switch(search_method) {
                                            case REPAIR: { // repair + static ordering
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_repair(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max, cliques_local);
                                                break;
                                            }
                                            case REPAIR_ONLY: {
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_repair_only(vs,es, P, ind, C, C_max, colors, pruned, max, cliques_local);
                                                break;
                                            }
                                            case STATIC_ORDER: { // static ordering only
                                                vector<Vertex> ColOrd = P;
                                                clique_enum_static_order(vs,es, P, ColOrd, ind, C, C_max, colors, pruned, max, cliques_local);
                                                break;
                                            }
                                            case DYNAMIC: {
                                                clique_enum_dynamic(vs,es,P,ind,C,C_max,colors,pruned,max,S,levels,steps,cliques_local);
                                                break;
                                            }
                                            case DYNAMIC_COLORING: {
                                                // only greedy coloring
                                                clique_enum(vs,es,P, ind, C, C_max, colors, pruned, max, cliques_local);
                                                break;
                                            }
                                            case BASIC: {
                                                // no greedy coloring
                                                clique_enum_basic(P, ind, C, C_max, pruned, max, cliques_local);
                                                break;
                                            }
                                            default: {
                                                clique_enum_basic(P, ind, C, C_max, pruned, max, cliques_local);
                                                break;
                                            }
                                        } // switch
                                    } // sparse
                                }
                            }
                        }
                        break;
                    } //end of CSC case
                }
                P = T;
            }
//            pruned[u] = 1;
//            if (G.is_dense_graph()) {
//                for (long long j = vs[u]; j < vs[u + 1]; j++) {
//                    adj[u][es[j]] = false;
//                    adj[es[j]][u] = false;
//                }
//            }
            // dynamically reduce graph in a thread-safe manner
//            if ((get_time() - induce_time[omp_get_thread_num()]) > wait_time) {
//                G.reduce_graph( vs, es, pruned, G, i+lb_idx, max);
//                G.graph_stats(G, vs, es, max, i+lb_idx, sec);
//                induce_time[omp_get_thread_num()] = get_time();
//            }
        }
    }

    if (pruned) delete[] pruned;

    if (search_method == DYNAMIC_KCORE || search_method == DYNAMIC) {
        int total_steps = 0;
        for (int t = 0; t < num_threads; ++t) { total_steps += steps[t]; };
        printf("total steps = %d \n", total_steps);

        if (G.verbose) {
            for (int i = 0; i < 10; ++i) {
                for (int t = 0; t < num_threads; ++t)
                    cout << "tid = " << t+1 << ", depth = " << i << ", i1 = " << S[t][i].get_i1() << ", i2 = " << S[t][i].get_i2() <<endl;
                printf("\n");
            }
        }
    }
    print_line(80);

    double s = tic();
    sol = cliques_local[0];
    for (int t = 1; t < num_threads; ++t) {
        add_all_cliques(cliques_local[t], sol);
        printf("[ pmc worker %d ]  found %lu cliques of size %d, total unique cliques = %lu \n",
                t, cliques_local[t].size(), max, sol.size());
    }
    toc(s);
    printf("took %lg seconds to find unique cliques \n", s);

    G.print_break();
    return max;
}

void graphpack_mc::branch_neighborhood_enum(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<short> &ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

    while (P.size() > 0) {
        // terminating condition
        if (C.size() + P.back().get_bound() >= max) {

            int v = P.back().get_id();
            C.push_back(N.vertex_lookup[v]);

            vector<Vertex> R;
            R.reserve(P.size());

            // intersection of N(v) and P - {v}
            for (int k = 0; k < P.size() - 1; k++) {
                if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                    if (!pruned[N.vertex_lookup[P[k].get_id()]])
                        if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                            R.push_back(P[k]);
            }


            // update only if ub has not been reached yet (for k-clique enumeration)
            if (C.size() > max && not_reached_ub) {
                // obtain lock
                #pragma omp critical (update_mc)
                if (C.size() > max && not_reached_ub) {

                    // reset all worker sets
                    for (int t = 0; t < num_threads; ++t)   cliques_local[t].clear();

                    // share global bound
                    max = C.size();
                    C_max = C;

                    print_mc_info(C_max,sec);
                    if (max >= k_clique_size) {
                        not_reached_ub = false;
                        cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                    }

                    // add clique to local worker set
                    add_clique(C,cliques_local[omp_get_thread_num()]);
                }
            }
            if (C.size() == max) {
                //                    #pragma omp critical (update_mc)
                if (C.size() == max) {
                    add_clique(C,cliques_local[omp_get_thread_num()]);
                }
            }


            if (R.size() > 0) {
                // color graph induced by R and sort for O(1)
                greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj, true);
                branch_neighborhood_enum(N, R, ind, C, C_max, colors, pruned, max, cliques_local);
            }
            // backtrack and search another branch
            R.clear();
            C.pop_back();
        }
        else return;
        P.pop_back();
    }
}

void graphpack_mc::branch_enum(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

    while (P.size() > 0) {
        // terminating condition
        if (C.size() + P.back().get_bound() >= max) {
            int v = P.back().get_id();   C.push_back(v);

            vector<Vertex> R;   R.reserve(P.size());
            for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

            // intersection of N(v) and P - {v}
            for (int k = 0; k < P.size() - 1; k++)
                if (ind[P[k].get_id()])
                    if (!pruned[P[k].get_id()])
                        if ((*bound)[P[k].get_id()] >= max)
                            R.push_back(P[k]);

            for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;


            // update only if ub has not been reached yet (for k-clique enumeration)
            if (C.size() > max && not_reached_ub) {
                // obtain lock
                #pragma omp critical (update_mc)
                if (C.size() > max && not_reached_ub) {

                    // reset all worker sets
                    for (int t = 0; t < num_threads; ++t)
                        cliques_local[t].clear();

                    // share global bound
                    max = C.size();
                    C_max = C;

                    // add clique to local worker set
                    add_clique(C,cliques_local[omp_get_thread_num()]);


                    print_mc_info(C_max,sec);
                    if (max >= k_clique_size) {
                        not_reached_ub = false;
                        cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                    }
                }
            }

            if (C.size() == max) {
//                #pragma omp critical (update_mc)
                if (C.size() == max) {
                    add_clique(C,cliques_local[omp_get_thread_num()]);
                }
            }


            if (R.size() > 0) {
                // color graph induced by R and sort for O(1)
                neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                branch_enum(vs, es, R, ind, C, C_max, colors, pruned, max, cliques_local);
            }

            // backtrack and search another branch
            R.clear();
            C.pop_back();
        }
        else return;
        P.pop_back();
    }

}

void graphpack_mc::branch_dense_enum(
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
        vector< set< vector<int> > > &cliques_local) {

    while (P.size() > 0) {
        // terminating condition
        if (C.size() + P.back().get_bound() >= max) {
            int v = P.back().get_id();   C.push_back(v);
            vector<Vertex> R;    R.reserve(P.size());

            for (int k = 0; k < P.size() - 1; k++)
                // indicates neighbor AND pruned, since threads dynamically update it
                if (adj[v][P[k].get_id()])
                    if ((*bound)[P[k].get_id()] >= max)
                        R.push_back(P[k]);


            // update only if ub has not been reached yet (for k-clique enumeration)
            if (C.size() > max && not_reached_ub) {
                // obtain lock
                #pragma omp critical (update_mc)
                if (C.size() > max && not_reached_ub) {

                    for (int t = 0; t < num_threads; ++t)
                        cliques_local[t].clear();

                    max = C.size();
                    C_max = C;

                    add_clique(C,cliques_local[omp_get_thread_num()]);


                    print_mc_info(C_max,sec);
                    if (max >= k_clique_size) {
                        not_reached_ub = false;
                        cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                    }
                }
            }

            if (C.size() == max) {
//                #pragma omp critical (update_mc)
                if (C.size() == max) {
                    add_clique(C,cliques_local[omp_get_thread_num()]);
                }
            }

            if (R.size() > 0) {
                // color graph induced by R and sort for O(1)
                neigh_coloring_dense_enum(vs, es, R, ind, C, C_max, colors, max, adj);
                branch_dense_enum(vs, es, R, ind, C, C_max, colors, pruned, max, adj, cliques_local);
            }

            // backtrack and search another branch
            R.clear();
            C.pop_back();
        }
        else return;
        P.pop_back();
    }
}



//
// Neighborhood is Explicitly Induced.
// Hybrid using both CSC and Adjacency Matrix
//
void graphpack_mc::clique_enum_neighborhood(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj, true);
                    clique_enum_neighborhood(N, R,  C, C_max, colors, pruned, max, cliques_local);
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_basic_neighborhood(
        graphpack_graph & N,
        vector<Vertex> &P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }
//                cout << "here" <<endl;

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (R.size() > 0) {
                    clique_enum_neighborhood(N, R,  C, C_max, colors, pruned, max, cliques_local);
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_neighborhood_static_order(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<Vertex> & ColOrd,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
                int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (N.adj[v][w])
                            if (!pruned[N.vertex_lookup[w]])
                                newColOrd.push_back(ColOrd[k]);
                    }
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {

                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_static(N.vertices, N.edges, R, newColOrd, C, C_max, colors, max, N.adj, is_enum);
                    clique_enum_neighborhood_static_order(N, R, newColOrd, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_neighborhood_repair(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<Vertex> & ColOrd,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }
//                cout << "here" << endl;

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (N.adj[v][w])
                            if (!pruned[N.vertex_lookup[w]])
                                if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {

                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair(N.vertices, N.edges, R, newColOrd, C, C_max, colors, max, N.adj, is_enum);
                    clique_enum_neighborhood_repair(N, R, newColOrd, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_neighborhood_repair_only(
        graphpack_graph & N,
        vector<Vertex> & P,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();

                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());


                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (N.adj[v][u])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[u]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair_only(N.vertices, N.edges, R, C, C_max, colors, max, N.adj, is_enum);
                    clique_enum_neighborhood_repair_only(N, R, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_dynamic_neighborhood(
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
        vector< set< vector<int> > > &cliques_local) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(N.vertex_lookup[v]);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (N.adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[N.vertex_lookup[P[k].get_id()]])
                            if ((*bound)[N.vertex_lookup[P[k].get_id()]] >= max)
                                R.push_back(P[k]);
                }


                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (C.size() == max) {
//                    #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {

                        switch (search_pruning) {
                            case KCORE: {
                                neigh_cores_enum(N.vertices,N.edges,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                                break;
                            }
                            case TRIANGLES: {
                                N.compute_vertex_triangles(R,pruned,max,false);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }

                        if (search_ordering != COLORING && search_ordering != KCORE && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,N,search_ordering,search_small_to_large,0,max);
                        }
                    }


                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(N.vertices, N.edges, R, C, C_max, colors, max, N.adj, is_enum);
                    }
                    else {
                        greedy_neighborhood_coloring(N.vertices, N.edges, R, C, C_max, colors, max, N.adj, is_enum);
                    }

                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    clique_enum_dynamic_neighborhood(N, R, C, C_max, ind, colors, pruned, max, S, levels, steps, cliques_local);
                    levels[tid]--;
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}



//
// CSC ONLY -- Neighborhood is not explicitly induced.
//
void graphpack_mc::clique_enum(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());
                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;


                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                    clique_enum(vs, es, R, ind, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_basic(
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();   C.push_back(v);

                vector<Vertex> R;   R.reserve(P.size());
                for (long long j = (*vertices)[v]; j < (*vertices)[v + 1]; j++)   ind[(*edges)[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);

                for (long long j = (*vertices)[v]; j < (*vertices)[v + 1]; j++)  ind[(*edges)[j]] = 0;

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                    clique_enum_basic(R, ind, C, C_max, pruned, max, cliques_local);
//                    clique_enum(vs, es, R, ind, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_static_order(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);


                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (ind[w])
                            if (!pruned[w])
                                if ((*bound)[w] > max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                    greedy_neighborhood_coloring_static(vs,es,R,newColOrd,C,C_max,colors,max, ind, true); // csc only -- sparse
                    clique_enum_static_order(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_repair(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);


                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (ind[w])
                            if (!pruned[w])
                                if ((*bound)[w] > max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;


                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                    greedy_neighborhood_coloring_repair(vs,es,R,newColOrd,C,C_max,colors,max, ind, true); // csc only -- sparse
                    clique_enum_repair(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_repair_only(
        vector<long long>& vs,
        vector<int>& es,
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        vector< vector<int> >& colors,
        int* &pruned,
        int& max,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++)
                    if (ind[P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);


                for (long long j = vs[v]; j < vs[v + 1]; j++)  ind[es[j]] = 0;

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1) bound check
                    greedy_neighborhood_coloring_repair_only(vs,es,R,C,C_max,colors,max, ind, is_enum); // csc only -- sparse
                    clique_enum_repair_only(vs, es, R, ind, C, C_max, colors, pruned, max, cliques_local);
                }


                else if (C.size() > max) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max) {
                        // ensure updated max is flushed
                        max = C.size();
                        C_max = C;
                        print_mc_info(C,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_dynamic(
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
        vector< set< vector<int> > > &cliques_local) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 1;

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (ind[P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);
                }

                for (long long j = vs[v]; j < vs[v + 1]; j++)   ind[es[j]] = 0;



                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {

                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {

                        switch (search_pruning) {
                            case KCORE: {
//                              // prunes verts from R and sorts automatically
                                neigh_cores_enum(vs,es,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                                break;
                            }
                            case TRIANGLES: {
                                // saves triangles in R, does not sort!
                                compute_vertex_triangles_dense(R,pruned,max,false, vs, es);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }

                        if (search_ordering != COLORING && search_ordering != KCORE && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,vs,es,search_ordering,search_small_to_large,0,max);
                        }
                    }

                    // color graph induced by R and sort for O(1) bound check
                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, ind, is_enum);
                    }
                    else {
                        neigh_coloring_enum(vs, es, R, ind, C, C_max, colors, pruned, max);
                    }



                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    clique_enum_dynamic(vs, es, R, ind, C, C_max, colors, pruned, max, S, levels, steps, cliques_local);
                    levels[tid]--;
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}


//
// CSC/ADJ ONLY -- Neighborhood is not explicitly induced.
//
void graphpack_mc::clique_enum_dense_static_order(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {

            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);
                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned, since threads dynamically update it
                    if (adj[v][P[k].get_id()])
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (adj[v][w])
                            if (!pruned[w])
                                if ((*bound)[P[k].get_id()] >= max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);

                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_static(vs, es, R, newColOrd, C, C_max, colors, max, adj, is_enum);
                    clique_enum_dense_static_order(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, adj, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_basic_dense(
        vector<Vertex> &P,
        vector<short>& ind,
        vector<int>& C,
        vector<int>& C_max,
        int* &pruned,
        int& max,
        bool** &adj,
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();   C.push_back(v);
                vector<Vertex> R;    R.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned
                    if (adj[v][P[k].get_id()])
                        if ((*bound)[P[k].get_id()] >= max)
                            R.push_back(P[k]);

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    clique_enum_basic_dense(R, ind, C, C_max, pruned, max, adj, cliques_local);
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_dense(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);
                vector<Vertex> R;
                R.reserve(P.size());

                for (int k = 0; k < P.size() - 1; k++)
                    // indicates neighbor AND pruned, since threads dynamically update it
                    if (adj[v][P[k].get_id()])
                        if ((*bound)[P[k].get_id()] >= max)
                            R.push_back(P[k]);

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }


                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
//                    neigh_coloring_dense(vs, es, R, ind, C, C_max, colors, max, adj, is_enum);
                    greedy_neighborhood_coloring(vs, es, R, C, C_max, colors, max, adj, is_enum);
                    clique_enum_dense(vs, es, R, ind, C, C_max, colors, pruned, max, adj, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

void graphpack_mc::clique_enum_dense_repair(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();

                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                vector<Vertex> newColOrd;
                newColOrd.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (adj[v][u] && !pruned[u])
                        if ((*bound)[u] >= max)
                            R.push_back(P[k]);
                }

                int v_pos = 0;
                for (int k = 0; k < ColOrd.size(); k++) {

                    int w = ColOrd[k].get_id();
                    if (w == v) { // mark it
                        v_pos = k;
                    }
                    else {
                        if (adj[v][w])
                            if (!pruned[w])
                                if ((*bound)[w] >= max)
                                    newColOrd.push_back(ColOrd[k]);
                    }
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair(vs, es, R, newColOrd, C, C_max, colors, max, adj, is_enum);
                    clique_enum_dense_repair(vs, es, R, newColOrd, ind, C, C_max, colors, pruned, max, adj, cliques_local);
                }

                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();

                // remove vertex from P from ColOrd
                bool found = false;
                for (int k = v_pos+1; k < ColOrd.size(); k++) {
                    ColOrd[k-1].set_id(ColOrd[k].get_id());
                }
                ColOrd.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_dense_repair_only(
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
        vector< set< vector<int> > > &cliques_local) {

        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {

                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    int u = P[k].get_id();
                    if (adj[v][u] && !pruned[u])
                        if ((*bound)[u] >= max)
                            R.push_back(P[k]);
                }


                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (R.size() > 0) {
                    // color graph induced by R and sort for O(1)
                    greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, adj, is_enum);
                    clique_enum_dense_repair_only(vs, es, R, ind, C, C_max, colors, pruned, max, adj, cliques_local);
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
                P.pop_back();
            }
            else return;
        }
}

void graphpack_mc::clique_enum_dense_dynamic(
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
        vector< set< vector<int> > > &cliques_local) {

    int tid = omp_get_thread_num();
    S[tid][levels[tid]].set_i1(S[tid][levels[tid]].get_i1() +
            S[tid][levels[tid] - 1].get_i1() - S[tid][levels[tid]].get_i2());
    S[tid][levels[tid]].set_i2(S[tid][levels[tid] - 1].get_i1());


        while (P.size() > 0) {
            // terminating condition
            if (C.size() + P.back().get_bound() >= max) {
                int v = P.back().get_id();
                C.push_back(v);

                vector<Vertex> R;
                R.reserve(P.size());

                // intersection of N(v) and P - {v}
                for (int k = 0; k < P.size() - 1; k++) {
                    if (adj[v][P[k].get_id()])   // check if neighbor of v is neighbor of w in P
                        if (!pruned[P[k].get_id()])
                            if ((*bound)[P[k].get_id()] >= max)
                                R.push_back(P[k]);
                }

                // update only if ub has not been reached yet (for k-clique enumeration)
                if (C.size() > max && not_reached_ub) {
                    // obtain lock
                    #pragma omp critical (update_mc)
                    if (C.size() > max && not_reached_ub) {

                        // reset all worker sets
                        for (int t = 0; t < num_threads; ++t)
                            cliques_local[t].clear();

                        // share global bound
                        max = C.size();
                        C_max = C;

                        // add clique to local worker set
                        add_clique(C,cliques_local[omp_get_thread_num()]);


                        print_mc_info(C_max,sec);
                        if (max >= k_clique_size) {
                            not_reached_ub = false;
                            cout << "[graphpack: upper bound reached]  omega = " << max <<endl;
                        }
                    }
                }

                if (C.size() == max) {
    //                #pragma omp critical (update_mc)
                    if (C.size() == max) {
                        add_clique(C,cliques_local[omp_get_thread_num()]);
                    }
                }

                if (R.size() > 0) {

                    if ((double)S[tid][levels[tid]].get_i1()/++steps[tid] < threshold) {
                        // color graph induced by R and sort for O(1)

                        switch (search_pruning) {
                            case KCORE: {
//                                printf("dynamic kcores \n");
                                neigh_cores_enum(vs,es,R,ind,max); // todo: use bound, C.size() + P.size() etc..
                                break;
                            }
                            case TRIANGLES: {
                                // creates t array in N, does not sort!
                                compute_vertex_triangles_dense(R,pruned,max,false, vs, es, adj);
                                if (search_ordering == TRIANGLES) {
                                    if (search_small_to_large)  std::sort(R.begin(), R.end(), incr_bound); // smallest to largest
                                    else                        std::sort(R.begin(), R.end(), decr_bound); // largest to smallest
                                }
                                break;
                            }
                            case TRIANGLE_CORES: {
                                break;
                            }
                        }

                        if (search_ordering != COLORING && search_ordering != KCORE && search_ordering != TRIANGLES) {
                            // degree, triangles, ...
                            compute_local_ordering(R,vs,es,search_ordering,search_small_to_large,0,max);
                        }
                    }

                    if (coloring_method == REPAIR_COLORING) {
                        greedy_neighborhood_coloring_repair_only(vs, es, R, C, C_max, colors, max, adj, is_enum);
                    }
                    else {
                        greedy_neighborhood_coloring(vs, es, R, C, C_max, colors, max, adj, is_enum);
                    }

                    S[tid][levels[tid]].inc_i1();
                    levels[tid]++;
                    clique_enum_dense_dynamic(vs, es, R, C, C_max, ind, colors, pruned, max, adj, S, levels, steps, cliques_local);
                    levels[tid]--;
                }
                // backtrack and search another branch
                R.clear();
                C.pop_back();
            }
            else return;
            P.pop_back();
        }
}

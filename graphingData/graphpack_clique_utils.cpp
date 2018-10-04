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

#include "graphpack_graph.h"
#include <algorithm>

using namespace std;
using namespace graphpack;

inline
void graphpack_graph::prune_vertex(int* &pruned, int i) {
    pruned[kcore_order[i]] = 1;
    for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
        adj[kcore_order[i]][edges[j]] = false;
        adj[edges[j]][kcore_order[i]] = false;
    }
}


int graphpack_graph::initial_pruning(graphpack_graph& G, int* &pruned, int lb, bool** &adj, string pruning_strat, bool is_enum) {

    printf("[graphpack: initial pruning] \n");
    int lb_idx = 0;
    if (pruning_strat == "kcore") {
        for (int i = G.num_vertices()-1; i >= 0; i--) {
            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
            printf("v=%d, K(v)=%d, lb=%d \n", i, kcore[kcore_order[i]], lb);
            if (is_enum) {
                if (kcore[kcore_order[i]] < lb) {
                    prune_vertex(pruned,i);
                }
            }
            else if (kcore[kcore_order[i]] <= lb) {
                prune_vertex(pruned,i);
            }
        }
    }
    else if (pruning_strat == "tcore") {
        vector<long long> V(vertices.size(),0);
        vector<int> E;
        E.reserve(edges.size());

        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (before)" <<endl;
        int start = 0;
        for (int i = 0; i < num_vertices(); i++) {
            start = E.size();
            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
                if (is_enum) {
                    if (tri_core[edges[j]]+2 >= lb)
                        E.push_back(edges[j]);
                }
                else if (tri_core[edges[j]]+2 > lb)
                    E.push_back(edges[j]);
            }
            V[i] = start;
            V[i + 1] = E.size();
        }
        vertices = V;
        edges = E;
        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (after)" <<endl;

        G.update_degrees();

        // now for k-core pruning
        G.compute_cores();
        cout << "K: " << get_max_core() <<endl;
        for (int i = G.num_vertices()-1; i >= 0; i--) {
            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
            if (is_enum) {
                if (kcore[kcore_order[i]] < lb) {
                    prune_vertex(pruned,i);
                }
            }
            else if (kcore[kcore_order[i]] <= lb) {
                prune_vertex(pruned,i);
            }
//            if (kcore[kcore_order[i]] <= lb) {
//                pruned[kcore_order[i]] = 1;
//                for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
//                    adj[kcore_order[i]][edges[j]] = false;
//                    adj[edges[j]][kcore_order[i]] = false;
//                }
//            }
        }

        //            G.degree_bucket_sort(true);

    }

    // reduces graph, updates degrees, and sorts edges
    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

int graphpack_graph::initial_pruning_induced(graphpack_graph& G, int* &pruned, int lb, string pruning_strat) {

    int lb_idx = 0;
    if (pruning_strat == "kcore") {
        for (int i = G.num_vertices()-1; i >= 0; i--) {
            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
            if (kcore[kcore_order[i]] <= lb)  pruned[kcore_order[i]] = 1;
        }
    }
    else if (pruning_strat == "tcore") {
        vector<long long> V(vertices.size(),0);
        vector<int> E;
        E.reserve(edges.size());

        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (before)" <<endl;
        int start = 0;
        for (int i = 0; i < num_vertices(); i++) {
            start = E.size();
            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
                if (tri_core[edges[j]]+2 >= lb)
                    E.push_back(edges[j]);
            }
            V[i] = start;
            V[i + 1] = E.size();
        }
        vertices = V;
        edges = E;
        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (after)" <<endl;

        G.update_degrees();

        // now for k-core pruning
//        G.compute_cores();
//        cout << "K: " << get_max_core() <<endl;
//        for (int i = G.num_vertices()-1; i >= 0; i--) {
//            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
//            if (kcore[kcore_order[i]] <= lb) {
//                pruned[kcore_order[i]] = 1;
//                for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
//                    adj[kcore_order[i]][edges[j]] = false;
//                    adj[edges[j]][kcore_order[i]] = false;
//                }
//            }
//        }
    }

    // reduces graph, updates degrees, and sorts edges
    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

int graphpack_graph::initial_pruning(graphpack_graph& G, int* &pruned, int lb, string pruning_strat, bool is_enum) {

    if (G.verbose)
        printf("[graphpack: initial pruning] \n");

    int lb_idx = 0;
    if (pruning_strat == "kcore") {
        for (int i = G.num_vertices()-1; i >= 0; i--) {
            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
            if (is_enum && kcore[kcore_order[i]] < lb)  pruned[kcore_order[i]] = 1;
            else if (kcore[kcore_order[i]] <= lb)  pruned[kcore_order[i]] = 1;
        }
    }
    else if (pruning_strat == "tcore") {
        vector<long long> V(vertices.size(),0);
        vector<int> E;
        E.reserve(edges.size());

        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (before)" <<endl;
        int start = 0;
        for (int i = 0; i < num_vertices(); i++) {
            start = E.size();
            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
                if (is_enum && tri_core[edges[j]]+2 >= lb)
                    E.push_back(edges[j]);
                else if (tri_core[edges[j]]+2 >= lb)
                    E.push_back(edges[j]);
            }
            V[i] = start;
            V[i + 1] = E.size();
        }
        vertices = V;
        edges = E;
        cout << "[TRIANGLE CORE PRUNING]  |E| = " << edges.size() << " (after)" <<endl;

        G.update_degrees();

        // now for k-core pruning
//        G.compute_cores();
//        cout << "K: " << get_max_core() <<endl;
//        for (int i = G.num_vertices()-1; i >= 0; i--) {
//            if (kcore[kcore_order[i]] == lb)  lb_idx = i;
//            if (kcore[kcore_order[i]] <= lb) {
//                pruned[kcore_order[i]] = 1;
//                for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
//                    adj[kcore_order[i]][edges[j]] = false;
//                    adj[edges[j]][kcore_order[i]] = false;
//                }
//            }
//        }
    }

    // reduces graph, updates degrees, and sorts edges
    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

int graphpack_graph::initial_pruning(graphpack_graph& G, int* &pruned, int lb) {
    printf("[graphpack: initial pruning] \n");

    int lb_idx = 0;
    for (int i = G.num_vertices()-1; i >= 0; i--) {
        if (kcore[kcore_order[i]] == lb)  lb_idx = i;
        if (kcore[kcore_order[i]] <= lb)  pruned[kcore_order[i]] = 1;
    }

    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

int graphpack_graph::initial_pruning(graphpack_graph& G, int* &pruned, int lb, bool** &adj) {
    int lb_idx = 0;
    for (int i = G.num_vertices()-1; i >= 0; i--) {
        if (kcore[kcore_order[i]] == lb)  lb_idx = i;
        if (kcore[kcore_order[i]] <= lb) {
            pruned[kcore_order[i]] = 1;
            for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
                adj[kcore_order[i]][edges[j]] = false;
                adj[edges[j]][kcore_order[i]] = false;
            }
        }
    }

    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

// TODO: remove these enum funcs, and simply call them above with lb-1
int graphpack_graph::initial_pruning_enum(graphpack_graph& G, int* &pruned, int lb) {
    int lb_idx = 0;
    for (int i = G.num_vertices()-1; i >= 0; i--) {
        if (kcore[kcore_order[i]] == lb)  lb_idx = i;
        if (kcore[kcore_order[i]] < lb)  pruned[kcore_order[i]] = 1;
    }

//    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

int graphpack_graph::initial_pruning_enum(graphpack_graph& G, int* &pruned, int lb, bool** &adj) {
    int lb_idx = 0;
    for (int i = G.num_vertices()-1; i >= 0; i--) {
        if (kcore[kcore_order[i]] == lb)  lb_idx = i;
        if (kcore[kcore_order[i]] < lb) {
            pruned[kcore_order[i]] = 1;
            for (long long j = vertices[kcore_order[i]]; j < vertices[kcore_order[i] + 1]; j++) {
                adj[kcore_order[i]][edges[j]] = false;
                adj[edges[j]][kcore_order[i]] = false;
            }
        }
    }

//    reduce_and_sort_edges(G, pruned, lb_idx);
    return lb_idx;
}

void graphpack_graph::reduce_and_sort_edges(graphpack_graph& G, int* &pruned, int lb_idx) {
    double sec = get_time();
    if (G.verbose) {
        cout << "[graphpack: initial k-core pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
    }
    G.reduce_graph(pruned);

    if (G.verbose) {
        cout << "[graphpack: initial k-core pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;
        cout << "[mcpack]  initial pruning took " << get_time()-sec << " sec" <<endl;
    }

    G.update_degrees();
//    G.degree_bucket_sort(true);
    G.degree_bucket_sort_parallel(true);
}


/*
 * Function modified when updating triangle cores
 */
void graphpack_graph::order_vertices_tcore(vector<Vertex> &V, graphpack_graph &G,
        int &lb_idx, int &lb, string vertex_ordering, bool decr_order,
        vector<int> & pruned_vertices) {

    cout << "[mcpack]  ORDERING VERTICES BY " << vertex_ordering <<endl;

    srand (time(NULL));
    int u = 0, val = 0;
    for (int u = lb_idx; u < G.num_vertices(); u++) {
        // edges and vertices from G may have been removed via kcore/triangle-core pruning
        if (degree[u] >= lb - 1) {

            if (vertex_ordering == "deg")
                val = vertices[u + 1] - vertices[u];
            else if (vertex_ordering == "kcore")
                val = kcore[u];
            else if (vertex_ordering == "kcore_deg")
                val = degree[u] * kcore[u];
            else if (vertex_ordering == "rand")
                val = rand() % vertices.size();
            // neighbor degrees
            else if (vertex_ordering == "dual_deg") {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    val = val + G.vertex_degree(edges[j]);
                }
            }
            // neighbor degrees
            else if (vertex_ordering == "dual_kcore") {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    val = val + kcore[edges[j]];
                }
            }
            else  val = vertices[u + 1] - vertices[u];
            V.push_back(Vertex(u,val));
        }
        else { // prune the vertex, note could use triangle core vertices here
            pruned_vertices[u] = 1;
        }
    }
    if (decr_order)
        std::sort(V.begin(), V.end(), decr_bound);
    else
        std::sort(V.begin(), V.end(), incr_bound);
}


/*
 * Original ordering function
 */
void graphpack_graph::order_vertices(vector<Vertex> &V, graphpack_graph &G,
        int &lb_idx, int &lb, string vertex_ordering, bool decr_order) {

    srand (time(NULL));
    int u = 0, val = 0;
    for (int k = lb_idx; k < G.num_vertices(); k++) {
        if (degree[kcore_order[k]] >= lb - 1) {
            u = kcore_order[k];

            if (vertex_ordering == "deg")
                val = vertices[u + 1] - vertices[u];
            else if (vertex_ordering == "kcore")
                val = kcore[u];
            else if (vertex_ordering == "kcore_deg")
                val = degree[u] * kcore[u];
            else if (vertex_ordering == "rand")
                val = rand() % vertices.size();
            // neighbor degrees
            else if (vertex_ordering == "dual_deg") {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    val = val + G.vertex_degree(edges[j]);
                }
            }
            // neighbor degrees
            else if (vertex_ordering == "dual_kcore") {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    val = val + kcore[edges[j]];
                }
            }
            else  val = vertices[u + 1] - vertices[u];
            V.push_back(Vertex(u,val));
        }
    }
    if (decr_order)
        std::sort(V.begin(), V.end(), decr_bound);
    else
        std::sort(V.begin(), V.end(), incr_bound);
}




/*
 * Simple but fast ordering function
 * Used in graphpack_heu_pmcx.cpp
 *
 * Last updated: 04/11/2013
 * Copyright Ryan A. Rossi
 */
void graphpack_graph::order_vertices(vector<Vertex> &V, graphpack_graph &G,
        string vertex_ordering, bool decr_order) {

    srand (time(NULL));
    int val = 0;
    for (int u = 0; u < G.num_vertices(); u++) {

        ordering_strategies(u,val,vertex_ordering);
        V.push_back(Vertex(u,val));
    }
    if (decr_order)
        std::sort(V.begin(), V.end(), decr_bound);
    else
        std::sort(V.begin(), V.end(), incr_bound);
}

inline
void graphpack_graph::ordering_strategies(int & u, int & val, string & vertex_ordering) {

    if (vertex_ordering == "deg")
        val = vertices[u + 1] - vertices[u];
    else if (vertex_ordering == "kcore")
        val = kcore[u];
    else if (vertex_ordering == "kcore_deg")
        val = degree[u] * kcore[u];
    else if (vertex_ordering == "rand")
        val = rand() % vertices.size();
    // neighbor degrees
    else if (vertex_ordering == "dual_deg") {
        val = 0;
        for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
            val = val + vertex_degree(edges[j]);
        }
    }
    // neighbor degrees
    else if (vertex_ordering == "dual_kcore") {
        val = 0;
        for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
            val = val + kcore[edges[j]];
        }
    }
    else  val = vertices[u + 1] - vertices[u];
}

void graphpack_graph::order_vertices(vector<Vertex> &V, string vertex_ordering, bool decr_order) {

    srand (time(NULL));
    int val = 0;
    for (int u = 0; u < num_vertices(); u++) {

        ordering_strategies(u,val,vertex_ordering);
        V.push_back(Vertex(u,val));
    }
    if (decr_order)
        std::sort(V.begin(), V.end(), decr_bound);
    else
        std::sort(V.begin(), V.end(), incr_bound);
}

void graphpack_graph::reduce_graph(
        vector<long long>& vs,
        vector<int>& es,
        int* &pruned,
        graphpack_graph& G,
        int id,
        int& mc) {

    int num_vs = vs.size();

    vector<long long> V(num_vs,0);
    vector<int> E;
    E.reserve(es.size());


    int start = 0;
    for (int i = 0; i < num_vs - 1; i++) {
        start = E.size();
        if (!pruned[i]) { //skip these V_local...
            for (long long j = vs[i]; j < vs[i + 1]; j++ ) {
                if (!pruned[es[j]])
                    E.push_back(es[j]);
            }
        }
        V[i] = start;
        V[i + 1] = E.size();
        G.degree[i] = V[i+1] - V[i];
    }
    vs = V;
    es = E;

    /**
     * compute k-cores and share bounds: ensure operation completed by single process
     * if bounds are currently being updated by another worker/processor, then skip it.
     */
    #pragma omp single nowait
    {
//        G.vertices = vs;
//        G.edges = es;

        cout << ">>> [graphpack: thread " << omp_get_thread_num() + 1 << "]  updating dynamic bounds" <<endl;

//        if (G.is_dense) {
        G.update_degrees(vs,es,mc,pruned);
//        }
        G.induced_cores_ordering(vs,es,pruned);
    }
    V.clear();
    E.clear();
}

void graphpack_graph::reduce_graph(
        vector<long long>& vs,
        vector<int>& es,
        vector<int> & pruned_vertices,
        graphpack_graph& G,
        int id,
        int& mc) {

    int num_vs = vs.size();
    vector<long long> V(num_vs,0);
    vector<int> E;
    E.reserve(es.size());

    int start = 0;
    for (int i = 0; i < num_vs - 1; i++) {
        start = E.size();
        if (!pruned_vertices[i]) { //skip these V_local...
            for (long long j = vs[i]; j < vs[i + 1]; j++ ) {
                if (!pruned_vertices[es[j]])
                    E.push_back(es[j]);
            }
        }
        V[i] = start;
        V[i + 1] = E.size();
    }
    vs = V;
    es = E;

//    // compute k-cores and share bounds: ensure operation completed by single process
    #pragma omp single nowait
    {
        G.edges = es;
        G.vertices = vs;
        cout << ">>> [graphpack: thread " << omp_get_thread_num() + 1 << "]" <<endl;
        G.induced_cores_ordering(vs,es,pruned_vertices);
    }
    V.clear();
    E.clear();
}

void graphpack_graph::reduce_graph(vector<int> & pruned_vertices, graphpack_graph& G) {

    int num_vs = G.vertices.size();

    vector<long long> V(num_vs,0);
    vector<int> E;
    E.reserve(G.edges.size());

    int start = 0;
    for (int i = 0; i < num_vs - 1; i++) {
        start = E.size();
        if (!pruned_vertices[i]) { //skip these V_local...
            for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {
                if (!pruned_vertices[G.edges[j]])
                    E.push_back(G.edges[j]);
            }
        }
        V[i] = start;
        V[i + 1] = E.size();
        G.degree[i] = V[i+1] - V[i];
    }
    G.vertices = V;
    G.edges = E;

//    // compute k-cores and share bounds: ensure operation completed by single process
//    #pragma omp single nowait
//    {
//        G.edges = es;
//        G.vertices = vs;
//        cout << ">>> [graphpack: thread " << omp_get_thread_num() + 1 << "]" <<endl;
//        G.induced_cores_ordering(vs,es,pruned_vertices);
//    }
//    V.clear();
//    E.clear();
}

void graphpack_graph::print_info(vector<int> &C_max, double &sec) {
    cout << "*** [graphpack: thread " << omp_get_thread_num() + 1;
    cout << "]   current max clique = " << C_max.size();
    cout << ",  time = " << get_time() - sec << " sec" <<endl;
}

void graphpack_graph::print_break() {
    if (verbose) {
    cout << "-----------------------------------------------------------------------" <<endl;
    }
}

bool graphpack_graph::time_left(vector<int> &C_max, double sec, double time_limit, bool &time_expired_msg) {
    if ((get_time() - sec) > time_limit) {
        if (time_expired_msg) {
            cout << "\n### Time limit expired, terminating search. ###" <<endl;
            cout << "Size: " << C_max.size() <<endl;
            print_max_clique(C_max);
            time_expired_msg = false;
        }
        return false;
    }
    return true;
}

void graphpack_graph::graph_stats(graphpack_graph& G, vector<long long> & vs, vector<int> & es,
        int& mc, int id, double &sec) {
    cout << "[graphpack: bounds updated - thread " << omp_get_thread_num() + 1 << "]  ";
    cout << "time = " << get_time() - sec << " sec, ";
    cout << "|V| = " << (G.num_vertices() - id);
    cout << " (" << id << " / " << G.num_vertices();
    cout << "), |E| = " << es.size()/2;
    cout << ", w = " << mc;
    cout << ", p = " << G.density(es.size()/2,G.num_vertices() - id);
//    cout << ", d_min = " << G.get_min_degree();
    cout << ", d_avg = " << G.get_avg_degree();
    cout << ", d_max = " << G.get_max_degree();
    cout << ", k_max = " << G.get_max_core();
    cout <<endl;
}

void graphpack_graph::graph_stats(graphpack_graph& G, int& mc, int id, double &sec) {
    if (G.verbose) {
        cout << "[graphpack: bounds updated - thread " << omp_get_thread_num() + 1 << "]  ";
        cout << "time = " << get_time() - sec << " sec, ";
        cout << "|V| = " << (G.num_vertices() - id);
        cout << " (" << id << " / " << G.num_vertices();
        cout << "), |E| = " << G.num_edges();
        cout << ", w = " << mc;
        cout << ", p = " << G.density();
        cout << ", d_min = " << G.get_min_degree();
        cout << ", d_avg = " << G.get_avg_degree();
        cout << ", d_max = " << G.get_max_degree();
        cout << ", k_max = " << G.get_max_core();
        cout <<endl;
    }
}

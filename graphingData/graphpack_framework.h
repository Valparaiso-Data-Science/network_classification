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

#ifndef MCPACK_FRAMEWORK_H_
#define MCPACK_FRAMEWORK_H_

#include "graphpack_vertex.h"
#include "graphpack_neigh_coloring.h"
#include "graphpack_neigh_cores.h"
#include <algorithm>
#include <set>

using namespace std;

namespace graphpack {

    enum  { // coloring methods
        REPAIR_COLORING,            // 0
        BASIC_COLORING,             // 1
        SPARSE_COLORING,            // 2
    };

    enum  { // order and pruning methods
        DEGREE,             // 0
        KCORE,              // 1
        KCORE_DEG,          // 2
        TRIANGLES,          // 3
        WEDGES,             // 4
        RAND,               // 5
        VAR,                // 6
        TRIANGLE_CORES,     // 7
        COLORING,           // 8
        DEGREE_VOL,         // 9
        KCORE_VOL,          // 10
        TRIANGLE_VOL,       // 11
        TRIANGLE_CORE_VOL,  // 12
        TRIANGLES_ONLY,     // 13 (for global pruning)
        NATURAL,            // 14
        DYNAMIC_LARGEST_FIRST,          // 15
        // orderings for greedy coloring
        TRIANGLE_CORE_MAX,              // 16
        DEGREE_TRIANGLES,               // 17
        KCORE_TRIANGLES,                // 18
        KCORE_DEG_TRI,                  // 19
        DEGREE_KCORE_VOL,               // 20
        KCORE_TRIANGLE_VOL,             // 21
        DEGREE_KCORE_TRIANGLE_VOL,      // 22
        DIST_TWO_LARGEST_FIRST,         // 23
        DIST_TWO_DYNAMIC_LARGEST_FIRST, // 24
        DIST_TWO_SMALLEST_LAST,         // 25
        INCIDENCE_DEGREE,               // 26
        DIST_TWO_INCIDENCE_DEGREE,      // 27
        TRIANGLE_CORE_DIST_TWO,         // 28
        TRIANGLE_CORE_MIN,              // 29
    };

    enum  { // representation
        CSC,         // 0
        HYBRID,      // 1
        ADJ,         // 2
    };

    enum  { // search methods
        REPAIR,             // 0
        STATIC_ORDER,       // 1
        REPAIR_ONLY,        // 2
        KCORE_SEARCH,       // 3
        DYNAMIC_COLORING,   // 4
        BASIC,              // 5
        DYNAMIC_KCORE,      // 6
        DYNAMIC,            // 7
    };

    static int set_enum_search_method(string &method) {
        if      (method == "repair" || method == "repair_static_order")                             return 0;
        else if (method == "static_order" || method == "static_ordering" || method == "static")     return 1;
        else if (method == "repair_only" || method == "repair_coloring_only")                       return 2;
        else if (method == "kcores" || method == "cores" || method == "kcore")                      return 3;
        else if (method == "coloring" || method == "greedy_coloring")                               return 4;
        else if (method == "basic" || method == "normal" || method == "original")                   return 5;
        else if (method == "dynamic_kcore" || method == "dynamic_cores")                            return 6;
        else if (method == "dynamic_bounds" || method == "dynamic")                                 return 7;
        return 0;
    }

    static int set_enum_representation(string &rep) {
        if      (rep == "csc" || rep == "sparse")                    return 0;
        else if (rep == "hybrid" || rep == "csc_adj")                  return 1;
        else if (rep == "adj" || rep == "adjcency")        return 2;
        return 0;
    }

    static int set_enum(string &ordering) {
        if      (ordering == "deg" || ordering == "degree")                                 return 0;
        else if (ordering == "kcore" || ordering == "kcores")                               return 1;
        else if (ordering == "kcore_deg" || ordering == "kcore_degree")                     return 2;
        else if (ordering == "triangles" || ordering == "tri" || ordering == "triangle")    return 3;
        else if (ordering == "wedges")                                                      return 4;
        else if (ordering == "rand" || ordering == "random")                                return 5;
        else if (ordering == "var")                                                         return 6;
        else if (ordering == "tcores" || ordering == "triangle_cores" ||
                ordering == "triangle_core" || ordering == "tcore")                         return 7;
        else if (ordering == "coloring" || ordering == "greedy_coloring")                   return 8;
        else if (ordering == "deg_vol" || ordering == "degree_vol")                         return 9;
        else if (ordering == "kcore_vol" || ordering == "kcore_volume")                     return 10;
        else if (ordering == "tri_vol" || ordering == "triangle_vol")                       return 11;
        else if (ordering == "tcore_vol" || ordering == "tcore_volume" ||
                ordering == "triangle_core_vol")                                            return 12;
        else if (ordering == "triangles_only" || ordering == "triangle_only")               return 13;
        else if (ordering == "natural")                                                     return 14;
        else if (ordering == "dynamic_largest_first" || ordering == "lfo" || ordering =="LFO")  return 15; // graph coloring
        else if (ordering == "triangle_core_max" || ordering == "tcore_max" ||
                ordering == "tcores_max")                                                   return 16;

        else if (ordering == "degree_triangles" || ordering == "degree_tri" ||
                ordering == "deg_tri" || ordering == "deg_triangles")                       return 17;

        else if (ordering == "kcore_tri" || ordering == "kcore_triangles" ||
                ordering =="kcore_triangle")                                                return 18; // graph coloring


        else if (ordering == "kcore_deg_tri" || ordering == "kcore_degree_triangle" ||
                ordering == "kcore_degree_tri")                                             return 19;

        else if (ordering == "kcore_deg_vol" || ordering == "kcore_degree_vol")             return 20;

        else if (ordering == "kcore_tri_vol" || ordering == "kcore_triangle_vol")           return 21;
        else if (ordering == "kcore_deg_tri_vol" || ordering == "kcore_degree_tri_vol")     return 22;
        else if (ordering == "dist_two_lfo" || ordering == "dist_two_largest_first")        return 23;
        else if (ordering == "dist_two_dlfo" || ordering == "dist_two_dynamic_lfo")         return 24;
        else if (ordering == "dist_two_slo" || ordering == "dist_two_smallest_last")        return 25;
        else if (ordering == "ido" || ordering == "incidence" ||
                ordering == "incidence_deg")                                                return 26;
        else if (ordering == "dist_two_ido" || ordering == "dist_two_incidence")            return 27;
        else if (ordering == "dist_two_tcore" || ordering == "dist_two_triangle_core")      return 28;


        else if (ordering == "tcore_min" || ordering == "triangle_core_min" ||
                ordering =="tri_core_min")                                                  return 29; // graph coloring
        else return 1; // default kcore
    }

    static int set_ordering_enum(string &ordering) {
        if      (ordering == "deg" || ordering == "degree")                    return 0;
        else if (ordering == "kcore" || ordering == "kcores")                  return 1;
        else if (ordering == "kcore_deg" || ordering == "kcore_degree")        return 2;
        else if (ordering == "triangles" || ordering == "tri")                 return 3;
        else if (ordering == "wedges")                                         return 4;
        else if (ordering == "rand" || ordering == "random")                   return 5;
        else if (ordering == "var" || ordering == "variance")                  return 6;
        else if (ordering == "tcores" || ordering == "triangle_cores")         return 7;
        else if (ordering == "coloring" || ordering == "greedy_coloring")      return 8;
        else if (ordering == "deg_vol" || ordering == "degree_vol")            return 9;
        else if (ordering == "kcore_vol" || ordering == "kcore_volume")        return 10;
        else if (ordering == "tri_vol" || ordering == "triangle_vol")          return 11;
        else return 1; // default kcore
    }

    static int set_enum_coloring_method(string &coloring_method) {
        if      (coloring_method == "repair" || coloring_method == "repair_colors")                                         return 0;
        else if (coloring_method == "normal" || coloring_method == "basic" || coloring_method == "original")                return 1;
        else return 1; // default kcore
    }


    inline
    static void ensure_graph_parameters(graphpack_graph &G, int ordering, int is_parallel = false) {
        if (ordering == DEGREE && G.degree.size() == 0) {
            G.vertex_degrees();
        }
        else if (ordering == KCORE && G.kcore.size() == 0) {
            printf("[graphpack: ordering]  kcores do not exist, computing...\n");
            G.compute_cores();
        }
        else if (ordering == TRIANGLES && G.t.size() == 0) {
            printf("[graphpack: ordering]  triangle counts do not exist, computing...\n");
            G.compute_vertex_triangles(0,is_parallel);
        }
        else if (ordering == TRIANGLE_CORES && G.tri_core.size() == 0) {
            printf("[graphpack: ordering]  triangle cores do not exist, computing...\n");
            G.compute_triangle_cores();
        }
    }

    // assumes V is empty, since vertices are pushed into V
    // 1,...,|V| vertices are in order, where |V| = G.num_vertices(),
    // thus if some vertices are implicitly pruned from G, then don't use this function
    // since G is used, assumes to be serial, unless G is a neighborhood local to a worker thread

    inline
    static void set_graph_parameters(vector<Vertex> &V, graphpack_graph &G, int ordering) {


        srand (time(NULL));
        int val = 0;
        for (int u = 0; u < G.num_vertices(); u++) {

            switch (ordering) {
                case NATURAL: { // natural order (read from disk)
                    val = u;
                    break;
                }
                case DEGREE: {
                    val = G.vertices[u + 1] - G.vertices[u];
                    break;
                }
                case WEDGES: {
                    val = G.vertices[u + 1] - G.vertices[u]; //store degree in val
                    val = (val*(val-1))/2;
                    break;
                }
                case KCORE: {
                    val = G.kcore[u];
                    break;
                }
                case TRIANGLES: {
                    val = G.t[u];
                    break;
                }
                case RAND: {
//                    val = drand() % G.vertices.size();
                    val = rand() % G.vertices.size();
                    break;
                }
                case VAR: {
                    val = G.kcore[u] * ((int)G.degree[u]/G.kcore[u]);
                    break;
                }
                case DEGREE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]];//G.vertex_degree(edges[j]);
                    }
                    break;
                }
                case KCORE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + G.kcore[G.edges[j]];
                    }
                    break;
                }
                case TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + G.t[G.edges[j]];
                    }
                    break;
                }
                case TRIANGLE_CORE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                    break;
                }
                case TRIANGLE_CORE_MAX: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        if (G.tri_core[G.eid[j]] > val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORES: { // sum
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                    val = val / ((double)(G.vertices[u+1] - G.vertices[u]));
                    break;
                }
                case TRIANGLE_CORE_MIN: {
                    val = G.num_vertices();
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        if (G.tri_core[G.eid[j]] < val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORE_DIST_TWO: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val += G.tri_core[G.eid[j]];
                        int w = G.edges[j];
                        for (long long k = G.vertices[w]; k < G.vertices[w + 1]; k++) {
                                val += G.tri_core[G.eid[k]]; //decay weight
                        }
                    }
                    break;
                }
                case KCORE_DEG: {
                    val = G.degree[u] * G.kcore[u];
                    break;
                }
                case DEGREE_TRIANGLES: {
                    val = G.degree[u] * G.t[u];
                    break;
                }
                case KCORE_TRIANGLES: {
                    val = G.kcore[u] * G.t[u];
                    break;
                }
                case KCORE_DEG_TRI: {
                    val = G.degree[u] * G.kcore[u] * G.t[u];
                    break;
                }
                case DEGREE_KCORE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + ((G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]]) * G.kcore[G.edges[j]]);
                    }
                    break;
                }
                case KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + (G.kcore[G.edges[j]] * G.t[G.edges[j]]);
                    }
                    break;
                }
                case DEGREE_KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                        val = val + (G.degree[G.edges[j]] * G.kcore[G.edges[j]] * G.t[G.edges[j]]);
                    }
                    break;
                }
                default: {
                    val = G.vertices[u + 1] - G.vertices[u];
                    break;
                }
            }
            V.push_back(Vertex(u,val));

        }
    }



    // vertices are processed 1 to |V|
    // u is the vertex at the ith position in V
    // val is added to the ith position in V
    inline
    static void set_graph_parameters(
            vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> & vertices,
            vector<int> & edges,
            int ordering) {


        srand (time(NULL));
        int val = 0;
        int u = 0;
        for (int i = 0; i < V.size(); i++) {

            u = V[i].get_id();

            switch (ordering) {
                case NATURAL: { // natural order (read from disk)
                    val = u;
                    break;
                }
                case DEGREE: {
                    val = vertices[u + 1] - vertices[u];
                    break;
                }
                case WEDGES: {
                    val = vertices[u + 1] - vertices[u]; //store degree in val
                    val = (val*(val-1))/2;
                    break;
                }
                case KCORE: {
                    val = G.kcore[u];
                    break;
                }
                case TRIANGLES: {
                    val = G.t[u];
                    break;
                }
                case RAND: {
//                    val = drand() % G.vertices.size();
                    val = rand() % vertices.size();
                    break;
                }
                case DEGREE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + vertices[edges[j] + 1] - vertices[edges[j]];//G.vertex_degree(edges[j]);
                    }
                    break;
                }
                case KCORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + G.kcore[edges[j]];
                    }
                    break;
                }
                case TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + G.t[edges[j]];
                    }
                    break;
                }
                case TRIANGLE_CORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                    break;
                }
                case TRIANGLE_CORE_MAX: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (G.tri_core[G.eid[j]] > val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORES: { // sum
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                    val = val / ((double)(vertices[u+1] - vertices[u]));
                    break;
                }
                case TRIANGLE_CORE_MIN: {
                    val = G.num_vertices();
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (G.tri_core[G.eid[j]] < val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORE_DIST_TWO: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val += G.tri_core[G.eid[j]];
                        int w = edges[j];
                        for (long long k = vertices[w]; k < vertices[w + 1]; k++) {
                                val += G.tri_core[G.eid[k]]; //decay weight
                        }
                    }
                    break;
                }
                case KCORE_DEG: {
                    val = G.degree[u] * G.kcore[u];
                    break;
                }
                case DEGREE_TRIANGLES: {
                    val = G.degree[u] * G.t[u];
                    break;
                }
                case KCORE_TRIANGLES: {
                    val = G.kcore[u] * G.t[u];
                    break;
                }
                case KCORE_DEG_TRI: {
                    val = G.degree[u] * G.kcore[u] * G.t[u];
                    break;
                }
                case DEGREE_KCORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + ((vertices[edges[j] + 1] - vertices[edges[j]]) * G.kcore[edges[j]]);
                    }
                    break;
                }
                case KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + (G.kcore[edges[j]] * G.t[edges[j]]);
                    }
                    break;
                }
                case DEGREE_KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        val = val + (G.degree[edges[j]] * G.kcore[edges[j]] * G.t[edges[j]]);
                    }
                    break;
                }
                default: {
                    val = vertices[u + 1] - vertices[u];
                    break;
                }
            }
            V[i].set_bound(val);
//            V.push_back(Vertex(u,val));

        }
    }



    template<typename T>
    static void neigh_cores_no_pruning(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            std::vector<T>& ind,
            int ordering) {

        int n = P.size() + 1;

        // lookup table
        vector<int> newids_to_actual(n, 0);
        vector<int> vert_order(n,0);
        vector<int> deg(n,0);
        vector<int> static_degree(n,0);
        vector<int> pos(n,0);

//        // lookup table for neighbors
//        for (int v = 1; v < n; v++) ind[P[v-1].get_id()] = 1;

        // compute degrees of induced neighborhood
        int md = 0, x, u;
        for (int v = 1; v < n; v++) {   // for each v in P
            u = P[v-1].get_id();
            x = 0;
            for (long long j=vs[u]; j<vs[u+1]; j++) { //induced degree
                if (ind[es[j]]) x++;
            }
            deg[v] = x;
            if (deg[v] > md)  md = deg[v];
        }
        static_degree = deg;

        int md_end = md+1;
        vector<int> bin(md_end,0);
        for (int v = 1; v < n; v++) bin[deg[v]]++;

        int start = 1, num = 0;
        for (int d=0; d < md_end; d++) { //for each deg, set bin to be the pos of the first vertex of that degree
            num = bin[d];
            bin[d] = start;
            start = start + num;
        }


        for (int v=1; v<n; v++) {
            pos[v] = bin[deg[v]];

            //view this step as relabeling the vertices
            vert_order[pos[v]] = v;
            ind[P[v-1].get_id()] = v;               // set bit for actual vertex id
            newids_to_actual[v] = P[v-1].get_id();
            bin[deg[v]]++;
        }

        for (int d=md; d > 1; d--)  bin[d] = bin[d-1];
        bin[0] = 1;


        int v_newid, v_actual, u_newid, du, pu, pw, w;
        long long j = 0;
        for (int i = 1; i < n; i++) {                           // neighborhood K-cores
            v_newid = vert_order[i];                            //relabeled id
            v_actual = newids_to_actual[v_newid];               // real id
            for (j = vs[v_actual]; j<vs[v_actual+1]; j++) {
                if (ind[es[j]] > 0) {                           // find common induced neighbors of k

                    u_newid = ind[es[j]];
                    if (deg[u_newid] > deg[v_newid]) {
                        du = deg[u_newid];
                        pu = pos[u_newid];
                        pw = bin[du];
                        w = vert_order[pw];
                        if (u_newid != w) {
                            pos[u_newid] = pw;
                            vert_order[pu] = w;
                            pos[w] = pu;
                            vert_order[pw] = u_newid;
                        }
                        bin[du]++;   deg[u_newid]--;
                    }
                }
            }
        }




        // reset index
//        for (int v=1; v<n; v++) ind[P[v-1].get_id()] = 0;


        // neighborhood core pruning and sorting
        int id = 0, prune_vert = 0;
        int v = 0;
        for (int i = n-1; i > 0; --i) {
            u = vert_order[i];
            v_actual = newids_to_actual[u];

            switch(ordering) {
                case KCORE_VOL: {
                    ind[v_actual] = deg[u];
                    P[id].set_bound(deg[u]);
                    P[id].set_id(v_actual);
                    break;
                }
                case DEGREE_KCORE_VOL: {
                    P[id].set_bound(deg[u]);
                    P[id].set_id(v_actual);
                    ind[v_actual] = deg[u] * static_degree[u];
                    break;
                }
                case KCORE_DEG: {
                    ind[v_actual] = deg[u] * static_degree[u];
                    P[id].set_bound(deg[u] * static_degree[u]);
                    P[id].set_id(v_actual);
                    break;
                }
                default: {
                    P[id].set_bound(deg[u]);
                    P[id].set_id(v_actual);
                    break;
                }
            }


            id++;
        }
//
//        // remove pruned verts from P
//        for (int i = 0; i < prune_vert; i++)
//            P.pop_back();
    }


    inline
    static void vertex_degree_among_set(
            vector<Vertex> &P,
            graphpack_graph &G,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind) {

        // if |P| > d_avg and adj matrix is stored, then faster to check P
//        if (G.avg_degree > P.size()) {
//
//        }
//        else {
//
//        }


        int u = 0;
        int deg = 0;
        for (int i = 0; i < P.size(); i++) {
            // get a vertex
            u = P[i].get_id();
            deg = 0;

            for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                if (ind[edges[j]]) deg++;
            }
            P[i].set_bound(i);
        }
    }


    inline
    static void local_vertex_degree(
            int & u,
            vector<Vertex> &P,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            int & val) {

        val = 0;
        for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
            if (ind[edges[j]]) val++;
        }

    }



    /**
     * Compute local triangle count for a particular vertex in an induced set P
     * @param v
     * @param P
     * @param vertices
     * @param edges
     * @param ind
     * @param marked
     * @param val
     */
    inline
    static void local_vertex_triangle_count(
            int & v,
            int & v_idx,
            vector<Vertex> &P,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            vector<int> & marked,
            int & val) {

        // create lookup table on neighbors of u
        for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
            // avoid reset by assigning marked[v] = u, then
            // use: marked[w]-1 == v, to check for neighbor
            marked[edges[j]] = v+1;
        }

        int u = 0, w = 0;
        int deg = 0;
        val = 0;

        // neighbors of v
        for (long long i = vertices[v]; i < vertices[v + 1]; i++) {

            // valid neighbor of v defined as u \in N(v), and u \in P
            if (ind[edges[i]]) {
                u = edges[i];
                deg++; // update degree

                // neighbors of u
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {

                    if (edges[j] == v) continue; // w is actually v
                    // is w a valid neighbor of u
                    if (ind[edges[j]]) {
//                        w = edges[j];

                        // triangle check: if w is a neighbor of v and valid (w must be valid from above)
                        if (marked[edges[j]]-1 == v) {
                            val++; // update triangle counts
                        }
                    } // if w is a valid neighbor of u
                } // neighbors of u
            } // if u is a valid neighbor of v
        } // neighbors of v

        // store degree back to P for free!
        P[v_idx].set_bound(val);
    }

    inline
    static void vertex_triangle_count_vertex_set(
//            graphpack_graph &G,
            vector<Vertex> &P,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            vector<int> & marked,
            int & ordering) {

        vector<int> N, T;
        N.reserve(P.size());
        T.reserve(P.size());


        int u = 0, w = 0;
        int v = 0;

        int num_triangles = 0;

        for (int i = 0; i < P.size(); i++) {


            v = P[i].get_id();

            // N consists of the neighbors of v that are also in P
            // create lookup table on neighbors of u
            for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
                // avoid reset by assigning marked[v] = u, then
                // use: marked[w]-1 == v, to check for neighbor
                if (ind[edges[j]]) {
                    marked[edges[j]] = v+1;
                    N.push_back(edges[j]);
                }
            }

            // degree of v
//            P[i].set_bound(N.size());

            num_triangles = 0;

            // neighbors of v
            //            for (long long i = vertices[v]; i < vertices[v + 1]; i++) {
            for (int k = 0; k < N.size(); ++k) {

                u = N[k]; // avoid checking ind for valid neighbor

                // neighbors of u
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {

                    if (edges[j] == v) continue; // w is actually v
                    // is w a valid neighbor of u
                    if (ind[edges[j]]) {
                        //                        w = edges[j];

                        // triangle check: if w is a neighbor of v and valid (w must be valid from above)
                        if (marked[edges[j]]-1 == v) {
                            num_triangles++; // update triangle counts
                        }
                    } // if w is a valid neighbor of u
                } // neighbors of u
            } // neighbors of v

            // store degree back to P for free!
            //            P[v_idx].set_bound(deg);



            switch (ordering) {

                case DEGREE_TRIANGLES: {
//                    ind[v] = N.size() * num_triangles;
                    P[i].set_bound(N.size() * num_triangles);
                    break;
                }

                case KCORE_TRIANGLES: {
                    // k-cores is assumed to be stored in P
//                    ind[v] = P[i].get_bound() * num_triangles;
                    P[i].set_bound(P[i].get_bound() * num_triangles);
                     break;
                 }
                case KCORE_DEG_TRI: {
//                    ind[v] = P[i].get_bound() * num_triangles * N.size();
                    P[i].set_bound(P[i].get_bound() * num_triangles * N.size());
//                       val = G.degree[u] * G.kcore[u] * G.t[u];
                       break;
                   }
                case TRIANGLE_VOL: {
                    ind[v] = num_triangles;
//                    P[i].set_bound(num_triangles);
                    break;
                }
                case KCORE_TRIANGLE_VOL: {
                    ind[v] = P[i].get_bound() * num_triangles;
//                    P[i].set_bound(P[i].get_bound() * num_triangles);
                    break;
                }
                case DEGREE_KCORE_TRIANGLE_VOL: {
                    ind[v] = P[i].get_bound() * num_triangles * N.size();
//                    P[i].set_bound(P[i].get_bound() * num_triangles * N.size());
                    break;
                }
            }

            P[i].set_bound(N.size());

            N = T; // reset array
        }
    }


    inline
    static void degree_vertex_set(
//            graphpack_graph &G,
            vector<Vertex> &P,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind) {

        int v = 0;
        int deg = 0;

        for (int i = 0; i < P.size(); i++) {

            v = P[i].get_id();

            // N consists of the neighbors of v that are also in P
            // create lookup table on neighbors of u
            for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
                // avoid reset by assigning marked[v] = u, then
                // use: marked[w]-1 == v, to check for neighbor
                if (ind[edges[j]]) {
                    deg++;
                }
            }
            ind[v] = deg;
            P[i].set_bound(deg);
        }
    }



    inline
    static void egonet_computation(
            vector<Vertex> &P,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind) {

        // assumes values are stored in ind

        int u = 0, val = 0;
        for (int i = 0; i < P.size(); ++i) {
            u = P[i].get_id();
            val = 0;
            for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                val += ind[edges[j]]; // if invalid neighbors, then zero!
            }
            P[i].set_bound(val);
        }

    }


    // used in exact
    inline
    static void set_graph_parameters_exact_simple(
            vector<Vertex> &P,
            graphpack_graph &G,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            vector<int> & marked,
            int ordering) {

        srand (time(NULL));
        int val = 0;
        int u = 0;
        for (int i = 0; i < P.size(); i++) {

            u = P[i].get_id();

            switch (ordering) {
                case NATURAL: { // natural order (read from disk)
                    val = u;
                    break;
                }
                case RAND: {
                    //                    val = drand() % G.vertices.size();
                    val = rand() % vertices.size();
                    break;
                }
                case DEGREE: {
                    // stores result in val
                    local_vertex_degree(u,P,vertices,edges,ind,val);
                    break;
                }
                case WEDGES: {
                    local_vertex_degree(u,P,vertices,edges,ind,val);
                    //                    val = vertices[u + 1] - vertices[u]; //store degree in val
                    val = (val*(val-1))/2;
                    break;
                }
                case TRIANGLES: {
                    //
                    local_vertex_triangle_count(u,i,P,vertices,edges,ind,marked,val);
                    //                    val = G.t[u];
                    break;
                }
                default: {
                    val = vertices[u + 1] - vertices[u];
                    break;
                }
            }
            P[i].set_bound(val);
        }
    }


    /**
     * see greedy_coloring in graphpack_mc.cpp
     *
     * Only difference from above, is that for the egonet/vol measures,
     * the measures are computed strictly over the vertices in P AND are EXACT***,
     * and not over G via vertices/edges array!
     *
     * The vertices in P are marked in the ind array.
     *
     * Note that for degree and those other parameters, we use the degree in the entire graph for efficientcy
     *
     * Degree properties can be computed from the vertices/edges array
     * For k-core, triangles, triangle-cores, these are stored in G.
     * However, another function may store them in P itself.
     *
     * 1. Compute properties explicitly over a vertex set P *exactly*, for this, we must
     * 2. Create induced graph H explicitly using P, then compute properties over it, and store them in H itself.
     * 3. Compute properties fast approximately over a vertex set P, using possibly stale values from G and vertices/edges for degree-related properties
     *
     * @param V
     * @param G
     * @param vertices
     * @param edges
     * @param ordering
     */
    inline
    static void set_graph_parameters_exact(
            vector<Vertex> &P,
            graphpack_graph &G,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            vector<int> & marked,
            int ordering) {

//        neigh_cores_bound(vs,es,P,ind,max);

//        if (ordering == KCORE || ordering == KCORE_VOL || ordering == KCORE_DEG || ordering == DEGREE_KCORE_VOL) {
//            neigh_cores_bound(vs,es,P,ind,marked);
//
//            if (ordering == KCORE) { // kcore number stored back into P, so can avoid
//                return;
//            }
//        }
//
//        if (ordering == DEGREE_VOL) {
//            // must first compute degree for each vertex and store it in P
//        }
//
//        if (ordering == TRIANGLE_VOL) {
//            // must first compute triangles for each vertex, store it in P
//        }

        int val = 0, u;

        switch(ordering) {

            case KCORE: {
                int max = 0;
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                break;
            }
            case KCORE_DEG: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                break;
            }
            case DEGREE_TRIANGLES: {
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
                break;
            }
            case KCORE_TRIANGLES: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
//                val = G.kcore[u] * G.t[u];
                break;
            }
            case KCORE_DEG_TRI: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
                break;
            }


            case DEGREE_VOL: {
                // puts degree in ind
                degree_vertex_set(P,vertices,edges,ind);
                egonet_computation(P,vertices,edges,ind);
                break;
            }
            case KCORE_VOL: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                egonet_computation(P,vertices,edges,ind);
                break;
            }
            case TRIANGLE_VOL: {
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
                egonet_computation(P,vertices,edges,ind);
                break;
            }
            case DEGREE_KCORE_VOL: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                egonet_computation(P,vertices,edges,ind);
                break;
            }
            case KCORE_TRIANGLE_VOL: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
                egonet_computation(P,vertices,edges,ind);

                break;
            }
            case DEGREE_KCORE_TRIANGLE_VOL: {
                neigh_cores_no_pruning(vertices,edges,P,ind,ordering);
                vertex_triangle_count_vertex_set(P,vertices,edges,ind,marked,ordering);
                egonet_computation(P,vertices,edges,ind);
                break;
            }



            /**
             * TRIANGLE CORE ORDERINGS
             **                  */
            case TRIANGLE_CORE_VOL: {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    if (ind[edges[j]]) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                }
                break;
            }
            case TRIANGLE_CORE_MAX: {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    if (ind[edges[j]]) {
                        if (G.tri_core[G.eid[j]] > val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                }
                break;
            }
            case TRIANGLE_CORES: { // sum
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    if (ind[edges[j]]) {
                        val = val + G.tri_core[G.eid[j]];
                    }
                }
                val = val / ((double)(vertices[u+1] - vertices[u]));
                break;
            }
            case TRIANGLE_CORE_MIN: {
                val = G.num_vertices();
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    if (ind[edges[j]]) {
                        if (G.tri_core[G.eid[j]] < val) {
                            val = G.tri_core[G.eid[j]];
                        }
                    }
                }
                break;
            }
            case TRIANGLE_CORE_DIST_TWO: {
                val = 0;
                for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                    if (ind[edges[j]]) {
                        val += G.tri_core[G.eid[j]];
                        int w = edges[j];
                        for (long long k = vertices[w]; k < vertices[w + 1]; k++) {
                            if (ind[edges[k]]) {
                                val += G.tri_core[G.eid[k]]; //decay weight
                            }
                        }
                    }
                }
                break;
            }
            default: {
                // sets all other orderings, in a similar fashion
                set_graph_parameters_exact_simple(P,G,vertices,edges,ind,marked,ordering);
                break;
            }
        } // end switch
    } // end function


    /**
     * see greedy_coloring in graphpack_mc.cpp
     *
     * Only difference from above, is that for the egonet/vol measures,
     * the measures are computed strictly over the vertices in P,
     * and not over G via vertices/edges array!
     *
     * The vertices in P are marked in the ind array.
     *
     * Note that for degree and those other parameters, we use the degree in the entire graph for efficientcy
     *
     * Degree properties can be computed from the vertices/edges array
     * For k-core, triangles, triangle-cores, these are stored in G.
     * However, another function may store them in P itself.
     *
     * 1. Compute properties explicitly over a vertex set P *exactly*, for this, we must
     * 2. Create induced graph H explicitly using P, then compute properties over it, and store them in H itself.
     * 3. Compute properties fast approximately over a vertex set P, using possibly stale values from G and vertices/edges for degree-related properties
     *
     * @param V
     * @param G
     * @param vertices
     * @param edges
     * @param ordering
     */
    inline
    static void set_graph_parameters_induced(
            vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> & vertices,
            vector<int> & edges,
            vector<int> & ind,
            int ordering) {


        srand (time(NULL));
        int val = 0;
        int u = 0;
        for (int i = 0; i < V.size(); i++) {

            u = V[i].get_id();

            switch (ordering) {
                case NATURAL: { // natural order (read from disk)
                    val = u;
                    break;
                }
                case RAND: {
//                    val = drand() % G.vertices.size();
                    val = rand() % vertices.size();
                    break;
                }
                case DEGREE: {
                    val = vertices[u + 1] - vertices[u];
                    break;
                }
                case WEDGES: {
                    val = vertices[u + 1] - vertices[u]; //store degree in val
                    val = (val*(val-1))/2;
                    break;
                }
                case KCORE: {
                    val = G.kcore[u];
                    break;
                }
                case TRIANGLES: {
                    val = G.t[u];
                    break;
                }

                case DEGREE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + vertices[edges[j] + 1] - vertices[edges[j]];//G.vertex_degree(edges[j]);
                        }
                    }
                    break;
                }
                case KCORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + G.kcore[edges[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + G.t[edges[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + G.tri_core[G.eid[j]];
                        }
                    }
                    break;
                }
                case TRIANGLE_CORE_MAX: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            if (G.tri_core[G.eid[j]] > val) {
                                val = G.tri_core[G.eid[j]];
                            }
                        }
                    }
                    break;
                }
                case TRIANGLE_CORES: { // sum
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + G.tri_core[G.eid[j]];
                        }
                    }
                    val = val / ((double)(vertices[u+1] - vertices[u]));
                    break;
                }
                case TRIANGLE_CORE_MIN: {
                    val = G.num_vertices();
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            if (G.tri_core[G.eid[j]] < val) {
                                val = G.tri_core[G.eid[j]];
                            }
                        }
                    }
                    break;
                }
                case TRIANGLE_CORE_DIST_TWO: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val += G.tri_core[G.eid[j]];
                            int w = edges[j];
                            for (long long k = vertices[w]; k < vertices[w + 1]; k++) {
                                if (ind[edges[k]]) {
                                    val += G.tri_core[G.eid[k]]; //decay weight
                                }
                            }
                        }
                    }
                    break;
                }
                case KCORE_DEG: {
                    val = G.degree[u] * G.kcore[u];
                    break;
                }
                case DEGREE_TRIANGLES: {
                    val = G.degree[u] * G.t[u];
                    break;
                }
                case KCORE_TRIANGLES: {
                    val = G.kcore[u] * G.t[u];
                    break;
                }
                case KCORE_DEG_TRI: {
                    val = G.degree[u] * G.kcore[u] * G.t[u];
                    break;
                }
                case DEGREE_KCORE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + ((vertices[edges[j] + 1] - vertices[edges[j]]) * G.kcore[edges[j]]);
                        }
                    }
                    break;
                }
                case KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + (G.kcore[edges[j]] * G.t[edges[j]]);
                        }
                    }
                    break;
                }
                case DEGREE_KCORE_TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vertices[u]; j < vertices[u + 1]; j++) {
                        if (ind[edges[j]]) {
                            val = val + (G.degree[edges[j]] * G.kcore[edges[j]] * G.t[edges[j]]);
                        }
                    }
                    break;
                }
                default: {
                    val = vertices[u + 1] - vertices[u];
                    break;
                }
            }
            V[i].set_bound(val);
//            V.push_back(Vertex(u,val));

        }
    }

    inline
    static void get_ordering(vector<Vertex> &V, bool smallest_to_largest = true) {
    if (smallest_to_largest) // smallest to largest
        std::sort(V.begin(), V.end(), incr_bound);
    else            // largest to smallest
        std::sort(V.begin(), V.end(), decr_bound);
}

    // if vertices are pruned implicitly, first induce the subgraph
    static void order_vertices(vector<Vertex> &V, graphpack_graph &G, int ordering, bool smallest_to_largest = true) {

        ensure_graph_parameters(G, ordering);   // check and compute graph parameters if needed
        if (V.size() != 0) {
            print_line(80);printf("ERROR: |V| > 0 : V is not empty, use local_order_vertices \n");print_line(80);
        }
        set_graph_parameters(V,G,ordering);     // add selected graph parameter to V
        get_ordering(V,smallest_to_largest);      // order vertices

    }


    /**
     * Fastest but least accurate
     *
     * @param V
     * @param G
     * @param vs
     * @param es
     * @param ordering
     * @param smallest_to_largest
     */
    inline
    static void order_local_vertices(
            vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> &vs,
            vector<int> &es,
            int ordering,
            bool smallest_to_largest = true) {

//        if (ordering == KCORE) {
//            neigh_cores_bound(vs,es,P,ind,max);
//        }
//        else if (ordering == TRIANGLES) {
//
//        }

//        ensure_graph_parameters(G, ordering);   // check and compute graph parameters if needed
        set_graph_parameters(V,G,vs,es,ordering);     // add selected graph parameter to V
        get_ordering(V,smallest_to_largest);      // order vertices

    }

    /**
     * Slower than order_local_vertices, but more accurate, since
     * ind array is used to ensure egonet/vol related properties use only the neighbors
     *
     * Also, this method is faster than exact, but less accurate
     *
     * @param V
     * @param G
     * @param vs
     * @param es
     * @param ind
     * @param ordering
     * @param smallest_to_largest
     */
    // uses ind to ensure egonet/vol related properties use only the neighbors
    inline
    static void order_local_vertices_induced(
            vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> &vs,
            vector<int> &es,
            vector<int> &ind,
            int ordering,
            bool smallest_to_largest = true) {

//        if (ordering == KCORE) {
//            neigh_cores_bound(vs,es,P,ind,max);
//        }
//        else if (ordering == TRIANGLES) {
//
//        }

//        ensure_graph_parameters(G, ordering);   // check and compute graph parameters if needed
        set_graph_parameters_induced(V,G,vs,es,ind,ordering);     // add selected graph parameter to V
        get_ordering(V,smallest_to_largest);      // order vertices

    }

    /**
     * Exact computation of local properties of the set of vertices P
     * Slowest but most accurate!
     *
     * @param V
     * @param G
     * @param vs
     * @param es
     * @param ind
     * @param ordering
     * @param smallest_to_largest
     */
    // uses ind/marked to ensure properties are computed exact, slowest, but most accurate to ensure egonet/vol related properties use only the neighbors
    inline
    static void order_local_vertices_exact(
            vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> &vs,
            vector<int> &es,
            vector<int> &ind,
            vector<int> &marked,
            int ordering,
            bool smallest_to_largest = true) {

//        if (ordering == KCORE) {
//            neigh_cores_bound(vs,es,P,ind,max);
//        }
//        else if (ordering == TRIANGLES) {
//
//        }

//        ensure_graph_parameters(G, ordering);   // check and compute graph parameters if needed
        set_graph_parameters_exact(V,G,vs,es,ind,marked,ordering);     // add selected graph parameter to V
        get_ordering(V,smallest_to_largest);      // order vertices

    }


    inline
    static void get_parameter(graphpack_graph & G, int & ordering, int & u, int & val) {

        switch (ordering) {
            case NATURAL: { // natural order (read from disk)
                val = u;
                break;
            }
            case DEGREE: {
                val = G.vertices[u + 1] - G.vertices[u];
                break;
            }
            case WEDGES: {
                val = G.vertices[u + 1] - G.vertices[u]; //store degree in val
                val = (val*(val-1))/2;
                break;
            }
            case KCORE: {
                val = G.kcore[u];
                break;
            }
            case TRIANGLES: {
                val = G.t[u];
                break;
            }
            case RAND: {
//                    val = drand() % G.vertices.size();
                val = rand() % G.vertices.size();
                break;
            }
            case VAR: {
                val = G.kcore[u] * ((int)G.degree[u]/G.kcore[u]);
                break;
            }
            case DEGREE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]];//G.vertex_degree(edges[j]);
                }
                break;
            }
            case KCORE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + G.kcore[G.edges[j]];
                }
                break;
            }
            case TRIANGLE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + G.t[G.edges[j]];
                }
                break;
            }
            case TRIANGLE_CORE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + G.tri_core[G.eid[j]];
                }
                break;
            }
            case TRIANGLE_CORE_MAX: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    if (G.tri_core[G.eid[j]] > val) {
                        val = G.tri_core[G.eid[j]];
                    }
                }
                break;
            }
            case TRIANGLE_CORES: { // sum
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + G.tri_core[G.eid[j]];
                }
                val = val / ((double)(G.vertices[u+1] - G.vertices[u]));
                break;
            }
            case TRIANGLE_CORE_MIN: {
                val = G.num_vertices();
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    if (G.tri_core[G.eid[j]] < val) {
                        val = G.tri_core[G.eid[j]];
                    }
                }
                break;
            }
            case TRIANGLE_CORE_DIST_TWO: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val += G.tri_core[G.eid[j]];
                    int w = G.edges[j];
                    for (long long k = G.vertices[w]; k < G.vertices[w + 1]; k++) {
                            val += G.tri_core[G.eid[k]]; //decay weight
                    }
                }
                break;
            }
            case KCORE_DEG: {
                val = G.degree[u] * G.kcore[u];
                break;
            }
            case DEGREE_TRIANGLES: {
                val = G.degree[u] * G.t[u];
                break;
            }
            case KCORE_TRIANGLES: {
                val = G.kcore[u] * G.t[u];
                break;
            }
            case KCORE_DEG_TRI: {
                val = G.degree[u] * G.kcore[u] * G.t[u];
                break;
            }
            case DEGREE_KCORE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + ((G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]]) * G.kcore[G.edges[j]]);
                }
                break;
            }
            case KCORE_TRIANGLE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + (G.kcore[G.edges[j]] * G.t[G.edges[j]]);
                }
                break;
            }
            case DEGREE_KCORE_TRIANGLE_VOL: {
                val = 0;
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                    val = val + (G.degree[G.edges[j]] * G.kcore[G.edges[j]] * G.t[G.edges[j]]);
                }
                break;
            }
            default: {
                val = G.vertices[u + 1] - G.vertices[u];
                break;
            }
        }
    }

    /**
     * Used for ordering vertices in G
     * V is assumed to be empty
     * Uses G directly!
     * Processes vertices in order of k-cores
     */
    static void compute_ordering(vector<Vertex> &V, graphpack_graph &G,
            int ordering, bool decr_order = true, int lb_idx = 0, int lb = 0) {

        if (decr_order)
            cout << "[computing ordering]  " << ordering << ", direction = largest to smallest" <<endl;
        else
            cout << "[computing ordering]  " << ordering << ", direction = smallest to largest" <<endl;


        srand (time(NULL));
        int u = 0, val = 0;
        for (int k = lb_idx; k < G.num_vertices(); k++) {
            if (G.degree[G.kcore_order[k]] >= lb - 1) {
                u = G.kcore_order[k];

                get_parameter(G,ordering,u,val); // val is returned by reference

                V.push_back(Vertex(u,val));
            }
        }
        if (decr_order) // smallest to largest
            std::sort(V.begin(), V.end(), incr_bound);
        else            // largest to smallest
            std::sort(V.begin(), V.end(), decr_bound);
    }


    static int optimize_block_size(graphpack_graph & G, int & num_workers, int n = 0) {
        if (n == 0)
            n = G.num_vertices();

        if (n < 10000) {
            return 1;
        }
        else { //if (n < 100000)
            return 64;
        }
    }

    /**
     * + Orders only vertices in V
     * + Assumes that the ordering strategy (triangles, kcores)
     *    has been computed and already exists in G
     * + Uses G directly, so this may be old
     */
    inline
    static void compute_local_ordering(vector<Vertex> &V, graphpack_graph &G,
            int ordering, bool decr_order = true, int lb_idx = 0, int lb = 0) {

        if (ordering == KCORE && G.kcore.size()) {
            G.compute_cores();
            printf("k-core numbers not found, computing them now... max k-core number = %d \n ", G.max_core);
            G.basic_stats();
        }

        if (ordering == TRIANGLES && G.t.size()) {
            G.compute_cores();
            printf("triangles not found, computing them now... \n");
            G.basic_stats();
        }

        srand (time(NULL));
        int u = 0, val = 0;
        for (int k = 0; k < V.size(); ++k) {
            u = V[k].get_id();

            get_parameter(G,ordering,u,val); // val is returned by reference



//
//            switch (ordering) {
//                case NATURAL: { // natural order (read from disk)
//                    val = u;
//                    break;
//                }
//                case DEGREE: {
//                    val = G.vertices[u + 1] - G.vertices[u];
//                    break;
//                }
//                case WEDGES: {
//                    val = G.vertices[u + 1] - G.vertices[u]; //store degree in val
//                    val = (val*(val-1))/2;
//                    break;
//                }
//                case KCORE: {
//                    val = G.kcore[u];
//                    break;
//                }
//                case TRIANGLES: {
//                    val = G.t[u];
//                    break;
//                }
//                case RAND: {
////                    val = drand() % G.vertices.size();
//                    val = rand() % G.vertices.size();
//                    break;
//                }
//                case VAR: {
//                    val = G.kcore[u] * ((int)G.degree[u]/G.kcore[u]);
//                    break;
//                }
//                case DEGREE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]];//G.vertex_degree(edges[j]);
//                    }
//                    break;
//                }
//                case KCORE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + G.kcore[G.edges[j]];
//                    }
//                    break;
//                }
//                case TRIANGLE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + G.t[G.edges[j]];
//                    }
//                    break;
//                }
//                case TRIANGLE_CORE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + G.tri_core[G.eid[j]];
//                    }
//                    break;
//                }
//                case TRIANGLE_CORE_MAX: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        if (G.tri_core[G.eid[j]] > val) {
//                            val = G.tri_core[G.eid[j]];
//                        }
//                    }
//                    break;
//                }
//                case TRIANGLE_CORES: { // sum
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + G.tri_core[G.eid[j]];
//                    }
//                    val = val / ((double)(G.vertices[u+1] - G.vertices[u]));
//                    break;
//                }
//                case TRIANGLE_CORE_MIN: {
//                    val = G.num_vertices();
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        if (G.tri_core[G.eid[j]] < val) {
//                            val = G.tri_core[G.eid[j]];
//                        }
//                    }
//                    break;
//                }
//                case TRIANGLE_CORE_DIST_TWO: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val += G.tri_core[G.eid[j]];
//                        int w = G.edges[j];
//                        for (long long k = G.vertices[w]; k < G.vertices[w + 1]; k++) {
//                                val += G.tri_core[G.eid[k]]; //decay weight
//                        }
//                    }
//                    break;
//                }
//                case KCORE_DEG: {
//                    val = G.degree[u] * G.kcore[u];
//                    break;
//                }
//                case DEGREE_TRIANGLES: {
//                    val = G.degree[u] * G.t[u];
//                    break;
//                }
//                case KCORE_TRIANGLES: {
//                    val = G.kcore[u] * G.t[u];
//                    break;
//                }
//                case KCORE_DEG_TRI: {
//                    val = G.degree[u] * G.kcore[u] * G.t[u];
//                    break;
//                }
//                case DEGREE_KCORE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + ((G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]]) * G.kcore[G.edges[j]]);
//                    }
//                    break;
//                }
//                case KCORE_TRIANGLE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + (G.kcore[G.edges[j]] * G.t[G.edges[j]]);
//                    }
//                    break;
//                }
//                case DEGREE_KCORE_TRIANGLE_VOL: {
//                    val = 0;
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        val = val + (G.degree[G.edges[j]] * G.kcore[G.edges[j]] * G.t[G.edges[j]]);
//                    }
//                    break;
//                }
//                default: {
//                    val = G.vertices[u + 1] - G.vertices[u];
//                    break;
//                }
//            }
//
            V[k].set_bound(val);
        }
        if (decr_order) // smallest to largest
            std::sort(V.begin(), V.end(), incr_bound);
        else            // largest to smallest
            std::sort(V.begin(), V.end(), decr_bound);

    }


    /**
     * + Orders only vertices in V
     * + Assumes that the ordering strategy (triangles, kcores)
     *    has been computed and already exists in G
     * + Uses G directly, to get bound!
     * + Uses V and E for volume and degree measures!
     *      For parallel, ensures these are always up to date,
     *      and thus the tightest possible.
     *
     * For explicit neighborhood induction branch methods, the previous one works fine where G is used
     *      Since G is the updated neighborhood
     * However, for CSC and Hybrid representations then the updated es and vs should be used.
     */
    inline
    static void compute_local_ordering(vector<Vertex> &V,
            graphpack_graph &G,
            vector<long long> &vs,
            vector<int> &es,
            int ordering,
            bool decr_order = true,
            int lb_idx = 0,
            int lb = 0) {

        if (ordering == KCORE && G.kcore.size()) {
            G.compute_cores();
            printf("k-core numbers NOT found, computing them now... max k-core number = %d \n ", G.max_core);
            G.basic_stats();
        }

        if (ordering == TRIANGLES && G.t.size()) {
            G.compute_cores();
            printf("triangles NOT found, computing them now... \n");
            G.basic_stats();
        }


        srand (time(NULL));
        int u = 0, val = 0;
        for (int k = 0; k < V.size(); k++) {
            u = V[k].get_id();

            switch (ordering) {
                case DEGREE: {
                    val = vs[u + 1] - vs[u];
                    break;
                }
                case KCORE: {
                    val = G.kcore[u];
                    break;
                }
                case TRIANGLES: {
                    val = G.t[u];
                    break;
                }
                case KCORE_DEG: {
                    val = (vs[u + 1] - vs[u]) * G.kcore[u];
                    break;
                }
                case RAND: {
                    val = rand() % vs.size()-1;
                    break;
                }
                case VAR: {
                    val = G.kcore[u] * ((int)(vs[u + 1] - vs[u])/G.kcore[u]);
                    break;
                }
                case DEGREE_VOL: {
                    val = 0;
                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                        val = val + vs[es[j] + 1] - vs[es[j]];//G.vertex_degree(edges[j]);
                    }
                    break;
                }
                case KCORE_VOL: {
                    val = 0;
                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                        val = val + G.kcore[es[j]];
                    }
                    break;
                }
                case TRIANGLE_VOL: {
                    val = 0;
                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                        val = val + G.t[es[j]];
                    }
                    break;
                }
                default: {
                    val = vs[u + 1] - vs[u];
                    break;
                }
            }
            V[k].set_bound(val);
        }
        if (decr_order) // smallest to largest
            std::sort(V.begin(), V.end(), incr_bound);
        else            // largest to smallest
            std::sort(V.begin(), V.end(), decr_bound);

    }




// dynamic search and local
    inline
    static void compute_local_ordering(
            vector<Vertex> &V,
//            graphpack_graph &G,
            vector<long long> &vs,
            vector<int> &es,
            int ordering,
            bool decr_order = true,
            int lb_idx = 0,
            int lb = 0) {




        // skip, since k-core and triangles assumed to be stored in the array of vertices
//        if (ordering != KCORE && ordering != TRIANGLES) {
            srand (time(NULL));
            int u = 0, val = 0;
            for (int k = 0; k < V.size(); k++) {
                u = V[k].get_id();

                switch (ordering) {
                    case DEGREE: {
                        val = vs[u + 1] - vs[u];
                        break;
                    }
                    case KCORE: {

                        val = V[k].get_bound();
//                        printf("local k-core K(%d) = %d \n", u, val);
                        break;
                    }
                    case TRIANGLES: {
                        //                    val = G.t[u];
                        val = V[k].get_bound();
//                        printf("local triangle count tr(%d) = %d \n", u, val);
                        break;
                    }
                    case KCORE_DEG: {
                        val = (vs[u + 1] - vs[u]) * V[k].get_bound();
                        break;
                    }
                    case RAND: {
                        val = rand() % vs.size()-1;
                        break;
                    }
                    case DEGREE_VOL: {
                        val = 0;
                        for (long long j = vs[u]; j < vs[u + 1]; j++) {
                            val = val + vs[es[j] + 1] - vs[es[j]];//G.vertex_degree(edges[j]);
                        }
                        break;
                    }
                    //                case KCORE_VOL: {
                    //                    val = 0;
                    //                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                    //                        val = val + G.kcore[es[j]];
                    //                    }
                    //                    break;
                    //                }
                    //                case TRIANGLE_VOL: {
                    //                    val = 0;
                    //                    for (long long j = vs[u]; j < vs[u + 1]; j++) {
                    //                        val = val + G.t[es[j]];
                    //                    }
                    //                    break;
                    //                }
                    default: {
                        val = vs[u + 1] - vs[u];
                        break;
                    }
                }
                V[k].set_bound(val);
            } // end for each vertex
//        }
//        else {
//            printf("LOCAL ORDERING: TRIANGLES OR KCORES DETECTED! \n");
//        }

        if (decr_order) // smallest to largest
            std::sort(V.begin(), V.end(), incr_bound);
        else            // largest to smallest
            std::sort(V.begin(), V.end(), decr_bound);

    }





//    graphpack_graph H = G.preprocess_triangle_core(G,p.lb);


    // uses pruned array, might also just use a vector<Vertex> array
    static int pruning(graphpack_graph& G, int max, string pruning_method, int* &pruned, bool is_explicit_induced = true, bool is_enum = false) {

        int lb_idx = 0;

        if (pruning_method == "kcore" || pruning_method == "kcores" || pruning_method == "cores") {
            if (G.kcore.size() == 0) {
                G.compute_cores();
            }

            if (is_enum) {
                if (G.is_dense_graph()) {
                    for (int i = G.num_vertices()-1; i >= 0; i--) {
                        if (G.kcore[G.kcore_order[i]] == max)  lb_idx = i;
                        if (G.kcore[G.kcore_order[i]] < max) {
                            pruned[G.kcore_order[i]] = 1;
                            for (long long j = G.vertices[G.kcore_order[i]]; j < G.vertices[G.kcore_order[i] + 1]; j++) {
                                G.adj[G.kcore_order[i]][G.edges[j]] = false;
                                G.adj[G.edges[j]][G.kcore_order[i]] = false;
                            }
                        }
                    }
                }
                else {
                    for (int i = G.num_vertices()-1; i >= 0; i--) {
                        if (G.kcore[G.kcore_order[i]] == max)  lb_idx = i;
                        if (G.kcore[G.kcore_order[i]] < max) pruned[G.kcore_order[i]] = 1;
                    }
                }
            }
            else {
                if (G.is_dense_graph()) {
                    for (int i = G.num_vertices()-1; i >= 0; i--) {
                        if (G.kcore[G.kcore_order[i]] == max)  lb_idx = i;
                        if (G.kcore[G.kcore_order[i]] <= max) {
                            pruned[G.kcore_order[i]] = 1;
                            for (long long j = G.vertices[G.kcore_order[i]]; j < G.vertices[G.kcore_order[i] + 1]; j++) {
                                G.adj[G.kcore_order[i]][G.edges[j]] = false;
                                G.adj[G.edges[j]][G.kcore_order[i]] = false;
                            }
                        }
                    }
                } // if adj struct
                else {
                    for (int i = G.num_vertices()-1; i >= 0; i--) {
                        if (G.kcore[G.kcore_order[i]] == max)  lb_idx = i;
                        if (G.kcore[G.kcore_order[i]] <= max)  pruned[G.kcore_order[i]] = 1;
                    }
                }
            }
        }// kcore

        if (is_explicit_induced) {
            cout << "[maxpack: initial k-core pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
            G.reduce_graph(pruned);
            G.update_degrees();
            cout << "[graphpack: initial k-core pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;
        }
//        reduce_and_sort_edges(G, pruned, lb_idx);
        return lb_idx;
    }


    static graphpack_graph pruning_vertex_set(graphpack_graph & G, int & max, string pruning_method, bool is_parallel = true) {

        cout << "pruning_vertex_set" <<endl;

        if (pruning_method == "kcore" || pruning_method == "kcores") {

            // this code relaxes the pruning using the lb
            vector<int> W;
            W.reserve(G.num_vertices());
            for (int i = G.num_vertices()-1; i >= 0; i--) {
                if (G.kcore[G.kcore_order[i]] >= max)
                    W.push_back(G.kcore_order[i]);
            }

            double sec = get_time();
            graphpack_graph H(W,G.vertices,G.edges);
            return H;

        }
        else if (pruning_method == "tcores" || pruning_method == "triangle_core") {

            // this code relaxes the pruning using the lb
            vector<int> W;
            W.reserve(G.num_vertices());
            for (int i = G.num_vertices()-1; i >= 0; i--) {
                if (G.kcore[G.kcore_order[i]] >= max)
                    W.push_back(G.kcore_order[i]);
            }

            double sec = get_time();
            graphpack_graph H(W,G.vertices,G.edges);


            /*
             * Compute Triangle Cores from H; the induced k-core graph
             */


            // given position in edge array (neighbor), eid returns the unique edge id
            H.map_eid_to_vid(is_parallel);

        //    double triangle_core_time;
        //    double fast_triangle_core_time;
        //    double edge_triangles;
        //    double fast_edge_triangles;

//            double sec = get_time();
            // number of triangles each edge participates
            vector<int> pruned_edges(H.num_edges() + 1,0);
            if (H.is_dense_graph()) // dense and/or small adj
                H.triangle_counting_adj(max, pruned_edges, is_parallel);
            else // sparse
                H.triangle_counting(max, pruned_edges, is_parallel);

            cout << "triangles: " << H.total_t << " (time: " << get_time() - sec << "), " <<endl;


            // compute the triangle cores
            sec = get_time();
            H.triangle_core_numbers();
            cout << "max triangle core: " << H.max_tri_core+2 << " (time: " << get_time() - sec << "), " <<endl;










//            int threshold = 15000;
////            double fast_triangle_core_perf = 0;
////            double fast_triangles_perf = 0;
////            double fast_triangle_cores = 0;
//            if (H.vertices.size() < threshold && H.adj != NULL) {
////                cout << "verts in adj: " << H.vertices.size() <<endl;
////                sec = get_time();
//                H.triangle_cores_approx_adj(fast_triangle_core_perf,fast_triangles_perf);
////                fast_triangle_cores = get_time() - sec;
////                cout << "[mcpack]  time to compute triangle cores from H: " << get_time() - sec << endl;
//            }
//            else {
////                sec = get_time();
//                H.triangle_cores_approx(fast_triangle_core_perf,fast_triangles_perf);
////                fast_triangle_cores = get_time() - sec;
//            }


            vector<long long> VS(H.vertices.size(),0);
            vector<int> E;
            E.reserve(H.edges.size());


            vector<Vertex> V;
            V.reserve(H.num_vertices());

            long long start = 0;
            long long edge_pos = 0;

            /*
             * Removes edges with Triangle Core less than K
             * After removing the edges with triangle core less than K,
             * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
             */
            int num_edges_pruned = 0;
            int num_verts_pruned = 0;
            if (max == G.max_tri_core+2) {
                cout << "inside tricore pruning, lb == max tri core..." <<endl;
                for (int i = 0; i < H.num_vertices(); i++) {
                    start = E.size();
                    for (long long j = H.vertices[i]; j < H.vertices[i + 1]; j++ ) {

                        edge_pos = H.eid[j]; // tri_core is m/2, so this finds the right tri_core value given an edge id
                        if (H.tri_core[edge_pos]+2 >= max) { // get correct tri_core for edge's eid
                            E.push_back(H.edges[j]);
                        }
                        else {
                            num_edges_pruned++;
                        }
                    }
                    if ((E.size() - start) >= max-1) { // basically the updated degree bound
                        V.push_back(Vertex(i,0));
                    } else num_verts_pruned++;
                    VS[i] = start;
                    VS[i + 1] = E.size();
                }
            }
            else {
                for (int i = 0; i < H.num_vertices(); i++) {
                    start = E.size();
                    for (long long j = H.vertices[i]; j < H.vertices[i + 1]; j++ ) {

                        edge_pos = H.eid[j]; // tri_core is m/2, so this finds the right tri_core value given an edge id
                        if (H.tri_core[edge_pos]+2 >= max) { // get correct tri_core for edge's eid
                            E.push_back(H.edges[j]);
                        }
                        else {
                            num_edges_pruned++;
                        }
                    }
                    if ((E.size() - start) >= max-1) { // basically the updated degree bound, but what about the edges???
                        V.push_back(Vertex(i,0));
                    }
                    else num_verts_pruned++;
                    VS[i] = start;
                    VS[i + 1] = E.size();
                }
            }



            //    cout << "after triangle core pruning..." <<endl;
            //
            //    cout << "num_verts pruned = " << num_verts_pruned << " out of " << G.num_vertices() <<endl;
            //    cout << "num edges pruned = " << num_edges_pruned << " out of " << G.edges.size() <<endl;
            //
            //    cout << "\n**********************************" <<endl;
            //    cout << "T: " << G.max_tri_core+2 <<endl;
            //    cout << "T(H): " << H.max_tri_core+2 <<endl;
            //    cout << "T_fast_time:\t " << fast_triangle_cores <<endl;
            //
            //
            //    G.basic_stats("G-");
            //    cout << "K-|T|: " << H.total_t <<endl;
            //    cout << "K-|T|_time: " << fast_triangles_perf <<endl;
            graphpack_graph F(V,VS,E);
            F.basic_stats("T-");
            return F;
        }
//        else { // degree
//            if (pruning_method == "kcore" || pruning_method == "kcores") {
//
//                // this code relaxes the pruning using the lb
//                vector<int> W;
//                W.reserve(G.num_vertices());
//                for (int i = G.num_vertices()-1; i >= 0; i--) {
//                    if (G.kcore[G.kcore_order[i]] >= lb)
//                        W.push_back(G.kcore_order[i]);
//                }
//
//                double sec = get_time();
//                graphpack_graph H(W,G.vertices,G.edges);
//                return H;
//
//            }
//        }
        return G;

    }



    static void explicit_reduce(graphpack_graph& G, int* &pruned, int lb_idx) {
        double sec = get_time();
        if (G.verbose) {
            cout << "[graphpack: initial k-core pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
        }

        G.reduce_graph(pruned);


        if (G.verbose) {
            cout << "[graphpack: initial k-core pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;
            cout << "[mcpack]  initial pruning took " << get_time()-sec << " sec" <<endl;
        }

//        G.update_degrees();
//    //    G.degree_bucket_sort(true);
////        G.degree_bucket_sort_parallel(true);
//        edge_bucket_sort_parallel(G, G.kcore, true);

//        edge_bucket_sort_parallel(graphpack_graph &G, vector<int> &bound, bool desc);
    }





    /*
     * Each worker sorts the edges of a vertex, and
     * afterwards, grabs the next available vertex to sort
     *
     * Note: sort neighbors by degree (largest to smallest)
     */
    static void edge_bucket_sort_parallel(graphpack_graph &G, vector<int> &bound, params &p) {

        cout << "ordering edges/neighbors of the vertices in parallel" <<endl;

    //  int v, u, n, md, md_end, start, d, num;

        vector<int> tmp_edges(G.edges.size(),0);
    //  tmp_edges.reserve(edges.size());
        int v = 0;
        double sec = get_time();
        #pragma omp parallel for schedule(dynamic,64) \
            shared(tmp_edges)
    //  firstprivate(colors,ind,vs,es) \
    //  private(u, P, C, C_max, mc_cur)
        for (v = 0; v < G.num_vertices(); ++v) {

            int u, n, md, md_end, start, d, num;

            n = G.vertices[v+1] - G.vertices[v] + 1;
            vector<int> vert(n);
            vector<int> pos(n);
            vector<int> deg(n);

            md = 0;
            for(u=1; u<n; u++) {
                deg[u] = bound[G.edges[G.vertices[v] + (u-1)]];
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
                vert[pos[u]] = G.edges[G.vertices[v] + (u-1)];
                bin[deg[u]]++;
            }

    //      int insert_pos = vertices[v];
            if (p.global_edge_small_to_large) {
                //from smallest degree to largest
                for (u=n-1; u>0; --u) {
                    tmp_edges[G.vertices[v] + (u-1)] = vert[u];
    //                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
                }
            }
            else {
                // largest to smallest, example: 5,3,2
                for (u = 1; u < n; ++u) {
    //                tmp_edges[vertices[v] + (u-1)] = vert[u];
                    tmp_edges[G.vertices[v] + (u-1)] = vert[n-u];
                }
            }

//    //            for (u=n-1; u>0; --u) {
//    ////                tmp_edges[vertices[v] + (u-1)] = vert[u];
//    //                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
//    //            }
//
//    //            print_break();
//    //            cout << "v = " << v <<endl;
//    //            for (int ii=0; ii < num_vertices(); ii++)
//    //                for (u=vertices[ii]; u < vertices[ii+1]; u++) {
//    //                    cout << "\t(" << ii << "," << tmp_edges[u];
//    //                    cout << "), d = " << degree[tmp_edges[u]] <<endl;
//    //                }
//
//                // largest to smallest
//    //          tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.rbegin(),vert.rend()-1);
//            }
//            else {
//                //from smallest degree to largest
//                for (u=n-1; u>0; --u) {
//                    tmp_edges[G.vertices[v] + (u-1)] = vert[u];
//    //                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
//                }
//
//    //            for (u=1; u<n; u++) {
//    //                tmp_edges[vertices[v] + (u-1)] = vert[u];
//    //            }
//
//                //from smallest degree to largest
//    //          tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.begin()+1,vert.end());
//            }
        }

        if (G.verbose) {
            cout << "[graphpack: sorting neighbors]  |E| = " << G.edges.size();
            cout << ", |E_sorted| = " << tmp_edges.size() <<endl;
            cout << "sorting neighbors took " << get_time() - sec <<endl;
        }
        G.edges = tmp_edges;
    }



    static void triangle_core_remove_edges(graphpack_graph &G, int* &pruned, int & max) {


        vector<long long> V(G.vertices.size(),0);
        vector<int> E;
        E.reserve(G.edges.size());

        long long start = 0;
        long long edge_pos = 0;


        /*
         * Removes edges with Triangle Core less than K
         * After removing the edges with triangle core less than K,
         * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
         */

        if (G.is_dense_graph()) {

            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;

                    if (G.tri_core[edge_pos+1]+2 > max) { // get correct tri_core for edge's eid
                        E.push_back(G.edges[j]);
                    }
                    else {
                        G.adj[G.e_u[edge_pos]][G.e_v[edge_pos]] = false;
                        G.adj[G.e_v[edge_pos]][G.e_u[edge_pos]] = false;
                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                }
            }
        }
        else {
            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;

                    if (G.tri_core[edge_pos+1]+2 > max) { // get correct tri_core for edge's eid
                        E.push_back(G.edges[j]);
                    }
                    //                    else {
                    //                        G.adj[G.e_u[edge_pos]][G.e_v[edge_pos]] = false;
                    //                        G.adj[G.e_v[edge_pos]][G.e_u[edge_pos]] = false;
                    //                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                }
            }
        }
    }



    static void triangle_core_remove_edges(graphpack_graph &G, vector<int> &pruned, vector<int> &pruned_edges, int & max) {


        vector<long long> V(G.vertices.size(),0);
        vector<int> E;
        E.reserve(G.edges.size());

        long long start = 0;
        long long edge_pos = 0;


        /*
         * Removes edges with Triangle Core less than K
         * After removing the edges with triangle core less than K,
         * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
         */

        if (G.is_dense_graph()) {

            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;
                    if (!pruned_edges[edge_pos]) {
                        if (G.tri_core[edge_pos]+2 > max) { // get correct tri_core for edge's eid
                            E.push_back(G.edges[j]);
                        }
                        else {
                            G.adj[G.e_u[edge_pos]][G.e_v[edge_pos]] = false;
                            G.adj[G.e_v[edge_pos]][G.e_u[edge_pos]] = false;
                        }
                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                    G.degree[i] = 0;
                }
            }
        }
        else {
            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;
                    if (!pruned_edges[edge_pos]) {
                        if (G.tri_core[edge_pos]+2 > max) { // get correct tri_core for edge's eid
                            E.push_back(G.edges[j]);
                        }
                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                }
            }
        }
    }




    static void triangle_core_remove_edges(graphpack_graph &G, int* &pruned, vector<int> &pruned_edges, int & max) {


        vector<long long> V(G.vertices.size(),0);
        vector<int> E;
        E.reserve(G.edges.size());

        long long start = 0;
        long long edge_pos = 0;

        int num_pruned = 0;

        /*
         * Removes edges with Triangle Core less than K
         * After removing the edges with triangle core less than K,
         * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
         */

        if (G.is_dense_graph()) {

            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;
                    if (!pruned_edges[edge_pos]) {
                        if (G.tri_core[edge_pos]+2 > max) { // get correct tri_core for edge's eid
                            E.push_back(G.edges[j]);
                        }
                        else {
                            G.adj[G.e_u[edge_pos]][G.e_v[edge_pos]] = false;
                            G.adj[G.e_v[edge_pos]][G.e_u[edge_pos]] = false;
                        }
                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                    G.degree[i] = 0;
                    num_pruned++;
                }
            }
        }
        else {
            for (int i = 0; i < G.num_vertices(); i++) {
                start = E.size();
                for (long long j = G.vertices[i]; j < G.vertices[i + 1]; j++ ) {

                    // add 1 to edge_pos for tri_core, but not e_u and e_v
                    edge_pos = G.eid[j];// + 1;
                    if (!pruned_edges[edge_pos]) {
                        if (G.tri_core[edge_pos]+2 > max) { // get correct tri_core for edge's eid
                            E.push_back(G.edges[j]);
                        }
                    }
                }

                V[i] = start;
                V[i + 1] = E.size();
                G.degree[i] = V[i+1] - V[i];

                if (G.degree[i] < max) {
                    pruned[i] = 1;
                    num_pruned++;
                }
            }
        }
        cout << "[graphpack: triangle-core pruning]  |V_prev| = " << G.vertices.size()-1 <<
                ", |V| = " << (G.vertices.size()-1)-num_pruned << ", " <<
                "pruned vertices = " << num_pruned <<endl;

        cout << "[graphpack: triangle-core pruning]  |E_prev| = " << G.num_edges()<<
                ", |E| = " << (E.size()/2) << ", " <<
                "pruned edges = " << (G.num_edges()) - (E.size()/2) <<endl;

        G.vertices = V;
        G.edges = E;
    }


    /**
     * num tris saved in P
     * Parallel Vertex Triangle Counting
     *
     * Handles both adj and csc representations
     *
     * @param max
     * @param is_parallel
     */
    static void compute_vertex_triangles_dense(vector<Vertex> &P, int* &pruned, int max, bool is_parallel,
            vector<long long>& vs, vector<int>& es) {
        long long v, u, w, i, j, k;
        long long n = vs.size()-1;
        //    long long n = P.size();
        double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;

        vector<long long> t;
        t.resize(n);
        //    kappa.resize(n);

        //    total_t = 0;
        vector<int> ind(n, 0);
        vector<int> tris(n, 0);
        //    vector<int> tris_max(n, 0);
        //    vector<int> tris_pruned(n, 0);

        //    if (is_dense_graph()) {
//        if (is_parallel) {
//#pragma omp parallel for schedule(dynamic,64) \
//        shared(tris) firstprivate(ind) \
//private(v,u,w,d,i,j)
//            for (k = 0; k < P.size(); k++) {
//
//                v = P[k].get_id();
//
//                if (pruned[v]) continue;
//
//                d = vs[v+1] - vs[v];
//
//                //                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[es[i]] = 1;
//                //                ind[v] = 1;
//
//                for (i=vs[v]; i<vs[v+1]; i++) {
//                    w = es[i];
//                    if (v==w) { d-=1.; continue; }
//                    for (j=vs[w]; j<vs[w+1]; j++) {
//                        u = es[j];
//                        if (u==w) { continue; }
//                        if (u==v) { continue; }
//                        if (adj[u][v])  tris[v] += 1.;
//                    }
//                }
//                t[v] = tris[v]/2.;
//                P[k].set_bound(t[v]);
//
//                if (t[v]+2 < max) {
//                    t[v] = 0;
//                    P[k].set_bound(0);
//                    //                    pruned[v] = 1;
//                }
//            }
//        }
//        else {
//            for (k = 0; k < P.size(); k++) {
//
//                v = P[k].get_id();
//
//                if (pruned[v]) continue;
//
//                for (i=vs[v]; i<vs[v+1]; i++) {
//                    w = es[i];
//                    if (v==w) { d-=1.; continue; }
//                    for (j=vs[w]; j<vs[w+1]; j++) {
//                        u = es[j];
//                        if (u==w) { continue; }
//                        if (u==v) { continue; }
//                        if (adj[u][v])  tris[v] += 1.;
//                    }
//                }
//
//                t[v] = tris[v]/2.;
//                P[k].set_bound(t[v]);
//
//                if (t[v]+2 < max) {
//                    t[v] = 0;
//                    P[k].set_bound(0);
//                    //                    pruned[v] = 1;
//                }
//            }
//        }


        //    }
        //    else { // csc only
        if (is_parallel) {
#pragma omp parallel for schedule(dynamic,64) \
        shared(tris) firstprivate(ind) \
private(v,u,w,d,i,j)
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[v]) continue;

                d = vs[v+1] - vs[v];

                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 1;
                ind[v] = 1;

                for (i=vs[v]; i<vs[v+1]; i++) {
                    w = es[i];
                    if (v==w) { d-=1.; continue; }
                    for (j=vs[w]; j<vs[w+1]; j++) {
                        u = es[j];
                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u])  tris[v] += 1.;
                    }
                }
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
                }

                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 0;
                ind[v] = 0;
            }
        }
        else {
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[v]) continue;

                d = vs[v+1] - vs[v];

                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 1;
                ind[v] = 1;

                for (i=vs[v]; i<vs[v+1]; i++) {
                    w = es[i];
                    if (v==w) { d-=1.; continue; }
                    for (j=vs[w]; j<vs[w+1]; j++) {
                        u = es[j];
                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u])  tris[v] += 1.;
                    }
                }
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
                    //                    pruned[v] = 1;
                }

                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 0;
                ind[v] = 0;
            }
        }


        vector<Vertex> R;
        R.reserve(P.size());
        for (int i = 0; i < P.size(); ++i) {
            if (P[i].get_bound()+2 >= max) {
                R.push_back(P[i]);
            }
        }
        P = R;



        //    for (v = 0; v < n; v++) {
        //        total_t += t[v];
        //        t_sum += tris[v];
        //        t_pot += tris_max[v];
        //    }
        //
        //    global_cc = t_sum / t_pot;
        ind.clear();
        tris.clear();
        //    tris_max.clear();

        //    triangle_stats();
        //  triangle_bound();
    }




    // for GRAPH COLORING, upper bound is relaxed!
    static void prune_graph(graphpack_graph & G, params & p, int & max) {

        int* pruned = new int[G.num_vertices()];
        memset(pruned, 0, G.num_vertices() * sizeof(int));

        int global_pruning = set_enum(p.global_pruning);

        bool is_enum = false;
        if (p.problem == "mce") is_enum = true;

        int lb_idx = 0;
        switch (global_pruning) {
            case KCORE: {
                printf("pruning method: k-core \n");

                vector<int> W;
                for (int i = G.num_vertices()-1; i >= 0; i--) {
                    if (G.kcore[G.kcore_order[i]]+1 >= max)
                        W.push_back(G.kcore_order[i]);
                }
                graphpack_graph H(W,G.vertices,G.edges);
                H.basic_stats();
//                return H;
                G = H;

//                lb_idx = pruning(G, max, p.global_pruning, pruned, false, is_enum);
                break;
            }
            case TRIANGLES: {
                // selects between parallel and sequential or csc only vs. hybrid of csc/adj
                printf("pruning method: k-core and triangles \n");

                // first prune using k-cores
//                lb_idx = pruning(G, max, "kcore", pruned, false, is_enum);

                vector<int> W;
                for (int i = G.num_vertices()-1; i >= 0; i--) {
                    if (G.kcore[G.kcore_order[i]] >= max)
                        W.push_back(G.kcore_order[i]);
                }
                graphpack_graph H(W,G.vertices,G.edges);
                H.basic_stats();

                // then prune via triangles
                H.compute_vertex_triangles(pruned, max, true);
                G = H;
//                return H;
                break;
            }
//            case TRIANGLES_ONLY: {
//                printf("pruning method: triangles only \n");
//                // selects between parallel and sequential
//                // or csc only vs. hybrid of csc/adj
//                G.compute_vertex_triangles(pruned, max, true, is_enum);
//                break;
//            }
            case TRIANGLE_CORES: {
                // || global_pruning == "triangle_cores") {

                // first prune using k-cores
                //        lb_idx = pruning(G, max, "kcore", pruned, false, is_enum); // do NOT explicitly induced graph
                //        lb_idx = pruning(G, max, "kcore", pruned, true, is_enum); // explicitly induce the graph

                vector<int> W;
                for (int i = G.num_vertices()-1; i >= 0; i--) {
                    if (G.kcore[G.kcore_order[i]]+1 >= max)
                        W.push_back(G.kcore_order[i]);
                }
                graphpack_graph H(W,G.vertices,G.edges);
                H.basic_stats();

                H.map_edgeid_to_vertex_pair(true);

                vector<int> pruned_edges(G.num_edges()+1,0);
                H.k_triangle_counting(max,pruned_edges, true);
                printf("[triangle cores]  argmax tr(u) = %lld,  sum_u tr(u) = %llu \n", H.max_t_edge, H.total_t);
                //
                //        // compute the triangle cores
                H.triangle_core_numbers_parallel();

                /**
                 * this function explicitly removes edges,
                 * but only marks vertices as pruned
                 */
                //        triangle_core_remove_edges(G,pruned,max);
                triangle_core_remove_edges(H,pruned,pruned_edges,max);

                //        G.vertices = H.vertices;
                //        G.edges = H.edges;
                //        G.kcore = H.kcore;
                G = H;

//                return H;
                break;
            }
        } // end of global_pruning switch

        cout << "[graphpack: initial pruning]  before pruning: |V| = " << G.num_vertices() << ", |E| = " << G.num_edges() <<endl;
//        G.reduce_graph(pruned);
        G.update_degrees();
        cout << "[graphpack: initial pruning]  after pruning:  |V| = " << G.num_vertices() - lb_idx << ", |E| = " << G.num_edges() <<endl;

        G.compute_cores();
        printf("updated max core = %d \n", G.max_core);



        if (p.global_edge_ordering == "kcore") {
            printf("ordering edges by k-cores \n");
            edge_bucket_sort_parallel(G, G.kcore, p);
        }
        else {
            printf("ordering edges by degree (default) \n");
            edge_bucket_sort_parallel(G, G.degree, p);
        }

        if (pruned) delete[] pruned;

    }




    /**
     * Parallel Vertex Triangle Counting
     *
     * Handles both adj and csc representations
     *
     * @param max
     * @param is_parallel
     */
    static void compute_vertex_triangles_dense(vector<Vertex> &P, int* &pruned, int max, bool is_parallel,
            vector<long long>& vs, vector<int>& es, bool** &adj) {
        long long v, u, w, i, j, k;
        long long n = vs.size()-1;
        //    long long n = P.size();
        double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;

        vector<long long> t;
        t.resize(n);
        //    kappa.resize(n);

        //    total_t = 0;
        vector<int> ind(n, 0);
        vector<int> tris(n, 0);
        //    vector<int> tris_max(n, 0);
        //    vector<int> tris_pruned(n, 0);

        //    if (is_dense_graph()) {
        if (is_parallel) {
#pragma omp parallel for schedule(dynamic,64) \
        shared(tris) firstprivate(ind) \
private(v,u,w,d,i,j)
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[v]) continue;

                d = vs[v+1] - vs[v];

                //                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[es[i]] = 1;
                //                ind[v] = 1;

                for (i=vs[v]; i<vs[v+1]; i++) {
                    w = es[i];
                    if (v==w) { d-=1.; continue; }
                    for (j=vs[w]; j<vs[w+1]; j++) {
                        u = es[j];
                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
                    //                    pruned[v] = 1;
                }
            }
        }
        else {
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[v]) continue;

                for (i=vs[v]; i<vs[v+1]; i++) {
                    w = es[i];
                    if (v==w) { d-=1.; continue; }
                    for (j=vs[w]; j<vs[w+1]; j++) {
                        u = es[j];
                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }

                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
                    //                    pruned[v] = 1;
                }
            }
        }


        //    }
        //    else { // csc only
        //        if (is_parallel) {
        //            #pragma omp parallel for schedule(dynamic,64) \
        //                shared(tris) firstprivate(ind) \
        //                private(v,u,w,d,i,j)
        //            for (k = 0; k < P.size(); k++) {
        //
        //                v = P[k].get_id();
        //
        //                if (pruned[v) continue;
        //
        //                d = vs[v+1] - vs[v];
        //
        //                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 1;
        //                ind[v] = 1;
        //
        //                for (i=vs[v]; i<vs[v+1]; i++) {
        //                    w = es[i];
        //                    if (v==w) { d-=1.; continue; }
        //                    for (j=vs[w]; j<vs[w+1]; j++) {
        //                        u = es[j];
        //                        if (u==w) { continue; }
        //                        if (u==v) { continue; }
        //                        if (ind[u])  tris[v] += 1.;
        //                    }
        //                }
        //                t[v] = tris[v]/2.;
        //                P[k].set_bound(t[v]);
        //
        //                if (t[v]+2 < max) {
        //                    t[v] = 0;
        //                    P[k].set_bound(0);
        //                }
        //
        //                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 0;
        //                ind[v] = 0;
        //            }
        //        }
        //        else {
        //            for (k = 0; k < P.size(); k++) {
        //
        //                v = P[k].get_id();
        //
        //                if (pruned[v]) continue;
        //
        //                d = vs[v+1] - vs[v];
        //
        //                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 1;
        //                ind[v] = 1;
        //
        //                for (i=vs[v]; i<vs[v+1]; i++) {
        //                    w = es[i];
        //                    if (v==w) { d-=1.; continue; }
        //                    for (j=vs[w]; j<vs[w+1]; j++) {
        //                        u = es[j];
        //                        if (u==w) { continue; }
        //                        if (u==v) { continue; }
        //                        if (ind[u])  tris[v] += 1.;
        //                    }
        //                }
        //                t[v] = tris[v]/2.;
        //                P[k].set_bound(t[v]);
        //
        //                if (t[v]+2 < max) {
        //                    t[v] = 0;
        //                    P[k].set_bound(0);
        ////                    pruned[v] = 1;
        //                }
        //
        //                for (i=vs[v]; i<vs[v+1]; i++)  ind[es[i]] = 0;
        //                ind[v] = 0;
        //            }
        //        }
        //    }

        vector<Vertex> R;
        R.reserve(P.size());
        for (int i = 0; i < P.size(); ++i) {
            if (P[i].get_bound()+2 >= max) {
                R.push_back(P[i]);
            }
        }
        P = R;



        //    for (v = 0; v < n; v++) {
        //        total_t += t[v];
        //        t_sum += tris[v];
        //        t_pot += tris_max[v];
        //    }
        //
        //    global_cc = t_sum / t_pot;
        ind.clear();
        tris.clear();
        //    tris_max.clear();

        //    triangle_stats();
        //  triangle_bound();
    }


    inline
    static graphpack_graph global_kcore_pruning(graphpack_graph &G, int &max) {
        vector<int> W;
        for (int i = G.num_vertices()-1; i >= 0; i--) {
            if (G.kcore[G.kcore_order[i]] >= max)
                W.push_back(G.kcore_order[i]);
        }
        graphpack_graph H(W,G.vertices,G.edges);
        return H;
    }

    inline
    static graphpack_graph global_triangle_pruning(graphpack_graph &G, int &max) {
        vector<int> W;
        for (int i = 0; i < G.num_vertices(); ++i) {
            if (G.t[i]+2 >= max)
                W.push_back(i);
        }
        graphpack_graph H(W,G.vertices,G.edges);
        return H;
    }

    inline
    static graphpack_graph global_pruning_methods(graphpack_graph& G, int &max, int global_pruning) {
//        graphpack_graph H;
        switch (global_pruning) {
            case KCORE: {
                printf("pruning method: k-core \n");
                graphpack_graph H = global_kcore_pruning(G,max);
                G = H;
//                lb_idx = pruning(G, max, p.global_pruning, pruned, false, is_enum);
                break;
            }
            case TRIANGLES: {
                // selects between parallel and sequential or csc only vs. hybrid of csc/adj
                printf("pruning method: k-core and triangles \n");

                graphpack_graph H = global_kcore_pruning(G,max);
                H.compute_vertex_triangles(max,true);
                H = global_triangle_pruning(H,max);
                G = H;


//                // first prune using k-cores
//                lb_idx = pruning(G, max, "kcore", pruned, false, is_enum);
//
//                // then prune via triangles
//                G.compute_vertex_triangles(pruned, max, true);
                break;
            }
//            case TRIANGLES_ONLY: {
//                printf("pruning method: triangles only \n");
//                // selects between parallel and sequential
//                // or csc only vs. hybrid of csc/adj
//                G.compute_vertex_triangles(pruned, max, true, is_enum);
//                break;
//            }
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

                 H.triangle_core_numbers_parallel();

                 /**
                  * this function explicitly removes edges,
                  * but only marks vertices as pruned
                  */
                 //        triangle_core_remove_edges(G,pruned,max);
                 vector<int> pruned(G.num_vertices(),0);
                 triangle_core_remove_edges(H,pruned,pruned_edges,max);



                 H.compute_cores();
                 printf("updated max core = %d \n", H.max_core);
         //        G.vertices = H.vertices;
         //        G.edges = H.edges;
         //        G.kcore = H.kcore;
                 G = H;
                 break;
            }
        } // end of global_pruning switch
    }


}
#endif /* MCPACK_FRAMEWORK_H_ */

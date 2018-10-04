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

#ifndef GRAPHPACK_NEIGH_COLORING_H_
#define GRAPHPACK_NEIGH_COLORING_H_

#include "graphpack_vertex.h"

using namespace std;

namespace graphpack {

    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     *
     * @param v
     * @param k
     * @param colors
     * @param adj
     * @return
     */
    static bool repair_coloring(graphpack_graph &G, vector<Vertex> &P,
            vector< int >& color, vector<int> &used, int & v, int c) {

//        int k = c;
        vector<int> conflicts(G.max_degree+1,0);


//        if (G.verbose) printf("v=%d, color(v)=%d, c=%d \n", v,color[v],c);

//        vector<int> used_colors(G.max_degree+1,0);
        vector<int> forbidden_colors(G.max_degree+1,0);

        // mark the colors of neighbors, store their id in that color
        for (long long j = G.vertices[v]; j < G.vertices[v + 1]; ++j) {
            int w = G.edges[j]; // w is a neighbor of v



            if (color[w] > 0) {

                conflicts[color[w]]++;// = conflicts[color[w]] + 1;
                // get the color of w, index into used array with it, and mark that position using the label of v
//                used_colors[color[w]] = w; //mark color of w in the used array with the label of v.

                used[color[w]] = w; //mark color of w in the used array with the label of v.

//                if (G.verbose)
//                    printf("\t w=%d, color(w) = %d,  adjacent to %d vertices \n",
//                        w,color[w],conflicts[color[w]]);
            }

        }
        // find C_i such that v is adjacent to a single vertex w \in C_i
        for (int k = 1; k < c-1; ++k) { //must start from 1, since 0 is not a color

            if (conflicts[k] == 1) { // if we have a single conflicting color class

                // w MUST BE the single conflicting vertex, hence, only need a single vertex, so only need to mark it in used_colors!
//                int w = used_colors[k];
                int w = used[k];



//                if (G.verbose) printf("v=%d (c=%d) is adj to w=%d, k=%d, |N(w) intersect C_i| = %d \n",v,c, w, k, conflicts[k]);

                // mark the colors of w's neighbors
                for (long long j = G.vertices[w]; j < G.vertices[w + 1]; ++j) { // neighbors of w

                    int u = G.edges[j]; // neigh of w (conflicting with v)
                    // set position of u's color to vertex w, so that if
//                    used[color[u]] = w; // get color of u, and use it to index into used array, set id of w to that position



                    if (color[u] > 0) {
//                        if (G.verbose)  printf("\t\t u=%d, color=%d, mark w=%d, color=%d \n",u,color[u],w,color[w]);

                        // mark color of adjacent vertex u
                        forbidden_colors[color[u]] = w;
                    }
                }

                // find a color class C_j for which w does not conflict!
                for (int i = k+1; i < c-1; ++i) {
                    if (forbidden_colors[i] != w) { // check if color i marked with w, if so, then neighbor of u
                        if (k < color[v]) {
                            if (G.verbose) printf("repairing!!! v=%d assigned to color=%d,  w=%d assigned to color=%d \n ",v,color[w],w,c);
                            color[v] = color[w];
                            color[w] = i;
                            return true;
                        }
                    }
                } // for loop searching over adj colors
            } // if single conflict/adj vertex in color class
        } // search over color classes
        return false;
    }

    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     *
     * @param v
     * @param k
     * @param colors
     * @param adj
     * @return
     */
    static bool repair_coloring_parallel(graphpack_graph &G, vector<Vertex> &P,
            vector< int >& color, int v, int & c) {

        // shared variables: color array, and max colors


        // must only _share_ the color array
        vector<int> conflicts(G.max_degree+1,0);

        // both are private arrays
        vector<int> used_colors(G.max_degree+1,0);
        vector<int> forbidden_colors(G.max_degree+1,0);

        // mark the colors of neighbors, store their id in that color
        for (long long j = G.vertices[v]; j < G.vertices[v + 1]; ++j) {
            int w = G.edges[j]; // w is a neighbor of v

            if (color[w] > 0) {

                conflicts[color[w]]++;// = conflicts[color[w]] + 1;
                // get the color of w, index into used array with it, and mark that position using the label of v
                used_colors[color[w]] = w; //mark color of w in the used array with the label of v.

//                if (G.verbose)
//                    printf("\t w=%d, color(w) = %d,  adjacent to %d vertices \n",
//                        w,color[w],conflicts[color[w]]);
            }

        }
        // find C_i such that v is adjacent to a single vertex w \in C_i
        for (int k = 1; k < c-1; ++k) { //must start from 1, since 0 is not a color

            if (conflicts[k] == 1) { // if we have a single conflicting color class

                // w MUST BE the single conflicting vertex, hence, only need a single vertex, so only need to mark it in used_colors!
                int w = used_colors[k];

//                if (G.verbose) printf("v=%d (c=%d) is adj to w=%d, k=%d, |N(w) intersect C_i| = %d \n",v,c, w, k, conflicts[k]);

                // mark the colors of w's neighbors
                for (long long j = G.vertices[w]; j < G.vertices[w + 1]; ++j) { // neighbors of w

                    int u = G.edges[j]; // neigh of w (conflicting with v)
                    if (color[u] > 0) {
//                        if (G.verbose)  printf("\t\t u=%d, color=%d, mark w=%d, color=%d \n",u,color[u],w,color[w]);

                        // mark color of adjacent vertex u
                        forbidden_colors[color[u]] = w;
                    }
                }

                // find a color class C_j for which w does not conflict!
                for (int i = k+1; i < c-1; ++i) {
                    if (forbidden_colors[i] != w) { // check if color i marked with w, if so, then neighbor of u
                        if (k < color[v]) {
                            if (G.verbose) printf("repairing!!! v=%d assigned to color=%d,  w=%d assigned to color=%d \n ",v,color[w],w,c);
                            color[v] = color[w];
                            color[w] = i;
                            return true;
                        }
                    }
                } // for loop searching over adj colors
            } // if single conflict/adj vertex in color class
        } // search over color classes
        return false;
    }


        static void greedy_coloring(
                graphpack_graph &G,
                vector<long long>& vs,
                vector<int>& es,
                vector<Vertex> &P,
                int & lb, // upper bound such as K
                int & ub, // lower bound such as maximum clique or a large clique
                int &num_colors,
                bool is_repair_coloring = false) {

            int max_color = 1;
//            int min_color
            int max_degree = G.max_degree+2;

            // must be initialized to a value that is not a vertex identifier in P
            vector<int> used(max_degree,-1);    // used array is for marking the colors of neighbors
            vector<int> color(G.num_vertices(),0);//+1,0);  // set all vertex colors to zero

    //        G.color.resize(G.num_vertices(),0);

            // assign initial vertex to the first color
            color[P[0].get_id()] = 1;
            printf("[graphpack: greedy coloring]  |P| = %lu, d_max = %d \n", P.size(), G.max_degree);

            if (is_repair_coloring) {
                printf("[graphpack: repair coloring]  \n");
                // color vertices in order of P
                for (int i = 1; i < P.size(); ++i) {
                    int v = P[i].get_id(); // get next vertex in ordering

                    // The color search for assigning a color to a single vertex is as follows:
                    // 1. mark colors of neighbors costing at most O(\dmax)
                    // 2. To find the minimum permissible color, we must iterate over each color c = 1,...,max, terminate
                    for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v+1; //mark color of w in the used array with the label of v.

                    }

                    // assign v to the minimum color class
                    for (int c = 1; c <= max_degree+1; ++c) { // check if color is marked with v or not
                        if (used[c] != v+1) {
                            color[v] = c;   // color of neighbor

//                            if (v < 10) printf("[graphpack: coloring graph! ]   v = %d, color = %d \n",v,color[v]);

                            if (c > max_color) {
                                max_color = c; // compute max

                                if (repair_coloring(G,P,color,used,v,c)) {
                                    max_color--;
                                    //                        printf("*** repaired the coloring *** \n");
                                }
                            }
                            break; // terminate above loop, since color was assigned to v
                        }
                    }
                }
            }
            else {
                printf("[graphpack: normal fast coloring, without recolor procedure]  \n");
                // color vertices in order of P
                for (int i = 1; i < P.size(); i++) {
                    int v = P[i].get_id(); // get next vertex in ordering

                    for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v+1; //mark color of w in the used array with the label of v.

                    }

                    // assign v to the minimum color class
                    for (int c = 1; c <= max_degree+1; ++c) { // check if color is marked with v or not
                        if (used[c] != v+1) {
                            color[v] = c;   // color of neighbor
//                            if (v < 10) printf("[graphpack: coloring graph! ]   v = %d, color = %d \n",v,color[v]);

                            if (c > max_color) max_color = c; // compute max
                            break;
                        }
                    }
                }
            } // is_repair_coloring
            num_colors = max_color;

            if (G.num_vertices() < 10) {
                print_line(80);
                for (int i = 0; i < color.size(); ++i) {
                    printf("[graphpack: color array after coloring] i = %d, color = %d \n",i,color[i]);
                }
                print_line(80);
            }


            double avg_color_class = 0;
            vector<int> color_class_sizes(num_colors+2,0);

            G.color_sizes.clear();
            G.color_sizes.resize(color.size(),0);

//            for (int i = 0; i < G.color.size(); ++i) {
//                color_class_sizes[G.color[i]]++;
//            }

            for (int i = 0; i < color.size(); ++i) {
                G.color_sizes[color[i]]++;
            }

            int max_color_size;
            int min_color_size = max_color_size = G.color_sizes[1];
            for (int i = 1; i < num_colors; ++i) {
                avg_color_class += G.color_sizes[i]; //color_class_sizes[i];
                if (G.color_sizes[i] > max_color_size)
                    max_color_size = G.color_sizes[i];
                if (G.color_sizes[i] < min_color_size)
                    min_color_size = G.color_sizes[i];
            }
            avg_color_class = avg_color_class / num_colors;
            G.max_independent_set_size = max_color_size;
            G.min_independent_set_size = min_color_size;
            G.avg_independent_set_size = avg_color_class;
            printf("\t\t max color class size (independent set) = %d, avg size = %lg \n", G.max_independent_set_size, G.avg_independent_set_size);

            G.color.resize(color.size());
            G.color = color;


            if (G.num_vertices() < 10) {
                print_line(80);
                for (int i = 0; i < G.color.size(); ++i) {
                    printf("[graphpack: G.color array after coloring and copying the colors into G]   i = %d, color = %d \n",i,G.color[i]);
                }
                print_line(80);
            }
        }




        static void greedy_distance_two_coloring(
                graphpack_graph &G,
                vector<long long>& vs,
                vector<int>& es,
                vector<Vertex> &P,
                int & lb, // upper bound such as K
                int & ub, // lower bound such as maximum clique or a large clique
                int &num_colors,
                bool is_repair_coloring = false) {

            int max_color = 1;
            int max_degree = G.max_degree+2;

            // must be initialized to a value that is not a vertex identifier in P
            vector<int> used(max_degree,-1);    // used array is for marking the colors of neighbors
            vector<int> color(G.num_vertices(),0);//+1,0);  // set all vertex colors to zero

    //        G.color.resize(G.num_vertices(),0);

            // assign initial vertex to the first color
            color[P[0].get_id()] = 1;
            printf("[graphpack: greedy coloring]  |P| = %lu, d_max = %d \n", P.size(), G.max_degree);

            if (is_repair_coloring) {
                // color vertices in order of P
                for (int i = 1; i < P.size(); ++i) {
                    int v = P[i].get_id(); // get next vertex in ordering

                    // neighbors of v
                    for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                        int u = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        // mark color of u as used! avoid it!
                        used[color[u]] = v; //mark color of w in the used array with the label of v.

                        // neighbors of u
                        for (long long jj = vs[u]; jj < vs[u + 1]; ++jj) {
                            int w = es[jj];

                            // mark the color of w as used!
                            used[color[w]] = v;
                        }

                    }

                    // assign v to the minimum color class
                    for (int c = 1; c <= max_degree+1; ++c) { // check if color is marked with v or not
                        if (used[c] != v) {
                            color[v] = c;   // color of neighbor

//                            if (v < 10) printf("[graphpack: coloring graph! ]   v = %d, color = %d \n",v,color[v]);

                            if (c > max_color) {
                                max_color = c; // compute max

                                if (repair_coloring(G,P,color,used,v,c)) {
                                    max_color--;
                                    //                        printf("*** repaired the coloring *** \n");
                                }
                            }
                            break; // terminate above loop, since color was assigned to v
                        }
                    }
                }
            }
            else {
                // color vertices in order of P
                for (int i = 1; i < P.size(); i++) {
                    int v = P[i].get_id(); // get next vertex in ordering

//                    for (long long j = vs[v]; j < vs[v + 1]; ++j) {
//                        int w = es[j]; // w is a neighbor of v
//
//                        /** get the color of w, index into used array with it,
//                         *  and mark that position using the label of v     */
//                        used[color[w]] = v; //mark color of w in the used array with the label of v.
//
//                    }


//                    mark_used_colors();
                    // neighbors of v
                    for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                        int u = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        // mark color of u as used! avoid it!
                        used[color[u]] = v; //mark color of w in the used array with the label of v.

                        // neighbors of u
                        for (long long jj = vs[u]; jj < vs[u + 1]; ++jj) {
                            int w = es[jj];

                            // mark the color of w as used!
                            used[color[w]] = v;
                        }

                    }

                    // assign v to the minimum color class
                    for (int c = 1; c <= max_degree+1; ++c) { // check if color is marked with v or not
                        if (used[c] != v) {
                            color[v] = c;   // color of neighbor
//                            if (v < 10) printf("[graphpack: coloring graph! ]   v = %d, color = %d \n",v,color[v]);

                            if (c > max_color) max_color = c; // compute max
                            break;
                        }
                    }
                }
            } // is_repair_coloring
            num_colors = max_color;


            double avg_color_class = 0;
            vector<int> color_class_sizes(num_colors+2,0);

            G.color_sizes.clear();
            G.color_sizes.resize(color.size(),0);

//            for (int i = 0; i < G.color.size(); ++i) {
//                color_class_sizes[G.color[i]]++;
//            }

            for (int i = 0; i < color.size(); ++i) {
                G.color_sizes[color[i]]++;
            }

            int max_color_size;
            int min_color_size = max_color_size = G.color_sizes[1];
            for (int i = 1; i < num_colors; ++i) {
                avg_color_class += G.color_sizes[i]; //color_class_sizes[i];
                if (G.color_sizes[i] > max_color_size)
                    max_color_size = G.color_sizes[i];
                if (G.color_sizes[i] < min_color_size)
                    min_color_size = G.color_sizes[i];
            }
            avg_color_class = avg_color_class / num_colors;
            G.max_independent_set_size = max_color_size;
            G.min_independent_set_size = min_color_size;
            G.avg_independent_set_size = avg_color_class;
            printf("\t\t max color class size (independent set) = %d, avg size = %lg \n", G.max_independent_set_size, G.avg_independent_set_size);

            G.color.resize(G.num_vertices());
            G.color = color;

        }


    static void greedy_coloring_maxsize(
            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            int & lb, // upper bound such as K
            int & ub, // lower bound such as maximum clique or a large clique
            int &num_colors) {

        int max_color = 0;
        int max_degree = G.max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(G.max_degree+1,-1);    // used array is for marking the colors of neighbors
        vector<int> color(G.num_vertices(),0);  // set all vertex colors to zero

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;

        // color vertices in order of P
        for (int i = 1; i < P.size(); i++) {
            int v = P[i].get_id(); // get next vertex in ordering

            for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                int w = es[j]; // w is a neighbor of v

                /** get the color of w, index into used array with it,
                 *  and mark that position using the label of v     */
                used[color[w]] = v; //mark color of w in the used array with the label of v.

            }

            // assign v to the minimum color class
            for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                if (used[c] != v) {
                    color[v] = c;   // color of neighbor
                    if (c > max_color) max_color = c; // compute max

                    break;
                }
            }


        }
        num_colors = max_color;





        double avg_color_class = 0;
        int max_color_size = 0;

        vector<int> color_class_sizes(num_colors+2,0);
        for (int i = 0; i < color.size(); ++i) {
            color_class_sizes[color[i]]++;
        }

        for (int i = 0; i < num_colors; ++i) {
            avg_color_class += color_class_sizes[i];
            if (color_class_sizes[i] > max_color_size)
                max_color_size = color_class_sizes[i];
        }
        avg_color_class = avg_color_class / num_colors;
        printf("max color class size (independent set) = %d \n", max_color_size);
        printf("avg color class size = %lg \n", avg_color_class);

    }


    /**
     * O(|E| + |V|) greedy coloring algorithm
     * Large sparse graphs
     */
    static void greedy_coloring(
            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            int & lb, // upper bound such as K
            int & ub, // lower bound such as maximum clique or a large clique
            int &num_colors,
            vector<int> &color) {

        int max_color = 0;
        int max_degree = G.max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(G.max_degree+1,-1);    // used array is for marking the colors of neighbors
//        vector<int> color(G.num_vertices(),0);  // set all vertex colors to zero

        color.resize(G.num_vertices(),0);

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;

        // color vertices in order of P
        for (int i = 1; i < P.size(); i++) {
            int v = P[i].get_id(); // get next vertex in ordering

            for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                int w = es[j]; // w is a neighbor of v

                /** get the color of w, index into used array with it,
                 *  and mark that position using the label of v     */
                used[color[w]] = v; //mark color of w in the used array with the label of v.

            }

            // assign v to the minimum color class
            for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                if (used[c] != v) {
                    color[v] = c;   // color of neighbor
                    if (c > max_color) max_color = c; // compute max

                    break;
                }
            }


        }
        num_colors = max_color;

    }


    /**
     * O(|E| + |V|) greedy coloring algorithm
     * Large sparse graphs
     */
    static void greedy_coloring_repair(
            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            int & lb, // upper bound such as K
            int & ub, // lower bound such as maximum clique or a large clique
            int &num_colors) {

        int max_color = 0;
#ifdef _DEBUG
        int steps = 0;
#endif

        int max_degree = G.max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(G.max_degree+1,-1);    // used array is for marking the colors of neighbors
        vector<int> color(G.num_vertices(),0);  // set all vertex colors to zero
        vector<int> S(max_degree,0); // num of verts assigned to each color class

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;

        // color vertices in order of P
        for (int i = 1; i < P.size(); i++) {
            int v = P[i].get_id(); // get next vertex in ordering

            for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                int w = es[j]; // w is a neighbor of v

                /** get the color of w, index into used array with it,
                 *  and mark that position using the label of v     */
                used[color[w]] = v; //mark color of w in the used array with the label of v.

            }

            // assign v to the minimum color class
            for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                if (used[c] != v) {
                    color[v] = c;   // color of neighbor
                    S[c]++; // increment size of color class
                    if (c > max_color) {
                        max_color = c; // compute max
                        if (repair_coloring(G,P,color,used,v,c)) {
                            max_color--;
                            printf("*** repaired the coloring *** \n");
                        }

                    }
#ifdef _DEBUG
                    steps++;
#endif
                    break;
                }
            }
        }
        num_colors = max_color;

#ifdef _DEBUG
        cout << "number of steps in greedy coloring = " << steps <<endl;
#endif

    }

    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     *
     * @param v
     * @param k
     * @param colors
     * @param adj
     * @return
     */
    inline
    static bool repair_coloring_fast_neighs(
//            graphpack_graph &G,
            vector<long long>& vertices,
            vector<int>& edges,
            vector<Vertex> &P,
            vector< int >& color,
            vector<int> &used,
            int & v,
            int c,
            int &max_degree,
            vector<int> & ind) {

//        int k = c;
        vector<int> conflicts(max_degree+1,0);


//        if (G.verbose) printf("v=%d, color(v)=%d, c=%d \n", v,color[v],c);

//        vector<int> used_colors(G.max_degree+1,0);
        vector<int> forbidden_colors(max_degree+1,0);

        // mark the colors of neighbors, store their id in that color
        for (long long j = vertices[v]; j < vertices[v + 1]; ++j) {
            if (ind[edges[j]]) { // ensure neighbor is in P

                int w = edges[j]; // w is a neighbor of v

                if (color[w] > 0) {

                    conflicts[color[w]]++;// = conflicts[color[w]] + 1;
                    // get the color of w, index into used array with it, and mark that position using the label of v
                    //                used_colors[color[w]] = w; //mark color of w in the used array with the label of v.

                    used[color[w]] = w; //mark color of w in the used array with the label of v.

                    //                if (G.verbose)
                    //                    printf("\t w=%d, color(w) = %d,  adjacent to %d vertices \n",
                    //                        w,color[w],conflicts[color[w]]);
                }
            }
        }
        // find C_i such that v is adjacent to a single vertex w \in C_i
        for (int k = 1; k < c-1; ++k) { //must start from 1, since 0 is not a color

            if (conflicts[k] == 1) { // if we have a single conflicting color class

                // w MUST BE the single conflicting vertex, hence, only need a single vertex, so only need to mark it in used_colors!
//                int w = used_colors[k];
                int w = used[k];



//                if (G.verbose) printf("v=%d (c=%d) is adj to w=%d, k=%d, |N(w) intersect C_i| = %d \n",v,c, w, k, conflicts[k]);

                // mark the colors of w's neighbors
                for (long long j = vertices[w]; j < vertices[w + 1]; ++j) { // neighbors of w
                    if (ind[edges[j]]) { // ensure neighbor is in P

                        int u = edges[j]; // neigh of w (conflicting with v)
                        // set position of u's color to vertex w, so that if
                        //                    used[color[u]] = w; // get color of u, and use it to index into used array, set id of w to that position



                        if (color[u] > 0) {
                            //                        if (G.verbose)  printf("\t\t u=%d, color=%d, mark w=%d, color=%d \n",u,color[u],w,color[w]);

                            // mark color of adjacent vertex u
                            forbidden_colors[color[u]] = w;
                        }
                    }
                }

                // find a color class C_j for which w does not conflict!
                for (int i = k+1; i < c-1; ++i) {
                    if (forbidden_colors[i] != w) { // check if color i marked with w, if so, then neighbor of u
                        if (k < color[v]) {
//                            if (verbose) printf("repairing!!! v=%d assigned to color=%d,  w=%d assigned to color=%d \n ",v,color[w],w,c);
                            color[v] = color[w];
                            color[w] = i;
                            return true;
                        }
                    }
                } // for loop searching over adj colors
            } // if single conflict/adj vertex in color class
        } // search over color classes
        return false;
    }

    static void greedy_coloring_fast_neighs(
            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
//            vector<int> & color,
//            int & lb, // upper bound such as K
//            int & ub, // lower bound such as maximum clique or a large clique
            int & num_colors,
            int & max_color_class,
            int & max_degree,
            vector<int> & ind,
            bool is_repair_coloring = false) {

        int max_color = 1;
//        max_degree = max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(max_degree+2,-1);    // used array is for marking the colors of neighbors
        vector<int> color(vs.size()-1,0);  // set all vertex colors to zero

//        color.resize(vs.size()-1,0);

//        G.color.resize(G.num_vertices(),0);

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;


        for (int i = 0; i < P.size(); ++i) { // mark vertices in P only!
            ind[P[i].get_id()] = 1;
        }

        if (is_repair_coloring) {
            // color vertices in order of P
            for (int i = 1; i < P.size(); ++i) {
                int v = P[i].get_id(); // get next vertex in ordering

                for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                    if (ind[es[j]]) {

                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v; //mark color of w in the used array with the label of v.
                    }

                }

                // assign v to the minimum color class
                for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                    if (used[c] != v) {
                        color[v] = c;   // color of neighbor
                        if (c > max_color) {
                            max_color = c; // compute max

                            if (repair_coloring_fast_neighs(vs,es,P,color,used,v,c,max_degree,ind)) {
                                max_color--;
                                //                        printf("*** repaired the coloring *** \n");
                            }
                        }
                        break; // terminate above loop, since color was assigned to v
                    }
                }
            }
        }
        else {
            // color vertices in order of P
            for (int i = 1; i < P.size(); i++) {
                int v = P[i].get_id(); // get next vertex in ordering

                for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                    if (ind[es[j]]) {

                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v; //mark color of w in the used array with the label of v.
                    }
                }

                // assign v to the minimum color class
                for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                    if (used[c] != v) {
                        color[v] = c;   // color of neighbor
                        if (c > max_color) max_color = c; // compute max
                        break;
                    }
                }
            }
        } // is_repair_coloring
        num_colors = max_color;

//        G.color.resize(color.size(),0);
//        G.color.swap(color);

        vector<int> color_class_sizes(num_colors+2,0);

        for (int i = 0; i < color.size(); ++i) {
            color_class_sizes[color[i]]++;
//            G.color[i] = color[i];
        }

        max_color_class = color_class_sizes[1];
        for (int i = 1; i < num_colors; ++i) {
            if (color_class_sizes[i] > max_color_class)
                max_color_class = color_class_sizes[i];
        }

//        return max_color_size;
    }

    static void greedy_coloring_fast_neighs(
//            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
//            vector<int> & color,
//            int & lb, // upper bound such as K
//            int & ub, // lower bound such as maximum clique or a large clique
            int & num_colors,
            int & max_color_class,
            int & max_degree,
            vector<int> & ind,
            bool is_repair_coloring = false) {

        int max_color = 1;
//        max_degree = max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(max_degree+2,-1);    // used array is for marking the colors of neighbors
        vector<int> color(vs.size()-1,0);  // set all vertex colors to zero

//        color.resize(vs.size()-1,0);

//        G.color.resize(G.num_vertices(),0);

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;


        for (int i = 0; i < P.size(); ++i) { // mark vertices in P only!
            ind[P[i].get_id()] = 1;
        }

        if (is_repair_coloring) {
            // color vertices in order of P
            for (int i = 1; i < P.size(); ++i) {
                int v = P[i].get_id(); // get next vertex in ordering

                for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                    if (ind[es[j]]) {

                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v; //mark color of w in the used array with the label of v.
                    }

                }

                // assign v to the minimum color class
                for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                    if (used[c] != v) {
                        color[v] = c;   // color of neighbor
                        if (c > max_color) {
                            max_color = c; // compute max

                            if (repair_coloring_fast_neighs(vs,es,P,color,used,v,c,max_degree,ind)) {
                                max_color--;
                                //                        printf("*** repaired the coloring *** \n");
                            }
                        }
                        break; // terminate above loop, since color was assigned to v
                    }
                }
            }
        }
        else {
            // color vertices in order of P
            for (int i = 1; i < P.size(); i++) {
                int v = P[i].get_id(); // get next vertex in ordering

                for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                    if (ind[es[j]]) {

                        int w = es[j]; // w is a neighbor of v

                        /** get the color of w, index into used array with it,
                         *  and mark that position using the label of v     */
                        used[color[w]] = v; //mark color of w in the used array with the label of v.
                    }
                }

                // assign v to the minimum color class
                for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                    if (used[c] != v) {
                        color[v] = c;   // color of neighbor
                        if (c > max_color) max_color = c; // compute max
                        break;
                    }
                }
            }
        } // is_repair_coloring
        num_colors = max_color;

//        G.color = color;

        vector<int> color_class_sizes(num_colors+2,0);

        for (int i = 0; i < color.size(); ++i)
            color_class_sizes[color[i]]++;

        max_color_class = color_class_sizes[1];
        for (int i = 1; i < num_colors; ++i) {
            if (color_class_sizes[i] > max_color_class)
                max_color_class = color_class_sizes[i];
        }

//        return max_color_size;
    }

    /**
     * O(|E| + |V|) greedy coloring algorithm
     * Large sparse graphs
     */
    static void greedy_coloring_desaturation(
            graphpack_graph &G,
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            int & lb, // upper bound such as K
            int & ub, // lower bound such as maximum clique or a large clique
            int &num_colors) {

        int max_color = 0;
        int max_degree = G.max_degree+1;

        // must be initialized to a value that is not a vertex identifier in P
        vector<int> used(G.max_degree+1,-1);    // used array is for marking the colors of neighbors
        vector<int> color(G.num_vertices(),0);  // set all vertex colors to zero

        //        Suppose we are given a graph and a partial coloring of this graph (some of the vertices are already colored).
        //    We define $$deg_s(v)$$ saturation degree of a vertex as the
        // number of different colors in the neighborhood of that vertex.


        // SDO: saturation degree of a vertex v is the number of its adjacent differently colored neighbors
        // IDO: incidence degree ordering is the number of adjacent colored vertices
        vector<int> degree_saturation(G.vertices.size(),0);
        vector<int> cand;
        cand.reserve(max_degree+1);
        vector<int> T;
        T.reserve(max_degree+1);

        // assign initial vertex to the first color
        color[P[0].get_id()] = 1;

        // color vertices in order of P
        for (int i = 1; i < P.size(); i++) {
            int v = P[i].get_id(); // get next vertex in ordering

            for (long long j = vs[v]; j < vs[v + 1]; ++j) {
                int w = es[j]; // w is a neighbor of v
                // color[w] is the color of a neighbor w if color[w] > 0,
                // otherwise w is uncolored (if color[w] = 0)


                if (color[w] > 0) // neighbor w is colored
                    if (used[color[w]] != v) // num of unique colors
                        degree_saturation[v]++;
                //                    else {
                //                        cand.push_back();
                //                    }

                /** get the color of w, index into used array with it,
                 *  and mark that position using the label of v     */
                used[color[w]] = v; //mark color of w in the used array with the label of v.


            }

            // degree saturation: Let deg_{s}(v) be the number of different colors in N(v) -- neighborhood of $v$.
            // number of different colors in the neighborhood of that vertex.


            // assign v to the minimum color class
            for (int c = 1; c < max_degree; ++c) { // check if color is marked with v or not
                if (used[c] != v) {
                    color[v] = c;   // color of neighbor
                    if (c > max_color) max_color = c; // compute max
                    break;
                }
            }

        }
        num_colors = max_color;

    }


    int static DynamicLargestFirstOrdering(graphpack_graph &G, vector<Vertex> &V) {


        int i, u, l;

        int i_HighestInducedVertexDegree;

        int i_VertexCount, i_InducedVertexDegree;

        int i_InducedVertexDegreeCount;

        int i_SelectedVertex, i_SelectedVertexCount;

        vector<int> vi_InducedVertexDegree;

        vector< vector <int> > vvi_GroupedInducedVertexDegree;

        vector< int > vi_VertexLocation;

        i_VertexCount = (signed) G.vertices.size() - 1;

        vi_InducedVertexDegree.clear();
        vi_InducedVertexDegree.reserve((unsigned) i_VertexCount);

        vvi_GroupedInducedVertexDegree.clear();
        vvi_GroupedInducedVertexDegree.resize((unsigned) i_VertexCount);

        vi_VertexLocation.clear();
        vi_VertexLocation.reserve((unsigned) i_VertexCount);

        i_SelectedVertex = -1;

        i_HighestInducedVertexDegree = false;

        for(i=0; i<i_VertexCount; i++)
        {
            //get vertex degree for each vertex
            i_InducedVertexDegree = G.vertices[i+1] - G.vertices[i];

            //vi_InducedVertexDegree[i] = vertex degree of vertex i
            vi_InducedVertexDegree.push_back(i_InducedVertexDegree);

            // vector vvi_GroupedInducedVertexDegree[i] = all the vertices with degree i
            // for every new vertex with degree i, it will be pushed to the back of vector vvi_GroupedInducedVertexDegree[i]
            vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].push_back(i);

            //vi_VertexLocation[i] = location of vertex i in vvi_GroupedInducedVertexDegree[i_InducedVertexDegree]
            vi_VertexLocation.push_back(vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].size() - 1);

            //get max degree (i_HighestInducedVertexDegree)
            if(i_HighestInducedVertexDegree < i_InducedVertexDegree)
            {
                i_HighestInducedVertexDegree = i_InducedVertexDegree;
            }
        }

        vector<long long> ordered;
        ordered.clear();
        ordered.reserve((unsigned) i_VertexCount);

        i_SelectedVertexCount = false;

        // just counting the number of vertices that we have worked with,
        // stop when i_SelectedVertexCount == i_VertexCount, i.e. we have looked through all the vertices
        while(i_SelectedVertexCount < i_VertexCount)
        {
            //pick the vertex with largest degree
            for(i = i_HighestInducedVertexDegree; i >= 0; i--)
            {
                i_InducedVertexDegreeCount = (signed) vvi_GroupedInducedVertexDegree[i].size();

                if(i_InducedVertexDegreeCount != false)
                {
                    i_SelectedVertex = vvi_GroupedInducedVertexDegree[i].back();
                    //remove the i_SelectedVertex from vvi_GroupedInducedVertexDegree
                    vvi_GroupedInducedVertexDegree[i].pop_back();
                    break;
                }
                else
                    i_HighestInducedVertexDegree--;
            }

            //for every D1 neighbor of the i_SelectedVertex, decrease their degree by one and then update their position in vvi_GroupedInducedVertexDegree
            // and vi_VertexLocation
            for(i=G.vertices[i_SelectedVertex]; i<G.vertices[i_SelectedVertex+1]; i++)
            {
                u = G.edges[i];

                if(vi_InducedVertexDegree[u] == -1)
                {
                    continue;
                }

                // move the last element in this bucket to u's position to get rid of expensive erase operation
                if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]].size() > 1)
                {
                    l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]].back();

                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]][vi_VertexLocation[u]] = l;


                    vi_VertexLocation[l] = vi_VertexLocation[u];
                }

                // remove last element from this bucket
                vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]].pop_back();

                // reduce degree of u by 1
                vi_InducedVertexDegree[u]--;

                // move u to appropriate bucket
                vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]].push_back(u);

                // update vi_VertexLocation[u] since it has now been changed
                vi_VertexLocation[u] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegree[u]].size() - 1;
            }

            //Mark the i_SelectedVertex as read (-1), so that we don't look at it again
            vi_InducedVertexDegree[i_SelectedVertex] = -1;

            //Select the vertex by pushing it to the end of ordered
            ordered.push_back(i_SelectedVertex);

            //increment i_SelectedVertexCount
            i_SelectedVertexCount = i_SelectedVertexCount + 1;
        }

        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }

        // clear the buffers
        vi_InducedVertexDegree.clear();
        vi_VertexLocation.clear();
        vvi_GroupedInducedVertexDegree.clear();

        return true;
    }


    //Public Function 1357
    static int DistanceTwoLargestFirstOrdering(graphpack_graph &G, vector<Vertex> &V)
    {

        int i, j, k;

        int i_VertexCount;

        int i_HighestDistanceTwoVertexDegree;

        int i_DistanceTwoVertexDegree, i_DistanceTwoVertexDegreeCount;

        vector<int> vi_IncludedVertices;

        vector< vector<int> > v2i_GroupedDistanceTwoVertexDegree;

        i_HighestDistanceTwoVertexDegree = false;

        i_VertexCount = ((signed) G.vertices.size() - 1);

        v2i_GroupedDistanceTwoVertexDegree.clear();
        v2i_GroupedDistanceTwoVertexDegree.resize((unsigned) i_VertexCount);

        vi_IncludedVertices.clear();
        vi_IncludedVertices.resize((unsigned) i_VertexCount, -1);

        for(i=0; i<i_VertexCount; i++)
        {
            vi_IncludedVertices[i] = i;

            i_DistanceTwoVertexDegree = false;

            for(j=G.vertices[i]; j<G.vertices[i + 1]; j++)
            {
                if(vi_IncludedVertices[G.edges[j]] != i)
                {
                    i_DistanceTwoVertexDegree++;

                    vi_IncludedVertices[G.edges[j]] = i;
                }

                for(k=G.vertices[G.edges[j]]; k<G.vertices[G.edges[j] + 1]; k++)
                {
                    if(vi_IncludedVertices[G.edges[k]] != i)
                    {
                        i_DistanceTwoVertexDegree++;

                        vi_IncludedVertices[G.edges[k]] = i;
                    }
                }
            }

            v2i_GroupedDistanceTwoVertexDegree[i_DistanceTwoVertexDegree].push_back(i);

            if(i_HighestDistanceTwoVertexDegree < i_DistanceTwoVertexDegree)
            {
                i_HighestDistanceTwoVertexDegree = i_DistanceTwoVertexDegree;
            }
        }
        vector<long long> ordered;
        ordered.clear();
        ordered.reserve((unsigned) i_VertexCount);

        for(i=i_HighestDistanceTwoVertexDegree; i>=0; i--)
        {
            i_DistanceTwoVertexDegreeCount = (signed) v2i_GroupedDistanceTwoVertexDegree[i].size();

            for(j=0; j<i_DistanceTwoVertexDegreeCount; j++)
            {
                ordered.push_back(v2i_GroupedDistanceTwoVertexDegree[i][j]);
            }
        }

        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }

        vi_IncludedVertices.clear();
        v2i_GroupedDistanceTwoVertexDegree.clear();

        return(true);
    }



    static int DistanceTwoDynamicLargestFirstOrdering(graphpack_graph &G, vector<Vertex> &V)
    {
        int i, j, k, l, u, v;

        int i_HighestInducedVertexDegree;

        int i_VertexCount, i_InducedVertexDegree;

        int i_InducedVertexDegreeCount;

        int i_SelectedVertex, i_SelectedVertexCount;

        vector < int > vi_IncludedVertices;

        vector < int > vi_InducedVertexDegrees;

        vector < vector < int > > vvi_GroupedInducedVertexDegree;

        vector < int > vi_VertexLocations;

        i_VertexCount = ((signed) G.vertices.size() - 1);

                vi_IncludedVertices.clear();
                vi_IncludedVertices.resize((unsigned) i_VertexCount, -1);

        vi_InducedVertexDegrees.clear();
        vi_InducedVertexDegrees.reserve((unsigned) i_VertexCount);

        vvi_GroupedInducedVertexDegree.clear();
        vvi_GroupedInducedVertexDegree.resize((unsigned) i_VertexCount);

        vi_VertexLocations.clear();
        vi_VertexLocations.reserve((unsigned) i_VertexCount);


        i_SelectedVertex = -1;

        i_HighestInducedVertexDegree = false;

        for(i=0; i<i_VertexCount; i++)
        {
            vi_IncludedVertices[i] = i;

            i_InducedVertexDegree = false;

            for(j=G.vertices[i]; j<G.vertices[i + 1]; j++)
            {
                if(vi_IncludedVertices[G.edges[j]] != i)
                {
                    i_InducedVertexDegree++;

                    vi_IncludedVertices[G.edges[j]] = i;
                }

                for(k=G.vertices[G.edges[j]]; k<G.vertices[G.edges[j] + 1]; k++)
                {
                    if (vi_IncludedVertices[G.edges[k]] != i)
                    {
                        i_InducedVertexDegree++;

                        vi_IncludedVertices[G.edges[k]] = i;
                    }
                }
            }

            vi_InducedVertexDegrees.push_back(i_InducedVertexDegree);

            vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].push_back(i);

            vi_VertexLocations.push_back(vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].size() - 1);

            if(i_HighestInducedVertexDegree < i_InducedVertexDegree)
            {
                i_HighestInducedVertexDegree = i_InducedVertexDegree;
            }
        }

        vector<long long> ordered;
        ordered.clear();
        ordered.reserve((unsigned) i_VertexCount);

        vi_IncludedVertices.assign((unsigned) i_VertexCount, -1);

        i_SelectedVertexCount = false;

        while(i_SelectedVertexCount < i_VertexCount)
        {
            for(i=i_HighestInducedVertexDegree; i >= 0; i--)
            {
                i_InducedVertexDegreeCount = (signed) vvi_GroupedInducedVertexDegree[i].size();

                if(i_InducedVertexDegreeCount != false)
                {
                    i_SelectedVertex = vvi_GroupedInducedVertexDegree[i].back();
                    vvi_GroupedInducedVertexDegree[i].pop_back();
                    break;
                }
                else
                    i_HighestInducedVertexDegree--;

            }

            vi_IncludedVertices[i_SelectedVertex] = i_SelectedVertex;

            for(i=G.vertices[i_SelectedVertex]; i<G.vertices[i_SelectedVertex + 1]; i++)
            {
                u = G.edges[i];

                if(vi_InducedVertexDegrees[u] == -1)
                {
                    continue;
                }

                if(vi_IncludedVertices[u] != i_SelectedVertex)
                {
                    // move the last element in this bucket to u's position to get rid of expensive erase operation
                    if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() > 1)
                    {
                        l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].back();
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]][vi_VertexLocations[u]] = l;
                        vi_VertexLocations[l] = vi_VertexLocations[u];
                    }

                    // remove last element from this bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].pop_back();

                    // reduce degree of u by 1
                    vi_InducedVertexDegrees[u]--;

                    // move u to appropriate bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].push_back(u);

                    // update vi_VertexLocation[u] since it has now been changed
                                    vi_VertexLocations[u] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() - 1;

                    // this neighbour has been visited
                    vi_IncludedVertices[u] = i_SelectedVertex;
                }

                for(j=G.vertices[G.edges[i]]; j<G.vertices[G.edges[i] + 1]; j++)
                {
                    v = G.edges[j];

                    if(vi_InducedVertexDegrees[v] == -1)
                    {
                        continue;
                    }

                    if(vi_IncludedVertices[v] != i_SelectedVertex)
                    {
                        // move the last element in this bucket to v's position to get rid of expensive erase operation
                        if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() > 1)
                        {
                            l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].back();
                            vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]][vi_VertexLocations[v]] = l;
                            vi_VertexLocations[l] = vi_VertexLocations[v];
                        }

                        // remove last element from this bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].pop_back();

                        // reduce degree of v by 1
                        vi_InducedVertexDegrees[v]--;

                        // move v to appropriate bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].push_back(v);

                        // update vi_VertexLocation[v] since it has now been changed
                                        vi_VertexLocations[v] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() - 1;

                        // this neighbour has been visited
                        vi_IncludedVertices[v] = i_SelectedVertex;
                    }
                }
            }

            vi_InducedVertexDegrees[i_SelectedVertex] = -1;
            ordered.push_back(i_SelectedVertex);
            i_SelectedVertexCount = i_SelectedVertexCount + 1;
        }




        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }

        vi_IncludedVertices.clear();
                vi_InducedVertexDegrees.clear();
                vvi_GroupedInducedVertexDegree.clear();
                vi_VertexLocations.clear();

        return(true);
    }


    //Public Function 1359
    static int DistanceTwoSmallestLastOrdering(graphpack_graph &G, vector<Vertex> &V)
    {

        int i, j, k, l, u, v;

        int i_HighestInducedVertexDegree;

        int i_VertexCount, i_InducedVertexDegree;

        int i_VertexCountMinus1;

        int i_InducedVertexDegreeCount;

        int i_SelectedVertex, i_SelectedVertexCount;

        vector < int > vi_IncludedVertices;

        vector < int > vi_InducedVertexDegrees;

        vector < vector < int > > vvi_GroupedInducedVertexDegree;

        vector < int > vi_VertexLocations;

        i_VertexCount = ((signed) G.vertices.size() - 1);
        i_VertexCountMinus1 = i_VertexCount - 1; // = i_VertexCount - 1, used when inserting selected vertices into ordered

                vi_IncludedVertices.clear();
                vi_IncludedVertices.resize((unsigned) i_VertexCount, -1);

        vi_InducedVertexDegrees.clear();
        vi_InducedVertexDegrees.reserve((unsigned) i_VertexCount);

        vvi_GroupedInducedVertexDegree.clear();
        vvi_GroupedInducedVertexDegree.resize((unsigned) i_VertexCount);

        vi_VertexLocations.clear();
        vi_VertexLocations.reserve((unsigned) i_VertexCount);


        i_SelectedVertex = -1;

        i_HighestInducedVertexDegree = false;

        for(i=0; i<i_VertexCount; i++)
        {
            vi_IncludedVertices[i] = i;

            i_InducedVertexDegree = false;

            for(j=G.vertices[i]; j<G.vertices[i + 1]; j++)
            {
                if(vi_IncludedVertices[G.edges[j]] != i)
                {
                    i_InducedVertexDegree++;

                    vi_IncludedVertices[G.edges[j]] = i;
                }

                for(k=G.vertices[G.edges[j]]; k<G.vertices[G.edges[j] + 1]; k++)
                {
                    if(vi_IncludedVertices[G.edges[k]] != i)
                    {
                        i_InducedVertexDegree++;

                        vi_IncludedVertices[G.edges[k]] = i;
                    }
                }
            }

            vi_InducedVertexDegrees.push_back(i_InducedVertexDegree);

            vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].push_back(i);

            vi_VertexLocations.push_back(vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].size() - 1);

            if(i_HighestInducedVertexDegree < i_InducedVertexDegree)
            {
                i_HighestInducedVertexDegree = i_InducedVertexDegree;
            }
        }

        vector<long long> ordered;
        ordered.clear();
        ordered.resize((unsigned) i_VertexCount, -1);

        vi_IncludedVertices.assign((unsigned) i_VertexCount, -1);

        i_SelectedVertexCount = false;

        int iMin = 1;

        while(i_SelectedVertexCount < i_VertexCount)
        {
            if(iMin != 0 && vvi_GroupedInducedVertexDegree[iMin -1].size() != false)
                iMin--;

            for(i= iMin; i < i_HighestInducedVertexDegree + 1; i++)
            {
                i_InducedVertexDegreeCount = (signed) vvi_GroupedInducedVertexDegree[i].size();

                if(i_InducedVertexDegreeCount != false)
                {
                    i_SelectedVertex = vvi_GroupedInducedVertexDegree[i].back();
                    vvi_GroupedInducedVertexDegree[i].pop_back();
                    break;
                }
                else
                    iMin++;
            }

            vi_IncludedVertices[i_SelectedVertex] = i_SelectedVertex;

            for(i=G.vertices[i_SelectedVertex]; i<G.vertices[i_SelectedVertex + 1]; i++)
            {
                u = G.edges[i];

                if(vi_InducedVertexDegrees[u] == -1)
                {
                    continue;
                }

                if(vi_IncludedVertices[u] != i_SelectedVertex)
                {
                    // move the last element in this bucket to u's position to get rid of expensive erase operation
                    if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() > 1)
                    {
                        l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].back();
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]][vi_VertexLocations[u]] = l;
                        vi_VertexLocations[l] = vi_VertexLocations[u];
                    }

                    // remove last element from this bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].pop_back();

                    // reduce degree of u by 1
                    vi_InducedVertexDegrees[u]--;

                    // move u to appropriate bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].push_back(u);

                    // update vi_VertexLocation[u] since it has now been changed
                                    vi_VertexLocations[u] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() - 1;

                    // this neighbour has been visited
                    vi_IncludedVertices[u] = i_SelectedVertex;
                }

                for(j=G.vertices[G.edges[i]]; j<G.vertices[G.edges[i] + 1]; j++)
                {
                    v = G.edges[j];

                    if(vi_InducedVertexDegrees[v] == -1)
                    {
                        continue;
                    }

                    if(vi_IncludedVertices[v] != i_SelectedVertex)
                    {
                        // move the last element in this bucket to v's position to get rid of expensive erase operation
                        if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() > 1)
                        {
                            l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].back();
                            vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]][vi_VertexLocations[v]] = l;
                            vi_VertexLocations[l] = vi_VertexLocations[v];
                        }

                        // remove last element from this bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].pop_back();

                        // reduce degree of v by 1
                        vi_InducedVertexDegrees[v]--;

                        // move v to appropriate bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].push_back(v);

                        // update vi_VertexLocation[v] since it has now been changed
                                        vi_VertexLocations[v] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() - 1;

                        // this neighbour has been visited
                        vi_IncludedVertices[v] = i_SelectedVertex;
                    }
                }
            }

            vi_InducedVertexDegrees[i_SelectedVertex] = -1;
            //ordered.push_back(i_SelectedVertex);
            ordered[i_VertexCountMinus1 - i_SelectedVertexCount] = i_SelectedVertex;
            i_SelectedVertexCount = i_SelectedVertexCount + 1;
        }


        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }



        vi_IncludedVertices.clear();
                vi_InducedVertexDegrees.clear();
                vvi_GroupedInducedVertexDegree.clear();
                vi_VertexLocations.clear();

        return(true);
    }


    //Public Function 1360
        static int IncidenceDegreeOrdering(graphpack_graph &G, vector<Vertex> &V)
        {
                int i, u, v, l;

                int i_HighestDegreeVertex, i_MaximumVertexDegree;

        int i_VertexCount, i_VertexDegree;

                int i_IncidenceVertexDegree, i_IncidenceVertexDegreeCount;

                int i_SelectedVertex, i_SelectedVertexCount;

                vector< int > vi_IncidenceVertexDegree;

                vector< vector< int > > vvi_GroupedIncidenceVertexDegree;

                vector< int > vi_VertexLocation;

                i_VertexCount = ((signed) G.vertices.size() - 1);

                vi_IncidenceVertexDegree.clear();
        vi_IncidenceVertexDegree.reserve((unsigned) i_VertexCount);

                vvi_GroupedIncidenceVertexDegree.clear();
                vvi_GroupedIncidenceVertexDegree.resize((unsigned) i_VertexCount);

                vi_VertexLocation.clear();
        vi_VertexLocation.reserve((unsigned) i_VertexCount);

                i_HighestDegreeVertex = i_MaximumVertexDegree = -1;

                i_SelectedVertex = -1;

                i_IncidenceVertexDegree = false;


        // initilly push all the vertices into the first bucket assuming that IncidenceVertexDegree is all 0
        vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree].reserve((unsigned) i_VertexCount); // ONLY FOR THE FIRST BUCKET SINCE WE KNOW in THIS case


                for(i=0; i<i_VertexCount; i++)
                {
            // i_IncidenceVertexDegree is 0 and insert that
            vi_IncidenceVertexDegree.push_back(i_IncidenceVertexDegree);

            // insert vertex i into bucket vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree]
                        vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree].push_back(i);

            // store the location
            vi_VertexLocation.push_back(vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree].size() - 1);

            // calculate the degree
            i_VertexDegree = G.vertices[i + 1] - G.vertices[i];

            // get the max degree vertex
                        if(i_MaximumVertexDegree < i_VertexDegree)
                        {
                                i_MaximumVertexDegree = i_VertexDegree;

                                i_HighestDegreeVertex = i;
                        }
                }

        // reserver memory for ordered
        vector<long long> ordered;
        ordered.clear();
        ordered.reserve((unsigned) i_VertexCount);

        i_SelectedVertexCount = false;

        // NOW SWAP THE MAX DEGREE VERTEX WITH THE LAST VERTEX IN THE FIRST BUCKET
        l = vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree].size() - 1;
        v = vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree][l];
        //u = vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree][i_HighestDegreeVertex];
        u = vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree][vi_VertexLocation[i_HighestDegreeVertex]];

        swap(vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree][vi_VertexLocation[i_HighestDegreeVertex]], vvi_GroupedIncidenceVertexDegree[i_IncidenceVertexDegree][l]);
        swap(vi_VertexLocation[v], vi_VertexLocation[u]);

        int iMax = i_MaximumVertexDegree - 1;
        // just counting the number of vertices that we have worked with,
        // stop when i_SelectedVertexCount == i_VertexCount, i.e. we have looked through all the vertices
        while(i_SelectedVertexCount < i_VertexCount)
                {

                        if(iMax != i_MaximumVertexDegree && vvi_GroupedIncidenceVertexDegree[iMax + 1].size() != false)
                                iMax++;

            //pick the vertex with maximum incidence degree
            for(i=iMax; i>=0; i--)
                        {
                            i_IncidenceVertexDegreeCount = (signed) vvi_GroupedIncidenceVertexDegree[i].size();

                                if(i_IncidenceVertexDegreeCount != false)
                                {
                                    i_SelectedVertex = vvi_GroupedIncidenceVertexDegree[i].back();
                    // remove i_SelectedVertex from  vvi_GroupedIncidenceVertexDegree[i]
                    vvi_GroupedIncidenceVertexDegree[i].pop_back();
                                        break;
                            }
                else
                    iMax--;
                        }

            //for every D1 neighbor of the i_SelectedVertex, decrease their degree by one and then update their position in vvi_GroupedInducedVertexDegree
            // and vi_VertexLocation
            for(i=G.vertices[i_SelectedVertex]; i<G.vertices[i_SelectedVertex + 1]; i++)
                        {
                                u = G.edges[i];

                                if(vi_IncidenceVertexDegree[u] == -1)
                                {
                                        continue;
                                }

                        // move the last element in this bucket to u's position to get rid of expensive erase operation
                if(vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]].size() > 1)
                                {
                                        l = vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]].back();

                                        vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]][vi_VertexLocation[u]] = l;

                                        vi_VertexLocation[l] = vi_VertexLocation[u];
                                }

                // remove the last element from vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]]
                                vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]].pop_back();

                // increase incidence degree of u
                vi_IncidenceVertexDegree[u]++;

                // insert u into appropriate bucket
                                vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]].push_back(u);

                // update location of u
                                vi_VertexLocation[u] = vvi_GroupedIncidenceVertexDegree[vi_IncidenceVertexDegree[u]].size() - 1;
                        }

            //Mark the i_SelectedVertex as read, so that we don't look at it again
            vi_IncidenceVertexDegree[i_SelectedVertex] = -1;
            // insert i_SelectedVertex into ordered
            ordered.push_back(i_SelectedVertex);
            // increament i_SelectedVertexCount
            i_SelectedVertexCount = i_SelectedVertexCount + 1;
        }




        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }




        // clear the buffer
                vi_IncidenceVertexDegree.clear();
                vi_VertexLocation.clear();
                vvi_GroupedIncidenceVertexDegree.clear();

        return(true);
    }


    //Public Function 1361
    static int DistanceTwoIncidenceDegreeOrdering(graphpack_graph & G, vector<Vertex> &V)
    {

        int i, j, k, l, u, v;

        //int i_HighestInducedVertexDegree;
        int i_DistanceTwoVertexDegree;

        int i_HighestDistanceTwoDegreeVertex, i_HighestDistanceTwoVertexDegree;

        int i_VertexCount, i_InducedVertexDegree;

        int i_InducedVertexDegreeCount;

        int i_SelectedVertex, i_SelectedVertexCount;

        vector < int > vi_IncludedVertices;

        vector < int > vi_InducedVertexDegrees;

        vector < vector < int > > vvi_GroupedInducedVertexDegree;

        vector < int > vi_VertexLocations;

        i_VertexCount = ((signed) G.vertices.size() - 1);

                vi_IncludedVertices.clear();
                vi_IncludedVertices.resize((unsigned) i_VertexCount, -1);

        vi_InducedVertexDegrees.clear();
        vi_InducedVertexDegrees.reserve((unsigned) i_VertexCount);

        vvi_GroupedInducedVertexDegree.clear();
        vvi_GroupedInducedVertexDegree.resize((unsigned) i_VertexCount);

        vi_VertexLocations.clear();
        vi_VertexLocations.reserve((unsigned) i_VertexCount);

        i_SelectedVertex = -1;

        i_HighestDistanceTwoDegreeVertex = i_HighestDistanceTwoVertexDegree = -1;
        i_InducedVertexDegree = false;

        // initilly push all the vertices into the first bucket assuming that IncidenceVertexDegree is all 0
        vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].reserve((unsigned) i_VertexCount); // ONLY FOR THE FIRST BUCKET SINCE WE KNOW in THIS case

        for(i=0; i<i_VertexCount; i++)
        {
                        vi_InducedVertexDegrees.push_back(i_InducedVertexDegree);

                        vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].push_back(i);

                        vi_VertexLocations.push_back(vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].size() - 1);

            vi_IncludedVertices[i] = i;

            i_DistanceTwoVertexDegree = false;

            for(j=G.vertices[i]; j<G.vertices[i + 1]; j++)
            {
                if(vi_IncludedVertices[G.edges[j]] != i)
                {
                    i_DistanceTwoVertexDegree++;

                    vi_IncludedVertices[G.edges[j]] = i;
                }

                for(k=G.vertices[G.edges[j]]; k<G.vertices[G.edges[j] + 1]; k++)
                {
                    if(vi_IncludedVertices[G.edges[k]] != i)
                    {
                        i_DistanceTwoVertexDegree++;

                        vi_IncludedVertices[G.edges[k]] = i;
                    }
                }
            }

            if(i_HighestDistanceTwoVertexDegree < i_DistanceTwoVertexDegree)
            {
                i_HighestDistanceTwoVertexDegree = i_DistanceTwoVertexDegree;
                i_HighestDistanceTwoDegreeVertex = i;
            }
        }

        vector<long long> ordered;
        ordered.clear();
        ordered.reserve((unsigned) i_VertexCount);

        vi_IncludedVertices.assign((unsigned) i_VertexCount, -1);


        // NOW SWAP THE MAX DEGREE VERTEX WITH THE LAST VERTEX IN THE FIRST BUCKET
        l = vvi_GroupedInducedVertexDegree[i_InducedVertexDegree].size() - 1;
        v = vvi_GroupedInducedVertexDegree[i_InducedVertexDegree][l];
        u = vvi_GroupedInducedVertexDegree[i_InducedVertexDegree][vi_VertexLocations[i_HighestDistanceTwoDegreeVertex]];
        swap(vvi_GroupedInducedVertexDegree[i_InducedVertexDegree][vi_VertexLocations[i_HighestDistanceTwoDegreeVertex]], vvi_GroupedInducedVertexDegree[i_InducedVertexDegree][l]);
        swap(vi_VertexLocations[v], vi_VertexLocations[u]);

        i_SelectedVertexCount = false;

        int iMax = i_HighestDistanceTwoVertexDegree - 1;

        while(i_SelectedVertexCount < i_VertexCount)
        {
                        if(iMax != i_HighestDistanceTwoVertexDegree && vvi_GroupedInducedVertexDegree[iMax + 1].size() != false)
                                iMax++;

            for(i= iMax; i>= 0; i--)
            {
                i_InducedVertexDegreeCount = (signed) vvi_GroupedInducedVertexDegree[i].size();

                if(i_InducedVertexDegreeCount != false)
                {
                    i_SelectedVertex = vvi_GroupedInducedVertexDegree[i].back();
                    vvi_GroupedInducedVertexDegree[i].pop_back();
                    break;
                }
                else
                    iMax--;
            }

            vi_IncludedVertices[i_SelectedVertex] = i_SelectedVertex;

            for(i=G.vertices[i_SelectedVertex]; i<G.vertices[i_SelectedVertex + 1]; i++)
            {
                u = G.edges[i];

                if(vi_InducedVertexDegrees[u] == -1)
                {
                    continue;
                }

                if(vi_IncludedVertices[u] != i_SelectedVertex)
                {
                    // move the last element in this bucket to u's position to get rid of expensive erase operation
                    if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() > 1)
                    {
                        l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].back();
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]][vi_VertexLocations[u]] = l;
                        vi_VertexLocations[l] = vi_VertexLocations[u];
                    }

                    // remove last element from this bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].pop_back();

                    // reduce degree of u by 1
                    vi_InducedVertexDegrees[u]++;

                    // move u to appropriate bucket
                    vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].push_back(u);

                    // update vi_VertexLocation[u] since it has now been changed
                                    vi_VertexLocations[u] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[u]].size() - 1;

                    // this neighbour has been visited
                    vi_IncludedVertices[u] = i_SelectedVertex;
                }

                for(j=G.vertices[G.edges[i]]; j<G.vertices[G.edges[i] + 1]; j++)
                {
                    v = G.edges[j];

                    if(vi_InducedVertexDegrees[v] == -1)
                    {
                        continue;
                    }

                    if(vi_IncludedVertices[v] != i_SelectedVertex)
                    {
                        // move the last element in this bucket to v's position to get rid of expensive erase operation
                        if(vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() > 1)
                        {
                            l = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].back();
                            vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]][vi_VertexLocations[v]] = l;
                            vi_VertexLocations[l] = vi_VertexLocations[v];
                        }

                        // remove last element from this bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].pop_back();

                        // reduce degree of v by 1
                        vi_InducedVertexDegrees[v]++;

                        // move v to appropriate bucket
                        vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].push_back(v);

                        // update vi_VertexLocation[v] since it has now been changed
                                        vi_VertexLocations[v] = vvi_GroupedInducedVertexDegree[vi_InducedVertexDegrees[v]].size() - 1;

                        // this neighbour has been visited
                        vi_IncludedVertices[v] = i_SelectedVertex;
                    }
                }
            }

            vi_InducedVertexDegrees[i_SelectedVertex] = -1;
            ordered.push_back(i_SelectedVertex);
            i_SelectedVertexCount = i_SelectedVertexCount + 1;
        }


        V.clear();
        V.reserve(G.num_vertices());
        // assumed to be sorted
        for (int i = 0; i < G.num_vertices(); ++i) {
            int u = ordered[i];
            V.push_back(Vertex(u,G.num_vertices()-i));
        }







        vi_IncludedVertices.clear();
                vi_InducedVertexDegrees.clear();
                vvi_GroupedInducedVertexDegree.clear();
                vi_VertexLocations.clear();

        return(true);
    }


    /**
     * @brief
 Searches a color class for a vertex w that only has a
     *              single conflict with v. The vertex w is exchanged with
     *              another suitable candidate, thereby reducing the number
     *              of colors by one.
     *
     * @param v
     * @param k             color class to search
     * @param colors        two dimensional array where colors[k][i]
     *                      represents the ith vertex in the k color class
     * @param adj           adj graph representation
     * @param w_pos         ptr to w in the k color class
     *
     * @return vertex to exchange
     */
    inline
    static int get_conflict(int v, int k, vector< vector<int> >& colors, bool** & adj, int & w_pos) {
        int conflictVar = -1;
        int count = 0;

        for (int i = 0; i < colors[k].size(); i++) { // check vertices in the k color class
            int w = colors[k][i];   // get a vertex
            if (adj[v][w]) {        // edge between w and v, i.e., conflict!
                conflictVar = w;
                w_pos = i;
                count++;            // count how many conflicts v has!
                if (count > 1)  return -1;  // if more than one, ignore
            }
        }
        return conflictVar;
    }


    /**
     * @brief       Find a color class for the vertex such that
     *              the vertex is not adjacent to any other vertex
     *              in that color class.
     */
    inline
    static bool uncolored(int v, int k, vector< vector<int> >& colors, bool** & adj) {

        for (int h = 0; h < colors[k].size(); h++)
            if (adj[v][colors[k][h]])  return true;
        return false;
    }


    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     *
     * @param v
     * @param k
     * @param colors
     * @param adj
     * @return
     */
    static bool repair(int v, int k, vector< vector<int> >& colors, bool** & adj) {

        //! check each color class
        //! todo: more likely to encounter suitable vertices from reverse order
        for (int color_class = k - 1; color_class >= 0; color_class--) {
            //        for (int color_class = 0; color_class < k - 1; color_class++) {

            /* find a color class where v has only a single conflict!
             * w is the vertex conflicting with u  */
            int w_pos = 0; // holds the position of w in colors[k]
            int w = get_conflict(v,color_class,colors,adj,w_pos);

            if (w >= 0) {
                //! check each color class
                for (int col_j = color_class + 1; col_j < k; col_j++) {
                    //! find a new color class for w, then swap w and v!
                    if (!uncolored(w,col_j,colors,adj)) {
                        colors[k].clear(); // remove v from kth color class
                        colors[color_class][w_pos] = v;
                        colors[col_j].push_back(w);
                        return true;
                    }
                }
            }
        }
        return false;
    }

    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     *
     * @param v
     * @param k
     * @param colors
     * @param adj
     * @return
     */
    static bool repair_forward(int v, int k, vector< vector<int> >& colors, bool** & adj) {

        //! check each color class
        //! todo: more likely to encounter suitable vertices from reverse order
        for (int color_class = 0; color_class < k-1; ++color_class) {
            //        for (int color_class = 0; color_class < k - 1; color_class++) {

            /* find a color class where v has only a single conflict!
             * w is the vertex conflicting with u  */
            int w_pos = 0; // holds the position of w in colors[k]
            int w = get_conflict(v,color_class,colors,adj,w_pos);

            if (w >= 0) {
                //! check each color class
                for (int col_j = color_class + 1; col_j < k; ++col_j) {
                    //! find a new color class for w, then swap w and v!
                    if (!uncolored(w,col_j,colors,adj)) {
                        colors[k].clear(); // remove v from kth color class
                        colors[color_class][w_pos] = v;
                        colors[col_j].push_back(w);
                        return true;
                    }
                }
            }
        }
        return false;
    }


    static void greedy_neighborhood_coloring_repair(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<Vertex> & ColOrd,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        for (int w = 0; w < ColOrd.size(); w++) {
            u = ColOrd[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below


            while (uncolored(u,k,colors,adj))  k++;



            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            // only try to repair if vertex caused a new color class to be added
            if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,adj))  max_k--;


            //            if (k < min_k) {
            //                P[j].set_id(u);
            //                j++;
            //            }
        }

        //        if (j > 0)  P[j-1].set_bound(0);
        //        if (min_k <= 0)  min_k = 1;

        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }

        // start at the min color class
        //        for (k = min_k; k <= max_k; k++)
        //            for (int w = 0; w < colors[k].size(); w++) {
        //                P[j].set_id(colors[k][w]);
        //                P[j].set_bound(k);
        //                j++;
        //            }
    }


    static void greedy_neighborhood_coloring_repair_only(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        for (int w = 0; w < P.size(); w++) {
            u = P[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below


            while (uncolored(u,k,colors,adj))  k++;

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            // only try to repair if vertex caused a new color class to be added
            if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,adj))  max_k--;
        }


        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }

    static void greedy_neighborhood_coloring_repair_only_forward(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& num_colors,
            bool** &adj,
            bool is_enum = false,
            bool is_repair = true) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

//        int min_k = mc - C.size() + 1;
//        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        if (is_repair) {
            for (int w = 0; w < P.size(); w++) {
                u = P[w].get_id(); // use initial order
                k = 1;                  // color class for vertex u, so far...
                k_prev = 0;             // used to check of k changes in the while loop below


                while (uncolored(u,k,colors,adj))  k++;

                // start new color class
                if (k > max_k) {
                    max_k = k;
                    colors[max_k+1].clear();
                }

                // add vertex to color class
                colors[k].push_back(u);

                // only try to repair if vertex caused a new color class to be added
//                if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,adj))  max_k--;
                if (colors[k].size() == 1 && repair(u,k,colors,adj))  max_k--;
            }
        }
        else {
            for (int w = 0; w < P.size(); w++) {
                u = P[w].get_id(); // use initial order
                k = 1;                  // color class for vertex u, so far...
                k_prev = 0;             // used to check of k changes in the while loop below


                while (uncolored(u,k,colors,adj))  k++;

                // start new color class
                if (k > max_k) {
                    max_k = k;
                    colors[max_k+1].clear();
                }

                // add vertex to color class
                colors[k].push_back(u);

                // only try to repair if vertex caused a new color class to be added
//                if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,adj))  max_k--;
            }
        }

        num_colors = max_k;


//        for (k = 1; k <= max_k; k++)
//            for (int w = 0; w < colors[k].size(); w++) {
//                P[j].set_id(colors[k][w]);
//                P[j].set_bound(k);
//                j++;
//            }
    }

    /**
     * Only STATIC Ordering
     */
    static void greedy_neighborhood_coloring_static(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<Vertex> & ColOrd,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();

        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        for (int w = 0; w < ColOrd.size(); w++) {
            u = ColOrd[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below

            while (k > k_prev) { // stop when we've found a color class
                k_prev = k;
                /*
                 * While a vertex w in C_k (in order) does not have an edge to u,
                 * If u does not have an edge to any vertex in C_k, then add u to C_k,
                 * Else try the next color class!
                 */
                for (int i = 0; i < colors[k].size(); i++) {
                    if (adj[u][colors[k][i]]) { // if u has an edge with a vertex in C_k
                        k++;                    // then set k=k+1 and break out of for loop!
                        break;
                    }
                }
            }

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            //            if (k < min_k) {
            //                P[j].set_id(u);
            //                j++;
            //            }
        }

        //        if (j > 0)  P[j-1].set_bound(0);
        //        if (min_k <= 0)  min_k = 1;

        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }

        // start at the min color class
        //        for (k = min_k; k <= max_k; k++)
        //            for (int w = 0; w < colors[k].size(); w++) {
        //                P[j].set_id(colors[k][w]);
        //                P[j].set_bound(k);
        //                j++;
        //            }
    }

    // sequential dynamic greedy coloring and sort
    static void greedy_neighborhood_coloring(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size() + 1;
        if (is_enum) { min_k = mc - C.size(); }

        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        for (int w=0; w < P.size(); w++) {
            // get the
            u = P[w].get_id();
            k = 1;      // color class for vertex u, so far...
            k_prev = 0; // used to check of k changes in the while loop below

            while (k > k_prev) { // stop when we've found a color class
                k_prev = k;
                /*
                 * While a vertex w in C_k (in order) does not have an edge to u,
                 * If u does not have an edge to any vertex in C_k, then add u to C_k,
                 * Else try the next color class!
                 */
                for (int i = 0; i < colors[k].size(); i++) {
                    if (adj[u][colors[k][i]]) { // if u has an edge with a vertex in C_k
                        k++;                    // then set k=k+1 and break out of for loop!
                        break;
                    }
                }
            }

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            // does this screw up the vertex_lookup table?
            // may need to adjust it as well.
            if (k < min_k) {
                P[j].set_id(u);
                j++;
            }
        }

        if (j > 0)  P[j-1].set_bound(0);
        if (min_k <= 0)  min_k = 1;

        // start at the min color class
        for (k = min_k; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }




    // sequential dynamic greedy coloring and sort
    static void neigh_coloring_bound(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<short>& ind,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int* pruned,
            int& mc) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size() + 1;
        colors[1].clear();   colors[2].clear();

        for (int w=0; w < P.size(); w++) {
            u = P[w].get_id();
            k = 1, k_prev = 0;

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

            while (k > k_prev) {
                k_prev = k;
                for (int i = 0; i < colors[k].size(); i++) {
                    if (ind[colors[k][i]]) {
                        k++;
                        break;
                    }
                }
            }

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;

            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            colors[k].push_back(u);
            if (k < min_k) {
                P[j].set_id(u);
                j++;
            }
        }

        if (j > 0)  P[j-1].set_bound(0);
        if (min_k <= 0)  min_k = 1;

        for (k = min_k; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }

    // sequential dynamic greedy coloring and sort
    static void neigh_coloring_dense(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<short>& ind,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size() + 1;

        colors[1].clear();   colors[2].clear();

        for (int w=0; w < P.size(); w++) {
            u = P[w].get_id();
            k = 1, k_prev = 0;

            while (k > k_prev) {
                k_prev = k;
                for (int i = 0; i < colors[k].size(); i++) { //use directly, sort makes it fast!
                    if (adj[u][colors[k][i]]) {
                        k++;
                        break;
                    }
                }
            }

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            //
            if (k < min_k) {
                P[j].set_id(u);
                j++;
            }
        }

        if (j > 0)  P[j-1].set_bound(0);
        if (min_k <= 0)  min_k = 1;

        // start at the min color class
        for (k = min_k; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }



    /**
     * MAXIMUM CLIQUE ENUMERATION
     */

    // sequential dynamic greedy coloring and sort
    static void neigh_coloring_enum(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<short>& ind,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int* pruned,
            int& mc) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size();// + 1;
        colors[1].clear();   colors[2].clear();

        for (int w=0; w < P.size(); w++) {
            u = P[w].get_id();
            k = 1, k_prev = 0;

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

            while (k > k_prev) {
                k_prev = k;
                for (int i = 0; i < colors[k].size(); i++) {
                    if (ind[colors[k][i]]) {
                        k++;
                        break;
                    }
                }
            }

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;

            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            colors[k].push_back(u);
            if (k < min_k) {
                P[j].set_id(u);
                j++;
            }
        }

        if (j > 0)  P[j-1].set_bound(0);
        if (min_k <= 0)  min_k = 1;

        for (k = min_k; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }

    // sequential dynamic greedy coloring and sort
    static void neigh_coloring_dense_enum(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<short>& ind,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            bool** &adj) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        int min_k = mc - C.size();// + 1;

        colors[1].clear();   colors[2].clear();

        for (int w=0; w < P.size(); w++) {
            u = P[w].get_id();
            k = 1, k_prev = 0;

            while (k > k_prev) {
                k_prev = k;
                for (int i = 0; i < colors[k].size(); i++) { //use directly, sort makes it fast!
                    if (adj[u][colors[k][i]]) {
                        k++;
                        break;
                    }
                }
            }

            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            colors[k].push_back(u);
            if (k < min_k) {
                P[j].set_id(u);
                j++;
            }
        }

        if (j > 0)  P[j-1].set_bound(0);
        if (min_k <= 0)  min_k = 1;

        for (k = min_k; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }



    /**
     * SPARSE GRAPHS - ONLY CSC
     */

    /**
     * @brief       Searches a color class for a vertex w that only has a
     *              single conflict with v. The vertex w is exchanged with
     *              another suitable candidate, thereby reducing the number
     *              of colors by one.
     *
     * @param v
     * @param k             color class to search
     * @param colors        two dimensional array where colors[k][i]
     *                      represents the ith vertex in the k color class
     * @param adj           adj graph representation
     * @param w_pos         ptr to w in the k color class
     *
     * @return vertex to exchange
     */
    inline
    static int get_conflict(int v, int k, vector< vector<int> >& colors, int & w_pos,
            vector<long long>& vs, vector<int>& es, vector<short> &ind) {
        int conflictVar = -1;
        int count = 0;

        for (int i = 0; i < colors[k].size(); i++) { // check vertices in the k color class
            int w = colors[k][i];   // get a vertex
            if (ind[w]) {        // edge between w and v, i.e., conflict!
                conflictVar = w;
                w_pos = i;
                count++;            // count how many conflicts v has!
                if (count > 1)  return -1;  // if more than one, ignore
            }
        }
        return conflictVar;
    }



    inline
    static bool uncolored(int v, int k, vector< vector<int> >& colors,
            vector<long long>& vs, vector<int>& es, vector<short> &ind) {

        for (int h = 0; h < colors[k].size(); h++)
            if (ind[colors[k][h]])  return true;
        return false;
    }


    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     */
    static bool repair(int v, int k, vector< vector<int> >& colors,
            vector<long long>& vs, vector<int>& es, vector<short> &ind, vector<short> &hash) {

        //! check each color class
        //! todo: more likely to encounter suitable vertices from reverse order
        for (int color_class = k - 1; color_class >= 0; color_class--) {
            //        for (int color_class = 0; color_class < k - 1; color_class++) {

            /* find a color class where v has only a single conflict!
             * w is the vertex conflicting with u  */
            int w_pos = 0; // holds the position of w in colors[k]
            int w = get_conflict(v,color_class,colors,w_pos,vs,es,ind);

            //            for (long long h = vs[v]; h < vs[v + 1]; h++)  ind[es[h]] = 0;

            if (w >= 0) {

                for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 1;

                //! check each color class
                for (int col_j = color_class + 1; col_j < k; col_j++) {
                    //! find a new color class for w, then swap w and v!
                    if (!uncolored(w,col_j,colors,vs,es,hash)) {
                        colors[k].clear(); // remove v from kth color class
                        colors[color_class][w_pos] = v;
                        colors[col_j].push_back(w);
                        for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 0;
                        return true;
                    }
                }

                for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 0;
            }
        }
        return false;
    }




    /**
     * @brief       Find a color class for the vertex such that
     *              the vertex is not adjacent to any other vertex
     *              in that color class.
     */
    inline
    static bool uncolored(int v, int k, vector< vector<int> >& colors,
            vector<long long>& vs, vector<int>& es, vector<int> &ind) {

        for (int h = 0; h < colors[k].size(); h++)
            if (ind[colors[k][h]])  return true;
        return false;
    }

    inline
    static int get_conflict(int v, int k, vector< vector<int> >& colors, int & w_pos,
            vector<long long>& vs, vector<int>& es, vector<int> &ind) {
        int conflictVar = -1;
        int count = 0;

        for (int i = 0; i < colors[k].size(); i++) { // check vertices in the k color class
            int w = colors[k][i];   // get a vertex
            if (ind[w]) {        // edge between w and v, i.e., conflict!
                conflictVar = w;
                w_pos = i;
                count++;            // count how many conflicts v has!
                if (count > 1)  return -1;  // if more than one, ignore
            }
        }
        return conflictVar;
    }


    /**
     * @brief   If a new color class is started, then we attempt to
     *          fix it, by swapping it with another suitable vertex.
     *          The idea is to find a vertex in another color class
     *          that is adjacent to only one other vertex.
     */
    static bool repair_forward(int v, int k, vector< vector<int> >& colors,
            vector<long long>& vs, vector<int>& es, vector<int> &ind, vector<int> &hash) {

        //! check each color class
        //! todo: more likely to encounter suitable vertices from reverse order
        for (int color_class = 1; color_class < k-1; ++color_class) {
            //        for (int color_class = 0; color_class < k - 1; color_class++) {

            /* find a color class where v has only a single conflict!
             * w is the vertex conflicting with u  */
            int w_pos = 0; // holds the position of w in colors[k]
            int w = get_conflict(v,color_class,colors,w_pos,vs,es,ind);

            //            for (long long h = vs[v]; h < vs[v + 1]; h++)  ind[es[h]] = 0;

            if (w >= 0) {

                for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 1;

                //! check each color class
                for (int col_j = color_class + 1; col_j < k; ++col_j) {
                    //! find a new color class for w, then swap w and v!
                    if (!uncolored(w,col_j,colors,vs,es,hash)) {
                        colors[k].clear(); // remove v from kth color class
                        colors[color_class][w_pos] = v;
                        colors[col_j].push_back(w);
                        for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 0;
                        return true;
                    }
                }

                for (long long h = vs[w]; h < vs[w + 1]; h++)  hash[es[h]] = 0;
            }
        }
        return false;
    }



    static void greedy_neighborhood_coloring_repair(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<Vertex> & ColOrd,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            vector<short> &ind,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        vector<short> hash(vs.size(),0); // hash table for repair function

        // color vertices from 0 to |P|
        for (int w = 0; w < ColOrd.size(); w++) {
            u = ColOrd[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below


            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

            while (uncolored(u,k,colors,vs,es,ind))  k++;

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            // ind is reset in repair!
            //            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;

            // only try to repair if vertex caused a new color class to be added
            if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,vs,es,ind,hash))  max_k--;

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;


        }
        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }




    static void greedy_neighborhood_coloring_repair_only(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            vector<short> &ind,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        vector<short> hash(vs.size(),0); // hash table for repair function

        // color vertices from 0 to |P|
        for (int w = 0; w < P.size(); w++) {
            u = P[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below


            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

            while (uncolored(u,k,colors,vs,es,ind))  k++;

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            // only try to repair if vertex caused a new color class to be added
            if (k+1 > min_k && colors[k].size() == 1 && repair(u,k,colors,vs,es,ind,hash))  max_k--;

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;

        }
        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }


    static void greedy_neighborhood_coloring_repair_only_forward(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& num_colors,
            vector<int> &ind,
            bool is_enum = false,
            bool is_repair = true) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;

//        int min_k = mc - C.size() + 1;
//        if (is_enum) min_k = mc - C.size();


        colors[1].clear();   colors[2].clear();

        vector<int> hash(vs.size(),0); // hash table for repair function


        if (is_repair) {
            // color vertices from 0 to |P|
            for (int w = 0; w < P.size(); w++) {
                u = P[w].get_id(); // use initial order
                k = 1;                  // color class for vertex u, so far...
                k_prev = 0;             // used to check of k changes in the while loop below


                for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

                while (uncolored(u,k,colors,vs,es,ind))  k++;

                // start new color class
                if (k > max_k) {
                    max_k = k;
                    colors[max_k+1].clear();
                }

                // add vertex to color class
                colors[k].push_back(u);

                // only try to repair if vertex caused a new color class to be added
//                if (k+1 > min_k && colors[k].size() == 1 && repair_forward(u,k,colors,vs,es,ind,hash))  max_k--;

                if (colors[k].size() == 1 && repair_forward(u,k,colors,vs,es,ind,hash))  max_k--;

                for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;
            }
        }
        else {
            for (int w = 0; w < P.size(); w++) {
                u = P[w].get_id(); // use initial order
                k = 1;                  // color class for vertex u, so far...
                k_prev = 0;             // used to check of k changes in the while loop below


                for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

                while (uncolored(u,k,colors,vs,es,ind))  k++;

                // start new color class
                if (k > max_k) {
                    max_k = k;
                    colors[max_k+1].clear();
                }

                // add vertex to color class
                colors[k].push_back(u);

                // only try to repair if vertex caused a new color class to be added
//                if (k+1 > min_k && colors[k].size() == 1 && repair_forward(u,k,colors,vs,es,ind,hash))  max_k--;

                for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;
            }
        }

        num_colors = max_k;

//        for (k = 1; k <= max_k; k++)
//            for (int w = 0; w < colors[k].size(); w++) {
//                P[j].set_id(colors[k][w]);
//                P[j].set_bound(k);
//                j++;
//            }
    }




    static void greedy_neighborhood_coloring_static(
            vector<long long>& vs,
            vector<int>& es,
            vector<Vertex> &P,
            vector<Vertex> & ColOrd,
            vector<int>& C,
            vector<int>& C_max,
            vector< vector<int> >& colors,
            int& mc,
            vector<short> &ind,
            bool is_enum = false) {

        int j = 0, u = 0, k = 1, k_prev = 0;
        int max_k = 1;
        //        int min_k = mc - C.size() + 1;
        int min_k = mc - C.size() + 1;
        if (is_enum) min_k = mc - C.size();

        colors[1].clear();   colors[2].clear();

        // color vertices from 0 to |P|
        for (int w = 0; w < ColOrd.size(); w++) {
            u = ColOrd[w].get_id(); // use initial order
            k = 1;                  // color class for vertex u, so far...
            k_prev = 0;             // used to check of k changes in the while loop below

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 1;

            while (uncolored(u,k,colors,vs,es,ind))  k++;

            //            while (k > k_prev) { // stop when we've found a color class
            //                k_prev = k;
            //                /*
            //                 * While a vertex w in C_k (in order) does not have an edge to u,
            //                 * If u does not have an edge to any vertex in C_k, then add u to C_k,
            //                 * Else try the next color class!
            //                 */
            //                for (int i = 0; i < colors[k].size(); i++) {
            //                    if (ind[colors[k][i]]) { // if u has an edge with a vertex in C_k
            //                        k++;                    // then set k=k+1 and break out of for loop!
            //                        break;
            //                    }
            //                }
            //            }

            for (long long h = vs[u]; h < vs[u + 1]; h++)  ind[es[h]] = 0;

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

        }

        for (k = 1; k <= max_k; k++)
            for (int w = 0; w < colors[k].size(); w++) {
                P[j].set_id(colors[k][w]);
                P[j].set_bound(k);
                j++;
            }
    }

}
#endif

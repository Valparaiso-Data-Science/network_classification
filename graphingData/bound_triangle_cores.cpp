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
#include "graphpack_headers.h"
#include <iostream>
#include <fstream>
#include "graphpack_vertex.h"
#include <algorithm>
#include "graphpack_utils.h"

using namespace graphpack;
using namespace std;

/**
 * set block_size via set_block_size() before calling this function
 * @param is_parallel
 */
void graphpack_graph::triangle_cores_in_parallel(bool is_parallel) {

    map_edgeid_to_vertex_pair(is_parallel);

    // number of triangles each edge participates
    if (is_dense)
        triangle_counting_in_parallel_adj(is_parallel);
    else
        triangle_counting_in_parallel(is_parallel);

    printf("[triangle cores]  argmax tr(u) = %lld,  sum_u tr(u) = %llu \n", max_t_edge, total_t);

    // compute the triangle cores
    triangle_core_numbers_parallel();

    // compute bound
//    triangle_max_core();
}

/**
 * @brief   Compute the maximum triangle core
 */
long long graphpack_graph::triangle_max_core() {
    long long m = num_edges() + 1;
    max_tri_core = -1;

    for (long long e=0; e<m; e++) {
        if (e < 20)
            cout << tri_core[e] << " ";
        if (tri_core[e] > max_tri_core)
            max_tri_core = tri_core[e];
    }
    return max_tri_core;
}





void graphpack_graph::map_edgeid_to_vertex_pair(bool is_parallel) {
    // edges (instead of vert ID store EID into e_u and e_v....
    long long i, j, u, v;
    long long m = edges.size();

    eid.resize(m);

    if (is_parallel) {

        #pragma omp parallel for schedule(dynamic,block_size) private(i,j,u,v)
        for (i=0; i < num_edges(); i++) {
            v = e_v[i];
            u = e_u[i];

            for (j=vertices[v]; j<vertices[v+1]; j++) {
                if (edges[j] == u) {
                    eid[j] = i;
                    break;
                }
            }

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                if (edges[j] == v) {
                    eid[j] = i;
                    break;
                }
            }
        }

    }
    else { // serial
        for (i=0; i < num_edges(); i++) {
            v = e_v[i];
            u = e_u[i];

            for (j=vertices[v]; j<vertices[v+1]; j++) {
                if (edges[j] == u) {
                    eid[j] = i;
                    break;
                }
            }

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                if (edges[j] == v) {
                    eid[j] = i;
                    break;
                }
            }
        }
    }
}








// FIXED OFFSET
void graphpack_graph::triangle_counting_in_parallel(bool is_parallel) {
    long long e, v, u, w, i, j;
    long long n = num_vertices();
    long long m = num_edges(); // this is |E|/2

//    vector<long long> triangles(m,0);
//    tri_core.resize(m,0);
    vector<long long> triangles_tmp; // = &tri_core;
    triangles_tmp.resize(m,0);
    vector<long long> ind(n, 0);

    if (is_parallel) {
        #pragma omp parallel for schedule(dynamic,block_size) \
        shared(triangles_tmp) firstprivate(ind) private(v,u,w,e,i,j)
        for (e=0; e<m; e++) {
            v = e_v[e];
            u = e_u[e];
            //        if (u==v) { continue; }

            for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = e+1;
//            ind[v] = 1;

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                w = edges[j];
//                if (w==u) { continue; }  //self-loop
//                if (w==v) { continue; }  //discard 1-hop self-loop
                if (ind[w] == e+1) {
                    triangles_tmp[e] += 1;
                }
            }

//            for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
//            ind[v] = 0;
        }
    }
    else {
        for (e=0; e<m; e++) {
            v = e_v[e];
            u = e_u[e];
            //        if (u==v) { continue; }
            for (i=vertices[v]; i<vertices[v+1]; i++)
                ind[edges[i]] = e+1;
//            ind[v] = 1;

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                w = edges[j];
//                if (w==u) { continue; }  //self-loop
//                if (w==v) { continue; }  //discard 1-hop self-loop
                if (ind[w] == e+1) {
                    triangles_tmp[e] += 1;
                }
            }

//            for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
            ind[v] = 0;
        }
    }

    max_t_edge = 0;
    total_t = 0;
    for (e=0; e<m; e++) {
        total_t += triangles_tmp[e];
        if (triangles_tmp[e] > max_t_edge) {
            max_t_edge = triangles_tmp[e];
        }
    }

    tri_core = triangles_tmp;

//    tmp_tri_core.clear();
    ind.clear();
}


/*
 * @brief   ADJ data structure basically removes the cost
 *          of creating the fast index structure!
 *
 */
void graphpack_graph::triangle_counting_in_parallel_adj(bool is_parallel) {
    long long e, v, u, w, i, j;
    long long n = num_vertices();
    long long m = num_edges(); // this is |E|/2

//    vector<long long> triangles(m,0);
//    triangles.resize(m,0);
    vector<long long> triangles_tmp;// = &tri_core;
//    tri_core.resize(m,0);
    triangles_tmp.resize(m,0);

    if (is_parallel) {
        #pragma omp parallel for schedule(dynamic,block_size) \
        shared(triangles_tmp) private(v,u,w,e,i,j)
        for (e=0; e<m; ++e) { // v--u
            v = e_v[e];
            u = e_u[e];
            //        if (u==v) { continue; }

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
                w = edges[j];
//                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (adj[v][w]) {
                    triangles_tmp[e] += 1;
                }
            }
        }
    }
    else {
        for (e=0; e<m; ++e) { // v--u
            v = e_v[e];   u = e_u[e];
            //        if (u==v) { continue; }

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
                w = edges[j];
//                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (adj[v][w]) {
                    triangles_tmp[e] += 1;
                }
            }
        }
    }

    max_t_edge = 0;
    total_t = 0;
    for (e=0; e<m; e++) {
        total_t += triangles_tmp[e];
        if (triangles_tmp[e] > max_t_edge) {
            max_t_edge = triangles_tmp[e];
        }
    }

    tri_core = triangles_tmp;
//    tmp_tri_core.clear();
//    ind.clear();
}




void graphpack_graph::triangle_core_numbers_parallel() {
    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges();// + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed


    for (e=0; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;


    /** compute starting positions for the edges in each bin */
//    start = 1;
    start = 0;
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    /** bucket sort the edges by their triangle counts
      *  place the edges into the correct buckets */
    for (e=0; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    /** fixup the bin for the main computation */
    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 0;
//    bin[0] = 1;

#ifndef DEBUG
    print_line(80);
    printf("printing the bin array -- edges are binned by their triangle counts \n");
    for (e=0; e<m; e++) {
        printf("e_idx = %lld, T[e_idx] = %lld (num triangles), pos[e_idx] = %lld, order[pos[e_idx]] = %lld \n",
                e, tri_core[e], pos[e], tris[pos[e]]);
//        pos[e] = bin[tri_core[e]];
//        tris[pos[e]] = e;
//        bin[tri_core[e]]++;
    }
    print_line(80);

    // print eid and e_v and e_u arrays!

    printf("printing the bin array -- edges are binned by their triangle counts \n");
    for (tri=max_t_edge; tri > 0; tri--) {
        printf("%lld ", bin[tri]);
    }
    printf("\n");
#endif


//    if (m < 20) {
//        for (e=0; e<m; e++) {
//            eid_vu = tris[e]; // @TODO: Could save the ordering..
//
//            /** get the vertex ids of each end of the edge */
//            v = e_v[eid_vu];
//            u = e_u[eid_vu];
//
//            eid[]
//        }
//    }

    /**
     * Process each edge in tris in the order of number of triangles
     */
    double sec = get_time();
    vector<long long> ind(n, 0);
    for (e=0; e<m; e++) {


        /** get the edge_id */
        eid_vu = tris[e]; // @TODO: Could save the ordering..

        /** get the vertex ids of each end of the edge */
        v = e_v[eid_vu];
        u = e_u[eid_vu];

        if (e % 100000 == 0) printf("(%lld/%lld),  T(%lld,%lld)=%lld,  %lg sec\n", e,m,v,u,tri_core[eid_vu], get_time()-sec);



        /**
         * Create lookup table for the neighbors of v.
         * The storage cost is O(|E|).
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, and therefore this indicates the exact neighbor.
         * --------------
         * For each neighbor of v denoted u \in N(v), we mark ind(u),
         * but instead of simply marking it, we store the ptr to that
         * position in the edges/neighbors array.
         *
         * This provides us with a way to check for a triangle in O(1) time, and also
         * if a triangle does exist from that vertex, it allows us to grab its edgeid in O(1) time.
         *
         * The only problem is if an edge has a ptr of 0, then when we check for a neighbor, it will appear as the neighbor does not exist, while its actual position is 0.
         * For this case, we simply add 1 to edge pos when marking that vertex, and after if we need it, then we subtract 1 from it to get the real edge position.
         *
         *
         * We store the edge ptr of each vertex, and then use
         * This works
         */
        for (i=vertices[v]; i<vertices[v+1]; i++) {
//            if (proc[eid[i]]) continue;
                ind[edges[i]] = i+1;
        }
//        ind[v] = v+1;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */

        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
            if (ind[edges[j]] && !proc[eid[j]] && !proc[eid[ind[edges[j]]-1]]) {//continue;
//            w = edges[j];
//            if (w==v) { continue; }  //self-loop
//            if (w==u) { continue; }  //discard 1-hop self-loop

            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
//            if (ind[w]) {
//                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
                edge_neigh[0] = eid[j];
                /**
                 * We subtract one here, since we added one above when marking v's neighbors.
                 * This was done precisely for the reason of getting the egdeid in O(1) time.
                 * Since the edge_ids begin at 0, if we had not added 1,
                 * then the vertex stored in the 0 position (edges[0]) would seem to not be connected to v.
                 * Hence, we subtract 1, again to retrieve the correct edge position, while avoiding the above case.
                 */
//                eid_vw = eid[ind[w]-1];   // get the edge_id of v--w
//                edge_neigh[1] = eid[ind[w]-1];
                edge_neigh[1] = eid[ind[edges[j]]-1];

                /**
                 * We should be able to perform these checks as soon as know the edge_id, then we can check in proc[], and if marked, then we just skip.
                 * This can also be done, when setting the ind, if it has already been processed, then skip it!
                 */

//                if (!proc[edge_neigh[0]] && !proc[edge_neigh[1]]) {
//                    edge_neigh[0] = eid_uw;
//                    edge_neigh[1] = eid_vw;

                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {

                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris[px];                 //starting edge id

                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris[pw] = e_x;
                                pos[e_x] = pw;
                                tris[px] = edge_neigh[k];
                            }
                            bin[tw]++;   // update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
//            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
            ind[edges[i]] = 0;
//        ind[v] = 0;

        proc[eid_vu] = 1; // set edge to be processed
    }

#ifndef DEBUG
        cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
#endif

        max_tri_core = tri_core[tris[m-2]];

//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos

        printf("max triangle core = %lld \n", max_tri_core+2);
//        if (validate) {
//            validate_triangle_cores(bin,tris,pos,ind);
//        }

        tris_ordering = tris;

        if (m < 20) {
            print_line(80);
            printf("[graphpack: triangle core debug info]  \n");
            print_line(80);
            for (long long e=0; e<m; e++) {
//                if (e < 20)
                eid_vu = tris[e]; // @TODO: Could save the ordering..

                /** get the vertex ids of each end of the edge */
                v = e_v[eid_vu];
                u = e_u[eid_vu];

                cout << "e = " << e << ", eid_vu: " << eid_vu << ", v = " << v << ", u = "
                        << u << ", T(" << v << "," << u << ") = " << tri_core[e] << " \n";
//                if (tri_core[e] > max_tri_core)
//                    max_tri_core = tri_core[e];
            }
            print_line(80);
        }

        bin.clear();
        pos.clear();
//        tris.clear();
        ind.clear();
}


/**
 * @brief           sparse parallel edge-centric triangle counting
 *                  these codes have been extremely optimized for
 *                  use within branch-and-bound algorithms. In this case,
 *                  mc can be simply treated as the bound.
 *
 * @param mc
 * @param pruned_edges
 */
void graphpack_graph::k_triangle_counting(int & max, vector<int> & pruned_edges, bool is_parallel) {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);

    if (is_dense_graph()) {
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,block_size) \
            shared(tmp_tri_core) private(v,u,w,e,i,j)
            for (e=0; e<m; ++e) { // v--u
                v = e_v[e];
                u = e_u[e];
                //        if (u==v) { continue; }

                for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
                    w = edges[j];
                    if (w==u) { continue; }  //self-loop
                    if (w==v) { continue; }  //discard 1-hop self-loop
                    if (adj[v][w]) {
                        tmp_tri_core[e] += 1;
                    }
                }
//                if (tmp_tri_core[e]+2 < max) {
//                    tmp_tri_core[e] = 0;
//                    pruned_edges[e] = 1;
//                }
            }
        }
        else {
            for (e=0; e<m; ++e) { // v--u
                v = e_v[e];
                u = e_u[e];
                //        if (u==v) { continue; }

                for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
                    w = edges[j];
                    if (w==u) { continue; }  //self-loop
                    if (w==v) { continue; }  //discard 1-hop self-loop
                    if (adj[v][w]) {
                        tmp_tri_core[e] += 1;
                    }
                }
//                if (tmp_tri_core[e]+2 < max) {
//                    tmp_tri_core[e] = 0;
//                    pruned_edges[e] = 1;
//                }
            }
        }
    }
    else {

        if (is_parallel) {
#pragma omp parallel for schedule(dynamic,block_size) \
        shared(tmp_tri_core) private(v,u,w,e,i,j)
            for (e=0; e<m; e++) { // v--u
                v = e_v[e];   u = e_u[e];
                if (u==v) { continue; }

                for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = i; // store the ptr location in the edge array
                ind[v] = v;

                for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                    /**
                     * @brief       if at any point, there are not enough neighbors to
                     *              to satisfy the bound, then we terminate the counting
                     *              early, rather than doing unnecessary work. Note that
                     *              (vertices[u+1] - j) is an upper-bound on the number
                     *              of possible triangles remaining for that vertex
                     *
                     *              Advantages
                     *              ----------
                     *              * Works well for graphs where the max clique is large, then
                     *              we have many early terminations, which can significantly
                     *              increase the performance (and cost of using such a tight bound)
                     */
                    //                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    //                    tmp_tri_core[e] = 0;
                    //                    pruned_edges[e] = 1;
                    //
                    //                    break;        // avoid computing the rest
                    //                }

                    w = edges[j];
                    if (w==u) { continue; }  //self-loop
                    if (w==v) { continue; }  //discard 1-hop self-loop
                    if (ind[w]) {
                        tmp_tri_core[e] += 1;
                    }
                }


//                if (tmp_tri_core[e]+2 < max) {
//                    tmp_tri_core[e] = 0;
//                    pruned_edges[e] = 1;
//
//                    // removing one edge, also impacts the triangle counts of other edges!
//                    // however, since we force an ordering, then we can check very quickly up to that ordering
//                    // hence, we only must check the previous edge triangle counts, since any future triangle counts would first check the pruned edges
//                    //            update_neighbors(tmp_tri_core);
//                }

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
                ind[v] = 0;
            }
        }
        else {
            for (e=0; e<m; e++) { // v--u
                v = e_v[e];   u = e_u[e];
                if (u==v) { continue; }

                //        for (i=vertices[v]; i<vertices[v+1]; i++)
                //            ind[edges[i]] = 1;
                //        ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = i; // store the ptr location in the edge array
                ind[v] = v;

                for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                    /**
                     * @brief       if at any point, there are not enough neighbors to
                     *              to satisfy the bound, then we terminate the counting
                     *              early, rather than doing unnecessary work. Note that
                     *              (vertices[u+1] - j) is an upper-bound on the number
                     *              of possible triangles remaining for that vertex
                     *
                     *              Advantages
                     *              ----------
                     *              * Works well for graphs where the max clique is large, then
                     *              we have many early terminations, which can significantly
                     *              increase the performance (and cost of using such a tight bound)
                     */
                    //                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    //                    tmp_tri_core[e] = 0;
                    //                    pruned_edges[e] = 1;
                    //
                    //                    break;        // avoid computing the rest
                    //                }

                    w = edges[j];
                    if (w==u) { continue; }  //self-loop
                    if (w==v) { continue; }  //discard 1-hop self-loop
                    if (ind[w]) {
                        tmp_tri_core[e] += 1;
                    }
                }


//                if (tmp_tri_core[e]+2 < max) {
//                    tmp_tri_core[e] = 0;
//                    pruned_edges[e] = 1;
//
//                    // removing one edge, also impacts the triangle counts of other edges!
//                    // however, since we force an ordering, then we can check very quickly up to that ordering
//                    // hence, we only must check the previous edge triangle counts, since any future triangle counts would first check the pruned edges
//                    //            update_neighbors(tmp_tri_core);
//                }

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
                ind[v] = 0;
            }
        }
    }


//    if (is_gstats) { // don't compute unless specified
        max_t_edge = 0;
        total_t = 0;
        for (e=0; e<m; e++) {
            total_t += tmp_tri_core[e];
            if (tmp_tri_core[e] > max_t_edge) {
                max_t_edge = tmp_tri_core[e];
            }
        }
//    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}















void graphpack_graph::validate_triangle_cores(vector<long long> bin, vector<long long> tris, vector<long long> pos, int mt_end) {
    int m = num_edges();
    for (uint64_t e = 0; e < m; e++) {

    }

}

//void graphpack_graph::triangle_core_numbers_parallel() {
//    triangle_core_numbers_parallel(false);
//}


//
////void triangle_cores_corrected();
//
//
//void graphpack_graph::triangle_cores_corrected() {
//    long long e, e_ux, e_x, tri, i, j, k;
//    long long start, num, tmp;
//    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
//    long long n = num_vertices(), m = num_edges() + 1;
//    long long mt_end = max_t_edge+1;
//
//    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
//    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
//    vector <long long> pos(m);        // position of the edges in tris
//
//    vector<long long> edge_neigh(2,0);  //exactly two neighboring edges of e form a triangle
//    vector<short> proc(m,0); //a triangle is processed, if any one of its edges are processed
//
//
//    for (e=1; e<m; e++) // bin edges by triangle counts
//        bin[tri_core[e]]++;
//
//    start = 1;  // bucket sort
//    for (tri=0; tri < mt_end; tri++) {
//        num = bin[tri];
//        bin[tri] = start;
//        start = start + num;
//    }
//
//    for (e=1; e<m; e++) {
//        pos[e] = bin[tri_core[e]];
//        tris[pos[e]] = e;
//        bin[tri_core[e]]++;
//    }
//
//    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
//    bin[0] = 1;
//
//
//    vector<long long> ind(n, 0);
//
//    // only allows two edges
//    for (e = 1; e < m; e++) { //for each edge in tris in the order of number of triangles
//        v = e_v[tris[e]-1];
//        u = e_u[tris[e]-1];
//
//        for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = i;
////      ind[v] = v;
//
//        for (j=vertices[u]; j<vertices[u+1]; j++) { // for neigh of u
//
//            if (ind[edges[j]]) { // if neigh of v
//                if (edges[j]==v) { continue; }  //self-loop
//                if (edges[j]==u) { continue; }  //discard 1-hop self-loop
//
//                // check if neighboring edges uw and vw have been marked
//                if (proc[eid[j]] == 0 && proc[eid[ind[edges[j]]]] == 0) {
//                    edge_neigh[0] = eid[j];
//                    edge_neigh[1] = eid[ind[edges[j]]];
//                    for (k=0; k<2; k++) { //each edge neighbors with two other edges to form a triangle
//                        if (tri_core[edge_neigh[k]] > tri_core[tris[e]]) {
//                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
//                            pw = pos[edge_neigh[k]];        //actual pos of w
//                            px = bin[tw];                   //starting position of tw
//                            e_x = tris[px];                 //starting edge id
//                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
//                                pos[edge_neigh[k]] = px;   //starting position of e_
//                                tris[pw] = e_x;
//                                pos[e_x] = pw;
//                                tris[px] = edge_neigh[k];
//                            }
//                            bin[tw]++;   //update starting pos
//                            tri_core[edge_neigh[k]]--;
//                        }
//                    }
//                }
//            }
//        }
//
//        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
//        ind[v] = 0;
//
//        proc[tris[e]] = 1; //set edge to be processed
//    }
//
////  cout << "median triangle core = " << tri_core[tris[ceil(tris.size()-2/2)]] <<endl;
////  cout << "min triangle core = " << tri_core[tris[1]] <<endl;
//    cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
//    max_tri_core = tri_core[tris[m-1]];
//
////  string::size_type result;
////  string::size_type beg;
////  string name = "";
////  result = fn.rfind('.', fn.size() - 1);
////  beg = fn.rfind('../data/output/', fn.size() - 1);
////  if(result != string::npos)
////      name = fn.substr(beg+1,result);
//
////  ofstream myfile;
////  fn.append("-tcore.txt");
////  char *fileName = (char*)fn.c_str();
////  myfile.open (fileName);
//
//
//    int sum = 0;
//    for (e=0; e<m-1; e++) {
////      cout << "e: " << e << ", id: " << tris[e+1] << ", T = " << tri_core[tris[e+1]] <<endl;
////      myfile << tri_core[tris[e]] <<endl;
////      tri_core[e] = tri_core[e+1];
////      tris[e] = tris[e];
////      cout << tri_core[tris[e]] <<endl;
//        sum += tri_core[e];
//
//    }
////  myfile.close();
//
////  cout << "avg triangle core = " << sum/m <<endl;
//
//    triangle_max_core();
//
//    if (m < 20) {
//        print_line(80);
//        printf("[graphpack: triangle core debug info]  \n");
//        print_line(80);
//        for (long long e=0; e<m; e++) {
////                if (e < 20)
//            eid_vu = tris[e]; // @TODO: Could save the ordering..
//
//            /** get the vertex ids of each end of the edge */
//            v = e_v[eid_vu];
//            u = e_u[eid_vu];
//
//            cout << "e = " << e << ", eid_vu: " << eid_vu << ", v = " << v << ", u = "
//                    << u << ", T(" << v << "," << u << ") = " << tri_core[e] << " \n";
////                if (tri_core[e] > max_tri_core)
////                    max_tri_core = tri_core[e];
//        }
//        print_line(80);
//    }
//
//    bin.clear();
//    pos.clear();
//    tris.clear();
//    ind.clear();
//}






/**
 * Assign a k-core number to each edge
 *    K_e(u,v) =  argmax K(u),K(v)
 */
void graphpack_graph::parallel_edge_kcores() {
//    long long e, v, u, w, i, j;
//    long long n = num_vertices(), m = num_edges();
//
//    long long tmp_max_t_edge = 0;
//    vector<long long> tri_core(m+1,0);
//    vector<long long> ind(n, 0);
//    int block_size = 64;
//
//    #pragma omp parallel for schedule(dynamic,block_size) \
//        shared(tri_core)
//    //firstprivate(ind) private(v,u,w,e,i,j)
//    for (long long e = 0; e < m; e++) {
//        int v = e_v[e];
//        int u = e_u[e];
//
//        if (kcore[v] >= kcore[u]) tri_core[e+1] = kcore[v]
//        else  tri_core[e+1] = kcore[u];
//    }
}







/*
 * @brief       Triangle Core Decomposition
 *              Runtime = O(m^1.5), majority of work is for
 *              computing the triangle counts
 *
 * Author: Ryan A. Rossi
 *
 */

void graphpack_graph::compute_triangle_cores() {
	//	map_edges();
	map_edges_parallel();

	double sec = get_time();
    // number of triangles each edge participates
    if (is_dense)
        parallel_triangle_counting_adj();
    else
        parallel_triangle_counting();
	cout << "triangles: " << total_t << " (time: " << get_time() - sec << "), " <<endl;

	// compute the triangle cores
	sec = get_time();
	parallel_triangle_cores();
	cout << "max triangle core: " << max_tri_core+2 << " (time: " << get_time() - sec << "), " <<endl;
}

/*
 * @brief Parallelized Triangle Core Decomposition
 * Runtime = O(m^1.5), majority of work is for computing the triangle counts
 *
 * Author: Ryan A. Rossi
 *
 */
void graphpack_graph::compute_parallel_triangle_cores() {
	map_edges_parallel();

	// number of triangles each edge participates
	if (is_dense)
	    parallel_triangle_counting_adj();
	else
	    parallel_triangle_counting();

	// compute the triangle cores
	parallel_triangle_cores();

	// compute bound
	triangle_max_core();
}




void graphpack_graph::test_total_triangles() {
	if (e_u.size() == 0)
		cout << "must use -s option" <<endl;
	else {
		double sec = get_time();
		map_edges_parallel();
		parallel_count_triangles();
		cout << "parallel counting triangles took " << get_time() - sec << " sec, ";
		cout << "total triangles = " << total_t << ", max edge triangles = " << max_t_edge <<endl;
		parallel_triangle_cores();
		triangle_max_core();
		cout << "parallel triangle cores took " << get_time() - sec << " sec, ";
		cout << "max triangle core = " << max_tri_core+2 <<endl;

		double sec2 = get_time();
		map_edges_serial(); //TODO -1 off
		count_triangles();
		cout << "serial counting triangles took " << get_time() - sec2 << " sec, ";
		cout << "total triangles = " << total_t << ", max edge triangles = " << max_t_edge <<endl;
		triangle_cores();
		triangle_max_core();
		cout << "serial triangle cores took " << get_time() - sec2 << " sec, ";
		cout << "max triangle core = " << max_tri_core+2 <<endl;
	}
}




void graphpack_graph::parallel_count_triangles() {
	long long e, v, u, w, i, j;
	long long n = num_vertices(), m = num_edges();


	long long tmp_max_t_edge = 0;
	vector<long long> tmp_tri_core(m+1,0);
	vector<long long> ind(n, 0);

#pragma omp parallel for schedule(dynamic) \
		shared(tmp_tri_core, tmp_max_t_edge) firstprivate(ind) \
private(v,u,w,e,i,j)
	for (e=0; e<m; e++) {
		v = e_v[e];   u = e_u[e];
		if (u==v) { continue; }
		for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 1;
		ind[v] = 1;

		for (j=vertices[u]; j<vertices[u+1]; j++) {
			w = edges[j];
			if (w==u) { continue; }  //self-loop
			if (w==v) { continue; }  //discard 1-hop self-loop
			if (ind[w]) {
				tmp_tri_core[e+1] += 1;
			}
		}

		for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
		ind[v] = 0;
	}

	max_t_edge = 0;
	total_t = 0;
	for (e=0; e<m; e++) {
		total_t += tmp_tri_core[e+1];
		if (tmp_tri_core[e+1] > max_t_edge) {
			max_t_edge = tmp_tri_core[e+1];
		}
	}

	tri_core = tmp_tri_core;

	tmp_tri_core.clear();
	ind.clear();
}







void graphpack_graph::parallel_triangle_cores() {
	long long e, e_ux, e_x, tri, i, j, k;
	long long start, num, tmp;
	long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
	long long n = num_vertices(), m = num_edges() + 1;
	long long mt_end = max_t_edge+1;

	vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
	vector <long long> tris(m);    	  // sort edges in increasing order of triangles: t[i] < t[j]
	vector <long long> pos(m);        // position of the edges in tris

	vector<long long> edge_neigh(2,0);  //exactly two neighboring edges of e form a triangle
	vector<short> proc(m,0); //a triangle is processed, if any one of its edges are processed


	for (e=1; e<m; e++) // bin edges by triangle counts
		bin[tri_core[e]]++;

	start = 1;  // bucket sort
	for (tri=0; tri < mt_end; tri++) {
		num = bin[tri];
		bin[tri] = start;
		start = start + num;
	}

	for (e=1; e<m; e++) {
		pos[e] = bin[tri_core[e]];
		tris[pos[e]] = e;
		bin[tri_core[e]]++;
	}

	for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
	bin[0] = 1;


	vector<long long> ind(n, 0);

	// only allows two edges
	for (e = 1; e < m; e++) { //for each edge in tris in the order of number of triangles
	    v = e_v[tris[e]-1];
	    u = e_u[tris[e]-1];

	    for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = i;
//	    ind[v] = v;

	    for (j=vertices[u]; j<vertices[u+1]; j++) { // for neigh of u

	        if (ind[edges[j]]) { // if neigh of v
	            if (edges[j]==v) { continue; }  //self-loop
	            if (edges[j]==u) { continue; }  //discard 1-hop self-loop

	            // check if neighboring edges uw and vw have been marked
	            if (proc[eid[j]] == 0 && proc[eid[ind[edges[j]]]] == 0) {
	                edge_neigh[0] = eid[j];
	                edge_neigh[1] = eid[ind[edges[j]]];
	                for (k=0; k<2; k++) { //each edge neighbors with two other edges to form a triangle
	                    if (tri_core[edge_neigh[k]] > tri_core[tris[e]]) {
	                        tw = tri_core[edge_neigh[k]];   //number of current tri for w
	                        pw = pos[edge_neigh[k]]; 		//actual pos of w
	                        px = bin[tw]; 					//starting position of tw
	                        e_x = tris[px];  				//starting edge id
	                        if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
	                            pos[edge_neigh[k]] = px;   //starting position of e_
	                            tris[pw] = e_x;
	                            pos[e_x] = pw;
	                            tris[px] = edge_neigh[k];
	                        }
	                        bin[tw]++;   //update starting pos
	                        tri_core[edge_neigh[k]]--;
	                    }
	                }
	            }
	        }
	    }

	    for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
	    ind[v] = 0;

	    proc[tris[e]] = 1; //set edge to be processed
	}

//	cout << "median triangle core = " << tri_core[tris[ceil(tris.size()-2/2)]] <<endl;
//	cout << "min triangle core = " << tri_core[tris[1]] <<endl;
	cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
	max_tri_core = tri_core[tris[m-1]];

//	string::size_type result;
//	string::size_type beg;
//	string name = "";
//	result = fn.rfind('.', fn.size() - 1);
//	beg = fn.rfind('../data/output/', fn.size() - 1);
//	if(result != string::npos)
//		name = fn.substr(beg+1,result);

//	ofstream myfile;
//	fn.append("-tcore.txt");
//	char *fileName = (char*)fn.c_str();
//	myfile.open (fileName);


	int sum = 0;
	for (e=0; e<m-1; e++) {
//		cout << "e: " << e << ", id: " << tris[e+1] << ", T = " << tri_core[tris[e+1]] <<endl;
//		myfile << tri_core[tris[e]] <<endl;
//		tri_core[e] = tri_core[e+1];
//		tris[e] = tris[e];
//		cout << tri_core[tris[e]] <<endl;
		sum += tri_core[e];

	}
//	myfile.close();

//	cout << "avg triangle core = " << sum/m <<endl;

	triangle_max_core();

    if (m < 20) {
        print_line(80);
        printf("[graphpack: triangle core debug info]  \n");
        print_line(80);
        for (long long e=0; e<m; e++) {
//                if (e < 20)
            eid_vu = tris[e]; // @TODO: Could save the ordering..

            /** get the vertex ids of each end of the edge */
            v = e_v[eid_vu];
            u = e_u[eid_vu];

            cout << "e = " << e << ", eid_vu: " << eid_vu << ", v = " << v << ", u = "
                    << u << ", T(" << v << "," << u << ") = " << tri_core[e] << " \n";
//                if (tri_core[e] > max_tri_core)
//                    max_tri_core = tri_core[e];
        }
        print_line(80);
    }

	bin.clear();
	pos.clear();
	tris.clear();
	ind.clear();
}

// fast compute of ONLY max triangle core
void graphpack_graph::compute_max_triangle_core(int lb) {

	//	vector<int> es = get_edges_array();
	//	vector<long long> vs = get_vertices_array();

	int* pruned = new int[num_vertices()];
	memset(pruned, 0, num_vertices() * sizeof(int));

	int lb_idx = 0;
	for (int i = num_vertices()-1; i >= 0; i--) {
		if (kcore[kcore_order[i]] == lb)  lb_idx = i;
		if (kcore[kcore_order[i]] <= lb)  pruned[kcore_order[i]] = 1;
	}


	cout << "[initial pruning: reduce graph]  before pruning: |V| = " << num_vertices() << ", |E| = " << num_edges() <<endl;
	reduce_graph(pruned);
	update_degrees(true);
	cout << "[initial pruning: reduce graph]  after pruning:  |V| = " << num_vertices() - lb_idx << ", |E| = " << num_edges() <<endl;

	double sec = get_time();
	compute_parallel_triangle_cores();

	cout << "[triangle cores]  computed in " << get_time() - sec << " sec, ";
	cout << "max triangle core = " << max_tri_core+2 <<endl;

	delete [] pruned;

}




/*
 * @brief   Serial Triangle Core Decomposition
 */

void graphpack_graph::count_triangles() {
	long long e, v, u, w, i, j;
	long long n = num_vertices(), m = num_edges();

	total_t = 0;
	max_t_edge = -1;
	tri_core.clear();
	tri_core.resize(m+1,0);
	vector<long long> ind(n, 0);

	for (e=0; e<m; e++) {
		v = e_v[e];   u = e_u[e];
		if (u==v) { continue; }
		for (i=vertices[v]; i<vertices[v+1]; i++)
			ind[edges[i]] = 1;
		ind[v] = 1;

		for (j=vertices[u]; j<vertices[u+1]; j++) {
			w = edges[j];
			if (w==u) { continue; }  //self-loop
			if (w==v) { continue; }  //discard 1-hop self-loop
			if (ind[w]) tri_core[e+1] += 1;
		}
		if (tri_core[e+1] > max_t_edge)  max_t_edge = tri_core[e+1];

		for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
		ind[v] = 0;

		//		total_t += tri_core[e+1];

#ifdef _DEBUG
		cout << "(v=" << v << ", u=" << u << ") | e=" << e+1 << " deg_v=" << degree[v] << " deg_u=" << degree[u] << " t_edge=" << t_edge[e] <<endl;
#endif
	}

	//	max_t_edge = 0;
	total_t = 0;
	for (e=0; e<m; e++) {
		total_t += tri_core[e+1];
		//		if (tmp_tri_core[e+1] > max_t_edge) {
		//			max_t_edge = tmp_tri_core[e+1];
		//		}
	}

	ind.clear();
}




void graphpack_graph::triangle_cores() {
	long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
	long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
	long long n = num_vertices(), m = num_edges() + 1;
	long long mt_end = max_t_edge+1;

	vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
	vector <long long> tris(m);    	  // sort edges in increasing order of triangles: t[i] < t[j]
	vector <long long> pos(m);        // position of the edges in tris

	vector<long long> edge_neigh(2);  //exactly two neighboring edges of e form a triangle
	vector<long long> proc(m,0); //a triangle is processed, if any one of its edges are processed


	for (e=1; e<m; e++) // bin edges by triangle counts
		bin[tri_core[e]]++;

	start = 1;  // bucket sort
	for (tri=0; tri < mt_end; tri++) {
		num = bin[tri];
		bin[tri] = start;
		start = start + num;
	}

	for (e=1; e<m; e++) {
		pos[e] = bin[tri_core[e]];
		tris[pos[e]] = e;
		bin[tri_core[e]]++;
	}

	for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
	bin[0] = 1;

#ifdef _DEBUG
	print_tcore_info(bin,tris,pos,mt_end);
#endif

	vector<long long> ind(n, 0);
	// note that m here is actually num_edges() + 1.
	for (e=1; e<m; e++) { //for each edge in tris in the order of number of triangles
		eid_vu = tris[e];
		v = e_v[eid_vu-1];
		u = e_u[eid_vu-1];

		//		sort tris
		//		if (degree[u] < degree[v]) v = u;

		for (i=vertices[v]; i<vertices[v+1]; i++)
			ind[edges[i]] = i; //i+1; //idx = neigh vertex, value = edge_id+1
		ind[v] = v;

#ifdef _DEBUG
		cout << "processing edge: " << eid_vu;
#endif

		for (j=vertices[u]; j<vertices[u+1]; j++) { // for neigh of u
			w = edges[j];
			if (w==v) { continue; }  //self-loop
			if (w==u) { continue; }  //discard 1-hop self-loop

			if (ind[w]) { // if neigh of v
				eid_uw = eid[j];// + 1;
				eid_vw = eid[ind[w]];// + 1; //-1] + 1; //get edge ids

#ifdef _DEBUG
				cout << "  vu-uw-vw: " << eid_vu << "-" << eid_uw << "-" << eid_vw;
				cout << "    Tvu: " << tri_core[eid_vu] << ", Tuw: " << tri_core[eid_uw] << ", Tvw: " << tri_core[eid_vw];
				cout << "\t eid_vu: " << (eid_vu-1) << ", eid_vw j: " << j << ", eid_uw ind[w]-1: " << ind[w]-1 <<endl;
#endif

				if (!proc[eid_uw] && !proc[eid_vw]) {
					edge_neigh[0] = eid_uw;
					edge_neigh[1] = eid_vw;
					for (k=0; k<2; k++) { //each edge neighbors with two other edges to form a triangle
						if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {
#ifdef _DEBUG
							print_tcore_info(bin,tris,pos,mt_end);
#endif
							tw = tri_core[edge_neigh[k]];   //number of current tri for w
							pw = pos[edge_neigh[k]]; 		//actual pos of w
							px = bin[tw]; 					//starting position of tw
							e_x = tris[px];  				//starting edge id
#ifdef _DEBUG
							cout << "\ttw: " << tw << ", pw: " << pw << ", px: " << px << ", e_x: " << e_x <<endl;
#endif
							if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
								pos[edge_neigh[k]] = px;   //starting position of e_
								tris[pw] = e_x;
								pos[e_x] = pw;
								tris[px] = edge_neigh[k];
							}
							bin[tw]++;   //update starting pos
							tri_core[edge_neigh[k]]--;
#ifdef _DEBUG
							print_tcore_info(bin,tris,pos,mt_end);
#endif
						}
					}
				}
			}
		}

		for (i=vertices[v]; i<vertices[v+1]; i++) //reset index struct
			ind[edges[i]] = 0;
		ind[v] = 0;

		proc[eid_vu] = 1; //set edge to be processed


#ifdef _DEBUG
		cout << "\tproc = ";
		for (i=1; i < m; i++) {
		    cout << proc[i] << " ";
		}
		cout <<endl;
#endif

	}

	//	int sum = 0;
	//finished Triangle Core Decomposition
//	for (e=0; e<m-1; e++) {
//		tri_core[e] = tri_core[e+1];
		//		sum += tri_core[e];
//#ifdef _DEBUG
//		cout << "e=" << e+1 << ", T = " << tri_core[e] <<endl;
//#endif
//	}

	//	cout << "avg triangle core = " << sum/m <<endl;
	//	cout << "median triangle core = " << tris[ceil(tris.size()-2/2)] <<endl;
	//	cout << "min triangle core = " << tris[1] <<endl;
	//	cout << "max triangle core = " << tris[tris.size()-3] <<endl;
	triangle_max_core();

	//	cout << "max triangle core: " << max_tri_core << ", ub is max_tri_core+2: " << max_tri_core+2 <<endl;

	// clear the buffers
	bin.clear();
	pos.clear();
	tris.clear();
	ind.clear();
}




void graphpack_graph::print_tcore_info(vector<long long> bin, vector<long long> tris, vector<long long> pos, int mt_end) {
	long long e, tri, m = num_edges()+1;

	cout << "\n\ttdeg = ";
	for (e=1; e<m; e++)
		cout << tri_core[e] << "\t";

	cout << "\n\tbins = ";
	for (tri=0; tri < mt_end; tri++)
		cout << bin[tri] << "\t";

	cout << "\n\ttris = ";
	for (e=1; e<m; e++)
		cout << tris[e] << "\t";

	cout << "\n\tpos  = ";
	for (e=1; e<m; e++)
		cout << pos[e] << "\t";
	cout << "\n" <<endl;
}

//
//
///**
// * @brief   Compute the maximum triangle core
// */
//long long graphpack_graph::triangle_max_core() {
//	long long m = num_edges();
//	max_tri_core = -1;
//
//	for (long long e=0; e<m; e++) {
//	    if (e < 20)
//	        cout << tri_core[e] << " ";
//		if (tri_core[e] > max_tri_core)
//			max_tri_core = tri_core[e];
//	}
//	return max_tri_core;
//}




/**
 * @brief   create a reverse lookup index on the edges
 *          given a ptr to a neighbor in the compressed
 *          sparse representation, what vertices are stored
 *          in that position. Namely, what neighbor does the ptr
 *          correspond to, and for which vertex is it a neighbor of..
 *          For this we design a fast data structure that gives us
 *          O(1) time for retrieving the vertices of an edge.
 */

void graphpack_graph::map_edges_serial() {
    // edges (instead of vert ID store EID into e_u and e_v....
    long long i, j, u, v;
    long long m = edges.size();

    eid.resize(m);

    for (i=0; i < m/2; i++) {
        v = e_v[i];    u = e_u[i];

        for (j=vertices[v]; j<vertices[v+1]; j++) {
            if (edges[j] == u) {
                eid[j] = i;
//                eid[j] = i + 1;
                break;
            }
        }

        for (j=vertices[u]; j<vertices[u+1]; j++) {
            if (edges[j] == v) {
                eid[j] = i;
//                eid[j] = i + 1;
                break;
            }
        }
    }
}

/**
 * @brief   create a reverse lookup index on the edges
 *          given a ptr to a neighbor in the compressed
 *          sparse representation, what vertices are stored
 *          in that position. Namely, what neighbor does the ptr
 *          correspond to, and for which vertex is it a neighbor of..
 *          For this we design a fast data structure that gives us
 *          O(1) time for retrieving the vertices of an edge.
 */
void graphpack_graph::map_edges_parallel() {
    // edges (instead of vert ID store EID into e_u and e_v....
    long long i, j, u, v;
    long long m = edges.size();

    eid.resize(m);
#pragma omp parallel for schedule(dynamic) private(i,j,u,v)
    for (i=0; i < m/2; i++) {
        v = e_v[i];    u = e_u[i];

        for (j=vertices[v]; j<vertices[v+1]; j++) {
            if (edges[j] == u) {
                eid[j] = i + 1;
                break;
            }
        }

        for (j=vertices[u]; j<vertices[u+1]; j++) {
            if (edges[j] == v) {
                eid[j] = i + 1;
                break;
            }
        }
    }
}


void graphpack_graph::map_eid_to_vid(bool is_parallel) {
    // edges (instead of vert ID store EID into e_u and e_v....
    long long i, j, u, v;
    long long m = edges.size();

    eid.resize(m);
    if (is_parallel) {
        #pragma omp parallel for schedule(dynamic) private(i,j,u,v)
        for (i=0; i < m/2; i++) {
            v = e_v[i];    u = e_u[i];

            for (j=vertices[v]; j<vertices[v+1]; j++) {
                if (edges[j] == u) {
                    eid[j] = i + 1;
                    break;
                }
            }

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                if (edges[j] == v) {
                    eid[j] = i + 1;
                    break;
                }
            }
        }
    }
    else {
        for (i=0; i < m/2; i++) {
            v = e_v[i];    u = e_u[i];

            for (j=vertices[v]; j<vertices[v+1]; j++) {
                if (edges[j] == u) {
                    eid[j] = i + 1;
                    break;
                }
            }

            for (j=vertices[u]; j<vertices[u+1]; j++) {
                if (edges[j] == v) {
                    eid[j] = i + 1;
                    break;
                }
            }
        }
    }
}





//
//void graphpack_graph::map_edges() {
//	// edges (instead of vert ID store EID into e_u and e_v....
//	long long i, j, u, v;
//	long long m = edges.size();
//
//	eid.resize(m);
//#pragma omp parallel for schedule(dynamic) private(i,j,u,v)
//	for (i=0; i < m/2; i++) {
//		v = e_v[i];
//		u = e_u[i];
//
//		for (j=vertices[v]; j<vertices[v+1]; j++) {
//			if (edges[j] == u) {
//				eid[j] = i;
//				break;
//			}
//		}
//
//		for (j=vertices[u]; j<vertices[u+1]; j++) {
//			if (edges[j] == v) {
//				eid[j] = i;
//				break;
//			}
//		}
//	}
//}












/**
 * Used as a subroutine in branch and bound algorithms
 */


/*
 * @brief   Parallelized Triangle Core Decomposition
 *          Runtime = O(m^1.5), majority of work is for computing the triangle counts
 *
 * Note the following:
 * Given an edge_id (or neighbor id of some vertex), e_v and e_u return the vertex v and u, respectively.
 * Also given the position of a neighbor in the edge array, eid returns the edge_id
 *
 * Note that e_v and e_u are constructed when the graph is initially read.
 * If these functions are called after the graph is reduced,
 * then e_u and e_v must be updated, and eid must be recomputed.
 *
 * @See also: mcbound_triangle_cores.cpp
 *
 * Author: Ryan A. Rossi
 *
 */

/** default to the largest degree strategy */
void graphpack_graph::triangle_cores_approx() {
    string strat = "degree";
    triangle_cores_approx(strat);
}

void graphpack_graph::triangle_cores_approx(string strategy) {
//void graphpack_graph::triangle_cores_approx() {

    double total_sec = get_time();
    // given position in edge array (neighbor), eid returns the unique edge id
    map_edges_parallel();

//    double triangle_core_time;
//    double fast_triangle_core_time;
//    double edge_triangles;
//    double fast_edge_triangles;

    double sec = get_time();
    // number of triangles each edge participates
    if (is_dense)
        parallel_triangle_counting_adj();
    else
        parallel_triangle_counting();

    cout << "triangles: " << total_t << " (time: " << get_time() - sec << "), " <<endl;

    sec = get_time();
    // compute the triangle cores
    triangle_core_numbers();
    cout << "max triangle core: " << max_tri_core+2 << " (time: " << get_time() - sec << "), " <<endl;

    cout << "total time of triangle cores: " << get_time() - total_sec;
    cout << ", max triangle core: " << max_tri_core+2 <<endl;
}

void graphpack_graph::triangle_cores_approx(double & triangle_core_perf, double & triangle_perf) {

    // given position in edge array (neighbor), eid returns the unique edge id
    map_edges_parallel();

     triangle_core_perf = get_time();
     triangle_perf = get_time();
    // number of triangles each edge participates
     // number of triangles each edge participates
      if (is_dense)
          parallel_triangle_counting_adj();
      else
          parallel_triangle_counting();

    triangle_perf = get_time() - triangle_perf;

    triangle_core_numbers();
    triangle_core_perf = get_time() - triangle_core_perf;
}



void graphpack_graph::triangle_core_numbers_ordering() {
    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges() + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
//    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    tris_ordering.reserve(m);
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed


    for (e=1; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;


    /** compute starting positions for the edges in each bin */
    start = 1;
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    /** bucket sort the edges by their triangle counts
      *  place the edges into the correct buckets */
    for (e=1; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris_ordering[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    /** fixup the bin for the main computation */
    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 1;

    /**
     * Process each edge in tris_ordering in the order of number of triangles
     * Note that m here is actually num_edges() + 1
     */
    vector<long long> ind(n, 0);
    for (e=1; e<m; e++) {

        /** get the edge_id */
        eid_vu = tris_ordering[e]; // @TODO: Could save the ordering..

        /** get the vertex ids of each end of the edge */
        v = e_v[eid_vu-1];
        u = e_u[eid_vu-1];


        /**
         * Create lookup table for the neighbors of v,
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, indicating the exact neighbor.
         */
        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i;
        ind[v] = v;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */
        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
            w = edges[j];
            if (w==v) { continue; }  //self-loop
            if (w==u) { continue; }  //discard 1-hop self-loop

            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
            if (ind[w]) {
                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
                eid_vw = eid[ind[w]];   // get the edge_id of v--w

                if (!proc[eid_uw] && !proc[eid_vw]) {
                    edge_neigh[0] = eid_uw;
                    edge_neigh[1] = eid_vw;
                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {

                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris_ordering[px];                 //starting edge id

                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris_ordering[pw] = e_x;
                                pos[e_x] = pw;
                                tris_ordering[px] = edge_neigh[k];
                            }
                            bin[tw]++;   // update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
            ind[edges[i]] = 0;
        ind[v] = 0;

        proc[eid_vu] = 1; // set edge to be processed
    }

#ifndef DEBUG
        cout << "max triangle core = " << tri_core[tris_ordering[m-1]] <<endl;
#endif

        max_tri_core = tri_core[tris_ordering[m-1]];

//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos

        bin.clear();
        pos.clear();
//        tris.clear();
        ind.clear();
}







void graphpack_graph::parallel_triangle_counting() {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);

    #pragma omp parallel for schedule(dynamic) \
        shared(tmp_tri_core, tmp_max_t_edge) firstprivate(ind) private(v,u,w,e,i,j)
    for (e=0; e<m; e++) {
        v = e_v[e];   u = e_u[e];
//        if (u==v) { continue; }
        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 1;
        ind[v] = 1;

        for (j=vertices[u]; j<vertices[u+1]; j++) {
            w = edges[j];
            if (w==u) { continue; }  //self-loop
            if (w==v) { continue; }  //discard 1-hop self-loop
            if (ind[w]) {
                tmp_tri_core[e+1] += 1;
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
        ind[v] = 0;
    }

    max_t_edge = 0;
    total_t = 0;
    for (e=0; e<m; e++) {
        total_t += tmp_tri_core[e+1];
        if (tmp_tri_core[e+1] > max_t_edge) {
            max_t_edge = tmp_tri_core[e+1];
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}


/*
 * @brief   ADJ data structure basically removes the cost
 *          of creating the fast index structure!
 *
 */
void graphpack_graph::parallel_triangle_counting_adj() {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
//    vector<long long> ind(n, 0);

    #pragma omp parallel for schedule(dynamic) \
        shared(tmp_tri_core, tmp_max_t_edge) private(v,u,w,e,i,j)
    for (e=0; e<m; e++) { // v--u
        v = e_v[e];   u = e_u[e];
//        if (u==v) { continue; }

        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
            w = edges[j];
            if (w==u) { continue; }  //self-loop
            if (w==v) { continue; }  //discard 1-hop self-loop
            if (adj[v][w]) {
                tmp_tri_core[e+1] += 1;
            }
        }
    }

    max_t_edge = 0;
    total_t = 0;

    for (e=0; e<m; e++) {
        total_t += tmp_tri_core[e+1];
        if (tmp_tri_core[e+1] > max_t_edge) {
            max_t_edge = tmp_tri_core[e+1];
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
//    ind.clear();
}




void graphpack_graph::triangle_core_numbers() {
    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges() + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed


    for (e=1; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;


    /** compute starting positions for the edges in each bin */
    start = 1;
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    /** bucket sort the edges by their triangle counts
      *  place the edges into the correct buckets */
    for (e=1; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    /** fixup the bin for the main computation */
    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 1;

    /**
     * Process each edge in tris in the order of number of triangles
     * Note that m here is actually num_edges() + 1
     */
    vector<long long> ind(n, 0);
    for (e=1; e<m; e++) {

        /** get the edge_id */
        eid_vu = tris[e]; // @TODO: Could save the ordering..

        /** get the vertex ids of each end of the edge */
        v = e_v[eid_vu-1];
        u = e_u[eid_vu-1];


        /**
         * Create lookup table for the neighbors of v,
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, indicating the exact neighbor.
         */
        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i;
        ind[v] = v;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */
        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
            w = edges[j];
            if (w==v) { continue; }  //self-loop
            if (w==u) { continue; }  //discard 1-hop self-loop

            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
            if (ind[w]) {
                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
                eid_vw = eid[ind[w]];   // get the edge_id of v--w

                if (!proc[eid_uw] && !proc[eid_vw]) {
                    edge_neigh[0] = eid_uw;
                    edge_neigh[1] = eid_vw;
                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {

                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris[px];                 //starting edge id

                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris[pw] = e_x;
                                pos[e_x] = pw;
                                tris[px] = edge_neigh[k];
                            }
                            bin[tw]++;   // update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
            ind[edges[i]] = 0;
        ind[v] = 0;

        proc[eid_vu] = 1; // set edge to be processed
    }

#ifndef DEBUG
        cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
#endif

        tris_ordering = tris;
        max_tri_core = tri_core[tris[m-1]];

//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos

        bin.clear();
        pos.clear();
//        tris.clear();
        ind.clear();
}






/**
 * Optimized Triangle Core Numbers (Slightly)
 */

void graphpack_graph::triangle_core_numbers_optimized() {
    long long e, e_ux, e_x, tri, i, j, k;
    long long start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges() + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2,0);  //exactly two neighboring edges of e form a triangle
    vector<short> proc(m,0); //a triangle is processed, if any one of its edges are processed


    for (e=1; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;

    start = 1;  // bucket sort
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    for (e=1; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 1;


    vector<long long> ind(n, 0);

    // only allows two edges
    for (e = 1; e < m; e++) { //for each edge in tris in the order of number of triangles

        v = e_v[tris[e]-1];
        u = e_u[tris[e]-1];

        /**
         * Create lookup table for the neighbors of v,
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, indicating the exact neighbor.
         */
        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i;
        ind[v] = v;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */
        for (j=vertices[u]; j<vertices[u+1]; j++) { // for neigh of u

            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
            if (ind[edges[j]]) { // if neigh of v
                if (edges[j]==v) { continue; }  //self-loop
                if (edges[j]==u) { continue; }  //discard 1-hop self-loop

                // check if neighboring edges uw and vw have been marked
                if (proc[eid[j]] == 0 && proc[eid[ind[edges[j]]]] == 0) {
                    edge_neigh[0] = eid[j];
                    edge_neigh[1] = eid[ind[edges[j]]];
                    for (k=0; k<2; k++) { //each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[tris[e]]) {
                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris[px];                 //starting edge id
                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris[pw] = e_x;
                                pos[e_x] = pw;
                                tris[px] = edge_neigh[k];
                            }
                            bin[tw]++;   //update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
        ind[v] = 0;

        proc[tris[e]] = 1; //set edge to be processed
    }

    cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
    max_tri_core = tri_core[tris[m-1]];

    write_vector(tri_core,"-tcore.txt");

    triangle_max_core();

    bin.clear();
    pos.clear();
    tris.clear();
    ind.clear();
}







/**
 * Triangle Cores using Adj
 */

/** default to the largest degree strategy */
void graphpack_graph::triangle_cores_approx_adj() {
    string strat = "degree";
    triangle_cores_approx_adj(strat);
}

void graphpack_graph::triangle_cores_approx_adj(string strategy) {

    // given position in edge array (neighbor), eid returns the unique edge id
    double total_sec = get_time();
    map_edges_parallel();

    double sec = get_time();
    // number of triangles each edge participates
    parallel_triangle_counting_adj();
    cout << "triangles: " << total_t << " (time: " << get_time() - sec << "), " <<endl;

    sec = get_time();
    // compute the triangle cores
    triangle_core_numbers();
    cout << "max triangle core: " << max_tri_core+2 << " (time: " << get_time() - sec << "), " <<endl;

    cout << "total time of triangle cores: " << get_time() - total_sec;
    cout << ", max triangle core: " << max_tri_core+2 <<endl;

}


void graphpack_graph::triangle_cores_approx_adj(double & triangle_core_perf, double & triangle_perf) {
//void graphpack_graph::triangle_cores_approx() {

    double total_sec = get_time();
    // given position in edge array (neighbor), eid returns the unique edge id
    map_edges_parallel();

//    double triangle_core_time;
//    double fast_triangle_core_time;
//    double edge_triangles;
//    double fast_edge_triangles;

     triangle_core_perf = get_time();
     triangle_perf = get_time();
    // number of triangles each edge participates
     // number of triangles each edge participates
     if (is_dense)
         parallel_triangle_counting_adj();
     else
         parallel_triangle_counting();
    triangle_perf = get_time() - triangle_perf;
//    cout << "triangles: " << total_t << " (time: " << get_time() - sec << "), " <<endl;

//    sec = get_time();
    // compute the triangle cores
    triangle_core_numbers();
    triangle_core_perf = get_time() - triangle_core_perf;

//    cout << "max triangle core: " << max_tri_core+2 << " (time: " << get_time() - sec << "), " <<endl;
//
//    cout << "total time of triangle cores: " << get_time() - total_sec;
//    cout << ", max triangle core: " << max_tri_core+2 <<endl;
}








/*
 * @brief    TRIANGLE NEIGHBORHOOD CORES (MAIN PURPOSE)
 * This function is customized for neighborhood cores
 * Assuming we can't compute triangles in parallel or compute the edge map in parallel
 */
void graphpack_graph::triangle_cores_serial_adj() {

    // given position in edge array (neighbor), eid returns the unique edge id
    double total_sec = get_time();
    map_edges_serial();

    double triangle_time = get_time();
    // number of triangles each edge participates
    triangle_counting_adj();
    triangle_time = get_time() - triangle_time;
    //cout << "triangles (serial): " << total_t << " (time: " << get_time() - triangle_time << "), " <<endl;

    double tcore_time = get_time();
    // compute the triangle cores
    triangle_core_numbers();
    tcore_time = get_time() - tcore_time;
    //cout << "max triangle core: " << max_tri_core+2 << " (time: " << tcore_time << "), " <<endl;

#ifndef DEBUG_TCORE
    cout << "max triangle core: " << max_tri_core+2 << ", ";
    cout << "total time of triangle cores (fully serial): " << get_time() - total_sec << " sec, ";
    cout << "triangle counting time: " << triangle_time << " sec, ";
    cout << "triangle core time: " << tcore_time << " sec" <<endl;
#endif

}


/*
 * This function's main purpose is for neighborhood triangle cores.
 *
 * Note: ADJ data structure basically removes the cost
 * of creating the fast index structure!
 */
void graphpack_graph::triangle_counting_adj() {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);

    for (e=0; e<m; e++) { // v--u
        v = e_v[e];   u = e_u[e];
        if (u==v) { continue; }

        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
            w = edges[j];
            if (w==u) { continue; }  //self-loop
            if (w==v) { continue; }  //discard 1-hop self-loop
            if (adj[v][w]) {
                tmp_tri_core[e+1] += 1;
            }
        }
    }

    max_t_edge = 0;
    total_t = 0;
    for (e=0; e<m; e++) {
        total_t += tmp_tri_core[e+1];
        if (tmp_tri_core[e+1] > max_t_edge) {
            max_t_edge = tmp_tri_core[e+1];
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}


/*
 * How many triangles must an edge participate in to contain a clique of size k?
 * For edges, without duplicate triangles, then: T(v)-2 triangles, so
 * T(v)-2 > mc
 *
 * @brief   Idea is to use adj, and what we already know from the number of triangles, to
 *          speedup the algorithm by removing vertices, etc...
 */

void graphpack_graph::triangle_core_numbers(int & lb) {
    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges() + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed


    for (e=1; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;


    /** compute starting positions for the edges in each bin */
    start = 1;
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    /** bucket sort the edges by their triangle counts
      *  place the edges into the correct buckets */
    for (e=1; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    /** fixup the bin for the main computation */
    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 1;

    /**
     * Process each edge in tris in the order of number of triangles
     * Note that m here is actually num_edges() + 1
     */
    vector<long long> ind(n, 0);
    for (e=1; e<m; e++) {

        /** get the edge_id */
        eid_vu = tris[e]; // @TODO: Could save the ordering..

        /** get the vertex ids of each end of the edge */
        v = e_v[eid_vu-1];
        u = e_u[eid_vu-1];


        /**
         * Create lookup table for the neighbors of v,
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, indicating the exact neighbor.
         */
        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i;
        ind[v] = v;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */
        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
            w = edges[j];
            if (w==v) { continue; }  //self-loop
            if (w==u) { continue; }  //discard 1-hop self-loop

            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
            if (ind[w]) {
                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
                eid_vw = eid[ind[w]];   // get the edge_id of v--w

                if (!proc[eid_uw] && !proc[eid_vw]) {
                    edge_neigh[0] = eid_uw;
                    edge_neigh[1] = eid_vw;
                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {

                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris[px];                 //starting edge id

                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris[pw] = e_x;
                                pos[e_x] = pw;
                                tris[px] = edge_neigh[k];
                            }
                            bin[tw]++;   // update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
            ind[edges[i]] = 0;
        ind[v] = 0;

        proc[eid_vu] = 1; // set edge to be processed
    }

#ifndef DEBUG
        cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
#endif

        max_tri_core = tri_core[tris[m-1]];

//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos

        bin.clear();
        pos.clear();
        tris.clear();
        ind.clear();
}










/**
 * START OF CODES FOR FAST MAXIMUM CLIQUE PRUNING VIA TRIANGLE CORES
 * @param mc
 */



/*
 * @brief TRIANGLE NEIGHBORHOOD CORES (MAIN PURPOSE)
 * This function is customized for neighborhood cores
 * Assuming we can't compute triangles in parallel or compute the edge map in parallel
 */
void graphpack_graph::triangle_cores_serial_adj(int & mc) {

    // given position in edge array (neighbor), eid returns the unique edge id
#ifndef DEBUG_TCORE
    double total_sec = get_time();
#endif
    map_edges_serial();


    vector<int> pruned_edges(num_edges()+1,0);

#ifndef DEBUG_TCORE
    double triangle_time = get_time();
#endif
    // number of triangles each edge participates
    triangle_counting_adj(mc,pruned_edges);
#ifndef DEBUG_TCORE
    triangle_time = get_time() - triangle_time;
#endif
    //cout << "triangles (serial): " << total_t << " (time: " << get_time() - triangle_time << "), " <<endl;

#ifndef DEBUG_TCORE
    double tcore_time = get_time();
#endif
    // compute the triangle cores
    triangle_core_numbers(mc,pruned_edges);

#ifndef DEBUG_TCORE
    tcore_time = get_time() - tcore_time;
#endif
    //cout << "max triangle core: " << max_tri_core+2 << " (time: " << tcore_time << "), " <<endl;

#ifndef DEBUG_TCORE
    cout << "max triangle core: " << max_tri_core+2 << ", ";
    cout << "total time of triangle cores (fully serial): " << get_time() - total_sec << " sec, ";
    cout << "triangle counting time: " << triangle_time << " sec, ";
    cout << "triangle core time: " << tcore_time << " sec" <<endl;
#endif

}

//
//
//
///*
// * This function's main purpose is for neighborhood triangle cores.
// *
// * Note: ADJ data structure basically removes the cost
// * of creating the fast index structure!
// */
///**
// *
// * @param mc            pointer to the current size of the max clique,
// *                      this value may be updated dynamically by another
// *                      worker (as larger cliques are found and communicated)
// *
// * @param pruned_edges  edges are set to be pruned dynamically after
// *                      computing this value. An edge is pruned if
// *                      it does not have at least omega+2 triangles.
// *                      Since if it has less, than that edge must not
// *                      take part in a clique larger than omega.
// *                      This computation ensures additional unnecessary
// *                      work is not carried out in triangle cores,
// *                      since these edges are removed.
// *                      This sort of dynamic computation
// */
//void graphpack_graph::triangle_counting_adj(int & mc, vector<int> & pruned_edges) {
//    long long e, v, u, w, i, j;
//    long long n = num_vertices(), m = num_edges();
//
//    long long tmp_max_t_edge = 0;
//    vector<long long> tmp_tri_core(m+1,0);
//    vector<long long> ind(n, 0);
//
//    #pragma omp parallel for schedule(dynamic) \
//        shared(tmp_tri_core, tmp_max_t_edge) private(v,u,w,e,i,j)
//    for (e=0; e<m; e++) { // v--u
//        v = e_v[e];   u = e_u[e];
//        if (u==v) { continue; }
//
//        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
//
//            /**
//             * @brief       if at any point, there are not enough neighbors to
//             *              to satisfy the bound, then we terminate the counting
//             *              early, rather than doing unnecessary work. Note that
//             *              (vertices[u+1] - j) is an upper-bound on the number
//             *              of possible triangles remaining for that vertex
//             *
//             *              Advantages
//             *              ----------
//             *              * Works well for graphs where the max clique is large, then
//             *              we have many early terminations, which can significantly
//             *              increase the performance (and cost of using such a tight bound)
//             */
//            if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
//                tmp_tri_core[e+1] = 0;
//                pruned_edges[e+1] = 1;
//                break;
//            }
//
//            if (pruned_edges[eid[j]+1])
//                continue;
//
//            w = edges[j];
//
//            if (w==u) { continue; }  //self-loop
//            if (w==v) { continue; }  //discard 1-hop self-loop
//            if (adj[v][w]) {
//                tmp_tri_core[e+1] += 1;
//            }
//        }
//        if (tmp_tri_core[e+1]+2 < mc) {
//            tmp_tri_core[e+1] = 0;
//            pruned_edges[e+1] = 1;
//            adj[v][u] = 0;
//            adj[u][v] = 0;
//        }
//    }
//
//
//
//    max_t_edge = 0;
//    total_t = 0;
//    for (e=0; e<m; e++) {
//        total_t += tmp_tri_core[e+1];
//        if (tmp_tri_core[e+1] > max_t_edge) {
//            max_t_edge = tmp_tri_core[e+1];
//        }
//    }
//
//    tri_core = tmp_tri_core;
//
//    tmp_tri_core.clear();
//    ind.clear();
//}






/*
 * This function's main purpose is for neighborhood triangle cores.
 *
 * Note: ADJ data structure basically removes the cost
 * of creating the fast index structure!
 */
/**
 *
 * @param mc            pointer to the current size of the max clique,
 *                      this value may be updated dynamically by another
 *                      worker (as larger cliques are found and communicated)
 *
 * @param pruned_edges  edges are set to be pruned dynamically after
 *                      computing this value. An edge is pruned if
 *                      it does not have at least omega+2 triangles.
 *                      Since if it has less, than that edge must not
 *                      take part in a clique larger than omega.
 *                      This computation ensures additional unnecessary
 *                      work is not carried out in triangle cores,
 *                      since these edges are removed.
 *                      This sort of dynamic computation
 */
void graphpack_graph::triangle_counting_adj(int & mc, vector<int> & pruned_edges, bool is_parallel) {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);

    if (is_parallel) {
        #pragma omp parallel for schedule(dynamic) \
        shared(tmp_tri_core, tmp_max_t_edge) private(v,u,w,e,i,j)
        for (e=0; e<m; e++) { // v--u
            v = e_v[e];   u = e_u[e];
            if (u==v) { continue; }

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                /**
                 * @brief       if at any point, there are not enough neighbors to
                 *              to satisfy the bound, then we terminate the counting
                 *              early, rather than doing unnecessary work. Note that
                 *              (vertices[u+1] - j) is an upper-bound on the number
                 *              of possible triangles remaining for that vertex
                 *
                 *              Advantages
                 *              ----------
                 *              * Works well for graphs where the max clique is large, then
                 *              we have many early terminations, which can significantly
                 *              increase the performance (and cost of using such a tight bound)
                 */
                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    tmp_tri_core[e+1] = 0;
                    pruned_edges[e+1] = 1;
                    break;
                }

                if (pruned_edges[eid[j]+1])
                    continue;

                w = edges[j];

                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (adj[v][w]) {
                    tmp_tri_core[e+1] += 1;
                }
            }
            if (tmp_tri_core[e+1]+2 < mc) {
                tmp_tri_core[e+1] = 0;
                pruned_edges[e+1] = 1;
                adj[v][u] = 0;
                adj[u][v] = 0;
            }
        }
    }
    else {
        for (e=0; e<m; e++) { // v--u
            v = e_v[e];   u = e_u[e];
            if (u==v) { continue; }

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                /**
                 * @brief       if at any point, there are not enough neighbors to
                 *              to satisfy the bound, then we terminate the counting
                 *              early, rather than doing unnecessary work. Note that
                 *              (vertices[u+1] - j) is an upper-bound on the number
                 *              of possible triangles remaining for that vertex
                 *
                 *              Advantages
                 *              ----------
                 *              * Works well for graphs where the max clique is large, then
                 *              we have many early terminations, which can significantly
                 *              increase the performance (and cost of using such a tight bound)
                 */
                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    tmp_tri_core[e+1] = 0;
                    pruned_edges[e+1] = 1;
                    break;
                }

                if (pruned_edges[eid[j]+1])
                    continue;

                w = edges[j];

                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (adj[v][w]) {
                    tmp_tri_core[e+1] += 1;
                }
            }
            if (tmp_tri_core[e+1]+2 < mc) {
                tmp_tri_core[e+1] = 0;
                pruned_edges[e+1] = 1;
                adj[v][u] = 0;
                adj[u][v] = 0;
            }
        }
    }


    max_t_edge = 0;
    total_t = 0;
    for (e=0; e<m; e++) {
        total_t += tmp_tri_core[e+1];
        if (tmp_tri_core[e+1] > max_t_edge) {
            max_t_edge = tmp_tri_core[e+1];
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}

void graphpack_graph::triangle_counting(int & mc, vector<int> & pruned_edges, bool is_parallel) {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);


    if (is_parallel) {
        #pragma omp parallel for schedule(dynamic) \
        shared(tmp_tri_core, tmp_max_t_edge) private(v,u,w,e,i,j)
        for (e=0; e<m; e++) { // v--u
            v = e_v[e];   u = e_u[e];
//            if (u==v) { continue; }

            //        for (i=vertices[v]; i<vertices[v+1]; i++)
            //            ind[edges[i]] = 1;
            //        ind[v] = 1;

            for (i=vertices[v]; i<vertices[v+1]; i++)
                ind[edges[i]] = i; // store the ptr location in the edge array
//            ind[v] = v;

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                /**
                 * @brief       if at any point, there are not enough neighbors to
                 *              to satisfy the bound, then we terminate the counting
                 *              early, rather than doing unnecessary work. Note that
                 *              (vertices[u+1] - j) is an upper-bound on the number
                 *              of possible triangles remaining for that vertex
                 *
                 *              Advantages
                 *              ----------
                 *              * Works well for graphs where the max clique is large, then
                 *              we have many early terminations, which can significantly
                 *              increase the performance (and cost of using such a tight bound)
                 */
                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    tmp_tri_core[e+1] = 0;
                    pruned_edges[e+1] = 1;

                    break;        // avoid computing the rest
                }

                w = edges[j];
//                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (ind[w]) {
                    tmp_tri_core[e+1] += 1;
                }
            }


            if (tmp_tri_core[e+1]+2 < mc) {
                tmp_tri_core[e+1] = 0;
                pruned_edges[e+1] = 1;

                //            update_neighbors(tmp_tri_core);
            }

            for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            ind[v] = 0;
        }
    }
    else {
        for (e=0; e<m; e++) { // v--u
            v = e_v[e];   u = e_u[e];
            if (u==v) { continue; }

            //        for (i=vertices[v]; i<vertices[v+1]; i++)
            //            ind[edges[i]] = 1;
            //        ind[v] = 1;

            for (i=vertices[v]; i<vertices[v+1]; i++)
                ind[edges[i]] = i; // store the ptr location in the edge array
            ind[v] = v;

            for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

                /**
                 * @brief       if at any point, there are not enough neighbors to
                 *              to satisfy the bound, then we terminate the counting
                 *              early, rather than doing unnecessary work. Note that
                 *              (vertices[u+1] - j) is an upper-bound on the number
                 *              of possible triangles remaining for that vertex
                 *
                 *              Advantages
                 *              ----------
                 *              * Works well for graphs where the max clique is large, then
                 *              we have many early terminations, which can significantly
                 *              increase the performance (and cost of using such a tight bound)
                 */
                if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                    tmp_tri_core[e+1] = 0;
                    pruned_edges[e+1] = 1;

                    break;        // avoid computing the rest
                }

                w = edges[j];
                if (w==u) { continue; }  //self-loop
                if (w==v) { continue; }  //discard 1-hop self-loop
                if (ind[w]) {
                    tmp_tri_core[e+1] += 1;
                }
            }


            if (tmp_tri_core[e+1]+2 < mc) {
                tmp_tri_core[e+1] = 0;
                pruned_edges[e+1] = 1;

                //            update_neighbors(tmp_tri_core);
            }

            for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            ind[v] = 0;
        }
    }


    if (is_gstats) { // don't compute unless specified
        max_t_edge = 0;
        total_t = 0;
        for (e=0; e<m; e++) {
            total_t += tmp_tri_core[e+1];
            if (tmp_tri_core[e+1] > max_t_edge) {
                max_t_edge = tmp_tri_core[e+1];
            }
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}




/**
 * @brief           sparse parallel edge-centric triangle counting
 *                  these codes have been extremely optimized for
 *                  use within branch-and-bound algorithms. In this case,
 *                  mc can be simply treated as the bound.
 *
 * @param mc
 * @param pruned_edges
 */
void graphpack_graph::triangle_counting_index(int & mc, vector<int> & pruned_edges) {
    long long e, v, u, w, i, j;
    long long n = num_vertices(), m = num_edges();

    long long tmp_max_t_edge = 0;
    vector<long long> tmp_tri_core(m+1,0);
    vector<long long> ind(n, 0);

    for (e=0; e<m; e++) { // v--u
        v = e_v[e];   u = e_u[e];
        if (u==v) { continue; }

//        for (i=vertices[v]; i<vertices[v+1]; i++)
//            ind[edges[i]] = 1;
//        ind[v] = 1;

        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i; // store the ptr location in the edge array
        ind[v] = v;

        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w

            /**
             * @brief       if at any point, there are not enough neighbors to
             *              to satisfy the bound, then we terminate the counting
             *              early, rather than doing unnecessary work. Note that
             *              (vertices[u+1] - j) is an upper-bound on the number
             *              of possible triangles remaining for that vertex
             *
             *              Advantages
             *              ----------
             *              * Works well for graphs where the max clique is large, then
             *              we have many early terminations, which can significantly
             *              increase the performance (and cost of using such a tight bound)
             */
            if ((tmp_tri_core[e+1]+2) + (vertices[u+1] - j) < mc) {
                tmp_tri_core[e+1] = 0;
                pruned_edges[e+1] = 1;

                break;        // avoid computing the rest
            }

            w = edges[j];
            if (w==u) { continue; }  //self-loop
            if (w==v) { continue; }  //discard 1-hop self-loop
            if (ind[w]) {
                tmp_tri_core[e+1] += 1;
            }
        }


        if (tmp_tri_core[e+1]+2 < mc) {
            tmp_tri_core[e+1] = 0;
            pruned_edges[e+1] = 1;

//            update_neighbors(tmp_tri_core);
        }

        for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
        ind[v] = 0;
    }


    if (is_gstats) { // don't compute unless specified
        max_t_edge = 0;
        total_t = 0;
        for (e=0; e<m; e++) {
            total_t += tmp_tri_core[e+1];
            if (tmp_tri_core[e+1] > max_t_edge) {
                max_t_edge = tmp_tri_core[e+1];
            }
        }
    }

    tri_core = tmp_tri_core;

    tmp_tri_core.clear();
    ind.clear();
}




/**
 * @brief               Applys a primitive operation (addition,
 *                      subtract, or a custom operation via a function)
 *                      to a vertices neighbors. Used mainly for dynamically updating
 *                      neighbor counts when a vertex (or edge) is removed from the graph.
 */
//void graphpack_graph::neighbor_operation(vector<long long> & tmp_tri_core, int & u, int & u) {
////    //! if a vertex (or edge) is pruned, then all my neighboring edges,
////    //! that use me to form a triangle should not be allowed to do so, or
////    //! we must subtract one.
////
////    //! from the edge perspective, how do we update this in linear time?
////    //! we must be able to index directly the neighboring edges?
////
//    //                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
//    //                eid_vw = eid[ind[w]];   // get the edge_id of v--w
//
//}


/*
 * How many triangles must an edge participate in to contain a clique of size k?
 * For edges, without duplicate triangles, then: T(v)-2 triangles, so
 * T(v)-2 > mc
 *
 * TODO: Idea is to use adj, and what we already know from the number of triangles, to
 *       speedup the algorithm by removing vertices, etc...
 */

void graphpack_graph::triangle_core_numbers(int & mc, vector<int> & pruned_edges) {
    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
    long long n = num_vertices(), m = num_edges() + 1;
    long long mt_end = max_t_edge+1;

    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
    vector <long long> pos(m);        // position of the edges in tris

    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed


    for (e=1; e<m; e++) // bin edges by triangle counts
        bin[tri_core[e]]++;


    /** compute starting positions for the edges in each bin */
    start = 1;
    for (tri=0; tri < mt_end; tri++) {
        num = bin[tri];
        bin[tri] = start;
        start = start + num;
    }

    /** bucket sort the edges by their triangle counts
      *  place the edges into the correct buckets */
    for (e=1; e<m; e++) {
        pos[e] = bin[tri_core[e]];
        tris[pos[e]] = e;
        bin[tri_core[e]]++;
    }

    /** fixup the bin for the main computation */
    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
    bin[0] = 1;

    /**
     * Process each edge in tris in the order of number of triangles
     * Note that m here is actually num_edges() + 1
     *
     * @todo:   make all funcs consistent, there doesn't seem to be a
     *          reason for the edge ids to start at 1
     */
    vector<long long> ind(n, 0);
    for (e=1; e<m; e++) {

        if (pruned_edges[e])
            continue;

        /** get the edge_id */
        eid_vu = tris[e]; // @TODO: Could save the ordering..

        /** get the vertex ids of each end of the edge */
        v = e_v[eid_vu-1];
        u = e_u[eid_vu-1];


        /**
         * Create lookup table for the neighbors of v,
         * But instead of only setting the lookup value to 1,
         * We set it to be the edge_id,
         * or the position in the edge list, indicating the exact neighbor.
         */
        for (i=vertices[v]; i<vertices[v+1]; i++)
            ind[edges[i]] = i;
        ind[v] = v;

        /**
         * Get a neighbor of u, so that we have a wedge: v--u--w
         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
         */
        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
            if (pruned_edges[eid[j]+1])
                continue;

            w = edges[j];
            if (w==v) { continue; }  //self-loop
            if (w==u) { continue; }  //discard 1-hop self-loop


            /**
             * v--w, if w is a neigh of v,
             * so if true, then we know: v--u, u--w, w--v forming a triangle.
             */
            if (ind[w]) {
                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
                eid_vw = eid[ind[w]];   // get the edge_id of v--w

                if (!proc[eid_uw] && !proc[eid_vw]) {
                    edge_neigh[0] = eid_uw;
                    edge_neigh[1] = eid_vw;
                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {

                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
                            pw = pos[edge_neigh[k]];        //actual pos of w
                            px = bin[tw];                   //starting position of tw
                            e_x = tris[px];                 //starting edge id

                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
                                pos[edge_neigh[k]] = px;   //starting position of e_
                                tris[pw] = e_x;
                                pos[e_x] = pw;
                                tris[px] = edge_neigh[k];
                            }
                            bin[tw]++;   // update starting pos
                            tri_core[edge_neigh[k]]--;
                        }
                    }
                }
            }
        }

        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
            ind[edges[i]] = 0;
        ind[v] = 0;

        proc[eid_vu] = 1; // set edge to be processed
    }

#ifndef DEBUG
        cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
#endif

        max_tri_core = tri_core[tris[m-1]];

//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos

        bin.clear();
        pos.clear();
        tris.clear();
        ind.clear();
}





//
///**
// * Purpose of function is to prune the edges and vertices via triangle cores
// * To make the termination of the branch function as fast as possible.
// * Note that for DIMACs, the price of this computation is insignificant (for the harder graphs)
// *
// * Function returns:
// *   P, which may be of smaller size.
// *   ADJ, smaller size
// *
// *   Less important:
// *   VS/E of smaller size
// *
// *
// * @param P
// * @param mc
// */
//void graphpack_graph::triangle_core_pruning(vector<Vertex> & P, int & mc) {
//
//    /*
//     * Pruning using Triangle Cores
//     * TODO: Revise, to understand best option, should we recompute to remove edges,
//     *  or just leave them and pay the price
//     */
//    vector<long long> VS(vertices.size(),0);
//    vector<int> E;
//    E.reserve(edges.size());
//
//    cout << "[TRIANGLE CORE PRUNING]  |E| = " << num_edges() << " (before)" <<endl;
//    int start = 0;
//
//    /* Vertex triangle cores
//     * a vertex max triangle core is the maximum among all neighbors
//     */
//    int vertex_max_tcore = 0;
//    vector<int> vertex_tcore(num_vertices(),0);
//
//    /*
//     * IMPORTANT NOTE: tri_core should only be m/2,
//     * so we must use e_u and e_v to map edges, or eid!
//     * tri_core and max_tri_core are not set to be proper bounds, we must add 2 to each
//     */
//    print_break();
//    cout << "tri_core.size(): " << tri_core.size();
//    cout << ", |E| = " << num_edges() << ", edges.size(): " << edges.size() <<endl;
//    print_break();
//
//    // NOTE: When using H.eid, use j (edge_id) instead of edges[j] (node_id)
//    long long edge_pos = 0;
//    int valid_vertex = 0;
//    for (int i = 0; i < num_vertices(); i++) {
//        start = E.size();
//        for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
//
//            edge_pos = eid[j]; // tri_core is m/2, so this finds the right tri_core value given an edge id
//            if (tri_core[edge_pos]+2 > mc) { // get correct tri_core for edge's eid
//                E.push_back(edges[j]);
//
//
//                if ((tri_core[edge_pos]+2) > vertex_tcore[i])
//                    vertex_tcore[i] = tri_core[edge_pos]+2;
//            }
//            else { // remove edge
//                adj[i][edges[j]] = false;
//            }
//        }
//        if ((E.size() - start) > mc || vertex_tcore[i] > mc) { // basically the updated degree bound, but what about the edges???
//            P[valid_vertex].set_id(i);
//            valid_vertex++;
//         }
//        VS[i] = start;
//        VS[i + 1] = E.size();
//    }
//
//    for (int i = 0; i < num_vertices()-valid_vertex; i++)
//        P.pop_back();
//}
//
//
//
//
//void graphpack_graph::triangle_core_pruning_only_adj_update_degree(vector<Vertex> & P, int & mc) {
//    long long edge_pos = 0;
//    int v = 0, u = 0;
//    for (int e = 0; e < num_edges(); e++) { // tri_core is m/2
//
//        /** finds the right pos in e_u and e_v for that given edge id (from edges) */
//        edge_pos = eid[e];
//
//        /** get the vertex ids of each end of the edge */
//        v = e_v[edge_pos];
//        u = e_u[edge_pos];
//
//        // keep edges that may lead us to find a clique of _strictly a larger size_
//        if (tri_core[edge_pos]+2 <= mc) { // remove edge
//            adj[v][u] = false;
//            adj[u][v] = false;
//            degree[v]--; // if we updated these, we could then construct/update P very quickly by checking degree
//            degree[u]--;
//        }
//    }
//
//    int valid = 0;
//    for (int v = 0; v < num_vertices(); v++) {
//        if (degree[v] < mc) {
//            P[valid].set_id(v);
//            valid++;
//        }
//    }
//
//    // remove invalid verts
//    for (int v = 0; v < num_vertices() - valid; v++) {
//        P.pop_back();
//    }
//
//}
//
//
//
//
///*
// * Could update degree on the fly.
// * If we updated degree, then we could
// *  construct/update P very quickly by checking degree
// */
//void graphpack_graph::triangle_core_pruning_only_adj(int & mc) {
//    long long edge_pos = 0;
//    int v = 0, u = 0;
//    for (int e = 0; e < num_edges(); e++) { // tri_core is m/2
//
//        /** finds the right pos in e_u and e_v for that given edge id (from edges) */
//        edge_pos = eid[e];
//
//        /** get the vertex ids of each end of the edge */
//        v = e_v[edge_pos];
//        u = e_u[edge_pos];
//
//        // keep edges that may lead us to find a clique of _strictly a larger size_
//        if (tri_core[edge_pos]+2 < mc) { // remove edge
//            adj[v][u] = false;
//            adj[u][v] = false;
//        }
//    }
//}



void graphpack_graph::test_triangle_cores() {

    if (tri_core.size() == 0) {
        cout << "triangle cores have not been computed." <<endl;
    }
    else {
        set<int> tcore;
        for (int i = 0; i < num_vertices(); i++) {
            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {

                // tri_core is m/2, so this finds the right tri_core value given an edge id
                int edge_pos = eid[j];

                tcore.insert(tri_core[edge_pos+1]);
            }
        }
        cout << "finished finding the unique triangle core numbers, total = " << tcore.size() <<endl;
        set< int >::iterator it;
        for( it = tcore.begin(); it != tcore.end(); it++) {
            //        const int& core_num = (*it);
            cout << *it << " " <<endl;
        }
    }

}

















// SEE MCPACK_GRAPH_ECSC.cpp

//
///**
// * Attempt at parallelization and correcting the annoying +1 offset
// * June 12, 2013
// */
//
//
//void graphpack_graph::triangle_cores_in_parallel(bool is_parallel = true, int block_size = 64) {
//
//    map_edgeid_to_vertex_pair(is_parallel);
//
//    // number of triangles each edge participates
//    if (is_dense)
//        triangle_counting_in_parallel_adj();
//    else
//        triangle_counting_in_parallel();
//
//    // compute the triangle cores
//    triangle_core_numbers_parallel();
//
//    // compute bound
//    triangle_max_core();
//}
//
//
//
//
////void graphpack_graph::map_edgeid_to_vertex_pair() {
////    // edges (instead of vert ID store EID into e_u and e_v....
////    long long i, j, u, v;
////    long long m = edges.size();
////
////    eid.resize(m);
////    #pragma omp parallel for schedule(dynamic,512) private(i,j,u,v)
////    for (i=0; i < m/2; i++) {
////        v = e_v[i];
////        u = e_u[i];
////
////        for (j=vertices[v]; j<vertices[v+1]; j++) {
////            if (edges[j] == u) {
////                eid[j] = i;
////                break;
////            }
////        }
////
////        for (j=vertices[u]; j<vertices[u+1]; j++) {
////            if (edges[j] == v) {
////                eid[j] = i;
////                break;
////            }
////        }
////    }
////
////
////
////}
////
////
////
////void graphpack_graph::map_edgeid_to_vertex_pair(bool is_validate) {
////    // edges (instead of vert ID store EID into e_u and e_v....
////    long long m = edges.size();
////
////    printf("printing i=0...m/2-1,  e_v and e_u\n");
////    eid.resize(m);
////    for (int i=0; i < m/2; i++) {
////        int v = e_v[i];
////        int u = e_u[i];
////
////        printf("i=%d v=%d u=%d \n", i, v, u);
////
////    }
////
////
////    printf("\nv u | j eid e_v e_u\n");
////    for (int i=0; i < vertices.size()-1; i++) {
////
////        int v = i;
////
////        printf("N(%d)\n",v);
////        for (long long j=vertices[v]; j<vertices[v+1]; j++) {
////            int u = edges[j];
////            int edge_idx = j; // this indexes into eid array
////            int ptr_to_edge_in_emap = eid[edge_idx]; // this is a ptr to edge in the e_v and e_u arrays
////            int x = e_v[ptr_to_edge_in_emap];
////            int y = e_u[ptr_to_edge_in_emap];
////
////            printf("\tv=%d u=%d | edge_idx=%lld  eid[edge_idx]=%d  x=%d  y=%d \n",v,u,edge_idx,ptr_to_edge_in_emap,x,y);
////
////            validate((x==v && y==u) || (x==u && y==v),"edge-csc is invalid!");
////
////        }
////    }
////
////
////
////}
//
//
//
//
//// FIXED OFFSET
//void graphpack_graph::triangle_counting_in_parallel() {
//    long long e, v, u, w, i, j;
//    long long n = num_vertices();
//    long long m = num_edges(); // this is |E|/2
//
//    vector<long long> triangles(m,0);
//    vector<long long> ind(n, 0);
//
//    #pragma omp parallel for schedule(dynamic,512) \
//        shared(triangles) firstprivate(ind) private(v,u,w,e,i,j)
//    for (e=0; e<m; e++) {
//        v = e_v[e];   u = e_u[e];
////        if (u==v) { continue; }
//        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 1;
//        ind[v] = 1;
//
//        for (j=vertices[u]; j<vertices[u+1]; j++) {
//            w = edges[j];
//            if (w==u) { continue; }  //self-loop
//            if (w==v) { continue; }  //discard 1-hop self-loop
//            if (ind[w]) {
//                triangles[e] += 1;
//            }
//        }
//
//        for (i=vertices[v]; i<vertices[v+1]; i++) ind[edges[i]] = 0;
//        ind[v] = 0;
//    }
//
//    max_t_edge = 0;
//    total_t = 0;
//    for (e=0; e<m; e++) {
//        total_t += triangles[e];
//        if (triangles[e] > max_t_edge) {
//            max_t_edge = triangles[e];
//        }
//    }
//
//    tri_core = triangles;
//
////    tmp_tri_core.clear();
//    ind.clear();
//}
//
//
///*
// * @brief   ADJ data structure basically removes the cost
// *          of creating the fast index structure!
// *
// */
//void graphpack_graph::triangle_counting_in_parallel_adj() {
//    long long e, v, u, w, i, j;
//    long long n = num_vertices();
//    long long m = num_edges(); // this is |E|/2
//
//    vector<long long> triangles(m,0);
//
//    #pragma omp parallel for schedule(dynamic,512) \
//        shared(triangles) private(v,u,w,e,i,j)
//    for (e=0; e<m; ++e) { // v--u
//        v = e_v[e];   u = e_u[e];
////        if (u==v) { continue; }
//
//        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w
//            w = edges[j];
//            if (w==u) { continue; }  //self-loop
//            if (w==v) { continue; }  //discard 1-hop self-loop
//            if (adj[v][w]) {
//                triangles[e] += 1;
//            }
//        }
//    }
//
//    max_t_edge = 0;
//    total_t = 0;
//    for (e=0; e<m; e++) {
//        total_t += triangles[e];
//        if (triangles[e] > max_t_edge) {
//            max_t_edge = triangles[e];
//        }
//    }
//
//    tri_core = triangles;
////    tmp_tri_core.clear();
////    ind.clear();
//}
//
//
//void graphpack_graph::validate_triangle_cores(vector<long long> bin, vector<long long> tris, vector<long long> pos, int mt_end) {
//    int m = num_edges();
//    for (uint64_t e = 0; e < m; e++) {
//
//    }
//
//}
//
////void graphpack_graph::triangle_core_numbers_parallel() {
////    triangle_core_numbers_parallel(false);
////}
//
//void graphpack_graph::triangle_core_numbers_parallel() {
//    long long e, e_ux, e_x, tri, i, j, k, start, num, tmp;
//    long long eid_vu, eid_uw, eid_vw, v, u, w, tw, pw, px;
//    long long n = num_vertices(), m = num_edges();// + 1;
//    long long mt_end = max_t_edge+1;
//
//    vector <long long> bin(mt_end,0); // 0,...,max_t : bin[i] is the first edge with i triangles
//    vector <long long> tris(m);       // sort edges in increasing order of triangles: t[i] < t[j]
//    vector <long long> pos(m);        // position of the edges in tris
//
//    vector<long long> edge_neigh(2);  // exactly two neighboring edges of e form a triangle
//    vector<long long> proc(m,0);      // a triangle is processed, if any one of its edges are processed
//
//
//    for (e=0; e<m; e++) // bin edges by triangle counts
//        bin[tri_core[e]]++;
//
//
//    /** compute starting positions for the edges in each bin */
////    start = 1;
//    start = 0;
//    for (tri=0; tri < mt_end; tri++) {
//        num = bin[tri];
//        bin[tri] = start;
//        start = start + num;
//    }
//
//    /** bucket sort the edges by their triangle counts
//      *  place the edges into the correct buckets */
//    for (e=0; e<m; e++) {
//        pos[e] = bin[tri_core[e]];
//        tris[pos[e]] = e;
//        bin[tri_core[e]]++;
//    }
//
//    /** fixup the bin for the main computation */
//    for (tri=max_t_edge; tri > 0; tri--)  bin[tri] = bin[tri-1];
//    bin[0] = 0;
////    bin[0] = 1;
//
//    /**
//     * Process each edge in tris in the order of number of triangles
//     */
//    double sec = get_time();
//    vector<long long> ind(n, 0);
//    for (e=0; e<m; e++) {
//
//        if (e % 100000 == 0) { printf("%ld/%ld %lg\n", e,m, get_time()-sec); }
//        /** get the edge_id */
//        eid_vu = tris[e]; // @TODO: Could save the ordering..
//
//        /** get the vertex ids of each end of the edge */
//        v = e_v[eid_vu];
//        u = e_u[eid_vu];
//
//
//        /**
//         * Create lookup table for the neighbors of v.
//         * The storage cost is O(|E|).
//         * But instead of only setting the lookup value to 1,
//         * We set it to be the edge_id,
//         * or the position in the edge list, and therefore this indicates the exact neighbor.
//         * --------------
//         * For each neighbor of v denoted u \in N(v), we mark ind(u),
//         * but instead of simply marking it, we store the ptr to that
//         * position in the edges/neighbors array.
//         *
//         * This provides us with a way to check for a triangle in O(1) time, and also
//         * if a triangle does exist from that vertex, it allows us to grab its edgeid in O(1) time.
//         *
//         * The only problem is if an edge has a ptr of 0, then when we check for a neighbor, it will appear as the neighbor does not exist, while its actual position is 0.
//         * For this case, we simply add 1 to edge pos when marking that vertex, and after if we need it, then we subtract 1 from it to get the real edge position.
//         *
//         *
//         * We store the edge ptr of each vertex, and then use
//         * This works
//         */
//        for (i=vertices[v]; i<vertices[v+1]; i++) {
////            if (proc[eid[i]]) continue;
//                ind[edges[i]] = i+1;
//        }
////        ind[v] = v+1;
//
//        /**
//         * Get a neighbor of u, so that we have a wedge: v--u--w
//         * Now, we must check if w is also a neighbor of v, if so then we have a triangle.
//         */
//
//        for (j=vertices[u]; j<vertices[u+1]; j++) { // u--w for neigh of u
//            if (ind[edges[j]] && !proc[eid[j]] && !proc[eid[ind[edges[j]]-1]]) {//continue;
////            w = edges[j];
////            if (w==v) { continue; }  //self-loop
////            if (w==u) { continue; }  //discard 1-hop self-loop
//
//            /**
//             * v--w, if w is a neigh of v,
//             * so if true, then we know: v--u, u--w, w--v forming a triangle.
//             */
////            if (ind[w]) {
////                eid_uw = eid[j];        // get edge_id of u--w: the neighbor of u, so the edge is u--w
//                edge_neigh[0] = eid[j];
//                /**
//                 * We subtract one here, since we added one above when marking v's neighbors.
//                 * This was done precisely for the reason of getting the egdeid in O(1) time.
//                 * Since the edge_ids begin at 0, if we had not added 1,
//                 * then the vertex stored in the 0 position (edges[0]) would seem to not be connected to v.
//                 * Hence, we subtract 1, again to retrieve the correct edge position, while avoiding the above case.
//                 */
////                eid_vw = eid[ind[w]-1];   // get the edge_id of v--w
////                edge_neigh[1] = eid[ind[w]-1];
//                edge_neigh[1] = eid[ind[edges[j]]-1];
//
//                /**
//                 * We should be able to perform these checks as soon as know the edge_id, then we can check in proc[], and if marked, then we just skip.
//                 * This can also be done, when setting the ind, if it has already been processed, then skip it!
//                 */
//
////                if (!proc[edge_neigh[0]] && !proc[edge_neigh[1]]) {
////                    edge_neigh[0] = eid_uw;
////                    edge_neigh[1] = eid_vw;
//
//                    for (k=0; k<2; k++) { // each edge neighbors with two other edges to form a triangle
//                        if (tri_core[edge_neigh[k]] > tri_core[eid_vu]) {
//
//                            tw = tri_core[edge_neigh[k]];   //number of current tri for w
//                            pw = pos[edge_neigh[k]];        //actual pos of w
//                            px = bin[tw];                   //starting position of tw
//                            e_x = tris[px];                 //starting edge id
//
//                            if (edge_neigh[k] != e_x) {    //swap & update edge_neigh position
//                                pos[edge_neigh[k]] = px;   //starting position of e_
//                                tris[pw] = e_x;
//                                pos[e_x] = pw;
//                                tris[px] = edge_neigh[k];
//                            }
//                            bin[tw]++;   // update starting pos
//                            tri_core[edge_neigh[k]]--;
//                        }
//                    }
//                }
////            }
//        }
//
//        for (i=vertices[v]; i<vertices[v+1]; i++) // reset lookup table
//            ind[edges[i]] = 0;
////        ind[v] = 0;
//
//        proc[eid_vu] = 1; // set edge to be processed
//    }
//
//#ifndef DEBUG
//        cout << "max triangle core = " << tri_core[tris[m-1]] <<endl;
//#endif
//
////        max_tri_core = tri_core[tris[m-1]];
//
//        triangle_max_core(); // probably can remove this, since max tri core must exist in last pos
//
////        if (validate) {
////            validate_triangle_cores(bin,tris,pos,ind);
////        }
//
//        bin.clear();
//        pos.clear();
//        tris.clear();
//        ind.clear();
//}
//
//



//void graphpack_graph::print_tcore_info(vector<long long> bin, vector<long long> tris, vector<long long> pos, int mt_end) {
//    long long e, tri, m = num_edges();
//
//    cout << "\n\ttdeg = ";
//    for (e=1; e<m; e++)
//        cout << tri_core[e] << "\t";
//
//    cout << "\n\tbins = ";
//    for (tri=0; tri < mt_end; tri++)
//        cout << bin[tri] << "\t";
//
//    cout << "\n\ttris = ";
//    for (e=1; e<m; e++)
//        cout << tris[e] << "\t";
//
//    cout << "\n\tpos  = ";
//    for (e=1; e<m; e++)
//        cout << pos[e] << "\t";
//    cout << "\n" <<endl;
//}


void graphpack_graph::bin_triangle_core_numbers() {
    vector <int> bin(max_t_edge,0);
    int m = num_edges();

    for (uint64_t e = 0; e < m; e++) { // bin edges by triangle counts
        if (tri_core[e] > max_tri_core) {
            printf("ERROR: FOUND NUMBER LARGER THAN MAX\n\n");
        }
        else
            bin[tri_core[e]]++;
    }


//    string filename = string("trianglecore/triangle_core_dist/")+fn;
//    string filename = "trianglecore/triangle_core_dist/";
//    string fn_only = fn;
//    string fn_new = extract_filename(fn);
//    filename.append(fn_new);
//    filename.append(".txt");
//    printf("%s\n",filename.c_str());
//    FILE *fp = fopen(filename.c_str(), "w");
    printf("max triangle core = %lld\n",get_triangle_core_bound());
    for (int i = 0; i <= max_tri_core; ++i) {
        printf("T(%d) = %d \n",i,bin[i]);
//        fprintf(fp,"%d \t %d\n",i,bin[i]); /*writes*/
    }
//    fclose(fp); /*done!*/

    // count the unique triangle cores
    int unique_triangle_cores = 0;
    int triangle_core_gap = 0;
    for (int i = 1; i <= max_tri_core; ++i) {
        if (bin[i] == 0)
            triangle_core_gap++;
        else
            unique_triangle_cores++;
    }

    printf("unique triangle cores = (%d / %lld)\ntriangle core gap: %d\n",
            unique_triangle_cores,max_tri_core,triangle_core_gap);

    double tcore_ratio = (double)unique_triangle_cores/ (double) max_tri_core;
    printf("unique triangle cores ratio: %lg\n", tcore_ratio);


}





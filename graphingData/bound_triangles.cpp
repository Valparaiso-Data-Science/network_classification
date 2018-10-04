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

using namespace graphpack;


/**
 * Fastest Triangle Counter
 *
 * Returns an array of triangle counts for each vertex in G
 *
 * Computes the following:
 *      local clustering-coeffient
 *      local clustering-coeffient
 *      max number of triangles
 *      avg number of triangles
 *      total triangles in G
 *
 * Number of triangles in k-core
 *
 * @param block_size
 */
void graphpack_graph::compute_triangles_simple(int block_size) {
    long long v, u, w, i, j;
    long long n = num_vertices();
    double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;
    double num_local_wedges = 0;

    t.resize(n);
    vector<int> ind(n, 0);
    vector<int> tris(n, 0);

    int num_workers = omp_get_max_threads()+1;
    vector<long long> max_triangles(num_workers,0);
    double sum_local_cc = 0;
    long long sum_triangles = 0;
//    long long sum_triangles_unnorm = 0;
    long long sum_wedges = 0;

    printf("computing fast triangles *** \n");

#pragma omp parallel for schedule(dynamic,block_size) \
    shared(tris) firstprivate(ind) \
    private(v,u,w,d,i,j,num_local_wedges) \
    reduction(+:sum_local_cc,sum_triangles,sum_wedges)
    for (v = 0; v < n; ++v) {
        d = vertices[v+1] - vertices[v];

        for (i=vertices[v]; i<vertices[v+1]; ++i)
            ind[edges[i]] = v + 1;

        // v--u
        for (i=vertices[v]; i<vertices[v+1]; ++i) {
            u = edges[i];
            // u--w
            for (j=vertices[u]; j<vertices[u+1]; ++j) {
                w = edges[j];
                if (w==v) { continue; }
                if (ind[w]-1 == v)  tris[v] += 1.;
            }
        }
        // d must be greater than 1 for triangle to exist
        if (d > 1.)  {
            // number of triangles
            t[v] = tris[v]/2.;

            // total triangles
            sum_triangles += t[v];
//            sum_triangles_unnorm = tris[v];

            // number of wedges
            num_local_wedges = (d*(d-1.))/2.;
            sum_wedges += num_local_wedges;

            // local clustering coeff
            sum_local_cc += double(t[v] / num_local_wedges);
//            sum_local_cc += double(tris[v] / num_local_wedges);
        }
//        t[v] = tris[v]/2.;
//        sum_triangles += t[v];

        // update worker-local max triangles
        if (t[v] > max_triangles[omp_get_thread_num()])
            max_triangles[omp_get_thread_num()] = t[v];
    }


    // max triangles
    max_t = max_triangles[0];
    for (int p = 1; p < num_workers; ++p) {
        if (max_triangles[p] > max_t)
            max_t = max_triangles[p];
    }

    // local clustering-coeffient
    mean_cc = sum_local_cc / n;

    // global clustering coeffient
    global_cc = sum_triangles / (double)sum_wedges;

    // save total triangles ( 3 times the number of triangles)
    total_t = sum_triangles;

    // average triangles per vertex
    mean_t = sum_triangles / (1.0*n);

    printf("cc_avg: %lg, cc_global: %lg, T_max: %lld, T_avg: %lg, |T|=%lld \n",
            mean_cc, global_cc, max_t, mean_t, total_t);

    ind.clear();
    tris.clear();
}


/**
 * Fastest Triangle Counter
 *
 * Returns an array of triangle counts for each vertex in G
 *
 * Computes the following:
 *      local clustering-coeffient
 *      local clustering-coeffient
 *      max number of triangles
 *      avg number of triangles
 *      total triangles in G
 *
 * Number of triangles in k-core
 *
 * @param block_size
 */
void graphpack_graph::compute_triangles(int block_size) {
    long long v, u, w, i, j;
    long long n = num_vertices();
    double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;
    double num_local_wedges = 0;

    t.resize(n);
    vector<int> ind(n, 0);
    vector<int> tris(n, 0);

    int num_workers = omp_get_max_threads()+1;
    vector<long long> max_triangles(num_workers,0);
    double sum_local_cc = 0;
    long long sum_triangles = 0;
//    long long sum_triangles_unnorm = 0;
    long long sum_wedges = 0;

    printf("computing fast triangles *** \n");

#pragma omp parallel for schedule(dynamic,block_size) \
    shared(tris) firstprivate(ind) \
    private(v,u,w,d,i,j,num_local_wedges) \
    reduction(+:sum_local_cc,sum_triangles,sum_wedges)
    for (v = 0; v < n; ++v) {
        d = vertices[v+1] - vertices[v];

        for (i=vertices[v]; i<vertices[v+1]; ++i)
            ind[edges[i]] = v + 1;

        // v--u
        for (i=vertices[v]; i<vertices[v+1]; ++i) {
            u = edges[i];
//            if (u < v) continue;
            // u--w
            for (j=vertices[u]; j<vertices[u+1]; ++j) {
                w = edges[j];
//                if (w < v || w < u) continue;
                if (w==v) { continue; }
                if (ind[w]-1 == v) {
//                    tris[u] += 1;
//                    tris[w] += 1;

                    tris[v] += 1.;
                }
            }
        }
        // d must be greater than 1 for triangle to exist
        if (d > 1.)  {
            // number of triangles
            t[v] = tris[v]/2.;

            // total triangles
            sum_triangles += t[v];
//            sum_triangles_unnorm = tris[v];

            // number of wedges
            num_local_wedges = (d*(d-1.))/2.;
            sum_wedges += num_local_wedges;

            // local clustering coeff
            sum_local_cc += double(t[v] / num_local_wedges);
//            sum_local_cc += double(tris[v] / num_local_wedges);
        }
//        t[v] = tris[v]/2.;
//        sum_triangles += t[v];

        // update worker-local max triangles
        if (t[v] > max_triangles[omp_get_thread_num()])
            max_triangles[omp_get_thread_num()] = t[v];
    }


    // max triangles
    max_t = max_triangles[0];
    for (int p = 1; p < num_workers; ++p) {
        if (max_triangles[p] > max_t)
            max_t = max_triangles[p];
    }

    // local clustering-coeffient
    mean_cc = sum_local_cc / n;

    // global clustering coeffient
    global_cc = sum_triangles / (double)sum_wedges;

    // save total triangles ( 3 times the number of triangles)
    total_t = sum_triangles;

    // average triangles per vertex
    mean_t = sum_triangles / n;

    printf("cc_avg: %lg, cc_global: %lg, T_max: %lld, T_avg: %lg, |T|=%lld \n",
            mean_cc, global_cc, max_t, mean_t, total_t);

    ind.clear();
    tris.clear();
}

// BIGraph
void graphpack_graph::compute_weighted_triangles(int block_size) {
    int v, u, w;
    long long i, j, k;
    long long n = num_vertices();


//    weighted_triangles.resize(n,0);
//    weighted_wedges.resize(n,0);

    // fast lookup table for determining a neighbor in O(1) time
    vector<long long> ind(n, 0);

    // thread local/shared arrays for summing triangles/wedges counted by each worker-thread
//    vector<double> sum_weighted_triangles(omp_get_max_threads()+1,0);
//    vector<double> sum_weighted_wedges(omp_get_max_threads()+1,0);

    vector<double> tmp_weighted_triangles(n,0);
    vector<double> tmp_weighted_wedges(n,0);

//    printf("number of workers = %d \n", omp_num_threads()); //omp_get_num_threads());

//#pragma omp parallel for schedule(dynamic,block_size) \
//    shared(sum_weighted_triangles, sum_weighted_wedges) \
//    firstprivate(ind) \
//    private(v,u,w,i,j,k)


#pragma omp parallel for schedule(dynamic,block_size) \
    shared(tmp_weighted_triangles, tmp_weighted_wedges) \
    firstprivate(ind) \
    private(v,u,w,i,j,k)
    for (v = 0; v < n; ++v) {

//        if (verbose) printf("v = %d \n", v);

        // mark the neighbors of v, set neighbors to 1
        for (i = vertices[v]; i < vertices[v + 1]; ++i) {
            // index into neighbors position given by edges[i], and at that position, store the ptr
            ind[edges[i]] = i+1; // i is a ptr into edges array, REMEMBER to subtract 1..
        }

        // for each u \in N(v)
        for (i = vertices[v]; i < vertices[v+1]; ++i) {
            u = edges[i]; //get_neighbor(i); // get neighbor using the ptr i

            // for each w \in N(u), so now we know that v--u, and we marked v's neighbors in ind, so now w
            for (j = vertices[u]; j < vertices[u + 1]; ++j) {
                w = edges[j]; //get_neighbor(j);

//                if (w==u) { continue; } // check if w and u is a SELF-LOOP
                if (w==v) { continue; } // required

                // only need to get these weights if a triangle v--u--w exists
                double uv_weight = wt[i]; //get_weight(i);
                double uw_weight = wt[j]; //get_weight(j);

                // check if the neighbor w of u, IS actually a neighbor of v, using ind gives O(1) time check
                if (ind[w] > 0)  { // add weighted triangle

                    // ptr is stored in ind[w], we added 1 to avoid the case where the ptr = 0, so now we fix it by subtracting 1
                    k = ind[w] - 1;
                    double vw_weight = wt[k]; //get_weight(k);

                    // compute the weight of the triangle
                    double triangle_weight = (uv_weight * uw_weight * vw_weight);

//                    if (verbose)
//                        printf("\t v = %d, triangle = %d--%d--%d,  weight = %lg  (uv_wt = %lg, uw_wt = %lg, vw_wt = %lg) \n",
//                                v, u,w,v, triangle_weight, uv_weight, uw_weight, vw_weight);

                    tmp_weighted_triangles[v] = tmp_weighted_triangles[v] + (1. / triangle_weight);
                }
                else { // add weighted wedges
                    double wedge_weight = (uv_weight * uw_weight);
                    tmp_weighted_wedges[v] = tmp_weighted_wedges[v] + (1. / wedge_weight);
                }
            }
        }

        tmp_weighted_triangles[v] = tmp_weighted_triangles[v] / 2.;

        // each thread saves its wedges/triangles counted thus far
//        sum_weighted_wedges[omp_get_thread_num()] += tmp_weighted_wedges[v];
//        sum_weighted_triangles[omp_get_thread_num()] += tmp_weighted_triangles[v];

        for (i=vertices[v]; i<vertices[v+1]; ++i)
            ind[edges[i]] = 0;

    } // end parallel code

    // store triangles
//    weighted_triangles = tmp_weighted_triangles;
//    weighted_wedges = tmp_weighted_wedges;

    total_weighted_triangles = 0.0;
    total_weighted_wedges = 0.0;

    for (int v = 0; v < n; ++v) {
        total_weighted_triangles += tmp_weighted_triangles[v];
        total_weighted_wedges += tmp_weighted_wedges[v];
    }


//    // reduce computation
//    for (int tid = 1; tid < omp_get_num_threads(); ++tid) {
//        sum_weighted_triangles[0] += sum_weighted_triangles[tid];
//        sum_weighted_wedges[0] += sum_weighted_wedges[tid];
//    }
//    total_weighted_triangles = sum_weighted_triangles[0];
//    total_weighted_wedges = sum_weighted_wedges[0];

    // clear tmp vars
    ind.clear();
//    tmp_weighted_triangles.clear();
//    tmp_weighted_wedges.clear();
}

// BIGraph
void graphpack_graph::compute_weighted_triangles_reduce(int block_size) {
    int v, u, w;
    long long i, j, k;
    long long n = num_vertices();


//    weighted_triangles.resize(n,0);
//    weighted_wedges.resize(n,0);

    // fast lookup table for determining a neighbor in O(1) time
    vector<long long> ind(n, 0);

    // thread local/shared arrays for summing triangles/wedges counted by each worker-thread
    vector<double> sum_weighted_triangles(omp_get_max_threads()+1,0);
    vector<double> sum_weighted_wedges(omp_get_max_threads()+1,0);

    vector<double> tmp_weighted_triangles(n,0);
    vector<double> tmp_weighted_wedges(n,0);

//    printf("number of workers = %d \n", omp_num_threads()); //omp_get_num_threads());


#pragma omp parallel for schedule(dynamic,block_size) \
    shared(tmp_weighted_triangles, tmp_weighted_wedges) \
    firstprivate(ind) \
    private(v,u,w,i,j,k)
    for (v = 0; v < n; ++v) {

//        if (verbose) printf("v = %d \n", v);

        // mark the neighbors of v, set neighbors to 1
        for (i = vertices[v]; i < vertices[v + 1]; ++i) {
            // index into neighbors position given by edges[i], and at that position, store the ptr
            ind[edges[i]] = i+1; // i is a ptr into edges array, REMEMBER to subtract 1..
        }

        // for each u \in N(v)
        for (i = vertices[v]; i < vertices[v+1]; ++i) {
            u = edges[i]; //get_neighbor(i); // get neighbor using the ptr i

            // for each w \in N(u), so now we know that v--u, and we marked v's neighbors in ind, so now w
            for (j = vertices[u]; j < vertices[u + 1]; ++j) {
                w = edges[j]; //get_neighbor(j);

//                if (w==u) { continue; } // check if w and u is a SELF-LOOP
                if (w==v) { continue; } // required

                // only need to get these weights if a triangle v--u--w exists
                double uv_weight = wt[i]; //get_weight(i);
                double uw_weight = wt[j]; //get_weight(j);

                // check if the neighbor w of u, IS actually a neighbor of v, using ind gives O(1) time check
                if (ind[w] > 0)  { // add weighted triangle

                    // ptr is stored in ind[w], we added 1 to avoid the case where the ptr = 0, so now we fix it by subtracting 1
                    k = ind[w] - 1;
                    double vw_weight = wt[k]; //get_weight(k);

                    // compute the weight of the triangle
                    double triangle_weight = (uv_weight * uw_weight * vw_weight);

//                    if (verbose)
//                        printf("\t v = %d, triangle = %d--%d--%d,  weight = %lg  (uv_wt = %lg, uw_wt = %lg, vw_wt = %lg) \n",
//                                v, u,w,v, triangle_weight, uv_weight, uw_weight, vw_weight);

                    tmp_weighted_triangles[v] = tmp_weighted_triangles[v] + (1. / triangle_weight);
                }
                else { // add weighted wedges
                    double wedge_weight = (uv_weight * uw_weight);
                    tmp_weighted_wedges[v] = tmp_weighted_wedges[v] + (1. / wedge_weight);
                }
            }
        }

        tmp_weighted_triangles[v] = tmp_weighted_triangles[v] / 2.;

        // each thread saves its wedges/triangles counted thus far
        sum_weighted_wedges[omp_get_thread_num()] += tmp_weighted_wedges[v];
        sum_weighted_triangles[omp_get_thread_num()] += tmp_weighted_triangles[v];

        for (i=vertices[v]; i<vertices[v+1]; ++i)
            ind[edges[i]] = 0;

    } // end parallel code

    // store triangles
//    weighted_triangles = tmp_weighted_triangles;
//    weighted_wedges = tmp_weighted_wedges;

    total_weighted_triangles = 0.0;
    total_weighted_wedges = 0.0;

//    // reduce computation
    for (int tid = 1; tid < omp_get_max_threads(); ++tid) {
        sum_weighted_triangles[0] += sum_weighted_triangles[tid];
        sum_weighted_wedges[0] += sum_weighted_wedges[tid];
    }
    total_weighted_triangles = sum_weighted_triangles[0];
    total_weighted_wedges = sum_weighted_wedges[0];

    // clear tmp vars
    ind.clear();
//    tmp_weighted_triangles.clear();
//    tmp_weighted_wedges.clear();
}

/**
 * Parallel Vertex Triangle Counting
 *
 * Handles both adj and csc representations
 *
 * @param max
 * @param is_parallel
 */
void graphpack_graph::compute_vertex_triangles(int* &pruned, int max, bool is_parallel, bool is_enum) {
    long long v, u, w, i, j;
    long long n = num_vertices();
    double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;

    t.resize(n);
//    kappa.resize(n);

    total_t = 0;
    vector<int> ind(n, 0);

    vector<int> tris(n, 0);
//    vector<int> tris_max(n, 0);


    if (is_dense_graph()) {
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
                private(v,u,w,d,i,j)
            for (v = 0; v < n; v++) {
                if (pruned[v]) continue;

                d = vertices[v+1] - vertices[v];

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
                    if (pruned[w]) continue;
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
                        if (pruned[u]) continue;
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v]+2 < max) {
                    t[v] = 0;
                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
        else {
            for (v = 0; v < n; v++) {
                if (pruned[v]) continue;

//                d = vertices[v+1] - vertices[v];

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
                    if (pruned[w]) continue;
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
                        if (pruned[u]) continue;
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v]+2 < max) {
                    t[v] = 0;
                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
    }
    else { // csc only
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
                private(v,u,w,d,i,j)
            for (v = 0; v < n; v++) {
                if (pruned[v]) continue;

//                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v + 1;


                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
                    if (pruned[w]) continue;
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
                        if (pruned[u]) continue;
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v]+2 < max) {
                    t[v] = 0;
                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            }
        }
        else {
            for (v = 0; v < n; v++) {
                if (pruned[v]) continue;

                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v + 1;


                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
                    if (pruned[w]) continue;
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
                        if (pruned[u]) continue;
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);

//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v]+2 < max) {
                    t[v] = 0;
                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;

            }
        }
    }




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

/**
 * Parallel Vertex Triangle Counting
 *
 * Handles both adj and csc representations
 *
 * @param max
 * @param is_parallel
 */
void graphpack_graph::compute_vertex_triangles(vector<Vertex> &P, int* &pruned, int max, bool is_parallel) {
    long long v, u, w, i, j, k;
    long long n = num_vertices();
//    long long n = P.size();
    double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;

    t.resize(n);
//    kappa.resize(n);

//    total_t = 0;
    vector<int> ind(n, 0);

    vector<int> tris(n, 0);
//    vector<int> tris_max(n, 0);
//    vector<int> tris_pruned(n, 0);

    if (is_dense_graph()) {
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
                private(v,u,w,d,i,j)
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[vertex_lookup[v]]) continue;

//                d = vertices[v+1] - vertices[v];

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
//                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
        else {
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[vertex_lookup[v]]) continue;

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
//                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
    }
    else { // csc only
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
                private(v,u,w,d,i,j)
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[vertex_lookup[v]]) continue;

//                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v+1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
//                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            }
        }
        else {
            for (k = 0; k < P.size(); k++) {

                v = P[k].get_id();

                if (pruned[vertex_lookup[v]]) continue;

//                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v + 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;
                P[k].set_bound(t[v]);

                if (t[v]+2 < max) {
                    t[v] = 0;
                    P[k].set_bound(0);
//                    pruned[v] = 1;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            }
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


/**
 * Parallel Vertex Triangle Counting
 *
 * Handles both adj and csc representations
 *
 * @param max
 * @param is_parallel
 */
void graphpack_graph::compute_vertex_triangles(int max, bool is_parallel) {
    long long v, u, w, i, j;
    long long n = num_vertices();
    double d = 0, tri = 0, t_sum = 0, t_pot = 0, local_max = 0;

    t.resize(n);
//    kappa.resize(n);

    total_t = 0;
    vector<int> ind(n, 0);

    vector<int> tris(n, 0);
//    vector<int> tris_max(n, 0);


    if (is_dense_graph()) {
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
	            private(v,u,w,d,i,j)
            for (v = 0; v < n; v++) {
                d = vertices[v+1] - vertices[v];

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v] < max) {
                    t[v] = 0;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
        else {
            for (v = 0; v < n; v++) {
//                d = vertices[v+1] - vertices[v];

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 1;
//                ind[v] = 1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (adj[u][v])  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v] < max) {
                    t[v] = 0;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
//                ind[v] = 0;
            }
        }
    }
    else { // csc only
        if (is_parallel) {
            #pragma omp parallel for schedule(dynamic,64) \
                shared(tris) firstprivate(ind) \
	            private(v,u,w,d,i,j)
            for (v = 0; v < n; v++) {
//                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v + 1;


                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v] < max) {
                    t[v] = 0;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            }
        }
        else {
            for (v = 0; v < n; v++) {
//                d = vertices[v+1] - vertices[v];

                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = v+1;

                for (i=vertices[v]; i<vertices[v+1]; i++) {
                    w = edges[i];
//                    if (v==w) { d-=1.; continue; }
                    for (j=vertices[w]; j<vertices[w+1]; j++) {
                        u = edges[j];
//                        if (u==w) { continue; }
                        if (u==v) { continue; }
                        if (ind[u]-1 == v)  tris[v] += 1.;
                    }
                }
//                if (d > 1.)  {
//                    tris_max[v] = d*(d-1.);
                    //            kappa[v] = double(tris[v] / tris_max[v]);
//                }
                //        else
                //            kappa[v] = 0.;
                t[v] = tris[v]/2.;

                if (t[v] < max) {
                    t[v] = 0;
                }

//                for (i=vertices[v]; i<vertices[v+1]; i++)  ind[edges[i]] = 0;
            }
        }
    }




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


void graphpack_graph::triangle_stats() {
	long long v;
	long long n = num_vertices();

	max_t = 0;
	mean_cc = 0;
	mean_t = 0;
	tri_ub = 0;
	for (v = 0; v < n; v++) {
		mean_cc += kappa[v];
		mean_t += t[v];
		if (t[v] > max_t)
			max_t = t[v];
	}
	tri_ub = sqrt(2*max_t) + 1;
	mean_t = mean_t / (1.0*n);
	mean_cc = mean_cc / (1.0*n);
}

void graphpack_graph::triangle_bound() {
	long long v;
	long long n = vertices.size() - 1;

	t_ub.resize(n);
	for (v = 0; v < n; v++) {
		t_ub[v] = sqrt(2*t[v]);
	}
}

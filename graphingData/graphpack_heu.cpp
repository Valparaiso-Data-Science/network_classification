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

#include "graphpack_heu.h"

using namespace graphpack;
using namespace std;

void graphpack_heu::branch(vector<Vertex>& P, int sz,
        int& mc, vector<int>& C, vector<int>& ind) {

    if (P.size() > 0) {

        int u = P.back().get_id();
        P.pop_back();

        for (long long j = (*V)[u]; j < (*V)[u + 1]; j++)  ind[(*E)[j]] = u;

        vector <Vertex> R;
        R.reserve(P.size());
        for (int i = 0; i < P.size(); i++)
            if (ind[P[i].get_id()] == u)
                if ((*K)[P[i].get_id()] > mc)
                    R.push_back(P[i]);

        for (long long j = (*V)[u]; j < (*V)[u + 1]; j++)  ind[(*E)[j]] = 0;

        int mc_prev = mc;
        branch(R, sz + 1, mc, C, ind);

        if (mc > mc_prev)  C.push_back(u);

        R.clear();  P.clear();
    }
    else if (sz > mc)
        mc = sz;
    return;
}

int graphpack_heu::search_bounds(graphpack_graph& G,
        vector<int>& C_max) {

    V = G.get_vertices();
    E = G.get_edges();
    degree = G.get_degree();
    vector <int> C, X;
    C.reserve(ub);
    C_max.reserve(ub);
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);
    vector<int> ind(G.num_vertices(),0);

    bool found_ub = false;
    int mc = 0, mc_prev, mc_cur, i, v, k, lb_idx = 0;


    printf("heuristic optimized for global search over k-cores\n");

    int start = G.num_vertices()-1;
    int end = 0;
    if (global_small_to_large) {
        start = 0;
        end = G.num_vertices()-1;
    }

    if (local_small_to_large) {

#pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, X, mc, C_max, lb_idx) private(i, v, P, mc_prev, mc_cur, C, k) firstprivate(ind)
        for (i = start; i >= end; --i) {
            if (found_ub) continue;

            v = (*order)[i];
            mc_prev = mc_cur = mc;

            if ((*K)[v] > mc) {
                for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)
                    if ((*K)[(*E)[j]] > mc)
                        P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                if (P.size() > mc_cur) {
                    std::sort(P.begin(), P.end(), incr_heur);
                    branch(P, 1 , mc_cur, C, ind);

                    if (mc_cur > mc_prev) {
                        if (mc < mc_cur) {
                            #pragma omp critical
                            if (mc < mc_cur) {
                                mc = mc_cur;
                                C.push_back(v);
                                C_max = C;
                                if (mc >= ub) found_ub = true;
                                if (G.verbose)  print_info(C_max);
                            }
                        }
                    }
                }
                C = X; P = T;
            }
        }
    }
    else {
#pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, X, mc, C_max, lb_idx) private(i, v, P, mc_prev, mc_cur, C, k) firstprivate(ind)
        for (i = start; i >= end; --i) {
            if (found_ub) continue;

            v = (*order)[i];
            mc_prev = mc_cur = mc;

            if ((*K)[v] > mc) {
                for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)
                    if ((*K)[(*E)[j]] > mc)
                        P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                if (P.size() > mc_cur) {
                    std::sort(P.begin(), P.end(), desc_heur);
                    branch(P, 1 , mc_cur, C, ind);

                    if (mc_cur > mc_prev) {
                        if (mc < mc_cur) {
                            #pragma omp critical
                            if (mc < mc_cur) {
                                mc = mc_cur;
                                C.push_back(v);
                                C_max = C;
                                if (mc >= ub) found_ub = true;
                                if (G.verbose)  print_info(C_max);
                            }
                        }
                    }
                }
                C = X; P = T;
            }
        }
    }
    if (G.verbose)
        cout << "[mcpack heuristic]\t mc = " << mc <<endl;
    return mc;
}

int graphpack_heu::search_bounds_ordering(graphpack_graph& G, vector<int>& C_max, params & p) {

    V = G.get_vertices();
    E = G.get_edges();
    degree = G.get_degree();
    vector <int> C, X;
    C.reserve(ub);
    C_max.reserve(ub);
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);
    vector<int> ind(G.num_vertices(),0);

    bool found_ub = false;
    int mc = 0, mc_prev, mc_cur, i, v, k, lb_idx = 0;

    // order verts for our search routine
    vector<Vertex> VS;
    VS.reserve(G.num_vertices());
    compute_ordering(VS,G,global_ordering,p.heu_global_s2l,0,0);
    cout << "|V| = " << VS.size() <<endl;



    if (local_small_to_large) {

            #pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, X, mc, C_max, lb_idx) private(i, v, P, mc_prev, mc_cur, C, k) firstprivate(ind)
        //    for (i = G.num_vertices()-1; i >= 0; --i) {
        for (i = 0; i < G.num_vertices(); ++i) {
            if (found_ub) continue;

            //        v = (*order)[i];

            mc_prev = mc_cur = mc;

            v = VS[i].get_id();
            if ((*K)[v] > mc) {

                for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)
                    if ((*K)[(*E)[j]] > mc)
                        P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                if (P.size() > mc_cur) {
                    std::sort(P.begin(), P.end(), incr_heur);
                    branch(P, 1 , mc_cur, C, ind);

                    if (mc_cur > mc_prev) {
                        if (mc < mc_cur) {
                            #pragma omp critical
                            if (mc < mc_cur) {
                                mc = mc_cur;
                                C.push_back(v);
                                C_max = C;
                                if (mc >= ub) found_ub = true;
                                if (G.verbose)  print_info(C_max);
                            }
                        }
                    }
                }
                C = X; P = T;
            }
        }
    }
    else {
        #pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, X, mc, C_max, lb_idx) private(i, v, P, mc_prev, mc_cur, C, k) firstprivate(ind)
        //    for (i = G.num_vertices()-1; i >= 0; --i) {
        for (i = 0; i < G.num_vertices(); ++i) {
            if (found_ub) continue;

            //        v = (*order)[i];

            mc_prev = mc_cur = mc;

            v = VS[i].get_id();
            if ((*K)[v] > mc) {

                for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)
                    if ((*K)[(*E)[j]] > mc)
                        P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                if (P.size() > mc_cur) {
                    std::sort(P.begin(), P.end(), desc_heur);
                    branch(P, 1 , mc_cur, C, ind);

                    if (mc_cur > mc_prev) {
                        if (mc < mc_cur) {
                            #pragma omp critical
                            if (mc < mc_cur) {
                                mc = mc_cur;
                                C.push_back(v);
                                C_max = C;
                                if (mc >= ub) found_ub = true;
                                if (G.verbose)  print_info(C_max);
                            }
                        }
                    }
                }
                C = X; P = T;
            }
        }
    }
    if (G.verbose)
        cout << "[mcpack heuristic]\t mc = " << mc <<endl;
    return mc;
}

int graphpack_heu::clique_compression(
        graphpack_graph& G,
        vector<int>& C_max,
        vector<int> & vertex_order_by_color,
        set< vector<int> >& cliques,
        params & p) {

//    V = G.get_vertices();
//    E = G.get_edges();
//    degree = G.get_degree();
    vector <int> X;
    X.reserve(G.get_max_degree()+1);
    C_max.reserve(G.get_max_degree()+1);
//    vector<Vertex> P, T;
//    P.reserve(G.get_max_degree()+1);
//    T.reserve(G.get_max_degree()+1);
    vector<int> ind(G.num_vertices(),0);

    bool found_ub = false;
    int max = 0; //, mc_prev, mc_cur,
    int i, v, k, lb_idx = 0;

    int* pruned = new int[G.num_vertices()];
    memset(pruned, 0, G.num_vertices() * sizeof(int));



    // assumes the vertices are sorted from smallest color index to largest, otherwise pass in number of colors
    int num_colors = G.color[vertex_order_by_color[vertex_order_by_color.size()-1]] + 1; // add 1 since colors start at 1
    printf("number of colors used in distance-two coloring is %d \n", num_colors-1);



    int n = G.num_vertices();


    // order verts for our search routine
    vector<Vertex> VS;
    VS.reserve(G.num_vertices());
    compute_ordering(VS,G,global_ordering,p.heu_global_s2l,0,0);

//    print_line(80);
    printf("vertices ordered globally by first ordering S_i, then the vertices in each color-class are ordered using kcore_sum \n");
    cout << "|V| = " << VS.size() << ", f(v_0)=" << VS[0].get_bound() << ",  f(v_n)=" << VS[VS.size()-1].get_bound() <<endl;

    vector< set< vector<int> > > cliques_local(p.threads);

    int cutoff = 1;
    if (n > 20)
        cutoff = 3;

    printf("cutoff = %d, number of threads = %d \n", cutoff, p.threads);

    int search_limit = 1000; // per iteration!
    int iter = 0;
    int num_removed = 0;
    while (VS.size() > 0 && num_removed < n) {

        //        for (i = 0; i < limit; ++i) {
        //
        //            v = W[i].get_id();
        //        }


        max = cutoff;

        C_max.clear();
        C_max = X;
        //        C_max.reserve(G.get_max_degree()+1);

        int limit = (VS.size()-1) - search_limit;
        if (limit < 0) limit = -1;

        #pragma omp parallel for schedule(dynamic,block_size) \
            shared(G, C_max, max, pruned) \
            private(v) \
            firstprivate(ind)
        for (i = (VS.size()-1); i > limit; --i) {
            //        for (i = G.num_vertices()-1; i >= 0; --i) {
            //        for (i = 0; i < G.num_vertices(); ++i) {

            v = VS[i].get_id();

            //            if (!pruned[])
            //            if (G.degree[v] <= cutoff) {
            //                pruned[v] = 1;
            //                continue;
            //            }



            if (G.kcore[v] > max) {

                vector<Vertex> P;
                P.reserve(G.get_max_degree()+1);

                for (long long j = G.vertices[v]; j < G.vertices[v + 1]; ++j)
                    if (!pruned[G.edges[j]])
                        if (G.kcore[G.edges[j]] > max)
                            P.push_back( Vertex(G.edges[j], compute_heuristic(G.edges[j],G)) );


                if (P.size() + 1 > max) {

                    vector<int> C;
                    C.reserve(G.get_max_degree()+1);
                    C.push_back(v);



                    if (local_small_to_large)
                        std::sort(P.begin(), P.end(), incr_heur);
                    else
                        std::sort(P.begin(), P.end(), desc_heur);
                    //                    G.vertex_bucket_sort(P,local_small_to_large);

                    // find largest clique at each vertex
                    while (P.size() + C.size() > max) {
                        int u = P.back().get_id();
                        P.pop_back();

                        C.push_back(u);


                        for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j)
                            ind[G.edges[j]] = v+1;

                        vector <Vertex> R;
                        R.reserve(P.size());

                        for (int w = 0; w < P.size(); ++w)
                            if (ind[P[w].get_id()] == v+1)
                                if (!pruned[P[w].get_id()])
                                    // no pruning, since we want to return the largest clique each iteration.
                                    if (G.kcore[P[i].get_id()] > max)
                                        R.push_back(P[w]);

                        for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j)
                            ind[G.edges[j]] = 0;

                        //                        if (R.size()) {
                        //                            //                    printf("[ mcpack worker %d ] BEFORE:  |P| = %lu, |R| = %lu \n",
                        //                            //                            omp_get_thread_num(), P.size(), R.size());
                        //                            P = R;
                        //                            //                    printf("[ mcpack worker %d ] AFTER:   |P| = %lu, |R| = %lu \n",
                        //                            //                            omp_get_thread_num(), P.size(), R.size());
                        //                        }
                        //
                        if (R.size() > 0) {
                            P = R;
                        }
                        if (C.size() > max) {
                            // obtain lock
#pragma omp critical (update_mc)
                            if (C.size() > max) {
                                printf("\t found larger clique of size = %d \n", max);

                                // ensure updated max is flushed
                                max = C.size();
                                C_max = C;
                                //                                print_mc_info(C,sec);
                            }
                        }
                        // backtrack and search another branch
                        R.clear();



                    } // end while

                } // end if B(v) > max
                //        pruned[v] = 1;
            }

        } // parallel for each vertex in order



        // sort clique C, then add to set
//        if (n <= 20)
//            printf("\t [mcpack worker %d]  searched v=%d, found clique of size %lu \n",
//                    omp_get_thread_num(), v, C.size());
//
//        if (C_max.size() == 0) {
//            printf("terminating early since |C_max| = 0, hence, no cliques of a larger size");
//            break;
//        }


        std::sort(C_max.begin(), C_max.end());
        cliques.insert(C_max);
//        cliques_local[omp_get_thread_num()].insert(C);

        // implicitly prune each vertex in C from the graph
        for (int h = 0; h < C_max.size(); ++h) {
            pruned[C_max[h]] = 1;
        }

//        print_mc_info(C_max,sec);

        // for each vertex in C, decrease the degree of each neighbor by 1
//        for (int h = 0; h < C_max.size(); ++h) {
//            int u = C_max[h];
//            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j) {
//                G.degree[G.edges[j]]--;
//            }
//        }

//        int last_vertex_idx = VS.size() - 1;

        int num_vertices = VS.size();
        if (num_vertices < search_limit) search_limit = num_vertices;

        for (int j = (num_vertices-1); j > ((num_vertices-1) - search_limit); --j) {
            int w = VS[j].get_id();
            if (pruned[w]) {

                VS[j].set_id(VS[VS.size()-1].get_id());
                VS[j].set_bound(VS[VS.size()-1].get_bound());
                VS.pop_back();
//                last_vertex_idx = VS.size()-1;

            }
        }

        num_removed += C_max.size();
        printf("iter=%d, |C| = %lu,  |V|=%lu  (%d),   ", iter, C_max.size(), VS.size(),G.num_vertices());
        for (int h = 0; h < C_max.size(); ++h) {
            printf("%d ", C_max[h]);
        }
        printf("\n");
        iter++;


    } // while more vertices

    printf("finished \n\n");



//    print_line(80);
//    double s = tic();
//
//    cliques = cliques_local[0];
//    for (int t = 1; t < p.threads; ++t) {
//
//        set< vector<int> >::iterator it;
//        for( it = cliques_local[t].begin(); it != cliques_local[t].end(); it++) {
//            const vector<int>& clq = (*it);
//
//            // find max clique found overall, useful to report
//            if (clq.size() > max)
//                max = clq.size();
//
//            cliques.insert(clq);
//        }
//
//        printf("[ mcpack worker %d ]  found %lu cliques of size %d, total unique cliques = %lu \n",
//                t, cliques_local[t].size(), max, cliques.size());
//    }
//    toc(s);
//    printf("took %lg seconds to find unique cliques \n\n", s);
//
//
//    if (cliques.size() < 10) {
//
//        set< vector<int> >::iterator it;
//        for( it = cliques.begin(); it != cliques.end(); it++) {
//            const vector<int>& clq = (*it);
//            printf("|C| = %lu,  C = {", clq.size());
//            for (int i = 0; i < clq.size(); ++i) {
//                if (i == clq.size()-1)  printf("%d} \n", clq[i]);
//                else    printf("%d, ", clq[i]);
//            }
//        }
//
//    }
//
//    print_line(80);

    if (pruned) delete[] pruned;
    return max;
}

int graphpack_heu::clique_compression_distance_two(
        graphpack_graph& G,
        vector<int>& C_max,
        vector<int> & vertex_order_by_color,
        set< vector<int> >& cliques,
        params & p) {

//    V = G.get_vertices();
//    E = G.get_edges();
//    degree = G.get_degree();
//    vector <int> C, X;
//    C.reserve(ub);
    C_max.reserve(G.get_max_degree()+1);
//    vector<Vertex> P, T;
//    P.reserve(G.get_max_degree()+1);
//    T.reserve(G.get_max_degree()+1);
    vector<int> ind(G.num_vertices(),0);

    bool found_ub = false;
    int max = 0; //, mc_prev, mc_cur,
    int i, v, k, lb_idx = 0;

    int* pruned = new int[G.num_vertices()];
    memset(pruned, 0, G.num_vertices() * sizeof(int));

    // assumes the vertices are sorted from smallest color index to largest, otherwise pass in number of colors
    int num_colors = G.color[vertex_order_by_color[vertex_order_by_color.size()-1]] + 1; // add 1 since colors start at 1
    printf("number of colors used in distance-two coloring is %d \n", num_colors-1);



    int n = vertex_order_by_color.size();

    vector<int> bin(num_colors,0);
    vector<int> kcore_sum_array(num_colors,0);
    vector<Vertex> kcore_sum;
    kcore_sum.reserve(num_colors);
    for (int i = 0; i < num_colors; ++i) {
        kcore_sum.push_back(Vertex(i,0));
    }


    // sorting the vertices now
//    vector<Vertex> VS;
//    VS.reserve(G.num_vertices());


    for (int i = 0; i < vertex_order_by_color.size(); ++i) {
        int v = vertex_order_by_color[i];
        int vertex_color = G.color[v];

        bin[vertex_color]++; // count vertices assigned to the kth color
        int sum = kcore_sum[vertex_color].get_bound() + G.kcore[v];
        kcore_sum[vertex_color].set_bound(sum);
        kcore_sum_array[vertex_color] = sum;


    }

    // sort color classes by kcore_sum
    // incr_heur = smallest to largest, desc_heur = largest to smallest
    std::sort(kcore_sum.begin(), kcore_sum.end(), decr_bound);

    // use old color as index, returns new color idx (from ordering)
    vector<int> color_to_newcolor(num_colors);
    int max_color_class = 0;
    print_line(80);
    for (int i = 0; i < num_colors-1; ++i) {
        int color = kcore_sum[i].get_id();

        color_to_newcolor[color] = i;

        if (bin[color] > max_color_class) {
            max_color_class = bin[color];
        }

        printf("\t color = %d, kcore_sum = %d, size of color class = %d, normalized score = %lg \n",
                color, kcore_sum[i].get_bound(), bin[color],
                (kcore_sum[i].get_bound()/(double)bin[color]));
    }
    print_line(80);
    printf("max size of distance-two color class = %d \n", max_color_class);


    // init the neigh coloring array
    vector< vector<Vertex> > color_classes(num_colors);
    for (int i = 0; i < num_colors; i++)  {
        // reserve space for each color class
        color_classes[i].reserve(max_color_class+1);
    }




    // color classes are ordered using kcore_sum
    // now add vertices to each color, not using the ordering

    // add vertices to their color classes (ordered)
    for (int i = 0; i < G.num_vertices(); ++i) {
        int v = vertex_order_by_color[i];
        int old_color_idx = G.color[v];
        int new_color_idx = color_to_newcolor[old_color_idx];

        color_classes[new_color_idx].push_back(Vertex(v,G.kcore[v]));
    }

    // sort vertices in each color class by kcore or any other function
    for (int c = 0; c < num_colors; ++c) {
        std::sort(color_classes[c].begin(), color_classes[c].end(), decr_bound);
    }

    // now copy



    vector<Vertex> VS;
    VS.reserve(G.num_vertices());
    for (int c = 0; c < num_colors; ++c) {
        VS.insert(VS.end(), color_classes[c].begin(), color_classes[c].end());
//        std::copy(color_classes[c].begin(), color_classes[c].end(), VS.end());
        printf("|V|=%lu, ", VS.size());
    }
    printf("\n");

    print_line(80);
    for (int i = 0; i < G.num_vertices(); ++i) {
        if (n <= 20) printf("[graphpack: vertices sorted by S_i and color-class]  i=%d, v = %d, K(v) = %d \n", i, VS[i].get_id(), VS[i].get_bound());
    }
    print_line(80);

//
//    int n = vertex_order_by_color.size();
//    if (n <= 20) print_line(80);
//    for (int i = 0; i < n; ++i) {
////        int v = vertex_order_by_color[i];
//        int color = G.color[i];
//
//        if (n <= 20)
//            printf("[mcpack adding vertex to set]  v = %d, color = %d, kcore sum = %d \n",
//                i, color, kcore_sum_array[color]);
//
//        VS.push_back(Vertex(i,kcore_sum_array[color]));
//    }
//    if (n <= 20)    print_line(80);
//
//    std::sort(VS.begin(), VS.end(), desc_heur);
//    for (int i = 0; i < n; ++i) {
//        int v = VS[i].get_id();
//        int sum_kcore = VS[i].get_bound();
//        if (n <= 20) printf("[mcpack after sorting]  v = %d, S_%d, kcore_sum = %d \n", v, G.color[v], sum_kcore);
//    }
//    if (n <= 20) print_line(80);
//
    // use -1 since last color is 0, using the above ordering
//    for (int i = 0; i < num_colors-1; ++i) {
//        int color = kcore_sum[i].get_id();
//
//
//
//    }

//    int start = 0;
//    for (int d = 0; d < num_colors; ++d) {
//        int num = bin[d];
//        bin[d] = start;
//        start = start + num;
//    }
//
//    int n = vertex_order_by_color.size();
//    for (v = 0; v < n; ++v) {
//        pos[v] = bin[bound[v]];
//        order[pos[v]] = v;
//        bin[bound[v]]++;
//    }



    // order verts for our search routine
//    vector<Vertex> VS;
//    VS.reserve(G.num_vertices());
//    compute_ordering(VS,G,global_ordering,p.heu_global_s2l,0,0);

//    print_line(80);
    printf("vertices ordered globally by first ordering S_i, then the vertices in each color-class are ordered using kcore_sum \n");
    cout << "|V| = " << VS.size() << ", f(v_0)=" << VS[0].get_bound() << ",  f(v_n)=" << VS[VS.size()-1].get_bound() <<endl;

    vector< set< vector<int> > > cliques_local(p.threads);

    int cutoff = 1;
    if (n > 20)
        cutoff = 3;

    printf("cutoff = %d, number of threads = %d \n", cutoff, p.threads);


    // 1. algorithm does not prune vertices since each worker must return largest clique each iteration
    // 2. degree is updated on the fly, then used to skip each vertex with degree <= cutoff.
    // 3. initial |P \cup \{v\}| > cutoff, where cutoff is usually relatively small such as 3.
    // 4. order color classes by \sum_{v \in S_i} \K(v), then order vertices within each S_i by k-cores

    #pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, pruned) \
        private(v, k) \
        firstprivate(ind)
    for (i = 0; i < G.num_vertices(); ++i) {

        v = VS[i].get_id();

        if (G.degree[v] <= cutoff) {
            pruned[v] = 1;
            continue;
        }


        vector<Vertex> P;
        P.reserve(G.get_max_degree()+1);

        for (long long j = G.vertices[v]; j < G.vertices[v + 1]; ++j)
            if (!pruned[G.edges[j]])
                P.push_back( Vertex(G.edges[j], compute_heuristic(G.edges[j],G)) );


        if (P.size() + 1 > cutoff) {

            vector<int> C;
            C.reserve(G.get_max_degree()+1);
            C.push_back(v);



            if (local_small_to_large)
                std::sort(P.begin(), P.end(), incr_heur);
            else
                std::sort(P.begin(), P.end(), desc_heur);
            //                    G.vertex_bucket_sort(P,local_small_to_large);

            // find largest clique at each vertex
            while (P.size() > 0) {
                int u = P.back().get_id();
                P.pop_back();

                C.push_back(u);


                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j)
                    ind[G.edges[j]] = v+1;

                vector <Vertex> R;
                R.reserve(P.size());

                for (int w = 0; w < P.size(); ++w)
                    if (ind[P[w].get_id()] == v+1)
                        if (!pruned[P[w].get_id()])
                            // no pruning, since we want to return the largest clique each iteration.
                            //                                if (G.kcore[P[i].get_id()] > max)
                            R.push_back(P[w]);

                if (R.size()) {
//                    printf("[ mcpack worker %d ] BEFORE:  |P| = %lu, |R| = %lu \n",
//                            omp_get_thread_num(), P.size(), R.size());
                    P = R;
//                    printf("[ mcpack worker %d ] AFTER:   |P| = %lu, |R| = %lu \n",
//                            omp_get_thread_num(), P.size(), R.size());
                }

                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j)
                    ind[G.edges[j]] = 0;

            } // end while


            // sort clique C, then add to set
            if (n <= 20)
                printf("\t [mcpack worker %d]  searched v=%d, found clique of size %lu \n",
                    omp_get_thread_num(), v, C.size());

            std::sort(C.begin(), C.end());
            cliques_local[omp_get_thread_num()].insert(C);

            // implicitly prune each vertex in C from the graph
            for (int h = 0; h < C.size(); ++h) {
                pruned[C[h]] = 1;
                //                    degree[C[h]] = 0;
            }

            // for each vertex in C, decrease the degree of each neighbor by 1
            for (int h = 0; h < C.size(); ++h) {
                int u = C[h];
                for (long long j = G.vertices[u]; j < G.vertices[u + 1]; ++j) {
                    G.degree[G.edges[j]]--;
                }
            }

        } // end if B(v) > max
//        pruned[v] = 1;

    } // for each vertex in order


    print_line(80);
    double s = tic();

    cliques = cliques_local[0];
    for (int t = 1; t < p.threads; ++t) {

        set< vector<int> >::iterator it;
        for( it = cliques_local[t].begin(); it != cliques_local[t].end(); it++) {
            const vector<int>& clq = (*it);

            // find max clique found overall, useful to report
            if (clq.size() > max)
                max = clq.size();

            cliques.insert(clq);
        }

        printf("[ mcpack worker %d ]  found %lu cliques of size %d, total unique cliques = %lu \n",
                t, cliques_local[t].size(), max, cliques.size());
    }
    toc(s);
    printf("took %lg seconds to find unique cliques \n\n", s);


    if (cliques.size() < 10) {

        set< vector<int> >::iterator it;
        for( it = cliques.begin(); it != cliques.end(); it++) {
            const vector<int>& clq = (*it);
            printf("|C| = %lu,  C = {", clq.size());
            for (int i = 0; i < clq.size(); ++i) {
                if (i == clq.size()-1)  printf("%d} \n", clq[i]);
                else    printf("%d, ", clq[i]);
            }
        }

    }

    print_line(80);

    if (pruned) delete[] pruned;
    return max;
}

int graphpack_heu::search_for_large_clique(graphpack_graph& G, vector<int>& sol, params & p) {
    V = G.get_vertices();
    E = G.get_edges();
    degree = G.get_degree();
    int mc = 0, i = 0, u = 0;
    if (p.k>0) mc=p.k-1;

    // set to worst case bound of cores/coloring
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);

    vector<int> C, C_max;
    C.reserve(G.get_max_degree()+1);
    C_max.reserve(G.get_max_degree()+1);

    int* pruned = new int[G.num_vertices()];
    memset(pruned, 0, G.num_vertices() * sizeof(int));

    // order verts for our search routine
    vector<Vertex> VS;
    VS.reserve(G.num_vertices());
    compute_ordering(VS,G,global_ordering,p.heu_global_s2l,0,0);
    if (G.verbose) cout << "|V| = " << VS.size() <<endl;

    vector<int> ind(G.num_vertices(),0);

    if (p.heuristic_type == "pruning" || p.heuristic_type == "pruned") {
        if (local_small_to_large) {
            printf("[graphpack: heuristic]  small to large, PRUNED \n");
            #pragma omp parallel for schedule(dynamic,block_size) shared(G, VS, mc, C_max, pruned) \
                firstprivate(ind) private(u, P, C)
            for (i = 0; i < (VS.size()) - (mc-1); ++i) {

                u = VS[i].get_id();
                if ((*K)[u] > mc) {
                    P.push_back(VS[i]);
                    for (long long j = (*V)[u]; j < (*V)[u + 1]; ++j)
                        if (!pruned[(*E)[j]])
                            if ((*K)[(*E)[j]] > mc)
                                P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );

                    if (P.size() > mc) {
                        std::sort(P.begin(), P.end(), incr_heur);
                        clique_search(P, ind, C, C_max, mc);
                    }
                    P = T;
                }
                pruned[u] = 1;
            }
        }
        else {
            printf("[graphpack: heuristic]  large to small, PRUNED \n");
            #pragma omp parallel for schedule(dynamic,block_size) shared(G, VS, mc, C_max, pruned) \
                firstprivate(ind) private(u, P, C)
            for (i = 0; i < (VS.size()) - (mc-1); ++i) {

                u = VS[i].get_id();
                if ((*K)[u] > mc) {
                    P.push_back(VS[i]);
                    for (long long j = (*V)[u]; j < (*V)[u + 1]; ++j)
                        if (!pruned[(*E)[j]])
                            if ((*K)[(*E)[j]] > mc)
                                P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                    if (P.size() > mc) {
                        std::sort(P.begin(), P.end(), desc_heur);
                        clique_search(P, ind, C, C_max, mc);
                    }
                    P = T;
                }
                pruned[u] = 1;
            }
        }
    }
    else {
        if (local_small_to_large) {

            printf("[graphpack: heuristic]  small to large, INPLACE \n");
            #pragma omp parallel for schedule(dynamic,block_size) shared(G, VS, mc, C_max, pruned) \
                firstprivate(ind) private(u, P, C)
            for (i = 0; i < (VS.size()) - (mc-1); ++i) {

                u = VS[i].get_id();
                if ((*K)[u] > mc) {
                    P.push_back(VS[i]);
                    for (long long j = (*V)[u]; j < (*V)[u + 1]; ++j)
//                        if (!pruned[(*E)[j]])
                            if ((*K)[(*E)[j]] > mc)
                                P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );


                    if (P.size() > mc) {
                        std::sort(P.begin(), P.end(), incr_heur);
                        clique_search(P, ind, C, C_max, mc);
                    }
                    P = T;
                }
//                pruned[u] = 1;
            }
        }
        else {

            printf("[graphpack: heuristic]  large to small, INPLACE \n");
            #pragma omp parallel for schedule(dynamic,block_size) shared(G, VS, mc, C_max, pruned) \
                firstprivate(ind) private(u, P, C)
            for (i = 0; i < (VS.size()) - (mc-1); ++i) {

                u = VS[i].get_id();
                if ((*K)[u] > mc) {
                    P.push_back(VS[i]);
                    for (long long j = (*V)[u]; j < (*V)[u + 1]; ++j)
//                        if (!pruned[(*E)[j]])
                            if ((*K)[(*E)[j]] > mc)
                                P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );

                    if (P.size() > mc) {
                        std::sort(P.begin(), P.end(), desc_heur);
                        clique_search(P, ind, C, C_max, mc);
                    }
                    P = T;
                }
//                pruned[u] = 1;
            }
        }
    }
    if (pruned) delete[] pruned;

    sol.resize(mc);
    for (int i = 0; i < C_max.size(); i++)  sol[i] = C_max[i];
    if (G.verbose)
        G.print_break();
    return sol.size();
}

inline
void graphpack_heu::clique_search(
        vector<Vertex> &P,
        vector<int>& ind,
        vector<int>& C,
        vector<int>& C_max,
        int& mc) {

    // terminating condition
    if (C.size() + P.size() > mc) {
        int v = P.back().get_id();
        P.pop_back();
        C.push_back(v);

        vector<Vertex> R;
        R.reserve(P.size());

        for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)   ind[(*E)[j]] = v+1;

        // intersection of N(v) and P - {v}
        for (int k = 0; k < P.size(); k++)
            if (ind[P[k].get_id()] == v+1)
                if ((*K)[P[k].get_id()] > mc)
                    R.push_back(P[k]);

        for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)  ind[(*E)[j]] = 0;


        if (R.size() > 0) {
            clique_search(R, ind, C, C_max, mc);
        }
        else if (C.size() > mc) {
            // obtain lock
            #pragma omp critical (update_mc)
            if (C.size() > mc) {
                // ensure updated max is flushed
                mc = C.size();
                C_max = C;
                if (verbose)    print_mc_info(C,sec);
            }
        }
        // backtrack and search another branch
        R.clear();
        C.pop_back();
    }
    else return;
}

inline
int graphpack_heu::compute_heuristic(int u, graphpack_graph & G) {
    switch (local_ordering) {
        case NATURAL: { // natural order (read from disk)
            return u;
            break;
        }
        case DEGREE: {
            return G.vertices[u + 1] - G.vertices[u];
            break;
        }
        case WEDGES: {
            int val = G.vertices[u + 1] - G.vertices[u]; //store degree in val
            return (val*(val-1))/2;
            break;
        }
        case KCORE: {
            return G.kcore[u];
            break;
        }
        case TRIANGLES: {
            return G.t[u];
            break;
        }
        case RAND: {
//                    return drand() % G.vertices.size();
            return rand() % G.vertices.size();
            break;
        }
        case VAR: {
            return G.kcore[u] * ((int)G.degree[u]/G.kcore[u]);
            break;
        }
        case DEGREE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]];//G.vertex_degree(edges[j]);
            }
            return val;
            break;
        }
        case KCORE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + G.kcore[G.edges[j]];
            }
            return val;
            break;
        }
        case TRIANGLE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + G.t[G.edges[j]];
            }
            return val;
            break;
        }
        case TRIANGLE_CORE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + G.tri_core[G.eid[j]];
            }
            return val;
            break;
        }
        case TRIANGLE_CORE_MAX: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                if (G.tri_core[G.eid[j]] > val) {
                    val = G.tri_core[G.eid[j]];
                }
            }
            return val;
            break;
        }
//                case TRIANGLE_CORE_MIN: {
//                    val = G.num_vertices();
//                    for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
//                        if (G.tri_core[G.eid[j]] < val) {
//                            val = G.tri_core[G.eid[j]];
//                        }
//                    }
//                    break;
//                }
        case KCORE_DEG: {
            return G.degree[u] * G.kcore[u];
            break;
        }
        case DEGREE_TRIANGLES: {
            return G.degree[u] * G.t[u];
            break;
        }
        case KCORE_TRIANGLES: {
            return G.kcore[u] * G.t[u];
            break;
        }
        case KCORE_DEG_TRI: {
            return G.degree[u] * G.kcore[u] * G.t[u];
            break;
        }
        case DEGREE_KCORE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + ((G.vertices[G.edges[j] + 1] - G.vertices[G.edges[j]]) * G.kcore[G.edges[j]]);
            }
            return val;
            break;
        }
        case KCORE_TRIANGLE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + (G.kcore[G.edges[j]] * G.t[G.edges[j]]);
            }
            return val;
            break;
        }
        case DEGREE_KCORE_TRIANGLE_VOL: {
            double val = 0;
            for (long long j = G.vertices[u]; j < G.vertices[u + 1]; j++) {
                val = val + (G.degree[G.edges[j]] * G.kcore[G.edges[j]] * G.t[G.edges[j]]);
            }
            return val;
            break;
        }
        default: {
            return G.vertices[u + 1] - G.vertices[u];
            break;
        }
    }
}

int graphpack_heu::search_cores(graphpack_graph& G, vector<int>& C_max, int lb) {

    vector <int> C, X;
    C.reserve(ub);
    C_max.reserve(ub);
    vector<Vertex> P, T;
    P.reserve(G.get_max_degree()+1);
    T.reserve(G.get_max_degree()+1);
    vector<int> ind(G.num_vertices(),0);

    int mc = lb, mc_prev, mc_cur, i;
    int lb_idx = 0, v = 0;
    for (i = G.num_vertices()-1; i >= 0; i--) {
        v = (*order)[i];
        if ((*K)[v] == lb)   lb_idx = i;
    }

    #pragma omp parallel for schedule(dynamic,block_size) \
        shared(G, X, mc, C_max) private(i, v, P, mc_prev, mc_cur, C) firstprivate(ind)
    for (i = lb_idx; i <= G.num_vertices()-1; i++) {

        v = (*order)[i];
        mc_prev = mc_cur = mc;

        if ((*K)[v] > mc_cur) {
            for (long long j = (*V)[v]; j < (*V)[v + 1]; j++)
                if ((*K)[(*E)[j]] > mc_cur)
                    P.push_back( Vertex((*E)[j], compute_heuristic((*E)[j],G)) );

            if (P.size() > mc_cur) {
                std::sort(P.begin(), P.end(), incr_heur);
                branch(P, 1 , mc_cur, C, ind);

                if (mc_cur > mc_prev) {
                    if (mc < mc_cur) {
                        #pragma omp critical
                        if (mc < mc_cur) {
                            mc = mc_cur;
                            C.push_back(v);
                            C_max = C;
                            if (G.verbose)
                                print_info(C_max);
                        }
                    }
                }
            }
        }
        C = X; P = T;
    }
    C.clear();
    if (G.verbose)
        cout << "[search_cores]\t mc = " << mc <<endl;
    return mc;
}

int graphpack_heu::search(graphpack_graph& G, vector<int>& C_max, params & p) {


    if (p.heu_global_ordering == "kcore" || p.heu_global_ordering == "kcores") {
        return search_bounds(G, C_max); // global ordering is k-core by default
    }
    else {
        return search_bounds_ordering(G,C_max,p);
    }
}

// global ordering is k-core by default
int graphpack_heu::search(graphpack_graph& G, vector<int>& C_max) {
        return search_bounds(G, C_max);
}

void graphpack_heu::print_info(vector<int> C_max) {
    if (verbose) {
        cout << "*** [mcpack heuristic: thread " << omp_get_thread_num() + 1;
        cout << "]   current max clique = " << C_max.size();
        cout << ",  time = " << get_time() - sec << " sec" <<endl;
    }
}

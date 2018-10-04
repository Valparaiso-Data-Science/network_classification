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

void graphpack_graph::induced_cores_ordering(
        vector<long long>& V,
        vector<int>& E,
        int* &pruned) {

    long long n, d, i, j, start, num, md;
    long long v, u, w, du, pu, pw, md_end;
    n = vertices.size();

    vector <int> pos_tmp(n);
    vector <int> core_tmp(n);
    vector <int> order_tmp(n);

    md = 0;
    for(v=1; v<n; v++) {
        core_tmp[v] = V[v] - V[v-1];
        if (core_tmp[v] > md)  md = core_tmp[v];
    }

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[core_tmp[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos_tmp[v] = bin[core_tmp[v]];
        order_tmp[pos_tmp[v]] = v;
        bin[core_tmp[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (i = 1; i < n; i++) {
        v=order_tmp[i];
        for (j = V[v-1]; j < V[v]; j++) {
            u = E[j] + 1;
            if (core_tmp[u] > core_tmp[v]) {
                du = core_tmp[u];   pu = pos_tmp[u];
                pw = bin[du];       w = order_tmp[pw];
                if (u != w) {
                    pos_tmp[u] = pw;   order_tmp[pu] = w;
                    pos_tmp[w] = pu;   order_tmp[pw] = u;
                }
                bin[du]++;   core_tmp[u]--;
            }
        }
    }

    for (v=0; v<n-1; v++) {
        core_tmp[v] = core_tmp[v+1]+1;
        order_tmp[v] = order_tmp[v+1]-1;
    }

    // added
    max_core = kcore[kcore_order[num_vertices()-1]] - 1;

    kcore = core_tmp;
    kcore_order = order_tmp;
    bin.clear();
}

void graphpack_graph::induced_cores_ordering(
        vector<long long>& V,
        vector<int>& E,
        vector<int> & pruned) {

    long long n, d, i, j, start, num, md;
    long long v, u, w, du, pu, pw, md_end;
    n = vertices.size();

    vector <int> pos_tmp(n);
    vector <int> core_tmp(n);
    vector <int> order_tmp(n);

    md = 0;
    for(v=1; v<n; v++) {
        core_tmp[v] = V[v] - V[v-1];
        if (core_tmp[v] > md)  md = core_tmp[v];
    }

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[core_tmp[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos_tmp[v] = bin[core_tmp[v]];
        order_tmp[pos_tmp[v]] = v;
        bin[core_tmp[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (i = 1; i < n; i++) {
        v=order_tmp[i];
        for (j = V[v-1]; j < V[v]; j++) {
            u = E[j] + 1;
            if (core_tmp[u] > core_tmp[v]) {
                du = core_tmp[u];   pu = pos_tmp[u];
                pw = bin[du];       w = order_tmp[pw];
                if (u != w) {
                    pos_tmp[u] = pw;   order_tmp[pu] = w;
                    pos_tmp[w] = pu;   order_tmp[pw] = u;
                }
                bin[du]++;   core_tmp[u]--;
            }
        }
    }

    for (v=0; v<n-1; v++) {
        core_tmp[v] = core_tmp[v+1]+1;
        order_tmp[v] = order_tmp[v+1]-1;
    }

    max_core = kcore[kcore_order[num_vertices()-1]] - 1;

    kcore = core_tmp;
    kcore_order = order_tmp;
    bin.clear();
}

// normalized for maximum clique upper bound
void graphpack_graph::compute_cores() {
    long long j;
    int n, d, i, start, num, md;
    int v, u, w, du, pu, pw, md_end;
    n = vertices.size();

    vector <int> pos(n);
    if (kcore_order.size() > 0) {
        vector<int> tmp(n,0);
        kcore = tmp;
        kcore_order = tmp;
    }
    else {
        kcore_order.resize(n);
        kcore.resize(n);
    }

    md = 0;
    for (v=1; v<n; v++) {
        kcore[v] = vertices[v] - vertices[v-1];
        if (kcore[v] > md)  md = kcore[v];
    }

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[kcore[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    // bucket sort
    for (v=1; v<n; v++) {
        pos[v] = bin[kcore[v]];
        kcore_order[pos[v]] = v;
        bin[kcore[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    // kcores
    for (i=1; i<n; i++) {
        v=kcore_order[i];
        for (j=vertices[v-1]; j<vertices[v]; j++) {
            u = edges[j] + 1;
            if (kcore[u] > kcore[v]) {
                du = kcore[u];   pu = pos[u];
                pw = bin[du];    w = kcore_order[pw];
                if (u != w) {
                    pos[u] = pw;   kcore_order[pu] = w;
                    pos[w] = pu;   kcore_order[pw] = u;
                }
                bin[du]++;   kcore[u]--;
            }
        }
    }

    for (v = 0; v < n-1; v++) {
        kcore[v] = kcore[v+1] + 1; // K + 1
        kcore_order[v] = kcore_order[v+1]-1;
    }
    max_core = kcore[kcore_order[num_vertices()-1]] - 1;

    bin.clear();
    pos.clear();
}

void graphpack_graph::compute_cores(vector<int> &pruned) {

    long long n, d, i, j, start, num, md;
    long long v, u, w, du, pu, pw, md_end;
    n = vertices.size();
//    kcore.resize(n);
//    fill(kcore.begin(), kcore.end(), 0);
    vector <int> pos_tmp(n);
//    vector <int> order_tmp(n,0);

    if (kcore_order.size() > 0) {
        vector<int> tmp(n,0);
        kcore = tmp;
        kcore_order = tmp;
    }
    else {
        kcore_order.resize(n);
        kcore.resize(n);
    }

    md = 0;
    for(v=1; v<n; v++) {
        if (!pruned[v-1]) {
            kcore[v] = degree[v-1];
            if (kcore[v] > md)  md = kcore[v];
        }
    }
    if (verbose)
        printf("[cores]  max degree = %lld \n", md);

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[kcore[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos_tmp[v] = bin[kcore[v]];
        kcore_order[pos_tmp[v]] = v;
        bin[kcore[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (i = 1; i < n; i++) {
        v=kcore_order[i];
        if (!pruned[v-1]) {
            for (j = vertices[v-1]; j < vertices[v]; j++) {
                if (!pruned[edges[j]]) {
                    u = edges[j] + 1;
                    if (kcore[u] > kcore[v]) {
                        du = kcore[u];   pu = pos_tmp[u];
                        pw = bin[du];  w = kcore_order[pw];
                        if (u != w) {
                            pos_tmp[u] = pw;   kcore_order[pu] = w;
                            pos_tmp[w] = pu;   kcore_order[pw] = u;
                        }
                        bin[du]++;   kcore[u]--;
                    }
                }
            }
        }
    }

    max_core = 0;
    for (v=0; v<n-1; v++) {
        if (!pruned[v]) {
            kcore[v] = kcore[v+1] + 1; // K+1
            kcore_order[v] = kcore_order[v+1]-1;
            if (kcore[v] > max_core)  max_core = kcore[v];
        }
        else kcore[v] = 0;
    }
    if (verbose)    cout << "[graphpack: updated cores]  K: " << max_core <<endl;

//    kcore_order = kcore_order;

    bin.clear();
    pos_tmp.clear();
//    kcore_order.clear();
}

void graphpack_graph::update_kcores(int* &pruned) {

    long long n, d, i, j, start, num, md;
    long long v, u, w, du, pu, pw, md_end;
    n = vertices.size();
    kcore.resize(n);
    fill(kcore.begin(), kcore.end(), 0);
    vector <int> pos_tmp(n);
    vector <int> order_tmp(n);

    md = 0;
    for(v=1; v<n; v++) {
        if (!pruned[v-1]) {
            kcore[v] = degree[v-1];
            if (kcore[v] > md)  md = kcore[v];
        }
    }

    md_end = md+1;
    vector < int > bin(md_end,0);

    for (v=1; v < n; v++)  bin[kcore[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos_tmp[v] = bin[kcore[v]];
        order_tmp[pos_tmp[v]] = v;
        bin[kcore[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (i = 1; i < n; i++) {
        v=order_tmp[i];
        if (!pruned[v-1]) {
            for (j = vertices[v-1]; j < vertices[v]; j++) {
                if (!pruned[edges[j]]) {
                    u = edges[j] + 1;
                    if (kcore[u] > kcore[v]) {
                        du = kcore[u];   pu = pos_tmp[u];
                        pw = bin[du];  w = order_tmp[pw];
                        if (u != w) {
                            pos_tmp[u] = pw;   order_tmp[pu] = w;
                            pos_tmp[w] = pu;   order_tmp[pw] = u;
                        }
                        bin[du]++;   kcore[u]--;
                    }
                }
            }
        }
    }

    max_core = 0;
    for (v=0; v<n-1; v++) {
        if (!pruned[v]) {
            kcore[v] = kcore[v+1] + 1; // K+1
            order_tmp[v] = order_tmp[v+1]-1;
            if (kcore[v] > max_core)  max_core = kcore[v];
        }
        else kcore[v] = 0;
    }
    cout << "[graphpack: updated cores]  K: " << max_core <<endl;

    bin.clear();
    pos_tmp.clear();
    order_tmp.clear();
}

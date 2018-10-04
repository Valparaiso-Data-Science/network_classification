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

#ifndef GRAPHPACK_UTILS_H_
#define GRAPHPACK_UTILS_H_

#include <cstddef>
#include <sys/time.h>
#include <unistd.h>
#include <iostream>
#include "assert.h"
#include <sys/types.h>
#include <dirent.h>
#include <errno.h>
#include <string>
#include <vector>
#include <set>
#include "graphpack_headers.h"


using namespace std;

bool fexists(const char *filename);
void usage(char *argv0);
double get_time();
double tic();
void toc(double & start);
string memory_usage();
void validate(bool condition, const string& msg);
void indent(int level);
void indent(int level, string str);
void print_max_clique(vector<int>& max_clique_data);
void print_n_maxcliques(set< vector<int> > C, int n);
void print_line(int n);
string extract_filename(string fn, bool remove_ext = true);
string remove_file_extension(string fn);
int getdir (string dir, vector<string> &files);

// must include templates in headers
template<typename T>
void write_results(std::vector<T> &data, string &filename, bool output_id = false) {

    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);

    if (output_id) {
        for (long long e = 0; e < data.size(); e++) {
            myfile << e << "," << data[e] << "\n";
        }
    }
    else {

        for (long long e = 0; e < data.size(); e++) {
            myfile << data[e] << "\n";
        }
    }
    myfile.close();
}

std::ifstream::pos_type get_filesize(const char* filename);

template<typename T>
void bin_values(std::vector<T> &data, vector<int> & bin) {
//    vector<int> bin_deg(G.max_degree,0);
    printf("[utils]  num verts = %lu \n", data.size());
    for (int v = 0; v < data.size(); ++v) {
        int val = data[v];
        bin[val]++;
    }
}

// use to save space if vector contains many zeros
template<typename T>
void write_sparse_results(std::vector<T> &data, string &filename) {

    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);

    for (long long i = 0; i < data.size(); i++) {
        if (data[i] != 0) {
            myfile << i << " " << data[i] << "\n";
        }
    }
    myfile.close();
}

template<typename T>
void write_vector(std::vector<T> &data, string suffix);

#endif

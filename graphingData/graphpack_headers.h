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

#ifndef GRAPHPACK_HEADERS_H_
#define GRAPHPACK_HEADERS_H_

#include <iostream>
#include <stdlib.h>
#include <vector>
#include <time.h>
#include <string.h>
#include <stdio.h>
#include <map>
#include <fstream>
#include <sstream>
#include <stdint.h>

#ifdef _OPENMP
#  include <omp.h>
#else
int omp_get_max_threads()       { return 1; }
void omp_set_num_threads(int)   {}
int omp_get_thread_num()        { return 1; }
#endif

using namespace std;

// common notations
typedef vector<double>  vec;                // vector of doubles
typedef vector<int>     vec_i;              // vector of int
typedef vector<short>   vec_s;              // vector of short ints
typedef vector<bool>    vec_b;              // vector of bools
typedef vector< vector<double> > dmat;  // dense matrix (or double matrix =])

//typedef vector<int>     index;           // vector of int

typedef vector<int>     hash;           // vector of int
typedef vector<double>  hash_d;        // vector of int


#ifndef LINE_LENGTH
#define LINE_LENGTH 256
#endif
#define NANOSECOND 1000000000

#define DEBUG
#define DEBUG_TCORE

//typedef uint64_t int;
//typedef long long         uint64_t;

#endif

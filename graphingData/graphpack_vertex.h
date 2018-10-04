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

#ifndef GRAPHPACK_VERTEX_H_
#define GRAPHPACK_VERTEX_H_

using namespace std;

namespace graphpack {

    class Vertex {
        private:
            int id, b;
        public:
            Vertex(int vertex_id, int bound): id(vertex_id), b(bound) {};

            void set_id(int vid)        { id = vid; }
            int get_id()                { return id; }

            void set_bound(int value)   { b = value; }
            int get_bound()             { return b; }
    };

    static bool decr_bound(Vertex v,  Vertex u) {
        return (v.get_bound() > u.get_bound() ||
                (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
    }
    static bool incr_bound(Vertex v,  Vertex u) {
        return (v.get_bound() < u.get_bound() ||
                (v.get_bound() == u.get_bound() && v.get_id() > u.get_id()));
    };

    typedef vector<Vertex>     VertexSet;

};
#endif

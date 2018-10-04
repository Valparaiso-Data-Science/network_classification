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

using namespace graphpack;
using namespace std;


void graphpack_graph::initialize() {
    max_degree = 0;
    min_degree = 0;
    avg_degree = 0;
    max_core = 0;
    max_tri_core = 0;
    max_t_edge = 0;
    total_t = 0;
    max_t = 0;
    block_size = 64;
    fn = "";
    is_gstats = false;
    adj = NULL;
    is_dense = false;
    verbose = false;
}

void graphpack_graph::read_graph(const string& filename) {
    // NOTE: flag is_gstats is set to true for approx triangle cores.
    // Need e_u and e_v, and these are computed if is_gstats is true, when -s option is used!
    is_gstats = true;
    fn = filename;
    double sec = get_time();
    string ext = get_file_extension(filename);

    cout << "[graphpack: graph reader]  All graphs are assumed to be undirected" <<endl;
    cout << "[graphpack: graph reader]  Self-loops and weights (if any) are discarded" <<endl;

    if (ext == "edges" || ext == "eg2" || ext == "txt" || ext == "csv") {
        cout << "[graphpack: general graph reader]  reading the edge list" <<endl;
        read_edge_list(filename);
//        read_edges(filename);
    }
    else if (ext == "mtx") {
        cout << "[graphpack: mtx graph reader]  Assuming matrix is undirected, and upper-triangular " <<endl;
        read_mtx(filename);
    }
    else if (ext == "gr") {
        read_metis(filename);
    }
    else {
        cout << "[graphpack: general graph reader] Unsupported graph format. Attempting to read the graph." <<endl;
        read_edge_list(filename);
        return;
    }
    cout << "Reading time " << get_time() - sec << endl;
    basic_stats();

    // if ()
    // create_adj()
}

bool graphpack_graph::detect_weighted_graph(string & line, string & delim) {

    int num_tokens = 0;

    string buf; // Have a buffer string
    stringstream ss(line); // Insert the string into a stream

    vector<string> tokens; // Create vector to hold our words

    while (ss >> buf) {
        tokens.push_back(buf);
    }
    printf("number of tokens in line = %lu \n", tokens.size());


    if (tokens.size() == 3) return true; // weighted graph (3rd column)
    return false;   // unweighted, only edge list with two columns
}

void graphpack_graph::detect_delim(string & line, string & delim) {

    // try to figure out the delimiter (first line of file)
    if (delim == "") {
        std::size_t prev = 0, pos;
        std::string tab_spaces("   ");
//        std::size_t found = line.find(tab_spaces);
        if ((pos = line.find_first_of(',', prev)) != std::string::npos) {
            delim = ',';
        }
        else if ((pos = line.find_first_of('\t', prev)) != std::string::npos) {
            delim = '\t';
        }
        else if ((pos = line.find(tab_spaces)) != std::string::npos) {
//            std::string tab_spaces("   ");
//            std::size_t found = line.find(tab_spaces);
            printf("found tab-spaces delimiter \n");
            delim = "   ";
        }
        else if ((pos = line.find_first_of(' ', prev)) != std::string::npos) {
            delim = ' ';
        }
    }

    if (delim == "") {
        if (get_file_extension(fn) == "csv")
            delim = ',';
        else if (get_file_extension(fn) == "tab")
            delim = '\t';
        else if (get_file_extension(fn) == "mtx")
            delim = ' ';
        else
            delim = ' '; // set as default
        cout << "[mcpack]  no delimiter recognized, using \"" << delim.c_str()[0] << "\" as delimiter!" <<endl;
    }
    else
        cout << "[mcpack]  detected \"" << delim << "\" as delimiter! " <<endl;
//        cout << "[mcpack]  detected \"" << delim.c_str()[0] << "\" as delimiter! " <<endl;
}

inline
void graphpack_graph::get_token(int & v, string & line, string & delim, size_t & pos, size_t & prev) {
    // check line for delimiter
//    if ((pos = line.find_first_of(delim, prev)) != std::string::npos) {
    if ((pos = line.find(delim, prev)) != std::string::npos) {
        if (pos > prev) {
            v = atoi(line.substr(prev, pos-prev).c_str());
        }
        prev = pos+1; // set current pos
    }
    else if (prev < line.length()) // last token in the line
        v = atoi(line.substr(prev, std::string::npos).c_str());
}


inline
void graphpack_graph::get_token(double & weight, string & line, string & delim, size_t & pos, size_t & prev, bool & is_weighted_graph) {
    // check line for delimiter
//    if ((pos = line.find_first_of(delim, prev)) != std::string::npos) {
    if ((pos = line.find(delim, prev)) != std::string::npos) {
        if (pos > prev) {
            weight = atof(line.substr(prev, pos-prev).c_str());
        }
        prev = pos+1; // set current pos
    }
    else if (prev < line.length()) // last token in the line
        weight = atof(line.substr(prev, std::string::npos).c_str());
}



/**
 * @brief:
 * reads a more general edge list, limited assumptions about the graph
 *
 * WEIGHTS: All weights are discarded, unless the graph is temporal
 * LABELS:  Vertices are relabeled, and the old ids are discarded, unless specified.
 * TODO: Input file has vertices with gaps
 * TODO: Detect if graph is temporal, and add timestamps to wt vector
 */
void graphpack_graph::read_edge_list(const string& filename) {
    string line = "";
    map< int, vector<int> > vert_list;
    map< int, vector<double> > wt_vert_list;
    int v = -1, u = -1, num_es = 0, self_edges = 0;
    string delimiters = " ,\t"; // spaces, commas, or tabs
    string delim = "";
//    int weights = 3;        // weights should be the 3rd column, which is ignored.
    double weight;
    string graph_exif = ""; // data about the graph (if any)

    ifstream file (filename.c_str());
    if (!file) { cout << filename << "File not found!" <<endl; return; }

    // check if vertex ids start at 0 or 1
    is_weighted = false;
    bool fix_start_idx = true;
    bool ignore_first_line = false;
    string token;
    stringstream iss;

    // save graph info/comments at top of file
    while(std::getline(file, line) && (line[0] == '#' || line[0] == '%')) {
        graph_exif += line;
        if (line.find("MatrixMarket matrix coordinate pattern symmetric") != std::string::npos) {
            // these files should be delimited by space and first line contains <num_nodes,num_nodes,num_edges>
            delim = ' ';
            ignore_first_line = true;
        }
    }

    int num_verts = 0, num_edges = 0;
    if (get_file_extension(filename) == ".mtx") {
        cout << "[graphpack: graph reader]  mtx file detected!" <<endl;
        iss << line;
        int cols = 0;
        iss >> num_verts >> cols >> num_edges;

        if(num_verts!=cols) {
            cout<<"[mcpack]  error: this is not a square matrix, attempting to proceed."<<endl;
        }
    }

    // detect the delim for reading the graph
    detect_delim(line,delim);

    // detect if line has three columns, third is assumed to be for weights
    is_weighted = detect_weighted_graph(line,delim);
    if (is_weighted)    printf("weighted graph detected \n");


    // handle the first line (find starting vertex id)
    if (line != "") {
        iss.clear();
        iss.str(line);
        iss >> v >> u; //>> weight;

        if (v == 0 || u == 0) {
            fix_start_idx = false;
        }
//        cout << "weight = " << weight <<endl;
    }

    int max = 0; // largest vertex id (assumed to be ints)
    print_line(80);
    cout << "[graphpack: graph reader]  reading a general edge list" <<endl;
    cout << "[graphpack: graph reader]     - multiple edges are removed" <<endl;
    cout << "[graphpack: graph reader]     - self-loops are removed" <<endl;
    cout << "[graphpack: graph reader]     - weight/time is stored in the third column, 4th col is igored if it exists" <<endl;
    cout << "[graphpack: graph reader]     - vertex ids are adjusted to start from 0 (if not already)" <<endl;
    cout << "[graphpack: graph reader]     - gaps in vertex ids are removed (on the fly)" <<endl;
    cout << "[graphpack: graph reader]     - vertex/edge wts, time, and the vertex lookup table are not stored, unless specified" <<endl;
    print_line(80);
    // find starting vertex id,
    // compute the number of vertices to expect (since gaps in vertex ids are possible)
    while(std::getline(file, line)) { // && fix_start_idx) {
        if (line != "") { // ensure line actually contains data
            iss << line;

            // ignore comments
            if (line[0] == '%' || line[0] == '#') continue;

            std::size_t prev = 0, pos;
            get_token(v,line,delim,pos,prev);
            get_token(u,line,delim,pos,prev);

            if (v == 0 || u == 0) {
                fix_start_idx = false;
            }
            // get max vertex id
            if (v > max) max = v;
            if (u > max) max = u;
        }
    }
    cout << "[graphpack: graph reader]  largest vertex id is " << max <<endl;
    //    }

    file.close();
    if (fix_start_idx) cout << "[graphpack: graph reader]  vertex ids from the file begin at 1" <<endl;
    else cout << "[graphpack: graph reader]  vertex ids begin at 0" <<endl;


    // mtx file had a zero vertex id, then we must perform some additional work
//    if (get_file_extension(filename) == ".mtx" && !fix_start_idx) {
//        cout << "[graphpack: graph reader]   mtx format error detected! vertex ids begin at zero! fixing it." <<endl;
//        map< int, vector<int> > vert_list_new;
//        // reset list
//        vert_list = vert_list_new;
//        max = num_verts;
//    }


    ifstream fin (filename.c_str());
    if (!fin) { cout << filename << "File not found!" <<endl; return; }

    int vertex_id = 0;
    vector<int> vert_lookup(max+1,-1); // add 1 since vertex ids may start at zero
//    vertex_lookup.resize(max+1,-1);

    if (is_weighted) {

        // uncomment for lookup
//        vertex_lookup.resize(max+1,-1);

//        wt_vert_list;
        while(std::getline(fin, line)) {
            if (line != "") { // ensure line actually contains data
                iss << line;

                // ignore comments
                if (line[0] == '%' || line[0] == '#') continue;

                std::size_t prev = 0, pos; // prev is last location in the line
                get_token(v,line,delim,pos,prev);
                get_token(u,line,delim,pos,prev);

                // get the weight (3rd column)
//                get_token(weight,line,delim,pos,prev,is_weighted);

//                cout << line << " |  (" << v << "," << u << ", weight = " << weight << ")" <<endl;

                if (fix_start_idx) {
                    v--;
                    u--;
                }
                if (v == u)  self_edges++;
                else {
                    if (vert_lookup[v] == -1) { // new vertex
                        // vertex_id is the new identifier we created (next)
                        vert_lookup[v] = vertex_id; // store the new id
//                        vertex_lookup[vertex_id] = v;
                        vertex_id++;
                    }
                    v = vert_lookup[v]; // get the consistent vertex id

                    if (vert_lookup[u] == -1) { // new vertex
                        vert_lookup[u] = vertex_id; // store the new id
//                        vertex_lookup[vertex_id] = u;
                        vertex_id++;
                    }
                    u = vert_lookup[u]; // get the consistent vertex id

    //                cout << line << " |  (" << v << "," << u << ")" <<endl;

                    // add edge
                    vert_list[v].push_back(u);
                    vert_list[u].push_back(v);

                    // add weight
                    wt_vert_list[v].push_back(weight);
                    wt_vert_list[u].push_back(weight);

    //                if (is_gstats) {
                        // e_v and e_u are used for triangle cores.
                        // edge-csc (for fast computations on the edges)
                        e_v.push_back(v);
                        e_u.push_back(u);
    //                }
                }
            }
        }
    }
    else { // unweighted graph (two columns)
        while(std::getline(fin, line)) {
            if (line != "") { // ensure line actually contains data
                iss << line;

                // ignore comments
                if (line[0] == '%' || line[0] == '#') continue;

                std::size_t prev = 0, pos; // prev is last location in the line
                get_token(v,line,delim,pos,prev);
                get_token(u,line,delim,pos,prev);

                //            cout << line << " |  (" << v << "," << u << ")" <<endl;

                if (fix_start_idx) {
                    v--;
                    u--;
                }
                if (v == u)  self_edges++;
                else {
                    if (vert_lookup[v] == -1) { // new vertex
                        vert_lookup[v] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    v = vert_lookup[v]; // get the consistent vertex id

                    if (vert_lookup[u] == -1) { // new vertex
                        vert_lookup[u] = vertex_id; // store the new id
                        vertex_id++;
                    }
                    u = vert_lookup[u]; // get the consistent vertex id

                    //                cout << line << " |  (" << v << "," << u << ")" <<endl;

                    // add edge
                    vert_list[v].push_back(u);
                    vert_list[u].push_back(v);
                    //                if (is_gstats) {
                    // e_v and e_u are used for triangle cores.
                    e_v.push_back(v);
                    e_u.push_back(u);
                    //                }
                }
            }
        }
    }
    fin.close();

    // clean up
    vector<int> tmp;
    vert_lookup.clear();
    vert_lookup = tmp;

    cout << "vert_list size: " << vert_list.size() <<endl;
    if (is_weighted) cout << "wt_vert_list size: " << wt_vert_list.size() <<endl;


    vertices.push_back(edges.size()); // inserts a 0 in first position of vertices array

    if (is_weighted) { // weighted graph (3rd column)
        for (int i=0; i < vert_list.size(); i++) {
            edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
            wt.insert(wt.end(),wt_vert_list[i].begin(),wt_vert_list[i].end());
            vertices.push_back(edges.size());
        }
        wt_vert_list.clear();
    }
    else { // unweighted
        for (int i=0; i < vert_list.size(); i++) {
            edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
            vertices.push_back(edges.size());
        }
    }
    vert_list.clear();

//    print_edges();


    /*
     * 1. temporal graphs tend to have many edges
     * 2. undirected graphs are assumed, note if edge list has either:
     *       - both (u,v) and (v,u), then the func below removes one copy
     *       - directed edge list, then it is made to be undirected
     */
    bool delete_multiple_edges = false;//true; // TODO: move this out
    if (delete_multiple_edges) {
        remove_multiple_edges();
    }

    vertex_degrees();
    cout << "self-loops: " << self_edges <<endl;

    if (vertices.size() < 10) // small graph
        print_weighted_graph();

}


/**
 * Specifically designed to be _fast_ for very large graphs
 * Impossible to store full adj of large sparse graphs, instead
 * we create a lookup table for each vertex, and build it on the fly,
 * using this info to mark and essentially remove the multiple edges
 */
void graphpack_graph::remove_multiple_edges() {

    // create index
    vector<short> ind(num_vertices(),0);

    // arrays for the processed graph
    vector<long long> vs(vertices.size(),0);
    vector<int> es;
    es.reserve(edges.size());

    int start = 0;
    for (int i = 0; i < num_vertices(); i++) { // check every vertex
        start = es.size();
        // for every neighbor, check if it has been indexed,
        // if not, add it to the new edge set and mark it, otherwise, ignore.
        for (long long j = vertices[i]; j < vertices[i + 1]; j++) {
            int u = edges[j];
            if (ind[u] == 0) {
                es.push_back(edges[j]);
                ind[u] = 1; // mark neighbor
            }
        }
        vs[i] = start;
        vs[i + 1] = es.size();

        // reset index fast using the smaller processed array
        for (long long j = vs[i]; j < vs[i + 1]; j++) {
            ind[es[j]] = 0;
        }

    }
    cout << "[graphpack: graph reader]  removed " << (edges.size() - es.size())/2 << " duplicate edges (multigraph)" <<endl;

    // set these to be the actual graph
    vertices = vs;
    edges = es;
    vs.clear();
    es.clear();
}




/**
 * Matlab Library Interface
 */
graphpack_graph::graphpack_graph(int nverts, int* heads, int* tails) {
    initialize();
    int num_vs = nverts;
    int num_es = heads[nverts];

    vector<long long> V(nverts,0);
    vector<int> E;
    E.reserve(num_es);

    int start = 0;
    for (int i = 0; i < nverts; i++) {
//        V[i] = heads[i];
        start = E.size();
        for (long long j = heads[i]; j < heads[i + 1]; j++ ) {

            E.push_back(tails[j]);
        }

        V[i] = start;
        V[i + 1] = E.size();
    }
    vertices = V;
    edges = E;

    vertex_degrees();
}


/**
 * graphpack library interface
 *
 * Basic assumptions:
 *  - for undirected graphs, we assume one unique pair for each edge
    - vertex ids are assumed to start at 0
 */
graphpack_graph::graphpack_graph(int nverts, int nedges,
        std::pair<int,int>* edge_pairs) {
    initialize();
    int num_vs = nverts;
    int num_es = nedges;

    map< int, vector<int> > vert_list;
    int v = 0, u = 0, self_edges = 0;

    for (int i = 0; i < num_es; ++i) {
        v = edge_pairs[i].first;
        u = edge_pairs[i].second;

        if (v == u)  self_edges++;
        else {
            vert_list[v].push_back(u);
            vert_list[u].push_back(v);
            if (is_gstats) {
                // e_v and e_u are used for triangle cores.
                e_v.push_back(v);
                e_u.push_back(u);
            }
        }
    }
    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
        edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
        vertices.push_back(edges.size());
    }
    vert_list.clear();
    vertex_degrees();
    cout << "self-loops: " << self_edges <<endl;

    // clean up
    vert_list.clear();
}

int graphpack_graph::read_edges(const string& filename) {
    istringstream in_stream;
    string line = "";
    map< int, vector<int> > vert_list;
    int v = 0, u = 0, num_es = 0, self_edges = 0;

    ifstream in_check (filename.c_str());
    if (!in_check) { cout << filename << "File not found!" <<endl; return 0; }

    bool fix_start_idx = true;
    while (!in_check.eof()) {
        getline(in_check,line);
        if (line[0] == '%' || line[0] == '#') continue;
        if (line != "") {
            in_stream.clear();
            in_stream.str(line);
            in_stream >> v >> u;
            if (v == 0 || u == 0) {
                fix_start_idx = false;
                break;
            }
        }
    }
    ifstream in (filename.c_str());
    if (!in) { cout << filename << "File not found!" <<endl; return 0; }

    while (!in.eof()) {
        getline(in,line);
        if (line[0] == '%' || line[0] == '#') continue;
        num_es++;
        if (line != "") {
            in_stream.clear();
            in_stream.str(line);
            in_stream >> v >> u;

            if (fix_start_idx) {
                v--;
                u--;
            }
            if (v == u)  self_edges++;
            else {
                vert_list[v].push_back(u);
                vert_list[u].push_back(v);
            }
        }
    }
    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
        edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
        vertices.push_back(edges.size());
    }
    vert_list.clear();
    vertex_degrees();
    cout << "self-loops: " << self_edges <<endl;

    return 1;
}

int graphpack_graph::read_mtx(const string& filename) {
    string line = "";
    map< int, vector<int> > vert_list;
    map< int, vector<double> > value_list;
    int v = -1, u = -1, num_es = 0, self_edges = 0;
    string delim = " ";
    string graph_exif = ""; // data about the graph (if any)

    ifstream file (filename.c_str());
    if (!file) { cout << filename << "File not found!" <<endl; return 0; }

    // check if vertex ids start at 0 or 1
    bool fix_start_idx = true;
    bool ignore_first_line = false;
    string token;
    stringstream iss;

    // save graph info/comments at top of file
    while(std::getline(file, line) && (line[0] == '#' || line[0] == '%')) {
        graph_exif += line;
    }

    int num_verts = 0, num_edges = 0;
    iss << line;
    int cols = 0;
    iss >> num_verts >> cols >> num_edges;

    if(num_verts!=cols)
        cout<<"[mcpack]  error: this is not a square matrix, attempting to proceed."<<endl;

    // handle the first line (find starting vertex id)
    if (line != "") {
        iss.clear();
        iss.str(line);
        iss >> v >> u;

        if (v == 0 || u == 0) {
            fix_start_idx = false;
        }
    }

    double value;
    cout << "[graphpack: graph reader]  reading mtx file" <<endl;
    // assume there are no gaps in vertex ids (mtx standard)
    // possibly faster, since we can break if a zero is encountered
    // find starting vertex id (break asap)
    while(std::getline(file, line) && fix_start_idx) {
        if (line != "") { // ensure line actually contains data

            // ignore comments
            if (line[0] == '%' || line[0] == '#') continue;

            iss.clear();
            iss.str(line);
            iss >> v >> u >> value;

            v--;
            u--;
            if (v == u)  self_edges++;
            else {
                vert_list[v].push_back(u);
                vert_list[u].push_back(v);
//                if (is_temporal() || is_weighted()) {
//                    value_list[v].push_back(value);
//                    value_list[u].push_back(value);
//                }
                if (is_gstats) {
                    // e_v and e_u are used for triangle cores.
                    // both must be |E|/2, since undirected. Must only have unique edges.
                    e_v.push_back(v);
                    e_u.push_back(u);
                }
            }
        }
    }
    vertices.push_back(edges.size());
    for (int i=0; i < vert_list.size(); i++) {
        edges.insert(edges.end(),vert_list[i].begin(),vert_list[i].end());
        wt.insert(wt.end(),value_list[i].begin(),value_list[i].end());
        vertices.push_back(edges.size());
    }
    vert_list.clear();
    vertex_degrees();
    cout << "self-loops: " << self_edges <<endl;
    return 1;
}

int graphpack_graph::read_mtx_file(const string& filename) {
    float connStrength = -DBL_MAX;
    istringstream in2;
    string line="";
    map<int,vector<int> > v_map;
    map<int,vector<double> > valueList;
    int col=0, row=0, ridx=0, cidx=0;
    int entry_counter = 0, num_of_entries = 0;
    double value;

    ifstream in (filename.c_str());
    if(!in) {
        cout<<filename<<" not Found!"<<endl;
        return 0;
    }
    char data[LINE_LENGTH], banner[LINE_LENGTH], mtx[LINE_LENGTH];
    char crd[LINE_LENGTH], data_type[LINE_LENGTH], storage_scheme[LINE_LENGTH];
    char* p;
    bool b_getValue = true;

    getline(in, line);
    strcpy(data, line.c_str());
    if (sscanf(data, "%s %s %s %s %s", banner, mtx, crd, data_type, storage_scheme) != 5) {
        cout << "ERROR: mtx header is missing" << endl;
        return 0;
    }

    for (p=data_type; *p!='\0'; *p=tolower(*p),p++);

    if (strcmp(data_type, "pattern") == 0)  b_getValue = false;

    getline(in, line);
    while(line.size()>0&&line[0]=='%') getline(in,line);
    in2.str(line);
    in2 >> row >> col >> num_of_entries;

    if(row!=col) {
        cout<<"* ERROR: This is not a square matrix."<<endl;
        return 0;
    }

    while(!in.eof() && entry_counter<num_of_entries) {
        getline(in,line);
        entry_counter++;

        if(line!="") {
            in2.clear();
            in2.str(line);
            in2 >> ridx >> cidx >> value;
            ridx--;
            cidx--;

            if (ridx < 0 || ridx >= row)  cout << "sym-mtx error: " << ridx << " row " << row << endl;
            if (cidx < 0 || cidx >= col)  cout << "sym-mtx error: " << cidx << " col " << col << endl;
            if (ridx == cidx)  continue;

            if (ridx > cidx) {
                if (b_getValue) {
                    if(value > connStrength) {
                        v_map[ridx].push_back(cidx);
                        v_map[cidx].push_back(ridx);
                        if (is_gstats) {
                            // e_v and e_u are used for triangle cores.
                            // both must be |E|/2, since undirected. Must only have unique edges.
                            e_v.push_back(ridx);
                            e_u.push_back(cidx);
                        }
                    }
                } else {
                    v_map[ridx].push_back(cidx);
                    v_map[cidx].push_back(ridx);
                    if (is_gstats) {
                        e_v.push_back(ridx);
                        e_u.push_back(cidx);
                    }
                }

                if (b_getValue && value > connStrength) {
                    valueList[ridx].push_back(value);
                    valueList[cidx].push_back(value);
                }
            } else {
                cout << "* WARNING: Found a nonzero in the upper triangular. ";
                break;
            }
        }
    }
    vertices.push_back(edges.size());
    for (int i=0;i < row; i++) {
        edges.insert(edges.end(),v_map[i].begin(),v_map[i].end());
        vertices.push_back(edges.size());
    }
    v_map.clear();
    valueList.clear();
    vertex_degrees();
    return 1;
}

int graphpack_graph::read_metis(const string& filename) { return 0; };






void graphpack_graph::create_adj() {
    is_dense = true;
    double sec = get_time();

    int size = num_vertices();
    adj = new bool*[size];
    for (int i = 0; i < size; i++) {
        adj[i] = new bool[size];
        memset(adj[i], 0, size * sizeof(bool));
    }

    for (int i = 0; i < num_vertices(); i++) {
        for (long long j = vertices[i]; j < vertices[i + 1]; j++ )
            adj[i][edges[j]] = true;
    }
    cout << "Created adjacency matrix in " << get_time() - sec << " seconds" <<endl;
}

graphpack_graph::~graphpack_graph() { }

graphpack_graph::graphpack_graph(const string& filename) {
    initialize();
    fn = filename;
    read_graph(filename);
}

graphpack_graph::graphpack_graph(bool graph_stats, const string& filename) {
    initialize();
    fn = filename;
    is_gstats = graph_stats; // construct e_u and e_v
    read_graph(filename);
}

graphpack_graph::graphpack_graph(const string& filename, bool make_adj) {
    initialize();
    fn = filename;
    read_graph(filename);
    if (make_adj && adj == NULL)
        create_adj();
}

/**
 * @brief       adapt graph representation based on some simple statistics
 *              and memory requirements. This is meant for performance and
 *              flexibility.
 *
 * @todo        remove this from the user, simply decide when graph is read
 *              or created via other graph operations
 *
 * @param adj_limit
 */
void graphpack_graph::optimize_graph_ops(int adj_limit) {
    if (num_vertices() < adj_limit || density() > 0.9) { // small or dense graph
        is_dense = true;
        create_adj();
        printf("[graphpack: DENSE graph detected]  optimizing internal data structures \n");
    }
    else {
        is_dense = false;
        printf("[graphpack: SPARSE graph detected]  optimizing internal data structures \n");
    }
}


/**
 *
 * @param  U is the vertex set to be induced from G
 * @param  G is the initial graph
 * @return H is the induced graph, containing all vertices in vs and all edges between them.
 */
graphpack_graph::graphpack_graph(
        vector<int>& U,
        graphpack_graph& G) {
    initialize();
    // induce H
    cout << "done" <<endl;
}

/*
 * We assume X is initially empty...
 */

graphpack_graph::graphpack_graph(
        vector<Vertex>& U,
        vector<int> & ind, // v --> v'  (vold to vnew), so ind[]
        vector<long long>& vs,
        vector<int>& es) {
    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);
#ifndef DEBUG
    cout << "[graphpack before induced graph]  |edges| = " << edges.size() << ", |es| = " << es.size() <<endl;
#endif
    /**
     * Suppose edge (u,v) where vertex u is in U,
     * then ind is used to check if vertex v is also in the vertex set U
     * If so, then we add that edge to the graph!
     */
//    int sum = 0;
    // maps new vertex_id in H to vertex_id in G
    // only O(|V(H)|)
    vertex_lookup.resize(U.size());
//    vector<int> ind(vs.size(),0);

    // @TODO, ind, reverse_lookup, may need to be saved!
    // perhaps, we should save ind, instead of vertex_lookup to save space
    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
//        sum += degree[U[i].get_id()];
    }
//    edges.reserve(sum);


    /* use U as an explicit relabeling of vertices, i --> U[i]
    /* i is the id in H, and U[i] is the vertex_id in G
     * Complexity: O(\sum_{w \in N(v)} |N(w)|) - has to go over all neighbors of v,
     *    and all their edges, to create the induced graph containing only the edges between N(v)
     * */
    int start = 0;
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id(); // v is a pointer into the edge array (neighbors)
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            // ind maps: old_vertex_id --> new_vertex_id,
            // so if 5 is old id, then ind[5] returns the new_id + 1..
            // note that e_u and e_v are mainly needed for computing triangle cores fast, see is_gstats flag.
            if (ind[es[j]]) {                       // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
//        cout << "vertex: " << v << ", edges:" << edges.size() <<endl;
    }
//    edges.swap(edges);
    vertex_degrees(); // O(|V|) initialize the degrees of H

#ifndef DEBUG
    cout << "[graphpack induced graph]  |E| = " << edges.size() <<endl;
    cout << "[graphpack induced graph]  d_avg(H) = " << avg_degree << ", d_max(H) = " << max_degree <<endl; // induced graph
#endif

    // setup e_u and e_v, also creates fast adj on the induced graph O(1)
    // idea is that the neighborhood will likely be dense anyway
    init_edge_data_structs();
    ind.clear();
}

/**
 * Used with triangle cores
 * @param U
 * @param vs
 * @param es
 */
graphpack_graph::graphpack_graph(
        vector<int>& U,
        vector<long long>& vs,
        vector<int>& es) {
    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);
#ifndef DEBUG
    cout << "[graphpack before induced graph]  |edges| = " << edges.size() << ", |es| = " << es.size() <<endl;
#endif

    /**
     * Suppose edge (u,v) where vertex u is in U,
     * then ind is used to check if vertex v is also in the vertex set U
     * If so, then we add that edge to the graph!
     */
    // maps new vertex_id in H to vertex_id in G
    // only O(|V(H)|)
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);
    for (int i = 0; i < U.size(); i++) {
        ind[U[i]] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i]; // v' --> v
    }


    /* use U as an explicit relabeling of vertices, i --> U[i]
    /* i is the id in H, and U[i] is the vertex_id in G
     * Complexity: O(\sum_{w \in N(v)} |N(w)|) - has to go over all neighbors of v,
     *    and all their edges, to create the induced graph containing only the edges between N(v)
     * */
    int start = 0;
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
//        v = U[i].get_id(); // v is a pointer into the edge array (neighbors)
        for (long long j = vs[U[i]]; j < vs[U[i] + 1]; j++ ) {
            // ind maps: old_vertex_id --> new_vertex_id,
            // so if 5 is old id, then ind[5] returns the new_id + 1..
            // note that e_u and e_v are mainly needed for computing triangle cores fast, see is_gstats flag.
            if (ind[es[j]]) {                       // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }

            // @TODO: can we fix e_u, e_v here...? or just compute afterwards? note eid is computed in tcores
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
//        cout << "vertex: " << v << ", edges:" << edges.size() <<endl;
    }
    vertex_degrees(); // O(|V|) initialize the degrees of H

#ifndef DEBUG
    cout << "[graphpack induced graph]  |E| = " << edges.size() <<endl;
    cout << "[graphpack induced graph]  d_avg(H) = " << avg_degree << ", d_max(H) = " << max_degree <<endl; // induced graph
#endif

    if (U.size() > 10000) {
        // setup e_u and e_v
        init_edge_data_structs_no_adj();
    }
    else {
        // setup e_u and e_v, also creates fast adj on the induced graph O(1)
        // idea is that the neighborhood will likely be dense anyway
        init_edge_data_structs();
    }
    ind.clear();
}


graphpack_graph::graphpack_graph(
        vector<Vertex>& U,
        vector<long long>& vs,
        vector<int>& es) {
    initialize();

    cout << "[graphpack - constructing subgraph]  induced subgraph has " << U.size() << " vertices" <<endl;

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);

    /**
     * Suppose edge (u,v) where vertex u is in U,
     * then ind is used to check if vertex v is also in the vertex set U
     * If so, then we add that edge to the graph!
     */
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);

    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
    }


    /* use U as an explicit relabeling of vertices, i --> U[i]
    /* i is the id in H, and U[i] is the vertex_id in G
     * Complexity: O(\sum_{w \in N(v)} |N(w)|) - has to go over all neighbors of v,
     *    and all their edges, to create the induced graph containing only the edges between N(v)
     * */
    int start = 0;
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id();
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            if (ind[es[j]]) {                       // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
    }
    vertex_degrees(); // O(|V|) initialize the degrees of H

#ifndef DEBUG
    cout << "[graphpack induced graph]  |E| = " << edges.size() <<endl;
    cout << "[graphpack induced graph]  d_avg(H) = " << avg_degree << ", d_max(H) = " << max_degree <<endl; // induced graph
#endif

    // setup e_u and e_v, also creates fast adj on the induced graph O(1)
    // idea is that the neighborhood will likely be dense anyway
    if (U.size() > 15000)
        init_edge_data_structs_no_adj();
    else
        init_edge_data_structs();
    ind.clear();
}



/*
 * Note that this version, updates the vertex IDs in U, so that _ind_ doesn't need to be saved.
 *
 */
graphpack_graph::graphpack_graph(
        vector<Vertex>& U,
        vector<long long>& vs,
        vector<int>& es,
        bool is_U_relabeled) {
    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);
#ifndef DEBUG
    cout << "[graphpack before induced graph]  |edges| = " << edges.size() << ", |es| = " << es.size() <<endl;
#endif
    /**
     * Suppose edge (u,v) where vertex u is in U,
     * then ind is used to check if vertex v is also in the vertex set U
     * If so, then we add that edge to the graph!
     */
//    int sum = 0;
    // maps new vertex_id in H to vertex_id in G
    // only O(|V(H)|)
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);

    // @TODO, ind, reverse_lookup, may need to be saved!
    // perhaps, we should save ind, instead of vertex_lookup to save space
    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
//        sum += degree[U[i].get_id()];
    }
//    edges.reserve(sum);


    /* use U as an explicit relabeling of vertices, i --> U[i]
    /* i is the id in H, and U[i] is the vertex_id in G
     * Complexity: O(\sum_{w \in N(v)} |N(w)|) - has to go over all neighbors of v,
     *    and all their edges, to create the induced graph containing only the edges between N(v)
     * */
    int start = 0;
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id(); // v is a pointer into the edge array (neighbors)
        // @TODO: CHECK THIS, MIGHT BE BETTER OPTION
        U[i].set_id(i);//ind[i]-1); // set P back to..
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            // ind maps: old_vertex_id --> new_vertex_id,
            // so if 5 is old id, then ind[5] returns the new_id + 1..
            // note that e_u and e_v are mainly needed for computing triangle cores fast, see is_gstats flag.
            if (ind[es[j]]) {                       // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
//        cout << "vertex: " << v << ", edges:" << edges.size() <<endl;
    }
//    edges.swap(edges);
    vertex_degrees(); // O(|V|) initialize the degrees of H

#ifndef DEBUG
    cout << "[graphpack induced graph]  |E| = " << edges.size() <<endl;
    cout << "[graphpack induced graph]  d_avg(H) = " << avg_degree << ", d_max(H) = " << max_degree <<endl; // induced graph
#endif

    // setup e_u and e_v, also creates fast adj on the induced graph O(1)
    // idea is that the neighborhood will likely be dense anyway
    init_edge_data_structs();
    ind.clear();
}



void graphpack_graph::init_edge_data_structs_no_adj() {
    // set to be edges.size()/2, unique edges..
    // use the fact that if u is smaller than v, then its a duplicate!
    // since we construct e_u and e_v in order...
    // @TODO: create its own function, and use adj matrix on neighborhood,
    //        since the neighborhood will likely be dense anyway.


//    if (adj == NULL) { // function should only be called if adj is not already set...
        int size = num_vertices();
//        adj = new bool*[size];
//        for (int i = 0; i < size; i++) {
//            adj[i] = new bool[size];
//            memset(adj[i], 0, size * sizeof(bool));
//        }

        // @TODO: we may also resize both e_v and e_u upfront,
        //      then use an edge_idx to keep track of edge pos so far.
        // Complexity: Only O(|E(H)|) as compared to the previous time.
        int u = 0;
        for (int v = 0; v < num_vertices(); v++) {
            for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
                u = edges[j];
                if (u > v) {
//                if (!adj[u][v]) { // check if u--v edge does NOT exists,
                    // note that if u > v, then edge was already added..
                    e_v.push_back(v);
                    e_u.push_back(u);
                    //                eid[j] = e_v.size(); // only works for half... ughh
//                }
                }
//                adj[v][u] = true;
            }
//        }
#ifndef DEBUG
        cout << "e_v.size() = " << e_v.size() << ", e_u.size() = " << e_u.size() <<endl;
#endif
    }
}


void graphpack_graph::init_edge_data_structs() {
    // set to be edges.size()/2, unique edges..
    // use the fact that if u is smaller than v, then its a duplicate!
    // since we construct e_u and e_v in order...
    // @TODO: create its own function, and use adj matrix on neighborhood,
    //        since the neighborhood will likely be dense anyway.


    if (adj == NULL) { // function should only be called if adj is not already set...
        int size = num_vertices();
        adj = new bool*[size];
        for (int i = 0; i < size; i++) {
            adj[i] = new bool[size];
            memset(adj[i], 0, size * sizeof(bool));
        }

        // @TODO: we may also resize both e_v and e_u upfront,
        //      then use an edge_idx to keep track of edge pos so far.
        // Complexity: Only O(|E(H)|) as compared to the previous time.
        int u = 0;
        for (int v = 0; v < num_vertices(); v++) {
            for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
                u = edges[j];
                if (!adj[u][v]) { // check if u--v edge does NOT exists,
                    // note that if u > v, then edge was already added..
                    e_v.push_back(v);
                    e_u.push_back(u);
                    //                eid[j] = e_v.size(); // only works for half... ughh
                }
                adj[v][u] = true;
            }
        }
#ifndef DEBUG
        cout << "e_v.size() = " << e_v.size() << ", e_u.size() = " << e_u.size() <<endl;
#endif
    }
}


/*
 * Afterwards, U may be less than adj, or vs, etc...
 * The vertices in U may also be in various order
 */
graphpack_graph::graphpack_graph(
        vector<Vertex>& U, vector<long long>& vs, vector<int>& es, int & mc, bool is_tcore = false) {

    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);

    // maps new vertex_id in H to vertex_id in G
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);
    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
    }

    degree.resize(U.size(),0);
    int start = 0, valid_verts = 0, n = U.size();
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id();                          // v is a pointer into the edge array (neighbors)
        U[i].set_id(i);
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            if (ind[es[j]]) {                      // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        if ((edges.size() - start) > mc) {
            degree[i] = edges.size() - start;
            // put valid verts upfront in same order as they appear
            U[valid_verts].set_id(i);
            U[valid_verts].set_bound(degree[i]);
            valid_verts++;
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
    }
    // remove pruned verts
    for (int i = valid_verts; i < n; i++)
        U.pop_back();

    /* setup e_u and e_v, also creates fast adj on the induced graph O(1)
     * idea is that the neighborhood will likely be dense anyway
     **/
    if (is_tcore) {
        init_data_structs(mc);
    }
    else {
        init_edge_data_structs_pruned(mc,U);
    }
    ind.clear();
}


/**
 * @brief: Given G=(V,E), induce H(V-C), let W = V - C,
 * where V-C = {u : u \in V \text{and} u \not \in C}
 *
 * C is the set of vertices to remove from vs and es
 *
 * Output is H,
 * Also creates,
 * vertices array
 * edges array
 * vertex_lookup
 * degrees
 */
graphpack_graph::graphpack_graph(graphpack_graph &G, vector<int>& C) { //vector<long long>& vs, vector<int>& es) {

    initialize();


    int n_prev = G.vertices.size()-1;
    int n = (G.vertices.size()-1) - C.size();
    vertices.resize(n); // |V_new| = |V| - |C|

    // maps new vertex_id in H to vertex_id in G
    vertex_lookup.resize(n);
    vector<int> ind(n_prev,0);
//    for (int i = 0; i < C.size(); i++) {
//        ind[C[i].get_id()] = i + 1;       // v  --> v'
//        vertex_lookup[i] = C[i].get_id(); // v' --> v
//    }

    // hash the vertices to be removed
    vector<int> pruned(n_prev,0);



//    if (G.vertex_lookup.size() == 0) {
//        printf("vertex_lookup does not exist! \n");
//        for (int i = 0; i < C.size(); ++i) {
//            pruned[C[i]] = 1;
//        }
//    }
//    else { // assuming vertices in C are not mapped back to their original ids
//        printf("vertex_lookup exists! \n");
//        for (int i = 0; i < C.size(); ++i) {
//
//            int u = G.vertex_lookup[C[i]];
//            pruned[u] = 1;
//        }
//    }

    print_line(100);
    printf("n = %d, n_prev = %d \n", n, n_prev);

    // vertices in C are from G.vertices/G.edges, hence, consistent! map ids later!
    for (int i = 0; i < C.size(); ++i) {
        pruned[C[i]+1] = 1;
        printf("\t pruned v=%d \n",C[i]+1);
    }

    printf("created the pruned array \n");

    int vertex_id = 0;
    for (int i = 0; i < n_prev; ++i) {
        if (pruned[i+1] == 0) { // not pruned or marked to be removed!
            // v_prev  --> v_new (new vertex id),
            // example: ind[vid] returns new vertex id
            ind[i] = vertex_id + 1;

            // v_new --> v_prev  (store original vertex id)
            // vertex_lookup[i] returns previous vertex id
//            if (G.vertex_lookup.size() == 0) {
//
//            }
//            else {
//
//            }
            vertex_lookup[vertex_id] = i;
            vertex_id++;
        }
    }

//    if (false) {
//        if (G.vertex_lookup.size() == 0) {
//            printf("\t graph has not been induced before \n");
//            // ids of the vertices are the indices (graph has not been induced before)
//            for (int i = 0; i < n_prev; ++i) {
//                if (pruned[i] == 0) { // not pruned or marked to be removed!
//
//                    // v_prev  --> v_new (new vertex id),
//                    // example: ind[vid] returns new vertex id
//                    //                ind[G.vertex_lookup[i]] = i + 1;
//                    ind[i] = i + 1;
//
//                    // v_new --> v_prev  (store original vertex id)
//                    // vertex_lookup[i] returns previous vertex id
//                    vertex_lookup[i] = i;
//                }//if pruned
//            }//for
//        }// no prev induced graph! (new)
//        else {
//            printf("\t induced previously! need to copy over the ids from the previous graphs vertex_lookup \n");
//            // need to copy over the ids from the previous graphs vertex_lookup (if they exist)
//            // graph given as params, was previously induced!
//            for (int i = 0; i < n_prev; ++i) {
//                if (pruned[i] == 0) { // not pruned or marked to be removed!
//
//                    // v_prev  --> v_new (new vertex id),
//                    // example: ind[vid] returns new vertex id
//                    ind[G.vertex_lookup[i]] = i + 1;
//
//                    // v_new --> v_prev  (store original vertex id)
//                    // vertex_lookup[i] returns previous vertex id
//                    vertex_lookup[i] = G.vertex_lookup[i];
//                } // if
//            }// for each vertex
//        } // prev induced graph
//    }

    printf("created ind and vertex_lookup \n");



    print_line(50);
    degree.resize(n,0);
    int start = 0;
    for (int i = 0; i < n; ++i) {
        start = edges.size();
        int v = vertex_lookup[i]; // returns previous vertex id from G

        printf("\t i=%d is v=%d, start = %d \n", i,v,start);

        // for each $u \in N(v)$, add u if and only if $u \not\in C$
        // this reduces the set of neighbors
        for (long long j = G.vertices[v]; j < G.vertices[v + 1]; j++ ) {
            int u = G.edges[j];
            if (ind[u]) {                      // returns new vertex id. (v has an edge with u, and u is also in U)
                edges.push_back(ind[u]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
    }


//    for (int v = 0; v < n_prev; ++v) {
//        // if pruned, then all v's edges are ignored/discarded!
//        if (pruned[v] == 0) { // add v to H if not removed (i.e., in the set C)
//            start = edges.size();
//
//            // for each $u \in N(v)$, add u if and only if $u \not\in C$
//            // this reduces the set of neighbors
//            for (long long j = G.vertices[v]; j < G.vertices[v + 1]; j++ ) {
//                int u = G.edges[j];
//                if (ind[u]) {                      // v has an edge with u, and u is also in U
//                    edges.push_back(ind[u]-1);     // map the id of u in G, to the new vertex id in U
//                }
//            }
//            int id = ind[v]-1;
//            vertices[id] = start;         // set pointers for vertex i to point to its edges */
//            vertices[id + 1] = edges.size();
//        }
//    }

    print_vertices_array();
    print_edges();

    printf("finished inducing the graph: removing the vertices in C from the graph! \n");

    if (G.is_dense) {
        // update adj
        // or recreate adj!
    }

    printf("finished creating new adj matrix \n");

    ind.clear();
}




graphpack_graph::graphpack_graph(
        vector<Vertex>& U, vector<long long>& vs, vector<int>& es, int & mc, bool is_tcore, bool is_enum) {

    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);

    // maps new vertex_id in H to vertex_id in G
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);
    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
    }

    degree.resize(U.size(),0);
    int start = 0, valid_verts = 0, n = U.size();
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id();                          // v is a pointer into the edge array (neighbors)
        U[i].set_id(i);
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            if (ind[es[j]]) {                      // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        if (is_enum) {
            if ((edges.size() - start) >= mc - 1) {
                degree[i] = edges.size() - start;
                // put valid verts upfront in same order as they appear
                U[valid_verts].set_id(i);
                U[valid_verts].set_bound(degree[i]);
                valid_verts++;
            }
        }
        else if ((edges.size() - start) > mc - 1) {
            degree[i] = edges.size() - start;
            // put valid verts upfront in same order as they appear
            U[valid_verts].set_id(i);
            U[valid_verts].set_bound(degree[i]);
            valid_verts++;
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
    }
    // remove pruned verts
    for (int i = valid_verts; i < n; i++)
        U.pop_back();

    /* setup e_u and e_v, also creates fast adj on the induced graph O(1)
     * idea is that the neighborhood will likely be dense anyway
     **/
    if (is_tcore) {
        init_data_structs(mc,is_enum);
    }
    else {
        init_edge_data_structs_pruned(mc,U,is_enum);
    }
    ind.clear();
}




void graphpack_graph::init_data_structs(int & mc, bool is_enum) { // = false) {

    if (adj == NULL) { // function should only be called if adj is not already set...
        int size = num_vertices();
        adj = new bool*[size];
        for (int i = 0; i < size; i++) {
            adj[i] = new bool[size];
            memset(adj[i], 0, size * sizeof(bool));
        }

        //! Complexity: Only O(|E(H)|) as compared to the previous time.
        int u = 0;
        for (int v = 0; v < num_vertices(); v++) {
            for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
                u = edges[j];
                if (!adj[u][v]) { // check if u--v edge does NOT exists,
                    e_v.push_back(v);
                    e_u.push_back(u);
                }
                if (is_enum && degree[u] >= mc - 1) // we use this in the if above
                    adj[v][u] = true;
                else if (degree[u] > mc - 1) // we use this in the if above
                    adj[v][u] = true;
            }
        }
    }
}



/*
 * Afterwards, U may be less than adj, or vs, etc...
 * The vertices in U may also be in various order
 */
graphpack_graph::graphpack_graph(
        vector<Vertex>& U, vector<long long>& vs, vector<int>& es, int & mc) {

    initialize();

    // data structures for the new induced graph H
    vertices.resize(U.size()+1,0);

    // maps new vertex_id in H to vertex_id in G
    vertex_lookup.resize(U.size());
    vector<int> ind(vs.size(),0);
    for (int i = 0; i < U.size(); i++) {
        ind[U[i].get_id()] = i + 1;       // v  --> v'
        vertex_lookup[i] = U[i].get_id(); // v' --> v
    }


    /* use U as an explicit relabeling of vertices, i --> U[i]
    /* i is the id in H, and U[i] is the vertex_id in G
     * Complexity: O(\sum_{w \in N(v)} |N(w)|) - has to go over all neighbors of v,
     *    and all their edges, to create the induced graph containing only the edges between N(v)
     * */
    degree.resize(U.size(),0); // compute degree on the fly
    int start = 0, valid_verts = 0, n = U.size();
    long long v = 0;
    for (int i = 0; i < U.size(); i++) {
        start = edges.size();
        v = U[i].get_id(); // v is a pointer into the edge array (neighbors)
        U[i].set_id(i);
        for (long long j = vs[v]; j < vs[v + 1]; j++ ) {
            if (ind[es[j]]) {                      // v has an edge with u, and u is also in U
                edges.push_back(ind[es[j]]-1);     // map the id of u in G, to the new vertex id in U
            }
        }
        if ((edges.size() - start) >= mc) {
            degree[i] = edges.size() - start;
            // put valid verts upfront in same order as they appear
            U[valid_verts].set_id(i);
//            U[valid_verts].set_bound(degree[i]);
            valid_verts++;
        }
        vertices[i] = start;         // set pointers for vertex i to point to its edges */
        vertices[i + 1] = edges.size();
    }
    // remove pruned verts
    for (int i = valid_verts; i < n; i++)
        U.pop_back();

    /* setup e_u and e_v, also creates fast adj on the induced graph O(1)
     * idea is that the neighborhood will likely be dense anyway
     **/
    init_edge_data_structs_pruned(mc,U);
    ind.clear();
}


void graphpack_graph::init_edge_data_structs_pruned(int & mc, vector<Vertex> & U, bool is_enum) {


    if (adj == NULL) { // function should only be called if adj is not already set...
        int size = num_vertices();
        adj = new bool*[size];
        for (int i = 0; i < size; i++) {
            adj[i] = new bool[size];
            memset(adj[i], 0, size * sizeof(bool));
        }

        // Complexity: Only O(|E(H)|) as compared to the previous time.
        int u = 0;
        for (int i = 0; i < U.size(); i++) {
            int v = U[i].get_id();
            int sum = 0;
            for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
                u = edges[j];
                sum += degree[u];

                if (is_enum && degree[u] >= mc - 1) // we use this in the if above
                    adj[v][u] = true;
                else if (degree[u] > mc - 1) // we use this in the if above
                    adj[v][u] = true;

//                if (degree[u] >= mc) // we use this in the if above
//                    adj[v][u] = true;
//                else degree[u] = 0;
            }
            U[i].set_bound(sum);
        }

//        for (int v = 0; v < num_vertices(); v++) {
//            int sum = 0;
//            for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
//                u = edges[j];
//                sum += degree[u];
//
//                if (degree[u] >= mc) // we use this in the if above
//                    adj[v][u] = true;
//            }
//        }
    }
}

/*
 * Remove adj from this func, and move it above... since we may want to ignore this step on some graphs
 */
void graphpack_graph::init_edge_data_structs_pruned(int & mc, bool is_enum) {
    // set to be edges.size()/2, unique edges..
    // use the fact that if u is smaller than v, then its a duplicate!
    // since we construct e_u and e_v in order...
    // @TODO: create its own function, and use adj matrix on neighborhood,
    //        since the neighborhood will likely be dense anyway.


    if (adj == NULL) { // function should only be called if adj is not already set...
        int size = num_vertices();
        adj = new bool*[size];
        for (int i = 0; i < size; i++) {
            adj[i] = new bool[size];
            memset(adj[i], 0, size * sizeof(bool));
        }

        // @TODO: we may also resize both e_v and e_u upfront,
        //      then use an edge_idx to keep track of edge pos so far.
        // Complexity: Only O(|E(H)|) as compared to the previous time.
        int u = 0;
        for (int v = 0; v < num_vertices(); v++) {
//            if (degree[v] > mc) { // @TODO: FIX THIS... WE NEED e_v.push_back() for tricores..
                for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
                    u = edges[j];
                    if (!adj[u][v]) { // check if u--v edge does NOT exists,
                        // note that if u > v, then edge was already added..
                        e_v.push_back(v);
                        e_u.push_back(u);
                        //                eid[j] = e_v.size(); // only works for half... ughh
                    }
                    if (degree[u] >= mc) // we use this in the if above
                        adj[v][u] = true;
                }
//            }
        }
#ifndef DEBUG
        cout << "e_v.size() = " << e_v.size() << ", e_u.size() = " << e_u.size() <<endl;
#endif
    }
}



// afterwards, order verts by dual_deg
void graphpack_graph::triangle_core_pruning(vector<Vertex> & P, int & lb) {


//    vector<long long> VS(vertices.size(),0);
    vector<int> E;
    E.reserve(edges.size());


    vector<Vertex> V;
    V.reserve(num_vertices());

    long long start = 0;
    long long edge_pos = 0;

    /*
     * Removes edges with Triangle Core less than K
     * After removing the edges with triangle core less than K,
     * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
     */

    /**
     * Note: do not need E or VS, just adj
     */
//    int num_edges_pruned = 0;
//    int num_verts_pruned = 0;
//    if (lb == max_tri_core+2) {
//        cout << "inside tricore pruning, lb == max tri core..." <<endl;
//        for (int i = 0; i < num_vertices(); i++) {
//            start = E.size();
//            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
//
//                edge_pos = eid[j];// + 1;              // tri_core is m/2, so this finds the right tri_core value given an edge id
//                if (tri_core[edge_pos+1]+2 >= lb) { // get correct tri_core for edge's eid
//                    E.push_back(edges[j]);
//                }
//                else {
//                    adj[e_u[edge_pos]][e_v[edge_pos]] = false;
//                    adj[e_v[edge_pos]][e_u[edge_pos]] = false;
//                    num_edges_pruned++;
//                }
//            }
//            if ((E.size() - start) >= lb-1) { // basically the updated degree bound
//                V.push_back(Vertex(i,(E.size() - start)));
//            } else num_verts_pruned++;
////            VS[i] = start;
////            VS[i + 1] = E.size();
//        }
//    }
//    else {
//    int sum_tcore = 0;
    for (int i = 0; i < num_vertices(); i++) {
        start = E.size();
        for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {

            // add 1 to edge_pos for tri_core, but not e_u and e_v
            edge_pos = eid[j];// + 1;
            if (tri_core[edge_pos+1]+2 > lb) { // get correct tri_core for edge's eid
                E.push_back(edges[j]);
//                sum_tcore += tri_core[edge_pos+1]+2;
            }
            else {
                adj[e_u[edge_pos]][e_v[edge_pos]] = false;
                adj[e_v[edge_pos]][e_u[edge_pos]] = false;
//                num_edges_pruned++;
            }
        }

        // add vert if deg >= lb
        if ((E.size() - start) >= lb) { // basically the updated degree bound
            V.push_back(Vertex(i,(E.size() - start)));
        }
//        else num_verts_pruned++;


        //            VS[i] = start;
        //            VS[i + 1] = E.size();
    }
//    }
//    cout << "after triangle core pruning..." <<endl;
//
//    cout << "num_verts pruned = " << num_verts_pruned << " out of " << num_vertices() <<endl;
//    cout << "num edges pruned = " << num_edges_pruned << " out of " << edges.size() <<endl;

//    std::sort(V.begin(), V.end(), decr_bound);
    P = V;

}



// afterwards, order verts by dual_deg
void graphpack_graph::triangle_core_pruning(vector<Vertex> & P, int & lb, string neigh_ordering) {


//    vector<long long> VS(vertices.size(),0);
    vector<int> E;
    E.reserve(edges.size());


    vector<Vertex> V;
    V.reserve(num_vertices());

    long long start = 0;
    long long edge_pos = 0;

    /*
     * @brief Removes edges with Triangle Core less than K
     * After removing the edges with triangle core less than K,
     * We check if the number of edges remaining is at least K, and if not, we remove the vertex.
     */

    /**
     * Note: do not need E or VS, just adj
     */
//    int num_edges_pruned = 0;
//    int num_verts_pruned = 0;

    int sum_tcore = 0;
    for (int i = 0; i < num_vertices(); i++) {

        sum_tcore = 0;
        start = E.size();
        for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {

            // add 1 to edge_pos for tri_core, but not e_u and e_v
            edge_pos = eid[j];// + 1;
            if (tri_core[edge_pos+1]+2 > lb) { // get correct tri_core for edge's eid
                E.push_back(edges[j]);
                sum_tcore += tri_core[edge_pos+1]+2;
            }
            else {
                adj[e_u[edge_pos]][e_v[edge_pos]] = false;
                adj[e_v[edge_pos]][e_u[edge_pos]] = false;
            }
        }

        // add vert if deg >= lb
        if ((E.size() - start) >= lb) { // basically the updated degree bound
            V.push_back(Vertex(i,sum_tcore));
//            V.push_back(Vertex(i,(E.size() - start)));
        }
    }
    P = V;

}





void graphpack_graph::tcore_pruner(graphpack_graph & H, int & lb, vector<int>& U_new) {
    vector<long long> VS(H.vertices.size(),0);
    vector<int> E;
    E.reserve(H.edges.size());

    // perhaps the best thing to do is to create a set of vertices to keep,
    // and then ship this off to the appropriate constructor
    // which would reduce the graph. even if a bit of work is overlapping.
    //
    // at this point, we don't know how many edges or vertices will be remaining
    // also if we prune with degree, then we can't add edges initially, since they may be removed..
    cout << "[TRIANGLE CORE PRUNING]  |E| = " << H.num_edges();
    cout << ", H.edges.size() = " << H.edges.size() << " (before)" <<endl;
    int vertex_id = 0;
    int start = 0;
    for (int i = 0; i < H.num_vertices(); i++) {
        start = E.size();
        // how can we  map the edges if we don't know the vertices?
        // keep an ind, and each time
        for (long long j = H.vertices[i]; j < H.vertices[i + 1]; j++ ) {
            if (H.tri_core[H.edges[j]]+2 >= lb)
                E.push_back(H.edges[j]); // at this point, we know
        }
        if (start != E.size()) { // no edges to vertex
            if ((start - E.size()) >= lb) { // basically the updated degree bound, but what about the edges???
                U_new.push_back(i);
                vertex_lookup.push_back(i);
                // can we add the vertex_lookup check here????
            }
        }
        VS[i] = start;
        VS[i + 1] = E.size();
    }
    vertex_degrees(); // O(|V|) initialize the degrees of H

 #ifndef DEBUG
     cout << "[graphpack induced graph]  |E| = " << edges.size() <<endl;
     cout << "[graphpack induced graph]  d_avg(H) = " << avg_degree;
     cout << ", d_max(H) = " << max_degree <<endl; // induced graph
 #endif

     // setup e_u and e_v, also creates fast adj on the induced graph O(1)
     // idea is that the neighborhood will likely be dense anyway
     init_edge_data_structs();
}



void graphpack_graph::triangle_core_pruning_only_adj_update_degree(vector<Vertex> & P, int & mc) {
    int valid = 0;
    int prune_vert = 0;
    long long edge_pos = 0;

    int n = P.size();
    int v = 0;
    for (int i = 0; i < n; i++) {

        v = P[i].get_id(); // get the actual vertex id in position P[i]
        for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {

            edge_pos = eid[j] + 1; // tri_core is m/2, so this finds the right tri_core value given an edge id

            // keep edges that may lead us to find a clique of _strictly a larger size_
            if (tri_core[edge_pos]+2 < mc) { // remove edge
                adj[v][edges[j]] = false;
                degree[v]--;
            }
        }

        if (degree[v] >= mc) { // clique strictly larger than mc
            P[valid].set_id(v);
            valid++;
        }
        else prune_vert++;
    }

    // remove pruned verts from P
    for (int i = 0; i < prune_vert; i++)
        P.pop_back();
}



/*
 * Could update degree on the fly.
 * If we updated degree, then we could
 *  construct/update P very quickly by checking degree
 */
void graphpack_graph::triangle_core_pruning_only_adj(int & mc) {

    long long edge_pos = 0;
    for (int v = 0; v < num_vertices(); v++) {
        for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {

            edge_pos = eid[j] + 1; // tri_core is m/2, so this finds the right tri_core value given an edge id

            // keep edges that may lead us to find a clique of _strictly a larger size_
            if (tri_core[edge_pos]+2 <= mc) { // remove edge
                adj[v][edges[j]] = false;
            }
        }
    }
}






void graphpack_graph::sum_vertex_degrees() {
    int n = vertices.size() - 1;

    uint64_t sum = 0;
    for (long long v = 0; v < n; v++) {
        degree[v] = vertices[v+1] - vertices[v];
        sum += (degree[v] * degree[v]-1) / 2;
    }
    cout << "sum of degrees: " << sum <<endl;
}



void graphpack_graph::vertex_degrees() {
    int n = vertices.size() - 1;
    degree.resize(n);

    // initialize min and max to degree of first vertex
    max_degree = min_degree = vertices[1] - vertices[0];
    for (long long v=0; v<n; v++) {
        degree[v] = vertices[v+1] - vertices[v];
        if (max_degree < degree[v])  max_degree = degree[v];
        if (degree[v] < min_degree)  min_degree = degree[v];
    }
    avg_degree = (double)edges.size()/n;
    return;
}




// fast update
void graphpack_graph::update_degrees() {
    for (long long v = 0; v < num_vertices(); v++)
        degree[v] = vertices[v+1] - vertices[v];
}

void graphpack_graph::update_degrees(bool flag) {

    int p = 0;
    max_degree = vertices[1] - vertices[0];
    for (long long v = 0; v < num_vertices(); v++) {
        degree[v] = vertices[v+1] - vertices[v];
        if (degree[v] > 0) {
            if (max_degree < degree[v])  max_degree = degree[v];
            p++;
        }
    }
    avg_degree = (double)edges.size() / p;
    return;
}


/**
 * Thread-safe dynamic degree update
 */
void graphpack_graph::update_degrees(vector<long long> & vs, vector<int> & es, int & mc, int* & pruned) {

    vector<int> deg(vs.size()-1,0);
    int num_vs = 0;
    int max_deg = vs[1] - vs[0];
    for (long long v = 0; v < num_vertices(); v++) {
        deg[v] = vs[v+1] - vs[v];
        if (deg[v] > mc) {
            if (max_deg < deg[v]) {
                max_deg = deg[v];
            }
            num_vs++;
        }
        else {
            pruned[v] = 1;
        }
    }
    avg_degree = (double)es.size() / num_vs;
    max_degree = max_deg;
    degree = deg;
    return;
}


void graphpack_graph::update_degrees(int* &pruned, int& mc) {
    max_degree = -1;
    min_degree = std::numeric_limits<int>::max();
    int p = 0;
    for (long long v=0; v < num_vertices(); v++) {
        degree[v] = vertices[v+1] - vertices[v];
        if (degree[v] < mc) {
            if (!pruned[v])  pruned[v] = 1;
            p++;
        }
        else {
            if (max_degree < degree[v])  max_degree = degree[v];
            if (degree[v] < min_degree)  min_degree = degree[v];
        }
    }
    avg_degree = (double)edges.size() / p;
    cout << ", pruned: " << p << endl;
}






void graphpack_graph::basic_stats() {
    cout << "|V|: " << num_vertices() <<endl;
    cout << "|E|: " << num_edges() <<endl;
    cout << "p: " << density() <<endl;
    cout << "d_max: " << get_max_degree() <<endl;
    cout << "d_avg: " << get_avg_degree() <<endl;
}

void graphpack_graph::basic_stats_line() {
    cout << "|V|: " << num_vertices();
    cout << ", |E|: " << num_edges();
    cout << ", p: " << density();
    cout << ", d_max: " << get_max_degree();
    cout << ", d_avg: " << get_avg_degree();
}

void graphpack_graph::basic_stats(string prefix) {
    cout << prefix << "|V|: " << num_vertices() <<endl;
    cout << prefix << "|E|: " << num_edges() <<endl;
    cout << prefix << "p: " << density() <<endl;
    cout << prefix << "d_max: " << get_max_degree() <<endl;
    cout << prefix << "d_avg: " << get_avg_degree() <<endl;
}











string graphpack_graph::get_file_extension(const string& filename) {
    string::size_type result;
    string fileExtension = "";
    result = filename.rfind('.', filename.size() - 1);
    if(result != string::npos)
        fileExtension = filename.substr(result+1);
    return fileExtension;
}



void graphpack_graph::reduce_graph(int* &pruned) {
    vector<long long> V(vertices.size(),0);
    vector<int> E;
    E.reserve(edges.size());

//    if (is_dense_graph()) {
//        int start = 0;
//        for (int i = 0; i < num_vertices(); i++) {
//            start = E.size();
//            if (!pruned[i]) {
//                for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
//                    if (!pruned[edges[j]])
//                        E.push_back(edges[j]);
//                    else {
//                        adj[i][edges[j]] = false;
//                        adj[edges[j]][i] = false;
//                    }
//                }
//            }
//            V[i] = start;
//            V[i + 1] = E.size();
//        }
//    }
//    else {
//
//    }
    int start = 0;
    for (int i = 0; i < num_vertices(); i++) {
        start = E.size();
        if (!pruned[i]) {
            for (long long j = vertices[i]; j < vertices[i + 1]; j++ ) {
                if (!pruned[edges[j]])
                    E.push_back(edges[j]);
            }
        }
        V[i] = start;
        V[i + 1] = E.size();
    }
    vertices = V;
    edges = E;
}


void graphpack_graph::reduce_graph(
        vector<long long>& vs,
        vector<int>& es,
        int* &pruned,
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
    }
    vs = V;
    es = E;
}


int graphpack_graph::coloring_number(string vertex_ordering, bool decr_order) {
    vector<Vertex> V;
    V.reserve(num_vertices());
    order_vertices(V,vertex_ordering,decr_order);
    cout << "|V| = " << V.size() << ", ";

    if (vertex_ordering == "deg" && degree.size() == 0) {
        cout << "[mcpack]  degrees have not been computed." <<endl;
        vertex_degrees(); // compute degrees
    }
    if (vertex_ordering == "kcore" && kcore.size() == 0) {
        cout << "[mcpack]  k-cores have not been computed." <<endl;
        compute_cores(); // compute k-cores
    }

    int coloring_number = color_vertices(V);
    return coloring_number;
}

bool graphpack_graph::conflicts(vector< vector<int> >& colors, int k, int u) {

    for (int i = 0; i < colors[k].size(); i++)
        if (adj[u][colors[k][i]])  return true;
    return false;
}

bool graphpack_graph::conflicts(vector< vector<int> >& colors, int k, vector<short> & ind) {

    for (int i = 0; i < colors[k].size(); i++)
        if (ind[colors[k][i]])  return true;
    return false;
}


int graphpack_graph::color_vertices(vector< Vertex > & V) {

    int j = 0, u = 0, k = 1;
    int max_k = 1;
    int min_k = 3;

    // init the neigh coloring array
    vector< vector<int> > colors(get_max_degree()+1);
    for (int i = 0; i < get_max_degree()+1; i++)
        colors[i].reserve(get_max_degree()+1);
    colors[1].clear(); colors[2].clear();


    if (adj == NULL) { // large sparse graph, use CSC
        cout << "large sparse graph detected, using CSC to color" <<endl;
        vector<short> ind(num_vertices(),0);
        for (int w = 0; w < V.size(); w++) {
            u = V[w].get_id();
            k = 1;

            for (long long h = vertices[u]; h < vertices[u + 1]; h++)  ind[edges[h]] = 1;

            while (conflicts(colors,k,ind)) {
                k++; // try the next color class
            }

            for (long long h = vertices[u]; h < vertices[u + 1]; h++)  ind[edges[h]] = 0;

            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            colors[k].push_back(u);
            if (k < min_k) {
                V[j].set_id(u);
                j++;
            }
        }
    }
    else { // dense graph, use adj for faster lookups

        cout << "dense graph detected, using adj for computing the coloring number" <<endl;
        for (int w=0; w < V.size(); w++) {
            u = V[w].get_id();
            k = 1;

            while (conflicts(colors,k,u)) {
                k++; // try the next color class
            }

            // start new color class
            if (k > max_k) {
                max_k = k;
                colors[max_k+1].clear();
            }

            // add vertex to color class
            colors[k].push_back(u);

            if (k < min_k) {
                V[j].set_id(u);
                j++;
            }
        }
    }
    /*
     * Sort vertices such that vertices with largest coloring number are in the back of V
     * O(1) time to check coloring number
     */
    if (j > 0)  V[j-1].set_bound(0);
    if (min_k <= 0)  min_k = 1;

    for (k = min_k; k <= max_k; k++)
        for (int w = 0; w < colors[k].size(); w++) {
            V[j].set_id(colors[k][w]);
            V[j].set_bound(k);
            j++;
        }
    return V.back().get_bound();
}





void graphpack_graph::bound_stats(int alg, int lb, graphpack_graph& G) {
    double seconds = get_time();
    compute_triangle_cores();
    double tcore_sec = get_time() - seconds;

    double sec = get_time();
    compute_triangles();
    sec = get_time() - sec;

    double sec_core = get_time();
    compute_cores();
    sec_core = get_time() - sec_core;

    cout << "graph: " << fn <<endl;
    cout << "alg: " << alg <<endl;
    cout << "-------------------------------" <<endl;
    cout << "Graph Stats for Max-Clique:" <<endl;
    cout << "-------------------------------" <<endl;
    cout << "|V|: " << num_vertices() <<endl;
    cout << "|E|: " << num_edges() <<endl;
    cout << "d_max: " << get_max_degree() <<endl;
    cout << "d_avg: " << get_avg_degree() <<endl;
    cout << "p: " << density() <<endl;
    cout << "|T|: " << num_triangles() <<endl;
    cout << "T_avg: " << get_avg_triangles() <<endl;
    cout << "T_max: " << get_max_triangles() <<endl;
    cout << "cc_avg: " << get_cc_avg() <<endl;
    cout << "cc_global: " << get_cc_global() <<endl;

    cout << "--------------------------" <<endl;
    cout << "Maximum Clique Bounds:" <<endl;
    cout << "--------------------------" <<endl;
    cout << "d_max: " << get_max_degree() <<endl;
    cout << "T_ub: " << get_triangle_bound() <<endl;
    cout << "sqrt(2Tmaxi) computed in " << sec << " seconds" <<endl;
    cout << "K: " << get_max_core()+1 <<endl;
    cout << "Kcores computed in " << sec_core << " seconds" <<endl;
    cout << "T: " << get_triangle_core_bound() <<endl;
    cout << "Tcores computed in " << tcore_sec << " seconds" <<endl;
}

template<typename T>
void graphpack_graph::compute_ordering(vector<T>& bound, vector<T>& order) {
    long long n, d, start, num, md;
    long long v, md_end;

    n = bound.size();
    order.resize(n,0);
    vector < long long > pos(n);

    md = 0;
    for(v=1; v<n; v++)
        if (bound[v] > md)  md = bound[v];

    md_end = md+1;
    vector < long long > bin(md_end,0);

    for (v=1; v < n; v++) bin[bound[v]]++;

    start = 1;
    for (d=0; d < md_end; d++) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v=1; v<n; v++) {
        pos[v] = bin[bound[v]];
        order[pos[v]] = v;
        bin[bound[v]]++;
    }

    for (d=md; d > 1; d--)  bin[d] = bin[d-1];
    bin[0] = 1;

    for (v=0; v<n-1; v++) {
        bound[v] = bound[v+1];
        order[v] = order[v+1]-1;
    }
}

void graphpack_graph::vertex_bucket_sort(vector<int>& bound, vector<int>& order) {
    // int max = 0;
    // if (max == 0), then compute max,
    // otherwise ignore!
    long long n, d, start, num, md;
    long long v, md_end;

    n = bound.size();
    order.resize(n,0);
    vector < long long > pos(n);

    md = 0;
    for(v = 0; v < n; ++v)
        if (bound[v] > md)  md = bound[v];

    printf("[graphpack: vertex bucket sort] \t num vertices = %lld, max bound = %lld \n", n, md);

    // since if bound starts at 1
    md_end = md+1;
    vector < long long > bin(md_end,0);

    for (v = 0; v < n; ++v) bin[bound[v]]++;

    start = 0;
    for (d = 0; d < md_end; ++d) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    for (v = 0; v < n; ++v) {
        pos[v] = bin[bound[v]];
        order[pos[v]] = v;
        bin[bound[v]]++;
    }

//    for (d=md; d >= 0; d--)  bin[d] = bin[d-1];
//    bin[0] = 1;

//    for (v=0; v<n-1; v++) {
//        bound[v] = bound[v+1];
//        order[v] = order[v+1]-1;
//    }
}




template<typename VertexType>
void graphpack_graph::vertex_bucket_sort(vector<VertexType>& P, bool smallest_to_largest) {// = true) {
    //, vector<VertexType>& bound, vector<VertexType>& order) {
    long long n, d, start, num, max_bound;
    long long v, max_end;


    n = P.size(); // num of vertices

//    vector<int> bound(n,0);
//    vector<int> order(n,0);
    vector<Vertex> R;
    R.reserve(n);

    vector < long long > pos(n);

    max_bound = 0;
    for(v=1; v<n; ++v) {
        if (P[v].get_bound() > max_bound) {
            max_bound = P[v].get_bound();
        }
    }

    max_end = max_bound + 1;
    vector < int > bin(max_end,0);

    // count the number of times a particular upper bound value appears
    for (v = 0; v < n; v++)
        bin[P[v].get_bound()]++;


    start = 0;
    for (d = 0; d < max_end; ++d) {
        num = bin[d];
        bin[d] = start;
        start = start + num;
    }

    if (smallest_to_largest) {
        for (v = 0; v < n; ++v) {
            pos[v] = bin[P[v].get_bound()];
            P[pos[v]].set_id(P[v].get_id());
            P[pos[v]].set_bound(P[v].get_bound());
    //        order[pos[v]] = v;
            bin[P[v].get_bound()]++;
        }
    }
    else { // largest to smallest
        for (v = n-1; v >= 0; --v) {
            pos[v] = bin[P[v].get_bound()];
            P[pos[v]].set_id(v);
            P[pos[v]].set_bound(P[v].get_bound());
    //        order[pos[v]] = v;
            bin[P[v].get_bound()]++;
        }
    }

//    for (v = 0; v < n; ++v) {
//        pos[v] = bin[P[v].get_bound()];
////        R.push_back()
//        P[pos[v]].set_id(v);
//        P[pos[v]].set_bound(P[v].get_bound());
////        order[pos[v]] = v;
//        bin[P[v].get_bound()]++;
//    }

    // bin can be used to get the buckets fast
    /**

    for (d = max_bound; d > 0; --d)
        bin[d] = bin[d-1];
    bin[0] = 0;

    **/

//    for (v=0; v<n-1; v++) {
//        bound[v] = bound[v+1];
//        order[v] = order[v+1]-1;
//    }
}




void graphpack_graph::degree_bucket_sort_parallel() {
    degree_bucket_sort(false);
}


/*
 * Each worker sorts the edges of a vertex, and
 * afterwards, grabs the next available vertex to sort
 *
 * Note: sort neighbors by degree (largest to smallest)
 */
void graphpack_graph::degree_bucket_sort_parallel(bool desc) {

//  int v, u, n, md, md_end, start, d, num;

    vector<int> tmp_edges(edges.size(),0);
//  tmp_edges.reserve(edges.size());
    int v = 0;
    double sec = get_time();
    #pragma omp parallel for schedule(dynamic) \
        shared(tmp_edges)
//  firstprivate(colors,ind,vs,es) \
//  private(u, P, C, C_max, mc_cur)
    for (v = 0; v < num_vertices(); ++v) {

        int u, n, md, md_end, start, d, num;

        n = vertices[v+1] - vertices[v] + 1;
        vector<int> vert(n);
        vector<int> pos(n);
        vector<int> deg(n);

        md = 0;
        for(u=1; u<n; u++) {
            deg[u] = degree[edges[vertices[v] + (u-1)]];
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
            vert[pos[u]] = edges[vertices[v] + (u-1)];
            bin[deg[u]]++;
        }

//      int insert_pos = vertices[v];
        if (desc) {

            // largest to smallest, example: 5,3,2
            for (u = 1; u < n; ++u) {
//                tmp_edges[vertices[v] + (u-1)] = vert[u];
                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
            }

//            for (u=n-1; u>0; --u) {
////                tmp_edges[vertices[v] + (u-1)] = vert[u];
//                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
//            }

//            print_break();
//            cout << "v = " << v <<endl;
//            for (int ii=0; ii < num_vertices(); ii++)
//                for (u=vertices[ii]; u < vertices[ii+1]; u++) {
//                    cout << "\t(" << ii << "," << tmp_edges[u];
//                    cout << "), d = " << degree[tmp_edges[u]] <<endl;
//                }

            // largest to smallest
//          tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.rbegin(),vert.rend()-1);
        }
        else {
            //from smallest degree to largest
            for (u=n-1; u>0; --u) {
                tmp_edges[vertices[v] + (u-1)] = vert[u];
//                tmp_edges[vertices[v] + (u-1)] = vert[n-u];
            }

//            for (u=1; u<n; u++) {
//                tmp_edges[vertices[v] + (u-1)] = vert[u];
//            }

            //from smallest degree to largest
//          tmp_edges.insert(tmp_edges.begin()+insert_pos, vert.begin()+1,vert.end());
        }
    }

    if (verbose) {
        cout << "[graphpack: sorting neighbors]  |E| = " << edges.size();
        cout << ", |E_sorted| = " << tmp_edges.size() <<endl;
        cout << "sorting neighbors took " << get_time() - sec <<endl;
    }
    edges = tmp_edges;
}




void graphpack_graph::degree_bucket_sort() {
    degree_bucket_sort(false);
}


// sort neighbors by degree (largest to smallest)
void graphpack_graph::degree_bucket_sort(bool desc) {

    int v, u, n, md, md_end, start, d, num;

    double sec = get_time();
    vector<int> tmp_edges;
    tmp_edges.reserve(edges.size());

    for (v = 0; v < num_vertices(); v++) {

        n = vertices[v+1] - vertices[v] + 1;
        vector<int> vert(n);
        vector<int> pos(n);
        vector<int> deg(n);

        md = 0;
        for(u=1; u<n; u++) {
            deg[u] = degree[edges[vertices[v] + (u-1)]];
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
            vert[pos[u]] = edges[vertices[v] + (u-1)];
            bin[deg[u]]++;
        }

        if (desc) {
            // largest to smallest
            tmp_edges.insert(tmp_edges.end(),vert.rbegin(),vert.rend()-1);
        }
        else {
            //from smallest degree to largest
            tmp_edges.insert(tmp_edges.end(),vert.begin()+1,vert.end());
        }
    }

    cout << "[graphpack: sorting neighbors]  |E| = " << edges.size();
    cout << ", |E_sorted| = " << tmp_edges.size() <<endl;
    cout << "sorting neighbors took " << get_time() - sec <<endl;
    edges = tmp_edges;
}



// note when reducing graph explicitly, then forced to read in or keep a copy
bool graphpack_graph::clique_test(graphpack_graph& G, vector<int> C) {
    int u = 0;
    vector<short> ind(G.num_vertices(),0);
    for (size_t i = 0; i < C.size(); i++) ind[C[i]] = 1;


    // ensure each vertex in C has |C|-1 edges between each other
    for (size_t i = 0; i < C.size(); i++) {
        u = C[i];
        int sz = 0;
        for (long long j = G.vertices[u]; j < G.vertices[u+1]; j++)
            if (ind[G.edges[j]])  sz++;

        // check if connected to |C|-1 vertices
        if (sz != C.size()-1)
            return false;
    }
    return true;
}

void graphpack_graph::print_base_stats() {
    cout << "|V| = " << num_vertices() << ", vertices.size() = " << vertices.size() <<endl;
    cout << "|E| = " << num_edges() << ", edges.size() = " << edges.size() <<endl;
}

void graphpack_graph::print_edges() {
    cout << "[graphpack debug]  printing edges..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
//            int u = edges[j];
            cout << "(" << v << "," << edges[j] << ")" <<endl;
        }
    }
}



void graphpack_graph::print_vertices_array() {
    cout << "[graphpack debug]  printing vertex array..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        cout << "v = " << v << ", ptr = " << vertices[v] <<endl;
    }
}

void graphpack_graph::print_weighted_graph() {
    cout << "[graphpack debug]  printing vertex array..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        cout << "v = " << v << ", ptr = " << vertices[v] <<endl;
    }

    print_line(80);

    for (int v = 0; v < num_vertices(); v++) {
        printf("v = %d, N(v) =  ", v);
        for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
            int u = edges[j];
            double weight = wt[j];
            printf("%d [wt=%lg],  ", u, weight);
        }
        printf("\n");
    }
    print_line(80);

    printf("mapping vertex identifiers back to the original given as params \n");
    print_line(80);
    for (int v = 0; v < num_vertices(); v++) {
        printf("v = %d, N(v) =  ", v);
        for (long long j = vertices[v]; j < vertices[v + 1]; j++) {
            int u = get_params_vertex_id(edges[j]);
            double weight = wt[j];
            printf("%d [wt=%lg],  ", u, weight);
        }
        printf("\n");
    }
    print_line(80);

}


void graphpack_graph::print_degrees() {
    cout << "[graphpack debug]  printing vertex degrees..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        cout << "v = " << v << ", d_v = " << degree[v] <<endl;
    }
}


void graphpack_graph::write_edge_list() {
    write_edge_list("edge_list/");
}

void graphpack_graph::write_edge_list(string dir) {

    string::size_type result;
    string::size_type beg;
    //string name = "edge_list/";
    string name = dir;
    cout << "fn: " << fn <<endl;

    string rawname, fn_str;
    int last_index = fn.find_last_of(".");
    if (last_index != string::npos) {
        rawname = fn.substr(0, last_index);
        cout << rawname <<endl;
        int beg_index = rawname.find_last_of("/");
        string fn_new = rawname.substr(beg_index,rawname.size()-1);
        cout << "rawname: " << fn_new <<endl;
//        fn_str = string("edge_list") + fn_new + ".txt"; //string("-tcore.txt");
        fn_str = string(dir) + fn_new + ".txt"; //string("-tcore.txt");
        cout << "fn_str: " << fn_str <<endl;
    }

    ofstream myfile;
    char *fileName = (char*)fn_str.c_str();
    myfile.open (fileName);


    cout << "[mcpack]  writing tab delimited edge list..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
//            int u = edges[j];
            myfile << v << "\t" << edges[j] << "\t" << tri_core[eid[j]+1] <<endl;
        }
    }
    myfile.close();
}


void graphpack_graph::write_only_edge_list(string dir) {

    string::size_type result;
    string::size_type beg;
    //string name = "edge_list/";
    string name = dir;
    cout << "fn: " << fn <<endl;

    string rawname, fn_str;
    int last_index = fn.find_last_of(".");
    if (last_index != string::npos) {
        rawname = fn.substr(0, last_index);
        cout << rawname <<endl;
        int beg_index = rawname.find_last_of("/");
        string fn_new = rawname.substr(beg_index,rawname.size()-1);
        cout << "rawname: " << fn_new <<endl;
//        fn_str = string("edge_list") + fn_new + ".txt"; //string("-tcore.txt");
        fn_str = string(dir) + fn_new + ".txt"; //string("-tcore.txt");
        cout << "fn_str: " << fn_str <<endl;
    }

    ofstream myfile;
    char *fileName = (char*)fn_str.c_str();
    myfile.open (fileName);


    cout << "[mcpack]  writing tab delimited edge list..." <<endl;
    print_break();
    for (int v = 0; v < num_vertices(); v++) {
        for (long long j = vertices[v]; j < vertices[v + 1]; j++ ) {
            myfile << v << "\t" << edges[j] << "\t" << 1 <<endl;
        }
    }
    myfile.close();
}



template<typename T>
void graphpack_graph::write_file(std::vector<T> &data, string &filename) {

    ofstream myfile;
    char *fn = (char*)filename.c_str();
    myfile.open(fn);

    for (long long e = 0; e < data.size(); e++) {
        myfile << data[e] << "\n";
    }
    myfile.close();
}


template<typename T>
void graphpack_graph::write_vector(std::vector<T> &data, string suffix) {
    string::size_type result;
    string::size_type beg;
    string name = "stats/";
    cout << "fn: " << fn <<endl;

    string rawname, fn_str;
    int last_index = fn.find_last_of(".");
    if (last_index != string::npos) {
        rawname = fn.substr(0, last_index);
        cout << rawname <<endl;
        int beg_index = rawname.find_last_of("/");
        string fn_new = rawname.substr(beg_index,rawname.size()-1);
        cout << "rawname: " << fn_new <<endl;
        fn_str = string("stats") + fn_new + suffix; //string("-tcore.txt");
        cout << "fn_str: " << fn_str <<endl;
    }

    ofstream myfile;
    char *fileName = (char*)fn_str.c_str();
    myfile.open (fileName);


//    int sum = 0;
    for (long long e = 0; e < data.size(); e++) {
        myfile << data[e] << "\n";// <<endl;
//        sum += data[e];
    }
    myfile.close();
}
template void graphpack_graph::write_vector(std::vector<long long> &data, string suffix);
template void graphpack_graph::write_vector(std::vector<int> &data, string suffix);
template void graphpack_graph::write_vector(std::vector<double> &data, string suffix);


/**
 * @brief write a string to a file and return
 * @param filepath
 * @param data
 */
void graphpack_graph::write_line(const char *filepath, const char *data)
{
    FILE *fp = fopen(filepath, "ab");
    if (fp != NULL)
    {
        fputs(data, fp);
        fclose(fp);
    }
}

/**
 * @brief Gets the filename from an arbitrary path
 * @param s
 * @return
 */
string graphpack_graph::get_filename_from_path(const string& s) {
   char sep = '/';
#ifdef _WIN32
   sep = '\\';
#endif
   size_t i = s.rfind(sep, s.length( ));
   if (i != string::npos) {
      return(s.substr(i+1, s.length( ) - i));
   }
   return("");
}


void graphpack_graph::greedy_coloring_perf() {

    int color_num_kcore_incr = coloring_number("kcore",false);
    int color_num_kcore_decr = coloring_number("kcore",true);
    int color_num_deg_incr = coloring_number("deg",false);
    int color_num_deg_decr = coloring_number("deg",true);
    int color_num_dual_deg_incr = coloring_number("dual_deg",false);
    int color_num_dual_deg_decr = coloring_number("dual_deg",true);
    int color_num_dual_kcore_incr = coloring_number("dual_kcore",false);
    int color_num_dual_kcore_decr = coloring_number("dual_kcore",true);
    int color_num_rand = coloring_number("rand",false);

    cout << "[graphpack: coloring number]  X(kcore, incr) = " << color_num_kcore_incr;
    cout << ",\t\t\t X(kcore, decr) = " << color_num_kcore_decr <<endl;
    cout << "[graphpack: coloring number]  X(deg, incr) = " << color_num_deg_incr;
    cout << ",\t\t\t X(deg, decr) = " << color_num_deg_decr <<endl;
    cout << "[graphpack: coloring number]  X(dual_deg, incr) = " << color_num_dual_deg_incr;
    cout << ",\t\t X(dual_deg, decr): " << color_num_dual_deg_decr <<endl;
    cout << "[graphpack: coloring number]  X(dual_kcore, incr) = " << color_num_dual_kcore_incr;
    cout << ",\t\t X(dual_kcore, decr) = " << color_num_dual_kcore_decr <<endl;
    cout << "[graphpack: coloring number]  X(rand) = " << color_num_rand <<endl;
}

graphpack_graph graphpack_graph::preprocess_triangle_core(graphpack_graph & G, int & lb) {

    // this code relaxes the pruning using the lb
    vector<int> W;
    for (int i = G.num_vertices()-1; i >= 0; i--) {
        if (G.kcore[G.kcore_order[i]] >= lb)
            W.push_back(G.kcore_order[i]);
    }

    double sec = get_time();
    graphpack_graph H(W,G.vertices,G.edges);

    /*
     * Compute Triangle Cores from H; the induced k-core graph
     */
    int threshold = 15000;
    double fast_triangle_core_perf = 0;
    double fast_triangles_perf = 0;
    double fast_triangle_cores = 0;
    if (H.vertices.size() < threshold && H.adj != NULL) {
        cout << "verts in adj: " << H.vertices.size() <<endl;
        sec = get_time();
        H.triangle_cores_approx_adj(fast_triangle_core_perf,fast_triangles_perf);
        fast_triangle_cores = get_time() - sec;
        cout << "[mcpack]  time to compute triangle cores from H: " << get_time() - sec << endl;
    }
    else {
        sec = get_time();
        H.triangle_cores_approx(fast_triangle_core_perf,fast_triangles_perf);
        fast_triangle_cores = get_time() - sec;
    }


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
    if (lb == G.max_tri_core+2) {
        cout << "inside tricore pruning, lb == max tri core..." <<endl;
        for (int i = 0; i < H.num_vertices(); i++) {
            start = E.size();
            for (long long j = H.vertices[i]; j < H.vertices[i + 1]; j++ ) {

                edge_pos = H.eid[j]; // tri_core is m/2, so this finds the right tri_core value given an edge id
                if (H.tri_core[edge_pos]+2 >= lb) { // get correct tri_core for edge's eid
                    E.push_back(H.edges[j]);
                }
                else {
                    num_edges_pruned++;
                }
            }
            if ((E.size() - start) >= lb-1) { // basically the updated degree bound
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
                if (H.tri_core[edge_pos]+2 >= lb) { // get correct tri_core for edge's eid
                    E.push_back(H.edges[j]);
                }
                else {
                    num_edges_pruned++;
                }
            }
            if ((E.size() - start) >= lb-1) { // basically the updated degree bound, but what about the edges???
                V.push_back(Vertex(i,0));
            }
            else num_verts_pruned++;
            VS[i] = start;
            VS[i + 1] = E.size();
        }
    }
    cout << "after triangle core pruning..." <<endl;

    cout << "num_verts pruned = " << num_verts_pruned << " out of " << G.num_vertices() <<endl;
    cout << "num edges pruned = " << num_edges_pruned << " out of " << G.edges.size() <<endl;

    cout << "\n**********************************" <<endl;
    cout << "T: " << G.max_tri_core+2 <<endl;
    cout << "T(H): " << H.max_tri_core+2 <<endl;
    cout << "T_fast_time:\t " << fast_triangle_cores <<endl;


    G.basic_stats("G-");
    cout << "K-|T|: " << H.total_t <<endl;
    cout << "K-|T|_time: " << fast_triangles_perf <<endl;

    graphpack_graph F(V,VS,E);
    F.basic_stats("T-");
    return F;
}



/**
 * Positive assortativity coefficient indicates that vertices tend
 * to have an edge to other vertices with the same or similar degree.
 */
void graphpack_graph::compute_assort() {
	long long n, m, v, e;
	n = vertices.size() - 1;   m = e_v.size();

//	vector <int> deg(n,0);
//	for(v=0; v<n; v++)
//	    deg[v] = vertices[v+1] - vertices[v];

	double mu = 0, jd = 0, psi = 0, tau = 0; r = 0;
	for (e=0; e<e_v.size(); e++) {
		mu += double(degree[e_v[e]] * degree[e_u[e]]);
		jd += (0.5 * double((degree[e_v[e]] + degree[e_u[e]])));
		psi += (0.5 * double((degree[e_v[e]]*degree[e_v[e]]) + (degree[e_u[e]]*degree[e_u[e]])));
	}
	tau = (jd/m)*(jd/m);
	r = (mu/m - tau) / (psi/m - tau);
//	deg.clear();
}


void graphpack_graph::compute_assort_parallel(int block_size) {
    long long n, m, v, e;
    n = vertices.size() - 1;
    m = e_v.size();

//  vector <int> deg(n,0);
//  for(v=0; v<n; v++)
//      deg[v] = vertices[v+1] - vertices[v];

    double mu = 0, jd = 0, psi = 0, tau = 0; r = 0;

#pragma omp parallel for schedule(dynamic,block_size) \
    shared(n,m) reduction(+:mu,jd,psi)
    for (e=0; e<e_v.size(); e++) {
        mu += double(degree[e_v[e]] * degree[e_u[e]]);
        jd += (0.5 * double((degree[e_v[e]] + degree[e_u[e]])));
        psi += (0.5 * double((degree[e_v[e]]*degree[e_v[e]]) + (degree[e_u[e]]*degree[e_u[e]])));
    }
    tau = (jd/m)*(jd/m);
    r = (mu/m - tau) / (psi/m - tau);
//  deg.clear();
}




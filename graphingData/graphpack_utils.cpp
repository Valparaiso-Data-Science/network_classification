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

#include "graphpack_utils.h"

using namespace std;

bool fexists(const char *filename) {
    ifstream ifile(filename);
    return (bool)ifile;
}

#if defined(_WIN32) || defined(_WIN64)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
double get_time() {
    LARGE_INTEGER t, freq;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&freq);
    return double(t.QuadPart)/double(freq.QuadPart);
}
#else
#include <sys/types.h>
#include <sys/timeb.h>
#include <sys/time.h>
double get_time() {
    struct timeval t;
    gettimeofday(&t, 0);
    return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
}
#endif

double tic() { return get_time(); }
void toc(double & start) {
    start = get_time() - start;
}

// works only for linux
string memory_usage() {
    ostringstream mem;
    ifstream proc("/proc/self/status");
    string s;
    while(getline(proc, s), !proc.fail()) {
        if(s.substr(0, 6) == "VmSize") {
            mem << s;
            return mem.str();
        }
    }
    return mem.str();
}


std::ifstream::pos_type get_filesize(const char* filename) {
    std::ifstream in(filename, std::ifstream::ate | std::ifstream::binary);
    return in.tellg();
}

void indent(int level = 0, string str = "") {
    for (int i = 0; i < level; i++)
        cout << "   ";
    cout << "(" << level << ") ";
}

void print_max_clique(vector<int>& C) {
    cout << "Maximum clique: ";
    for(int i = 0; i < C.size(); i++)
        cout << C[i] + 1 << " ";
    cout << endl;
}

void print_n_maxcliques(set< vector<int> > C, int n) {
    set< vector<int> >::iterator it;
    int mc = 0;
    for( it = C.begin(); it != C.end(); it++) {
        if (mc < n) {
            cout << "Maximum clique: ";
            const vector<int>& clq = (*it);
            for (int j = 0; j < clq.size(); j++)
                cout << clq[j] << " ";
            printf("\n\n");
            ++mc;
        }
        else break;
    }
}

void validate(bool condition, const string& msg) {
    if (!condition) {
        cerr << msg << endl;
        assert(condition);
    }
}

int getdir (string dir, vector<string> &files) {
    DIR *dp;
    struct dirent *dirp;
    if((dp  = opendir(dir.c_str())) == NULL) {
        cout << "Error(" << errno << ") opening " << dir << endl;
        return errno;
    }

    while ((dirp = readdir(dp)) != NULL) {
        if (dirp->d_name != ".")
            files.push_back(string(dirp->d_name));
    }
    closedir(dp);
    return 0;
}

void print_line(int n) {
    for (int i = 0; i < n; ++i) printf("-");
    printf("\n");
}

template<typename T>
void write_vector(std::vector<T> &data, string suffix, string fn) {
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


    int sum = 0;
    for (long long e = 0; e < data.size(); e++) {
        myfile << data[e] << "\n";// <<endl;
        sum += data[e];
    }
    myfile.close();
}
template void write_vector(std::vector<long long> &data, string suffix, string fn);
template void write_vector(std::vector<int> &data, string suffix, string fn);
template void write_vector(std::vector<double> &data, string suffix, string fn);

string extract_filename(string fn, bool remove_ext) { //, bool remove_ext = true) {
    if (remove_ext) {
        fn = remove_file_extension(fn);
    }
    int idx = 0;
    idx = fn.find_last_of("/")+1;
    return fn.substr(idx, idx - fn.size());//(fn.size() - fn.find_last_of(".")) );
}

string remove_file_extension(string fn) {
    return fn.substr(0,fn.size() - (fn.size()-fn.find_last_of(".")));
}

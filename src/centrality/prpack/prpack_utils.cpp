/**
 * @file prpack_utils.cpp
 * An assortment of utility functions for reporting errors, checking time,
 * and working with vectors.
 */

#include <stdlib.h>
#include "prpack_utils.h"
#include <cassert>
#include <iostream>
#include <string>
using namespace prpack;
using namespace std;

#ifdef PRPACK_IGRAPH_SUPPORT
#include "igraph_error.h"
#endif

#if defined(_WIN32)
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#endif
double prpack_utils::get_time() {
    LARGE_INTEGER t, freq;
    QueryPerformanceCounter(&t);
    QueryPerformanceFrequency(&freq);
    return double(t.QuadPart)/double(freq.QuadPart);
}
#else
#include <sys/types.h>
#include <sys/time.h>
double prpack_utils::get_time() {
    struct timeval t;
    gettimeofday(&t, 0);
    return (t.tv_sec*1.0 + t.tv_usec/1000000.0);
}
#endif

// Fails and outputs 'msg' if 'condition' is false.
void prpack_utils::validate(const bool condition, const string& msg) {
    if (!condition) {
#ifdef PRPACK_IGRAPH_SUPPORT
        igraph_error("Internal error in PRPACK", IGRAPH_FILE_BASENAME, __LINE__, 
	             IGRAPH_EINTERNAL);
#else
        cerr << msg << endl;
        exit(-1);
#endif
    }
}

// Permute a vector.
double* prpack_utils::permute(const int length, const double* a, const int* coding) {
    double* ret = new double[length];
    for (int i = 0; i < length; ++i)
        ret[coding[i]] = a[i];
    return ret;
}


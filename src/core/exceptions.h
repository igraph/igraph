#ifndef IGRAPH_HANDLE_EXCEPTIONS_H
#define IGRAPH_HANDLE_EXCEPTIONS_H

#include <exception>
#include <new>

/* igraph functions which may be called from C code must not throw C++ exceptions.
 * This includes all public functions. This macro is meant to handle exceptions thrown
 * by C++ libraries used by igraph (such as bliss). Wrap the entire body
 * of public functions implemented in C++ in IGRAPH_HANDLE_EXCEPTIONS().
 */
#define IGRAPH_HANDLE_EXCEPTIONS(code) \
    try { code; } \
    catch (const std::bad_alloc &e) { IGRAPH_ERROR(e.what(), IGRAPH_ENOMEM); } \
    catch (const std::exception &e) { IGRAPH_ERROR(e.what(), IGRAPH_FAILURE); } \
    catch (...) { IGRAPH_ERROR("Unknown exception caught.", IGRAPH_FAILURE); }


#endif // IGRAPH_HANDLE_EXCEPTIONS_H

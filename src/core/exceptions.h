#ifndef IGRAPH_HANDLE_EXCEPTIONS_H
#define IGRAPH_HANDLE_EXCEPTIONS_H

#include <exception>
#include <new>

#define IGRAPH_HANDLE_EXCEPTIONS(code) \
    try { code; } \
    catch (const std::bad_alloc &e) { IGRAPH_ERROR(e.what(), IGRAPH_ENOMEM); } \
    catch (const std::exception &e) { IGRAPH_ERROR(e.what(), IGRAPH_FAILURE); } \
    catch (...) { IGRAPH_ERROR("Unknown exception caught.", IGRAPH_FAILURE); }


#endif // IGRAPH_HANDLE_EXCEPTIONS_H

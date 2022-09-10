#ifndef IGRAPH_HANDLE_EXCEPTIONS_H
#define IGRAPH_HANDLE_EXCEPTIONS_H

#include "igraph_error.h"

#include <exception>
#include <new>
#include <stdexcept>

/* igraph functions which may be called from C code must not throw C++ exceptions.
 * This includes all public functions. This macro is meant to handle exceptions thrown
 * by C++ libraries used by igraph (such as bliss). Wrap the entire body
 * of public functions implemented in C++ in IGRAPH_HANDLE_EXCEPTIONS().
 *
 * In some cases IGRAPH_HANDLE_EXCEPTIONS() won't work because the
 * C preprocessor gets confused by the code block. In that case one can use
 * IGRAPH_HANDLE_EXCEPTIONS_BEGIN; and IGRAPH_HANDLE_EXCEPTIONS_END at the
 * beginning and end of the code block.
 */
#define IGRAPH_HANDLE_EXCEPTIONS_BEGIN \
    try {
#define IGRAPH_HANDLE_EXCEPTIONS_END \
    } \
    catch (const std::bad_alloc &e) { IGRAPH_ERROR(e.what(), IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */ } \
    catch (const std::range_error &e) { IGRAPH_ERROR(e.what(), IGRAPH_EOVERFLOW); /* LCOV_EXCL_LINE */ } \
    catch (const std::exception &e) { IGRAPH_ERROR(e.what(), IGRAPH_FAILURE); /* LCOV_EXCL_LINE */ } \
    catch (...) { IGRAPH_ERROR("Unknown exception caught.", IGRAPH_FAILURE); /* LCOV_EXCL_LINE */ }

#define IGRAPH_HANDLE_EXCEPTIONS(code) \
    IGRAPH_HANDLE_EXCEPTIONS_BEGIN; \
    code; \
    IGRAPH_HANDLE_EXCEPTIONS_END;

#endif // IGRAPH_HANDLE_EXCEPTIONS_H

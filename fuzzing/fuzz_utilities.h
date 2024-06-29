#ifndef FUZZ_UTILITIES_H
#define FUZZ_UTILITIES_H

#include <igraph.h>

#define CHECK_ERROR(funcall, expected_err) \
    do { \
            igraph_error_handler_t *handler; \
            handler = igraph_set_error_handler(igraph_error_handler_ignore); \
            IGRAPH_ASSERT(funcall == expected_err); \
            igraph_set_error_handler(handler); \
    } while (0)

#endif // FUZZ_UTILITIES_H

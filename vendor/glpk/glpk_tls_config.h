
#include "igraph_threading.h" /* IGRAPH_THREAD_SAFE */

/* This includes igraph's config.h.
 * The vendored GLPK must not have a config.h. */
#include "config.h"

#if IGRAPH_THREAD_SAFE
#define TLS IGRAPH_THREAD_LOCAL
#endif

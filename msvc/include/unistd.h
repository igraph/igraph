
/*
 * unistd.h replacement for MSVC
 *
 * Provide the minimum that igraph needs.
 * At presents this is:
 *  - isatty() for f2c and flex-generated sources
 *  - unlink() for examples/simple/graphml.c
 */

#include <io.h>

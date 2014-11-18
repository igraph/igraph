#include "uuid.h"

#include <Rinternals.h>

SEXP UUID_gen(SEXP sTime) {
    uuid_t u;
    char c[40];
    int use_time = asInteger(sTime);
    if (use_time == TRUE)
	uuid_generate_time(u);
    else if (use_time == FALSE)
	uuid_generate_random(u);
    else
	uuid_generate(u);
    uuid_unparse_lower(u, c);
    return mkString(c);
}


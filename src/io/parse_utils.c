
#include "parse_utils.h"

#include "igraph_memory.h"

#include <errno.h>
#include <limits.h>
#include <stdlib.h>
#include <string.h>

/* TODO: Support for reporting line number where parse error occurred. */

/* Converts a string to an integer. Throws an error if the result is not representable.
 *
 * The input is a not-necesarily-null-terminated string that must contain only the number.
 * Any additional characters at the end of the string, such as whitespace, will trigger
 * a parsing error.
 */
igraph_error_t igraph_i_parse_integer(const char *str, size_t length, igraph_integer_t *value) {
    char buffer[128];
    char *tmp, *end;
    char last_char;
    igraph_bool_t out_of_range, dynamic_alloc;
    long long val;

    dynamic_alloc = length+1 > sizeof(buffer) / sizeof(buffer[0]);

    if (dynamic_alloc) {
        tmp = IGRAPH_CALLOC(length+1, char);
        if (tmp == NULL) {
            IGRAPH_ERROR("Failed to parse integer.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
    } else {
        tmp = buffer;
    }

    strncpy(tmp, str, length);
    tmp[length]='\0';

    /* To avoid having to choose the appropriate strto?() function based on
     * the definition of igraph_integer_t, we first use a long long variable
     * which should be at least as large as igraph_integer_t on any platform. */
    errno = 0;
    val = strtoll(tmp, &end, 10);
    out_of_range = errno == ERANGE;
    *value = (igraph_integer_t) val;
    last_char = *end;
    if (*value != val) {
        out_of_range = 1;
    }

    /* Free memory before raising any errors. */
    if (dynamic_alloc) {
        IGRAPH_FREE(tmp);
    }

    if (out_of_range) {
        IGRAPH_ERROR("Failed to parse integer.", val > 0 ? IGRAPH_EOVERFLOW : IGRAPH_EUNDERFLOW);
    }

    /* Did we parse to the end of the string? */
    if (last_char) {
        IGRAPH_ERRORF("Unexpected character '%c' while parsing integer.", IGRAPH_PARSEERROR, last_char);
    }

    return IGRAPH_SUCCESS;
}


/* Converts a string to a real number. Throws an error if the result is not representable.
 *
 * The input is a not-necesarily-null-terminated string that must contain only the number.
 * Any additional characters at the end of the string, such as whitespace, will trigger
 * a parsing error.
 */
igraph_error_t igraph_i_parse_real(const char *str, size_t length, igraph_real_t *value) {
    char buffer[128];
    char *tmp, *end;
    char last_char;
    igraph_bool_t out_of_range, dynamic_alloc;

    dynamic_alloc = length+1 > sizeof(buffer) / sizeof(buffer[0]);

    if (dynamic_alloc) {
        tmp = IGRAPH_CALLOC(length+1, char);
        if (tmp == NULL) {
            IGRAPH_ERROR("Failed to parse real number.", IGRAPH_ENOMEM); /* LCOV_EXCL_LINE */
        }
    } else {
        tmp = buffer;
    }

    strncpy(tmp, str, length);
    tmp[length]='\0';

    errno = 0;
    *value = strtod(tmp, &end);
    out_of_range = errno == ERANGE; /* This does not trigger when reading +-Inf. */
    last_char = *end;

    /* Free memory before raising any errors. */
    if (dynamic_alloc) {
        IGRAPH_FREE(tmp);
    }

    if (out_of_range) {
        IGRAPH_ERROR("Failed to parse real number.", *value > 0 ? IGRAPH_EOVERFLOW : IGRAPH_EUNDERFLOW);
    }

    /* Did we parse to the end of the string? */
    if (last_char) {
        IGRAPH_ERRORF("Unexpected character '%c' while parsing real number.", IGRAPH_PARSEERROR, last_char);
    }

    return IGRAPH_SUCCESS;
}

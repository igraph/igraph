/* error.c
 *
 * Copyright (C) 2010-2011 Tamas Nepusz
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */

#include <stdio.h>
#include <stdlib.h>
#include "error.h"
#include "platform.h"

static char *plfit_i_error_strings[] = {
    "No error",
    "Failed",
    "Invalid value",
    "Underflow",
    "Overflow",
    "Not enough memory"
};

#ifndef USING_R
static plfit_error_handler_t* plfit_error_handler = plfit_error_handler_abort;
#else
/* This is overwritten, anyway */
static plfit_error_handler_t* plfit_error_handler = plfit_error_handler_ignore;
#endif

const char* plfit_strerror(const int plfit_errno) {
  return plfit_i_error_strings[plfit_errno];
}

plfit_error_handler_t* plfit_set_error_handler(plfit_error_handler_t* new_handler) {
    plfit_error_handler_t* old_handler = plfit_error_handler;
    plfit_error_handler = new_handler;
    return old_handler;
}

void plfit_error(const char *reason, const char *file, int line,
        int plfit_errno) {
    plfit_error_handler(reason, file, line, plfit_errno);
}

#ifndef USING_R
void plfit_error_handler_abort(const char *reason, const char *file, int line,
        int plfit_errno) {
    fprintf(stderr, "Error at %s:%i : %s, %s\n", file, line, reason,
            plfit_strerror(plfit_errno));
    abort();
}
#endif

#ifndef USING_R
void plfit_error_handler_printignore(const char *reason, const char *file, int line,
        int plfit_errno) {
    fprintf(stderr, "Error at %s:%i : %s, %s\n", file, line, reason,
            plfit_strerror(plfit_errno));
}
#endif

void plfit_error_handler_ignore(const char* reason, const char* file, int line,
		int plfit_errno) {
}

/*
   IGraph library.
   Copyright (C) 2021  The igraph development team <igraph@igraph.org>

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/

#include <igraph_error.h>
#include "test_utilities.inc"

int cause_error() {
    IGRAPH_ERRORF("%d %f %ld %c", IGRAPH_EINVAL, 1, 1.0, 1, 'a');
    return IGRAPH_SUCCESS;
}

int cause_warning() {
    IGRAPH_WARNINGF("%d %f %ld %c", 1, 1.0, 1, 'a');
    return IGRAPH_SUCCESS;
}

int cause_fatal() {
    IGRAPH_FATALF("%d %f %ld %c", 1, 1.0, 1, 'a');
}

void error_handler(const char *reason, const char *file, int line, int igraph_errno) {
    printf("Error. Reason: %s\nErrno: %d\n", reason, igraph_errno);
}

void warning_handler(const char *reason, const char *file, int line, int igraph_errno) {
    printf("Warning. Reason: %s\nErrno: %d\n", reason, igraph_errno);
}

void fatal_handler(const char *reason, const char *file, int line) {
    printf("Fatal. Reason: %s\n", reason);
    exit(0);
}

int main() {
    igraph_set_error_handler(&error_handler);
    igraph_set_warning_handler(&warning_handler);
    igraph_set_fatal_handler(&fatal_handler);

    IGRAPH_ASSERT(cause_error() == IGRAPH_EINVAL);
    IGRAPH_ASSERT(cause_warning() == IGRAPH_SUCCESS);
    cause_fatal();

    /* The igraph_fatal() call must not return, so the following lines should not run. */

    printf("This should not be printed.");

    return 1;
}

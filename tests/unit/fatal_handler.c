/*
   IGraph library.
   Copyright (C) 2021-2022  The igraph development team <igraph@igraph.org>

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

#include <igraph.h>
#include <stdio.h>

igraph_fatal_handler_t hanlder;

void handler(const char *reason, const char *file, int line) {
    printf("Reason: %s\nFile: %s\nLine: %d\n", reason, file, line);
    exit(0); /* We use exit(0) instead of abort() to allow the test to succeed. */
}

int main(void) {
    igraph_set_fatal_handler(&handler);

    igraph_fatal("REASON", "FILENAME", 123);

    /* The igraph_fatal() call must not return, so the following lines should not run. */

    printf("This should not be printed.");

    return 0;
}

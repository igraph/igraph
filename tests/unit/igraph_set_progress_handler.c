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

#include <igraph.h>
#include "test_utilities.inc"

int handler(const char* message, igraph_real_t percent, void*data) {
    printf("handler, %s, %f, %d\n", message, percent, *(int*)data);
    return IGRAPH_SUCCESS;
}

int main() {
    igraph_set_progress_handler(handler);
    int data = 10;

    printf("progress with set progress handler:\n");
    IGRAPH_PROGRESS("message", 100.0, &data);

    igraph_progress_handler_t *previous = igraph_set_progress_handler(NULL);

    printf("\nprogress with no handler:\n");
    IGRAPH_PROGRESS("message", 100.0, &data);

    igraph_set_progress_handler(previous);

    printf("\nprogress with previous handler:\n");
    IGRAPH_PROGRESS("message", 100.0, &data);

    VERIFY_FINALLY_STACK();
    return 0;
}

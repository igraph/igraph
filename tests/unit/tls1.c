/* -*- mode: C -*-  */
/*
   IGraph library.
   Copyright (C) 2011-2012  Gabor Csardi <csardi.gabor@gmail.com>
   334 Harvard street, Cambridge MA, 02139 USA

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
   02110-1301 USA

*/

#include <igraph.h>
#include <pthread.h>

void *thread_function(void *arg) {
    IGRAPH_FINALLY(igraph_free, NULL);
    return 0;
}

int main() {
    pthread_t thread_id;
    void *exit_status;

    /* Skip if igraph is not thread-safe */
    if (!IGRAPH_THREAD_SAFE) {
        return 77;
    }

    /* Run a thread that leaves some junk in the error stack */
    pthread_create(&thread_id, NULL, thread_function, 0);
    pthread_join(thread_id, &exit_status);

    /* Check that the error stack is not common */
    if (!IGRAPH_FINALLY_STACK_EMPTY) {
        printf("Foobar\n");
        return 1;
    }

    return 0;
}

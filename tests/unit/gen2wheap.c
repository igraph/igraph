
#include <igraph.h>

#include "core/genheap.h"

#include "test_utilities.h"

typedef struct intpair_t {
    int x, y;
} intpair_t;

int cmp(const void *a, const void *b) {
    const intpair_t *sa = a, *sb = b;
    if (sa->x < sb->x) return -1;
    if (sa->x > sb->x) return 1;
    if (sa->y < sb->y) return -1;
    if (sa->y > sb->y) return 1;
    return 0;
}

void intpair_print(const intpair_t *p) {
    printf("(%d, %d)", p->x, p->y);
}

int main(void) {
    igraph_integer_t n;
    igraph_vector_int_t idx;
    igraph_gen2wheap_t h;
    igraph_integer_t i;
    intpair_t p;

    igraph_rng_seed(igraph_rng_default(), 42);

    n = 10;

    igraph_gen2wheap_init(&h, cmp, sizeof(intpair_t), n);
    IGRAPH_ASSERT(igraph_gen2wheap_max_size(&h) == n);

    igraph_vector_int_init_range(&idx, 0, n);
    igraph_vector_int_shuffle(&idx);

    for (i=0; i < n; i++) {
        igraph_integer_t j = VECTOR(idx)[i];
        p.x = RNG_INTEGER(1, 10);
        p.y = RNG_INTEGER(1, 10);
        printf("Adding %2" IGRAPH_PRId ": ", j);
        intpair_print(&p);
        printf("\n");
        IGRAPH_ASSERT(! igraph_gen2wheap_has_elem(&h, j));
        igraph_gen2wheap_push_with_index(&h, j, &p);
        igraph_gen2wheap_check(&h);
        IGRAPH_ASSERT(igraph_gen2wheap_has_elem(&h, j));
        IGRAPH_ASSERT(igraph_gen2wheap_size(&h) == i+1);
    }

    printf("\n");

    i = igraph_gen2wheap_size(&h);
    while (!igraph_gen2wheap_empty(&h)) {
        printf("Removing %2" IGRAPH_PRId ": ", igraph_gen2wheap_max_index(&h));
        intpair_print(igraph_gen2wheap_max(&h));
        printf("\n");
        igraph_gen2wheap_delete_max(&h);
        igraph_gen2wheap_check(&h);
        i--;
        IGRAPH_ASSERT(igraph_gen2wheap_size(&h) == i);
    }

    printf("\n");

    n=5;
    for (igraph_integer_t i=0; i < n; i++) {
        p.x = RNG_INTEGER(1, 10);
        p.y = RNG_INTEGER(1, 10);
        printf("Adding %2" IGRAPH_PRId ": ", i);
        intpair_print(&p);
        printf("\n");
        IGRAPH_ASSERT(! igraph_gen2wheap_has_elem(&h, i));
        igraph_gen2wheap_push_with_index(&h, i, &p);
        igraph_gen2wheap_check(&h);
        IGRAPH_ASSERT(igraph_gen2wheap_has_elem(&h, i));
        IGRAPH_ASSERT(igraph_gen2wheap_size(&h) == i+1);
    }

    i = 1;
    p.x = -1; p.y = -1;
    printf("\nModifying %" IGRAPH_PRId " to ", i);
    intpair_print(&p);
    printf("\n");
    igraph_gen2wheap_modify(&h, i, &p);
    igraph_gen2wheap_check(&h);

    i = igraph_gen2wheap_max_index(&h);
    p = * (intpair_t *) igraph_gen2wheap_max(&h);
    p.x -= 3;
    printf("Modifying max to ");
    intpair_print(&p);
    printf("\n");
    igraph_gen2wheap_modify(&h, i, &p);
    igraph_gen2wheap_check(&h);

    printf("\n");

    i = igraph_gen2wheap_size(&h);
    while (!igraph_gen2wheap_empty(&h)) {
        printf("Removing %2" IGRAPH_PRId ": ", igraph_gen2wheap_max_index(&h));
        intpair_print(igraph_gen2wheap_max(&h));
        printf("\n");
        igraph_gen2wheap_delete_max(&h);
        igraph_gen2wheap_check(&h);
        i--;
        IGRAPH_ASSERT(igraph_gen2wheap_size(&h) == i);
    }

    printf("\n");

    igraph_vector_int_destroy(&idx);
    igraph_gen2wheap_destroy(&h);

    VERIFY_FINALLY_STACK();

    return 0;
}

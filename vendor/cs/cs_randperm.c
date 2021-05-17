#include "cs.h"

#include "igraph_random.h"

/* return a random permutation vector, the identity perm, or p = n-1:-1:0.
 * seed = -1 means p = n-1:-1:0.  seed = 0 means p = identity.  otherwise
 * p = random permutation.  */
CS_INT *cs_randperm (CS_INT n, CS_INT seed)
{
    CS_INT *p, k, j, t ;
    if (seed == 0) return (NULL) ;      /* return p = NULL (identity) */
    p = cs_malloc (n, sizeof (CS_INT)) ;   /* allocate result */
    if (!p) return (NULL) ;             /* out of memory */
    for (k = 0 ; k < n ; k++) p [k] = n-k-1 ;
    if (seed == -1) return (p) ;        /* return reverse permutation */
    /* srand (seed) ;                      /\* get new random number seed *\/ */
    RNG_BEGIN();
    for (k = 0 ; k < n ; k++)
    {
        /* j = k + (rand ( ) % (n-k)) ;    /\* j = rand CS_INT in range k to n-1 *\/ */
        j = RNG_INTEGER(k, n-1) ;
        t = p [j] ;                     /* swap p[k] and p[j] */
        p [j] = p [k] ;
        p [k] = t ;
    }
    RNG_END();
    return (p) ;
}

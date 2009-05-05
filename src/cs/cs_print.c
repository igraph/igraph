#include "cs.h"
/* print a sparse matrix */
CS_INT cs_print (const cs *A, CS_INT brief)
{
    CS_INT p, j, m, n, nzmax, nz, *Ap, *Ai ;
    CS_ENTRY *Ax ;
    if (!A) { printf ("(null)\n") ; return (0) ; }
    m = A->m ; n = A->n ; Ap = A->p ; Ai = A->i ; Ax = A->x ;
    nzmax = A->nzmax ; nz = A->nz ;
    printf ("CXSparse Version %d.%d.%d, %s.  %s\n", CS_VER, CS_SUBVER,
        CS_SUBSUB, CS_DATE, CS_COPYRIGHT) ;
    if (nz < 0)
    {
        printf (""CS_ID"-by-"CS_ID", nzmax: "CS_ID" nnz: "CS_ID", 1-norm: %g\n", m, n, nzmax,
                Ap [n], cs_norm (A)) ;
        for (j = 0 ; j < n ; j++)
        {
            printf ("    col "CS_ID" : locations "CS_ID" to "CS_ID"\n", j, Ap [j], Ap [j+1]-1);
            for (p = Ap [j] ; p < Ap [j+1] ; p++)
            {
#ifdef CS_COMPLEX
                printf ("      "CS_ID" : (%g, %g)\n", Ai [p], 
		    Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
                printf ("      "CS_ID" : %g\n", Ai [p], Ax ? Ax [p] : 1) ;
#endif
                if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
            }
        }
    }
    else
    {
        printf ("triplet: "CS_ID"-by-"CS_ID", nzmax: "CS_ID" nnz: "CS_ID"\n", m, n, nzmax, nz) ;
        for (p = 0 ; p < nz ; p++)
        {
#ifdef CS_COMPLEX
            printf ("    "CS_ID" "CS_ID" : (%g, %g)\n", Ai [p], Ap [p], 
		    Ax ? CS_REAL (Ax [p]) : 1, Ax ? CS_IMAG (Ax [p]) : 0) ;
#else
            printf ("    "CS_ID" "CS_ID" : %g\n", Ai [p], Ap [p], Ax ? Ax [p] : 1) ;
#endif
            if (brief && p > 20) { printf ("  ...\n") ; return (1) ; }
        }
    }
    return (1) ;
}

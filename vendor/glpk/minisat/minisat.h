/* minisat.h */

/* Modified by Andrew Makhorin <mao@gnu.org>, August 2011 */

/***********************************************************************
*  MiniSat -- Copyright (c) 2005, Niklas Sorensson
*  http://www.cs.chalmers.se/Cs/Research/FormalMethods/MiniSat/
*
*  Permission is hereby granted, free of charge, to any person
*  obtaining a copy of this software and associated documentation files
*  (the "Software"), to deal in the Software without restriction,
*  including without limitation the rights to use, copy, modify, merge,
*  publish, distribute, sublicense, and/or sell copies of the Software,
*  and to permit persons to whom the Software is furnished to do so,
*  subject to the following conditions:
*
*  The above copyright notice and this permission notice shall be
*  included in all copies or substantial portions of the Software.
*
*  THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
*  EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
*  MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
*  NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
*  BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
*  ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
*  CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
*  SOFTWARE.
***********************************************************************/
/* Modified to compile with MS Visual Studio 6.0 by Alan Mishchenko */

#ifndef MINISAT_H
#define MINISAT_H

/*====================================================================*/
/* Simple types: */

typedef int bool;

#define true  1
#define false 0

typedef int  lit;
#if 0 /* by mao */
typedef char lbool;
#else
typedef int lbool;
#endif

#define var_Undef (int)(-1)
#define lit_Undef (lit)(-2)

#define l_Undef (lbool)0
#define l_True  (lbool)1
#define l_False (lbool)(-1)

#define toLit(v) (lit)((v) + (v))
#define lit_neg(l) (lit)((l) ^ 1)
#define lit_var(l) (int)((l) >> 1)
#define lit_sign(l) (int)((l) & 1)

/*====================================================================*/
/* Vectors: */

/* vector of 32-bit intergers (added for 64-bit portability) */
typedef struct /* veci_t */ {
    int    size;
    int    cap;
    int*   ptr;
} veci;

#define veci_new(v) \
{   (v)->size = 0; \
    (v)->cap  = 4; \
    (v)->ptr  = (int*)malloc(sizeof(int)*(v)->cap); \
}

#define veci_delete(v) free((v)->ptr)

#define veci_begin(v) ((v)->ptr)

#define veci_size(v) ((v)->size)

#define veci_resize(v, k) (void)((v)->size = (k))
/* only safe to shrink !! */

#define veci_push(v, e) \
{   if ((v)->size == (v)->cap) \
    {   int newsize = (v)->cap * 2+1; \
        (v)->ptr = (int*)realloc((v)->ptr,sizeof(int)*newsize); \
        (v)->cap = newsize; \
    } \
    (v)->ptr[(v)->size++] = (e); \
}

/* vector of 32- or 64-bit pointers */
typedef struct /* vecp_t */ {
    int    size;
    int    cap;
    void** ptr;
} vecp;

#define vecp_new(v) \
{   (v)->size = 0; \
    (v)->cap  = 4; \
    (v)->ptr  = (void**)malloc(sizeof(void*)*(v)->cap); \
}

#define vecp_delete(v) free((v)->ptr)

#define vecp_begin(v) ((v)->ptr)

#define vecp_size(v) ((v)->size)

#define vecp_resize(v, k) (void)((v)->size = (k))
/* only safe to shrink !! */

#define vecp_push(v, e) \
{   if ((v)->size == (v)->cap) \
    {   int newsize = (v)->cap * 2+1; \
        (v)->ptr = (void**)realloc((v)->ptr,sizeof(void*)*newsize); \
        (v)->cap = newsize; \
    } \
    (v)->ptr[(v)->size++] = (e); \
}

/*====================================================================*/
/* Solver representation: */

typedef struct /* clause_t */
{
    int size_learnt;
    lit lits[1];
} clause;

typedef struct /* stats_t */
{
    double   starts, decisions, propagations, inspects, conflicts;
    double   clauses, clauses_literals, learnts, learnts_literals,
             max_literals, tot_literals;
} stats;

typedef struct /* solver_t */
{
    int      size;          /* nof variables */
    int      cap;           /* size of varmaps */
    int      qhead;         /* Head index of queue. */
    int      qtail;         /* Tail index of queue. */

    /* clauses */
    vecp     clauses;       /* List of problem constraints.
                               (contains: clause*) */
    vecp     learnts;       /* List of learnt clauses.
                               (contains: clause*) */

    /* activities */
    double   var_inc;       /* Amount to bump next variable with. */
    double   var_decay;     /* INVERSE decay factor for variable
                               activity: stores 1/decay. */
    float    cla_inc;       /* Amount to bump next clause with. */
    float    cla_decay;     /* INVERSE decay factor for clause
                               activity: stores 1/decay. */

    vecp*    wlists;
    double*  activity;      /* A heuristic measurement of the activity
                               of a variable. */
    lbool*   assigns;       /* Current values of variables. */
    int*     orderpos;      /* Index in variable order. */
    clause** reasons;
    int*     levels;
    lit*     trail;

    clause*  binary;        /* A temporary binary clause */
    lbool*   tags;
    veci     tagged;        /* (contains: var) */
    veci     stack;         /* (contains: var) */

    veci     order;         /* Variable order. (heap) (contains: var) */
    veci     trail_lim;     /* Separator indices for different decision
                               levels in 'trail'. (contains: int) */
    veci     model;         /* If problem is solved, this vector
                               contains the model (contains: lbool). */

    int      root_level;    /* Level of first proper decision. */
    int      simpdb_assigns;/* Number of top-level assignments at last
                               'simplifyDB()'. */
    int      simpdb_props;  /* Number of propagations before next
                               'simplifyDB()'. */
    double   random_seed;
    double   progress_estimate;
    int      verbosity;     /* Verbosity level.
                               0=silent,
                               1=some progress report,
                               2=everything */

    stats    stats;
} solver;

/*====================================================================*/
/* Public interface: */

#if 1 /* by mao; to keep namespace clean */
#define solver_new        _glp_minisat_new
#define solver_delete     _glp_minisat_delete
#define solver_addclause  _glp_minisat_addclause
#define solver_simplify   _glp_minisat_simplify
#define solver_solve      _glp_minisat_solve
#define solver_nvars      _glp_minisat_nvars
#define solver_nclauses   _glp_minisat_nclauses
#define solver_nconflicts _glp_minisat_nconflicts
#define solver_setnvars   _glp_minisat_setnvars
#define solver_propagate  _glp_minisat_propagate
#define solver_reducedb   _glp_minisat_reducedb
#endif

solver* solver_new(void);
void    solver_delete(solver* s);

bool    solver_addclause(solver* s, lit* begin, lit* end);
bool    solver_simplify(solver* s);
bool    solver_solve(solver* s, lit* begin, lit* end);

int     solver_nvars(solver* s);
int     solver_nclauses(solver* s);
int     solver_nconflicts(solver* s);

void    solver_setnvars(solver* s,int n);

#endif

/* eof */

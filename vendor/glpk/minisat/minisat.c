/* minisat.c */

/* Modified by Andrew Makhorin <mao@gnu.org>, August 2011 */
/* May 2017: Changes were made to provide 64-bit portability; thanks to
 * Chris Matrakidis <cmatraki@gmail.com> for patch */

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

#include "env.h"
#include "minisat.h"

#if 1 /* by mao */
static void *ymalloc(int size)
{     void *ptr;
      xassert(size > 0);
      ptr = malloc(size);
      if (ptr == NULL)
         xerror("MiniSat: no memory available\n");
      return ptr;
}

static void *yrealloc(void *ptr, int size)
{     xassert(size > 0);
      if (ptr == NULL)
         ptr = malloc(size);
      else
         ptr = realloc(ptr, size);
      if (ptr == NULL)
         xerror("MiniSat: no memory available\n");
      return ptr;
}

static void yfree(void *ptr)
{     xassert(ptr != NULL);
      free(ptr);
      return;
}

#define assert  xassert
#define printf  xprintf
#define fflush(f) /* nop */
#define malloc  ymalloc
#define realloc yrealloc
#define free    yfree
#define inline  /* empty */
#endif

/*====================================================================*/
/* Debug: */

#if 0
#define VERBOSEDEBUG 1
#endif

/* For derivation output (verbosity level 2) */
#define L_IND    "%-*d"
#define L_ind    solver_dlevel(s)*3+3,solver_dlevel(s)
#define L_LIT    "%sx%d"
#define L_lit(p) lit_sign(p)?"~":"", (lit_var(p))

#if 0 /* by mao */
/* Just like 'assert()' but expression will be evaluated in the release
   version as well. */
static inline void check(int expr) { assert(expr); }
#endif

#if 0 /* by mao */
static void printlits(lit* begin, lit* end)
{
    int i;
    for (i = 0; i < end - begin; i++)
        printf(L_LIT" ",L_lit(begin[i]));
}
#endif

/*====================================================================*/
/* Random numbers: */

/* Returns a random float 0 <= x < 1. Seed must never be 0. */
static inline double drand(double* seed) {
    int q;
    *seed *= 1389796;
    q = (int)(*seed / 2147483647);
    *seed -= (double)q * 2147483647;
    return *seed / 2147483647; }

/* Returns a random integer 0 <= x < size. Seed must never be 0. */
static inline int irand(double* seed, int size) {
    return (int)(drand(seed) * size); }

/*====================================================================*/
/* Predeclarations: */

static void sort(void** array, int size,
    int(*comp)(const void *, const void *));

/*====================================================================*/
/* Clause datatype + minor functions: */

#if 0 /* by mao; see minisat.h */
struct clause_t
{
    int size_learnt;
    lit lits[0];
};
#endif

#define clause_size(c) ((c)->size_learnt >> 1)

#define clause_begin(c) ((c)->lits)

#define clause_learnt(c) ((c)->size_learnt & 1)

#define clause_activity(c) \
    (*((float*)&(c)->lits[(c)->size_learnt>>1]))

#define clause_setactivity(c, a) \
    (void)(*((float*)&(c)->lits[(c)->size_learnt>>1]) = (a))

/*====================================================================*/
/* Encode literals in clause pointers: */

#if 0 /* 8/I-2017 by cmatraki (64-bit portability) */
#define clause_from_lit(l) \
      (clause*)((unsigned long)(l) + (unsigned long)(l) + 1)

#define clause_is_lit(c) \
      ((unsigned long)(c) & 1)

#define clause_read_lit(c) \
      (lit)((unsigned long)(c) >> 1)
#else
#define clause_from_lit(l) \
      (clause*)((size_t)(l) + (size_t)(l) + 1)

#define clause_is_lit(c) \
      ((size_t)(c) & 1)

#define clause_read_lit(c) \
      (lit)((size_t)(c) >> 1)
#endif

/*====================================================================*/
/* Simple helpers: */

#define solver_dlevel(s) \
    (int)veci_size(&(s)->trail_lim)

#define solver_read_wlist(s, l) \
    (vecp *)(&(s)->wlists[l])

static inline void    vecp_remove(vecp* v, void* e)
{
    void** ws = vecp_begin(v);
    int    j  = 0;

    for (; ws[j] != e  ; j++);
    assert(j < vecp_size(v));
    for (; j < vecp_size(v)-1; j++) ws[j] = ws[j+1];
    vecp_resize(v,vecp_size(v)-1);
}

/*====================================================================*/
/* Variable order functions: */

static inline void order_update(solver* s, int v)
{   /* updateorder */
    int*    orderpos = s->orderpos;
    double* activity = s->activity;
    int*    heap     = veci_begin(&s->order);
    int     i        = orderpos[v];
    int     x        = heap[i];
    int     parent   = (i - 1) / 2;

    assert(s->orderpos[v] != -1);

    while (i != 0 && activity[x] > activity[heap[parent]]){
        heap[i]           = heap[parent];
        orderpos[heap[i]] = i;
        i                 = parent;
        parent            = (i - 1) / 2;
    }
    heap[i]     = x;
    orderpos[x] = i;
}

#define order_assigned(s, v) /* nop */

static inline void order_unassigned(solver* s, int v)
{   /* undoorder */
    int* orderpos = s->orderpos;
    if (orderpos[v] == -1){
        orderpos[v] = veci_size(&s->order);
        veci_push(&s->order,v);
        order_update(s,v);
    }
}

static int order_select(solver* s, float random_var_freq)
{   /* selectvar */
    int*    heap;
    double* activity;
    int*    orderpos;

    lbool* values = s->assigns;

    /* Random decision: */
    if (drand(&s->random_seed) < random_var_freq){
        int next = irand(&s->random_seed,s->size);
        assert(next >= 0 && next < s->size);
        if (values[next] == l_Undef)
            return next;
    }

    /* Activity based decision: */

    heap     = veci_begin(&s->order);
    activity = s->activity;
    orderpos = s->orderpos;

    while (veci_size(&s->order) > 0){
        int    next  = heap[0];
        int    size  = veci_size(&s->order)-1;
        int    x     = heap[size];

        veci_resize(&s->order,size);

        orderpos[next] = -1;

        if (size > 0){
            double act   = activity[x];

            int    i     = 0;
            int    child = 1;

            while (child < size){
                if (child+1 < size
                    && activity[heap[child]] < activity[heap[child+1]])
                    child++;

                assert(child < size);

                if (act >= activity[heap[child]])
                    break;

                heap[i]           = heap[child];
                orderpos[heap[i]] = i;
                i                 = child;
                child             = 2 * child + 1;
            }
            heap[i]           = x;
            orderpos[heap[i]] = i;
        }

        if (values[next] == l_Undef)
            return next;
    }

    return var_Undef;
}

/*====================================================================*/
/* Activity functions: */

static inline void act_var_rescale(solver* s) {
    double* activity = s->activity;
    int i;
    for (i = 0; i < s->size; i++)
        activity[i] *= 1e-100;
    s->var_inc *= 1e-100;
}

static inline void act_var_bump(solver* s, int v) {
    double* activity = s->activity;
    if ((activity[v] += s->var_inc) > 1e100)
        act_var_rescale(s);

    /* printf("bump %d %f\n", v-1, activity[v]); */

    if (s->orderpos[v] != -1)
        order_update(s,v);

}

static inline void act_var_decay(solver* s)
    { s->var_inc *= s->var_decay; }

static inline void act_clause_rescale(solver* s) {
    clause** cs = (clause**)vecp_begin(&s->learnts);
    int i;
    for (i = 0; i < vecp_size(&s->learnts); i++){
        float a = clause_activity(cs[i]);
        clause_setactivity(cs[i], a * (float)1e-20);
    }
    s->cla_inc *= (float)1e-20;
}

static inline void act_clause_bump(solver* s, clause *c) {
    float a = clause_activity(c) + s->cla_inc;
    clause_setactivity(c,a);
    if (a > 1e20) act_clause_rescale(s);
}

static inline void act_clause_decay(solver* s)
    { s->cla_inc *= s->cla_decay; }

/*====================================================================*/
/* Clause functions: */

/* pre: size > 1 && no variable occurs twice
 */
static clause* clause_new(solver* s, lit* begin, lit* end, int learnt)
{
    int size;
    clause* c;
    int i;

    assert(end - begin > 1);
    assert(learnt >= 0 && learnt < 2);
    size           = end - begin;
    c              = (clause*)malloc(sizeof(clause)
                     + sizeof(lit) * size + learnt * sizeof(float));
    c->size_learnt = (size << 1) | learnt;
#if 1 /* by mao & cmatraki; non-portable check that is a fundamental  \
       * assumption of minisat code: bit 0 is used as a flag (zero    \
       * for pointer, one for shifted int) so allocated memory should \
       * be at least 16-bit aligned */
    assert(((size_t)c & 1) == 0);
#endif

    for (i = 0; i < size; i++)
        c->lits[i] = begin[i];

    if (learnt)
        *((float*)&c->lits[size]) = 0.0;

    assert(begin[0] >= 0);
    assert(begin[0] < s->size*2);
    assert(begin[1] >= 0);
    assert(begin[1] < s->size*2);

    assert(lit_neg(begin[0]) < s->size*2);
    assert(lit_neg(begin[1]) < s->size*2);

    /* vecp_push(solver_read_wlist(s,lit_neg(begin[0])),(void*)c); */
    /* vecp_push(solver_read_wlist(s,lit_neg(begin[1])),(void*)c); */

    vecp_push(solver_read_wlist(s,lit_neg(begin[0])),
        (void*)(size > 2 ? c : clause_from_lit(begin[1])));
    vecp_push(solver_read_wlist(s,lit_neg(begin[1])),
        (void*)(size > 2 ? c : clause_from_lit(begin[0])));

    return c;
}

static void clause_remove(solver* s, clause* c)
{
    lit* lits = clause_begin(c);
    assert(lit_neg(lits[0]) < s->size*2);
    assert(lit_neg(lits[1]) < s->size*2);

    /* vecp_remove(solver_read_wlist(s,lit_neg(lits[0])),(void*)c); */
    /* vecp_remove(solver_read_wlist(s,lit_neg(lits[1])),(void*)c); */

    assert(lits[0] < s->size*2);
    vecp_remove(solver_read_wlist(s,lit_neg(lits[0])),
        (void*)(clause_size(c) > 2 ? c : clause_from_lit(lits[1])));
    vecp_remove(solver_read_wlist(s,lit_neg(lits[1])),
        (void*)(clause_size(c) > 2 ? c : clause_from_lit(lits[0])));

    if (clause_learnt(c)){
        s->stats.learnts--;
        s->stats.learnts_literals -= clause_size(c);
    }else{
        s->stats.clauses--;
        s->stats.clauses_literals -= clause_size(c);
    }

    free(c);
}

static lbool clause_simplify(solver* s, clause* c)
{
    lit*   lits   = clause_begin(c);
    lbool* values = s->assigns;
    int i;

    assert(solver_dlevel(s) == 0);

    for (i = 0; i < clause_size(c); i++){
        lbool sig = !lit_sign(lits[i]); sig += sig - 1;
        if (values[lit_var(lits[i])] == sig)
            return l_True;
    }
    return l_False;
}

/*====================================================================*/
/* Minor (solver) functions: */

void solver_setnvars(solver* s,int n)
{
    int var;

    if (s->cap < n){

        while (s->cap < n) s->cap = s->cap*2+1;

        s->wlists    = (vecp*)   realloc(s->wlists,
                                 sizeof(vecp)*s->cap*2);
        s->activity  = (double*) realloc(s->activity,
                                 sizeof(double)*s->cap);
        s->assigns   = (lbool*)  realloc(s->assigns,
                                 sizeof(lbool)*s->cap);
        s->orderpos  = (int*)    realloc(s->orderpos,
                                 sizeof(int)*s->cap);
        s->reasons   = (clause**)realloc(s->reasons,
                                 sizeof(clause*)*s->cap);
        s->levels    = (int*)    realloc(s->levels,
                                 sizeof(int)*s->cap);
        s->tags      = (lbool*)  realloc(s->tags,
                                 sizeof(lbool)*s->cap);
        s->trail     = (lit*)    realloc(s->trail,
                                 sizeof(lit)*s->cap);
    }

    for (var = s->size; var < n; var++){
        vecp_new(&s->wlists[2*var]);
        vecp_new(&s->wlists[2*var+1]);
        s->activity [var] = 0;
        s->assigns  [var] = l_Undef;
        s->orderpos [var] = veci_size(&s->order);
        s->reasons  [var] = (clause*)0;
        s->levels   [var] = 0;
        s->tags     [var] = l_Undef;

        /* does not hold because variables enqueued at top level will
           not be reinserted in the heap
           assert(veci_size(&s->order) == var);
         */
        veci_push(&s->order,var);
        order_update(s, var);
    }

    s->size = n > s->size ? n : s->size;
}

static inline bool enqueue(solver* s, lit l, clause* from)
{
    lbool* values = s->assigns;
    int    v      = lit_var(l);
    lbool  val    = values[v];
    lbool  sig;
#ifdef VERBOSEDEBUG
    printf(L_IND"enqueue("L_LIT")\n", L_ind, L_lit(l));
#endif

    /* lbool */ sig = !lit_sign(l); sig += sig - 1;
    if (val != l_Undef){
        return val == sig;
    }else{
        /* New fact -- store it. */
        int*     levels;
        clause** reasons;
#ifdef VERBOSEDEBUG
        printf(L_IND"bind("L_LIT")\n", L_ind, L_lit(l));
#endif
        /* int*     */ levels  = s->levels;
        /* clause** */ reasons = s->reasons;

        values [v] = sig;
        levels [v] = solver_dlevel(s);
        reasons[v] = from;
        s->trail[s->qtail++] = l;

        order_assigned(s, v);
        return true;
    }
}

static inline void assume(solver* s, lit l){
    assert(s->qtail == s->qhead);
    assert(s->assigns[lit_var(l)] == l_Undef);
#ifdef VERBOSEDEBUG
    printf(L_IND"assume("L_LIT")\n", L_ind, L_lit(l));
#endif
    veci_push(&s->trail_lim,s->qtail);
    enqueue(s,l,(clause*)0);
}

static inline void solver_canceluntil(solver* s, int level) {
    lit*     trail;
    lbool*   values;
    clause** reasons;
    int      bound;
    int      c;

    if (solver_dlevel(s) <= level)
        return;

    trail   = s->trail;
    values  = s->assigns;
    reasons = s->reasons;
    bound   = (veci_begin(&s->trail_lim))[level];

    for (c = s->qtail-1; c >= bound; c--) {
        int     x  = lit_var(trail[c]);
        values [x] = l_Undef;
        reasons[x] = (clause*)0;
    }

    for (c = s->qhead-1; c >= bound; c--)
        order_unassigned(s,lit_var(trail[c]));

    s->qhead = s->qtail = bound;
    veci_resize(&s->trail_lim,level);
}

static void solver_record(solver* s, veci* cls)
{
    lit*    begin = veci_begin(cls);
    lit*    end   = begin + veci_size(cls);
    clause* c     = (veci_size(cls) > 1) ? clause_new(s,begin,end,1)
                                         : (clause*)0;
    enqueue(s,*begin,c);

    assert(veci_size(cls) > 0);

    if (c != 0) {
        vecp_push(&s->learnts,c);
        act_clause_bump(s,c);
        s->stats.learnts++;
        s->stats.learnts_literals += veci_size(cls);
    }
}

static double solver_progress(solver* s)
{
    lbool*  values = s->assigns;
    int*    levels = s->levels;
    int     i;

    double  progress = 0;
    double  F        = 1.0 / s->size;
    for (i = 0; i < s->size; i++)
        if (values[i] != l_Undef)
            progress += pow(F, levels[i]);
    return progress / s->size;
}

/*====================================================================*/
/* Major methods: */

static bool solver_lit_removable(solver* s, lit l, int minl)
{
    lbool*   tags    = s->tags;
    clause** reasons = s->reasons;
    int*     levels  = s->levels;
    int      top     = veci_size(&s->tagged);

    assert(lit_var(l) >= 0 && lit_var(l) < s->size);
    assert(reasons[lit_var(l)] != 0);
    veci_resize(&s->stack,0);
    veci_push(&s->stack,lit_var(l));

    while (veci_size(&s->stack) > 0){
        clause* c;
        int v = veci_begin(&s->stack)[veci_size(&s->stack)-1];
        assert(v >= 0 && v < s->size);
        veci_resize(&s->stack,veci_size(&s->stack)-1);
        assert(reasons[v] != 0);
        c    = reasons[v];

        if (clause_is_lit(c)){
            int v = lit_var(clause_read_lit(c));
            if (tags[v] == l_Undef && levels[v] != 0){
                if (reasons[v] != 0
                    && ((1 << (levels[v] & 31)) & minl)){
                    veci_push(&s->stack,v);
                    tags[v] = l_True;
                    veci_push(&s->tagged,v);
                }else{
                    int* tagged = veci_begin(&s->tagged);
                    int j;
                    for (j = top; j < veci_size(&s->tagged); j++)
                        tags[tagged[j]] = l_Undef;
                    veci_resize(&s->tagged,top);
                    return false;
                }
            }
        }else{
            lit*    lits = clause_begin(c);
            int     i, j;

            for (i = 1; i < clause_size(c); i++){
                int v = lit_var(lits[i]);
                if (tags[v] == l_Undef && levels[v] != 0){
                    if (reasons[v] != 0
                        && ((1 << (levels[v] & 31)) & minl)){

                        veci_push(&s->stack,lit_var(lits[i]));
                        tags[v] = l_True;
                        veci_push(&s->tagged,v);
                    }else{
                        int* tagged = veci_begin(&s->tagged);
                        for (j = top; j < veci_size(&s->tagged); j++)
                            tags[tagged[j]] = l_Undef;
                        veci_resize(&s->tagged,top);
                        return false;
                    }
                }
            }
        }
    }

    return true;
}

static void solver_analyze(solver* s, clause* c, veci* learnt)
{
    lit*     trail   = s->trail;
    lbool*   tags    = s->tags;
    clause** reasons = s->reasons;
    int*     levels  = s->levels;
    int      cnt     = 0;
    lit      p       = lit_Undef;
    int      ind     = s->qtail-1;
    lit*     lits;
    int      i, j, minl;
    int*     tagged;

    veci_push(learnt,lit_Undef);

    do{
        assert(c != 0);

        if (clause_is_lit(c)){
            lit q = clause_read_lit(c);
            assert(lit_var(q) >= 0 && lit_var(q) < s->size);
            if (tags[lit_var(q)] == l_Undef && levels[lit_var(q)] > 0){
                tags[lit_var(q)] = l_True;
                veci_push(&s->tagged,lit_var(q));
                act_var_bump(s,lit_var(q));
                if (levels[lit_var(q)] == solver_dlevel(s))
                    cnt++;
                else
                    veci_push(learnt,q);
            }
        }else{

            if (clause_learnt(c))
                act_clause_bump(s,c);

            lits = clause_begin(c);
            /* printlits(lits,lits+clause_size(c)); printf("\n"); */
            for (j = (p == lit_Undef ? 0 : 1); j < clause_size(c); j++){
                lit q = lits[j];
                assert(lit_var(q) >= 0 && lit_var(q) < s->size);
                if (tags[lit_var(q)] == l_Undef
                    && levels[lit_var(q)] > 0){
                    tags[lit_var(q)] = l_True;
                    veci_push(&s->tagged,lit_var(q));
                    act_var_bump(s,lit_var(q));
                    if (levels[lit_var(q)] == solver_dlevel(s))
                        cnt++;
                    else
                        veci_push(learnt,q);
                }
            }
        }

        while (tags[lit_var(trail[ind--])] == l_Undef);

        p = trail[ind+1];
        c = reasons[lit_var(p)];
        cnt--;

    }while (cnt > 0);

    *veci_begin(learnt) = lit_neg(p);

    lits = veci_begin(learnt);
    minl = 0;
    for (i = 1; i < veci_size(learnt); i++){
        int lev = levels[lit_var(lits[i])];
        minl    |= 1 << (lev & 31);
    }

    /* simplify (full) */
    for (i = j = 1; i < veci_size(learnt); i++){
        if (reasons[lit_var(lits[i])] == 0
            || !solver_lit_removable(s,lits[i],minl))
            lits[j++] = lits[i];
    }

    /* update size of learnt + statistics */
    s->stats.max_literals += veci_size(learnt);
    veci_resize(learnt,j);
    s->stats.tot_literals += j;

    /* clear tags */
    tagged = veci_begin(&s->tagged);
    for (i = 0; i < veci_size(&s->tagged); i++)
        tags[tagged[i]] = l_Undef;
    veci_resize(&s->tagged,0);

#ifdef DEBUG
    for (i = 0; i < s->size; i++)
        assert(tags[i] == l_Undef);
#endif

#ifdef VERBOSEDEBUG
    printf(L_IND"Learnt {", L_ind);
    for (i = 0; i < veci_size(learnt); i++)
        printf(" "L_LIT, L_lit(lits[i]));
#endif
    if (veci_size(learnt) > 1){
        int max_i = 1;
        int max   = levels[lit_var(lits[1])];
        lit tmp;

        for (i = 2; i < veci_size(learnt); i++)
            if (levels[lit_var(lits[i])] > max){
                max   = levels[lit_var(lits[i])];
                max_i = i;
            }

        tmp         = lits[1];
        lits[1]     = lits[max_i];
        lits[max_i] = tmp;
    }
#ifdef VERBOSEDEBUG
    {
        int lev = veci_size(learnt) > 1 ? levels[lit_var(lits[1])] : 0;
        printf(" } at level %d\n", lev);
    }
#endif
}

clause* solver_propagate(solver* s)
{
    lbool*  values = s->assigns;
    clause* confl  = (clause*)0;
    lit*    lits;

    /* printf("solver_propagate\n"); */
    while (confl == 0 && s->qtail - s->qhead > 0){
        lit  p  = s->trail[s->qhead++];
        vecp* ws = solver_read_wlist(s,p);
        clause **begin = (clause**)vecp_begin(ws);
        clause **end   = begin + vecp_size(ws);
        clause **i, **j;

        s->stats.propagations++;
        s->simpdb_props--;

        /* printf("checking lit %d: "L_LIT"\n", veci_size(ws),
               L_lit(p)); */
        for (i = j = begin; i < end; ){
            if (clause_is_lit(*i)){
                *j++ = *i;
                if (!enqueue(s,clause_read_lit(*i),clause_from_lit(p))){
                    confl = s->binary;
                    (clause_begin(confl))[1] = lit_neg(p);
                    (clause_begin(confl))[0] = clause_read_lit(*i++);

                    /* Copy the remaining watches: */
                    while (i < end)
                        *j++ = *i++;
                }
            }else{
                lit false_lit;
                lbool sig;

                lits = clause_begin(*i);

                /* Make sure the false literal is data[1]: */
                false_lit = lit_neg(p);
                if (lits[0] == false_lit){
                    lits[0] = lits[1];
                    lits[1] = false_lit;
                }
                assert(lits[1] == false_lit);
                /* printf("checking clause: ");
                   printlits(lits, lits+clause_size(*i));
                   printf("\n"); */

                /* If 0th watch is true, then clause is already
                   satisfied. */
                sig = !lit_sign(lits[0]); sig += sig - 1;
                if (values[lit_var(lits[0])] == sig){
                    *j++ = *i;
                }else{
                    /* Look for new watch: */
                    lit* stop = lits + clause_size(*i);
                    lit* k;
                    for (k = lits + 2; k < stop; k++){
                        lbool sig = lit_sign(*k); sig += sig - 1;
                        if (values[lit_var(*k)] != sig){
                            lits[1] = *k;
                            *k = false_lit;
                            vecp_push(solver_read_wlist(s,
                                lit_neg(lits[1])),*i);
                            goto next; }
                    }

                    *j++ = *i;
                    /* Clause is unit under assignment: */
                    if (!enqueue(s,lits[0], *i)){
                        confl = *i++;
                        /* Copy the remaining watches: */
                        while (i < end)
                            *j++ = *i++;
                    }
                }
            }
        next:
            i++;
        }

        s->stats.inspects += j - (clause**)vecp_begin(ws);
        vecp_resize(ws,j - (clause**)vecp_begin(ws));
    }

    return confl;
}

static inline int clause_cmp (const void* x, const void* y) {
    return clause_size((clause*)x) > 2
           && (clause_size((clause*)y) == 2
               || clause_activity((clause*)x)
                  < clause_activity((clause*)y)) ? -1 : 1; }

void solver_reducedb(solver* s)
{
    int      i, j;
    double   extra_lim = s->cla_inc / vecp_size(&s->learnts);
             /* Remove any clause below this activity */
    clause** learnts = (clause**)vecp_begin(&s->learnts);
    clause** reasons = s->reasons;

    sort(vecp_begin(&s->learnts), vecp_size(&s->learnts), clause_cmp);

    for (i = j = 0; i < vecp_size(&s->learnts) / 2; i++){
        if (clause_size(learnts[i]) > 2
            && reasons[lit_var(*clause_begin(learnts[i]))]
               != learnts[i])
            clause_remove(s,learnts[i]);
        else
            learnts[j++] = learnts[i];
    }
    for (; i < vecp_size(&s->learnts); i++){
        if (clause_size(learnts[i]) > 2
            && reasons[lit_var(*clause_begin(learnts[i]))]
               != learnts[i]
            && clause_activity(learnts[i]) < extra_lim)
            clause_remove(s,learnts[i]);
        else
            learnts[j++] = learnts[i];
    }

    /* printf("reducedb deleted %d\n", vecp_size(&s->learnts) - j); */

    vecp_resize(&s->learnts,j);
}

static lbool solver_search(solver* s, int nof_conflicts,
                           int nof_learnts)
{
    int*    levels          = s->levels;
    double  var_decay       = 0.95;
    double  clause_decay    = 0.999;
    double  random_var_freq = 0.02;

    int     conflictC       = 0;
    veci    learnt_clause;

    assert(s->root_level == solver_dlevel(s));

    s->stats.starts++;
    s->var_decay = (float)(1 / var_decay   );
    s->cla_decay = (float)(1 / clause_decay);
    veci_resize(&s->model,0);
    veci_new(&learnt_clause);

    for (;;){
        clause* confl = solver_propagate(s);
        if (confl != 0){
            /* CONFLICT */
            int blevel;

#ifdef VERBOSEDEBUG
            printf(L_IND"**CONFLICT**\n", L_ind);
#endif
            s->stats.conflicts++; conflictC++;
            if (solver_dlevel(s) == s->root_level){
                veci_delete(&learnt_clause);
                return l_False;
            }

            veci_resize(&learnt_clause,0);
            solver_analyze(s, confl, &learnt_clause);
            blevel = veci_size(&learnt_clause) > 1
                     ? levels[lit_var(veci_begin(&learnt_clause)[1])]
                     : s->root_level;
            blevel = s->root_level > blevel ? s->root_level : blevel;
            solver_canceluntil(s,blevel);
            solver_record(s,&learnt_clause);
            act_var_decay(s);
            act_clause_decay(s);

        }else{
            /* NO CONFLICT */
            int next;

            if (nof_conflicts >= 0 && conflictC >= nof_conflicts){
                /* Reached bound on number of conflicts: */
                s->progress_estimate = solver_progress(s);
                solver_canceluntil(s,s->root_level);
                veci_delete(&learnt_clause);
                return l_Undef; }

            if (solver_dlevel(s) == 0)
                /* Simplify the set of problem clauses: */
                solver_simplify(s);

            if (nof_learnts >= 0
                && vecp_size(&s->learnts) - s->qtail >= nof_learnts)
                /* Reduce the set of learnt clauses: */
                solver_reducedb(s);

            /* New variable decision: */
            s->stats.decisions++;
            next = order_select(s,(float)random_var_freq);

            if (next == var_Undef){
                /* Model found: */
                lbool* values = s->assigns;
                int i;
                for (i = 0; i < s->size; i++)
                    veci_push(&s->model,(int)values[i]);
                solver_canceluntil(s,s->root_level);
                veci_delete(&learnt_clause);

                /*
                veci apa; veci_new(&apa);
                for (i = 0; i < s->size; i++)
                    veci_push(&apa,(int)(s->model.ptr[i] == l_True
                        ? toLit(i) : lit_neg(toLit(i))));
                printf("model: ");
                printlits((lit*)apa.ptr,
                    (lit*)apa.ptr + veci_size(&apa)); printf("\n");
                veci_delete(&apa);
                */

                return l_True;
            }

            assume(s,lit_neg(toLit(next)));
        }
    }

#if 0 /* by mao; unreachable code */
    return l_Undef; /* cannot happen */
#endif
}

/*====================================================================*/
/* External solver functions: */

solver* solver_new(void)
{
    solver* s = (solver*)malloc(sizeof(solver));

    /* initialize vectors */
    vecp_new(&s->clauses);
    vecp_new(&s->learnts);
    veci_new(&s->order);
    veci_new(&s->trail_lim);
    veci_new(&s->tagged);
    veci_new(&s->stack);
    veci_new(&s->model);

    /* initialize arrays */
    s->wlists    = 0;
    s->activity  = 0;
    s->assigns   = 0;
    s->orderpos  = 0;
    s->reasons   = 0;
    s->levels    = 0;
    s->tags      = 0;
    s->trail     = 0;

    /* initialize other vars */
    s->size                   = 0;
    s->cap                    = 0;
    s->qhead                  = 0;
    s->qtail                  = 0;
    s->cla_inc                = 1;
    s->cla_decay              = 1;
    s->var_inc                = 1;
    s->var_decay              = 1;
    s->root_level             = 0;
    s->simpdb_assigns         = 0;
    s->simpdb_props           = 0;
    s->random_seed            = 91648253;
    s->progress_estimate      = 0;
    s->binary                 = (clause*)malloc(sizeof(clause)
                                                + sizeof(lit)*2);
    s->binary->size_learnt    = (2 << 1);
    s->verbosity              = 0;

    s->stats.starts           = 0;
    s->stats.decisions        = 0;
    s->stats.propagations     = 0;
    s->stats.inspects         = 0;
    s->stats.conflicts        = 0;
    s->stats.clauses          = 0;
    s->stats.clauses_literals = 0;
    s->stats.learnts          = 0;
    s->stats.learnts_literals = 0;
    s->stats.max_literals     = 0;
    s->stats.tot_literals     = 0;

    return s;
}

void solver_delete(solver* s)
{
    int i;
    for (i = 0; i < vecp_size(&s->clauses); i++)
        free(vecp_begin(&s->clauses)[i]);

    for (i = 0; i < vecp_size(&s->learnts); i++)
        free(vecp_begin(&s->learnts)[i]);

    /* delete vectors */
    vecp_delete(&s->clauses);
    vecp_delete(&s->learnts);
    veci_delete(&s->order);
    veci_delete(&s->trail_lim);
    veci_delete(&s->tagged);
    veci_delete(&s->stack);
    veci_delete(&s->model);
    free(s->binary);

    /* delete arrays */
    if (s->wlists != 0){
        int i;
        for (i = 0; i < s->size*2; i++)
            vecp_delete(&s->wlists[i]);

        /* if one is different from null, all are */
        free(s->wlists);
        free(s->activity );
        free(s->assigns  );
        free(s->orderpos );
        free(s->reasons  );
        free(s->levels   );
        free(s->trail    );
        free(s->tags     );
    }

    free(s);
}

bool solver_addclause(solver* s, lit* begin, lit* end)
{
    lit *i,*j;
    int maxvar;
    lbool* values;
    lit last;

    if (begin == end) return false;

    /* printlits(begin,end); printf("\n"); */
    /* insertion sort */
    maxvar = lit_var(*begin);
    for (i = begin + 1; i < end; i++){
        lit l = *i;
        maxvar = lit_var(l) > maxvar ? lit_var(l) : maxvar;
        for (j = i; j > begin && *(j-1) > l; j--)
            *j = *(j-1);
        *j = l;
    }
    solver_setnvars(s,maxvar+1);

    /* printlits(begin,end); printf("\n"); */
    values = s->assigns;

    /* delete duplicates */
    last = lit_Undef;
    for (i = j = begin; i < end; i++){
        /* printf("lit: "L_LIT", value = %d\n", L_lit(*i),
        (lit_sign(*i) ? -values[lit_var(*i)] : values[lit_var(*i)])); */
        lbool sig = !lit_sign(*i); sig += sig - 1;
        if (*i == lit_neg(last) || sig == values[lit_var(*i)])
            return true;   /* tautology */
        else if (*i != last && values[lit_var(*i)] == l_Undef)
            last = *j++ = *i;
    }

    /* printf("final: "); printlits(begin,j); printf("\n"); */

    if (j == begin)          /* empty clause */
        return false;
    else if (j - begin == 1) /* unit clause */
        return enqueue(s,*begin,(clause*)0);

    /* create new clause */
    vecp_push(&s->clauses,clause_new(s,begin,j,0));

    s->stats.clauses++;
    s->stats.clauses_literals += j - begin;

    return true;
}

bool   solver_simplify(solver* s)
{
    clause** reasons;
    int type;

    assert(solver_dlevel(s) == 0);

    if (solver_propagate(s) != 0)
        return false;

    if (s->qhead == s->simpdb_assigns || s->simpdb_props > 0)
        return true;

    reasons = s->reasons;
    for (type = 0; type < 2; type++){
        vecp*    cs  = type ? &s->learnts : &s->clauses;
        clause** cls = (clause**)vecp_begin(cs);

        int i, j;
        for (j = i = 0; i < vecp_size(cs); i++){
            if (reasons[lit_var(*clause_begin(cls[i]))] != cls[i] &&
                clause_simplify(s,cls[i]) == l_True)
                clause_remove(s,cls[i]);
            else
                cls[j++] = cls[i];
        }
        vecp_resize(cs,j);
    }

    s->simpdb_assigns = s->qhead;
    /* (shouldn't depend on 'stats' really, but it will do for now) */
    s->simpdb_props   = (int)(s->stats.clauses_literals
                              + s->stats.learnts_literals);

    return true;
}

bool   solver_solve(solver* s, lit* begin, lit* end)
{
    double  nof_conflicts = 100;
    double  nof_learnts   = solver_nclauses(s) / 3;
    lbool   status        = l_Undef;
    lbool*  values        = s->assigns;
    lit*    i;

    /* printf("solve: "); printlits(begin, end); printf("\n"); */
    for (i = begin; i < end; i++){
        switch (lit_sign(*i) ? -values[lit_var(*i)]
                             : values[lit_var(*i)]){
        case 1: /* l_True: */
            break;
        case 0: /* l_Undef */
            assume(s, *i);
            if (solver_propagate(s) == NULL)
                break;
            /* falltrough */
        case -1: /* l_False */
            solver_canceluntil(s, 0);
            return false;
        }
    }

    s->root_level = solver_dlevel(s);

    if (s->verbosity >= 1){
        printf("==================================[MINISAT]============"
               "=======================\n");
        printf("| Conflicts |     ORIGINAL     |              LEARNT   "
               "           | Progress |\n");
        printf("|           | Clauses Literals |   Limit Clauses Litera"
               "ls  Lit/Cl |          |\n");
        printf("======================================================="
               "=======================\n");
    }

    while (status == l_Undef){
        double Ratio = (s->stats.learnts == 0)? 0.0 :
            s->stats.learnts_literals / (double)s->stats.learnts;

        if (s->verbosity >= 1){
            printf("| %9.0f | %7.0f %8.0f | %7.0f %7.0f %8.0f %7.1f | %"
                   "6.3f %% |\n",
                (double)s->stats.conflicts,
                (double)s->stats.clauses,
                (double)s->stats.clauses_literals,
                (double)nof_learnts,
                (double)s->stats.learnts,
                (double)s->stats.learnts_literals,
                Ratio,
                s->progress_estimate*100);
            fflush(stdout);
        }
        status = solver_search(s,(int)nof_conflicts, (int)nof_learnts);
        nof_conflicts *= 1.5;
        nof_learnts   *= 1.1;
    }
    if (s->verbosity >= 1)
        printf("======================================================="
               "=======================\n");

    solver_canceluntil(s,0);
    return status != l_False;
}

int solver_nvars(solver* s)
{
    return s->size;
}

int solver_nclauses(solver* s)
{
    return vecp_size(&s->clauses);
}

int solver_nconflicts(solver* s)
{
    return (int)s->stats.conflicts;
}

/*====================================================================*/
/* Sorting functions (sigh): */

static inline void selectionsort(void** array, int size,
                                 int(*comp)(const void *, const void *))
{
    int     i, j, best_i;
    void*   tmp;

    for (i = 0; i < size-1; i++){
        best_i = i;
        for (j = i+1; j < size; j++){
            if (comp(array[j], array[best_i]) < 0)
                best_i = j;
        }
        tmp = array[i]; array[i] = array[best_i]; array[best_i] = tmp;
    }
}

static void sortrnd(void** array, int size,
                    int(*comp)(const void *, const void *),
                    double* seed)
{
    if (size <= 15)
        selectionsort(array, size, comp);

    else{
        void*       pivot = array[irand(seed, size)];
        void*       tmp;
        int         i = -1;
        int         j = size;

        for(;;){
            do i++; while(comp(array[i], pivot)<0);
            do j--; while(comp(pivot, array[j])<0);

            if (i >= j) break;

            tmp = array[i]; array[i] = array[j]; array[j] = tmp;
        }

        sortrnd(array    , i     , comp, seed);
        sortrnd(&array[i], size-i, comp, seed);
    }
}

static void sort(void** array, int size,
          int(*comp)(const void *, const void *))
{
    double seed = 91648253;
    sortrnd(array,size,comp,&seed);
}

/* eof */

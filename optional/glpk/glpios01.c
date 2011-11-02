/* glpios01.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*
*  Copyright (C) 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008,
*  2009, 2010 Andrew Makhorin, Department for Applied Informatics,
*  Moscow Aviation Institute, Moscow, Russia. All rights reserved.
*  E-mail: <mao@gnu.org>.
*
*  GLPK is free software: you can redistribute it and/or modify it
*  under the terms of the GNU General Public License as published by
*  the Free Software Foundation, either version 3 of the License, or
*  (at your option) any later version.
*
*  GLPK is distributed in the hope that it will be useful, but WITHOUT
*  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
*  or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public
*  License for more details.
*
*  You should have received a copy of the GNU General Public License
*  along with GLPK. If not, see <http://www.gnu.org/licenses/>.
***********************************************************************/

#include "glpios.h"

/***********************************************************************
*  NAME
*
*  ios_create_tree - create branch-and-bound tree
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  glp_tree *ios_create_tree(glp_prob *mip, const glp_iocp *parm);
*
*  DESCRIPTION
*
*  The routine ios_create_tree creates the branch-and-bound tree.
*
*  Being created the tree consists of the only root subproblem whose
*  reference number is 1. Note that initially the root subproblem is in
*  frozen state and therefore needs to be revived.
*
*  RETURNS
*
*  The routine returns a pointer to the tree created. */

static IOSNPD *new_node(glp_tree *tree, IOSNPD *parent);

glp_tree *ios_create_tree(glp_prob *mip, const glp_iocp *parm)
{     int m = mip->m;
      int n = mip->n;
      glp_tree *tree;
      int i, j;
      xassert(mip->tree == NULL);
      mip->tree = tree = xmalloc(sizeof(glp_tree));
      tree->pool = dmp_create_pool();
      tree->n = n;
      /* save original problem components */
      tree->orig_m = m;
      tree->orig_type = xcalloc(1+m+n, sizeof(char));
      tree->orig_lb = xcalloc(1+m+n, sizeof(double));
      tree->orig_ub = xcalloc(1+m+n, sizeof(double));
      tree->orig_stat = xcalloc(1+m+n, sizeof(char));
      tree->orig_prim = xcalloc(1+m+n, sizeof(double));
      tree->orig_dual = xcalloc(1+m+n, sizeof(double));
      for (i = 1; i <= m; i++)
      {  GLPROW *row = mip->row[i];
         tree->orig_type[i] = (char)row->type;
         tree->orig_lb[i] = row->lb;
         tree->orig_ub[i] = row->ub;
         tree->orig_stat[i] = (char)row->stat;
         tree->orig_prim[i] = row->prim;
         tree->orig_dual[i] = row->dual;
      }
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = mip->col[j];
         tree->orig_type[m+j] = (char)col->type;
         tree->orig_lb[m+j] = col->lb;
         tree->orig_ub[m+j] = col->ub;
         tree->orig_stat[m+j] = (char)col->stat;
         tree->orig_prim[m+j] = col->prim;
         tree->orig_dual[m+j] = col->dual;
      }
      tree->orig_obj = mip->obj_val;
      /* initialize the branch-and-bound tree */
      tree->nslots = 0;
      tree->avail = 0;
      tree->slot = NULL;
      tree->head = tree->tail = NULL;
      tree->a_cnt = tree->n_cnt = tree->t_cnt = 0;
      /* the root subproblem is not solved yet, so its final components
         are unknown so far */
      tree->root_m = 0;
      tree->root_type = NULL;
      tree->root_lb = tree->root_ub = NULL;
      tree->root_stat = NULL;
      /* the current subproblem does not exist yet */
      tree->curr = NULL;
      tree->mip = mip;
      /*tree->solved = 0;*/
      tree->non_int = xcalloc(1+n, sizeof(char));
      memset(&tree->non_int[1], 0, n);
      /* arrays to save parent subproblem components will be allocated
         later */
      tree->pred_m = tree->pred_max = 0;
      tree->pred_type = NULL;
      tree->pred_lb = tree->pred_ub = NULL;
      tree->pred_stat = NULL;
      /* cut generator */
      tree->local = ios_create_pool(tree);
      /*tree->first_attempt = 1;*/
      /*tree->max_added_cuts = 0;*/
      /*tree->min_eff = 0.0;*/
      /*tree->miss = 0;*/
      /*tree->just_selected = 0;*/
      tree->mir_gen = NULL;
      tree->clq_gen = NULL;
      /*tree->round = 0;*/
#if 0
      /* create the conflict graph */
      tree->n_ref = xcalloc(1+n, sizeof(int));
      memset(&tree->n_ref[1], 0, n * sizeof(int));
      tree->c_ref = xcalloc(1+n, sizeof(int));
      memset(&tree->c_ref[1], 0, n * sizeof(int));
      tree->g = scg_create_graph(0);
      tree->j_ref = xcalloc(1+tree->g->n_max, sizeof(int));
#endif
      /* pseudocost branching */
      tree->pcost = NULL;
      tree->iwrk = xcalloc(1+n, sizeof(int));
      tree->dwrk = xcalloc(1+n, sizeof(double));
      /* initialize control parameters */
      tree->parm = parm;
      tree->tm_beg = xtime();
      tree->tm_lag = xlset(0);
      tree->sol_cnt = 0;
      /* initialize advanced solver interface */
      tree->reason = 0;
      tree->reopt = 0;
      tree->reinv = 0;
      tree->br_var = 0;
      tree->br_sel = 0;
      tree->child = 0;
      tree->next_p = 0;
      /*tree->btrack = NULL;*/
      tree->stop = 0;
      /* create the root subproblem, which initially is identical to
         the original MIP */
      new_node(tree, NULL);
      return tree;
}

/***********************************************************************
*  NAME
*
*  ios_revive_node - revive specified subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_revive_node(glp_tree *tree, int p);
*
*  DESCRIPTION
*
*  The routine ios_revive_node revives the specified subproblem, whose
*  reference number is p, and thereby makes it the current subproblem.
*  Note that the specified subproblem must be active. Besides, if the
*  current subproblem already exists, it must be frozen before reviving
*  another subproblem. */

void ios_revive_node(glp_tree *tree, int p)
{     glp_prob *mip = tree->mip;
      IOSNPD *node, *root;
      /* obtain pointer to the specified subproblem */
      xassert(1 <= p && p <= tree->nslots);
      node = tree->slot[p].node;
      xassert(node != NULL);
      /* the specified subproblem must be active */
      xassert(node->count == 0);
      /* the current subproblem must not exist */
      xassert(tree->curr == NULL);
      /* the specified subproblem becomes current */
      tree->curr = node;
      /*tree->solved = 0;*/
      /* obtain pointer to the root subproblem */
      root = tree->slot[1].node;
      xassert(root != NULL);
      /* at this point problem object components correspond to the root
         subproblem, so if the root subproblem should be revived, there
         is nothing more to do */
      if (node == root) goto done;
      xassert(mip->m == tree->root_m);
      /* build path from the root to the current node */
      node->temp = NULL;
      for (node = node; node != NULL; node = node->up)
      {  if (node->up == NULL)
            xassert(node == root);
         else
            node->up->temp = node;
      }
      /* go down from the root to the current node and make necessary
         changes to restore components of the current subproblem */
      for (node = root; node != NULL; node = node->temp)
      {  int m = mip->m;
         int n = mip->n;
         /* if the current node is reached, the problem object at this
            point corresponds to its parent, so save attributes of rows
            and columns for the parent subproblem */
         if (node->temp == NULL)
         {  int i, j;
            tree->pred_m = m;
            /* allocate/reallocate arrays, if necessary */
            if (tree->pred_max < m + n)
            {  int new_size = m + n + 100;
               if (tree->pred_type != NULL) xfree(tree->pred_type);
               if (tree->pred_lb != NULL) xfree(tree->pred_lb);
               if (tree->pred_ub != NULL) xfree(tree->pred_ub);
               if (tree->pred_stat != NULL) xfree(tree->pred_stat);
               tree->pred_max = new_size;
               tree->pred_type = xcalloc(1+new_size, sizeof(char));
               tree->pred_lb = xcalloc(1+new_size, sizeof(double));
               tree->pred_ub = xcalloc(1+new_size, sizeof(double));
               tree->pred_stat = xcalloc(1+new_size, sizeof(char));
            }
            /* save row attributes */
            for (i = 1; i <= m; i++)
            {  GLPROW *row = mip->row[i];
               tree->pred_type[i] = (char)row->type;
               tree->pred_lb[i] = row->lb;
               tree->pred_ub[i] = row->ub;
               tree->pred_stat[i] = (char)row->stat;
            }
            /* save column attributes */
            for (j = 1; j <= n; j++)
            {  GLPCOL *col = mip->col[j];
               tree->pred_type[mip->m+j] = (char)col->type;
               tree->pred_lb[mip->m+j] = col->lb;
               tree->pred_ub[mip->m+j] = col->ub;
               tree->pred_stat[mip->m+j] = (char)col->stat;
            }
         }
         /* change bounds of rows and columns */
         {  IOSBND *b;
            for (b = node->b_ptr; b != NULL; b = b->next)
            {  if (b->k <= m)
                  glp_set_row_bnds(mip, b->k, b->type, b->lb, b->ub);
               else
                  glp_set_col_bnds(mip, b->k-m, b->type, b->lb, b->ub);
            }
         }
         /* change statuses of rows and columns */
         {  IOSTAT *s;
            for (s = node->s_ptr; s != NULL; s = s->next)
            {  if (s->k <= m)
                  glp_set_row_stat(mip, s->k, s->stat);
               else
                  glp_set_col_stat(mip, s->k-m, s->stat);
            }
         }
         /* add new rows */
         if (node->r_ptr != NULL)
         {  IOSROW *r;
            IOSAIJ *a;
            int i, len, *ind;
            double *val;
            ind = xcalloc(1+n, sizeof(int));
            val = xcalloc(1+n, sizeof(double));
            for (r = node->r_ptr; r != NULL; r = r->next)
            {  i = glp_add_rows(mip, 1);
               glp_set_row_name(mip, i, r->name);
#if 1 /* 20/IX-2008 */
               xassert(mip->row[i]->level == 0);
               mip->row[i]->level = node->level;
               mip->row[i]->origin = r->origin;
               mip->row[i]->klass = r->klass;
#endif
               glp_set_row_bnds(mip, i, r->type, r->lb, r->ub);
               len = 0;
               for (a = r->ptr; a != NULL; a = a->next)
                  len++, ind[len] = a->j, val[len] = a->val;
               glp_set_mat_row(mip, i, len, ind, val);
               glp_set_rii(mip, i, r->rii);
               glp_set_row_stat(mip, i, r->stat);
            }
            xfree(ind);
            xfree(val);
         }
#if 0
         /* add new edges to the conflict graph */
         /* add new cliques to the conflict graph */
         /* (not implemented yet) */
         xassert(node->own_nn == 0);
         xassert(node->own_nc == 0);
         xassert(node->e_ptr == NULL);
#endif
      }
      /* the specified subproblem has been revived */
      node = tree->curr;
      /* delete its bound change list */
      while (node->b_ptr != NULL)
      {  IOSBND *b;
         b = node->b_ptr;
         node->b_ptr = b->next;
         dmp_free_atom(tree->pool, b, sizeof(IOSBND));
      }
      /* delete its status change list */
      while (node->s_ptr != NULL)
      {  IOSTAT *s;
         s = node->s_ptr;
         node->s_ptr = s->next;
         dmp_free_atom(tree->pool, s, sizeof(IOSTAT));
      }
#if 1 /* 20/XI-2009 */
      /* delete its row addition list (additional rows may appear, for
         example, due to branching on GUB constraints */
      while (node->r_ptr != NULL)
      {  IOSROW *r;
         r = node->r_ptr;
         node->r_ptr = r->next;
         xassert(r->name == NULL);
         while (r->ptr != NULL)
         {  IOSAIJ *a;
            a = r->ptr;
            r->ptr = a->next;
            dmp_free_atom(tree->pool, a, sizeof(IOSAIJ));
         }
         dmp_free_atom(tree->pool, r, sizeof(IOSROW));
      }
#endif
done: return;
}

/***********************************************************************
*  NAME
*
*  ios_freeze_node - freeze current subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_freeze_node(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_freeze_node freezes the current subproblem. */

void ios_freeze_node(glp_tree *tree)
{     glp_prob *mip = tree->mip;
      int m = mip->m;
      int n = mip->n;
      IOSNPD *node;
      /* obtain pointer to the current subproblem */
      node = tree->curr;
      xassert(node != NULL);
      if (node->up == NULL)
      {  /* freeze the root subproblem */
         int k;
         xassert(node->p == 1);
         xassert(tree->root_m == 0);
         xassert(tree->root_type == NULL);
         xassert(tree->root_lb == NULL);
         xassert(tree->root_ub == NULL);
         xassert(tree->root_stat == NULL);
         tree->root_m = m;
         tree->root_type = xcalloc(1+m+n, sizeof(char));
         tree->root_lb = xcalloc(1+m+n, sizeof(double));
         tree->root_ub = xcalloc(1+m+n, sizeof(double));
         tree->root_stat = xcalloc(1+m+n, sizeof(char));
         for (k = 1; k <= m+n; k++)
         {  if (k <= m)
            {  GLPROW *row = mip->row[k];
               tree->root_type[k] = (char)row->type;
               tree->root_lb[k] = row->lb;
               tree->root_ub[k] = row->ub;
               tree->root_stat[k] = (char)row->stat;
            }
            else
            {  GLPCOL *col = mip->col[k-m];
               tree->root_type[k] = (char)col->type;
               tree->root_lb[k] = col->lb;
               tree->root_ub[k] = col->ub;
               tree->root_stat[k] = (char)col->stat;
            }
         }
      }
      else
      {  /* freeze non-root subproblem */
         int root_m = tree->root_m;
         int pred_m = tree->pred_m;
         int i, j, k;
         xassert(pred_m <= m);
         /* build change lists for rows and columns which exist in the
            parent subproblem */
         xassert(node->b_ptr == NULL);
         xassert(node->s_ptr == NULL);
         for (k = 1; k <= pred_m + n; k++)
         {  int pred_type, pred_stat, type, stat;
            double pred_lb, pred_ub, lb, ub;
            /* determine attributes in the parent subproblem */
            pred_type = tree->pred_type[k];
            pred_lb = tree->pred_lb[k];
            pred_ub = tree->pred_ub[k];
            pred_stat = tree->pred_stat[k];
            /* determine attributes in the current subproblem */
            if (k <= pred_m)
            {  GLPROW *row = mip->row[k];
               type = row->type;
               lb = row->lb;
               ub = row->ub;
               stat = row->stat;
            }
            else
            {  GLPCOL *col = mip->col[k - pred_m];
               type = col->type;
               lb = col->lb;
               ub = col->ub;
               stat = col->stat;
            }
            /* save type and bounds of a row/column, if changed */
            if (!(pred_type == type && pred_lb == lb && pred_ub == ub))
            {  IOSBND *b;
               b = dmp_get_atom(tree->pool, sizeof(IOSBND));
               b->k = k;
               b->type = (unsigned char)type;
               b->lb = lb;
               b->ub = ub;
               b->next = node->b_ptr;
               node->b_ptr = b;
            }
            /* save status of a row/column, if changed */
            if (pred_stat != stat)
            {  IOSTAT *s;
               s = dmp_get_atom(tree->pool, sizeof(IOSTAT));
               s->k = k;
               s->stat = (unsigned char)stat;
               s->next = node->s_ptr;
               node->s_ptr = s;
            }
         }
         /* save new rows added to the current subproblem */
         xassert(node->r_ptr == NULL);
         if (pred_m < m)
         {  int i, len, *ind;
            double *val;
            ind = xcalloc(1+n, sizeof(int));
            val = xcalloc(1+n, sizeof(double));
            for (i = m; i > pred_m; i--)
            {  GLPROW *row = mip->row[i];
               IOSROW *r;
               const char *name;
               r = dmp_get_atom(tree->pool, sizeof(IOSROW));
               name = glp_get_row_name(mip, i);
               if (name == NULL)
                  r->name = NULL;
               else
               {  r->name = dmp_get_atom(tree->pool, strlen(name)+1);
                  strcpy(r->name, name);
               }
#if 1 /* 20/IX-2008 */
               r->origin = row->origin;
               r->klass = row->klass;
#endif
               r->type = (unsigned char)row->type;
               r->lb = row->lb;
               r->ub = row->ub;
               r->ptr = NULL;
               len = glp_get_mat_row(mip, i, ind, val);
               for (k = 1; k <= len; k++)
               {  IOSAIJ *a;
                  a = dmp_get_atom(tree->pool, sizeof(IOSAIJ));
                  a->j = ind[k];
                  a->val = val[k];
                  a->next = r->ptr;
                  r->ptr = a;
               }
               r->rii = row->rii;
               r->stat = (unsigned char)row->stat;
               r->next = node->r_ptr;
               node->r_ptr = r;
            }
            xfree(ind);
            xfree(val);
         }
         /* remove all rows missing in the root subproblem */
         if (m != root_m)
         {  int nrs, *num;
            nrs = m - root_m;
            xassert(nrs > 0);
            num = xcalloc(1+nrs, sizeof(int));
            for (i = 1; i <= nrs; i++) num[i] = root_m + i;
            glp_del_rows(mip, nrs, num);
            xfree(num);
         }
         m = mip->m;
         /* and restore attributes of all rows and columns for the root
            subproblem */
         xassert(m == root_m);
         for (i = 1; i <= m; i++)
         {  glp_set_row_bnds(mip, i, tree->root_type[i],
               tree->root_lb[i], tree->root_ub[i]);
            glp_set_row_stat(mip, i, tree->root_stat[i]);
         }
         for (j = 1; j <= n; j++)
         {  glp_set_col_bnds(mip, j, tree->root_type[m+j],
               tree->root_lb[m+j], tree->root_ub[m+j]);
            glp_set_col_stat(mip, j, tree->root_stat[m+j]);
         }
#if 1
         /* remove all edges and cliques missing in the conflict graph
            for the root subproblem */
         /* (not implemented yet) */
#endif
      }
      /* the current subproblem has been frozen */
      tree->curr = NULL;
      return;
}

/***********************************************************************
*  NAME
*
*  ios_clone_node - clone specified subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_clone_node(glp_tree *tree, int p, int nnn, int ref[]);
*
*  DESCRIPTION
*
*  The routine ios_clone_node clones the specified subproblem, whose
*  reference number is p, creating its nnn exact copies. Note that the
*  specified subproblem must be active and must be in the frozen state
*  (i.e. it must not be the current subproblem).
*
*  Each clone, an exact copy of the specified subproblem, becomes a new
*  active subproblem added to the end of the active list. After cloning
*  the specified subproblem becomes inactive.
*
*  The reference numbers of clone subproblems are stored to locations
*  ref[1], ..., ref[nnn]. */

static int get_slot(glp_tree *tree)
{     int p;
      /* if no free slots are available, increase the room */
      if (tree->avail == 0)
      {  int nslots = tree->nslots;
         IOSLOT *save = tree->slot;
         if (nslots == 0)
            tree->nslots = 20;
         else
         {  tree->nslots = nslots + nslots;
            xassert(tree->nslots > nslots);
         }
         tree->slot = xcalloc(1+tree->nslots, sizeof(IOSLOT));
         if (save != NULL)
         {  memcpy(&tree->slot[1], &save[1], nslots * sizeof(IOSLOT));
            xfree(save);
         }
         /* push more free slots into the stack */
         for (p = tree->nslots; p > nslots; p--)
         {  tree->slot[p].node = NULL;
            tree->slot[p].next = tree->avail;
            tree->avail = p;
         }
      }
      /* pull a free slot from the stack */
      p = tree->avail;
      tree->avail = tree->slot[p].next;
      xassert(tree->slot[p].node == NULL);
      tree->slot[p].next = 0;
      return p;
}

static IOSNPD *new_node(glp_tree *tree, IOSNPD *parent)
{     IOSNPD *node;
      int p;
      /* pull a free slot for the new node */
      p = get_slot(tree);
      /* create descriptor of the new subproblem */
      node = dmp_get_atom(tree->pool, sizeof(IOSNPD));
      tree->slot[p].node = node;
      node->p = p;
      node->up = parent;
      node->level = (parent == NULL ? 0 : parent->level + 1);
      node->count = 0;
      node->b_ptr = NULL;
      node->s_ptr = NULL;
      node->r_ptr = NULL;
      node->solved = 0;
#if 0
      node->own_nn = node->own_nc = 0;
      node->e_ptr = NULL;
#endif
#if 1 /* 04/X-2008 */
      node->lp_obj = (parent == NULL ? (tree->mip->dir == GLP_MIN ?
         -DBL_MAX : +DBL_MAX) : parent->lp_obj);
#endif
      node->bound = (parent == NULL ? (tree->mip->dir == GLP_MIN ?
         -DBL_MAX : +DBL_MAX) : parent->bound);
      node->br_var = 0;
      node->br_val = 0.0;
      node->ii_cnt = 0;
      node->ii_sum = 0.0;
#if 1 /* 30/XI-2009 */
      node->changed = 0;
#endif
      if (tree->parm->cb_size == 0)
         node->data = NULL;
      else
      {  node->data = dmp_get_atom(tree->pool, tree->parm->cb_size);
         memset(node->data, 0, tree->parm->cb_size);
      }
      node->temp = NULL;
      node->prev = tree->tail;
      node->next = NULL;
      /* add the new subproblem to the end of the active list */
      if (tree->head == NULL)
         tree->head = node;
      else
         tree->tail->next = node;
      tree->tail = node;
      tree->a_cnt++;
      tree->n_cnt++;
      tree->t_cnt++;
      /* increase the number of child subproblems */
      if (parent == NULL)
         xassert(p == 1);
      else
         parent->count++;
      return node;
}

void ios_clone_node(glp_tree *tree, int p, int nnn, int ref[])
{     IOSNPD *node;
      int k;
      /* obtain pointer to the subproblem to be cloned */
      xassert(1 <= p && p <= tree->nslots);
      node = tree->slot[p].node;
      xassert(node != NULL);
      /* the specified subproblem must be active */
      xassert(node->count == 0);
      /* and must be in the frozen state */
      xassert(tree->curr != node);
      /* remove the specified subproblem from the active list, because
         it becomes inactive */
      if (node->prev == NULL)
         tree->head = node->next;
      else
         node->prev->next = node->next;
      if (node->next == NULL)
         tree->tail = node->prev;
      else
         node->next->prev = node->prev;
      node->prev = node->next = NULL;
      tree->a_cnt--;
      /* create clone subproblems */
      xassert(nnn > 0);
      for (k = 1; k <= nnn; k++)
         ref[k] = new_node(tree, node)->p;
      return;
}

/***********************************************************************
*  NAME
*
*  ios_delete_node - delete specified subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_delete_node(glp_tree *tree, int p);
*
*  DESCRIPTION
*
*  The routine ios_delete_node deletes the specified subproblem, whose
*  reference number is p. The subproblem must be active and must be in
*  the frozen state (i.e. it must not be the current subproblem).
*
*  Note that deletion is performed recursively, i.e. if a subproblem to
*  be deleted is the only child of its parent, the parent subproblem is
*  also deleted, etc. */

void ios_delete_node(glp_tree *tree, int p)
{     IOSNPD *node, *temp;
      /* obtain pointer to the subproblem to be deleted */
      xassert(1 <= p && p <= tree->nslots);
      node = tree->slot[p].node;
      xassert(node != NULL);
      /* the specified subproblem must be active */
      xassert(node->count == 0);
      /* and must be in the frozen state */
      xassert(tree->curr != node);
      /* remove the specified subproblem from the active list, because
         it is gone from the tree */
      if (node->prev == NULL)
         tree->head = node->next;
      else
         node->prev->next = node->next;
      if (node->next == NULL)
         tree->tail = node->prev;
      else
         node->next->prev = node->prev;
      node->prev = node->next = NULL;
      tree->a_cnt--;
loop: /* recursive deletion starts here */
      /* delete the bound change list */
      {  IOSBND *b;
         while (node->b_ptr != NULL)
         {  b = node->b_ptr;
            node->b_ptr = b->next;
            dmp_free_atom(tree->pool, b, sizeof(IOSBND));
         }
      }
      /* delete the status change list */
      {  IOSTAT *s;
         while (node->s_ptr != NULL)
         {  s = node->s_ptr;
            node->s_ptr = s->next;
            dmp_free_atom(tree->pool, s, sizeof(IOSTAT));
         }
      }
      /* delete the row addition list */
      while (node->r_ptr != NULL)
      {  IOSROW *r;
         r = node->r_ptr;
         if (r->name != NULL)
            dmp_free_atom(tree->pool, r->name, strlen(r->name)+1);
         while (r->ptr != NULL)
         {  IOSAIJ *a;
            a = r->ptr;
            r->ptr = a->next;
            dmp_free_atom(tree->pool, a, sizeof(IOSAIJ));
         }
         node->r_ptr = r->next;
         dmp_free_atom(tree->pool, r, sizeof(IOSROW));
      }
#if 0
      /* delete the edge addition list */
      /* delete the clique addition list */
      /* (not implemented yet) */
      xassert(node->own_nn == 0);
      xassert(node->own_nc == 0);
      xassert(node->e_ptr == NULL);
#endif
      /* free application-specific data */
      if (tree->parm->cb_size == 0)
         xassert(node->data == NULL);
      else
         dmp_free_atom(tree->pool, node->data, tree->parm->cb_size);
      /* free the corresponding node slot */
      p = node->p;
      xassert(tree->slot[p].node == node);
      tree->slot[p].node = NULL;
      tree->slot[p].next = tree->avail;
      tree->avail = p;
      /* save pointer to the parent subproblem */
      temp = node->up;
      /* delete the subproblem descriptor */
      dmp_free_atom(tree->pool, node, sizeof(IOSNPD));
      tree->n_cnt--;
      /* take pointer to the parent subproblem */
      node = temp;
      if (node != NULL)
      {  /* the parent subproblem exists; decrease the number of its
            child subproblems */
         xassert(node->count > 0);
         node->count--;
         /* if now the parent subproblem has no childs, it also must be
            deleted */
         if (node->count == 0) goto loop;
      }
      return;
}

/***********************************************************************
*  NAME
*
*  ios_delete_tree - delete branch-and-bound tree
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_delete_tree(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_delete_tree deletes the branch-and-bound tree, which
*  the parameter tree points to, and frees all the memory allocated to
*  this program object.
*
*  On exit components of the problem object are restored to correspond
*  to the original MIP passed to the routine ios_create_tree. */

void ios_delete_tree(glp_tree *tree)
{     glp_prob *mip = tree->mip;
      int i, j;
      int m = mip->m;
      int n = mip->n;
      xassert(mip->tree == tree);
      /* remove all additional rows */
      if (m != tree->orig_m)
      {  int nrs, *num;
         nrs = m - tree->orig_m;
         xassert(nrs > 0);
         num = xcalloc(1+nrs, sizeof(int));
         for (i = 1; i <= nrs; i++) num[i] = tree->orig_m + i;
         glp_del_rows(mip, nrs, num);
         xfree(num);
      }
      m = tree->orig_m;
      /* restore original attributes of rows and columns */
      xassert(m == tree->orig_m);
      xassert(n == tree->n);
      for (i = 1; i <= m; i++)
      {  glp_set_row_bnds(mip, i, tree->orig_type[i],
            tree->orig_lb[i], tree->orig_ub[i]);
         glp_set_row_stat(mip, i, tree->orig_stat[i]);
         mip->row[i]->prim = tree->orig_prim[i];
         mip->row[i]->dual = tree->orig_dual[i];
      }
      for (j = 1; j <= n; j++)
      {  glp_set_col_bnds(mip, j, tree->orig_type[m+j],
            tree->orig_lb[m+j], tree->orig_ub[m+j]);
         glp_set_col_stat(mip, j, tree->orig_stat[m+j]);
         mip->col[j]->prim = tree->orig_prim[m+j];
         mip->col[j]->dual = tree->orig_dual[m+j];
      }
      mip->pbs_stat = mip->dbs_stat = GLP_FEAS;
      mip->obj_val = tree->orig_obj;
      /* delete the branch-and-bound tree */
      xassert(tree->local != NULL);
      ios_delete_pool(tree, tree->local);
      dmp_delete_pool(tree->pool);
      xfree(tree->orig_type);
      xfree(tree->orig_lb);
      xfree(tree->orig_ub);
      xfree(tree->orig_stat);
      xfree(tree->orig_prim);
      xfree(tree->orig_dual);
      xfree(tree->slot);
      if (tree->root_type != NULL) xfree(tree->root_type);
      if (tree->root_lb != NULL) xfree(tree->root_lb);
      if (tree->root_ub != NULL) xfree(tree->root_ub);
      if (tree->root_stat != NULL) xfree(tree->root_stat);
      xfree(tree->non_int);
#if 0
      xfree(tree->n_ref);
      xfree(tree->c_ref);
      xfree(tree->j_ref);
#endif
      if (tree->pcost != NULL) ios_pcost_free(tree);
      xfree(tree->iwrk);
      xfree(tree->dwrk);
#if 0
      scg_delete_graph(tree->g);
#endif
      if (tree->pred_type != NULL) xfree(tree->pred_type);
      if (tree->pred_lb != NULL) xfree(tree->pred_lb);
      if (tree->pred_ub != NULL) xfree(tree->pred_ub);
      if (tree->pred_stat != NULL) xfree(tree->pred_stat);
#if 0
      xassert(tree->cut_gen == NULL);
#endif
      xassert(tree->mir_gen == NULL);
      xassert(tree->clq_gen == NULL);
      xfree(tree);
      mip->tree = NULL;
      return;
}

/***********************************************************************
*  NAME
*
*  ios_eval_degrad - estimate obj. degrad. for down- and up-branches
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  void ios_eval_degrad(glp_tree *tree, int j, double *dn, double *up);
*
*  DESCRIPTION
*
*  Given optimal basis to LP relaxation of the current subproblem the
*  routine ios_eval_degrad performs the dual ratio test to compute the
*  objective values in the adjacent basis for down- and up-branches,
*  which are stored in locations *dn and *up, assuming that x[j] is a
*  variable chosen to branch upon. */

void ios_eval_degrad(glp_tree *tree, int j, double *dn, double *up)
{     glp_prob *mip = tree->mip;
      int m = mip->m, n = mip->n;
      int len, kase, k, t, stat;
      double alfa, beta, gamma, delta, dz;
      int *ind = tree->iwrk;
      double *val = tree->dwrk;
      /* current basis must be optimal */
      xassert(glp_get_status(mip) == GLP_OPT);
      /* basis factorization must exist */
      xassert(glp_bf_exists(mip));
      /* obtain (fractional) value of x[j] in optimal basic solution
         to LP relaxation of the current subproblem */
      xassert(1 <= j && j <= n);
      beta = mip->col[j]->prim;
      /* since the value of x[j] is fractional, it is basic; compute
         corresponding row of the simplex table */
      len = lpx_eval_tab_row(mip, m+j, ind, val);
      /* kase < 0 means down-branch; kase > 0 means up-branch */
      for (kase = -1; kase <= +1; kase += 2)
      {  /* for down-branch we introduce new upper bound floor(beta)
            for x[j]; similarly, for up-branch we introduce new lower
            bound ceil(beta) for x[j]; in the current basis this new
            upper/lower bound is violated, so in the adjacent basis
            x[j] will leave the basis and go to its new upper/lower
            bound; we need to know which non-basic variable x[k] should
            enter the basis to keep dual feasibility */
#if 0 /* 23/XI-2009 */
         k = lpx_dual_ratio_test(mip, len, ind, val, kase, 1e-7);
#else
         k = lpx_dual_ratio_test(mip, len, ind, val, kase, 1e-9);
#endif
         /* if no variable has been chosen, current basis being primal
            infeasible due to the new upper/lower bound of x[j] is dual
            unbounded, therefore, LP relaxation to corresponding branch
            has no primal feasible solution */
         if (k == 0)
         {  if (mip->dir == GLP_MIN)
            {  if (kase < 0)
                  *dn = +DBL_MAX;
               else
                  *up = +DBL_MAX;
            }
            else if (mip->dir == GLP_MAX)
            {  if (kase < 0)
                  *dn = -DBL_MAX;
               else
                  *up = -DBL_MAX;
            }
            else
               xassert(mip != mip);
            continue;
         }
         xassert(1 <= k && k <= m+n);
         /* row of the simplex table corresponding to specified basic
            variable x[j] is the following:
               x[j] = ... + alfa * x[k] + ... ;
            we need to know influence coefficient, alfa, at non-basic
            variable x[k] chosen with the dual ratio test */
         for (t = 1; t <= len; t++)
            if (ind[t] == k) break;
         xassert(1 <= t && t <= len);
         alfa = val[t];
         /* determine status and reduced cost of variable x[k] */
         if (k <= m)
         {  stat = mip->row[k]->stat;
            gamma = mip->row[k]->dual;
         }
         else
         {  stat = mip->col[k-m]->stat;
            gamma = mip->col[k-m]->dual;
         }
         /* x[k] cannot be basic or fixed non-basic */
         xassert(stat == GLP_NL || stat == GLP_NU || stat == GLP_NF);
         /* if the current basis is dual degenerative, some reduced
            costs, which are close to zero, may have wrong sign due to
            round-off errors, so correct the sign of gamma */
         if (mip->dir == GLP_MIN)
         {  if (stat == GLP_NL && gamma < 0.0 ||
                stat == GLP_NU && gamma > 0.0 ||
                stat == GLP_NF) gamma = 0.0;
         }
         else if (mip->dir == GLP_MAX)
         {  if (stat == GLP_NL && gamma > 0.0 ||
                stat == GLP_NU && gamma < 0.0 ||
                stat == GLP_NF) gamma = 0.0;
         }
         else
            xassert(mip != mip);
         /* determine the change of x[j] in the adjacent basis:
            delta x[j] = new x[j] - old x[j] */
         delta = (kase < 0 ? floor(beta) : ceil(beta)) - beta;
         /* compute the change of x[k] in the adjacent basis:
            delta x[k] = new x[k] - old x[k] = delta x[j] / alfa */
         delta /= alfa;
         /* compute the change of the objective in the adjacent basis:
            delta z = new z - old z = gamma * delta x[k] */
         dz = gamma * delta;
         if (mip->dir == GLP_MIN)
            xassert(dz >= 0.0);
         else if (mip->dir == GLP_MAX)
            xassert(dz <= 0.0);
         else
            xassert(mip != mip);
         /* compute the new objective value in the adjacent basis:
            new z = old z + delta z */
         if (kase < 0)
            *dn = mip->obj_val + dz;
         else
            *up = mip->obj_val + dz;
      }
      /*xprintf("obj = %g; dn = %g; up = %g\n",
         mip->obj_val, *dn, *up);*/
      return;
}

/***********************************************************************
*  NAME
*
*  ios_round_bound - improve local bound by rounding
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  double ios_round_bound(glp_tree *tree, double bound);
*
*  RETURNS
*
*  For the given local bound for any integer feasible solution to the
*  current subproblem the routine ios_round_bound returns an improved
*  local bound for the same integer feasible solution.
*
*  BACKGROUND
*
*  Let the current subproblem has the following objective function:
*
*     z =   sum  c[j] * x[j] + s >= b,                               (1)
*         j in J
*
*  where J = {j: c[j] is non-zero and integer, x[j] is integer}, s is
*  the sum of terms corresponding to fixed variables, b is an initial
*  local bound (minimization).
*
*  From (1) it follows that:
*
*     d *  sum  (c[j] / d) * x[j] + s >= b,                          (2)
*        j in J
*
*  or, equivalently,
*
*     sum  (c[j] / d) * x[j] >= (b - s) / d = h,                     (3)
*   j in J
*
*  where d = gcd(c[j]). Since the left-hand side of (3) is integer,
*  h = (b - s) / d can be rounded up to the nearest integer:
*
*     h' = ceil(h) = (b' - s) / d,                                   (4)
*
*  that gives an rounded, improved local bound:
*
*     b' = d * h' + s.                                               (5)
*
*  In case of maximization '>=' in (1) should be replaced by '<=' that
*  leads to the following formula:
*
*     h' = floor(h) = (b' - s) / d,                                  (6)
*
*  which should used in the same way as (4).
*
*  NOTE: If b is a valid local bound for a child of the current
*        subproblem, b' is also valid for that child subproblem. */

double ios_round_bound(glp_tree *tree, double bound)
{     glp_prob *mip = tree->mip;
      int n = mip->n;
      int d, j, nn, *c = tree->iwrk;
      double s, h;
      /* determine c[j] and compute s */
      nn = 0, s = mip->c0, d = 0;
      for (j = 1; j <= n; j++)
      {  GLPCOL *col = mip->col[j];
         if (col->coef == 0.0) continue;
         if (col->type == GLP_FX)
         {  /* fixed variable */
            s += col->coef * col->prim;
         }
         else
         {  /* non-fixed variable */
            if (col->kind != GLP_IV) goto skip;
            if (col->coef != floor(col->coef)) goto skip;
            if (fabs(col->coef) <= (double)INT_MAX)
               c[++nn] = (int)fabs(col->coef);
            else
               d = 1;
         }
      }
      /* compute d = gcd(c[1],...c[nn]) */
      if (d == 0)
      {  if (nn == 0) goto skip;
         d = gcdn(nn, c);
      }
      xassert(d > 0);
      /* compute new local bound */
      if (mip->dir == GLP_MIN)
      {  if (bound != +DBL_MAX)
         {  h = (bound - s) / (double)d;
            if (h >= floor(h) + 0.001)
            {  /* round up */
               h = ceil(h);
               /*xprintf("d = %d; old = %g; ", d, bound);*/
               bound = (double)d * h + s;
               /*xprintf("new = %g\n", bound);*/
            }
         }
      }
      else if (mip->dir == GLP_MAX)
      {  if (bound != -DBL_MAX)
         {  h = (bound - s) / (double)d;
            if (h <= ceil(h) - 0.001)
            {  /* round down */
               h = floor(h);
               bound = (double)d * h + s;
            }
         }
      }
      else
         xassert(mip != mip);
skip: return bound;
}

/***********************************************************************
*  NAME
*
*  ios_is_hopeful - check if subproblem is hopeful
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_is_hopeful(glp_tree *tree, double bound);
*
*  DESCRIPTION
*
*  Given the local bound of a subproblem the routine ios_is_hopeful
*  checks if the subproblem can have an integer optimal solution which
*  is better than the best one currently known.
*
*  RETURNS
*
*  If the subproblem can have a better integer optimal solution, the
*  routine returns non-zero; otherwise, if the corresponding branch can
*  be pruned, the routine returns zero. */

int ios_is_hopeful(glp_tree *tree, double bound)
{     glp_prob *mip = tree->mip;
      int ret = 1;
      double eps;
      if (mip->mip_stat == GLP_FEAS)
      {  eps = tree->parm->tol_obj * (1.0 + fabs(mip->mip_obj));
         switch (mip->dir)
         {  case GLP_MIN:
               if (bound >= mip->mip_obj - eps) ret = 0;
               break;
            case GLP_MAX:
               if (bound <= mip->mip_obj + eps) ret = 0;
               break;
            default:
               xassert(mip != mip);
         }
      }
      else
      {  switch (mip->dir)
         {  case GLP_MIN:
               if (bound == +DBL_MAX) ret = 0;
               break;
            case GLP_MAX:
               if (bound == -DBL_MAX) ret = 0;
               break;
            default:
               xassert(mip != mip);
         }
      }
      return ret;
}

/***********************************************************************
*  NAME
*
*  ios_best_node - find active node with best local bound
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_best_node(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_best_node finds an active node whose local bound is
*  best among other active nodes.
*
*  It is understood that the integer optimal solution of the original
*  mip problem cannot be better than the best bound, so the best bound
*  is an lower (minimization) or upper (maximization) global bound for
*  the original problem.
*
*  RETURNS
*
*  The routine ios_best_node returns the subproblem reference number
*  for the best node. However, if the tree is empty, it returns zero. */

int ios_best_node(glp_tree *tree)
{     IOSNPD *node, *best = NULL;
      switch (tree->mip->dir)
      {  case GLP_MIN:
            /* minimization */
            for (node = tree->head; node != NULL; node = node->next)
               if (best == NULL || best->bound > node->bound)
                  best = node;
            break;
         case GLP_MAX:
            /* maximization */
            for (node = tree->head; node != NULL; node = node->next)
               if (best == NULL || best->bound < node->bound)
                  best = node;
            break;
         default:
            xassert(tree != tree);
      }
      return best == NULL ? 0 : best->p;
}

/***********************************************************************
*  NAME
*
*  ios_relative_gap - compute relative mip gap
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  double ios_relative_gap(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_relative_gap computes the relative mip gap using the
*  formula:
*
*     gap = |best_mip - best_bnd| / (|best_mip| + DBL_EPSILON),
*
*  where best_mip is the best integer feasible solution found so far,
*  best_bnd is the best (global) bound. If no integer feasible solution
*  has been found yet, rel_gap is set to DBL_MAX.
*
*  RETURNS
*
*  The routine ios_relative_gap returns the relative mip gap. */

double ios_relative_gap(glp_tree *tree)
{     glp_prob *mip = tree->mip;
      int p;
      double best_mip, best_bnd, gap;
      if (mip->mip_stat == GLP_FEAS)
      {  best_mip = mip->mip_obj;
         p = ios_best_node(tree);
         if (p == 0)
         {  /* the tree is empty */
            gap = 0.0;
         }
         else
         {  best_bnd = tree->slot[p].node->bound;
            gap = fabs(best_mip - best_bnd) / (fabs(best_mip) +
               DBL_EPSILON);
         }
      }
      else
      {  /* no integer feasible solution has been found yet */
         gap = DBL_MAX;
      }
      return gap;
}

/***********************************************************************
*  NAME
*
*  ios_solve_node - solve LP relaxation of current subproblem
*
*  SYNOPSIS
*
*  #include "glpios.h"
*  int ios_solve_node(glp_tree *tree);
*
*  DESCRIPTION
*
*  The routine ios_solve_node re-optimizes LP relaxation of the current
*  subproblem using the dual simplex method.
*
*  RETURNS
*
*  The routine returns the code which is reported by glp_simplex. */

int ios_solve_node(glp_tree *tree)
{     glp_prob *mip = tree->mip;
      glp_smcp parm;
      int ret;
      /* the current subproblem must exist */
      xassert(tree->curr != NULL);
      /* set some control parameters */
      glp_init_smcp(&parm);
      switch (tree->parm->msg_lev)
      {  case GLP_MSG_OFF:
            parm.msg_lev = GLP_MSG_OFF; break;
         case GLP_MSG_ERR:
            parm.msg_lev = GLP_MSG_ERR; break;
         case GLP_MSG_ON:
         case GLP_MSG_ALL:
            parm.msg_lev = GLP_MSG_ON; break;
         case GLP_MSG_DBG:
            parm.msg_lev = GLP_MSG_ALL; break;
         default:
            xassert(tree != tree);
      }
      parm.meth = GLP_DUALP;
      if (tree->parm->msg_lev < GLP_MSG_DBG)
         parm.out_dly = tree->parm->out_dly;
      else
         parm.out_dly = 0;
      /* if the incumbent objective value is already known, use it to
         prematurely terminate the dual simplex search */
      if (mip->mip_stat == GLP_FEAS)
      {  switch (tree->mip->dir)
         {  case GLP_MIN:
               parm.obj_ul = mip->mip_obj;
               break;
            case GLP_MAX:
               parm.obj_ll = mip->mip_obj;
               break;
            default:
               xassert(mip != mip);
         }
      }
      /* try to solve/re-optimize the LP relaxation */
      ret = glp_simplex(mip, &parm);
      tree->curr->solved++;
#if 0
      xprintf("ret = %d; status = %d; pbs = %d; dbs = %d; some = %d\n",
         ret, glp_get_status(mip), mip->pbs_stat, mip->dbs_stat,
         mip->some);
      lpx_print_sol(mip, "sol");
#endif
      return ret;
}

/**********************************************************************/

IOSPOOL *ios_create_pool(glp_tree *tree)
{     /* create cut pool */
      IOSPOOL *pool;
#if 0
      pool = dmp_get_atom(tree->pool, sizeof(IOSPOOL));
#else
      xassert(tree == tree);
      pool = xmalloc(sizeof(IOSPOOL));
#endif
      pool->size = 0;
      pool->head = pool->tail = NULL;
      pool->ord = 0, pool->curr = NULL;
      return pool;
}

int ios_add_row(glp_tree *tree, IOSPOOL *pool,
      const char *name, int klass, int flags, int len, const int ind[],
      const double val[], int type, double rhs)
{     /* add row (constraint) to the cut pool */
      IOSCUT *cut;
      IOSAIJ *aij;
      int k;
      xassert(pool != NULL);
      cut = dmp_get_atom(tree->pool, sizeof(IOSCUT));
      if (name == NULL || name[0] == '\0')
         cut->name = NULL;
      else
      {  for (k = 0; name[k] != '\0'; k++)
         {  if (k == 256)
               xerror("glp_ios_add_row: cut name too long\n");
            if (iscntrl((unsigned char)name[k]))
               xerror("glp_ios_add_row: cut name contains invalid chara"
                  "cter(s)\n");
         }
         cut->name = dmp_get_atom(tree->pool, strlen(name)+1);
         strcpy(cut->name, name);
      }
      if (!(0 <= klass && klass <= 255))
         xerror("glp_ios_add_row: klass = %d; invalid cut class\n",
            klass);
      cut->klass = (unsigned char)klass;
      if (flags != 0)
         xerror("glp_ios_add_row: flags = %d; invalid cut flags\n",
            flags);
      cut->ptr = NULL;
      if (!(0 <= len && len <= tree->n))
         xerror("glp_ios_add_row: len = %d; invalid cut length\n",
            len);
      for (k = 1; k <= len; k++)
      {  aij = dmp_get_atom(tree->pool, sizeof(IOSAIJ));
         if (!(1 <= ind[k] && ind[k] <= tree->n))
            xerror("glp_ios_add_row: ind[%d] = %d; column index out of "
               "range\n", k, ind[k]);
         aij->j = ind[k];
         aij->val = val[k];
         aij->next = cut->ptr;
         cut->ptr = aij;
      }
      if (!(type == GLP_LO || type == GLP_UP || type == GLP_FX))
         xerror("glp_ios_add_row: type = %d; invalid cut type\n",
            type);
      cut->type = (unsigned char)type;
      cut->rhs = rhs;
      cut->prev = pool->tail;
      cut->next = NULL;
      if (cut->prev == NULL)
         pool->head = cut;
      else
         cut->prev->next = cut;
      pool->tail = cut;
      pool->size++;
      return pool->size;
}

IOSCUT *ios_find_row(IOSPOOL *pool, int i)
{     /* find row (constraint) in the cut pool */
      /* (smart linear search) */
      xassert(pool != NULL);
      xassert(1 <= i && i <= pool->size);
      if (pool->ord == 0)
      {  xassert(pool->curr == NULL);
         pool->ord = 1;
         pool->curr = pool->head;
      }
      xassert(pool->curr != NULL);
      if (i < pool->ord)
      {  if (i < pool->ord - i)
         {  pool->ord = 1;
            pool->curr = pool->head;
            while (pool->ord != i)
            {  pool->ord++;
               xassert(pool->curr != NULL);
               pool->curr = pool->curr->next;
            }
         }
         else
         {  while (pool->ord != i)
            {  pool->ord--;
               xassert(pool->curr != NULL);
               pool->curr = pool->curr->prev;
            }
         }
      }
      else if (i > pool->ord)
      {  if (i - pool->ord < pool->size - i)
         {  while (pool->ord != i)
            {  pool->ord++;
               xassert(pool->curr != NULL);
               pool->curr = pool->curr->next;
            }
         }
         else
         {  pool->ord = pool->size;
            pool->curr = pool->tail;
            while (pool->ord != i)
            {  pool->ord--;
               xassert(pool->curr != NULL);
               pool->curr = pool->curr->prev;
            }
         }
      }
      xassert(pool->ord == i);
      xassert(pool->curr != NULL);
      return pool->curr;
}

void ios_del_row(glp_tree *tree, IOSPOOL *pool, int i)
{     /* remove row (constraint) from the cut pool */
      IOSCUT *cut;
      IOSAIJ *aij;
      xassert(pool != NULL);
      if (!(1 <= i && i <= pool->size))
         xerror("glp_ios_del_row: i = %d; cut number out of range\n",
            i);
      cut = ios_find_row(pool, i);
      xassert(pool->curr == cut);
      if (cut->next != NULL)
         pool->curr = cut->next;
      else if (cut->prev != NULL)
         pool->ord--, pool->curr = cut->prev;
      else
         pool->ord = 0, pool->curr = NULL;
      if (cut->name != NULL)
         dmp_free_atom(tree->pool, cut->name, strlen(cut->name)+1);
      if (cut->prev == NULL)
      {  xassert(pool->head == cut);
         pool->head = cut->next;
      }
      else
      {  xassert(cut->prev->next == cut);
         cut->prev->next = cut->next;
      }
      if (cut->next == NULL)
      {  xassert(pool->tail == cut);
         pool->tail = cut->prev;
      }
      else
      {  xassert(cut->next->prev == cut);
         cut->next->prev = cut->prev;
      }
      while (cut->ptr != NULL)
      {  aij = cut->ptr;
         cut->ptr = aij->next;
         dmp_free_atom(tree->pool, aij, sizeof(IOSAIJ));
      }
      dmp_free_atom(tree->pool, cut, sizeof(IOSCUT));
      pool->size--;
      return;
}

void ios_clear_pool(glp_tree *tree, IOSPOOL *pool)
{     /* remove all rows (constraints) from the cut pool */
      xassert(pool != NULL);
      while (pool->head != NULL)
      {  IOSCUT *cut = pool->head;
         pool->head = cut->next;
         if (cut->name != NULL)
            dmp_free_atom(tree->pool, cut->name, strlen(cut->name)+1);
         while (cut->ptr != NULL)
         {  IOSAIJ *aij = cut->ptr;
            cut->ptr = aij->next;
            dmp_free_atom(tree->pool, aij, sizeof(IOSAIJ));
         }
         dmp_free_atom(tree->pool, cut, sizeof(IOSCUT));
      }
      pool->size = 0;
      pool->head = pool->tail = NULL;
      pool->ord = 0, pool->curr = NULL;
      return;
}

void ios_delete_pool(glp_tree *tree, IOSPOOL *pool)
{     /* delete cut pool */
      xassert(pool != NULL);
      ios_clear_pool(tree, pool);
      xfree(pool);
      return;
}

/**********************************************************************/

#if 0
static int refer_to_node(glp_tree *tree, int j)
{     /* determine node number corresponding to binary variable x[j] or
         its complement */
      glp_prob *mip = tree->mip;
      int n = mip->n;
      int *ref;
      if (j > 0)
         ref = tree->n_ref;
      else
         ref = tree->c_ref, j = - j;
      xassert(1 <= j && j <= n);
      if (ref[j] == 0)
      {  /* new node is needed */
         SCG *g = tree->g;
         int n_max = g->n_max;
         ref[j] = scg_add_nodes(g, 1);
         if (g->n_max > n_max)
         {  int *save = tree->j_ref;
            tree->j_ref = xcalloc(1+g->n_max, sizeof(int));
            memcpy(&tree->j_ref[1], &save[1], g->n * sizeof(int));
            xfree(save);
         }
         xassert(ref[j] == g->n);
         tree->j_ref[ref[j]] = j;
         xassert(tree->curr != NULL);
         if (tree->curr->level > 0) tree->curr->own_nn++;
      }
      return ref[j];
}
#endif

#if 0
void ios_add_edge(glp_tree *tree, int j1, int j2)
{     /* add new edge to the conflict graph */
      glp_prob *mip = tree->mip;
      int n = mip->n;
      SCGRIB *e;
      int first, i1, i2;
      xassert(-n <= j1 && j1 <= +n && j1 != 0);
      xassert(-n <= j2 && j2 <= +n && j2 != 0);
      xassert(j1 != j2);
      /* determine number of the first node, which was added for the
         current subproblem */
      xassert(tree->curr != NULL);
      first = tree->g->n - tree->curr->own_nn + 1;
      /* determine node numbers for both endpoints */
      i1 = refer_to_node(tree, j1);
      i2 = refer_to_node(tree, j2);
      /* add edge (i1,i2) to the conflict graph */
      e = scg_add_edge(tree->g, i1, i2);
      /* if the current subproblem is not the root and both endpoints
         were created on some previous levels, save the edge */
      if (tree->curr->level > 0 && i1 < first && i2 < first)
      {  IOSRIB *rib;
         rib = dmp_get_atom(tree->pool, sizeof(IOSRIB));
         rib->j1 = j1;
         rib->j2 = j2;
         rib->e = e;
         rib->next = tree->curr->e_ptr;
         tree->curr->e_ptr = rib;
      }
      return;
}
#endif

/* eof */

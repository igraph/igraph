/* mpl2.c */

/***********************************************************************
*  This code is part of GLPK (GNU Linear Programming Kit).
*  Copyright (C) 2003-2016 Free Software Foundation, Inc.
*  Written by Andrew Makhorin <mao@gnu.org>.
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

#include "mpl.h"

/**********************************************************************/
/* * *                  PROCESSING DATA SECTION                   * * */
/**********************************************************************/

/*----------------------------------------------------------------------
-- create_slice - create slice.
--
-- This routine creates a slice, which initially has no components. */

SLICE *create_slice(MPL *mpl)
{     SLICE *slice;
      xassert(mpl == mpl);
      slice = NULL;
      return slice;
}

/*----------------------------------------------------------------------
-- expand_slice - append new component to slice.
--
-- This routine expands slice appending to it either a given symbol or
-- null component, which becomes the last component of the slice. */

SLICE *expand_slice
(     MPL *mpl,
      SLICE *slice,           /* destroyed */
      SYMBOL *sym             /* destroyed */
)
{     SLICE *tail, *temp;
      /* create a new component */
      tail = dmp_get_atom(mpl->tuples, sizeof(SLICE));
      tail->sym = sym;
      tail->next = NULL;
      /* and append it to the component list */
      if (slice == NULL)
         slice = tail;
      else
      {  for (temp = slice; temp->next != NULL; temp = temp->next);
         temp->next = tail;
      }
      return slice;
}

/*----------------------------------------------------------------------
-- slice_dimen - determine dimension of slice.
--
-- This routine returns dimension of slice, which is number of all its
-- components including null ones. */

int slice_dimen
(     MPL *mpl,
      SLICE *slice            /* not changed */
)
{     SLICE *temp;
      int dim;
      xassert(mpl == mpl);
      dim = 0;
      for (temp = slice; temp != NULL; temp = temp->next) dim++;
      return dim;
}

/*----------------------------------------------------------------------
-- slice_arity - determine arity of slice.
--
-- This routine returns arity of slice, i.e. number of null components
-- (indicated by asterisks) in the slice. */

int slice_arity
(     MPL *mpl,
      SLICE *slice            /* not changed */
)
{     SLICE *temp;
      int arity;
      xassert(mpl == mpl);
      arity = 0;
      for (temp = slice; temp != NULL; temp = temp->next)
         if (temp->sym == NULL) arity++;
      return arity;
}

/*----------------------------------------------------------------------
-- fake_slice - create fake slice of all asterisks.
--
-- This routine creates a fake slice of given dimension, which contains
-- asterisks in all components. Zero dimension is allowed. */

SLICE *fake_slice(MPL *mpl, int dim)
{     SLICE *slice;
      slice = create_slice(mpl);
      while (dim-- > 0) slice = expand_slice(mpl, slice, NULL);
      return slice;
}

/*----------------------------------------------------------------------
-- delete_slice - delete slice.
--
-- This routine deletes specified slice. */

void delete_slice
(     MPL *mpl,
      SLICE *slice            /* destroyed */
)
{     SLICE *temp;
      while (slice != NULL)
      {  temp = slice;
         slice = temp->next;
         if (temp->sym != NULL) delete_symbol(mpl, temp->sym);
xassert(sizeof(SLICE) == sizeof(TUPLE));
         dmp_free_atom(mpl->tuples, temp, sizeof(TUPLE));
      }
      return;
}

/*----------------------------------------------------------------------
-- is_number - check if current token is number.
--
-- If the current token is a number, this routine returns non-zero.
-- Otherwise zero is returned. */

int is_number(MPL *mpl)
{     return
         mpl->token == T_NUMBER;
}

/*----------------------------------------------------------------------
-- is_symbol - check if current token is symbol.
--
-- If the current token is suitable to be a symbol, the routine returns
-- non-zero. Otherwise zero is returned. */

int is_symbol(MPL *mpl)
{     return
         mpl->token == T_NUMBER ||
         mpl->token == T_SYMBOL ||
         mpl->token == T_STRING;
}

/*----------------------------------------------------------------------
-- is_literal - check if current token is given symbolic literal.
--
-- If the current token is given symbolic literal, this routine returns
-- non-zero. Otherwise zero is returned.
--
-- This routine is used on processing the data section in the same way
-- as the routine is_keyword on processing the model section. */

int is_literal(MPL *mpl, char *literal)
{     return
         is_symbol(mpl) && strcmp(mpl->image, literal) == 0;
}

/*----------------------------------------------------------------------
-- read_number - read number.
--
-- This routine reads the current token, which must be a number, and
-- returns its numeric value. */

double read_number(MPL *mpl)
{     double num;
      xassert(is_number(mpl));
      num = mpl->value;
      get_token(mpl /* <number> */);
      return num;
}

/*----------------------------------------------------------------------
-- read_symbol - read symbol.
--
-- This routine reads the current token, which must be a symbol, and
-- returns its symbolic value. */

SYMBOL *read_symbol(MPL *mpl)
{     SYMBOL *sym;
      xassert(is_symbol(mpl));
      if (is_number(mpl))
         sym = create_symbol_num(mpl, mpl->value);
      else
         sym = create_symbol_str(mpl, create_string(mpl, mpl->image));
      get_token(mpl /* <symbol> */);
      return sym;
}

/*----------------------------------------------------------------------
-- read_slice - read slice.
--
-- This routine reads slice using the syntax:
--
-- <slice> ::= [ <symbol list> ]
-- <slice> ::= ( <symbol list> )
-- <symbol list> ::= <symbol or star>
-- <symbol list> ::= <symbol list> , <symbol or star>
-- <symbol or star> ::= <symbol>
-- <symbol or star> ::= *
--
-- The bracketed form of slice is used for members of multi-dimensional
-- objects while the parenthesized form is used for elemental sets. */

SLICE *read_slice
(     MPL *mpl,
      char *name,             /* not changed */
      int dim
)
{     SLICE *slice;
      int close;
      xassert(name != NULL);
      switch (mpl->token)
      {  case T_LBRACKET:
            close = T_RBRACKET;
            break;
         case T_LEFT:
            xassert(dim > 0);
            close = T_RIGHT;
            break;
         default:
            xassert(mpl != mpl);
      }
      if (dim == 0)
         error(mpl, "%s cannot be subscripted", name);
      get_token(mpl /* ( | [ */);
      /* read slice components */
      slice = create_slice(mpl);
      for (;;)
      {  /* the current token must be a symbol or asterisk */
         if (is_symbol(mpl))
            slice = expand_slice(mpl, slice, read_symbol(mpl));
         else if (mpl->token == T_ASTERISK)
         {  slice = expand_slice(mpl, slice, NULL);
            get_token(mpl /* * */);
         }
         else
            error(mpl, "number, symbol, or asterisk missing where expec"
               "ted");
         /* check a token that follows the symbol */
         if (mpl->token == T_COMMA)
            get_token(mpl /* , */);
         else if (mpl->token == close)
            break;
         else
            error(mpl, "syntax error in slice");
      }
      /* number of slice components must be the same as the appropriate
         dimension */
      if (slice_dimen(mpl, slice) != dim)
      {  switch (close)
         {  case T_RBRACKET:
               error(mpl, "%s must have %d subscript%s, not %d", name,
                  dim, dim == 1 ? "" : "s", slice_dimen(mpl, slice));
               break;
            case T_RIGHT:
               error(mpl, "%s has dimension %d, not %d", name, dim,
                  slice_dimen(mpl, slice));
               break;
            default:
               xassert(close != close);
         }
      }
      get_token(mpl /* ) | ] */);
      return slice;
}

/*----------------------------------------------------------------------
-- select_set - select set to saturate it with elemental sets.
--
-- This routine selects set to saturate it with elemental sets provided
-- in the data section. */

SET *select_set
(     MPL *mpl,
      char *name              /* not changed */
)
{     SET *set;
      AVLNODE *node;
      xassert(name != NULL);
      node = avl_find_node(mpl->tree, name);
      if (node == NULL || avl_get_node_type(node) != A_SET)
         error(mpl, "%s not a set", name);
      set = (SET *)avl_get_node_link(node);
      if (set->assign != NULL || set->gadget != NULL)
         error(mpl, "%s needs no data", name);
      set->data = 1;
      return set;
}

/*----------------------------------------------------------------------
-- simple_format - read set data block in simple format.
--
-- This routine reads set data block using the syntax:
--
-- <simple format> ::= <symbol> , <symbol> , ... , <symbol>
--
-- where <symbols> are used to construct a complete n-tuple, which is
-- included in elemental set assigned to the set member. Commae between
-- symbols are optional and may be omitted anywhere.
--
-- Number of components in the slice must be the same as dimension of
-- n-tuples in elemental sets assigned to the set members. To construct
-- complete n-tuple the routine replaces null positions in the slice by
-- corresponding <symbols>.
--
-- If the slice contains at least one null position, the current token
-- must be symbol. Otherwise, the routine reads no symbols to construct
-- the n-tuple, so the current token is not checked. */

void simple_format
(     MPL *mpl,
      SET *set,               /* not changed */
      MEMBER *memb,           /* modified */
      SLICE *slice            /* not changed */
)
{     TUPLE *tuple;
      SLICE *temp;
      SYMBOL *sym, *with = NULL;
      xassert(set != NULL);
      xassert(memb != NULL);
      xassert(slice != NULL);
      xassert(set->dimen == slice_dimen(mpl, slice));
      xassert(memb->value.set->dim == set->dimen);
      if (slice_arity(mpl, slice) > 0) xassert(is_symbol(mpl));
      /* read symbols and construct complete n-tuple */
      tuple = create_tuple(mpl);
      for (temp = slice; temp != NULL; temp = temp->next)
      {  if (temp->sym == NULL)
         {  /* substitution is needed; read symbol */
            if (!is_symbol(mpl))
            {  int lack = slice_arity(mpl, temp);
               /* with cannot be null due to assertion above */
               xassert(with != NULL);
               if (lack == 1)
                  error(mpl, "one item missing in data group beginning "
                     "with %s", format_symbol(mpl, with));
               else
                  error(mpl, "%d items missing in data group beginning "
                     "with %s", lack, format_symbol(mpl, with));
            }
            sym = read_symbol(mpl);
            if (with == NULL) with = sym;
         }
         else
         {  /* copy symbol from the slice */
            sym = copy_symbol(mpl, temp->sym);
         }
         /* append the symbol to the n-tuple */
         tuple = expand_tuple(mpl, tuple, sym);
         /* skip optional comma *between* <symbols> */
         if (temp->next != NULL && mpl->token == T_COMMA)
            get_token(mpl /* , */);
      }
      /* add constructed n-tuple to elemental set */
      check_then_add(mpl, memb->value.set, tuple);
      return;
}

/*----------------------------------------------------------------------
-- matrix_format - read set data block in matrix format.
--
-- This routine reads set data block using the syntax:
--
-- <matrix format> ::= <column> <column> ... <column> :=
--               <row>   +/-      +/-    ...   +/-
--               <row>   +/-      +/-    ...   +/-
--                 .  .  .  .  .  .  .  .  .  .  .
--               <row>   +/-      +/-    ...   +/-
--
-- where <rows> are symbols that denote rows of the matrix, <columns>
-- are symbols that denote columns of the matrix, "+" and "-" indicate
-- whether corresponding n-tuple needs to be included in the elemental
-- set or not, respectively.
--
-- Number of the slice components must be the same as dimension of the
-- elemental set. The slice must have two null positions. To construct
-- complete n-tuple for particular element of the matrix the routine
-- replaces first null position of the slice by the corresponding <row>
-- (or <column>, if the flag tr is on) and second null position by the
-- corresponding <column> (or by <row>, if the flag tr is on). */

void matrix_format
(     MPL *mpl,
      SET *set,               /* not changed */
      MEMBER *memb,           /* modified */
      SLICE *slice,           /* not changed */
      int tr
)
{     SLICE *list, *col, *temp;
      TUPLE *tuple;
      SYMBOL *row;
      xassert(set != NULL);
      xassert(memb != NULL);
      xassert(slice != NULL);
      xassert(set->dimen == slice_dimen(mpl, slice));
      xassert(memb->value.set->dim == set->dimen);
      xassert(slice_arity(mpl, slice) == 2);
      /* read the matrix heading that contains column symbols (there
         may be no columns at all) */
      list = create_slice(mpl);
      while (mpl->token != T_ASSIGN)
      {  /* read column symbol and append it to the column list */
         if (!is_symbol(mpl))
            error(mpl, "number, symbol, or := missing where expected");
         list = expand_slice(mpl, list, read_symbol(mpl));
      }
      get_token(mpl /* := */);
      /* read zero or more rows that contain matrix data */
      while (is_symbol(mpl))
      {  /* read row symbol (if the matrix has no columns, row symbols
            are just ignored) */
         row = read_symbol(mpl);
         /* read the matrix row accordingly to the column list */
         for (col = list; col != NULL; col = col->next)
         {  int which = 0;
            /* check indicator */
            if (is_literal(mpl, "+"))
               ;
            else if (is_literal(mpl, "-"))
            {  get_token(mpl /* - */);
               continue;
            }
            else
            {  int lack = slice_dimen(mpl, col);
               if (lack == 1)
                  error(mpl, "one item missing in data group beginning "
                     "with %s", format_symbol(mpl, row));
               else
                  error(mpl, "%d items missing in data group beginning "
                     "with %s", lack, format_symbol(mpl, row));
            }
            /* construct complete n-tuple */
            tuple = create_tuple(mpl);
            for (temp = slice; temp != NULL; temp = temp->next)
            {  if (temp->sym == NULL)
               {  /* substitution is needed */
                  switch (++which)
                  {  case 1:
                        /* substitute in the first null position */
                        tuple = expand_tuple(mpl, tuple,
                           copy_symbol(mpl, tr ? col->sym : row));
                        break;
                     case 2:
                        /* substitute in the second null position */
                        tuple = expand_tuple(mpl, tuple,
                           copy_symbol(mpl, tr ? row : col->sym));
                        break;
                     default:
                        xassert(which != which);
                  }
               }
               else
               {  /* copy symbol from the slice */
                  tuple = expand_tuple(mpl, tuple, copy_symbol(mpl,
                     temp->sym));
               }
            }
            xassert(which == 2);
            /* add constructed n-tuple to elemental set */
            check_then_add(mpl, memb->value.set, tuple);
            get_token(mpl /* + */);
         }
         /* delete the row symbol */
         delete_symbol(mpl, row);
      }
      /* delete the column list */
      delete_slice(mpl, list);
      return;
}

/*----------------------------------------------------------------------
-- set_data - read set data.
--
-- This routine reads set data using the syntax:
--
-- <set data> ::= set <set name> <assignments> ;
-- <set data> ::= set <set name> [ <symbol list> ] <assignments> ;
-- <set name> ::= <symbolic name>
-- <assignments> ::= <empty>
-- <assignments> ::= <assignments> , :=
-- <assignments> ::= <assignments> , ( <symbol list> )
-- <assignments> ::= <assignments> , <simple format>
-- <assignments> ::= <assignments> , : <matrix format>
-- <assignments> ::= <assignments> , (tr) <matrix format>
-- <assignments> ::= <assignments> , (tr) : <matrix format>
--
-- Commae in <assignments> are optional and may be omitted anywhere. */

void set_data(MPL *mpl)
{     SET *set;
      TUPLE *tuple;
      MEMBER *memb;
      SLICE *slice;
      int tr = 0;
      xassert(is_literal(mpl, "set"));
      get_token(mpl /* set */);
      /* symbolic name of set must follows the keyword 'set' */
      if (!is_symbol(mpl))
         error(mpl, "set name missing where expected");
      /* select the set to saturate it with data */
      set = select_set(mpl, mpl->image);
      get_token(mpl /* <symbolic name> */);
      /* read optional subscript list, which identifies member of the
         set to be read */
      tuple = create_tuple(mpl);
      if (mpl->token == T_LBRACKET)
      {  /* subscript list is specified */
         if (set->dim == 0)
            error(mpl, "%s cannot be subscripted", set->name);
         get_token(mpl /* [ */);
         /* read symbols and construct subscript list */
         for (;;)
         {  if (!is_symbol(mpl))
               error(mpl, "number or symbol missing where expected");
            tuple = expand_tuple(mpl, tuple, read_symbol(mpl));
            if (mpl->token == T_COMMA)
               get_token(mpl /* , */);
            else if (mpl->token == T_RBRACKET)
               break;
            else
               error(mpl, "syntax error in subscript list");
         }
         if (set->dim != tuple_dimen(mpl, tuple))
            error(mpl, "%s must have %d subscript%s rather than %d",
               set->name, set->dim, set->dim == 1 ? "" : "s",
               tuple_dimen(mpl, tuple));
         get_token(mpl /* ] */);
      }
      else
      {  /* subscript list is not specified */
         if (set->dim != 0)
            error(mpl, "%s must be subscripted", set->name);
      }
      /* there must be no member with the same subscript list */
      if (find_member(mpl, set->array, tuple) != NULL)
         error(mpl, "%s%s already defined",
            set->name, format_tuple(mpl, '[', tuple));
      /* add new member to the set and assign it empty elemental set */
      memb = add_member(mpl, set->array, tuple);
      memb->value.set = create_elemset(mpl, set->dimen);
      /* create an initial fake slice of all asterisks */
      slice = fake_slice(mpl, set->dimen);
      /* read zero or more data assignments */
      for (;;)
      {  /* skip optional comma */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
         /* process assignment element */
         if (mpl->token == T_ASSIGN)
         {  /* assignment ligature is non-significant element */
            get_token(mpl /* := */);
         }
         else if (mpl->token == T_LEFT)
         {  /* left parenthesis begins either new slice or "transpose"
               indicator */
            int is_tr;
            get_token(mpl /* ( */);
            is_tr = is_literal(mpl, "tr");
            unget_token(mpl /* ( */);
            if (is_tr) goto left;
            /* delete the current slice and read new one */
            delete_slice(mpl, slice);
            slice = read_slice(mpl, set->name, set->dimen);
            /* each new slice resets the "transpose" indicator */
            tr = 0;
            /* if the new slice is 0-ary, formally there is one 0-tuple
               (in the simple format) that follows it */
            if (slice_arity(mpl, slice) == 0)
               simple_format(mpl, set, memb, slice);
         }
         else if (is_symbol(mpl))
         {  /* number or symbol begins data in the simple format */
            simple_format(mpl, set, memb, slice);
         }
         else if (mpl->token == T_COLON)
         {  /* colon begins data in the matrix format */
            if (slice_arity(mpl, slice) != 2)
err1:          error(mpl, "slice currently used must specify 2 asterisk"
                  "s, not %d", slice_arity(mpl, slice));
            get_token(mpl /* : */);
            /* read elemental set data in the matrix format */
            matrix_format(mpl, set, memb, slice, tr);
         }
         else if (mpl->token == T_LEFT)
left:    {  /* left parenthesis begins the "transpose" indicator, which
               is followed by data in the matrix format */
            get_token(mpl /* ( */);
            if (!is_literal(mpl, "tr"))
err2:          error(mpl, "transpose indicator (tr) incomplete");
            if (slice_arity(mpl, slice) != 2) goto err1;
            get_token(mpl /* tr */);
            if (mpl->token != T_RIGHT) goto err2;
            get_token(mpl /* ) */);
            /* in this case the colon is optional */
            if (mpl->token == T_COLON) get_token(mpl /* : */);
            /* set the "transpose" indicator */
            tr = 1;
            /* read elemental set data in the matrix format */
            matrix_format(mpl, set, memb, slice, tr);
         }
         else if (mpl->token == T_SEMICOLON)
         {  /* semicolon terminates the data block */
            get_token(mpl /* ; */);
            break;
         }
         else
            error(mpl, "syntax error in set data block");
      }
      /* delete the current slice */
      delete_slice(mpl, slice);
      return;
}

/*----------------------------------------------------------------------
-- select_parameter - select parameter to saturate it with data.
--
-- This routine selects parameter to saturate it with data provided in
-- the data section. */

PARAMETER *select_parameter
(     MPL *mpl,
      char *name              /* not changed */
)
{     PARAMETER *par;
      AVLNODE *node;
      xassert(name != NULL);
      node = avl_find_node(mpl->tree, name);
      if (node == NULL || avl_get_node_type(node) != A_PARAMETER)
         error(mpl, "%s not a parameter", name);
      par = (PARAMETER *)avl_get_node_link(node);
      if (par->assign != NULL)
         error(mpl, "%s needs no data", name);
      if (par->data)
         error(mpl, "%s already provided with data", name);
      par->data = 1;
      return par;
}

/*----------------------------------------------------------------------
-- set_default - set default parameter value.
--
-- This routine sets default value for specified parameter. */

void set_default
(     MPL *mpl,
      PARAMETER *par,         /* not changed */
      SYMBOL *altval          /* destroyed */
)
{     xassert(par != NULL);
      xassert(altval != NULL);
      if (par->option != NULL)
         error(mpl, "default value for %s already specified in model se"
            "ction", par->name);
      xassert(par->defval == NULL);
      par->defval = altval;
      return;
}

/*----------------------------------------------------------------------
-- read_value - read value and assign it to parameter member.
--
-- This routine reads numeric or symbolic value from the input stream
-- and assigns to new parameter member specified by its n-tuple, which
-- (the member) is created and added to the parameter array. */

MEMBER *read_value
(     MPL *mpl,
      PARAMETER *par,         /* not changed */
      TUPLE *tuple            /* destroyed */
)
{     MEMBER *memb;
      xassert(par != NULL);
      xassert(is_symbol(mpl));
      /* there must be no member with the same n-tuple */
      if (find_member(mpl, par->array, tuple) != NULL)
         error(mpl, "%s%s already defined",
            par->name, format_tuple(mpl, '[', tuple));
      /* create new parameter member with given n-tuple */
      memb = add_member(mpl, par->array, tuple);
      /* read value and assigns it to the new parameter member */
      switch (par->type)
      {  case A_NUMERIC:
         case A_INTEGER:
         case A_BINARY:
            if (!is_number(mpl))
               error(mpl, "%s requires numeric data", par->name);
            memb->value.num = read_number(mpl);
            break;
         case A_SYMBOLIC:
            memb->value.sym = read_symbol(mpl);
            break;
         default:
            xassert(par != par);
      }
      return memb;
}

/*----------------------------------------------------------------------
-- plain_format - read parameter data block in plain format.
--
-- This routine reads parameter data block using the syntax:
--
-- <plain format> ::= <symbol> , <symbol> , ... , <symbol> , <value>
--
-- where <symbols> are used to determine a complete subscript list for
-- parameter member, <value> is a numeric or symbolic value assigned to
-- the parameter member. Commae between data items are optional and may
-- be omitted anywhere.
--
-- Number of components in the slice must be the same as dimension of
-- the parameter. To construct the complete subscript list the routine
-- replaces null positions in the slice by corresponding <symbols>. */

void plain_format
(     MPL *mpl,
      PARAMETER *par,         /* not changed */
      SLICE *slice            /* not changed */
)
{     TUPLE *tuple;
      SLICE *temp;
      SYMBOL *sym, *with = NULL;
      xassert(par != NULL);
      xassert(par->dim == slice_dimen(mpl, slice));
      xassert(is_symbol(mpl));
      /* read symbols and construct complete subscript list */
      tuple = create_tuple(mpl);
      for (temp = slice; temp != NULL; temp = temp->next)
      {  if (temp->sym == NULL)
         {  /* substitution is needed; read symbol */
            if (!is_symbol(mpl))
            {  int lack = slice_arity(mpl, temp) + 1;
               xassert(with != NULL);
               xassert(lack > 1);
               error(mpl, "%d items missing in data group beginning wit"
                  "h %s", lack, format_symbol(mpl, with));
            }
            sym = read_symbol(mpl);
            if (with == NULL) with = sym;
         }
         else
         {  /* copy symbol from the slice */
            sym = copy_symbol(mpl, temp->sym);
         }
         /* append the symbol to the subscript list */
         tuple = expand_tuple(mpl, tuple, sym);
         /* skip optional comma */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
      }
      /* read value and assign it to new parameter member */
      if (!is_symbol(mpl))
      {  xassert(with != NULL);
         error(mpl, "one item missing in data group beginning with %s",
            format_symbol(mpl, with));
      }
      read_value(mpl, par, tuple);
      return;
}

/*----------------------------------------------------------------------
-- tabular_format - read parameter data block in tabular format.
--
-- This routine reads parameter data block using the syntax:
--
-- <tabular format> ::= <column> <column> ... <column> :=
--                <row> <value>  <value>  ... <value>
--                <row> <value>  <value>  ... <value>
--                  .  .  .  .  .  .  .  .  .  .  .
--                <row> <value>  <value>  ... <value>
--
-- where <rows> are symbols that denote rows of the table, <columns>
-- are symbols that denote columns of the table, <values> are numeric
-- or symbolic values assigned to the corresponding parameter members.
-- If <value> is specified as single point, no value is provided.
--
-- Number of components in the slice must be the same as dimension of
-- the parameter. The slice must have two null positions. To construct
-- complete subscript list for particular <value> the routine replaces
-- the first null position of the slice by the corresponding <row> (or
-- <column>, if the flag tr is on) and the second null position by the
-- corresponding <column> (or by <row>, if the flag tr is on). */

void tabular_format
(     MPL *mpl,
      PARAMETER *par,         /* not changed */
      SLICE *slice,           /* not changed */
      int tr
)
{     SLICE *list, *col, *temp;
      TUPLE *tuple;
      SYMBOL *row;
      xassert(par != NULL);
      xassert(par->dim == slice_dimen(mpl, slice));
      xassert(slice_arity(mpl, slice) == 2);
      /* read the table heading that contains column symbols (the table
         may have no columns) */
      list = create_slice(mpl);
      while (mpl->token != T_ASSIGN)
      {  /* read column symbol and append it to the column list */
         if (!is_symbol(mpl))
            error(mpl, "number, symbol, or := missing where expected");
         list = expand_slice(mpl, list, read_symbol(mpl));
      }
      get_token(mpl /* := */);
      /* read zero or more rows that contain tabular data */
      while (is_symbol(mpl))
      {  /* read row symbol (if the table has no columns, these symbols
            are just ignored) */
         row = read_symbol(mpl);
         /* read values accordingly to the column list */
         for (col = list; col != NULL; col = col->next)
         {  int which = 0;
            /* if the token is single point, no value is provided */
            if (is_literal(mpl, "."))
            {  get_token(mpl /* . */);
               continue;
            }
            /* construct complete subscript list */
            tuple = create_tuple(mpl);
            for (temp = slice; temp != NULL; temp = temp->next)
            {  if (temp->sym == NULL)
               {  /* substitution is needed */
                  switch (++which)
                  {  case 1:
                        /* substitute in the first null position */
                        tuple = expand_tuple(mpl, tuple,
                           copy_symbol(mpl, tr ? col->sym : row));
                        break;
                     case 2:
                        /* substitute in the second null position */
                        tuple = expand_tuple(mpl, tuple,
                           copy_symbol(mpl, tr ? row : col->sym));
                        break;
                     default:
                        xassert(which != which);
                  }
               }
               else
               {  /* copy symbol from the slice */
                  tuple = expand_tuple(mpl, tuple, copy_symbol(mpl,
                     temp->sym));
               }
            }
            xassert(which == 2);
            /* read value and assign it to new parameter member */
            if (!is_symbol(mpl))
            {  int lack = slice_dimen(mpl, col);
               if (lack == 1)
                  error(mpl, "one item missing in data group beginning "
                     "with %s", format_symbol(mpl, row));
               else
                  error(mpl, "%d items missing in data group beginning "
                     "with %s", lack, format_symbol(mpl, row));
            }
            read_value(mpl, par, tuple);
         }
         /* delete the row symbol */
         delete_symbol(mpl, row);
      }
      /* delete the column list */
      delete_slice(mpl, list);
      return;
}

/*----------------------------------------------------------------------
-- tabbing_format - read parameter data block in tabbing format.
--
-- This routine reads parameter data block using the syntax:
--
-- <tabbing format> ::=  <prefix> <name>  , ... , <name>  , := ,
--    <symbol> , ... , <symbol> , <value> , ... , <value> ,
--    <symbol> , ... , <symbol> , <value> , ... , <value> ,
--     .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .  .
--    <symbol> , ... , <symbol> , <value> , ... , <value>
-- <prefix> ::= <empty>
-- <prefix> ::= <set name> :
--
-- where <names> are names of parameters (all the parameters must be
-- subscripted and have identical dimensions), <symbols> are symbols
-- used to define subscripts of parameter members, <values> are numeric
-- or symbolic values assigned to the corresponding parameter members.
-- Optional <prefix> may specify a simple set, in which case n-tuples
-- built of <symbols> for each row of the data table (i.e. subscripts
-- of parameter members) are added to the specified set. Commae between
-- data items are optional and may be omitted anywhere.
--
-- If the parameter altval is not NULL, it specifies a default value
-- provided for all the parameters specified in the data block.  */

void tabbing_format
(     MPL *mpl,
      SYMBOL *altval          /* not changed */
)
{     SET *set = NULL;
      PARAMETER *par;
      SLICE *list, *col;
      TUPLE *tuple;
      int next_token, j, dim = 0;
      char *last_name = NULL;
      /* read the optional <prefix> */
      if (is_symbol(mpl))
      {  get_token(mpl /* <symbol> */);
         next_token = mpl->token;
         unget_token(mpl /* <symbol> */);
         if (next_token == T_COLON)
         {  /* select the set to saturate it with data */
            set = select_set(mpl, mpl->image);
            /* the set must be simple (i.e. not set of sets) */
            if (set->dim != 0)
               error(mpl, "%s must be a simple set", set->name);
            /* and must not be defined yet */
            if (set->array->head != NULL)
               error(mpl, "%s already defined", set->name);
            /* add new (the only) member to the set and assign it empty
               elemental set */
            add_member(mpl, set->array, NULL)->value.set =
               create_elemset(mpl, set->dimen);
            last_name = set->name, dim = set->dimen;
            get_token(mpl /* <symbol> */);
            xassert(mpl->token == T_COLON);
            get_token(mpl /* : */);
         }
      }
      /* read the table heading that contains parameter names */
      list = create_slice(mpl);
      while (mpl->token != T_ASSIGN)
      {  /* there must be symbolic name of parameter */
         if (!is_symbol(mpl))
            error(mpl, "parameter name or := missing where expected");
         /* select the parameter to saturate it with data */
         par = select_parameter(mpl, mpl->image);
         /* the parameter must be subscripted */
         if (par->dim == 0)
            error(mpl, "%s not a subscripted parameter", mpl->image);
         /* the set (if specified) and all the parameters in the data
            block must have identical dimension */
         if (dim != 0 && par->dim != dim)
         {  xassert(last_name != NULL);
            error(mpl, "%s has dimension %d while %s has dimension %d",
               last_name, dim, par->name, par->dim);
         }
         /* set default value for the parameter (if specified) */
         if (altval != NULL)
            set_default(mpl, par, copy_symbol(mpl, altval));
         /* append the parameter to the column list */
         list = expand_slice(mpl, list, (SYMBOL *)par);
         last_name = par->name, dim = par->dim;
         get_token(mpl /* <symbol> */);
         /* skip optional comma */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
      }
      if (slice_dimen(mpl, list) == 0)
         error(mpl, "at least one parameter name required");
      get_token(mpl /* := */);
      /* skip optional comma */
      if (mpl->token == T_COMMA) get_token(mpl /* , */);
      /* read rows that contain tabbing data */
      while (is_symbol(mpl))
      {  /* read subscript list */
         tuple = create_tuple(mpl);
         for (j = 1; j <= dim; j++)
         {  /* read j-th subscript */
            if (!is_symbol(mpl))
            {  int lack = slice_dimen(mpl, list) + dim - j + 1;
               xassert(tuple != NULL);
               xassert(lack > 1);
               error(mpl, "%d items missing in data group beginning wit"
                  "h %s", lack, format_symbol(mpl, tuple->sym));
            }
            /* read and append j-th subscript to the n-tuple */
            tuple = expand_tuple(mpl, tuple, read_symbol(mpl));
            /* skip optional comma *between* <symbols> */
            if (j < dim && mpl->token == T_COMMA)
               get_token(mpl /* , */);
         }
         /* if the set is specified, add to it new n-tuple, which is a
            copy of the subscript list just read */
         if (set != NULL)
            check_then_add(mpl, set->array->head->value.set,
               copy_tuple(mpl, tuple));
         /* skip optional comma between <symbol> and <value> */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
         /* read values accordingly to the column list */
         for (col = list; col != NULL; col = col->next)
         {  /* if the token is single point, no value is provided */
            if (is_literal(mpl, "."))
            {  get_token(mpl /* . */);
               continue;
            }
            /* read value and assign it to new parameter member */
            if (!is_symbol(mpl))
            {  int lack = slice_dimen(mpl, col);
               xassert(tuple != NULL);
               if (lack == 1)
                  error(mpl, "one item missing in data group beginning "
                     "with %s", format_symbol(mpl, tuple->sym));
               else
                  error(mpl, "%d items missing in data group beginning "
                     "with %s", lack, format_symbol(mpl, tuple->sym));
            }
            read_value(mpl, (PARAMETER *)col->sym, copy_tuple(mpl,
               tuple));
            /* skip optional comma preceding the next value */
            if (col->next != NULL && mpl->token == T_COMMA)
               get_token(mpl /* , */);
         }
         /* delete the original subscript list */
         delete_tuple(mpl, tuple);
         /* skip optional comma (only if there is next data group) */
         if (mpl->token == T_COMMA)
         {  get_token(mpl /* , */);
            if (!is_symbol(mpl)) unget_token(mpl /* , */);
         }
      }
      /* delete the column list (it contains parameters, not symbols,
         so nullify it before) */
      for (col = list; col != NULL; col = col->next) col->sym = NULL;
      delete_slice(mpl, list);
      return;
}

/*----------------------------------------------------------------------
-- parameter_data - read parameter data.
--
-- This routine reads parameter data using the syntax:
--
-- <parameter data> ::= param <default value> : <tabbing format> ;
-- <parameter data> ::= param <parameter name> <default value>
--                      <assignments> ;
-- <parameter name> ::= <symbolic name>
-- <default value> ::= <empty>
-- <default value> ::= default <symbol>
-- <assignments> ::= <empty>
-- <assignments> ::= <assignments> , :=
-- <assignments> ::= <assignments> , [ <symbol list> ]
-- <assignments> ::= <assignments> , <plain format>
-- <assignemnts> ::= <assignments> , : <tabular format>
-- <assignments> ::= <assignments> , (tr) <tabular format>
-- <assignments> ::= <assignments> , (tr) : <tabular format>
--
-- Commae in <assignments> are optional and may be omitted anywhere. */

void parameter_data(MPL *mpl)
{     PARAMETER *par;
      SYMBOL *altval = NULL;
      SLICE *slice;
      int tr = 0;
      xassert(is_literal(mpl, "param"));
      get_token(mpl /* param */);
      /* read optional default value */
      if (is_literal(mpl, "default"))
      {  get_token(mpl /* default */);
         if (!is_symbol(mpl))
            error(mpl, "default value missing where expected");
         altval = read_symbol(mpl);
         /* if the default value follows the keyword 'param', the next
            token must be only the colon */
         if (mpl->token != T_COLON)
            error(mpl, "colon missing where expected");
      }
      /* being used after the keyword 'param' or the optional default
         value the colon begins data in the tabbing format */
      if (mpl->token == T_COLON)
      {  get_token(mpl /* : */);
         /* skip optional comma */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
         /* read parameter data in the tabbing format */
         tabbing_format(mpl, altval);
         /* on reading data in the tabbing format the default value is
            always copied, so delete the original symbol */
         if (altval != NULL) delete_symbol(mpl, altval);
         /* the next token must be only semicolon */
         if (mpl->token != T_SEMICOLON)
            error(mpl, "symbol, number, or semicolon missing where expe"
               "cted");
         get_token(mpl /* ; */);
         goto done;
      }
      /* in other cases there must be symbolic name of parameter, which
         follows the keyword 'param' */
      if (!is_symbol(mpl))
         error(mpl, "parameter name missing where expected");
      /* select the parameter to saturate it with data */
      par = select_parameter(mpl, mpl->image);
      get_token(mpl /* <symbol> */);
      /* read optional default value */
      if (is_literal(mpl, "default"))
      {  get_token(mpl /* default */);
         if (!is_symbol(mpl))
            error(mpl, "default value missing where expected");
         altval = read_symbol(mpl);
         /* set default value for the parameter */
         set_default(mpl, par, altval);
      }
      /* create initial fake slice of all asterisks */
      slice = fake_slice(mpl, par->dim);
      /* read zero or more data assignments */
      for (;;)
      {  /* skip optional comma */
         if (mpl->token == T_COMMA) get_token(mpl /* , */);
         /* process current assignment */
         if (mpl->token == T_ASSIGN)
         {  /* assignment ligature is non-significant element */
            get_token(mpl /* := */);
         }
         else if (mpl->token == T_LBRACKET)
         {  /* left bracket begins new slice; delete the current slice
               and read new one */
            delete_slice(mpl, slice);
            slice = read_slice(mpl, par->name, par->dim);
            /* each new slice resets the "transpose" indicator */
            tr = 0;
         }
         else if (is_symbol(mpl))
         {  /* number or symbol begins data in the plain format */
            plain_format(mpl, par, slice);
         }
         else if (mpl->token == T_COLON)
         {  /* colon begins data in the tabular format */
            if (par->dim == 0)
err1:          error(mpl, "%s not a subscripted parameter",
                  par->name);
            if (slice_arity(mpl, slice) != 2)
err2:          error(mpl, "slice currently used must specify 2 asterisk"
                  "s, not %d", slice_arity(mpl, slice));
            get_token(mpl /* : */);
            /* read parameter data in the tabular format */
            tabular_format(mpl, par, slice, tr);
         }
         else if (mpl->token == T_LEFT)
         {  /* left parenthesis begins the "transpose" indicator, which
               is followed by data in the tabular format */
            get_token(mpl /* ( */);
            if (!is_literal(mpl, "tr"))
err3:          error(mpl, "transpose indicator (tr) incomplete");
            if (par->dim == 0) goto err1;
            if (slice_arity(mpl, slice) != 2) goto err2;
            get_token(mpl /* tr */);
            if (mpl->token != T_RIGHT) goto err3;
            get_token(mpl /* ) */);
            /* in this case the colon is optional */
            if (mpl->token == T_COLON) get_token(mpl /* : */);
            /* set the "transpose" indicator */
            tr = 1;
            /* read parameter data in the tabular format */
            tabular_format(mpl, par, slice, tr);
         }
         else if (mpl->token == T_SEMICOLON)
         {  /* semicolon terminates the data block */
            get_token(mpl /* ; */);
            break;
         }
         else
            error(mpl, "syntax error in parameter data block");
      }
      /* delete the current slice */
      delete_slice(mpl, slice);
done: return;
}

/*----------------------------------------------------------------------
-- data_section - read data section.
--
-- This routine reads data section using the syntax:
--
-- <data section> ::= <empty>
-- <data section> ::= <data section> <data block> ;
-- <data block> ::= <set data>
-- <data block> ::= <parameter data>
--
-- Reading data section is terminated by either the keyword 'end' or
-- the end of file. */

void data_section(MPL *mpl)
{     while (!(mpl->token == T_EOF || is_literal(mpl, "end")))
      {  if (is_literal(mpl, "set"))
            set_data(mpl);
         else if (is_literal(mpl, "param"))
            parameter_data(mpl);
         else
            error(mpl, "syntax error in data section");
      }
      return;
}

/* eof */

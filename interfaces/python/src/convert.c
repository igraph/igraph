/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
   MTA RMKI, Konkoly-Thege Miklos st. 29-33, Budapest 1121, Hungary
   
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
   Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
   02110-1301 USA

*/

/************************ Miscellaneous functions *************************/

#include "graphobject.h"
#include "convert.h"
#include "error.h"
#include "memory.h"

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python integer list to an igraph \c igraph_vector_t
 * The incoming \c igraph_vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c igraph_vector_t containing the result
 * \param need_non_negative if true, checks whether all elements are non-negative
 * \param pairs if true, assumes that every list element is a pair of integers
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyList_to_vector_t(PyObject *list, igraph_vector_t *v, igraph_bool_t need_non_negative, igraph_bool_t pairs)
{
   PyObject *item, *i1, *i2;
   int i, j, k, ok;
   long idx=0, idx2=0;

   if (!PyList_Check(list)) 
     {
	ok=1;
	
	if (pairs)
	  {
	     // a pair was given instead of a list
	     // Let's assume that the user meant a list consisting of this single pair
	     if (PyTuple_Check(list) && PyTuple_Size(list)==2)
	       {
		  i1=i2=NULL;
		  i1=PyTuple_GetItem(list, 0);
		  if (i1) i2=PyTuple_GetItem(list, 1);
		  if (i1 && i2) 
		    {
		       if (PyInt_Check(i1) && PyInt_Check(i2)) 
			 {
			    idx=PyInt_AsLong(i1); idx2=PyInt_AsLong(i2);
			    // domain checking
			    if (need_non_negative && (idx<0 || idx2<0)) ok=0;
			 }
		       else ok=0;
		    }
		  else
		    {
		       // should not ever get here
		       // Exception was set by PyTuple_GetItem, so return immediately
		       return 1;
		    }
		  
	       }
	     else 
	       {
		  PyErr_SetString(PyExc_TypeError, "List of pairs or a single pair of integer expected");
		  return 1;
	       }
	     
	     if (ok) 
	       {
		  igraph_vector_init(v, 2);
		  VECTOR(*v)[0]=(igraph_real_t)idx;
		  VECTOR(*v)[1]=(igraph_real_t)idx2;
	       }
	     else
	       {
		  if (need_non_negative)
		    PyErr_SetString(PyExc_TypeError, "List elements must be non-negative integers");
		  else
		    PyErr_SetString(PyExc_TypeError, "List elements must be integers");
		  return 1;
	       }
	  }
	else if (PyInt_Check(list)) 
	  {
	     // a single integer was given instead of a list
	     // Let's assume that the user meant a list consisting of this single item
	     igraph_vector_init(v, 1);
	     VECTOR(*v)[0]=(igraph_real_t)PyInt_AsLong(list);
	  }
	else 
	  {
	     PyErr_SetString(PyExc_TypeError, "List or single integer expected");
	     return 1;
	  }
     }
   else 
     {
	// we received a list, so loop through all elements and add them
	// to a vector.
	// Non-integer or negative elements raise an exception
	j=PyList_Size(list);
	if (pairs)
	  igraph_vector_init(v, 2*j);
	else
	  igraph_vector_init(v, j);
	for (i=0, k=0; i<j; i++)
	  {
	     item=PyList_GetItem(list, i);
	     if (item) 
	       {
		  ok=1;
		  
		  if (pairs) 
		    {
		       // item is a borrowed reference, so no need to hassle with reference counting
		       if (PyTuple_Check(item) && PyTuple_Size(item)==2)
			 {
			    i1=NULL; i2=NULL;
			    i1=PyTuple_GetItem(item, 0);
			    if (i1) i2=PyTuple_GetItem(item, 1);
			    if (i1 && i2) 
			      {
				 if (PyInt_Check(i1) && PyInt_Check(i2)) 
				   {
				      // domain checking
				      idx=PyInt_AsLong(i1);
				      idx2=PyInt_AsLong(i2);
				      if (need_non_negative && (idx<0 || idx2<0)) ok=0;
				   }
				 else ok=0;
			      }
			    else 
			      {
				 // this should not happen, but we return anyway.
				 // an IndexError exception was set by PyTuple_GetItem
				 // at this point
				 igraph_vector_destroy(v);
				 return 1;
			      }
			    
			 } else ok=0;
		    }
		  else 
		    {
		       // item is a borrowed reference, so no need to hassle with reference counting
		       if (PyInt_Check(item)) 
			 {
			    // domain checking
			    idx=PyInt_AsLong(item);
			    if (need_non_negative && idx<0) ok=0;
			 }
		       else ok=0;
		    }
	     
		  if (!ok)
		    {
		      if (pairs) {
			// item is not a non-negative integer pair, so throw an exception
			if (need_non_negative)
			  PyErr_SetString(PyExc_TypeError, "List elements must be non-negative integer pairs");
			else
			  PyErr_SetString(PyExc_TypeError, "List elements must be integer pairs");
			igraph_vector_destroy(v);
			return 1;
		      } else {
			// item is not a non-negative integer, so throw an exception
			if (need_non_negative)
			  PyErr_SetString(PyExc_TypeError, "List elements must be non-negative integers");
			else
			  PyErr_SetString(PyExc_TypeError, "List elements must be integers");
			igraph_vector_destroy(v);
			return 1;
		      }
		    }
		  
		 // add idx into index vector
		  VECTOR(*v)[k]=(igraph_real_t)idx;
		  k++;
		  if (pairs) 
		    {
		       // if we are working on pairs, add idx and idx2 as well
		       VECTOR(*v)[k]=(igraph_real_t)idx2;
		       k++;
		    }
	       }
	     else 
	       {
		  // this should not happen, but we return anyway.
		  // an IndexError exception was set by PyList_GetItem
		  // at this point
		  igraph_vector_destroy(v);
		  return 1;
	       }
	  }
     }
   
   return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python integer list
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList(igraph_vector_t *v) {
  PyObject* list;
  int n, i;
   
  n=igraph_vector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();

  list=PyList_New(n);
  for (i=0; i<n; i++) {
    PyObject *item=PyInt_FromLong((long)VECTOR(*v)[i]);
    if (!item) {
      Py_DECREF(list);
      return NULL;
    }
    PyList_SET_ITEM(list, i, item);
  }

  return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python integer tuple
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer tuple as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyTuple(igraph_vector_t *v) {
  PyObject* tuple;
  int n, i;
  
  n=igraph_vector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();
  
  tuple=PyTuple_New(n);
  for (i=0; i<n; i++) {
    PyObject *item=PyInt_FromLong((long)VECTOR(*v)[i]);
    if (!item) {
      Py_DECREF(tuple);
      return NULL;
    }
    PyTuple_SET_ITEM(tuple, i, item);
  }

  return tuple;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python list of integer pairs
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList_pairs(igraph_vector_t *v) {
   PyObject *list, *pair;
   long n, i, j;
   
   n=igraph_vector_size(v);
   if (n<0) return igraphmodule_handle_igraph_error();
   if (n%2) return igraphmodule_handle_igraph_error();
   
   /* create a new Python list */
   n>>=1;
   list=PyList_New(n);
   
   /* populate the list with data */
   for (i=0, j=0; i<n; i++, j+=2) {
     pair=Py_BuildValue("(ll)", (long)VECTOR(*v)[j], (long)VECTOR(*v)[j+1]);
     if (pair==NULL || PyList_SetItem(list, i, pair)) {
       /* error occurred while populating the list, return immediately */
       Py_DECREF(pair);
       Py_DECREF(list);
       return NULL;
     }
   }

   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts two igraph \c igraph_vector_t vectors to a Python list of integer pairs
 * 
 * \param v1 the \c igraph_vector_t containing the 1st vector to be converted
 * \param v2 the \c igraph_vector_t containing the 2nd vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_pair_to_PyList(igraph_vector_t *v1,
					       igraph_vector_t *v2) {
   PyObject *list, *pair;
   long n, i;
   
   n=igraph_vector_size(v1);
   if (n<0) return igraphmodule_handle_igraph_error();
   if (igraph_vector_size(v2) != n) return igraphmodule_handle_igraph_error();

   /* create a new Python list */
   list=PyList_New(n);
   
   /* populate the list with data */
   for (i=0; i<n; i++) {
     pair=Py_BuildValue("(ll)", (long)VECTOR(*v1)[i], (long)VECTOR(*v2)[i]);
     if (pair==NULL || PyList_SetItem(list, i, pair)) {
       /* error occurred while populating the list, return immediately */
       Py_DECREF(pair);
       Py_DECREF(list);
       return NULL;
     }
   }

   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_vector_t to a Python float list
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python float list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_float_PyList(igraph_vector_t *v) {
   PyObject* list;
   int n, i;
   
   n=igraph_vector_size(v);
   if (n<0) return igraphmodule_handle_igraph_error();

   // create a new Python list
   list=PyList_New(igraph_vector_size(v));
   // populate the list with data
   for (i=0; i<n; i++) 
     {
	if (PyList_SetItem(list, i, PyFloat_FromDouble(VECTOR(*v)[i]))) 
	  {
	     // error occurred while populating the list, return immediately
	     Py_DECREF(list);
	     return NULL;
	  }
     }
   // return the list
   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c igraph_matrix_t to a Python list of lists
 * 
 * \param m the \c igraph_matrix_t containing the matrix to be converted
 * \param type the type of conversion. If equals to IGRAPHMODULE_TYPE_INT,
 *        returns an integer matrix, else returns a float matrix.
 * \return the Python list of lists as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_matrix_t_to_PyList(igraph_matrix_t *m,
						 igraphmodule_conv_t type) {
   PyObject *list, *row, *item;
   int nr, nc, i, j;
   
   
   nr=igraph_matrix_nrow(m); nc=igraph_matrix_ncol(m);
   if (nr<0 || nc<0) return igraphmodule_handle_igraph_error();

   // create a new Python list
   list=PyList_New(nr);
   // populate the list with data
   for (i=0; i<nr; i++) 
     {
	row=PyList_New(nc);
	for (j=0; j<nc; j++) 
	  {
	     if (type==IGRAPHMODULE_TYPE_INT)
	       item=PyInt_FromLong(MATRIX(*m,i,j));
	     else
	       item=PyFloat_FromDouble(MATRIX(*m,i,j));
	       
	     if (PyList_SetItem(row, j, item))
	       {
		  // error occurred while populating the list, return immediately
		  Py_DECREF(row);
		  Py_DECREF(list);
		  return NULL;
	       }
	  }
	if (PyList_SetItem(list, i, row)) 
	  {
	     Py_DECREF(row);
	     Py_DECREF(list);
	     return NULL;
	  }
     }
   // return the list
   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python list of lists to an \c igraph_matrix_t
 * 
 * \param o the Python object representing the list of lists
 * \param m the address of an uninitialized \c igraph_matrix_t
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyList_to_matrix_t(PyObject* o, igraph_matrix_t *m) {
  int nr, nc, n, i, j;
  PyObject *row, *item;
  
  /* calculate the matrix dimensions */
  if (!PyList_Check(o)) return 1;
  nr = PyList_Size(o);
  nc = 0;
  for (i=0; i<nr; i++) {
    row=PyList_GetItem(o, i);
    if (!PyList_Check(row)) return 1;
    n=PyList_Size(row);
    if (n>nc) nc=n;
  }
  
  igraph_matrix_init(m, nr, nc);
  for (i=0; i<nr; i++) {
    row=PyList_GetItem(o, i);
    n=PyList_Size(row);
    for (j=0; j<n; j++) {
      item=PyList_GetItem(row, j);
      if (PyInt_Check(item)) {
	MATRIX(*m, i, j) = (igraph_real_t)PyInt_AsLong(item);
      } else if (PyLong_Check(item)) {
	MATRIX(*m, i, j) = (igraph_real_t)PyLong_AsLong(item);
      } else if (PyFloat_Check(item)) {
	MATRIX(*m, i, j) = (igraph_real_t)PyFloat_AsDouble(item);
      } else {
	/* warning? */
      }
    }
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an \c igraph_strvector_t to a Python string list
 * 
 * \param v the \c igraph_strvector_t containing the vector to be converted
 * \return the Python string list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_strvector_t_to_PyList(igraph_strvector_t *v) {
  PyObject* list;
  int n, i;
  char* ptr;
  
  n=igraph_strvector_size(v);
  if (n<0) return igraphmodule_handle_igraph_error();
  
  // create a new Python list
  list=PyList_New(n);
  /* populate the list with data */
  for (i=0; i<n; i++) {
    igraph_strvector_get(v, i, &ptr);
    if (PyList_SetItem(list, i, PyString_FromString(ptr))) {
      /* error occurred while populating the list, return immediately */
      Py_DECREF(list);
      return NULL;
    }
  }
  
   /* return the list */
   return list;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python string list to an \c igraph_strvector_t
 * 
 * \param v the Python list as a \c PyObject*
 * \param result an \c igraph_strvector_t containing the result
 * The incoming \c igraph_strvector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyList_to_strvector_t(PyObject* v, igraph_strvector_t *result) {
  int n, i;
  static char* empty = "";
  
  if (!PyList_Check(v)) {
    PyErr_SetString(PyExc_TypeError, "expected list");
    return 1;
  }
  
  n=PyList_Size(v);

  /* initialize the string vector */
  if (igraph_strvector_init(result, n)) return 1;
  
  /* populate the vector with data */
  for (i=0; i<n; i++) {
    PyObject *item = PyList_GetItem(v, i);
    char* ptr;
    if (PyString_Check(item))
      ptr=PyString_AS_STRING(item);
    else
      ptr=empty;
    if (igraph_strvector_set(result, i, ptr)) {
      igraph_strvector_destroy(result);
      return 1;
    }
  }

  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python iterator returning graphs to an \c igraph_vectorptr_t
 * The incoming \c igraph_vector_ptr_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param it the Python iterator
 * \param v the \c igraph_vector_ptr_t which will contain the result
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyIter_to_vector_ptr_t(PyObject *it, igraph_vector_ptr_t *v) {
  PyObject *t;
  
  if (igraph_vector_ptr_init(v, 0)) {
    igraphmodule_handle_igraph_error();
    return 1;
  }
    
  while ((t=PyIter_Next(it))) {
    if (!PyObject_TypeCheck(t, &igraphmodule_GraphType)) {
      PyErr_SetString(PyExc_TypeError, "iterable argument must contain graphs");
      igraph_vector_ptr_destroy(v);
      Py_DECREF(t);
      Py_DECREF(it);
      return 1;
    }
    igraph_vector_ptr_push_back(v, &((igraphmodule_GraphObject*)t)->g);
    Py_DECREF(t);
  }  
  
  return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a vertex selector
 * 
 * \param o the Python object
 * \param vs the \c igraph_vs_t which will contain the result
 * \param return_single will be 1 if the selector selected only a single vertex,
 * 0 otherwise
 * \return 0 if everything was OK, 1 otherwise
 */
int igraphmodule_PyObject_to_vs_t(PyObject *o, igraph_vs_t *vs,
				  igraph_bool_t *return_single) {
  *return_single=0;
  if (o==NULL || o == Py_None) {
    /* Returns a vertex sequence for all vertices */
    igraph_vs_all(vs);
  } else if (PyInt_Check(o)) {
    /* Returns a vertex sequence for a single vertex ID */
    igraph_vs_1(vs, PyInt_AsLong(o));
    *return_single=1;
  } else if (PyLong_Check(o)) {
    /* Returns a vertex sequence for a single vertex ID */
    igraph_vs_1(vs, PyLong_AsLong(o));
    *return_single=1;
  } else {
    /* Returns a vertex sequence with the IDs returned by the iterator */
    PyObject *iterator = PyObject_GetIter(o);
    PyObject *item;
    igraph_vector_t vector;

    if (iterator == NULL) {
      PyErr_SetString(PyExc_TypeError, "integer, long, iterable or None expected");
      return 1;
    }

    IGRAPH_CHECK(igraph_vector_init(&vector, 0));
    IGRAPH_FINALLY(igraph_vector_destroy, &vector);
    IGRAPH_CHECK(igraph_vector_reserve(&vector, 20));

    while ((item = PyIter_Next(iterator))) {
      long val=-1;
      if (PyInt_Check(item)) val=PyInt_AsLong(item);
      else if (PyLong_Check(item)) val=PyLong_AsLong(item);
      Py_DECREF(item);

      if (val >= 0) igraph_vector_push_back(&vector, val);
      else {
	PyErr_SetString(PyExc_TypeError, "integer or long expected");
      }

      if (PyErr_Occurred()) break;
    }
    Py_DECREF(iterator);

    if (PyErr_Occurred()) {
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
      return 1;
    } else {
      igraph_vs_vector_copy(vs, &vector);
      igraph_vector_destroy(&vector);
      IGRAPH_FINALLY_CLEAN(1);
    }
  }
  return 0;
}


/**
 * \ingroup python_interface_conversion
 * \brief Tries to interpret a Python object as a numeric attribute value list
 * 
 * \param o the Python object
 * \param v the \c igraph_vector_t which will contain the result
 * \param g a \c igraphmodule_GraphObject object or \c NULL - used when the
 * provided Python object is not a list and we're trying to interpret it as
 * an attribute name.
 * \param type the attribute type (graph = 0, vertex = 1, edge = 2) to be used
 * \param def default value if the attribute name supplied is \c None
 * if \c o is not a list.
 * \return 0 if everything was OK, 1 otherwise
 *
 * If the Python object is not a list, tries to interpret it as an attribute
 * name.
 */
int igraphmodule_PyObject_to_attribute_values(PyObject *o,
					      igraph_vector_t *v,
					      igraphmodule_GraphObject* g,
					      int type, igraph_real_t def) {
  PyObject* list = o;
  long i, n;

  if (o==NULL) return 1;
  
  if (o == Py_None) {
    if (type == ATTRHASH_IDX_VERTEX) n=igraph_vcount(&g->g);
    else if (type == ATTRHASH_IDX_EDGE) n=igraph_ecount(&g->g);
    else n=1;

    if (igraph_vector_init(v, n)) return 1;
    for (i=0; i<n; i++) VECTOR(*v)[i] = def;
    return 0;
  }

  if (!PyList_Check(o)) {
    list = PyDict_GetItem(((PyObject**)g->g.attr)[type], o);
    if (!list) {
      if (!PyErr_Occurred())
	PyErr_SetString(PyExc_KeyError, "Attribute does not exist");
      return 1;
    }
  }

  n=PyList_Size(list);
  if (igraph_vector_init(v, n)) return 1;

  for (i=0; i<n; i++) {
    PyObject *item = PyList_GetItem(list, i);
    if (!item) {
      igraph_vector_destroy(v);
      return 1;
    }

    if (PyInt_Check(item))
      VECTOR(*v)[i] = PyInt_AsLong(item);
    else if (PyLong_Check(item))
      VECTOR(*v)[i] = PyLong_AsLong(item);
    else if (PyFloat_Check(item))
      VECTOR(*v)[i] = PyFloat_AsDouble(item);
    else
      VECTOR(*v)[i] = def;
  }

  return 0;
}

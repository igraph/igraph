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
int igraphmodule_PyList_to_vector_t(PyObject *list, igraph_vector_t *v, bool_t need_non_negative, bool_t pairs)
{
   PyObject *item, *i1, *i2;
   int i, j, k, ok;
   long idx, idx2=0;

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
		  VECTOR(*v)[0]=(real_t)idx;
		  VECTOR(*v)[1]=(real_t)idx2;
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
	     VECTOR(*v)[0]=(real_t)PyInt_AsLong(list);
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
		  VECTOR(*v)[k]=(real_t)idx;
		  k++;
		  if (pairs) 
		    {
		       // if we are working on pairs, add idx and idx2 as well
		       VECTOR(*v)[k]=(real_t)idx2;
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

   // create a new Python list
   list=PyList_New(n);
   // populate the list with data
   for (i=0; i<n; i++) 
     {
	if (PyList_SetItem(list, i, PyInt_FromLong(VECTOR(*v)[i]))) 
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
 * \brief Converts an igraph \c igraph_vector_t to a Python list of integer pairs
 * 
 * \param v the \c igraph_vector_t containing the vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
PyObject* igraphmodule_vector_t_to_PyList_pairs(igraph_vector_t *v) {
   PyObject *list, *pair;
   int n, i, j;
   
   n=igraph_vector_size(v);
   if (n<0) return igraphmodule_handle_igraph_error();
   if (n%2) return igraphmodule_handle_igraph_error();
   
   // create a new Python list
   n>>=1;
   list=PyList_New(n);
   
   // populate the list with data
   for (i=0, j=0; i<n; i++, j+=2)
     {
	pair=Py_BuildValue("(ll)", (long)VECTOR(*v)[j], (long)VECTOR(*v)[j+1]);
	if (pair==NULL || PyList_SetItem(list, i, pair)) 
	  {
	     // error occurred while populating the list, return immediately
	     Py_DECREF(pair);
	     Py_DECREF(list);
	     return NULL;
	  }
     }
   // return the list
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
  // populate the list with data
  for (i=0; i<n; i++) {
    igraph_strvector_get(v, i, &ptr);
    if (PyList_SetItem(list, i, PyString_FromString(ptr))) {
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

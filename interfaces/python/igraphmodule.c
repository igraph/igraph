#include <Python.h>
#include "igraph.h"

/**
 * \defgroup python_interface Python module implementation
 * \brief Functions implementing a Python interface to \a igraph
 * 
 * These functions provide a way to access \a igraph functions from Python.
 * It should be of interest of \a igraph developers only. Classes, functions
 * and methods exposed to Python are still to be documented. Until it is done,
 * just type the following to get help about \a igraph functions in Python
 * (assuming you have \c igraph.so somewhere in your Python library path):
 * 
 * \verbatim
import igraph
help(igraph)
help(igraph.Graph)
\endverbatim
 * 
 * Most of the functions provided here share the same calling conventions
 * (which are determined by the Python/C API). Since the role of the
 * arguments are the same across many functions, I won't explain them
 * everywhere, just give a quick overview of the common argument names here.
 * 
 * \param self the Python igraph.Graph object the method is working on
 * \param args pointer to the Python tuple containing the arguments
 * \param kwds pointer to the Python hash containing the keyword parameters
 * \param type the type object of a Python igraph.Graph object. Used usually
 * in constructors and class methods.
 * 
 * Any arguments not documented here should be mentioned at the documentation
 * of the appropriate method.
 * 
 * The functions which implement a Python method always return a pointer to
 * a \c PyObject. According to Python conventions, this is \c NULL if and
 * only if an exception was thrown by the method (or any of the functions
 * it has called). When I explain the return value of a function which
 * provides interface to an \a igraph function, I won't cover the case of
 * returning a \c NULL value, because this is the same for every such method.
 * The conclusion is that a method can return \c NULL even if I don't state
 * it explicitly.
 * 
 * Also please take into consideration that I'm documenting the C calls
 * with the abovementioned parameters here, and \em not the Python methods
 * which are presented to the user using the Python interface of \a igraph.
 * If you are looking for the documentation of the classes, methods and
 * functions exposed to Python, please use the \c help calls from Python
 * or use \c pydoc to generate a formatted version.
 */

/********************** Internal types *********************************/
typedef enum { IGRAPHMODULE_TYPE_INT=0, IGRAPHMODULE_TYPE_FLOAT }
igraphmodule_conv_t;

/********************** Error handling functions **************************/

/** \defgroup python_interface_errors Error handling
 * \ingroup python_interface */

/** \ingroup python_interface_errors
 * \brief Exception type to be returned when an internal \c igraph error occurs.
 */
static PyObject* igraphmodule_InternalError;

/**
 * \ingroup python_interface_errors
 * \brief Generic error handler for internal \c igraph errors.
 * 
 * Since the \c igraph functions don't provide any clue regarding what went
 * wrong (at least at the time of writing), this function simply generates
 * a generic exception. It should be extended as soon as some "real" error
 * handling is provided by \c igraph.
 * 
 * \return Always returns \c NULL, and all callers are advised to pass this
 * \c NULL value to their callers until it is propagated to the Python
 * interpreter.
 */
PyObject* igraphmodule_handle_igraph_error() 
{
   PyErr_SetString(igraphmodule_InternalError,
		   "Internal igraph error. Please contact the author!");
   return NULL;
}

/************************ Miscellaneous functions *************************/

/** \defgroup python_interface_conversion Converting between Python and igraph data types
 * \ingroup python_interface */

/**
 * \ingroup python_interface_conversion
 * \brief Converts a Python integer list to an igraph \c vector_t
 * The incoming \c vector_t should be uninitialized. Raises suitable
 * Python exceptions when needed.
 * 
 * \param list the Python list to be converted
 * \param v the \c vector_t containing the result
 * \param need_non_negative if true, checks whether all elements are non-negative
 * \param pairs if true, assumes that every list element is a pair of integers
 * \return 0 if everything was OK, 1 otherwise
 */
static int igraphmodule_PyList_to_vector_t(PyObject *list, vector_t *v, bool_t need_non_negative, bool_t pairs)
{
   PyObject *item, *i1, *i2;
   int i, j, k, ok;
   long idx, idx2;

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
		  i1=PyTuple_GetItem(list, 1);
		  if (i1) i2=PyTuple_GetItem(list, 2);
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
		  vector_init(v, 2);
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
	     vector_init(v, 1);
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
	  vector_init(v, 2*j);
	else
	  vector_init(v, j);
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
				 vector_destroy(v);
				 return 1;
			      }
			    
			 }
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
		       // item is not a non-negative integer, so throw an exception
		       if (need_non_negative)
			 PyErr_SetString(PyExc_TypeError, "List elements must be non-negative integers");
		       else
			 PyErr_SetString(PyExc_TypeError, "List elements must be integers");
		       vector_destroy(v);
		       return 1;
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
		  vector_destroy(v);
		  return 1;
	       }
	  }
     }
   
   return 0;
}

/**
 * \ingroup python_interface_conversion
 * \brief Converts an igraph \c vector_t to a Python integer list
 * 
 * \param v the \c vector_t containing the vector to be converted
 * \return the Python integer list as a \c PyObject*, or \c NULL if an error occurred
 */
static PyObject* igraphmodule_vector_t_to_PyList(vector_t *v) {
   PyObject* list;
   int n, i;
   
   n=vector_size(v);
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
 * \brief Converts an igraph \c vector_t to a Python list of integer pairs
 * 
 * \param v the \c vector_t containing the vector to be converted
 * \return the Python integer pair list as a \c PyObject*, or \c NULL if an error occurred
 */
static PyObject* igraphmodule_vector_t_to_PyList_pairs(vector_t *v) {
   PyObject *list, *pair;
   int n, i, j;
   
   n=vector_size(v);
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
 * \brief Converts an igraph \c vector_t to a Python float list
 * 
 * \param v the \c vector_t containing the vector to be converted
 * \return the Python float list as a \c PyObject*, or \c NULL if an error occurred
 */
static PyObject* igraphmodule_vector_t_to_float_PyList(vector_t *v) {
   PyObject* list;
   int n, i;
   
   n=vector_size(v);
   if (n<0) return igraphmodule_handle_igraph_error();

   // create a new Python list
   list=PyList_New(vector_size(v));
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
 * \brief Converts an igraph \c matrix_t to a Python list of lists
 * 
 * \param m the \c matrix_t containing the matrix to be converted
 * \param type the type of conversion. If equals to IGRAPHMODULE_TYPE_INT,
 *        returns an integer matrix, else returns a float matrix.
 * \return the Python list of lists as a \c PyObject*, or \c NULL if an error occurred
 */
static PyObject* igraphmodule_matrix_t_to_PyList(matrix_t *m,
						 igraphmodule_conv_t type) {
   PyObject *list, *row, *item;
   int nr, nc, i, j;
   
   
   nr=matrix_nrow(m); nc=matrix_ncol(m);
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
 * \ingroup python_interface
 * \brief Handler function for all unimplemented \c igraph.Graph methods
 * 
 * This function is called whenever an unimplemented \c igraph.Graph method
 * is called ("unimplemented" meaning that there is a method name in the
 * method table of \c igraph.Graph , but there isn't any working implementation
 * either because the underlying \c igraph API might be subject to change
 * or because the calling format from Python is not decided yet (or maybe
 * because of laziness or lack of time ;))
 * 
 * All of the parameters are ignored, they are here just to make the
 * function satisfy the requirements of \c PyCFunction, thus allowing it
 * to be included in a method table.
 * 
 * \return NULL
 */
static PyObject* igraphmodule_unimplemented(PyObject* self, PyObject* args, PyObject* kwds)
{
   PyErr_SetString(PyExc_NotImplementedError, "This method is unimplemented.");
   return NULL;
}

/**************************** igraph object *******************************/

/**
 * \ingroup python_interface
 * \brief A structure containing all the fields required to access an igraph from Python
 */
typedef struct 
{
  PyObject_HEAD
  igraph_t g;             // The graph object
  PyObject* destructor;   // Python object to be called upon destruction
} igraphmodule_GraphObject;

/**
 * \ingroup python_interface
 * \brief Deallocates a Python representation of a given igraph object
 */
static void igraphmodule_Graph_dealloc(igraphmodule_GraphObject* self) 
{
  PyObject* r;
  
  igraph_destroy(&self->g);
  if (PyCallable_Check(self->destructor)) {
    r=PyObject_CallObject(self->destructor, NULL);
    if (r) Py_DECREF(r);
  }
  self->ob_type->tp_free((PyObject*)self);
}

/**
 * \ingroup python_interface
 * \brief Creates a new igraph object in Python
 * 
 * This function is called whenever a new \c igraph.Graph object is created in
 * Python. An optional \c n parameter can be passed from Python,
 * representing the number of vertices in the graph. If it is omitted,
 * the default value is 1.
 * 
 * <b>Example call from Python:</b>
\verbatim
g = igraph.Graph(5);
\endverbatim
 *
 * In fact, the parameters are processed by \c igraphmodule_Graph_init
 * 
 * \return the new \c igraph.Graph object or NULL if an error occurred.
 * 
 * \sa igraphmodule_Graph_init
 * \sa igraph_empty
 */
static PyObject* igraphmodule_Graph_new(PyTypeObject *type,
					PyObject *args,
					PyObject *kwds) 
{
   igraphmodule_GraphObject *self;
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	igraph_empty(&self->g, 1, 0);
     }
   
   return (PyObject*)self;
}

/**
 * \ingroup python_interface
 * \brief Initializes a new \c igraph object in Python
 * 
 * This function is called whenever a new \c igraph.Graph object is initialized in
 * Python (note that initializing is not equal to creating: an object might
 * be created but not initialized when it is being recovered from a serialized
 * state).
 * 
 * Throws \c AssertionError in Python if \c vcount is less than or equal to zero.
 * \return the new \c igraph.Graph object or NULL if an error occurred.
 * 
 * \sa igraphmodule_Graph_new
 * \sa igraph_empty
 * \sa igraph_create
 */
static int igraphmodule_Graph_init(igraphmodule_GraphObject *self,
				   PyObject *args,
				   PyObject *kwds)
{
   static char *kwlist[] = 
     {
	"n", "edges", "directed", NULL
     }
   ;
   int n=1;
   PyObject *edges=NULL, *dir=Py_False;
   vector_t edges_vector;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|iO!O!", kwlist,
				    &n, &PyList_Type, &edges,
				    &PyBool_Type, &dir))
     return -1;

   if (n<0) 
     {
	// Throw an exception
	PyErr_SetString(PyExc_AssertionError, "Number of vertices must be greater than zero.");
	return -1;
     }
   
   if (edges && PyList_Check(edges))
     {
	// Caller specified an edge list, so we use igraph_create
	// We have to convert the Python list to a vector_t
	if (igraphmodule_PyList_to_vector_t(edges, &edges_vector, 1, 1)) 
	  {
	     igraphmodule_handle_igraph_error();
	     return -1;
	  }
	
	igraph_destroy(&self->g);
	igraph_create(&self->g, &edges_vector, (integer_t)n, (dir==Py_True));
	
	vector_destroy(&edges_vector);
     }
   else
     {
	// No edge list was specified, let's use igraph_empty
	igraph_destroy(&self->g);
	igraph_empty(&self->g, n, (dir==Py_True));
     }
   
   return 0;
}

/** \ingroup python_interface
 * \brief Formats an \c igraph.Graph object in a human-consumable format.
 * 
 * This function is rather simple now, it returns the number of vertices
 * and edges in a string.
 * 
 * \return the formatted textual representation as a \c PyObject
 */
static PyObject* igraphmodule_Graph_str(igraphmodule_GraphObject *self)
{
   if (igraph_is_directed(&self->g))
     return PyString_FromFormat("Directed graph (|V| = %ld, |E| = %ld)",
				(long)igraph_vcount(&self->g),
				(long)igraph_ecount(&self->g));
   else
     return PyString_FromFormat("Undirected graph (|V| = %ld, |E| = %ld)",
				(long)igraph_vcount(&self->g),
				(long)igraph_ecount(&self->g));
}

/** \ingroup python_interface
 * \brief Returns the number of vertices in an \c igraph.Graph object.
 * \return the number of vertices as a \c PyObject
 * \sa igraph_vcount
 */
static PyObject* igraphmodule_Graph_vcount(igraphmodule_GraphObject *self) 
{
   PyObject *result;
   result=Py_BuildValue("l", (long)igraph_vcount(&self->g));
   return result;
}

/** \ingroup python_interface
 * \brief Returns the number of edges in an \c igraph.Graph object.
 * \return the number of edges as a \c PyObject
 * \sa igraph_ecount
 */
static PyObject* igraphmodule_Graph_ecount(igraphmodule_GraphObject *self) 
{
   PyObject *result;
   result=Py_BuildValue("l", (long)igraph_ecount(&self->g));
   return result;
}

/** \ingroup python_interface
 * \brief Checks whether an \c igraph.Graph object is directed.
 * \return \c True if the graph is directed, \c False otherwise.
 * \sa igraph_is_directed
 */
static PyObject* igraphmodule_Graph_is_directed(igraphmodule_GraphObject *self) 
{
   if (igraph_is_directed(&self->g))
     Py_RETURN_TRUE;
   else
     Py_RETURN_FALSE;
}

/** \ingroup python_interface
 * \brief Adds vertices to an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * \sa igraph_add_vertices
 */
static PyObject* igraphmodule_Graph_add_vertices(igraphmodule_GraphObject *self,
						 PyObject *args,
						 PyObject *kwds) 
{
   long n;
   
   if (!PyArg_ParseTuple(args, "l", &n)) return NULL;
   if (n<0)
     {
	// Throw an exception
	PyErr_SetString(PyExc_AssertionError, "Number of vertices to be added can't be negative.");
	return NULL;
     }
   if (igraph_add_vertices(&self->g, (integer_t)n)) {
      igraphmodule_handle_igraph_error();
      return NULL;
   }
   
   Py_INCREF(self);
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Removes vertices from an \c igraph.Graph
 * \return the modified \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_delete_vertices
 */
static PyObject* igraphmodule_Graph_delete_vertices(igraphmodule_GraphObject *self,
						    PyObject *args,
						    PyObject *kwds)
{
   PyObject *list;
   vector_t v;
   
   if (!PyArg_ParseTuple(args, "O", &list)) return NULL;
   Py_INCREF(list);
   
   if (igraphmodule_PyList_to_vector_t(list, &v, 1, 0))
     {
	// something bad happened during conversion, release the
	// list reference and return immediately
	Py_DECREF(list);
	return NULL;
     }
   
   // do the hard work :)
   if (igraph_delete_vertices(&self->g, &v)) 
     {
	igraphmodule_handle_igraph_error();
	Py_DECREF(list);
	vector_destroy(&v);
	return NULL;
     }

   vector_destroy(&v);
   
   Py_DECREF(list);
   
   Py_INCREF(self);
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Adds edges to an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_add_edges
 */
static PyObject* igraphmodule_Graph_add_edges(igraphmodule_GraphObject *self,
					      PyObject *args,
					      PyObject *kwds) 
{
   PyObject *list;
   vector_t v;

   if (!PyArg_ParseTuple(args, "O", &list)) return NULL;
   Py_INCREF(list);
   
   if (igraphmodule_PyList_to_vector_t(list, &v, 1, 1))
     {
	// something bad happened during conversion, release the
	// list reference and return immediately
	Py_DECREF(list);
	return NULL;
     }
   
   // do the hard work :)
   if (igraph_add_edges(&self->g, &v)) 
     {
	igraphmodule_handle_igraph_error();
	Py_DECREF(list);
	vector_destroy(&v);
	return NULL;
     }
   
   Py_DECREF(list);
   
   Py_INCREF(self);
   
   vector_destroy(&v);
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Deletes edges from an \c igraph.Graph
 * \return the extended \c igraph.Graph object
 * 
 * \todo Need more error checking on vertex IDs. (igraph fails when an
 * invalid vertex ID is given)
 * \sa igraph_delete_edges
 */
static PyObject* igraphmodule_Graph_delete_edges(igraphmodule_GraphObject *self,
						 PyObject *args,
						 PyObject *kwds) 
{
   PyObject *list;
   vector_t v;

   if (!PyArg_ParseTuple(args, "O", &list)) return NULL;
   Py_INCREF(list);
   
   if (igraphmodule_PyList_to_vector_t(list, &v, 1, 1))
     {
	// something bad happened during conversion, release the
	// list reference and return immediately
	Py_DECREF(list);
	return NULL;
     }
   
   // do the hard work :)
   if (igraph_delete_edges(&self->g, &v)) 
     {
	igraphmodule_handle_igraph_error();
	Py_DECREF(list);
	vector_destroy(&v);
	return NULL;
     }
   
   Py_DECREF(list);
   
   Py_INCREF(self);
   
   vector_destroy(&v);
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief The degree of some vertices in an \c igraph.Graph
 * \return the degree list as a Python object
 * \sa igraph_degree
 */
static PyObject* igraphmodule_Graph_degree(igraphmodule_GraphObject *self,
					   PyObject *args,
					   PyObject *kwds) 
{
   PyObject *list=NULL;
   int dtype=IGRAPH_ALL;
   bool_t input_was_list;
   PyObject *loops = Py_False;
   vector_t vids, result;
   
   static char *kwlist[] = 
     {
	"vertices", "type", "loops", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OiO!", kwlist,
				    &list, &dtype, &PyBool_Type, &loops))
	return NULL;
   
   if (dtype!=IGRAPH_ALL && dtype!=IGRAPH_OUT && dtype!=IGRAPH_IN) 
     {
	PyErr_SetString(PyExc_ValueError, "dtype should be either ALL or IN or OUT");
	return NULL;
     }
   
   Py_INCREF(loops);
   if (list)
     {
	Py_INCREF(list);
	if (igraphmodule_PyList_to_vector_t(list, &vids, 1, 0)) 
	  {
	     // something bad happened during conversion, release the
	     // list reference and return immediately
	     Py_DECREF(list);
	     Py_DECREF(loops);
	     return NULL;
	  }
	input_was_list=PyList_Check(list);
     }
   else
     {
	// no list was given, so we assume that the user wanted to
	// retrieve the degrees of all vertices
	if (igraph_vcount(&self->g)>0)
	  vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1);
	else
	  vector_init(&vids, 0);
	input_was_list=1;
     }

   vector_init(&result, 0);
   if (igraph_degree(&self->g, &result, &vids,
		     (igraph_neimode_t)dtype, (bool_t)(loops == Py_True))) 
     {
	igraphmodule_handle_igraph_error();
	Py_DECREF(list);
	Py_DECREF(loops);
	vector_destroy(&result);
	return NULL;
     }
   
   Py_DECREF(loops);
   if (list) Py_DECREF(list);
   
   if (input_was_list) 
     list=igraphmodule_vector_t_to_PyList(&result);
   else
     list=PyInt_FromLong(VECTOR(result)[0]);
   
   vector_destroy(&result);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief The neighbors of a given vertex in an \c igraph.Graph
 * This method accepts a single vertex ID as a parameter, and returns the
 * neighbors of the given vertex in the form of an integer list. A
 * second argument may be passed as well, meaning the type of neighbors to
 * be returned (\c OUT for successors, \c IN for predecessors or \c ALL
 * for both of them). This argument is ignored for undirected graphs.
 * 
 * \return the neighbor list as a Python list object
 * \sa igraph_neighbors
 */
static PyObject* igraphmodule_Graph_neighbors(igraphmodule_GraphObject *self,
					   PyObject *args,
					   PyObject *kwds) 
{
   PyObject *list;
   int dtype=IGRAPH_ALL;
   long idx;
   vector_t result;
   
   static char *kwlist[] = 
     {
	"vertex", "type", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|i", kwlist,
				    &idx, &dtype))
     return NULL;
   
   if (dtype!=IGRAPH_ALL && dtype!=IGRAPH_OUT && dtype!=IGRAPH_IN) 
     {
	PyErr_SetString(PyExc_ValueError, "dtype should be either ALL or IN or OUT");
	return NULL;
     }
   
   vector_init(&result, 1);
   if (igraph_neighbors(&self->g, &result, idx, (igraph_neimode_t)dtype))
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&result);
	return NULL;
     }
   
   list=igraphmodule_vector_t_to_PyList(&result);
   vector_destroy(&result);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the diameter of an \c igraph.Graph
 * This method accepts two optional parameters: the first one is
 * a boolean meaning whether to consider directed paths (and is
 * ignored for undirected graphs). The second one is only meaningful
 * in unconnected graphs: it is \c True if the longest geodesic
 * within a component should be returned and \c False if the number of
 * vertices should be returned. They both have a default value of \c False.
 * 
 * \return the diameter as a Python integer
 * \sa igraph_diameter
 */
static PyObject* igraphmodule_Graph_diameter(igraphmodule_GraphObject *self,
					     PyObject *args,
					     PyObject *kwds) 
{
  PyObject *dir=NULL, *vcount_if_unconnected=NULL;
  integer_t i;
  int r;
   
  static char *kwlist[] = {
    "directed", "unconn", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!", kwlist,
				   &PyBool_Type, &dir,
				   &PyBool_Type, &vcount_if_unconnected))
    return NULL;
  
  Py_BEGIN_ALLOW_THREADS
  r=igraph_diameter(&self->g, &i, (bool_t)(dir == Py_True),
		    (bool_t)(vcount_if_unconnected == Py_True));
  Py_END_ALLOW_THREADS
  if (r) 
     {
	igraphmodule_handle_igraph_error();
	return NULL;
     }
   
   return PyInt_FromLong((long)i);
}

/** \ingroup python_interface
 * \brief Generates a graph based on the Barabási-Albert model
 * This is intended to be a class method in Python, so the first argument
 * is the type object and not the Python igraph object (because we have
 * to allocate that in this method).
 * 
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_barabasi_game
 */
static PyObject* igraphmodule_Graph_Barabasi(PyTypeObject *type,
					     PyObject *args,
					     PyObject *kwds) 
{
   igraphmodule_GraphObject *self;
   long n, m;
   vector_t outseq;
   PyObject *m_obj, *outpref=NULL, *directed=NULL;
   
   static char *kwlist[] = 
     {
	"n", "m", "outpref", "directed", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "lO|O!O!", kwlist,
				    &n, &m_obj,
				    &PyBool_Type, &outpref,
				    &PyBool_Type, &directed))
     return NULL;
   
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   // let's check whether we have a constant out-degree or a list
   if (PyInt_Check(m_obj)) 
     {
	m=PyInt_AsLong(m_obj);
	vector_init(&outseq, 0);
     }
   else if (PyList_Check(m_obj))
     {
	if (igraphmodule_PyList_to_vector_t(m_obj, &outseq, 1, 0)) 
	  {
	     // something bad happened during conversion
	     return NULL;
	  }
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_barabasi_game(&self->g, (integer_t)n, (integer_t)m,
				  &outseq, (outpref == Py_True),
				  (directed == Py_True)))
	  {
	     igraphmodule_handle_igraph_error();
	     vector_destroy(&outseq);
	     return NULL;
	  }
     }
   
   vector_destroy(&outseq);
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a graph based on the Erdõs-Rényi model
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_erdos_renyi_game
 */
static PyObject* igraphmodule_Graph_Erdos_Renyi(PyTypeObject *type,
						PyObject *args,
						PyObject *kwds) 
{
   igraphmodule_GraphObject *self;
   long n, m=-1;
   double p=-1.0;
   igraph_erdos_renyi_t t;
   PyObject *loops=NULL, *directed=NULL;
   
   static char *kwlist[] =
     {
	"n", "p", "m", "directed", "loops", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|dlO!O!", kwlist,
				    &n, &p, &m,
				    &PyBool_Type, &directed,
				    &PyBool_Type, &loops))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   if (m==-1 && p==-1.0)
     {
	// no density parameters were given, throw exception
	PyErr_SetString(PyExc_TypeError, "Either m or p must be given.");
	return NULL;
     }
   if (m!=-1 && p!=-1.0)
     {
	// both density parameters were given, throw exception
	PyErr_SetString(PyExc_TypeError, "Only one must be given from m and p.");
	return NULL;
     }
   
   t=(m==-1)?IGRAPH_ERDOS_RENYI_GNP:IGRAPH_ERDOS_RENYI_GNM;
   
   if (t==IGRAPH_ERDOS_RENYI_GNP) 
     {
	if (p<0.0 || p>1.0)
	  {
	     // Invalid probability was given, throw exception
	     PyErr_SetString(PyExc_ValueError, "p must be between 0 and 1.");
	     return NULL;
	  }	
     }
   else
     {
	if (m<0 || m>n*n)
	  {
	     // Invalid edge count was given, throw exception
	     PyErr_SetString(PyExc_ValueError, "m must be between 0 and n^2.");
	     return NULL;
	  }	
     }
      
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_erdos_renyi_game(&self->g, t, (integer_t)n,
				    (real_t)((t==IGRAPH_ERDOS_RENYI_GNM)?m:p),
				    (directed == Py_True),
				    (loops == Py_True)))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a full graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_full
 */
static PyObject* igraphmodule_Graph_Full(PyTypeObject *type,
					 PyObject *args,
					 PyObject *kwds) 
{
   igraphmodule_GraphObject *self;
   long n;
   PyObject *loops=NULL, *directed=NULL;
   
   static char *kwlist[] =
     {
	"n", "directed", "loops", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|O!O!", kwlist, &n,
				    &PyBool_Type, &directed,
				    &PyBool_Type, &loops))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_full(&self->g, (integer_t)n,
			(directed == Py_True), (loops == Py_True)))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a growing random graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_growing_random_game
 */
static PyObject* igraphmodule_Graph_Growing_Random(PyTypeObject *type,
						   PyObject *args,
						   PyObject *kwds) 
{
   long n, m;
   PyObject *directed=NULL, *citation=NULL;
   igraphmodule_GraphObject *self;
   
   static char *kwlist[] =
     {
	"n", "m", "directed", "citation", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|O!O!", kwlist, &n, &m,
				    &PyBool_Type, &directed,
				    &PyBool_Type, &citation))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   if (m<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of new edges per iteration must be positive.");
	return NULL;
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_growing_random_game(&self->g, (integer_t)n,
				       (integer_t)m, (directed == Py_True),
				       (citation == Py_True)))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a star graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_star
 */
static PyObject* igraphmodule_Graph_Star(PyTypeObject *type,
					 PyObject *args,
					 PyObject *kwds) 
{
   long n, center=0;
   igraph_star_mode_t mode=IGRAPH_STAR_UNDIRECTED;
   igraphmodule_GraphObject *self;
   
   static char *kwlist[] =
     {
	"n", "mode", "center", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|ll", kwlist,
				    &n, &mode, &center))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   if (center>=n || center<0)
     {
	PyErr_SetString(PyExc_ValueError, "Central vertex ID should be between 0 and n-1");
	return NULL;
     }
   
   if (mode!=IGRAPH_STAR_UNDIRECTED && mode!=IGRAPH_STAR_IN &&
       mode!=IGRAPH_STAR_OUT) 
     {
	PyErr_SetString(PyExc_ValueError, "Mode should be either STAR_IN, STAR_OUT or STAR_UNDIRECTED.");
	return NULL;
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_star(&self->g, (integer_t)n, mode, (integer_t)center))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a ring-shaped graph
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_ring
 */
static PyObject* igraphmodule_Graph_Ring(PyTypeObject *type,
					 PyObject *args,
					 PyObject *kwds) 
{
   long n, m;
   PyObject *directed=Py_False, *mutual=Py_False, *circular=Py_True;
   igraphmodule_GraphObject *self;
   
   static char *kwlist[] =
     {
	"n", "directed", "mutual", "circular", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|O!O!O!", kwlist, &n,
				    &PyBool_Type, &directed,
				    &PyBool_Type, &mutual,
				    &PyBool_Type, &circular))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_ring(&self->g, (integer_t)n, (directed == Py_True),
			(mutual == Py_True), (circular == Py_True)))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Generates a tree graph where almost all vertices have an equal number of children
 * \return a reference to the newly generated Python igraph object
 * \sa igraph_tree
 */
static PyObject* igraphmodule_Graph_Tree(PyTypeObject *type,
					 PyObject *args,
					 PyObject *kwds) 
{
   long n, children;
   igraph_tree_mode_t mode=IGRAPH_TREE_UNDIRECTED;
   igraphmodule_GraphObject *self;
   
   static char *kwlist[] =
     {
	"n", "children", "type", NULL
     }
   ;

   if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll|l", kwlist,
				    &n, &children, &mode))
     return NULL;
      
   if (n<0)
     {
	PyErr_SetString(PyExc_ValueError, "Number of vertices must be positive.");
	return NULL;
     }
   
   if (mode!=IGRAPH_TREE_UNDIRECTED && mode!=IGRAPH_TREE_IN &&
       mode!=IGRAPH_TREE_OUT) 
     {
	PyErr_SetString(PyExc_ValueError, "Mode should be either TREE_IN, TREE_OUT or TREE_UNDIRECTED.");
	return NULL;
     }
   
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
	if (igraph_tree(&self->g, (integer_t)n, (integer_t)children, mode))
	  {
	     igraphmodule_handle_igraph_error();
	     return NULL;
	  }
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Decides whether a graph is connected.
 * \return Py_True if the graph is connected, Py_False otherwise
 * \sa igraph_is_connected
 */
static PyObject* igraphmodule_Graph_is_connected(igraphmodule_GraphObject *self,
						 PyObject *args,
						 PyObject *kwds) 
{
   static char *kwlist[] = {"mode", NULL};
   igraph_connectedness_t mode=IGRAPH_STRONG;
   bool_t res;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &mode))
     return NULL;

   if (mode != IGRAPH_STRONG && mode != IGRAPH_WEAK) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be either STRONG or WEAK");
	return NULL;
     }
   
   if (igraph_is_connected(&self->g, &res, mode)) 
     {
	igraphmodule_handle_igraph_error();
	return NULL;
     }
   if (res) Py_RETURN_TRUE; else Py_RETURN_FALSE;
}

/** \ingroup python_interface
 * \brief Decides whether there is an edge from a given vertex to an other one.
 * \return Py_True if the vertices are directly connected, Py_False otherwise
 * \sa igraph_are_connected
 */
static PyObject* igraphmodule_Graph_are_connected(igraphmodule_GraphObject *self,
						  PyObject *args,
						  PyObject *kwds) 
{
   static char *kwlist[] = {"v1", "v2", NULL};
   long v1, v2;
   bool_t res;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "ll", kwlist, &v1, &v2))
     return NULL;

   res=igraph_are_connected(&self->g, (integer_t)v1, (integer_t)v2);
   if (res) Py_RETURN_TRUE; else Py_RETURN_FALSE;
}

/** \ingroup python_interface
 * \brief Calculates the average path length in a graph.
 * \return the average path length as a PyObject
 * \sa igraph_average_path_length
 */
static PyObject* igraphmodule_Graph_average_path_length(igraphmodule_GraphObject *self,
							PyObject *args,
							PyObject *kwds) 
{
   static char *kwlist[] = {"directed", "unconn", NULL};
   long v1, v2;
   PyObject *directed=Py_True, *unconn=Py_True;
   real_t res;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!O!", kwlist,
				    &PyBool_Type, &directed,
				    &PyBool_Type, &unconn))
     return NULL;

   if (igraph_average_path_length(&self->g, &res, (directed==Py_True),
				  (unconn==Py_True))) 
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   return PyFloat_FromDouble(res);
}

/** \ingroup python_interface
 * \brief Calculates the betweennesses of some nodes in a graph.
 * \return the betweennesses as a list (or a single float)
 * \sa igraph_betweenness
 */
static PyObject* igraphmodule_Graph_betweenness(igraphmodule_GraphObject *self,
						PyObject *args,
						PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", "directed", NULL};
   PyObject *directed=Py_True;
   PyObject *vobj=NULL, *list=NULL;
   vector_t vids;
   vector_t res;
   int return_single=0;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,
				    &vobj, &directed))
     return NULL;

   if (vobj == NULL) 
     {
	// no vertex list was supplied
	if (igraph_vcount(&self->g)>0) 
	  {
	     if (vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1)) 
	       return igraphmodule_handle_igraph_error();
	  }
	else
	  {
	     if (vector_init(&vids, 0)) return igraphmodule_handle_igraph_error();
	  }
     }
   else
     {
	if (PyInt_Check(vobj)) return_single=1;
	
	Py_INCREF(vobj);
	// vertex list was specified, convert to vector_t
	if (igraphmodule_PyList_to_vector_t(vobj, &vids, 1, 0)) {
	   Py_DECREF(vobj);
	   return NULL;
	}
	Py_DECREF(vobj);
     }

   
   if (vector_init(&res, vector_size(&vids))) return igraphmodule_handle_igraph_error();
   
   if (igraph_betweenness(&self->g, &res, &vids, PyObject_IsTrue(directed)))
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   if (!return_single)
     list=igraphmodule_vector_t_to_float_PyList(&res);
   else
     list=PyFloat_FromDouble(VECTOR(res)[0]);
   
   vector_destroy(&res);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the bibliographic coupling of some nodes in a graph.
 * \return the bibliographic coupling values in a matrix
 * \sa igraph_bibcoupling
 */
static PyObject* igraphmodule_Graph_bibcoupling(igraphmodule_GraphObject *self,
						PyObject *args,
						PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", NULL};
   PyObject *vobj=NULL, *list=NULL;
   vector_t vids;
   matrix_t res;
   int return_single=0;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &vobj))
     return NULL;

   if (vobj == NULL) 
     {
	// no vertex list was supplied
	if (igraph_vcount(&self->g)>0) 
	  {
	     if (vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1)) 
	       return igraphmodule_handle_igraph_error();
	  }
	else
	  {
	     if (vector_init(&vids, 0)) return igraphmodule_handle_igraph_error();
	  }
     }
   else
     {
	if (PyInt_Check(vobj)) return_single=1;
	
	Py_INCREF(vobj);
	// vertex list was specified, convert to vector_t
	if (igraphmodule_PyList_to_vector_t(vobj, &vids, 1, 0)) {
	   Py_DECREF(vobj);
	   return NULL;
	}
	Py_DECREF(vobj);
     }

   
   if (matrix_init(&res, vector_size(&vids), igraph_vcount(&self->g)))
     return igraphmodule_handle_igraph_error();
   
   if (igraph_bibcoupling(&self->g, &res, &vids))
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }

   /// \todo Return a single list instead of a matrix if only one vertex was given
   list=igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
   
   matrix_destroy(&res);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the closeness centrality of some nodes in a graph.
 * \return the closeness centralities as a list (or a single float)
 * \sa igraph_betweenness
 */
static PyObject* igraphmodule_Graph_closeness(igraphmodule_GraphObject *self,
					      PyObject *args,
					      PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", "mode", NULL};
   PyObject *vobj=NULL, *list=NULL;
   vector_t vids;
   vector_t res;
   igraph_neimode_t mode=IGRAPH_ALL;
   int return_single=0;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Ol", kwlist,
				    &vobj, &mode))
     return NULL;

   if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be one of IN, OUT or ALL");
	return NULL;
     }
   
   if (vobj == NULL) 
     {
	// no vertex list was supplied
	if (igraph_vcount(&self->g)>0) 
	  {
	     if (vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1)) 
	       return igraphmodule_handle_igraph_error();
	  }
	else
	  {
	     if (vector_init(&vids, 0)) return igraphmodule_handle_igraph_error();
	  }
     }
   else
     {
	if (PyInt_Check(vobj)) return_single=1;
	
	Py_INCREF(vobj);
	// vertex list was specified, convert to vector_t
	if (igraphmodule_PyList_to_vector_t(vobj, &vids, 1, 0)) {
	   Py_DECREF(vobj);
	   return NULL;
	}
	Py_DECREF(vobj);
     }

   
   if (vector_init(&res, vector_size(&vids))) return igraphmodule_handle_igraph_error();
   
   if (igraph_closeness(&self->g, &res, &vids, mode))
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   if (!return_single)
     list=igraphmodule_vector_t_to_float_PyList(&res);
   else
     list=PyFloat_FromDouble(VECTOR(res)[0]);
   
   vector_destroy(&res);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the (weakly or strongly) connected components in a graph.
 * \return a list containing the cluster ID for every vertex in the graph
 * \sa igraph_clusters
 */
static PyObject* igraphmodule_Graph_clusters(igraphmodule_GraphObject *self,
					     PyObject *args,
					     PyObject *kwds) 
{
   static char *kwlist[] = {"mode", NULL};
   igraph_connectedness_t mode=IGRAPH_STRONG;
   vector_t res1, res2;
   PyObject *list;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|l", kwlist, &mode))
     return NULL;

   if (mode != IGRAPH_STRONG && mode != IGRAPH_WEAK) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be either STRONG or WEAK");
	return NULL;
     }
   
   vector_init(&res1, igraph_vcount(&self->g));
   vector_init(&res2, 10);
   
   if (igraph_clusters(&self->g, &res1, &res2, mode)) 
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&res1);
	vector_destroy(&res2);
	return NULL;
     }
   
   list=igraphmodule_vector_t_to_PyList(&res1);
   vector_destroy(&res1);
   vector_destroy(&res2);
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the cocitation scores of some nodes in a graph.
 * \return the cocitation scores in a matrix
 * \sa igraph_cocitation
 */
static PyObject* igraphmodule_Graph_cocitation(igraphmodule_GraphObject *self,
					       PyObject *args,
					       PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", NULL};
   PyObject *vobj=NULL, *list=NULL;
   vector_t vids;
   matrix_t res;
   int return_single=0;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &vobj))
     return NULL;

   if (vobj == NULL) 
     {
	// no vertex list was supplied
	if (igraph_vcount(&self->g)>0) 
	  {
	     if (vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1)) 
	       return igraphmodule_handle_igraph_error();
	  }
	else
	  {
	     if (vector_init(&vids, 0)) return igraphmodule_handle_igraph_error();
	  }
     }
   else
     {
	if (PyInt_Check(vobj)) return_single=1;
	
	Py_INCREF(vobj);
	// vertex list was specified, convert to vector_t
	if (igraphmodule_PyList_to_vector_t(vobj, &vids, 1, 0)) {
	   Py_DECREF(vobj);
	   return NULL;
	}
	Py_DECREF(vobj);
     }

   
   if (matrix_init(&res, vector_size(&vids), igraph_vcount(&self->g)))
     return igraphmodule_handle_igraph_error();
   
   if (igraph_cocitation(&self->g, &res, &vids))
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }

   /// \todo Return a single list instead of a matrix if only one vertex was given
   list=igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
   
   matrix_destroy(&res);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the edge betweennesses in the graph
 * \return a list containing the edge betweenness for every edge
 * \sa igraph_edge_betweenness
 */
static PyObject* igraphmodule_Graph_edge_betweenness(igraphmodule_GraphObject *self,
						     PyObject *args,
						     PyObject *kwds) 
{
   static char *kwlist[] = {"directed", NULL};
   vector_t res;
   PyObject *list, *directed=Py_True;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist,
				    &PyBool_Type, &directed))
     return NULL;

   vector_init(&res, igraph_ecount(&self->g));
   
   if (igraph_edge_betweenness(&self->g, &res, (directed==Py_True)))
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&res);
	return NULL;
     }
   
   list=igraphmodule_vector_t_to_float_PyList(&res);
   vector_destroy(&res);
   return list;
}

/** \ingroup python_interface
 * \brief Calculates the shortest paths from/to a given node in the graph
 * \return a list containing shortest paths from/to the given node
 * \sa igraph_get_shortest_paths
 */
static PyObject* igraphmodule_Graph_get_shortest_paths(igraphmodule_GraphObject *self,
						       PyObject *args,
						       PyObject *kwds) 
{
   static char *kwlist[] = {"v", "mode", NULL};
   vector_t *res;
   igraph_neimode_t mode=IGRAPH_ALL;
   long from0, i, j;
   integer_t from;
   PyObject *list, *item;
   long int no_of_nodes=igraph_vcount(&self->g);
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|l", kwlist,
				    &from0, &mode))
     return NULL;

   if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be either IN, OUT or ALL");
	return NULL;
     }
   
   if (from0<0 || from0>=igraph_vcount(&self->g)) 
     {
	PyErr_SetString(PyExc_ValueError, "vertex ID must be non-negative and less than the number of edges");
	return NULL;
     }
   from=(integer_t)from0;
   
   res=(vector_t*)calloc(no_of_nodes, sizeof(vector_t));
   if (!res) 
     {
	PyErr_SetString(PyExc_MemoryError, "");
	return NULL;
     }
   
   for (i=0; i<no_of_nodes; i++) vector_init(&res[i], 5);
   
   if (igraph_get_shortest_paths(&self->g, res, from, mode))
     {
	igraphmodule_handle_igraph_error();
	for (j=0; j<no_of_nodes; j++) vector_destroy(&res[j]);
	free(res);
	return NULL;
     }

   list=PyList_New(no_of_nodes);
   if (!list) {
      for (j=0; j<no_of_nodes; j++) vector_destroy(&res[j]);
      free(res);
      return NULL;
   }
   
   for (i=0; i<no_of_nodes; i++) 
     {
	item=igraphmodule_vector_t_to_PyList(&res[i]);
	if (!item) 
	  {
	     Py_DECREF(list);
	     for (j=0; j<no_of_nodes; j++) vector_destroy(&res[j]);
	     free(res);
	     return NULL;
	  }
	if (PyList_SetItem(list, i, item)) 
	  {
	     Py_DECREF(list);
	     for (j=0; j<no_of_nodes; j++) vector_destroy(&res[j]);
	     free(res);
	     return NULL;
	  }
     }
   
   for (j=0; j<no_of_nodes; j++) vector_destroy(&res[j]);
   free(res);
   return list;
}

/** \ingroup python_interface
 * \brief Calculates shortest paths in a graph.
 * \return the shortest path lengths for the given vertices
 * \sa igraph_shortest_paths
 */
static PyObject* igraphmodule_Graph_shortest_paths(igraphmodule_GraphObject *self,
						   PyObject *args,
						   PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", "mode", NULL};
   PyObject *vobj=NULL, *list=NULL;
   vector_t vids;
   matrix_t res;
   igraph_neimode_t mode=IGRAPH_ALL;
   int return_single=0;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|Ol", kwlist, &vobj, &mode))
     return NULL;

   if (mode!=IGRAPH_IN && mode!=IGRAPH_OUT && mode!=IGRAPH_ALL) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be either IN or OUT or ALL");
	return NULL;
     }
   
   if (vobj == NULL) 
     {
	// no vertex list was supplied
	if (igraph_vcount(&self->g)>0) 
	  {
	     if (vector_init_seq(&vids, 0, igraph_vcount(&self->g)-1)) 
	       return igraphmodule_handle_igraph_error();
	  }
	else
	  {
	     if (vector_init(&vids, 0)) return igraphmodule_handle_igraph_error();
	  }
     }
   else
     {
	if (PyInt_Check(vobj)) return_single=1;
	
	Py_INCREF(vobj);
	// vertex list was specified, convert to vector_t
	if (igraphmodule_PyList_to_vector_t(vobj, &vids, 1, 0)) {
	   Py_DECREF(vobj);
	   return NULL;
	}
	Py_DECREF(vobj);
     }

   
   if (matrix_init(&res, vector_size(&vids), igraph_vcount(&self->g)))
     return igraphmodule_handle_igraph_error();
   
   if (igraph_shortest_paths(&self->g, &res, &vids, mode))
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }

   /// \todo Return a single list instead of a matrix if only one vertex was given
   list=igraphmodule_matrix_t_to_PyList(&res, IGRAPHMODULE_TYPE_INT);
   
   matrix_destroy(&res);
   vector_destroy(&vids);
   
   return list;
}

/** \ingroup python_interface
 * \brief Calculates a spanning tree for a graph
 * \return a list containing the edge betweenness for every edge
 * \sa igraph_minimum_spanning_tree_unweighted
 * \sa igraph_minimum_spanning_tree_prim
 */
static PyObject* igraphmodule_Graph_spanning_tree(igraphmodule_GraphObject *self,
						  PyObject *args,
						  PyObject *kwds) 
{
   static char *kwlist[] = {"weights", NULL};
   igraph_t mst;
   int err;
   vector_t ws;
   PyObject *weights=NULL;
   igraphmodule_GraphObject *result=NULL;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O!", kwlist,
				    &PyList_Type, &weights))
     return NULL;

   if (weights && (PyList_Size(weights)<igraph_vcount(&self->g))) 
     {
	PyErr_SetString(PyExc_ValueError, "Weight list must have at least |V| elements (|V| = node count in the graph)");
	return NULL;
     }

   if (!weights)
     err=igraph_minimum_spanning_tree_unweighted(&self->g, &mst);
   else 
     {
	if (igraphmodule_PyList_to_vector_t(weights, &ws, 1, 0)) return NULL;
	err=igraph_minimum_spanning_tree_prim(&self->g, &mst, &ws);
     }
   
   if (err)
     {
	igraphmodule_handle_igraph_error();
	if (weights) vector_destroy(&ws);
	return NULL;
     }

   result = (igraphmodule_GraphObject*)(self->ob_type->tp_alloc(self->ob_type, 0));
   if (result != NULL) result->g=mst;
   
   if (weights) vector_destroy(&ws);
   
   return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Simplifies a graph by removing loops and/or multiple edges
 * \return the simplified graph.
 * \sa igraph_simplify
 */
static PyObject* igraphmodule_Graph_simplify(igraphmodule_GraphObject *self,
					     PyObject *args,
					     PyObject *kwds) 
{
   static char *kwlist[] = {"multiple", "loops", NULL};
   PyObject *multiple=Py_True, *loops=Py_True;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|OO", kwlist,
				    &multiple, &loops))
     return NULL;

   if (igraph_simplify(&self->g, PyObject_IsTrue(multiple),
		       PyObject_IsTrue(loops)))
     {
	igraphmodule_handle_igraph_error();
	return NULL;
     }

   Py_INCREF(self);
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Calculates the vertex indices within the same component as a given vertex
 * \return the vertex indices in a list
 * \sa igraph_subcomponent
 */
static PyObject* igraphmodule_Graph_subcomponent(igraphmodule_GraphObject *self,
					     PyObject *args,
					     PyObject *kwds) 
{
   static char *kwlist[] = {"v", "mode", NULL};
   vector_t res;
   igraph_neimode_t mode=IGRAPH_ALL;
   long from0;
   real_t from;
   PyObject *multiple=Py_True, *loops=Py_True, *list=NULL;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "l|l", kwlist,
				    &from0, &mode))
     return NULL;

   if (mode != IGRAPH_OUT && mode != IGRAPH_IN && mode != IGRAPH_ALL) 
     {
	PyErr_SetString(PyExc_ValueError, "mode must be either IN, OUT or ALL");
	return NULL;
     }
   
   if (from0<0 || from0>=igraph_vcount(&self->g)) 
     {
	PyErr_SetString(PyExc_ValueError, "vertex ID must be non-negative and less than the number of edges");
	return NULL;
     }
   from=(real_t)from0;

   vector_init(&res, 0);
   if (igraph_subcomponent(&self->g, &res, from, mode))
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&res);
	return NULL;
     }

   list=igraphmodule_vector_t_to_PyList(&res);
   vector_destroy(&res);
   
   return list;
}

/** \ingroup python_interface
 * \brief Returns a subgraph of the graph based on the given vertices
 * \return the subgraph as a new igraph object
 * \sa igraph_subgraph
 */
static PyObject* igraphmodule_Graph_subgraph(igraphmodule_GraphObject *self,
					     PyObject *args,
					     PyObject *kwds) 
{
   static char *kwlist[] = {"vertices", NULL};
   vector_t vertices;
   igraph_t sg;
   igraphmodule_GraphObject *result;
   PyObject *list;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &list))
     return NULL;

   if (igraphmodule_PyList_to_vector_t(list, &vertices, 1, 0))
     return NULL;
   
   if (igraph_subgraph(&self->g, &sg, &vertices))
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&vertices);
	return NULL;
     }

   result = (igraphmodule_GraphObject*)(self->ob_type->tp_alloc(self->ob_type, 0));
   if (result != NULL) result->g=sg;
   
   vector_destroy(&vertices);
   
   return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Places the vertices of a graph uniformly on a circle.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_circle
 */
static PyObject* igraphmodule_Graph_layout_circle(igraphmodule_GraphObject *self,
						  PyObject *args,
						  PyObject *kwds) 
{
   matrix_t m;
   PyObject *result;
   
   if (matrix_init(&m, 1, 1)) 
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   if (igraph_layout_circle(&self->g, &m))
     {
	matrix_destroy(&m);
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   result=igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
   
   matrix_destroy(&m);
   
   return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Places the vertices of a graph randomly.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_circle
 */
static PyObject* igraphmodule_Graph_layout_random(igraphmodule_GraphObject *self,
						  PyObject *args,
						  PyObject *kwds) 
{
   matrix_t m;
   PyObject *result;
   
   if (matrix_init(&m, 1, 1)) 
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   if (igraph_layout_random(&self->g, &m))
     {
	matrix_destroy(&m);
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   result=igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
   
   matrix_destroy(&m);
   
   return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Places the vertices on a plane according to the Kamada-Kawai algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_kamada_kawai
 */
static PyObject* igraphmodule_Graph_layout_kamada_kawai(igraphmodule_GraphObject *self,
							PyObject *args,
							PyObject *kwds) 
{
  static char *kwlist[] = {"n", "sigma", "initemp", "coolexp", "kkconst", NULL};
  matrix_t m;
  long niter=1000;
  double sigma, initemp, coolexp, kkconst;
  PyObject *result;
   
  sigma=igraph_vcount(&self->g);
  kkconst=sigma*sigma; sigma=sigma/4.0;
  initemp=10.0; coolexp=0.99;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ldddd", kwlist,
				   &niter, &sigma, &initemp, &coolexp, &kkconst))
    return NULL;
  
  if (matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
   
  if (igraph_layout_kamada_kawai(&self->g, &m, niter, sigma, initemp, coolexp, kkconst)) {
    matrix_destroy(&m);
    igraphmodule_handle_igraph_error(); return NULL;
  }
   
  result=igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
   
  matrix_destroy(&m);
   
  return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Places the vertices on a plane according to the Fruchterman-Reingold algorithm.
 * \return the calculated coordinates as a Python list of lists
 * \sa igraph_layout_kamada_kawai
 */
static PyObject* igraphmodule_Graph_layout_fruchterman_reingold(igraphmodule_GraphObject *self,
								PyObject *args,
								PyObject *kwds) 
{
  static char *kwlist[] = {"n", "maxdelta", "area", "coolexp", "repulserad", NULL};
  matrix_t m;
  long niter=500;
  double maxdelta, area, coolexp, repulserad;
  PyObject *result;
   
  maxdelta=igraph_vcount(&self->g);
  area=maxdelta*maxdelta; coolexp=1.5;
  repulserad=area*maxdelta;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|ldddd", kwlist,
				   &niter, &maxdelta, &area, &coolexp, &repulserad))
    return NULL;
  
  if (matrix_init(&m, 1, 1)) {
    igraphmodule_handle_igraph_error(); return NULL;
  }
   
  if (igraph_layout_fruchterman_reingold(&self->g, &m, niter, maxdelta, area, coolexp, repulserad, 0)) {
    matrix_destroy(&m);
    igraphmodule_handle_igraph_error(); return NULL;
  }
   
  result=igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_FLOAT);
   
  matrix_destroy(&m);
   
  return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Returns the adjacency matrix of a graph.
 * \return the adjacency matrix as a Python list of lists
 * \sa igraph_get_adjacency
 */
static PyObject* igraphmodule_Graph_get_adjacency(igraphmodule_GraphObject *self,
						  PyObject *args,
						  PyObject *kwds) 
{
   static char *kwlist[] = {"type", NULL};
   igraph_get_adjacency_t t=IGRAPH_GET_ADJACENCY_BOTH;
   matrix_t m;
   PyObject *result;
   
   if (!PyArg_ParseTupleAndKeywords(args, kwds, "|i", kwlist, &t)) return NULL;
   
   if (t!=IGRAPH_GET_ADJACENCY_UPPER && t!=IGRAPH_GET_ADJACENCY_LOWER &&
       t!=IGRAPH_GET_ADJACENCY_BOTH)
     {
	PyErr_SetString(PyExc_ValueError, "type must be either GET_ADJACENCY_LOWER or GET_ADJACENCY_UPPER or GET_ADJACENCY_BOTH");
	return NULL;
     }
   
   if (matrix_init(&m, igraph_vcount(&self->g), igraph_vcount(&self->g))) 
     {
	igraphmodule_handle_igraph_error(); return NULL;
     }
   
   if (igraph_get_adjacency(&self->g, &m, t)) 
     {
	igraphmodule_handle_igraph_error();
	matrix_destroy(&m);
	return NULL;
     }
   
   result=igraphmodule_matrix_t_to_PyList(&m, IGRAPHMODULE_TYPE_INT);
   matrix_destroy(&m);
   return result;
}

/** \ingroup python_interface
 * \brief Returns the list of edges in a graph.
 * \return the list of edges, every edge is represented by a pair
 * \sa igraph_get_edgelist
 */
static PyObject* igraphmodule_Graph_get_edgelist(igraphmodule_GraphObject *self,
						 PyObject *args,
						 PyObject *kwds) 
{
   vector_t edgelist;
   PyObject *result;
   
   vector_init(&edgelist, igraph_ecount(&self->g));
   if (igraph_get_edgelist(&self->g, &edgelist, 0))
     {
	igraphmodule_handle_igraph_error();
	vector_destroy(&edgelist);
	return NULL;
     }
   
   result=igraphmodule_vector_t_to_PyList_pairs(&edgelist);
   vector_destroy(&edgelist);
   
   return (PyObject*)result;
}

/** \ingroup python_interface
 * \brief Reads an edge list from a file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_edgelist
 */
static PyObject* igraphmodule_Graph_Read_Edgelist(PyTypeObject *type,
						  PyObject *args,
						  PyObject *kwds)
{
  igraphmodule_GraphObject *self;
  char* fname=NULL;
  FILE* f;
  PyObject *directed=Py_True;
  
  static char *kwlist[] =
  {
    "f", "directed", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|O", kwlist,
				   &fname, &directed))
     return NULL;
      
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
       f=fopen(fname, "r");
       if (!f) {
	 PyErr_SetString(PyExc_IOError, strerror(errno));
	 return NULL;
       }
       if (igraph_read_graph_edgelist(&self->g, f, 0, PyObject_IsTrue(directed)))
       {
	 igraphmodule_handle_igraph_error();
	 fclose(f);
	 return NULL;
       }
       fclose(f);
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Reads an edge list from an NCOL file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_ncol
 */
static PyObject* igraphmodule_Graph_Read_Ncol(PyTypeObject *type,
					      PyObject *args,
					      PyObject *kwds)
{
  igraphmodule_GraphObject *self;
  char* fname=NULL;
  FILE* f;
  PyObject *names=Py_True, *weights=Py_True;
  
  static char *kwlist[] =
  {
    "f", "names", "weights", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|OO", kwlist,
				   &fname, &names, &weights))
     return NULL;
      
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
       f=fopen(fname, "r");
       if (!f) {
	 PyErr_SetString(PyExc_IOError, strerror(errno));
	 return NULL;
       }
       if (igraph_read_graph_ncol(&self->g, f, PyObject_IsTrue(names), PyObject_IsTrue(weights)))
       {
	 igraphmodule_handle_igraph_error();
	 fclose(f);
	 return NULL;
       }
       fclose(f);
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Reads an edge list from an LGL file and creates a graph from it.
 * \return the graph
 * \sa igraph_read_graph_lgl
 */
static PyObject* igraphmodule_Graph_Read_Lgl(PyTypeObject *type,
					     PyObject *args,
					     PyObject *kwds)
{
  igraphmodule_GraphObject *self;
  char* fname=NULL;
  FILE* f;
  PyObject *names=Py_True, *weights=Py_True;
  
  static char *kwlist[] =
  {
    "f", "names", "weights", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|OO", kwlist,
				   &fname, &names, &weights))
     return NULL;
      
   self = (igraphmodule_GraphObject*)type->tp_alloc(type, 0);
   if (self != NULL) 
     {
       f=fopen(fname, "r");
       if (!f) {
	 PyErr_SetString(PyExc_IOError, strerror(errno));
	 return NULL;
       }
       if (igraph_read_graph_lgl(&self->g, f, PyObject_IsTrue(names), PyObject_IsTrue(weights)))
       {
	 igraphmodule_handle_igraph_error();
	 fclose(f);
	 return NULL;
       }
       fclose(f);
     }
   
   return (PyObject*)self;
}

/** \ingroup python_interface
 * \brief Writes the edge list to a file
 * \return none
 * \sa igraph_write_graph_edgelist
 */
static PyObject* igraphmodule_Graph_write_edgelist(igraphmodule_GraphObject *self,
						   PyObject *args,
						   PyObject *kwds)
{
  char* fname=NULL;
  FILE* f;
  
  static char *kwlist[] =
  {
    "f", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s", kwlist, &fname))
     return NULL;
      
  f=fopen(fname, "w");
  if (!f) {
    PyErr_SetString(PyExc_IOError, strerror(errno));
    return NULL;
  }
  if (igraph_write_graph_edgelist(&self->g, f))
  {
    igraphmodule_handle_igraph_error();
    fclose(f);
    return NULL;
  }
  fclose(f);
  
  Py_RETURN_NONE;
}

/** \ingroup python_interface
 * \brief Writes the edge list to a file in .ncol format
 * \return none
 * \sa igraph_write_graph_ncol
 */
static PyObject* igraphmodule_Graph_write_ncol(igraphmodule_GraphObject *self,
					       PyObject *args,
					       PyObject *kwds)
{
  char* fname=NULL;
  char* names="name";
  char* weights="weight";
  FILE* f;
  
  static char *kwlist[] =
  {
    "f", "names", "weights", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|zz", kwlist,
				   &fname, &names, &weights))
     return NULL;

  f=fopen(fname, "w");
  if (!f) {
    PyErr_SetString(PyExc_IOError, strerror(errno));
    return NULL;
  }
  if (igraph_write_graph_ncol(&self->g, f, names, weights))
  {
    igraphmodule_handle_igraph_error();
    fclose(f);
    return NULL;
  }
  fclose(f);
  
  Py_RETURN_NONE;
}

/** \ingroup python_interface
 * \brief Writes the edge list to a file in .lgl format
 * \return none
 * \sa igraph_write_graph_lgl
 */
static PyObject* igraphmodule_Graph_write_lgl(igraphmodule_GraphObject *self,
					      PyObject *args,
					      PyObject *kwds)
{
  char* fname=NULL;
  char* names="name";
  char* weights="weight";
  PyObject* isolates=Py_True;
  FILE* f;
  
  static char *kwlist[] =
  {
    "f", "names", "weights", "isolates", NULL
  };

  if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|zzO", kwlist,
				   &fname, &names, &weights, &isolates))
     return NULL;

  f=fopen(fname, "w");
  if (!f) {
    PyErr_SetString(PyExc_IOError, strerror(errno));
    return NULL;
  }
  if (igraph_write_graph_lgl(&self->g, f, names, weights,
			     PyObject_IsTrue(isolates)))
  {
    igraphmodule_handle_igraph_error();
    fclose(f);
    return NULL;
  }
  fclose(f);
  
  Py_RETURN_NONE;
}

/** \defgroup python_interface_internal Internal functions
 * \ingroup python_interface */

/** \ingroup python_interface_internal
 * \brief Returns the encapsulated igraph graph as a PyCObject
 * \return a new PyCObject
 */
static PyObject* igraphmodule_Graph___graph_as_cobject__(igraphmodule_GraphObject *self,
						 PyObject *args,
						 PyObject *kwds) 
{
  /*
  static char *kwlist[] = {"ref", NULL};
  PyObject *incref=Py_True;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "|O", kwlist, &incref)) return NULL;
  if (PyObject_IsTrue(incref)) Py_INCREF(self);
  */
  
  return PyCObject_FromVoidPtr((void*)&self->g, NULL);
}

/** \ingroup python_interface_internal
 * \brief Registers a destructor to be called when the object is destroyed
 * \return the previous destructor (if any)
 * Unimplemented.
 */
static PyObject* igraphmodule_Graph___register_destructor__(igraphmodule_GraphObject *self,
							    PyObject *args,
							    PyObject *kwds) 
{
  static char *kwlist[] = {"destructor", NULL};
  PyObject *destructor = NULL, *result;
  
  if (!PyArg_ParseTupleAndKeywords(args, kwds, "O", kwlist, &destructor)) return NULL;
  
  if (!PyCallable_Check(destructor)) {
    PyErr_SetString(PyExc_TypeError, "The destructor must be callable!");
    return NULL;
  }
  
  result=self->destructor;
  self->destructor=destructor;
  Py_INCREF(self->destructor);

  if (!result) Py_RETURN_NONE;
  
  return result;
}

/** \ingroup python_interface
 * Python type object referencing the methods Python calls when it performs various operations on an igraph (creating, printing and so on)
 */
static PyTypeObject igraphmodule_GraphType =
{
     PyObject_HEAD_INIT(NULL)                  // 
     0,                                        // ob_size
     "igraph.Graph",                           // tp_name
     sizeof(igraphmodule_GraphObject),         // tp_basicsize
     0,                                        // tp_itemsize
     (destructor)igraphmodule_Graph_dealloc,   // tp_dealloc
     0,                                        // tp_print
     0,                                        // tp_getattr
     0,                                        // tp_setattr
     0,                                        // tp_compare
     0,                                        // tp_repr
     0,                                        // tp_as_number
     0,                                        // tp_as_sequence
     0,                                        // tp_as_mapping
     0,                                        // tp_hash
     0,                                        // tp_call
     (reprfunc)igraphmodule_Graph_str,         // tp_str
     0,                                        // tp_getattro
     0,                                        // tp_setattro
     0,                                        // tp_as_buffer
     Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE, // tp_flags
     "igraph graph object"                     // tp_doc
}
;

/** \ingroup python_interface
 * \brief Method list of the \c igraph.Graph object type
 */
static PyMethodDef igraphmodule_Graph_methods[] = 
{
  ////////////////////////////
  // BASIC IGRAPH INTERFACE //
  ////////////////////////////
  
  // interface to igraph_vcount
  {"vcount", (PyCFunction)igraphmodule_Graph_vcount,
      METH_NOARGS,
      "Returns the number of vertices in the graph"
  },
  // interface to igraph_ecount
  {"ecount", (PyCFunction)igraphmodule_Graph_ecount,
      METH_NOARGS,
      "Returns the number of edges in the graph"
  },
  // interface to igraph_is_directed
  {"is_directed", (PyCFunction)igraphmodule_Graph_is_directed,
      METH_NOARGS,
      "Checks whether the graph is directed"
  },
  // interface to igraph_add_vertices
  {"add_vertices", (PyCFunction)igraphmodule_Graph_add_vertices,
      METH_VARARGS,
      "Adds vertices to the graph. The only parameter is the number of "
      "vertices to be added"
  },
  // interface to igraph_delete_vertices
  {"delete_vertices", (PyCFunction)igraphmodule_Graph_delete_vertices,
      METH_VARARGS,
      "Deletes vertices and all its edges from the graph. The only "
      "parameter is a list of the vertices to be added. It is allowed "
      "to provide a single integer instead of a list consisting of only "
      "one integer."
  },
  // interface to igraph_add_edges
  {"add_edges", (PyCFunction)igraphmodule_Graph_add_edges,
      METH_VARARGS,
      "Adds edges to the graph. The only parameter is a list of "
      "edges to be added. Every edge is represented with a tuple, "
      "containing the vertex IDs of the two endpoints. Vertices are "
      "enumerated from zero. It is allowed to provide a single pair "
      "instead of a list consisting of only one pair."
  },
  // interface to igraph_delete_edges
  {"delete_edges", (PyCFunction)igraphmodule_Graph_delete_edges,
      METH_VARARGS,
      "Removes edges from the graph. The only parameter is a list of "
      "edges to be removed. Every edge is represented with a tuple, "
      "containing the vertex IDs of the two endpoints. Vertices are "
      "enumerated from zero. It is allowed to provide a single pair "
      "instead of a list consisting of only one pair. Nonexistent "
      "edges will be silently ignored. All vertices will be kept, even "
      "if they lose all their edges."
  },
  // interface to igraph_degree
  {"degree", (PyCFunction)igraphmodule_Graph_degree,
      METH_VARARGS | METH_KEYWORDS,
      "Returns some vertex degrees from the graph.\n"
      "This method accepts a single vertex ID or a list of vertex IDs as a "
      "parameter, and returns the degree of the given vertices (in the form of "
      "a single integer or a list, depending on the input parameter). A "
      "second and a third argument may be passed as well, the second one "
      "meaning the type of degree to be returned (OUT for out-degrees, "
      "IN for in-degrees or ALL for the sum of them) and the third one "
      "meaning whether self-loops should be counted. The default for them is "
      "ALL and False. The type of degree is ignored for undirected graphs."
  },
  // interface to igraph_neighbors
  {"neighbors", (PyCFunction)igraphmodule_Graph_neighbors,
      METH_VARARGS | METH_KEYWORDS,
      "Returns adjacent vertices to a given vertex.\n"
      "This method accepts a single vertex ID as an argument, "
      "and returns the adjacent vertices of that vertex. An optional "
      "second argument allows the user to limit the result to only "
      "predecessors (IN), only successors (OUT) or both of them (ALL). "
      "The default behaviour is the latter one. "
      "This argument is ignored for undirected graphs."
  },
  
  //////////////////////
  // GRAPH GENERATORS //
  //////////////////////
  
  // interface to igraph_barabasi_game
  {"Barabasi", (PyCFunction)igraphmodule_Graph_Barabasi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph based on the Barabási-Albert model.\n"
      "The first two arguments are mandatory: the first one is the "
      "number of vertices, the second one is either the number of "
      "outgoing edges generated for each vertex or a list containing the "
      "number of outgoing edges for each vertex explicitly. The third "
      "argument is True if the out-degree of a given vertex should also "
      "increase its citation probability (as well as its in-degree), but "
      "it defaults to False. The fourth argument is True if the generated "
      "graph should be directed (default: False).\n\n"
      "Keywords for the arguments: n, m, outpref, directed"
  },
  
  // interface to igraph_erdos_renyi_game
  {"Erdos_Renyi", (PyCFunction)igraphmodule_Graph_Erdos_Renyi,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph based on the Erdõs-Rényi model.\n"
      "There are a total of five possible arguments, two of them are "
      "mutually exclusive. The first argument (keyword: n) is the number "
      "of vertices. The second and the third (keywords: p and m) define "
      "the density of the graph: if p is missing, there will be m edges; "
      "if m is missing, every edge will be present with a probability of "
      "p. These two arguments influence the same graph property (the "
      "number of edges) in two different ways, so either p or m must be "
      "present (but not both of them). The remaining two arguments are "
      "optional. The fourth argument (keyword: directed) is True if the "
      "generated graph should be directed (default: False), the fifth "
      "(keyword: loops) is True if self-loops are allowed (default: False)."
  },
  
  // interface to igraph_full_game
  {"Full", (PyCFunction)igraphmodule_Graph_Full,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a full graph (directed or undirected, with or without loops).\n"
      "The only mandatory argument (keyword: n) is the number "
      "of vertices. The remaining two arguments are optional. "
      "The second argument (keyword: directed) is True if the "
      "generated graph should be directed (default: False), the third "
      "(keyword: loops) is True if self-loops are allowed (default: False)."
  },
  
  // interface to igraph_growing_random_game
  {"Growing_Random", (PyCFunction)igraphmodule_Graph_Growing_Random,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a growing random graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "m -- The number of edges to add in each step (after adding a new vertex)\n"
      "directed -- whether the graph should be directed.\n"
      "            Optional, defaults to False.\n"
      "citation -- whether the new edges should originate from the most\n"
      "            recently added vertex.\n"
      "            Optional, defaults to False."
  },
  
  // interface to igraph_star
  {"Star", (PyCFunction)igraphmodule_Graph_Star,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a star graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "mode -- Gives the type of the star graph to create. Should be\n"
      "        one of the constants STAR_OUT, STAR_IN and STAR_UNDIRECTED.\n"
      "        Optional, defaults to STAR_UNDIRECTED.\n"
      "center -- Vertex ID for the central vertex in the star.\n"
      "          Optional, defaults to zero.\n"
  },
  
  // interface to igraph_lattice
  {"Lattice", (PyCFunction)igraphmodule_unimplemented,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a lattice. This function is yet unimplemented.\n\n"
      "Throws a NotImplementedError."
  },
  
  // interface to igraph_ring
  {"Ring", (PyCFunction)igraphmodule_Graph_Ring,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a ring graph.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the ring\n"
      "directed -- whether to create a directed ring.\n"
      "            Optional, defaults to False.\n"
      "mutual -- whether to create mutual edges in a directed ring.\n"
      "          Optional, defaults to False.\n"
      "          Ignored for undirected graphs.\n"
      "circular -- whether to create a closed ring.\n"
      "            Optional, defaults to True."
  },
  
  // interface to igraph_tree
  {"Tree", (PyCFunction)igraphmodule_Graph_Tree,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a tree in which almost all vertices have the same number of children.\n\n"
      "Keyword arguments:\n"
      "n -- The number of vertices in the graph\n"
      "children -- The number of children of a vertex in the graph\n"
      "type -- determines whether the tree should be directed, and if\n"
      "        this is the case, also its orientation. Must be one of\n"
      "        TREE_IN, TREE_OUT and TREE_UNDIRECTED.\n"
      "        Optional, defaults to TREE_UNDIRECTED\n"
  },
  
  // interface to igraph_adjacency
  {"Adjacency", (PyCFunction)igraphmodule_unimplemented,
      METH_VARARGS | METH_CLASS | METH_KEYWORDS,
      "Generates a graph from an adjacency matrix. This function is yet unimplemented.\n\n"
      "Throws a NotImplementedError."
  },
  
  /////////////////////////////////////
  // STRUCTURAL PROPERTIES OF GRAPHS //
  /////////////////////////////////////
  
  // interface to igraph_are_connected
  {"are_connected", (PyCFunction)igraphmodule_Graph_are_connected,
      METH_VARARGS | METH_KEYWORDS,
      "Decides whether two given vertices are directly connected.\n\n"
      "Keyword arguments:\n"
      "v1 -- the first vertex\n"
      "v2 -- the second vertex\n"
      "Returns true if there exists an edge from v1 to v2."
  },
  
  // interface to igraph_average_path_length
  {"average_path_length", (PyCFunction)igraphmodule_Graph_average_path_length,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the average path length in a graph.\n\n"
      "Keyword arguments:\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
      "unconn -- what to do when the graph is unconnected. If True, the\n"
      "          average of the geodesic lengths in the components is\n"
      "          calculated. Otherwise for all unconnected vertex pairs,\n"
      "          a path length equal to the number of vertices is used.\n"
  },
  
  // interface to igraph_betweenness
  {"betweenness", (PyCFunction)igraphmodule_Graph_betweenness,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates the betweenness of nodes in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the vertices for which the betweennesses must be returned.\n"
      "            Optional, defaults to all of the vertices in the graph.\n"
      "directed -- whether to consider directed paths.\n"
      "            Ignored for undirected graphs. Optional, defaults to True.\n"
  },
  
  // interface to igraph_bibcoupling
  {"bibcoupling", (PyCFunction)igraphmodule_Graph_bibcoupling,
      METH_VARARGS | METH_KEYWORDS,
      "Calculates bibliographic coupling values for given vertices in a graph.\n\n"
      "Keyword arguments:\n"
      "vertices -- the vertices to be analysed.\n"
      "Returns bibliographic coupling values for all given vertices in a matrix."
  },
  
   // interface to igraph_closeness
     {"closeness", (PyCFunction)igraphmodule_Graph_closeness,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates the closeness centralities of given nodes in a graph.\n\n"
	  "The closeness centerality of a vertex measures how easily other\n"
	  "vertices can be reached from it (or the other way: how easily it\n"
	  "can be reached from the other vertices). It is defined as the\n"
	  "number of the number of vertices minus one divided by the sum of\n"
	  "the lengths of all geodesics from/to the given vertex.\n\n"
	  "If the graph is not connected, and there is no path between two\n"
	  "vertices, the number of vertices is used instead the length of\n"
	  "the geodesic. This is always longer than the longest possible\n"
	  "geodesic.\n\n"
	  "Keyword arguments:\n"
	  "vertices -- the vertices for which the betweennesses must be returned.\n"
	  "            Optional, defaults to all of the vertices in the graph.\n"
	  "mode -- must be one of IN, OUT and ALL. IN means that the length of\n"
	  "        incoming paths, OUT means that the length of the outgoing\n"
	  "        paths must be calculated. ALL means that both of them must\n"
	  "        be calculated. Optional, defaults to ALL.\n"
     },   	    

   // interface to igraph_clusters
     {"clusters", (PyCFunction)igraphmodule_Graph_clusters,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates the (strong or weak) clusters for a given graph.\n\n"
	  "Keyword arguments:\n"
	  "mode -- must be either STRONG or WEAK, depending on the clusters\n"
	  "        being sought. Optional, defaults to STRONG.\n"
	  "Returns the component index for every node in the graph."
     },
     {"components", (PyCFunction)igraphmodule_Graph_clusters,
	  METH_VARARGS | METH_KEYWORDS,
	  "Alias for 'clusters'.\n\n"
	  "See the documentation of 'clusters' for details."
     },
   
   // interface to igraph_cocitation
     {"cocitation", (PyCFunction)igraphmodule_Graph_cocitation,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates cocitation scores for given vertices in a graph.\n\n"
	  "Keyword arguments:\n"
	  "vertices -- the vertices to be analysed.\n"
	  "Returns cocitation scores for all given vertices in a matrix."
     },
   
   // interface to igraph_diameter
     {"diameter", (PyCFunction)igraphmodule_Graph_diameter,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates the diameter of the graph.\n\n"
	  "Keyword arguments:\n"
	  "directed -- whether to consider directed paths.\n"
	  "            Ignored for undirected graphs. Optional, defaults to True.\n"
	  "unconn -- if True and the graph is undirected, the longest geodesic "
	  "          within a component will be returned. If False and the "
	  "          graph is undirected, the result is the number of vertices."
     },
   
   // interface to igraph_edge_betweenness
     {"edge_betweenness", (PyCFunction)igraphmodule_Graph_edge_betweenness,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates the edge betweennesses in a graph.\n\n"
	  "Keyword arguments:\n"
	  "directed -- whether to consider directed paths.\n"
	  "            Ignored for undirected graphs. Optional, defaults to True.\n"
	  "Returns a list with the edge betweennesses of all the edges.\n"
     },

   // interface to igraph_get_shortest_paths
     {"get_shortest_paths", (PyCFunction)igraphmodule_Graph_get_shortest_paths,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates the shortest paths from/to a given node in a graph.\n\n"
	  "Keyword arguments:\n"
	  "v    -- the source/destination for the calculated paths\n"
	  "mode -- the directionality of the paths. IN means to calculate\n"
	  "        incoming paths, OUT means to calculate outgoing paths,\n"
	  "        ALL means to calculate both ones. Ignored for undirected\n"
	  "        graphs. Optional, defaults to ALL\n"
	  "Returns at most one shortest path for every node in the graph in a\n"
	  "list. For unconnected graphs, some of the list elements will be\n"
	  "an empty list. Note that in case of mode=IN, the nodes in a path\n"
	  "are returned in reversed order!"
     },
   
   // interface to igraph_is_connected
     {"is_connected", (PyCFunction)igraphmodule_Graph_is_connected,
	  METH_VARARGS | METH_KEYWORDS,
	  "Decides whether a graph is connected.\n\n"
	  "Keyword arguments:\n"
	  "mode -- whether we should calculate strong or weak connectivity.\n"
	  "        Ignored for undirected graphs. Optional, defaults to\n"
	  "        STRONG."
     },
   
   // interface to igraph_shortest_paths
     {"shortest_paths", (PyCFunction)igraphmodule_Graph_shortest_paths,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates shortest path lengths for given nodes in a graph.\n\n"
	  "Keyword arguments:\n"
	  "vertices -- a list containing the vertex IDs which should be included in the result.\n"
	  "mode -- the type of shortest paths to be used for the calculation in directed graphs.\n"
	  "        OUT -- outgoing paths, IN -- incoming paths, ALL -- the directed graph\n"
	  "        is considered as an undirected one."
     },

   // interface to igraph_simplify
     {"simplify", (PyCFunction)igraphmodule_Graph_simplify,
	  METH_VARARGS | METH_KEYWORDS,
	  "Simplifies a graph by removing selp-loops and/or multiple edges."
	  "Keywords arguments:\n"
	  "multiple -- whether to remove multiple edges. Optional, defaults\n"
	  "            to True.\n"
	  "loops -- whether to remove loops. Optional, defaults to True.\n"
     },
   
   // interface to igraph_minimum_spanning_tree_unweighted and
   // igraph_minimum_spanning_tree_prim
     {"spanning_tree", (PyCFunction)igraphmodule_Graph_spanning_tree,
	  METH_VARARGS | METH_KEYWORDS,
	  "Calculates a minimum spanning tree for a graph (weighted or unweighted)\n\n"
	  "Keyword arguments:\n"
	  "weights -- a vector containing weights for every edge in the graph.\n"
	  "           If omitted, every edge is assumed to have an equal weight.\n"
	  "Returns the spanning tree as an igraph.Graph object."
     },

   // interface to igraph_subcomponent
     {"subcomponent", (PyCFunction)igraphmodule_Graph_subcomponent,
	  METH_VARARGS | METH_KEYWORDS,
	  "Returns the indices of vertices which are in the same component as a given vertex.\n\n"
	  "Keyword arguments:\n"
	  "v -- the index of the vertex used as the source/destination\n"
	  "mode -- if equals to IN, returns the vertex IDs from where the\n"
	  "        given vertex can be reached. If equals to OUT, returns the\n"
	  "        vertex IDs which are reachable from the given vertex. If\n"
	  "        equals to ALL, returns all vertices within the same component\n"
	  "        as the given vertex, ignoring edge directions. Note that this\n"
	  "        not equals to calculating the union of the results of IN and OUT.\n"
     },

   // interface to igraph_subgraph
     {"subgraph", (PyCFunction)igraphmodule_Graph_subgraph,
	  METH_VARARGS | METH_KEYWORDS,
	  "Returns a subgraph based on the given vertices.\n\n"
	  "Keyword arguments:\n"
	  "vertices -- a list containing the vertex IDs which should be included in the result.\n"
     },

   //////////////////////
   // LAYOUT FUNCTIONS //
   //////////////////////
   
  // interface to igraph_layout_circle
  {"layout_circle", (PyCFunction)igraphmodule_Graph_layout_circle,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph uniformly on a circle.\n\n"
      "Returns the calculated coordinate pairs in a vector."
  },

  // interface to igraph_layout_kamada_kawai
  {"layout_kamada_kawai", (PyCFunction)igraphmodule_Graph_layout_kamada_kawai,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a plane according the Kamada-Kawai algorithm.\n\n"
      "This is a force directed layout, see Kamada, T. and Kawai, S.:\n"
      "An Algorithm for Drawing General Undirected Graphs.\n"
      "Information Processing Letters, 31/1, 7--15, 1989.\n\n"
      "Keyword arguments:\n"
      "n       -- the number of iterations to perform. Optional, defaults to 1000.\n"
      "sigma   -- the standard base deviation of the position change proposals.\n"
      "           Optional, defaults to the number of vertices * 0.25\n"
      "initemp -- initial temperature of the simulated annealing.\n"
      "           Optional, defaults to 10.\n"
      "coolexp -- cooling exponent of the simulated annealing.\n"
      "           Optional, defaults to 0.99\n"
      "kkconst -- the Kamada-Kawai vertex attraction constant. Optional,\n"
      "           defaults to the square of the number of vertices.\n"
  },
  
  // interface to igraph_layout_fruchterman_reingold
  {"layout_fruchterman_reingold", (PyCFunction)igraphmodule_Graph_layout_fruchterman_reingold,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices on a plane according the Fruchterman-Reingold algorithm.\n\n"
      "This is a force directed layout, see Fruchterman, T. M. J. and Reingold, E. M.:\n"
      "Graph Drawing by Force-directed Placement.\n"
      "Software -- Practice and Experience, 21/11, 1129--1164, 1991\n\n"
      "Keyword arguments:\n"
      "n          -- the number of iterations to perform. Optional, defaults to 500.\n"
      "maxdelta   -- the maximum distance to move a vertex in an iteration\n"
      "              Optional, defaults to the number of vertices\n"
      "area       -- the area parameter of the algorithm. Optional, defaults\n"
      "              to the square of the number of vertices\n"
      "coolexp    -- the cooling exponent of the simulated annealing.\n"
      "              Optional, defaults to 0.99\n"
      "repulserad -- Determines the radius at which vertex-vertex repulsion\n"
      "              cancels out attraction of adjacent vertices.\n"
  },
  
  // interface to igraph_layout_random
  {"layout_random", (PyCFunction)igraphmodule_Graph_layout_random,
      METH_VARARGS | METH_KEYWORDS,
      "Places the vertices of the graph randomly.\n\n"
      "Returns the \"calculated\" coordinate pairs in a vector."
  },
   
   //////////////////////////////////////////////////////
   // CONVERT A GRAPH TO EDGE LIST OR ADJACENCY MATRIX //
   //////////////////////////////////////////////////////
   
   // interface to igraph_get_edgelist
     {"get_adjacency", (PyCFunction)igraphmodule_Graph_get_adjacency,
	  METH_VARARGS | METH_KEYWORDS,
	  "Returns the edge list of a graph.\n\n"
	  "Keyword arguments:\n"
	  "type -- either GET_ADJACENCY_LOWER (uses the lower triangle\n"
	  "        of the matrix) or GET_ADJACENCY_UPPER (uses the upper\n"
	  "        triangle) or GET_ADJACENCY_BOTH (uses both parts).\n"
	  "        Optional, defaults to GET_ADJACENCY_BOTH, ignored for\n"
	  "        directed graphs."
     },
   
   // interface to igraph_get_edgelist
     {"get_edgelist", (PyCFunction)igraphmodule_Graph_get_edgelist,
	  METH_NOARGS,
	  "Returns the edge list of a graph."
     },

  ///////////////////////////////
  // LOADING AND SAVING GRAPHS //
  ///////////////////////////////
  
  // interface to igraph_read_graph_edgelist
  {"Read_Edgelist", (PyCFunction)igraphmodule_Graph_Read_Edgelist,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an edge list from a file and creates a graph based on it.\n"
      "Please note that the vertex indices are zero-based.\n\n"
      "Keyword arguments:\n"
      "f        -- the name of the file\n"
      "directed -- whether the generated graph should be directed.\n"
      "            Optional, defaults to True.\n\n"
  },
  // interface to igraph_read_graph_ncol
  {"Read_Ncol", (PyCFunction)igraphmodule_Graph_Read_Ncol,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an .ncol file used by LGL, also useful for creating graphs\n"
      "from \"named\" (and optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "documentation of LGL regarding the exact format description:\n"
      "http://bioinformatics.icmb.utexas.edu/bgl\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file\n"
      "names   -- logical value. If True, the vertex names are added as a\n"
      "           vertex attribute called 'name'. Optional, defaults to\n"
      "           True.\n"
      "weights -- logical value. If True, the edge weights are added as an\n"
      "           edge attribute called 'weight'. Optional, defaults to\n"
      "           True.\n"
  },
  // interface to igraph_read_graph_lgl
  {"Read_Lgl", (PyCFunction)igraphmodule_Graph_Read_Lgl,
      METH_VARARGS | METH_KEYWORDS | METH_CLASS,
      "Reads an .lgl file used by LGL, also useful for creating graphs\n"
      "from \"named\" (and optionally weighted) edge lists.\n\n"
      "This format is used by the Large Graph Layout program. See the\n"
      "documentation of LGL regarding the exact format description:\n"
      "http://bioinformatics.icmb.utexas.edu/bgl\n\n"
      "LGL originally cannot deal with graphs containing multiple or loop\n"
      "edges, but this condition is not checked here, as igraph is happy\n"
      "with these.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file\n"
      "names   -- logical value. If True, the vertex names are added as a\n"
      "           vertex attribute called 'name'. Optional, defaults to\n"
      "           True.\n"
      "weights -- logical value. If True, the edge weights are added as an\n"
      "           edge attribute called 'weight'. Optional, defaults to\n"
      "           True.\n"
  },
  // interface to igraph_write_graph_edgelist
  {"write_edgelist", (PyCFunction)igraphmodule_Graph_write_edgelist,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file. Directed edges are\n"
      "written in (from, to) order.\n\n"
      "Keyword arguments:\n"
      "f -- the name of the file to be written\n"
  },
  // interface to igraph_write_graph_ncol
  {"write_ncol", (PyCFunction)igraphmodule_Graph_write_ncol,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file in .ncol format.\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call simplify() before saving.\n\n"
      "Keyword arguments:\n"
      "f       -- the name of the file to be written\n"
      "names   -- the name of the vertex attribute containing the name\n"
      "           of the vertices. Optional, defaults to 'name'. If you\n"
      "           don't want to store vertex names, supply None here.\n"
      "weights -- the name of the edge attribute containing the weight\n"
      "           of the vertices. Optional, defaults to 'weight'. If you\n"
      "           don't want to store weights, supply None here.\n\n"
      "PLEASE NOTE THAT IT IS VITAL TO SUPPLY CORRECT 'names' AND 'weights'\n"
      "PARAMETERS, otherwise the underlying igraph library will segfault\n"
      "when it tries to reach nonexistent attributes. This issue will be\n"
      "corrected soon.\n"
  },
  // interface to igraph_write_graph_lgl
  {"write_lgl", (PyCFunction)igraphmodule_Graph_write_lgl,
      METH_VARARGS | METH_KEYWORDS,
      "Writes the edge list of a graph to a file in .lgl format.\n"
      "Note that multiple edges and/or loops break the LGL software,\n"
      "but igraph does not check for this condition. Unless you know\n"
      "that the graph does not have multiple edges and/or loops, it\n"
      "is wise to call simplify() before saving.\n\n"
      "Keyword arguments:\n"
      "f        -- the name of the file to be written\n"
      "names    -- the name of the vertex attribute containing the name\n"
      "            of the vertices. Optional, defaults to 'name'. If you\n"
      "            don't want to store vertex names, supply None here.\n"
      "weights  -- the name of the edge attribute containing the weight\n"
      "            of the vertices. Optional, defaults to 'weight'. If you\n"
      "            don't want to store weights, supply None here.\n"
      "isolates -- whether to include isolated vertices in the output.\n"
      "            Optional, defaults to True.\n\n"
      "PLEASE NOTE THAT IT IS VITAL TO SUPPLY CORRECT 'names' AND 'weights'\n"
      "PARAMETERS, otherwise the underlying igraph library will segfault\n"
      "when it tries to reach nonexistent attributes. This issue will be\n"
      "corrected soon.\n"
  },
  
  ////////////////////////////////////
  // INTERNAL/DEVELOPMENT FUNCTIONS //
  ////////////////////////////////////
  {"__graph_as_cobject__", (PyCFunction)igraphmodule_Graph___graph_as_cobject__,
      METH_VARARGS | METH_KEYWORDS,
      "Returns the igraph graph encapsulated by the Python object as\n"
      "a PyCObject (which is barely a regular C pointer). This function\n"
      "should not be used directly by igraph users, it is useful only\n"
      "in the case when the underlying igraph object must be passed to\n"
      "another C code through Python.\n\n"
      /*"Keyword arguments:\n"
      "ref -- increases the reference count of the graph when True.\n"
      "       Optional, defaults to False.\n"*/
  },
  {"__register_destructor__", (PyCFunction)igraphmodule_Graph___register_destructor__,
      METH_VARARGS | METH_KEYWORDS,
      "Registers a destructor to be called when the object is freed by "
      "Python. This function should not be used directly by igraph users."
  },
  {NULL}

}
;

/** \ingroup python_interface
 * \brief Method table for the igraph Python module
 */
static PyMethodDef igraphmodule_methods[] = 
{
   {NULL, NULL, 0, NULL}
};

#ifndef PyMODINIT_FUNC
#define PyMODINIT_FUNC void
#endif

PyMODINIT_FUNC
initigraph(void)
{
   PyObject* m;  ///< igraph module object
   
   igraphmodule_GraphType.tp_new = igraphmodule_Graph_new;
   igraphmodule_GraphType.tp_init = (initproc)igraphmodule_Graph_init;
   igraphmodule_GraphType.tp_methods = igraphmodule_Graph_methods;
   
   if (PyType_Ready(&igraphmodule_GraphType) < 0) return;

   igraphmodule_InternalError =
     PyErr_NewException("igraph.InternalError", NULL, NULL);
   
   Py_INCREF(igraphmodule_InternalError);
   
   m = Py_InitModule3("igraph", igraphmodule_methods,
		      "Python interface for the igraph library");
   Py_INCREF(&igraphmodule_GraphType);
   
   PyModule_AddObject(m, "Graph", (PyObject*)&igraphmodule_GraphType);
   PyModule_AddObject(m, "InternalError", igraphmodule_InternalError);
   PyModule_AddIntConstant(m, "OUT", IGRAPH_OUT);
   PyModule_AddIntConstant(m, "IN", IGRAPH_IN);
   PyModule_AddIntConstant(m, "ALL", IGRAPH_ALL);
   PyModule_AddIntConstant(m, "STAR_OUT", IGRAPH_STAR_OUT);
   PyModule_AddIntConstant(m, "STAR_IN", IGRAPH_STAR_IN);
   PyModule_AddIntConstant(m, "STAR_UNDIRECTED", IGRAPH_STAR_UNDIRECTED);
   PyModule_AddIntConstant(m, "TREE_OUT", IGRAPH_TREE_OUT);
   PyModule_AddIntConstant(m, "TREE_IN", IGRAPH_TREE_IN);
   PyModule_AddIntConstant(m, "TREE_UNDIRECTED", IGRAPH_TREE_UNDIRECTED);
   PyModule_AddIntConstant(m, "STRONG", IGRAPH_STRONG);
   PyModule_AddIntConstant(m, "WEAK", IGRAPH_WEAK);
   PyModule_AddIntConstant(m, "GET_ADJACENCY_UPPER", IGRAPH_GET_ADJACENCY_UPPER);
   PyModule_AddIntConstant(m, "GET_ADJACENCY_LOWER", IGRAPH_GET_ADJACENCY_LOWER);
   PyModule_AddIntConstant(m, "GET_ADJACENCY_BOTH", IGRAPH_GET_ADJACENCY_BOTH);
}

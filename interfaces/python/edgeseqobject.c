#include "edgeseqobject.h"
// #include "edgeobject.h"
#include "common.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_edgeseq Edge sequence object
 */

/**
 * \ingroup python_interface_edgeseq
 * \brief Allocate a new edge sequence object for a given graph
 * \param g the graph object being referenced
 * \return the allocated PyObject
 */
PyObject* igraphmodule_EdgeSeq_New(igraphmodule_GraphObject *g) {
  igraphmodule_EdgeSeqObject* o;
  
  o=PyObject_GC_New(igraphmodule_EdgeSeqObject, &igraphmodule_EdgeSeqType);
  o->gref=PyWeakref_NewRef((PyObject*)g, NULL);
  //Py_INCREF(g);
  PyObject_GC_Track(o);
  
  RC_ALLOC("EdgeSeq", o);
  
  return (PyObject*)o;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.EdgeSeq object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_EdgeSeq_traverse(igraphmodule_EdgeSeqObject *self,
				  visitproc visit, void *arg) {
  int vret;

  RC_TRAVERSE("EdgeSeq", self);
  
  if (self->gref) {
    vret=visit(self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_EdgeSeq_clear(igraphmodule_EdgeSeqObject *self) {
  PyObject *tmp;

  PyObject_GC_UnTrack(self);
  
  tmp=self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Deallocates a Python representation of a given edge sequence object
 */
void igraphmodule_EdgeSeq_dealloc(igraphmodule_EdgeSeqObject* self) {
  igraphmodule_EdgeSeq_clear(self);

  RC_DEALLOC("EdgeSeq", self);
  
  PyObject_GC_Del(self);
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Returns the length of the sequence (i.e. the number of edges in the graph)
 */
int igraphmodule_EdgeSeq_sq_length(igraphmodule_EdgeSeqObject* self) {
  igraph_t *g;
  
  g=&((igraphmodule_GraphObject*)PyWeakref_GetObject(self->gref))->g;
  if ((PyObject*)g == Py_None) return 0;
  
  return (int)igraph_ecount(g);
}

/**
 * \ingroup python_interface_edgeseq
 * \brief Returns the item at the given index in the sequence
 */
PyObject* igraphmodule_EdgeSeq_sq_item(igraphmodule_EdgeSeqObject* self,
				       int i) {
  igraphmodule_GraphObject *o;
  igraph_t *g;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  g=&o->g;
  if (i<0 || i>=(int)igraph_ecount(g)) {
    PyErr_SetString(PyExc_IndexError, "edge index out of range");
    return NULL;
  }
  /// @todo caching
  return igraphmodule_Edge_New(self->gref, i);
}

/** \ingroup python_interface_edge
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_EdgeSeq_attributes(igraphmodule_EdgeSeqObject* self) {
  vector_t t;
  vector_ptr_t ns;
  long result;
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_edge_attributes(o, NULL, NULL);
}

/**
 * \ingroup python_interface_edgeseq
 * Method table for the \c igraph.EdgeSeq object
 */
PyMethodDef igraphmodule_EdgeSeq_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_EdgeSeq_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph's edges\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_edgeseq
 * This is the collection of functions necessary to implement the
 * edge sequence as a real sequence (e.g. allowing to reference
 * edges by indices)
 */
static PySequenceMethods igraphmodule_EdgeSeq_as_sequence = {
  (inquiry)igraphmodule_EdgeSeq_sq_length,
  0,               /* sq_concat */
  0,               /* sq_repeat */
  (intargfunc)igraphmodule_EdgeSeq_sq_item, /* sq_item */
  0,                                          /* sq_slice */
  0,                                          /* sq_ass_item */
  0,                                          /* sq_ass_slice */
  0,                                          /* sq_contains */
  0,                                          /* sq_inplace_concat */
  0,                                          /* sq_inplace_repeat */
};

/** \ingroup python_interface_edgeseq
 * Python type object referencing the methods Python calls when it performs various operations on
 * an edge sequence of a graph
 */
PyTypeObject igraphmodule_EdgeSeqType =
{
  PyObject_HEAD_INIT(NULL)                  //
  0,                                        // ob_size
  "igraph.EdgeSeq",                         // tp_name
  sizeof(igraphmodule_EdgeSeqObject),       // tp_basicsize
  0,                                        // tp_itemsize
  (destructor)igraphmodule_EdgeSeq_dealloc, // tp_dealloc
  0,                                        // tp_print
  0,                                        // tp_getattr
  0,                                        // tp_setattr
  0,                                        // tp_compare
  0,                                        // tp_repr
  0,                                        // tp_as_number
  &igraphmodule_EdgeSeq_as_sequence,        // tp_as_sequence
  0,                                        // tp_as_mapping
  0,                                        // tp_hash
  0,                                        // tp_call
  0,                                        // tp_str
  0,                                        // tp_getattro
  0,                                        // tp_setattro
  0,                                        // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "igraph edge sequence object",            // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  0,                                          // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  igraphmodule_EdgeSeq_methods,             // tp_methods
};


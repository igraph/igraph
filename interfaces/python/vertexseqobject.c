#include "vertexseqobject.h"
#include "vertexobject.h"
#include "common.h"

/**
 * \ingroup python_interface
 * \defgroup python_interface_vertexseq Vertex sequence object
 */

/**
 * \ingroup python_interface_vertexseq
 * \brief Allocate a new vertex sequence object for a given graph
 * \param g the graph object being referenced
 * \return the allocated PyObject
 */
PyObject* igraphmodule_VertexSeq_New(igraphmodule_GraphObject *g) {
  igraphmodule_VertexSeqObject* o;
  
  o=PyObject_GC_New(igraphmodule_VertexSeqObject, &igraphmodule_VertexSeqType);
  o->gref=PyWeakref_NewRef((PyObject*)g, NULL);
  //Py_INCREF(g);
  PyObject_GC_Track(o);
  
  RC_ALLOC("VertexSeq", o);
  
  return (PyObject*)o;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Support for cyclic garbage collection in Python
 * 
 * This is necessary because the \c igraph.VertexSeq object contains several
 * other \c PyObject pointers and they might point back to itself.
 */
int igraphmodule_VertexSeq_traverse(igraphmodule_VertexSeqObject *self,
				    visitproc visit, void *arg) {
  int vret;

  RC_TRAVERSE("VertexSeq", self);
  
  if (self->gref) {
    vret=visit(self->gref, arg);
    if (vret != 0) return vret;
  }
  
  return 0;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Clears the graph object's subobject (before deallocation)
 */
int igraphmodule_VertexSeq_clear(igraphmodule_VertexSeqObject *self) {
  PyObject *tmp;

  PyObject_GC_UnTrack(self);
  
  tmp=self->gref;
  self->gref=NULL;
  Py_XDECREF(tmp);

  return 0;
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Deallocates a Python representation of a given vertex sequence object
 */
void igraphmodule_VertexSeq_dealloc(igraphmodule_VertexSeqObject* self) {
  igraphmodule_VertexSeq_clear(self);

  RC_DEALLOC("VertexSeq", self);
  
  PyObject_GC_Del(self);
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Returns the length of the sequence (i.e. the number of vertices in the graph)
 */
int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject* self) {
  igraph_t *g;
  
  g=&((igraphmodule_GraphObject*)PyWeakref_GetObject(self->gref))->g;
  if ((PyObject*)g == Py_None) return 0;
  
  return (int)igraph_vcount(g);
}

/**
 * \ingroup python_interface_vertexseq
 * \brief Returns the item at the given index in the sequence
 */
PyObject* igraphmodule_VertexSeq_sq_item(igraphmodule_VertexSeqObject* self,
					 int i) {
  igraphmodule_GraphObject *o;
  igraph_t *g;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;
  
  g=&o->g;
  if (i<0 || i>=(int)igraph_vcount(g)) {
    PyErr_SetString(PyExc_IndexError, "vertex index out of range");
    return NULL;
  }
  /// @todo caching
  return igraphmodule_Vertex_New(self->gref, i);
}

/** \ingroup python_interface_vertex
 * \brief Returns the list of attribute names
 */
PyObject* igraphmodule_VertexSeq_attributes(igraphmodule_VertexObject* self) {
  vector_t t;
  vector_ptr_t ns;
  long result;
  igraphmodule_GraphObject *o;
  
  o=(igraphmodule_GraphObject*)igraphmodule_resolve_graph_weakref(self->gref);
  if (!o) return NULL;

  return igraphmodule_Graph_vertex_attributes(o, NULL, NULL);
}

/**
 * \ingroup python_interface_vertexseq
 * Method table for the \c igraph.VertexSeq object
 */
PyMethodDef igraphmodule_VertexSeq_methods[] = {
  {"attributes", (PyCFunction)igraphmodule_VertexSeq_attributes,
      METH_NOARGS,
      "Returns the attribute list of the graph's vertices\n"
  },
  {NULL}
};

/**
 * \ingroup python_interface_vertexseq
 * This is the collection of functions necessary to implement the
 * vertex sequence as a real sequence (e.g. allowing to reference
 * vertices by indices)
 */
static PySequenceMethods igraphmodule_VertexSeq_as_sequence = {
  (inquiry)igraphmodule_VertexSeq_sq_length,
  0,               /* sq_concat */
  0,               /* sq_repeat */
  (intargfunc)igraphmodule_VertexSeq_sq_item, /* sq_item */
  0,                                          /* sq_slice */
  0,                                          /* sq_ass_item */
  0,                                          /* sq_ass_slice */
  0,                                          /* sq_contains */
  0,                                          /* sq_inplace_concat */
  0,                                          /* sq_inplace_repeat */
};

/** \ingroup python_interface_vertexseq
 * Python type object referencing the methods Python calls when it performs various operations on
 * a vertex sequence of a graph
 */
PyTypeObject igraphmodule_VertexSeqType =
{
  PyObject_HEAD_INIT(NULL)                    //
  0,                                          // ob_size
  "igraph.VertexSeq",                         // tp_name
  sizeof(igraphmodule_VertexSeqObject),       // tp_basicsize
  0,                                          // tp_itemsize
  (destructor)igraphmodule_VertexSeq_dealloc,  // tp_dealloc
  0,                                          // tp_print
  0,                                          // tp_getattr
  0,                                          // tp_setattr
  0,                                          // tp_compare
  0,                                          // tp_repr
  0,                                          // tp_as_number
  &igraphmodule_VertexSeq_as_sequence,        // tp_as_sequence
  0,                                          // tp_as_mapping
  0,                                          // tp_hash
  0,                                          // tp_call
  0,                                          // tp_str
  0,                                          // tp_getattro
  0,                                          // tp_setattro
  0,                                          // tp_as_buffer
  Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_GC, // tp_flags
  "igraph vertex sequence object",            // tp_doc
  0,                                          // tp_traverse
  0,                                          // tp_clear
  0,                                          // tp_richcompare
  0,                                          // tp_weaklistoffset
  0,                                          // tp_iter
  0,                                          // tp_iternext
  igraphmodule_VertexSeq_methods,             // tp_methods
};


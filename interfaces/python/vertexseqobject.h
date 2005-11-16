#ifndef PYTHON_VERTEXSEQOBJECT_H
#define PYTHON_VERTEXSEQOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_vertexseq
 * \brief A structure representing the vertex sequence of a graph
 */
typedef struct
{
  PyObject_HEAD
  PyObject* gref;
} igraphmodule_VertexSeqObject;

PyObject* igraphmodule_VertexSeq_New(igraphmodule_GraphObject *g);
int igraphmodule_VertexSeq_traverse(igraphmodule_VertexSeqObject *self,
					   visitproc visit, void *arg);
int igraphmodule_VertexSeq_clear(igraphmodule_VertexSeqObject *self);
void igraphmodule_VertexSeq_dealloc(igraphmodule_VertexSeqObject* self);

int igraphmodule_VertexSeq_sq_length(igraphmodule_VertexSeqObject *self);

PyTypeObject igraphmodule_VertexSeqType;

#endif

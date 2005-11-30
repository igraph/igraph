#ifndef PYTHON_EDGESEQOBJECT_H
#define PYTHON_EDGESEQOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_edgeseq
 * \brief A structure representing the edge sequence of a graph
 */
typedef struct
{
  PyObject_HEAD
  PyObject* gref;
} igraphmodule_EdgeSeqObject;

PyObject* igraphmodule_EdgeSeq_New(igraphmodule_GraphObject *g);
int igraphmodule_EdgeSeq_traverse(igraphmodule_EdgeSeqObject *self,
				  visitproc visit, void *arg);
int igraphmodule_EdgeSeq_clear(igraphmodule_EdgeSeqObject *self);
void igraphmodule_EdgeSeq_dealloc(igraphmodule_EdgeSeqObject* self);

int igraphmodule_EdgeSeq_sq_length(igraphmodule_EdgeSeqObject *self);

PyTypeObject igraphmodule_EdgeSeqType;

#endif

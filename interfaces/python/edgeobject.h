#ifndef PYTHON_EDGEOBJECT_H
#define PYTHON_EDGEOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_edge
 * \brief A structure representing an edge of a graph
 */
typedef struct
{
  PyObject_HEAD
  PyObject* gref;
  long idx;
} igraphmodule_EdgeObject;

int igraphmodule_Edge_traverse(igraphmodule_EdgeObject *self,
			       visitproc visit, void *arg);
int igraphmodule_Edge_clear(igraphmodule_EdgeObject *self);
void igraphmodule_Edge_dealloc(igraphmodule_EdgeObject* self);

PyObject* igraphmodule_Edge_New(PyObject *gref, long idx);
PyObject* igraphmodule_Edge_str(igraphmodule_EdgeObject *self);
//int igraphmodule_Edge_attribute_count(igraphmodule_EdgeObject* self);
PyObject* igraphmodule_Edge_attributes(igraphmodule_EdgeObject* self);
/*PyObject* igraphmodule_Edge_get_attribute(igraphmodule_EdgeObject* self,
					   PyObject* s);
int igraphmodule_Edge_set_attribute(igraphmodule_EdgeObject* self, PyObject* k, PyObject* v);*/

PyTypeObject igraphmodule_EdgeType;

#endif

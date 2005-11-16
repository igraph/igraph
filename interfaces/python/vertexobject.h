#ifndef PYTHON_VERTEXOBJECT_H
#define PYTHON_VERTEXOBJECT_H

#include <Python.h>
#include "graphobject.h"

/**
 * \ingroup python_interface_vertex
 * \brief A structure representing a vertex of a graph
 */
typedef struct
{
  PyObject_HEAD
  PyObject* gref;
  long idx;
} igraphmodule_VertexObject;

int igraphmodule_Vertex_traverse(igraphmodule_VertexObject *self,
					   visitproc visit, void *arg);
int igraphmodule_Vertex_clear(igraphmodule_VertexObject *self);
void igraphmodule_Vertex_dealloc(igraphmodule_VertexObject* self);

PyObject* igraphmodule_Vertex_New(PyObject *gref, long idx);
PyObject* igraphmodule_Vertex_str(igraphmodule_VertexObject *self);
//int igraphmodule_Vertex_attribute_count(igraphmodule_VertexObject* self);
PyObject* igraphmodule_Vertex_attributes(igraphmodule_VertexObject* self);
/*PyObject* igraphmodule_Vertex_get_attribute(igraphmodule_VertexObject* self,
					   PyObject* s);
int igraphmodule_Vertex_set_attribute(igraphmodule_VertexObject* self, PyObject* k, PyObject* v);*/

PyTypeObject igraphmodule_VertexType;

#endif

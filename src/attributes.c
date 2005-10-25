/* -*- mode: C -*-  */
/* 
   IGraph library.
   Copyright (C) 2005  Gabor Csardi <csardi@rmki.kfki.hu>
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
   Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.

*/

#include "igraph.h"
#include "memory.h"

/**
 * \ingroup attributes
 * \brief Adds a graph attribute.
 * 
 * Attributes have to be added by calling this function before setting
 * or getting them.
 * @param graph A graph object.
 * @param name The name of the attribute to install.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of graph
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_add_graph_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_add(&graph->gal, name);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Removes a graph attribute.
 * 
 * @param graph A graph object.
 * @param name The name of the attribute to remove.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of graph
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_remove_graph_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_remove(&graph->gal, name);  
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the value of a graph attribute.
 * 
 * @param graph A graph object.
 * @param name The name of the attribute to query.
 * @param value This will be set to the value of the attribute.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of graph
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_get_graph_attribute(igraph_t *graph, const char *name, 
			       real_t *value) {
  igraph_attribute_list_get(&graph->gal, name, 0, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Sets the value of a graph attribute.
 * 
 * @param graph A graph object.
 * @param name The name of the attribute to set.
 * @param value The new value of the attribute.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of graph
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_set_graph_attribute(igraph_t *graph, const char *name, 
			       real_t value) {
  igraph_attribute_list_set(&graph->gal, name, 0, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the list of installed graph attributes.
 * 
 * @param graph A graph object.
 * @param l This string array will contain the names of the 
 *        attributes. It should be initialized and will be resized.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of graph
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_list_graph_attributes(igraph_t *graph, igraph_strarray_t *l) {
  igraph_attribute_list_list(&graph->gal, l);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Adds a vertex attribute.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to install.
 * @return Error code.
 * 
 * Time complexity: <code>O(|V|)</code>, the number of vertices in the
 * graph.
 */

int igraph_add_vertex_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_add(&graph->val, name);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Removes a vertex attribute.
 *
 * @param graph A graph object.
 * @param name The name of the attribute to remove.
 * @return Error code.
 *
 * Time complexity: <code>O(|V|)</code>, assuming that the graph has
 * <code>O(1)</code> vertex attributes. <code>|V|</code> is the number
 * of vertices.
 */

int igraph_remove_vertex_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_remove(&graph->val, name);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the value of a vertex attribute for a single vertex
 * 
 * @param graph The graph object.
 * @param name The name of the vertex attribute.
 * @param v The id of the vertex of which the attribute is requested.
 * @param value Pointer to a real number, the result will be stored
 *        here.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>, assuming that the graph has
 * <code>O(1)</code> vertex attributes installed.
 */

int igraph_get_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, real_t *value) {
  igraph_attribute_list_get(&graph->val, name, v, value);
  return 0;  
}

/**
 * \ingroup attributes
 * \brief Set the value of a vertex attribute for a single vertex.
 * 
 * @param graph The graph object.
 * @param name Name of the vertex attribute.
 * @param v The id of the vertex of which the attribute is set.
 * @param value The new value of the attribute.
 * @return Error code.
 *
 * Time complexity: <code>O(1)</code>, assuming that the graph has
 * <code>O(1)</code> vertex attributes installed.
 */

int igraph_set_vertex_attribute(igraph_t *graph, const char *name, 
				long int v, real_t value) {
  igraph_attribute_list_set(&graph->val, name, v, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Query the value of a vertex attribute for many vertices.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to get.
 * @param v Vector with the vertex ids of the vertices of which the
 *        attribute will be returned.
 * @param value Pointer to an initialized vector, the result will be
 *        stored here. It will be resized if needed.
 * @return Error code.
 * 
 * Time complexity: <code>O(|v|)</code>, the number of queried
 * vertices, assuming the graph has <code>O(1)</code> vertex
 * attributes. 
 */

int igraph_get_vertex_attributes(igraph_t *graph, const char *name, 
				 vector_t *v, vector_t *value) {
  igraph_attribute_list_gets(&graph->val, name, v, value);
  return 0;  
}

/**
 * \ingroup attributes
 * \brief Set the value of a vertex attribute for many vertices.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to set.
 * @param v Vector with the vertex ids of the vertices of which the
 *        attribute will be set.
 * @param value The new value(s) of the attribute. This vector may be of
 *        different length than <code>v</code>, if it is shorter it
 *        will be recycled (ie. after the last element the first one
 *        is used again), if it is longer the unneeded values are
 *        ignored. Thus it is easy to set an attribute to a single
 *        constant value for many vertices, just give a vector of
 *        length 1 here.
 * @return Error code.
 * 
 * Time complexity: <code>O(|v|)</code>, the number of affected
 * vertices, assuming the graph has <code>O(1)</code> vertex
 * attributes. 
 */

int igraph_set_vertex_attributes(igraph_t *graph, const char *name, 
				 vector_t *v, vector_t *value) {
  igraph_attribute_list_sets(&graph->val, name, v, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the list of installed vertex attributes.
 * 
 * @param graph A graph object.
 * @param l This string array will contain the names of the 
 *        attributes. It should be initialized and will be resized.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of vertex
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_list_vertex_attributes(igraph_t *graph, igraph_strarray_t *l) {
  igraph_attribute_list_list(&graph->val, l);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Adds an edge attribute.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to install.
 * @return Error code.
 * 
 * Time complexity: <code>O(|E|)</code>, the number of edges in the
 * graph.
 */

int igraph_add_edge_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_add(&graph->eal, name);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Removes an edge attribute.
 *
 * @param graph A graph object.
 * @param name The name of the attribute to remove.
 * @return Error code.
 *
 * Time complexity: <code>O(|E|)</code>, assuming that the graph has
 * <code>O(1)</code> edge attributes. <code>|E|</code> is the number
 * of edges.
 */

int igraph_remove_edge_attribute(igraph_t *graph, const char *name) {
  igraph_attribute_list_remove(&graph->eal, name);  
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the value of an edge attribute for a single edge
 * 
 * It is easy to get the id of an edge by using an edge iterator.
 * @param graph The graph object.
 * @param name The name of the edge attribute.
 * @param e The id of the edge of which the attribute is requested.
 * @param value Pointer to a real number, the result will be stored
 *        here.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>, assuming that the graph has
 * <code>O(1)</code> edge attributes installed.
 */

int igraph_get_edge_attribute(igraph_t *graph, const char *name, 
			      long int e, real_t *value) {
  igraph_attribute_list_get(&graph->eal, name, e, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Set the value of an edge attribute for a single edge.
 * 
 * @param graph The graph object.
 * @param name Name of the edge attribute.
 * @param e The id of the edge of which the attribute is set.
 * @param value The new value of the attribute.
 * @return Error code.
 *
 * Time complexity: <code>O(1)</code>, assuming that the graph has
 * <code>O(1)</code> edge attributes installed.
 */

int igraph_set_edge_attribute(igraph_t *graph, const char *name, 
			      long int e, real_t value) {
  igraph_attribute_list_set(&graph->eal, name, e, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Query the value of an edge attribute for many edges.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to get.
 * @param e Vector with the edge ids of the edges of which the
 *        attribute will be returned.
 * @param value Pointer to an initialized vector, the result will be
 *        stored here. It will be resized if needed.
 * @return Error code.
 * 
 * Time complexity: <code>O(|e|)</code>, the number of queried
 * edges, assuming the graph has <code>O(1)</code> edge
 * attributes. 
 */

int igraph_get_edge_attributes(igraph_t *graph, const char *name, 
			       vector_t *e, vector_t *value) {
  igraph_attribute_list_gets(&graph->eal, name, e, value);
  return 0;  
}

/**
 * \ingroup attributes
 * \brief Set the value of an edge attribute for many edges.
 *
 * @param graph The graph object.
 * @param name The name of the attribute to set.
 * @param e Vector with the edge ids of the edges of which the
 *        attribute will be set.
 * @param value The new value(s) of the attribute. This vector may be of
 *        different length than <code>e</code>, if it is shorter it
 *        will be recycled (ie. after the last element the first one
 *        is used again), if it is longer the unneeded values are
 *        ignored. Thus it is easy to set an attribute to a single
 *        constant value for many edges, just give a vector of
 *        length 1 here.
 * @return Error code.
 * 
 * Time complexity: <code>O(|v|)</code>, the number of affected
 * edges, assuming the graph has <code>O(1)</code> edge
 * attributes. 
 */

int igraph_set_edge_attributes(igraph_t *graph, const char *name, 
			       vector_t *e, vector_t *value) {
  igraph_attribute_list_sets(&graph->eal, name, e, value);
  return 0;
}

/**
 * \ingroup attributes
 * \brief Queries the list of installed edge attributes.
 * 
 * @param graph A graph object.
 * @param l This string array will contain the names of the 
 *        attributes. It should be initialized and will be resized.
 * @return Error code.
 * 
 * Time complexity: <code>O(1)</code>. (Assuming the number of edge
 * attributes of <code>graph</code> is <code>O(1)</code>.)
 */

int igraph_list_edge_attributes(igraph_t *graph, igraph_strarray_t *l) {
  igraph_attribute_list_list(&graph->eal, l);
  return 0;
}


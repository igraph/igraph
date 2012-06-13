# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""
Implementation of `igraph.Graph.Formula()`

You should use this module directly only if you have a very strong reason
to do so. In almost all cases, you are better off with calling
`igraph.Graph.Formula()`.
"""

from cStringIO import StringIO
from igraph.datatypes import UniqueIdGenerator

import tokenize
import token

__all__ = ["construct_graph_from_formula"]

__license__ = u"""\
Copyright (C) 2006-2012  Tamás Nepusz <ntamas@gmail.com>
Pázmány Péter sétány 1/a, 1117 Budapest, Hungary

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
"""

def generate_edges(formula):
    """Parses an edge specification from the head of the given
    formula part and yields the following:
    
      - startpoint(s) of the edge by vertex names
      - endpoint(s) of the edge by names or an empty list if the vertices are isolated
      - a pair of bools to denote whether we had arrowheads at the start and end vertices 
    """
    if formula == "":
        yield [], [""], [False, False]
        return

    name_tokens = set([token.NAME, token.NUMBER, token.STRING])
    edge_chars = "<>-+"
    start_names, end_names, arrowheads = [], [], [False, False]
    parsing_vertices = True

    # Tokenize the formula
    token_gen = tokenize.generate_tokens(StringIO(formula).next)
    for token_info in token_gen:
        # Do the state transitions
        token_type, tok, _, _, _ = token_info
        if parsing_vertices:
            if all(ch in edge_chars for ch in tok) and token_type == token.OP:
                parsing_vertices = False
                # Check the edge we currently have
                if start_names and end_names:
                    # We have a whole edge
                    yield start_names, end_names, arrowheads
                start_names, end_names = end_names, []
                arrowheads = [False, False]
        else:
            if any(ch not in edge_chars for ch in tok):
                parsing_vertices = True

        if parsing_vertices:
            # We are parsing vertex names at the moment
            if token_type in name_tokens:
                # We found a vertex name
                if token_type == token.STRING:
                    end_names.append(eval(tok))
                else:
                    end_names.append(str(tok))
            elif tok == ":" and token_type == token.OP:
                # Separating semicolon between vertex names, we just go on
                continue
            elif token_type == token.ENDMARKER:
                # End markers are fine
                pass
            else:
                msg = "invalid token found in edge specification: %s" % formula
                raise SyntaxError(msg)
        else:
            # We are parsing an edge operator
            if "<" in tok:
                if ">" in tok:
                    arrowheads = [True, True]
                else:
                    arrowheads[0] = True
            elif ">" in tok:
                arrowheads[1] = True
            elif "+" in tok:
                if tok[0] == "+":
                    arrowheads[0] = True
                if len(tok) > 1 and tok[-1] == "+":
                    arrowheads[1] = True

    # The final edge
    yield start_names, end_names, arrowheads


def construct_graph_from_formula(cls, formula = None, attr = "name",
        simplify = True):
    """Graph.Formula(formula = None, attr = "name", simplify = True)
    
    Generates a graph from a graph formula

    A graph formula is a simple string representation of a graph.
    It is very handy for creating small graphs quickly. The string
    consists of vertex names separated by edge operators. An edge
    operator is a sequence of dashes (C{-}) that may or may not
    start with an arrowhead (C{<} at the beginning of the sequence
    or C{>} at the end of the sequence). The edge operators can
    be arbitrarily long, i.e., you may use as many dashes to draw
    them as you like. This makes a total of four different edge
    operators:

      - C{-----} makes an undirected edge
      - C{<----} makes a directed edge pointing from the vertex
        on the right hand side of the operator to the vertex on
        the left hand side
      - C{---->} is the opposite of C{<----}
      - C{<--->} creates a mutual directed edge pair between
        the two vertices

    If you only use the undirected edge operator (C{-----}),
    the graph will be undirected. Otherwise it will be directed.
    Vertex names used in the formula will be assigned to the
    C{name} vertex attribute of the graph.

    Some simple examples:

      >>> from igraph import Graph
      >>> print Graph.Formula()           # empty graph
      IGRAPH UN-- 0 0 --
      + attr: name (v)
      >>> g = Graph.Formula("A-B")        # undirected graph
      >>> g.vs["name"]
      ['A', 'B']
      >>> print g
      IGRAPH UN-- 2 1 --
      + attr: name (v)
      + edges (vertex names):
      A--B
      >>> g.get_edgelist()
      [(0, 1)]
      >>> g2 = Graph.Formula("A-----------B")
      >>> g2.isomorphic(g)
      True
      >>> g = Graph.Formula("A  --->  B") # directed graph
      >>> g.vs["name"]
      ['A', 'B']
      >>> print g
      IGRAPH DN-- 2 1 --
      + attr: name (v)
      + edges (vertex names):
      A->B
      
    If you have may disconnected componnets, you can separate them
    with commas. You can also specify isolated vertices:

      >>> g = Graph.Formula("A--B, C--D, E--F, G--H, I, J, K")
      >>> print ", ".join(g.vs["name"])
      A, B, C, D, E, F, G, H, I, J, K
      >>> g.clusters().membership
      [0, 0, 1, 1, 2, 2, 3, 3, 4, 5, 6]

    The colon (C{:}) operator can be used to specify vertex sets.
    If an edge operator connects two vertex sets, then every vertex
    from the first vertex set will be connected to every vertex in
    the second set:

      >>> g = Graph.Formula("A:B:C:D --- E:F:G")
      >>> g.isomorphic(Graph.Full_Bipartite(4, 3))
      True

    Note that you have to quote vertex names if they include spaces
    or special characters:

      >>> g = Graph.Formula('"this is" +- "a silly" -+ "graph here"')
      >>> g.vs["name"]
      ['this is', 'a silly', 'graph here']

    @param formula: the formula itself
    @param attr: name of the vertex attribute where the vertex names
                 will be stored
    @param simplify: whether the simplify the constructed graph
    @return: the constructed graph:
    """
    
    # If we have no formula, return an empty graph
    if formula is None:
        return cls(0, vertex_attrs = {attr: []})

    vertex_ids, edges, directed = UniqueIdGenerator(), [], False
    # Loop over each part in the formula
    for part in formula.split(","):
        # Drop newlines from the part
        part = part.strip().replace("\n", "").replace("\t", "")
        # Parse the first vertex specification from the formula
        for start_names, end_names, arrowheads in generate_edges(part):
            start_ids = [vertex_ids[name] for name in start_names]
            end_ids   = [vertex_ids[name] for name in end_names]
            if not arrowheads[0] and not arrowheads[1]:
                # This is an undirected edge. Do we have a directed graph?
                if not directed:
                    # Nope, add the edge
                    edges.extend((id1, id2) for id1 in start_ids \
                                 for id2 in end_ids)
            else:
                # This is a directed edge
                directed = True
                if arrowheads[1]:
                    edges.extend((id1, id2) for id1 in start_ids \
                                 for id2 in end_ids)
                if arrowheads[0]:
                    edges.extend((id2, id1) for id1 in start_ids \
                                 for id2 in end_ids)

    # Grab the vertex names into a list
    vertex_attrs = {}
    vertex_attrs[attr] = vertex_ids.values()
    # Construct and return the graph
    result = cls(len(vertex_ids), edges, directed, vertex_attrs=vertex_attrs)
    if simplify:
        result.simplify()
    return result

#!/usr/bin/env python
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

__all__ = ["Formula"]

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
    for token_type, tok, (_, start_char), (_, end_char), _ in token_gen:
        # Do the state transitions
        if parsing_vertices:
            if tok in edge_chars and token_type == token.OP:
                parsing_vertices = False
                # Check the edge we currently have
                if start_names and end_names:
                    # We have a whole edge
                    yield start_names, end_names, arrowheads
                start_names, end_names, arrowheads = end_names, [], [False, False]
        else:
            if tok not in edge_chars:
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
                raise SyntaxError, "invalid token found in edge specification: %s" % formula
        else:
            # We are parsing an edge operator
            if   tok == "<":  arrowheads[0] = True
            elif tok == ">":  arrowheads[1] = True
            elif tok == "<>": arrowheads = [True, True]
            elif tok == "+":
                pass   # TODO!

    # The final edge
    yield start_names, end_names, arrowheads


def Formula(klass, formula = None, attr = "name"):
    """Generates a graph from a graph formula

    TODO"""
    
    # If we have no formula, return an empty graph
    if formula is None: return klass(0, vertex_attrs = {attr: []})

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
                    edges.extend([(id1, id2) for id1 in start_ids for id2 in end_ids])
            else:
                # This is a directed edge
                directed = True
                if arrowheads[1]:
                    edges.extend([(id1, id2) for id1 in start_ids for id2 in end_ids])
                if arrowheads[0]:
                    edges.extend([(id2, id1) for id1 in start_ids for id2 in end_ids])

    # Grab the vertex names into a list
    names = sorted(((v, k) for k, v in vertex_ids._ids.iteritems()))
    names = [k for _, k in names]
    # Construct and return the graph
    return klass(len(names), edges, directed, vertex_attrs={attr: names})

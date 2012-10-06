# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""Classes representing cuts and flows on graphs."""

from igraph.clustering import VertexClustering

__license__ = """\
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

class Cut(VertexClustering):
    """A cut of a given graph.

    This is a simple class used to represent cuts returned by
    L{Graph.mincut()}, L{Graph.all_st_cuts()} and other functions
    that calculate cuts.
    
    A cut is a special vertex clustering with only two clusters.
    Besides the usual L{VertexClustering} methods, it also has the
    following attributes:

      - C{value} - the value (capacity) of the cut. It is equal to
        the number of edges if there are no capacities on the
        edges.

      - C{partition} - vertex IDs in the parts created
        after removing edges in the cut

      - C{cut} - edge IDs in the cut

      - C{es} - an edge selector restricted to the edges
        in the cut.

    You can use indexing on this object to obtain lists of vertex IDs
    for both sides of the partition.

    This class is usually not instantiated directly, everything
    is taken care of by the functions that return cuts.

    Examples:

      >>> from igraph import Graph
      >>> g = Graph.Ring(20)
      >>> mc = g.mincut()
      >>> print mc.value
      2.0
      >>> print min(map(len, mc))
      1
      >>> mc.es["color"] = "red"
    """

    # pylint: disable-msg=R0913
    def __init__(self, graph, value=None, cut=None, partition=None,
            partition2=None):
        """Initializes the cut.

        This should not be called directly, everything is taken care of by
        the functions that return cuts.
        """
        # Input validation
        if partition is None or cut is None:
            raise ValueError("partition and cut must be given")

        # Set up a membership vector, initialize parent class
        membership = [1] * graph.vcount()
        for vid in partition:
            membership[vid] = 0
        VertexClustering.__init__(self, graph, membership)

        if value is None:
            # Value of the cut not given, count the number of edges
            value = len(cut)
        self._value = float(value)

        self._partition = sorted(partition)
        self._cut = cut

    def __repr__(self):
        return "%s(%r, %r, %r, %r)" % \
          (self.__class__.__name__, self._graph, \
           self._value, self._cut, self._partition)

    def __str__(self):
        return "Graph cut (%d edges, %d vs %d vertices, value=%.4f)" % \
          (len(self._cut), len(self._partition),
          self._graph.vcount() - len(self._partition), self._value)

    # pylint: disable-msg=C0103
    @property
    def es(self):
        """Returns an edge selector restricted to the cut"""
        return self._graph.es.select(self.cut)

    @property
    def partition(self):
        """Returns the vertex IDs partitioned according to the cut"""
        return list(self)

    @property
    def cut(self):
        """Returns the edge IDs in the cut"""
        return self._cut

    @property
    def value(self):
        """Returns the sum of edge capacities in the cut"""
        return self._value


class Flow(Cut):
    """A flow of a given graph.

    This is a simple class used to represent flows returned by
    L{Graph.maxflow}. It has the following attributes:

      - C{graph} - the graph on which this flow is defined

      - C{value} - the value (capacity) of the flow 

      - C{flow} - the flow values on each edge. For directed graphs,
        this is simply a list where element M{i} corresponds to the
        flow on edge M{i}. For undirected graphs, the direction of
        the flow is not constrained (since the edges are undirected),
        hence positive flow always means a flow from the smaller vertex
        ID to the larger, while negative flow means a flow from the
        larger vertex ID to the smaller.

      - C{cut} - edge IDs in the minimal cut corresponding to
        the flow.

      - C{partition} - vertex IDs in the parts created
        after removing edges in the cut

      - C{es} - an edge selector restricted to the edges
        in the cut.

    This class is usually not instantiated directly, everything
    is taken care of by L{Graph.maxflow}.

    Examples:

      >>> from igraph import Graph
      >>> g = Graph.Ring(20)
      >>> mf = g.maxflow(0, 10)
      >>> print mf.value
      2.0
      >>> mf.es["color"] = "red"
    """

    # pylint: disable-msg=R0913
    def __init__(self, graph, value, flow, cut, partition):
        """Initializes the flow.

        This should not be called directly, everything is
        taken care of by L{Graph.maxflow}.
        """
        super(Flow, self).__init__(graph, value, cut, partition)
        self._flow = flow

    def __repr__(self):
        return "%s(%r, %r, %r, %r, %r)" % \
          (self.__class__.__name__, self._graph, \
           self._value, self._flow, self._cut, self._partition)

    def __str__(self):
        return "Graph flow (%d edges, %d vs %d vertices, value=%.4f)" % \
          (len(self._cut), len(self._partition),
          self._graph.vcount() - len(self._partition), self._value)

    @property
    def flow(self):
        """Returns the flow values for each edge.
        
        For directed graphs, this is simply a list where element M{i}
        corresponds to the flow on edge M{i}. For undirected graphs, the
        direction of the flow is not constrained (since the edges are
        undirected), hence positive flow always means a flow from the smaller
        vertex ID to the larger, while negative flow means a flow from the
        larger vertex ID to the smaller.
        """
        return self._flow




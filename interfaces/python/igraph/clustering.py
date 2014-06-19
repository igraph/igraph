# vim:ts=4:sw=4:sts=4:et
# -*- coding: utf-8 -*-
"""Classes related to graph clustering.

@undocumented: _handle_mark_groups_arg_for_clustering, _prepare_community_comparison"""

__license__ = u"""
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

from copy import deepcopy
from itertools import izip
from math import pi
from cStringIO import StringIO

from igraph import community_to_membership
from igraph.compat import property
from igraph.configuration import Configuration
from igraph.datatypes import UniqueIdGenerator
from igraph.drawing.colors import ClusterColoringPalette
from igraph.statistics import Histogram
from igraph.summary import _get_wrapper_for_width
from igraph.utils import str_to_orientation

class Clustering(object):
    """Class representing a clustering of an arbitrary ordered set.
    
    This is now used as a base for L{VertexClustering}, but it might be
    useful for other purposes as well.
    
    Members of an individual cluster can be accessed by the C{[]} operator:
    
      >>> cl = Clustering([0,0,0,0,1,1,1,2,2,2,2])
      >>> cl[0]
      [0, 1, 2, 3]
    
    The membership vector can be accessed by the C{membership} property:

      >>> cl.membership
      [0, 0, 0, 0, 1, 1, 1, 2, 2, 2, 2]

    The number of clusters can be retrieved by the C{len} function:

      >>> len(cl)
      3

    You can iterate over the clustering object as if it were a regular list
    of clusters:

      >>> for cluster in cl:
      ...     print " ".join(str(idx) for idx in cluster)
      ...
      0 1 2 3
      4 5 6
      7 8 9 10

    If you need all the clusters at once as lists, you can simply convert
    the clustering object to a list:

      >>> cluster_list = list(cl)
      >>> print cluster_list
      [[0, 1, 2, 3], [4, 5, 6], [7, 8, 9, 10]]

    @undocumented: _formatted_cluster_iterator
    """

    def __init__(self, membership, params = None):
        """Constructor.

        @param membership: the membership list -- that is, the cluster
          index in which each element of the set belongs to.
        @param params: additional parameters to be stored in this
          object's dictionary."""
        self._membership = list(membership)
        if len(self._membership)>0:
            self._len = max(m for m in self._membership if m is not None)+1
        else:
            self._len = 0

        if params:
            self.__dict__.update(params)
    
    def __getitem__(self, idx):
        """Returns the members of the specified cluster.

        @param idx: the index of the cluster
        @return: the members of the specified cluster as a list
        @raise IndexError: if the index is out of bounds"""
        if idx < 0 or idx >= self._len:
            raise IndexError("cluster index out of range")
        return [i for i, e in enumerate(self._membership) if e == idx]

    def __iter__(self):
        """Iterates over the clusters in this clustering.

        This method will return a generator that generates the clusters
        one by one."""
        clusters = [[] for _ in xrange(self._len)]
        for idx, cluster in enumerate(self._membership):
            clusters[cluster].append(idx)
        return iter(clusters)

    def __len__(self):
        """Returns the number of clusters.

        @return: the number of clusters
        """
        return self._len

    def __str__(self):
        return self.summary(verbosity=1, width=78)

    def as_cover(self):
        """Returns a L{Cover} that contains the same clusters as this clustering."""
        return Cover(self._graph, self)

    def compare_to(self, other, *args, **kwds):
        """Compares this clustering to another one using some similarity or
        distance metric.

        This is a convenience method that simply calls L{compare_communities}
        with the two clusterings as arguments. Any extra positional or keyword
        argument is also forwarded to L{compare_communities}."""
        return compare_communities(self, other, *args, **kwds)

    @property
    def membership(self):
        """Returns the membership vector."""
        return self._membership[:]

    @property
    def n(self):
        """Returns the number of elements covered by this clustering."""
        return len(self._membership)

    def size(self, idx):
        """Returns the size of a given cluster.

        @param idx: the cluster in which we are interested.
        """
        return len(self[idx])

    def sizes(self, *args):
        """Returns the size of given clusters.
        
        The indices are given as positional arguments. If there are no
        positional arguments, the function will return the sizes of all clusters.
        """
        counts = [0] * len(self)
        for x in self._membership:
            counts[x] += 1

        if args:
            return [counts[idx] for idx in args]

        return counts
    
    def size_histogram(self, bin_width = 1):
        """Returns the histogram of cluster sizes.

        @param bin_width: the bin width of the histogram
        @return: a L{Histogram} object
        """
        return Histogram(bin_width, self.sizes())

    def summary(self, verbosity=0, width=None):
        """Returns the summary of the clustering.
        
        The summary includes the number of items and clusters, and also the
        list of members for each of the clusters if the verbosity is nonzero.
        
        @param verbosity: determines whether the cluster members should be
          printed. Zero verbosity prints the number of items and clusters only.
        @return: the summary of the clustering as a string.
        """
        out = StringIO()
        print >>out, "Clustering with %d elements and %d clusters" % \
                (len(self._membership), len(self))

        if verbosity < 1:
            return out.getvalue().strip()

        ndigits = len(str(len(self)))
        wrapper = _get_wrapper_for_width(width,
                subsequent_indent = " " * (ndigits+3))

        for idx, cluster in enumerate(self._formatted_cluster_iterator()):
            wrapper.initial_indent = "[%*d] " % (ndigits, idx)
            print >>out, "\n".join(wrapper.wrap(cluster))

        return out.getvalue().strip()

    def _formatted_cluster_iterator(self):
        """Iterates over the clusters and formats them into a string to be
        presented in the summary."""
        for cluster in self:
            yield ", ".join(str(member) for member in cluster)


class VertexClustering(Clustering):
    """The clustering of the vertex set of a graph.

    This class extends L{Clustering} by linking it to a specific L{Graph} object
    and by optionally storing the modularity score of the clustering.
    It also provides some handy methods like getting the subgraph corresponding
    to a cluster and such.

    @note: since this class is linked to a L{Graph}, destroying the graph by the
      C{del} operator does not free the memory occupied by the graph if there
      exists a L{VertexClustering} that references the L{Graph}.

    @undocumented: _formatted_cluster_iterator
    """

    # Allow None to be passed to __plot__ as the "palette" keyword argument
    _default_palette = None

    def __init__(self, graph, membership = None, modularity = None, \
                 params = None, modularity_params = None):
        """Creates a clustering object for a given graph.

        @param graph: the graph that will be associated to the clustering
        @param membership: the membership list. The length of the list must
          be equal to the number of vertices in the graph. If C{None}, every
          vertex is assumed to belong to the same cluster.
        @param modularity: the modularity score of the clustering. If C{None},
          it will be calculated when needed.
        @param params: additional parameters to be stored in this object.
        @param modularity_params: arguments that should be passed to
          L{Graph.modularity} when the modularity is (re)calculated. If the
          original graph was weighted, you should pass a dictionary
          containing a C{weight} key with the appropriate value here.
        """
        if membership is None:
            Clustering.__init__(self, [0]*graph.vcount(), params)
        else:
            if len(membership) != graph.vcount():
                raise ValueError("membership list has invalid length")
            Clustering.__init__(self, membership, params)

        self._graph = graph
        self._modularity = modularity
        if modularity_params is None:
            self._modularity_params = {}
        else:
            self._modularity_params = dict(modularity_params)

    # pylint: disable-msg=C0103
    @classmethod
    def FromAttribute(cls, graph, attribute, intervals=None, params=None):
        """Creates a vertex clustering based on the value of a vertex attribute.

        Vertices having the same attribute will correspond to the same cluster.

        @param graph: the graph on which we are working
        @param attribute: name of the attribute on which the clustering
            is based.
        @param intervals: for numeric attributes, you can either pass a single
            number or a list of numbers here. A single number means that the
            vertices will be put in bins of that width and vertices ending up
            in the same bin will be in the same cluster. A list of numbers
            specify the bin positions explicitly; e.g., C{[10, 20, 30]} means
            that there will be four categories: vertices with the attribute
            value less than 10, between 10 and 20, between 20 and 30 and over 30.
            Intervals are closed from the left and open from the right.
        @param params: additional parameters to be stored in this object.

        @return: a new VertexClustering object
        """
        from bisect import bisect

        def safeintdiv(x, y):
            """Safe integer division that handles None gracefully"""
            if x is None:
                return None
            return int(x / y)

        def safebisect(intervals, x):
            """Safe list bisection that handles None gracefully"""
            if x is None:
                return None
            return bisect(intervals, x)

        try:
            _ = iter(intervals)
            iterable = True
        except TypeError:
            iterable = False
        if intervals is None:
            vec = graph.vs[attribute]
        elif iterable:
            intervals = list(intervals)
            vec = [safebisect(intervals, x) for x in graph.vs[attribute]]
        else:
            intervals = float(intervals)
            vec = [safeintdiv(x, intervals) for x in graph.vs[attribute]]

        idgen = UniqueIdGenerator()
        idgen[None] = None
        vec = [idgen[i] for i in vec]
        return cls(graph, vec, None, params)

    def as_cover(self):
        """Returns a L{VertexCover} that contains the same clusters as this
        clustering."""
        return VertexCover(self._graph, self)

    def cluster_graph(self, combine_vertices=None, combine_edges=None):
        """Returns a graph where each cluster is contracted into a single
        vertex.

        In the resulting graph, vertex M{i} represents cluster M{i} in this
        clustering. Vertex M{i} and M{j} will be connected if there was
        at least one connected vertex pair M{(a, b)} in the original graph such
        that vertex M{a} was in cluster M{i} and vertex M{b} was in cluster
        M{j}.

        @param combine_vertices: specifies how to derive the attributes of
          the vertices in the new graph from the attributes of the old ones.
          See L{Graph.contract_vertices()} for more details.
        @param combine_edges: specifies how to derive the attributes of the
          edges in the new graph from the attributes of the old ones. See
          L{Graph.simplify()} for more details. If you specify C{False}
          here, edges will not be combined, and the number of edges between
          the vertices representing the original clusters will be equal to
          the number of edges between the members of those clusters in the
          original graph.

        @return: the new graph.
        """
        result = self.graph.copy()
        result.contract_vertices(self.membership, combine_vertices)
        if combine_edges != False:
            result.simplify(combine_edges=combine_edges)
        return result

    def crossing(self):
        """Returns a boolean vector where element M{i} is C{True} iff edge
        M{i} lies between clusters, C{False} otherwise."""
        membership = self.membership
        return [membership[v1] != membership[v2] \
                for v1, v2 in self.graph.get_edgelist()]

    @property
    def modularity(self):
        """Returns the modularity score"""
        if self._modularity is None:
            return self.recalculate_modularity()
        return self._modularity
    q = modularity

    @property
    def graph(self):
        """Returns the graph belonging to this object"""
        return self._graph

    def recalculate_modularity(self):
        """Recalculates the stored modularity value.

        This method must be called before querying the modularity score of the
        clustering through the class member C{modularity} or C{q} if the
        graph has been modified (edges have been added or removed) since the
        creation of the L{VertexClustering} object.
        
        @return: the new modularity score
        """
        self._modularity = self._graph.modularity(self._membership,
                **self._modularity_params)
        return self._modularity


    def subgraph(self, idx):
        """Get the subgraph belonging to a given cluster.

        @param idx: the cluster index
        @return: a copy of the subgraph
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the clustering was constructed.
        """
        return self._graph.subgraph(self[idx])


    def subgraphs(self):
        """Gets all the subgraphs belonging to each of the clusters.

        @return: a list containing copies of the subgraphs
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the clustering was constructed.
        """
        return [self._graph.subgraph(cl) for cl in self]


    def giant(self):
        """Returns the giant community of the clustered graph.

        The giant component a community for which no larger community exists.
        @note: there can be multiple giant communities, this method will return
          the copy of an arbitrary one if there are multiple giant communities.

        @return: a copy of the giant community.
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the clustering was constructed.
        """
        ss = self.sizes()
        max_size = max(ss)
        return self.subgraph(ss.index(max_size))

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Plots the clustering to the given Cairo context in the given
        bounding box.

        This is done by calling L{Graph.__plot__()} with the same arguments, but
        coloring the graph vertices according to the current clustering (unless
        overridden by the C{vertex_color} argument explicitly).

        This method understands all the positional and keyword arguments that
        are understood by L{Graph.__plot__()}, only the differences will be
        highlighted here:

          - C{mark_groups}: whether to highlight some of the vertex groups by
            colored polygons. Besides the values accepted by L{Graph.__plot__}
            (i.e., a dict mapping colors to vertex indices, a list containing
            lists of vertex indices, or C{False}), the following are also
            accepted:

              - C{True}: all the groups will be highlighted, the colors matching
                the corresponding color indices from the current palette
                (see the C{palette} keyword argument of L{Graph.__plot__}.

              - A dict mapping cluster indices or tuples of vertex indices to
                color names.  The given clusters or vertex groups will be
                highlighted by the given colors.

              - A list of cluster indices. This is equivalent to passing a
                dict mapping numeric color indices from the current palette
                to cluster indices; therefore, the cluster referred to by element
                I{i} of the list will be highlighted by color I{i} from the
                palette.

            The value of the C{plotting.mark_groups} configuration key is also
            taken into account here; if that configuration key is C{True} and
            C{mark_groups} is not given explicitly, it will automatically be set
            to C{True}.

            In place of lists of vertex indices, you may also use L{VertexSeq}
            instances.

            In place of color names, you may also use color indices into the
            current palette. C{None} as a color name will mean that the
            corresponding group is ignored.

          - C{palette}: the palette used to resolve numeric color indices to RGBA
            values. By default, this is an instance of L{ClusterColoringPalette}.

        @see: L{Graph.__plot__()} for more supported keyword arguments.
        """
        if "edge_color" not in kwds and "color" not in self.graph.edge_attributes():
            # Set up a default edge coloring based on internal vs external edges
            colors = ["grey20", "grey80"]
            kwds["edge_color"] = [colors[is_crossing]
                                  for is_crossing in self.crossing()]

        if palette is None:
            palette = ClusterColoringPalette(len(self))

        if "mark_groups" not in kwds:
            if Configuration.instance()["plotting.mark_groups"]:
                kwds["mark_groups"] = (
                    (group, color) for color, group in enumerate(self)
                )
        else:
            kwds["mark_groups"] = _handle_mark_groups_arg_for_clustering(
                    kwds["mark_groups"], self)

        if "vertex_color" not in kwds:
            kwds["vertex_color"] = self.membership

        return self._graph.__plot__(context, bbox, palette, *args, **kwds)

    def _formatted_cluster_iterator(self):
        """Iterates over the clusters and formats them into a string to be
        presented in the summary."""
        if self._graph.is_named():
            names = self._graph.vs["name"]
            for cluster in self:
                yield ", ".join(str(names[member]) for member in cluster)
        else:
            for cluster in self:
                yield ", ".join(str(member) for member in cluster)


###############################################################################

class Dendrogram(object):
    """The hierarchical clustering (dendrogram) of some dataset.

    A hierarchical clustering means that we know not only the way the
    elements are separated into groups, but also the exact history of
    how individual elements were joined into larger subgroups.

    This class internally represents the hierarchy by a matrix with n rows
    and 2 columns -- or more precisely, a list of lists of size 2. This is
    exactly the same as the original format used by C{igraph}'s C core.
    The M{i}th row of the matrix contains the indices of the two clusters
    being joined in time step M{i}. The joint group will be represented by
    the ID M{n+i}, with M{i} starting from one. The ID of the joint group
    will be referenced in the upcoming steps instead of any of its individual
    members. So, IDs less than or equal to M{n} (where M{n} is the number of
    rows in the matrix) mean the original members of the dataset (with ID
    from 0 to M{n}), while IDs up from M{n+1} mean joint groups. As an
    example, take a look at the dendrogram and the internal representation of
    a given clustering of five nodes::

      0 -+
         |
      1 -+-+
           |
      2 ---+-+        <====>   [[0, 1], [3, 4], [2, 5], [6, 7]]
             |
      3 -+   |
         |   |
      4 -+---+---

    @undocumented: _item_box_size, _plot_item, _traverse_inorder
    """

    def __init__(self, merges):
        """Creates a hierarchical clustering.

        @param merges: the merge history either in matrix or tuple format"""
        self._merges = [tuple(pair) for pair in merges]
        self._nmerges = len(self._merges)
        if self._nmerges:
            self._nitems = max(self._merges[-1])-self._nmerges+2
        else:
            self._nitems = 0
        self._names = None

    @staticmethod
    def _convert_matrix_to_tuple_repr(merges, n=None):
        """Converts the matrix representation of a clustering to a tuple
        representation.
        
        @param merges: the matrix representation of the clustering
        @return: the tuple representation of the clustering
        """
        if n is None:
            n = len(merges)+1
        tuple_repr = range(n)
        idxs = range(n)
        for rowidx, row in enumerate(merges):
            i, j = row
            try:
                idxi, idxj = idxs[i], idxs[j]
                tuple_repr[idxi] = (tuple_repr[idxi], tuple_repr[idxj])
                tuple_repr[idxj] = None
            except IndexError:
                raise ValueError("malformed matrix, subgroup referenced "+
                                 "before being created in step %d" % rowidx)
            idxs.append(j)
        return [x for x in tuple_repr if x is not None]

    def _traverse_inorder(self):
        """Conducts an inorder traversal of the merge tree.

        The inorder traversal returns the nodes on the last level in the order
        they should be drawn so that no edges cross each other.

        @return: the result of the inorder traversal in a list."""
        stack = [self._merges[-1]]
        result = []
        while len(stack)>0:
            last = stack[-1]
            if len(last) == 0:
                stack.pop()
                continue
            elif len(last) == 1:       # Right child
                stack[-1] = ()
                last = last[0]
            else:                      # Left child
                stack[-1] = (last[1],)
                last = last[0]
            if last < self._nitems: # This will be a regular node
                result.append(last)
            else:        # This is a merge node, proceed towards left
                stack.append(self._merges[last-self._nitems])

        return result

    def __str__(self):
        return self.summary(verbosity=1)

    def format(self, format="newick"):
        """Formats the dendrogram in a foreign format.

        Currently only the Newick format is supported.
        
        Example:
            
            >>> d = Dendrogram([(2, 3), (0, 1), (4, 5)])
            >>> d.format()
            '((2,3)4,(0,1)5)6;'
            >>> d.names = list("ABCDEFG")
            >>> d.format()
            '((C,D)E,(A,B)F)G;'
        """
        if format == "newick":
            n = self._nitems + self._nmerges
            if self._names is None:
                nodes = range(n)
            else:
                nodes = list(self._names)
            if len(nodes) < n:
                nodes.extend("" for _ in xrange(n - len(nodes)))
            for k, (i, j) in enumerate(self._merges, self._nitems):
                nodes[k] = "(%s,%s)%s" % (nodes[i], nodes[j], nodes[k])
                nodes[i] = nodes[j] = None
            return nodes[-1] + ";"
        raise ValueError("unsupported format: %r" % format)

    def summary(self, verbosity=0, max_leaf_count=40):
        """Returns the summary of the dendrogram.
        
        The summary includes the number of leafs and branches, and also an
        ASCII art representation of the dendrogram unless it is too large.
        
        @param verbosity: determines whether the ASCII representation of the
          dendrogram should be printed. Zero verbosity prints only the number
          of leafs and branches.
        @param max_leaf_count: the maximal number of leafs to print in the
          ASCII representation. If the dendrogram has more leafs than this
          limit, the ASCII representation will not be printed even if the
          verbosity is larger than or equal to 1.
        @return: the summary of the dendrogram as a string.
        """
        from array import array

        out = StringIO()
        print >>out, "Dendrogram, %d elements, %d merges" % \
                (self._nitems, self._nmerges)

        if self._nitems == 0 or verbosity < 1 or self._nitems > max_leaf_count:
            return out.getvalue().strip()
            
        print >>out

        positions = [None] * self._nitems
        inorder = self._traverse_inorder()
        distance = 2
        level_distance = 2
        nextp = 0
        for idx, element in enumerate(inorder):
            positions[element] = nextp
            inorder[idx] = str(element)
            nextp += max(distance, len(inorder[idx])+1)

        width = max(positions)+1

        # Print the nodes on the lowest level
        print >>out, (" " * (distance-1)).join(inorder)
        midx = 0
        max_community_idx = self._nitems
        while midx < self._nmerges:
            char_array = array("c", " "*width)
            for position in positions:
                if position >= 0:
                    char_array[position] = "|"
            char_str = char_array.tostring()
            for _ in xrange(level_distance-1):
                print >>out, char_str # Print the lines
            
            cidx_incr = 0
            while midx < self._nmerges:
                id1, id2 = self._merges[midx]
                if id1 >= max_community_idx or id2 >= max_community_idx:
                    break
                midx += 1

                pos1, pos2 = positions[id1], positions[id2]
                positions[id1], positions[id2] = -1, -1
                positions.append((pos1+pos2)/2)

                dashes = "-" * (pos2 - pos1 - 1)
                char_array[pos1:(pos2+1)] = array("c", "`%s'" % dashes)

                cidx_incr += 1
            
            max_community_idx += cidx_incr

            print >>out, char_array.tostring()

        return out.getvalue().strip()

    def _item_box_size(self, context, horiz, idx):
        """Calculates the amount of space needed for drawing an
        individual vertex at the bottom of the dendrogram."""
        if self._names is None or self._names[idx] is None:
            x_bearing, _, _, height, x_advance, _ = context.text_extents("")
        else:
            x_bearing, _, _, height, x_advance, _ = context.text_extents(str(self._names[idx]))

        if horiz:
            return x_advance - x_bearing, height
        return height, x_advance - x_bearing

    # pylint: disable-msg=R0913
    def _plot_item(self, context, horiz, idx, x, y):
        """Plots a dendrogram item to the given Cairo context

        @param context: the Cairo context we are plotting on
        @param horiz: whether the dendrogram is horizontally oriented
        @param idx: the index of the item
        @param x: the X position of the item
        @param y: the Y position of the item
        """
        if self._names is None or self._names[idx] is None:
            return

        height = self._item_box_size(context, True, idx)[1]
        if horiz:
            context.move_to(x, y+height)
            context.show_text(str(self._names[idx]))
        else:
            context.save()
            context.translate(x, y)
            context.rotate(-pi/2.)
            context.move_to(0, height)
            context.show_text(str(self._names[idx]))
            context.restore()

    # pylint: disable-msg=C0103,W0613
    # W0613 = unused argument 'palette'
    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Draws the dendrogram on the given Cairo context

        Supported keyword arguments are:

          - C{orientation}: the orientation of the dendrogram. Must be one of
            the following values: C{left-right}, C{bottom-top}, C{right-left}
            or C{top-bottom}. Individual elements are always placed at the
            former edge and merges are performed towards the latter edge.
            Possible aliases: C{horizontal} = C{left-right},
            C{vertical} = C{bottom-top}, C{lr} = C{left-right},
            C{rl} = C{right-left}, C{tb} = C{top-bottom}, C{bt} = C{bottom-top}.
            The default is C{left-right}.

        """
        from igraph.layout import Layout

        if self._names is None:
            self._names = [str(x) for x in xrange(self._nitems)]

        orientation = str_to_orientation(kwds.get("orientation", "lr"),
                reversed_vertical=True)
        horiz = orientation in ("lr", "rl")

        # Get the font height
        font_height = context.font_extents()[2]

        # Calculate space needed for individual items at the
        # bottom of the dendrogram
        item_boxes = [self._item_box_size(context, horiz, idx) \
          for idx in xrange(self._nitems)]

        # Small correction for cases when the right edge of the labels is
        # aligned with the tips of the dendrogram branches
        ygap = 2 if orientation == "bt" else 0
        xgap = 2 if orientation == "lr" else 0
        item_boxes = [(x+xgap, y+ygap) for x, y in item_boxes]

        # Calculate coordinates
        layout = Layout([(0, 0)] * self._nitems, dim=2)
        inorder = self._traverse_inorder()
        if not horiz:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                layout[element] = (x, 0)
                x += max(font_height, item_boxes[element][0])

            for id1, id2 in self._merges:
                y += 1
                layout.append(((layout[id1][0]+layout[id2][0])/2., y))

            # Mirror or rotate the layout if necessary
            if orientation == "bt":
                layout.mirror(1)
        else:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                layout[element] = (0, y)
                y += max(font_height, item_boxes[element][1])

            for id1, id2 in self._merges:
                x += 1
                layout.append((x, (layout[id1][1]+layout[id2][1])/2.))

            # Mirror or rotate the layout if necessary
            if orientation == "rl":
                layout.mirror(0)
        
        # Rescale layout to the bounding box
        maxw = max(e[0] for e in item_boxes)
        maxh = max(e[1] for e in item_boxes)

        # w, h: width and height of the area containing the dendrogram
        # tree without the items.
        # delta_x, delta_y: displacement of the dendrogram tree
        width, height = float(bbox.width), float(bbox.height)
        delta_x, delta_y = 0, 0
        if horiz:
            width -= maxw
            if orientation == "lr":
                delta_x = maxw
        else:
            height -= maxh
            if orientation == "tb":
                delta_y = maxh

        if horiz:
            delta_y += font_height / 2.
        else:
            delta_x += font_height / 2.
        layout.fit_into((delta_x, delta_y, width - delta_x, height - delta_y),
                        keep_aspect_ratio=False)

        context.save()

        context.translate(bbox.left, bbox.top)
        context.set_source_rgb(0., 0., 0.)
        context.set_line_width(1)
        
        # Draw items
        if horiz:
            sgn = 0 if orientation == "rl" else -1
            for idx in xrange(self._nitems):
                x = layout[idx][0] + sgn * item_boxes[idx][0]
                y = layout[idx][1] - item_boxes[idx][1]/2.
                self._plot_item(context, horiz, idx, x, y)
        else:
            sgn = 1 if orientation == "bt" else 0
            for idx in xrange(self._nitems):
                x = layout[idx][0] - item_boxes[idx][0]/2.
                y = layout[idx][1] + sgn * item_boxes[idx][1]
                self._plot_item(context, horiz, idx, x, y)

        # Draw dendrogram lines
        if not horiz:
            for idx, (id1, id2) in enumerate(self._merges):
                x0, y0 = layout[id1]
                x1, y1 = layout[id2]
                x2, y2 = layout[idx + self._nitems]
                context.move_to(x0, y0)
                context.line_to(x0, y2)
                context.line_to(x1, y2)
                context.line_to(x1, y1)
                context.stroke()
        else:
            for idx, (id1, id2) in enumerate(self._merges):
                x0, y0 = layout[id1]
                x1, y1 = layout[id2]
                x2, y2 = layout[idx + self._nitems]
                context.move_to(x0, y0)
                context.line_to(x2, y0)
                context.line_to(x2, y1)
                context.line_to(x1, y1)
                context.stroke()

        context.restore()

    @property
    def merges(self):
        """Returns the performed merges in matrix format"""
        return deepcopy(self._merges)

    @property
    def names(self):
        """Returns the names of the nodes in the dendrogram"""
        return self._names

    @names.setter
    def names(self, items):
        """Sets the names of the nodes in the dendrogram"""
        if items is None:
            self._names = None
            return

        items = list(items)
        if len(items) < self._nitems:
            raise ValueError("must specify at least %d names" % self._nitems)

        n = self._nitems + self._nmerges
        self._names = items[:n]
        if len(self._names) < n:
            self._names.extend("" for _ in xrange(n-len(self._names)))


class VertexDendrogram(Dendrogram):
    """The dendrogram resulting from the hierarchical clustering of the
    vertex set of a graph."""

    def __init__(self, graph, merges, optimal_count = None, params = None,
            modularity_params = None):
        """Creates a dendrogram object for a given graph.

        @param graph: the graph that will be associated to the clustering
        @param merges: the merges performed given in matrix form.
        @param optimal_count: the optimal number of clusters where the
          dendrogram should be cut. This is a hint usually provided by the
          clustering algorithm that produces the dendrogram. C{None} means
          that such a hint is not available; the optimal count will then be
          selected based on the modularity in such a case.
        @param params: additional parameters to be stored in this object.
        @param modularity_params: arguments that should be passed to
          L{Graph.modularity} when the modularity is (re)calculated. If the
          original graph was weighted, you should pass a dictionary
          containing a C{weight} key with the appropriate value here.
        """
        Dendrogram.__init__(self, merges)
        self._graph = graph
        self._optimal_count = optimal_count
        if modularity_params is None:
            self._modularity_params = {}
        else:
            self._modularity_params = dict(modularity_params)

    def as_clustering(self, n=None):
        """Cuts the dendrogram at the given level and returns a corresponding
        L{VertexClustering} object.

        @param n: the desired number of clusters. Merges are replayed from the
          beginning until the membership vector has exactly M{n} distinct elements
          or until there are no more recorded merges, whichever happens first.
          If C{None}, the optimal count hint given by the clustering algorithm
          will be used If the optimal count was not given either, it will be
          calculated by selecting the level where the modularity is maximal.
        @return: a new L{VertexClustering} object.
        """
        if n is None:
            n = self.optimal_count
        num_elts = self._graph.vcount()
        idgen = UniqueIdGenerator()
        membership = community_to_membership(self._merges, num_elts, \
                                             num_elts - n)
        membership = [idgen[m] for m in membership]
        return VertexClustering(self._graph, membership,
                modularity_params=self._modularity_params)

    @property
    def optimal_count(self):
        """Returns the optimal number of clusters for this dendrogram.

        If an optimal count hint was given at construction time, this
        property simply returns the hint. If such a count was not given,
        this method calculates the optimal number of clusters by maximizing
        the modularity along all the possible cuts in the dendrogram.
        """
        if self._optimal_count is not None:
            return self._optimal_count

        n = self._graph.vcount()
        max_q, optimal_count = 0, 1
        for step in xrange(min(n-1, len(self._merges))):
            membs = community_to_membership(self._merges, n, step)
            q = self._graph.modularity(membs, **self._modularity_params)
            if q > max_q:
                optimal_count = n-step
                max_q = q
        self._optimal_count = optimal_count
        return self._optimal_count

    @optimal_count.setter
    def optimal_count(self, value):
        self._optimal_count = max(int(value), 1)

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Draws the vertex dendrogram on the given Cairo context

        See L{Dendrogram.__plot__} for the list of supported keyword
        arguments."""
        from igraph.drawing.metamagic import AttributeCollectorBase

        class VisualVertexBuilder(AttributeCollectorBase):
            _kwds_prefix = "vertex_"
            label = None

        builder = VisualVertexBuilder(self._graph.vs, kwds)
        self._names = [vertex.label for vertex in builder]
        self._names = [name if name is not None else str(idx)
                       for idx, name in enumerate(self._names)]
        result = Dendrogram.__plot__(self, context, bbox, palette, \
                *args, **kwds)
        del self._names

        return result

###############################################################################

class Cover(object):
    """Class representing a cover of an arbitrary ordered set.

    Covers are similar to clusterings, but each element of the set may
    belong to more than one cluster in a cover, and elements not belonging
    to any cluster are also allowed.

    L{Cover} instances provide a similar API as L{Clustering} instances;
    for instance, iterating over a L{Cover} will iterate over the clusters
    just like with a regular L{Clustering} instance. However, they are not
    derived from each other or from a common superclass, and there might
    be functions that exist only in one of them or the other.

    Clusters of an individual cover can be accessed by the C{[]} operator:
    
      >>> cl = Cover([[0,1,2,3], [2,3,4], [0,1,6]])
      >>> cl[0]
      [0, 1, 2, 3]
    
    The membership vector can be accessed by the C{membership} property.
    Note that contrary to L{Clustering} instances, the membership vector
    will contain lists that contain the cluster indices each item belongs
    to:

      >>> cl.membership
      [[0, 2], [0, 2], [0, 1], [0, 1], [1], [], [2]]

    The number of clusters can be retrieved by the C{len} function:

      >>> len(cl)
      3

    You can iterate over the cover as if it were a regular list of
    clusters:

      >>> for cluster in cl:
      ...     print " ".join(str(idx) for idx in cluster)
      ...
      0 1 2 3
      2 3 4
      0 1 6

    If you need all the clusters at once as lists, you can simply convert
    the cover to a list:

      >>> cluster_list = list(cl)
      >>> print cluster_list
      [[0, 1, 2, 3], [2, 3, 4], [0, 1, 6]]

    L{Clustering} objects can readily be converted to L{Cover} objects
    using the constructor:

      >>> clustering = Clustering([0, 0, 0, 0, 1, 1, 1, 2, 2, 2])
      >>> cover = Cover(clustering)
      >>> list(clustering) == list(cover)
      True

    @undocumented: _formatted_cluster_iterator
    """

    def __init__(self, clusters, n=0):
        """Constructs a cover with the given clusters.

        @param clusters: the clusters in this cover, as a list or iterable.
          Each cluster is specified by a list or tuple that contains the
          IDs of the items in this cluster. IDs start from zero.

        @param n: the total number of elements in the set that is covered
          by this cover. If it is less than the number of unique elements
          found in all the clusters, we will simply use the number of unique
          elements, so it is safe to leave this at zero. You only have to
          specify this parameter if there are some elements that are covered
          by none of the clusters.
        """

        self._clusters = [list(cluster) for cluster in clusters]
        try:
            self._n = max(max(cluster)+1 for cluster in self._clusters if cluster)
        except ValueError:
            self._n = 0
        self._n = max(n, self._n)

    def __getitem__(self, index):
        """Returns the cluster with the given index."""
        return self._clusters[index]

    def __iter__(self):
        """Iterates over the clusters in this cover."""
        return iter(self._clusters)

    def __len__(self):
        """Returns the number of clusters in this cover."""
        return len(self._clusters)

    def __str__(self):
        """Returns a string representation of the cover."""
        return self.summary(verbosity=1, width=78)

    @property
    def membership(self):
        """Returns the membership vector of this cover.

        The membership vector of a cover covering I{n} elements is a list of
        length I{n}, where element I{i} contains the cluster indices of the
        I{i}th item.
        """
        result = [[] for _ in xrange(self._n)]
        for idx, cluster in enumerate(self):
            for item in cluster:
                result[item].append(idx)
        return result

    @property
    def n(self):
        """Returns the number of elements in the set covered by this cover."""
        return self._n

    def size(self, idx):
        """Returns the size of a given cluster.

        @param idx: the cluster in which we are interested.
        """
        return len(self[idx])

    def sizes(self, *args):
        """Returns the size of given clusters.
        
        The indices are given as positional arguments. If there are no
        positional arguments, the function will return the sizes of all clusters.
        """
        if args:
            return [len(self._clusters[idx]) for idx in args]
        return [len(cluster) for cluster in self]

    def size_histogram(self, bin_width = 1):
        """Returns the histogram of cluster sizes.

        @param bin_width: the bin width of the histogram
        @return: a L{Histogram} object
        """
        return Histogram(bin_width, self.sizes())

    def summary(self, verbosity=0, width=None):
        """Returns the summary of the cover.
        
        The summary includes the number of items and clusters, and also the
        list of members for each of the clusters if the verbosity is nonzero.
        
        @param verbosity: determines whether the cluster members should be
          printed. Zero verbosity prints the number of items and clusters only.
        @return: the summary of the cover as a string.
        """
        out = StringIO()
        print >>out, "Cover with %d clusters" % len(self)

        if verbosity < 1:
            return out.getvalue().strip()

        ndigits = len(str(len(self)))
        wrapper = _get_wrapper_for_width(width,
                subsequent_indent = " " * (ndigits+3))

        for idx, cluster in enumerate(self._formatted_cluster_iterator()):
            wrapper.initial_indent = "[%*d] " % (ndigits, idx)
            print >>out, "\n".join(wrapper.wrap(cluster))

        return out.getvalue().strip()

    def _formatted_cluster_iterator(self):
        """Iterates over the clusters and formats them into a string to be
        presented in the summary."""
        for cluster in self:
            yield ", ".join(str(member) for member in cluster)


class VertexCover(Cover):
    """The cover of the vertex set of a graph.

    This class extends L{Cover} by linking it to a specific L{Graph} object.
    It also provides some handy methods like getting the subgraph corresponding
    to a cluster and such.

    @note: since this class is linked to a L{Graph}, destroying the graph by the
      C{del} operator does not free the memory occupied by the graph if there
      exists a L{VertexCover} that references the L{Graph}.

    @undocumented: _formatted_cluster_iterator
    """

    def __init__(self, graph, clusters = None):
        """Creates a cover object for a given graph.

        @param graph: the graph that will be associated to the cover
        @param clusters: the list of clusters. If C{None}, it is assumed
          that there is only a single cluster that covers the whole graph.
        """
        if clusters is None:
            clusters = [range(graph.vcount())]

        Cover.__init__(self, clusters, n = graph.vcount())
        if self._n > graph.vcount():
            raise ValueError("cluster list contains vertex ID larger than the "
                             "number of vertices in the graph")

        self._graph = graph

    def crossing(self):
        """Returns a boolean vector where element M{i} is C{True} iff edge
        M{i} lies between clusters, C{False} otherwise."""
        membership = [frozenset(cluster) for cluster in self.membership]
        return [membership[v1].isdisjoint(membership[v2]) \
                for v1, v2 in self.graph.get_edgelist()]

    @property
    def graph(self):
        """Returns the graph belonging to this object"""
        return self._graph

    def subgraph(self, idx):
        """Get the subgraph belonging to a given cluster.

        @param idx: the cluster index
        @return: a copy of the subgraph
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the cover was constructed.
        """
        return self._graph.subgraph(self[idx])

    def subgraphs(self):
        """Gets all the subgraphs belonging to each of the clusters.

        @return: a list containing copies of the subgraphs
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the cover was constructed.
        """
        return [self._graph.subgraph(cl) for cl in self]

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Plots the cover to the given Cairo context in the given
        bounding box.

        This is done by calling L{Graph.__plot__()} with the same arguments, but
        drawing nice colored blobs around the vertex groups.

        This method understands all the positional and keyword arguments that
        are understood by L{Graph.__plot__()}, only the differences will be
        highlighted here:

          - C{mark_groups}: whether to highlight the vertex clusters by
            colored polygons. Besides the values accepted by L{Graph.__plot__}
            (i.e., a dict mapping colors to vertex indices, a list containing
            lists of vertex indices, or C{False}), the following are also
            accepted:

              - C{True}: all the clusters will be highlighted, the colors matching
                the corresponding color indices from the current palette
                (see the C{palette} keyword argument of L{Graph.__plot__}.

              - A dict mapping cluster indices or tuples of vertex indices to
                color names.  The given clusters or vertex groups will be
                highlighted by the given colors.

              - A list of cluster indices. This is equivalent to passing a
                dict mapping numeric color indices from the current palette
                to cluster indices; therefore, the cluster referred to by element
                I{i} of the list will be highlighted by color I{i} from the
                palette.

            The value of the C{plotting.mark_groups} configuration key is also
            taken into account here; if that configuration key is C{True} and
            C{mark_groups} is not given explicitly, it will automatically be set
            to C{True}.

            In place of lists of vertex indices, you may also use L{VertexSeq}
            instances.

            In place of color names, you may also use color indices into the
            current palette. C{None} as a color name will mean that the
            corresponding group is ignored.

          - C{palette}: the palette used to resolve numeric color indices to RGBA
            values. By default, this is an instance of L{ClusterColoringPalette}.

        @see: L{Graph.__plot__()} for more supported keyword arguments.
        """
        if "edge_color" not in kwds and "color" not in self.graph.edge_attributes():
            # Set up a default edge coloring based on internal vs external edges
            colors = ["grey20", "grey80"]
            kwds["edge_color"] = [colors[is_crossing]
                                  for is_crossing in self.crossing()]

        if "palette" in kwds:
            palette = kwds["palette"]
        else:
            palette = ClusterColoringPalette(len(self))

        if "mark_groups" not in kwds:
            if Configuration.instance()["plotting.mark_groups"]:
                kwds["mark_groups"] = enumerate(self) 
        else:
            kwds["mark_groups"] = _handle_mark_groups_arg_for_clustering(
                    kwds["mark_groups"], self)

        return self._graph.__plot__(context, bbox, palette, *args, **kwds)

    def _formatted_cluster_iterator(self):
        """Iterates over the clusters and formats them into a string to be
        presented in the summary."""
        if self._graph.is_named():
            names = self._graph.vs["name"]
            for cluster in self:
                yield ", ".join(str(names[member]) for member in cluster)
        else:
            for cluster in self:
                yield ", ".join(str(member) for member in cluster)


class CohesiveBlocks(VertexCover):
    """The cohesive block structure of a graph.

    Instances of this type are created by L{Graph.cohesive_blocks()}. See
    the documentation of L{Graph.cohesive_blocks()} for an explanation of
    what cohesive blocks are.

    This class provides a few more methods that make handling of cohesive
    block structures easier.
    """

    def __init__(self, graph, blocks = None, cohesion = None, parent = None):
        """Constructs a new cohesive block structure for the given graph.

        If any of I{blocks}, I{cohesion} or I{parent} is C{None}, all the
        arguments will be ignored and L{Graph.cohesive_blocks()} will be
        called to calculate the cohesive blocks. Otherwise, these three
        variables should describe the *result* of a cohesive block structure
        calculation. Chances are that you never have to construct L{CohesiveBlocks}
        instances directly, just use L{Graph.cohesive_blocks()}.

        @param graph: the graph itself
        @param blocks: a list containing the blocks; each block is described
          as a list containing vertex IDs.
        @param cohesion: the cohesion of each block. The length of this list
          must be equal to the length of I{blocks}.
        @param parent: the parent block of each block. Negative values or
          C{None} mean that there is no parent block for that block. There
          should be only one parent block, which covers the entire graph.
        @see: Graph.cohesive_blocks()
        """
        if blocks is None or cohesion is None or parent is None:
            blocks, cohesion, parent = graph.cohesive_blocks()

        VertexCover.__init__(self, graph, blocks)

        self._cohesion = cohesion
        self._parent = parent
        for idx, p in enumerate(self._parent):
            if p < 0:
                self._parent[idx] = None

    def cohesion(self, idx):
        """Returns the cohesion of the group with the given index."""
        return self._cohesion[idx]

    def cohesions(self):
        """Returns the list of cohesion values for each group."""
        return self._cohesion[:]

    def hierarchy(self):
        """Returns a new graph that describes the hierarchical relationships
        between the groups.

        The new graph will be a directed tree; an edge will point from
        vertex M{i} to vertex M{j} if group M{i} is a superset of group M{j}.
        In other words, the edges point downwards.
        """
        from igraph import Graph
        edges = [pair for pair in izip(self._parent, xrange(len(self)))
                 if pair[0] is not None]
        return Graph(edges, directed=True)

    def max_cohesion(self, idx):
        """Finds the maximum cohesion score among all the groups that contain
        the given vertex."""
        result = 0
        for cohesion, cluster in izip(self._cohesion, self._clusters):
            if idx in cluster:
                result = max(result, cohesion)
        return result

    def max_cohesions(self):
        """For each vertex in the graph, returns the maximum cohesion score
        among all the groups that contain the vertex."""
        result = [0] * self._graph.vcount()
        for cohesion, cluster in izip(self._cohesion, self._clusters):
            for idx in cluster:
                result[idx] = max(result[idx], cohesion)
        return result

    def parent(self, idx):
        """Returns the parent group index of the group with the given index
        or C{None} if the given group is the root."""
        return self._parent[idx]

    def parents(self):
        """Returns the list of parent group indices for each group or C{None}
        if the given group is the root."""
        return self._parent[:]

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Plots the cohesive block structure to the given Cairo context in
        the given bounding box.

        Since a L{CohesiveBlocks} instance is also a L{VertexCover}, keyword
        arguments accepted by L{VertexCover.__plot__()} are also accepted here.
        The only difference is that the vertices are colored according to their
        maximal cohesions by default, and groups are marked by colored blobs
        except the last group which encapsulates the whole graph.

        See the documentation of L{VertexCover.__plot__()} for more details.
        """
        prepare_groups = False
        if "mark_groups" not in kwds:
            if Configuration.instance()["plotting.mark_groups"]:
                prepare_groups = True
        elif kwds["mark_groups"] == True:
            prepare_groups = True

        if prepare_groups:
            colors = [pair for pair in enumerate(self.cohesions())
                if pair[1] > 1]
            kwds["mark_groups"] = colors

        if "vertex_color" not in kwds:
            kwds["vertex_color"] = self.max_cohesions()

        return VertexCover.__plot__(self, context, bbox, palette, *args, **kwds)

def _handle_mark_groups_arg_for_clustering(mark_groups, clustering):
    """Handles the mark_groups=... keyword argument in plotting methods of
    clusterings.

    This is an internal method, you shouldn't need to mess around with it.
    Its purpose is to handle the extended semantics of the mark_groups=...
    keyword argument in the C{__plot__} method of L{VertexClustering} and
    L{VertexCover} instances, namely the feature that numeric IDs are resolved
    to clusters automatically.
    """
    # Handle the case of mark_groups = True, mark_groups containing a list or
    # tuple of cluster IDs, and and mark_groups yielding (cluster ID, color)
    # pairs
    if mark_groups is True:
        group_iter = ((group, color) for color, group in enumerate(clustering))
    elif isinstance(mark_groups, dict):
        group_iter = mark_groups.iteritems()
    elif hasattr(mark_groups, "__getitem__") and hasattr(mark_groups, "__len__"):
        # Lists, tuples
        try:
            first = mark_groups[0]
        except:
            # Hmm. Maybe not a list or tuple?
            first = None
        if first is not None:
            # Okay. Is the first element of the list a single number?
            if isinstance(first, (int, long)):
                # Yes. Seems like we have a list of cluster indices.
                # Assign color indices automatically.
                group_iter = ((group, color)
                        for color, group in enumerate(mark_groups))
            else:
                # No. Seems like we have good ol' group-color pairs.
                group_iter = mark_groups
        else:
            group_iter = mark_groups
    elif hasattr(mark_groups, "__iter__"):
        # Iterators etc
        group_iter = mark_groups
    else:
        group_iter = {}.iteritems()

    def cluster_index_resolver():
        for group, color in group_iter:
            if isinstance(group, (int, long)):
                group = clustering[group]
            yield group, color

    return cluster_index_resolver()

##############################################################

def _prepare_community_comparison(comm1, comm2, remove_none=False):
    """Auxiliary method that takes two community structures either as
    membership lists or instances of L{Clustering}, and returns a
    tuple whose two elements are membership lists.

    This is used by L{compare_communities} and L{split_join_distance}.

    @param comm1: the first community structure as a membership list or
      as a L{Clustering} object.
    @param comm2: the second community structure as a membership list or
      as a L{Clustering} object.
    @param remove_none: whether to remove C{None} entries from the membership
      lists. If C{remove_none} is C{False}, a C{None} entry in either C{comm1}
      or C{comm2} will result in an exception. If C{remove_none} is C{True},
      C{None} values are filtered away and only the remaining lists are
      compared.
    """
    def _ensure_list(obj):
        if isinstance(obj, Clustering):
            return obj.membership
        return list(obj)

    vec1, vec2 = _ensure_list(comm1), _ensure_list(comm2)
    if len(vec1) != len(vec2):
        raise ValueError("the two membership vectors must be equal in length")

    if remove_none and (None in vec1 or None in vec2):
        idxs_to_remove = [i for i in xrange(len(vec1)) \
                if vec1[i] is None or vec2[i] is None]
        idxs_to_remove.reverse()
        n = len(vec1)
        for i in idxs_to_remove:
            n -= 1
            vec1[i], vec1[n] = vec1[n], vec1[i]
            vec2[i], vec2[n] = vec2[n], vec2[i]
        del vec1[n:]
        del vec2[n:]

    return vec1, vec2


def compare_communities(comm1, comm2, method="vi", remove_none=False):
    """Compares two community structures using various distance measures.

    @param comm1: the first community structure as a membership list or
      as a L{Clustering} object.
    @param comm2: the second community structure as a membership list or
      as a L{Clustering} object.
    @param method: the measure to use. C{"vi"} or C{"meila"} means the
      variation of information metric of Meila (2003), C{"nmi"} or C{"danon"}
      means the normalized mutual information as defined by Danon et al (2005),
      C{"split-join"} means the split-join distance of van Dongen (2000),
      C{"rand"} means the Rand index of Rand (1971), C{"adjusted_rand"}
      means the adjusted Rand index of Hubert and Arabie (1985).
    @param remove_none: whether to remove C{None} entries from the membership
      lists. This is handy if your L{Clustering} object was constructed using
      L{VertexClustering.FromAttribute} using an attribute which was not defined
      for all the vertices. If C{remove_none} is C{False}, a C{None} entry in
      either C{comm1} or C{comm2} will result in an exception. If C{remove_none}
      is C{True}, C{None} values are filtered away and only the remaining lists
      are compared.

    @return: the calculated measure.
    @newfield ref: Reference
    @ref: Meila M: Comparing clusterings by the variation of information.
          In: Scholkopf B, Warmuth MK (eds). Learning Theory and Kernel
          Machines: 16th Annual Conference on Computational Learning Theory
          and 7th Kernel Workship, COLT/Kernel 2003, Washington, DC, USA.
          Lecture Notes in Computer Science, vol. 2777, Springer, 2003.
          ISBN: 978-3-540-40720-1.
    @ref: Danon L, Diaz-Guilera A, Duch J, Arenas A: Comparing community
          structure identification. J Stat Mech P09008, 2005.
    @ref: van Dongen D: Performance criteria for graph clustering and Markov
          cluster experiments. Technical Report INS-R0012, National Research
          Institute for Mathematics and Computer Science in the Netherlands,
          Amsterdam, May 2000.
    @ref: Rand WM: Objective criteria for the evaluation of clustering
          methods. J Am Stat Assoc 66(336):846-850, 1971.
    @ref: Hubert L and Arabie P: Comparing partitions. Journal of
          Classification 2:193-218, 1985.
    """
    import igraph._igraph
    vec1, vec2 = _prepare_community_comparison(comm1, comm2, remove_none)
    return igraph._igraph._compare_communities(vec1, vec2, method)


def split_join_distance(comm1, comm2, remove_none=False):
    """Calculates the split-join distance between two community structures.

    The split-join distance is a distance measure defined on the space of
    partitions of a given set. It is the sum of the projection distance of
    one partition from the other and vice versa, where the projection
    number of A from B is if calculated as follows:

      1. For each set in A, find the set in B with which it has the
         maximal overlap, and take note of the size of the overlap.

      2. Take the sum of the maximal overlap sizes for each set in A.

      3. Subtract the sum from M{n}, the number of elements in the
         partition.

    Note that the projection distance is asymmetric, that's why it has to be
    calculated in both directions and then added together.  This function
    returns the projection distance of C{comm1} from C{comm2} and the
    projection distance of C{comm2} from C{comm1}, and returns them in a pair.
    The actual split-join distance is the sum of the two distances. The reason
    why it is presented this way is that one of the elements being zero then
    implies that one of the partitions is a subpartition of the other (and if
    it is close to zero, then one of the partitions is close to being a
    subpartition of the other).

    @param comm1: the first community structure as a membership list or
      as a L{Clustering} object.
    @param comm2: the second community structure as a membership list or
      as a L{Clustering} object.
    @param remove_none: whether to remove C{None} entries from the membership
      lists. This is handy if your L{Clustering} object was constructed using
      L{VertexClustering.FromAttribute} using an attribute which was not defined
      for all the vertices. If C{remove_none} is C{False}, a C{None} entry in
      either C{comm1} or C{comm2} will result in an exception. If C{remove_none}
      is C{True}, C{None} values are filtered away and only the remaining lists
      are compared.

    @return: the projection distance of C{comm1} from C{comm2} and vice versa
      in a tuple. The split-join distance is the sum of the two.
    @newfield ref: Reference
    @ref: van Dongen D: Performance criteria for graph clustering and Markov
          cluster experiments. Technical Report INS-R0012, National Research
          Institute for Mathematics and Computer Science in the Netherlands,
          Amsterdam, May 2000.

    @see: L{compare_communities()} with C{method = "split-join"} if you are
      not interested in the individual projection distances but only the
      sum of them.
    """
    import igraph._igraph
    vec1, vec2 = _prepare_community_comparison(comm1, comm2, remove_none)
    return igraph._igraph._split_join_distance(vec1, vec2)



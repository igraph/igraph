"""Objects related to graph clustering"""

__license__ = """
Copyright (C) 2006-2007  Gabor Csardi <csardi@rmki.kfki.hu>,
Tamas Nepusz <ntamas@rmki.kfki.hu>

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
Foundation, Inc.,  51 Franklin Street, Fifth Floor, Boston, MA 
02110-1301 USA
"""

from copy import deepcopy
from StringIO import StringIO

from igraph import community_to_membership
from igraph.configuration import Configuration
from igraph.datatypes import UniqueIdGenerator
from igraph.drawing.colors import ClusterColoringPalette
from igraph.statistics import Histogram

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

    @property
    def membership(self):
        """Returns the membership vector."""
        return self._membership[:]

    def size(self, idx):
        """Returns the size of a given cluster.

        @param idx: the cluster in which we are interested.
        """
        return len(self[idx])

    def sizes(self, *args):
        """Returns the size of given clusters.

        @keyword idxs: the cluster indices in which we are interested. If C{None},
          defaults to all clusters
        """
        idxs = args
        if len(idxs) == 0:
            idxs = None

        counts = [0] * len(self)
        for x in self._membership:
            counts[x] += 1

        if idxs is None:
            return counts

        return [counts[idx] for idx in idxs]
    
    def size_histogram(self, bin_width = 1):
        """Returns the histogram of cluster sizes.

        @param bin_width: the bin width of the histogram
        @return: a L{Histogram} object
        """
        return Histogram(bin_width, self.sizes())


class OverlappingClustering(Clustering):
    """Extension of L{Clustering} that allows for overlapping clusters.

    With overlapping clusters, a single vertex can be the member of multiple
    clusters (or even none of them). Therefore, each item of the membership
    vector is a list, set or tuple containing the cluster indices for all
    vertices.
    
    Members of an individual cluster can be accessed by the C{[]} operator:
    
      >>> cl = OverlappingClustering([(0,), (0,), (0,1), (1,), (1,), ()])
      >>> cl[0]
      [0, 1, 2]
      >>> cl[1]
      [2, 3, 4]
    
    The membership vector can be accessed by the C{membership} property:

      >>> cl.membership
      [frozenset([0]), frozenset([0]), frozenset([0, 1]), frozenset([1]), frozenset([1]), frozenset([])]

    The number of clusters can be retrieved by the C{len} function:

      >>> len(cl)
      2
    """
    def __init__(self, membership, params = None):
        """Constructor.

        @param membership: the membership list -- that is, the cluster
          index in which each element of the set belongs to.
        @param params: additional parameters to be stored in this
          object's dictionary."""
        Clustering.__init__(self, [], params)

        self._membership = [frozenset(clusters) for clusters in membership]
        try:
            self._len = max(max(m) for m in self._membership if m) + 1
        except ValueError:
            self._len = 0
    
    def __getitem__(self, idx):
        """Returns the members of the specified cluster.

        @param idx: the index of the cluster
        @return: the members of the specified cluster as a list
        @raise IndexError: if the index is out of bounds"""
        if idx < 0 or idx >= self._len:
            raise IndexError("cluster index out of range")
        return [i for i, cl in enumerate(self._membership) if idx in cl]
   
    def sizes(self, *args):
        """Returns the size of given clusters.

        @keyword idxs: the cluster indices in which we are interested. If C{None},
          defaults to all clusters
        """
        idxs = args
        if len(idxs) == 0:
            idxs = None

        counts = [0] * len(self)
        for members in self._membership:
            for member in members:
                counts[member] += 1

        if idxs is None:
            return counts

        return [counts[idx] for idx in idxs]
    

class VertexClustering(Clustering):
    """The clustering of the vertex set of a graph.

    This class extends L{Clustering} by linking it to a specific L{Graph} object
    and by optionally storing the modularity score of the clustering.

    @note: since this class is linked to a L{Graph}, destroying the graph by the
      C{del} operator does not free the memory occupied by the graph if there
      exists a L{VertexClustering} that references the L{Graph}.
    """

    def __init__(self, graph, membership = None, modularity = None, \
                 params = None):
        """Creates a clustering object for a given graph.

        @param graph: the graph that will be associated to the clustering
        @param membership: the membership list. The length of the list must
          be equal to the number of vertices in the graph. If C{None}, every
          vertex is assumed to belong to the same cluster.
        @param modularity: the modularity score of the clustering. If C{None},
          it will be calculated when needed.
        @param params: additional parameters to be stored in this object.
        """
        self._graph = graph

        if membership is None:
            Clustering.__init__(self, [0]*graph.vcount(), params)
        else:
            if len(membership) != graph.vcount():
                raise ValueError("membership list has invalid length")
            Clustering.__init__(self, membership, params)

        self._modularity = modularity

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
        self._modularity = self._graph.modularity(self._membership)
        return self._modularity


    def subgraph(self, idx):
        """Get the subgraph belonging to a given cluster.

        @param idx: the cluster index
        @return: a copy of the subgraph
        @precondition: the vertex set of the graph hasn't been modified since
          the moment the clustering was constructed.
        """
        return self._graph.subgraph(self[idx])


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
        coloring the graph vertices according to the current clustering.

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

              - A dict mapping color names to cluster indices. The given clusters
                will be highlighted by the given colors.

              - A list of cluster indices. This is equivalent to passing a
                dict mapping numeric color indices from the current palette
                to cluster indices; therefore, the cluster referred to by element
                I{i} of the list will be highlighted by color I{i} from the
                palette.

            The value of the C{plotting.mark_groups} configuration key is also
            taken into account here; if that configuration key is C{True} and
            C{mark_groups} is not given explicitly, it will automatically be set
            to C{True}.

          - C{palette}: the palette used to resolve numeric color indices to RGBA
            values. By default, this is an instance of L{ClusterColoringPalette}.

          - C{vertex_color}: this keyword argument is not allowed as it would override
            the coloring.

        @see: L{Graph.__plot__()} for more supported keyword arguments.
        """
        if "vertex_color" in kwds:
            raise ValueError("you are not allowed to define vertex colors "+
                             "when plotting a clustering")

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
            mark_groups = kwds["mark_groups"]
            del kwds["mark_groups"]

            # Handle the case of mark_groups = True and mark_groups yielding
            # cluster IDs
            if mark_groups is True:
                group_iter = enumerate(self)
            elif isinstance(mark_groups, dict):
                group_iter = mark_groups.iteritems()
            elif hasattr(mark_groups, "__iter__"):
                if hasattr(mark_groups, "next"):
                    # Already an iterator, let's hope it works
                    group_iter = mark_groups
                else:
                    # Lists, tuples etc
                    group_iter = enumerate(mark_groups)
            else:
                group_iter = {}.iteritems()

            def cluster_index_resolver():
                for color_id, group in group_iter:
                    if isinstance(group, (int, long)):
                        group = self[group]
                    yield color_id, group

            kwds["mark_groups"] = cluster_index_resolver()

        kwds["vertex_color"] = self.membership
        return self._graph.__plot__(context, bbox, palette, *args, **kwds)


class OverlappingVertexClustering(OverlappingClustering, VertexClustering):
    """Overlapping clustering of the vertex set of a graph.

    This class extends L{OverlappingClustering} by linking it to a specific
    L{Graph} object and by optionally storing the modularity score of the
    clustering.

    Modularity in the case of overlapping communities is defined similarly
    to the nonoverlapping case, but the statement that ``vertex M{i} and
    M{j} is in the same community'' (expressed by the Kronecker-delta at
    the end of the formula in the original paper) is replaced by the
    statement that ``vertex M{i} and M{j} are both contained by at least
    one of the communities''.

    @note: since this class is linked to a L{Graph}, destroying the graph by the
      C{del} operator does not free the memory occupied by the graph if there
      exists an L{OverlappingVertexClustering} that references the L{Graph}.
    """

    def __init__(self, graph, membership = None, modularity = None, \
            params = None):
        """Creates an overlapping clustering object for a given graph.

        @param graph: the graph that will be associated to the clustering
        @param membership: the membership list. The length of the list must
          be equal to the number of vertices in the graph. If C{None}, every
          vertex is assumed to belong to the same cluster.
        @param modularity: the modularity score of the clustering. If C{None},
          it will be calculated.
        @param params: additional parameters to be stored in this object.
        """
        self._graph = graph

        if membership is None:
            membership = [set(0)] * graph.vcount()
        if len(membership) != graph.vcount():
            raise ValueError("membership list is too short")
        OverlappingClustering.__init__(self, membership, params)

        if modularity is None:
            self._modularity = self.recalculate_modularity()
        else:
            self._modularity = modularity

    def recalculate_modularity(self):
        """Recalculates the stored modularity value.

        This method must be called before querying the modularity score of the
        clustering through the class member C{modularity} or C{q} if the
        graph has been modified (edges have been added or removed) since the
        creation of the L{OverlappingVertexClustering} object.
        
        @return: the new modularity score
        @todo: this is pretty slow now, it should eventually be moved to
          the C layer.
        """
        degrees = self._graph.degree()
        edge_set = set(self._graph.get_edgelist())
        if not self._graph.is_directed():
            mirrored_edge_set = set((v2, v1) for v1, v2 in edge_set)
            edge_set = edge_set.union(mirrored_edge_set)
        ecount = float(len(edge_set))
        result = 0.0
        for source, cl1 in enumerate(self._membership):
            source_degree = degrees[source]
            for target, cl2 in enumerate(self._membership):
                if len(cl1.intersection(cl2))>0:
                    if (source, target) in edge_set:
                        result += 1.0
                    result -= source_degree*degrees[target] / ecount

        self._modularity = result / ecount
        return self._modularity


class Dendrogram(Clustering):
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
    """

    def __init__(self, merges):
        """Creates a hierarchical clustering.

        @param merges: the merge history either in matrix or tuple format"""
        Clustering.__init__(self, [0]*(len(merges)+1))
        self._merges = [tuple(pair) for pair in merges]
        self._nmerges = len(self._merges)
        self._nitems = max(self._merges[-1])-self._nmerges+2
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
        return "Dendrogram, %d elements, %d merges" % \
                (self._nitems, self._nmerges)

    def summary(self):
        """Draws the dendrogram of the hierarchical clustering in a string"""
        from array import array

        out = StringIO()
        print >>out, str(self)
        if self._nitems == 0:
            return out.getvalue()
            
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
                char_array[pos1:(pos2+1)] = array("c", "+%s+" % dashes)

                cidx_incr += 1
            
            max_community_idx += cidx_incr

            print >>out, char_array.tostring()


        return out.getvalue()

    def _item_box_size(self, context, horiz, idx):
        """Calculates the amount of space needed for drawing an
        individual vertex at the bottom of the dendrogram."""
        if self._names[idx] is None:
            _, _, _, height, x_ascent, _ = context.text_extents("")
        else:
            _, _, _, height, x_ascent, _ = context.text_extents(str(self._names[idx]))
        if horiz:
            return x_ascent, height
        return height, x_ascent

    # pylint: disable-msg=R0913
    def _plot_item(self, context, horiz, idx, x, y):
        """Plots a dendrogram item to the given Cairo context

        @param context: the Cairo context we are plotting on
        @param horiz: whether the dendrogram is horizontally oriented
        @param idx: the index of the item
        @param x: the X position of the item
        @param y: the Y position of the item
        """
        if not self._names[idx]:
            return

        height, _ = self._item_box_size(context, horiz, idx)
        if horiz:
            context.move_to(x, y+height)
            context.show_text(str(self._names[idx]))
        else:
            context.save()
            context.translate(x, y)
            context.rotate(-1.5707963285)    # pi/2
            context.move_to(0, height)
            context.show_text(self._names[idx])
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

        orientation = kwds.get("orientation", "lr")
        
        orientation_aliases = {
            "left-right": "lr", "right-left": "rl",
            "top-bottom": "tb", "bottom-top": "bt",
            "horizontal": "lr", "horiz": "lr", "h": "lr",
            "vertical": "bt", "vert": "bt", "v": "bt"
        }
        orientation = orientation_aliases.get(orientation, orientation)
        if orientation not in ("lr", "rl", "tb", "bt"):
            raise ValueError("unknown orientation: %s" % orientation)
        horiz = orientation in ("lr", "rl")

        # Calculate space needed for individual items at the
        # bottom of the dendrogram
        item_boxes = [self._item_box_size(context, horiz, idx) \
          for idx in xrange(self._nitems)]

        # Calculate coordinates
        layout = Layout([(0, 0)] * self._nitems, dim=2)
        inorder = self._traverse_inorder()
        if not horiz:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                layout[element] = (x + item_boxes[element][0]/2., 0)
                x += item_boxes[element][0]

            for id1, id2 in self._merges:
                y += 1
                layout.append(((layout[id1][0]+layout[id2][0])/2., y))

            # Mirror or rotate the layout if necessary
            if orientation == "bt":
                layout.mirror(1)
        else:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                layout[element] = (0, y + item_boxes[element][1]/2.)
                y += item_boxes[element][1]

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

        bbox = layout.bounding_box()
        rx, ry = width / max(bbox.width, 1), height / max(bbox.height, 1)
        delta_x -= (bbox.left * rx - bbox.left)
        delta_y -= (bbox.top * ry - bbox.top)
        layout.scale(rx, ry)
        layout.translate(delta_x, delta_y)

        context.set_source_rgb(0., 0., 0.)
        context.set_line_width(1)
        
        # Draw items
        if horiz:
            sgn = -1
            if orientation == "rl":
                sgn = 0
            for idx in xrange(self._nitems):
                x = layout[idx][0] + sgn * item_boxes[idx][0]
                y = layout[idx][1] - item_boxes[idx][1]/2.
                self._plot_item(context, horiz, idx, x, y)
        else:
            sgn = 0
            if orientation == "bt":
                sgn = 1
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

    @property
    def merges(self):
        """Returns the performed merges in matrix format"""
        return deepcopy(self._merges)

class VertexDendrogram(VertexClustering, Dendrogram):
    """The dendrogram resulting from the hierarchical clustering of the
    vertex set of a graph."""

    def __init__(self, graph, merges, membership = None, modularity = None, \
            params = None):
        """Creates a dendrogram object for a given graph.

        @param graph: the graph that will be associated to the clustering
        @param merges: the merges performed given in matrix form.
        @param membership: the membership list. The length of the list must
          be equal to the number of vertices in the graph. If C{None}, the
          dendrogram will be cut at the level where the modularity is maximized
          and the membership list will represent this state.
        @param modularity: the modularity score of the clustering on each
          level of the dendrogram starting from the fully decomposed state.
          If C{None}, it will be calculated.
        @param params: additional parameters to be stored in this object.
        """
        if modularity is None:
            # TODO: this is a fairly simple way to calculate the modularity
            membs = range(graph.vcount())
            modularity = []
            n = graph.vcount()
            for step in xrange(min(n-1, len(merges))):
                membs = community_to_membership(merges, n, step)
                modularity.append(graph.modularity(membs))

        if membership is None:
            maxmod = max(modularity)
            maxidx = modularity.index(maxmod)
            membership = community_to_membership(merges, graph.vcount(), maxidx)
            
            recoding, n = {}, 0
            for idx, m in enumerate(membership):
                try:
                    membership[idx] = recoding[m]
                except KeyError:
                    recoding[m], membership[idx] = n, n
                    n += 1

        else:
            maxmod = None

        Dendrogram.__init__(self, merges)
        VertexClustering.__init__(self, graph, membership, maxmod, params)

    def as_clustering(self, n=None):
        """Cuts the dendrogram at the given level and returns a corresponding
        L{VertexClustering} object.

        @param n: the desired number of clusters. Merges are replayed from the
          beginning until the membership vector has exactly M{n} distinct elements
          or until there are no more recorded merges, whichever happens first.
          If C{None}, the current membership vector will be used.
        @return: a new L{VertexClustering} object.
        """
        if n is not None:
            num_elts = self._graph.vcount()
            membership = community_to_membership(self._merges, num_elts, \
                                                 num_elts - n)
            idgen = UniqueIdGenerator()
            membership = [idgen[m] for m in membership]
        else:
            membership = self.membership
        return VertexClustering(self._graph, membership)


    def cut(self, n):
        """Cuts the dendrogram at a given level.

        @param n: the desired number of clusters. Merges are replayed from the
          beginning until the membership vector has exactly M{n} distinct elements
          or until there are no more recorded merges, whichever happens first.
        @return: the membership vector
        """
        num_elts = self._graph.vcount()
        membership = community_to_membership(self._merges, num_elts, num_elts-n)
        idgen = UniqueIdGenerator()
        self._membership = [idgen[m] for m in membership]
        self._len = max(self._membership) + 1
        return self._membership[:]

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Draws the vertex dendrogram on the given Cairo context

        See L{Dendrogram.__plot__} for the list of supported keyword
        arguments."""
        from igraph.drawing.metamagic import AttributeCollectorBase

        class VisualVertexBuilder(AttributeCollectorBase):
            label = None

        builder = VisualVertexBuilder(self._graph.vs, kwds)
        self._names = [vertex.label for vertex in builder]
        result = Dendrogram.__plot__(self, context, bbox, palette, \
                *args, **kwds)
        del self._names

        return result


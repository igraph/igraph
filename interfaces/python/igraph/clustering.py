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

from statistics import Histogram
from copy import copy, deepcopy
from StringIO import StringIO
from igraph import community_to_membership
try:
    set, frozenset
except NameError:
    import sets
    set, frozenset = sets.Set, sets.ImmutableSet

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
    """

    def __init__(self, membership, params = {}):
        """Constructor.

        @param membership: the membership list -- that is, the cluster
          index in which each element of the set belongs to.
        @param params: additional parameters to be stored in this
          object's dictionary."""
        self._membership = list(membership)
        if len(self._membership)>0:
            self._len = max(self._membership)+1
        else:
            self._len = 0
        self.__dict__.update(params)
    
    def __getitem__(self, idx):
        """Returns the members of the specified cluster.

        @param idx: the index of the cluster
        @return: the members of the specified cluster as a list
        @raise IndexError: if the index is out of bounds"""
        if idx<0 or idx>=self._len:
            raise IndexError, "cluster index out of range"
        return [i for i,e in enumerate(self._membership) if e==idx]

    def __len__(self):
        """Returns the number of clusters.

        @return: the number of clusters
        """
        return self._len

    def _get_membership(self): return deepcopy(self._membership)
    membership = property(_get_membership, doc = "The membership vector (read only)")

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
        if len(idxs) == 0: idxs = None
        counts = [0] * len(self)
        m = min(self._membership)
        for x in self._membership: counts[x+m] += 1
        if idxs is None: return counts
        result = []
        for idx in idxs: result.append(counts[idx-m])
        return result
    
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
    
      >>> cl = Clustering([(0,), (0,), (0,1), (1,), (1,), ()])
      >>> cl[0]
      [0, 1, 2]
      >>> cl[1]
      [2, 3, 4]
    
    The membership vector can be accessed by the C{membership} property:

      >>> cl.membership
      [frozenset([0]), frozenset([0]), frozenset([0,1]), frozenset([1]), frozenset([1]), frozenset([])]

    The number of clusters can be retrieved by the C{len} function:

      >>> len(cl)
      2
    """
    def __init__(self, membership, params = {}):
        """Constructor.

        @param membership: the membership list -- that is, the cluster
          index in which each element of the set belongs to.
        @param params: additional parameters to be stored in this
          object's dictionary."""
        self._membership = map(frozenset, list(membership))
        self._len = -1
        for m in self._membership:
            if len(m) == 0: continue
            self._len = max(self._len, max(m))
        self._len += 1
        self.__dict__.update(params)
    
    def __getitem__(self, idx):
        """Returns the members of the specified cluster.

        @param idx: the index of the cluster
        @return: the members of the specified cluster as a list
        @raise IndexError: if the index is out of bounds"""
        if idx<0 or idx>=self._len:
            raise IndexError, "cluster index out of range"
        return [i for i,cl in enumerate(self._membership) if idx in cl]
   
    def sizes(self, *args):
        """Returns the size of given clusters.

        @keyword idxs: the cluster indices in which we are interested. If C{None},
          defaults to all clusters
        """
        idxs = args
        if len(idxs) == 0: idxs = None
        counts = [0] * len(self)
        for members in self._membership:
            for member in members: counts[member] += 1
        if idxs is None: return counts
        result = []
        for idx in idxs: result.append(counts[idx])
        return result
    

class VertexClustering(Clustering):
    """The clustering of the vertex set of a graph.

    This class extends L{Clustering} by linking it to a specific L{Graph} object
    and by optionally storing the modularity score of the clustering.

    @note: since this class is linked to a L{Graph}, destroying the graph by the
      C{del} operator does not free the memory occupied by the graph if there
      exists a L{VertexClustering} that references the L{Graph}.
    """

    def __init__(self, graph, membership = None, modularity = None, params = {}):
        """Creates a clustering object for a given graph.

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
            Clustering.__init__(self, [0]*graph.vcount(), params)
        else:
            if len(membership) != graph.vcount():
                raise ValueError, "membership list has invalid length"
            Clustering.__init__(self, membership, params)

        if modularity is None:
            self._q = graph.modularity(membership)
        else:
            self._q = modularity

    def _get_modularity(self): return self._q
    modularity = property(_get_modularity, doc = "The modularity score")
    q = modularity

    def _get_graph(self): return self._graph
    graph = property(_get_graph, doc = "The graph belonging to this object")

    def recalculate_modularity(self):
        """Recalculates the stored modularity value.

        This method must be called before querying the modularity score of the
        clustering through the class member C{modularity} or C{q} if the
        graph has been modified (edges have been added or removed) since the
        creation of the L{VertexClustering} object.
        
        @return: the new modularity score
        """
        self._q = self._graph.modularity(self._membership)
        return self._q


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

    def __init__(self, graph, membership = None, modularity = None, params = {}):
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
            OverlappingClustering.__init__(self, [set(0)]*graph.vcount(), params)
        else:
            if len(membership) != graph.vcount():
                raise ValueError, "membership list is too short"
            OverlappingClustering.__init__(self, membership, params)

        if modularity is None:
            self._q = self.recalculate_modularity()
        else:
            self._q = modularity

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
        n = self._graph.vcount()
        el = set(self._graph.get_edgelist())
        if not self._graph.is_directed():
            el2 = set([(v2,v1) for v1,v2 in self._graph.get_edgelist()])
            el = el.union(el2)
        ecount = float(len(el))
        result = 0.0
        for v1 in xrange(n):
            for v2 in xrange(n):
                if len(self._membership[v1].intersection(self._membership[v2]))>0:
                    if (v1,v2) in el: result += 1.0
                    result -= degrees[v1]*degrees[v2] / ecount

        self._q = result / ecount
        return self._q


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
        self._n = max(self._merges[-1])-self._nmerges+2
    
    def _convert_matrix_to_tuple_repr(merges, n=None):
        if n is None: n = len(merges)+1
        t = range(n)
        idxs = range(n)
        for rowidx, row in enumerate(merges):
            i, j = row
            try:
                idxi, idxj = idxs[i], idxs[j]
                t[idxi] = (t[idxi], t[idxj])
                t[idxj] = None
            except IndexError:
                raise ValueError, "malformed matrix, subgroup referenced before being created in step %d" % rowidx
            idxs.append(j)
        return [x for x in t if x is not None]

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
            if last < self._n: # This will be a regular node
                result.append(last)
            else:        # This is a merge node, proceed towards left
                stack.append(self._merges[last-self._n])

        return result

    def __str__(self):
        return "Dendrogram, %d elements, %d merges" % (self._n, self._nmerges)

    def summary(self):
        """Draws the dendrogram of the hierarchical clustering in a string"""
        out = StringIO()
        print >>out, str(self)
        if self._n == 0: return out.getvalue()
            
        print >>out

        positions = [None] * self._n
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
        print >>out, (" "*(distance-1)).join(inorder)
        midx = 0
        max_community_idx = self._n
        while midx < self._nmerges:
            s = [" "]*width
            for p in positions:
                if p >= 0: s[p] = "|"
            for i in xrange(level_distance-1): print >>out, "".join(s) # Print the lines
            
            cidx_incr = 0
            while midx < self._nmerges:
                v1, v2 = self._merges[midx]
                if v1 >= max_community_idx or v2 >= max_community_idx: break
                midx += 1
                p1 = positions[v1]
                p2 = positions[v2]
                positions[v1] = -1
                positions[v2] = -1
                positions.append((p1+p2)/2)
                s[p1:(p2+1)] = "+%s+" % ("-" * (p2-p1-1))
                cidx_incr += 1
            
            max_community_idx += cidx_incr

            print >>out, "".join(s)


        return out.getvalue()

    def _item_box_size(self, context, horiz, idx):
        """Calculates the amount of space needed for drawing an individual vertex
        at the bottom of the dendrogram."""
        xb, yb, w, h, xa, ya = context.text_extents(self._names[idx])
        if horiz: return (xa, h)
        return (h, xa)

    def _plot_item(self, context, horiz, idx, x, y):
        xb, yb, w, h, xa, ya = context.text_extents(self._names[idx])
        if horiz:
            context.move_to(x, y+h)
            context.show_text(self._names[idx])
        else:
            context.save()
            context.translate(x, y)
            context.rotate(-1.5707963285)    # pi/2
            context.move_to(0, h)
            context.show_text(self._names[idx])
            context.restore()

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

        if not hasattr(self, "_names"): self._names = map(str, xrange(self._n))

        orientation = kwds.get("orientation", "lr")
        
        orientation_aliases = {
            "lr": "left-right", "rl": "right-left",
            "tb": "top-bottom", "bt": "bottom-top",
            "horizontal": "left-right", "horiz": "left-right", "h": "left-right",
            "vertical": "bottom-top", "vert": "bottom-top", "v": "bottom-top"
        }
        orientation = orientation_aliases.get(orientation, orientation)
        if orientation not in ("left-right", "right-left", "top-bottom", "bottom-top"):
            raise ValueError, "unknown orientation: %s" % orientation
        horiz = orientation in ("left-right", "right-left")

        # Calculate space needed for individual items at the bottom of the dendrogram
        item_boxes = [self._item_box_size(context, horiz, idx) \
          for idx in xrange(self._n)]

        # Calculate coordinates
        w, h = bbox.width, bbox.height
        lo = Layout([(0,0)]*self._n, dim=2)
        inorder = self._traverse_inorder()
        if not horiz:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                lo[element] = (x + item_boxes[element][0]/2., 0)
                x += item_boxes[element][0]

            for c1, c2 in self._merges:
                y += 1
                lo.append(((lo[c1][0]+lo[c2][0])/2., y))

            # Mirror or rotate the layout if necessary
            if orientation == "bottom-top": lo.mirror(1)
        else:
            x, y = 0, 0
            for idx, element in enumerate(inorder):
                lo[element] = (0, y + item_boxes[element][1]/2.)
                y += item_boxes[element][1]

            for c1, c2 in self._merges:
                x += 1
                lo.append((x, (lo[c1][1]+lo[c2][1])/2.))

            # Mirror or rotate the layout if necessary
            if orientation == "right-left": lo.mirror(0)
        
        # Rescale layout to the bounding box
        maxw, maxh = max([e[0] for e in item_boxes]), max([e[1] for e in item_boxes])
        # w, h: width and height of the area containing the dendrogram tree without
        # the items. dx, dy: displacement of the dendrogram tree
        w, h, dx, dy = float(bbox.width), float(bbox.height), 0, 0
        if horiz:
            w -= maxw
            if orientation == "left-right": dx = maxw
        else:
            h -= maxh
            if orientation == "top-bottom": dy = maxh
        sl, st, sr, sb = lo.bounding_box()
        sw, sh = max(sr-sl, 1), max(sb-st, 1)
        rx, ry = w/sw, h/sh
        lo.scale(rx, ry)
        lo.translate(dx-sl*rx+bbox.coords[0], dy-st*ry+bbox.coords[1])

        context.set_source_rgb(0.,0.,0.)
        context.set_line_width(1)
        
        # Draw items
        if horiz:
            sgn = -1
            if orientation == "right-left": sgn = 0
            for idx in xrange(self._n):
                x = lo[idx][0] + sgn * item_boxes[idx][0]
                y = lo[idx][1] - item_boxes[idx][1]/2.
                self._plot_item(context, horiz, idx, x, y)
        else:
            sgn = 0
            if orientation == "bottom-top": sgn = 1
            for idx in xrange(self._n):
                x = lo[idx][0] - item_boxes[idx][0]/2.
                y = lo[idx][1] + sgn * item_boxes[idx][1]
                self._plot_item(context, horiz, idx, x, y)

        # Draw dendrogram lines
        if not horiz:
            for idx, (c1, c2) in enumerate(self._merges):
                x0, y0 = lo[c1]
                x1, y1 = lo[c2]
                x2, y2 = lo[idx + self._n]
                context.move_to(x0, y0)
                context.line_to(x0, y2)
                context.line_to(x1, y2)
                context.line_to(x1, y1)
                context.stroke()
        else:
            for idx, (c1, c2) in enumerate(self._merges):
                x0, y0 = lo[c1]
                x1, y1 = lo[c2]
                x2, y2 = lo[idx + self._n]
                context.move_to(x0, y0)
                context.line_to(x2, y0)
                context.line_to(x2, y1)
                context.line_to(x1, y1)
                context.stroke()

    def _get_merges(self): return copy(self._merges)
    merges = property(_get_merges, doc = "The performed merges in matrix format")

class VertexDendrogram(VertexClustering, Dendrogram):
    """The dendrogram resulting from the hierarchical clustering of the
    vertex set of a graph."""

    def __init__(self, graph, merges, membership = None, modularity = None, params = {}):
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
            ms = range(graph.vcount())
            communities = range(graph.vcount())
            modularity = []
            n = graph.vcount()
            for step in xrange(min(n-1, len(merges))):
                ms = community_to_membership(merges, graph.vcount(), step)
                modularity.append(graph.modularity(ms))

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

    def cut(self, n):
        """Cuts the dendrogram at a given level.

        @param n: the desired number of clusters. Merges are replayed from the
          beginning until the membership vector has exactly M{n} distinct elements
          or until there are no more recorded merges, whichever happens first.
        @return: the membership vector
        """
        num_elts = self._graph.vcount()
        membership = community_to_membership(self._merges, num_elts, num_elts-n)
        
        recoding = {}
        n=0
        for idx, m in enumerate(membership):
            try:
                v = recoding[m]
            except KeyError:
                recoding[m], v = n, n
                n += 1
            membership[idx] = v
 
        self._membership = membership
        self._len = max(membership) + 1
        return copy(membership)

    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Draws the vertex dendrogram on the given Cairo context

        See L{Dendrogram.__plot__} for the list of supported keyword arguments."""
        from igraph.drawing import collect_attributes
        from igraph import config

        if not kwds.has_key("vertex_label") and \
            "label" not in self._graph.vs.attribute_names():
            self._names = map(str, xrange(self._graph.vcount()))
        elif kwds.get("vertex_label", []) is None:
            self._names = map(str, xrange(self._graph.vcount()))
        else:
            self._names = collect_attributes(self._graph.vcount(), \
                "vertex_label", "label", kwds, self._graph.vs, config, None)

        result = Dendrogram.__plot__(self, context, bbox, palette, *args, **kwds)
        del self._names
        return result


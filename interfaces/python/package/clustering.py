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
        @raise: C{IndexError} if the index is out of bounds"""
        if idx<0 or idx>self._len:
            raise IndexError, "cluster index out of range"
        return [i for i,e in enumerate(self._membership) if e==idx]

    def __len__(self):
        """Returns the number of clusters.

        @return: the number of clusters
        """
        return self._len

    def _get_membership(self): return copy(self._membership)
    membership = property(_get_membership, doc = "The membership vector (read only)")

    def size(self, idx):
        """Returns the size of a given cluster.

        @param idx: the cluster in which we are interested.
        """
        if idx<0 or idx>self._len:
            raise IndexError, "cluster index out of range"
        return len([1 for x in self._membership if x == idx])

    def sizes(self, *args):
        """Returns the size of given clusters.

        @param idxs: the cluster indices in which we are interested. If C{None},
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
                raise ValueError, "membership list is too short"
            Clustering.__init__(self, membership, params)

        if modularity is None:
            self._q = graph.modularity(membership)
        else:
            self._q = modularity

    def _get_modularity(self): return self._q
    modularity = property(_get_modularity, doc = "The modularity score")
    q = modularity

    def _get_graph(self): return self._graph
    modularity = property(_get_graph, doc = "The graph belonging to this object")

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
      2 ---+-+        <====>   [[0, 1], [4, 5], [2, 3], [6, 7]]
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
            modularity.append(graph.modularity(ms))
            for c1, c2 in merges:
                try:
                    cidx1 = communities[c1]
                    cidx2 = communities[c2]
                except IndexError:
                    raise ValueError, "invalid merge matrix, referencing nonexisting community"
                if cidx1 == -1 or cidx2 == -1:
                    raise ValueError, "invalid merge matrix, referencing already joined community"
                for idx, m in enumerate(ms):
                    if m == cidx2: ms[idx] = cidx1
                communities.append(communities[c1])
                communities[c1] = -1
                communities[c2] = -1
                modularity.append(graph.modularity(ms))

        if membership is None:
            maxmod = max(modularity)
            maxidx = modularity.index(maxmod)
            membership = range(graph.vcount())
            communities = range(graph.vcount())
            midx = 0
            while maxidx>0:
                maxidx -= 1
                c1, c2 = merges[midx]
                midx += 1
                try:
                    cidx1 = communities[c1]
                    cidx2 = communities[c2]
                except IndexError:
                    raise ValueError, "invalid merge matrix, referencing nonexisting community"

                if cidx1 == -1 or cidx2 == -1:
                    raise ValueError, "invalid merge matrix, referencing already joined community"

                for idx, m in enumerate(membership):
                    if m == cidx2: membership[idx] = cidx1

                communities.append(communities[c1])
                communities[c1] = -1
                communities[c2] = -1

            recoding = {}
            n=0
            for idx, m in enumerate(membership):
                try:
                    v = recoding[m]
                except KeyError:
                    recoding[m], v = n, n
                    n += 1
                membership[idx] = v
        
        Dendrogram.__init__(self, merges)
        VertexClustering.__init__(self, graph, membership, None, params)


    def cut(self, n):
        """Cuts the dendrogram at a given level.

        @param n: the desired number of clusters. Merges are replayed from the
          beginning until the membership vector has exactly M{n} distinct elements
          or until there are no more recorded merges, whichever happens first.
        @return: the membership vector
        """
        num_elts = self._graph.vcount()
        membership = range(self._graph.vcount())
        communities = range(self._graph.vcount())
        midx = 0
        while num_elts>n:
            num_elts -= 1
            c1, c2 = self._merges[midx]
            midx += 1
            try:
                cidx1 = communities[c1]
                cidx2 = communities[c2]
            except IndexError:
                raise ValueError, "invalid merge matrix, referencing nonexisting community in row %d" % idx

            if cidx1 == -1 or cidx2 == -1:
                raise ValueError, "invalid merge matrix, referencing already joined community in row %d" % idx

            for idx, m in enumerate(membership):
                if m == cidx2: membership[idx] = cidx1

            communities.append(communities[c1])
            communities[c1] = -1
            communities[c2] = -1

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
        self._len = max(membership) - min(membership) + 1
        return copy(membership)

# vim:ts=4:sw=4:sts=4:et
"""
IGraph library.

@undocumented: igraph.formula, igraph.test
"""

__license__ = """
Copyright (C) 2006-2009  Gabor Csardi <csardi@rmki.kfki.hu>,
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

from igraph.core import *
from igraph.core import __version__, __build_date__
from igraph.clustering import *
from igraph.configuration import Configuration
from igraph.colors import *
from igraph.drawing import *
from igraph.datatypes import *
from igraph.formula import *
from igraph.layout import *

import os
import math
import gzip
import sys
import operator
from tempfile import mkstemp
from warnings import warn

class Graph(core.GraphBase):
    """Generic graph.
    
    This class is built on top of L{GraphBase}, so the order of the
    methods in the Epydoc documentation is a little bit obscure:
    inherited methods come after the ones implemented directly in the
    subclass. L{Graph} provides many functions that L{GraphBase} does not,
    mostly because these functions are not speed critical and they were
    easier to implement in Python than in pure C. An example is the
    attribute handling in the constructor: the constructor of L{Graph}
    accepts three dictionaries corresponding to the graph, vertex and edge
    attributes while the constructor of L{GraphBase} does not. This extension
    was needed to make L{Graph} serializable through the C{pickle} module.
    """

    # Some useful aliases
    omega = core.GraphBase.clique_number
    alpha = core.GraphBase.independence_number
    shell_index = core.GraphBase.coreness
    cut_vertices = core.GraphBase.articulation_points
    blocks = core.GraphBase.biconnected_components
    evcent = core.GraphBase.eigenvector_centrality
    vertex_disjoint_paths = core.GraphBase.vertex_connectivity
    edge_disjoint_paths = core.GraphBase.edge_connectivity
    cohesion = core.GraphBase.vertex_connectivity
    adhesion = core.GraphBase.edge_connectivity

    # Compatibility aliases
    shortest_paths_dijkstra = core.GraphBase.shortest_paths

    def __init__(self, *args, **kwds):
        """__init__(n=None, edges=None, directed=None, graph_attrs=None, vertex_attrs=None, edge_attrs=None)
        
        Constructs a graph from scratch.

        @keyword n: the number of vertices. Can be omitted.
        @keyword edges: the edge list where every list item is a pair of integers.
          If any of the integers is larger than M{n-1}, the number of vertices
          is adjusted accordingly.
        @keyword directed: whether the graph should be directed
        @keyword graph_attrs: the attributes of the graph as a dictionary.
        @keyword vertex_attrs: the attributes of the vertices as a dictionary.
          Every dictionary value must be an iterable with exactly M{n} items.
        @keyword edge_attrs: the attributes of the edges as a dictionary. Every
          dictionary value must be an iterable with exactly M{m} items where
          M{m} is the number of edges.
        """
        # Set up default values for the parameters. This should match the order in *args
        kwd_order = ["n", "edges", "directed", "graph_attrs", "vertex_attrs", "edge_attrs"]
        params = [1, [], False, {}, {}, {}]
        # If the first argument is a list, assume that the number of vertices were omitted
        args = list(args)
        if len(args) > 0:
            if isinstance(args[0], list) or isinstance(args[0], tuple):
                args.insert(0, params[0])
        # Override default parameters from args
        params[:len(args)] = args
        # Override default parameters from keywords
        for idx, k in enumerate(kwd_order):
            if k in kwds: params[idx] = kwds[k]
        # Now, translate the params list to argument names
        n, edges, directed, graph_attrs, vertex_attrs, edge_attrs = params

        # Initialize the graph
        core.GraphBase.__init__(self, n, edges, directed)
        # Set the graph attributes
        for k, v in graph_attrs.iteritems():
            if isinstance(k, int) or isinstance(k, long): k=str(k)
            self[k]=v
        # Set the vertex attributes
        for k, v in vertex_attrs.iteritems():
            if isinstance(k, int) or isinstance(k, long): k=str(k)
            self.vs[k]=v
        # Set the edge attributes
        for k, v in edge_attrs.iteritems():
            if isinstance(k, int) or isinstance(k, long): k=str(k)
            self.es[k]=v

    def delete_edges(self, *args, **kwds):
        """Deletes some edges from the graph.

        The set of edges to be deleted is determined by the positional and
        keyword arguments. If any keyword argument is present, or the
        first positional argument is callable, an edge
        sequence is derived by calling L{EdgeSeq.select} with the same
        positional and keyword arguments. Edges in the derived edge sequence
        will be removed. Otherwise the first positional argument is considered
        as follows:

          - C{None} - deletes all edges
          - a single integer - deletes the edge with the given ID
          - a list of integers - deletes the edges denoted by the given IDs
          - a list of 2-tuples - deletes the edges denoted by the given
            source-target vertex pairs. When multiple edges are present
            between a given source-target vertex pair, only one is removed.
        """
        if len(args) == 0 and len(kwds) == 0:
            raise ValueError, "expected at least one argument"
        if len(kwds)>0 or (callable(args[0]) and not isinstance(args[0], core.EdgeSeq)):
            es = self.es(*args, **kwds)
        else:
            es = args[0]
        return core.GraphBase.delete_edges(self, es)


    def indegree(self, *args, **kwds):
        """Returns the in-degrees in a list.
        
        See L{degree} for possible arguments.
        """
        kwds['type']=IN
        return self.degree(*args, **kwds)

    def outdegree(self, *args, **kwds):
        """Returns the out-degrees in a list.
        
        See L{degree} for possible arguments.
        """
        kwds['type']=OUT
        return self.degree(*args, **kwds)

    def biconnected_components(self, return_articulation_points=False):
        """biconnected_components(return_articulation_points=False)

        Calculates the biconnected components of the graph.
        
        @param return_articulation_points: whether to return the articulation
          points as well
        @return: an L{OverlappingVertexClustering} object describing the
          biconnected components and optionally the list of articulation points
        """
        if return_articulation_points:
            trees, aps = GraphBase.biconnected_components(self, True)
        else:
            trees = GraphBase.biconnected_components(self, False)

        membership = [set() for _ in xrange(self.vcount())]
        for idx, tree in enumerate(trees):
            for e in tree:
                membership[self.es[e].source].add(idx)
                membership[self.es[e].target].add(idx)
        cl = OverlappingVertexClustering(self, membership)

        if return_articulation_points:
            return cl, aps
        else:
            return cl
    blocks = biconnected_components

    def clusters(self, mode=STRONG):
        """clusters(mode=STRONG)

        Calculates the (strong or weak) clusters (connected components) for
        a given graph.

        @param mode: must be either C{STRONG} or C{WEAK}, depending on the
          clusters being sought. Optional, defaults to C{STRONG}.
        @return: a L{VertexClustering} object"""
        return VertexClustering(self, GraphBase.clusters(self, mode))
    components = clusters

    def degree_distribution(self, bin_width = 1, *args, **kwds):
        """degree_distribution(bin_width=1, ...)

        Calculates the degree distribution of the graph.

        Unknown keyword arguments are directly passed to L{degree()}.

        @param bin_width: the bin width of the histogram
        @return: a histogram representing the degree distribution of the
          graph.
        """
        result = Histogram(bin_width, self.degree(*args, **kwds))
        return result

    def dyad_census(self, *args, **kwds):
        """dyad_census()

        Calculates the dyad census of the graph.

        Dyad census means classifying each pair of vertices of a directed
        graph into three categories: mutual (there is an edge from I{a} to
        I{b} and also from I{b} to I{a}), asymmetric (there is an edge
        from I{a} to I{b} or from I{b} to I{a} but not the other way round)
        and null (there is no connection between I{a} and I{b}).

        @return: a L{DyadCensus} object.
        @newfield ref: Reference
        @ref: Holland, P.W. and Leinhardt, S.  (1970).  A Method for Detecting
          Structure in Sociometric Data.  American Journal of Sociology, 70,
          492-513.
        """
        return DyadCensus(GraphBase.dyad_census(self, *args, **kwds))

    def eccentricity(self, vertices=None):
        """Calculates eccentricities for vertices with the given indices.
        
        Eccentricity is given as the reciprocal of the greatest distance
        between the vertex being considered and any other vertex in the
        graph.

        Please note that for any unconnected graph, eccentricities will
        all be equal to 1 over the number of vertices, since for all vertices
        the greatest distance will be equal to the number of vertices (this
        is how L{shortest_paths} denotes vertex pairs where it is impossible
        to reach one from the other).

        @param vertices: the vertices to consider. If C{None}, all
          vertices are considered.
        @return: the eccentricities in a list
        """
        if self.vcount() == 0: return []
        if self.vcount() == 1: return [1.0]
        distance_matrix = self.shortest_paths(mode=OUT)
        distance_maxs = map(max, distance_matrix)
        
        if vertices is None:
            result = [1.0/x for x in distance_maxs]
        else:
            result = [1.0/distance_maxs[idx] for idx in vertices]

        return result

    def get_adjacency(self, type=GET_ADJACENCY_BOTH, attribute=None, default=None):
        """Returns the adjacency matrix of a graph.

        @param type: either C{GET_ADJACENCY_LOWER} (uses the lower
          triangle of the matrix) or C{GET_ADJACENCY_UPPER}
          (uses the upper triangle) or C{GET_ADJACENCY_BOTH}
          (uses both parts). Ignored for directed graphs.
        @param attribute: if C{None}, returns the ordinary adjacency
          matrix. When the name of a valid edge attribute is given
          here, the matrix returned will contain the default value 
          at the places where there is no edge or the value of the
          given attribute where there is an edge. Multiple edges are
          not supported, the value written in the matrix in this case
          will be unpredictable.
        @param default: the default value written to the cells in the
          case of adjacency matrices with attributes.
        @return: the adjacency matrix as a L{Matrix}.
        """
        if type != GET_ADJACENCY_LOWER and type != GET_ADJACENCY_UPPER and \
          type != GET_ADJACENCY_BOTH:
            # Maybe it was called with the first argument as the attribute name
            type, attribute = attribute, type
            if type is None: type = GET_ADJACENCY_BOTH

        if attribute is None: return Matrix(GraphBase.get_adjacency(self, type))
        if attribute not in self.es.attribute_names():
            raise ValueError, "Attribute does not exist"
        data = [[default] * self.vcount() for _ in xrange(self.vcount())]

        if not self.is_directed():
            if type == GET_ADJACENCY_BOTH:
                for e in self.es:
                    data[e.tuple[0]][e.tuple[1]] = e[attribute]
                    data[e.tuple[1]][e.tuple[0]] = e[attribute]
            elif type == GET_ADJACENCY_UPPER:
                for e in self.es: data[min(e.tuple)][max(e.tuple)] = e[attribute]
            else:
                for e in self.es: data[max(e.tuple)][min(e.tuple)] = e[attribute]
        else:
            for e in self.es: data[e.tuple[0]][e.tuple[1]] = e[attribute]

        return Matrix(data)


    def get_adjlist(self, type=OUT):
        """get_adjlist(type=OUT)

        Returns the adjacency list representation of the graph.

        The adjacency list representation is a list of lists. Each item of the
        outer list belongs to a single vertex of the graph. The inner list
        contains the neighbors of the given vertex.

        @param type: if L{OUT}, returns the successors of the vertex. If
          L{IN}, returns the predecessors of the vertex. If L{ALL}, both
          the predecessors and the successors will be returned. Ignored
          for undirected graphs.
        """
        return [self.neighbors(idx, type) for idx in xrange(self.vcount())]


    def get_adjedgelist(self, type=OUT):
        """get_adjedgelist(type=OUT)

        Returns the adjacency edge list representation of the graph.

        The adjacency edge list representation is a list of lists. Each
        item of the outer list belongs to a single vertex of the graph.
        The inner list contains the IDs of the adjacent edges of the
        given vertex.

        @param type: if L{OUT}, returns the successors of the vertex. If
          L{IN}, returns the predecessors of the vertex. If L{ALL}, both
          the predecessors and the successors will be returned. Ignored
          for undirected graphs.
        """
        return [self.adjacent(idx, type) for idx in xrange(self.vcount())]


    def mincut(self, capacity=None):
        """mincut(capacity=None)

        Returns a minimum cut in an undirected graph.

        @param capacity: the edge capacities (weights). If C{None}, all
          edges have equal weight.
        @return: a L{Cut} object describing the minimum cut
        """
        return Cut(self, *GraphBase.mincut(self, capacity))

    def modularity(self, membership, weights=None):
        """modularity(membership, weights=None)

        Calculates the modularity score of the graph with respect to a given
        clustering.
        
        The modularity of a graph w.r.t. some division measures how good the
        division is, or how separated are the different vertex types from each
        other. It is defined as M{Q=1/(2m)*sum(Aij-ki*kj/(2m)delta(ci,cj),i,j)}.
        M{m} is the number of edges, M{Aij} is the element of the M{A} adjacency
        matrix in row M{i} and column M{j}, M{ki} is the degree of node M{i},
        M{kj} is the degree of node M{j}, and M{Ci} and C{cj} are the types of
        the two vertices (M{i} and M{j}). M{delta(x,y)} is one iff M{x=y}, 0
        otherwise.
        
        If edge weights are given, the definition of modularity is modified as
        follows: M{Aij} becomes the weight of the corresponding edge, M{ki}
        is the total weight of edges adjacent to vertex M{i}, M{kj} is the
        total weight of edges adjacent to vertex M{j} and M{m} is the total
        edge weight in the graph.
        
        @param membership: a membership list or a L{VertexClustering} object
        @param weights: optional edge weights or C{None} if all edges are weighed
          equally. Attribute names are also allowed.
        @return: the modularity score
        
        @newfield ref: Reference
        @ref: MEJ Newman and M Girvan: Finding and evaluating community
          structure in networks. Phys Rev E 69 026113, 2004.
        """
        if isinstance(membership, VertexClustering):
            if membership.graph != self:
                raise ValueError, "clustering object belongs to a different graph"
            return GraphBase.modularity(self, membership.membership, weights)
        else:
            return GraphBase.modularity(self, membership, weights)

    def path_length_hist(self, directed=True):
        """path_length_hist(directed=True)

        Returns the path length histogram of the graph

        @param directed: whether to consider directed paths. Ignored for
          undirected graphs.
        @return: a L{Histogram} object. The object will also have an
          C{unconnected} attribute that stores the number of unconnected
          vertex pairs (where the second vertex can not be reached from
          the first one). The latter one will be of type long (and not
          a simple integer), since this can be I{very} large.
        """
        data, unconn = GraphBase.path_length_hist(self, directed)
        hist = Histogram(bin_width=1)
        for i, n in enumerate(data): hist.add(i+1, n)
        hist.unconnected = long(unconn)
        return hist

    def triad_census(self, *args, **kwds):
        """triad_census()

        Calculates the triad census of the graph.

        @return: a L{TriadCensus} object.
        @newfield ref: Reference
        @ref: Davis, J.A. and Leinhardt, S.  (1972).  The Structure of
          Positive Interpersonal Relations in Small Groups.  In:
          J. Berger (Ed.), Sociological Theories in Progress, Volume 2,
          218-251. Boston: Houghton Mifflin.
        """
        return TriadCensus(GraphBase.triad_census(self, *args, **kwds))

    # Automorphisms
    def count_automorphisms_vf2(self):
        """Returns the number of automorphisms of the graph"""
        return self.count_isomorphisms_vf2(self)

    def get_automorphisms_vf2(self):
        """Returns all automorphisms of the graph
        
        @return: a list of lists, each item containing a possible mapping
          of the graph vertices to itself according to the automorphism"""
        return self.get_isomorphisms_vf2(self)

    # Various clustering algorithms -- mostly wrappers around GraphBase
    def community_fastgreedy(self, weights=None):
        """Community structure based on the greedy optimization of modularity.

        This algorithm merges individual nodes into communities in a way that
        greedily maximizes the modularity score of the graph. It can be proven
        that if no merge can increase the current modularity score, the algorithm
        can be stopped since no further increase can be achieved.

        This algorithm is said to run almost in linear time on sparse graphs.

        @param weights: edge attribute name or a list containing edge
          weights
        @return: an appropriate L{VertexDendrogram} object.

        @newfield ref: Reference
        @ref: A Clauset, MEJ Newman and C Moore: Finding community structure
          in very large networks. Phys Rev E 70, 066111 (2004).
        """
        merges, qs = GraphBase.community_fastgreedy(self, weights, True)
        return VertexDendrogram(self, merges, None, qs)


    def community_leading_eigenvector_naive(self, clusters=None, return_merges = False):
        """community_leading_eigenvector_naive(clusters=None, return_merges=False)
        A naive implementation of Newman's eigenvector community structure
        detection. This function splits the network into two components
        according to the leading eigenvector of the modularity matrix and
        then recursively takes the given number of steps by splitting the
        communities as individual networks. This is not the correct way,
        however, see the reference for explanation. Consider using the
        correct L{community_leading_eigenvector} method instead.

        @param clusters: the desired number of communities. If C{None}, the algorithm
          tries to do as many splits as possible. Note that the algorithm
          won't split a community further if the signs of the leading eigenvector
          are all the same, so the actual number of discovered communities can be
          less than the desired one.
        @param return_merges: whether the returned object should be a
          dendrogram instead of a single clustering.
        @return: an appropriate L{VertexClustering} or L{VertexDendrogram}
          object.
        
        @newfield ref: Reference
        @ref: MEJ Newman: Finding community structure in networks using the
        eigenvectors of matrices, arXiv:physics/0605087"""
        if clusters is None: clusters=-1
        cl, merges = GraphBase.community_leading_eigenvector_naive(self, clusters, return_merges)
        if merges is None:
            return VertexClustering(self, cl)
        else:
            return VertexDendrogram(self, merges, cl)


    def community_leading_eigenvector(self, clusters=None, return_merges = False):
        """community_leading_eigenvector(clusters=None, return_merges=False)
        
        Newman's leading eigenvector method for detecting community structure.
        This is the proper implementation of the recursive, divisive algorithm:
        each split is done by maximizing the modularity regarding the
        original network.
        
        @param clusters: the desired number of communities. If C{None}, the algorithm
          tries to do as many splits as possible. Note that the algorithm
          won't split a community further if the signs of the leading eigenvector
          are all the same, so the actual number of discovered communities can be
          less than the desired one.
        @param return_merges: whether the returned object should be a
          dendrogram instead of a single clustering.
        @return: an appropriate L{VertexClustering} or L{VertexDendrogram}
          object.
        
        @newfield ref: Reference
        @ref: MEJ Newman: Finding community structure in networks using the
        eigenvectors of matrices, arXiv:physics/0605087"""
        if clusters is None: clusters=-1
        cl, merges = GraphBase.community_leading_eigenvector(self, clusters, return_merges)
        if merges is None:
            return VertexClustering(self, cl)
        else:
            return VertexDendrogram(self, merges, cl)


    def community_label_propagation(self, weights = None, initial = None, fixed = None):
        """community_label_propagation(weights=None, initial=None, fixed=None)

        Finds the community structure of the graph according to the label
        propagation method of Raghavan et al.
        Initially, each vertex is assigned a different label. After that,
        each vertex chooses the dominant label in its neighbourhood in each
        iteration. Ties are broken randomly and the order in which the
        vertices are updated is randomized before every iteration. The
        algorithm ends when vertices reach a consensus.
        Note that since ties are broken randomly, there is no guarantee that
        the algorithm returns the same community structure after each run.
        In fact, they frequently differ. See the paper of Raghavan et al
        on how to come up with an aggregated community structure.
        @param weights: name of an edge attribute or a list containing
          edge weights
        @param initial: name of a vertex attribute or a list containing
          the initial vertex labels. Labels are identified by integers from
          zero to M{n-1} where M{n} is the number of vertices. Negative
          numbers may also be present in this vector, they represent unlabeled
          vertices.
        @param fixed: a list of booleans for each vertex. C{True} corresponds
          to vertices whose labeling should not change during the algorithm.
          It only makes sense if initial labels are also given. Unlabeled
          vertices cannot be fixed.
        @return: the resulting membership vector

        @newfield ref: Reference
        @ref: Raghavan, U.N. and Albert, R. and Kumara, S. Near linear
          time algorithm to detect community structures in large-scale
          networks. Phys Rev E 76:036106, 2007.
          U{http://arxiv.org/abs/0709.2938}.
        """
        if isinstance(fixed, (str, unicode)):
            fixed = [bool(o) for o in g.vs[fixed]]
        cl = GraphBase.community_label_propagation(self, weights, initial, fixed)
        return VertexClustering(self, cl)


    def community_edge_betweenness(self, clusters = None, directed = True):
        """Community structure based on the betweenness of the edges in the network.

        The idea is that the betweenness of the edges connecting two communities
        is typically high, as many of the shortest paths between nodes in separate
        communities go through them. So we gradually remove the edge with the
        highest betweenness and recalculate the betweennesses after every
        removal. This way sooner or later the network falls of to separate
        components. The result of the clustering will be represented by a
        dendrogram.

        @param clusters: the number of clusters we would like to see. This
          practically defines the "level" where we "cut" the dendrogram to
          get the membership vector of the vertices. If C{None}, the dendrogram
          is cut at the level which maximizes the modularity.
        @param directed: whether the directionality of the edges should be taken
          into account or not.
        @return: a L{VertexDendrogram} object, initally cut at the maximum
          modularity.
        """
        d = VertexDendrogram(self, GraphBase.community_edge_betweenness(self, directed));
        if clusters is not None: d.cut(clusters)
        return d

    def community_walktrap(self, weights=None, steps=4):
        """Community detection algorithm of Latapy & Pons, based on random walks.
        
        The basic idea of the algorithm is that short random walks tend to stay
        in the same community. The result of the clustering will be represented
        as a dendrogram.
        
        @param weights: name of an edge attribute or a list containing
          edge weights
        @param steps: length of random walks to perform
        
        @return: a L{VertexDendrogram} object, initially cut at the maximum
          modularity.
          
        @newfield ref: Reference
        @ref: Pascal Pons, Matthieu Latapy: Computing communities in large networks
          using random walks, U{http://arxiv.org/abs/physics/0512106}.
        """
        merges, qs = GraphBase.community_walktrap(self, weights, steps, True)
        d = VertexDendrogram(self, merges, modularity=qs)
        return d


    def k_core(self, *args):
        """Returns some k-cores of the graph.

        The method accepts an arbitrary number of arguments representing
        the desired indices of the M{k}-cores to be returned. The arguments
        can also be lists or tuples. The result is a single L{Graph} object
        if an only integer argument was given, otherwise the result is a
        list of L{Graph} objects representing the desired k-cores in the
        order the arguments were specified. If no argument is given, returns
        all M{k}-cores in increasing order of M{k}.
        """
        if len(args) == 0:
            indices = xrange(self.vcount())
            return_single = False
        else:
            return_single = True
            indices = []
            for arg in args:
                try:
                    indices.extend(arg)
                except:
                    indices.append(arg)

        if len(indices)>1 or hasattr(args[0], "__iter__"):
            return_single = False

        corenesses = self.coreness()
        result = []
        vidxs = xrange(self.vcount())
        for idx in indices:
            core_idxs = [vidx for vidx in vidxs if corenesses[vidx] >= idx]
            result.append(self.subgraph(core_idxs))

        if return_single: return result[0]
        return result


    def layout(self, layout=None, *args, **kwds):
        """Returns the layout of the graph according to a layout algorithm.

        Parameters and keyword arguments not specified here are passed to the
        layout algorithm directly. See the documentation of the layout
        algorithms for the explanation of these parameters.

        Registered layout names understood by this method are:

          - C{circle}, C{circular}: circular layout (see L{Graph.layout_circle})

          - C{drl}: DrL layout for large graphs (see L{Graph.layout_drl})

          - C{fr}, C{fruchterman_reingold}: Fruchterman-Reingold layout
            (see L{Graph.layout_fruchterman_reingold}).

          - C{fr_3d}, C{fr3d}, C{fruchterman_reingold_3d}: 3D Fruchterman-Reingold
            layout (see L{Graph.layout_fruchterman_reingold_3d}).

          - C{graphopt}: the graphopt algorithm (see L{Graph.layout_graphopt})

          - C{gfr}, C{grid_fr}, C{grid_fruchterman_reingold}: grid-based
            Fruchterman-Reingold layout
            (see L{Graph.layout_grid_fruchterman_reingold})

          - C{kk}, C{kamada_kawai}: Kamada-Kawai layout
            (see L{Graph.layout_kamada_kawai})

          - C{kk_3d}, C{kk3d}, C{kamada_kawai_3d}: 3D Kamada-Kawai layout
            (see L{Graph.layout_kamada_kawai_3d})

          - C{lgl}, C{large}, C{large_graph}: Large Graph Layout
            (see L{Graph.layout_lgl})

          - C{random}: random layout (see L{Graph.layout_random})

          - C{random_3d}: random 3D layout (see L{Graph.layout_random_3d})

          - C{rt}, C{tree}, C{reingold_tilford}: Reingold-Tilford tree
            layout (see L{Graph.layout_reingold_tilford})

          - C{rt_circular}, C{reingold_tilford_circular}: circular
            Reingold-Tilford tree layout
            (see L{Graph.layout_reingold_tilford_circular})

          - C{sphere}, C{spherical}, C{circle_3d}, C{circular_3d}: spherical
            layout (see L{Graph.layout_sphere})

        @param layout: the layout to use. This can be one of the registered
          layout names or a callable which returns either a L{Layout} object or
          a list of lists containing the coordinates. If C{None}, uses the
          value of the C{plotting.layout} configuration key.
        @return: a L{Layout} object.
        """
        if layout is None: layout = config["plotting.layout"]
        if callable(layout):
            method = layout
        else:
            method = getattr(self.__class__, self._layout_mapping[layout.lower()])
        if not callable(method): raise ValueError, "layout method must be callable"
        l=method(self, *args, **kwds)
        if not isinstance(l, Layout): l=Layout(l)
        return l

    # Auxiliary I/O functions

    def write_adjacency(self, f, sep=" ", eol="\n", *args, **kwds):
        """Writes the adjacency matrix of the graph to the given file

        All the remaining arguments not mentioned here are passed intact
        to L{Graph.get_adjacency}.

        @param f: the name of the file to be written.
        @param sep: the string that separates the matrix elements in a row
        @param eol: the string that separates the rows of the matrix. Please
          note that igraph is able to read back the written adjacency matrix
          if and only if this is a single newline character
        """
        if not isinstance(f, file): f = file(f, "w")
        matrix = self.get_adjacency(*args, **kwds)
        for row in matrix:
            f.write(sep.join(map(str, row)))
            f.write(eol)
        f.close()

    @classmethod
    def Read_Adjacency(klass, f, sep=None, comment_char = "#", attribute=None,
        *args, **kwds):
        """Constructs a graph based on an adjacency matrix from the given file

        Additional positional and keyword arguments not mentioned here are
        passed intact to L{Graph.Adjacency}.

        @param f: the name of the file to be read or a file object
        @param sep: the string that separates the matrix elements in a row.
          C{None} means an arbitrary sequence of whitespace characters.
        @param comment_char: lines starting with this string are treated
          as comments.
        @param attribute: an edge attribute name where the edge weights are
          stored in the case of a weighted adjacency matrix. If C{None},
          no weights are stored, values larger than 1 are considered as
          edge multiplicities.
        @return: the created graph"""
        if not isinstance(f, file): f = file(f)
        matrix, ri, weights = [], 0, {} 
        for line in f:
            line = line.strip()
            if len(line) == 0: continue
            if line.startswith(comment_char): continue
            row = map(float, line.split(sep))
            matrix.append(row)
            ri += 1

        f.close()

        if attribute is None:
            graph=klass.Adjacency(matrix, *args, **kwds)
        else:
            kwds["attr"] = attribute
            graph=klass.Weighted_Adjacency(matrix, *args, **kwds)

        return graph

    def write_graphmlz(self, f, compresslevel=9):
        """Writes the graph to a zipped GraphML file.

        The library uses the gzip compression algorithm, so the resulting
        file can be unzipped with regular gzip uncompression (like
        C{gunzip} or C{zcat} from Unix command line) or the Python C{gzip}
        module.

        Uses a temporary file to store intermediate GraphML data, so
        make sure you have enough free space to store the unzipped
        GraphML file as well.

        @param f: the name of the file to be written.
        @param compresslevel: the level of compression. 1 is fastest and
          produces the least compression, and 9 is slowest and produces
          the most compression."""
        tmpfilename=None
        try:
            tmpfileno, tmpfilename = mkstemp(text=True)
            os.close(tmpfileno)
            self.write_graphml(tmpfilename)
            outf = gzip.GzipFile(f, "wb", compresslevel)
            inf = open(tmpfilename)
            for line in inf: outf.write(line)
            inf.close()
            outf.close()
        finally:
            if tmpfilename is not None: os.unlink(tmpfilename)

    @classmethod
    def Read_GraphMLz(cls, f, *params, **kwds):
        """Read_GraphMLz(f, directed=True, index=0)

        Reads a graph from a zipped GraphML file.

        @param f: the name of the file
        @param index: if the GraphML file contains multiple graphs,
          specified the one that should be loaded. Graph indices
          start from zero, so if you want to load the first graph,
          specify 0 here.
        @return: the loaded graph object"""
        tmpfilename=None
        try:
            tmpfileno, tmpfilename = mkstemp(text=True)
            os.close(tmpfileno)
            inf = gzip.GzipFile(f, "rb")
            outf = open(tmpfilename, "wt")
            for line in inf: outf.write(line)
            inf.close()
            outf.close()
            result=cls.Read_GraphML(tmpfilename)
        finally:
            if tmpfilename is not None: os.unlink(tmpfilename)
        return result


    def write_pickle(self, fname=None, version=-1):
        """Saves the graph in Python pickled format

        @param fname: the name of the file or a stream to save to. If
          C{None}, saves the graph to a string and returns the string.
        @param version: pickle protocol version to be used. If -1, uses
          the highest protocol available
        @return: C{None} if the graph was saved successfully to the
          file given, or a string if C{fname} was C{None}.
        """
        import cPickle
        if fname is None: return cPickle.dumps(self, version)
        if not isinstance(fname, file):
            file_was_opened=True
            fname=open(fname, 'wb')
        result=cPickle.dump(self, fname, version)
        if file_was_opened: fname.close()
        return result

    @classmethod
    def Read_Pickle(klass, fname=None):
        """Reads a graph from Python pickled format

        @param fname: the name of the file, a stream to read from, or
          a string containing the pickled data. The string is assumed to
          hold pickled data if it is longer than 40 characters and
          contains a substring that's peculiar to pickled versions
          of an C{igraph} Graph object.
        @return: the created graph object.
        """
        import cPickle
        if len(fname)>40 and "cigraph\nGraph\n" in fname:
            return cPickle.loads(fname)
        if not isinstance(fname, file):
            file_was_opened=True
            fname=open(fname, 'rb')
        result=cPickle.load(fname)
        if file_was_opened: fname.close()
        if not isinstance(result, klass):
            raise TypeError, "unpickled object is not a %s" % klass.__name__
        return result

    def write_svg(self, fname, layout, width = None, height = None, \
                  labels = "label", colors = "color", shapes = "shape", \
                  vertex_size = 10, edge_colors = "color", \
                  font_size = 16, *args, **kwds):
        """Saves the graph as an SVG (Scalable Vector Graphics) file
        
        @param fname: the name of the file
        @param layout: the layout of the graph. Can be either an
          explicitly specified layout (using a list of coordinate
          pairs) or the name of a layout algorithm (which should
          refer to a method in the L{Graph} object, but without
          the C{layout_} prefix.
        @param width: the preferred width in pixels (default: 400)
        @param height: the preferred height in pixels (default: 400)
        @param labels: the vertex labels. Either it is the name of
          a vertex attribute to use, or a list explicitly specifying
          the labels. It can also be C{None}.
        @param colors: the vertex colors. Either it is the name of
          a vertex attribute to use, or a list explicitly specifying
          the colors. A color can be anything acceptable in an SVG
          file.
        @param shapes: the vertex shapes. Either it is the name of
          a vertex attribute to use, or a list explicitly specifying
          the shapes as integers. Shape 0 means hidden (nothing is drawn),
          shape 1 is a circle, shape 2 is a rectangle.
        @param vertex_size: vertex size in pixels
        @param edge_colors: the edge colors. Either it is the name
          of an edge attribute to use, or a list explicitly specifying
          the colors. A color can be anything acceptable in an SVG
          file.
        @param font_size: font size. If it is a string, it is written into
          the SVG file as-is (so you can specify anything which is valid
          as the value of the C{font-size} style). If it is a number, it
          is interpreted as pixel size and converted to the proper attribute
          value accordingly.
        """
        if width is None and height is None:
            width = 400
            height = 400
        elif width is None:
            width = height
        elif height is None:
            height = width
                
        if width<=0 or height<=0:
            raise ValueError, "width and height must be positive"

        if isinstance(layout, str):
            f=getattr(Graph, "layout_"+layout);
            layout=f(self, *args)

        if isinstance(labels, str):
            try:
                labels = self.vs.get_attribute_values(labels)
            except KeyError:
                labels = [x+1 for x in xrange(self.vcount())]
        elif labels is None:
            labels = [""] * self.vcount()

        if isinstance(colors, str):
            try:
                colors = self.vs.get_attribute_values(colors)
            except KeyError:
                colors = ["red" for x in xrange(self.vcount())]

        if isinstance(shapes, str):
            try:
                shapes = self.vs.get_attribute_values(shapes)
            except KeyError:
                shapes = [1]*self.vcount()
        if isinstance(edge_colors, str):
            try:
                edge_colors = self.es.get_attribute_values(edge_colors)
            except KeyError:
                edge_colors = ["black" for x in xrange(self.ecount())]

        if not isinstance(font_size, str):
            font_size = "%spx" % str(font_size)
        else:
            if ";" in font_size:
                raise ValueError, "font size can't contain a semicolon"

        vc = self.vcount()
        while len(labels)<vc: labels.append(len(labels)+1)
        while len(colors)<vc: colors.append("red")

        f=open(fname, "w")
                
        maxs=[layout[0][dim] for dim in range(2)]
        mins=[layout[0][dim] for dim in range(2)]
                
        for rowidx in range(1, len(layout)):
            row = layout[rowidx]
            for dim in range(0, 2):
                if maxs[dim]<row[dim]: maxs[dim]=row[dim]
                if mins[dim]>row[dim]: mins[dim]=row[dim]
                
        sizes=[width-2*vertex_size, height-2*vertex_size]
        halfsizes=[(maxs[dim]+mins[dim])/2.0 for dim in range(2)]
        ratios=[sizes[dim]/(maxs[dim]-mins[dim]) for dim in range(2)]
        layout=[[(row[0]-halfsizes[0])*ratios[0], \
                 (row[1]-halfsizes[1])*ratios[1]] \
                for row in layout]
                
        directed=self.is_directed()

        print >>f, "<?xml version=\"1.0\" standalone=\"no\"?>"
        print >>f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\""
        print >>f, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
        
        print >>f, "<svg width=\"%d\" height=\"%d\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\" xmlns:xlink=\"http://www.w3.org/1999/xlink\">" % (width, height)
        print >>f, "<!-- Created by igraph -->"
        print >>f
        print >>f, "<defs>"
        print >>f, "  <symbol id=\"Triangle\" overflow=\"visible\">"
        print >>f, "    <path d=\"M 0 0 L 10 -5 L 10 5 z\"/>"
        print >>f, "  </symbol>"
        print >>f, "  <style type=\"text/css\">"
        print >>f, "    <![CDATA["
        print >>f, "#vertices circle { stroke: black; stroke-width: 1 }"
        print >>f, "#vertices rect { stroke: black; stroke-width: 1 }"
        print >>f, "#vertices text { text-anchor: middle; font-size: %s; font-family: sans-serif; font-weight: normal }" % font_size
        print >>f, "#edges line { stroke-width: 1 }"
        print >>f, "    ]]>"
        print >>f, "  </style>"
        print >>f, "</defs>"
        print >>f
        print >>f, "<g transform=\"translate(%.4f,%.4f)\">" % (width/2.0, height/2.0)
        print >>f, "  <g id=\"edges\">"
        print >>f, "  <!-- Edges -->"

        has_edge_opacities = "opacity" in self.edge_attributes()
        for eidx, edge in enumerate(self.es):
            vidxs = edge.tuple
            x1 = layout[vidxs[0]][0]
            y1 = layout[vidxs[0]][1]
            x2 = layout[vidxs[1]][0]
            y2 = layout[vidxs[1]][1]
            angle = math.atan2(y2-y1, x2-x1)
            x2 = x2 - vertex_size*math.cos(angle)
            y2 = y2 - vertex_size*math.sin(angle)
            if directed:
                # Dirty hack because of the SVG specification:
                # markers do not inherit stroke colors
                print >>f, "    <g transform=\"translate(%.4f,%.4f)\" fill=\"%s\" stroke=\"%s\">" % (x2, y2, edge_colors[eidx], edge_colors[eidx]) 
                print >>f, "      <line x1=\"%.4f\" y1=\"%.4f\" x2=\"0\" y2=\"0\"/>" % (x1-x2, y1-y2)
                print >>f, "      <use x=\"0\" y=\"0\" xlink:href=\"#Triangle\" transform=\"rotate(%.4f)\"/>" % (180+angle*180/math.pi,)
                print >>f, "    </g>\n"
            else:
                print >>f, "    <line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\" style=\"stroke: %s\"/>" % (x1, y1, x2, y2, edge_colors[eidx])

        print >>f, "  </g>"
        print >>f

        print >>f, "  <g id=\"vertices\">"
        print >>f, "  <!-- Vertices -->"
        for vidx in range(self.vcount()):
            print >>f, "    <g transform=\"translate(%.4f %.4f)\">" % (layout[vidx][0], layout[vidx][1])
            if shapes[vidx] == 1:
                # Undocumented feature: can handle two colors
                c = str(colors[vidx])
                if " " in c:
                    c = c.split(" ")
                    vs = str(vertex_size)
                    print >>f, "      <path d=\"M -%s,0 A%s,%s 0 0,0 %s,0 L -%s,0\" fill=\"%s\"/>" % (vs,vs,vs,vs,vs,c[0])
                    print >>f, "      <path d=\"M -%s,0 A%s,%s 0 0,1 %s,0 L -%s,0\" fill=\"%s\"/>" % (vs,vs,vs,vs,vs,c[1])
                    print >>f, "      <circle cx=\"0\" cy=\"0\" r=\"%s\" fill=\"none\"/>" % vs
                else:
                    print >>f, "      <circle cx=\"0\" cy=\"0\" r=\"%s\" fill=\"%s\"/>" % (str(vertex_size), str(colors[vidx]))
            elif shapes[vidx] == 2:
                print >>f, "      <rect x=\"-%s\" y=\"-%s\" width=\"%s\" height=\"%s\" fill=\"%s\"/>" % (str(vertex_size), str(vertex_size), str(2*vertex_size), str(2*vertex_size), str(colors[vidx]))
            print >>f, "      <text x=\"0\" y=\"5\">%s</text>" % str(labels[vidx])
            print >>f, "    </g>"

        print >>f, "  </g>"
        print >>f, "</g>"
        print >>f
        print >>f, "</svg>"
                
        f.close()


    @classmethod
    def _identify_format(klass, filename):
        """_identify_format(filename)

        Tries to identify the format of the graph stored in the file with the
        given filename. It identifies most file formats based on the extension
        of the file (and not on syntactic evaluation). The only exception is
        the adjacency matrix format and the edge list format: the first few
        lines of the file are evaluated to decide between the two.

        @note: Internal function, should not be called directly.

        @param filename: the name of the file or a file object whose C{name}
          attribute is set.
        @return: the format of the file as a string.
        """
        import os.path
        if isinstance(filename, file):
            try:
                filename=filename.name
            except:
                return None

        root, ext = os.path.splitext(filename)
        ext = ext.lower()
        
        if ext in [".graphml", ".graphmlz", ".lgl", ".ncol", ".pajek",
            ".gml", ".dimacs", ".edgelist", ".edges", ".edge", ".net",
            ".pickle", ".dot"]:
            return ext[1:]

        if ext == ".txt" or ext == ".dat":
            # Most probably an adjacency matrix or an edge list
            f = open(filename, "r")
            line = f.readline()
            if line is None: return "edges"
            parts = line.strip().split()
            if len(parts) == 2:
                line = f.readline()
                if line is None: return "edges"
                parts = line.strip().split()
                if len(parts) == 2:
                    line = f.readline()
                    if line is None:
                        # This is a 2x2 matrix, it can be a matrix or an edge list
                        return None
                    else:
                        parts = line.strip().split()
                        if len(parts) == 0:
                            return None
                    return "edges"
                else:
                    # Not a matrix
                    return None
            else:
                return "adjacency"

    @classmethod
    def Read(klass, f, format=None, *args, **kwds):
        """Unified reading function for graphs.

        This method tries to identify the format of the graph given in
        the first parameter and calls the corresponding reader method.

        The remaining arguments are passed to the reader method without
        any changes.

        @param f: the file containing the graph to be loaded
        @param format: the format of the file (if known in advance).
          C{None} means auto-detection. Possible values are: C{"ncol"}
          (NCOL format), C{"lgl"} (LGL format), C{"graphml"}, C{"graphmlz"}
          (GraphML and gzipped GraphML format), C{"gml"} (GML format),
          C{"net"}, C{"pajek"} (Pajek format), C{"dimacs"} (DIMACS format),
          C{"edgelist"}, C{"edges"} or C{"edge"} (edge list), C{"adjacency"}
          (adjacency matrix), C{"pickle"} (Python pickled format).
        @raises IOError: if the file format can't be identified and
          none was given.
        """
        if format is None: format = klass._identify_format(f)
        try:
            reader = klass._format_mapping[format][0]
        except KeyError, IndexError:
            raise IOError, "unknown file format: %s" % str(format)
        if reader is None:
            raise IOError, "no reader method for file format: %s" % str(format)
        reader = getattr(klass, reader)
        return reader(f, *args, **kwds)
    Load = Read

    
    def write(self, f, format=None, *args, **kwds):
        """Unified writing function for graphs.

        This method tries to identify the format of the graph given in
        the first parameter (based on extension) and calls the corresponding
        writer method.

        The remaining arguments are passed to the writer method without
        any changes.

        @param f: the file containing the graph to be saved
        @param format: the format of the file (if one wants to override the
          format determined from the filename extension, or the filename itself
          is a stream). C{None} means auto-detection. Possible values are: C{"ncol"}
          (NCOL format), C{"lgl"} (LGL format), C{"graphml"}, C{"graphmlz"}
          (GraphML and gzipped GraphML format), C{"gml"} (GML format),
          C{"dot"}, C{"graphviz"} (DOT format, used by GraphViz), C{"net"},
          C{"pajek"} (Pajek format), C{"dimacs"} (DIMACS format),
          C{"edgelist"}, C{"edges"} or C{"edge"} (edge list), C{"adjacency"}
          (adjacency matrix), C{"pickle"} (Python pickled format),
          C{"svg"} (Scalable Vector Graphics).
        @raises IOError: if the file format can't be identified and
          none was given.
        """
        if format is None: format = self._identify_format(f)
        try:
            writer = self._format_mapping[format][1]
        except KeyError, IndexError:
            raise IOError, "unknown file format: %s" % str(format)
        if writer is None:
            raise IOError, "no writer method for file format: %s" % str(format)
        writer = getattr(self, writer)
        return writer(f, *args, **kwds)
    save = write

    #####################################################
    # Constructor for dict-like representation of graphs

    @classmethod
    def DictList(klass, vertices, edges, directed=False, vertex_name_attr="name", \
            edge_foreign_keys=("source", "target"), iterative=False):
        """Constructs a graph from a list-of-dictionaries representation.

        This representation assumes that vertices and edges are encoded in
        two lists, each list containing a Python dict for each vertex and
        each edge, respectively. A distinguished element of the vertex dicts
        contain a vertex ID which is used in the edge dicts to refer to
        source and target vertices. All the remaining elements of the dict
        are considered vertex and edge attributes. Note that the implementation
        does not assume that the objects passed to this method are indeed
        lists of dicts, but they should be iterable and they should yield
        objects that behave as dicts. So, for instance, a database query
        result is likely to be fit as long as it's iterable and yields
        dict-like objects with every iteration.

        @param vertices: the data source for the vertices or C{None} if
          there are no special attributes assigned to vertices and we
          should simply use the edge list of dicts to infer vertex names.
        @param edges: the data source for the edges.
        @param directed: whether the constructed graph will be directed
        @param vertex_name_attr: the name of the distinguished key in the
          dicts in the vertex data source that contains the vertex names.
          Ignored if C{vertices} is C{None}.
        @param edge_foreign_keys: the name of the attributes in the dicts
          in the edge data source that contain the source and target
          vertex names.
        @param iterative: whether to add the edges to the graph one by one,
          iteratively, or to build a large edge list first and use that to
          construct the graph. The latter approach is faster but it may
          not be suitable if your dataset is large. The default is to
          add the edges in a batch from an edge list.
        @return: the graph that was constructed
        """
        def create_list_from_indices(l, n):
            result = [None] * n
            for i, v in l: result[i] = v
            return result

        # Construct the vertices
        vertex_attrs, n = {}, 0
        if vertices:
            for idx, vertex_data in enumerate(vertices):
                for k, v in vertex_data.iteritems():
                    try:
                        vertex_attrs[k].append((idx, v))
                    except KeyError:
                        vertex_attrs[k] = [(idx, v)]
                n += 1
            for k, v in vertex_attrs.iteritems():
                vertex_attrs[k] = create_list_from_indices(v, n)
        else:
            vertex_attrs[vertex_name_attr] = []

        vertex_names = vertex_attrs[vertex_name_attr]
        # Check for duplicates in vertex_names
        if len(vertex_names) != len(set(vertex_names)):
            raise ValueError, "vertex names are not unique"
        # Create a reverse mapping from vertex names to indices
        vertex_name_map = UniqueIdGenerator(initial = vertex_names)

        # Construct the edges
        efk_src, efk_dest = edge_foreign_keys
        if iterative:
            g = klass(n, [], directed, {}, vertex_attrs)
            for idx, edge_data in enumerate(edges):
                src_name, dst_name = edge_data[efk_src], edge_data[efk_dest]
                v1 = vertex_name_map[src_name]
                if v1 == n:
                    g.add_vertices(1)
                    g.vs[n][vertex_name_attr] = src_name
                    n += 1
                v2 = vertex_name_map[dst_name]
                if v2 == n:
                    g.add_vertices(1)
                    g.vs[n][vertex_name_attr] = dst_name
                    n += 1
                g.add_edges((v1, v2))
                for k, v in edge_data.iteritems():
                    g.es[idx][k] = v

            return g
        else:
            edge_list, edge_attrs, m = [], {}, 0
            for idx, edge_data in enumerate(edges):
                v1, v2 = vertex_name_map[edge_data[efk_src]], vertex_name_map[edge_data[efk_dest]]
                edge_list.append((v1, v2))
                for k, v in edge_data.iteritems():
                    try:
                        edge_attrs[k].append((idx, v))
                    except KeyError:
                        edge_attrs[k] = [(idx, v)]
                m += 1
            for k, v in edge_attrs.iteritems():
                edge_attrs[k] = create_list_from_indices(v, m)

            # It may have happened that some vertices were added during the process
            if len(vertex_name_map) > n:
                diff = len(vertex_name_map) - n
                more = [None] * diff
                for k, v in vertex_attrs.iteritems(): v.extend(more)
                vertex_attrs[vertex_name_attr] = vertex_name_map.values()
                n = len(vertex_name_map)

            # Create the graph
            return klass(n, edge_list, directed, {}, vertex_attrs, edge_attrs)

    #################################
    # Constructor for graph formulae
    Formula=classmethod(Formula)

    ###########################
    # Vertex and edge sequence
    def _get_vs(self): return VertexSeq(self)
    def _get_es(self): return EdgeSeq(self)
    vs=property(_get_vs, doc="The vertex sequence of the graph")
    es=property(_get_es, doc="The edge sequence of the graph")

    #############################################
    # Friendlier interface for bipartite methods
    @classmethod
    def Bipartite(klass, types, *args, **kwds):
        """Bipartite(types, edges, directed=False)

        Creates a bipartite graph with the given vertex types and edges.
        This is similar to the default constructor of the graph, the
        only difference is that it checks whether all the edges go
        between the two vertex classes and it assigns the type vector
        to a C{type} attribute afterwards.

        Examples:

        >>> g = Graph.Bipartite([0, 1, 0, 1], [(0, 1), (2, 3), (0, 3)])
        >>> g.is_bipartite()
        True
        >>> g.vs["type"]
        [False, True, False, True]

        @param types: the vertex types as a boolean list. Anything that
          evaluates to C{False} will denote a vertex of the first kind,
          anything that evaluates to C{True} will denote a vertex of the
          second kind.
        @param edges: the edges as a list of tuples.
        @param directed: whether to create a directed graph. Bipartite
          networks are usually undirected, so the default is C{False}

        @return: the graph with a binary vertex attribute named C{"type"} that
          stores the vertex classes.
        """
        result = klass._Bipartite(types, *args, **kwds)
        result.vs["type"] = [bool(x) for x in types]
        return result

    @classmethod
    def Full_Bipartite(klass, *args, **kwds):
        """Full_Bipartite(n1, n2, directed=False, mode=ALL)

        Generates a full bipartite graph (directed or undirected, with or
        without loops).

        >>> g = Graph.Full_Bipartite(2, 3)
        >>> g.is_bipartite()
        True
        >>> g.vs["type"]
        [False, False, True, True, True]

        @param n1: the number of vertices of the first kind.
        @param n2: the number of vertices of the second kind.
        @param directed: whether tp generate a directed graph.
        @param mode: if C{OUT}, then all vertices of the first kind are
          connected to the others; C{IN} specifies the opposite direction,
          C{ALL} creates mutual edges. Ignored for undirected graphs.

        @return: the graph with a binary vertex attribute named C{"type"} that
          stores the vertex classes.
        """
        result, types = klass._Full_Bipartite(*args, **kwds)
        result.vs["type"] = types
        return result

    @classmethod
    def Incidence(klass, *args, **kwds):
        """Incidence(matrix, directed=False, mode=ALL, multiple=False)

        Creates a bipartite graph from an incidence matrix.

        Example:

        >>> g = Graph.Incidence([[0, 1, 1], [1, 1, 0]])

        @param matrix: the incidence matrix.
        @param directed: whether to create a directed graph.
        @param mode: defines the direction of edges in the graph. If
          C{OUT}, then edges go from vertices of the first kind
          (corresponding to rows of the matrix) to vertices of the
          second kind (the columns of the matrix). If C{IN}, the
          opposite direction is used. C{ALL} creates mutual edges.
          Ignored for undirected graphs.
        @param multiple: defines what to do with non-zero entries in the
          matrix. If C{False}, non-zero entries will create an edge no matter
          what the value is. If C{True}, non-zero entries are rounded up to
          the nearest integer and this will be the number of multiple edges
          created.

        @return: the graph with a binary vertex attribute named C{"type"} that
          stores the vertex classes.
        """
        result, types = klass._Incidence(*args, **kwds)
        result.vs["type"] = types
        return result

    def bipartite_projection(self, types="type", *args, **kwds):
        """bipartite_projection(types="type", probe1=-1)

        Projects a bipartite graph into two one-mode graphs. Edge directions
        are ignored while projecting.

        Examples:

        >>> g = Graph.Full_Bipartite(10, 5)
        >>> g1, g2 = g.bipartite_projection()
        >>> g1.isomorphic(Graph.Full(10))
        True
        >>> g2.isomorphic(Graph.Full(5))
        True
        
        @param types: an igraph vector containing the vertex types, or an attribute
          name. Anything that evalulates to C{False} corresponds to vertices of the
          first kind, everything else to the second kind.
        @param probe1: this argument can be used to specify the order of the
          projections in the resulting list. If given and non-negative, then it is
          considered as a vertex ID; the projection containing the vertex will be
          the first one in the result.
        @return: a tuple containing the two projected one-mode graphs.
        """
        return super(Graph, self).bipartite_projection(types, *args, **kwds)

    def bipartite_projection_size(self, types="type", *args, **kwds):
        """bipartite_projection(types="type")

        Calculates the number of vertices and edges in the bipartite projections
        of this graph according to the specified vertex types. This is useful if
        you have a bipartite graph and you want to estimate the amount of memory
        you would need to calculate the projections themselves.
        
        @param types: an igraph vector containing the vertex types, or an attribute
          name. Anything that evalulates to C{False} corresponds to vertices of the
          first kind, everything else to the second kind.
        @return: a 4-tuple containing the number of vertices and edges in the
          first projection, followed by the number of vertices and edges in the
          second projection.
        """
        return super(Graph, self).bipartite_projection_size(types, *args, **kwds)

    def get_incidence(self, types="type", *args, **kwds):
        """get_incidence(self, types="type")

        Returns the incidence matrix of a bipartite graph. The incidence matrix
        is an M{n} times M{m} matrix, where M{n} and M{m} are the number of vertices
        in the two vertex classes.

        @param types: an igraph vector containing the vertex types, or an attribute
          name. Anything that evalulates to C{False} corresponds to vertices of the
          first kind, everything else to the second kind.
        @return: the incidence matrix and two lists in a triplet. The first list
          define the mapping between row indices of the matrix and the
          original vertex IDs. The second list is the same for the column indices.
        """
        return super(Graph, self).get_incidence(types, *args, **kwds)

    ###################
    # Custom operators

    def __iadd__(self, other):
        """In-place addition (disjoint union).

        @see: L{__add__}
        """
        if isinstance(other, int):
            return self.add_vertices(other)
        elif isinstance(other, tuple) and len(other) == 2:
            return self.add_edges([other])
        elif isinstance(other, list):
            if len(other)>0:
                if isinstance(other[0], tuple):
                    return self.add_edges(other)
            else:
                return self

        return NotImplemented


    def __add__(self, other):
        """Copies the graph and extends the copy depending on the type of
        the other object given.

        @param other: if it is an integer, the copy is extended by the given
          number of vertices. If it is a tuple with two elements, the copy
          is extended by a single edge. If it is a list of tuples, the copy
          is extended by multiple edges. If it is a L{Graph}, a disjoint
          union is performed.
        """
        if isinstance(other, int):
            g = self.copy()
            g.add_vertices(other)
        elif isinstance(other, tuple) and len(other) == 2:
            g = self.copy()
            g.add_edges([other])
        elif isinstance(other, list):
            if len(other)>0:
                if isinstance(other[0], tuple):
                    g = self.copy()
                    g.add_edges(other)
                elif isinstance(other[0], Graph):
                    return self.disjoint_union(other)
                else:
                    return NotImplemented
            else:
                return self.copy()

        elif isinstance(other, Graph):
            return self.disjoint_union(other)
        else:
            return NotImplemented

        return g


    def __isub__(self, other):
        """In-place subtraction (difference).

        @see: L{__sub__}"""
        if isinstance(other, int):
            return self.delete_vertices(other)
        elif isinstance(other, tuple) and len(other) == 2:
            return self.delete_edges([other])
        elif isinstance(other, list):
            if len(other)>0:
                if isinstance(other[0], tuple):
                    return self.delete_edges(other)
                elif isinstance(other[0], int):
                    return self.delete_vertices(other)
            else:
                return self
        elif isinstance(other, core.Vertex):
            return self.delete_vertices(other)
        elif isinstance(other, core.VertexSeq):
            return self.delete_vertices(other)
        elif isinstance(other, core.Edge):
            return self.delete_edges(other)
        elif isinstance(other, core.EdgeSeq):
            return self.delete_edges(other)
        return NotImplemented


    def __sub__(self, other):
        """Removes the given object(s) from the graph

        @param other: if it is an integer, removes the vertex with the given
          ID from the graph (note that the remaining vertices will get
          re-indexed!). If it is a tuple, removes the given edge. If it is
          a graph, takes the difference of the two graphs. Accepts
          lists of integers or lists of tuples as well, but they can't be
          mixed! Also accepts L{Edge} and L{EdgeSeq} objects.
        """
        if isinstance(other, int):
            return self.copy().delete_vertices(other)
        elif isinstance(other, tuple) and len(other) == 2:
            return self.copy().delete_edges(other)
        elif isinstance(other, list):
            if len(other)>0:
                if isinstance(other[0], tuple):
                    return self.copy().delete_edges(other)
                elif isinstance(other[0], int):
                    return self.copy().delete_vertices(other)
            else:
                return self.copy()
        elif isinstance(other, core.Vertex):
            return self.copy().delete_vertices(other)
        elif isinstance(other, core.VertexSeq):
            return self.copy().delete_vertices(other)
        elif isinstance(other, core.Edge):
            return self.copy().delete_edges(other)
        elif isinstance(other, core.EdgeSeq):
            return self.copy().delete_edges(other)
        elif isinstance(other, Graph):
            return self.difference(other)

        return NotImplemented

    def __mul__(self, other):
        """Copies exact replicas of the original graph an arbitrary number of times.

        @param other: if it is an integer, multiplies the graph by creating the
          given number of identical copies and taking the disjoint union of
          them.
        """
        if isinstance(other, int):
            if other == 0:
                return Graph()
            elif other == 1:
                return self
            elif other > 1:
                # TODO: should make it more efficient - powers of 2?
                return self.disjoint_union([self]*(other-1))
            else:
                return NotImplemented

        return NotImplemented
    
    def __coerce__(self, other):
        """Coercion rules.

        This method is needed to allow the graph to react to additions
        with lists, tuples, integers, vertices, edges and so on.
        """
        if type(other) in [int, tuple, list]:
            return self, other
        if isinstance(other, core.Vertex):
            return self, other
        if isinstance(other, core.VertexSeq):
            return self, other
        if isinstance(other, core.Edge):
            return self, other
        if isinstance(other, core.EdgeSeq):
            return self, other


    def __reduce__(self):
        """Support for pickling."""
        import warnings
        constructor = self.__class__
        graph_attr_names = self.attributes()
        vertex_attr_names = self.vs.attribute_names()
        edge_attr_names = self.es.attribute_names()
        gattrs, vattrs, eattrs = {}, {}, {}
        for a in graph_attr_names: gattrs[a] = self[a]
        for a in vertex_attr_names: vattrs[a] = self.vs[a]
        for a in edge_attr_names: eattrs[a] = self.es[a]
        parameters = (self.vcount(), self.get_edgelist(), self.is_directed(), \
            gattrs, vattrs, eattrs)
        return (constructor, parameters, {})


    def __plot__(self, context, bbox, palette, *args, **kwds):
        """Plots the graph to the given Cairo context in the given bounding box
       
        The visual style of vertices and edges can be modified at three
        places in the following order of precedence (lower indices override
        higher indices):

          1. Keyword arguments of this function (or of L{plot()} which is
             passed intact to C{Graph.__plot__()}.

          2. Vertex or edge attributes, specified later in the list of
             keyword arguments.

          3. igraph-wide plotting defaults (see
             L{igraph.config.Configuration})

          4. Built-in defaults.

        E.g., if the C{vertex_size} keyword attribute is not present,
        but there exists a vertex attribute named C{size}, the sizes of
        the vertices will be specified by that attribute.

        Besides the usual self-explanatory plotting parameters (C{context},
        C{bbox}, C{palette}), it accepts the following keyword arguments:
        
          - C{layout}: the layout to be used. If not an instance of
            L{Layout}, it will be passed to L{Graph.layout} to calculate
            the layout. Note that if you want a deterministic layout that
            does not change with every plot, you must either use a deterministic
            layout function (like L{Graph.layout_circle}) or calculate the
            layout in advance and pass a L{Layout} object here.
            
          - C{margin}: the top, right, bottom, left margins as a 4-tuple.
            If it has less than 4 elements or is a single float, the elements
            will be re-used until the length is at least 4.

          - C{vertex_size}: size of the vertices. The corresponding vertex
            attribute is called C{size}. The default is 10. Vertex sizes
            are measured in the unit of the Cairo context on which igraph
            is drawing.

          - C{vertex_color}: color of the vertices. The corresponding vertex
            attribute is C{color}, the default is red.  Colors can be
            specified either by common X11 color names (see the source
            code of L{igraph.colors} for a list of known colors), by
            3-tuples of floats (ranging between 0 and 1 for the R, G and
            B components), by CSS-style string specifications (C{#rrggbb})
            or by integer color indices of the specified palette.

          - C{vertex_shape}: shape of the vertices. Alternatively it can
            be specified by the C{shape} vertex attribute. Possibilities
            are: C{square}, {circle}, {triangle}, {triangle-down} or
            C{hidden}. See the source code of L{igraph.drawing} for a
            list of alternative shape names that are also accepted and
            mapped to these.

          - C{vertex_label}: labels drawn next to the vertices.
            The corresponding vertex attribute is C{label}.

          - C{vertex_label_dist}: distance of the midpoint of the vertex
            label from the center of the corresponding vertex.
            The corresponding vertex attribute is C{label_dist}.

          - C{vertex_label_color}: color of the label. Corresponding
            vertex attribute: C{label_color}. See C{vertex_color} for
            color specification syntax.

          - C{vertex_label_size}: font size of the label, specified
            in the unit of the Cairo context on which we are drawing.
            Corresponding vertex attribute: C{label_size}.

          - C{vertex_label_angle}: the direction of the line connecting
            the midpoint of the vertex with the midpoint of the label.
            This can be used to position the labels relative to the
            vertices themselves in conjunction with C{vertex_label_dist}.
            Corresponding vertex attribute: C{label_angle}. The
            default is C{-math.pi/4}.

          - C{edge_color}: color of the edges. The corresponding edge
            attribute is C{color}, the default is red. See C{vertex_color}
            for color specification syntax.

          - C{edge_width}: width of the edges in the default unit of the
            Cairo context on which we are drawing. The corresponding
            edge attribute is C{width}, the default is 1.

          - C{edge_arrow_size}: arrow size of the edges. The
            corresponding edge attribute is C{arrow_size}, the default
            is 1.

          - C{edge_arrow_width}: width of the arrowhead on the edge. The
            corresponding edge attribute is C{arrow_width}, the default
            is 1.
        """
        import colors
        import cairo

        directed = self.is_directed()
        margin = kwds.get("margin", [0., 0., 0., 0.])
        try:
            margin = list(margin)
        except TypeError:
            margin = [margin]
        while len(margin)<4: margin.extend(margin)
        margin = tuple(map(float, margin[:4]))

        vertex_colors = drawing.collect_attributes(self.vcount(), "vertex_color", \
            "color", kwds, self.vs, config, "red", palette.get)
        vertex_sizes = drawing.collect_attributes(self.vcount(), "vertex_size", \
            "size", kwds, self.vs, config, 10, float)
        vertex_shapes = [drawing.known_shapes.get(x, drawing.NullDrawer) \
            for x in drawing.collect_attributes(self.vcount(), "vertex_shape", \
            "shape", kwds, self.vs, config, "circle")]

        max_vertex_size = max(vertex_sizes)

        layout = kwds.get("layout", None)
        if isinstance(layout, Layout):
            layout = Layout(layout.coords)
        elif isinstance(layout, str) or layout is None:
            layout = self.layout(layout)
        else:
            layout = Layout(layout)

        sl, st, sr, sb = layout.bounding_box()
        sw, sh = sr-sl, sb-st
        if sw == 0 and sh == 0: sw, sh = 1, 1
        if sw == 0: sw = sh
        if sh == 0: sh = sw
        rx, ry = float(bbox.width-max_vertex_size-margin[1]-margin[3])/sw, \
          float(bbox.height-max_vertex_size-margin[0]-margin[2])/sh
        layout.scale(rx, ry)
        layout.translate(-sl*rx+max_vertex_size/2.+margin[1]+bbox.coords[0], \
          -st*ry+max_vertex_size/2.+margin[0]+bbox.coords[1])
        context.set_line_width(1)

        edge_colors = drawing.collect_attributes(self.ecount(), "edge_color", \
            "color", kwds, self.es, config, "black", palette.get)
        edge_widths = drawing.collect_attributes(self.ecount(), "edge_width", \
            "width", kwds, self.es, config, 1, float)
        edge_arrow_sizes = drawing.collect_attributes(self.ecount(), \
            "edge_arrow_size", "arrow_size", kwds, self.es, config, 1, float)
        edge_arrow_widths = drawing.collect_attributes(self.ecount(), \
            "edge_arrow_width", "arrow_width", kwds, self.es, config, 1, float)

        # Draw the edges
        for idx, e in enumerate(self.es):
            context.set_source_rgb(*edge_colors[idx])
            context.set_line_width(edge_widths[idx])

            src, tgt = e.tuple
            if src == tgt:
                # Loop edge
                r = vertex_sizes[src]*2
                cx, cy = layout[src][0]+math.cos(math.pi/4)*r/2, \
                  layout[src][1]-math.sin(math.pi/4)*r/2
                context.arc(cx, cy, r/2., 0, math.pi*2)
            else:
                # Determine where the edge intersects the circumference of the
                # vertex shape. TODO: theoretically this need not to be done
                # if there are no arrowheads on the edge, but maybe it's not worth
                # testing for
                p1 = vertex_shapes[src].intersection_point( \
                    layout[src][0], layout[src][1], layout[tgt][0], layout[tgt][1],
                    vertex_sizes[src])
                p2 = vertex_shapes[tgt].intersection_point( \
                    layout[tgt][0], layout[tgt][1], layout[src][0], layout[src][1],
                    vertex_sizes[tgt])
                context.move_to(*p1)
                context.line_to(*p2)
            context.stroke()

            if directed and src != tgt:
                # Draw an arrowhead
                angle = math.atan2(p2[1]-p1[1], p2[0]-p1[0])
                arrow_size = 15.*edge_arrow_sizes[idx]
                arrow_width = 10./edge_arrow_widths[idx]
                a1 = (p2[0]-arrow_size*math.cos(angle-math.pi/arrow_width),
                  p2[1]-arrow_size*math.sin(angle-math.pi/arrow_width))
                a2 = (p2[0]-arrow_size*math.cos(angle+math.pi/arrow_width),
                  p2[1]-arrow_size*math.sin(angle+math.pi/arrow_width))
                context.move_to(*p2)
                context.line_to(*a1)
                context.line_to(*a2)
                context.line_to(*p2)
                context.fill()

        del edge_colors
        del edge_widths

        # Draw the vertices
        context.set_line_width(1)
        for idx, v in enumerate(self.vs):
            vertex_shapes[idx].draw_path(context, layout[idx][0], layout[idx][1], vertex_sizes[idx])
            context.set_source_rgb(*vertex_colors[idx])
            context.fill_preserve()
            context.set_source_rgb(0., 0., 0.)
            context.stroke()
        del vertex_colors
        del vertex_shapes

        # Draw the vertex labels
        if not kwds.has_key("vertex_label") and "label" not in self.vs.attribute_names():
            vertex_labels = map(str, xrange(self.vcount()))
        elif kwds.has_key("vertex_label") and kwds["vertex_label"] is None:
            vertex_labels = [""] * self.vcount()
        else:
            vertex_labels = drawing.collect_attributes(self.vcount(), "vertex_label", \
                "label", kwds, self.vs, config, None)
        vertex_dists = drawing.collect_attributes(self.vcount(), "vertex_label_dist", \
            "label_dist", kwds, self.vs, config, 1, float)
        vertex_degrees = drawing.collect_attributes(self.vcount(), \
            "vertex_label_angle", "label_angle", kwds, self.vs, config, \
            -math.pi/4, float)
        vertex_label_colors = drawing.collect_attributes(self.vcount(), \
            "vertex_label_color", "label_color", kwds, self.vs, config, \
            "black", palette.get)
        vertex_label_sizes = drawing.collect_attributes(self.vcount(), \
            "vertex_label_size", "label_size", kwds, self.vs, \
            config, 14, float)

        context.select_font_face("sans-serif", cairo.FONT_SLANT_NORMAL, \
            cairo.FONT_WEIGHT_BOLD)
        
        for idx, v in enumerate(self.vs):
            xb, yb, w, h = context.text_extents(vertex_labels[idx])[:4]
            cx, cy = layout[idx]
            cx += math.cos(vertex_degrees[idx]) * vertex_dists[idx] * vertex_sizes[idx]
            cy += math.sin(vertex_degrees[idx]) * vertex_dists[idx] * vertex_sizes[idx]
            cx -= w/2. + xb
            cy -= h/2. + yb
            context.move_to(cx, cy)
            context.set_font_size(vertex_label_sizes[idx])
            context.set_source_rgb(*vertex_label_colors[idx])
            context.text_path(vertex_labels[idx])
            context.fill()

    def summary(self, verbosity=0):
        """Returns basic statistics about the graph in a string
        
        @param verbosity: the amount of statistics to be returned. 0 returns
          the usual statistics (node, edge count, directedness, number of
          strong components, density, reciprocity, average path length,
          diameter). 1 also returns the detailed degree distributions."""
        output=[]
        output.append("%d nodes, %d edges, %sdirected" % \
            (self.vcount(), self.ecount(), ["un", ""][self.is_directed()]))
        output.append("")
        output.append("Number of components: %d" % len(self.clusters()))
        output.append("Diameter: %d" % self.diameter(unconn=True))
        output.append("Density: %.4f" % self.density())
        # output.append("Transitivity: %.4f" % self.transitivity())
        if self.is_directed():
            output.append("Reciprocity: %.4f" % self.reciprocity())
        output.append("Average path length: %.4f" % self.average_path_length())

        if verbosity>=1:
            maxdegree=self.maxdegree()
            binwidth=max(1, maxdegree/20)
            output.append("")
            output.append("Degree distribution:")
            output.append(str(self.degree_distribution(binwidth)))

            if self.is_directed():
                output.append("")
                output.append("Degree distribution (only in-degrees):")
                output.append(str(self.degree_distribution(binwidth, type=IN)))
                output.append("")
                output.append("Degree distribution (only out-degrees):")
                output.append(str(self.degree_distribution(binwidth, type=OUT)))

        return "\n".join(output)

    _format_mapping = {
          "ncol":       ("Read_Ncol", "write_ncol"),
          "lgl":        ("Read_Lgl", "write_lgl"),
          "graphmlz":   ("Read_GraphMLz", "write_graphmlz"),
          "graphml":    ("Read_GraphML", "write_graphml"),
          "gml":        ("Read_GML", "write_gml"),
          "dot":		(None, "write_dot"),
          "graphviz":	(None, "write_dot"),
          "net":        ("Read_Pajek", "write_pajek"),
          "pajek":      ("Read_Pajek", "write_pajek"),
          "dimacs":     ("Read_DIMACS", "write_dimacs"),
          "adjacency":  ("Read_Adjacency", "write_adjacency"),
          "adj":        ("Read_Adjacency", "write_adjacency"),
          "edgelist":   ("Read_Edgelist", "write_edgelist"),
          "edge":       ("Read_Edgelist", "write_edgelist"),
          "edges":      ("Read_Edgelist", "write_edgelist"),
          "pickle":     ("Read_Pickle", "write_pickle"),
          "svg":        (None, "write_svg")
    }

    _layout_mapping = {
        "circle": "layout_circle",
        "circular": "layout_circle",
        "drl": "layout_drl",
        "fr": "layout_fruchterman_reingold",
        "fruchterman_reingold": "layout_fruchterman_reingold",
        "fr3d": "layout_fruchterman_reingold_3d",
        "fr_3d": "layout_fruchterman_reingold_3d",
        "fruchterman_reingold_3d": "layout_fruchterman_reingold_3d",
        "gfr": "layout_grid_fruchterman_reingold",
        "graphopt": "layout_graphopt",
        "grid_fr": "layout_grid_fruchterman_reingold",
        "grid_fruchterman_reingold": "layout_grid_fruchterman_reingold",
        "kk": "layout_kamada_kawai",
        "kamada_kawai": "layout_kamada_kawai",
        "kk3d": "layout_kamada_kawai_3d",
        "kk_3d": "layout_kamada_kawai_3d",
        "kamada_kawai_3d": "layout_kamada_kawai_3d",
        "lgl": "layout_lgl",
        "large": "layout_lgl",
        "large_graph": "layout_lgl",
        "random": "layout_random",
        "random_3d": "layout_random_3d",
        "rt": "layout_reingold_tilford",
        "tree": "layout_reingold_tilford",
        "reingold_tilford": "layout_reingold_tilford",
        "rt_circular": "layout_reingold_tilford_circular",
        "reingold_tilford_circular": "layout_reingold_tilford_circular",
        "sphere": "layout_sphere",
        "spherical": "layout_sphere",
        "circle_3d": "layout_sphere",
        "circular_3d": "layout_sphere",
    }

    # After adjusting something here, don't forget to update the docstring
    # of Graph.layout if necessary!

##############################################################

class VertexSeq(core.VertexSeq):
    """Class representing a sequence of vertices in the graph.
    
    This class is most easily accessed by the C{vs} field of the
    L{Graph} object, which returns an ordered sequence of all vertices in
    the graph. The vertex sequence can be refined by invoking the
    L{VertexSeq.select()} method. L{VertexSeq.select()} can also be
    accessed by simply calling the L{VertexSeq} object.

    An alternative way to create a vertex sequence referring to a given
    graph is to use the constructor directly:
    
      >>> g = Graph.Full(3)
      >>> vs = VertexSeq(g)
      >>> restricted_vs = VertexSeq(g, [0, 1])

    The individual vertices can be accessed by indexing the vertex sequence
    object. It can be used as an iterable as well, or even in a list
    comprehension:
    
      >>> g=Graph.Full(3)
      >>> for v in g.vs:
      ...   v["value"] = v.index ** 2
      ...
      >>> [v["value"] ** 0.5 for v in g.vs]
      [0.0, 1.0, 2.0]
      
    The vertex set can also be used as a dictionary where the keys are the
    attribute names. The values corresponding to the keys are the values
    of the given attribute for every vertex selected by the sequence.
    
      >>> g=Graph.Full(3)
      >>> for idx, v in enumerate(g.vs):
      ...   v["weight"] = idx*(idx+1)
      ...
      >>> g.vs["weight"]
      [0, 2, 6]
      >>> g.vs.select(1,2)["weight"] = [10, 20]
      >>> g.vs["weight"]
      [0, 10, 20]

    Some methods of the vertex sequences are simply proxy methods to the
    corresponding methods in the L{Graph} object. One such example is
    L{VertexSeq.degree()}:

      >>> g=Graph.Tree(7, 2)
      >>> g.vs.degree()
      [2, 3, 3, 1, 1, 1, 1]
      >>> g.vs.degree() == g.degree()
      True
    """

    def select(self, *args, **kwds):
        """Selects a subset of the vertex sequence based on some criteria
        
        The selection criteria can be specified by the positional and the keyword
        arguments. Positional arguments are always processed before keyword
        arguments.
        
          - If the first positional argument is C{None}, an empty sequence is
            returned.
            
          - If the first positional argument is a callable object, the object
            will be called for every vertex in the sequence. If it returns
            C{True}, the vertex will be included, otherwise it will
            be excluded.
            
          - If the first positional argument is an iterable, it must return
            integers and they will be considered as indices of the current
            vertex set (NOT the whole vertex set of the graph -- the
            difference matters when one filters a vertex set that has
            already been filtered by a previous invocation of
            L{VertexSeq.select()}. In this case, the indices do not refer
            directly to the vertices of the graph but to the elements of
            the filtered vertex sequence.
            
          - If the first positional argument is an integer, all remaining
            arguments are expected to be integers. They are considered as
            indices of the current vertex set again.

        Keyword arguments can be used to filter the vertices based on their
        attributes. The name of the keyword specifies the name of the attribute
        and the filtering operator, they should be concatenated by an
        underscore (C{_}) character. Attribute names can also contain
        underscores, but operator names don't, so the operator is always the
        largest trailing substring of the keyword name that does not contain
        an underscore. Possible operators are:

          - C{eq}: equal to

          - C{ne}: not equal to

          - C{lt}: less than
          
          - C{gt}: greater than

          - C{le}: less than or equal to

          - C{ge}: greater than or equal to

          - C{in}: checks if the value of an attribute is in a given list

          - C{notin}: checks if the value of an attribute is not in a given
            list

        For instance, if you want to filter vertices with a numeric C{age}
        property larger than 200, you have to write:

          >>> g.vs.select(age_gt=200)

        Similarly, to filter vertices whose C{type} is in a list of predefined
        types:

          >>> list_of_types = ["HR", "Finance", "Management"]
          >>> g.vs.select(type_in=list_of_types)

        If the operator is omitted, it defaults to C{eq}. For instance, the
        following selector selects vertices whose C{cluster} property equals
        to 2:

          >>> g.vs.select(cluster=2)

        In the case of an unknown operator, it is assumed that the
        recognized operator is part of the attribute name and the actual
        operator is C{eq}.

        Attribute names inferred from keyword arguments are treated specially
        if they start with an underscore (C{_}). These are not real attributes
        but refer to specific properties of the vertices, e.g., its degree.
        The rule is as follows: if an attribute name starts with an underscore,
        the rest of the name is interpreted as a method of the L{Graph} object.
        This method is called with the vertex sequence as its first argument
        (all others left at default values) and vertices are filtered
        according to the value returned by the method. For instance, if you
        want to exclude isolated vertices:

          >>> non_isolated = g.vs.select(_degree_gt=0)

        For properties that take a long time to be computed (e.g., betweenness
        centrality for large graphs), it is advised to calculate the values
        in advance and store it in a graph attribute. The same applies when
        you are selecting based on the same property more than once in the
        same C{select()} call to avoid calculating it twice unnecessarily.
        For instance, the following would calculate betweenness centralities
        twice:

          >>> g.vs.select(_betweenness_gt=10, _betweenness_lt=30)

        It is advised to use this instead:

          >>> g.vs["bs"] = g.betwenness()
          >>> g.vs.select(bs_gt=10, bs_lt=30)

        @return: the new, filtered vertex sequence"""
        vs = core.VertexSeq.select(self, *args)

        operators = {
            "lt": operator.lt, \
            "gt": operator.gt, \
            "le": operator.le, \
            "ge": operator.ge, \
            "eq": operator.eq, \
            "ne": operator.ne, \
            "in": lambda a, b: a in b, \
            "notin": lambda a, b: a not in b }
        for keyword, value in kwds.iteritems():
            if "_" not in keyword or keyword.rindex("_") == 0: keyword = keyword+"_eq"
            pos = keyword.rindex("_")
            attr, op = keyword[0:pos], keyword[pos+1:]
            try:
                func = operators[op]
            except KeyError:
                # No such operator, assume that it's part of the attribute name
                attr = "%s_%s" % (attr,op)
                func = operators["eq"]

            if attr[0] == '_':
                # Method call, not an attribute
                values = getattr(vs.graph, attr[1:])(vs) 
            else:
                values = vs[attr]
            filtered_idxs=[i for i, v in enumerate(values) if func(v, value)]
            vs = vs.select(filtered_idxs)

        return vs

    def __call__(self, *args, **kwds):
        """Shorthand notation to select()

        This method simply passes all its arguments to L{VertexSeq.select()}.
        """
        return self.select(*args, **kwds)

##############################################################

class EdgeSeq(core.EdgeSeq):
    """Class representing a sequence of edges in the graph.
    
    This class is most easily accessed by the C{es} field of the
    L{Graph} object, which returns an ordered sequence of all edges in
    the graph. The edge sequence can be refined by invoking the
    L{EdgeSeq.select()} method. L{EdgeSeq.select()} can also be
    accessed by simply calling the L{EdgeSeq} object.

    An alternative way to create an edge sequence referring to a given
    graph is to use the constructor directly:
    
      >>> g = Graph.Full(3)
      >>> es = EdgeSeq(g)
      >>> restricted_es = EdgeSeq(g, [0, 1])

    The individual edges can be accessed by indexing the edge sequence
    object. It can be used as an iterable as well, or even in a list
    comprehension:
    
      >>> g=Graph.Full(3)
      >>> for e in g.es:
      ...   print e.tuple
      ...
      (0, 1)
      (0, 2)
      (1, 2)
      >>> [max(e.tuple) for e in g.es]
      [1, 2, 2]
      
    The edge sequence can also be used as a dictionary where the keys are the
    attribute names. The values corresponding to the keys are the values
    of the given attribute of every edge in the graph:
    
      >>> g=Graph.Full(3)
      >>> for idx, e in enumerate(g.es):
      ...   e["weight"] = idx*(idx+1)
      ...
      >>> g.es["weight"]
      [0, 2, 6]
      >>> g.es["weight"] = range(3)
      >>> g.es["weight"]
      [0, 1, 2]

    Some methods of the edge sequences are simply proxy methods to the
    corresponding methods in the L{Graph} object. One such example is
    L{EdgeSeq.is_multiple()}:

      >>> g=Graph(3, [(0,1), (1,0), (1,2)])
      >>> g.es.is_multiple()
      [False, True, False]
      >>> g.es.is_multiple() == g.is_multiple()
      True
    """

    def select(self, *args, **kwds):
        """Selects a subset of the edge sequence based on some criteria
        
        The selection criteria can be specified by the positional and the keyword
        arguments. Positional arguments are always processed before keyword
        arguments.
        
          - If the first positional argument is C{None}, an empty sequence is
            returned.
            
          - If the first positional argument is a callable object, the object
            will be called for every edge in the sequence. If it returns
            C{True}, the edge will be included, otherwise it will
            be excluded.
            
          - If the first positional argument is an iterable, it must return
            integers and they will be considered as indices of the current
            edge set (NOT the whole edge set of the graph -- the
            difference matters when one filters an edge set that has
            already been filtered by a previous invocation of
            L{EdgeSeq.select()}. In this case, the indices do not refer
            directly to the edges of the graph but to the elements of
            the filtered edge sequence.
            
          - If the first positional argument is an integer, all remaining
            arguments are expected to be integers. They are considered as
            indices of the current edge set again.

        Keyword arguments can be used to filter the edges based on their
        attributes. The name of the keyword specifies the name of the attribute
        and the filtering operator, they should be concatenated by an
        underscore (C{_}) character. Attribute names can also contain
        underscores, but operator names don't, so the operator is always the
        largest trailing substring of the keyword name that does not contain
        an underscore. Possible operators are:

          - C{eq}: equal to

          - C{ne}: not equal to

          - C{lt}: less than
          
          - C{gt}: greater than

          - C{le}: less than or equal to

          - C{ge}: greater than or equal to

          - C{in}: checks if the value of an attribute is in a given list

          - C{notin}: checks if the value of an attribute is not in a given
            list

        For instance, if you want to filter edges with a numeric C{weight}
        property larger than 50, you have to write:

          >>> g.es.select(weight_gt=50)

        Similarly, to filter edges whose C{type} is in a list of predefined
        types:

          >>> list_of_types = ["inhibitory", "excitatory"]
          >>> g.es.select(type_in=list_of_types)

        If the operator is omitted, it defaults to C{eq}. For instance, the
        following selector selects edges whose C{type} property is
        C{intracluster}:

          >>> g.es.select(type="intracluster")

        In the case of an unknown operator, it is assumed that the
        recognized operator is part of the attribute name and the actual
        operator is C{eq}.

        Attribute names inferred from keyword arguments are treated specially
        if they start with an underscore (C{_}). These are not real attributes
        but refer to specific properties of the edges, e.g., their centrality.
        The rule is as follows: if an attribute name starts with an underscore,
        the rest of the name is interpreted as a method of the L{Graph} object.
        This method is called with the edge sequence as its first argument
        (all others left at default values) and edges are filtered
        according to the value returned by the method. For instance, if you
        want to exclude edges with a betweenness centrality less than 2:

          >>> excl = g.es.select(_edge_betweenness_ge = 2)

        For properties that take a long time to be computed (e.g., betweenness
        centrality for large graphs), it is advised to calculate the values
        in advance and store it in a graph attribute. The same applies when
        you are selecting based on the same property more than once in the
        same C{select()} call to avoid calculating it twice unnecessarily.
        For instance, the following would calculate betweenness centralities
        twice:

          >>> g.es.select(_edge_betweenness_gt=10, _edge_betweenness_lt=30)

        It is advised to use this instead:

          >>> g.es["bs"] = g.edge_betwenness()
          >>> g.es.select(bs_gt=10, bs_lt=30)

        @return: the new, filtered edge sequence"""
        es = core.EdgeSeq.select(self, *args)

        operators = {
            "lt": operator.lt, \
            "gt": operator.gt, \
            "le": operator.le, \
            "ge": operator.ge, \
            "eq": operator.eq, \
            "ne": operator.ne, \
            "in": lambda a, b: a in b, \
            "notin": lambda a, b: a not in b }
        for keyword, value in kwds.iteritems():
            if "_" not in keyword or keyword.rindex("_") == 0: keyword = keyword+"_eq"
            pos = keyword.rindex("_")
            attr, op = keyword[0:pos], keyword[pos+1:]
            try:
                func = operators[op]
            except KeyError:
                # No such operator, assume that it's part of the attribute name
                attr = "%s_%s" % (attr,op)
                func = operators["eq"]

            if attr[0] == '_':
                # Method call, not an attribute
                values = getattr(es.graph, attr[1:])(es) 
            else:
                values = es[attr]
            filtered_idxs=[i for i, v in enumerate(values) if func(v, value)]
            es = es.select(filtered_idxs)

        return es


    def __call__(self, *args, **kwds):
        """Shorthand notation to select()

        This method simply passes all its arguments to L{EdgeSeq.select()}.
        """
        return self.select(*args, **kwds)

##############################################################
# Additional methods of VertexSeq and EdgeSeq that call Graph methods

def _graphmethod(func=None, name=None):
    """Auxiliary decorator
    
    This decorator allows some methods of L{VertexSeq} and L{EdgeSeq} to
    call their respective counterparts in L{Graph} to avoid code duplication.

    @param func: the function being decorated. This function will be
      called on the results of the original L{Graph} method.
      If C{None}, defaults to the identity function.
    @param name: the name of the corresponding method in L{Graph}. If
      C{None}, it defaults to the name of the decorated function.
    @return: the decorated function
    """
    if name is None: name = func.__name__
    method = getattr(Graph, name)

    if callable(func):
        def decorated(*args, **kwds):
            self = args[0].graph
            return func(args[0], method(self, *args, **kwds))
    else:
        def decorated(*args, **kwds):
            self = args[0].graph
            return method(self, *args, **kwds)

    decorated.__name__ = name
    decorated.__doc__ = """Proxy method to L{Graph.%(name)s()}

This method calls the C{%(name)s()} method of the L{Graph} class
restricted to this sequence, and returns the result.

@see: Graph.%(name)s() for details.
""" % { "name": name }

    return decorated

def _add_proxy_methods():
    decorated_methods = {}
    decorated_methods[VertexSeq] = \
        ["degree", "betweenness", "bibcoupling", "closeness", "cocitation",
        "constraint", "eccentricity", "get_shortest_paths", "maxdegree",
        "pagerank", "shortest_paths", "similarity_dice", "similarity_jaccard",
        "subgraph", "indegree", "outdegree", "isoclass"]
    decorated_methods[EdgeSeq] = \
        ["count_multiple", "is_loop", "is_multiple"]

    for klass, methods in decorated_methods.iteritems():
        for method in methods:
            setattr(klass, method, _graphmethod(None, method))

    setattr(EdgeSeq, "edge_betweenness", _graphmethod( \
      lambda self, result: [result[i] for i in self.indices], "edge_betweenness"))

_add_proxy_methods()

##############################################################

def read(filename, *args, **kwds):
    """Loads a graph from the given filename.

    This is just a convenience function, calls L{Graph.Read} directly.
    All arguments are passed unchanged to L{Graph.Read}
    
    @param filename: the name of the file to be loaded
    """
    return Graph.Read(filename, *args, **kwds)
load=read

def write(graph, filename, *args, **kwds):
    """Saves a graph to the given file.

    This is just a convenience function, calls L{Graph.write} directly.
    All arguments are passed unchanged to L{Graph.write}

    @param graph: the graph to be saved
    @param filename: the name of the file to be written
    """
    return graph.write(filename, *args, **kwds)
save=write

def summary(o, f=sys.stdout):
    """Prints a summary of object o to a given stream

    @param o: the object about which a human-readable summary is requested.
    @param f: the stream to be used
    """
    if hasattr(o, "summary"):
        print >>f, o.summary()
    else:
        print >>f, str(o)

config = configuration.init()

"""
IGraph library.

@undocumented: test
"""

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

from core import *
from core import __version__, __build_date__
from clustering import *
from layout import *

import os
import math
import gzip
import sys
from tempfile import mkstemp
from warnings import warn

def summary(o, f=sys.stdout):
    """Prints a summary of object o to a given stream

    @param o: the object about which a human-readable summary is requested.
    @param f: the stream to be used
    """
    if hasattr(o, "summary"):
        print >>f, o.summary()
    else:
        print >>f, str(o)


class Graph(core.GraphBase):
    """Generic graph.
    
    This class is built on top of L{GraphBase}, so the order of the
    methods in the Epydoc documentation is a little bit obscure:
    inherited methods come after the ones implemented directly in the
    subclass.
    """

    # Some useful aliases
    omega = core.GraphBase.clique_number
    alpha = core.GraphBase.independence_number
    shell_index = core.GraphBase.coreness
    cut_vertices = core.GraphBase.articulation_points

    def indegree(self, *args, **kwds):
        """Returns the in-degrees in a list.
        
        See L{degree} for possible arguments.
        """
        kwds['degree']=_igraph.IN
        return self.degree(*args, **kwds)

    def outdegree(self, *args, **kwds):
        """Returns the out-degrees in a list.
        
        See L{degree} for possible arguments.
        """
        kwds['degree']=_igraph.OUT
        return self.degree(*args, **kwds)

    def eccentricity(self, nodes=None):
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
        
        if nodes is None:
            result = [1.0/x for x in distance_maxs]
        else:
            result = [1.0/distance_maxs[idx] for idx in nodes]

        return result

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

    def modularity(self, membership):
        """modularity(membership)

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

        @param membership: a membership list or a L{VertexClustering} object
        @return: the modularity score
        
        @newfield ref: Reference
        @ref: MEJ Newman and M Girvan: Finding and evaluating community
          structure in networks. Phys Rev E 69 026113, 2004.
        """
        if isinstance(membership, VertexClustering):
            if membership.graph != self:
                raise ValueError, "clustering object belongs to a different graph"
            return GraphBase.modularity(self, membership.membership)
        else:
            return GraphBase.modularity(self, membership)

    # Various clustering algorithms -- mostly wrappers around GraphBase

    def community_fastgreedy(self):
        """Community structure based on the greedy optimization of modularity.

        This algorithm merges individual nodes into communities in a way that
        greedily maximizes the modularity score of the graph. It can be proven
        that if no merge can increase the current modularity score, the algorithm
        can be stopped since no further increase can be achieved.

        This algorithm is said to run almost in linear time on sparse graphs.

        @return: an appropriate L{VertexClustering} object.

        @newfield ref: Reference
        @ref: A Clauset, MEJ Newman and C Moore: Finding community structure
          in very large networks. Phys Rev E 70, 066111 (2004).
        """
        merges, qs = GraphBase.community_fastgreedy(self, True)
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
        @param return_merges: whether the returned L{VertexClustering} object
          should contain information about the merges performed on the graph.
        @return: an appropriate L{VertexClustering} object.
        @param return_merges: whether 
        
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
        @param return_merges: whether the returned L{VertexClustering} object
          should contain information about the merges performed on the graph.
        @return: an appropriate L{VertexClustering} object.
        @param return_merges: whether 
        
        @newfield ref: Reference
        @ref: MEJ Newman: Finding community structure in networks using the
        eigenvectors of matrices, arXiv:physics/0605087"""
        if clusters is None: clusters=-1
        cl, merges = GraphBase.community_leading_eigenvector(self, clusters, return_merges)
        if merges is None:
            return VertexClustering(self, cl)
        else:
            return VertexDendrogram(self, merges, cl)


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
        @return: a L{VertexClustering} object. Besides the usual methods and members,
          this object will have a member called C{merges} which records information
          used to produce the dendrogram. It is practically a list of tuples where
          each tuple defines two nodes which will be joined in a step. Node IDs
          from 0 to M{n-1} (where M{n} is the number of vertices) correspond to
          the individual vertices, while node IDs up from M{n} correspond to
          merged communities. M{n} means the community created after the first
          merge, M{n+1} means the community created after the second merge and
          so on...
        """
        d = VertexDendrogram(self, GraphBase.community_edge_betweenness(self, directed));
        if clusters is not None: d.cut(clusters)
        return d
    
    def edge_betweenness_clustering(self, clusters = None, steps = None):
        """Newman's edge betweenness clustering.

        Iterative removal of edges with the largest edge betweenness until
        the given number of steps is reached or until the graph is decomposed
        to the given number of clusters. Edge betweennesses are recalculated
        after every run.

        @param clusters: the desired number of clusters.
        @param steps: the number of tests to take.

        @return: an appropriate L{VertexClustering} object.

        @newfield ref: Reference
        @ref: Girvan, M and Newman, MEJ: Community structure in social and
          biological networks. Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002)
        """
        warn("Graph.edge_betweenness_clustering is deprecated and will be removed soon. Use Graph.community_edge_betweenness instead", DeprecationWarning)
        g = self.copy()
        number_of_steps = 0

        directed = g.is_directed()

        while True:
            if clusters is not None:
                cl = g.clusters()
                if max(cl)+1 >= clusters: break

            if steps is not None:
                if number_of_steps > steps: break

            ebs = g.edge_betweenness(directed)

            if len(ebs) == 0: break

            eb_max = max(ebs)
            eb_max_index = ebs.index(eb_max)

            g.delete_edges(eb_max_index, by_index=True)
            number_of_steps += 1

        return VertexClustering(self, g.clusters()) 


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
    Read_GraphMLz = classmethod(Read_GraphMLz)

    def write_svg(self, fname, layout, width = None, height = None, \
                  labels = "label", colors = "color", shapes = "shape", \
                  vertex_size = 10, font_size = 16, *args, **kwds):
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

        print >>f, "<?xml version=\"1.0\" standalone=\"no\"?>"
        print >>f, "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\""
        print >>f, "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
        
        print >>f, "<svg width=\"%d\" height=\"%d\" version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">" % (width, height)
        print >>f, "<!-- Created by igraph -->"
        print >>f
        print >>f, "<defs>"
        print >>f, "  <marker id=\"Triangle\" viewBox=\"0 0 10 10\""
        print >>f, "    refX=\"10\" refY=\"5\" orient=\"auto\""
        print >>f, "    markerUnits=\"strokeWidth\""
        print >>f, "    markerWidth=\"16\" markerHeight=\"8\">"
        print >>f, "    <path d=\"M 0 0 L 10 5 L 0 10 z\"/>"
        print >>f, "  </marker>"
        print >>f, "  <style type=\"text/css\">"
        print >>f, "    <![CDATA["
        print >>f, "#vertices circle { stroke: black; stroke-width: 1 }"
        print >>f, "#vertices rect { stroke: black; stroke-width: 1 }"
        print >>f, "#vertices text { text-anchor: middle; font-size: %s; font-family: sans-serif; font-weight: normal }" % font_size
        if self.is_directed():
            print >>f, "#edges line { stroke: black; stroke-width: 1; marker-end: url(#Triangle) }"
        else:
            print >>f, "#edges line { stroke: black; stroke-width: 1 }"
        print >>f, "    ]]>"
        print >>f, "  </style>"
        print >>f, "</defs>"
        print >>f
        print >>f, "<g transform=\"translate(%.4f,%.4f)\">" % (width/2.0, height/2.0)
        print >>f, "  <g id=\"edges\">"
        print >>f, "  <!-- Edges -->"

        has_edge_opacities = "opacity" in self.edge_attributes()
        for edge in self.es:
            vidxs = edge.tuple
            x1 = layout[vidxs[0]][0]
            y1 = layout[vidxs[0]][1]
            x2 = layout[vidxs[1]][0]
            y2 = layout[vidxs[1]][1]
            angle = math.atan2(y2-y1, x2-x1)
            x2 = x2 - vertex_size*math.cos(angle)
            y2 = y2 - vertex_size*math.sin(angle)
            if has_edge_opacities:
                print >>f, "    <line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\" style=\"stroke-opacity:%.2f\"/>" % (x1, y1, x2, y2, float(edge["opacity"]))
            else:
                print >>f, "    <line x1=\"%.4f\" y1=\"%.4f\" x2=\"%.4f\" y2=\"%.4f\"/>" % (x1, y1, x2, y2)

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
          C{"edgelist"}, C{"edges"} or C{"edge"} (edge list).
        @raises IOError: if the file format can't be identified and
          none was given.
        """
        if format is None: format = klass._identify_format(f)
        try:
            reader = klass._format_mapping[format][0]
        except KeyError, IndexError:
            raise IOError, "unknown file format: %s" % str(format)
        if reader is None:
            raise IOError, "no loader method for file format: %s" % str(format)
        reader = getattr(klass, reader)
        return reader(f, *args, **kwds)
    Read = classmethod(Read)
    Load = Read

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
            return self.delete_edges(other)
        elif isinstance(other, list):
            if len(other)>0:
                if isinstance(other[0], tuple):
                    return self.delete_edges(other)
                elif isinstance(other[0], int):
                    return self.delete_vertices(other)
            else:
                return self

        return NotImplemented


    def __sub__(self, other):
        """Removes the given object(s) from the graph

        @param other: if it is an integer, removes the vertex with the given
          ID from the graph (note that the remaining vertices will get
          re-indexed!). If it is a tuple, removes the given edge. If it is
          a graph, takes the difference of the two graphs. Accepts
          lists of integers or lists of tuples as well, but they can't be
          mixed!
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
        with lists, tuples or integers.
        """
        if type(other) in [int, tuple, list]:
            return self, other


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
        output.append("Number of components: %d" % (len(self.clusters())+1))
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
          "net":        ("Read_Pajek", None),
          "pajek":      ("Read_Pajek", None),
          "dimacs":     ("Read_DIMACS", "write_dimacs"),
          #"adjacency":  ("Read_Adjacency", "write_adjacency"),
          #"adj":        ("Read_Adjacency", "write_adjacency),
          "edgelist":   ("Read_Edgelist", "write_edgelist"),
          "edge":       ("Read_Edgelist", "write_edgelist"),
          "edges":      ("Read_Edgelist", "write_edgelist")
    }


"""
IGraph library.
Copyright (C) 2006  Gabor Csardi <csardi@rmki.kfki.hu>
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

from igraph._igraph import *
from igraph._igraph import __version__, __build_date__

import os
import math
import gzip
from tempfile import mkstemp


class Interval:
    """A class representing an interval over the real numbers"""
    def __init__(self, left, right):
	"""Constructor.

	@param left: the left edge of the interval
	@param right: the right edge of the interval"""
	if left == right:
	    raise ArgumentError, "Interval can't have zero length"""
	self._left = min(left, right)
	self._right = max(left, right)

    def __contains__(self, x):
	return (x >= self._left and x < self._right)

    def __eq__(self, x):
	return (isinstance(x, Interval) and \
		x.left == self._left and x.right == self._right)
 
    def __cmp__(self, x): return cmp(self._left, x.left)
    def __getattr__(self, attr):
	if attr == "left": return self._left
	if attr == "right": return self._right
	return object.__getattr__(self, attr)

    def __hash__(self): return hash(self._left) | hash(self._right)

    def __str__(self): return "[%f, %f)" % (self._left, self._right)

class Histogram:
    """Generic histogram class for real numbers"""

    def __init__(self, bin_width = 1, data = []):
	"""Initializes the histogram with the given data set.

	@param bin_width: the bin width of the histogram.
	@param data: the data set to be used. Must contain real numbers.
	"""
	self._bin_width = float(bin_width)
	self.clear()
	self.add_many(data)

    def _get_bin(self, num, create = False):
	"""Returns the bin corresponding to the given number.

	@param num: the number for which the bin is being sought
	@param create: whether to create a new bin if no bin exists yet.
	@return: the range of the bin or C{None} if no bin exists yet and
	  {create} is C{False}."""
	for bin in self._bins.keys():
	    if num in bin: return bin

	if create:
	    left = (num // self._bin_width) * self._bin_width
	    right = left + self._bin_width
	    bin = Interval(left, right)
	    self._bins[bin] = 0

	    if self._min is None or left < self._min: self._min = left
	    if self._max is None or right > self._max: self._max = right

	    return bin

	return None

    def add(self, num):
	"""Adds a single number to the histogram.
	
	@param num: the number to be added"""
	bin = self._get_bin(num, True)
	self._bins[bin] += 1

    def add_many(self, data):
	"""Adds a single number or elements of an iterable to the histogram.

	@param data: the data to be added"""
	try:
	    it = iter(data)
	except:
	    it = iter([data])
	for x in it: self.add(x)
    __lshift__ = add_many

    def clear(self):
	"""Clears the collected data"""
	self._bins = {}
	self._min = None
	self._max = None
	
    def __str__(self):
	"""Returns the string representation of the histogram"""
	if self._min is None or self._max is None: return str()
	num_length = max(len("%.3f" % self._min), \
			 len("%.3f" % self._max))
	format_string = "[%%%d.3f, %%%d.3f): %%s" % (num_length, num_length)
	#bins = self._bins.keys()
	#bins.sort()
	maxval = max(self._bins.itervalues())
	scale = maxval // (70-2*num_length)
	if scale<1: scale = 1

	result = []

	if scale>1: result.append("Each * represents %d items" % scale)

	x = self._min
	while x < self._max:
	    bin = self._get_bin(x)
	    if bin is None:
		cnt = 0
	    else:
		cnt = self._bins[bin] // scale
	    result.append(format_string % (x, x+self._bin_width, '*'*cnt))
	    x += self._bin_width

	return "\n".join(result)

class Graph(_igraph.Graph):
    def indegree(self, *args, **kwds):
	"""indegree(...)

	Returns the in-degrees in a list. See L{degree} for possible
	arguments.
	"""
	kwds['degree']=_igraph.IN
	return self.degree(*args, **kwds)

    def outdegree(self, *args, **kwds):
	"""outdegree(...)

	Returns the out-degrees in a list. See L{degree} for possible
	arguments.
	"""
	kwds['degree']=_igraph.OUT
	return self.degree(*args, **kwds)

    def eccentricity(self, nodes=None):
	"""eccentricity(self, vertices=None)

	Calculates eccentricities for vertices with the given indices.
	Eccentricity is given as the reciprocal of the greatest distance
	between the vertex being considered and any other vertex in the
	graph.

	Please note that for any unconnected graph, eccentricities will
	all be equal to 1 over the number of vertices, since for all vertices
	the greatest distance will be equal to the number of vertices (this
	is how C{shortest_paths} denotes vertex pairs where it is impossible
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

    def edge_betweenness_clustering(self, clusters = None, steps = None, \
				    return_graph = False, \
				    return_removed_edges = False):
	"""edge_betweenness_clustering(clusters = None, steps = None, return_graph = False, return_removed_edges = False)

	Newman's edge betweenness clustering.

	Iterative removal of edges with the largest edge betweenness until
	the given number of steps is reached or until the graph is decomposed
	to the given number of clusters. Edge betweennesses are recalculated
	after every run.

	For more information, see the original paper of Girvan and Newman:
	Girvan, M and Newman, MEJ: Community structure in social and
	biological networks. Proc. Natl. Acad. Sci. USA 99, 7821-7826 (2002)

	@param clusters: the desired number of clusters.
	@param steps: the number of tests to take.
	@param return_graph: whether to return the state of the graph when
	  the iteration ended.

        @return: the cluster index for every node. If C{return_graph} is
	  C{True}, the clustered graph is also returned. If
	  C{return_removed_edges} is C{True}, the list of removed edges is
	  also returned."""
	g = self.copy()
	number_of_steps = 0
	removed_edges = [None, []][return_removed_edges]

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

	    if return_removed_edges:
		removed_edges.append(g.es[eb_max_index].tuple)

	    g.delete_edges(eb_max_index, by_index=True)
	    number_of_steps += 1

	if return_graph and return_removed_edges:
	    return g.clusters(), g, removed_edges
	elif return_graph:
	    return g.clusters(), g
	elif return_removed_edges:
	    return g.clusters(), removed_edges
	return g.clusters()



    def write_graphmlz(self, f, compresslevel=9):
	"""write_graphmlz(f, compresslevel=9)
    
	Writes the graph to a zipped GraphML file.

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

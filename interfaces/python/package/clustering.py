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

	def __init__(self, membership):
		"""Constructor.

		@param membership: the membership list -- that is, the cluster
		  index in which each element of the set belongs to.
		"""
		self._membership = list(membership)
	
	def __getitem__(self, idx):
		"""Returns the members of the specified cluster.

		@param idx: the index of the cluster
		@return: the members of the specified cluster as a list
		"""
		return [i for i,e in enumerate(self._membership) if e==idx]

	def __len__(self):
		"""Returns the number of clusters.

		@note: the result is calculated by taking the difference between
		the minimum and maximum elements of the membership vector and adding
		it to 1. So if there are empty clusters, they are still counted in the
		number of clusters -- except when the clusters with the largest or smallest
		indices are empty!

		@return: the number of clusters
		"""
		if len(self._membership) == 0: return 0
		return max(self._membership) - min(self._membership) + 1

	def _get_membership(self): return self._membership
	membership = property(_get_membership, doc = "The membership vector (read only)")

	def size(self, idx):
		"""Returns the size of a given cluster.

		@param idx: the cluster in which we are interested.
		"""
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

	def __init__(self, graph, membership = None, modularity = None, *args, **kwds):
		"""Creates a clustering object for a given graph.

		@param graph: the graph that will be associated to the clustering
		@param membership: the membership list. The length of the list must
		  be equal to the number of vertices in the graph. If C{None}, every
		  vertex is assumed to belong to the same cluster.
		@param modularity: the modularity score of the clustering. If C{None},
		  it will be calculated.
		"""
		self._graph = graph

		if membership is None:
			Clustering.__init__(self, [0]*graph.vcount())
		else:
			if len(membership) != graph.vcount():
				raise ValueError, "membership list is too short"
			Clustering.__init__(self, membership)

		if modularity is None:
			self._q = 0 # TODO
		else:
			self._q = modularity

	def _get_modularity(self): return self._q
	modularity = property(_get_modularity, doc = "The modularity score")
	q = modularity


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

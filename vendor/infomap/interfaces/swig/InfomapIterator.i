%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/iterators/InfomapIterator.h"
%}

%include "std_deque.i"

namespace std {
    %template(deque_uint) std::deque<unsigned int>;
}

%include "InfoNode.i"

/* Parse the header file to generate wrappers */
%include "src/core/iterators/InfomapIterator.h"


#ifdef SWIGPYTHON
%extend infomap::InfomapIterator
{

	// Make the class iterable, and wait until
	// first call to next() to yield first element
	%insert("python") %{
		def __iter__(self):
			self._firstYielded = False
			return self

		def __next__(self):
			if not self._firstYielded:
				self._firstYielded = True
			else:
				self.stepForward()

			if self.isEnd():
				raise StopIteration

			return self


		@property
		def module_id(self):
			"""Get the module id of the node.

			Returns
			-------
			int
				The module id
			"""
			return self.moduleId()

		_path = path
		@property
		def path(self):
			"""Get the path to the node in the tree.

			Returns
			-------
			tuple of ints
				The path
			"""
			return self._path()

		_depth = depth
		@property
		def depth(self):
			"""Get the depth.

			Returns
			-------
			int
				The depth
			"""
			return self._depth()

		@property
		def modular_centrality(self):
			"""Get the modular centrality of the node.

			Returns
			-------
			float
				The modular centrality
			"""
			return self.modularCentrality()

		@property
		def child_index(self):
			"""Get the child index.

			Returns
			-------
			int
				The child index
			"""
			return self.childIndex()

		# Forward to the node it currently points to
		def __getattr__(self, name):
			return getattr(self.current(), name)

	%}
}



%extend infomap::InfomapLeafModuleIterator
{

	// Make the class iterable, and wait until
	// first call to next() to yield first element
	%insert("python") %{
		def __iter__(self):
			self._firstYielded = False
			return self

		def __next__(self):
			if not self._firstYielded:
				self._firstYielded = True
			else:
				self.stepForward()

			if self.isEnd():
				raise StopIteration

			return self

		@property
		def module_id(self):
			"""Get the module id of the node.

			Returns
			-------
			int
				The module id
			"""
			return self.moduleId()

		_path = path
		@property
		def path(self):
			"""Get the path to the node in the tree.

			Returns
			-------
			tuple of ints
				The path
			"""
			return self._path()

		_depth = depth
		@property
		def depth(self):
			"""Get the depth.

			Returns
			-------
			int
				The depth
			"""
			return self._depth()

		@property
		def modular_centrality(self):
			"""Get the modular centrality of the node.

			Returns
			-------
			float
				The modular centrality
			"""
			return self.modularCentrality()

		@property
		def child_index(self):
			"""Get the child index.

			Returns
			-------
			int
				The child index
			"""
			return self.childIndex()

		# Forward to the node it currently points to
		def __getattr__(self, name):
			return getattr(self.current(), name)

	%}
}



%extend infomap::InfomapLeafIterator
{

	// Make the class iterable, and wait until
	// first call to next() to yield first element
	%insert("python") %{
		def __iter__(self):
			self._firstYielded = False
			return self

		def __next__(self):
			if not self._firstYielded:
				self._firstYielded = True
			else:
				self.stepForward()

			if self.isEnd():
				raise StopIteration

			return self


		@property
		def module_id(self):
			"""Get the module id of the node.

			Returns
			-------
			int
				The module id
			"""
			return self.moduleId()

		_path = path
		@property
		def path(self):
			"""Get the path to the node in the tree.

			Returns
			-------
			tuple of ints
				The path
			"""
			return self._path()

		_depth = depth
		@property
		def depth(self):
			"""Get the depth.

			Returns
			-------
			int
				The depth
			"""
			return self._depth()

		@property
		def modular_centrality(self):
			"""Get the modular centrality of the node.

			Returns
			-------
			float
				The modular centrality
			"""
			return self.modularCentrality()

		@property
		def child_index(self):
			"""Get the child index.

			Returns
			-------
			int
				The child index
			"""
			return self.childIndex()

		# Forward to the node it currently points to
		def __getattr__(self, name):
			return getattr(self.current(), name)

	%}
}



%extend infomap::InfomapIteratorPhysical
{

	// Make the class iterable, and wait until
	// first call to next() to yield first element
	%insert("python") %{
		def __iter__(self):
			self._firstYielded = False
			return self

		def __next__(self):
			if not self._firstYielded:
				self._firstYielded = True
			else:
				self.stepForward()

			if self.isEnd():
				raise StopIteration

			return self


		@property
		def module_id(self):
			"""Get the module id of the node.

			Returns
			-------
			int
				The module id
			"""
			return self.moduleId()

		_path = path
		@property
		def path(self):
			"""Get the path to the node in the tree.

			Returns
			-------
			tuple of ints
				The path
			"""
			return self._path()

		_depth = depth
		@property
		def depth(self):
			"""Get the depth.

			Returns
			-------
			int
				The depth
			"""
			return self._depth()

		@property
		def modular_centrality(self):
			"""Get the modular centrality of the node.

			Returns
			-------
			float
				The modular centrality
			"""
			return self.modularCentrality()

		@property
		def child_index(self):
			"""Get the child index.

			Returns
			-------
			int
				The child index
			"""
			return self.childIndex()

		# Forward to the node it currently points to
		def __getattr__(self, name):
			return getattr(self.current(), name)

	%}
}



%extend infomap::InfomapLeafIteratorPhysical
{

	// Make the class iterable, and wait until
	// first call to next() to yield first element
	%insert("python") %{
		def __iter__(self):
			self._firstYielded = False
			return self

		def __next__(self):
			if not self._firstYielded:
				self._firstYielded = True
			else:
				self.stepForward()

			if self.isEnd():
				raise StopIteration

			return self


		@property
		def module_id(self):
			"""Get the module id of the node.

			Returns
			-------
			int
				The module id
			"""
			return self.moduleId()

		_path = path
		@property
		def path(self):
			"""Get the path to the node in the tree.

			Returns
			-------
			tuple of ints
				The path
			"""
			return self._path()

		_depth = depth
		@property
		def depth(self):
			"""Get the depth.

			Returns
			-------
			int
				The depth
			"""
			return self._depth()

		@property
		def modular_centrality(self):
			"""Get the modular centrality of the node.

			Returns
			-------
			float
				The modular centrality
			"""
			return self.modularCentrality()

		@property
		def child_index(self):
			"""Get the child index.

			Returns
			-------
			int
				The child index
			"""
			return self.childIndex()

		# Forward to the node it currently points to
		def __getattr__(self, name):
			return getattr(self.current(), name)

	%}
}
#endif
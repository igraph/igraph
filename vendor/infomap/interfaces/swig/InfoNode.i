%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/InfoNode.h"
%}

%include "FlowData.i"
%include "InfoEdge.i"
%include "std_string.i"
%include "std_vector.i"

namespace std {
    %template(vector_uint) std::vector<unsigned int>;
}

/* Parse the header file to generate wrappers */
%include "src/core/InfoNode.h"


#ifdef SWIGPYTHON
%extend infomap::InfoNode
{
	%insert("python") %{

		@property
		def node_id(self):
			"""Get the physical node id.

			Returns
			-------
			int
				The node id
			"""
			return self.physicalId

		@property
		def state_id(self):
			"""Get the state id of the node.

			Returns
			-------
			int
				The state id
			"""
			return self.stateId

		@property
		def flow(self):
			"""Get the flow of the node.

			Returns
			-------
			float
				The flow
			"""
			return self.data.flow

		@property
		def layer_id(self):
			"""Get the layer id of a multilayer node.

			Returns
			-------
			int
				The layer id
			"""
			return self.layerId

		@property
		def child_degree(self):
			"""The number of children.

			Returns
			-------
			int
				Number of children
			"""
			return self.childDegree()

		@property
		def is_leaf(self):
			"""True if the node has no children.

			Returns
			-------
			bool
				Is leaf node
			"""
			return self.isLeaf()

		@property
		def is_leaf_module(self):
			"""True if the node has children but no grandchildren.

			Returns
			-------
			bool
				Is leaf module
			"""
			return self.isLeafModule()

		@property
		def is_root(self):
			"""True if the node has no parent.

			Returns
			-------
			bool
				Is root
			"""
			return self.isRoot()

		@property
		def meta_data(self):
			"""Meta data (on first dimension if more).

			Returns
			-------
			int
				The meta data
			"""
			return self.getMetaData()

		def get_meta_data(self, dimension = 0):
			"""Get meta data on a specific dimension.

			Parameters
			----------
			dimension : int
				The dimension (default 0)

			Returns
			-------
			int
				The meta data
			"""
			return self.getMetaData(dimension)
	%}
}
#endif

%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/FlowData.h"
%}

/* Parse the header file to generate wrappers */
%include "src/core/FlowData.h"

#ifdef SWIGPYTHON
%extend infomap::FlowData
{
	%insert("python") %{

		@property
		def enter_flow(self):
			"""Get the flow entering the node."""
			return self.enterFlow
		
		@property
		def exit_flow(self):
			"""Get the flow exiting the node."""
			return self.exitFlow
		
	%}
}
#endif
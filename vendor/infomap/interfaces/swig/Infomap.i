%module infomap
#pragma SWIG nowarn=302,325,362,383,389,503,508,509

%include "exceptions.i"

%{
/* Includes the header in the wrapper code */
#include "src/Infomap.h"
// SWIG strips namespaces, so include infomap in global namespace in wrapper code
using namespace infomap;
%}

// Makes it possible to use for example numpy.int64 as link weights
%begin %{
#define SWIG_PYTHON_CAST_MODE
%}

%include "std_string.i"
%include "Config.i"
%include "InfomapBase.i"

%include "std_map.i"
%include "std_vector.i"
%include "std_pair.i"

namespace std {
    %template(vector_uint) std::vector<unsigned int>;
    %template(map_uint_uint) std::map<unsigned int, unsigned int>;
    %template(map_uint_vector_uint) std::map<unsigned int, std::vector<unsigned int>>;
    %template(map_uint_string) std::map<unsigned int, std::string>;
    %template(pair_uint_uint) std::pair<unsigned int, unsigned int>;
    %template(map_pair_uint_uint_double) std::map<std::pair<unsigned int, unsigned int>, double>;
}

#ifdef SWIGPYTHON
%extend infomap::InfomapWrapper
{
	// overrides to convert swig map proxy object to pure dict
	%insert("python") %{
        def getModules(self, level=1, states=False):
            return dict(_infomap.InfomapWrapper_getModules(self, level, states))

        def getMultilevelModules(self, states=False):
            return dict(_infomap.InfomapWrapper_getMultilevelModules(self, states))

        def getNames(self):
            return dict(_infomap.InfomapWrapper_getNames(self))

        def getLinks(self, flow=False):
            return dict(_infomap.InfomapWrapper_getLinks(self, flow))
    %}
}
#endif

/* Parse the header file to generate wrappers */
%include "src/Infomap.h"

#ifdef SWIGPYTHON
%pythoncode "interfaces/python/infomap.py"
#endif

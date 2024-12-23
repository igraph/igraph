%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/InfomapBase.h"
%}

%include "InfomapConfig.i"
%include "Config.i"
%include "InfoNode.i"
%include "InfomapIterator.i"
%include "Network.i"

%include "std_map.i"

namespace std {
    %template(map_uint_uint) std::map<unsigned int, unsigned int>;
    %template(vector_double) std::vector<double>;
}
namespace infomap {
    %template(InfomapConfigInfomapBase) InfomapConfig<InfomapBase>;
}

/* Parse the header file to generate wrappers */
%include "src/core/InfomapBase.h"

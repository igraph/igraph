%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/InfomapOptimizer.h"
%}

%include "InfomapBase.i"
%include "InfoNode.i"
%include "Config.i"

/* Parse the header file to generate wrappers */
%include "src/core/InfomapOptimizer.h"

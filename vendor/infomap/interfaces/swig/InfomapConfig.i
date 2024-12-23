%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/InfomapConfig.h"
%}

%include "Config.i"
%include "std_string.i"

/* Parse the header file to generate wrappers */
%include "src/core/InfomapConfig.h"

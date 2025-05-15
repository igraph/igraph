%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/InfoEdge.h"
%}

/* Parse the header file to generate wrappers */
%include "src/core/InfoEdge.h"

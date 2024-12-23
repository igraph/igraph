%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/core/MapEquation.h"
%}

%include "InfoNode.i"

/* Parse the header file to generate wrappers */
%include "src/core/MapEquation.h"

%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/io/Network.h"
%}

%include "std_string.i"
%include "StateNetwork.i"

/* Parse the header file to generate wrappers */
%include "src/io/Network.h"
%module infomap

%{
/* Includes the header in the wrapper code */
#include "src/utils/MetaCollection.h"
%}

/* Parse the header file to generate wrappers */
%include "src/utils/MetaCollection.h"
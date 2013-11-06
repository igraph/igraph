#!/bin/sh

## Convert HTML files created by R CMD INSTALL, to 
## file that can be used as Jekyll inputs

if [ "$#" -ne 2 ] || ! [ -d "$1" ] || ! [ -d "$2" ]; then
  echo "Usage: $0 SOURCE-DIR TARGET-DIR" >&2
  exit 1
fi

indir=$1
outdir=$2

inhtml=`ls $indir/*.html`

header='---
layout: manual
title: igraph R manual pages
mainheader: R igraph manual pages
lead: Use this if you are using igraph from R
---
'

for ih in ${inhtml}; do
    oh=${outdir}/`basename ${ih}`
    
    # YAML header
    echo "${header}" > ${oh}
    
    # Begin raw
    echo '{% raw %}' >> ${oh}
    
    # Main text
    cat ${ih} |
 
    # Strip header and footer    
    sed -n '1,/<body[ >]/!{ /<\/body>/,/<body[ >]/!p; }' |

    # Remove top table
    grep -v '^<table.*summary="page for ' | 

    # Done 
    cat >> ${oh}
    
    # End raw
    echo '{% endraw %}' >> ${oh}
done;

# Index file is treated separately, 
# here we also want to remove more from the top
ih=`dirname ${ih}`/00Index.html
oh=${outdir}/00Index.html
# YAML header
echo "${header}" > ${oh}
echo '{% raw %}' >> ${oh}
cat ${ih} |
sed -n '1,/<body[ >]/!{ /<\/body>/,/<body[ >]/!p; }' |
grep -v '^<table.*summary="page for ' |
sed -n '/<h2>Help Pages<\/h2>/,$p' | tail +2 |
cat >> ${oh}
echo '{% endraw %}' >> ${oh}


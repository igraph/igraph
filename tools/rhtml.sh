#!/bin/bash

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

basepkg="base boot class cluster codetools compiler datasets               \
    foreign graphics grDevices grid KernSmooth lattice MASS Matrix methods \
    mgcv nlme nnet parallel rpart spatial splines stats stats4 survival    \
    tcltk tools utils"

sedfile=`mktemp /tmp/XXXXXX`
echo -n > $sedfile
for i in $basepkg; do 
    printf "s/href=\\\"\\.\\.\\/\\.\\.\\/${i}\\/html\/\\([^\\\"]*\)\\.html/\\\"" >> $sedfile
    printf "href=\\\"http:\\/\\/www.inside-r.org\\/r-doc\\/${i}\\/\\\1/g\n" >> $sedfile
done

# echo foo | sed -f $sedfile 

for ih in ${inhtml}; do
    ihf=`basename ${ih}`
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

    # From index remove the top
    if [ ${ihf} = "00Index.html" ]; then 
	sed -n '/<h2>Help Pages<\/h2>/,$p' | tail +2 
    else
	cat
    fi |

    # Rewrite links to other packages, looks like 
    # we are only referring to base and recommended
    # packages, maybe only these are allowed by R CMD check
    sed -f $sedfile |
    
    # Done 
    cat >> ${oh}
    
    # End raw
    echo '{% endraw %}' >> ${oh}
done;

rm -f ${sedfile}

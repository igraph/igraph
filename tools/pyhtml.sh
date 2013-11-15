#!/bin/bash

## Convert HTML files created by epydoc, to files to be used as 
## Jekyll inputs

if [ "$#" -ne 1 ] || ! [ -d "$1" ]; then
  echo "Usage: $0 DIR" >&2
  exit 1
fi

dir=$1

inhtml=`ls $dir/*.html`

header='---
layout: epydoc
title: python-igraph manual
mainheader: python-igraph manual
lead: For using igraph from Python
---
'

tmpfile=`mktemp /tmp/XXXXXX`

for ih in ${inhtml}; do
    hf=`basename ${ih}`

    printf "%b" "Converting ${hf}..."
    
    # YAML header
    echo "${header}" > ${tmpfile}
    
    # Begin raw
    echo '{% raw %}' >> ${tmpfile}

    # Main text
    cat ${ih} |

    # Strip header and footer
    sed -n '1h;1!H;${;g;s/.*<body\([^>]*\)>/<div\1/g;p;}' |
    sed -n '1h;1!H;${;g;s/<\/body>.*/<\/div>/;p;}' |

    # Take out links to frames / no frames
    sed -n '1h;1!H;${;g;s/\[[^[]*frames.html[^[]*]//;p;}' |

    # Repair a bug in epydoc
    sed 's/^<\/div>\(<a name=\)/\1/' |

    # Remove link from navbar
    sed 's/<a class="navbar" [^>]*_top[^>]*igraph\.org[^>]*>[^<]*<\/a>//' |

    # Remove navbar class, replace with epynavbar
    sed 's/class="navbar"/class="epynavbar"/g' |

    # Done 
    cat >> ${tmpfile}
    
    # End raw
    echo '{% endraw %}' >> ${tmpfile}
    
    # Replace
    cp ${tmpfile} ${ih}
    
    printf "Done.\n"
    
done;

rm -f ${tmpfile}


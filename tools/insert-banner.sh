#!/bin/bash

## Insert a banner into a html file, right at the start of <body>

if [ $# != "2" -a $# != "3" ] || [ ! -d $1 ] || [ ! -f $2 ]; then
    printf "Usage: $0 <directory> <banner-file> [<exclude-file>]\n"
    exit 1
fi

banner=$2
exclude=/dev/null
if [ -n "$3" ]; then exclude=$3; fi

tmpfile=`mktemp -t XXXXXX`

function insert {
    printf "%b" "Doing $1..."

    if [ -n "$exclude" ] && (echo $1 | grep -q -f $exclude); then
        printf "%b" " excluded\n"
    else
        insert2 $1 > "$tmpfile"
        cp $tmpfile $1
        printf "%b" " DONE\n"
    fi
}

function insert2  {
    cat $1 |
    sed -n '1h;1!H;${;g;s/\(<body[^>]*>\)/\1\n<!--endbanner-->/g;p;}' |
    sed "/<!--endbanner-->/ {
        r $banner
        N
    }"
}

find $1 -name "*.html" |
while read; do
    insert $REPLY
done

rm "$tmpfile"


#! /bin/sh

vstring=`cat ../../VERSION`
if echo "$vstring" | grep -q '\+'; then
    pvstring=`echo $vstring | 
         sed 's/\(^[0-9\.]*\)[^+]*[+]*\([0-9][0-9]*\)\.\([0-9a-f]*\)$/\1|\2|\3/'`
    version=`echo $pvstring | cut -f1 -d"|"`
    dist=-`echo $pvstring | cut -f2 -d"|"`
    hash=`echo $pvstring | cut -f3 -d"|"`

    ## Erase trailing zeros from version
    while echo $version | grep -q '\.0'; do
	version=`echo $version | sed 's/\.0$//'`
    done

    ## Decrease last number
    if echo $version | grep -q '\.'; then
	lastnum=`echo $version | sed 's/.*\.\([0-9]*\)$/\1/'`
	rest=`echo $version | sed 's/\(.*\.\)[0-9]*$/\1/'`
    else 
	lastnum=$version
	rest=""
    fi
    
    ## Add .999 and distance
    echo ${rest}$((lastnum-1)).999${dist}

else 
    echo $vstring
fi

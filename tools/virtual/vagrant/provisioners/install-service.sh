#! /bin/bash

## Quit immediately on error
set -e

## Want to run this as 'vagrant', so rerun if root
if [ "$(id -u)" = "0" ]; then
    sudo -u vagrant bash $0 $@
    exit 0
fi

scriptdir="/tekton"

## Check arguments, we need two of them, 
## first specifies the script to run, the second
## is the crontab specification for when to run it.
if [ ! "x"$# = "x2" ]; then exit 2; fi
script=$1
if [ ! -e ${scriptdir}/${script} ]; then exit 3; fi




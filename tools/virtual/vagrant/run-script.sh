#! /bin/sh

## Quit immediately on error
set -e

## Check arguments, at least the script to run is needed.
## Additional arguments will be passed to the script
if [ $# -lt 1 ]; then 
    echo "Error: not enough arguments, need script to run at least"
    exit 1
fi
script=$1 
shift

if [ ! -f "../scripts/$script" ]; then
    echo "Script '$script' does not exist"
    exit 3
fi

if [ ! -x "../scripts/$script" ]; then
    echo "Script '$script' is not executable"
    exit 4
fi

vagrant up
vagrant ssh -- /tekton/$script $@

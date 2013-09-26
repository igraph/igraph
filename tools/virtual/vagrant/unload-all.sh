#! /bin/bash

for a in `launchctl list | cut -f3 | grep "^org\.igraph\.tekton\."`; do
    launchctl remove $a
done

#! /bin/bash

launchctl list | cut -f3 | grep "^org\.igraph\.tekton\." | xargs launchctl remove


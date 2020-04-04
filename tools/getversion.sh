#! /bin/bash

thistag=$(git describe --exact-match --tags HEAD 2>/dev/null || true)

if [ -z "${thistag}" ]; then
    script_dir=$( dirname "${BASH_SOURCE[0]}" )
    if [ -f $script_dir/NEXT_VERSION ]; then
        # taghash=$(git rev-list --tags --max-count=1)
        # tag=$(git describe --tags "$taghash")
        next_version=$( cat $script_dir/NEXT_VERSION )
        current=$(git rev-parse --short HEAD)
        echo "${next_version}-pre+${current}"
    else
        latest_version=$(git describe --abbrev=0 HEAD)
        current=$(git rev-parse --short HEAD)
        echo "${latest_version}-post+${current}"
    fi
else
    echo "${thistag}"
fi

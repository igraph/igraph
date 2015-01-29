#! /bin/bash

thistag=$(git describe --exact-match --tags HEAD 2>/dev/null || true)

if [ -z "${thistag}" ]; then
    # taghash=$(git rev-list --tags --max-count=1)
    # tag=$(git describe --tags "$taghash")
    next_version=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && cat NEXT_VERSION )
    current=$(git rev-parse --short HEAD)
    echo "${next_version}-pre+${current}"
else
    echo "${thistag}"
fi

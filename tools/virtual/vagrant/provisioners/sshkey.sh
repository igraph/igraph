#! /bin/sh

#! /bin/sh

## Quit immediately on error
set -e

## Want to run this as 'vagrant', so rerun if root
if [ "$(id -u)" = "0" ]; then
    sudo -u vagrant bash $0 $@
    exit 0
fi

mkdir -p ~/.ssh
cp /tekton/key/id_rsa ~/.ssh/id_rsa
cp /tekton/key/known_hosts ~/.ssh/

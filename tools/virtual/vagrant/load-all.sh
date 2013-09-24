#! /bin/sh

agentdir=`dirname $0`/agents
agents=`ls $agentdir/*.plist`
tmpdir=`mktemp -d -t tekton`
mywd=$(pwd | sed -e 's/[\/&]/\\&/g')
trap "rm -rf $tmpdir" EXIT

for agent in $agents; do
    agfile=$tmpdir/`basename $agent`    
    basename $agent
    cat "$agent" | sed 's/\$PWD/'$mywd'/g' >> "$agfile"
    launchctl load "$agfile"
done

rm -rf $tmpdir


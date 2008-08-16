#! /bin/bash

export repohost=cneurocvs.rmki.kfki.hu
export repodir=igraph-new
export debianrepodir=/var/www/packages

# Get current version
version="`head -1 configure.in | cut -f2 -d, | tr -d ' '`"
nextversion=`echo $version+0.1 | bc | sed s/^\./0\./` 

# Create directories on the server
ssh $repohost mkdir -p ${repodir}/download
ssh $repohost mkdir -p ${repodir}/doc-${version}/{R,html,python}
ssh $repohost mkdir -p ${repodir}/images/screenshots
ssh $repohost rm -f ${repodir}/doc
ssh $repohost ln -s doc-${version} ${repodir}/doc

# Make everything 
./configure &&
make &&
make check &&
cd doc && make html && make pdf && make info && cd .. || exit 1

#################################################
# start jed to write the NEWS section
jed NEWS

#################################################
# make a tar.gz source distribution and upload it to the 
#       igraph homepage
make dist || exit 1
scp igraph-${version}.tar.gz ${repohost}:${repodir}/download/ || exit 1

#################################################
# make documentation and upload it to the igraph homepage
scp -rC doc/html ${repohost}:${repodir}/doc/ && 
cd doc && mv html igraph-docs-${version} &&
tar czf igraph-docs-${version}.tar.gz igraph-docs-${version} && 
mv igraph-docs-${version} html && cd - && 
scp doc/igraph-docs-${version}.tar.gz ${repohost}:${repodir}/doc/ && 
scp -C doc/igraph-docs.info doc/igraph-docs.pdf ${repohost}:${repodir}/doc/ ||
exit 1

#################################################
# update the igraph homepage (html) itself

cd doc/homepage
./generate.py $version
scp -r * ${repohost}:${repodir}/
cd ../..

#################################################
# symlink the uploaded .tar.gz to the proper .orig.tar.gz name
# for debian source package creation
mkdir ../debian-package &&
ln igraph-${version}.tar.gz ../debian-package/igraph_${version}.orig.tar.gz || exit 1

#################################################
# extract the original source .tar.gz to the proper directory
# (debian needs this -- we don't want our arch specific files to
# be included in the diff file created by dpkg-buildpackage
cd ../debian-package &&
tar -xvvzf igraph_${version}.orig.tar.gz &&
cd igraph-${version} || exit 1

#################################################
# make debian packages and upload them to the igraph homepage
make deb &&
scp ../libigraph_${version}_i386.deb ${repohost}:${debianrepodir}/binary/ &&
scp ../libigraph-dev_${version}_i386.deb ${repohost}:${debianrepodir}/binary/ &&
scp ../igraph_${version}.orig.tar.gz ${repohost}:${debianrepodir}/source/ &&
scp ../igraph_${version}.diff.gz ${repohost}:${debianrepodir}/source/ &&
scp ../igraph_${version}.dsc ${repohost}:${debianrepodir}/source/ || exit 1

######################
# clean up some stuff
cd ../.. && rm -rf debian-package || exit 1

#################################################
# make an R source package and upload that too

cd interfaces/R
make
echo Please upload the package to the windows build server...

scp igraph_${version}.tar.gz ${repohost}:${repodir}/download/ &&
scp igraph_${version}.zip ${repohost}:${repodir}/download/ &&
scp igraph_${version}.tgz ${repohost}:${repodir}/download/ ||
exit 1
cd ../..

#################################################
# Upload R documentation

cp interfaces/R/igraph_${version}.tar.gz /tmp
cd /tmp
tar xzf igraph_${version}.tar.gz
R CMD check -l ~/.R/library igraph
scp -rC ~/.R/library/igraph/html/* ${repohost}:${repodir}/doc/R/
R CMD Rd2dvi --pdf igraph
scp igraph.pdf ${repohost}:${repodir}/doc/R/igraph.pdf
cd -

#################################################
# mirror the tla repo

# tla archive-mirror csardi@rmki.kfki.hu--2004-public

cat <<EOF
DONE. Don't forget to 
 * upload the Python package to the Python Package Index
 * create a tla tag for the release:
   tla tag igraph--main--${version} igraph--release--${version}
 * create a new tla version:
   tla tag igraph--main--${version} igraph--main--$nextversion
 * increase the version number:
   in configure.in replace ${version} with $nextversion
   tla commit
EOF

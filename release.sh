#! /bin/bash

repohost=cneurocvs.rmki.kfki.hu
repodir=/var/www/igraph.bak

# Get current version
version="`head -1 configure.in | cut -f2 -d, | tr -d ' '`"
nextversion=`echo $version+0.1 | bc | sed s/^\./0\./` 

# Make everything 
./configure --enable-python --enable-docs --enable-rpackage &&
make &&
make html && make pdf && make info || exit 1

#################################################
# start jed to write the NEWS section
jed NEWS

#################################################
# make a tar.gz source distribution and upload it to the 
#       igraph homepage
make dist || exit 1
scp igraph-${version}.tar.gz ${repohost}:${repodir}/download || exit 1

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

# insert the news
python > doc/igraph2.html -c '
import re
newsfile=open("NEWS", "r")
news=newsfile.read()
newsfile.close()

reg=re.compile(r"^=====*\n(?P<title>[^\n]*)\n====*\n", re.MULTILINE|re.DOTALL)
news=reg.sub("<table class=\"navigation\"><tr><td><h3>\g<title></h3></td></tr></table>\n", news)

reg=re.compile(r"^- (?P<entry>.*?)(?=(^$)|(^-[ ])|(\Z))", 
  re.MULTILINE|re.DOTALL)
news=reg.sub("<li>\g<entry></li>\n", news)
reg=re.compile(r"^[\s]*$", re.MULTILINE)
news=reg.sub("<p></p>", news)
reg=re.compile(r"</li>\s*<p>", re.MULTILINE)
news=reg.sub("</li></ul><p>", news)
reg=re.compile(r"</p>\s*<li>", re.MULTILINE)
news=reg.sub("</p><ul><li>", news)

htmlfile=open("doc/igraph.html")
html=htmlfile.read()
htmlfile.close()

reg=re.compile(r"(?P<head><!-- NEWS.*?-->)")
html=reg.sub("\g<head>"+news, html)

print html
' || exit 1

# replace the things to be replaced
sed 's/@VERSION@/'$version'/g' doc/igraph2.html >doc/igraph3.html && 
scp doc/igraph3.html ${repohost}:${repodir}/igraph.html && 
scp doc/*.png doc/*.jpg ${repohost}:${repodir}/ || exit 1

#################################################
make debian packages and upload them to the igraph homepage
make deb &&
scp ../*igraph*.deb ${repohost}:${repodir}/debian

#################################################
# make an R source package and upload that too

scp interfaces/R/igraph_${version}.tar.gz ${repohost}:${repodir}/download &&
scp interfaces/R/igraph_${version}.zip ${repohost}:${repodir}/download ||
exit 1

#################################################
# Upload R documentation

cp interfaces/R/igraph_${version}.tar.gz /tmp
cd /tmp
tar xzf igraph_${version}.tar.gz
R CMD check igraph
scp -rC igraph.Rcheck/igraph/html/* ${repohost}:${repodir}/doc/R/
dvipdf igraph.Rcheck/igraph-manual.dvi
scp igraph-manual.pdf ${repohost}:${repodir}/doc/R/igraph.pdf

#################################################
# mirror the tla repo

# tla archive-mirror csardi@rmki.kfki.hu--2004-public

cat <<EOF
DONE. Don't forget to 
 * create a tla tag for the release:
   tla tag igraph--main--${version} igraph--release--${version}
 * create a new tla version:
   tla tag igraph--main--${version} igraph--main--$nextversion
 * increase the version number:
   in configure.in replace ${version} with $nextversion
   tla commit
EOF

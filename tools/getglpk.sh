#! /bin/sh

if [ -d ../optional/glpk ]; then 
    echo "GLPK directory '../optional/glpk' already exists, remove it first"
#    exit 1
fi

THIS=`pwd`
IDIR=${THIS}/../optional/glpk/
mkdir $IDIR

GLPK="http://ftp.gnu.org/gnu/glpk/glpk-4.45.tar.gz"
TARGZ=`echo $GLPK | sed 's/^.*\///'`
DIR=`echo $TARGZ | sed 's/\.tar\.gz$//'`

cd /tmp
if [ ! -f $TARGZ ]; then curl -O $GLPK; fi
tar xzf $TARGZ

#cp -R $DIR/include/*.h $DIR/src/*.{c,h} $DIR/src/amd $DIR/src/colamd \
#    $DIR/{README,COPYING} $IDIR

cd $THIS

SRC=`ls ../optional/glpk/*.h ../optional/glpk/*.c`
SRC2=`ls ../optional/glpk/amd/*.h ../optional/glpk/amd/*.c`
SRC3=`ls ../optional/glpk/colamd/*.h ../optional/glpk/colamd/*.c`

INC=$IDIR/glpk.inc

/bin/echo -n "GLPK = " > $INC
for i in $SRC; do /bin/echo -n "$i " >>$INC; done
for i in $SRC2; do /bin/echo -n "$i " >>$INC; done
for i in $SRC3; do /bin/echo -n "$i " >>$INC; done


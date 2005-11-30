#!/bin/bash
if [ x$1 == x-d -a -x /usr/bin/valgrind ]; then
  # Checking memory leaks with Valgrind
  echo "Valgrind memory leak debugging enabled"
  FNAME=/tmp/igraph_${RANDOM}.supp
  cat /usr/lib/valgrind/python.supp >$FNAME
  cat $0|awk 'BEGIN { ok=0 } /[S]UPPRESSIONS/ { first=1 } { if (first) { ok=1; first=0; } else if (ok) { print; } }' >>$FNAME
  PRE="valgrind --tool=memcheck --leak-check=yes --trace-children=yes --suppressions=$FNAME"
  shift
else
  PRE=""
fi
LD_LIBRARY_PATH=src/.libs $PRE python $1 interfaces/python/test.py
if [ x$FNAME != x ]; then rm -f $FNAME; fi

exit 0
################## SUPPRESSIONS for Valgrind ########################
{
   <insert a suppression name here>
   Memcheck:Cond
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
}
{
   <insert a suppression name here>
   Memcheck:Cond
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/tls/i686/cmov/libc-2.3.5.so
   obj:/lib/ld-2.3.5.so
   fun:_dl_open
   obj:/lib/tls/i686/cmov/libdl-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/tls/i686/cmov/libdl-2.3.5.so
   fun:dlopen
   fun:_PyImport_GetDynLoadFunc
   fun:_PyImport_LoadDynamicModule
}
{
   <insert a suppression name here>
   Memcheck:Cond
   obj:/lib/ld-2.3.5.so
   obj:/lib/tls/i686/cmov/libc-2.3.5.so
   obj:/lib/ld-2.3.5.so
   fun:_dl_open
   obj:/lib/tls/i686/cmov/libdl-2.3.5.so
   obj:/lib/ld-2.3.5.so
   obj:/lib/tls/i686/cmov/libdl-2.3.5.so
   fun:dlopen
   fun:_PyImport_GetDynLoadFunc
   fun:_PyImport_LoadDynamicModule
   obj:/usr/bin/python2.4
   obj:/usr/bin/python2.4
}
					  

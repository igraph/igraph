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

# Need a patch to get rid of an abort() call. We call igraph_error()
# instead.

patch -p1 -d ../optional/glpk <<-EOF
--- glpk.old/glpenv01.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpenv01.c	2012-03-30 12:03:54.000000000 -0400
@@ -23,6 +23,7 @@
 ***********************************************************************/
 
 #include "glpapi.h"
+#include "igraph_error.h"
 
 /***********************************************************************
 *  NAME
@@ -126,19 +127,15 @@
       {  /* not initialized yet; perform initialization */
          if (glp_init_env() != 0)
          {  /* initialization failed; display an error message */
-            fprintf(stderr, "GLPK initialization failed\n");
-            fflush(stderr);
-            /* and abnormally terminate the program */
-            abort();
+	   IGRAPH_ERROR("GLPK initialization failed", IGRAPH_EGLP);
          }
          /* initialization successful; retrieve the pointer */
          env = tls_get_ptr();
       }
       /* check if the environment block is valid */
       if (env->magic != ENV_MAGIC)
-      {  fprintf(stderr, "Invalid GLPK environment\n");
-         fflush(stderr);
-         abort();
+      {  
+	IGRAPH_ERROR("Invalid GLPK environment", IGRAPH_EGLP);
       }
       return env;
 }
@@ -200,9 +197,8 @@
       if (env == NULL) return 1;
       /* check if the environment block is valid */
       if (env->magic != ENV_MAGIC)
-      {  fprintf(stderr, "Invalid GLPK environment\n");
-         fflush(stderr);
-         abort();
+      {  
+	 IGRAPH_ERROR("Invalid GLPK environment", IGRAPH_EGLP);
       }
       /* close handles to shared libraries */
       if (env->h_odbc != NULL)
--- glpk.old/glpenv04.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpenv04.c	2012-03-30 11:48:28.000000000 -0400
@@ -23,6 +23,7 @@
 ***********************************************************************/
 
 #include "glpapi.h"
+#include "../../include/igraph_error.h"
 
 /***********************************************************************
 *  NAME
@@ -44,14 +45,7 @@
       va_list arg;
       env->term_out = GLP_ON;
       va_start(arg, fmt);
-      xvprintf(fmt, arg);
-      va_end(arg);
-      xprintf("Error detected in file %s at line %d\n", env->err_file,
-         env->err_line);
-      if (env->err_hook != NULL)
-         env->err_hook(env->err_info);
-      abort();
-      exit(EXIT_FAILURE);
+      igraph_errorvf(fmt, env->err_file, env->err_line, IGRAPH_EGLP, arg);
       /* no return */
 }
EOF

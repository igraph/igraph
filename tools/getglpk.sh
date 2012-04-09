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
diff -ru glpk.old/glpenv01.c glpk/glpenv01.c
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
diff -ru glpk.old/glpenv03.c glpk/glpenv03.c
--- glpk.old/glpenv03.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpenv03.c	2012-04-02 11:18:42.000000000 -0400
@@ -40,9 +40,9 @@
 
 void glp_printf(const char *fmt, ...)
 {     va_list arg;
-      va_start(arg, fmt);
-      xvprintf(fmt, arg);
-      va_end(arg);
+      /* va_start(arg, fmt); */
+      /* xvprintf(fmt, arg); */
+      /* va_end(arg); */
       return;
 }
 
@@ -64,22 +64,22 @@
 void glp_vprintf(const char *fmt, va_list arg)
 {     ENV *env = get_env_ptr();
       /* if terminal output is disabled, do nothing */
-      if (!env->term_out) goto skip;
-      /* format the output */
-      vsprintf(env->term_buf, fmt, arg);
-      /* pass the output to the user-defined routine */
-      if (env->term_hook != NULL)
-      {  if (env->term_hook(env->term_info, env->term_buf) != 0)
-            goto skip;
-      }
-      /* send the output to the terminal */
-      fputs(env->term_buf, stdout);
-      fflush(stdout);
-      /* copy the output to the text file */
-      if (env->tee_file != NULL)
-      {  fputs(env->term_buf, env->tee_file);
-         fflush(env->tee_file);
-      }
+      /* if (!env->term_out) goto skip; */
+      /* /\* format the output *\/ */
+      /* vsprintf(env->term_buf, fmt, arg); */
+      /* /\* pass the output to the user-defined routine *\/ */
+      /* if (env->term_hook != NULL) */
+      /* {  if (env->term_hook(env->term_info, env->term_buf) != 0) */
+      /*       goto skip; */
+      /* } */
+      /* /\* send the output to the terminal *\/ */
+      /* fputs(env->term_buf, stdout); */
+      /* fflush(stdout); */
+      /* /\* copy the output to the text file *\/ */
+      /* if (env->tee_file != NULL) */
+      /* {  fputs(env->term_buf, env->tee_file); */
+      /*    fflush(env->tee_file); */
+      /* } */
 skip: return;
 }
 
diff -ru glpk.old/glpenv04.c glpk/glpenv04.c
--- glpk.old/glpenv04.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpenv04.c	2012-03-30 11:56:41.000000000 -0400
@@ -23,6 +23,7 @@
 ***********************************************************************/
 
 #include "glpapi.h"
+#include "igraph_error.h"
 
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
 
diff -ru glpk.old/glpenv07.c glpk/glpenv07.c
--- glpk.old/glpenv07.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpenv07.c	2012-03-31 13:21:03.000000000 -0400
@@ -413,13 +413,13 @@
 
 static void *c_fopen(const char *fname, const char *mode)
 {     FILE *fh;
-      if (strcmp(fname, "/dev/stdin") == 0)
-         fh = stdin;
-      else if (strcmp(fname, "/dev/stdout") == 0)
-         fh = stdout;
-      else if (strcmp(fname, "/dev/stderr") == 0)
-         fh = stderr;
-      else
+      /* if (strcmp(fname, "/dev/stdin") == 0) */
+      /*    fh = stdin; */
+      /* else if (strcmp(fname, "/dev/stdout") == 0) */
+      /*    fh = stdout; */
+      /* else if (strcmp(fname, "/dev/stderr") == 0) */
+      /*    fh = stderr; */
+      /* else */
          fh = fopen(fname, mode);
       if (fh == NULL)
          lib_err_msg(strerror(errno));
@@ -484,11 +484,11 @@
 static int c_fclose(void *_fh)
 {     FILE *fh = _fh;
       int ret;
-      if (fh == stdin)
-         ret = 0;
-      else if (fh == stdout || fh == stderr)
-         fflush(fh), ret = 0;
-      else
+      /* if (fh == stdin) */
+      /*    ret = 0; */
+      /* else if (fh == stdout || fh == stderr) */
+      /*    fflush(fh), ret = 0; */
+      /* else */
          ret = fclose(fh);
       if (ret != 0)
       {  lib_err_msg(strerror(errno));
diff -ru glpk.old/glpgmp.c glpk/glpgmp.c
--- glpk.old/glpgmp.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpgmp.c	2012-04-01 00:05:13.000000000 -0400
@@ -860,7 +860,7 @@
          d[j] = (unsigned char)r->val;
       }
       /* output the integer to the stream */
-      if (fp == NULL) fp = stdout;
+      /* if (fp == NULL) fp = stdout; */
       if (mpz_sgn(x) < 0)
          fputc('-', fp), nwr++;
       for (j = n-1; j >= 0; j--)
@@ -1091,7 +1091,7 @@
       int nwr;
       if (!(2 <= base && base <= 36))
          xfault("mpq_out_str: base = %d; invalid base\n", base);
-      if (fp == NULL) fp = stdout;
+      /* if (fp == NULL) fp = stdout; */
       nwr = mpz_out_str(fp, base, &x->p);
       if (x->q.val == 1 && x->q.ptr == NULL)
          ;
diff -ru glpk.old/glpmpl04.c glpk/glpmpl04.c
--- glpk.old/glpmpl04.c	2012-03-30 11:30:58.000000000 -0400
+++ glpk/glpmpl04.c	2012-04-01 00:07:09.000000000 -0400
@@ -341,11 +341,11 @@
 
 void open_output(MPL *mpl, char *file)
 {     xassert(mpl->out_fp == NULL);
-      if (file == NULL)
-      {  file = "<stdout>";
-         mpl->out_fp = (void *)stdout;
-      }
-      else
+      /* if (file == NULL) */
+      /* {  file = "<stdout>"; */
+      /*    mpl->out_fp = (void *)stdout; */
+      /* } */
+      /* else */
       {  mpl->out_fp = xfopen(file, "w");
          if (mpl->out_fp == NULL)
             error(mpl, "unable to create %s - %s", file, xerrmsg());
@@ -362,9 +362,9 @@
 
 void write_char(MPL *mpl, int c)
 {     xassert(mpl->out_fp != NULL);
-      if (mpl->out_fp == (void *)stdout)
-         xprintf("%c", c);
-      else
+      /* if (mpl->out_fp == (void *)stdout) */
+      /*    xprintf("%c", c); */
+      /* else */
          xfprintf(mpl->out_fp, "%c", c);
       return;
 }
@@ -393,7 +393,7 @@
 
 void flush_output(MPL *mpl)
 {     xassert(mpl->out_fp != NULL);
-      if (mpl->out_fp != (void *)stdout)
+      /* if (mpl->out_fp != (void *)stdout) */
       {  xfflush(mpl->out_fp);
          if (xferror(mpl->out_fp))
             error(mpl, "write error on %s - %s", mpl->out_file,
@@ -1410,7 +1410,7 @@
       if (mpl->row != NULL) xfree(mpl->row);
       if (mpl->col != NULL) xfree(mpl->col);
       if (mpl->in_fp != NULL) xfclose(mpl->in_fp);
-      if (mpl->out_fp != NULL && mpl->out_fp != (void *)stdout)
+      if (mpl->out_fp != NULL /* && mpl->out_fp != (void *)stdout */)
          xfclose(mpl->out_fp);
       if (mpl->out_file != NULL) xfree(mpl->out_file);
       if (mpl->prt_fp != NULL) xfclose(mpl->prt_fp);
EOF

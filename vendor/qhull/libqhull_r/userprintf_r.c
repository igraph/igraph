/*<html><pre>  -<a                             href="qh-user_r.htm"
  >-------------------------------</a><a name="TOP">-</a>

  userprintf_r.c
  user redefinable function -- qh_fprintf

  see README.txt  see COPYING.txt for copyright information.

  If you recompile and load this file, then userprintf_r.o will not be loaded
  from qhull_r.a or qhull_r.lib

  See libqhull_r.h for data structures, macros, and user-callable functions.
  See user_r.c for qhull-related, redefinable functions
  see user_r.h for user-definable constants
  See usermem_r.c for qh_exit(), qh_free(), and qh_malloc()
  see Qhull.cpp and RboxPoints.cpp for examples.

  qh_printf is a good location for debugging traps, checked on each log line

  Please report any errors that you fix to qhull@qhull.org
*/

#include "libqhull_r.h"
#include "poly_r.h" /* for qh.tracefacet */

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>

/*-<a                             href="qh-user_r.htm#TOC"
  >-------------------------------</a><a name="qh_fprintf">-</a>

  qh_fprintf(qh, fp, msgcode, format, list of args )
    print arguments to *fp according to format
    Use qh_fprintf_rbox() for rboxlib_r.c

  notes:
    sets qh.last_errcode if msgcode is error 6000..6999
    same as fprintf()
    fgets() is not trapped like fprintf()
    exit qh_fprintf via qh_errexit()
    may be called for errors in qh_initstatistics and qh_meminit
*/

void qh_fprintf(qhT *qh, FILE *fp, int msgcode, const char *fmt, ... ) {
    va_list args;
    if (qh) {
        if (msgcode >= MSG_ERROR && msgcode < MSG_WARNING) {
            qh->last_errcode = msgcode;
        }
        if (qh->FLUSHprint) {
            fflush(fp);
        }
    }
    if (!fp) {
        return;
    }
    if ((qh && qh->ANNOTATEoutput) || msgcode < MSG_TRACE4) {
        fprintf(fp, "[QH%.4d]", msgcode);
    } else if (msgcode >= MSG_ERROR && msgcode < MSG_STDERR ) {
        fprintf(fp, "QH%.4d ", msgcode);
    }
    va_start(args, fmt);
    vfprintf(fp, fmt, args);
    va_end(args);

} /* qh_fprintf */


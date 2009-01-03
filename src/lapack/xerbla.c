/* xerbla.f -- translated by f2c (version 19991025).
   You must link the resulting object file with the libraries:
	-lf2c -lm   (in that order)
*/

#include "f2c.h"
#include "config.h"
#include "arpack_internal.h"

/* Table of constant values */

static integer c__1 = 1;

/* Subroutine */ int igraphxerbla_(srname, info, srname_len)
char *srname;
integer *info;
ftnlen srname_len;
{
    /* Format strings */
    static char fmt_9999[] = "(\002 ** On entry to \002,a6,\002 parameter nu\
mber \002,i2,\002 had \002,\002an illegal value\002)";

    /* Builtin functions */
    integer s_wsfe(), do_fio(), e_wsfe();
    /* Subroutine */ int s_stop();

    /* Fortran I/O blocks */
    static cilist io___1 = { 0, 6, 0, fmt_9999, 0 };



/*  -- LAPACK auxiliary routine (version 2.0) -- */
/*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd., */
/*     Courant Institute, Argonne National Lab, and Rice University */
/*     September 30, 1994 */

/*     .. Scalar Arguments .. */
/*     .. */

/*  Purpose */
/*  ======= */

/*  XERBLA  is an error handler for the LAPACK routines. */
/*  It is called by an LAPACK routine if an input parameter has an */
/*  invalid value.  A message is printed and execution stops. */

/*  Installers may consider modifying the STOP statement in order to */
/*  call system-specific exception-handling facilities. */

/*  Arguments */
/*  ========= */

/*  SRNAME  (input) CHARACTER*6 */
/*          The name of the routine which called XERBLA. */

/*  INFO    (input) INTEGER */
/*          The position of the invalid parameter in the parameter list */
/*          of the calling routine. */

/* ===================================================================== */

/*     .. Executable Statements .. */

/*     s_wsfe(&io___1); */
/*     do_fio(&c__1, srname, (ftnlen)6); */
/*     do_fio(&c__1, (char *)&(*info), (ftnlen)sizeof(integer)); */
/*     e_wsfe(); */

/*     s_stop("", (ftnlen)0); */


/*     End of XERBLA */

} /* igraphxerbla_ */


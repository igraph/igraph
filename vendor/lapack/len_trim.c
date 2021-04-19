/*  -- translated by f2c (version 20191129).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://www.netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"


/*  -- LEN_TRIM is Fortran 95, so we use a replacement here */

integer igraphlen_trim__(char *s, ftnlen s_len)
{
    /* System generated locals */
    integer ret_val;

    /* Builtin functions */
    integer i_len(char *, ftnlen);




    for (ret_val = i_len(s, s_len); ret_val >= 1; --ret_val) {
	if (*(unsigned char *)&s[ret_val - 1] != ' ') {
	    return ret_val;
	}
    }
    return ret_val;
} /* igraphlen_trim__ */


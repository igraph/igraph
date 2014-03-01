/* ========================================================================= */
/* === AMD_control ========================================================= */
/* ========================================================================= */

/* ------------------------------------------------------------------------- */
/* AMD, Copyright (c) Timothy A. Davis,					     */
/* Patrick R. Amestoy, and Iain S. Duff.  See ../README.txt for License.     */
/* email: DrTimothyAldenDavis@gmail.com                                      */
/* ------------------------------------------------------------------------- */

/* User-callable.  Prints the control parameters for AMD.  See amd.h
 * for details.  If the Control array is not present, the defaults are
 * printed instead.
 */

#include "amd_internal.h"

GLOBAL void AMD_control
(
    double Control [ ]
)
{
    double alpha ;
    Int aggressive ;

    if (Control != (double *) NULL)
    {
	alpha = Control [AMD_DENSE] ;
	aggressive = Control [AMD_AGGRESSIVE] != 0 ;
    }
    else
    {
	alpha = AMD_DEFAULT_DENSE ;
	aggressive = AMD_DEFAULT_AGGRESSIVE ;
    }

    PRINTF (("\nAMD version %d.%d.%d, %s: approximate minimum degree ordering\n"
	"    dense row parameter: %g\n", AMD_MAIN_VERSION, AMD_SUB_VERSION,
	AMD_SUBSUB_VERSION, AMD_DATE, alpha)) ;

    if (alpha < 0)
    {
	PRINTF (("    no rows treated as dense\n")) ;
    }
    else
    {
	PRINTF ((
	"    (rows with more than max (%g * sqrt (n), 16) entries are\n"
	"    considered \"dense\", and placed last in output permutation)\n",
	alpha)) ;
    }

    if (aggressive)
    {
	PRINTF (("    aggressive absorption:  yes\n")) ;
    }
    else
    {
	PRINTF (("    aggressive absorption:  no\n")) ;
    }

    PRINTF (("    size of AMD integer: %d\n\n", sizeof (Int))) ;
}

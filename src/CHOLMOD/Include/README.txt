CHOLMOD: a sparse Cholesky factorization package.  http://www.suitesparse.com

The Include/*.h files in this directory provide a basic documentation of all
user-callable routines and user-visible data structures in the CHOLMOD
package.  Start with cholmod.h, which describes the general structure of
the parameter lists of CHOLMOD routines.  cholmod_core.h describes the
data structures and basic operations on them (creating and deleting them).

cholmod.h		single include file for all user programs
cholmod_config.h	CHOLMOD compile-time configuration

cholmod_core.h		Core module: data structures and basic support routines
cholmod_check.h		Check module: check/print CHOLMOD data structures
cholmod_cholesky.h	Cholesky module: LL' and LDL' factorization
cholmod_matrixops.h	MatrixOps module: sparse matrix operators (add, mult,..)
cholmod_modify.h	Modify module: update/downdate/...
cholmod_partition.h	Partition module: nested dissection ordering
cholmod_supernodal.h	Supernodal module: supernodal Cholesky

These include files are not used in user programs, but in CHOLMOD only:

cholmod_blas.h		BLAS definitions
cholmod_complexity.h	complex arithmetic
cholmod_template.h	complex arithmetic for template routines
cholmod_internal.h	internal definitions, not visible to user program

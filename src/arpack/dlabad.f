      SUBROUTINE IGRAPHARPACKDLABAD( SMALL, LARGE )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   LARGE, SMALL
*     ..
*
*  Purpose
*  =======
*
*  IGRAPHARPACKDLABAD takes as input the values computed by SLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by IGRAPHARPACKDLAMCH.  This subroutine is needed because
*  IGRAPHARPACKDLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by IGRAPHARPACKDLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by IGRAPHARPACKDLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*
*  =====================================================================
*
*     .. Intrinsic Functions ..
      INTRINSIC          LOG10, SQRT
*     ..
*     .. Executable Statements ..
*
*     If it looks like we're on a Cray, take the square root of
*     SMALL and LARGE to avoid overflow and underflow problems.
*
      IF( LOG10( LARGE ).GT.2000.D0 ) THEN
         SMALL = SQRT( SMALL )
         LARGE = SQRT( LARGE )
      END IF
*
      RETURN
*
*     End of IGRAPHARPACKDLABAD
*
      END

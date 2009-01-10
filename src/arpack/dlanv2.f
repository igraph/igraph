      SUBROUTINE IGRAPHARPACKDLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994 
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
*     ..
*
*  Purpose
*  =======
*
*  IGRAPHARPACKDLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
*  matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
*
*  where either
*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
*  conjugate eigenvalues.
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION
*  B       (input/output) DOUBLE PRECISION
*  C       (input/output) DOUBLE PRECISION
*  D       (input/output) DOUBLE PRECISION
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) DOUBLE PRECISION
*  RT1I    (output) DOUBLE PRECISION
*  RT2R    (output) DOUBLE PRECISION
*  RT2I    (output) DOUBLE PRECISION
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) DOUBLE PRECISION
*          Parameters of the rotation matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   AA, BB, CC, CS1, DD, P, SAB, SAC, SIGMA, SN1,
     $                   TAU, TEMP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   IGRAPHARPACKDLAPY2
      EXTERNAL           IGRAPHARPACKDLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Initialize CS and SN
*
      CS = ONE
      SN = ZERO
*
      IF( C.EQ.ZERO ) THEN
         GO TO 10
*
      ELSE IF( B.EQ.ZERO ) THEN
*
*        Swap rows and columns
*
         CS = ZERO
         SN = ONE
         TEMP = D
         D = A
         A = TEMP
         B = -C
         C = ZERO
         GO TO 10
      ELSE IF( (A-D).EQ.ZERO .AND. SIGN( ONE, B ).NE.
     $   SIGN( ONE, C ) ) THEN
         GO TO 10
      ELSE
*
*        Make diagonal elements equal
*
         TEMP = A - D
         P = HALF*TEMP
         SIGMA = B + C
         TAU = IGRAPHARPACKDLAPY2( SIGMA, TEMP )
         CS1 = SQRT( HALF*( ONE+ABS( SIGMA ) / TAU ) )
         SN1 = -( P / ( TAU*CS1 ) )*SIGN( ONE, SIGMA )
*
*        Compute [ AA  BB ] = [ A  B ] [ CS1 -SN1 ]
*                [ CC  DD ]   [ C  D ] [ SN1  CS1 ]
*
         AA = A*CS1 + B*SN1
         BB = -A*SN1 + B*CS1
         CC = C*CS1 + D*SN1
         DD = -C*SN1 + D*CS1
*
*        Compute [ A  B ] = [ CS1  SN1 ] [ AA  BB ]
*                [ C  D ]   [-SN1  CS1 ] [ CC  DD ]
*
         A = AA*CS1 + CC*SN1
         B = BB*CS1 + DD*SN1
         C = -AA*SN1 + CC*CS1
         D = -BB*SN1 + DD*CS1
*
*        Accumulate transformation
*
         TEMP = CS*CS1 - SN*SN1
         SN = CS*SN1 + SN*CS1
         CS = TEMP
*
         TEMP = HALF*( A+D )
         A = TEMP
         D = TEMP
*
         IF( C.NE.ZERO ) THEN
            IF( B.NE.ZERO ) THEN
               IF( SIGN( ONE, B ).EQ.SIGN( ONE, C ) ) THEN
*
*                 Real eigenvalues: reduce to upper triangular form
*
                  SAB = SQRT( ABS( B ) )
                  SAC = SQRT( ABS( C ) )
                  P = SIGN( SAB*SAC, C )
                  TAU = ONE / SQRT( ABS( B+C ) )
                  A = TEMP + P
                  D = TEMP - P
                  B = B - C
                  C = ZERO
                  CS1 = SAB*TAU
                  SN1 = SAC*TAU
                  TEMP = CS*CS1 - SN*SN1
                  SN = CS*SN1 + SN*CS1
                  CS = TEMP
               END IF
            ELSE
               B = -C
               C = ZERO
               TEMP = CS
               CS = -SN
               SN = TEMP
            END IF
         END IF
      END IF
*
   10 CONTINUE
*
*     Store eigenvalues in (RT1R,RT1I) and (RT2R,RT2I).
*
      RT1R = A
      RT2R = D
      IF( C.EQ.ZERO ) THEN
         RT1I = ZERO
         RT2I = ZERO
      ELSE
         RT1I = SQRT( ABS( B ) )*SQRT( ABS( C ) )
         RT2I = -RT1I
      END IF
      RETURN
*
*     End of IGRAPHARPACKDLANV2
*
      END

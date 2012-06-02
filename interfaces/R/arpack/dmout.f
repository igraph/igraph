*-----------------------------------------------------------------------
*  Routine:    DMOUT
*
*  Purpose:    Real matrix output routine.
*
*  Usage:      CALL DMOUT (LOUT, M, N, A, LDA, IDIGIT, IFMT)
*
*  Arguments
*     M      - Number of rows of A.  (Input)
*     N      - Number of columns of A.  (Input)
*     A      - Real M by N matrix to be printed.  (Input)
*     LDA    - Leading dimension of A exactly as specified in the
*              dimension statement of the calling program.  (Input)
*     IFMT   - Format to be used in printing matrix A.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE IGRAPHDMOUT( LOUT, M, N, A, LDA, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
*     .. Scalar Arguments ..
      CHARACTER*( * )    IFMT
      INTEGER            IDIGIT, LDA, LOUT, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       LINE
      INTEGER            I, J, K1, K2, LLL, NDIGIT
*     ..
*     .. Local Arrays ..
      CHARACTER          ICOL( 3 )
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN, MIN0
*     ..
*     .. Data statements ..
      DATA               ICOL( 1 ), ICOL( 2 ), ICOL( 3 ) / 'C', 'o',
     $                   'l' /
*     ..
*     .. Executable Statements ..
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
c$$$      LLL = MIN( LEN( IFMT ), 80 )
c$$$      DO 10 I = 1, LLL
c$$$         LINE( I: I ) = '-'
c$$$   10 CONTINUE
c$$$*
c$$$      DO 20 I = LLL + 1, 80
c$$$         LINE( I: I ) = ' '
c$$$   20 CONTINUE
c$$$*
c$$$      WRITE( LOUT, FMT = 9999 )IFMT, LINE( 1: LLL )
c$$$ 9999 FORMAT( / 1X, A, / 1X, A )
c$$$*
c$$$      IF( M.LE.0 .OR. N.LE.0 .OR. LDA.LE.0 )
c$$$     $   RETURN
c$$$      NDIGIT = IDIGIT
c$$$      IF( IDIGIT.EQ.0 )
c$$$     $   NDIGIT = 4
c$$$*
c$$$*=======================================================================
c$$$*             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
c$$$*=======================================================================
c$$$*
c$$$      IF( IDIGIT.LT.0 ) THEN
c$$$         NDIGIT = -IDIGIT
c$$$         IF( NDIGIT.LE.4 ) THEN
c$$$            DO 40 K1 = 1, N, 5
c$$$               K2 = MIN0( N, K1+4 )
c$$$               WRITE( LOUT, FMT = 9998 )( ICOL, I, I = K1, K2 )
c$$$               DO 30 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9994 )I, ( A( I, J ), J = K1, K2 )
c$$$   30          CONTINUE
c$$$   40       CONTINUE
c$$$*
c$$$         ELSE IF( NDIGIT.LE.6 ) THEN
c$$$            DO 60 K1 = 1, N, 4
c$$$               K2 = MIN0( N, K1+3 )
c$$$               WRITE( LOUT, FMT = 9997 )( ICOL, I, I = K1, K2 )
c$$$               DO 50 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9993 )I, ( A( I, J ), J = K1, K2 )
c$$$   50          CONTINUE
c$$$   60       CONTINUE
c$$$*
c$$$         ELSE IF( NDIGIT.LE.10 ) THEN
c$$$            DO 80 K1 = 1, N, 3
c$$$               K2 = MIN0( N, K1+2 )
c$$$               WRITE( LOUT, FMT = 9996 )( ICOL, I, I = K1, K2 )
c$$$               DO 70 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9992 )I, ( A( I, J ), J = K1, K2 )
c$$$   70          CONTINUE
c$$$   80       CONTINUE
c$$$*
c$$$         ELSE
c$$$            DO 100 K1 = 1, N, 2
c$$$               K2 = MIN0( N, K1+1 )
c$$$               WRITE( LOUT, FMT = 9995 )( ICOL, I, I = K1, K2 )
c$$$               DO 90 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9991 )I, ( A( I, J ), J = K1, K2 )
c$$$   90          CONTINUE
c$$$  100       CONTINUE
c$$$         END IF
c$$$*
c$$$*=======================================================================
c$$$*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
c$$$*=======================================================================
c$$$*
c$$$      ELSE
c$$$         IF( NDIGIT.LE.4 ) THEN
c$$$            DO 120 K1 = 1, N, 10
c$$$               K2 = MIN0( N, K1+9 )
c$$$               WRITE( LOUT, FMT = 9998 )( ICOL, I, I = K1, K2 )
c$$$               DO 110 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9994 )I, ( A( I, J ), J = K1, K2 )
c$$$  110          CONTINUE
c$$$  120       CONTINUE
c$$$*
c$$$         ELSE IF( NDIGIT.LE.6 ) THEN
c$$$            DO 140 K1 = 1, N, 8
c$$$               K2 = MIN0( N, K1+7 )
c$$$               WRITE( LOUT, FMT = 9997 )( ICOL, I, I = K1, K2 )
c$$$               DO 130 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9993 )I, ( A( I, J ), J = K1, K2 )
c$$$  130          CONTINUE
c$$$  140       CONTINUE
c$$$*
c$$$         ELSE IF( NDIGIT.LE.10 ) THEN
c$$$            DO 160 K1 = 1, N, 6
c$$$               K2 = MIN0( N, K1+5 )
c$$$               WRITE( LOUT, FMT = 9996 )( ICOL, I, I = K1, K2 )
c$$$               DO 150 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9992 )I, ( A( I, J ), J = K1, K2 )
c$$$  150          CONTINUE
c$$$  160       CONTINUE
c$$$*
c$$$         ELSE
c$$$            DO 180 K1 = 1, N, 5
c$$$               K2 = MIN0( N, K1+4 )
c$$$               WRITE( LOUT, FMT = 9995 )( ICOL, I, I = K1, K2 )
c$$$               DO 170 I = 1, M
c$$$                  WRITE( LOUT, FMT = 9991 )I, ( A( I, J ), J = K1, K2 )
c$$$  170          CONTINUE
c$$$  180       CONTINUE
c$$$         END IF
c$$$      END IF
c$$$      WRITE( LOUT, FMT = 9990 )
c$$$*
c$$$ 9998 FORMAT( 10X, 10( 4X, 3A1, I4, 1X ) )
c$$$ 9997 FORMAT( 10X, 8( 5X, 3A1, I4, 2X ) )
c$$$ 9996 FORMAT( 10X, 6( 7X, 3A1, I4, 4X ) )
c$$$ 9995 FORMAT( 10X, 5( 9X, 3A1, I4, 6X ) )
c$$$ 9994 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 10D12.3 )
c$$$ 9993 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 8D14.5 )
c$$$ 9992 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 6D18.9 )
c$$$ 9991 FORMAT( 1X, ' Row', I4, ':', 1X, 1P, 5D22.13 )
c$$$ 9990 FORMAT( 1X, ' ' )
*
      RETURN
      END

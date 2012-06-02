*-----------------------------------------------------------------------
*  Routine:    DVOUT
*
*  Purpose:    Real vector output routine.
*
*  Usage:      CALL DVOUT (LOUT, N, SX, IDIGIT, IFMT)
*
*  Arguments
*     N      - Length of array SX.  (Input)
*     SX     - Real array to be printed.  (Input)
*     IFMT   - Format to be used in printing array SX.  (Input)
*     IDIGIT - Print up to IABS(IDIGIT) decimal digits per number.  (In)
*              If IDIGIT .LT. 0, printing is done with 72 columns.
*              If IDIGIT .GT. 0, printing is done with 132 columns.
*
*-----------------------------------------------------------------------
*
      SUBROUTINE IGRAPHDVOUT( LOUT, N, SX, IDIGIT, IFMT )
*     ...
*     ... SPECIFICATIONS FOR ARGUMENTS
*     ...
*     ... SPECIFICATIONS FOR LOCAL VARIABLES
*     .. Scalar Arguments ..
      CHARACTER*( * )    IFMT
      INTEGER            IDIGIT, LOUT, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   SX( * )
*     ..
*     .. Local Scalars ..
      CHARACTER*80       LINE
      INTEGER            I, K1, K2, LLL, NDIGIT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN, MIN, MIN0
*     ..
*     .. Executable Statements ..
*     ...
*     ... FIRST EXECUTABLE STATEMENT
*
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
c$$$      IF( N.LE.0 )
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
c$$$            DO 30 K1 = 1, N, 5
c$$$               K2 = MIN0( N, K1+4 )
c$$$               WRITE( LOUT, FMT = 9998 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   30       CONTINUE
c$$$         ELSE IF( NDIGIT.LE.6 ) THEN
c$$$            DO 40 K1 = 1, N, 4
c$$$               K2 = MIN0( N, K1+3 )
c$$$               WRITE( LOUT, FMT = 9997 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   40       CONTINUE
c$$$         ELSE IF( NDIGIT.LE.10 ) THEN
c$$$            DO 50 K1 = 1, N, 3
c$$$               K2 = MIN0( N, K1+2 )
c$$$               WRITE( LOUT, FMT = 9996 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   50       CONTINUE
c$$$         ELSE
c$$$            DO 60 K1 = 1, N, 2
c$$$               K2 = MIN0( N, K1+1 )
c$$$               WRITE( LOUT, FMT = 9995 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   60       CONTINUE
c$$$         END IF
c$$$*
c$$$*=======================================================================
c$$$*             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
c$$$*=======================================================================
c$$$*
c$$$      ELSE
c$$$         IF( NDIGIT.LE.4 ) THEN
c$$$            DO 70 K1 = 1, N, 10
c$$$               K2 = MIN0( N, K1+9 )
c$$$               WRITE( LOUT, FMT = 9998 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   70       CONTINUE
c$$$         ELSE IF( NDIGIT.LE.6 ) THEN
c$$$            DO 80 K1 = 1, N, 8
c$$$               K2 = MIN0( N, K1+7 )
c$$$               WRITE( LOUT, FMT = 9997 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   80       CONTINUE
c$$$         ELSE IF( NDIGIT.LE.10 ) THEN
c$$$            DO 90 K1 = 1, N, 6
c$$$               K2 = MIN0( N, K1+5 )
c$$$               WRITE( LOUT, FMT = 9996 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$   90       CONTINUE
c$$$         ELSE
c$$$            DO 100 K1 = 1, N, 5
c$$$               K2 = MIN0( N, K1+4 )
c$$$               WRITE( LOUT, FMT = 9995 )K1, K2, ( SX( I ), I = K1, K2 )
c$$$  100       CONTINUE
c$$$         END IF
c$$$      END IF
c$$$      WRITE( LOUT, FMT = 9994 )
c$$$      RETURN
c$$$ 9998 FORMAT( 1X, I4, ' - ', I4, ':', 1P, 10D12.3 )
c$$$ 9997 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 8D14.5 )
c$$$ 9996 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 6D18.9 )
c$$$ 9995 FORMAT( 1X, I4, ' - ', I4, ':', 1X, 1P, 5D24.13 )
c$$$ 9994 FORMAT( 1X, ' ' )
      END

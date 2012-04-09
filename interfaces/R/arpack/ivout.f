C-----------------------------------------------------------------------
C  Routine:    IVOUT
C
C  Purpose:    Integer vector output routine.
C
C  Usage:      CALL IVOUT (LOUT, N, IX, IDIGIT, IFMT)
C
C  Arguments
C     N      - Length of array IX. (Input)
C     IX     - Integer array to be printed. (Input)
C     IFMT   - Format to be used in printing array IX. (Input)
C     IDIGIT - Print up to ABS(IDIGIT) decimal digits / number. (Input)
C              If IDIGIT .LT. 0, printing is done with 72 columns.
C              If IDIGIT .GT. 0, printing is done with 132 columns.
C
C-----------------------------------------------------------------------
C
      SUBROUTINE IGRAPHIVOUT (LOUT, N, IX, IDIGIT, IFMT)
C     ...
C     ... SPECIFICATIONS FOR ARGUMENTS
      INTEGER    IX(*), N, IDIGIT, LOUT
      CHARACTER  IFMT*(*)
C     ...
C     ... SPECIFICATIONS FOR LOCAL VARIABLES
      INTEGER    I, NDIGIT, K1, K2, LLL
      CHARACTER*80 LINE
*     ...
*     ... SPECIFICATIONS INTRINSICS
      INTRINSIC          MIN
*
C
c$$$      LLL = MIN ( LEN ( IFMT ), 80 )
c$$$      DO 1 I = 1, LLL
c$$$          LINE(I:I) = '-'
c$$$    1 CONTINUE
c$$$C
c$$$      DO 2 I = LLL+1, 80
c$$$          LINE(I:I) = ' '
c$$$    2 CONTINUE
c$$$C
c$$$      WRITE ( LOUT, 2000 ) IFMT, LINE(1:LLL)
c$$$ 2000 FORMAT ( /1X, A  /1X, A )
c$$$C
c$$$      IF (N .LE. 0) RETURN
c$$$      NDIGIT = IDIGIT
c$$$      IF (IDIGIT .EQ. 0) NDIGIT = 4
c$$$C
c$$$C=======================================================================
c$$$C             CODE FOR OUTPUT USING 72 COLUMNS FORMAT
c$$$C=======================================================================
c$$$C
c$$$      IF (IDIGIT .LT. 0) THEN
c$$$C
c$$$      NDIGIT = -IDIGIT
c$$$      IF (NDIGIT .LE. 4) THEN
c$$$         DO 10 K1 = 1, N, 10
c$$$            K2 = MIN0(N,K1+9)
c$$$            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
c$$$   10    CONTINUE
c$$$C
c$$$      ELSE IF (NDIGIT .LE. 6) THEN
c$$$         DO 30 K1 = 1, N, 7
c$$$            K2 = MIN0(N,K1+6)
c$$$            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
c$$$   30    CONTINUE
c$$$C
c$$$      ELSE IF (NDIGIT .LE. 10) THEN
c$$$         DO 50 K1 = 1, N, 5
c$$$            K2 = MIN0(N,K1+4)
c$$$            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
c$$$   50    CONTINUE
c$$$C
c$$$      ELSE
c$$$         DO 70 K1 = 1, N, 3
c$$$            K2 = MIN0(N,K1+2)
c$$$            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
c$$$   70    CONTINUE
c$$$      END IF
c$$$C
c$$$C=======================================================================
c$$$C             CODE FOR OUTPUT USING 132 COLUMNS FORMAT
c$$$C=======================================================================
c$$$C
c$$$      ELSE
c$$$C
c$$$      IF (NDIGIT .LE. 4) THEN
c$$$         DO 90 K1 = 1, N, 20
c$$$            K2 = MIN0(N,K1+19)
c$$$            WRITE(LOUT,1000) K1,K2,(IX(I),I=K1,K2)
c$$$   90    CONTINUE
c$$$C
c$$$      ELSE IF (NDIGIT .LE. 6) THEN
c$$$         DO 110 K1 = 1, N, 15
c$$$            K2 = MIN0(N,K1+14)
c$$$            WRITE(LOUT,1001) K1,K2,(IX(I),I=K1,K2)
c$$$  110    CONTINUE
c$$$C
c$$$      ELSE IF (NDIGIT .LE. 10) THEN
c$$$         DO 130 K1 = 1, N, 10
c$$$            K2 = MIN0(N,K1+9)
c$$$            WRITE(LOUT,1002) K1,K2,(IX(I),I=K1,K2)
c$$$  130    CONTINUE
c$$$C
c$$$      ELSE
c$$$         DO 150 K1 = 1, N, 7
c$$$            K2 = MIN0(N,K1+6)
c$$$            WRITE(LOUT,1003) K1,K2,(IX(I),I=K1,K2)
c$$$  150    CONTINUE
c$$$      END IF
c$$$      END IF
c$$$      WRITE (LOUT,1004)
c$$$C
c$$$ 1000 FORMAT(1X,I4,' - ',I4,':',20(1X,I5))
c$$$ 1001 FORMAT(1X,I4,' - ',I4,':',15(1X,I7))
c$$$ 1002 FORMAT(1X,I4,' - ',I4,':',10(1X,I11))
c$$$ 1003 FORMAT(1X,I4,' - ',I4,':',7(1X,I15))
c$$$ 1004 FORMAT(1X,' ')
c$$$C
      RETURN
      END

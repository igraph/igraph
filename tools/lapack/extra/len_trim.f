*
*  -- LEN_TRIM is Fortran 95, so we use a replacement here
*
      FUNCTION LEN_TRIM(S)
*
      CHARACTER*(*)  S
      INTEGER LEN_TRIM
*
      INTRINSIC LEN
*
      DO LEN_TRIM = LEN(S), 1, -1
         IF (s(LEN_TRIM:LEN_TRIM) .NE. ' ') RETURN
      END DO
      END

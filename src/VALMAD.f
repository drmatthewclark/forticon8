C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      REAL*8 FUNCTION VALMAD(A,B,R)
C
C     FUNCTION ROUTINE FOR CALCULATING MADELUNG PARAMETERS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      IF(A.LT.0.01D0.OR.B.LT.0.01D0) GO TO 1
      AB=(A+B)/(2.0D0*A*B)
      VALMAD=1.0D0/DSQRT(R*R+AB*AB)
      GO TO 2
1     VALMAD=0.0D0
      IF(R.GT.0.001D0) VALMAD=1.0D0/R
2     RETURN
      END

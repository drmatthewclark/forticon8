C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      REAL*8 FUNCTION DSUM(B,A,IP1,LIMIT)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION B(*),A(*)
C
C    FUNCTION FOR USE IN GIVENS DIAGNOLIZATION PACKAGE.
C
      JJ=1
      DSUM=0.D0
      DO 180 II=IP1,LIMIT
      DSUM=DSUM+B(II+1)*A(JJ)
 180  JJ=JJ+II
      RETURN
      END

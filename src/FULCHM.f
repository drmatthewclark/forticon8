C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE FULCHM(E,U,C,H,N5)
      IMPLICIT REAL*8 (A-H,O-Z)
      REAL*8 E(N5),C(N5*N5),U(N5,N5),H(N5,N5)
      LOGICAL*1 JGO
C
C    SUBROUTINE TO CALCULATE COMPLETE CHARGE MATRIX.
C
      KJ=1
      DO 31 I=1,N5
      JGO=.FALSE.
      IJ=KJ
      DO 32 K=1,N5
      IF(JGO) GO TO 34
      IF(I.NE.K) GO TO 35
      JGO=.TRUE.
      E(K)=1.0D0
      GO TO 32
 35   E(K)=C(IJ)
      IJ=IJ+1
      GO TO 32
 34   IJ=IJ+K-2
      E(K)=C(IJ)
 32   CONTINUE
      KJ=KJ+I-1
      DO 31 J=1,N5
      UB=0.0D0
      DO 36 K=1,N5
 36   UB=UB+H(K,J)*E(K)
 31   U(I,J)=2.0D0*H(I,J)*UB
      RETURN
      END

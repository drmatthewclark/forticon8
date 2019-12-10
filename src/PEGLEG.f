C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE PEGLEG(A,N,NL)
C
C  SUBROUTINE TO PRINT OUT MATRICES IN READABLE FORMAT.
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NL,NL)
      NROW=N
      NCOL=N
      GO TO 10
      ENTRY OUTMAT(A,NL,NR,NC)
      NROW=NR
      NCOL=NC
10    KITE=0
20    LOW=KITE+1
      KITE=KITE+14
      IF(KITE.GT.NCOL) KITE=NCOL
      WRITE(6,1000) (I,I=LOW,KITE)
1000  FORMAT(/5X,14I8,//)
      DO 30 I=1,NROW
30    WRITE(6,1001) I,(A(I,J),J=LOW,KITE)
1001  FORMAT(I5,2X,14F8.4)
      IF(KITE.LT.NCOL) GO TO 20
      RETURN
      END

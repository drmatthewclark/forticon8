C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE PEGLEG3(AA,NN,IUNIT)
C
C     ROUTINE FOR WRITING OUT THE ORBITAL OCCUPANCIES TO DISK FILE 13 FOR
C     PLOTTING ROUTINES.   JJN  8-8-90
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION AA(NN)

      WRITE(IUNIT,9998) (AA(J), J=1,NN)
9998  FORMAT(8F10.6)
      CLOSE(IUNIT)

      RETURN
      END

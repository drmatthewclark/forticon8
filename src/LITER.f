C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE LITER(NH,NA,E,W,NCYC)
C
C  SUBROUTINE FOR CALCULATING Q*SENSE WHEN USING CHARGE ITERATION
C  OPTION ( METH = 1 ).
C
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARAMETERS'
      INCLUDE 'ATOMCOMMON'

      DIMENSION E(1),W(1)
 
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*1 PARTIT,PRINTX,ITABLE
      NCYC=ICYCLE+1-(ICYCLE/4)*4
      ICYCLE=ICYCLE+1
      DELTA=0.D0
      NATOM=NH+NA
      INDEX=1
      DO 60 I=1,NATOM
      INCR=1
      IF(I.GT.NH) GO TO 100
      CHG=1.0D0-E(INDEX)
      GO TO 150
100   KEYI=KEY(I-NH)
      CHG=E(INDEX)
      IF(ND(KEYI).EQ.0) GO TO 120
      INCR=9
      CHG=CHG+E(INDEX+4)+E(INDEX+5)+E(INDEX+6)+E(INDEX+7)+E(INDEX+8)
      GO TO 130
120   IF(NP(KEYI).EQ.0) GO TO 140
      INCR=4
130   CHG=CHG+E(INDEX+1)+E(INDEX+2)+E(INDEX+3)
140   CHG=VELEC(KEYI)-CHG
150   INDEX=INDEX+INCR
      GO TO (10,20,30,40),NCYC
10    CHG=0.25D0*(CHG+Y(I)+Z(I)+W(I))
      X(I)=CHG
      Y(I)=CHG
      Z(I)=CHG
      W(I)=CHG
      GO TO 60
20    DELTA=DELTA+DABS(CHG-X(I))
      Y(I)=CHG
      GO TO 60
30    DELTA=DELTA+DABS(CHG-Y(I))
      Z(I)=CHG
      GO TO 60
40    DELTA=DELTA+DABS(CHG-Z(I))
      W(I)=CHG
60    E(I)=CHG*SENSE
      IF(DELTAC.LT.DELTA.OR.NCYC.EQ.1) RETURN
      ICYCLE=15000
      RETURN
      END
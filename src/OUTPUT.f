C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
C     NDIM - NUMBER OF ELECTRONS
C
      SUBROUTINE OUTPUT(H,U,MAD,C,E,W,IOCC,HDG,NDIM,NTYPE,NC,NHDG)
C
C  SUBROUTINE TO ANALYSE AND PRINT OUT RESULTS.
C
      IMPLICIT REAL*8(A-H,O-Z)
      INCLUDE 'PARAMETERS'

      DIMENSION H(NDIM,NDIM),U(NDIM,NDIM),MAD(NTYPE,NTYPE),C(NDIM*NDIM),
     1E(NDIM),W(NDIM),IOCC(NDIM),HDG(NHDG)

      REAL*8 MAD
      REAL*4 IOCC

      COMMON/TITLE/AB(10)
      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,IPRINT,
     1IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      LOGICAL*1 PRT,PUN
      INTEGER*2 IOVPOP,IENRGY
 
      INCLUDE 'ATOMCOMMON'

      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,
     1ICYCLE,NCON,PARTIT,PRINTX,ITABLE(20)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*1 PARTIT,PRINTX,ITABLE
      CHARACTER*2 SYMB,HYDROG
      DATA HYDROG/' H'/
      NMIN=0

      DO 5 I=1,NDIM
      NMIN=I
      IF(IOCC(I).GT.0.0001) GO TO 7
5     CONTINUE
7     IF(ITERAT) GO TO 38
C
C  CALCULATE AND PRINT OUT SUM OF ONE-ELECTRON ENERGIES.
C
      SUM=0.0D0
      DO 13 I=1,NDIM
      W(I)=DBLE(IOCC(I))
13    SUM=SUM+E(I)*W(I)
      IF(.NOT.PRT(9)) WRITE(6,2001) SUM
2001  FORMAT(T10,'SUM OF ONE-ELECTRON ENERGIES =',F16.8,' EV.',///)
      IF(PUN(9)) WRITE(7,2002) SUM
2002  FORMAT(F20.8)
      IF(ONEMAT) GO TO 9999
C
C  PRINT OUT WAVE FUNCTIONS.
C
      IF(PRT(10)) GO TO 1003
      WRITE(6,1002)
1002  FORMAT('WAVE FUNCTIONS'/'MO''S IN COLUMNS, AO''S IN ROWS')
      CALL PEGLEG(H,NDIM,NDIM)
C
C  ** THIS CALL IS TO A SIMILAR ROUTINE TO WRITE THE MO COEFFICIENTS
C     TO DISK FILE 13. THE FORMAT IS DIFFERENT.
C     FIRST NDIM IS MO
C     SECOND NDIM is AO

      CALL PEGLEG2(H,NDIM,NDIM,13)
C
C  **
C
1003  IF(PUN(10)) WRITE(7,2003) H
2003  FORMAT(8F9.6)
C
C  CALCULATE AND PRINT OUT DENSITY MATRIX.
C
      IF(PRT(11).AND..NOT.PUN(11)) GO TO 38
      DO 300 I=1,NDIM
      DO 300 J=1,I
      U(I,J)=0.0D0
      DO 310 K=1,NDIM
310   U(I,J)=U(I,J)+H(I,K)*H(J,K)*W(K)
300   U(J,I)=U(I,J)
      IF(PRT(11)) GO TO 360
      WRITE(6,350)
350   FORMAT(/,'DENSITY MATRIX')
      CALL PEGLEG(U,NDIM,NDIM)
360   IF(PUN(11)) WRITE(7,2003) U
C
C  CALCULATE ATOMIC ORBITAL OCCUPATIONS AND STORE IN E(I).
C  CALCULATE OVERLAP POPULATION MATRIX.
C
38    IJ=1
      DO 60 I=1,NDIM
      E(I)=0.0D0
      DO 60 J=1,I
      UB=0.0D0
      DO 41 K=NMIN,NDIM
41    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))
      UB=UB*0.5D0
      IF(I.EQ.J) GO TO 50
      UB=(UB+UB)*C(IJ)
      IJ=IJ+1
50    E(I)=E(I)+UB
      E(J)=E(J)+UB
      IF(ITERAT) GO TO 60
      UB=UB+UB
      U(I,J)=UB
      U(J,I)=UB
60    CONTINUE
C
C  IF DOING CHARGE ITERATION ( METH=1 ) CALL LITER.
C
      IF(.NOT.ITERAT) GO TO 80
      K=ICYCLE
      PRINTX=((ICYCLE/PRTCYC)*PRTCYC.EQ.ICYCLE)
      IF(ICYCLE.EQ.MAXCYC) PRINTX=.TRUE.
      CALL LITER(NH,NA,E,W,J)
      IF(ICYCLE.EQ.15000) PRINTX=.TRUE.
      IF(K.EQ.1) WRITE(6,355)
355   FORMAT('      ATOMIC CHARGES',//)
      GO TO (81,82,83,84),J
81    WRITE(6,375) K,(X(I),I=1,NATM)
375   FORMAT(/,T3,'CYCLE NO.',I3,(T20,10F10.5))
      GO TO 9000
82    WRITE(6,375) K,(Y(I),I=1,NATM)
      GO TO 9000
83    WRITE(6,375) K,(Z(I),I=1,NATM)
      GO TO 9000
84    WRITE(6,375) K,(W(I),I=1,NATM)
9000  RETURN
C
C  PRINT OUT OVERLAP POPULATION MATRIX.
C
80    IF(PRT(12)) GO TO 1009
      WRITE(6,1006) NELEC
1006  FORMAT(//,'OVERLAP POPULATION MATRIX FOR',I4,' ELECTRONS')
      CALL PEGLEG(U,NDIM,NDIM)
1009  IF(PUN(12)) WRITE(7,2003) U
C
C  CALL REDUCE TO CALCULATE REDUCED OVERLAP MATRIX. PRINT OUT
C  REDUCED OVERLAP MATRIX.
C
      IF(PRT(13).AND..NOT.PUN(13)) GO TO 2005
      CALL REDUCE(U,NDIM,NA,NH)
      DO 100 I=2,NATM
      K=I-1
      DO 100 J=1,K
100   U(J,I)=U(I,J)
      IF(PRT(13)) GO TO 1015
      WRITE(6,1007)
1007  FORMAT(/,'REDUCED OVERLAP POPULATION MATRIX, ATOM BY ATOM')
      CALL PEGLEG(U,NATM,NDIM)
1015  IF(PUN(13)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NATM)
C
C  IF L3 IS TRUE, CALCULATE AND PRINT OUT OVERLAP POPULATION ANALYSIS,
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.
C
2005  IF(.NOT.L3) GO TO 21
      DO 600 N=1,23,2
      IF(IOVPOP(N).EQ.0) GO TO 25
      KMIN=IOVPOP(N)
      KMAX=IOVPOP(N+1)
      DO 600 K=KMIN,KMAX
      WRITE(6,2004) K,W(K)
2004  FORMAT(///'OVERLAP POPULATION MATRIX, ORBITAL BY ORBITAL, FOR MOL
     1ECULAR ORBITAL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')
      IF(IOCC(K).LT.0.0001) GO TO 600
      SUM=0.5D0*W(K)
      IJ=1
      DO 460 I=1,NDIM
      DO 460 J=1,I
      UB=H(I,K)*H(J,K)
      IF(I.EQ.J) GO TO 450
      UB=(UB+UB)*C(IJ)
      IJ=IJ+1
450   UB=(UB+UB)*SUM
      U(J,I)=UB
      U(I,J)=UB
460   CONTINUE
      CALL PEGLEG(U,NDIM,NDIM)
600   CONTINUE
25    WRITE(6,7003)
7003  FORMAT(///)
C
C  CALL FULCHM TO CALCULATE COMPLETE CHARGE MATRIX. PRINT OUT
C  COMPLETE CHARGE MATRIX.
C
21    L1=PRT(14).AND..NOT.PUN(14)
      L2=PRT(15).AND..NOT.PUN(15)
      IF(L1.AND.L2) GO TO 1020
      CALL FULCHM(W,U,C,H,NDIM)
      IF(PRT(14)) GO TO 2021
      WRITE(6,1008)
1008  FORMAT('COMPLETE CHARGE MATRIX FOR EACH MO, NORMALIZED TO TWO ELE
     1CTRONS REGARDLESS OF OCCUPATION')
      CALL PEGLEG(U,NDIM,NDIM)
2021  IF(PUN(14)) WRITE(7,2003) U
C
C  CALL REDCHM TO CALCULATE REDUCED CHARGE MATRIX. PRINT OUT
C  REDUCED CHARGE MATRIX.
C
      IF(L2) GO TO 1020
      CALL REDCHM(U,NDIM)
      IF(PRT(15)) GO TO 1022
      WRITE(6,1019)
1019  FORMAT(/,'REDUCED CHARGE MATRIX, MO''S IN COLUMNS, ATOMS IN ROWS')
      CALL OUTMAT(U,NDIM,NATM,NDIM)
1022  IF(PUN(15)) WRITE(7,2003) ((U(I,J),I=1,NATM),J=1,NDIM)
C
C  PRINT OUT ATOMIC CHARGES AND ORBITAL OCCUPATIONS.
C
1020  IF(PRT(16)) GO TO 40
      WRITE(6,1010)
1010  FORMAT(/'ATOM',T12,'NET CHG.',T35,'ATOMIC ORBITAL OCCUPATION FOR
     *GIVEN MO OCCUPATION'/T35,'S',T45,'X',T55,'Y',T65,'Z',T75,'X2-Y2',
     *T85,'Z2',T95,'XY',T105,'XZ',T115,'YZ'/)
      SYMB=HYDROG
      J=1
      DO 140 I=1,NATM
      IF(I.GT.NH) GO TO 120
      UB=1.0-E(I)
      N=I
      GO TO 130
120   KITE=I-NH
      KEYI=KEY(KITE)
      SYMB=SYMBOL(KEYI)
      UB=VELEC(KEYI)-E(J)
      N=J
      IF(NP(KEYI).EQ.0) GO TO 130
      UB=UB-E(J+1)-E(J+2)-E(J+3)
      N=J+3
      IF(ND(KEYI).EQ.0) GO TO 130
      UB=UB-E(J+4)-E(J+5)-E(J+6)-E(J+7)-E(J+8)
      N=J+8
130   WRITE(6,1011) SYMB,I,UB,(E(K),K=J,N)
1011  FORMAT(1X,A2,I3,T10,F10.5,T30,9F10.5)
140   J=N+1
C
C     WRITE MO OCCUPANCIES.  E IS THE TOTAL COEFFICIENT
C     FOR EACH ATOMIC ORBITAL. THIS IS FOR WRITING THE
C     TOTAL DENSITY
C
C      CALL PEGLEG3(E,NDIM,13)
C
C  IF CALCULATING ENERGY MATRIX, PUT DIAGONAL ELEMENTS OF HUCKEL
C  MATRIX IN W(I).
C
40    PRT(1)=PRT(17).AND..NOT.PUN(17)
      PRT(2)=PRT(18).AND..NOT.PUN(18)
      PRT(3)=PRT(19).AND..NOT.PUN(19)
      PRT(4)=PRT(20).AND..NOT.PUN(20)
      L1=PRT(1).AND.PRT(2)
      L2=PRT(3).AND.PRT(4)
      IF(L1.AND.L2.AND..NOT.L4) GO TO 9999
      IF(NH.EQ.0) GO TO 23
      DO 22 I=1,NH
22    W(I)=X(I)
23    J=NH+1
      K=NH+1
      DO 24 I=1,NA
      KEYI=KEY(I)
      W(J)=X(K)
      J=J+1
      IF(NP(KEYI).EQ.0) GO TO 24
      UB=Y(K)
      W(J)=UB
      W(J+1)=UB
      W(J+2)=UB
      J=J+3
      IF(ND(KEYI).EQ.0) GO TO 24
      UB=Z(K)
      W(J)=UB
      W(J+1)=UB
      W(J+2)=UB
      W(J+3)=UB
      W(J+4)=UB
      J=J+5
24    K=K+1
C
C  IF DOING CHARGE ITERATION WITH MADELUNG CORRECTION ON MOLECULE
C  WITH NON-ZERO CHARGE, PUT MADELUNG TERMS IN E(I).
C
      QON= DFLOAT(KA)/DFLOAT(NELEC)
      IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 180
      IF(NH.EQ.0) GO TO 181
      DO 182 I=1,NH
182   E(I)=MAD(I,I)
181   J=NH+1
      K=NH+1
      DO 183 I=1,NA
      KEYI=KEY(I)
      E(J)=MAD(K,K)
      J=J+1
      K=K+1
      IF(NP(KEYI).EQ.0) GO TO 183
      UB=MAD(K,K)
      E(J)=UB
      E(J+1)=UB
      E(J+2)=UB
      J=J+3
      K=K+1
      IF(ND(KEYI).EQ.0) GO TO 183
      UB=MAD(K,K)
      E(J)=UB
      E(J+1)=UB
      E(J+2)=UB
      E(J+3)=UB
      E(J+4)=UB
      J=J+5
      K=K+1
183   CONTINUE
C
C  CALCULATE AND PRINT OUT ENERGY MATRIX.
C
180   ONEMAT=.TRUE.
      IF(L1) GO TO 7000
170   SUM=1.0D0
      IF(ONEMAT) SUM=0.0D0
      IJ=1
      CNST=2.0D0*CON
      CN2=CNST
      DO 28 I=1,NDIM
      U(I,I)=0.0D0
      DO 28 J=1,I
      UB=0.0D0
      DO 26 K=NMIN,NDIM
26    UB=UB+H(I,K)*H(J,K)*DBLE(IOCC(K))
      IF(I.EQ.J) GO TO 28
      IF(.NOT.L5) GO TO 35
      UC=W(I)
      ET=W(J)
      IF(NHDG.EQ.1) GO TO 36
      UC=HDG(I)
      ET=HDG(J)
36    UC=(UC-ET)/(UC+ET)
      UC=UC*UC
      CNST=CN2+UC+UC*UC*(1.0D0-CN2)
35    IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 31
      UC=UB*C(IJ)*((CNST-SUM)*(W(I)+W(J))-QON*(1.0D0-CNST)*(E(I)+E(J)))
      GO TO 32
31    UC=UB*(CNST-SUM)*C(IJ)*(W(I)+W(J))
32    UB=SUM*UB*C(IJ)
      IJ=IJ+1
      U(J,I)=UC
      U(I,J)=UC
      U(J,J)=U(J,J)+UB*W(J)
28    U(I,I)=U(I,I)+UB*W(I)
      IF(PRT(17)) GO TO 3020
      IF(ONEMAT) WRITE(6,3005)
3005  FORMAT(///,'ENERGY MATRIX')
      IF(.NOT.ONEMAT) WRITE(6,3006)
3006  FORMAT('ENERGY PARTITIONING')
      CALL PEGLEG(U,NDIM,NDIM)
3020  IF(PUN(17)) WRITE(7,2800) U
2800  FORMAT(8F9.5)
C
C  CALL REDUCE TO CALCULATE REDUCED ENERGY MATRIX. PRINT OUT
C  REDUCED ENERGY MATRIX.
C
      IF(PRT(2)) GO TO 7000
      CALL REDUCE(U,NDIM,NA,NH)
      DO 700 I=2,NATM
      K=I-1
      DO 700 J=1,K
700   U(J,I)=U(I,J)
      IF(PRT(18)) GO TO 3021
      IF(ONEMAT) WRITE(6,3007)
3007  FORMAT(/,'REDUCED ENERGY MATRIX, ATOM BY ATOM')
      IF(.NOT.ONEMAT) WRITE(6,3008)
3008  FORMAT(/,'REDUCED ENERGY PARTITIONING, ATOM BY ATOM')
      KITE=0
701   LOW=KITE+1
      KITE=KITE+13
      IF(KITE.GT.NATM) KITE=NATM
      WRITE(6,702) (I,I=LOW,KITE)
702   FORMAT(/5X,13I9,//)
      DO 703 I=1,NATM
703   WRITE(6,704) I,(U(I,J),J=LOW,KITE)
704   FORMAT(I5,2X,13F9.4)
      IF(KITE.LT.NATM) GO TO 701
3021  IF(PUN(18)) WRITE(7,2111) ((U(I,J),I=1,NATM),J=1,NATM)
2111  FORMAT(7F10.5)
C
C  IF L4 IS TRUE, CALCULATE AND PRINT OUT ENERGY MATRIX ANALYSIS,
C  ORBITAL BY ORBITAL, FOR EACH MOLECULAR ORBITAL SPECIFIED.
C
7000  IF(.NOT.ONEMAT) GO TO 9999
      IF(.NOT.L4) GO TO 71
      DO 246 N=1,23,2
      IF(IENRGY(N).EQ.0) GO TO 72
      KMIN=IENRGY(N)
      KMAX=IENRGY(N+1)
      DO 246 K=KMIN,KMAX
      WRITE(6,3004) K,IOCC(K)
3004  FORMAT(///'ENERGY MATRIX, ORBITAL BY ORBITAL, FOR MOLECULAR ORBIT
     1AL',I4,5X,'OCCUPATION IS',F7.4,' ELECTRONS')
      IF(IOCC(K).LT.0.0001) GO TO 246
      IJ=1
      CNST=CON
      EX=DBLE(IOCC(K))
      DO 244 I=1,NDIM
      DO 244 J=1,I
      UB=H(I,K)*H(J,K)
      IF(I.EQ.J) GO TO 242
      IF(.NOT.L5) GO TO 236
      UC=W(I)
      ET=W(J)
      IF(NHDG.EQ.1) GO TO 237
      UC=HDG(I)
      ET=HDG(J)
237   UC=(UC-ET)/(UC+ET)
      UC=UC*UC
      CNST=CON+UC/2.0D0+UC*UC*(0.5D0-CON)
236   IF(METH.LT.3.OR.DABS(QON).LT.0.0001D0) GO TO 247
      UB=UB*C(IJ)*(CNST*(W(I)+W(J))-QON*(0.5D0-CNST)*(E(I)+E(J)))
      GO TO 248
247   UB=UB*CNST*C(IJ)*(W(I)+W(J))
248   IJ=IJ+1
      UB=2.0D0*UB*EX
      U(I,J)=UB
      U(J,I)=UB
      GO TO 244
242   U(I,I)=UB*W(I)*EX
244   CONTINUE
      CALL PEGLEG(U,NDIM,NDIM)
246   CONTINUE
72    WRITE(6,7003)
C
C  CALCULATE AND PRINT OUT ENERGY PARTITIONING AND REDUCED
C  ENERGY PARTITIONING.
C
71    IF(L2) GO TO 9999
      ONEMAT=.FALSE.
      PRT(17)=PRT(19)
      PUN(17)=PUN(19)
      PRT(18)=PRT(20)
      PUN(18)=PUN(20)
      PRT(2)=PRT(4)
      GO TO 170
9999  RETURN
      END

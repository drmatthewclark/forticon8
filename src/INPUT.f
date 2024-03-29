C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE INPUT(NATOM,NDIM,NTYPE)
C
C     INPUT FORMAT
C
C     TITLE  FORMAT(8A8,A6,A2)
C       IF COL 71 IS '*' READ TITLE CONTINUATION
C     PARAMETERS FORMAT(6I3,5L1,F5.2,2F6.3,40L1)
C        NH,NA,KA,METH,IPRINT,IPUNCH,L1,L2,L3,L4,L5,CON, PEEP,COULH,(PRT(I),I=1,20),(PUN(J),J=1,20)

C        NH   NUMBER OF HYDROGENS
C        NA   NUMBER OF NON-HYDROGENS
C        KA    CHARGE
C        METH   0 - EXTENDED HUCKEL
C               1 - EXTENDED HUCKEL WITH CHARGE ITERATION LINEAR
C                   CHARGE DEPENDENCE FOR SENSE*CHARGE
C               2 - EXTENDED HUCKEL WITH CHARGE INTERATION
C               3 - 2 WITH MADELING CORRECTION
C
C        IPRINT   -2 SET PRT 6,7,11,17,19,12,14,20,13,15,16,18,10 TRUE
C                 -1 SET PRT 13,15,16,18 TRUE
C                  0 SET PRT 6,7,11,17,19,12,14,20 TRUE
C                  1 SET PRT 6,7,11,17,18 TRUE
C                  2 PRINT EVERYTHING, ALL PRT FALSE
C        IPUNCH
C        L1   SET BY OUTPUT.F  L1=PRT(14).AND..NOT.PUN(14)
C        L2   SET BY OUTPUT.F  L2=PRT(15).AND..NOT.PUN(15)
C        L3   CALC AND PRINT OVERLAP POPULATIONS IF TRUE
C        L4   CALC AND PRINT ENERGY MATRIX ANALYSIS IF TRUE
C        L5   USE T/F WEIGHTED HIJ FORMULA
C        CON   HUCKEL CONSTANT DEFAULT IS 1.75, MIN 1E-5
C        PEEP  HYDROGEN ORBITAL EXPONENT DEFAULT is 1.3, MIN 1E-5
C        COULH HYDROGEN H(I,I)  VALUE DEFAULT -13.6, MAX -1E5
C
C        (PRT(I),I=1,20)
C        THESE ARE COMPLICATED BY CONNECTION TO THE PUN FLAGS AND
C        SOME VALUES ARE CHANGED IN THE OUTPUT SUBROUTINE
C
C                        PRT(1) PRINT OUT ATOMIC COORDS AND PARAMS
C                        PRT(6) DO NOT PRINT RESULTS ORHUCKEL MATRIX
C                        PRT(7) DO NOT PRINT HUCKEL MATRIX
C                        PRT(8) DO NOT PRINT ENERGY LEVELS
C                        PRT(9) DO NO PRINT SUM OF 1-ELECTRON ENERGIES
C                        PRT(10) DO NOT PRINT WAVE FUNCTIONS
C                        PRT(11) DO NOT PRINT DENSITY MATRIX
C                        PRT(12) DON NOT PRINT OVERLAP MATRIX
C                        PRT(13) DO NOT PRINTREDUCED OVERLAP
C                        PRT(14) DO NOT PRINT CHARGE MATRIX FOR MO
C                        PRT(15) DO NOT PRINT REDUCED CHARGE MATRIX
C                        PRT(16) DO NOT PRING ORBITAL OCCUPATION
C                        PRT(17) DO NOT PRINT E PARTITIONING
C                        PRT(18) DO NOT PRINT REDUCED E MATRIX BY ATOM
C                        PRT(19) DO NOT PRINT REDUCED E PARTITION
C                        PRT(20) DO NOT PRINT PARTITION OR REDUCED E

C        (PUN(J),J=1,20)
C                        PUNCH OUTPUT CORRESPONDING TO PRT
C
C
C ATOM COORDINATES  X,Y,Z, NA+NH lines  FORMAT(3F15.6)
C	IT IS ASSUMED THAT THE HYDROGENS ARE FIRST, THEN THE HEAVY ATOMS
C HEAVY ATOMS FORMAT(40A2)
C           LIST OF HEAVY ATOM SYMBOLS 
C
C IF A HEAVY ATOM SYMBOL IS * THEN READ DEFINITION, IF ** RE-USE PREVIOUS DEFINITION
C USER ATOMS  SYMBOL,VELEC,NS,EXPS,COULS,NP.EXPP,COULP,ND,EXPD,COULD,C1,EXPD2,C2
C                    FORMAT(A2,I3,3(I3,2F6.3),F6.4,F6.3,F6.4)
C        USED TO DEFINE PARAMETERS NOT STORED IN THE PROGRAM
C
C IF CHARGE ITERATION READ THIS LINE
C READ DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,PRTCYC,NCON,ITEMP  FORMAT(6F10.5,3I5,4X,I1)
C ITERATION PARAMETERS
C
C READ SYMBOLS FOR CHARGE ITERATION
C   CHANGE FORMAT(20A2)
C
C READ VSIE AND MADELUNG PARAMETERS
C AS, BS1, CS1 MADS  FORMAT(4F10.8)
C AP1,BP1, CP1, MADP FORMAT(4F10.8)
C AP2, BP2,CP2, AP3,BP3,CP3   FORMAT(3F10.8)
C AD1, BD1 CD1,MADD FORMAT(4F10.8)
C AD2, BD2, CD2, AD3,BD3, CD3   FORMAT(3F10.8)
C
C
 
      IMPLICIT REAL*8(A-H,O-Z)
C
C     SUBROUTINE FOR READING IN AND PRINTING OUT INPUT DATA.
C
      INCLUDE 'PARAMETERS'
      INCLUDE 'ATOMCOMMON'

      COMMON/CNTRL/CON,PEEP,COULH,NH,NA,NATM,KA,NELEC,METH,
     .IPRINT,IPUNCH,L1,L2,L3,L4,L5,ONEMAT,ITERAT
      LOGICAL*1 L1,L2,L3,L4,L5,ONEMAT,ITERAT
      COMMON/OUT/PRT(20),PUN(20),IOVPOP(24),IENRGY(24)
      LOGICAL*1 PRT,PUN
      INTEGER*2 IOVPOP,IENRGY
      CHARACTER*32 UNIT13
      CHARACTER*29 AB
      
      COMMON/ITPARM/DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,
     .PRTCYC,ICYCLE,NCON,PARTIT,PRINTX,ITABLE(230)
      REAL*8 LAMPRI
      INTEGER*4 PRTCYC
      LOGICAL*1 PARTIT,PRINTX,ITABLE
C
C     SINCE THE ITERNAL ATOMIC PARAMETERS ( EXPS, EXPP, ETC. ) ARE
C     NOT USED WHEN DOING CHARGE ITERATION ( METH >1 ) THE SPACE
C     ALLOCATED TO THEM CAN BE USED FOR THE VSIE CHARGE ITERATION
C     PARAMETERS.
C
      DIMENSION AS1(MXUSER),BS1(MXUSER),CS1(MXUSER),AP1(MXUSER),
     .BP1(MXUSER),CP1(MXUSER),AD1(MXUSER),BD1(MXUSER),CD1(MXUSER)         
      EQUIVALENCE (AS1(1),EXPS(MXUSR2)),(BS1(1),EXPP(MXUSR2)),
     .(CS1(1),EXPD(MXUSR2)),(AP1(1),EXPD2(MXUSR2)),
     .(BP1(1),C1(MXUSR2)),(CP1(1),C2(MXUSR2)),
     .(AD1(1),COULS(MXUSR2)),(BD1(1),COULP(MXUSR2)), 
     .(CD1(1),COULD(MXUSR2))

      COMMON/ABC/AS2(5),BS2(5),CS2(5),AP2(5),BP2(5),CP2(5),AD2(5),
     .BD2(5),CD2(5),AS3(5),BS3(5),CS3(5),AP3(5),BP3(5),CP3(5),
     .AD3(5),BD3(5),CD3(5)

      CHARACTER*2 CHANGE(MXUSER)
      EQUIVALENCE (CHANGE(1),AS2(1))
      REAL*4 MADS(MXUSER),MADP(MXUSER),MADD(MXUSER)
      EQUIVALENCE (MADS(1),NS(MXUSR2)),(MADP(1),NP(MXUSR2)),
     .(MADD(1),ND(MXUSR2))                              
      
      COMMON/STARS/STAR,STAR2
      CHARACTER*2 STAR,STAR2

      COMMON/START/NUSER
      DIMENSION EXTRA(9)
      EQUIVALENCE (X(1),EXTRA(1))
      CHARACTER*2 CONTIN
      EQUIVALENCE (AB,CONTIN)
      CHARACTER*2 HYDROG
      DATA HYDROG/' H'/
 
C
C     READ AND WRITE TITLE.
C     IF CONTIN IS EQUAL TO STAR THEN ANOTHER TITLE CARD WILL
C     BE READ AND PRINTED. HOWEVER ONLY THE FIRST IS STORED
C     FOR PRINTING LATER ON. GET YOUR GOODIES ON THE FIRST.
C
      READ(5,1,END=115) AB
1     FORMAT(A29)
C

2730  FORMAT(8A8,A6,A2)
C
      WRITE(6,2) AB
2     FORMAT(T10,A32)
11    IF(CONTIN.NE.STAR) GO TO 9
      READ(5,1) EXTRA,CONTIN
      WRITE(6,12) EXTRA,CONTIN
12    FORMAT(T10,8A8,A6,A2)
      GO TO 11
9     CONTINUE
C
C
C     DEFINE FILENAME FOR PUNCH OUTPUT UNIT 13
      UNIT13=TRIM(AB) // '.MO'
      OPEN(13,FILE=UNIT13)
      WRITE(6,*) 'UNIT13: ', UNIT13
C     WRITE TITLE TO DISK FILE 13.   JJN  9-3-90
      WRITE(13,*) AB

C
C     READ PARAMETER CARD.
C
      READ(5,3,ERR=333) NH,NA,KA,METH,IPRINT,IPUNCH,L1,L2,L3,L4,L5,CON,
     .PEEP,COULH,(PRT(I),I=1,20),(PUN(J),J=1,20)
3     FORMAT(6I3,5L1,F5.2,2F6.3,40L1)
C
C     INSERT DEFAULT PARAMETERS.
C
333   IF(CON.LT.1.E-05) CON=1.75D0
      IF(PEEP.LT.1.E-05) PEEP=1.3D0
      IF(COULH.GT.-1.E-05) COULH=-13.6D0
      ITERAT=METH.NE.0
      NATOM=NH+NA
      NATM=NATOM
C
C     SET IPRINT OPTION.
C
      IF(IPRINT.GT.1) GO TO 250
      PRT(6)=.TRUE.
      PRT(7)=.TRUE.
      PRT(11)=.TRUE.
      PRT(17)=.TRUE.
      PRT(19)=.TRUE.
      IF(IPRINT.GT.0) GO TO 250
      PRT(12)=.TRUE.
      PRT(14)=.TRUE.
      PRT(20)=.TRUE.
      IF(IPRINT.GT.-1) GO TO 250
      PRT(13)=.TRUE.
      PRT(15)=.TRUE.
      PRT(16)=.TRUE.
      PRT(18)=.TRUE.
      IF(IPRINT.GT.-2) GO TO 250
      PRT(10)=.TRUE.
C
C     READ COORDINATES AND HEAVY ATOM CARD.
C
250   READ(5,5) (X(I),Y(I),Z(I),I=1,NATOM)
5     FORMAT(3F15.6)
      READ(5,8) (AC(I),I=1,NA)
8     FORMAT(40A2)
 
C
C     ORIGINAL CODE MADE CHAR*2 AC and INT*2 KEY SHARE THE SAME ARRAY
C     THAT HASN'T WORKED SINCE FORTRAN IV SO WE DO IT MANUALLY
C     HERE BY LOOKING UP INDEX OF SYMBOL.
C     IT IS NOT CLEAR HOW THIS EVER COULD HAVE WORKED EVEN USING THE
C     ASCII CODES FOR THE LETTERS   MC
 
      DO 888 II=1,NA
       DO 888 JJ=1,IB
         IF(AC(II).EQ.SYMBOL(JJ)) THEN
           KEY(II) = JJ
           EXIT
         END IF
888   CONTINUE
C
C
C     READ AND DECODE ATOM DEFINITION CARDS.
C
      JOHN=0
      NDIM=NH
      NTYPE=NH
      NELEC=NH-KA
C
C    NUSER DEFINED AS 231 IN BLOCKDATA
C
      K=NUSER
      NUSER2=IB
      IF(METH.GE.2) NUSER2=MXUSER
      DO 100 I=1,NA
      IF(NUSER.GT.NUSER2) GO TO 103
 
      DO 102 J=NUSER,IB
      JSAVE=J
      IF(AC(I) .EQ. SYMBOL(J)) GO TO 101
102   CONTINUE
C
C     PROVISION FOR USER SPECIFIED DATA.
C
      IF(AC(I).EQ.STAR) GO TO 103
      IF(AC(I).EQ.STAR2) GO TO 105
      WRITE(6,6) I,AC(I)
6     FORMAT(//,T10,'HEAVY ATOM',I3,' NOT RECOGNIZED. SYMBOL0',A2)
      IF(METH.GE.2) WRITE(6,13)
13    FORMAT(/,T10,'REMEMBER IF USING METH > 1 ALL ATOMIC',
     .' PARAMETERS MUST BE DEFINED BY THE USER.')
115   REWIND 7
      STOP

103   NUSER=NUSER-1
105   READ(5,7) SYMBOL(NUSER),VELEC(NUSER),NS(NUSER),EXPS(NUSER),
     .COULS(NUSER),NP(NUSER),EXPP(NUSER),COULP(NUSER),ND(NUSER),
     .EXPD(NUSER),COULD(NUSER),C1(NUSER),EXPD2(NUSER),C2(NUSER)
7     FORMAT(A2,I3,3(I3,2F6.3),F6.4,F6.3,F6.4)
C
C    THIS SECTION WRITE OUT THE USER-DEFINED PARAMETERS TO DISK
C    FILE 18. THESE WILL BE USED FOR CONSTRUCTION OF THE PSI1
C    INPUT FILE AS THESE PARAMETERS ARE REQUIRED FOR ELEMENTS
C    GREATER THAN ATOMIC NUMBER 18.    JJN  9-8-90
C
      JOHN=JOHN+1
      WRITE(18,767) SYMBOL(NUSER),VELEC(NUSER),NS(NUSER),
     .EXPS(NUSER),NP(NUSER),EXPP(NUSER),ND(NUSER),EXPD(NUSER),
     .C1(NUSER),EXPD2(NUSER),C2(NUSER)
767   FORMAT(A2,I3,I3,F6.3,I3,F6.3,I3,F6.3,F6.4,F6.3,F6.4)
C
C
      JSAVE=NUSER
C
C     NORMALIZE USER SPECIFIED CONTRACTED D ORBITAL.
C
      IF(C2(NUSER).EQ.0.) GO TO 101
      S=(4.D0*EXPD(NUSER)*EXPD2(NUSER)/(EXPD(NUSER)+EXPD2(NUSER)
     .)**2)**(ND(NUSER)+.5D0)
      S=1.D0/DSQRT(C1(NUSER)**2+C2(NUSER)**2+(S+S)*C1(NUSER)
     .*C2(NUSER))
      C1(NUSER)=S*C1(NUSER)
      C2(NUSER)=S*C2(NUSER)
101   NELEC=NELEC+VELEC(JSAVE)
C
C     AC, LATER REFERENCED AS KEY, IS A POINTER TO THE PARAMETER TABLES.
C
      KEY(I)=JSAVE
      NDIM=NDIM+4
      IF(NP(JSAVE).EQ.0) NDIM=NDIM-3
      IF(ND(JSAVE).NE.0) NDIM=NDIM+5
      NTYPE=NTYPE+2
      IF(NP(JSAVE).EQ.0) NTYPE=NTYPE-1
      IF(ND(JSAVE).NE.0) NTYPE=NTYPE+1
100   CONTINUE
C
C     READ IN CHARGE ITERATION PARAMETERS. SET DEFAULT VALUES.
C
      IF(.NOT.ITERAT) GO TO 60
      IF(K.NE.MXUSR2) GO TO 60
      PARTIT = .FALSE.
      READ(5,61) DAMP1,DAMP2,DAMP3,LAMPRI,DELTAC,SENSE,MAXCYC,
     .PRTCYC,NCON,ITEMP
61    FORMAT(6F10.5,3I5,4X,I1)
C       
C          READINT INT IS LESS ERROR PRONE
C
      IF (ITEMP .NE. 0 ) PARTIT = .TRUE.
 
      IF(.NOT.PARTIT.OR.METH.EQ.2) GO TO 65
      WRITE(6,66)
66    FORMAT(///,T10,'PARTIAL ITERATION ( PARTIT = TRUE ) MAY',
     .' ONLY BE USED IF METH = 2.')
      STOP
65    IF(DELTAC.EQ.0.0D0) DELTAC=0.0001D0
      IF(SENSE.EQ.0.0D0) SENSE=2.0D0
      IF(MAXCYC.EQ.0) MAXCYC=100
      IF(PRTCYC.EQ.0) PRTCYC=MAXCYC
      IF(NCON.EQ.0) NCON=3
      IF(DAMP1.EQ.0.0D0) DAMP1=0.1D0
      IF(METH.GE.3) GO TO 62
      IF(DAMP2.EQ.0.0D0) DAMP2=0.25D0
      IF(LAMPRI.EQ.0.0D0) LAMPRI=0.25D0
      GO TO 63
62    IF(DAMP2.EQ.0.0D0) DAMP2=0.75D0
      IF(LAMPRI.EQ.0.0D0) LAMPRI=0.75D0
63    IF(METH.LT.2) GO TO 60
      DO 32 I=1,20
32    ITABLE(I)=.FALSE.
      NUSER2=MXUSR2-NUSER
C
C     READ IN SYMBOLS OF ATOMS ON WHICH CHARGE ITERATION
C     IS TO BE PERFORMED.
C
      IF(.NOT.PARTIT) GO TO 30
      READ(5,31) CHANGE
31    FORMAT(20A2)
      DO 33 I=1,NUSER2
      J=MXUSR2-I
      DO 33 K=1,NUSER2
      IF(SYMBOL(J).EQ.CHANGE(K)) ITABLE(J)=.TRUE.
33    CONTINUE
      GO TO 34
30    DO 35 I=1,NUSER2
      J=MXUSR2-I
35    ITABLE(J)=.TRUE.
C
C     READ IN VSIE AND MADELUNG PARAMETERS.
C
34    DO 36 I=1,NUSER2
      J=MXUSR2-I
      IF(.NOT.ITABLE(J)) GO TO 36
      READ(5,37) AS1(I),BS1(I),CS1(I),MADS(I)
37    FORMAT(4F10.8)
      IF(NP(J).EQ.0) GO TO 36
      IF(ND(J).NE.0) GO TO 38
      READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)
      GO TO 36
38    IF(NCON.EQ.3) GO TO 39
      READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I),AD1(I),BD1(I),
     .CD1(I),MADD(I)
      GO TO 36
39    READ(5,40) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I)
40    FORMAT(3F10.8)
      READ(5,37) AP1(I),BP1(I),CP1(I),MADP(I)
      READ(5,40) AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I)
      READ(5,37) AD1(I),BD1(I),CD1(I),MADD(I)
      READ(5,40) AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),CD3(I)
36    CONTINUE
C
C     READ IN IOVPOP(I) AND IENRGY(I). INDIVIDULAL OVERLAP POPULATION
C     ANALYSES ARE PERFORMED FROM ORBITAL IOVPOP(N) TO ORBITAL
C     IOVPOP(N+1). INDIVIDUAL ENERGY MATRIX ANALYSES ARE PERFORMED
C     FROM ORBITAL IENRGY(N) TO ORBITAL IENRGY(N+1).
C
60    IF(L3) READ(5,67) IOVPOP
67    FORMAT(24I3)
      IF(L4) READ(5,67) IENRGY
C
C     PRINT OUT TYPE OF CALCULATION.
C
      IF(METH.EQ.0) WRITE(6,90)
      IF(METH.EQ.1) WRITE(6,91)
      IF(METH.GE.2) WRITE(6,92)
      IF(METH.GT.2) WRITE(6,93)
      IF(L5) WRITE(6,94)
90    FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION.')
91    FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',
     .' ITERATION.',/,T10,'LINEAR CHARGE DEPENDENCE OF SENSE*CHARG'
     .,'E FOR H(I,I)''S.')
92    FORMAT(///,T10,'EXTENDED HUCKEL CALCULATION WITH CHARGE',
     .' ITERATION.')
93    FORMAT(T10,'MADELUNG CORRECTION INCLUDED.')
94    FORMAT(T10,'WEIGHTED HIJ FORMULA USED.')
C
C     PRINT OUT ATOMIC COORDINATES AND PARAMETERS.
C
      IF(PRT(1)) GO TO 80
      WRITE(6,74)
74    FORMAT(///,T5,'ATOM',T17,'X',T29,'Y',T41,'Z',T56,'S',T76,'P',
     .T96,'D',T113,'CONTRACTED D'/T47,'N',T50,'EXP',T59,'COUL',
     .T67,'N',T70,'EXP',T79,'COUL',T87,'N',T90,'EXPD1',T99,'COUL',
     .T109,'C1',T118,'C2',T125,'EXPD2')
      IF(NH.EQ.0) GO TO 72
      J=1
      DO 76 I=1,NH
76    WRITE(6,53) HYDROG,I,X(I),Y(I),Z(I),J,PEEP,COULH
53    FORMAT(T4,A2,I3,3F12.5,3(I3,F8.4,F9.4),2F9.5,F8.4)
C
C  ** THIS CHUNK WRITES ANY HYDROGEN ATOM COORDINATES TO DISK FILE 13
C     THAT ARE DEFINED SPECIFICALLY AS HYDROGENS (I.E., NOT DEFINED AS
C     HEAVY ATOMS). THE HYDROGEN AS A HEAVY ATOM CASE COMES LATER.
C     IT ALSO WRITES THE NUMBER OF ATOMS AND THE NUMBER OF VALENCE
C     ORBITALS (ALSO MOLECULAR ORBITALS) FOR USE IN SETTING UP THE
C     PLOTTING FILES.       JJN   8-28-90
C
      IF (NH.EQ.0) GO TO 72
C     NDIM HERE IS THE NUMBER OF ATOMIC ORBITALS
C     JOHN IS THE NUMBER OF USER PARAMETER CARDS READ
C     KA IS TOTAL CHARGE
      WRITE(13,9998) NATOM,NDIM,NELEC,KA,JOHN
9998  FORMAT(I3,2I4,4I3)
      DO 9991 I=1,NH
9991  WRITE(13,9992) X(I),Y(I),Z(I)
9992  FORMAT(' 1',3F12.6)
C
C  **
C
72    CONTINUE
      DO 151 I=1,NA
      KEYI=KEY(I)
      INH=I+NH
      IF(NP(KEYI).NE.0) GO TO 152
      WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .EXPS(KEYI),COULS(KEYI)
      GO TO 151
152   IF(ND(KEYI).NE.0) GO TO 153
      WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI)
      GO TO 151
153   IF(C2(KEYI).NE.0.0D0) GO TO 154
      WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
     .ND(KEYI),EXPD(KEYI),COULD(KEYI)
      GO TO 151
154   WRITE(6,53) SYMBOL(KEYI),INH,X(INH),Y(INH),Z(INH),NS(KEYI),
     .EXPS(KEYI),COULS(KEYI),NP(KEYI),EXPP(KEYI),COULP(KEYI),
     .ND(KEYI),EXPD(KEYI),COULD(KEYI),C1(KEYI),C2(KEYI),EXPD2(KEYI)
151   CONTINUE
C
C     THIS CHUNK WRITES THE COORDINATES AND ATOMIC NUMBER 
C     (FOR ELEMENTS LESS THAN OR EQUAL TO 118) TO DISK FILE
C     13 FOR NON-HYDROGEN ATOMS.     JJN  9-8-90
C
C     THE FIRST ITEMS WRITTEN OUT ARE THE NUMBER OF ATOMS AND THE
C     NUMBER OF VALENCE ORBITALS (I.E., MOLECULAR ORBITALS) FOR USE
C     IN CONSTRUCTING THE PLOTTING ROUTINE INPUT FILES.
C
      IF (NH.NE.0) GO TO 10051
      WRITE(13,10050) NATOM,NDIM,NELEC,KA,JOHN
10050 FORMAT(I3,2I4,4I3)
10051 DO 9996 I=1,NA
      KEYI=KEY(I)
      INH=I+NH
      DO 9995 J=1,118
      IF (SYMBOL(KEYI).EQ.ATS(J)) THEN
      SYMBL(KEYI)=J
      END IF
9995  CONTINUE
      WRITE(13,9994) SYMBL(KEYI),X(INH),Y(INH),Z(INH)
9994  FORMAT(I2,3F12.6)
9996  CONTINUE
C
C  
C
      WRITE(6,160) KA,IPRINT,IPUNCH,CON
160   FORMAT(///,T10,'CHARGE =',I3,8X,'IPRINT =',I3,8X,'IPUNCH =',
     .I3,8X,'HUCKEL CONSTANT =',F7.3)
C
C     PRINT OUT ITERATION PARAMETERS.
C
      IF(.NOT.ITERAT) GO TO 80
      WRITE(6,81) DAMP1,DAMP2,DAMP3,LAMPRI,MAXCYC,PRTCYC,
     .SENSE,DELTAC
81    FORMAT(/,T10,'DAMP1 =',F6.3,6X,'DAMP2 =',F6.3,6X,'DAMP3 =',
     .F6.3,6X,'LAMPRI =',F6.3,//,T10,'MAXCYC =',I3,8X,'PRTCYC =',
     .I3,8X,'SENSE =',F6.3,6X,'DELTAC =',F10.7)
C
C     PRINT OUT VSIE PARAMETERS.
C
      IF(METH.LT.2) GO TO 80
      WRITE(6,82)
82    FORMAT(///,' VSIE PARAMETERS',//,T10,'ATOM',T26,'A',T39,'B',
     .T52,'C')
      NUSER2=MXUSR2-NUSER
      DO 83 I=1,NUSER2
      J=MXUSR2-I
      IF(.NOT.ITABLE(J)) GO TO 83
      WRITE(6,84) SYMBOL(J),AS1(I),BS1(I),CS1(I)
84    FORMAT(/,T11,A2,4X,3F13.5)
      IF(NP(J).EQ.0) GO TO 83
      IF(ND(J).NE.0) GO TO 85
      WRITE(6,86) AP1(I),BP1(I),CP1(I)
86    FORMAT(T17,3F13.5)
      GO TO 83
85    IF(NCON.EQ.3) GO TO 87
      WRITE(6,86) AP1(I),BP1(I),CP1(I),AD1(I),BD1(I),CD1(I)
      GO TO 83
87    WRITE(6,86) AS2(I),BS2(I),CS2(I),AS3(I),BS3(I),CS3(I),AP1(I),
     .BP1(I),CP1(I),AP2(I),BP2(I),CP2(I),AP3(I),BP3(I),CP3(I),
     .AD1(I),BD1(I),CD1(I),AD2(I),BD2(I),CD2(I),AD3(I),BD3(I),
     .CD3(I)
83    CONTINUE
80    WRITE(6,99)
99    FORMAT(///)
      RETURN
      END

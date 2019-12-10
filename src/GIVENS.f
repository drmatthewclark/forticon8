C
C     FORTICON8  JANUARY 2018
C     MATTHEW CLARK
C
      SUBROUTINE GIVENS(NX,NROOTX,NJX,A,B,ROOT,VECT)
 
c     ---------------------------------------------------
 
      IMPLICIT REAL*8 (A-H,O-Z)
 
C     corrected code for overflows fetched from here:
C     http://scs.illinois.edu/~makri/New-Web-Site/550-web-site/givens.f
 
      DIMENSION B(NX,6),A(NX*(NX+1)/2),VECT(NX,NROOTX),ROOT(NX)
 
C      CALCULATES EIGENVALUES AND EIGENVECTORS OF REAL SYMMETRIC MATRIX 
C      STORED IN PACKED UPPER TRIANGULAR FORM.                          
C     62.3  GIVENS  -EIGENVALUES AND EIGENVECTORS BY THE GIVENS METHOD. 
C      BY FRANKLIN PROSSER, INDIANA UNIVERSITY.                         
C      SEPTEMBER, 1967                                                  
C                                                                       
C      THANKS ARE DUE TO F. E. HARRIS (STANFORD UNIVERSITY) AND H. H.   
C      MICHELS (UNITED AIRCRAFT RESEARCH LABORATORIES) FOR EXCELLENT    
C      WORK ON NUMERICAL DIFFICULTIES WITH EARLIER VERSIONS OF THIS     
C      PROGRAM.                                                         
C                                                                       
C      THE PARAMETERS FOR THE ROUTINE ARE...                            
C          NX     ORDER OF MATRIX                                       
C          NROOTX NUMBER OF ROOTS WANTED.  THE NROOTX SMALLEST (MOST    
C                  NEGATIVE) ROOTS WILL BE CALCULATED.  IF NO VECTORS   
C                  ARE WANTED, MAKE THIS NUMBER NEGATIVE.               
C          NJX    ROW DIMENSION OF VECT ARRAY.  SEE :VECT: BELOW.       
C                  NJX MUST BE NOT LESS THAN NX.                        
C          A      MATRIX STORED BY COLUMNS IN PACKED UPPER TRIANGULAR   
C                 FORM, I.E. OCCUPYING NX*(NX+1)/2 CONSECUTIVE          
C                 LOCATIONS.                                            
C          B      SCRATCH ARRAY USED BY GIVENS.  MUST BE AT LEAST       
C                  NX*5 CELLS.                                          
C          ROOT   ARRAY TO HOLD THE EIGENVALUES.  MUST BE AT LEAST      
C                 NROOTX CELLS LONG.  THE NROOTX SMALLEST ROOTS ARE     
C                  ORDERED LARGEST FIRST IN THIS ARRAY.                 
C          VECT   EIGENVECTOR ARRAY.  EACH COLUMN WILL HOLD AN          
C                  EIGENVECTOR FOR THE CORRESPONDING ROOT.  MUST BE     
C                  DIMENSIONED WITH :NJX: ROWS AND AT LEAST :NROOTX:    
C                  COLUMNS, UNLESS NO VECTORS                           
C                  ARE REQUESTED (NEGATIVE NROOTX).  IN THIS LATTER     
C                  CASE, THE ARGUMENT VECT IS JUST A DUMMY, AND THE     
C                  STORAGE IS NOT USED.                                 
C                                                                       
C      THE ARRAYS A AND B ARE DESTROYED BY THE COMPUTATION.  THE RESULTS
C      APPEAR IN ROOT AND VECT.                                         
C      FOR PROPER FUNCTIONING OF THIS ROUTINE, THE RESULT OF A FLOATING 
C      POINT UNDERFLOW SHOULD BE A ZERO.                                
C                                                                       
C      THE ORIGINAL REFERENCE TO THE GIVENS TECHNIQUE IS IN OAK RIDGE   
C      REPORT NUMBER ORNL 1574 (PHYSICS), BY WALLACE GIVENS.            
C      THE METHOD AS PRESENTED IN THIS PROGRAM CONSISTS OF FOUR STEPS,  
C      ALL MODIFICATIONS OF THE ORIGINAL METHOD...                      
C      FIRST, THE INPUT MATRIX IS REDUCED TO TRIDIAGONAL FORM BY THE    
C      HOUSEHOLDER TECHNIQUE (J. H. WILKINSON, COMP. J. 3, 23 (1960)).  
C      THE ROOTS ARE THEN LOCATED BY THE STURM SEQUENCE METHOD (J. M.   
C      ORTEGA (SEE REFERENCE BELOW).  THE VECTORS OF THE TRIDIAGONAL    
C      FORM ARE THEN EVALUATED (J. H. WILKINSON, COMP. J. 1, 90 (1958)),
C      AND LAST THE TRIDIAGONAL VECTORS ARE ROTATED TO VECTORS OF THE   
C      ORIGINAL ARRAY (FIRST REFERENCE).                                
C      VECTORS FOR DEGENERATE (OR NEAR-DEGENERATE) ROOTS ARE FORCED     
C      TO BE ORTHOGONAL, USING A METHOD SUGGESTED BY B. GARBOW, ARGONNE 
C      NATIONAL LABS (PRIVATE COMMUNICATION, 1964).  THE GRAM-SCHMIDT   
C      PROCESS IS USED FOR THE ORTHOGONALIZATION.                       
C                                                                       
C      AN EXCELLENT PRESENTATION OF THE GIVENS TECHNIQUE IS FOUND IN    
C      J. M. ORTEGA:S ARTICLE IN :MATHEMATICS FOR DIGITAL COMPUTERS,:   
C      VOLUME 2, ED. BY RALSTON AND WILF, WILEY (1967), PAGE 94.        
C                                                                       
C                                                                       
C      ALL REAL LIBRARY FUNCTIONS AND INTERNAL FUNCTIONS ARE DEFINED    
C      FOLLOWING STATEMENT FUNCTIONS.  THIS IS TO FACILITATE CONVERSION 
C      OF THIS ROUTINE TO DOUBLE PRECISION ON IBM 360 MACHINES.         
C      TO ACCOMPLISH THIS, CHANGE THE FUNCTION DEFINITIONS TO THEIR     
C      DOUBLE PRECISION COUNTERPARTS, AND ADD AN :IMPLICIT: STATEMENT   
C      OF THE FORM...      IMPLICIT REAL*8(A-H),REAL*8(O-Z)             
C     ZSQRT(X)=SQRT(X)                                                  
C     ZABS(X)=ABS(X)                                                    
C     ZMAX1(X,Y)=AMAX1(X,Y)                                             
C     ZMIN1(X,Y)=AMIN1(X,Y)                                             
C     ZSIGN(X,Y)=SIGN(X,Y)                                              
C                                                                       
C      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
C      USERS PLEASE NOTE...                                             
C      THE FOLLOWING TWO PARAMETERS, ETA AND THETA, SHOULD BE ADJUSTED  
C      BY THE USER FOR HIS PARTICULAR MACHINE.                          
C      ETA IS AN INDICATION OF THE PRECISION OF THE FLOATING POINT      
C      REPRESENTATION ON THE COMPUTER BEING USED (ROUGHLY 10**(-M),     
C      WHERE M IS THE NUMBER OF DECIMALS OF PRECISION ).                
C      THETA IS AN INDICATION OF THE RANGE OF NUMBERS THAT CAN BE       
C      EXPRESSED IN THE FLOATING POINT REPRESENTATION (ROUGHLY THE      
C      LARGEST NUMBER).                                                 
C      SOME RECOMMENDED VALUES FOLLOW.                                  
C      FOR CONTROL DATA 3600 (36-BIT BINARY FRACTION, 11-BIT BINARY     
C      EXPONENT), ETA=1.E-11, THETA=1.E307.                             
C      FOR CONTROL DATA 6600 (48-BIT BINARY FRACTION, 11-BIT BINARY     
C      EXPONENT), ETA=1.E-14, THETA=1.E307.                             
C      FOR IBM 7094, UNIVAC 1108, ETC. (27-BIT BINARY FRACTION, 8-BIT   
C      BINARY EXPONENT), ETA=1.E-8, THETA=1.E37.                        
C      FOR IBM 360/50 AND 360/65 DOUBLE PRECISION (56-BIT HEXADECIMAL   
C      FRACTION, 7-BIT HEXADECIMAL EXPONENT), ETA=1.E-16, THETA=1.E75.  
C                                                                       
C      ETA = 1.D-15                                                      
C      THETA = 1.D+65                                                    
C     IEEE values (INTEL) THETA REDUCED TO PREVENT OVERFLOW
 
      ETA = 2.22045D-16
      THETA =  1.7976931348623157D304
      TEMP=0.D+00
C      * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      DO 20 I=1,NJX
      DO 10 J=1,NROOTX
      VECT(I,J)=0.0D+00
   10 CONTINUE
   20 CONTINUE
      DO 40 I=1,NX
      DO 30 J=1,5
      B(I,J)=0.0D+00
   30 CONTINUE
   40 CONTINUE
C                                                                       
      DEL1=ETA/100.D+00
      DELTA=ETA**2*100.D+00
      SMALL=ETA**2/100.D+00
      DELBIG=THETA*DELTA/1000.D+00
      THETA1=1.D-03/THETA
C      TOLER  IS A FACTOR USED TO DETERMINE IF TWO ROOTS ARE CLOSE      
C      ENOUGH TO BE CONSIDERED DEGENERATE FOR PURPOSES OF ORTHOGONALI-  
C      ZING THEIR VECTORS.  FOR THE MATRIX NORMED TO UNITY, IF THE      
C      DIFFERENCE BETWEEN TWO ROOTS IS LESS THAN TOLER, THEN            
C      ORTHOGONALIZATION WILL OCCUR.                                    
      TOLER=ETA*10000.0D+00
C                                                                       
C      INITIAL VALUE FOR PSEUDORANDOM NUMBER GENERATOR... (2**23)-3     
      RAND1=8388605.0D+00
C                                                                       
      N=NX
      NROOT=IABS(NROOTX)
      IF (NROOT.EQ.0) GO TO 700
      IF (N-1) 700,50,60
   50 ROOT(1)=A(1)
      IF (NROOTX.GT.0) VECT(1,1)=1.0D+00
      GO TO 700
   60 CONTINUE
C     NSIZE    NUMBER OF ELEMENTS IN THE PACKED ARRAY                   
      NSIZE=(N*(N+1))/2
      NM1=N-1
      NM2=N-2
C                                                                       
C     SCALE MATRIX TO EUCLIDEAN NORM OF 1.  SCALE FACTOR IS ANORM.      
      FACTOR=0.D+00
      DO 70 I=1,NSIZE
      ZZ=A(I)
      ZZ = ABS(ZZ)
   70 FACTOR=MAX(FACTOR,ZZ)
      IF (FACTOR.NE.0.0D+00) GO TO 100
C     NULL MATRIX.  FIX UP ROOTS AND VECTORS, THEN EXIT.                
      DO 90 I=1,NROOT
      IF (NROOTX.LT.0) GO TO 90
      DO 80 J=1,N
   80 VECT(J,I)=0.D+00
      VECT(I,I)=1.0D+00
   90 ROOT(I)=0.D+00
      GO TO 700
C                                                                       
  100 ANORM=0.D+00
      J=1
      K=1
      DO 120 I=1,NSIZE
      IF (I.NE.J) GO TO 110
      ANORM=ANORM+(A(I)/FACTOR)**2/2.D+00
      K=K+1
      J=J+K
      GO TO 120
  110 ANORM=ANORM+(A(I)/FACTOR)**2
  120 CONTINUE
      ANORM=SQRT(ANORM*2.0D+00)*FACTOR
      DO 130 I=1,NSIZE
  130 A(I)=A(I)/ANORM
      ALIMIT=1.0D+00
C                                                                       
C      TRIDIA SECTION.                                                  
C      TRIDIAGONALIZATION OF SYMMETRIC MATRIX                           
      ID=0
      IA=1
      IF (NM2.EQ.0) GO TO 240
      DO 230 J=1,NM2
C      J       COUNTS ROW  OF A-MATRIX TO BE DIAGONALIZED               
C      IA      START OF NON-CODIAGONAL ELEMENTS IN THE ROW              
C      ID      INDEX OF CODIAGONAL ELEMENT ON ROW BEING CODIAGONALIZED. 
      IA=IA+J+2
      ID=ID+J+1
      JP2=J+2
C      SUM SQUARES OF NON-CODIAGONAL ELEMENTS IN ROW J                  
      II=IA
      SUM=0.0D+00
      DO 140 I=JP2,N
      SUM=SUM+A(II)**2
  140 II=II+I
      TEMP=A(ID)
      IF (SUM.GT.SMALL) GO TO 150
C      NO TRANSFORMATION NECESSARY IF ALL THE NON-CODIAGONAL            
C      ELEMENTS ARE TINY.                                               
      B(J,1)=TEMP
      A(ID)=0.D+00
      GO TO 230
C      NOW COMPLETE THE SUM OF OFF-DIAGONAL SQUARES                     
  150 SUM=SQRT(SUM+TEMP**2)
C      NEW CODIAGONAL ELEMENT                                           
      B(J,1)=-SIGN(SUM,TEMP)
C      FIRST NON-ZERO ELEMENT OF THIS W-VECTOR                          
      B(J+1,2)=SQRT((1.0D+00+ABS(TEMP)/SUM)/2.0D+00)
C      FORM REST OF THE W-VECTOR ELEMENTS                               
      TEMP=SIGN(0.5/(B(J+1,2)*SUM),TEMP)
      II=IA
      DO 160 I=JP2,N
      B(I,2)=A(II)*TEMP
  160 II=II+I
C      FORM P-VECTOR AND SCALAR.  P-VECTOR = A-MATRIX*W-VECTOR.         
C      SCALAR = W-VECTOR*P-VECTOR.                                      
      AK=0.0D+00
C      IC      LOCATION OF NEXT DIAGONAL ELEMENT                        
      IC=ID+1
      J1=J+1
      DO 190 I=J1,N
      JJ=IC
      TEMP=0.D+00
      DO 180 II=J1,N
C      I       RUNS OVER THE NON-ZERO P-ELEMENTS                        
C      II      RUNS OVER ELEMENTS OF W-VECTOR                           
      TEMP=TEMP+B(II,2)*A(JJ)
C      CHANGE INCREMENTING MODE AT THE DIAGONAL ELEMENTS.               
      IF (II.LT.I) GO TO 170
      JJ=JJ+II
      GO TO 180
  170 JJ=JJ+1
  180 CONTINUE
C      BUILD UP THE K-SCALAR (AK)                                       
      AK=AK+TEMP*B(I,2)
      B(I,1)=TEMP
C      MOVE IC TO TOP OF NEXT A-MATRIX :ROW:                            
  190 IC=IC+I
C      FORM THE Q-VECTOR                                                
      DO 200 I=J1,N
  200 B(I,1)=B(I,1)-AK*B(I,2)
C      TRANSFORM THE REST OF THE A-MATRIX                               
C      JJ      START-1 OF THE REST OF THE A-MATRIX                      
      JJ=ID
C      MOVE W-VECTOR INTO THE OLD A-MATRIX LOCATIONS TO SAVE SPACE      
C      I       RUNS OVER THE SIGNIFICANT ELEMENTS OF THE W-VECTOR       
      DO 220 I=J1,N
      A(JJ)=B(I,2)
      DO 210 II=J1,I
      JJ=JJ+1
  210 A(JJ)=A(JJ)-2.0D+00*(B(I,1)*B(II,2)+B(I,2)*B(II,1))
  220 JJ=JJ+J
  230 CONTINUE
C      MOVE LAST CODIAGONAL ELEMENT OUT INTO ITS PROPER PLACE           
  240 CONTINUE
      B(NM1,1)=A(NSIZE-1)
      A(NSIZE-1)=0.D+00
C                                                                       
C     STURM SECTION.                                                    
C     STURM SEQUENCE ITERATION TO OBTAIN ROOTS OF TRIDIAGONAL FORM.     
C     MOVE DIAGONAL ELEMENTS INTO SECOND N ELEMENTS OF B-VECTOR.        
C     THIS IS A MORE CONVENIENT INDEXING POSITION.                      
C     ALSO, PUT SQUARE OF CODIAGONAL ELEMENTS IN THIRD N ELEMENTS.      
      JUMP=1
      DO 250 J=1,N
      B(J,2)=A(JUMP)
      B(J,3)=B(J,1)**2
  250 JUMP=JUMP+J+1
      DO 260 I=1,NROOT
  260 ROOT(I)=+ALIMIT
      ROOTL=-ALIMIT
C     ISOLATE THE ROOTS.  THE NROOT LOWEST ROOTS ARE FOUND, LOWEST FIRST
      DO 340 I=1,NROOT
C     FIND CURRENT :BEST: UPPER BOUND                                   
      ROOTX=+ALIMIT
      DO 270 J=I,NROOT
  270 ROOTX=MIN(ROOTX,ROOT(J))
      ROOT(I)=ROOTX
C     GET IMPROVED TRIAL ROOT                                           
  280 TRIAL=(ROOTL+ROOT(I))*0.5D+00
      IF (TRIAL.EQ.ROOTL.OR.TRIAL.EQ.ROOT(I)) GO TO 340
C     FORM STURM SEQUENCE RATIOS, USING ORTEGA:S ALGORITHM (MODIFIED).  
C     NOMTCH IS THE NUMBER OF ROOTS LESS THAN THE TRIAL VALUE.          
      NOMTCH=N
      J=1
  290 F0=B(J,2)-TRIAL
  300 CONTINUE
      IF (ABS(F0).LT.THETA1) GO TO 310
      IF (F0.GE.0.0D+00) NOMTCH=NOMTCH-1
      J=J+1
      IF (J.GT.N) GO TO 320
C     SINCE MATRIX IS NORMED TO UNITY, MAGNITUDE OF B(J,3) IS LESS THAN 
C     ONE, SO OVERFLOW IS NOT POSSIBLE AT THE DIVISION STEP, SINCE      
C     F0 IS GREATER THAN THETA1.                                        
      F0=B(J,2)-TRIAL-B(J-1,3)/F0
      GO TO 300
  310 J=J+2
      NOMTCH=NOMTCH-1
      IF (J.LE.N) GO TO 290
  320 CONTINUE
C     FIX NEW BOUNDS ON ROOTS                                           
      IF (NOMTCH.GE.I) GO TO 330
      ROOTL=TRIAL
      GO TO 280
  330 ROOT(I)=TRIAL
      NOM=MIN0(NROOT,NOMTCH)
      ROOT(NOM)=TRIAL
      GO TO 280
  340 CONTINUE
C     REVERSE THE ORDER OF THE EIGENVALUES, SINCE CUSTOM DICTATES       
C     :LARGEST FIRST:.  THIS SECTION MAY BE REMOVED IF DESIRED WITHOUT  
C     AFFECTING THE REMAINDER OF THE ROUTINE.                           
      NRT=NROOT/2
      DO 350 I=1,NRT
      SAVE=ROOT(I)
      NMIP1=NROOT-I+1
      ROOT(I)=ROOT(NMIP1)
  350 ROOT(NMIP1)=SAVE
C     TRIVEC SECTION.                                                   
C     EIGENVECTORS OF CODIAGONAL FORM                                   
C     QUIT NOW IF NO VECTORS WERE REQUESTED.                            
      IF (NROOTX.LT.0) GO TO 680
C     INITIALIZE VECTOR ARRAY.                                          
      DO 360 I=1,N
      DO 360 J=1,NROOT
  360 VECT(I,J)=1.0D+00
      DO 620 I=1,NROOT
      AROOT=ROOT(I)
C     ORTHOGONALIZE IF ROOTS ARE CLOSE.                                 
      IF (I.EQ.1) GO TO 370
C     THE ABSOLUTE VALUE IN THE NEXT TEST IS TO ASSURE THAT THE TRIVEC  
C     SECTION IS INDEPENDENT OF THE ORDER OF THE EIGENVALUES.           
      IF (ABS(ROOT(I-1)-AROOT).LT.TOLER) GO TO 380
  370 IA=-1
  380 IA=IA+1
      ELIM1=A(1)-AROOT
      ELIM2=B(1,1)
      JUMP=1
      DO 410 J=1,NM1
      JUMP=JUMP+J+1
C     GET THE CORRECT PIVOT EQUATION FOR THIS STEP.                     
      IF (ABS(ELIM1).LE.ABS(B(J,1))) GO TO 390
C     FIRST (ELIM1) EQUATION IS THE PIVOT THIS TIME.  CASE 1.           
      B(J,2)=ELIM1
      B(J,3)=ELIM2
      B(J,4)=0.D+00
      TEMP=B(J,1)/ELIM1
      ELIM1=A(JUMP)-AROOT-TEMP*ELIM2
      ELIM2=B(J+1,1)
      GO TO 400
C     SECOND EQUATION IS THE PIVOT THIS TIME.  CASE 2.                  
  390 B(J,2)=B(J,1)
      B(J,3)=A(JUMP)-AROOT
      B(J,4)=B(J+1,1)
      TEMP=1.0D+00
      IF (ABS(B(J,1)).GT.THETA1) TEMP=ELIM1/B(J,1)
      ELIM1=ELIM2-TEMP*B(J,3)
      ELIM2=-TEMP*B(J+1,1)
C     SAVE FACTOR FOR SECOND ITERATION.                                 
  400 B(J,5)=TEMP
  410 CONTINUE
      B(N,2)=ELIM1
      B(N,3)=0.D+00
      B(N,4)=0.D+00
      B(NM1,4)=0.D+00
      ITER=1
      IF (IA.NE.0) GO TO 510
C     BACK SUBSTITUTE TO GET THIS VECTOR.                               
  420 L=N+1
      DO 460 J=1,N
      L=L-1
  430 CONTINUE
      ELIM1=VECT(L,I)
      IF (L.LT.N) ELIM1=ELIM1-VECT(L+1,I)*B(L,3)
      IF (L.LT.N-1) ELIM1=ELIM1-VECT(L+2,I)*B(L,4)
C     IF OVERFLOW IS CONCEIVABLE, SCALE THE VECTOR DOWN.                
C     THIS APPROACH IS USED TO AVOID MACHINE-DEPENDENT AND SYSTEM-      
C     DEPENDENT CALLS TO OVERFLOW ROUTINES.                             
      IF (ABS(ELIM1).GT.DELBIG) GO TO 440
      TEMP=B(L,2)
      IF (ABS(B(L,2)).LT.DELTA) TEMP=DELTA
      VECT(L,I)=ELIM1/TEMP
      GO TO 460
C     VECTOR IS TOO BIG.  SCALE IT DOWN.                                
  440 DO 450 K=1,N
  450 VECT(K,I)=VECT(K,I)/DELBIG
      GO TO 430
  460 CONTINUE
      GO TO (470,530), ITER
C     SECOND ITERATION.  (BOTH ITERATIONS FOR REPEATED-ROOT VECTORS).   
  470 ITER=ITER+1
  480 ELIM1=VECT(1,I)
      DO 500 J=1,NM1
      IF (B(J,2).EQ.B(J,1)) GO TO 490
C     CASE ONE.                                                         
      VECT(J,I)=ELIM1
      ELIM1=VECT(J+1,I)-ELIM1*B(J,5)
      GO TO 500
C     CASE TWO.                                                         
  490 VECT(J,I)=VECT(J+1,I)
      ELIM1=ELIM1-VECT(J+1,I)*TEMP
  500 CONTINUE
      VECT(N,I)=ELIM1
      GO TO 420
C     PRODUCE A RANDOM VECTOR                                           
  510 CONTINUE
      DO 520 J=1,N
C     GENERATE PSEUDORANDOM NUMBERS WITH UNIFORM DISTRIBUTION IN (-1,1).
C     THIS RANDOM NUMBER SCHEME IS OF THE FORM...                       
C     RAND1 = AMOD((2**12+3)*RAND1,2**23)                               
C     IT HAS A PERIOD OF 2**21 NUMBERS.                                 
      AAAA=4099.0D+00*RAND1
      BBBB=8388608.0D+00
      IIII=AAAA/BBBB
      CCCC=IIII
      RAND1=AAAA-CCCC*BBBB
  520 VECT(J,I)=RAND1/4194304.0D+00-1.0D+00
      GO TO 420
C                                                                       
C     ORTHOGONALIZE THIS REPEATED-ROOT VECTOR TO OTHERS WITH THIS ROOT. 
  530 IF (IA.EQ.0) GO TO 570
      DO 560 J1=1,IA
      K=I-J1
      TEMP=0.D+00
      DO 540 J=1,N
  540 TEMP=TEMP+VECT(J,I)*VECT(J,K)
      DO 550 J=1,N
  550 VECT(J,I)=VECT(J,I)-TEMP*VECT(J,K)
  560 CONTINUE
  570 GO TO (480,580), ITER
C     NORMALIZE THE VECTOR                                              
  580 ELIM1=0.D+00
      DO 590 J=1,N
  590 ELIM1=MAX(ABS(VECT(J,I)),ELIM1)
      TEMP=0.D+00
      DO 600 J=1,N
      ELIM2=VECT(J,I)/ELIM1
  600 TEMP=TEMP+ELIM2**2
      TEMP=1.0D+00/(SQRT(TEMP)*ELIM1)
      DO 610 J=1,N
      VECT(J,I)=VECT(J,I)*TEMP
      IF (ABS(VECT(J,I)).LT.DEL1) VECT(J,I)=0.D+00
  610 CONTINUE
  620 CONTINUE
C                                                                       
C      SIMVEC SECTION.                                                  
C      ROTATE CODIAGONAL VECTORS INTO VECTORS OF ORIGINAL ARRAY         
C      LOOP OVER ALL THE TRANSFORMATION VECTORS                         
      IF (NM2.EQ.0) GO TO 680
      JUMP=NSIZE-(N+1)
      IM=NM1
      DO 670 I=1,NM2
      J1=JUMP
C      MOVE A TRANSFORMATION VECTOR OUT INTO BETTER INDEXING POSITION.  
      DO 630 J=IM,N
      B(J,2)=A(J1)
  630 J1=J1+J
C      MODIFY ALL REQUESTED VECTORS.                                    
      DO 660 K=1,NROOT
      TEMP=0.D+00
C      FORM SCALAR PRODUCT OF TRANSFORMATION VECTOR WITH EIGENVECTOR    
      DO 640 J=IM,N
  640 TEMP=TEMP+B(J,2)*VECT(J,K)
      TEMP=TEMP+TEMP
      DO 650 J=IM,N
  650 VECT(J,K)=VECT(J,K)-TEMP*B(J,2)
  660 CONTINUE
      JUMP=JUMP-IM
  670 IM=IM-1
  680 CONTINUE
C      RESTORE ROOTS TO THEIR PROPER SIZE.                              
      DO 690 I=1,NROOT
  690 ROOT(I)=ROOT(I)*ANORM
  700 CONTINUE
      RETURN
      END
 

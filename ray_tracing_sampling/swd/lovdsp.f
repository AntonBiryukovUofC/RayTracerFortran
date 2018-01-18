C                                                                       LOV02690
C     LOVE WAVE DISPERSION INTERPOLATION                                LOV02700
C                                                                       LOV02710
      SUBROUTINE  LOVDSP(DIFFEQ,H,RHO,VS,AX,L,W,CMN,CMX,                LOV02720
     *                   DC,TOL,ITR,IA,C,U,EK,Y0,IER)                   LOV02730
C                                                                       LOV02740
C     INPUT                                                             LOV02750
C       DIFFEQ : SUBROUTINE TO INTEGRATE EQ. OF MOTION                  LOV02760
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      LOV02770
C       RHO  : DENSITY                                                  LOV02780
C       VS   : SHEAR WAVE VELOCITY                                      LOV02790
C       AX   : ANISOTROPY FACTOR XI                                     LOV02800
C       L    : NUMBER OF LAYERS INCLUDING THE BOTTOM                    LOV02810
C              HALF-SPACE                                               LOV02820
C       W    : ANGULAR FREQUENCY                                        LOV02830
C       CMN  : LOWER LIMIT OF PHASE VELOCITY                            LOV02840
C       CMX  : UPPER LIMIT OF PHASE VELOCITY                            LOV02850
C       DC   : INCREMENT   OF PHASE VELOCITY                            LOV02860
C       TOL  : RELATIVE ACCURACY OF PHASE VELOCITY                      LOV02870
C       ITR  : MAXIMUM NUMBER OF ITERATIONS                             LOV02880
C       IA   : = 0 ;   ISOTROPIC MODEL                                  LOV02890
C              = 1 ; ANISOTROPIC MODEL                                  LOV02900
C     OUTPUT                                                            LOV02910
C       C    : PHASE VELOCITY                                           LOV02920
C       U    : GROUP VELOCITY BY DIFFERENTIATION                        LOV02930
C       EK   : KINETIC ENERGY 2*K**2*I3 BY DIFFERENTIATION              LOV02940
C       Y0   : SURFACE VALUE OF Y, C*DY/DC, AND W*DY/DW                 LOV02950
C              Y0(1) ; Y1                                               LOV02960
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  LOV02970
C       IER  : RETURN CODE                                              LOV02980
C              < 0 ; INPUT ERROR                                        LOV02990
C              = 0 ; NO ERROR                                           LOV03000
C              = 1 ; SLOW CONVERGENCE                                   LOV03010
C              = 2 ; ROOT NOT FOUND                                     LOV03020
C                                                                       LOV03030
C     SUBROUTINE : DIFFEQ (LOVMRX OR LOVRKG)                            LOV03040
C                                                                       LOV03050
C     DISPER-80,  VER-86                                                LOV03060
C                                                                       LOV03070
C     M. SAITO  23/VI/79                                                LOV03080
C     REVISED   24/IX/85                                                LOV03090
C                                                                       LOV03100
      DIMENSION  H(L),RHO(L),VS(L),AX(L),Y0(6)                          LOV03110
C                                                                       LOV03120
C     INITIALIZATION                                                    LOV03130
C                                                                       LOV03140
      IF( L.LE.0 .OR. W.LE.0 .OR. DC.EQ.0 )  GO TO  90                  LOV03150
      ONE = 1                                                           LOV03160
C                                                                       LOV03170
C     MACHINE EPSILON                                                   LOV03180
C                                                                       LOV03190
      EPS = 1                                                           LOV03200
    1 IF( (1+EPS).LE.1 )  GO TO  2                                      LOV03210
        EPS = EPS/2                                                     LOV03220
        GO TO  1                                                        LOV03230
    2 EPS = EPS*2                                                       LOV03240
C                                                                       LOV03250
      TOL1 = TOL                                                        LOV03260
      IF( TOL1.GT.0 )  TOL1 = MAX(TOL1,EPS)                             LOV03270
C     WRITE(6,3)  W                                                     LOV03280
    3   FORMAT(/7X,'W',8X,'LOVE WAVE'/                                  LOV03290
     *         1PE18.6//7X,'C',17X,'Y2')                                LOV03300
C    *         1PD18.6//7X,'C',17X,'Y2')                                LOV03310
      C3 = CMN                                                          LOV03320
      IF( C3.LE.0 )  GO TO  92                                          LOV03330
      CALL  DIFFEQ(H,RHO,VS,AX,L,W,C3,IA,1,U,EK,Y0,IER)                 LOV03340
      IF( IER.NE.0 )  RETURN                                            LOV03350
      F3 = Y0(2)                                                        LOV03360
C     WRITE(6,4)  C3,F3                                                 LOV03370
    4   FORMAT(1P2E18.6)                                                LOV03380
C   4   FORMAT(1P2D18.6)                                                LOV03390
      IF( F3.EQ.0 .AND. TOL1.GT.0 )  GO TO  9                           LOV03400
C                                                                       LOV03410
C     FIND A ZERO-CROSS                                                 LOV03420
C                                                                       LOV03430
      KX = (CMX - CMN)/DC + 0.5                                         LOV03440
      KX = MAX(1,KX)                                                    LOV03450
      DO  5  K=1,KX                                                     LOV03460
        CC = CMN + DC*K                                                 LOV03470
        C1 = C3                                                         LOV03480
        F1 = F3                                                         LOV03490
        C3 = CC                                                         LOV03500
        IF( C3.LE.0 )  GO TO  92                                        LOV03510
        CALL  DIFFEQ(H,RHO,VS,AX,L,W,C3,IA,1,U,EK,Y0,IER)               LOV03520
        IF( IER.NE.0 )  RETURN                                          LOV03530
        F3 = Y0(2)                                                      LOV03540
C       WRITE(6,4)  C3,F3                                               LOV03550
        IF(  TOL1.LE.0 )  GO TO  5                                      LOV03560
        IF( F3*SIGN(ONE,F1).LE.0 )  GO TO  6                            LOV03570
    5 CONTINUE                                                          LOV03580
C                                                                       LOV03590
      IER = 2                                                           LOV03600
      IF( TOL1.LE.0 )  RETURN                                           LOV03610
      GO TO  94                                                         LOV03620
C                                                                       LOV03630
C     INTERPOLATION                                                     LOV03640
C                                                                       LOV03650
    6 IF( F3.EQ.0 )  GO TO  9                                           LOV03660
      C2 = C3                                                           LOV03670
      F2 = F3                                                           LOV03680
      E  = C1 - C2                                                      LOV03690
      D  = E/2                                                          LOV03700
      C3 = C2 + D                                                       LOV03710
      KX = MAX(1,ITR)                                                   LOV03720
C                                                                       LOV03730
      DO  7  K=1,KX                                                     LOV03740
        CALL  DIFFEQ(H,RHO,VS,AX,L,W,C3,IA,1,U,EK,Y0,IER)               LOV03750
        F3 = Y0(2)                                                      LOV03760
C       WRITE(6,4)  C3,F3                                               LOV03770
        IF( F3*SIGN(ONE,F2).GT.0 )  THEN                                LOV03780
          FF = C1                                                       LOV03790
          C1 = C2                                                       LOV03800
          C2 = FF                                                       LOV03810
          FF = F1                                                       LOV03820
          F1 = F2                                                       LOV03830
          F2 = FF                                                       LOV03840
        END IF                                                          LOV03850
        IF( ABS(F3).GT.ABS(F2) )  THEN                                  LOV03860
          FF = C2                                                       LOV03870
          C2 = C3                                                       LOV03880
          C3 = FF                                                       LOV03890
          FF = F2                                                       LOV03900
          F2 = F3                                                       LOV03910
          F3 = FF                                                       LOV03920
        END IF                                                          LOV03930
        E = C2 - C3                                                     LOV03940
        IF( F3.EQ.0 )  GO TO  9                                         LOV03950
        TOLC = C3*TOL1                                                  LOV03960
        DD = D                                                          LOV03970
        F32 = F3/F2                                                     LOV03980
        F31 = F3/F1                                                     LOV03990
        F21 = F2/F1                                                     LOV04000
        Q = F32*(E*(1 - F31) + F21*(F31 - F21)*(C1 - C3))               LOV04010
        S = (F21 - 1)*(F32 - 1)*(F31 - 1)                               LOV04020
C                                                                       LOV04030
C       TEST RANGE                                                      LOV04040
C                                                                       LOV04050
        IF( Q.LT.0 )  S =-S                                             LOV04060
        Q = ABS(Q)                                                      LOV04070
        IF( Q.GE.(E*S-ABS(TOLC*S)) )  THEN                              LOV04080
C                                                                       LOV04090
C       LINEAR INTERPOLATION                                            LOV04100
C                                                                       LOV04110
          D = E*F32/(F32 - 1)                                           LOV04120
C                                                                       LOV04130
C       INVERSE QUADRATIC INTERPOLATION                                 LOV04140
C                                                                       LOV04150
        ELSE                                                            LOV04160
          D = Q/S                                                       LOV04170
        END IF                                                          LOV04180
C                                                                       LOV04190
C       TEST CONVERGENCE                                                LOV04200
C                                                                       LOV04210
        C1 = C2                                                         LOV04220
        F1 = F2                                                         LOV04230
        C2 = C3                                                         LOV04240
        F2 = F3                                                         LOV04250
        C3 = C2 + D                                                     LOV04260
        IF( ABS(E).LE.TOLC )  GO TO  9                                  LOV04270
        IF( ABS(D).LE.TOLC )  THEN                                      LOV04280
C                                                                       LOV04290
C       BISECTION                                                       LOV04300
C                                                                       LOV04310
          IF( ABS(DD).LE.TOLC )  THEN                                   LOV04320
            D  = E/2                                                    LOV04330
            C3 = C2 + D                                                 LOV04340
          ELSE                                                          LOV04350
            C3 = C2 + SIGN(TOLC,D)                                      LOV04360
          END IF                                                        LOV04370
        END IF                                                          LOV04380
C                                                                       LOV04390
    7 CONTINUE                                                          LOV04400
C                                                                       LOV04410
C     SLOW CONVERGENCE                                                  LOV04420
C                                                                       LOV04430
C     WRITE(6,8)  KX                                                    LOV04440
    8   FORMAT(20X,5('?'),3X,'(LOVDSP)   SLOW CONV. ',                  LOV04450
     *         'AFTER',I5,' ITERATIONS',3X,5('?'))                      LOV04460
      IER = 1                                                           LOV04470
C                                                                       LOV04480
C     ROOT IS FOUND                                                     LOV04490
C                                                                       LOV04500
    9 CALL  DIFFEQ(H,RHO,VS,AX,L,W,C3,IA,3,U,EK,Y0,IER1)                LOV04510
C     WRITE(6,4)  C3,Y0(2)                                              LOV04520
      C = C3                                                            LOV04530
      RETURN                                                            LOV04540
C                                                                       LOV04550
C     INPUT ERROR                                                       LOV04560
   90 CONTINUE                                                          LOV04570
C  90 WRITE(6,91)  L,W,DC                                               LOV04580
   91   FORMAT(20X,5('?'),3X,'(LOVDSP)   INPUT ERROR',                  LOV04590
     *         'L  =',I5,3X,'W =',1PE13.6,3X,                           LOV04600
     1         'DC =',E13.6,3X,5('?'))                                  LOV04610
      IER =-1                                                           LOV04620
      RETURN                                                            LOV04630
   92 CONTINUE                                                          LOV04640
C  92 WRITE(6,93)  C3                                                   LOV04650
   93   FORMAT(20X,5('?'),3X,'(LOVDSP)   INPUT ERROR',                  LOV04660
     *         'C  =',1PE13.6,3X,5('?'))                                LOV04670
      IER =-1                                                           LOV04680
      RETURN                                                            LOV04690
C                                                                       LOV04700
C     NO ROOT                                                           LOV04710
   94 CONTINUE                                                          LOV04720
C  94 WRITE(6,95)                                                       LOV04730
   95   FORMAT(20X,5('?'),3X,'(LOVDSP)',3X,                             LOV04740
     *         'ROOT NOT FOUND',3X,5('?'))                              LOV04750
      RETURN                                                            LOV04760
      END                                                               LOV04770

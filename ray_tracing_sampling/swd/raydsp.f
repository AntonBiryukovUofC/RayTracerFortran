C                                                                       RAY02740
C     RAYLEIGH WAVE EIGENVALUE (R-K-G VERSION)                          RAY02750
C                                                                       RAY02760
      SUBROUTINE  RAYDSP(DIFFEQ,H,RHO,VP,VS,AP,AE,L,                    RAY02770
     *                   W,CMN,CMX,DC,TOL,ITR,IA,C,U,EK,                RAY02780
     1                   Y0,YIJ,IER)                                    RAY02790
C                                                                       RAY02800
C     INPUT                                                             RAY02810
C       DIFFEQ : SUBROUTINE TO INTEGRATE EQ. OF MOTION                  RAY02820
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      RAY02830
C       RHO  : DENSITY                                                  RAY02840
C       VP   : COMPRESSIONAL WAVE VELOCITY                              RAY02850
C       VS   : SHEAR WAVE VELOCITY                                      RAY02860
C       AP   : ANISOTROPY FACTOR PHI                                    RAY02870
C       AE   : ANISOTROPY FACTOR ETA                                    RAY02880
C       L    : NUMBER OF LAYERS INCLUDING THE BOTTOM                    RAY02890
C              HALF-SPACE                                               RAY02900
C       W    : ANGULAR FREQUENCY                                        RAY02910
C       CMN  : LOWER LIMIT OF PHASE VELOCITY                            RAY02920
C       CMX  : UPPER LIMIT OF PHASE VELOCITY                            RAY02930
C       DC   : INCREMENT   OF PHASE VELOCITY                            RAY02940
C       TOL  : RELATIVE ACCURACY OF PHASE VELOCITY                      RAY02950
C       ITR  : MAXIMUM NUMBER OF ITERATIONS                             RAY02960
C       IA   : = 0 ;   ISOTROPIC MODEL                                  RAY02970
C              = 1 ; ANISOTROPIC MODEL                                  RAY02980
C     OUTPUT                                                            RAY02990
C       C    : PHASE VELOCITY                                           RAY03000
C       U    : GROUP VELOCITY BY DIFFERENTIATION                        RAY03010
C       EKD  : ENERGY INTEGRAL BY DIFFERENTIATION                       RAY03020
C              2*K**2*I3                                                RAY03030
C       Y0   : SURFACE VALUES OF EIGENFUNCTION                          RAY03040
C              Y0(1) ; Y1 (SCALE FACTOR)                                RAY03050
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  RAY03060
C              Y0(3) ; Y3/Y1                                            RAY03070
C       YIJ  : SURFACE VALUES OF YIJ (COMPOUNDS),                       RAY03080
C              C*DYIJ/DC, AND W*DYIJ/DW                                 RAY03090
C       IER  : RETURN CODE                                              RAY03100
C              < 0 ; INPUT ERROR                                        RAY03110
C              = 0 ; NO ERROR                                           RAY03120
C              = 1 ; SLOW CONVERGENCE                                   RAY03130
C              = 2 ; ROOT NOT FOUND                                     RAY03140
C                                                                       RAY03150
C     SUBROUTINE : DIFFEQ (RAYMRX OR RAYRKG)                            RAY03160
C                                                                       RAY03170
C     DISPER-80,  VER-86                                                RAY03180
C                                                                       RAY03190
C     M. SAITO  23/VI/79                                                RAY03200
C     REVISED   25/IX/85                                                RAY03210
C                                                                       RAY03220
      DIMENSION  H(L),RHO(L),VP(L),VS(L),AP(L),AE(L),                   RAY03230
     *           Y0(3),YIJ(15)                                          RAY03240
C                                                                       RAY03250
C     INITIALIZATION                                                    RAY03260
C                                                                       RAY03270
      IF( L.LE.2 .OR. W.LE.0 .OR. DC.EQ.0 )  GO TO  90                  RAY03280
      ONE = 1                                                           RAY03290
C                                                                       RAY03300
C     MACHINE EPSILON                                                   RAY03310
C                                                                       RAY03320
      EPS = 1                                                           RAY03330
    1 IF( (1+EPS).LE.1 )  GO TO  2                                      RAY03340
        EPS = EPS/2                                                     RAY03350
        GO TO  1                                                        RAY03360
    2 EPS = EPS/2                                                       RAY03370
C                                                                       RAY03380
      TOL1 = TOL                                                        RAY03390
      IF( TOL1.GT.0 )  TOL1 = MAX(TOL1,EPS)                             RAY03400
      IER1 = 0                                                          RAY03410
C     WRITE(6,3)  W                                                     RAY03420
    3   FORMAT(/7X,'W',15X,'RAYLEIGH WAVE'/                             RAY03430
     *         1PE18.6//7X,'C',17X,'Y2',16X,'Y3')                       RAY03440
C    *         1PD18.6//7X,'C',17X,'Y2',16X,'Y3')                       RAY03450
      C3 = CMN                                                          RAY03460
      IF( C3.LE.0 )  GO TO  92                                          RAY03470
      CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                     RAY03480
     *             EK,Y0,YIJ,IER)                                       RAY03490
      IF( IER.NE.0 )  RETURN                                            RAY03500
      F3 = Y0(2)                                                        RAY03510
C     WRITE(6,4)  C3,F3,Y0(3)                                           RAY03520
    4   FORMAT(1P3E18.6)                                                RAY03530
C   4   FORMAT(1P3D18.6)                                                RAY03540
      IF( F3.EQ.0 .AND. TOL1.GT.0 )  GO TO  9                           RAY03550
C                                                                       RAY03560
C     FIND A ZERO-CROSS                                                 RAY03570
C                                                                       RAY03580
      KX = (CMX - CMN)/DC + 0.5                                         RAY03590
      KX = MAX(1,KX)                                                    RAY03600
      DO  5  K=1,KX                                                     RAY03610
        CC = CMN + K*DC                                                 RAY03620
        C1 = C3                                                         RAY03630
        F1 = F3                                                         RAY03640
        C3 = CC                                                         RAY03650
        IF( C3.LE.0 )  GO TO  92                                        RAY03660
        CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                   RAY03670
     *               EK,Y0,YIJ,IER)                                     RAY03680
        IF( IER.NE.0 )  RETURN                                          RAY03690
        F3 = Y0(2)                                                      RAY03700
C       WRITE(6,4)  C3,F3,Y0(3)                                         RAY03710
        IF(  TOL1.LE.0 )  GO TO  5                                      RAY03720
        IF( F3*SIGN(ONE,F1).LE.0 )  GO TO  6                            RAY03730
    5 CONTINUE                                                          RAY03740
C                                                                       RAY03750
      IER = 2                                                           RAY03760
      IF( TOL1.LE.0 )  RETURN                                           RAY03770
      GO TO  94                                                         RAY03780
C                                                                       RAY03790
C     INTERPOLATION                                                     RAY03800
C                                                                       RAY03810
    6 IF( F3.EQ.0 )  GO TO  9                                           RAY03820
      C2 = C3                                                           RAY03830
      F2 = F3                                                           RAY03840
      E  = C1 - C2                                                      RAY03850
      D  = E/2                                                          RAY03860
      C3 = C2 + D                                                       RAY03870
      KX = MAX(1,ITR)                                                   RAY03880
C                                                                       RAY03890
      DO  7  K=1,KX                                                     RAY03900
        CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,1,U,                   RAY03910
     *               EK,Y0,YIJ,IER)                                     RAY03920
        F3 = Y0(2)                                                      RAY03930
C       WRITE(6,4) C3,F3,Y0(3)                                          RAY03940
        IF( F3*SIGN(ONE,F2).GT.0 )  THEN                                RAY03950
          FF = C1                                                       RAY03960
          C1 = C2                                                       RAY03970
          C2 = FF                                                       RAY03980
          FF = F1                                                       RAY03990
          F1 = F2                                                       RAY04000
          F2 = FF                                                       RAY04010
        END IF                                                          RAY04020
        IF( ABS(F3).GT.ABS(F2) )  THEN                                  RAY04030
          FF = C2                                                       RAY04040
          C2 = C3                                                       RAY04050
          C3 = FF                                                       RAY04060
          FF = F2                                                       RAY04070
          F2 = F3                                                       RAY04080
          F3 = FF                                                       RAY04090
        END IF                                                          RAY04100
        E = C2 - C3                                                     RAY04110
        IF( F3.EQ.0 )  GO TO  9                                         RAY04120
        TOLC = C3*TOL1                                                  RAY04130
        DD = D                                                          RAY04140
        F32 = F3/F2                                                     RAY04150
        F31 = F3/F1                                                     RAY04160
        F21 = F2/F1                                                     RAY04170
        Q = F32*(E*(1 - F31) + F21*(F31 - F21)*(C1 - C3))               RAY04180
        S = (F21 - 1)*(F32 - 1)*(F31 - 1)                               RAY04190
C                                                                       RAY04200
C       TEST RANGE                                                      RAY04210
C                                                                       RAY04220
        IF( Q.LT.0 )  S =-S                                             RAY04230
        Q = ABS(Q)                                                      RAY04240
        IF( Q.GE.(E*S-ABS(TOLC*S)) )  THEN                              RAY04250
C                                                                       RAY04260
C       LINEAR INTERPOLATION                                            RAY04270
C                                                                       RAY04280
          D = E*F32/(F32 - 1)                                           RAY04290
C                                                                       RAY04300
C       INVERSE QUADRATIC INTERPOLATION                                 RAY04310
C                                                                       RAY04320
        ELSE                                                            RAY04330
          D = Q/S                                                       RAY04340
        END IF                                                          RAY04350
C                                                                       RAY04360
C       TEST CONVERGENCE                                                RAY04370
C                                                                       RAY04380
        C1 = C2                                                         RAY04390
        F1 = F2                                                         RAY04400
        C2 = C3                                                         RAY04410
        F2 = F3                                                         RAY04420
        C3 = C2 + D                                                     RAY04430
        IF( ABS(E).LE.TOLC )  GO TO  9                                  RAY04440
        IF( ABS(D).LE.TOLC )  THEN                                      RAY04450
C                                                                       RAY04460
C       BISECTION                                                       RAY04470
C                                                                       RAY04480
          IF( ABS(DD).LE.TOLC )  THEN                                   RAY04490
            D = E/2                                                     RAY04500
            C3 = C2 + D                                                 RAY04510
          ELSE                                                          RAY04520
            C3 = C2 + SIGN(TOLC,D)                                      RAY04530
          END IF                                                        RAY04540
        END IF                                                          RAY04550
C                                                                       RAY04560
    7 CONTINUE                                                          RAY04570
C                                                                       RAY04580
C     SLOW CONVERGENCE                                                  RAY04590
C                                                                       RAY04600
C     WRITE(6,8)  KX                                                    RAY04610
    8   FORMAT(20X,5('?'),3X,'(RAYDSP)   SLOW CONV. ',                  RAY04620
     *         'AFTER',I5,' ITERATIONS',3X,5('?'))                      RAY04630
      IER = 1                                                           RAY04640
C                                                                       RAY04650
C     ROOT IS FOUND                                                     RAY04660
C                                                                       RAY04670
    9 CALL  DIFFEQ(H,RHO,VP,VS,AP,AE,L,W,C3,IA,3,U,                     RAY04680
     *             EK,Y0,YIJ,IER1)                                      RAY04690
C     WRITE(6,4)  C3,Y0(2),Y0(3)                                        RAY04700
      C   = C3                                                          RAY04710
      RETURN                                                            RAY04720
C                                                                       RAY04730
C     INPUT ERROR                                                       RAY04740
90    CONTINUE                                                          RAY04750
C  90 WRITE(6,91)  L,W,DC                                               RAY04760
   91   FORMAT(20X,5('?'),3X,'(RAYDSP)   INPUT ERROR',3X,               RAY04770
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          RAY04780
     1         'DC =',E13.6,3X,5('?'))                                  RAY04790
      IER =-1                                                           RAY04800
      RETURN                                                            RAY04810
92    CONTINUE                                                          RAY04820
C  92 WRITE(6,93)  C3                                                   RAY04830
   93   FORMAT(20X,5('?'),3X,'(RAYDSP)   INPUT ERROR',3X,               RAY04840
     *         'C  =',1PE13.6,3X,5('?'))                                RAY04850
      IER =-1                                                           RAY04860
      RETURN                                                            RAY04870
C                                                                       RAY04880
C     NO ROOT                                                           RAY04890
94    CONTINUE                                                          RAY04900
C  94 WRITE(6,95)                                                       RAY04910
   95   FORMAT(20X,5('?'),3X,'(RAYDSP)',3X,                             RAY04920
     *         'ROOT NOT FOUND',3X,5('?'))                              RAY04930
      RETURN                                                            RAY04940
      END                                                               RAY04950

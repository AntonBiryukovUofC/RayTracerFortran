C                                                                       LOV06850
C     LOVE WAVE MATRIX METHOD INTEGRATION                               LOV06860
C                                                                       LOV06870
      SUBROUTINE  LOVMRX(H,RHO,VS,AX,L,W,C,                             LOV06880
     *                   IA,IG,U,EK,Y0,IER)                             LOV06890
C                                                                       LOV06900
C     INPUT                                                             LOV06910
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      LOV06920
C       RHO  : DENSITY                                                  LOV06930
C       VS   : SHEAR WAVE VELOCITY                                      LOV06940
C       AX   : ANISOTROPY FACTOR XI                                     LOV06950
C       L    : NUMBER OF LAYERS INCLUDING THE BOTTOM                    LOV06960
C              HALF-SPACE                                               LOV06970
C       W    : ANGULAR FREQUENCY                                        LOV06980
C       C    : PHASE VELOCITY                                           LOV06990
C       IA   : = 0 ;   ISOTROPIC MODEL                                  LOV07000
C              = 1 ; ANISOTROPIC MODEL                                  LOV07010
C       IG   : = 1 ; TO COMPUTE Y ONLY                                  LOV07020
C              = 2 ; TO COMPUTE Y AND C*DY/DC                           LOV07030
C              = 3 ; TO COMPUTE Y, C*DY/DC AND W*DY/DW                  LOV07040
C     OUTPUT                                                            LOV07050
C       U    : GROUP VELOCITY BY DIFFERENTIATION (IG=3)                 LOV07060
C       EK   : ENERGY INTEGRAL 2*K**2*I3 BY                             LOV07070
C              DIFFERENTIATION  (IG>=2)                                 LOV07080
C       Y0   : SURFACE VALUE OF Y                                       LOV07090
C              SURFACE VALUE OF C*DY/DC (WHEN IG>=2)                    LOV07100
C              SURFACE VALUE OF W*DY/DW (WHEN IG=3)                     LOV07110
C              Y0(1) ; Y1                                               LOV07120
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  LOV07130
C                 Y2 HERE IS Y2/K IN THE CONVENTIONAL                   LOV07140
C                 DEFINITION                                            LOV07150
C       IER  : RETURN CODE                                              LOV07160
C              < 0 ; INPUT ERROR                                        LOV07170
C              = 0 ; NO ERROR                                           LOV07180
C                                                                       LOV07190
C     DISPER-80,  VER-86                                                LOV07200
C                                                                       LOV07210
C     M. SAITO  23/VI/79                                                LOV07220
C     REVISED    3/XII/86                                               LOV07230
C                                                                       LOV07240
      DIMENSION  H(L),RHO(L),VS(L),AX(L),Y0(6)                          LOV07250
      DATA  EPS/1.E-10/,BIG/1.E+10/                                     LOV07260
C                                                                       LOV07270
C     DEFINE  SH0(X**2)= SINH(X)/X AND                                  LOV07280
C             SH1(X**2)=(COSH(X)-SINH(X)/X)/X**2                        LOV07290
C                                                                       LOV07300
      SH0(X) = 0.9999997 + X*(0.1666667                                 LOV07310
     *       + X*(0.0083361 + X*0.0001984))                             LOV07320
      SH1(X) = 0.3333333 + X*(0.0333333                                 LOV07330
     *       + X*(0.0011907 + X*0.0000220))                             LOV07340
C                                                                       LOV07350
C     FOR DOUBLE PRECISION USE                                          LOV07360
C                                                                       LOV07370
C     SH0(X) = 1.0D0 + X*(1.6666 66666 66666 7D-1                       LOV07380
C    *       + X*(8.33 33333 33334 0D-3                                 LOV07390
C    1       + X*( 1.9 84126 98412 7D-4                                 LOV07400
C    2       + X*(2.7557 31918 9D-6 + X*(2.50 12108 4D-8                LOV07410
C    3       + X*(1.60596 1D-10 + X*7.64 7D-13))))))                    LOV07420
C     SH1(X) = 3.3333 33333 33333 3D-1                                  LOV07430
C    *       + X*(3.333 33333 33333 3D-2                                LOV07440
C    1       + X*(1.19 04761 90476 2D-3                                 LOV07450
C    2       + X*(2.20458 55379 2D-5                                    LOV07460
C    3       + X*(2.505 21083 7D-7 + X*(1.9 27085 3D-9                  LOV07470
C    4       + X*(1.0706 3D-11 + X*4.50D-14))))))                       LOV07480
C                                                                       LOV07490
C     INITIAL VALUE                                                     LOV07500
C                                                                       LOV07510
      IF( L.LE.1 .OR. W.LE.0 .OR. C.LE.0 )  GO TO  90                   LOV07520
      I   = L                                                           LOV07530
      RO  = RHO(I)                                                      LOV07540
      IF( RO.LE.0 )  GO TO  92                                          LOV07550
      DO  1  J=2,L                                                      LOV07560
        IF( VS(I-1).GT.0 )  GO TO  2                                    LOV07570
        I = I - 1                                                       LOV07580
    1 CONTINUE                                                          LOV07590
      GO TO  92                                                         LOV07600
C                                                                       LOV07610
    2 IER = 0                                                           LOV07620
      IGG = MAX(1,MIN(3,IG))                                            LOV07630
      WN  = W/C                                                         LOV07640
      CC  = C*C                                                         LOV07650
      ROC = RO*CC                                                       LOV07660
      Y1 = EPS                                                          LOV07670
      Y2 = 0                                                            LOV07680
      Y3 = 0                                                            LOV07690
      Y4 = 0                                                            LOV07700
      Y5 = 0                                                            LOV07710
      Y6 = 0                                                            LOV07720
      EL  = RO*VS(I)**2                                                 LOV07730
      IF( EL.GT.0 )  THEN                                               LOV07740
        XI  = 1                                                         LOV07750
        IF( IA.GT.0 )  XI = AX(I)                                       LOV07760
        RBB = XI - (C/VS(I))**2                                         LOV07770
C       WRITE (*,*),"XI C VS(I) RBB= ",RBB,XI,C,VS(I)                   LOV07775
        IF( RBB.LE.0 )  GO TO  92                                       LOV07780
        RB  = SQRT(RBB)                                                 LOV07790
        Y2  = Y1*EL*RB                                                  LOV07800
        IF( IGG.GE.2 )  Y4 =-Y1*ROC/RB                                  LOV07810
      ENDIF                                                             LOV07820
C                                                                       LOV07830
C     INTEGRATE UPWARD                                                  LOV07840
C                                                                       LOV07850
      DO  3  II=2,L                                                     LOV07860
        I   = I - 1                                                     LOV07870
        RO  = RHO(I)                                                    LOV07880
        ROC = RO*CC                                                     LOV07890
        XI  = 1                                                         LOV07900
        IF( IA.GT.0 )  XI = AX(I)                                       LOV07910
        EL  = RO*VS(I)**2                                               LOV07920
        IF( EL.LE.0 .OR. H(I).LE.0 )  GO TO  92                         LOV07930
        RBB = XI - (C/VS(I))**2                                         LOV07940
        HK  = H(I)*WN                                                   LOV07950
        X   = RBB*HK**2                                                 LOV07960
C                                                                       LOV07970
C       SINH(X)/X                                                       LOV07980
C                                                                       LOV07990
        AA = ABS(X)                                                     LOV08000
        IF( AA.LE.1 )  THEN                                             LOV08010
          SHB = SH0(X)                                                  LOV08020
          CHB = 1 + X*SH0(X/4)**2/2                                     LOV08030
          IF( IGG.GE.2 )  DHB = HK**2*SH1(X)                            LOV08040
        ELSE                                                            LOV08050
          AA = SQRT(AA)                                                 LOV08060
          IF( X.GE.0 )  THEN                                            LOV08070
            CHB = 1                                                     LOV08080
            SHB = TANH(AA)/AA                                           LOV08090
          ELSE                                                          LOV08100
            SHB = SIN(AA)/AA                                            LOV08110
            CHB = COS(AA)                                               LOV08120
          ENDIF                                                         LOV08130
          IF( IGG.GE.2 )  DHB = HK**2*(CHB - SHB)/X                     LOV08140
        ENDIF                                                           LOV08150
        SHB = HK*SHB                                                    LOV08160
C                                                                       LOV08170
C       LAYER MATRIX                                                    LOV08180
C                                                                       LOV08190
        B11 = CHB                                                       LOV08200
        B12 = SHB/EL                                                    LOV08210
        B21 = EL*RBB*SHB                                                LOV08220
C                                                                       LOV08230
        Z1  = Y1                                                        LOV08240
        Z2  = Y2                                                        LOV08250
        Y1  = B11*Z1 + B12*Z2                                           LOV08260
        Y2  = B21*Z1 + B11*Z2                                           LOV08270
C                                                                       LOV08280
        IF( IGG.GE.2 )  THEN                                            LOV08290
          C11 =-HK*SHB                                                  LOV08300
          C12 =-(HK*DHB + SHB)/EL                                       LOV08310
          C21 =-EL*HK*CHB - ROC*SHB                                     LOV08320
C                                                                       LOV08330
          Z3  = Y3                                                      LOV08340
          Z4  = Y4                                                      LOV08350
          Y3  = B11*Z3 + B12*Z4 + C11*Z1 + C12*Z2                       LOV08360
          Y4  = B21*Z3 + B11*Z4 + C21*Z1 + C11*Z2                       LOV08370
C                                                                       LOV08380
          IF( IGG.GT.2 )  THEN                                          LOV08390
            W11 = HK*RBB*SHB                                            LOV08400
            W12 = HK*CHB/EL                                             LOV08410
            W21 = HK*EL*RBB*CHB                                         LOV08420
C                                                                       LOV08430
            Z5  = Y5                                                    LOV08440
            Z6  = Y6                                                    LOV08450
            Y5  = B11*Z5 + B12*Z6 + W11*Z1 + W12*Z2                     LOV08460
            Y6  = B21*Z5 + B11*Z6 + W21*Z1 + W11*Z2                     LOV08470
          ENDIF                                                         LOV08480
        ENDIF                                                           LOV08490
C                                                                       LOV08500
C     NORMALIZATION                                                     LOV08510
C                                                                       LOV08520
        Z1 = MAX(ABS(Y1), ABS(Y2))                                      LOV08530
        IF( Z1.GT.BIG )  THEN                                           LOV08540
          Y1 = EPS*Y1                                                   LOV08550
          Y2 = EPS*Y2                                                   LOV08560
          IF( IGG.GE.2 )  THEN                                          LOV08570
            Y3 = EPS*Y3                                                 LOV08580
            Y4 = EPS*Y4                                                 LOV08590
            IF( IGG.GT.2 )  THEN                                        LOV08600
              Y5 = EPS*Y5                                               LOV08610
              Y6 = EPS*Y6                                               LOV08620
            ENDIF                                                       LOV08630
          ENDIF                                                         LOV08640
        ENDIF                                                           LOV08650
C                                                                       LOV08660
    3 CONTINUE                                                          LOV08670
C                                                                       LOV08680
C     NORMAL EXIT                                                       LOV08690
C                                                                       LOV08700
      Y0(1) = Y1                                                        LOV08710
      IF( ABS(Y2)*EPS.LE.ABS(Y1) )  THEN                                LOV08720
        Y0(2) = Y2/ABS(Y1)                                              LOV08730
      ELSE                                                              LOV08740
        Y0(2) = SIGN(BIG,Y2)                                            LOV08750
      ENDIF                                                             LOV08760
      IF( IGG.GE.2 )  THEN                                              LOV08770
        Y0(3) = Y3                                                      LOV08780
        Y0(4) = Y4                                                      LOV08790
        EK =-WN*Y4/Y1                                                   LOV08800
        IF( IGG.GT.2 )  THEN                                            LOV08810
          U =C*Y4/(Y4 + Y6)                                             LOV08820
          Y0(5) = Y5                                                    LOV08830
          Y0(6) = Y6                                                    LOV08840
        ENDIF                                                           LOV08850
      ENDIF                                                             LOV08860
C                                                                       LOV08870
      RETURN                                                            LOV08880
C                                                                       LOV08890
C     INPUT ERROR                                                       LOV08900
C                                                                       LOV08910
   90 WRITE(6,91)  L,W,C                                                LOV08920
   91   FORMAT(20X,5('?'),3X,'(LOVMRX)   INPUT ERROR',3X,               LOV08930
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          LOV08940
     1         'C  =',E13.6,3X,5('?'))                                  LOV08950
      IER =-1                                                           LOV08960
      RETURN                                                            LOV08970
C                                                                       LOV08980
   92 WRITE(6,93)  I                                                    LOV08990
   93   FORMAT(20X,5('?'),3X,'(LOVMRX)   INPUT ERROR',3X,               LOV09000
     *         'AT LAYER',I5,3X,5('?'))                                 LOV09010
      IER =-MAX(1,I)                                                    LOV09020
      RETURN                                                            LOV09030
      END                                                               LOV09040

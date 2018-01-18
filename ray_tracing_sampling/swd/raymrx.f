C                                                                       RAY08420
C     RAYLEIGH WAVE MATRIX METHOD INTEGRATION                           RAY08430
C                                                                       RAY08440
      SUBROUTINE  RAYMRX(H,RHO,VP,VS,AP,AE,L,W,C,                       RAY08450
     *                   IA,IG,U,EK,Y0,Y,IER)                           RAY08460
C                                                                       RAY08470
C     INPUT                                                             RAY08480
C       H    : LAYER THICKNESS.  H(L) IS ARBITRARY                      RAY08490
C       RHO  : DENSITY                                                  RAY08500
C       VP   : P WAVE VELOCITY                                          RAY08510
C       VS   : S WAVE VELOCITY                                          RAY08520
C       L    : NUMBER OF LAYERS INCLDING THE BOTTOM                     RAY08530
C              HALF-SPACE                                               RAY08540
C       W    : ANGULAR FREQUENCY                                        RAY08550
C       C    : PHASE VELOCITY                                           RAY08560
C       IG   : = 1 ; TO COMPUTE Y ONLY                                  RAY08570
C              = 2 ; TO COMPUTE Y AND C*CY/DC                           RAY08580
C              = 3 ; TO COMPUTE Y, C*DY/DC, AND W*DY/DW                 RAY08590
C     OUTPUT                                                            RAY08600
C       U    : GROUP VELOCITY BY DIFFERENTIATION (IG=3)                 RAY08610
C       EK   : ENERGY INTEGRAL 2*K**2*I3  (IG>=2)                       RAY08620
C       Y0   : SURFACE VALUES OF EIGENFUNCTIONS                         RAY08630
C              Y0(1) ; Y1 (SCALE FACTOR)                                RAY08640
C              Y0(2) ; Y2/ABS(Y1), DISPERSION FUNCTION                  RAY08650
C              Y0(3) ; Y3/Y1                                            RAY08660
C       Y    : SURFACE VALUES OF YIJ (COMPOUNDS),                       RAY08670
C              C*DYIJ/DC, AND W*DYIJ/DW                                 RAY08680
C              FOR SOLID LAYER                                          RAY08690
C                (IJ) = (12),(13),(14),(23),(24)                        RAY08700
C       IER  : RETURN CODE                                              RAY08710
C              < 0 ; INPUT ERROR                                        RAY08720
C              = 0 ; NO ERROR                                           RAY08730
C                                                                       RAY08740
C     ISOTROPIC MODEL ONLY.  AP, AE, IA ARE DUMMY.                      RAY08750
C                                                                       RAY08760
C     DISPER-80,  VER-86                                                RAY08770
C                                                                       RAY08780
C     M. SAITO  30/VII/79                                               RAY08790
C     REVISED   10/XII/86                                               RAY08800
C                                                                       RAY08810
      DIMENSION  H(L),RHO(L),VP(L),VS(L),Y0(3),Y(15),Z(15)              RAY08820
     *          ,AP(1),AE(1)                                            RAY08830
      DATA  EPS/1.E-10/,BIG/1.E+10/                                     RAY08840
C                                                                       RAY08850
C     DEFINE SINH(X)/X AND (COSH(X)-SINH(X)/X)/X**2                     RAY08860
C                                                                       RAY08870
      SH0(X) = 0.9999997 + X*(0.1666667 + X*(0.0083361                  RAY08880
     *       + X*0.0001984))                                            RAY08890
      SH1(X) = 0.3333333 + X*(0.0333333 + X*(0.0011907                  RAY08900
     *       + X*0.0000220))                                            RAY08910
C                                                                       RAY08920
C     FOR DOUBLE PRECISION USE                                          RAY08930
C                                                                       RAY08940
C     SH0(X) = 1.0D0 + X*(1.6666 66666 66666 7D-1                       RAY08950
C    *       + X*(8.33 33333 33334 0D-3                                 RAY08960
C    1       + X*(1.9 84126 98412 7D-4                                  RAY08970
C    2       + X*(2.7557 31918 9D-6 + X*(2.50 12108 4D-8                RAY08980
C    3       + X*(1.60596 1D-10 + X*7.64 7D-13))))))                    RAY08990
C     SH1(X) = 3.3333 33333 33333 3D-1                                  RAY09000
C    *       + X*(3.333 33333 33333 3D-2                                RAY09010
C    1       + X*(1.19 04761 90476 2D-3                                 RAY09020
C    2       + X*(2.20458 55379 2D-5                                    RAY09030
C    3       + X*(2.505 21083 7D-7 + X*(1.9 27085 3D-9                  RAY09040
C    4       + X*(1.0706 3D-11 + X*4.50D-14))))))                       RAY09050
C                                                                       RAY09060
C     INITIAL VALUE                                                     RAY09070
C                                                                       RAY09080
      IF( L.LE.0 .OR. W.LE.0 .OR. C.LE.0 )  GO TO  90                   RAY09090
      I   = L                                                           RAY09100
      RO  = RHO(I)                                                      RAY09110
      IF( RO.LE.0 .OR. C.GE.VP(I) )  GO TO  92                          RAY09120
      IER = 0                                                           RAY09130
      CC  = C*C                                                         RAY09140
      WN  = W/C                                                         RAY09150
      IGG = MAX(1,MIN(3,IG))                                            RAY09160
      ROC = RO*CC                                                       RAY09170
      SV  = VS(I)                                                       RAY09180
      CP  = C/VP(I)                                                     RAY09190
      RAA = (1 + CP)*(1 - CP)                                           RAY09200
      RA  = SQRT(RAA)                                                   RAY09210
      DO  1  J=1,15                                                     RAY09220
    1   Y(J) = 0                                                        RAY09230
C                                                                       RAY09240
C     LIQUID BOTTOM                                                     RAY09250
C                                                                       RAY09260
      IF( SV.LE.0 )  THEN                                               RAY09270
        Y(1) = RA*EPS                                                   RAY09280
        Y(2) =-ROC*EPS                                                  RAY09290
        JX = 2                                                          RAY09300
C                                                                       RAY09310
        IF( IGG.GE.2 )  THEN                                            RAY09320
          Y(3) =-CP**2*EPS/RA                                           RAY09330
          Y(4) =-2*ROC*EPS                                              RAY09340
          JX = 4                                                        RAY09350
C                                                                       RAY09360
          IF( IGG.GT.2 )  JX = 6                                        RAY09370
        ENDIF                                                           RAY09380
C                                                                       RAY09390
C     SOLID BOTTOM                                                      RAY09400
C                                                                       RAY09410
      ELSE                                                              RAY09420
        IF( C.GE.SV )  GO TO  92                                        RAY09430
        CS  = C/SV                                                      RAY09440
        RBB = (1 + CS)*(1 - CS)                                         RAY09450
        RB  = SQRT(RBB)                                                 RAY09460
        RG  = 2*RO*VS(I)**2                                             RAY09470
        Y(3) =-RA*EPS                                                   RAY09480
        Y(4) =-RB*EPS                                                   RAY09490
        Y(2) =-EPS*(CP**2*RBB + CS**2)/(ROC*(RA*RB + 1))                RAY09500
        Y(1) = RG*Y(2) + EPS                                            RAY09510
        Y(5) =-RG*(Y(1) + EPS) + ROC*EPS                                RAY09520
        JX = 5                                                          RAY09530
C                                                                       RAY09540
        IF( IGG.GE.2 )  THEN                                            RAY09550
          Y(8) = EPS*CP**2/RA                                           RAY09560
          Y(9) = EPS*CS**2/RB                                           RAY09570
          Y(7) =-(RB*Y(8) + RA*Y(9))/ROC - 2*Y(2)                       RAY09580
          Y(6) = RG*Y(7)                                                RAY09590
          Y(10)=-RG*Y(6) + EPS*ROC*2                                    RAY09600
          JX = 10                                                       RAY09610
C                                                                       RAY09620
          IF( IGG.GT.2 )  JX = 15                                       RAY09630
        ENDIF                                                           RAY09640
      ENDIF                                                             RAY09650
C                                                                       RAY09660
C     INTEGRATE UPWARD                                                  RAY09670
C                                                                       RAY09680
      IF( L.LE.1 )  GO TO  8                                            RAY09690
      DO  7  II=2,L                                                     RAY09700
        I  = I - 1                                                      RAY09710
        RO = RHO(I)                                                     RAY09720
        ROC= RO*CC                                                      RAY09730
        PV = VP(I)                                                      RAY09740
        SV = VS(I)                                                      RAY09750
        IF( PV.LE.0 .OR. RO.LE.0 .OR.                                   RAY09760
     *      H(I).LE.0 )  GO TO  92                                      RAY09770
        DO  2  J=1,15                                                   RAY09780
    2     Z(J) = Y(J)                                                   RAY09790
        IF( (SV.LE.0 .AND. VS(I+1).LE.0) .OR.                           RAY09800
     *      (SV.GT.0 .AND. VS(I+1).GT.0) )  GO TO  3                    RAY09810
C                                                                       RAY09820
C       SOLID  -----> LIQUID BOUNDARY                                   RAY09830
C                                                                       RAY09840
        IF( SV.LE.0 )  THEN                                             RAY09850
          Y0(3) =-Y(1)/Y(3)                                             RAY09860
          Z(1) = Y(3)                                                   RAY09870
          Z(2) = Y(5)                                                   RAY09880
          JX = 2                                                        RAY09890
C                                                                       RAY09900
          IF( IGG.GE.2 )  THEN                                          RAY09910
            Z(3) = Y(8)                                                 RAY09920
            Z(4) = Y(10)                                                RAY09930
            JX = 4                                                      RAY09940
C                                                                       RAY09950
            IF( IGG.GT.2 )  THEN                                        RAY09960
              Z(5) = Y(13)                                              RAY09970
              Z(6) = Y(15)                                              RAY09980
              JX = 6                                                    RAY09990
            ENDIF                                                       RAY10000
          ENDIF                                                         RAY10010
C                                                                       RAY10020
C       LIQUID -----> SOLID BOUNDARY                                    RAY10030
C                                                                       RAY10040
        ELSE                                                            RAY10050
          Z( 2) = Y(1)                                                  RAY10060
          Z( 4) = Y(2)                                                  RAY10070
          Z( 1) = 0                                                     RAY10080
          Z( 3) = 0                                                     RAY10090
          Z( 5) = 0                                                     RAY10100
          JX = 5                                                        RAY10110
C                                                                       RAY10120
          IF( IGG.GE.2 )  THEN                                          RAY10130
            Z( 7) = Y(3)                                                RAY10140
            Z( 9) = Y(4)                                                RAY10150
            Z( 6) = 0                                                   RAY10160
            Z( 8) = 0                                                   RAY10170
            Z(10) = 0                                                   RAY10180
            JX = 10                                                     RAY10190
C                                                                       RAY10200
            IF( IGG.GT.2 )  THEN                                        RAY10210
              Z(12) = Y(5)                                              RAY10220
              Z(14) = Y(6)                                              RAY10230
              Z(11) = 0                                                 RAY10240
              Z(13) = 0                                                 RAY10250
              Z(15) = 0                                                 RAY10260
              JX = 15                                                   RAY10270
            ENDIF                                                       RAY10280
          ENDIF                                                         RAY10290
        ENDIF                                                           RAY10300
C                                                                       RAY10310
    3   R2  = 1/ROC                                                     RAY10320
        CP  = C/PV                                                      RAY10330
        RAA = (1 + CP)*(1 - CP)                                         RAY10340
        CS  = 0                                                         RAY10350
        IF( SV.GT.0 )  CS = C/SV                                        RAY10360
        RBB = (1 + CS)*(1 - CS)                                         RAY10370
        HK  = H(I)*WN                                                   RAY10380
        HKK = HK**2                                                     RAY10390
        XX  = RAA*HKK                                                   RAY10400
        ONE = 1                                                         RAY10410
C                                                                       RAY10420
C       SINH(X)/X                                                       RAY10430
C                                                                       RAY10440
        DO  4  K=1,2                                                    RAY10450
          CHA = CHB                                                     RAY10460
          SHA = SHB                                                     RAY10470
          DHA = DHB                                                     RAY10480
          AA  = ABS(XX)                                                 RAY10490
          IF( AA.LE.1 )  THEN                                           RAY10500
            SHB = SH0(XX)                                               RAY10510
            CHB = 1 + XX*SH0(XX/4)**2/2                                 RAY10520
            IF( IGG.GE.2 )  DHB = SH1(XX)*HKK                           RAY10530
          ELSE                                                          RAY10540
            AA  = SQRT(AA)                                              RAY10550
            IF( XX.LE.0 )  THEN                                         RAY10560
              CHB = COS(AA)                                             RAY10570
              SHB = SIN(AA)/AA                                          RAY10580
            ELSE                                                        RAY10590
              IF( AA.GT.100 )  ONE = 0                                  RAY10600
              IF( AA.LE.100 )  ONE = ONE/COSH(AA)                       RAY10610
              CHB = 1                                                   RAY10620
              SHB = TANH(AA)/AA                                         RAY10630
            ENDIF                                                       RAY10640
            IF( IGG.GE.2 )  DHB = (HKK/XX)*(CHB - SHB)                  RAY10650
          ENDIF                                                         RAY10660
          XX  = HKK*RBB                                                 RAY10670
          SHB = HK*SHB                                                  RAY10680
    4   CONTINUE                                                        RAY10690
C                                                                       RAY10700
C     LAYER MATRICES                                                    RAY10710
C                                                                       RAY10720
C       LIQUID LAYER                                                    RAY10730
C                                                                       RAY10740
        IF( SV.LE.0 )  THEN                                             RAY10750
          B11 = CHA                                                     RAY10760
          B12 =-RAA*SHA*R2                                              RAY10770
          B21 =-ROC*SHA                                                 RAY10780
C                                                                       RAY10790
          Y(1) = B11*Z(1) + B12*Z(2)                                    RAY10800
          Y(2) = B21*Z(1) + B11*Z(2)                                    RAY10810
C                                                                       RAY10820
          IF( IGG.GE.2 )  THEN                                          RAY10830
            C11 =-HK*SHA                                                RAY10840
            C12 = (HK*CHA + (1 + RAA)*SHA)*R2                           RAY10850
            C21 = (HK*DHA - SHA)*ROC                                    RAY10860
C                                                                       RAY10870
            Y(3) = B11*Z(3) + B12*Z(4) + C11*Z(1)                       RAY10880
     *           + C12*Z(2)                                             RAY10890
            Y(4) = B21*Z(3) + B11*Z(4) + C21*Z(1)                       RAY10900
     *           + C11*Z(2)                                             RAY10910
C                                                                       RAY10920
            IF( IGG.GT.2 )  THEN                                        RAY10930
              W11 = HK*RAA*SHA                                          RAY10940
              W12 =-HK*RAA*CHA*R2                                       RAY10950
              W21 =-HK*ROC*CHA                                          RAY10960
C                                                                       RAY10970
              Y(5) = B11*Z(5) + B12*Z(6) + W11*Z(1)                     RAY10980
     *             + W12*Z(2)                                           RAY10990
              Y(6) = B21*Z(5) + B11*Z(6) + W21*Z(1)                     RAY11000
     *             + W11*Z(2)                                           RAY11010
            ENDIF                                                       RAY11020
          ENDIF                                                         RAY11030
C                                                                       RAY11040
C       SOLID LAYER                                                     RAY11050
C                                                                       RAY11060
        ELSE                                                            RAY11070
          G1  = 2/CS**2                                                 RAY11080
          RG  = G1*ROC                                                  RAY11090
          R4  = RG - ROC                                                RAY11100
          E1  = CHA*CHB                                                 RAY11110
          E2  = E1 - ONE                                                RAY11120
          E3  = SHA*SHB                                                 RAY11130
          E5  = SHA*CHB                                                 RAY11140
          E6  = SHB*CHA                                                 RAY11150
          F1  = E2 - E3                                                 RAY11160
          F2  = R2*F1                                                   RAY11170
          F3  = G1*F1 + E3                                              RAY11180
          B33 = E1                                                      RAY11190
          B34 = RAA*E3                                                  RAY11200
          B43 = RBB*E3                                                  RAY11210
          B25 =-R2*(F2 + R2*(E2 - RAA*B43))                             RAY11220
          B15 = RG*B25 + F2                                             RAY11230
          B16 =-RG*B15 - F3                                             RAY11240
          B22 = B16 + E1                                                RAY11250
          B12 = RG*B16 - R4*F3                                          RAY11260
          B52 =-RG*B12 + R4*(RG*F3 + R4*E3)                             RAY11270
          B23 = R2*(E5 - RBB*E6)                                        RAY11280
          B13 = RG*B23 - E5                                             RAY11290
          B42 =-RG*B13 + R4*E5                                          RAY11300
          B24 = R2*(E6 - RAA*E5)                                        RAY11310
          B14 = RG*B24 - E6                                             RAY11320
          B32 =-RG*B14 + R4*E6                                          RAY11330
          B11 = ONE - B16-B16                                           RAY11340
          B21 = B15+B15                                                 RAY11350
          B31 = B14+B14                                                 RAY11360
          B41 = B13+B13                                                 RAY11370
          B51 = B12+B12                                                 RAY11380
C                                                                       RAY11390
          Y(1) = B11*Z(1) + B12*Z(2) + B13*Z(3)                         RAY11400
     *         + B14*Z(4) + B15*Z(5)                                    RAY11410
          Y(2) = B21*Z(1) + B22*Z(2) + B23*Z(3)                         RAY11420
     *         + B24*Z(4) + B25*Z(5)                                    RAY11430
          Y(3) = B31*Z(1) + B32*Z(2) + B33*Z(3)                         RAY11440
     *         + B34*Z(4) + B24*Z(5)                                    RAY11450
          Y(4) = B41*Z(1) + B42*Z(2) + B43*Z(3)                         RAY11460
     *         + B33*Z(4) + B23*Z(5)                                    RAY11470
          Y(5) = B51*Z(1) + B52*Z(2) + B42*Z(3)                         RAY11480
     *         + B32*Z(4) + B22*Z(5)                                    RAY11490
C                                                                       RAY11500
          IF( IGG.GE.2 )  THEN                                          RAY11510
            RAAC =-2*CP*CP                                              RAY11520
            RBBC =-2*CS*CS                                              RAY11530
            R1C = ROC+ROC                                               RAY11540
            E1C =-HK*(E5 + E6)                                          RAY11550
            E3C =-E3-E3 - HK*(DHA*SHB + DHB*SHA)                        RAY11560
            E5C =-E5 - HK*(DHA*CHB + E3)                                RAY11570
            E6C =-E6 - HK*(DHB*CHA + E3)                                RAY11580
            F1C = E1C - E3C                                             RAY11590
            F2C = R2*(F1C - F1-F1)                                      RAY11600
            F3C = G1*(F1C - F1-F1) + E3C                                RAY11610
            C33 = E1C                                                   RAY11620
            C34 = RAA*E3C + RAAC*E3                                     RAY11630
            C43 = RBB*E3C + RBBC*E3                                     RAY11640
            C25 =-R2*(F2C + R2*(E1C - RAA*C43 - RAAC*B43))              RAY11650
     *          - 2*(B25+B25 + R2*F2)                                   RAY11660
            C15 = RG*C25 + F2C                                          RAY11670
            C16 =-RG*C15 - F3C                                          RAY11680
            C22 = C16 + E1C                                             RAY11690
            C12 = RG*C16 + R1C*F3 - R4*F3C                              RAY11700
            C52 =-RG*C12 + R4*(RG*F3C + R4*E3C)                         RAY11710
     *          - R1C*(RG*F3 + 2*R4*E3)                                 RAY11720
            C23 = R2*(E5C - RBB*E6C - RBBC*E6) - B23-B23                RAY11730
            C13 = RG*C23 - E5C                                          RAY11740
            C42 =-RG*C13 + R4*E5C - R1C*E5                              RAY11750
            C24 = R2*(E6C - RAA*E5C - RAAC*E5) - B24-B24                RAY11760
            C14 = RG*C24 - E6C                                          RAY11770
            C32 =-RG*C14 + R4*E6C - R1C*E6                              RAY11780
            C11 =-C16-C16                                               RAY11790
            C21 = C15+C15                                               RAY11800
            C31 = C14+C14                                               RAY11810
            C41 = C13+C13                                               RAY11820
            C51 = C12+C12                                               RAY11830
C                                                                       RAY11840
            Y( 6) = B11*Z(6) + B12*Z(7) + B13*Z(8)                      RAY11850
     *            + B14*Z(9) + B15*Z(10)                                RAY11860
     *            + C11*Z(1) + C12*Z(2) + C13*Z(3)                      RAY11870
     *            + C14*Z(4) + C15*Z(5)                                 RAY11880
            Y( 7) = B21*Z(6) + B22*Z(7) + B23*Z(8)                      RAY11890
     *            + B24*Z(9) + B25*Z(10)                                RAY11900
     *            + C21*Z(1) + C22*Z(2) + C23*Z(3)                      RAY11910
     *            + C24*Z(4) + C25*Z(5)                                 RAY11920
            Y( 8) = B31*Z(6) + B32*Z(7) + B33*Z(8)                      RAY11930
     *            + B34*Z(9) + B24*Z(10)                                RAY11940
     *            + C31*Z(1) + C32*Z(2) + C33*Z(3)                      RAY11950
     *            + C34*Z(4) + C24*Z(5)                                 RAY11960
            Y( 9) = B41*Z(6) + B42*Z(7) + B43*Z(8)                      RAY11970
     *            + B33*Z(9) + B23*Z(10)                                RAY11980
     *            + C41*Z(1) + C42*Z(2) + C43*Z(3)                      RAY11990
     *            + C33*Z(4) + C23*Z(5)                                 RAY12000
            Y(10) = B51*Z(6) + B52*Z(7) + B42*Z(8)                      RAY12010
     *            + B32*Z(9) + B22*Z(10)                                RAY12020
     *            + C51*Z(1) + C52*Z(2) + C42*Z(3)                      RAY12030
     *            + C32*Z(4) + C22*Z(5)                                 RAY12040
C                                                                       RAY12050
            IF( IGG.GT.2 )  THEN                                        RAY12060
              E1W = HK*(RAA*E5 + RBB*E6)                                RAY12070
              E3W = HK*(E5 + E6)                                        RAY12080
              E5W = HK*(E1 + B43)                                       RAY12090
              E6W = HK*(E1 + B34)                                       RAY12100
              F1W = E1W - E3W                                           RAY12110
              F2W = R2*F1W                                              RAY12120
              F3W = G1*F1W + E3W                                        RAY12130
              W33 = E1W                                                 RAY12140
              W34 = RAA*E3W                                             RAY12150
              W43 = RBB*E3W                                             RAY12160
              W25 =-R2*(F2W + R2*(E1W - RAA*W43))                       RAY12170
              W15 = RG*W25 + F2W                                        RAY12180
              W16 =-RG*W15 - F3W                                        RAY12190
              W22 = W16 + E1W                                           RAY12200
              W12 = RG*W16 - R4*F3W                                     RAY12210
              W52 =-RG*W12 + R4*(RG*F3W + R4*E3W)                       RAY12220
              W23 = R2*(E5W - RBB*E6W)                                  RAY12230
              W13 = RG*W23 - E5W                                        RAY12240
              W42 =-RG*W13 + R4*E5W                                     RAY12250
              W24 = R2*(E6W - RAA*E5W)                                  RAY12260
              W14 = RG*W24 - E6W                                        RAY12270
              W32 =-RG*W14 + R4*E6W                                     RAY12280
              W11 =-W16-W16                                             RAY12290
              W21 = W15+W15                                             RAY12300
              W31 = W14+W14                                             RAY12310
              W41 = W13+W13                                             RAY12320
              W51 = W12+W12                                             RAY12330
C                                                                       RAY12340
              Y(11) = B11*Z(11) + B12*Z(12) + B13*Z(13)                 RAY12350
     *              + B14*Z(14) + B15*Z(15)                             RAY12360
     1              + W11*Z( 1) + W12*Z( 2) + W13*Z( 3)                 RAY12370
     2              + W14*Z( 4) + W15*Z( 5)                             RAY12380
              Y(12) = B21*Z(11) + B22*Z(12) + B23*Z(13)                 RAY12390
     *              + B24*Z(14) + B25*Z(15)                             RAY12400
     1              + W21*Z( 1) + W22*Z( 2) + W23*Z( 3)                 RAY12410
     2              + W24*Z( 4) + W25*Z( 5)                             RAY12420
              Y(13) = B31*Z(11) + B32*Z(12) + B33*Z(13)                 RAY12430
     *              + B34*Z(14) + B24*Z(15)                             RAY12440
     1              + W31*Z( 1) + W32*Z( 2) + W33*Z( 3)                 RAY12450
     2              + W34*Z( 4) + W24*Z( 5)                             RAY12460
              Y(14) = B41*Z(11) + B42*Z(12) + B43*Z(13)                 RAY12470
     *              + B33*Z(14) + B23*Z(15)                             RAY12480
     1              + W41*Z( 1) + W42*Z( 2) + W43*Z( 3)                 RAY12490
     2              + W33*Z( 4) + W23*Z( 5)                             RAY12500
              Y(15) = B51*Z(11) + B52*Z(12) + B42*Z(13)                 RAY12510
     *              + B32*Z(14) + B22*Z(15)                             RAY12520
     1              + W51*Z( 1) + W52*Z( 2) + W42*Z( 3)                 RAY12530
     2              + W32*Z( 4) + W22*Z( 5)                             RAY12540
            ENDIF                                                       RAY12550
          ENDIF                                                         RAY12560
        ENDIF                                                           RAY12570
C                                                                       RAY12580
C     NORMALIZATION                                                     RAY12590
C                                                                       RAY12600
        Z1 = 0                                                          RAY12610
        DO  5  J=1,JX                                                   RAY12620
    5     Z1 = MAX(Z1,ABS(Y(J)))                                        RAY12630
        IF( Z1.GT.BIG )  THEN                                           RAY12640
          DO  6  J=1,JX                                                 RAY12650
    6       Y(J) = EPS*Y(J)                                             RAY12660
        ENDIF                                                           RAY12670
C                                                                       RAY12680
    7 CONTINUE                                                          RAY12690
C                                                                       RAY12700
C     NORMAL EXIT                                                       RAY12710
C                                                                       RAY12720
C     LIQUID SURFACE                                                    RAY12730
C                                                                       RAY12740
    8 IF( SV.LE.0 )  THEN                                               RAY12750
        Y0(1) = Y(1)                                                    RAY12760
        IF( ABS(Y(2))*EPS.LE.ABS(Y(1)) )  THEN                          RAY12770
          Y0(2) = Y(2)/ABS(Y(1))                                        RAY12780
        ELSE                                                            RAY12790
          Y0(2) = SIGN(BIG,Y(2))                                        RAY12800
        ENDIF                                                           RAY12810
C                                                                       RAY12820
        IF( IGG.GE.2 )  THEN                                            RAY12830
          EK =-WN*Y(4)/Y(1)                                             RAY12840
          IF( IGG.GT.2 )  U = C*Y(4)/(Y(4) + Y(6))                      RAY12850
        ENDIF                                                           RAY12860
C                                                                       RAY12870
C     SOLID SURFACE                                                     RAY12880
C                                                                       RAY12890
      ELSE                                                              RAY12900
        Y0(1) = Y(3)                                                    RAY12910
        IF( ABS(Y(5))*EPS.LE.ABS(Y(3)) )  THEN                          RAY12920
          Y0(2) = Y(5)/ABS(Y(3))                                        RAY12930
          Y0(3) =-Y(1)/Y(3)                                             RAY12940
        ELSE                                                            RAY12950
          Y0(2) = SIGN(BIG,Y(5))                                        RAY12960
        ENDIF                                                           RAY12970
C                                                                       RAY12980
        IF( IGG.GE.2 )  THEN                                            RAY12990
          EK =-WN*Y(10)/Y(3)                                            RAY13000
          IF( IGG.GT.2 )  U = C*Y(10)/(Y(10) + Y(15))                   RAY13010
        ENDIF                                                           RAY13020
      ENDIF                                                             RAY13030
      RETURN                                                            RAY13040
C                                                                       RAY13050
C     INPUT ERROR                                                       RAY13060
C                                                                       RAY13070
   90 WRITE(6,91)  L,W,C                                                RAY13080
   91   FORMAT(20X,5('?'),3X,'(RAYMRX)   INPUT ERROR',3X,               RAY13090
     *         'L  =',I5,3X,'W  =',1PE13.6,3X,                          RAY13100
     1         'C  =',E13.6,3X,5('?'))                                  RAY13110
      IER =-1                                                           RAY13120
C                                                                       RAY13130
C     DUMMY STATEMENTS                                                  RAY13140
C                                                                       RAY13150
      AP(1) = AP(1)                                                     RAY13160
      AE(1) = AE(1)                                                     RAY13170
      IA = IA                                                           RAY13180
      RETURN                                                            RAY13190
C                                                                       RAY13200		
C   92 WRITE(6,93)  I                                                    RAY13210
C   93   FORMAT(20X,5('?'),3X,'(RAYMRX)   INPUT ERROR ',                 RAY13220
C     *         'AT LAYER',I5,3X,5('?'))                                 RAY13230
C	WRITE(*,*) 'SAMERE',I,RO,C,VP(I),SV
   92  IER =-MAX(1,I)                                                   RAY13240
      RETURN                                                            RAY13250
      END                                                               RAY13260

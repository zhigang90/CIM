      BLOCK DATA
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      SAVE
      COMMON /NATORB/ NATORB(54)
***********************************************************************
*
*     COMMON BLOCKS FOR AM1
*
***********************************************************************
     $       /ALPHA / ALP(54)
     $       /CORES / CORE(54)
     $       /MULTIP/ DD(54),QQ(54),AM(54),AD(54),AQ(54)
     $       /EXPONT/ ZS(54),ZP(54),ZD(54)
     $       /ONELEC/ USS(54),UPP(54),UDD(54)
     $       /BETAS / BETAS(54),BETAP(54),BETAD(54)
     $       /TWOELE/ GSS(54),GSP(54),GPP(54),GP2(54),HSP(54),
     $                GSD(54),GPD(54),GDD(54)
     $       /ATOMIC/ EISOL(54),EHEAT(54)
     $       /AM1REF/ AM1REF(54)
     $       /VSIPS / VS(54),VP(54),VD(54)
     $       /IDEAS / GUESS1(54,10),GUESS2(54,10),GUESS3(54,10)
     $       /GAUSS / FN1(54),FN2(54)
***********************************************************************
*
*     COMMON BLOCKS FOR MNDO
*
***********************************************************************
      COMMON /MNDO/  USSM(54), UPPM(54), UDDM(54), ZSM(54), ZPM(54),
     $ZDM(54), BETASM(54), BETAPM(54), BETADM(54), ALPM(54),
     $EISOLM(54), DDM(54), QQM(54), AMM(54), ADM(54), AQM(54)
     $,GSSM(54),GSPM(54),GPPM(54),GP2M(54),HSPM(54), POLVOM(54)
C**********************************************************************
C
C     COMMON BLOCKS FOR PM3
C
C**********************************************************************
      COMMON /PM3 /  USSPM3(54),UPPPM3(54),UDDPM3(54),ZSPM3(54),
     $      ZPPM3(54),ZDPM3(54),BETASP(54),BETAPP(54),BETADP(54),
     $      ALPPM3(54),EISOLP(54),DDPM3(54),QQPM3(54),AMPM3(54),
     $      ADPM3(54),AQPM3(54),GSSPM3(54),GSPPM3(54),GPPPM3(54),
     $      GP2PM3(54),HSPPM3(54),POLVOP(54)
      COMMON /IDEAP / GuesP1(54,10),GuesP2(54,10),GuesP3(54,10)
***********************************************************************
*
*  COMMON BLOCKS FOR MINDO/3
*
***********************************************************************
      COMMON /ONELE3 /  USS3(18),UPP3(18)
     $       /TWOEL3 /  F03(54)
     $       /ATOMI3 /  EISOL3(18),EHEAT3(18)
     $       /BETA3  /  BETA3(153)
     $       /ALPHA3 /  ALP3(153)
     $       /EXPON3 /  ZS3(18),ZP3(18)
C
C
C   NATORB IS THE NUMBER OF ATOMIC ORBITALS PER ATOM.
C
      DATA NATORB/2*1,8*4,8*4,2*4,9*9,7*4,2*4,9*9,7*4/
      DATA AM1REF /54*1.D0/
      DATA      POLVOM(1) /0.2287D0/
      DATA      POLVOM(6) /0.2647D0/
      DATA      POLVOM(7) /0.3584D0/
      DATA      POLVOM(8) /0.2324D0/
      DATA      POLVOM(9) /0.1982D0/
      DATA      POLVOM(17)/1.3236D0/
      DATA      POLVOM(35)/2.2583D0/
      DATA      POLVOM(53)/4.0930D0/
C
C   CORE IS THE CHARGE ON THE ATOM AS SEEN BY THE ELECTRONS
C
      DATA CORE/1.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0,
     $ 1.D0,2.D0,3.D0,4.D0,5.D0,6.D0,7.D0,8.D0,9.D0,10.D0,11.D0,2.D0,
     $ 3.D0,4.D0,5.D0,6.D0,7.D0,0.D0/
C
C     ENTHALPIES OF FORMATION OF GASEOUS ATOMS ARE TAKEN FROM "ANNUAL
C     REPORTS,1974,71B,P 117".  THERE ARE SOME SIGNIFICANT DIFFERENCES
C     BETWEEN THE VALUES REPORTED THERE AND THE VALUES PREVIOUSLY IN
C     THE BLOCK DATA OF THIS PROGRAM.  ONLY THE THIRD ROW ELEMENTS
C     HAVE BEEN UPDATED.
C
      DATA EHEAT(1)/52.102D0/
      data eheat(2)/ 0.0d0/
      DATA EHEAT(3)/38.41D0/
      DATA EHEAT(4)/76.96D0/
      DATA EHEAT(5)/135.7D0/
      DATA EHEAT(6)/170.89D0/
      DATA EHEAT(7)/113.0D0/
      DATA EHEAT(8)/59.559D0/
      DATA EHEAT(9)/18.89D0/
      data eheat(10)/ 0.0d0/
      data eheat(11)/25.85d0/
      data eheat(12)/35.00d0/
      DATA EHEAT(13)/79.49D0/
      DATA EHEAT(14)/108.39D0/
      DATA EHEAT(15)/75.57D0/
      DATA EHEAT(16)/66.40D0/
      DATA EHEAT(17)/28.99D0/
      data eheat(18)/ 0.0d0/
      data eheat(19)/21.42d0/
      data eheat(20)/42.60d0/
      data eheat(30)/31.17d0/
      data eheat(31)/65.40d0/
      DATA EHEAT(32)/89.5D0/
      data eheat(33)/72.30d0/
      data eheat(34)/54.30d0/
      DATA EHEAT(35)/26.74D0/
      data eheat(36)/ 0.0d0/
      data eheat(37)/19.60d0/
      data eheat(38)/39.10d0/
      data eheat(48)/26.72d0/
      data eheat(49)/58.00d0/
      DATA EHEAT(50)/72.2D0/
      data eheat(51)/63.20d0/
      data eheat(52)/47.00d0/
      DATA EHEAT(53)/25.517D0/
      data eheat(54)/0.0d0/
C *** VS AND VP ARE THE VALENCE STATE IONIZATION POTENTIAL OF S AND P
C     ELECTRONS IN E.V. : USED IN THE EHT RESONANCE INTEGRALS.
C
c -- missing data set to zero - NEEDS TO BE CHECKED  (JB)
      DATA VS(1)/ -13.605D0/
      DATA VS(3)/-5.39D0/
      data vs(4)/ 0.0d0/
      DATA VS(5)/-15.16D0/
      DATA VS(6)/-21.34D0/
      DATA VS(7)/-27.51D0/
      DATA VS(8)/-35.30D0/
      DATA VS(9)/-43.70D0/
      data vs(13)/ 0.0d0/
      DATA VS(14)/-17.82D0/
      DATA VS(15)/-21.10D0/
      DATA VS(16)/-23.84D0/
      DATA VS(17)/-25.26D0/
      data vs(30)/ 0.0d0/
      data vs(32)/ 0.0d0/
      data vs(35)/ 0.0d0/
      data vs(50)/ 0.0d0/
      data vs(53)/ 0.0d0/
      DATA VP(1)/ 0.0D0/
      DATA VP(3)/-3.54D0/
      data vp(4)/ 0.0d0/
      DATA VP(5)/-8.52D0/
      DATA VP(6)/-11.54D0/
      DATA VP(7)/-14.34D0/
      DATA VP(8)/-17.91D0/
      DATA VP(9)/-20.89D0/
      data vp(13)/ 0.0d0/
      DATA VP(14)/-8.51D0/
      DATA VP(15)/-10.29D0/
      DATA VP(16)/-12.41D0/
      DATA VP(17)/-15.09D0/
      data vp(30)/ 0.0d0/
      data vp(32)/ 0.0d0/
      data vp(35)/ 0.0d0/
      data vp(50)/ 0.0d0/
      data vp(53)/ 0.0d0/
C      DATA NPQ/1,1, 2,2,2,2,2,2,2,2, 3,3,3,3,3,3,3,3, 4,4,4,4,4,4,4,4,
C     +4,4,4,4,4,4,4,4,4,4, 5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5/
C
C *** ONE CENTER REPULSION INTEGRALS
C     GSS ::= (SS,SS)
C     GPP ::= (PP,PP)
C     GSP ::= (SS,PP)
C     GP2 ::= (PP,P*P*)
C     HSP ::= (SP,SP)
      DATA GSSM(1)/ 12.848D00 /
cc      data gssm(3)/ 0.0d0/
      DATA GSSM(4)/9.00D00/
      DATA GSSM(5)/10.59D00/
      DATA GSSM(6)/12.23D00 /
      DATA GSSM(7)/13.59D00/
      DATA GSSM(8)/15.42D00/
      DATA GSSM(9)/16.92D00/
      DATA GSSM(13)/8.09D00/
      DATA GSSM(14)/9.82D00/
      DATA GSSM(15)/11.56D00/
      DATA GSSM(16)/12.88D00/
      DATA GSSM(17)/15.03D00/
cc      data gssm(30)/ 0.0d0/
cc      data gssm(32)/ 0.0d0/
      DATA GSSM(35)/15.03643948D0/
cc      data gssm(50)/ 0.0d0/
      DATA GSSM(53)/15.04044855D0/
c
      data gppm(1)/ 0.0d0/
cc      data gppm(3)/ 0.0d0/
      DATA GPPM(4)/6.97D00/
      DATA GPPM(5)/8.86D00/
      DATA GPPM(6)/11.08D00 /
      DATA GPPM(7)/12.98D00/
      DATA GPPM(8)/14.52D00/
      DATA GPPM(9)/16.71D00/
      DATA GPPM(13)/5.98D00/
      DATA GPPM(14)/7.31D00/
      DATA GPPM(15)/8.64D00/
      DATA GPPM(16)/9.90D00/
      DATA GPPM(17)/11.30D00/
cc      data gppm(30)/ 0.0d0/
cc      data gppm(32)/ 0.0d0/
      DATA GPPM(35)/11.27632539D0/
cc      data gppm(50)/ 0.0d0/
      DATA GPPM(53)/11.14778369D0/
c
      data gspm(1)/ 0.0d0/
cc      data gspm(3)/ 0.0d0/
      DATA GSPM(4)/7.43D00/
      DATA GSPM(5)/9.56D00/
      DATA GSPM(6) / 11.47D00 /
      DATA GSPM(7)/12.66D00/
      DATA GSPM(8)/14.48D00/
      DATA GSPM(9)/17.25D00/
      DATA GSPM(13)/6.63D00/
      DATA GSPM(14)/8.36D00/
      DATA GSPM(15)/10.08D00/
      DATA GSPM(16)/11.26D00/
      DATA GSPM(17)/13.16D00/
cc      data gspm(30)/ 0.0d0/
cc      data gspm(32)/ 0.0d0/
      DATA GSPM(35)/13.03468242D0/
cc      data gspm(50)/ 0.0d0/
      DATA GSPM(53)/13.05655798D0/
c
      data gp2m(1)/ 0.0d0/
cc      data gp2m(3)/ 0.0d0/
      DATA GP2M(4)/6.22D00/
      DATA GP2M(5)/7.86D00/
      DATA GP2M(6) / 9.84D00 /
      DATA GP2M(7)/11.59D00/
      DATA GP2M(8)/12.98D00/
      DATA GP2M(9)/14.91D00/
      DATA GP2M(13)/5.40D00/
      DATA GP2M(14)/6.54D00/
      DATA GP2M(15)/7.68D00/
      DATA GP2M(16)/8.83D00/
      DATA GP2M(17)/9.97D00/
cc      data gp2m(30)/ 0.0d0/
cc      data gp2m(32)/ 0.0d0/
      DATA GP2M(35)/9.85442552D0/
cc      data gp2m(50)/ 0.0d0/
      DATA GP2M(53)/9.91409071D0/
c
      data hspm(1)/ 0.0d0/
cc      data hspm(3)/ 0.0d0/
      DATA HSPM(4)/1.28D00/
      DATA HSPM(5)/1.81D00/
      DATA HSPM(6) / 2.43D00 /
      DATA HSPM(7)/3.14D00/
      DATA HSPM(8)/3.94D00/
      DATA HSPM(9)/4.83D00/
      DATA HSPM(13)/0.70D00/
      DATA HSPM(14)/1.32D00/
      DATA HSPM(15)/1.92D00/
      DATA HSPM(16)/2.26D00/
      DATA HSPM(17)/2.42D00/
cc      data hspm(30)/ 0.0d0/
cc      data hspm(32)/ 0.0d0/
      DATA HSPM(35)/2.45586832D0/
cc      data hspm(50)/ 0.0d0/
      DATA HSPM(53)/2.45638202D0/
C
C     THE MONOCENTRIC INTEGRALS HSP AND GSP FOR ALUMINIUM ARE ONLY
C     ESTIMATES. A VALUE OF G1 FOR AL IS NEEDED TO RESOLVE OLEARIS
C     INTEGRALS.
C
C     OPTIMIZED MNDO PARAMETERS FOR H, BE, B, C, N, O, F
C                                                     CL
C     ESTIMATED MNDO PARAMETERS FOR       AL,SI, P, S
C
C     ELEMENTS H ,LI, C, N, O WERE PARAMETERIZED BY WALTER THIEL
C     ELEMENTS B,SI,P,S       WERE      ..          MICHAEL MCKEE
C     ELEMENTS BE,F,AL,CL     WERE      ..          HENRY RZEPA
C
***********************************************************************
*
*    START OF MINDO/3 PARAMETERS
*
***********************************************************************
C *** F03 IS THE ONE CENTER AVERAGED REPULSION INTEGRAL FOR USE IN THE
C        TWO CENTER ELECTRONIC REPULSION INTEGRAL EVALUATION.
      DATA F03              /  12.848D0, 10.0D0, 0.0D0, 0.0D0,
     .  8.958D0, 10.833D0, 12.377D0, 13.985D0, 16.250D0,
     .         10.000D0, 0.000D0, 0.000D0, 0.000D0,
     .    7.57D0 ,  9.00D0 , 10.20D0 , 11.73,10.0D0,35*0.D0,10.D0/
C *** USS AND UPP ARE THE ONE-CENTER CORE ELECTRON ATTRACTION AND KINETI
C     ENERGY INTEGRALS FOR S AND P ELECTRONS RESPECTIVELY IN E.V.
      DATA USS3             / -12.505D0, 0.000D0, 0.000D0, 0.000D0,
     .                       -33.61D0, -51.79D0, -66.06D0, -91.73D0 ,
     .                       -129.86D0,
     .                        0.0000D0 , 0.000 D0 ,0.000D0 , 0.000D0 ,
     .          -39.82D0 , -56.23D0 , -73.39D0 , -98.99D0 ,.0D0/
      DATA UPP3             /   0.0D0, 0.0D0, 0.0D0, 0.0D0,
     .     -25.11D0 , -39.18D0 , -56.40D0 , -78.80D0 , -105.93D0 ,
     .                        0.000D0 , 0.000D0 , 0.000D0 , 0.000D0 ,
     .         -29.15D0 , -42.31D0 , -57.25D0 , -76.43D0 ,.0D0/
C *** EISOL3 AND EHEAT3 ARE THE GS ELECTRONIC ENERGY OF THE NEUTRAL ATOM
C     (IN E.V.) AND THE HEAT OF FORMATION IF THE FREE ATOM (IN KCAL/MOL)
      DATA EISOL3             /-12.505D0 , 0.0D0 , 0.0D0 ,0.0D0 ,
     .        -61.70D0 ,-119.47D0 , -187.51D0 , -307.07D0 , -475.00D0 ,
     .                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     .          -90.98D0 , -150.81D0 , -229.15D0 , -345.93D0 , 0.0D0/
      DATA EHEAT3             / 52.102D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     .     135.7 D0 , 170.89D0 ,  113.0 D0 ,  59.559D0 ,  18.86D0 ,
     .                         0.0D0 , 0.0D0 , 0.0D0 , 0.0D0 ,
     .     106.0D0 ,   79.8D0 ,  65.65D0 ,  28.95D0 , 0.0D0 /
C *** BETA3 AND ALP3 ARE THE BOND PARAMETERS USED IN THE
C     RESONANCE INTEGRAL AND THE CORE CORE REPULSION INTEGRAL RESPECTIVE
C     THAT IS ACCORDING TO THE FOLLOWING CONVENTION
C
C     HERE IS THE
C     BOND TYPE DESIGNATION
C
C
C         H   B   C   N   O   F  SI   P   S  CL
C       -----------------------------------------
C      H  1  11  16  22  29  37  92 106 121 137
C      B     15  20  26  33  41
C      C         21  27  34  42  97 111 126 142
C      N             28  35  43         127 143
C      O                 36  44         128 144
C      F                     45         129
C     SI                        105
C      P                            120     151
C      S                                136 152
C     CL                                    153
      DATA BETA3(1),ALP3(1)   /  0.244770D0 ,  1.489450D0 /
      DATA BETA3(11),ALP3(11)   /  0.185347D0 ,  2.090352D0 /
      DATA BETA3(15),ALP3(15)   /  0.151324D0 ,  2.280544D0 /
      DATA BETA3(16),ALP3(16)   /  0.315011D0 ,  1.475836D0 /
      DATA BETA3(20),ALP3(20)   /  0.250031D0 ,  2.138291D0 /
      DATA BETA3(21),ALP3(21)   /  0.419907D0 ,  1.371208D0 /
      DATA BETA3(22),ALP3(22)   /  0.360776D0 ,  0.589380D0 /
      DATA BETA3(26),ALP3(26)   /  0.310959D0 ,  1.909763D0 /
      DATA BETA3(27),ALP3(27)   /  0.410886D0 ,  1.635259D0 /
      DATA BETA3(28),ALP3(28) /  0.377342D0 ,  2.029618D0 /
      DATA BETA3(29),ALP3(29) /  0.417759D0 ,  0.478901D0 /
      DATA BETA3(33),ALP3(33) /  0.349745D0 ,  2.484827D0 /
      DATA BETA3(34),ALP3(34) /  0.464514D0 ,  1.820975D0 /
      DATA BETA3(35),ALP3(35) /  0.458110D0 ,  1.873859D0 /
      DATA BETA3(36),ALP3(36) /  0.659407D0 ,  1.537190D0 /
      DATA BETA3(37),ALP3(37) /  0.195242D0 ,  3.771362D0 /
      DATA BETA3(41),ALP3(41) /  0.219591D0 ,  2.862183D0 /
      DATA BETA3(42),ALP3(42) /  0.247494D0 ,  2.725913D0 /
      DATA BETA3(43),ALP3(43) /  0.205347D0 ,  2.861667D0 /
      DATA BETA3(44),ALP3(44) /  0.334044D0 ,  2.266949D0 /
      DATA BETA3(45),ALP3(45) /  0.197464D0 ,  3.864997D0 /
      DATA BETA3(92),ALP3(92) /  0.289647D0 ,  0.940789D0 /
      DATA BETA3(97),ALP3(97) /  0.411377D0 ,  1.101382D0 /
      DATA BETA3(105),ALP3(105) /  0.291703D0 ,  0.918432D0 /
      DATA BETA3(106),ALP3(106) /  0.320118D0 ,  0.923170D0 /
      DATA BETA3(111),ALP3(111) /  0.457816D0 ,  1.029693D0 /
      DATA BETA3(120),ALP3(120) /  0.311790D0 ,  1.186652D0 /
      DATA BETA3(121),ALP3(121) /  0.220654D0 ,  1.700698D0 /
      DATA BETA3(126),ALP3(126) /  0.284620D0 ,  1.761370D0 /
      DATA BETA3(127),ALP3(127) /  0.313170D0 ,  1.878176D0/
      DATA BETA3(128),ALP3(128) /  0.422890D0 ,  2.077240D0 /
      DATA BETA3(129),ALP3(129)  /  0.000000D0 ,  0.000000D0 /
      DATA BETA3(136),ALP3(136) /  0.202489D0 ,  1.751617D0 /
      DATA BETA3(137),ALP3(137) /  0.231653D0 ,  2.089404D0 /
      DATA BETA3(142),ALP3(142) /  0.315480D0 ,  1.676222D0 /
      DATA BETA3(143),ALP3(143) /  0.302298D0 ,  1.817064D0 /
      DATA BETA3(144),ALP3(144) /  0.000000D0 ,  0.000000D0 /
      DATA BETA3(151),ALP3(151) /  0.277322D0 ,  1.543720D0 /
      DATA BETA3(152),ALP3(152) /  0.221764D0 ,  1.950318D0 /
      DATA BETA3(153),ALP3(153) /  0.258969D0 ,  1.792125D0 /
C *** HERE COMES THE OPTIMIZED SLATER_S EXPONENTS FOR THE EVALUATION
C     OF THE OVERLAP INTEGRALS AND MOLECULAR DIPOLE MOMENTS.
      DATA ZS3(1),ZP3(1)      /  1.3D0       ,  0.0D0      /
      DATA ZS3(5),ZP3(5)      /  1.211156D0 ,  0.972826D0 /
      DATA ZS3(6),ZP3(6)      /  1.739391D0 ,  1.709645D0 /
      DATA ZS3(7),ZP3(7)      /  2.704546D0 ,  1.870839D0 /
      DATA ZS3(8),ZP3(8)      /  3.640575D0 ,  2.168448D0 /
      DATA ZS3(9),ZP3(9)      /  3.111270D0 ,  1.41986D0 /
      DATA ZS3(14),ZP3(14)    /  1.629173D0 ,  1.381721D0 /
      DATA ZS3(15),ZP3(15)    /  1.926108D0 ,  1.590665D0 /
      DATA ZS3(16),ZP3(16)    /  1.719480D0 ,  1.403205D0 /
      DATA ZS3(17),ZP3(17)    /  3.430887D0 ,  1.627017D0 /
*************************************************************
*                                                           *
*               DATA FOR THE SPARKLES                       *
*                                                           *
*************************************************************
*                               DATA FOR THE " ++ " SPARKLE
cc      DATA NATORB(2)   / 0    /
cc      DATA F03(2)      /10.0D5/
cc      DATA CORE(2)     / 2.0D0/
cc      DATA EHEAT(2)    / 0.0D0/
cc      DATA VS(2)       /10.0D0/
cc      DATA ALP(2)      / 1.5D0/
cc      DATA EISOL(2)    / 0.0D0/
cc      DATA AM(2)       / 0.5D0/
cc      DATA ALPM(2)      / 1.5D0/
cc      DATA EISOLM(2)    / 0.0D0/
cc      DATA AMM(2)       / 0.5D0/
*                               DATA FOR THE " + " SPARKLE
cc      DATA NATORB(10)   / 0    /
cc      DATA F03(10)      /10.0D5/
cc      DATA CORE(10)     / 1.0D0/
cc      DATA EHEAT(10)    / 0.0D0/
cc      DATA VS(10)       /10.0D0/
cc      DATA ALP(10)      / 1.5D0/
cc      DATA EISOL(10)    / 0.0D0/
cc      DATA AM(10)       / 0.5D0/
cc      DATA ALPM(10)      / 1.5D0/
cc      DATA EISOLM(10)    / 0.0D0/
cc      DATA AMM(10)       / 0.5D0/
*                               DATA FOR THE " -- " SPARKLE
cc      DATA NATORB(18)   / 0    /
cc      DATA F03(18)      /10.0D5/
cc      DATA CORE(18)     /-2.0D0/
cc      DATA EHEAT(18)    / 0.0D0/
cc      DATA VS(18)       /10.0D0/
cc      DATA ALP(18)      / 1.5D0/
cc      DATA EISOL(18)    / 0.0D0/
cc      DATA AM(18)       / 0.5D0/
cc      DATA ALPM(18)      / 1.5D0/
cc      DATA EISOLM(18)    / 0.0D0/
cc      DATA AMM(18)       / 0.5D0/
*                               DATA FOR THE " - " SPARKLE
cc      DATA NATORB(36)   / 0    /
cc      DATA F03(36)      /10.0D5/
cc      DATA CORE(36)     /-1.0D0/
cc      DATA EHEAT(36)    / 0.0D0/
cc      DATA VS(36)       /10.0D0/
cc      DATA ALP(36)      / 1.5D0/
cc      DATA EISOL(36)    / 0.0D0/
cc      DATA AM(36)       / 0.5D0/
cc      DATA ALPM(36)      / 1.5D0/
cc      DATA EISOLM(36)    / 0.0D0/
cc      DATA AMM(36)       / 0.5D0/
***********************************************************************
*
*    START OF MNDO PARAMETERS
*
***********************************************************************
C                    DATA FOR ELEMENT  1        HYDROGEN
      DATA USSM   ( 1)/     -11.9062760D0/
      DATA BETASM ( 1)/      -6.9890640D0/
      DATA ZSM    ( 1)/       1.3319670D0/
      DATA ALPM   ( 1)/       2.5441341D0/
      DATA EISOLM ( 1)/     -11.9062760D0/
      DATA AMM    ( 1)/       0.4721793D0/
      DATA ADM    ( 1)/       0.4721793D0/
      DATA AQM    ( 1)/       0.4721793D0/
C                    DATA FOR ELEMENT  3        LITHIUM
      DATA USSM   ( 3)/      -5.1280000D0/
      DATA UPPM   ( 3)/      -2.7212000D0/
      DATA BETASM ( 3)/      -1.3500400D0/
      DATA BETAPM ( 3)/      -1.3500400D0/
      DATA ZSM    ( 3)/       0.7023800D0/
      DATA ZPM    ( 3)/       0.7023800D0/
      DATA ALPM   ( 3)/       1.2501400D0/
      DATA EISOLM ( 3)/      -5.1280000D0/
      DATA GSSM   ( 3)/       7.3D00     /
      DATA GPPM   ( 3)/       5.0D00     /
      DATA GSPM   ( 3)/       5.42D00    /
      DATA GP2M   ( 3)/       4.52D00    /
      DATA HSPM   ( 3)/       0.83D00    /
      DATA DDM    ( 3)/       2.05497832D0/
      DATA QQM    ( 3)/       1.74370693D0/
      DATA AMM    ( 3)/       0.26828372D0/
      DATA ADM    ( 3)/       0.22697935D0/
      DATA AQM    ( 3)/       0.26145812D0/
C                    DATA FOR ELEMENT  4        BERYLLIUM
      DATA USSM   ( 4)/     -16.6023780D0/
      DATA UPPM   ( 4)/     -10.7037710D0/
      DATA BETASM ( 4)/      -4.0170960D0/
      DATA BETAPM ( 4)/      -4.0170960D0/
      DATA ZSM    ( 4)/       1.0042100D0/
      DATA ZPM    ( 4)/       1.0042100D0/
      DATA ALPM   ( 4)/       1.6694340D0/
      DATA EISOLM ( 4)/     -24.2047560D0/
      DATA DDM    ( 4)/       1.4373245D0/
      DATA QQM    ( 4)/       1.2196103D0/
      DATA AMM    ( 4)/       0.3307607D0/
      DATA ADM    ( 4)/       0.3356142D0/
      DATA AQM    ( 4)/       0.3846373D0/
C                    DATA FOR ELEMENT  5        BORON
      DATA USSM   ( 5)/     -34.5471300D0/
      DATA UPPM   ( 5)/     -23.1216900D0/
      DATA BETASM ( 5)/      -8.2520540D0/
      DATA BETAPM ( 5)/      -8.2520540D0/
      DATA ZSM    ( 5)/       1.5068010D0/
      DATA ZPM    ( 5)/       1.5068010D0/
      DATA ALPM   ( 5)/       2.1349930D0/
      DATA EISOLM ( 5)/     -64.3159500D0/
      DATA DDM    ( 5)/       0.9579073D0/
      DATA QQM    ( 5)/       0.8128113D0/
      DATA AMM    ( 5)/       0.3891951D0/
      DATA ADM    ( 5)/       0.4904730D0/
      DATA AQM    ( 5)/       0.5556979D0/
C                    DATA FOR ELEMENT  6        CARBON
      DATA USSM   ( 6)/     -52.2797450D0/
      DATA UPPM   ( 6)/     -39.2055580D0/
      DATA BETASM ( 6)/     -18.9850440D0/
      DATA BETAPM ( 6)/      -7.9341220D0/
      DATA ZSM    ( 6)/       1.7875370D0/
      DATA ZPM    ( 6)/       1.7875370D0/
      DATA ALPM   ( 6)/       2.5463800D0/
      DATA EISOLM ( 6)/    -120.5006060D0/
      DATA DDM    ( 6)/       0.8074662D0/
      DATA QQM    ( 6)/       0.6851578D0/
      DATA AMM    ( 6)/       0.4494671D0/
      DATA ADM    ( 6)/       0.6149474D0/
      DATA AQM    ( 6)/       0.6685897D0/
C                    DATA FOR ELEMENT  7        NITROGEN
      DATA USSM   ( 7)/     -71.9321220D0/
      DATA UPPM   ( 7)/     -57.1723190D0/
      DATA BETASM ( 7)/     -20.4957580D0/
      DATA BETAPM ( 7)/     -20.4957580D0/
      DATA ZSM    ( 7)/       2.2556140D0/
      DATA ZPM    ( 7)/       2.2556140D0/
      DATA ALPM   ( 7)/       2.8613420D0/
      DATA EISOLM ( 7)/    -202.5812010D0/
      DATA DDM    ( 7)/       0.6399037D0/
      DATA QQM    ( 7)/       0.5429763D0/
      DATA AMM    ( 7)/       0.4994487D0/
      DATA ADM    ( 7)/       0.7843643D0/
      DATA AQM    ( 7)/       0.8144720D0/
C                    DATA FOR ELEMENT  8        OXYGEN
      DATA USSM   ( 8)/     -99.6443090D0/
      DATA UPPM   ( 8)/     -77.7974720D0/
      DATA BETASM ( 8)/     -32.6880820D0/
      DATA BETAPM ( 8)/     -32.6880820D0/
      DATA ZSM    ( 8)/       2.6999050D0/
      DATA ZPM    ( 8)/       2.6999050D0/
      DATA ALPM   ( 8)/       3.1606040D0/
      DATA EISOLM ( 8)/    -317.8685060D0/
      DATA DDM    ( 8)/       0.5346024D0/
      DATA QQM    ( 8)/       0.4536252D0/
      DATA AMM    ( 8)/       0.5667034D0/
      DATA ADM    ( 8)/       0.9592562D0/
      DATA AQM    ( 8)/       0.9495934D0/
C                    DATA FOR ELEMENT  9        FLUORINE
      DATA USSM   ( 9)/    -131.0715480D0/
      DATA UPPM   ( 9)/    -105.7821370D0/
      DATA BETASM ( 9)/     -48.2904660D0/
      DATA BETAPM ( 9)/     -36.5085400D0/
      DATA ZSM    ( 9)/       2.8484870D0/
      DATA ZPM    ( 9)/       2.8484870D0/
      DATA ALPM   ( 9)/       3.4196606D0/
      DATA EISOLM ( 9)/    -476.6837810D0/
      DATA DDM    ( 9)/       0.5067166D0/
      DATA QQM    ( 9)/       0.4299633D0/
      DATA AMM    ( 9)/       0.6218302D0/
      DATA ADM    ( 9)/       1.0850301D0/
      DATA AQM    ( 9)/       1.0343643D0/
C                    DATA FOR ELEMENT 13        ALUMINIUM
      DATA USSM   (13)/     -23.8070970D0/
      DATA UPPM   (13)/     -17.5198780D0/
      DATA BETASM (13)/      -2.6702840D0/
      DATA BETAPM (13)/      -2.6702840D0/
      DATA ZSM    (13)/       1.4441610D0/
      DATA ZPM    (13)/       1.4441610D0/
      DATA ZDM    (13)/       1.0000000D0/
      DATA ALPM   (13)/       1.8688394D0/
**************************************************************************
*
*   THE FOLLOWING PARAMETER IS IN ERROR, DO NOT CORRECT IT AS THE ERROR
*   WAS INTRODUCED IN THE ORIGINAL PARAMETRIZATION, AND HAS BEEN ABSORBED
*   BY THE OTHER PARAMETERS.
*
      DATA EISOLM (13)/     -44.4840711D0/
*
*   THE CORRECT VALUE SHOULD BE -43.0840720D0
*
*************************************************************************
      DATA DDM    (13)/       1.3992387D0/
      DATA QQM    (13)/       1.1586797D0/
      DATA AMM    (13)/       0.2973172D0/
      DATA ADM    (13)/       0.2635574D0/
      DATA AQM    (13)/       0.3673560D0/
C                    DATA FOR ELEMENT 14        SILICON
      DATA USSM   (14)/     -37.0375330D0/
      DATA UPPM   (14)/     -27.7696780D0/
      DATA BETASM (14)/      -9.0868040D0/
      DATA BETAPM (14)/      -1.0758270D0/
      DATA ZSM    (14)/       1.3159860D0/
      DATA ZPM    (14)/       1.7099430D0/
      DATA ZDM    (14)/       1.0000000D0/
      DATA ALPM   (14)/       2.2053160D0/
      DATA EISOLM (14)/     -82.8394220D0/
      DATA DDM    (14)/       1.2580349D0/
      DATA QQM    (14)/       0.9785824D0/
      DATA AMM    (14)/       0.3608967D0/
      DATA ADM    (14)/       0.3664244D0/
      DATA AQM    (14)/       0.4506740D0/
C                    DATA FOR ELEMENT 15        PHOSPHORUS
      DATA USSM   (15)/     -56.1433600D0/
      DATA UPPM   (15)/     -42.8510800D0/
      DATA BETASM (15)/      -6.7916000D0/
      DATA BETAPM (15)/      -6.7916000D0/
      DATA ZSM    (15)/       2.1087200D0/
      DATA ZPM    (15)/       1.7858100D0/
      DATA ZDM    (15)/       1.0000000D0/
      DATA ALPM   (15)/       2.4152800D0/
      DATA EISOLM (15)/    -152.9599600D0/
      DATA DDM    (15)/       1.0129699D0/
      DATA QQM    (15)/       0.9370090D0/
      DATA AMM    (15)/       0.4248438D0/
      DATA ADM    (15)/       0.4882420D0/
      DATA AQM    (15)/       0.4979406D0/
C                    DATA FOR ELEMENT 16        SULPHUR
      DATA USSM   (16)/     -72.2422810D0/
      DATA UPPM   (16)/     -56.9732070D0/
      DATA BETASM (16)/     -10.7616700D0/
      DATA BETAPM (16)/     -10.1084330D0/
      DATA ZSM    (16)/       2.3129620D0/
      DATA ZPM    (16)/       2.0091460D0/
      DATA ZDM    (16)/       1.0000000D0/
      DATA ALPM   (16)/       2.4780260D0/
      DATA EISOLM (16)/    -226.0123900D0/
      DATA DDM    (16)/       0.9189935D0/
      DATA QQM    (16)/       0.8328514D0/
      DATA AMM    (16)/       0.4733554D0/
      DATA ADM    (16)/       0.5544502D0/
      DATA AQM    (16)/       0.5585244D0/
C                    DATA FOR ELEMENT 17        CHLORINE
      DATA USSM   (17)/    -100.2271660D0/
      DATA UPPM   (17)/     -77.3786670D0/
      DATA BETASM (17)/     -14.2623200D0/
      DATA BETAPM (17)/     -14.2623200D0/
      DATA ZSM    (17)/       3.7846450D0/
      DATA ZPM    (17)/       2.0362630D0/
      DATA ZDM    (17)/       1.0000000D0/
      DATA ALPM   (17)/       2.5422010D0/
      DATA EISOLM (17)/    -353.1176670D0/
      DATA DDM    (17)/       0.4986870D0/
      DATA QQM    (17)/       0.8217603D0/
      DATA AMM    (17)/       0.5523705D0/
      DATA ADM    (17)/       0.8061220D0/
      DATA AQM    (17)/       0.6053435D0/
C                    DATA FOR ELEMENT 30       ZINC
      DATA USSM   (30)/     -20.8397160D0/
      DATA UPPM   (30)/     -19.6252240D0/
      DATA BETASM (30)/      -1.0000000D0/
      DATA BETAPM (30)/      -2.0000000D0/
      DATA ZSM    (30)/       2.0473590D0/
      DATA ZPM    (30)/       1.4609460D0/
      DATA ZDM    (30)/       1.0000000D0/
      DATA ALPM   (30)/       1.5064570D0/
      DATA EISOLM (30)/     -29.8794320D0/
      DATA GSSM   (30)/      11.8000000D0/
      DATA GSPM   (30)/      11.1820180D0/
      DATA GPPM   (30)/      13.3000000D0/
      DATA GP2M   (30)/      12.9305200D0/
      DATA HSPM   (30)/       0.4846060D0/
      DATA DDM    (30)/       1.3037826D0/
      DATA QQM    (30)/       1.4520183D0/
      DATA AMM    (30)/       0.4336641D0/
      DATA ADM    (30)/       0.2375912D0/
      DATA AQM    (30)/       0.2738858D0/
C                    DATA FOR ELEMENT 32       GERMANIUM
      DATA USSM   (32)/     -33.9493670D0/
      DATA UPPM   (32)/     -27.4251050D0/
      DATA BETASM (32)/      -4.5164790D0/
      DATA BETAPM (32)/      -1.7555170D0/
      DATA ZSM    (32)/       1.2931800D0/
      DATA ZPM    (32)/       2.0205640D0/
      DATA ZDM    (32)/       1.0000000D0/
      DATA ALPM   (32)/       1.9784980D0/
      DATA EISOLM (32)/     -76.2489440D0/
      DATA GSSM   (32)/       9.8000000D0/
      DATA GSPM   (32)/       8.3000000D0/
      DATA GPPM   (32)/       7.3000000D0/
      DATA GP2M   (32)/       6.5000000D0/
      DATA HSPM   (32)/       1.3000000D0/
      DATA DDM    (32)/       1.2556091D0/
      DATA QQM    (32)/       1.0498655D0/
      DATA AMM    (32)/       0.3601617D0/
      DATA ADM    (32)/       0.3643722D0/
      DATA AQM    (32)/       0.4347337D0/
C                    DATA FOR ELEMENT 35        BROMINE
      DATA USSM   (35)/     -99.9864405D0/
      DATA UPPM   (35)/     -75.6713075D0/
      DATA BETASM (35)/      -8.9171070D0/
      DATA BETAPM (35)/      -9.9437400D0/
      DATA ZSM    (35)/       3.8543019D0/
      DATA ZPM    (35)/       2.1992091D0/
      DATA ZDM    (35)/       1.0000000D0/
      DATA ALPM   (35)/       2.4457051D0/
      DATA EISOLM (35)/    -346.6812500D0/
      DATA DDM    (35)/       0.6051074D0/
      DATA QQM    (35)/       0.9645873D0/
      DATA AMM    (35)/       0.5526068D0/
      DATA ADM    (35)/       0.7258330D0/
      DATA AQM    (35)/       0.5574589D0/
C                    DATA FOR ELEMENT 50        TIN
      DATA USSM   (50)/     -40.8518020D0/
      DATA UPPM   (50)/     -28.5602490D0/
      DATA BETASM (50)/      -3.2351470D0/
      DATA BETAPM (50)/      -4.2904160D0/
      DATA ZSM    (50)/       2.0803800D0/
      DATA ZPM    (50)/       1.9371060D0/
      DATA ALPM   (50)/       1.8008140D0/
      DATA EISOLM (50)/     -92.3241020D0/
      DATA GSSM   (50)/       9.8000000D0/
      DATA GSPM   (50)/       8.3000000D0/
      DATA GPPM   (50)/       7.3000000D0/
      DATA GP2M   (50)/       6.5000000D0/
      DATA HSPM   (50)/       1.3000000D0/
      DATA DDM    (50)/       1.5697766D0/
      DATA QQM    (50)/       1.3262292D0/
      DATA AMM    (50)/       0.3601617D0/
      DATA ADM    (50)/       0.3219998D0/
      DATA AQM    (50)/       0.3713827D0/
C                    DATA FOR ELEMENT 53        IODINE
      DATA USSM   (53)/    -100.0030538D0/
      DATA UPPM   (53)/     -74.6114692D0/
      DATA BETASM (53)/      -7.4144510D0/
      DATA BETAPM (53)/      -6.1967810D0/
      DATA ZSM    (53)/       2.2729610D0/
      DATA ZPM    (53)/       2.1694980D0/
      DATA ZDM    (53)/       1.0000000D0/
      DATA ALPM   (53)/       2.2073200D0/
      DATA EISOLM (53)/    -340.5983600D0/
      DATA DDM    (53)/       1.4253233D0/
      DATA QQM    (53)/       1.1841707D0/
      DATA AMM    (53)/       0.5527541D0/
      DATA ADM    (53)/       0.4593451D0/
      DATA AQM    (53)/       0.4585376D0/
***********************************************************************
*
*    START OF AM1 PARAMETERS
*
***********************************************************************
C                    DATA FOR ELEMENT  1
      DATA USS   ( 1)/     -11.3964270D0/
      DATA BETAS ( 1)/      -6.1737870D0/
      DATA ZS    ( 1)/       1.1880780D0/
      DATA ALP   ( 1)/       2.8823240D0/
      DATA EISOL ( 1)/     -11.3964270D0/
      DATA GSS   ( 1)/      12.8480000D0/
      DATA AM    ( 1)/       0.4721793D0/
      DATA AD    ( 1)/       0.4721793D0/
      DATA AQ    ( 1)/       0.4721793D0/
      DATA GUESS1( 1,1)/       0.1227960D0/
      DATA GUESS2( 1,1)/       5.0000000D0/
      DATA GUESS3( 1,1)/       1.2000000D0/
      DATA GUESS1( 1,2)/       0.0050900D0/
      DATA GUESS2( 1,2)/       5.0000000D0/
      DATA GUESS3( 1,2)/       1.8000000D0/
      DATA GUESS1( 1,3)/      -0.0183360D0/
      DATA GUESS2( 1,3)/       2.0000000D0/
      DATA GUESS3( 1,3)/       2.1000000D0/
      Data (Guess1(1,J),J=4,10)/7*0.0d0/
      Data (Guess2(1,J),J=4,10)/7*0.0d0/
      Data (Guess3(1,J),J=4,10)/7*0.0d0/
C                    DATA FOR ELEMENT  4
      DATA USS   ( 4)/     -16.6023780D0/
      DATA UPP   ( 4)/     -10.7037710D0/
      DATA BETAS ( 4)/      -4.0170960D0/
      DATA BETAP ( 4)/      -4.0170960D0/
      DATA ZS    ( 4)/       1.0042100D0/
      DATA ZP    ( 4)/       1.0042100D0/
      DATA ALP   ( 4)/       1.6694340D0/
      DATA EISOL ( 4)/     -24.2047560D0/
      DATA GSS   ( 4)/       9.0000000D0/
      DATA GSP   ( 4)/       7.4300000D0/
      DATA GPP   ( 4)/       6.9700000D0/
      DATA GP2   ( 4)/       6.2200000D0/
      DATA HSP   ( 4)/       1.2800000D0/
      DATA DD    ( 4)/       1.4373245D0/
      DATA QQ    ( 4)/       1.2196103D0/
      DATA AM    ( 4)/       0.3307607D0/
      DATA AD    ( 4)/       0.3356142D0/
      DATA AQ    ( 4)/       0.3846373D0/
      Data (Guess1(4,J),J=1,10)/10*0.0d0/
      Data (Guess2(4,J),J=1,10)/10*0.0d0/
      Data (Guess3(4,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT  5  (MAY 1990)
      DATA USS   ( 5)/    -43.4928700D0/
      DATA UPP   ( 5)/    -22.6315250D0/
      DATA BETAS ( 5)/     -9.5991140D0/
      DATA BETAP ( 5)/     -6.2737570D0/
      DATA ZS    ( 5)/      1.6117090D0/
      DATA ZP    ( 5)/      1.5553850D0/
      DATA ALP   ( 5)/      2.4469090D0/
      DATA EISOL ( 5)/    -63.7172650D0/
      DATA GSS   ( 5)/     10.5900000D0/
      DATA GSP   ( 5)/      9.5600000D0/
      DATA GPP   ( 5)/      8.8600000D0/
      DATA GP2   ( 5)/      7.8600000D0/
      DATA HSP   ( 5)/      1.8100000D0/
      DATA DD    ( 5)/      0.9107622D0/
      DATA QQ    ( 5)/      0.7874223D0/
      DATA AM    ( 5)/      0.3891951D0/
      DATA AD    ( 5)/      0.5045152D0/
      DATA AQ    ( 5)/      0.5678856D0/
      DATA GUESS1( 5,1)/      0.1826130D0/
      DATA GUESS2( 5,1)/      6.0000000D0/
      DATA GUESS3( 5,1)/      0.7275920D0/
      DATA GUESS1( 5,2)/      0.1185870D0/
      DATA GUESS2( 5,2)/      6.0000000D0/
      DATA GUESS3( 5,2)/      1.4666390D0/
      DATA GUESS1( 5,3)/     -0.0732800D0/
      DATA GUESS2( 5,3)/      5.0000000D0/
      DATA GUESS3( 5,3)/      1.5709750D0/
      DATA GUESS1( 5,4)/      0.4122530D0/
      DATA GUESS2( 5,4)/     10.0000000D0/
      DATA GUESS3( 5,4)/      0.8325860D0/
      DATA GUESS1( 5,5)/     -0.1499170D0/
      DATA GUESS2( 5,5)/      6.0000000D0/
      DATA GUESS3( 5,5)/      1.1862200D0/
      DATA GUESS1( 5,6)/      0.2617510D0/
      DATA GUESS2( 5,6)/      8.0000000D0/
      DATA GUESS3( 5,6)/      1.0639950D0/
      DATA GUESS1( 5,7)/      0.0502750D0/
      DATA GUESS2( 5,7)/      5.0000000D0/
      DATA GUESS3( 5,7)/      1.9364920D0/
      DATA GUESS1( 5,8)/      0.3592440D0/
      DATA GUESS2( 5,8)/      9.0000000D0/
      DATA GUESS3( 5,8)/      0.8193510D0/
      DATA GUESS1( 5,9)/      0.0747290D0/
      DATA GUESS2( 5,9)/      9.0000000D0/
      DATA GUESS3( 5,9)/      1.5744140D0/
      data guess1(5,10)/ 0.0d0/
      data guess2(5,10)/ 0.0d0/
      data guess3(5,10)/ 0.0d0/
C                    DATA FOR ELEMENT  6
      DATA USS   ( 6)/     -52.0286580D0/
      DATA UPP   ( 6)/     -39.6142390D0/
      DATA BETAS ( 6)/     -15.7157830D0/
      DATA BETAP ( 6)/      -7.7192830D0/
      DATA ZS    ( 6)/       1.8086650D0/
      DATA ZP    ( 6)/       1.6851160D0/
      DATA ALP   ( 6)/       2.6482740D0/
      DATA EISOL ( 6)/    -120.8157940D0/
      DATA GSS   ( 6)/      12.2300000D0/
      DATA GSP   ( 6)/      11.4700000D0/
      DATA GPP   ( 6)/      11.0800000D0/
      DATA GP2   ( 6)/       9.8400000D0/
      DATA HSP   ( 6)/       2.4300000D0/
      DATA DD    ( 6)/       0.8236736D0/
      DATA QQ    ( 6)/       0.7268015D0/
      DATA AM    ( 6)/       0.4494671D0/
      DATA AD    ( 6)/       0.6082946D0/
      DATA AQ    ( 6)/       0.6423492D0/
      DATA GUESS1( 6,1)/       0.0113550D0/
      DATA GUESS2( 6,1)/       5.0000000D0/
      DATA GUESS3( 6,1)/       1.6000000D0/
      DATA GUESS1( 6,2)/       0.0459240D0/
      DATA GUESS2( 6,2)/       5.0000000D0/
      DATA GUESS3( 6,2)/       1.8500000D0/
      DATA GUESS1( 6,3)/      -0.0200610D0/
      DATA GUESS2( 6,3)/       5.0000000D0/
      DATA GUESS3( 6,3)/       2.0500000D0/
      DATA GUESS1( 6,4)/      -0.0012600D0/
      DATA GUESS2( 6,4)/       5.0000000D0/
      DATA GUESS3( 6,4)/       2.6500000D0/
      Data (Guess1(6,J),J=5,10)/6*0.0d0/
      Data (Guess2(6,J),J=5,10)/6*0.0d0/
      Data (Guess3(6,J),J=5,10)/6*0.0d0/
C                    DATA FOR ELEMENT  7
      DATA USS   ( 7)/     -71.8600000D0/
      DATA UPP   ( 7)/     -57.1675810D0/
      DATA BETAS ( 7)/     -20.2991100D0/
      DATA BETAP ( 7)/     -18.2386660D0/
      DATA ZS    ( 7)/       2.3154100D0/
      DATA ZP    ( 7)/       2.1579400D0/
      DATA ALP   ( 7)/       2.9472860D0/
      DATA EISOL ( 7)/    -202.4077430D0/
      DATA GSS   ( 7)/      13.5900000D0/
      DATA GSP   ( 7)/      12.6600000D0/
      DATA GPP   ( 7)/      12.9800000D0/
      DATA GP2   ( 7)/      11.5900000D0/
      DATA HSP   ( 7)/       3.1400000D0/
      DATA DD    ( 7)/       0.6433247D0/
      DATA QQ    ( 7)/       0.5675528D0/
      DATA AM    ( 7)/       0.4994487D0/
      DATA AD    ( 7)/       0.7820840D0/
      DATA AQ    ( 7)/       0.7883498D0/
      DATA GUESS1( 7,1)/       0.0252510D0/
      DATA GUESS2( 7,1)/       5.0000000D0/
      DATA GUESS3( 7,1)/       1.5000000D0/
      DATA GUESS1( 7,2)/       0.0289530D0/
      DATA GUESS2( 7,2)/       5.0000000D0/
      DATA GUESS3( 7,2)/       2.1000000D0/
      DATA GUESS1( 7,3)/      -0.0058060D0/
      DATA GUESS2( 7,3)/       2.0000000D0/
      DATA GUESS3( 7,3)/       2.4000000D0/
      Data (Guess1(7,J),J=4,10)/7*0.0d0/
      Data (Guess2(7,J),J=4,10)/7*0.0d0/
      Data (Guess3(7,J),J=4,10)/7*0.0d0/
C                    DATA FOR ELEMENT  8
      DATA USS   ( 8)/     -97.8300000D0/
      DATA UPP   ( 8)/     -78.2623800D0/
      DATA BETAS ( 8)/     -29.2727730D0/
      DATA BETAP ( 8)/     -29.2727730D0/
      DATA ZS    ( 8)/       3.1080320D0/
      DATA ZP    ( 8)/       2.5240390D0/
      DATA ALP   ( 8)/       4.4553710D0/
      DATA EISOL ( 8)/    -316.0995200D0/
      DATA GSS   ( 8)/      15.4200000D0/
      DATA GSP   ( 8)/      14.4800000D0/
      DATA GPP   ( 8)/      14.5200000D0/
      DATA GP2   ( 8)/      12.9800000D0/
      DATA HSP   ( 8)/       3.9400000D0/
      DATA DD    ( 8)/       0.4988896D0/
      DATA QQ    ( 8)/       0.4852322D0/
      DATA AM    ( 8)/       0.5667034D0/
      DATA AD    ( 8)/       0.9961066D0/
      DATA AQ    ( 8)/       0.9065223D0/
      DATA GUESS1( 8,1)/       0.2809620D0/
      DATA GUESS2( 8,1)/       5.0000000D0/
      DATA GUESS3( 8,1)/       0.8479180D0/
      DATA GUESS1( 8,2)/       0.0814300D0/
      DATA GUESS2( 8,2)/       7.0000000D0/
      DATA GUESS3( 8,2)/       1.4450710D0/
      Data (Guess1(8,J),J=3,10)/8*0.0d0/
      Data (Guess2(8,J),J=3,10)/8*0.0d0/
      Data (Guess3(8,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  9
      DATA USS   ( 9)/    -136.1055790D0/
      DATA UPP   ( 9)/    -104.8898850D0/
      DATA BETAS ( 9)/     -69.5902770D0/
      DATA BETAP ( 9)/     -27.9223600D0/
      DATA ZS    ( 9)/       3.7700820D0/
      DATA ZP    ( 9)/       2.4946700D0/
      DATA ALP   ( 9)/       5.5178000D0/
      DATA EISOL ( 9)/    -482.2905830D0/
      DATA GSS   ( 9)/      16.9200000D0/
      DATA GSP   ( 9)/      17.2500000D0/
      DATA GPP   ( 9)/      16.7100000D0/
      DATA GP2   ( 9)/      14.9100000D0/
      DATA HSP   ( 9)/       4.8300000D0/
      DATA DD    ( 9)/       0.4145203D0/
      DATA QQ    ( 9)/       0.4909446D0/
      DATA AM    ( 9)/       0.6218302D0/
      DATA AD    ( 9)/       1.2088792D0/
      DATA AQ    ( 9)/       0.9449355D0/
      DATA GUESS1( 9,1)/       0.2420790D0/
      DATA GUESS2( 9,1)/       4.8000000D0/
      DATA GUESS3( 9,1)/       0.9300000D0/
      DATA GUESS1( 9,2)/       0.0036070D0/
      DATA GUESS2( 9,2)/       4.6000000D0/
      DATA GUESS3( 9,2)/       1.6600000D0/
      Data (Guess1(9,J),J=3,10)/8*0.0d0/
      Data (Guess2(9,J),J=3,10)/8*0.0d0/
      Data (Guess3(9,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 12 (taken from J.Phys.Chem.B 102 (1998) 8080)
      DATA USS   (12)/    -14.96959313D0/
      DATA UPP   (12)/    -11.56229248D0/
      DATA BETAS (12)/     -1.25974355D0/
      DATA BETAP (12)/     -0.77836604D0/
      DATA ZS    (12)/      1.22339270D0/
      DATA ZP    (12)/      1.02030798D0/
      DATA ZD    (12)/      1.0000000D0/
      DATA ALP   (12)/      1.67049799D0/
      DATA EISOL (12)/    -00.00000000D0/
      DATA GSS   (12)/      7.50132277D0/
      DATA GSP   (12)/      6.34591536D0/
      DATA GPP   (12)/      4.77534467D0/
      DATA GP2   (12)/      4.34017279D0/
      DATA HSP   (12)/      0.48930466D0/
      DATA DD    (12)/      0.00000000D0/
      DATA QQ    (12)/      0.0000000D0/
      DATA AM    (12)/      0.0000000D0/
      DATA AD    (12)/      0.0000000D0/
      DATA AQ    (12)/      0.0000000D0/
      DATA GUESS1(12,1)/       2.55017735D0/
      DATA GUESS2(12,1)/       4.29397225D0/
      DATA GUESS3(12,1)/       0.79989601D0/
      DATA GUESS1(12,2)/      -0.00565806D0/
      DATA GUESS2(12,2)/       2.96053910D0/
      DATA GUESS3(12,2)/       1.47499983D0/
      DATA GUESS1(12,3)/      -0.00610286D0/
      DATA GUESS2(12,3)/       2.61416919D0/
      DATA GUESS3(12,3)/       2.42604040D0/
      Data (Guess1(12,J),J=4,10)/7*0.0d0/
      Data (Guess2(12,J),J=4,10)/7*0.0d0/
      Data (Guess3(12,J),J=4,10)/7*0.0d0/
C                    DATA FOR ELEMENT 13 (MAY 1990)
      DATA USS   (13)/    -24.3535850D0/
      DATA UPP   (13)/    -18.3636450D0/
      DATA BETAS (13)/     -3.8668220D0/
      DATA BETAP (13)/     -2.3171460D0/
      DATA ZS    (13)/      1.5165930D0/
      DATA ZP    (13)/      1.3063470D0/
      DATA ZD    (13)/      1.0000000D0/
      DATA ALP   (13)/      1.9765860D0/
      DATA EISOL (13)/    -46.4208150D0/
      DATA GSS   (13)/      8.0900000D0/
      DATA GSP   (13)/      6.6300000D0/
      DATA GPP   (13)/      5.9800000D0/
      DATA GP2   (13)/      5.4000000D0/
      DATA HSP   (13)/      0.7000000D0/
      DATA DD    (13)/      1.4040443D0/
      DATA QQ    (13)/      1.2809154D0/
      DATA AM    (13)/      0.2973172D0/
      DATA AD    (13)/      0.2630229D0/
      DATA AQ    (13)/      0.3427832D0/
      DATA GUESS1(13,1)/      0.0900000D0/
      DATA GUESS2(13,1)/     12.3924430D0/
      DATA GUESS3(13,1)/      2.0503940D0/
      Data (Guess1(13,J),J=2,10)/9*0.0d0/
      Data (Guess2(13,J),J=2,10)/9*0.0d0/
      Data (Guess3(13,J),J=2,10)/9*0.0d0/
C                    DATA FOR ELEMENT 14
      DATA USS   (14)/     -33.9536220D0/
      DATA UPP   (14)/     -28.9347490D0/
      DATA BETAS (14)/      -3.7848520D0/
      DATA BETAP (14)/      -1.9681230D0/
      DATA ZS    (14)/       1.8306970D0/
      DATA ZP    (14)/       1.2849530D0/
      DATA ZD    (14)/       1.0000000D0/
      DATA ALP   (14)/       2.2578160D0/
      DATA EISOL (14)/     -79.0017420D0/
      DATA GSS   (14)/       9.8200000D0/
      DATA GSP   (14)/       8.3600000D0/
      DATA GPP   (14)/       7.3100000D0/
      DATA GP2   (14)/       6.5400000D0/
      DATA HSP   (14)/       1.3200000D0/
      DATA DD    (14)/       1.1631107D0/
      DATA QQ    (14)/       1.3022422D0/
      DATA AM    (14)/       0.3608967D0/
      DATA AD    (14)/       0.3829813D0/
      DATA AQ    (14)/       0.3712106D0/
      DATA GUESS1(14,1)/       0.2500000D0/
      DATA GUESS2(14,1)/       9.0000000D0/
      DATA GUESS3(14,1)/       0.9114530D0/
      DATA GUESS1(14,2)/       0.0615130D0/
      DATA GUESS2(14,2)/       5.0000000D0/
      DATA GUESS3(14,2)/       1.9955690D0/
      DATA GUESS1(14,3)/       0.0207890D0/
      DATA GUESS2(14,3)/       5.0000000D0/
      DATA GUESS3(14,3)/       2.9906100D0/
      Data (Guess1(14,J),J=4,10)/7*0.0d0/
      Data (Guess2(14,J),J=4,10)/7*0.0d0/
      Data (Guess3(14,J),J=4,10)/7*0.0d0/
C                    DATA FOR ELEMENT 15
      DATA USS   (15)/     -56.1433600D0/
      DATA UPP   (15)/     -42.8510800D0/
      DATA BETAS (15)/      -6.7916000D0/
      DATA BETAP (15)/      -6.7916000D0/
      DATA ZS    (15)/       2.1087200D0/
      DATA ZP    (15)/       1.7858100D0/
      DATA ZD    (15)/       1.0000000D0/
      DATA ALP   (15)/       2.4152800D0/
      DATA EISOL (15)/    -152.9599600D0/
      DATA GSS   (15)/      11.5600000D0/
      DATA GSP   (15)/      10.0800000D0/
      DATA GPP   (15)/       8.6400000D0/
      DATA GP2   (15)/       7.6800000D0/
      DATA HSP   (15)/       1.9200000D0/
      DATA DD    (15)/       1.0129699D0/
      DATA QQ    (15)/       0.9370090D0/
      DATA AM    (15)/       0.4248438D0/
      DATA AD    (15)/       0.4882420D0/
      DATA AQ    (15)/       0.4979406D0/
      Data (Guess1(15,J),J=1,10)/10*0.0d0/
      Data (Guess2(15,J),J=1,10)/10*0.0d0/
      Data (Guess3(15,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT 16  (JAN. 1990)
      DATA USS   (16)/     -56.6940560D0/
      DATA UPP   (16)/     -48.7170490D0/
      DATA BETAS (16)/      -3.9205660D0/
      DATA BETAP (16)/      -7.9052780D0/
      DATA ZS    (16)/       2.3665150D0/
      DATA ZP    (16)/       1.6672630D0/
      DATA ZD    (16)/       1.0000000D0/
      DATA ALP   (16)/       2.4616480D0/
      DATA EISOL (16)/    -191.7321930D0/
      DATA GSS   (16)/      11.7863290D0/
      DATA GSP   (16)/       8.6631270D0/
      DATA GPP   (16)/      10.0393080D0/
      DATA GP2   (16)/       7.7816880D0/
      DATA HSP   (16)/       2.5321370D0/
      DATA DD    (16)/       0.9004265D0/
      DATA QQ    (16)/       1.0036329D0/
      DATA AM    (16)/       0.4331617D0/
      DATA AD    (16)/       0.5907115D0/
      DATA AQ    (16)/       0.6454943D0/
      DATA GUESS1(16,1)/      -0.5091950D0/
      DATA GUESS2(16,1)/       4.5936910D0/
      DATA GUESS3(16,1)/       0.7706650D0/
      DATA GUESS1(16,2)/      -0.0118630D0/
      DATA GUESS2(16,2)/       5.8657310D0/
      DATA GUESS3(16,2)/       1.5033130D0/
      DATA GUESS1(16,3)/       0.0123340D0/
      DATA GUESS2(16,3)/      13.5573360D0/
      DATA GUESS3(16,3)/       2.0091730D0/
      Data (Guess1(16,J),J=4,10)/7*0.0d0/
      Data (Guess2(16,J),J=4,10)/7*0.0d0/
      Data (Guess3(16,J),J=4,10)/7*0.0d0/
C                    DATA FOR ELEMENT 17
      DATA USS   (17)/    -111.6139480D0/
      DATA UPP   (17)/     -76.6401070D0/
      DATA BETAS (17)/     -24.5946700D0/
      DATA BETAP (17)/     -14.6372160D0/
      DATA ZS    (17)/       3.6313760D0/
      DATA ZP    (17)/       2.0767990D0/
      DATA ZD    (17)/       1.0000000D0/
      DATA ALP   (17)/       2.9193680D0/
      DATA EISOL (17)/    -372.1984310D0/
      DATA GSS   (17)/      15.0300000D0/
      DATA GSP   (17)/      13.1600000D0/
      DATA GPP   (17)/      11.3000000D0/
      DATA GP2   (17)/       9.9700000D0/
      DATA HSP   (17)/       2.4200000D0/
      DATA DD    (17)/       0.5406286D0/
      DATA QQ    (17)/       0.8057208D0/
      DATA AM    (17)/       0.5523705D0/
      DATA AD    (17)/       0.7693200D0/
      DATA AQ    (17)/       0.6133369D0/
      DATA GUESS1(17,1)/       0.0942430D0/
      DATA GUESS2(17,1)/       4.0000000D0/
      DATA GUESS3(17,1)/       1.3000000D0/
      DATA GUESS1(17,2)/       0.0271680D0/
      DATA GUESS2(17,2)/       4.0000000D0/
      DATA GUESS3(17,2)/       2.1000000D0/
      Data (Guess1(17,J),J=3,10)/8*0.0d0/
      Data (Guess2(17,J),J=3,10)/8*0.0d0/
      Data (Guess3(17,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  30   (MAY 1990)
      DATA USS   (30)/    -21.0400080D0/
      DATA UPP   (30)/    -17.6555740D0/
      DATA BETAS (30)/     -1.9974290D0/
      DATA BETAP (30)/     -4.7581190D0/
      DATA ZS    (30)/      1.9542990D0/
      DATA ZP    (30)/      1.3723650D0/
      DATA ZD    (30)/      1.0000000D0/
      DATA ALP   (30)/      1.4845630D0/
      DATA EISOL (30)/    -30.2800160D0/
      DATA GSS   (30)/     11.8000000D0/
      DATA GSP   (30)/     11.1820180D0/
      DATA GPP   (30)/     13.3000000D0/
      DATA GP2   (30)/     12.9305200D0/
      DATA HSP   (30)/      0.4846060D0/
      DATA DD    (30)/      1.3581113D0/
      DATA QQ    (30)/      1.5457406D0/
      DATA AM    (30)/      0.4336641D0/
      DATA AD    (30)/      0.2317423D0/
      DATA AQ    (30)/      0.2621165D0/
      Data (Guess1(30,J),J=1,10)/10*0.0d0/
      Data (Guess2(30,J),J=1,10)/10*0.0d0/
      Data (Guess3(30,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT  32  (MAY 1990)
      DATA USS   (32)/    -34.1838890D0/
      DATA UPP   (32)/    -28.6408110D0/
      DATA BETAS (32)/     -4.3566070D0/
      DATA BETAP (32)/     -0.9910910D0/
      DATA ZS    (32)/      1.2196310D0/
      DATA ZP    (32)/      1.9827940D0/
      DATA ZD    (32)/      1.0000000D0/
      DATA ALP   (32)/      2.1364050D0/
      DATA EISOL (32)/    -78.7084810D0/
      DATA GSS   (32)/     10.1686050D0/
      DATA GSP   (32)/      8.1444730D0/
      DATA GPP   (32)/      6.6719020D0/
      DATA GP2   (32)/      6.2697060D0/
      DATA HSP   (32)/      0.9370930D0/
      DATA DD    (32)/      1.2472095D0/
      DATA QQ    (32)/      1.0698642D0/
      DATA AM    (32)/      0.3737084D0/
      DATA AD    (32)/      0.3180310D0/
      DATA AQ    (32)/      0.3485612D0/
      Data (Guess1(32,J),J=1,10)/10*0.0d0/
      Data (Guess2(32,J),J=1,10)/10*0.0d0/
      Data (Guess3(32,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT 35
      DATA USS   (35)/    -104.6560630D0/
      DATA UPP   (35)/     -74.9300520D0/
      DATA BETAS (35)/     -19.3998800D0/
      DATA BETAP (35)/      -8.9571950D0/
      DATA ZS    (35)/       3.0641330D0/
      DATA ZP    (35)/       2.0383330D0/
      DATA ZD    (35)/       1.0000000D0/
      DATA ALP   (35)/       2.5765460D0/
      DATA EISOL (35)/    -352.3142087D0/
      DATA GSS   (35)/      15.0364395D0/
      DATA GSP   (35)/      13.0346824D0/
      DATA GPP   (35)/      11.2763254D0/
      DATA GP2   (35)/       9.8544255D0/
      DATA HSP   (35)/       2.4558683D0/
      DATA DD    (35)/       0.8458104D0/
      DATA QQ    (35)/       1.0407133D0/
      DATA AM    (35)/       0.5526071D0/
      DATA AD    (35)/       0.6024598D0/
      DATA AQ    (35)/       0.5307555D0/
      DATA GUESS1(35,1)/       0.0666850D0/
      DATA GUESS2(35,1)/       4.0000000D0/
      DATA GUESS3(35,1)/       1.5000000D0/
      DATA GUESS1(35,2)/       0.0255680D0/
      DATA GUESS2(35,2)/       4.0000000D0/
      DATA GUESS3(35,2)/       2.3000000D0/
      Data (Guess1(35,J),J=3,10)/8*0.0d0/
      Data (Guess2(35,J),J=3,10)/8*0.0d0/
      Data (Guess3(35,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 50
      DATA USS   (50)/     -35.4967410D0/
      DATA UPP   (50)/     -28.0976360D0/
      DATA BETAS (50)/      -3.2350000D0/
      DATA BETAP (50)/      -2.5778900D0/
      DATA ZS    (50)/       2.5993760D0/
      DATA ZP   (50)/        1.6959620D0/
      DATA ZD    (50)/       1.0000000D0/
      DATA ALP   (50)/       1.8369360D0/
      DATA EISOL (50)/     -80.6887540D0/
      DATA GSS   (50)/       9.8000000D0/
      DATA GSP   (50)/       8.3000000D0/
      DATA GPP   (50)/       7.3000000D0/
      DATA GP2   (50)/       6.5000000D0/
      DATA HSP   (50)/       1.3000000D0/
      DATA DD    (50)/       1.1528229D0/
      DATA QQ    (50)/       1.5148019D0/
      DATA AM    (50)/       0.3601617D0/
      DATA AD    (50)/       0.3823775D0/
      DATA AQ    (50)/       0.3400755D0/
      Data (Guess1(50,J),J=1,10)/10*0.0d0/
      Data (Guess2(50,J),J=1,10)/10*0.0d0/
      Data (Guess3(50,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT 53
      DATA USS   (53)/    -103.5896630D0/
      DATA UPP   (53)/     -74.4299970D0/
      DATA BETAS (53)/      -8.4433270D0/
      DATA BETAP (53)/      -6.3234050D0/
      DATA ZS    (53)/       2.1028580D0/
      DATA ZP    (53)/       2.1611530D0/
      DATA ZD    (53)/       1.0000000D0/
      DATA ALP   (53)/       2.2994240D0/
      DATA EISOL (53)/    -346.8642857D0/
      DATA GSS   (53)/      15.0404486D0/
      DATA GSP   (53)/      13.0565580D0/
      DATA GPP   (53)/      11.1477837D0/
      DATA GP2   (53)/       9.9140907D0/
      DATA HSP   (53)/       2.4563820D0/
      DATA DD    (53)/       1.4878778D0/
      DATA QQ    (53)/       1.1887388D0/
      DATA AM    (53)/       0.5527544D0/
      DATA AD    (53)/       0.44975223D0/
      DATA AQ    (53)/       0.4631775D0/
      DATA GUESS1(53,1)/       0.0043610D0/
      DATA GUESS2(53,1)/       2.3000000D0/
      DATA GUESS3(53,1)/       1.8000000D0/
      DATA GUESS1(53,2)/       0.0157060D0/
      DATA GUESS2(53,2)/       3.0000000D0/
      DATA GUESS3(53,2)/       2.2400000D0/
      Data (Guess1(53,J),J=3,10)/8*0.0d0/
      Data (Guess2(53,J),J=3,10)/8*0.0d0/
      Data (Guess3(53,J),J=3,10)/8*0.0d0/
***********************************************************************
*
*    START OF PM3 PARAMETERS
*
***********************************************************************
C                    DATA FOR ELEMENT  1        HYDROGEN
      DATA USSPM3( 1)/     -13.0733210D0/
      DATA BETASP( 1)/      -5.6265120D0/
      DATA ZSPM3 ( 1)/       0.9678070D0/
      DATA ALPPM3( 1)/       3.3563860D0/
      DATA EISOLP( 1)/     -13.0733210D0/
      DATA GSSPM3( 1)/      14.7942080D0/
      DATA AMPM3 ( 1)/       0.5437048D0/
      DATA ADPM3 ( 1)/       0.5437048D0/
      DATA AQPM3 ( 1)/       0.5437048D0/
      DATA GuesP1( 1,1)/       1.1287500D0/
      DATA GuesP2( 1,1)/       5.0962820D0/
      DATA GuesP3( 1,1)/       1.5374650D0/
      DATA GuesP1( 1,2)/      -1.0603290D0/
      DATA GuesP2( 1,2)/       6.0037880D0/
      DATA GuesP3( 1,2)/       1.5701890D0/
      Data (GuesP1(1,J),J=3,10)/8*0.0d0/
      Data (GuesP2(1,J),J=3,10)/8*0.0d0/
      Data (GuesP3(1,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  3        LITHIUM
      DATA USSPM3( 3)/      -4.9414130D0/
      DATA UPPPM3( 3)/      -4.8996300D0/
      DATA BETASP( 3)/      -0.2783390D0/
      DATA BETAPP( 3)/      -2.8900310D0/
      DATA ZSPM3 ( 3)/       1.5000000D0/
      DATA ZPPM3 ( 3)/       0.7000000D0/
      DATA ALPPM3( 3)/       1.2303780D0/
      DATA EISOLP( 3)/      -4.9414130D0/
      DATA GSSPM3( 3)/       4.2382110D0/
      DATA GSPPM3( 3)/      16.2028110D0/
      DATA GPPPM3( 3)/       5.8425930D0/
      DATA GP2PM3( 3)/      13.2881150D0/
      DATA HSPPM3( 3)/       0.0000001D0/
      DATA DDPM3 ( 3)/       0.9204422D0/
      DATA QQPM3 ( 3)/       1.7496355D0/
      DATA AMPM3 ( 3)/       0.1557593D0/
      DATA ADPM3 ( 3)/       0.0016309D0/
      DATA AQPM3 ( 3)/       0.2010255D0/
      DATA GuesP1( 3,1)/       1.7757820D0/
      DATA GuesP2( 3,1)/       5.7070980D0/
      DATA GuesP3( 3,1)/       1.1333310D0/
      DATA GuesP1( 3,2)/      -1.8605030D0/
      DATA GuesP2( 3,2)/       3.8427580D0/
      DATA GuesP3( 3,2)/       1.0267740D0/
      Data (GuesP1(3,J),J=3,10)/8*0.0d0/
      Data (GuesP2(3,J),J=3,10)/8*0.0d0/
      Data (GuesP3(3,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  4        BERYLIUM
      DATA UPPPM3( 4)/     -11.3042430D0/
      DATA BETASP( 4)/      -3.9620530D0/
      DATA BETAPP( 4)/      -2.7806840D0/
      DATA ZSPM3 ( 4)/       0.8774390D0/
      DATA ZPPM3 ( 4)/       1.5087550D0/
      DATA ALPPM3( 4)/       1.5935360D0/
      DATA EISOLP( 4)/     -25.5166530D0/
      DATA GSSPM3( 4)/       9.0128510D0/
      DATA GSPPM3( 4)/       6.5761990D0/
      DATA GPPPM3( 4)/       6.0571820D0/
      DATA GP2PM3( 4)/       9.0052190D0/
      DATA HSPPM3( 4)/       0.5446790D0/
      DATA DDPM3 ( 4)/       1.0090531D0/
      DATA QQPM3 ( 4)/       0.8117586D0/
      DATA AMPM3 ( 4)/       0.3312330D0/
      DATA ADPM3 ( 4)/       0.2908996D0/
      DATA AQPM3 ( 4)/       0.3530008D0/
      DATA GuesP1( 4,1)/       1.6315720D0/
      DATA GuesP2( 4,1)/       2.6729620D0/
      DATA GuesP3( 4,1)/       1.7916860D0/
      DATA GuesP1( 4,2)/      -2.1109590D0/
      DATA GuesP2( 4,2)/       1.9685940D0/
      DATA GuesP3( 4,2)/       1.7558710D0/
      Data (GuesP1(4,J),J=3,10)/8*0.0d0/
      Data (GuesP2(4,J),J=3,10)/8*0.0d0/
      Data (GuesP3(4,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  5        BORON
      DATA USSPM3( 5)/     -28.3539160D0/
      DATA UPPPM3( 5)/     -23.9738780D0/
      DATA BETASP( 5)/      -5.5215820D0/
      DATA BETAPP( 5)/      -5.9611770D0/
      DATA ZSPM3 ( 5)/       0.9885410D0/
      DATA ZPPM3 ( 5)/       1.9094410D0/
      DATA ALPPM3( 5)/       2.4981070D0/
      DATA EISOLP( 5)/     -52.1098330D0/
      DATA GSSPM3( 5)/       8.8104530D0/
      DATA GSPPM3( 5)/      10.2992280D0/
      DATA GPPPM3( 5)/       4.9891680D0/
      DATA GP2PM3( 5)/       7.2899770D0/
      DATA HSPPM3( 5)/       0.8370320D0/
      DATA DDPM3 ( 5)/       0.7633744D0/
      DATA QQPM3 ( 5)/       0.6414154D0/
      DATA AMPM3 ( 5)/       0.3237947D0/
      DATA ADPM3 ( 5)/       0.4075660D0/
      DATA AQPM3 ( 5)/       0.4208854D0/
      DATA GuesP1( 5,1)/       3.0046640D0/
      DATA GuesP2( 5,1)/       6.0288690D0/
      DATA GuesP3( 5,1)/       0.3472230D0/
      DATA GuesP1( 5,2)/       0.0201950D0/
      DATA GuesP2( 5,2)/       6.0106580D0/
      DATA GuesP3( 5,2)/       1.2140670D0/
      Data (GuesP1(5,J),J=3,10)/8*0.0d0/
      Data (GuesP2(5,J),J=3,10)/8*0.0d0/
      Data (GuesP3(5,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  6        CARBON
      DATA USSPM3( 6)/     -47.2703200D0/
      DATA UPPPM3( 6)/     -36.2669180D0/
      DATA BETASP( 6)/     -11.9100150D0/
      DATA BETAPP( 6)/      -9.8027550D0/
      DATA ZSPM3 ( 6)/       1.5650850D0/
      DATA ZPPM3 ( 6)/       1.8423450D0/
      DATA ALPPM3( 6)/       2.7078070D0/
      DATA EISOLP( 6)/    -111.2299170D0/
      DATA GSSPM3( 6)/      11.2007080D0/
      DATA GSPPM3( 6)/      10.2650270D0/
      DATA GPPPM3( 6)/      10.7962920D0/
      DATA GP2PM3( 6)/       9.0425660D0/
      DATA HSPPM3( 6)/       2.2909800D0/
      DATA DDPM3 ( 6)/       0.8332396D0/
      DATA QQPM3 ( 6)/       0.6647750D0/
      DATA AMPM3 ( 6)/       0.4116394D0/
      DATA ADPM3 ( 6)/       0.5885862D0/
      DATA AQPM3 ( 6)/       0.7647667D0/
      DATA GuesP1( 6,1)/       0.0501070D0/
      DATA GuesP2( 6,1)/       6.0031650D0/
      DATA GuesP3( 6,1)/       1.6422140D0/
      DATA GuesP1( 6,2)/       0.0507330D0/
      DATA GuesP2( 6,2)/       6.0029790D0/
      DATA GuesP3( 6,2)/       0.8924880D0/
      Data (GuesP1(6,J),J=3,10)/8*0.0d0/
      Data (GuesP2(6,J),J=3,10)/8*0.0d0/
      Data (GuesP3(6,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  7        NITROGEN
      DATA USSPM3( 7)/     -49.3356720D0/
      DATA UPPPM3( 7)/     -47.5097360D0/
      DATA BETASP( 7)/     -14.0625210D0/
      DATA BETAPP( 7)/     -20.0438480D0/
      DATA ZSPM3 ( 7)/       2.0280940D0/
      DATA ZPPM3 ( 7)/       2.3137280D0/
      DATA ALPPM3( 7)/       2.8305450D0/
      DATA EISOLP( 7)/    -157.6137755D0/
      DATA GSSPM3( 7)/      11.9047870D0/
      DATA GSPPM3( 7)/       7.3485650D0/
      DATA GPPPM3( 7)/      11.7546720D0/
      DATA GP2PM3( 7)/      10.8072770D0/
      DATA HSPPM3( 7)/       1.1367130D0/
      DATA DDPM3 ( 7)/       0.6577006D0/
      DATA QQPM3 ( 7)/       0.5293383D0/
      DATA AMPM3 ( 7)/       0.4375151D0/
      DATA ADPM3 ( 7)/       0.5030995D0/
      DATA AQPM3 ( 7)/       0.7364933D0/
      DATA GuesP1( 7,1)/       1.5016740D0/
      DATA GuesP2( 7,1)/       5.9011480D0/
      DATA GuesP3( 7,1)/       1.7107400D0/
      DATA GuesP1( 7,2)/      -1.5057720D0/
      DATA GuesP2( 7,2)/       6.0046580D0/
      DATA GuesP3( 7,2)/       1.7161490D0/
      Data (GuesP1(7,J),J=3,10)/8*0.0d0/
      Data (GuesP2(7,J),J=3,10)/8*0.0d0/
      Data (GuesP3(7,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  8        OXYGEN
      DATA USSPM3( 8)/     -86.9930020D0/
      DATA UPPPM3( 8)/     -71.8795800D0/
      DATA BETASP( 8)/     -45.2026510D0/
      DATA BETAPP( 8)/     -24.7525150D0/
      DATA ZSPM3 ( 8)/       3.7965440D0/
      DATA ZPPM3 ( 8)/       2.3894020D0/
      DATA ALPPM3( 8)/       3.2171020D0/
      DATA EISOLP( 8)/    -289.3422065D0/
      DATA GSSPM3( 8)/      15.7557600D0/
      DATA GSPPM3( 8)/      10.6211600D0/
      DATA GPPPM3( 8)/      13.6540160D0/
      DATA GP2PM3( 8)/      12.4060950D0/
      DATA HSPPM3( 8)/       0.5938830D0/
      DATA DDPM3 ( 8)/       0.4086173D0/
      DATA QQPM3 ( 8)/       0.5125738D0/
      DATA AMPM3 ( 8)/       0.5790430D0/
      DATA ADPM3 ( 8)/       0.5299630D0/
      DATA AQPM3 ( 8)/       0.8179630D0/
      DATA GuesP1( 8,1)/      -1.1311280D0/
      DATA GuesP2( 8,1)/       6.0024770D0/
      DATA GuesP3( 8,1)/       1.6073110D0/
      DATA GuesP1( 8,2)/       1.1378910D0/
      DATA GuesP2( 8,2)/       5.9505120D0/
      DATA GuesP3( 8,2)/       1.5983950D0/
      Data (GuesP1(8,J),J=3,10)/8*0.0d0/
      Data (GuesP2(8,J),J=3,10)/8*0.0d0/
      Data (GuesP3(8,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT  9        FLUORINE
      DATA USSPM3( 9)/    -110.4353030D0/
      DATA UPPPM3( 9)/    -105.6850470D0/
      DATA BETASP( 9)/     -48.4059390D0/
      DATA BETAPP( 9)/     -27.7446600D0/
      DATA ZSPM3 ( 9)/       4.7085550D0/
      DATA ZPPM3 ( 9)/       2.4911780D0/
      DATA ALPPM3( 9)/       3.3589210D0/
      DATA EISOLP( 9)/    -437.5171690D0/
      DATA GSSPM3( 9)/      10.4966670D0/
      DATA GSPPM3( 9)/      16.0736890D0/
      DATA GPPPM3( 9)/      14.8172560D0/
      DATA GP2PM3( 9)/      14.4183930D0/
      DATA HSPPM3( 9)/       0.7277630D0/
      DATA DDPM3 ( 9)/       0.3125302D0/
      DATA QQPM3 ( 9)/       0.4916328D0/
      DATA AMPM3 ( 9)/       0.3857650D0/
      DATA ADPM3 ( 9)/       0.6768503D0/
      DATA AQPM3 ( 9)/       0.6120047D0/
      DATA GuesP1( 9,1)/      -0.0121660D0/
      DATA GuesP2( 9,1)/       6.0235740D0/
      DATA GuesP3( 9,1)/       1.8568590D0/
      DATA GuesP1( 9,2)/      -0.0028520D0/
      DATA GuesP2( 9,2)/       6.0037170D0/
      DATA GuesP3( 9,2)/       2.6361580D0/
      Data (GuesP1(9,J),J=3,10)/8*0.0d0/
      Data (GuesP2(9,J),J=3,10)/8*0.0d0/
      Data (GuesP3(9,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 11        SODIUM
      DATA USSPM3(11)/      -4.7683450D0/
      DATA UPPPM3(11)/      -4.7703090D0/
      DATA BETASP(11)/      -0.5354170D0/
      DATA BETAPP(11)/       0.0548760D0/
      DATA ZSPM3 (11)/       0.6767940D0/
      DATA ZPPM3 (11)/       0.9876530D0/
      DATA ALPPM3(11)/       1.5133320D0/
      DATA EISOLP(11)/      -4.7683450D0/
      DATA GSSPM3(11)/       5.1531420D0/
      DATA GSPPM3(11)/       3.1778570D0/
      DATA GPPPM3(11)/       3.7335330D0/
      DATA GP2PM3(11)/      11.5386740D0/
      DATA HSPPM3(11)/       0.5455510D0/
      DATA DDPM3 (11)/       2.1443749D0/
      DATA QQPM3 (11)/       1.6942388D0/
      DATA AMPM3 (11)/       0.1893841D0/
      DATA ADPM3 (11)/       0.1850465D0/
      DATA AQPM3 (11)/       0.2057464D0/
      DATA GuesP1(11,1)/       0.4175960D0/
      DATA GuesP2(11,1)/       8.5969780D0/
      DATA GuesP3(11,1)/       1.6731780D0/
      DATA GuesP1(11,2)/       0.4432600D0/
      DATA GuesP2(11,2)/       3.3594960D0/
      DATA GuesP3(11,2)/       2.2315110D0/
      Data (GuesP1(11,J),J=3,10)/8*0.0d0/
      Data (GuesP2(11,J),J=3,10)/8*0.0d0/
      Data (GuesP3(11,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 12        MAGNESIUM
      DATA USSPM3(12)/     -14.6236880D0/
      DATA UPPPM3(12)/     -14.1734600D0/
      DATA BETASP(12)/      -2.0716910D0/
      DATA BETAPP(12)/      -0.5695810D0/
      DATA ZSPM3 (12)/       0.6985520D0/
      DATA ZPPM3 (12)/       1.4834530D0/
      DATA ALPPM3(12)/       1.3291470D0/
      DATA EISOLP(12)/     -22.5530760D0/
      DATA GSSPM3(12)/       6.6943000D0/
      DATA GSPPM3(12)/       6.7939950D0/
      DATA GPPPM3(12)/       6.9104460D0/
      DATA GP2PM3(12)/       7.0908230D0/
      DATA HSPPM3(12)/       0.5433000D0/
      DATA DDPM3 (12)/       1.1403950D0/
      DATA QQPM3 (12)/       1.1279899D0/
      DATA AMPM3 (12)/       0.2460235D0/
      DATA ADPM3 (12)/       0.2695751D0/
      DATA AQPM3 (12)/       0.2767522D0/
      DATA GuesP1(12,1)/       2.1170500D0/
      DATA GuesP2(12,1)/       6.0094770D0/
      DATA GuesP3(12,1)/       2.0844060D0/
      DATA GuesP1(12,2)/      -2.5477670D0/
      DATA GuesP2(12,2)/       4.3953700D0/
      DATA GuesP3(12,2)/       2.0636740D0/
      Data (GuesP1(12,J),J=3,10)/8*0.0d0/
      Data (GuesP2(12,J),J=3,10)/8*0.0d0/
      Data (GuesP3(12,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 13        ALUMINIUM
      DATA USSPM3(13)/     -24.8454040D0/
      DATA UPPPM3(13)/     -22.2641590D0/
      DATA BETASP(13)/      -0.5943010D0/
      DATA BETAPP(13)/      -0.9565500D0/
      DATA ZSPM3 (13)/       1.7028880D0/
      DATA ZPPM3 (13)/       1.0736290D0/
      DATA ZDPM3 (13)/       1.0000000D0/
      DATA ALPPM3(13)/       1.5217030D0/
      DATA EISOLP(13)/     -46.8647630D0/
      DATA GSSPM3(13)/       5.7767370D0/
      DATA GSPPM3(13)/      11.6598560D0/
      DATA GPPPM3(13)/       6.3477900D0/
      DATA GP2PM3(13)/       6.1210770D0/
      DATA HSPPM3(13)/       4.0062450D0/
      DATA DDPM3 (13)/       1.2102799D0/
      DATA QQPM3 (13)/       1.5585645D0/
      DATA AMPM3 (13)/       0.2123020D0/
      DATA ADPM3 (13)/       0.6418584D0/
      DATA AQPM3 (13)/       0.2262838D0/
      DATA GuesP1(13,1)/      -0.4730900D0/
      DATA GuesP2(13,1)/       1.9158250D0/
      DATA GuesP3(13,1)/       1.4517280D0/
      DATA GuesP1(13,2)/      -0.1540510D0/
      DATA GuesP2(13,2)/       6.0050860D0/
      DATA GuesP3(13,2)/       2.5199970D0/
      Data (GuesP1(13,J),J=3,10)/8*0.0d0/
      Data (GuesP2(13,J),J=3,10)/8*0.0d0/
      Data (GuesP3(13,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 14        SILICON
      DATA USSPM3(14)/     -26.7634830D0/
      DATA UPPPM3(14)/     -22.8136350D0/
      DATA BETASP(14)/      -2.8621450D0/
      DATA BETAPP(14)/      -3.9331480D0/
      DATA ZSPM3 (14)/       1.6350750D0/
      DATA ZPPM3 (14)/       1.3130880D0/
      DATA ZDPM3 (14)/       1.0000000D0/
      DATA ALPPM3(14)/       2.1358090D0/
      DATA EISOLP(14)/     -67.7882140D0/
      DATA GSSPM3(14)/       5.0471960D0/
      DATA GSPPM3(14)/       5.9490570D0/
      DATA GPPPM3(14)/       6.7593670D0/
      DATA GP2PM3(14)/       5.1612970D0/
      DATA HSPPM3(14)/       0.9198320D0/
      DATA DDPM3 (14)/       1.3144550D0/
      DATA QQPM3 (14)/       1.2743396D0/
      DATA AMPM3 (14)/       0.1854905D0/
      DATA ADPM3 (14)/       0.3060715D0/
      DATA AQPM3 (14)/       0.4877432D0/
      DATA GuesP1(14,1)/      -0.3906000D0/
      DATA GuesP2(14,1)/       6.0000540D0/
      DATA GuesP3(14,1)/       0.6322620D0/
      DATA GuesP1(14,2)/       0.0572590D0/
      DATA GuesP2(14,2)/       6.0071830D0/
      DATA GuesP3(14,2)/       2.0199870D0/
      Data (GuesP1(14,J),J=3,10)/8*0.0d0/
      Data (GuesP2(14,J),J=3,10)/8*0.0d0/
      Data (GuesP3(14,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 15        PHOSPHORUS
      DATA USSPM3(15)/     -40.4130960D0/
      DATA UPPPM3(15)/     -29.5930520D0/
      DATA BETASP(15)/     -12.6158790D0/
      DATA BETAPP(15)/      -4.1600400D0/
      DATA ZSPM3 (15)/       2.0175630D0/
      DATA ZPPM3 (15)/       1.5047320D0/
      DATA ZDPM3 (15)/       1.0000000D0/
      DATA ALPPM3(15)/       1.9405340D0/
      DATA EISOLP(15)/    -117.9591740D0/
      DATA GSSPM3(15)/       7.8016150D0/
      DATA GSPPM3(15)/       5.1869490D0/
      DATA GPPPM3(15)/       6.6184780D0/
      DATA GP2PM3(15)/       6.0620020D0/
      DATA HSPPM3(15)/       1.5428090D0/
      DATA DDPM3 (15)/       1.0644947D0/
      DATA QQPM3 (15)/       1.1120386D0/
      DATA AMPM3 (15)/       0.2867187D0/
      DATA ADPM3 (15)/       0.4309446D0/
      DATA AQPM3 (15)/       0.3732517D0/
      DATA GuesP1(15,1)/      -0.6114210D0/
      DATA GuesP2(15,1)/       1.9972720D0/
      DATA GuesP3(15,1)/       0.7946240D0/
      DATA GuesP1(15,2)/      -0.0939350D0/
      DATA GuesP2(15,2)/       1.9983600D0/
      DATA GuesP3(15,2)/       1.9106770D0/
      Data (GuesP1(15,J),J=3,10)/8*0.0d0/
      Data (GuesP2(15,J),J=3,10)/8*0.0d0/
      Data (GuesP3(15,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 16        SULPHUR
      DATA USSPM3(16)/     -49.8953710D0/
      DATA UPPPM3(16)/     -44.3925830D0/
      DATA BETASP(16)/      -8.8274650D0/
      DATA BETAPP(16)/      -8.0914150D0/
      DATA ZSPM3 (16)/       1.8911850D0/
      DATA ZPPM3 (16)/       1.6589720D0/
      DATA ZDPM3 (16)/       1.0000000D0/
      DATA ALPPM3(16)/       2.2697060D0/
      DATA EISOLP(16)/    -183.4537395D0/
      DATA GSSPM3(16)/       8.9646670D0/
      DATA GSPPM3(16)/       6.7859360D0/
      DATA GPPPM3(16)/       9.9681640D0/
      DATA GP2PM3(16)/       7.9702470D0/
      DATA HSPPM3(16)/       4.0418360D0/
      DATA DDPM3 (16)/       1.1214313D0/
      DATA QQPM3 (16)/       1.0086488D0/
      DATA AMPM3 (16)/       0.3294622D0/
      DATA ADPM3 (16)/       0.6679118D0/
      DATA AQPM3 (16)/       0.6137472D0/
      DATA GuesP1(16,1)/      -0.3991910D0/
      DATA GuesP2(16,1)/       6.0006690D0/
      DATA GuesP3(16,1)/       0.9621230D0/
      DATA GuesP1(16,2)/      -0.0548990D0/
      DATA GuesP2(16,2)/       6.0018450D0/
      DATA GuesP3(16,2)/       1.5799440D0/
      Data (GuesP1(16,J),J=3,10)/8*0.0d0/
      Data (GuesP2(16,J),J=3,10)/8*0.0d0/
      Data (GuesP3(16,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 17        CHLORINE
      DATA USSPM3(17)/    -100.6267470D0/
      DATA UPPPM3(17)/     -53.6143960D0/
      DATA BETASP(17)/     -27.5285600D0/
      DATA BETAPP(17)/     -11.5939220D0/
      DATA ZSPM3 (17)/       2.2462100D0/
      DATA ZPPM3 (17)/       2.1510100D0/
      DATA ZDPM3 (17)/       1.0000000D0/
      DATA ALPPM3(17)/       2.5172960D0/
      DATA EISOLP(17)/    -315.1949480D0/
      DATA GSSPM3(17)/      16.0136010D0/
      DATA GSPPM3(17)/       8.0481150D0/
      DATA GPPPM3(17)/       7.5222150D0/
      DATA GP2PM3(17)/       7.5041540D0/
      DATA HSPPM3(17)/       3.4811530D0/
      DATA DDPM3 (17)/       0.9175856D0/
      DATA QQPM3 (17)/       0.7779230D0/
      DATA AMPM3 (17)/       0.5885190D0/
      DATA ADPM3 (17)/       0.6814522D0/
      DATA AQPM3 (17)/       0.3643694D0/
      DATA GuesP1(17,1)/      -0.1715910D0/
      DATA GuesP2(17,1)/       6.0008020D0/
      DATA GuesP3(17,1)/       1.0875020D0/
      DATA GuesP1(17,2)/      -0.0134580D0/
      DATA GuesP2(17,2)/       1.9666180D0/
      DATA GuesP3(17,2)/       2.2928910D0/
      Data (GuesP1(17,J),J=3,10)/8*0.0d0/
      Data (GuesP2(17,J),J=3,10)/8*0.0d0/
      Data (GuesP3(17,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 19        POTASSIUM
      DATA USSPM3(19)/      -4.1079950D0/
      DATA UPPPM3(19)/      -1.9126300D0/
      DATA BETASP(19)/      -0.7235020D0/
      DATA BETAPP(19)/       0.0109380D0/
      DATA ZSPM3 (19)/       0.9158730D0/
      DATA ZPPM3 (19)/       1.3430990D0/
      DATA ALPPM3(19)/       4.4452300D0/
      DATA EISOLP(19)/      -4.1079950D0/
      DATA GSSPM3(19)/       4.0502660D0/
      DATA GSPPM3(19)/       2.9879670D0/
      DATA GPPPM3(19)/       2.9983050D0/
      DATA GP2PM3(19)/       4.9936140D0/
      DATA HSPPM3(19)/       0.5653070D0/
      DATA DDPM3 (19)/       1.9524869D0/
      DATA QQPM3 (19)/       1.5794222D0/
      DATA AMPM3 (19)/       0.1488521D0/
      DATA ADPM3 (19)/       0.1983782D0/
      DATA AQPM3 (19)/       0.2164587D0/
      DATA GuesP1(19,1)/       6.0409170D0/
      DATA GuesP2(19,1)/      10.0465540D0/
      DATA GuesP3(19,1)/       1.6005120D0/
      DATA GuesP1(19,2)/       0.7242240D0/
      DATA GuesP2(19,2)/       3.0082750D0/
      DATA GuesP3(19,2)/       2.4643520D0/
      Data (GuesP1(19,J),J=3,10)/8*0.0d0/
      Data (GuesP2(19,J),J=3,10)/8*0.0d0/
      Data (GuesP3(19,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 20        CALCIUM
      DATA USSPM3(20)/     -11.6466770D0/
      DATA UPPPM3(20)/     -11.7003040D0/
      DATA BETASP(20)/      -1.3663020D0/
      DATA BETAPP(20)/      -0.6800840D0/
      DATA ZSPM3 (20)/       0.8224500D0/
      DATA ZPPM3 (20)/       1.2791170D0/
      DATA ALPPM3(20)/       3.2251430D0/
      DATA EISOLP(20)/     -17.7272050D0/
      DATA GSSPM3(20)/       5.5661490D0/
      DATA GSPPM3(20)/      15.8559960D0/
      DATA GPPPM3(20)/       6.7057840D0/
      DATA GP2PM3(20)/       4.8642700D0/
      DATA HSPPM3(20)/       1.0053200D0/
      DATA DDPM3 (20)/       1.9888793D0/
      DATA QQPM3 (20)/       1.6584256D0/
      DATA AMPM3 (20)/       0.2045626D0/
      DATA ADPM3 (20)/       0.2518009D0/
      DATA AQPM3 (20)/       0.4412688D0/
      DATA GuesP1(20,1)/       0.7940250D0/
      DATA GuesP2(20,1)/       3.8548560D0/
      DATA GuesP3(20,1)/       1.5284410D0/
      DATA GuesP1(20,2)/       0.3400220D0/
      DATA GuesP2(20,2)/       5.7455560D0/
      DATA GuesP3(20,2)/       2.3571110D0/
      Data (GuesP1(20,J),J=3,10)/8*0.0d0/
      Data (GuesP2(20,J),J=3,10)/8*0.0d0/
      Data (GuesP3(20,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 24        CHROMIUM (?)
      DATA USSPM3(24)/     -17.5170270D0/
      DATA UPPPM3(24)/     -12.5337290D0/
      DATA UDDPM3(24)/     -44.1249280D0/
      DATA BETASP(24)/      -0.1000000D0/
      DATA BETAPP(24)/      -0.1000000D0/
      DATA BETADP(24)/      -8.7766360D0/
      DATA ZSPM3 (24)/       1.5000000D0/
      DATA ZPPM3 (24)/       1.5000000D0/
      DATA ZDPM3 (24)/       2.8845490D0/
      DATA ALPPM3(24)/       3.0683070D0/
      DATA EISOLP(24)/    -134.8187920D0/
      DATA GSSPM3(24)/       6.0000000D0/
      DATA GSPPM3(24)/       4.1500000D0/
      DATA GPPPM3(24)/       5.0000000D0/
      DATA GP2PM3(24)/       3.5000000D0/
      DATA HSPPM3(24)/       1.0000000D0/
      DATA DDPM3 (24)/       1.7320508D0/
      DATA QQPM3 (24)/       1.4142136D0/
      DATA AMPM3 (24)/       0.2205072D0/
      DATA ADPM3 (24)/       0.2711332D0/
      DATA AQPM3 (24)/       0.4464656D0/
      Data (GuesP1(24,J),J=1,10)/10*0.0d0/
      Data (GuesP2(24,J),J=1,10)/10*0.0d0/
      Data (GuesP3(24,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT 30        ZINC
      DATA USSPM3(30)/     -18.5321980D0/
      DATA UPPPM3(30)/     -11.0474090D0/
      DATA BETASP(30)/      -0.7155780D0/
      DATA BETAPP(30)/      -6.3518640D0/
      DATA ZSPM3 (30)/       1.8199890D0/
      DATA ZPPM3 (30)/       1.5069220D0/
      DATA ZDPM3 (30)/       1.0000000D0/
      DATA ALPPM3(30)/       1.3501260D0/
      DATA EISOLP(30)/     -27.3872000D0/
      DATA GSSPM3(30)/       9.6771960D0/
      DATA GSPPM3(30)/       7.7362040D0/
      DATA GPPPM3(30)/       4.9801740D0/
      DATA GP2PM3(30)/       4.6696560D0/
      DATA HSPPM3(30)/       0.6004130D0/
      DATA DDPM3 (30)/       1.5005758D0/
      DATA QQPM3 (30)/       1.4077174D0/
      DATA AMPM3 (30)/       0.3556485D0/
      DATA ADPM3 (30)/       0.2375689D0/
      DATA AQPM3 (30)/       0.2661069D0/
      DATA GuesP1(30,1)/      -0.1112340D0/
      DATA GuesP2(30,1)/       6.0014780D0/
      DATA GuesP3(30,1)/       1.5160320D0/
      DATA GuesP1(30,2)/      -0.1323700D0/
      DATA GuesP2(30,2)/       1.9958390D0/
      DATA GuesP3(30,2)/       2.5196420D0/
      Data (GuesP1(30,J),J=3,10)/8*0.0d0/
      Data (GuesP2(30,J),J=3,10)/8*0.0d0/
      Data (GuesP3(30,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 31        GALLIUM
      DATA USSPM3(31)/     -29.8555930D0/
      DATA UPPPM3(31)/     -21.8753710D0/
      DATA BETASP(31)/      -4.9456180D0/
      DATA BETAPP(31)/      -0.4070530D0/
      DATA ZSPM3 (31)/       1.8470400D0/
      DATA ZPPM3 (31)/       0.8394110D0/
      DATA ALPPM3(31)/       1.6051150D0/
      DATA EISOLP(31)/     -57.3280250D0/
      DATA GSSPM3(31)/       8.4585540D0/
      DATA GSPPM3(31)/       8.9256190D0/
      DATA GPPPM3(31)/       5.0868550D0/
      DATA GP2PM3(31)/       4.9830450D0/
      DATA HSPPM3(31)/       2.0512600D0/
      DATA DDPM3 (31)/       0.9776692D0/
      DATA QQPM3 (31)/       2.5271534D0/
      DATA AMPM3 (31)/       0.3108620D0/
      DATA ADPM3 (31)/       0.5129360D0/
      DATA AQPM3 (31)/       0.1546208D0/
      DATA GuesP1(31,1)/      -0.5601790D0/
      DATA GuesP2(31,1)/       5.6232730D0/
      DATA GuesP3(31,1)/       1.5317800D0/
      DATA GuesP1(31,2)/      -0.2727310D0/
      DATA GuesP2(31,2)/       1.9918430D0/
      DATA GuesP3(31,2)/       2.1838640D0/
      Data (GuesP1(31,J),J=3,10)/8*0.0d0/
      Data (GuesP2(31,J),J=3,10)/8*0.0d0/
      Data (GuesP3(31,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 32        GERMANIUM
      DATA USSPM3(32)/     -35.4671955D0/
      DATA UPPPM3(32)/     -31.5863583D0/
      DATA BETASP(32)/      -5.3250024D0/
      DATA BETAPP(32)/      -2.2501567D0/
      DATA ZSPM3 (32)/       2.2373526D0/
      DATA ZPPM3 (32)/       1.5924319D0/
      DATA ALPPM3(32)/       1.9723370D0/
      DATA EISOLP(32)/     -84.0156006D0/
      DATA GSSPM3(32)/       5.3769635D0/
      DATA GSPPM3(32)/      10.2095293D0/
      DATA GPPPM3(32)/       7.6718647D0/
      DATA GP2PM3(32)/       6.9242663D0/
      DATA HSPPM3(32)/       1.3370204D0/
      DATA DDPM3 (32)/       1.1920304D0/
      DATA QQPM3 (32)/       1.3321263D0/
      DATA AMPM3 (32)/       0.1976098D0/
      DATA ADPM3 (32)/       0.3798182D0/
      DATA AQPM3 (32)/       0.3620669D0/
      DATA GuesP1(32,1)/       0.9631726D0/
      DATA GuesP2(32,1)/       6.0120134D0/
      DATA GuesP3(32,1)/       2.1633655D0/
      DATA GuesP1(32,2)/      -0.9593891D0/
      DATA GuesP2(32,2)/       5.7491802D0/
      DATA GuesP3(32,2)/       2.1693724D0/
      Data (GuesP1(32,J),J=3,10)/8*0.0d0/
      Data (GuesP2(32,J),J=3,10)/8*0.0d0/
      Data (GuesP3(32,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 33        ARSENIC
      DATA USSPM3(33)/     -38.5074240D0/
      DATA UPPPM3(33)/     -35.1524150D0/
      DATA BETASP(33)/      -8.2321650D0/
      DATA BETAPP(33)/      -5.0173860D0/
      DATA ZSPM3 (33)/       2.6361770D0/
      DATA ZPPM3 (33)/       1.7038890D0/
      DATA ALPPM3(33)/       1.7944770D0/
      DATA EISOLP(33)/    -122.6326140D0/
      DATA GSSPM3(33)/       8.7890010D0/
      DATA GSPPM3(33)/       5.3979830D0/
      DATA GPPPM3(33)/       8.2872500D0/
      DATA GP2PM3(33)/       8.2103460D0/
      DATA HSPPM3(33)/       1.9510340D0/
      DATA DDPM3 (33)/       0.9679655D0/
      DATA QQPM3 (33)/       1.2449874D0/
      DATA AMPM3 (33)/       0.3230063D0/
      DATA ADPM3 (33)/       0.5042239D0/
      DATA AQPM3 (33)/       0.2574219D0/
      DATA GuesP1(33,1)/      -0.4600950D0/
      DATA GuesP2(33,1)/       1.9831150D0/
      DATA GuesP3(33,1)/       1.0867930D0/
      DATA GuesP1(33,2)/      -0.0889960D0/
      DATA GuesP2(33,2)/       1.9929440D0/
      DATA GuesP3(33,2)/       2.1400580D0/
      Data (GuesP1(33,J),J=3,10)/8*0.0d0/
      Data (GuesP2(33,J),J=3,10)/8*0.0d0/
      Data (GuesP3(33,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 34        SELENIUM
      DATA USSPM3(34)/     -55.3781350D0/
      DATA UPPPM3(34)/     -49.8230760D0/
      DATA BETASP(34)/      -6.1578220D0/
      DATA BETAPP(34)/      -5.4930390D0/
      DATA ZSPM3 (34)/       2.8280510D0/
      DATA ZPPM3 (34)/       1.7325360D0/
      DATA ALPPM3(34)/       3.0439570D0/
      DATA EISOLP(34)/    -192.7748115D0/
      DATA GSSPM3(34)/       7.4325910D0/
      DATA GSPPM3(34)/      10.0604610D0/
      DATA GPPPM3(34)/       9.5683260D0/
      DATA GP2PM3(34)/       7.7242890D0/
      DATA HSPPM3(34)/       4.0165580D0/
      DATA DDPM3 (34)/       0.8719813D0/
      DATA QQPM3 (34)/       1.2244019D0/
      DATA AMPM3 (34)/       0.2731566D0/
      DATA ADPM3 (34)/       0.7509697D0/
      DATA AQPM3 (34)/       0.5283737D0/
      DATA GuesP1(34,1)/       0.0478730D0/
      DATA GuesP2(34,1)/       6.0074000D0/
      DATA GuesP3(34,1)/       2.0817170D0/
      DATA GuesP1(34,2)/       0.1147200D0/
      DATA GuesP2(34,2)/       6.0086720D0/
      DATA GuesP3(34,2)/       1.5164230D0/
      Data (GuesP1(34,J),J=3,10)/8*0.0d0/
      Data (GuesP2(34,J),J=3,10)/8*0.0d0/
      Data (GuesP3(34,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 35        BROMINE
      DATA USSPM3(35)/    -116.6193110D0/
      DATA UPPPM3(35)/     -74.2271290D0/
      DATA BETASP(35)/     -31.1713420D0/
      DATA BETAPP(35)/      -6.8140130D0/
      DATA ZSPM3 (35)/       5.3484570D0/
      DATA ZPPM3 (35)/       2.1275900D0/
      DATA ZDPM3 (35)/       1.0000000D0/
      DATA ALPPM3(35)/       2.5118420D0/
      DATA EISOLP(35)/    -352.5398970D0/
      DATA GSSPM3(35)/      15.9434250D0/
      DATA GSPPM3(35)/      16.0616800D0/
      DATA GPPPM3(35)/       8.2827630D0/
      DATA GP2PM3(35)/       7.8168490D0/
      DATA HSPPM3(35)/       0.5788690D0/
      DATA DDPM3 (35)/       0.2759025D0/
      DATA QQPM3 (35)/       0.9970532D0/
      DATA AMPM3 (35)/       0.5859399D0/
      DATA ADPM3 (35)/       0.6755383D0/
      DATA AQPM3 (35)/       0.3823719D0/
      DATA GuesP1(35,1)/       0.9604580D0/
      DATA GuesP2(35,1)/       5.9765080D0/
      DATA GuesP3(35,1)/       2.3216540D0/
      DATA GuesP1(35,2)/      -0.9549160D0/
      DATA GuesP2(35,2)/       5.9447030D0/
      DATA GuesP3(35,2)/       2.3281420D0/
      Data (GuesP1(35,J),J=3,10)/8*0.0d0/
      Data (GuesP2(35,J),J=3,10)/8*0.0d0/
      Data (GuesP3(35,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 37        RUBIDIUM
      DATA USSPM3(37)/      -4.4841970D0/
      DATA UPPPM3(37)/      -3.0722520D0/
      DATA BETASP(37)/      -0.7442240D0/
      DATA BETAPP(37)/      -0.5677060D0/
      DATA ZSPM3 (37)/       0.6801340D0/
      DATA ZPPM3 (37)/       0.9658140D0/
      DATA ALPPM3(37)/       5.7668440D0/
      DATA EISOLP(37)/      -4.4841970D0/
      DATA GSSPM3(37)/       4.9848710D0/
      DATA GSPPM3(37)/       4.9983590D0/
      DATA GPPPM3(37)/       6.9857170D0/
      DATA GP2PM3(37)/       5.5079160D0/
      DATA HSPPM3(37)/       0.6211430D0/
      DATA DDPM3 (37)/       3.2610134D0/
      DATA QQPM3 (37)/       2.6599806D0/
      DATA AMPM3 (37)/       0.1832000D0/
      DATA ADPM3 (37)/       0.1544789D0/
      DATA AQPM3 (37)/       0.3073402D0/
      DATA GuesP1(37,1)/       4.7724480D0/
      DATA GuesP2(37,1)/       5.8039230D0/
      DATA GuesP3(37,1)/       1.4442230D0/
      DATA GuesP1(37,2)/       0.5062010D0/
      DATA GuesP2(37,2)/       2.8303820D0/
      DATA GuesP3(37,2)/       2.6433260D0/
      Data (GuesP1(37,J),J=3,10)/8*0.0d0/
      Data (GuesP2(37,J),J=3,10)/8*0.0d0/
      Data (GuesP3(37,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 38        STRONTIUM
      DATA USSPM3(38)/     -11.4870450D0/
      DATA UPPPM3(38)/      -7.7152160D0/
      DATA BETASP(38)/      -0.7533420D0/
      DATA BETAPP(38)/      -0.5760720D0/
      DATA ZSPM3 (38)/       1.5592340D0/
      DATA ZPPM3 (38)/       1.2129820D0/
      DATA ALPPM3(38)/       4.5905860D0/
      DATA EISOLP(38)/     -16.9990050D0/
      DATA GSSPM3(38)/       5.9750850D0/
      DATA GSPPM3(38)/       4.8869600D0/
      DATA GPPPM3(38)/       7.3407630D0/
      DATA GP2PM3(38)/      15.8140850D0/
      DATA HSPPM3(38)/       4.0512180D0/
      DATA DDPM3 (38)/       2.1011077D0/
      DATA QQPM3 (38)/       2.1179593D0/
      DATA AMPM3 (38)/       0.2195915D0/
      DATA ADPM3 (38)/       0.5137005D0/
      DATA AQPM3 (38)/       0.1752812D0/
      DATA GuesP1(38,1)/       0.4398870D0/
      DATA GuesP2(38,1)/       5.7510950D0/
      DATA GuesP3(38,1)/       1.6731940D0/
      DATA GuesP1(38,2)/       0.3957700D0/
      DATA GuesP2(38,2)/       4.4993190D0/
      DATA GuesP3(38,2)/       2.3281550D0/
      Data (GuesP1(38,J),J=3,10)/8*0.0d0/
      Data (GuesP2(38,J),J=3,10)/8*0.0d0/
      Data (GuesP3(38,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 48        CADMIUM
      DATA USSPM3(48)/     -15.8285840D0/
      DATA UPPPM3(48)/       8.7497950D0/
      DATA BETASP(48)/      -8.5819440D0/
      DATA BETAPP(48)/      -0.6010340D0/
      DATA ZSPM3 (48)/       1.6793510D0/
      DATA ZPPM3 (48)/       2.0664120D0/
      DATA ALPPM3(48)/       1.5253820D0/
      DATA EISOLP(48)/     -22.4502080D0/
      DATA GSSPM3(48)/       9.2069600D0/
      DATA GSPPM3(48)/       8.2315390D0/
      DATA GPPPM3(48)/       4.9481040D0/
      DATA GP2PM3(48)/       4.6696560D0/
      DATA HSPPM3(48)/       1.6562340D0/
      DATA DDPM3 (48)/       1.5982681D0/
      DATA QQPM3 (48)/       1.2432402D0/
      DATA AMPM3 (48)/       0.3383668D0/
      DATA ADPM3 (48)/       0.3570290D0/
      DATA AQPM3 (48)/       0.2820582D0/
      Data (GuesP1(48,J),J=1,10)/10*0.0d0/
      Data (GuesP2(48,J),J=1,10)/10*0.0d0/
      Data (GuesP3(48,J),J=1,10)/10*0.0d0/
C                    DATA FOR ELEMENT 49        INDIUM
      DATA USSPM3(49)/     -26.1762050D0/
      DATA UPPPM3(49)/     -20.0058220D0/
      DATA BETASP(49)/      -2.9933190D0/
      DATA BETAPP(49)/      -1.8289080D0/
      DATA ZSPM3 (49)/       2.0161160D0/
      DATA ZPPM3 (49)/       1.4453500D0/
      DATA ALPPM3(49)/       1.4183850D0/
      DATA EISOLP(49)/     -51.9750470D0/
      DATA GSSPM3(49)/       6.5549000D0/
      DATA GSPPM3(49)/       8.2298730D0/
      DATA GPPPM3(49)/       6.2992690D0/
      DATA GP2PM3(49)/       4.9842110D0/
      DATA HSPPM3(49)/       2.6314610D0/
      DATA DDPM3 (49)/       1.5766241D0/
      DATA QQPM3 (49)/       1.7774563D0/
      DATA AMPM3 (49)/       0.2409004D0/
      DATA ADPM3 (49)/       0.4532655D0/
      DATA AQPM3 (49)/       0.3689812D0/
      DATA GuesP1(49,1)/      -0.3431380D0/
      DATA GuesP2(49,1)/       1.9940340D0/
      DATA GuesP3(49,1)/       1.6255160D0/
      DATA GuesP1(49,2)/      -0.1095320D0/
      DATA GuesP2(49,2)/       5.6832170D0/
      DATA GuesP3(49,2)/       2.8670090D0/
      Data (GuesP1(49,J),J=3,10)/8*0.0d0/
      Data (GuesP2(49,J),J=3,10)/8*0.0d0/
      Data (GuesP3(49,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 50        TIN
      DATA USSPM3(50)/     -34.5501920D0/
      DATA UPPPM3(50)/     -25.8944190D0/
      DATA BETASP(50)/      -2.7858020D0/
      DATA BETAPP(50)/      -2.0059990D0/
      DATA ZSPM3 (50)/       2.3733280D0/
      DATA ZPPM3 (50)/       1.6382330D0/
      DATA ALPPM3(50)/       1.6996500D0/
      DATA EISOLP(50)/     -78.8877790D0/
      DATA GSSPM3(50)/      10.1900330D0/
      DATA GSPPM3(50)/       7.2353270D0/
      DATA GPPPM3(50)/       5.6738100D0/
      DATA GP2PM3(50)/       5.1822140D0/
      DATA HSPPM3(50)/       1.0331570D0/
      DATA DDPM3 (50)/       1.3120038D0/
      DATA QQPM3 (50)/       1.5681814D0/
      DATA AMPM3 (50)/       0.3744959D0/
      DATA ADPM3 (50)/       0.3218163D0/
      DATA AQPM3 (50)/       0.2832529D0/
      DATA GuesP1(50,1)/      -0.1503530D0/
      DATA GuesP2(50,1)/       6.0056940D0/
      DATA GuesP3(50,1)/       1.7046420D0/
      DATA GuesP1(50,2)/      -0.0444170D0/
      DATA GuesP2(50,2)/       2.2573810D0/
      DATA GuesP3(50,2)/       2.4698690D0/
      Data (GuesP1(50,J),J=3,10)/8*0.0d0/
      Data (GuesP2(50,J),J=3,10)/8*0.0d0/
      Data (GuesP3(50,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 51        ANTIMONY
      DATA USSPM3(51)/     -56.4321960D0/
      DATA UPPPM3(51)/     -29.4349540D0/
      DATA BETASP(51)/     -14.7942170D0/
      DATA BETAPP(51)/      -2.8179480D0/
      DATA ZSPM3 (51)/       2.3430390D0/
      DATA ZPPM3 (51)/       1.8999920D0/
      DATA ALPPM3(51)/       2.0343010D0/
      DATA EISOLP(51)/    -148.9382890D0/
      DATA GSSPM3(51)/       9.2382770D0/
      DATA GSPPM3(51)/       5.2776800D0/
      DATA GPPPM3(51)/       6.3500000D0/
      DATA GP2PM3(51)/       6.2500000D0/
      DATA HSPPM3(51)/       2.4244640D0/
      DATA DDPM3 (51)/       1.4091903D0/
      DATA QQPM3 (51)/       1.3521354D0/
      DATA AMPM3 (51)/       0.3395177D0/
      DATA ADPM3 (51)/       0.4589010D0/
      DATA AQPM3 (51)/       0.2423472D0/
      DATA GuesP1(51,1)/       3.0020280D0/
      DATA GuesP2(51,1)/       6.0053420D0/
      DATA GuesP3(51,1)/       0.8530600D0/
      DATA GuesP1(51,2)/      -0.0188920D0/
      DATA GuesP2(51,2)/       6.0114780D0/
      DATA GuesP3(51,2)/       2.7933110D0/
      Data (GuesP1(51,J),J=3,10)/8*0.0d0/
      Data (GuesP2(51,J),J=3,10)/8*0.0d0/
      Data (GuesP3(51,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 52      TELLURIUM
      DATA USSPM3(52)/     -44.9380360D0/
      DATA UPPPM3(52)/     -46.3140990D0/
      DATA BETASP(52)/      -2.6651460D0/
      DATA BETAPP(52)/      -3.8954300D0/
      DATA ZSPM3 (52)/       4.1654920D0/
      DATA ZPPM3 (52)/       1.6475550D0/
      DATA ALPPM3(52)/       2.4850190D0/
      DATA EISOLP(52)/    -168.0945925D0/
      DATA GSSPM3(52)/      10.2550730D0/
      DATA GSPPM3(52)/       8.1691450D0/
      DATA GPPPM3(52)/       7.7775920D0/
      DATA GP2PM3(52)/       7.7551210D0/
      DATA HSPPM3(52)/       3.7724620D0/
      DATA DDPM3 (52)/       0.3484177D0/
      DATA QQPM3 (52)/       1.5593085D0/
      DATA AMPM3 (52)/       0.3768862D0/
      DATA ADPM3 (52)/       1.1960743D0/
      DATA AQPM3 (52)/       0.2184786D0/
      DATA GuesP1(52,1)/       0.0333910D0/
      DATA GuesP2(52,1)/       5.9563790D0/
      DATA GuesP3(52,1)/       2.2775750D0/
      DATA GuesP1(52,2)/      -1.9218670D0/
      DATA GuesP2(52,2)/       4.9732190D0/
      DATA GuesP3(52,2)/       0.5242430D0/
      Data (GuesP1(52,J),J=3,10)/8*0.0d0/
      Data (GuesP2(52,J),J=3,10)/8*0.0d0/
      Data (GuesP3(52,J),J=3,10)/8*0.0d0/
C                    DATA FOR ELEMENT 53        IODINE
      DATA USSPM3(53)/     -96.4540370D0/
      DATA UPPPM3(53)/     -61.0915820D0/
      DATA BETASP(53)/     -14.4942340D0/
      DATA BETAPP(53)/      -5.8947030D0/
      DATA ZSPM3 (53)/       7.0010130D0/
      DATA ZPPM3 (53)/       2.4543540D0/
      DATA ZDPM3 (53)/       1.0000000D0/
      DATA ALPPM3(53)/       1.9901850D0/
      DATA EISOLP(53)/    -288.3160860D0/
      DATA GSSPM3(53)/      13.6319430D0/
      DATA GSPPM3(53)/      14.9904060D0/
      DATA GPPPM3(53)/       7.2883300D0/
      DATA GP2PM3(53)/       5.9664070D0/
      DATA HSPPM3(53)/       2.6300350D0/
      DATA DDPM3 (53)/       0.1581469D0/
      DATA QQPM3 (53)/       1.0467302D0/
      DATA AMPM3 (53)/       0.5009902D0/
      DATA ADPM3 (53)/       1.6699104D0/
      DATA AQPM3 (53)/       0.5153082D0/
      DATA GuesP1(53,1)/      -0.1314810D0/
      DATA GuesP2(53,1)/       5.2064170D0/
      DATA GuesP3(53,1)/       1.7488240D0/
      DATA GuesP1(53,2)/      -0.0368970D0/
      DATA GuesP2(53,2)/       6.0101170D0/
      DATA GuesP3(53,2)/       2.7103730D0/
      Data (GuesP1(53,J),J=3,10)/8*0.0d0/
      Data (GuesP2(53,J),J=3,10)/8*0.0d0/
      Data (GuesP3(53,J),J=3,10)/8*0.0d0/
      END

c =====================================================================
c  This file contains PPs original MINDO program as modified for
c  use in SCF GUESS module in PQS package by JB in June 1997
c
c  ** WARNING **
c  subroutines OSINV (inversion) SDIAG2 (diagonalization) and EIG
c  eliminated as these routines are essentially identical to same
c  routines in <service.f>.
c  subroutine SPUR renamed to SPUR1 to avoid name conflict
c =====================================================================
c
      BLOCK DATA OLDMINDO
      IMPLICIT REAL*8(A-H,O-Z)
C.....BASIC DATA
C.....F0 IS THE AVERAGE OF THE ONE-CENTER REPULSION INTEGRALS (EV)
C.....USS AND UPP ARE THE ONE-CENTER CORE INTEGRALS (EV)
C.....VS AND VP ARE THE VALENCE STATE IONIZATION POTENTIALS (EV)
C.....EISOL AND EHEAT ARE THE ENERGIES AND THE HEATS OF FORMATION OF
C        ISOLATED ATOMS (EV AND KCAL/G-ATOM)
C.....GSS,GPP,GSP,GP2,HSP AND HP2 ARE THE ONE-CENTER COULOMB AND EXCHAN
C        GE INTEGRALS (EV)
C.....ZS AND ZP ARE THE OPTIMIZED ORBITAL EXPONENTS
C.....BETA IS THE PARAMETER FOR RESONANCE INTEGRALS
C.....ALFA IS THE PARAMETER FOR CORE REPULSION FUNCTION
C.....ZET IS  (10./SQRT(3.))*((4.*ZS*ZP)**2.5)/((ZS+ZP)**6)
      COMMON/IKB/ IK,NP,NPL,NRE,NDF,NDE
      COMMON/UNIT/ BOHR,HARTRE,DEBYE,CBM,CKALM,AJOULE,AMDYNE
      COMMON/GG/ F0(9),USS(9),UPP(9),VS(9),VP(9),EISOL(9),EHEAT(9),
     1   GSS(9),GPP(9),GSP(9),GP2(9),HSP(9),HP2(9),ZS(9),ZP(9),
     2   BETA(9,9),ALFA(9,9),ZET(9)
      DATA F0(1),F0(2),F0(3),F0(4),F0(5),F0(6),F0(7),F0(8),F0(9)/12.8480
     1   ,0.01,0.01,0.01,8.9580,10.8330,12.3770,13.9850,16.250/
c ...........................................................................
c  note:  above values 0.01 (formerly undefined) to prevent division by zero
c ...........................................................................
      DATA USS(1),USS(2),USS(3),USS(4),USS(5),USS(6),USS(7),USS(8),USS(
     1  9)/-12.5050,0.0,0.0,0.0,-33.610,-51.790,-66.060,-91.730,-129.860
     2   /
      DATA UPP(1),UPP(2),UPP(3),UPP(4),UPP(5),UPP(6),UPP(7),UPP(8),UPP(
     1  9)/0.0,0.0,0.0,0.0,-25.110,-39.180,-56.40,-78.80,-105.930/
      DATA VS(1),VS(2),VS(3),VS(4),VS(5),VS(6),VS(7),VS(8),VS(9)/
     1   -13.6050,0.0,0.0,0.0,-15.160,-21.340,-27.510,-35.30,-43.70/
      DATA VP(1),VP(2),VP(3),VP(4),VP(5),VP(6),VP(7),VP(8),VP(9)/0.0,
     1   0.0,0.0,0.0,-8.520,-11.540,-14.340,-17.910,-20.890/
      DATA EISOL(1),EISOL(2),EISOL(3),EISOL(4),EISOL(5),EISOL(6),EISOL(
     1   7),EISOL(8),EISOL(9)/-12.5050,0.0,0.0,0.0,-61.70,-119.470,
     2   -187.510,-307.070,-475.0/
      DATA EHEAT(1),EHEAT(2),EHEAT(3),EHEAT(4),EHEAT(5),EHEAT(6),EHEAT(
     1   7),EHEAT(8),EHEAT(9)/52.1020,0.0,0.0,0.0,135.70,170.890,113.0,
     2   59.559,18.860/
      DATA GSS(1),GSS(2),GSS(3),GSS(4),GSS(5),GSS(6),GSS(7),GSS(8),GSS(
     1   9)/12.8480,0.0,0.0,0.0,10.590,12.230,13.590,15.420,16.920/
      DATA GPP(1),GPP(2),GPP(3),GPP(4),GPP(5),GPP(6),GPP(7),GPP(8),GPP(
     1   9)/0.0,0.0,0.0,0.0,8.860,11.080,12.980,14.520,16.710/
      DATA GSP(1),GSP(2),GSP(3),GSP(4),GSP(5),GSP(6),GSP(7),GSP(8),GSP(
     1   9)/0.0,0.0,0.0,0.0,9.560,11.470,12.660,14.480,17.250/
      DATA GP2(1),GP2(2),GP2(3),GP2(4),GP2(5),GP2(6),GP2(7),GP2(8),GP2(
     1   9)/0.0,0.0,0.0,0.0,7.860,9.840,11.590,12.980,14.910/
      DATA HSP(1),HSP(2),HSP(3),HSP(4),HSP(5),HSP(6),HSP(7),HSP(8),HSP(
     1   9)/0.0,0.0,0.0,0.0,1.810,2.430,3.140,3.940,4.830/
      DATA HP2(1),HP2(2),HP2(3),HP2(4),HP2(5),HP2(6),HP2(7),HP2(8),HP2(
     1   9)/0.0,0.0,0.0,0.0,0.50,0.620,0.70,0.770,0.90/
      DATA ZS(1),ZS(2),ZS(3),ZS(4),ZS(5),ZS(6),ZS(7),ZS(8),ZS(9)/1.30,
     1   0.0,0.0,0.0,1.2111560,1.7393910,2.7045460,3.640575,3.111270/
      DATA ZP(1),ZP(2),ZP(3),ZP(4),ZP(5),ZP(6),ZP(7),ZP(8),ZP(9)/0.0,
     1   0.0,0.0,0.0,0.9728260,1.7096450,1.8708390,2.1684480,1.419860/
      DATA BETA(1,1),BETA(1,2),BETA(1,3),BETA(1,4),BETA(1,5),BETA(1,6),
     1   BETA(1,7),BETA(1,8),BETA(1,9)/0.244770,3*0.0,0.185347,0.315011,
     2   0.360776,0.417759,0.195242/
      DATA BETA(2,1),BETA(2,2),BETA(2,3),BETA(2,4),BETA(2,5),BETA(2,6),
     1   BETA(2,7),BETA(2,8),BETA(2,9)/9*0.0/
      DATA BETA(3,1),BETA(3,2),BETA(3,3),BETA(3,4),BETA(3,5),BETA(3,6),
     1   BETA(3,7),BETA(3,8),BETA(3,9)/9*0.0/
      DATA BETA(4,1),BETA(4,2),BETA(4,3),BETA(4,4),BETA(4,5),BETA(4,6),
     1   BETA(4,7),BETA(4,8),BETA(4,9)/9*0.0/
      DATA BETA(5,1),BETA(5,2),BETA(5,3),BETA(5,4),BETA(5,5),BETA(5,6),
     1   BETA(5,7),BETA(5,8),BETA(5,9)/4*0.0,0.151324,0.250031,0.310959,
     2   0.349745,0.219591/
      DATA BETA(6,1),BETA(6,2),BETA(6,3),BETA(6,4),BETA(6,5),BETA(6,6),
     1   BETA(6,7),BETA(6,8),BETA(6,9)/5*0.0,0.419907,0.410886,0.464514,
     2   0.247494/
      DATA BETA(7,1),BETA(7,2),BETA(7,3),BETA(7,4),BETA(7,5),BETA(7,6),
     1   BETA(7,7),BETA(7,8),BETA(7,9)/6*0.0,0.377342,0.458110,0.205347/
      DATA BETA(8,1),BETA(8,2),BETA(8,3),BETA(8,4),BETA(8,5),BETA(8,6),
     1   BETA(8,7),BETA(8,8),BETA(8,9)/7*0.0,0.659407,0.334044/
      DATA BETA(9,1),BETA(9,2),BETA(9,3),BETA(9,4),BETA(9,5),BETA(9,6),
     1   BETA(9,7),BETA(9,8),BETA(9,9)/8*0.0,0.197464/
      DATA ALFA(1,1),ALFA(1,2),ALFA(1,3),ALFA(1,4),ALFA(1,5),ALFA(1,6),
     1   ALFA(1,7),ALFA(1,8),ALFA(1,9)/1.489450,3*0.0,2.090352,1.475836,
     2   0.589380,0.478901,3.771362/
      DATA ALFA(2,1),ALFA(2,2),ALFA(2,3),ALFA(2,4),ALFA(2,5),ALFA(2,6),
     1   ALFA(2,7),ALFA(2,8),ALFA(2,9)/9*0.0/
      DATA ALFA(3,1),ALFA(3,2),ALFA(3,3),ALFA(3,4),ALFA(3,5),ALFA(3,6),
     1   ALFA(3,7),ALFA(3,8),ALFA(3,9)/9*0.0/
      DATA ALFA(4,1),ALFA(4,2),ALFA(4,3),ALFA(4,4),ALFA(4,5),ALFA(4,6),
     1   ALFA(4,7),ALFA(4,8),ALFA(4,9)/9*0.0/
      DATA ALFA(5,1),ALFA(5,2),ALFA(5,3),ALFA(5,4),ALFA(5,5),ALFA(5,6),
     1   ALFA(5,7),ALFA(5,8),ALFA(5,9)/4*0.0,2.280544,2.138291,1.909763,
     2   2.484827,2.862183/
      DATA ALFA(6,1),ALFA(6,2),ALFA(6,3),ALFA(6,4),ALFA(6,5),ALFA(6,6),
     1   ALFA(6,7),ALFA(6,8),ALFA(6,9)/5*0.0,1.371208,1.635259,1.820975,
     2   2.725913/
      DATA ALFA(7,1),ALFA(7,2),ALFA(7,3),ALFA(7,4),ALFA(7,5),ALFA(7,6),
     1   ALFA(7,7),ALFA(7,8),ALFA(7,9)/6*0.0,2.029618,1.873859,2.861667/
      DATA ALFA(8,1),ALFA(8,2),ALFA(8,3),ALFA(8,4),ALFA(8,5),ALFA(8,6),
     1   ALFA(8,7),ALFA(8,8),ALFA(8,9)/7*0.0,1.537190,2.266949/
      DATA ALFA(9,1),ALFA(9,2),ALFA(9,3),ALFA(9,4),ALFA(9,5),ALFA(9,6),
     1   ALFA(9,7),ALFA(9,8),ALFA(9,9)/8*0.0,3.864997/
      DATA ZET(1),ZET(2),ZET(3),ZET(4),ZET(5),ZET(6),ZET(7),ZET(8),ZET(
     1   9)/4*0.0,2.5655658,1.6736355,1.1597128,0.8419165,0.8756051/
      DATA BOHR/0.52917715/
      DATA HARTRE/27.21160555/
      DATA DEBYE/0.39342658/
      DATA CBM/0.117946332E+30/
      DATA CKALM/627.5093314/
      DATA AJOULE/0.22936758/
      DATA AMDYNE/8.238855474/
C  IO Unit Nos.
      DATA IK,NP,NPL,NRE,NDF,NDE /6,1,2,3,20,21/
      END
*/Deck RAD
      DOUBLE PRECISION FUNCTION RAD(I,J,X)
      REAL*8 X(3,*)
      RAD=SQRT((X(1,I)-X(1,J))**2+(X(2,I)-X(2,J))**2+(X(3,I)-X(3,J))**2)
      RETURN
      END
*/Deck IDIA
      INTEGER FUNCTION IDIA(K,L)
      IF(K.GE.L) THEN
        IDIA=(K*(K-1))/2+L
      ELSE
        IDIA=(L*(L-1))/2+K
      ENDIF
      RETURN
      END
*/Deck DRUMB
      SUBROUTINE DRUMB(H,N,TEXT)
      REAL*8 H(*)
      CHARACTER*8 TEXT
      COMMON/IKB/ IK,NP,NPL,NRE,NDF,NDE
      WRITE(IK,1001) TEXT
 1001 FORMAT(//,10X,'***',A8,'***',/)
      MN=1
      DO 900 I=1,N
      MN1=MN+I-1
      WRITE(IK,1000) (H(IJ),IJ=MN,MN1)
      MN=MN1+1
  900 CONTINUE
 1000 FORMAT(/,1X,(10F12.5))
      RETURN
      END
*/Deck HMAT
      SUBROUTINE HMAT(NAtoms,IAN,CH,XC,IFUN,GAM,N,SZOR,H)
      IMPLICIT REAL*8(A-H,O-Z)
cc  ...................................................................
cc  Arguments appear to be
cc  NAtoms  -  number of atoms          (formerly IND)
cc  IAN     -  atomic numbers
cc  CH      -  effective atomic charges
cc  XC      -  Cartesian coordinates    (formerly in common/CX/)
cc  IFUN    -  basis function pointer array
cc  GAM     -  storage for 2-centre 2-electron integrals
cc  N       -  number of basis functions
cc  SZOR    -   ???
cc  H       -  on exit contains H matrix
cc  ...................................................................
      DIMENSION IAN(NAtoms),CH(NAtoms),XC(3,NAtoms),IFUN(NAtoms)
      DIMENSION GAM(NAtoms*(NAtoms+1)/2),H(N*(N+1)/2)
      COMMON/GG/ F0(9),USS(9),UPP(9),VS(9),VP(9),EISOL(9),EHEAT(9),
     1   GSS(9),GPP(9),GSP(9),GP2(9),HSP(9),HP2(9),ZS(9),ZP(9),
     2   BETA(9,9),ALFA(9,9),ZET(9)
      COMMON/UNIT/ BOHR,HARTRE,DEBYE,CBM,CKALM,AJOULE,AMDYNE
      COMMON/IKB/ IK,NP,NPL,NRE,NDF,NDE
      DIMENSION RHO(9),X(3)
      LOGICAL TEST
      COMMON /TE/ TEST
      HALFH=HARTRE/2.0
      DO 7 I=1,9
    7 RHO(I)=HALFH/F0(I)
      IFX=1
      DO 100 I=1,NAtoms
      IFUN(I)=IFX
      IFX=IFX+1
      II=IAN(I)
      CH(I)=II
      IF(II.EQ.1) GOTO 100
      CH(I)=II-2
      IFX=IFX+3
  100 CONTINUE
C.... THE ARRAY IFUN(I) CONTAINS THE SERIAL NUMBER  OF THE FIRST BASIS
C.... FUNCTION (1S OR 2S) FOR ATOM I
      MN=0
C.....TWO-CENTER, TWO-ELECTRON INTEGRALS
      DO 600 I=1,NAtoms
      II=IAN(I)
      DO 600 J=1,I
      JJ=IAN(J)
      MN=MN+1
      IF(I.EQ.J) GOTO 610
      RIJ2=RAD(I,J,XC)**2
      GAM(MN)=1.0/SQRT(RIJ2+(RHO(II)+RHO(JJ))**2)
      GOTO 600
  610 GAM(MN)=0.0d0
  600 CONTINUE
      MN=N*(N+1)/2
      DO 199 I=1,MN
  199 H(I)=0.0d0
      IF(TEST) CALL DRUMB(GAM,NAtoms,' COULOMB')
C.....CORE HAMILTONIAN
C.....H  -- DIAGONAL VALUES
      DO 200 I=1,NAtoms
      II=IAN(I)
      IFI=IFUN(I)
      ZSI=ZS(II)
      ZPI=ZP(II)
      DO 300 J=1,I
      JJ=IAN(J)
      IFJ=IFUN(J)
      ZSJ=ZS(JJ)
      ZPJ=ZP(JJ)
      BETAIJ=BETA(II,JJ)
      IF(II.GT.JJ) BETAIJ=BETA(JJ,II)
      X(1) = XC(1,I) - XC(1,J)
      X(2) = XC(2,I) - XC(2,J)
      X(3) = XC(3,I) - XC(3,J)
      R=RAD(I,J,XC)
      IF(I.NE.J) THEN
        RR=1.0d0/R
        GOTO 320
      ENDIF
C.....THE ATTRACTION OF ALL OTHER CORES FOR THE ORBITALS OF ATOM I
      OTHER=0.0d0
      DO 230 L=1,NAtoms
      IF(L.EQ.I) GOTO 230
      IL=IDIA(I,L)
      OTHER=OTHER-CH(L)*GAM(IL)*SZOR
  230 CONTINUE
      ME=1
      IF(II.GT.1) ME=4
      DO 210 MU=1,ME
      MI=IFI+MU-1
      IF(MU.EQ.1) HM=OTHER+USS(II)/HARTRE
      IF(MU.NE.1) HM=OTHER+UPP(II)/HARTRE
      MN=MI*(MI+1)/2
      H(MN)=HM
  210 CONTINUE
      GOTO 300
C.....H  -- OFF-DIAGONAL VALUES
  320 IF(II.NE.1.OR.JJ.NE.1) GOTO 330
C.....BOTH ATOMS HYDROGEN
      CALL OVL(SX,SD,1,.TRUE.,R,ZSI,ZSJ)
      MN=IDIA(IFI,IFJ)
      TERMI=VS(1)/HARTRE
      TERMJ=TERMI
      H(MN)=BETAIJ*(TERMI+TERMJ)*SX
      GOTO 300
  330 IF(II.NE.1) GOTO 340
C.... I=H, J=B,C,N,O,F
      CALL OVL(SX,SD,2,.TRUE.,R,ZSI,ZSJ)
      CALL OVL(SSZ,SD,3,.TRUE.,R,ZSI,ZPJ)
      TERMI=VS(1)/HARTRE
      DO 331 NU=1,4
      TERMJ=VS(JJ)/HARTRE
      IF(NU.EQ.1) GOTO 332
      SX=RR*X(NU-1)*SSZ
      TERMJ=VP(JJ)/HARTRE
  332 MN=IDIA(IFI,IFJ+NU-1)
      H(MN)=BETAIJ*(TERMI+TERMJ)*SX
  331 CONTINUE
      GOTO 300
  340 IF(JJ.NE.1) GOTO 350
C.... J=H, I=B,C,N,O,F
      CALL OVL(SX,SD,2,.TRUE.,R,ZSJ,ZSI)
      CALL OVL(SSZ,SD,3,.TRUE.,R,ZSJ,ZPI)
      TERMJ=VS(1)/HARTRE
      DO 341 MU=1,4
      TERMI=VS(II)/HARTRE
      IF(MU.EQ.1) GOTO 342
      SX=-RR*X(MU-1)*SSZ
      TERMI=VP(II)/HARTRE
  342 MN=IDIA(IFJ,IFI+MU-1)
      H(MN)=BETAIJ*(TERMI+TERMJ)*SX
  341 CONTINUE
      GOTO 300
C.... BOTH ATOMS SECOND ROW
  350 CALL OVL(SX,SD,4,.TRUE.,R,ZSI,ZSJ)
      CALL OVL(SSZ,SD,5,.TRUE.,R,ZSI,ZPJ)
      CALL OVL(SZS,SD,5,.TRUE.,R,ZSJ,ZPI)
      CALL OVL(SZSZ,SD,6,.TRUE.,R,ZPI,ZPJ)
      CALL OVL(PP,SD,7,.TRUE.,R,ZPI,ZPJ)
      DO 351 MU=1,4
      DO 351 NU=1,4
      FF1=0.0d0
      IF(NU.EQ.MU) FF1=1.0d0
      IF(NU.EQ.1.AND.MU.GT.1) SX=-RR*X(MU-1)*SZS
      IF(MU.EQ.1.AND.NU.GT.1) SX=RR*X(NU-1)*SSZ
      IF(NU.GT.1.AND.MU.GT.1) SX=-X(MU-1)*X(NU-1)*RR**2*(PP+SZSZ)+FF1*PP
      I1=IFI+MU-1
      J1=IFJ+NU-1
      MN=IDIA(I1,J1)
      TERMI=VS(II)/HARTRE
      TERMJ=VS(JJ)/HARTRE
      IF(MU.GT.1) TERMI=VP(II)/HARTRE
      IF(NU.GT.1) TERMJ=VP(JJ)/HARTRE
      H(MN)=BETAIJ*(TERMI+TERMJ)*SX
  351 CONTINUE
  300 CONTINUE
  200 CONTINUE
      IF(TEST) WRITE(IK,400) SZOR
  400 FORMAT(//,6H SZOR=,F4.1,/)
      IF(TEST) CALL DRUMB(H,N,'H-MATRIX')
      RETURN
      END
*/Deck A
      DOUBLE PRECISION FUNCTION A(K)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ARG1/X
      C=EXP(-X)
      A=C/X
      IF(K.EQ.0) GOTO 11
      DO 10 I=1,K
   10 A=(A*FLOAT(I)+C)/X
   11 CONTINUE
      RETURN
      END
*/Deck B
      DOUBLE PRECISION FUNCTION B(K)
      IMPLICIT REAL*8(A-H,O-Z)
      COMMON/ARG2/X
      IO=0
      B=0.0
      ABSX=ABS(X)
      IF(ABSX.GT.3.0) GOTO 15
      IF(ABSX.LE.2.0) GOTO 10
      IF(K.LE.10) GOTO 15
      LAST=15
      GOTO 17
   10 IF(ABSX.LE.1.0) GOTO 11
      IF(K.LE.7) GOTO 15
      LAST=12
      GOTO 17
   11 IF(ABSX.LE.0.5) GOTO 12
      IF(K.LE.5) GOTO 15
      LAST=7
      GOTO 17
   12 IF(ABSX.LE.0.000001) GOTO 22
      LAST=6
      GOTO 17
   15 EXPX=EXP(X)
      EXPMX=1.0/EXPX
      B=(EXPX-EXPMX)/X
      IF(K.EQ.0) GOTO 30
      DO 16 I=1,K
   16 B=(FLOAT(I)*B+(-1.0)**I*EXPX-EXPMX)/X
      GOTO 30
   17 I=K
      DO 20 M=IO,LAST
   20 B=B+(-X)**M*FLOAT(2*MOD(M+I+1,2))/(FACT(M)*FLOAT(M+I+1))
      GOTO 30
   22 I=K
      B=FLOAT(2*MOD(I+1,2))/FLOAT(I+1)
   30 CONTINUE
      RETURN
      END
*/Deck FACT
      DOUBLE PRECISION FUNCTION FACT(M)
      FACT=1.0d0
      DO 30 I=1,M
      FACT=FACT*FLOAT(I)
   30 CONTINUE
      RETURN
      END
*/Deck OVL
      SUBROUTINE OVL(S,SD,I,NDER,R,ZA,ZB)
      IMPLICIT REAL*8(A-H,O-Z)
      LOGICAL NDER
      COMMON/ARG1/P
      COMMON/ARG2/Q
      Z1=0.5*(ZA+ZB)
      Z2=0.5*(ZA-ZB)
      P=Z1*R
      Q=Z2*R
      T=Z2/Z1
      IF(ZA.GE.ZB.OR.I.EQ.2.OR.I.EQ.3.OR.I.EQ.5) GOTO 27
      T=-T
      Q=-Q
      Z2=-Z2
   27 EP=EXP(-P)
      RP=1.0/P
      EQ=EXP(-Q)
      RQ=1.0
      IF(Q.NE.0.0) RQ=1.0/Q
      GOTO(28,30,32,34,36,38,40),I
   28 P1=(P*SQRT(1.0-T**2))**3*0.25
      S=P1*(A(2)*B(0)-A(0)*B(2))
      IF(NDER) GOTO 42
      SD=3.0*S/R-P1*(Z 1 *(A(3)*B(0)-A(1)*B(2))+Z2*(A(2)*B(1)-A(0)*B(3)
     1))
      GOTO 42
   30 P1=P*(1.0-T)*(P*SQRT(1.0-T**2))**3/13.8564064
      S=P1*(A(3)*B(0)-A(2)*B(1)-A(1)*B(2)+A(0)*B(3))
      IF(NDER) GOTO 42
      SD=4.0*S/R-P1*(Z1*(A(4)*B(0)-A(3)*B(1)-A(2)*B(2)+A(1)*B(3))
     1 +Z2*(A(3)*B(1)-A(2)*B(2)-A(1)*B(3)+A(0)*B(4)))
      GOTO 42
   32 P1=P*(1.0-T)*(P*SQRT(1.0-T**2))**3*0.125
      S=P1*(-A(3)*B(1)+A(2)*B(0)+A(1)*B(3)-A(0)*B(2))
      IF(NDER) GOTO 42
      SD=4.0*S/R+P1*(Z1*(A(4)*B(1)-A(3)*B(0)-A(2)*B(3)+A(1)*B(2))
     1 +Z2*(A(3)*B(2)-A(2)*B(1)-A(1)*B(4)+A(0)*B(3)))
      GOTO 42
   34 P1=(P*SQRT(1.0-T**2))**5/48.0
      S=P1*(A(4)*B(0)-2.0*A(2)*B(2)+A(0)*B(4))
      IF(NDER) GOTO 42
      SD=5.0*S/R-P1*(Z1*(A(5)*B(0)-2.0*A(3)*B(2)+A(1)*B(4))+
     1 Z2*(A(4)*B(1)-2.0*A(2)*B(3)+A(0)*B(5)))
      GOTO 42
   36 P1=(P*SQRT(1.0-T**2))**5/27.7128128
      S=P1*(A(3)*(B(0)-B(2))+A(1)*(B(4)-B(2))+B(1)*(A(2)-A(4))+
     1 B(3)*(A(2)-A(0)))
      IF(NDER) GOTO 42
      SD=5.0*S/R-P1*(Z1*(A(4)*(B(0)-B(2))+A(2)*(B(4)-B(2))+
     1 B(1)*(A(3)-A(5))+B(3)*(A(3)-A(1)))+Z2*(A(3)*(B(1)-B(3))+
     2 A(1)*(B(5)-B(3))+B(2)*(A(2)-A(4))+B(4)*(A(2)-A(0))))
      GOTO 42
   38 P1=-(P*SQRT(1.0-T**2))**5/16.0
      S=P1*(B(2)*(A(0)+A(4))-A(2)*(B(0)+B(4)))
      IF(NDER) GOTO 42
      SD=5.0*S/R-P1*(Z1*(B(2)*(A(1)+A(5))-A(3)*(B(0)+B(4)))+
     1 Z2*(B(3)*(A(0)+A(4))-A(2)*(B(1)+B(5))))
      GOTO 42
   40 P1=(P*SQRT(1.0-T**2))**5/32.0
      S=P1*(A(4)*(B(0)-B(2))+A(2)*(B(4)-B(0))+A(0)*(B(2)-B(4)))
      IF(NDER) GOTO 42
      SD=5.0*S/R-P1*(Z1*(A(5)*(B(0)-B(2))+A(3) *(B(4)-B(0))+
     1 A(1)*(B(2)-B(4)))+Z2*(A(4)*(B(1)-B(3))+A(2)*(B(5)-B(1))+
     2 A(0)*(B(3)-B(5))))
   42 IF(.NOT.NDER) SD=SD/R
  100 FORMAT(1X,3I5,3F12.5,2F14.5)
      RETURN
      END
*/Deck HFCLOS
      SUBROUTINE HFCLOS(NAtoms, AtSymb, IAN,    CH,     XC,
     $                  IFUN,   GAM,    IPRNT,  N,      P,
     $                  H,      DIAG,   NOC,    IT,     VANF,
     $                  PAR,    CONLIM, PAR1,   REFPR,  CONA,
     $                  CONB,   KEZD,   MFOCK,  LL,     MM,
     $                  NOSD,   A,      IXCLU,  force,  E,
     $                  VEV,    IErr)
      IMPLICIT REAL*8(A-H,O-Z)
cc  ...................................................................
cc  Arguments appear to be
cc  NAtoms  -  number of atoms          (formerly IND)
cc  AtSymb  -  atomic symbols
cc  IAN     -  atomic numbers
cc  CH      -  effective atomic charges
cc  XC      -  Cartesian coordinates    (formerly in common/CX/)
cc  IFUN    -  basis function pointer array
cc  GAM     -  storage for 2-centre 2-electron integrals
cc  IPRNT   -  print flag
cc             <5 - no print out
cc              5 - print out energy etc.. each SCF cycle
cc             >5 -  ditto + final summary
cc  N       -  number of basis functions
cc  P       -  density matrix/MOs
cc  H       -  MINDO core Hamiltonian
cc             (NOTE: Although P and H are dimensioned as lower triangle
cc                    there should be storage for a full N**2 array)
cc  DIAG    -  orbital energies
cc  NOC     -  number of doubly occupied MOs
cc  IT      -  maximum number SCF iterations
cc  VANF    -  logical flag (indicating presence of F matrix)
cc  PAR     -  real - This is probably the level shift - it should be
cc              about 1-2
cc  CONLIM  -  convergence criterion
cc             (maximum difference in 1-electron energy?)
cc  PAR1    -  real - something to do with damping?   ** UNKNOWN **
cc  REFPR   -  logical flag                           ** UNKNOWN **
cc  CONA    -  real - Fock matrix damping?            ** UNKNOWN **
cc  CONB    -    ditto?                               ** UNKNOWN **
cc  KEZD    -  integer (starts in Hungarian)
cc  MFOCK   -  maximum size of DIIS subspace
cc   (HOWEVER, CAREFUL AS SEEMS TO BE NO DISCARDING OF OLD FOCK MATRICES)
cc  LL      -  integer scratch storage
cc  MM      -    ditto
cc  NOSD    -  logical flag (true if No DIIS is allowed?) ** UNKNOWN **
cc  A       -  Diis matrix
cc  IXCLU   -  integer array - something to do with excited states?
cc  force   -  logical flag for calculation of forces
cc
cc  on exit
cc
cc  E       -  total energy
cc  VEV     -  dipole moment vector
cc  IErr    -  error flag
cc               0 - scf converged
cc              -1 - did not converge
cc  ...................................................................
      LOGICAL VANF,REFPR,NOSD,force
      COMMON/IKB/ IK,NP,NPL,NRE,NDF,NDE
      DIMENSION IAN(NAtoms),CH(NAtoms),XC(3,NAtoms),IFUN(NAtoms),
     $          GAM(NAtoms*(NAtoms+1)/2),IXCLU(7)
      DIMENSION P((N*(N+1))/2),H((N*(N+1))/2)
      DIMENSION DIAG(N),LL(*),MM(*),A(MFOCK+1,MFOCK+1)
      DIMENSION VE(3),VEV(3),VEV1(3),VEV2(3)
      CHARACTER*8 AtSymb(NAtoms)
      Character*256 jobname,MOS
      LOGICAL TEST,SCF
      COMMON /TE/ TEST
      COMMON/GG/ XXX(45),EISOL(9),EHEAT(9),XXXX(153),ALFA(9,9),ZET(9)
      COMMON/UNIT/ BOHR,HARTRE,DEBYE,CBM,CKALM,AJOULE,AMDYNE
      Common /job/jobname,lenJ
C
      ND = (N*(N+1))/2
      SCF=.FALSE.
      NIT=0
      M=0
      MX=0
      IF(VANF) NIT=2
C.....CORE REPULSION -- EMAG
      MN=0
      EMAG=0.0d0
      DO 10 I=1,NAtoms
      II=IAN(I)
      DO 10 J=1,I
      JJ=IAN(J)
      MN=MN+1
      IF(I.EQ.J) GOTO 10
      ALFAIJ=ALFA(II,JJ)
      IF(II.GT.JJ) ALFAIJ=ALFA(JJ,II)
      R = RAD(I,J,XC)
      RIJ=R*BOHR
      FUN3=EXP(-ALFAIJ*RIJ)
      IF((((II.EQ.7).OR.(II.EQ.8)).AND.(JJ.EQ.1)).OR.
     1   ((II.EQ.1).AND.((JJ.EQ.7).OR.(JJ.EQ.8)))) FUN3=ALFAIJ*EXP(-RIJ)
      DJI=GAM(MN)+(1.0d0/R-GAM(MN))*FUN3
      EMAG=EMAG+CH(I)*CH(J)*DJI
   10 CONTINUE
      E1O=0.0d0
      REWIND NDF
      REWIND NDE
      IF(VANF) GOTO 1800
      REWIND NPL
      WRITE(NPL) H
C ========================================================================
C -- for semiempirical proper, see if MOs from previous cycle available
cc      IF(force) THEN
cc        call tstchval('mos-file',iyes)
cc        If(iyes.eq.1) Then
cc          call getchval('mos-file',MOS)
cc        Else
cc          MOS = jobname(1:lenJ)//'.mos'
cc        EndIf
cc        call rmblan2(MOS,256,lenM)
cc        itype = 2
cc        call ReadMOS(N,P,lenM,MOS,itype,IErr)
cc        If(IErr.EQ.0) Then
cc          PAR = 0.0d0
cc          PAR1 = 0.0d0
cc          If(IPRNT.GT.0) write(IK,950)
cc 950      Format(' MINDO MOs read from previous cycle')
cc        EndIf
cc      ENDIF
C ========================================================================
cc      If(IErr.NE.0) CALL EIG(H,P,DIAG,N)
      ncalc=min0(noc+10,n)
      CALL EIG1(H,P,DIAG,N,noc)
      CALL DENS(N,NOC,P,H,P)
  200 NIT=NIT+1
C.....DENSITY MATRIX IN P
      REWIND NP
      READ(NP) H
      E1=0.5*SPUR1(H,P,N)
      CALL FOCK(NAtoms,IAN,DIAG,IFUN,N,GAM,1.0d0,P,H)
C.....NEW FOCK MATRIX IN H
      E2=0.5*SPUR1(H,P,N)
      E=E1+E2+EMAG
      DEL=Abs(E1-E1O)
      E1O=E1
      IF(.NOT.NOSD) GOTO 20
      IJ=0
      DO 204 I=1,N
      DO 204 J=1,I
      IJ=IJ+1
      H(IJ)=H(IJ)-PAR*P(IJ)
      IF(I.EQ.J) H(IJ)=H(IJ)+2.0*PAR
  204 CONTINUE
   20 CONTINUE
      IF(NIT.GE.KEZD.AND.(.NOT.NOSD)) M=M+1
      IF(NOSD.OR.M.LE.0) GOTO 70
      IF(MFOCK.EQ.1) GOTO 65
      IF(M.LE.MFOCK) GOTO 70
      MF1=MFOCK-1
      NFILE=NDF
   30 REWIND NFILE
      DO 60 II=1,MF1
      DO 40 I=1,2
CPP it dies at this point
CPP before label 60 it writes and then tries to read immediately
      READ(NFILE) P
   40 CONTINUE
      DO 50 I=1,2
      BACKSPACE NFILE
   50 CONTINUE
      WRITE(NFILE) P
   60 CONTINUE
      IF(NFILE.EQ.NDE) GOTO 70
      NFILE=NDE
      GOTO 30
   65 REWIND NDF
      REWIND NDE
C.....FILES NDF AND NDE ARE IN THE CORRECT POSITION
   70 CONTINUE
      REWIND NPL
      READ(NPL) P
C.....OLD FOCK MATRIX IN P
      XNORM=0.0d0
      DO 205 K=1,ND
      P(K)=H(K)-P(K)
  205 XNORM=XNORM+P(K)*P(K)
      XNORM=SQRT(XNORM)
C.....NEW ERROR VECTOR IN P, ITS NORM IN XNORM
      If(IPRNT.gt.4) WRITE(IK,201) NIT,E,DEL,XNORM
  201 Format(' Cycle: ',I3,2X,' E = ',F13.7,2X,' Delta: ',D13.5,2X,
     $       ' Norm: ',D13.5)
CPP Try to use a high level shift and no DIIS at first
      IF(NIt.GT.8) THEN
        NOSD=.FALSE.
        PAR1=0.0D0
      END IF
CPP
      IF(NIT.GT.IT.OR.DEL.LE.CONLIM) SCF=.TRUE.
      IF(NOSD.OR.M.LE.0) GOTO 80
      WRITE(NDF) H
      WRITE(NDE) P
C.....NEW FOCK MATRIX AND ERROR VECTOR ARE SAVED
   80 IF(NOSD) GOTO 260
C.....GENERATION OF DIIS MATRIX IN A
      MX=M
      IF(M.GT.MFOCK) MX=MFOCK
      IF(MFOCK.EQ.1) GOTO 260
      MXP1=MX+1
      MXM1=MX-1
      IF(M.GT.0) GOTO 90
      A(1,1)=0.0d0
      GOTO 260
   90 IF(M.LE.MFOCK) GOTO 110
      DO 100 I=2,MX
      II=I+1
      DO 100 J=I,MX
      JJ=J+1
      A(I,J)=A(II,JJ)
  100 A(J,I)=A(I,J)
  110 DO 150 II=1,MX
      BACKSPACE NDE
      IF(II.GT.1) BACKSPACE NDE
      I=MXP1-II+1
      IF(II.GT.1) GOTO 120
      READ(NDE) P
      A(MXP1,MXP1)=XNORM*XNORM
      IF(MX.EQ.1) GOTO 140
      GOTO 150
  120 READ(NDE) H
      S=0.0d0
      DO 130 K=1,ND
  130 S=S+P(K)*H(K)
      A(MXP1,I)=S
      A(I,MXP1)=S
      IF(II.LT.MX) GOTO 150
  140 A(MXP1,1)=-1.0
      A(1,MXP1)=-1.0
  150 CONTINUE
      IF(MX.EQ.1) GOTO 260
      DO 160 I=1,MXM1
      READ(NDE)
  160 CONTINUE
C.....FILE NDE IS IN THE ORIGINAL POSITION
C.....COPY MATRIX A INTO ARRAY P (SAVE DIIS MATRIX)
      IJ=0
      DO 170 I=1,MXP1
      DO 170 J=1,MXP1
      IJ=IJ+1
  170 P(IJ)=A(J,I)
      IF(.NOT.TEST) GOTO 177
      WRITE(IK,174)
  174 FORMAT(/1X,11HDIIS MATRIX)
      DO 175 I=1,MXP1
      WRITE(IK,176) (A(I,J),J=1,I)
  175 CONTINUE
  176 FORMAT(/1X,(5E13.5))
  177 CONTINUE
C.....SCALING DIIS MATRIX BEFORE INVERSION
      H(1)=1.0d0
      DO 180 I=2,MXP1
  180 H(I)=1.0/SQRT(A(I,I))
      DO 190 I=1,MXP1
      DO 190 J=I,MXP1
      IJ=(J-1)*MXP1+I
      JI=(I-1)*MXP1+J
      P(IJ)=H(I)*P(IJ)*H(J)
  190 P(JI)=P(IJ)
C.....DIRECT INVERSION IN THE ITERATIVE SUBSPACE
      TOL=1.0E-11
      CALL OSINV(P,MXP1,DET,TOL,LL,MM)
      IF(DET.NE.0.0) GOTO 210
      If(IPRNT.gt.5) WRITE(IK,206)
  206 FORMAT(1X,23HDIIS MATRIX IS SINGULAR)
      NOSD=.TRUE.
      M=0
      MX=0
      BACKSPACE NDF
      READ(NDF) H
      GOTO 260
  210 XLAM=-P(1)
      DO 220 I=2,MXP1
  220 DIAG(I)=-P(I)*H(I)
      IF(TEST) WRITE(IK,225) (DIAG(I),I=2,MXP1)
  225 FORMAT(1X,17HDIIS COEFFICIENTS,/1X,1(10F10.4,/))
      DO 230 K=1,ND
  230 H(K)=0.0
C...GENERATE THE BEST FOCK MATRIX OF THE SUBSPACE
      REWIND NDF
      DO 250 I=2,MXP1
      READ(NDF) P
      DO 240 K=1,ND
  240 H(K)=H(K)+DIAG(I)*P(K)
  250 CONTINUE
C.....FILE NDF IS IN THE ORIGINAL POSITION
C.....DIIS IS READY, BEST FOCK MATRIX IN H
  260 IF((.NOT.NOSD).OR.NIT.LE.2) GOTO 390
      REWIND NPL
      READ(NPL) P
C.....OLD FOCK MATRIX IN P
      P2=1.0-PAR1
      DO 300 IJ=1,ND
  300 H(IJ)=P2*H(IJ)+PAR1*P(IJ)
  390 CONTINUE
      REWIND NPL
      WRITE(NPL) H
C.....SAVE ACTUAL FOCK MATRIX
c  during SCF iteration calculate only the occupied orbitals +2 for swapping
      if(scf) then
        ncalc=n
      else
        ncalc=min0(noc+2,n)
      end if
 1900 CALL EIG1(H,P,DIAG,N,ncalc)
      If(TEST.OR.(SCF.and.IPRNT.gt.5)) CALL PRINWF(P,DIAG,IXCLU,N,NOC)
      IF(SCF) GOTO 2200
      CALL DENS(N,NOC,P,H,P)
      GOTO 200
 1800 REWIND NPL
      REWIND NRE
      DO 1810 IJ=1,ND
 1810 H(IJ)=0.0d0
      IF(CONA.EQ.0.0) GOTO 1830
      READ(NRE) H
      DO 1820 IJ=1,ND
 1820 H(IJ)=CONA*H(IJ)
 1830 IF(CONB.EQ.0.0) GOTO 1850
      READ(NPL) P
      DO 1840 IJ=1,ND
 1840 H(IJ)=H(IJ)+CONB*P(IJ)
 1850 REWIND NPL
      WRITE(NPL) H
      IF(TEST) CALL DRUMB(H,N,'F ESTIMA')
      GOTO 1900
 2200 IF(.NOT.REFPR) GOTO 2300
      REWIND NPL
      REWIND NRE
      READ(NPL) H
      WRITE(NRE) H
C ---- finished ---------
 2300 CONTINUE
C
c  ..............................................................
      If(.NOT.force) Then
        rewind NDF
        write(NDF) (P(I),I=1,N*N)            ! save all MOs
cc      Else
cc        call WriteMOS(N,N,P,lenM,MOS,itype)
      EndIf
c  ..............................................................
c
      CALL DENS(N,NOC,P,H,P)                 ! form density
      If(IPRNT.GT.5) Then
        CALL DRUMB(P,N,' DENSITY')
        CALL CHARG(NAtoms,P,H,IAN,IFUN)
        WRITE(IK,605)
  605   FORMAT(/,1X,41HATOMIC ELECTRON DENSITIES AND NET CHARGES,/)
      EndIF
C
C...calculate dipole moment
      DO 106 L=1,3
 106  VE(L)=0.0
      J=0
      DO 104 I=1,NAtoms
      IE=4
      IF(IAN(I).EQ.1) IE=1
      BEX=0.0
      DO 105 NU=1,IE
      J=J+1
      JJ=IDIA(J,J)
      BEX=BEX+P(JJ)
 105  CONTINUE
      TO=CH(I)-BEX
cc
      If(IPRNT.GT.5) Then
        WRITE(IK,606) I,AtSymb(I),IAN(I),BEX,TO
 606    FORMAT(1X,I3,1H.,3X,A4,2X,I3,1H.,2F12.5)
      EndIf
cc
      DO 104 L=1,3
      VE(L)=VE(L)+TO*XC(L,I)
      VEV(L)=VE(L)
 104  CONTINUE
      DO 89 L=1,3
      DO 89 I=1,NAtoms
      IF(IAN(I).LE.2) GO TO 89
      II=IAN(I)
      LM=IDIA(IFUN(I),IFUN(I)+L)
      VEV(L)=VEV(L)-P(LM)*ZET(II)
  89  CONTINUE
c
      If(IPRNT.GT.5) Then
        WRITE(IK,607) VE
  607   FORMAT(/,1X,28HDIPOLE FROM ATOMIC DENSITIES,/,3X,3F11.7,3H AU,/)
        DO 620 I=1,3
        VEV1(I)=VEV(I)/DEBYE
        VEV2(I)=VEV(I)/CBM
  620   CONTINUE
        WRITE(IK,621) VEV,VEV1,VEV2
  621   FORMAT(/,1X,15HCOMPLETE DIPOLE,/,3X,3F11.7,3H AU,/,3X,3F11.7,
     1     6H DEBYE,/,3X,3E14.6,5H CB.M,/)
        EAT=0.0
        ATHEAT=0.0
        DO 622 I=1,NAtoms
        II=IAN(I)
        EAT=EAT+EISOL(II)
  622   ATHEAT=ATHEAT+EHEAT(II)
        EAT=EAT/HARTRE
        ATHEAT=ATHEAT/CKALM
        ELEC=E-EMAG
        EBON=E-EAT
        HFORM=EBON+ATHEAT
        EX=E*HARTRE
        ELECX=ELEC*HARTRE
        EMAGX=EMAG*HARTRE
        EBONX=EBON*HARTRE
        HFORMX=HFORM*HARTRE
        EY=E*CKALM
        ELECY=ELEC*CKALM
        EMAGY=EMAG*CKALM
        EBONY=EBON*CKALM
        HFORMY=HFORM*CKALM
        EZ=E/AJOULE
        ELECZ=ELEC/AJOULE
        EMAGZ=EMAG/AJOULE
        EBONZ=EBON/AJOULE
        HFORMZ=HFORM/AJOULE
      WRITE(IK,625) E,EX,EY,EZ,ELEC,ELECX,ELECY,ELECZ,EMAG,EMAGX,EMAGY,
     1     EMAGZ,EBON,EBONX,EBONY,EBONZ,HFORM,HFORMX,HFORMY,HFORMZ
  625 FORMAT(/,19H TOTAL ENERGY     =,F15.8,2X,5HA.U.=,F15.6,2X,3HEV=,
     1   F15.3,2X,10HKCAL/MOLE=,F15.8,2X,6HAJOULE,/,
     2   19H ELECTRONIC ENERGY=,F15.8,2X,5HA.U.=,F15.6,2X,3HEV=,F15.3,
     3   2X,10HKCAL/MOLE=,F15.8,2X,6HAJOULE,/,19H CORE-CORE REPULS.=,
     4   F15.8,2X,5HA.U.=,F15.6,2X,3HEV=,F15.3,2X,10HKCAL/MOLE=,F15.8,
     5   2X,6HAJOULE,/,19H BONDING ENERGY   =,F15.8,2X,5HA.U.=,F15.6,
     6   2X,3HEV=,F15.3,2X,10HKCAL/MOLE=,F15.8,2X,6HAJOULE,/,
     7   19H HEAT OF FORMATION=,F15.8,2X,5HA.U.=,F15.6,2X,3HEV=,F15.3,
     8   2X,10HKCAL/MOLE=,F15.8,2X,6HAJOULE,/)
      EndIf
C  ................................................................
      If(.NOT.force) Then
        rewind NDF
        read(NDF) (P(I),I=1,N*N)        ! restore MOs
      EndIf
C  .................................................................
      IErr = NIT
      If(DEL.LE.CONLIM) IErr = 0        ! converged
      RETURN
      END
*/Deck DENS
      SUBROUTINE DENS(N,NOC,C,F,P)
      IMPLICIT REAL*8(A-H,O-Z)
C -- reformulated March 99    JB   (IXCLU ignored)
C.... RESULT IN P IS THE DENSITY MATRIX.
C.... SIMULATED EQUIVALENCE(P,C) SHOULD BE OBSERVED AT THE CALL
      DIMENSION P(1),C(N,NOC),F(NOC,N)
      LOGICAL TEST
      COMMON /TE/ TEST
c
c -- transpose MO coefficients
      DO 2000 I=1,N
      DO 2000 J=1,NOC
      F(J,I)=C(I,J)
 2000 CONTINUE
c
c -- form density
      IJ = 0
      DO 2100 I=1,N
      DO 2100 J=1,I
      IJ = IJ+1
      VAL = 0.0d0
      DO 2200 K=1,NOC
      VAL = VAL + F(K,I)*F(K,J)
 2200 CONTINUE
      P(IJ) = VAL+VAL
 2100 CONTINUE
c
      IF(TEST) CALL DRUMB(P,N,' DENSITY')
      RETURN
      END
*/Deck PRINWF
      SUBROUTINE PRINWF(C,DIAG,IXCLU,N,NOC)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION C(N,N),DIAG(1),IXCLU(7)
      COMMON/IKB/ IK,NP,NPL,NRE,NDF,NDE
      IEX=0
      DO 10 K=1,7
      IF(IXCLU(K).GT.0) IEX=IEX+1
   10 CONTINUE
      NX=NOC+IEX+4
      IF(NX.GT.N) NX=N
      WRITE(IK,603)
  603 FORMAT(/,1X,37HSCF ORBITAL ENERGIES AND EIGENVECTORS,/,1X,
     1 3HNR.,6X,3HEPS,6X,4HOCC.,8X,12HCOEFFICIENTS,/)
      DO 103 I=1,NX
      IOCC=11
      DO 20 K=1,7
      IF(I.EQ.IXCLU(K)) IOCC=0
   20 CONTINUE
      IF(I.GT.NOC+IEX) IOCC=0
      WRITE(IK,604) I,DIAG(I),IOCC,(C(J,I),J=1,N)
  604 FORMAT(/,1X,I3,F13.6,I5,4X,10F11.6,/,(26X,10F11.6))
  103 CONTINUE
      RETURN
      END
*/Deck SPUR1
      DOUBLE PRECISION FUNCTION SPUR1(A,B,N)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(1),B(1)
C...  SPUR1(A*B)
      IJ=0
      SPUR1=0.0d0
      DO 100 I=1,N
      DO 100 J=1,I
      IJ=IJ+1
      X=A(IJ)*B(IJ)
      IF(I.EQ.J) GOTO 101
      X=X+X
  101 SPUR1=SPUR1+X
  100 CONTINUE
      RETURN
      END
*/Deck CHARG
      SUBROUTINE CHARG(NAtoms,P,CHARGE,IAN,IFUN)
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION P(1),CHARGE(1),IAN(NAtoms),IFUN(NAtoms)
      DO 100 I=1,NAtoms
      IE=IFUN(I)
      IK=IDIA(IE,IE)
      CHARGE(I)=P(IK)
      IF(IAN(I).EQ.1) GOTO 100
      DO 101 J=1,3
      IK=IK+IE+J
      CHARGE(I)=CHARGE(I)+P(IK)
  101 CONTINUE
  100 CONTINUE
      RETURN
      END
*/Deck FOCK
      SUBROUTINE FOCK(NAtoms,IAN,CHARGE,IFUN,N,GAM,SZOR,P,H)
      IMPLICIT REAL*8(A-H,O-Z)
cc  ...................................................................
cc  Arguments appear to be
cc  NAtoms  -  number of atoms          (formerly IND)
cc  IAN     -  atomic numbers
cc  CHARGE  -  effective atomic charges
cc  IFUN    -  basis function pointer array
cc  N       -  number of basis functions
cc  GAM     -  storage for 2-centre 2-electron integrals
cc  SZOR    -   ???
cc  P       -  density matrix
cc  H       -  on input contains 1-electron Hamiltonian
cc             on exit contains Fock matrix
cc  ...................................................................
      DIMENSION IAN(NAtoms),CHARGE(NAtoms),IFUN(NAtoms)
      DIMENSION GAM(NAtoms*(NAtoms+1)/2),P(*),H(*)
      COMMON/GG/ F0(9),USS(9),UPP(9),VS(9),VP(9),EISOL(9),EHEAT(9),GSS(
     1   9),GPP(9),GSP(9),GP2(9),HSP(9),HP2(9),ZS(9),ZP(9),BETA(9,9),
     2   ALFA(9,9),Z(9)
      COMMON/UNIT/ BOHR,HARTRE,DEBYE,CBM,CKALM,AJOULE,AMDYNE
      LOGICAL TEST
      COMMON /TE/ TEST
C.... P IS THE DENSITY MATRIX, H IS THE ONE-ELECRON HAMILTONIAN
C.... RESULT IN H IS THE FOCK MATRIX
      CALL CHARG(NAtoms,P,CHARGE,IAN,IFUN)
      DO 200 I=1,NAtoms
      II=IAN(I)
      ID=IDIA(I,I)
      IE=1
      IF(II.GT.1) IE=4
      IFI=IFUN(I)
      GUG=0.0
      DO 209 J=1,NAtoms
      IF(I.EQ.J) GOTO 209
      IJ=IDIA(I,J)
      GUG=CHARGE(J)*GAM(IJ)+GUG
  209 CONTINUE
C.....F -- DIAGONAL VALUES
      DO 210 MU=1,IE
      MI=IFI+MU-1
      ME=IDIA(MI,MI)
      IF(MU.GT.1) GOTO 10
C.....F(S,S)
      FME=(H(ME)+GUG*SZOR)*HARTRE+0.5*P(ME)*GSS(II)+(CHARGE(I)-P(ME))*
     1   (GSP(II)-0.5*HSP(II))
      H(ME)=FME/HARTRE
      MES=ME
      GOTO 210
C.....F(P,P)
   10 FME=(H(ME)+GUG*SZOR)*HARTRE+0.5*P(ME)*GPP(II)+P(MES)*(GSP(II)-0.5*
     1   HSP(II))+(CHARGE(I)-P(ME)-P(MES))*(GP2(II)-0.5*HP2(II))
      H(ME)=FME/HARTRE
  210 CONTINUE
C.....F -- OFF-DIAGONAL VALUES
      DO 220 J=1,I
      IF(I.NE.J) GOTO 100
      IF(II.LT.2) GOTO 220
      IE1=IE-1
      DO 30 MU=1,IE1
      MI=IFI+MU-1
      MU1=MU+1
      DO 30 NU=MU1,IE
      NI=IFI+NU-1
      ME=IDIA(MI,NI)
      IF(MU.NE.1) GOTO 20
C.....F(S,P)
      H(ME)=(0.5*P(ME)*(3.0*HSP(II)-GSP(II)))/HARTRE
      GOTO 30
C.....F(P,PP)
   20 H(ME)=(0.5*P(ME)*(3.0*HP2(II)-GP2(II)))/HARTRE
   30 CONTINUE
      GOTO 220
  100 JJ=IAN(J)
      IJ=ID-I+J
      JE=1
      IF(JJ.GT.1) JE=4
      IFJ=IFUN(J)
      DO 230 MU=1,IE
      MI=IFI+MU-1
      DO 240 NU=1,JE
      NI=IFJ+NU-1
      IF(NI.GT.MI) GOTO 240
      ME=IDIA(MI,NI)
      H(ME)=H(ME)-0.5*P(ME)*GAM(IJ)*SZOR
  240 CONTINUE
  230 CONTINUE
  220 CONTINUE
  200 CONTINUE
      IF(TEST) CALL DRUMB(H,N,'FOCK MAT')
      RETURN
      END

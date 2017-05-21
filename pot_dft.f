      SUBROUTINE SetDFTblock(inp,iout,dft,aXX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  sets up the coefficients for the various DFT functionals
C  stores them in the DFT common block
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All standard methods include Slater exchange)
C              0 - no dft contribution (this routine should NOT be called)
C              1 - Slater exchange only
C              2 - VWN (Slater exchange + local RPA correlation)
C              3 - VWN5 (ditto but RECOMMENDED Ceperley-Alder fit)
C              4 - Becke 88 exchange
C              5 - Becke 88 + VWN (local correlation)
C              6 - Becke 88 + VWN5 (local correlation)
C              7 - Becke 88 + Perdew 86
C              8 - Becke 88 + Perdew-Wang 91
C              9 - Becke 88 + LYP (correlation)
C             10 - Becke 88 + VWN5 (local) + Perdew 86 (nonlocal)
C             11 - Handy/Cohen exchange
C             12 - Handy/Cohen + VWN
C             13 - Handy/Cohen + VWN5
C             14 - Handy/Cohen + Perdew 86
C             15 - Handy/Cohen + Perdew-Wang 91
C             16 - Handy/Cohen + LYP
C             17 - PW91 exchange + correlation
C             18 - PBE exchange + correlation
C             19 - O3LYP   (VWN5 + OPTX     + LYP  + HF exchange)
C             20 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             21 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             22 - WAH     (VWN5 + Becke 88 + LYP  + HF exchange;
C                           modified B3LYP specifically for NMR)
C             23 - user-defined functional
C -- The following functionals are self-contained and were downloaded from the web
C             24 - B97    (Becke 97 exchange-correlation functional)
C             25 - B97-1  ( ditto, but reparametrized)
C             26 - B97-2  ( ditto, but reparametrized again)
C             27 - HCTH   (Hamprecht-Cohen-Tozer-Handy)
C                  (this is similar in form to B97, but has no HF exchange)
C  aXX     -  amount of Hartree-Fock exchange
C             (passed back for use in integral routines)
C
      INTEGER dft
c
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
C
      Parameter (Zero=0.0d0,One=1.0d0)
C
C  coefficients are
C  aX      -  amount of exact (Hartree-Fock) exchange
C  XS      -  amount of Slater (Xalpha) exchange
C  XVWN    -  amount of VWN local correlation (RPA)
C  XVWN5   -  amount of VWN local correlation (Ceperley-Alder fit)
C  XB88    -  amount of Becke 88 nonlocal exchange
C  XP86L   -  amount of Perdew 86 local correlation
C  XP86NL  -   ditto nonlocal correlation
C  XPW91X  -  amount of Perdew-Wang 91 nonlocal exchange
C  XPW91L  -   ditto local correlation
C  XPW91NL -   ditto nonlocal correlation
C  XLYP    -  amount of LYP nonlocal correlation
C  XOptX   -  amount of Handy-Cohen nonlocal exchange
C  XPBEX   -  amount of PBE nonlocal exchange
C  XPBEC   -   ditto correlation (local & nonlocal)
C  EFS     -  factor for Slater-exchange
C  VFS     -  factor for Slater-potential
C  EFSC    -  factor for Slater-exchange closed-shell
C  VFSC    -  factor for Slater-potential closed-shell
C
C
C  initialize Slater factors
cc    XAlpha method  (with Alpha = 2/3)
cc    EFS = -(9.0d0*Alpha/Four)*(Three/(Four*PI))**Third
cc    VFS = -(Three*Alpha)*(Three/(Four*PI))**Third
cc    EFSC = EFS*(Half**Third)                  ! closed shell
cc    VFSC = VFS*(Half**Third)                  ! closed shell
c
      IF(dft.GT.23) THEN
C
C  SELF-CONTAINED FUNCTIONALS
C
        aX = 0.21d0
        If(dft.EQ.24) Then
          aX = 0.19430d0
        Else If(dft.EQ.27) Then
          aX = Zero
        EndIf
cc
      ELSE
cc
        EFS = -0.930525736349100185d0
        VFS = -1.24070098179880017d0
        EFSC= -0.738558766382022558d0
        VFSC= -0.984745021842696744d0
c
cc        If((dft.GT.10.AND.dft.LT.17).OR.dft.EQ.19) Then   ! JB  Oct 2011
        If(dft.GT.10.AND.dft.LT.17) Then
c -- using OPTX non-local exchange
          EFS = 1.05151d0*EFS
          VFS = 1.05151d0*VFS
          EFSC = 1.05151d0*EFSC
          VFSC = 1.05151d0*VFSC
        EndIf
C
C  initialize XS to one and all other coefficients to zero
C
        aX = Zero
        XS = One
        XVWN = Zero
        XVWN5 = Zero
        XB88 = Zero
        XP86L = Zero
        XP86NL = Zero
        XPW91X = Zero
        XPW91L = Zero
        XPW91NL = Zero
        XLYP = Zero
        XOptX = Zero
        XPBEX = Zero
        XPBEC = Zero
C
C  now initialize coefficients depending on the functional
C
        IF(dft.EQ.2) THEN
c -- svwn
          XVWN = One
        ELSE IF(dft.EQ.3) THEN
c -- svwn5
          XVWN5 = One
        ELSE IF(dft.EQ.4) THEN
c -- b88
          XB88 = One
        ELSE IF(dft.EQ.5) THEN
c -- bvwn
          XB88 = One
          XVWN = One
        ELSE IF(dft.EQ.6) THEN
c -- bvwn5
          XB88 = One
          XVWN5 = One
        ELSE IF(dft.EQ.7) THEN
c -- bp86
          XB88 = One
          XP86L = One
          XP86NL = One
        ELSE IF(dft.EQ.8) THEN
c -- bpw91
          XB88 = One
          XPW91L = One
          XPW91NL = One
        ELSE IF(dft.EQ.9) THEN
c -- blyp
          XB88 = One
          XLYP = One
        ELSE IF(dft.EQ.10) THEN
c -- bvp86 (for COSMO)
          XB88 = One
          XVWN5 = One
          XP86NL = One
        ELSE IF(dft.EQ.11) THEN
c -- optx
          XOptX = One
        ELSE IF(dft.EQ.12) THEN
c -- ovwn
          XOptX = One
          XVWN = One
        ELSE IF(dft.EQ.13) THEN
c -- ovwn5
          XOptX = One
          XVWN5 = One
        ELSE IF(dft.EQ.14) THEN
c -- op86
          XOptX = One
          XP86L = One
          XP86NL = One
        ELSE IF(dft.EQ.15) THEN
c -- opw91
          XOptX = One
          XPW91L = One
          XPW91NL = One
        ELSE IF(dft.EQ.16) THEN
c -- olyp
          XOPtX = One
          XLYP = One
        ELSE IF(dft.EQ.17) THEN
c -- pw91
          XPW91X = One
          XPW91L = One
          XPW91NL = One
        ELSE IF(dft.EQ.18) THEN
c -- pbe
          XPBEX = One
          XPBEC = One
        ELSE IF(dft.EQ.19) THEN
c -- o3lyp
          aX = 0.1161d0
          XS = 0.9262d0
          XVWN5 = 0.19d0
          XOptX = 0.8133d0
          XLYP = 0.81d0
        ELSE IF(dft.EQ.20) THEN
c -- b3lyp
          aX = 0.2d0
          XS = 0.8d0
          XVWN = 0.19d0
          XB88 = 0.72d0
          XLYP = 0.81d0
        ELSE IF(dft.EQ.21) THEN
c -- b3pw91
          ax = 0.2d0
          XS = 0.8d0
          XB88 = 0.72d0
          XPW91L = One
          XPW91NL = 0.81d0
        ELSE IF(dft.EQ.22) THEN
c -- wah
          aX = 0.05d0
          XS = 0.95d0
          XVWN5 = 0.19d0
          XB88 = 0.72d0
          XLYP = 0.81d0
        ELSE IF(dft.EQ.23) THEN
c -- user-defined functional (need to read in coefficients)
          call RdDFTUser(inp)
c ....................................................................
          WRITE(iout,*) '  User-defined density functional:'
          If(aX.NE.Zero)    WRITE(iout,*) ' HF exchange:      ',aX
          If(XS.NE.Zero)    WRITE(iout,*) ' Slater exchange:  ',XS
          If(XB88.NE.Zero)  WRITE(iout,*) ' B88 exchange:     ',XB88
          If(XOptX.NE.Zero) WRITE(iout,*) ' OPTX exchange:    ',XOptX
          If(XVWN.NE.Zero)  WRITE(iout,*) ' VWN correlation:  ',XVWN
          If(XVWN5.NE.Zero) WRITE(iout,*) ' VWN5 correlation: ',XVWN5
          If(XP86L.NE.Zero)
     $        WRITE(iout,*) ' P86 Local correlation: ',XP86L
          If(XP86NL.NE.Zero)
     $        WRITE(iout,*) ' P86 Nonlocal correlation: ',XP86NL
          If(XPW91X.NE.Zero)
     $        WRITE(iout,*) ' PW91 Nonlocal exchange: ',XPW91X
          If(XPW91L.NE.Zero)
     $        WRITE(iout,*) ' PW91 Local correlation: ',XPW91L
          If(XPW91NL.NE.Zero)
     $        WRITE(iout,*) ' PW91 Nonlocal correlation: ',XPW91NL
          If(XPBEX.NE.Zero)
     $        WRITE(iout,*) ' PBE Nonlocal exchange: ',XPBEX
          If(XPBEC.NE.Zero) WRITE(iout,*) ' PBE correlation: ',XPBEC
          If(XLYP.NE.Zero)  WRITE(iout,*) ' LYP correlation:  ',XLYP
c ....................................................................
        ELSE IF(dft.NE.1) THEN
c -- undefined functional
          Call nerror(1,'Setting DFT functional',
     $      'Undefined Density Functional - Contact PQS',0,0)
        ENDIF
cc
      ENDIF
c
      aXX = aX   ! pass back to calling routine
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE RdDFTUser(inp)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Read from very next line in input file (following SCF command)
C  coefficients for user-defined density functional
C  (linear combination of existing functionals)
C
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
C
      Character*90 Char
C
C
C  Read next line in input file
C
      READ(inp,900) Char
      call lowerca2(Char,80)
C
C  first 5 characters MUST be $user
C
      If(Char(1:5).NE.'$user') call nerror(13,'SCF Module',
     $   'user-defined functional and no coefficients found',0,0)
C
C  start looking for coefficients
C  expected form    $user ax=0.17 xs=0.83 xvwn=0.81
C  must be at least 1 blank space between entries
C
      I = 5
 10   CONTINUE
      I = I+1
      If(I.GE.80) RETURN
c
      IF(Char(I:I+1).EQ.'ax') THEN
        DO J=I+3,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+3:J),*) aX
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+1).EQ.'xs') THEN
        DO J=I+3,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+3:J),*) XS
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+4).EQ.'xvwn5') THEN
        DO J=I+6,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+6:J),*) XVWN5
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+3).EQ.'xvwn') THEN
        DO J=I+5,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+5:J),*) XVWN
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+3).EQ.'xb88') THEN
        DO J=I+5,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+5:J),*) XB88
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+4).EQ.'xp86l') THEN
        DO J=I+6,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+6:J),*) XP86L
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+5).EQ.'xp86nl') THEN
        DO J=I+7,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+7:J),*) XP86NL
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+5).EQ.'xpw91x') THEN
        DO J=I+7,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+7:J),*) XPW91X
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+5).EQ.'xpw91l') THEN
        DO J=I+7,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+7:J),*) XPW91L
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+6).EQ.'xpw91nl') THEN
        DO J=I+8,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+8:J),*) XPW91NL
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+3).EQ.'xlyp') THEN
        DO J=I+5,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+5:J),*) XLYP
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+4).EQ.'xoptx') THEN
        DO J=I+6,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+6:J),*) XOptX
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+4).EQ.'xpbex') THEN
        DO J=I+6,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+6:J),*) XPBEX
          GO TO 20
        EndIf
        EndDO
      ELSE IF(Char(I:I+4).EQ.'xpbec') THEN
        DO J=I+6,80
        If(Char(J:J).EQ.' ') Then
          Read(Char(I+6:J),*) XPBEC
          GO TO 20
        EndIf
        EndDO
      ELSE
c --nothing found (?)
        J = I
      ENDIF
c
 20   CONTINUE
      I = J
      GO TO 10
c
  900 Format(A80)
c
      END
c =====================================================================
c
      SUBROUTINE vxcpot(dft,    NP,     thrsh,  IdWt,   Den,
     $                  DenX,   WGHT,   FPot,   PotX,   ETot,
     $                  EL,     EVec)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Evaluates the exchange-correlation potential and determines
C  the exchange-correlation energy for the current batch of grid points
C  ** SPECTRAL VERSION OF FUNCTIONALS **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All standard methods include Slater exchange)
C              0 - no dft contribution (this routine should NOT be called)
C              1 - Slater exchange only
C              2 - VWN (Slater exchange + local RPA correlation)
C              3 - VWN5 (ditto but RECOMMENDED Ceperley-Alder fit)
C              4 - Becke 88 exchange
C              5 - Becke 88 + VWN (local correlation)
C              6 - Becke 88 + VWN5 (local correlation)
C              7 - Becke 88 + Perdew 86
C              8 - Becke 88 + Perdew-Wang 91
C              9 - Becke 88 + LYP (correlation)
C             10 - Becke 88 + VWN5 (local) + Perdew 86 (nonlocal)
C             11 - Handy/Cohen exchange
C             12 - Handy/Cohen + VWN
C             13 - Handy/Cohen + VWN5
C             14 - Handy/Cohen + Perdew 86
C             15 - Handy/Cohen + Perdew-Wang 91
C             16 - Handy/Cohen + LYP
C             17 - PW91 exchange + correlation
C             18 - PBE exchange + correlation
C             19 - O3LYP   (VWN5 + OPTX     + LYP  + HF exchange)
C             20 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             21 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             22 - WAH     (VWN5 + Becke 88 + LYP  + HF exchange;
C                           modified B3LYP specifically for NMR)
C             23 - user-defined functional
C -- The following functionals are self-contained and were downloaded from the web
C             24 - B97    (Becke 97 exchange-correlation functional)
C             25 - B97-1  ( ditto, but reparametrized)
C             26 - B97-2  ( ditto, but reparametrized again)
C             27 - HCTH   (Hamprecht-Cohen-Tozer-Handy)
C                  (this is similar in form to B97, but has no HF exchange)
C  NP      -  number of grid points in this batch
C  thrsh   -  threshold for neglect of contribution
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C  Den     -  density
C  DenX    -  density derivatives
C             ** IMPORTANT - These have been predivided by 2 **
C  WGHT    -  grid quadrature weights
C
C  on exit
C
C  FPot    -  potential at each grid point
C  PotX    -  potential derivatives at each grid point
C  ETot    -  total exchange-correlation energy contribution
C  EL      -  number of electrons over grid
C  EVec    -  vector of exchange-correlation energy contributions
C             per grid point
C             (needed for quadrature weight derivatives in DFT gradient)
C
C
      DIMENSION Den(NP),DenX(3,NP),FPot(NP),PotX(3,NP)
      DIMENSION WGHT(NP),EVec(NP)
      INTEGER dft
C
      parameter (Zero=0.0d0,Third=1.0d0/3.0d0,Half=0.5d0,Four=4.0d0)
cc      parameter (PI=3.14159 26535 89793d0)
C
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
C
C
      DO 20 I=1,NP
      Dens = Den(I)
c
      If(Dens.LT.thrsh) Then
       FPot(I) = Zero
       If(dft.GT.3) Then
         PotX(1,I) = Zero
         PotX(2,I) = Zero
         PotX(3,I) = Zero
       EndIf
       GO TO 20
      EndIf
C
C  Do the Self-Contained functionals first
C
      IF(dft.GT.23) THEN
c
       HalfDen = Half*Dens
       g2 = DenX(1,I)**2 + DenX(2,I)**2 + DenX(3,I)**2
c
       If(dft.EQ.24) Then
c...
c...  Becke 97.
c...
        call b97r(halfden,g2,exb97,dra,dsaa)
        EXC = EXB97
        Pot = dra
        VPot = dsaa+dsaa
cc
       Else If(dft.EQ.25) Then
c...
c...  Becke 97_1.
c...
        call b971r(halfden,g2,exb97_1,dra,dsaa)
        EXC = EXB97_1
        Pot = dra
        VPot = dsaa+dsaa
cc
       Else If(dft.EQ.26) Then
c...
c...  Becke 97_2.
c...
        call b972r(halfden,g2,exb97_2,dra,dsaa)
        EXC = EXB97_2
        Pot = dra
        VPot = dsaa+dsaa
cc
       Else
c...
c...  HCTH.
c...
        call hcth407r(halfden,g2,ehcth,dra,dsaa)
        EXC = ehcth
        Pot = dra
        VPot = dsaa+dsaa
cc
       EndIf
ccc
      ELSE
C
C  get the local Slater exchange
C
       Dens13 = Dens**Third
       EXC = XS*EFSC*Dens13
       Pot = XS*VFSC*Dens13
       If(dft.EQ.1) GO TO 10
C
C  VWN local Correlation
C
       IF(XVWN.NE.Zero) THEN
         CALL VWN(Dens,Jnk,.true.,ECVWN,VCVWN,Jnk)
         EXC = EXC + ECVWN*XVWN
         Pot = Pot + VCVWN*XVWN
         If(dft.EQ.2) GO TO 10
       ELSE IF(XVWN5.NE.Zero) THEN
         CALL VWN5(Dens,Jnk,.true.,ECVWN,VCVWN,Jnk)
         EXC = EXC + ECVWN*XVWN5
         Pot = Pot + VCVWN*XVWN5
         If(dft.EQ.3) GO TO 10
       ENDIF
c
       HalfDen = Half*Dens
       Dens13 = HalfDen**Third
       Dens43 = HalfDen*Dens13
c
       g2 = DenX(1,I)**2 + DenX(2,I)**2 + DenX(3,I)**2
       g1 = SQRT(g2)
c
       VPot = Zero
C
C  Becke non-local exchange
C
       IF(XB88.NE.Zero) THEN
         CALL Becke88(HalfDen,Dens13, Dens43, g1,    EXB88,
     $                V1A,    V2A)
         EXC = EXC + EXB88*XB88
         Pot = Pot + V1A*XB88
         VPot = VPot + V2A*XB88
       ENDIF
C
C  Handy/Cohen optimized non-local exchange
C
       IF(XOptX.NE.Zero) THEN
         CALL OPTX(HalfDen,Dens13, Dens43, g1,    EXHC,
     $             V1A,    V2A)
         EXC = EXC + EXHC*XOptX
         Pot = Pot + V1A*XOptX
         VPot = VPot + V2A*XOptX
       ENDIF
c
c  Perdew 86 non-local correlation
c
       IF(XP86L.NE.Zero.OR.XP86NL.NE.Zero) THEN
         call p86nlcr(halfden,g2,xp86l,xp86nl,ecp86,dra,dsaa)
         EXC = EXC + ECP86
         Pot = Pot + dra
         VPot = VPot + four*dsaa
       ENDIF
C
C  Perdew-Wang 91 non-local exchange
C
       IF(XPW91X.NE.Zero) THEN
         call pw91x(halfden,g2,expw91,dra,dsaa)
         EXC = EXC + expw91*xpw91x
         Pot = Pot + dra*xpw91x
         VPot = VPot + (dsaa+dsaa)*xpw91x
       ENDIF
C
C  Perdew-Wang 91 non-local correlation
C
       IF(XPW91L.NE.Zero.OR.XPW91NL.NE.Zero) THEN
c        g2 = Four*g2
c        CALL dftacg_pw91c(.true., .false.,Dens,   Zero,  g2,
c    $                     Zero,   Zero,  ECPW91, VC1A,   jnk,
c    $                     VC2A,   jnk,   jnk)
c        VC2A = Four*VC2A
c        EXC = EXC + ECPW91*XPW91
c        Pot = Pot + VC1A*XPW91
c        VPot = VPot + VC2A*XPW91
cc         write(6,*) ' About to call <pw91nlcr>'
cc         write(6,*) ' halfden:',halfden,' g2:',g2,
cc     $              ' XPW91L:',xpw91l,' XPW91NL:',xpw91nl
         call pw91nlcr(halfden,g2,XPW91L,XPW91NL,ecpw91,dra,dsaa)
cc         write(6,*) ' ecpw91:',ecpw91,' dra:',dra,' dsaa:',dsaa
         EXC = EXC + ECPW91
         Pot = Pot + dra
         VPot = VPot + four*dsaa
       ENDIF
C
C  Lee, Yang & Parr non-local correlation
C
       IF(XLYP.NE.Zero) THEN
         CALL LYP(HalfDen,Dens13, Dens43, g2,     ECLYP,
     $            VC1A,   VC2A)
         EXC = EXC + ECLYP*XLYP
         Pot = Pot + VC1A*XLYP
         VPot = VPot + VC2A*XLYP
       ENDIF
C
C  Perdew, Burke and Ernzerhof non local exchange
C
       IF(XPBEX.NE.Zero) THEN
         call pbexu(halfden,halfden,g2,g2,excpbe,dra,drb,dsaa,dsbb)
         EXC = EXC + excpbe*xpbex
         Pot = Pot + dra*xpbex
         VPot = VPot + (dsaa+dsaa)*xpbex
       ENDIF
C
C  Perdew, Burke and Ernzerhof non local correlation
C
       IF(XPBEC.NE.Zero) THEN
         call pbecu(halfden,halfden,g2,g2,g2,ecpbe,dra,drb,
     +              dsaa,dsbb,dsab)
         EXC = EXC + ecpbe*xpbec
         Pot = Pot + dra*xpbec
         VPot = VPot + (dsaa+dsaa+dsab)*xpbec
       ENDIF
      ENDIF
C
      PotX(1,I) = VPot*DenX(1,I)
      PotX(2,I) = VPot*DenX(2,I)
      PotX(3,I) = VPot*DenX(3,I)
C
 10   CONTINUE
C
C  Accumulate into total Exchange-Correlation energy
C  and total number of electrons
C
      DVal = Den(I)*WGHT(I)
      ETot = ETot + EXC*DVal
      EL   = EL + DVal
C
C  store full potential and exchange-correlation energy
C
      FPot(I) = Pot
      EVec(I) = EXC                           ! for Malkin
      If(IdWt.EQ.1) EVec(I) = EXC*DEN(I)      ! for weight derivatives
cc
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE vxcpotu(dft,    NP,     thrsh,  DenA,   DenB,
     $                   DenAX,  DenBX,  WGHT,   FPotA,  FPotB,
     $                   PotAX,  PotBX,  ETot,   EL,     EVec)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Evaluates the exchange-correlation potential and determines
C  the exchange-correlation energy for the current batch of grid points
C  ** SPECTRAL VERSION OF FUNCTIONALS **
C
C  ARGUMENTS
C
C  dft     -  method flag (NOTE: All standard methods include Slater exchange)
C              0 - no dft contribution (this routine should NOT be called)
C              1 - Slater exchange only
C              2 - VWN (Slater exchange + local RPA correlation)
C              3 - VWN5 (ditto but RECOMMENDED Ceperley-Alder fit)
C              4 - Becke 88 exchange
C              5 - Becke 88 + VWN (local correlation)
C              6 - Becke 88 + VWN5 (local correlation)
C              7 - Becke 88 + Perdew 86
C              8 - Becke 88 + Perdew-Wang 91
C              9 - Becke 88 + LYP (correlation)
C             10 - Becke 88 + VWN5 (local) + Perdew 86 (nonlocal)
C             11 - Handy/Cohen exchange
C             12 - Handy/Cohen + VWN
C             13 - Handy/Cohen + VWN5
C             14 - Handy/Cohen + Perdew 86
C             15 - Handy/Cohen + Perdew-Wang 91
C             16 - Handy/Cohen + LYP
C             17 - PW91 exchange + correlation
C             18 - PBE exchange + correlation
C             19 - O3LYP   (VWN5 + OPTX     + LYP  + HF exchange)
C             20 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             21 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             22 - WAH     (VWN5 + Becke 88 + LYP  + HF exchange;
C                           modified B3LYP specifically for NMR)
C             23 - user-defined functional
C -- The following functionals are self-contained and were downloaded from the web
C             24 - B97    (Becke 97 exchange-correlation functional)
C             25 - B97-1  ( ditto, but reparametrized)
C             26 - B97-2  ( ditto, but reparametrized again)
C             27 - HCTH   (Hamprecht-Cohen-Tozer-Handy)
C                  (this is similar in form to B97, but has no HF exchange)
C  NP      -  number of grid points in this batch
C  thrsh   -  threshold for neglecting contribution
C  DenA    -  alpha density
C  DenB    -  beta density
C  DenAX   -  alpha density derivatives
C  DenBX   -  beta density derivatives
C             ** IMPORTANT - These have been predivided by 2 **
C  WGHT    -  grid quadrature weights
C
C  on exit
C
C  FPotA   -  alpha potential at each grid point
C  FPotB   -  beta  potential at each grid point
C  PotAX   -  alpha potential derivatives at each grid point
C  PotBX   -  beta  potential derivatives at each grid point
C  ETot    -  total exchange-correlation energy contribution
C  EL      -  number of electrons over grid
C  EVec    -  vector of exchange-correlation energy contributions
C             per grid point
C             (needed for quadrature weight derivatives in DFT gradient)
C
C
      DIMENSION DenA(NP),DenAX(3,NP),FPotA(NP),PotAX(3,NP)
      DIMENSION DenB(NP),DenBX(3,NP),FPotB(NP),PotBX(3,NP)
      DIMENSION WGHT(NP),EVec(NP)
      INTEGER dft
C
      parameter (Zero=0.0d0,Third=1.0d0/3.0d0,Half=0.5d0,Four=4.0d0)
cc      parameter (PI=3.14159 26535 89793d0)
C
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
C
C
      DO 20 I=1,NP
      DA = DenA(I)
      DB = DenB(I)
c
      If(DA+DB.LT.thrsh) Then
       FPotA(I) = Zero
       FPotB(I) = Zero
       If(dft.GT.3) Then
         PotAX(1,I) = Zero
         PotAX(2,I) = Zero
         PotAX(3,I) = Zero
         PotBX(1,I) = Zero
         PotBX(2,I) = Zero
         PotBX(3,I) = Zero
       EndIf
       GO TO 20
      EndIf
C
C  Do the Self-Contained functionals first
C
      IF(dft.GT.23) THEN
c
       VPotAB = Zero
       g2A = DenAX(1,I)**2 + DenAX(2,I)**2 + DenAX(3,I)**2
       g2B = DenBX(1,I)**2 + DenBX(2,I)**2 + DenBX(3,I)**2
c
       If(dft.EQ.24) Then
c...
c...  Becke 97.
c...
        call b97u(da,db,g2a,g2b,exb97,dra,drb,dsaa,dsbb)
        EXC = EXB97
        PotA = dra
        PotB = drb
        VPotA = dsaa+dsaa
        VPotB = dsbb+dsbb
cc
       Else If(dft.EQ.25) Then
c...
c...  Becke 97_1.
c...
        call b971u(da,db,g2a,g2b,exb97_1,dra,drb,dsaa,dsbb)
        EXC = EXB97_1
        PotA = dra
        PotB = drb
        VPotA = dsaa+dsaa
        VPotB = dsbb+dsbb
cc
       Else If(dft.EQ.26) Then
c...
c...  Becke 97_2.
c...
        call b972u(da,db,g2a,g2b,exb97_2,dra,drb,dsaa,dsbb)
        EXC = EXB97_2
        PotA = dra
        PotB = drb
        VPotA = dsaa+dsaa
        VPotB = dsbb+dsbb
cc
       Else
c...
c...  HCTH.
c...
        call hcth407u(da,db,g2a,g2b,ehcth,dra,drb,dsaa,dsbb)
        EXC = ehcth
        PotA = dra
        PotB = drb
        VPotA = dsaa+dsaa
        VPotB = dsbb+dsbb
cc
       EndIf
ccc
      ELSE
C
C  get the local Slater exchange
C
       DA13 = DA**Third
       DB13 = DB**Third
       DA43 = DA*DA13
       DB43 = DB*DB13
       EXC = XS*EFS*(DA43+DB43)/(DA+DB)
       PotA = XS*VFS*DA13
       PotB = XS*VFS*DB13
       If(dft.EQ.1) GO TO 10
C
C  VWN local Correlation
C
       IF(XVWN.NE.Zero) THEN
         CALL VWN(DA,DB,.false.,ECVWN,VAVWN,VBVWN)
         EXC = EXC + ECVWN*XVWN
         PotA = PotA + VAVWN*XVWN
         PotB = PotB + VBVWN*XVWN
         If(dft.EQ.2) GO TO 10
       ELSE IF(XVWN5.NE.Zero) THEN
         CALL VWN5(DA,DB,.false.,ECVWN,VAVWN,VBVWN)
         EXC = EXC + ECVWN*XVWN5
         PotA = PotA + VAVWN*XVWN5
         PotB = PotB + VBVWN*XVWN5
         If(dft.EQ.3) GO TO 10
       ENDIF
c
       g2A = DenAX(1,I)**2 + DenAX(2,I)**2 + DenAX(3,I)**2
       g1A = SQRT(g2A)
       g2B = DenBX(1,I)**2 + DenBX(2,I)**2 + DenBX(3,I)**2
       g1B = SQRT(g2B)
c
       VPotA = Zero
       VPotB = Zero
       VPotAB = Zero
C
C  Becke non-local exchange
C
       IF(XB88.NE.Zero) THEN
         CALL Becke88U(DA,     DA13,   DA43,   g1A,    DB,
     $                 DB13,   DB43,   g1B,    EXB88,  V1A,
     $                 V2A,    V1B,    V2B)
         EXC = EXC + EXB88*XB88
         PotA = PotA + V1A*XB88
         PotB = PotB + V1B*XB88
         VPotA = VPotA + V2A*XB88
         VPotB = VPotB + V2B*XB88
       ENDIF
C
C  Handy/Cohen optimized non-local exchange
C
       IF(XOptX.NE.Zero) THEN
         CALL OPTXU(DA,     DA13,   DA43,   g1A,    DB,
     $              DB13,   DB43,   g1B,    EXHC,   V1A,
     $              V2A,    V1B,    V2B)
         EXC = EXC + EXHC*XOptX
         PotA = PotA + V1A*XOptX
         PotB = PotB + V1B*XOptX
         VPotA = VPotA + V2A*XOptX
         VPotB = VPotB + V2B*XOptX
       ENDIF
c
c  Perdew 86 non-local correlation
c
       IF(XP86L.NE.Zero.OR.XP86NL.NE.Zero) THEN
         g2AB = DenAX(1,I)*DenBX(1,I) + DenAX(2,I)*DenBX(2,I)
     $              + DenAX(3,I)*DenBX(3,I)
         call p86nlcu(da,db,g2a,g2b,g2ab,xp86l,xp86nl,ecp86,
     $                dra,drb,dsaa,dsbb,dsab)
         EXC = EXC + ECP86
         PotA = PotA + dra
         PotB = PotB + drb
         VPotA = VPotA + (dsaa+dsaa)
         VPotB = VPotB + (dsbb+dsbb)
         VPotAB = VPotAB + dsab
       ENDIF
C
C  Perdew-Wang 91 non-local exchange
C
       IF(XPW91X.NE.Zero) THEN
         call pw91xu(da,db,g2a,g2b,expw91,dra,drb,dsaa,dsbb)
         EXC = EXC + expw91*xpw91x
         PotA = PotA + dra*xpw91x
         PotB = PotB + drb*xpw91x
         VPotA = VPotA + (dsaa+dsaa)*xpw91x
         VPotB = VPotB + (dsbb+dsbb)*xpw91x
       ENDIF
C
C  Perdew-Wang non-local correlation
       IF(XPW91L.NE.Zero.OR.XPW91NL.NE.Zero) THEN
c        DenC = DA+DB
c        DenO = DA-DB
c        sigCC = (DenAX(1,I)+DenBX(1,I))**2 +
c    $           (DenAX(2,I)+DenBX(2,I))**2 +
c    $           (DenAX(3,I)+DenBX(3,I))**2
c        sigOO = (DenAX(1,I)-DenBX(1,I))**2 +
c    $           (DenAX(2,I)-DenBX(2,I))**2 +
c    $           (DenAX(3,I)-DenBX(3,I))**2
c        sigCO = (DenAX(1,I)+DenBX(1,I))*(DenAX(1,I)-DenBX(1,I))
c    $         + (DenAX(2,I)+DenBX(2,I))*(DenAX(2,I)-DenBX(2,I))
c    $         + (DenAX(3,I)+DenBX(3,I))*(DenAX(3,I)-DenBX(3,I))
c        CALL dftacg_pw91c(.true., .true., DenC,  DenO,  sigCC,
c    $                     sigCO,  sigOO,  EPW91, VPWC,  VPWO,
c    $                     VCC,    VCO,    VOO)
c        EPW91 = EPW91/DenC
c        VC1A = (VPWC+VPWO)
c        VC1B = (VPWC-VPWO)
c        VC2A = 2.0d0*(VCC+VOO+VCO+VCO)
c        VC2B = 2.0d0*(VCC+VOO-VCO-VCO)
c        VCAB = 2.0d0*(VCC-VOO)
c        EXC = EXC + EPW91*XPW91
c        PotA = PotA + VC1A*XPW91
c        PotB = PotB + VC1B*XPW91
c        VPotA = VPotA + VC2A*XPW91
c        VPotB = VPotB + VC2B*XPW91
c        VPotAB = VPotAB + VCAB*XPW91
         g2AB = DenAX(1,I)*DenBX(1,I) + DenAX(2,I)*DenBX(2,I)
     $              + DenAX(3,I)*DenBX(3,I)
         call pw91nlcu(da,db,g2a,g2b,g2ab,XPW91L,XPW91NL,ecpw91,
     $                 dra,drb,dsaa,dsbb,dsab)
         EXC = EXC + ECPW91
         PotA = PotA + dra
         PotB = PotB + drb
         VPotA = VPotA + (dsaa+dsaa)
         VPotB = VPotB + (dsbb+dsbb)
         VPotAB = VPotAB + dsab
       ENDIF
C
C  Lee, Yang & Parr non-local correlation
C
       IF(XLYP.NE.Zero) THEN
         g2AB = DenAX(1,I)*DenBX(1,I) + DenAX(2,I)*DenBX(2,I)
     $              + DenAX(3,I)*DenBX(3,I)
         CALL ULYP(DA,     DB,     g2A,    g2B,    g2AB,
     $             ECLYP,  VC1A,   VC2A,   VC1B,   VC2B,
     $             VCAB)
         EXC = EXC + ECLYP*XLYP
         PotA = PotA + VC1A*XLYP
         PotB = PotB + VC1B*XLYP
         VPotA = VPotA + VC2A*XLYP
         VPotB = VPotB + VC2B*XLYP
         VPotAB = VPotAB + VCAB*XLYP
       ENDIF
C
C  Perdew, Burke and Ernzerhof non-local exchange
C
       IF(XPBEX.NE.Zero) THEN
         call pbexu(da,db,g2a,g2b,excpbe,dra,drb,dsaa,dsbb)
         EXC = EXC + excpbe*xpbex
         PotA = PotA + dra*xpbex
         PotB = PotB + drb*xpbex
         VPotA = VPotA + (dsaa+dsaa)*xpbex
         VPotB = VPotB + (dsbb+dsbb)*xpbex
       ENDIF
C
C  Perdew, Burke and Ernzerhof non-local correlation
C
       IF(XPBEC.NE.Zero) THEN
         g2AB = DenAX(1,I)*DenBX(1,I) + DenAX(2,I)*DenBX(2,I)
     $              + DenAX(3,I)*DenBX(3,I)
         call pbecu(da,db,g2a,g2b,g2ab,ecpbe,dra,drb,dsaa,dsbb,dsab)
         EXC = EXC + ecpbe*xpbec
         PotA = PotA + dra*xpbec
         PotB = PotB + drb*xpbec
         VPotA = VPotA + (dsaa+dsaa)*xpbec
         VPotB = VPotB + (dsbb+dsbb)*xpbec
         VPotAB = VPotAB + dsab*xpbec
       ENDIF
C
cc
      ENDIF
C
      IF(VPotAB.EQ.Zero) THEN
        PotAX(1,I) = VPotA*DenAX(1,I)
        PotAX(2,I) = VPotA*DenAX(2,I)
        PotAX(3,I) = VPotA*DenAX(3,I)
        PotBX(1,I) = VPotB*DenBX(1,I)
        PotBX(2,I) = VPotB*DenBX(2,I)
        PotBX(3,I) = VPotB*DenBX(3,I)
      ELSE
        PotAX(1,I) = VPotA*DenAX(1,I) + VPotAB*DenBX(1,I)
        PotAX(2,I) = VPotA*DenAX(2,I) + VPotAB*DenBX(2,I)
        PotAX(3,I) = VPotA*DenAX(3,I) + VPotAB*DenBX(3,I)
        PotBX(1,I) = VPotB*DenBX(1,I) + VPotAB*DenAX(1,I)
        PotBX(2,I) = VPotB*DenBX(2,I) + VPotAB*DenAX(2,I)
        PotBX(3,I) = VPotB*DenBX(3,I) + VPOtAB*DenAX(3,I)
      ENDIF
C
 10   CONTINUE
C
C  Accumulate into total Exchange-Correlation energy
C  and total number of electrons
C
      DVal = DenA(I)+DenB(I)
      ETot = ETot + EXC*DVal*WGHT(I)
      EL   = EL + DVal*WGHT(I)
C
C  store full potential and exchange-correlation energy
C
      FPotA(I) = PotA
      FPotB(I) = PotB
      EVec(I) = EXC*DVal    ! no need for this without weight derivatives
cc
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE VWN5(DenA,DenB,rhf,EC,PotA,PotB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Vosko, Wilk & Nusair (VWN)  Local correlation
C  modified from version used in DMOL
C  (This is the fit to the electron gas results of Ceperley &
C   Alder, and is the version recommended in the paper)
C
C  reference
C  S.H.Vosko, L.Wilk and M.Nusair,  Can J. Phys. 58 (1980) 1200
C
C  ARGUMENTS
C
C  DenA    -  alpha (or closed shell) density over grid
C  DenB    -  beta density over grid
C  rhf     -  logical flag indicating open or closed-shell
C
C  on exit
C
C  EC      -  contribution to local correlation energy
C  PotA    -  alpha (or closed shell) potential
C  PotB    -  beta potential
C
C
      REAL*8 A(3),B(3),C(3),X0(3),Q(3)
      LOGICAL rhf
C
      parameter (Half=0.5d0,One=1.0d0,Two=2.0d0,Four=4.0d0)
      parameter (Third=One/3.0d0,Sixth=One/6.0d0)
      parameter (PI=3.14159 26535 89793d0)
c
      PARAMETER (aa=2.519842099789746d0,ab=0.2599210498948732d0,
     $           f2=1.709920934161365d0)
C
      DATA A /-0.0337737d0,  0.0621814d0, 0.0310907d0/
      DATA B / 1.13107d0,    3.72744d0,   7.06042d0/
      DATA C /13.00450d0,   12.93520d0,  18.05780d0/
      DATA X0/-0.0047584d0, -0.10498d0,  -0.32500d0/
      DATA Q/7.12310891781812d0,6.15199081975908d0,4.73092690956011d0/
cc
cc
c  function statements
      f(z) = ((One+z)*(One+z)**Third+(One-z)*(One-z)**Third-Two)
     $       /(aa-Two)
      rs(z) = (0.75d0/(PI*z))**Third
      cx(x,m) = x**2 + b(m)*x + c(m)
      fx(x,m) = atan(q(m)/(Two*x+b(m)))
      e1(x,m) = log(x**2/cx(x,m)) + Two*b(m)*fx(x,m)/q(m)
      e2(x,m) = (x-x0(m))**2/cx(x,m)
      e3(m) = log(e2(x,m)) + (Two*b(m)+Four*x0(m))*fx(x,m)/q(m)
      epsilc(m) = a(m)*Half*(e1(x,m)-b(m)*x0(m)*e3(m)/cx(x0(m),m))
      g1(x,m) = Two/x-(One-b(m)*x0(m)/cx(x0(m),m))*(two*x+b(m))/cx(x,m)
      g2(x,m) = Two*b(m)*x0(m)/(x-x0(m))/cx(x0(m),m)
      g3(x,m) = Four*b(m)*(One-(b(m)+Two*x0(m))*x0(m)/cx(x0(m),m))
     $          /(q(m)**2+(Two*x+b(m))**2)
      dgdx(m) = a(m)*(g1(x,m)-g2(x,m)-g3(x,m))*Half
C
C
      IF(rhf) THEN
C
C  ** Closed Shell **
C
        rs1 = rs(DenA)
        x = sqrt(rs1)
cc  energy
        EC = epsilc(2)
cc  potential
        PotA = EC - x*sixth*dgdx(2)
C
      ELSE
C
C  ** Unrestricted Open Shell **
C     This is coded in terms of the total and spin (difference) densities
C
        DVal = DenA + DenB
        DDif = DenA - DenB
c
        rs1 = rs(DVal)
        x = sqrt(rs1)
        ecp = epsilc(2)
        decp = dgdx(2)
c
        ac = epsilc(1)
        ecf = epsilc(3)
        dac = dgdx(1)
        decf = dgdx(3)
        beta = f2*(ecf-ecp)/ac - One
        z = min( max ((DDif/DVal),-One),One)
        d1 = z**4
        d2 = f(z)
        d4 = Four*ac/f2
        gu1 = x*Sixth*((One-d1*d2)*decp + d1*d2*decf +
     $                 (One-d1)*d2*dac/f2)
        gu2 = d4*(z**3*d2*beta + (One+beta*d1)*Sixth/ab *
     $       ((One+z)**Third - (One-z)**Third))
        d5 = (ecf-ecp-ac/f2)*d1
cc  energy
        EC = ecp + d2*(ac/f2+d5)
cc  potential
        PotA = EC - gu1 + gu2*(One-z)
        PotB = EC - gu1 - gu2*(One+z)
C
      ENDIF
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE VWN(DenA,DenB,rhf,EC,PotA,PotB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Vosko, Wilk & Nusair (VWN)  Local correlation
C  (This is the RPA fit)
C
C  reference
C  S.H.Vosko, L.Wilk and M.Nusair,  Can J. Phys. 58 (1980) 1200
C
C  ARGUMENTS
C
C  DenA    -  alpha (or closed shell) density over grid
C  DenB    -  beta density over grid
C  rhf        logical flag indicating open or closed shell
C
C  on exit
C
C  EC      -  contribution to local correlation energy
C  PotA    -  alpha (or closed shell) potential
C  PotB    -  beta potential
C
C
      DIMENSION AD(2),bD(2),cD(2),x0D(2),A2D(2),QD(2),bXx0D(2),
     $          dXx0D(2)
      LOGICAL rhf
      SAVE AD,bD,cD,x0D,A2D,QD,bXx0D,dXx0D
C
      parameter (Zero=0.0d0,One=1.0d0,Two=2.0d0)
      parameter (Third=One/3.0d0)
      parameter (PI=3.14159 26535 89793d0)
      parameter (FP3=3.0d0/(4.0d0*PI))
c
      DATA AD /0.0621814d0, 0.0310907d0/
      DATA bD /13.0720d0, 20.1231d0/
      DATA cD /42.7198d0, 101.578d0/
      DATA x0D /-0.409286d0, -0.743294d0/
c
      DATA A2D /0.0310907d0, 0.01554535d0/
      DATA QD  /0.0448998886415767975d0, 1.17168527770897146d0/
      DATA bXx0D /-0.142530524167983924d0, -0.171582499414507567d0/
      DATA dXx0D /12.253428d0, 18.636512d0/
C
C .................................................................
C  NOTE  (as worked out by PP)
C  This routine is coded in terms of the Wigner-Seitz radius
C   rs = (3/(4*Pi*density))
C  We have  (R = density)
C   Exc = Integral{E(R)*R}
C   Vxc = Integral{(dE/dR)*R + E(R)}
C  Since this routine forms the "potential" in terms of rs
C  the final potential is given by
C   Vxc = Integral{(dE/drs)*(drs/dR)*R + E(rs)}
C .................................................................
C
C 1. Closed-Shell/paramagnetic term
C
      If(rhf) Then
        rs = (FP3/DenA)**Third
      Else
        DVal = DenA+DenB
        rs = (FP3/DVal)**Third
      EndIf
      x = SQRT(rs)
c .................................................................
      A2 = A2D(1)
      b = bD(1)
      c = cD(1)
      x0 = x0D(1)
      bXx0 = bXx0D(1)
      dXx0 = dXx0D(1)
      Q = QD(1)
c ..................................................................
c
      U = rs + b*x + c
      S1 = LOG(rs/U)
      xdif = x - x0
      S2 = Log(xdif*xdif/U)
      V = Two*x + b
      T = (Two/Q)*ATAN(Q/V)
c
      dUdr = One + b/(Two*x)
      dLogU = dUdr/U
      dS1dr = One/rs - dLogU
      dS2dr = One/(x*xdif) - dLogU
      QV = Q*Q + V*V
      dTdr = -Two/(x*QV)
c
      ECP = A2*(S1 + b*T - bXx0*(S2 + dXx0*T))
      PotA = A2*(dS1dr + b*dTdr - bXx0*(dS2dr + dXx0*dTdr))
c
      If(rhf) Then
       EC = ECP
c -- now multiply by (drs/dR)
       PotA = -Third*PotA*rs + ECP
       RETURN
      EndIf
C
C
C 2. Open-Shell/ferromagnetic term
C
c .................................................................
      A2 = A2D(2)
      b = bD(2)
      c = cD(2)
      x0 = x0D(2)
      bXx0 = bXx0D(2)
      dXx0 = dXx0D(2)
      Q = QD(2)
c ..................................................................
c
      U = rs + b*x + c
      S1 = LOG(rs/U)
      xdif = x - x0
      S2 = Log(xdif*xdif/U)
      V = Two*x + b
      T = (Two/Q)*ATAN(Q/V)
c
      dUdr = One + b/(Two*x)
      dLogU = dUdr/U
      dS1dr = One/rs - dLogU
      dS2dr = One/(x*xdif) - dLogU
      QV = Q*Q + V*V
      dTdr = -Two/(x*QV)
c
      ECF = A2*(S1 + b*T - bXx0*(S2 + dXx0*T))
      PotB = A2*(dS1dr + b*dTdr - bXx0*(dS2dr + dXx0*dTdr))
c .................................................................
c
      TwoCon = Two*(Two**Third - One)
      TwoCI = One/TwoCon
      TwoCID = TwoCI + Third*TwoCI
c
      Zeta = (DenA-DenB)/DVal
      fZet = -Two
      OMZ = One - Zeta
      OPZ = One + Zeta
      OMZ3 = OMZ**Third
      fZet = fZet + OMZ*OMZ3
      OPZ3 = OPZ**Third
      fZet = fZet + OPZ*OPZ3
      fZet = fZet*TwoCI
      EC = ECP + fZet*(ECF - ECP)
c
      dZdRA =  OMZ/DVal
      dZdRB = -OPZ/DVal
      dfZdz = TwoCID*(-OMZ3 + OPZ3)
      dXdRS = PotA + fZet*(PotB - PotA)
      dXTemp = -Third*rs*dXdRS/DVal
      dXdz = dfZdz*(ECF - ECP)
      dXdRA = dXTemp + dXdz*dZdRA
      dXdRB = dXTemp + dXdz*dZdRB
c
      PotA = DVal*dXdRA + EC
      PotB = DVal*dXdRB + EC
C
      RETURN
      END
c =======================================================================
c
      SUBROUTINE Becke88(Dens,   Dens13, Dens43, g1,     EXB88,
     $                   Pot,    PotX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Becke 1988 nonlocal exchange
C  gradient-corrected exchange term
C  implemented as per Johnson, Gill & Pople   J.Chem.Phys. 98 (1993) 5612
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  Dens    -  Half closed-shell density at current grid point
C  Dens13  -  Dens**(1/3)
C  Dens43  -  Dens**(4/3)
C  g1      -  absolute value of density derivative
C             (WARNING: evaluated using half closed-shell density)
C
C  on exit
C
C  EXB88   -  exchange energy contribution at current grid point
C  Pot     -  potential
C  PotX    -  potential derivative
C
      parameter(Zero=0.0d0,One=1.0d0,Two=2.0d0)
      parameter(FThird=4.0d0/3.0d0,thrsh=1.0d-12)
c
      PARAMETER (b=0.0042d0,b6=0.0252d0)
      PARAMETER (b6b=b6*b)
C
C
      XA = g1/Dens43
      XA2 = XA*XA
c
      SS = SQRT(XA2+One)
      SXA = LOG(XA+SS)
      fac = (One + b6*XA*SXA)
      fXA = -b*XA2/fac
cc energy
      EXB88 = Dens13*fXA
cc potential
      gXA = ( b6b*XA2*( XA/SS - SXA) - Two*b*XA ) / (fac*fac)
      Pot = FThird*Dens13*(fXA - XA*gXA)
      If(XA.LT.thrsh) Then
        PotX = Zero
      Else
        PotX = gXA/(Dens43*XA)     ! factor of 1/2 already included
      EndIf                        ! in density derivatives
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE Becke88U(DA,    DA13,   DA43,   g1A,    DB,
     $                    DB13,  DB43,   g1B,    EXB88,  PotA,
     $                    PotAX, PotB,   PotBX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Becke 1988 nonlocal exchange
C  gradient-corrected exchange term
C  implemented as per Johnson, Gill & Pople   J.Chem.Phys. 98 (1993) 5612
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  DA      -  alpha density at current grid point
C  DA13    -  DA**(1/3)
C  DA43    -  DA**(4/3)
C  g1A     -  absolute value of alpha density derivative
C  DB      -  beta density at current grid point
C  DB13    -  DB**(1/3)
C  DB43    -  DB**(4/3)
C  g1B     -  absolute value of beta density derivative
C
C  on exit
C
C  EXB88   -  exchange energy contribution at current grid point
C  PotA    -  alpha potential
C  PotAX   -  alpha potential derivative
C  PotB    -  beta potential
C  PotBX   -  beta potential derivative
C
      parameter(Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0)
      parameter(FThird=4.0d0/3.0d0,thrsh=1.0d-12)
c
      PARAMETER (b=0.0042d0,b6=0.0252d0)
      PARAMETER (b6b=b6*b)
C
C
C  alpha spin
      IF(DA.GT.thrsh) THEN
       XA = g1A/DA43
       XA2 = XA*XA
c
       SS = SQRT(XA2+One)
       SXA = LOG(XA+SS)
       fac = (One + b6*XA*SXA)
       fXA = -b*XA2/fac
C  potential
       gXA = ( b6b*XA2*( XA/SS - SXA) - Two*b*XA ) / (fac*fac)
       PotA = FThird*DA13*(fXA - XA*gXA)
       If(XA.LT.thrsh) Then
         PotAX = Zero
       Else
         PotAX = gXA/(DA43*XA)     ! factor of 1/2 already included
       EndIf                       ! in density derivatives
      ELSE
       fXA = Zero
       PotA = Zero
       PotAX = Zero
      ENDIF
C
C  beta spin
      IF(DB.GT.thrsh) THEN
       XA = g1B/DB43
       XA2 = XA*XA
c
       SS = SQRT(XA2+One)
       SXA = LOG(XA+SS)
       fac = (One + b6*XA*SXA)
       fXB = -b*XA2/fac
C  potential
       gXA = ( b6b*XA2*( XA/SS - SXA) - Two*b*XA ) / (fac*fac)
       PotB = FThird*DB13*(fXB - XA*gXA)
       If(XA.LT.thrsh) Then
         PotBX = Zero
       Else
         PotBX = gXA/(DB43*XA)     ! factor of 1/2 already included
       EndIf                       ! in density derivatives
      ELSE
       fXB = Zero
       PotB = Zero
       PotBX = Zero
      ENDIF
C
C  energy
      EXB88 = (DA43*fXA + DB43*fXB)/(DA+DB)
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE LYP(Dens,   Dens13, Dens43, g2,    ECLYP,
     $               Pot,    PotX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Lee-Yang-Parr non-local correlation term, as
C  reformulated by A.Savin
C  ** CLOSED SHELL **
C
C  reference
C  C.Lee, W.Yang and R.G.Parr,  Phys.Rev. B37 (1988) 785
C  B.Miehlich, A.Savin, H.Stoll and H.Preuss, Chem.Phys.Lett.
C
C  ARGUMENTS
C
C  Dens    -  half of the closed-shell density at current grid point
C  Dens13  -  Dens**(1/3)
C  Dens43  -  Dens**(4/3)
C  g2      -  absolute value of the density gradient squared
C
C  on exit
C
C  ECLYP   -  exchange energy contribution at current grid point
C  V1A     -  potential
C  V2A     -  potential gradient
C
      parameter(half=0.5d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0,
     1 five=5.0d0,six=6.0d0,seven=7.0d0,eight=8.0d0,ten=10.0d0)
      parameter(third=one/three,a23=two*third,a14p9=14.0d0/9.0d0,
     1 a14p27=14.0d0/27.0d0,a43=four*third,a53=five*third,
     2 a73=seven*third,egy9=one/9.0d0,egy24=one/24.0d0,
     3 half13=0.7937005259841d0)
      parameter(alyp=0.04918d0,blyp=0.132d0,clyp=0.2533d0,dlyp=0.349d0)
      parameter(ablyp=alyp*blyp,clyp1=clyp*half13,dlyp1=dlyp*half13)
c    2(-11/3) and (2*3/10)*(3*pi**2)**(2/3)
      parameter(twom113=0.07874506561842d0,const=5.74246800037638d0)
      parameter(alyp2=two*alyp)
c
      Densm13=One/Dens13
      Densm43=One/Dens43
      Densm53=Dens13/(Dens*Dens)
      denom=One/(Densm13*DLYP1+One)
c
c    first term
      XLYP1 = -ALYP2*DENOM*Dens
      DERLYP1 = -ALYP*(DENOM*DENOM*Densm13*DLYP1*Third+DENOM)
c
c   Second term
      excterm=EXP(-CLYP1*Densm13)*denom
c   exct is just the exponential part of excterm - replace ir w/ excterm
      XLYP2 = -ABLYP*CONST*EXCTERM*Dens
c  take half of the MACSYMA output because of der. w.r.t. densalpha
      derlyp2=-ablyp*const*excterm*(Densm13*Third*(denom*dlyp1+clyp1)
     1   +One)*Half
c
c  Third term
c  newest - directly out of the machine
      GLYP3 = ABLYP*(A14P9*CLYP1*Densm13-A14P9*DENOM+A23+A14P9)
     1   *Densm53*EXCTERM*TWOM113
      XLYP3=GLYP3*G2
c
      FACTOR = EXCTERM*TWOM113*g2
      Dens23=Dens13*Dens13
      DERLYP3 = (-ABLYP*DENOM**2*Densm13*
     1   (((-Seven*CLYP1**2/Dens23-Ten*CLYP1*Densm13-Three)*DLYP1**2
     2   +(-14.0d0*CLYP1**2*Densm13-Six*CLYP1+Four*Dens13)*DLYP1
     3    -Seven*CLYP1**2+Four*CLYP1*Dens13)
     x    +((35.0d0*CLYP1*Densm13
     4    +50.0d0)*DLYP1**2+(70.0d0*CLYP1+65.0d0*Dens13)*DLYP1
     5   +35.0d0*CLYP1*Dens13+15.0d0*Dens23)
     6    -Seven*DLYP1**2)*FACTOR)/(27.0d0*Dens**3)
c
      ECLYP = (xlyp1+xlyp2+xlyp3)/(Dens+Dens)
      Pot = derlyp1+derlyp2+derlyp3
      PotX = glyp3                   ! factor of 1/2 already included
c                                    ! in density derivatives
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE ULYP(DA,     DB,     g2A,    g2B,   g2AB,
     $                ECLYP,  PotA,   PotAX,  PotB,  PotBX,
     $                PotABX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Lee-Yang-Parr non-local correlation term
C  ** OPEN SHELL **
C
C  reference
C  C.Lee, W.Yang and R.G.Parr,  Phys.Rev. B37 (1988) 785
C  B.Miehlich, A.Savin, H.Stoll and H.Preuss, Chem.Phys.Lett.
C
C  ARGUMENTS
C
C  DA      -  alpha density at current grid point
C  DB      -  beta density at current grid point
C  g2A     -  absolute value of alpha density derivative squared
C  g2B     -  absolute value of beta density derivative squared
C  g2AB    -  mixed alpha-beta density derivative square
C
C  on exit
C
C  ECLYP   -  exchange energy contribution at current grid point
C  PotA    -  alpha potential
C  PotAX   -  alpha potential derivative
C  PotB    -  beta potential
C  PotBX   -  beta potential derivative
C  PotABX  -  mixed alpha-beta potential derivative
C
C
      Parameter (a=0.04918d0,b=0.132d0,c=0.2533d0,d=0.349d0)
      Parameter (o83=8.0d0/3.0d0,o23=2.0d0/3.0d0,o14=0.25d0,
     $           o19=1.0d0/9.0d0,o113=11.0d0/3.0d0)
      Parameter (Zero=0.0d0,One=1.0d0,Two=2.0d0,Three=3.0d0,
     $           Four=4.0d0,Seven=7.0d0,Eleven=11.0d0,
     $           TwtyFour=24.0d0,FtySeven=47.0d0)
      Parameter (Third=One/Three,FThird=Four/Three)
      Parameter (PI=3.14159 26535 89793d0)
c
cc    fac2 = (Two**o113)*0.3d0*(Three*PI*PI)**o23
      Parameter (fac2=36.4623989787647673d0)
C
      DATA epsma/1.0d-20/, thrsh/1.0d-12/
c
      DVal = DA + DB
c
      ro13 = DVal**(-Third)
      ro43 = ro13**Four
      dg = One/(One + d*ro13)
      wr = exp(-c*ro13)*dg*DVal**(-o113)
c
      IF(wr.gt.epsma) THEN
cc
        delr = (c + d*dg)*ro13
        ror = One/DVal
        abwr = a*b*wr
        abwud = abwr*DA*DB
        rudo = DA*DB*ror
        ru83 = DA**o83
        rd83 = DB**o83
c
        fp1 = -Four*a*dg*rudo
        fp2 = -fac2*abwud*(ru83+rd83)
        fpg = -abwud*o19
        fpi = One - Three*delr
        fpj = (delr-Eleven)*ror
        f47 = FtySeven - Seven*delr
        dlguu = fpg*(fpi - fpj*DA) + abwr*DB*DB
        dlgud = fpg*f47 + FThird*abwr*DVal*DVal
        dlgdd = fpg*(fpi - fpj*DB) + abwr*DA*DA

        vc1 = fp1 + fp2 + dlguu*g2A + dlgud*g2AB + dlgdd*g2B
c energy
        ECLYP = vc1*ror
c potential
        wprim = -Third*ro43*wr*(Eleven/ro13 - c - d*dg)
        delprim = Third*(d*d*ro43*ro13*dg*dg - delr*ror)
        wprw = wprim/wr
c
        fpju = fpj*ror*DB
        fpjd = fpj*ror*DA
        fpkd =  - abwr*o19*(One - Three*delr - fpj*DA)
        fpku =  - abwr*o19*(One - Three*delr - fpj*DB)
        fplu =  - fpg*(Three+DA*ror)*delprim
        fpld =  - fpg*(Three+DB*ror)*delprim
        fp9 = o19*DA*DB
        fpn = wprw*dlgud
     $         + abwr*o19*(Seven*DA*DB*delprim + TwtyFour*DVal)
c
        dlu1 = fp1*(Third*d*ro43*dg - ror) + fp2*wprw
c
c potential alpha (u)
c
        IF(DA.GT.thrsh) THEN
          dluguu = wprw*dlguu + fpkd*DB + fplu - fpg*fpju
          dlugud = fpn - abwr*o19*DB*f47
          dlugdd = wprw*dlgdd + fpku*DB +fpld + fpg*fpju
     $               + abwr*Two*DA
c
          dlu2u = fp1/DA
     $            - fac2*abwr*DB*(o113*ru83+rd83)
     $            + dluguu*g2A + dlugud*g2AB + dlugdd*g2B
          dlu = dlu1 + dlu2u
c
          PotA  = dlu
          PotAX = dlguu + dlguu
        ELSE
          PotA  = Zero
          PotAX = Zero
        ENDIF
c
c potential beta(d)
c
        IF(DB.GT.thrsh) THEN
          dldguu = wprw*dlguu + fpkd*DA + fplu + fpg*fpjd
     $             + abwr*Two*DB
          dldgud = fpn - abwr*o19*DA*f47
          dldgdd = wprw*dlgdd + fpku*DA + fpld - fpg*fpjd

          dlu2d = fp1/DB
     $            - fac2*abwr*DA*(o113*rd83+ru83)
     $            + dldguu*g2A + dldgud*g2AB + dldgdd*g2B
          dld = dlu1 + dlu2d
c
          PotB  = dld
          PotBX = dlgdd + dlgdd
        ELSE
          PotB  = Zero
          PotBX = Zero
        ENDIF
c
c  alpha-beta term
        PotABX = dlgud
cc
      ELSE
cc
        ECLYP = Zero
        PotA = Zero
        PotB = Zero
        PotAX = Zero
        PotBX = Zero
        PotABX = Zero
      ENDIF
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE OPTX(Dens,   Dens13, Dens43, g1,    EXHC,
     $                Pot,    PotX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Handy/Cohen 2001 nonlocal exchange
C  gradient-corrected exchange term
C  N. C. Handy & A. J. Cohen, Mol.Phys. 99 (2001) 403
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  Dens    -  Half closed-shell density at current grid point
C  Dens13  -  Dens**(1/3)
C  Dens43  -  Dens**(4/3)
C  g1      -  absolute value of density derivative
C             (WARNING: evaluated using half closed-shell density)
C
C  on exit
C
C  EXHC    -  exchange energy contribution at current grid point
C  Pot     -  potential
C  PotX    -  potential derivative
C
      parameter(Zero=0.0d0,One=1.0d0,Four=4.0d0)
      parameter(FThird=4.0d0/3.0d0,thrsh=1.0d-12)
c
      PARAMETER (A=-1.43169d0,b=0.006d0)
C
C
      XA = g1/Dens43
      XA2 = XA**2
c
      fac = (One + b*XA2)
      rnum = b*XA2/fac
      fXA = A*rnum**2
cc energy
      EXHC = Dens13*fXA
cc potential
      gXA = (Four*A*b*XA*rnum)/(fac**2)
      Pot = FThird*Dens13*(fXA - XA*gXA)
      If(XA.LT.thrsh) Then
        PotX = Zero
      Else
        PotX = gXA/(Dens43*XA)     ! factor of 1/2 already included
      EndIf                        ! in density derivatives
C
      RETURN
      END
c =====================================================================
c
      SUBROUTINE OPTXU(DA,    DA13,   DA43,   g1A,    DB,
     $                 DB13,  DB43,   g1B,    EXHC,   PotA,
     $                 PotAX, PotB,   PotBX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Handy/Cohen 2001 nonlocal exchange
C  gradient-corrected exchange term
C  N. C. Handy & A. J. Cohen, Mol.Phys. 99 (2001) 403
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  DA      -  alpha density at current grid point
C  DA13    -  DA**(1/3)
C  DA43    -  DA**(4/3)
C  g1A     -  absolute value of alpha density derivative
C  DB      -  beta density at current grid point
C  DB13    -  DB**(1/3)
C  DB43    -  DB**(4/3)
C  g1B     -  absolute value of beta density derivative
C
C  on exit
C
C  EXHC    -  exchange energy contribution at current grid point
C  PotA    -  alpha potential
C  PotAX   -  alpha potential derivative
C  PotB    -  beta potential
C  PotBX   -  beta potential derivative
C
      parameter(Zero=0.0d0,One=1.0d0,Four=4.0d0)
      parameter(FThird=4.0d0/3.0d0,thrsh=1.0d-12)
c
      PARAMETER (A=-1.43169d0,b=0.006d0)
C
C
C  alpha spin
      IF(DA.GT.thrsh) THEN
       XA = g1A/DA43
       XA2 = XA**2
c
       fac = (One + b*XA2)
       rnum = b*XA2/fac
       fXA = A*rnum**2
C  potential
       gXA = (Four*A*b*XA*rnum)/(fac**2)
       PotA = FThird*DA13*(fXA - XA*gXA)
       If(XA.LT.thrsh) Then
         PotAX = Zero
       Else
         PotAX = gXA/(DA43*XA)     ! factor of 1/2 already included
       EndIf                       ! in density derivatives
      ELSE
       fXA = Zero
       PotA = Zero
       PotAX = Zero
      ENDIF
C
C  beta spin
      IF(DB.GT.thrsh) THEN
       XA = g1B/DB43
       XA2 = XA**2
c
       fac = (One + b*XA2)
       rnum = b*XA2/fac
       fXB = A*rnum**2
C  potential
       gXA = (Four*A*b*XA*rnum)/(fac**2)
       PotB = FThird*DB13*(fXB - XA*gXA)
       If(XA.LT.thrsh) Then
         PotBX = Zero
       Else
         PotBX = gXA/(DB43*XA)     ! factor of 1/2 already included
       EndIf                       ! in density derivatives
      ELSE
       fXB = Zero
       PotB = Zero
       PotBX = Zero
      ENDIF
C
C  energy
      EXHC = (DA43*fXA + DB43*fXB)/(DA+DB)
C
      RETURN
      END
c=======================================================================
      subroutine pw91lcr(r,ec,dra)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/21/2003) generated with maxima, file pw91lc.mc
c...
c...  Perdew and Wang local correlation functional
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...              closed shell version
c...
c...  Formulation taken from the Molpro manual (PW92C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first
c...  derivative.
c...
c...  Input parameters:
c...
c...  r     closed shell total density
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...
      real*8 r,ec,dra
      real*8 rm13
      real*8 x,x2,x12,x32,polx1
      real*8 eps11,eps1p1,eps1,epsp1
      real*8 mthird,th
      parameter (mthird=-1.0d0/3.0d0,th=1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      rm13 = r**mthird
c...
c...  the auxiliary variable x and some of its powers.
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x12 = SQRT(x)
      x32 = x**th
c...
c...  the polynomial of x appearing in eps1.
c...
      polx1 = 1.6382D0*x32+4.9294D-1*x2+7.5957D0*x12+3.5876D0*x
c...
c...  the logarithm part of eps1.
c...
      eps11 = LOG(1.608197949869254D1/polx1+1.0D0)
c...
c...  and finally eps1.
c...
      eps1 = -6.21814D-2*eps11*(2.137D-1*x+1.0D0)
c...
c...  Functional value (divided by the density).
c...
      ec = eps1
c...
c...  First derivative (potential).
c...
      eps1p1 = -5.0D-1*(x*(1.97176D0*x12+4.9146D0)+7.1752D0*x12+7.5957D0
     1   )/(polx1*(6.21814D-2*polx1+1.0D0)*x12)
      epsp1 = -6.21814D-2*eps1p1*(2.137D-1*x+1.0D0)-1.328816518D-2*eps11
      dra = eps1-2.067834969664667D-1*epsp1*rm13
c...
      return
      end
c=======================================================================
      subroutine pw91nlcr(ra,saa,scl,scnl,ec,dra,dsaa)
      implicit none
c...
c...  MM (05/30/2003) generated with maxima, file pw91nlcr.mc
c...
c...  Perdew and Wang non local correlation functional
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...              closed shell version
c...
c...  Formulation taken from the Molpro manual (PW91C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first
c...  derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  saa   alpha density gradient invariant
c...  scl   scaling factor for local part
c...  scnl  scaling factor for non local part
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    functional derivative with respect to alpha gradient
c...
      real*8 ra,saa,scl,scnl,ec,dra,dsaa
      real*8 pw91lc,dpw91lcdra
      real*8 j,jpd,jpx,jpr,djdra,djdsaa,expj,polj
      real*8 l,lpe,lpd,dldra,dldsaa
      real*8 r,r13,rm13,rm43,rm76,rm136
      real*8 s,s12
      real*8 x,x2,x3
      real*8 d,d2,d4
      real*8 e,epra
      real*8 numt,dent,dentp,t,tp
      real*8 expa,a,ape,loga,logape,logapd,numlog,denlog
      real*8 zero,one,third,ssm,ttsm
      parameter (zero=0.0d0,one=1.0d0,third=1.0d0/3.0d0)
      parameter (ssm=-7.0d0/6.0d0,ttsm=-13.0d0/6.0d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      r = 2.0D0*ra
      call pw91lcr(r,pw91lc,dpw91lcdra)
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      s = 4.0D0*saa
      if(s.lt.epsi.or.scnl.eq.zero)then
      j = ZERO
      djdra = ZERO
      djdsaa = ZERO
      l = ZERO
      dldra = ZERO
      dldsaa = ZERO
      else
c...
c...  Non local part.
c...  we compute some powers of r and s
c...
      r13 = r**THIRD
      rm13 = one/r13
      rm43 = rm13/r
      rm76 = r**ssm
      rm136 = r**ttsm
      s12 = SQRT(s)
c...
c...  the auxiliary variables x some of its powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x3 = x*x2
c...
c...  the variable d.
c...
      d = 2.519289703422449D-1*rm76*s12
      d2 = d**2
      d4 = d2**2
c...
c...  and the last auxiliary variable, e.
c...
      e = pw91lc
c...
c...  we start constructing the term j.
c...
c...  function t.
c...
      numt = 7.389D-3*x2+2.3266D1*x+2.568D0
      dent = 7.389D-2*x3+4.72D-1*x2+8.723D0*x+1.0D0
c...
c...  we build j.
c...
      t = 1.0D-3*numt/dent
      expj = EXP(-4.11563120990411453D1*d2*rm13)
      polj = t-1.853268571428571D-3
      j = 1.575592034948315D1*d2*expj*polj
c...
c...  now the term l
c...
c...  function a.
c...
      expa = EXP(-4.04276151175637226D1*e)
      a = 2.697586091519874D0/(expa-1.0D0)
c...
c...  we build l.
c...
      numlog = a*d4+d2
      denlog = a**2*d4+a*d2+1.0D0
      loga = LOG(2.697586091519874D0*numlog/denlog+1.0D0)
      l = 2.473556743557577D-2*loga
c...
c...  First derivative of e
c...
      epra = -1.0D0*(e-1.0D0*dpw91lcdra)/r
c...
c...  First derivatives of j w.r.t. auxiliary variables
c...
      dentp = 2.2167D-1*x2+9.44D-1*x+8.723D0
      tp = 1.0D-3*(dent*(1.4778D-2*x+2.3266D1)-dentp*numt)/dent**2
      jpd = 2.057815604952057D-1*j*(9.719043801529514D0*r13-4.0D2*d2)*rm
     1   13/d
      jpx = 1.575592034948315D1*d2*expj*tp
      jpr = 1.371877069968038D1*d2*j*rm43
c...
c...  First derivatives of l w.r.t. auxiliary variables
c...
      ape = 1.090569722544585D2*expa/(expa-1.0D0)**2
      logape = -1.8D-1*ape*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*denlog
     1   )/(denlog*(1.8D-1*numlog+6.672632268006112D-2*denlog))
      logapd = -3.6D-1*d*(2.0D0*a*d2+1.0D0)*(a*numlog-denlog)/(den
     1   log*(1.8D-1*numlog+6.672632268006112D-2*denlog))
      lpe = 2.473556743557577D-2*logape
      lpd = 2.473556743557577D-2*logapd
c...
c...  Total first derivatives of the non local terms:
c...   1) j
c...   2) l
c...
      djdra = -2.939171320659524D-1*jpd*rm136*s12-2.067834969664667D-1*j
     1   px*rm43+jpr
      djdsaa = 1.259644851711224D-1*jpd*rm76/s12
      dldra = epra*lpe-2.939171320659524D-1*lpd*rm136*s12
      dldsaa = 1.259644851711224D-1*lpd*rm76/s12
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = (l+j)*scnl+pw91lc*scl
c...
c...  First derivatives of the functional (potential).
c...
      dra = ((dldra+djdra)*r+l+j)*scnl+dpw91lcdra*scl
      dsaa = (dldsaa+djdsaa)*r*scnl
c...
      return
      end
c=======================================================================
      subroutine pw91lcu(ra,rb,ec,dra,drb)
      implicit none !do not comment this line, no implicit convention here
c...  MM (05/21/2003) generated with maxima, file pw91lcu.mc
c...
c...  Perdew and Wang local correlation functional
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...              open shell version
c...
c...  Formulation taken from the Molpro manual (PW92C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dra     funct. deriv. w.r.t. beta density
c...
      real*8 ra,rb,ec,dra,drb
      real*8 z,z2,z3,z4
      real*8 zp1,zp113,zp143
      real*8 zm1,zm113,zm143
      real*8 x,x2,x12,x32,polx1,polx2,polx3
      real*8 r,r13,rm13,rm43
      real*8 eps11,eps12,eps13,eps1p1,eps1p2,eps1p3
      real*8 eps1,eps2,eps3,epsp1,epsp2,epsp3
      real*8 om,omp
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
      real*8 third,tthird,etm,th
      parameter (third=0.33333333333333333333d0)
      parameter (tthird=0.66666666666666666667d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      r = rb+ra
      r13 = r**THIRD
      rm13 = one/r13
      rm43 = rm13/r
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.2035049089940001667D-1*rm13
      x2 = x**2
      x12 = SQRT(x)
      x32 = x**th
      z = min(max((ra-rb)/r,-one),one)
      z2 = z**2
      z3 = z*z2
      z4 = z*z3
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp143 = zp1*zp113
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm143 = zm1*zm113
c...
c...  the polynomials of x appearing in eps1, eps2, and eps3.
c...
      polx1 = 1.6382D0*x32+4.9294D-1*x2+7.5957D0*x12+3.5876D0*x
      polx2 = 3.3662D0*x32+6.2517D-1*x2+1.41189D1*x12+6.1977D0*x
      polx3 = 8.8026D-1*x32+4.9671D-1*x2+1.0357D1*x12+3.6231D0*x
c...
c...  the logarithm parts of eps1, eps2, and eps3.
c...  if we are executing this code, it means that r (density) has already
c...  been screened, hence x, (and polx) are significant. Thus, no need to test
c...  the argument of LOG
c...
      eps11 = LOG(1.6081979498692535067D1/polx1+1.0D0)
      eps12 = LOG(3.2163958997385070134D1/polx2+1.0D0)
      eps13 = LOG(2.9608749977793437517D1/polx3+1.0D0)
c...
c...  and finally eps1, eps2, and eps3.
c...
      eps1 = -6.21814D-2*eps11*(2.137D-1*x+1.0D0)
      eps2 = -3.10907D-2*eps12*(2.0548D-1*x+1.0D0)
      eps3 = -3.37738D-2*eps13*(1.1125D-1*x+1.0D0)
c...
c...  the function omega.
c...  the constant is 1/(2^(4/3)-2)
c...
      om = 1.9236610509315363198D0*(zp143+zm143-2.0D0)
c...
c...  Functional value (divided by the density).
c...
      ec = 5.8482233974552040708D-1*(-1.709921D0*eps1*(om*z4-1.0D0)+eps3
     1   *om*(z4-1.0D0)+1.709921D0*eps2*om*z4)
c...
c...  First derivative (potential).
c...
      eps1p1 = -5.0D-1*(x*(1.97176D0*x12+4.9146D0)+7.1752D0*x12+7.5957D0
     1   )/(polx1*(6.21814D-2*polx1+1.0D0)*x12)
      epsp1 = -6.21814D-2*eps1p1*(2.137D-1*x+1.0D0)-1.328816518D-2*eps11
      eps1p2 = -5.0D-1*(x*(2.50068D0*x12+1.00986D1)+1.23954D1*x12+1.4118
     1   9D1)/(polx2*(3.10907D-2*polx2+1.0D0)*x12)
      epsp2 = -3.10907D-2*eps1p2*(2.0548D-1*x+1.0D0)-6.388517036D-3*eps1
     1   2
      eps1p3 = -5.0D-1*(x*(1.98684D0*x12+2.64078D0)+7.2462D0*x12+1.0357D
     1   1)/(polx3*(3.37738D-2*polx3+1.0D0)*x12)
      epsp3 = -3.37738D-2*eps1p3*(1.1125D-1*x+1.0D0)-3.75733525D-3*eps13
      omp = 2.5648814012420484263D0*(zp113-zm113)
      dra = 8.3849294190404081567D-2*rm43*(1.464591887561523263D0*eps3*r
     1   13*(9.5244063118091968485D0*rb*(omp*z4+4.0D0*om*z3-omp)+4
     2   .7622031559045984243D0*om*r*(z4-1.0D0))+1.464591887561523263D0*
     3   eps1*r13*(-1.6285982365095093684D1*rb*(omp*z4+4.0D0*om*z3)-8.14
     4   29911825475468422D0*r*(om*z4-1.0D0))+2.3852317652888304153D1*ep
     5   s2*r13*rb*(omp*z4+4.0D0*om*z3)+2.4661328275096140485D0*epsp1*r*
     6   (om*z4-1.0D0)-1.4422495703074083823D0*epsp3*om*r*(z4-1.0D0)+1.7
     7   09921D0*om*r*(6.9746841090577588534D0*eps2*r13-1.44224957030740
     8   83823D0*epsp2)*z4)
      drb = -8.3849294190404081567D-2*rm43*(1.464591887561523263D0*eps3*
     1   r13*(9.5244063118091968485D0*ra*(omp*z4+4.0D0*om*z3-omp)-
     2   4.7622031559045984243D0*om*r*(z4-1.0D0))+1.464591887561523263D0
     3   *eps1*r13*(8.1429911825475468422D0*r*(om*z4-1.0D0)-1.6285982365
     4   095093684D1*ra*(omp*z4+4.0D0*om*z3))+2.3852317652888304153D1*ep
     5   s2*r13*ra*(omp*z4+4.0D0*om*z3)-2.4661328275096140485D0*epsp1*r*
     6   (om*z4-1.0D0)+1.4422495703074083823D0*epsp3*om*r*(z4-1.0D0)-1.7
     7   09921D0*om*r*(6.9746841090577588534D0*eps2*r13-1.44224957030740
     8   83823D0*epsp2)*z4)
c...
      return
      end
c=======================================================================
      subroutine pw91nlcu(ra,rb,saa,sbb,sab,scl,scnl,ec,
     $                     dra,drb,dsaa,dsbb,dsab)
      implicit none !do not comment this line, no implicit convention here
c...  MM (05/21/2003) generated with maxima, file pw91nlcu.mc
c...
c...  *** IMPORTANT: remember to edit this file to:
c...           2) remove the coefficient 1 of some negative terms
c...
c...
c...  Perdew and Wang non local correlation functional
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...              open shell version
c...
c...  Formulation taken from the Molpro manual (PW91C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...  saa   alpha density gradient invariant
c...  sbb   beta density gradient invariant
c...  sab   alpha beta beta density gradient invariant
c...  scl   scaling factor for local part
c...  scnl  scaling factor for non local part
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dra     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...  dsab    funct. deriv. w.r.t. alpha beta gradient invariant
c...
      real*8 ra,rb,saa,sbb,sab,scl,scnl,ec,dra,drb,dsaa,dsbb,dsab
      real*8 pw91lc,dpw91lcdra,dpw91lcdrb
      real*8 r,r2,r32,r72,rm13,r13,rm76,r76,rm43,r136,rm196
      real*8 s,s12
      real*8 z
      real*8 zp1,zp113,zp123
      real*8 zm1,zm113,zm123
      real*8 x,x2,x3
      real*8 d,d2,d4
      real*8 polj,expj
      real*8 numlog,denlog
      real*8 om,omp,om2,om3,om4
      real*8 e,epra,eprb
      real*8 numt,dent,dentp,t,tp
      real*8 j,jpd,jpom,jpx,jpr
      real*8 djdra,djdrb,djdsaa
      real*8 expa,a,ape,apom
      real*8 loga,logape,logapd,logapom
      real*8 l,lpe,lpom,lpd
      real*8 dldra,dldrb,dldsaa
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,sixth,etm,th
      parameter (third=1.0d0/3.0d0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      call pw91lcu(ra,rb,pw91lc,dpw91lcdra,dpw91lcdrb)
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      r = rb+ra
      s = sbb+2.0D0*sab+saa
      if(s.lt.epsi.or.scnl.eq.zero)then
      j = ZERO
      djdra = ZERO
      djdrb = ZERO
      djdsaa = ZERO
      l = ZERO
      dldra = ZERO
      dldrb = ZERO
      dldsaa = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
      r2 = r**2
      r32 = r**th
      r72 = r2*r32
      r13 = r**THIRD
      rm13 = one/r13
      r76 = r72**THIRD
      rm76 = one/r76
      rm43 = rm13/r
      r136 = r*r76
      rm196 = one/(r*r136)
c...
      s12 = SQRT(s)
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x3 = x*x2
      z = MIN(MAX((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp123 = zp113**2
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm123 = zm113**2
c...
c...  the variables omega and d.
c...
      om = 5.0D-1*(zp123+zm123)
      om2 = om**2
      om3 = om*om2
      om4 = om*om3
      d = 2.519289703422449D-1*rm76*s12/om
      d2 = d**2
      d4 = d2**2
c...
c...  and the last auxiliary variable, e.
c...
      e = pw91lc
c...
c...  we start constructing the term j.
c...
c...  function t.
c...
      numt = 7.389D-3*x2+2.3266D1*x+2.568D0
      dent = 7.389D-2*x3+4.72D-1*x2+8.723D0*x+1.0D0
c...
c...  we build j.
c...
      t = 1.0D-3*numt/dent
      expj = EXP(-4.11563120990411453D1*d2*om4*rm13)
      polj = t-1.853268571428571D-3
      j = 1.575592034948315D1*d2*expj*om3*polj
c...
c...  now the term l
c...
c...  function a.
c...
      expa = EXP(-4.04276151175637226D1*e/om3)
      a = 2.697586091519874D0/(expa-1.0D0)
c...
c...  we build l.
c...
      numlog = a*d4+d2
      denlog = a**2*d4+a*d2+1.0D0
      loga = LOG(2.697586091519874D0*numlog/denlog+1.0D0)
      l = 2.473556743557577D-2*loga*om3
c...
c...  First derivative of om
c...
      if(zm113*zp113.gt.epsi)then
         omp = -3.333333333333333D-1*(zp113-zm113)/(zm113*zp113)
      else
        omp=zero
      endif
c...
c...  First derivatives of e
c...
      epra = (-e+dpw91lcdra)/r
      eprb = (-e+dpw91lcdrb)/r
c...
c...  First derivatives of j w.r.t. auxiliary variables
c...
      dentp = 2.2167D-1*x2+9.44D-1*x+8.723D0
      tp = 1.0D-3*(dent*(1.4778D-2*x+2.3266D1)-dentp*numt)/dent**2
      jpd = 2.057815604952057D-1*j*(9.719043801529514D0*r13-4.0D2*d2*om4
     1   )*rm13/d
      jpom = 1.028907802476029D-1*j*(2.915713140458854D1*r13-1.6D3*d2*om
     1   4)*rm13/om
      jpx = 1.575592034948315D1*d2*expj*om3*tp
      jpr = 1.371877069968038D1*d2*j*om4*rm43
c...
c...  First derivatives of l w.r.t. auxiliary variables
c...
      ape = 1.090569722544585D2*expa/((expa-1.0D0)**2*om3)
      apom = -3.271709167633755D2*e*expa/((expa-1.0D0)**2*om4)
      logape = -1.8D-1*ape*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*denlog
     1   )/(denlog*(1.8D-1*numlog+6.672632268006112D-2*denlog))
      logapom = -1.8D-1*apom*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*denl
     1   og)/(denlog*(1.8D-1*numlog+6.672632268006112D-2*denlog))
      logapd = -3.6D-1*d*(2.0D0*a*d2+1.0D0)*(a*numlog-denlog)/(den
     1   log*(1.8D-1*numlog+6.672632268006112D-2*denlog))
      lpe = 2.473556743557577D-2*logape*om3
      lpom = 2.473556743557577D-2*(logapom*om+3.0D0*loga)*om2
      lpd = 2.473556743557577D-2*logapd*om3
c...
c...  Total first derivatives of the non local terms:
c...  1) j
c...  2) l
c...
      djdra = -1.492331345469521D-2*(2.813595107492502D0*jpd*r13*(1.2D1*
     1   omp*rb+7.0D0*om*r)*s12-6.700924717797021D1*om2*r32*(2.0D0*jpom*
     2   omp*rb+jpr*r2)+1.385640646055102D1*jpx*om2*r136)/(om2*r72)
      djdrb = 1.492331345469521D-2*(2.813595107492502D0*jpd*r13*(1.2D1*o
     1   mp*ra-7.0D0*om*r)*s12-6.700924717797021D1*om2*r32*(2.0D0*jpom*o
     2   mp*ra-jpr*r2)-1.385640646055102D1*jpx*om2*r136)/(om2*r72)
      djdsaa = 1.259644851711224D-1*jpd*rm76/(om*s12)
      dldra = -3.469513240231685D-2*rm196*(1.210203242253764D0*lpd*(1.2D
     1   1*omp*rb+7.0D0*om*r)*s12-2.882248692422407D1*om2*r76*(2.0D0*lpo
     2   m*omp*rb+epra*lpe*r2))/om2
      dldrb = 3.469513240231685D-2*rm196*(1.210203242253764D0*lpd*(1.2D1
     1   *omp*ra-7.0D0*om*r)*s12-2.882248692422407D1*om2*r76*(2.0D0*lpom
     2   *omp*ra-eprb*lpe*r2))/om2
      dldsaa = 1.259644851711224D-1*lpd*rm76/(om*s12)
c...
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = (l+j)*scnl+pw91lc*scl
c...
c...  First derivatives of the functional (potential).
c...
      dra = ((dldra+djdra)*r+l+j)*scnl+dpw91lcdra*scl
      drb = ((dldrb+djdrb)*r+l+j)*scnl+dpw91lcdrb*scl
      dsaa = (dldsaa+djdsaa)*r*scnl
      dsbb = dsaa
      dsab = 2.0D0*dsaa
c...
      return
      end
c=======================================================================
      subroutine pw91xu(ra,rb,saa,sbb,exc,dra,drb,dsaa,dsbb)
      implicit none
c...  MM (01/27/2004) generated with maxima, file pw91xu.mc
c...
c...                   Perdew-Wang 1991 (PW91)
c...                     Exchange Fuctional
c...                    open-shell version
c...
c...  see Kieron Burke, John P Perdew, and Yuc Wang, in Electronic
c...  Density Functional Theory: Recent Progress and New Directions,
c...  eds. J.F. Dobson, G. Vignale, and M.P. Das (Plenum, NY, 1997),
c...  page 81.
c...
c...       Formulation taken from the Molpro manual (PW91X)
c...         http://www.molpro.net/current/molpro_manual
c...
c...  This subroutine computes the functional and its first derivatives
c...  (minus the Slater exchange contribution, which is computed separately).
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...  saa   alpha density gradient invariant
c...  sbb   beta density gradient invariant
c...
c...  Output parameters:
c...
c...  exc     contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     functional derivative with respect to beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. beta gradient invariant
c...
c...  Note that this functional does not have mixed alpha-beta
c...  derivatives
c...
      real*8 ra,rb,saa,sbb,exc,eca,esla,ecb,eslb,dra,d1ra,dslra,drb,d1rb
     $  ,dslrb,dsaa,dsbb
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 rbm1,rb13,rb43,rbm43,sbb12,sbbm12
      real*8 essea,essea2,essea3,essea4
      real*8 esseb,esseb2,esseb3,esseb4
      real*8 cessea1,cessea12,cesseam12
      real*8 cesseb1,cesseb12,cessebm12
      real*8 nfa,dfa,nfap,dfap,dfa2,dfam2
      real*8 nfb,dfb,nfbp,dfbp,dfb2,dfbm2
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third
      parameter (third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...
c...  Alpha part
c...
      if(ra.gt.epsi)then
c...
c...  Some povers of ra and saa
c...
      ram1 = one/ra
      ra13 = ra**THIRD
      ra43 = ra*ra13
      ram43 = one/ra43
      saa12 = SQRT(saa)
c...
c...  We compute essea, and some of its powers
c...
      essea = 1.282782438530422D-1*ram43*saa12
      essea2 = essea**2
      essea3 = essea*essea2
      essea4 = essea*essea3
c...
c...  Function effe
c...
      cessea1 = 6.077137936D1*essea2+1.0D0
      cessea12 = SQRT(cessea1)
      nfa = essea2*(2.743D-1-1.508D-1*EXP(-1.0D2*essea2))+1.9645D-1*esse
     $  a*LOG(7.7956D0*essea+cessea12)+1.0D0
      dfa = 4.0D-3*essea4+1.9645D-1*essea*LOG(7.7956D0*essea+cessea12)+1
     $  .0D0
c...
c...  Functional value
c...
      eca = -9.305257363491D-1*nfa*ra43/dfa
      esla = -9.305257363491D-1*ra43
c...
c...  First derivatives
c...
      cesseam12 = one/cessea12
      nfap = cesseam12*(3.016D-1*cessea12*essea*(1.0D1*essea-1.0D0)*(1.0
     $  D1*essea+1.0D0)*EXP(-1.0D2*essea2)+cessea12*(1.9645D-1*LOG(7.795
     $  6D0*essea+cessea12)+5.486D-1*essea)+1.53144562D0*essea)
      dfa2 = dfa**2
      dfam2 = one/dfa2
      dfap = 1.6D-2*essea3+1.9645D-1*LOG(7.7956D0*essea+cessea12)+1.5314
     $  4562D0*cesseam12*essea
c...
c...  First derivatives of the functional
c...
      d1ra = 1.086684587314497D-1*dfam2*ram1*(1.464591887561523D0*dfa*nf
     $  ap*saa12-1.464591887561523D0*dfap*nfa*saa12-1.141730541025636D1*
     $  dfa*nfa*ra43)
      dslra = -1.2407009817988D0*ra13
      if(saa.gt.epsi)then
      saam12 = one/saa12
      dsaa = -5.968310365946075D-2*dfam2*(dfa*nfap-dfap*nfa)*saam12
      else
      dsaa = ZERO
      endif
      else
      eca = ZERO
      esla = ZERO
      d1ra = ZERO
      dslra = ZERO
      dsaa = ZERO
      endif
c...
c...  Beta part
c...
      if(rb.gt.epsi)then
c...
c...  Some povers of rb and sbb
c...
      rbm1 = one/rb
      rb13 = rb**THIRD
      rb43 = rb*rb13
      rbm43 = one/rb43
      sbb12 = SQRT(sbb)
c...
c...  We compute esseb and some of its powers
c...
      esseb = 1.282782438530422D-1*rbm43*sbb12
      esseb2 = esseb**2
      esseb3 = esseb*esseb2
      esseb4 = esseb*esseb3
c...
c...  Function effe
c...
      cesseb1 = 6.077137936D1*esseb2+1.0D0
      cesseb12 = SQRT(cesseb1)
      nfb = esseb2*(2.743D-1-1.508D-1*EXP(-1.0D2*esseb2))+1.9645D-1*esse
     $  b*LOG(7.7956D0*esseb+cesseb12)+1.0D0
      dfb = 4.0D-3*esseb4+1.9645D-1*esseb*LOG(7.7956D0*esseb+cesseb12)+1
     $  .0D0
c...
c...  Functional value
c...
      ecb = -9.305257363491D-1*nfb*rb43/dfb
      eslb = -9.305257363491D-1*rb43
c...
c...  First derivatives
c...
      cessebm12 = one/cesseb12
      nfbp = cessebm12*(3.016D-1*cesseb12*esseb*(1.0D1*esseb-1.0D0)*(1.0
     $  D1*esseb+1.0D0)*EXP(-1.0D2*esseb2)+cesseb12*(1.9645D-1*LOG(7.795
     $  6D0*esseb+cesseb12)+5.486D-1*esseb)+1.53144562D0*esseb)
      dfb2 = dfb**2
      dfbm2 = one/dfb2
      dfbp = 1.6D-2*esseb3+1.9645D-1*LOG(7.7956D0*esseb+cesseb12)+1.5314
     $  4562D0*cessebm12*esseb
c...
c...  First derivatives of the functional
c...
      d1rb = 1.086684587314497D-1*dfbm2*rbm1*(1.464591887561523D0*dfb*nf
     $  bp*sbb12-1.464591887561523D0*dfbp*nfb*sbb12-1.141730541025636D1*
     $  dfb*nfb*rb43)
      dslrb = -1.2407009817988D0*rb13
      if(sbb.gt.epsi)then
      sbbm12 = one/sbb12
      dsbb = -5.968310365946075D-2*dfbm2*(dfb*nfbp-dfbp*nfb)*sbbm12
      else
      dsbb = ZERO
      endif
      else
      ecb = ZERO
      eslb = ZERO
      d1rb = ZERO
      dslrb = ZERO
      dsbb = ZERO
      endif
c...
c...  Assemble the output values, subtracting the Slater term
c...
      exc = (-eslb-esla+ecb+eca)/(rb+ra)
      dra = d1ra-dslra
      drb = d1rb-dslrb
c...
      return
      end
c=======================================================================
      subroutine pw91x(ra,saa,exc,dra,dsaa)
      implicit none
c...  MM (01/22/2004) generated with maxima, file pw91xc.mc
c...
c...                   Perdew-Wang 1991 (PW91)
c...                     Exchange Fuctional
c...                    closed shell version
c...
c...  see Kieron Burke, John P Perdew, and Yuc Wang, in Electronic
c...  Density Functional Theory: Recent Progress and New Directions,
c...  eds. J.F. Dobson, G. Vignale, and M.P. Das (Plenum, NY, 1997),
c...  page 81.
c...
c...       Formulation taken from the Molpro manual (PW91X)
c...         http://www.molpro.net/current/molpro_manual
c...
c...  This subroutine computes the functional and its first derivatives
c...  (minus the Slater exchange contribution, which is computed separately).
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  saa   alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  exc     contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,exc,dra,dsaa,eca,esla,d1ra,dslra
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 essea,essea2,essea3,essea4
      real*8 cessea1,cessea12,cesseam12
      real*8 nf,df,nfp,dfp,df2,dfm2
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third
      parameter (third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...
c...  Some povers of ra and saa
c...
      ram1 = one/ra
      ra13 = ra**THIRD
      ra43 = ra*ra13
      ram43 = one/ra43
      saa12 = SQRT(saa)
c...
c...  We compute esse and some of its powers
c...
      essea = 1.282782438530422D-1*ram43*saa12
      essea2 = essea**2
      essea3 = essea*essea2
      essea4 = essea*essea3
c...
c...  Function effe
c...
      cessea1 = 6.077137936D1*essea2+1.0D0
      cessea12 = SQRT(cessea1)
      nf = essea2*(2.743D-1-1.508D-1*EXP(-1.0D2*essea2))+1.9645D-1*essea
     $  *LOG(7.7956D0*essea+cessea12)+1.0D0
      df = 4.0D-3*essea4+1.9645D-1*essea*LOG(7.7956D0*essea+cessea12)+1.
     $  0D0
c...
c...  Functional value
c...
      eca = -9.305257363491D-1*nf*ra43/df
      esla = -9.305257363491D-1*ra43
c...
c...  First derivatives
c...
      cesseam12 = one/cessea12
      nfp = cesseam12*(3.016D-1*cessea12*essea*(1.0D1*essea-1.0D0)*(1.0D
     $  1*essea+1.0D0)*EXP(-1.0D2*essea2)+cessea12*(1.9645D-1*LOG(7.7956
     $  D0*essea+cessea12)+5.486D-1*essea)+1.53144562D0*essea)
      df2 = df**2
      dfm2 = one/df2
      dfp = 1.6D-2*essea3+1.9645D-1*LOG(7.7956D0*essea+cessea12)+1.53144
     $  562D0*cesseam12*essea
c...
c...  First derivatives of the functional
c...
      d1ra = 1.086684587314497D-1*dfm2*ram1*(1.464591887561523D0*df*nfp*
     $  saa12-1.464591887561523D0*dfp*nf*saa12-1.141730541025636D1*df*nf
     $  *ra43)
      dslra = -1.2407009817988D0*ra13
      if(saa.gt.epsi)then
      saam12 = one/saa12
      dsaa = -5.968310365946075D-2*dfm2*(df*nfp-dfp*nf)*saam12
      else
      dsaa = ZERO
      endif
c...
c...  Assemble the output values, subtracting the Slater term
c...
      exc = (eca-esla)/ra
      dra = d1ra-dslra
c...
      return
      end
c=======================================================================
      subroutine pbexu(ra,rb,saa,sbb,exc,dra,drb,dsaa,dsbb)
      implicit none
c...
c...  MM (01/28/2004) generated with maxima, file pbexu.mc
c...
c...  Perdew-Burke-Ernzerhof non-local exchange functional (PBEX)
c...
c...       J.P. Perdew, Kieron Burke, Matthias Ernzerhof,
c...              Phys. Rev. Lett. 77, 3865 (1996).
c...
c...                   open shell version
c...
c...  This subroutine computes the functional and its first derivatives
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...  saa   alpha density gradient invariant
c...  sbb   beta density gradient invariant
c...
c...  Output parameters:
c...
c...  exc     contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     functional derivative with respect to beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. beta gradient invariant
c...
c...  Note that this functional does not have mixed alpha-beta
c...  derivatives
c...
      real*8 ra,rb,saa,sbb,exc,exa,exb,dra,drb,dsaa,dsbb
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 rbm1,rb13,rb43,rbm43,sbb12,sbbm12
      real*8 dfam1,dfam2
      real*8 dfbm1,dfbm2
      real*8 essea,essea2
      real*8 esseb,esseb2
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third
      parameter (third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...
c...  Alpha part
c...
      if(ra.gt.epsi)then
c...
c...  Some povers of ra and saa
c...
      ram1 = one/ra
      ra13 = ra**THIRD
      ra43 = ra*ra13
      ram43 = one/ra43
      saa12 = SQRT(saa)
c...
c...  We compute essea, and some of its powers
c...
      essea = 1.282782438530422D-1*ram43*saa12
      essea2 = essea**2
c...
c...  Functional value
c...
      dfam1 = 1/(2.730304119662883D-1*essea2+1.0D0)
      exa = -9.305257363491D-1*(8.04D-1-8.04D-1*dfam1)*ra43
c...
c...  First derivatives
c...
      dfam2 = dfam1**2
      dra = 6.987425660359299D-2*dfam2*essea*ram1*saa12-1.2407009817988D
     $  0*(8.04D-1-8.04D-1*dfam1)*ra13
      if(saa.gt.epsi)then
      saam12 = one/saa12
      dsaa = -2.620284622634737D-2*dfam2*essea*saam12
      else
      dsaa = ZERO
      endif
      else
      exa = ZERO
      dra = ZERO
      dsaa = ZERO
      endif
c...
c...  Beta part
c...
      if(rb.gt.epsi)then
c...
c...  Some povers of rb and sbb
c...
      rbm1 = one/rb
      rb13 = rb**THIRD
      rb43 = rb*rb13
      rbm43 = one/rb43
      sbb12 = SQRT(sbb)
c...
c...  We compute esseb and some of its powers
c...
      esseb = 1.282782438530422D-1*rbm43*sbb12
      esseb2 = esseb**2
c...
c...  Functional value
c...
      dfbm1 = 1/(2.730304119662883D-1*esseb2+1.0D0)
      exb = -9.305257363491D-1*(8.04D-1-8.04D-1*dfbm1)*rb43
c...
c...  First derivatives
c...
      dfbm2 = dfbm1**2
      drb = 6.987425660359299D-2*dfbm2*esseb*rbm1*sbb12-1.2407009817988D
     $  0*(8.04D-1-8.04D-1*dfbm1)*rb13
      if(sbb.gt.epsi)then
      sbbm12 = one/sbb12
      dsbb = -2.620284622634737D-2*dfbm2*esseb*sbbm12
      else
      dsbb = ZERO
      endif
      else
      exb = ZERO
      drb = ZERO
      dsbb = ZERO
      endif
c...
c...  functional output value
c...
      exc = (exb+exa)/(rb+ra)
c...
      return
      end
c=======================================================================
      subroutine pbecu(ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab)
      implicit none
c...  MM (01/28/2004) generated with maxima, file pbecu.mc
c...
c...  Perdew-Burke-Ernzerhof non-local correlation functional (PBEC)
c...
c...       J.P. Perdew, Kieron Burke, Matthias Ernzerhof,
c...              Phys. Rev. Lett. 77, 3865 (1996).
c...
c...                   open shell version
c...
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...  saa   alpha density gradient invariant
c...  sbb   beta density gradient invariant
c...  sab   alpha beta beta density gradient invariant
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dra     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...  dsab    funct. deriv. w.r.t. alpha beta gradient invariant
c...
      real*8 ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab
      real*8 pw91lc,dpw91lcdra,dpw91lcdrb
      real*8 r,r2,r32,r72,rm13,r13,rm76,r76,rm43,r136,rm196
      real*8 s,s12
      real*8 z
      real*8 zp1,zp113,zp123
      real*8 zm1,zm113,zm123
      real*8 d,d2,d4
      real*8 numlog,denlog
      real*8 om,omp,om2,om3,om4
      real*8 e,epra,eprb
      real*8 expa,a,ape,apom
      real*8 loga,logape,logapd,logapom
      real*8 h,hpe,hpom,hpd
      real*8 dhdra,dhdrb,dhdsaa
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,sixth,etm,th
      parameter (third=1.0d0/3.0d0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      call pw91lcu(ra,rb,pw91lc,dpw91lcdra,dpw91lcdrb)
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      r = rb+ra
      s = sbb+2.0D0*sab+saa
      if(s.lt.epsi)then
      h = ZERO
      dhdra = ZERO
      dhdrb = ZERO
      dhdsaa = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
      r = rb+ra
      r2 = r**2
      r32 = r**th
      r72 = r2*r32
      r13 = r**THIRD
      rm13 = one/r13
      r76 = r72**THIRD
      rm76 = one/r76
      rm43 = rm13/r
      r136 = r*r76
      rm196 = one/(r*r136)
c...
      s12 = SQRT(s)
c...
c...  the auxiliary variables z and some of its powers.
c...
      z = MIN(MAX((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp123 = zp113**2
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm123 = zm113**2
c...
c...  the variables omega and d.
c...
      om = 5.0D-1*(zp123+zm123)
      om2 = om**2
      om3 = om*om2
      om4 = om*om3
      d = 2.519289703422449D-1*rm76*s12/om
      d2 = d**2
      d4 = d2**2
c...
c...  and the last auxiliary variable, e.
c...
      e = pw91lc
c...
c...  the term h
c...
c...  function a.
c...
      expa = EXP(-3.21639684429148211D1*e/om3)
      a = 2.146140794353492D0/(expa-1.0D0)
c...
c...  we build h.
c...
      numlog = a*d4+d2
      denlog = a**2*d4+a*d2+1.0D0
      loga = LOG(2.146140794353492D0*numlog/denlog+1.0D0)
      h = 3.10906908696549D-2*loga*om3
c...
c...  First derivative of om
c...
      if(zm113*zp113.gt.epsi)then
      omp = -3.333333333333333D-1*(zp113-zm113)/(zm113*zp113)
      else
      omp = ZERO
      endif
c...
c...  First derivatives of e
c...
      epra = -(e-dpw91lcdra)/r
      eprb = -(e-dpw91lcdrb)/r
c...
c...  First derivatives of h w.r.t. auxiliary variables
c...
      ape = 6.902840478363784D1*expa/((expa-1.0D0)**2*om3)
      apom = -2.070852143509135D2*e*expa/((expa-1.0D0)**2*om4)
      logape = -6.6725D-2*ape*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*den
     $  log)/(denlog*(6.6725D-2*numlog+3.10906908696549D-2*denlog))
      logapom = -6.6725D-2*apom*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*d
     $  enlog)/(denlog*(6.6725D-2*numlog+3.10906908696549D-2*denlog))
      logapd = -1.3345D-1*d*(2.0D0*a*d2+1.0D0)*(a*numlog-denlog)/(
     $  denlog*(6.6725D-2*numlog+3.10906908696549D-2*denlog))
      hpe = 3.10906908696549D-2*logape*om3
      hpom = 3.10906908696549D-2*(logapom*om+3.0D0*loga)*om2
      hpd = 3.10906908696549D-2*logapd*om3
c...
c...  Total first derivatives of the non local term h
c...
      dhdra = -3.469513240231685D-2*rm196*(1.210203242253764D0*hpd*(1.2D
     $  1*omp*rb+7.0D0*om*r)*s12-2.882248692422407D1*om2*r76*(2.0D0*hpom
     $  *omp*rb+epra*hpe*r2))/om2
      dhdrb = 3.469513240231685D-2*rm196*(1.452243890704517D1*hpd*omp*ra
     $  *s12-8.47142269577635D0*hpd*om*r*s12-5.764497384844813D1*hpom*om
     $  2*omp*r76*ra+2.882248692422407D1*eprb*hpe*om2*r**(19.0/6.0))/om2
      dhdsaa = 1.259644851711224D-1*hpd*rm76/(om*s12)
c...
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = pw91lc+h
c...
c...  First derivatives of the functional (potential).
c...
      dra = dhdra*r+h+dpw91lcdra
      drb = dhdrb*r+h+dpw91lcdrb
      dsaa = dhdsaa*r
      dsbb = dsaa
      dsab = 2.0D0*dsaa
c...
      return
      end

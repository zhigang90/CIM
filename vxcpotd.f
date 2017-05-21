      SUBROUTINE vxcpotD(dft,    NP,     thrsh,  Den,    DenX,
     $                   WGHT,   pra,    prara,  prarb,  pga,
     $                   pgc,    praga,  pragb,  pragc,  pgaga,
     $                   pgagb,  pgagc,  pgcgc,  ETot,   EL,
     $                   EVec,   IdWt)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Evaluates the exchange-correlation potential and
C  potential derivative, and determines
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
C  Den     -  density
C  DenX    -  density gradient (dft > 3 only)
C             ** IMPORTANT - The latter has been predivided by 2 **
C  WGHT    -  grid quadrature weights
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C
C  on exit
C
C  pra     -  Functional derivative w.r.t. alpha density at grid points
C  prara   -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  ETot    -  total exchange-correlation energy contribution
C  EL      -  number of electrons over grid
C  EVec    -  vector of exchange-correlation energy contributions
C             per grid point (computed only if (IdWt.eq.1)
C             (needed for quadrature weight derivatives in DFT hessian)
C
C  For nonlocal functionals only ( dft > 3)
C
C  prarb -   Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -   Funct. deriv. w.r.t. alpha gradient
C  pgc   -   Funct. deriv. w.r.t. alpha beta gradient
C  praga -   Funct. 2nd. deriv. w.r.t. alpha dens. and  grad.
C  pragb -   Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C  pragc -   Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C  pgaga -   Funct. 2nd. deriv. w.r.t. alpha grad.
C  pgagb -   Funct. 2nd. deriv. w.r.t. alpha and beta grad.
C  pgagc -   Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C  pgcgc -   Funct. 2nd. deriv. w.r.t. alpha beta grad.
C
C  NOTE: For the closed shell case, Den contains the total
C        density, while DenX contains half of the total
C        density gradient. The potentials (and potential derivatives)
C        returned by this subroutine are computed with respect
C        to half the density and half the density gradient invariant
C
C  NOTE ALSO: For non local functionals, the potentials (and
C             derivatives) returned by this subroutine are computed
C             with respect to the gradient invariant.
C             This is a difference with the PQS implementation of
C             the SCF and force calculation, where subroutines
C             vxcpot and vxcpotu compute potentials with respect to
C             the density gradient.
C
      DIMENSION Den(NP),DenX(3,NP)
      DIMENSION WGHT(NP),EVec(NP)
      INTEGER dft
      dimension pra(np),pga(np),pgc(np)
      dimension prara(np),prarb(np)
      dimension praga(np),pragb(np),pragc(np)
      dimension pgaga(np),pgagb(np),pgagc(np),pgcgc(np)
C
      parameter (Zero=0.0d0,Third=1.0d0/3.0d0,Half=0.5d0,Four=4.0d0)
cc      parameter (PI=3.14159 26535 89793d0)
      parameter (two=2.0d0,TThird=2.0d0/3.0d0)
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
       pra(I) = Zero
       prara(I) = Zero
       If(dft.GT.3) Then
         prarb(I) = Zero
         pga(I) = Zero
         pgc(I) = Zero
         praga(I) = Zero
         pragb(I) = Zero
         pragc(I) = Zero
         pgaga(I) = Zero
         pgagb(I) = Zero
         pgagc(I) = Zero
         pgcgc(I) = Zero
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
        ra=zero
        rara=zero
        ga = zero
        gc = zero
        rarb= zero
        raga = zero
        ragb = zero
        ragc = zero
        gaga = zero
        gagb = zero
        gagc = zero
        gcgc = zero
c
       If(dft.EQ.24) Then
c...
c...  Becke 97.
c...
        call b97dr(halfden,g2,ec,dra,dsaa,drara,drasaa,dsaasaa)
        exc=ec
        ra=dra
        rara=drara
        ga = dsaa
        raga = drasaa
        gaga = dsaasaa
cc
       Else If(dft.EQ.25) Then
c...
c...  Becke 97_1.
c...
        call b971dr(halfden,g2,ec,dra,dsaa,drara,drasaa,dsaasaa)
        exc=ec
        ra=dra
        rara=drara
        ga = dsaa
        raga = drasaa
        gaga = dsaasaa
cc
       Else If(dft.EQ.26) Then
c...
c...  Becke 97_2.
c...
        call b972dr(halfden,g2,ec,dra,dsaa,drara,drasaa,dsaasaa)
        exc=ec
        ra=dra
        rara=drara
        ga = dsaa
        raga = drasaa
        gaga = dsaasaa
cc
       Else
c...
c...  HCTH.
c...
        call hcth407dr(halfden,g2,ec,dra,dsaa,drara,drasaa,dsaasaa)
        exc=ec
        ra=dra
        rara=drara
        ga = dsaa
        raga = drasaa
        gaga = dsaasaa
cc
       EndIf
      ELSE
C
C  get the local Slater exchange
C
        Dens13 = Dens**Third
        EXC = XS*EFSC*Dens13
        ra = XS*VFSC*Dens13
        rara = TThird*XS*VFSC/(Dens13**2)
        If(dft.EQ.1) GO TO 10
C
C  VWN local Correlation
C
        IF(XVWN.NE.Zero) THEN
          CALL VWND(Dens,Jnk,.true.,ECVWN,VCVWN,Jnk,VCVWND,Jnk,Jnk)
          EXC = EXC + ECVWN*XVWN
          ra = ra + VCVWN*XVWN
          rara= rara + VCVWND*XVWN
          If(dft.EQ.2) GO TO 10
        ELSE IF(XVWN5.NE.Zero) THEN
          CALL VWN5D(Dens,Jnk,.true.,ECVWN,VCVWN,Jnk,VCVWND,Jnk,Jnk)
          EXC = EXC + ECVWN*XVWN5
          ra = ra + VCVWN*XVWN5
          rara= rara + VCVWND*XVWN5
          If(dft.EQ.3) GO TO 10
        ENDIF
c
        HalfDen = Half*Dens
        Dens13 = HalfDen**Third
        Dens43 = HalfDen*Dens13
c
        g2 = (DenX(1,I)**2 + DenX(2,I)**2 + DenX(3,I)**2)
        g1 = SQRT(g2)
        ga = zero
        gc = zero
        rarb= zero
        raga = zero
        ragb = zero
        ragc = zero
        gaga = zero
        gagb = zero
        gagc = zero
        gcgc = zero
c
        IF(XB88.NE.Zero) THEN
C
C  Becke non-local exchange
C
          Call Becke88D(HalfDen,Dens13, Dens43, g1,    g2,
     $                  unk,    unk,    unk,    unk,   unk,
     $    .true., EXB88,V1A,    V2A,    VV1A,   VV2A,  VV3A,
     $                  unk,    unk,    unk,    unk,   unk)
          EXC = EXC + EXB88*XB88
          ra = ra + V1A*XB88
          ga = ga + V2A*XB88
          rara = rara + VV1A*XB88
          gaga = gaga + VV2A*XB88
          raga = raga + VV3A*XB88
        ENDIF
C
C  Handy/Cohen optimized non-local exchange
C
        IF(XOptX.NE.Zero) THEN
          Call OPTXD(  HalfDen,Dens13, Dens43, g1,    g2,
     $                 unk,    unk,    unk,    unk,   unk,
     $   .true.,EXOPTX,V1A,    V2A,    VV1A,   VV2A,  VV3A,
     $                 unk,    unk,    unk,    unk,   unk)
          EXC = EXC + EXOPTX*XOptX
          ra = ra + V1A*XOptX
          ga = ga + V2A*XOptX
          rara = rara + VV1A*XOptX
          gaga = gaga + VV2A*XOptX
          raga = raga + VV3A*XOptX
        ENDIF
c
c  Perdew 86 non-local correlation
c
       IF(XP86L.NE.Zero.OR.XP86NL.NE.Zero) THEN
         call p86nlcdu(halfden,halfden,g2,g2,g2,xp86l,xp86nl,ecp86,
     $                    dra,drb,dsaa,dsbb,dsab,
     $                    drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                    drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $                    dsabsab,dsaasbb,dsaasab,dsbbsab)
         EXC = EXC + ECP86
          ra=ra+dra
          ga=ga+dsaa
          gc=gc+dsab
          rara=rara+drara
          rarb=rarb+drarb
          raga=raga+drasaa
          ragb=ragb+drasbb
          ragc=ragc+drasab
          gaga=gaga+dsaasaa
          gagb=gagb+dsaasbb
          gagc=gagc+dsaasab
          gcgc=gcgc+dsabsab
       ENDIF
C
C  Perdew-Wang 91 non-local exchange
C
       IF(XPW91X.NE.Zero) THEN
         call pw91xd(halfden,g2,expw91,dra,dsaa,drara,drasaa,dsaasaa)
         EXC = EXC + expw91*xpw91x
         ra = ra + dra*xpw91x
         ga = ga + dsaa*xpw91x
         rara = rara + drara*xpw91x
         raga = raga + drasaa*xpw91x
         gaga = gaga + dsaasaa*xpw91x
       ENDIF
C
C  Perdew-Wang non-local correlation
        IF(XPW91L.NE.Zero.OR.XPW91NL.NE.Zero) THEN
          call pw91nlcdu(halfden,halfden,g2,g2,g2,XPW91L,XPW91NL,ecpw91,
     $                         dra,drb,dsaa,dsbb,dsab,
     $                         drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                         drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $                         dsabsab,dsaasbb,dsaasab,dsbbsab)
          EXC = EXC + ecpw91
          ra=ra+dra
          ga=ga+dsaa
          gc=gc+dsab
          rara=rara+drara
          rarb=rarb+drarb
          raga=raga+drasaa
          ragb=ragb+drasbb
          ragc=ragc+drasab
          gaga=gaga+dsaasaa
          gagb=gagb+dsaasbb
          gagc=gagc+dsaasab
          gcgc=gcgc+dsabsab
        ENDIF
C
C  Lee, Yang & Parr non-local correlation
C
C   this functional is linear on the density gradient, thus
C   its second derivatives with respect to the density gradients
C   two times are zero
C
        IF(XLYP.NE.Zero) THEN
          CALL LYPCD(HalfDen,Dens13,Dens43,g2,ECLYP,dra,dgaa,dgab,
     $                  dra2,drab,dragaa,dragbb,dragab)
          EXC = EXC + ECLYP*XLYP
          ra=ra+dra*XLYP
          ga=ga+dgaa*XLYP
          gc=gc+dgab*XLYP
          rara=rara+dra2*XLYP
          rarb=rarb+drab*XLYP
          raga=raga+dragaa*XLYP
          ragb=ragb+dragbb*XLYP
          ragc=ragc+dragab*XLYP
        ENDIF
C
C  Perdew, Burke and Ernzerhof non local exchange
C
        IF(XPBEX.NE.Zero) THEN
          call pbexdu(halfden,halfden,g2,g2,excpbe,dra,drb,dsaa,dsbb,
     $                drara,drbrb,drasaa,drbsbb,dsaasaa,dsbbsbb)
          EXC=EXC+excpbe*xpbex
          ra = ra + dra*xpbex
          ga = ga + dsaa*xpbex
          rara = rara + drara*xpbex
          raga = raga + drasaa*xpbex
          gaga = gaga + dsaasaa*xpbex
        ENDIF
C
C  Perdew, Burke and Ernzerhof non local correlation
C
        IF(XPBEC.NE.Zero) THEN
          call pbecdu(halfden,halfden,g2,g2,g2,ecpbe,dra,drb,dsaa,dsbb,
     $                dsab,drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,dsabsab,
     $                dsaasbb,dsaasab,dsbbsab)
          EXC=EXC+ecpbe*xpbec
          ra=ra+dra*xpbec
          ga=ga+dsaa*xpbec
          gc=gc+dsab*xpbec
          rara=rara+drara*xpbec
          rarb=rarb+drarb*xpbec
          raga=raga+drasaa*xpbec
          ragb=ragb+drasbb*xpbec
          ragc=ragc+drasab*xpbec
          gaga=gaga+dsaasaa*xpbec
          gagb=gagb+dsaasbb*xpbec
          gagc=gagc+dsaasab*xpbec
          gcgc=gcgc+dsabsab*xpbec
        ENDIF
C
C ..................................................
C  other functional density derivative code here
C ..................................................
C
C
      ENDIF
      pga(i) = ga
      pgc(i) = gc
      prarb(I) = rarb
      praga(I) = raga
      pragb(I) = ragb
      pragc(I) = ragc
      pgaga(I) = gaga
      pgagb(I) = gagb
      pgagc(I) = gagc
      pgcgc(I) = gcgc
C
 10   CONTINUE
C
C  Accumulate into total Exchange-Correlation energy
C  and total number of electrons
C
      DVal = Dens*WGHT(I)
      ETot = ETot + EXC*DVal
      EL   = EL + DVal
C
C  store full potential and exchange-correlation energy
C
      pra(I) = ra
      prara(I) = rara
      if(IdWt.eq.1) EVec(I) = EXC*dens ! for weight derivatives
cc      EVec(I) = EXC           ! for Malkin
cc
 20   CONTINUE
C
      RETURN
      END
C======================================================================
      SUBROUTINE vxcpotdu(dft,    NP,     thrsh,  DenA,   DenB,
     $                    DenXA,  DenXB,  WGHT,   pra,    prb,
     $                    prara,  prbrb,  prarb,  pga,    pgb,
     $                    pgc,    praga,  pragb,  pragc,  prbga,
     $                    prbgb,  prbgc,  pgaga,  pgagb,  pgagc,
     $                    pgbgb,  pgbgc,  pgcgc,  ETot,   EL,
     $                    EVec,   IdWt)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Evaluates the exchange-correlation potential and
C  potential derivative, and determines
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
C  DenA    -  alpha density
C  DenB    -  beta density
C  DenXA   -  alpha density gradient (dft > 3 only)
C  DenXB   -  beta density gradient (dft > 3 only)
C  WGHT    -  grid quadrature weights
C  IdWt    -  flag controlling whether to include weight derivatives
C               0 - do not take weight derivatives
C               1 - take quadrature weight derivatives
C
C  on exit
C
C  pra     -  Functional derivative w.r.t. alpha density at grid points
C  prb     -  Functional derivative w.r.t. beta density at grid points
C  prara   -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb   -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb   -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  ETot    -  total exchange-correlation energy contribution
C  EL      -  number of electrons over grid
C  EVec    -  vector of exchange-correlation energy contributions
C             per grid point (computed only if (IdWt.eq.1)
C             (needed for quadrature weight derivatives in DFT hessian)
C
C  For nonlocal functionals only ( dft > 3)
C
C  pga   -   Funct. deriv. w.r.t. alpha gradient
C  pgb   -   Funct. deriv. w.r.t. beta gradient
C  pgc   -   Funct. deriv. w.r.t. alpha beta gradient
C  praga -   Funct. 2nd. deriv. w.r.t. alpha dens. and  grad.
C  pragb -   Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C  pragc -   Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C  prbga -   Funct. 2nd. deriv. w.r.t. beta dens. and  alpha grad.
C  prbgb -   Funct. 2nd. deriv. w.r.t. beta dens. and  grad.
C  prbgc -   Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C  pgaga -   Funct. 2nd. deriv. w.r.t. alpha grad.
C  pgagb -   Funct. 2nd. deriv. w.r.t. alpha and beta grad.
C  pgagc -   Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C  pgbgb -   Funct. 2nd. deriv. w.r.t. beta grad.
C  pgbgc -   Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C  pgcgc -   Funct. 2nd. deriv. w.r.t. alpha beta grad.
C
C
C  NOTE: For non local functionals, the potentials (and
C        derivatives) returned by this subroutine are computed
C        with respect to the gradient invariant.
C        This is a difference with the PQS implementation of
C        the SCF and force calculation, where subroutines
C        vxcpot and vxcpotu compute potentials with respect to
C        the density gradient.
C
      DIMENSION DenA(NP),DenB(NP),DenXA(3,NP),DenXB(3,NP)
      DIMENSION WGHT(NP),EVec(NP)
      INTEGER dft
      dimension pra(np),prb(np),pga(np),pgb(np),pgc(np)
      dimension prara(np),prarb(np),prbrb(np)
      dimension praga(np),pragb(np),pragc(np)
      dimension prbga(np),prbgb(np),prbgc(np)
      dimension pgaga(np),pgagb(np),pgagc(np)
      dimension pgbgb(np),pgbgc(np),pgcgc(np)
C
      parameter (Zero=0.0d0,Third=1.0d0/3.0d0,Half=0.5d0,Four=4.0d0)
cc      parameter (PI=3.14159 26535 89793d0)
      parameter (two=2.0d0,TThird=2.0d0/3.0d0)
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
       pra(I) = Zero
       prb(I) = Zero
       prara(I) = Zero
       prbrb(I) = Zero
       prarb(I) = Zero
       If(dft.GT.3) Then
         pga(I) = Zero
         pgb(I) = Zero
         pgc(I) = Zero
         praga(I) = Zero
         pragb(I) = Zero
         pragc(I) = Zero
         prbga(I) = Zero
         prbgb(I) = Zero
         prbgc(I) = Zero
         pgaga(I) = Zero
         pgagb(I) = Zero
         pgagc(I) = Zero
         pgbgb(I) = Zero
         pgbgc(I) = Zero
         pgcgc(I) = Zero
       EndIf
       GO TO 20
      EndIf
C
C  Do the Self-Contained functionals first
C
      IF(dft.GT.23) THEN

        saa = DenXA(1,I)**2 + DenXA(2,I)**2 + DenXA(3,I)**2
        sbb = DenXB(1,I)**2 + DenXB(2,I)**2 + DenXB(3,I)**2
        ra=zero
        rb=zero
        rara=zero
        rbrb=zero
        rarb=zero
        ga = zero
        gb = zero
        gc = zero
        raga = zero
        ragb = zero
        ragc = zero
        rbga = zero
        rbgb = zero
        rbgc = zero
        gaga = zero
        gagb = zero
        gagc = zero
        gbgb = zero
        gbgc = zero
        gcgc = zero
c
       If(dft.EQ.24) Then
c...
c...  Becke 97.
c...
        call b97du(da,db,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $             drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $             dsaasaa,dsbbsbb,dsaasbb)
        exc=ec
        ra=dra
        rb=drb
        rara=drara
        rbrb=drbrb
        rarb=drarb
        ga = dsaa
        gb = dsbb
        raga = drasaa
        ragb = drasbb
        rbga = drbsaa
        rbgb = drbsbb
        gaga = dsaasaa
        gbgb = dsbbsbb
        gagb = dsaasbb
cc
       Else If(dft.EQ.25) Then
c...
c...  Becke 97_1.
c...
        call b971du(da,db,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $              drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $              dsaasaa,dsbbsbb,dsaasbb)
        exc=ec
        ra=dra
        rb=drb
        rara=drara
        rbrb=drbrb
        rarb=drarb
        ga = dsaa
        gb = dsbb
        raga = drasaa
        ragb = drasbb
        rbga = drbsaa
        rbgb = drbsbb
        gaga = dsaasaa
        gbgb = dsbbsbb
        gagb = dsaasbb
cc
       Else If(dft.EQ.26) Then
c...
c...  Becke 97_2.
c...
        call b972du(da,db,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $              drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $              dsaasaa,dsbbsbb,dsaasbb)
        exc=ec
        ra=dra
        rb=drb
        rara=drara
        rbrb=drbrb
        rarb=drarb
        ga = dsaa
        gb = dsbb
        raga = drasaa
        ragb = drasbb
        rbga = drbsaa
        rbgb = drbsbb
        gaga = dsaasaa
        gbgb = dsbbsbb
        gagb = dsaasbb
cc
       Else
c...
c...  HCTH.
c...
        call hcth407du(da,db,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                 drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                 dsaasaa,dsbbsbb,dsaasbb)
        exc=ec
        ra=dra
        rb=drb
        rara=drara
        rbrb=drbrb
        rarb=drarb
        ga = dsaa
        gb = dsbb
        raga = drasaa
        ragb = drasbb
        rbga = drbsaa
        rbgb = drbsbb
        gaga = dsaasaa
        gbgb = dsbbsbb
        gagb = dsaasbb
cc
       EndIf
      ELSE
C
C  get the local Slater exchange
C
        DA13 = DA**Third
        DB13 = DB**Third
        DA43 = DA*DA13
        DB43 = DB*DB13
        EXC = XS*EFS*(DA43+DB43)/(DA+DB)
        ra = XS*VFS*DA13
        rb = XS*VFS*DB13
        rara = Third*XS*VFS/(DA13*DA13)
        rbrb = Third*XS*VFS/(DB13*DB13)
        rarb=zero
        If(dft.EQ.1) GO TO 10
C
C  VWN local Correlation
C
        IF(XVWN.NE.Zero) THEN
          CALL VWND(DA,DB,.false.,ECVWN,dra,drb,drara,drbrb,drarb)
          EXC = EXC + ECVWN*XVWN
          ra = ra + dra*XVWN
          rb = rb + drb*XVWN
          rara= rara + drara*XVWN
          rbrb= rbrb + drbrb*XVWN
          rarb= rarb + drarb*XVWN
          If(dft.EQ.2) GO TO 10
        ELSE IF(XVWN5.NE.Zero) THEN
          CALL VWN5D(DA,DB,.false.,ECVWN,dra,drb,drara,drbrb,drarb)
          EXC = EXC + ECVWN*XVWN5
          ra = ra + dra*XVWN5
          rb = rb + drb*XVWN5
          rara= rara + drara*XVWN5
          rbrb= rbrb + drbrb*XVWN5
          rarb= rarb + drarb*XVWN5
          If(dft.EQ.3) GO TO 10
        ENDIF
c
        g2a = (DenXA(1,I)**2 + DenXA(2,I)**2 + DenXA(3,I)**2)
        g1a = SQRT(g2a)
        g2b = (DenXB(1,I)**2 + DenXB(2,I)**2 + DenXB(3,I)**2)
        g1b = SQRT(g2b)
        g2ab = DenXA(1,I)*DenXB(1,I) +DenXA(2,I)*DenXB(2,I) +
     $         DenXA(3,I)*DenXB(3,I)
        ga = zero
        gb = zero
        gc = zero
        raga = zero
        ragb = zero
        ragc = zero
        rbga = zero
        rbgb = zero
        rbgc = zero
        gaga = zero
        gagb = zero
        gagc = zero
        gbgb = zero
        gbgc = zero
        gcgc = zero
        IF(XB88.NE.Zero) THEN
C
C  Becke non-local exchange
C
          Call Becke88D(DA,     DA13,   DA43,   g1a,   g2a,
     $                  DB,     DB13,   DB43,   g1b,   g2b,
     $    .false.,EXB88,dra,    dga,    drara,  dgaga, draga,
     $                  drb,    dgb,    drbrb,  dgbgb, drbgb)
          EXC = EXC + EXB88*XB88
          ra = ra + dra*XB88
          rb = rb + drb*XB88
          ga = ga + dga*XB88
          gb = gb + dgb*XB88
          rara = rara + drara*XB88
          rbrb = rbrb + drbrb*XB88
          gaga = gaga + dgaga*XB88
          gbgb = gbgb + dgbgb*XB88
          raga = raga + draga*XB88
          rbgb = rbgb + drbgb*XB88
        ENDIF
C
C  Handy/Cohen optimized non-local exchange
C
        IF(XOptX.NE.Zero) THEN
          Call    OPTXD(DA,     DA13,   DA43,   g1a,   g2a,
     $                  DB,     DB13,   DB43,   g1b,   g2b,
     $   .false.,EXOPTX,dra,    dga,    drara,  dgaga, draga,
     $                  drb,    dgb,    drbrb,  dgbgb, drbgb)
          EXC = EXC + EXOPTX*XOptX
          ra = ra + dra*XOptX
          rb = rb + drb*XOptX
          ga = ga + dga*XOptX
          gb = gb + dgb*XOptX
          rara = rara + drara*XOptX
          rbrb = rbrb + drbrb*XOptX
          gaga = gaga + dgaga*XOptX
          gbgb = gbgb + dgbgb*XOptX
          raga = raga + draga*XOptX
          rbgb = rbgb + drbgb*XOptX
        ENDIF
c
c  Perdew 86 non-local correlation
c
       IF(XP86L.NE.Zero.OR.XP86NL.NE.Zero) THEN
         call p86nlcdu(da,db,g2a,g2b,g2ab,xp86l,xp86nl,ecp86,
     $                    dra,drb,dsaa,dsbb,dsab,
     $                    drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                    drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $                    dsabsab,dsaasbb,dsaasab,dsbbsab)
         EXC = EXC + ECP86
          ra=ra+dra
          rb=rb+drb
          ga=ga+dsaa
          gb=gb+dsbb
          gc=gc+dsab
          rara=rara+drara
          rbrb=rbrb+drbrb
          rarb=rarb+drarb
          raga=raga+drasaa
          ragb=ragb+drasbb
          ragc=ragc+drasab
          rbga=rbga+drbsaa
          rbgb=rbgb+drbsbb
          rbgc=rbgc+drbsab
          gaga=gaga+dsaasaa
          gagb=gagb+dsaasbb
          gagc=gagc+dsaasab
          gbgb=gbgb+dsbbsbb
          gbgc=gbgc+dsbbsab
          gcgc=gcgc+dsabsab
       ENDIF
C
C  Perdew-Wang 91 non-local exchange
C
       IF(XPW91X.NE.Zero) THEN
        call pw91xdu(da,db,g2a,g2b,expw91,dra,drb,dsaa,dsbb,drara,drbrb,
     $               drasaa,drbsbb,dsaasaa,dsbbsbb)
         EXC = EXC + expw91*xpw91x
         ra = ra + dra*xpw91x
         rb = rb + drb*xpw91x
         ga = ga + dsaa*xpw91x
         gb = gb + dsbb*xpw91x
         rara = rara + drara*xpw91x
         rbrb = rbrb + drbrb*xpw91x
         raga = raga + drasaa*xpw91x
         rbgb = rbgb + drbsbb*xpw91x
         gaga = gaga + dsaasaa*xpw91x
         gbgb = gbgb + dsbbsbb*xpw91x
       ENDIF
C
C  Perdew-Wang non-local correlation
C
        IF(XPW91L.NE.Zero.OR.XPW91NL.NE.Zero) THEN
          call pw91nlcdu(da,db,g2a,g2b,g2ab,XPW91L,XPW91NL,ecpw91,
     $                         dra,drb,dsaa,dsbb,dsab,
     $                         drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                         drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $                         dsabsab,dsaasbb,dsaasab,dsbbsab)
          EXC = EXC + ecpw91
          ra=ra+dra
          rb=rb+drb
          ga=ga+dsaa
          gb=gb+dsbb
          gc=gc+dsab
          rara=rara+drara
          rbrb=rbrb+drbrb
          rarb=rarb+drarb
          raga=raga+drasaa
          ragb=ragb+drasbb
          ragc=ragc+drasab
          rbga=rbga+drbsaa
          rbgb=rbgb+drbsbb
          rbgc=rbgc+drbsab
          gaga=gaga+dsaasaa
          gagb=gagb+dsaasbb
          gagc=gagc+dsaasab
          gbgb=gbgb+dsbbsbb
          gbgc=gbgc+dsbbsab
          gcgc=gcgc+dsabsab
        ENDIF
C
C  Lee, Yang & Parr non-local correlation
C
C   this functional is linear on the density gradient, thus
C   its second derivatives with respect to the density gradients
C   two times are zero
C
        IF(XLYP.NE.Zero) THEN
          call lypdu(da,db,g2a,g2b,g2ab,ec,dra,drb,dsaa,dsbb,dsab,drara
     $         ,drbrb,drarb,drasaa,drasbb,drasab,drbsaa,drbsbb,drbsab)
          EXC = EXC + EC*XLYP
          ra=ra+dra*XLYP
          rb=rb+drb*XLYP
          ga=ga+dsaa*XLYP
          gb=gb+dsbb*XLYP
          gc=gc+dsab*XLYP
          rara=rara+drara*XLYP
          rbrb=rbrb+drbrb*XLYP
          rarb=rarb+drarb*XLYP
          raga=raga+drasaa*XLYP
          ragb=ragb+drasbb*XLYP
          ragc=ragc+drasab*XLYP
          rbga=rbga+drbsaa*XLYP
          rbgb=rbgb+drbsbb*XLYP
          rbgc=rbgc+drbsab*XLYP
        ENDIF
C
C  Perdew, Burke and Ernzerhof non-local exchange
C
        IF(XPBEX.NE.Zero) THEN
          call pbexdu(da,db,g2a,g2b,excpbe,dra,drb,dsaa,dsbb,
     $                drara,drbrb,drasaa,drbsbb,dsaasaa,dsbbsbb)
          EXC=EXC+excpbe*xpbex
          ra = ra + dra*xpbex
          rb = rb + drb*xpbex
          ga = ga + dsaa*xpbex
          gb = gb + dsbb*xpbex
          rara = rara + drara*xpbex
          rbrb = rbrb + drbrb*xpbex
          raga = raga + drasaa*xpbex
          rbgb = rbgb + drbsbb*xpbex
          gaga = gaga + dsaasaa*xpbex
          gbgb = gbgb + dsbbsbb*xpbex
        ENDIF
C
C  Perdew, Burke and Ernzerhof non-local correlation
C
        IF(XPBEC.NE.Zero) THEN
          call pbecdu(da,db,g2a,g2b,g2ab,ecpbe,dra,drb,dsaa,dsbb,
     $                dsab,drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,dsabsab,
     $                dsaasbb,dsaasab,dsbbsab)
          EXC=EXC+ecpbe*xpbec
          ra=ra+dra*xpbec
          rb=rb+drb*xpbec
          ga=ga+dsaa*xpbec
          gb=gb+dsbb*xpbec
          gc=gc+dsab*xpbec
          rara=rara+drara*xpbec
          rbrb=rbrb+drbrb*xpbec
          rarb=rarb+drarb*xpbec
          raga=raga+drasaa*xpbec
          ragb=ragb+drasbb*xpbec
          ragc=ragc+drasab*xpbec
          rbga=rbga+drbsaa*xpbec
          rbgb=rbgb+drbsbb*xpbec
          rbgc=rbgc+drbsab*xpbec
          gaga=gaga+dsaasaa*xpbec
          gagb=gagb+dsaasbb*xpbec
          gagc=gagc+dsaasab*xpbec
          gbgb=gbgb+dsbbsbb*xpbec
          gbgc=gbgc+dsbbsab*xpbec
          gcgc=gcgc+dsabsab*xpbec
        ENDIF
C
C ..................................................
C  other functional density derivative code here
C ..................................................
C
C
      ENDIF
      pga(i) = ga
      pgb(i) = gb
      pgc(i) = gc
      praga(I) = raga
      pragb(I) = ragb
      pragc(I) = ragc
      prbga(I) = rbga
      prbgb(I) = rbgb
      prbgc(I) = rbgc
      pgaga(I) = gaga
      pgagb(I) = gagb
      pgagc(I) = gagc
      pgbgb(I) = gbgb
      pgbgc(I) = gbgc
      pgcgc(I) = gcgc
C
 10   CONTINUE
C
C  Accumulate into total Exchange-Correlation energy
C  and total number of electrons
C
      DVal = (DA+DB)*WGHT(I)
      ETot = ETot + EXC*DVal
      EL   = EL + DVal
C
C  store full potential and exchange-correlation energy
C
      pra(I) = ra
      prb(I) = rb
      prara(I) = rara
      prbrb(I) = rbrb
      prarb(I) = rarb
      if(IdWt.eq.1) EVec(I) = EXC*(DA+DB) ! for weight derivatives
 20   CONTINUE
C
      RETURN
      END

c =====================================================================
      SUBROUTINE VWN5D(DenA,DenB,rhf,EC,PotA,PotB,PotDA,PotDB,
     +                 PotDAB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     MM (05/13/2003)
C
C  Vosko, Wilk & Nusair (VWN)  Local correlation
C  this is the parametrization number 5 in the Vosko paper
C  (fit of the Cepeley-Alder monte carlo values of correlation
C  energy, using equations [3.2] and [4.8] to define Delta epsilon_c)
C  The subroutine is coded according to the expressions given in
C  Poples' paper.
C
C  reference
C  S.H.Vosko, L.Wilk and M.Nusair,  Can J. Phys. 58 (1980) 1200
C  B.G.Johnson, P.M.W.Gill, ans J.A.Pople, J. Chem. Phys. 98 (1993) 5612
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
C  PotDA   -  derivative of the alpha (or closed shell) potential
C  PotDB   -  derivative of the beta potential
C  PotDAB  -  cross derivative for open-shell case
C
C
      REAL*8 A(3),B(3),C(3),X0(3),Q(3)
      LOGICAL rhf
C
      parameter (zero=0.0d0,epsi=1.0d-15)
      parameter (Half=0.5d0,One=1.0d0,Two=2.0d0,Three=3.0d+0,Four=4.0d0)
      parameter (Five=5.0d+0,Eight=8.0d+0,Nine=9.0d+0)
      parameter (Twelve=12.0d+0,Onesix=16.0d+0)
      parameter (Third=One/3.0d0,Sixth=One/6.0d0)
      parameter (PI=3.14159 26535 89793d0)
c
c...  aa= 2*2**(1/3)
c...  ab= 2**(1/3)-1
c...  f2= 4/(9*(2**(1/3)-1))
c...
      PARAMETER (aa=2.519842099789746d0,ab=0.2599210498948732d0,
     $           f2=1.709920934161365d0)
C
c...
c...  in the notation of Poples' paper, the first element of
c...  the following constants corresponds to superscript A,
c...  the second to P and the third to F
c...
c...  the set of data is for parametrization 5
c...  (Table 5 and Section 4 of Volsko paper)
c...
      DATA A /-0.0337737d0,  0.0621814d0, 0.0310907d0/
      DATA B / 1.13107d0,    3.72744d0,   7.06042d0/
      DATA C /13.00450d0,   12.93520d0,  18.05780d0/
      DATA X0/-0.0047584d0, -0.10498d0,  -0.32500d0/
      DATA Q/7.12310891781812d0,6.15199081975908d0,4.73092690956011d0/
cc
cc
c  function statements
c...
c...  function g(z), equation (A11)
c...
      g(z) = nine/eight*((One+z)*(One+z)**Third+(One-z)*(One-z)**Third-
     +       Two)
c...
c...   g'(z), (A20)
c...
      gp(z) = three/two*((One+z)**Third-(One-z)**Third)
c...
c...   g''(z)
c...
      gpp(z)=half*(one/((one+z)**Third*(one+z)**Third)+
     +             one/((one-z)**Third*(one-z)**Third))
c...
c...  rs corresponds to x**2 (see eq. A10)
c...
      rs(z) = (0.75d0/(PI*z))**Third
c...
c...  capital X(x), (A14)
c...
      cx(x,m) = x**2 + b(m)*x + c(m)
c...
c...  epsilonc(x), (A13)
c...
      e0(x,m)=two*x+b(m)
      e1(m)=one/Q(m)
      e2(x,m)=ATAN(Q(m)/e0(x,m))
      e3(x,m)=one/cx(x,m)
      epsc(x,m)= Half*A(m)*(-b(m)*x0(m)*(two*e1(m)*e2(x,m)*e0(x0(m),m)
     +   + LOG(e3(x,m)*(x-x0(m))**2))/cx(x0(m),m)+LOG(x**2*e3(x,m))+
     +     two*b(m)*e1(m)*e2(x,m))
c...
c...  h(x), (A15)
c...
      h(x)=f2*(epsc(x,3)-epsc(x,2))/epsc(x,1)-one
c...
c...  epsilonc'(x), (A18)
c...
      e4(x,m)=one/(e0(x,m)**2+Q(m)**2)
     +
      epscp(x,m)=half*A(m)*((-b(m)*x0(m)*(two/(x-x0(m))-four*e4(x,m)*
     +     e0(x0(m),m))+b(m)*e0(x,m)*e3(x,m)*x0(m))*e3(x0(m),m)-
     +     four*b(m)*e4(x,m)-e0(x,m)*e3(x,m)+two/x)
c...
c...  h'(x), (A19)
c...
      hp(x)=(f2*(epscp(x,3)-epscp(x,2))-epscp(x,1)*(h(x)+one))/epsc(x,1)
c...
c...  epsilonc''(x)
c...
      epscpp(x,m)=half*a(m)*(onesix*b(m)*e4(x,m)*(one-e4(x,m)*Q(m)**2)*
     +    (one-x0(m)*e0(x0(m),m)*e3(x0(m),m))/e0(x,m)+(two-e0(x,m)**2*
     +    e3(x,m))*(b(m)*x0(m)*e3(x0(m),m)-one)*e3(x,m)+two*b(m)*
     +    x0(m)/(x-x0(m))**2*e3(x0(m),m)-two/x**2)
c...
c...  h''(x)
c...
      hpp(x)=(f2*(epscpp(x,3)-epscpp(x,2))-two*epscp(x,1)*hp(x)-
     +       epscpp(x,1)*(h(x)+one))/epsc(x,1)
c...
c...  end of function statements
c...
      IF(rhf) THEN
C
C  ** Closed Shell **
C
        rs1 = rs(DenA)
        x = sqrt(rs1)
cc  energy
        EC = epsc(x,2)
cc  potential
        ecpp = epscp(x,2)
        PotA = EC - x*sixth*ecpp
c...
c...  derivative of the potential
c...  Note: This is multiplied by a factor two
c...
        PotDA=Third*sixth*x/DenA*(x*epscpp(x,2)-five*ecpp)
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
c
        eca=epsc(x,1)
        ecap=epscp(x,1)
        ecp=epsc(x,2)
        ecpp=epscp(x,2)
        z = min( max ((DDif/DVal),-One),One)
        z4 = z**4
        z3 = z**3
        gz = g(z)
        gpz= gp(z)
        hx=h(x)
        hpx=hp(x)
        z4h1=z4*hx+one
c...
        gu1=ecpp+eca*z4*gz*hpx+ecap*gz*z4h1
        guu2=gpz*z4h1+four*z3*gz*hx
        gu2=eca*guu2
        sx=sixth*x
        omz=one-z
        opz=one+z
cc  energy
        EC = ecp+eca*gz*z4h1
cc  potential
        PotA = EC - sx*gu1 + gu2*omz
        PotB = EC - sx*gu1 - gu2*opz
c...
c...  derivativatives of the potential
c...
        ecapp=epscpp(x,1)
        ecppp=epscpp(x,2)
        hppx=hpp(x)
        if(omz.gt.epsi.and.opz.gt.epsi)then
          gppz= gpp(z)
        else
          gppz=zero
        endif
        ggu2=ecppp+eca*z4*gz*hppx+two*ecap*z4*gz*hpx+ecapp*gz*z4h1
        ggu3=eca*(two*z4*gpz+eight*z3*gz)*hpx+ecap*two*guu2
        ggu4=eca*(gppz*z4h1+eight*z3*gpz*hx+twelve*z*z*gz*hx)
c...
        PotDA=(sx*(-five*sixth*gu1+sx*ggu2-omz*ggu3)+omz**2*ggu4)/DVal
        PotDB=(sx*(-five*sixth*gu1+sx*ggu2+opz*ggu3)+opz**2*ggu4)/DVal
        PotDAB=(sx*(-five*sixth*gu1+sx*ggu2+z*ggu3)-(one-z*z)*ggu4)/Dval
C
      ENDIF
C
      RETURN
      END
c =====================================================================
      SUBROUTINE VWND(DenA,DenB,rhf,EC,PotA,PotB,PotDA,PotDB,
     +                 PotDAB)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     MM (05/13/2003)
C
C  Vosko, Wilk & Nusair (VWN)  Local correlation.
C  This is the fit of the RPA correlation energy,
C  using expression [3.1] of Volsko paper for Delta epsilon_c.
C  That is, the difference with VWN5 is not only in the parameters,
C  but also in the form of the functional.
C  This is the functional used by Gaussian with the LDA keyword.
C
C  reference
C  S.H.Vosko, L.Wilk and M.Nusair,  Can J. Phys. 58 (1980) 1200
C  B.G.Johnson, P.M.W.Gill, ans J.A.Pople, J. Chem. Phys. 98 (1993) 5612
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
C  PotDA   -  derivative of the alpha (or closed shell) potential
C  PotDB   -  derivative of the beta potential
C  PotDAB  -  cross derivative for open-shell case
C
C
      REAL*8 A(3),B(3),C(3),X0(3),Q(3)
      LOGICAL rhf
C
      parameter (zero=0.0d0,epsi=1.0d-15)
      parameter (Half=0.5d0,One=1.0d0,Two=2.0d0,Three=3.0d+0,Four=4.0d0)
      parameter (Five=5.0d+0,Eight=8.0d+0,Nine=9.0d+0)
      parameter (Twelve=12.0d+0,Onesix=16.0d+0)
      parameter (Third=One/3.0d0,Sixth=One/6.0d0)
      parameter (PI=3.14159 26535 89793d0)
c
c...  aa= 2*2**(1/3)
c...  ab= 2**(1/3)-1
c...  f2= 4/(9*(2**(1/3)-1))
c...
      PARAMETER (aa=2.519842099789746d0,ab=0.2599210498948732d0,
     $           f2=1.709920934161365d0)
C
c...
c...  In the notation of Poples' paper, the first element of
c...  the following constants corresponds to superscript A,
c...  the second to P and the third to F
c...
c...  this set of data is for the fit of the RPA correlation energy
c...  (Section 4 of Volsko paper). Only columns 2 and 3 are actually
c...  used in the calculation, due to the form of delta epsilon_c
c...  (equation [3.1] of Volsko paper).
c...
      DATA A /-0.0337737d0,  0.0621814d0, 0.0310907d0/
      DATA B / 1.06835d0,    13.0720d0,  20.1231d0/
      DATA C /11.4813d0,   42.7198d0,  101.578d0/
      DATA X0/-0.228344d0, -0.409286d0,  -0.743294d0/
      DATA Q/6.692072046645942d0,0.0448998886415768d0,
     +       1.171685277708971d0/
cc
cc
c  function statements
c...
c...  function g(z)
c...
      g(z) =((One+z)*(One+z)**Third+(One-z)*(One-z)**Third-Two)/two/ab
c...
c...   g'(z)
c...
      gp(z) = two*third/ab*((One+z)**Third-(One-z)**Third)
c...
c...   g''(z)
c...
      gpp(z)=two/(nine*ab)*(one/((one+z)**Third*(one+z)**Third)+
     +             one/((one-z)**Third*(one-z)**Third))
c...
c...  rs corresponds to x**2 (see eq. A10)
c...
      rs(z) = (0.75d0/(PI*z))**Third
c...
c...  capital X(x), (A14)
c...
      cx(x,m) = x**2 + b(m)*x + c(m)
c...
c...  epsilonc(x), (A13)
c...
      e0(x,m)=two*x+b(m)
      e1(m)=one/Q(m)
      e2(x,m)=ATAN(Q(m)/e0(x,m))
      e3(x,m)=one/cx(x,m)
      epsc(x,m)= Half*A(m)*(-b(m)*x0(m)*(two*e1(m)*e2(x,m)*e0(x0(m),m)
     +   + LOG(e3(x,m)*(x-x0(m))**2))/cx(x0(m),m)+LOG(x**2*e3(x,m))+
     +     two*b(m)*e1(m)*e2(x,m))
c...
c...  epsilonc'(x), (A18)
c...
      e4(x,m)=one/(e0(x,m)**2+Q(m)**2)
     +
      epscp(x,m)=half*A(m)*((-b(m)*x0(m)*(two/(x-x0(m))-four*e4(x,m)*
     +     e0(x0(m),m))+b(m)*e0(x,m)*e3(x,m)*x0(m))*e3(x0(m),m)-
     +     four*b(m)*e4(x,m)-e0(x,m)*e3(x,m)+two/x)
c...
c...  epsilonc''(x)
c...
      epscpp(x,m)=half*a(m)*(onesix*b(m)*e4(x,m)*(one-e4(x,m)*Q(m)**2)*
     +    (one-x0(m)*e0(x0(m),m)*e3(x0(m),m))/e0(x,m)+(two-e0(x,m)**2*
     +    e3(x,m))*(b(m)*x0(m)*e3(x0(m),m)-one)*e3(x,m)+two*b(m)*
     +    x0(m)/(x-x0(m))**2*e3(x0(m),m)-two/x**2)
c...
c...  end of function statements
c...
      IF(rhf) THEN
C
C  ** Closed Shell **
C
        rs1 = rs(DenA)
        x = sqrt(rs1)
cc  energy
        EC = epsc(x,2)
cc  potential
        ecpp = epscp(x,2)
        PotA = EC - x*sixth*ecpp
c...
c...  derivative of the potential
c...  Note: This is multiplied by a factor two
c...
        PotDA=Third*sixth*x/DenA*(x*epscpp(x,2)-five*ecpp)
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
c
        ecp=epsc(x,2)
        ecpp=epscp(x,2)
        ecf=epsc(x,3)
        ecfp=epscp(x,3)
        z = min( max ((DDif/DVal),-One),One)
        gz = g(z)
        gpz = gp(z)
c...
        dec=ecf-ecp
        decp=ecfp-ecpp
        sx=sixth*x
        gu1=sx*(ecpp+gz*decp)
        gu2=gpz*dec
        omz=one-z
        opz=one+z
cc  energy
        EC = ecp+gz*dec
cc  potential
        PotA = EC - gu1 + gu2*omz
        PotB = EC - gu1 - gu2*opz
c...
c...  derivativatives of the potential
c...
        ecppp=epscpp(x,2)
        ecfpp=epscpp(x,3)
        if(omz.gt.epsi.and.opz.gt.epsi)then
          gppz= gpp(z)
        else
          gppz=zero
        endif
        ggu1=sixth*sixth*(gz*((ecfpp-ecppp)*rs1-five*decp*x)+ecppp*rs1-
     +        five*ecpp*x)
        ggu2=third*decp*gpz*x
        ggu3=gppz*dec
c...
        PotDA=(ggu1-ggu2*omz+ggu3*omz*omz)/DVal
        PotDB=(ggu1+ggu2*opz+ggu3*opz*opz)/Dval
        PotDAB=(ggu1+ggu2*z-ggu3*opz*omz)/DVal
C
      ENDIF
C
      RETURN
      END
c=======================================================================

      SUBROUTINE Becke88D(DA,    DA13,   DA43,   g1A,    gama,
     $                    DB,    DB13,   DB43,   g1B,    gamb,
     $     rhf,   EXB88,  dfdda, dfdga,  dfdda2, dfdga2, dfddaga,
     $                    dfddb, dfdgb,  dfddb2, dfdgb2, dfddbgb)
      IMPLICIT REAL*8(A-H,O-Z)
C
C   MM (02/13/2003)
C
C  Becke 1988 nonlocal exchange
C  gradient-corrected exchange term
C
C  see  B.G. Johnson, P.M.W. Gill and J.A. Pople,
C       J. Chem. Phys. 98, 5612 (1993) eqs (A2)--(A8).
C
C    !!! attention to the misprints: in the definition of the
C        auxiliary variable x (eq. A5), gamma should actually be
C        elevated at 1/2. Note also that the B88 definition
C        includes a Slater Exchange term, which has
C        been subtracted from the functional here defined.
C
C  this subroutine computes the functional and its first and
C  second derivatives for both closed and open shell
C
C  ARGUMENTS
C
C  DA      -  alpha density at current grid point
C  DA13    -  DA**(1/3)
C  DA43    -  DA**(4/3)
C  g1A     -  absolute value of alpha density derivative
C  gama    -  square of g1a
C  DB      -  beta density at current grid point
C  DB13    -  DB**(1/3)
C  DB43    -  DB**(4/3)
C  g1B     -  absolute value of beta density derivative
C  gamb    -  square of g1b
C  rhf     -  logical flag for rhf/uhf
C
C  on exit
C
C  EXB88   -  contribution to the exchange energy, that is,
C             the functional values divided by the total density
C             at the current grid point.
C  dfdda   -  functional derivative with respect the alpha density
C  dfdga   -  derivative with respect the alpha density gradient
C  dfdda2  -  second derivative with respect the alpha density
C  dfdga2  -  second derivative with respect the alpha gradient
C  dfddaga -  second derivative with respect density and gradient
C
C  dfddb   -  functional derivative with respect the beta density
C  dfdgb   -  derivative with respect the beta density gradient
C  dfddb2  -  second derivative with respect the beta density
C  dfdgb2  -  second derivative with respect the beta gradient
C  dfddbgb -  second derivative with respect density and gradient
C
      logical rhf
      parameter(Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0,Four=4.0d0)
      parameter(TThird=2.0d0/3.0d0,FThird=4.0d0/3.0d0,thrsh=1.0d-12)
      parameter(Fourth=1.0d0/4.0d0,FNine=4.0d0/9.0d0)
c
      PARAMETER (b=0.0042d0,b6=0.0252d0)
      PARAMETER (b6b=b6*b)
c
c...
      if(rhf)then
c...
c...  Closed Shell
c...
        xa = g1a/da43
        x21a = xa**2+One
        xs21a = SQRT(x21a)
        xl21a = LOG(xs21a+xa)
        dgfa = b6*xa*xl21a+One
        gfa = -b*xa**2/dgfa
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
        exb88 = da13*gfa
c...
c...  First derivatives of the functional (potential)
c...
        dgfpa = b6*(xl21a*xs21a+xa)/xs21a
        gfpa = b*xa*(dgfpa*xa-Two*dgfa)/dgfa**2
        dfdda = -FThird*(g1a*gfpa-da43*gfa)/da
        if(g1a.gt.thrsh)then
          dfdga = Half*gfpa/g1a
        else
          dfdga = ZERO
        endif
c...
c...  Second derivatives of the functional
c...
        dgfsa = -b6*(xa**2-Two*x21a)/xs21a**3
        gfsa = b*(-Two*dgfpa**2*xa**2+dgfa*xa*(dgfsa*xa+Four*dgfpa)-
     1     Two*dgfa**2)/dgfa**3
        dfdda2 = FNine*(Four*da13*gama*gfsa-
     1     da43*da13*g1a*gfpa+da**3*gfa)/(da43*da43*da)
        if(g1a.gt.thrsh)then
          dfdga2 = Fourth*(g1a*gfsa-da43*gfpa)/(da43*g1a**3)
        else
          dfdga2 = ZERO
        endif
        dfddaga = -TThird*gfsa/(da43*da)
c...
      else
c...
c...  Open Shell - Alpha Density
c...
        if(da.gt.thrsh)then
          xa = g1a/da43
          x21a = xa**2+One
          xs21a = SQRT(x21a)
          xl21a = LOG(xs21a+xa)
          dgfa = b6*xa*xl21a+One
          gfa = -b*xa**2/dgfa
c...
c...  First derivatives of the functional (potential)
c...
          dgfpa = b6*(xl21a*xs21a+xa)/xs21a
          gfpa = b*xa*(dgfpa*xa-Two*dgfa)/dgfa**2
          dfdda = -FThird*(g1a*gfpa-da43*gfa)/da
          if(g1a.gt.thrsh)then
            dfdga = Half*gfpa/g1a
          else
            dfdga = ZERO
          endif
c...
c...  Second derivatives of the functional
c...
          dgfsa = -b6*(xa**2-Two*x21a)/xs21a**3
          gfsa = b*(-Two*dgfpa**2*xa**2+dgfa*xa*(dgfsa*xa+Four*dgfpa)-
     1       Two*dgfa**2)/dgfa**3
          dfdda2 = FNine*(Four*da13*gama*gfsa-
     1       da43*da13*g1a*gfpa+da**3*gfa)/(da43*da43*da)
          if(g1a.gt.thrsh)then
            dfdga2 = Fourth*(g1a*gfsa-da43*gfpa)/(da43*g1a**3)
          else
            dfdga2 = ZERO
          endif
          dfddaga = -TThird*gfsa/(da43*da)
        else
          gfa= ZERO
          dfdda = ZERO
          dfdda2 = ZERO
          dfdga = ZERO
          dfdga2 = ZERO
          dfddaga = ZERO
        endif
c...
c...  Open Shell - Beta Density
c...
        if(db.gt.thrsh)then
          xb = g1b/db43
          x21b = xb**2+One
          xs21b = SQRT(x21b)
          xl21b = LOG(xs21b+xb)
          dgfb = b6*xb*xl21b+One
          gfb = -b*xb**2/dgfb
c...
c...  First derivatives of the functional (potential)
c...
          dgfpb = b6*(xl21b*xs21b+xb)/xs21b
          gfpb = b*xb*(dgfpb*xb-Two*dgfb)/dgfb**2
          dfddb = -Fthird*(g1b*gfpb-db43*gfb)/db
          if(g1b.gt.thrsh)then
            dfdgb = Half*gfpb/g1b
          else
            dfdgb = ZERO
          endif
c...
c...  Second derivatives of the functional
c...
          dgfsb = -b6*(xb**2-Two*x21b)/xs21b**3
          gfsb = b*(-Two*dgfpb**2*xb**2+dgfb*xb*(dgfsb*xb+Four*dgfpb)-
     1       Two*dgfb**2)/dgfb**3
          dfddb2 = Fnine*(Four*db13*gamb*gfsb-
     1       db43*db13*g1b*gfpb+db**3*gfb)/(db43*db43*db)
          if(g1b.gt.thrsh)then
            dfdgb2 = Fourth*(g1b*gfsb-db43*gfpb)/(db43*g1b**3)
          else
            dfdgb2 = ZERO
          endif
          dfddbgb = -TThird*gfsb/(db43*db)
        else
          gfb = ZERO
          dfddb = ZERO
          dfddb2 = ZERO
          dfdgb = ZERO
          dfdgb2 = ZERO
          dfddbgb = ZERO
        endif
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
        exb88 = (db43*gfb+da43*gfa)/(db+da)
      endif
c...
      return
      end
c=====================================================================
      SUBROUTINE LYPCD(DENA,DENA13,DENA43,GAMA,ECLYP,DRA,DGAA,DGAB,
     $                 DRA2,DRAB,DRAGAA,DRAGBB,DRAGAB)
      IMPLICIT NONE    ! trying to avoid typing errors
C
C   MM (04/23/2003)  (partly generated with maxima)
C
C  Lee-Yang-Parr non-local correlation term, as
C  reformulated by A.Savin
C  ** CLOSED SHELL **
C
C  reference
C  C.Lee, W.Yang and R.G.Parr,  Phys.Rev. B37 (1988) 785
C
C  B.Miehlich, A.Savin, H.Stoll and H.Preuss, Chem.Phys.Lett.
C
C  B.G. Johnson, P.M.W. Gill and J.A. Pople,
C  J. Chem. Phys. 98, 5612 (1993) eqs (A22)--(A33).
C
C  this subroutine computes the functional and its first and
C  second derivatives. The LYP functional is linear in the
C  density gradients, thus the second derivatives with respect
C  the gradients two times are zero.
C
C  ARGUMENTS
C
C  DENA    Half of the total closed shell density at the current
C          grid point
C  DENA13  DENA**(1/3)
C  DENA43  DENA**(4/3)
C  GAMA    Density gradient invariant at the current grid point
C
C  on exit
C
C  ECLYP   exchange energy contribution at current grid point
C  DRA functional derivative with respect to the density
C  DGAA functional derivative w.r.t. the alpha gradient invariant
C  DGAB functional derivative w.r.t. the alpha beta gradient invariant
C  DRA2 second derivative w.r.t. the alpha density
C  DRAB mixed second deriv. w.r.t. the alpha and beta densities
C  DRAGAA mixed second derivative w.r.t. alpha density and gradient
C  DRAGBB mixed second deriv. w.r.t. alpha density and beta gradient
C  DRAGAB mixed second deriv. w.r.t. alpha dens. and alpha beta grad.
C
      REAL*8 DENA13,DENA43,DENA,GAMA,eclyp,dra,dgaa,dgab,
     $       dra2,drab,dragaa,dragbb,dragab
      REAL*8 z,zm1,ze
      REAL*8 w,wp,ws,delt,deltp,delts,sdelt,sdelt1,sdelt5
      REAL*8 eclyp1,eclyp1pa,eclyp1saa,eclyp1sab
      REAL*8 eclyp2,eclyp2pa,eclyp2saa,eclyp2sab
      REAL*8 pgaa,pgaapa,pgaasaa,pgaasab
      REAL*8 pgab,pgabpa,pgabsaa,pgabsab
      REAL*8 pgbbpa,pgbbsaa
      REAL*8 DENA2,DENA83,DENA113
      REAL*8 DEN,DEN2,DENM1,DEN13,DENM13,DENM53,DENM83,DENM113,DENM173
c
      REAL*8 A,B,C,D
      parameter (a=0.04918d0,b=0.132d0,c=0.2533d0,d=0.349d0)
c
      REAL*8 CF  !  CF= 2**(2/3)*(3/10)*(3*pi**2)**(2/3)
      parameter (CF=4.5577998723455959d0)
c
      REAL*8 one,two,three
      parameter (one=1.0d0,two=2.0d0,three=3.0d0)
      REAL*8 third,tt
      parameter (third=0.333333333333333333333d0)
      parameter (tt=0.66666666666666666667d0)
      REAL*8 on,tn,fn,stt,stn
      parameter (on=one/9.0d0,tn=two/9.0d0,fn=4.0d0/9.0d0)
      parameter (stt=16.0d0/3.0d0,stn=16.0d0/9.0d0)
      REAL*8 oet,ots
      parameter (oet=one/18.0d0,ots=one/36.0d0)
      real*8 zero,epsi
      parameter (zero=0.0d00,epsi=1.0d-15)
c...
c...  First  a check to make sure that the density is
c...  significant. This should always be passed on
c...  normal runs, as the density is checked by the
c...  calling routine, but sometimes I like to reduce
c...  the thresholds, and in that case we could run into
c...  trouble when we divide by the density.
c...
      if(dena.lt.epsi)then
         eclyp=zero
         dra=zero
         dgaa=zero
         dgab=zero
         dra2=zero
         drab=zero
         dragaa=zero
         dragbb=zero
         dragab=zero
         return
      endif
c...
c...  generates some additional powers of the density
c...
      DENA2=DENA*DENA
      DENA83=DENA43*DENA43
      DENA113=DENA83*DENA
      DEN = two*DENA
      DEN2=DEN*DEN
      DENM1=one/DEN
      DEN13=DEN**(third)
      DENM13=one/DEN13
      DENM53=DENM1*DENM13*DENM13
      DENM83=DENM53*DENM1
      DENM113=DENM83*DENM1
      DENM173=DENM113*DENM1*DENM1
c...
      z = d*DENM13+one
      zm1 = one/z
      ze = EXP(-C*DENM13)
c...
      w = DENM113*ze*zm1
      delt = DENM13*(d*zm1+C)
      sdelt= 7.0d0*delt
      sdelt1=sdelt+one
      sdelt5=sdelt+5.0d0
c...
      eclyp1 = -two*a*DENA*zm1
      eclyp2 = -16.0d0*a*b*CF*DENA113*DENA*w
      pgaa = a*b*sdelt5*DENA2*w*oet
      pgab = a*b*sdelt1*DENA2*w*on
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
      eclyp = (GAMA*(pgab+pgaa+pgaa)+eclyp2+eclyp1)*DENM1
c...
c...  First derivatives of the functional
c...
      wp = -DENM173*ze*zm1*(DEN*(zm1+10.0d0)-C*DEN13*DEN13)*third
      deltp = (d*d*DENM53*zm1*zm1-delt*DENM1)*third
      eclyp1pa = a*(zm1-4.0d0)*zm1*third
      eclyp2pa = -stt*a*b*CF*DENA113*(three*DENA*wp+7.0d0*w)
      pgaapa = a*b*DENA*(two*sdelt5*DENA*wp+(14.0d0*deltp*DENA
     $         +15.0d0*delt-37.0d0)*w)*ots
      pgabpa = a*b*DENA*(sdelt1*DENA*wp+(7.0d0*deltp*DENA+sdelt1)*w)*on
      pgbbpa = a*b*DENA*(two*sdelt5*DENA*wp+(14.0d0*deltp*DENA+
     $         13.0d0*delt+57.0d0)*w)*ots
      dra = GAMA*(pgbbpa+pgabpa+pgaapa)+eclyp2pa+eclyp1pa
      dgaa = pgaa
      dgab = pgab
c...
c...  Second derivatives. The LYP functional is linear on the
c...  density gradient, thus the second derivative with respect
c...  to two times the density gradients are zero.
c...
      ws = DENM173*DENM1*ze*zm1*(two*DEN*(zm1*zm1+11.0d0*zm1+65.0d0)
     $     -C*DEN13*(two*DEN13*(zm1+12.0d0)-C))*on
      delts = tn*d**3*zm1**3/(DEN*DEN2)-tt*d*d*DENM83
     $        *zm1*zm1+fn*delt/DEN2
      eclyp1saa = -a*zm1*(zm1*zm1-three*zm1-7.0d0)/DENA*on
      eclyp1sab = -a*zm1*(zm1*zm1-three*zm1+11.0d0)/DENA*on
      eclyp2saa = -stn*a*b*CF*DENA83*(9.0d0*DENA2*ws+42.0d0*DENA*w
     $             p+44.0d0*w)
      eclyp2sab = -stt*a*b*CF*DENA83*(3.0d0*DENA2*ws+14.0d0*DENA*w
     $             p+11.0d0*w)
      pgaasaa = a*b*(two*sdelt5*DENA2*ws+DEN*(14.0d0*deltp*DENA+
     $   15.0d0*delt-37.0d0)*wp+(14.0d0*delts*DENA2+30.0d0*deltp*DENA+
     $   delt-11.0d0)*w)*ots
      pgaasab = a*b*(sdelt5*DENA2*ws+DEN*(7.0d0*deltp*DENA+sdelt5)
     1   *wp+(7.0d0*delts*DENA2+14.0d0*deltp*DENA+sdelt-13.0d0)*w)*oet
      pgabsaa = a*b*(sdelt1*DENA2*ws+DEN*(7.0d0*deltp*DENA+sdelt1)*wp+
     $          (7.0d0*delts*DENA2+14.0d0*deltp*DENA+24.0d0)*w)*on
      pgabsab = a*b*(sdelt1*DENA2*ws+DEN*(7.0d0*deltp*DENA+sdelt1)*wp+
     1      (7.0d0*delts*DENA2+14.0d0*deltp*DENA+sdelt-23.0d0)*w)*on
      pgbbsaa = a*b*(two*sdelt5*DENA2*ws+DEN*(14.0d0*deltp*DENA+
     $   13.0d0*delt+57.0d0)*wp+(14.0d0*delts*DENA2+26.0d0*deltp*DENA
     $   -delt+83.0d0)*w)*ots
      dra2 = GAMA*(pgbbsaa+pgabsaa+pgaasaa)+eclyp2saa+eclyp1saa
      drab = GAMA*(pgabsab+pgaasab+pgaasab)+eclyp2sab+eclyp1sab
      dragaa = pgaapa
      dragbb = pgbbpa
      dragab = pgabpa
c...
      return
      end
c=======================================================================
      subroutine lypdu(ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab,drara
     $  ,drbrb,drarb,drasaa,drasbb,drasab,drbsaa,drbsbb,drbsab)
      implicit none
c...  MM (09/19/2003) generated with maxima, file lyp.mc
c...
c...                  Lee, Yang, and Parr (LYP)
c...                   Correlation Fuctional
c...                     open shell version
c...
c...  see  B.G. Johnson, P.M.W. Gill and J.A. Pople,
c...     J. Chem. Phys. 98, 5612 (1993) eqs (A22)--(A33).
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  rb    beta density
c...  saa   alpha density gradient invariant
c...  sbb   beta density gradient invariant
c...  sab   alpha beta beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...  dsab    funct. deriv. w.r.t. alpha beta gradient invariant
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drasab  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  drbsab  funct. mixed 2nd deriv. w.r.t. beta den. and alpha beta grad.
c...
c...  Note that the LYP fuctional is linear in the density gradient,
c...  thus the second derivatives involving the density gradient
c...  two times are zero.
c...
      real*8 ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab
      real*8 drara,drbrb,drarb,drasaa,drasbb,drasab
      real*8 drbsaa,drbsbb,drbsab
      real*8 r,z,zp,zs,ze,w,wp,ws,delt,deltp,delts
      real*8 rm1,r2,rm2,r3,rm3,r13,rm13,r23,rm43,rm53,rm73,r83,rm103,rm1
     $  13
      real*8 z2,zm1,zm2,zm3,zp2
      real*8 ra2,ra53,ra83,ra113
      real*8 rb2,rb53,rb83,rb113
      real*8 eclyp1,eclyp1pa,eclyp1pb,eclyp1saa,eclyp1sbb,eclyp1sab
      real*8 eclyp2,eclyp2pa,eclyp2pb,eclyp2saa,eclyp2sbb,eclyp2sab
      real*8 pgaa,pgaapa,pgaapb,pgaasaa,pgaasbb,pgaasab
      real*8 pgbb,pgbbpa,pgbbpb,pgbbsaa,pgbbsbb,pgbbsab
      real*8 pgab,pgabpa,pgabpb,pgabsaa,pgabsbb,pgabsab
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,thirdm,tthird,sixth,et,th,ssm,eht,ft
      parameter (third=1.0d0/3.0d0,thirdm=-1.0d0/3.0d0,tthird=2.0d0/3.0d
     $  0,sixth=1.0d0/6.0d0)
      parameter (et=11.0d0/3.0d0,th=1.5d0,ssm=-7.0d0/6.0d0,eht=-8.0d0/3.
     $  0d0,ft=5.0d0/3.0d0)
c...
      r = rb+ra
c...
c...  Some povers of r, ra and rb
c...
      rm1 = one/r
      r2 = r**2
      r3 = r*r2
      rm2 = rm1**2
      rm3 = rm1*rm2
      r13 = r**THIRD
      rm13 = one/r13
      r23 = r13**2
      rm43 = rm1*rm13
      rm53 = rm13*rm43
      rm73 = rm1*rm43
      r83 = r2*r23
      rm103 = rm1*rm73
      rm113 = rm103*rm13
      ra2 = ra**2
      rb2 = rb**2
      ra53 = ra**ft
      rb53 = rb**ft
      ra83 = ra*ra53
      rb83 = rb*rb53
      ra113 = ra*ra83
      rb113 = rb*rb83
c...
c...  Auxiliary functions z,ze,w, and delt
c...
      z = 3.49D-1*rm13+1.0D0
      z2 = z**2
      zm1 = one/z
      zm2 = zm1**2
      zm3 = zm1*zm2
      ze = EXP(-2.53299999999999997D-1*rm13)
      w = rm113*ze*zm1
      delt = rm13*(2.533D-1*z+3.49D-1)*zm1
c...
c...  Functional value
c...
      eclyp1 = -1.9672D-1*ra*rb*rm1*zm1
      eclyp2 = -2.36705143194386D-1*ra*rb*(rb83+ra83)*w
      pgaa = 7.213066666666667D-4*rb*(r*(9.0D0*rb+(3.0D0*delt-1.0D0)*ra)
     $  +(delt-1.1D1)*ra2)*rm1*w
      pgbb = 7.213066666666667D-4*ra*((delt-1.1D1)*rb2+r*((3.0D0*delt-1.
     $  0D0)*rb+9.0D0*ra))*rm1*w
      pgab = 7.213066666666667D-4*((7.0D0*delt-4.7D1)*ra*rb+1.2D1*r2)*w
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
      ec = rm1*(pgbb*sbb+pgab*sab+pgaa*saa+eclyp2+eclyp1)
c...
c...
c...  First derivatives
c...
      zp = -1.163333333333333D-1*rm43
      wp = 3.333333333333333D-1*rm2*w*(r23*(2.533D-1*z+3.49D-1)-1.1D1*r*
     $  z)*zm1
      deltp = -3.333333333333333D-1*rm53*(delt*r23*z2-1.21801D-1)*zm2
c...
      eclyp1pa = 1.9672D-1*rb*rm2*zm2*(ra*(r*zp+z)-r*z)
      eclyp1pb = 1.9672D-1*ra*rm2*zm2*(rb*(r*zp+z)-r*z)
      eclyp2pa = -7.890171439812868D-2*rb*(ra83*(3.0D0*ra*wp+1.1D1*w)+3.
     $  0D0*rb83*(ra*wp+w))
      eclyp2pb = -7.890171439812868D-2*ra*(rb83*(3.0D0*rb*wp+1.1D1*w)+3.
     $  0D0*ra83*(rb*wp+w))
      pgaapa = 7.213066666666667D-4*rb*rm2*((r2*(9.0D0*rb+3.0D0*delt*ra-
     $  ra)+(delt-1.1D1)*r*ra2)*wp+(-1.0D0*(delt-1.1D1)*ra2+r2*(3.
     $  0D0*deltp*ra+3.0D0*delt-1.0D0)+r*ra*(deltp*ra+2.0D0*delt-2.2D1))
     $  *w)
      pgaapb = 7.213066666666667D-4*rm2*((r2*rb*(9.0D0*rb+3.0D0*delt*ra-
     $  ra)+(delt-1.1D1)*r*ra2*rb)*wp+(r2*(3.0D0*deltp*ra*rb+1.8D1
     $  *rb+3.0D0*delt*ra-ra)+r*ra2*(deltp*rb+delt-1.1D1)-(d
     $  elt-1.1D1)*ra2*rb)*w)
      pgbbpa = 7.213066666666667D-4*rm2*(((delt-1.1D1)*r*ra*rb2+r2*ra*(3
     $  .0D0*delt*rb-rb+9.0D0*ra))*wp+(r*(deltp*ra+delt-1.1D1)*rb2
     $  -(delt-1.1D1)*ra*rb2+r2*(3.0D0*deltp*ra*rb+3.0D0*delt*rb-
     $  rb+1.8D1*ra))*w)
      pgbbpb = 7.213066666666667D-4*ra*rm2*(((delt-1.1D1)*r*rb2+r2*(3.0D
     $  0*delt*rb-rb+9.0D0*ra))*wp+(-1.0D0*(delt-1.1D1)*rb2+r2*(3.
     $  0D0*deltp*rb+3.0D0*delt-1.0D0)+r*rb*(deltp*rb+2.0D0*delt-2.2D1))
     $  *w)
      pgabpa = 7.213066666666667D-4*(((7.0D0*delt-4.7D1)*ra*rb+1.2D1*r2)
     $  *wp+((7.0D0*deltp*ra+7.0D0*delt-4.7D1)*rb+2.4D1*r)*w)
      pgabpb = 7.213066666666667D-4*(((7.0D0*delt-4.7D1)*ra*rb+1.2D1*r2)
     $  *wp+(ra*(7.0D0*deltp*rb+7.0D0*delt-4.7D1)+2.4D1*r)*w)
c...
c...  First derivatives of the functional
c...
      dra = pgbbpa*sbb+pgabpa*sab+pgaapa*saa+eclyp2pa+eclyp1pa
      drb = pgbbpb*sbb+pgabpb*sab+pgaapb*saa+eclyp2pb+eclyp1pb
      dsaa = pgaa
      dsbb = pgbb
      dsab = pgab
c...
c...  Second derivatives
c...
      zp2 = zp**2
      zs = 1.551111111111111D-1*rm73
      ws = 1.111111111111111D-1*w*(1.54D2*r23*z2+6.416089D-2*z2-2.6D1*r1
     $  3*z*(2.533D-1*z+3.49D-1)+1.768034D-1*z+2.43602D-1)*zm2/r83
      delts = 1.111111111111111D-1*rm103*(4.0D0*r*(2.533D-1*z+3.49D-1)*z
     $  2-7.30806D-1*r23*z+8.5017098D-2*r13)*zm3
c...
      eclyp1saa = 1.9672D-1*rb*rm3*zm3*(ra*(r2*z*zs-2.0D0*r2*zp2-2.0D0*r
     $  *z*zp-2.0D0*z2)+2.0D0*r*z*(r*zp+z))
      eclyp1sbb = 1.9672D-1*ra*rm3*zm3*(rb*(r2*z*zs-2.0D0*r2*zp2-2.0D0*r
     $  *z*zp-2.0D0*z2)+2.0D0*r*z*(r*zp+z))
      eclyp1sab = 1.9672D-1*rm3*zm3*(ra*(rb*(r2*(z*zs-2.0D0*zp2)-2.0D0*r
     $  *z*zp-2.0D0*z2)+r2*z*zp+r*z2)+rb*(r2*z*zp+r*z2)-r2*z2)
      eclyp2saa = -2.630057146604289D-2*rb*(ra53*(9.0D0*ra2*ws+6.6D1*ra*
     $  wp+8.8D1*w)+9.0D0*rb83*(ra*ws+2.0D0*wp))
      eclyp2sbb = -2.630057146604289D-2*ra*(rb53*(9.0D0*rb2*ws+6.6D1*rb*
     $  wp+8.8D1*w)+9.0D0*ra83*(rb*ws+2.0D0*wp))
      eclyp2sab = -7.890171439812868D-2*(3.0D0*ra*rb*(rb83+ra83)*ws+(1.1
     $  D1*ra*rb83+3.0D0*rb113+1.1D1*ra83*rb+3.0D0*ra113)*wp+1.1D1*(rb83
     $  +ra83)*w)
      pgaasaa = 7.213066666666667D-4*rb*rm3*((r3*(9.0D0*rb+3.0D0*delt*ra
     $  -ra)+(delt-1.1D1)*r2*ra2)*ws+(-2.0D0*(delt-1.1D1)*r*ra2+2.
     $  0D0*r3*(3.0D0*deltp*ra+3.0D0*delt-1.0D0)+2.0D0*r2*ra*(deltp*ra+2
     $  .0D0*delt-2.2D1))*wp+(r2*(delts*ra2+4.0D0*deltp*ra+2.0D0*delt-2.
     $  2D1)+2.0D0*(delt-1.1D1)*ra2+3.0D0*r3*(delts*ra+2.0D0*deltp)-2.0D
     $  0*r*ra*(deltp*ra+2.0D0*delt-2.2D1))*w)
      pgaasbb = 7.213066666666667D-4*rm3*((r3*rb*(9.0D0*rb+3.0D0*delt*ra
     $  -ra)+(delt-1.1D1)*r2*ra2*rb)*ws+(2.0D0*r3*(3.0D0*deltp*ra*
     $  rb+1.8D1*rb+3.0D0*delt*ra-ra)+2.0D0*r2*ra2*(deltp*rb+delt-
     $  1.1D1)-2.0D0*(delt-1.1D1)*r*ra2*rb)*wp+(3.0D0*r3*(delts*ra*rb+2.
     $  0D0*deltp*ra+6.0D0)+r2*ra2*(delts*rb+2.0D0*deltp)-2.0D0*r*ra2*(d
     $  eltp*rb+delt-1.1D1)+2.0D0*(delt-1.1D1)*ra2*rb)*w)
      pgaasab = 7.213066666666667D-4*rm3*((r3*(rb*(9.0D0*rb-ra)+3.
     $  0D0*delt*ra*rb)+(delt-1.1D1)*r2*ra2*rb)*ws+(r2*(delt*ra*(2.0D0*r
     $  b+ra)-1.1D1*ra*(2.0D0*rb+ra)+2.0D0*deltp*ra2*rb)+r3*(3.0D0*delt*
     $  (rb+ra)+6.0D0*deltp*ra*rb+1.7D1*rb-ra)+(2.2D1-2.0D0*delt)*
     $  r*ra2*rb)*wp+(r2*(ra*(delts*ra*rb-2.2D1)+deltp*ra*(2.0D0*rb+ra)+
     $  2.0D0*delt*ra)+r*(-delt*ra*(2.0D0*rb+ra)+1.1D1*ra*(2.0D0*r
     $  b+ra)-2.0D0*deltp*ra2*rb)+r3*(3.0D0*deltp*(rb+ra)+3.0D0*delts*ra
     $  *rb+3.0D0*delt-1.0D0)+(2.0D0*delt-2.2D1)*ra2*rb)*w)
      pgbbsaa = 7.213066666666667D-4*rm3*(((delt-1.1D1)*r2*ra*rb2+r3*ra*
     $  (3.0D0*delt*rb-rb+9.0D0*ra))*ws+(2.0D0*r2*(deltp*ra+delt-1
     $  .1D1)*rb2-2.0D0*(delt-1.1D1)*r*ra*rb2+2.0D0*r3*(3.0D0*deltp*ra*r
     $  b+3.0D0*delt*rb-1.0D0*rb+1.8D1*ra))*wp+(r2*(delts*ra+2.0D0*deltp
     $  )*rb2-2.0D0*r*(deltp*ra+delt-1.1D1)*rb2+2.0D0*(delt-1.1D1)*ra*rb
     $  2+3.0D0*r3*(delts*ra*rb+2.0D0*deltp*rb+6.0D0))*w)
      pgbbsbb = 7.213066666666667D-4*ra*rm3*(((delt-1.1D1)*r2*rb2+r3*(3.
     $  0D0*delt*rb-rb+9.0D0*ra))*ws+(-2.0D0*(delt-1.1D1)*r*rb2+2.
     $  0D0*r3*(3.0D0*deltp*rb+3.0D0*delt-1.0D0)+2.0D0*r2*rb*(deltp*rb+2
     $  .0D0*delt-2.2D1))*wp+(r2*(delts*rb2+4.0D0*deltp*rb+2.0D0*delt-2.
     $  2D1)+2.0D0*(delt-1.1D1)*rb2+3.0D0*r3*(delts*rb+2.0D0*deltp)-2.0D
     $  0*r*rb*(deltp*rb+2.0D0*delt-2.2D1))*w)
      pgbbsab = 7.213066666666667D-4*rm3*(((delt-1.1D1)*r2*ra*rb2+r3*(3.
     $  0D0*delt*ra*rb-ra*(rb-9.0D0*ra)))*ws+(r2*(2.0D0*deltp*ra*r
     $  b2+delt*rb*(rb+2.0D0*ra)-1.1D1*rb*(rb+2.0D0*ra))+(2.2D1-2.0D0*de
     $  lt)*r*ra*rb2+r3*(3.0D0*delt*(rb+ra)+6.0D0*deltp*ra*rb-rb+1
     $  .7D1*ra))*wp+(r*(-2.0D0*deltp*ra*rb2-delt*rb*(rb+2.0D0*ra)
     $  +1.1D1*rb*(rb+2.0D0*ra))+(2.0D0*delt-2.2D1)*ra*rb2+r2*(rb*(delts
     $  *ra*rb-2.2D1)+deltp*rb*(rb+2.0D0*ra)+2.0D0*delt*rb)+r3*(3.0D0*de
     $  ltp*(rb+ra)+3.0D0*delts*ra*rb+3.0D0*delt-1.0D0))*w)
      pgabsaa = 7.213066666666667D-4*(((7.0D0*delt-4.7D1)*ra*rb+1.2D1*r2
     $  )*ws+(2.0D0*(7.0D0*deltp*ra+7.0D0*delt-4.7D1)*rb+4.8D1*r)*wp+(7.
     $  0D0*(delts*ra+2.0D0*deltp)*rb+2.4D1)*w)
      pgabsbb = 7.213066666666667D-4*(((7.0D0*delt-4.7D1)*ra*rb+1.2D1*r2
     $  )*ws+(2.0D0*ra*(7.0D0*deltp*rb+7.0D0*delt-4.7D1)+4.8D1*r)*wp+(7.
     $  0D0*ra*(delts*rb+2.0D0*deltp)+2.4D1)*w)
      pgabsab = 7.213066666666667D-4*((7.0D0*delt*ra*rb-4.7D1*ra*rb+1.2D
     $  1*r2)*ws+(1.4D1*deltp*ra*rb-4.7D1*rb-4.7D1*ra+7.0D0*delt*r+4.8D1
     $  *r)*wp+(7.0D0*delts*ra*rb+7.0D0*deltp*r+7.0D0*delt-2.3D1)*w)
c...
c...  Second derivatives of the functional
c...
      drara = pgbbsaa*sbb+pgabsaa*sab+pgaasaa*saa+eclyp2saa+eclyp1saa
      drbrb = pgbbsbb*sbb+pgabsbb*sab+pgaasbb*saa+eclyp2sbb+eclyp1sbb
      drarb = pgbbsab*sbb+pgabsab*sab+pgaasab*saa+eclyp2sab+eclyp1sab
      drasaa = pgaapa
      drasbb = pgbbpa
      drasab = pgabpa
      drbsaa = pgaapb
      drbsbb = pgbbpb
      drbsab = pgabpb
c...
      return
      end
c=======================================================================

      SUBROUTINE OPTXD(DA,    DA13,   DA43,   g1A,    gama,
     $                 DB,    DB13,   DB43,   g1B,    gamb,
     $  rhf,   EXOPTX, dfdda, dfdga,  dfdda2, dfdga2, dfddaga,
     $                 dfddb, dfdgb,  dfddb2, dfdgb2, dfddbgb)
      IMPLICIT REAL*8(A-H,O-Z)
C
C   MM (04/23/2003)
C
C  Handy and Cohen optimized nonlocal exchange
C
C  this subroutine computes the functional and its first and
C  second derivatives for both closed and open shell
C
C  ARGUMENTS
C
C  DA      -  alpha density at current grid point
C  DA13    -  DA**(1/3)
C  DA43    -  DA**(4/3)
C  g1A     -  absolute value of alpha density derivative
C  gama    -  square of g1a
C  DB      -  beta density at current grid point
C  DB13    -  DB**(1/3)
C  DB43    -  DB**(4/3)
C  g1B     -  absolute value of beta density derivative
C  gamb    -  square of g1b
C  rhf     -  logical flag for rhf/uhf
C
C  on exit
C
C  EXOPTX  -  contribution to the exchange energy, that is,
C             the functional values divided by the total density
C             at the current grid point.
C  dfdda   -  functional derivative with respect the alpha density
C  dfdga   -  derivative with respect the alpha density gradient
C  dfdda2  -  second derivative with respect the alpha density
C  dfdga2  -  second derivative with respect the alpha gradient
C  dfddaga -  second derivative with respect density and gradient
C
C  dfddb   -  functional derivative with respect the beta density
C  dfdgb   -  derivative with respect the beta density gradient
C  dfddb2  -  second derivative with respect the beta density
C  dfdgb2  -  second derivative with respect the beta gradient
C  dfddbgb -  second derivative with respect density and gradient
C
      logical rhf
      parameter(Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0,Four=4.0d0)
      parameter(TThird=2.0d0/3.0d0,FThird=4.0d0/3.0d0,thrsh=1.0d-12)
      parameter(Fourth=1.0d0/4.0d0,FNine=4.0d0/9.0d0)
c
      PARAMETER (a=-1.43169d0,b=0.006d0)
c
c...
      if(rhf)then
c...
c...  Closed Shell
c...
        xa = g1a/da43
        xa2 = xa*xa
        dgfa = b*xa2+One
        gfa = a*b**2*xa2*xa2/dgfa**2
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
        exoptx = da13*gfa
c...
c...  First derivatives of the functional (potential)
c...
        gfpa = -four*a*b**2*xa2*xa*(b*xa2-dgfa)/dgfa**3
        dfdda = -FThird*(g1a*gfpa-da43*gfa)/da
        if(g1a.gt.thrsh)then
          dfdga = Half*gfpa/g1a
        else
          dfdga = ZERO
        endif
c...
c...  Second derivatives of the functional
c...
        gfsa = 12.0d0*a*b**2*xa2*(b*xa2-dgfa)*(Two*b*xa2-dgfa)/dgfa**4
        dfdda2 = FNine*(Four*da13*gama*gfsa-
     $     da43*da13*g1a*gfpa+da**3*gfa)/(da43*da43*da)
        if(g1a.gt.thrsh)then
          dfdga2 = Fourth*(g1a*gfsa-da43*gfpa)/(da43*g1a**3)
        else
          dfdga2 = ZERO
        endif
        dfddaga = -TThird*gfsa/(da43*da)
c...
      else
c...
c...  Open Shell - Alpha Density
c...
        if(da.gt.thrsh)then
          xa = g1a/da43
          xa2 = xa*xa
          dgfa = b*xa2+One
          gfa = a*b**2*xa2*xa2/dgfa**2
c...
c...  First derivatives of the functional (potential)
c...
          gfpa = -four*a*b**2*xa2*xa*(b*xa2-dgfa)/dgfa**3
          dfdda = -FThird*(g1a*gfpa-da43*gfa)/da
          if(g1a.gt.thrsh)then
            dfdga = Half*gfpa/g1a
          else
            dfdga = ZERO
          endif
c...
c...  Second derivatives of the functional
c...
          gfsa = 12.0d0*a*b**2*xa2*(b*xa2-dgfa)*(Two*b*xa2-dgfa)/dgfa**4
          dfdda2 = FNine*(Four*da13*gama*gfsa-
     1       da43*da13*g1a*gfpa+da**3*gfa)/(da43*da43*da)
          if(g1a.gt.thrsh)then
            dfdga2 = Fourth*(g1a*gfsa-da43*gfpa)/(da43*g1a**3)
          else
            dfdga2 = ZERO
          endif
          dfddaga = -TThird*gfsa/(da43*da)
        else
          gfa = ZERO
          dfdda = ZERO
          dfdda2 = ZERO
          dfdga = ZERO
          dfdga2 = ZERO
          dfddaga = ZERO
        endif
c...
c...  Open Shell - Beta Density
c...
        if(db.gt.thrsh)then
          xb = g1b/db43
          xb2 = xb*xb
          dgfb = b*xb2+One
          gfb = a*b**2*xb2*xb2/dgfb**2
c...
c...  First derivatives of the functional (potential)
c...
          gfpb = -four*a*b**2*xb2*xb*(b*xb2-dgfb)/dgfb**3
          dfddb = -Fthird*(g1b*gfpb-db43*gfb)/db
          if(g1b.gt.thrsh)then
            dfdgb = Half*gfpb/g1b
          else
            dfdgb = ZERO
          endif
c...
c...  Second derivatives of the functional
c...
          gfsb = 12.0d0*a*b**2*xb2*(b*xb2-dgfb)*(Two*b*xb2-dgfb)/dgfb**4
          dfddb2 = Fnine*(Four*db13*gamb*gfsb-
     1       db43*db13*g1b*gfpb+db**3*gfb)/(db43*db43*db)
          if(g1b.gt.thrsh)then
            dfdgb2 = Fourth*(g1b*gfsb-db43*gfpb)/(db43*g1b**3)
          else
            dfdgb2 = ZERO
          endif
          dfddbgb = -TThird*gfsb/(db43*db)
        else
          gfb = ZERO
          dfddb = ZERO
          dfddb2 = ZERO
          dfdgb = ZERO
          dfdgb2 = ZERO
          dfddbgb = ZERO
        endif
c...
c...  Functional value (divided by the density, to be consistent
c...  with the PQS implementation)
c...
        exoptx = (db43*gfb+da43*gfa)/(db+da)
      endif
c...
      return
      end
c=======================================================================
      subroutine pw91nlcdu(ra,rb,saa,sbb,sab,scl,scnl,ec,
     $                     dra,drb,dsaa,dsbb,dsab,
     $                     drara,drbrb,drarb,drasaa,drasbb,drasab,
     $                     drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $                     dsabsab,dsaasbb,dsaasab,dsbbsab)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (04/29/2003) generated with maxima, file pw91lcu.mc
c...
c...  Perdew and Wang nonlocal correlation functional
c...  *** open shell version
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...
c...  Formulation taken from the Molpro manual (PW91C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
c...
c...  Input parameters:
c...
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  saa     beta density gradient invariant
c...  sab     alpha beta density gradient invariant
c...  scl     scaling factor for local part
c...  scnl    scaling factor for non local part
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drasab  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  drbsab  funct. mixed 2nd deriv. w.r.t. beta den. and alpha beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...  dsaasab funct. mixed 2nd deriv. w.r.t. alpha and alpha beta grad.
c...  dsbbsab funct. mixed 2nd deriv. w.r.t. beta and alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,sab,scl,scnl,ec,dra,drb,dsaa,dsbb,dsab,
     $       drara,drbrb,drarb,drasaa,drasbb,drasab,
     $       drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb,
     $       dsabsab,dsaasbb,dsaasab,dsbbsab
      real*8 pw91lc,dpw91lcdra,dpw91lcdrb,dpw91lcdrara,
     $       dpw91lcdrbrb,dpw91lcdrarb
      real*8 r,r2,r3,r4,r5,r6,r7,r12,r32,r52,r72,r16,rm13,r13,r23,r56,
     $       rm76,r76,rm43,r43,r53,r73,r83,r136,r113,r133,rm193,
     $       rm196,rm296
      real*8 s,s12,sm32
      real*8 x,x2,x3
      real*8 z
      real*8 zp1,zp113,zp123,zp143
      real*8 zm1,zm113,zm123,zm143
      real*8 om,omp,oms,om2,om3,om4,om6,om7,om8
      real*8 d,d2,d4
      real*8 e,epra,eprb,esra,esrb,esrarb
      real*8 numt,dent,dentp,t,tp,ts
      real*8 polj,expj
      real*8 j,jpd,jpom,jpx,jpr,jsd,jsom,jsx,jsr,jsdom,jsdx,jsdr,
     $       jsxom,jsxr,jsomr
      real*8 djdra,djdrb,djdsaa,djdrara,djdrbrb,djdrarb,djdrasaa,
     $       djdrbsaa,djdsaasaa
      real*8 expa,a,ape,apom,ase,asom,aseom
      real*8 numlog,denlog
      real*8 loga,logape,logapd,logapom,logase,logasd,logasom,
     $       logased,logaseom,logasomd
      real*8 l,lpe,lpom,lpd,lse,lsom,lsd,lseom,lsed,lsomd
      real*8 dldra,dldrb,dldsaa,dldrara,dldrbrb,dldrarb,dldrasaa,
     $       dldrbsaa,dldsaasaa
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,sixth,etm
      parameter (third=1.0d0/3.0d0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0)
c...
c... threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      call pw91lcdu(ra,rb,pw91lc,dpw91lcdra,dpw91lcdrb,dpw91lcdrara,
     $              dpw91lcdrbrb,dpw91lcdrarb)
c...
c...   Now, if the total gradient is bigger than the threshold,
c...   we compute the non local part, otherwise we are basically done
c...
      r = rb+ra
      s=saa+sbb+sab+sab
      if(s.lt.epsi.or.scnl.eq.zero)then
        j = ZERO
        djdra = ZERO
        djdrb = ZERO
        djdsaa = ZERO
        djdrara = ZERO
        djdrbrb = ZERO
        djdrarb = ZERO
        djdrasaa = ZERO
        djdrbsaa = ZERO
        djdsaasaa = ZERO
        l = ZERO
        dldra = ZERO
        dldrb = ZERO
        dldsaa = ZERO
        dldrara = ZERO
        dldrbrb = ZERO
        dldrarb = ZERO
        dldrasaa = ZERO
        dldrbsaa = ZERO
        dldsaasaa = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
        r2 = r*r
        r3 = r*r2
        r4 = r*r3
        r5 = r*r4
        r6 = r*r5
        r7 = r*r6
        r12 = SQRT(r)
        r32 = r*r12
        r52 = r12*r2
        r72 = r12*r3
        r16 = r**SIXTH
        r13 = r16*r16
        rm13 = one/r13
        r23 = r13*r13
        r56 = r16*r23
        r76 = r*r16
        rm76 = one/r76
        r43 = r*r13
        rm43 = one/r43
        r53 = r13*r43
        r73 = r13*r2
        r83 = r13*r73
        r136 = r16*r2
        r113 = r23*r3
        r133 = r13*r4
        rm193 = one/(r13*r6)
        rm196 = one/(r16*r3)
        rm296 = one/(r4*r56)
c...
        s12 = SQRT(s)
        sm32 = one/(s*s12)
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
        x = 6.2035049089940001667D-1*rm13
        x2 = x*x
        x3 = x2*x
        z = min(max((ra-rb)/r,-one),one)
        zp1 = z+1.0D0
        zp113 = zp1**THIRD
        zp123 = zp113*zp113
        zp143 = zp123*zp123
        zm1 = 1.0D0-z
        zm113 = zm1**THIRD
        zm123 = zm113*zm113
        zm143 = zm123*zm123
c...
c...  the variables omega and d.
c...
        om = 5.0D-1*(zp123+zm123)
        om2 = om*om
        om3 = om*om2
        om4 = om*om3
        om6 = om2*om4
        om7 = om*om6
        om8 = om*om7
        d = 2.5192897034224488769D-1*rm76*s12/om
        d2 = d*d
        d4 = d2*d2
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
        expj = EXP(-4.115631209904114531594D1*d2*om4*rm13)
        polj = t-1.8532685714285714286D-3
        j = 1.5755920349483144659D1*d2*expj*om3*polj
c...
c...  now the term l
c...
c...  function a.
c...
        expa = EXP(-4.042761511756372495028D1*e/om3)
        a = 2.6975860915198740865D0/(expa-1.0D0)
c...
c...  we build l.
c...
        numlog = a*d4+d2
        denlog = a**2*d4+a*d2+1.0D0
        loga = LOG(2.6975860915198740865D0*numlog/denlog+1.0D0)
        l = 2.473556743557577051D-2*loga*om3
c...
c...  First derivative of om
c...
        if(zm113*zp113.gt.epsi)then
          omp = -3.3333333333333333D-1*(zp113-zm113)/(zm113*zp113)
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
        tp = 1.0D-3*(dent*(1.4778D-2*x+2.3266D1)-dentp*numt)/dent*
     1     *2
        jpd = 2.0578156049520572658D-1*j*(9.7190438015295143676D0*r13-4.
     1     0D2*d2*om4)*rm13/d
        jpom = 1.0289078024760286329D-1*j*(2.9157131404588543103D1*r13-1
     1     .6D3*d2*om4)*rm13/om
        jpx = 1.5755920349483144659D1*d2*expj*om3*tp
        jpr = 1.3718770699680381772D1*d2*j*om4*rm43
c...
c...  First derivatives of l w.r.t. auxiliary variables
c...
        ape = 1.0905697225445850371D2*expa/((expa-1.0D0)**2*om3)
        apom = -3.2717091676337551113D2*e*expa/((expa-1.0D0)**2*om4)
        logape = -1.8D-1*ape*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*denl
     1     og)/(denlog*(1.8D-1*numlog+6.672632268006111763D-2*denlog))
        logapom = -1.8D-1*apom*d2*((2.0D0*a*d2+1.0D0)*numlog-d2*de
     1     nlog)/(denlog*(1.8D-1*numlog+6.672632268006111763D-2*denlog))
        logapd = -3.6D-1*d*(2.0D0*a*d2+1.0D0)*(a*numlog-denlog)/(d
     1     enlog*(1.8D-1*numlog+6.672632268006111763D-2*denlog))
        lpe = 2.473556743557577051D-2*logape*om3
        lpom = 2.473556743557577051D-2*(logapom*om+3.0D0*loga)*om2
        lpd = 2.473556743557577051D-2*logapd*om3
c...
c... Total first derivatives of the non local terms:
c...  1) j
c...  2) l
c...
        djdra = -1.4923313454695211986D-2*(2.8135951074925023309D0*jpd*r
     1     13*(1.2D1*omp*rb+7.0D0*om*r)*s12-6.700924717797021006D1*om2*r
     2     32*(2.0D0*jpom*omp*rb+jpr*r2)+1.3856406460551018348D1*jpx*om2
     3     *r136)/(om2*r72)
        djdrb = 1.4923313454695211986D-2*(2.8135951074925023309D0*jpd*r1
     1     3*(1.2D1*omp*ra-7.0D0*om*r)*s12-6.700924717797021006D1*om2*r3
     2     2*(2.0D0*jpom*omp*ra-jpr*r2)-1.3856406460551018348D1*jp
     3     x*om2*r136)/(om2*r72)
        djdsaa = 1.2596448517112244384D-1*jpd*rm76/(om*s12)
        dldra = -3.469513240231684665D-2*rm196*(1.210203242253764276D0*l
     1     pd*(1.2D1*omp*rb+7.0D0*om*r)*s12-2.882248692422406544D1*om2*r
     2     76*(2.0D0*lpom*omp*rb+epra*lpe*r2))/om2
        dldrb = 3.469513240231684665D-2*rm196*(1.210203242253764276D0*lp
     1     d*(1.2D1*omp*ra-7.0D0*om*r)*s12-2.882248692422406544D1*om2*r7
     2     6*(2.0D0*lpom*omp*ra-eprb*lpe*r2))/om2
        dldsaa = 1.2596448517112244384D-1*lpd*rm76/(om*s12)
c...
c...  second derivative of om
c...
        if(zm143*zp143.gt.epsi)then
          oms = -1.1111111111111111111D-1*(z*zp113+zp113-z*zm113+zm1
     1       13)/(zm143*zp143)
        else
          oms=zero
        endif
c...
c...  second derivatives of e
c...
        esra = -(2*epra-dpw91lcdrara)/r
        esrb = -(2*eprb-dpw91lcdrbrb)/r
        esrarb = -(eprb+epra-dpw91lcdrarb)/r
c...
c...  Second derivatives of j w.r.t. auxiliary variables
c...
        ts = -2.0D-3*(numt*(dent*(2.2167D-1*x+4.72D-1)-dentp**2)+d
     1     ent*dentp*(1.4778D-2*x+2.3266D1)-7.389D-3*dent**2)/dent**3
        jsd = 1.4456604259342921967D-2*j*(1.3834507496512901777D2*r23-2.
     1     8468865413150468245D4*d2*om4*r13+4.6866940401968744417D5*d4*o
     2     m8)/(d2*r23)
        jsom = 1.4456604259342921967D-2*j*(4.1503522489538705332D2*r23-1
     1     .0248791548734168568D5*d2*om4*r13+1.8746776160787497767D6*d4*
     2     om8)/(om2*r23)
        jsx = 1.5755920349483144659D1*d2*expj*om3*ts
        jsr = -1.2850314897193708415D0*d2*j*om4*(1.4234432706575234123D1
     1     *r13-1.464591887561523263D2*d2*om4)/r83
        jsdom = 1.4456604259342921967D-2*j*(4.1503522489538705332D2*r23-
     1     6.263150390893103014D4*d2*om4*r13+9.3733880803937488833D5*d4*
     2     om8)/(d*om*r23)
        jsdx = 3.2422778765548086862D0*d*expj*om3*(9.7190438015295143676
     1     D0*r13-4.0D2*d2*om4)*rm13*tp
        jsdr = 3.8550944691581125245D0*d*j*om4*(1.4234432706575234123D1*
     1     r13-2.929183775123046526D2*d2*om4)/r53
        jsxom = 1.6211389382774043431D0*d2*expj*om2*(2.91571314045885431
     1     03D1*r13-1.6D3*d2*om4)*rm13*tp
        jsxr = 2.1615185843698724575D2*d4*expj*om7*rm43*tp
        jsomr = 9.6377361728952813112D-1*d2*j*om3*(9.9641028946026638859
     1     D1*r13-2.3433470200984372208D3*d2*om4)/r53
c...
c...  Second derivatives of l w.r.t. auxiliary variables
c...
        ase = 4.4089133001900743116D3*expa*(expa+1.0D0)/((expa-1.0D0)**3
     1     *om6)
        asom = 1.4696377667300247705D5*e*expa*(expa*(8.90480427680727738
     1     37D-3*om3+2.7D-1*e)-8.9048042768072773837D-3*om3+2.7D-1*e)/((
     2     expa-1.0D0)**3*om8)
        aseom = -7.3481888336501238526D4*expa*(expa*(4.45240213840363869
     1     19D-3*om3+1.8D-1*e)-4.4524021384036386919D-3*om3+1.8D-1*e)/((
     2     expa-1.0D0)**3*om7)
        logase = -1.8D-1*d2*(1.8D-1*(2.0D0*a*ase*d2+2.0D0*ape**2*d2+ase)
     1     *denlog*numlog**2-1.8D-1*ape**2*d2*(2.0D0*a*d2+1.0D0)**2*numl
     2     og**2+(1.3345264536012223526D-1*a*ase*d2-1.8D-1*ase*d2+1.3345
     3     264536012223526D-1*ape**2*d2+6.672632268006111763D-2*ase)*den
     4     log**2*numlog-1.3345264536012223526D-1*ape**2*d2*(2.0D0*a*d2+
     5     1.0D0)**2*denlog*numlog-6.672632268006111763D-2*ase*d2*denlog
     6     **3+2.0D0*ape**2*(1.3345264536012223526D-1*a*d2+9.0D-2*d2+6.6
     7     72632268006111763D-2)*d4*denlog**2)/(denlog**2*(1.8D-1*numlog
     8     +6.672632268006111763D-2*denlog)**2)
        logasom = -1.8D-1*d2*(1.8D-1*(2.0D0*a*asom*d2+2.0D0*apom**2*d2+a
     1     som)*denlog*numlog**2-1.8D-1*apom**2*d2*(2.0D0*a*d2+1.0D0)**2
     2     *numlog**2+(1.3345264536012223526D-1*a*asom*d2-1.8D-1*asom*d2
     3     +1.3345264536012223526D-1*apom**2*d2+6.672632268006111763D-2*
     4     asom)*denlog**2*numlog-1.3345264536012223526D-1*apom**2*d2*(2
     5     .0D0*a*d2+1.0D0)**2*denlog*numlog-6.672632268006111763D-2*aso
     6     m*d2*denlog**3+2.0D0*apom**2*(1.3345264536012223526D-1*a*d2+9
     7     .0D-2*d2+6.672632268006111763D-2)*d4*denlog**2)/(denlog**2*(1
     8     .8D-1*numlog+6.672632268006111763D-2*denlog)**2)
        logasd = -3.6D-1*(a*numlog-denlog)*(denlog*(1.8D-1*(6.0D0*
     1     a*d2+1.0D0)*numlog-4.0D0*(6.672632268006111763D-2*a+9.0D-2)*d
     2     2*(2.0D0*a*d2+1.0D0)**2)-3.6D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*nu
     3     mlog+6.672632268006111763D-2*(6.0D0*a*d2+1.0D0)*denlog**2)/(d
     4     enlog**2*(1.8D-1*numlog+6.672632268006111763D-2*denlog)**2)
        logaseom = -1.8D-1*d2*(denlog*((1.8D-1*aseom*(2.0D0*a*d2+1.0D0)+
     1     3.6D-1*ape*apom*d2)*numlog**2-1.3345264536012223526D-1*ape*ap
     2     om*d2*(2.0D0*a*d2+1.0D0)**2*numlog)-1.8D-1*ape*apom*d2*(2.0D0
     3     *a*d2+1.0D0)**2*numlog**2+denlog**2*((aseom*(1.33452645360122
     4     23526D-1*a*d2-1.8D-1*d2+6.672632268006111763D-2)+1.3345264536
     5     012223526D-1*ape*apom*d2)*numlog+2.0D0*ape*apom*(1.3345264536
     6     012223526D-1*a*d2+9.0D-2*d2+6.672632268006111763D-2)*d4)-6.67
     7     2632268006111763D-2*aseom*d2*denlog**3)/(denlog**2*(1.8D-1*nu
     8     mlog+6.672632268006111763D-2*denlog)**2)
        logased = -3.6D-1*ape*d*(denlog*(1.8D-1*(4.0D0*a*d2+1.0D0)*numlo
     1     g**2-1.3345264536012223526D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*numl
     2     og)-1.8D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*numlog**2+denlog**2*((2
     3     .6690529072024447052D-1*a*d2-3.6D-1*d2+6.672632268006111763D-
     4     2)*numlog+d2*(2.0017896804018335289D-1*a*d2+1.8D-1*d2+6.67263
     5     2268006111763D-2)*(2.0D0*a*d2+1.0D0))-1.3345264536012223526D-
     6     1*d2*denlog**3)/(denlog**2*(1.8D-1*numlog+6.67263226800611176
     7     3D-2*denlog)**2)
        logasomd = -3.6D-1*apom*d*(denlog*(1.8D-1*(4.0D0*a*d2+1.0D0)*num
     1     log**2-1.3345264536012223526D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*nu
     2     mlog)-1.8D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*numlog**2+denlog**2*(
     3     (2.6690529072024447052D-1*a*d2-3.6D-1*d2+6.672632268006111763
     4     D-2)*numlog+d2*(2.0017896804018335289D-1*a*d2+1.8D-1*d2+6.672
     5     632268006111763D-2)*(2.0D0*a*d2+1.0D0))-1.3345264536012223526
     6     D-1*d2*denlog**3)/(denlog**2*(1.8D-1*numlog+6.672632268006111
     7     763D-2*denlog)**2)
        lse = 2.473556743557577051D-2*logase*om3
        lsom = 2.473556743557577051D-2*om*(logasom*om2+6.0D0*logapom*om+
     1     6.0D0*loga)
        lsd = 2.473556743557577051D-2*logasd*om3
        lseom = 2.473556743557577051D-2*(logaseom*om+3.0D0*logape)*om2
        lsed = 2.473556743557577051D-2*logased*om3
        lsomd = 2.473556743557577051D-2*(logasomd*om+3.0D0*logapd)*om2
c...
c... Total second derivatives of the non local terms:
c...  1) j
c...  2) l
c...
        djdrara = 2.2270528446708714288D-4*((2.5959215311134004471D0*r56
     1     *(-1.210468599538321271D1*jpd*om*r*(1.44D2*om*oms*rb**2-2.88D
     2     2*omp**2*rb**2-3.12D2*om*omp*r*rb-9.1D1*om2*r2)-1.45256231944
     3     59855252D2*om2*r*(1.2D1*omp*rb+7.0D0*om*r)*(2.0D0*jsdom*omp*r
     4     b+jsdr*r2))+7.7972634849667692456D1*jsdx*om2*r52*(1.2D1*omp*r
     5     b+7.0D0*om*r))*s12+7.916317428905745746D0*jsd*r23*(1.2D1*omp*
     6     rb+7.0D0*om*r)**2*s+2.1450293971110256001D0*(8.37329169177049
     7     95482D3*om4*r3*rb*(jpom*oms*rb+jsom*omp**2*rb+jsomr*omp*r2-
     8     jpom*omp*r)+2.093322922942624887D3*jsr*om4*r7)-6.19004910
     9     34232434203D2*om4*r113*(6.0D0*jsxom*omp*rb+3.0D0*jsxr*r2-2.0D
     :     0*jpx*r)+1.92D2*jsx*om4*r133)/(om4*r7)
        djdrbrb = 2.2270528446708714288D-4*((2.5959215311134004471D0*r56
     1     *(-1.210468599538321271D1*jpd*om*r*(1.44D2*om*oms*ra**2-2.88D
     2     2*omp**2*ra**2+3.12D2*om*omp*r*ra-9.1D1*om2*r2)-1.45256231944
     3     59855252D2*om2*r*(1.2D1*omp*ra-7.0D0*om*r)*(2.0D0*jsdom*omp*r
     4     a-jsdr*r2))-7.7972634849667692456D1*jsdx*om2*r52*(1.2D1
     5     *omp*ra-7.0D0*om*r))*s12+7.916317428905745746D0*jsd*r23*(1.2D
     6     1*omp*ra-7.0D0*om*r)**2*s+2.1450293971110256001D0*(8.37329169
     7     17704995482D3*om4*r3*ra*(jpom*oms*ra+jsom*omp**2*ra-jso
     8     mr*omp*r2+jpom*omp*r)+2.093322922942624887D3*jsr*om4*r7)+6.19
     9     00491034232434203D2*om4*r113*(6.0D0*jsxom*omp*ra-3.0D0*jsxr*r
     :     2+2.0D0*jpx*r)+1.92D2*jsx*om4*r133)/(om4*r7)
        djdrarb = -2.2270528446708714288D-4*((2.5959215311134004471D0*r5
     1     6*(-1.210468599538321271D1*jpd*om*r*(1.44D2*om*oms*ra*rb-2.88
     2     D2*omp**2*ra*rb+1.56D2*om*omp*r*rb-1.56D2*om*omp*r*ra+9.1D1*o
     3     m2*r2)-1.4525623194459855252D2*jsdom*om2*omp*r*(2.4D1*omp*ra*
     4     rb-7.0D0*om*r*rb+7.0D0*om*r*ra)+1.4525623194459855252D2*jsdr*
     5     om2*r3*(6.0D0*omp*rb-6.0D0*omp*ra+7.0D0*om*r))-7.797263484966
     6     7692456D1*jsdx*om2*r52*(6.0D0*omp*rb-6.0D0*omp*ra+7.0D0*om*r)
     7     )*s12+7.916317428905745746D0*jsd*r23*(1.2D1*omp*ra-7.0D0*om*r
     8     )*(1.2D1*omp*rb+7.0D0*om*r)*s+2.1450293971110256001D0*(4.1866
     9     458458852497741D3*jpom*om4*r3*(2.0D0*oms*ra*rb+omp*r*rb-
     :     omp*r*ra)+2.093322922942624887D3*om4*r3*(4.0D0*jsom*omp**2*r
     ;     a*rb-jsr*r4)-4.1866458458852497741D3*jsomr*om4*omp*r5*(
     <     rb-ra))+1.464591887561523263D0*r23*(1.26794006357553664
     =     66D3*jsxom*om4*omp*r3*(rb-ra)+4.2264668785851221554D2*o
     >     m4*(3.0D0*jsxr*r-2.0D0*jpx)*r4)-1.92D2*jsx*om4*r133)/(om4*r7)
        djdrasaa = -1.8798074963679670784D-3*rm296*(2.813595107492502330
     1     9D0*jsd*r12*(1.2D1*omp*rb+7.0D0*om*r)*s12+1.46459188756152326
     2     3D0*r23*(7.6254743439754925901D0*jpd*om*r*(1.2D1*omp*rb+7.0D0
     3     *om*r)-4.575284606385295554D1*om2*r*(2.0D0*jsdom*omp*rb+jsdr*
     4     r2))+1.3856406460551018348D1*jsdx*om2*r73)/(om3*s12)
        djdrbsaa = 1.8798074963679670784D-3*rm296*(2.8135951074925023309
     1     D0*jsd*r12*(1.2D1*omp*ra-7.0D0*om*r)*s12+1.464591887561523263
     2     D0*r23*(7.6254743439754925901D0*jpd*om*r*(1.2D1*omp*ra-7.0D0*
     3     om*r)-4.575284606385295554D1*om2*r*(2.0D0*jsdom*omp*ra-
     4     jsdr*r2))-1.3856406460551018348D1*jsdx*om2*r73)/(om3*s12)
        djdsaasaa = 9.0210979560879025705D-3*(1.7588825220236102883D0*js
     1     d*r16*s12-6.9816604245004959155D0*jpd*om*r43)*sm32/(om2*r52)
c...
        dldrara = 1.2037522124142963626D-3*rm193*(1.210203242253764276D0
     1     *r16*(-4.8037478207040109067D0*lpd*om*r*(1.44D2*om*oms*rb**2-
     2     2.88D2*omp**2*rb**2-3.12D2*om*omp*r*rb-9.1D1*om2*r2)-5.764497
     3     384844813088D1*om2*r*(1.2D1*omp*rb+7.0D0*om*r)*(2.0D0*lsomd*o
     4     mp*rb+epra*lsed*r2))*s12+1.464591887561523263D0*lsd*(1.2D1*om
     5     p*rb+7.0D0*om*r)**2*s+8.3073575249706722822D2*om4*r73*(4.0D0*
     6     lpom*oms*rb**2+4.0D0*lsom*omp**2*rb**2+4.0D0*epra*lseom*omp*r
     7     2*rb-4.0D0*lpom*omp*r*rb+epra**2*lse*r4+esra*lpe*r4))/om4
        dldrbrb = 1.2037522124142963626D-3*rm193*(1.210203242253764276D0
     1     *r16*(-4.8037478207040109067D0*lpd*om*r*(1.44D2*om*oms*ra**2-
     2     2.88D2*omp**2*ra**2+3.12D2*om*omp*r*ra-9.1D1*om2*r2)-5.764497
     3     384844813088D1*om2*r*(1.2D1*omp*ra-7.0D0*om*r)*(2.0D0*lsomd*o
     4     mp*ra-eprb*lsed*r2))*s12+1.464591887561523263D0*lsd*(1.
     5     2D1*omp*ra-7.0D0*om*r)**2*s+8.3073575249706722822D2*om4*r73*(
     6     4.0D0*lpom*oms*ra**2+4.0D0*lsom*omp**2*ra**2-4.0D0*eprb*lseom
     7     *omp*r2*ra+4.0D0*lpom*omp*r*ra+eprb**2*lse*r4+esrb*lpe*r4))/o
     8     m4
        dldrarb = -1.2037522124142963626D-3*rm193*(1.210203242253764276D
     1     0*r16*(lpd*(-4.8037478207040109067D0*om2*r*(1.44D2*oms*ra*rb+
     2     9.1D1*om*r2)-7.4938466002982570145D2*om2*omp*r2*(rb-ra)
     3     +1.3834793723627551411D3*om*omp**2*r*ra*rb)+lsed*(3.458698430
     4     9068878528D2*om2*omp*r3*(eprb*rb-epra*ra)+2.01757408469
     5     56845808D2*(eprb+epra)*om3*r4)+4.0351481693913691616D2*lsomd*
     6     om3*omp*r2*(rb-ra)-1.3834793723627551411D3*lsomd*om2*om
     7     p**2*r*ra*rb)*s12+1.464591887561523263D0*(-8.4D1*lsd*om*omp*r
     8     *(rb-ra)+1.44D2*lsd*omp**2*ra*rb-4.9D1*lsd*om2*r2)*s+r1
     9     3*(lpom*(1.6614715049941344564D3*om4*omp*r3*(rb-ra)+3.3
     :     229430099882689129D3*om4*oms*r2*ra*rb)-1.6614715049941344564D
     ;     3*lseom*om4*omp*r4*(eprb*rb-epra*ra)+3.3229430099882689
     <     129D3*lsom*om4*omp**2*r2*ra*rb-8.3073575249706722822D2*(epra*
     =     eprb*lse+esrarb*lpe)*om4*r6))/om4
        dldrasaa = -4.3703544910017702413D-3*(1.210203242253764276D0*lsd
     1     *(1.2D1*omp*rb+7.0D0*om*r)*s12+r16*(5.764497384844813088D1*lp
     2     d*om*omp*r*rb-4.8037478207040109067D0*om2*r*(1.2D1*lsomd*omp*
     3     rb+6.0D0*epra*lsed*r2-7.0D0*lpd*r)))/(om3*r133*s12)
        dldrbsaa = 4.3703544910017702413D-3*(1.210203242253764276D0*lsd*
     1     (1.2D1*omp*ra-7.0D0*om*r)*s12+r16*(5.764497384844813088D1*lpd
     2     *om*omp*r*ra-4.8037478207040109067D0*om2*r*(1.2D1*lsomd*omp*r
     3     a-6.0D0*eprb*lsed*r2+7.0D0*lpd*r)))/(om3*r133*s12)
        dldsaasaa = 9.0210979560879025705D-3*(1.7588825220236102883D0*ls
     1     d*r16*s12-6.9816604245004959155D0*lpd*om*r43)*sm32/(om2*r52)
c...
      endif
c...
c...  We assemble the output values.
c...
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
      dsab=dsaa+dsaa
c...
c...  Second derivatives of the functional.
c...
      drara = scnl*((dldra+djdra)*two+(dldrara+djdrara)*r)+dpw91lcdrara*
     1   scl
      drbrb = scnl*((dldrb+djdrb)*two+(dldrbrb+djdrbrb)*r)+dpw91lcdrbrb*
     1   scl
      drarb = ((dldrarb+djdrarb)*r+dldrb+dldra+djdrb+djdra)*scnl+dpw91lc
     1   drarb*scl
      drasaa = ((dldrasaa+djdrasaa)*r+dldsaa+djdsaa)*scnl
      drasbb = drasaa
      drasab=drasaa+drasaa
      drbsaa = ((dldrbsaa+djdrbsaa)*r+dldsaa+djdsaa)*scnl
      drbsbb = drbsaa
      drbsab=drbsaa+drbsaa
      dsaasaa = (dldsaasaa+djdsaasaa)*r*scnl
      dsaasbb = dsaasaa
      dsaasab=dsaasaa+dsaasaa
      dsbbsbb = dsaasaa
      dsbbsab = dsaasab
      dsabsab=dsaasab+dsaasab
c...
      return
      end
c=======================================================================
      subroutine pw91lcdu(ra,rb,ec,dra,drb,drara,drbrb,drarb)
      implicit none   !trying to avoid typing errors
c...
c...  MM (04/25/2003) generated with maxima, file pw91lcu.mc
c...
c...  Perdew and Wang local correlation functional
c...  *** open shell version
c...
c...  J.P. Perdew and Y. Wang, Phys. Rev. B 45, 13244 (1992)
c...
c...  Formulation taken from the Molpro manual (PW92C)
c...  http://www.molpro.net/current/molpro_manual
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dra     funct. deriv. w.r.t. beta density
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...
      real*8 ra,rb,ec,dra,drb,drara,drbrb,drarb
c...
c...  Numerical parameters of the functional.
c...  these instructions are no longer needed, but i leave them as
c...  comments for reference.
c...
c...  real*8 ct1,ct2,ct3,cu1,cu2,cu3,cv1,cv2,cv3,cw1,cw2,cw3,cx1,cx2,cx3
c...  real*8 cy1,cy2,cy3,c
c...  parameter (ct1=0.031091D0,  ct2=0.015545D0,   ct3=0.016887D0)
c...  parameter (cu1=0.21370D0,   cu2=0.20548D0,    cu3=0.11125D0)
c...  parameter (cv1=7.5957D0,    cv2=14.1189D0,    cv3=10.357D0)
c...  parameter (cw1=3.5876D0,    cw2=6.1977D0,     cw3=3.6231D0)
c...  parameter (cx1=1.6382D0,    cx2=3.3662D0,     cx3=0.88026D0)
c...  parameter (cy1=0.49294D0,   cy2=0.62517D0,    cy3=0.49671D0)
c...  parameter (  c=1.709921D0)
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
      real*8 third,tthird,etm
      parameter (third=0.33333333333333333333d0)
      parameter (tthird=0.66666666666666666667d0)
      parameter (etm=-11.0d0/3.0d0)
      real*8 r,r2,rm13,r13,r23,rm43,r43,rm113
      real*8 x,x2,x3,x12,x32,polx1,polx2,polx3
      real*8 z,z2,z3,z4
      real*8 zp1,zp113,zp123,zp143
      real*8 zm1,zm113,zm123,zm143
      real*8 eps11,eps12,eps13,eps1p1,eps1p2,eps1p3,eps1s1,eps1s2,eps1s3
      real*8 eps1,eps2,eps3,epsp1,epsp2,epsp3,epss1,epss2,epss3
      real*8 om,omp,oms
c...
c...  first we compute some powers of r that will be used later
c...
      r = rb+ra
      r2 = r*r
      r13 = r**THIRD
      rm13 = one/r13
      r23 = r13*r13
      r43 = r*r13
      rm43 = one/r43
      rm113 = r**etm
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.2035049089940001667D-1*rm13
      x2 = x*x
      x3 = x2*x
      x12 = sqrt(x)
      x32 = x**1.5d0
      z = min(max((ra-rb)/r,-one),one)
      z2 = z*z
      z3 = z2*z
      z4 = z3*z
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp123 = zp113*zp113
      zp143 = zp123*zp123
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm123 = zm113*zm113
      zm143 = zm123*zm123
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
c...  Second derivatives.
c...
      eps1s1 = 2.5D-1*(((4.8350235714652928D-1*polx1+3.8878374976D0)*x12
     1   +2.4102534633346176D0*polx1+1.9380823392D1)*x3+((6.522684333513
     2   0592D0*polx1+5.2448837864D1)*x12+1.24960078095074112D1*polx1+1.
     3   00480270704D2)*x2+x*((-2.45213594528D-1*polx1**2+1.174399527015
     4   1344D1*polx1+1.2614314948D2)*x12-3.0559670844D-1*polx1**2+8.641
     5   111010433984D0*polx1+1.0900133328D2)+(7.175069274860172D0*polx1
     6   +5.769465849D1)*x12+4.7231125998D-1*polx1**2+7.5957D0*polx1)/(p
     7   olx1**2*(6.21814D-2*polx1+1.0D0)**2*x32)
      epss1 = -6.21814D-2*eps1s1*(2.137D-1*x+1.0D0)-2.657633036D-2*eps1p
     1   1
      eps1s2 = 2.5D-1*(((3.8884519551267936D-1*polx2+6.2534004624D0)*x12
     1   +3.1405794355170144D0*polx2+5.0506734096D1)*x3+((1.019623111180
     2   63056D1*polx2+1.63975579704D2)*x12+1.99581198701375376D1*polx2+
     3   3.20966074584D2)*x2+x*((-1.55495783352D-1*polx2**2+2.2284347476
     4   226736D1*polx2+4.3880818824D2)*x12-3.1397254302D-1*polx2**2+1.1
     5   666060634498168D1*polx2+3.5001882612D2)+(1.2395447788389894D1*p
     6   olx2+1.9934333721D2)*x12+4.3896648423D-1*polx2**2+1.41189D1*pol
     7   x2)/(polx2**2*(3.10907D-2*polx2+1.0D0)**2*x32)
      epss2 = -3.10907D-2*eps1s2*(2.0548D-1*x+1.0D0)-1.2777034072D-2*eps
     1   1p2
      eps1s3 = 2.5D-1*(((2.6664639260763456D-1*polx3+3.9475331856D0)*x12
     1   +7.0881848631031104D-1*polx3+1.04936146704D1)*x3+((2.4160289813
     2   8056144D0*polx3+3.57677990244D1)*x12+5.3650791669064032D0*polx3
     3   +7.9426643832D1)*x2+x*((-1.34206273584D-1*polx3**2+3.2679989928
     4   92736D0*polx3+1.0720853136D2)*x12-8.9189175564D-2*polx3**2+7.49
     5   796526365168D0*polx3+1.500977868D2)+(7.2456587380724D0*polx3+1.
     6   07267449D2)*x12+3.497952466D-1*polx3**2+1.0357D1*polx3)/(polx3*
     7   *2*(3.37738D-2*polx3+1.0D0)**2*x32)
      epss3 = -3.37738D-2*eps1s3*(1.1125D-1*x+1.0D0)-7.5146705D-3*eps1p3
      if(zm123*zp123.gt.epsi)then
        oms = 8.5496046708068280878D-1*(zp123+zm123)/(zm123*zp123)
      else
        oms=zero
      endif
      drara = 1.2021948647324711061D-2*rm113*(1.464591887561523263D0*eps
     1   p3*r13*(-2.7473141821279964827D1*r*rb*(omp*z4+4.0D0*om*z3-
     2   omp)-4.5788569702133274712D0*om*r2*(z4-1.0D0))+1.4645918875615
     3   23263D0*epsp1*r13*(4.6976902136184858738D1*r*rb*(omp*z4+4.0D0*o
     4   m*z3)+7.8294836893641431229D0*r2*(om*z4-1.0D0))+1.9458487368457
     5   129358D2*eps3*r23*rb**2*(oms*z4+8.0D0*omp*z3+1.2D1*om*z2-
     6   oms)+3.3272476179559583089D2*eps2*r23*rb**2*(oms*z4+8.0D0*omp*z
     7   3+1.2D1*om*z2)-3.3272476179559583089D2*eps1*r23*rb**2*(oms*z4+8
     8   .0D0*omp*z3+1.2D1*om*z2)-6.8801989771427936613D1*epsp2*r43*rb*(
     9   omp*z4+4.0D0*om*z3)-3.5567790107967349354D0*epss1*r2*(om*z4-1.0
     :   D0)+2.0800838230519041145D0*epss3*om*r2*(z4-1.0D0)-2.4661328275
     ;   096140485D0*om*(4.6497894060385059023D0*epsp2*r13-1.44224957030
     <   74083823D0*epss2)*r2*z4)
      drbrb = 1.2021948647324711061D-2*rm113*(1.464591887561523263D0*eps
     1   p3*r13*(2.7473141821279964827D1*r*ra*(omp*z4+4.0D0*om*z3-
     2   omp)-4.5788569702133274712D0*om*r2*(z4-1.0D0))+1.46459188756152
     3   3263D0*epsp1*r13*(7.8294836893641431229D0*r2*(om*z4-1.0D0)-4.69
     4   76902136184858738D1*r*ra*(omp*z4+4.0D0*om*z3))+1.94584873684571
     5   29358D2*eps3*r23*ra**2*(oms*z4+8.0D0*omp*z3+1.2D1*om*z2-o
     6   ms)+3.3272476179559583089D2*eps2*r23*ra**2*(oms*z4+8.0D0*omp*z3
     7   +1.2D1*om*z2)-3.3272476179559583089D2*eps1*r23*ra**2*(oms*z4+8.
     8   0D0*omp*z3+1.2D1*om*z2)+6.8801989771427936613D1*epsp2*r43*ra*(o
     9   mp*z4+4.0D0*om*z3)-3.5567790107967349354D0*epss1*r2*(om*z4-1.0D
     :   0)+2.0800838230519041145D0*epss3*om*r2*(z4-1.0D0)-2.46613282750
     ;   96140485D0*om*(4.6497894060385059023D0*epsp2*r13-1.442249570307
     <   4083823D0*epss2)*r2*z4)
      drarb = -1.2021948647324711061D-2*rm113*(1.464591887561523263D0*ep
     1   sp3*r13*(1.3736570910639982414D1*r*(rb-1.0D0*ra)*(omp*z4+4.0D0*
     2   om*z3-omp)+4.5788569702133274712D0*om*r2*(z4-1.0D0))+1.46
     3   4591887561523263D0*epsp1*r13*(-2.3488451068092429369D1*r*(rb-
     4   ra)*(omp*z4+4.0D0*om*z3)-7.8294836893641431229D0*r2*(om*z4-
     5   1.0D0))+1.9458487368457129358D2*eps3*r23*ra*rb*(oms*z4+8.0D0*om
     6   p*z3+1.2D1*om*z2-oms)+3.3272476179559583089D2*eps2*r23*ra
     7   *rb*(oms*z4+8.0D0*omp*z3+1.2D1*om*z2)-3.3272476179559583089D2*e
     8   ps1*r23*ra*rb*(oms*z4+8.0D0*omp*z3+1.2D1*om*z2)+3.4400994885713
     9   968307D1*epsp2*r43*(rb-ra)*(omp*z4+4.0D0*om*z3)+3.5567790
     :   107967349354D0*epss1*r2*(om*z4-1.0D0)-2.0800838230519041145D0*e
     ;   pss3*om*r2*(z4-1.0D0)+2.4661328275096140485D0*om*(4.64978940603
     <   85059023D0*epsp2*r13-1.4422495703074083823D0*epss2)*r2*z4)
c...
      return
      end
c=======================================================================
      subroutine pw91xdu(ra,rb,saa,sbb,exc,dra,drb,dsaa,dsbb,drara,drbrb
     $  ,drasaa,drbsbb,dsaasaa,dsbbsbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...
c...  Note that this functional does not have mixed alpha-beta
c...  derivatives
c...
      real*8 ra,rb,saa,sbb,exc,eca,esla,ecb,eslb,dra,d1ra,dslra,drb,d1rb
     $  ,dslrb,dsaa,dsbb
      real*8 drara,d1rara,dslrara,drasaa,dsaasaa
      real*8 drbrb,d1rbrb,dslrbrb,drbsbb,dsbbsbb
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 rbm1,rb13,rb43,rbm43,sbb12,sbbm12
      real*8 ram23,ra53,ram73,ram113,saam32
      real*8 rbm23,rb53,rbm73,rbm113,sbbm32
      real*8 essea,essea2,essea3,essea4
      real*8 esseb,esseb2,esseb3,esseb4
      real*8 cessea1,cessea12,cesseam12
      real*8 cesseb1,cesseb12,cessebm12
      real*8 cessea32,cesseam32
      real*8 cesseb32,cessebm32
      real*8 nfa,dfa,nfap,dfap,dfa2,dfam2
      real*8 nfb,dfb,nfbp,dfbp,dfb2,dfbm2
      real*8 nfas,dfas,dfam3
      real*8 nfbs,dfbs,dfbm3
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
      ram23 = SQRT(ram43)
      ra53 = ra13*ra43
      ram73 = ram1*ram43
      ram113 = ram43*ram73
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
c...
c...  Second derivatives
c...
      cessea32 = cessea1*cessea12
      cesseam32 = one/cessea32
      dfam3 = dfam2/dfa
      nfas = cesseam32*(-3.016D-1*cessea32*EXP(-1.0D2*essea2)*(2.0D4*ess
     $  ea4-5.0D2*essea2+1.0D0)-9.306806274223041D1*essea2+5.486D-1*cess
     $  ea32+3.06289124D0*cessea1)
      dfas = -9.306806274223041D1*cesseam32*essea2+4.8D-2*essea2+3.06289
     $  124D0*cesseam12
c...
c...  Second derivatives of the functional
c...
      d1rara = -1.361074440023947D-2*dfam3*ram113*(-3.897777089720754D0*
     $  dfa*(dfa*nfap-dfap*nfa)*ra53*saa12+2.0D0*(dfa2*nfas-2.0D0*
     $  dfa*dfap*nfap-dfa*dfas*nfa+2.0D0*dfap**2*nfa)*ra13*saa+3.0
     $  38533248230398D1*dfa2*nfa*ra**3)
      dslrara = -4.135669939329334D-1*ram23
      drasaa = 1.02080583001796D-2*dfam3*(dfa2*nfas-2.0D0*dfa*dfap*nfap-
     $  dfa*dfas*nfa+2.0D0*dfap**2*nfa)*ram73
      if(saa.gt.epsi)then
      saam32 = saam12/saa
      dsaasaa = -3.828021862567351D-3*dfam3*ram43*((dfa2*nfas-2.0D0*dfa*
     $  dfap*nfap-dfa*dfas*nfa+2.0D0*dfap**2*nfa)*saa12-7.79555417
     $  9441508D0*dfa*(dfa*nfap-dfap*nfa)*ra43)*saam32
      else
      dsaasaa = ZERO
      endif
      else
      eca = ZERO
      esla = ZERO
      d1ra = ZERO
      dslra = ZERO
      dsaa = ZERO
      d1rara = ZERO
      dslrara = ZERO
      drasaa = ZERO
      dsaasaa = ZERO
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
      rbm23 = SQRT(rbm43)
      rb53 = rb13*rb43
      rbm73 = rbm1*rbm43
      rbm113 = rbm43*rbm73
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
c...
c...  Second derivatives
c...
      cesseb32 = cesseb1*cesseb12
      cessebm32 = one/cesseb32
      dfbm3 = dfbm2/dfb
      nfbs = cessebm32*(-3.016D-1*cesseb32*EXP(-1.0D2*esseb2)*(2.0D4*ess
     $  eb4-5.0D2*esseb2+1.0D0)-9.306806274223041D1*esseb2+5.486D-1*cess
     $  eb32+3.06289124D0*cesseb1)
      dfbs = -9.306806274223041D1*cessebm32*esseb2+4.8D-2*esseb2+3.06289
     $  124D0*cessebm12
c...
c...  Second derivatives of the functional
c...
      d1rbrb = -1.361074440023947D-2*dfbm3*rbm113*(-3.897777089720754D0*
     $  dfb*(dfb*nfbp-dfbp*nfb)*rb53*sbb12+2.0D0*(dfb2*nfbs-2.0D0*
     $  dfb*dfbp*nfbp-dfb*dfbs*nfb+2.0D0*dfbp**2*nfb)*rb13*sbb+3.0
     $  38533248230398D1*dfb2*nfb*rb**3)
      dslrbrb = -4.135669939329334D-1*rbm23
      drbsbb = 1.02080583001796D-2*dfbm3*(dfb2*nfbs-2.0D0*dfb*dfbp*nfbp-
     $  dfb*dfbs*nfb+2.0D0*dfbp**2*nfb)*rbm73
      if(sbb.gt.epsi)then
      sbbm32 = sbbm12/sbb
      dsbbsbb = -3.828021862567351D-3*dfbm3*rbm43*((dfb2*nfbs-2.0D0*dfb*
     $  dfbp*nfbp-dfb*dfbs*nfb+2.0D0*dfbp**2*nfb)*sbb12-7.79555417
     $  9441508D0*dfb*(dfb*nfbp-dfbp*nfb)*rb43)*sbbm32
      else
      dsbbsbb = ZERO
      endif
      else
      ecb = ZERO
      eslb = ZERO
      d1rb = ZERO
      dslrb = ZERO
      dsbb = ZERO
      d1rbrb = ZERO
      dslrbrb = ZERO
      drbsbb = ZERO
      dsbbsbb = ZERO
      endif
c...
c...  Assemble the output values, subtracting the Slater term
c...
      exc = (-eslb-esla+ecb+eca)/(rb+ra)
      dra = d1ra-dslra
      drb = d1rb-dslrb
      drara = d1rara-dslrara
      drbrb = d1rbrb-dslrbrb
c...
      return
      end
c=======================================================================
      subroutine pw91xd(ra,saa,exc,dra,dsaa,drara,drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,exc,dra,dsaa,eca,esla,d1ra,dslra
      real*8 drara,d1rara,dslrara,drasaa,dsaasaa
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 ram23,ra53,ram73,ram113,saam32
      real*8 essea,essea2,essea3,essea4
      real*8 cessea1,cessea12,cesseam12
      real*8 cessea32,cesseam32
      real*8 nf,df,nfp,dfp,df2,dfm2
      real*8 nfs,dfs,dfm3
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
      ram23 = SQRT(ram43)
      ra53 = ra13*ra43
      ram73 = ram1*ram43
      ram113 = ram43*ram73
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
c...  Second derivatives
c...
      cessea32 = cessea1*cessea12
      cesseam32 = one/cessea32
      dfm3 = dfm2/df
      nfs = cesseam32*(-3.016D-1*cessea32*EXP(-1.0D2*essea2)*(2.0D4*esse
     $  a4-5.0D2*essea2+1.0D0)-9.306806274223041D1*essea2+5.486D-1*cesse
     $  a32+3.06289124D0*cessea1)
      dfs = -9.306806274223041D1*cesseam32*essea2+4.8D-2*essea2+3.062891
     $  24D0*cesseam12
c...
c...  Second derivatives of the functional
c...
      d1rara = -1.361074440023947D-2*dfm3*ram113*(-3.897777089720754D0*d
     $  f*(df*nfp-dfp*nf)*ra53*saa12+2.0D0*(df2*nfs-2.0D0*df*dfp*n
     $  fp-df*dfs*nf+2.0D0*dfp**2*nf)*ra13*saa+3.038533248230398D1
     $  *df2*nf*ra**3)
      dslrara = -4.135669939329334D-1*ram23
      drasaa = 1.02080583001796D-2*dfm3*(df2*nfs-2.0D0*df*dfp*nfp-
     $  df*dfs*nf+2.0D0*dfp**2*nf)*ram73
      if(saa.gt.epsi)then
      saam32 = saam12/saa
      dsaasaa = -3.828021862567351D-3*dfm3*ram43*((df2*nfs-2.0D0*df*dfp*
     $  nfp-df*dfs*nf+2.0D0*dfp**2*nf)*saa12-7.795554179441508D0*d
     $  f*(df*nfp-dfp*nf)*ra43)*saam32
      else
      dsaasaa = ZERO
      endif
c...
c...  Assemble the output values, subtracting the Slater term
c...
      exc = (eca-esla)/ra
      dra = d1ra-dslra
      drara = d1rara-dslrara
c...
      return
      end
c=======================================================================
      subroutine pbecdu(ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab,drar
     $  a,drbrb,drarb,drasaa,drasbb,drasab,drbsaa,drbsbb,drbsab,dsaasaa,
     $  dsbbsbb,dsabsab,dsaasbb,dsaasab,dsbbsab)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drasab  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  drbsab  funct. mixed 2nd deriv. w.r.t. beta den. and alpha beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...  dsaasab funct. mixed 2nd deriv. w.r.t. alpha and alpha beta grad.
c...  dsbbsab funct. mixed 2nd deriv. w.r.t. beta and alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,sab,ec,dra,drb,dsaa,dsbb,dsab
      real* 8 drara,drbrb,drarb,drasaa,drasbb,drasab
      real*8 drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb
      real*8 dsabsab,dsaasbb,dsaasab,dsbbsab
      real*8 pw91lc,dpw91lcdra,dpw91lcdrb,dpw91lcdrara,dpw91lcdrbrb,dpw9
     $  1lcdrarb
      real*8 r,r2,r3,r4,r5,r6,r7,r12,r32,r52,r72,r16,rm13,r13,r23,r56,rm
     $  76,r76,rm43,r43,r53,r73,r83,r136,r113,r133,rm193,rm196,rm296
      real*8 s,s12,sm32
      real*8 z
      real*8 zp1,zp113,zp123,zp143
      real*8 zm1,zm113,zm123,zm143
      real*8 d,d2,d4
      real*8 numlog,denlog
      real*8 om,omp,oms,om2,om3,om4,om6,om7,om8
      real*8 e,epra,eprb,esra,esrb,esrarb
      real*8 expa,a,ape,apom,ase,asom,aseom
      real*8 loga,logape,logapd,logapom,logase,logasd,logasom,logased,lo
     $  gaseom,logasomd
      real*8 h,hpe,hpom,hpd,hse,hsom,hsd,hseom,hsed,hsomd
      real*8 dhdra,dhdrb,dhdsaa,dhdrara,dhdrbrb,dhdrarb,dhdrasaa,dhdrbsa
     $  a,dhdsaasaa
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
      call pw91lcdu(ra,rb,pw91lc,dpw91lcdra,dpw91lcdrb,dpw91lcdrara,dpw9
     $  1lcdrbrb,dpw91lcdrarb)
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
      dhdrara = ZERO
      dhdrbrb = ZERO
      dhdrarb = ZERO
      dhdrasaa = ZERO
      dhdrbsaa = ZERO
      dhdsaasaa = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
      r2 = r**2
      r3 = r*r2
      r4 = r*r3
      r5 = r*r4
      r6 = r*r5
      r7 = r*r6
      r12 = SQRT(r)
      r32 = r*r12
      r52 = r12*r2
      r72 = r12*r3
      r16 = r**SIXTH
      r13 = r16**2
      rm13 = one/r13
      r23 = r13**2
      r56 = r16*r23
      r76 = r*r16
      rm76 = one/r76
      r43 = r*r13
      rm43 = one/r43
      r53 = r13*r43
      r73 = r13*r2
      r83 = r13*r73
      r136 = r16*r2
      r113 = r23*r3
      r133 = r13*r4
      rm193 = one/(r13*r6)
      rm196 = one/(r16*r3)
      rm296 = one/(r4*r56)
c...
      s12 = SQRT(s)
      sm32 = one/(s*s12)
c...
c...  the auxiliary variables z and some of its powers.
c...
      z = MIN(MAX((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp123 = zp113**2
      zp143 = zp123**2
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm123 = zm113**2
      zm143 = zm123**2
c...
c...  the variables omega and d.
c...
      om = 5.0D-1*(zp123+zm123)
      om2 = om**2
      om3 = om*om2
      om4 = om*om3
      om6 = om2*om4
      om7 = om*om6
      om8 = om*om7
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
c...  second derivative of om
c...
      if(zm143*zp143.gt.epsi)then
      oms = -1.111111111111111D-1*(z*zp113+zp113-z*zm113+zm113)/(z
     $  m143*zp143)
      else
      oms = ZERO
      endif
c...
c...  second derivatives of e
c...
      esra = -(2.0D0*epra-dpw91lcdrara)/r
      esrb = -(2.0D0*eprb-dpw91lcdrbrb)/r
      esrarb = -(eprb+epra-dpw91lcdrarb)/r
c...
c...  Second derivatives of h w.r.t. auxiliary variables
c...
      ase = 2.220227433125678D3*expa*(expa+1.0D0)/((expa-1.0D0)**3*om6)
      asom = 6.660682299377035D3*e*expa*(expa*(1.243627634786196D-1*om3+
     $  3.0D0*e)-1.243627634786196D-1*om3+3.0D0*e)/((expa-1.0D0)**3*om8)
      aseom = -6.660682299377035D3*expa*(expa*(3.10906908696549D-2*om3+e
     $  )-3.10906908696549D-2*om3+e)/((expa-1.0D0)**3*om7)
      logase = -6.6725D-2*d2*(6.6725D-2*(2.0D0*a*ase*d2+2.0D0*ape**2*d2+
     $  ase)*denlog*numlog**2-6.6725D-2*ape**2*d2*(2.0D0*a*d2+1.0D0)**2*
     $  numlog**2+(6.218138173930979D-2*a*ase*d2-6.6725D-2*ase*d2+6.2181
     $  38173930979D-2*ape**2*d2+3.10906908696549D-2*ase)*denlog**2*numl
     $  og-6.218138173930979D-2*ape**2*d2*(2.0D0*a*d2+1.0D0)**2*denlog*n
     $  umlog-3.10906908696549D-2*ase*d2*denlog**3+ape**2*(1.24362763478
     $  6196D-1*a*d2+6.6725D-2*d2+6.218138173930979D-2)*d4*denlog**2)/(d
     $  enlog**2*(6.6725D-2*numlog+3.10906908696549D-2*denlog)**2)
      logasom = -6.6725D-2*d2*(6.6725D-2*(2.0D0*a*asom*d2+2.0D0*apom**2*
     $  d2+asom)*denlog*numlog**2-6.6725D-2*apom**2*d2*(2.0D0*a*d2+1.0D0
     $  )**2*numlog**2+(6.218138173930979D-2*a*asom*d2-6.6725D-2*asom*d2
     $  +6.218138173930979D-2*apom**2*d2+3.10906908696549D-2*asom)*denlo
     $  g**2*numlog-6.218138173930979D-2*apom**2*d2*(2.0D0*a*d2+1.0D0)**
     $  2*denlog*numlog-3.10906908696549D-2*asom*d2*denlog**3+apom**2*(1
     $  .243627634786196D-1*a*d2+6.6725D-2*d2+6.218138173930979D-2)*d4*d
     $  enlog**2)/(denlog**2*(6.6725D-2*numlog+3.10906908696549D-2*denlo
     $  g)**2)
      logasd = -1.3345D-1*(a*numlog-1.0D0*denlog)*(denlog*(6.6725D-2*(6.
     $  0D0*a*d2+1.0D0)*numlog-2.0D0*(6.218138173930979D-2*a+6.6725D-2)*
     $  d2*(2.0D0*a*d2+1.0D0)**2)-1.3345D-1*a*d2*(2.0D0*a*d2+1.0D0)**2*n
     $  umlog+3.10906908696549D-2*(6.0D0*a*d2+1.0D0)*denlog**2)/(denlog*
     $  *2*(6.6725D-2*numlog+3.10906908696549D-2*denlog)**2)
      logaseom = -6.6725D-2*d2*(denlog*((6.6725D-2*aseom*(2.0D0*a*d2+1.0
     $  D0)+1.3345D-1*ape*apom*d2)*numlog**2-6.218138173930979D-2*ape*ap
     $  om*d2*(2.0D0*a*d2+1.0D0)**2*numlog)-6.6725D-2*ape*apom*d2*(2.0D0
     $  *a*d2+1.0D0)**2*numlog**2+denlog**2*((aseom*(6.218138173930979D-
     $  2*a*d2-6.6725D-2*d2+3.10906908696549D-2)+6.218138173930979D-2*ap
     $  e*apom*d2)*numlog+ape*apom*(1.243627634786196D-1*a*d2+6.6725D-2*
     $  d2+6.218138173930979D-2)*d4)-3.10906908696549D-2*aseom*d2*denlog
     $  **3)/(denlog**2*(6.6725D-2*numlog+3.10906908696549D-2*denlog)**2
     $  )
      logased = -1.3345D-1*ape*d*(denlog*(6.6725D-2*(4.0D0*a*d2+1.0D0)*n
     $  umlog**2-6.218138173930979D-2*a*d2*(2.0D0*a*d2+1.0D0)**2*numlog)
     $  -6.6725D-2*a*d2*(2.0D0*a*d2+1.0D0)**2*numlog**2+denlog**2*((1.24
     $  3627634786196D-1*a*d2-1.3345D-1*d2+3.10906908696549D-2)*numlog+d
     $  2*(9.327207260896469D-2*a*d2+6.6725D-2*d2+3.10906908696549D-2)*(
     $  2.0D0*a*d2+1.0D0))-6.218138173930979D-2*d2*denlog**3)/(denlog**2
     $  *(6.6725D-2*numlog+3.10906908696549D-2*denlog)**2)
      logasomd = -1.3345D-1*apom*d*(denlog*(6.6725D-2*(4.0D0*a*d2+1.0D0)
     $  *numlog**2-6.218138173930979D-2*a*d2*(2.0D0*a*d2+1.0D0)**2*numlo
     $  g)-6.6725D-2*a*d2*(2.0D0*a*d2+1.0D0)**2*numlog**2+denlog**2*((1.
     $  243627634786196D-1*a*d2-1.3345D-1*d2+3.10906908696549D-2)*numlog
     $  +d2*(9.327207260896469D-2*a*d2+6.6725D-2*d2+3.10906908696549D-2)
     $  *(2.0D0*a*d2+1.0D0))-6.218138173930979D-2*d2*denlog**3)/(denlog*
     $  *2*(6.6725D-2*numlog+3.10906908696549D-2*denlog)**2)
      hse = 3.10906908696549D-2*logase*om3
      hsom = 3.10906908696549D-2*om*(logasom*om2+6.0D0*logapom*om+6.0D0*
     $  loga)
      hsd = 3.10906908696549D-2*logasd*om3
      hseom = 3.10906908696549D-2*(logaseom*om+3.0D0*logape)*om2
      hsed = 3.10906908696549D-2*logased*om3
      hsomd = 3.10906908696549D-2*(logasomd*om+3.0D0*logapd)*om2
c...
c...  Total second derivatives of the non local term h
c...
c...
      dhdrara = 1.203752212414296D-3*rm193*(1.210203242253764D0*r16*(-4.
     $  803747820704011D0*hpd*om*r*(1.44D2*om*oms*rb**2-2.88D2*omp**2*rb
     $  **2-3.12D2*om*omp*r*rb-9.1D1*om2*r2)-5.764497384844813D1*om2*r*(
     $  1.2D1*omp*rb+7.0D0*om*r)*(2.0D0*hsomd*omp*rb+epra*hsed*r2))*s12+
     $  1.464591887561523D0*hsd*(1.2D1*omp*rb+7.0D0*om*r)**2*s+8.3073575
     $  24970672D2*om4*r73*(4.0D0*hpom*oms*rb**2+4.0D0*hsom*omp**2*rb**2
     $  +4.0D0*epra*hseom*omp*r2*rb-4.0D0*hpom*omp*r*rb+epra**2*hse*r4+e
     $  sra*hpe*r4))/om4
      dhdrbrb = 1.203752212414296D-3*rm193*(1.210203242253764D0*r16*(-4.
     $  803747820704011D0*hpd*om*r*(1.44D2*om*oms*ra**2-2.88D2*omp**2*ra
     $  **2+3.12D2*om*omp*r*ra-9.1D1*om2*r2)-5.764497384844813D1*om2*r*(
     $  1.2D1*omp*ra-7.0D0*om*r)*(2.0D0*hsomd*omp*ra-1.0D0*eprb*hsed*r2)
     $  )*s12+1.464591887561523D0*hsd*(1.2D1*omp*ra-7.0D0*om*r)**2*s+8.3
     $  07357524970672D2*om4*r73*(4.0D0*hpom*oms*ra**2+4.0D0*hsom*omp**2
     $  *ra**2-4.0D0*eprb*hseom*omp*r2*ra+4.0D0*hpom*omp*r*ra+eprb**2*hs
     $  e*r4+esrb*hpe*r4))/om4
      dhdrarb = -1.203752212414296D-3*rm193*(1.210203242253764D0*r16*(hp
     $  d*(-4.803747820704011D0*om2*r*(1.44D2*oms*ra*rb+9.1D1*om*r2)-7.4
     $  93846600298257D2*om2*omp*r2*(rb-ra)+1.383479372362755D3*om
     $  *omp**2*r*ra*rb)+hsed*(3.458698430906888D2*om2*omp*r3*(eprb*rb-
     $  epra*ra)+2.017574084695685D2*(eprb+epra)*om3*r4)+4.03514816
     $  9391369D2*hsomd*om3*omp*r2*(rb-ra)-1.383479372362755D3*hso
     $  md*om2*omp**2*r*ra*rb)*s12+1.464591887561523D0*(-8.4D1*hsd*om*om
     $  p*r*(rb-1.0D0*ra)+1.44D2*hsd*omp**2*ra*rb-4.9D1*hsd*om2*r2)*s+r1
     $  3*(hpom*(1.661471504994134D3*om4*omp*r3*(rb-ra)+3.32294300
     $  9988269D3*om4*oms*r2*ra*rb)-1.661471504994134D3*hseom*om4*omp*r4
     $  *(eprb*rb-epra*ra)+3.322943009988269D3*hsom*om4*omp**2*r2*
     $  ra*rb-8.307357524970672D2*(epra*eprb*hse+esrarb*hpe)*om4*r6))/om
     $  4
      dhdrasaa = -4.37035449100177D-3*(1.210203242253764D0*hsd*(1.2D1*om
     $  p*rb+7.0D0*om*r)*s12+r16*(5.764497384844813D1*hpd*om*omp*r*rb-4.
     $  803747820704011D0*om2*r*(1.2D1*hsomd*omp*rb+6.0D0*epra*hsed*r2-7
     $  .0D0*hpd*r)))/(om3*r133*s12)
      dhdrbsaa = 4.37035449100177D-3*(1.210203242253764D0*hsd*(1.2D1*omp
     $  *ra-7.0D0*om*r)*s12+r16*(5.764497384844813D1*hpd*om*omp*r*ra-4.8
     $  03747820704011D0*om2*r*(1.2D1*hsomd*omp*ra-6.0D0*eprb*hsed*r2+7.
     $  0D0*hpd*r)))/(om3*r133*s12)
      dhdsaasaa = 9.021097956087903D-3*(1.75888252202361D0*hsd*r16*s12-6
     $  .981660424500496D0*hpd*om*r43)*sm32/(om2*r52)
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
c...  Second derivatives of the functional.
c...
      drara = dhdra*two+dhdrara*r+dpw91lcdrara
      drbrb = dhdrb*two+dhdrbrb*r+dpw91lcdrbrb
      drarb = dhdrarb*r+dpw91lcdrarb+dhdrb+dhdra
      drasaa = dhdrasaa*r+dhdsaa
      drasbb = drasaa
      drasab = 2.0D0*drasaa
      drbsaa = dhdrbsaa*r+dhdsaa
      drbsbb = drbsaa
      drbsab = 2.0D0*drbsaa
      dsaasaa = dhdsaasaa*r
      dsaasbb = dsaasaa
      dsaasab = 2.0D0*dsaasaa
      dsbbsbb = dsaasaa
      dsbbsab = dsaasab
      dsabsab = 2.0D0*dsaasab
c...
      return
      end
c=======================================================================
      subroutine pbexdu(ra,rb,saa,sbb,exc,dra,drb,dsaa,dsbb,drara,drbrb,
     $  drasaa,drbsbb,dsaasaa,dsbbsbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...
c...  Note that this functional does not have mixed alpha-beta
c...  derivatives
c...
      real*8 ra,rb,saa,sbb,exc,exa,exb,dra,drb,dsaa,dsbb
      real*8 drara,drasaa,dsaasaa
      real*8 drbrb,drbsbb,dsbbsbb
      real*8 ram1,ra13,ra43,ram43,saa12,saam12
      real*8 rbm1,rb13,rb43,rbm43,sbb12,sbbm12
      real*8 dfam1,dfam2
      real*8 dfbm1,dfbm2
      real*8 ram2,ram23,ram73,ram103,saam1,saam32
      real*8 rbm2,rbm23,rbm73,rbm103,sbbm1,sbbm32
      real*8 dfam3,dfbm3
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
      ram2 = ram1**2
      ram23 = SQRT(ram43)
      ram73 = ram1*ram43
      ram103 = ram1*ram73
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
c...
c...  Second derivatives
c...
      dfam3 = dfam1*dfam2
      drara = 2.329141886786433D-2*dfam2*essea*ram2*saa12+1.305208695601
     $  004D-2*dfam3*essea2*ram103*saa-1.195112923686099D-2*dfam2*ram103
     $  *saa-4.135669939329334D-1*(8.04D-1-8.04D-1*dfam1)*ram23
      drasaa = 4.481673463822872D-3*dfam2*ram73-4.894532608503765D-3*dfa
     $  m3*essea2*ram73
      if(saa.gt.epsi)then
      saam1 = one/saa
      saam32 = saam1*saam12
      dsaasaa = 1.310142311317369D-2*dfam2*essea*saam32+1.83544972818891
     $  2D-3*dfam3*essea2*ram43*saam1-1.680627548933577D-3*dfam2*ram43*s
     $  aam1
      else
      dsaasaa = ZERO
      endif
      else
      exa = ZERO
      dra = ZERO
      dsaa = ZERO
      drara = ZERO
      drasaa = ZERO
      dsaasaa = ZERO
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
      rbm2 = rbm1**2
      rbm23 = SQRT(rbm43)
      rbm73 = rbm1*rbm43
      rbm103 = rbm1*rbm73
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
c...
c...  Second derivatives
c...
      dfbm3 = dfbm1*dfbm2
      drbrb = 2.329141886786433D-2*dfbm2*esseb*rbm2*sbb12+1.305208695601
     $  004D-2*dfbm3*esseb2*rbm103*sbb-1.195112923686099D-2*dfbm2*rbm103
     $  *sbb-4.135669939329334D-1*(8.04D-1-8.04D-1*dfbm1)*rbm23
      drbsbb = 4.481673463822872D-3*dfbm2*rbm73-4.894532608503765D-3*dfb
     $  m3*esseb2*rbm73
      if(sbb.gt.epsi)then
      sbbm1 = one/sbb
      sbbm32 = sbbm1*sbbm12
      dsbbsbb = 1.310142311317369D-2*dfbm2*esseb*sbbm32+1.83544972818891
     $  2D-3*dfbm3*esseb2*rbm43*sbbm1-1.680627548933577D-3*dfbm2*rbm43*s
     $  bbm1
      else
      dsbbsbb = ZERO
      endif
      else
      exb = ZERO
      drb = ZERO
      dsbb = ZERO
      drbrb = ZERO
      drbsbb = ZERO
      dsbbsbb = ZERO
      endif
c...
c...  functional output value
c...
      exc = (exb+exa)/(rb+ra)
c...
      return
      end

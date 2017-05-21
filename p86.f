      subroutine p86nlcr(ra,saa,scl,scnl,ec,dra,dsaa)
      implicit none
c...  MM (07/10/2003) generated with maxima, file p86nlcr.mc
c...
c...           Perdew 86 correlation functional
c...                   non local part
c...
c...     J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...                   closed shell version
c...
c...     Formulation taken from the Molpro manual (P86)
c...     http://www.molpro.net/current/molpro_manual
c...  (note that the local part is different here than in Molpro)
c...
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra    alpha density
c...  saa   alpha density gradient invariant
c...  scl   scaling factor for local part
c...  scnl  scaling factor for non local part
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,scl,scnl,ec,dra,dsaa
      real*8 eloc,plocra
      real*8 enloc,pnlocra,pnlocsaa
      real*8 r,rm13,r13,rm73,rm76,rm43,rm83,rm136,rm113
      real*8 s,s12,sm12
      real*8 x,x2,x3
      real*8 numc,numcp,dencm1,dencm2,dencp,c,cm1,cm2,cp
      real*8 phi,expphi,phipc,phipr,phips
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,thirdm,tthird,sixth,etm,th,ssm
      parameter (third=1.0d0/3.0d0,thirdm=-1.0d0/3.0d0,tthird=2.0d0/3.0d
     $  0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0,ssm=-7.0d0/6.0d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      if(scl.ne.zero)then
      call pz81lcr(ra,eloc,plocra)
      else
      eloc = ZERO
      plocra = ZERO
      endif
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      r = 2.0D0*ra
      s = 4.0D0*saa
      if(s.lt.epsi.or.scnl.eq.zero)then
      enloc = ZERO
      pnlocra = ZERO
      pnlocsaa= ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
      r13 = r**THIRD
      rm13 = one/r13
      rm76 = r**ssm
      rm43 = rm13/r
      rm73 = rm43/r
      rm83 = rm13*rm73
      rm113 = r**etm
      rm136 = rm76/r
c...
      s12 = SQRT(s)
      sm12 = one/s12
c...
c...  the auxiliary variables x and  some of its powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x3 = x*x2
c...
c...  the function c.
c...
      numc = 7.389D-6*x2+2.3266D-2*x+2.568D-3
      dencm1 = 1/(7.389D-2*x3+4.72D-1*x2+8.723D0*x+1.0D0)
      C = dencm1*numc+1.667D-3
      cm1 = one/C
c...
c...  the function phi.
c...
      PHI = 8.13101627188389D-4*cm1*rm76*s12
      expphi = EXP(-PHI)
c...
c...  Non local contribution to the functional (dived by the total density)
c...
      enloc = C*expphi*rm73*s
c...
c...  First derivative of c
c...
      numcp = 1.4778D-5*x+2.3266D-2
      dencm2 = dencm1**2
      dencp = 2.2167D-1*x2+9.44D-1*x+8.723D0
      cp = dencm1*numcp-dencm2*dencp*numc
      cm2 = cm1**2
c...
c...  First derivative of phi
c...
      phipc = -8.13101627188389D-4*cm2*rm76*s12
      phipr = -9.486185650531205D-4*cm1*rm136*s12
      phips = 4.065508135941945D-4*cm1*rm76*sm12
c...
c...  First derivatives of the non local contribution
c...
      pnlocra = -1.433756689713499D-1*expphi*(2.324894703019253D0*C*(3.0
     $  D0*phipr*r+4.0D0)*r13-1.442249570307408D0*cp*(C*phipc-1.0D0))*rm
     $  83*s
      pnlocsaa = -C*expphi*rm43*(phips*s-1.0D0)
c...
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = enloc*scnl+eloc*scl
c...
c...  First derivatives of the functional (potential).
c...
      dra = pnlocra*scnl+plocra*scl
      dsaa = pnlocsaa*scnl
c...
      return
      end
c=======================================================================
      subroutine p86nlcdu(ra,rb,saa,sbb,sab,scl,scnl,ec,dra,drb,dsaa,dsb
     $  b,dsab,drara,drbrb,drarb,drasaa,drasbb,drasab,drbsaa,drbsbb,drbs
     $  ab,dsaasaa,dsbbsbb,dsabsab,dsaasbb,dsaasab,dsbbsab)
      implicit none
c...  MM (07/08/2003) generated with maxima, file p86nlcu.mc
c...
c...           Perdew 86 correlation functional
c...                   non local part
c...
c...     J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...                   open shell version
c...
c...     Formulation taken from the Molpro manual (P86)
c...     http://www.molpro.net/current/molpro_manual
c...  (note that the local part is different here than in Molpro)
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
      real*8 ra,rb,saa,sbb,sab,scl,scnl,ec,dra,drb,dsaa,dsbb,dsab
      real* 8 drara,drbrb,drarb,drasaa,drasbb,drasab
      real*8 drbsaa,drbsbb,drbsab,dsaasaa,dsbbsbb
      real*8 dsabsab,dsaasbb,dsaasab,dsbbsab
      real*8 eloc,plocra,plocrb,plocrara,plocrbrb,plocrarb
      real*8 enloc,pnlocra,pnlocrb,pnlocsaa,pnlocsbb,pnlocsab,pnlocrara,
     $  pnlocrbrb,pnlocrarb
      real*8 pnlocrasaa,pnlocrasbb,pnlocrasab,pnlocrbsaa,pnlocrbsbb,pnlo
     $  crbsab
      real*8 pnlocsaasaa,pnlocsaasbb,pnlocsaasab,pnlocsbbsbb,pnlocsbbsab
     $  ,pnlocsabsab
      real*8 r,r2,r3,r4,rm13,r13,r53,rm73,rm76,rm43,rm136,rm196,rm113,rm
     $  173
      real*8 s,s12,sm12,sm32
      real*8 z
      real*8 zp1,zp123,zp153,zp1m13,szpm53
      real*8 zm1,zm123,zm153,zm1m13
      real*8 x,x2,x3
      real*8 d,dp,dm1,dm2,d2,dm3,ds
      real*8 numc,numcp,dencm1,dencm2,dencm3,dencp,c,cm1,cm2,cm3,cp,cs
      real*8 phi,expphi,phipc,phipr,phips,phisc,phisr,phiss,phiscr,phisc
     $  s,phisrs
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,thirdm,tthird,sixth,etm,th,ssm
      parameter (third=1.0d0/3.0d0,thirdm=-1.0d0/3.0d0,tthird=2.0d0/3.0d
     $  0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0,ssm=-7.0d0/6.0d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      if(scl.ne.zero)then
      call pz81lcdu(ra,rb,eloc,plocra,plocrb,plocrara,plocrbrb,plocrarb)
      else
      eloc = ZERO
      plocra = ZERO
      plocrb = ZERO
      plocrara = ZERO
      plocrbrb = ZERO
      plocrarb = ZERO
      endif
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      r = rb+ra
      s = sbb+2.0D0*sab+saa
      if(s.lt.epsi.or.scnl.eq.zero)then
      enloc = ZERO
      pnlocra = ZERO
      pnlocrb = ZERO
      pnlocsaa = ZERO
      pnlocrara = ZERO
      pnlocrbrb = ZERO
      pnlocrarb = ZERO
      pnlocrasaa = ZERO
      pnlocrasbb = ZERO
      pnlocrasab = ZERO
      pnlocrbsaa = ZERO
      pnlocrbsbb = ZERO
      pnlocrbsab = ZERO
      pnlocsaasaa = ZERO
      pnlocsaasbb = ZERO
      pnlocsaasab = ZERO
      pnlocsbbsbb = ZERO
      pnlocsbbsab = ZERO
      pnlocsabsab = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
c...
      r2 = r**2
      r3 = r*r2
      r4 = r*r3
      r13 = r**THIRD
      r53 = r*r13**2
      rm13 = one/r13
      rm73 = rm13/r2
      rm76 = r**ssm
      rm43 = rm13/r
      rm113 = r**etm
      rm136 = rm76/r
      rm196 = rm136/r
      rm173 = rm113/r2
c...
      s12 = SQRT(s)
      sm12 = one/s12
      sm32 = sm12/s
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x3 = x*x2
      z = MIN(MAX((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp123 = zp1**tthird
      zp153 = zp1*zp123
      zp1m13 = zp1**thirdm
      zm1 = 1.0D0-z
      zm123 = zm1**tthird
      zm153 = zm1*zm123
      zm1m13 = zm1**thirdm
c...
c...  the function  d.
c...
      szpm53 = SQRT(zp153+zm153)
      d = 7.071067811865475D-1*szpm53
      dm1 = one/d
      dm2 = dm1**2
      d2 = d**2
      dm3 = dm2/d
c...
c...  the function c.
c...
      numc = 7.389D-6*x2+2.3266D-2*x+2.568D-3
      dencm1 = 1/(7.389D-2*x3+4.72D-1*x2+8.723D0*x+1.0D0)
      C = dencm1*numc+1.667D-3
      cm1 = one/C
c...
c...  the function phi.
c...
      PHI = 8.13101627188389D-4*cm1*rm76*s12
      expphi = EXP(-PHI)
c...
c...  Non local contribution to the functional (dived by the total density)
c...
      enloc = C*expphi*rm73*s/d
c...
c...  First derivative of d
c...
      if(abs(szpm53).gt.epsi)then
      dp = 5.892556509887896D-1*(zp123-zm123)/szpm53
      else
      dp = ZERO
      endif
c...
c...  First derivative of c
c...
      numcp = 1.4778D-5*x+2.3266D-2
      dencm2 = dencm1**2
      dencp = 2.2167D-1*x2+9.44D-1*x+8.723D0
      cp = dencm1*numcp-dencm2*dencp*numc
      cm2 = cm1**2
c...
c...  First derivative of phi
c...
      phipc = -8.13101627188389D-4*cm2*rm76*s12
      phipr = -9.486185650531205D-4*cm1*rm136*s12
      phips = 4.065508135941945D-4*cm1*rm76*sm12
c...
c...  First derivatives of the non local contribution
c...
      pnlocra = -1.433756689713499D-1*dm2*expphi*(2.324894703019253D0*C*
     $  r13*(6.0D0*dp*rb+3.0D0*d*phipr*r2+4.0D0*d*r)-1.442249570307408D0
     $  *cp*d*(C*phipc-1.0D0)*r)*rm113*s
      pnlocrb = 1.433756689713499D-1*dm2*expphi*(2.324894703019253D0*C*r
     $  13*(6.0D0*dp*ra-3.0D0*d*phipr*r2-4.0D0*d*r)+1.442249570307408D0*
     $  cp*d*(C*phipc-1.0D0)*r)*rm113*s
      pnlocsaa = -C*dm1*expphi*rm43*(phips*s-1.0D0)
c...
c...  second derivative of d
c...
      if(abs(szpm53).gt.epsi)then
      ds = 3.928371006591931D-1*(zp1m13+zm1m13)/szpm53-4.910463758239913
     $  D-1*(zp123-zm123)**2/szpm53**3
      else
      ds = ZERO
      endif
c...
c...  second derivative of c
c...
      dencm3 = dencm1*dencm2
      cs = -2.0D0*(numc*(dencm2*(2.2167D-1*x+4.72D-1)-dencm3*dencp
     $  **2)+dencm2*dencp*numcp-7.389D-6*dencm1)
      cm3 = cm2/C
c...
c...  Second derivative of phi
c...
      phisc = 1.626203254376778D-3*cm3*rm76*s12
      phisr = 2.055340224281761D-3*cm1*rm196*s12
      phiss = -2.032754067970973D-4*cm1*rm76*sm32
      phiscr = 9.486185650531205D-4*cm2*rm136*s12
      phiscs = -4.065508135941945D-4*cm2*rm76*sm12
      phisrs = -4.743092825265603D-4*cm1*rm136*sm12
c...
c...  Second derivatives of the non local contribution
c...
      pnlocrara = -2.055658245298212D-2*dm3*expphi*(2.145029397111026D0*
     $  C*r13*(9.071431559243087D1*(d*ds-2.0D0*dp**2)*rb**2-3.0238105197
     $  47696D1*d*dp*r*(3.0D0*phipr*r+7.0D0)*rb+2.519842099789746D0*d2*r
     $  2*(9.0D0*phisr*r2-9.0D0*phipr**2*r2-2.4D1*phipr*r-2.8D1))+C*(1.4
     $  64591887561523D0*cp*(2.747314182127997D1*d*dp*phipc*r*rb-1.37365
     $  7091063998D1*d2*(phiscr*r-phipc*phipr*r-2.0D0*phipc)*r2)+2
     $  .080083823051904D0*cp**2*d2*(phisc-phipc**2)*r53+2.0800838
     $  23051904D0*cs*d2*phipc*r53)+1.464591887561523D0*cp*(-2.747314182
     $  127997D1*d*dp*r*rb-1.373657091063998D1*d2*(phipr*r+2.0D0)*r2)+4.
     $  160167646103808D0*cp**2*d2*phipc*r53-2.080083823051904D0*cs*d2*r
     $  53)*rm173*s
      pnlocrbrb = -2.055658245298212D-2*dm3*expphi*(2.145029397111026D0*
     $  C*r13*(9.071431559243087D1*(d*ds-2.0D0*dp**2)*ra**2+3.0238105197
     $  47696D1*d*dp*r*(3.0D0*phipr*r+7.0D0)*ra+2.519842099789746D0*d2*r
     $  2*(9.0D0*phisr*r2-9.0D0*phipr**2*r2-2.4D1*phipr*r-2.8D1))+C*(1.4
     $  64591887561523D0*cp*(-2.747314182127997D1*d*dp*phipc*r*ra-1.3736
     $  57091063998D1*d2*(phiscr*r-phipc*phipr*r-2.0D0*phipc)*r2)+
     $  2.080083823051904D0*cp**2*d2*(phisc-phipc**2)*r53+2.080083
     $  823051904D0*cs*d2*phipc*r53)+1.464591887561523D0*cp*(2.747314182
     $  127997D1*d*dp*r*ra-1.373657091063998D1*d2*(phipr*r+2.0D0)*r2)+4.
     $  160167646103808D0*cp**2*d2*phipc*r53-2.080083823051904D0*cs*d2*r
     $  53)*rm173*s
      pnlocrarb = 2.055658245298212D-2*dm3*expphi*(C*(-2.011847031863692
     $  D1*cp*d*r*(dp*phipc*rb-dp*phipc*ra-d*phiscr*r2+d*phi
     $  pc*phipr*r2+2.0D0*d*phipc*r)-2.080083823051904D0*cp**2*d2*(phisc
     $  -phipc**2)*r53-2.080083823051904D0*cs*d2*phipc*r53)+5.4051
     $  3538012698D0*C*r13*(3.6D1*d*ds*ra*rb-7.2D1*dp**2*ra*rb+1.8D1*d*d
     $  p*phipr*r2*rb+4.2D1*d*dp*r*rb-1.8D1*d*dp*phipr*r2*ra-4.2D1*d*dp*
     $  r*ra-9.0D0*d2*phisr*r4+9.0D0*d2*phipr**2*r4+2.4D1*d2*phipr*r3+2.
     $  8D1*d2*r2)+2.011847031863692D1*cp*d*r*(dp*rb-dp*ra+d*phipr
     $  *r2+2.0D0*d*r)-4.160167646103808D0*cp**2*d2*phipc*r53+2.08008382
     $  3051904D0*cs*d2*r53)*rm173*s
      pnlocrasaa = 1.433756689713499D-1*dm2*expphi*rm113*(1.464591887561
     $  523D0*r13*(-4.762203155904598D0*C*d*r2*(phisrs*s-phipr*phi
     $  ps*s+phipr)+9.524406311809197D0*C*dp*rb*(phips*s-1.0D0)+6.349604
     $  207872798D0*C*d*r*(phips*s-1.0D0))+cp*d*r*(1.442249570307408D0*C
     $  *(phiscs*s-phipc*phips*s+phipc)+1.442249570307408D0*(phips
     $  *s-1.0D0)))
      pnlocrbsaa = -1.433756689713499D-1*dm2*expphi*rm113*(1.46459188756
     $  1523D0*r13*(4.762203155904598D0*C*d*r2*(phisrs*s-phipr*phi
     $  ps*s+phipr)+9.524406311809197D0*C*dp*ra*(phips*s-1.0D0)-6.349604
     $  207872798D0*C*d*r*(phips*s-1.0D0))+cp*d*r*(-1.442249570307408D0*
     $  C*(phiscs*s-phipc*phips*s+phipc)-1.442249570307408D0*(phip
     $  s*s-1.0D0)))
      pnlocsaasaa = -C*dm1*expphi*rm43*((phiss-phips**2)*s+2
     $  .0D0*phips)
c...
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = enloc*scnl+eloc*scl
c...
c...  First derivatives of the functional (potential).
c...
      dra = pnlocra*scnl+plocra*scl
      drb = pnlocrb*scnl+plocrb*scl
      dsaa = pnlocsaa*scnl
      dsbb = dsaa
      dsab = 2.0D0*dsaa
c...
c...  Second derivatives of the functional.
c...
      drara = pnlocrara*scnl+plocrara*scl
      drbrb = pnlocrbrb*scnl+plocrbrb*scl
      drarb = pnlocrarb*scnl+plocrarb*scl
      drasaa = pnlocrasaa*scnl
      drasbb = drasaa
      drasab = 2.0D0*drasaa
      drbsaa = pnlocrbsaa*scnl
      drbsbb = drbsaa
      drbsab = 2.0D0*drbsaa
      dsaasaa = pnlocsaasaa*scnl
      dsaasbb = dsaasaa
      dsaasab = 2.0D0*dsaasaa
      dsbbsbb = dsaasaa
      dsbbsab = dsaasab
      dsabsab = 2.0D0*dsaasab
c...
      return
      end
c=======================================================================
      subroutine p86nlcu(ra,rb,saa,sbb,sab,scl,scnl,ec,dra,drb,dsaa,dsbb
     $  ,dsab)
      implicit none
c...  MM (07/08/2003) generated with maxima, file p86nlcu.mc
c...
c...           Perdew 86 correlation functional
c...                   non local part
c...
c...     J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...                   open shell version
c...
c...     Formulation taken from the Molpro manual (P86)
c...     http://www.molpro.net/current/molpro_manual
c...  (note that the local part is different here than in Molpro)
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
      real*8 eloc,plocra,plocrb
      real*8 enloc,pnlocra,pnlocrb,pnlocsaa,pnlocsbb,pnlocsab
      real*8 r,r2,rm13,r13,rm73,rm76,rm43,rm136,rm113
      real*8 s,s12,sm12
      real*8 z
      real*8 zp1,zp123,zp153,szpm53
      real*8 zm1,zm123,zm153
      real*8 x,x2,x3
      real*8 d,dp,dm1,dm2
      real*8 numc,numcp,dencm1,dencm2,dencp,c,cm1,cm2,cp
      real*8 phi,expphi,phipc,phipr,phips
c...
      real*8 zero,one,two
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      real*8 third,thirdm,tthird,sixth,etm,th,ssm
      parameter (third=1.0d0/3.0d0,thirdm=-1.0d0/3.0d0,tthird=2.0d0/3.0d
     $  0,sixth=1.0d0/6.0d0)
      parameter (etm=-11.0d0/3.0d0,th=1.5d0,ssm=-7.0d0/6.0d0)
c...
c...  threshold for gradient
c...
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  first, lets get the local part out of the way
c...
      if(scl.ne.zero)then
      call pz81lcu(ra,rb,eloc,plocra,plocrb)
      else
      eloc = ZERO
      plocra = ZERO
      plocrb = ZERO
      endif
c...
c...  Now, if the total gradient is bigger than the threshold,
c...  we compute the non local part, otherwise we are basically done
c...
      r = rb+ra
      s = sbb+2.0D0*sab+saa
      if(s.lt.epsi.or.scnl.eq.zero)then
      enloc = ZERO
      pnlocra = ZERO
      pnlocrb = ZERO
      pnlocsaa = ZERO
      else
c...
c...  Non local part.
c...  we compute the total density r and some powers of s
c...
      r2 = r**2
      r13 = r**THIRD
      rm13 = one/r13
      rm73 = rm13/r2
      rm76 = r**ssm
      rm43 = rm13/r
      rm113 = r**etm
      rm136 = rm76/r
c...
      s12 = SQRT(s)
      sm12 = one/s12
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      x2 = x**2
      x3 = x*x2
      z = MIN(MAX((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp123 = zp1**tthird
      zp153 = zp1*zp123
      zm1 = 1.0D0-z
      zm123 = zm1**tthird
      zm153 = zm1*zm123
c...
c...  the function  d.
c...
      szpm53 = SQRT(zp153+zm153)
      d = 7.071067811865475D-1*szpm53
      dm1 = one/d
      dm2 = dm1**2
c...
c...  the function c.
c...
      numc = 7.389D-6*x2+2.3266D-2*x+2.568D-3
      dencm1 = 1/(7.389D-2*x3+4.72D-1*x2+8.723D0*x+1.0D0)
      C = dencm1*numc+1.667D-3
      cm1 = one/C
c...
c...  the function phi.
c...
      PHI = 8.13101627188389D-4*cm1*rm76*s12
      expphi = EXP(-PHI)
c...
c...  Non local contribution to the functional (dived by the total density)
c...
      enloc = C*expphi*rm73*s/d
c...
c...  First derivative of d
c...
      if(abs(szpm53).gt.epsi)then
      dp = 5.892556509887896D-1*(zp123-zm123)/szpm53
      else
      dp = ZERO
      endif
c...
c...  First derivative of c
c...
      numcp = 1.4778D-5*x+2.3266D-2
      dencm2 = dencm1**2
      dencp = 2.2167D-1*x2+9.44D-1*x+8.723D0
      cp = dencm1*numcp-dencm2*dencp*numc
      cm2 = cm1**2
c...
c...  First derivative of phi
c...
      phipc = -8.13101627188389D-4*cm2*rm76*s12
      phipr = -9.486185650531205D-4*cm1*rm136*s12
      phips = 4.065508135941945D-4*cm1*rm76*sm12
c...
c...  First derivatives of the non local contribution
c...
      pnlocra = -1.433756689713499D-1*dm2*expphi*(2.324894703019253D0*C*
     $  r13*(6.0D0*dp*rb+3.0D0*d*phipr*r2+4.0D0*d*r)-1.442249570307408D0
     $  *cp*d*(C*phipc-1.0D0)*r)*rm113*s
      pnlocrb = 1.433756689713499D-1*dm2*expphi*(2.324894703019253D0*C*r
     $  13*(6.0D0*dp*ra-3.0D0*d*phipr*r2-4.0D0*d*r)+1.442249570307408D0*
     $  cp*d*(C*phipc-1.0D0)*r)*rm113*s
      pnlocsaa = -C*dm1*expphi*rm43*(phips*s-1.0D0)
c...
      endif
c...
c...  We assemble the output values.
c...
c...  Functional value (divided by the density).
c...
      ec = enloc*scnl+eloc*scl
c...
c...  First derivatives of the functional (potential).
c...
      dra = pnlocra*scnl+plocra*scl
      drb = pnlocrb*scnl+plocrb*scl
      dsaa = pnlocsaa*scnl
      dsbb = dsaa
      dsab = 2.0D0*dsaa
c...
      return
      end
c=======================================================================
      subroutine pz81lcr(ra,ec,dr)
      implicit none
c...
c...  MM (07/08/2003) generated with maxima, file p86lcr.mc
c...
c...                Perdew 86 correlation functional
c...                local part (aka Perdew-Zunger 81)
c...
c...          J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...     J.P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
c...                        closed shell version
c...
c...         Formulation taken from the Density Functional Depository
c...     http://www.dl.ac.uk/TCSC/QuantumChem/dft_lib/lib/def_c_p86.html
c...  (for the actual definition it is better to look directly at the code)
c...
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra    half the closed-shell density
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dr     functional derivative with respect to density
c...
      real*8 ra,ec,dr
      real*8 z
      real*8 x,xm12,x12,xm1,xln
      real*8 r,rm13
      real*8 denepsp,denepsp2,denepspp,epsp,epspp
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 thirdm,thm
      parameter (thirdm=-0.33333333333333333333d0,thm=-1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      r = 2.0D0*ra
      rm13 = r**thirdm
c...
c...  the auxiliary variable x and some of its powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
c...
c...  epsp.
c...  There are two definitions, one for x < 1, the other elsewhere.
c...
      if(x.lt.one)then
      xln = LOG(x)
      epsp = 2.0D-3*x*xln+3.11D-2*xln-1.16D-2*x-4.8D-2
      else
      x12 = SQRT(x)
      denepsp = 1/(1.0529D0*x12+3.334D-1*x+1.0D0)
      epsp = -1.423D-1*denepsp
      endif
c...
c...  Functional value (divided by the density).
c...
      ec = epsp
c...
c...  First derivative (potential).
c...
      if(x.lt.one)then
      xm1 = one/x
      epspp = 3.11D-2*xm1+2.0D-3*xln-9.6D-3
      else
      xm12 = one/x12
      denepsp2 = denepsp**2
      denepspp = 5.2645D-1*xm12+3.334D-1
      epspp = 1.423D-1*denepsp2*denepspp
      endif
      dr = epsp-2.067834969664667D-1*epspp*rm13
c...
      return
      end
c=======================================================================
      subroutine pz81lcdr(ra,ec,dr,drr)
      implicit none
c...
c...  MM (07/08/2003) generated with maxima, file p86lcr.mc
c...
c...                Perdew 86 correlation functional
c...                local part (aka Perdew-Zunger 81)
c...
c...          J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...     J.P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
c...                        closed shell version
c...
c...         Formulation taken from the Density Functional Depository
c...     http://www.dl.ac.uk/TCSC/QuantumChem/dft_lib/lib/def_c_p86.html
c...  (for the actual definition it is better to look directly at the code)
c...
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
c...
c...  Input parameters:
c...
c...  ra    half the closed-shell density
c...
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dr     functional derivative with respect to density
c...  drr   funct. 2nd deriv. w.r.t. density
c...
      real*8 ra,ec,dr
      real*8 drr
      real*8 z
      real*8 x,xm12,x12,xm32,xm1,xm2,xln
      real*8 r,rm13,rm43,rm53
      real*8 denepsp,denepsp2,denepsp3,denepspp,epsp,epspp,epsps
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 thirdm,thm
      parameter (thirdm=-0.33333333333333333333d0,thm=-1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      r = 2.0D0*ra
      rm13 = r**thirdm
      rm43 = rm13/r
      rm53 = rm13*rm43
c...
c...  the auxiliary variable x and some of its powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
c...
c...  epsp.
c...  There are two definitions, one for x < 1, the other elsewhere.
c...
      if(x.lt.one)then
      xln = LOG(x)
      epsp = 2.0D-3*x*xln+3.11D-2*xln-1.16D-2*x-4.8D-2
      else
      x12 = SQRT(x)
      denepsp = 1/(1.0529D0*x12+3.334D-1*x+1.0D0)
      epsp = -1.423D-1*denepsp
      endif
c...
c...  Functional value (divided by the density).
c...
      ec = epsp
c...
c...  First derivative (potential).
c...
      if(x.lt.one)then
      xm1 = one/x
      epspp = 3.11D-2*xm1+2.0D-3*xln-9.6D-3
      else
      xm12 = one/x12
      denepsp2 = denepsp**2
      denepspp = 5.2645D-1*xm12+3.334D-1
      epspp = 1.423D-1*denepsp2*denepspp
      endif
      dr = epsp-2.067834969664667D-1*epspp*rm13
c...
c...  Second derivatives.
c...
      if(x.lt.one)then
      xm2 = xm1**2
      epsps = 2.0D-3*xm1-3.11D-2*xm2
      else
      xm32 = x**thm
      denepsp3 = denepsp*denepsp2
      epsps = -3.74569175D-2*denepsp2*xm32-2.846D-1*denepsp3*denepspp**2
      endif
      drr = 4.275941461768073D-2*epsps*rm53-1.378556646443111D-1*epspp*r
     $  m43
c...
      return
      end
c=======================================================================
      subroutine pz81lcu(ra,rb,ec,dra,drb)
      implicit none
c...
c...  MM (07/02/2003) generated with maxima, file p86lcu.mc
c...
c...                Perdew 86 correlation functional
c...                local part (aka Perdew-Zunger 81)
c...
c...          J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...     J.P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
c...                        open shell version
c...
c...         Formulation taken from the Density Functional Depository
c...     http://www.dl.ac.uk/TCSC/QuantumChem/dft_lib/lib/def_c_p86.html
c...  (for the actual definition it is better to look directly at the code)
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
      real*8 z
      real*8 zp1,zp113,zp143
      real*8 zm1,zm113,zm143
      real*8 x,xm12,x12,xm1,xln
      real*8 r,rm13,r13,rm43,r43
      real*8 denepsp,denepsp2,denepspp,epsp,denepsf,denepsf2,denepsfp,ep
     $  sf,epspp,epsfp
      real*8 om,omp
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
      real*8 third,tthird,etm,thm
      parameter (third=0.33333333333333333333d0)
      parameter (tthird=0.66666666666666666667d0)
      parameter (etm=-11.0d0/3.0d0,thm=-1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      r = rb+ra
      r13 = r**THIRD
      rm13 = one/r13
      r43 = r*r13
      rm43 = one/r43
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      z = min(max((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp143 = zp1*zp113
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm143 = zm1*zm113
c...
c...  epsp and epsf.
c...  There are two definitions, one for x < 1, the other elsewhere.
c...
      if(x.lt.one)then
      xln = LOG(x)
      epsp = 2.0D-3*x*xln+3.11D-2*xln-1.16D-2*x-4.8D-2
      epsf = 7.0D-4*x*xln+1.555D-2*xln-4.8D-3*x-2.69D-2
      else
      x12 = SQRT(x)
      denepsp = 1/(1.0529D0*x12+3.334D-1*x+1.0D0)
      epsp = -1.423D-1*denepsp
      denepsf = 1/(1.3981D0*x12+2.611D-1*x+1.0D0)
      epsf = -8.43D-2*denepsf
      endif
c...
c...  the function omega.
c...  the constant is 1/(2^(4/3)-2)
c...
      om = 1.923661050931536D0*(zp143+zm143-2.0D0)
c...
c...  Functional value (divided by the density).
c...
      ec = (epsf-epsp)*om+epsp
c...
c...  First derivative (potential).
c...
      if(x.lt.one)then
      xm1 = one/x
      epspp = 3.11D-2*xm1+2.0D-3*xln-9.6D-3
      epsfp = 1.555D-2*xm1+7.0D-4*xln-4.1D-3
      else
      xm12 = one/x12
      denepsp2 = denepsp**2
      denepspp = 5.2645D-1*xm12+3.334D-1
      epspp = 1.423D-1*denepsp2*denepspp
      denepsf2 = denepsf**2
      denepsfp = 6.9905D-1*xm12+2.611D-1
      epsfp = 8.43D-2*denepsf2*denepsfp
      endif
      omp = 2.564881401242048D0*(zp113-zm113)
      dra = -1.433756689713499D-1*(6.974684109057759D0*epsp*r13*(2.0D0*o
     $  mp*rb+om*r-r)-6.974684109057759D0*epsf*r13*(2.0D0*omp*rb+o
     $  m*r)-1.442249570307408D0*epspp*(om-1.0D0)*r+1.442249570307408D0*
     $  epsfp*om*r)*rm43
      drb = 1.433756689713499D-1*(6.974684109057759D0*epsp*r13*(2.0D0*om
     $  p*ra-1.0D0*om*r+r)-6.974684109057759D0*epsf*r13*(2.0D0*omp*ra-
     $  om*r)+1.442249570307408D0*epspp*(om-1.0D0)*r-1.4422495703074
     $  08D0*epsfp*om*r)*rm43
c...
      return
      end
c=======================================================================
      subroutine pz81lcdu(ra,rb,ec,dra,drb,drara,drbrb,drarb)
      implicit none
c...
c...  MM (07/02/2003) generated with maxima, file p86lcu.mc
c...
c...                Perdew 86 correlation functional
c...                local part (aka Perdew-Zunger 81)
c...
c...          J.P. Perdew, Phys. Rev. B 33, 8822 (1986)
c...     J.P. Perdew and A. Zunger, Phys. Rev. B 23, 5048 (1981)
c...                        open shell version
c...
c...         Formulation taken from the Density Functional Depository
c...     http://www.dl.ac.uk/TCSC/QuantumChem/dft_lib/lib/def_c_p86.html
c...  (for the actual definition it is better to look directly at the code)
c...
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...
      real*8 ra,rb,ec,dra,drb
      real*8 drara,drbrb,drarb
      real*8 z
      real*8 zp1,zp113,zp123,zp143
      real*8 zm1,zm113,zm123,zm143
      real*8 x,xm12,x12,xm32,xm1,xm2,xln
      real*8 r,ra2,rb2,r2,rm13,r13,r23,rm43,r43,rm113
      real*8 denepsp,denepsp2,denepsp3,denepspp,epsp,denepsf,denepsf2,de
     $  nepsf3,denepsfp,epsf,epspp,epsfp,epsps,epsfs
      real*8 om,omp,oms
c...
      real*8 zero,one
      parameter (zero=0.0d0,one=1.0d0)
      real*8 epsi,epsim1
      parameter (epsi=1.0d-15,epsim1=1.0d+15)
      real*8 third,tthird,etm,thm
      parameter (third=0.33333333333333333333d0)
      parameter (tthird=0.66666666666666666667d0)
      parameter (etm=-11.0d0/3.0d0,thm=-1.5d0)
c...
c...  first we compute some powers of r that will be used later
c...
      r = rb+ra
      ra2 = ra**2
      rb2 = rb**2
      r2 = r**2
      r13 = r**THIRD
      rm13 = one/r13
      r23 = r13**2
      r43 = r*r13
      rm43 = one/r43
      rm113 = r**etm
c...
c...  the auxiliary variables x and z and some of their powers.
c...  the constant is (3/(4*%pi))^(1/3)
c...
      x = 6.203504908994D-1*rm13
      z = min(max((ra-rb)/r,-one),one)
      zp1 = z+1.0D0
      zp113 = zp1**THIRD
      zp123 = zp113**2
      zp143 = zp1*zp113
      zm1 = 1.0D0-z
      zm113 = zm1**THIRD
      zm123 = zm113**2
      zm143 = zm1*zm113
c...
c...  epsp and epsf.
c...  There are two definitions, one for x < 1, the other elsewhere.
c...
      if(x.lt.one)then
      xln = LOG(x)
      epsp = 2.0D-3*x*xln+3.11D-2*xln-1.16D-2*x-4.8D-2
      epsf = 7.0D-4*x*xln+1.555D-2*xln-4.8D-3*x-2.69D-2
      else
      x12 = SQRT(x)
      denepsp = 1/(1.0529D0*x12+3.334D-1*x+1.0D0)
      epsp = -1.423D-1*denepsp
      denepsf = 1/(1.3981D0*x12+2.611D-1*x+1.0D0)
      epsf = -8.43D-2*denepsf
      endif
c...
c...  the function omega.
c...  the constant is 1/(2^(4/3)-2)
c...
      om = 1.923661050931536D0*(zp143+zm143-2.0D0)
c...
c...  Functional value (divided by the density).
c...
      ec = (epsf-epsp)*om+epsp
c...
c...  First derivative (potential).
c...
      if(x.lt.one)then
      xm1 = one/x
      epspp = 3.11D-2*xm1+2.0D-3*xln-9.6D-3
      epsfp = 1.555D-2*xm1+7.0D-4*xln-4.1D-3
      else
      xm12 = one/x12
      denepsp2 = denepsp**2
      denepspp = 5.2645D-1*xm12+3.334D-1
      epspp = 1.423D-1*denepsp2*denepspp
      denepsf2 = denepsf**2
      denepsfp = 6.9905D-1*xm12+2.611D-1
      epsfp = 8.43D-2*denepsf2*denepsfp
      endif
      omp = 2.564881401242048D0*(zp113-zm113)
      dra = -1.433756689713499D-1*(6.974684109057759D0*epsp*r13*(2.0D0*o
     $  mp*rb+om*r-r)-6.974684109057759D0*epsf*r13*(2.0D0*omp*rb+o
     $  m*r)-1.442249570307408D0*epspp*(om-1.0D0)*r+1.442249570307408D0*
     $  epsfp*om*r)*rm43
      drb = 1.433756689713499D-1*(6.974684109057759D0*epsp*r13*(2.0D0*om
     $  p*ra-om*r+r)-6.974684109057759D0*epsf*r13*(2.0D0*omp*ra-
     $  om*r)+1.442249570307408D0*epspp*(om-1.0D0)*r-1.4422495703074
     $  08D0*epsfp*om*r)*rm43
c...
c...  Second derivatives.
c...
      if(x.lt.one)then
      xm2 = xm1**2
      epsps = 2.0D-3*xm1-3.11D-2*xm2
      epsfs = 7.0D-4*xm1-1.555D-2*xm2
      else
      xm32 = x**thm
      denepsp3 = denepsp*denepsp2
      epsps = -3.74569175D-2*denepsp2*xm32-2.846D-1*denepsp3*denepspp**2
      denepsf3 = denepsf*denepsf2
      epsfs = -2.94649575D-2*denepsf2*xm32-1.686D-1*denepsf3*denepsfp**2
      endif
      if(zm123*zp123.gt.epsi)then
      oms = 8.549604670806828D-1*(zp123+zm123)/(zm123*zp123)
      else
      oms = ZERO
      endif
      drara = -2.055658245298212D-2*(1.945848736845713D2*epsp*oms*r23*rb
     $  2-1.945848736845713D2*epsf*oms*r23*rb2-6.706156772878975D0*epspp
     $  *r43*(6.0D0*omp*rb+om*r-r)+6.706156772878975D0*epsfp*r43*(
     $  6.0D0*omp*rb+om*r)+2.080083823051904D0*epsps*(om-1.0D0)*r2-2.080
     $  083823051904D0*epsfs*om*r2)*rm113
      drbrb = -2.055658245298212D-2*(1.945848736845713D2*epsp*oms*r23*ra
     $  2-1.945848736845713D2*epsf*oms*r23*ra2+6.706156772878975D0*epspp
     $  *r43*(6.0D0*omp*ra-om*r+r)-6.706156772878975D0*epsfp*r43*(
     $  6.0D0*omp*ra-om*r)+2.080083823051904D0*epsps*(om-1.0D0)*r2
     $  -2.080083823051904D0*epsfs*om*r2)*rm113
      drarb = 2.055658245298212D-2*(6.706156772878975D0*epspp*r43*(3.0D0
     $  *omp*rb-3.0D0*omp*ra+om*r-r)-6.706156772878975D0*epsfp*r43
     $  *(3.0D0*omp*rb-3.0D0*omp*ra+om*r)+1.945848736845713D2*epsp*oms*r
     $  23*ra*rb-1.945848736845713D2*epsf*oms*r23*ra*rb-2.08008382305190
     $  4D0*epsps*(om-1.0D0)*r2+2.080083823051904D0*epsfs*om*r2)*rm113
c...
      return
      end

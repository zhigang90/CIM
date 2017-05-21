c=======================================================================
      subroutine hcth407r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 407
c...
c...               A. D. Doese and  N. C. Handy
c...              J. Chem. Phys. 114, 5497 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.2005D1*ETAcvl1chi2a4+4.2572D1*ETAcvl1chi2a3-1.9222D1
     1   *ETAcvl1chi2a2+4.4237D0*ETAcvl1chi2a1+5.8908D-1
      sumachiapchi2a = -4.2005D1*ETAcvl1chi2a4pchi2a+4.2572D1*ETAcvl1chi
     1   2a3pchi2a-1.9222D1*ETAcvl1chi2a2pchi2a+4.4237D0*ETAcvl1chi2a1pc
     2   hi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 6.248D0*ETAcvl2chi2a4-9.1792D0*ETAcvl2chi2a3+5.6174D0*E
     1   TAcvl2chi2a2-2.4029D0*ETAcvl2chi2a1+1.18777D0
      sumbchiapchi2a = 6.248D0*ETAcvl2chi2a4pchi2a-9.1792D0*ETAcvl2chi2a
     1   3pchi2a+5.6174D0*ETAcvl2chi2a2pchi2a-2.4029D0*ETAcvl2chi2a1pchi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 2.2886D0*ETAcvl3chi2a4-2.629D0*ETAcvl3chi2a3+3.4256D0*E
     1   TAcvl3chi2a2-5.183D-1*ETAcvl3chi2a1+1.08184D0
      sumcchiapchi2a = 2.2886D0*ETAcvl3chi2a4pchi2a-2.629D0*ETAcvl3chi2a
     1   3pchi2a+3.4256D0*ETAcvl3chi2a2pchi2a-5.183D-1*ETAcvl3chi2a1pchi
     2   2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth407u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 407
c...
c...               A. D. Doese and  N. C. Handy
c...              J. Chem. Phys. 114, 5497 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 6.248D0*ETAcvl2chi2a4-9.1792D0*ETAcvl2chi2a3+5.6174D0*E
     1   TAcvl2chi2a2-2.4029D0*ETAcvl2chi2a1+1.18777D0
      sumbchiapchi2a = 6.248D0*ETAcvl2chi2a4pchi2a-9.1792D0*ETAcvl2chi2a
     1   3pchi2a+5.6174D0*ETAcvl2chi2a2pchi2a-2.4029D0*ETAcvl2chi2a1pchi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 2.2886D0*ETAcvl3chi2a4-2.629D0*ETAcvl3chi2a3+3.4256D0*E
     1   TAcvl3chi2a2-5.183D-1*ETAcvl3chi2a1+1.08184D0
      sumcchiapchi2a = 2.2886D0*ETAcvl3chi2a4pchi2a-2.629D0*ETAcvl3chi2a
     1   3pchi2a+3.4256D0*ETAcvl3chi2a2pchi2a-5.183D-1*ETAcvl3chi2a1pchi
     2   2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      chi2b3 = chi2b*chi2b2
      chi2b4 = chi2b*chi2b3
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 6.248D0*ETAcvl2chi2b4-9.1792D0*ETAcvl2chi2b3+5.6174D0*E
     1   TAcvl2chi2b2-2.4029D0*ETAcvl2chi2b1+1.18777D0
      sumbchibpchi2b = 6.248D0*ETAcvl2chi2b4pchi2b-9.1792D0*ETAcvl2chi2b
     1   3pchi2b+5.6174D0*ETAcvl2chi2b2pchi2b-2.4029D0*ETAcvl2chi2b1pchi
     2   2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 2.2886D0*ETAcvl3chi2b4-2.629D0*ETAcvl3chi2b3+3.4256D0*E
     1   TAcvl3chi2b2-5.183D-1*ETAcvl3chi2b1+1.08184D0
      sumcchibpchi2b = 2.2886D0*ETAcvl3chi2b4pchi2b-2.629D0*ETAcvl3chi2b
     1   3pchi2b+3.4256D0*ETAcvl3chi2b2pchi2b-5.183D-1*ETAcvl3chi2b1pchi
     2   2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      d23 = d2*d22
      d24 = d2*d23
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
c...
c...  sumad and its derivatives
c...
      sumad = -4.2005D1*ETAcvl1d24+4.2572D1*ETAcvl1d23-1.9222D1*ETAcvl1d
     1   22+4.4237D0*ETAcvl1d21+5.8908D-1
      sumadpd2 = -4.2005D1*ETAcvl1d24pd2+4.2572D1*ETAcvl1d23pd2-1.9222D1
     1   *ETAcvl1d22pd2+4.4237D0*ETAcvl1d21pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine b97r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...           Becke 97 (B97) functional
c...
c...  A. D. Becke, J. Chem. Phys. 107, 8554 (1997)
c...
c...  Formulation taken from the Molpro manual (B97)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.5961D0*ETAcvl1chi2a2+7.471D-1*ETAcvl1chi2a1+9.454D-1
      sumachiapchi2a = 7.471D-1*ETAcvl1chi2a1pchi2a-4.5961D0*ETAcvl1chi2
     1   a2pchi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.4868D0*ETAcvl2chi2a2+2.3487D0*ETAcvl2chi2a1+1.737D-1
      sumbchiapchi2a = 2.3487D0*ETAcvl2chi2a1pchi2a-2.4868D0*ETAcvl2chi2
     1   a2pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 7.481D-1*ETAcvl3chi2a2+5.073D-1*ETAcvl3chi2a1+8.094D-1
      sumcchiapchi2a = 7.481D-1*ETAcvl3chi2a2pchi2a+5.073D-1*ETAcvl3chi2
     1   a1pchi2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine b97u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...           Becke 97 (B97) functional
c...
c...  A. D. Becke, J. Chem. Phys. 107, 8554 (1997)
c...
c...  Formulation taken from the Molpro manual (B97)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.4868D0*ETAcvl2chi2a2+2.3487D0*ETAcvl2chi2a1+1.737D-1
      sumbchiapchi2a = 2.3487D0*ETAcvl2chi2a1pchi2a-2.4868D0*ETAcvl2chi2
     1   a2pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 7.481D-1*ETAcvl3chi2a2+5.073D-1*ETAcvl3chi2a1+8.094D-1
      sumcchiapchi2a = 7.481D-1*ETAcvl3chi2a2pchi2a+5.073D-1*ETAcvl3chi2
     1   a1pchi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
c...
c...  sumbchib and its derivatives
c...
      sumbchib = -2.4868D0*ETAcvl2chi2b2+2.3487D0*ETAcvl2chi2b1+1.737D-1
      sumbchibpchi2b = 2.3487D0*ETAcvl2chi2b1pchi2b-2.4868D0*ETAcvl2chi2
     1   b2pchi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 7.481D-1*ETAcvl3chi2b2+5.073D-1*ETAcvl3chi2b1+8.094D-1
      sumcchibpchi2b = 7.481D-1*ETAcvl3chi2b2pchi2b+5.073D-1*ETAcvl3chi2
     1   b1pchi2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
c...
c...  sumad and its derivatives
c...
      sumad = -4.5961D0*ETAcvl1d22+7.471D-1*ETAcvl1d21+9.454D-1
      sumadpd2 = 7.471D-1*ETAcvl1d21pd2-4.5961D0*ETAcvl1d22pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine b971r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...             Becke 97-1  functional
c...  This is the Becke-97 functional reparametrized
c...
c...  F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N. C. Handy,
c...              J. Chem. Phys. 109, 6264 (1998)
c...
c...  Formulation taken from the Molpro manual (B97R)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumachia and its derivatives
c...
      sumachia = -5.47869D0*ETAcvl1chi2a2+7.88552D-1*ETAcvl1chi2a1+9.556
     1   89D-1
      sumachiapchi2a = 7.88552D-1*ETAcvl1chi2a1pchi2a-5.47869D0*ETAcvl1c
     1   hi2a2pchi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.87103D0*ETAcvl2chi2a2+2.71681D0*ETAcvl2chi2a1+8.2001
     1   1D-2
      sumbchiapchi2a = 2.71681D0*ETAcvl2chi2a1pchi2a-2.87103D0*ETAcvl2ch
     1   i2a2pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 6.60975D-1*ETAcvl3chi2a2+5.73805D-1*ETAcvl3chi2a1+7.895
     1   18D-1
      sumcchiapchi2a = 6.60975D-1*ETAcvl3chi2a2pchi2a+5.73805D-1*ETAcvl3
     1   chi2a1pchi2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra1
     5   13
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine b971u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...             Becke 97-1  functional
c...  This is the Becke-97 functional reparametrized
c...
c...  F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N. C. Handy,
c...              J. Chem. Phys. 109, 6264 (1998)
c...
c...  Formulation taken from the Molpro manual (D97R)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.87103D0*ETAcvl2chi2a2+2.71681D0*ETAcvl2chi2a1+8.2001
     1   1D-2
      sumbchiapchi2a = 2.71681D0*ETAcvl2chi2a1pchi2a-2.87103D0*ETAcvl2ch
     1   i2a2pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 6.60975D-1*ETAcvl3chi2a2+5.73805D-1*ETAcvl3chi2a1+7.895
     1   18D-1
      sumcchiapchi2a = 6.60975D-1*ETAcvl3chi2a2pchi2a+5.73805D-1*ETAcvl3
     1   chi2a1pchi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
c...
c...  sumbchib and its derivatives
c...
      sumbchib = -2.87103D0*ETAcvl2chi2b2+2.71681D0*ETAcvl2chi2b1+8.2001
     1   1D-2
      sumbchibpchi2b = 2.71681D0*ETAcvl2chi2b1pchi2b-2.87103D0*ETAcvl2ch
     1   i2b2pchi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 6.60975D-1*ETAcvl3chi2b2+5.73805D-1*ETAcvl3chi2b1+7.895
     1   18D-1
      sumcchibpchi2b = 6.60975D-1*ETAcvl3chi2b2pchi2b+5.73805D-1*ETAcvl3
     1   chi2b1pchi2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
c...
c...  sumad and its derivatives
c...
      sumad = -5.47869D0*ETAcvl1d22+7.88552D-1*ETAcvl1d21+9.55689D-1
      sumadpd2 = 7.88552D-1*ETAcvl1d21pd2-5.47869D0*ETAcvl1d22pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine b972r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...             Becke 97-2  functional
c...  This is the Becke-97 functional reparametrized
c...
c...  P.J. Wilson, T.J. Bradley, and D.J.Tozer,
c...              J. Chem. Phys. 115, 9233 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumachia and its derivatives
c...
      sumachia = -7.4406D0*ETAcvl1chi2a2+1.40626D0*ETAcvl1chi2a1+9.99849
     1   D-1
      sumachiapchi2a = 1.40626D0*ETAcvl1chi2a1pchi2a-7.4406D0*ETAcvl1chi
     1   2a2pchi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 3.94796D-1*ETAcvl2chi2a2-6.91682D-1*ETAcvl2chi2a1+5.858
     1   08D-1
      sumbchiapchi2a = 3.94796D-1*ETAcvl2chi2a2pchi2a-6.91682D-1*ETAcvl2
     1   chi2a1pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.76125D0*ETAcvl3chi2a2+4.784D-2*ETAcvl3chi2a1+8.27642D
     1   -1
      sumcchiapchi2a = 1.76125D0*ETAcvl3chi2a2pchi2a+4.784D-2*ETAcvl3chi
     1   2a1pchi2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra1
     5   13
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine b972u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...             Becke 97-2  functional
c...  This is the Becke-97 functional reparametrized
c...
c...  P.J. Wilson, T.J. Bradley, and D.J.Tozer,
c...              J. Chem. Phys. 115, 9233 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 3.94796D-1*ETAcvl2chi2a2-6.91682D-1*ETAcvl2chi2a1+5.858
     1   08D-1
      sumbchiapchi2a = 3.94796D-1*ETAcvl2chi2a2pchi2a-6.91682D-1*ETAcvl2
     1   chi2a1pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.76125D0*ETAcvl3chi2a2+4.784D-2*ETAcvl3chi2a1+8.27642D
     1   -1
      sumcchiapchi2a = 1.76125D0*ETAcvl3chi2a2pchi2a+4.784D-2*ETAcvl3chi
     1   2a1pchi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 3.94796D-1*ETAcvl2chi2b2-6.91682D-1*ETAcvl2chi2b1+5.858
     1   08D-1
      sumbchibpchi2b = 3.94796D-1*ETAcvl2chi2b2pchi2b-6.91682D-1*ETAcvl2
     1   chi2b1pchi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 1.76125D0*ETAcvl3chi2b2+4.784D-2*ETAcvl3chi2b1+8.27642D
     1   -1
      sumcchibpchi2b = 1.76125D0*ETAcvl3chi2b2pchi2b+4.784D-2*ETAcvl3chi
     1   2b1pchi2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
c...
c...  sumad and its derivatives
c...
      sumad = -7.4406D0*ETAcvl1d22+1.40626D0*ETAcvl1d21+9.99849D-1
      sumadpd2 = 1.40626D0*ETAcvl1d21pd2-7.4406D0*ETAcvl1d22pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth93r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 93
c...
c...  F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N. C. Handy,
c...              J. Chem. Phys. 109, 6264 (1998)
c...
c...  Formulation taken from the Molpro manual (HCTH93)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.47857D0*ETAcvl1chi2a4+8.08564D0*ETAcvl1chi2a3-1.1543
     1   D1*ETAcvl1chi2a2+3.35287D0*ETAcvl1chi2a1+7.29974D-1
      sumachiapchi2a = -4.47857D0*ETAcvl1chi2a4pchi2a+8.08564D0*ETAcvl1c
     1   hi2a3pchi2a-1.1543D1*ETAcvl1chi2a2pchi2a+3.35287D0*ETAcvl1chi2a
     2   1pchi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 1.55396D0*ETAcvl2chi2a4-8.02496D-1*ETAcvl2chi2a3-1.2517
     1   D-2*ETAcvl2chi2a2-3.38622D-2*ETAcvl2chi2a1+2.22601D-1
      sumbchiapchi2a = 1.55396D0*ETAcvl2chi2a4pchi2a-8.02496D-1*ETAcvl2c
     1   hi2a3pchi2a-1.2517D-2*ETAcvl2chi2a2pchi2a-3.38622D-2*ETAcvl2chi
     2   2a1pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 4.49357D0*ETAcvl3chi2a4-6.78549D0*ETAcvl3chi2a3+5.5992D
     1   0*ETAcvl3chi2a2-7.44056D-1*ETAcvl3chi2a1+1.0932D0
      sumcchiapchi2a = 4.49357D0*ETAcvl3chi2a4pchi2a-6.78549D0*ETAcvl3ch
     1   i2a3pchi2a+5.5992D0*ETAcvl3chi2a2pchi2a-7.44056D-1*ETAcvl3chi2a
     2   1pchi2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra1
     5   13
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth93u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 93
c...
c...  F.A. Hamprecht, A.J. Cohen, D.J. Tozer, and N. C. Handy,
c...              J. Chem. Phys. 109, 6264 (1998)
c...
c...  Formulation taken from the Molpro manual (HCTH93)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 1.55396D0*ETAcvl2chi2a4-8.02496D-1*ETAcvl2chi2a3-1.2517
     1   D-2*ETAcvl2chi2a2-3.38622D-2*ETAcvl2chi2a1+2.22601D-1
      sumbchiapchi2a = 1.55396D0*ETAcvl2chi2a4pchi2a-8.02496D-1*ETAcvl2c
     1   hi2a3pchi2a-1.2517D-2*ETAcvl2chi2a2pchi2a-3.38622D-2*ETAcvl2chi
     2   2a1pchi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 4.49357D0*ETAcvl3chi2a4-6.78549D0*ETAcvl3chi2a3+5.5992D
     1   0*ETAcvl3chi2a2-7.44056D-1*ETAcvl3chi2a1+1.0932D0
      sumcchiapchi2a = 4.49357D0*ETAcvl3chi2a4pchi2a-6.78549D0*ETAcvl3ch
     1   i2a3pchi2a+5.5992D0*ETAcvl3chi2a2pchi2a-7.44056D-1*ETAcvl3chi2a
     2   1pchi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      chi2b3 = chi2b*chi2b2
      chi2b4 = chi2b*chi2b3
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 1.55396D0*ETAcvl2chi2b4-8.02496D-1*ETAcvl2chi2b3-1.2517
     1   D-2*ETAcvl2chi2b2-3.38622D-2*ETAcvl2chi2b1+2.22601D-1
      sumbchibpchi2b = 1.55396D0*ETAcvl2chi2b4pchi2b-8.02496D-1*ETAcvl2c
     1   hi2b3pchi2b-1.2517D-2*ETAcvl2chi2b2pchi2b-3.38622D-2*ETAcvl2chi
     2   2b1pchi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 4.49357D0*ETAcvl3chi2b4-6.78549D0*ETAcvl3chi2b3+5.5992D
     1   0*ETAcvl3chi2b2-7.44056D-1*ETAcvl3chi2b1+1.0932D0
      sumcchibpchi2b = 4.49357D0*ETAcvl3chi2b4pchi2b-6.78549D0*ETAcvl3ch
     1   i2b3pchi2b+5.5992D0*ETAcvl3chi2b2pchi2b-7.44056D-1*ETAcvl3chi2b
     2   1pchi2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      d23 = d2*d22
      d24 = d2*d23
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
c...
c...  sumad and its derivatives
c...
      sumad = -4.47857D0*ETAcvl1d24+8.08564D0*ETAcvl1d23-1.1543D1*ETAcvl
     1   1d22+3.35287D0*ETAcvl1d21+7.29974D-1
      sumadpd2 = -4.47857D0*ETAcvl1d24pd2+8.08564D0*ETAcvl1d23pd2-1.1543
     1   D1*ETAcvl1d22pd2+3.35287D0*ETAcvl1d21pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth120r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 120
c...
c...  A. D. Boese, N. L. Doltsinis, N. C. Handy and M. Sprik,
c...              J. Chem. Phys. 112, 1670 (2000)
c...
c...  Formulation taken from the Molpro manual (HCTH120)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumachia and its derivatives
c...
      sumachia = -1.1323D1*ETAcvl1chi2a4+2.311D1*ETAcvl1chi2a3-2.4707D1*
     1   ETAcvl1chi2a2+6.9298D0*ETAcvl1chi2a1+5.1473D-1
      sumachiapchi2a = -1.1323D1*ETAcvl1chi2a4pchi2a+2.311D1*ETAcvl1chi2
     1   a3pchi2a-2.4707D1*ETAcvl1chi2a2pchi2a+6.9298D0*ETAcvl1chi2a1pch
     2   i2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 2.4853D0*ETAcvl2chi2a4-1.9925D0*ETAcvl2chi2a3+4.329D-1*
     1   ETAcvl2chi2a2-2.607D-1*ETAcvl2chi2a1+4.8951D-1
      sumbchiapchi2a = 2.4853D0*ETAcvl2chi2a4pchi2a-1.9925D0*ETAcvl2chi2
     1   a3pchi2a+4.329D-1*ETAcvl2chi2a2pchi2a-2.607D-1*ETAcvl2chi2a1pch
     2   i2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.1717D0*ETAcvl3chi2a4-4.1075D0*ETAcvl3chi2a3+5.0783D0*
     1   ETAcvl3chi2a2-7.472D-1*ETAcvl3chi2a1+1.09163D0
      sumcchiapchi2a = 1.1717D0*ETAcvl3chi2a4pchi2a-4.1075D0*ETAcvl3chi2
     1   a3pchi2a+5.0783D0*ETAcvl3chi2a2pchi2a-7.472D-1*ETAcvl3chi2a1pch
     2   i2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra1
     5   13
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth120u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 120
c...
c...  A. D. Boese, N. L. Doltsinis, N. C. Handy and M. Sprik,
c...              J. Chem. Phys. 112, 1670 (2000)
c...
c...  Formulation taken from the Molpro manual (HCTH120)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 2.4853D0*ETAcvl2chi2a4-1.9925D0*ETAcvl2chi2a3+4.329D-1*
     1   ETAcvl2chi2a2-2.607D-1*ETAcvl2chi2a1+4.8951D-1
      sumbchiapchi2a = 2.4853D0*ETAcvl2chi2a4pchi2a-1.9925D0*ETAcvl2chi2
     1   a3pchi2a+4.329D-1*ETAcvl2chi2a2pchi2a-2.607D-1*ETAcvl2chi2a1pch
     2   i2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.1717D0*ETAcvl3chi2a4-4.1075D0*ETAcvl3chi2a3+5.0783D0*
     1   ETAcvl3chi2a2-7.472D-1*ETAcvl3chi2a1+1.09163D0
      sumcchiapchi2a = 1.1717D0*ETAcvl3chi2a4pchi2a-4.1075D0*ETAcvl3chi2
     1   a3pchi2a+5.0783D0*ETAcvl3chi2a2pchi2a-7.472D-1*ETAcvl3chi2a1pch
     2   i2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      chi2b3 = chi2b*chi2b2
      chi2b4 = chi2b*chi2b3
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 2.4853D0*ETAcvl2chi2b4-1.9925D0*ETAcvl2chi2b3+4.329D-1*
     1   ETAcvl2chi2b2-2.607D-1*ETAcvl2chi2b1+4.8951D-1
      sumbchibpchi2b = 2.4853D0*ETAcvl2chi2b4pchi2b-1.9925D0*ETAcvl2chi2
     1   b3pchi2b+4.329D-1*ETAcvl2chi2b2pchi2b-2.607D-1*ETAcvl2chi2b1pch
     2   i2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 1.1717D0*ETAcvl3chi2b4-4.1075D0*ETAcvl3chi2b3+5.0783D0*
     1   ETAcvl3chi2b2-7.472D-1*ETAcvl3chi2b1+1.09163D0
      sumcchibpchi2b = 1.1717D0*ETAcvl3chi2b4pchi2b-4.1075D0*ETAcvl3chi2
     1   b3pchi2b+5.0783D0*ETAcvl3chi2b2pchi2b-7.472D-1*ETAcvl3chi2b1pch
     2   i2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      d23 = d2*d22
      d24 = d2*d23
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
c...
c...  sumad and its derivatives
c...
      sumad = -1.1323D1*ETAcvl1d24+2.311D1*ETAcvl1d23-2.4707D1*ETAcvl1d2
     1   2+6.9298D0*ETAcvl1d21+5.1473D-1
      sumadpd2 = -1.1323D1*ETAcvl1d24pd2+2.311D1*ETAcvl1d23pd2-2.4707D1*
     1   ETAcvl1d22pd2+6.9298D0*ETAcvl1d21pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth147r(ra,saa,ec,dra,dsaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 147
c...
c...  A. D. Boese, N. L. Doltsinis, N. C. Handy and M. Sprik,
c...              J. Chem. Phys. 112, 1670 (2000)
c...
c...  Formulation taken from the Molpro manual (HCTH147)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  saa     alpha density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 eablc
      real*8 eaa,eaapra
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 junk
      real*8 r
      real*8 ra3,ra23,ra43,ra83,ram83,ra113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcu(ra,ra,eablc,eaapra,junk)
      eaa = eablc*ra
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      ra3 = ra**3
      ra23 = ra**tt
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl1chi2ap1 = 6.0D-3*chi2a+1.0D0
      cvl1chi2ap12 = cvl1chi2ap1**2
      cvl1chi2ap13 = cvl1chi2ap1*cvl1chi2ap12
      ETAcvl1chi2a1 = 6.0D-3*chi2a/cvl1chi2ap1
      ETAcvl1chi2a1pchi2a = 6.0D-3/cvl1chi2ap12
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumachia and its derivatives
c...
      sumachia = -2.0428D1*ETAcvl1chi2a4+3.5033D1*ETAcvl1chi2a3-2.8382D1
     1   *ETAcvl1chi2a2+7.0146D0*ETAcvl1chi2a1+5.4235D-1
      sumachiapchi2a = -2.0428D1*ETAcvl1chi2a4pchi2a+3.5033D1*ETAcvl1chi
     1   2a3pchi2a-2.8382D1*ETAcvl1chi2a2pchi2a+7.0146D0*ETAcvl1chi2a1pc
     2   hi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 8.854D-1*ETAcvl2chi2a4+1.0575D0*ETAcvl2chi2a3-1.3064D0*
     1   ETAcvl2chi2a2-1.71D-2*ETAcvl2chi2a1+5.6258D-1
      sumbchiapchi2a = 8.854D-1*ETAcvl2chi2a4pchi2a+1.0575D0*ETAcvl2chi2
     1   a3pchi2a-1.3064D0*ETAcvl2chi2a2pchi2a-1.71D-2*ETAcvl2chi2a1pchi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 3.0454D0*ETAcvl3chi2a4-5.8676D0*ETAcvl3chi2a3+5.5721D0*
     1   ETAcvl3chi2a2-7.992D-1*ETAcvl3chi2a1+1.09025D0
      sumcchiapchi2a = 3.0454D0*ETAcvl3chi2a4pchi2a-5.8676D0*ETAcvl3chi2
     1   a3pchi2a+5.5721D0*ETAcvl3chi2a2pchi2a-7.992D-1*ETAcvl3chi2a1pch
     2   i2a
c...
c...  We compute the output quantities
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*ra43*sumcchia+ea0*sumbchia+(eaa-ea0)
     1   *sumachia)/ra
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-1.8599157
     2   62415402D1*(eaa-ea0)*sumachiapchi2a)-3.0D0*ra3*(2.8844991
     3   40614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra1
     5   13
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth147u(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 147
c...
c...  A. D. Boese, N. L. Doltsinis, N. C. Handy and M. Sprik,
c...              J. Chem. Phys. 112, 1670 (2000)
c...
c...  Formulation taken from the Molpro manual (HCTH147)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first derivatives.
c...
c...  Input parameters:
c...
c...  ra      alpha density
c...  rb      beta density
c...  saa     alpha density gradient invariant
c...  sbb     beta density gradient invariant
c...
c...  Output parameters:
c...
c...  ec      contribution to exchange-correlation energy
c...          (i.e.: functional value dived by the total density)
c...  dra     functional derivative with respect to alpha density
c...  drb     funct. deriv. w.r.t. beta density
c...  dsaa    funct. deriv. w.r.t. alpha gradient invariant
c...  dsbb    funct. deriv. w.r.t. alpha gradient invariant
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 eablc,deablcdra,deablcdrb
      real*8 eab,eabpra,eabprb
      real*8 ea0lc,dea0lcdra
      real*8 ea0,ea0pra
      real*8 eb0lc,deb0lcdrb
      real*8 eb0,eb0prb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra23,ra43,ra83,ram83,ra113
      real*8 rb2,rb3,rb23,rb43,rb83,rbm83,rb113
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2
      real*8 sumad,sumadpd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a
      real*8 sumbchia,sumbchiapchi2a
      real*8 sumcchia,sumcchiapchi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b
      real*8 sumbchib,sumbchibpchi2b
      real*8 sumcchib,sumcchibpchi2b
      real*8 zero,one,third
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0)
      real*8 epsi
      parameter (epsi=1.0d-15)
c...
c...  First, let's get the local terms out of the way.
c...
c...  We compute only one term (either ea0 or eb0) or
c...  all three, according to the values of ra and rb
c...
      r = rb+ra
c...
c...  ea0 only
c...
      if(rb.lt.epsi)then
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcu(ra,rb,eablc,eabpra,eabprb)
      eab = eablc*r
      call pw91lcu(ra,zero,ea0lc,ea0pra,junk)
      ea0 = ea0lc*ra
      call pw91lcu(rb,zero,eb0lc,eb0prb,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra23 = ra2**THIRD
      ra43 = ra23**2
      ra83 = ra43**2
      ram83 = one/ra83
      ra113 = ra*ra83
c...
c...  chi2a
c...
      chi2a = ram83*saa
c...
c...  the functions eta of chi2a and their powers
c...
      chi2a2 = chi2a**2
      chi2a3 = chi2a*chi2a2
      chi2a4 = chi2a*chi2a3
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 8.854D-1*ETAcvl2chi2a4+1.0575D0*ETAcvl2chi2a3-1.3064D0*
     1   ETAcvl2chi2a2-1.71D-2*ETAcvl2chi2a1+5.6258D-1
      sumbchiapchi2a = 8.854D-1*ETAcvl2chi2a4pchi2a+1.0575D0*ETAcvl2chi2
     1   a3pchi2a-1.3064D0*ETAcvl2chi2a2pchi2a-1.71D-2*ETAcvl2chi2a1pchi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 3.0454D0*ETAcvl3chi2a4-5.8676D0*ETAcvl3chi2a3+5.5721D0*
     1   ETAcvl3chi2a2-7.992D-1*ETAcvl3chi2a1+1.09025D0
      sumcchiapchi2a = 3.0454D0*ETAcvl3chi2a4pchi2a-5.8676D0*ETAcvl3chi2
     1   a3pchi2a+5.5721D0*ETAcvl3chi2a2pchi2a-7.992D-1*ETAcvl3chi2a1pch
     2   i2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb23 = rb2**THIRD
      rb43 = rb23**2
      rb83 = rb43**2
      rbm83 = one/rb83
      rb113 = rb*rb83
c...
c...  chi2b
c...
      chi2b = rbm83*sbb
c...
c...  the functions eta of chi2b and their powers
c...
      chi2b2 = chi2b**2
      chi2b3 = chi2b*chi2b2
      chi2b4 = chi2b*chi2b3
      cvl2chi2bp1 = 2.0D-1*chi2b+1.0D0
      cvl2chi2bp12 = cvl2chi2bp1**2
      cvl2chi2bp13 = cvl2chi2bp1*cvl2chi2bp12
      ETAcvl2chi2b1 = 2.0D-1*chi2b/cvl2chi2bp1
      ETAcvl2chi2b1pchi2b = 2.0D-1/cvl2chi2bp12
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 8.854D-1*ETAcvl2chi2b4+1.0575D0*ETAcvl2chi2b3-1.3064D0*
     1   ETAcvl2chi2b2-1.71D-2*ETAcvl2chi2b1+5.6258D-1
      sumbchibpchi2b = 8.854D-1*ETAcvl2chi2b4pchi2b+1.0575D0*ETAcvl2chi2
     1   b3pchi2b-1.3064D0*ETAcvl2chi2b2pchi2b-1.71D-2*ETAcvl2chi2b1pchi
     2   2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 3.0454D0*ETAcvl3chi2b4-5.8676D0*ETAcvl3chi2b3+5.5721D0*
     1   ETAcvl3chi2b2-7.992D-1*ETAcvl3chi2b1+1.09025D0
      sumcchibpchi2b = 3.0454D0*ETAcvl3chi2b4pchi2b-5.8676D0*ETAcvl3chi2
     1   b3pchi2b+5.5721D0*ETAcvl3chi2b2pchi2b-7.992D-1*ETAcvl3chi2b1pch
     2   i2b
      endif
c...
c...  these quantities depend on bot ra and rb, thus no need
c...  for testing, as one of the two must be non zero,
c...  otherwise we would not be executing this code
c...
c...  auxiliary variable d2
c...
      d2 = 5.0D-1*(chi2b+chi2a)
c...
c...  the function eta of d2 and its powers
c...
      d22 = d2**2
      d23 = d2*d22
      d24 = d2*d23
      cvl1d2p1 = 6.0D-3*d2+1.0D0
      cvl1d2p12 = cvl1d2p1**2
      cvl1d2p13 = cvl1d2p1*cvl1d2p12
      ETAcvl1d21 = 6.0D-3*d2/cvl1d2p1
      ETAcvl1d21pd2 = 6.0D-3/cvl1d2p12
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
c...
c...  sumad and its derivatives
c...
      sumad = -2.0428D1*ETAcvl1d24+3.5033D1*ETAcvl1d23-2.8382D1*ETAcvl1d
     1   22+7.0146D0*ETAcvl1d21+5.4235D-1
      sumadpd2 = -2.0428D1*ETAcvl1d24pd2+3.5033D1*ETAcvl1d23pd2-2.8382D1
     1   *ETAcvl1d22pd2+7.0146D0*ETAcvl1d21pd2
c...
c...  We compute the output quantities, using different
c...  functional forms according to the values of ra and rb
c...
c...  Case 1:  rb is zero
c...
      if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (ea0*sumbchia-9.305257363491D-1*ra43*sumcchia)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(4.0D0*saa*(4.326748710922225D0*ra43*su
     1   mcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)-3.0D0*ra3*
     2   (2.884499140614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra2
     3   3*sumbchia))/ra113
      drb = ZERO
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)
      dsbb = ZERO
c...
c...  Case 2:  ra is zero
c...
      else if(rb.lt.epsi)then
c...
c...  functional value (divided by the total density)
c...
      ec = (eb0*sumbchib-9.305257363491D-1*rb43*sumcchib)/r
c...
c...  First derivatives of the functional
c...
      dra = ZERO
      drb = 1.433756689713499D-1*(4.0D0*sbb*(4.326748710922225D0*rb43*su
     1   mcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)-3.0D0*rb3*
     2   (2.884499140614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb2
     3   3*sumbchib))/rb113
      dsaa = ZERO
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)
c...
c...  Case 3:  both ra and rb are nonzero
c...
      else
c...
c...  functional value (divided by the total density)
c...
      ec = (-9.305257363491D-1*rb43*sumcchib-9.305257363491D-1*ra43*sumc
     1   chia+eb0*sumbchib+ea0*sumbchia+(-eb0+eab-ea0)*sumad)/r
c...
c...  First derivatives of the functional
c...
      dra = 1.433756689713499D-1*(saa*(4.0D0*(4.326748710922225D0*ra43*s
     1   umcchiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*ra3*(2.884499140
     3   614817D0*ra*sumcchia-2.324894703019253D0*ea0pra*ra23*sumbchia)+
     4   6.974684109057759D0*(eabpra-ea0pra)*ra113*sumad)/ra113
      drb = 1.433756689713499D-1*(sbb*(4.0D0*(4.326748710922225D0*rb43*s
     1   umcchibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b)+9.2995788
     2   12077012D0*(eb0-eab+ea0)*sumadpd2)-3.0D0*rb3*(2.884499140
     3   614817D0*rb*sumcchib-2.324894703019253D0*eb0prb*rb23*sumbchib)-
     4   6.974684109057759D0*(eb0prb-eabprb)*rb113*sumad)/rb113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      dsbb = -2.150635034570249D-1*rbm83*(4.326748710922225D0*rb43*sumcc
     1   hibpchi2b-4.649789406038506D0*eb0*sumbchibpchi2b+2.324894703019
     2   253D0*(eb0-eab+ea0)*sumadpd2)
      endif
c...
      return
      end

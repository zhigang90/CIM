c=======================================================================
      subroutine b97dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.5961D0*ETAcvl1chi2a2+7.471D-1*ETAcvl1chi2a1+9.454D-1
      sumachiapchi2a = 7.471D-1*ETAcvl1chi2a1pchi2a-4.5961D0*ETAcvl1chi2
     1   a2pchi2a
      sumachiaschi2a = 7.471D-1*ETAcvl1chi2a1schi2a-4.5961D0*ETAcvl1chi2
     1   a2schi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.4868D0*ETAcvl2chi2a2+2.3487D0*ETAcvl2chi2a1+1.737D-1
      sumbchiapchi2a = 2.3487D0*ETAcvl2chi2a1pchi2a-2.4868D0*ETAcvl2chi2
     1   a2pchi2a
      sumbchiaschi2a = 2.3487D0*ETAcvl2chi2a1schi2a-2.4868D0*ETAcvl2chi2
     1   a2schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 7.481D-1*ETAcvl3chi2a2+5.073D-1*ETAcvl3chi2a1+8.094D-1
      sumcchiapchi2a = 7.481D-1*ETAcvl3chi2a2pchi2a+5.073D-1*ETAcvl3chi2
     1   a1pchi2a
      sumcchiaschi2a = 7.481D-1*ETAcvl3chi2a2schi2a+5.073D-1*ETAcvl3chi2
     1   a1schi2a
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-1.0D0*ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine b97du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...           Becke 97 (B97) functional
c...
c...  A. D. Becke, J. Chem. Phys. 107, 8554 (1997)
c...
c...  Formulation taken from the Molpro manual (D97)
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.4868D0*ETAcvl2chi2a2+2.3487D0*ETAcvl2chi2a1+1.737D-1
      sumbchiapchi2a = 2.3487D0*ETAcvl2chi2a1pchi2a-2.4868D0*ETAcvl2chi2
     1   a2pchi2a
      sumbchiaschi2a = 2.3487D0*ETAcvl2chi2a1schi2a-2.4868D0*ETAcvl2chi2
     1   a2schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 7.481D-1*ETAcvl3chi2a2+5.073D-1*ETAcvl3chi2a1+8.094D-1
      sumcchiapchi2a = 7.481D-1*ETAcvl3chi2a2pchi2a+5.073D-1*ETAcvl3chi2
     1   a1pchi2a
      sumcchiaschi2a = 7.481D-1*ETAcvl3chi2a2schi2a+5.073D-1*ETAcvl3chi2
     1   a1schi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
c...
c...  sumbchib and its derivatives
c...
      sumbchib = -2.4868D0*ETAcvl2chi2b2+2.3487D0*ETAcvl2chi2b1+1.737D-1
      sumbchibpchi2b = 2.3487D0*ETAcvl2chi2b1pchi2b-2.4868D0*ETAcvl2chi2
     1   b2pchi2b
      sumbchibschi2b = 2.3487D0*ETAcvl2chi2b1schi2b-2.4868D0*ETAcvl2chi2
     1   b2schi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 7.481D-1*ETAcvl3chi2b2+5.073D-1*ETAcvl3chi2b1+8.094D-1
      sumcchibpchi2b = 7.481D-1*ETAcvl3chi2b2pchi2b+5.073D-1*ETAcvl3chi2
     1   b1pchi2b
      sumcchibschi2b = 7.481D-1*ETAcvl3chi2b2schi2b+5.073D-1*ETAcvl3chi2
     1   b1schi2b
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
c...
c...  sumad and its derivatives
c...
      sumad = -4.5961D0*ETAcvl1d22+7.471D-1*ETAcvl1d21+9.454D-1
      sumadpd2 = 7.471D-1*ETAcvl1d21pd2-4.5961D0*ETAcvl1d22pd2
      sumadsd2 = 7.471D-1*ETAcvl1d21sd2-4.5961D0*ETAcvl1d22sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine b971dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumachia and its derivatives
c...
      sumachia = -5.47869D0*ETAcvl1chi2a2+7.88552D-1*ETAcvl1chi2a1+9.556
     1   89D-1
      sumachiapchi2a = 7.88552D-1*ETAcvl1chi2a1pchi2a-5.47869D0*ETAcvl1c
     1   hi2a2pchi2a
      sumachiaschi2a = 7.88552D-1*ETAcvl1chi2a1schi2a-5.47869D0*ETAcvl1c
     1   hi2a2schi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.87103D0*ETAcvl2chi2a2+2.71681D0*ETAcvl2chi2a1+8.2001
     1   1D-2
      sumbchiapchi2a = 2.71681D0*ETAcvl2chi2a1pchi2a-2.87103D0*ETAcvl2ch
     1   i2a2pchi2a
      sumbchiaschi2a = 2.71681D0*ETAcvl2chi2a1schi2a-2.87103D0*ETAcvl2ch
     1   i2a2schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 6.60975D-1*ETAcvl3chi2a2+5.73805D-1*ETAcvl3chi2a1+7.895
     1   18D-1
      sumcchiapchi2a = 6.60975D-1*ETAcvl3chi2a2pchi2a+5.73805D-1*ETAcvl3
     1   chi2a1pchi2a
      sumcchiaschi2a = 6.60975D-1*ETAcvl3chi2a2schi2a+5.73805D-1*ETAcvl3
     1   chi2a1schi2a
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-1.0D0*ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine b971du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumbchia and its derivatives
c...
      sumbchia = -2.87103D0*ETAcvl2chi2a2+2.71681D0*ETAcvl2chi2a1+8.2001
     1   1D-2
      sumbchiapchi2a = 2.71681D0*ETAcvl2chi2a1pchi2a-2.87103D0*ETAcvl2ch
     1   i2a2pchi2a
      sumbchiaschi2a = 2.71681D0*ETAcvl2chi2a1schi2a-2.87103D0*ETAcvl2ch
     1   i2a2schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 6.60975D-1*ETAcvl3chi2a2+5.73805D-1*ETAcvl3chi2a1+7.895
     1   18D-1
      sumcchiapchi2a = 6.60975D-1*ETAcvl3chi2a2pchi2a+5.73805D-1*ETAcvl3
     1   chi2a1pchi2a
      sumcchiaschi2a = 6.60975D-1*ETAcvl3chi2a2schi2a+5.73805D-1*ETAcvl3
     1   chi2a1schi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
c...
c...  sumbchib and its derivatives
c...
      sumbchib = -2.87103D0*ETAcvl2chi2b2+2.71681D0*ETAcvl2chi2b1+8.2001
     1   1D-2
      sumbchibpchi2b = 2.71681D0*ETAcvl2chi2b1pchi2b-2.87103D0*ETAcvl2ch
     1   i2b2pchi2b
      sumbchibschi2b = 2.71681D0*ETAcvl2chi2b1schi2b-2.87103D0*ETAcvl2ch
     1   i2b2schi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 6.60975D-1*ETAcvl3chi2b2+5.73805D-1*ETAcvl3chi2b1+7.895
     1   18D-1
      sumcchibpchi2b = 6.60975D-1*ETAcvl3chi2b2pchi2b+5.73805D-1*ETAcvl3
     1   chi2b1pchi2b
      sumcchibschi2b = 6.60975D-1*ETAcvl3chi2b2schi2b+5.73805D-1*ETAcvl3
     1   chi2b1schi2b
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
c...
c...  sumad and its derivatives
c...
      sumad = -5.47869D0*ETAcvl1d22+7.88552D-1*ETAcvl1d21+9.55689D-1
      sumadpd2 = 7.88552D-1*ETAcvl1d21pd2-5.47869D0*ETAcvl1d22pd2
      sumadsd2 = 7.88552D-1*ETAcvl1d21sd2-5.47869D0*ETAcvl1d22sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-1.0D0*eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-1.0D0*ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-1.0D0*eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine b972dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...             Becke 97-2  functional
c...  This is the Becke-97 functional reparametrized
c...
c...  P.J. Wilson, T.J. Dradley, and D.J.Tozer,
c...              J. Chem. Phys. 115, 9233 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumachia and its derivatives
c...
      sumachia = -7.4406D0*ETAcvl1chi2a2+1.40626D0*ETAcvl1chi2a1+9.99849
     1   D-1
      sumachiapchi2a = 1.40626D0*ETAcvl1chi2a1pchi2a-7.4406D0*ETAcvl1chi
     1   2a2pchi2a
      sumachiaschi2a = 1.40626D0*ETAcvl1chi2a1schi2a-7.4406D0*ETAcvl1chi
     1   2a2schi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 3.94796D-1*ETAcvl2chi2a2-6.91682D-1*ETAcvl2chi2a1+5.858
     1   08D-1
      sumbchiapchi2a = 3.94796D-1*ETAcvl2chi2a2pchi2a-6.91682D-1*ETAcvl2
     1   chi2a1pchi2a
      sumbchiaschi2a = 3.94796D-1*ETAcvl2chi2a2schi2a-6.91682D-1*ETAcvl2
     1   chi2a1schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.76125D0*ETAcvl3chi2a2+4.784D-2*ETAcvl3chi2a1+8.27642D
     1   -1
      sumcchiapchi2a = 1.76125D0*ETAcvl3chi2a2pchi2a+4.784D-2*ETAcvl3chi
     1   2a1pchi2a
      sumcchiaschi2a = 1.76125D0*ETAcvl3chi2a2schi2a+4.784D-2*ETAcvl3chi
     1   2a1schi2a
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine b972du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 3.94796D-1*ETAcvl2chi2a2-6.91682D-1*ETAcvl2chi2a1+5.858
     1   08D-1
      sumbchiapchi2a = 3.94796D-1*ETAcvl2chi2a2pchi2a-6.91682D-1*ETAcvl2
     1   chi2a1pchi2a
      sumbchiaschi2a = 3.94796D-1*ETAcvl2chi2a2schi2a-6.91682D-1*ETAcvl2
     1   chi2a1schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.76125D0*ETAcvl3chi2a2+4.784D-2*ETAcvl3chi2a1+8.27642D
     1   -1
      sumcchiapchi2a = 1.76125D0*ETAcvl3chi2a2pchi2a+4.784D-2*ETAcvl3chi
     1   2a1pchi2a
      sumcchiaschi2a = 1.76125D0*ETAcvl3chi2a2schi2a+4.784D-2*ETAcvl3chi
     1   2a1schi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 3.94796D-1*ETAcvl2chi2b2-6.91682D-1*ETAcvl2chi2b1+5.858
     1   08D-1
      sumbchibpchi2b = 3.94796D-1*ETAcvl2chi2b2pchi2b-6.91682D-1*ETAcvl2
     1   chi2b1pchi2b
      sumbchibschi2b = 3.94796D-1*ETAcvl2chi2b2schi2b-6.91682D-1*ETAcvl2
     1   chi2b1schi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 1.76125D0*ETAcvl3chi2b2+4.784D-2*ETAcvl3chi2b1+8.27642D
     1   -1
      sumcchibpchi2b = 1.76125D0*ETAcvl3chi2b2pchi2b+4.784D-2*ETAcvl3chi
     1   2b1pchi2b
      sumcchibschi2b = 1.76125D0*ETAcvl3chi2b2schi2b+4.784D-2*ETAcvl3chi
     1   2b1schi2b
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
c...
c...  sumad and its derivatives
c...
      sumad = -7.4406D0*ETAcvl1d22+1.40626D0*ETAcvl1d21+9.99849D-1
      sumadpd2 = 1.40626D0*ETAcvl1d21pd2-7.4406D0*ETAcvl1d22pd2
      sumadsd2 = 1.40626D0*ETAcvl1d21sd2-7.4406D0*ETAcvl1d22sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth93dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a,ETAcvl1chi2a3schi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a,ETAcvl1chi2a4schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      ETAcvl1chi2a3schi2a = -1.296D-6*chi2a*(6.0D-3*chi2a-1.0D0)/cvl1chi
     1   2ap15
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      ETAcvl1chi2a4schi2a = -5.184D-9*(1.2D-2*chi2a-3.0D0)*chi2a2/cvl1ch
     1   i2ap16
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.47857D0*ETAcvl1chi2a4+8.08564D0*ETAcvl1chi2a3-1.1543
     1   D1*ETAcvl1chi2a2+3.35287D0*ETAcvl1chi2a1+7.29974D-1
      sumachiapchi2a = -4.47857D0*ETAcvl1chi2a4pchi2a+8.08564D0*ETAcvl1c
     1   hi2a3pchi2a-1.1543D1*ETAcvl1chi2a2pchi2a+3.35287D0*ETAcvl1chi2a
     2   1pchi2a
      sumachiaschi2a = -4.47857D0*ETAcvl1chi2a4schi2a+8.08564D0*ETAcvl1c
     1   hi2a3schi2a-1.1543D1*ETAcvl1chi2a2schi2a+3.35287D0*ETAcvl1chi2a
     2   1schi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 1.55396D0*ETAcvl2chi2a4-8.02496D-1*ETAcvl2chi2a3-1.2517
     1   D-2*ETAcvl2chi2a2-3.38622D-2*ETAcvl2chi2a1+2.22601D-1
      sumbchiapchi2a = 1.55396D0*ETAcvl2chi2a4pchi2a-8.02496D-1*ETAcvl2c
     1   hi2a3pchi2a-1.2517D-2*ETAcvl2chi2a2pchi2a-3.38622D-2*ETAcvl2chi
     2   2a1pchi2a
      sumbchiaschi2a = 1.55396D0*ETAcvl2chi2a4schi2a-8.02496D-1*ETAcvl2c
     1   hi2a3schi2a-1.2517D-2*ETAcvl2chi2a2schi2a-3.38622D-2*ETAcvl2chi
     2   2a1schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 4.49357D0*ETAcvl3chi2a4-6.78549D0*ETAcvl3chi2a3+5.5992D
     1   0*ETAcvl3chi2a2-7.44056D-1*ETAcvl3chi2a1+1.0932D0
      sumcchiapchi2a = 4.49357D0*ETAcvl3chi2a4pchi2a-6.78549D0*ETAcvl3ch
     1   i2a3pchi2a+5.5992D0*ETAcvl3chi2a2pchi2a-7.44056D-1*ETAcvl3chi2a
     2   1pchi2a
      sumcchiaschi2a = 4.49357D0*ETAcvl3chi2a4schi2a-6.78549D0*ETAcvl3ch
     1   i2a3schi2a+5.5992D0*ETAcvl3chi2a2schi2a-7.44056D-1*ETAcvl3chi2a
     2   1schi2a
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth93du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2,ETAcvl1d23sd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2,ETAcvl1d24sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b,ETAcvl2chi2b3schi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b,ETAcvl2chi2b4schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b,ETAcvl3chi2b3schi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b,ETAcvl3chi2b4schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 1.55396D0*ETAcvl2chi2a4-8.02496D-1*ETAcvl2chi2a3-1.2517
     1   D-2*ETAcvl2chi2a2-3.38622D-2*ETAcvl2chi2a1+2.22601D-1
      sumbchiapchi2a = 1.55396D0*ETAcvl2chi2a4pchi2a-8.02496D-1*ETAcvl2c
     1   hi2a3pchi2a-1.2517D-2*ETAcvl2chi2a2pchi2a-3.38622D-2*ETAcvl2chi
     2   2a1pchi2a
      sumbchiaschi2a = 1.55396D0*ETAcvl2chi2a4schi2a-8.02496D-1*ETAcvl2c
     1   hi2a3schi2a-1.2517D-2*ETAcvl2chi2a2schi2a-3.38622D-2*ETAcvl2chi
     2   2a1schi2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 4.49357D0*ETAcvl3chi2a4-6.78549D0*ETAcvl3chi2a3+5.5992D
     1   0*ETAcvl3chi2a2-7.44056D-1*ETAcvl3chi2a1+1.0932D0
      sumcchiapchi2a = 4.49357D0*ETAcvl3chi2a4pchi2a-6.78549D0*ETAcvl3ch
     1   i2a3pchi2a+5.5992D0*ETAcvl3chi2a2pchi2a-7.44056D-1*ETAcvl3chi2a
     2   1pchi2a
      sumcchiaschi2a = 4.49357D0*ETAcvl3chi2a4schi2a-6.78549D0*ETAcvl3ch
     1   i2a3schi2a+5.5992D0*ETAcvl3chi2a2schi2a-7.44056D-1*ETAcvl3chi2a
     2   1schi2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      ETAcvl2chi2b3schi2b = -4.8D-2*chi2b*(2.0D-1*chi2b-1.0D0)/cvl2chi2b
     1   p15
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      ETAcvl2chi2b4schi2b = -6.4D-3*(4.0D-1*chi2b-3.0D0)*chi2b2/cvl2chi2
     1   bp16
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      ETAcvl3chi2b3schi2b = -3.84D-7*chi2b*(4.0D-3*chi2b-1.0D0)/cvl3chi2
     1   bp15
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
      ETAcvl3chi2b4schi2b = -1.024D-9*(8.0D-3*chi2b-3.0D0)*chi2b2/cvl3ch
     1   i2bp16
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 1.55396D0*ETAcvl2chi2b4-8.02496D-1*ETAcvl2chi2b3-1.2517
     1   D-2*ETAcvl2chi2b2-3.38622D-2*ETAcvl2chi2b1+2.22601D-1
      sumbchibpchi2b = 1.55396D0*ETAcvl2chi2b4pchi2b-8.02496D-1*ETAcvl2c
     1   hi2b3pchi2b-1.2517D-2*ETAcvl2chi2b2pchi2b-3.38622D-2*ETAcvl2chi
     2   2b1pchi2b
      sumbchibschi2b = 1.55396D0*ETAcvl2chi2b4schi2b-8.02496D-1*ETAcvl2c
     1   hi2b3schi2b-1.2517D-2*ETAcvl2chi2b2schi2b-3.38622D-2*ETAcvl2chi
     2   2b1schi2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 4.49357D0*ETAcvl3chi2b4-6.78549D0*ETAcvl3chi2b3+5.5992D
     1   0*ETAcvl3chi2b2-7.44056D-1*ETAcvl3chi2b1+1.0932D0
      sumcchibpchi2b = 4.49357D0*ETAcvl3chi2b4pchi2b-6.78549D0*ETAcvl3ch
     1   i2b3pchi2b+5.5992D0*ETAcvl3chi2b2pchi2b-7.44056D-1*ETAcvl3chi2b
     2   1pchi2b
      sumcchibschi2b = 4.49357D0*ETAcvl3chi2b4schi2b-6.78549D0*ETAcvl3ch
     1   i2b3schi2b+5.5992D0*ETAcvl3chi2b2schi2b-7.44056D-1*ETAcvl3chi2b
     2   1schi2b
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      ETAcvl1d23sd2 = -1.296D-6*d2*(6.0D-3*d2-1.0D0)/cvl1d2p15
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
      ETAcvl1d24sd2 = -5.184D-9*(1.2D-2*d2-3.0D0)*d22/cvl1d2p16
c...
c...  sumad and its derivatives
c...
      sumad = -4.47857D0*ETAcvl1d24+8.08564D0*ETAcvl1d23-1.1543D1*ETAcvl
     1   1d22+3.35287D0*ETAcvl1d21+7.29974D-1
      sumadpd2 = -4.47857D0*ETAcvl1d24pd2+8.08564D0*ETAcvl1d23pd2-1.1543
     1   D1*ETAcvl1d22pd2+3.35287D0*ETAcvl1d21pd2
      sumadsd2 = -4.47857D0*ETAcvl1d24sd2+8.08564D0*ETAcvl1d23sd2-1.1543
     1   D1*ETAcvl1d22sd2+3.35287D0*ETAcvl1d21sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-1.0D0*eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth120dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a,ETAcvl1chi2a3schi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a,ETAcvl1chi2a4schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      ETAcvl1chi2a3schi2a = -1.296D-6*chi2a*(6.0D-3*chi2a-1.0D0)/cvl1chi
     1   2ap15
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      ETAcvl1chi2a4schi2a = -5.184D-9*(1.2D-2*chi2a-3.0D0)*chi2a2/cvl1ch
     1   i2ap16
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumachia and its derivatives
c...
      sumachia = -1.1323D1*ETAcvl1chi2a4+2.311D1*ETAcvl1chi2a3-2.4707D1*
     1   ETAcvl1chi2a2+6.9298D0*ETAcvl1chi2a1+5.1473D-1
      sumachiapchi2a = -1.1323D1*ETAcvl1chi2a4pchi2a+2.311D1*ETAcvl1chi2
     1   a3pchi2a-2.4707D1*ETAcvl1chi2a2pchi2a+6.9298D0*ETAcvl1chi2a1pch
     2   i2a
      sumachiaschi2a = -1.1323D1*ETAcvl1chi2a4schi2a+2.311D1*ETAcvl1chi2
     1   a3schi2a-2.4707D1*ETAcvl1chi2a2schi2a+6.9298D0*ETAcvl1chi2a1sch
     2   i2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 2.4853D0*ETAcvl2chi2a4-1.9925D0*ETAcvl2chi2a3+4.329D-1*
     1   ETAcvl2chi2a2-2.607D-1*ETAcvl2chi2a1+4.8951D-1
      sumbchiapchi2a = 2.4853D0*ETAcvl2chi2a4pchi2a-1.9925D0*ETAcvl2chi2
     1   a3pchi2a+4.329D-1*ETAcvl2chi2a2pchi2a-2.607D-1*ETAcvl2chi2a1pch
     2   i2a
      sumbchiaschi2a = 2.4853D0*ETAcvl2chi2a4schi2a-1.9925D0*ETAcvl2chi2
     1   a3schi2a+4.329D-1*ETAcvl2chi2a2schi2a-2.607D-1*ETAcvl2chi2a1sch
     2   i2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.1717D0*ETAcvl3chi2a4-4.1075D0*ETAcvl3chi2a3+5.0783D0*
     1   ETAcvl3chi2a2-7.472D-1*ETAcvl3chi2a1+1.09163D0
      sumcchiapchi2a = 1.1717D0*ETAcvl3chi2a4pchi2a-4.1075D0*ETAcvl3chi2
     1   a3pchi2a+5.0783D0*ETAcvl3chi2a2pchi2a-7.472D-1*ETAcvl3chi2a1pch
     2   i2a
      sumcchiaschi2a = 1.1717D0*ETAcvl3chi2a4schi2a-4.1075D0*ETAcvl3chi2
     1   a3schi2a+5.0783D0*ETAcvl3chi2a2schi2a-7.472D-1*ETAcvl3chi2a1sch
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
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth120du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2,ETAcvl1d23sd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2,ETAcvl1d24sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b,ETAcvl2chi2b3schi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b,ETAcvl2chi2b4schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b,ETAcvl3chi2b3schi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b,ETAcvl3chi2b4schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 2.4853D0*ETAcvl2chi2a4-1.9925D0*ETAcvl2chi2a3+4.329D-1*
     1   ETAcvl2chi2a2-2.607D-1*ETAcvl2chi2a1+4.8951D-1
      sumbchiapchi2a = 2.4853D0*ETAcvl2chi2a4pchi2a-1.9925D0*ETAcvl2chi2
     1   a3pchi2a+4.329D-1*ETAcvl2chi2a2pchi2a-2.607D-1*ETAcvl2chi2a1pch
     2   i2a
      sumbchiaschi2a = 2.4853D0*ETAcvl2chi2a4schi2a-1.9925D0*ETAcvl2chi2
     1   a3schi2a+4.329D-1*ETAcvl2chi2a2schi2a-2.607D-1*ETAcvl2chi2a1sch
     2   i2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 1.1717D0*ETAcvl3chi2a4-4.1075D0*ETAcvl3chi2a3+5.0783D0*
     1   ETAcvl3chi2a2-7.472D-1*ETAcvl3chi2a1+1.09163D0
      sumcchiapchi2a = 1.1717D0*ETAcvl3chi2a4pchi2a-4.1075D0*ETAcvl3chi2
     1   a3pchi2a+5.0783D0*ETAcvl3chi2a2pchi2a-7.472D-1*ETAcvl3chi2a1pch
     2   i2a
      sumcchiaschi2a = 1.1717D0*ETAcvl3chi2a4schi2a-4.1075D0*ETAcvl3chi2
     1   a3schi2a+5.0783D0*ETAcvl3chi2a2schi2a-7.472D-1*ETAcvl3chi2a1sch
     2   i2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      ETAcvl2chi2b3schi2b = -4.8D-2*chi2b*(2.0D-1*chi2b-1.0D0)/cvl2chi2b
     1   p15
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      ETAcvl2chi2b4schi2b = -6.4D-3*(4.0D-1*chi2b-3.0D0)*chi2b2/cvl2chi2
     1   bp16
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      ETAcvl3chi2b3schi2b = -3.84D-7*chi2b*(4.0D-3*chi2b-1.0D0)/cvl3chi2
     1   bp15
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
      ETAcvl3chi2b4schi2b = -1.024D-9*(8.0D-3*chi2b-3.0D0)*chi2b2/cvl3ch
     1   i2bp16
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 2.4853D0*ETAcvl2chi2b4-1.9925D0*ETAcvl2chi2b3+4.329D-1*
     1   ETAcvl2chi2b2-2.607D-1*ETAcvl2chi2b1+4.8951D-1
      sumbchibpchi2b = 2.4853D0*ETAcvl2chi2b4pchi2b-1.9925D0*ETAcvl2chi2
     1   b3pchi2b+4.329D-1*ETAcvl2chi2b2pchi2b-2.607D-1*ETAcvl2chi2b1pch
     2   i2b
      sumbchibschi2b = 2.4853D0*ETAcvl2chi2b4schi2b-1.9925D0*ETAcvl2chi2
     1   b3schi2b+4.329D-1*ETAcvl2chi2b2schi2b-2.607D-1*ETAcvl2chi2b1sch
     2   i2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 1.1717D0*ETAcvl3chi2b4-4.1075D0*ETAcvl3chi2b3+5.0783D0*
     1   ETAcvl3chi2b2-7.472D-1*ETAcvl3chi2b1+1.09163D0
      sumcchibpchi2b = 1.1717D0*ETAcvl3chi2b4pchi2b-4.1075D0*ETAcvl3chi2
     1   b3pchi2b+5.0783D0*ETAcvl3chi2b2pchi2b-7.472D-1*ETAcvl3chi2b1pch
     2   i2b
      sumcchibschi2b = 1.1717D0*ETAcvl3chi2b4schi2b-4.1075D0*ETAcvl3chi2
     1   b3schi2b+5.0783D0*ETAcvl3chi2b2schi2b-7.472D-1*ETAcvl3chi2b1sch
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      ETAcvl1d23sd2 = -1.296D-6*d2*(6.0D-3*d2-1.0D0)/cvl1d2p15
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
      ETAcvl1d24sd2 = -5.184D-9*(1.2D-2*d2-3.0D0)*d22/cvl1d2p16
c...
c...  sumad and its derivatives
c...
      sumad = -1.1323D1*ETAcvl1d24+2.311D1*ETAcvl1d23-2.4707D1*ETAcvl1d2
     1   2+6.9298D0*ETAcvl1d21+5.1473D-1
      sumadpd2 = -1.1323D1*ETAcvl1d24pd2+2.311D1*ETAcvl1d23pd2-2.4707D1*
     1   ETAcvl1d22pd2+6.9298D0*ETAcvl1d21pd2
      sumadsd2 = -1.1323D1*ETAcvl1d24sd2+2.311D1*ETAcvl1d23sd2-2.4707D1*
     1   ETAcvl1d22sd2+6.9298D0*ETAcvl1d21sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-1.0D0*eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth147dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a,ETAcvl1chi2a3schi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a,ETAcvl1chi2a4schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      ETAcvl1chi2a3schi2a = -1.296D-6*chi2a*(6.0D-3*chi2a-1.0D0)/cvl1chi
     1   2ap15
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      ETAcvl1chi2a4schi2a = -5.184D-9*(1.2D-2*chi2a-3.0D0)*chi2a2/cvl1ch
     1   i2ap16
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumachia and its derivatives
c...
      sumachia = -2.0428D1*ETAcvl1chi2a4+3.5033D1*ETAcvl1chi2a3-2.8382D1
     1   *ETAcvl1chi2a2+7.0146D0*ETAcvl1chi2a1+5.4235D-1
      sumachiapchi2a = -2.0428D1*ETAcvl1chi2a4pchi2a+3.5033D1*ETAcvl1chi
     1   2a3pchi2a-2.8382D1*ETAcvl1chi2a2pchi2a+7.0146D0*ETAcvl1chi2a1pc
     2   hi2a
      sumachiaschi2a = -2.0428D1*ETAcvl1chi2a4schi2a+3.5033D1*ETAcvl1chi
     1   2a3schi2a-2.8382D1*ETAcvl1chi2a2schi2a+7.0146D0*ETAcvl1chi2a1sc
     2   hi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 8.854D-1*ETAcvl2chi2a4+1.0575D0*ETAcvl2chi2a3-1.3064D0*
     1   ETAcvl2chi2a2-1.71D-2*ETAcvl2chi2a1+5.6258D-1
      sumbchiapchi2a = 8.854D-1*ETAcvl2chi2a4pchi2a+1.0575D0*ETAcvl2chi2
     1   a3pchi2a-1.3064D0*ETAcvl2chi2a2pchi2a-1.71D-2*ETAcvl2chi2a1pchi
     2   2a
      sumbchiaschi2a = 8.854D-1*ETAcvl2chi2a4schi2a+1.0575D0*ETAcvl2chi2
     1   a3schi2a-1.3064D0*ETAcvl2chi2a2schi2a-1.71D-2*ETAcvl2chi2a1schi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 3.0454D0*ETAcvl3chi2a4-5.8676D0*ETAcvl3chi2a3+5.5721D0*
     1   ETAcvl3chi2a2-7.992D-1*ETAcvl3chi2a1+1.09025D0
      sumcchiapchi2a = 3.0454D0*ETAcvl3chi2a4pchi2a-5.8676D0*ETAcvl3chi2
     1   a3pchi2a+5.5721D0*ETAcvl3chi2a2pchi2a-7.992D-1*ETAcvl3chi2a1pch
     2   i2a
      sumcchiaschi2a = 3.0454D0*ETAcvl3chi2a4schi2a-5.8676D0*ETAcvl3chi2
     1   a3schi2a+5.5721D0*ETAcvl3chi2a2schi2a-7.992D-1*ETAcvl3chi2a1sch
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
     4   )+6.974684109057759D0*(eaapra-ea0pra)*ra113*sumachia)/ra113
      dsaa = -2.150635034570249D-1*ram83*(4.326748710922225D0*ra43*sumcc
     1   hiapchi2a-4.649789406038506D0*ea0*sumbchiapchi2a-4.649789406038
     2   506D0*(eaa-ea0)*sumachiapchi2a)
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-1.0D0*ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth147du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
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
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2,ETAcvl1d23sd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2,ETAcvl1d24sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b,ETAcvl2chi2b3schi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b,ETAcvl2chi2b4schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b,ETAcvl3chi2b3schi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b,ETAcvl3chi2b4schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 8.854D-1*ETAcvl2chi2a4+1.0575D0*ETAcvl2chi2a3-1.3064D0*
     1   ETAcvl2chi2a2-1.71D-2*ETAcvl2chi2a1+5.6258D-1
      sumbchiapchi2a = 8.854D-1*ETAcvl2chi2a4pchi2a+1.0575D0*ETAcvl2chi2
     1   a3pchi2a-1.3064D0*ETAcvl2chi2a2pchi2a-1.71D-2*ETAcvl2chi2a1pchi
     2   2a
      sumbchiaschi2a = 8.854D-1*ETAcvl2chi2a4schi2a+1.0575D0*ETAcvl2chi2
     1   a3schi2a-1.3064D0*ETAcvl2chi2a2schi2a-1.71D-2*ETAcvl2chi2a1schi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 3.0454D0*ETAcvl3chi2a4-5.8676D0*ETAcvl3chi2a3+5.5721D0*
     1   ETAcvl3chi2a2-7.992D-1*ETAcvl3chi2a1+1.09025D0
      sumcchiapchi2a = 3.0454D0*ETAcvl3chi2a4pchi2a-5.8676D0*ETAcvl3chi2
     1   a3pchi2a+5.5721D0*ETAcvl3chi2a2pchi2a-7.992D-1*ETAcvl3chi2a1pch
     2   i2a
      sumcchiaschi2a = 3.0454D0*ETAcvl3chi2a4schi2a-5.8676D0*ETAcvl3chi2
     1   a3schi2a+5.5721D0*ETAcvl3chi2a2schi2a-7.992D-1*ETAcvl3chi2a1sch
     2   i2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      ETAcvl2chi2b3schi2b = -4.8D-2*chi2b*(2.0D-1*chi2b-1.0D0)/cvl2chi2b
     1   p15
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      ETAcvl2chi2b4schi2b = -6.4D-3*(4.0D-1*chi2b-3.0D0)*chi2b2/cvl2chi2
     1   bp16
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      ETAcvl3chi2b3schi2b = -3.84D-7*chi2b*(4.0D-3*chi2b-1.0D0)/cvl3chi2
     1   bp15
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
      ETAcvl3chi2b4schi2b = -1.024D-9*(8.0D-3*chi2b-3.0D0)*chi2b2/cvl3ch
     1   i2bp16
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 8.854D-1*ETAcvl2chi2b4+1.0575D0*ETAcvl2chi2b3-1.3064D0*
     1   ETAcvl2chi2b2-1.71D-2*ETAcvl2chi2b1+5.6258D-1
      sumbchibpchi2b = 8.854D-1*ETAcvl2chi2b4pchi2b+1.0575D0*ETAcvl2chi2
     1   b3pchi2b-1.3064D0*ETAcvl2chi2b2pchi2b-1.71D-2*ETAcvl2chi2b1pchi
     2   2b
      sumbchibschi2b = 8.854D-1*ETAcvl2chi2b4schi2b+1.0575D0*ETAcvl2chi2
     1   b3schi2b-1.3064D0*ETAcvl2chi2b2schi2b-1.71D-2*ETAcvl2chi2b1schi
     2   2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 3.0454D0*ETAcvl3chi2b4-5.8676D0*ETAcvl3chi2b3+5.5721D0*
     1   ETAcvl3chi2b2-7.992D-1*ETAcvl3chi2b1+1.09025D0
      sumcchibpchi2b = 3.0454D0*ETAcvl3chi2b4pchi2b-5.8676D0*ETAcvl3chi2
     1   b3pchi2b+5.5721D0*ETAcvl3chi2b2pchi2b-7.992D-1*ETAcvl3chi2b1pch
     2   i2b
      sumcchibschi2b = 3.0454D0*ETAcvl3chi2b4schi2b-5.8676D0*ETAcvl3chi2
     1   b3schi2b+5.5721D0*ETAcvl3chi2b2schi2b-7.992D-1*ETAcvl3chi2b1sch
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      ETAcvl1d23sd2 = -1.296D-6*d2*(6.0D-3*d2-1.0D0)/cvl1d2p15
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
      ETAcvl1d24sd2 = -5.184D-9*(1.2D-2*d2-3.0D0)*d22/cvl1d2p16
c...
c...  sumad and its derivatives
c...
      sumad = -2.0428D1*ETAcvl1d24+3.5033D1*ETAcvl1d23-2.8382D1*ETAcvl1d
     1   22+7.0146D0*ETAcvl1d21+5.4235D-1
      sumadpd2 = -2.0428D1*ETAcvl1d24pd2+3.5033D1*ETAcvl1d23pd2-2.8382D1
     1   *ETAcvl1d22pd2+7.0146D0*ETAcvl1d21pd2
      sumadsd2 = -2.0428D1*ETAcvl1d24sd2+3.5033D1*ETAcvl1d23sd2-2.8382D1
     1   *ETAcvl1d22sd2+7.0146D0*ETAcvl1d21sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end
c=======================================================================
      subroutine hcth407dr(ra,saa,ec,dra,dsaa,drara,
     $                    drasaa,dsaasaa)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthr.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 407
c...
c...               A. D. Boese and  N. C. Handy
c...              J. Chem. Phys. 114, 5497 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Closed Shell
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...
      real*8 ra,saa,ec,dra,dsaa
      real*8 drara,drasaa,dsaasaa
      real*8 eablc
      real*8 eaa,eaapra,eaasra,eaasrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 saa2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 cvl1chi2ap1,cvl1chi2ap12,cvl1chi2ap13,cvl1chi2ap14
      real*8 cvl1chi2ap15,cvl1chi2ap16
      real*8 ETAcvl1chi2a1,ETAcvl1chi2a1pchi2a,ETAcvl1chi2a1schi2a
      real*8 ETAcvl1chi2a2,ETAcvl1chi2a2pchi2a,ETAcvl1chi2a2schi2a
      real*8 ETAcvl1chi2a3,ETAcvl1chi2a3pchi2a,ETAcvl1chi2a3schi2a
      real*8 ETAcvl1chi2a4,ETAcvl1chi2a4pchi2a,ETAcvl1chi2a4schi2a
      real*8 sumachia,sumachiapchi2a,sumachiaschi2a
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 zero,one,third,tt
      parameter (zero=0.0d0,one=1.0d0,Third=1.0d0/3.0d0,tt=2.0d0/3.0d0)
c...
c...  First, let's get the local terms out of the way.
c...
      r = 2.0D0*ra
      call pw91lcdu(ra,ra,eablc,eaapra,junk,eaasra,
     $              junk,eaasrarb)
      eaa = eablc*ra
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl1chi2a1schi2a = -7.2D-5/cvl1chi2ap13
      cvl1chi2ap14 = cvl1chi2ap1*cvl1chi2ap13
      ETAcvl1chi2a2 = 3.6D-5*chi2a2/cvl1chi2ap12
      ETAcvl1chi2a2pchi2a = 7.2D-5*chi2a/cvl1chi2ap13
      ETAcvl1chi2a2schi2a = -7.2D-5*(1.2D-2*chi2a-1.0D0)/cvl1chi2ap14
      cvl1chi2ap15 = cvl1chi2ap1*cvl1chi2ap14
      ETAcvl1chi2a3 = 2.16D-7*chi2a3/cvl1chi2ap13
      ETAcvl1chi2a3pchi2a = 6.48D-7*chi2a2/cvl1chi2ap14
      ETAcvl1chi2a3schi2a = -1.296D-6*chi2a*(6.0D-3*chi2a-1.0D0)/cvl1chi
     1   2ap15
      cvl1chi2ap16 = cvl1chi2ap1*cvl1chi2ap15
      ETAcvl1chi2a4 = 1.296D-9*chi2a4/cvl1chi2ap14
      ETAcvl1chi2a4pchi2a = 5.184D-9*chi2a3/cvl1chi2ap15
      ETAcvl1chi2a4schi2a = -5.184D-9*(1.2D-2*chi2a-3.0D0)*chi2a2/cvl1ch
     1   i2ap16
      cvl2chi2ap1 = 2.0D-1*chi2a+1.0D0
      cvl2chi2ap12 = cvl2chi2ap1**2
      cvl2chi2ap13 = cvl2chi2ap1*cvl2chi2ap12
      ETAcvl2chi2a1 = 2.0D-1*chi2a/cvl2chi2ap1
      ETAcvl2chi2a1pchi2a = 2.0D-1/cvl2chi2ap12
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumachia and its derivatives
c...
      sumachia = -4.2005D1*ETAcvl1chi2a4+4.2572D1*ETAcvl1chi2a3-1.9222D1
     1   *ETAcvl1chi2a2+4.4237D0*ETAcvl1chi2a1+5.8908D-1
      sumachiapchi2a = -4.2005D1*ETAcvl1chi2a4pchi2a+4.2572D1*ETAcvl1chi
     1   2a3pchi2a-1.9222D1*ETAcvl1chi2a2pchi2a+4.4237D0*ETAcvl1chi2a1pc
     2   hi2a
      sumachiaschi2a = -4.2005D1*ETAcvl1chi2a4schi2a+4.2572D1*ETAcvl1chi
     1   2a3schi2a-1.9222D1*ETAcvl1chi2a2schi2a+4.4237D0*ETAcvl1chi2a1sc
     2   hi2a
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 6.248D0*ETAcvl2chi2a4-9.1792D0*ETAcvl2chi2a3+5.6174D0*E
     1   TAcvl2chi2a2-2.4029D0*ETAcvl2chi2a1+1.18777D0
      sumbchiapchi2a = 6.248D0*ETAcvl2chi2a4pchi2a-9.1792D0*ETAcvl2chi2a
     1   3pchi2a+5.6174D0*ETAcvl2chi2a2pchi2a-2.4029D0*ETAcvl2chi2a1pchi
     2   2a
      sumbchiaschi2a = 6.248D0*ETAcvl2chi2a4schi2a-9.1792D0*ETAcvl2chi2a
     1   3schi2a+5.6174D0*ETAcvl2chi2a2schi2a-2.4029D0*ETAcvl2chi2a1schi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 2.2886D0*ETAcvl3chi2a4-2.629D0*ETAcvl3chi2a3+3.4256D0*E
     1   TAcvl3chi2a2-5.183D-1*ETAcvl3chi2a1+1.08184D0
      sumcchiapchi2a = 2.2886D0*ETAcvl3chi2a4pchi2a-2.629D0*ETAcvl3chi2a
     1   3pchi2a+3.4256D0*ETAcvl3chi2a2pchi2a-5.183D-1*ETAcvl3chi2a1pchi
     2   2a
      sumcchiaschi2a = 2.2886D0*ETAcvl3chi2a4schi2a-2.629D0*ETAcvl3chi2a
     1   3schi2a+3.4256D0*ETAcvl3chi2a2schi2a-5.183D-1*ETAcvl3chi2a1schi
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+1.11594945
     3   7449241D2*(eaapra-ea0pra)*sumachiapchi2a)+ra7*(3.0D0*(2.8
     4   84499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbc
     5   hia)-2.092405232717328D1*(eaasrarb+eaasra-ea0sra)*ra23*su
     6   machia)+1.464591887561523D0*ra13*saa2*(-1.015936673259648D2*ea0
     7   *sumbchiaschi2a-1.015936673259648D2*(eaa-ea0)*sumachiasch
     8   i2a)+1.464591887561523D0*ra3*saa*(-1.396912925732016D2*ea0*sumb
     9   chiapchi2a-1.396912925732016D2*(eaa-ea0)*sumachiapchi2a))
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(3.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eaapra-ea0pra)*sumachiapchi2a)+1.464591887561523D0*r
     4   a13*saa*(-1.26992084157456D1*ea0*sumbchiaschi2a-1.2699208415745
     5   6D1*(eaa-ea0)*sumachiaschi2a)+1.464591887561523D0*ra3*(-1
     6   .26992084157456D1*ea0*sumbchiapchi2a-1.26992084157456D1*(eaa-
     7   ea0)*sumachiapchi2a))
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a-4.64978940
     2   6038506D0*(eaa-ea0)*sumachiaschi2a)
c...
      return
      end
c=======================================================================
      subroutine hcth407du(ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb,drara,
     $                    drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb,
     $                    dsaasaa,dsbbsbb,dsaasbb)
      implicit none !do not comment this line, no implicit convention here
c...
c...  MM (05/14/2003) generated with maxima, file hcthu.mc
c...
c...  Hamprecht, Cohen, Tozer and Handy (HCTH) functional,
c...                    parametrization 407
c...
c...               A. D. Boese and  N. C. Handy
c...              J. Chem. Phys. 114, 5497 (2001)
c...
c...  Formulation taken from the Molpro manual
c...  http://www.molpro.net/current/molpro_manual
c...
c...
c...   *** Open Shell
c...
c...  This subroutine computes the functional and its first
c...  and second derivatives.
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
c...  drara   funct. 2nd deriv. w.r.t. alpha density
c...  drbrb   funct. 2nd deriv. w.r.t. beta density
c...  drarb   funct. mixed 2nd deriv. w.r.t. alpha and beta density
c...  drasaa  funct. mixed 2nd deriv. w.r.t. alpha den. and alpha grad.
c...  drasbb  funct. mixed 2nd deriv. w.r.t. alpha den. and beta grad.
c...  drbsaa  funct. mixed 2nd deriv. w.r.t. beta den. and alpha grad.
c...  drbsbb  funct. mixed 2nd deriv. w.r.t. beta den. and beta grad.
c...  dsaasaa funct. 2nd deriv. w.r.t. alpha grad.
c...  dsbbsbb funct. 2nd deriv. w.r.t. beta grad.
c...  dsabsab funct. 2nd deriv. w.r.t. alpha beta grad.
c...
      real*8 ra,rb,saa,sbb,ec,dra,drb,dsaa,dsbb
      real*8 drara,drbrb,drarb,drasaa,drasbb,drbsaa,drbsbb
      real*8 dsaasaa,dsbbsbb,dsaasbb
      real*8 eablc,deablcdra,deablcdrb,deablcdrara,deablcdrbrb,
     $       deablcdrarb
      real*8 eab,eabpra,eabprb,eabsra,eabsrb,eabsrarb
      real*8 ea0lc,dea0lcdra,dea0lcdrara
      real*8 ea0,ea0pra,ea0sra
      real*8 eb0lc,deb0lcdrb,deb0lcdrbrb
      real*8 eb0,eb0prb,eb0srb
      real*8 junk
      real*8 r
      real*8 ra2,ra3,ra4,ra7,ra13,ra23,ra43,ra53,ra83,ram83,ra113,
     $       ram113,ram163,ram203,ram233
      real*8 rb2,rb3,rb4,rb7,rb13,rb23,rb43,rb53,rb83,rbm83,rb113,
     $       rbm113,rbm163,rbm203,rbm233
      real*8 saa2,sbb2
      real*8 chi2a,chi2a2,chi2a3,chi2a4
      real*8 chi2b,chi2b2,chi2b3,chi2b4
      real*8 d2,d22,d23,d24
      real*8 cvl1d2p1,cvl1d2p12,cvl1d2p13,cvl1d2p14
      real*8 cvl1d2p15,cvl1d2p16
      real*8 ETAcvl1d21,ETAcvl1d21pd2,ETAcvl1d21sd2
      real*8 ETAcvl1d22,ETAcvl1d22pd2,ETAcvl1d22sd2
      real*8 ETAcvl1d23,ETAcvl1d23pd2,ETAcvl1d23sd2
      real*8 ETAcvl1d24,ETAcvl1d24pd2,ETAcvl1d24sd2
      real*8 sumad,sumadpd2,sumadsd2
      real*8 cvl2chi2ap1,cvl2chi2ap12,cvl2chi2ap13,cvl2chi2ap14
      real*8 cvl2chi2ap15,cvl2chi2ap16
      real*8 ETAcvl2chi2a1,ETAcvl2chi2a1pchi2a,ETAcvl2chi2a1schi2a
      real*8 ETAcvl2chi2a2,ETAcvl2chi2a2pchi2a,ETAcvl2chi2a2schi2a
      real*8 ETAcvl2chi2a3,ETAcvl2chi2a3pchi2a,ETAcvl2chi2a3schi2a
      real*8 ETAcvl2chi2a4,ETAcvl2chi2a4pchi2a,ETAcvl2chi2a4schi2a
      real*8 cvl3chi2ap1,cvl3chi2ap12,cvl3chi2ap13,cvl3chi2ap14
      real*8 cvl3chi2ap15,cvl3chi2ap16
      real*8 ETAcvl3chi2a1,ETAcvl3chi2a1pchi2a,ETAcvl3chi2a1schi2a
      real*8 ETAcvl3chi2a2,ETAcvl3chi2a2pchi2a,ETAcvl3chi2a2schi2a
      real*8 ETAcvl3chi2a3,ETAcvl3chi2a3pchi2a,ETAcvl3chi2a3schi2a
      real*8 ETAcvl3chi2a4,ETAcvl3chi2a4pchi2a,ETAcvl3chi2a4schi2a
      real*8 sumbchia,sumbchiapchi2a,sumbchiaschi2a
      real*8 sumcchia,sumcchiapchi2a,sumcchiaschi2a
      real*8 cvl2chi2bp1,cvl2chi2bp12,cvl2chi2bp13,cvl2chi2bp14
      real*8 cvl2chi2bp15,cvl2chi2bp16
      real*8 ETAcvl2chi2b1,ETAcvl2chi2b1pchi2b,ETAcvl2chi2b1schi2b
      real*8 ETAcvl2chi2b2,ETAcvl2chi2b2pchi2b,ETAcvl2chi2b2schi2b
      real*8 ETAcvl2chi2b3,ETAcvl2chi2b3pchi2b,ETAcvl2chi2b3schi2b
      real*8 ETAcvl2chi2b4,ETAcvl2chi2b4pchi2b,ETAcvl2chi2b4schi2b
      real*8 cvl3chi2bp1,cvl3chi2bp12,cvl3chi2bp13,cvl3chi2bp14
      real*8 cvl3chi2bp15,cvl3chi2bp16
      real*8 ETAcvl3chi2b1,ETAcvl3chi2b1pchi2b,ETAcvl3chi2b1schi2b
      real*8 ETAcvl3chi2b2,ETAcvl3chi2b2pchi2b,ETAcvl3chi2b2schi2b
      real*8 ETAcvl3chi2b3,ETAcvl3chi2b3pchi2b,ETAcvl3chi2b3schi2b
      real*8 ETAcvl3chi2b4,ETAcvl3chi2b4pchi2b,ETAcvl3chi2b4schi2b
      real*8 sumbchib,sumbchibpchi2b,sumbchibschi2b
      real*8 sumcchib,sumcchibpchi2b,sumcchibschi2b
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
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
c...
c...  eb0 only
c...
      else if(ra.lt.epsi)then
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
c...
c...  eab, ea0, and eb0
c...
      else
      call pw91lcdu(ra,rb,eablc,eabpra,eabprb,eabsra,
     $              eabsrb,eabsrarb)
      eab = eablc*r
      call pw91lcdu(ra,zero,ea0lc,ea0pra,junk,ea0sra,
     $              junk,junk)
      ea0 = ea0lc*ra
      call pw91lcdu(rb,zero,eb0lc,eb0prb,junk,eb0srb,
     $              junk,junk)
      eb0 = eb0lc*rb
      endif
c...
c...  Nonlocal part
c...
      saa2 = saa**2
      sbb2 = sbb**2
c...
c...  Tests ra, and computes the related quantities
c...
      if(ra.gt.epsi)then
      ra2 = ra**2
      ra3 = ra*ra2
      ra4 = ra*ra3
      ra7 = ra3*ra4
      ra13 = ra**THIRD
      ra23 = ra13**2
      ra43 = ra23**2
      ra53 = ra13*ra43
      ra83 = ra*ra53
      ram83 = one/ra83
      ra113 = ra*ra83
      ram113 = ram83/ra
      ram163 = ram113/ra53
      ram203 = ram163/ra43
      ram233 = ram203/ra
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
      ETAcvl2chi2a1schi2a = -8.0D-2/cvl2chi2ap13
      cvl2chi2ap14 = cvl2chi2ap1*cvl2chi2ap13
      ETAcvl2chi2a2 = 4.0D-2*chi2a2/cvl2chi2ap12
      ETAcvl2chi2a2pchi2a = 8.0D-2*chi2a/cvl2chi2ap13
      ETAcvl2chi2a2schi2a = -8.0D-2*(4.0D-1*chi2a-1.0D0)/cvl2chi2ap14
      cvl2chi2ap15 = cvl2chi2ap1*cvl2chi2ap14
      ETAcvl2chi2a3 = 8.0D-3*chi2a3/cvl2chi2ap13
      ETAcvl2chi2a3pchi2a = 2.4D-2*chi2a2/cvl2chi2ap14
      ETAcvl2chi2a3schi2a = -4.8D-2*chi2a*(2.0D-1*chi2a-1.0D0)/cvl2chi2a
     1   p15
      cvl2chi2ap16 = cvl2chi2ap1*cvl2chi2ap15
      ETAcvl2chi2a4 = 1.6D-3*chi2a4/cvl2chi2ap14
      ETAcvl2chi2a4pchi2a = 6.4D-3*chi2a3/cvl2chi2ap15
      ETAcvl2chi2a4schi2a = -6.4D-3*(4.0D-1*chi2a-3.0D0)*chi2a2/cvl2chi2
     1   ap16
      cvl3chi2ap1 = 4.0D-3*chi2a+1.0D0
      cvl3chi2ap12 = cvl3chi2ap1**2
      cvl3chi2ap13 = cvl3chi2ap1*cvl3chi2ap12
      ETAcvl3chi2a1 = 4.0D-3*chi2a/cvl3chi2ap1
      ETAcvl3chi2a1pchi2a = 4.0D-3/cvl3chi2ap12
      ETAcvl3chi2a1schi2a = -3.2D-5/cvl3chi2ap13
      cvl3chi2ap14 = cvl3chi2ap1*cvl3chi2ap13
      ETAcvl3chi2a2 = 1.6D-5*chi2a2/cvl3chi2ap12
      ETAcvl3chi2a2pchi2a = 3.2D-5*chi2a/cvl3chi2ap13
      ETAcvl3chi2a2schi2a = -3.2D-5*(8.0D-3*chi2a-1.0D0)/cvl3chi2ap14
      cvl3chi2ap15 = cvl3chi2ap1*cvl3chi2ap14
      ETAcvl3chi2a3 = 6.4D-8*chi2a3/cvl3chi2ap13
      ETAcvl3chi2a3pchi2a = 1.92D-7*chi2a2/cvl3chi2ap14
      ETAcvl3chi2a3schi2a = -3.84D-7*chi2a*(4.0D-3*chi2a-1.0D0)/cvl3chi2
     1   ap15
      cvl3chi2ap16 = cvl3chi2ap1*cvl3chi2ap15
      ETAcvl3chi2a4 = 2.56D-10*chi2a4/cvl3chi2ap14
      ETAcvl3chi2a4pchi2a = 1.024D-9*chi2a3/cvl3chi2ap15
      ETAcvl3chi2a4schi2a = -1.024D-9*(8.0D-3*chi2a-3.0D0)*chi2a2/cvl3ch
     1   i2ap16
c...
c...  sumbchia and its derivatives
c...
      sumbchia = 6.248D0*ETAcvl2chi2a4-9.1792D0*ETAcvl2chi2a3+5.6174D0*E
     1   TAcvl2chi2a2-2.4029D0*ETAcvl2chi2a1+1.18777D0
      sumbchiapchi2a = 6.248D0*ETAcvl2chi2a4pchi2a-9.1792D0*ETAcvl2chi2a
     1   3pchi2a+5.6174D0*ETAcvl2chi2a2pchi2a-2.4029D0*ETAcvl2chi2a1pchi
     2   2a
      sumbchiaschi2a = 6.248D0*ETAcvl2chi2a4schi2a-9.1792D0*ETAcvl2chi2a
     1   3schi2a+5.6174D0*ETAcvl2chi2a2schi2a-2.4029D0*ETAcvl2chi2a1schi
     2   2a
c...
c...  sumcchia and its derivatives
c...
      sumcchia = 2.2886D0*ETAcvl3chi2a4-2.629D0*ETAcvl3chi2a3+3.4256D0*E
     1   TAcvl3chi2a2-5.183D-1*ETAcvl3chi2a1+1.08184D0
      sumcchiapchi2a = 2.2886D0*ETAcvl3chi2a4pchi2a-2.629D0*ETAcvl3chi2a
     1   3pchi2a+3.4256D0*ETAcvl3chi2a2pchi2a-5.183D-1*ETAcvl3chi2a1pchi
     2   2a
      sumcchiaschi2a = 2.2886D0*ETAcvl3chi2a4schi2a-2.629D0*ETAcvl3chi2a
     1   3schi2a+3.4256D0*ETAcvl3chi2a2schi2a-5.183D-1*ETAcvl3chi2a1schi
     2   2a
      endif
c...
c...  Tests rb, and computes the related quantities
c...
      if(rb.gt.epsi)then
      rb2 = rb**2
      rb3 = rb*rb2
      rb4 = rb*rb3
      rb7 = rb3*rb4
      rb13 = rb**THIRD
      rb23 = rb13**2
      rb43 = rb23**2
      rb53 = rb13*rb43
      rb83 = rb*rb53
      rbm83 = one/rb83
      rb113 = rb*rb83
      rbm113 = rbm83/rb
      rbm163 = rbm113/rb53
      rbm203 = rbm163/rb43
      rbm233 = rbm203/rb
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
      ETAcvl2chi2b1schi2b = -8.0D-2/cvl2chi2bp13
      cvl2chi2bp14 = cvl2chi2bp1*cvl2chi2bp13
      ETAcvl2chi2b2 = 4.0D-2*chi2b2/cvl2chi2bp12
      ETAcvl2chi2b2pchi2b = 8.0D-2*chi2b/cvl2chi2bp13
      ETAcvl2chi2b2schi2b = -8.0D-2*(4.0D-1*chi2b-1.0D0)/cvl2chi2bp14
      cvl2chi2bp15 = cvl2chi2bp1*cvl2chi2bp14
      ETAcvl2chi2b3 = 8.0D-3*chi2b3/cvl2chi2bp13
      ETAcvl2chi2b3pchi2b = 2.4D-2*chi2b2/cvl2chi2bp14
      ETAcvl2chi2b3schi2b = -4.8D-2*chi2b*(2.0D-1*chi2b-1.0D0)/cvl2chi2b
     1   p15
      cvl2chi2bp16 = cvl2chi2bp1*cvl2chi2bp15
      ETAcvl2chi2b4 = 1.6D-3*chi2b4/cvl2chi2bp14
      ETAcvl2chi2b4pchi2b = 6.4D-3*chi2b3/cvl2chi2bp15
      ETAcvl2chi2b4schi2b = -6.4D-3*(4.0D-1*chi2b-3.0D0)*chi2b2/cvl2chi2
     1   bp16
      cvl3chi2bp1 = 4.0D-3*chi2b+1.0D0
      cvl3chi2bp12 = cvl3chi2bp1**2
      cvl3chi2bp13 = cvl3chi2bp1*cvl3chi2bp12
      ETAcvl3chi2b1 = 4.0D-3*chi2b/cvl3chi2bp1
      ETAcvl3chi2b1pchi2b = 4.0D-3/cvl3chi2bp12
      ETAcvl3chi2b1schi2b = -3.2D-5/cvl3chi2bp13
      cvl3chi2bp14 = cvl3chi2bp1*cvl3chi2bp13
      ETAcvl3chi2b2 = 1.6D-5*chi2b2/cvl3chi2bp12
      ETAcvl3chi2b2pchi2b = 3.2D-5*chi2b/cvl3chi2bp13
      ETAcvl3chi2b2schi2b = -3.2D-5*(8.0D-3*chi2b-1.0D0)/cvl3chi2bp14
      cvl3chi2bp15 = cvl3chi2bp1*cvl3chi2bp14
      ETAcvl3chi2b3 = 6.4D-8*chi2b3/cvl3chi2bp13
      ETAcvl3chi2b3pchi2b = 1.92D-7*chi2b2/cvl3chi2bp14
      ETAcvl3chi2b3schi2b = -3.84D-7*chi2b*(4.0D-3*chi2b-1.0D0)/cvl3chi2
     1   bp15
      cvl3chi2bp16 = cvl3chi2bp1*cvl3chi2bp15
      ETAcvl3chi2b4 = 2.56D-10*chi2b4/cvl3chi2bp14
      ETAcvl3chi2b4pchi2b = 1.024D-9*chi2b3/cvl3chi2bp15
      ETAcvl3chi2b4schi2b = -1.024D-9*(8.0D-3*chi2b-3.0D0)*chi2b2/cvl3ch
     1   i2bp16
c...
c...  sumbchib and its derivatives
c...
      sumbchib = 6.248D0*ETAcvl2chi2b4-9.1792D0*ETAcvl2chi2b3+5.6174D0*E
     1   TAcvl2chi2b2-2.4029D0*ETAcvl2chi2b1+1.18777D0
      sumbchibpchi2b = 6.248D0*ETAcvl2chi2b4pchi2b-9.1792D0*ETAcvl2chi2b
     1   3pchi2b+5.6174D0*ETAcvl2chi2b2pchi2b-2.4029D0*ETAcvl2chi2b1pchi
     2   2b
      sumbchibschi2b = 6.248D0*ETAcvl2chi2b4schi2b-9.1792D0*ETAcvl2chi2b
     1   3schi2b+5.6174D0*ETAcvl2chi2b2schi2b-2.4029D0*ETAcvl2chi2b1schi
     2   2b
c...
c...  sumcchib and its derivatives
c...
      sumcchib = 2.2886D0*ETAcvl3chi2b4-2.629D0*ETAcvl3chi2b3+3.4256D0*E
     1   TAcvl3chi2b2-5.183D-1*ETAcvl3chi2b1+1.08184D0
      sumcchibpchi2b = 2.2886D0*ETAcvl3chi2b4pchi2b-2.629D0*ETAcvl3chi2b
     1   3pchi2b+3.4256D0*ETAcvl3chi2b2pchi2b-5.183D-1*ETAcvl3chi2b1pchi
     2   2b
      sumcchibschi2b = 2.2886D0*ETAcvl3chi2b4schi2b-2.629D0*ETAcvl3chi2b
     1   3schi2b+3.4256D0*ETAcvl3chi2b2schi2b-5.183D-1*ETAcvl3chi2b1schi
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
      ETAcvl1d21sd2 = -7.2D-5/cvl1d2p13
      cvl1d2p14 = cvl1d2p1*cvl1d2p13
      ETAcvl1d22 = 3.6D-5*d22/cvl1d2p12
      ETAcvl1d22pd2 = 7.2D-5*d2/cvl1d2p13
      ETAcvl1d22sd2 = -7.2D-5*(1.2D-2*d2-1.0D0)/cvl1d2p14
      cvl1d2p15 = cvl1d2p1*cvl1d2p14
      ETAcvl1d23 = 2.16D-7*d23/cvl1d2p13
      ETAcvl1d23pd2 = 6.48D-7*d22/cvl1d2p14
      ETAcvl1d23sd2 = -1.296D-6*d2*(6.0D-3*d2-1.0D0)/cvl1d2p15
      cvl1d2p16 = cvl1d2p1*cvl1d2p15
      ETAcvl1d24 = 1.296D-9*d24/cvl1d2p14
      ETAcvl1d24pd2 = 5.184D-9*d23/cvl1d2p15
      ETAcvl1d24sd2 = -5.184D-9*(1.2D-2*d2-3.0D0)*d22/cvl1d2p16
c...
c...  sumad and its derivatives
c...
      sumad = -4.2005D1*ETAcvl1d24+4.2572D1*ETAcvl1d23-1.9222D1*ETAcvl1d
     1   22+4.4237D0*ETAcvl1d21+5.8908D-1
      sumadpd2 = -4.2005D1*ETAcvl1d24pd2+4.2572D1*ETAcvl1d23pd2-1.9222D1
     1   *ETAcvl1d22pd2+4.4237D0*ETAcvl1d21pd2
      sumadsd2 = -4.2005D1*ETAcvl1d24sd2+4.2572D1*ETAcvl1d23sd2-1.9222D1
     1   *ETAcvl1d22sd2+4.4237D0*ETAcvl1d21sd2
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
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+1.2D1*ra4*saa*(4.326748710922225D0*ra13*sumcch
     2   iapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+3.0D0*ra7*(
     3   2.884499140614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*su
     4   mbchia)-1.487932609932322D2*ea0*ra13*saa2*sumbchiaschi2a-2.0459
     5   07338656943D2*ea0*ra3*saa*sumbchiapchi2a)
      drbrb = ZERO
      drarb = ZERO
      drasaa = 1.433756689713499D-1*ram203*(1.73069948436889D1*ra53*saa*
     1   sumcchiaschi2a+3.0D0*ra4*(2.884499140614817D0*ra13*sumcchiapchi
     2   2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)-1.859915762415402
     3   D1*ea0*ra13*saa*sumbchiaschi2a-1.859915762415402D1*ea0*ra3*sumb
     4   chiapchi2a)
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = ZERO
      dsaasaa = -2.150635034570249D-1*ram163*(4.326748710922225D0*ra43*s
     1   umcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)
      dsbbsbb = ZERO
      dsaasbb = ZERO
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
c...  Second derivatives of the functional
c...
      drara = ZERO
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+1.2D1*rb4*sbb*(4.326748710922225D0*rb13*sumcch
     2   ibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)+3.0D0*rb7*(
     3   2.884499140614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*su
     4   mbchib)-1.487932609932322D2*eb0*rb13*sbb2*sumbchibschi2b-2.0459
     5   07338656943D2*eb0*rb3*sbb*sumbchibpchi2b)
      drarb = ZERO
      drasaa = ZERO
      drasbb = ZERO
      drbsaa = ZERO
      drbsbb = 1.433756689713499D-1*rbm203*(1.73069948436889D1*rb53*sbb*
     1   sumcchibschi2b+3.0D0*rb4*(2.884499140614817D0*rb13*sumcchibpchi
     2   2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-1.859915762415402
     3   D1*eb0*rb13*sbb*sumbchibschi2b-1.859915762415402D1*eb0*rb3*sumb
     4   chibpchi2b)
      dsaasaa = ZERO
      dsbbsbb = -2.150635034570249D-1*rbm163*(4.326748710922225D0*rb43*s
     1   umcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)
      dsaasbb = ZERO
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
c...
c...  Second derivatives of the functional
c...
      drara = -4.779188965711665D-2*ram233*(1.384559587495112D2*ra53*saa
     1   2*sumcchiaschi2a+ra4*saa*(1.2D1*(4.326748710922225D0*ra13*sumcc
     2   hiapchi2a+9.299578812077012D0*ea0pra*sumbchiapchi2a)+5.57974728
     3   7246207D1*(eabpra-ea0pra)*sumadpd2)+ra7*(3.0D0*(2.8844991
     4   40614817D0*sumcchia-6.974684109057759D0*ea0sra*ra23*sumbchia)-2
     5   .092405232717328D1*(eabsra-ea0sra)*ra23*sumad)+1.46459188
     6   7561523D0*ra13*saa2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*ea0*sumbchiaschi2a)+1.46459188756152
     8   3D0*ra3*saa*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*ea0*sumbchiapchi2a))
      drbrb = -4.779188965711665D-2*rbm233*(1.384559587495112D2*rb53*sbb
     1   2*sumcchibschi2b+rb4*sbb*(1.2D1*(4.326748710922225D0*rb13*sumcc
     2   hibpchi2b+9.299578812077012D0*eb0prb*sumbchibpchi2b)-5.57974728
     3   7246207D1*(eb0prb-eabprb)*sumadpd2)+rb7*(3.0D0*(2.8844991
     4   40614817D0*sumcchib-6.974684109057759D0*eb0srb*rb23*sumbchib)+2
     5   .092405232717328D1*(eb0srb-eabsrb)*rb23*sumad)+1.46459188
     6   7561523D0*rb13*sbb2*(2.539841683149119D1*(eb0-eab+ea0)*su
     7   madsd2-1.015936673259648D2*eb0*sumbchibschi2b)+1.46459188756152
     8   3D0*rb3*sbb*(6.984564628660078D1*(eb0-eab+ea0)*sumadpd2-1
     9   .396912925732016D2*eb0*sumbchibpchi2b))
      drarb = -1.111111111111111D-1*(1.6D1*(eb0-eab+ea0)*saa*sbb*s
     1   umadsd2+1.2D1*(eabpra-1.0D0*ea0pra)*ra113*sbb*sumadpd2-1.2D1*(e
     2   b0prb-1.0D0*eabprb)*rb113*saa*sumadpd2-9.0D0*eabsrarb*ra113*rb1
     3   13*sumad)/(ra113*rb113)
      drasaa = 7.168783448567497D-2*ram203*(3.46139896873778D1*ra53*saa*
     1   sumcchiaschi2a+ra4*(6.0D0*(2.884499140614817D0*ra13*sumcchiapch
     2   i2a+2.324894703019253D0*ea0pra*sumbchiapchi2a)+6.97468410905775
     3   9D0*(eabpra-ea0pra)*sumadpd2)+1.464591887561523D0*ra13*sa
     4   a*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*ea0*sumbchiaschi2a)+1.464591887561523D0*ra3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*ea0*
     7   sumbchiapchi2a))
      drasbb = 1.666666666666667D-1*rbm83*(4.0D0*(eb0-eab+ea0)*saa
     1   *sumadsd2+3.0D0*(eabpra-ea0pra)*ra113*sumadpd2)/ra113
      drbsaa = 1.666666666666667D-1*ram83*(4.0D0*(eb0-eab+ea0)*sbb
     1   *sumadsd2-3.0D0*(eb0prb-eabprb)*rb113*sumadpd2)/rb113
      drbsbb = 7.168783448567497D-2*rbm203*(3.46139896873778D1*rb53*sbb*
     1   sumcchibschi2b+rb4*(6.0D0*(2.884499140614817D0*rb13*sumcchibpch
     2   i2b+2.324894703019253D0*eb0prb*sumbchibpchi2b)-6.97468410905775
     3   9D0*(eb0prb-eabprb)*sumadpd2)+1.464591887561523D0*rb13*sb
     4   b*(6.349604207872798D0*(eb0-eab+ea0)*sumadsd2-2.539841683
     5   149119D1*eb0*sumbchibschi2b)+1.464591887561523D0*rb3*(1.2699208
     6   4157456D1*(eb0-eab+ea0)*sumadpd2-2.539841683149119D1*eb0*
     7   sumbchibpchi2b))
      dsaasaa = -1.075317517285125D-1*ram163*(2.0D0*(4.326748710922225D0
     1   *ra43*sumcchiaschi2a-4.649789406038506D0*ea0*sumbchiaschi2a)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsbbsbb = -1.075317517285125D-1*rbm163*(2.0D0*(4.326748710922225D0
     1   *rb43*sumcchibschi2b-4.649789406038506D0*eb0*sumbchibschi2b)+2.
     2   324894703019253D0*(eb0-eab+ea0)*sumadsd2)
      dsaasbb = -2.5D-1*(eb0-eab+ea0)*ram83*rbm83*sumadsd2
      endif
c...
      return
      end

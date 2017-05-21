c...
c...  MM March-September 2004
c...
c...  thys file contains routines specific to the calculation
c...  of the angular part of integrals over ab initio
c...  pseudopotentials.
c...  The method of calculation is based on
c...
c...  L. E. McMurchie and E. R. Davidson, J. Comp. Phys. 44, 289 (1981)
c...
c...  Altough all the code has been originally written,
c..   this work has been grately facilitated by the the availability
c...  of Fortran routines from the ARGOS code that have been used as
c...  model and for testing.
c...
c=======================================================================
      subroutine preomega(mijk,xyz)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 9/15/2004 15:14:15
c...
c...  this function computes:
c...
c...          / /
c...     1    [ [   i   j   k               (i-1)!! (j-1)!! (k-1)!!
c...   -----  I I xs  ys  zs  dPhi dTheta = -----------------------
c...   4 %pi  ] ]                                (i+j+k+1)!!
c...          / /
c...
c...  if i, j and k are all even, or zero otherwise
c...
c...
      integer maxn
      parameter (maxn=40)
c...
      integer mijk,i,j,k
      real*8 xyz(0:mijk,0:mijk,0:mijk)
c...
c...
      real*8 dbf(-1:40)
      save dbf
      data dbf/1.0D0,1.0D0,1.0D0,2.0D0,3.0D0,8.0D0,1.5D1,4.8D1,1.05D2,3.
     $  84D2,9.45D2,3.84D3,1.0395D4,4.608D4,1.35135D5,6.4512D5,2.027025D
     $  6,1.032192D7,3.4459425D7,1.8579456D8,6.54729075D8,3.7158912D9,1.
     $  3749310575D10,8.17496064D10,3.16234143225D11,1.9619905536D12,7.9
     $  05853580625D12,5.10117543936D13,2.13458046676875D14,1.4283291230
     $  208D15,6.190283353629375D15,4.2849873690624D16,1.918987839625106
     $  D17,1.371195958099968D18,6.332659870762851D18,4.662066257539891D
     $  19,2.216430954766998D20,1.678343852714361D21,8.200794532637891D2
     $  1,6.377706640314571D22,3.198309867728778D23,2.551082656125828D24
     $  /
c...
      if(mijk.ge.maxn)
     $  call nerror(1,'Preomega','maximum value exceeded',mijk,maxn)
      do i=0,mijk,2
        do j=0,mijk,2
          do k=0,mijk,2
            if(k+j+i.le.mijk)then
              xyz(i,j,k)=dbf(i-1)*dbf(j-1)*dbf(k-1)/dbf(k+j+i+1)
            endif
          enddo
        enddo
      enddo
c...
      end
c=======================================================================
      subroutine yrleg0(l,THETA,yr)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 6/28/2004 18:12:38
c...
c...  this function computes only the YR(0,l) real spherical harmonics,
c...  based on the recurrence relation of the Legendre polynomials:
c...
c...  yr(0,l,theta)=  sqrt((2*l+1)/(4*%pi)) * P(l,cos(theta))
c...
c...  where P(l,cos(theta)) is a Legendre polynomial of degree l.
c...
c...  basically, this version is to be used when the angle phi is zero.
c...
c...  The legendre polynomials obey the recurrence:
c...
c...  P(l,x)=(2*l-1)/l*x*P(l-1,x)-(l-1)/l*P(l-2,x)
c...
c...
      real*8 spi2,stsp2
      parameter (spi2=2.820947917738782D-1,stsp2=4.886025119029199D-1)
c...
      integer l,nl,IND,ind1,ind2
      real*8 THETA,x,yr(*)
c...
      x=COS(THETA)
c...
c...  First generates the Legendre polynomials
c...
      yr(1)=1.0D0
      if(l.gt.0)then
        yr(3)=x
      endif
      if(l.gt.1)then
        ind2=1
        ind1=3
        do nl=2,l,1
          IND=nl*(nl+1)+1
          yr(IND)=yr(ind1)*(2.0D0*nl-1.0D0)*x/nl-1.0D0*yr(ind2)*(nl-1.0D
     $      0)/nl
          ind2=ind1
          ind1=IND
        enddo
      endif
c...
c...  Now transforms to real spherical harmonics
c...
      yr(1)=spi2
      if(l.gt.0)then
        yr(3)=yr(3)*stsp2
      endif
      if(l.gt.1)then
        ind2=1
        ind1=3
        do nl=2,l,1
          IND=nl*(nl+1)+1
          yr(IND)=yr(IND)*SQRT(2.0D0*nl+1.0D0)*spi2
          ind2=ind1
          ind1=IND
        enddo
      endif
c...
      end
c=======================================================================
      subroutine yrleg(l,theta,phi,yr)
      implicit none
c...
c...  MM (06/09/2004) generated with maxima, file ecpcode.mc
c...
c...  This subroutine computes the real spherical harmonics up to
c...  a degree l, based on the recurrence of the associate
c...  Legendre polynomials
c...
c...  the relation between real spherical harmonics and associate
c...  Legendre polynomials is as follow:
c...
c...  yr(0,l,theta,phi)= sqrt((2*l+1)/(4*%pi))* P(0,l,cos(theta))
c...
c...  yr(|m|,l,theta,phi)= sqrt((2*l+1)/(4*%pi)*(l-|m|)!/(l+|m|)!)*
c...               P(|m|,l,cos(theta)) * 2/sqrt(2) * cos(|m|*phi)
c...
c...  yr(-|m|,l,theta,phi)=sqrt((2*l+1)/(4*%pi)*(l-|m|)!/(l+|m|)!)*
c...               P(|m|,l,cos(theta)) * 2/sqrt(2) * sin(|m|*phi)
c...
c...  where P are the associate Legendre Polynomials
c...
c...  Press, Teukolosky, Vetterling, Flannery, Numerical Recipes,
c...  Second edition, (Cambridge, 1992), page 246-248
c...
c...  Input parameters:
c...
c...  l           maximum value of quantum number l
c...  theta, phi  spherical polar angles (Radians)
c...
c...  Output parameter:
c...
c...  yr          array values of real spherical harmonics
c...
      integer l
      real*8 yr(*),theta,phi
c...
      integer lm,l2p1,ind,ind1,ind2,m,i
      real*8 x,somx2,ff,rl2m1,rllm1,rllmm,c,cc,cf,fm
c...
c...  real numerical parameters
c...
      real*8 ZERO,one,two,three,spi4m1,s3pi4m1,s22
      parameter (ZERO=0.0D0,one=1.0D0,two=2.0D0,three=3.0D0,spi4m1=2.820
     $  947917738782D-1,s3pi4m1=4.886025119029199D-1,s22=1.4142135623730
     $  95D0)
c...
c...  integer numerical parameters
c...
      integer lmax
      parameter (lmax=50)
c...
c...  local arrays
c...
      real*8 cmp(lmax),smp(lmax)
c...
      if(l.gt.lmax)
     $  call nerror(1,'yrleg','maximum l value exceeded',l,lmax)
c...
c...  First we compute the associate Legendre polynomials
c...  by recurrence, and temporarily store them in yr
c...
c...
c...  l = 0
c...
      yr(1) = one
c...
c...  l = 1
c...
      if(l.gt.0)then
        x = COS(THETA)
        somx2 = SQRT(one-x**2)
        yr(3) = yr(1)*x
        yr(4) = -yr(1)*somx2
      endif
c...
c...  l > 1
c...
      if(l.gt.1)then
        ff = three
        ind=4
        ind1=1
        do lm=2,l
          rl2m1=dfloat(2*lm-1)
          l2p1=2*lm+1
          ind2=ind1
          ind1=ind
          ind=ind+l2p1
          yr(IND) = -ff*yr(ind1)*somx2
          yr(IND-1) = yr(ind1)*rl2m1*x
          ff = two+ff
          i=0
          do m=lm-2,0,-1
            rllm1=dfloat(lm+m-1)
            rllmm=dfloat(lm-m)
            yr(IND-i-2) = (yr(ind1-i-1)*rl2m1*x-yr(ind2-i)*rllm1)/
     $        rllmm
            i=i+1
          enddo
        enddo
      endif
c...
c...  now we transform to the real spherical harmonics
c...
      yr(1) = yr(1)*spi4m1
c...
c...  l = 1
c...
      if(l.gt.0)then
        cmp(1) = COS(PHI)
        smp(1) = SIN(PHI)
        yr(3) = yr(3)*s3pi4m1
        cc = yr(4)*s3pi4m1
        yr(2) = smp(1)*cc
        yr(4) = cmp(1)*cc
      endif
c...
c...  l > 1
c...
      if(l.gt.1)then
        do m=2,l
          fm=dfloat(m)
          cmp(m) = COS(fm*PHI)
          smp(m) = SIN(fm*PHI)
        enddo
        ind=4
        do lm=2,l
          c=sqrt(dfloat(2*lm+1))*spi4m1
          ind=ind+lm+1
          yr(IND) = C*yr(IND)
          CF = one
          do m=1,lm
            cf=cf/dfloat((lm+m)*(lm-m+1))
            cc = C*SQRT(CF)*yr(m+IND)*s22
            yr(IND-m) = cc*smp(m)
            yr(m+IND) = cc*cmp(m)
          enddo
          ind=ind+lm
        enddo
      endif
c...
      return
      end
c=======================================================================
      subroutine omega0(i,j,k,yrk,mijk,xyzi,OMEGA)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 9/15/2004 15:21:40
c...
c...  This subroutine computes type 1 angular integrals over pseudopotentials
c...
c...
c...   omega(i,j,k,l)=
c...
c...                 / /
c...                 [ [                     i   j   k
c...  yr(0,l,thetak) I I yr(0,l,THETA,PHI) xs  ys  zs  dPHi dTHETA
c...                 ] ]
c...                 / /
c...
c...  where yr are real spherical harmonic functions, and xs, ys and zs
c...  are cartesian coordinates costrained over the surface of the unit
c...  sphere:
c...              xs=sin(THETA) cos(PHI)
c...              ys=sin(THETA) sin(PHI)
c...              zs=cos(THETA)
c...
c...  this version is to be used when the angle phik cannot be defined
c...  (both kx and ky are zero), and thus only the real spherical
c...  harmonics with m=0 survive
c...
c...  The integral is computed by expanding the real sperical harmonic
c...  as a polynomial in xs, ys and zs, obtaining integrals of the type:
c...
c...                      / /
c...                 1    [ [   a   b   c         (a-1)!! (b-1)!! (c-1)!!
c...  xyzi(a,b,c)= -----  I I xs  ys  zs  dP dT = -----------------------
c...               4 %pi  ] ]                          (a+b+c+1)!!
c...                      / /
c...
c...  if a, b and c are all even, or zero otherwise
c...
c...  This subroutine computes all integrals with 0 <= l <= i+j+k,
c...  with the additional constraint that l-(i+j+k) must be even
c...
c...  **** This version works for l up to 16 ****
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (18),(27-30)
c...
c...  David B. Cook, Handbook of Computational Quantum Chemistry,
c...  (Oxford, 1998), page 596-601
c...
c...  Input parameters:
c...
c...  i       exponent of xs
c...  j       exponent of ys
c...  k       exponent of zs
c...  yrk     array of real spherical harmonics evaluated at k
c...  mijk    maximum exponent value, for dimensioning
c...  xyzi    table of double factorial products
c...
c...  Output parameter:
c...
c...  omega   array of type one angular integrals
c...
      real*8 ZERO
      parameter (ZERO=0.0D0)
c...
      integer i,j,k,lmax,mijk,ii
      real*8 yrk(*),OMEGA(0:k+j+i),xyzi(0:mijk,0:mijk,0:mijk)
      logical ieven,jeven,keven
c...
c...
c...  determine whether i,j and k are even or odd
c...
      ieven=mod(i,2).eq.0
      jeven=mod(j,2).eq.0
      keven=mod(k,2).eq.0
c...
c...  compute the upper limit for l
c...  and whether we must do the even or odd series
c...
      lmax=k+j+i
      if(lmax.gt.16)
     $  call nerror(1,'omega0','maximum l value exceeded',lmax,16)
c...
c...  zero out omega
c...
      do ii=0,lmax,1
        OMEGA(ii)=ZERO
      enddo
c...
c...  integrals survive only if i and j are both even
c...
      if((ieven.and.jeven))then
        if(keven)then
c...
c...  l=0
c...
          if(lmax.ge.0)then
            OMEGA(0)=3.544907701811032D0*yrk(1)*xyzi(i,j,k)
          endif
c...
c...  l=2
c...
          if(lmax.ge.2)then
            OMEGA(2)=yrk(7)*(1.188998189281803D1*xyzi(i,j,k+2)-3.9633272
     $        97606011D0*xyzi(i,j,k))
          endif
c...
c...  l=4
c...
          if(lmax.ge.4)then
            OMEGA(4)=yrk(21)*(4.65269135862698D1*xyzi(i,j,k+4)-3.9880211
     $        64537411D1*xyzi(i,j,k+2)+3.988021164537411D0*xyzi(i,j,k))
          endif
c...
c...  l=6
c...
          if(lmax.ge.6)then
            OMEGA(6)=yrk(43)*(1.845306898868157D2*xyzi(i,j,k+6)-2.516327
     $        589365668D2*xyzi(i,j,k+4)+8.387758631218894D1*xyzi(i,j,k+2
     $        )-3.994170776770902D0*xyzi(i,j,k))
          endif
c...
c...  l=8
c...
          if(lmax.ge.8)then
            OMEGA(8)=yrk(73)*(7.347980147805839D2*xyzi(i,j,k+8)-1.371622
     $        960923757D3*xyzi(i,j,k+6)+7.91320938994475D2*xyzi(i,j,k+4)
     $        -1.438765343626318D2*xyzi(i,j,k+2)+3.996570398961995D0*xyz
     $        i(i,j,k))
          endif
c...
c...  l=10
c...
          if(lmax.ge.10)then
            OMEGA(10)=yrk(111)*(2.930982152135684D3*xyzi(i,j,k+10)-6.941
     $        799834005568D3*xyzi(i,j,k+8)+5.716776333886939D3*xyzi(i,j,
     $        k+6)-1.905592111295646D3*xyzi(i,j,k+4)+2.198760128418053D2
     $        *xyzi(i,j,k+2)-3.997745688032824D0*xyzi(i,j,k))
          endif
c...
c...  l=12
c...
          if(lmax.ge.12)then
            OMEGA(12)=yrk(157)*(1.170163993078432D4*xyzi(i,j,k+12)-3.357
     $        861893181587D4*xyzi(i,j,k+10)+3.597709171265986D4*xyzi(i,j
     $        ,k+8)-1.767295733253467D4*xyzi(i,j,k+6)+3.898446470412059D
     $        3*xyzi(i,j,k+4)-3.118757176329647D2*xyzi(i,j,k+2)+3.998406
     $        636320061D0*xyzi(i,j,k))
          endif
c...
c...  l=14
c...
          if(lmax.ge.14)then
            OMEGA(14)=yrk(211)*(4.674208812108592D4*xyzi(i,j,k+14)-1.575
     $        381488525489D5*xyzi(i,j,k+12)+2.079503564853645D5*xyzi(i,j
     $        ,k+10)-1.356197977078464D5*xyzi(i,j,k+8)+4.52065992359488D
     $        4*xyzi(i,j,k+6)-7.137884089886653D3*xyzi(i,j,k+4)+4.198755
     $        346992149D2*xyzi(i,j,k+2)-3.998814616182999D0*xyzi(i,j,k))
          endif
c...
c...  l=16
c...
          if(lmax.ge.16)then
            OMEGA(16)=yrk(273)*(1.867731876130107D5*xyzi(i,j,k+16)-7.229
     $        929843084283D5*xyzi(i,j,k+14)+1.134351061587362D6*xyzi(i,j
     $        ,k+12)-9.242860501822947D5*xyzi(i,j,k+10)+4.15928722582032
     $        6D5*xyzi(i,j,k+8)-1.012696020199732D5*xyzi(i,j,k+6)+1.2055
     $        90500237776D4*xyzi(i,j,k+4)-5.438754136411018D2*xyzi(i,j,k
     $        +2)+3.999083923831631D0*xyzi(i,j,k))
          endif
        else
c...
c...  l=1
c...
          if(lmax.ge.1)then
            OMEGA(1)=6.139960247678931D0*yrk(3)*xyzi(i,j,k+1)
          endif
c...
c...  l=3
c...
          if(lmax.ge.3)then
            OMEGA(3)=yrk(13)*(2.344736049917376D1*xyzi(i,j,k+3)-1.406841
     $        629950425D1*xyzi(i,j,k+1))
          endif
c...
c...  l=5
c...
          if(lmax.ge.5)then
            OMEGA(5)=yrk(31)*(9.258738901136752D1*xyzi(i,j,k+5)-1.028748
     $        766792972D2*xyzi(i,j,k+3)+2.204461643127798D1*xyzi(i,j,k+1
     $        ))
          endif
c...
c...  l=7
c...
          if(lmax.ge.7)then
            OMEGA(7)=yrk(57)*(3.681186927173971D2*xyzi(i,j,k+7)-5.946532
     $        728511799D2*xyzi(i,j,k+5)+2.702969422050818D2*xyzi(i,j,k+3
     $        )-3.003299357834242D1*xyzi(i,j,k+1))
          endif
c...
c...  l=9
c...
          if(lmax.ge.9)then
            OMEGA(9)=yrk(91)*(1.467326381829043D3*xyzi(i,j,k+9)-3.107279
     $        396814444D3*xyzi(i,j,k+7)+2.175095577770111D3*xyzi(i,j,k+5
     $        )-5.57716814812849D2*xyzi(i,j,k+3)+3.802614646451243D1*xyz
     $        i(i,j,k+1))
          endif
c...
c...  l=11
c...
          if(lmax.ge.11)then
            OMEGA(11)=yrk(133)*(5.855905424818466D3*xyzi(i,j,k+11)-1.533
     $        689516023884D4*xyzi(i,j,k+9)+1.452969015180522D4*xyzi(i,j,
     $        k+7)-5.982813591919795D3*xyzi(i,j,k+5)+9.971355986532992D2
     $        *xyzi(i,j,k+3)-4.602164301476766D1*xyzi(i,j,k+1))
          endif
c...
c...  l=13
c...
          if(lmax.ge.13)then
            OMEGA(13)=yrk(183)*(2.338596333691754D4*xyzi(i,j,k+13)-7.296
     $        420561118272D4*xyzi(i,j,k+11)+8.72398110568489D4*xyzi(i,j,
     $        k+9)-4.985132060391366D4*xyzi(i,j,k+7)+1.377470700897614D4
     $        *xyzi(i,j,k+5)-1.620553765761899D3*xyzi(i,j,k+3)+5.4018458
     $        85872997D1*xyzi(i,j,k+1))
          endif
c...
c...  l=15
c...
          if(lmax.ge.15)then
            OMEGA(15)=yrk(241)*(9.343222615411324D4*xyzi(i,j,k+15)-3.382
     $        890946959272D5*xyzi(i,j,k+13)+4.886398034496727D5*xyzi(i,j
     $        ,k+11)-3.583358558630933D5*xyzi(i,j,k+9)+1.402183783812104
     $        D5*xyzi(i,j,k+7)-2.804367567624208D4*xyzi(i,j,k+5)+2.45997
     $        1550547551D3*xyzi(i,j,k+3)-6.201608950960213D1*xyzi(i,j,k+
     $        1))
          endif
        endif
      endif
c...
      end
c=======================================================================
      subroutine omega1(i,j,k,yrk,mijk,xyzi,OMEGA)
      implicit none
c...
c...  Max2fort (MM) version 0.0 (June 2004)
c...  translated on 9/15/2004 15:43:35
c...
c...  This subroutine computes type 1 angular integrals over pseudopotentials
c...
c...
c...   omega(i,j,k,l)=
c...
c...   l
c...  ====                    / /
c...  \                       [ [                     i   j   k
c...   >  yr(m,l,thetak,phik) I I yr(m,l,THETA,PHI) xs  ys  zs  dPHi dTHETA
c...  /                       ] ]
c...  ====                    / /
c...  m = - l
c...
c...  where yr are real spherical harmonic functions, and xs, ys and zs
c...  are cartesian coordinates costrained over the surface of the unit
c...  sphere:
c...              xs=sin(THETA) cos(PHI)
c...              ys=sin(THETA) sin(PHI)
c...              zs=cos(THETA)
c...
c...  The integral is computed by expanding the real sperical harmonic
c...  as a polynomial in xs, ys and zs, obtaining integrals of the type:
c...
c...                      / /
c...                 1    [ [   a   b   c         (a-1)!! (b-1)!! (c-1)!!
c...  xyzi(a,b,c)= -----  I I xs  ys  zs  dP dT = -----------------------
c...               4 %pi  ] ]                          (a+b+c+1)!!
c...                      / /
c...
c...  if a, b and c are all even, or zero otherwise.
c...
c...  This subroutine computes all integrals with 0 <= l <= i+j+k,
c...  with the additional constraint that l-(i+j+k) must be even
c...
c...  **** This version works for l up to 16 ****
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (18),(27-30)
c...
c...  David B. Cook, Handbook of Computational Quantum Chemistry,
c...  (Oxford, 1998), page 596-601
c...
c...  Input parameters:
c...
c...  i       exponent of xs
c...  j       exponent of ys
c...  k       exponent of zs
c...  yrk     array of real spherical harmonics evaluated at k
c...  mijk    maximum exponent value, for dimensioning
c...  xyzi    table of double factorial products
c...
c...  Output parameter:
c...
c...  omega   array of type one angular integrals
c...
      real*8 ZERO
      parameter (ZERO=0.0D0)
c...
      integer i,j,k,lmax,mijk,ii
      real*8 yrk(*),OMEGA(0:k+j+i),xyzi(0:mijk,0:mijk,0:mijk),angi
      logical leven,ieven,jeven,keven
c...
c...
c...  determine whether i,j and k are even or odd
c...
      ieven=mod(i,2).eq.0
      jeven=mod(j,2).eq.0
      keven=mod(k,2).eq.0
c...
c...  compute the upper limit for l
c...  and whether we must do the even or odd series
c...
      lmax=k+j+i
      if(lmax.gt.16)
     $  call nerror(1,'omega1','maximum l value exceeded',lmax,16)
      leven=mod(lmax,2).eq.0
c...
c...  zero out omega
c...
      do ii=0,lmax,1
        OMEGA(ii)=ZERO
      enddo
c...
c...  chose even or odd series
c...
      if(leven)then
c...
c...  l=0
c...
        if(lmax.ge.0)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=3.544907701811032D0*yrk(1)*xyzi(i,j,k)+angi
          endif
          OMEGA(0)=angi
        endif
c...
c...  l=2
c...
        if(lmax.ge.2)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(7)*(1.188998189281803D1*xyzi(i,j,k+2)-3.96332729760
     $        6011D0*xyzi(i,j,k))+angi
            angi=yrk(9)*(6.864684246478268D0*xyzi(i+2,j,k)-6.86468424647
     $        8268D0*xyzi(i,j+2,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=1.372936849295654D1*yrk(5)*xyzi(i+1,j+1,k)+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=angi-1.372936849295654D1*yrk(8)*xyzi(i+1,j,k+1)
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=angi-1.372936849295654D1*yrk(6)*xyzi(i,j+1,k+1)
          endif
          OMEGA(2)=angi
        endif
c...
c...  l=4
c...
        if(lmax.ge.4)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(21)*(4.65269135862698D1*xyzi(i,j,k+4)-3.98802116453
     $        7411D1*xyzi(i,j,k+2)+3.988021164537411D0*xyzi(i,j,k))+angi
            angi=yrk(23)*(-4.161493662486312D1*xyzi(i+4,j,k)+3.566994567
     $        84541D1*xyzi(i+2,j,k)+4.161493662486312D1*xyzi(i,j+4,k)-3.
     $        56699456784541D1*xyzi(i,j+2,k))+angi
            angi=yrk(25)*(7.864483795364389D0*xyzi(i+4,j,k)-4.7186902772
     $        18633D1*xyzi(i+2,j+2,k)+7.864483795364389D0*xyzi(i,j+4,k))
     $        +angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(17)*(3.145793518145755D1*xyzi(i+3,j+1,k)-3.14579351
     $        8145755D1*xyzi(i+1,j+3,k))+angi
            angi=yrk(19)*(-8.322987324972623D1*xyzi(i+3,j+1,k)-8.3229873
     $        24972623D1*xyzi(i+1,j+3,k)+7.13398913569082D1*xyzi(i+1,j+1
     $        ,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(22)*(5.885240777217825D1*xyzi(i+3,j,k+1)+5.88524077
     $        7217825D1*xyzi(i+1,j+2,k+1)-3.362994729838757D1*xyzi(i+1,j
     $        ,k+1))+angi
            angi=yrk(24)*(6.67323578668065D1*xyzi(i+1,j+2,k+1)-2.2244119
     $        2889355D1*xyzi(i+3,j,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(18)*(2.22441192889355D1*xyzi(i,j+3,k+1)-6.673235786
     $        68065D1*xyzi(i+2,j+1,k+1))+angi
            angi=yrk(20)*(5.885240777217825D1*xyzi(i+2,j+1,k+1)+5.885240
     $        777217825D1*xyzi(i,j+3,k+1)-3.362994729838757D1*xyzi(i,j+1
     $        ,k+1))+angi
          endif
          OMEGA(4)=angi
        endif
c...
c...  l=6
c...
        if(lmax.ge.6)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(43)*(1.845306898868157D2*xyzi(i,j,k+6)-2.5163275893
     $        65668D2*xyzi(i,j,k+4)+8.387758631218894D1*xyzi(i,j,k+2)-3.
     $        994170776770902D0*xyzi(i,j,k))+angi
            angi=yrk(45)*(1.910074105988639D2*xyzi(i+6,j,k)+1.9100741059
     $        88639D2*xyzi(i+4,j+2,k)-2.778289608710748D2*xyzi(i+4,j,k)-
     $        1.910074105988639D2*xyzi(i+2,j+4,k)+9.260965362369161D1*xy
     $        zi(i+2,j,k)-1.910074105988639D2*xyzi(i,j+6,k)+2.7782896087
     $        10748D2*xyzi(i,j+4,k)-9.260965362369161D1*xyzi(i,j+2,k))+a
     $        ngi
            angi=yrk(47)*(-6.974604495709942D1*xyzi(i+6,j,k)+3.487302247
     $        854971D2*xyzi(i+4,j+2,k)+6.340549541554493D1*xyzi(i+4,j,k)
     $        +3.487302247854971D2*xyzi(i+2,j+4,k)-3.804329724932696D2*x
     $        yzi(i+2,j+2,k)-6.974604495709942D1*xyzi(i,j+6,k)+6.3405495
     $        41554493D1*xyzi(i,j+4,k))+angi
            angi=yrk(49)*(8.585144663680938D0*xyzi(i+6,j,k)-1.2877716995
     $        52141D2*xyzi(i+4,j+2,k)+1.287771699552141D2*xyzi(i+2,j+4,k
     $        )-8.585144663680938D0*xyzi(i,j+6,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(37)*(5.151086798208563D1*xyzi(i+5,j+1,k)-1.71702893
     $        2736188D2*xyzi(i+3,j+3,k)+5.151086798208563D1*xyzi(i+1,j+5
     $        ,k))+angi
            angi=yrk(39)*(-2.789841798283977D2*xyzi(i+5,j+1,k)+2.5362198
     $        16621797D2*xyzi(i+3,j+1,k)+2.789841798283977D2*xyzi(i+1,j+
     $        5,k)-2.536219816621797D2*xyzi(i+1,j+3,k))+angi
            angi=yrk(41)*(3.820148211977279D2*xyzi(i+5,j+1,k)+7.64029642
     $        3954557D2*xyzi(i+3,j+3,k)-5.556579217421496D2*xyzi(i+3,j+1
     $        ,k)+3.820148211977279D2*xyzi(i+1,j+5,k)-5.556579217421496D
     $        2*xyzi(i+1,j+3,k)+1.852193072473832D2*xyzi(i+1,j+1,k))+ang
     $        i
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(44)*(-2.416073869853585D2*xyzi(i+5,j,k+1)-4.8321477
     $        39707171D2*xyzi(i+3,j+2,k+1)+2.635716948931184D2*xyzi(i+3,
     $        j,k+1)-2.416073869853585D2*xyzi(i+1,j+4,k+1)+2.63571694893
     $        1184D2*xyzi(i+1,j+2,k+1)-5.857148775402631D1*xyzi(i+1,j,k+
     $        1))+angi
            angi=yrk(46)*(1.27338273732576D2*xyzi(i+5,j,k+1)-2.546765474
     $        651519D2*xyzi(i+3,j+2,k+1)-9.260965362369161D1*xyzi(i+3,j,
     $        k+1)-3.820148211977279D2*xyzi(i+1,j+4,k+1)+2.7782896087107
     $        48D2*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(48)*(-2.973981349564841D1*xyzi(i+5,j,k+1)+2.9739813
     $        49564841D2*xyzi(i+3,j+2,k+1)-1.486990674782421D2*xyzi(i+1,
     $        j+4,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(38)*(-1.486990674782421D2*xyzi(i+4,j+1,k+1)+2.97398
     $        1349564841D2*xyzi(i+2,j+3,k+1)-2.973981349564841D1*xyzi(i,
     $        j+5,k+1))+angi
            angi=yrk(40)*(3.820148211977279D2*xyzi(i+4,j+1,k+1)+2.546765
     $        474651519D2*xyzi(i+2,j+3,k+1)-2.778289608710748D2*xyzi(i+2
     $        ,j+1,k+1)-1.27338273732576D2*xyzi(i,j+5,k+1)+9.26096536236
     $        9161D1*xyzi(i,j+3,k+1))+angi
            angi=yrk(42)*(-2.416073869853585D2*xyzi(i+4,j+1,k+1)-4.83214
     $        7739707171D2*xyzi(i+2,j+3,k+1)+2.635716948931184D2*xyzi(i+
     $        2,j+1,k+1)-2.416073869853585D2*xyzi(i,j+5,k+1)+2.635716948
     $        931184D2*xyzi(i,j+3,k+1)-5.857148775402631D1*xyzi(i,j+1,k+
     $        1))+angi
          endif
          OMEGA(6)=angi
        endif
c...
c...  l=8
c...
        if(lmax.ge.8)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(73)*(7.347980147805839D2*xyzi(i,j,k+8)-1.3716229609
     $        23757D3*xyzi(i,j,k+6)+7.91320938994475D2*xyzi(i,j,k+4)-1.4
     $        38765343626318D2*xyzi(i,j,k+2)+3.996570398961995D0*xyzi(i,
     $        j,k))+angi
            angi=yrk(75)*(-8.197015020580125D2*xyzi(i+8,j,k)-1.639403004
     $        116025D3*xyzi(i+6,j+2,k)+1.639403004116025D3*xyzi(i+6,j,k)
     $        +1.639403004116025D3*xyzi(i+4,j+2,k)-1.008863387148323D3*x
     $        yzi(i+4,j,k)+1.639403004116025D3*xyzi(i+2,j+6,k)-1.6394030
     $        04116025D3*xyzi(i+2,j+4,k)+1.834297067542406D2*xyzi(i+2,j,
     $        k)+8.197015020580125D2*xyzi(i,j+8,k)-1.639403004116025D3*x
     $        yzi(i,j+6,k)+1.008863387148323D3*xyzi(i,j+4,k)-1.834297067
     $        542406D2*xyzi(i,j+2,k))+angi
            angi=yrk(77)*(3.907773582803669D2*xyzi(i+8,j,k)-1.5631094331
     $        21468D3*xyzi(i+6,j+2,k)-6.252437732485871D2*xyzi(i+6,j,k)-
     $        3.907773582803669D3*xyzi(i+4,j+4,k)+3.126218866242935D3*xy
     $        zi(i+4,j+2,k)+2.404783743263796D2*xyzi(i+4,j,k)-1.56310943
     $        3121468D3*xyzi(i+2,j+6,k)+3.126218866242935D3*xyzi(i+2,j+4
     $        ,k)-1.442870245958278D3*xyzi(i+2,j+2,k)+3.907773582803669D
     $        2*xyzi(i,j+8,k)-6.252437732485871D2*xyzi(i,j+6,k)+2.404783
     $        743263796D2*xyzi(i,j+4,k))+angi
            angi=yrk(79)*(-1.003423624270676D2*xyzi(i+8,j,k)+1.404793073
     $        978946D3*xyzi(i+6,j+2,k)+9.365287159859641D1*xyzi(i+6,j,k)
     $        -1.404793073978946D3*xyzi(i+4,j+2,k)-1.404793073978946D3*x
     $        yzi(i+2,j+6,k)+1.404793073978946D3*xyzi(i+2,j+4,k)+1.00342
     $        3624270676D2*xyzi(i,j+8,k)-9.365287159859641D1*xyzi(i,j+6,
     $        k))+angi
            angi=yrk(81)*(9.159962562443957D0*xyzi(i+8,j,k)-2.5647895174
     $        84308D2*xyzi(i+6,j+2,k)+6.41197379371077D2*xyzi(i+4,j+4,k)
     $        -2.564789517484308D2*xyzi(i+2,j+6,k)+9.159962562443957D0*x
     $        yzi(i,j+8,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(65)*(7.327970049955166D1*xyzi(i+7,j+1,k)-5.12957903
     $        4968616D2*xyzi(i+5,j+3,k)+5.129579034968616D2*xyzi(i+3,j+5
     $        ,k)-7.327970049955166D1*xyzi(i+1,j+7,k))+angi
            angi=yrk(67)*(-6.020541745624055D2*xyzi(i+7,j+1,k)+1.4047930
     $        73978946D3*xyzi(i+5,j+3,k)+5.619172295915784D2*xyzi(i+5,j+
     $        1,k)+1.404793073978946D3*xyzi(i+3,j+5,k)-1.873057431971928
     $        D3*xyzi(i+3,j+3,k)-6.020541745624055D2*xyzi(i+1,j+7,k)+5.6
     $        19172295915784D2*xyzi(i+1,j+5,k))+angi
            angi=yrk(69)*(1.563109433121468D3*xyzi(i+7,j+1,k)+1.56310943
     $        3121468D3*xyzi(i+5,j+3,k)-2.500975092994348D3*xyzi(i+5,j+1
     $        ,k)-1.563109433121468D3*xyzi(i+3,j+5,k)+9.619134973055186D
     $        2*xyzi(i+3,j+1,k)-1.563109433121468D3*xyzi(i+1,j+7,k)+2.50
     $        0975092994348D3*xyzi(i+1,j+5,k)-9.619134973055186D2*xyzi(i
     $        +1,j+3,k))+angi
            angi=yrk(71)*(-1.639403004116025D3*xyzi(i+7,j+1,k)-4.9182090
     $        12348075D3*xyzi(i+5,j+3,k)+3.27880600823205D3*xyzi(i+5,j+1
     $        ,k)-4.918209012348075D3*xyzi(i+3,j+5,k)+6.5576120164641D3*
     $        xyzi(i+3,j+3,k)-2.017726774296646D3*xyzi(i+3,j+1,k)-1.6394
     $        03004116025D3*xyzi(i+1,j+7,k)+3.27880600823205D3*xyzi(i+1,
     $        j+5,k)-2.017726774296646D3*xyzi(i+1,j+3,k)+3.6685941350848
     $        11D2*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(74)*(9.797306863741119D2*xyzi(i+7,j,k+1)+2.93919205
     $        9122336D3*xyzi(i+5,j+2,k+1)-1.567569098198579D3*xyzi(i+5,j
     $        ,k+1)+2.939192059122336D3*xyzi(i+3,j+4,k+1)-3.135138196397
     $        158D3*xyzi(i+3,j+2,k+1)+7.234934299378057D2*xyzi(i+3,j,k+1
     $        )+9.797306863741119D2*xyzi(i+1,j+6,k+1)-1.567569098198579D
     $        3*xyzi(i+1,j+4,k+1)+7.234934299378057D2*xyzi(i+1,j+2,k+1)-
     $        8.769617332579463D1*xyzi(i+1,j,k+1))+angi
            angi=yrk(76)*(-6.05389680277916D2*xyzi(i+7,j,k+1)+6.05389680
     $        277916D2*xyzi(i+5,j+2,k+1)+8.071862403705547D2*xyzi(i+5,j,
     $        k+1)+3.02694840138958D3*xyzi(i+3,j+4,k+1)-1.61437248074110
     $        9D3*xyzi(i+3,j+2,k+1)-2.483649970370938D2*xyzi(i+3,j,k+1)+
     $        1.816169040833748D3*xyzi(i+1,j+6,k+1)-2.421558721111664D3*
     $        xyzi(i+1,j+4,k+1)+7.450949911112813D2*xyzi(i+1,j+2,k+1))+a
     $        ngi
            angi=yrk(78)*(2.167642773184962D2*xyzi(i+7,j,k+1)-1.95087849
     $        5866466D3*xyzi(i+5,j+2,k+1)-1.73411421854797D2*xyzi(i+5,j,
     $        k+1)-1.083821386592481D3*xyzi(i+3,j+4,k+1)+1.7341142185479
     $        7D3*xyzi(i+3,j+2,k+1)+1.083821386592481D3*xyzi(i+1,j+6,k+1
     $        )-8.670571092739848D2*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(80)*(-3.663985024977583D1*xyzi(i+7,j,k+1)+7.6943685
     $        52452924D2*xyzi(i+5,j+2,k+1)-1.282394758742154D3*xyzi(i+3,
     $        j+4,k+1)+2.564789517484308D2*xyzi(i+1,j+6,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(66)*(-2.564789517484308D2*xyzi(i+6,j+1,k+1)+1.28239
     $        4758742154D3*xyzi(i+4,j+3,k+1)-7.694368552452924D2*xyzi(i+
     $        2,j+5,k+1)+3.663985024977583D1*xyzi(i,j+7,k+1))+angi
            angi=yrk(68)*(1.083821386592481D3*xyzi(i+6,j+1,k+1)-1.083821
     $        386592481D3*xyzi(i+4,j+3,k+1)-8.670571092739848D2*xyzi(i+4
     $        ,j+1,k+1)-1.950878495866466D3*xyzi(i+2,j+5,k+1)+1.73411421
     $        854797D3*xyzi(i+2,j+3,k+1)+2.167642773184962D2*xyzi(i,j+7,
     $        k+1)-1.73411421854797D2*xyzi(i,j+5,k+1))+angi
            angi=yrk(70)*(-1.816169040833748D3*xyzi(i+6,j+1,k+1)-3.02694
     $        840138958D3*xyzi(i+4,j+3,k+1)+2.421558721111664D3*xyzi(i+4
     $        ,j+1,k+1)-6.05389680277916D2*xyzi(i+2,j+5,k+1)+1.614372480
     $        741109D3*xyzi(i+2,j+3,k+1)-7.450949911112813D2*xyzi(i+2,j+
     $        1,k+1)+6.05389680277916D2*xyzi(i,j+7,k+1)-8.07186240370554
     $        7D2*xyzi(i,j+5,k+1)+2.483649970370938D2*xyzi(i,j+3,k+1))+a
     $        ngi
            angi=yrk(72)*(9.797306863741119D2*xyzi(i+6,j+1,k+1)+2.939192
     $        059122336D3*xyzi(i+4,j+3,k+1)-1.567569098198579D3*xyzi(i+4
     $        ,j+1,k+1)+2.939192059122336D3*xyzi(i+2,j+5,k+1)-3.13513819
     $        6397158D3*xyzi(i+2,j+3,k+1)+7.234934299378057D2*xyzi(i+2,j
     $        +1,k+1)+9.797306863741119D2*xyzi(i,j+7,k+1)-1.567569098198
     $        579D3*xyzi(i,j+5,k+1)+7.234934299378057D2*xyzi(i,j+3,k+1)-
     $        8.769617332579463D1*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA(8)=angi
        endif
c...
c...  l=10
c...
        if(lmax.ge.10)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(111)*(2.930982152135684D3*xyzi(i,j,k+10)-6.94179983
     $        4005568D3*xyzi(i,j,k+8)+5.716776333886939D3*xyzi(i,j,k+6)-
     $        1.905592111295646D3*xyzi(i,j,k+4)+2.198760128418053D2*xyzi
     $        (i,j,k+2)-3.997745688032824D0*xyzi(i,j,k))+angi
            angi=yrk(113)*(3.422649766190848D3*xyzi(i+10,j,k)+1.02679492
     $        9857254D4*xyzi(i+8,j+2,k)-8.646694146166353D3*xyzi(i+8,j,k
     $        )+6.845299532381696D3*xyzi(i+6,j+4,k)-1.729338829233271D4*
     $        xyzi(i+6,j+2,k)+7.629436011323252D3*xyzi(i+6,j,k)-6.845299
     $        532381696D3*xyzi(i+4,j+6,k)+7.629436011323252D3*xyzi(i+4,j
     $        +2,k)-2.712688359581601D3*xyzi(i+4,j,k)-1.026794929857254D
     $        4*xyzi(i+2,j+8,k)+1.729338829233271D4*xyzi(i+2,j+6,k)-7.62
     $        9436011323252D3*xyzi(i+2,j+4,k)+3.130025030286462D2*xyzi(i
     $        +2,j,k)-3.422649766190848D3*xyzi(i,j+10,k)+8.6466941461663
     $        53D3*xyzi(i,j+8,k)-7.629436011323252D3*xyzi(i,j+6,k)+2.712
     $        688359581601D3*xyzi(i,j+4,k)-3.130025030286462D2*xyzi(i,j+
     $        2,k))+angi
            angi=yrk(115)*(-1.898544496916298D3*xyzi(i+10,j,k)+5.6956334
     $        90748894D3*xyzi(i+8,j+2,k)+4.196782572130764D3*xyzi(i+8,j,
     $        k)+2.657962295682817D4*xyzi(i+6,j+4,k)-1.678713028852305D4
     $        *xyzi(i+6,j+2,k)-2.962434756798186D3*xyzi(i+6,j,k)+2.65796
     $        2295682817D4*xyzi(i+4,j+6,k)-4.196782572130764D4*xyzi(i+4,
     $        j+4,k)+1.481217378399093D4*xyzi(i+4,j+2,k)+6.5831883484404
     $        14D2*xyzi(i+4,j,k)+5.695633490748894D3*xyzi(i+2,j+8,k)-1.6
     $        78713028852305D4*xyzi(i+2,j+6,k)+1.481217378399093D4*xyzi(
     $        i+2,j+4,k)-3.949913009064248D3*xyzi(i+2,j+2,k)-1.898544496
     $        916298D3*xyzi(i,j+10,k)+4.196782572130764D3*xyzi(i,j+8,k)-
     $        2.962434756798186D3*xyzi(i,j+6,k)+6.583188348440414D2*xyzi
     $        (i,j+4,k))+angi
            angi=yrk(117)*(6.712368440769583D2*xyzi(i+10,j,k)-8.72607897
     $        3000458D3*xyzi(i+8,j+2,k)-1.130504158445403D3*xyzi(i+8,j,k
     $        )-9.397315817077416D3*xyzi(i+6,j+4,k)+1.582705821823565D4*
     $        xyzi(i+6,j+2,k)+4.655017123010485D2*xyzi(i+6,j,k)+9.397315
     $        817077416D3*xyzi(i+4,j+6,k)-6.982525684515727D3*xyzi(i+4,j
     $        +2,k)+8.726078973000458D3*xyzi(i+2,j+8,k)-1.58270582182356
     $        5D4*xyzi(i+2,j+6,k)+6.982525684515727D3*xyzi(i+2,j+4,k)-6.
     $        712368440769583D2*xyzi(i,j+10,k)+1.130504158445403D3*xyzi(
     $        i,j+8,k)-4.655017123010485D2*xyzi(i,j+6,k))+angi
            angi=yrk(119)*(-1.329247023836435D2*xyzi(i+10,j,k)+3.5889669
     $        64358373D3*xyzi(i+8,j+2,k)+1.259286654160833D2*xyzi(i+8,j,
     $        k)-5.582837500113025D3*xyzi(i+6,j+4,k)-3.526002631650332D3
     $        *xyzi(i+6,j+2,k)-5.582837500113025D3*xyzi(i+4,j+6,k)+8.815
     $        006579125829D3*xyzi(i+4,j+4,k)+3.588966964358373D3*xyzi(i+
     $        2,j+8,k)-3.526002631650332D3*xyzi(i+2,j+6,k)-1.32924702383
     $        6435D2*xyzi(i,j+10,k)+1.259286654160833D2*xyzi(i,j+8,k))+a
     $        ngi
            angi=yrk(121)*(9.643371463227499D0*xyzi(i+10,j,k)-4.33951715
     $        8452375D2*xyzi(i+8,j+2,k)+2.025108007277775D3*xyzi(i+6,j+4
     $        ,k)-2.025108007277775D3*xyzi(i+4,j+6,k)+4.339517158452375D
     $        2*xyzi(i+2,j+8,k)-9.643371463227499D0*xyzi(i,j+10,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(101)*(9.643371463227499D1*xyzi(i+9,j+1,k)-1.1572045
     $        755873D3*xyzi(i+7,j+3,k)+2.43012960873333D3*xyzi(i+5,j+5,k
     $        )-1.1572045755873D3*xyzi(i+3,j+7,k)+9.643371463227499D1*xy
     $        zi(i+1,j+9,k))+angi
            angi=yrk(103)*(-1.063397619069148D3*xyzi(i+9,j+1,k)+6.380385
     $        714414886D3*xyzi(i+7,j+3,k)+1.007429323328666D3*xyzi(i+7,j
     $        +1,k)-7.052005263300663D3*xyzi(i+5,j+3,k)-6.38038571441488
     $        6D3*xyzi(i+3,j+7,k)+7.052005263300663D3*xyzi(i+3,j+5,k)+1.
     $        063397619069148D3*xyzi(i+1,j+9,k)-1.007429323328666D3*xyzi
     $        (i+1,j+7,k))+angi
            angi=yrk(105)*(4.02742106446175D3*xyzi(i+9,j+1,k)-5.36989475
     $        2615666D3*xyzi(i+7,j+3,k)-6.783024950672421D3*xyzi(i+7,j+1
     $        ,k)-1.879463163415483D4*xyzi(i+5,j+5,k)+1.582705821823565D
     $        4*xyzi(i+5,j+3,k)+2.793010273806291D3*xyzi(i+5,j+1,k)-5.36
     $        9894752615666D3*xyzi(i+3,j+7,k)+1.582705821823565D4*xyzi(i
     $        +3,j+5,k)-9.31003424602097D3*xyzi(i+3,j+3,k)+4.02742106446
     $        175D3*xyzi(i+1,j+9,k)-6.783024950672421D3*xyzi(i+1,j+7,k)+
     $        2.793010273806291D3*xyzi(i+1,j+5,k))+angi
            angi=yrk(107)*(-7.594177987665191D3*xyzi(i+9,j+1,k)-1.518835
     $        597533038D4*xyzi(i+7,j+3,k)+1.678713028852305D4*xyzi(i+7,j
     $        +1,k)+1.678713028852305D4*xyzi(i+5,j+3,k)-1.18497390271927
     $        4D4*xyzi(i+5,j+1,k)+1.518835597533038D4*xyzi(i+3,j+7,k)-1.
     $        678713028852305D4*xyzi(i+3,j+5,k)+2.633275339376166D3*xyzi
     $        (i+3,j+1,k)+7.594177987665191D3*xyzi(i+1,j+9,k)-1.67871302
     $        8852305D4*xyzi(i+1,j+7,k)+1.184973902719274D4*xyzi(i+1,j+5
     $        ,k)-2.633275339376166D3*xyzi(i+1,j+3,k))+angi
            angi=yrk(109)*(6.845299532381696D3*xyzi(i+9,j+1,k)+2.7381198
     $        12952678D4*xyzi(i+7,j+3,k)-1.729338829233271D4*xyzi(i+7,j+
     $        1,k)+4.107179719429018D4*xyzi(i+5,j+5,k)-5.188016487699811
     $        D4*xyzi(i+5,j+3,k)+1.52588720226465D4*xyzi(i+5,j+1,k)+2.73
     $        8119812952678D4*xyzi(i+3,j+7,k)-5.188016487699811D4*xyzi(i
     $        +3,j+5,k)+3.051774404529301D4*xyzi(i+3,j+3,k)-5.4253767191
     $        63202D3*xyzi(i+3,j+1,k)+6.845299532381696D3*xyzi(i+1,j+9,k
     $        )-1.729338829233271D4*xyzi(i+1,j+7,k)+1.52588720226465D4*x
     $        yzi(i+1,j+5,k)-5.425376719163202D3*xyzi(i+1,j+3,k)+6.26005
     $        0060572925D2*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(112)*(-3.952135527704191D3*xyzi(i+9,j,k+1)-1.580854
     $        211081677D4*xyzi(i+7,j+2,k+1)+8.320285321482508D3*xyzi(i+7
     $        ,j,k+1)-2.371281316622515D4*xyzi(i+5,j+4,k+1)+2.4960855964
     $        44752D4*xyzi(i+5,j+2,k+1)-5.873142579870006D3*xyzi(i+5,j,k
     $        +1)-1.580854211081677D4*xyzi(i+3,j+6,k+1)+2.49608559644475
     $        2D4*xyzi(i+3,j+4,k+1)-1.174628515974001D4*xyzi(i+3,j+2,k+1
     $        )+1.566171354632002D3*xyzi(i+3,j,k+1)-3.952135527704191D3*
     $        xyzi(i+1,j+8,k+1)+8.320285321482508D3*xyzi(i+1,j+6,k+1)-5.
     $        873142579870006D3*xyzi(i+1,j+4,k+1)+1.566171354632002D3*xy
     $        zi(i+1,j+2,k+1)-1.20474719587077D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(114)*(2.684947376307833D3*xyzi(i+9,j,k+1)-5.0872687
     $        13004316D3*xyzi(i+7,j,k+1)-1.6109684257847D4*xyzi(i+5,j+4,
     $        k+1)+5.087268713004316D3*xyzi(i+5,j+2,k+1)+2.9925110076495
     $        97D3*xyzi(i+5,j,k+1)-2.147957901046267D4*xyzi(i+3,j+6,k+1)
     $        +2.543634356502158D4*xyzi(i+3,j+4,k+1)-5.985022015299195D3
     $        *xyzi(i+3,j+2,k+1)-5.32001956915484D2*xyzi(i+3,j,k+1)-8.05
     $        48421289235D3*xyzi(i+1,j+8,k+1)+1.526180613901295D4*xyzi(i
     $        +1,j+6,k+1)-8.977533022948792D3*xyzi(i+1,j+4,k+1)+1.596005
     $        870746452D3*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(116)*(-1.200744969886805D3*xyzi(i+9,j,k+1)+9.605959
     $        759094437D3*xyzi(i+7,j+2,k+1)+1.769518902991081D3*xyzi(i+7
     $        ,j,k+1)+1.681042957841527D4*xyzi(i+5,j+4,k+1)-1.5925670126
     $        91973D4*xyzi(i+5,j+2,k+1)-6.245360834086167D2*xyzi(i+5,j,k
     $        +1)-8.847594514955403D3*xyzi(i+3,j+4,k+1)+6.24536083408616
     $        7D3*xyzi(i+3,j+2,k+1)-6.003724849434023D3*xyzi(i+1,j+8,k+1
     $        )+8.847594514955403D3*xyzi(i+1,j+6,k+1)-3.122680417043083D
     $        3*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(118)*(3.255976950512413D2*xyzi(i+9,j,k+1)-6.5119539
     $        01024826D3*xyzi(i+7,j+2,k+1)-2.741875326747295D2*xyzi(i+7,
     $        j,k+1)+4.558367730717379D3*xyzi(i+5,j+4,k+1)+5.75793818616
     $        932D3*xyzi(i+5,j+2,k+1)+9.116735461434757D3*xyzi(i+3,j+6,k
     $        +1)-9.596563643615534D3*xyzi(i+3,j+4,k+1)-2.27918386535868
     $        9D3*xyzi(i+1,j+8,k+1)+1.919312728723107D3*xyzi(i+1,j+6,k+1
     $        ))+angi
            angi=yrk(120)*(-4.31264682481166D1*xyzi(i+9,j,k+1)+1.5525528
     $        56932198D3*xyzi(i+7,j+2,k+1)-5.433934999262692D3*xyzi(i+5,
     $        j+4,k+1)+3.622623332841795D3*xyzi(i+3,j+6,k+1)-3.881382142
     $        330494D2*xyzi(i+1,j+8,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(102)*(-3.881382142330494D2*xyzi(i+8,j+1,k+1)+3.6226
     $        23332841795D3*xyzi(i+6,j+3,k+1)-5.433934999262692D3*xyzi(i
     $        +4,j+5,k+1)+1.552552856932198D3*xyzi(i+2,j+7,k+1)-4.312646
     $        82481166D1*xyzi(i,j+9,k+1))+angi
            angi=yrk(104)*(2.279183865358689D3*xyzi(i+8,j+1,k+1)-9.11673
     $        5461434757D3*xyzi(i+6,j+3,k+1)-1.919312728723107D3*xyzi(i+
     $        6,j+1,k+1)-4.558367730717379D3*xyzi(i+4,j+5,k+1)+9.5965636
     $        43615534D3*xyzi(i+4,j+3,k+1)+6.511953901024826D3*xyzi(i+2,
     $        j+7,k+1)-5.75793818616932D3*xyzi(i+2,j+5,k+1)-3.2559769505
     $        12413D2*xyzi(i,j+9,k+1)+2.741875326747295D2*xyzi(i,j+7,k+1
     $        ))+angi
            angi=yrk(106)*(-6.003724849434023D3*xyzi(i+8,j+1,k+1)+8.8475
     $        94514955403D3*xyzi(i+6,j+1,k+1)+1.681042957841527D4*xyzi(i
     $        +4,j+5,k+1)-8.847594514955403D3*xyzi(i+4,j+3,k+1)-3.122680
     $        417043083D3*xyzi(i+4,j+1,k+1)+9.605959759094437D3*xyzi(i+2
     $        ,j+7,k+1)-1.592567012691973D4*xyzi(i+2,j+5,k+1)+6.24536083
     $        4086167D3*xyzi(i+2,j+3,k+1)-1.200744969886805D3*xyzi(i,j+9
     $        ,k+1)+1.769518902991081D3*xyzi(i,j+7,k+1)-6.24536083408616
     $        7D2*xyzi(i,j+5,k+1))+angi
            angi=yrk(108)*(8.0548421289235D3*xyzi(i+8,j+1,k+1)+2.1479579
     $        01046267D4*xyzi(i+6,j+3,k+1)-1.526180613901295D4*xyzi(i+6,
     $        j+1,k+1)+1.6109684257847D4*xyzi(i+4,j+5,k+1)-2.54363435650
     $        2158D4*xyzi(i+4,j+3,k+1)+8.977533022948792D3*xyzi(i+4,j+1,
     $        k+1)-5.087268713004316D3*xyzi(i+2,j+5,k+1)+5.9850220152991
     $        95D3*xyzi(i+2,j+3,k+1)-1.596005870746452D3*xyzi(i+2,j+1,k+
     $        1)-2.684947376307833D3*xyzi(i,j+9,k+1)+5.087268713004316D3
     $        *xyzi(i,j+7,k+1)-2.992511007649597D3*xyzi(i,j+5,k+1)+5.320
     $        01956915484D2*xyzi(i,j+3,k+1))+angi
            angi=yrk(110)*(-3.952135527704191D3*xyzi(i+8,j+1,k+1)-1.5808
     $        54211081677D4*xyzi(i+6,j+3,k+1)+8.320285321482508D3*xyzi(i
     $        +6,j+1,k+1)-2.371281316622515D4*xyzi(i+4,j+5,k+1)+2.496085
     $        596444752D4*xyzi(i+4,j+3,k+1)-5.873142579870006D3*xyzi(i+4
     $        ,j+1,k+1)-1.580854211081677D4*xyzi(i+2,j+7,k+1)+2.49608559
     $        6444752D4*xyzi(i+2,j+5,k+1)-1.174628515974001D4*xyzi(i+2,j
     $        +3,k+1)+1.566171354632002D3*xyzi(i+2,j+1,k+1)-3.9521355277
     $        04191D3*xyzi(i,j+9,k+1)+8.320285321482508D3*xyzi(i,j+7,k+1
     $        )-5.873142579870006D3*xyzi(i,j+5,k+1)+1.566171354632002D3*
     $        xyzi(i,j+3,k+1)-1.20474719587077D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA(10)=angi
        endif
c...
c...  l=12
c...
        if(lmax.ge.12)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(157)*(1.170163993078432D4*xyzi(i,j,k+12)-3.35786189
     $        3181587D4*xyzi(i,j,k+10)+3.597709171265986D4*xyzi(i,j,k+8)
     $        -1.767295733253467D4*xyzi(i,j,k+6)+3.898446470412059D3*xyz
     $        i(i,j,k+4)-3.118757176329647D2*xyzi(i,j,k+2)+3.99840663632
     $        0061D0*xyzi(i,j,k))+angi
            angi=yrk(159)*(-1.409330983563191D4*xyzi(i+12,j,k)-5.6373239
     $        34252763D4*xyzi(i+10,j+2,k)+4.289268210844493D4*xyzi(i+10,
     $        j,k)-7.046654917815953D4*xyzi(i+8,j+4,k)+1.286780463253348
     $        D5*xyzi(i+8,j+2,k)-4.902020812393707D4*xyzi(i+8,j,k)+8.578
     $        536421688987D4*xyzi(i+6,j+4,k)-9.804041624787413D4*xyzi(i+
     $        6,j+2,k)+2.580010953891425D4*xyzi(i+6,j,k)+7.0466549178159
     $        53D4*xyzi(i+4,j+8,k)-8.578536421688987D4*xyzi(i+4,j+6,k)+2
     $        .580010953891425D4*xyzi(i+4,j+2,k)-6.070614009156293D3*xyz
     $        i(i+4,j,k)+5.637323934252763D4*xyzi(i+2,j+10,k)-1.28678046
     $        3253348D5*xyzi(i+2,j+8,k)+9.804041624787413D4*xyzi(i+2,j+6
     $        ,k)-2.580010953891425D4*xyzi(i+2,j+4,k)+4.856491207325035D
     $        2*xyzi(i+2,j,k)+1.409330983563191D4*xyzi(i,j+12,k)-4.28926
     $        8210844493D4*xyzi(i,j+10,k)+4.902020812393707D4*xyzi(i,j+8
     $        ,k)-2.580010953891425D4*xyzi(i,j+6,k)+6.070614009156293D3*
     $        xyzi(i,j+4,k)-4.856491207325035D2*xyzi(i,j+2,k))+angi
            angi=yrk(161)*(8.630354471061409D3*xyzi(i+12,j,k)-1.72607089
     $        4212282D4*xyzi(i+10,j+2,k)-2.401489939773609D4*xyzi(i+10,j
     $        ,k)-1.467160260080439D5*xyzi(i+8,j+4,k)+7.204469819320828D
     $        4*xyzi(i+8,j+2,k)+2.401489939773609D4*xyzi(i+8,j,k)-2.4164
     $        99251897194D5*xyzi(i+6,j+6,k)+3.362085915683053D5*xyzi(i+6
     $        ,j+4,k)-9.605959759094437D4*xyzi(i+6,j+2,k)-1.011153658852
     $        046D4*xyzi(i+6,j,k)-1.467160260080439D5*xyzi(i+4,j+8,k)+3.
     $        362085915683053D5*xyzi(i+4,j+6,k)-2.401489939773609D5*xyzi
     $        (i+4,j+4,k)+5.05576829426023D4*xyzi(i+4,j+2,k)+1.486990674
     $        782421D3*xyzi(i+4,j,k)-1.726070894212282D4*xyzi(i+2,j+10,k
     $        )+7.204469819320828D4*xyzi(i+2,j+8,k)-9.605959759094437D4*
     $        xyzi(i+2,j+6,k)+5.05576829426023D4*xyzi(i+2,j+4,k)-8.92194
     $        4048694524D3*xyzi(i+2,j+2,k)+8.630354471061409D3*xyzi(i,j+
     $        12,k)-2.401489939773609D4*xyzi(i,j+10,k)+2.401489939773609
     $        D4*xyzi(i,j+8,k)-1.011153658852046D4*xyzi(i,j+6,k)+1.48699
     $        0674782421D3*xyzi(i,j+4,k))+angi
            angi=yrk(163)*(-3.692002053806592D3*xyzi(i+12,j,k)+4.4304024
     $        6456791D4*xyzi(i+10,j+2,k)+8.668178735024172D3*xyzi(i+10,j
     $        ,k)+9.968405545277798D4*xyzi(i+8,j+4,k)-1.126863235553142D
     $        5*xyzi(i+8,j+2,k)-6.604326655256512D3*xyzi(i+8,j,k)-1.2135
     $        45022903384D5*xyzi(i+6,j+4,k)+9.246057317359117D4*xyzi(i+6
     $        ,j+2,k)+1.622115318834933D3*xyzi(i+6,j,k)-9.96840554527779
     $        8D4*xyzi(i+4,j+8,k)+1.213545022903384D5*xyzi(i+4,j+6,k)-2.
     $        433172978252399D4*xyzi(i+4,j+2,k)-4.43040246456791D4*xyzi(
     $        i+2,j+10,k)+1.126863235553142D5*xyzi(i+2,j+8,k)-9.24605731
     $        7359117D4*xyzi(i+2,j+6,k)+2.433172978252399D4*xyzi(i+2,j+4
     $        ,k)+3.692002053806592D3*xyzi(i,j+12,k)-8.668178735024172D3
     $        *xyzi(i,j+10,k)+6.604326655256512D3*xyzi(i,j+8,k)-1.622115
     $        318834933D3*xyzi(i,j+6,k))+angi
            angi=yrk(165)*(1.037363021977718D3*xyzi(i+12,j,k)-2.69714385
     $        7142068D4*xyzi(i+10,j+2,k)-1.80410960343951D3*xyzi(i+10,j,
     $        k)+1.556044532966577D4*xyzi(i+8,j+4,k)+4.871095929286677D4
     $        *xyzi(i+8,j+2,k)+7.731898300455043D2*xyzi(i+8,j,k)+8.71384
     $        9384612834D4*xyzi(i+6,j+6,k)-7.577260334445942D4*xyzi(i+6,
     $        j+4,k)-2.164931524127412D4*xyzi(i+6,j+2,k)+1.5560445329665
     $        77D4*xyzi(i+4,j+8,k)-7.577260334445942D4*xyzi(i+4,j+6,k)+5
     $        .41232881031853D4*xyzi(i+4,j+4,k)-2.697143857142068D4*xyzi
     $        (i+2,j+10,k)+4.871095929286677D4*xyzi(i+2,j+8,k)-2.1649315
     $        24127412D4*xyzi(i+2,j+6,k)+1.037363021977718D3*xyzi(i,j+12
     $        ,k)-1.80410960343951D3*xyzi(i,j+10,k)+7.731898300455043D2*
     $        xyzi(i,j+8,k))+angi
            angi=yrk(167)*(-1.671861890280821D2*xyzi(i+12,j,k)+7.3561923
     $        17235614D3*xyzi(i+10,j+2,k)+1.599172242877307D2*xyzi(i+10,
     $        j,k)-2.758572118963355D4*xyzi(i+8,j+4,k)-7.196275092947883
     $        D3*xyzi(i+8,j+2,k)+3.358261710042346D4*xyzi(i+6,j+4,k)+2.7
     $        58572118963355D4*xyzi(i+4,j+8,k)-3.358261710042346D4*xyzi(
     $        i+4,j+6,k)-7.356192317235614D3*xyzi(i+2,j+10,k)+7.19627509
     $        2947883D3*xyzi(i+2,j+8,k)+1.671861890280821D2*xyzi(i,j+12,
     $        k)-1.599172242877307D2*xyzi(i,j+10,k))+angi
            angi=yrk(169)*(1.006342599515217D1*xyzi(i+12,j,k)-6.64186115
     $        6800431D2*xyzi(i+10,j+2,k)+4.981395867600323D3*xyzi(i+8,j+
     $        4,k)-9.298605619520603D3*xyzi(i+6,j+6,k)+4.981395867600323
     $        D3*xyzi(i+4,j+8,k)-6.641861156800431D2*xyzi(i+2,j+10,k)+1.
     $        006342599515217D1*xyzi(i,j+12,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(145)*(1.20761111941826D2*xyzi(i+11,j+1,k)-2.2139537
     $        18933477D3*xyzi(i+9,j+3,k)+7.970233388160517D3*xyzi(i+7,j+
     $        5,k)-7.970233388160517D3*xyzi(i+5,j+7,k)+2.213953718933477
     $        D3*xyzi(i+3,j+9,k)-1.20761111941826D2*xyzi(i+1,j+11,k))+an
     $        gi
            angi=yrk(147)*(-1.671861890280821D3*xyzi(i+11,j+1,k)+1.83904
     $        8079308904D4*xyzi(i+9,j+3,k)+1.599172242877307D3*xyzi(i+9,
     $        j+1,k)-2.206857695170684D4*xyzi(i+7,j+5,k)-1.9190066914527
     $        69D4*xyzi(i+7,j+3,k)-2.206857695170684D4*xyzi(i+5,j+7,k)+4
     $        .029914052050815D4*xyzi(i+5,j+5,k)+1.839048079308904D4*xyz
     $        i(i+3,j+9,k)-1.919006691452769D4*xyzi(i+3,j+7,k)-1.6718618
     $        90280821D3*xyzi(i+1,j+11,k)+1.599172242877307D3*xyzi(i+1,j
     $        +9,k))+angi
            angi=yrk(149)*(8.298904175821746D3*xyzi(i+11,j+1,k)-4.149452
     $        087910873D4*xyzi(i+9,j+3,k)-1.443287682751608D4*xyzi(i+9,j
     $        +1,k)-4.979342505493048D4*xyzi(i+7,j+5,k)+8.65972609650964
     $        8D4*xyzi(i+7,j+3,k)+6.185518640364035D3*xyzi(i+7,j+1,k)+4.
     $        979342505493048D4*xyzi(i+5,j+7,k)-4.329863048254824D4*xyzi
     $        (i+5,j+3,k)+4.149452087910873D4*xyzi(i+3,j+9,k)-8.65972609
     $        6509648D4*xyzi(i+3,j+7,k)+4.329863048254824D4*xyzi(i+3,j+5
     $        ,k)-8.298904175821746D3*xyzi(i+1,j+11,k)+1.443287682751608
     $        D4*xyzi(i+1,j+9,k)-6.185518640364035D3*xyzi(i+1,j+7,k))+an
     $        gi
            angi=yrk(151)*(-2.215201232283955D4*xyzi(i+11,j+1,k)+7.38400
     $        4107613184D3*xyzi(i+9,j+3,k)+5.200907241014503D4*xyzi(i+9,
     $        j+1,k)+1.329120739370373D5*xyzi(i+7,j+5,k)-6.9345429880193
     $        38D4*xyzi(i+7,j+3,k)-3.962595993153907D4*xyzi(i+7,j+1,k)+1
     $        .329120739370373D5*xyzi(i+5,j+7,k)-2.427090045806768D5*xyz
     $        i(i+5,j+5,k)+9.246057317359117D4*xyzi(i+5,j+3,k)+9.7326919
     $        13009597D3*xyzi(i+5,j+1,k)+7.384004107613184D3*xyzi(i+3,j+
     $        9,k)-6.934542988019338D4*xyzi(i+3,j+7,k)+9.246057317359117
     $        D4*xyzi(i+3,j+5,k)-3.244230637669866D4*xyzi(i+3,j+3,k)-2.2
     $        15201232283955D4*xyzi(i+1,j+11,k)+5.200907241014503D4*xyzi
     $        (i+1,j+9,k)-3.962595993153907D4*xyzi(i+1,j+7,k)+9.73269191
     $        3009597D3*xyzi(i+1,j+5,k))+angi
            angi=yrk(153)*(3.452141788424563D4*xyzi(i+11,j+1,k)+1.035642
     $        536527369D5*xyzi(i+9,j+3,k)-9.605959759094437D4*xyzi(i+9,j
     $        +1,k)+6.904283576849127D4*xyzi(i+7,j+5,k)-1.92119195181888
     $        7D5*xyzi(i+7,j+3,k)+9.605959759094437D4*xyzi(i+7,j+1,k)-6.
     $        904283576849127D4*xyzi(i+5,j+7,k)+9.605959759094437D4*xyzi
     $        (i+5,j+3,k)-4.044614635408184D4*xyzi(i+5,j+1,k)-1.03564253
     $        6527369D5*xyzi(i+3,j+9,k)+1.921191951818887D5*xyzi(i+3,j+7
     $        ,k)-9.605959759094437D4*xyzi(i+3,j+5,k)+5.947962699129683D
     $        3*xyzi(i+3,j+1,k)-3.452141788424563D4*xyzi(i+1,j+11,k)+9.6
     $        05959759094437D4*xyzi(i+1,j+9,k)-9.605959759094437D4*xyzi(
     $        i+1,j+7,k)+4.044614635408184D4*xyzi(i+1,j+5,k)-5.947962699
     $        129683D3*xyzi(i+1,j+3,k))+angi
            angi=yrk(155)*(-2.818661967126381D4*xyzi(i+11,j+1,k)-1.40933
     $        0983563191D5*xyzi(i+9,j+3,k)+8.578536421688987D4*xyzi(i+9,
     $        j+1,k)-2.818661967126381D5*xyzi(i+7,j+5,k)+3.4314145686755
     $        95D5*xyzi(i+7,j+3,k)-9.804041624787413D4*xyzi(i+7,j+1,k)-2
     $        .818661967126381D5*xyzi(i+5,j+7,k)+5.147121853013392D5*xyz
     $        i(i+5,j+5,k)-2.941212487436224D5*xyzi(i+5,j+3,k)+5.1600219
     $        07782849D4*xyzi(i+5,j+1,k)-1.409330983563191D5*xyzi(i+3,j+
     $        9,k)+3.431414568675595D5*xyzi(i+3,j+7,k)-2.941212487436224
     $        D5*xyzi(i+3,j+5,k)+1.03200438155657D5*xyzi(i+3,j+3,k)-1.21
     $        4122801831259D4*xyzi(i+3,j+1,k)-2.818661967126381D4*xyzi(i
     $        +1,j+11,k)+8.578536421688987D4*xyzi(i+1,j+9,k)-9.804041624
     $        787413D4*xyzi(i+1,j+7,k)+5.160021907782849D4*xyzi(i+1,j+5,
     $        k)-1.214122801831259D4*xyzi(i+1,j+3,k)+9.712982414650069D2
     $        *xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(158)*(1.589939778654773D4*xyzi(i+11,j,k+1)+7.949698
     $        893273864D4*xyzi(i+9,j+2,k+1)-4.147668987795059D4*xyzi(i+9
     $        ,j,k+1)+1.589939778654773D5*xyzi(i+7,j+4,k+1)-1.6590675951
     $        18024D5*xyzi(i+7,j+2,k+1)+3.9501609407572D4*xyzi(i+7,j,k+1
     $        )+1.589939778654773D5*xyzi(i+5,j+6,k+1)-2.488601392677036D
     $        5*xyzi(i+5,j+4,k+1)+1.18504828222716D5*xyzi(i+5,j+2,k+1)-1
     $        .663225659266189D4*xyzi(i+5,j,k+1)+7.949698893273864D4*xyz
     $        i(i+3,j+8,k+1)-1.659067595118024D5*xyzi(i+3,j+6,k+1)+1.185
     $        04828222716D5*xyzi(i+3,j+4,k+1)-3.326451318532379D4*xyzi(i
     $        +3,j+2,k+1)+2.935104104587393D3*xyzi(i+3,j,k+1)+1.58993977
     $        8654773D4*xyzi(i+1,j+10,k+1)-4.147668987795059D4*xyzi(i+1,
     $        j+8,k+1)+3.9501609407572D4*xyzi(i+1,j+6,k+1)-1.66322565926
     $        6189D4*xyzi(i+1,j+4,k+1)+2.935104104587393D3*xyzi(i+1,j+2,
     $        k+1)-1.565388855779943D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(160)*(-1.150713929474854D4*xyzi(i+11,j,k+1)-1.15071
     $        3929474854D4*xyzi(i+9,j+2,k+1)+2.801738263069211D4*xyzi(i+
     $        9,j,k+1)+6.904283576849127D4*xyzi(i+7,j+4,k+1)-2.401489939
     $        773609D4*xyzi(i+7,j,k+1)+1.610999501264796D5*xyzi(i+5,j+6,
     $        k+1)-1.681042957841527D5*xyzi(i+5,j+4,k+1)+2.4014899397736
     $        09D4*xyzi(i+5,j+2,k+1)+8.426280490433717D3*xyzi(i+5,j,k+1)
     $        +1.26578532242234D5*xyzi(i+3,j+8,k+1)-2.241390610455369D5*
     $        xyzi(i+3,j+6,k+1)+1.200744969886805D5*xyzi(i+3,j+4,k+1)-1.
     $        685256098086743D4*xyzi(i+3,j+2,k+1)-9.913271165216138D2*xy
     $        zi(i+3,j,k+1)+3.452141788424563D4*xyzi(i+1,j+10,k+1)-8.405
     $        214789207633D4*xyzi(i+1,j+8,k+1)+7.204469819320828D4*xyzi(
     $        i+1,j+6,k+1)-2.527884147130115D4*xyzi(i+1,j+4,k+1)+2.97398
     $        1349564841D3*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(162)*(5.920374324261427D3*xyzi(i+11,j,k+1)-4.144262
     $        026982999D4*xyzi(i+9,j+2,k+1)-1.235556380715428D4*xyzi(i+9
     $        ,j,k+1)-1.302482351337514D5*xyzi(i+7,j+4,k+1)+9.8844510457
     $        23426D4*xyzi(i+7,j+2,k+1)+8.237042538102855D3*xyzi(i+7,j,k
     $        +1)-8.288524053965998D4*xyzi(i+5,j+6,k+1)+1.7297789330016D
     $        5*xyzi(i+5,j+4,k+1)-7.41333828429257D4*xyzi(i+5,j+2,k+1)-1
     $        .73411421854797D3*xyzi(i+5,j,k+1)+2.960187162130714D4*xyzi
     $        (i+3,j+8,k+1)-4.118521269051428D4*xyzi(i+3,j+4,k+1)+1.7341
     $        14218547969D4*xyzi(i+3,j+2,k+1)+2.960187162130714D4*xyzi(i
     $        +1,j+10,k+1)-6.177781903577141D4*xyzi(i+1,j+8,k+1)+4.11852
     $        1269051428D4*xyzi(i+1,j+6,k+1)-8.670571092739847D3*xyzi(i+
     $        1,j+4,k+1))+angi
            angi=yrk(164)*(-2.074726043955437D3*xyzi(i+11,j,k+1)+3.94197
     $        948351533D4*xyzi(i+9,j+2,k+1)+3.247397286191118D3*xyzi(i+9
     $        ,j,k+1)+1.244835626373262D4*xyzi(i+7,j+4,k+1)-6.4947945723
     $        82236D4*xyzi(i+7,j+2,k+1)-1.237103728072807D3*xyzi(i+7,j,k
     $        +1)-8.713849384612834D4*xyzi(i+5,j+6,k+1)+4.54635620066756
     $        5D4*xyzi(i+5,j+4,k+1)+2.597917828952894D4*xyzi(i+5,j+2,k+1
     $        )-4.356924692306417D4*xyzi(i+3,j+8,k+1)+9.092712401335131D
     $        4*xyzi(i+3,j+6,k+1)-4.329863048254824D4*xyzi(i+3,j+4,k+1)+
     $        1.452308230768806D4*xyzi(i+1,j+10,k+1)-2.273178100333783D4
     $        *xyzi(i+1,j+8,k+1)+8.659726096509648D3*xyzi(i+1,j+6,k+1))+
     $        angi
            angi=yrk(166)*(4.527423401296222D2*xyzi(i+11,j,k+1)-1.584598
     $        190453678D4*xyzi(i+9,j+2,k+1)-3.936889914170628D2*xyzi(i+9
     $        ,j,k+1)+4.0746810611666D4*xyzi(i+7,j+4,k+1)+1.417280369101
     $        426D4*xyzi(i+7,j+2,k+1)+1.901517828544413D4*xyzi(i+5,j+6,k
     $        +1)-4.960481291854991D4*xyzi(i+5,j+4,k+1)-3.39556755097216
     $        6D4*xyzi(i+3,j+8,k+1)+3.306987527903327D4*xyzi(i+3,j+6,k+1
     $        )+4.0746810611666D3*xyzi(i+1,j+10,k+1)-3.543200922753565D3
     $        *xyzi(i+1,j+8,k+1))+angi
            angi=yrk(168)*(-4.930051750476566D1*xyzi(i+11,j,k+1)+2.71152
     $        8462762111D3*xyzi(i+9,j+2,k+1)-1.626917077657267D4*xyzi(i+
     $        7,j+4,k+1)+2.277683908720174D4*xyzi(i+5,j+6,k+1)-8.1345853
     $        88286334D3*xyzi(i+3,j+8,k+1)+5.423056925524223D2*xyzi(i+1,
     $        j+10,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(146)*(-5.423056925524223D2*xyzi(i+10,j+1,k+1)+8.134
     $        585388286334D3*xyzi(i+8,j+3,k+1)-2.277683908720174D4*xyzi(
     $        i+6,j+5,k+1)+1.626917077657267D4*xyzi(i+4,j+7,k+1)-2.71152
     $        8462762111D3*xyzi(i+2,j+9,k+1)+4.930051750476566D1*xyzi(i,
     $        j+11,k+1))+angi
            angi=yrk(148)*(4.0746810611666D3*xyzi(i+10,j+1,k+1)-3.395567
     $        550972166D4*xyzi(i+8,j+3,k+1)-3.543200922753565D3*xyzi(i+8
     $        ,j+1,k+1)+1.901517828544413D4*xyzi(i+6,j+5,k+1)+3.30698752
     $        7903327D4*xyzi(i+6,j+3,k+1)+4.0746810611666D4*xyzi(i+4,j+7
     $        ,k+1)-4.960481291854991D4*xyzi(i+4,j+5,k+1)-1.584598190453
     $        678D4*xyzi(i+2,j+9,k+1)+1.417280369101426D4*xyzi(i+2,j+7,k
     $        +1)+4.527423401296222D2*xyzi(i,j+11,k+1)-3.936889914170628
     $        D2*xyzi(i,j+9,k+1))+angi
            angi=yrk(150)*(-1.452308230768806D4*xyzi(i+10,j+1,k+1)+4.356
     $        924692306417D4*xyzi(i+8,j+3,k+1)+2.273178100333783D4*xyzi(
     $        i+8,j+1,k+1)+8.713849384612834D4*xyzi(i+6,j+5,k+1)-9.09271
     $        2401335131D4*xyzi(i+6,j+3,k+1)-8.659726096509648D3*xyzi(i+
     $        6,j+1,k+1)-1.244835626373262D4*xyzi(i+4,j+7,k+1)-4.5463562
     $        00667565D4*xyzi(i+4,j+5,k+1)+4.329863048254824D4*xyzi(i+4,
     $        j+3,k+1)-3.94197948351533D4*xyzi(i+2,j+9,k+1)+6.4947945723
     $        82236D4*xyzi(i+2,j+7,k+1)-2.597917828952894D4*xyzi(i+2,j+5
     $        ,k+1)+2.074726043955437D3*xyzi(i,j+11,k+1)-3.2473972861911
     $        18D3*xyzi(i,j+9,k+1)+1.237103728072807D3*xyzi(i,j+7,k+1))+
     $        angi
            angi=yrk(152)*(2.960187162130714D4*xyzi(i+10,j+1,k+1)+2.9601
     $        87162130714D4*xyzi(i+8,j+3,k+1)-6.177781903577141D4*xyzi(i
     $        +8,j+1,k+1)-8.288524053965998D4*xyzi(i+6,j+5,k+1)+4.118521
     $        269051428D4*xyzi(i+6,j+1,k+1)-1.302482351337514D5*xyzi(i+4
     $        ,j+7,k+1)+1.7297789330016D5*xyzi(i+4,j+5,k+1)-4.1185212690
     $        51428D4*xyzi(i+4,j+3,k+1)-8.670571092739847D3*xyzi(i+4,j+1
     $        ,k+1)-4.144262026982999D4*xyzi(i+2,j+9,k+1)+9.884451045723
     $        426D4*xyzi(i+2,j+7,k+1)-7.41333828429257D4*xyzi(i+2,j+5,k+
     $        1)+1.734114218547969D4*xyzi(i+2,j+3,k+1)+5.920374324261427
     $        D3*xyzi(i,j+11,k+1)-1.235556380715428D4*xyzi(i,j+9,k+1)+8.
     $        237042538102855D3*xyzi(i,j+7,k+1)-1.73411421854797D3*xyzi(
     $        i,j+5,k+1))+angi
            angi=yrk(154)*(-3.452141788424563D4*xyzi(i+10,j+1,k+1)-1.265
     $        78532242234D5*xyzi(i+8,j+3,k+1)+8.405214789207633D4*xyzi(i
     $        +8,j+1,k+1)-1.610999501264796D5*xyzi(i+6,j+5,k+1)+2.241390
     $        610455369D5*xyzi(i+6,j+3,k+1)-7.204469819320828D4*xyzi(i+6
     $        ,j+1,k+1)-6.904283576849127D4*xyzi(i+4,j+7,k+1)+1.68104295
     $        7841527D5*xyzi(i+4,j+5,k+1)-1.200744969886805D5*xyzi(i+4,j
     $        +3,k+1)+2.527884147130115D4*xyzi(i+4,j+1,k+1)+1.1507139294
     $        74854D4*xyzi(i+2,j+9,k+1)-2.401489939773609D4*xyzi(i+2,j+5
     $        ,k+1)+1.685256098086743D4*xyzi(i+2,j+3,k+1)-2.973981349564
     $        841D3*xyzi(i+2,j+1,k+1)+1.150713929474854D4*xyzi(i,j+11,k+
     $        1)-2.801738263069211D4*xyzi(i,j+9,k+1)+2.401489939773609D4
     $        *xyzi(i,j+7,k+1)-8.426280490433717D3*xyzi(i,j+5,k+1)+9.913
     $        271165216138D2*xyzi(i,j+3,k+1))+angi
            angi=yrk(156)*(1.589939778654773D4*xyzi(i+10,j+1,k+1)+7.9496
     $        98893273864D4*xyzi(i+8,j+3,k+1)-4.147668987795059D4*xyzi(i
     $        +8,j+1,k+1)+1.589939778654773D5*xyzi(i+6,j+5,k+1)-1.659067
     $        595118024D5*xyzi(i+6,j+3,k+1)+3.9501609407572D4*xyzi(i+6,j
     $        +1,k+1)+1.589939778654773D5*xyzi(i+4,j+7,k+1)-2.4886013926
     $        77036D5*xyzi(i+4,j+5,k+1)+1.18504828222716D5*xyzi(i+4,j+3,
     $        k+1)-1.663225659266189D4*xyzi(i+4,j+1,k+1)+7.9496988932738
     $        64D4*xyzi(i+2,j+9,k+1)-1.659067595118024D5*xyzi(i+2,j+7,k+
     $        1)+1.18504828222716D5*xyzi(i+2,j+5,k+1)-3.326451318532379D
     $        4*xyzi(i+2,j+3,k+1)+2.935104104587393D3*xyzi(i+2,j+1,k+1)+
     $        1.589939778654773D4*xyzi(i,j+11,k+1)-4.147668987795059D4*x
     $        yzi(i,j+9,k+1)+3.9501609407572D4*xyzi(i,j+7,k+1)-1.6632256
     $        59266189D4*xyzi(i,j+5,k+1)+2.935104104587393D3*xyzi(i,j+3,
     $        k+1)-1.565388855779943D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA(12)=angi
        endif
c...
c...  l=14
c...
        if(lmax.ge.14)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(211)*(4.674208812108592D4*xyzi(i,j,k+14)-1.57538148
     $        8525489D5*xyzi(i,j,k+12)+2.079503564853645D5*xyzi(i,j,k+10
     $        )-1.356197977078464D5*xyzi(i,j,k+8)+4.52065992359488D4*xyz
     $        i(i,j,k+6)-7.137884089886653D3*xyzi(i,j,k+4)+4.19875534699
     $        2149D2*xyzi(i,j,k+2)-3.998814616182999D0*xyzi(i,j,k))+angi
            angi=yrk(213)*(5.756429376136187D4*xyzi(i+14,j,k)+2.87821468
     $        8068094D5*xyzi(i+12,j+2,k)-2.046730444848422D5*xyzi(i+12,j
     $        ,k)+5.180786438522568D5*xyzi(i+10,j+4,k)-8.186921779393688
     $        D5*xyzi(i+10,j+2,k)+2.865422622787791D5*xyzi(i+10,j,k)+2.8
     $        78214688068094D5*xyzi(i+8,j+6,k)-1.023365222424211D6*xyzi(
     $        i+8,j+4,k)+8.596267868363373D5*xyzi(i+8,j+2,k)-1.993337476
     $        721942D5*xyzi(i+8,j,k)-2.878214688068094D5*xyzi(i+6,j+8,k)
     $        +5.730845245575582D5*xyzi(i+6,j+4,k)-3.986674953443883D5*x
     $        yzi(i+6,j+2,k)+7.119062416864077D4*xyzi(i+6,j,k)-5.1807864
     $        38522568D5*xyzi(i+4,j+10,k)+1.023365222424211D6*xyzi(i+4,j
     $        +8,k)-5.730845245575582D5*xyzi(i+4,j+6,k)+7.11906241686407
     $        7D4*xyzi(i+4,j+2,k)-1.19899998599816D4*xyzi(i+4,j,k)-2.878
     $        214688068094D5*xyzi(i+2,j+12,k)+8.186921779393688D5*xyzi(i
     $        +2,j+10,k)-8.596267868363373D5*xyzi(i+2,j+8,k)+3.986674953
     $        443883D5*xyzi(i+2,j+6,k)-7.119062416864077D4*xyzi(i+2,j+4,
     $        k)+7.052941094106825D2*xyzi(i+2,j,k)-5.756429376136187D4*x
     $        yzi(i,j+14,k)+2.046730444848422D5*xyzi(i,j+12,k)-2.8654226
     $        22787791D5*xyzi(i,j+10,k)+1.993337476721942D5*xyzi(i,j+8,k
     $        )-7.119062416864077D4*xyzi(i,j+6,k)+1.19899998599816D4*xyz
     $        i(i,j+4,k)-7.052941094106825D2*xyzi(i,j+2,k))+angi
            angi=yrk(215)*(-3.780762817453435D4*xyzi(i+14,j,k)+3.7807628
     $        17453435D4*xyzi(i+12,j+2,k)+1.260254272484478D5*xyzi(i+12,
     $        j,k)+7.183449353161527D5*xyzi(i+10,j+4,k)-2.52050854496895
     $        7D5*xyzi(i+10,j+2,k)-1.613125468780132D5*xyzi(i+10,j,k)+1.
     $        701343267854046D6*xyzi(i+8,j+6,k)-2.142432263223613D6*xyzi
     $        (i+8,j+4,k)+4.839376406340397D5*xyzi(i+8,j+2,k)+9.81902459
     $        2574719D4*xyzi(i+8,j,k)+1.701343267854046D6*xyzi(i+6,j+8,k
     $        )-3.52871196295654D6*xyzi(i+6,j+6,k)+2.258375656292185D6*x
     $        yzi(i+6,j+4,k)-3.927609837029888D5*xyzi(i+6,j+2,k)-2.80543
     $        5597878491D4*xyzi(i+6,j,k)+7.183449353161527D5*xyzi(i+4,j+
     $        10,k)-2.142432263223613D6*xyzi(i+4,j+8,k)+2.25837565629218
     $        5D6*xyzi(i+4,j+6,k)-9.819024592574719D5*xyzi(i+4,j+4,k)+1.
     $        402717798939246D5*xyzi(i+4,j+2,k)+2.953090103029991D3*xyzi
     $        (i+4,j,k)+3.780762817453435D4*xyzi(i+2,j+12,k)-2.520508544
     $        968957D5*xyzi(i+2,j+10,k)+4.839376406340397D5*xyzi(i+2,j+8
     $        ,k)-3.927609837029888D5*xyzi(i+2,j+6,k)+1.402717798939246D
     $        5*xyzi(i+2,j+4,k)-1.771854061817994D4*xyzi(i+2,j+2,k)-3.78
     $        0762817453435D4*xyzi(i,j+14,k)+1.260254272484478D5*xyzi(i,
     $        j+12,k)-1.613125468780132D5*xyzi(i,j+10,k)+9.8190245925747
     $        19D4*xyzi(i,j+8,k)-2.805435597878491D4*xyzi(i,j+6,k)+2.953
     $        090103029991D3*xyzi(i,j+4,k))+angi
            angi=yrk(217)*(1.839962151616926D4*xyzi(i+14,j,k)-2.02395836
     $        6778619D5*xyzi(i+12,j+2,k)-5.451739708494596D4*xyzi(i+12,j
     $        ,k)-7.175852391306012D5*xyzi(i+10,j+4,k)+6.542087650193515
     $        D5*xyzi(i+10,j+2,k)+5.887878885174163D4*xyzi(i+10,j,k)-4.9
     $        678978093657D5*xyzi(i+8,j+6,k)+1.471969721293541D6*xyzi(i+
     $        8,j+4,k)-7.654242550726412D5*xyzi(i+8,j+2,k)-2.73061049747
     $        2076D4*xyzi(i+8,j,k)+4.9678978093657D5*xyzi(i+6,j+8,k)-8.2
     $        43030439243829D5*xyzi(i+6,j+4,k)+3.822854696460906D5*xyzi(
     $        i+6,j+2,k)+4.551017495786793D3*xyzi(i+6,j,k)+7.17585239130
     $        6012D5*xyzi(i+4,j+10,k)-1.471969721293541D6*xyzi(i+4,j+8,k
     $        )+8.243030439243829D5*xyzi(i+4,j+6,k)-6.826526243680189D4*
     $        xyzi(i+4,j+2,k)+2.023958366778619D5*xyzi(i+2,j+12,k)-6.542
     $        087650193515D5*xyzi(i+2,j+10,k)+7.654242550726412D5*xyzi(i
     $        +2,j+8,k)-3.822854696460906D5*xyzi(i+2,j+6,k)+6.8265262436
     $        80189D4*xyzi(i+2,j+4,k)-1.839962151616926D4*xyzi(i,j+14,k)
     $        +5.451739708494596D4*xyzi(i,j+12,k)-5.887878885174163D4*xy
     $        zi(i,j+10,k)+2.730610497472076D4*xyzi(i,j+8,k)-4.551017495
     $        786793D3*xyzi(i,j+6,k))+angi
            angi=yrk(219)*(-6.405925968013536D3*xyzi(i+14,j,k)+1.6014814
     $        92003384D5*xyzi(i+12,j+2,k)+1.565893014403309D4*xyzi(i+12,
     $        j,k)+7.04651856481489D4*xyzi(i+10,j+4,k)-4.071321837448603
     $        D5*xyzi(i+10,j+2,k)-1.252714411522647D4*xyzi(i+10,j,k)-6.3
     $        41866708333401D5*xyzi(i+8,j+6,k)+2.348839521604963D5*xyzi(
     $        i+8,j+4,k)+3.382328911111147D5*xyzi(i+8,j+2,k)+3.267950638
     $        754731D3*xyzi(i+8,j,k)-6.341866708333401D5*xyzi(i+6,j+8,k)
     $        +1.315350132098779D6*xyzi(i+6,j+6,k)-5.261400528395117D5*x
     $        yzi(i+6,j+4,k)-9.150261788513248D4*xyzi(i+6,j+2,k)+7.04651
     $        856481489D4*xyzi(i+4,j+10,k)+2.348839521604963D5*xyzi(i+4,
     $        j+8,k)-5.261400528395117D5*xyzi(i+4,j+6,k)+2.2875654471283
     $        12D5*xyzi(i+4,j+4,k)+1.601481492003384D5*xyzi(i+2,j+12,k)-
     $        4.071321837448603D5*xyzi(i+2,j+10,k)+3.382328911111147D5*x
     $        yzi(i+2,j+8,k)-9.150261788513248D4*xyzi(i+2,j+6,k)-6.40592
     $        5968013536D3*xyzi(i,j+14,k)+1.565893014403309D4*xyzi(i,j+1
     $        2,k)-1.252714411522647D4*xyzi(i,j+10,k)+3.267950638754731D
     $        3*xyzi(i,j+8,k))+angi
            angi=yrk(221)*(1.493389191601027D3*xyzi(i+14,j,k)-6.42157352
     $        3884417D4*xyzi(i+12,j+2,k)-2.654914118401826D3*xyzi(i+12,j
     $        ,k)+1.807000921837243D5*xyzi(i+10,j+4,k)+1.168162212096803
     $        D5*xyzi(i+10,j+2,k)+1.168162212096803D3*xyzi(i+10,j,k)+2.4
     $        64092166141695D5*xyzi(i+8,j+6,k)-4.380608295363013D5*xyzi(
     $        i+8,j+4,k)-5.256729954435616D4*xyzi(i+8,j+2,k)-2.464092166
     $        141695D5*xyzi(i+6,j+8,k)+2.453140645403287D5*xyzi(i+6,j+4,
     $        k)-1.807000921837243D5*xyzi(i+4,j+10,k)+4.380608295363013D
     $        5*xyzi(i+4,j+8,k)-2.453140645403287D5*xyzi(i+4,j+6,k)+6.42
     $        1573523884417D4*xyzi(i+2,j+12,k)-1.168162212096803D5*xyzi(
     $        i+2,j+10,k)+5.256729954435616D4*xyzi(i+2,j+8,k)-1.49338919
     $        1601027D3*xyzi(i,j+14,k)+2.654914118401826D3*xyzi(i,j+12,k
     $        )-1.168162212096803D3*xyzi(i,j+10,k))+angi
            angi=yrk(223)*(-2.029116341627528D2*xyzi(i+14,j,k)+1.3189256
     $        22057893D4*xyzi(i+12,j+2,k)+1.953963884530212D2*xyzi(i+12,
     $        j,k)-8.704909105582095D4*xyzi(i+10,j+4,k)-1.28961616378994
     $        D4*xyzi(i+10,j+2,k)+8.704909105582095D4*xyzi(i+8,j+6,k)+9.
     $        672121228424549D4*xyzi(i+8,j+4,k)+8.704909105582095D4*xyzi
     $        (i+6,j+8,k)-1.805462629305916D5*xyzi(i+6,j+6,k)-8.70490910
     $        5582095D4*xyzi(i+4,j+10,k)+9.672121228424549D4*xyzi(i+4,j+
     $        8,k)+1.318925622057893D4*xyzi(i+2,j+12,k)-1.28961616378994
     $        D4*xyzi(i+2,j+10,k)-2.029116341627528D2*xyzi(i,j+14,k)+1.9
     $        53963884530212D2*xyzi(i,j+12,k))+angi
            angi=yrk(225)*(1.04366482991984D1*xyzi(i+14,j,k)-9.497349952
     $        270546D2*xyzi(i+12,j+2,k)+1.04470849474976D4*xyzi(i+10,j+4
     $        ,k)-3.13412548424928D4*xyzi(i+8,j+6,k)+3.13412548424928D4*
     $        xyzi(i+6,j+8,k)-1.04470849474976D4*xyzi(i+4,j+10,k)+9.4973
     $        49952270546D2*xyzi(i+2,j+12,k)-1.04366482991984D1*xyzi(i,j
     $        +14,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(197)*(1.461130761887776D2*xyzi(i+13,j+1,k)-3.798939
     $        980908219D3*xyzi(i+11,j+3,k)+2.08941698949952D4*xyzi(i+9,j
     $        +5,k)-3.581857696284892D4*xyzi(i+7,j+7,k)+2.08941698949952
     $        D4*xyzi(i+5,j+9,k)-3.798939980908219D3*xyzi(i+3,j+11,k)+1.
     $        461130761887776D2*xyzi(i+1,j+13,k))+angi
            angi=yrk(199)*(-2.434939609953033D3*xyzi(i+13,j+1,k)+4.22056
     $        1990585258D4*xyzi(i+11,j+3,k)+2.344756661436254D3*xyzi(i+1
     $        1,j+1,k)-1.160654547410946D5*xyzi(i+9,j+5,k)-4.29872054596
     $        6466D4*xyzi(i+9,j+3,k)+1.547539396547928D5*xyzi(i+7,j+5,k)
     $        +1.160654547410946D5*xyzi(i+5,j+9,k)-1.547539396547928D5*x
     $        yzi(i+5,j+7,k)-4.220561990585258D4*xyzi(i+3,j+11,k)+4.2987
     $        20545966466D4*xyzi(i+3,j+9,k)+2.434939609953033D3*xyzi(i+1
     $        ,j+13,k)-2.344756661436254D3*xyzi(i+1,j+11,k))+angi
            angi=yrk(201)*(1.493389191601027D4*xyzi(i+13,j+1,k)-1.493389
     $        191601027D5*xyzi(i+11,j+3,k)-2.654914118401826D4*xyzi(i+11
     $        ,j+1,k)+3.28545622152226D4*xyzi(i+9,j+5,k)+2.9204055302420
     $        09D5*xyzi(i+9,j+3,k)+1.168162212096803D4*xyzi(i+9,j+1,k)+3
     $        .942547465826712D5*xyzi(i+7,j+7,k)-3.50448663629041D5*xyzi
     $        (i+7,j+5,k)-1.401794654516164D5*xyzi(i+7,j+3,k)+3.28545622
     $        152226D4*xyzi(i+5,j+9,k)-3.50448663629041D5*xyzi(i+5,j+7,k
     $        )+2.943768774483945D5*xyzi(i+5,j+5,k)-1.493389191601027D5*
     $        xyzi(i+3,j+11,k)+2.920405530242009D5*xyzi(i+3,j+9,k)-1.401
     $        794654516164D5*xyzi(i+3,j+7,k)+1.493389191601027D4*xyzi(i+
     $        1,j+13,k)-2.654914118401826D4*xyzi(i+1,j+11,k)+1.168162212
     $        096803D4*xyzi(i+1,j+9,k))+angi
            angi=yrk(203)*(-5.124740774410829D4*xyzi(i+13,j+1,k)+2.04989
     $        6309764331D5*xyzi(i+11,j+3,k)+1.252714411522647D5*xyzi(i+1
     $        1,j+1,k)+5.637214851851912D5*xyzi(i+9,j+5,k)-6.26357205761
     $        3235D5*xyzi(i+9,j+3,k)-1.002171529218118D5*xyzi(i+9,j+1,k)
     $        -7.516286469135882D5*xyzi(i+7,j+5,k)+6.013029175308706D5*x
     $        yzi(i+7,j+3,k)+2.614360511003785D4*xyzi(i+7,j+1,k)-5.63721
     $        4851851912D5*xyzi(i+5,j+9,k)+7.516286469135882D5*xyzi(i+5,
     $        j+7,k)-1.83005235770265D5*xyzi(i+5,j+3,k)-2.04989630976433
     $        1D5*xyzi(i+3,j+11,k)+6.263572057613235D5*xyzi(i+3,j+9,k)-6
     $        .013029175308706D5*xyzi(i+3,j+7,k)+1.83005235770265D5*xyzi
     $        (i+3,j+5,k)+5.124740774410829D4*xyzi(i+1,j+13,k)-1.2527144
     $        11522647D5*xyzi(i+1,j+11,k)+1.002171529218118D5*xyzi(i+1,j
     $        +9,k)-2.614360511003785D4*xyzi(i+1,j+7,k))+angi
            angi=yrk(205)*(1.103977290970156D5*xyzi(i+13,j+1,k)+7.359848
     $        606467704D4*xyzi(i+11,j+3,k)-3.271043825096757D5*xyzi(i+11
     $        ,j+1,k)-6.991856176144319D5*xyzi(i+9,j+5,k)+1.090347941698
     $        919D5*xyzi(i+9,j+3,k)+3.532727331104498D5*xyzi(i+9,j+1,k)-
     $        1.324772749164187D6*xyzi(i+7,j+7,k)+1.962626295058055D6*xy
     $        zi(i+7,j+5,k)-4.710303108139331D5*xyzi(i+7,j+3,k)-1.638366
     $        298483245D5*xyzi(i+7,j+1,k)-6.991856176144319D5*xyzi(i+5,j
     $        +9,k)+1.962626295058055D6*xyzi(i+5,j+7,k)-1.64860608784876
     $        6D6*xyzi(i+5,j+5,k)+3.822854696460906D5*xyzi(i+5,j+3,k)+2.
     $        730610497472076D4*xyzi(i+5,j+1,k)+7.359848606467704D4*xyzi
     $        (i+3,j+11,k)+1.090347941698919D5*xyzi(i+3,j+9,k)-4.7103031
     $        08139331D5*xyzi(i+3,j+7,k)+3.822854696460906D5*xyzi(i+3,j+
     $        5,k)-9.102034991573586D4*xyzi(i+3,j+3,k)+1.103977290970156
     $        D5*xyzi(i+1,j+13,k)-3.271043825096757D5*xyzi(i+1,j+11,k)+3
     $        .532727331104498D5*xyzi(i+1,j+9,k)-1.638366298483245D5*xyz
     $        i(i+1,j+7,k)+2.730610497472076D4*xyzi(i+1,j+5,k))+angi
            angi=yrk(207)*(-1.512305126981374D5*xyzi(i+13,j+1,k)-6.04922
     $        0507925497D5*xyzi(i+11,j+3,k)+5.041017089937914D5*xyzi(i+1
     $        1,j+1,k)-7.561525634906871D5*xyzi(i+9,j+5,k)+1.51230512698
     $        1374D6*xyzi(i+9,j+3,k)-6.45250187512053D5*xyzi(i+9,j+1,k)+
     $        1.008203417987583D6*xyzi(i+7,j+5,k)-1.290500375024106D6*xy
     $        zi(i+7,j+3,k)+3.927609837029888D5*xyzi(i+7,j+1,k)+7.561525
     $        634906871D5*xyzi(i+5,j+9,k)-1.008203417987583D6*xyzi(i+5,j
     $        +7,k)+3.927609837029888D5*xyzi(i+5,j+3,k)-1.12217423915139
     $        6D5*xyzi(i+5,j+1,k)+6.049220507925497D5*xyzi(i+3,j+11,k)-1
     $        .512305126981374D6*xyzi(i+3,j+9,k)+1.290500375024106D6*xyz
     $        i(i+3,j+7,k)-3.927609837029888D5*xyzi(i+3,j+5,k)+1.1812360
     $        41211996D4*xyzi(i+3,j+1,k)+1.512305126981374D5*xyzi(i+1,j+
     $        13,k)-5.041017089937914D5*xyzi(i+1,j+11,k)+6.4525018751205
     $        3D5*xyzi(i+1,j+9,k)-3.927609837029888D5*xyzi(i+1,j+7,k)+1.
     $        122174239151396D5*xyzi(i+1,j+5,k)-1.181236041211996D4*xyzi
     $        (i+1,j+3,k))+angi
            angi=yrk(209)*(1.151285875227237D5*xyzi(i+13,j+1,k)+6.907715
     $        251363424D5*xyzi(i+11,j+3,k)-4.093460889696844D5*xyzi(i+11
     $        ,j+1,k)+1.726928812840856D6*xyzi(i+9,j+5,k)-2.046730444848
     $        422D6*xyzi(i+9,j+3,k)+5.730845245575582D5*xyzi(i+9,j+1,k)+
     $        2.302571750454475D6*xyzi(i+7,j+7,k)-4.093460889696844D6*xy
     $        zi(i+7,j+5,k)+2.292338098230233D6*xyzi(i+7,j+3,k)-3.986674
     $        953443883D5*xyzi(i+7,j+1,k)+1.726928812840856D6*xyzi(i+5,j
     $        +9,k)-4.093460889696844D6*xyzi(i+5,j+7,k)+3.43850714734534
     $        9D6*xyzi(i+5,j+5,k)-1.196002486033165D6*xyzi(i+5,j+3,k)+1.
     $        423812483372815D5*xyzi(i+5,j+1,k)+6.907715251363424D5*xyzi
     $        (i+3,j+11,k)-2.046730444848422D6*xyzi(i+3,j+9,k)+2.2923380
     $        98230233D6*xyzi(i+3,j+7,k)-1.196002486033165D6*xyzi(i+3,j+
     $        5,k)+2.847624966745631D5*xyzi(i+3,j+3,k)-2.397999971996321
     $        D4*xyzi(i+3,j+1,k)+1.151285875227237D5*xyzi(i+1,j+13,k)-4.
     $        093460889696844D5*xyzi(i+1,j+11,k)+5.730845245575582D5*xyz
     $        i(i+1,j+9,k)-3.986674953443883D5*xyzi(i+1,j+7,k)+1.4238124
     $        83372815D5*xyzi(i+1,j+5,k)-2.397999971996321D4*xyzi(i+1,j+
     $        3,k)+1.410588218821365D3*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(212)*(-6.386185008998833D4*xyzi(i+13,j,k+1)-3.83171
     $        10053993D5*xyzi(i+11,j+2,k+1)+1.986813113910748D5*xyzi(i+1
     $        1,j,k+1)-9.579277513498249D5*xyzi(i+9,j+4,k+1)+9.934065569
     $        553739D5*xyzi(i+9,j+2,k+1)-2.384175736692897D5*xyzi(i+9,j,
     $        k+1)-1.277237001799767D6*xyzi(i+7,j+6,k+1)+1.9868131139107
     $        48D6*xyzi(i+7,j+4,k+1)-9.53670294677159D5*xyzi(i+7,j+2,k+1
     $        )+1.382130861850955D5*xyzi(i+7,j,k+1)-9.579277513498249D5*
     $        xyzi(i+5,j+8,k+1)+1.986813113910748D6*xyzi(i+5,j+6,k+1)-1.
     $        430505442015738D6*xyzi(i+5,j+4,k+1)+4.146392585552865D5*xy
     $        zi(i+5,j+2,k+1)-3.948945319574157D4*xyzi(i+5,j,k+1)-3.8317
     $        110053993D5*xyzi(i+3,j+10,k+1)+9.934065569553739D5*xyzi(i+
     $        3,j+8,k+1)-9.53670294677159D5*xyzi(i+3,j+6,k+1)+4.14639258
     $        5552865D5*xyzi(i+3,j+4,k+1)-7.897890639148315D4*xyzi(i+3,j
     $        +2,k+1)+4.988141456304199D3*xyzi(i+3,j,k+1)-6.386185008998
     $        833D4*xyzi(i+1,j+12,k+1)+1.986813113910748D5*xyzi(i+1,j+10
     $        ,k+1)-2.384175736692897D5*xyzi(i+1,j+8,k+1)+1.382130861850
     $        955D5*xyzi(i+1,j+6,k+1)-3.948945319574157D4*xyzi(i+1,j+4,k
     $        +1)+4.988141456304199D3*xyzi(i+1,j+2,k+1)-1.95613390443301
     $        9D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(214)*(4.83636804631037D4*xyzi(i+13,j,k+1)+9.6727360
     $        92620741D4*xyzi(i+11,j+2,k+1)-1.432997939647517D5*xyzi(i+1
     $        1,j,k+1)-2.418184023155185D5*xyzi(i+9,j+4,k+1)-1.432997939
     $        647517D5*xyzi(i+9,j+2,k+1)+1.604957692405219D5*xyzi(i+9,j,
     $        k+1)-9.672736092620741D5*xyzi(i+7,j+6,k+1)+8.5979876378851
     $        03D5*xyzi(i+7,j+4,k+1)-8.373692308201144D4*xyzi(i+7,j,k+1)
     $        -1.209092011577593D6*xyzi(i+5,j+8,k+1)+2.006197115506524D6
     $        *xyzi(i+5,j+6,k+1)-9.629746154431316D5*xyzi(i+5,j+4,k+1)+8
     $        .373692308201144D4*xyzi(i+5,j+2,k+1)+1.993736263857415D4*x
     $        yzi(i+5,j,k+1)-6.770915264834519D5*xyzi(i+3,j+10,k+1)+1.57
     $        6297733612269D6*xyzi(i+3,j+8,k+1)-1.283966153924175D6*xyzi
     $        (i+3,j+6,k+1)+4.186846154100572D5*xyzi(i+3,j+4,k+1)-3.9874
     $        7252771483D4*xyzi(i+3,j+2,k+1)-1.678935801143087D3*xyzi(i+
     $        3,j,k+1)-1.450910413893111D5*xyzi(i+1,j+12,k+1)+4.29899381
     $        8942552D5*xyzi(i+1,j+10,k+1)-4.814873077215658D5*xyzi(i+1,
     $        j+8,k+1)+2.512107692460343D5*xyzi(i+1,j+6,k+1)-5.981208791
     $        572246D4*xyzi(i+1,j+4,k+1)+5.036807403429259D3*xyzi(i+1,j+
     $        2,k+1))+angi
            angi=yrk(216)*(-2.742853631361481D4*xyzi(i+13,j,k+1)+1.64571
     $        2178816889D5*xyzi(i+11,j+2,k+1)+7.314276350297282D4*xyzi(i
     $        +11,j,k+1)+7.954275530948294D5*xyzi(i+9,j+4,k+1)-5.1199934
     $        45208098D5*xyzi(i+9,j+2,k+1)-7.021705296285391D4*xyzi(i+9,
     $        j,k+1)+9.874273072901331D5*xyzi(i+7,j+6,k+1)-1.60914079706
     $        5402D6*xyzi(i+7,j+4,k+1)+5.617364237028313D5*xyzi(i+7,j+2,
     $        k+1)+2.849387656463637D4*xyzi(i+7,j,k+1)+2.468568268225333
     $        D5*xyzi(i+5,j+8,k+1)-1.02399868904162D6*xyzi(i+5,j+6,k+1)+
     $        9.830387414799548D5*xyzi(i+5,j+4,k+1)-2.564448890817273D5*
     $        xyzi(i+5,j+2,k+1)-4.070553794948053D3*xyzi(i+5,j,k+1)-2.74
     $        2853631361481D5*xyzi(i+3,j+10,k+1)+3.657138175148641D5*xyz
     $        i(i+3,j+8,k+1)-1.424693828231818D5*xyzi(i+3,j+4,k+1)+4.070
     $        553794948053D4*xyzi(i+3,j+2,k+1)-1.37142681568074D5*xyzi(i
     $        +1,j+12,k+1)+3.657138175148641D5*xyzi(i+1,j+10,k+1)-3.5108
     $        52648142696D5*xyzi(i+1,j+8,k+1)+1.424693828231818D5*xyzi(i
     $        +1,j+6,k+1)-2.035276897474026D4*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(218)*(1.135649295191801D4*xyzi(i+13,j,k+1)-2.044168
     $        731345241D5*xyzi(i+11,j+2,k+1)-2.523665100426224D4*xyzi(i+
     $        11,j,k+1)-2.839123237979502D5*xyzi(i+9,j+4,k+1)+4.79496369
     $        0809825D5*xyzi(i+9,j+2,k+1)+1.817038872306881D4*xyzi(i+9,j
     $        ,k+1)+4.088337462690482D5*xyzi(i+7,j+6,k+1)+1.514199060255
     $        734D5*xyzi(i+7,j+4,k+1)-3.634077744613762D5*xyzi(i+7,j+2,k
     $        +1)-4.213423472015956D3*xyzi(i+7,j,k+1)+7.154590559708344D
     $        5*xyzi(i+5,j+8,k+1)-1.059939342179014D6*xyzi(i+5,j+6,k+1)+
     $        2.543854421229634D5*xyzi(i+5,j+4,k+1)+8.848189291233507D4*
     $        xyzi(i+5,j+2,k+1)+1.589909013268521D5*xyzi(i+3,j+10,k+1)-5
     $        .29969671089507D5*xyzi(i+3,j+8,k+1)+5.087708842459267D5*xy
     $        zi(i+3,j+6,k+1)-1.474698215205585D5*xyzi(i+3,j+4,k+1)-7.94
     $        9545066342605D4*xyzi(i+1,j+12,k+1)+1.766565570298357D5*xyz
     $        i(i+1,j+10,k+1)-1.271927210614817D5*xyzi(i+1,j+8,k+1)+2.94
     $        9396430411169D4*xyzi(i+1,j+6,k+1))+angi
            angi=yrk(220)*(-3.271851789497149D3*xyzi(i+13,j,k+1)+1.11242
     $        9608429031D5*xyzi(i+11,j+2,k+1)+5.331906619921279D3*xyzi(i
     $        +11,j,k+1)-1.799518484223432D5*xyzi(i+9,j+4,k+1)-1.8661673
     $        16972448D5*xyzi(i+9,j+2,k+1)-2.132762647968512D3*xyzi(i+9,
     $        j,k+1)-4.318844362136236D5*xyzi(i+7,j+6,k+1)+4.79871595792
     $        9151D5*xyzi(i+7,j+4,k+1)+7.677945532686642D4*xyzi(i+7,j+2,
     $        k+1)+1.079711090534059D5*xyzi(i+5,j+8,k+1)+2.2394007803669
     $        37D5*xyzi(i+5,j+6,k+1)-2.687280936440325D5*xyzi(i+5,j+4,k+
     $        1)+2.159422181068118D5*xyzi(i+3,j+10,k+1)-3.99892996494095
     $        9D5*xyzi(i+3,j+8,k+1)+1.79152062429355D5*xyzi(i+3,j+6,k+1)
     $        -2.944666610547434D4*xyzi(i+1,j+12,k+1)+4.798715957929151D
     $        4*xyzi(i+1,j+10,k+1)-1.919486383171661D4*xyzi(i+1,j+8,k+1)
     $        )+angi
            angi=yrk(222)*(5.973556766404109D2*xyzi(i+13,j,k+1)-3.225720
     $        653858219D4*xyzi(i+11,j+2,k+1)-5.309828236803652D2*xyzi(i+
     $        11,j,k+1)+1.64272811076113D5*xyzi(i+9,j+4,k+1)+2.920405530
     $        242009D4*xyzi(i+9,j+2,k+1)-7.885094931653424D4*xyzi(i+7,j+
     $        6,k+1)-1.752243318145205D5*xyzi(i+7,j+4,k+1)-1.77414635962
     $        202D5*xyzi(i+5,j+8,k+1)+2.453140645403287D5*xyzi(i+5,j+6,k
     $        +1)+9.199277420262327D4*xyzi(i+3,j+10,k+1)-8.7612165907260
     $        26D4*xyzi(i+3,j+8,k+1)-6.570912443044519D3*xyzi(i+1,j+12,k
     $        +1)+5.840811060484017D3*xyzi(i+1,j+10,k+1))+angi
            angi=yrk(224)*(-5.522555184144841D1*xyzi(i+13,j,k+1)+4.30759
     $        3043632976D3*xyzi(i+11,j+2,k+1)-3.948626956663561D4*xyzi(i
     $        +9,j+4,k+1)+9.476704695992546D4*xyzi(i+7,j+6,k+1)-7.107528
     $        52199441D4*xyzi(i+5,j+8,k+1)+1.579450782665424D4*xyzi(i+3,
     $        j+10,k+1)-7.179321739388293D2*xyzi(i+1,j+12,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(198)*(-7.179321739388293D2*xyzi(i+12,j+1,k+1)+1.579
     $        450782665424D4*xyzi(i+10,j+3,k+1)-7.10752852199441D4*xyzi(
     $        i+8,j+5,k+1)+9.476704695992546D4*xyzi(i+6,j+7,k+1)-3.94862
     $        6956663561D4*xyzi(i+4,j+9,k+1)+4.307593043632976D3*xyzi(i+
     $        2,j+11,k+1)-5.522555184144841D1*xyzi(i,j+13,k+1))+angi
            angi=yrk(200)*(6.570912443044519D3*xyzi(i+12,j+1,k+1)-9.1992
     $        77420262327D4*xyzi(i+10,j+3,k+1)-5.840811060484017D3*xyzi(
     $        i+10,j+1,k+1)+1.77414635962202D5*xyzi(i+8,j+5,k+1)+8.76121
     $        6590726026D4*xyzi(i+8,j+3,k+1)+7.885094931653424D4*xyzi(i+
     $        6,j+7,k+1)-2.453140645403287D5*xyzi(i+6,j+5,k+1)-1.6427281
     $        1076113D5*xyzi(i+4,j+9,k+1)+1.752243318145205D5*xyzi(i+4,j
     $        +7,k+1)+3.225720653858219D4*xyzi(i+2,j+11,k+1)-2.920405530
     $        242009D4*xyzi(i+2,j+9,k+1)-5.973556766404109D2*xyzi(i,j+13
     $        ,k+1)+5.309828236803652D2*xyzi(i,j+11,k+1))+angi
            angi=yrk(202)*(-2.944666610547434D4*xyzi(i+12,j+1,k+1)+2.159
     $        422181068118D5*xyzi(i+10,j+3,k+1)+4.798715957929151D4*xyzi
     $        (i+10,j+1,k+1)+1.079711090534059D5*xyzi(i+8,j+5,k+1)-3.998
     $        929964940959D5*xyzi(i+8,j+3,k+1)-1.919486383171661D4*xyzi(
     $        i+8,j+1,k+1)-4.318844362136236D5*xyzi(i+6,j+7,k+1)+2.23940
     $        0780366937D5*xyzi(i+6,j+5,k+1)+1.79152062429355D5*xyzi(i+6
     $        ,j+3,k+1)-1.799518484223432D5*xyzi(i+4,j+9,k+1)+4.79871595
     $        7929151D5*xyzi(i+4,j+7,k+1)-2.687280936440325D5*xyzi(i+4,j
     $        +5,k+1)+1.112429608429031D5*xyzi(i+2,j+11,k+1)-1.866167316
     $        972448D5*xyzi(i+2,j+9,k+1)+7.677945532686642D4*xyzi(i+2,j+
     $        7,k+1)-3.271851789497149D3*xyzi(i,j+13,k+1)+5.331906619921
     $        279D3*xyzi(i,j+11,k+1)-2.132762647968512D3*xyzi(i,j+9,k+1)
     $        )+angi
            angi=yrk(204)*(7.949545066342605D4*xyzi(i+12,j+1,k+1)-1.5899
     $        09013268521D5*xyzi(i+10,j+3,k+1)-1.766565570298357D5*xyzi(
     $        i+10,j+1,k+1)-7.154590559708344D5*xyzi(i+8,j+5,k+1)+5.2996
     $        9671089507D5*xyzi(i+8,j+3,k+1)+1.271927210614817D5*xyzi(i+
     $        8,j+1,k+1)-4.088337462690482D5*xyzi(i+6,j+7,k+1)+1.0599393
     $        42179014D6*xyzi(i+6,j+5,k+1)-5.087708842459267D5*xyzi(i+6,
     $        j+3,k+1)-2.949396430411169D4*xyzi(i+6,j+1,k+1)+2.839123237
     $        979502D5*xyzi(i+4,j+9,k+1)-1.514199060255734D5*xyzi(i+4,j+
     $        7,k+1)-2.543854421229634D5*xyzi(i+4,j+5,k+1)+1.47469821520
     $        5585D5*xyzi(i+4,j+3,k+1)+2.044168731345241D5*xyzi(i+2,j+11
     $        ,k+1)-4.794963690809825D5*xyzi(i+2,j+9,k+1)+3.634077744613
     $        762D5*xyzi(i+2,j+7,k+1)-8.848189291233507D4*xyzi(i+2,j+5,k
     $        +1)-1.135649295191801D4*xyzi(i,j+13,k+1)+2.523665100426224
     $        D4*xyzi(i,j+11,k+1)-1.817038872306881D4*xyzi(i,j+9,k+1)+4.
     $        213423472015956D3*xyzi(i,j+7,k+1))+angi
            angi=yrk(206)*(-1.37142681568074D5*xyzi(i+12,j+1,k+1)-2.7428
     $        53631361481D5*xyzi(i+10,j+3,k+1)+3.657138175148641D5*xyzi(
     $        i+10,j+1,k+1)+2.468568268225333D5*xyzi(i+8,j+5,k+1)+3.6571
     $        38175148641D5*xyzi(i+8,j+3,k+1)-3.510852648142696D5*xyzi(i
     $        +8,j+1,k+1)+9.874273072901331D5*xyzi(i+6,j+7,k+1)-1.023998
     $        68904162D6*xyzi(i+6,j+5,k+1)+1.424693828231818D5*xyzi(i+6,
     $        j+1,k+1)+7.954275530948294D5*xyzi(i+4,j+9,k+1)-1.609140797
     $        065402D6*xyzi(i+4,j+7,k+1)+9.830387414799548D5*xyzi(i+4,j+
     $        5,k+1)-1.424693828231818D5*xyzi(i+4,j+3,k+1)-2.03527689747
     $        4026D4*xyzi(i+4,j+1,k+1)+1.645712178816889D5*xyzi(i+2,j+11
     $        ,k+1)-5.119993445208098D5*xyzi(i+2,j+9,k+1)+5.617364237028
     $        313D5*xyzi(i+2,j+7,k+1)-2.564448890817273D5*xyzi(i+2,j+5,k
     $        +1)+4.070553794948053D4*xyzi(i+2,j+3,k+1)-2.74285363136148
     $        1D4*xyzi(i,j+13,k+1)+7.314276350297282D4*xyzi(i,j+11,k+1)-
     $        7.021705296285391D4*xyzi(i,j+9,k+1)+2.849387656463637D4*xy
     $        zi(i,j+7,k+1)-4.070553794948053D3*xyzi(i,j+5,k+1))+angi
            angi=yrk(208)*(1.450910413893111D5*xyzi(i+12,j+1,k+1)+6.7709
     $        15264834519D5*xyzi(i+10,j+3,k+1)-4.298993818942552D5*xyzi(
     $        i+10,j+1,k+1)+1.209092011577593D6*xyzi(i+8,j+5,k+1)-1.5762
     $        97733612269D6*xyzi(i+8,j+3,k+1)+4.814873077215658D5*xyzi(i
     $        +8,j+1,k+1)+9.672736092620741D5*xyzi(i+6,j+7,k+1)-2.006197
     $        115506524D6*xyzi(i+6,j+5,k+1)+1.283966153924175D6*xyzi(i+6
     $        ,j+3,k+1)-2.512107692460343D5*xyzi(i+6,j+1,k+1)+2.41818402
     $        3155185D5*xyzi(i+4,j+9,k+1)-8.597987637885103D5*xyzi(i+4,j
     $        +7,k+1)+9.629746154431316D5*xyzi(i+4,j+5,k+1)-4.1868461541
     $        00572D5*xyzi(i+4,j+3,k+1)+5.981208791572246D4*xyzi(i+4,j+1
     $        ,k+1)-9.672736092620741D4*xyzi(i+2,j+11,k+1)+1.43299793964
     $        7517D5*xyzi(i+2,j+9,k+1)-8.373692308201144D4*xyzi(i+2,j+5,
     $        k+1)+3.98747252771483D4*xyzi(i+2,j+3,k+1)-5.03680740342925
     $        9D3*xyzi(i+2,j+1,k+1)-4.83636804631037D4*xyzi(i,j+13,k+1)+
     $        1.432997939647517D5*xyzi(i,j+11,k+1)-1.604957692405219D5*x
     $        yzi(i,j+9,k+1)+8.373692308201144D4*xyzi(i,j+7,k+1)-1.99373
     $        6263857415D4*xyzi(i,j+5,k+1)+1.678935801143087D3*xyzi(i,j+
     $        3,k+1))+angi
            angi=yrk(210)*(-6.386185008998833D4*xyzi(i+12,j+1,k+1)-3.831
     $        7110053993D5*xyzi(i+10,j+3,k+1)+1.986813113910748D5*xyzi(i
     $        +10,j+1,k+1)-9.579277513498249D5*xyzi(i+8,j+5,k+1)+9.93406
     $        5569553739D5*xyzi(i+8,j+3,k+1)-2.384175736692897D5*xyzi(i+
     $        8,j+1,k+1)-1.277237001799767D6*xyzi(i+6,j+7,k+1)+1.9868131
     $        13910748D6*xyzi(i+6,j+5,k+1)-9.53670294677159D5*xyzi(i+6,j
     $        +3,k+1)+1.382130861850955D5*xyzi(i+6,j+1,k+1)-9.5792775134
     $        98249D5*xyzi(i+4,j+9,k+1)+1.986813113910748D6*xyzi(i+4,j+7
     $        ,k+1)-1.430505442015738D6*xyzi(i+4,j+5,k+1)+4.146392585552
     $        865D5*xyzi(i+4,j+3,k+1)-3.948945319574157D4*xyzi(i+4,j+1,k
     $        +1)-3.8317110053993D5*xyzi(i+2,j+11,k+1)+9.934065569553739
     $        D5*xyzi(i+2,j+9,k+1)-9.53670294677159D5*xyzi(i+2,j+7,k+1)+
     $        4.146392585552865D5*xyzi(i+2,j+5,k+1)-7.897890639148315D4*
     $        xyzi(i+2,j+3,k+1)+4.988141456304199D3*xyzi(i+2,j+1,k+1)-6.
     $        386185008998833D4*xyzi(i,j+13,k+1)+1.986813113910748D5*xyz
     $        i(i,j+11,k+1)-2.384175736692897D5*xyzi(i,j+9,k+1)+1.382130
     $        861850955D5*xyzi(i,j+7,k+1)-3.948945319574157D4*xyzi(i,j+5
     $        ,k+1)+4.988141456304199D3*xyzi(i,j+3,k+1)-1.95613390443301
     $        9D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA(14)=angi
        endif
c...
c...  l=16
c...
        if(lmax.ge.16)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(273)*(1.867731876130107D5*xyzi(i,j,k+16)-7.22992984
     $        3084283D5*xyzi(i,j,k+14)+1.134351061587362D6*xyzi(i,j,k+12
     $        )-9.242860501822947D5*xyzi(i,j,k+10)+4.159287225820326D5*x
     $        yzi(i,j,k+8)-1.012696020199732D5*xyzi(i,j,k+6)+1.205590500
     $        237776D4*xyzi(i,j,k+4)-5.438754136411018D2*xyzi(i,j,k+2)+3
     $        .999083923831631D0*xyzi(i,j,k))+angi
            angi=yrk(275)*(-2.339238140133424D5*xyzi(i+16,j,k)-1.4035428
     $        84080055D6*xyzi(i+14,j+2,k)+9.507871150219725D5*xyzi(i+14,
     $        j,k)-3.274933396186794D6*xyzi(i+12,j+4,k)+4.75393557510986
     $        2D6*xyzi(i+12,j+2,k)-1.573716604174299D6*xyzi(i+12,j,k)-3.
     $        274933396186794D6*xyzi(i+10,j+6,k)+8.557084035197753D6*xyz
     $        i(i+10,j+4,k)-6.294866416697197D6*xyzi(i+10,j+2,k)+1.36000
     $        2003607419D6*xyzi(i+10,j,k)+4.753935575109862D6*xyzi(i+8,j
     $        +6,k)-7.868583020871497D6*xyzi(i+8,j+4,k)+4.08000601082225
     $        8D6*xyzi(i+8,j+2,k)-6.528009617315612D5*xyzi(i+8,j,k)+3.27
     $        4933396186794D6*xyzi(i+6,j+10,k)-4.753935575109862D6*xyzi(
     $        i+6,j+8,k)+2.720004007214838D6*xyzi(i+6,j+4,k)-1.305601923
     $        463122D6*xyzi(i+6,j+2,k)+1.702959030604073D5*xyzi(i+6,j,k)
     $        +3.274933396186794D6*xyzi(i+4,j+12,k)-8.557084035197753D6*
     $        xyzi(i+4,j+10,k)+7.868583020871497D6*xyzi(i+4,j+8,k)-2.720
     $        004007214838D6*xyzi(i+4,j+6,k)+1.702959030604073D5*xyzi(i+
     $        4,j+2,k)-2.162487657909934D4*xyzi(i+4,j,k)+1.4035428840800
     $        55D6*xyzi(i+2,j+14,k)-4.753935575109862D6*xyzi(i+2,j+12,k)
     $        +6.294866416697197D6*xyzi(i+2,j+10,k)-4.080006010822258D6*
     $        xyzi(i+2,j+8,k)+1.305601923463122D6*xyzi(i+2,j+6,k)-1.7029
     $        59030604073D5*xyzi(i+2,j+4,k)+9.755583419142558D2*xyzi(i+2
     $        ,j,k)+2.339238140133424D5*xyzi(i,j+16,k)-9.507871150219725
     $        D5*xyzi(i,j+14,k)+1.573716604174299D6*xyzi(i,j+12,k)-1.360
     $        002003607419D6*xyzi(i,j+10,k)+6.528009617315612D5*xyzi(i,j
     $        +8,k)-1.702959030604073D5*xyzi(i,j+6,k)+2.162487657909934D
     $        4*xyzi(i,j+4,k)-9.755583419142558D2*xyzi(i,j+2,k))+angi
            angi=yrk(277)*(1.618893696225904D5*xyzi(i+16,j,k)-6.26668527
     $        5713177D5*xyzi(i+14,j,k)-3.237787392451808D6*xyzi(i+12,j+4
     $        ,k)+6.266685275713177D5*xyzi(i+12,j+2,k)+9.724166807141138
     $        D5*xyzi(i+12,j,k)-1.036091965584579D7*xyzi(i+10,j+6,k)+1.1
     $        90670202385504D7*xyzi(i+10,j+4,k)-1.944833361428228D6*xyzi
     $        (i+10,j+2,k)-7.683292292062133D5*xyzi(i+10,j,k)-1.45700432
     $        6603314D7*xyzi(i+8,j+8,k)+2.82000837407093D7*xyzi(i+8,j+6,
     $        k)-1.653108357213993D7*xyzi(i+8,j+4,k)+2.30498768761864D6*
     $        xyzi(i+8,j+2,k)+3.226982762666096D5*xyzi(i+8,j,k)-1.036091
     $        965584579D7*xyzi(i+6,j+10,k)+2.82000837407093D7*xyzi(i+6,j
     $        +8,k)-2.722766705999518D7*xyzi(i+6,j+6,k)+1.07566092088869
     $        9D7*xyzi(i+6,j+4,k)-1.290793105066438D6*xyzi(i+6,j+2,k)-6.
     $        734572722085765D4*xyzi(i+6,j,k)-3.237787392451808D6*xyzi(i
     $        +4,j+12,k)+1.190670202385504D7*xyzi(i+4,j+10,k)-1.65310835
     $        7213993D7*xyzi(i+4,j+8,k)+1.075660920888699D7*xyzi(i+4,j+6
     $        ,k)-3.226982762666096D6*xyzi(i+4,j+4,k)+3.367286361042883D
     $        5*xyzi(i+4,j+2,k)+5.344898985782354D3*xyzi(i+4,j,k)+6.2666
     $        85275713177D5*xyzi(i+2,j+12,k)-1.944833361428228D6*xyzi(i+
     $        2,j+10,k)+2.30498768761864D6*xyzi(i+2,j+8,k)-1.29079310506
     $        6438D6*xyzi(i+2,j+6,k)+3.367286361042883D5*xyzi(i+2,j+4,k)
     $        -3.206939391469412D4*xyzi(i+2,j+2,k)+1.618893696225904D5*x
     $        yzi(i,j+16,k)-6.266685275713177D5*xyzi(i,j+14,k)+9.7241668
     $        07141138D5*xyzi(i,j+12,k)-7.683292292062133D5*xyzi(i,j+10,
     $        k)+3.226982762666096D5*xyzi(i,j+8,k)-6.734572722085765D4*x
     $        yzi(i,j+6,k)+5.344898985782354D3*xyzi(i,j+4,k))+angi
            angi=yrk(279)*(-8.653350795550605D4*xyzi(i+16,j,k)+8.6533507
     $        95550605D5*xyzi(i+14,j+2,k)+3.070543830679247D5*xyzi(i+14,
     $        j,k)+4.326675397775303D6*xyzi(i+12,j+4,k)-3.37759821374717
     $        2D6*xyzi(i+12,j+2,k)-4.23523286990241D5*xyzi(i+12,j,k)+5.7
     $        11211525063399D6*xyzi(i+10,j+6,k)-1.197512093964906D7*xyzi
     $        (i+10,j+4,k)+5.082279443882891D6*xyzi(i+10,j+2,k)+2.823488
     $        57993494D5*xyzi(i+10,j,k)-8.290468342833967D6*xyzi(i+8,j+6
     $        ,k)+1.143512874873651D7*xyzi(i+8,j+4,k)-3.670535153915422D
     $        6*xyzi(i+8,j+2,k)-9.035163455791807D4*xyzi(i+8,j,k)-5.7112
     $        11525063399D6*xyzi(i+6,j+10,k)+8.290468342833967D6*xyzi(i+
     $        6,j+8,k)-3.952884011908916D6*xyzi(i+6,j+4,k)+1.26492288381
     $        0853D6*xyzi(i+6,j+2,k)+1.09993294244422D4*xyzi(i+6,j,k)-4.
     $        326675397775303D6*xyzi(i+4,j+12,k)+1.197512093964906D7*xyz
     $        i(i+4,j+10,k)-1.143512874873651D7*xyzi(i+4,j+8,k)+3.952884
     $        011908916D6*xyzi(i+4,j+6,k)-1.64989941366633D5*xyzi(i+4,j+
     $        2,k)-8.653350795550605D5*xyzi(i+2,j+14,k)+3.37759821374717
     $        2D6*xyzi(i+2,j+12,k)-5.082279443882891D6*xyzi(i+2,j+10,k)+
     $        3.670535153915422D6*xyzi(i+2,j+8,k)-1.264922883810853D6*xy
     $        zi(i+2,j+6,k)+1.64989941366633D5*xyzi(i+2,j+4,k)+8.6533507
     $        95550605D4*xyzi(i,j+16,k)-3.070543830679247D5*xyzi(i,j+14,
     $        k)+4.23523286990241D5*xyzi(i,j+12,k)-2.82348857993494D5*xy
     $        zi(i,j+10,k)+9.035163455791807D4*xyzi(i,j+8,k)-1.099932942
     $        44422D4*xyzi(i,j+6,k))+angi
            angi=yrk(281)*(3.494105595363806D4*xyzi(i+16,j,k)-8.38585342
     $        8873134D5*xyzi(i+14,j+2,k)-1.082045603725566D5*xyzi(i+14,j
     $        ,k)-1.25787801433097D6*xyzi(i+12,j+4,k)+2.705114009313914D
     $        6*xyzi(i+12,j+2,k)+1.231293273204954D5*xyzi(i+12,j,k)+3.07
     $        4812923920149D6*xyzi(i+10,j+6,k)+1.190250164098122D6*xyzi(
     $        i+10,j+4,k)-3.20136251033288D6*xyzi(i+10,j+2,k)-6.08046060
     $        8419526D4*xyzi(i+10,j,k)+6.918329078820336D6*xyzi(i+8,j+8,
     $        k)-1.07122514768831D7*xyzi(i+8,j+6,k)+1.846939909807431D6*
     $        xyzi(i+8,j+4,k)+1.641724364273272D6*xyzi(i+8,j+2,k)+1.0944
     $        82909515515D4*xyzi(i+8,j,k)+3.074812923920149D6*xyzi(i+6,j
     $        +10,k)-1.07122514768831D7*xyzi(i+6,j+8,k)+1.03428634949216
     $        1D7*xyzi(i+6,j+6,k)-2.553793455536201D6*xyzi(i+6,j+4,k)-3.
     $        064552146643441D5*xyzi(i+6,j+2,k)-1.25787801433097D6*xyzi(
     $        i+4,j+12,k)+1.190250164098122D6*xyzi(i+4,j+10,k)+1.8469399
     $        09807431D6*xyzi(i+4,j+8,k)-2.553793455536201D6*xyzi(i+4,j+
     $        6,k)+7.661380366608603D5*xyzi(i+4,j+4,k)-8.385853428873134
     $        D5*xyzi(i+2,j+14,k)+2.705114009313914D6*xyzi(i+2,j+12,k)-3
     $        .20136251033288D6*xyzi(i+2,j+10,k)+1.641724364273272D6*xyz
     $        i(i+2,j+8,k)-3.064552146643441D5*xyzi(i+2,j+6,k)+3.4941055
     $        95363806D4*xyzi(i,j+16,k)-1.082045603725566D5*xyzi(i,j+14,
     $        k)+1.231293273204954D5*xyzi(i,j+12,k)-6.080460608419526D4*
     $        xyzi(i,j+10,k)+1.094482909515515D4*xyzi(i,j+8,k))+angi
            angi=yrk(283)*(-1.025589015787025D4*xyzi(i+16,j,k)+4.3074738
     $        66305504D5*xyzi(i+14,j+2,k)+2.580514297786707D4*xyzi(i+14,
     $        j,k)-7.999594323138793D5*xyzi(i+12,j+4,k)-1.10962114804828
     $        4D6*xyzi(i+12,j+2,k)-2.13559803954762D4*xyzi(i+12,j,k)-2.9
     $        33184585150891D6*xyzi(i+10,j+6,k)+3.122422300321916D6*xyzi
     $        (i+10,j+4,k)+9.396631374009528D5*xyzi(i+10,j+2,k)+5.800389
     $        737042918D3*xyzi(i+10,j,k)+4.257848591348067D6*xyzi(i+8,j+
     $        6,k)-3.523736765253573D6*xyzi(i+8,j+4,k)-2.610175381669313
     $        D5*xyzi(i+8,j+2,k)+2.933184585150891D6*xyzi(i+6,j+10,k)-4.
     $        257848591348067D6*xyzi(i+6,j+8,k)+1.218081844779013D6*xyzi
     $        (i+6,j+4,k)+7.999594323138793D5*xyzi(i+4,j+12,k)-3.1224223
     $        00321916D6*xyzi(i+4,j+10,k)+3.523736765253573D6*xyzi(i+4,j
     $        +8,k)-1.218081844779013D6*xyzi(i+4,j+6,k)-4.30747386630550
     $        4D5*xyzi(i+2,j+14,k)+1.109621148048284D6*xyzi(i+2,j+12,k)-
     $        9.396631374009528D5*xyzi(i+2,j+10,k)+2.610175381669313D5*x
     $        yzi(i+2,j+8,k)+1.025589015787025D4*xyzi(i,j+16,k)-2.580514
     $        297786707D4*xyzi(i,j+14,k)+2.13559803954762D4*xyzi(i,j+12,
     $        k)-5.800389737042918D3*xyzi(i,j+10,k))+angi
            angi=yrk(285)*(2.043022221812925D3*xyzi(i+16,j,k)-1.30753422
     $        1960272D5*xyzi(i+14,j+2,k)-3.69062078779109D3*xyzi(i+14,j,
     $        k)+7.436600887399047D5*xyzi(i+12,j+4,k)+2.398903512064209D
     $        5*xyzi(i+12,j+2,k)+1.654416215216696D3*xyzi(i+12,j,k)-1.58
     $        3276317962378D6*xyzi(i+10,j+4,k)-1.091914702043019D5*xyzi(
     $        i+10,j+2,k)-1.75291306631549D6*xyzi(i+8,j+8,k)+1.583276317
     $        962378D6*xyzi(i+8,j+6,k)+8.189360265322644D5*xyzi(i+8,j+4,
     $        k)+1.583276317962378D6*xyzi(i+6,j+8,k)-1.528680582860227D6
     $        *xyzi(i+6,j+6,k)+7.436600887399047D5*xyzi(i+4,j+12,k)-1.58
     $        3276317962378D6*xyzi(i+4,j+10,k)+8.189360265322644D5*xyzi(
     $        i+4,j+8,k)-1.307534221960272D5*xyzi(i+2,j+14,k)+2.39890351
     $        2064209D5*xyzi(i+2,j+12,k)-1.091914702043019D5*xyzi(i+2,j+
     $        10,k)+2.043022221812925D3*xyzi(i,j+16,k)-3.69062078779109D
     $        3*xyzi(i,j+14,k)+1.654416215216696D3*xyzi(i,j+12,k))+angi
            angi=yrk(287)*(-2.399407915132806D2*xyzi(i+16,j,k)+2.1594671
     $        23619526D4*xyzi(i+14,j+2,k)+2.322007659805941D2*xyzi(i+14,
     $        j,k)-2.183461202770854D5*xyzi(i+12,j+4,k)-2.11302697042340
     $        7D4*xyzi(i+12,j+2,k)+4.803614646095878D5*xyzi(i+10,j+6,k)+
     $        2.324329667465747D5*xyzi(i+10,j+4,k)-6.972989002397242D5*x
     $        yzi(i+8,j+6,k)-4.803614646095878D5*xyzi(i+6,j+10,k)+6.9729
     $        89002397242D5*xyzi(i+6,j+8,k)+2.183461202770854D5*xyzi(i+4
     $        ,j+12,k)-2.324329667465747D5*xyzi(i+4,j+10,k)-2.1594671236
     $        19526D4*xyzi(i+2,j+14,k)+2.113026970423407D4*xyzi(i+2,j+12
     $        ,k)+2.399407915132806D2*xyzi(i,j+16,k)-2.322007659805941D2
     $        *xyzi(i,j+14,k))+angi
            angi=yrk(289)*(1.077365958207155D1*xyzi(i+16,j,k)-1.29283914
     $        9848586D3*xyzi(i+14,j+2,k)+1.960806043937022D4*xyzi(i+12,j
     $        +4,k)-8.627546593322896D4*xyzi(i+10,j+6,k)+1.3865699882126
     $        08D5*xyzi(i+8,j+8,k)-8.627546593322896D4*xyzi(i+6,j+10,k)+
     $        1.960806043937022D4*xyzi(i+4,j+12,k)-1.292839149848586D3*x
     $        yzi(i+2,j+14,k)+1.077365958207155D1*xyzi(i,j+16,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(257)*(1.723785533131448D2*xyzi(i+15,j+1,k)-6.033249
     $        365960067D3*xyzi(i+13,j+3,k)+4.705934505448852D4*xyzi(i+11
     $        ,j+5,k)-1.232506656188985D5*xyzi(i+9,j+7,k)+1.232506656188
     $        985D5*xyzi(i+7,j+9,k)-4.705934505448852D4*xyzi(i+5,j+11,k)
     $        +6.033249365960067D3*xyzi(i+3,j+13,k)-1.723785533131448D2*
     $        xyzi(i+1,j+15,k))+angi
            angi=yrk(259)*(-3.359171081185929D3*xyzi(i+15,j+1,k)+8.39792
     $        7702964821D4*xyzi(i+13,j+3,k)+3.250810723728318D3*xyzi(i+1
     $        3,j+1,k)-3.930230164987537D5*xyzi(i+11,j+5,k)-8.4521078816
     $        93627D4*xyzi(i+11,j+3,k)+3.431153318639913D5*xyzi(i+9,j+7,
     $        k)+4.648659334931495D5*xyzi(i+9,j+5,k)+3.431153318639913D5
     $        *xyzi(i+7,j+9,k)-7.969130288453991D5*xyzi(i+7,j+7,k)-3.930
     $        230164987537D5*xyzi(i+5,j+11,k)+4.648659334931495D5*xyzi(i
     $        +5,j+9,k)+8.397927702964821D4*xyzi(i+3,j+13,k)-8.452107881
     $        693627D4*xyzi(i+3,j+11,k)-3.359171081185929D3*xyzi(i+1,j+1
     $        5,k)+3.250810723728318D3*xyzi(i+1,j+13,k))+angi
            angi=yrk(261)*(2.45162666617551D4*xyzi(i+15,j+1,k)-4.0043235
     $        54753333D5*xyzi(i+13,j+3,k)-4.428744945349308D4*xyzi(i+13,
     $        j+1,k)+7.436600887399047D5*xyzi(i+11,j+5,k)+7.676491238605
     $        468D5*xyzi(i+11,j+3,k)+1.985299458260035D4*xyzi(i+11,j+1,k
     $        )+1.168608710876993D6*xyzi(i+9,j+7,k)-2.111035090616504D6*
     $        xyzi(i+9,j+5,k)-3.639715673476731D5*xyzi(i+9,j+3,k)-1.1686
     $        08710876993D6*xyzi(i+7,j+9,k)+1.310297642451623D6*xyzi(i+7
     $        ,j+5,k)-7.436600887399047D5*xyzi(i+5,j+11,k)+2.11103509061
     $        6504D6*xyzi(i+5,j+9,k)-1.310297642451623D6*xyzi(i+5,j+7,k)
     $        +4.004323554753333D5*xyzi(i+3,j+13,k)-7.676491238605468D5*
     $        xyzi(i+3,j+11,k)+3.639715673476731D5*xyzi(i+3,j+9,k)-2.451
     $        62666617551D4*xyzi(i+1,j+15,k)+4.428744945349308D4*xyzi(i+
     $        1,j+13,k)-1.985299458260035D4*xyzi(i+1,j+11,k))+angi
            angi=yrk(263)*(-1.025589015787025D5*xyzi(i+15,j+1,k)+9.23030
     $        1142083223D5*xyzi(i+13,j+3,k)+2.580514297786707D5*xyzi(i+1
     $        3,j+1,k)+7.999594323138793D5*xyzi(i+11,j+5,k)-2.5805142977
     $        86707D6*xyzi(i+11,j+3,k)-2.13559803954762D5*xyzi(i+11,j+1,
     $        k)-2.933184585150891D6*xyzi(i+9,j+7,k)+5.677131455130756D5
     $        *xyzi(i+9,j+5,k)+2.349157843502382D6*xyzi(i+9,j+3,k)+5.800
     $        389737042918D4*xyzi(i+9,j+1,k)-2.933184585150891D6*xyzi(i+
     $        7,j+9,k)+6.812557746156907D6*xyzi(i+7,j+7,k)-2.81898941220
     $        2858D6*xyzi(i+7,j+5,k)-6.960467684451502D5*xyzi(i+7,j+3,k)
     $        +7.999594323138793D5*xyzi(i+5,j+11,k)+5.677131455130756D5*
     $        xyzi(i+5,j+9,k)-2.818989412202858D6*xyzi(i+5,j+7,k)+1.4616
     $        98213734815D6*xyzi(i+5,j+5,k)+9.230301142083223D5*xyzi(i+3
     $        ,j+13,k)-2.580514297786707D6*xyzi(i+3,j+11,k)+2.3491578435
     $        02382D6*xyzi(i+3,j+9,k)-6.960467684451502D5*xyzi(i+3,j+7,k
     $        )-1.025589015787025D5*xyzi(i+1,j+15,k)+2.580514297786707D5
     $        *xyzi(i+1,j+13,k)-2.13559803954762D5*xyzi(i+1,j+11,k)+5.80
     $        0389737042918D4*xyzi(i+1,j+9,k))+angi
            angi=yrk(265)*(2.795284476291045D5*xyzi(i+15,j+1,k)-8.385853
     $        428873134D5*xyzi(i+13,j+3,k)-8.656364829804525D5*xyzi(i+13
     $        ,j+1,k)-4.192926714436567D6*xyzi(i+11,j+5,k)+3.46254593192
     $        181D6*xyzi(i+11,j+3,k)+9.850346185639632D5*xyzi(i+11,j+1,k
     $        )-3.074812923920149D6*xyzi(i+9,j+7,k)+9.522001312784978D6*
     $        xyzi(i+9,j+5,k)-4.925173092819816D6*xyzi(i+9,j+3,k)-4.8643
     $        68486735621D5*xyzi(i+9,j+1,k)+3.074812923920149D6*xyzi(i+7
     $        ,j+9,k)-5.910207711383779D6*xyzi(i+7,j+5,k)+2.918621092041
     $        373D6*xyzi(i+7,j+3,k)+8.755863276124118D4*xyzi(i+7,j+1,k)+
     $        4.192926714436567D6*xyzi(i+5,j+11,k)-9.522001312784978D6*x
     $        yzi(i+5,j+9,k)+5.910207711383779D6*xyzi(i+5,j+7,k)-6.12910
     $        4293286882D5*xyzi(i+5,j+3,k)+8.385853428873134D5*xyzi(i+3,
     $        j+13,k)-3.46254593192181D6*xyzi(i+3,j+11,k)+4.925173092819
     $        816D6*xyzi(i+3,j+9,k)-2.918621092041373D6*xyzi(i+3,j+7,k)+
     $        6.129104293286882D5*xyzi(i+3,j+5,k)-2.795284476291045D5*xy
     $        zi(i+1,j+15,k)+8.656364829804525D5*xyzi(i+1,j+13,k)-9.8503
     $        46185639632D5*xyzi(i+1,j+11,k)+4.864368486735621D5*xyzi(i+
     $        1,j+9,k)-8.755863276124118D4*xyzi(i+1,j+7,k))+angi
            angi=yrk(267)*(-5.192010477330363D5*xyzi(i+15,j+1,k)-8.65335
     $        0795550605D5*xyzi(i+13,j+3,k)+1.842326298407548D6*xyzi(i+1
     $        3,j+1,k)+2.942139270487206D6*xyzi(i+11,j+5,k)+1.2282175322
     $        71699D6*xyzi(i+11,j+3,k)-2.541139721941446D6*xyzi(i+11,j+1
     $        ,k)+9.518685875105666D6*xyzi(i+9,j+7,k)-1.166806655658114D
     $        7*xyzi(i+9,j+5,k)+8.470465739804819D5*xyzi(i+9,j+3,k)+1.69
     $        4093147960964D6*xyzi(i+9,j+1,k)+9.518685875105666D6*xyzi(i
     $        +7,j+9,k)-2.210791558089058D7*xyzi(i+7,j+7,k)+1.5246838331
     $        64867D7*xyzi(i+7,j+5,k)-2.258790863947952D6*xyzi(i+7,j+3,k
     $        )-5.421098073475084D5*xyzi(i+7,j+1,k)+2.942139270487206D6*
     $        xyzi(i+5,j+11,k)-1.166806655658114D7*xyzi(i+5,j+9,k)+1.524
     $        683833164867D7*xyzi(i+5,j+7,k)-7.905768023817831D6*xyzi(i+
     $        5,j+5,k)+1.264922883810853D6*xyzi(i+5,j+3,k)+6.59959765466
     $        532D4*xyzi(i+5,j+1,k)-8.653350795550605D5*xyzi(i+3,j+13,k)
     $        +1.228217532271699D6*xyzi(i+3,j+11,k)+8.470465739804819D5*
     $        xyzi(i+3,j+9,k)-2.258790863947952D6*xyzi(i+3,j+7,k)+1.2649
     $        22883810853D6*xyzi(i+3,j+5,k)-2.19986588488844D5*xyzi(i+3,
     $        j+3,k)-5.192010477330363D5*xyzi(i+1,j+15,k)+1.842326298407
     $        548D6*xyzi(i+1,j+13,k)-2.541139721941446D6*xyzi(i+1,j+11,k
     $        )+1.694093147960964D6*xyzi(i+1,j+9,k)-5.421098073475084D5*
     $        xyzi(i+1,j+7,k)+6.59959765466532D4*xyzi(i+1,j+5,k))+angi
            angi=yrk(269)*(6.475574784903617D5*xyzi(i+15,j+1,k)+3.237787
     $        392451808D6*xyzi(i+13,j+3,k)-2.506674110285271D6*xyzi(i+13
     $        ,j+1,k)+5.828017306413255D6*xyzi(i+11,j+5,k)-1.00266964411
     $        4108D7*xyzi(i+11,j+3,k)+3.889666722856455D6*xyzi(i+11,j+1,
     $        k)+3.237787392451808D6*xyzi(i+9,j+7,k)-1.253337055142636D7
     $        *xyzi(i+9,j+5,k)+1.166900016856936D7*xyzi(i+9,j+3,k)-3.073
     $        316916824853D6*xyzi(i+9,j+1,k)-3.237787392451808D6*xyzi(i+
     $        7,j+9,k)+7.77933344571291D6*xyzi(i+7,j+5,k)-6.146633833649
     $        707D6*xyzi(i+7,j+3,k)+1.290793105066438D6*xyzi(i+7,j+1,k)-
     $        5.828017306413255D6*xyzi(i+5,j+11,k)+1.253337055142636D7*x
     $        yzi(i+5,j+9,k)-7.77933344571291D6*xyzi(i+5,j+7,k)+1.290793
     $        105066438D6*xyzi(i+5,j+3,k)-2.693829088834306D5*xyzi(i+5,j
     $        +1,k)-3.237787392451808D6*xyzi(i+3,j+13,k)+1.0026696441141
     $        08D7*xyzi(i+3,j+11,k)-1.166900016856936D7*xyzi(i+3,j+9,k)+
     $        6.146633833649707D6*xyzi(i+3,j+7,k)-1.290793105066438D6*xy
     $        zi(i+3,j+5,k)+2.137959594312941D4*xyzi(i+3,j+1,k)-6.475574
     $        784903617D5*xyzi(i+1,j+15,k)+2.506674110285271D6*xyzi(i+1,
     $        j+13,k)-3.889666722856455D6*xyzi(i+1,j+11,k)+3.07331691682
     $        4853D6*xyzi(i+1,j+9,k)-1.290793105066438D6*xyzi(i+1,j+7,k)
     $        +2.693829088834306D5*xyzi(i+1,j+5,k)-2.137959594312941D4*x
     $        yzi(i+1,j+3,k))+angi
            angi=yrk(271)*(-4.678476280266849D5*xyzi(i+15,j+1,k)-3.27493
     $        3396186794D6*xyzi(i+13,j+3,k)+1.901574230043945D6*xyzi(i+1
     $        3,j+1,k)-9.824800188560383D6*xyzi(i+11,j+5,k)+1.1409445380
     $        26367D7*xyzi(i+11,j+3,k)-3.147433208348599D6*xyzi(i+11,j+1
     $        ,k)-1.637466698093397D7*xyzi(i+9,j+7,k)+2.852361345065918D
     $        7*xyzi(i+9,j+5,k)-1.573716604174299D7*xyzi(i+9,j+3,k)+2.72
     $        0004007214838D6*xyzi(i+9,j+1,k)-1.637466698093397D7*xyzi(i
     $        +7,j+9,k)+3.80314846008789D7*xyzi(i+7,j+7,k)-3.14743320834
     $        8599D7*xyzi(i+7,j+5,k)+1.088001602885935D7*xyzi(i+7,j+3,k)
     $        -1.305601923463122D6*xyzi(i+7,j+1,k)-9.824800188560383D6*x
     $        yzi(i+5,j+11,k)+2.852361345065918D7*xyzi(i+5,j+9,k)-3.1474
     $        33208348599D7*xyzi(i+5,j+7,k)+1.632002404328903D7*xyzi(i+5
     $        ,j+5,k)-3.916805770389367D6*xyzi(i+5,j+3,k)+3.405918061208
     $        145D5*xyzi(i+5,j+1,k)-3.274933396186794D6*xyzi(i+3,j+13,k)
     $        +1.140944538026367D7*xyzi(i+3,j+11,k)-1.573716604174299D7*
     $        xyzi(i+3,j+9,k)+1.088001602885935D7*xyzi(i+3,j+7,k)-3.9168
     $        05770389367D6*xyzi(i+3,j+5,k)+6.811836122416291D5*xyzi(i+3
     $        ,j+3,k)-4.324975315819867D4*xyzi(i+3,j+1,k)-4.678476280266
     $        849D5*xyzi(i+1,j+15,k)+1.901574230043945D6*xyzi(i+1,j+13,k
     $        )-3.147433208348599D6*xyzi(i+1,j+11,k)+2.720004007214838D6
     $        *xyzi(i+1,j+9,k)-1.305601923463122D6*xyzi(i+1,j+7,k)+3.405
     $        918061208145D5*xyzi(i+1,j+5,k)-4.324975315819867D4*xyzi(i+
     $        1,j+3,k)+1.951116683828512D3*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(274)*(2.562506993455015D5*xyzi(i+15,j,k+1)+1.793754
     $        89541851D6*xyzi(i+13,j+2,k+1)-9.258089782805215D5*xyzi(i+1
     $        3,j,k+1)+5.381264686255531D6*xyzi(i+11,j+4,k+1)-5.55485386
     $        9683129D6*xyzi(i+11,j+2,k+1)+1.34082679613041D6*xyzi(i+11,
     $        j,k+1)+8.968774477092552D6*xyzi(i+9,j+6,k+1)-1.38871346742
     $        0782D7*xyzi(i+9,j+4,k+1)+6.704133980652052D6*xyzi(i+9,j+2,
     $        k+1)-9.932050341706744D5*xyzi(i+9,j,k+1)+8.968774477092552
     $        D6*xyzi(i+7,j+8,k+1)-1.851617956561043D7*xyzi(i+7,j+6,k+1)
     $        +1.34082679613041D7*xyzi(i+7,j+4,k+1)-3.972820136682698D6*
     $        xyzi(i+7,j+2,k+1)+3.972820136682698D5*xyzi(i+7,j,k+1)+5.38
     $        1264686255531D6*xyzi(i+5,j+10,k+1)-1.388713467420782D7*xyz
     $        i(i+5,j+8,k+1)+1.34082679613041D7*xyzi(i+5,j+6,k+1)-5.9592
     $        30205024046D6*xyzi(i+5,j+4,k+1)+1.191846041004809D6*xyzi(i
     $        +5,j+2,k+1)-8.291102893946499D4*xyzi(i+5,j,k+1)+1.79375489
     $        541851D6*xyzi(i+3,j+12,k+1)-5.554853869683129D6*xyzi(i+3,j
     $        +10,k+1)+6.704133980652052D6*xyzi(i+3,j+8,k+1)-3.972820136
     $        682698D6*xyzi(i+3,j+6,k+1)+1.191846041004809D6*xyzi(i+3,j+
     $        4,k+1)-1.6582205787893D5*xyzi(i+3,j+2,k+1)+7.8962884704252
     $        37D3*xyzi(i+3,j,k+1)+2.562506993455015D5*xyzi(i+1,j+14,k+1
     $        )-9.258089782805215D5*xyzi(i+1,j+12,k+1)+1.34082679613041D
     $        6*xyzi(i+1,j+10,k+1)-9.932050341706744D5*xyzi(i+1,j+8,k+1)
     $        +3.972820136682698D5*xyzi(i+1,j+6,k+1)-8.291102893946499D4
     $        *xyzi(i+1,j+4,k+1)+7.896288470425237D3*xyzi(i+1,j+2,k+1)-2
     $        .374823600127891D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(276)*(-2.007990499395227D5*xyzi(i+15,j,k+1)-6.02397
     $        1498185681D5*xyzi(i+13,j+2,k+1)+6.995579804344662D5*xyzi(i
     $        +13,j,k+1)+6.023971498185681D5*xyzi(i+11,j+4,k+1)+1.399115
     $        960868932D6*xyzi(i+11,j+2,k+1)-9.649075592199533D5*xyzi(i+
     $        11,j,k+1)+5.019976248488067D6*xyzi(i+9,j+6,k+1)-3.49778990
     $        2172331D6*xyzi(i+9,j+4,k+1)-9.649075592199533D5*xyzi(i+9,j
     $        +2,k+1)+6.670965841520665D5*xyzi(i+9,j,k+1)+9.035957247278
     $        521D6*xyzi(i+7,j+8,k+1)-1.399115960868932D7*xyzi(i+7,j+6,k
     $        +1)+5.78944535531972D6*xyzi(i+7,j+4,k+1)-2.401547702947439
     $        D5*xyzi(i+7,j,k+1)+7.831162947641385D6*xyzi(i+5,j+10,k+1)-
     $        1.748894951086165D7*xyzi(i+5,j+8,k+1)+1.350870582907935D7*
     $        xyzi(i+5,j+6,k+1)-4.002579504912399D6*xyzi(i+5,j+4,k+1)+2.
     $        401547702947439D5*xyzi(i+5,j+2,k+1)+4.176604700778155D4*xy
     $        zi(i+5,j,k+1)+3.413583848971886D6*xyzi(i+3,j+12,k+1)-9.793
     $        811726082526D6*xyzi(i+3,j+10,k+1)+1.061398315141949D7*xyzi
     $        (i+3,j+8,k+1)-5.336772673216532D6*xyzi(i+3,j+6,k+1)+1.2007
     $        7385147372D6*xyzi(i+3,j+4,k+1)-8.353209401556311D4*xyzi(i+
     $        3,j+2,k+1)-2.651812508430575D3*xyzi(i+3,j,k+1)+6.023971498
     $        185681D5*xyzi(i+1,j+14,k+1)-2.098673941303398D6*xyzi(i+1,j
     $        +12,k+1)+2.89472267765986D6*xyzi(i+1,j+10,k+1)-2.001289752
     $        456199D6*xyzi(i+1,j+8,k+1)+7.204643108842318D5*xyzi(i+1,j+
     $        6,k+1)-1.252981410233447D5*xyzi(i+1,j+4,k+1)+7.95543752529
     $        1725D3*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(278)*(1.223768605503968D5*xyzi(i+15,j,k+1)-6.118843
     $        027519839D5*xyzi(i+13,j+2,k+1)-3.947640662916025D5*xyzi(i+
     $        13,j,k+1)-4.283190119263887D6*xyzi(i+11,j+4,k+1)+2.3685843
     $        97749615D6*xyzi(i+11,j+2,k+1)+4.900519443619893D5*xyzi(i+1
     $        1,j,k+1)-7.95449593577579D6*xyzi(i+9,j+6,k+1)+1.1448157922
     $        45647D7*xyzi(i+9,j+4,k+1)-3.430363610533925D6*xyzi(i+9,j+2
     $        ,k+1)-2.904011522145122D5*xyzi(i+9,j,k+1)-5.50695872476785
     $        5D6*xyzi(i+7,j+8,k+1)+1.421150638649769D7*xyzi(i+7,j+6,k+1
     $        )-1.078114277596376D7*xyzi(i+7,j+4,k+1)+2.323209217716097D
     $        6*xyzi(i+7,j+2,k+1)+8.131232262006341D4*xyzi(i+7,j,k+1)+1.
     $        223768605503968D5*xyzi(i+5,j+10,k+1)+3.552876596624422D6*x
     $        yzi(i+5,j+8,k+1)-6.86072722106785D6*xyzi(i+5,j+6,k+1)+4.06
     $        561613100317D6*xyzi(i+5,j+4,k+1)-7.318109035805707D5*xyzi(
     $        i+5,j+2,k+1)-8.484764099484878D3*xyzi(i+5,j,k+1)+1.8356529
     $        08255952D6*xyzi(i+3,j+12,k+1)-3.947640662916025D6*xyzi(i+3
     $        ,j+10,k+1)+2.450259721809947D6*xyzi(i+3,j+8,k+1)-4.0656161
     $        31003171D5*xyzi(i+3,j+4,k+1)+8.484764099484878D4*xyzi(i+3,
     $        j+2,k+1)+6.118843027519839D5*xyzi(i+1,j+14,k+1)-1.97382033
     $        1458012D6*xyzi(i+1,j+12,k+1)+2.450259721809947D6*xyzi(i+1,
     $        j+10,k+1)-1.452005761072561D6*xyzi(i+1,j+8,k+1)+4.06561613
     $        1003171D5*xyzi(i+1,j+6,k+1)-4.242382049742439D4*xyzi(i+1,j
     $        +4,k+1))+angi
            angi=yrk(280)*(-5.705850544029968D4*xyzi(i+15,j,k+1)+9.69994
     $        5924850946D5*xyzi(i+13,j+2,k+1)+1.619725315724636D5*xyzi(i
     $        +13,j,k+1)+2.453515733932886D6*xyzi(i+11,j+4,k+1)-2.915505
     $        568304345D6*xyzi(i+11,j+2,k+1)-1.675577912818589D5*xyzi(i+
     $        11,j,k+1)-6.276435598432965D5*xyzi(i+9,j+6,k+1)-4.04931328
     $        931159D6*xyzi(i+9,j+4,k+1)+3.183598034355319D6*xyzi(i+9,j+
     $        2,k+1)+7.447012945860396D4*xyzi(i+9,j,k+1)-5.6487920385896
     $        68D6*xyzi(i+7,j+8,k+1)+5.83101113660869D6*xyzi(i+7,j+6,k+1
     $        )+1.005346747691153D6*xyzi(i+7,j+4,k+1)-1.489402589172079D
     $        6*xyzi(i+7,j+2,k+1)-1.191522071337663D4*xyzi(i+7,j,k+1)-4.
     $        393504918903075D6*xyzi(i+5,j+10,k+1)+1.020426948906521D7*x
     $        yzi(i+5,j+8,k+1)-7.037427233838074D6*xyzi(i+5,j+6,k+1)+1.0
     $        42581812420455D6*xyzi(i+5,j+4,k+1)+2.502196349809093D5*xyz
     $        i(i+5,j+2,k+1)-3.994095380820978D5*xyzi(i+3,j+12,k+1)+2.26
     $        7615442014491D6*xyzi(i+3,j+10,k+1)-3.518713616919037D6*xyz
     $        i(i+3,j+8,k+1)+2.085163624840911D6*xyzi(i+3,j+6,k+1)-4.170
     $        327249681822D5*xyzi(i+3,j+4,k+1)+3.994095380820978D5*xyzi(
     $        i+1,j+14,k+1)-1.133807721007245D6*xyzi(i+1,j+12,k+1)+1.172
     $        904538973012D6*xyzi(i+1,j+10,k+1)-5.212909062102277D5*xyzi
     $        (i+1,j+8,k+1)+8.340654499363643D4*xyzi(i+1,j+6,k+1))+angi
            angi=yrk(282)*(1.976564608530885D4*xyzi(i+15,j,k+1)-6.522663
     $        20815192D5*xyzi(i+13,j+2,k+1)-4.590730703684636D4*xyzi(i+1
     $        3,j,k+1)+4.150785677914858D5*xyzi(i+11,j+4,k+1)+1.56084843
     $        9252776D6*xyzi(i+11,j+2,k+1)+3.482623292450413D4*xyzi(i+11
     $        ,j,k+1)+3.696175817952755D6*xyzi(i+9,j+6,k+1)-2.5249018870
     $        2655D6*xyzi(i+9,j+4,k+1)-1.218918152357645D6*xyzi(i+9,j+2,
     $        k+1)-8.599069857902255D3*xyzi(i+9,j,k+1)+1.956798962445576
     $        D6*xyzi(i+7,j+8,k+1)-6.059764528863719D6*xyzi(i+7,j+6,k+1)
     $        +3.134360963205372D6*xyzi(i+7,j+4,k+1)+3.095665148844812D5
     $        *xyzi(i+7,j+2,k+1)-1.956798962445576D6*xyzi(i+5,j+10,k+1)+
     $        1.51494113221593D6*xyzi(i+5,j+8,k+1)+1.462701782829174D6*x
     $        yzi(i+5,j+6,k+1)-1.083482802095684D6*xyzi(i+5,j+4,k+1)-1.1
     $        26641826862604D6*xyzi(i+3,j+12,k+1)+3.029882264431859D6*xy
     $        zi(i+3,j+10,k+1)-2.61196746933781D6*xyzi(i+3,j+8,k+1)+7.22
     $        3218680637894D5*xyzi(i+3,j+6,k+1)+1.778908147677796D5*xyzi
     $        (i+1,j+14,k+1)-4.131657633316172D5*xyzi(i+1,j+12,k+1)+3.13
     $        4360963205372D5*xyzi(i+1,j+10,k+1)-7.73916287211203D4*xyzi
     $        (i+1,j+8,k+1))+angi
            angi=yrk(284)*(-4.834672985156282D3*xyzi(i+15,j,k+1)+2.56237
     $        668213283D5*xyzi(i+13,j+2,k+1)+8.109774039616989D3*xyzi(i+
     $        13,j,k+1)-1.068462729719538D6*xyzi(i+11,j+4,k+1)-4.3792779
     $        81393174D5*xyzi(i+11,j+2,k+1)-3.355768568117375D3*xyzi(i+1
     $        1,j,k+1)-6.913582368773484D5*xyzi(i+9,j+6,k+1)+2.230187860
     $        894672D6*xyzi(i+9,j+4,k+1)+1.845672712464556D5*xyzi(i+9,j+
     $        2,k+1)+2.074074710632045D6*xyzi(i+7,j+8,k+1)-1.07049017322
     $        9443D6*xyzi(i+7,j+6,k+1)-1.107403627478734D6*xyzi(i+7,j+4,
     $        k+1)+6.913582368773484D5*xyzi(i+5,j+10,k+1)-2.408602889766
     $        246D6*xyzi(i+5,j+8,k+1)+1.550365078470227D6*xyzi(i+5,j+6,k
     $        +1)-6.913582368773484D5*xyzi(i+3,j+12,k+1)+1.2489052021010
     $        16D6*xyzi(i+3,j+10,k+1)-5.537018137393669D5*xyzi(i+3,j+8,k
     $        +1)+5.31814028367191D4*xyzi(i+1,j+14,k+1)-8.92075144357868
     $        8D4*xyzi(i+1,j+12,k+1)+3.691345424929112D4*xyzi(i+1,j+10,k
     $        +1))+angi
            angi=yrk(286)*(7.58759404765566D2*xyzi(i+15,j,k+1)-5.8424474
     $        16694858D4*xyzi(i+13,j+2,k+1)-6.853310752721241D2*xyzi(i+1
     $        3,j,k+1)+4.833297408356655D5*xyzi(i+11,j+4,k+1)+5.34558238
     $        7122568D4*xyzi(i+11,j+2,k+1)-7.595181641703315D5*xyzi(i+9,
     $        j+6,k+1)-4.900117188195687D5*xyzi(i+9,j+4,k+1)-3.255077846
     $        444278D5*xyzi(i+7,j+8,k+1)+1.176028125166965D6*xyzi(i+7,j+
     $        6,k+1)+7.595181641703315D5*xyzi(i+5,j+10,k+1)-8.8202109387
     $        52237D5*xyzi(i+5,j+8,k+1)-2.071413175009995D5*xyzi(i+3,j+1
     $        2,k+1)+1.960046875278275D5*xyzi(i+3,j+10,k+1)+9.8638722619
     $        52358D3*xyzi(i+1,j+14,k+1)-8.909303978537614D3*xyzi(i+1,j+
     $        12,k+1))+angi
            angi=yrk(288)*(-6.094502198942574D1*xyzi(i+15,j,k+1)+6.39922
     $        7308889703D3*xyzi(i+13,j+2,k+1)-8.318995501556614D4*xyzi(i
     $        +11,j+4,k+1)+3.050298350570758D5*xyzi(i+9,j+6,k+1)-3.92181
     $        2165019546D5*xyzi(i+7,j+8,k+1)+1.830179010342455D5*xyzi(i+
     $        5,j+10,k+1)-2.772998500518871D4*xyzi(i+3,j+12,k+1)+9.14175
     $        3298413861D2*xyzi(i+1,j+14,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(258)*(-9.141753298413861D2*xyzi(i+14,j+1,k+1)+2.772
     $        998500518871D4*xyzi(i+12,j+3,k+1)-1.830179010342455D5*xyzi
     $        (i+10,j+5,k+1)+3.921812165019546D5*xyzi(i+8,j+7,k+1)-3.050
     $        298350570758D5*xyzi(i+6,j+9,k+1)+8.318995501556614D4*xyzi(
     $        i+4,j+11,k+1)-6.399227308889703D3*xyzi(i+2,j+13,k+1)+6.094
     $        502198942574D1*xyzi(i,j+15,k+1))+angi
            angi=yrk(260)*(9.863872261952358D3*xyzi(i+14,j+1,k+1)-2.0714
     $        13175009995D5*xyzi(i+12,j+3,k+1)-8.909303978537614D3*xyzi(
     $        i+12,j+1,k+1)+7.595181641703315D5*xyzi(i+10,j+5,k+1)+1.960
     $        046875278275D5*xyzi(i+10,j+3,k+1)-3.255077846444278D5*xyzi
     $        (i+8,j+7,k+1)-8.820210938752237D5*xyzi(i+8,j+5,k+1)-7.5951
     $        81641703315D5*xyzi(i+6,j+9,k+1)+1.176028125166965D6*xyzi(i
     $        +6,j+7,k+1)+4.833297408356655D5*xyzi(i+4,j+11,k+1)-4.90011
     $        7188195687D5*xyzi(i+4,j+9,k+1)-5.842447416694858D4*xyzi(i+
     $        2,j+13,k+1)+5.345582387122568D4*xyzi(i+2,j+11,k+1)+7.58759
     $        404765566D2*xyzi(i,j+15,k+1)-6.853310752721241D2*xyzi(i,j+
     $        13,k+1))+angi
            angi=yrk(262)*(-5.31814028367191D4*xyzi(i+14,j+1,k+1)+6.9135
     $        82368773484D5*xyzi(i+12,j+3,k+1)+8.920751443578688D4*xyzi(
     $        i+12,j+1,k+1)-6.913582368773484D5*xyzi(i+10,j+5,k+1)-1.248
     $        905202101016D6*xyzi(i+10,j+3,k+1)-3.691345424929112D4*xyzi
     $        (i+10,j+1,k+1)-2.074074710632045D6*xyzi(i+8,j+7,k+1)+2.408
     $        602889766246D6*xyzi(i+8,j+5,k+1)+5.537018137393669D5*xyzi(
     $        i+8,j+3,k+1)+6.913582368773484D5*xyzi(i+6,j+9,k+1)+1.07049
     $        0173229443D6*xyzi(i+6,j+7,k+1)-1.550365078470227D6*xyzi(i+
     $        6,j+5,k+1)+1.068462729719538D6*xyzi(i+4,j+11,k+1)-2.230187
     $        860894672D6*xyzi(i+4,j+9,k+1)+1.107403627478734D6*xyzi(i+4
     $        ,j+7,k+1)-2.56237668213283D5*xyzi(i+2,j+13,k+1)+4.37927798
     $        1393174D5*xyzi(i+2,j+11,k+1)-1.845672712464556D5*xyzi(i+2,
     $        j+9,k+1)+4.834672985156282D3*xyzi(i,j+15,k+1)-8.1097740396
     $        16989D3*xyzi(i,j+13,k+1)+3.355768568117375D3*xyzi(i,j+11,k
     $        +1))+angi
            angi=yrk(264)*(1.778908147677796D5*xyzi(i+14,j+1,k+1)-1.1266
     $        41826862604D6*xyzi(i+12,j+3,k+1)-4.131657633316172D5*xyzi(
     $        i+12,j+1,k+1)-1.956798962445576D6*xyzi(i+10,j+5,k+1)+3.029
     $        882264431859D6*xyzi(i+10,j+3,k+1)+3.134360963205372D5*xyzi
     $        (i+10,j+1,k+1)+1.956798962445576D6*xyzi(i+8,j+7,k+1)+1.514
     $        94113221593D6*xyzi(i+8,j+5,k+1)-2.61196746933781D6*xyzi(i+
     $        8,j+3,k+1)-7.73916287211203D4*xyzi(i+8,j+1,k+1)+3.69617581
     $        7952755D6*xyzi(i+6,j+9,k+1)-6.059764528863719D6*xyzi(i+6,j
     $        +7,k+1)+1.462701782829174D6*xyzi(i+6,j+5,k+1)+7.2232186806
     $        37894D5*xyzi(i+6,j+3,k+1)+4.150785677914858D5*xyzi(i+4,j+1
     $        1,k+1)-2.52490188702655D6*xyzi(i+4,j+9,k+1)+3.134360963205
     $        372D6*xyzi(i+4,j+7,k+1)-1.083482802095684D6*xyzi(i+4,j+5,k
     $        +1)-6.52266320815192D5*xyzi(i+2,j+13,k+1)+1.56084843925277
     $        6D6*xyzi(i+2,j+11,k+1)-1.218918152357645D6*xyzi(i+2,j+9,k+
     $        1)+3.095665148844812D5*xyzi(i+2,j+7,k+1)+1.976564608530885
     $        D4*xyzi(i,j+15,k+1)-4.590730703684636D4*xyzi(i,j+13,k+1)+3
     $        .482623292450413D4*xyzi(i,j+11,k+1)-8.599069857902255D3*xy
     $        zi(i,j+9,k+1))+angi
            angi=yrk(266)*(-3.994095380820978D5*xyzi(i+14,j+1,k+1)+3.994
     $        095380820978D5*xyzi(i+12,j+3,k+1)+1.133807721007245D6*xyzi
     $        (i+12,j+1,k+1)+4.393504918903075D6*xyzi(i+10,j+5,k+1)-2.26
     $        7615442014491D6*xyzi(i+10,j+3,k+1)-1.172904538973012D6*xyz
     $        i(i+10,j+1,k+1)+5.648792038589668D6*xyzi(i+8,j+7,k+1)-1.02
     $        0426948906521D7*xyzi(i+8,j+5,k+1)+3.518713616919037D6*xyzi
     $        (i+8,j+3,k+1)+5.212909062102277D5*xyzi(i+8,j+1,k+1)+6.2764
     $        35598432965D5*xyzi(i+6,j+9,k+1)-5.83101113660869D6*xyzi(i+
     $        6,j+7,k+1)+7.037427233838074D6*xyzi(i+6,j+5,k+1)-2.0851636
     $        24840911D6*xyzi(i+6,j+3,k+1)-8.340654499363643D4*xyzi(i+6,
     $        j+1,k+1)-2.453515733932886D6*xyzi(i+4,j+11,k+1)+4.04931328
     $        931159D6*xyzi(i+4,j+9,k+1)-1.005346747691153D6*xyzi(i+4,j+
     $        7,k+1)-1.042581812420455D6*xyzi(i+4,j+5,k+1)+4.17032724968
     $        1822D5*xyzi(i+4,j+3,k+1)-9.699945924850946D5*xyzi(i+2,j+13
     $        ,k+1)+2.915505568304345D6*xyzi(i+2,j+11,k+1)-3.18359803435
     $        5319D6*xyzi(i+2,j+9,k+1)+1.489402589172079D6*xyzi(i+2,j+7,
     $        k+1)-2.502196349809093D5*xyzi(i+2,j+5,k+1)+5.7058505440299
     $        68D4*xyzi(i,j+15,k+1)-1.619725315724636D5*xyzi(i,j+13,k+1)
     $        +1.675577912818589D5*xyzi(i,j+11,k+1)-7.447012945860396D4*
     $        xyzi(i,j+9,k+1)+1.191522071337663D4*xyzi(i,j+7,k+1))+angi
            angi=yrk(268)*(6.118843027519839D5*xyzi(i+14,j+1,k+1)+1.8356
     $        52908255952D6*xyzi(i+12,j+3,k+1)-1.973820331458012D6*xyzi(
     $        i+12,j+1,k+1)+1.223768605503968D5*xyzi(i+10,j+5,k+1)-3.947
     $        640662916025D6*xyzi(i+10,j+3,k+1)+2.450259721809947D6*xyzi
     $        (i+10,j+1,k+1)-5.506958724767855D6*xyzi(i+8,j+7,k+1)+3.552
     $        876596624422D6*xyzi(i+8,j+5,k+1)+2.450259721809947D6*xyzi(
     $        i+8,j+3,k+1)-1.452005761072561D6*xyzi(i+8,j+1,k+1)-7.95449
     $        593577579D6*xyzi(i+6,j+9,k+1)+1.421150638649769D7*xyzi(i+6
     $        ,j+7,k+1)-6.86072722106785D6*xyzi(i+6,j+5,k+1)+4.065616131
     $        003171D5*xyzi(i+6,j+1,k+1)-4.283190119263887D6*xyzi(i+4,j+
     $        11,k+1)+1.144815792245647D7*xyzi(i+4,j+9,k+1)-1.0781142775
     $        96376D7*xyzi(i+4,j+7,k+1)+4.06561613100317D6*xyzi(i+4,j+5,
     $        k+1)-4.065616131003171D5*xyzi(i+4,j+3,k+1)-4.2423820497424
     $        39D4*xyzi(i+4,j+1,k+1)-6.118843027519839D5*xyzi(i+2,j+13,k
     $        +1)+2.368584397749615D6*xyzi(i+2,j+11,k+1)-3.4303636105339
     $        25D6*xyzi(i+2,j+9,k+1)+2.323209217716097D6*xyzi(i+2,j+7,k+
     $        1)-7.318109035805707D5*xyzi(i+2,j+5,k+1)+8.484764099484878
     $        D4*xyzi(i+2,j+3,k+1)+1.223768605503968D5*xyzi(i,j+15,k+1)-
     $        3.947640662916025D5*xyzi(i,j+13,k+1)+4.900519443619893D5*x
     $        yzi(i,j+11,k+1)-2.904011522145122D5*xyzi(i,j+9,k+1)+8.1312
     $        32262006341D4*xyzi(i,j+7,k+1)-8.484764099484878D3*xyzi(i,j
     $        +5,k+1))+angi
            angi=yrk(270)*(-6.023971498185681D5*xyzi(i+14,j+1,k+1)-3.413
     $        583848971886D6*xyzi(i+12,j+3,k+1)+2.098673941303398D6*xyzi
     $        (i+12,j+1,k+1)-7.831162947641385D6*xyzi(i+10,j+5,k+1)+9.79
     $        3811726082526D6*xyzi(i+10,j+3,k+1)-2.89472267765986D6*xyzi
     $        (i+10,j+1,k+1)-9.035957247278521D6*xyzi(i+8,j+7,k+1)+1.748
     $        894951086165D7*xyzi(i+8,j+5,k+1)-1.061398315141949D7*xyzi(
     $        i+8,j+3,k+1)+2.001289752456199D6*xyzi(i+8,j+1,k+1)-5.01997
     $        6248488067D6*xyzi(i+6,j+9,k+1)+1.399115960868932D7*xyzi(i+
     $        6,j+7,k+1)-1.350870582907935D7*xyzi(i+6,j+5,k+1)+5.3367726
     $        73216532D6*xyzi(i+6,j+3,k+1)-7.204643108842318D5*xyzi(i+6,
     $        j+1,k+1)-6.023971498185681D5*xyzi(i+4,j+11,k+1)+3.49778990
     $        2172331D6*xyzi(i+4,j+9,k+1)-5.78944535531972D6*xyzi(i+4,j+
     $        7,k+1)+4.002579504912399D6*xyzi(i+4,j+5,k+1)-1.20077385147
     $        372D6*xyzi(i+4,j+3,k+1)+1.252981410233447D5*xyzi(i+4,j+1,k
     $        +1)+6.023971498185681D5*xyzi(i+2,j+13,k+1)-1.3991159608689
     $        32D6*xyzi(i+2,j+11,k+1)+9.649075592199533D5*xyzi(i+2,j+9,k
     $        +1)-2.401547702947439D5*xyzi(i+2,j+5,k+1)+8.35320940155631
     $        1D4*xyzi(i+2,j+3,k+1)-7.955437525291725D3*xyzi(i+2,j+1,k+1
     $        )+2.007990499395227D5*xyzi(i,j+15,k+1)-6.995579804344662D5
     $        *xyzi(i,j+13,k+1)+9.649075592199533D5*xyzi(i,j+11,k+1)-6.6
     $        70965841520665D5*xyzi(i,j+9,k+1)+2.401547702947439D5*xyzi(
     $        i,j+7,k+1)-4.176604700778155D4*xyzi(i,j+5,k+1)+2.651812508
     $        430575D3*xyzi(i,j+3,k+1))+angi
            angi=yrk(272)*(2.562506993455015D5*xyzi(i+14,j+1,k+1)+1.7937
     $        5489541851D6*xyzi(i+12,j+3,k+1)-9.258089782805215D5*xyzi(i
     $        +12,j+1,k+1)+5.381264686255531D6*xyzi(i+10,j+5,k+1)-5.5548
     $        53869683129D6*xyzi(i+10,j+3,k+1)+1.34082679613041D6*xyzi(i
     $        +10,j+1,k+1)+8.968774477092552D6*xyzi(i+8,j+7,k+1)-1.38871
     $        3467420782D7*xyzi(i+8,j+5,k+1)+6.704133980652052D6*xyzi(i+
     $        8,j+3,k+1)-9.932050341706744D5*xyzi(i+8,j+1,k+1)+8.9687744
     $        77092552D6*xyzi(i+6,j+9,k+1)-1.851617956561043D7*xyzi(i+6,
     $        j+7,k+1)+1.34082679613041D7*xyzi(i+6,j+5,k+1)-3.9728201366
     $        82698D6*xyzi(i+6,j+3,k+1)+3.972820136682698D5*xyzi(i+6,j+1
     $        ,k+1)+5.381264686255531D6*xyzi(i+4,j+11,k+1)-1.38871346742
     $        0782D7*xyzi(i+4,j+9,k+1)+1.34082679613041D7*xyzi(i+4,j+7,k
     $        +1)-5.959230205024046D6*xyzi(i+4,j+5,k+1)+1.19184604100480
     $        9D6*xyzi(i+4,j+3,k+1)-8.291102893946499D4*xyzi(i+4,j+1,k+1
     $        )+1.79375489541851D6*xyzi(i+2,j+13,k+1)-5.554853869683129D
     $        6*xyzi(i+2,j+11,k+1)+6.704133980652052D6*xyzi(i+2,j+9,k+1)
     $        -3.972820136682698D6*xyzi(i+2,j+7,k+1)+1.191846041004809D6
     $        *xyzi(i+2,j+5,k+1)-1.6582205787893D5*xyzi(i+2,j+3,k+1)+7.8
     $        96288470425237D3*xyzi(i+2,j+1,k+1)+2.562506993455015D5*xyz
     $        i(i,j+15,k+1)-9.258089782805215D5*xyzi(i,j+13,k+1)+1.34082
     $        679613041D6*xyzi(i,j+11,k+1)-9.932050341706744D5*xyzi(i,j+
     $        9,k+1)+3.972820136682698D5*xyzi(i,j+7,k+1)-8.2911028939464
     $        99D4*xyzi(i,j+5,k+1)+7.896288470425237D3*xyzi(i,j+3,k+1)-2
     $        .374823600127891D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA(16)=angi
        endif
      else
c...
c...  l=1
c...
        if(lmax.ge.1)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=angi-6.139960247678931D0*yrk(4)*xyzi(i+1,j,k)
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=angi-6.139960247678931D0*yrk(2)*xyzi(i,j+1,k)
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=6.139960247678931D0*yrk(3)*xyzi(i,j,k+1)+angi
          endif
          OMEGA(1)=angi
        endif
c...
c...  l=3
c...
        if(lmax.ge.3)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(14)*(2.871703451903279D1*xyzi(i+3,j,k)+2.8717034519
     $        03279D1*xyzi(i+1,j+2,k)-2.297362761522623D1*xyzi(i+1,j,k))
     $        +angi
            angi=yrk(16)*(2.22441192889355D1*xyzi(i+1,j+2,k)-7.414706429
     $        645167D0*xyzi(i+3,j,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(10)*(7.414706429645167D0*xyzi(i,j+3,k)-2.2244119288
     $        9355D1*xyzi(i+2,j+1,k))+angi
            angi=yrk(12)*(2.871703451903279D1*xyzi(i+2,j+1,k)+2.87170345
     $        1903279D1*xyzi(i,j+3,k)-2.297362761522623D1*xyzi(i,j+1,k))
     $        +angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(13)*(2.344736049917376D1*xyzi(i,j,k+3)-1.4068416299
     $        50425D1*xyzi(i,j,k+1))+angi
            angi=yrk(15)*(1.816224734516432D1*xyzi(i+2,j,k+1)-1.81622473
     $        4516432D1*xyzi(i,j+2,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=3.632449469032864D1*yrk(11)*xyzi(i+1,j+1,k+1)+angi
          endif
          OMEGA(3)=angi
        endif
c...
c...  l=5
c...
        if(lmax.ge.5)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(32)*(-1.19529805236618D2*xyzi(i+5,j,k)-2.3905961047
     $        3236D2*xyzi(i+3,j+2,k)+1.59373073648824D2*xyzi(i+3,j,k)-1.
     $        19529805236618D2*xyzi(i+1,j+4,k)+1.59373073648824D2*xyzi(i
     $        +1,j+2,k)-4.5535163899664D1*xyzi(i+1,j,k))+angi
            angi=yrk(34)*(5.533154810497966D1*xyzi(i+5,j,k)-1.1066309620
     $        99593D2*xyzi(i+3,j+2,k)-4.918359831553747D1*xyzi(i+3,j,k)-
     $        1.65994644314939D2*xyzi(i+1,j+4,k)+1.475507949466124D2*xyz
     $        i(i+1,j+2,k))+angi
            angi=yrk(36)*(-8.248340190868946D0*xyzi(i+5,j,k)+8.248340190
     $        868946D1*xyzi(i+3,j+2,k)-4.124170095434473D1*xyzi(i+1,j+4,
     $        k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(26)*(-4.124170095434473D1*xyzi(i+4,j+1,k)+8.2483401
     $        90868946D1*xyzi(i+2,j+3,k)-8.248340190868946D0*xyzi(i,j+5,
     $        k))+angi
            angi=yrk(28)*(1.65994644314939D2*xyzi(i+4,j+1,k)+1.106630962
     $        099593D2*xyzi(i+2,j+3,k)-1.475507949466124D2*xyzi(i+2,j+1,
     $        k)-5.533154810497966D1*xyzi(i,j+5,k)+4.918359831553747D1*x
     $        yzi(i,j+3,k))+angi
            angi=yrk(30)*(-1.19529805236618D2*xyzi(i+4,j+1,k)-2.39059610
     $        473236D2*xyzi(i+2,j+3,k)+1.59373073648824D2*xyzi(i+2,j+1,k
     $        )-1.19529805236618D2*xyzi(i,j+5,k)+1.59373073648824D2*xyzi
     $        (i,j+3,k)-4.5535163899664D1*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(31)*(9.258738901136752D1*xyzi(i,j,k+5)-1.0287487667
     $        92972D2*xyzi(i,j,k+3)+2.204461643127798D1*xyzi(i,j,k+1))+a
     $        ngi
            angi=yrk(33)*(-9.035603969030778D1*xyzi(i+4,j,k+1)+6.0237359
     $        79353852D1*xyzi(i+2,j,k+1)+9.035603969030778D1*xyzi(i,j+4,
     $        k+1)-6.023735979353852D1*xyzi(i,j+2,k+1))+angi
            angi=yrk(35)*(2.608354191905385D1*xyzi(i+4,j,k+1)-1.56501251
     $        5143231D2*xyzi(i+2,j+2,k+1)+2.608354191905385D1*xyzi(i,j+4
     $        ,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(27)*(1.043341676762154D2*xyzi(i+3,j+1,k+1)-1.043341
     $        676762154D2*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(29)*(-1.807120793806156D2*xyzi(i+3,j+1,k+1)-1.80712
     $        0793806156D2*xyzi(i+1,j+3,k+1)+1.20474719587077D2*xyzi(i+1
     $        ,j+1,k+1))+angi
          endif
          OMEGA(5)=angi
        endif
c...
c...  l=7
c...
        if(lmax.ge.7)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(58)*(4.869752569422183D2*xyzi(i+7,j,k)+1.4609257708
     $        26655D3*xyzi(i+5,j+2,k)-8.990312435856337D2*xyzi(i+5,j,k)+
     $        1.460925770826655D3*xyzi(i+3,j+4,k)-1.798062487171267D3*xy
     $        zi(i+3,j+2,k)+4.903806783194366D2*xyzi(i+3,j,k)+4.86975256
     $        9422183D2*xyzi(i+1,j+6,k)-8.990312435856337D2*xyzi(i+1,j+4
     $        ,k)+4.903806783194366D2*xyzi(i+1,j+2,k)-7.264898938065727D
     $        1*xyzi(i+1,j,k))+angi
            angi=yrk(60)*(-2.811552956842769D2*xyzi(i+7,j,k)+2.811552956
     $        842769D2*xyzi(i+5,j+2,k)+4.325466087450414D2*xyzi(i+5,j,k)
     $        +1.405776478421384D3*xyzi(i+3,j+4,k)-8.650932174900827D2*x
     $        yzi(i+3,j+2,k)-1.572896759072878D2*xyzi(i+3,j,k)+8.4346588
     $        70528306D2*xyzi(i+1,j+6,k)-1.297639826235124D3*xyzi(i+1,j+
     $        4,k)+4.718690277218633D2*xyzi(i+1,j+2,k))+angi
            angi=yrk(62)*(8.477151123692503D1*xyzi(i+7,j,k)-7.6294360113
     $        23252D2*xyzi(i+5,j+2,k)-7.825062575716156D1*xyzi(i+5,j,k)-
     $        4.238575561846251D2*xyzi(i+3,j+4,k)+7.825062575716156D2*xy
     $        zi(i+3,j+2,k)+4.238575561846251D2*xyzi(i+1,j+6,k)-3.912531
     $        287858078D2*xyzi(i+1,j+4,k))+angi
            angi=yrk(64)*(-8.886468981567021D0*xyzi(i+7,j,k)+1.866158486
     $        129074D2*xyzi(i+5,j+2,k)-3.110264143548457D2*xyzi(i+3,j+4,
     $        k)+6.220528287096915D1*xyzi(i+1,j+6,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(50)*(-6.220528287096915D1*xyzi(i+6,j+1,k)+3.1102641
     $        43548457D2*xyzi(i+4,j+3,k)-1.866158486129074D2*xyzi(i+2,j+
     $        5,k)+8.886468981567021D0*xyzi(i,j+7,k))+angi
            angi=yrk(52)*(4.238575561846251D2*xyzi(i+6,j+1,k)-4.23857556
     $        1846251D2*xyzi(i+4,j+3,k)-3.912531287858078D2*xyzi(i+4,j+1
     $        ,k)-7.629436011323252D2*xyzi(i+2,j+5,k)+7.825062575716156D
     $        2*xyzi(i+2,j+3,k)+8.477151123692503D1*xyzi(i,j+7,k)-7.8250
     $        62575716156D1*xyzi(i,j+5,k))+angi
            angi=yrk(54)*(-8.434658870528306D2*xyzi(i+6,j+1,k)-1.4057764
     $        78421384D3*xyzi(i+4,j+3,k)+1.297639826235124D3*xyzi(i+4,j+
     $        1,k)-2.811552956842769D2*xyzi(i+2,j+5,k)+8.650932174900827
     $        D2*xyzi(i+2,j+3,k)-4.718690277218633D2*xyzi(i+2,j+1,k)+2.8
     $        11552956842769D2*xyzi(i,j+7,k)-4.325466087450414D2*xyzi(i,
     $        j+5,k)+1.572896759072878D2*xyzi(i,j+3,k))+angi
            angi=yrk(56)*(4.869752569422183D2*xyzi(i+6,j+1,k)+1.46092577
     $        0826655D3*xyzi(i+4,j+3,k)-8.990312435856337D2*xyzi(i+4,j+1
     $        ,k)+1.460925770826655D3*xyzi(i+2,j+5,k)-1.798062487171267D
     $        3*xyzi(i+2,j+3,k)+4.903806783194366D2*xyzi(i+2,j+1,k)+4.86
     $        9752569422183D2*xyzi(i,j+7,k)-8.990312435856337D2*xyzi(i,j
     $        +5,k)+4.903806783194366D2*xyzi(i,j+3,k)-7.264898938065727D
     $        1*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(57)*(3.681186927173971D2*xyzi(i,j,k+7)-5.9465327285
     $        11799D2*xyzi(i,j,k+5)+2.702969422050818D2*xyzi(i,j,k+3)-3.
     $        003299357834242D1*xyzi(i,j,k+1))+angi
            angi=yrk(59)*(3.976136322897221D2*xyzi(i+6,j,k+1)+3.97613632
     $        2897221D2*xyzi(i+4,j+2,k+1)-4.89370624356581D2*xyzi(i+4,j,
     $        k+1)-3.976136322897221D2*xyzi(i+2,j+4,k+1)+1.3346471573361
     $        3D2*xyzi(i+2,j,k+1)-3.976136322897221D2*xyzi(i,j+6,k+1)+4.
     $        89370624356581D2*xyzi(i,j+4,k+1)-1.33464715733613D2*xyzi(i
     $        ,j+2,k+1))+angi
            angi=yrk(61)*(-1.695430224738501D2*xyzi(i+6,j,k+1)+8.4771511
     $        23692502D2*xyzi(i+4,j+2,k+1)+1.304177095952693D2*xyzi(i+4,
     $        j,k+1)+8.477151123692502D2*xyzi(i+2,j+4,k+1)-7.82506257571
     $        6156D2*xyzi(i+2,j+2,k+1)-1.695430224738501D2*xyzi(i,j+6,k+
     $        1)+1.304177095952693D2*xyzi(i,j+4,k+1))+angi
            angi=yrk(63)*(3.325012230721775D1*xyzi(i+6,j,k+1)-4.98751834
     $        6082662D2*xyzi(i+4,j+2,k+1)+4.987518346082662D2*xyzi(i+2,j
     $        +4,k+1)-3.325012230721775D1*xyzi(i,j+6,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(51)*(1.995007338433065D2*xyzi(i+5,j+1,k+1)-6.650024
     $        46144355D2*xyzi(i+3,j+3,k+1)+1.995007338433065D2*xyzi(i+1,
     $        j+5,k+1))+angi
            angi=yrk(53)*(-6.781720898954002D2*xyzi(i+5,j+1,k+1)+5.21670
     $        8383810771D2*xyzi(i+3,j+1,k+1)+6.781720898954002D2*xyzi(i+
     $        1,j+5,k+1)-5.216708383810771D2*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(55)*(7.952272645794442D2*xyzi(i+5,j+1,k+1)+1.590454
     $        529158888D3*xyzi(i+3,j+3,k+1)-9.787412487131621D2*xyzi(i+3
     $        ,j+1,k+1)+7.952272645794442D2*xyzi(i+1,j+5,k+1)-9.78741248
     $        7131621D2*xyzi(i+1,j+3,k+1)+2.66929431467226D2*xyzi(i+1,j+
     $        1,k+1))+angi
          endif
          OMEGA(7)=angi
        endif
c...
c...  l=9
c...
        if(lmax.ge.9)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(92)*(-1.968624920969132D3*xyzi(i+9,j,k)-7.874499683
     $        876526D3*xyzi(i+7,j+2,k)+4.632058637574427D3*xyzi(i+7,j,k)
     $        -1.181174952581479D4*xyzi(i+5,j+4,k)+1.389617591272328D4*x
     $        yzi(i+5,j+2,k)-3.705646910059542D3*xyzi(i+5,j,k)-7.8744996
     $        83876526D3*xyzi(i+3,j+6,k)+1.389617591272328D4*xyzi(i+3,j+
     $        4,k)-7.411293820119083D3*xyzi(i+3,j+2,k)+1.14019904924909D
     $        3*xyzi(i+3,j,k)-1.968624920969132D3*xyzi(i+1,j+8,k)+4.6320
     $        58637574427D3*xyzi(i+1,j+6,k)-3.705646910059542D3*xyzi(i+1
     $        ,j+4,k)+1.14019904924909D3*xyzi(i+1,j+2,k)-1.0365445902264
     $        45D2*xyzi(i+1,j,k))+angi
            angi=yrk(94)*(1.2822420836111D3*xyzi(i+9,j,k)-2.715336177058
     $        8D3*xyzi(i+7,j,k)-7.693452501666601D3*xyzi(i+5,j+4,k)+2.71
     $        53361770588D3*xyzi(i+5,j+2,k)+1.8102241180392D3*xyzi(i+5,j
     $        ,k)-1.02579366688888D4*xyzi(i+3,j+6,k)+1.3576680885294D4*x
     $        yzi(i+3,j+4,k)-3.6204482360784D3*xyzi(i+3,j+2,k)-3.7132802
     $        42131693D2*xyzi(i+3,j,k)-3.8467262508333D3*xyzi(i+1,j+8,k)
     $        +8.146008531176401D3*xyzi(i+1,j+6,k)-5.430672354117601D3*x
     $        yzi(i+1,j+4,k)+1.113984072639508D3*xyzi(i+1,j+2,k))+angi
            angi=yrk(96)*(-5.205889671223938D2*xyzi(i+9,j,k)+4.164711736
     $        97915D3*xyzi(i+7,j+2,k)+8.574406517310015D2*xyzi(i+7,j,k)+
     $        7.288245539713513D3*xyzi(i+5,j+4,k)-7.716965865579014D3*xy
     $        zi(i+5,j+2,k)-3.429762606924006D2*xyzi(i+5,j,k)-4.28720325
     $        8655008D3*xyzi(i+3,j+4,k)+3.429762606924006D3*xyzi(i+3,j+2
     $        ,k)-2.602944835611969D3*xyzi(i+1,j+8,k)+4.287203258655008D
     $        3*xyzi(i+1,j+6,k)-1.714881303462003D3*xyzi(i+1,j+4,k))+ang
     $        i
            angi=yrk(98)*(1.164072318822076D2*xyzi(i+9,j,k)-2.3281446376
     $        44151D3*xyzi(i+7,j+2,k)-1.095597476538424D2*xyzi(i+7,j,k)+
     $        1.629701246350906D3*xyzi(i+5,j+4,k)+2.300754700730691D3*xy
     $        zi(i+5,j+2,k)+3.259402492701812D3*xyzi(i+3,j+6,k)-3.834591
     $        167884484D3*xyzi(i+3,j+4,k)-8.148506231754529D2*xyzi(i+1,j
     $        +8,k)+7.669182335768969D2*xyzi(i+1,j+6,k))+angi
            angi=yrk(100)*(-9.410966914433519D0*xyzi(i+9,j,k)+3.38794808
     $        9196067D2*xyzi(i+7,j+2,k)-1.185781831218623D3*xyzi(i+5,j+4
     $        ,k)+7.905212208124156D2*xyzi(i+3,j+6,k)-8.469870222990167D
     $        1*xyzi(i+1,j+8,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(82)*(-8.469870222990167D1*xyzi(i+8,j+1,k)+7.9052122
     $        08124156D2*xyzi(i+6,j+3,k)-1.185781831218623D3*xyzi(i+4,j+
     $        5,k)+3.387948089196067D2*xyzi(i+2,j+7,k)-9.410966914433519
     $        D0*xyzi(i,j+9,k))+angi
            angi=yrk(84)*(8.148506231754529D2*xyzi(i+8,j+1,k)-3.25940249
     $        2701812D3*xyzi(i+6,j+3,k)-7.669182335768969D2*xyzi(i+6,j+1
     $        ,k)-1.629701246350906D3*xyzi(i+4,j+5,k)+3.834591167884484D
     $        3*xyzi(i+4,j+3,k)+2.328144637644151D3*xyzi(i+2,j+7,k)-2.30
     $        0754700730691D3*xyzi(i+2,j+5,k)-1.164072318822076D2*xyzi(i
     $        ,j+9,k)+1.095597476538424D2*xyzi(i,j+7,k))+angi
            angi=yrk(86)*(-2.602944835611969D3*xyzi(i+8,j+1,k)+4.2872032
     $        58655008D3*xyzi(i+6,j+1,k)+7.288245539713513D3*xyzi(i+4,j+
     $        5,k)-4.287203258655008D3*xyzi(i+4,j+3,k)-1.714881303462003
     $        D3*xyzi(i+4,j+1,k)+4.16471173697915D3*xyzi(i+2,j+7,k)-7.71
     $        6965865579014D3*xyzi(i+2,j+5,k)+3.429762606924006D3*xyzi(i
     $        +2,j+3,k)-5.205889671223938D2*xyzi(i,j+9,k)+8.574406517310
     $        015D2*xyzi(i,j+7,k)-3.429762606924006D2*xyzi(i,j+5,k))+ang
     $        i
            angi=yrk(88)*(3.8467262508333D3*xyzi(i+8,j+1,k)+1.0257936668
     $        8888D4*xyzi(i+6,j+3,k)-8.146008531176401D3*xyzi(i+6,j+1,k)
     $        +7.693452501666601D3*xyzi(i+4,j+5,k)-1.3576680885294D4*xyz
     $        i(i+4,j+3,k)+5.430672354117601D3*xyzi(i+4,j+1,k)-2.7153361
     $        770588D3*xyzi(i+2,j+5,k)+3.6204482360784D3*xyzi(i+2,j+3,k)
     $        -1.113984072639508D3*xyzi(i+2,j+1,k)-1.2822420836111D3*xyz
     $        i(i,j+9,k)+2.7153361770588D3*xyzi(i,j+7,k)-1.8102241180392
     $        D3*xyzi(i,j+5,k)+3.713280242131693D2*xyzi(i,j+3,k))+angi
            angi=yrk(90)*(-1.968624920969132D3*xyzi(i+8,j+1,k)-7.8744996
     $        83876526D3*xyzi(i+6,j+3,k)+4.632058637574427D3*xyzi(i+6,j+
     $        1,k)-1.181174952581479D4*xyzi(i+4,j+5,k)+1.389617591272328
     $        D4*xyzi(i+4,j+3,k)-3.705646910059542D3*xyzi(i+4,j+1,k)-7.8
     $        74499683876526D3*xyzi(i+2,j+7,k)+1.389617591272328D4*xyzi(
     $        i+2,j+5,k)-7.411293820119083D3*xyzi(i+2,j+3,k)+1.140199049
     $        24909D3*xyzi(i+2,j+1,k)-1.968624920969132D3*xyzi(i,j+9,k)+
     $        4.632058637574427D3*xyzi(i,j+7,k)-3.705646910059542D3*xyzi
     $        (i,j+5,k)+1.14019904924909D3*xyzi(i,j+3,k)-1.0365445902264
     $        45D2*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(91)*(1.467326381829043D3*xyzi(i,j,k+9)-3.1072793968
     $        14444D3*xyzi(i,j,k+7)+2.175095577770111D3*xyzi(i,j,k+5)-5.
     $        57716814812849D2*xyzi(i,j,k+3)+3.802614646451243D1*xyzi(i,
     $        j,k+1))+angi
            angi=yrk(93)*(-1.678848973544503D3*xyzi(i+8,j,k+1)-3.3576979
     $        47089007D3*xyzi(i+6,j+2,k+1)+2.962674659196182D3*xyzi(i+6,
     $        j,k+1)+2.962674659196182D3*xyzi(i+4,j+2,k+1)-1.58009315157
     $        1297D3*xyzi(i+4,j,k+1)+3.357697947089007D3*xyzi(i+2,j+6,k+
     $        1)-2.962674659196182D3*xyzi(i+2,j+4,k+1)+2.430912540878919
     $        D2*xyzi(i+2,j,k+1)+1.678848973544503D3*xyzi(i,j+8,k+1)-2.9
     $        62674659196182D3*xyzi(i,j+6,k+1)+1.580093151571297D3*xyzi(
     $        i,j+4,k+1)-2.430912540878919D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(95)*(8.711119580919379D2*xyzi(i+8,j,k+1)-3.48444783
     $        2367752D3*xyzi(i+6,j+2,k+1)-1.229805117306265D3*xyzi(i+6,j
     $        ,k+1)-8.711119580919379D3*xyzi(i+4,j+4,k+1)+6.149025586531
     $        326D3*xyzi(i+4,j+2,k+1)+4.099350391020884D2*xyzi(i+4,j,k+1
     $        )-3.484447832367752D3*xyzi(i+2,j+6,k+1)+6.149025586531326D
     $        3*xyzi(i+2,j+4,k+1)-2.459610234612531D3*xyzi(i+2,j+2,k+1)+
     $        8.711119580919379D2*xyzi(i,j+8,k+1)-1.229805117306265D3*xy
     $        zi(i,j+6,k+1)+4.099350391020884D2*xyzi(i,j+4,k+1))+angi
            angi=yrk(97)*(-2.688309866512469D2*xyzi(i+8,j,k+1)+3.7636338
     $        13117456D3*xyzi(i+6,j+2,k+1)+2.213902243010268D2*xyzi(i+6,
     $        j,k+1)-3.320853364515403D3*xyzi(i+4,j+2,k+1)-3.76363381311
     $        7456D3*xyzi(i+2,j+6,k+1)+3.320853364515403D3*xyzi(i+2,j+4,
     $        k+1)+2.688309866512469D2*xyzi(i,j+8,k+1)-2.213902243010268
     $        D2*xyzi(i,j+6,k+1))+angi
            angi=yrk(99)*(3.992735113630909D1*xyzi(i+8,j,k+1)-1.11796583
     $        1816654D3*xyzi(i+6,j+2,k+1)+2.794914579541636D3*xyzi(i+4,j
     $        +4,k+1)-1.117965831816654D3*xyzi(i+2,j+6,k+1)+3.9927351136
     $        30909D1*xyzi(i,j+8,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(83)*(3.194188090904727D2*xyzi(i+7,j+1,k+1)-2.235931
     $        663633309D3*xyzi(i+5,j+3,k+1)+2.235931663633309D3*xyzi(i+3
     $        ,j+5,k+1)-3.194188090904727D2*xyzi(i+1,j+7,k+1))+angi
            angi=yrk(85)*(-1.612985919907481D3*xyzi(i+7,j+1,k+1)+3.76363
     $        3813117456D3*xyzi(i+5,j+3,k+1)+1.328341345806161D3*xyzi(i+
     $        5,j+1,k+1)+3.763633813117456D3*xyzi(i+3,j+5,k+1)-4.4278044
     $        86020537D3*xyzi(i+3,j+3,k+1)-1.612985919907481D3*xyzi(i+1,
     $        j+7,k+1)+1.328341345806161D3*xyzi(i+1,j+5,k+1))+angi
            angi=yrk(87)*(3.484447832367752D3*xyzi(i+7,j+1,k+1)+3.484447
     $        832367752D3*xyzi(i+5,j+3,k+1)-4.919220469225061D3*xyzi(i+5
     $        ,j+1,k+1)-3.484447832367752D3*xyzi(i+3,j+5,k+1)+1.63974015
     $        6408354D3*xyzi(i+3,j+1,k+1)-3.484447832367752D3*xyzi(i+1,j
     $        +7,k+1)+4.919220469225061D3*xyzi(i+1,j+5,k+1)-1.6397401564
     $        08354D3*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(89)*(-3.357697947089007D3*xyzi(i+7,j+1,k+1)-1.00730
     $        9384126702D4*xyzi(i+5,j+3,k+1)+5.925349318392365D3*xyzi(i+
     $        5,j+1,k+1)-1.007309384126702D4*xyzi(i+3,j+5,k+1)+1.1850698
     $        63678473D4*xyzi(i+3,j+3,k+1)-3.160186303142594D3*xyzi(i+3,
     $        j+1,k+1)-3.357697947089007D3*xyzi(i+1,j+7,k+1)+5.925349318
     $        392365D3*xyzi(i+1,j+5,k+1)-3.160186303142594D3*xyzi(i+1,j+
     $        3,k+1)+4.861825081757838D2*xyzi(i+1,j+1,k+1))+angi
          endif
          OMEGA(9)=angi
        endif
c...
c...  l=11
c...
        if(lmax.ge.11)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(134)*(7.928933427523546D3*xyzi(i+11,j,k)+3.96446671
     $        3761773D4*xyzi(i+9,j+2,k)-2.265409550721013D4*xyzi(i+9,j,k
     $        )+7.928933427523546D4*xyzi(i+7,j+4,k)-9.061638202884053D4*
     $        xyzi(i+7,j+2,k)+2.384641632337909D4*xyzi(i+7,j,k)+7.928933
     $        427523546D4*xyzi(i+5,j+6,k)-1.359245730432608D5*xyzi(i+5,j
     $        +4,k)+7.153924897013726D4*xyzi(i+5,j+2,k)-1.12218429757078
     $        1D4*xyzi(i+5,j,k)+3.964466713761773D4*xyzi(i+3,j+8,k)-9.06
     $        1638202884053D4*xyzi(i+3,j+6,k)+7.153924897013726D4*xyzi(i
     $        +3,j+4,k)-2.244368595141561D4*xyzi(i+3,j+2,k)+2.2443685951
     $        41561D3*xyzi(i+3,j,k)+7.928933427523546D3*xyzi(i+1,j+10,k)
     $        -2.265409550721013D4*xyzi(i+1,j+8,k)+2.384641632337909D4*x
     $        yzi(i+1,j+6,k)-1.122184297570781D4*xyzi(i+1,j+4,k)+2.24436
     $        8595141561D3*xyzi(i+1,j+2,k)-1.381149904702499D2*xyzi(i+1,
     $        j,k))+angi
            angi=yrk(136)*(-5.575711986679481D3*xyzi(i+11,j,k)-5.5757119
     $        86679481D3*xyzi(i+9,j+2,k)+1.486856529781195D4*xyzi(i+9,j,
     $        k)+3.345427192007688D4*xyzi(i+7,j+4,k)-1.408600922950606D4
     $        *xyzi(i+7,j,k)+7.805996781351273D4*xyzi(i+5,j+6,k)-8.92113
     $        9178687169D4*xyzi(i+5,j+4,k)+1.408600922950606D4*xyzi(i+5,
     $        j+2,k)+5.523925188041591D3*xyzi(i+5,j,k)+6.133283185347429
     $        D4*xyzi(i+3,j+8,k)-1.189485223824956D5*xyzi(i+3,j+6,k)+7.0
     $        43004614753028D4*xyzi(i+3,j+4,k)-1.104785037608318D4*xyzi(
     $        i+3,j+2,k)-7.365233584055455D2*xyzi(i+3,j,k)+1.67271359600
     $        3844D4*xyzi(i+1,j+10,k)-4.460569589343585D4*xyzi(i+1,j+8,k
     $        )+4.225802768851817D4*xyzi(i+1,j+6,k)-1.657177556412477D4*
     $        xyzi(i+1,j+4,k)+2.209570075216636D3*xyzi(i+1,j+2,k))+angi
            angi=yrk(138)*(2.693324767573892D3*xyzi(i+11,j,k)-1.88532733
     $        7301724D4*xyzi(i+9,j+2,k)-6.156170897311752D3*xyzi(i+9,j,k
     $        )-5.925314488662561D4*xyzi(i+7,j+4,k)+4.924936717849402D4*
     $        xyzi(i+7,j+2,k)+4.536125924334975D3*xyzi(i+7,j,k)-3.770654
     $        674603448D4*xyzi(i+5,j+6,k)+8.618639256236453D4*xyzi(i+5,j
     $        +4,k)-4.082513331901478D4*xyzi(i+5,j+2,k)-1.06732374690234
     $        7D3*xyzi(i+5,j,k)+1.346662383786946D4*xyzi(i+3,j+8,k)-2.26
     $        8062962167488D4*xyzi(i+3,j+4,k)+1.067323746902347D4*xyzi(i
     $        +3,j+2,k)+1.346662383786946D4*xyzi(i+1,j+10,k)-3.078085448
     $        655876D4*xyzi(i+1,j+8,k)+2.268062962167488D4*xyzi(i+1,j+6,
     $        k)-5.336618734511736D3*xyzi(i+1,j+4,k))+angi
            angi=yrk(140)*(-8.433126966180176D2*xyzi(i+11,j,k)+1.6022941
     $        23574233D4*xyzi(i+9,j+2,k)+1.44567890848803D3*xyzi(i+9,j,k
     $        )+5.059876179708106D3*xyzi(i+7,j+4,k)-2.89135781697606D4*x
     $        yzi(i+7,j+2,k)-6.087069088370653D2*xyzi(i+7,j,k)-3.5419133
     $        25795674D4*xyzi(i+5,j+6,k)+2.023950471883242D4*xyzi(i+5,j+
     $        4,k)+1.278284508557837D4*xyzi(i+5,j+2,k)-1.770956662897837
     $        D4*xyzi(i+3,j+8,k)+4.047900943766485D4*xyzi(i+3,j+6,k)-2.1
     $        30474180929729D4*xyzi(i+3,j+4,k)+5.903188876326123D3*xyzi(
     $        i+1,j+10,k)-1.011975235941621D4*xyzi(i+1,j+8,k)+4.26094836
     $        1859457D3*xyzi(i+1,j+6,k))+angi
            angi=yrk(142)*(1.49860598832503D2*xyzi(i+11,j,k)-5.245120959
     $        137605D3*xyzi(i+9,j+2,k)-1.42724379840479D2*xyzi(i+9,j,k)+
     $        1.348745389492527D4*xyzi(i+7,j+4,k)+5.138077674257246D3*xy
     $        zi(i+7,j+2,k)+6.294145150965126D3*xyzi(i+5,j+6,k)-1.798327
     $        185990036D4*xyzi(i+5,j+4,k)-1.123954491243772D4*xyzi(i+3,j
     $        +8,k)+1.198884790660024D4*xyzi(i+3,j+6,k)+1.34874538949252
     $        7D3*xyzi(i+1,j+10,k)-1.284519418564311D3*xyzi(i+1,j+8,k))+
     $        angi
            angi=yrk(144)*(-9.860103500953133D0*xyzi(i+11,j,k)+5.4230569
     $        25524223D2*xyzi(i+9,j+2,k)-3.253834155314534D3*xyzi(i+7,j+
     $        4,k)+4.555367817440347D3*xyzi(i+5,j+6,k)-1.626917077657267
     $        D3*xyzi(i+3,j+8,k)+1.084611385104845D2*xyzi(i+1,j+10,k))+a
     $        ngi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(122)*(-1.084611385104845D2*xyzi(i+10,j+1,k)+1.62691
     $        7077657267D3*xyzi(i+8,j+3,k)-4.555367817440347D3*xyzi(i+6,
     $        j+5,k)+3.253834155314534D3*xyzi(i+4,j+7,k)-5.4230569255242
     $        23D2*xyzi(i+2,j+9,k)+9.860103500953133D0*xyzi(i,j+11,k))+a
     $        ngi
            angi=yrk(124)*(1.348745389492527D3*xyzi(i+10,j+1,k)-1.123954
     $        491243772D4*xyzi(i+8,j+3,k)-1.284519418564311D3*xyzi(i+8,j
     $        +1,k)+6.294145150965126D3*xyzi(i+6,j+5,k)+1.19888479066002
     $        4D4*xyzi(i+6,j+3,k)+1.348745389492527D4*xyzi(i+4,j+7,k)-1.
     $        798327185990036D4*xyzi(i+4,j+5,k)-5.245120959137605D3*xyzi
     $        (i+2,j+9,k)+5.138077674257246D3*xyzi(i+2,j+7,k)+1.49860598
     $        832503D2*xyzi(i,j+11,k)-1.42724379840479D2*xyzi(i,j+9,k))+
     $        angi
            angi=yrk(126)*(-5.903188876326123D3*xyzi(i+10,j+1,k)+1.77095
     $        6662897837D4*xyzi(i+8,j+3,k)+1.011975235941621D4*xyzi(i+8,
     $        j+1,k)+3.541913325795674D4*xyzi(i+6,j+5,k)-4.0479009437664
     $        85D4*xyzi(i+6,j+3,k)-4.260948361859457D3*xyzi(i+6,j+1,k)-5
     $        .059876179708106D3*xyzi(i+4,j+7,k)-2.023950471883242D4*xyz
     $        i(i+4,j+5,k)+2.130474180929729D4*xyzi(i+4,j+3,k)-1.6022941
     $        23574233D4*xyzi(i+2,j+9,k)+2.89135781697606D4*xyzi(i+2,j+7
     $        ,k)-1.278284508557837D4*xyzi(i+2,j+5,k)+8.433126966180176D
     $        2*xyzi(i,j+11,k)-1.44567890848803D3*xyzi(i,j+9,k)+6.087069
     $        088370653D2*xyzi(i,j+7,k))+angi
            angi=yrk(128)*(1.346662383786946D4*xyzi(i+10,j+1,k)+1.346662
     $        383786946D4*xyzi(i+8,j+3,k)-3.078085448655876D4*xyzi(i+8,j
     $        +1,k)-3.770654674603448D4*xyzi(i+6,j+5,k)+2.26806296216748
     $        8D4*xyzi(i+6,j+1,k)-5.925314488662561D4*xyzi(i+4,j+7,k)+8.
     $        618639256236453D4*xyzi(i+4,j+5,k)-2.268062962167488D4*xyzi
     $        (i+4,j+3,k)-5.336618734511736D3*xyzi(i+4,j+1,k)-1.88532733
     $        7301724D4*xyzi(i+2,j+9,k)+4.924936717849402D4*xyzi(i+2,j+7
     $        ,k)-4.082513331901478D4*xyzi(i+2,j+5,k)+1.067323746902347D
     $        4*xyzi(i+2,j+3,k)+2.693324767573892D3*xyzi(i,j+11,k)-6.156
     $        170897311752D3*xyzi(i,j+9,k)+4.536125924334975D3*xyzi(i,j+
     $        7,k)-1.067323746902347D3*xyzi(i,j+5,k))+angi
            angi=yrk(130)*(-1.672713596003844D4*xyzi(i+10,j+1,k)-6.13328
     $        3185347429D4*xyzi(i+8,j+3,k)+4.460569589343585D4*xyzi(i+8,
     $        j+1,k)-7.805996781351273D4*xyzi(i+6,j+5,k)+1.1894852238249
     $        56D5*xyzi(i+6,j+3,k)-4.225802768851817D4*xyzi(i+6,j+1,k)-3
     $        .345427192007688D4*xyzi(i+4,j+7,k)+8.921139178687169D4*xyz
     $        i(i+4,j+5,k)-7.043004614753028D4*xyzi(i+4,j+3,k)+1.6571775
     $        56412477D4*xyzi(i+4,j+1,k)+5.575711986679481D3*xyzi(i+2,j+
     $        9,k)-1.408600922950606D4*xyzi(i+2,j+5,k)+1.104785037608318
     $        D4*xyzi(i+2,j+3,k)-2.209570075216636D3*xyzi(i+2,j+1,k)+5.5
     $        75711986679481D3*xyzi(i,j+11,k)-1.486856529781195D4*xyzi(i
     $        ,j+9,k)+1.408600922950606D4*xyzi(i,j+7,k)-5.52392518804159
     $        1D3*xyzi(i,j+5,k)+7.365233584055455D2*xyzi(i,j+3,k))+angi
            angi=yrk(132)*(7.928933427523546D3*xyzi(i+10,j+1,k)+3.964466
     $        713761773D4*xyzi(i+8,j+3,k)-2.265409550721013D4*xyzi(i+8,j
     $        +1,k)+7.928933427523546D4*xyzi(i+6,j+5,k)-9.06163820288405
     $        3D4*xyzi(i+6,j+3,k)+2.384641632337909D4*xyzi(i+6,j+1,k)+7.
     $        928933427523546D4*xyzi(i+4,j+7,k)-1.359245730432608D5*xyzi
     $        (i+4,j+5,k)+7.153924897013726D4*xyzi(i+4,j+3,k)-1.12218429
     $        7570781D4*xyzi(i+4,j+1,k)+3.964466713761773D4*xyzi(i+2,j+9
     $        ,k)-9.061638202884053D4*xyzi(i+2,j+7,k)+7.153924897013726D
     $        4*xyzi(i+2,j+5,k)-2.244368595141561D4*xyzi(i+2,j+3,k)+2.24
     $        4368595141561D3*xyzi(i+2,j+1,k)+7.928933427523546D3*xyzi(i
     $        ,j+11,k)-2.265409550721013D4*xyzi(i,j+9,k)+2.3846416323379
     $        09D4*xyzi(i,j+7,k)-1.122184297570781D4*xyzi(i,j+5,k)+2.244
     $        368595141561D3*xyzi(i,j+3,k)-1.381149904702499D2*xyzi(i,j+
     $        1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(133)*(5.855905424818466D3*xyzi(i,j,k+11)-1.53368951
     $        6023884D4*xyzi(i,j,k+9)+1.452969015180522D4*xyzi(i,j,k+7)-
     $        5.982813591919795D3*xyzi(i,j,k+5)+9.971355986532992D2*xyzi
     $        (i,j,k+3)-4.602164301476766D1*xyzi(i,j,k+1))+angi
            angi=yrk(135)*(6.954134647161096D3*xyzi(i+10,j,k+1)+2.086240
     $        394148329D4*xyzi(i+8,j+2,k+1)-1.589516490779679D4*xyzi(i+8
     $        ,j,k+1)+1.390826929432219D4*xyzi(i+6,j+4,k+1)-3.1790329815
     $        59358D4*xyzi(i+6,j+2,k+1)+1.25488144008922D4*xyzi(i+6,j,k+
     $        1)-1.390826929432219D4*xyzi(i+4,j+6,k+1)+1.25488144008922D
     $        4*xyzi(i+4,j+2,k+1)-3.936882949299515D3*xyzi(i+4,j,k+1)-2.
     $        086240394148329D4*xyzi(i+2,j+8,k+1)+3.179032981559358D4*xy
     $        zi(i+2,j+6,k+1)-1.25488144008922D4*xyzi(i+2,j+4,k+1)+3.936
     $        882949299515D2*xyzi(i+2,j,k+1)-6.954134647161096D3*xyzi(i,
     $        j+10,k+1)+1.589516490779679D4*xyzi(i,j+8,k+1)-1.2548814400
     $        8922D4*xyzi(i,j+6,k+1)+3.936882949299515D3*xyzi(i,j+4,k+1)
     $        -3.936882949299515D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(137)*(-4.071924305675061D3*xyzi(i+10,j,k+1)+1.22157
     $        7291702518D4*xyzi(i+8,j+2,k+1)+8.143848611350123D3*xyzi(i+
     $        8,j,k+1)+5.700694027945086D4*xyzi(i+6,j+4,k+1)-3.257539444
     $        540049D4*xyzi(i+6,j+2,k+1)-5.143483333484288D3*xyzi(i+6,j,
     $        k+1)+5.700694027945086D4*xyzi(i+4,j+6,k+1)-8.1438486113501
     $        23D4*xyzi(i+4,j+4,k+1)+2.571741666742144D4*xyzi(i+4,j+2,k+
     $        1)+1.008526143820449D3*xyzi(i+4,j,k+1)+1.221577291702518D4
     $        *xyzi(i+2,j+8,k+1)-3.257539444540049D4*xyzi(i+2,j+6,k+1)+2
     $        .571741666742144D4*xyzi(i+2,j+4,k+1)-6.051156862922692D3*x
     $        yzi(i+2,j+2,k+1)-4.071924305675061D3*xyzi(i,j+10,k+1)+8.14
     $        3848611350123D3*xyzi(i,j+8,k+1)-5.143483333484288D3*xyzi(i
     $        ,j+6,k+1)+1.008526143820449D3*xyzi(i,j+4,k+1))+angi
            angi=yrk(139)*(1.600073340630907D3*xyzi(i+10,j,k+1)-2.080095
     $        342820179D4*xyzi(i+8,j+2,k+1)-2.438206995247096D3*xyzi(i+8
     $        ,j,k+1)-2.240102676883269D4*xyzi(i+6,j+4,k+1)+3.4134897933
     $        45934D4*xyzi(i+6,j+2,k+1)+8.982867877226143D2*xyzi(i+6,j,k
     $        +1)+2.240102676883269D4*xyzi(i+4,j+6,k+1)-1.34743018158392
     $        1D4*xyzi(i+4,j+2,k+1)+2.080095342820179D4*xyzi(i+2,j+8,k+1
     $        )-3.413489793345934D4*xyzi(i+2,j+6,k+1)+1.347430181583921D
     $        4*xyzi(i+2,j+4,k+1)-1.600073340630907D3*xyzi(i,j+10,k+1)+2
     $        .438206995247096D3*xyzi(i,j+8,k+1)-8.982867877226143D2*xyz
     $        i(i,j+6,k+1))+angi
            angi=yrk(141)*(-3.869384023539698D2*xyzi(i+10,j,k+1)+1.04473
     $        3686355719D4*xyzi(i+8,j+2,k+1)+3.316614877319741D2*xyzi(i+
     $        8,j,k+1)-1.625141289886673D4*xyzi(i+6,j+4,k+1)-9.286521656
     $        495276D3*xyzi(i+6,j+2,k+1)-1.625141289886673D4*xyzi(i+4,j+
     $        6,k+1)+2.321630414123819D4*xyzi(i+4,j+4,k+1)+1.04473368635
     $        5719D4*xyzi(i+2,j+8,k+1)-9.286521656495276D3*xyzi(i+2,j+6,
     $        k+1)-3.869384023539698D2*xyzi(i,j+10,k+1)+3.31661487731974
     $        1D2*xyzi(i,j+8,k+1))+angi
            angi=yrk(143)*(4.624798485436075D1*xyzi(i+10,j,k+1)-2.081159
     $        318446234D3*xyzi(i+8,j+2,k+1)+9.712076819415756D3*xyzi(i+6
     $        ,j+4,k+1)-9.712076819415756D3*xyzi(i+4,j+6,k+1)+2.08115931
     $        8446234D3*xyzi(i+2,j+8,k+1)-4.624798485436075D1*xyzi(i,j+1
     $        0,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(123)*(4.624798485436074D2*xyzi(i+9,j+1,k+1)-5.54975
     $        8182523289D3*xyzi(i+7,j+3,k+1)+1.165449218329891D4*xyzi(i+
     $        5,j+5,k+1)-5.549758182523289D3*xyzi(i+3,j+7,k+1)+4.6247984
     $        85436074D2*xyzi(i+1,j+9,k+1))+angi
            angi=yrk(125)*(-3.095507218831759D3*xyzi(i+9,j+1,k+1)+1.8573
     $        04331299055D4*xyzi(i+7,j+3,k+1)+2.653291901855793D3*xyzi(i
     $        +7,j+1,k+1)-1.857304331299055D4*xyzi(i+5,j+3,k+1)-1.857304
     $        331299055D4*xyzi(i+3,j+7,k+1)+1.857304331299055D4*xyzi(i+3
     $        ,j+5,k+1)+3.095507218831759D3*xyzi(i+1,j+9,k+1)-2.65329190
     $        1855793D3*xyzi(i+1,j+7,k+1))+angi
            angi=yrk(127)*(9.60044004378544D3*xyzi(i+9,j+1,k+1)-1.280058
     $        672504725D4*xyzi(i+7,j+3,k+1)-1.462924197148258D4*xyzi(i+7
     $        ,j+1,k+1)-4.480205353766539D4*xyzi(i+5,j+5,k+1)+3.41348979
     $        3345934D4*xyzi(i+5,j+3,k+1)+5.389720726335685D3*xyzi(i+5,j
     $        +1,k+1)-1.280058672504725D4*xyzi(i+3,j+7,k+1)+3.4134897933
     $        45934D4*xyzi(i+3,j+5,k+1)-1.796573575445229D4*xyzi(i+3,j+3
     $        ,k+1)+9.60044004378544D3*xyzi(i+1,j+9,k+1)-1.4629241971482
     $        58D4*xyzi(i+1,j+7,k+1)+5.389720726335685D3*xyzi(i+1,j+5,k+
     $        1))+angi
            angi=yrk(129)*(-1.628769722270025D4*xyzi(i+9,j+1,k+1)-3.2575
     $        39444540049D4*xyzi(i+7,j+3,k+1)+3.257539444540049D4*xyzi(i
     $        +7,j+1,k+1)+3.257539444540049D4*xyzi(i+5,j+3,k+1)-2.057393
     $        333393715D4*xyzi(i+5,j+1,k+1)+3.257539444540049D4*xyzi(i+3
     $        ,j+7,k+1)-3.257539444540049D4*xyzi(i+3,j+5,k+1)+4.03410457
     $        5281795D3*xyzi(i+3,j+1,k+1)+1.628769722270025D4*xyzi(i+1,j
     $        +9,k+1)-3.257539444540049D4*xyzi(i+1,j+7,k+1)+2.0573933333
     $        93715D4*xyzi(i+1,j+5,k+1)-4.034104575281795D3*xyzi(i+1,j+3
     $        ,k+1))+angi
            angi=yrk(131)*(1.390826929432219D4*xyzi(i+9,j+1,k+1)+5.56330
     $        7717728877D4*xyzi(i+7,j+3,k+1)-3.179032981559358D4*xyzi(i+
     $        7,j+1,k+1)+8.344961576593315D4*xyzi(i+5,j+5,k+1)-9.5370989
     $        44678074D4*xyzi(i+5,j+3,k+1)+2.509762880178441D4*xyzi(i+5,
     $        j+1,k+1)+5.563307717728877D4*xyzi(i+3,j+7,k+1)-9.537098944
     $        678074D4*xyzi(i+3,j+5,k+1)+5.019525760356881D4*xyzi(i+3,j+
     $        3,k+1)-7.873765898599029D3*xyzi(i+3,j+1,k+1)+1.39082692943
     $        2219D4*xyzi(i+1,j+9,k+1)-3.179032981559358D4*xyzi(i+1,j+7,
     $        k+1)+2.509762880178441D4*xyzi(i+1,j+5,k+1)-7.8737658985990
     $        29D3*xyzi(i+1,j+3,k+1)+7.873765898599029D2*xyzi(i+1,j+1,k+
     $        1))+angi
          endif
          OMEGA(11)=angi
        endif
c...
c...  l=13
c...
        if(lmax.ge.13)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(184)*(-3.186969598569298D4*xyzi(i+13,j,k)-1.9121817
     $        59141579D5*xyzi(i+11,j+2,k)+1.070821785119284D5*xyzi(i+11,
     $        j,k)-4.780454397853947D5*xyzi(i+9,j+4,k)+5.35410892559642D
     $        5*xyzi(i+9,j+2,k)-1.396724067546892D5*xyzi(i+9,j,k)-6.3739
     $        39197138596D5*xyzi(i+7,j+6,k)+1.070821785119284D6*xyzi(i+7
     $        ,j+4,k)-5.586896270187569D5*xyzi(i+7,j+2,k)+8.868089317758
     $        046D4*xyzi(i+7,j,k)-4.780454397853947D5*xyzi(i+5,j+8,k)+1.
     $        070821785119284D6*xyzi(i+5,j+6,k)-8.380344405281354D5*xyzi
     $        (i+5,j+4,k)+2.660426795327414D5*xyzi(i+5,j+2,k)-2.80044925
     $        8239383D4*xyzi(i+5,j,k)-1.912181759141579D5*xyzi(i+3,j+10,
     $        k)+5.35410892559642D5*xyzi(i+3,j+8,k)-5.586896270187569D5*
     $        xyzi(i+3,j+6,k)+2.660426795327414D5*xyzi(i+3,j+4,k)-5.6008
     $        98516478766D4*xyzi(i+3,j+2,k)+3.953575423396776D3*xyzi(i+3
     $        ,j,k)-3.186969598569298D4*xyzi(i+1,j+12,k)+1.0708217851192
     $        84D5*xyzi(i+1,j+10,k)-1.396724067546892D5*xyzi(i+1,j+8,k)+
     $        8.868089317758046D4*xyzi(i+1,j+6,k)-2.800449258239383D4*xy
     $        zi(i+1,j+4,k)+3.953575423396776D3*xyzi(i+1,j+2,k)-1.757144
     $        632620789D2*xyzi(i+1,j,k))+angi
            angi=yrk(186)*(2.36351991153295D4*xyzi(i+13,j,k)+4.727039823
     $        0659D4*xyzi(i+11,j+2,k)-7.56326371690544D4*xyzi(i+11,j,k)-
     $        1.181759955766475D5*xyzi(i+9,j+4,k)-7.56326371690544D4*xyz
     $        i(i+9,j+2,k)+9.2074514814501D4*xyzi(i+9,j,k)-4.72703982306
     $        59D5*xyzi(i+7,j+6,k)+4.537958230143264D5*xyzi(i+7,j+4,k)-5
     $        .261400846542915D4*xyzi(i+7,j,k)-5.908799778832375D5*xyzi(
     $        i+5,j+8,k)+1.058856920366762D6*xyzi(i+5,j+6,k)-5.524470888
     $        87006D5*xyzi(i+5,j+4,k)+5.261400846542915D4*xyzi(i+5,j+2,k
     $        )+1.384579170142872D4*xyzi(i+5,j,k)-3.30892787614613D5*xyz
     $        i(i+3,j+10,k)+8.319590088595984D5*xyzi(i+3,j+8,k)-7.365961
     $        18516008D5*xyzi(i+3,j+6,k)+2.630700423271457D5*xyzi(i+3,j+
     $        4,k)-2.769158340285745D4*xyzi(i+3,j+2,k)-1.303133336605056
     $        D3*xyzi(i+3,j,k)-7.09055973459885D4*xyzi(i+1,j+12,k)+2.268
     $        979115071632D5*xyzi(i+1,j+10,k)-2.76223544443503D5*xyzi(i+
     $        1,j+8,k)+1.578420253962874D5*xyzi(i+1,j+6,k)-4.15373751042
     $        8617D4*xyzi(i+1,j+4,k)+3.909400009815169D3*xyzi(i+1,j+2,k)
     $        )+angi
            angi=yrk(188)*(-1.281798641180881D4*xyzi(i+13,j,k)+7.6907918
     $        47085288D4*xyzi(i+11,j+2,k)+3.691580086600938D4*xyzi(i+11,
     $        j,k)+3.717216059424556D5*xyzi(i+9,j+4,k)-2.584106060620657
     $        D5*xyzi(i+9,j+2,k)-3.852083568627066D4*xyzi(i+9,j,k)+4.614
     $        475108251173D5*xyzi(i+7,j+6,k)-8.121476190522064D5*xyzi(i+
     $        7,j+4,k)+3.081666854901653D5*xyzi(i+7,j+2,k)+1.71203714161
     $        2029D4*xyzi(i+7,j,k)+1.153618777062793D5*xyzi(i+5,j+8,k)-5
     $        .168212121241314D5*xyzi(i+5,j+6,k)+5.392916996077893D5*xyz
     $        i(i+5,j+4,k)-1.540833427450826D5*xyzi(i+5,j+2,k)-2.7032165
     $        39387415D3*xyzi(i+5,j,k)-1.281798641180881D5*xyzi(i+3,j+10
     $        ,k)+1.845790043300469D5*xyzi(i+3,j+8,k)-8.560185708060147D
     $        4*xyzi(i+3,j+4,k)+2.703216539387415D4*xyzi(i+3,j+2,k)-6.40
     $        8993205904407D4*xyzi(i+1,j+12,k)+1.845790043300469D5*xyzi(
     $        i+1,j+10,k)-1.926041784313533D5*xyzi(i+1,j+8,k)+8.56018570
     $        8060147D4*xyzi(i+1,j+6,k)-1.351608269693707D4*xyzi(i+1,j+4
     $        ,k))+angi
            angi=yrk(190)*(4.920644864827347D3*xyzi(i+13,j,k)-8.85716075
     $        6689225D4*xyzi(i+11,j+2,k)-1.180954767558563D4*xyzi(i+11,j
     $        ,k)-1.230161216206837D5*xyzi(i+9,j+4,k)+2.24381405836127D5
     $        *xyzi(i+9,j+2,k)+9.242254702632235D3*xyzi(i+9,j,k)+1.77143
     $        2151337845D5*xyzi(i+7,j+6,k)+7.08572860535138D4*xyzi(i+7,j
     $        +4,k)-1.848450940526447D5*xyzi(i+7,j+2,k)-2.34723928955739
     $        3D3*xyzi(i+7,j,k)+3.100006264841229D5*xyzi(i+5,j+8,k)-4.96
     $        0010023745966D5*xyzi(i+5,j+6,k)+1.293915658368513D5*xyzi(i
     $        +5,j+4,k)+4.929202508070525D4*xyzi(i+5,j+2,k)+6.8889028107
     $        58286D4*xyzi(i+3,j+10,k)-2.480005011872983D5*xyzi(i+3,j+8,
     $        k)+2.587831316737026D5*xyzi(i+3,j+6,k)-8.215337513450875D4
     $        *xyzi(i+3,j+4,k)-3.444451405379143D4*xyzi(i+1,j+12,k)+8.26
     $        6683372909943D4*xyzi(i+1,j+10,k)-6.469578291842564D4*xyzi(
     $        i+1,j+8,k)+1.643067502690175D4*xyzi(i+1,j+6,k))+angi
            angi=yrk(192)*(-1.253896417710616D3*xyzi(i+13,j,k)+4.2632478
     $        20216095D4*xyzi(i+11,j+2,k)+2.206857695170684D3*xyzi(i+11,
     $        j,k)-6.896430297408388D4*xyzi(i+9,j+4,k)-7.724001933097395
     $        D4*xyzi(i+9,j+2,k)-9.595033457263844D2*xyzi(i+9,j,k)-1.655
     $        143271378013D5*xyzi(i+7,j+6,k)+1.986171925653616D5*xyzi(i+
     $        7,j+4,k)+3.454212044614984D4*xyzi(i+7,j+2,k)+4.13785817844
     $        5033D4*xyzi(i+5,j+8,k)+9.268802319716874D4*xyzi(i+5,j+6,k)
     $        -1.208974215615244D5*xyzi(i+5,j+4,k)+8.275716356890066D4*x
     $        yzi(i+3,j+10,k)-1.655143271378013D5*xyzi(i+3,j+8,k)+8.0598
     $        28104101629D4*xyzi(i+3,j+6,k)-1.128506775939554D4*xyzi(i+1
     $        ,j+12,k)+1.986171925653616D4*xyzi(i+1,j+10,k)-8.6355301115
     $        3746D3*xyzi(i+1,j+8,k))+angi
            angi=yrk(194)*(1.848769406428712D2*xyzi(i+13,j,k)-9.98335479
     $        4715047D3*xyzi(i+11,j+2,k)-1.774818630171564D2*xyzi(i+11,j
     $        ,k)+5.084115867678959D4*xyzi(i+9,j+4,k)+9.761502465943601D
     $        3*xyzi(i+9,j+2,k)-2.4403756164859D4*xyzi(i+7,j+6,k)-5.8569
     $        01479566161D4*xyzi(i+7,j+4,k)-5.490845137093276D4*xyzi(i+5
     $        ,j+8,k)+8.199662071392625D4*xyzi(i+5,j+6,k)+2.847104885900
     $        217D4*xyzi(i+3,j+10,k)-2.92845073978308D4*xyzi(i+3,j+8,k)-
     $        2.033646347071584D3*xyzi(i+1,j+12,k)+1.95230049318872D3*xy
     $        zi(i+1,j+10,k))+angi
            angi=yrk(196)*(-1.025512752521207D1*xyzi(i+13,j,k)+7.9989994
     $        69665415D2*xyzi(i+11,j+2,k)-7.332416180526631D3*xyzi(i+9,j
     $        +4,k)+1.759779883326391D4*xyzi(i+7,j+6,k)-1.31983491249479
     $        4D4*xyzi(i+5,j+8,k)+2.932966472210652D3*xyzi(i+3,j+10,k)-1
     $        .333166578277569D2*xyzi(i+1,j+12,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(170)*(-1.333166578277569D2*xyzi(i+12,j+1,k)+2.93296
     $        6472210652D3*xyzi(i+10,j+3,k)-1.319834912494794D4*xyzi(i+8
     $        ,j+5,k)+1.759779883326391D4*xyzi(i+6,j+7,k)-7.332416180526
     $        631D3*xyzi(i+4,j+9,k)+7.998999469665415D2*xyzi(i+2,j+11,k)
     $        -1.025512752521207D1*xyzi(i,j+13,k))+angi
            angi=yrk(172)*(2.033646347071584D3*xyzi(i+12,j+1,k)-2.847104
     $        885900217D4*xyzi(i+10,j+3,k)-1.95230049318872D3*xyzi(i+10,
     $        j+1,k)+5.490845137093276D4*xyzi(i+8,j+5,k)+2.9284507397830
     $        8D4*xyzi(i+8,j+3,k)+2.4403756164859D4*xyzi(i+6,j+7,k)-8.19
     $        9662071392625D4*xyzi(i+6,j+5,k)-5.084115867678959D4*xyzi(i
     $        +4,j+9,k)+5.856901479566161D4*xyzi(i+4,j+7,k)+9.9833547947
     $        15047D3*xyzi(i+2,j+11,k)-9.761502465943601D3*xyzi(i+2,j+9,
     $        k)-1.848769406428712D2*xyzi(i,j+13,k)+1.774818630171564D2*
     $        xyzi(i,j+11,k))+angi
            angi=yrk(174)*(-1.128506775939554D4*xyzi(i+12,j+1,k)+8.27571
     $        6356890066D4*xyzi(i+10,j+3,k)+1.986171925653616D4*xyzi(i+1
     $        0,j+1,k)+4.137858178445033D4*xyzi(i+8,j+5,k)-1.65514327137
     $        8013D5*xyzi(i+8,j+3,k)-8.63553011153746D3*xyzi(i+8,j+1,k)-
     $        1.655143271378013D5*xyzi(i+6,j+7,k)+9.268802319716874D4*xy
     $        zi(i+6,j+5,k)+8.059828104101629D4*xyzi(i+6,j+3,k)-6.896430
     $        297408388D4*xyzi(i+4,j+9,k)+1.986171925653616D5*xyzi(i+4,j
     $        +7,k)-1.208974215615244D5*xyzi(i+4,j+5,k)+4.26324782021609
     $        5D4*xyzi(i+2,j+11,k)-7.724001933097395D4*xyzi(i+2,j+9,k)+3
     $        .454212044614984D4*xyzi(i+2,j+7,k)-1.253896417710616D3*xyz
     $        i(i,j+13,k)+2.206857695170684D3*xyzi(i,j+11,k)-9.595033457
     $        263844D2*xyzi(i,j+9,k))+angi
            angi=yrk(176)*(3.444451405379143D4*xyzi(i+12,j+1,k)-6.888902
     $        810758286D4*xyzi(i+10,j+3,k)-8.266683372909943D4*xyzi(i+10
     $        ,j+1,k)-3.100006264841229D5*xyzi(i+8,j+5,k)+2.480005011872
     $        983D5*xyzi(i+8,j+3,k)+6.469578291842564D4*xyzi(i+8,j+1,k)-
     $        1.771432151337845D5*xyzi(i+6,j+7,k)+4.960010023745966D5*xy
     $        zi(i+6,j+5,k)-2.587831316737026D5*xyzi(i+6,j+3,k)-1.643067
     $        502690175D4*xyzi(i+6,j+1,k)+1.230161216206837D5*xyzi(i+4,j
     $        +9,k)-7.08572860535138D4*xyzi(i+4,j+7,k)-1.293915658368513
     $        D5*xyzi(i+4,j+5,k)+8.215337513450875D4*xyzi(i+4,j+3,k)+8.8
     $        57160756689225D4*xyzi(i+2,j+11,k)-2.24381405836127D5*xyzi(
     $        i+2,j+9,k)+1.848450940526447D5*xyzi(i+2,j+7,k)-4.929202508
     $        070525D4*xyzi(i+2,j+5,k)-4.920644864827347D3*xyzi(i,j+13,k
     $        )+1.180954767558563D4*xyzi(i,j+11,k)-9.242254702632235D3*x
     $        yzi(i,j+9,k)+2.347239289557393D3*xyzi(i,j+7,k))+angi
            angi=yrk(178)*(-6.408993205904407D4*xyzi(i+12,j+1,k)-1.28179
     $        8641180881D5*xyzi(i+10,j+3,k)+1.845790043300469D5*xyzi(i+1
     $        0,j+1,k)+1.153618777062793D5*xyzi(i+8,j+5,k)+1.84579004330
     $        0469D5*xyzi(i+8,j+3,k)-1.926041784313533D5*xyzi(i+8,j+1,k)
     $        +4.614475108251173D5*xyzi(i+6,j+7,k)-5.168212121241314D5*x
     $        yzi(i+6,j+5,k)+8.560185708060147D4*xyzi(i+6,j+1,k)+3.71721
     $        6059424556D5*xyzi(i+4,j+9,k)-8.121476190522064D5*xyzi(i+4,
     $        j+7,k)+5.392916996077893D5*xyzi(i+4,j+5,k)-8.5601857080601
     $        47D4*xyzi(i+4,j+3,k)-1.351608269693707D4*xyzi(i+4,j+1,k)+7
     $        .690791847085288D4*xyzi(i+2,j+11,k)-2.584106060620657D5*xy
     $        zi(i+2,j+9,k)+3.081666854901653D5*xyzi(i+2,j+7,k)-1.540833
     $        427450826D5*xyzi(i+2,j+5,k)+2.703216539387415D4*xyzi(i+2,j
     $        +3,k)-1.281798641180881D4*xyzi(i,j+13,k)+3.691580086600938
     $        D4*xyzi(i,j+11,k)-3.852083568627066D4*xyzi(i,j+9,k)+1.7120
     $        37141612029D4*xyzi(i,j+7,k)-2.703216539387415D3*xyzi(i,j+5
     $        ,k))+angi
            angi=yrk(180)*(7.09055973459885D4*xyzi(i+12,j+1,k)+3.3089278
     $        7614613D5*xyzi(i+10,j+3,k)-2.268979115071632D5*xyzi(i+10,j
     $        +1,k)+5.908799778832375D5*xyzi(i+8,j+5,k)-8.31959008859598
     $        4D5*xyzi(i+8,j+3,k)+2.76223544443503D5*xyzi(i+8,j+1,k)+4.7
     $        270398230659D5*xyzi(i+6,j+7,k)-1.058856920366762D6*xyzi(i+
     $        6,j+5,k)+7.36596118516008D5*xyzi(i+6,j+3,k)-1.578420253962
     $        874D5*xyzi(i+6,j+1,k)+1.181759955766475D5*xyzi(i+4,j+9,k)-
     $        4.537958230143264D5*xyzi(i+4,j+7,k)+5.52447088887006D5*xyz
     $        i(i+4,j+5,k)-2.630700423271457D5*xyzi(i+4,j+3,k)+4.1537375
     $        10428617D4*xyzi(i+4,j+1,k)-4.7270398230659D4*xyzi(i+2,j+11
     $        ,k)+7.56326371690544D4*xyzi(i+2,j+9,k)-5.261400846542915D4
     $        *xyzi(i+2,j+5,k)+2.769158340285745D4*xyzi(i+2,j+3,k)-3.909
     $        400009815169D3*xyzi(i+2,j+1,k)-2.36351991153295D4*xyzi(i,j
     $        +13,k)+7.56326371690544D4*xyzi(i,j+11,k)-9.2074514814501D4
     $        *xyzi(i,j+9,k)+5.261400846542915D4*xyzi(i,j+7,k)-1.3845791
     $        70142872D4*xyzi(i,j+5,k)+1.303133336605056D3*xyzi(i,j+3,k)
     $        )+angi
            angi=yrk(182)*(-3.186969598569298D4*xyzi(i+12,j+1,k)-1.91218
     $        1759141579D5*xyzi(i+10,j+3,k)+1.070821785119284D5*xyzi(i+1
     $        0,j+1,k)-4.780454397853947D5*xyzi(i+8,j+5,k)+5.35410892559
     $        642D5*xyzi(i+8,j+3,k)-1.396724067546892D5*xyzi(i+8,j+1,k)-
     $        6.373939197138596D5*xyzi(i+6,j+7,k)+1.070821785119284D6*xy
     $        zi(i+6,j+5,k)-5.586896270187569D5*xyzi(i+6,j+3,k)+8.868089
     $        317758046D4*xyzi(i+6,j+1,k)-4.780454397853947D5*xyzi(i+4,j
     $        +9,k)+1.070821785119284D6*xyzi(i+4,j+7,k)-8.38034440528135
     $        4D5*xyzi(i+4,j+5,k)+2.660426795327414D5*xyzi(i+4,j+3,k)-2.
     $        800449258239383D4*xyzi(i+4,j+1,k)-1.912181759141579D5*xyzi
     $        (i+2,j+11,k)+5.35410892559642D5*xyzi(i+2,j+9,k)-5.58689627
     $        0187569D5*xyzi(i+2,j+7,k)+2.660426795327414D5*xyzi(i+2,j+5
     $        ,k)-5.600898516478766D4*xyzi(i+2,j+3,k)+3.953575423396776D
     $        3*xyzi(i+2,j+1,k)-3.186969598569298D4*xyzi(i,j+13,k)+1.070
     $        821785119284D5*xyzi(i,j+11,k)-1.396724067546892D5*xyzi(i,j
     $        +9,k)+8.868089317758046D4*xyzi(i,j+7,k)-2.800449258239383D
     $        4*xyzi(i,j+5,k)+3.953575423396776D3*xyzi(i,j+3,k)-1.757144
     $        632620789D2*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(183)*(2.338596333691754D4*xyzi(i,j,k+13)-7.29642056
     $        1118272D4*xyzi(i,j,k+11)+8.72398110568489D4*xyzi(i,j,k+9)-
     $        4.985132060391366D4*xyzi(i,j,k+7)+1.377470700897614D4*xyzi
     $        (i,j,k+5)-1.620553765761899D3*xyzi(i,j,k+3)+5.401845885872
     $        997D1*xyzi(i,j,k+1))+angi
            angi=yrk(185)*(-2.850512265850467D4*xyzi(i+12,j,k+1)-1.14020
     $        4906340187D5*xyzi(i+10,j+2,k+1)+7.981434344381306D4*xyzi(i
     $        +10,j,k+1)-1.425256132925233D5*xyzi(i+8,j+4,k+1)+2.3944303
     $        03314392D5*xyzi(i+8,j+2,k+1)-8.328453228919624D4*xyzi(i+8,
     $        j,k+1)+1.596286868876261D5*xyzi(i+6,j+4,k+1)-1.66569064578
     $        3925D5*xyzi(i+6,j+2,k+1)+3.965930109009345D4*xyzi(i+6,j,k+
     $        1)+1.425256132925233D5*xyzi(i+4,j+8,k+1)-1.596286868876261
     $        D5*xyzi(i+4,j+6,k+1)+3.965930109009345D4*xyzi(i+4,j+2,k+1)
     $        -8.349326545282831D3*xyzi(i+4,j,k+1)+1.140204906340187D5*x
     $        yzi(i+2,j+10,k+1)-2.394430303314392D5*xyzi(i+2,j+8,k+1)+1.
     $        665690645783925D5*xyzi(i+2,j+6,k+1)-3.965930109009345D4*xy
     $        zi(i+2,j+4,k+1)+5.893642267258469D2*xyzi(i+2,j,k+1)+2.8505
     $        12265850467D4*xyzi(i,j+12,k+1)-7.981434344381306D4*xyzi(i,
     $        j+10,k+1)+8.328453228919624D4*xyzi(i,j+8,k+1)-3.9659301090
     $        09345D4*xyzi(i,j+6,k+1)+8.349326545282831D3*xyzi(i,j+4,k+1
     $        )-5.893642267258469D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(187)*(1.812737022589407D4*xyzi(i+12,j,k+1)-3.625474
     $        045178814D4*xyzi(i+10,j+2,k+1)-4.640606777828882D4*xyzi(i+
     $        10,j,k+1)-3.081652938401992D5*xyzi(i+8,j+4,k+1)+1.39218203
     $        3348665D5*xyzi(i+8,j+2,k+1)+4.237075753669848D4*xyzi(i+8,j
     $        ,k+1)-5.075663663250339D5*xyzi(i+6,j+6,k+1)+6.496849488960
     $        434D5*xyzi(i+6,j+4,k+1)-1.694830301467939D5*xyzi(i+6,j+2,k
     $        +1)-1.614124096636133D4*xyzi(i+6,j,k+1)-3.081652938401992D
     $        5*xyzi(i+4,j+8,k+1)+6.496849488960434D5*xyzi(i+4,j+6,k+1)-
     $        4.237075753669848D5*xyzi(i+4,j+4,k+1)+8.070620483180664D4*
     $        xyzi(i+4,j+2,k+1)+2.123847495573859D3*xyzi(i+4,j,k+1)-3.62
     $        5474045178814D4*xyzi(i+2,j+10,k+1)+1.392182033348665D5*xyz
     $        i(i+2,j+8,k+1)-1.694830301467939D5*xyzi(i+2,j+6,k+1)+8.070
     $        620483180664D4*xyzi(i+2,j+4,k+1)-1.274308497344315D4*xyzi(
     $        i+2,j+2,k+1)+1.812737022589407D4*xyzi(i,j+12,k+1)-4.640606
     $        777828882D4*xyzi(i,j+10,k+1)+4.237075753669848D4*xyzi(i,j+
     $        8,k+1)-1.614124096636133D4*xyzi(i,j+6,k+1)+2.1238474955738
     $        59D3*xyzi(i,j+4,k+1))+angi
            angi=yrk(189)*(-8.317407887033718D3*xyzi(i+12,j,k+1)+9.98088
     $        9464440461D4*xyzi(i+10,j+2,k+1)+1.796560103599283D4*xyzi(i
     $        +10,j,k+1)+2.245700129499104D5*xyzi(i+8,j+4,k+1)-2.3355281
     $        34679068D5*xyzi(i+8,j+2,k+1)-1.249780941634284D4*xyzi(i+8,
     $        j,k+1)-2.515184145038996D5*xyzi(i+6,j+4,k+1)+1.74969331828
     $        7997D5*xyzi(i+6,j+2,k+1)+2.777290981409519D3*xyzi(i+6,j,k+
     $        1)-2.245700129499104D5*xyzi(i+4,j+8,k+1)+2.515184145038996
     $        D5*xyzi(i+4,j+6,k+1)-4.165936472114279D4*xyzi(i+4,j+2,k+1)
     $        -9.980889464440461D4*xyzi(i+2,j+10,k+1)+2.335528134679068D
     $        5*xyzi(i+2,j+8,k+1)-1.749693318287997D5*xyzi(i+2,j+6,k+1)+
     $        4.165936472114279D4*xyzi(i+2,j+4,k+1)+8.317407887033718D3*
     $        xyzi(i,j+12,k+1)-1.796560103599283D4*xyzi(i,j+10,k+1)+1.24
     $        9780941634284D4*xyzi(i,j+8,k+1)-2.777290981409519D3*xyzi(i
     $        ,j+6,k+1))+angi
            angi=yrk(191)*(2.630195315167501D3*xyzi(i+12,j,k+1)-6.838507
     $        819435502D4*xyzi(i+10,j+2,k+1)-4.208312504268001D3*xyzi(i+
     $        10,j,k+1)+3.945292972751251D4*xyzi(i+8,j+4,k+1)+1.13624437
     $        615236D5*xyzi(i+8,j+2,k+1)+1.646730979930957D3*xyzi(i+8,j,
     $        k+1)+2.209364064740701D5*xyzi(i+6,j+6,k+1)-1.7674912517925
     $        61D5*xyzi(i+6,j+4,k+1)-4.61084674380668D4*xyzi(i+6,j+2,k+1
     $        )+3.945292972751251D4*xyzi(i+4,j+8,k+1)-1.767491251792561D
     $        5*xyzi(i+4,j+6,k+1)+1.15271168595167D5*xyzi(i+4,j+4,k+1)-6
     $        .838507819435502D4*xyzi(i+2,j+10,k+1)+1.13624437615236D5*x
     $        yzi(i+2,j+8,k+1)-4.61084674380668D4*xyzi(i+2,j+6,k+1)+2.63
     $        0195315167501D3*xyzi(i,j+12,k+1)-4.208312504268001D3*xyzi(
     $        i,j+10,k+1)+1.646730979930957D3*xyzi(i,j+8,k+1))+angi
            angi=yrk(193)*(-5.229109536543883D2*xyzi(i+12,j,k+1)+2.30080
     $        8196079309D4*xyzi(i+10,j+2,k+1)+4.601616392158617D2*xyzi(i
     $        +10,j,k+1)-8.628030735297407D4*xyzi(i+8,j+4,k+1)-2.0707273
     $        76471378D4*xyzi(i+8,j+2,k+1)+9.663394423533096D4*xyzi(i+6,
     $        j+4,k+1)+8.628030735297407D4*xyzi(i+4,j+8,k+1)-9.663394423
     $        533096D4*xyzi(i+4,j+6,k+1)-2.300808196079309D4*xyzi(i+2,j+
     $        10,k+1)+2.070727376471378D4*xyzi(i+2,j+8,k+1)+5.2291095365
     $        43883D2*xyzi(i,j+12,k+1)-4.601616392158617D2*xyzi(i,j+10,k
     $        +1))+angi
            angi=yrk(195)*(5.229109536543883D1*xyzi(i+12,j,k+1)-3.451212
     $        294118963D3*xyzi(i+10,j+2,k+1)+2.588409220589222D4*xyzi(i+
     $        8,j+4,k+1)-4.831697211766548D4*xyzi(i+6,j+6,k+1)+2.5884092
     $        20589222D4*xyzi(i+4,j+8,k+1)-3.451212294118963D3*xyzi(i+2,
     $        j+10,k+1)+5.229109536543883D1*xyzi(i,j+12,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(171)*(6.27493144385266D2*xyzi(i+11,j+1,k+1)-1.15040
     $        4098039654D4*xyzi(i+9,j+3,k+1)+4.141454752942756D4*xyzi(i+
     $        7,j+5,k+1)-4.141454752942756D4*xyzi(i+5,j+7,k+1)+1.1504040
     $        98039654D4*xyzi(i+3,j+9,k+1)-6.27493144385266D2*xyzi(i+1,j
     $        +11,k+1))+angi
            angi=yrk(173)*(-5.229109536543883D3*xyzi(i+11,j+1,k+1)+5.752
     $        020490198272D4*xyzi(i+9,j+3,k+1)+4.601616392158617D3*xyzi(
     $        i+9,j+1,k+1)-6.902424588237926D4*xyzi(i+7,j+5,k+1)-5.52193
     $        9670590341D4*xyzi(i+7,j+3,k+1)-6.902424588237926D4*xyzi(i+
     $        5,j+7,k+1)+1.159607330823972D5*xyzi(i+5,j+5,k+1)+5.7520204
     $        90198272D4*xyzi(i+3,j+9,k+1)-5.521939670590341D4*xyzi(i+3,
     $        j+7,k+1)-5.229109536543883D3*xyzi(i+1,j+11,k+1)+4.60161639
     $        2158617D3*xyzi(i+1,j+9,k+1))+angi
            angi=yrk(175)*(2.104156252134001D4*xyzi(i+11,j+1,k+1)-1.0520
     $        78126067D5*xyzi(i+9,j+3,k+1)-3.366650003414401D4*xyzi(i+9,
     $        j+1,k+1)-1.2624937512804D5*xyzi(i+7,j+5,k+1)+2.01999000204
     $        8641D5*xyzi(i+7,j+3,k+1)+1.317384783944766D4*xyzi(i+7,j+1,
     $        k+1)+1.2624937512804D5*xyzi(i+5,j+7,k+1)-9.221693487613359
     $        D4*xyzi(i+5,j+3,k+1)+1.052078126067D5*xyzi(i+3,j+9,k+1)-2.
     $        019990002048641D5*xyzi(i+3,j+7,k+1)+9.221693487613359D4*xy
     $        zi(i+3,j+5,k+1)-2.104156252134001D4*xyzi(i+1,j+11,k+1)+3.3
     $        66650003414401D4*xyzi(i+1,j+9,k+1)-1.317384783944766D4*xyz
     $        i(i+1,j+7,k+1))+angi
            angi=yrk(177)*(-4.990444732220231D4*xyzi(i+11,j+1,k+1)+1.663
     $        481577406744D4*xyzi(i+9,j+3,k+1)+1.07793606215957D5*xyzi(i
     $        +9,j+1,k+1)+2.994266839332138D5*xyzi(i+7,j+5,k+1)-1.437248
     $        082879426D5*xyzi(i+7,j+3,k+1)-7.498685649805703D4*xyzi(i+7
     $        ,j+1,k+1)+2.994266839332138D5*xyzi(i+5,j+7,k+1)-5.03036829
     $        0077993D5*xyzi(i+5,j+5,k+1)+1.749693318287997D5*xyzi(i+5,j
     $        +3,k+1)+1.666374588845712D4*xyzi(i+5,j+1,k+1)+1.6634815774
     $        06744D4*xyzi(i+3,j+9,k+1)-1.437248082879426D5*xyzi(i+3,j+7
     $        ,k+1)+1.749693318287997D5*xyzi(i+3,j+5,k+1)-5.554581962819
     $        039D4*xyzi(i+3,j+3,k+1)-4.990444732220231D4*xyzi(i+1,j+11,
     $        k+1)+1.07793606215957D5*xyzi(i+1,j+9,k+1)-7.49868564980570
     $        3D4*xyzi(i+1,j+7,k+1)+1.666374588845712D4*xyzi(i+1,j+5,k+1
     $        ))+angi
            angi=yrk(179)*(7.250948090357628D4*xyzi(i+11,j+1,k+1)+2.1752
     $        84427107288D5*xyzi(i+9,j+3,k+1)-1.856242711131553D5*xyzi(i
     $        +9,j+1,k+1)+1.450189618071526D5*xyzi(i+7,j+5,k+1)-3.712485
     $        422263105D5*xyzi(i+7,j+3,k+1)+1.694830301467939D5*xyzi(i+7
     $        ,j+1,k+1)-1.450189618071526D5*xyzi(i+5,j+7,k+1)+1.69483030
     $        1467939D5*xyzi(i+5,j+3,k+1)-6.456496386544531D4*xyzi(i+5,j
     $        +1,k+1)-2.175284427107288D5*xyzi(i+3,j+9,k+1)+3.7124854222
     $        63105D5*xyzi(i+3,j+7,k+1)-1.694830301467939D5*xyzi(i+3,j+5
     $        ,k+1)+8.495389982295436D3*xyzi(i+3,j+1,k+1)-7.250948090357
     $        628D4*xyzi(i+1,j+11,k+1)+1.856242711131553D5*xyzi(i+1,j+9,
     $        k+1)-1.694830301467939D5*xyzi(i+1,j+7,k+1)+6.4564963865445
     $        31D4*xyzi(i+1,j+5,k+1)-8.495389982295436D3*xyzi(i+1,j+3,k+
     $        1))+angi
            angi=yrk(181)*(-5.701024531700933D4*xyzi(i+11,j+1,k+1)-2.850
     $        512265850467D5*xyzi(i+9,j+3,k+1)+1.596286868876261D5*xyzi(
     $        i+9,j+1,k+1)-5.701024531700933D5*xyzi(i+7,j+5,k+1)+6.38514
     $        7475505045D5*xyzi(i+7,j+3,k+1)-1.665690645783925D5*xyzi(i+
     $        7,j+1,k+1)-5.701024531700933D5*xyzi(i+5,j+7,k+1)+9.5777212
     $        13257568D5*xyzi(i+5,j+5,k+1)-4.997071937351775D5*xyzi(i+5,
     $        j+3,k+1)+7.93186021801869D4*xyzi(i+5,j+1,k+1)-2.8505122658
     $        50467D5*xyzi(i+3,j+9,k+1)+6.385147475505045D5*xyzi(i+3,j+7
     $        ,k+1)-4.997071937351775D5*xyzi(i+3,j+5,k+1)+1.586372043603
     $        738D5*xyzi(i+3,j+3,k+1)-1.669865309056566D4*xyzi(i+3,j+1,k
     $        +1)-5.701024531700933D4*xyzi(i+1,j+11,k+1)+1.5962868688762
     $        61D5*xyzi(i+1,j+9,k+1)-1.665690645783925D5*xyzi(i+1,j+7,k+
     $        1)+7.93186021801869D4*xyzi(i+1,j+5,k+1)-1.669865309056566D
     $        4*xyzi(i+1,j+3,k+1)+1.178728453451694D3*xyzi(i+1,j+1,k+1))
     $        +angi
          endif
          OMEGA(13)=angi
        endif
c...
c...  l=15
c...
        if(lmax.ge.15)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(242)*(1.279373446563299D5*xyzi(i+15,j,k)+8.95561412
     $        5943095D5*xyzi(i+13,j+2,k)-4.941028483278949D5*xyzi(i+13,j
     $        ,k)+2.686684237782929D6*xyzi(i+11,j+4,k)-2.964617089967369
     $        D6*xyzi(i+11,j+2,k)+7.68604430732281D5*xyzi(i+11,j,k)+4.47
     $        7807062971548D6*xyzi(i+9,j+6,k)-7.411542724918424D6*xyzi(i
     $        +9,j+4,k)+3.843022153661405D6*xyzi(i+9,j+2,k)-6.1488354458
     $        58248D5*xyzi(i+9,j,k)+4.477807062971548D6*xyzi(i+7,j+8,k)-
     $        9.882056966557898D6*xyzi(i+7,j+6,k)+7.68604430732281D6*xyz
     $        i(i+7,j+4,k)-2.459534178343299D6*xyzi(i+7,j+2,k)+2.6734067
     $        15590543D5*xyzi(i+7,j,k)+2.686684237782929D6*xyzi(i+5,j+10
     $        ,k)-7.411542724918424D6*xyzi(i+5,j+8,k)+7.68604430732281D6
     $        *xyzi(i+5,j+6,k)-3.689301267514949D6*xyzi(i+5,j+4,k)+8.020
     $        220146771627D5*xyzi(i+5,j+2,k)-6.110643921349812D4*xyzi(i+
     $        5,j,k)+8.955614125943095D5*xyzi(i+3,j+12,k)-2.964617089967
     $        369D6*xyzi(i+3,j+10,k)+3.843022153661405D6*xyzi(i+3,j+8,k)
     $        -2.459534178343299D6*xyzi(i+3,j+6,k)+8.020220146771627D5*x
     $        yzi(i+3,j+4,k)-1.222128784269962D5*xyzi(i+3,j+2,k)+6.43225
     $        6759315591D3*xyzi(i+3,j,k)+1.279373446563299D5*xyzi(i+1,j+
     $        14,k)-4.941028483278949D5*xyzi(i+1,j+12,k)+7.6860443073228
     $        1D5*xyzi(i+1,j+10,k)-6.148835445858248D5*xyzi(i+1,j+8,k)+2
     $        .673406715590543D5*xyzi(i+1,j+6,k)-6.110643921349812D4*xyz
     $        i(i+1,j+4,k)+6.432256759315591D3*xyzi(i+1,j+2,k)-2.1621031
     $        12374989D2*xyzi(i+1,j,k))+angi
            angi=yrk(244)*(-9.866708857725911D4*xyzi(i+15,j,k)-2.9600126
     $        57317773D5*xyzi(i+13,j+2,k)+3.674498471153098D5*xyzi(i+13,
     $        j,k)+2.960012657317773D5*xyzi(i+11,j+4,k)+7.34899694230619
     $        6D5*xyzi(i+11,j+2,k)-5.44370143874533D5*xyzi(i+11,j,k)+2.4
     $        66677214431478D6*xyzi(i+9,j+6,k)-1.837249235576549D6*xyzi(
     $        i+9,j+4,k)-5.44370143874533D5*xyzi(i+9,j+2,k)+4.0646304075
     $        96513D5*xyzi(i+9,j,k)+4.44001898597666D6*xyzi(i+7,j+8,k)-7
     $        .348996942306196D6*xyzi(i+7,j+6,k)+3.266220863247198D6*xyz
     $        i(i+7,j+4,k)-1.590507550798636D5*xyzi(i+7,j,k)+3.848016454
     $        513105D6*xyzi(i+5,j+10,k)-9.186246177882745D6*xyzi(i+5,j+8
     $        ,k)+7.621182014243462D6*xyzi(i+5,j+6,k)-2.438778244557908D
     $        6*xyzi(i+5,j+4,k)+1.590507550798636D5*xyzi(i+5,j+2,k)+3.02
     $        9538191997401D4*xyzi(i+5,j,k)+1.677340505813405D6*xyzi(i+3
     $        ,j+12,k)-5.144297859614337D6*xyzi(i+3,j+10,k)+5.9880715826
     $        19863D6*xyzi(i+3,j+8,k)-3.25170432607721D6*xyzi(i+3,j+6,k)
     $        +7.952537753993178D5*xyzi(i+3,j+4,k)-6.059076383994802D4*x
     $        yzi(i+3,j+2,k)-2.125991713682387D3*xyzi(i+3,j,k)+2.9600126
     $        57317773D5*xyzi(i+1,j+14,k)-1.102349541345929D6*xyzi(i+1,j
     $        +12,k)+1.633110431623599D6*xyzi(i+1,j+10,k)-1.219389122278
     $        954D6*xyzi(i+1,j+8,k)+4.771522652395907D5*xyzi(i+1,j+6,k)-
     $        9.088614575992203D4*xyzi(i+1,j+4,k)+6.37797514104716D3*xyz
     $        i(i+1,j+2,k))+angi
            angi=yrk(246)*(5.81523782519791D4*xyzi(i+15,j,k)-2.907618912
     $        598955D5*xyzi(i+13,j+2,k)-2.005254422482038D5*xyzi(i+13,j,
     $        k)-2.035333238819269D6*xyzi(i+11,j+4,k)+1.203152653489223D
     $        6*xyzi(i+11,j+2,k)+2.673672563309384D5*xyzi(i+11,j,k)-3.77
     $        9904586378641D6*xyzi(i+9,j+6,k)+5.81523782519791D6*xyzi(i+
     $        9,j+4,k)-1.871570794316569D6*xyzi(i+9,j+2,k)-1.71115044051
     $        8006D5*xyzi(i+9,j,k)-2.61685702133906D6*xyzi(i+7,j+8,k)+7.
     $        218915920935337D6*xyzi(i+7,j+6,k)-5.882079639280644D6*xyzi
     $        (i+7,j+4,k)+1.368920352414405D6*xyzi(i+7,j+2,k)+5.20784916
     $        679393D4*xyzi(i+7,j,k)+5.81523782519791D4*xyzi(i+5,j+10,k)
     $        +1.804728980233834D6*xyzi(i+5,j+8,k)-3.743141588633137D6*x
     $        yzi(i+5,j+6,k)+2.395610616725208D6*xyzi(i+5,j+4,k)-4.68706
     $        4250114537D5*xyzi(i+5,j+2,k)-5.951827619193063D3*xyzi(i+5,
     $        j,k)+8.722856737796865D5*xyzi(i+3,j+12,k)-2.00525442248203
     $        8D6*xyzi(i+3,j+10,k)+1.336836281654692D6*xyzi(i+3,j+8,k)-2
     $        .603924583396965D5*xyzi(i+3,j+4,k)+5.951827619193063D4*xyz
     $        i(i+3,j+2,k)+2.907618912598955D5*xyzi(i+1,j+14,k)-1.002627
     $        211241019D6*xyzi(i+1,j+12,k)+1.336836281654692D6*xyzi(i+1,
     $        j+10,k)-8.555752202590029D5*xyzi(i+1,j+8,k)+2.603924583396
     $        965D5*xyzi(i+1,j+6,k)-2.975913809596532D4*xyzi(i+1,j+4,k))
     $        +angi
            angi=yrk(248)*(-2.566656485077824D4*xyzi(i+15,j,k)+4.3633160
     $        246323D5*xyzi(i+13,j+2,k)+7.788474851270638D4*xyzi(i+13,j,
     $        k)+1.103662288583464D6*xyzi(i+11,j+4,k)-1.401925473228715D
     $        6*xyzi(i+11,j+2,k)-8.653860945856264D4*xyzi(i+11,j,k)-2.82
     $        3322133585606D5*xyzi(i+9,j+6,k)-1.947118712817659D6*xyzi(i
     $        +9,j+4,k)+1.64423357971269D6*xyzi(i+9,j+2,k)+4.15385325401
     $        1007D4*xyzi(i+9,j,k)-2.540989920227045D6*xyzi(i+7,j+8,k)+2
     $        .80385094645743D6*xyzi(i+7,j+6,k)+5.192316567513758D5*xyzi
     $        (i+7,j+4,k)-8.307706508022013D5*xyzi(i+7,j+2,k)-7.22409261
     $        5671316D3*xyzi(i+7,j,k)-1.976325493509924D6*xyzi(i+5,j+10,
     $        k)+4.906739156300502D6*xyzi(i+5,j+8,k)-3.634621597259631D6
     $        *xyzi(i+5,j+6,k)+5.815394555615409D5*xyzi(i+5,j+4,k)+1.517
     $        059449290976D5*xyzi(i+5,j+2,k)-1.796659539554477D5*xyzi(i+
     $        3,j+12,k)+1.090386479177889D6*xyzi(i+3,j+10,k)-1.817310798
     $        629815D6*xyzi(i+3,j+8,k)+1.163078911123082D6*xyzi(i+3,j+6,
     $        k)-2.528432415484961D5*xyzi(i+3,j+4,k)+1.796659539554477D5
     $        *xyzi(i+1,j+14,k)-5.451932395889446D5*xyzi(i+1,j+12,k)+6.0
     $        57702662099385D5*xyzi(i+1,j+10,k)-2.907697277807705D5*xyzi
     $        (i+1,j+8,k)+5.056864830969921D4*xyzi(i+1,j+6,k))+angi
            angi=yrk(250)*(8.17508397215609D3*xyzi(i+15,j,k)-2.697777710
     $        81151D5*xyzi(i+13,j+2,k)-2.029676020673236D4*xyzi(i+13,j,k
     $        )+1.716767634152779D5*xyzi(i+11,j+4,k)+6.900898470289003D5
     $        *xyzi(i+11,j+2,k)+1.653810090918933D4*xyzi(i+11,j,k)+1.528
     $        740702793189D6*xyzi(i+9,j+6,k)-1.11632181137028D6*xyzi(i+9
     $        ,j+4,k)-5.788335318216266D5*xyzi(i+9,j+2,k)-4.410160242450
     $        489D3*xyzi(i+9,j,k)+8.093333132434529D5*xyzi(i+7,j+8,k)-2.
     $        679172347288672D6*xyzi(i+7,j+6,k)+1.48842908182704D6*xyzi(
     $        i+7,j+4,k)+1.587657687282176D5*xyzi(i+7,j+2,k)-8.093333132
     $        434529D5*xyzi(i+5,j+10,k)+6.697930868221679D5*xyzi(i+5,j+8
     $        ,k)+6.946002381859519D5*xyzi(i+5,j+6,k)-5.556801905487616D
     $        5*xyzi(i+5,j+4,k)-4.659797864128971D5*xyzi(i+3,j+12,k)+1.3
     $        39586173644336D6*xyzi(i+3,j+10,k)-1.2403575681892D6*xyzi(i
     $        +3,j+8,k)+3.70453460365841D5*xyzi(i+3,j+6,k)+7.35757557494
     $        0481D4*xyzi(i+1,j+14,k)-1.826708418605913D5*xyzi(i+1,j+12,
     $        k)+1.48842908182704D5*xyzi(i+1,j+10,k)-3.96914421820544D4*
     $        xyzi(i+1,j+8,k))+angi
            angi=yrk(252)*(-1.756289768694704D3*xyzi(i+15,j,k)+9.3083357
     $        74081929D4*xyzi(i+13,j+2,k)+3.149209240418089D3*xyzi(i+13,
     $        j,k)-3.881400388815295D5*xyzi(i+11,j+4,k)-1.70057298982576
     $        8D5*xyzi(i+11,j+2,k)-1.399648551296928D3*xyzi(i+11,j,k)-2.
     $        511494369233426D5*xyzi(i+9,j+6,k)+8.660325411149745D5*xyzi
     $        (i+9,j+4,k)+7.698067032133107D4*xyzi(i+9,j+2,k)+7.53448310
     $        7700278D5*xyzi(i+7,j+8,k)-4.156956197351878D5*xyzi(i+7,j+6
     $        ,k)-4.618840219279864D5*xyzi(i+7,j+4,k)+2.511494369233426D
     $        5*xyzi(i+5,j+10,k)-9.353151444041725D5*xyzi(i+5,j+8,k)+6.4
     $        6637630699181D5*xyzi(i+5,j+6,k)-2.511494369233426D5*xyzi(i
     $        +3,j+12,k)+4.849782230243857D5*xyzi(i+3,j+10,k)-2.30942010
     $        9639932D5*xyzi(i+3,j+8,k)+1.931918745564174D4*xyzi(i+1,j+1
     $        4,k)-3.464130164459898D4*xyzi(i+1,j+12,k)+1.53961340642662
     $        1D4*xyzi(i+1,j+10,k))+angi
            angi=yrk(254)*(2.212717122920637D2*xyzi(i+15,j,k)-1.70379218
     $        4648891D4*xyzi(i+13,j+2,k)-2.136416532475098D2*xyzi(i+13,j
     $        ,k)+1.409500807300446D5*xyzi(i+11,j+4,k)+1.666404895330577
     $        D4*xyzi(i+11,j+2,k)-2.214929840043558D5*xyzi(i+9,j+6,k)-1.
     $        527537820719695D5*xyzi(i+9,j+4,k)-9.492556457329534D4*xyzi
     $        (i+7,j+8,k)+3.666090769727268D5*xyzi(i+7,j+6,k)+2.21492984
     $        0043558D5*xyzi(i+5,j+10,k)-2.749568077295451D5*xyzi(i+5,j+
     $        8,k)-6.04071774557334D4*xyzi(i+3,j+12,k)+6.110151282878781
     $        D4*xyzi(i+3,j+10,k)+2.876532259796829D3*xyzi(i+1,j+14,k)-2
     $        .777341492217627D3*xyzi(i+1,j+12,k))+angi
            angi=yrk(256)*(-1.060916657008769D1*xyzi(i+15,j,k)+1.1139624
     $        89859208D3*xyzi(i+13,j+2,k)-1.44815123681697D4*xyzi(i+11,j
     $        +4,k)+5.309887868328891D4*xyzi(i+9,j+6,k)-6.82699868785143
     $        1D4*xyzi(i+7,j+8,k)+3.185932720997334D4*xyzi(i+5,j+10,k)-4
     $        .827170789389901D3*xyzi(i+3,j+12,k)+1.591374985513154D2*xy
     $        zi(i+1,j+14,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(226)*(-1.591374985513154D2*xyzi(i+14,j+1,k)+4.82717
     $        0789389901D3*xyzi(i+12,j+3,k)-3.185932720997334D4*xyzi(i+1
     $        0,j+5,k)+6.826998687851431D4*xyzi(i+8,j+7,k)-5.30988786832
     $        8891D4*xyzi(i+6,j+9,k)+1.44815123681697D4*xyzi(i+4,j+11,k)
     $        -1.113962489859208D3*xyzi(i+2,j+13,k)+1.060916657008769D1*
     $        xyzi(i,j+15,k))+angi
            angi=yrk(228)*(2.876532259796829D3*xyzi(i+14,j+1,k)-6.040717
     $        74557334D4*xyzi(i+12,j+3,k)-2.777341492217627D3*xyzi(i+12,
     $        j+1,k)+2.214929840043558D5*xyzi(i+10,j+5,k)+6.110151282878
     $        781D4*xyzi(i+10,j+3,k)-9.492556457329534D4*xyzi(i+8,j+7,k)
     $        -2.749568077295451D5*xyzi(i+8,j+5,k)-2.214929840043558D5*x
     $        yzi(i+6,j+9,k)+3.666090769727268D5*xyzi(i+6,j+7,k)+1.40950
     $        0807300446D5*xyzi(i+4,j+11,k)-1.527537820719695D5*xyzi(i+4
     $        ,j+9,k)-1.703792184648891D4*xyzi(i+2,j+13,k)+1.66640489533
     $        0577D4*xyzi(i+2,j+11,k)+2.212717122920637D2*xyzi(i,j+15,k)
     $        -2.136416532475098D2*xyzi(i,j+13,k))+angi
            angi=yrk(230)*(-1.931918745564174D4*xyzi(i+14,j+1,k)+2.51149
     $        4369233426D5*xyzi(i+12,j+3,k)+3.464130164459898D4*xyzi(i+1
     $        2,j+1,k)-2.511494369233426D5*xyzi(i+10,j+5,k)-4.8497822302
     $        43857D5*xyzi(i+10,j+3,k)-1.539613406426621D4*xyzi(i+10,j+1
     $        ,k)-7.534483107700278D5*xyzi(i+8,j+7,k)+9.353151444041725D
     $        5*xyzi(i+8,j+5,k)+2.309420109639932D5*xyzi(i+8,j+3,k)+2.51
     $        1494369233426D5*xyzi(i+6,j+9,k)+4.156956197351878D5*xyzi(i
     $        +6,j+7,k)-6.46637630699181D5*xyzi(i+6,j+5,k)+3.88140038881
     $        5295D5*xyzi(i+4,j+11,k)-8.660325411149745D5*xyzi(i+4,j+9,k
     $        )+4.618840219279864D5*xyzi(i+4,j+7,k)-9.308335774081929D4*
     $        xyzi(i+2,j+13,k)+1.700572989825768D5*xyzi(i+2,j+11,k)-7.69
     $        8067032133107D4*xyzi(i+2,j+9,k)+1.756289768694704D3*xyzi(i
     $        ,j+15,k)-3.149209240418089D3*xyzi(i,j+13,k)+1.399648551296
     $        928D3*xyzi(i,j+11,k))+angi
            angi=yrk(232)*(7.357575574940481D4*xyzi(i+14,j+1,k)-4.659797
     $        864128971D5*xyzi(i+12,j+3,k)-1.826708418605913D5*xyzi(i+12
     $        ,j+1,k)-8.093333132434529D5*xyzi(i+10,j+5,k)+1.33958617364
     $        4336D6*xyzi(i+10,j+3,k)+1.48842908182704D5*xyzi(i+10,j+1,k
     $        )+8.093333132434529D5*xyzi(i+8,j+7,k)+6.697930868221679D5*
     $        xyzi(i+8,j+5,k)-1.2403575681892D6*xyzi(i+8,j+3,k)-3.969144
     $        21820544D4*xyzi(i+8,j+1,k)+1.528740702793189D6*xyzi(i+6,j+
     $        9,k)-2.679172347288672D6*xyzi(i+6,j+7,k)+6.946002381859519
     $        D5*xyzi(i+6,j+5,k)+3.70453460365841D5*xyzi(i+6,j+3,k)+1.71
     $        6767634152779D5*xyzi(i+4,j+11,k)-1.11632181137028D6*xyzi(i
     $        +4,j+9,k)+1.48842908182704D6*xyzi(i+4,j+7,k)-5.55680190548
     $        7616D5*xyzi(i+4,j+5,k)-2.69777771081151D5*xyzi(i+2,j+13,k)
     $        +6.900898470289003D5*xyzi(i+2,j+11,k)-5.788335318216266D5*
     $        xyzi(i+2,j+9,k)+1.587657687282176D5*xyzi(i+2,j+7,k)+8.1750
     $        8397215609D3*xyzi(i,j+15,k)-2.029676020673236D4*xyzi(i,j+1
     $        3,k)+1.653810090918933D4*xyzi(i,j+11,k)-4.410160242450489D
     $        3*xyzi(i,j+9,k))+angi
            angi=yrk(234)*(-1.796659539554477D5*xyzi(i+14,j+1,k)+1.79665
     $        9539554477D5*xyzi(i+12,j+3,k)+5.451932395889446D5*xyzi(i+1
     $        2,j+1,k)+1.976325493509924D6*xyzi(i+10,j+5,k)-1.0903864791
     $        77889D6*xyzi(i+10,j+3,k)-6.057702662099385D5*xyzi(i+10,j+1
     $        ,k)+2.540989920227045D6*xyzi(i+8,j+7,k)-4.906739156300502D
     $        6*xyzi(i+8,j+5,k)+1.817310798629815D6*xyzi(i+8,j+3,k)+2.90
     $        7697277807705D5*xyzi(i+8,j+1,k)+2.823322133585606D5*xyzi(i
     $        +6,j+9,k)-2.80385094645743D6*xyzi(i+6,j+7,k)+3.63462159725
     $        9631D6*xyzi(i+6,j+5,k)-1.163078911123082D6*xyzi(i+6,j+3,k)
     $        -5.056864830969921D4*xyzi(i+6,j+1,k)-1.103662288583464D6*x
     $        yzi(i+4,j+11,k)+1.947118712817659D6*xyzi(i+4,j+9,k)-5.1923
     $        16567513758D5*xyzi(i+4,j+7,k)-5.815394555615409D5*xyzi(i+4
     $        ,j+5,k)+2.528432415484961D5*xyzi(i+4,j+3,k)-4.363316024632
     $        3D5*xyzi(i+2,j+13,k)+1.401925473228715D6*xyzi(i+2,j+11,k)-
     $        1.64423357971269D6*xyzi(i+2,j+9,k)+8.307706508022013D5*xyz
     $        i(i+2,j+7,k)-1.517059449290976D5*xyzi(i+2,j+5,k)+2.5666564
     $        85077824D4*xyzi(i,j+15,k)-7.788474851270638D4*xyzi(i,j+13,
     $        k)+8.653860945856264D4*xyzi(i,j+11,k)-4.153853254011007D4*
     $        xyzi(i,j+9,k)+7.224092615671316D3*xyzi(i,j+7,k))+angi
            angi=yrk(236)*(2.907618912598955D5*xyzi(i+14,j+1,k)+8.722856
     $        737796865D5*xyzi(i+12,j+3,k)-1.002627211241019D6*xyzi(i+12
     $        ,j+1,k)+5.81523782519791D4*xyzi(i+10,j+5,k)-2.005254422482
     $        038D6*xyzi(i+10,j+3,k)+1.336836281654692D6*xyzi(i+10,j+1,k
     $        )-2.61685702133906D6*xyzi(i+8,j+7,k)+1.804728980233834D6*x
     $        yzi(i+8,j+5,k)+1.336836281654692D6*xyzi(i+8,j+3,k)-8.55575
     $        2202590029D5*xyzi(i+8,j+1,k)-3.779904586378641D6*xyzi(i+6,
     $        j+9,k)+7.218915920935337D6*xyzi(i+6,j+7,k)-3.7431415886331
     $        37D6*xyzi(i+6,j+5,k)+2.603924583396965D5*xyzi(i+6,j+1,k)-2
     $        .035333238819269D6*xyzi(i+4,j+11,k)+5.81523782519791D6*xyz
     $        i(i+4,j+9,k)-5.882079639280644D6*xyzi(i+4,j+7,k)+2.3956106
     $        16725208D6*xyzi(i+4,j+5,k)-2.603924583396965D5*xyzi(i+4,j+
     $        3,k)-2.975913809596532D4*xyzi(i+4,j+1,k)-2.907618912598955
     $        D5*xyzi(i+2,j+13,k)+1.203152653489223D6*xyzi(i+2,j+11,k)-1
     $        .871570794316569D6*xyzi(i+2,j+9,k)+1.368920352414405D6*xyz
     $        i(i+2,j+7,k)-4.687064250114537D5*xyzi(i+2,j+5,k)+5.9518276
     $        19193063D4*xyzi(i+2,j+3,k)+5.81523782519791D4*xyzi(i,j+15,
     $        k)-2.005254422482038D5*xyzi(i,j+13,k)+2.673672563309384D5*
     $        xyzi(i,j+11,k)-1.711150440518006D5*xyzi(i,j+9,k)+5.2078491
     $        6679393D4*xyzi(i,j+7,k)-5.951827619193063D3*xyzi(i,j+5,k))
     $        +angi
            angi=yrk(238)*(-2.960012657317773D5*xyzi(i+14,j+1,k)-1.67734
     $        0505813405D6*xyzi(i+12,j+3,k)+1.102349541345929D6*xyzi(i+1
     $        2,j+1,k)-3.848016454513105D6*xyzi(i+10,j+5,k)+5.1442978596
     $        14337D6*xyzi(i+10,j+3,k)-1.633110431623599D6*xyzi(i+10,j+1
     $        ,k)-4.44001898597666D6*xyzi(i+8,j+7,k)+9.186246177882745D6
     $        *xyzi(i+8,j+5,k)-5.988071582619863D6*xyzi(i+8,j+3,k)+1.219
     $        389122278954D6*xyzi(i+8,j+1,k)-2.466677214431478D6*xyzi(i+
     $        6,j+9,k)+7.348996942306196D6*xyzi(i+6,j+7,k)-7.62118201424
     $        3462D6*xyzi(i+6,j+5,k)+3.25170432607721D6*xyzi(i+6,j+3,k)-
     $        4.771522652395907D5*xyzi(i+6,j+1,k)-2.960012657317773D5*xy
     $        zi(i+4,j+11,k)+1.837249235576549D6*xyzi(i+4,j+9,k)-3.26622
     $        0863247198D6*xyzi(i+4,j+7,k)+2.438778244557908D6*xyzi(i+4,
     $        j+5,k)-7.952537753993178D5*xyzi(i+4,j+3,k)+9.0886145759922
     $        03D4*xyzi(i+4,j+1,k)+2.960012657317773D5*xyzi(i+2,j+13,k)-
     $        7.348996942306196D5*xyzi(i+2,j+11,k)+5.44370143874533D5*xy
     $        zi(i+2,j+9,k)-1.590507550798636D5*xyzi(i+2,j+5,k)+6.059076
     $        383994802D4*xyzi(i+2,j+3,k)-6.37797514104716D3*xyzi(i+2,j+
     $        1,k)+9.866708857725911D4*xyzi(i,j+15,k)-3.674498471153098D
     $        5*xyzi(i,j+13,k)+5.44370143874533D5*xyzi(i,j+11,k)-4.06463
     $        0407596513D5*xyzi(i,j+9,k)+1.590507550798636D5*xyzi(i,j+7,
     $        k)-3.029538191997401D4*xyzi(i,j+5,k)+2.125991713682387D3*x
     $        yzi(i,j+3,k))+angi
            angi=yrk(240)*(1.279373446563299D5*xyzi(i+14,j+1,k)+8.955614
     $        125943095D5*xyzi(i+12,j+3,k)-4.941028483278949D5*xyzi(i+12
     $        ,j+1,k)+2.686684237782929D6*xyzi(i+10,j+5,k)-2.96461708996
     $        7369D6*xyzi(i+10,j+3,k)+7.68604430732281D5*xyzi(i+10,j+1,k
     $        )+4.477807062971548D6*xyzi(i+8,j+7,k)-7.411542724918424D6*
     $        xyzi(i+8,j+5,k)+3.843022153661405D6*xyzi(i+8,j+3,k)-6.1488
     $        35445858248D5*xyzi(i+8,j+1,k)+4.477807062971548D6*xyzi(i+6
     $        ,j+9,k)-9.882056966557898D6*xyzi(i+6,j+7,k)+7.686044307322
     $        81D6*xyzi(i+6,j+5,k)-2.459534178343299D6*xyzi(i+6,j+3,k)+2
     $        .673406715590543D5*xyzi(i+6,j+1,k)+2.686684237782929D6*xyz
     $        i(i+4,j+11,k)-7.411542724918424D6*xyzi(i+4,j+9,k)+7.686044
     $        30732281D6*xyzi(i+4,j+7,k)-3.689301267514949D6*xyzi(i+4,j+
     $        5,k)+8.020220146771627D5*xyzi(i+4,j+3,k)-6.110643921349812
     $        D4*xyzi(i+4,j+1,k)+8.955614125943095D5*xyzi(i+2,j+13,k)-2.
     $        964617089967369D6*xyzi(i+2,j+11,k)+3.843022153661405D6*xyz
     $        i(i+2,j+9,k)-2.459534178343299D6*xyzi(i+2,j+7,k)+8.0202201
     $        46771627D5*xyzi(i+2,j+5,k)-1.222128784269962D5*xyzi(i+2,j+
     $        3,k)+6.432256759315591D3*xyzi(i+2,j+1,k)+1.279373446563299
     $        D5*xyzi(i,j+15,k)-4.941028483278949D5*xyzi(i,j+13,k)+7.686
     $        04430732281D5*xyzi(i,j+11,k)-6.148835445858248D5*xyzi(i,j+
     $        9,k)+2.673406715590543D5*xyzi(i,j+7,k)-6.110643921349812D4
     $        *xyzi(i,j+5,k)+6.432256759315591D3*xyzi(i,j+3,k)-2.1621031
     $        12374989D2*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(241)*(9.343222615411324D4*xyzi(i,j,k+15)-3.38289094
     $        6959272D5*xyzi(i,j,k+13)+4.886398034496727D5*xyzi(i,j,k+11
     $        )-3.583358558630933D5*xyzi(i,j,k+9)+1.402183783812104D5*xy
     $        zi(i,j,k+7)-2.804367567624208D4*xyzi(i,j,k+5)+2.4599715505
     $        47551D3*xyzi(i,j,k+3)-6.201608950960213D1*xyzi(i,j,k+1))+a
     $        ngi
            angi=yrk(243)*(1.161012484626535D5*xyzi(i+14,j,k+1)+5.805062
     $        423132676D5*xyzi(i+12,j+2,k+1)-3.843351673246462D5*xyzi(i+
     $        12,j,k+1)+1.044911236163882D6*xyzi(i+10,j+4,k+1)-1.5373406
     $        69298585D6*xyzi(i+10,j+2,k+1)+4.982122539393561D5*xyzi(i+1
     $        0,j,k+1)+5.805062423132676D5*xyzi(i+8,j+6,k+1)-1.921675836
     $        623231D6*xyzi(i+8,j+4,k+1)+1.494636761818068D6*xyzi(i+8,j+
     $        2,k+1)-3.188558425211879D5*xyzi(i+8,j,k+1)-5.8050624231326
     $        76D5*xyzi(i+6,j+8,k+1)+9.964245078787123D5*xyzi(i+6,j+4,k+
     $        1)-6.377116850423759D5*xyzi(i+6,j+2,k+1)+1.039747312569091
     $        D5*xyzi(i+6,j,k+1)-1.044911236163882D6*xyzi(i+4,j+10,k+1)+
     $        1.921675836623231D6*xyzi(i+4,j+8,k+1)-9.964245078787123D5*
     $        xyzi(i+4,j+6,k+1)+1.039747312569091D5*xyzi(i+4,j+2,k+1)-1.
     $        584376857248139D4*xyzi(i+4,j,k+1)-5.805062423132676D5*xyzi
     $        (i+2,j+12,k+1)+1.537340669298585D6*xyzi(i+2,j+10,k+1)-1.49
     $        4636761818068D6*xyzi(i+2,j+8,k+1)+6.377116850423759D5*xyzi
     $        (i+2,j+6,k+1)-1.039747312569091D5*xyzi(i+2,j+4,k+1)+8.3388
     $        25564463888D2*xyzi(i+2,j,k+1)-1.161012484626535D5*xyzi(i,j
     $        +14,k+1)+3.843351673246462D5*xyzi(i,j+12,k+1)-4.9821225393
     $        93561D5*xyzi(i,j+10,k+1)+3.188558425211879D5*xyzi(i,j+8,k+
     $        1)-1.039747312569091D5*xyzi(i,j+6,k+1)+1.584376857248139D4
     $        *xyzi(i,j+4,k+1)-8.338825564463888D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(245)*(-7.841265083878948D4*xyzi(i+14,j,k+1)+7.84126
     $        5083878948D4*xyzi(i+12,j+2,k+1)+2.433496060514156D5*xyzi(i
     $        +12,j,k+1)+1.489840365937D6*xyzi(i+10,j+4,k+1)-4.866992121
     $        028313D5*xyzi(i+10,j+2,k+1)-2.884143479127889D5*xyzi(i+10,
     $        j,k+1)+3.528569287745527D6*xyzi(i+8,j+6,k+1)-4.13694330287
     $        4065D6*xyzi(i+8,j+4,k+1)+8.652430437383667D5*xyzi(i+8,j+2,
     $        k+1)+1.615120348311618D5*xyzi(i+8,j,k+1)+3.528569287745527
     $        D6*xyzi(i+6,j+8,k+1)-6.813788969439638D6*xyzi(i+6,j+6,k+1)
     $        +4.037800870779044D6*xyzi(i+6,j+4,k+1)-6.460481393246471D5
     $        *xyzi(i+6,j+2,k+1)-4.213357430378134D4*xyzi(i+6,j,k+1)+1.4
     $        89840365937D6*xyzi(i+4,j+10,k+1)-4.136943302874065D6*xyzi(
     $        i+4,j+8,k+1)+4.037800870779044D6*xyzi(i+4,j+6,k+1)-1.61512
     $        0348311618D6*xyzi(i+4,j+4,k+1)+2.106678715189067D5*xyzi(i+
     $        4,j+2,k+1)+4.012721362264889D3*xyzi(i+4,j,k+1)+7.841265083
     $        878948D4*xyzi(i+2,j+12,k+1)-4.866992121028313D5*xyzi(i+2,j
     $        +10,k+1)+8.652430437383667D5*xyzi(i+2,j+8,k+1)-6.460481393
     $        246471D5*xyzi(i+2,j+6,k+1)+2.106678715189067D5*xyzi(i+2,j+
     $        4,k+1)-2.407632817358933D4*xyzi(i+2,j+2,k+1)-7.84126508387
     $        8948D4*xyzi(i,j+14,k+1)+2.433496060514156D5*xyzi(i,j+12,k+
     $        1)-2.884143479127889D5*xyzi(i,j+10,k+1)+1.615120348311618D
     $        5*xyzi(i,j+8,k+1)-4.213357430378134D4*xyzi(i,j+6,k+1)+4.01
     $        2721362264889D3*xyzi(i,j+4,k+1))+angi
            angi=yrk(247)*(4.012895342554011D4*xyzi(i+14,j,k+1)-4.414184
     $        876809412D5*xyzi(i+12,j+2,k+1)-1.107005611739038D5*xyzi(i+
     $        12,j,k+1)-1.565029183596064D6*xyzi(i+10,j+4,k+1)+1.3284067
     $        34086845D6*xyzi(i+10,j+2,k+1)+1.107005611739038D5*xyzi(i+1
     $        0,j,k+1)-1.083481742489583D6*xyzi(i+8,j+6,k+1)+2.988915151
     $        695401D6*xyzi(i+8,j+4,k+1)-1.439107295260749D6*xyzi(i+8,j+
     $        2,k+1)-4.723223943419894D4*xyzi(i+8,j,k+1)+1.0834817424895
     $        83D6*xyzi(i+6,j+8,k+1)-1.549807856434653D6*xyzi(i+6,j+4,k+
     $        1)+6.612513520787851D5*xyzi(i+6,j+2,k+1)+7.187514696508534
     $        D3*xyzi(i+6,j,k+1)+1.565029183596064D6*xyzi(i+4,j+10,k+1)-
     $        2.988915151695401D6*xyzi(i+4,j+8,k+1)+1.549807856434653D6*
     $        xyzi(i+4,j+6,k+1)-1.07812720447628D5*xyzi(i+4,j+2,k+1)+4.4
     $        14184876809412D5*xyzi(i+2,j+12,k+1)-1.328406734086845D6*xy
     $        zi(i+2,j+10,k+1)+1.439107295260749D6*xyzi(i+2,j+8,k+1)-6.6
     $        12513520787851D5*xyzi(i+2,j+6,k+1)+1.07812720447628D5*xyzi
     $        (i+2,j+4,k+1)-4.012895342554011D4*xyzi(i,j+14,k+1)+1.10700
     $        5611739038D5*xyzi(i,j+12,k+1)-1.107005611739038D5*xyzi(i,j
     $        +10,k+1)+4.723223943419894D4*xyzi(i,j+8,k+1)-7.18751469650
     $        8534D3*xyzi(i,j+6,k+1))+angi
            angi=yrk(249)*(-1.513731411750108D4*xyzi(i+14,j,k+1)+3.78432
     $        8529375269D5*xyzi(i+12,j+2,k+1)+3.445043902603693D4*xyzi(i
     $        +12,j,k+1)+1.665104552925118D5*xyzi(i+10,j+4,k+1)-8.957114
     $        146769602D5*xyzi(i+10,j+2,k+1)-2.551884372299032D4*xyzi(i+
     $        10,j,k+1)-1.498594097632606D6*xyzi(i+8,j+6,k+1)+5.16756585
     $        390554D5*xyzi(i+8,j+4,k+1)+6.890087805207386D5*xyzi(i+8,j+
     $        2,k+1)+6.124522493517676D3*xyzi(i+8,j,k+1)-1.4985940976326
     $        06D6*xyzi(i+6,j+8,k+1)+2.893836878187102D6*xyzi(i+6,j+6,k+
     $        1)-1.071791436365593D6*xyzi(i+6,j+4,k+1)-1.714866298184949
     $        D5*xyzi(i+6,j+2,k+1)+1.665104552925118D5*xyzi(i+4,j+10,k+1
     $        )+5.16756585390554D5*xyzi(i+4,j+8,k+1)-1.071791436365593D6
     $        *xyzi(i+4,j+6,k+1)+4.287165745462373D5*xyzi(i+4,j+4,k+1)+3
     $        .784328529375269D5*xyzi(i+2,j+12,k+1)-8.957114146769602D5*
     $        xyzi(i+2,j+10,k+1)+6.890087805207386D5*xyzi(i+2,j+8,k+1)-1
     $        .714866298184949D5*xyzi(i+2,j+6,k+1)-1.513731411750108D4*x
     $        yzi(i,j+14,k+1)+3.445043902603693D4*xyzi(i,j+12,k+1)-2.551
     $        884372299032D4*xyzi(i,j+10,k+1)+6.124522493517676D3*xyzi(i
     $        ,j+8,k+1))+angi
            angi=yrk(251)*(4.004956867237501D3*xyzi(i+14,j,k+1)-1.722131
     $        452912125D5*xyzi(i+12,j+2,k+1)-6.62889412508276D3*xyzi(i+1
     $        2,j,k+1)+4.845997809357376D5*xyzi(i+10,j+4,k+1)+2.91671341
     $        5036414D5*xyzi(i+10,j+2,k+1)+2.700660569478161D3*xyzi(i+10
     $        ,j,k+1)+6.608178830941876D5*xyzi(i+8,j+6,k+1)-1.0937675306
     $        38655D6*xyzi(i+8,j+4,k+1)-1.215297256265173D5*xyzi(i+8,j+2
     $        ,k+1)-6.608178830941876D5*xyzi(i+6,j+8,k+1)+5.671387195904
     $        139D5*xyzi(i+6,j+4,k+1)-4.845997809357376D5*xyzi(i+4,j+10,
     $        k+1)+1.093767530638655D6*xyzi(i+4,j+8,k+1)-5.6713871959041
     $        39D5*xyzi(i+4,j+6,k+1)+1.722131452912125D5*xyzi(i+2,j+12,k
     $        +1)-2.916713415036414D5*xyzi(i+2,j+10,k+1)+1.2152972562651
     $        73D5*xyzi(i+2,j+8,k+1)-4.004956867237501D3*xyzi(i,j+14,k+1
     $        )+6.62889412508276D3*xyzi(i,j+12,k+1)-2.700660569478161D3*
     $        xyzi(i,j+10,k+1))+angi
            angi=yrk(253)*(-6.759962471539151D2*xyzi(i+14,j,k+1)+4.39397
     $        5606500448D4*xyzi(i+12,j+2,k+1)+6.060656008966135D2*xyzi(i
     $        +12,j,k+1)-2.900023900290296D5*xyzi(i+10,j+4,k+1)-4.000032
     $        965917649D4*xyzi(i+10,j+2,k+1)+2.900023900290296D5*xyzi(i+
     $        8,j+6,k+1)+3.000024724438237D5*xyzi(i+8,j+4,k+1)+2.9000239
     $        00290296D5*xyzi(i+6,j+8,k+1)-5.600046152284709D5*xyzi(i+6,
     $        j+6,k+1)-2.900023900290296D5*xyzi(i+4,j+10,k+1)+3.00002472
     $        4438237D5*xyzi(i+4,j+8,k+1)+4.393975606500448D4*xyzi(i+2,j
     $        +12,k+1)-4.000032965917649D4*xyzi(i+2,j+10,k+1)-6.75996247
     $        1539151D2*xyzi(i,j+14,k+1)+6.060656008966135D2*xyzi(i,j+12
     $        ,k+1))+angi
            angi=yrk(255)*(5.810879846766743D1*xyzi(i+14,j,k+1)-5.287900
     $        660557736D3*xyzi(i+12,j+2,k+1)+5.816690726613509D4*xyzi(i+
     $        10,j+4,k+1)-1.745007217984053D5*xyzi(i+8,j+6,k+1)+1.745007
     $        217984053D5*xyzi(i+6,j+8,k+1)-5.816690726613509D4*xyzi(i+4
     $        ,j+10,k+1)+5.287900660557736D3*xyzi(i+2,j+12,k+1)-5.810879
     $        846766743D1*xyzi(i,j+14,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(227)*(8.13523178547344D2*xyzi(i+13,j+1,k+1)-2.11516
     $        0264223094D4*xyzi(i+11,j+3,k+1)+1.163338145322702D5*xyzi(i
     $        +9,j+5,k+1)-1.994293963410346D5*xyzi(i+7,j+7,k+1)+1.163338
     $        145322702D5*xyzi(i+5,j+9,k+1)-2.115160264223094D4*xyzi(i+3
     $        ,j+11,k+1)+8.13523178547344D2*xyzi(i+1,j+13,k+1))+angi
            angi=yrk(229)*(-8.111954965846981D3*xyzi(i+13,j+1,k+1)+1.406
     $        072194080143D5*xyzi(i+11,j+3,k+1)+7.272787210759363D3*xyzi
     $        (i+11,j+1,k+1)-3.866698533720394D5*xyzi(i+9,j+5,k+1)-1.333
     $        34432197255D5*xyzi(i+9,j+3,k+1)+4.800039559101179D5*xyzi(i
     $        +7,j+5,k+1)+3.866698533720394D5*xyzi(i+5,j+9,k+1)-4.800039
     $        559101179D5*xyzi(i+5,j+7,k+1)-1.406072194080143D5*xyzi(i+3
     $        ,j+11,k+1)+1.33334432197255D5*xyzi(i+3,j+9,k+1)+8.11195496
     $        5846981D3*xyzi(i+1,j+13,k+1)-7.272787210759363D3*xyzi(i+1,
     $        j+11,k+1))+angi
            angi=yrk(231)*(4.004956867237501D4*xyzi(i+13,j+1,k+1)-4.0049
     $        56867237501D5*xyzi(i+11,j+3,k+1)-6.628894125082759D4*xyzi(
     $        i+11,j+1,k+1)+8.810905107922501D4*xyzi(i+9,j+5,k+1)+7.2917
     $        83537591036D5*xyzi(i+9,j+3,k+1)+2.700660569478161D4*xyzi(i
     $        +9,j+1,k+1)+1.0573086129507D6*xyzi(i+7,j+7,k+1)-8.75014024
     $        5109242D5*xyzi(i+7,j+5,k+1)-3.240792683373793D5*xyzi(i+7,j
     $        +3,k+1)+8.810905107922501D4*xyzi(i+5,j+9,k+1)-8.7501402451
     $        09242D5*xyzi(i+5,j+7,k+1)+6.805664635084966D5*xyzi(i+5,j+5
     $        ,k+1)-4.004956867237501D5*xyzi(i+3,j+11,k+1)+7.29178353759
     $        1036D5*xyzi(i+3,j+9,k+1)-3.240792683373793D5*xyzi(i+3,j+7,
     $        k+1)+4.004956867237501D4*xyzi(i+1,j+13,k+1)-6.628894125082
     $        759D4*xyzi(i+1,j+11,k+1)+2.700660569478161D4*xyzi(i+1,j+9,
     $        k+1))+angi
            angi=yrk(233)*(-1.210985129400086D5*xyzi(i+13,j+1,k+1)+4.843
     $        940517600344D5*xyzi(i+11,j+3,k+1)+2.756035122082954D5*xyzi
     $        (i+11,j+1,k+1)+1.332083642340095D6*xyzi(i+9,j+5,k+1)-1.378
     $        017561041477D6*xyzi(i+9,j+3,k+1)-2.041507497839225D5*xyzi(
     $        i+9,j+1,k+1)-1.653621073249773D6*xyzi(i+7,j+5,k+1)+1.22490
     $        4498703535D6*xyzi(i+7,j+3,k+1)+4.899617994814141D4*xyzi(i+
     $        7,j+1,k+1)-1.332083642340095D6*xyzi(i+5,j+9,k+1)+1.6536210
     $        73249773D6*xyzi(i+5,j+7,k+1)-3.429732596369899D5*xyzi(i+5,
     $        j+3,k+1)-4.843940517600344D5*xyzi(i+3,j+11,k+1)+1.37801756
     $        1041477D6*xyzi(i+3,j+9,k+1)-1.224904498703535D6*xyzi(i+3,j
     $        +7,k+1)+3.429732596369899D5*xyzi(i+3,j+5,k+1)+1.2109851294
     $        00086D5*xyzi(i+1,j+13,k+1)-2.756035122082954D5*xyzi(i+1,j+
     $        11,k+1)+2.041507497839225D5*xyzi(i+1,j+9,k+1)-4.8996179948
     $        14141D4*xyzi(i+1,j+7,k+1))+angi
            angi=yrk(235)*(2.407737205532407D5*xyzi(i+13,j+1,k+1)+1.6051
     $        58137021604D5*xyzi(i+11,j+3,k+1)-6.642033670434226D5*xyzi(
     $        i+11,j+1,k+1)-1.524900230170524D6*xyzi(i+9,j+5,k+1)+2.2140
     $        11223478075D5*xyzi(i+9,j+3,k+1)+6.642033670434226D5*xyzi(i
     $        +9,j+1,k+1)-2.889284646638888D6*xyzi(i+7,j+7,k+1)+3.985220
     $        202260535D6*xyzi(i+7,j+5,k+1)-8.856044893912301D5*xyzi(i+7
     $        ,j+3,k+1)-2.833934366051936D5*xyzi(i+7,j+1,k+1)-1.52490023
     $        0170524D6*xyzi(i+5,j+9,k+1)+3.985220202260535D6*xyzi(i+5,j
     $        +7,k+1)-3.099615712869305D6*xyzi(i+5,j+5,k+1)+6.6125135207
     $        87851D5*xyzi(i+5,j+3,k+1)+4.31250881790512D4*xyzi(i+5,j+1,
     $        k+1)+1.605158137021604D5*xyzi(i+3,j+11,k+1)+2.214011223478
     $        075D5*xyzi(i+3,j+9,k+1)-8.856044893912301D5*xyzi(i+3,j+7,k
     $        +1)+6.612513520787851D5*xyzi(i+3,j+5,k+1)-1.43750293930170
     $        7D5*xyzi(i+3,j+3,k+1)+2.407737205532407D5*xyzi(i+1,j+13,k+
     $        1)-6.642033670434226D5*xyzi(i+1,j+11,k+1)+6.64203367043422
     $        6D5*xyzi(i+1,j+9,k+1)-2.833934366051936D5*xyzi(i+1,j+7,k+1
     $        )+4.31250881790512D4*xyzi(i+1,j+5,k+1))+angi
            angi=yrk(237)*(-3.136506033551579D5*xyzi(i+13,j+1,k+1)-1.254
     $        602413420632D6*xyzi(i+11,j+3,k+1)+9.733984242056625D5*xyzi
     $        (i+11,j+1,k+1)-1.56825301677579D6*xyzi(i+9,j+5,k+1)+2.9201
     $        95272616987D6*xyzi(i+9,j+3,k+1)-1.153657391651156D6*xyzi(i
     $        +9,j+1,k+1)+1.946796848411325D6*xyzi(i+7,j+5,k+1)-2.307314
     $        783302311D6*xyzi(i+7,j+3,k+1)+6.460481393246471D5*xyzi(i+7
     $        ,j+1,k+1)+1.56825301677579D6*xyzi(i+5,j+9,k+1)-1.946796848
     $        411325D6*xyzi(i+5,j+7,k+1)+6.460481393246471D5*xyzi(i+5,j+
     $        3,k+1)-1.685342972151253D5*xyzi(i+5,j+1,k+1)+1.25460241342
     $        0632D6*xyzi(i+3,j+11,k+1)-2.920195272616987D6*xyzi(i+3,j+9
     $        ,k+1)+2.307314783302311D6*xyzi(i+3,j+7,k+1)-6.460481393246
     $        471D5*xyzi(i+3,j+5,k+1)+1.605088544905956D4*xyzi(i+3,j+1,k
     $        +1)+3.136506033551579D5*xyzi(i+1,j+13,k+1)-9.7339842420566
     $        25D5*xyzi(i+1,j+11,k+1)+1.153657391651156D6*xyzi(i+1,j+9,k
     $        +1)-6.460481393246471D5*xyzi(i+1,j+7,k+1)+1.68534297215125
     $        3D5*xyzi(i+1,j+5,k+1)-1.605088544905956D4*xyzi(i+1,j+3,k+1
     $        ))+angi
            angi=yrk(239)*(2.322024969253071D5*xyzi(i+13,j+1,k+1)+1.3932
     $        14981551842D6*xyzi(i+11,j+3,k+1)-7.686703346492923D5*xyzi(
     $        i+11,j+1,k+1)+3.483037453879606D6*xyzi(i+9,j+5,k+1)-3.8433
     $        51673246461D6*xyzi(i+9,j+3,k+1)+9.964245078787123D5*xyzi(i
     $        +9,j+1,k+1)+4.644049938506141D6*xyzi(i+7,j+7,k+1)-7.686703
     $        346492923D6*xyzi(i+7,j+5,k+1)+3.985698031514849D6*xyzi(i+7
     $        ,j+3,k+1)-6.377116850423759D5*xyzi(i+7,j+1,k+1)+3.48303745
     $        3879606D6*xyzi(i+5,j+9,k+1)-7.686703346492923D6*xyzi(i+5,j
     $        +7,k+1)+5.978547047272274D6*xyzi(i+5,j+5,k+1)-1.9131350551
     $        27128D6*xyzi(i+5,j+3,k+1)+2.079494625138182D5*xyzi(i+5,j+1
     $        ,k+1)+1.393214981551842D6*xyzi(i+3,j+11,k+1)-3.84335167324
     $        6461D6*xyzi(i+3,j+9,k+1)+3.985698031514849D6*xyzi(i+3,j+7,
     $        k+1)-1.913135055127128D6*xyzi(i+3,j+5,k+1)+4.1589892502763
     $        64D5*xyzi(i+3,j+3,k+1)-3.168753714496277D4*xyzi(i+3,j+1,k+
     $        1)+2.322024969253071D5*xyzi(i+1,j+13,k+1)-7.68670334649292
     $        3D5*xyzi(i+1,j+11,k+1)+9.964245078787123D5*xyzi(i+1,j+9,k+
     $        1)-6.377116850423759D5*xyzi(i+1,j+7,k+1)+2.079494625138182
     $        D5*xyzi(i+1,j+5,k+1)-3.168753714496277D4*xyzi(i+1,j+3,k+1)
     $        +1.667765112892778D3*xyzi(i+1,j+1,k+1))+angi
          endif
          OMEGA(15)=angi
        endif
      endif
c...
      end
c=======================================================================
      subroutine omega2one(i,j,k,lp,ii,lpmax,ikang,yrk,mijk,xyzi,o1,OMEG
     $  A)
      implicit none
c...
c...  Max2fort (MM) version 0.1 (September 2004)
c...  translated on 12/9/2004 17:52:23
c...
c...  This subroutine computes type 2 angular integrals over pseudopotentials
c...
c...
c...   omega(i,j,k,l,mp,lp) =
c...
c...   l
c...  ====              / /
c...  \                 [ [                             i   j   k
c...   >  yr(m,l,tk,pk) I I yr(m,l,T,P) yr(mp,lp,T,P) xs  ys  zs  dP dT
c...  /                 ] ]
c...  ====              / /
c...  m = - l
c...
c...  where yr are real spherical harmonic functions, and xs, ys and zs
c...  are cartesian coordinates costrained over the surface of the unit
c...  sphere:
c...              xs=sin(T) cos(P)
c...              ys=sin(T) sin(P)
c...              zs=cos(T)
c...
c...  The integral is computed by expanding the real sperical harmonics
c...  as polynomials in xs, ys and zs, obtaining integrals of the type:
c...
c...                      / /
c...                 1    [ [   a   b   c         (a-1)!! (b-1)!! (c-1)!!
c...  xyzi(a,b,c)= -----  I I xs  ys  zs  dP dT = -----------------------
c...               4 %pi  ] ]                          (a+b+c+1)!!
c...                      / /
c...
c...  if a, b and c are all even, or zero otherwise
c...
c...  This subroutine computes a batch of integrals
c...  with l=ii,  and mp=-lp,lp
c...  with the additional constraint that l-(i+j+k) must be even
c...
c...  **** This version works for lp up to 6 ****
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (24),(29-30)
c...
c...  David B. Cook, Handbook of Computational Quantum Chemistry,
c...  (Oxford, 1998), page 603
c...
c...  Input parameters:
c...
c...  i       exponent of xs
c...  j       exponent of ys
c...  k       exponent of zs
c...  lp      angular momentum of the current projector
c...  ii      value of l to compute
c...  lpmax   maximum angular momentum of projectors,
c...          for dimensioning
c...  ikang   flag for polar angles of vector k:
c...          0 vector k is zero (thetak and phik are undefined)
c...            only integrals with l=0 survive
c...          1 angle phik is undefined
c...            only integrals with m=0 survive
c...          > 1 general case (both thetak and phik are defined)
c...  yrk     array of real spherical harmonics evaluated at k
c...  mijk    maximum exponent value, for dimensioning
c...  xyzi    table of double factorial products
c...  o1      scratch array for type 1 angular integrals
c...
c...  Output parameter:
c...
c...  omega   array of type two angular integrals
c...
      real*8 ZERO
      parameter (ZERO=0.0D0)
c...
      integer i,j,k,lp,lpmax,ikang,mijk,lmax,ijk,mp,ii
      real*8 yrk(*),OMEGA(-lpmax:lpmax),o1,xyzi(0:mijk,0:mijk,0:mijk)
c...
      if(lp.gt.6)then
        call nerror(1,'omega2one','maximum lp value exceeded',lp,6)
      endif
      ijk=k+j+i
      lmax=lp+ijk
      if(ii.gt.lmax)then
        call nerror(2,'omega2one','inconsistent call',ii,lmax)
      endif
c...
c...  zero out omega
c...
      do mp=-lp,lp,1
        OMEGA(mp)=ZERO
      enddo
      if(ikang.eq.0.and.ii.ne.0)return
c...
c...  lp =0
c...
      if(lp.eq.0)then
        if(ikang.eq.0)then
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(0)=2.820947917738782D-1*xyzi(i,j,k)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=2.820947917738782D-1*o1
          else
            call omega1one(i,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=2.820947917738782D-1*o1
          endif
        endif
      endif
c...
c...  lp =1
c...
      if(lp.eq.1)then
        if(ikang.eq.0)then
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=-4.886025119029199D-1*xyzi(i,j+1,k)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))then
            OMEGA(0)=4.886025119029199D-1*xyzi(i,j,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=-4.886025119029199D-1*xyzi(i+1,j,k)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=-4.886025119029199D-1*o1
            call omega0one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=4.886025119029199D-1*o1
            call omega0one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=-4.886025119029199D-1*o1
          else
            call omega1one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=-4.886025119029199D-1*o1
            call omega1one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=4.886025119029199D-1*o1
            call omega1one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=-4.886025119029199D-1*o1
          endif
        endif
      endif
c...
c...  lp =2
c...
      if(lp.eq.2)then
        if(ikang.eq.0)then
          if((mod(i+1,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=1.092548430592079D0*xyzi(i+1,j+1,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=-1.092548430592079D0*xyzi(i,j+1,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+2,2).eq.0))then
            OMEGA(0)=9.4617469575756D-1*xyzi(i,j,k+2)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(0)=OMEGA(0)-3.1539156525252D-1*xyzi(i,j,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=-1.092548430592079D0*xyzi(i+1,j,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=5.462742152960396D-1*xyzi(i+2,j,k)+OMEGA(2)
          endif
          if((mod(i,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-5.462742152960396D-1*xyzi(i,j+2,k)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=1.092548430592079D0*o1
            call omega0one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=-1.092548430592079D0*o1
            call omega0one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=9.4617469575756D-1*o1+OMEGA(0)
            if(ii.le.ijk)then
              call omega0one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.1539156525252D-1*o1
            endif
            call omega0one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=-1.092548430592079D0*o1
            call omega0one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=5.462742152960396D-1*o1+OMEGA(2)
            call omega0one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-5.462742152960396D-1*o1
          else
            call omega1one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=1.092548430592079D0*o1
            call omega1one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=-1.092548430592079D0*o1
            call omega1one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=9.4617469575756D-1*o1+OMEGA(0)
            if(ii.le.ijk)then
              call omega1one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.1539156525252D-1*o1
            endif
            call omega1one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=-1.092548430592079D0*o1
            call omega1one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=5.462742152960396D-1*o1+OMEGA(2)
            call omega1one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-5.462742152960396D-1*o1
          endif
        endif
      endif
c...
c...  lp =3
c...
      if(lp.eq.3)then
        if(ikang.eq.0)then
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-3)=5.900435899266435D-1*xyzi(i,j+3,k)+OMEGA(-3)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-3)=OMEGA(-3)-1.770130769779931D0*xyzi(i+2,j+1,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-2)=2.890611442640554D0*xyzi(i+1,j+1,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=2.285228997322329D0*xyzi(i,j+3,k)+OMEGA(-1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-1)=2.285228997322329D0*xyzi(i+2,j+1,k)+OMEGA(-1)
          endif
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=OMEGA(-1)-1.828183197857863D0*xyzi(i,j+1,k)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+3,2).eq.0))then
            OMEGA(0)=1.865881662950577D0*xyzi(i,j,k+3)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))then
            OMEGA(0)=OMEGA(0)-1.119528997770346D0*xyzi(i,j,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(1)=2.285228997322329D0*xyzi(i+1,j+2,k)+OMEGA(1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=2.285228997322329D0*xyzi(i+3,j,k)+OMEGA(1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=OMEGA(1)-1.828183197857863D0*xyzi(i+1,j,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=1.445305721320277D0*xyzi(i+2,j,k+1)+OMEGA(2)
          endif
          if((mod(i,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=OMEGA(2)-1.445305721320277D0*xyzi(i,j+2,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(3)=1.770130769779931D0*xyzi(i+1,j+2,k)+OMEGA(3)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(3)=OMEGA(3)-5.900435899266435D-1*xyzi(i+3,j,k)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=5.900435899266435D-1*o1+OMEGA(-3)
            call omega0one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-1.770130769779931D0*o1
            call omega0one(i+1,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=2.890611442640554D0*o1
            call omega0one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=2.285228997322329D0*o1+OMEGA(-1)
            call omega0one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=2.285228997322329D0*o1+OMEGA(-1)
            if(ii.le.ijk+1)then
              call omega0one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-1.828183197857863D0*o1
            endif
            call omega0one(i,j,k+3,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=1.865881662950577D0*o1+OMEGA(0)
            if(ii.le.ijk+1)then
              call omega0one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-1.119528997770346D0*o1
            endif
            call omega0one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=2.285228997322329D0*o1+OMEGA(1)
            call omega0one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=2.285228997322329D0*o1+OMEGA(1)
            if(ii.le.ijk+1)then
              call omega0one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-1.828183197857863D0*o1
            endif
            call omega0one(i+2,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.445305721320277D0*o1+OMEGA(2)
            call omega0one(i,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.445305721320277D0*o1
            call omega0one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=1.770130769779931D0*o1+OMEGA(3)
            call omega0one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-5.900435899266435D-1*o1
          else
            call omega1one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=5.900435899266435D-1*o1+OMEGA(-3)
            call omega1one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-1.770130769779931D0*o1
            call omega1one(i+1,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=2.890611442640554D0*o1
            call omega1one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=2.285228997322329D0*o1+OMEGA(-1)
            call omega1one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=2.285228997322329D0*o1+OMEGA(-1)
            if(ii.le.ijk+1)then
              call omega1one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-1.828183197857863D0*o1
            endif
            call omega1one(i,j,k+3,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=1.865881662950577D0*o1+OMEGA(0)
            if(ii.le.ijk+1)then
              call omega1one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-1.119528997770346D0*o1
            endif
            call omega1one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=2.285228997322329D0*o1+OMEGA(1)
            call omega1one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=2.285228997322329D0*o1+OMEGA(1)
            if(ii.le.ijk+1)then
              call omega1one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-1.828183197857863D0*o1
            endif
            call omega1one(i+2,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.445305721320277D0*o1+OMEGA(2)
            call omega1one(i,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.445305721320277D0*o1
            call omega1one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=1.770130769779931D0*o1+OMEGA(3)
            call omega1one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-5.900435899266435D-1*o1
          endif
        endif
      endif
c...
c...  lp =4
c...
      if(lp.eq.4)then
        if(ikang.eq.0)then
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=2.503342941796705D0*xyzi(i+3,j+1,k)+OMEGA(-4)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=OMEGA(-4)-2.503342941796705D0*xyzi(i+1,j+3,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-3)=1.770130769779931D0*xyzi(i,j+3,k+1)+OMEGA(-3)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-3)=OMEGA(-3)-5.310392309339792D0*xyzi(i+2,j+1,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*xyzi(i+1,j+3,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*xyzi(i+3,j+1,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=5.67704817454536D0*xyzi(i+1,j+1,k)+OMEGA(-2)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=4.683325804901024D0*xyzi(i,j+3,k+1)+OMEGA(-1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-1)=4.683325804901024D0*xyzi(i+2,j+1,k+1)+OMEGA(-1)
          endif
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=OMEGA(-1)-2.676186174229157D0*xyzi(i,j+1,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+4,2).eq.0))then
            OMEGA(0)=3.702494142032151D0*xyzi(i,j,k+4)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+2,2).eq.0))then
            OMEGA(0)=OMEGA(0)-3.173566407456129D0*xyzi(i,j,k+2)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(0)=3.173566407456129D-1*xyzi(i,j,k)+OMEGA(0)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(1)=4.683325804901024D0*xyzi(i+1,j+2,k+1)+OMEGA(1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=4.683325804901024D0*xyzi(i+3,j,k+1)+OMEGA(1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=OMEGA(1)-2.676186174229157D0*xyzi(i+1,j,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=3.31161143515146D0*xyzi(i,j+4,k)+OMEGA(2)
          endif
          if((mod(i,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-2.83852408727268D0*xyzi(i,j+2,k)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-3.31161143515146D0*xyzi(i+4,j,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=2.83852408727268D0*xyzi(i+2,j,k)+OMEGA(2)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(3)=5.310392309339792D0*xyzi(i+1,j+2,k+1)+OMEGA(3)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(3)=OMEGA(3)-1.770130769779931D0*xyzi(i+3,j,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=6.258357354491762D-1*xyzi(i,j+4,k)+OMEGA(4)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(4)=OMEGA(4)-3.755014412695057D0*xyzi(i+2,j+2,k)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=6.258357354491762D-1*xyzi(i+4,j,k)+OMEGA(4)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=2.503342941796705D0*o1+OMEGA(-4)
            call omega0one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-2.503342941796705D0*o1
            call omega0one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=1.770130769779931D0*o1+OMEGA(-3)
            call omega0one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-5.310392309339792D0*o1
            call omega0one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*o1
            call omega0one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*o1
            if(ii.le.ijk+2)then
              call omega0one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=5.67704817454536D0*o1+OMEGA(-2)
            endif
            call omega0one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=4.683325804901024D0*o1+OMEGA(-1)
            call omega0one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=4.683325804901024D0*o1+OMEGA(-1)
            if(ii.le.ijk+2)then
              call omega0one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-2.676186174229157D0*o1
            endif
            call omega0one(i,j,k+4,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=3.702494142032151D0*o1+OMEGA(0)
            if(ii.le.ijk+2)then
              call omega0one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.173566407456129D0*o1
            endif
            if(ii.le.ijk)then
              call omega0one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=3.173566407456129D-1*o1+OMEGA(0)
            endif
            call omega0one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=4.683325804901024D0*o1+OMEGA(1)
            call omega0one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=4.683325804901024D0*o1+OMEGA(1)
            if(ii.le.ijk+2)then
              call omega0one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-2.676186174229157D0*o1
            endif
            call omega0one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=3.31161143515146D0*o1+OMEGA(2)
            if(ii.le.ijk+2)then
              call omega0one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-2.83852408727268D0*o1
            endif
            call omega0one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-3.31161143515146D0*o1
            if(ii.le.ijk+2)then
              call omega0one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=2.83852408727268D0*o1+OMEGA(2)
            endif
            call omega0one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=5.310392309339792D0*o1+OMEGA(3)
            call omega0one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-1.770130769779931D0*o1
            call omega0one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=6.258357354491762D-1*o1+OMEGA(4)
            call omega0one(i+2,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-3.755014412695057D0*o1
            call omega0one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=6.258357354491762D-1*o1+OMEGA(4)
          else
            call omega1one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=2.503342941796705D0*o1+OMEGA(-4)
            call omega1one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-2.503342941796705D0*o1
            call omega1one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=1.770130769779931D0*o1+OMEGA(-3)
            call omega1one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-5.310392309339792D0*o1
            call omega1one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*o1
            call omega1one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-6.62322287030292D0*o1
            if(ii.le.ijk+2)then
              call omega1one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=5.67704817454536D0*o1+OMEGA(-2)
            endif
            call omega1one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=4.683325804901024D0*o1+OMEGA(-1)
            call omega1one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=4.683325804901024D0*o1+OMEGA(-1)
            if(ii.le.ijk+2)then
              call omega1one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-2.676186174229157D0*o1
            endif
            call omega1one(i,j,k+4,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=3.702494142032151D0*o1+OMEGA(0)
            if(ii.le.ijk+2)then
              call omega1one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.173566407456129D0*o1
            endif
            if(ii.le.ijk)then
              call omega1one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=3.173566407456129D-1*o1+OMEGA(0)
            endif
            call omega1one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=4.683325804901024D0*o1+OMEGA(1)
            call omega1one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=4.683325804901024D0*o1+OMEGA(1)
            if(ii.le.ijk+2)then
              call omega1one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-2.676186174229157D0*o1
            endif
            call omega1one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=3.31161143515146D0*o1+OMEGA(2)
            if(ii.le.ijk+2)then
              call omega1one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-2.83852408727268D0*o1
            endif
            call omega1one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-3.31161143515146D0*o1
            if(ii.le.ijk+2)then
              call omega1one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=2.83852408727268D0*o1+OMEGA(2)
            endif
            call omega1one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=5.310392309339792D0*o1+OMEGA(3)
            call omega1one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-1.770130769779931D0*o1
            call omega1one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=6.258357354491762D-1*o1+OMEGA(4)
            call omega1one(i+2,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-3.755014412695057D0*o1
            call omega1one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=6.258357354491762D-1*o1+OMEGA(4)
          endif
        endif
      endif
c...
c...  lp =5
c...
      if(lp.eq.5)then
        if(ikang.eq.0)then
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-5)=OMEGA(-5)-6.563820568401701D-1*xyzi(i,j+5,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-5)=6.563820568401701D0*xyzi(i+2,j+3,k)+OMEGA(-5)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-5)=OMEGA(-5)-3.281910284200851D0*xyzi(i+4,j+1,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-4)=8.302649259524165D0*xyzi(i+3,j+1,k+1)+OMEGA(-4)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-4)=OMEGA(-4)-8.302649259524165D0*xyzi(i+1,j+3,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-3)=OMEGA(-3)-4.403144694917254D0*xyzi(i,j+5,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-3)=8.806289389834507D0*xyzi(i+2,j+3,k)+OMEGA(-3)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-3)=3.913906395482003D0*xyzi(i,j+3,k)+OMEGA(-3)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-3)=1.320943408475176D1*xyzi(i+4,j+1,k)+OMEGA(-3)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-3)=OMEGA(-3)-1.174171918644601D1*xyzi(i+2,j+1,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*xyzi(i+1,j+3,k+1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*xyzi(i+3,j+1,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-2)=9.587073569946648D0*xyzi(i+1,j+1,k+1)+OMEGA(-2)
          endif
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*xyzi(i,j+5,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-1)=OMEGA(-1)-1.902375935021927D1*xyzi(i+2,j+3,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=1.268250623347951D1*xyzi(i,j+3,k)+OMEGA(-1)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*xyzi(i+4,j+1,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-1)=1.268250623347951D1*xyzi(i+2,j+1,k)+OMEGA(-1)
          endif
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(-1)=OMEGA(-1)-3.623573209565575D0*xyzi(i,j+1,k)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+5,2).eq.0))then
            OMEGA(0)=7.367870314565687D0*xyzi(i,j,k+5)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+3,2).eq.0))then
            OMEGA(0)=OMEGA(0)-8.186522571739652D0*xyzi(i,j,k+3)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))then
            OMEGA(0)=1.754254836801354D0*xyzi(i,j,k+1)+OMEGA(0)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*xyzi(i+1,j+4,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(1)=OMEGA(1)-1.902375935021927D1*xyzi(i+3,j+2,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(1)=1.268250623347951D1*xyzi(i+1,j+2,k)+OMEGA(1)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*xyzi(i+5,j,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=1.268250623347951D1*xyzi(i+3,j,k)+OMEGA(1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(1)=OMEGA(1)-3.623573209565575D0*xyzi(i+1,j,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=7.190305177459986D0*xyzi(i,j+4,k+1)+OMEGA(2)
          endif
          if((mod(i,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=OMEGA(2)-4.793536784973324D0*xyzi(i,j+2,k+1)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=OMEGA(2)-7.190305177459986D0*xyzi(i+4,j,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(2)=4.793536784973324D0*xyzi(i+2,j,k+1)+OMEGA(2)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(3)=OMEGA(3)-1.320943408475176D1*xyzi(i+1,j+4,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(3)=OMEGA(3)-8.806289389834507D0*xyzi(i+3,j+2,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(3)=1.174171918644601D1*xyzi(i+1,j+2,k)+OMEGA(3)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(3)=4.403144694917254D0*xyzi(i+5,j,k)+OMEGA(3)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(3)=OMEGA(3)-3.913906395482003D0*xyzi(i+3,j,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(4)=2.075662314881041D0*xyzi(i,j+4,k+1)+OMEGA(4)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(4)=OMEGA(4)-1.245397388928625D1*xyzi(i+2,j+2,k+1)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(4)=2.075662314881041D0*xyzi(i+4,j,k+1)+OMEGA(4)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(5)=OMEGA(5)-3.281910284200851D0*xyzi(i+1,j+4,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(5)=6.563820568401701D0*xyzi(i+3,j+2,k)+OMEGA(5)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(5)=OMEGA(5)-6.563820568401701D-1*xyzi(i+5,j,k)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-6.563820568401701D-1*o1
            call omega0one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=6.563820568401701D0*o1+OMEGA(-5)
            call omega0one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-3.281910284200851D0*o1
            call omega0one(i+3,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=8.302649259524165D0*o1+OMEGA(-4)
            call omega0one(i+1,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-8.302649259524165D0*o1
            call omega0one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-4.403144694917254D0*o1
            call omega0one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=8.806289389834507D0*o1+OMEGA(-3)
            if(ii.le.ijk+3)then
              call omega0one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=3.913906395482003D0*o1+OMEGA(-3)
            endif
            call omega0one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=1.320943408475176D1*o1+OMEGA(-3)
            if(ii.le.ijk+3)then
              call omega0one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=OMEGA(-3)-1.174171918644601D1*o1
            endif
            call omega0one(i+1,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*o1
            call omega0one(i+3,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*o1
            if(ii.le.ijk+3)then
              call omega0one(i+1,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=9.587073569946648D0*o1+OMEGA(-2)
            endif
            call omega0one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*o1
            call omega0one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.902375935021927D1*o1
            if(ii.le.ijk+3)then
              call omega0one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=1.268250623347951D1*o1+OMEGA(-1)
            endif
            call omega0one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*o1
            if(ii.le.ijk+3)then
              call omega0one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=1.268250623347951D1*o1+OMEGA(-1)
            endif
            if(ii.le.ijk+1)then
              call omega0one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-3.623573209565575D0*o1
            endif
            call omega0one(i,j,k+5,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=7.367870314565687D0*o1+OMEGA(0)
            if(ii.le.ijk+3)then
              call omega0one(i,j,k+3,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-8.186522571739652D0*o1
            endif
            if(ii.le.ijk+1)then
              call omega0one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=1.754254836801354D0*o1+OMEGA(0)
            endif
            call omega0one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*o1
            call omega0one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.902375935021927D1*o1
            if(ii.le.ijk+3)then
              call omega0one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=1.268250623347951D1*o1+OMEGA(1)
            endif
            call omega0one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*o1
            if(ii.le.ijk+3)then
              call omega0one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=1.268250623347951D1*o1+OMEGA(1)
            endif
            if(ii.le.ijk+1)then
              call omega0one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-3.623573209565575D0*o1
            endif
            call omega0one(i,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=7.190305177459986D0*o1+OMEGA(2)
            if(ii.le.ijk+3)then
              call omega0one(i,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-4.793536784973324D0*o1
            endif
            call omega0one(i+4,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-7.190305177459986D0*o1
            if(ii.le.ijk+3)then
              call omega0one(i+2,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=4.793536784973324D0*o1+OMEGA(2)
            endif
            call omega0one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-1.320943408475176D1*o1
            call omega0one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-8.806289389834507D0*o1
            if(ii.le.ijk+3)then
              call omega0one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=1.174171918644601D1*o1+OMEGA(3)
            endif
            call omega0one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=4.403144694917254D0*o1+OMEGA(3)
            if(ii.le.ijk+3)then
              call omega0one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=OMEGA(3)-3.913906395482003D0*o1
            endif
            call omega0one(i,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.075662314881041D0*o1+OMEGA(4)
            call omega0one(i+2,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-1.245397388928625D1*o1
            call omega0one(i+4,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.075662314881041D0*o1+OMEGA(4)
            call omega0one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-3.281910284200851D0*o1
            call omega0one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=6.563820568401701D0*o1+OMEGA(5)
            call omega0one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-6.563820568401701D-1*o1
          else
            call omega1one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-6.563820568401701D-1*o1
            call omega1one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=6.563820568401701D0*o1+OMEGA(-5)
            call omega1one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-3.281910284200851D0*o1
            call omega1one(i+3,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=8.302649259524165D0*o1+OMEGA(-4)
            call omega1one(i+1,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-8.302649259524165D0*o1
            call omega1one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-4.403144694917254D0*o1
            call omega1one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=8.806289389834507D0*o1+OMEGA(-3)
            if(ii.le.ijk+3)then
              call omega1one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=3.913906395482003D0*o1+OMEGA(-3)
            endif
            call omega1one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=1.320943408475176D1*o1+OMEGA(-3)
            if(ii.le.ijk+3)then
              call omega1one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=OMEGA(-3)-1.174171918644601D1*o1
            endif
            call omega1one(i+1,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*o1
            call omega1one(i+3,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=OMEGA(-2)-1.438061035491997D1*o1
            if(ii.le.ijk+3)then
              call omega1one(i+1,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=9.587073569946648D0*o1+OMEGA(-2)
            endif
            call omega1one(i,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*o1
            call omega1one(i+2,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.902375935021927D1*o1
            if(ii.le.ijk+3)then
              call omega1one(i,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=1.268250623347951D1*o1+OMEGA(-1)
            endif
            call omega1one(i+4,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-9.511879675109636D0*o1
            if(ii.le.ijk+3)then
              call omega1one(i+2,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=1.268250623347951D1*o1+OMEGA(-1)
            endif
            if(ii.le.ijk+1)then
              call omega1one(i,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-3.623573209565575D0*o1
            endif
            call omega1one(i,j,k+5,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=7.367870314565687D0*o1+OMEGA(0)
            if(ii.le.ijk+3)then
              call omega1one(i,j,k+3,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-8.186522571739652D0*o1
            endif
            if(ii.le.ijk+1)then
              call omega1one(i,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=1.754254836801354D0*o1+OMEGA(0)
            endif
            call omega1one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*o1
            call omega1one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.902375935021927D1*o1
            if(ii.le.ijk+3)then
              call omega1one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=1.268250623347951D1*o1+OMEGA(1)
            endif
            call omega1one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-9.511879675109636D0*o1
            if(ii.le.ijk+3)then
              call omega1one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=1.268250623347951D1*o1+OMEGA(1)
            endif
            if(ii.le.ijk+1)then
              call omega1one(i+1,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-3.623573209565575D0*o1
            endif
            call omega1one(i,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=7.190305177459986D0*o1+OMEGA(2)
            if(ii.le.ijk+3)then
              call omega1one(i,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-4.793536784973324D0*o1
            endif
            call omega1one(i+4,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-7.190305177459986D0*o1
            if(ii.le.ijk+3)then
              call omega1one(i+2,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=4.793536784973324D0*o1+OMEGA(2)
            endif
            call omega1one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-1.320943408475176D1*o1
            call omega1one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-8.806289389834507D0*o1
            if(ii.le.ijk+3)then
              call omega1one(i+1,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=1.174171918644601D1*o1+OMEGA(3)
            endif
            call omega1one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=4.403144694917254D0*o1+OMEGA(3)
            if(ii.le.ijk+3)then
              call omega1one(i+3,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=OMEGA(3)-3.913906395482003D0*o1
            endif
            call omega1one(i,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.075662314881041D0*o1+OMEGA(4)
            call omega1one(i+2,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-1.245397388928625D1*o1
            call omega1one(i+4,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.075662314881041D0*o1+OMEGA(4)
            call omega1one(i+1,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-3.281910284200851D0*o1
            call omega1one(i+3,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=6.563820568401701D0*o1+OMEGA(5)
            call omega1one(i+5,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-6.563820568401701D-1*o1
          endif
        endif
      endif
c...
c...  lp =6
c...
      if(lp.eq.6)then
        if(ikang.eq.0)then
          if((mod(i+1,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-6)=4.099104631151486D0*xyzi(i+1,j+5,k)+OMEGA(-6)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-6)=OMEGA(-6)-1.366368210383829D1*xyzi(i+3,j+3,k)
          endif
          if((mod(i+5,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-6)=4.099104631151486D0*xyzi(i+5,j+1,k)+OMEGA(-6)
          endif
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-5)=OMEGA(-5)-2.366619162231752D0*xyzi(i,j+5,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-5)=2.366619162231752D1*xyzi(i+2,j+3,k+1)+OMEGA(-5)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-5)=OMEGA(-5)-1.183309581115876D1*xyzi(i+4,j+1,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=2.220085563206386D1*xyzi(i+1,j+5,k)+OMEGA(-4)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=OMEGA(-4)-2.018259602914897D1*xyzi(i+1,j+3,k)
          endif
          if((mod(i+5,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=OMEGA(-4)-2.220085563206386D1*xyzi(i+5,j+1,k)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-4)=2.018259602914897D1*xyzi(i+3,j+1,k)+OMEGA(-4)
          endif
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-3)=OMEGA(-3)-1.013325785466416D1*xyzi(i,j+5,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-3)=2.026651570932832D1*xyzi(i+2,j+3,k+1)+OMEGA(-3)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-3)=7.369642076119388D0*xyzi(i,j+3,k+1)+OMEGA(-3)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-3)=3.039977356399248D1*xyzi(i+4,j+1,k+1)+OMEGA(-3)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-3)=OMEGA(-3)-2.210892622835816D1*xyzi(i+2,j+1,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=3.039977356399248D1*xyzi(i+1,j+5,k)+OMEGA(-2)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=6.079954712798495D1*xyzi(i+3,j+3,k)+OMEGA(-2)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*xyzi(i+1,j+3,k)
          endif
          if((mod(i+5,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=3.039977356399248D1*xyzi(i+5,j+1,k)+OMEGA(-2)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*xyzi(i+3,j+1,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(-2)=1.473928415223878D1*xyzi(i+1,j+1,k)+OMEGA(-2)
          endif
          if((mod(i,2).eq.0.and.mod(j+5,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*xyzi(i,j+5,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-1)=OMEGA(-1)-3.845300992623627D1*xyzi(i+2,j+3,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+3,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=2.097436905067433D1*xyzi(i,j+3,k+1)+OMEGA(-1)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*xyzi(i+4,j+1,k+1)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(-1)=2.097436905067433D1*xyzi(i+2,j+1,k+1)+OMEGA(-1)
          endif
          if((mod(i,2).eq.0.and.mod(j+1,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(-1)=OMEGA(-1)-4.660970900149851D0*xyzi(i,j+1,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+6,2).eq.0))then
            OMEGA(0)=1.468448572382217D1*xyzi(i,j,k+6)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+4,2).eq.0))then
            OMEGA(0)=OMEGA(0)-2.002429871430295D1*xyzi(i,j,k+4)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k+2,2).eq.0))then
            OMEGA(0)=6.674766238100985D0*xyzi(i,j,k+2)+OMEGA(0)
          endif
          if((mod(i,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(0)=OMEGA(0)-3.178460113381421D-1*xyzi(i,j,k)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*xyzi(i+1,j+4,k+1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(1)=OMEGA(1)-3.845300992623627D1*xyzi(i+3,j+2,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(1)=2.097436905067433D1*xyzi(i+1,j+2,k+1)+OMEGA(1)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*xyzi(i+5,j,k+1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=2.097436905067433D1*xyzi(i+3,j,k+1)+OMEGA(1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(1)=OMEGA(1)-4.660970900149851D0*xyzi(i+1,j,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+6,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*xyzi(i,j+6,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*xyzi(i+2,j+4,k)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=2.210892622835816D1*xyzi(i,j+4,k)+OMEGA(2)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(2)=1.519988678199624D1*xyzi(i+4,j+2,k)+OMEGA(2)
          endif
          if((mod(i,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-7.369642076119388D0*xyzi(i,j+2,k)
          endif
          if((mod(i+6,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=1.519988678199624D1*xyzi(i+6,j,k)+OMEGA(2)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=OMEGA(2)-2.210892622835816D1*xyzi(i+4,j,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(2)=7.369642076119388D0*xyzi(i+2,j,k)+OMEGA(2)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(3)=OMEGA(3)-3.039977356399248D1*xyzi(i+1,j+4,k+1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(3)=OMEGA(3)-2.026651570932832D1*xyzi(i+3,j+2,k+1)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(3)=2.210892622835816D1*xyzi(i+1,j+2,k+1)+OMEGA(3)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(3)=1.013325785466416D1*xyzi(i+5,j,k+1)+OMEGA(3)
          endif
          if((mod(i+3,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(3)=OMEGA(3)-7.369642076119388D0*xyzi(i+3,j,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+6,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*xyzi(i,j+6,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(4)=2.775106954007983D1*xyzi(i+2,j+4,k)+OMEGA(4)
          endif
          if((mod(i,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=5.045649007287242D0*xyzi(i,j+4,k)+OMEGA(4)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(4)=2.775106954007983D1*xyzi(i+4,j+2,k)+OMEGA(4)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(4)=OMEGA(4)-3.027389404372345D1*xyzi(i+2,j+2,k)
          endif
          if((mod(i+6,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*xyzi(i+6,j,k)
          endif
          if((mod(i+4,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(4)=5.045649007287242D0*xyzi(i+4,j,k)+OMEGA(4)
          endif
          if((mod(i+1,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(5)=OMEGA(5)-1.183309581115876D1*xyzi(i+1,j+4,k+1)
          endif
          if((mod(i+3,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k+1,2).eq.0))t
     $      hen
            OMEGA(5)=2.366619162231752D1*xyzi(i+3,j+2,k+1)+OMEGA(5)
          endif
          if((mod(i+5,2).eq.0.and.mod(j,2).eq.0.and.mod(k+1,2).eq.0))the
     $      n
            OMEGA(5)=OMEGA(5)-2.366619162231752D0*xyzi(i+5,j,k+1)
          endif
          if((mod(i,2).eq.0.and.mod(j+6,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(6)=OMEGA(6)-6.831841051919143D-1*xyzi(i,j+6,k)
          endif
          if((mod(i+2,2).eq.0.and.mod(j+4,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(6)=1.024776157787872D1*xyzi(i+2,j+4,k)+OMEGA(6)
          endif
          if((mod(i+4,2).eq.0.and.mod(j+2,2).eq.0.and.mod(k,2).eq.0))the
     $      n
            OMEGA(6)=OMEGA(6)-1.024776157787872D1*xyzi(i+4,j+2,k)
          endif
          if((mod(i+6,2).eq.0.and.mod(j,2).eq.0.and.mod(k,2).eq.0))then
            OMEGA(6)=6.831841051919143D-1*xyzi(i+6,j,k)+OMEGA(6)
          endif
        else
          if(ikang.eq.1)then
            call omega0one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=4.099104631151486D0*o1+OMEGA(-6)
            call omega0one(i+3,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=OMEGA(-6)-1.366368210383829D1*o1
            call omega0one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=4.099104631151486D0*o1+OMEGA(-6)
            call omega0one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-2.366619162231752D0*o1
            call omega0one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=2.366619162231752D1*o1+OMEGA(-5)
            call omega0one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-1.183309581115876D1*o1
            call omega0one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=2.220085563206386D1*o1+OMEGA(-4)
            if(ii.le.ijk+4)then
              call omega0one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-4)=OMEGA(-4)-2.018259602914897D1*o1
            endif
            call omega0one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-2.220085563206386D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-4)=2.018259602914897D1*o1+OMEGA(-4)
            endif
            call omega0one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-1.013325785466416D1*o1
            call omega0one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=2.026651570932832D1*o1+OMEGA(-3)
            if(ii.le.ijk+4)then
              call omega0one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=7.369642076119388D0*o1+OMEGA(-3)
            endif
            call omega0one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=3.039977356399248D1*o1+OMEGA(-3)
            if(ii.le.ijk+4)then
              call omega0one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=OMEGA(-3)-2.210892622835816D1*o1
            endif
            call omega0one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=3.039977356399248D1*o1+OMEGA(-2)
            call omega0one(i+3,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=6.079954712798495D1*o1+OMEGA(-2)
            if(ii.le.ijk+4)then
              call omega0one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*o1
            endif
            call omega0one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=3.039977356399248D1*o1+OMEGA(-2)
            if(ii.le.ijk+4)then
              call omega0one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega0one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=1.473928415223878D1*o1+OMEGA(-2)
            endif
            call omega0one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*o1
            call omega0one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-3.845300992623627D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=2.097436905067433D1*o1+OMEGA(-1)
            endif
            call omega0one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=2.097436905067433D1*o1+OMEGA(-1)
            endif
            if(ii.le.ijk+2)then
              call omega0one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-4.660970900149851D0*o1
            endif
            call omega0one(i,j,k+6,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=1.468448572382217D1*o1+OMEGA(0)
            if(ii.le.ijk+4)then
              call omega0one(i,j,k+4,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-2.002429871430295D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega0one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=6.674766238100985D0*o1+OMEGA(0)
            endif
            if(ii.le.ijk)then
              call omega0one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.178460113381421D-1*o1
            endif
            call omega0one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*o1
            call omega0one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-3.845300992623627D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=2.097436905067433D1*o1+OMEGA(1)
            endif
            call omega0one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=2.097436905067433D1*o1+OMEGA(1)
            endif
            if(ii.le.ijk+2)then
              call omega0one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-4.660970900149851D0*o1
            endif
            call omega0one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*o1
            call omega0one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=2.210892622835816D1*o1+OMEGA(2)
            endif
            call omega0one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.519988678199624D1*o1+OMEGA(2)
            if(ii.le.ijk+2)then
              call omega0one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-7.369642076119388D0*o1
            endif
            call omega0one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.519988678199624D1*o1+OMEGA(2)
            if(ii.le.ijk+4)then
              call omega0one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-2.210892622835816D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega0one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=7.369642076119388D0*o1+OMEGA(2)
            endif
            call omega0one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-3.039977356399248D1*o1
            call omega0one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-2.026651570932832D1*o1
            if(ii.le.ijk+4)then
              call omega0one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=2.210892622835816D1*o1+OMEGA(3)
            endif
            call omega0one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=1.013325785466416D1*o1+OMEGA(3)
            if(ii.le.ijk+4)then
              call omega0one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=OMEGA(3)-7.369642076119388D0*o1
            endif
            call omega0one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*o1
            call omega0one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.775106954007983D1*o1+OMEGA(4)
            if(ii.le.ijk+4)then
              call omega0one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=5.045649007287242D0*o1+OMEGA(4)
            endif
            call omega0one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.775106954007983D1*o1+OMEGA(4)
            if(ii.le.ijk+4)then
              call omega0one(i+2,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=OMEGA(4)-3.027389404372345D1*o1
            endif
            call omega0one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*o1
            if(ii.le.ijk+4)then
              call omega0one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=5.045649007287242D0*o1+OMEGA(4)
            endif
            call omega0one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-1.183309581115876D1*o1
            call omega0one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=2.366619162231752D1*o1+OMEGA(5)
            call omega0one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-2.366619162231752D0*o1
            call omega0one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=OMEGA(6)-6.831841051919143D-1*o1
            call omega0one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=1.024776157787872D1*o1+OMEGA(6)
            call omega0one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=OMEGA(6)-1.024776157787872D1*o1
            call omega0one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=6.831841051919143D-1*o1+OMEGA(6)
          else
            call omega1one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=4.099104631151486D0*o1+OMEGA(-6)
            call omega1one(i+3,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=OMEGA(-6)-1.366368210383829D1*o1
            call omega1one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-6)=4.099104631151486D0*o1+OMEGA(-6)
            call omega1one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-2.366619162231752D0*o1
            call omega1one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=2.366619162231752D1*o1+OMEGA(-5)
            call omega1one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-5)=OMEGA(-5)-1.183309581115876D1*o1
            call omega1one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=2.220085563206386D1*o1+OMEGA(-4)
            if(ii.le.ijk+4)then
              call omega1one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-4)=OMEGA(-4)-2.018259602914897D1*o1
            endif
            call omega1one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-4)=OMEGA(-4)-2.220085563206386D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-4)=2.018259602914897D1*o1+OMEGA(-4)
            endif
            call omega1one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=OMEGA(-3)-1.013325785466416D1*o1
            call omega1one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=2.026651570932832D1*o1+OMEGA(-3)
            if(ii.le.ijk+4)then
              call omega1one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=7.369642076119388D0*o1+OMEGA(-3)
            endif
            call omega1one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-3)=3.039977356399248D1*o1+OMEGA(-3)
            if(ii.le.ijk+4)then
              call omega1one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-3)=OMEGA(-3)-2.210892622835816D1*o1
            endif
            call omega1one(i+1,j+5,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=3.039977356399248D1*o1+OMEGA(-2)
            call omega1one(i+3,j+3,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=6.079954712798495D1*o1+OMEGA(-2)
            if(ii.le.ijk+4)then
              call omega1one(i+1,j+3,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*o1
            endif
            call omega1one(i+5,j+1,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(-2)=3.039977356399248D1*o1+OMEGA(-2)
            if(ii.le.ijk+4)then
              call omega1one(i+3,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=OMEGA(-2)-4.421785245671633D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega1one(i+1,j+1,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(-2)=1.473928415223878D1*o1+OMEGA(-2)
            endif
            call omega1one(i,j+5,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*o1
            call omega1one(i+2,j+3,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-3.845300992623627D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i,j+3,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=2.097436905067433D1*o1+OMEGA(-1)
            endif
            call omega1one(i+4,j+1,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(-1)=OMEGA(-1)-1.922650496311814D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i+2,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=2.097436905067433D1*o1+OMEGA(-1)
            endif
            if(ii.le.ijk+2)then
              call omega1one(i,j+1,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(-1)=OMEGA(-1)-4.660970900149851D0*o1
            endif
            call omega1one(i,j,k+6,ii,yrk,mijk,xyzi,o1)
            OMEGA(0)=1.468448572382217D1*o1+OMEGA(0)
            if(ii.le.ijk+4)then
              call omega1one(i,j,k+4,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-2.002429871430295D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega1one(i,j,k+2,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=6.674766238100985D0*o1+OMEGA(0)
            endif
            if(ii.le.ijk)then
              call omega1one(i,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(0)=OMEGA(0)-3.178460113381421D-1*o1
            endif
            call omega1one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*o1
            call omega1one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-3.845300992623627D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=2.097436905067433D1*o1+OMEGA(1)
            endif
            call omega1one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(1)=OMEGA(1)-1.922650496311814D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=2.097436905067433D1*o1+OMEGA(1)
            endif
            if(ii.le.ijk+2)then
              call omega1one(i+1,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(1)=OMEGA(1)-4.660970900149851D0*o1
            endif
            call omega1one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*o1
            call omega1one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=OMEGA(2)-1.519988678199624D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=2.210892622835816D1*o1+OMEGA(2)
            endif
            call omega1one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.519988678199624D1*o1+OMEGA(2)
            if(ii.le.ijk+2)then
              call omega1one(i,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-7.369642076119388D0*o1
            endif
            call omega1one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(2)=1.519988678199624D1*o1+OMEGA(2)
            if(ii.le.ijk+4)then
              call omega1one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=OMEGA(2)-2.210892622835816D1*o1
            endif
            if(ii.le.ijk+2)then
              call omega1one(i+2,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(2)=7.369642076119388D0*o1+OMEGA(2)
            endif
            call omega1one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-3.039977356399248D1*o1
            call omega1one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=OMEGA(3)-2.026651570932832D1*o1
            if(ii.le.ijk+4)then
              call omega1one(i+1,j+2,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=2.210892622835816D1*o1+OMEGA(3)
            endif
            call omega1one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(3)=1.013325785466416D1*o1+OMEGA(3)
            if(ii.le.ijk+4)then
              call omega1one(i+3,j,k+1,ii,yrk,mijk,xyzi,o1)
              OMEGA(3)=OMEGA(3)-7.369642076119388D0*o1
            endif
            call omega1one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*o1
            call omega1one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.775106954007983D1*o1+OMEGA(4)
            if(ii.le.ijk+4)then
              call omega1one(i,j+4,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=5.045649007287242D0*o1+OMEGA(4)
            endif
            call omega1one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=2.775106954007983D1*o1+OMEGA(4)
            if(ii.le.ijk+4)then
              call omega1one(i+2,j+2,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=OMEGA(4)-3.027389404372345D1*o1
            endif
            call omega1one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(4)=OMEGA(4)-5.550213908015966D0*o1
            if(ii.le.ijk+4)then
              call omega1one(i+4,j,k,ii,yrk,mijk,xyzi,o1)
              OMEGA(4)=5.045649007287242D0*o1+OMEGA(4)
            endif
            call omega1one(i+1,j+4,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-1.183309581115876D1*o1
            call omega1one(i+3,j+2,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=2.366619162231752D1*o1+OMEGA(5)
            call omega1one(i+5,j,k+1,ii,yrk,mijk,xyzi,o1)
            OMEGA(5)=OMEGA(5)-2.366619162231752D0*o1
            call omega1one(i,j+6,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=OMEGA(6)-6.831841051919143D-1*o1
            call omega1one(i+2,j+4,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=1.024776157787872D1*o1+OMEGA(6)
            call omega1one(i+4,j+2,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=OMEGA(6)-1.024776157787872D1*o1
            call omega1one(i+6,j,k,ii,yrk,mijk,xyzi,o1)
            OMEGA(6)=6.831841051919143D-1*o1+OMEGA(6)
          endif
        endif
      endif
c...
      end
c=======================================================================
      subroutine omega1one(i,j,k,lone,yrk,mijk,xyzi,OMEGA)
      implicit none
c...
c...  Max2fort (MM) version 0.1 (September 2004)
c...  translated on 12/9/2004 18:13:32
c...
c...  This subroutine computes type 1 angular integrals over pseudopotentials
c...
c...
c...   omega(i,j,k,l)=
c...
c...   l
c...  ====                    / /
c...  \                       [ [                     i   j   k
c...   >  yr(m,l,thetak,phik) I I yr(m,l,THETA,PHI) xs  ys  zs  dPHi dTHETA
c...  /                       ] ]
c...  ====                    / /
c...  m = - l
c...
c...  where yr are real spherical harmonic functions, and xs, ys and zs
c...  are cartesian coordinates costrained over the surface of the unit
c...  sphere:
c...              xs=sin(THETA) cos(PHI)
c...              ys=sin(THETA) sin(PHI)
c...              zs=cos(THETA)
c...
c...  The integral is computed by expanding the real sperical harmonic
c...  as a polynomial in xs, ys and zs, obtaining integrals of the type:
c...
c...                      / /
c...                 1    [ [   a   b   c         (a-1)!! (b-1)!! (c-1)!!
c...  xyzi(a,b,c)= -----  I I xs  ys  zs  dP dT = -----------------------
c...               4 %pi  ] ]                          (a+b+c+1)!!
c...                      / /
c...
c...  if a, b and c are all even, or zero otherwise.
c...
c...  This subroutine computes the integral with l = lone,
c...  with the additional constraint that l-(i+j+k) must be even
c...
c...  **** This version works for l up to 14 ****
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (18),(27-30)
c...
c...  David B. Cook, Handbook of Computational Quantum Chemistry,
c...  (Oxford, 1998), page 596-601
c...
c...  Input parameters:
c...
c...  i       exponent of xs
c...  j       exponent of ys
c...  k       exponent of zs
c...  lone    l value to compute
c...  yrk     array of real spherical harmonics evaluated at k
c...  mijk    maximum exponent value, for dimensioning
c...  xyzi    table of double factorial products
c...
c...  Output parameter:
c...
c...  omega   array of type one angular integrals
c...
      real*8 ZERO
      parameter (ZERO=0.0D0)
c...
      integer i,j,k,lmax,mijk,lone
      real*8 yrk(*),OMEGA,xyzi(0:mijk,0:mijk,0:mijk),angi
      logical leven,ieven,jeven,keven
c...
c...
c...  determine whether i,j and k are even or odd
c...
      ieven=mod(i,2).eq.0
      jeven=mod(j,2).eq.0
      keven=mod(k,2).eq.0
c...
c...  compute the upper limit for l
c...  and whether we must do the even or odd series
c...
      lmax=k+j+i
      if(lmax.gt.14)then
        write(6,*)'omega1one: maximum l value exceeded',lmax,14
        call nerror(1,'omega1one','maximum l value exceeded',lmax,14)
      endif
      if(lone.gt.14)then
        call nerror(2,'omega1one','inconsistent call',lone,14)
      endif
      leven=mod(lmax,2).eq.0
c...
c...  zero out omega
c...
      OMEGA=ZERO
c...
c...  chose even or odd series
c...
      if(leven)then
c...
c...  l=0
c...
        if(lone.eq.0)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=3.544907701811032D0*yrk(1)*xyzi(i,j,k)+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=2
c...
        if(lone.eq.2)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(7)*(1.188998189281803D1*xyzi(i,j,k+2)-3.96332729760
     $        6011D0*xyzi(i,j,k))+angi
            angi=yrk(9)*(6.864684246478268D0*xyzi(i+2,j,k)-6.86468424647
     $        8268D0*xyzi(i,j+2,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=1.372936849295654D1*yrk(5)*xyzi(i+1,j+1,k)+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=angi-1.372936849295654D1*yrk(8)*xyzi(i+1,j,k+1)
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=angi-1.372936849295654D1*yrk(6)*xyzi(i,j+1,k+1)
          endif
          OMEGA=angi
        endif
c...
c...  l=4
c...
        if(lone.eq.4)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(21)*(4.65269135862698D1*xyzi(i,j,k+4)-3.98802116453
     $        7411D1*xyzi(i,j,k+2)+3.988021164537411D0*xyzi(i,j,k))+angi
            angi=yrk(23)*(-4.161493662486312D1*xyzi(i+4,j,k)+3.566994567
     $        84541D1*xyzi(i+2,j,k)+4.161493662486312D1*xyzi(i,j+4,k)-3.
     $        56699456784541D1*xyzi(i,j+2,k))+angi
            angi=yrk(25)*(7.864483795364389D0*xyzi(i+4,j,k)-4.7186902772
     $        18633D1*xyzi(i+2,j+2,k)+7.864483795364389D0*xyzi(i,j+4,k))
     $        +angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(17)*(3.145793518145755D1*xyzi(i+3,j+1,k)-3.14579351
     $        8145755D1*xyzi(i+1,j+3,k))+angi
            angi=yrk(19)*(-8.322987324972623D1*xyzi(i+3,j+1,k)-8.3229873
     $        24972623D1*xyzi(i+1,j+3,k)+7.13398913569082D1*xyzi(i+1,j+1
     $        ,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(22)*(5.885240777217825D1*xyzi(i+3,j,k+1)+5.88524077
     $        7217825D1*xyzi(i+1,j+2,k+1)-3.362994729838757D1*xyzi(i+1,j
     $        ,k+1))+angi
            angi=yrk(24)*(6.67323578668065D1*xyzi(i+1,j+2,k+1)-2.2244119
     $        2889355D1*xyzi(i+3,j,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(18)*(2.22441192889355D1*xyzi(i,j+3,k+1)-6.673235786
     $        68065D1*xyzi(i+2,j+1,k+1))+angi
            angi=yrk(20)*(5.885240777217825D1*xyzi(i+2,j+1,k+1)+5.885240
     $        777217825D1*xyzi(i,j+3,k+1)-3.362994729838757D1*xyzi(i,j+1
     $        ,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=6
c...
        if(lone.eq.6)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(43)*(1.845306898868157D2*xyzi(i,j,k+6)-2.5163275893
     $        65668D2*xyzi(i,j,k+4)+8.387758631218894D1*xyzi(i,j,k+2)-3.
     $        994170776770902D0*xyzi(i,j,k))+angi
            angi=yrk(45)*(1.910074105988639D2*xyzi(i+6,j,k)+1.9100741059
     $        88639D2*xyzi(i+4,j+2,k)-2.778289608710748D2*xyzi(i+4,j,k)-
     $        1.910074105988639D2*xyzi(i+2,j+4,k)+9.260965362369161D1*xy
     $        zi(i+2,j,k)-1.910074105988639D2*xyzi(i,j+6,k)+2.7782896087
     $        10748D2*xyzi(i,j+4,k)-9.260965362369161D1*xyzi(i,j+2,k))+a
     $        ngi
            angi=yrk(47)*(-6.974604495709942D1*xyzi(i+6,j,k)+3.487302247
     $        854971D2*xyzi(i+4,j+2,k)+6.340549541554493D1*xyzi(i+4,j,k)
     $        +3.487302247854971D2*xyzi(i+2,j+4,k)-3.804329724932696D2*x
     $        yzi(i+2,j+2,k)-6.974604495709942D1*xyzi(i,j+6,k)+6.3405495
     $        41554493D1*xyzi(i,j+4,k))+angi
            angi=yrk(49)*(8.585144663680938D0*xyzi(i+6,j,k)-1.2877716995
     $        52141D2*xyzi(i+4,j+2,k)+1.287771699552141D2*xyzi(i+2,j+4,k
     $        )-8.585144663680938D0*xyzi(i,j+6,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(37)*(5.151086798208563D1*xyzi(i+5,j+1,k)-1.71702893
     $        2736188D2*xyzi(i+3,j+3,k)+5.151086798208563D1*xyzi(i+1,j+5
     $        ,k))+angi
            angi=yrk(39)*(-2.789841798283977D2*xyzi(i+5,j+1,k)+2.5362198
     $        16621797D2*xyzi(i+3,j+1,k)+2.789841798283977D2*xyzi(i+1,j+
     $        5,k)-2.536219816621797D2*xyzi(i+1,j+3,k))+angi
            angi=yrk(41)*(3.820148211977279D2*xyzi(i+5,j+1,k)+7.64029642
     $        3954557D2*xyzi(i+3,j+3,k)-5.556579217421496D2*xyzi(i+3,j+1
     $        ,k)+3.820148211977279D2*xyzi(i+1,j+5,k)-5.556579217421496D
     $        2*xyzi(i+1,j+3,k)+1.852193072473832D2*xyzi(i+1,j+1,k))+ang
     $        i
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(44)*(-2.416073869853585D2*xyzi(i+5,j,k+1)-4.8321477
     $        39707171D2*xyzi(i+3,j+2,k+1)+2.635716948931184D2*xyzi(i+3,
     $        j,k+1)-2.416073869853585D2*xyzi(i+1,j+4,k+1)+2.63571694893
     $        1184D2*xyzi(i+1,j+2,k+1)-5.857148775402631D1*xyzi(i+1,j,k+
     $        1))+angi
            angi=yrk(46)*(1.27338273732576D2*xyzi(i+5,j,k+1)-2.546765474
     $        651519D2*xyzi(i+3,j+2,k+1)-9.260965362369161D1*xyzi(i+3,j,
     $        k+1)-3.820148211977279D2*xyzi(i+1,j+4,k+1)+2.7782896087107
     $        48D2*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(48)*(-2.973981349564841D1*xyzi(i+5,j,k+1)+2.9739813
     $        49564841D2*xyzi(i+3,j+2,k+1)-1.486990674782421D2*xyzi(i+1,
     $        j+4,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(38)*(-1.486990674782421D2*xyzi(i+4,j+1,k+1)+2.97398
     $        1349564841D2*xyzi(i+2,j+3,k+1)-2.973981349564841D1*xyzi(i,
     $        j+5,k+1))+angi
            angi=yrk(40)*(3.820148211977279D2*xyzi(i+4,j+1,k+1)+2.546765
     $        474651519D2*xyzi(i+2,j+3,k+1)-2.778289608710748D2*xyzi(i+2
     $        ,j+1,k+1)-1.27338273732576D2*xyzi(i,j+5,k+1)+9.26096536236
     $        9161D1*xyzi(i,j+3,k+1))+angi
            angi=yrk(42)*(-2.416073869853585D2*xyzi(i+4,j+1,k+1)-4.83214
     $        7739707171D2*xyzi(i+2,j+3,k+1)+2.635716948931184D2*xyzi(i+
     $        2,j+1,k+1)-2.416073869853585D2*xyzi(i,j+5,k+1)+2.635716948
     $        931184D2*xyzi(i,j+3,k+1)-5.857148775402631D1*xyzi(i,j+1,k+
     $        1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=8
c...
        if(lone.eq.8)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(73)*(7.347980147805839D2*xyzi(i,j,k+8)-1.3716229609
     $        23757D3*xyzi(i,j,k+6)+7.91320938994475D2*xyzi(i,j,k+4)-1.4
     $        38765343626318D2*xyzi(i,j,k+2)+3.996570398961995D0*xyzi(i,
     $        j,k))+angi
            angi=yrk(75)*(-8.197015020580125D2*xyzi(i+8,j,k)-1.639403004
     $        116025D3*xyzi(i+6,j+2,k)+1.639403004116025D3*xyzi(i+6,j,k)
     $        +1.639403004116025D3*xyzi(i+4,j+2,k)-1.008863387148323D3*x
     $        yzi(i+4,j,k)+1.639403004116025D3*xyzi(i+2,j+6,k)-1.6394030
     $        04116025D3*xyzi(i+2,j+4,k)+1.834297067542406D2*xyzi(i+2,j,
     $        k)+8.197015020580125D2*xyzi(i,j+8,k)-1.639403004116025D3*x
     $        yzi(i,j+6,k)+1.008863387148323D3*xyzi(i,j+4,k)-1.834297067
     $        542406D2*xyzi(i,j+2,k))+angi
            angi=yrk(77)*(3.907773582803669D2*xyzi(i+8,j,k)-1.5631094331
     $        21468D3*xyzi(i+6,j+2,k)-6.252437732485871D2*xyzi(i+6,j,k)-
     $        3.907773582803669D3*xyzi(i+4,j+4,k)+3.126218866242935D3*xy
     $        zi(i+4,j+2,k)+2.404783743263796D2*xyzi(i+4,j,k)-1.56310943
     $        3121468D3*xyzi(i+2,j+6,k)+3.126218866242935D3*xyzi(i+2,j+4
     $        ,k)-1.442870245958278D3*xyzi(i+2,j+2,k)+3.907773582803669D
     $        2*xyzi(i,j+8,k)-6.252437732485871D2*xyzi(i,j+6,k)+2.404783
     $        743263796D2*xyzi(i,j+4,k))+angi
            angi=yrk(79)*(-1.003423624270676D2*xyzi(i+8,j,k)+1.404793073
     $        978946D3*xyzi(i+6,j+2,k)+9.365287159859641D1*xyzi(i+6,j,k)
     $        -1.404793073978946D3*xyzi(i+4,j+2,k)-1.404793073978946D3*x
     $        yzi(i+2,j+6,k)+1.404793073978946D3*xyzi(i+2,j+4,k)+1.00342
     $        3624270676D2*xyzi(i,j+8,k)-9.365287159859641D1*xyzi(i,j+6,
     $        k))+angi
            angi=yrk(81)*(9.159962562443957D0*xyzi(i+8,j,k)-2.5647895174
     $        84308D2*xyzi(i+6,j+2,k)+6.41197379371077D2*xyzi(i+4,j+4,k)
     $        -2.564789517484308D2*xyzi(i+2,j+6,k)+9.159962562443957D0*x
     $        yzi(i,j+8,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(65)*(7.327970049955166D1*xyzi(i+7,j+1,k)-5.12957903
     $        4968616D2*xyzi(i+5,j+3,k)+5.129579034968616D2*xyzi(i+3,j+5
     $        ,k)-7.327970049955166D1*xyzi(i+1,j+7,k))+angi
            angi=yrk(67)*(-6.020541745624055D2*xyzi(i+7,j+1,k)+1.4047930
     $        73978946D3*xyzi(i+5,j+3,k)+5.619172295915784D2*xyzi(i+5,j+
     $        1,k)+1.404793073978946D3*xyzi(i+3,j+5,k)-1.873057431971928
     $        D3*xyzi(i+3,j+3,k)-6.020541745624055D2*xyzi(i+1,j+7,k)+5.6
     $        19172295915784D2*xyzi(i+1,j+5,k))+angi
            angi=yrk(69)*(1.563109433121468D3*xyzi(i+7,j+1,k)+1.56310943
     $        3121468D3*xyzi(i+5,j+3,k)-2.500975092994348D3*xyzi(i+5,j+1
     $        ,k)-1.563109433121468D3*xyzi(i+3,j+5,k)+9.619134973055186D
     $        2*xyzi(i+3,j+1,k)-1.563109433121468D3*xyzi(i+1,j+7,k)+2.50
     $        0975092994348D3*xyzi(i+1,j+5,k)-9.619134973055186D2*xyzi(i
     $        +1,j+3,k))+angi
            angi=yrk(71)*(-1.639403004116025D3*xyzi(i+7,j+1,k)-4.9182090
     $        12348075D3*xyzi(i+5,j+3,k)+3.27880600823205D3*xyzi(i+5,j+1
     $        ,k)-4.918209012348075D3*xyzi(i+3,j+5,k)+6.5576120164641D3*
     $        xyzi(i+3,j+3,k)-2.017726774296646D3*xyzi(i+3,j+1,k)-1.6394
     $        03004116025D3*xyzi(i+1,j+7,k)+3.27880600823205D3*xyzi(i+1,
     $        j+5,k)-2.017726774296646D3*xyzi(i+1,j+3,k)+3.6685941350848
     $        11D2*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(74)*(9.797306863741119D2*xyzi(i+7,j,k+1)+2.93919205
     $        9122336D3*xyzi(i+5,j+2,k+1)-1.567569098198579D3*xyzi(i+5,j
     $        ,k+1)+2.939192059122336D3*xyzi(i+3,j+4,k+1)-3.135138196397
     $        158D3*xyzi(i+3,j+2,k+1)+7.234934299378057D2*xyzi(i+3,j,k+1
     $        )+9.797306863741119D2*xyzi(i+1,j+6,k+1)-1.567569098198579D
     $        3*xyzi(i+1,j+4,k+1)+7.234934299378057D2*xyzi(i+1,j+2,k+1)-
     $        8.769617332579463D1*xyzi(i+1,j,k+1))+angi
            angi=yrk(76)*(-6.05389680277916D2*xyzi(i+7,j,k+1)+6.05389680
     $        277916D2*xyzi(i+5,j+2,k+1)+8.071862403705547D2*xyzi(i+5,j,
     $        k+1)+3.02694840138958D3*xyzi(i+3,j+4,k+1)-1.61437248074110
     $        9D3*xyzi(i+3,j+2,k+1)-2.483649970370938D2*xyzi(i+3,j,k+1)+
     $        1.816169040833748D3*xyzi(i+1,j+6,k+1)-2.421558721111664D3*
     $        xyzi(i+1,j+4,k+1)+7.450949911112813D2*xyzi(i+1,j+2,k+1))+a
     $        ngi
            angi=yrk(78)*(2.167642773184962D2*xyzi(i+7,j,k+1)-1.95087849
     $        5866466D3*xyzi(i+5,j+2,k+1)-1.73411421854797D2*xyzi(i+5,j,
     $        k+1)-1.083821386592481D3*xyzi(i+3,j+4,k+1)+1.7341142185479
     $        7D3*xyzi(i+3,j+2,k+1)+1.083821386592481D3*xyzi(i+1,j+6,k+1
     $        )-8.670571092739848D2*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(80)*(-3.663985024977583D1*xyzi(i+7,j,k+1)+7.6943685
     $        52452924D2*xyzi(i+5,j+2,k+1)-1.282394758742154D3*xyzi(i+3,
     $        j+4,k+1)+2.564789517484308D2*xyzi(i+1,j+6,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(66)*(-2.564789517484308D2*xyzi(i+6,j+1,k+1)+1.28239
     $        4758742154D3*xyzi(i+4,j+3,k+1)-7.694368552452924D2*xyzi(i+
     $        2,j+5,k+1)+3.663985024977583D1*xyzi(i,j+7,k+1))+angi
            angi=yrk(68)*(1.083821386592481D3*xyzi(i+6,j+1,k+1)-1.083821
     $        386592481D3*xyzi(i+4,j+3,k+1)-8.670571092739848D2*xyzi(i+4
     $        ,j+1,k+1)-1.950878495866466D3*xyzi(i+2,j+5,k+1)+1.73411421
     $        854797D3*xyzi(i+2,j+3,k+1)+2.167642773184962D2*xyzi(i,j+7,
     $        k+1)-1.73411421854797D2*xyzi(i,j+5,k+1))+angi
            angi=yrk(70)*(-1.816169040833748D3*xyzi(i+6,j+1,k+1)-3.02694
     $        840138958D3*xyzi(i+4,j+3,k+1)+2.421558721111664D3*xyzi(i+4
     $        ,j+1,k+1)-6.05389680277916D2*xyzi(i+2,j+5,k+1)+1.614372480
     $        741109D3*xyzi(i+2,j+3,k+1)-7.450949911112813D2*xyzi(i+2,j+
     $        1,k+1)+6.05389680277916D2*xyzi(i,j+7,k+1)-8.07186240370554
     $        7D2*xyzi(i,j+5,k+1)+2.483649970370938D2*xyzi(i,j+3,k+1))+a
     $        ngi
            angi=yrk(72)*(9.797306863741119D2*xyzi(i+6,j+1,k+1)+2.939192
     $        059122336D3*xyzi(i+4,j+3,k+1)-1.567569098198579D3*xyzi(i+4
     $        ,j+1,k+1)+2.939192059122336D3*xyzi(i+2,j+5,k+1)-3.13513819
     $        6397158D3*xyzi(i+2,j+3,k+1)+7.234934299378057D2*xyzi(i+2,j
     $        +1,k+1)+9.797306863741119D2*xyzi(i,j+7,k+1)-1.567569098198
     $        579D3*xyzi(i,j+5,k+1)+7.234934299378057D2*xyzi(i,j+3,k+1)-
     $        8.769617332579463D1*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=10
c...
        if(lone.eq.10)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(111)*(2.930982152135684D3*xyzi(i,j,k+10)-6.94179983
     $        4005568D3*xyzi(i,j,k+8)+5.716776333886939D3*xyzi(i,j,k+6)-
     $        1.905592111295646D3*xyzi(i,j,k+4)+2.198760128418053D2*xyzi
     $        (i,j,k+2)-3.997745688032824D0*xyzi(i,j,k))+angi
            angi=yrk(113)*(3.422649766190848D3*xyzi(i+10,j,k)+1.02679492
     $        9857254D4*xyzi(i+8,j+2,k)-8.646694146166353D3*xyzi(i+8,j,k
     $        )+6.845299532381696D3*xyzi(i+6,j+4,k)-1.729338829233271D4*
     $        xyzi(i+6,j+2,k)+7.629436011323252D3*xyzi(i+6,j,k)-6.845299
     $        532381696D3*xyzi(i+4,j+6,k)+7.629436011323252D3*xyzi(i+4,j
     $        +2,k)-2.712688359581601D3*xyzi(i+4,j,k)-1.026794929857254D
     $        4*xyzi(i+2,j+8,k)+1.729338829233271D4*xyzi(i+2,j+6,k)-7.62
     $        9436011323252D3*xyzi(i+2,j+4,k)+3.130025030286462D2*xyzi(i
     $        +2,j,k)-3.422649766190848D3*xyzi(i,j+10,k)+8.6466941461663
     $        53D3*xyzi(i,j+8,k)-7.629436011323252D3*xyzi(i,j+6,k)+2.712
     $        688359581601D3*xyzi(i,j+4,k)-3.130025030286462D2*xyzi(i,j+
     $        2,k))+angi
            angi=yrk(115)*(-1.898544496916298D3*xyzi(i+10,j,k)+5.6956334
     $        90748894D3*xyzi(i+8,j+2,k)+4.196782572130764D3*xyzi(i+8,j,
     $        k)+2.657962295682817D4*xyzi(i+6,j+4,k)-1.678713028852305D4
     $        *xyzi(i+6,j+2,k)-2.962434756798186D3*xyzi(i+6,j,k)+2.65796
     $        2295682817D4*xyzi(i+4,j+6,k)-4.196782572130764D4*xyzi(i+4,
     $        j+4,k)+1.481217378399093D4*xyzi(i+4,j+2,k)+6.5831883484404
     $        14D2*xyzi(i+4,j,k)+5.695633490748894D3*xyzi(i+2,j+8,k)-1.6
     $        78713028852305D4*xyzi(i+2,j+6,k)+1.481217378399093D4*xyzi(
     $        i+2,j+4,k)-3.949913009064248D3*xyzi(i+2,j+2,k)-1.898544496
     $        916298D3*xyzi(i,j+10,k)+4.196782572130764D3*xyzi(i,j+8,k)-
     $        2.962434756798186D3*xyzi(i,j+6,k)+6.583188348440414D2*xyzi
     $        (i,j+4,k))+angi
            angi=yrk(117)*(6.712368440769583D2*xyzi(i+10,j,k)-8.72607897
     $        3000458D3*xyzi(i+8,j+2,k)-1.130504158445403D3*xyzi(i+8,j,k
     $        )-9.397315817077416D3*xyzi(i+6,j+4,k)+1.582705821823565D4*
     $        xyzi(i+6,j+2,k)+4.655017123010485D2*xyzi(i+6,j,k)+9.397315
     $        817077416D3*xyzi(i+4,j+6,k)-6.982525684515727D3*xyzi(i+4,j
     $        +2,k)+8.726078973000458D3*xyzi(i+2,j+8,k)-1.58270582182356
     $        5D4*xyzi(i+2,j+6,k)+6.982525684515727D3*xyzi(i+2,j+4,k)-6.
     $        712368440769583D2*xyzi(i,j+10,k)+1.130504158445403D3*xyzi(
     $        i,j+8,k)-4.655017123010485D2*xyzi(i,j+6,k))+angi
            angi=yrk(119)*(-1.329247023836435D2*xyzi(i+10,j,k)+3.5889669
     $        64358373D3*xyzi(i+8,j+2,k)+1.259286654160833D2*xyzi(i+8,j,
     $        k)-5.582837500113025D3*xyzi(i+6,j+4,k)-3.526002631650332D3
     $        *xyzi(i+6,j+2,k)-5.582837500113025D3*xyzi(i+4,j+6,k)+8.815
     $        006579125829D3*xyzi(i+4,j+4,k)+3.588966964358373D3*xyzi(i+
     $        2,j+8,k)-3.526002631650332D3*xyzi(i+2,j+6,k)-1.32924702383
     $        6435D2*xyzi(i,j+10,k)+1.259286654160833D2*xyzi(i,j+8,k))+a
     $        ngi
            angi=yrk(121)*(9.643371463227499D0*xyzi(i+10,j,k)-4.33951715
     $        8452375D2*xyzi(i+8,j+2,k)+2.025108007277775D3*xyzi(i+6,j+4
     $        ,k)-2.025108007277775D3*xyzi(i+4,j+6,k)+4.339517158452375D
     $        2*xyzi(i+2,j+8,k)-9.643371463227499D0*xyzi(i,j+10,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(101)*(9.643371463227499D1*xyzi(i+9,j+1,k)-1.1572045
     $        755873D3*xyzi(i+7,j+3,k)+2.43012960873333D3*xyzi(i+5,j+5,k
     $        )-1.1572045755873D3*xyzi(i+3,j+7,k)+9.643371463227499D1*xy
     $        zi(i+1,j+9,k))+angi
            angi=yrk(103)*(-1.063397619069148D3*xyzi(i+9,j+1,k)+6.380385
     $        714414886D3*xyzi(i+7,j+3,k)+1.007429323328666D3*xyzi(i+7,j
     $        +1,k)-7.052005263300663D3*xyzi(i+5,j+3,k)-6.38038571441488
     $        6D3*xyzi(i+3,j+7,k)+7.052005263300663D3*xyzi(i+3,j+5,k)+1.
     $        063397619069148D3*xyzi(i+1,j+9,k)-1.007429323328666D3*xyzi
     $        (i+1,j+7,k))+angi
            angi=yrk(105)*(4.02742106446175D3*xyzi(i+9,j+1,k)-5.36989475
     $        2615666D3*xyzi(i+7,j+3,k)-6.783024950672421D3*xyzi(i+7,j+1
     $        ,k)-1.879463163415483D4*xyzi(i+5,j+5,k)+1.582705821823565D
     $        4*xyzi(i+5,j+3,k)+2.793010273806291D3*xyzi(i+5,j+1,k)-5.36
     $        9894752615666D3*xyzi(i+3,j+7,k)+1.582705821823565D4*xyzi(i
     $        +3,j+5,k)-9.31003424602097D3*xyzi(i+3,j+3,k)+4.02742106446
     $        175D3*xyzi(i+1,j+9,k)-6.783024950672421D3*xyzi(i+1,j+7,k)+
     $        2.793010273806291D3*xyzi(i+1,j+5,k))+angi
            angi=yrk(107)*(-7.594177987665191D3*xyzi(i+9,j+1,k)-1.518835
     $        597533038D4*xyzi(i+7,j+3,k)+1.678713028852305D4*xyzi(i+7,j
     $        +1,k)+1.678713028852305D4*xyzi(i+5,j+3,k)-1.18497390271927
     $        4D4*xyzi(i+5,j+1,k)+1.518835597533038D4*xyzi(i+3,j+7,k)-1.
     $        678713028852305D4*xyzi(i+3,j+5,k)+2.633275339376166D3*xyzi
     $        (i+3,j+1,k)+7.594177987665191D3*xyzi(i+1,j+9,k)-1.67871302
     $        8852305D4*xyzi(i+1,j+7,k)+1.184973902719274D4*xyzi(i+1,j+5
     $        ,k)-2.633275339376166D3*xyzi(i+1,j+3,k))+angi
            angi=yrk(109)*(6.845299532381696D3*xyzi(i+9,j+1,k)+2.7381198
     $        12952678D4*xyzi(i+7,j+3,k)-1.729338829233271D4*xyzi(i+7,j+
     $        1,k)+4.107179719429018D4*xyzi(i+5,j+5,k)-5.188016487699811
     $        D4*xyzi(i+5,j+3,k)+1.52588720226465D4*xyzi(i+5,j+1,k)+2.73
     $        8119812952678D4*xyzi(i+3,j+7,k)-5.188016487699811D4*xyzi(i
     $        +3,j+5,k)+3.051774404529301D4*xyzi(i+3,j+3,k)-5.4253767191
     $        63202D3*xyzi(i+3,j+1,k)+6.845299532381696D3*xyzi(i+1,j+9,k
     $        )-1.729338829233271D4*xyzi(i+1,j+7,k)+1.52588720226465D4*x
     $        yzi(i+1,j+5,k)-5.425376719163202D3*xyzi(i+1,j+3,k)+6.26005
     $        0060572925D2*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(112)*(-3.952135527704191D3*xyzi(i+9,j,k+1)-1.580854
     $        211081677D4*xyzi(i+7,j+2,k+1)+8.320285321482508D3*xyzi(i+7
     $        ,j,k+1)-2.371281316622515D4*xyzi(i+5,j+4,k+1)+2.4960855964
     $        44752D4*xyzi(i+5,j+2,k+1)-5.873142579870006D3*xyzi(i+5,j,k
     $        +1)-1.580854211081677D4*xyzi(i+3,j+6,k+1)+2.49608559644475
     $        2D4*xyzi(i+3,j+4,k+1)-1.174628515974001D4*xyzi(i+3,j+2,k+1
     $        )+1.566171354632002D3*xyzi(i+3,j,k+1)-3.952135527704191D3*
     $        xyzi(i+1,j+8,k+1)+8.320285321482508D3*xyzi(i+1,j+6,k+1)-5.
     $        873142579870006D3*xyzi(i+1,j+4,k+1)+1.566171354632002D3*xy
     $        zi(i+1,j+2,k+1)-1.20474719587077D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(114)*(2.684947376307833D3*xyzi(i+9,j,k+1)-5.0872687
     $        13004316D3*xyzi(i+7,j,k+1)-1.6109684257847D4*xyzi(i+5,j+4,
     $        k+1)+5.087268713004316D3*xyzi(i+5,j+2,k+1)+2.9925110076495
     $        97D3*xyzi(i+5,j,k+1)-2.147957901046267D4*xyzi(i+3,j+6,k+1)
     $        +2.543634356502158D4*xyzi(i+3,j+4,k+1)-5.985022015299195D3
     $        *xyzi(i+3,j+2,k+1)-5.32001956915484D2*xyzi(i+3,j,k+1)-8.05
     $        48421289235D3*xyzi(i+1,j+8,k+1)+1.526180613901295D4*xyzi(i
     $        +1,j+6,k+1)-8.977533022948792D3*xyzi(i+1,j+4,k+1)+1.596005
     $        870746452D3*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(116)*(-1.200744969886805D3*xyzi(i+9,j,k+1)+9.605959
     $        759094437D3*xyzi(i+7,j+2,k+1)+1.769518902991081D3*xyzi(i+7
     $        ,j,k+1)+1.681042957841527D4*xyzi(i+5,j+4,k+1)-1.5925670126
     $        91973D4*xyzi(i+5,j+2,k+1)-6.245360834086167D2*xyzi(i+5,j,k
     $        +1)-8.847594514955403D3*xyzi(i+3,j+4,k+1)+6.24536083408616
     $        7D3*xyzi(i+3,j+2,k+1)-6.003724849434023D3*xyzi(i+1,j+8,k+1
     $        )+8.847594514955403D3*xyzi(i+1,j+6,k+1)-3.122680417043083D
     $        3*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(118)*(3.255976950512413D2*xyzi(i+9,j,k+1)-6.5119539
     $        01024826D3*xyzi(i+7,j+2,k+1)-2.741875326747295D2*xyzi(i+7,
     $        j,k+1)+4.558367730717379D3*xyzi(i+5,j+4,k+1)+5.75793818616
     $        932D3*xyzi(i+5,j+2,k+1)+9.116735461434757D3*xyzi(i+3,j+6,k
     $        +1)-9.596563643615534D3*xyzi(i+3,j+4,k+1)-2.27918386535868
     $        9D3*xyzi(i+1,j+8,k+1)+1.919312728723107D3*xyzi(i+1,j+6,k+1
     $        ))+angi
            angi=yrk(120)*(-4.31264682481166D1*xyzi(i+9,j,k+1)+1.5525528
     $        56932198D3*xyzi(i+7,j+2,k+1)-5.433934999262692D3*xyzi(i+5,
     $        j+4,k+1)+3.622623332841795D3*xyzi(i+3,j+6,k+1)-3.881382142
     $        330494D2*xyzi(i+1,j+8,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(102)*(-3.881382142330494D2*xyzi(i+8,j+1,k+1)+3.6226
     $        23332841795D3*xyzi(i+6,j+3,k+1)-5.433934999262692D3*xyzi(i
     $        +4,j+5,k+1)+1.552552856932198D3*xyzi(i+2,j+7,k+1)-4.312646
     $        82481166D1*xyzi(i,j+9,k+1))+angi
            angi=yrk(104)*(2.279183865358689D3*xyzi(i+8,j+1,k+1)-9.11673
     $        5461434757D3*xyzi(i+6,j+3,k+1)-1.919312728723107D3*xyzi(i+
     $        6,j+1,k+1)-4.558367730717379D3*xyzi(i+4,j+5,k+1)+9.5965636
     $        43615534D3*xyzi(i+4,j+3,k+1)+6.511953901024826D3*xyzi(i+2,
     $        j+7,k+1)-5.75793818616932D3*xyzi(i+2,j+5,k+1)-3.2559769505
     $        12413D2*xyzi(i,j+9,k+1)+2.741875326747295D2*xyzi(i,j+7,k+1
     $        ))+angi
            angi=yrk(106)*(-6.003724849434023D3*xyzi(i+8,j+1,k+1)+8.8475
     $        94514955403D3*xyzi(i+6,j+1,k+1)+1.681042957841527D4*xyzi(i
     $        +4,j+5,k+1)-8.847594514955403D3*xyzi(i+4,j+3,k+1)-3.122680
     $        417043083D3*xyzi(i+4,j+1,k+1)+9.605959759094437D3*xyzi(i+2
     $        ,j+7,k+1)-1.592567012691973D4*xyzi(i+2,j+5,k+1)+6.24536083
     $        4086167D3*xyzi(i+2,j+3,k+1)-1.200744969886805D3*xyzi(i,j+9
     $        ,k+1)+1.769518902991081D3*xyzi(i,j+7,k+1)-6.24536083408616
     $        7D2*xyzi(i,j+5,k+1))+angi
            angi=yrk(108)*(8.0548421289235D3*xyzi(i+8,j+1,k+1)+2.1479579
     $        01046267D4*xyzi(i+6,j+3,k+1)-1.526180613901295D4*xyzi(i+6,
     $        j+1,k+1)+1.6109684257847D4*xyzi(i+4,j+5,k+1)-2.54363435650
     $        2158D4*xyzi(i+4,j+3,k+1)+8.977533022948792D3*xyzi(i+4,j+1,
     $        k+1)-5.087268713004316D3*xyzi(i+2,j+5,k+1)+5.9850220152991
     $        95D3*xyzi(i+2,j+3,k+1)-1.596005870746452D3*xyzi(i+2,j+1,k+
     $        1)-2.684947376307833D3*xyzi(i,j+9,k+1)+5.087268713004316D3
     $        *xyzi(i,j+7,k+1)-2.992511007649597D3*xyzi(i,j+5,k+1)+5.320
     $        01956915484D2*xyzi(i,j+3,k+1))+angi
            angi=yrk(110)*(-3.952135527704191D3*xyzi(i+8,j+1,k+1)-1.5808
     $        54211081677D4*xyzi(i+6,j+3,k+1)+8.320285321482508D3*xyzi(i
     $        +6,j+1,k+1)-2.371281316622515D4*xyzi(i+4,j+5,k+1)+2.496085
     $        596444752D4*xyzi(i+4,j+3,k+1)-5.873142579870006D3*xyzi(i+4
     $        ,j+1,k+1)-1.580854211081677D4*xyzi(i+2,j+7,k+1)+2.49608559
     $        6444752D4*xyzi(i+2,j+5,k+1)-1.174628515974001D4*xyzi(i+2,j
     $        +3,k+1)+1.566171354632002D3*xyzi(i+2,j+1,k+1)-3.9521355277
     $        04191D3*xyzi(i,j+9,k+1)+8.320285321482508D3*xyzi(i,j+7,k+1
     $        )-5.873142579870006D3*xyzi(i,j+5,k+1)+1.566171354632002D3*
     $        xyzi(i,j+3,k+1)-1.20474719587077D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=12
c...
        if(lone.eq.12)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(157)*(1.170163993078432D4*xyzi(i,j,k+12)-3.35786189
     $        3181587D4*xyzi(i,j,k+10)+3.597709171265986D4*xyzi(i,j,k+8)
     $        -1.767295733253467D4*xyzi(i,j,k+6)+3.898446470412059D3*xyz
     $        i(i,j,k+4)-3.118757176329647D2*xyzi(i,j,k+2)+3.99840663632
     $        0061D0*xyzi(i,j,k))+angi
            angi=yrk(159)*(-1.409330983563191D4*xyzi(i+12,j,k)-5.6373239
     $        34252763D4*xyzi(i+10,j+2,k)+4.289268210844493D4*xyzi(i+10,
     $        j,k)-7.046654917815953D4*xyzi(i+8,j+4,k)+1.286780463253348
     $        D5*xyzi(i+8,j+2,k)-4.902020812393707D4*xyzi(i+8,j,k)+8.578
     $        536421688987D4*xyzi(i+6,j+4,k)-9.804041624787413D4*xyzi(i+
     $        6,j+2,k)+2.580010953891425D4*xyzi(i+6,j,k)+7.0466549178159
     $        53D4*xyzi(i+4,j+8,k)-8.578536421688987D4*xyzi(i+4,j+6,k)+2
     $        .580010953891425D4*xyzi(i+4,j+2,k)-6.070614009156293D3*xyz
     $        i(i+4,j,k)+5.637323934252763D4*xyzi(i+2,j+10,k)-1.28678046
     $        3253348D5*xyzi(i+2,j+8,k)+9.804041624787413D4*xyzi(i+2,j+6
     $        ,k)-2.580010953891425D4*xyzi(i+2,j+4,k)+4.856491207325035D
     $        2*xyzi(i+2,j,k)+1.409330983563191D4*xyzi(i,j+12,k)-4.28926
     $        8210844493D4*xyzi(i,j+10,k)+4.902020812393707D4*xyzi(i,j+8
     $        ,k)-2.580010953891425D4*xyzi(i,j+6,k)+6.070614009156293D3*
     $        xyzi(i,j+4,k)-4.856491207325035D2*xyzi(i,j+2,k))+angi
            angi=yrk(161)*(8.630354471061409D3*xyzi(i+12,j,k)-1.72607089
     $        4212282D4*xyzi(i+10,j+2,k)-2.401489939773609D4*xyzi(i+10,j
     $        ,k)-1.467160260080439D5*xyzi(i+8,j+4,k)+7.204469819320828D
     $        4*xyzi(i+8,j+2,k)+2.401489939773609D4*xyzi(i+8,j,k)-2.4164
     $        99251897194D5*xyzi(i+6,j+6,k)+3.362085915683053D5*xyzi(i+6
     $        ,j+4,k)-9.605959759094437D4*xyzi(i+6,j+2,k)-1.011153658852
     $        046D4*xyzi(i+6,j,k)-1.467160260080439D5*xyzi(i+4,j+8,k)+3.
     $        362085915683053D5*xyzi(i+4,j+6,k)-2.401489939773609D5*xyzi
     $        (i+4,j+4,k)+5.05576829426023D4*xyzi(i+4,j+2,k)+1.486990674
     $        782421D3*xyzi(i+4,j,k)-1.726070894212282D4*xyzi(i+2,j+10,k
     $        )+7.204469819320828D4*xyzi(i+2,j+8,k)-9.605959759094437D4*
     $        xyzi(i+2,j+6,k)+5.05576829426023D4*xyzi(i+2,j+4,k)-8.92194
     $        4048694524D3*xyzi(i+2,j+2,k)+8.630354471061409D3*xyzi(i,j+
     $        12,k)-2.401489939773609D4*xyzi(i,j+10,k)+2.401489939773609
     $        D4*xyzi(i,j+8,k)-1.011153658852046D4*xyzi(i,j+6,k)+1.48699
     $        0674782421D3*xyzi(i,j+4,k))+angi
            angi=yrk(163)*(-3.692002053806592D3*xyzi(i+12,j,k)+4.4304024
     $        6456791D4*xyzi(i+10,j+2,k)+8.668178735024172D3*xyzi(i+10,j
     $        ,k)+9.968405545277798D4*xyzi(i+8,j+4,k)-1.126863235553142D
     $        5*xyzi(i+8,j+2,k)-6.604326655256512D3*xyzi(i+8,j,k)-1.2135
     $        45022903384D5*xyzi(i+6,j+4,k)+9.246057317359117D4*xyzi(i+6
     $        ,j+2,k)+1.622115318834933D3*xyzi(i+6,j,k)-9.96840554527779
     $        8D4*xyzi(i+4,j+8,k)+1.213545022903384D5*xyzi(i+4,j+6,k)-2.
     $        433172978252399D4*xyzi(i+4,j+2,k)-4.43040246456791D4*xyzi(
     $        i+2,j+10,k)+1.126863235553142D5*xyzi(i+2,j+8,k)-9.24605731
     $        7359117D4*xyzi(i+2,j+6,k)+2.433172978252399D4*xyzi(i+2,j+4
     $        ,k)+3.692002053806592D3*xyzi(i,j+12,k)-8.668178735024172D3
     $        *xyzi(i,j+10,k)+6.604326655256512D3*xyzi(i,j+8,k)-1.622115
     $        318834933D3*xyzi(i,j+6,k))+angi
            angi=yrk(165)*(1.037363021977718D3*xyzi(i+12,j,k)-2.69714385
     $        7142068D4*xyzi(i+10,j+2,k)-1.80410960343951D3*xyzi(i+10,j,
     $        k)+1.556044532966577D4*xyzi(i+8,j+4,k)+4.871095929286677D4
     $        *xyzi(i+8,j+2,k)+7.731898300455043D2*xyzi(i+8,j,k)+8.71384
     $        9384612834D4*xyzi(i+6,j+6,k)-7.577260334445942D4*xyzi(i+6,
     $        j+4,k)-2.164931524127412D4*xyzi(i+6,j+2,k)+1.5560445329665
     $        77D4*xyzi(i+4,j+8,k)-7.577260334445942D4*xyzi(i+4,j+6,k)+5
     $        .41232881031853D4*xyzi(i+4,j+4,k)-2.697143857142068D4*xyzi
     $        (i+2,j+10,k)+4.871095929286677D4*xyzi(i+2,j+8,k)-2.1649315
     $        24127412D4*xyzi(i+2,j+6,k)+1.037363021977718D3*xyzi(i,j+12
     $        ,k)-1.80410960343951D3*xyzi(i,j+10,k)+7.731898300455043D2*
     $        xyzi(i,j+8,k))+angi
            angi=yrk(167)*(-1.671861890280821D2*xyzi(i+12,j,k)+7.3561923
     $        17235614D3*xyzi(i+10,j+2,k)+1.599172242877307D2*xyzi(i+10,
     $        j,k)-2.758572118963355D4*xyzi(i+8,j+4,k)-7.196275092947883
     $        D3*xyzi(i+8,j+2,k)+3.358261710042346D4*xyzi(i+6,j+4,k)+2.7
     $        58572118963355D4*xyzi(i+4,j+8,k)-3.358261710042346D4*xyzi(
     $        i+4,j+6,k)-7.356192317235614D3*xyzi(i+2,j+10,k)+7.19627509
     $        2947883D3*xyzi(i+2,j+8,k)+1.671861890280821D2*xyzi(i,j+12,
     $        k)-1.599172242877307D2*xyzi(i,j+10,k))+angi
            angi=yrk(169)*(1.006342599515217D1*xyzi(i+12,j,k)-6.64186115
     $        6800431D2*xyzi(i+10,j+2,k)+4.981395867600323D3*xyzi(i+8,j+
     $        4,k)-9.298605619520603D3*xyzi(i+6,j+6,k)+4.981395867600323
     $        D3*xyzi(i+4,j+8,k)-6.641861156800431D2*xyzi(i+2,j+10,k)+1.
     $        006342599515217D1*xyzi(i,j+12,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(145)*(1.20761111941826D2*xyzi(i+11,j+1,k)-2.2139537
     $        18933477D3*xyzi(i+9,j+3,k)+7.970233388160517D3*xyzi(i+7,j+
     $        5,k)-7.970233388160517D3*xyzi(i+5,j+7,k)+2.213953718933477
     $        D3*xyzi(i+3,j+9,k)-1.20761111941826D2*xyzi(i+1,j+11,k))+an
     $        gi
            angi=yrk(147)*(-1.671861890280821D3*xyzi(i+11,j+1,k)+1.83904
     $        8079308904D4*xyzi(i+9,j+3,k)+1.599172242877307D3*xyzi(i+9,
     $        j+1,k)-2.206857695170684D4*xyzi(i+7,j+5,k)-1.9190066914527
     $        69D4*xyzi(i+7,j+3,k)-2.206857695170684D4*xyzi(i+5,j+7,k)+4
     $        .029914052050815D4*xyzi(i+5,j+5,k)+1.839048079308904D4*xyz
     $        i(i+3,j+9,k)-1.919006691452769D4*xyzi(i+3,j+7,k)-1.6718618
     $        90280821D3*xyzi(i+1,j+11,k)+1.599172242877307D3*xyzi(i+1,j
     $        +9,k))+angi
            angi=yrk(149)*(8.298904175821746D3*xyzi(i+11,j+1,k)-4.149452
     $        087910873D4*xyzi(i+9,j+3,k)-1.443287682751608D4*xyzi(i+9,j
     $        +1,k)-4.979342505493048D4*xyzi(i+7,j+5,k)+8.65972609650964
     $        8D4*xyzi(i+7,j+3,k)+6.185518640364035D3*xyzi(i+7,j+1,k)+4.
     $        979342505493048D4*xyzi(i+5,j+7,k)-4.329863048254824D4*xyzi
     $        (i+5,j+3,k)+4.149452087910873D4*xyzi(i+3,j+9,k)-8.65972609
     $        6509648D4*xyzi(i+3,j+7,k)+4.329863048254824D4*xyzi(i+3,j+5
     $        ,k)-8.298904175821746D3*xyzi(i+1,j+11,k)+1.443287682751608
     $        D4*xyzi(i+1,j+9,k)-6.185518640364035D3*xyzi(i+1,j+7,k))+an
     $        gi
            angi=yrk(151)*(-2.215201232283955D4*xyzi(i+11,j+1,k)+7.38400
     $        4107613184D3*xyzi(i+9,j+3,k)+5.200907241014503D4*xyzi(i+9,
     $        j+1,k)+1.329120739370373D5*xyzi(i+7,j+5,k)-6.9345429880193
     $        38D4*xyzi(i+7,j+3,k)-3.962595993153907D4*xyzi(i+7,j+1,k)+1
     $        .329120739370373D5*xyzi(i+5,j+7,k)-2.427090045806768D5*xyz
     $        i(i+5,j+5,k)+9.246057317359117D4*xyzi(i+5,j+3,k)+9.7326919
     $        13009597D3*xyzi(i+5,j+1,k)+7.384004107613184D3*xyzi(i+3,j+
     $        9,k)-6.934542988019338D4*xyzi(i+3,j+7,k)+9.246057317359117
     $        D4*xyzi(i+3,j+5,k)-3.244230637669866D4*xyzi(i+3,j+3,k)-2.2
     $        15201232283955D4*xyzi(i+1,j+11,k)+5.200907241014503D4*xyzi
     $        (i+1,j+9,k)-3.962595993153907D4*xyzi(i+1,j+7,k)+9.73269191
     $        3009597D3*xyzi(i+1,j+5,k))+angi
            angi=yrk(153)*(3.452141788424563D4*xyzi(i+11,j+1,k)+1.035642
     $        536527369D5*xyzi(i+9,j+3,k)-9.605959759094437D4*xyzi(i+9,j
     $        +1,k)+6.904283576849127D4*xyzi(i+7,j+5,k)-1.92119195181888
     $        7D5*xyzi(i+7,j+3,k)+9.605959759094437D4*xyzi(i+7,j+1,k)-6.
     $        904283576849127D4*xyzi(i+5,j+7,k)+9.605959759094437D4*xyzi
     $        (i+5,j+3,k)-4.044614635408184D4*xyzi(i+5,j+1,k)-1.03564253
     $        6527369D5*xyzi(i+3,j+9,k)+1.921191951818887D5*xyzi(i+3,j+7
     $        ,k)-9.605959759094437D4*xyzi(i+3,j+5,k)+5.947962699129683D
     $        3*xyzi(i+3,j+1,k)-3.452141788424563D4*xyzi(i+1,j+11,k)+9.6
     $        05959759094437D4*xyzi(i+1,j+9,k)-9.605959759094437D4*xyzi(
     $        i+1,j+7,k)+4.044614635408184D4*xyzi(i+1,j+5,k)-5.947962699
     $        129683D3*xyzi(i+1,j+3,k))+angi
            angi=yrk(155)*(-2.818661967126381D4*xyzi(i+11,j+1,k)-1.40933
     $        0983563191D5*xyzi(i+9,j+3,k)+8.578536421688987D4*xyzi(i+9,
     $        j+1,k)-2.818661967126381D5*xyzi(i+7,j+5,k)+3.4314145686755
     $        95D5*xyzi(i+7,j+3,k)-9.804041624787413D4*xyzi(i+7,j+1,k)-2
     $        .818661967126381D5*xyzi(i+5,j+7,k)+5.147121853013392D5*xyz
     $        i(i+5,j+5,k)-2.941212487436224D5*xyzi(i+5,j+3,k)+5.1600219
     $        07782849D4*xyzi(i+5,j+1,k)-1.409330983563191D5*xyzi(i+3,j+
     $        9,k)+3.431414568675595D5*xyzi(i+3,j+7,k)-2.941212487436224
     $        D5*xyzi(i+3,j+5,k)+1.03200438155657D5*xyzi(i+3,j+3,k)-1.21
     $        4122801831259D4*xyzi(i+3,j+1,k)-2.818661967126381D4*xyzi(i
     $        +1,j+11,k)+8.578536421688987D4*xyzi(i+1,j+9,k)-9.804041624
     $        787413D4*xyzi(i+1,j+7,k)+5.160021907782849D4*xyzi(i+1,j+5,
     $        k)-1.214122801831259D4*xyzi(i+1,j+3,k)+9.712982414650069D2
     $        *xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(158)*(1.589939778654773D4*xyzi(i+11,j,k+1)+7.949698
     $        893273864D4*xyzi(i+9,j+2,k+1)-4.147668987795059D4*xyzi(i+9
     $        ,j,k+1)+1.589939778654773D5*xyzi(i+7,j+4,k+1)-1.6590675951
     $        18024D5*xyzi(i+7,j+2,k+1)+3.9501609407572D4*xyzi(i+7,j,k+1
     $        )+1.589939778654773D5*xyzi(i+5,j+6,k+1)-2.488601392677036D
     $        5*xyzi(i+5,j+4,k+1)+1.18504828222716D5*xyzi(i+5,j+2,k+1)-1
     $        .663225659266189D4*xyzi(i+5,j,k+1)+7.949698893273864D4*xyz
     $        i(i+3,j+8,k+1)-1.659067595118024D5*xyzi(i+3,j+6,k+1)+1.185
     $        04828222716D5*xyzi(i+3,j+4,k+1)-3.326451318532379D4*xyzi(i
     $        +3,j+2,k+1)+2.935104104587393D3*xyzi(i+3,j,k+1)+1.58993977
     $        8654773D4*xyzi(i+1,j+10,k+1)-4.147668987795059D4*xyzi(i+1,
     $        j+8,k+1)+3.9501609407572D4*xyzi(i+1,j+6,k+1)-1.66322565926
     $        6189D4*xyzi(i+1,j+4,k+1)+2.935104104587393D3*xyzi(i+1,j+2,
     $        k+1)-1.565388855779943D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(160)*(-1.150713929474854D4*xyzi(i+11,j,k+1)-1.15071
     $        3929474854D4*xyzi(i+9,j+2,k+1)+2.801738263069211D4*xyzi(i+
     $        9,j,k+1)+6.904283576849127D4*xyzi(i+7,j+4,k+1)-2.401489939
     $        773609D4*xyzi(i+7,j,k+1)+1.610999501264796D5*xyzi(i+5,j+6,
     $        k+1)-1.681042957841527D5*xyzi(i+5,j+4,k+1)+2.4014899397736
     $        09D4*xyzi(i+5,j+2,k+1)+8.426280490433717D3*xyzi(i+5,j,k+1)
     $        +1.26578532242234D5*xyzi(i+3,j+8,k+1)-2.241390610455369D5*
     $        xyzi(i+3,j+6,k+1)+1.200744969886805D5*xyzi(i+3,j+4,k+1)-1.
     $        685256098086743D4*xyzi(i+3,j+2,k+1)-9.913271165216138D2*xy
     $        zi(i+3,j,k+1)+3.452141788424563D4*xyzi(i+1,j+10,k+1)-8.405
     $        214789207633D4*xyzi(i+1,j+8,k+1)+7.204469819320828D4*xyzi(
     $        i+1,j+6,k+1)-2.527884147130115D4*xyzi(i+1,j+4,k+1)+2.97398
     $        1349564841D3*xyzi(i+1,j+2,k+1))+angi
            angi=yrk(162)*(5.920374324261427D3*xyzi(i+11,j,k+1)-4.144262
     $        026982999D4*xyzi(i+9,j+2,k+1)-1.235556380715428D4*xyzi(i+9
     $        ,j,k+1)-1.302482351337514D5*xyzi(i+7,j+4,k+1)+9.8844510457
     $        23426D4*xyzi(i+7,j+2,k+1)+8.237042538102855D3*xyzi(i+7,j,k
     $        +1)-8.288524053965998D4*xyzi(i+5,j+6,k+1)+1.7297789330016D
     $        5*xyzi(i+5,j+4,k+1)-7.41333828429257D4*xyzi(i+5,j+2,k+1)-1
     $        .73411421854797D3*xyzi(i+5,j,k+1)+2.960187162130714D4*xyzi
     $        (i+3,j+8,k+1)-4.118521269051428D4*xyzi(i+3,j+4,k+1)+1.7341
     $        14218547969D4*xyzi(i+3,j+2,k+1)+2.960187162130714D4*xyzi(i
     $        +1,j+10,k+1)-6.177781903577141D4*xyzi(i+1,j+8,k+1)+4.11852
     $        1269051428D4*xyzi(i+1,j+6,k+1)-8.670571092739847D3*xyzi(i+
     $        1,j+4,k+1))+angi
            angi=yrk(164)*(-2.074726043955437D3*xyzi(i+11,j,k+1)+3.94197
     $        948351533D4*xyzi(i+9,j+2,k+1)+3.247397286191118D3*xyzi(i+9
     $        ,j,k+1)+1.244835626373262D4*xyzi(i+7,j+4,k+1)-6.4947945723
     $        82236D4*xyzi(i+7,j+2,k+1)-1.237103728072807D3*xyzi(i+7,j,k
     $        +1)-8.713849384612834D4*xyzi(i+5,j+6,k+1)+4.54635620066756
     $        5D4*xyzi(i+5,j+4,k+1)+2.597917828952894D4*xyzi(i+5,j+2,k+1
     $        )-4.356924692306417D4*xyzi(i+3,j+8,k+1)+9.092712401335131D
     $        4*xyzi(i+3,j+6,k+1)-4.329863048254824D4*xyzi(i+3,j+4,k+1)+
     $        1.452308230768806D4*xyzi(i+1,j+10,k+1)-2.273178100333783D4
     $        *xyzi(i+1,j+8,k+1)+8.659726096509648D3*xyzi(i+1,j+6,k+1))+
     $        angi
            angi=yrk(166)*(4.527423401296222D2*xyzi(i+11,j,k+1)-1.584598
     $        190453678D4*xyzi(i+9,j+2,k+1)-3.936889914170628D2*xyzi(i+9
     $        ,j,k+1)+4.0746810611666D4*xyzi(i+7,j+4,k+1)+1.417280369101
     $        426D4*xyzi(i+7,j+2,k+1)+1.901517828544413D4*xyzi(i+5,j+6,k
     $        +1)-4.960481291854991D4*xyzi(i+5,j+4,k+1)-3.39556755097216
     $        6D4*xyzi(i+3,j+8,k+1)+3.306987527903327D4*xyzi(i+3,j+6,k+1
     $        )+4.0746810611666D3*xyzi(i+1,j+10,k+1)-3.543200922753565D3
     $        *xyzi(i+1,j+8,k+1))+angi
            angi=yrk(168)*(-4.930051750476566D1*xyzi(i+11,j,k+1)+2.71152
     $        8462762111D3*xyzi(i+9,j+2,k+1)-1.626917077657267D4*xyzi(i+
     $        7,j+4,k+1)+2.277683908720174D4*xyzi(i+5,j+6,k+1)-8.1345853
     $        88286334D3*xyzi(i+3,j+8,k+1)+5.423056925524223D2*xyzi(i+1,
     $        j+10,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(146)*(-5.423056925524223D2*xyzi(i+10,j+1,k+1)+8.134
     $        585388286334D3*xyzi(i+8,j+3,k+1)-2.277683908720174D4*xyzi(
     $        i+6,j+5,k+1)+1.626917077657267D4*xyzi(i+4,j+7,k+1)-2.71152
     $        8462762111D3*xyzi(i+2,j+9,k+1)+4.930051750476566D1*xyzi(i,
     $        j+11,k+1))+angi
            angi=yrk(148)*(4.0746810611666D3*xyzi(i+10,j+1,k+1)-3.395567
     $        550972166D4*xyzi(i+8,j+3,k+1)-3.543200922753565D3*xyzi(i+8
     $        ,j+1,k+1)+1.901517828544413D4*xyzi(i+6,j+5,k+1)+3.30698752
     $        7903327D4*xyzi(i+6,j+3,k+1)+4.0746810611666D4*xyzi(i+4,j+7
     $        ,k+1)-4.960481291854991D4*xyzi(i+4,j+5,k+1)-1.584598190453
     $        678D4*xyzi(i+2,j+9,k+1)+1.417280369101426D4*xyzi(i+2,j+7,k
     $        +1)+4.527423401296222D2*xyzi(i,j+11,k+1)-3.936889914170628
     $        D2*xyzi(i,j+9,k+1))+angi
            angi=yrk(150)*(-1.452308230768806D4*xyzi(i+10,j+1,k+1)+4.356
     $        924692306417D4*xyzi(i+8,j+3,k+1)+2.273178100333783D4*xyzi(
     $        i+8,j+1,k+1)+8.713849384612834D4*xyzi(i+6,j+5,k+1)-9.09271
     $        2401335131D4*xyzi(i+6,j+3,k+1)-8.659726096509648D3*xyzi(i+
     $        6,j+1,k+1)-1.244835626373262D4*xyzi(i+4,j+7,k+1)-4.5463562
     $        00667565D4*xyzi(i+4,j+5,k+1)+4.329863048254824D4*xyzi(i+4,
     $        j+3,k+1)-3.94197948351533D4*xyzi(i+2,j+9,k+1)+6.4947945723
     $        82236D4*xyzi(i+2,j+7,k+1)-2.597917828952894D4*xyzi(i+2,j+5
     $        ,k+1)+2.074726043955437D3*xyzi(i,j+11,k+1)-3.2473972861911
     $        18D3*xyzi(i,j+9,k+1)+1.237103728072807D3*xyzi(i,j+7,k+1))+
     $        angi
            angi=yrk(152)*(2.960187162130714D4*xyzi(i+10,j+1,k+1)+2.9601
     $        87162130714D4*xyzi(i+8,j+3,k+1)-6.177781903577141D4*xyzi(i
     $        +8,j+1,k+1)-8.288524053965998D4*xyzi(i+6,j+5,k+1)+4.118521
     $        269051428D4*xyzi(i+6,j+1,k+1)-1.302482351337514D5*xyzi(i+4
     $        ,j+7,k+1)+1.7297789330016D5*xyzi(i+4,j+5,k+1)-4.1185212690
     $        51428D4*xyzi(i+4,j+3,k+1)-8.670571092739847D3*xyzi(i+4,j+1
     $        ,k+1)-4.144262026982999D4*xyzi(i+2,j+9,k+1)+9.884451045723
     $        426D4*xyzi(i+2,j+7,k+1)-7.41333828429257D4*xyzi(i+2,j+5,k+
     $        1)+1.734114218547969D4*xyzi(i+2,j+3,k+1)+5.920374324261427
     $        D3*xyzi(i,j+11,k+1)-1.235556380715428D4*xyzi(i,j+9,k+1)+8.
     $        237042538102855D3*xyzi(i,j+7,k+1)-1.73411421854797D3*xyzi(
     $        i,j+5,k+1))+angi
            angi=yrk(154)*(-3.452141788424563D4*xyzi(i+10,j+1,k+1)-1.265
     $        78532242234D5*xyzi(i+8,j+3,k+1)+8.405214789207633D4*xyzi(i
     $        +8,j+1,k+1)-1.610999501264796D5*xyzi(i+6,j+5,k+1)+2.241390
     $        610455369D5*xyzi(i+6,j+3,k+1)-7.204469819320828D4*xyzi(i+6
     $        ,j+1,k+1)-6.904283576849127D4*xyzi(i+4,j+7,k+1)+1.68104295
     $        7841527D5*xyzi(i+4,j+5,k+1)-1.200744969886805D5*xyzi(i+4,j
     $        +3,k+1)+2.527884147130115D4*xyzi(i+4,j+1,k+1)+1.1507139294
     $        74854D4*xyzi(i+2,j+9,k+1)-2.401489939773609D4*xyzi(i+2,j+5
     $        ,k+1)+1.685256098086743D4*xyzi(i+2,j+3,k+1)-2.973981349564
     $        841D3*xyzi(i+2,j+1,k+1)+1.150713929474854D4*xyzi(i,j+11,k+
     $        1)-2.801738263069211D4*xyzi(i,j+9,k+1)+2.401489939773609D4
     $        *xyzi(i,j+7,k+1)-8.426280490433717D3*xyzi(i,j+5,k+1)+9.913
     $        271165216138D2*xyzi(i,j+3,k+1))+angi
            angi=yrk(156)*(1.589939778654773D4*xyzi(i+10,j+1,k+1)+7.9496
     $        98893273864D4*xyzi(i+8,j+3,k+1)-4.147668987795059D4*xyzi(i
     $        +8,j+1,k+1)+1.589939778654773D5*xyzi(i+6,j+5,k+1)-1.659067
     $        595118024D5*xyzi(i+6,j+3,k+1)+3.9501609407572D4*xyzi(i+6,j
     $        +1,k+1)+1.589939778654773D5*xyzi(i+4,j+7,k+1)-2.4886013926
     $        77036D5*xyzi(i+4,j+5,k+1)+1.18504828222716D5*xyzi(i+4,j+3,
     $        k+1)-1.663225659266189D4*xyzi(i+4,j+1,k+1)+7.9496988932738
     $        64D4*xyzi(i+2,j+9,k+1)-1.659067595118024D5*xyzi(i+2,j+7,k+
     $        1)+1.18504828222716D5*xyzi(i+2,j+5,k+1)-3.326451318532379D
     $        4*xyzi(i+2,j+3,k+1)+2.935104104587393D3*xyzi(i+2,j+1,k+1)+
     $        1.589939778654773D4*xyzi(i,j+11,k+1)-4.147668987795059D4*x
     $        yzi(i,j+9,k+1)+3.9501609407572D4*xyzi(i,j+7,k+1)-1.6632256
     $        59266189D4*xyzi(i,j+5,k+1)+2.935104104587393D3*xyzi(i,j+3,
     $        k+1)-1.565388855779943D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=14
c...
        if(lone.eq.14)then
          angi=ZERO
          if((ieven.and.jeven.and.keven))then
            angi=yrk(211)*(4.674208812108592D4*xyzi(i,j,k+14)-1.57538148
     $        8525489D5*xyzi(i,j,k+12)+2.079503564853645D5*xyzi(i,j,k+10
     $        )-1.356197977078464D5*xyzi(i,j,k+8)+4.52065992359488D4*xyz
     $        i(i,j,k+6)-7.137884089886653D3*xyzi(i,j,k+4)+4.19875534699
     $        2149D2*xyzi(i,j,k+2)-3.998814616182999D0*xyzi(i,j,k))+angi
            angi=yrk(213)*(5.756429376136187D4*xyzi(i+14,j,k)+2.87821468
     $        8068094D5*xyzi(i+12,j+2,k)-2.046730444848422D5*xyzi(i+12,j
     $        ,k)+5.180786438522568D5*xyzi(i+10,j+4,k)-8.186921779393688
     $        D5*xyzi(i+10,j+2,k)+2.865422622787791D5*xyzi(i+10,j,k)+2.8
     $        78214688068094D5*xyzi(i+8,j+6,k)-1.023365222424211D6*xyzi(
     $        i+8,j+4,k)+8.596267868363373D5*xyzi(i+8,j+2,k)-1.993337476
     $        721942D5*xyzi(i+8,j,k)-2.878214688068094D5*xyzi(i+6,j+8,k)
     $        +5.730845245575582D5*xyzi(i+6,j+4,k)-3.986674953443883D5*x
     $        yzi(i+6,j+2,k)+7.119062416864077D4*xyzi(i+6,j,k)-5.1807864
     $        38522568D5*xyzi(i+4,j+10,k)+1.023365222424211D6*xyzi(i+4,j
     $        +8,k)-5.730845245575582D5*xyzi(i+4,j+6,k)+7.11906241686407
     $        7D4*xyzi(i+4,j+2,k)-1.19899998599816D4*xyzi(i+4,j,k)-2.878
     $        214688068094D5*xyzi(i+2,j+12,k)+8.186921779393688D5*xyzi(i
     $        +2,j+10,k)-8.596267868363373D5*xyzi(i+2,j+8,k)+3.986674953
     $        443883D5*xyzi(i+2,j+6,k)-7.119062416864077D4*xyzi(i+2,j+4,
     $        k)+7.052941094106825D2*xyzi(i+2,j,k)-5.756429376136187D4*x
     $        yzi(i,j+14,k)+2.046730444848422D5*xyzi(i,j+12,k)-2.8654226
     $        22787791D5*xyzi(i,j+10,k)+1.993337476721942D5*xyzi(i,j+8,k
     $        )-7.119062416864077D4*xyzi(i,j+6,k)+1.19899998599816D4*xyz
     $        i(i,j+4,k)-7.052941094106825D2*xyzi(i,j+2,k))+angi
            angi=yrk(215)*(-3.780762817453435D4*xyzi(i+14,j,k)+3.7807628
     $        17453435D4*xyzi(i+12,j+2,k)+1.260254272484478D5*xyzi(i+12,
     $        j,k)+7.183449353161527D5*xyzi(i+10,j+4,k)-2.52050854496895
     $        7D5*xyzi(i+10,j+2,k)-1.613125468780132D5*xyzi(i+10,j,k)+1.
     $        701343267854046D6*xyzi(i+8,j+6,k)-2.142432263223613D6*xyzi
     $        (i+8,j+4,k)+4.839376406340397D5*xyzi(i+8,j+2,k)+9.81902459
     $        2574719D4*xyzi(i+8,j,k)+1.701343267854046D6*xyzi(i+6,j+8,k
     $        )-3.52871196295654D6*xyzi(i+6,j+6,k)+2.258375656292185D6*x
     $        yzi(i+6,j+4,k)-3.927609837029888D5*xyzi(i+6,j+2,k)-2.80543
     $        5597878491D4*xyzi(i+6,j,k)+7.183449353161527D5*xyzi(i+4,j+
     $        10,k)-2.142432263223613D6*xyzi(i+4,j+8,k)+2.25837565629218
     $        5D6*xyzi(i+4,j+6,k)-9.819024592574719D5*xyzi(i+4,j+4,k)+1.
     $        402717798939246D5*xyzi(i+4,j+2,k)+2.953090103029991D3*xyzi
     $        (i+4,j,k)+3.780762817453435D4*xyzi(i+2,j+12,k)-2.520508544
     $        968957D5*xyzi(i+2,j+10,k)+4.839376406340397D5*xyzi(i+2,j+8
     $        ,k)-3.927609837029888D5*xyzi(i+2,j+6,k)+1.402717798939246D
     $        5*xyzi(i+2,j+4,k)-1.771854061817994D4*xyzi(i+2,j+2,k)-3.78
     $        0762817453435D4*xyzi(i,j+14,k)+1.260254272484478D5*xyzi(i,
     $        j+12,k)-1.613125468780132D5*xyzi(i,j+10,k)+9.8190245925747
     $        19D4*xyzi(i,j+8,k)-2.805435597878491D4*xyzi(i,j+6,k)+2.953
     $        090103029991D3*xyzi(i,j+4,k))+angi
            angi=yrk(217)*(1.839962151616926D4*xyzi(i+14,j,k)-2.02395836
     $        6778619D5*xyzi(i+12,j+2,k)-5.451739708494596D4*xyzi(i+12,j
     $        ,k)-7.175852391306012D5*xyzi(i+10,j+4,k)+6.542087650193515
     $        D5*xyzi(i+10,j+2,k)+5.887878885174163D4*xyzi(i+10,j,k)-4.9
     $        678978093657D5*xyzi(i+8,j+6,k)+1.471969721293541D6*xyzi(i+
     $        8,j+4,k)-7.654242550726412D5*xyzi(i+8,j+2,k)-2.73061049747
     $        2076D4*xyzi(i+8,j,k)+4.9678978093657D5*xyzi(i+6,j+8,k)-8.2
     $        43030439243829D5*xyzi(i+6,j+4,k)+3.822854696460906D5*xyzi(
     $        i+6,j+2,k)+4.551017495786793D3*xyzi(i+6,j,k)+7.17585239130
     $        6012D5*xyzi(i+4,j+10,k)-1.471969721293541D6*xyzi(i+4,j+8,k
     $        )+8.243030439243829D5*xyzi(i+4,j+6,k)-6.826526243680189D4*
     $        xyzi(i+4,j+2,k)+2.023958366778619D5*xyzi(i+2,j+12,k)-6.542
     $        087650193515D5*xyzi(i+2,j+10,k)+7.654242550726412D5*xyzi(i
     $        +2,j+8,k)-3.822854696460906D5*xyzi(i+2,j+6,k)+6.8265262436
     $        80189D4*xyzi(i+2,j+4,k)-1.839962151616926D4*xyzi(i,j+14,k)
     $        +5.451739708494596D4*xyzi(i,j+12,k)-5.887878885174163D4*xy
     $        zi(i,j+10,k)+2.730610497472076D4*xyzi(i,j+8,k)-4.551017495
     $        786793D3*xyzi(i,j+6,k))+angi
            angi=yrk(219)*(-6.405925968013536D3*xyzi(i+14,j,k)+1.6014814
     $        92003384D5*xyzi(i+12,j+2,k)+1.565893014403309D4*xyzi(i+12,
     $        j,k)+7.04651856481489D4*xyzi(i+10,j+4,k)-4.071321837448603
     $        D5*xyzi(i+10,j+2,k)-1.252714411522647D4*xyzi(i+10,j,k)-6.3
     $        41866708333401D5*xyzi(i+8,j+6,k)+2.348839521604963D5*xyzi(
     $        i+8,j+4,k)+3.382328911111147D5*xyzi(i+8,j+2,k)+3.267950638
     $        754731D3*xyzi(i+8,j,k)-6.341866708333401D5*xyzi(i+6,j+8,k)
     $        +1.315350132098779D6*xyzi(i+6,j+6,k)-5.261400528395117D5*x
     $        yzi(i+6,j+4,k)-9.150261788513248D4*xyzi(i+6,j+2,k)+7.04651
     $        856481489D4*xyzi(i+4,j+10,k)+2.348839521604963D5*xyzi(i+4,
     $        j+8,k)-5.261400528395117D5*xyzi(i+4,j+6,k)+2.2875654471283
     $        12D5*xyzi(i+4,j+4,k)+1.601481492003384D5*xyzi(i+2,j+12,k)-
     $        4.071321837448603D5*xyzi(i+2,j+10,k)+3.382328911111147D5*x
     $        yzi(i+2,j+8,k)-9.150261788513248D4*xyzi(i+2,j+6,k)-6.40592
     $        5968013536D3*xyzi(i,j+14,k)+1.565893014403309D4*xyzi(i,j+1
     $        2,k)-1.252714411522647D4*xyzi(i,j+10,k)+3.267950638754731D
     $        3*xyzi(i,j+8,k))+angi
            angi=yrk(221)*(1.493389191601027D3*xyzi(i+14,j,k)-6.42157352
     $        3884417D4*xyzi(i+12,j+2,k)-2.654914118401826D3*xyzi(i+12,j
     $        ,k)+1.807000921837243D5*xyzi(i+10,j+4,k)+1.168162212096803
     $        D5*xyzi(i+10,j+2,k)+1.168162212096803D3*xyzi(i+10,j,k)+2.4
     $        64092166141695D5*xyzi(i+8,j+6,k)-4.380608295363013D5*xyzi(
     $        i+8,j+4,k)-5.256729954435616D4*xyzi(i+8,j+2,k)-2.464092166
     $        141695D5*xyzi(i+6,j+8,k)+2.453140645403287D5*xyzi(i+6,j+4,
     $        k)-1.807000921837243D5*xyzi(i+4,j+10,k)+4.380608295363013D
     $        5*xyzi(i+4,j+8,k)-2.453140645403287D5*xyzi(i+4,j+6,k)+6.42
     $        1573523884417D4*xyzi(i+2,j+12,k)-1.168162212096803D5*xyzi(
     $        i+2,j+10,k)+5.256729954435616D4*xyzi(i+2,j+8,k)-1.49338919
     $        1601027D3*xyzi(i,j+14,k)+2.654914118401826D3*xyzi(i,j+12,k
     $        )-1.168162212096803D3*xyzi(i,j+10,k))+angi
            angi=yrk(223)*(-2.029116341627528D2*xyzi(i+14,j,k)+1.3189256
     $        22057893D4*xyzi(i+12,j+2,k)+1.953963884530212D2*xyzi(i+12,
     $        j,k)-8.704909105582095D4*xyzi(i+10,j+4,k)-1.28961616378994
     $        D4*xyzi(i+10,j+2,k)+8.704909105582095D4*xyzi(i+8,j+6,k)+9.
     $        672121228424549D4*xyzi(i+8,j+4,k)+8.704909105582095D4*xyzi
     $        (i+6,j+8,k)-1.805462629305916D5*xyzi(i+6,j+6,k)-8.70490910
     $        5582095D4*xyzi(i+4,j+10,k)+9.672121228424549D4*xyzi(i+4,j+
     $        8,k)+1.318925622057893D4*xyzi(i+2,j+12,k)-1.28961616378994
     $        D4*xyzi(i+2,j+10,k)-2.029116341627528D2*xyzi(i,j+14,k)+1.9
     $        53963884530212D2*xyzi(i,j+12,k))+angi
            angi=yrk(225)*(1.04366482991984D1*xyzi(i+14,j,k)-9.497349952
     $        270546D2*xyzi(i+12,j+2,k)+1.04470849474976D4*xyzi(i+10,j+4
     $        ,k)-3.13412548424928D4*xyzi(i+8,j+6,k)+3.13412548424928D4*
     $        xyzi(i+6,j+8,k)-1.04470849474976D4*xyzi(i+4,j+10,k)+9.4973
     $        49952270546D2*xyzi(i+2,j+12,k)-1.04366482991984D1*xyzi(i,j
     $        +14,k))+angi
          endif
          if((.not.ieven.and..not.jeven.and.keven))then
            angi=yrk(197)*(1.461130761887776D2*xyzi(i+13,j+1,k)-3.798939
     $        980908219D3*xyzi(i+11,j+3,k)+2.08941698949952D4*xyzi(i+9,j
     $        +5,k)-3.581857696284892D4*xyzi(i+7,j+7,k)+2.08941698949952
     $        D4*xyzi(i+5,j+9,k)-3.798939980908219D3*xyzi(i+3,j+11,k)+1.
     $        461130761887776D2*xyzi(i+1,j+13,k))+angi
            angi=yrk(199)*(-2.434939609953033D3*xyzi(i+13,j+1,k)+4.22056
     $        1990585258D4*xyzi(i+11,j+3,k)+2.344756661436254D3*xyzi(i+1
     $        1,j+1,k)-1.160654547410946D5*xyzi(i+9,j+5,k)-4.29872054596
     $        6466D4*xyzi(i+9,j+3,k)+1.547539396547928D5*xyzi(i+7,j+5,k)
     $        +1.160654547410946D5*xyzi(i+5,j+9,k)-1.547539396547928D5*x
     $        yzi(i+5,j+7,k)-4.220561990585258D4*xyzi(i+3,j+11,k)+4.2987
     $        20545966466D4*xyzi(i+3,j+9,k)+2.434939609953033D3*xyzi(i+1
     $        ,j+13,k)-2.344756661436254D3*xyzi(i+1,j+11,k))+angi
            angi=yrk(201)*(1.493389191601027D4*xyzi(i+13,j+1,k)-1.493389
     $        191601027D5*xyzi(i+11,j+3,k)-2.654914118401826D4*xyzi(i+11
     $        ,j+1,k)+3.28545622152226D4*xyzi(i+9,j+5,k)+2.9204055302420
     $        09D5*xyzi(i+9,j+3,k)+1.168162212096803D4*xyzi(i+9,j+1,k)+3
     $        .942547465826712D5*xyzi(i+7,j+7,k)-3.50448663629041D5*xyzi
     $        (i+7,j+5,k)-1.401794654516164D5*xyzi(i+7,j+3,k)+3.28545622
     $        152226D4*xyzi(i+5,j+9,k)-3.50448663629041D5*xyzi(i+5,j+7,k
     $        )+2.943768774483945D5*xyzi(i+5,j+5,k)-1.493389191601027D5*
     $        xyzi(i+3,j+11,k)+2.920405530242009D5*xyzi(i+3,j+9,k)-1.401
     $        794654516164D5*xyzi(i+3,j+7,k)+1.493389191601027D4*xyzi(i+
     $        1,j+13,k)-2.654914118401826D4*xyzi(i+1,j+11,k)+1.168162212
     $        096803D4*xyzi(i+1,j+9,k))+angi
            angi=yrk(203)*(-5.124740774410829D4*xyzi(i+13,j+1,k)+2.04989
     $        6309764331D5*xyzi(i+11,j+3,k)+1.252714411522647D5*xyzi(i+1
     $        1,j+1,k)+5.637214851851912D5*xyzi(i+9,j+5,k)-6.26357205761
     $        3235D5*xyzi(i+9,j+3,k)-1.002171529218118D5*xyzi(i+9,j+1,k)
     $        -7.516286469135882D5*xyzi(i+7,j+5,k)+6.013029175308706D5*x
     $        yzi(i+7,j+3,k)+2.614360511003785D4*xyzi(i+7,j+1,k)-5.63721
     $        4851851912D5*xyzi(i+5,j+9,k)+7.516286469135882D5*xyzi(i+5,
     $        j+7,k)-1.83005235770265D5*xyzi(i+5,j+3,k)-2.04989630976433
     $        1D5*xyzi(i+3,j+11,k)+6.263572057613235D5*xyzi(i+3,j+9,k)-6
     $        .013029175308706D5*xyzi(i+3,j+7,k)+1.83005235770265D5*xyzi
     $        (i+3,j+5,k)+5.124740774410829D4*xyzi(i+1,j+13,k)-1.2527144
     $        11522647D5*xyzi(i+1,j+11,k)+1.002171529218118D5*xyzi(i+1,j
     $        +9,k)-2.614360511003785D4*xyzi(i+1,j+7,k))+angi
            angi=yrk(205)*(1.103977290970156D5*xyzi(i+13,j+1,k)+7.359848
     $        606467704D4*xyzi(i+11,j+3,k)-3.271043825096757D5*xyzi(i+11
     $        ,j+1,k)-6.991856176144319D5*xyzi(i+9,j+5,k)+1.090347941698
     $        919D5*xyzi(i+9,j+3,k)+3.532727331104498D5*xyzi(i+9,j+1,k)-
     $        1.324772749164187D6*xyzi(i+7,j+7,k)+1.962626295058055D6*xy
     $        zi(i+7,j+5,k)-4.710303108139331D5*xyzi(i+7,j+3,k)-1.638366
     $        298483245D5*xyzi(i+7,j+1,k)-6.991856176144319D5*xyzi(i+5,j
     $        +9,k)+1.962626295058055D6*xyzi(i+5,j+7,k)-1.64860608784876
     $        6D6*xyzi(i+5,j+5,k)+3.822854696460906D5*xyzi(i+5,j+3,k)+2.
     $        730610497472076D4*xyzi(i+5,j+1,k)+7.359848606467704D4*xyzi
     $        (i+3,j+11,k)+1.090347941698919D5*xyzi(i+3,j+9,k)-4.7103031
     $        08139331D5*xyzi(i+3,j+7,k)+3.822854696460906D5*xyzi(i+3,j+
     $        5,k)-9.102034991573586D4*xyzi(i+3,j+3,k)+1.103977290970156
     $        D5*xyzi(i+1,j+13,k)-3.271043825096757D5*xyzi(i+1,j+11,k)+3
     $        .532727331104498D5*xyzi(i+1,j+9,k)-1.638366298483245D5*xyz
     $        i(i+1,j+7,k)+2.730610497472076D4*xyzi(i+1,j+5,k))+angi
            angi=yrk(207)*(-1.512305126981374D5*xyzi(i+13,j+1,k)-6.04922
     $        0507925497D5*xyzi(i+11,j+3,k)+5.041017089937914D5*xyzi(i+1
     $        1,j+1,k)-7.561525634906871D5*xyzi(i+9,j+5,k)+1.51230512698
     $        1374D6*xyzi(i+9,j+3,k)-6.45250187512053D5*xyzi(i+9,j+1,k)+
     $        1.008203417987583D6*xyzi(i+7,j+5,k)-1.290500375024106D6*xy
     $        zi(i+7,j+3,k)+3.927609837029888D5*xyzi(i+7,j+1,k)+7.561525
     $        634906871D5*xyzi(i+5,j+9,k)-1.008203417987583D6*xyzi(i+5,j
     $        +7,k)+3.927609837029888D5*xyzi(i+5,j+3,k)-1.12217423915139
     $        6D5*xyzi(i+5,j+1,k)+6.049220507925497D5*xyzi(i+3,j+11,k)-1
     $        .512305126981374D6*xyzi(i+3,j+9,k)+1.290500375024106D6*xyz
     $        i(i+3,j+7,k)-3.927609837029888D5*xyzi(i+3,j+5,k)+1.1812360
     $        41211996D4*xyzi(i+3,j+1,k)+1.512305126981374D5*xyzi(i+1,j+
     $        13,k)-5.041017089937914D5*xyzi(i+1,j+11,k)+6.4525018751205
     $        3D5*xyzi(i+1,j+9,k)-3.927609837029888D5*xyzi(i+1,j+7,k)+1.
     $        122174239151396D5*xyzi(i+1,j+5,k)-1.181236041211996D4*xyzi
     $        (i+1,j+3,k))+angi
            angi=yrk(209)*(1.151285875227237D5*xyzi(i+13,j+1,k)+6.907715
     $        251363424D5*xyzi(i+11,j+3,k)-4.093460889696844D5*xyzi(i+11
     $        ,j+1,k)+1.726928812840856D6*xyzi(i+9,j+5,k)-2.046730444848
     $        422D6*xyzi(i+9,j+3,k)+5.730845245575582D5*xyzi(i+9,j+1,k)+
     $        2.302571750454475D6*xyzi(i+7,j+7,k)-4.093460889696844D6*xy
     $        zi(i+7,j+5,k)+2.292338098230233D6*xyzi(i+7,j+3,k)-3.986674
     $        953443883D5*xyzi(i+7,j+1,k)+1.726928812840856D6*xyzi(i+5,j
     $        +9,k)-4.093460889696844D6*xyzi(i+5,j+7,k)+3.43850714734534
     $        9D6*xyzi(i+5,j+5,k)-1.196002486033165D6*xyzi(i+5,j+3,k)+1.
     $        423812483372815D5*xyzi(i+5,j+1,k)+6.907715251363424D5*xyzi
     $        (i+3,j+11,k)-2.046730444848422D6*xyzi(i+3,j+9,k)+2.2923380
     $        98230233D6*xyzi(i+3,j+7,k)-1.196002486033165D6*xyzi(i+3,j+
     $        5,k)+2.847624966745631D5*xyzi(i+3,j+3,k)-2.397999971996321
     $        D4*xyzi(i+3,j+1,k)+1.151285875227237D5*xyzi(i+1,j+13,k)-4.
     $        093460889696844D5*xyzi(i+1,j+11,k)+5.730845245575582D5*xyz
     $        i(i+1,j+9,k)-3.986674953443883D5*xyzi(i+1,j+7,k)+1.4238124
     $        83372815D5*xyzi(i+1,j+5,k)-2.397999971996321D4*xyzi(i+1,j+
     $        3,k)+1.410588218821365D3*xyzi(i+1,j+1,k))+angi
          endif
          if((.not.ieven.and.jeven.and..not.keven))then
            angi=yrk(212)*(-6.386185008998833D4*xyzi(i+13,j,k+1)-3.83171
     $        10053993D5*xyzi(i+11,j+2,k+1)+1.986813113910748D5*xyzi(i+1
     $        1,j,k+1)-9.579277513498249D5*xyzi(i+9,j+4,k+1)+9.934065569
     $        553739D5*xyzi(i+9,j+2,k+1)-2.384175736692897D5*xyzi(i+9,j,
     $        k+1)-1.277237001799767D6*xyzi(i+7,j+6,k+1)+1.9868131139107
     $        48D6*xyzi(i+7,j+4,k+1)-9.53670294677159D5*xyzi(i+7,j+2,k+1
     $        )+1.382130861850955D5*xyzi(i+7,j,k+1)-9.579277513498249D5*
     $        xyzi(i+5,j+8,k+1)+1.986813113910748D6*xyzi(i+5,j+6,k+1)-1.
     $        430505442015738D6*xyzi(i+5,j+4,k+1)+4.146392585552865D5*xy
     $        zi(i+5,j+2,k+1)-3.948945319574157D4*xyzi(i+5,j,k+1)-3.8317
     $        110053993D5*xyzi(i+3,j+10,k+1)+9.934065569553739D5*xyzi(i+
     $        3,j+8,k+1)-9.53670294677159D5*xyzi(i+3,j+6,k+1)+4.14639258
     $        5552865D5*xyzi(i+3,j+4,k+1)-7.897890639148315D4*xyzi(i+3,j
     $        +2,k+1)+4.988141456304199D3*xyzi(i+3,j,k+1)-6.386185008998
     $        833D4*xyzi(i+1,j+12,k+1)+1.986813113910748D5*xyzi(i+1,j+10
     $        ,k+1)-2.384175736692897D5*xyzi(i+1,j+8,k+1)+1.382130861850
     $        955D5*xyzi(i+1,j+6,k+1)-3.948945319574157D4*xyzi(i+1,j+4,k
     $        +1)+4.988141456304199D3*xyzi(i+1,j+2,k+1)-1.95613390443301
     $        9D2*xyzi(i+1,j,k+1))+angi
            angi=yrk(214)*(4.83636804631037D4*xyzi(i+13,j,k+1)+9.6727360
     $        92620741D4*xyzi(i+11,j+2,k+1)-1.432997939647517D5*xyzi(i+1
     $        1,j,k+1)-2.418184023155185D5*xyzi(i+9,j+4,k+1)-1.432997939
     $        647517D5*xyzi(i+9,j+2,k+1)+1.604957692405219D5*xyzi(i+9,j,
     $        k+1)-9.672736092620741D5*xyzi(i+7,j+6,k+1)+8.5979876378851
     $        03D5*xyzi(i+7,j+4,k+1)-8.373692308201144D4*xyzi(i+7,j,k+1)
     $        -1.209092011577593D6*xyzi(i+5,j+8,k+1)+2.006197115506524D6
     $        *xyzi(i+5,j+6,k+1)-9.629746154431316D5*xyzi(i+5,j+4,k+1)+8
     $        .373692308201144D4*xyzi(i+5,j+2,k+1)+1.993736263857415D4*x
     $        yzi(i+5,j,k+1)-6.770915264834519D5*xyzi(i+3,j+10,k+1)+1.57
     $        6297733612269D6*xyzi(i+3,j+8,k+1)-1.283966153924175D6*xyzi
     $        (i+3,j+6,k+1)+4.186846154100572D5*xyzi(i+3,j+4,k+1)-3.9874
     $        7252771483D4*xyzi(i+3,j+2,k+1)-1.678935801143087D3*xyzi(i+
     $        3,j,k+1)-1.450910413893111D5*xyzi(i+1,j+12,k+1)+4.29899381
     $        8942552D5*xyzi(i+1,j+10,k+1)-4.814873077215658D5*xyzi(i+1,
     $        j+8,k+1)+2.512107692460343D5*xyzi(i+1,j+6,k+1)-5.981208791
     $        572246D4*xyzi(i+1,j+4,k+1)+5.036807403429259D3*xyzi(i+1,j+
     $        2,k+1))+angi
            angi=yrk(216)*(-2.742853631361481D4*xyzi(i+13,j,k+1)+1.64571
     $        2178816889D5*xyzi(i+11,j+2,k+1)+7.314276350297282D4*xyzi(i
     $        +11,j,k+1)+7.954275530948294D5*xyzi(i+9,j+4,k+1)-5.1199934
     $        45208098D5*xyzi(i+9,j+2,k+1)-7.021705296285391D4*xyzi(i+9,
     $        j,k+1)+9.874273072901331D5*xyzi(i+7,j+6,k+1)-1.60914079706
     $        5402D6*xyzi(i+7,j+4,k+1)+5.617364237028313D5*xyzi(i+7,j+2,
     $        k+1)+2.849387656463637D4*xyzi(i+7,j,k+1)+2.468568268225333
     $        D5*xyzi(i+5,j+8,k+1)-1.02399868904162D6*xyzi(i+5,j+6,k+1)+
     $        9.830387414799548D5*xyzi(i+5,j+4,k+1)-2.564448890817273D5*
     $        xyzi(i+5,j+2,k+1)-4.070553794948053D3*xyzi(i+5,j,k+1)-2.74
     $        2853631361481D5*xyzi(i+3,j+10,k+1)+3.657138175148641D5*xyz
     $        i(i+3,j+8,k+1)-1.424693828231818D5*xyzi(i+3,j+4,k+1)+4.070
     $        553794948053D4*xyzi(i+3,j+2,k+1)-1.37142681568074D5*xyzi(i
     $        +1,j+12,k+1)+3.657138175148641D5*xyzi(i+1,j+10,k+1)-3.5108
     $        52648142696D5*xyzi(i+1,j+8,k+1)+1.424693828231818D5*xyzi(i
     $        +1,j+6,k+1)-2.035276897474026D4*xyzi(i+1,j+4,k+1))+angi
            angi=yrk(218)*(1.135649295191801D4*xyzi(i+13,j,k+1)-2.044168
     $        731345241D5*xyzi(i+11,j+2,k+1)-2.523665100426224D4*xyzi(i+
     $        11,j,k+1)-2.839123237979502D5*xyzi(i+9,j+4,k+1)+4.79496369
     $        0809825D5*xyzi(i+9,j+2,k+1)+1.817038872306881D4*xyzi(i+9,j
     $        ,k+1)+4.088337462690482D5*xyzi(i+7,j+6,k+1)+1.514199060255
     $        734D5*xyzi(i+7,j+4,k+1)-3.634077744613762D5*xyzi(i+7,j+2,k
     $        +1)-4.213423472015956D3*xyzi(i+7,j,k+1)+7.154590559708344D
     $        5*xyzi(i+5,j+8,k+1)-1.059939342179014D6*xyzi(i+5,j+6,k+1)+
     $        2.543854421229634D5*xyzi(i+5,j+4,k+1)+8.848189291233507D4*
     $        xyzi(i+5,j+2,k+1)+1.589909013268521D5*xyzi(i+3,j+10,k+1)-5
     $        .29969671089507D5*xyzi(i+3,j+8,k+1)+5.087708842459267D5*xy
     $        zi(i+3,j+6,k+1)-1.474698215205585D5*xyzi(i+3,j+4,k+1)-7.94
     $        9545066342605D4*xyzi(i+1,j+12,k+1)+1.766565570298357D5*xyz
     $        i(i+1,j+10,k+1)-1.271927210614817D5*xyzi(i+1,j+8,k+1)+2.94
     $        9396430411169D4*xyzi(i+1,j+6,k+1))+angi
            angi=yrk(220)*(-3.271851789497149D3*xyzi(i+13,j,k+1)+1.11242
     $        9608429031D5*xyzi(i+11,j+2,k+1)+5.331906619921279D3*xyzi(i
     $        +11,j,k+1)-1.799518484223432D5*xyzi(i+9,j+4,k+1)-1.8661673
     $        16972448D5*xyzi(i+9,j+2,k+1)-2.132762647968512D3*xyzi(i+9,
     $        j,k+1)-4.318844362136236D5*xyzi(i+7,j+6,k+1)+4.79871595792
     $        9151D5*xyzi(i+7,j+4,k+1)+7.677945532686642D4*xyzi(i+7,j+2,
     $        k+1)+1.079711090534059D5*xyzi(i+5,j+8,k+1)+2.2394007803669
     $        37D5*xyzi(i+5,j+6,k+1)-2.687280936440325D5*xyzi(i+5,j+4,k+
     $        1)+2.159422181068118D5*xyzi(i+3,j+10,k+1)-3.99892996494095
     $        9D5*xyzi(i+3,j+8,k+1)+1.79152062429355D5*xyzi(i+3,j+6,k+1)
     $        -2.944666610547434D4*xyzi(i+1,j+12,k+1)+4.798715957929151D
     $        4*xyzi(i+1,j+10,k+1)-1.919486383171661D4*xyzi(i+1,j+8,k+1)
     $        )+angi
            angi=yrk(222)*(5.973556766404109D2*xyzi(i+13,j,k+1)-3.225720
     $        653858219D4*xyzi(i+11,j+2,k+1)-5.309828236803652D2*xyzi(i+
     $        11,j,k+1)+1.64272811076113D5*xyzi(i+9,j+4,k+1)+2.920405530
     $        242009D4*xyzi(i+9,j+2,k+1)-7.885094931653424D4*xyzi(i+7,j+
     $        6,k+1)-1.752243318145205D5*xyzi(i+7,j+4,k+1)-1.77414635962
     $        202D5*xyzi(i+5,j+8,k+1)+2.453140645403287D5*xyzi(i+5,j+6,k
     $        +1)+9.199277420262327D4*xyzi(i+3,j+10,k+1)-8.7612165907260
     $        26D4*xyzi(i+3,j+8,k+1)-6.570912443044519D3*xyzi(i+1,j+12,k
     $        +1)+5.840811060484017D3*xyzi(i+1,j+10,k+1))+angi
            angi=yrk(224)*(-5.522555184144841D1*xyzi(i+13,j,k+1)+4.30759
     $        3043632976D3*xyzi(i+11,j+2,k+1)-3.948626956663561D4*xyzi(i
     $        +9,j+4,k+1)+9.476704695992546D4*xyzi(i+7,j+6,k+1)-7.107528
     $        52199441D4*xyzi(i+5,j+8,k+1)+1.579450782665424D4*xyzi(i+3,
     $        j+10,k+1)-7.179321739388293D2*xyzi(i+1,j+12,k+1))+angi
          endif
          if((ieven.and..not.jeven.and..not.keven))then
            angi=yrk(198)*(-7.179321739388293D2*xyzi(i+12,j+1,k+1)+1.579
     $        450782665424D4*xyzi(i+10,j+3,k+1)-7.10752852199441D4*xyzi(
     $        i+8,j+5,k+1)+9.476704695992546D4*xyzi(i+6,j+7,k+1)-3.94862
     $        6956663561D4*xyzi(i+4,j+9,k+1)+4.307593043632976D3*xyzi(i+
     $        2,j+11,k+1)-5.522555184144841D1*xyzi(i,j+13,k+1))+angi
            angi=yrk(200)*(6.570912443044519D3*xyzi(i+12,j+1,k+1)-9.1992
     $        77420262327D4*xyzi(i+10,j+3,k+1)-5.840811060484017D3*xyzi(
     $        i+10,j+1,k+1)+1.77414635962202D5*xyzi(i+8,j+5,k+1)+8.76121
     $        6590726026D4*xyzi(i+8,j+3,k+1)+7.885094931653424D4*xyzi(i+
     $        6,j+7,k+1)-2.453140645403287D5*xyzi(i+6,j+5,k+1)-1.6427281
     $        1076113D5*xyzi(i+4,j+9,k+1)+1.752243318145205D5*xyzi(i+4,j
     $        +7,k+1)+3.225720653858219D4*xyzi(i+2,j+11,k+1)-2.920405530
     $        242009D4*xyzi(i+2,j+9,k+1)-5.973556766404109D2*xyzi(i,j+13
     $        ,k+1)+5.309828236803652D2*xyzi(i,j+11,k+1))+angi
            angi=yrk(202)*(-2.944666610547434D4*xyzi(i+12,j+1,k+1)+2.159
     $        422181068118D5*xyzi(i+10,j+3,k+1)+4.798715957929151D4*xyzi
     $        (i+10,j+1,k+1)+1.079711090534059D5*xyzi(i+8,j+5,k+1)-3.998
     $        929964940959D5*xyzi(i+8,j+3,k+1)-1.919486383171661D4*xyzi(
     $        i+8,j+1,k+1)-4.318844362136236D5*xyzi(i+6,j+7,k+1)+2.23940
     $        0780366937D5*xyzi(i+6,j+5,k+1)+1.79152062429355D5*xyzi(i+6
     $        ,j+3,k+1)-1.799518484223432D5*xyzi(i+4,j+9,k+1)+4.79871595
     $        7929151D5*xyzi(i+4,j+7,k+1)-2.687280936440325D5*xyzi(i+4,j
     $        +5,k+1)+1.112429608429031D5*xyzi(i+2,j+11,k+1)-1.866167316
     $        972448D5*xyzi(i+2,j+9,k+1)+7.677945532686642D4*xyzi(i+2,j+
     $        7,k+1)-3.271851789497149D3*xyzi(i,j+13,k+1)+5.331906619921
     $        279D3*xyzi(i,j+11,k+1)-2.132762647968512D3*xyzi(i,j+9,k+1)
     $        )+angi
            angi=yrk(204)*(7.949545066342605D4*xyzi(i+12,j+1,k+1)-1.5899
     $        09013268521D5*xyzi(i+10,j+3,k+1)-1.766565570298357D5*xyzi(
     $        i+10,j+1,k+1)-7.154590559708344D5*xyzi(i+8,j+5,k+1)+5.2996
     $        9671089507D5*xyzi(i+8,j+3,k+1)+1.271927210614817D5*xyzi(i+
     $        8,j+1,k+1)-4.088337462690482D5*xyzi(i+6,j+7,k+1)+1.0599393
     $        42179014D6*xyzi(i+6,j+5,k+1)-5.087708842459267D5*xyzi(i+6,
     $        j+3,k+1)-2.949396430411169D4*xyzi(i+6,j+1,k+1)+2.839123237
     $        979502D5*xyzi(i+4,j+9,k+1)-1.514199060255734D5*xyzi(i+4,j+
     $        7,k+1)-2.543854421229634D5*xyzi(i+4,j+5,k+1)+1.47469821520
     $        5585D5*xyzi(i+4,j+3,k+1)+2.044168731345241D5*xyzi(i+2,j+11
     $        ,k+1)-4.794963690809825D5*xyzi(i+2,j+9,k+1)+3.634077744613
     $        762D5*xyzi(i+2,j+7,k+1)-8.848189291233507D4*xyzi(i+2,j+5,k
     $        +1)-1.135649295191801D4*xyzi(i,j+13,k+1)+2.523665100426224
     $        D4*xyzi(i,j+11,k+1)-1.817038872306881D4*xyzi(i,j+9,k+1)+4.
     $        213423472015956D3*xyzi(i,j+7,k+1))+angi
            angi=yrk(206)*(-1.37142681568074D5*xyzi(i+12,j+1,k+1)-2.7428
     $        53631361481D5*xyzi(i+10,j+3,k+1)+3.657138175148641D5*xyzi(
     $        i+10,j+1,k+1)+2.468568268225333D5*xyzi(i+8,j+5,k+1)+3.6571
     $        38175148641D5*xyzi(i+8,j+3,k+1)-3.510852648142696D5*xyzi(i
     $        +8,j+1,k+1)+9.874273072901331D5*xyzi(i+6,j+7,k+1)-1.023998
     $        68904162D6*xyzi(i+6,j+5,k+1)+1.424693828231818D5*xyzi(i+6,
     $        j+1,k+1)+7.954275530948294D5*xyzi(i+4,j+9,k+1)-1.609140797
     $        065402D6*xyzi(i+4,j+7,k+1)+9.830387414799548D5*xyzi(i+4,j+
     $        5,k+1)-1.424693828231818D5*xyzi(i+4,j+3,k+1)-2.03527689747
     $        4026D4*xyzi(i+4,j+1,k+1)+1.645712178816889D5*xyzi(i+2,j+11
     $        ,k+1)-5.119993445208098D5*xyzi(i+2,j+9,k+1)+5.617364237028
     $        313D5*xyzi(i+2,j+7,k+1)-2.564448890817273D5*xyzi(i+2,j+5,k
     $        +1)+4.070553794948053D4*xyzi(i+2,j+3,k+1)-2.74285363136148
     $        1D4*xyzi(i,j+13,k+1)+7.314276350297282D4*xyzi(i,j+11,k+1)-
     $        7.021705296285391D4*xyzi(i,j+9,k+1)+2.849387656463637D4*xy
     $        zi(i,j+7,k+1)-4.070553794948053D3*xyzi(i,j+5,k+1))+angi
            angi=yrk(208)*(1.450910413893111D5*xyzi(i+12,j+1,k+1)+6.7709
     $        15264834519D5*xyzi(i+10,j+3,k+1)-4.298993818942552D5*xyzi(
     $        i+10,j+1,k+1)+1.209092011577593D6*xyzi(i+8,j+5,k+1)-1.5762
     $        97733612269D6*xyzi(i+8,j+3,k+1)+4.814873077215658D5*xyzi(i
     $        +8,j+1,k+1)+9.672736092620741D5*xyzi(i+6,j+7,k+1)-2.006197
     $        115506524D6*xyzi(i+6,j+5,k+1)+1.283966153924175D6*xyzi(i+6
     $        ,j+3,k+1)-2.512107692460343D5*xyzi(i+6,j+1,k+1)+2.41818402
     $        3155185D5*xyzi(i+4,j+9,k+1)-8.597987637885103D5*xyzi(i+4,j
     $        +7,k+1)+9.629746154431316D5*xyzi(i+4,j+5,k+1)-4.1868461541
     $        00572D5*xyzi(i+4,j+3,k+1)+5.981208791572246D4*xyzi(i+4,j+1
     $        ,k+1)-9.672736092620741D4*xyzi(i+2,j+11,k+1)+1.43299793964
     $        7517D5*xyzi(i+2,j+9,k+1)-8.373692308201144D4*xyzi(i+2,j+5,
     $        k+1)+3.98747252771483D4*xyzi(i+2,j+3,k+1)-5.03680740342925
     $        9D3*xyzi(i+2,j+1,k+1)-4.83636804631037D4*xyzi(i,j+13,k+1)+
     $        1.432997939647517D5*xyzi(i,j+11,k+1)-1.604957692405219D5*x
     $        yzi(i,j+9,k+1)+8.373692308201144D4*xyzi(i,j+7,k+1)-1.99373
     $        6263857415D4*xyzi(i,j+5,k+1)+1.678935801143087D3*xyzi(i,j+
     $        3,k+1))+angi
            angi=yrk(210)*(-6.386185008998833D4*xyzi(i+12,j+1,k+1)-3.831
     $        7110053993D5*xyzi(i+10,j+3,k+1)+1.986813113910748D5*xyzi(i
     $        +10,j+1,k+1)-9.579277513498249D5*xyzi(i+8,j+5,k+1)+9.93406
     $        5569553739D5*xyzi(i+8,j+3,k+1)-2.384175736692897D5*xyzi(i+
     $        8,j+1,k+1)-1.277237001799767D6*xyzi(i+6,j+7,k+1)+1.9868131
     $        13910748D6*xyzi(i+6,j+5,k+1)-9.53670294677159D5*xyzi(i+6,j
     $        +3,k+1)+1.382130861850955D5*xyzi(i+6,j+1,k+1)-9.5792775134
     $        98249D5*xyzi(i+4,j+9,k+1)+1.986813113910748D6*xyzi(i+4,j+7
     $        ,k+1)-1.430505442015738D6*xyzi(i+4,j+5,k+1)+4.146392585552
     $        865D5*xyzi(i+4,j+3,k+1)-3.948945319574157D4*xyzi(i+4,j+1,k
     $        +1)-3.8317110053993D5*xyzi(i+2,j+11,k+1)+9.934065569553739
     $        D5*xyzi(i+2,j+9,k+1)-9.53670294677159D5*xyzi(i+2,j+7,k+1)+
     $        4.146392585552865D5*xyzi(i+2,j+5,k+1)-7.897890639148315D4*
     $        xyzi(i+2,j+3,k+1)+4.988141456304199D3*xyzi(i+2,j+1,k+1)-6.
     $        386185008998833D4*xyzi(i,j+13,k+1)+1.986813113910748D5*xyz
     $        i(i,j+11,k+1)-2.384175736692897D5*xyzi(i,j+9,k+1)+1.382130
     $        861850955D5*xyzi(i,j+7,k+1)-3.948945319574157D4*xyzi(i,j+5
     $        ,k+1)+4.988141456304199D3*xyzi(i,j+3,k+1)-1.95613390443301
     $        9D2*xyzi(i,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
      else
c...
c...  l=1
c...
        if(lone.eq.1)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=angi-6.139960247678931D0*yrk(4)*xyzi(i+1,j,k)
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=angi-6.139960247678931D0*yrk(2)*xyzi(i,j+1,k)
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=6.139960247678931D0*yrk(3)*xyzi(i,j,k+1)+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=3
c...
        if(lone.eq.3)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(14)*(2.871703451903279D1*xyzi(i+3,j,k)+2.8717034519
     $        03279D1*xyzi(i+1,j+2,k)-2.297362761522623D1*xyzi(i+1,j,k))
     $        +angi
            angi=yrk(16)*(2.22441192889355D1*xyzi(i+1,j+2,k)-7.414706429
     $        645167D0*xyzi(i+3,j,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(10)*(7.414706429645167D0*xyzi(i,j+3,k)-2.2244119288
     $        9355D1*xyzi(i+2,j+1,k))+angi
            angi=yrk(12)*(2.871703451903279D1*xyzi(i+2,j+1,k)+2.87170345
     $        1903279D1*xyzi(i,j+3,k)-2.297362761522623D1*xyzi(i,j+1,k))
     $        +angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(13)*(2.344736049917376D1*xyzi(i,j,k+3)-1.4068416299
     $        50425D1*xyzi(i,j,k+1))+angi
            angi=yrk(15)*(1.816224734516432D1*xyzi(i+2,j,k+1)-1.81622473
     $        4516432D1*xyzi(i,j+2,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=3.632449469032864D1*yrk(11)*xyzi(i+1,j+1,k+1)+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=5
c...
        if(lone.eq.5)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(32)*(-1.19529805236618D2*xyzi(i+5,j,k)-2.3905961047
     $        3236D2*xyzi(i+3,j+2,k)+1.59373073648824D2*xyzi(i+3,j,k)-1.
     $        19529805236618D2*xyzi(i+1,j+4,k)+1.59373073648824D2*xyzi(i
     $        +1,j+2,k)-4.5535163899664D1*xyzi(i+1,j,k))+angi
            angi=yrk(34)*(5.533154810497966D1*xyzi(i+5,j,k)-1.1066309620
     $        99593D2*xyzi(i+3,j+2,k)-4.918359831553747D1*xyzi(i+3,j,k)-
     $        1.65994644314939D2*xyzi(i+1,j+4,k)+1.475507949466124D2*xyz
     $        i(i+1,j+2,k))+angi
            angi=yrk(36)*(-8.248340190868946D0*xyzi(i+5,j,k)+8.248340190
     $        868946D1*xyzi(i+3,j+2,k)-4.124170095434473D1*xyzi(i+1,j+4,
     $        k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(26)*(-4.124170095434473D1*xyzi(i+4,j+1,k)+8.2483401
     $        90868946D1*xyzi(i+2,j+3,k)-8.248340190868946D0*xyzi(i,j+5,
     $        k))+angi
            angi=yrk(28)*(1.65994644314939D2*xyzi(i+4,j+1,k)+1.106630962
     $        099593D2*xyzi(i+2,j+3,k)-1.475507949466124D2*xyzi(i+2,j+1,
     $        k)-5.533154810497966D1*xyzi(i,j+5,k)+4.918359831553747D1*x
     $        yzi(i,j+3,k))+angi
            angi=yrk(30)*(-1.19529805236618D2*xyzi(i+4,j+1,k)-2.39059610
     $        473236D2*xyzi(i+2,j+3,k)+1.59373073648824D2*xyzi(i+2,j+1,k
     $        )-1.19529805236618D2*xyzi(i,j+5,k)+1.59373073648824D2*xyzi
     $        (i,j+3,k)-4.5535163899664D1*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(31)*(9.258738901136752D1*xyzi(i,j,k+5)-1.0287487667
     $        92972D2*xyzi(i,j,k+3)+2.204461643127798D1*xyzi(i,j,k+1))+a
     $        ngi
            angi=yrk(33)*(-9.035603969030778D1*xyzi(i+4,j,k+1)+6.0237359
     $        79353852D1*xyzi(i+2,j,k+1)+9.035603969030778D1*xyzi(i,j+4,
     $        k+1)-6.023735979353852D1*xyzi(i,j+2,k+1))+angi
            angi=yrk(35)*(2.608354191905385D1*xyzi(i+4,j,k+1)-1.56501251
     $        5143231D2*xyzi(i+2,j+2,k+1)+2.608354191905385D1*xyzi(i,j+4
     $        ,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(27)*(1.043341676762154D2*xyzi(i+3,j+1,k+1)-1.043341
     $        676762154D2*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(29)*(-1.807120793806156D2*xyzi(i+3,j+1,k+1)-1.80712
     $        0793806156D2*xyzi(i+1,j+3,k+1)+1.20474719587077D2*xyzi(i+1
     $        ,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=7
c...
        if(lone.eq.7)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(58)*(4.869752569422183D2*xyzi(i+7,j,k)+1.4609257708
     $        26655D3*xyzi(i+5,j+2,k)-8.990312435856337D2*xyzi(i+5,j,k)+
     $        1.460925770826655D3*xyzi(i+3,j+4,k)-1.798062487171267D3*xy
     $        zi(i+3,j+2,k)+4.903806783194366D2*xyzi(i+3,j,k)+4.86975256
     $        9422183D2*xyzi(i+1,j+6,k)-8.990312435856337D2*xyzi(i+1,j+4
     $        ,k)+4.903806783194366D2*xyzi(i+1,j+2,k)-7.264898938065727D
     $        1*xyzi(i+1,j,k))+angi
            angi=yrk(60)*(-2.811552956842769D2*xyzi(i+7,j,k)+2.811552956
     $        842769D2*xyzi(i+5,j+2,k)+4.325466087450414D2*xyzi(i+5,j,k)
     $        +1.405776478421384D3*xyzi(i+3,j+4,k)-8.650932174900827D2*x
     $        yzi(i+3,j+2,k)-1.572896759072878D2*xyzi(i+3,j,k)+8.4346588
     $        70528306D2*xyzi(i+1,j+6,k)-1.297639826235124D3*xyzi(i+1,j+
     $        4,k)+4.718690277218633D2*xyzi(i+1,j+2,k))+angi
            angi=yrk(62)*(8.477151123692503D1*xyzi(i+7,j,k)-7.6294360113
     $        23252D2*xyzi(i+5,j+2,k)-7.825062575716156D1*xyzi(i+5,j,k)-
     $        4.238575561846251D2*xyzi(i+3,j+4,k)+7.825062575716156D2*xy
     $        zi(i+3,j+2,k)+4.238575561846251D2*xyzi(i+1,j+6,k)-3.912531
     $        287858078D2*xyzi(i+1,j+4,k))+angi
            angi=yrk(64)*(-8.886468981567021D0*xyzi(i+7,j,k)+1.866158486
     $        129074D2*xyzi(i+5,j+2,k)-3.110264143548457D2*xyzi(i+3,j+4,
     $        k)+6.220528287096915D1*xyzi(i+1,j+6,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(50)*(-6.220528287096915D1*xyzi(i+6,j+1,k)+3.1102641
     $        43548457D2*xyzi(i+4,j+3,k)-1.866158486129074D2*xyzi(i+2,j+
     $        5,k)+8.886468981567021D0*xyzi(i,j+7,k))+angi
            angi=yrk(52)*(4.238575561846251D2*xyzi(i+6,j+1,k)-4.23857556
     $        1846251D2*xyzi(i+4,j+3,k)-3.912531287858078D2*xyzi(i+4,j+1
     $        ,k)-7.629436011323252D2*xyzi(i+2,j+5,k)+7.825062575716156D
     $        2*xyzi(i+2,j+3,k)+8.477151123692503D1*xyzi(i,j+7,k)-7.8250
     $        62575716156D1*xyzi(i,j+5,k))+angi
            angi=yrk(54)*(-8.434658870528306D2*xyzi(i+6,j+1,k)-1.4057764
     $        78421384D3*xyzi(i+4,j+3,k)+1.297639826235124D3*xyzi(i+4,j+
     $        1,k)-2.811552956842769D2*xyzi(i+2,j+5,k)+8.650932174900827
     $        D2*xyzi(i+2,j+3,k)-4.718690277218633D2*xyzi(i+2,j+1,k)+2.8
     $        11552956842769D2*xyzi(i,j+7,k)-4.325466087450414D2*xyzi(i,
     $        j+5,k)+1.572896759072878D2*xyzi(i,j+3,k))+angi
            angi=yrk(56)*(4.869752569422183D2*xyzi(i+6,j+1,k)+1.46092577
     $        0826655D3*xyzi(i+4,j+3,k)-8.990312435856337D2*xyzi(i+4,j+1
     $        ,k)+1.460925770826655D3*xyzi(i+2,j+5,k)-1.798062487171267D
     $        3*xyzi(i+2,j+3,k)+4.903806783194366D2*xyzi(i+2,j+1,k)+4.86
     $        9752569422183D2*xyzi(i,j+7,k)-8.990312435856337D2*xyzi(i,j
     $        +5,k)+4.903806783194366D2*xyzi(i,j+3,k)-7.264898938065727D
     $        1*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(57)*(3.681186927173971D2*xyzi(i,j,k+7)-5.9465327285
     $        11799D2*xyzi(i,j,k+5)+2.702969422050818D2*xyzi(i,j,k+3)-3.
     $        003299357834242D1*xyzi(i,j,k+1))+angi
            angi=yrk(59)*(3.976136322897221D2*xyzi(i+6,j,k+1)+3.97613632
     $        2897221D2*xyzi(i+4,j+2,k+1)-4.89370624356581D2*xyzi(i+4,j,
     $        k+1)-3.976136322897221D2*xyzi(i+2,j+4,k+1)+1.3346471573361
     $        3D2*xyzi(i+2,j,k+1)-3.976136322897221D2*xyzi(i,j+6,k+1)+4.
     $        89370624356581D2*xyzi(i,j+4,k+1)-1.33464715733613D2*xyzi(i
     $        ,j+2,k+1))+angi
            angi=yrk(61)*(-1.695430224738501D2*xyzi(i+6,j,k+1)+8.4771511
     $        23692502D2*xyzi(i+4,j+2,k+1)+1.304177095952693D2*xyzi(i+4,
     $        j,k+1)+8.477151123692502D2*xyzi(i+2,j+4,k+1)-7.82506257571
     $        6156D2*xyzi(i+2,j+2,k+1)-1.695430224738501D2*xyzi(i,j+6,k+
     $        1)+1.304177095952693D2*xyzi(i,j+4,k+1))+angi
            angi=yrk(63)*(3.325012230721775D1*xyzi(i+6,j,k+1)-4.98751834
     $        6082662D2*xyzi(i+4,j+2,k+1)+4.987518346082662D2*xyzi(i+2,j
     $        +4,k+1)-3.325012230721775D1*xyzi(i,j+6,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(51)*(1.995007338433065D2*xyzi(i+5,j+1,k+1)-6.650024
     $        46144355D2*xyzi(i+3,j+3,k+1)+1.995007338433065D2*xyzi(i+1,
     $        j+5,k+1))+angi
            angi=yrk(53)*(-6.781720898954002D2*xyzi(i+5,j+1,k+1)+5.21670
     $        8383810771D2*xyzi(i+3,j+1,k+1)+6.781720898954002D2*xyzi(i+
     $        1,j+5,k+1)-5.216708383810771D2*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(55)*(7.952272645794442D2*xyzi(i+5,j+1,k+1)+1.590454
     $        529158888D3*xyzi(i+3,j+3,k+1)-9.787412487131621D2*xyzi(i+3
     $        ,j+1,k+1)+7.952272645794442D2*xyzi(i+1,j+5,k+1)-9.78741248
     $        7131621D2*xyzi(i+1,j+3,k+1)+2.66929431467226D2*xyzi(i+1,j+
     $        1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=9
c...
        if(lone.eq.9)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(92)*(-1.968624920969132D3*xyzi(i+9,j,k)-7.874499683
     $        876526D3*xyzi(i+7,j+2,k)+4.632058637574427D3*xyzi(i+7,j,k)
     $        -1.181174952581479D4*xyzi(i+5,j+4,k)+1.389617591272328D4*x
     $        yzi(i+5,j+2,k)-3.705646910059542D3*xyzi(i+5,j,k)-7.8744996
     $        83876526D3*xyzi(i+3,j+6,k)+1.389617591272328D4*xyzi(i+3,j+
     $        4,k)-7.411293820119083D3*xyzi(i+3,j+2,k)+1.14019904924909D
     $        3*xyzi(i+3,j,k)-1.968624920969132D3*xyzi(i+1,j+8,k)+4.6320
     $        58637574427D3*xyzi(i+1,j+6,k)-3.705646910059542D3*xyzi(i+1
     $        ,j+4,k)+1.14019904924909D3*xyzi(i+1,j+2,k)-1.0365445902264
     $        45D2*xyzi(i+1,j,k))+angi
            angi=yrk(94)*(1.2822420836111D3*xyzi(i+9,j,k)-2.715336177058
     $        8D3*xyzi(i+7,j,k)-7.693452501666601D3*xyzi(i+5,j+4,k)+2.71
     $        53361770588D3*xyzi(i+5,j+2,k)+1.8102241180392D3*xyzi(i+5,j
     $        ,k)-1.02579366688888D4*xyzi(i+3,j+6,k)+1.3576680885294D4*x
     $        yzi(i+3,j+4,k)-3.6204482360784D3*xyzi(i+3,j+2,k)-3.7132802
     $        42131693D2*xyzi(i+3,j,k)-3.8467262508333D3*xyzi(i+1,j+8,k)
     $        +8.146008531176401D3*xyzi(i+1,j+6,k)-5.430672354117601D3*x
     $        yzi(i+1,j+4,k)+1.113984072639508D3*xyzi(i+1,j+2,k))+angi
            angi=yrk(96)*(-5.205889671223938D2*xyzi(i+9,j,k)+4.164711736
     $        97915D3*xyzi(i+7,j+2,k)+8.574406517310015D2*xyzi(i+7,j,k)+
     $        7.288245539713513D3*xyzi(i+5,j+4,k)-7.716965865579014D3*xy
     $        zi(i+5,j+2,k)-3.429762606924006D2*xyzi(i+5,j,k)-4.28720325
     $        8655008D3*xyzi(i+3,j+4,k)+3.429762606924006D3*xyzi(i+3,j+2
     $        ,k)-2.602944835611969D3*xyzi(i+1,j+8,k)+4.287203258655008D
     $        3*xyzi(i+1,j+6,k)-1.714881303462003D3*xyzi(i+1,j+4,k))+ang
     $        i
            angi=yrk(98)*(1.164072318822076D2*xyzi(i+9,j,k)-2.3281446376
     $        44151D3*xyzi(i+7,j+2,k)-1.095597476538424D2*xyzi(i+7,j,k)+
     $        1.629701246350906D3*xyzi(i+5,j+4,k)+2.300754700730691D3*xy
     $        zi(i+5,j+2,k)+3.259402492701812D3*xyzi(i+3,j+6,k)-3.834591
     $        167884484D3*xyzi(i+3,j+4,k)-8.148506231754529D2*xyzi(i+1,j
     $        +8,k)+7.669182335768969D2*xyzi(i+1,j+6,k))+angi
            angi=yrk(100)*(-9.410966914433519D0*xyzi(i+9,j,k)+3.38794808
     $        9196067D2*xyzi(i+7,j+2,k)-1.185781831218623D3*xyzi(i+5,j+4
     $        ,k)+7.905212208124156D2*xyzi(i+3,j+6,k)-8.469870222990167D
     $        1*xyzi(i+1,j+8,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(82)*(-8.469870222990167D1*xyzi(i+8,j+1,k)+7.9052122
     $        08124156D2*xyzi(i+6,j+3,k)-1.185781831218623D3*xyzi(i+4,j+
     $        5,k)+3.387948089196067D2*xyzi(i+2,j+7,k)-9.410966914433519
     $        D0*xyzi(i,j+9,k))+angi
            angi=yrk(84)*(8.148506231754529D2*xyzi(i+8,j+1,k)-3.25940249
     $        2701812D3*xyzi(i+6,j+3,k)-7.669182335768969D2*xyzi(i+6,j+1
     $        ,k)-1.629701246350906D3*xyzi(i+4,j+5,k)+3.834591167884484D
     $        3*xyzi(i+4,j+3,k)+2.328144637644151D3*xyzi(i+2,j+7,k)-2.30
     $        0754700730691D3*xyzi(i+2,j+5,k)-1.164072318822076D2*xyzi(i
     $        ,j+9,k)+1.095597476538424D2*xyzi(i,j+7,k))+angi
            angi=yrk(86)*(-2.602944835611969D3*xyzi(i+8,j+1,k)+4.2872032
     $        58655008D3*xyzi(i+6,j+1,k)+7.288245539713513D3*xyzi(i+4,j+
     $        5,k)-4.287203258655008D3*xyzi(i+4,j+3,k)-1.714881303462003
     $        D3*xyzi(i+4,j+1,k)+4.16471173697915D3*xyzi(i+2,j+7,k)-7.71
     $        6965865579014D3*xyzi(i+2,j+5,k)+3.429762606924006D3*xyzi(i
     $        +2,j+3,k)-5.205889671223938D2*xyzi(i,j+9,k)+8.574406517310
     $        015D2*xyzi(i,j+7,k)-3.429762606924006D2*xyzi(i,j+5,k))+ang
     $        i
            angi=yrk(88)*(3.8467262508333D3*xyzi(i+8,j+1,k)+1.0257936668
     $        8888D4*xyzi(i+6,j+3,k)-8.146008531176401D3*xyzi(i+6,j+1,k)
     $        +7.693452501666601D3*xyzi(i+4,j+5,k)-1.3576680885294D4*xyz
     $        i(i+4,j+3,k)+5.430672354117601D3*xyzi(i+4,j+1,k)-2.7153361
     $        770588D3*xyzi(i+2,j+5,k)+3.6204482360784D3*xyzi(i+2,j+3,k)
     $        -1.113984072639508D3*xyzi(i+2,j+1,k)-1.2822420836111D3*xyz
     $        i(i,j+9,k)+2.7153361770588D3*xyzi(i,j+7,k)-1.8102241180392
     $        D3*xyzi(i,j+5,k)+3.713280242131693D2*xyzi(i,j+3,k))+angi
            angi=yrk(90)*(-1.968624920969132D3*xyzi(i+8,j+1,k)-7.8744996
     $        83876526D3*xyzi(i+6,j+3,k)+4.632058637574427D3*xyzi(i+6,j+
     $        1,k)-1.181174952581479D4*xyzi(i+4,j+5,k)+1.389617591272328
     $        D4*xyzi(i+4,j+3,k)-3.705646910059542D3*xyzi(i+4,j+1,k)-7.8
     $        74499683876526D3*xyzi(i+2,j+7,k)+1.389617591272328D4*xyzi(
     $        i+2,j+5,k)-7.411293820119083D3*xyzi(i+2,j+3,k)+1.140199049
     $        24909D3*xyzi(i+2,j+1,k)-1.968624920969132D3*xyzi(i,j+9,k)+
     $        4.632058637574427D3*xyzi(i,j+7,k)-3.705646910059542D3*xyzi
     $        (i,j+5,k)+1.14019904924909D3*xyzi(i,j+3,k)-1.0365445902264
     $        45D2*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(91)*(1.467326381829043D3*xyzi(i,j,k+9)-3.1072793968
     $        14444D3*xyzi(i,j,k+7)+2.175095577770111D3*xyzi(i,j,k+5)-5.
     $        57716814812849D2*xyzi(i,j,k+3)+3.802614646451243D1*xyzi(i,
     $        j,k+1))+angi
            angi=yrk(93)*(-1.678848973544503D3*xyzi(i+8,j,k+1)-3.3576979
     $        47089007D3*xyzi(i+6,j+2,k+1)+2.962674659196182D3*xyzi(i+6,
     $        j,k+1)+2.962674659196182D3*xyzi(i+4,j+2,k+1)-1.58009315157
     $        1297D3*xyzi(i+4,j,k+1)+3.357697947089007D3*xyzi(i+2,j+6,k+
     $        1)-2.962674659196182D3*xyzi(i+2,j+4,k+1)+2.430912540878919
     $        D2*xyzi(i+2,j,k+1)+1.678848973544503D3*xyzi(i,j+8,k+1)-2.9
     $        62674659196182D3*xyzi(i,j+6,k+1)+1.580093151571297D3*xyzi(
     $        i,j+4,k+1)-2.430912540878919D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(95)*(8.711119580919379D2*xyzi(i+8,j,k+1)-3.48444783
     $        2367752D3*xyzi(i+6,j+2,k+1)-1.229805117306265D3*xyzi(i+6,j
     $        ,k+1)-8.711119580919379D3*xyzi(i+4,j+4,k+1)+6.149025586531
     $        326D3*xyzi(i+4,j+2,k+1)+4.099350391020884D2*xyzi(i+4,j,k+1
     $        )-3.484447832367752D3*xyzi(i+2,j+6,k+1)+6.149025586531326D
     $        3*xyzi(i+2,j+4,k+1)-2.459610234612531D3*xyzi(i+2,j+2,k+1)+
     $        8.711119580919379D2*xyzi(i,j+8,k+1)-1.229805117306265D3*xy
     $        zi(i,j+6,k+1)+4.099350391020884D2*xyzi(i,j+4,k+1))+angi
            angi=yrk(97)*(-2.688309866512469D2*xyzi(i+8,j,k+1)+3.7636338
     $        13117456D3*xyzi(i+6,j+2,k+1)+2.213902243010268D2*xyzi(i+6,
     $        j,k+1)-3.320853364515403D3*xyzi(i+4,j+2,k+1)-3.76363381311
     $        7456D3*xyzi(i+2,j+6,k+1)+3.320853364515403D3*xyzi(i+2,j+4,
     $        k+1)+2.688309866512469D2*xyzi(i,j+8,k+1)-2.213902243010268
     $        D2*xyzi(i,j+6,k+1))+angi
            angi=yrk(99)*(3.992735113630909D1*xyzi(i+8,j,k+1)-1.11796583
     $        1816654D3*xyzi(i+6,j+2,k+1)+2.794914579541636D3*xyzi(i+4,j
     $        +4,k+1)-1.117965831816654D3*xyzi(i+2,j+6,k+1)+3.9927351136
     $        30909D1*xyzi(i,j+8,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(83)*(3.194188090904727D2*xyzi(i+7,j+1,k+1)-2.235931
     $        663633309D3*xyzi(i+5,j+3,k+1)+2.235931663633309D3*xyzi(i+3
     $        ,j+5,k+1)-3.194188090904727D2*xyzi(i+1,j+7,k+1))+angi
            angi=yrk(85)*(-1.612985919907481D3*xyzi(i+7,j+1,k+1)+3.76363
     $        3813117456D3*xyzi(i+5,j+3,k+1)+1.328341345806161D3*xyzi(i+
     $        5,j+1,k+1)+3.763633813117456D3*xyzi(i+3,j+5,k+1)-4.4278044
     $        86020537D3*xyzi(i+3,j+3,k+1)-1.612985919907481D3*xyzi(i+1,
     $        j+7,k+1)+1.328341345806161D3*xyzi(i+1,j+5,k+1))+angi
            angi=yrk(87)*(3.484447832367752D3*xyzi(i+7,j+1,k+1)+3.484447
     $        832367752D3*xyzi(i+5,j+3,k+1)-4.919220469225061D3*xyzi(i+5
     $        ,j+1,k+1)-3.484447832367752D3*xyzi(i+3,j+5,k+1)+1.63974015
     $        6408354D3*xyzi(i+3,j+1,k+1)-3.484447832367752D3*xyzi(i+1,j
     $        +7,k+1)+4.919220469225061D3*xyzi(i+1,j+5,k+1)-1.6397401564
     $        08354D3*xyzi(i+1,j+3,k+1))+angi
            angi=yrk(89)*(-3.357697947089007D3*xyzi(i+7,j+1,k+1)-1.00730
     $        9384126702D4*xyzi(i+5,j+3,k+1)+5.925349318392365D3*xyzi(i+
     $        5,j+1,k+1)-1.007309384126702D4*xyzi(i+3,j+5,k+1)+1.1850698
     $        63678473D4*xyzi(i+3,j+3,k+1)-3.160186303142594D3*xyzi(i+3,
     $        j+1,k+1)-3.357697947089007D3*xyzi(i+1,j+7,k+1)+5.925349318
     $        392365D3*xyzi(i+1,j+5,k+1)-3.160186303142594D3*xyzi(i+1,j+
     $        3,k+1)+4.861825081757838D2*xyzi(i+1,j+1,k+1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=11
c...
        if(lone.eq.11)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(134)*(7.928933427523546D3*xyzi(i+11,j,k)+3.96446671
     $        3761773D4*xyzi(i+9,j+2,k)-2.265409550721013D4*xyzi(i+9,j,k
     $        )+7.928933427523546D4*xyzi(i+7,j+4,k)-9.061638202884053D4*
     $        xyzi(i+7,j+2,k)+2.384641632337909D4*xyzi(i+7,j,k)+7.928933
     $        427523546D4*xyzi(i+5,j+6,k)-1.359245730432608D5*xyzi(i+5,j
     $        +4,k)+7.153924897013726D4*xyzi(i+5,j+2,k)-1.12218429757078
     $        1D4*xyzi(i+5,j,k)+3.964466713761773D4*xyzi(i+3,j+8,k)-9.06
     $        1638202884053D4*xyzi(i+3,j+6,k)+7.153924897013726D4*xyzi(i
     $        +3,j+4,k)-2.244368595141561D4*xyzi(i+3,j+2,k)+2.2443685951
     $        41561D3*xyzi(i+3,j,k)+7.928933427523546D3*xyzi(i+1,j+10,k)
     $        -2.265409550721013D4*xyzi(i+1,j+8,k)+2.384641632337909D4*x
     $        yzi(i+1,j+6,k)-1.122184297570781D4*xyzi(i+1,j+4,k)+2.24436
     $        8595141561D3*xyzi(i+1,j+2,k)-1.381149904702499D2*xyzi(i+1,
     $        j,k))+angi
            angi=yrk(136)*(-5.575711986679481D3*xyzi(i+11,j,k)-5.5757119
     $        86679481D3*xyzi(i+9,j+2,k)+1.486856529781195D4*xyzi(i+9,j,
     $        k)+3.345427192007688D4*xyzi(i+7,j+4,k)-1.408600922950606D4
     $        *xyzi(i+7,j,k)+7.805996781351273D4*xyzi(i+5,j+6,k)-8.92113
     $        9178687169D4*xyzi(i+5,j+4,k)+1.408600922950606D4*xyzi(i+5,
     $        j+2,k)+5.523925188041591D3*xyzi(i+5,j,k)+6.133283185347429
     $        D4*xyzi(i+3,j+8,k)-1.189485223824956D5*xyzi(i+3,j+6,k)+7.0
     $        43004614753028D4*xyzi(i+3,j+4,k)-1.104785037608318D4*xyzi(
     $        i+3,j+2,k)-7.365233584055455D2*xyzi(i+3,j,k)+1.67271359600
     $        3844D4*xyzi(i+1,j+10,k)-4.460569589343585D4*xyzi(i+1,j+8,k
     $        )+4.225802768851817D4*xyzi(i+1,j+6,k)-1.657177556412477D4*
     $        xyzi(i+1,j+4,k)+2.209570075216636D3*xyzi(i+1,j+2,k))+angi
            angi=yrk(138)*(2.693324767573892D3*xyzi(i+11,j,k)-1.88532733
     $        7301724D4*xyzi(i+9,j+2,k)-6.156170897311752D3*xyzi(i+9,j,k
     $        )-5.925314488662561D4*xyzi(i+7,j+4,k)+4.924936717849402D4*
     $        xyzi(i+7,j+2,k)+4.536125924334975D3*xyzi(i+7,j,k)-3.770654
     $        674603448D4*xyzi(i+5,j+6,k)+8.618639256236453D4*xyzi(i+5,j
     $        +4,k)-4.082513331901478D4*xyzi(i+5,j+2,k)-1.06732374690234
     $        7D3*xyzi(i+5,j,k)+1.346662383786946D4*xyzi(i+3,j+8,k)-2.26
     $        8062962167488D4*xyzi(i+3,j+4,k)+1.067323746902347D4*xyzi(i
     $        +3,j+2,k)+1.346662383786946D4*xyzi(i+1,j+10,k)-3.078085448
     $        655876D4*xyzi(i+1,j+8,k)+2.268062962167488D4*xyzi(i+1,j+6,
     $        k)-5.336618734511736D3*xyzi(i+1,j+4,k))+angi
            angi=yrk(140)*(-8.433126966180176D2*xyzi(i+11,j,k)+1.6022941
     $        23574233D4*xyzi(i+9,j+2,k)+1.44567890848803D3*xyzi(i+9,j,k
     $        )+5.059876179708106D3*xyzi(i+7,j+4,k)-2.89135781697606D4*x
     $        yzi(i+7,j+2,k)-6.087069088370653D2*xyzi(i+7,j,k)-3.5419133
     $        25795674D4*xyzi(i+5,j+6,k)+2.023950471883242D4*xyzi(i+5,j+
     $        4,k)+1.278284508557837D4*xyzi(i+5,j+2,k)-1.770956662897837
     $        D4*xyzi(i+3,j+8,k)+4.047900943766485D4*xyzi(i+3,j+6,k)-2.1
     $        30474180929729D4*xyzi(i+3,j+4,k)+5.903188876326123D3*xyzi(
     $        i+1,j+10,k)-1.011975235941621D4*xyzi(i+1,j+8,k)+4.26094836
     $        1859457D3*xyzi(i+1,j+6,k))+angi
            angi=yrk(142)*(1.49860598832503D2*xyzi(i+11,j,k)-5.245120959
     $        137605D3*xyzi(i+9,j+2,k)-1.42724379840479D2*xyzi(i+9,j,k)+
     $        1.348745389492527D4*xyzi(i+7,j+4,k)+5.138077674257246D3*xy
     $        zi(i+7,j+2,k)+6.294145150965126D3*xyzi(i+5,j+6,k)-1.798327
     $        185990036D4*xyzi(i+5,j+4,k)-1.123954491243772D4*xyzi(i+3,j
     $        +8,k)+1.198884790660024D4*xyzi(i+3,j+6,k)+1.34874538949252
     $        7D3*xyzi(i+1,j+10,k)-1.284519418564311D3*xyzi(i+1,j+8,k))+
     $        angi
            angi=yrk(144)*(-9.860103500953133D0*xyzi(i+11,j,k)+5.4230569
     $        25524223D2*xyzi(i+9,j+2,k)-3.253834155314534D3*xyzi(i+7,j+
     $        4,k)+4.555367817440347D3*xyzi(i+5,j+6,k)-1.626917077657267
     $        D3*xyzi(i+3,j+8,k)+1.084611385104845D2*xyzi(i+1,j+10,k))+a
     $        ngi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(122)*(-1.084611385104845D2*xyzi(i+10,j+1,k)+1.62691
     $        7077657267D3*xyzi(i+8,j+3,k)-4.555367817440347D3*xyzi(i+6,
     $        j+5,k)+3.253834155314534D3*xyzi(i+4,j+7,k)-5.4230569255242
     $        23D2*xyzi(i+2,j+9,k)+9.860103500953133D0*xyzi(i,j+11,k))+a
     $        ngi
            angi=yrk(124)*(1.348745389492527D3*xyzi(i+10,j+1,k)-1.123954
     $        491243772D4*xyzi(i+8,j+3,k)-1.284519418564311D3*xyzi(i+8,j
     $        +1,k)+6.294145150965126D3*xyzi(i+6,j+5,k)+1.19888479066002
     $        4D4*xyzi(i+6,j+3,k)+1.348745389492527D4*xyzi(i+4,j+7,k)-1.
     $        798327185990036D4*xyzi(i+4,j+5,k)-5.245120959137605D3*xyzi
     $        (i+2,j+9,k)+5.138077674257246D3*xyzi(i+2,j+7,k)+1.49860598
     $        832503D2*xyzi(i,j+11,k)-1.42724379840479D2*xyzi(i,j+9,k))+
     $        angi
            angi=yrk(126)*(-5.903188876326123D3*xyzi(i+10,j+1,k)+1.77095
     $        6662897837D4*xyzi(i+8,j+3,k)+1.011975235941621D4*xyzi(i+8,
     $        j+1,k)+3.541913325795674D4*xyzi(i+6,j+5,k)-4.0479009437664
     $        85D4*xyzi(i+6,j+3,k)-4.260948361859457D3*xyzi(i+6,j+1,k)-5
     $        .059876179708106D3*xyzi(i+4,j+7,k)-2.023950471883242D4*xyz
     $        i(i+4,j+5,k)+2.130474180929729D4*xyzi(i+4,j+3,k)-1.6022941
     $        23574233D4*xyzi(i+2,j+9,k)+2.89135781697606D4*xyzi(i+2,j+7
     $        ,k)-1.278284508557837D4*xyzi(i+2,j+5,k)+8.433126966180176D
     $        2*xyzi(i,j+11,k)-1.44567890848803D3*xyzi(i,j+9,k)+6.087069
     $        088370653D2*xyzi(i,j+7,k))+angi
            angi=yrk(128)*(1.346662383786946D4*xyzi(i+10,j+1,k)+1.346662
     $        383786946D4*xyzi(i+8,j+3,k)-3.078085448655876D4*xyzi(i+8,j
     $        +1,k)-3.770654674603448D4*xyzi(i+6,j+5,k)+2.26806296216748
     $        8D4*xyzi(i+6,j+1,k)-5.925314488662561D4*xyzi(i+4,j+7,k)+8.
     $        618639256236453D4*xyzi(i+4,j+5,k)-2.268062962167488D4*xyzi
     $        (i+4,j+3,k)-5.336618734511736D3*xyzi(i+4,j+1,k)-1.88532733
     $        7301724D4*xyzi(i+2,j+9,k)+4.924936717849402D4*xyzi(i+2,j+7
     $        ,k)-4.082513331901478D4*xyzi(i+2,j+5,k)+1.067323746902347D
     $        4*xyzi(i+2,j+3,k)+2.693324767573892D3*xyzi(i,j+11,k)-6.156
     $        170897311752D3*xyzi(i,j+9,k)+4.536125924334975D3*xyzi(i,j+
     $        7,k)-1.067323746902347D3*xyzi(i,j+5,k))+angi
            angi=yrk(130)*(-1.672713596003844D4*xyzi(i+10,j+1,k)-6.13328
     $        3185347429D4*xyzi(i+8,j+3,k)+4.460569589343585D4*xyzi(i+8,
     $        j+1,k)-7.805996781351273D4*xyzi(i+6,j+5,k)+1.1894852238249
     $        56D5*xyzi(i+6,j+3,k)-4.225802768851817D4*xyzi(i+6,j+1,k)-3
     $        .345427192007688D4*xyzi(i+4,j+7,k)+8.921139178687169D4*xyz
     $        i(i+4,j+5,k)-7.043004614753028D4*xyzi(i+4,j+3,k)+1.6571775
     $        56412477D4*xyzi(i+4,j+1,k)+5.575711986679481D3*xyzi(i+2,j+
     $        9,k)-1.408600922950606D4*xyzi(i+2,j+5,k)+1.104785037608318
     $        D4*xyzi(i+2,j+3,k)-2.209570075216636D3*xyzi(i+2,j+1,k)+5.5
     $        75711986679481D3*xyzi(i,j+11,k)-1.486856529781195D4*xyzi(i
     $        ,j+9,k)+1.408600922950606D4*xyzi(i,j+7,k)-5.52392518804159
     $        1D3*xyzi(i,j+5,k)+7.365233584055455D2*xyzi(i,j+3,k))+angi
            angi=yrk(132)*(7.928933427523546D3*xyzi(i+10,j+1,k)+3.964466
     $        713761773D4*xyzi(i+8,j+3,k)-2.265409550721013D4*xyzi(i+8,j
     $        +1,k)+7.928933427523546D4*xyzi(i+6,j+5,k)-9.06163820288405
     $        3D4*xyzi(i+6,j+3,k)+2.384641632337909D4*xyzi(i+6,j+1,k)+7.
     $        928933427523546D4*xyzi(i+4,j+7,k)-1.359245730432608D5*xyzi
     $        (i+4,j+5,k)+7.153924897013726D4*xyzi(i+4,j+3,k)-1.12218429
     $        7570781D4*xyzi(i+4,j+1,k)+3.964466713761773D4*xyzi(i+2,j+9
     $        ,k)-9.061638202884053D4*xyzi(i+2,j+7,k)+7.153924897013726D
     $        4*xyzi(i+2,j+5,k)-2.244368595141561D4*xyzi(i+2,j+3,k)+2.24
     $        4368595141561D3*xyzi(i+2,j+1,k)+7.928933427523546D3*xyzi(i
     $        ,j+11,k)-2.265409550721013D4*xyzi(i,j+9,k)+2.3846416323379
     $        09D4*xyzi(i,j+7,k)-1.122184297570781D4*xyzi(i,j+5,k)+2.244
     $        368595141561D3*xyzi(i,j+3,k)-1.381149904702499D2*xyzi(i,j+
     $        1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(133)*(5.855905424818466D3*xyzi(i,j,k+11)-1.53368951
     $        6023884D4*xyzi(i,j,k+9)+1.452969015180522D4*xyzi(i,j,k+7)-
     $        5.982813591919795D3*xyzi(i,j,k+5)+9.971355986532992D2*xyzi
     $        (i,j,k+3)-4.602164301476766D1*xyzi(i,j,k+1))+angi
            angi=yrk(135)*(6.954134647161096D3*xyzi(i+10,j,k+1)+2.086240
     $        394148329D4*xyzi(i+8,j+2,k+1)-1.589516490779679D4*xyzi(i+8
     $        ,j,k+1)+1.390826929432219D4*xyzi(i+6,j+4,k+1)-3.1790329815
     $        59358D4*xyzi(i+6,j+2,k+1)+1.25488144008922D4*xyzi(i+6,j,k+
     $        1)-1.390826929432219D4*xyzi(i+4,j+6,k+1)+1.25488144008922D
     $        4*xyzi(i+4,j+2,k+1)-3.936882949299515D3*xyzi(i+4,j,k+1)-2.
     $        086240394148329D4*xyzi(i+2,j+8,k+1)+3.179032981559358D4*xy
     $        zi(i+2,j+6,k+1)-1.25488144008922D4*xyzi(i+2,j+4,k+1)+3.936
     $        882949299515D2*xyzi(i+2,j,k+1)-6.954134647161096D3*xyzi(i,
     $        j+10,k+1)+1.589516490779679D4*xyzi(i,j+8,k+1)-1.2548814400
     $        8922D4*xyzi(i,j+6,k+1)+3.936882949299515D3*xyzi(i,j+4,k+1)
     $        -3.936882949299515D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(137)*(-4.071924305675061D3*xyzi(i+10,j,k+1)+1.22157
     $        7291702518D4*xyzi(i+8,j+2,k+1)+8.143848611350123D3*xyzi(i+
     $        8,j,k+1)+5.700694027945086D4*xyzi(i+6,j+4,k+1)-3.257539444
     $        540049D4*xyzi(i+6,j+2,k+1)-5.143483333484288D3*xyzi(i+6,j,
     $        k+1)+5.700694027945086D4*xyzi(i+4,j+6,k+1)-8.1438486113501
     $        23D4*xyzi(i+4,j+4,k+1)+2.571741666742144D4*xyzi(i+4,j+2,k+
     $        1)+1.008526143820449D3*xyzi(i+4,j,k+1)+1.221577291702518D4
     $        *xyzi(i+2,j+8,k+1)-3.257539444540049D4*xyzi(i+2,j+6,k+1)+2
     $        .571741666742144D4*xyzi(i+2,j+4,k+1)-6.051156862922692D3*x
     $        yzi(i+2,j+2,k+1)-4.071924305675061D3*xyzi(i,j+10,k+1)+8.14
     $        3848611350123D3*xyzi(i,j+8,k+1)-5.143483333484288D3*xyzi(i
     $        ,j+6,k+1)+1.008526143820449D3*xyzi(i,j+4,k+1))+angi
            angi=yrk(139)*(1.600073340630907D3*xyzi(i+10,j,k+1)-2.080095
     $        342820179D4*xyzi(i+8,j+2,k+1)-2.438206995247096D3*xyzi(i+8
     $        ,j,k+1)-2.240102676883269D4*xyzi(i+6,j+4,k+1)+3.4134897933
     $        45934D4*xyzi(i+6,j+2,k+1)+8.982867877226143D2*xyzi(i+6,j,k
     $        +1)+2.240102676883269D4*xyzi(i+4,j+6,k+1)-1.34743018158392
     $        1D4*xyzi(i+4,j+2,k+1)+2.080095342820179D4*xyzi(i+2,j+8,k+1
     $        )-3.413489793345934D4*xyzi(i+2,j+6,k+1)+1.347430181583921D
     $        4*xyzi(i+2,j+4,k+1)-1.600073340630907D3*xyzi(i,j+10,k+1)+2
     $        .438206995247096D3*xyzi(i,j+8,k+1)-8.982867877226143D2*xyz
     $        i(i,j+6,k+1))+angi
            angi=yrk(141)*(-3.869384023539698D2*xyzi(i+10,j,k+1)+1.04473
     $        3686355719D4*xyzi(i+8,j+2,k+1)+3.316614877319741D2*xyzi(i+
     $        8,j,k+1)-1.625141289886673D4*xyzi(i+6,j+4,k+1)-9.286521656
     $        495276D3*xyzi(i+6,j+2,k+1)-1.625141289886673D4*xyzi(i+4,j+
     $        6,k+1)+2.321630414123819D4*xyzi(i+4,j+4,k+1)+1.04473368635
     $        5719D4*xyzi(i+2,j+8,k+1)-9.286521656495276D3*xyzi(i+2,j+6,
     $        k+1)-3.869384023539698D2*xyzi(i,j+10,k+1)+3.31661487731974
     $        1D2*xyzi(i,j+8,k+1))+angi
            angi=yrk(143)*(4.624798485436075D1*xyzi(i+10,j,k+1)-2.081159
     $        318446234D3*xyzi(i+8,j+2,k+1)+9.712076819415756D3*xyzi(i+6
     $        ,j+4,k+1)-9.712076819415756D3*xyzi(i+4,j+6,k+1)+2.08115931
     $        8446234D3*xyzi(i+2,j+8,k+1)-4.624798485436075D1*xyzi(i,j+1
     $        0,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(123)*(4.624798485436074D2*xyzi(i+9,j+1,k+1)-5.54975
     $        8182523289D3*xyzi(i+7,j+3,k+1)+1.165449218329891D4*xyzi(i+
     $        5,j+5,k+1)-5.549758182523289D3*xyzi(i+3,j+7,k+1)+4.6247984
     $        85436074D2*xyzi(i+1,j+9,k+1))+angi
            angi=yrk(125)*(-3.095507218831759D3*xyzi(i+9,j+1,k+1)+1.8573
     $        04331299055D4*xyzi(i+7,j+3,k+1)+2.653291901855793D3*xyzi(i
     $        +7,j+1,k+1)-1.857304331299055D4*xyzi(i+5,j+3,k+1)-1.857304
     $        331299055D4*xyzi(i+3,j+7,k+1)+1.857304331299055D4*xyzi(i+3
     $        ,j+5,k+1)+3.095507218831759D3*xyzi(i+1,j+9,k+1)-2.65329190
     $        1855793D3*xyzi(i+1,j+7,k+1))+angi
            angi=yrk(127)*(9.60044004378544D3*xyzi(i+9,j+1,k+1)-1.280058
     $        672504725D4*xyzi(i+7,j+3,k+1)-1.462924197148258D4*xyzi(i+7
     $        ,j+1,k+1)-4.480205353766539D4*xyzi(i+5,j+5,k+1)+3.41348979
     $        3345934D4*xyzi(i+5,j+3,k+1)+5.389720726335685D3*xyzi(i+5,j
     $        +1,k+1)-1.280058672504725D4*xyzi(i+3,j+7,k+1)+3.4134897933
     $        45934D4*xyzi(i+3,j+5,k+1)-1.796573575445229D4*xyzi(i+3,j+3
     $        ,k+1)+9.60044004378544D3*xyzi(i+1,j+9,k+1)-1.4629241971482
     $        58D4*xyzi(i+1,j+7,k+1)+5.389720726335685D3*xyzi(i+1,j+5,k+
     $        1))+angi
            angi=yrk(129)*(-1.628769722270025D4*xyzi(i+9,j+1,k+1)-3.2575
     $        39444540049D4*xyzi(i+7,j+3,k+1)+3.257539444540049D4*xyzi(i
     $        +7,j+1,k+1)+3.257539444540049D4*xyzi(i+5,j+3,k+1)-2.057393
     $        333393715D4*xyzi(i+5,j+1,k+1)+3.257539444540049D4*xyzi(i+3
     $        ,j+7,k+1)-3.257539444540049D4*xyzi(i+3,j+5,k+1)+4.03410457
     $        5281795D3*xyzi(i+3,j+1,k+1)+1.628769722270025D4*xyzi(i+1,j
     $        +9,k+1)-3.257539444540049D4*xyzi(i+1,j+7,k+1)+2.0573933333
     $        93715D4*xyzi(i+1,j+5,k+1)-4.034104575281795D3*xyzi(i+1,j+3
     $        ,k+1))+angi
            angi=yrk(131)*(1.390826929432219D4*xyzi(i+9,j+1,k+1)+5.56330
     $        7717728877D4*xyzi(i+7,j+3,k+1)-3.179032981559358D4*xyzi(i+
     $        7,j+1,k+1)+8.344961576593315D4*xyzi(i+5,j+5,k+1)-9.5370989
     $        44678074D4*xyzi(i+5,j+3,k+1)+2.509762880178441D4*xyzi(i+5,
     $        j+1,k+1)+5.563307717728877D4*xyzi(i+3,j+7,k+1)-9.537098944
     $        678074D4*xyzi(i+3,j+5,k+1)+5.019525760356881D4*xyzi(i+3,j+
     $        3,k+1)-7.873765898599029D3*xyzi(i+3,j+1,k+1)+1.39082692943
     $        2219D4*xyzi(i+1,j+9,k+1)-3.179032981559358D4*xyzi(i+1,j+7,
     $        k+1)+2.509762880178441D4*xyzi(i+1,j+5,k+1)-7.8737658985990
     $        29D3*xyzi(i+1,j+3,k+1)+7.873765898599029D2*xyzi(i+1,j+1,k+
     $        1))+angi
          endif
          OMEGA=angi
        endif
c...
c...  l=13
c...
        if(lone.eq.13)then
          angi=ZERO
          if((.not.ieven.and.jeven.and.keven))then
            angi=yrk(184)*(-3.186969598569298D4*xyzi(i+13,j,k)-1.9121817
     $        59141579D5*xyzi(i+11,j+2,k)+1.070821785119284D5*xyzi(i+11,
     $        j,k)-4.780454397853947D5*xyzi(i+9,j+4,k)+5.35410892559642D
     $        5*xyzi(i+9,j+2,k)-1.396724067546892D5*xyzi(i+9,j,k)-6.3739
     $        39197138596D5*xyzi(i+7,j+6,k)+1.070821785119284D6*xyzi(i+7
     $        ,j+4,k)-5.586896270187569D5*xyzi(i+7,j+2,k)+8.868089317758
     $        046D4*xyzi(i+7,j,k)-4.780454397853947D5*xyzi(i+5,j+8,k)+1.
     $        070821785119284D6*xyzi(i+5,j+6,k)-8.380344405281354D5*xyzi
     $        (i+5,j+4,k)+2.660426795327414D5*xyzi(i+5,j+2,k)-2.80044925
     $        8239383D4*xyzi(i+5,j,k)-1.912181759141579D5*xyzi(i+3,j+10,
     $        k)+5.35410892559642D5*xyzi(i+3,j+8,k)-5.586896270187569D5*
     $        xyzi(i+3,j+6,k)+2.660426795327414D5*xyzi(i+3,j+4,k)-5.6008
     $        98516478766D4*xyzi(i+3,j+2,k)+3.953575423396776D3*xyzi(i+3
     $        ,j,k)-3.186969598569298D4*xyzi(i+1,j+12,k)+1.0708217851192
     $        84D5*xyzi(i+1,j+10,k)-1.396724067546892D5*xyzi(i+1,j+8,k)+
     $        8.868089317758046D4*xyzi(i+1,j+6,k)-2.800449258239383D4*xy
     $        zi(i+1,j+4,k)+3.953575423396776D3*xyzi(i+1,j+2,k)-1.757144
     $        632620789D2*xyzi(i+1,j,k))+angi
            angi=yrk(186)*(2.36351991153295D4*xyzi(i+13,j,k)+4.727039823
     $        0659D4*xyzi(i+11,j+2,k)-7.56326371690544D4*xyzi(i+11,j,k)-
     $        1.181759955766475D5*xyzi(i+9,j+4,k)-7.56326371690544D4*xyz
     $        i(i+9,j+2,k)+9.2074514814501D4*xyzi(i+9,j,k)-4.72703982306
     $        59D5*xyzi(i+7,j+6,k)+4.537958230143264D5*xyzi(i+7,j+4,k)-5
     $        .261400846542915D4*xyzi(i+7,j,k)-5.908799778832375D5*xyzi(
     $        i+5,j+8,k)+1.058856920366762D6*xyzi(i+5,j+6,k)-5.524470888
     $        87006D5*xyzi(i+5,j+4,k)+5.261400846542915D4*xyzi(i+5,j+2,k
     $        )+1.384579170142872D4*xyzi(i+5,j,k)-3.30892787614613D5*xyz
     $        i(i+3,j+10,k)+8.319590088595984D5*xyzi(i+3,j+8,k)-7.365961
     $        18516008D5*xyzi(i+3,j+6,k)+2.630700423271457D5*xyzi(i+3,j+
     $        4,k)-2.769158340285745D4*xyzi(i+3,j+2,k)-1.303133336605056
     $        D3*xyzi(i+3,j,k)-7.09055973459885D4*xyzi(i+1,j+12,k)+2.268
     $        979115071632D5*xyzi(i+1,j+10,k)-2.76223544443503D5*xyzi(i+
     $        1,j+8,k)+1.578420253962874D5*xyzi(i+1,j+6,k)-4.15373751042
     $        8617D4*xyzi(i+1,j+4,k)+3.909400009815169D3*xyzi(i+1,j+2,k)
     $        )+angi
            angi=yrk(188)*(-1.281798641180881D4*xyzi(i+13,j,k)+7.6907918
     $        47085288D4*xyzi(i+11,j+2,k)+3.691580086600938D4*xyzi(i+11,
     $        j,k)+3.717216059424556D5*xyzi(i+9,j+4,k)-2.584106060620657
     $        D5*xyzi(i+9,j+2,k)-3.852083568627066D4*xyzi(i+9,j,k)+4.614
     $        475108251173D5*xyzi(i+7,j+6,k)-8.121476190522064D5*xyzi(i+
     $        7,j+4,k)+3.081666854901653D5*xyzi(i+7,j+2,k)+1.71203714161
     $        2029D4*xyzi(i+7,j,k)+1.153618777062793D5*xyzi(i+5,j+8,k)-5
     $        .168212121241314D5*xyzi(i+5,j+6,k)+5.392916996077893D5*xyz
     $        i(i+5,j+4,k)-1.540833427450826D5*xyzi(i+5,j+2,k)-2.7032165
     $        39387415D3*xyzi(i+5,j,k)-1.281798641180881D5*xyzi(i+3,j+10
     $        ,k)+1.845790043300469D5*xyzi(i+3,j+8,k)-8.560185708060147D
     $        4*xyzi(i+3,j+4,k)+2.703216539387415D4*xyzi(i+3,j+2,k)-6.40
     $        8993205904407D4*xyzi(i+1,j+12,k)+1.845790043300469D5*xyzi(
     $        i+1,j+10,k)-1.926041784313533D5*xyzi(i+1,j+8,k)+8.56018570
     $        8060147D4*xyzi(i+1,j+6,k)-1.351608269693707D4*xyzi(i+1,j+4
     $        ,k))+angi
            angi=yrk(190)*(4.920644864827347D3*xyzi(i+13,j,k)-8.85716075
     $        6689225D4*xyzi(i+11,j+2,k)-1.180954767558563D4*xyzi(i+11,j
     $        ,k)-1.230161216206837D5*xyzi(i+9,j+4,k)+2.24381405836127D5
     $        *xyzi(i+9,j+2,k)+9.242254702632235D3*xyzi(i+9,j,k)+1.77143
     $        2151337845D5*xyzi(i+7,j+6,k)+7.08572860535138D4*xyzi(i+7,j
     $        +4,k)-1.848450940526447D5*xyzi(i+7,j+2,k)-2.34723928955739
     $        3D3*xyzi(i+7,j,k)+3.100006264841229D5*xyzi(i+5,j+8,k)-4.96
     $        0010023745966D5*xyzi(i+5,j+6,k)+1.293915658368513D5*xyzi(i
     $        +5,j+4,k)+4.929202508070525D4*xyzi(i+5,j+2,k)+6.8889028107
     $        58286D4*xyzi(i+3,j+10,k)-2.480005011872983D5*xyzi(i+3,j+8,
     $        k)+2.587831316737026D5*xyzi(i+3,j+6,k)-8.215337513450875D4
     $        *xyzi(i+3,j+4,k)-3.444451405379143D4*xyzi(i+1,j+12,k)+8.26
     $        6683372909943D4*xyzi(i+1,j+10,k)-6.469578291842564D4*xyzi(
     $        i+1,j+8,k)+1.643067502690175D4*xyzi(i+1,j+6,k))+angi
            angi=yrk(192)*(-1.253896417710616D3*xyzi(i+13,j,k)+4.2632478
     $        20216095D4*xyzi(i+11,j+2,k)+2.206857695170684D3*xyzi(i+11,
     $        j,k)-6.896430297408388D4*xyzi(i+9,j+4,k)-7.724001933097395
     $        D4*xyzi(i+9,j+2,k)-9.595033457263844D2*xyzi(i+9,j,k)-1.655
     $        143271378013D5*xyzi(i+7,j+6,k)+1.986171925653616D5*xyzi(i+
     $        7,j+4,k)+3.454212044614984D4*xyzi(i+7,j+2,k)+4.13785817844
     $        5033D4*xyzi(i+5,j+8,k)+9.268802319716874D4*xyzi(i+5,j+6,k)
     $        -1.208974215615244D5*xyzi(i+5,j+4,k)+8.275716356890066D4*x
     $        yzi(i+3,j+10,k)-1.655143271378013D5*xyzi(i+3,j+8,k)+8.0598
     $        28104101629D4*xyzi(i+3,j+6,k)-1.128506775939554D4*xyzi(i+1
     $        ,j+12,k)+1.986171925653616D4*xyzi(i+1,j+10,k)-8.6355301115
     $        3746D3*xyzi(i+1,j+8,k))+angi
            angi=yrk(194)*(1.848769406428712D2*xyzi(i+13,j,k)-9.98335479
     $        4715047D3*xyzi(i+11,j+2,k)-1.774818630171564D2*xyzi(i+11,j
     $        ,k)+5.084115867678959D4*xyzi(i+9,j+4,k)+9.761502465943601D
     $        3*xyzi(i+9,j+2,k)-2.4403756164859D4*xyzi(i+7,j+6,k)-5.8569
     $        01479566161D4*xyzi(i+7,j+4,k)-5.490845137093276D4*xyzi(i+5
     $        ,j+8,k)+8.199662071392625D4*xyzi(i+5,j+6,k)+2.847104885900
     $        217D4*xyzi(i+3,j+10,k)-2.92845073978308D4*xyzi(i+3,j+8,k)-
     $        2.033646347071584D3*xyzi(i+1,j+12,k)+1.95230049318872D3*xy
     $        zi(i+1,j+10,k))+angi
            angi=yrk(196)*(-1.025512752521207D1*xyzi(i+13,j,k)+7.9989994
     $        69665415D2*xyzi(i+11,j+2,k)-7.332416180526631D3*xyzi(i+9,j
     $        +4,k)+1.759779883326391D4*xyzi(i+7,j+6,k)-1.31983491249479
     $        4D4*xyzi(i+5,j+8,k)+2.932966472210652D3*xyzi(i+3,j+10,k)-1
     $        .333166578277569D2*xyzi(i+1,j+12,k))+angi
          endif
          if((ieven.and..not.jeven.and.keven))then
            angi=yrk(170)*(-1.333166578277569D2*xyzi(i+12,j+1,k)+2.93296
     $        6472210652D3*xyzi(i+10,j+3,k)-1.319834912494794D4*xyzi(i+8
     $        ,j+5,k)+1.759779883326391D4*xyzi(i+6,j+7,k)-7.332416180526
     $        631D3*xyzi(i+4,j+9,k)+7.998999469665415D2*xyzi(i+2,j+11,k)
     $        -1.025512752521207D1*xyzi(i,j+13,k))+angi
            angi=yrk(172)*(2.033646347071584D3*xyzi(i+12,j+1,k)-2.847104
     $        885900217D4*xyzi(i+10,j+3,k)-1.95230049318872D3*xyzi(i+10,
     $        j+1,k)+5.490845137093276D4*xyzi(i+8,j+5,k)+2.9284507397830
     $        8D4*xyzi(i+8,j+3,k)+2.4403756164859D4*xyzi(i+6,j+7,k)-8.19
     $        9662071392625D4*xyzi(i+6,j+5,k)-5.084115867678959D4*xyzi(i
     $        +4,j+9,k)+5.856901479566161D4*xyzi(i+4,j+7,k)+9.9833547947
     $        15047D3*xyzi(i+2,j+11,k)-9.761502465943601D3*xyzi(i+2,j+9,
     $        k)-1.848769406428712D2*xyzi(i,j+13,k)+1.774818630171564D2*
     $        xyzi(i,j+11,k))+angi
            angi=yrk(174)*(-1.128506775939554D4*xyzi(i+12,j+1,k)+8.27571
     $        6356890066D4*xyzi(i+10,j+3,k)+1.986171925653616D4*xyzi(i+1
     $        0,j+1,k)+4.137858178445033D4*xyzi(i+8,j+5,k)-1.65514327137
     $        8013D5*xyzi(i+8,j+3,k)-8.63553011153746D3*xyzi(i+8,j+1,k)-
     $        1.655143271378013D5*xyzi(i+6,j+7,k)+9.268802319716874D4*xy
     $        zi(i+6,j+5,k)+8.059828104101629D4*xyzi(i+6,j+3,k)-6.896430
     $        297408388D4*xyzi(i+4,j+9,k)+1.986171925653616D5*xyzi(i+4,j
     $        +7,k)-1.208974215615244D5*xyzi(i+4,j+5,k)+4.26324782021609
     $        5D4*xyzi(i+2,j+11,k)-7.724001933097395D4*xyzi(i+2,j+9,k)+3
     $        .454212044614984D4*xyzi(i+2,j+7,k)-1.253896417710616D3*xyz
     $        i(i,j+13,k)+2.206857695170684D3*xyzi(i,j+11,k)-9.595033457
     $        263844D2*xyzi(i,j+9,k))+angi
            angi=yrk(176)*(3.444451405379143D4*xyzi(i+12,j+1,k)-6.888902
     $        810758286D4*xyzi(i+10,j+3,k)-8.266683372909943D4*xyzi(i+10
     $        ,j+1,k)-3.100006264841229D5*xyzi(i+8,j+5,k)+2.480005011872
     $        983D5*xyzi(i+8,j+3,k)+6.469578291842564D4*xyzi(i+8,j+1,k)-
     $        1.771432151337845D5*xyzi(i+6,j+7,k)+4.960010023745966D5*xy
     $        zi(i+6,j+5,k)-2.587831316737026D5*xyzi(i+6,j+3,k)-1.643067
     $        502690175D4*xyzi(i+6,j+1,k)+1.230161216206837D5*xyzi(i+4,j
     $        +9,k)-7.08572860535138D4*xyzi(i+4,j+7,k)-1.293915658368513
     $        D5*xyzi(i+4,j+5,k)+8.215337513450875D4*xyzi(i+4,j+3,k)+8.8
     $        57160756689225D4*xyzi(i+2,j+11,k)-2.24381405836127D5*xyzi(
     $        i+2,j+9,k)+1.848450940526447D5*xyzi(i+2,j+7,k)-4.929202508
     $        070525D4*xyzi(i+2,j+5,k)-4.920644864827347D3*xyzi(i,j+13,k
     $        )+1.180954767558563D4*xyzi(i,j+11,k)-9.242254702632235D3*x
     $        yzi(i,j+9,k)+2.347239289557393D3*xyzi(i,j+7,k))+angi
            angi=yrk(178)*(-6.408993205904407D4*xyzi(i+12,j+1,k)-1.28179
     $        8641180881D5*xyzi(i+10,j+3,k)+1.845790043300469D5*xyzi(i+1
     $        0,j+1,k)+1.153618777062793D5*xyzi(i+8,j+5,k)+1.84579004330
     $        0469D5*xyzi(i+8,j+3,k)-1.926041784313533D5*xyzi(i+8,j+1,k)
     $        +4.614475108251173D5*xyzi(i+6,j+7,k)-5.168212121241314D5*x
     $        yzi(i+6,j+5,k)+8.560185708060147D4*xyzi(i+6,j+1,k)+3.71721
     $        6059424556D5*xyzi(i+4,j+9,k)-8.121476190522064D5*xyzi(i+4,
     $        j+7,k)+5.392916996077893D5*xyzi(i+4,j+5,k)-8.5601857080601
     $        47D4*xyzi(i+4,j+3,k)-1.351608269693707D4*xyzi(i+4,j+1,k)+7
     $        .690791847085288D4*xyzi(i+2,j+11,k)-2.584106060620657D5*xy
     $        zi(i+2,j+9,k)+3.081666854901653D5*xyzi(i+2,j+7,k)-1.540833
     $        427450826D5*xyzi(i+2,j+5,k)+2.703216539387415D4*xyzi(i+2,j
     $        +3,k)-1.281798641180881D4*xyzi(i,j+13,k)+3.691580086600938
     $        D4*xyzi(i,j+11,k)-3.852083568627066D4*xyzi(i,j+9,k)+1.7120
     $        37141612029D4*xyzi(i,j+7,k)-2.703216539387415D3*xyzi(i,j+5
     $        ,k))+angi
            angi=yrk(180)*(7.09055973459885D4*xyzi(i+12,j+1,k)+3.3089278
     $        7614613D5*xyzi(i+10,j+3,k)-2.268979115071632D5*xyzi(i+10,j
     $        +1,k)+5.908799778832375D5*xyzi(i+8,j+5,k)-8.31959008859598
     $        4D5*xyzi(i+8,j+3,k)+2.76223544443503D5*xyzi(i+8,j+1,k)+4.7
     $        270398230659D5*xyzi(i+6,j+7,k)-1.058856920366762D6*xyzi(i+
     $        6,j+5,k)+7.36596118516008D5*xyzi(i+6,j+3,k)-1.578420253962
     $        874D5*xyzi(i+6,j+1,k)+1.181759955766475D5*xyzi(i+4,j+9,k)-
     $        4.537958230143264D5*xyzi(i+4,j+7,k)+5.52447088887006D5*xyz
     $        i(i+4,j+5,k)-2.630700423271457D5*xyzi(i+4,j+3,k)+4.1537375
     $        10428617D4*xyzi(i+4,j+1,k)-4.7270398230659D4*xyzi(i+2,j+11
     $        ,k)+7.56326371690544D4*xyzi(i+2,j+9,k)-5.261400846542915D4
     $        *xyzi(i+2,j+5,k)+2.769158340285745D4*xyzi(i+2,j+3,k)-3.909
     $        400009815169D3*xyzi(i+2,j+1,k)-2.36351991153295D4*xyzi(i,j
     $        +13,k)+7.56326371690544D4*xyzi(i,j+11,k)-9.2074514814501D4
     $        *xyzi(i,j+9,k)+5.261400846542915D4*xyzi(i,j+7,k)-1.3845791
     $        70142872D4*xyzi(i,j+5,k)+1.303133336605056D3*xyzi(i,j+3,k)
     $        )+angi
            angi=yrk(182)*(-3.186969598569298D4*xyzi(i+12,j+1,k)-1.91218
     $        1759141579D5*xyzi(i+10,j+3,k)+1.070821785119284D5*xyzi(i+1
     $        0,j+1,k)-4.780454397853947D5*xyzi(i+8,j+5,k)+5.35410892559
     $        642D5*xyzi(i+8,j+3,k)-1.396724067546892D5*xyzi(i+8,j+1,k)-
     $        6.373939197138596D5*xyzi(i+6,j+7,k)+1.070821785119284D6*xy
     $        zi(i+6,j+5,k)-5.586896270187569D5*xyzi(i+6,j+3,k)+8.868089
     $        317758046D4*xyzi(i+6,j+1,k)-4.780454397853947D5*xyzi(i+4,j
     $        +9,k)+1.070821785119284D6*xyzi(i+4,j+7,k)-8.38034440528135
     $        4D5*xyzi(i+4,j+5,k)+2.660426795327414D5*xyzi(i+4,j+3,k)-2.
     $        800449258239383D4*xyzi(i+4,j+1,k)-1.912181759141579D5*xyzi
     $        (i+2,j+11,k)+5.35410892559642D5*xyzi(i+2,j+9,k)-5.58689627
     $        0187569D5*xyzi(i+2,j+7,k)+2.660426795327414D5*xyzi(i+2,j+5
     $        ,k)-5.600898516478766D4*xyzi(i+2,j+3,k)+3.953575423396776D
     $        3*xyzi(i+2,j+1,k)-3.186969598569298D4*xyzi(i,j+13,k)+1.070
     $        821785119284D5*xyzi(i,j+11,k)-1.396724067546892D5*xyzi(i,j
     $        +9,k)+8.868089317758046D4*xyzi(i,j+7,k)-2.800449258239383D
     $        4*xyzi(i,j+5,k)+3.953575423396776D3*xyzi(i,j+3,k)-1.757144
     $        632620789D2*xyzi(i,j+1,k))+angi
          endif
          if((ieven.and.jeven.and..not.keven))then
            angi=yrk(183)*(2.338596333691754D4*xyzi(i,j,k+13)-7.29642056
     $        1118272D4*xyzi(i,j,k+11)+8.72398110568489D4*xyzi(i,j,k+9)-
     $        4.985132060391366D4*xyzi(i,j,k+7)+1.377470700897614D4*xyzi
     $        (i,j,k+5)-1.620553765761899D3*xyzi(i,j,k+3)+5.401845885872
     $        997D1*xyzi(i,j,k+1))+angi
            angi=yrk(185)*(-2.850512265850467D4*xyzi(i+12,j,k+1)-1.14020
     $        4906340187D5*xyzi(i+10,j+2,k+1)+7.981434344381306D4*xyzi(i
     $        +10,j,k+1)-1.425256132925233D5*xyzi(i+8,j+4,k+1)+2.3944303
     $        03314392D5*xyzi(i+8,j+2,k+1)-8.328453228919624D4*xyzi(i+8,
     $        j,k+1)+1.596286868876261D5*xyzi(i+6,j+4,k+1)-1.66569064578
     $        3925D5*xyzi(i+6,j+2,k+1)+3.965930109009345D4*xyzi(i+6,j,k+
     $        1)+1.425256132925233D5*xyzi(i+4,j+8,k+1)-1.596286868876261
     $        D5*xyzi(i+4,j+6,k+1)+3.965930109009345D4*xyzi(i+4,j+2,k+1)
     $        -8.349326545282831D3*xyzi(i+4,j,k+1)+1.140204906340187D5*x
     $        yzi(i+2,j+10,k+1)-2.394430303314392D5*xyzi(i+2,j+8,k+1)+1.
     $        665690645783925D5*xyzi(i+2,j+6,k+1)-3.965930109009345D4*xy
     $        zi(i+2,j+4,k+1)+5.893642267258469D2*xyzi(i+2,j,k+1)+2.8505
     $        12265850467D4*xyzi(i,j+12,k+1)-7.981434344381306D4*xyzi(i,
     $        j+10,k+1)+8.328453228919624D4*xyzi(i,j+8,k+1)-3.9659301090
     $        09345D4*xyzi(i,j+6,k+1)+8.349326545282831D3*xyzi(i,j+4,k+1
     $        )-5.893642267258469D2*xyzi(i,j+2,k+1))+angi
            angi=yrk(187)*(1.812737022589407D4*xyzi(i+12,j,k+1)-3.625474
     $        045178814D4*xyzi(i+10,j+2,k+1)-4.640606777828882D4*xyzi(i+
     $        10,j,k+1)-3.081652938401992D5*xyzi(i+8,j+4,k+1)+1.39218203
     $        3348665D5*xyzi(i+8,j+2,k+1)+4.237075753669848D4*xyzi(i+8,j
     $        ,k+1)-5.075663663250339D5*xyzi(i+6,j+6,k+1)+6.496849488960
     $        434D5*xyzi(i+6,j+4,k+1)-1.694830301467939D5*xyzi(i+6,j+2,k
     $        +1)-1.614124096636133D4*xyzi(i+6,j,k+1)-3.081652938401992D
     $        5*xyzi(i+4,j+8,k+1)+6.496849488960434D5*xyzi(i+4,j+6,k+1)-
     $        4.237075753669848D5*xyzi(i+4,j+4,k+1)+8.070620483180664D4*
     $        xyzi(i+4,j+2,k+1)+2.123847495573859D3*xyzi(i+4,j,k+1)-3.62
     $        5474045178814D4*xyzi(i+2,j+10,k+1)+1.392182033348665D5*xyz
     $        i(i+2,j+8,k+1)-1.694830301467939D5*xyzi(i+2,j+6,k+1)+8.070
     $        620483180664D4*xyzi(i+2,j+4,k+1)-1.274308497344315D4*xyzi(
     $        i+2,j+2,k+1)+1.812737022589407D4*xyzi(i,j+12,k+1)-4.640606
     $        777828882D4*xyzi(i,j+10,k+1)+4.237075753669848D4*xyzi(i,j+
     $        8,k+1)-1.614124096636133D4*xyzi(i,j+6,k+1)+2.1238474955738
     $        59D3*xyzi(i,j+4,k+1))+angi
            angi=yrk(189)*(-8.317407887033718D3*xyzi(i+12,j,k+1)+9.98088
     $        9464440461D4*xyzi(i+10,j+2,k+1)+1.796560103599283D4*xyzi(i
     $        +10,j,k+1)+2.245700129499104D5*xyzi(i+8,j+4,k+1)-2.3355281
     $        34679068D5*xyzi(i+8,j+2,k+1)-1.249780941634284D4*xyzi(i+8,
     $        j,k+1)-2.515184145038996D5*xyzi(i+6,j+4,k+1)+1.74969331828
     $        7997D5*xyzi(i+6,j+2,k+1)+2.777290981409519D3*xyzi(i+6,j,k+
     $        1)-2.245700129499104D5*xyzi(i+4,j+8,k+1)+2.515184145038996
     $        D5*xyzi(i+4,j+6,k+1)-4.165936472114279D4*xyzi(i+4,j+2,k+1)
     $        -9.980889464440461D4*xyzi(i+2,j+10,k+1)+2.335528134679068D
     $        5*xyzi(i+2,j+8,k+1)-1.749693318287997D5*xyzi(i+2,j+6,k+1)+
     $        4.165936472114279D4*xyzi(i+2,j+4,k+1)+8.317407887033718D3*
     $        xyzi(i,j+12,k+1)-1.796560103599283D4*xyzi(i,j+10,k+1)+1.24
     $        9780941634284D4*xyzi(i,j+8,k+1)-2.777290981409519D3*xyzi(i
     $        ,j+6,k+1))+angi
            angi=yrk(191)*(2.630195315167501D3*xyzi(i+12,j,k+1)-6.838507
     $        819435502D4*xyzi(i+10,j+2,k+1)-4.208312504268001D3*xyzi(i+
     $        10,j,k+1)+3.945292972751251D4*xyzi(i+8,j+4,k+1)+1.13624437
     $        615236D5*xyzi(i+8,j+2,k+1)+1.646730979930957D3*xyzi(i+8,j,
     $        k+1)+2.209364064740701D5*xyzi(i+6,j+6,k+1)-1.7674912517925
     $        61D5*xyzi(i+6,j+4,k+1)-4.61084674380668D4*xyzi(i+6,j+2,k+1
     $        )+3.945292972751251D4*xyzi(i+4,j+8,k+1)-1.767491251792561D
     $        5*xyzi(i+4,j+6,k+1)+1.15271168595167D5*xyzi(i+4,j+4,k+1)-6
     $        .838507819435502D4*xyzi(i+2,j+10,k+1)+1.13624437615236D5*x
     $        yzi(i+2,j+8,k+1)-4.61084674380668D4*xyzi(i+2,j+6,k+1)+2.63
     $        0195315167501D3*xyzi(i,j+12,k+1)-4.208312504268001D3*xyzi(
     $        i,j+10,k+1)+1.646730979930957D3*xyzi(i,j+8,k+1))+angi
            angi=yrk(193)*(-5.229109536543883D2*xyzi(i+12,j,k+1)+2.30080
     $        8196079309D4*xyzi(i+10,j+2,k+1)+4.601616392158617D2*xyzi(i
     $        +10,j,k+1)-8.628030735297407D4*xyzi(i+8,j+4,k+1)-2.0707273
     $        76471378D4*xyzi(i+8,j+2,k+1)+9.663394423533096D4*xyzi(i+6,
     $        j+4,k+1)+8.628030735297407D4*xyzi(i+4,j+8,k+1)-9.663394423
     $        533096D4*xyzi(i+4,j+6,k+1)-2.300808196079309D4*xyzi(i+2,j+
     $        10,k+1)+2.070727376471378D4*xyzi(i+2,j+8,k+1)+5.2291095365
     $        43883D2*xyzi(i,j+12,k+1)-4.601616392158617D2*xyzi(i,j+10,k
     $        +1))+angi
            angi=yrk(195)*(5.229109536543883D1*xyzi(i+12,j,k+1)-3.451212
     $        294118963D3*xyzi(i+10,j+2,k+1)+2.588409220589222D4*xyzi(i+
     $        8,j+4,k+1)-4.831697211766548D4*xyzi(i+6,j+6,k+1)+2.5884092
     $        20589222D4*xyzi(i+4,j+8,k+1)-3.451212294118963D3*xyzi(i+2,
     $        j+10,k+1)+5.229109536543883D1*xyzi(i,j+12,k+1))+angi
          endif
          if((.not.ieven.and..not.jeven.and..not.keven))then
            angi=yrk(171)*(6.27493144385266D2*xyzi(i+11,j+1,k+1)-1.15040
     $        4098039654D4*xyzi(i+9,j+3,k+1)+4.141454752942756D4*xyzi(i+
     $        7,j+5,k+1)-4.141454752942756D4*xyzi(i+5,j+7,k+1)+1.1504040
     $        98039654D4*xyzi(i+3,j+9,k+1)-6.27493144385266D2*xyzi(i+1,j
     $        +11,k+1))+angi
            angi=yrk(173)*(-5.229109536543883D3*xyzi(i+11,j+1,k+1)+5.752
     $        020490198272D4*xyzi(i+9,j+3,k+1)+4.601616392158617D3*xyzi(
     $        i+9,j+1,k+1)-6.902424588237926D4*xyzi(i+7,j+5,k+1)-5.52193
     $        9670590341D4*xyzi(i+7,j+3,k+1)-6.902424588237926D4*xyzi(i+
     $        5,j+7,k+1)+1.159607330823972D5*xyzi(i+5,j+5,k+1)+5.7520204
     $        90198272D4*xyzi(i+3,j+9,k+1)-5.521939670590341D4*xyzi(i+3,
     $        j+7,k+1)-5.229109536543883D3*xyzi(i+1,j+11,k+1)+4.60161639
     $        2158617D3*xyzi(i+1,j+9,k+1))+angi
            angi=yrk(175)*(2.104156252134001D4*xyzi(i+11,j+1,k+1)-1.0520
     $        78126067D5*xyzi(i+9,j+3,k+1)-3.366650003414401D4*xyzi(i+9,
     $        j+1,k+1)-1.2624937512804D5*xyzi(i+7,j+5,k+1)+2.01999000204
     $        8641D5*xyzi(i+7,j+3,k+1)+1.317384783944766D4*xyzi(i+7,j+1,
     $        k+1)+1.2624937512804D5*xyzi(i+5,j+7,k+1)-9.221693487613359
     $        D4*xyzi(i+5,j+3,k+1)+1.052078126067D5*xyzi(i+3,j+9,k+1)-2.
     $        019990002048641D5*xyzi(i+3,j+7,k+1)+9.221693487613359D4*xy
     $        zi(i+3,j+5,k+1)-2.104156252134001D4*xyzi(i+1,j+11,k+1)+3.3
     $        66650003414401D4*xyzi(i+1,j+9,k+1)-1.317384783944766D4*xyz
     $        i(i+1,j+7,k+1))+angi
            angi=yrk(177)*(-4.990444732220231D4*xyzi(i+11,j+1,k+1)+1.663
     $        481577406744D4*xyzi(i+9,j+3,k+1)+1.07793606215957D5*xyzi(i
     $        +9,j+1,k+1)+2.994266839332138D5*xyzi(i+7,j+5,k+1)-1.437248
     $        082879426D5*xyzi(i+7,j+3,k+1)-7.498685649805703D4*xyzi(i+7
     $        ,j+1,k+1)+2.994266839332138D5*xyzi(i+5,j+7,k+1)-5.03036829
     $        0077993D5*xyzi(i+5,j+5,k+1)+1.749693318287997D5*xyzi(i+5,j
     $        +3,k+1)+1.666374588845712D4*xyzi(i+5,j+1,k+1)+1.6634815774
     $        06744D4*xyzi(i+3,j+9,k+1)-1.437248082879426D5*xyzi(i+3,j+7
     $        ,k+1)+1.749693318287997D5*xyzi(i+3,j+5,k+1)-5.554581962819
     $        039D4*xyzi(i+3,j+3,k+1)-4.990444732220231D4*xyzi(i+1,j+11,
     $        k+1)+1.07793606215957D5*xyzi(i+1,j+9,k+1)-7.49868564980570
     $        3D4*xyzi(i+1,j+7,k+1)+1.666374588845712D4*xyzi(i+1,j+5,k+1
     $        ))+angi
            angi=yrk(179)*(7.250948090357628D4*xyzi(i+11,j+1,k+1)+2.1752
     $        84427107288D5*xyzi(i+9,j+3,k+1)-1.856242711131553D5*xyzi(i
     $        +9,j+1,k+1)+1.450189618071526D5*xyzi(i+7,j+5,k+1)-3.712485
     $        422263105D5*xyzi(i+7,j+3,k+1)+1.694830301467939D5*xyzi(i+7
     $        ,j+1,k+1)-1.450189618071526D5*xyzi(i+5,j+7,k+1)+1.69483030
     $        1467939D5*xyzi(i+5,j+3,k+1)-6.456496386544531D4*xyzi(i+5,j
     $        +1,k+1)-2.175284427107288D5*xyzi(i+3,j+9,k+1)+3.7124854222
     $        63105D5*xyzi(i+3,j+7,k+1)-1.694830301467939D5*xyzi(i+3,j+5
     $        ,k+1)+8.495389982295436D3*xyzi(i+3,j+1,k+1)-7.250948090357
     $        628D4*xyzi(i+1,j+11,k+1)+1.856242711131553D5*xyzi(i+1,j+9,
     $        k+1)-1.694830301467939D5*xyzi(i+1,j+7,k+1)+6.4564963865445
     $        31D4*xyzi(i+1,j+5,k+1)-8.495389982295436D3*xyzi(i+1,j+3,k+
     $        1))+angi
            angi=yrk(181)*(-5.701024531700933D4*xyzi(i+11,j+1,k+1)-2.850
     $        512265850467D5*xyzi(i+9,j+3,k+1)+1.596286868876261D5*xyzi(
     $        i+9,j+1,k+1)-5.701024531700933D5*xyzi(i+7,j+5,k+1)+6.38514
     $        7475505045D5*xyzi(i+7,j+3,k+1)-1.665690645783925D5*xyzi(i+
     $        7,j+1,k+1)-5.701024531700933D5*xyzi(i+5,j+7,k+1)+9.5777212
     $        13257568D5*xyzi(i+5,j+5,k+1)-4.997071937351775D5*xyzi(i+5,
     $        j+3,k+1)+7.93186021801869D4*xyzi(i+5,j+1,k+1)-2.8505122658
     $        50467D5*xyzi(i+3,j+9,k+1)+6.385147475505045D5*xyzi(i+3,j+7
     $        ,k+1)-4.997071937351775D5*xyzi(i+3,j+5,k+1)+1.586372043603
     $        738D5*xyzi(i+3,j+3,k+1)-1.669865309056566D4*xyzi(i+3,j+1,k
     $        +1)-5.701024531700933D4*xyzi(i+1,j+11,k+1)+1.5962868688762
     $        61D5*xyzi(i+1,j+9,k+1)-1.665690645783925D5*xyzi(i+1,j+7,k+
     $        1)+7.93186021801869D4*xyzi(i+1,j+5,k+1)-1.669865309056566D
     $        4*xyzi(i+1,j+3,k+1)+1.178728453451694D3*xyzi(i+1,j+1,k+1))
     $        +angi
          endif
          OMEGA=angi
        endif
      endif
c...
      end
c=======================================================================
      subroutine omega0one(i,j,k,lone,yrk,mijk,xyzi,OMEGA)
      implicit none
c...
c...  Max2fort (MM) version 0.1 (September 2004)
c...  translated on 12/9/2004 18:18:40
c...
c...  This subroutine computes type 1 angular integrals over pseudopotentials
c...
c...
c...   omega(i,j,k,l)=
c...
c...                 / /
c...                 [ [                     i   j   k
c...  yr(0,l,thetak) I I yr(0,l,THETA,PHI) xs  ys  zs  dPHi dTHETA
c...                 ] ]
c...                 / /
c...
c...  where yr are real spherical harmonic functions, and xs, ys and zs
c...  are cartesian coordinates costrained over the surface of the unit
c...  sphere:
c...              xs=sin(THETA) cos(PHI)
c...              ys=sin(THETA) sin(PHI)
c...              zs=cos(THETA)
c...
c...  this version is to be used when the angle phik cannot be defined
c...  (both kx and ky are zero), and thus only the real spherical
c...  harmonics with m=0 survive
c...
c...  The integral is computed by expanding the real sperical harmonic
c...  as a polynomial in xs, ys and zs, obtaining integrals of the type:
c...
c...                      / /
c...                 1    [ [   a   b   c         (a-1)!! (b-1)!! (c-1)!!
c...  xyzi(a,b,c)= -----  I I xs  ys  zs  dP dT = -----------------------
c...               4 %pi  ] ]                          (a+b+c+1)!!
c...                      / /
c...
c...  if a, b and c are all even, or zero otherwise
c...
c...  This subroutine computes the integral with l = lone,
c...  with the additional constraint that l-(i+j+k) must be even
c...
c...  **** This version works for l up to 14 ****
c...
c...  McMurchie and Davidson, J. Comp. Phys. 44, 289 (1981),
c...  eqs. (18),(27-30)
c...
c...  David B. Cook, Handbook of Computational Quantum Chemistry,
c...  (Oxford, 1998), page 596-601
c...
c...  Input parameters:
c...
c...  i       exponent of xs
c...  j       exponent of ys
c...  k       exponent of zs
c...  lone    lvalue to compute
c...  yrk     array of real spherical harmonics evaluated at k
c...  mijk    maximum exponent value, for dimensioning
c...  xyzi    table of double factorial products
c...
c...  Output parameter:
c...
c...  omega   array of type one angular integrals
c...
      real*8 ZERO
      parameter (ZERO=0.0D0)
c...
      integer i,j,k,lmax,mijk,lone
      real*8 yrk(*),OMEGA,xyzi(0:mijk,0:mijk,0:mijk)
      logical ieven,jeven,keven
c...
c...
c...  determine whether i,j and k are even or odd
c...
      ieven=mod(i,2).eq.0
      jeven=mod(j,2).eq.0
      keven=mod(k,2).eq.0
c...
c...  compute the upper limit for l
c...  and whether we must do the even or odd series
c...
      lmax=k+j+i
      if(lmax.gt.14)then
        write(6,*)'omega0one: maximum l value exceeded',lmax,14
        call nerror(1,'omega0one','maximum l value exceeded',lmax,14)
      endif
      if(lone.gt.14)then
        call nerror(2,'omega0one','inconsistent call',lone,14)
      endif
c...
c...  zero out omega
c...
      OMEGA=ZERO
c...
c...  integrals survive only if i and j are both even
c...
      if((ieven.and.jeven))then
        if(keven)then
c...
c...  l=0
c...
          if(lone.eq.0)then
            OMEGA=3.544907701811032D0*yrk(1)*xyzi(i,j,k)
          endif
c...
c...  l=2
c...
          if(lone.eq.2)then
            OMEGA=yrk(7)*(1.188998189281803D1*xyzi(i,j,k+2)-3.9633272976
     $        06011D0*xyzi(i,j,k))
          endif
c...
c...  l=4
c...
          if(lone.eq.4)then
            OMEGA=yrk(21)*(4.65269135862698D1*xyzi(i,j,k+4)-3.9880211645
     $        37411D1*xyzi(i,j,k+2)+3.988021164537411D0*xyzi(i,j,k))
          endif
c...
c...  l=6
c...
          if(lone.eq.6)then
            OMEGA=yrk(43)*(1.845306898868157D2*xyzi(i,j,k+6)-2.516327589
     $        365668D2*xyzi(i,j,k+4)+8.387758631218894D1*xyzi(i,j,k+2)-3
     $        .994170776770902D0*xyzi(i,j,k))
          endif
c...
c...  l=8
c...
          if(lone.eq.8)then
            OMEGA=yrk(73)*(7.347980147805839D2*xyzi(i,j,k+8)-1.371622960
     $        923757D3*xyzi(i,j,k+6)+7.91320938994475D2*xyzi(i,j,k+4)-1.
     $        438765343626318D2*xyzi(i,j,k+2)+3.996570398961995D0*xyzi(i
     $        ,j,k))
          endif
c...
c...  l=10
c...
          if(lone.eq.10)then
            OMEGA=yrk(111)*(2.930982152135684D3*xyzi(i,j,k+10)-6.9417998
     $        34005568D3*xyzi(i,j,k+8)+5.716776333886939D3*xyzi(i,j,k+6)
     $        -1.905592111295646D3*xyzi(i,j,k+4)+2.198760128418053D2*xyz
     $        i(i,j,k+2)-3.997745688032824D0*xyzi(i,j,k))
          endif
c...
c...  l=12
c...
          if(lone.eq.12)then
            OMEGA=yrk(157)*(1.170163993078432D4*xyzi(i,j,k+12)-3.3578618
     $        93181587D4*xyzi(i,j,k+10)+3.597709171265986D4*xyzi(i,j,k+8
     $        )-1.767295733253467D4*xyzi(i,j,k+6)+3.898446470412059D3*xy
     $        zi(i,j,k+4)-3.118757176329647D2*xyzi(i,j,k+2)+3.9984066363
     $        20061D0*xyzi(i,j,k))
          endif
c...
c...  l=14
c...
          if(lone.eq.14)then
            OMEGA=yrk(211)*(4.674208812108592D4*xyzi(i,j,k+14)-1.5753814
     $        88525489D5*xyzi(i,j,k+12)+2.079503564853645D5*xyzi(i,j,k+1
     $        0)-1.356197977078464D5*xyzi(i,j,k+8)+4.52065992359488D4*xy
     $        zi(i,j,k+6)-7.137884089886653D3*xyzi(i,j,k+4)+4.1987553469
     $        92149D2*xyzi(i,j,k+2)-3.998814616182999D0*xyzi(i,j,k))
          endif
        else
c...
c...  l=1
c...
          if(lone.eq.1)then
            OMEGA=6.139960247678931D0*yrk(3)*xyzi(i,j,k+1)
          endif
c...
c...  l=3
c...
          if(lone.eq.3)then
            OMEGA=yrk(13)*(2.344736049917376D1*xyzi(i,j,k+3)-1.406841629
     $        950425D1*xyzi(i,j,k+1))
          endif
c...
c...  l=5
c...
          if(lone.eq.5)then
            OMEGA=yrk(31)*(9.258738901136752D1*xyzi(i,j,k+5)-1.028748766
     $        792972D2*xyzi(i,j,k+3)+2.204461643127798D1*xyzi(i,j,k+1))
          endif
c...
c...  l=7
c...
          if(lone.eq.7)then
            OMEGA=yrk(57)*(3.681186927173971D2*xyzi(i,j,k+7)-5.946532728
     $        511799D2*xyzi(i,j,k+5)+2.702969422050818D2*xyzi(i,j,k+3)-3
     $        .003299357834242D1*xyzi(i,j,k+1))
          endif
c...
c...  l=9
c...
          if(lone.eq.9)then
            OMEGA=yrk(91)*(1.467326381829043D3*xyzi(i,j,k+9)-3.107279396
     $        814444D3*xyzi(i,j,k+7)+2.175095577770111D3*xyzi(i,j,k+5)-5
     $        .57716814812849D2*xyzi(i,j,k+3)+3.802614646451243D1*xyzi(i
     $        ,j,k+1))
          endif
c...
c...  l=11
c...
          if(lone.eq.11)then
            OMEGA=yrk(133)*(5.855905424818466D3*xyzi(i,j,k+11)-1.5336895
     $        16023884D4*xyzi(i,j,k+9)+1.452969015180522D4*xyzi(i,j,k+7)
     $        -5.982813591919795D3*xyzi(i,j,k+5)+9.971355986532992D2*xyz
     $        i(i,j,k+3)-4.602164301476766D1*xyzi(i,j,k+1))
          endif
c...
c...  l=13
c...
          if(lone.eq.13)then
            OMEGA=yrk(183)*(2.338596333691754D4*xyzi(i,j,k+13)-7.2964205
     $        61118272D4*xyzi(i,j,k+11)+8.72398110568489D4*xyzi(i,j,k+9)
     $        -4.985132060391366D4*xyzi(i,j,k+7)+1.377470700897614D4*xyz
     $        i(i,j,k+5)-1.620553765761899D3*xyzi(i,j,k+3)+5.40184588587
     $        2997D1*xyzi(i,j,k+1))
          endif
        endif
      endif
c...
      end

cc  ................................................................
cc  DFT now working with general basis
cc  (modified routines  <AOGrad> and <AOVal>)
cc  JB  12 July 1997
cc
cc  F-functions added
cc  JB   4 Nov 1999
cc .................................................................
      SUBROUTINE AODer3(NP,     ncs,    ncf,    XGRID,
     $                  xnuc,   bl,     basdat, inx,    thrsh,
     $                  ExpMIN, VAO,    VAOX,   VAOXX,  VAOXXX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C     MM (02/27/2003)  Generated with maxima (file aoder.mc)
C
C  calculates vector of basis function (AO) values and their
C  first, second, and third derivatives at the current grid point
C  the structure of the subroutine is derived from aoder2. The
C  expressions for the derivatives have been generated with the
C  maxima program.
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  ncs     -  number of shells
C  ncf     -  number of basis functions
C  XGRID   -  grid points
C  xnuc    -  nuclear coordinates
C  bl      -  array of precomputed normalization factors
C             (from routine <AOInit>)
C  basdat  -  basis set data for TEXAS
C  inx     -  more basis set data
C  thrsh   -  exponent threshold for neglecting contribution
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  VAO     -  on exit contains basis function values
C  VAOX    -  on exit contains basis function derivatives
C  VAOXX   -  on exit contains basis function 2nd derivatives
C  VAOXXX  -  on exit contains basis function 3rd derivatives
C
C
      DIMENSION xnuc(3,*),bl(*),basdat(13,*),inx(12,*),ExpMIN(*),
     $          VAO(ncf,NP),VAOX(3,ncf,NP),VAOXX(6,ncf,NP),
     $          VAOXXX(10,ncf,NP),XGRID(3,NP)
      dimension x0(3)
      parameter (Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0,Three=3.0d0)
      parameter (Four=4.0d0,Five=5.0d0,Six=6.0d0,Seven=7.0d0)
      parameter (Eight=8.0d0,Nine=9.0d0,Ten=10.0d0,Eleven=11.0d0)
      parameter (Twelve=12.0d0)
      parameter (TSeven=17.0d0,Teight=18.0d0,Twenty=20.0d0)
      parameter (TwFour=24.0d0,TwSeven=27.0d0,ThSix=36.0d0)
      parameter (FEight=48.0d0)
      parameter (PI=3.14159 26535 89793d0)
C
C
C  get some values
C
      sq3i = One/SQRT(Three)
      sq10 = SQRT(Ten)
c
      DO 950 IP=1,NP
        xp = XGRID(1,IP)
        yp = XGRID(2,IP)
        zp = XGRID(3,IP)
        do 900 ics=1,ncs
c len=number of ang. mom. components, e.g. 3 for P, 4 for L, 5 for D;
c ia=beginning, ie=ending of the contraction
c iat=atom number, xnuc=coordinates of the center(nucleus)
c   ityp: 1=S, 2=P, 3=SP(L), 4=5D, 5=6D, 6=7F, 7=10F
          iat=inx(2,ics)
          x=xp-xnuc(1,iat)
          y=yp-xnuc(2,iat)
          z=zp-xnuc(3,iat)
c
c  precomputes the first three powers of x, y, and z
c
          x2 = x**2
          x3 = x*x2
          Y2 = y**2
          Y3 = y*Y2
          z2 = z**2
          z3 = z*z2
          r02 = z2+Y2+x2
cc  ...............................................................
cc    global shell exponent test                ! jb march 97
cc  ...............................................................
          rrexp = r02*ExpMIN(ics)
c         write(6,*) ' ics is:',ics
c         write(6,*) ' rrexp is:',rrexp,' thrsh is:',thrsh
          If(rrexp.GT.thrsh) go to 900
cc  ...............................................................
          ityp=inx(12,ics)                   ! shell type
          ia=inx(1,ics)+1                    ! start of contraction
          ie=inx(5,ics)                      ! end of contraction
          ngr = inx(4,ics)                   ! no. of general contractions
          ifun=inx(11,ics)
c -- general contraction loop
          DO 800 igr=0,ngr
            sa0ex1 = ZERO
            sa1ex1 = ZERO
            sa2ex1 = ZERO
            sa3ex1 = ZERO
            if(ityp.eq.3) then
              sa0ex2 = ZERO
              sa1ex2 = ZERO
              sa2ex2 = ZERO
              sa3ex2 = ZERO
            end if
c -- contraction loop
            do 700 ish=ia,ie
c  this the Gaussian exponent times the distance squared, a*r02
              aexp=basdat(1,ish)         ! gaussian exponent
              ar02=aexp*r02
c  ..........................................................
c  do not calculate very small values
              If(ar02.GT.thrsh) go to 700
c  ..........................................................
c  bl(ish) contains the precomputed value (2*aexp/pi)**0.75
c   aexp is the exponent of the Gaussian
              xx=bl(ish)*exp(-ar02)
              a0ex1=basdat(2+igr,ish)*xx
              a1ex1 = a0ex1*aexp
              a2ex1 = a1ex1*aexp
              a3ex1 = a2ex1*aexp
              sa0ex1 = sa0ex1+a0ex1
              sa1ex1 = sa1ex1+a1ex1
              sa2ex1 = sa2ex1+a2ex1
              sa3ex1 = sa3ex1+a3ex1
c.................
              if(ityp.eq.3) then
                a0ex2=basdat(3,ish)*xx
                a1ex2 = a0ex2*aexp
                a2ex2 = a1ex2*aexp
                a3ex2 = a2ex2*aexp
                sa0ex2 = sa0ex2+a0ex2
                sa1ex2 = sa1ex2+a1ex2
                sa2ex2 = sa2ex2+a2ex2
                sa3ex2 = sa3ex2+a3ex2
              end if
c..................
 700        continue
c
c  precomputes the most used expressions
c
            sa1ex1mt = -sa1ex1*Two
            sa2ex1f = Four*sa2ex1
            sa3ex1me = -Eight*sa3ex1
c
c  The order of the first derivatives is x, y, z
c  The order of the second derivatives is xx, xy, yy, xz, yz, zz
c  The order of the third derivatives is xxx, xxy, xyy, yyy, xxz,
c                                        xyz, yyz, xzz, yzz, zzz
c
            IF(ityp.EQ.1) THEN
c
c -- s function
c
c
c   component  [1,[0,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1
              VAOX(1,ifun,IP) = sa1ex1mt*x
              VAOX(2,ifun,IP) = sa1ex1mt*y
              VAOX(3,ifun,IP) = sa1ex1mt*z
              VAOXX(1,ifun,IP) = sa2ex1f*x2+sa1ex1mt
              VAOXX(2,ifun,IP) = sa2ex1f*x*y
              VAOXX(3,ifun,IP) = sa2ex1f*Y2+sa1ex1mt
              VAOXX(4,ifun,IP) = sa2ex1f*x*z
              VAOXX(5,ifun,IP) = sa2ex1f*y*z
              VAOXX(6,ifun,IP) = sa2ex1f*z2+sa1ex1mt
              VAOXXX(1,ifun,IP) = sa3ex1me*x3+sa2ex1f*Three*x
              VAOXXX(2,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)
              VAOXXX(4,ifun,IP) = sa3ex1me*Y3+sa2ex1f*Three*y
              VAOXXX(5,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*z
              VAOXXX(6,ifun,IP) = sa3ex1me*x*y*z
              VAOXXX(7,ifun,IP) = (sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(8,ifun,IP) = x*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(9,ifun,IP) = y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(10,ifun,IP) = sa3ex1me*z3+sa2ex1f*Three*z
c
            ELSE IF(ityp.EQ.2) THEN
c
c -- p function
c
              x4 = x*x3
              y4 = y*y3
              z4 = z*z3
c
c   component  [1,[1,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x
              VAOX(1,ifun,IP) = sa1ex1mt*x2+sa0ex1
              VAOX(2,ifun,IP) = sa1ex1mt*x*y
              VAOX(3,ifun,IP) = sa1ex1mt*x*z
              VAOXX(1,ifun,IP) = sa2ex1f*x3+sa1ex1mt*Three*x
              VAOXX(2,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)
              VAOXX(4,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*z
              VAOXX(5,ifun,IP) = sa2ex1f*x*y*z
              VAOXX(6,ifun,IP) = x*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*
     1           Three
              VAOXXX(2,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y
              VAOXXX(3,ifun,IP) = sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+
     1           sa1ex1mt
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*z
              VAOXXX(6,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y*z
              VAOXXX(7,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(8,ifun,IP) = sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+
     1           sa1ex1mt
              VAOXXX(9,ifun,IP) = x*y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(10,ifun,IP) = x*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[0,1,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*y
              VAOX(1,ifun,IP) = sa1ex1mt*x*y
              VAOX(2,ifun,IP) = sa1ex1mt*Y2+sa0ex1
              VAOX(3,ifun,IP) = sa1ex1mt*y*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)
              VAOXX(3,ifun,IP) = sa2ex1f*Y3+sa1ex1mt*Three*y
              VAOXX(4,ifun,IP) = sa2ex1f*x*y*z
              VAOXX(5,ifun,IP) = (sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(6,ifun,IP) = y*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y
              VAOXXX(2,ifun,IP) = sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+
     1            sa1ex1mt
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)
              VAOXXX(4,ifun,IP) = sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt*
     1           Three
              VAOXXX(5,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y*z
              VAOXXX(6,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(7,ifun,IP) = (sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(8,ifun,IP) = x*y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(9,ifun,IP) = sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+
     1           sa1ex1mt
              VAOXXX(10,ifun,IP) = y*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[0,0,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*z
              VAOX(1,ifun,IP) = sa1ex1mt*x*z
              VAOX(2,ifun,IP) = sa1ex1mt*y*z
              VAOX(3,ifun,IP) = sa1ex1mt*z2+sa0ex1
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*z
              VAOXX(2,ifun,IP) = sa2ex1f*x*y*z
              VAOXX(3,ifun,IP) = (sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(4,ifun,IP) = x*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(5,ifun,IP) = y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(6,ifun,IP) = sa2ex1f*z3+sa1ex1mt*Three*z
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*z
              VAOXXX(2,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y*z
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(4,ifun,IP) = (sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(5,ifun,IP) = sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+
     1           sa1ex1mt
              VAOXXX(6,ifun,IP) = x*y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(7,ifun,IP) = sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+
     1           sa1ex1mt
              VAOXXX(8,ifun,IP) = x*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(9,ifun,IP) = y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(10,ifun,IP) = sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1mt*
     1           Three
c
            ELSE IF(ityp.EQ.3) THEN
c
c -- l function (sp)
c
c
c   component  [1,[0,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1
              VAOX(1,ifun,IP) = sa1ex1mt*x
              VAOX(2,ifun,IP) = sa1ex1mt*y
              VAOX(3,ifun,IP) = sa1ex1mt*z
              VAOXX(1,ifun,IP) = sa2ex1f*x2+sa1ex1mt
              VAOXX(2,ifun,IP) = sa2ex1f*x*y
              VAOXX(3,ifun,IP) = sa2ex1f*Y2+sa1ex1mt
              VAOXX(4,ifun,IP) = sa2ex1f*x*z
              VAOXX(5,ifun,IP) = sa2ex1f*y*z
              VAOXX(6,ifun,IP) = sa2ex1f*z2+sa1ex1mt
              VAOXXX(1,ifun,IP) = sa3ex1me*x3+sa2ex1f*Three*x
              VAOXXX(2,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)
              VAOXXX(4,ifun,IP) = sa3ex1me*Y3+sa2ex1f*Three*y
              VAOXXX(5,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*z
              VAOXXX(6,ifun,IP) = sa3ex1me*x*y*z
              VAOXXX(7,ifun,IP) = (sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(8,ifun,IP) = x*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(9,ifun,IP) = y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(10,ifun,IP) = sa3ex1me*z3+sa2ex1f*Three*z
c
c
c  precomputes the most used expressions
c
              sa1ex2mt = -sa1ex2*Two
              sa2ex2f = Four*sa2ex2
              sa3ex2me = -Eight*sa3ex2
              x4 = x*x3
              y4 = y*y3
              z4 = z*z3
c
c   component  [1,[1,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex2*x
              VAOX(1,ifun,IP) = sa1ex2mt*x2+sa0ex2
              VAOX(2,ifun,IP) = sa1ex2mt*x*y
              VAOX(3,ifun,IP) = sa1ex2mt*x*z
              VAOXX(1,ifun,IP) = sa2ex2f*x3+sa1ex2mt*Three*x
              VAOXX(2,ifun,IP) = (sa2ex2f*x2+sa1ex2mt)*y
              VAOXX(3,ifun,IP) = x*(sa2ex2f*Y2+sa1ex2mt)
              VAOXX(4,ifun,IP) = (sa2ex2f*x2+sa1ex2mt)*z
              VAOXX(5,ifun,IP) = sa2ex2f*x*y*z
              VAOXX(6,ifun,IP) = x*(sa2ex2f*z2+sa1ex2mt)
              VAOXXX(1,ifun,IP) = sa3ex2me*x4+sa2ex2f*Six*x2+sa1ex2mt*
     1           Three
              VAOXXX(2,ifun,IP) = (sa3ex2me*x3+sa2ex2f*Three*x)*y
              VAOXXX(3,ifun,IP) = sa2ex2f*(Y2+x2)+sa3ex2me*x2*Y2+
     1           sa1ex2mt
              VAOXXX(4,ifun,IP) = x*(sa3ex2me*Y3+sa2ex2f*Three*y)
              VAOXXX(5,ifun,IP) = (sa3ex2me*x3+sa2ex2f*Three*x)*z
              VAOXXX(6,ifun,IP) = (sa3ex2me*x2+sa2ex2f)*y*z
              VAOXXX(7,ifun,IP) = x*(sa3ex2me*Y2+sa2ex2f)*z
              VAOXXX(8,ifun,IP) = sa2ex2f*(z2+x2)+sa3ex2me*x2*z2+
     1           sa1ex2mt
              VAOXXX(9,ifun,IP) = x*y*(sa3ex2me*z2+sa2ex2f)
              VAOXXX(10,ifun,IP) = x*(sa3ex2me*z3+sa2ex2f*Three*z)
c
c   component  [1,[0,1,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex2*y
              VAOX(1,ifun,IP) = sa1ex2mt*x*y
              VAOX(2,ifun,IP) = sa1ex2mt*Y2+sa0ex2
              VAOX(3,ifun,IP) = sa1ex2mt*y*z
              VAOXX(1,ifun,IP) = (sa2ex2f*x2+sa1ex2mt)*y
              VAOXX(2,ifun,IP) = x*(sa2ex2f*Y2+sa1ex2mt)
              VAOXX(3,ifun,IP) = sa2ex2f*Y3+sa1ex2mt*Three*y
              VAOXX(4,ifun,IP) = sa2ex2f*x*y*z
              VAOXX(5,ifun,IP) = (sa2ex2f*Y2+sa1ex2mt)*z
              VAOXX(6,ifun,IP) = y*(sa2ex2f*z2+sa1ex2mt)
              VAOXXX(1,ifun,IP) = (sa3ex2me*x3+sa2ex2f*Three*x)*y
              VAOXXX(2,ifun,IP) = sa2ex2f*(Y2+x2)+sa3ex2me*x2*Y2+
     1           sa1ex2mt
              VAOXXX(3,ifun,IP) = x*(sa3ex2me*Y3+sa2ex2f*Three*y)
              VAOXXX(4,ifun,IP) = sa3ex2me*y4+sa2ex2f*Six*Y2+sa1ex2mt*
     1           Three
              VAOXXX(5,ifun,IP) = (sa3ex2me*x2+sa2ex2f)*y*z
              VAOXXX(6,ifun,IP) = x*(sa3ex2me*Y2+sa2ex2f)*z
              VAOXXX(7,ifun,IP) = (sa3ex2me*Y3+sa2ex2f*Three*y)*z
              VAOXXX(8,ifun,IP) = x*y*(sa3ex2me*z2+sa2ex2f)
              VAOXXX(9,ifun,IP) = sa2ex2f*(z2+Y2)+sa3ex2me*Y2*z2+
     1           sa1ex2mt
              VAOXXX(10,ifun,IP) = y*(sa3ex2me*z3+sa2ex2f*Three*z)
c
c   component  [1,[0,0,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex2*z
              VAOX(1,ifun,IP) = sa1ex2mt*x*z
              VAOX(2,ifun,IP) = sa1ex2mt*y*z
              VAOX(3,ifun,IP) = sa1ex2mt*z2+sa0ex2
              VAOXX(1,ifun,IP) = (sa2ex2f*x2+sa1ex2mt)*z
              VAOXX(2,ifun,IP) = sa2ex2f*x*y*z
              VAOXX(3,ifun,IP) = (sa2ex2f*Y2+sa1ex2mt)*z
              VAOXX(4,ifun,IP) = x*(sa2ex2f*z2+sa1ex2mt)
              VAOXX(5,ifun,IP) = y*(sa2ex2f*z2+sa1ex2mt)
              VAOXX(6,ifun,IP) = sa2ex2f*z3+sa1ex2mt*Three*z
              VAOXXX(1,ifun,IP) = (sa3ex2me*x3+sa2ex2f*Three*x)*z
              VAOXXX(2,ifun,IP) = (sa3ex2me*x2+sa2ex2f)*y*z
              VAOXXX(3,ifun,IP) = x*(sa3ex2me*Y2+sa2ex2f)*z
              VAOXXX(4,ifun,IP) = (sa3ex2me*Y3+sa2ex2f*Three*y)*z
              VAOXXX(5,ifun,IP) = sa2ex2f*(z2+x2)+sa3ex2me*x2*z2+
     1           sa1ex2mt
              VAOXXX(6,ifun,IP) = x*y*(sa3ex2me*z2+sa2ex2f)
              VAOXXX(7,ifun,IP) = sa2ex2f*(z2+Y2)+sa3ex2me*Y2*z2+
     1           sa1ex2mt
              VAOXXX(8,ifun,IP) = x*(sa3ex2me*z3+sa2ex2f*Three*z)
              VAOXXX(9,ifun,IP) = y*(sa3ex2me*z3+sa2ex2f*Three*z)
              VAOXXX(10,ifun,IP) = sa3ex2me*z4+sa2ex2f*Six*z2+sa1ex2mt*
     1           Three
c
            ELSE IF(ityp.EQ.4) THEN
c
c -- d function (5 components)    z**2,x**2-y**2,xy,xz,yz
c
              x4 = x*x3
              x5 = x*x4
              y4 = y*y3
              y5 = y*y4
              z4 = z*z3
              z5 = z*z4
c
c   component  [sq3i,[0,0,2],-sq3i/2,[2,0,0],-sq3i/2,[0,2,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = Half*sa0ex1*sq3i*(Two*z2-Y2-x2)
              VAOX(1,ifun,IP) = sq3i*(sa1ex1mt*x*z2+sa1ex1*(x*Y2+x3)-sa0
     1           ex1*x)
              VAOX(2,ifun,IP) = sq3i*(sa1ex1mt*y*z2+sa1ex1*(Y3+x2*y)-sa0
     1           ex1*y)
              VAOX(3,ifun,IP) = sq3i*(sa1ex1mt*z3+sa1ex1*(Y2+x2)*z+sa0ex
     1           1*Two*z)
              VAOXX(1,ifun,IP) = sq3i*(sa1ex1mt*(z2-Two*x2)+sa2ex1f*x2*z
     1           2-sa2ex1*Two*(x2*Y2+x4)+sa1ex1*(Y2+x2)-sa0ex1)
              VAOXX(2,ifun,IP) = sq3i*(sa2ex1f*x*y*z2-sa2ex1*Two*(x*Y3+x
     1           3*y)-sa1ex1mt*Two*x*y)
              VAOXX(3,ifun,IP) = sq3i*(sa1ex1mt*(z2-Two*Y2)+sa2ex1f*Y2*z
     1           2-sa2ex1*Two*(y4+x2*Y2)+sa1ex1*(Y2+x2)-sa0ex1)
              VAOXX(4,ifun,IP) = sq3i*(sa2ex1f*x*z3-sa2ex1*Two*(x*Y2+x3)
     1           *z+sa1ex1mt*x*z)
              VAOXX(5,ifun,IP) = sq3i*(sa2ex1f*y*z3-sa2ex1*Two*(Y3+x2*y)
     1           *z+sa1ex1mt*y*z)
              VAOXX(6,ifun,IP) = sq3i*(sa2ex1f*z4-sa2ex1*Two*(Y2+x2)*z2+
     1           Five*sa1ex1mt*z2+sa1ex1*(Y2+x2)+sa0ex1*Two)
              VAOXXX(1,ifun,IP) = sq3i*(sa2ex1f*(Three*x*z2-x*Y2-Four*x3
     1           )+sa3ex1me*x3*z2+Four*sa3ex1*(x3*Y2+x5)-sa2ex1*Two*(x*Y
     2           2+x3)-sa1ex1mt*Six*x)
              VAOXXX(2,ifun,IP) = sq3i*(sa2ex1f*y*(z2-Three*x2)+sa3ex1me
     1           *x2*y*z2+Four*sa3ex1*(x2*Y3+x4*y)-sa2ex1*Two*(Y3+x2*y)-
     2           sa1ex1mt*Two*y)
              VAOXXX(3,ifun,IP) = sq3i*(sa2ex1f*x*(z2-Three*Y2)+sa3ex1me
     1           *x*Y2*z2+Four*sa3ex1*(x*y4+x3*Y2)-sa2ex1*Two*(x*Y2+x3)-
     2           sa1ex1mt*Two*x)
              VAOXXX(4,ifun,IP) = sq3i*(sa2ex1f*(Three*y*z2-Four*Y3-x2*y
     1           )+sa3ex1me*Y3*z2+Four*sa3ex1*(y5+x2*Y3)-sa2ex1*Two*(Y3+
     2           x2*y)-sa1ex1mt*Six*y)
              VAOXXX(5,ifun,IP) = sq3i*(sa3ex1me*x2*z3+sa2ex1f*z3+Four*s
     1           a3ex1*(x2*Y2+x4)*z-sa2ex1*Two*(Y2+x2)*z+sa1ex1mt*z)
              VAOXXX(6,ifun,IP) = sq3i*(sa3ex1me*x*y*z3+Four*sa3ex1*(x*Y
     1           3+x3*y)*z)
              VAOXXX(7,ifun,IP) = sq3i*(sa3ex1me*Y2*z3+sa2ex1f*z3+Four*s
     1           a3ex1*(y4+x2*Y2)*z-sa2ex1*Two*(Y2+x2)*z+sa1ex1mt*z)
              VAOXXX(8,ifun,IP) = sq3i*(sa3ex1me*x*z4+Four*sa3ex1*(x*Y2+
     1           x3)*z2+Four*sa2ex1f*x*z2-sa2ex1*Two*(x*Y2+x3)+sa1ex1mt*
     2           x)
              VAOXXX(9,ifun,IP) = sq3i*(sa3ex1me*y*z4+Four*sa3ex1*(Y3+x2
     1           *y)*z2+Four*sa2ex1f*y*z2-sa2ex1*Two*(Y3+x2*y)+sa1ex1mt*
     2           y)
              VAOXXX(10,ifun,IP) = sq3i*(sa3ex1me*z5+sa2ex1f*(Nine*z3-Y2
     1           *z-x2*z)+Four*sa3ex1*(Y2+x2)*z3-sa2ex1*Two*(Y2+x2)*z+sa
     2           1ex1mt*Twelve*z)
c
c   component  [1/2,[2,0,0],-1/2,[0,2,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = -Half*sa0ex1*(Y2-x2)
              VAOX(1,ifun,IP) = sa1ex1*(x*Y2-x3)+sa0ex1*x
              VAOX(2,ifun,IP) = sa1ex1*(Y3-x2*y)-sa0ex1*y
              VAOX(3,ifun,IP) = sa1ex1*(Y2-x2)*z
              VAOXX(1,ifun,IP) = -sa2ex1*Two*(x2*Y2-x4)+sa1ex1*(Y2-x2)+s
     1           a1ex1mt*Two*x2+sa0ex1
              VAOXX(2,ifun,IP) = -sa2ex1*Two*(x*Y3-x3*y)
              VAOXX(3,ifun,IP) = -sa2ex1*Two*(y4-x2*Y2)+sa1ex1*(Y2-x2)-s
     1           a1ex1mt*Two*Y2-sa0ex1
              VAOXX(4,ifun,IP) = -(sa2ex1*Two*(x*Y2-x3)-sa1ex1mt*x)*z
              VAOXX(5,ifun,IP) = -(sa2ex1*Two*(Y3-x2*y)+sa1ex1mt*y)*z
              VAOXX(6,ifun,IP) = -(Y2-x2)*(sa2ex1*Two*z2-sa1ex1)
              VAOXXX(1,ifun,IP) = Four*sa3ex1*(x3*Y2-x5)-sa2ex1f*(x*Y2-F
     1           our*x3)-sa2ex1*Two*(x*Y2-x3)+sa1ex1mt*Six*x
              VAOXXX(2,ifun,IP) = Four*sa3ex1*(x2*Y3-x4*y)-sa2ex1*Two*(Y
     1           3-x2*y)+sa2ex1f*x2*y
              VAOXXX(3,ifun,IP) = Four*sa3ex1*(x*y4-x3*Y2)-sa2ex1*Two*(x
     1           *Y2-x3)-sa2ex1f*x*Y2
              VAOXXX(4,ifun,IP) = Four*sa3ex1*(y5-x2*Y3)-sa2ex1f*(Four*Y
     1           3-x2*y)-sa2ex1*Two*(Y3-x2*y)-sa1ex1mt*Six*y
              VAOXXX(5,ifun,IP) = (Four*sa3ex1*(x2*Y2-x4)-sa2ex1*Two*(Y2
     1           -x2)+sa2ex1f*Two*x2+sa1ex1mt)*z
              VAOXXX(6,ifun,IP) = Four*sa3ex1*(x*Y3-x3*y)*z
              VAOXXX(7,ifun,IP) = (Four*sa3ex1*(y4-x2*Y2)-sa2ex1*Two*(Y2
     1           -x2)-sa2ex1f*Two*Y2-sa1ex1mt)*z
              VAOXXX(8,ifun,IP) = Four*sa3ex1*(x*Y2-x3)*z2+sa2ex1f*x*z2-
     1           sa2ex1*Two*(x*Y2-x3)+sa1ex1mt*x
              VAOXXX(9,ifun,IP) = Four*sa3ex1*(Y3-x2*y)*z2-sa2ex1f*y*z2-
     1           sa2ex1*Two*(Y3-x2*y)-sa1ex1mt*y
              VAOXXX(10,ifun,IP) = (Y2-x2)*(Four*sa3ex1*z3-sa2ex1*Two*z-
     1           sa2ex1f*z)
c
c   component  [1,[1,1,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*y
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*y
              VAOX(2,ifun,IP) = x*(sa1ex1mt*Y2+sa0ex1)
              VAOX(3,ifun,IP) = sa1ex1mt*x*y*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*y
              VAOXX(2,ifun,IP) = sa1ex1mt*(Y2+x2)+sa2ex1f*x2*Y2+sa0ex1
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Three*y)
              VAOXX(4,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(5,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(6,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*y
              VAOXXX(2,ifun,IP) = sa2ex1f*(Three*x*Y2+x3)+sa3ex1me*x3*Y2
     1           +sa1ex1mt*Three*x
              VAOXXX(3,ifun,IP) = sa2ex1f*(Y3+Three*x2*y)+sa3ex1me*x2*Y3
     1           +sa1ex1mt*Three*y
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt
     1           *Three)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(6,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(7,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(8,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(9,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(10,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[1,0,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*z
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*z
              VAOX(2,ifun,IP) = sa1ex1mt*x*y*z
              VAOX(3,ifun,IP) = x*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*z
              VAOXX(2,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(4,ifun,IP) = sa1ex1mt*(z2+x2)+sa2ex1f*x2*z2+sa0ex1
              VAOXX(5,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(6,ifun,IP) = x*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*z
              VAOXXX(2,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(3,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(5,ifun,IP) = sa2ex1f*(Three*x*z2+x3)+sa3ex1me*x3*z2
     1           +sa1ex1mt*Three*x
              VAOXXX(6,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(7,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(8,ifun,IP) = sa2ex1f*(z3+Three*x2*z)+sa3ex1me*x2*z3
     1           +sa1ex1mt*Three*z
              VAOXXX(9,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(10,ifun,IP) = x*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1m
     1           t*Three)
c
c   component  [1,[0,1,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*y*z
              VAOX(1,ifun,IP) = sa1ex1mt*x*y*z
              VAOX(2,ifun,IP) = (sa1ex1mt*Y2+sa0ex1)*z
              VAOX(3,ifun,IP) = y*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(3,ifun,IP) = (sa2ex1f*Y3+sa1ex1mt*Three*y)*z
              VAOXX(4,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(5,ifun,IP) = sa1ex1mt*(z2+Y2)+sa2ex1f*Y2*z2+sa0ex1
              VAOXX(6,ifun,IP) = y*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(2,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(4,ifun,IP) = (sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt*T
     1           hree)*z
              VAOXXX(5,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(6,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(7,ifun,IP) = sa2ex1f*(Three*y*z2+Y3)+sa3ex1me*Y3*z2
     1           +sa1ex1mt*Three*y
              VAOXXX(8,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(9,ifun,IP) = sa2ex1f*(z3+Three*Y2*z)+sa3ex1me*Y2*z3
     1           +sa1ex1mt*Three*z
              VAOXXX(10,ifun,IP) = y*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1m
     1           t*Three)
c
            ELSE IF(ityp.EQ.5) THEN
c
c -- d function (6 components)    xx,yy,zz,xy,xz,yz
c
              x4 = x*x3
              x5 = x*x4
              y4 = y*y3
              y5 = y*y4
              z4 = z*z3
              z5 = z*z4
c
c   component  [1,[2,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x2
              VAOX(1,ifun,IP) = sa1ex1mt*x3+sa0ex1*Two*x
              VAOX(2,ifun,IP) = sa1ex1mt*x2*y
              VAOX(3,ifun,IP) = sa1ex1mt*x2*z
              VAOXX(1,ifun,IP) = sa2ex1f*x4+Five*sa1ex1mt*x2+sa0ex1*Two
              VAOXX(2,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Two*x)*y
              VAOXX(3,ifun,IP) = x2*(sa2ex1f*Y2+sa1ex1mt)
              VAOXX(4,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Two*x)*z
              VAOXX(5,ifun,IP) = sa2ex1f*x2*y*z
              VAOXX(6,ifun,IP) = x2*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = sa3ex1me*x5+Nine*sa2ex1f*x3+sa1ex1mt*T
     1           welve*x
              VAOXXX(2,ifun,IP) = (sa3ex1me*x4+Five*sa2ex1f*x2+sa1ex1mt*
     1           Two)*y
              VAOXXX(3,ifun,IP) = sa2ex1f*(Two*x*Y2+x3)+sa3ex1me*x3*Y2+s
     1           a1ex1mt*Two*x
              VAOXXX(4,ifun,IP) = x2*(sa3ex1me*Y3+sa2ex1f*Three*y)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x4+Five*sa2ex1f*x2+sa1ex1mt*
     1           Two)*z
              VAOXXX(6,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Two*x)*y*z
              VAOXXX(7,ifun,IP) = x2*(sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(8,ifun,IP) = sa2ex1f*(Two*x*z2+x3)+sa3ex1me*x3*z2+s
     1           a1ex1mt*Two*x
              VAOXXX(9,ifun,IP) = x2*y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(10,ifun,IP) = x2*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[0,2,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*Y2
              VAOX(1,ifun,IP) = sa1ex1mt*x*Y2
              VAOX(2,ifun,IP) = sa1ex1mt*Y3+sa0ex1*Two*y
              VAOX(3,ifun,IP) = sa1ex1mt*Y2*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*Y2
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Two*y)
              VAOXX(3,ifun,IP) = sa2ex1f*y4+Five*sa1ex1mt*Y2+sa0ex1*Two
              VAOXX(4,ifun,IP) = sa2ex1f*x*Y2*z
              VAOXX(5,ifun,IP) = (sa2ex1f*Y3+sa1ex1mt*Two*y)*z
              VAOXX(6,ifun,IP) = Y2*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*Y2
              VAOXXX(2,ifun,IP) = sa2ex1f*(Y3+Two*x2*y)+sa3ex1me*x2*Y3+s
     1           a1ex1mt*Two*y
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*y4+Five*sa2ex1f*Y2+sa1ex1m
     1           t*Two)
              VAOXXX(4,ifun,IP) = sa3ex1me*y5+Nine*sa2ex1f*Y3+sa1ex1mt*T
     1           welve*y
              VAOXXX(5,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*Y2*z
              VAOXXX(6,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Two*y)*z
              VAOXXX(7,ifun,IP) = (sa3ex1me*y4+Five*sa2ex1f*Y2+sa1ex1mt*
     1           Two)*z
              VAOXXX(8,ifun,IP) = x*Y2*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(9,ifun,IP) = sa2ex1f*(Two*y*z2+Y3)+sa3ex1me*Y3*z2+s
     1           a1ex1mt*Two*y
              VAOXXX(10,ifun,IP) = Y2*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[0,0,2]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*z2
              VAOX(1,ifun,IP) = sa1ex1mt*x*z2
              VAOX(2,ifun,IP) = sa1ex1mt*y*z2
              VAOX(3,ifun,IP) = sa1ex1mt*z3+sa0ex1*Two*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*z2
              VAOXX(2,ifun,IP) = sa2ex1f*x*y*z2
              VAOXX(3,ifun,IP) = (sa2ex1f*Y2+sa1ex1mt)*z2
              VAOXX(4,ifun,IP) = x*(sa2ex1f*z3+sa1ex1mt*Two*z)
              VAOXX(5,ifun,IP) = y*(sa2ex1f*z3+sa1ex1mt*Two*z)
              VAOXX(6,ifun,IP) = sa2ex1f*z4+Five*sa1ex1mt*z2+sa0ex1*Two
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*z2
              VAOXXX(2,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y*z2
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)*z2
              VAOXXX(4,ifun,IP) = (sa3ex1me*Y3+sa2ex1f*Three*y)*z2
              VAOXXX(5,ifun,IP) = sa2ex1f*(z3+Two*x2*z)+sa3ex1me*x2*z3+s
     1           a1ex1mt*Two*z
              VAOXXX(6,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Two*z)
              VAOXXX(7,ifun,IP) = sa2ex1f*(z3+Two*Y2*z)+sa3ex1me*Y2*z3+s
     1           a1ex1mt*Two*z
              VAOXXX(8,ifun,IP) = x*(sa3ex1me*z4+Five*sa2ex1f*z2+sa1ex1m
     1           t*Two)
              VAOXXX(9,ifun,IP) = y*(sa3ex1me*z4+Five*sa2ex1f*z2+sa1ex1m
     1           t*Two)
              VAOXXX(10,ifun,IP) = sa3ex1me*z5+Nine*sa2ex1f*z3+sa1ex1mt*
     1           Twelve*z
c
c   component  [1,[1,1,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*y
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*y
              VAOX(2,ifun,IP) = x*(sa1ex1mt*Y2+sa0ex1)
              VAOX(3,ifun,IP) = sa1ex1mt*x*y*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*y
              VAOXX(2,ifun,IP) = sa1ex1mt*(Y2+x2)+sa2ex1f*x2*Y2+sa0ex1
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Three*y)
              VAOXX(4,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(5,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(6,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*y
              VAOXXX(2,ifun,IP) = sa2ex1f*(Three*x*Y2+x3)+sa3ex1me*x3*Y2
     1           +sa1ex1mt*Three*x
              VAOXXX(3,ifun,IP) = sa2ex1f*(Y3+Three*x2*y)+sa3ex1me*x2*Y3
     1           +sa1ex1mt*Three*y
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt
     1           *Three)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(6,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(7,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(8,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(9,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(10,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[1,0,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*z
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*z
              VAOX(2,ifun,IP) = sa1ex1mt*x*y*z
              VAOX(3,ifun,IP) = x*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*z
              VAOXX(2,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(4,ifun,IP) = sa1ex1mt*(z2+x2)+sa2ex1f*x2*z2+sa0ex1
              VAOXX(5,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(6,ifun,IP) = x*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*z
              VAOXXX(2,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(3,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(5,ifun,IP) = sa2ex1f*(Three*x*z2+x3)+sa3ex1me*x3*z2
     1           +sa1ex1mt*Three*x
              VAOXXX(6,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(7,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(8,ifun,IP) = sa2ex1f*(z3+Three*x2*z)+sa3ex1me*x2*z3
     1           +sa1ex1mt*Three*z
              VAOXXX(9,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(10,ifun,IP) = x*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1m
     1           t*Three)
c
c   component  [1,[0,1,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*y*z
              VAOX(1,ifun,IP) = sa1ex1mt*x*y*z
              VAOX(2,ifun,IP) = (sa1ex1mt*Y2+sa0ex1)*z
              VAOX(3,ifun,IP) = y*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(3,ifun,IP) = (sa2ex1f*Y3+sa1ex1mt*Three*y)*z
              VAOXX(4,ifun,IP) = x*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(5,ifun,IP) = sa1ex1mt*(z2+Y2)+sa2ex1f*Y2*z2+sa0ex1
              VAOXX(6,ifun,IP) = y*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z
              VAOXXX(2,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(4,ifun,IP) = (sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt*T
     1           hree)*z
              VAOXXX(5,ifun,IP) = y*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1e
     1           x1mt)
              VAOXXX(6,ifun,IP) = x*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1e
     1           x1mt)
              VAOXXX(7,ifun,IP) = sa2ex1f*(Three*y*z2+Y3)+sa3ex1me*Y3*z2
     1           +sa1ex1mt*Three*y
              VAOXXX(8,ifun,IP) = x*y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(9,ifun,IP) = sa2ex1f*(z3+Three*Y2*z)+sa3ex1me*Y2*z3
     1           +sa1ex1mt*Three*z
              VAOXXX(10,ifun,IP) = y*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1m
     1           t*Three)
c
            ELSE IF(ityp.EQ.6) THEN
c
c -- f function (7 components)    (5xxy-rry),(5xxz-rrz),(5yyx-rrx),
c                                 (5yyz-rrz),(5zzx-rrx),(5zzy-rry),
c
              x4 = x*x3
              x5 = x*x4
              x6 = x*x5
              y4 = y*y3
              y5 = y*y4
              y6 = y*y5
              z4 = z*z3
              z5 = z*z4
              z6 = z*z5
c
c   component  [4,[2,1,0],-1,[0,3,0],-1,[0,1,2]]
c
              ifun = ifun+1
              VAO(ifun,IP) = -sa0ex1*(y*z2+Y3-Four*x2*y)
              VAOX(1,ifun,IP) = Eight*sa0ex1*x*y-sa1ex1mt*(x*y*z2+x*Y3-F
     1           our*x3*y)
              VAOX(2,ifun,IP) = -sa1ex1mt*(Y2*z2+y4-Four*x2*Y2)-sa0ex1*(
     1           z2+Three*Y2-Four*x2)
              VAOX(3,ifun,IP) = -sa1ex1mt*(y*z3+Y3*z-Four*x2*y*z)-sa0ex1
     1           *Two*y*z
              VAOXX(1,ifun,IP) = -sa2ex1f*(x2*y*z2+x2*Y3-Four*x4*y)-sa1e
     1           x1mt*(y*z2+Y3-Twenty*x2*y)+Eight*sa0ex1*y
              VAOXX(2,ifun,IP) = -sa2ex1f*(x*Y2*z2+x*y4-Four*x3*Y2)-sa1e
     1           x1mt*(x*z2-Five*x*Y2-Four*x3)+Eight*sa0ex1*x
              VAOXX(3,ifun,IP) = -sa2ex1f*(Y3*z2+y5-Four*x2*Y3)-sa1ex1mt
     1           *(Three*y*z2+Seven*Y3-Twelve*x2*y)-sa0ex1*Six*y
              VAOXX(4,ifun,IP) = sa1ex1mt*Six*x*y*z-sa2ex1f*(x*y*z3+x*Y3
     1           *z-Four*x3*y*z)
              VAOXX(5,ifun,IP) = -sa2ex1f*(Y2*z3+y4*z-Four*x2*Y2*z)-sa1e
     1           x1mt*(z3+Five*Y2*z-Four*x2*z)-sa0ex1*Two*z
              VAOXX(6,ifun,IP) = -sa2ex1f*(y*z4+Y3*z2-Four*x2*y*z2)-sa1e
     1           x1mt*(Five*y*z2+Y3-Four*x2*y)-sa0ex1*Two*y
              VAOXXX(1,ifun,IP) = -sa3ex1me*(x3*y*z2+x3*Y3-Four*x5*y)-sa
     1           2ex1f*(Three*x*y*z2+Three*x*Y3-ThSix*x3*y)+FEight*sa1ex
     2           1mt*x*y
              VAOXXX(2,ifun,IP) = -sa3ex1me*(x2*Y2*z2+x2*y4-Four*x4*Y2)-
     1           sa2ex1f*(Y2*z2+x2*z2+y4-TSeven*x2*Y2-Four*x4)-sa1ex1mt*
     2           (z2-Five*Y2-Twenty*x2)+Eight*sa0ex1
              VAOXXX(3,ifun,IP) = -sa3ex1me*(x*Y3*z2+x*y5-Four*x3*Y3)-sa
     1           2ex1f*(Three*x*y*z2-x*Y3-Twelve*x3*y)+sa1ex1mt*TEight*x
     2           *y
              VAOXXX(4,ifun,IP) = -sa3ex1me*(y4*z2+y6-Four*x2*y4)-sa2ex1
     1           f*(Six*Y2*z2+Twelve*y4-TwFour*x2*Y2)-sa1ex1mt*(Three*z2
     2           +TwSeven*Y2-Twelve*x2)-sa0ex1*Six
              VAOXXX(5,ifun,IP) = -sa3ex1me*(x2*y*z3+x2*Y3*z-Four*x4*y*z
     1           )-sa2ex1f*(y*z3+Y3*z-TEight*x2*y*z)+sa1ex1mt*Six*y*z
              VAOXXX(6,ifun,IP) = -sa3ex1me*(x*Y2*z3+x*y4*z-Four*x3*Y2*z
     1           )-sa2ex1f*(x*z3-Three*x*Y2*z-Four*x3*z)+sa1ex1mt*Six*x*
     2           z
              VAOXXX(7,ifun,IP) = -sa3ex1me*(Y3*z3+y5*z-Four*x2*Y3*z)-sa
     1           2ex1f*(Three*y*z3+Nine*Y3*z-Twelve*x2*y*z)-sa1ex1mt*Twe
     2           lve*y*z
              VAOXXX(8,ifun,IP) = -sa3ex1me*(x*y*z4+x*Y3*z2-Four*x3*y*z2
     1           )+sa2ex1f*(Three*x*y*z2-x*Y3+Four*x3*y)+sa1ex1mt*Six*x*
     2           y
              VAOXXX(9,ifun,IP) = -sa3ex1me*(Y2*z4+y4*z2-Four*x2*Y2*z2)-
     1           sa2ex1f*(z4+Eight*Y2*z2-Four*x2*z2+y4-Four*x2*Y2)-sa1ex
     2           1mt*(Five*z2+Five*Y2-Four*x2)-sa0ex1*Two
              VAOXXX(10,ifun,IP) = -sa3ex1me*(y*z5+Y3*z3-Four*x2*y*z3)-s
     1           a2ex1f*(Nine*y*z3+Three*Y3*z-Twelve*x2*y*z)-sa1ex1mt*Tw
     2           elve*y*z
c
c   component  [4,[2,0,1],-1,[0,2,1],-1,[0,0,3]]
c
              ifun = ifun+1
              VAO(ifun,IP) = -sa0ex1*(z3+Y2*z-Four*x2*z)
              VAOX(1,ifun,IP) = Eight*sa0ex1*x*z-sa1ex1mt*(x*z3+x*Y2*z-F
     1           our*x3*z)
              VAOX(2,ifun,IP) = -sa1ex1mt*(y*z3+Y3*z-Four*x2*y*z)-sa0ex1
     1           *Two*y*z
              VAOX(3,ifun,IP) = -sa1ex1mt*(z4+Y2*z2-Four*x2*z2)-sa0ex1*(
     1           Three*z2+Y2-Four*x2)
              VAOXX(1,ifun,IP) = -sa2ex1f*(x2*z3+x2*Y2*z-Four*x4*z)-sa1e
     1           x1mt*(z3+Y2*z-Twenty*x2*z)+Eight*sa0ex1*z
              VAOXX(2,ifun,IP) = sa1ex1mt*Six*x*y*z-sa2ex1f*(x*y*z3+x*Y3
     1           *z-Four*x3*y*z)
              VAOXX(3,ifun,IP) = -sa2ex1f*(Y2*z3+y4*z-Four*x2*Y2*z)-sa1e
     1           x1mt*(z3+Five*Y2*z-Four*x2*z)-sa0ex1*Two*z
              VAOXX(4,ifun,IP) = -sa2ex1f*(x*z4+x*Y2*z2-Four*x3*z2)+sa1e
     1           x1mt*(Five*x*z2-x*Y2+Four*x3)+Eight*sa0ex1*x
              VAOXX(5,ifun,IP) = -sa2ex1f*(y*z4+Y3*z2-Four*x2*y*z2)-sa1e
     1           x1mt*(Five*y*z2+Y3-Four*x2*y)-sa0ex1*Two*y
              VAOXX(6,ifun,IP) = -sa2ex1f*(z5+Y2*z3-Four*x2*z3)-sa1ex1mt
     1           *(Seven*z3+Three*Y2*z-Twelve*x2*z)-sa0ex1*Six*z
              VAOXXX(1,ifun,IP) = -sa3ex1me*(x3*z3+x3*Y2*z-Four*x5*z)-sa
     1           2ex1f*(Three*x*z3+Three*x*Y2*z-ThSix*x3*z)+FEight*sa1ex
     2           1mt*x*z
              VAOXXX(2,ifun,IP) = -sa3ex1me*(x2*y*z3+x2*Y3*z-Four*x4*y*z
     1           )-sa2ex1f*(y*z3+Y3*z-TEight*x2*y*z)+sa1ex1mt*Six*y*z
              VAOXXX(3,ifun,IP) = -sa3ex1me*(x*Y2*z3+x*y4*z-Four*x3*Y2*z
     1           )-sa2ex1f*(x*z3-Three*x*Y2*z-Four*x3*z)+sa1ex1mt*Six*x*
     2           z
              VAOXXX(4,ifun,IP) = -sa3ex1me*(Y3*z3+y5*z-Four*x2*Y3*z)-sa
     1           2ex1f*(Three*y*z3+Nine*Y3*z-Twelve*x2*y*z)-sa1ex1mt*Twe
     2           lve*y*z
              VAOXXX(5,ifun,IP) = -sa3ex1me*(x2*z4+x2*Y2*z2-Four*x4*z2)-
     1           sa2ex1f*(z4+Y2*z2-TSeven*x2*z2+x2*Y2-Four*x4)+sa1ex1mt*
     2           (Five*z2-Y2+Twenty*x2)+Eight*sa0ex1
              VAOXXX(6,ifun,IP) = -sa3ex1me*(x*y*z4+x*Y3*z2-Four*x3*y*z2
     1           )+sa2ex1f*(Three*x*y*z2-x*Y3+Four*x3*y)+sa1ex1mt*Six*x*
     2           y
              VAOXXX(7,ifun,IP) = -sa3ex1me*(Y2*z4+y4*z2-Four*x2*Y2*z2)-
     1           sa2ex1f*(z4+Eight*Y2*z2-Four*x2*z2+y4-Four*x2*Y2)-sa1ex
     2           1mt*(Five*z2+Five*Y2-Four*x2)-sa0ex1*Two
              VAOXXX(8,ifun,IP) = -sa3ex1me*(x*z5+x*Y2*z3-Four*x3*z3)+sa
     1           2ex1f*(x*z3-Three*x*Y2*z+Twelve*x3*z)+sa1ex1mt*TEight*x
     2           *z
              VAOXXX(9,ifun,IP) = -sa3ex1me*(y*z5+Y3*z3-Four*x2*y*z3)-sa
     1           2ex1f*(Nine*y*z3+Three*Y3*z-Twelve*x2*y*z)-sa1ex1mt*Twe
     2           lve*y*z
              VAOXXX(10,ifun,IP) = -sa3ex1me*(z6+Y2*z4-Four*x2*z4)-sa2ex
     1           1f*(Twelve*z4+Six*Y2*z2-TwFour*x2*z2)-sa1ex1mt*(TwSeven
     2           *z2+Three*Y2-Twelve*x2)-sa0ex1*Six
c
c   component  [4,[1,2,0],-1,[3,0,0],-1,[1,0,2]]
c
              ifun = ifun+1
              VAO(ifun,IP) = -sa0ex1*(x*z2-Four*x*Y2+x3)
              VAOX(1,ifun,IP) = -sa1ex1mt*(x2*z2-Four*x2*Y2+x4)-sa0ex1*(
     1           z2-Four*Y2+Three*x2)
              VAOX(2,ifun,IP) = Eight*sa0ex1*x*y-sa1ex1mt*(x*y*z2-Four*x
     1           *Y3+x3*y)
              VAOX(3,ifun,IP) = -sa1ex1mt*(x*z3-Four*x*Y2*z+x3*z)-sa0ex1
     1           *Two*x*z
              VAOXX(1,ifun,IP) = -sa2ex1f*(x3*z2-Four*x3*Y2+x5)-sa1ex1mt
     1           *(Three*x*z2-Twelve*x*Y2+Seven*x3)-sa0ex1*Six*x
              VAOXX(2,ifun,IP) = -sa2ex1f*(x2*y*z2-Four*x2*Y3+x4*y)-sa1e
     1           x1mt*(y*z2-Four*Y3-Five*x2*y)+Eight*sa0ex1*y
              VAOXX(3,ifun,IP) = -sa2ex1f*(x*Y2*z2-Four*x*y4+x3*Y2)-sa1e
     1           x1mt*(x*z2-Twenty*x*Y2+x3)+Eight*sa0ex1*x
              VAOXX(4,ifun,IP) = -sa2ex1f*(x2*z3-Four*x2*Y2*z+x4*z)-sa1e
     1           x1mt*(z3-Four*Y2*z+Five*x2*z)-sa0ex1*Two*z
              VAOXX(5,ifun,IP) = sa1ex1mt*Six*x*y*z-sa2ex1f*(x*y*z3-Four
     1           *x*Y3*z+x3*y*z)
              VAOXX(6,ifun,IP) = -sa2ex1f*(x*z4-Four*x*Y2*z2+x3*z2)-sa1e
     1           x1mt*(Five*x*z2-Four*x*Y2+x3)-sa0ex1*Two*x
              VAOXXX(1,ifun,IP) = -sa3ex1me*(x4*z2-Four*x4*Y2+x6)-sa2ex1
     1           f*(Six*x2*z2-TwFour*x2*Y2+Twelve*x4)-sa1ex1mt*(Three*z2
     2           -Twelve*Y2+TwSeven*x2)-sa0ex1*Six
              VAOXXX(2,ifun,IP) = -sa3ex1me*(x3*y*z2-Four*x3*Y3+x5*y)-sa
     1           2ex1f*(Three*x*y*z2-Twelve*x*Y3-x3*y)+sa1ex1mt*TEight*x
     2           *y
              VAOXXX(3,ifun,IP) = -sa3ex1me*(x2*Y2*z2-Four*x2*y4+x4*Y2)-
     1           sa2ex1f*(Y2*z2+x2*z2-Four*y4-TSeven*x2*Y2+x4)-sa1ex1mt*
     2           (z2-Twenty*Y2-Five*x2)+Eight*sa0ex1
              VAOXXX(4,ifun,IP) = -sa3ex1me*(x*Y3*z2-Four*x*y5+x3*Y3)-sa
     1           2ex1f*(Three*x*y*z2-ThSix*x*Y3+Three*x3*y)+FEight*sa1ex
     2           1mt*x*y
              VAOXXX(5,ifun,IP) = -sa3ex1me*(x3*z3-Four*x3*Y2*z+x5*z)-sa
     1           2ex1f*(Three*x*z3-Twelve*x*Y2*z+Nine*x3*z)-sa1ex1mt*Twe
     2           lve*x*z
              VAOXXX(6,ifun,IP) = -sa3ex1me*(x2*y*z3-Four*x2*Y3*z+x4*y*z
     1           )-sa2ex1f*(y*z3-Four*Y3*z-Three*x2*y*z)+sa1ex1mt*Six*y*
     2           z
              VAOXXX(7,ifun,IP) = -sa3ex1me*(x*Y2*z3-Four*x*y4*z+x3*Y2*z
     1           )-sa2ex1f*(x*z3-TEight*x*Y2*z+x3*z)+sa1ex1mt*Six*x*z
              VAOXXX(8,ifun,IP) = -sa3ex1me*(x2*z4-Four*x2*Y2*z2+x4*z2)-
     1           sa2ex1f*(z4-Four*Y2*z2+Eight*x2*z2-Four*x2*Y2+x4)-sa1ex
     2           1mt*(Five*z2-Four*Y2+Five*x2)-sa0ex1*Two
              VAOXXX(9,ifun,IP) = -sa3ex1me*(x*y*z4-Four*x*Y3*z2+x3*y*z2
     1           )+sa2ex1f*(Three*x*y*z2+Four*x*Y3-x3*y)+sa1ex1mt*Six*x*
     2           y
              VAOXXX(10,ifun,IP) = -sa3ex1me*(x*z5-Four*x*Y2*z3+x3*z3)-s
     1           a2ex1f*(Nine*x*z3-Twelve*x*Y2*z+Three*x3*z)-sa1ex1mt*Tw
     2           elve*x*z
c
c   component  [4,[0,2,1],-1,[2,0,1],-1,[0,0,3]]
c
              ifun = ifun+1
              VAO(ifun,IP) = -sa0ex1*(z3-Four*Y2*z+x2*z)
              VAOX(1,ifun,IP) = -sa1ex1mt*(x*z3-Four*x*Y2*z+x3*z)-sa0ex1
     1           *Two*x*z
              VAOX(2,ifun,IP) = Eight*sa0ex1*y*z-sa1ex1mt*(y*z3-Four*Y3*
     1           z+x2*y*z)
              VAOX(3,ifun,IP) = -sa1ex1mt*(z4-Four*Y2*z2+x2*z2)-sa0ex1*(
     1           Three*z2-Four*Y2+x2)
              VAOXX(1,ifun,IP) = -sa2ex1f*(x2*z3-Four*x2*Y2*z+x4*z)-sa1e
     1           x1mt*(z3-Four*Y2*z+Five*x2*z)-sa0ex1*Two*z
              VAOXX(2,ifun,IP) = sa1ex1mt*Six*x*y*z-sa2ex1f*(x*y*z3-Four
     1           *x*Y3*z+x3*y*z)
              VAOXX(3,ifun,IP) = -sa2ex1f*(Y2*z3-Four*y4*z+x2*Y2*z)-sa1e
     1           x1mt*(z3-Twenty*Y2*z+x2*z)+Eight*sa0ex1*z
              VAOXX(4,ifun,IP) = -sa2ex1f*(x*z4-Four*x*Y2*z2+x3*z2)-sa1e
     1           x1mt*(Five*x*z2-Four*x*Y2+x3)-sa0ex1*Two*x
              VAOXX(5,ifun,IP) = -sa2ex1f*(y*z4-Four*Y3*z2+x2*y*z2)+sa1e
     1           x1mt*(Five*y*z2+Four*Y3-x2*y)+Eight*sa0ex1*y
              VAOXX(6,ifun,IP) = -sa2ex1f*(z5-Four*Y2*z3+x2*z3)-sa1ex1mt
     1           *(Seven*z3-Twelve*Y2*z+Three*x2*z)-sa0ex1*Six*z
              VAOXXX(1,ifun,IP) = -sa3ex1me*(x3*z3-Four*x3*Y2*z+x5*z)-sa
     1           2ex1f*(Three*x*z3-Twelve*x*Y2*z+Nine*x3*z)-sa1ex1mt*Twe
     2           lve*x*z
              VAOXXX(2,ifun,IP) = -sa3ex1me*(x2*y*z3-Four*x2*Y3*z+x4*y*z
     1           )-sa2ex1f*(y*z3-Four*Y3*z-Three*x2*y*z)+sa1ex1mt*Six*y*
     2           z
              VAOXXX(3,ifun,IP) = -sa3ex1me*(x*Y2*z3-Four*x*y4*z+x3*Y2*z
     1           )-sa2ex1f*(x*z3-TEight*x*Y2*z+x3*z)+sa1ex1mt*Six*x*z
              VAOXXX(4,ifun,IP) = -sa3ex1me*(Y3*z3-Four*y5*z+x2*Y3*z)-sa
     1           2ex1f*(Three*y*z3-ThSix*Y3*z+Three*x2*y*z)+FEight*sa1ex
     2           1mt*y*z
              VAOXXX(5,ifun,IP) = -sa3ex1me*(x2*z4-Four*x2*Y2*z2+x4*z2)-
     1           sa2ex1f*(z4-Four*Y2*z2+Eight*x2*z2-Four*x2*Y2+x4)-sa1ex
     2           1mt*(Five*z2-Four*Y2+Five*x2)-sa0ex1*Two
              VAOXXX(6,ifun,IP) = -sa3ex1me*(x*y*z4-Four*x*Y3*z2+x3*y*z2
     1           )+sa2ex1f*(Three*x*y*z2+Four*x*Y3-x3*y)+sa1ex1mt*Six*x*
     2           y
              VAOXXX(7,ifun,IP) = -sa3ex1me*(Y2*z4-Four*y4*z2+x2*Y2*z2)-
     1           sa2ex1f*(z4-TSeven*Y2*z2+x2*z2-Four*y4+x2*Y2)+sa1ex1mt*
     2           (Five*z2+Twenty*Y2-x2)+Eight*sa0ex1
              VAOXXX(8,ifun,IP) = -sa3ex1me*(x*z5-Four*x*Y2*z3+x3*z3)-sa
     1           2ex1f*(Nine*x*z3-Twelve*x*Y2*z+Three*x3*z)-sa1ex1mt*Twe
     2           lve*x*z
              VAOXXX(9,ifun,IP) = -sa3ex1me*(y*z5-Four*Y3*z3+x2*y*z3)+sa
     1           2ex1f*(y*z3+Twelve*Y3*z-Three*x2*y*z)+sa1ex1mt*TEight*y
     2           *z
              VAOXXX(10,ifun,IP) = -sa3ex1me*(z6-Four*Y2*z4+x2*z4)-sa2ex
     1           1f*(Twelve*z4-TwFour*Y2*z2+Six*x2*z2)-sa1ex1mt*(TwSeven
     2           *z2-Twelve*Y2+Three*x2)-sa0ex1*Six
c
c   component  [4,[1,0,2],-1,[3,0,0],-1,[1,2,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*(Four*x*z2-x*Y2-x3)
              VAOX(1,ifun,IP) = sa1ex1mt*(Four*x2*z2-x2*Y2-x4)+sa0ex1*(F
     1           our*z2-Y2-Three*x2)
              VAOX(2,ifun,IP) = sa1ex1mt*(Four*x*y*z2-x*Y3-x3*y)-sa0ex1*
     1           Two*x*y
              VAOX(3,ifun,IP) = sa1ex1mt*(Four*x*z3-x*Y2*z-x3*z)+Eight*s
     1           a0ex1*x*z
              VAOXX(1,ifun,IP) = sa2ex1f*(Four*x3*z2-x3*Y2-x5)+sa1ex1mt*
     1           (Twelve*x*z2-Three*x*Y2-Seven*x3)-sa0ex1*Six*x
              VAOXX(2,ifun,IP) = sa2ex1f*(Four*x2*y*z2-x2*Y3-x4*y)+sa1ex
     1           1mt*(Four*y*z2-Y3-Five*x2*y)-sa0ex1*Two*y
              VAOXX(3,ifun,IP) = sa2ex1f*(Four*x*Y2*z2-x*y4-x3*Y2)+sa1ex
     1           1mt*(Four*x*z2-Five*x*Y2-x3)-sa0ex1*Two*x
              VAOXX(4,ifun,IP) = sa2ex1f*(Four*x2*z3-x2*Y2*z-x4*z)+sa1ex
     1           1mt*(Four*z3-Y2*z+Five*x2*z)+Eight*sa0ex1*z
              VAOXX(5,ifun,IP) = sa2ex1f*(Four*x*y*z3-x*Y3*z-x3*y*z)+sa1
     1           ex1mt*Six*x*y*z
              VAOXX(6,ifun,IP) = sa2ex1f*(Four*x*z4-x*Y2*z2-x3*z2)+sa1ex
     1           1mt*(Twenty*x*z2-x*Y2-x3)+Eight*sa0ex1*x
              VAOXXX(1,ifun,IP) = sa3ex1me*(Four*x4*z2-x4*Y2-x6)+sa2ex1f
     1           *(TwFour*x2*z2-Six*x2*Y2-Twelve*x4)+sa1ex1mt*(Twelve*z2
     2           -Three*Y2-TwSeven*x2)-sa0ex1*Six
              VAOXXX(2,ifun,IP) = sa3ex1me*(Four*x3*y*z2-x3*Y3-x5*y)+sa2
     1           ex1f*(Twelve*x*y*z2-Three*x*Y3-Nine*x3*y)-sa1ex1mt*Twel
     2           ve*x*y
              VAOXXX(3,ifun,IP) = sa3ex1me*(Four*x2*Y2*z2-x2*y4-x4*Y2)+s
     1           a2ex1f*(Four*Y2*z2+Four*x2*z2-y4-Eight*x2*Y2-x4)+sa1ex1
     2           mt*(Four*z2-Five*Y2-Five*x2)-sa0ex1*Two
              VAOXXX(4,ifun,IP) = sa3ex1me*(Four*x*Y3*z2-x*y5-x3*Y3)+sa2
     1           ex1f*(Twelve*x*y*z2-Nine*x*Y3-Three*x3*y)-sa1ex1mt*Twel
     2           ve*x*y
              VAOXXX(5,ifun,IP) = sa3ex1me*(Four*x3*z3-x3*Y2*z-x5*z)+sa2
     1           ex1f*(Twelve*x*z3-Three*x*Y2*z+x3*z)+sa1ex1mt*TEight*x*
     2           z
              VAOXXX(6,ifun,IP) = sa3ex1me*(Four*x2*y*z3-x2*Y3*z-x4*y*z)
     1           +sa2ex1f*(Four*y*z3-Y3*z+Three*x2*y*z)+sa1ex1mt*Six*y*z
              VAOXXX(7,ifun,IP) = sa3ex1me*(Four*x*Y2*z3-x*y4*z-x3*Y2*z)
     1           +sa2ex1f*(Four*x*z3+Three*x*Y2*z-x3*z)+sa1ex1mt*Six*x*z
              VAOXXX(8,ifun,IP) = sa3ex1me*(Four*x2*z4-x2*Y2*z2-x4*z2)+s
     1           a2ex1f*(Four*z4-Y2*z2+TSeven*x2*z2-x2*Y2-x4)+sa1ex1mt*(
     2           Twenty*z2-Y2+Five*x2)+Eight*sa0ex1
              VAOXXX(9,ifun,IP) = sa3ex1me*(Four*x*y*z4-x*Y3*z2-x3*y*z2)
     1           +sa2ex1f*(TEight*x*y*z2-x*Y3-x3*y)+sa1ex1mt*Six*x*y
              VAOXXX(10,ifun,IP) = sa3ex1me*(Four*x*z5-x*Y2*z3-x3*z3)+sa
     1           2ex1f*(ThSix*x*z3-Three*x*Y2*z-Three*x3*z)+FEight*sa1ex
     2           1mt*x*z
c
c   component  [4,[0,1,2],-1,[2,1,0],-1,[0,3,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*(Four*y*z2-Y3-x2*y)
              VAOX(1,ifun,IP) = sa1ex1mt*(Four*x*y*z2-x*Y3-x3*y)-sa0ex1*
     1           Two*x*y
              VAOX(2,ifun,IP) = sa1ex1mt*(Four*Y2*z2-y4-x2*Y2)+sa0ex1*(F
     1           our*z2-Three*Y2-x2)
              VAOX(3,ifun,IP) = sa1ex1mt*(Four*y*z3-Y3*z-x2*y*z)+Eight*s
     1           a0ex1*y*z
              VAOXX(1,ifun,IP) = sa2ex1f*(Four*x2*y*z2-x2*Y3-x4*y)+sa1ex
     1           1mt*(Four*y*z2-Y3-Five*x2*y)-sa0ex1*Two*y
              VAOXX(2,ifun,IP) = sa2ex1f*(Four*x*Y2*z2-x*y4-x3*Y2)+sa1ex
     1           1mt*(Four*x*z2-Five*x*Y2-x3)-sa0ex1*Two*x
              VAOXX(3,ifun,IP) = sa2ex1f*(Four*Y3*z2-y5-x2*Y3)+sa1ex1mt*
     1           (Twelve*y*z2-Seven*Y3-Three*x2*y)-sa0ex1*Six*y
              VAOXX(4,ifun,IP) = sa2ex1f*(Four*x*y*z3-x*Y3*z-x3*y*z)+sa1
     1           ex1mt*Six*x*y*z
              VAOXX(5,ifun,IP) = sa2ex1f*(Four*Y2*z3-y4*z-x2*Y2*z)+sa1ex
     1           1mt*(Four*z3+Five*Y2*z-x2*z)+Eight*sa0ex1*z
              VAOXX(6,ifun,IP) = sa2ex1f*(Four*y*z4-Y3*z2-x2*y*z2)+sa1ex
     1           1mt*(Twenty*y*z2-Y3-x2*y)+Eight*sa0ex1*y
              VAOXXX(1,ifun,IP) = sa3ex1me*(Four*x3*y*z2-x3*Y3-x5*y)+sa2
     1           ex1f*(Twelve*x*y*z2-Three*x*Y3-Nine*x3*y)-sa1ex1mt*Twel
     2           ve*x*y
              VAOXXX(2,ifun,IP) = sa3ex1me*(Four*x2*Y2*z2-x2*y4-x4*Y2)+s
     1           a2ex1f*(Four*Y2*z2+Four*x2*z2-y4-Eight*x2*Y2-x4)+sa1ex1
     2           mt*(Four*z2-Five*Y2-Five*x2)-sa0ex1*Two
              VAOXXX(3,ifun,IP) = sa3ex1me*(Four*x*Y3*z2-x*y5-x3*Y3)+sa2
     1           ex1f*(Twelve*x*y*z2-Nine*x*Y3-Three*x3*y)-sa1ex1mt*Twel
     2           ve*x*y
              VAOXXX(4,ifun,IP) = sa3ex1me*(Four*y4*z2-y6-x2*y4)+sa2ex1f
     1           *(TwFour*Y2*z2-Twelve*y4-Six*x2*Y2)+sa1ex1mt*(Twelve*z2
     2           -TwSeven*Y2-Three*x2)-sa0ex1*Six
              VAOXXX(5,ifun,IP) = sa3ex1me*(Four*x2*y*z3-x2*Y3*z-x4*y*z)
     1           +sa2ex1f*(Four*y*z3-Y3*z+Three*x2*y*z)+sa1ex1mt*Six*y*z
              VAOXXX(6,ifun,IP) = sa3ex1me*(Four*x*Y2*z3-x*y4*z-x3*Y2*z)
     1           +sa2ex1f*(Four*x*z3+Three*x*Y2*z-x3*z)+sa1ex1mt*Six*x*z
              VAOXXX(7,ifun,IP) = sa3ex1me*(Four*Y3*z3-y5*z-x2*Y3*z)+sa2
     1           ex1f*(Twelve*y*z3+Y3*z-Three*x2*y*z)+sa1ex1mt*TEight*y*
     2           z
              VAOXXX(8,ifun,IP) = sa3ex1me*(Four*x*y*z4-x*Y3*z2-x3*y*z2)
     1           +sa2ex1f*(TEight*x*y*z2-x*Y3-x3*y)+sa1ex1mt*Six*x*y
              VAOXXX(9,ifun,IP) = sa3ex1me*(Four*Y2*z4-y4*z2-x2*Y2*z2)+s
     1           a2ex1f*(Four*z4+TSeven*Y2*z2-x2*z2-y4-x2*Y2)+sa1ex1mt*(
     2           Twenty*z2+Five*Y2-x2)+Eight*sa0ex1
              VAOXXX(10,ifun,IP) = sa3ex1me*(Four*y*z5-Y3*z3-x2*y*z3)+sa
     1           2ex1f*(ThSix*y*z3-Three*Y3*z-Three*x2*y*z)+FEight*sa1ex
     2           1mt*y*z
c
c   component  [2*sq10,[1,1,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*sq10*Two*x*y*z
              VAOX(1,ifun,IP) = sq10*Two*(sa1ex1mt*x2+sa0ex1)*y*z
              VAOX(2,ifun,IP) = sq10*Two*x*(sa1ex1mt*Y2+sa0ex1)*z
              VAOX(3,ifun,IP) = sq10*Two*x*y*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = sq10*(sa2ex1f*Two*x3+sa1ex1mt*Six*x)*y*
     1           z
              VAOXX(2,ifun,IP) = sq10*Two*(sa1ex1mt*(Y2+x2)+sa2ex1f*x2*Y
     1           2+sa0ex1)*z
              VAOXX(3,ifun,IP) = sq10*x*(sa2ex1f*Two*Y3+sa1ex1mt*Six*y)*
     1           z
              VAOXX(4,ifun,IP) = sq10*Two*y*(sa1ex1mt*(z2+x2)+sa2ex1f*x2
     1           *z2+sa0ex1)
              VAOXX(5,ifun,IP) = sq10*Two*x*(sa1ex1mt*(z2+Y2)+sa2ex1f*Y2
     1           *z2+sa0ex1)
              VAOXX(6,ifun,IP) = sq10*x*y*(sa2ex1f*Two*z3+sa1ex1mt*Six*z
     1           )
              VAOXXX(1,ifun,IP) = sq10*(sa3ex1me*Two*x4+sa2ex1f*Twelve*x
     1           2+sa1ex1mt*Six)*y*z
              VAOXXX(2,ifun,IP) = sq10*(sa2ex1f*(Six*x*Y2+Two*x3)+sa3ex1
     1           me*Two*x3*Y2+sa1ex1mt*Six*x)*z
              VAOXXX(3,ifun,IP) = sq10*(sa2ex1f*(Two*Y3+Six*x2*y)+sa3ex1
     1           me*Two*x2*Y3+sa1ex1mt*Six*y)*z
              VAOXXX(4,ifun,IP) = sq10*x*(sa3ex1me*Two*y4+sa2ex1f*Twelve
     1           *Y2+sa1ex1mt*Six)*z
              VAOXXX(5,ifun,IP) = sq10*y*(sa2ex1f*(Six*x*z2+Two*x3)+sa3e
     1           x1me*Two*x3*z2+sa1ex1mt*Six*x)
              VAOXXX(6,ifun,IP) = sq10*Two*(sa2ex1f*(Y2*z2+x2*z2+x2*Y2)+
     1           sa1ex1mt*(z2+Y2+x2)+sa3ex1me*x2*Y2*z2+sa0ex1)
              VAOXXX(7,ifun,IP) = sq10*x*(sa2ex1f*(Six*y*z2+Two*Y3)+sa3e
     1           x1me*Two*Y3*z2+sa1ex1mt*Six*y)
              VAOXXX(8,ifun,IP) = sq10*y*(sa2ex1f*(Two*z3+Six*x2*z)+sa3e
     1           x1me*Two*x2*z3+sa1ex1mt*Six*z)
              VAOXXX(9,ifun,IP) = sq10*x*(sa2ex1f*(Two*z3+Six*Y2*z)+sa3e
     1           x1me*Two*Y2*z3+sa1ex1mt*Six*z)
              VAOXXX(10,ifun,IP) = sq10*x*y*(sa3ex1me*Two*z4+sa2ex1f*Twe
     1           lve*z2+sa1ex1mt*Six)
c
            ELSE IF(ityp.EQ.7) THEN
c -- f function (10 components)    xxx, xxy, xxz, xyy, xyz,
c                                  xzz, yyy, yyz, yzz, zzz
c
              x4 = x*x3
              x5 = x*x4
              x6 = x*x5
              y4 = y*y3
              y5 = y*y4
              y6 = y*y5
              z4 = z*z3
              z5 = z*z4
              z6 = z*z5
c
c   component  [1,[3,0,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x3
              VAOX(1,ifun,IP) = sa1ex1mt*x4+sa0ex1*Three*x2
              VAOX(2,ifun,IP) = sa1ex1mt*x3*y
              VAOX(3,ifun,IP) = sa1ex1mt*x3*z
              VAOXX(1,ifun,IP) = sa2ex1f*x5+sa1ex1mt*Seven*x3+sa0ex1*Six
     1           *x
              VAOXX(2,ifun,IP) = (sa2ex1f*x4+sa1ex1mt*Three*x2)*y
              VAOXX(3,ifun,IP) = x3*(sa2ex1f*Y2+sa1ex1mt)
              VAOXX(4,ifun,IP) = (sa2ex1f*x4+sa1ex1mt*Three*x2)*z
              VAOXX(5,ifun,IP) = sa2ex1f*x3*y*z
              VAOXX(6,ifun,IP) = x3*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = sa3ex1me*x6+sa2ex1f*Twelve*x4+sa1ex1mt
     1           *TwSeven*x2+sa0ex1*Six
              VAOXXX(2,ifun,IP) = (sa3ex1me*x5+sa2ex1f*Seven*x3+sa1ex1mt
     1           *Six*x)*y
              VAOXXX(3,ifun,IP) = sa2ex1f*(Three*x2*Y2+x4)+sa3ex1me*x4*Y
     1           2+sa1ex1mt*Three*x2
              VAOXXX(4,ifun,IP) = x3*(sa3ex1me*Y3+sa2ex1f*Three*y)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x5+sa2ex1f*Seven*x3+sa1ex1mt
     1           *Six*x)*z
              VAOXXX(6,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Three*x2)*y*z
              VAOXXX(7,ifun,IP) = x3*(sa3ex1me*Y2+sa2ex1f)*z
              VAOXXX(8,ifun,IP) = sa2ex1f*(Three*x2*z2+x4)+sa3ex1me*x4*z
     1           2+sa1ex1mt*Three*x2
              VAOXXX(9,ifun,IP) = x3*y*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(10,ifun,IP) = x3*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[2,1,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x2*y
              VAOX(1,ifun,IP) = (sa1ex1mt*x3+sa0ex1*Two*x)*y
              VAOX(2,ifun,IP) = x2*(sa1ex1mt*Y2+sa0ex1)
              VAOX(3,ifun,IP) = sa1ex1mt*x2*y*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x4+Five*sa1ex1mt*x2+sa0ex1*Two
     1           )*y
              VAOXX(2,ifun,IP) = sa1ex1mt*(Two*x*Y2+x3)+sa2ex1f*x3*Y2+sa
     1           0ex1*Two*x
              VAOXX(3,ifun,IP) = x2*(sa2ex1f*Y3+sa1ex1mt*Three*y)
              VAOXX(4,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Two*x)*y*z
              VAOXX(5,ifun,IP) = x2*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(6,ifun,IP) = x2*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x5+Nine*sa2ex1f*x3+sa1ex1mt*
     1           Twelve*x)*y
              VAOXXX(2,ifun,IP) = sa2ex1f*(Five*x2*Y2+x4)+sa1ex1mt*(Two*
     1           Y2+Five*x2)+sa3ex1me*x4*Y2+sa0ex1*Two
              VAOXXX(3,ifun,IP) = sa2ex1f*(Two*x*Y3+Three*x3*y)+sa3ex1me
     1           *x3*Y3+sa1ex1mt*Six*x*y
              VAOXXX(4,ifun,IP) = x2*(sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1m
     1           t*Three)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x4+Five*sa2ex1f*x2+sa1ex1mt*
     1           Two)*y*z
              VAOXXX(6,ifun,IP) = (sa2ex1f*(Two*x*Y2+x3)+sa3ex1me*x3*Y2+
     1           sa1ex1mt*Two*x)*z
              VAOXXX(7,ifun,IP) = x2*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(8,ifun,IP) = y*(sa2ex1f*(Two*x*z2+x3)+sa3ex1me*x3*z
     1           2+sa1ex1mt*Two*x)
              VAOXXX(9,ifun,IP) = x2*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1
     1           ex1mt)
              VAOXXX(10,ifun,IP) = x2*y*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[2,0,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x2*z
              VAOX(1,ifun,IP) = (sa1ex1mt*x3+sa0ex1*Two*x)*z
              VAOX(2,ifun,IP) = sa1ex1mt*x2*y*z
              VAOX(3,ifun,IP) = x2*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x4+Five*sa1ex1mt*x2+sa0ex1*Two
     1           )*z
              VAOXX(2,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Two*x)*y*z
              VAOXX(3,ifun,IP) = x2*(sa2ex1f*Y2+sa1ex1mt)*z
              VAOXX(4,ifun,IP) = sa1ex1mt*(Two*x*z2+x3)+sa2ex1f*x3*z2+sa
     1           0ex1*Two*x
              VAOXX(5,ifun,IP) = x2*y*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(6,ifun,IP) = x2*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x5+Nine*sa2ex1f*x3+sa1ex1mt*
     1           Twelve*x)*z
              VAOXXX(2,ifun,IP) = (sa3ex1me*x4+Five*sa2ex1f*x2+sa1ex1mt*
     1           Two)*y*z
              VAOXXX(3,ifun,IP) = (sa2ex1f*(Two*x*Y2+x3)+sa3ex1me*x3*Y2+
     1           sa1ex1mt*Two*x)*z
              VAOXXX(4,ifun,IP) = x2*(sa3ex1me*Y3+sa2ex1f*Three*y)*z
              VAOXXX(5,ifun,IP) = sa2ex1f*(Five*x2*z2+x4)+sa1ex1mt*(Two*
     1           z2+Five*x2)+sa3ex1me*x4*z2+sa0ex1*Two
              VAOXXX(6,ifun,IP) = y*(sa2ex1f*(Two*x*z2+x3)+sa3ex1me*x3*z
     1           2+sa1ex1mt*Two*x)
              VAOXXX(7,ifun,IP) = x2*(sa2ex1f*(z2+Y2)+sa3ex1me*Y2*z2+sa1
     1           ex1mt)
              VAOXXX(8,ifun,IP) = sa2ex1f*(Two*x*z3+Three*x3*z)+sa3ex1me
     1           *x3*z3+sa1ex1mt*Six*x*z
              VAOXXX(9,ifun,IP) = x2*y*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(10,ifun,IP) = x2*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1
     1           mt*Three)
c
c   component  [1,[1,2,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*Y2
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*Y2
              VAOX(2,ifun,IP) = x*(sa1ex1mt*Y3+sa0ex1*Two*y)
              VAOX(3,ifun,IP) = sa1ex1mt*x*Y2*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*Y2
              VAOXX(2,ifun,IP) = sa1ex1mt*(Y3+Two*x2*y)+sa2ex1f*x2*Y3+sa
     1           0ex1*Two*y
              VAOXX(3,ifun,IP) = x*(sa2ex1f*y4+Five*sa1ex1mt*Y2+sa0ex1*T
     1           wo)
              VAOXX(4,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*Y2*z
              VAOXX(5,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Two*y)*z
              VAOXX(6,ifun,IP) = x*Y2*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*Y2
              VAOXXX(2,ifun,IP) = sa2ex1f*(Three*x*Y3+Two*x3*y)+sa3ex1me
     1           *x3*Y3+sa1ex1mt*Six*x*y
              VAOXXX(3,ifun,IP) = sa2ex1f*(y4+Five*x2*Y2)+sa3ex1me*x2*y4
     1           +sa1ex1mt*(Five*Y2+Two*x2)+sa0ex1*Two
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*y5+Nine*sa2ex1f*Y3+sa1ex1m
     1           t*Twelve*y)
              VAOXXX(5,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*Y2*z
              VAOXXX(6,ifun,IP) = (sa2ex1f*(Y3+Two*x2*y)+sa3ex1me*x2*Y3+
     1           sa1ex1mt*Two*y)*z
              VAOXXX(7,ifun,IP) = x*(sa3ex1me*y4+Five*sa2ex1f*Y2+sa1ex1m
     1           t*Two)*z
              VAOXXX(8,ifun,IP) = Y2*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1
     1           ex1mt)
              VAOXXX(9,ifun,IP) = x*(sa2ex1f*(Two*y*z2+Y3)+sa3ex1me*Y3*z
     1           2+sa1ex1mt*Two*y)
              VAOXXX(10,ifun,IP) = x*Y2*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[1,1,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*y*z
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*y*z
              VAOX(2,ifun,IP) = x*(sa1ex1mt*Y2+sa0ex1)*z
              VAOX(3,ifun,IP) = x*y*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*y*z
              VAOXX(2,ifun,IP) = (sa1ex1mt*(Y2+x2)+sa2ex1f*x2*Y2+sa0ex1)
     1           *z
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Three*y)*z
              VAOXX(4,ifun,IP) = y*(sa1ex1mt*(z2+x2)+sa2ex1f*x2*z2+sa0ex
     1           1)
              VAOXX(5,ifun,IP) = x*(sa1ex1mt*(z2+Y2)+sa2ex1f*Y2*z2+sa0ex
     1           1)
              VAOXX(6,ifun,IP) = x*y*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*y*z
              VAOXXX(2,ifun,IP) = (sa2ex1f*(Three*x*Y2+x3)+sa3ex1me*x3*Y
     1           2+sa1ex1mt*Three*x)*z
              VAOXXX(3,ifun,IP) = (sa2ex1f*(Y3+Three*x2*y)+sa3ex1me*x2*Y
     1           3+sa1ex1mt*Three*y)*z
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt
     1           *Three)*z
              VAOXXX(5,ifun,IP) = y*(sa2ex1f*(Three*x*z2+x3)+sa3ex1me*x3
     1           *z2+sa1ex1mt*Three*x)
              VAOXXX(6,ifun,IP) = sa2ex1f*(Y2*z2+x2*z2+x2*Y2)+sa1ex1mt*(
     1           z2+Y2+x2)+sa3ex1me*x2*Y2*z2+sa0ex1
              VAOXXX(7,ifun,IP) = x*(sa2ex1f*(Three*y*z2+Y3)+sa3ex1me*Y3
     1           *z2+sa1ex1mt*Three*y)
              VAOXXX(8,ifun,IP) = y*(sa2ex1f*(z3+Three*x2*z)+sa3ex1me*x2
     1           *z3+sa1ex1mt*Three*z)
              VAOXXX(9,ifun,IP) = x*(sa2ex1f*(z3+Three*Y2*z)+sa3ex1me*Y2
     1           *z3+sa1ex1mt*Three*z)
              VAOXXX(10,ifun,IP) = x*y*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex
     1           1mt*Three)
c
c   component  [1,[1,0,2]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*x*z2
              VAOX(1,ifun,IP) = (sa1ex1mt*x2+sa0ex1)*z2
              VAOX(2,ifun,IP) = sa1ex1mt*x*y*z2
              VAOX(3,ifun,IP) = x*(sa1ex1mt*z3+sa0ex1*Two*z)
              VAOXX(1,ifun,IP) = (sa2ex1f*x3+sa1ex1mt*Three*x)*z2
              VAOXX(2,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z2
              VAOXX(3,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z2
              VAOXX(4,ifun,IP) = sa1ex1mt*(z3+Two*x2*z)+sa2ex1f*x2*z3+sa
     1           0ex1*Two*z
              VAOXX(5,ifun,IP) = x*y*(sa2ex1f*z3+sa1ex1mt*Two*z)
              VAOXX(6,ifun,IP) = x*(sa2ex1f*z4+Five*sa1ex1mt*z2+sa0ex1*T
     1           wo)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x4+sa2ex1f*Six*x2+sa1ex1mt*T
     1           hree)*z2
              VAOXXX(2,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z2
              VAOXXX(3,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z2
              VAOXXX(4,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z2
              VAOXXX(5,ifun,IP) = sa2ex1f*(Three*x*z3+Two*x3*z)+sa3ex1me
     1           *x3*z3+sa1ex1mt*Six*x*z
              VAOXXX(6,ifun,IP) = y*(sa2ex1f*(z3+Two*x2*z)+sa3ex1me*x2*z
     1           3+sa1ex1mt*Two*z)
              VAOXXX(7,ifun,IP) = x*(sa2ex1f*(z3+Two*Y2*z)+sa3ex1me*Y2*z
     1           3+sa1ex1mt*Two*z)
              VAOXXX(8,ifun,IP) = sa2ex1f*(z4+Five*x2*z2)+sa3ex1me*x2*z4
     1           +sa1ex1mt*(Five*z2+Two*x2)+sa0ex1*Two
              VAOXXX(9,ifun,IP) = x*y*(sa3ex1me*z4+Five*sa2ex1f*z2+sa1ex
     1           1mt*Two)
              VAOXXX(10,ifun,IP) = x*(sa3ex1me*z5+Nine*sa2ex1f*z3+sa1ex1
     1           mt*Twelve*z)
c
c   component  [1,[0,3,0]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*Y3
              VAOX(1,ifun,IP) = sa1ex1mt*x*Y3
              VAOX(2,ifun,IP) = sa1ex1mt*y4+sa0ex1*Three*Y2
              VAOX(3,ifun,IP) = sa1ex1mt*Y3*z
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*Y3
              VAOXX(2,ifun,IP) = x*(sa2ex1f*y4+sa1ex1mt*Three*Y2)
              VAOXX(3,ifun,IP) = sa2ex1f*y5+sa1ex1mt*Seven*Y3+sa0ex1*Six
     1           *y
              VAOXX(4,ifun,IP) = sa2ex1f*x*Y3*z
              VAOXX(5,ifun,IP) = (sa2ex1f*y4+sa1ex1mt*Three*Y2)*z
              VAOXX(6,ifun,IP) = Y3*(sa2ex1f*z2+sa1ex1mt)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*Y3
              VAOXXX(2,ifun,IP) = sa2ex1f*(y4+Three*x2*Y2)+sa3ex1me*x2*y
     1           4+sa1ex1mt*Three*Y2
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*y5+sa2ex1f*Seven*Y3+sa1ex1
     1           mt*Six*y)
              VAOXXX(4,ifun,IP) = sa3ex1me*y6+sa2ex1f*Twelve*y4+sa1ex1mt
     1           *TwSeven*Y2+sa0ex1*Six
              VAOXXX(5,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*Y3*z
              VAOXXX(6,ifun,IP) = x*(sa3ex1me*y4+sa2ex1f*Three*Y2)*z
              VAOXXX(7,ifun,IP) = (sa3ex1me*y5+sa2ex1f*Seven*Y3+sa1ex1mt
     1           *Six*y)*z
              VAOXXX(8,ifun,IP) = x*Y3*(sa3ex1me*z2+sa2ex1f)
              VAOXXX(9,ifun,IP) = sa2ex1f*(Three*Y2*z2+y4)+sa3ex1me*y4*z
     1           2+sa1ex1mt*Three*Y2
              VAOXXX(10,ifun,IP) = Y3*(sa3ex1me*z3+sa2ex1f*Three*z)
c
c   component  [1,[0,2,1]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*Y2*z
              VAOX(1,ifun,IP) = sa1ex1mt*x*Y2*z
              VAOX(2,ifun,IP) = (sa1ex1mt*Y3+sa0ex1*Two*y)*z
              VAOX(3,ifun,IP) = Y2*(sa1ex1mt*z2+sa0ex1)
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*Y2*z
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y3+sa1ex1mt*Two*y)*z
              VAOXX(3,ifun,IP) = (sa2ex1f*y4+Five*sa1ex1mt*Y2+sa0ex1*Two
     1           )*z
              VAOXX(4,ifun,IP) = x*Y2*(sa2ex1f*z2+sa1ex1mt)
              VAOXX(5,ifun,IP) = sa1ex1mt*(Two*y*z2+Y3)+sa2ex1f*Y3*z2+sa
     1           0ex1*Two*y
              VAOXX(6,ifun,IP) = Y2*(sa2ex1f*z3+sa1ex1mt*Three*z)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*Y2*z
              VAOXXX(2,ifun,IP) = (sa2ex1f*(Y3+Two*x2*y)+sa3ex1me*x2*Y3+
     1           sa1ex1mt*Two*y)*z
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*y4+Five*sa2ex1f*Y2+sa1ex1m
     1           t*Two)*z
              VAOXXX(4,ifun,IP) = (sa3ex1me*y5+Nine*sa2ex1f*Y3+sa1ex1mt*
     1           Twelve*y)*z
              VAOXXX(5,ifun,IP) = Y2*(sa2ex1f*(z2+x2)+sa3ex1me*x2*z2+sa1
     1           ex1mt)
              VAOXXX(6,ifun,IP) = x*(sa2ex1f*(Two*y*z2+Y3)+sa3ex1me*Y3*z
     1           2+sa1ex1mt*Two*y)
              VAOXXX(7,ifun,IP) = sa2ex1f*(Five*Y2*z2+y4)+sa1ex1mt*(Two*
     1           z2+Five*Y2)+sa3ex1me*y4*z2+sa0ex1*Two
              VAOXXX(8,ifun,IP) = x*Y2*(sa3ex1me*z3+sa2ex1f*Three*z)
              VAOXXX(9,ifun,IP) = sa2ex1f*(Two*y*z3+Three*Y3*z)+sa3ex1me
     1           *Y3*z3+sa1ex1mt*Six*y*z
              VAOXXX(10,ifun,IP) = Y2*(sa3ex1me*z4+sa2ex1f*Six*z2+sa1ex1
     1           mt*Three)
c
c   component  [1,[0,1,2]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*y*z2
              VAOX(1,ifun,IP) = sa1ex1mt*x*y*z2
              VAOX(2,ifun,IP) = (sa1ex1mt*Y2+sa0ex1)*z2
              VAOX(3,ifun,IP) = y*(sa1ex1mt*z3+sa0ex1*Two*z)
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*y*z2
              VAOXX(2,ifun,IP) = x*(sa2ex1f*Y2+sa1ex1mt)*z2
              VAOXX(3,ifun,IP) = (sa2ex1f*Y3+sa1ex1mt*Three*y)*z2
              VAOXX(4,ifun,IP) = x*y*(sa2ex1f*z3+sa1ex1mt*Two*z)
              VAOXX(5,ifun,IP) = sa1ex1mt*(z3+Two*Y2*z)+sa2ex1f*Y2*z3+sa
     1           0ex1*Two*z
              VAOXX(6,ifun,IP) = y*(sa2ex1f*z4+Five*sa1ex1mt*z2+sa0ex1*T
     1           wo)
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*y*z2
              VAOXXX(2,ifun,IP) = (sa2ex1f*(Y2+x2)+sa3ex1me*x2*Y2+sa1ex1
     1           mt)*z2
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y3+sa2ex1f*Three*y)*z2
              VAOXXX(4,ifun,IP) = (sa3ex1me*y4+sa2ex1f*Six*Y2+sa1ex1mt*T
     1           hree)*z2
              VAOXXX(5,ifun,IP) = y*(sa2ex1f*(z3+Two*x2*z)+sa3ex1me*x2*z
     1           3+sa1ex1mt*Two*z)
              VAOXXX(6,ifun,IP) = x*(sa2ex1f*(z3+Two*Y2*z)+sa3ex1me*Y2*z
     1           3+sa1ex1mt*Two*z)
              VAOXXX(7,ifun,IP) = sa2ex1f*(Three*y*z3+Two*Y3*z)+sa3ex1me
     1           *Y3*z3+sa1ex1mt*Six*y*z
              VAOXXX(8,ifun,IP) = x*y*(sa3ex1me*z4+Five*sa2ex1f*z2+sa1ex
     1           1mt*Two)
              VAOXXX(9,ifun,IP) = sa2ex1f*(z4+Five*Y2*z2)+sa3ex1me*Y2*z4
     1           +sa1ex1mt*(Five*z2+Two*Y2)+sa0ex1*Two
              VAOXXX(10,ifun,IP) = y*(sa3ex1me*z5+Nine*sa2ex1f*z3+sa1ex1
     1           mt*Twelve*z)
c
c   component  [1,[0,0,3]]
c
              ifun = ifun+1
              VAO(ifun,IP) = sa0ex1*z3
              VAOX(1,ifun,IP) = sa1ex1mt*x*z3
              VAOX(2,ifun,IP) = sa1ex1mt*y*z3
              VAOX(3,ifun,IP) = sa1ex1mt*z4+sa0ex1*Three*z2
              VAOXX(1,ifun,IP) = (sa2ex1f*x2+sa1ex1mt)*z3
              VAOXX(2,ifun,IP) = sa2ex1f*x*y*z3
              VAOXX(3,ifun,IP) = (sa2ex1f*Y2+sa1ex1mt)*z3
              VAOXX(4,ifun,IP) = x*(sa2ex1f*z4+sa1ex1mt*Three*z2)
              VAOXX(5,ifun,IP) = y*(sa2ex1f*z4+sa1ex1mt*Three*z2)
              VAOXX(6,ifun,IP) = sa2ex1f*z5+sa1ex1mt*Seven*z3+sa0ex1*Six
     1           *z
              VAOXXX(1,ifun,IP) = (sa3ex1me*x3+sa2ex1f*Three*x)*z3
              VAOXXX(2,ifun,IP) = (sa3ex1me*x2+sa2ex1f)*y*z3
              VAOXXX(3,ifun,IP) = x*(sa3ex1me*Y2+sa2ex1f)*z3
              VAOXXX(4,ifun,IP) = (sa3ex1me*Y3+sa2ex1f*Three*y)*z3
              VAOXXX(5,ifun,IP) = sa2ex1f*(z4+Three*x2*z2)+sa3ex1me*x2*z
     1           4+sa1ex1mt*Three*z2
              VAOXXX(6,ifun,IP) = x*y*(sa3ex1me*z4+sa2ex1f*Three*z2)
              VAOXXX(7,ifun,IP) = sa2ex1f*(z4+Three*Y2*z2)+sa3ex1me*Y2*z
     1           4+sa1ex1mt*Three*z2
              VAOXXX(8,ifun,IP) = x*(sa3ex1me*z5+sa2ex1f*Seven*z3+sa1ex1
     1           mt*Six*z)
              VAOXXX(9,ifun,IP) = y*(sa3ex1me*z5+sa2ex1f*Seven*z3+sa1ex1
     1           mt*Six*z)
              VAOXXX(10,ifun,IP) = sa3ex1me*z6+sa2ex1f*Twelve*z4+sa1ex1m
     1           t*TwSeven*z2+sa0ex1*Six
c
            ELSE
c  -- g functions -  currently not active
              call nerror(1,'AODer2:',
     $             'g-functions and above not yet coded in DFT',0,0)
            ENDIF
c -- end of loop over general contractions   (igr loop)
 800        Continue
c -- end of loop over contracted shells      (ics loop)
 900      continue
c -- end of loop over grid points            (IP loop)
 950    CONTINUE
C
      return
      end
c ======================================================================
      SUBROUTINE AODer2(NP,     ncs,    ncf,    XGRID,
     $                  xnuc,   bl,     basdat, inx,    thrsh,
     $                  ExpMIN, VAO,    VAOX,   VAOXX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates vector of basis function (AO) values and their
C  first and second derivatives at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  ncs     -  number of shells
C  ncf     -  number of basis functions
C  XGRID   -  grid points
C  xnuc    -  nuclear coordinates
C  bl      -  array of precomputed normalization factors
C             (from routine <AOInit>)
C  basdat  -  basis set data for TEXAS
C  inx     -  more basis set data
C  thrsh   -  exponent threshold for neglecting contribution
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  VAO     -  on exit contains basis function values
C  VAOX    -  on exit contains basis function derivatives
C  VAOXX   -  on exit contains basis function 2nd derivatives
C
C
      DIMENSION xnuc(3,*),bl(*),basdat(13,*),inx(12,*),ExpMIN(*),
     $          VAO(ncf,NP),VAOX(3,ncf,NP),VAOXX(6,ncf,NP),XGRID(3,NP)
      dimension x0(3)
      parameter (Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0,Three=3.0d0,
     $           Four=4.0d0,Five=5.0d0,Six=6.0d0,Eight=8.0d0)
      parameter (PI=3.14159 26535 89793d0)
C
C
C  get some values
C
        sqtw=One/sqrt(12.0d0)
        sqtw2=two*sqtw
        sqtw4=four*sqtw
        sqr40=sqrt(40.0d0)
c
        DO 950 IP=1,NP
        xp = XGRID(1,IP)
        yp = XGRID(2,IP)
        zp = XGRID(3,IP)
        do 900 ics=1,ncs
c len=number of ang. mom. components, e.g. 3 for P, 4 for L, 5 for D;
c ia=beginning, ie=ending of the contraction
c iat=atom number, xnuc=coordinates of the center(nucleus)
c   ityp: 1=S, 2=P, 3=SP(L), 4=5D, 5=6D, 6=7F, 7=10F
          iat=inx(2,ics)
          x0(1)=xp-xnuc(1,iat)
          x0(2)=yp-xnuc(2,iat)
          x0(3)=zp-xnuc(3,iat)
          r02=x0(1)**2+x0(2)**2+x0(3)**2
cc  ...............................................................
cc    global shell exponent test                ! jb march 97
cc  ...............................................................
          rrexp = r02*ExpMIN(ics)
c          write(6,*) ' ics is:',ics
c          write(6,*) ' rrexp is:',rrexp,' thrsh is:',thrsh
          If(rrexp.GT.thrsh) go to 900
cc  ...............................................................
          ityp=inx(12,ics)                       ! shell type
          ia=inx(1,ics)+1                        ! start of contraction
          ie=inx(5,ics)                          ! end of contraction
          ngr = inx(4,ics)                       ! no. of general contractions
          ifun=inx(11,ics)
c -- general contraction loop
            DO 800 igr=0,ngr
            xnorm=Zero
            xnorm1=Zero
            derx=Zero
            dery=Zero
            derz=Zero
            derxx=zero
            deryy=zero
            derzz=zero
            derxy=zero
            derxz=zero
            deryz=zero
            if(ityp.eq.3) then
                derx1=Zero
                dery1=Zero
                derz1=Zero
                derxx1=zero
                deryy1=zero
                derzz1=zero
                derxy1=zero
                derxz1=zero
                deryz1=zero
              end if
c -- contraction loop
            do 700 ish=ia,ie
c  this the Gaussian exponent times the distance squared, a*r02
              aexp=basdat(1,ish)                 ! gaussian exponent
              ar02=aexp*r02
c  ..........................................................
c  do not calculate very small values
              If(ar02.GT.thrsh) go to 700
c  ..........................................................
c  bl(ish) contains the precomputed value (2*aexp/pi)**0.75
c   aexp is the exponent of the Gaussian
              xx=bl(ish)*exp(-ar02)
              xx1=xx
              xx=basdat(2+igr,ish)*xx
              atwo=-aexp-aexp
              xatwo=xx*atwo
              dx0=atwo*x0(1)
              dy0=atwo*x0(2)
              dz0=atwo*x0(3)
              dx=dx0*xx
              dy=dy0*xx
              dz=dz0*xx
c.................
              dxx=xatwo+dx*dx0
              dyy=xatwo+dy*dy0
              dzz=xatwo+dz*dz0
              dxy=dx*dy0
              dxz=dx*dz0
              dyz=dy*dz0
              if(ityp.eq.3) then
                xx1=xx1*basdat(3,ish)
                xatwo1=xx1*atwo
                dx1=dx0*xx1
                dy1=dy0*xx1
                dz1=dz0*xx1
                dxx1=xatwo1+dx1*dx0
                dyy1=xatwo1+dy1*dy0
                dzz1=xatwo1+dz1*dz0
                dxy1=dx1*dy0
                dxz1=dx1*dz0
                dyz1=dy1*dz0
              end if
c..................
c  d/dx[exp(-a(x-xn)**2-a(y-yn)**2-a(z-zn)**2]= -2a*(x-xn)*exp(...)
c  for L type, calculate another for the p component
c basdat(2,ish) is the main contraction coefficient
c basdat(3,ish) is the P contraction coefficient
              xnorm=xnorm+xx
              derx=derx+dx
              dery=dery+dy
              derz=derz+dz
              derxx=derxx+dxx
              deryy=deryy+dyy
              derzz=derzz+dzz
              derxy=derxy+dxy
              derxz=derxz+dxz
              deryz=deryz+dyz
              if(ityp.eq.3) then         ! segmented basis only
                xnorm1=xnorm1+xx1
                derx1=derx1+dx1
                dery1=dery1+dy1
                derz1=derz1+dz1
                derxx1=derxx1+dxx1
                deryy1=deryy1+dyy1
                derzz1=derzz1+dzz1
                derxy1=derxy1+dxy1
                derxz1=derxz1+dxz1
                deryz1=deryz1+dyz1
              end if
 700        continue
c
c  the spherically symmetrical part of the AO at this point is in
c  xnorm. For L shells, the s component value is in xnorm, the P
c  component in xnorm1. These are the same for all angular momentum
c  components (e.g. X,Y,Z or XX,YY,ZZ,XY,XZ,YZ  etc.)
c  loop over the function components in one contracted shell
c  the radial part is identical
c
c  The spherical part of the gradient is in derx,dery,derz.
c  For L shells, the S component is in derx,dery,derz, the P is in
c  derx1,dery1,derz1
cc
c  The order of the second derivatives is xx, xy, yy, xz, yz, zz
              IF(ityp.EQ.1) THEN
c -- s function
                ifun = ifun+1
                VAO(ifun,IP)    = xnorm
                VAOX(1,ifun,IP) = derx
                VAOX(2,ifun,IP) = dery
                VAOX(3,ifun,IP) = derz
                VAOXX(1,ifun,IP)=derxx
                VAOXX(2,ifun,IP)=derxy
                VAOXX(3,ifun,IP)=deryy
                VAOXX(4,ifun,IP)=derxz
                VAOXX(5,ifun,IP)=deryz
                VAOXX(6,ifun,IP)=derzz
              ELSE IF(ityp.EQ.2) THEN
c -- p function
cc  Loop over the 3 components of a p function
                Do icomp=1,3
                ifun = ifun+1
                VAO(ifun,IP)    = xnorm*x0(icomp)
                VAOX(1,ifun,IP) = derx*x0(icomp)
                VAOX(2,ifun,IP) = dery*x0(icomp)
                VAOX(3,ifun,IP) = derz*x0(icomp)
                VAOXX(1,ifun,IP)=derxx*x0(icomp)
                VAOXX(2,ifun,IP)=derxy*x0(icomp)
                VAOXX(3,ifun,IP)=deryy*x0(icomp)
                VAOXX(4,ifun,IP)=derxz*x0(icomp)
                VAOXX(5,ifun,IP)=deryz*x0(icomp)
                VAOXX(6,ifun,IP)=derzz*x0(icomp)
c  add the derivative of the polynomial term - non-zero only if
c  icomp = derivative direction for the first derivatives
                if(icomp.eq.1) then
c d2/dx2[x*R]=2*dR/dx+x*d2R/dx2
                  VAOXX(1,ifun,IP)=VAOXX(1,ifun,IP)+two*derx
c d2/dxdy[x*R]=dR/dy+x*d2R/dxdy
                  VAOXX(2,ifun,IP)=VAOXX(2,ifun,IP)+dery
c d2/dxdz[x*R]=dR/dz+x*d2R/dxdz
                  VAOXX(4,ifun,IP)=VAOXX(4,ifun,IP)+derz
c rest of the derivatives(d2/dy2, d2/dydz. d2/dz2) are d2R/dy2 etc.
                else if(icomp.eq.2) then
c d2/dxdy[y*R]=dR/dx+y*d2R/dxdy
                  VAOXX(2,ifun,IP)=VAOXX(2,ifun,IP)+derx
c d2/dy2[y*R]=2*dR/dy+y*d2R/dy2
                  VAOXX(3,ifun,IP)=VAOXX(3,ifun,IP)+two*dery
c d2/dydz[y*R]=dR/dz+y*d2R/dydz
                  VAOXX(5,ifun,IP)=VAOXX(5,ifun,IP)+derz
                else if (icomp.eq.3) then
c d2[z*R]/dxdz=dR/dx+z*d2R/dxdz
                  VAOXX(4,ifun,IP)=VAOXX(4,ifun,IP)+derx
c d2[z*R]/dydz=dR/dy+z*d2R/dydz
                  VAOXX(5,ifun,IP)=VAOXX(5,ifun,IP)+dery
c d2[z*R]/dz2=2*dR/dz+z*d2R/dz2
                  VAOXX(6,ifun,IP)=VAOXX(6,ifun,IP)+two*derz
                end if
                VAOX(icomp,ifun,IP) = VAOX(icomp,ifun,IP) + xnorm
                EndDo
              ELSE IF(ityp.EQ.3) THEN
c -- l function (sp)
                ifun = ifun+1
                VAO(ifun,IP)    = xnorm
                VAOX(1,ifun,IP) = derx
                VAOX(2,ifun,IP) = dery
                VAOX(3,ifun,IP) = derz
c
                VAOXX(1,ifun,IP)=derxx
                VAOXX(2,ifun,IP)=derxy
                VAOXX(3,ifun,IP)=deryy
                VAOXX(4,ifun,IP)=derxz
                VAOXX(5,ifun,IP)=deryz
                VAOXX(6,ifun,IP)=derzz
                Do icomp=1,3
                  ifun = ifun+1
                  VAO(ifun,IP)    = xnorm1*x0(icomp)
                  VAOX(1,ifun,IP) = derx1*x0(icomp)
                  VAOX(2,ifun,IP) = dery1*x0(icomp)
                  VAOX(3,ifun,IP) = derz1*x0(icomp)
                  VAOX(icomp,ifun,IP) = VAOX(icomp,ifun,IP) + xnorm1
c
                  VAOXX(1,ifun,IP)=derxx1*x0(icomp)
                  VAOXX(2,ifun,IP)=derxy1*x0(icomp)
                  VAOXX(3,ifun,IP)=deryy1*x0(icomp)
                  VAOXX(4,ifun,IP)=derxz1*x0(icomp)
                  VAOXX(5,ifun,IP)=deryz1*x0(icomp)
                  VAOXX(6,ifun,IP)=derzz1*x0(icomp)
                if(icomp.eq.1) then
                  VAOXX(1,ifun,IP)=VAOXX(1,ifun,IP)+two*derx1
                  VAOXX(2,ifun,IP)=VAOXX(2,ifun,IP)+dery1
                  VAOXX(4,ifun,IP)=VAOXX(4,ifun,IP)+derz1
                else if(icomp.eq.2) then
                  VAOXX(2,ifun,IP)=VAOXX(2,ifun,IP)+derx1
                  VAOXX(3,ifun,IP)=VAOXX(3,ifun,IP)+two*dery1
                  VAOXX(5,ifun,IP)=VAOXX(5,ifun,IP)+derz1
                else if (icomp.eq.3) then
                  VAOXX(4,ifun,IP)=VAOXX(4,ifun,IP)+derx1
                  VAOXX(5,ifun,IP)=VAOXX(5,ifun,IP)+dery1
                  VAOXX(6,ifun,IP)=VAOXX(6,ifun,IP)+two*derz1
                end if
                EndDo
              ELSE IF(ityp.EQ.4) THEN
c -- d function (5 components)    z**2,x**2-y**2,xy,xz,yz
                ifun = ifun+1
                xnorm1=xnorm*x0(1)
                xnorm2=xnorm*x0(2)
                xnorm3=xnorm*x0(3)
c  (2*z**2 -x**2-y**2)/SQRT(12)  = dz2
                polynom = (Two*x0(3)**2-x0(1)**2-x0(2)**2)*sqtw
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx - sqtw2*xnorm1
                VAOX(2,ifun,IP) = polynom*dery - sqtw2*xnorm2
                VAOX(3,ifun,IP) = polynom*derz + sqtw4*xnorm3
c
                VAOXX(1,ifun,IP)=polynom*derxx - sqtw4*x0(1)*derx -
     1            sqtw2*xnorm
                VAOXX(2,ifun,IP) = polynom*derxy -
     1            sqtw2*(x0(2)*derx + x0(1)*dery)
                VAOXX(3,ifun,IP) = polynom*deryy -sqtw4*x0(2)*dery -
     1            sqtw2*xnorm
                VAOXX(4,ifun,IP) = polynom*derxz + sqtw4*x0(3)*derx -
     1            sqtw2*x0(1)*derz
                VAOXX(5,ifun,IP) = polynom*deryz + sqtw4*x0(3)*dery -
     1            sqtw2*x0(2)*derz
                VAOXX(6,ifun,IP) = polynom*derzz +
     1            sqtw4*(two*x0(3)*derz + xnorm)
c
                ifun = ifun+1
c   (x**2-y**2)/2 = d(x2-y2)
                polynom = (x0(1)**2-x0(2)**2)*Half
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm1
                VAOX(2,ifun,IP) = polynom*dery - xnorm2
                VAOX(3,ifun,IP) = polynom*derz
c
                VAOXX(1,ifun,IP)=polynom*derxx + two*x0(1)*derx +
     1            xnorm
                VAOXX(2,ifun,IP) = polynom*derxy -
     1            x0(2)*derx + x0(1)*dery
                VAOXX(3,ifun,IP) = polynom*deryy -two*x0(2)*dery -
     1            xnorm
                VAOXX(4,ifun,IP) = polynom*derxz + x0(1)*derz
                VAOXX(5,ifun,IP) = polynom*deryz - x0(2)*derz
                VAOXX(6,ifun,IP) = polynom*derzz
c
                ifun = ifun+1
c    xy
                polynom = x0(1)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(2)*xnorm
                VAOX(2,ifun,IP) = polynom*dery + x0(1)*xnorm
                VAOX(3,ifun,IP) = polynom*derz
c
                VAOXX(1,ifun,IP)=polynom*derxx+two*x0(2)*derx
                VAOXX(2,ifun,IP)=polynom*derxy+x0(1)*derx+x0(2)*dery +
     1            xnorm
                VAOXX(3,ifun,IP)=polynom*deryy+two*x0(1)*dery
                VAOXX(4,ifun,IP)=polynom*derxz+x0(2)*derz
                VAOXX(5,ifun,IP)=polynom*deryz+x0(1)*derz
                VAOXX(6,ifun,IP)=polynom*derzz
c
                ifun = ifun+1
c    xz
                polynom = x0(1)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(3)*xnorm
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + x0(1)*xnorm
c
                VAOXX(1,ifun,IP)=polynom*derxx+two*x0(3)*derx
                VAOXX(2,ifun,IP)=polynom*derxy+x0(3)*dery
                VAOXX(3,ifun,IP)=polynom*deryy
                VAOXX(4,ifun,IP)=polynom*derxz + x0(3)*derz +
     1             x0(1)*derx + xnorm
                VAOXX(5,ifun,IP)=polynom*deryz+x0(1)*dery   ! this one
                VAOXX(6,ifun,IP)=polynom*derzz+two*x0(1)*derz
c
                ifun = ifun+1
c    yz
                polynom = x0(2)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + x0(3)*xnorm
                VAOX(3,ifun,IP) = polynom*derz + x0(2)*xnorm
c
                VAOXX(1,ifun,IP)=polynom*derxx
                VAOXX(2,ifun,IP)=polynom*derxy+x0(3)*derx
                VAOXX(3,ifun,IP)=polynom*deryy+two*x0(3)*dery
                VAOXX(4,ifun,IP)=polynom*derxz+x0(2)*derx
                VAOXX(5,ifun,IP)=polynom*deryz + x0(2)*dery +
     1            x0(3)*derz + xnorm
                VAOXX(6,ifun,IP)=polynom*derzz+two*x0(2)*derz
c
              ELSE IF(ityp.EQ.5) THEN
c -- d function (6 components)    xx,yy,zz,xy,xz,yz
                ifun = ifun+1
c    xx
                polynom = x0(1)**2
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Two*x0(1)*xnorm
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz
c
                VAOXX(1,ifun,IP) = polynom*derxx +
     1            four*x0(1)*derx+two*xnorm
                VAOXX(2,ifun,IP)=polynom*derxy+two*x0(1)*dery
                VAOXX(3,ifun,IP)=polynom*deryy
                VAOXX(4,ifun,IP)=polynom*derxz+two*x0(1)*derz
                VAOXX(5,ifun,IP)=polynom*deryz
                VAOXX(6,ifun,IP)=polynom*derzz
                ifun = ifun+1
c    yy
                polynom = x0(2)**2
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + Two*x0(2)*xnorm
                VAOX(3,ifun,IP) = polynom*derz
c
                VAOXX(1,ifun,IP)=polynom*derxx
                VAOXX(2,ifun,IP)=polynom*derxy+two*x0(2)*derx
                VAOXX(3,ifun,IP)=polynom*deryy+
     1            four*x0(2)*dery+two*xnorm
                VAOXX(4,ifun,IP)=polynom*derxz
                VAOXX(5,ifun,IP)=polynom*deryz+two*x0(2)*derz
                VAOXX(6,ifun,IP)=polynom*derzz
                ifun = ifun+1
c    zz
                polynom = x0(3)**2
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + Two*x0(3)*xnorm
c
                VAOXX(1,ifun,IP) = polynom*derxx
                VAOXX(2,ifun,IP) = polynom*derxy
                VAOXX(3,ifun,IP) = polynom*deryy
                VAOXX(4,ifun,IP) = polynom*derxz + two*x0(3)*derx
                VAOXX(5,ifun,IP) = polynom*deryz + two*x0(3)*dery
                VAOXX(6,ifun,IP) = polynom*derzz + four*x0(3)*derz +
     1            two*xnorm
                ifun = ifun+1
c    xy
                polynom = x0(1)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(2)*xnorm
                VAOX(2,ifun,IP) = polynom*dery + x0(1)*xnorm
                VAOX(3,ifun,IP) = polynom*derz
c
                VAOXX(1,ifun,IP) = polynom*derxx + two*x0(2)*derx
                VAOXX(2,ifun,IP) = polynom*derxy +
     1             x0(1)*derx + x0(2)*dery + xnorm
                VAOXX(3,ifun,IP)=polynom*deryy+two*x0(1)*dery
                VAOXX(4,ifun,IP) = polynom*derxz + x0(2)*derz
                VAOXX(5,ifun,IP) = polynom*deryz + x0(1)*derz
                VAOXX(6,ifun,IP) = polynom*derzz
                ifun = ifun+1
c    xz
                polynom = x0(1)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(3)*xnorm
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + x0(1)*xnorm
c
                VAOXX(1,ifun,IP)=polynom*derxx+two*x0(3)*derx
                VAOXX(2,ifun,IP)=polynom*derxy+x0(3)*dery
                VAOXX(3,ifun,IP)=polynom*deryy
                VAOXX(4,ifun,IP)=polynom*derxz + x0(3)*derz +
     1             x0(1)*derx + xnorm
                VAOXX(5,ifun,IP)=polynom*deryz+x0(1)*dery   ! this one
                VAOXX(6,ifun,IP)=polynom*derzz+two*x0(1)*derz
                ifun = ifun+1
c    yz
                polynom = x0(2)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + x0(3)*xnorm
                VAOX(3,ifun,IP) = polynom*derz + x0(2)*xnorm
c
                VAOXX(1,ifun,IP)=polynom*derxx
                VAOXX(2,ifun,IP)=polynom*derxy+x0(3)*derx
                VAOXX(3,ifun,IP)=polynom*deryy+two*x0(3)*dery
                VAOXX(4,ifun,IP)=polynom*derxz+x0(2)*derx
                VAOXX(5,ifun,IP)=polynom*deryz + x0(2)*dery +
     1            x0(3)*derz + xnorm
                VAOXX(6,ifun,IP)=polynom*derzz+two*x0(2)*derz
c
              ELSE IF(ityp.EQ.6) THEN
c -- f function (7 components)    (5xxy-rry),(5xxz-rrz),(5yyx-rrx),
c                                 (5yyz-rrz),(5zzx-rrx),(5zzy-rry),
c                                    xyx*SQRT(40)
c
c  precalculate values
                xy=x0(1)*x0(2)*Two
                xz=x0(1)*x0(3)*Two
                yz=x0(2)*x0(3)*Two
                xy8=Four*xy
                xz8=xz*Four
                yz8=yz*Four
                xx=x0(1)**2
                yy=x0(2)**2
                zz=x0(3)**2
                xx4=Four*xx
                yy4=Four*yy
                zz4=Four*zz
                xx3=Three*xx
                yy3=Three*yy
                zz3=Three*zz
c
                ifun = ifun+1
c    5xxy
                polynom = (Five*x0(1)**2-r02)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xy8*xnorm
                VAOX(2,ifun,IP) = polynom*dery + (xx4-yy3-zz)*xnorm
                VAOX(3,ifun,IP) = polynom*derz - yz*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                      Eight*x0(2)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy + derx*(xx4-yy3-zz)
     $                            + Eight*x0(1)*(xnorm + dery*x0(2))
                VAOXX(3,ifun,IP) =  polynom*deryy +
     $                      Two*dery*(xx4-yy3-zz) - Six*xnorm*x0(2)
                VAOXX(4,ifun,IP) = polynom*derxz + Six*derz*x0(1)*x0(2)
                VAOXX(5,ifun,IP) = polynom*deryz + derz*(xx4-yy3-zz)
     $                            - Two*x0(3)*(xnorm + dery*x0(2))
                VAOXX(6,ifun,IP) = polynom*derzz -
     $                              Two*x0(2)*(xnorm + Two*derz*x0(3))
c .....................................................................
                ifun = ifun+1
c    5xxz
                polynom = (Five*x0(1)**2-r02)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xz8*xnorm
                VAOX(2,ifun,IP) = polynom*dery - yz*xnorm
                VAOX(3,ifun,IP) = polynom*derz + (xx4-zz3-yy)*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                      Eight*x0(3)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy + Six*dery*x0(1)*x0(3)
                VAOXX(3,ifun,IP) = polynom*deryy -
     $                              Two*x0(3)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derx*(xx4-zz3-yy)
     $                            + Eight*x0(1)*(xnorm + derz*x0(3))
                VAOXX(5,ifun,IP) = polynom*deryz + dery*(xx4-zz3-yy)
     $                            - Two*x0(2)*(xnorm + derz*x0(3))
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                      Two*derz*(xx4-zz3-yy) - Six*xnorm*x0(3)
c .....................................................................
                ifun = ifun+1
c    5yyx
                polynom = (Five*x0(2)**2-r02)*x0(1)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + (yy4-xx3-zz)*xnorm
                VAOX(2,ifun,IP) = polynom*dery + xy8*xnorm
                VAOX(3,ifun,IP) = polynom*derz - xz*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) =  polynom*derxx +
     $                      Two*derx*(yy4-xx3-zz) - Six*xnorm*x0(1)
                VAOXX(2,ifun,IP) = polynom*derxy + dery*(yy4-xx3-zz)
     $                            + Eight*x0(2)*(xnorm + derx*x0(1))
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                      Eight*x0(1)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derz*(yy4-xx3-zz)
     $                            - Two*x0(3)*(xnorm + derx*x0(1))
                VAOXX(5,ifun,IP) = polynom*deryz + Six*derz*x0(1)*x0(2)
                VAOXX(6,ifun,IP) = polynom*derzz -
     $                              Two*x0(1)*(xnorm + Two*derz*x0(3))
c .....................................................................
                ifun = ifun+1
c    5yyz
                polynom = (Five*x0(2)**2-r02)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx - xz*xnorm
                VAOX(2,ifun,IP) = polynom*dery + yz8*xnorm
                VAOX(3,ifun,IP) = polynom*derz + (yy4-zz3-xx)*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) = polynom*derxx -
     $                              Two*x0(3)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy + Six*derx*x0(2)*x0(3)
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                      Eight*x0(3)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derx*(yy4-zz3-xx)
     $                            - Two*x0(1)*(xnorm + derz*x0(3))
                VAOXX(5,ifun,IP) = polynom*deryz + dery*(yy4-zz3-xx)
     $                            + Eight*x0(2)*(xnorm + derz*x0(3))
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                      Two*derz*(yy4-zz3-xx) - Six*xnorm*x0(3)
c .....................................................................
                ifun = ifun+1
c    5zzx
                polynom = (Five*x0(3)**2-r02)*x0(1)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + (zz4-xx3-yy)*xnorm
                VAOX(2,ifun,IP) = polynom*dery - xy*xnorm
                VAOX(3,ifun,IP) = polynom*derz + xz8*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) =  polynom*derxx +
     $                      Two*derx*(zz4-xx3-yy) - Six*xnorm*x0(1)
                VAOXX(2,ifun,IP) = polynom*derxy + dery*(zz4-xx3-yy)
     $                            - Two*x0(2)*(xnorm + derx*x0(1))
                VAOXX(3,ifun,IP) = polynom*deryy -
     $                              Two*x0(1)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derz*(zz4-xx3-yy)
     $                            + Eight*x0(3)*(xnorm + derx*x0(1))
                VAOXX(5,ifun,IP) = polynom*deryz + Six*dery*x0(1)*x0(3)
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                      Eight*x0(1)*(xnorm + Two*derz*x0(3))
c .....................................................................
                ifun = ifun+1
c    5zzy
                polynom = (Five*x0(3)**2-r02)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx - xy*xnorm
                VAOX(2,ifun,IP) = polynom*dery + (zz4-yy3-xx)*xnorm
                VAOX(3,ifun,IP) = polynom*derz + yz8*xnorm
c .....................................................................
                VAOXX(1,ifun,IP) = polynom*derxx -
     $                              Two*x0(2)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy + derx*(zz4-yy3-xx)
     $                            - Two*x0(1)*(xnorm + dery*x0(2))
                VAOXX(3,ifun,IP) =  polynom*deryy +
     $                      Two*dery*(zz4-yy3-xx) - Six*xnorm*x0(2)
                VAOXX(4,ifun,IP) = polynom*derxz + Six*derx*x0(2)*x0(3)
                VAOXX(5,ifun,IP) = polynom*deryz + derz*(zz4-yy3-xx)
     $                            + Eight*x0(3)*(xnorm + dery*x0(2))
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                      Eight*x0(2)*(xnorm + Two*derz*x0(3))
c .....................................................................
                ifun = ifun+1
c    xyz
                polynom = x0(1)*x0(2)*x0(3)*sqr40
                xnorm1 = xnorm*sqr40
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(2)*x0(3)*xnorm1
                VAOX(2,ifun,IP) = polynom*dery + x0(1)*x0(3)*xnorm1
                VAOX(3,ifun,IP) = polynom*derz + x0(1)*x0(2)*xnorm1
c .....................................................................
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                              sqr40*Two*derx*x0(2)*x0(3)
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                   sqr40*x0(3)*(xnorm + derx*x0(1) + dery*x0(2))
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                              sqr40*Two*dery*x0(1)*x0(3)
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                   sqr40*x0(2)*(xnorm + derx*x0(1) + derz*x0(3))
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                   sqr40*x0(1)*(xnorm + dery*x0(2) + derz*x0(3))
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                              sqr40*Two*derz*x0(1)*x0(2)
c .....................................................................
              ELSE IF(ityp.EQ.7) THEN
c -- f function (10 components)    xxx, xxy, xxz, xyy, xyz,
c                                  xzz, yyy, yyz, yzz, zzz
                ifun = ifun+1
c    xxx
                polynom = x0(1)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Three*xnorm*x0(1)**2
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                            Six*x0(1)*(xnorm + derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                              Three*derx*x0(1)*x0(2)
                VAOXX(3,ifun,IP) = polynom*deryy
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                              Three*derx*x0(1)*x0(3)
                VAOXX(5,ifun,IP) = polynom*deryz
                VAOXX(6,ifun,IP) = polynom*derzz
c ......................................................................
                ifun = ifun+1
c    xxy
                polynom = x0(2)*x0(1)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Two*xnorm*x0(1)*x0(2)
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(1)**2
                VAOX(3,ifun,IP) = polynom*derz
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                              Two*x0(2)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                      x0(1)*(Two*(xnorm+dery*x0(2)) + derx*x0(1))
                VAOXX(3,ifun,IP) = polynom*deryy + Two*dery*x0(1)**2
                VAOXX(4,ifun,IP) = polynom*derxz + Two*derz*x0(1)*x0(2)
                VAOXX(5,ifun,IP) = polynom*deryz + derz*x0(1)**2
                VAOXX(6,ifun,IP) = polynom*derzz
c ........................................................................
                ifun = ifun+1
c    xxz
                polynom = x0(3)*x0(1)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Two*xnorm*x0(1)*x0(3)
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(1)**2
                VAOXX(1,ifun,IP) = polynom*derxx +
     $                              Two*x0(3)*(xnorm + Two*derx*x0(1))
                VAOXX(2,ifun,IP) = polynom*derxy + Two*dery*x0(1)*x0(3)
                VAOXX(3,ifun,IP) = polynom*deryy
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                      x0(1)*(Two*(xnorm+derz*x0(3)) + derx*x0(1))
                VAOXX(5,ifun,IP) = polynom*deryz + dery*x0(1)**2
                VAOXX(6,ifun,IP) = polynom*derzz + Two*derz*x0(1)**2
c .......................................................................
                ifun = ifun+1
c    xyy
                polynom = x0(1)*x0(2)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(2)**2
                VAOX(2,ifun,IP) = polynom*dery + Two*xnorm*x0(1)*x0(2)
                VAOX(3,ifun,IP) = polynom*derz
                VAOXX(1,ifun,IP) = polynom*derxx + Two*derx*x0(2)**2
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                      x0(2)*(Two*(xnorm+derx*x0(1)) + dery*x0(2))
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                              Two*x0(1)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derz*x0(2)**2
                VAOXX(5,ifun,IP) = polynom*deryz + Two*derz*x0(1)*x0(2)
                VAOXX(6,ifun,IP) = polynom*derzz
c .......................................................................
                ifun = ifun+1
c    xyz
                polynom = x0(1)*x0(2)*x0(3)
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(2)*x0(3)
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(1)*x0(3)
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(1)*x0(2)
                VAOXX(1,ifun,IP) = polynom*derxx + Two*derx*x0(2)*x0(3)
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                       x0(3)*(xnorm + derx*x0(1) + dery*x0(2))
                VAOXX(3,ifun,IP) = polynom*deryy + Two*dery*x0(1)*x0(3)
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                       x0(2)*(xnorm + derx*x0(1) + derz*x0(3))
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                       x0(1)*(xnorm + dery*x0(2) + derz*x0(3))
                VAOXX(6,ifun,IP) = polynom*derzz + Two*derz*x0(1)*x0(2)
c .......................................................................
                ifun = ifun+1
c    xzz
                polynom = x0(1)*x0(3)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(3)**2
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + Two*xnorm*x0(1)*x0(3)
                VAOXX(1,ifun,IP) = polynom*derxx + Two*derx*x0(3)**2
                VAOXX(2,ifun,IP) = polynom*derxy + dery*x0(3)**2
                VAOXX(3,ifun,IP) = polynom*deryy
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                      x0(3)*(Two*(xnorm+derx*x0(1)) + derz*x0(3))
                VAOXX(5,ifun,IP) = polynom*deryz + Two*dery*x0(1)*x0(3)
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                              Two*x0(1)*(xnorm + Two*derz*x0(3))
c ........................................................................
                ifun = ifun+1
c    yyy
                polynom = x0(2)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + Three*xnorm*x0(2)**2
                VAOX(3,ifun,IP) = polynom*derz
                VAOXX(1,ifun,IP) = polynom*derxx
                VAOXX(2,ifun,IP) = polynom*derxy +
     $                              Three*dery*x0(1)*x0(2)
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                            Six*x0(2)*(xnorm + dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                              Three*dery*x0(2)*x0(3)
                VAOXX(6,ifun,IP) = polynom*derzz
c .......................................................................
                ifun = ifun+1
c    yyz
                polynom = x0(3)*x0(2)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + Two*xnorm*x0(2)*x0(3)
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(2)**2
                VAOXX(1,ifun,IP) = polynom*derxx
                VAOXX(2,ifun,IP) = polynom*derxy + Two*derx*x0(2)*x0(3)
                VAOXX(3,ifun,IP) = polynom*deryy +
     $                              Two*x0(3)*(xnorm + Two*dery*x0(2))
                VAOXX(4,ifun,IP) = polynom*derxz + derx*x0(2)**2
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                      x0(2)*(Two*(xnorm+derz*x0(3)) + dery*x0(2))
                VAOXX(6,ifun,IP) = polynom*derzz + Two*derz*x0(2)**2
c ........................................................................
                ifun = ifun+1
c    yzz
                polynom = x0(2)*x0(3)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(3)**2
                VAOX(3,ifun,IP) = polynom*derz + Two*xnorm*x0(2)*x0(3)
                VAOXX(1,ifun,IP) = polynom*derxx
                VAOXX(2,ifun,IP) = polynom*derxy + derx*x0(3)**2
                VAOXX(3,ifun,IP) = polynom*deryy + Two*dery*x0(3)**2
                VAOXX(4,ifun,IP) = polynom*derxz + Two*derx*x0(2)*x0(3)
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                      x0(3)*(Two*(xnorm+dery*x0(2)) + derz*x0(3))
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                              Two*x0(2)*(xnorm + Two*derz*x0(3))
c .........................................................................
                ifun = ifun+1
c    zzz
                polynom = x0(3)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + Three*xnorm*x0(3)**2
                VAOXX(1,ifun,IP) = polynom*derxx
                VAOXX(2,ifun,IP) = polynom*derxy
                VAOXX(3,ifun,IP) = polynom*deryy
                VAOXX(4,ifun,IP) = polynom*derxz +
     $                              Three*derz*x0(1)*x0(3)
                VAOXX(5,ifun,IP) = polynom*deryz +
     $                              Three*derz*x0(2)*x0(3)
                VAOXX(6,ifun,IP) = polynom*derzz +
     $                            Six*x0(3)*(xnorm + derz*x0(3))
c
              ELSE
c  -- g functions -  currently not active
              call nerror(1,'AODer2:',
     $             'g-functions and above not yet coded in DFT',0,0)
              ENDIF
c -- end of loop over general contractions   (igr loop)
 800          Continue
c -- end of loop over contracted shells      (ics loop)
 900          continue
c -- end of loop over grid points            (IP loop)
 950      CONTINUE
C
      return
      end
c ======================================================================
      SUBROUTINE AOGrad(NP,     ncs,    ncf,    XGRID,
     $                  xnuc,   bl,     basdat, inx,    thrsh,
     $                  ExpMIN, VAO,    VAOX)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates vector of basis function (AO) values and their
C  first derivatives at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  ncs     -  number of shells
C  ncf     -  number of basis functions
C  XGRID   -  grid points
C  xnuc    -  nuclear coordinates
C  bl      -  array of precomputed normalization factors
C             (from routine <AOInit>)
C  basdat  -  basis set data for TEXAS
C  inx     -  more basis set data
C  thrsh   -  exponent threshold for neglecting contribution
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  VAO     -  on exit contains basis function values
C  VAOX    -  on exit contains basis function derivatives
C
C
      DIMENSION xnuc(3,*),bl(*),basdat(13,*),inx(12,*),ExpMIN(*),
     $          VAO(ncf,NP),VAOX(3,ncf,NP),XGRID(3,NP)
cc
      dimension x0(3)
      parameter (Zero=0.0d0,One=1.0d0,Two=2.0d0,Three=3.0d0,
     $           Four=4.0d0,Five=5.0d0,Half=0.5d0,Eight=8.0d0)
      parameter (PI=3.14159 26535 89793d0)
C
C
C  get some values
C
        sqtw=One/sqrt(12.0d0)
        sqr40=sqrt(40.0d0)
c
        DO 950 IP=1,NP
        xp = XGRID(1,IP)
        yp = XGRID(2,IP)
        zp = XGRID(3,IP)
        do 900 ics=1,ncs
c len=number of ang. mom. components, e.g. 3 for P, 4 for L, 5 for D;
c ia=beginning, ie=ending of the contraction
c iat=atom number, xnuc=coordinates of the center(nucleus)
c   ityp: 1=S, 2=P, 3=SP(L), 4=5D, 5=6D, 6=7F, 7=10F
          iat=inx(2,ics)
          x0(1)=xp-xnuc(1,iat)
          x0(2)=yp-xnuc(2,iat)
          x0(3)=zp-xnuc(3,iat)
          r02=x0(1)**2+x0(2)**2+x0(3)**2
cc  ...............................................................
cc    global shell exponent test                ! jb march 97
cc  ...............................................................
          rrexp = r02*ExpMIN(ics)
c          write(6,*) ' ics is:',ics
c          write(6,*) ' rrexp is:',rrexp,' thrsh is:',thrsh
          If(rrexp.GT.thrsh) go to 900
cc  ...............................................................
          ityp=inx(12,ics)                       ! shell type
          ia=inx(1,ics)+1                        ! start of contraction
          ie=inx(5,ics)                          ! end of contraction
          ngr = inx(4,ics)                       ! no. of general contractions
          ifun=inx(11,ics)
c -- general contraction loop
            DO 800 igr=0,ngr
            xnorm=Zero
            xnorm1=Zero
            derx=Zero
            dery=Zero
            derz=Zero
            derx1=Zero
            dery1=Zero
            derz1=Zero
c -- contraction loop
            do 700 ish=ia,ie
c  this the Gaussian exponent times the distance squared, a*r02
              aexp=basdat(1,ish)                 ! gaussian exponent
              ar02=aexp*r02
c  ..........................................................
c  do not calculate very small values
              If(ar02.GT.thrsh) go to 700
c  ..........................................................
c  bl(ish) contains the precomputed value (2*aexp/pi)**0.75
c   aexp is the exponent of the Gaussian
              xx=bl(ish)*exp(-ar02)
              atwo=-aexp-aexp
              dx=xx*atwo*x0(1)
              dy=xx*atwo*x0(2)
              dz=xx*atwo*x0(3)
c  d/dx[exp(-a(x-xn)**2-a(y-yn)**2-a(z-zn)**2]= -2a*(x-xn)*exp(...)
c  for L type, calculate another for the p component
c basdat(2,ish) is the main contraction coefficient
c basdat(3,ish) is the P contraction coefficient
              contrs=basdat(2+igr,ish)
              xnorm=xnorm+xx*contrs
              derx=derx+dx*contrs
              dery=dery+dy*contrs
              derz=derz+dz*contrs
              if(ityp.eq.3) then                 ! segmented basis only
                contrp=basdat(3,ish)
                xnorm1=xnorm1+xx*contrp
                derx1=derx1+dx*contrp
                dery1=dery1+dy*contrp
                derz1=derz1+dz*contrp
              end if
 700        continue
c
c  the spherically symmetrical part of the AO at this point is in
c  xnorm. For L shells, the s component value is in xnorm, the P
c  component in xnorm1. These are the same for all angular momentum
c  components (e.g. X,Y,Z or XX,YY,ZZ,XY,XZ,YZ  etc.)
c  loop over the function components in one contracted shell
c  the radial part is identical
c
c  The spherical part of the gradient is in derx,dery,derz.
c  For L shells, the S component is in derx,dery,derz, the P is in
c  derx1,dery1,derz1
cc
              IF(ityp.EQ.1) THEN
c -- s function
                ifun = ifun+1
                VAO(ifun,IP)   = xnorm
                VAOX(1,ifun,IP) = derx
                VAOX(2,ifun,IP) = dery
                VAOX(3,ifun,IP) = derz
              ELSE IF(ityp.EQ.2) THEN
c -- p function
                Do icomp=1,3
                ifun = ifun+1
                VAO(ifun,IP)   = xnorm*x0(icomp)
                VAOX(1,ifun,IP) = derx*x0(icomp)
                VAOX(2,ifun,IP) = dery*x0(icomp)
                VAOX(3,ifun,IP) = derz*x0(icomp)
c  add the derivative of the polynomial term - non-zero only if
c  icomp = derivative direction
                VAOX(icomp,ifun,IP) = VAOX(icomp,ifun,IP) + xnorm
                EndDo
              ELSE IF(ityp.EQ.3) THEN
c -- l function (sp)
                ifun = ifun+1
                VAO(ifun,IP)   = xnorm
                VAOX(1,ifun,IP) = derx
                VAOX(2,ifun,IP) = dery
                VAOX(3,ifun,IP) = derz
                Do icomp=1,3
                ifun = ifun+1
                VAO(ifun,IP)   = xnorm1*x0(icomp)
                VAOX(1,ifun,IP) = derx1*x0(icomp)
                VAOX(2,ifun,IP) = dery1*x0(icomp)
                VAOX(3,ifun,IP) = derz1*x0(icomp)
                VAOX(icomp,ifun,IP) = VAOX(icomp,ifun,IP) + xnorm1
                EndDo
              ELSE IF(ityp.EQ.4) THEN
c -- d function (5 components)    z**2,x**2-y**2,xy,xz,yz
                ifun = ifun+1
c  (2*z**2 -x**2-y**2)/SQRT(12)  = dz2
                polinom = (Two*x0(3)**2-x0(1)**2-x0(2)**2)*sqtw
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx - Two*x0(1)*sqtw*xnorm
                VAOX(2,ifun,IP) = polinom*dery - Two*x0(2)*sqtw*xnorm
                VAOX(3,ifun,IP) = polinom*derz + Four*x0(3)*sqtw*xnorm
                ifun = ifun+1
c   (x**2-y**2)/2 = d(x2-y2)
                polinom = (x0(1)**2-x0(2)**2)*Half
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + x0(1)*xnorm
                VAOX(2,ifun,IP) = polinom*dery - x0(2)*xnorm
                VAOX(3,ifun,IP) = polinom*derz
                ifun = ifun+1
c    xy
                polinom = x0(1)*x0(2)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + x0(2)*xnorm
                VAOX(2,ifun,IP) = polinom*dery + x0(1)*xnorm
                VAOX(3,ifun,IP) = polinom*derz
                ifun = ifun+1
c    xz
                polinom = x0(1)*x0(3)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + x0(3)*xnorm
                VAOX(2,ifun,IP) = polinom*dery
                VAOX(3,ifun,IP) = polinom*derz + x0(1)*xnorm
                ifun = ifun+1
c    yz
                polinom = x0(2)*x0(3)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx
                VAOX(2,ifun,IP) = polinom*dery + x0(3)*xnorm
                VAOX(3,ifun,IP) = polinom*derz + x0(2)*xnorm
              ELSE IF(ityp.EQ.5) THEN
c -- d function (6 components)    xx,yy,zz,xy,xz,yz
                ifun = ifun+1
c    xx
                polinom = x0(1)*x0(1)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + Two*x0(1)*xnorm
                VAOX(2,ifun,IP) = polinom*dery
                VAOX(3,ifun,IP) = polinom*derz
                ifun = ifun+1
c    yy
                polinom = x0(2)*x0(2)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx
                VAOX(2,ifun,IP) = polinom*dery + Two*x0(2)*xnorm
                VAOX(3,ifun,IP) = polinom*derz
                ifun = ifun+1
c    zz
                polinom = x0(3)*x0(3)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx
                VAOX(2,ifun,IP) = polinom*dery
                VAOX(3,ifun,IP) = polinom*derz + Two*x0(3)*xnorm
                ifun = ifun+1
c    xy
                polinom = x0(1)*x0(2)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + x0(2)*xnorm
                VAOX(2,ifun,IP) = polinom*dery + x0(1)*xnorm
                VAOX(3,ifun,IP) = polinom*derz
                ifun = ifun+1
c    xz
                polinom = x0(1)*x0(3)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx + x0(3)*xnorm
                VAOX(2,ifun,IP) = polinom*dery
                VAOX(3,ifun,IP) = polinom*derz + x0(1)*xnorm
                ifun = ifun+1
c    yz
                polinom = x0(2)*x0(3)
                VAO(ifun,IP)   = xnorm*polinom
                VAOX(1,ifun,IP) = polinom*derx
                VAOX(2,ifun,IP) = polinom*dery + x0(3)*xnorm
                VAOX(3,ifun,IP) = polinom*derz + x0(2)*xnorm
              ELSE IF(ityp.EQ.6) THEN
c -- f function (7 components)    (5xxy-rry),(5xxz-rrz),(5yyx-rrx),
c                                 (5yyz-rrz),(5zzx-rrx),(5zzy-rry),
c                                    xyx*SQRT(40)
c
c  precalculate values
                xy=x0(1)*x0(2)*Two
                xz=x0(1)*x0(3)*Two
                yz=x0(2)*x0(3)*Two
                xy8=Four*xy
                xz8=xz*Four
                yz8=yz*Four
                xx=x0(1)**2
                yy=x0(2)**2
                zz=x0(3)**2
                xx4=Four*xx
                yy4=Four*yy
                zz4=Four*zz
                xx3=Three*xx
                yy3=Three*yy
                zz3=Three*zz
c
                ifun = ifun+1
c    5xxy
                polynom = (Five*x0(1)**2-r02)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xy8*xnorm
                VAOX(2,ifun,IP) = polynom*dery + (xx4-yy3-zz)*xnorm
                VAOX(3,ifun,IP) = polynom*derz - yz*xnorm
                ifun = ifun+1
c    5xxz
                polynom = (Five*x0(1)**2-r02)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xz8*xnorm
                VAOX(2,ifun,IP) = polynom*dery - yz*xnorm
                VAOX(3,ifun,IP) = polynom*derz + (xx4-zz3-yy)*xnorm
                ifun = ifun+1
c    5yyx
                polynom = (Five*x0(2)**2-r02)*x0(1)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + (yy4-xx3-zz)*xnorm
                VAOX(2,ifun,IP) = polynom*dery + xy8*xnorm
                VAOX(3,ifun,IP) = polynom*derz - xz*xnorm
                ifun = ifun+1
c    5yyz
                polynom = (Five*x0(2)**2-r02)*x0(3)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx - xz*xnorm
                VAOX(2,ifun,IP) = polynom*dery + yz8*xnorm
                VAOX(3,ifun,IP) = polynom*derz + (yy4-zz3-xx)*xnorm
                ifun = ifun+1
c    5zzx
                polynom = (Five*x0(3)**2-r02)*x0(1)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + (zz4-xx3-yy)*xnorm
                VAOX(2,ifun,IP) = polynom*dery - xy*xnorm
                VAOX(3,ifun,IP) = polynom*derz + xz8*xnorm
                ifun = ifun+1
c    5zzy
                polynom = (Five*x0(3)**2-r02)*x0(2)
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx - xy*xnorm
                VAOX(2,ifun,IP) = polynom*dery + (zz4-yy3-xx)*xnorm
                VAOX(3,ifun,IP) = polynom*derz + yz8*xnorm
                ifun = ifun+1
c    xyz
                polynom = x0(1)*x0(2)*x0(3)*sqr40
                xnorm1 = xnorm*sqr40
                VAO(ifun,IP)   = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + x0(2)*x0(3)*xnorm1
                VAOX(2,ifun,IP) = polynom*dery + x0(1)*x0(3)*xnorm1
                VAOX(3,ifun,IP) = polynom*derz + x0(1)*x0(2)*xnorm1
c
              ELSE IF(ityp.EQ.7) THEN
c -- f function (10 components)    xxx, xxy, xxz, xyy, xyz,
c                                  xzz, yyy, yyz, yzz, zzz
                ifun = ifun+1
c    xxx
                polynom = x0(1)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Three*xnorm*x0(1)**2
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz
                ifun = ifun+1
c    xxy
                polynom = x0(2)*x0(1)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Two*xnorm*x0(1)*x0(2)
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(1)**2
                VAOX(3,ifun,IP) = polynom*derz
                ifun = ifun+1
c    xxz
                polynom = x0(3)*x0(1)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + Two*xnorm*x0(1)*x0(3)
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(1)**2
                ifun = ifun+1
c    xyy
                polynom = x0(1)*x0(2)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(2)**2
                VAOX(2,ifun,IP) = polynom*dery + Two*xnorm*x0(1)*x0(2)
                VAOX(3,ifun,IP) = polynom*derz
                ifun = ifun+1
c    xyz
                polynom = x0(1)*x0(2)*x0(3)
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(2)*x0(3)
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(1)*x0(3)
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(1)*x0(2)
                ifun = ifun+1
c    xzz
                polynom = x0(1)*x0(3)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx + xnorm*x0(3)**2
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + Two*xnorm*x0(1)*x0(3)
                ifun = ifun+1
c    yyy
                polynom = x0(2)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + Three*xnorm*x0(2)**2
                VAOX(3,ifun,IP) = polynom*derz
                ifun = ifun+1
c    yyz
                polynom = x0(3)*x0(2)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + Two*xnorm*x0(2)*x0(3)
                VAOX(3,ifun,IP) = polynom*derz + xnorm*x0(2)**2
                ifun = ifun+1
c    yzz
                polynom = x0(2)*x0(3)**2
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery + xnorm*x0(3)**2
                VAOX(3,ifun,IP) = polynom*derz + Two*xnorm*x0(2)*x0(3)
                ifun = ifun+1
c    zzz
                polynom = x0(3)**3
                VAO(ifun,IP) = xnorm*polynom
                VAOX(1,ifun,IP) = polynom*derx
                VAOX(2,ifun,IP) = polynom*dery
                VAOX(3,ifun,IP) = polynom*derz + Three*xnorm*x0(3)**2
c
              ELSE
c  -- g functions -  currently not active
              call nerror(1,'AOGrad:',
     $             'g-functions and above not yet coded in DFT',0,0)
              ENDIF
c -- end of loop over general contractions   (igr loop)
 800          Continue
c -- end of loop over contracted shells      (ics loop)
 900          continue
c -- end of loop over grid points            (IP loop)
 950      CONTINUE
C
      return
      end
c ======================================================================
      SUBROUTINE AOInit(ncs,    basdat, inx,    bl)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  precomputes (2*exponent/pi)**0.75 for each contracted shell
C
C  ARGUMENTS
C
C  ncs     -  number of shells
C  basdat  -  basis set data for TEXAS
C  inx     -  more basis set data
C  bl      -  on exit contains precomputed factors
C
C
      DIMENSION basdat(13,*),inx(12,*),bl(*)
      parameter (PI=3.14159 26535 89793d0,Two=2.0d00)
C
C
      xcons=sqrt(Two*sqrt(Two/PI)/PI)
c
      DO 20 ics=1,ncs
        ia=inx(1,ics)+1
        ie=inx(5,ics)
c  contraction loop
        DO 10 ish=ia,ie
          expnt=basdat(1,ish)
          sqa=sqrt(expnt)
          bl(ish)=xcons*sqrt(expnt*sqa)
cc      write(6,*) ' shell: ',ish,' normalization factor: ',bl(ish)
 10     CONTINUE
 20   CONTINUE
C
      RETURN
      END
c ===================================================================
      SUBROUTINE AOSort(NP,     NBas,   thrsh,  VAO,    nbf,
     $                  INB,    VAOO,   VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  determines number of "non-zero" basis function values
C  at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglecting contribution
C  VAO     -  array of basis function values
C
C  on exit
C
C  nbf     -  indexing array to indices of "non-zero" AOs
C  INB     -  indexing array for non-zero entries to VAO
C  VAOO    -  stores values for all "non-zero" AOs
C  VM      -  maximum magnitude AO value per grid point
C
C  ** WARNING:   VAO and VAOO may share storage **
C
C
      DIMENSION VAO(NBas,NP),VAOO(NBas*NP),VM(NP)
      INTEGER   INB(NBas*NP),nbf(NP+1)
C
      nf = 0
      nbf(1) = 0
      DO 20 IP=1,NP
      ValM = 0.0d0
      DO 10 I=1,NBas
      If(Abs(VAO(I,IP)).GT.thrsh) Then
       nf = nf+1
       INB(nf) = I
       VAOO(nf) = VAO(I,IP)
       If(Abs(VAO(I,IP)).GT.ValM) ValM = Abs(VAO(I,IP))
      EndIf
 10   CONTINUE
      nbf(IP+1) = nf
      VM(IP) = ValM
 20   CONTINUE
C
C  on exit, number of non-zero AOs is nf which is stored
C  as the last entry in the nbf array
C
      RETURN
      END
c =====================================================================
      SUBROUTINE AOVal(NP,     ncs,    ncf,    XGRID,
     $                 xnuc,   bl,     basdat, inx,    thrsh,
     $                 ExpMIN, VAO)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates vector of basis function (AO) values and their
C  derivatives at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  ncs     -  number of shells
C  ncf     -  number of basis functions
C  XGRID   -  grid points
C  xnuc    -  nuclear coordinates
C  bl      -  array of precomputed normalization factors
C             (from routine <AOInit>)
C  basdat  -  basis set data for TEXAS
C  inx     -  more basis set data
C  thrsh   -  exponent threshold for neglecting contribution
C  ExpMIN  -  array of smallest exponents for all basis
C             functions per shell
C  VAO     -  on exit contains basis function values
C
C
      DIMENSION xnuc(3,*),bl(*),basdat(13,*),inx(12,*),ExpMIN(ncs),
     $          VAO(ncf,NP),XGRID(3,NP)
cc
      dimension x0(3)
      parameter (Zero=0.0d0,One=1.0d0,Two=2.0d0,Three=3.0d0,
     $           Four=4.0d0,Five=5.0d0,Half=0.5d0,Eight=8.0d0)
      parameter (PI=3.14159 26535 89793d0)
C
C
C  get some values
C
        sqtw=One/sqrt(12.0d0)
        sqr40=sqrt(40.0d0)
c
        DO 950 IP=1,NP
        xp = XGRID(1,IP)
        yp = XGRID(2,IP)
        zp = XGRID(3,IP)
        do 900 ics=1,ncs
c len=number of ang. mom. components, e.g. 3 for P, 4 for L, 5 for D;
c ia=beginning, ie=ending of the contraction
c iat=atom number, xnuc=coordinates of the center(nucleus)
c   ityp: 1=S, 2=P, 3=SP(L), 4=5D, 5=6D, 6=7F, 7=10F
          iat=inx(2,ics)                        ! atomic centre
          x0(1)=xp-xnuc(1,iat)
          x0(2)=yp-xnuc(2,iat)
          x0(3)=zp-xnuc(3,iat)
          r02=x0(1)**2+x0(2)**2+x0(3)**2
cc  ...............................................................
cc    global shell exponent test                ! jb march 97
cc  ...............................................................
          rrexp = r02*ExpMIN(ics)
c          write(6,*) ' ics is:',ics
c          write(6,*) ' rrexp is:',rrexp,' thrsh is:',thrsh
          If(rrexp.GT.thrsh) go to 900
cc  ...............................................................
          ityp=inx(12,ics)                       ! shell type
          ia=inx(1,ics)+1                        ! start of contraction
          ie=inx(5,ics)                          ! end of contraction
          ngr = inx(4,ics)                       ! no. of general contractions
          ifun=inx(11,ics)
c -- general contraction loop
            DO 800 igr=0,ngr
            xnorm=Zero
            xnorm1=Zero
c -- contraction loop
            do 700 ish=ia,ie
c  this the Gaussian exponent times the distance squared, a*r02
              aexp=basdat(1,ish)                 ! gaussian exponent
              ar02=aexp*r02
c  ..........................................................
c  do not calculate very small values
              If(ar02.GT.thrsh) go to 700
c  ..........................................................
c  bl(ish) contains the precomputed value (2*aexp/pi)**0.75
c   aexp is the exponent of the Gaussian
              xx=bl(ish)*exp(-ar02)
c basdat(2,ish) is the main contraction coefficient
c basdat(3,ish) is the P contraction coefficient
              contrs=basdat(2+igr,ish)
              xnorm=xnorm+xx*contrs
              if(ityp.eq.3) then                  ! segmented basis only
                contrp=basdat(3,ish)
                xnorm1=xnorm1+xx*contrp
              end if
 700        continue
c
c  the spherically symmetrical part of the AO at this point is in
c  xnorm. For L shells, the s component value is in xnorm, the P
c  component in xnorm1. These are the same for all angular momentum
c  components (e.g. X,Y,Z or XX,YY,ZZ,XY,XZ,YZ  etc.)
c loop over the function components in one contracted shell
c the radial part is identical
c
              IF(ityp.EQ.1) THEN
c -- s function
                ifun = ifun+1
                VAO(ifun,IP) = xnorm
              ELSE IF(ityp.EQ.2) THEN
c -- p function
                Do icomp=1,3
                ifun = ifun+1
                VAO(ifun,IP) = xnorm*x0(icomp)
                EndDo
              ELSE IF(ityp.EQ.3) THEN
c -- l function (sp)
                ifun = ifun+1
                VAO(ifun,IP) = xnorm
                Do icomp=1,3
                ifun = ifun+1
                VAO(ifun,IP) = xnorm1*x0(icomp)
                EndDo
              ELSE IF(ityp.EQ.4) THEN
c -- d function (5 components)    z**2,x**2-y**2,xy,xz,yz
                ifun = ifun+1
c  (2*z**2 -x**2-y**2)/SQRT(12)  = dz2
             VAO(ifun,IP) = xnorm*(Two*x0(3)**2-x0(1)**2-x0(2)**2)*sqtw
                ifun = ifun+1
c   (x**2-y**2)/2 = d(x2-y2)
                VAO(ifun,IP) = xnorm*(x0(1)**2-x0(2)**2)*Half
                ifun = ifun+1
c    xy
                VAO(ifun,IP) = xnorm*x0(1)*x0(2)
                ifun = ifun+1
c    xz
                VAO(ifun,IP) = xnorm*x0(1)*x0(3)
                ifun = ifun+1
c    yz
                VAO(ifun,IP) = xnorm*x0(2)*x0(3)
c
              ELSE IF(ityp.EQ.5) THEN
c -- d function (6 components)    xx,yy,zz,xy,xz,yz
                ifun = ifun+1
c    xx
                VAO(ifun,IP)   = xnorm*x0(1)**2
                ifun = ifun+1
c    yy
                VAO(ifun,IP)   = xnorm*x0(2)**2
                ifun = ifun+1
c    zz
                VAO(ifun,IP)   = xnorm*x0(3)**2
                ifun = ifun+1
c    xy
                VAO(ifun,IP)   = xnorm*x0(1)*x0(2)
                ifun = ifun+1
c    xz
                VAO(ifun,IP)   = xnorm*x0(1)*x0(3)
                ifun = ifun+1
c    yz
                VAO(ifun,IP)   = xnorm*x0(2)*x0(3)
c
              ELSE IF(ityp.EQ.6) THEN
c -- f function (7 components)    (5xxy-rry),(5xxz-rrz),(5yyx-rrx),
c                                 (5yyz-rrz),(5zzx-rrx),(5zzy-rry),
c                                    xyx*SQRT(40)
                ifun = ifun+1
c    5xxy
                VAO(ifun,IP)   = xnorm*(Five*x0(1)**2-r02)*x0(2)
                ifun = ifun+1
c    5xxz
                VAO(ifun,IP)   = xnorm*(Five*x0(1)**2-r02)*x0(3)
                ifun = ifun+1
c    5yyx
                VAO(ifun,IP)   = xnorm*(Five*x0(2)**2-r02)*x0(1)
                ifun = ifun+1
c    5yyz
                VAO(ifun,IP)   = xnorm*(Five*x0(2)**2-r02)*x0(3)
                ifun = ifun+1
c    5zzx
                VAO(ifun,IP)   = xnorm*(Five*x0(3)**2-r02)*x0(1)
                ifun = ifun+1
c    5zzy
                VAO(ifun,IP)   = xnorm*(Five*x0(3)**2-r02)*x0(2)
                ifun = ifun+1
c    xyz
                VAO(ifun,IP)   = xnorm*x0(1)*x0(2)*x0(3)*sqr40
c
              ELSE IF(ityp.EQ.7) THEN
c -- f function (10 components)    xxx, xxy, xxz, xyy, xyz,
c                                  xzz, yyy, yyz, yzz, zzz
                ifun = ifun+1
c    xxx
                VAO(ifun,IP) = xnorm*x0(1)**3
                ifun = ifun+1
c    xxy
                VAO(ifun,IP) = xnorm*x0(2)*x0(1)**2
                ifun = ifun+1
c    xxz
                VAO(ifun,IP) = xnorm*x0(3)*x0(1)**2
                ifun = ifun+1
c    xyy
                VAO(ifun,IP) = xnorm*x0(1)*x0(2)**2
                ifun = ifun+1
c    xyz
                VAO(ifun,IP) = xnorm*x0(1)*x0(2)*x0(3)
                ifun = ifun+1
c    xzz
                VAO(ifun,IP) = xnorm*x0(1)*x0(3)**2
                ifun = ifun+1
c    yyy
                VAO(ifun,IP) = xnorm*x0(2)**3
                ifun = ifun+1
c    yyz
                VAO(ifun,IP) = xnorm*x0(3)*x0(2)**2
                ifun = ifun+1
c    yzz
                VAO(ifun,IP) = xnorm*x0(2)*x0(3)**2
                ifun = ifun+1
c    zzz
                VAO(ifun,IP) = xnorm*x0(3)**3
c
              ELSE
c  -- g functions -  currently not active
              call nerror(1,'AOVal:',
     $             'g-functions and above not yet coded in DFT',0,0)
              ENDIF
c -- end of loop over general contractions   (igr loop)
 800          Continue
c -- end of loop over contracted shells      (ics loop)
 900          continue
c -- end of loop over grid points            (IP loop)
 950      CONTINUE
C
      return
      end
c ===================================================================
      SUBROUTINE AOXSort(NP,     NBas,   thrsh,  VAO,    VAOX,
     $                   nbf,    INB,    VAOO,   VAOOX,  VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  determines number of "non-zero" basis function or derivative
C  values at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglecting contribution
C  VAO     -  array of basis function values
C  VAOX    -  array of basis function derivatives
C
C  on exit
C
C  nbf     -  number of "non-zero" AOs and/or derivatives
C  INB     -  indexing array for non-zero entries to VAO/VAOX
C  VAOO    -  stores values for all "non-zero" AOs
C  VAOOX   -  stores values for all "non-zero" AO derivatives
C  VM      -  maximum magnitude AO value per grid point
C
C  ** WARNING:   VAO/VAOO and VAOX/VAOOX may share storage **
C
C
      DIMENSION VAO(NBas,NP),VAOX(3,NBas,NP),VAOO(NBas*NP),VAOOX(3,*),
     $          VM(NP)
      INTEGER   INB(NBas*NP),nbf(NP+1)
C
      nf = 0
      nbf(1) = 0
      DO 20 IP=1,NP
      ValM = 0.0d0
      DO 10 I=1,NBas
      If(Abs(VAO(I,IP)).GT.thrsh.OR.Abs(VAOX(1,I,IP)).GT.thrsh.OR.
     $   Abs(VAOX(2,I,IP)).GT.thrsh.OR.Abs(VAOX(3,I,IP)).GT.thrsh) Then
       nf = nf+1
       INB(nf) = I
       VAOO(nf) = VAO(I,IP)
       VAOOX(1,nf) = VAOX(1,I,IP)
       VAOOX(2,nf) = VAOX(2,I,IP)
       VAOOX(3,nf) = VAOX(3,I,IP)
       If(Abs(VAO(I,IP)).GT.ValM) ValM = Abs(VAO(I,IP))
      EndIf
 10   CONTINUE
      nbf(IP+1) = nf
      VM(IP) = ValM
 20   CONTINUE
C
      RETURN
      END
c ===================================================================
      SUBROUTINE AOXXSort(NP,     NBas,   thrsh,  VAO,    VAOX,
     $                    VAOXX,  nbf,    INB,    VAOO,   VAOOX,
     $                    VAOOXX, VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  determines number of "non-zero" basis function or derivative
C  values at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglecting contribution
C  VAO     -  array of basis function values
C  VAOX    -  array of basis function derivatives
C  VAOXX   -  array of basis function second derivatives
C
C  on exit
C
C  nbf     -  number of "non-zero" AOs and/or derivatives
C  INB     -  indexing array for non-zero entries to VAO/VAOX
C  VAOO    -  stores values for all "non-zero" AOs
C  VAOOX   -  stores values for all "non-zero" AO first derivatives
C  VAOOXX  -  stores values for all "non-zero" AO second derivatives
C  VM      -  maximum magnitude AO value per grid point
C
C  ** WARNING:   VAO/VAOO, VAOX/VAOOX and VAOXX/VAOOXX may share storage **
C
C
      DIMENSION VAO(NBas,NP),VAOX(3,NBas,NP),VAOXX(6,NBas,NP),
     $          VAOO(NBas*NP),VAOOX(3,*),VAOOXX(6,*),VM(NP)
      INTEGER   INB(NBas*NP),nbf(NP+1)
C
      nf = 0
      nbf(1) = 0
      DO 20 IP=1,NP
      ValM = 0.0d0
      DO 10 I=1,NBas
      ValX = MAX(Abs(VAO(I,IP)),Abs(VAOX(1,I,IP)),Abs(VAOX(2,I,IP)),
     $         Abs(VAOX(3,I,IP)),Abs(VAOXX(1,I,IP)),Abs(VAOXX(2,I,IP)),
     $        Abs(VAOXX(3,I,IP)),Abs(VAOXX(4,I,IP)),Abs(VAOXX(5,I,IP)),
     $        Abs(VAOXX(6,I,IP)))
      If(ValX.GT.thrsh) Then
cc      If(Abs(VAO(I,IP)).GT.thrsh.OR.Abs(VAOX(1,I,IP)).GT.thrsh.OR.
cc     $ Abs(VAOX(2,I,IP)).GT.thrsh.OR.Abs(VAOX(3,I,IP)).GT.thrsh.OR.
cc     $ Abs(VAOXX(1,I,IP)).GT.thrsh.OR.Abs(VAOXX(2,I,IP)).GT.thrsh.OR.
cc     $ Abs(VAOXX(3,I,IP)).GT.thrsh.OR.Abs(VAOXX(4,I,IP)).GT.thrsh.OR.
cc     $ Abs(VAOXX(5,I,IP)).GT.thrsh.OR.Abs(VAOXX(6,I,IP)).GT.thrsh) Then
       nf = nf+1
       INB(nf) = I
       VAOO(nf) = VAO(I,IP)
       VAOOX(1,nf) = VAOX(1,I,IP)
       VAOOX(2,nf) = VAOX(2,I,IP)
       VAOOX(3,nf) = VAOX(3,I,IP)
       VAOOXX(1,nf) = VAOXX(1,I,IP)
       VAOOXX(2,nf) = VAOXX(2,I,IP)
       VAOOXX(3,nf) = VAOXX(3,I,IP)
       VAOOXX(4,nf) = VAOXX(4,I,IP)
       VAOOXX(5,nf) = VAOXX(5,I,IP)
       VAOOXX(6,nf) = VAOXX(6,I,IP)
cc       If(Abs(VAO(I,IP)).GT.ValM) ValM = Abs(VAO(I,IP))
       If(ValX.GT.ValM) ValM = ValX
      EndIf
 10   CONTINUE
      nbf(IP+1) = nf
      VM(IP) = ValM
 20   CONTINUE
C
      RETURN
      END
c ===================================================================
      SUBROUTINE AOXXXSort(NP,     NBas,   thrsh,  VAO,    VAOX,
     $                    VAOXX,   VAOXXX, nbf,    INB,    VAOO,
     $                    VAOOX,   VAOOXX, VAOOXXX,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  determines number of "non-zero" basis function or derivative
C  values at the current grid point
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch
C  NBas    -  total number of basis functions
C  thrsh   -  threshold for neglecting contribution
C  VAO     -  array of basis function values
C  VAOX    -  array of basis function derivatives
C  VAOXX   -  array of basis function second derivatives
C  VAOXXX  -  array of basis function third derivatives
C
C  on exit
C
C  nbf     -  number of "non-zero" AOs and/or derivatives
C  INB     -  indexing array for non-zero entries to VAO/VAOX
C  VAOO    -  stores values for all "non-zero" AOs
C  VAOOX   -  stores values for all "non-zero" AO first derivatives
C  VAOOXX  -  stores values for all "non-zero" AO second derivatives
C  VAOOXXX -  stores values for all "non-zero" AO third derivatives
C  VM      -  maximum magnitude AO value per grid point
C
C  ** WARNING:   VAO/VAOO, VAOX/VAOOX, VAOXX/VAOOXX,
C                and VAOXXX/VAOOXXX may share storage **
C
C
      DIMENSION VAO(NBas,NP),VAOX(3,NBas,NP),VAOXX(6,NBas,NP),
     $          VAOO(NBas*NP),VAOOX(3,*),VAOOXX(6,*),VM(NP)
      DIMENSION VAOXXX(10,NBas,NP),VAOOXXX(10,*)
      INTEGER   INB(NBas*NP),nbf(NP+1)
C
      nf = 0
      nbf(1) = 0
      DO 20 IP=1,NP
      ValM = 0.0d0
      DO 10 I=1,NBas
      ValX = MAX(Abs(VAO(I,IP)),Abs(VAOX(1,I,IP)),Abs(VAOX(2,I,IP)),
     $     Abs(VAOX(3,I,IP)),Abs(VAOXX(1,I,IP)),Abs(VAOXX(2,I,IP)),
     $     Abs(VAOXX(3,I,IP)),Abs(VAOXX(4,I,IP)),Abs(VAOXX(5,I,IP)),
     $     Abs(VAOXX(6,I,IP)),Abs(VAOXXX(1,I,IP)),Abs(VAOXXX(2,I,IP)),
     $     Abs(VAOXXX(3,I,IP)),Abs(VAOXXX(4,I,IP)),Abs(VAOXXX(5,I,IP)),
     $     Abs(VAOXXX(6,I,IP)),Abs(VAOXXX(7,I,IP)),Abs(VAOXXX(8,I,IP)),
     $     Abs(VAOXXX(9,I,IP)),Abs(VAOXXX(10,I,IP)))
      If(ValX.GT.thrsh) Then
       nf = nf+1
       INB(nf) = I
       VAOO(nf) = VAO(I,IP)
       VAOOX(1,nf) = VAOX(1,I,IP)
       VAOOX(2,nf) = VAOX(2,I,IP)
       VAOOX(3,nf) = VAOX(3,I,IP)
       VAOOXX(1,nf) = VAOXX(1,I,IP)
       VAOOXX(2,nf) = VAOXX(2,I,IP)
       VAOOXX(3,nf) = VAOXX(3,I,IP)
       VAOOXX(4,nf) = VAOXX(4,I,IP)
       VAOOXX(5,nf) = VAOXX(5,I,IP)
       VAOOXX(6,nf) = VAOXX(6,I,IP)
       VAOOXXX(1,nf)  = VAOXXX(1,I,IP)
       VAOOXXX(2,nf)  = VAOXXX(2,I,IP)
       VAOOXXX(3,nf)  = VAOXXX(3,I,IP)
       VAOOXXX(4,nf)  = VAOXXX(4,I,IP)
       VAOOXXX(5,nf)  = VAOXXX(5,I,IP)
       VAOOXXX(6,nf)  = VAOXXX(6,I,IP)
       VAOOXXX(7,nf)  = VAOXXX(7,I,IP)
       VAOOXXX(8,nf)  = VAOXXX(8,I,IP)
       VAOOXXX(9,nf)  = VAOXXX(9,I,IP)
       VAOOXXX(10,nf) = VAOXXX(10,I,IP)
cc       If(Abs(VAO(I,IP)).GT.ValM) ValM = Abs(VAO(I,IP))
       If(ValX.GT.ValM) ValM = ValX
      EndIf
 10   CONTINUE
      nbf(IP+1) = nf
      VM(IP) = ValM
 20   CONTINUE
C
      RETURN
      END
c =====================================================================
      SUBROUTINE RdAOVal(NP,nbf,VAO,INB,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads AO values over grid for current atom from disk
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch + 1
C  nbf     -  number of "non-zero" basis functions per grid point
C  VAO     -  basis function values at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  VM      -  maximum magnitude AO value per grid point
C
C
      DIMENSION VAO(*),VM(NP-1)
      INTEGER nbf(NP),INB(*)
c
      READ(42,ERR=95) NMax,nbf,VM
      If(NMax.GT.0) CALL ReadBinaryAO(42,NMax,INB,VAO)
      RETURN
c
 95   CONTINUE
      Call nerror(13,'File IO routine <RdAOVal>',
     $  'Error Reading AO grid file in DFT Energy Calculation',0,0)
C
      END
c =====================================================================
      subroutine ReadBinaryAO(IUnit,NMax,INB,VAO)
      real*8 INB(NMax),VAO(NMax)
c
      READ(IUnit) INB,VAO
c
      RETURN
      END
c =====================================================================
      SUBROUTINE RdAOGrad(NP,nbf,VAO,VAOX,INB,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads AO values over grid for current atom from disk
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch + 1
C  nbf     -  number of "non-zero" basis functions per grid point
C  VAO     -  basis function values at grid point
C  VAOX    -  basis function derivatives at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  VM      -  maximum magnitude AO value per grid point
C
C
      DIMENSION VAO(*),VAOX(3,*),VM(NP-1)
      INTEGER nbf(NP),INB(*)
c
      READ(42,ERR=95) NMax,nbf,VM
      If(NMax.GT.0) CALL ReadBinaryGD(42,NMax,INB,VAO,VAOX)
      RETURN
c
 95   CONTINUE
      Call nerror(13,'File IO routine <RdAOGrad>',
     $  'Error Reading AO grid file in DFT Energy Calculation',0,0)
C
      END
c =====================================================================
      subroutine ReadBinaryGD(IUnit,NMax,INB,VAO,VAOX)
      real*8 INB(NMax),VAO(NMax),VAOX(3,NMax)
c
      READ(IUnit) INB,VAO,VAOX
c
      RETURN
      END
c =====================================================================
      SUBROUTINE WrAOVal(NP,nbf,VAO,INB,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Writes AO values over grid for current atom to disk
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch + 1
C  nbf     -  number of "non-zero" basis functions per grid point
C  VAO     -  basis function values at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  VM      -  maximum magnitude AO value per grid point
C
C
      DIMENSION VAO(*),VM(NP-1)
      INTEGER nbf(NP),INB(*)
C
      nf = nbf(NP)
      WRITE(42,ERR=95) nf,nbf,VM
      If(nf.GT.0) CALL WriteBinaryAO(42,nf,INB,VAO)
      RETURN
c
 95   CONTINUE
      Call nerror(13,'File IO routine <WrAOVal>',
     $     'Unable to write to grid AO file  unit =',42,0)
C
      END
c =====================================================================
      subroutine WriteBinaryAO(IUnit,NMax,INB,VAO)
      real*8 INB(NMax),VAO(NMax)
c
      WRITE(IUnit) INB,VAO
c
      RETURN
      END
c =====================================================================
      SUBROUTINE WrAOGrad(NP,nbf,VAO,VAOX,INB,VM)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Writes AO values over grid for current atom to disk
C
C  ARGUMENTS
C
C  NP      -  number of grid points this batch + 1
C  nbf     -  number of "non-zero" basis functions per grid point
C  VAO     -  basis function values at grid point
C  VAOX    -  basis function derivatives at grid point
C  INB     -  indexing array for non-zero entries to VAO
C  VM      -  maximum magnitude AO value per grid point
C
C
      DIMENSION VAO(*),VAOX(3,*),VM(NP-1)
      INTEGER nbf(NP),INB(*)
c
      nf = nbf(NP)
      WRITE(42,ERR=95) nf,nbf,VM
      If(nf.GT.0) CALL WriteBinaryGD(42,nf,INB,VAO,VAOX)
      RETURN
c
 95   CONTINUE
      Call nerror(13,'File IO routine <WrAOVal>',
     $     'Unable to write to grid AO file  unit =',42,0)
C
      END
c =====================================================================
      subroutine WriteBinaryGD(IUnit,NMax,INB,VAO,VAOX)
      real*8 INB(NMax),VAO(NMax),VAOX(3,NMax)
c
      WRITE(IUnit) INB,VAO,VAOX
c
      RETURN
      END

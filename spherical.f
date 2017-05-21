      SUBROUTINE spherical_harmonics(x,y,z,LMax,Ylm)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Calculates values of the spherical harmonics, Y(L.M)
C  for a given grid point (x,y,z) and maximum L value
C
C  ARGUMENTS
C
C  x       -  x-coordinate of grid point
C  y       -  y-coordinate of grid point
C  z       -  z-coordinate of grid point
C  LMax    -  maximum L value
C
C  on exit
C
C  Ylm     -  spherical harmonics
C
C
      DIMENSION Plm(0:LMax,0:LMax),Ylm(0:LMax,-LMax:LMax)
      PARAMETER (One=1.0d0,Two=2.0d0,Four=4.0d0)
C
C
      pi = Four*ATAN(One)
      sqtwopi = One/SQRT(Two*pi)
      sqrt2 = One/SQRT(Two)
c
c -- convert Cartesian to Spherical coordinates
      call xyz_ang(x,y,z,th,phi)
c
c -- compute associated Legendre polynomials
      call AssocLegndr(LMax,Plm,th)
c
      do l=0,LMax
        fact=sqtwopi*sqrt(dble(2*l+1))
        do m=0,l
          Plm(l,m)=fact*Plm(l,m)
          if(m.eq.0) Plm(l,m) = Plm(l,m)*sqrt2
          fact=fact/sqrt(dble((l-m)*(l+m+1)))
        enddo
      enddo
c
      do l=0,Lmax
        Ylm(l,0)=Plm(l,0)
        do m=1,l
          Ylm(l,m)=Plm(l,m)*cos(m*phi)
          Ylm(l,-m)=Plm(l,m)*sin(m*phi)
        enddo
      enddo
C
      return
      end
c ============================================================
c
      subroutine AssocLegndr(LMax,Plm,z)
      implicit real*8 (a-h,o-z)
c  this routine calculates all associated Legendre polynomials
c  up to LMax, stored in the quadratic array Plm(l,m)
      dimension Plm(0:LMax,0:LMax)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      Plm(0,0)=one
      Plm(1,0)=z
      sqz=sqrt(one-z**2)
      pmm=one
      fact=one
      do i=1,LMax
        pmm=-pmm*fact*sqz
        fact=fact+two
        Plm(i,i)=pmm
        if(i.lt.LMax) Plm(i+1,i)=pmm*z*fact
      end do
      fact=one
      do l=2,LMax
        fact=fact+two
        fact1=l
        do m=0,l-2
          Plm(l,m)=(z*fact*Plm(l-1,m)-dble(l+m-1)*Plm(l-2,m))/fact1
          fact1=fact1-one
        end do
      end do
      return
      end
c ============================================================
c
      subroutine xyz_ang(x,y,z,th,phi)
      implicit real*8(a-h,o-z)
c  This routine converts cartesian coordinates to
c  spherical coordinates
      parameter(half=0.5D0,zero=0.0d0,thrshd=1.0D-10)
      pi=4.0d0*atan(1.0d0)
      r2=sqrt(x*x+y*y+z*z)
      if ( r2.gt.thrshd ) then
        th=z/r2
      else if ( abs(z).gt.thrshd ) then
        th=z/abs(z)
      else
        th=1.0d0
      endif
      if ( abs(x).gt.thrshd ) then
        phi=atan(y/x)
        if (x.lt.zero) phi=phi+pi
      else if (y.gt.zero) then
        phi=half*pi
      else
        phi=-half*pi
      endif
      return
      end
c ============================================================
c
      SUBROUTINE RadGRID(N,r1,r2,x,w)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Radial integration by Gauss-Legendre quadrature
C
C  ARGUMENTS
C
C  N       -  number of radial integration points
C  r1      -  lower integration limit
C  r2      -  upper integration limit
C
C  on exit
C
C  x       -  abscissas of the Gauss-Legendre quadrature
C  w       -  radial quadrature weights
C
      Dimension x(N),w(N)
C
      PARAMETER (PI=3.14159 26535 89793d0)
      PARAMETER (FourPI=4.0d0*PI)
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0)
      PARAMETER (EPS=3.0d-14)
C
      M = (N+1)/2
      XM = 0.5*(r2-r1)
      XL = 0.5*(r2+r1)
c
      DO I=1,M
      z = cos(PI*(I-0.25d0)/(N+0.5d0))
  10  continue
      p1 = One
      p2 = Zero
      do J = 1, N
        p3 = p2
        p2 = p1
        p1 = ((Two*J-One)*z*p2-(J-One)*p3)/J
      enddo
      pp = N*(z*p1-p2)/(z*z-One)
      z1 = z
      z  = z1-p1/pp
      If(abs(z-z1).gt.EPS) goto 10
      x(I) = XM-XL*z
      x(N+1-I) = XM+XL*z
      w(I) = Two*XL*FourPI*x(I)*x(I)/((One-z*z)*pp*pp)
      w(N+1-I) = w(I)*FourPI*x(N+1-I)**2
      EndDO
C
      RETURN
      END

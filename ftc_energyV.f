      Subroutine make_gridfunctions_f(
     &    itype,         ccs,           expvalx,       xexpvalx,
     &    xxexpvalx,     Px,            Py,            Pz,
     &    fockcontr,     iymin,         iymax,         izmin,
     &    izmax,         Lxmin,         Lymin,         Lzmin,
     &    npwy,          npwz,          griddens,      yexpstore,
     &    npwyexp,       zexpstore,     ro6,           ro4,
     &    npwzexp)
c
c THIS SUBROUTINE CALCULATES THE 2D COULOMB MATRIX CONTRIBUTION OF THE
C g(i)*g(j) SMOOTH PAIR USING THE PRECALCULATED SMOOTH+MIX COULOMB
C POTENTIAL.
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer itype,iymin,iymax,izmin,izmax,npwy,npwz
      integer npwyexp,npwzexp
      real*8 ro6(npwz,npwy) !g(i) basis function on the grid , input
      real*8 ro4(npwz,npwy) !smooth potential, input
      real*8 Lxmin,Lymin,Lzmin,ccs
      real*8 expvalx,xexpvalx,xxexpvalx,Px,Py,Pz,griddens
      real*8 yexpstore(npwyexp) !Pre calc exp functions
      real*8 zexpstore(npwzexp)
      real*8 fockcontr(10)  !fock mx contribution, output
c
      data pi/3.1415926535897932384626433d0/
c
      dx=1.0d0/griddens
      xinit=Lxmin - Px - dx
      yinit=Lymin - Py - dx
      y0=yinit+float(iymin)*dx
      y=y0
      zinit=Lzmin - Pz - dx
      z0=zinit+float(izmin)*dx
      z=z0
      ix=1
      x=xinit+dx
      indexy=1
c
      focksumm1=0.0d0
      focksumm2=0.0d0
      focksumm3=0.0d0
      focksumm4=0.0d0
      focksumm5=0.0d0
      focksumm6=0.0d0
      focksumm7=0.0d0
      focksumm8=0.0d0
      focksumm9=0.0d0
      focksumm10=0.0d0
c
      if (itype .eq. 1) then
c
c  s type function
c
        do iy=iymin,iymax
c
          valxy=expvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            pot=ro4(iz,iy)
            fsm=valxy*zexpstore(indexz)
            indexz=indexz+1
            focksumm1=focksumm1+pot*fc*fsm
c
          end do
        end do
        fockcontr(1)=focksumm1
c
      else if (itype .eq. 2) then
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxypx=xexpvalx*expy
          valxypy=expvalx*y*yexpstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fctpot=ro6(iz,iy)*ro4(iz,iy)*zexpstore(indexz)
c
            focksumm1=focksumm1+fctpot*valxypx
c now px is done
            focksumm2=focksumm2+fctpot*valxypy
c now py is done
            focksumm3=focksumm3+fctpot*valxypz*z
c now pz is done
c
            indexz=indexz+1
            z=z+dx
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
c
      else if (itype .eq. 3) then
c
        valxs=ccs*expvalx
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxys=valxs*expy
          valxypx=xexpvalx*expy
          valxypy=expvalx*y*yexpstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fctpot=ro6(iz,iy)*ro4(iz,iy)*zexpstore(indexz)
c
            focksumm1=focksumm1+fctpot*valxys
c now s is done
            focksumm2=focksumm2+fctpot*valxypx
c now px is done
            focksumm3=focksumm3+fctpot*valxypy
c now py is done
            focksumm4=focksumm4+fctpot*z*valxypz
c now pz is done
c
            indexz=indexz+1
            z=z+dx
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
        fockcontr(4)=focksumm4
c
      else if (itype .eq. 4) then
c
c D functions
c
        valx1=expvalx/dsqrt(12.0d0)
        valx2=expvalx/2.0d0
        valx34=expvalx
        r2x=x**2
        do iy=iymin,iymax
c
          y2=y**2
          r2xy=r2x+y2 !x^2+y^2
          r2xy2=r2xy-2.0d0*y2 !x^2-y^2
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*r2xy2*expy
          valxy5=valx34*expy
          valxy4=valxy5*x
          valxy3=valxy4*y
          valxy5=valxy5*y
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fctpot=ro6(iz,iy)*ro4(iz,iy)*zexpstore(indexz)
            a = 2.0d0*z**2 - r2xy
c
            focksumm1=focksumm1+fctpot*a*valxy1
c now D1 is done
            focksumm2=focksumm2+fctpot*valxy2
c now D2 is done
            focksumm3=focksumm3+fctpot*valxy3
c now D3 is done
            focksumm4=focksumm4+fctpot*z*valxy4
c now D4 is done
            focksumm5=focksumm5+fctpot*z*valxy5
c now D5 is done
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
        fockcontr(4)=focksumm4
        fockcontr(5)=focksumm5
c
      else if (itype .eq. 5) then
c
c D6 functions
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxy1=xxexpvalx*expy
          valxy2=expvalx*expy*y**2
          valxy3=expvalx*expy
          valxy5=xexpvalx*expy
          valxy4=valxy5*y
          valxy6=y*valxy3
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            pot=ro4(iz,iy)
            fctpot=fc*pot*zexpstore(indexz)
            zfctpot=z*fctpot
c
            focksumm1=focksumm1+fctpot*valxy1
c now D61 is done
            focksumm2=focksumm2+fctpot*valxy2
c now D62 is done
            focksumm3=focksumm3+zfctpot*z*valxy3
c now D63 is done
            focksumm4=focksumm4+fctpot*valxy4
c now D64 is done
            focksumm5=focksumm5+zfctpot*valxy5
c now D65 is done
            focksumm6=focksumm6+zfctpot*valxy6
c now D66 is done
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
        fockcontr(4)=focksumm4
        fockcontr(5)=focksumm5
        fockcontr(6)=focksumm6
c
      else if (itype .eq. 6) then
c
c F functions
c
        x=xinit+float(ix)*dx
        fourtx=4.0d0*x
        ax46=x**2
        ax1=4.0d0*ax46
        ax35=x*ax46
        ax7=dsqrt(40.0d0)
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          fourty=4.0d0*y
          y2=y**2
          y3=y2*y
          axy1=ax1*y-y3
          axy2=ax1-y2
          axy3=4.0d0*x*y2-ax35
          axy4=4.0d0*y2-ax46
          axy5=x*y2+ax35
          axy6=ax46*y+y3
          axy7=ax7*x*y
          valxy=expvalx*expy
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            vals=valxy*zexpstore(indexz)*ro6(iz,iy)*ro4(iz,iy)
c
            z2=z**2
            z3=z2*z
            a1=axy1-y*z2
            a2=axy2*z-z3
            a3=axy3-x*z2
            a4=axy4*z-z3
            a5=fourtx*z2-axy5
            a6=fourty*z2-axy6 !must be 4*y*z^2-x^2*y-y^3
            a7=axy7*z
c
            focksumm1=focksumm1+vals*a1
c now F1 is done
            focksumm2=focksumm2+vals*a2
c now F2 is done
            focksumm3=focksumm3+vals*a3
c now F3 is done
            focksumm4=focksumm4+vals*a4
c now F4 is done
            focksumm5=focksumm5+vals*a5
c now F5 is done
            focksumm6=focksumm6+vals*a6
c now F6 is done
            focksumm7=focksumm7+vals*a7
c now F7 is done
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
        fockcontr(4)=focksumm4
        fockcontr(5)=focksumm5
        fockcontr(6)=focksumm6
        fockcontr(7)=focksumm7
c
      else if (itype .eq. 7) then
c
c F10 functions
c
        x=xinit+float(ix)*dx
        x2=x**2
        x3=x*x2
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          y2=y**2
          y3=y2*y
          xy=x*y
          x2y=x*xy
          xy2=y*xy
          valxy=expvalx*expy
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            vals=valxy*zexpstore(indexz)*ro6(iz,iy)*ro4(iz,iy)
c
            z2=z**2
            a1=x3
            a2=x2y
            a3=x2*z
            a4=xy2
            a5=xy*z
            a6=x*z2
            a7=y3
            a8=y2*z
            a9=y*z2
            a10=z*z2
c
            focksumm1=focksumm1+vals*a1
c now F1 is done
            focksumm2=focksumm2+vals*a2
c now F2 is done
            focksumm3=focksumm3+vals*a3
c now F3 is done
            focksumm4=focksumm4+vals*a4
c now F4 is done
            focksumm5=focksumm5+vals*a5
c now F5 is done
            focksumm6=focksumm6+vals*a6
c now F6 is done
            focksumm7=focksumm7+vals*a7
c now F7 is done
            focksumm8=focksumm8+vals*a8
c now F8 is done
            focksumm9=focksumm9+vals*a9
c now F9 is done
            focksumm10=focksumm10+vals*a10
c now F10 is done
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
        fockcontr(4)=focksumm4
        fockcontr(5)=focksumm5
        fockcontr(6)=focksumm6
        fockcontr(7)=focksumm7
        fockcontr(8)=focksumm8
        fockcontr(9)=focksumm9
        fockcontr(10)=focksumm10
c
c
      else
        write(6,*) ' FTC: Can handle only S, P, D and F Shells'
        stop
      end if
c
c
      return
      end
c
c
c***********************************************************************
c
c
      Subroutine make_gridfunctions_d(
     &    itype,         ccs,           expvalx,       xexpvalx,
     &    xxexpvalx,     Px,            Py,            Pz,
     &    iymin,         iymax,         izmin,         izmax,
     &    x0,            Lymin,         Lzmin,         npwy,
     &    npwz,          griddens,      yexpstore,     npwyexp,
     &    zexpstore,     ro5,           ro6,           densmx,
     &    npwzexp)
c
c THIS SUBROUTINE BUILDS UP THE SMOOTH DENSITY CONTRIBUTION OF THE
c g(i)*summ{dij*g(j)} PRODUCT IN 2D AT A GIVEN ix GRIDPOINT.
c THE SUMMATION ABOVE IS OVER THE PRIMITIVE BASIS FUNCTIONS WITHIN
c A GIVEN SHELL.
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer itype !type of the j shell
      integer iymin,iymax,izmin,izmax,npwy,npwz,npwyexp !grid dims
                                                        !and ranges
      real*8 Px,Py,Pz !Center of the j shell
      real*8 x0,Lymin,Lzmin ! x value and left corner of the box
      real*8 expvalx,xexpvalx,xxexpvalx !Precalculated expx,
                                        !x*expx and x^2*expx
      real*8 ccs !spec factor for the l shells
      real*8 ro6(npwz,npwy)     ! grid function of g(i) basis function
      real*8 yexpstore(npwyexp) ! Precalculated 1d exp functions
      real*8 zexpstore(npwzexp)
      real*8 densmx(10)         ! the given denity mx elements
      real*8 ro5(npwz,npwy)     ! smooth density contribution
c
      data pi/3.1415926535897932384626433d0/
c
c
      dx=1.0d0/griddens
      x=x0 - Px
      yinit=Lymin - Py - dx
      y0=yinit+float(iymin)*dx
      y=y0
      zinit=Lzmin - Pz - dx
      z0=zinit+float(izmin)*dx
      z=z0
      indexy=1
c
      if (itype .eq. 1) then
c
c  s type function
c
        d1=densmx(1)
        do iy=iymin,iymax
c
          valxy=expvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            fsm=valxy*zexpstore(indexz)
            indexz=indexz+1
            ro5(iz,iy)=ro5(iz,iy)+d1*fc*fsm
c
          end do
        end do
c
      else if (itype .eq. 2) then
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxypx=xexpvalx*expy
          valxypy=expvalx*y*yexpstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            expz=zexpstore(indexz)
c
            fsm=valxypx*expz
            denssumm1=d1*fc*fsm
c now px is done
            fsm=valxypy*expz
            denssumm2=d2*fc*fsm
c now py is done
            fsm=z*valxypz*expz
            denssumm3=d3*fc*fsm
c now pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if (itype .eq. 3) then
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
        valxs=ccs*expvalx
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxys=valxs*expy
          valxypx=xexpvalx*expy
          valxypy=expvalx*y*yexpstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            expz=zexpstore(indexz)
c
            fsm=valxys*expz
            denssumm1=d1*fc*fsm
c now s is done
            fsm=valxypx*expz
            denssumm2=d2*fc*fsm
c now px is done
            fsm=valxypy*expz
            denssumm3=d3*fc*fsm
c now py is done
            fsm=z*valxypz*expz
            denssumm4=d4*fc*fsm
c now pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if (itype .eq. 4) then
c
c D functions
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
        d5=densmx(5)
c
        valx1=expvalx/dsqrt(12.0d0)
        valx2=expvalx/2.0d0
        valx34=expvalx
        r2x=x**2
        do iy=iymin,iymax
c
          y2=y**2
          r2xy=r2x+y2 !x^2+y^2
          r2xy2=r2xy-2.0d0*y2 !x^2-y^2
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*r2xy2*expy
          valxy5=valx34*expy
          valxy4=valxy5*x
          valxy3=valxy4*y
          valxy5=valxy5*y
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            a = 2.0d0*z**2 - r2xy
c
c
            fsm1=zexpstore(indexz)
            fsm=a*valxy1*fsm1
            denssumm1=d1*fc*fsm
c now D1 is done
            fsm=valxy2*fsm1
            denssumm2=d2*fc*fsm
c now D2 is done
            fsm=valxy3*fsm1
            denssumm3=d3*fc*fsm
c now D3 is done
            fsm=z*valxy4*fsm1
            denssumm4=d4*fc*fsm
c now D4 is done
            fsm=z*valxy5*fsm1
            denssumm5=d5*fc*fsm
c now D5 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if (itype .eq. 5) then
c
c D6 functions
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
        d5=densmx(5)
        d6=densmx(6)
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxy1=xxexpvalx*expy
          valxy2=expvalx*expy*y**2
          valxy3=expvalx*expy
          valxy5=xexpvalx*expy
          valxy4=valxy5*y
          valxy6=y*valxy3
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            fsm1=zexpstore(indexz)
            ztfsm1=z*fsm1
c
            fsm=valxy1*fsm1
            denssumm1=d1*fc*fsm
c now D61 is done
            fsm=valxy2*fsm1
            denssumm2=d2*fc*fsm
c now D62 is done
            fsm=valxy3*ztfsm1*z
            denssumm3=d3*fc*fsm
c now D63 is done
            fsm=valxy4*fsm1
            denssumm4=d4*fc*fsm
c now D64 is done
            fsm=valxy5*ztfsm1
            denssumm5=d5*fc*fsm
c now D65 is done
            fsm=valxy6*ztfsm1
            denssumm6=d6*fc*fsm
c now D66 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5+denssumm6
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if (itype .eq. 6) then
c
c F functions
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
        d5=densmx(5)
        d6=densmx(6)
        d7=densmx(7)
c
        fourtx=4.0d0*x
        ax46=x**2
        ax1=4.0d0*ax46
        ax35=x*ax46
        ax7=dsqrt(40.0d0)
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          fourty=4.0d0*y
          y2=y**2
          y3=y2*y
          axy1=ax1*y-y3
          axy2=ax1-y2
          axy3=4.0d0*x*y2-ax35
          axy4=4.0d0*y2-ax46
          axy5=x*y2+ax35
          axy6=ax46*y+y3
          axy7=ax7*x*y
          valxy=expvalx*expy
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            z2=z**2
            z3=z2*z
            a1=axy1-y*z2
            a2=axy2*z-z3
            a3=axy3-x*z2
            a4=axy4*z-z3
            a5=fourtx*z2-axy5
            a6=fourty*z2-axy6 !must be 4*y*z^2-x^2*y-y^3
            a7=axy7*z
c
            fsm1=zexpstore(indexz)
            vals=valxy*fsm1
c
            fsm=a1*vals
            denssumm1=d1*fc*fsm
c now F1 is done
            fsm=a2*vals
            denssumm2=d2*fc*fsm
c now F2 is done
            fsm=a3*vals
            denssumm3=d3*fc*fsm
c now F3 is done
            fsm=a4*vals
            denssumm4=d4*fc*fsm
c now F4 is done
            fsm=a5*vals
            denssumm5=d5*fc*fsm
c now F5 is done
            fsm=a6*vals
            denssumm6=d6*fc*fsm
c now F6 is done
            fsm=a7*vals
            denssumm7=d7*fc*fsm
c now F7 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                  denssumm4+denssumm5+denssumm6+denssumm7
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if (itype .eq. 7) then
c
c F10 functions
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
        d5=densmx(5)
        d6=densmx(6)
        d7=densmx(7)
        d8=densmx(8)
        d9=densmx(9)
        d10=densmx(10)
c
        x2=x**2
        x3=x*x2
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          y2=y**2
          y3=y2*y
          xy=x*y
          x2y=xy*x
          xy2=xy*y
          valxy=expvalx*expy
c
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz,iy)
            z2=z**2
            a1=x3
            a2=x2y
            a3=x2*z
            a4=xy2
            a5=xy*z
            a6=x*z2
            a7=y3
            a8=y2*z
            a9=y*z2
            a10=z*z2
c
            fsm1=zexpstore(indexz)
            vals=valxy*fsm1
c
            fsm=a1*vals
            denssumm1=d1*fc*fsm
c now F1 is done
            fsm=a2*vals
            denssumm2=d2*fc*fsm
c now F2 is done
            fsm=a3*vals
            denssumm3=d3*fc*fsm
c now F3 is done
            fsm=a4*vals
            denssumm4=d4*fc*fsm
c now F4 is done
            fsm=a5*vals
            denssumm5=d5*fc*fsm
c now F5 is done
            fsm=a6*vals
            denssumm6=d6*fc*fsm
c now F6 is done
            fsm=a7*vals
            denssumm7=d7*fc*fsm
c now F7 is done
            fsm=a8*vals
            denssumm8=d8*fc*fsm
c now F8 is done
            fsm=a9*vals
            denssumm9=d9*fc*fsm
c now F9 is done
            fsm=a10*vals
            denssumm10=d10*fc*fsm
c now F10 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3
     &                           +denssumm4+denssumm5+denssumm6
     &                           +denssumm7+denssumm8+denssumm9
     &                           +denssumm10
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else
        write(6,*) ' FTC: Can handle only S, P, D and F Shells'
        stop
      end if
c
c
      return
      end

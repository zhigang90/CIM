      Subroutine make_coulomb_smoothc(Lx,Ly,Lz,npwx,npwy,
     &npwz,D,ro1)
c
      IMPLICIT REAL*8(A-H,O-Z)
      complex*16 ro1(npwz,npwy,npwx)
      real*8 Lx,Ly,Lz
      data pi/3.1415926535897932384626433d0/
      fourpi=4.0d0*pi
c
      cx=2.0d0*pi/Lx
      cx=cx**2
      cy=2.0d0*pi/Ly
      cy=cy**2
      cz=2.0d0*pi/Lz
      cz=cz**2
c
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      do jkx=1,npwx
        kx=(jkx-1)-npwx*int(jkx/nx)
        valksqx=cx*dfloat(kx)**2
        do jky=1,npwy
          ky=(jky-1)-npwy*int(jky/ny)
          valksqxy=valksqx+cy*dfloat(ky)**2
          do jkz=1,npwz
            kz=(jkz-1)-npwz*int(jkz/nz)
c
            valksq=valksqxy + cz*dfloat(kz)**2
            valkabs=sqrt(valksq)
c
            if((kx .eq. 0) .and. (ky .eq. 0) .and. (kz .eq. 0)) then
c
              ro1(jkz,jky,jkx)= 2.0d0 * pi * D**2 *
     &        ro1(jkz,jky,jkx)
c
            else
c
              ro1(jkz,jky,jkx)=(fourpi/valksq)*(1.0d0-cos(valkabs*D))*
     &        ro1(jkz,jky,jkx)
c
            end if
c
          end do
        end do
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine make_coulomb_smoothc_2D(
     &    cx,         cy,         cz,         ix,         npwx,
     &    npwy,       npwz,       nx,         ny,         nz,
     &    D,          ro1)
c
      implicit real*8(A-H,O-Z)
      real*8 ro1(2,npwz,npwy)
c
      kx=(ix-1)-npwx*int(ix/nx)
      valksqx=cx*dfloat(kx)**2
      do jky=1,npwy
        ky=(jky-1)-npwy*int(jky/ny)
        valksqxy=valksqx+cy*dfloat(ky)**2
        do jkz=1,npwz
          kz=(jkz-1)-npwz*int(jkz/nz)
c
          valksq=valksqxy + cz*dfloat(kz)**2
          valkabs=sqrt(valksq)
c                       4*%pi
          cm= (1.2566370614359172954d1/valksq) * (1.0d0-cos(valkabs*D))
          ro1(1,jkz,jky)= cm * ro1(1,jkz,jky)
          ro1(2,jkz,jky)= cm * ro1(2,jkz,jky)
c
        end do
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine make_coulomb_smoothc_2D_k0(
     &    cx,         cy,         cz,         ix,         npwx,
     &    npwy,       npwz,       nx,         ny,         nz,
     &    D,          ro1)
c
      implicit real*8(A-H,O-Z)
      real*8 ro1(2,npwz,npwy)
c
c     kx=0
c     valksqx=cx*dfloat(kx)**2
      do jky=1,npwy
        ky=(jky-1)-npwy*int(jky/ny)
c       valksqxy=valksqx+cy*dfloat(ky)**2
        valksqxy=cy*dfloat(ky)**2
        do jkz=1,npwz
          kz=(jkz-1)-npwz*int(jkz/nz)
c
          if((ky .eq. 0) .and. (kz .eq. 0)) then
            cm = 6.2831853071795864769d0 * D**2
          else
            valksq=valksqxy + cz*dfloat(kz)**2
            valkabs=sqrt(valksq)
c
c                       4*%pi
            cm=(1.2566370614359172954d1/valksq) * (1.0d0-cos(valkabs*D))
          end if
c
            ro1(1,jkz,jky)= cm * ro1(1,jkz,jky)
            ro1(2,jkz,jky)= cm * ro1(2,jkz,jky)
c
        end do
      end do
c
      return
      end
c
c******************************************************************
c
      Subroutine make_COULOMB_on_FTGRID(
     &    npwxo,      npwyo,      npwzo,      npwxe,      npwye,
     &    npwze,      ftwork,     D,          Lxe,        Lye,
     &    Lze,        bfk,        w1,         vyze,       w3,
     &    w4,         vky,        vkz)
c
c This subroutine calculates the exact Coulomb mx. elements on the
c Fourier grid using the density and the value of D which is
c the maximum distance in the Coulomb interaction.
c
c The routine uses the real and complex one dimensional
c FFT subroutines
c
c This subroutine is written by Laszlo Fusti-Molnar at U of A
c 2001 August. Modified by MM, March-April 2005 (now the fftpack
c routines are used for fft. Also got rid of the complex arrays and
c incorporated the code for the Coulomb operator calculation)
c
c Input:
c     ftwork(npwzo,npwyo,npwxo) - array with the density
c     npwzo,npwyo,npwxo         - grid dimensions of the original box
c                                 ALL MUST BE EVEN NUMBER !!!
c     npwxe,npwxe,npwxe         - expanded grid dimensions
c                                 ALL MUST BE EVEN NUMBER !!!
c     D                         - max. possible  distance in the Coulomb
c                                 interaction
c     Lxe,Lye,Lze               - expanded box sizes
c
c Working arrays:
c
c       bfk(npwzo,npwyo,npwxe),w1(*),w3(*),w4(*),
c       vky(npwye),vkz(npwze)vyze(2,npwze,npwye)
c
c Output:
c       ftwork(npwzo,npwyo,npwxo) - array with the Coulomb mx element
c                                   on the grid points multiplied by
c                                   npwxe*npwxe*npwxe
c
      implicit none
c
      integer npwxo,npwyo,npwzo,npwxe,npwye,npwze
      integer ixstart,ixend,iystart,iyend,izstart,izend
      integer ix,iy,iz,ixe,iye,ize,nx,ny,nz
      real*8 D,Lxe,Lye,Lze,cx,cy,cz
c
c arrays
c
      real*8 ftwork(npwzo,npwyo,npwxo)
      real*8 bfk(npwzo,npwyo,npwxe),w1(*),w3(*),w4(*)
      real*8 vyze(2,npwze,npwye)
      real*8 vky(npwye),vkz(npwze)
c
      real*8 twopi2
      parameter(twopi2=3.9478417604357434475d1) ! (2*%pi)^2
cmax
      real*8 t1,t2,tex,tin,txf,tzero,tfill,tyf,tzf,tc,tzb,tyb
      real*8 tor,txb,tfb
      real*8 vkx,vkxy,valksq,cm,cmd2
      integer kx,jky,ky,jkz,kz
cmax
c
c   precompute some quantities for the Coulomb potential calculation
c
      cmd2 = 6.2831853071795864769d0 * D * D
      cx=twopi2/(Lxe*Lxe)
      cy=twopi2/(Lye*Lye)
      cz=twopi2/(Lze*Lze)
      nx=npwxe/2+2
      ny=npwye/2+2
      nz=npwze/2+2
      do jky=1,npwye
        ky=(jky-1)-npwye*int(jky/ny)
        vky(jky)=cy*dfloat(ky*ky)
      enddo
      do jkz=1,npwze
        kz=(jkz-1)-npwze*int(jkz/nz)
        vkz(jkz)=cz*dfloat(kz*kz)
      enddo
c
      ixstart=(npwxe-npwxo)/2 + 1
      ixend=ixstart+npwxo - 1
      iystart=(npwye-npwyo)/2 + 1
      iyend=iystart+npwyo - 1
      izstart=(npwze-npwzo)/2 + 1
      izend=izstart+npwzo - 1
c
c  expand the density on the extended grid
c
cc      call secund(t1)
      do ix=1,npwxo
        ixe=ixstart+ix-1
        do iy=1,npwyo
          do iz=1,npwzo
            bfk(iz,iy,ixe)=ftwork(iz,iy,ix)   !this should be
                                              !deleted if they overlap
          end do
        end do
      end do
cc      call secund(t2)
cc      tex=t2-t1
c
c    initialize the FFT
c
cc      call secund(t1)
      call vrffti(npwxe,w1)
      call vzffti(npwye,w3)
      call vzffti(npwze,w4)
      call secund(t2)
cc      tin=t2-t1
c
c for the new development: real_to_complex trafo of the real array
c
cc      call secund(t1)
      call vrfftf(1,npwzo*npwyo,1,npwxe,bfk,npwzo*npwyo,w1)
      call secund(t2)
cc      txf=t2-t1
c
c
c use vyze to create an expanded 2D function for every ix
c
cc      tzero=0.0d0
cc      tfill=0.0d0
cc      tyf=0.0d0
cc      tzf=0.0d0
cc      tc=0.0d0
cc      tzb=0.0d0
cc      tyb=0.0d0
cc      tor=0.0d0
c
      DO ix=1,npwxe/2+1 ! main loop over kx
        kx=(ix-1)-npwxe*int(ix/nx)
        vkx=cx*dfloat(kx*kx)
c
c set to zero the elements of vyze belonging to the expanded grid only
c
cc      call secund(t1)
       do iye=1,iystart-1
          call zeroit(vyze(1,1,iye),2*npwze)
       enddo
       do iye=iystart,iyend
          do ize=1,izstart-1
             vyze(1,ize,iye)=0.0d0
             vyze(2,ize,iye)=0.0d0
          enddo
          do ize=izend+1,npwze
             vyze(1,ize,iye)=0.0d0
             vyze(2,ize,iye)=0.0d0
          enddo
       enddo
       do iye=iyend+1,npwye
          call zeroit(vyze(1,1,iye),2*npwze)
       enddo
cc      call secund(t2)
cc      tzero=tzero+t2-t1
c
c  now fill the elements belonging to both original and expanded grids
c
cc      call secund(t1)
       if(ix.eq.1)then
          do iy=1,npwyo
            iye=iystart+iy-1
            do iz=1,npwzo
              ize=izstart+iz-1
              vyze(1,ize,iye)=bfk(iz,iy,1)
              vyze(2,ize,iye)=0.0d0
            end do
          end do
       else if(ix.eq.npwxe/2+1)then
          do iy=1,npwyo
            iye=iystart+iy-1
            do iz=1,npwzo
              ize=izstart+iz-1
              vyze(1,ize,iye)=bfk(iz,iy,npwxe)
              vyze(2,ize,iye)=0.0d0
            end do
          end do
       else
          do iy=1,npwyo
            iye=iystart+iy-1
            do iz=1,npwzo
              ize=izstart+iz-1
              vyze(1,ize,iye)=bfk(iz,iy,2*ix-2)
              vyze(2,ize,iye)=bfk(iz,iy,2*ix-1)
            end do
          end do
       endif
cc      call secund(t2)
cc      tfill=tfill+t2-t1
c
c FT along y for every npwzo sequence
c
cc        call secund(t1)
        call vzfftf(izstart,izend,1,npwye,vyze,npwze,w3)
cc        call secund(t2)
cc        tyf=tyf+t2-t1
c
c FT along z
c
cc        call secund(t1)
        do iy=1,npwye
          call zfftf(npwze,vyze(1,1,iy),w4)
        enddo
cc        call secund(t2)
cc        tzf=tzf+t2-t1
c
c Now we have everything in momentum space at this point
c Build up the Coulomb op. in momentum space
c
cc      call secund(t1)
        if(ix .eq. 1) then
          do jky=1,npwye
            do jkz=1,npwze
              if((vky(jky) .eq. 0.0d0) .and.
     &           (vkz(jkz) .eq. 0.0d0)) then
                vyze(1,jkz,jky)= cmd2 * vyze(1,jkz,jky)
                vyze(2,jkz,jky)= cmd2 * vyze(2,jkz,jky)
              else
                valksq=vky(jky)+vkz(jkz)
c                           4*%pi
                cm=(1.2566370614359172954d1/valksq)
     &                     * (1.0d0-cos(sqrt(valksq)*D))
                vyze(1,jkz,jky)= cm * vyze(1,jkz,jky)
                vyze(2,jkz,jky)= cm * vyze(2,jkz,jky)
              end if
            enddo
          enddo
        else
          do jky=1,npwye
            vkxy=vkx+vky(jky)
            do jkz=1,npwze
              valksq=vkxy+vkz(jkz)
c                          4*%pi
              cm= (1.2566370614359172954d1/valksq)
     &                    * (1.0d0-cos(sqrt(valksq)*D))
              vyze(1,jkz,jky)= cm * vyze(1,jkz,jky)
              vyze(2,jkz,jky)= cm * vyze(2,jkz,jky)
            enddo
          enddo
        end if
cc      call secund(t2)
cc      tc=tc+t2-t1
c
c FTBACK along z
c
cc        call secund(t1)
        do iy=1,npwye
          call zfftb(npwze,vyze(1,1,iy),w4)
        enddo
cc        call secund(t2)
cc        tzb=tzb+t2-t1
c
c FTBACK along y
c
cc        call secund(t1)
        call vzfftb(izstart,izend,1,npwye,vyze,npwze,w3)
cc        call secund(t2)
cc        tyb=tyb+t2-t1
c
c Now order back for complex_to_real trafo
c
cc        call secund(t1)
        if(ix.eq.1)then
           do iy=1,npwyo
             iye=iystart+iy-1
             do iz=1,npwzo
               ize=izstart+iz-1
               bfk(iz,iy,1)=vyze(1,ize,iye)
             end do
           end do
        else if(ix.eq.npwxe/2+1)then
           do iy=1,npwyo
             iye=iystart+iy-1
             do iz=1,npwzo
               ize=izstart+iz-1
               bfk(iz,iy,npwxe)=vyze(1,ize,iye)
             end do
           end do
        else
           do iy=1,npwyo
             iye=iystart+iy-1
             do iz=1,npwzo
               ize=izstart+iz-1
               bfk(iz,iy,2*ix-2)=vyze(1,ize,iye)
               bfk(iz,iy,2*ix-1)=vyze(2,ize,iye)
             end do
           end do
        endif
cc        call secund(t2)
cc        tor=tor+t2-t1
c
      ENDDO !loop over kx end
c
c Now back trafo along x
c
cc      call secund(t1)
      call vrfftb(1,npwzo*npwyo,1,npwxe,bfk,npwzo*npwyo,w1)
cc      call secund(t2)
cc      txb=t2-t1
c
cc      call secund(t1)
      do ix=1,npwxo
        ixe=ixstart+ix-1
        do iy=1,npwyo
          do iz=1,npwzo
            ftwork(iz,iy,ix)=bfk(iz,iy,ixe) !this line should
                                            !be deleted if they overlap
          end do
        end do
      end do
cc      call secund(t2)
cc      tfb=t2-t1
cc      write(*,*)
cc      write(*,*)'     Timing for Coulomb potential routine'
cc      write(*,*)
cc      write(*,*)' original grid dimensions: ',npwxo,npwyo,npwzo
cc      write(*,*)' extended grid dimensions: ',npwxe,npwye,npwze
cc      write(*,*)
cc      write(*,*)'Density expansion          =',tex
cc      write(*,*)'FFT initialization         =',tin
cc      write(*,*)'Forward FFT along x        =',txf
cc      write(*,*)'Zeroing 2D array           =',tzero
cc      write(*,*)'Filling 2D array           =',tfill
cc      write(*,*)'Forward FFT along y        =',tyf
cc      write(*,*)'Forward FFT along z        =',tzf
cc      write(*,*)'Coulomb potential          =',tc
cc      write(*,*)'Backward FFT along z       =',tzb
cc      write(*,*)'Backward FFT along y       =',tyb
cc      write(*,*)'Ordering Back              =',tor
cc      write(*,*)'Backward FFT along x       =',txb
cc      write(*,*)'Filling back the potential =',tfb
cc      write(*,*)
cc      write(*,*)'Total for FFT              =',
cc     &                                       tin+txf+tyf+tzf+tzb+tyb+txb
cc      write(*,*)'Total                      =',tex+tin+txf+tzero+tfill+
cc     &                                tyf+tzf+tc+tzb+tyb+tor+txb+tfb
cc      write(*,*)
c
c
      return
      end
c
c******************************************************************
c
      Subroutine anal_sharp_ongrid_new(
     &    itype,      eta,        concoef,    icomp,      cji,
     &    Rx,         Ry,         Rz,         npwx,       npwy,
     &    npwz,       Lx,         Ly,         Lz,         ftwork)
c
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      complex*16 ftwork(npwz,npwy,npwx)
      real*8 Lx,Ly,Lz
      integer itype,icomp,kx,ky,kz
      data pi/3.1415926535897932384626433d0/
c
      ax = 2.0d0 * pi / Lx
      ay = 2.0d0 * pi / Ly
      az = 2.0d0 * pi / Lz
c
      Rx=Rx+Lx/2.0d0 !In FFTW the origo is different !
      Ry=Ry+Ly/2.0d0
      Rz=Rz+Lz/2.0d0
      Rxax=Rx*ax
      Ryay=Ry*ay
      Rzaz=Rz*az
c
ct      const1=cc*sqrt(sqrt((2.0d0 * eta / pi)))**3  * cji * concoef
      const1=sqrt(sqrt((2.0d0 * eta / pi)))**3  * cji * concoef
      tfml1=-1.0d0/eta/4.0d0
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      if(itype .eq. 1) then
c
c  s functions
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = cos(arg)
              val=const1*t22
              valre = val*t30
              t30 = sin(arg)
              valim = -val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c    End of s type functions
c
      else if ((itype .eq. 2) .and. (icomp .eq. 1)) then
c
c Px function
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*ax/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          const1v=const1*kx !this is different from the s
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = sin(arg) !this is different
              val=-const1v*t22 !this is different
              valre = val*t30 !this is different
              t30 = cos(arg) !this is different
              valim = val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 2)) then
c
c Py function
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*ay/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            const1v=const1*ky !this is different from the s
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = sin(arg) !this is different
              val=-const1v*t22 !this is different
              valre = val*t30 !this is different
              t30 = cos(arg) !this is different
              valim = val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 3)) then
c
c Pz function
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*az/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
            kz=(jkz-1)-npwz*int(jkz/nz)
            const1v=const1*kz !this is different from the s
            t17 = kz**2
            t22 = tfml3*exp(tfml6*t17)
c            t22 = tfml3*exptfml6**t17
            arg=argxy+Rzaz*kz
            t30 = sin(arg) !this is different
            val=-const1v*t22 !this is different
            valre = val*t30 !this is different
            t30 = cos(arg) !this is different
            valim = val*t30
c
            ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 1)) then
c
c  s part of l
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = cos(arg)
              val=const1*t22
              valre = val*t30
              t30 = sin(arg)
              valim = -val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 2)) then
c
c Px part of l
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*ax/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          const1v=const1*kx !this is different from the s
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = sin(arg) !this is different
              val=-const1v*t22 !this is different
              valre = val*t30 !this is different
              t30 = cos(arg) !this is different
              valim = val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 3)) then
c
c Py part of l
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*ay/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            const1v=const1*ky !this is different from the s
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = sin(arg) !this is different
              val=-const1v*t22 !this is different
              valre = val*t30 !this is different
              t30 = cos(arg) !this is different
              valim = val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                            dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 4)) then
c
c Pz part of l
c
        t1 = sqrt(eta)
        t2 = t1**2
        t5 = sqrt(pi)
        t6 = t5**2
        t10 = ax**2
        t13 = ay**2
        t16 = az**2
c
        const1=const1/t2/t1*t6*t5
        tfml4=tfml1*t10
        tfml5=tfml1*t13
        tfml6=tfml1*t16
        exptfml6=exp(tfml6)
        const1=const1*az/2.0d0/eta !this is different from the s
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          t11 = kx**2
          tfml2=exp(tfml4*t11)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            t14 = ky**2
            tfml3=tfml2*exp(tfml5*t14)
            do jkz=1,npwz
              kz=(jkz-1)-npwz*int(jkz/nz)
              const1v=const1*kz !this is different from the s
              t17 = kz**2
              t22 = tfml3*exp(tfml6*t17)
c              t22 = tfml3*exptfml6**t17
              arg=argxy+Rzaz*kz
              t30 = sin(arg) !this is different
              val=-const1v*t22 !this is different
              valre = val*t30 !this is different
              t30 = cos(arg) !this is different
              valim = val*t30
c
              ftwork(jkz,jky,jkx)=ftwork(jkz,jky,jkx)+
     &                             dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else
        call nerror(1,'anal_sharp_ongrid_new',
     $    'only S, P or l',0,0)
      end if
      return
      end
c
c***********************************************************************
c
      Subroutine anal_sharp_derivatives_ongrid(
     &    itype,         etas,          concoefs,      icomp,
     &    Rx,            Ry,            Rz,            npwx,
     &    npwy,          npwz,          Lx,            Ly,
     &    Lz,            ftwork,        ncontr,        expfuncx,
     &    expfuncy,      expfuncz,      nzeffdim,      ftwork2)
c
c This subroutine buid up the analytical Fourier transformation of the
c derivative of the basis functions and put it into the ftwork array.
c
c Inputs:
c itype            type of the given basis function (s,p,l,d etc.)
c etas             precalculated 1/eta values for tha whole contraction
c concoefs         modified contr. coefs for the whole contraction
c icomp            component nr.
c Rx,Ry,Rz         modified (for fftw) centers of the given shell
c npwx,npwy,npwz   grid dimensions
c nzeffdim         effective dimension along z (npwz/2+1)
c Lx,Ly,Lz         the given box dimensions
c ncontr           dimension of the given contraction
c expfuncx,expfuncy,expfuncz    pre calculated exponential functions
c
c In/Outs
c ftwork2          array with the helping function
c
c Outputs
c ftwork           array for the derivatives:
c                       ftwork(1,n) -> d/dRx;
c                       ftwork(2,n) -> d/dRy;
c                       ftwork(3,n) -> d/dRz;
c
c
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      complex*16 ftwork(3,nzeffdim,npwy,npwx) ! x,y,z comp. of the
                                              ! derivatives
      complex*16 ftwork2(nzeffdim,npwy,npwx)
      real*8 etas(ncontr)
      real*8 concoefs(ncontr)
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 Lx,Ly,Lz
      complex*16 vals
      integer itype,icomp,kx,ky,kz,i
      data pi/3.1415926535897932384626433d0/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' Data on entry to <anal_sharp_derivatives_ongrid>'
cc      write(6,*) ' itype:',itype,' icomp:',icomp,' ncontr:',ncontr
cc      write(6,*) ' Rx:',rx,' Ry:',ry,' Rz:',rz
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' Lx:',lx,' Ly:',ly,' Lz:',lz
cc      write(6,*) ' etas and concoefs array are:'
cc      do i=1,ncontr
cc      write(6,*) i,'  ',etas(i),concoefs(i)
cc      enddo
cc      write(6,*) ' expfuncx array is:'
cc      do i=1,npwx
cc      write(6,*) i,'  ',expfuncx(i)
cc      enddo
cc      write(6,*) ' expfuncy array is:'
cc      do i=1,npwy
cc      write(6,*) i,'  ',expfuncy(i)
cc      enddo
cc      write(6,*) ' expfuncz array is:'
cc      do i=1,npwz
cc      write(6,*) i,'  ',expfuncz(i)
cc      enddo
cc      write(6,*) ' nzeffdim is:',nzeffdim
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
      ax = 2.0d0 * pi / Lx
      ay = 2.0d0 * pi / Ly
      az = 2.0d0 * pi / Lz
c
      Rx=Rx+Lx/2.0d0 !In FFTW the origo is different !
      Ry=Ry+Ly/2.0d0
      Rz=Rz+Lz/2.0d0
      Rxax=Rx*ax
      Ryay=Ry*ay
      Rzaz=Rz*az
c
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      if(itype .eq. 1) then
c
c  s function derivatives: first calculate a helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = -sin(arg)
              t31 = -cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork(1,jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c now calculate the derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
            end do
          end do
        end do
cc
      else if((itype .eq. 2) .and. (icomp .eq. 1)) then
c
c  Px function derivatives: first calculate  a helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = -cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c now copy the helping function (in ftwork2) to the first comp.
c of ftwork
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the Px derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          cx2=cx**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            cxtcy=cx*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cxtcz=cx*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx2*vals
              ftwork(2,jkz,jky,jkx)=cxtcy*vals
              ftwork(3,jkz,jky,jkx)=cxtcz*vals
            end do
          end do
        end do
c
      else if((itype .eq. 2) .and. (icomp .eq. 2)) then
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the Py derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            cy2=cy**2
            cxtcy=cx*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cytcz=cy*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cxtcy*vals
              ftwork(2,jkz,jky,jkx)=cy2*vals
              ftwork(3,jkz,jky,jkx)=cytcz*vals
            end do
          end do
        end do
c
      else if((itype .eq. 2) .and. (icomp .eq. 3)) then
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the Pz derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cz2=cz**2
              cytcz=cy*cz
              cxtcz=cx*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cxtcz*vals
              ftwork(2,jkz,jky,jkx)=cytcz*vals
              ftwork(3,jkz,jky,jkx)=cz2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 3) .and. (icomp .eq. 1)) then
c
c  s of L derivatives: first calculate the helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = -sin(arg)
              t31 = -cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork(1,jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c now   calculate the derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
            end do
          end do
        end do
c
      else if((itype .eq. 3) .and. (icomp .eq. 2)) then
c
c  Px of L derivatives: first calculate  a helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = -cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the derivatives
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          cx2=cx**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            cxtcy=cx*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cxtcz=cx*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx2*vals
              ftwork(2,jkz,jky,jkx)=cxtcy*vals
              ftwork(3,jkz,jky,jkx)=cxtcz*vals
            end do
          end do
        end do
c
      else if((itype .eq. 3) .and. (icomp .eq. 3)) then
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the derivatives of Py of L
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            cy2=cy**2
            cxtcy=cx*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cytcz=cy*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cxtcy*vals
              ftwork(2,jkz,jky,jkx)=cy2*vals
              ftwork(3,jkz,jky,jkx)=cytcz*vals
            end do
          end do
        end do
c
      else if((itype .eq. 3) .and. (icomp .eq. 4)) then
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
c now calculate the derivatives of Pz of L
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*dfloat(kx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*dfloat(ky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*dfloat(kz)
              cz2=cz**2
              cytcz=cy*cz
              cxtcz=cx*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cxtcz*vals
              ftwork(2,jkz,jky,jkx)=cytcz*vals
              ftwork(3,jkz,jky,jkx)=cz2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 4) .and. (icomp .eq. 1)) then
c
c now calculate the derivatives of the D1 component (2dzz-dxx-dyy)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=-cx*cx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2xy=c2x-cy*cy
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2=c2xy+2.0d0*cz*cz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 4) .and. (icomp .eq. 2)) then
c
c now calculate the derivatives of the D2 component (dxx-dyy)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2xy=c2x-cy*cy
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*c2xy*vals
              ftwork(2,jkz,jky,jkx)=cy*c2xy*vals
              ftwork(3,jkz,jky,jkx)=cz*c2xy*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 4) .and. (icomp .eq. 3)) then
c
c now calculate the derivatives of the D3 component (dxy)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            cxy=cx*cy
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*cxy*vals
              ftwork(2,jkz,jky,jkx)=cy*cxy*vals
              ftwork(3,jkz,jky,jkx)=cz*cxy*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 4) .and. (icomp .eq. 4)) then
c
c now calculate the derivatives of the D4 component (dxz)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              cxz=cx*cz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*cxz*vals
              ftwork(2,jkz,jky,jkx)=cy*cxz*vals
              ftwork(3,jkz,jky,jkx)=cz*cxz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 4) .and. (icomp .eq. 5)) then
c
c now calculate the derivatives of the D5 component (dyz)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              cyz=cy*cz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*cyz*vals
              ftwork(2,jkz,jky,jkx)=cy*cyz*vals
              ftwork(3,jkz,jky,jkx)=cz*cyz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 1)) then
c
c now calculate the derivatives of the DXX component
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+(c2x-2.0d0/etas(icontr))*
     &            concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 2)) then
c
c now calculate the derivatives of the DYY component
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+(c2y-2.0d0/etas(icontr))*
     &            concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 3)) then
c
c now calculate the derivatives of the DZZ component
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+(c2z-2.0d0/etas(icontr))*
     &            concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              vals=dcmplx(valre,valim)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 4)) then
c
c  DXY, DXZ, DYZ function derivatives: first calculate  a helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c now calculate the derivatives of the DXY component
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            cxy=cx*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
c
              vals=cxy*ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 5)) then
c
c now calculate the derivatives of the DXZ component
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              cxz=cx*cz
c
              vals=cxz*ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 5) .and. (icomp .eq. 6)) then
c
c now calculate the derivatives of the DYZ component
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              cyz=cy*cz
c
              vals=cyz*ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*vals
              ftwork(2,jkz,jky,jkx)=cy*vals
              ftwork(3,jkz,jky,jkx)=cz*vals
c
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 1)) then
c
c  F function derivatives: first calculate  a helping function
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = -cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = val*t31
c
              ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c now calculate the derivatives of the F1 component (-4fxxy+fyyy+fyzz)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=-cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2x=-cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cy*(-4.0d0*c2x+c2y+c2z)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 2)) then
c
c now calculate the derivatives of the F2 component (-4fxxz+fyyz+fzzz)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cz*(-4.0d0*c2x+c2y+c2z)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 3)) then
c
c now calculate the derivatives of the F3 component (-4fxyy+fxxx+fxzz)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cx*(-4.0d0*c2y+c2x+c2z)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 4)) then
c
c now calculate the derivatives of the F4 component (-4fyyz+fxxz+fzzz)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cz*(-4.0d0*c2y+c2x+c2z)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 5)) then
c
c now calculate the derivatives of the F5 component (-4fxzz+fxxx+fxyy)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cx*(-4.0d0*c2z+c2x+c2y)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 6)) then
c
c now calculate the derivatives of the F6 component (-4fyzz+fxxy+fyyy)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          c2x=cx*cx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            c2y=cy*cy
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2z=cz*cz
              c2=cy*(-4.0d0*c2z+c2x+c2y)
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else if((itype .eq. 6) .and. (icomp .eq. 7)) then
c
c now calculate the derivatives of the F7 component (fxyz)
c
        call zcopy(npwx*npwy*nzeffdim,ftwork2,1,ftwork,3)
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          cx=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            cy=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              cz=az*kz
              c2=-cx*cy*cz
              vals=ftwork(1,jkz,jky,jkx)
              ftwork(1,jkz,jky,jkx)=cx*c2*vals
              ftwork(2,jkz,jky,jkx)=cy*c2*vals
              ftwork(3,jkz,jky,jkx)=cz*c2*vals
            end do
          end do
        end do
c
      else
        call nerror(1,'anal_sharp_derivatives_ongrid',
     $    'only S,P,l,D,d6,F',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
      Subroutine anal_sharp_ongrid_new2(
     &    itype,      etas,       concoefs,   icomp,      Rx,
     &    Ry,         Rz,         npwx,       npwy,       npwz,
     &    Lx,         Ly,         Lz,         ftwork,     ncontr,
     &    expfuncx,   expfuncy,   expfuncz,   nzeffdim)
c
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      complex*16 ftwork(nzeffdim,npwy,npwx)
      real*8 etas(ncontr)
      real*8 concoefs(ncontr)
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 Lx,Ly,Lz
      integer itype,icomp,kx,ky,kz
      data pi/3.1415926535897932384626433d0/
c
      ax = 2.0d0 * pi / Lx
      ay = 2.0d0 * pi / Ly
      az = 2.0d0 * pi / Lz
c
      Rx=Rx+Lx/2.0d0 !In FFTW the origo is different !
      Ry=Ry+Ly/2.0d0
      Rz=Rz+Lz/2.0d0
      Rxax=Rx*ax
      Ryay=Ry*ay
      Rzaz=Rz*az
c
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      if(itype .eq. 1) then
c
c  s functions
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = -val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else
        call nerror(1,'anal_sharp_ongrid_new2',
     $    'only S here!!!',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
      Subroutine anal_sharp_ongrid_new3(
     &    itype,         etas,          concoefs,      icomp,
     &    Rx,            Ry,            Rz,            npwx,
     &    npwy,          npwz,          Lx,            Ly,
     &    Lz,            ftwork,        ncontr,        expfuncx,
     &    expfuncy,      expfuncz,      ftwork2,       nzeffdim)
c
c
c This subroutine calculates the analytical Fourier coefficients
c [int(f(r-R)*exp(-I*k*r))] of a given basis fuction.
c It can handle s,p,l,d,d6,f functions.
c
c Some additional usefull information for these basis functions:
c
c P: The order is x,y,z and 2*sqrt(eta) is inside the concoef.
c
c D: The order is: (3z^2-r^2),(x^2-y^2),xy,xz,yz.
c    4*eta inside the concoef
c
c D6: The order is: x^2,y^2,z^2,xy,xz,yz. The first three is normalized
c     to 3, the others to 1. 4*eta inside the concoef.
c
c F:  The order is: (5x^2*y-r^2*y),(5x^2*z-r^2*z),(5y^2*x-r^2*x),
c     (5y^2*z-r^2*z),(5z^2*x-r^2*x),(5z^2*y-r^2*y),x*y*z.
c     {(2*eta)^(3/2)}/sqrt(5) inside the concoef.
c
c Input:
c        itype:        type of the shell
c        etas:         1/eta (1/orbital exponent) for all contractions
c        icomp:        the given component of the shell
c        concoefs:     contraction coefs for all contractions
c                      times some constants
c        Rx,Ry,Rx      center of the given basis function
c        npwx,npwy,npwz Number of plane-wave in each direction
c        Lx,Ly,Lz      Box lengths for each direction
c        ncontr        Number of the contracted function in the
c                      given contraction
c        expfuncx,
c        expfuncy,     Precalculated exp functions
c        expfuncz
c
c  InOut:
c        ftwork2       helping array
c
c  Output:
c        ftwork:       The result array with the analitical
c                      fourier components
c
c
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      complex*16 ftwork(nzeffdim,npwy,npwx)
      complex*16 ftwork2(nzeffdim,npwy,npwx)
      real*8 etas(ncontr)
      real*8 concoefs(ncontr)
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 Lx,Ly,Lz
      integer itype,icomp,kx,ky,kz
      data pi/3.1415926535897932384626433d0/
c
      ax = 2.0d0 * pi / Lx
      ay = 2.0d0 * pi / Ly
      az = 2.0d0 * pi / Lz
c
      Rx=Rx+Lx/2.0d0   ! In FFTW the origin is different
      Ry=Ry+Ly/2.0d0
      Rz=Rz+Lz/2.0d0
      Rxax=Rx*ax
      Ryay=Ry*ay
      Rzaz=Rz*az
c
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      if(itype .eq. 1) then
c
c  s functions
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = -val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  Px functions
c
      else if ((itype .eq. 2) .and. (icomp .eq. 1)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
        valexpx=expfuncx(jkx)
          c2=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       call zcopy(npwx*npwy*nzeffdim,ftwork,1,ftwork2,1) !save Px into
                                                          ! ftwork2
c
c now recalculate the Kx=0 without ax*kx
c
c       do jkx=1,1
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         argx=Rxax*kx
c         valexpx=expfuncx(jkx)
c        c2=1.0d0
c         do jky=1,npwy
c           ky=(jky-1)-npwy*int(jky/ny)
c           argxy=argx+Ryay*ky
c           valexpxy=valexpx*expfuncy(jky)
c           do jkz=1,nzeffdim
c             kz=(jkz-1)-npwz*int(jkz/nz)
c             arg=argxy+Rzaz*kz
c             valexp=valexpxy*expfuncz(jkz)
c
c             val=0.0d0
c             t30 = sin(arg)
c             t31 = cos(arg)
c
c             do icontr=1,ncontr
c               val=val+concoefs(icontr)*valexp**etas(icontr)
c             end do
c
c             valre = -val*t30
c             valim = -val*t31
c
c             ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
c           end do
c         end do
c       end do
c
c
c  Py functions
c
      else if ((itype .eq. 2) .and. (icomp .eq. 2)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       do jkx=1,npwx
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         if(jkx .gt. 1) then
c           c2x=1.0d0/ax/kx
c         else
c           c2x=1.0d0
c         end if
c         do jky=1,npwy
c           ky=(jky-1)-npwy*int(jky/ny)
c           c2=c2x*ay*ky
c           do jkz=1,nzeffdim
c
c             ftwork(jkz,jky,jkx)=c2*ftwork2(jkz,jky,jkx)
c
c           end do
c         end do
c       end do
c
c  Pz functions
c
      else if ((itype .eq. 2) .and. (icomp .eq. 3)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=az*kz
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       do jkx=1,npwx
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         if(jkx .gt. 1) then
c           c2x=1.0d0/ax/kx
c         else
c           c2x=1.0d0
c         end if
c         do jky=1,npwy
c           do jkz=1,nzeffdim
c             kz=(jkz-1)-npwz*int(jkz/nz)
c             c2=c2x*az*kz
c
c             ftwork(jkz,jky,jkx)=c2*ftwork2(jkz,jky,jkx)
c
c           end do
c         end do
c       end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 1)) then

c
c  s part of L
c
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = val*t30
              valim = -val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  Px part of L
c
      else if ((itype .eq. 3) .and. (icomp .eq. 2)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       call dcopy(2*npwx*npwy*nzeffdim,ftwork,1,ftwork2,1) !save Px into ftwork2
c     write(79,*)'ftwork2 after copy'
c     do iz=1,nzeffdim
c     do iy=1,npwy
c     do ix=1,npwx
c     write(79,'(2e20.5)')real(ftwork2(iz,iy,ix)),
c    &                    imag(ftwork2(iz,iy,ix))
c     enddo
c     enddo
c     enddo
c
c now   recalculate the Kx=0 without ax*kx
c
c       do jkx=1,1
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         argx=Rxax*kx
c         valexpx=expfuncx(jkx)
c        c2=1.0d0
c         do jky=1,npwy
c           ky=(jky-1)-npwy*int(jky/ny)
c           argxy=argx+Ryay*ky
c           valexpxy=valexpx*expfuncy(jky)
c           do jkz=1,nzeffdim
c             kz=(jkz-1)-npwz*int(jkz/nz)
c             arg=argxy+Rzaz*kz
c             valexp=valexpxy*expfuncz(jkz)
c
c             val=0.0d0
c             t30 = sin(arg)
c             t31 = cos(arg)
c
c             do icontr=1,ncontr
c               val=val+concoefs(icontr)*valexp**etas(icontr)
c             end do
c
c             valre = -val*t30
c             valim = -val*t31
c
c             ftwork2(jkz,jky,jkx)=dcmplx(valre,valim)
c
c           end do
c         end do
c       end do
c     write(79,*)'ftwork2 after recalculating kx=0'
c     do iz=1,nzeffdim
c     do iy=1,npwy
c     do ix=1,npwx
c     write(79,'(2e20.5)')real(ftwork2(iz,iy,ix)),
c    &                    imag(ftwork2(iz,iy,ix))
c     enddo
c     enddo
c     enddo
c
c  Py part of L
c
      else if ((itype .eq. 3) .and. (icomp .eq. 3)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       do jkx=1,npwx
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         if(jkx .gt. 1) then
c           c2x=1.0d0/ax/kx
c         else
c           c2x=1.0d0
c         end if
c         do jky=1,npwy
c           ky=(jky-1)-npwy*int(jky/ny)
c           c2=c2x*ay*ky
c           do jkz=1,nzeffdim
c
c             ftwork(jkz,jky,jkx)=c2*ftwork2(jkz,jky,jkx)
c
c           end do
c         end do
c       end do
c
c  Pz part of L
c
      else if ((itype .eq. 3) .and. (icomp .eq. 4)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=az*kz
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c       do jkx=1,npwx
c         kx=(jkx-1)-npwx*int(jkx/nx)
c         if(jkx .gt. 1) then
c           c2x=1.0d0/ax/kx
c         else
c           c2x=1.0d0
c         end if
c         do jky=1,npwy
c           do jkz=1,nzeffdim
c             kz=(jkz-1)-npwz*int(jkz/nz)
c             c2=c2x*az*kz
c
c             ftwork(jkz,jky,jkx)=c2*ftwork2(jkz,jky,jkx)
c
c           end do
c         end do
c       end do
c
c  D1 function
c
      else if ((itype .eq. 4) .and. (icomp .eq. 1)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=-ax**2*kx**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2xy=c2x-ay**2*ky**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2xy+2.0d0*az**2*kz**2
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D2 function
c
      else if ((itype .eq. 4) .and. (icomp .eq. 2)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=ax**2*kx**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=c2x-ay**2*ky**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D3 function
c
      else if ((itype .eq. 4) .and. (icomp .eq. 3)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=c2x*ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D4 function
c
      else if ((itype .eq. 4) .and. (icomp .eq. 4)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2x*az*kz
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D5 function
c
      else if ((itype .eq. 4) .and. (icomp .eq. 5)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2y=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2y*az*kz
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 first function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 1)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2=ax**2*kx**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+(c2-2.0d0/etas(icontr))*
     &            concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -val*t30
              valim = val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 second function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 2)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=ay**2*ky**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+(c2-2.0d0/etas(icontr))*
     &            concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -val*t30
              valim = val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 third function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 3)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=az**2*kz**2
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+(c2-2.0d0/etas(icontr))*
     &          concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -val*t30
              valim = val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 4-th function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 4)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2=c2x*ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 5-th function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 5)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2x*az*kz
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  D6 6-th function
c
      else if ((itype .eq. 5) .and. (icomp .eq. 6)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2y=ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2y*az*kz
c
              val=0.0d0
              t30 = cos(arg)
              t31 = sin(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: first functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 1)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=-4.0d0*c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x+c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3y*(c2xy+c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: second functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 2)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=-4.0d0*c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x+c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3z*(c2xy+c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: 3-th functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 3)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x-4.0d0*c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3x*(c2xy+c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: 4-th functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 4)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x-4.0d0*c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3z*(c2xy+c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: 5-th functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 5)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x+c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3x*(c2xy-4.0d0*c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
                end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: 6-th functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 6)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c3x=ax*kx
          c2x=c3x**2
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c3y=ay*ky
            c2xy=c2x+c3y**2
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c3z=az*kz
              c2=c3y*(c2xy-4.0d0*c3z**2)
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
c  F: 7-th functions
c
      else if ((itype .eq. 6) .and. (icomp .eq. 7)) then
c
        do jkx=1,npwx
          kx=(jkx-1)-npwx*int(jkx/nx)
          argx=Rxax*kx
          valexpx=expfuncx(jkx)
          c2x=-ax*kx
          do jky=1,npwy
            ky=(jky-1)-npwy*int(jky/ny)
            argxy=argx+Ryay*ky
            valexpxy=valexpx*expfuncy(jky)
            c2xy=c2x*ay*ky
            do jkz=1,nzeffdim
              kz=(jkz-1)-npwz*int(jkz/nz)
              arg=argxy+Rzaz*kz
              valexp=valexpxy*expfuncz(jkz)
              c2=c2xy*az*kz
c
              val=0.0d0
              t30 = sin(arg)
              t31 = cos(arg)
c
              do icontr=1,ncontr
                val=val+concoefs(icontr)*valexp**etas(icontr)
              end do
c
              valre = -c2*val*t30
              valim = -c2*val*t31
c
              ftwork(jkz,jky,jkx)=dcmplx(valre,valim)
c
            end do
          end do
        end do
c
      else
        call nerror(1,'anal_sharp_ongrid_new3',
     $ 'only s,p,l,d,d6,f sharps right now !!!',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
c
      Subroutine make_sharpgrids_work(
     &    itype,         icomp,         npwx,          npwy,
     &    npwz,          Lx,            Ly,            Lz,
     &    cftwork,       isharpgrd,     ncontr,        icfpl,
     &    plbasdat,      ncfpl,         xp,            yp,
     &    zp,            i4ftwork,      expfuncx,      expfuncy,
     &    expfuncz,      cftwork2,      comprcoefs,    ifcount,
     &    nzeffdim,      fbk,           w1,            w3,
     &    w4)
c
c This subroutine calculates the (exact) projection of the given
c compact basis function into the given Fourier grid space
c (coordinate space) using analytical FT followed by numerical back FFT.
c The subroutine writes this function into a file on the disk using
c 2D sequences.
c
c Inputs:
c          itype              type of the shell
c          icomp              the given component of the shell
c          npwx,npwy,npwz     Number of plane-wave in each direction
c          Lx,Ly,Lz           Box lengths for each direction
c          isharpgrd          unit number
c          ncontr             number of primitive shells in the given
c                             contraction
c          plbasdat           pw basdat: see before
c          ncfpl              number of primitive core-like shells
c          xp,yp,zp           coordinates of the origin of the actual
c                             box
c          expfuncx,
c          expfuncy,          Precalculated exp functions
c          expfuncz
c
c scratch arrays:
c
c          cftwork,cftwork2,i4ftwork,fbk,w1,w3,w4
c
c InOut:
c          icfpl              the starting primitive shell number
c
c Output: The values on the disk
c
c
c

      use memory
      use kinds     ! this is needed for i4ftwork

      implicit none
c
      integer itype,icomp,npwx,npwy,npwz,ix,iy,iz
      integer isharpgrd,ncontr,icfpl,ncfpl,icontr,icfpl2,iconcoefs
      integer iexponents
      integer ifcount,nzeffdim
      integer ncspl
      real*8 comprcoefs(*)
      real*8 vmax,value
      real*8 exponent,concoef,cji,Px,Py,Pz,Lx,Ly,Lz,cc,xp,yp,zp,summ
      real*8 plbasdat(6,ncfpl)
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 t1,t2,tfftb
      save tfftb
c     common /big/bl(300)
      real*8 pi
      data pi/3.1415926535897932384626433d0/
      real*8 ar1,ar2,ai1,ai2,diffr,diffi
      real*8 cftwork2(2,npwx,npwy)
      real*8 cftwork(2,nzeffdim,npwy,npwx)
      integer*4 i4ftwork(npwz,npwy)
      real*8 fbk(npwx,npwy,npwz),w1(*),w3(*),w4(*)
c
c Perform the analytical FT of the given function
c
      cji=1.0d0
c
      call mmark
      call getmem(ncontr,iconcoefs)
      call getmem(ncontr,iexponents)
c
      icfpl2=icfpl
      do icontr=1,ncontr
        icfpl2=icfpl2+1
        concoef=plbasdat(2,icfpl2)
        if((itype .eq. 3) .and. (icomp .gt. 1)) then
          concoef=plbasdat(6,icfpl2)
        end if
        exponent=plbasdat(1,icfpl2)
c
        bl(iexponents+icontr-1)=1.0d0/exponent !1/eta
c
        if(itype .eq. 1) then
c
          bl(iconcoefs+icontr-1)=
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                         !(2pi/eta)^(3/4)*concoef
c
        else if((itype .eq. 3) .and. (icomp .eq. 1)) then !s part of L
c
          bl(iconcoefs+icontr-1)=
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                          !(2pi/eta)^(3/4)*concoef
c
        else if((itype .eq. 2) .or. (itype .eq. 3)) then
                                              !P or P part of L
c
          bl(iconcoefs+icontr-1)=0.5d0 / exponent *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                      !0.5*(2pi/eta)^(3/4)*concoef/eta
c
c   *** D functions ***
c
        else if((itype .eq. 4) .and. (icomp .eq. 1)) then !D1
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef /sqrt(12.0d0)
                           !0.25*(2pi/eta)^(3/4)*concoef/eta^2/sqrt(12)
c
        else if((itype .eq. 4) .and. (icomp .eq. 2)) then !D2
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef /2.0d0
                                  !0.25*(2pi/eta)^(3/4)*concoef/eta^2/2
c
        else if(itype .eq. 4) then !D3,D4,D5
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                    !0.25*(2pi/eta)^(3/4)*concoef/eta^2
c
c   *** D6 functions ***
c
        else if(itype .eq. 5) then !D6
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                    !0.25*(2pi/eta)^(3/4)*concoef/eta^2
c
c   ***  F function  ***
c
        else if((itype .eq. 6) .and. (icomp .eq. 7)) then
                                                 ! 7-th component of f
c
          bl(iconcoefs+icontr-1)=0.125d0 * sqrt(40.0d0) / exponent**3 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                         !0.125*sqrt(40)*(2pi/eta)^(3/4)*concoef/eta^3
c
        else if(itype .eq. 6) then ! other f functions
c
          bl(iconcoefs+icontr-1)=0.125d0 / exponent**3 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                  !0.125*(2pi/eta)^(3/4)*concoef/eta^3
c
        end if
      end do
c
      Px=plbasdat(3,icfpl+1)-xp
      Py=plbasdat(4,icfpl+1)-yp
      Pz=plbasdat(5,icfpl+1)-zp
c
      if(itype .eq. 1) then
c
        call anal_sharp_ongrid_new2(
     &    itype,         bl(iexponents),bl(iconcoefs), icomp,
     &    Px,            Py,            Pz,            npwx,
     &    npwy,          npwz,          Lx,            Ly,
     &    Lz,            cftwork,       ncontr,        expfuncx,
     &    expfuncy,      expfuncz,      nzeffdim)
c
      else
c
        call anal_sharp_ongrid_new3(
     &    itype,         bl(iexponents),bl(iconcoefs), icomp,
     &    Px,            Py,            Pz,            npwx,
     &    npwy,          npwz,          Lx,            Ly,
     &    Lz,            cftwork,       ncontr,        expfuncx,
     &    expfuncy,      expfuncz,      cftwork2,       nzeffdim)
c
      end if
c
      call retmark
      icfpl=icfpl+ncontr
c
c  initialize the FFT
c
      call vrffti(npwz,w1)
      call vzffti(npwy,w3)
      call vzffti(npwx,w4)
c
      do iz=1,nzeffdim
        do iy=1,npwy
          do ix=1,npwx
            cftwork2(1,ix,iy)=cftwork(1,iz,iy,ix)
            cftwork2(2,ix,iy)=cftwork(2,iz,iy,ix)
          enddo
        enddo
c
c  back fft transform along x
c
        do iy=1,npwy
          call zfftb(npwx,cftwork2(1,1,iy),w4)
        enddo
c
c  back fft transform along y
c
        call vzfftb(1,npwx,1,npwy,cftwork2,npwx,w3)
        if(iz.eq.1)then
          do iy=1,npwy
            do ix=1,npwx
              fbk(ix,iy,1)=cftwork2(1,ix,iy)
            enddo
          enddo
        else if(mod(npwz,2).eq.0.and.iz.eq.npwz/2+1)then
          do iy=1,npwy
            do ix=1,npwx
              fbk(ix,iy,npwz)=cftwork2(1,ix,iy)
            enddo
          enddo
        else
          do iy=1,npwy
            do ix=1,npwx
              fbk(ix,iy,2*iz-2)=cftwork2(1,ix,iy)
              fbk(ix,iy,2*iz-1)=cftwork2(2,ix,iy)
            enddo
          enddo
        endif
      enddo
c
c  back fft along z
c
      call vrfftb(1,npwx*npwy,1,npwz,fbk,npwx*npwy,w1)
c
c   save the projected function into a file
c
      cc=1.0d0/Lx**2/Ly**2/Lz**2
      cc=sqrt(cc)
      do ix=1,npwx
        ifcount=ifcount+1
        vmax=0.0d0
        do iy=1,npwy
          do iz=1,npwz
            vmax=max(vmax,abs(fbk(ix,iy,iz)))
          enddo
        enddo
        vmax=cc*vmax
        if(vmax .gt. 0.0d0) then
          value=2147483647.0d0/vmax
        else
          value=0.0d0
        end if
        comprcoefs(ifcount)=value
        do iy=1,npwy
          do iz=1,npwz
            i4ftwork(iz,iy)=int(cc*value*fbk(ix,iy,iz),kind=i_4)
          end do
        end do
c
        write(isharpgrd) i4ftwork
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine make_sharpgrid_derivatives_work(
     &    itype,         icomp,         npwx,          npwy,
     &    npwz,          Lx,            Ly,            Lz,
     &    cftwork,       ncontr,        icfpl,         plbasdat,
     &    ncfpl,         xp,            yp,            zp,
     &    ftwork,        expfuncx,      expfuncy,      expfuncz,
     &    cftwork2,      nzeffdim,      ftwork2)
c
c This subroutine calculates the (exact) projection of the first
c derivatives (derivatives with respect to Rx,Ry,Rz which are the
c center of the given function) of the given compact contracted basis
c function to the given Fourier grid space (coordinate space) using
c analytical FT followed by numerical back FT.
c
c Inputs:
c          itype              type of the shell
c          icomp              the given component of the shell
c          npwx,npwy,npwz     Number of plane-wave in each direction
c          Lx,Ly,Lz           Box lengths for each direction
c          ncontr             number of primitive shells in the given
c                             contraction
c          plbasdat           pw basdat: see before
c          ncfpl              number of primitive core-like shells
c          xp,yp,zp           coordinates of the origin of the
c                             actual box
c          expfuncx,
c          expfuncy,          Precalculated exp functions
c          expfuncz
c          cftwork,           arrays for FT
c          ftwork2,cftwork2
c
c InOut:
c          icfpl              the starting primitive shell number
c
c Output:
c         ftwork              array including the projection of the
c                             three derivatives  (all together) of
c                             the given basis function in coordinate
c                             space
c                            ftwork(1,iz,iy,ix) -> d/dRx of the function
c                            ftwork(2,iz,iy,ix) -> d/dRy of the function
c                            ftwork(3,iz,iy,ix) -> d/dRz of the function
c **WARNING**
c   arrays cftwork and ftwork2 share storage        ! JB Feb 2006
c

      use memory

      implicit none
c
      integer itype,icomp,npwx,npwy,npwz,ix,iy,iz,ndim(3),plan
      integer howmany, istride, idist, iout, ostride, odist
      integer irec,isharpgrd,ncontr,icfpl,ncfpl,icontr,icfpl2,iconcoefs
      integer iexponents
      integer ifcount,nymax,nzmax,nzeffdim,nzeffdim2
      integer ncspl,iodd
      real*8 vmax,value
      real*8 exponent,concoef,cji,Px,Py,Pz,Lx,Ly,Lz,cc,xp,yp,zp,summ
      complex*16 cftwork(3,nzeffdim,npwy,npwx),czero
      real*8 ftwork2(6,nzeffdim,npwy,npwx)
      real*8 ftwork(3,npwz,npwy,npwx)
      complex*16 cftwork2(nzeffdim,npwy,npwx)
      real*8 plbasdat(6,ncfpl)
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 t1,t2,tfftb
      save tfftb
      integer ifw1,ifw2,i
c     common /big/bl(5000000)
      real*8 pi
      data pi/3.1415926535897932384626433d0/
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' On entry to <make_sharpgrid_derivatives_work>'
cc      write(6,*) ' itype:',itype,' icomp:',icomp,' ncontr:',ncontr
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' Lx:',lx,' Ly:',ly,' Lz:',lz
cc      write(6,*) ' icfpl:',icfpl,' ncfpl:',ncfpl
cc      write(6,*) ' xp:',xp,' yp:',yp,' zp:',zp
cc      write(6,*) ' expfuncx array is:'
cc      do i=1,npwx
cc      write(6,*) i,'  ',expfuncx(i)
cc      enddo
cc      write(6,*) ' expfuncy array is:'
cc      do i=1,npwy
cc      write(6,*) i,'  ',expfuncy(i)
cc      enddo
cc      write(6,*) ' expfuncz array is:'
cc      do i=1,npwz
cc      write(6,*) i,'  ',expfuncz(i)
cc      enddo
cc      call f_lush(6)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      czero=dcmplx(0.0d0,0.0d0)
      cji=1.0d0
      nzeffdim2=nzeffdim-1 !this is npwz/2 for both odd and even cases
c
      call mmark
      call getmem(ncontr,iconcoefs)
      call getmem(ncontr,iexponents)
c
      icfpl2=icfpl
      do icontr=1,ncontr
        icfpl2=icfpl2+1
        concoef=plbasdat(2,icfpl2)
        if((itype .eq. 3) .and. (icomp .gt. 1)) then
          concoef=plbasdat(6,icfpl2)
        end if
        exponent=plbasdat(1,icfpl2)
c
        bl(iexponents+icontr-1)=1.0d0/exponent !1/eta
c
        if(itype .eq. 1) then
c
          bl(iconcoefs+icontr-1)=
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                         !(2pi/eta)^(3/4)*concoef
c
        else if((itype .eq. 3) .and. (icomp .eq. 1)) then !s part of L
c
          bl(iconcoefs+icontr-1)=
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                          !(2pi/eta)^(3/4)*concoef
c
        else if((itype .eq. 2) .or. (itype .eq. 3)) then
                                              !P or P part of L
c
          bl(iconcoefs+icontr-1)=0.5d0 / exponent *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                      !0.5*(2pi/eta)^(3/4)*concoef/eta
c
c   *** D functions ***
c
        else if((itype .eq. 4) .and. (icomp .eq. 1)) then !D1
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef /sqrt(12.0d0)
                           !0.25*(2pi/eta)^(3/4)*concoef/eta^2/sqrt(12)
c
        else if((itype .eq. 4) .and. (icomp .eq. 2)) then !D2
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef /2.0d0
                                  !0.25*(2pi/eta)^(3/4)*concoef/eta^2/2
c
        else if(itype .eq. 4) then !D3,D4,D5
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                    !0.25*(2pi/eta)^(3/4)*concoef/eta^2
c
c   *** D6 functions ***
c
        else if(itype .eq. 5) then !D6
c
          bl(iconcoefs+icontr-1)=0.25d0 / exponent**2 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                    !0.25*(2pi/eta)^(3/4)*concoef/eta^2
c
c   ***  F function  ***
c
        else if((itype .eq. 6) .and. (icomp .eq. 7)) then ! 7-th component of f
c
          bl(iconcoefs+icontr-1)=0.125d0 * sqrt(40.0d0) / exponent**3 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                          !0.125*sqrt(40)*(2pi/eta)^(3/4)*concoef/eta^3
c
        else if(itype .eq. 6) then ! other f functions
c
          bl(iconcoefs+icontr-1)=0.125d0 / exponent**3 *
     &    (sqrt(sqrt((2.0d0*pi/exponent)))**3)*concoef
                                   !0.125*(2pi/eta)^(3/4)*concoef/eta^3
c
        end if
      end do
c
      Px=plbasdat(3,icfpl+1)-xp
      Py=plbasdat(4,icfpl+1)-yp
      Pz=plbasdat(5,icfpl+1)-zp
c
      call anal_sharp_derivatives_ongrid(
     &    itype,         bl(iexponents),bl(iconcoefs), icomp,
     &    Px,            Py,            Pz,            npwx,
     &    npwy,          npwz,          Lx,            Ly,
     &    Lz,            cftwork,       ncontr,        expfuncx,
     &    expfuncy,      expfuncz,      nzeffdim,      cftwork2)
c
      call retmark
      icfpl=icfpl+ncontr
c
      call getmem(npwx*npwy*npwz,ifw1)
      call getmem(2*npwx*npwy,ifw2)
      call backfft_grad(cftwork,bl(ifw2),bl(ifw1),
     &                  npwx,npwy,npwz,nzeffdim)
      call retmem(2)
c
      cc=1.0d0/Lx**2/Ly**2/Lz**2
      cc=sqrt(cc)
c
      if((npwz/2)*2 .ne. npwz) then ! npwz is odd
        iodd=1
      else
        iodd=0
      end if
c
      if(iodd .eq. 0) then
c
        do ix=1,npwx
          do iy=1,npwy
            do iz=1,nzeffdim2
              ftwork(1,2*iz-1,iy,ix)=ftwork2(1,iz,iy,ix)*cc
              ftwork(2,2*iz-1,iy,ix)=ftwork2(2,iz,iy,ix)*cc
              ftwork(3,2*iz-1,iy,ix)=ftwork2(3,iz,iy,ix)*cc
              ftwork(1,2*iz,iy,ix)=ftwork2(4,iz,iy,ix)*cc
              ftwork(2,2*iz,iy,ix)=ftwork2(5,iz,iy,ix)*cc
              ftwork(3,2*iz,iy,ix)=ftwork2(6,iz,iy,ix)*cc
            end do
          end do
        end do
c
      else
c
        do ix=1,npwx
          do iy=1,npwy
            do iz=1,nzeffdim2
              ftwork(1,2*iz-1,iy,ix)=ftwork2(1,iz,iy,ix)*cc
              ftwork(2,2*iz-1,iy,ix)=ftwork2(2,iz,iy,ix)*cc
              ftwork(3,2*iz-1,iy,ix)=ftwork2(3,iz,iy,ix)*cc
              ftwork(1,2*iz,iy,ix)=ftwork2(4,iz,iy,ix)*cc
              ftwork(2,2*iz,iy,ix)=ftwork2(5,iz,iy,ix)*cc
              ftwork(3,2*iz,iy,ix)=ftwork2(6,iz,iy,ix)*cc
            end do
           ftwork(1,2*nzeffdim2+1,iy,ix)=ftwork2(1,nzeffdim2+1,iy,ix)*cc
           ftwork(2,2*nzeffdim2+1,iy,ix)=ftwork2(2,nzeffdim2+1,iy,ix)*cc
           ftwork(3,2*nzeffdim2+1,iy,ix)=ftwork2(3,nzeffdim2+1,iy,ix)*cc
          end do
        end do
c
      end if
c
      return
      end
c
c***********************************************************************
c
      subroutine backfft_grad(cftw,cftw2,fw,npwx,npwy,npwz,nzeffdim)
c
c backward fft for gradient calculation
c

      use memory

      implicit none
      integer npwx,npwy,npwz,nzeffdim
      complex*16 cftw(3,nzeffdim,npwy,npwx)
      real*8 cftw2(2,npwx,npwy)
      real*8 fw(npwx,npwy,npwz)
      real*8 val1,val2,val3,val4,val5,val6
c     common /big/bl(1)
      integer iwx,iwy,iwz
      integer ic,ix,iy,iz
c
c  allocate memory for work arrays for fft
c
      call mmark
      call getmem(2*npwx+15,iwx)
      call getmem(2*npwy+15,iwy)
      call getmem(npwz+15,iwz)
c
c  fft initialization
c
      call vzffti(npwx,bl(iwx))
      call vzffti(npwy,bl(iwy))
      call vrffti(npwz,bl(iwz))
      do ic=1,3 ! main loop over gradient components

        do iz=1,nzeffdim
          do ix=1,npwx
            do iy=1,npwy
              cftw2(1,ix,iy)=real(cftw(ic,iz,iy,ix))
              cftw2(2,ix,iy)=aimag(cftw(ic,iz,iy,ix))
            enddo
          enddo
c
c FTBACK along x
c
          do iy=1,npwy
            call zfftb(npwx,cftw2(1,1,iy),bl(iwx))
          enddo
c
c FTBACK along y
c
          call vzfftb(1,npwx,1,npwy,cftw2,npwx,bl(iwy))
          if(iz.eq.1)then
            do ix=1,npwx
              do iy=1,npwy
                fw(ix,iy,1)=cftw2(1,ix,iy)
              enddo
            enddo
          else if(mod(npwz,2).eq.0.and.iz.eq.npwz/2+1)then
            do ix=1,npwx
              do iy=1,npwy
                fw(ix,iy,npwz)=cftw2(1,ix,iy)
              enddo
            enddo
          else
            do ix=1,npwx
              do iy=1,npwy
                fw(ix,iy,2*iz-2)=cftw2(1,ix,iy)
                fw(ix,iy,2*iz-1)=cftw2(2,ix,iy)
              enddo
            enddo
          endif
        enddo
c
c FTBACK along z
c
        call vrfftb(1,npwx*npwy,1,npwz,fw,npwx*npwy,bl(iwz))
c
        do iz=1,nzeffdim-1
          do iy=1,npwy
            do ix=1,npwx
              cftw(ic,iz,iy,ix)=dcmplx(fw(ix,iy,2*iz-1),fw(ix,iy,2*iz))
            enddo
          enddo
        enddo
        if(nzeffdim.eq.npwz/2+1)then
          do iy=1,npwy
            do ix=1,npwx
             cftw(ic,nzeffdim,iy,ix)=dcmplx(fw(ix,iy,npwz),0.0d0)
            enddo
          enddo
        endif
      enddo
c
c   OK, now looks like i need to transpose some indexes, or
c   something like that. I have absolutely no clue of what I am
c   doing here and why this is necessary.
c
      do iz=1,nzeffdim
        do iy=1,npwy
          do ix=1,npwx
            val1=real(cftw(1,iz,iy,ix))
            val4=aimag(cftw(1,iz,iy,ix))
            val2=real(cftw(2,iz,iy,ix))
            val5=aimag(cftw(2,iz,iy,ix))
            val3=real(cftw(3,iz,iy,ix))
            val6=aimag(cftw(3,iz,iy,ix))
            cftw(1,iz,iy,ix)=dcmplx(val1,val2)
            cftw(2,iz,iy,ix)=dcmplx(val3,val4)
            cftw(3,iz,iy,ix)=dcmplx(val5,val6)
          enddo
        enddo
      enddo
c
      call retmark
      end
cccccccc
      subroutine print_ifac(ifac)
      integer ifac(15)
c
      write(6,*) ' ifac: ',(ifac(i),i=1,15)
c
      return
      end

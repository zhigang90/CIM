      subroutine fillin_rangeconst(rangeconst)
      IMPLICIT REAL*8(A-H,O-Z)
c
c This subroutine fill in my so called range const matrix.
c The values are distances in atomic unit. These approximate values
c mean what is the shortest distances where the values of a given shell
c (given type with eta=1 exponents) is already negligible.
c I has been determined these values for 8 different accuracy and for
c s,p,d,f shells separately.
c These values were determined with
c my mapple program using the following process (example for the s
c function):
c
c f=exp(-r^2)
c N1=int(r^2*f^2,r=0..infinity)
c N2=int(r^2*f^2,r=R..infinity)
c error=N2/N1.
c Here R is the shortest distance, which is the question,
c if the error is bellow then the given accuracy limit.
c
      real*8 rangeconst(4,8)
c
c accuracy limit (error) 10^-7
c
      rangeconst(1,1)=2.98d0
      rangeconst(2,1)=3.20d0
      rangeconst(3,1)=3.39d0
      rangeconst(4,1)=3.55d0
c
c accuracy limit 10^-8
c
      rangeconst(1,2)=3.17d0
      rangeconst(2,2)=3.39d0
      rangeconst(3,2)=3.57d0
      rangeconst(4,2)=3.73d0
c
c accuracy limit 10^-9
c
      rangeconst(1,3)=3.35d0
      rangeconst(2,3)=3.56d0
      rangeconst(3,3)=3.74d0
      rangeconst(4,3)=3.90d0
c
c accuracy limit 10^-10
c
      rangeconst(1,4)=3.52d0
      rangeconst(2,4)=3.73d0
      rangeconst(3,4)=3.91d0
      rangeconst(4,4)=4.06d0
c
c accuracy limit 10^-11
c
      rangeconst(1,5)=3.69d0
      rangeconst(2,5)=3.89d0
      rangeconst(3,5)=4.06d0
      rangeconst(4,5)=4.22d0
c
c accuracy limit 10^-12
c
      rangeconst(1,6)=3.84d0
      rangeconst(2,6)=4.04d0
      rangeconst(3,6)=4.21d0
      rangeconst(4,6)=4.36d0
c
c accuracy limit 10^-13
c
      rangeconst(1,7)=3.99d0
      rangeconst(2,7)=4.19d0
      rangeconst(3,7)=4.36d0
      rangeconst(4,7)=4.51d0
c
c accuracy limit 10^-14
c
      rangeconst(1,8)=4.14d0
      rangeconst(2,8)=4.33d0
      rangeconst(3,8)=4.50d0
      rangeconst(4,8)=4.64d0
c
      return
      end
c***********************************************************************
c
      subroutine fillin_expcuts(expcuts)
c
c This subroutine fill in the expcuts array. These numbers
c were determined the followings:
c 1. analytical FT of the given basis function
c 2. calculate the norm of the function in momentum space
c 3. determine  the limit of the exponent for the given
c    basis function, accuracy limit and grid density.
c  The first index belongs to the type of the basis function
c  (s=1,p=2,d=3,f=4), the second index belongs to the grid density
c  (d=2 -> 1, d=4 -> 2), and the third index belongs to the given
c  accuracy limit beginning from 10^-4 to 10^-13. Currently only d=4
c  is used in the program !
c ........................................................................
c IMPORTANT: modified so that S and P shells have same upper limit
c to prevent problems when L shell is segmented     ! JB Feb 2006
c ........................................................................
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 expcuts(4,2,10)
c
c accuracy limit 10^-4
c
cc      expcuts(1,1,1)=1.80d0
      expcuts(1,1,1)=1.52d0
      expcuts(2,1,1)=1.52d0
      expcuts(3,1,1)=1.24d0
      expcuts(4,1,1)=1.31d0
c
cc      expcuts(1,2,1)=7.18d0
      expcuts(1,2,1)=6.06d0
      expcuts(2,2,1)=6.06d0
      expcuts(3,2,1)=4.93d0
      expcuts(4,2,1)=5.23d0
c
c accuracy limit 10^-5
c
cc      expcuts(1,1,2)=1.49d0
      expcuts(1,1,2)=1.28d0
      expcuts(2,1,2)=1.28d0
      expcuts(3,1,2)=1.06d0
      expcuts(4,1,2)=1.08d0
c
cc      expcuts(1,2,2)=5.94d0
      expcuts(1,2,2)=5.09d0
      expcuts(2,2,2)=5.09d0
      expcuts(3,2,2)=4.24d0
      expcuts(4,2,2)=4.31d0
c
c accuracy limit 10^-6
c
cc      expcuts(1,1,3)=1.27d0
      expcuts(1,1,3)=1.10d0
      expcuts(2,1,3)=1.10d0
      expcuts(3,1,3)=0.94d0
      expcuts(4,1,3)=0.93d0
c
cc      expcuts(1,2,3)=5.06d0
      expcuts(1,2,3)=4.39d0
      expcuts(2,2,3)=4.39d0
      expcuts(3,2,3)=3.73d0
      expcuts(4,2,3)=3.71d0
c
c accuracy limit 10^-7
c
cc      expcuts(1,1,4)=1.11d0
      expcuts(1,1,4)=0.97d0
      expcuts(2,1,4)=0.97d0
      expcuts(3,1,4)=0.84d0
      expcuts(4,1,4)=0.82d0
c
cc      expcuts(1,2,4)=4.41d0
      expcuts(1,2,4)=3.87d0
      expcuts(2,2,4)=3.87d0
      expcuts(3,2,4)=3.33d0
      expcuts(4,2,4)=3.28d0
c
c accuracy limit 10^-8
c
cc      expcuts(1,1,5)=0.98d0
      expcuts(1,1,5)=0.87d0
      expcuts(2,1,5)=0.87d0
      expcuts(3,1,5)=0.76d0
      expcuts(4,1,5)=0.74d0
c
cc      expcuts(1,2,5)=3.91d0
      expcuts(1,2,5)=3.46d0
      expcuts(2,2,5)=3.46d0
      expcuts(3,2,5)=3.01d0
      expcuts(4,2,5)=2.95d0
c
c accuracy limit 10^-9
c
cc      expcuts(1,1,6)=0.88d0
      expcuts(1,1,6)=0.79d0
      expcuts(2,1,6)=0.79d0
      expcuts(3,1,6)=0.69d0
      expcuts(4,1,6)=0.68d0
c
cc      expcuts(1,2,6)=3.51d0
      expcuts(1,2,6)=3.13d0
      expcuts(2,2,6)=3.13d0
      expcuts(3,2,6)=2.75d0
      expcuts(4,2,6)=2.69d0
c
c accuracy limit 10^-10
c
cc      expcuts(1,1,7)=0.80d0
      expcuts(1,1,7)=0.72d0
      expcuts(2,1,7)=0.72d0
      expcuts(3,1,7)=0.64d0
      expcuts(4,1,7)=0.62d0
c
cc      expcuts(1,2,7)=3.19d0
      expcuts(1,2,7)=2.86d0
      expcuts(2,2,7)=2.86d0
      expcuts(3,2,7)=2.53d0
      expcuts(4,2,7)=2.47d0
c
c accuracy limit 10^-11
c
cc      expcuts(1,1,8)=0.73d0
      expcuts(1,1,8)=0.66d0
      expcuts(2,1,8)=0.66d0
      expcuts(3,1,8)=0.59d0
      expcuts(4,1,8)=0.58d0
c
cc      expcuts(1,2,8)=2.92d0
      expcuts(1,2,8)=2.63d0
      expcuts(2,2,8)=2.63d0
      expcuts(3,2,8)=2.35d0
      expcuts(4,2,8)=2.29d0
c
c accuracy limit 10^-12
c
cc      expcuts(1,1,9)=0.68d0
      expcuts(1,1,9)=0.61d0
      expcuts(2,1,9)=0.61d0
      expcuts(3,1,9)=0.55d0
      expcuts(4,1,9)=0.54d0
c
cc      expcuts(1,2,9)=2.69d0
      expcuts(1,2,9)=2.44d0
      expcuts(2,2,9)=2.44d0
      expcuts(3,2,9)=2.19d0
      expcuts(4,2,9)=2.13d0
c
c accuracy limit 10^-13
c
cc      expcuts(1,1,10)=0.63d0
      expcuts(1,1,10)=0.57d0
      expcuts(2,1,10)=0.57d0
      expcuts(3,1,10)=0.52d0
      expcuts(4,1,10)=0.50d0
c
cc      expcuts(1,2,10)=2.50d0
      expcuts(1,2,10)=2.26d0
      expcuts(2,2,10)=2.26d0
      expcuts(3,2,10)=2.05d0
      expcuts(4,2,10)=1.99d0
c
      return
      end
c
c***********************************************************************
c
      Subroutine determine_sharpdim(ncs,ncspl,ncfpl,basdat,inx,
     &       expcuts,sharpness,expacclimit)
      IMPLICIT REAL*8(A-H,O-Z)
c
c This small subroutine goes over the whole basis set and determines
c the sharpness of the given ics contracted shell and the number of the
c core-like contracted shells.
c
c sharpness=1 -> can be described with ro=2 grid density
c sharpness=2 -> can be described with ro=4 grid density
c sharpness=3 -> core like shells
c
c
c Input: ncs,basdat,inx  - as always
c     expacclimit - accuracy limit for the exponent cuts (possible 1-10)
c     expcuts - array for the exponent cutoffs
c
c Output:
c     ncspl  -  number of contracted core-like shells
c     ncfpl  -  number of primitive core-like shells
c     sharpness - array to give the sharpness for every ics
c                 (contracted shells)
c
c
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer sharpness(ncs)
      integer expacclimit
      real*8 expcuts(4,2,10)
c
      ncspl=0
      ncfpl=0
c
      do ics=1,ncs
c
        itype=inx(12,ics)
        ibeg=inx(1,ics)+1
        iend=inx(5,ics)
        isize=iend-ibeg+1
c
c determine the actual limits for the exponent
c
        if(itype .eq. 1) then
c s type
          explimit1=expcuts(1,1,expacclimit)
          explimit2=expcuts(1,2,expacclimit)
c
        else if((itype .eq. 2) .or. (itype .eq. 3)) then
c p and l type
          explimit1=expcuts(2,1,expacclimit)
          explimit2=expcuts(2,2,expacclimit)
c
        else if((itype .eq. 4) .or. (itype .eq. 5)) then
c d and d6 type
          explimit1=expcuts(3,1,expacclimit)
          explimit2=expcuts(3,2,expacclimit)
c
        else if((itype .eq. 6) .or. (itype .eq. 7)) then
c f and f10 type
          explimit1=expcuts(4,1,expacclimit)
          explimit2=expcuts(4,2,expacclimit)
c
        else
          call nerror(1,'determine_sharpdim',
     $    'unknown basis function',0,0)
        end if
c
        exponent=basdat(1,ibeg)
        if(exponent .gt. explimit2) then
c this function is a sharp one
          ncspl=ncspl+1
          ncfpl=ncfpl+isize
          sharpness(ics)=3
        else if(exponent .gt. explimit1) then
          sharpness(ics)=2
        else
          sharpness(ics)=1
        end if
c
      end do
C
      return
      end
c***********************************************************************
c
      Subroutine make_sharp_arrays(
     &    ncs,           ncfpl,         ncspl,         basdat,
     &    inx,           plbasdat,      icsplsize,     icspltype,
     &    sharpness,     cssharps,      cssharps2,     ncsdiff,
     &    icsdiff)
      IMPLICIT REAL*8(A-H,O-Z)
c
c This subroutine builds up some arrays for the core like gaussians for
c GAPW calculations.
c
c Input:
c         ncs,basdat,inx  - as always
c         sharpness - array for the sharpness for every ics
c                     (contracted shells)
c         ncspl  -  number of contracted core-like shells
c         ncfpl  -  number of primitive core-like shells
c         ncsdiff  -  number of contracted diffuse-like shells
c
c
c Outputs:
c         plbasdat  -  array for every primitive shells (icfpl)
c         plbasdat(1,icfpl) - exponent
c         plbasdat(2,icfpl) - contr. coef
c         plbasdat(3,icfpl) - Rx
c         plbasdat(4,icfpl) - Ry
c         plbasdat(5,icfpl) - Rz
c         plbasdat(6,icfpl) - contr. coef for the p part of l functions
c
c         icsplsize - integer array for the size of every
c                     contracted core-like shells
c
c         icspltype - integer array for the type of every contracted
c                     core-like shells (s=1,p=2,l=3,d=4 etc.)
c         cssharps - integer array to determine which ics are the sharps
c         icsdiff - integer array to determine which ics are the diffuse
c         cssharps2 - integer array: hanyadik sharp basisfugvenynel
c                     kezdodik az adott ics sharp
c
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 plbasdat(6,ncfpl)
      integer icsplsize(ncspl)
      integer icspltype(ncspl)
      integer sharpness(ncs)
      integer cssharps(ncspl)
      integer cssharps2(ncs)
      integer icsdiff(ncsdiff)
c
      Dimension NFunc(7)
      data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
c
      icfpl=0
      icspl=0
      icountcsf=1
      idiff=0
c
      do ics=1,ncs
        itype=inx(12,ics)
        ibeg=inx(1,ics)+1
        iend=inx(5,ics)
        isize=iend-ibeg+1
        exponent=basdat(1,ibeg)
        isharpness=sharpness(ics)
c check the shell whether it is sharp or not
        if(isharpness .eq. 3) then
          icspl=icspl+1
          icsplsize(icspl)=isize
          icspltype(icspl)=itype
          cssharps(icspl)=ics
          cssharps2(ics)=icountcsf
c
          ncomp = NFunc(itype)
c
          icountcsf=icountcsf+ncomp
c
          do ish=ibeg,iend
c
            icfpl=icfpl+1
            plbasdat(1,icfpl)=basdat(1,ish)
            plbasdat(2,icfpl)=basdat(2,ish)
            plbasdat(3,icfpl)=basdat(11,ish)
            plbasdat(4,icfpl)=basdat(12,ish)
            plbasdat(5,icfpl)=basdat(13,ish)
            plbasdat(6,icfpl)=basdat(3,ish)
c
          end do
        else
          idiff=idiff+1
          icsdiff(idiff)=ics
        end if
c
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine make_listsd(ncs, ncspl, cssharp, listsd)
      IMPLICIT INTEGER(A-Z)
c
c  This subroutine builds up the listsd integer array.
c  It sets 0 for diffuse shells and 1 for compact (sharp).
c
c  ARGUMENTS
c
c  ncs     -  number of shells
c  ncspl   -  number of core-like shells
c  cssharp -  integer array listing sharp shells
c  listsd  -  on exit, integer array listing all shells
c
c
      Dimension cssharp(ncspl),listsd(ncs)
c
c  Set all to 0
c
      call IZeroIT(listsd,ncs)
      do i=1,ncspl
        ii = cssharp(i)
        listsd(ii) = 1
      enddo
c
      return
      end
c***********************************************************************
c
      Subroutine calc_ranges(
     &    ncs,           basdat,        inx,           rangeconst,
     &    ranges,        Lxmin,         Lxmax,         Lymin,
     &    Lymax,         Lzmin,         Lzmax,         gridranges,
     &    irangeacc,     griddens,      ranges2,       gridranges2,
     &    irangeacc2,    npwx,          npwy,          npwz)
c
c This subroutine does 3 main things.
c 1. Determines a spherical range (r) for every contracted shell where
c    it is not negligible.
c 2. Determines the necessary size of the box along each dimension
c    for the whole molecule.
c 3. Determines the gridranges integer arrays for every cs using the
c    given size of the box.
c    These integer values in these arrays mean the first and the last
c    grid points (in every direction) where the value of the contr.
c    shell is not negligible.
c
c Input:
c     ncs, basdat,inx         as always
c     irangeacc,irangeacc2    range accuracy parameter for smooth-smooth
c                             and for mixed density. They must be equal
c                             in the current version !!!
c     rangeconst              already filled in mx with the range const
c                             elements; first dimension -> type(s,p,d,f)
c                             second dim -> accuracy class
c     griddens                grid density (Rydberg  is about
c                                           9.8696*griddens**2)
c
c Output:
c       ranges(ncs)           real*8 array for the spherical range for
c                             every ics
c       ranges2(ncs)          real*8 array for the spherical range for
c                             every ics (using different accuracy class
c                             for mixed density)
c       Lxmin,Lxmax,
c       Lymin,Lymax,          Box dimensions in au. (start,end)
c       Lzmin,Lzmax           (real*8 !!!)
c       gridranges(6,ncs)     Integers for the grid point ranges
c                             (ixmin,ixmax,iymin etc.) for every ics
c       gridranges2(6,ncs)     Same using  irangeacc2 accuracy class
c       npwx,npwy,npwz         Total number of grid points along x,y
c                              and z
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 rangeconst(4,8)
      real*8 ranges(ncs)
      real*8 ranges2(ncs)
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      integer gridranges(6,ncs)
      integer gridranges2(6,ncs)
c
      boxstepp=1.0d0/griddens
c
      Lxmin=0.0d0
      Lymin=0.0d0
      Lzmin=0.0d0
      Lxmax=0.0d0
      Lymax=0.0d0
      Lzmax=0.0d0
c
c
c Now determine the spherical range for every ics (contracted shells)
c
      do ics=1,ncs
        itype=inx(12,ics)
        iend=inx(5,ics)
c the last part of the contraction determines the range because it has
c the smallest exponent
        eta=basdat(1,iend)
        call range_determine2(itype,eta,irangeacc,rangeconst,r)
        ranges(ics)=r
c
        call range_determine2(itype,eta,irangeacc2,rangeconst,r2)
        ranges2(ics)=r2
c
        xmin=basdat(11,iend)-ranges(ics)
        xmax=basdat(11,iend)+ranges(ics)
        ymin=basdat(12,iend)-ranges(ics)
        ymax=basdat(12,iend)+ranges(ics)
        zmin=basdat(13,iend)-ranges(ics)
        zmax=basdat(13,iend)+ranges(ics)
c
        if(xmin .lt. Lxmin) then
          Lxmin=xmin
        end if
        if(ymin .lt. Lymin) then
          Lymin=ymin
        end if
        if(zmin .lt. Lzmin) then
          Lzmin=zmin
        end if
c
        if(xmax .gt. Lxmax) then
          Lxmax=xmax
        end if
        if(ymax .gt. Lymax) then
          Lymax=ymax
        end if
        if(zmax .gt. Lzmax) then
          Lzmax=zmax
        end if
c
      end do
c
c Now modify the Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax to have integer
c number of grid points along every region.
c
      if(Lxmin .le. 0.0d0) then
        np=int(Lxmin*griddens)-1 !np negative
        Lxmin=float(np)/griddens !Lxmin - 0 region has integer number
                                 !of grid points
      else
        np=int(Lxmin*griddens)+1 !np positive
        Lxmin=float(np)/griddens !0 - Lxmin region has integer number
                                 !of grid points
      end if
c
      if(Lxmax .le. 0.0d0) then
        np=int(Lxmax*griddens)-1 !np negative
        Lxmax=float(np)/griddens !Lxmax - 0 region has integer number
                                 !of grid points
      else
        np=int(Lxmax*griddens)+1 !np positive
        Lxmax=float(np)/griddens !0 - Lxmax region has integer number
                                 !of grid points
      end if
c
      if(Lymin .le. 0.0d0) then
        np=int(Lymin*griddens)-1 !np negative
        Lymin=float(np)/griddens !Lymin - 0 region has integer number
                                 !of grid points
      else
        np=int(Lymin*griddens)+1 !np positive
        Lymin=float(np)/griddens !0 - Lymin region has integer number
                                 !of grid points
      end if
c
      if(Lymax .le. 0.0d0) then
        np=int(Lymax*griddens)-1 !np negative
        Lymax=float(np)/griddens !Lymax - 0 region has integer number
                                 !of grid points
      else
        np=int(Lymax*griddens)+1 !np positive
        Lymax=float(np)/griddens !0 - Lymax region has integer number
                                 !of grid points
      end if
c
      if(Lzmin .le. 0.0d0) then
        np=int(Lzmin*griddens)-1 !np negative
        Lzmin=float(np)/griddens !Lzmin - 0 region has integer number
                                 !of grid points
      else
        np=int(Lzmin*griddens)+1 !np positive
        Lzmin=float(np)/griddens !0 - Lzmin region has integer number
                                 !of grid points
      end if
c
      if(Lzmax .le. 0.0d0) then
        np=int(Lzmax*griddens)-1 !np negative
        Lzmax=float(np)/griddens !Lzmax - 0 region has integer number
                                 !of grid points
      else
        np=int(Lzmax*griddens)+1 !np positive
        Lzmax=float(np)/griddens !0 - Lzmax region has integer number
                                 !of grid points
      end if
c
cc
 100  continue
c
c Now make the gridranges
c
        do ics=1,ncs
c
          npointadd=int(ranges(ics)*griddens)+1
          npointadd2=int(ranges2(ics)*griddens)+1
c
          iend=inx(5,ics)
c
          ixcenter=int((basdat(11,iend)-Lxmin)*griddens)
          gridranges(1,ics)=ixcenter-npointadd
          gridranges(2,ics)=gridranges(1,ics)+2*npointadd
          gridranges2(1,ics)=ixcenter-npointadd2
          gridranges2(2,ics)=gridranges2(1,ics)+2*npointadd2
c
          iycenter=int((basdat(12,iend)-Lymin)*griddens)
          gridranges(3,ics)=iycenter-npointadd
          gridranges(4,ics)=gridranges(3,ics)+2*npointadd
          gridranges2(3,ics)=iycenter-npointadd2
          gridranges2(4,ics)=gridranges2(3,ics)+2*npointadd2
c
          izcenter=int((basdat(13,iend)-Lzmin)*griddens)
          gridranges(5,ics)=izcenter-npointadd
          gridranges(6,ics)=gridranges(5,ics)+2*npointadd
          gridranges2(5,ics)=izcenter-npointadd2
          gridranges2(6,ics)=gridranges2(5,ics)+2*npointadd2
c
        end do
c
c checking
c
        do ics=1,ncs
c
          n=gridranges(1,ics)
          if(n .lt. 1) then
            Lxmin=Lxmin-boxstepp
            goto 100
          end if
c
          n=gridranges(3,ics)
          if(n .lt. 1) then
            Lymin=Lymin-boxstepp
            goto 100
          end if
c
          n=gridranges(5,ics)
          if(n .lt. 1) then
            Lzmin=Lzmin-boxstepp
            goto 100
          end if
c
          n=gridranges(2,ics)
          nn=int((Lxmax-Lxmin)*griddens-0.01d0)+1 !because the
                                                  !conversion is not
                                                  !exact !
          if(n .gt. nn) then
            Lxmax=Lxmax+boxstepp
            goto 100
          end if
c
          n=gridranges(4,ics)
          nn=int((Lymax-Lymin)*griddens-0.01d0)+1
          if(n .gt. nn) then
            Lymax=Lymax+boxstepp
            goto 100
          end if
c
          n=gridranges(6,ics)
          nn=int((Lzmax-Lzmin)*griddens-0.01d0)+1
          if(n .gt. nn) then
            Lzmax=Lzmax+boxstepp
            goto 100
          end if
c
        end do

c
c      write(6,*)
c      write(6,*)"Lxmin,Lxmax",Lxmin,Lxmax
c      write(6,*)"Lymin,Lymax",Lymin,Lymax
c      write(6,*)"Lzmin,Lzmax",Lzmin,Lzmax
c      write(6,*)
c for test
c      do ics=1,ncs
c      write(6,*) ics
c      write(6,*) gridranges(1,ics),gridranges(2,ics)
c      write(6,*) gridranges(3,ics),gridranges(4,ics)
c      write(6,*) gridranges(5,ics),gridranges(6,ics)
c     write(6,*)
c      end do
      npwx=int((Lxmax-Lxmin)*griddens-0.01d0)+1 !because the conversion
                                                ! is not exact !
      npwy=int((Lymax-Lymin)*griddens-0.01d0)+1
      npwz=int((Lzmax-Lzmin)*griddens-0.01d0)+1
c      write(6,*)"npwx,npwy,npwz:",npwx,npwy,npwz
c      write(6,*)"boxstepp=",boxstepp
c      stop
c
      return
      end
c***********************************************************************
c
      Subroutine calc_Dmax(ncs,basdat,inx,ranges,Dmax)
c
c This subroutine calculates the cutoff distance for the Coulomb
c operator (Dmax) from the spherical reagions of all csf.
c Inputs:
c   ncs,basdat,inx          as always
c   ranges                  shperical ranges for all ics
c Output:
c   Dmax                    cutoff distance for the Coulomb operator
c                           (maximum possible distance between two
c                            electrons in the molecule)
c
c
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 ranges(ncs)
c
      Dmax=0.0d0
      do ics=1,ncs
        iend=inx(5,ics)
        xi=basdat(11,iend)
        yi=basdat(12,iend)
        zi=basdat(13,iend)
        do jcs=1,ics
          jend=inx(5,jcs)
          xj=basdat(11,jend)
          yj=basdat(12,jend)
          zj=basdat(13,jend)
c
          D=(xi-xj)**2 + (yi-yj)**2 + (zi-zj)**2
          D=sqrt(D) + ranges(ics) + ranges(jcs)
          if(D .gt. Dmax) then
            Dmax=D
          end if
c
        end do
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine determine_expanded_dims(
     &    npwx,          npwy,          npwz,          npwxe,
     &    npwye,         npwze,         Dmax,          
     &                griddens,      Lxo,           Lyo,
     &    Lzo,           Lxe,           Lye,           Lze)
c
c This program determines the expanded box dimensions and box sizes
c along each dimension using small prime (2,3,5,7) factor
c decompositions for the number of grid points of the extended box
c (fftw efficient only in these cases).
c
c Inputs:
c  npwx,npwy,npwz             nr. of grid points in the orig. box
c  Lxo,Lyo,Lzo                orig box lengths
c  Dmax                       Coulomb cutoff distance
c  griddens                   grid density
c Outputs:
c  npwxe,npwye,npwze          nr. of grid points in the expanded box
c  Lxe,Lye,Lze                expanded box lengths
c Temporary working arrays:
c
c
      implicit real*8(a-h,o-z)
c
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
c
      npwadd=int(Dmax*griddens) + 1
c
 50   continue
        npwxe=npwx+npwadd
        if(npwxe/2*2 .ne. npwxe) then ! I would like to have an even
                                      ! number
          npwadd=npwadd+1
          npwxe=npwx+npwadd
        end if
        Lxe=Lxo + dfloat(npwadd)/griddens
        call prime(npwxe,n1)
        if(n1 .ne. npwxe) then
          npwadd=n1-npwx
          goto 50
        end if
c
      npwadd=int(Dmax*griddens) + 1
c
 60   continue
        npwye=npwy+npwadd
        if(npwye/2*2 .ne. npwye) then ! I would like to have an even
                                      ! number
          npwadd=npwadd+1
          npwye=npwy+npwadd
        end if
        Lye=Lyo + dfloat(npwadd)/griddens
        call prime(npwye,n1)
        if(n1 .ne. npwye) then
          npwadd=n1-npwy
          goto 60
        end if
c
      npwadd=int(Dmax*griddens) + 1
c
 70   continue
        npwze=npwz+npwadd
        if(npwze/2*2 .ne. npwze) then ! I would like to have an even
                                      ! number
          npwadd=npwadd+1
          npwze=npwz+npwadd
        end if
        Lze=Lzo + dfloat(npwadd)/griddens
        call prime(npwze,n1)
        if(n1 .ne. npwze) then
          npwadd=n1-npwz
          goto 70
        end if
c
c
c      write(6,*)
c      write(6,*)"Expanded dimensions:"
c      write(6,*)"npwxe,npwye,npwze:",npwxe,npwye,npwze
c      write(6,*)"Lxe:",Lxe
c      write(6,*)"Lye:",Lye
c      write(6,*)"Lze:",Lze
c      write(6,*)
c
c
      return
      end
c***********************************************************************
c
      Subroutine calc_overlapdim(ncs,ncspl,basdat,inx,ranges,
     &   cssharps,nsharpovs,max)
c
c This subroutine calculates what is the number of csh-s which
c are overlap with a given csh. It also calculates their maximum value.
c Input:
c      ncs,inx,basdat      As always
c      ncspl               Nr. of sharp contracted shells
c      cssharps            integer array to determine which ics
c                          are the sharps
c      ranges              Spherical ranges for every ics
c
c Output:
c      nsharpovs(ncspl)    integer array for the number of overlapping
c                          cs for every sharp type cs.
c      max                 their maximum value (for memory allocation
c                          later)
c
c

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer cssharps(ncspl)
      integer nsharpovs(ncspl)
      real*8 ranges(ncs)
      integer fsize,f1size
c
      Dimension NFunc(7)
      data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
      max=0
      nint=0
c
      do icspl=1,ncspl
        icount=0
        isharp=cssharps(icspl)
c
        jtype=inx(12,isharp)
        f1size = NFunc(jtype)
c
        isbeg=inx(1,isharp)+1
        isend=inx(5,isharp)
        jsize=isend-isbeg+1
        R1x=basdat(11,isend)
        R1y=basdat(12,isend)
        R1z=basdat(13,isend)
        range1=ranges(isharp)
c
        do ics=1,ncs
c
          itype=inx(12,ics)
          fsize = NFunc(itype)
c
          ibeg=inx(1,ics)+1
          iend=inx(5,ics)
          isize=iend-ibeg+1
c
          R2x=basdat(11,iend)
          R2y=basdat(12,iend)
          R2z=basdat(13,iend)
          range2=ranges(ics)
          if((R1x .eq. R2x) .and. (R1y .eq. R2y)
     &        .and. (R1z .eq. R2z)) then
            icount=icount+1
            nint=nint+fsize*f1size
          else
            d1=(R1x-R2x)**2 + (R1y-R2y)**2 + (R1z-R2z)**2
            distance=sqrt(d1)
            r12=range1+range2
            if(distance .lt. r12) then
              icount=icount+1
              nint=nint+fsize*f1size
            end if
          end if
        end do
c
        nsharpovs(icspl)=icount
c      write(6,*)"icspl,icount:",icspl,nsharpovs(icspl)
c
        if(icount .gt. max) max=icount

      end do
c
c      write(6,*)"npairs=",nint
      nint=nint*(nint+1)/2
c      write(6,*)"Nr. of traditional integrals=",nint
      store=float(nint)*4.0d0/1.0d9
c      write(6,*)"Storage in real*4 and in Gbyte is about",store
c
      return
      end
c***********************************************************************
c
      Subroutine calc_sharpovs(
     &    ncs,           ncspl,         basdat,        inx,
     &    ranges,        cssharps,      maxovs,        sharpovs,
     &    gridranges,    ifilesplit,    icspltype,     maxfilesize,
     &  isharpgridrange1,isharpgridrange2,icorecut)
c
c This subroutine calculates the matrix "sharpovs", which ics is
c overlapps with a given sharp type jcs. After that it determines
c the new gridregion for the sharps which is the union of the region
c of the all overlapping shells.
c
c This subroutine also determines which file will contain the
c sharp basisfunction belonging to a given sharp shell.
c
c Input:
c      ncs,basdat,inx     As always
c      ncspl              Nr. of sharp contracted shells
c      cssharps           integer array to determine which ics are
c                         the sharps
c      ranges             Spherical ranges for every ics
c      maxovs             Maximum possible number of overlapping
c                         shells with any ics
c      maxfilesize        max file size in word
c InOut:
c      gridranges         Cartesian grid point ranges for every ics.
c                         Modified for sharps.
c Output:
c sharpovs(maxovs,ncspl)  an integer matrix. The first element are
c                         the indexes of those ics whics are overlaps
c                         with the given sharp type jcs (second
c                         dimension).
c   ifilesplit(ncspl)     the given sharp shell is in the file 1 or 2
c                         or 3 etc.
c
c

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(mincorerange=15)
      parameter(mincorecut=12)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer cssharps(ncspl)
      integer ifilesplit(ncspl)
      real*8 ranges(ncs)
      integer sharpovs(maxovs,ncspl)
      integer gridranges(6,ncs)
      integer icspltype(ncspl)
      integer isharpgridrange1(6,ncspl)
      integer isharpgridrange2(6,ncspl)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! P6 polinom coefs
      Dimension NFunc(7)
      Data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
      str=dfloat(mincorecut)
      str3=str**3
      str4=str3*str
      str5=str4*str
      str6=str5*str
c
      c1r=-10.0d0/str6
      c2r=24.0d0/str5
      c3r=-15.0d0/str4
c
      c1l=-10.0d0/str6
      c2l=36.0d0/str5
      c3l=-45.0d0/str4
      c4l=20.0d0/str3

c
      ifilesizecount=0
      ifiles=1
c
      do icspl=1,ncspl
        icount=0
        isharp=cssharps(icspl)
        isend=inx(5,isharp)
        R1x=basdat(11,isend)
        R1y=basdat(12,isend)
        R1z=basdat(13,isend)
        range1=ranges(isharp)
c
        ixtotmin=gridranges(1,isharp) !these are the original grid ranges
        ixtotmax=gridranges(2,isharp)
        iytotmin=gridranges(3,isharp)
        iytotmax=gridranges(4,isharp)
        iztotmin=gridranges(5,isharp)
        iztotmax=gridranges(6,isharp)
c
        ixtotmino=ixtotmin !Save the original regions
        ixtotmaxo=ixtotmax
        iytotmino=iytotmin
        iytotmaxo=iytotmax
        iztotmino=iztotmin
        iztotmaxo=iztotmax
c
        do ics=1,ncs
          iend=inx(5,ics)
          R2x=basdat(11,iend)
          R2y=basdat(12,iend)
          R2z=basdat(13,iend)
          range2=ranges(ics)
          if((R1x .eq. R2x) .and. (R1y .eq. R2y)
     &        .and. (R1z .eq. R2z)) then
            icount=icount+1
            sharpovs(icount,icspl)=ics
c
            ixmin=gridranges(1,ics)
            ixmax=gridranges(2,ics)
            iymin=gridranges(3,ics)
            iymax=gridranges(4,ics)
            izmin=gridranges(5,ics)
            izmax=gridranges(6,ics)
            if(ixmin .lt. ixtotmin) ixtotmin=ixmin
            if(iymin .lt. iytotmin) iytotmin=iymin
            if(izmin .lt. iztotmin) iztotmin=izmin
            if(ixmax .gt. ixtotmax) ixtotmax=ixmax
            if(iymax .gt. iytotmax) iytotmax=iymax
            if(izmax .gt. iztotmax) iztotmax=izmax
c
          else
            d1=(R1x-R2x)**2 + (R1y-R2y)**2 + (R1z-R2z)**2
            distance=sqrt(d1)
            r12=range1+range2
            if(distance .lt. r12) then !there is overlap
              icount=icount+1
              sharpovs(icount,icspl)=ics
c
              ixmin=gridranges(1,ics)
              ixmax=gridranges(2,ics)
              iymin=gridranges(3,ics)
              iymax=gridranges(4,ics)
              izmin=gridranges(5,ics)
              izmax=gridranges(6,ics)
              if(ixmin .lt. ixtotmin) ixtotmin=ixmin
              if(iymin .lt. iytotmin) iytotmin=iymin
              if(izmin .lt. iztotmin) iztotmin=izmin
              if(ixmax .gt. ixtotmax) ixtotmax=ixmax
              if(iymax .gt. iytotmax) iytotmax=iymax
              if(izmax .gt. iztotmax) iztotmax=izmax
c
            end if
          end if
        end do
c
        if (icorecut .eq. 1) then
c
          irangeo=(ixtotmaxo-ixtotmino)/2
          irmin=max(mincorerange,irangeo)
          ixcenter=ixtotmino+irangeo
          ixcutmin1=max(ixcenter-irmin,ixtotmin)
          ixcutmin2=max(ixcenter-irmin-mincorecut,ixtotmin)
          ixcutmax1=min(ixcenter+irmin,ixtotmax)
          ixcutmax2=min(ixcenter+irmin+mincorecut,ixtotmax)
c
          iycenter=iytotmino+irangeo
          iycutmin1=max(iycenter-irmin,iytotmin)
          iycutmin2=max(iycenter-irmin-mincorecut,iytotmin)
          iycutmax1=min(iycenter+irmin,iytotmax)
          iycutmax2=min(iycenter+irmin+mincorecut,iytotmax)
c
          izcenter=iztotmino+irangeo
          izcutmin1=max(izcenter-irmin,iztotmin)
          izcutmin2=max(izcenter-irmin-mincorecut,iztotmin)
          izcutmax1=min(izcenter+irmin,iztotmax)
          izcutmax2=min(izcenter+irmin+mincorecut,iztotmax)
c
          isharpgridrange1(1,icspl)=ixcutmin1
          isharpgridrange1(2,icspl)=ixcutmax1
          isharpgridrange1(3,icspl)=iycutmin1
          isharpgridrange1(4,icspl)=iycutmax1
          isharpgridrange1(5,icspl)=izcutmin1
          isharpgridrange1(6,icspl)=izcutmax1
c
          isharpgridrange2(1,icspl)=ixcutmin2
          isharpgridrange2(2,icspl)=ixcutmax2
          isharpgridrange2(3,icspl)=iycutmin2
          isharpgridrange2(4,icspl)=iycutmax2
          isharpgridrange2(5,icspl)=izcutmin2
          isharpgridrange2(6,icspl)=izcutmax2
c
          gridranges(1,isharp)=ixtotmin ! These will be not needed later
          gridranges(2,isharp)=ixtotmax ! but the transition program
                                        ! needs them !
          gridranges(3,isharp)=iytotmin
          gridranges(4,isharp)=iytotmax
          gridranges(5,isharp)=iztotmin
          gridranges(6,isharp)=iztotmax
c
        else
          gridranges(1,isharp)=ixtotmin !these are the modified grid
          gridranges(2,isharp)=ixtotmax ! ranges for the delocalized
                                        ! sharps (no cut version)
          gridranges(3,isharp)=iytotmin
          gridranges(4,isharp)=iytotmax
          gridranges(5,isharp)=iztotmin
          gridranges(6,isharp)=iztotmax
        end if
c
c Now   for the file splitting staff
c
        itype=icspltype(icspl)
        ncomp = NFunc(itype)
c
        nsdim=(ixtotmax-ixtotmin+1) *(iytotmax-iytotmin+1) *
     &        (iztotmax-iztotmin+1) * ncomp
c
        ifilesizecount=ifilesizecount+nsdim
        if(ifilesizecount .gt. maxfilesize) then
          ifilesizecount=nsdim
          ifiles=ifiles+1
        end if
        ifilesplit(icspl)=ifiles
c        write(6,*)"icspl,ifilesplit:",icspl,ifilesplit(icspl)
c
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine calc_sharpovs_new(
     &    ncs,           ncspl,         basdat,        inx,
     &    ranges,        cssharps,      maxovs,        sharpovs,
     &    gridranges,    ifilesplit,    icspltype,     maxfilesize,
     &  isharpgridrange1,isharpgridrange2,icorecut,    iteration)
c
c This subroutine calculates the matrix "sharpovs", which ics is
c overlapps with a given sharp type jcs. After that it determines
c the new gridregion for the sharps which is the union of the region
c of the all overlapping shells.
c
c This subroutine also determines which file will contain the sharp
c basisfunction belonging to a given sharp shell.
c
c Input:
c      ncs,basdat,inx       As always
c      ncspl                Nr. of sharp contracted shells
c      cssharps             integer array to determine which ics are
c                           the sharps
c      ranges               Spherical ranges for every ics
c      maxovs               Maximum possible number of overlapping
c                           shells with any ics
c      maxfilesize          max file size in word
c InOut:
c      gridranges           Cartesian grid point ranges for every ics.
c                           Modified for sharps.
c Output:
c  sharpovs(maxovs,ncspl)   an integer matrix. The first element are
c                           the indexes of those ics whics are overlaps
c                           with the given sharp type jcs (second
c                           dimension).
c  ifilesplit(ncspl)        the given sharp shell is in the file 1
c                           or 2 or 3 etc.
c
c

      IMPLICIT REAL*8(A-H,O-Z)
      parameter(mincorerangei=10)
      parameter(mincorecuti=4)
      parameter(mincorerangef=16)
      parameter(mincorecutf=10)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer cssharps(ncspl)
      integer ifilesplit(ncspl)
      real*8 ranges(ncs)
      integer sharpovs(maxovs,ncspl)
      integer gridranges(6,ncs)
      integer icspltype(ncspl)
      integer isharpgridrange1(6,ncspl)
      integer isharpgridrange2(6,ncspl)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l !P6 polinom coefs
      Dimension NFunc(7)
      Data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
      if(iteration .eq. 0) then
        mincorerange=mincorerangei
        mincorecut=mincorecuti
      else
        mincorecut=mincorecutf
        mincorerange=mincorerangef
      end if
c
c        write(6,*)"Used mincorerange=",mincorerange
c        write(6,*)"Used mincorecut=",mincorecut
c        call f_lush(6)
c
      str=dfloat(mincorecut)
      str3=str**3
      str4=str3*str
      str5=str4*str
      str6=str5*str
c
      c1r=-10.0d0/str6
      c2r=24.0d0/str5
      c3r=-15.0d0/str4
c
      c1l=-10.0d0/str6
      c2l=36.0d0/str5
      c3l=-45.0d0/str4
      c4l=20.0d0/str3

c
      ifilesizecount=0
      ifiles=1
c
      do icspl=1,ncspl
        icount=0
        isharp=cssharps(icspl)
        isend=inx(5,isharp)
        R1x=basdat(11,isend)
        R1y=basdat(12,isend)
        R1z=basdat(13,isend)
        range1=ranges(isharp)
c
        ixtotmin=gridranges(1,isharp)!these are the original grid ranges
        ixtotmax=gridranges(2,isharp)
        iytotmin=gridranges(3,isharp)
        iytotmax=gridranges(4,isharp)
        iztotmin=gridranges(5,isharp)
        iztotmax=gridranges(6,isharp)
c
        ixtotmino=ixtotmin !Save the original regions
        ixtotmaxo=ixtotmax
        iytotmino=iytotmin
        iytotmaxo=iytotmax
        iztotmino=iztotmin
        iztotmaxo=iztotmax
c
        do ics=1,ncs
          iend=inx(5,ics)
          R2x=basdat(11,iend)
          R2y=basdat(12,iend)
          R2z=basdat(13,iend)
          range2=ranges(ics)
          if((R1x .eq. R2x) .and. (R1y .eq. R2y)
     &        .and. (R1z .eq. R2z)) then
            icount=icount+1
            sharpovs(icount,icspl)=ics
c
            ixmin=gridranges(1,ics)
            ixmax=gridranges(2,ics)
            iymin=gridranges(3,ics)
            iymax=gridranges(4,ics)
            izmin=gridranges(5,ics)
            izmax=gridranges(6,ics)
            if(ixmin .lt. ixtotmin) ixtotmin=ixmin
            if(iymin .lt. iytotmin) iytotmin=iymin
            if(izmin .lt. iztotmin) iztotmin=izmin
            if(ixmax .gt. ixtotmax) ixtotmax=ixmax
            if(iymax .gt. iytotmax) iytotmax=iymax
            if(izmax .gt. iztotmax) iztotmax=izmax
c
          else
            d1=(R1x-R2x)**2 + (R1y-R2y)**2 + (R1z-R2z)**2
            distance=sqrt(d1)
            r12=range1+range2
            if(distance .lt. r12) then !there is overlap
              icount=icount+1
              sharpovs(icount,icspl)=ics
c
              ixmin=gridranges(1,ics)
              ixmax=gridranges(2,ics)
              iymin=gridranges(3,ics)
              iymax=gridranges(4,ics)
              izmin=gridranges(5,ics)
              izmax=gridranges(6,ics)
              if(ixmin .lt. ixtotmin) ixtotmin=ixmin
              if(iymin .lt. iytotmin) iytotmin=iymin
              if(izmin .lt. iztotmin) iztotmin=izmin
              if(ixmax .gt. ixtotmax) ixtotmax=ixmax
              if(iymax .gt. iytotmax) iytotmax=iymax
              if(izmax .gt. iztotmax) iztotmax=izmax
c
            end if
          end if
        end do
c
        if (icorecut .eq. 1) then
c
          irangeo=(ixtotmaxo-ixtotmino)/2
          mincorerange=1.5d0*irangeo
c            write(6,*)"mincorerange,irangeo: ",mincorerange,irangeo
c            call f_lush(6)
          irmin=max(mincorerange,irangeo)
          ixcenter=ixtotmino+irangeo
          ixcutmin1=max(ixcenter-irmin,ixtotmin)
          ixcutmin2=max(ixcenter-irmin-mincorecut,ixtotmin)
          ixcutmax1=min(ixcenter+irmin,ixtotmax)
          ixcutmax2=min(ixcenter+irmin+mincorecut,ixtotmax)
c
          iycenter=iytotmino+irangeo
          iycutmin1=max(iycenter-irmin,iytotmin)
          iycutmin2=max(iycenter-irmin-mincorecut,iytotmin)
          iycutmax1=min(iycenter+irmin,iytotmax)
          iycutmax2=min(iycenter+irmin+mincorecut,iytotmax)
c
          izcenter=iztotmino+irangeo
          izcutmin1=max(izcenter-irmin,iztotmin)
          izcutmin2=max(izcenter-irmin-mincorecut,iztotmin)
          izcutmax1=min(izcenter+irmin,iztotmax)
          izcutmax2=min(izcenter+irmin+mincorecut,iztotmax)
c
          isharpgridrange1(1,icspl)=ixcutmin1
          isharpgridrange1(2,icspl)=ixcutmax1
          isharpgridrange1(3,icspl)=iycutmin1
          isharpgridrange1(4,icspl)=iycutmax1
          isharpgridrange1(5,icspl)=izcutmin1
          isharpgridrange1(6,icspl)=izcutmax1
c
          isharpgridrange2(1,icspl)=ixcutmin2
          isharpgridrange2(2,icspl)=ixcutmax2
          isharpgridrange2(3,icspl)=iycutmin2
          isharpgridrange2(4,icspl)=iycutmax2
          isharpgridrange2(5,icspl)=izcutmin2
          isharpgridrange2(6,icspl)=izcutmax2
c
          gridranges(1,isharp)=ixtotmin ! These will be not needed later
          gridranges(2,isharp)=ixtotmax ! but the transition program
                                        ! needs them !
          gridranges(3,isharp)=iytotmin
          gridranges(4,isharp)=iytotmax
          gridranges(5,isharp)=iztotmin
          gridranges(6,isharp)=iztotmax
c
        else
          gridranges(1,isharp)=ixtotmin !these are the modified grid
          gridranges(2,isharp)=ixtotmax ! ranges for the delocalized
                                        ! sharps (no cut version)
          gridranges(3,isharp)=iytotmin
          gridranges(4,isharp)=iytotmax
          gridranges(5,isharp)=iztotmin
          gridranges(6,isharp)=iztotmax
        end if
c
c Now   for the file splitting staff
c
        itype=icspltype(icspl)
        ncomp = NFunc(itype)
c
        nsdim=(ixtotmax-ixtotmin+1) *(iytotmax-iytotmin+1) *
     &        (iztotmax-iztotmin+1) * ncomp
c
        ifilesizecount=ifilesizecount+nsdim
        if(ifilesizecount .gt. maxfilesize) then
          ifilesizecount=nsdim
          ifiles=ifiles+1
        end if
        ifilesplit(icspl)=ifiles
c      write(6,*)"icspl,ifilesplit:",icspl,ifilesplit(icspl)
c
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine make_sharpgrids(
     &    ncs,           ncspl,         cssharps,      gridranges,
     &    plbasdat,      ncfpl,         nfsh,          icsplsize,
     &    icspltype,     Lxmin,         Lxmax,         Lymin,
     &    Lymax,         Lzmin,         Lzmax,         griddens,
     &    isharpgrd,     comprcoefs,  isharpgridrange2,icorecut)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c
c This subroutine goes over the sharp shells makes some necessary
c arrangements and calls  the make_sharpgrids_work subroutine.
c That subroutine calculates the (exact) projection of the given
c compact basis function into the given Fourier grid space
c (coordinate space) using analytical FT followed by
c numerical back FT and writes the results onto the disk.
c
c Inputs:
c   ncs                number of contracted shells
c   ncspl              number of compact contracted shells
c   cssharps           integer array to determine which ics are
c                      the sharps
c   gridranges         grid ranges for every contracted shell
c   plbasdat           pw basdat, see before
c   ncfpl              number of primitive core-like shells
c   nfsh               number of contracted sharp functions
c   icsplsize          integer array for the size of every contracted
c                      core-like shells
c
c   icspltype          integer array for the type of every contracted
c                      core-like shells (s=1,p=2,l=3,d=4 etc.)
c  Lxmin,Lxmax,Lymin,  Box dimensions
c  Lymax,Lzmin,Lzmax
c  griddens            grid density
c  isharpgrd           unit number
c
      integer cssharps(ncspl)
      integer gridranges(6,ncs)
      real*8 Lx,Ly,Lz
      real*8 plbasdat(6,ncfpl)
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      integer icsplsize(ncspl)
      integer icspltype(ncspl)
      real*8 comprcoefs(*)
      integer isharpgridrange2(6,ncspl)
c     common /big/bl(300)
c
      Dimension NFunc(7)
      data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
c
      dx=1.0d0/griddens
      icfpl=0
      ifcount=0
c
      do icspl=1,ncspl
c
        ics=cssharps(icspl)
        itype=icspltype(icspl)
        ncontr=icsplsize(icspl)
c
        if (icorecut .eq. 1) then
c
          ixdmin=gridranges(1,ics)
          ixdmax=gridranges(2,ics)
          iydmin=gridranges(3,ics)
          iydmax=gridranges(4,ics)
          izdmin=gridranges(5,ics)
          izdmax=gridranges(6,ics)
c
          ixcmin=isharpgridrange2(1,icspl)
          ixcmax=isharpgridrange2(2,icspl)
          iycmin=isharpgridrange2(3,icspl)
          iycmax=isharpgridrange2(4,icspl)
          izcmin=isharpgridrange2(5,icspl)
          izcmax=isharpgridrange2(6,icspl)
c
          ixmin=max(ixdmin,ixcmin)
          iymin=max(iydmin,iycmin)
          izmin=max(izdmin,izcmin)
          ixmax=min(ixdmax,ixcmax)
          iymax=min(iydmax,iycmax)
          izmax=min(izdmax,izcmax)
c
          npwx=ixmax-ixmin+1
          npwy=iymax-iymin+1
          npwz=izmax-izmin+1
          Lx=float(npwx)/griddens
          Ly=float(npwy)/griddens
          Lz=float(npwz)/griddens
c
c determine the coordinates of the origin of the actual box
c
          xp=float(ixmin-1)*dx + Lxmin + Lx/2.0d0
          yp=float(iymin-1)*dx + Lymin + Ly/2.0d0
          zp=float(izmin-1)*dx + Lzmin + Lz/2.0d0
c
        else
c
          npwx=gridranges(2,ics)-gridranges(1,ics)+1
          npwy=gridranges(4,ics)-gridranges(3,ics)+1
          npwz=gridranges(6,ics)-gridranges(5,ics)+1
          Lx=float(npwx)/griddens
          Ly=float(npwy)/griddens
          Lz=float(npwz)/griddens
c
c determine the coordinates of the origin of the actual box
c
          xp=float(gridranges(1,ics)-1)*dx + Lxmin + Lx/2.0d0
          yp=float(gridranges(3,ics)-1)*dx + Lymin + Ly/2.0d0
          zp=float(gridranges(5,ics)-1)*dx + Lzmin + Lz/2.0d0
        end if
c
        ncomp = NFunc(itype)
c
        call mmark
        call getmem(npwx,iexpfuncx)
        call getmem(npwy,iexpfuncy)
        call getmem(npwz,iexpfuncz)
c
c Now   precalculate the exp functions to speed up the analytical FT
c
        call Precalc_expfuncs(
     &       bl(iexpfuncx),bl(iexpfuncy),bl(iexpfuncz),npwx,npwy,
     &       npwz,Lx,Ly,Lz)
c
c  allocate some scratch arrays
c
        nzeffdim=npwz/2+1 !for complex to real back Fourier trafo
        call getmem(2*npwx*npwy*nzeffdim,icftwork)
        call getmem(2*npwx*npwy*npwz,icftwork2)
        call getint_4(npwy*npwz,iftwork)
        call getmem(npwx*npwy*npwz,ifbk)
        call getmem(npwz+15,iw1)
        call getmem(2*npwy+15,iw3)
        call getmem(2*npwx+15,iw4)
c
        do icomp=1,ncomp
c
          call make_sharpgrids_work(
     &    itype,         icomp,         npwx,          npwy,
     &    npwz,          Lx,            Ly,            Lz,
     &    bl(icftwork),  isharpgrd,     ncontr,        icfpl,
     &    plbasdat,      ncfpl,         xp,            yp,
     &    zp,            bl(iftwork),   bl(iexpfuncx), bl(iexpfuncy),
     &    bl(iexpfuncz), bl(icftwork2), comprcoefs,    ifcount,
     &    nzeffdim,      bl(ifbk),      bl(iw1),       bl(iw3),
     &    bl(iw4))
c
          if(icomp .ne. ncomp) icfpl=icfpl-ncontr
c
        end do
c
        call retmark
c
      end do
c .....................................................................
c -- save the comprcoefs array to disk 
cc      WRITE(isharpgrd+1) ifcount
cc      CALL WriteBinary1(isharpgrd+1,ifcount,comprcoefs)
C
      return
      end
c***********************************************************************
c
      Subroutine readin_mocoefs(ncf,nmo,mocoef)
c
c This subroutine reads in the mocoefs into a mx from a file.
c
c Input:
c      ncf,nmo     number of contracted functions and occupied orbitals
c output:
c   mocoef         mx with the mo coefs
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 mocoef(ncf,nmo)
      character*256 jobname,MOS
      common /job/jobname,lenJ
c
c
c read in Mo coefs from MOS file
c
c      write(6,*)"ncf,nmo:",ncf,nmo
c
      call tstchval('mos-file',iyes)
      If(iyes.eq.1) Then
        call getchval('mos-file',MOS)
      Else
        MOS = jobname(1:lenJ)//'.mos'
      EndIf
c
      call rmblan2(MOS,256,lenM)
      itype=1
      call ReadMOS(ncf,mocoef(1,1),lenM,MOS,itype,IErr)
c
c
      return
      end
c***********************************************************************
c
      Subroutine calc_densmx(ncf,  nmo,  mocoef, dens)
c
c This subroutine builds up the density mx. using a simple
c (and expensive) way from the mo coefs.
c
c Input:
c   ncf,nmo       number of contracted basis functions and occ. orbitals
c   mocoef        mx for the mo coefs.
c Output
c   dens          density mx.
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 mocoef(ncf,nmo)
      real*8 dens(ncf,ncf)
c
c      write(6,*)"ncf,nmo:",ncf,nmo
c
c
      do icf=1,ncf
        do jcf=1,ncf
          do imo=1,nmo
            dens(jcf,icf)=dens(jcf,icf)+
     &        2.0d0*mocoef(icf,imo)*mocoef(jcf,imo)
          end do
c            write(6,*)"jcf,icf,dens(jcf,icf):",jcf,icf,dens(jcf,icf)
        end do
      end do
c
c
      return
      end
c***********************************************************************
c
      Subroutine makemapf2spw(inx,ncs,ncf,map_fs)
c
c This routine makes the mapping from a basis function to its shell.
c Input:
c   inx,ncs,ncf          as always
c Output:
c   map_fs       integer array: a given function belongs to which shell.
c
      implicit real*8(a-h,o-z)
      dimension inx(12,ncs)
      dimension map_fs(ncf)
c
      do  ics=1,ncs
        icf_b=inx(11,ics)+1
        icf_e=inx(10,ics)
        do icf=icf_b,icf_e
          map_fs(icf)=ics
        end do
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine range_determine2(itype,eta,irangeacc,rangeconst,r)
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 rangeconst(4,8)
      data pi/3.1415926535897932384626433d0/
c
c
      if (itype .eq. 1) then
c
c s type
c
        r=rangeconst(1,irangeacc)/sqrt(eta)
c
      else if((itype .eq. 2) .or. (itype .eq. 3)) then
c
c P or L type
c
        r=rangeconst(2,irangeacc)/sqrt(eta)
c
      else if((itype .eq. 4) .or. (itype .eq. 5)) then
c
c D or D6 type
c
        r=rangeconst(3,irangeacc)/sqrt(eta)
c
      else if((itype .eq. 6) .or. (itype .eq. 7)) then
c
c F or F10 type
c
        r=rangeconst(4,irangeacc)/sqrt(eta)
c
      end if
c
      return
      end
c***********************************************************************
c
      Subroutine prime(norig,n1)
c
c This subroutine performs the prime factor decomposition.
c Input:
c   norig    original grid points
c Output:
c   n1       closest bigger integer to norig having only 2,3,4,5 primes
c
c
c
      implicit real*8(a-h,o-z)
c
      integer iprimes(4)
      integer ipower(4)
c
      n=norig
c
      iprimes(1)=2
      iprimes(2)=3
      iprimes(3)=4
      iprimes(4)=5
c
      do i=1,4
       ipower(i)=0
      end do
c
 50   continue
      nold=n
      do i=1,4
       iprime=iprimes(i)
       ipower(i)=0
c
 100  continue
        if((n/iprime)*iprime .eq. n) then
         ipower(i)=ipower(i)+1
         n=n/iprime
         goto 100
        end if
c
      end do
c
      n1=2**ipower(1)*3**ipower(2)*4**ipower(3)*5**ipower(4)
      if(n1 .ne. nold) then
        n=nold+2
        do i=1,4
         ipower(i)=0
        end do
       goto 50
      else
c
      end if
c
      return
      end

      subroutine slave_sharpgrid(
     &   griddens,      Lxmin,
     &   Lymin,         Lzmin,         ncfpl,         isharpgrd,
     &   icorecut,      filname1,      len1,          nslave,
     &   ntot,          cssharps,      icspltype,     icsplsize,
     &   gridrange2,  isharpgridrange2,plbasdat,      comprcoefs)

      use memory
      use newpara


      implicit real*8(a-h,o-z)
c
c   slave code for driving the calculation of the sharp grid
c
c   in output the array nslave will keep track of which shells have
c   been computed
c
c     include 'fpvm3.h'
c     include 'txmsg.h'
      integer nslave(*),cssharps(*),icspltype(*),icsplsize(*)
      integer gridrange2(6,*),isharpgridrange2(6,*)
      real*8 plbasdat(6,*),comprcoefs(*)
      character*256 filename,filname1
      real*8 Lxmin,Lymin,Lzmin
      real*8 Lx,Ly,Lz
c     common /big/ bl(300)
      parameter (half=0.5d0)
      integer NFunc(7)
      data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
c  open file
c
      filename = filname1(1:len1)//'.sharpgrd'
      OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c -- now start main computation
c    we are going to keep track of which shells and the total number
c
      NTot = 0
      dx = 1.0d0/griddens
      ifcount = 0
c
      do
c
c -- get another index
c
        call para_send(MY_GID,0,TxDftReq)
        call para_recv_pack(ifrom,TxFTCAssign)
        call para_unpack_int(ish,1)
        call para_unpack_int(icfpl,1)
c
        IF(ish.eq.0) exit
c
c -- store and calculate current shell
c
        NTot = NTot+1
c
c
        NSLAVE(ntot)=ish
        ics=cssharps(ish)
        itype=icspltype(ish)
        ncontr=icsplsize(ish)
c
        IF(icorecut.eq.1) THEN
c
          ixdmin=gridrange2(1,ics)
          ixdmax=gridrange2(2,ics)
          iydmin=gridrange2(3,ics)
          iydmax=gridrange2(4,ics)
          izdmin=gridrange2(5,ics)
          izdmax=gridrange2(6,ics)
c
          ixcmin=isharpgridrange2(1,ish)
          ixcmax=isharpgridrange2(2,ish)
          iycmin=isharpgridrange2(3,ish)
          iycmax=isharpgridrange2(4,ish)
          izcmin=isharpgridrange2(5,ish)
          izcmax=isharpgridrange2(6,ish)
c
          ixmin=max(ixdmin,ixcmin)
          iymin=max(iydmin,iycmin)
          izmin=max(izdmin,izcmin)
          ixmax=min(ixdmax,ixcmax)
          iymax=min(iydmax,iycmax)
          izmax=min(izdmax,izcmax)
c
          npx=ixmax-ixmin+1
          npy=iymax-iymin+1
          npz=izmax-izmin+1
          Lx=float(npx)*dx
          Ly=float(npy)*dx
          Lz=float(npz)*dx
c
c determine the coordinates of the origin of the actual box
c
          xp=float(ixmin-1)*dx + Lxmin + Half*Lx
          yp=float(iymin-1)*dx + Lymin + Half*Ly
          zp=float(izmin-1)*dx + Lzmin + Half*Lz
        ELSE
          npx=gridrange2(2,ics)-gridrange2(1,ics)+1
          npy=gridrange2(4,ics)-gridrange2(3,ics)+1
          npz=gridrange2(6,ics)-gridrange2(5,ics)+1
          Lx=float(npx)*dx
          Ly=float(npy)*dx
          Lz=float(npz)*dx
c
c determine the coordinates of the origin of the actual box
c
          xp=float(gridrange2(1,ics)-1)*dx + Lxmin + half*Lx
          yp=float(gridrange2(3,ics)-1)*dx + Lymin + half*Ly
          zp=float(gridrange2(5,ics)-1)*dx + Lzmin + half*Lz
        ENDIF
c
        ncomp = NFunc(itype)
c
c -- precalculate the EXP function to speedup analytical FT
        call mmark
c
c   allocate arrays expfuncx(npx),expfuncy(npy),expfuncz(npz)
c
        call getmem(npx,iexpfuncx)
        call getmem(npy,iexpfuncy)
        call getmem(npz,iexpfuncz)
        call Precalc_expfuncs(
     &             bl(iexpfuncx),bl(iexpfuncy),bl(iexpfuncz),
     &             npx, npy, npz, Lx, Ly, Lz)
c
        nzeffdim=npz/2+1         !for complex to real back FT
c
c   allocate complex arrays cftwork(nzeffdim,npy,npx),
c   and cftwork2(npy,npx), integer array ftwork(npz,npy)
c   and arrays fbk(npx,npy,npz), w1(npz+15), w3(2*npy+15), w4(2*npx+15)
c
        call getmem(2*npx*npy*nzeffdim,icftwork)
        call getmem(2*npx*npy*npz,icftwork2)
        call getint_4(npy*npz,iftwork)
        call getmem(npx*npy*npz,ifbk)
        call getmem(npz+15,iw1)
        call getmem(2*npy+15,iw3)
        call getmem(2*npx+15,iw4)
c
        do icomp=1,ncomp
c
          call make_sharpgrids_work(
     &    itype,         icomp,        npx,            npy,
     &    npz,           Lx,           Ly,             Lz,
     &    bl(icftwork),  isharpgrd,    ncontr,         icfpl,
     &    plbasdat,      ncfpl,        xp,             yp,
     &    zp,            bl(iftwork),  bl(iexpfuncx),  bl(iexpfuncy),
     &    bl(iexpfuncz), bl(icftwork2),comprcoefs,     ifcount,
     &    nzeffdim,      bl(ifbk),     bl(iw1),        bl(iw3),
     &    bl(iw4))
c
          if(icomp .ne. ncomp) icfpl=icfpl-ncontr
c
        end do
c
        call retmark   ! deallocate expfuncx,expfuncy,expfuncz,
                       ! cftwork,ctfwork2,ftwork,fbk,w1,w3,w4
c
      enddo
c
      CLOSE (UNIT=isharpgrd)
c
      end
c=======================================================================
c
      subroutine slave_smoothd(
     &   ncf,           ncs,
     &   ncspl,         npwy,          npwz,          griddens,
     &   Lxmin,         Lymin,         Lzmin,         gridranges,
     &   ro1,           ro2,           sharpness,     da,
     &   inx,           bss,           iystoredim,    izstoredim,
     &   ipexpdim,      yexpstore,     zexpstore,     ipyexpstore,
     &   ipzexpstore)
c
c   slave code for driving the calculation of the smooth density
c

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
c     include 'fpvm3.h'
c     include 'txmsg.h'
      integer gridranges(*),sharpness(*),inx(*)
      integer ipyexpstore(*),ipzexpstore(*)
      real*8 ro1(*),ro2(*),da(*),bss(*),yexpstore(*),zexpstore(*)
      real*8 lxmin,lymin,lzmin
c     common /big/ bl(300)
c
      dx=1.0d0/griddens
c
c -- first get one block of data
c
      call para_recv_pack(ifrom,TxFTCAssign)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c
c -- now start main computation
c
  130 CONTINUE
c
c    allocate array ro4(npwz,npwy,ixstart:ixend)
c
      call getmem(npwz*npwy*(ixend-ixstart+1),iro4)
      call zeroit(bl(iro4),npwz*npwy*(ixend-ixstart+1))
c
c -- calculate current block
c
      ixcounter=0
      do ix=ixstart,ixend
c
      x=Lxmin+float(ix-1)*dx
      ixcounter=ixcounter+1
      iro4p=npwy*npwz*(ixcounter-1)
c
      call build_smooth_density(
     &   ncf,           ro1,           npwy,          npwz,
     &   gridranges,    ro2,           ncs,           sharpness,
     &   ncspl,         dA,            inx,           bss,
     &   x,             Lymin,         Lzmin,         bl(iro4+iro4p),
     &   ix,            1,             npwy,          griddens,
     &   iystoredim,    yexpstore,     izstoredim,    zexpstore,
     &   ipexpdim,      ipyexpstore,   ipzexpstore)
c
      end do
c
c  Send back to the master
c
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_pack_real(bl(iro4),npwy*npwz*(ixend-ixstart+1))
      call para_send_pack(0,TxFTCBTM)
c
c -- now receive new data from master
c
      call retmem(1)     !  deallocate ro4
c
      call para_recv_pack(ifrom,TxFTCAssign)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c -- test for termination
      If(ixstart.NE.0.AND.ixend.NE.0) GO TO 130
c
c  finished this section
c
      end
c=======================================================================
c
      subroutine slave_mixed_dpf(
     &   ncf,           ncs,
     &   ncspl,         npwx,          npwy,          npwz,
     &   ntot,          maxovs,        icorecut,      isharpgrd,
     &   filname1,      len1,          griddens,      Lxmin,
     &   Lymin,         Lzmin,         nslave,        icspltype,
     &   cssharps,      gridrange2,    inx,           bss,
     & isharpgridrange1,sharpness,     sharpovs,      nsharpovs,
     & isharpgridrange2,ro1,           ro2,           da,
     &   fa,            comprcoefs,    iystoredim,    izstoredim,
     &   ipexpdim,      yexpstore,     zexpstore,     ipyexpstore,
     &   ipzexpstore)
c
c  slave code for driving the calculation of mixed (smooth+sharp)
c  density and fock matrix
c

      use memory
      use newpara


      implicit real*8 (a-h,o-z)
c     include 'fpvm3.h'
c     include 'txmsg.h'
      integer gridrange2(*),sharpness(*),sharpovs(*),nsharpovs(*),inx(*)
      integer isharpgridrange1(*),isharpgridrange2(*),nslave(*)
      integer ipyexpstore(*),ipzexpstore(*),icspltype(*),cssharps(*)
      real*8 ro1(*),ro2(*),da(*),fa(*),bss(*)
      real*8 comprcoefs(*),yexpstore(*),zexpstore(*)
      real*8 Lxmin,Lymin,Lzmin
      character*256 filename,filname1
c     common /big/ bl(300)
c
      dx=1.0d0/griddens
c -- note: greatest memory usage here
c
c   allocate arrays ro4(npwz,npwy,npwx),ro5(npwz,npwy,npwx)
c
      call getmem(npwx*npwy*npwz,iro4)
      call getmem(npwx*npwy*npwz,iro5)
      call zeroit(bl(iro5),npwx*npwy*npwz)
c
c  zero out ro1 and ro2
c
      call zeroit(ro1,npwy*npwz)
      call zeroit(ro2,npwy*npwz)
c
c -- receive the full smooth Coulomb operator from the master
      call para_recv_bcastreal(bl(iro4),npwz*npwy*npwx,TxFTCDen)
c
c  open file
      filename = filname1(1:len1)//'.sharpgrd'
      OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
c
c -- calculate density for shells previously computed
c    we already know which shells we have data for
c
      xmin0 = Lxmin-dx
      indxcomprcoef=0
c
      DO ish=1,NTot
      ics = NSLAVE(ish)
c
      call mixed_density_help1(
     &   ics,           icspltype,     ncspl,         cssharps,
     &   ncomp,         gridrange2,    ncs,           isxmin,
     &   isxmax,        isymin,        isymax,        iszmin,
     &   iszmax,        inx,           ifirstf,       ilastf,
     &   icorecut,      isharpgridrange2)
c
      npx=isxmax-isxmin+1
      npy=isymax-isymin+1
      npz=iszmax-iszmin+1
      ymin=Lymin+float(isymin-1)*dx
      npwregy=isymax-isymin+1
c
c   allocate array ro62d(npz,npy)
c
      call getmem(npz*npy,iro62d)
c
      do ifunc=ifirstf,ilastf
      xmin=xmin0+(isxmin-1)*dx
c
      do ix=isxmin,isxmax
      indxcomprcoef=indxcomprcoef+1
      comprc = comprcoefs(indxcomprcoef)
c
      call readin_sharps2d(isharpgrd,npy,npz,bl(iro62d))
c
      xmin=xmin+dx
      iro4p=npwz*(isymin-1)+npwy*npwz*(ix-1)
c
      call make_mixed_dpf(
     &   ncf,           bl(iro62d),    npwregy,       npwz,
     &   gridrange2,    ro1,           ncs,           sharpness,
     &   ncspl,         nfsh,          igrdsh,        iro2,
     &   dA,            inx,           bss,           xmin,
     &   ymin,          Lzmin,         bl(iro5+iro4p),ix,
     &   isymin,        isymax,        ics,           ifunc,
     &   sharpovs,      maxovs,        nsharpovs,     npz,
     &   iszmin,        bl(iro4+iro4p),fA,            griddens,
     &   iystoredim,    yexpstore,     izstoredim,    zexpstore,
     &   ipexpdim,      ipyexpstore,   ipzexpstore,   comprc,
     &   icorecut,    isharpgridrange1,isharpgridrange2)
c
      end do
      end do
c
      call retmem(1) !  deallocate ro62d
c
      end do
c
      CLOSE(UNIT=isharpgrd)
c
c -- send back the mixed density to the master
        call para_reduce(bl(iro5),npwx*npwy*npwz,TxFTCBTM)
c
      call retmem(2)   ! deallocate ro5,ro4
      end
c=======================================================================
c
      subroutine slave_smooth_f(
     &   ncf,           ncs,           ncspl,
     &   npwx,          npwy,          npwz,          griddens,
     &   Lxmin,         Lymin,         Lzmin,         inx,
     &   bss,           ro1,           ro2,           ro3,
     &   fa,            gridranges,    sharpness,     iystoredim,
     &   izstoredim,    ipexpdim,      yexpstore,     zexpstore,
     &   ipyexpstore,   ipzexpstore)
c
c  slave code for driving the calculation of the smooth+mixed
c  fock matrix
c

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
c     include 'fpvm3.h'
c     include 'txmsg.h'
      integer gridranges(*),sharpness(*),inx(*)
      integer ipyexpstore(*),ipzexpstore(*)
      real*8 ro1(*),ro2(*),ro3(*),fa(*),bss(*)
      real*8 yexpstore(*),zexpstore(*)
      real*8 Lxmin,Lymin,Lzmin
c     common /big/ bl(300)
c
      dx=1.0d0/griddens
c
c  zero out ro1, ro2, ro3
c
      call zeroit(ro1,npwy*npwz)
      call zeroit(ro2,npwy*npwz)
      call zeroit(ro3,npwy*npwz)
c
c -- get block of data
c
      call para_recv_pack(ifrom,TxFTCAssign)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c
c -- now start main computation
c
  140 CONTINUE
c
c    allocate array ro4(npwz,npwy,ixstart:ixend)
c
      call getmem(npwy*npwz*(ixend-ixstart+1),iro4)
c
c -- get density block from master, separate b/ mem. allocation
c
      call para_recv_real(bl(iro4),npwy*npwz*(ixend-ixstart+1),
     &                    0,TxFTCAssign)
c
c -- calculate Fock matrix elements from current block
c
      ixcounter=0
      do ix=ixstart,ixend
c
        x=Lxmin+float(ix-1)*dx
        ixcounter=ixcounter+1
        iro4p=npwy*npwz*(ixcounter-1)
c
        call calc_FTC_fock_smooth(
     &   ncf,           ro1,           1,             npwy,
     &   npwz,          gridranges,    ro2,           ncs,
     &   sharpness,     ncspl,         igrdsh,
     &   ro3,           inx,           bss,           x,
     &   Lymin,         Lzmin,         bl(iro4+iro4p),ix,
     &   1,             npwy,          fA,            griddens,
     &   iystoredim,    yexpstore,     izstoredim,    zexpstore,
     &   ipexpdim,      ipyexpstore,   ipzexpstore)
c
      enddo
c
      call retmem(1) ! deallocate ro4
c
c -- now receive new data from master
c
      call para_send(MY_GID,0,TxDftReq)
      call para_recv_pack(ifrom,TxFTCAssign)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c
c -- test for termination
c
      If(ixstart.NE.0.AND.ixend.NE.0) GO TO 140
c
      end
c=======================================================================
c
      Subroutine Slave_COULOMB_on_FTGRID(
     &                 npwx,    npwy,    npwz,
     &                 npwxe,   npwye,   npwze,   PLDmax,  Lxe,
     &                 Lye,     Lze)
c
c  slave code for driving the calculation of the Coulomb potential
c

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
c     include 'fpvm3.h'
c     include 'txmsg.h'
      real*8 Lxe,Lye,Lze
      Logical firstcall
c     common /big/ bl(300)
c
c
c -- first allocate work arrays
      call mmark
      call getmem(2*npwye+15,iw3)
      call getmem(2*npwze+15,iw4)
      call getmem(npwye,ivky)
      call getmem(npwze,ivkz)
      call getmem(2*npwze*npwye,ivyze)
c
c -- initialize the FFT
      call vzffti(npwye,bl(iw3))
      call vzffti(npwze,bl(iw4))
      firstcall = .TRUE.

c -- get data from master
      call para_recv_pack(ifrom,TxFTCFFT)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c
 130  CONTINUE
c
      call para_unpack_int(lstart,1)
      call para_unpack_int(ndata,1)
cc      write(6,*) ' SLAVE:',MY_GID,' Receiving <TxFTCFFT>'
cc      write(6,*) ' SLAVE:',MY_GID,' ixstart:',ixstart,' ixend:',ixend
cc      write(6,*) ' SLAVE:',MY_GID,' lstart:',lstart,' ndata:',ndata
cc      call f_lush(6)
c
c -- allocate memory for bfk array
      call getmem(ndata,ibfk)
c
c -- received separately to avoid mem. allocation conflict      
      call para_recv_real(bl(ibfk),ndata,0,TxFTCFFT)
c
c -- compute ending index for bfk array
      lend   = MIN(npwxe,2*ixend-1)
c
      call Slave_FFT(npwx,    npwy,    npwz,    npwxe,   npwye,
     $               npwze,   PLDmax,  Lxe,     Lye,     Lze,
     $               ixstart, ixend,   lstart,  lend,    firstcall,
     $               bl(ibfk),bl(ivyze),bl(iw3),bl(iw4), bl(ivky),
     $               bl(ivkz))
cc      write(6,*) ' SLAVE:',MY_GID,' Back from <slave_fft>'
cc      call f_lush(6)
c
c -- send data back to master
      call para_initsend
      call para_pack_int(lstart,1)
      call para_pack_int(ndata,1)
cc      write(6,*) ' SLAVE:',MY_GID,' Sending <TxFTCBTM>'
cc      write(6,*) ' SLAVE:',MY_GID,' lstart:',lstart,' ndata:',ndata
cc      call f_lush(6)
      call para_pack_real(bl(ibfk),ndata)
      call para_send_pack(0,TxFTCBTM)
c
c -- release memory for bfk
      call retmem(1)
c
c -- now receive new data from master
      call para_recv_pack(ifrom,TxFTCFFT)
      call para_unpack_int(ixstart,1)
      call para_unpack_int(ixend,1)
c -- test for termination
      If(ixstart.NE.0.AND.ixend.NE.0) GO TO 130
c
c -- we're finished
c
      call retmark
c
      return
      end
c=======================================================================
c
      Subroutine Slave_FFT(
     &    npwxo,      npwyo,      npwzo,      npwxe,      npwye,
     &    npwze,      D,          Lxe,        Lye,        Lze,
     &    ixstart,    ixend,      lstart,     lend,     firstcall,
     &    bfk,        vyze,       w3,         w4,         vky,
     &    vkz)
c
c This subroutine calculates the exact Coulomb mx. elements on the
c Fourier grid using the density and the value of D which is
c the maximum distance in the Coulomb interaction.
c
c The routine uses the real and complex one dimensional
c FFT subroutines
c
c This subroutine is a modified version of <make_COULOMB_on_FTGRID>
c used in parallel on the slave nodes
c
c Input:
c     npwzo,npwyo,npwxo         - grid dimensions of the original box
c                                 ALL MUST BE EVEN NUMBER !!!
c     npwxe,npwxe,npwxe         - expanded grid dimensions
c                                 ALL MUST BE EVEN NUMBER !!!
c     D                         - max. possible  distance in the Coulomb
c                                 interaction
c     Lxe,Lye,Lze               - expanded box sizes
c     ixstart                   - starting index for FFT
c     ixend                     - ending index for FFT
c     lstart                    - starting index of bfk array
c     lend                      - ending index of bfk array
c     firstcall                 - set to true on first call to this
c                                 routine; thereafter false
c     bfk                       - partially FFTed density
c
c Working arrays:
c
c       bfk(npwzo,npwyo,lstart:lend),w3(*),w4(*),
c       vky(npwye),vkz(npwze)vyze(2,npwze,npwye)
c
c Output:
c      bfk                        - array with the Coulomb mx element
c                                   on the grid points multiplied by
c                                   npwxe*npwxe*npwxe
c
      implicit none
c
      integer npwxo,npwyo,npwzo,npwxe,npwye,npwze
      integer ixstart,ixend,iystart,iyend,izstart,izend
      integer lstart,lend
      integer ix,iy,iz,ixe,iye,ize,nx,ny,nz
      real*8 D,Lxe,Lye,Lze,cx,cy,cz
      logical firstcall
c
c arrays
c
      real*8 bfk(npwzo,npwyo,lstart:lend),w3(*),w4(*)
      real*8 vyze(2,npwze,npwye)
      real*8 vky(npwye),vkz(npwze)
c
      real*8 twopi2
      parameter(twopi2=3.9478417604357434475d1) ! (2*%pi)^2
      real*8 vkx,vkxy,valksq,cm,cmd2
      integer kx,jky,ky,jkz,kz
c
c   precompute some quantities for the Coulomb potential calculation
c   SHOULD TAKE SOME OF THIS STUFF OUT - JB
c
      cmd2 = 6.2831853071795864769d0 * D * D
      cx=twopi2/(Lxe*Lxe)
      cy=twopi2/(Lye*Lye)
      cz=twopi2/(Lze*Lze)
      nx=npwxe/2+2
      ny=npwye/2+2
      nz=npwze/2+2
c
      If(firstcall) Then
        do jky=1,npwye
          ky=(jky-1)-npwye*int(jky/ny)
          vky(jky)=cy*dfloat(ky*ky)
        enddo
        do jkz=1,npwze
          kz=(jkz-1)-npwze*int(jkz/nz)
          vkz(jkz)=cz*dfloat(kz*kz)
        enddo
        firstcall = .False.
      EndIf
c
      iystart=(npwye-npwyo)/2 + 1
      iyend=iystart+npwyo - 1
      izstart=(npwze-npwzo)/2 + 1
      izend=izstart+npwzo - 1
c
      DO ix=ixstart,ixend       ! main loop over kx
        kx=(ix-1)-npwxe*int(ix/nx)
        vkx=cx*dfloat(kx*kx)
c
c set to zero the elements of vyze belonging to the expanded grid only
c
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
c
c  now fill the elements belonging to both original and expanded grids
c
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
c
c FT along y for every npwzo sequence
c
        call vzfftf(izstart,izend,1,npwye,vyze,npwze,w3)
c
c FT along z
c
        do iy=1,npwye
          call zfftf(npwze,vyze(1,1,iy),w4)
        enddo
c
c Now we have everything in momentum space at this point
c Build up the Coulomb op. in momentum space
c
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
c
c FTBACK along z
c
        do iy=1,npwye
          call zfftb(npwze,vyze(1,1,iy),w4)
        enddo
c
c FTBACK along y
c
        call vzfftb(izstart,izend,1,npwye,vyze,npwze,w3)
c
c Now order back for complex_to_real trafo
c
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
c
      ENDDO    !loop over kx end
c
      return
      end
c=======================================================================
c
      Subroutine Slave_Multipoles_Main(
     &                 ncs,     inx,     basdat,
     &                 ncspl,  cssharps, natom,   ncf,     DA,
     &               sharpovs, maxovs, nsharpovs, dist2mp, sharpness,
     &                 XC,    multipolecut,  fskal0,  FA)
c
c  slave code for driving the Multipole contribution
c

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
c     include 'fpvm3.h'
c     include 'txmsg.h'
c     common /big/ bl(300)
c
      Dimension inx(12,ncs),basdat(13,*),XC(3,natom),
     $          DA(ncf*(ncf+1)/2),FA(ncf*(ncf+1)/2)
      Integer cssharps(ncspl),sharpovs(maxovs,ncspl),
     $        sharpness(ncs),nsharpovs(ncspl)
      Integer*1 multipolecut(natom,natom)
c
c
c --  allocate memory
      call mmark
      call getmem(ncf*ncf,idens)
      call getint(2*ncs,iishellinfo)
      call getmem(35*natom,iatomicpoles)
      call zeroit(bl(iatomicpoles),35*natom)
c
c -- expand out density
      call fillin_dens(bl(idens),DA,ncf)
c
      call find_ncfp_descartes(ncs,inx,ncfp,bl(iishellinfo))
      If(ncfp.GT.ncf) Then
c
        call getmem(ncfp*ncfp,idensp)
        call getmem(ncfp,iw1)
        call getmem(ncfp,iw2)
        call getmem(ncfp,iw3)
        call getmem(ncfp,iw4)
        call getmem(ncfp,iw5)
        call getmem(ncfp,iw6)
        call getmem(ncfp,iw7)
c
        call denstrafo_to_descartes(
     &    ncs,        inx,        ncf,        ncfp,       bl(idens),
     &    bl(idensp), bl(iw1),    bl(iw2),    bl(iw3),    bl(iw4),
     &    bl(iw5),    bl(iw6),    bl(iw7))
        iden = idensp
        nbas = ncfp
        call retmem(7)
      Else
        iden = idens
        nbas = ncf
      EndIf
c
c -- get data from master
      call para_recv_pack(ifrom,TxFTCMult)
      call para_unpack_int(istart,1)
      call para_unpack_int(iend,1)
c
c -- calculate atomicpoles for this range of indices
c
      call Slave_Collect_Multipoles(
     $               ncs,     inx,     basdat,  ncspl,   istart,
     $               iend,   cssharps, natom,   nbas,    bl(iden),
     $             sharpovs,  maxovs, nsharpovs,bl(iatomicpoles),
     $             sharpness, bl(iishellinfo))
c
c -- send data back to master
      call para_reduce(bl(iatomicpoles),35*natom,TxFTCBTM)
c
c -- allocate array v(8,natom)
      call getmem(8*natom,iv)
      call zeroit(bl(iv),8*natom)
c
c -- now collect fully formed atomicpoles array from master
      call para_recv_bcastreal(bl(iatomicpoles),35*natom,TxFTCDen)
c
c -- calculate Fock matrix contribution from selected indices
      call Slave_Multipoles_Fock(
     $               ncs,     inx,     basdat,  ncspl,   istart,
     $               iend,   cssharps, natom,   ncf,     XC,
     $             sharpovs,  maxovs, nsharpovs,bl(iatomicpoles),
     $             sharpness, bl(iv),multipolecut, fskal0,  FA)
c
c -- we're finished
c
      call retmark
c
      return
      end
c=======================================================================
c
      Subroutine Slave_Collect_Multipoles(
     $               ncs,     inx,     basdat,  ncspl,   istart,
     $               iend,   cssharps, natom,   ncf,     DENS,
     $             sharpovs,  maxovs, nsharpovs, atomicpoles,
     $             sharpness, ishellinfo)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer cssharps(ncspl)
      integer sharpness(ncs)
      real*8 DENS(ncf,ncf)
      real*8 atomicpoles(35,natom)   !atomic poles times density mx
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
      integer ishellinfo(2,ncs)
c     common /big/bl(300)
c
      do ish=istart,iend
      ics=cssharps(ish)
      iatom=inx(2,ics)
      idim=inx(3,ics)
      itype=inx(12,ics)
      if(itype .eq. 4) then !d->d6
      idim=6
      else if(itype .eq. 6) then !f->f10
      idim=10
      end if
      ifirstf=ishellinfo(1,ics)
      ilastf=ishellinfo(2,ics)
        do jsh=1,nsharpovs(ish)
        jcs=sharpovs(jsh,ish)
        if(sharpness(jcs) .eq. 3) goto 150 !the function is sharp
        jdim=inx(3,jcs)
        jtype=inx(12,jcs)
        if(jtype .eq. 4) then !d->d6
        jdim=6
        else if(jtype .eq. 6) then !f->f10
        jdim=10
        end if
        jfirstf=ishellinfo(1,jcs)
        jlastf=ishellinfo(2,jcs)
c
c      allocate arrays poles(jdim,idim,35),poleshelp(jdim,idim,35)
c
        call mmark
        call getmem(35*idim*jdim,ipoles)
        call zeroit(bl(ipoles),35*idim*jdim)
        call getmem(35*idim*jdim,ipoleshelp)
c
        call mutipoles_jcs_ics(jcs,ics,jdim,idim,jtype,
     &                         itype,ncs,inx,basdat,
     &                         bl(ipoleshelp),bl(ipoles))
c
          do ipole=1,35
            i3=idim*jdim*(ipole-1)
            icomp=0
            do ifunc=ifirstf,ilastf
              icomp=icomp+1
              i2=jdim*(icomp-1)
              jcomp=0
              do jfunc=jfirstf,jlastf
              jcomp=jcomp+1
              atomicpoles(ipole,iatom)=atomicpoles(ipole,iatom)+
     &              DENS(jfunc,ifunc)*bl(ipoles+jcomp-1+i2+i3)
              end do
            end do
          end do
c
        call retmark       !  deallocate poles,poleshelp
 150    continue
        end do
      end do
c
      return
      end
c=======================================================================
c
      Subroutine Slave_Multipoles_Fock(
     $               ncs,     inx,     basdat,  ncspl,   istart,
     $               iend,   cssharps, natom,   ncf,     XC,
     $             sharpovs,  maxovs, nsharpovs, atomicpoles,
     $             sharpness, V,   multipolecut, fskal0,  Fockmx)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      dimension XC(3,natom),atomicpoles(35,natom),V(8,natom)
      integer cssharps(ncspl)
      integer sharpness(ncs)
      real*8 Fockmx(ncf*(ncf+1)/2)
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
      integer*1 multipolecut(natom,natom)
c     common /big/bl(300)
c
c
cccccccccc
cc      write(6,*) ' In <Slave_Multipoles_Fock>  Data on input is:'
cc      write(6,*) ' ncs:',ncs,' INX array:'
cc      do i=1,ncs
cc      write(6,*) (inx(j,i),j=1,12)
cc      enddo
cc      call f_lush(6)
cc      write(6,*) ' BASDAT array:'
cc      call prntmat(246,13,246,basdat)
cc      call f_lush(6)
cc      write(6,*) ' natom:',natom,' Atomic coordinates are:'
cc      do i=1,natom
cc      write(6,*) xc(1,i),xc(2,i),xc(3,i)
cc      enddo
cc      call f_lush(6)
cc      write(6,*) ' ncspl:',ncspl,' cssharps array:'
cc      write(6,*) (cssharps(i),i=1,ncspl)
cc      call f_lush(6)
cc      write(6,*) ' maxovs:',maxovs,' sharpovs array:'
cc      do i=1,ncspl
cc      write(6,*) (sharpovs(j,i),j=1,maxovs)
cc      enddo
cc      write(6,*) ' nsharpovs array:'
cc      write(6,*) (nsharpovs(i),i=1,ncspl)
cc      write(6,*) ' sharpness array:'
cc      write(6,*) (sharpness(i),i=1,ncs)
cc      write(6,*) ' multipolecut array:'
cc      do i=1,natom
cc      write(6,*) (multipolecut(i,j),j=1,natom)
cc      enddo
cc      write(6,*) ' fskal0:',fskal0,' Fock Matrix:'
cc        do kk=1,ncf
cc        kkk = kk*(kk-1)/2
cc        write(6,*) (Fockmx(ll),ll=kkk+1,kkk+kk)
cc        enddo
cc      call f_lush(6)
ccccccccccccccccccccc
      scalfock=4.0d0/fskal0
c
      DO ish=istart,iend
        ics=cssharps(ish)
        jiatom=inx(2,ics)
        idim=inx(3,ics)
        if(idim .eq. 5) then
          idimp=6
        else if(idim .eq. 7) then
          idimp=10
        else
          idimp=idim
        end if
        ilastf=inx(10,ics)
        ifirstf=inx(11,ics)+1
        itype=inx(12,ics)
        Rxj = XC(1,jiatom)
        Ryj = XC(2,jiatom)
        Rzj = XC(3,jiatom)
        do iatom=1,natom
          if(multipolecut(iatom,jiatom) .eq. 0) then  !they are far enough
            Rxi = XC(1,iatom)
            Ryi = XC(2,iatom)
            Rzi = XC(3,iatom)
            Rx=Rxj-Rxi
            Ry=Ryj-Ryi
            Rz=Rzj-Rzi
            call multipole_intermediates(atomicpoles(1,iatom),
     &                                     Rx,Ry,Rz,V(1,iatom))
          end if
        end do
        do jsh=1,nsharpovs(ish)
          jcs=sharpovs(jsh,ish)
          if(sharpness(jcs) .eq. 3) goto 150 !the function is sharp
          jdim=inx(3,jcs)
          if(jdim .eq. 5) then
            jdimp=6
          else if(jdim .eq. 7) then
            jdimp=10
          else
            jdimp=jdim
          end if
          jlastf=inx(10,jcs)
          jfirstf=inx(11,jcs)+1
          jtype=inx(12,jcs)
c
c    allocate arrays polesi(jdimp,idimp,35),polesit(35,jdimp,idimp)
c                    fcollect(jdimp,idimp),fcollect2(jdim,idim),
c                    poleshelp(jdimp,idimp,35)
c
          lenmp=idimp*jdimp
          lenmp35=35*lenmp
          call mmark
          call getmem(lenmp35,ipolesi)
          call zeroit(bl(ipolesi),lenmp35)
          call getmem(lenmp35,ipolesit)
          call getmem(lenmp,ifcollect)
          call zeroit(bl(ifcollect),lenmp)
          call getmem(jdim*idim,ifcollect2)
          call getmem(lenmp35,ipoleshelp)
c
          call mutipoles_jcs_ics(
     &    jcs,           ics,           jdimp,         idimp,
     &    jtype,         itype,         ncs,           inx,
     &    basdat,        bl(ipoleshelp),bl(ipolesi))
c
c   now transpose the array polesi. Do not ask me why
c
          call transpole(idimp,jdimp,bl(ipolesi),bl(ipolesit))
c
          do jatom=1,natom
            if(multipolecut(jatom,jiatom) .eq. 0) then  !they are far enough
              Rxi = XC(1,jatom)
              Ryi = XC(2,jatom)
              Rzi = XC(3,jatom)
              Rx=Rxj-Rxi
              Ry=Ryj-Ryi
              Rz=Rzj-Rzi
              R=dsqrt(Rx**2+Ry**2+Rz**2)
c
c Now the main things are coming
c
              do icomp=1,idimp
                ifcoll=jdimp*(icomp-1)
                i3=35*ifcoll
                do jcomp=1,jdimp
                  i2=35*(jcomp-1)
                  call multipoles_interaction_R5(
     &        bl(ipolesit+i2+i3),atomicpoles(1,jatom),V(1,jatom),
     &                                 Rx,Ry,Rz,R,final)
c
                  bl(ifcollect+ifcoll+jcomp-1)=
     &               bl(ifcollect+ifcoll+jcomp-1)+final
                end do
              end do
c
            end if
          end do !for jatom
c
c Now transform it back if necessary
c
          if((idim .ne. idimp) .or. (jdim .ne. jdimp)) then
c
c  allocate array work(jdim,idimp)
c
            call getmem(jdim*idimp,iwork)
c
            call shellpair_trafo_from_descartes(idim,jdim,idimp,
     &           jdimp,bl(ifcollect2),bl(ifcollect),bl(iwork))
            call retmem(1) ! deallocate work
c
            if(ics .gt. jcs)then
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc*(ifunc-1)/2
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect2+jdim*(icomp-1)+jcomp-1)
                end do
              end do
            else
              icomp=0
              do ifunc=ifirstf,ilastf
                icomp=icomp+1
                jcomp=0
                indxi=ifunc
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc*(jfunc-1)/2
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect2+jdim*(icomp-1)+jcomp-1)
                end do
              end do
            end if
c
          else
c
            if(ics .gt. jcs)then
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc*(ifunc-1)/2
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect+jdimp*(icomp-1)+jcomp-1)
                end do
              end do
            else
              icomp=0
              do ifunc=ifirstf,ilastf
                indxi=ifunc
                icomp=icomp+1
                jcomp=0
                do jfunc=jfirstf,jlastf
                  jcomp=jcomp+1
                  indxij=indxi+jfunc*(jfunc-1)/2
                  Fockmx(indxij)=Fockmx(indxij)+
     &                  scalfock*bl(ifcollect+jdimp*(icomp-1)+jcomp-1)
                end do
              end do
            end if
c
          end if
c
          call retmark  ! deallocate polesi,polesit,fcollect,fcollect2,
                        ! poleshelp
 150      continue
        end do
      END DO
C
      return
      end
c=======================================================================
c
      SUBROUTINE slave_sharp1_grad(
     $       ncs,           ncspl,         icssharps, igridrange2,
     $       ICSPLSIZE,     iicspltype,    Lxmin,     Lymin,
     $       Lzmin,         griddens,      iplbasdat, iisharpgridrange2,
     $       iro4,          isharpness,    idens,     iisharpgridrange1,
     $       maxovs,        insharpovs,    isharpovs, npwx,
     $       npwy,          npwz,          icorecut,  INX,
     $       iystoredim,    iyexpstore,  iy2expstore, iy3expstore,
     $       izstoredim,    izexpstore,  iz2expstore, iz3expstore,
     $       ipexpdim,    iipyexpstore, iipzexpstore, ibas,
     $       ncf,           ncfpl,         natom,     NTot,
     $       NSLAVE,        GRADX,         GRADY,     GRADZ)
c
c  slave code for calculating the first part of the gradient
c  using the mixed density
c

      use memory

      implicit real*8 (a-h,o-z)
      real*8 Lxmin,Lymin,Lzmin
      real*8 GRADX(natom),GRADY(natom),GRADZ(natom)
      integer NSLAVE(NTot),INX(12,*),ICSPLSIZE(ncspl)
c .................................................
c -- dynamically allocated array
      integer ICFPLT(ncspl)
c .................................................
c     common /big/ bl(300)
c
cc      write(6,*) ' SLAVE:  In <slave_sharp1_grad>  NTot:',ntot
cc      write(6,*) ' NSLAVE: ',(nslave(i),i=1,ntot)
cc      write(6,*) ' ICSPLSIZE:',(icsplsize(i),i=1,ncspl)
cc      call f_lush(6)
      npwyz = npwy*npwz
      call getmem(npwyz,iro3)
      call ZeroIT(bl(iro3),npwyz)
      ICFPLT(1) = 0
      do i=1,ncspl-1
      ICFPLT(i+1) = ICFPLT(i) + ICSPLSIZE(i)
      enddo
cc      write(6,*) ' ICFPLT:',(icfplt(i),i=1,ncspl)
cc      call f_lush(6)
c
      dx = 1.0d0/griddens
      npwregx=1
      xmin0=Lxmin-dx
c
      DO ish=1,NTot
      call mmark
      icspl = NSLAVE(ish)
      icfpl = ICFPLT(icspl)
cc      write(6,*) ' SLAVE:  ish:',ish,' icspl:',icspl,' icfpl:',icfpl
cc      call f_lush(6)
c
      call AFT_sharpgrid_derivatives_arrange(
     &    ncs,           ncspl,         bl(icssharps), bl(igridrange2),
     &    ICSPLSIZE,   bl(iicspltype),  Lxmin,         Lymin,
     &    Lzmin,         griddens,      icspl,         ncontr,
     &    rLx,           rLy,           rLz,           xp,
     &    yp,            zp,            npwxs,         npwys,
     &    npwzs,         ics,           itype,         nzeffdim,
     &    iexpfuncx,     iexpfuncy,     iexpfuncz,     icftwork,
     &    icftwork2,     iftwork,       isxmin,        isxmax,
     &    isymin,        isymax,        iszmin,        iszmax,
     &    INX,           ifirstf,       ilastf,  bl(iisharpgridrange2),
     &    icorecut)
c
      ymin=Lymin+float(isymin-1)*dx
      npwregy=npwys
      ishift1help=3*npwys*npwzs
      ishift2help=npwyz
      valgradx=0.0d0
      valgrady=0.0d0
      valgradz=0.0d0
      ncomp=ilastf-ifirstf+1
cc      write(6,*) ' SLAVE: ifirstf' ,ifirstf,' ilastf',ilastf
cc      call f_lush(6)
c
      do ifunc=ifirstf,ilastf
      icomp=ifunc-ifirstf+1
c
      call make_sharpgrid_derivatives_work(
     &    itype,         icomp,         npwxs,         npwys,
     &    npwzs,         rLx,           rLy,           rLz,
     &    bl(icftwork),  ncontr,        icfpl,         bl(iplbasdat),
     &    ncfpl,         xp,            yp,            zp,
     &    bl(iftwork),   bl(iexpfuncx), bl(iexpfuncy), bl(iexpfuncz),
     &    bl(icftwork2), nzeffdim,      bl(icftwork))
c
c -- now the derivatives (3D) for the given ifunc are in the ftwork array
c
      if(icomp .ne. ncomp) icfpl=icfpl-ncontr
c
      ishift1=0
      ishift2=(isxmin-1)*npwyz+(isymin-1)*npwz
      xmin=xmin0+(isxmin-1)*dx
c
      do ix=isxmin,isxmax
      xmin = xmin+dx
c
      call calc_mix1_gradients(
     &    ncf,      bl(iftwork+ishift1),npwregx,       npwregy,
     &    npwz,        bl(igridrange2), ncs,         bl(isharpness),
     &    ncspl,         bl(iro3),      bl(idens),     INX,
     &    bl(ibas),      xmin,          ymin,          Lzmin,
     &    ix,            isymin,        isymax,        icspl,
     &    ifunc,         bl(isharpovs), maxovs,      bl(insharpovs),
     &    npwzs,      bl(iro4+ishift2), griddens,      iystoredim,
     &  bl(iyexpstore),bl(iy2expstore),bl(iy3expstore),izstoredim,
     &  bl(izexpstore),bl(iz2expstore),bl(iz3expstore),ipexpdim,
     &  bl(iipyexpstore),bl(iipzexpstore),iszmin,      iszmax,
     &    valgradx,     valgrady,        valgradz,     icorecut,
     &  bl(iisharpgridrange1), bl(iisharpgridrange2))
c
      ishift1=ishift1+ishift1help
      ishift2=ishift2+ishift2help
      enddo
      enddo
c
      call add_to_gradients(ncs,INX,ics,natom,valgradx,
     &           valgrady,valgradz,GRADX,GRADY,GRADZ)
c
      call retmark
c
      EndDO
c
      call retmem(1)      ! deallocate ro3
C
      RETURN
      END
c=======================================================================
c
      SUBROUTINE slave_sharp2_grad(
     $       ncs,           ncspl,         icssharps, igridrange2,
     $       ii4sharps,     iicspltype,    Lxmin,     Lymin,
     $       Lzmin,         griddens,      isharpgrd, iisharpgridrange2,
     $       iro4,          isharpness,    idens,     iisharpgridrange1,
     $       maxovs,        insharpovs,    isharpovs, npwx,
     $       npwy,          npwz,          icorecut,  INX,
     $       iystoredim,    iyexpstore,  iy2expstore, iy3expstore,
     $       izstoredim,    izexpstore,  iz2expstore, iz3expstore,
     $       ipexpdim,    iipyexpstore, iipzexpstore, ibas,
     $       ncf,           natom,         NTot,      NSLAVE,
     $       GRADX,         GRADY,         GRADZ)
c
c  slave code for calculating the second part of the gradient
c  using the mixed density
c

      use memory

      implicit real*8 (a-h,o-z)
      real*8 Lxmin,Lymin,Lzmin
      real*8 GRADX(natom),GRADY(natom),GRADZ(natom)
      integer NSLAVE(NTot),INX(12,*)
c     common /big/ bl(300)
c
      npwyz = npwy*npwz
      call getmem(3*npwyz,iro5v)
      call ZeroIT(bl(iro5v),3*npwyz)
c
      dx = 1.0d0/griddens
      npwregx=1
      xmin0=Lxmin-dx
      ifunc=0
      indxcomprcoef=0
c
      DO ish=1,NTot
      call mmark
      icspl = NSLAVE(ish)
c
      call mixed_density_help1(
     &    icspl,         bl(iicspltype),ncspl,         bl(icssharps),
     &    ncomp,       bl(igridrange2), ncs,           isxmin,
     &    isxmax,        isymin,        isymax,        iszmin,
     &    iszmax,        INX,           ifirstf,       ilastf,
     &    icorecut,bl(iisharpgridrange2))
c
      npwsx=isxmax-isxmin+1
      npwsy=isymax-isymin+1
      npwsz=iszmax-iszmin+1
      ymin=Lymin+float(isymin-1)*dx
      npwregy=isymax-isymin+1
      ishifthelp1=(isymin-1)*npwz
      ishifthelp2=npwregy*npwz+(npwy-isymax)*npwz
      ishift0=(isxmin-1)*npwy*npwz
      ishift2help1=npwregy*npwsz
c
      call getint(npwsy*npwsz,iro62d)
c
      do ifunc=ifirstf,ilastf
      ishift2=0
      ishift=ishift0
      xmin=xmin0+(isxmin-1)*dx
c
      do ix=isxmin,isxmax ! loop over the x region of the given sharp bf
      indxcomprcoef=indxcomprcoef+1
      comprcoef=bl(ii4sharps+indxcomprcoef-1)
c
      call readin_sharps2d(isharpgrd,npwsy,npwsz,bl(iro62d))
c
      xmin=xmin+dx
      if(isymin .gt. 1) ishift=ishift+ishifthelp1
c
      call calc_mix2_gradients(
     &    ncf,           bl(iro62d),    npwregx,       npwregy,
     &    npwz,        bl(igridrange2), ncs,           bl(isharpness),
     &    ncspl,         bl(idens),     INX,           bl(ibas),
     &    xmin,          ymin,          Lzmin,         bl(iro5v),
     &    ix,            isymin,        isymax,        icspl,
     &    ifunc,         bl(isharpovs), maxovs,        bl(insharpovs),
     &    npwsz,         iszmin,      bl(iro4+ishift), griddens,
     &    iystoredim,  bl(iyexpstore),bl(iy2expstore),bl(iy3expstore),
     &    izstoredim,  bl(izexpstore),bl(iz2expstore),bl(iz3expstore),
     &    ipexpdim,  bl(iipyexpstore),bl(iipzexpstore),comprcoef,
     &    natom,         GRADX,         GRADY,         GRADZ,
     &    icorecut, bl(iisharpgridrange1),bl(iisharpgridrange2))
c
      ishift=ishift+ishifthelp2
      ishift2=ishift2+ishift2help1
      enddo
      enddo
c
      call retmark        ! deallocate ro62d
c
      EndDO
c
      call retmem(1)      ! deallocate ro5v
C
      RETURN
      END

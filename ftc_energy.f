      Subroutine FTC_ENERGY(
     &     griddens,   Lxmin,      Lymin,      Lzmin,      Lxo,
     &     Lyo,        Lzo,        Lxe,        Lye,        Lze,
     &     PLDmax,     npwx,       npwy,       npwz,       npwxe,
     &     npwye,      npwze,      isharpness, iinit,      ncspl,
     &     igridranges,igridranges2,icssharps, iicsdiff,   iicspltype,
     &     maxovs,     isharpovs,  insharpovs, ii4sharps,  icorecut,
     &iisharpgridrange1,iisharpgridrange2,denA,denB,       fockA,
     &     fockB,      thint,      dist2mp,  imultipolecut,nomultipole,
     &     conver)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      real*8 denA(*),denB(*),fockA(*),fockB(*)
      character*256 scrfile,filename
      logical conver
      parameter (isharpgrd=51)            ! unit number for FTC I/O
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(5000000)
cccccc
      common/ftctime/tinit,einit,tsmooth,esmooth,tftc1,eftc1,
     $               tmix,emix,tftc2,eftc2,tfock1,efock1,tmpole,empole
cccccc
c
c
      call elapsec(ee1)
      call secund(tt1)
c
      call getival('ncs ',ncs)
      call getival('ncf ',ncf)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('inuc',inuc)
c
      ntri = (ncf*(ncf+1))/2
      ncsdiff=ncs-ncspl
      dx=1.0d0/griddens
c
c -- determine value of Fock matrix scaling factor
      dv=(Lxo/dfloat(npwx)) * (Lyo/dfloat(npwy)) * (Lzo/dfloat(npwz))
      constfftw=1.0d0/npwxe/npwye/npwze
      call setrval('cfftw',constfftw)
      fskal0 = dv*constfftw
c
c -- scale the input Fock matrix
c    (contains Coulomb contribution from classical integrals)
      fskal = thint/fskal0
      call vscal(ntri,fskal,fockA)
c
      call elapsec(ee2)
      call secund(tt2)
      tinit = tinit + tt2-tt1
      einit = einit + ee2-ee1
c
c=======================================================================
c Now build up the smooth density
c
      call elapsec(ee1)
      call secund(tt1)
c
      call mmark
c
c    allocate  the array ro4(npwz,npwy,npwx)
c
      npwxyz=npwz*npwy*npwx
      call getmem(npwxyz,iro4)
      call zeroit(bl(iro4),npwxyz)
c
      npwyz=npwz*npwy
c
      call mmark
c
c    allocate arrays ro1(npwz,npwy) and ro2(npwz,npwy)
c
      call getmem(npwyz,iro1)
      call getmem(npwyz,iro2)
      call zeroit(bl(iro1),npwyz)
      call zeroit(bl(iro2),npwyz)
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges),iystoredim,izstoredim,ipexpdim)
c
c    allocate arrays yexpstore(iystoredim) and zexpstore(izstoredim)
c
      call getmem(iystoredim,iyexpstore)
      call getmem(izstoredim,izexpstore)
      call zeroit(bl(iyexpstore),iystoredim)
      call zeroit(bl(izexpstore),izstoredim)
c
c     allocate integer arrays ipyexpstore(ipexpdim) and
c                             ipzexpstore(ipexpdim)
c
      call getint(ipexpdim,iipyexpstore)
      call getint(ipexpdim,iipzexpstore)
      call izeroit(bl(iipyexpstore),ipexpdim)
      call izeroit(bl(iipzexpstore),ipexpdim)
c
      call calc_exp_store_yz_2vec(
     &    ncsdiff,         bl(iicsdiff),   ncs,       bl(ictr),
     &    bl(ibas),        bl(igridranges),iystoredim,bl(iyexpstore),
     &    izstoredim,      bl(izexpstore), ipexpdim,  bl(iipyexpstore),
     &    bl(iipzexpstore),Lymin,          Lzmin,     griddens)
c
      do ix=1,npwx
c
        x=Lxmin+float(ix-1)*dx
        iro4p=npwyz*(ix-1)
c
        call build_smooth_density(
     &    ncf,            bl(iro1),       npwy,          npwz,
     &    bl(igridranges),bl(iro2),       ncs,           bl(isharpness),
     &    ncspl,          denA,           bl(ictr),     bl(ibas),
     &    x,              Lymin,          Lzmin,         bl(iro4+iro4p),
     &    ix,             1,              npwy,          griddens,
     &    iystoredim,     bl(iyexpstore), izstoredim,    bl(izexpstore),
     &    ipexpdim,       bl(iipyexpstore),bl(iipzexpstore))
c
      end do
c
      call retmark    ! deallocate ro1,ro2,yexpstore,zexpstore,
                      ! ipyexpstore,ipzexpstore
c
      call elapsec(ee2)
      call secund(tt2)
      tsmooth = tsmooth + tt2-tt1
      esmooth = esmooth + ee2-ee1
c
c Now the smooth density is in ro4
c Now make the Coulomb potential from the smooth density !
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE SMOOTH DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      call elapsec(ee1)
      call secund(tt1)
c
      call mmark
c
c     allocate the real arrays roe(npwz,npwy,npwxe),w1(npwxe+15)
c       w3(2*npwye+15),w4(2*npwze+15),vky(npwye),vkz(npwze)
c     and the complex array vyze(npwze,npwye)
c
      call getmem(npwyz*npwxe,iroe)
      call zeroit(bl(iroe),npwyz*npwxe)
      call getmem(npwxe+15,iw1)
      call getmem(2*npwye+15,iw3)
      call getmem(2*npwze+15,iw4)
      call getmem(npwye,ivky)
      call getmem(npwze,ivkz)
      call getmem(2*npwze*npwye,ivyze)
c
c Now create the Coulomb potential from the smooth density in ro4
c
      call make_COULOMB_on_FTGRID(
     &    npwx,       npwy,       npwz,      npwxe,        npwye,
     &    npwze,      bl(iro4),   PLDmax,    Lxe,          Lye,
     &    Lze,        bl(iroe),   bl(iw1),   bl(ivyze),    bl(iw3),
     &    bl(iw4),    bl(ivky),   bl(ivkz))
c
      call retmark     ! deallocate vxyze,vkz,vky,w4,w3,w1,roe
c
c Now the Coulomb potential of the smooth density is in ro4
c Save the smooth potential if this is the last cycle
      If(conver) CALL SAVE_POT(npwxyz,'potS',bl(iro4))
c
      call secund(tt2)
      call elapsec(ee2)
      tftc1 = tftc1 + tt2-tt1
      eftc1 = eftc1 + ee2-ee1
c
c=======================================================================
c Now make the mixed (sharp*smooth) density (it will go to ro5)
c and the corresponding Fock matrix elements
c
      call elapsec(ee1)
      call secund(tt1)
c
c    allocate  the array ro5(npwz,npwy,npwx)
c
      call getmem(npwxyz,iro5)
      call zeroit(bl(iro5),npwxyz)
c
      call mmark
c
c    allocate arrays ro1(npwz,npwy) and ro3(npwz,npwy)
c
      call getmem(npwyz,iro1)
      call getmem(npwyz,iro3)
      call zeroit(bl(iro1),npwyz)
      call zeroit(bl(iro3),npwyz)
c
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges2),iystoredim,izstoredim,ipexpdim)
c
c    allocate arrays yexpstore(iystoredim) and zexpstore(izstoredim)
c
      call getmem(iystoredim,iyexpstore)
      call getmem(izstoredim,izexpstore)
      call zeroit(bl(iyexpstore),iystoredim)
      call zeroit(bl(izexpstore),izstoredim)
c
c     allocate integer arrays ipyexpstore(ipexpdim) and
c                             ipzexpstore(ipexpdim)
c
      call getint(ipexpdim,iipyexpstore)
      call getint(ipexpdim,iipzexpstore)
      call izeroit(bl(iipyexpstore),ipexpdim)
      call izeroit(bl(iipzexpstore),ipexpdim)
c
      call calc_exp_store_yz_2vec(
     &    ncsdiff,        bl(iicsdiff),    ncs,        bl(ictr),
     &    bl(ibas),       bl(igridranges2),iystoredim, bl(iyexpstore),
     &    izstoredim,     bl(izexpstore),  ipexpdim,   bl(iipyexpstore),
     &    bl(iipzexpstore),Lymin,          Lzmin,      griddens)
c
c  open file
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len1)
      filename = scrfile(1:len1)//'.sharpgrd'
      OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
c
      xmin0=Lxmin-dx
      indxcomprcoef=0
c
      do icspl=1,ncspl ! Loop over the sharp shells
c
        call mixed_density_help1(
     &    icspl,         bl(iicspltype),  ncspl,       bl(icssharps),
     &    ncomp,         bl(igridranges2),ncs,         isxmin,
     &    isxmax,        isymin,          isymax,      iszmin,
     &    iszmax,        bl(ictr),       ifirstf,     ilastf,
     &    icorecut,      bl(iisharpgridrange2))
c
        npwsx=isxmax-isxmin+1
        npwsy=isymax-isymin+1
        npwsz=iszmax-iszmin+1
        ymin=Lymin+float(isymin-1)*dx
        npwregy=isymax -isymin+1
c
c    allocate  the integer*4 array ro62d(npwsz,npwsy)
c
        call getint_4(npwsz*npwsy,iro62d)
c
        do ifunc=ifirstf,ilastf
          xmin=xmin0+(isxmin-1)*dx
          do ix=isxmin,isxmax ! loop over x region of the given sharp bf
            indxcomprcoef=indxcomprcoef+1
            comprcoef=bl(ii4sharps+indxcomprcoef-1)
c
            call readin_sharps2d(isharpgrd,npwsy,npwsz,bl(iro62d))
c
            xmin=xmin+dx
            iro4p=npwz*(isymin-1)+npwyz*(ix-1)
c
c           write(*,*)'ifunc,ix',ifunc,ix
            call make_mixed_dpf(
     & ncf,          bl(iro62d),        npwregy,         npwz,
     & bl(igridranges2),bl(iro1),       ncs,             bl(isharpness),
     & ncspl,        nfsh,              igrdsh,          bl(iro3),
     & denA,         bl(ictr),         bl(ibas),        xmin,
     & ymin,         Lzmin,             bl(iro5+iro4p),  ix,
     & isymin,       isymax,            icspl,           ifunc,
     & bl(isharpovs),maxovs,            bl(insharpovs),  npwsz,
     & iszmin,       bl(iro4+iro4p),    fockA,           griddens,
     & iystoredim,   bl(iyexpstore),    izstoredim,      bl(izexpstore),
     & ipexpdim,     bl(iipyexpstore),  bl(iipzexpstore),comprcoef,
     & icorecut,     bl(iisharpgridrange1),bl(iisharpgridrange2))
c
          end do !end of the loop over ix
c
        end do ! end of the loop over the component
c
        call retmem(1)          ! deallocate ro62d
c
      end do ! end of the loop over the sharp shells
c
      CLOSE(UNIT=isharpgrd)
c     stop
c
      call retmark    ! deallocate ro1,ro3,yexpstore,zexpstore,
                      ! ipyexpstore,ipzexpstore
c
      call elapsec(ee2)
      call secund(tt2)
      tmix = tmix + tt2-tt1
      emix = emix + ee2-ee1
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE MIXED DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      call elapsec(ee1)
      call secund(tt1)
c
      call mmark
c
c     allocate the real arrays roe(npwz,npwy,npwxe),w1(npwxe+15)
c       w3(2*npwye+15),w4(2*npwze+15),vky(npwye),vkz(npwze)
c     and the complex array vyze(npwze,npwye)
c
      call getmem(npwyz*npwxe,iroe)
      call zeroit(bl(iroe),npwyz*npwxe)
      call getmem(npwxe+15,iw1)
      call getmem(2*npwye+15,iw3)
      call getmem(2*npwze+15,iw4)
      call getmem(npwye,ivky)
      call getmem(npwze,ivkz)
      call getmem(2*npwze*npwye,ivyze)
c
c Now create the Coulomb potential from the mixed density in ro5
c
      call make_COULOMB_on_FTGRID(
     &    npwx,       npwy,       npwz,      npwxe,        npwye,
     &    npwze,      bl(iro5),   PLDmax,    Lxe,          Lye,
     &    Lze,        bl(iroe),   bl(iw1),   bl(ivyze),    bl(iw3),
     &    bl(iw4),    bl(ivky),   bl(ivkz))
c
      call retmark     ! deallocate vxyze,vkz,vky,w4,w3,w1,roe
c
c Now the Coulomb potential of the mixed density is in ro5
c Save the mixed potential if this is the last cycle
cc      If(conver) CALL SAVE_POT(npwxyz,'potM',bl(iro5))
c
c    sum smooth and mixed potentials. The result will be in ro5
c
      call daxpy(npwxyz,1.0d0,bl(iro4),1,bl(iro5),1)
      If(conver) CALL SAVE_POT(npwxyz,'potM',bl(iro5))
c
c
      call elapsec(ee2)
      call secund(tt2)
      tftc2 = tftc2 + tt2-tt1
      eftc2 = eftc2 + ee2-ee1
c
c
c At this point the mixed+smooth potential is in the ro5 array
c Now calculate the Fock matrix elements from the smooth+mixed density
c
      call elapsec(ee1)
      call secund(tt1)
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges),iystoredim,izstoredim,ipexpdim)
c
      call mmark
c
c    allocate arrays ro1(npwz,npwy), ro2(npwz,npwy) and ro3(npwz,npwy)
c
      call getmem(npwyz,iro1)
      call getmem(npwyz,iro2)
      call getmem(npwyz,iro3)
      call zeroit(bl(iro1),npwyz)
      call zeroit(bl(iro2),npwyz)
      call zeroit(bl(iro3),npwyz)
c
c    allocate arrays yexpstore(iystoredim) and zexpstore(izstoredim)
c
      call getmem(iystoredim,iyexpstore)
      call getmem(izstoredim,izexpstore)
      call zeroit(bl(iyexpstore),iystoredim)
      call zeroit(bl(izexpstore),izstoredim)
c
c     allocate integer arrays ipyexpstore(ipexpdim) and
c                             ipzexpstore(ipexpdim)
c
      call getint(ipexpdim,iipyexpstore)
      call getint(ipexpdim,iipzexpstore)
      call izeroit(bl(iipyexpstore),ipexpdim)
      call izeroit(bl(iipzexpstore),ipexpdim)
c
      call calc_exp_store_yz_2vec(
     &    ncsdiff,       bl(iicsdiff),   ncs,         bl(ictr),
     &    bl(ibas),      bl(igridranges),iystoredim, bl(iyexpstore),
     &    izstoredim,    bl(izexpstore), ipexpdim,    bl(iipyexpstore),
     &  bl(iipzexpstore),Lymin,          Lzmin,       griddens)
cc
      do ix=1,npwx
c
        x=Lxmin+float(ix-1)*dx
        iro5p=npwyz*(ix-1)
c
        call calc_FTC_fock_smooth(
     &    ncf,           bl(iro1),      1,             npwy,
     &    npwz,          bl(igridranges),bl(iro2),     ncs,
     &    bl(isharpness),ncspl,         igrdsh,
     &    bl(iro3),      bl(ictr),      bl(ibas),       x,
     &    Lymin,         Lzmin,         bl(iro5+iro5p),ix,
     &    1,             npwy,          fockA,         griddens,
     &    iystoredim,    bl(iyexpstore),izstoredim,    bl(izexpstore),
     &    ipexpdim,      bl(iipyexpstore),bl(iipzexpstore))
c
      end do
c
      call retmark    ! deallocate ro1,ro2,ro3,yexpstore,zexpstore,
                      ! ipyexpstore,ipzexpstore
c
      call retmark    ! deallocate ro4,ro5
c
      call elapsec(ee2)
      call secund(tt2)
      tfock1 = tfock1 + tt2-tt1
      efock1 = efock1 + ee2-ee1
c
cc        write(6,*) ' At start of multipoles section'
cc        write(6,*) ' Fock Matrix is:'
cc        do kk=1,ncf
cc        kkk = kk*(kk-1)/2
cc        write(6,*) (fockA(ll),ll=kkk+1,kkk+kk)
cc        enddo
      if(nomultipole .eq. 0) then
ct      call getival('iroute',iroute)
ct        if(iroute .ne. 1) then
ct        write(6,*)"FTC MULTIPOLE ERROR !"
ct        write(6,*)"THE MULTIPOLE CODE WORKS ONLY WITH ROUTE=1 !"
ct        write(6,*)" PUT THIS KEYWORD INTO THE inte LINE"
ct        write(6,*)"THE PROGRAM WILL STOP NOW !"
ct        stop
ct        end if
c
c Now the multipole stuff is coming
        call getival('na  ',na)
        call getival('ndum',ndum)
        natom = na-ndum
c
c     allocate array denstmp(ncf,ncf)
c
        call getmem(ncf*ncf,idenstmp)
c
        call fillin_dens(bl(idenstmp),denA,ncf)
c
        call elapsec(ee1)
        call secund(tt1)
c
        call multipoles_main(
     &    ncs,           bl(ictr),     bl(ibas),      ncspl,
     &    bl(icssharps), natom,         ncf,           bl(idenstmp),
     &    bl(isharpovs), maxovs,        bl(insharpovs),dist2mp,
     &    bl(isharpness),bl(inuc),   bl(imultipolecut),fskal0,
     &    fockA)
c
        call retmem(1)     !     deallocate denstmp
c
        call elapsec(ee2)
        call secund(tt2)
c
        tmpole = tmpole + tt2-tt1
        empole = empole + ee2-ee1
c
      else
        tmpole = 0.0d0
        empole = 0.0d0
      end if
cc        write(6,*) ' At end of multipoles section'
cc        write(6,*) ' Fock Matrix is:'
cc        do kk=1,ncf
cc        kkk = kk*(kk-1)/2
cc        write(6,*) (fockA(ll),ll=kkk+1,kkk+kk)
cc        enddo
c
c Multipole stuff ends
c
c -- scale final Fock matrix
      call Make_FockL(ncf,fskal0,fockA)
C
      return
      end
c***********************************************************************
c
      Subroutine calc_exp_dims(ncsdiff,icsdiff,ncs,inx,basdat,
     &           gridranges,iystoredim,izstoredim,ipexpdim)
      IMPLICIT REAL*8(A-H,O-Z)
c
c This small subroutine calculates the grid dimensions along y and z
c for pre calculating the exp functions for smooth basis functions.
c
c Input:
c
c ncsdiff               nr. of smooth contracted shell
c gridranges            coordinate space grid ranges for the
c                       given contr. shell
c inx,basdat            as always
c icsdiff               smooth contracted shell indexes
c
c Output:
c
c iystoredim,izstoredim sum of the grid ranges along y and z for the all
c                               smooth prim. shells
c ipexpdim              max among the inx(5,ics) (primitive shell index)
c                       ics means a smooth contr. shell here.
c
      integer icsdiff(ncsdiff)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer gridranges(6,ncs)
c
      iystoredim=1
      izstoredim=1
      ipexpdim=0
c
      do ish=1,ncsdiff ! loop over the smooth shells
        ics=icsdiff(ish)
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iydim=iymax-iymin+1
        izmin=gridranges(5,ics)
        izmax=gridranges(6,ics)
        izdim=izmax-izmin+1
        ibeg=inx(1,ics)+1
        iend=inx(5,ics)
        do icont=ibeg,iend  ! loop over the contractions
          iystoredim=iystoredim+iydim
          izstoredim=izstoredim+izdim
        end do
        ipexpdim=max(ipexpdim,iend)
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine calc_exp_store_yz(
     &    ncsdiff,       icsdiff,       ncs,           inx,
     &    basdat,        gridranges,    iystoredim,    yexpstore,
     &    y2expstore,    y3expstore,    izstoredim,    zexpstore,
     &    z2expstore,    z3expstore,    ipexpdim,      ipyexpstore,
     &    ipzexpstore,   Lymin,         Lzmin,         griddens)
c
c This subroutine fills in the array for  the 1d exponential functions
c (along y and z) for all smooth prim. shells.
c
c Input:
c
c ncsdiff               nr. of smooth contracted shell
c gridranges            coordinate space grid ranges for the given
c                       contr. shell
c inx,basdat,ncs        as always
c icsdiff               smooth contracted shell indexes
c iystoredim,izstoredim summ of the grid ranges along y and z for
c                       the all smooth prim. shells
c ipexpdim              max among the inx(5,ics) (primitive shell index)
c                       ics-> smooth
c Lymin,Lzmin           real box starting positions along y and z
c griddens              grid density
c
c Output:
c
c yexpstore,y2expstore, exp(-eta*y**2);y*exp(-eta*y**2);
c                       y**2*exp(-eta*y**2) along
c y3expstore            the grid for all prim. shell
c zexpstore,z2expstore, exp(-eta*z**2);z*exp(-eta*z**2);
c                       z**2*exp(-eta*z**2) along
c z3expstore            the grid for all prim. shell
c ipyexpstore,ipzexpstore   pointers to the first positions of the
c                           given 1D grid function of the given
c                           prim. shell
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer icsdiff(ncsdiff)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer gridranges(6,ncs)
      real*8 yexpstore(iystoredim)
      real*8 y2expstore(iystoredim)
      real*8 y3expstore(iystoredim)
      real*8 zexpstore(izstoredim)
      real*8 z2expstore(izstoredim)
      real*8 z3expstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      real*8 Lymin,Lzmin
c
      ipy=1
      ipz=1
c
      do ish=1,ncsdiff
        ics=icsdiff(ish)
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iydim=iymax-iymin+1
        izmin=gridranges(5,ics)
        izmax=gridranges(6,ics)
        izdim=izmax-izmin+1
        ibeg=inx(1,ics)+1
        Py=basdat(12,ibeg)
        Pz=basdat(13,ibeg)
        iend=inx(5,ics)
        do icont=ibeg,iend  ! loop over the contractions
          ipyexpstore(icont)=ipy
          ipzexpstore(icont)=ipz
          eta=basdat(1,icont)
c
          call precalc_expval_yz(
     &    yexpstore(ipy),y2expstore(ipy),y3expstore(ipy),zexpstore(ipz),
     &    z2expstore(ipz),z3expstore(ipz),Py,          Pz,
     &    iymin,         iymax,         izmin,         izmax,
     &    Lymin,         Lzmin,         eta,           griddens,
     &    iydim,         izdim)
c
          ipy=ipy+iydim
          ipz=ipz+izdim
        end do
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine calc_exp_store_yz_2vec(
     &    ncsdiff,       icsdiff,       ncs,           inx,
     &    basdat,        gridranges,    iystoredim,    yexpstore,
     &    izstoredim,    zexpstore,     ipexpdim,      ipyexpstore,
     &    ipzexpstore,   Lymin,         Lzmin,         griddens)
c
c This subroutine fills in the array for the 1d exponential functions
c (along y and z) for all smooth prim. shells.
c
c Input:
c
c ncsdiff               nr. of smooth contracted shell
c gridranges            coordinate space grid ranges for the
c                       given contr. shell
c inx,basdat,ncs        as always
c icsdiff               smooth contracted shell indexes
c iystoredim,izstoredim sum of the grid ranges along y and z for
C                       the all smooth prim. shells
c ipexpdim              max among the inx(5,ics) (primitive shell index)
c                       ics-> smooth
c Lymin,Lzmin           real box starting positions along y and z
c griddens              grid density
c
c Output:
c
c yexpstore,y2expstore, exp(-eta*y**2);y*exp(-eta*y**2);
c                       y**2*exp(-eta*y**2) along
c y3expstore            the grid for all prim. shell
c zexpstore,z2expstore, exp(-eta*z**2);z*exp(-eta*z**2);
c                       z**2*exp(-eta*z**2) along
c z3expstore            the grid for all prim. shell
c ipyexpstore,ipzexpstore pointers to the first positions of the
c                         given 1D grid function of the given
c                         prim. shell
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer icsdiff(ncsdiff)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      integer gridranges(6,ncs)
      real*8 yexpstore(iystoredim)
      real*8 zexpstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      real*8 Lymin,Lzmin
c
      ipy=1
      ipz=1
c
      do ish=1,ncsdiff
        ics=icsdiff(ish)
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iydim=iymax-iymin+1
        izmin=gridranges(5,ics)
        izmax=gridranges(6,ics)
        izdim=izmax-izmin+1
        ibeg=inx(1,ics)+1
        Py=basdat(12,ibeg)
        Pz=basdat(13,ibeg)
        iend=inx(5,ics)
        do icont=ibeg,iend  ! loop over the contractions
          ipyexpstore(icont)=ipy
          ipzexpstore(icont)=ipz
          eta=basdat(1,icont)
c
          call precalc_expval_yz_2vec(
     &    yexpstore(ipy),zexpstore(ipz),Py,            Pz,
     &    iymin,         iymax,         izmin,         izmax,
     &    Lymin,         Lzmin,         eta,           griddens,
     &    iydim,         izdim)
c
          ipy=ipy+iydim
          ipz=ipz+izdim
        end do
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine precalc_expval_yz(
     &    expvaly,    yexpvaly,   yyexpvaly,  expvalz,    zexpvalz,
     &    zzexpvalz,  Py,         Pz,         iymin,      iymax,
     &    izmin,      izmax,      Lymin,      Lzmin,      eta,
     &    griddens,   npwy,       npwz)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin,Lymin,Lzmin
      real*8 expvaly(npwy)
      real*8 yexpvaly(npwy)
      real*8 yyexpvaly(npwy)
      real*8 expvalz(npwz)
      real*8 zexpvalz(npwz)
      real*8 zzexpvalz(npwz)
c
      dx=1.0d0/griddens
      yinit=Lymin - Py - dx
      zinit=Lzmin - Pz - dx
c
      index=1
      do iy=iymin,iymax
        y=yinit+float(iy)*dx
        r2=-eta*y**2
        valexp=exp(r2)
        expvaly(index)=valexp
        yexpvaly(index)=y*valexp
        yyexpvaly(index)=y**2*valexp
        index=index+1
      end do
c
      index=1
c
      do iz=izmin,izmax
        z=zinit+float(iz)*dx
        r2=-eta*z**2
        valexp=exp(r2)
        expvalz(index)=valexp
        zexpvalz(index)=z*valexp
        zzexpvalz(index)=z**2*valexp
        index=index+1
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine precalc_expval_yz_2vec(
     &    expvaly,    expvalz,    Py,         Pz,         iymin,
     &    iymax,      izmin,      izmax,      Lymin,      Lzmin,
     &    eta,        griddens,   npwy,       npwz)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin,Lymin,Lzmin
      real*8 expvaly(npwy)
      real*8 expvalz(npwz)
c
      dx=1.0d0/griddens
      yinit=Lymin - Py - dx
      zinit=Lzmin - Pz - dx
c
      index=1
      do iy=iymin,iymax
        y=yinit+float(iy)*dx
        r2=-eta*y**2
        valexp=exp(r2)
        expvaly(index)=valexp
        index=index+1
      end do
c
      index=1
c
      do iz=izmin,izmax
        z=zinit+float(iz)*dx
        r2=-eta*z**2
        valexp=exp(r2)
        expvalz(index)=valexp
        index=index+1
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine build_smooth_density(
     &    ncf,        ro1,        npwy,       npwz,       gridranges,
     &    ro2,        ncs,        sharpness,  ncspl,      dens,
     &    inx,        basdat,     x,          Lymin,      Lzmin,
     &    ro4,        ixreg,      iyregmin,   iyregmax,   griddens,
     &    iystoredim, yexpstore,  izstoredim, zexpstore,  ipexpdim,
     &    ipyexpstore,ipzexpstore)
c
c THIS SUBROUTINE BUILDS UP THE SMOOTH DENSITY IN 2D AND PUT IT
c INTO THE RO4 ARRAY.
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      integer ncf,ncs,inx(12,ncs)
      integer npwy,npwz,ixreg,iyregmin,iyregmax
      integer iystoredim,izstoredim,ipexpdim
      integer ipyexpstore(ipexpdim),ipzexpstore(ipexpdim)
      integer gridranges(6,ncs),sharpness(ncs)
      real*8 x,Lymin,Lzmin
      real*8 basdat(13,*),dens(ncf*(ncf+1)/2)
      real*8 yexpstore(iystoredim)
      real*8 zexpstore(izstoredim)
c
      real*8 ro4(npwz,npwy)
      real*8 ro1(npwz,npwy),ro2(npwz,npwy)
      real*8 denspart(10)
c
      do ics=1,ncs  ! loop over the contracted shell ics
c
        isharpness=sharpness(ics)
        if(isharpness .eq. 3) goto 150 !the function is sharp
        ixmin=gridranges(1,ics)
        ixmax=gridranges(2,ics)
        ifirstf=inx(11,ics)+1
        ilastf=inx(10,ics)
c
        if((ixmin .gt. ixreg) .or. (ixmax .lt. ixreg)) then
          goto 150 !the function is zero
        else
          ixminold=ixmin
          ixmin=ixreg
          ixmax=ixreg
        end if
c
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iyminold=iymin
        iymaxold=iymax
c
        call get_y_gridrange(iymin,iymax,ixminold,ixreg,icsrad,
     &                       icsradx2,iycenter)
c
        if(iymax-iymin .eq. 1) goto 150
        if((iymin .ge. iyregmax) .or. (iymax .le. iyregmin)) then
          goto 150
        else
          iymin=max(iymin,iyregmin)
          iymax=min(iymax,iyregmax)
        end if
c
        iyindexmin=iymin-iyregmin+1
        iyindexmax=npwy-(iyregmax-iymax)
c
        izmin=gridranges(5,ics)
        izmax=gridranges(6,ics)
c
        ibeg=inx(1,ics)+1
        Px=basdat(11,ibeg)
        Py=basdat(12,ibeg)
        Pz=basdat(13,ibeg)
        iend=inx(5,ics)
        icssize=iend-ibeg+1
        itype=inx(12,ics)
c
c
        do ifunc=ifirstf,ilastf !loop over  contracted basis function i
          ii = (ifunc*(ifunc-1))/2
          icomp=ifunc-ifirstf+1
          cji=1.0d0
c
          if(ibeg .ne. iend) then
            do iy=iyindexmin,iyindexmax
              do iz=izmin,izmax
                ro2(iz,iy)=0.0d0
              end do
            end do
          end if
c
          do ish=ibeg,iend  ! loop over the contractions
            eta=basdat(1,ish)
            if((itype .eq. 3) .and. (icomp .gt. 1)) then
              concoef=basdat(3,ish)
            else
              concoef=basdat(2,ish)
            end if
c
            call calc_expvaluesv(expvalxv,Px,x,eta,
     &                      cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
            ipyexp=ipyexpstore(ish)+(iymin-iyminold)
            ipzexp=ipzexpstore(ish)
            npwyexp=npwy-(iymin-iyminold)
            if(icssize .eq. 1) then
cc
              call build_a_gridfunction(
     &    itype,  expvalxv,   icomp,      Px,         Py,
     &    Pz,     ro2,        iyindexmin, iyindexmax, izmin,
     &    izmax,  x,          Lymin,      Lzmin,      npwy,
     &    npwz,   griddens,   xexpvalxv,  xxexpvalxv, yexpstore(ipyexp),
     &    npwyexp,zexpstore(ipzexp))
c
            else
c
              call build_a_gridfunction(
     &    itype,  expvalxv,   icomp,      Px,         Py,
     &    Pz,     ro1,        iyindexmin, iyindexmax, izmin,
     &    izmax,  x,          Lymin,      Lzmin,      npwy,
     &    npwz,   griddens,   xexpvalxv,  xxexpvalxv, yexpstore(ipyexp),
     &    npwyexp,zexpstore(ipzexp))
c
            end if
c
            if(icssize .gt. 1) then
              do iy=iyindexmin,iyindexmax
                do iz=izmin,izmax
                  ro2(iz,iy)=ro2(iz,iy)+ro1(iz,iy)
                end do
              end do
            end if
c
          end do  ! end of the loop over the contractions ish
c
c fml
          do jcs=1,ics   ! loop over the contracted shell jcs
c
            jsharpness=sharpness(jcs)
            if(jsharpness .eq. 3) goto 120 !no sharp is allowed
            jxmin=gridranges(1,jcs)
            jxmax=gridranges(2,jcs)
            jymin=gridranges(3,jcs)
            jyminold=jymin
            jymax=gridranges(4,jcs)
c new part
            call get_y_gridrange(jymin,jymax,jxmin,ixreg,jcsrad,
     &                         jcsradx2,jycenter)
c
            if(jymax-jymin .eq. 1) goto 120
            jymin=max(jymin,iyregmin)
            jymax=min(jymax,iyregmax)
c
            jzmin=gridranges(5,jcs)
            jzmax=gridranges(6,jcs)
c
c now looking for the common region
c
            ijxmin=max(ixmin,jxmin)
            ijymin=max(iymin,jymin)
            ijzmin=max(izmin,jzmin)
            ijxmax=min(ixmax,jxmax)
            ijymax=min(iymax,jymax)
            ijzmax=min(izmax,jzmax)
c
            ijyindexmin=iyindexmin+ijymin-iymin
            ijyindexmax=iyindexmax-iymax+ijymax
            jydmin=ijyindexmin-(jymin-iyregmin+1)
            iydmin=ijyindexmin-(iymin-iyregmin+1)
            jfirstf=inx(11,jcs)+1
            jlastf=inx(10,jcs)
c
            if((ijxmin .eq. ijxmax) .and. (ijymin .lt. ijymax) .and.
     &         (ijzmin .lt. ijzmax)) then
c
              jbeg=inx(1,jcs)+1
              Px=basdat(11,jbeg)
              Py=basdat(12,jbeg)
              Pz=basdat(13,jbeg)
c
              jend=inx(5,jcs)
              jcssize=jend-jbeg+1
              jtype=inx(12,jcs)
              index=1
              do jfunc=jfirstf,jlastf
                if(ifunc.ge.jfunc) then
                  ij = ii+jfunc
                else
                  ij = (jfunc*(jfunc-1))/2 + ifunc
                endif
                denspart(index)=dens(ij)
                index=index+1
              end do
c
              if(jcs .ne. ics) then
                cji=2.0d0
              else
                cji=1.0d0
              end if
c
              do jsh=jbeg,jend
c
                eta=basdat(1,jsh)
                if(jtype .eq. 3) then
                  concoef=basdat(3,jsh)
                  ccs=basdat(2,jsh)/concoef
                else
                  concoef=basdat(2,jsh)
                end if
c
                call calc_expvaluesv(expvalxv,Px,x,eta,
     &               cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
                ipyexp=ipyexpstore(jsh)+(ijymin-jyminold)
                ipzexp=ipzexpstore(jsh)+(ijzmin-jzmin)
                npwyexp=npwy-(ijymin-jyminold)
                npwzexp=npwz-(ijzmin-jzmin)
cc
                call make_gridfunctions_d(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Px,            Py,            Pz,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    x,             Lymin,          Lzmin,        npwy,
     &    npwz,          griddens,   yexpstore(ipyexp),npwyexp,
     & zexpstore(ipzexp),ro4,            ro2,          denspart,
     &    npwzexp)
c
              end do !over the contractions
c
            end if
c
 120        continue
          end do  ! end of the loop over the contracted shell jcs
c
        end do  ! end of the loop over contracted basis function i
c
 150    continue
      end do  ! end of the loop over the contracted shell ics
c
      return
      end
c***********************************************************************
c
      Subroutine get_y_gridrange(iymin,iymax,ixmin,ix,icsrad,
     &                          icsradx2,iycenter)
      implicit real*8(a-h,o-z)
c
ct      icsrad=(iymax-iymin)/2
      icsrad=(iymax-iymin)/2 + 1
      ixcenter=ixmin+icsrad
      iycenter=iymin+icsrad
      icsradx2=icsrad**2-(ix-ixcenter)**2
      if(icsradx2 .gt. 0) then
        maxyrad=int(sqrt(float(icsradx2)))+1
      else
        maxyrad=0
      end if
      iymin=max(iymin,iycenter-maxyrad)
      iymax=min(iymax,iycenter+maxyrad)
c
      return
      end
c***********************************************************************
c
      Subroutine calc_expvaluesv(expvalx,Px,Lxmin,eta,
     &                          cji,concoef,griddens,xexpvalx,xxexpvalx)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin
      data pi/3.1415926535897932384626433d0/
c
      const1=sqrt(sqrt((2.0d0 * eta / pi)**3)) * cji * concoef
c
      dx=1.0d0/griddens
      x=Lxmin - Px
c
      r2=-eta*x**2
      expvalx=const1*exp(r2) ! So this includes the const1 also !!!
      xexpvalx=x*expvalx ! So this includes the const1 also !!!
      xxexpvalx=x*xexpvalx ! So this includes the const1 also !!!
c
      return
      end
c***********************************************************************
c
      Subroutine mixed_density_help1(
     &    icspl,      icspltype,  ncspl,      cssharps,   ncomp,
     &    gridranges, ncs,        ixmin,      ixmax,      iymin,
     &    iymax,      izmin,      izmax,      inx,        ifirstf,
     &    ilastf,     icorecut,   isharpgridrange2)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      dimension inx(12,ncs)
      integer icspltype(ncspl)
      integer gridranges(6,ncs)
      integer cssharps(ncspl)
      integer isharpgridrange2(6,ncspl)
      Dimension NFunc(7)
      Data NFunc/ 1, 3, 4, 5, 6, 7,10/
c
      itype=icspltype(icspl)
      isharp=cssharps(icspl)
      ifirstf=inx(11,isharp)+1
      ilastf=inx(10,isharp)
c
      ncomp = NFunc(itype)
c
      if(icorecut .eq. 1) then
        ixmin=max(gridranges(1,isharp),isharpgridrange2(1,icspl))
        ixmax=min(gridranges(2,isharp),isharpgridrange2(2,icspl))
        iymin=max(gridranges(3,isharp),isharpgridrange2(3,icspl))
        iymax=min(gridranges(4,isharp),isharpgridrange2(4,icspl))
        izmin=max(gridranges(5,isharp),isharpgridrange2(5,icspl))
        izmax=min(gridranges(6,isharp),isharpgridrange2(6,icspl))
      else
        ixmin=gridranges(1,isharp)
        ixmax=gridranges(2,isharp)
        iymin=gridranges(3,isharp)
        iymax=gridranges(4,isharp)
        izmin=gridranges(5,isharp)
        izmax=gridranges(6,isharp)
      end if
c
      return
      end
c***********************************************************************
c
      Subroutine readin_sharps2D(isharpgrd,npwy,npwz,ro)
      IMPLICIT REAL*8(A-H,O-Z)
ct      real*8 ro(npwz,npwy)
      integer*4 ro(npwz,npwy)
c
      read(isharpgrd) ro
c
      return
      end
c***********************************************************************
c
      Subroutine make_mixed_dpf(
     &    ncf,           ro6,           npwy,          npwz,
     &    gridranges,    ro1,           ncs,           sharpness,
     &    ncspl,         nfsh,          igrdsh,        ro3,
     &    dens,          inx,           basdat,        Lxmin,
     &    Lymin,         Lzmin,         ro5,           ixreg,
     &    iyregmin,      iyregmax,      icspl,         ifunc,
     &    sharpovs,      maxovs,        nsharpovs,     npwsz,
     &    iszmin,        ro4,           Fockmx,        griddens,
     &    iystoredim,    yexpstore,     izstoredim,    zexpstore,
     &    ipexpdim,      ipyexpstore,   ipzexpstore,   comprcoef,
     &    icorecut,    isharpgridrange1,isharpgridrange2)

      use memory

c
      IMPLICIT REAL*8(A-H,O-Z)
c
c So npwx=1 in this subroutine and npwy=npwregy in this subroutine
c Actual values of the Lxmin and Lymin are coming from outside
c
      integer ncf,npwy,npwz,ncs,ncspl,nfsh,igrdsh,ixreg
      integer iyregmin,iyregmax,icspl,ifunc,maxovs,npwsz
      integer iszmin,iystoredim,izstoredim,ipexpdim,icorecut
      integer*4 ro6(npwsz,npwy) !sharp function
      real*8 ro1(npwz,npwy) !working array if neded
      real*8 ro3(npwz,npwy) !smooth function
      real*8 ro5(npwz,npwy) !mixed density
      real*8 ro4(npwz,npwy) !Coulomb potential from the smooth density
      real*8 Fockmx(ncf*(ncf+1)/2)
      real*8 fockcontr(10)  !fock mx contribution
      real*8 denspart(10)   !density mx contribution
      integer gridranges(6,ncs)
      integer sharpness(ncs)
      integer inx(12,ncs)
      real*8 basdat(13,*)
      real*8 dens(ncf*(ncf+1)/2)
      real*8 Lxmin,Lymin,Lzmin
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
c
      real*8 yexpstore(iystoredim)
      real*8 zexpstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      integer isharpgridrange1(6,ncspl)
      integer isharpgridrange2(6,ncspl)
c
c     common /big/bl(5000000)
c
      if(icorecut .eq. 1) then
        icorexmin2=isharpgridrange2(1,icspl)
        icorexmax2=isharpgridrange2(2,icspl)
        if(icorexmin2 .gt. ixreg) return !the given x grid does not
                                         !contribute to the core
        if(icorexmax2 .lt. ixreg) return !the given x grid does not
                                         !contribute to the core
        icoreymin2=isharpgridrange2(3,icspl)
        icoreymax2=isharpgridrange2(4,icspl)
        icorezmin2=isharpgridrange2(5,icspl)
        icorezmax2=isharpgridrange2(6,icspl)
c
        icorexmin1=isharpgridrange1(1,icspl)
        icorexmax1=isharpgridrange1(2,icspl)
        icoreymin1=isharpgridrange1(3,icspl)
        icoreymax1=isharpgridrange1(4,icspl)
        icorezmin1=isharpgridrange1(5,icspl)
        icorezmax1=isharpgridrange1(6,icspl)
      end if
c
      ii = (ifunc*(ifunc-1))/2
      npwzdiff=iszmin-1
      istep=1
      ip=0
      cji=2.0d0
c
      do jcssmooth=1,nsharpovs(icspl)
        jcs=sharpovs(jcssmooth,icspl) ! Loop over the contracted shells
c
        jsharpness=sharpness(jcs)
        if(jsharpness .eq. 3) goto 120 !no sharp is allowed
        jxmin=gridranges(1,jcs)
        if(jxmin .gt. ixreg) goto 120 !the given x grid does not
                                      !contribute to the function
        jxmax=gridranges(2,jcs)
        if(jxmax .lt. ixreg) goto 120 !the given x grid does not
                                      !contribute to the function
        jymin=gridranges(3,jcs)
        jyminold=jymin
        jymax=gridranges(4,jcs)
        jzmin=gridranges(5,jcs)
        jzmax=gridranges(6,jcs)
c
cc new part to calculate spherical regions
c
        call get_y_gridrange(jymin,jymax,jxmin,ixreg,jcsrad,
     &                       jcsradx2,jycenter)
c
cc
c now looking for the common region
c
        if(icorecut .eq. 1) then
          ijxmin=max(ixreg,jxmin,icorexmin2)
          ijymin=max(iyregmin,jymin,icoreymin2)
          ijzmin=max(jzmin,icorezmin2)
          ijxmax=min(ixreg,jxmax,icorexmax2)
          ijymax=min(iyregmax,jymax,icoreymax2)
          ijzmax=min(jzmax,icorezmax2)
        else
          ijxmin=max(ixreg,jxmin)
          ijymin=max(iyregmin,jymin)
          ijzmin=jzmin
          ijxmax=min(ixreg,jxmax)
          ijymax=min(iyregmax,jymax)
          ijzmax=jzmax
        end if
c
        if((ijxmin .eq. ijxmax) .and. (ijymin .lt. ijymax) .and.
     &     (ijzmin .lt. ijzmax)) then
c        write(6,*)"There is common region !"
          ijyindexmin=1+ijymin-iyregmin
          ijyindexmax=npwy-iyregmax+ijymax
c
          jfirstf=inx(11,jcs)+1
          jlastf=inx(10,jcs)
c
c Now ijxmin=ijxmax=ixreg
c
          jbeg=inx(1,jcs)+1
          Px=basdat(11,jbeg)
          Py=basdat(12,jbeg)
          Pz=basdat(13,jbeg)
c
          jend=inx(5,jcs)
          jcssize=jend-jbeg+1
          jtype=inx(12,jcs)
c
          index=1
          do jfunc=jfirstf,jlastf
            if(ifunc.ge.jfunc) then
              ij = ii + jfunc
            else
              ij = (jfunc*(jfunc-1))/2 + ifunc
            endif
            denspart(index)=dens(ij)
            index=index+1
          end do
c
          if(comprcoef .gt. 0) then
            cji_new=cji/comprcoef    !for the compression to i4
          else
            cji_new=0.0d0
          end if
c
          do jsh=jbeg,jend
c
            eta=basdat(1,jsh)
            if(jtype .eq. 3) then
              concoef=basdat(3,jsh)
              ccs=basdat(2,jsh)/concoef
            else
              concoef=basdat(2,jsh)
            end if
c
            ipyexp=ipyexpstore(jsh)+(ijymin-jyminold)
            ipzexp=ipzexpstore(jsh)+(ijzmin-jzmin)
            npwyexp=ijyindexmax-ijyindexmin+1
            npwzexp=ijzmax-ijzmin+1
c
c Here must be the new trafo to cut the cores !
c
            if(icorecut .eq. 1) then
              call calc_expvalx_corecut(
     &    expvalxv,   Px,         Lxmin,      eta,        cji_new,
     &    concoef,    griddens,   xexpvalxv,  xxexpvalxv, icorexmin1,
     &    icorexmax1, icorexmin2, icorexmax2, ixreg)
c
c    allocate arrays core1y(npwyexp) and core1z(npwzexp)
c
              call mmark
              call getmem(npwyexp,icore1y)
              call getmem(npwzexp,icore1z)
c
              call corecut_1d_2vec(
     &    bl(icore1y),bl(icore1z),      npwyexp,   yexpstore(ipyexp),
     &    npwzexp,    zexpstore(ipzexp),icoreymin1,icoreymax1,
     &    icorezmin1, icorezmax1,       ijzmin,    ijzmax,
     &    ijymin,     ijymax)
c
              call make_gridfunctions_dpf(
     &    jtype,      ccs,        expvalxv,   xexpvalxv,  xxexpvalxv,
     &    Px,         Py,         Pz,         fockcontr,  ijyindexmin,
     &    ijyindexmax,ijzmin,     ijzmax,     Lxmin,      Lymin,
     &    Lzmin,      npwy,       npwz,       griddens,   bl(icore1y),
     &    npwyexp,    bl(icore1z),ro5,        ro6,        ro4,
     &    denspart,   npwsz,      npwzdiff,   npwzexp)
c
              call retmark    !  deallocate core1y,core1z
            else
              call calc_expvaluesv(expvalxv,Px,Lxmin,eta,
     &             cji_new,concoef,griddens,xexpvalxv,xxexpvalxv)
c
              call make_gridfunctions_dpf(
     &      jtype,              ccs,expvalxv,xexpvalxv,      xxexpvalxv,
     &         Px,               Py,      Pz,fockcontr,     ijyindexmin,
     &ijyindexmax,           ijzmin,  ijzmax,    Lxmin,           Lymin,
     &      Lzmin,             npwy,    npwz,griddens,yexpstore(ipyexp),
     &    npwyexp,zexpstore(ipzexp),     ro5,      ro6,             ro4,
     &   denspart,            npwsz, npwzdiff, npwzexp)
            end if
c
            index=1
            do jfunc=jfirstf,jlastf
              value_fock=fockcontr(index)
              if(ifunc.ge.jfunc) then
                ij = ii + jfunc
                if(ifunc.eq.jfunc) value_fock = value_fock + value_fock
              else
                ij = (jfunc*(jfunc-1))/2 + ifunc
              endif
              Fockmx(ij)=Fockmx(ij) + value_fock
              index=index+1
            end do
c
          end do ! over the contraction
c
        end if
c
 120    continue
c
      end do
c
      return
      end
c***********************************************************************
c
      Subroutine calc_expvalx_corecut(
     &    expvalx,    Px,         Lxmin,      eta,        cji,
     &    concoef,    griddens,   xexpvalx,   xxexpvalx,  ixmin1,
     &    ixmax1,     ixmin2,     ixmax2,     ix)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! precalculated
                                                       !P6 polinom coefs
      data pi/3.1415926535897932384626433d0/
c
      const1=sqrt(sqrt((2.0d0 * eta / pi)**3)) * cji * concoef
      dx=1.0d0/griddens
      x=Lxmin - Px
c
c
      if((ix .ge. ixmin1) .and. (ix .le. ixmax1)) then !no cut at all
        const2=const1
      else if(ix .lt. ixmin1) then !cut on left
        xleft=dfloat(ixmin1)-str
        vixv=dfloat(ix)-xleft
        vixv3=vixv**3
        fmult=(c4l+(c3l+(c2l+c1l*vixv)*vixv)*vixv)*vixv3
        const2=fmult*const1
      else if(ix .gt. ixmax1) then !cut on right
        ixv=ix-ixmax1
        vixv=dfloat(ixv)
        vixv4=vixv**4
        fmult=1.D0+(c3r+(c2r+c1r*vixv)*vixv)*vixv4
        const2=fmult*const1
      end if
c
      r2=-eta*x**2
      expvalx=const2*exp(r2) ! So this includes the const1 also !!!
      xexpvalx=x*expvalx
      xxexpvalx=x*xexpvalx
c
      return
      end
c***********************************************************************
c
      Subroutine corecut_1d(
     &    core1y,     core2y,     core1z,     npwyexp,    yexpstore,
     &    y2expstore, npwzexp,    zexpstore,  iymin1,     iymax1,
     &    izmin1,     izmax1,     izmin2,     izmax2,      iymin3,
     &    iymax3)
c
c     This subroutine performs the smoothing on a given diffuse
c     primitive function along 1d using a six-th order polinomial
c     and put the results into 1d output arrays.
c
c    iymin1,iymax1,izmin1,izmax1     soft limits of the sharp functions
c    izmin2,izmax2,iymin3,iymax3     hard limits of the sharp functions
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      integer npwyexp,npwzexp,iymin1,iymax1,izmin1
      integer izmax1,izmin2,izmax2,iymin3,iymax3
      real*8 yexpstore(npwyexp) ! original 1d exp(-eta*y^2)
      real*8 y2expstore(npwyexp) ! original 1d y*exp(-eta*y^2)
      real*8 zexpstore(npwzexp) ! original 1d exp(-eta*z^2)
      real*8 core1y(npwyexp),core2y(npwyexp)
      real*8 core1z(npwzexp)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! precalculated
                                                       !P6 polinom coefs
c
c FOR Y
c
      if((iymin1 .le. iymin3) .and. (iymax1 .ge. iymax3)) then
        index=1             ! no cut at all
        do iy=iymin3,iymax3 ! inner region
          core1y(index)=yexpstore(index)
          core2y(index)=y2expstore(index)
          index=index+1
        end do
      else if((iymin1 .le. iymin3) .and. (iymax1 .lt. iymax3)) then
        if(iymax1 .gt. iymin3) then      !cut on right side
          index=1
          do iy=iymin3,iymax1 ! inner region
            core1y(index)=yexpstore(index)
            core2y(index)=y2expstore(index)
            index=index+1
          end do
          do iy=iymax1+1,iymax3 !right side with cut
            iyv=iy-iymax1
            viyv=dfloat(iyv)
            viyv4=viyv**4
            fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
            core1y(index)=fmult*yexpstore(index)
            core2y(index)=fmult*y2expstore(index)
            index=index+1
          end do
        else
          index=1
          do iy=iymin3,iymax3 !right side with cut
            iyv=iy-iymax1
            viyv=dfloat(iyv)
            viyv4=viyv**4
            fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
            core1y(index)=fmult*yexpstore(index)
            core2y(index)=fmult*y2expstore(index)
            index=index+1
          end do
        end if
      else if((iymin1 .gt. iymin3) .and. (iymax1 .ge. iymax3)) then
                                                 !cut on left side
        if(iymin1 .lt. iymax3) then !ordinary case
          index=1
          yleft=dfloat(iymin1)-str
          do iy=iymin3,iymin1 !left side cut
            viyv=dfloat(iy)-yleft
            viyv3=viyv**3
            fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
            core1y(index)=fmult*yexpstore(index)
          core2y(index)=fmult*y2expstore(index)
          index=index+1
          end do
          do iy=iymin1+1,iymax3 ! inner region
            core1y(index)=yexpstore(index)
            core2y(index)=y2expstore(index)
            index=index+1
          end do
        else
          index=1
          yleft=dfloat(iymin1)-str
          do iy=iymin3,iymax3 !left side cut
            viyv=dfloat(iy)-yleft
            viyv3=viyv**3
            fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
            core1y(index)=fmult*yexpstore(index)
            core2y(index)=fmult*y2expstore(index)
            index=index+1
          end do
        end if
      else if((iymin1 .gt. iymin3) .and. (iymax1 .lt. iymax3)) then
                                                !cut on both sides
        index=1
        yleft=dfloat(iymin1)-str
        do iy=iymin3,iymin1 !left side cut
          viyv=dfloat(iy)-yleft
          viyv3=viyv**3
          fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
          core1y(index)=fmult*yexpstore(index)
          core2y(index)=fmult*y2expstore(index)
          index=index+1
        end do
        do iy=iymin1+1,iymax1 ! inner region
          core1y(index)=yexpstore(index)
          core2y(index)=y2expstore(index)
          index=index+1
        end do
        do iy=iymax1+1,iymax3 !right side with cut
          iyv=iy-iymax1
          viyv=dfloat(iyv)
          viyv4=viyv**4
          fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
          core1y(index)=fmult*yexpstore(index)
          core2y(index)=fmult*y2expstore(index)
          index=index+1
        end do
      end if
c
c
c FOR Z
c
      if((izmin1 .le. izmin2) .and. (izmax1 .ge. izmax2)) then
        index=1                                !no cut at all
        do iz=izmin2,izmax2 ! inner region
          core1z(index)=zexpstore(index)
          index=index+1
        end do
      else if((izmin1 .le. izmin2) .and. (izmax1 .lt. izmax2)) then
                                               !cut on right side
        if(izmax1 .gt. izmin2) then
          index=1
          do iz=izmin2,izmax1 ! inner region
            core1z(index)=zexpstore(index)
            index=index+1
          end do
          do iz=izmax1+1,izmax2 !right side with cut
            izv=iz-izmax1
            vizv=dfloat(izv)
            vizv4=vizv**4
            fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        else
          index=1
          do iz=izmin2,izmax2 !right side with cut
            izv=iz-izmax1
            vizv=dfloat(izv)
            vizv4=vizv**4
            fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        end if
      else if((izmin1 .gt. izmin2) .and. (izmax1 .ge. izmax2)) then
                                                   !cut on left side
        if(izmin1 .lt. izmax2) then !ordinary case
          index=1
          zleft=dfloat(izmin1)-str
          do iz=izmin2,izmin1 !left side cut
            vizv=dfloat(iz)-zleft
            vizv3=vizv**3
            fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
          do iz=izmin1+1,izmax2 ! inner region
            core1z(index)=zexpstore(index)
            index=index+1
          end do
        else
          index=1
          zleft=dfloat(izmin1)-str
          do iz=izmin2,izmax2 !left side cut
            vizv=dfloat(iz)-zleft
            vizv3=vizv**3
            fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        end if
      else if((izmin1 .gt. izmin2) .and. (izmax1 .lt. izmax2)) then
                                                !cut on both sides
        index=1
        zleft=dfloat(izmin1)-str
        do iz=izmin2,izmin1 !left side cut
          vizv=dfloat(iz)-zleft
          vizv3=vizv**3
          fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
          core1z(index)=fmult*zexpstore(index)
          index=index+1
        end do
        do iz=izmin1+1,izmax1 ! inner region
          core1z(index)=zexpstore(index)
          index=index+1
        end do
        do iz=izmax1+1,izmax2 !right side with cut
          izv=iz-izmax1
          vizv=dfloat(izv)
          vizv4=vizv**4
          fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
          core1z(index)=fmult*zexpstore(index)
          index=index+1
        end do
      end if
c
      return
      end
c***********************************************************************
c
      Subroutine corecut_1d_2vec(
     &    core1y,     core1z,     npwyexp,    yexpstore,  npwzexp,
     &    zexpstore,  iymin1,     iymax1,     izmin1,     izmax1,
     &    izmin2,     izmax2,     iymin3,     iymax3)
c
c     This subroutine performs the smoothing on a given diffuse
c     primitive function along 1d using a six-th order polinomial
c     and put the results into 1d output arrays.
c
c    iymin1,iymax1,izmin1,izmax1     soft limits of the sharp functions
c    izmin2,izmax2,iymin3,iymax3     hard limits of the sharp functions
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      integer npwyexp,npwzexp,iymin1,iymax1,izmin1
      integer izmax1,izmin2,izmax2,iymin3,iymax3
      real*8 yexpstore(npwyexp) ! original 1d exp(-eta*y^2)
      real*8 zexpstore(npwzexp) ! original 1d exp(-eta*z^2)
      real*8 core1y(npwyexp)
      real*8 core1z(npwzexp)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! precalculated
                                                       !P6 polinom coefs
c
c FOR Y
c
      if((iymin1 .le. iymin3) .and. (iymax1 .ge. iymax3)) then
        index=1                                      !no cut at all
        do iy=iymin3,iymax3 ! inner region
          core1y(index)=yexpstore(index)
          index=index+1
        end do
      else if((iymin1 .le. iymin3) .and. (iymax1 .lt. iymax3)) then
                                                     !cut on right side
        if(iymax1 .gt. iymin3) then
          index=1
          do iy=iymin3,iymax1 ! inner region
            core1y(index)=yexpstore(index)
            index=index+1
          end do
          do iy=iymax1+1,iymax3 !right side with cut
            iyv=iy-iymax1
            viyv=dfloat(iyv)
            viyv4=viyv**4
            fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
            core1y(index)=fmult*yexpstore(index)
            index=index+1
          end do
        else
          index=1
          do iy=iymin3,iymax3 !right side with cut
            iyv=iy-iymax1
            viyv=dfloat(iyv)
            viyv4=viyv**4
            fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
            core1y(index)=fmult*yexpstore(index)
            index=index+1
          end do
        end if
      else if((iymin1 .gt. iymin3) .and. (iymax1 .ge. iymax3)) then
                                                      !cut on left side
        if(iymin1 .lt. iymax3) then !ordinary case
          index=1
          yleft=dfloat(iymin1)-str
          do iy=iymin3,iymin1 !left side cut
            viyv=dfloat(iy)-yleft
            viyv3=viyv**3
            fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
            core1y(index)=fmult*yexpstore(index)
            index=index+1
          end do
          do iy=iymin1+1,iymax3 ! inner region
            core1y(index)=yexpstore(index)
            index=index+1
          end do
        else
          index=1
          yleft=dfloat(iymin1)-str
          do iy=iymin3,iymax3 !left side cut
            viyv=dfloat(iy)-yleft
            viyv3=viyv**3
            fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
            core1y(index)=fmult*yexpstore(index)
            index=index+1
          end do
        end if
      else if((iymin1 .gt. iymin3) .and. (iymax1 .lt. iymax3)) then
                                                 !cut on both sides
        index=1
        yleft=dfloat(iymin1)-str
        do iy=iymin3,iymin1 !left side cut
          viyv=dfloat(iy)-yleft
          viyv3=viyv**3
          fmult=(c4l+(c3l+(c2l+c1l*viyv)*viyv)*viyv)*viyv3
          core1y(index)=fmult*yexpstore(index)
          index=index+1
        end do
        do iy=iymin1+1,iymax1 ! inner region
          core1y(index)=yexpstore(index)
          index=index+1
        end do
        do iy=iymax1+1,iymax3 !right side with cut
          iyv=iy-iymax1
          viyv=dfloat(iyv)
          viyv4=viyv**4
          fmult=1.D0+(c3r+(c2r+c1r*viyv)*viyv)*viyv4
          core1y(index)=fmult*yexpstore(index)
          index=index+1
        end do
      end if
c
c
c FOR Z
c
      if((izmin1 .le. izmin2) .and. (izmax1 .ge. izmax2)) then
        index=1                                       !no cut at all
        do iz=izmin2,izmax2 ! inner region
          core1z(index)=zexpstore(index)
          index=index+1
        end do
      else if((izmin1 .le. izmin2) .and. (izmax1 .lt. izmax2)) then
                                                    !cut on right side
        if(izmax1 .gt. izmin2) then
          index=1
          do iz=izmin2,izmax1 ! inner region
            core1z(index)=zexpstore(index)
            index=index+1
          end do
          do iz=izmax1+1,izmax2 !right side with cut
            izv=iz-izmax1
            vizv=dfloat(izv)
            vizv4=vizv**4
            fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        else
          index=1
          do iz=izmin2,izmax2 !right side with cut
            izv=iz-izmax1
            vizv=dfloat(izv)
            vizv4=vizv**4
            fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        end if
      else if((izmin1 .gt. izmin2) .and. (izmax1 .ge. izmax2)) then
                                                     !cut on left side
        if(izmin1 .lt. izmax2) then !ordinary case
          index=1
          zleft=dfloat(izmin1)-str
          do iz=izmin2,izmin1 !left side cut
            vizv=dfloat(iz)-zleft
            vizv3=vizv**3
            fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
          do iz=izmin1+1,izmax2 ! inner region
            core1z(index)=zexpstore(index)
            index=index+1
          end do
        else
          index=1
          zleft=dfloat(izmin1)-str
          do iz=izmin2,izmax2 !left side cut
            vizv=dfloat(iz)-zleft
            vizv3=vizv**3
            fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
            core1z(index)=fmult*zexpstore(index)
            index=index+1
          end do
        end if
      else if((izmin1 .gt. izmin2) .and. (izmax1 .lt. izmax2)) then
                                                      !cut on both side
        index=1
        zleft=dfloat(izmin1)-str
        do iz=izmin2,izmin1 !left side cut
          vizv=dfloat(iz)-zleft
          vizv3=vizv**3
          fmult=(c4l+(c3l+(c2l+c1l*vizv)*vizv)*vizv)*vizv3
          core1z(index)=fmult*zexpstore(index)
          index=index+1
        end do
        do iz=izmin1+1,izmax1 ! inner region
          core1z(index)=zexpstore(index)
          index=index+1
        end do
        do iz=izmax1+1,izmax2 !right side with cut
          izv=iz-izmax1
          vizv=dfloat(izv)
          vizv4=vizv**4
          fmult=1.D0+(c3r+(c2r+c1r*vizv)*vizv)*vizv4
          core1z(index)=fmult*zexpstore(index)
          index=index+1
        end do
      end if
c
      return
      end
c***********************************************************************
c
      Subroutine calc_FTC_fock_smooth(
     &    ncf,        ro1,        npwx,       npwy,
     &    npwz,       gridranges, ro2,        ncs,
     &    sharpness,  ncspl,      igrdsh,
     &    ro3,        inx,        basdat,     Lxmin,
     &    Lymin,      Lzmin,      ro4,        ixreg,
     &    iyregmin,   iyregmax,   Fockmx,     griddens,
     &    iystoredim, yexpstore,  izstoredim, zexpstore,
     &    ipexpdim,   ipyexpstore,ipzexpstore)
      IMPLICIT REAL*8(A-H,O-Z)
c
c THIS SUBROUTINE CALCULATES THE 2D COULOMB MATRIX CONTRIBUTION
c FOR EVERY g(i),g(j) SMOOTH PAIRS USING THE PRECALCULATED SMOOTH+MIX
c COULOMB POTENTIAL.
C
c
      integer ncf,ncs,inx(12,ncs)
      integer npwy,npwz,ixreg,iyregmin,iyregmax !grid dims and intervals
      integer iystoredim,izstoredim,ipexpdim
      integer ipyexpstore(ipexpdim),ipzexpstore(ipexpdim)!pointer arrays
                                                         !for preexp
                                                         !functions
      integer gridranges(6,ncs),sharpness(ncs)
      real*8 Lxmin,Lymin,Lzmin
      real*8 basdat(13,*)
      real*8 yexpstore(iystoredim) !Precalc. exp functions
      real*8 zexpstore(izstoredim)
c
      real*8 ro4(npwz,npwy) ! Coulomb pot of smooth+mix density
      real*8 Fockmx(ncf*(ncf+1)/2)
      real*8 fockcontr(10) ! Local working array for Fock mx contribution
      real*8 ro1(npwz,npwy),ro2(npwz,npwy),ro3(npwz,npwy)!working arrays
c
      do ics=1,ncs  ! loop over the contracted shell ics
c
        isharpness=sharpness(ics)
        if(isharpness .eq. 3) goto 150 !the function is sharp
        ixmin=gridranges(1,ics)
        ixmax=gridranges(2,ics)
        ifirstf=inx(11,ics)+1
        ilastf=inx(10,ics)
c
        if((ixmin .gt. ixreg) .or. (ixmax .lt. ixreg)) then
          goto 150  !the function is zero
        else
          ixminold=ixmin
          ixmin=ixreg
          ixmax=ixreg
        end if
c
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iyminold=iymin
        iymaxold=iymax
        call get_y_gridrange(iymin,iymax,ixminold,ixreg,icsrad,
     &                       icsradx2,iycenter)
c
        if((iymin .gt. iyregmax) .or. (iymax .lt. iyregmin)) then
          goto 150
        else
          iymin=max(iymin,iyregmin)
          iymax=min(iymax,iyregmax)
        end if
        iyindexmin=iymin-iyregmin+1
        iyindexmax=npwy-(iyregmax-iymax)
c
        izmin=gridranges(5,ics)
        izmax=gridranges(6,ics)
c
cc
c      call fillin_y_storage(ics,ncs,ipystorage,iystorage,izranges,
c     &                   iyindexmin,iyindexmax,iycenter,icsrad,icsradx2,
c     &                   izmin,izmax,ip)
cc
        ibeg=inx(1,ics)+1
        Px=basdat(11,ibeg)
        Py=basdat(12,ibeg)
        Pz=basdat(13,ibeg)
        iend=inx(5,ics)
        icssize=iend-ibeg+1
        itype=inx(12,ics)
c
c
        do ifunc=ifirstf,ilastf !loop over the contracted basis function i
          ii = (ifunc*(ifunc-1))/2
          icomp=ifunc-ifirstf+1
          cji=1.0d0
c
          if(ibeg .ne. iend) then
            do iy=iyindexmin,iyindexmax
              do iz=izmin,izmax
c
                ro2(iz,iy)=0.0d0
c
              end do
            end do
          end if
c
          do ish=ibeg,iend  ! loop over the contractions
            eta=basdat(1,ish)
            if((itype .eq. 3) .and. (icomp .gt. 1)) then
              concoef=basdat(3,ish)
            else
              concoef=basdat(2,ish)
            end if
cc
c
            call calc_expvaluesv(expvalxv,Px,Lxmin,eta,
     &                      cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
            ipyexp=ipyexpstore(ish)+(iymin-iyminold)
            ipzexp=ipzexpstore(ish)
            npwyexp=npwy-(iymin-iyminold)
            if(icssize .eq. 1) then
cc
              call build_a_gridfunction(
     &    itype,      expvalxv,icomp,     Px,        Py,
     &    Pz,         ro2,     iyindexmin,iyindexmax,izmin,
     &    izmax,      Lxmin,   Lymin,     Lzmin,     npwy,
     &    npwz,       griddens,xexpvalxv, xxexpvalxv,yexpstore(ipyexp),
     &    npwyexp,    zexpstore(ipzexp))
cc
            else
cc
              call build_a_gridfunction(
     &    itype,      expvalxv,icomp,     Px,        Py,
     &    Pz,         ro1,     iyindexmin,iyindexmax,izmin,
     &    izmax,      Lxmin,   Lymin,     Lzmin,     npwy,
     &    npwz,       griddens,xexpvalxv, xxexpvalxv,yexpstore(ipyexp),
     &    npwyexp,    zexpstore(ipzexp))
cc
            end if
c
            if(icssize .gt. 1) then
              do iy=iyindexmin,iyindexmax
                do iz=izmin,izmax
c
                  ro2(iz,iy)=ro2(iz,iy)+ro1(iz,iy)
c
                end do
              end do
            end if
c
          end do  ! end of the loop over the contractions ish
c
c fml
          do jcs=1,ics   ! loop over the contracted shell jcs
c
            jsharpness=sharpness(jcs)
            if(jsharpness .eq. 3) goto 120 !no sharp is allowed
            jxmin=gridranges(1,jcs)
            jxmax=gridranges(2,jcs)
            jymin=gridranges(3,jcs)
            jyminold=jymin
            jymax=gridranges(4,jcs)
c
            call get_y_gridrange(jymin,jymax,jxmin,ixreg,jcsrad,
     &                     jcsradx2,jycenter)
c
            jymin=max(jymin,iyregmin)
            jymax=min(jymax,iyregmax)
c
            jzmin=gridranges(5,jcs)
            jzmax=gridranges(6,jcs)
c
c now looking for the common region
c
            ijxmin=max(ixmin,jxmin)
            ijymin=max(iymin,jymin)
            ijzmin=max(izmin,jzmin)
            ijxmax=min(ixmax,jxmax)
            ijymax=min(iymax,jymax)
            ijzmax=min(izmax,jzmax)
c
            ijyindexmin=iyindexmin+ijymin-iymin
            ijyindexmax=iyindexmax-iymax+ijymax
            jydmin=ijyindexmin-(jymin-iyregmin+1)
            iydmin=ijyindexmin-(iymin-iyregmin+1)
c
            jfirstf=inx(11,jcs)+1
            jlastf=inx(10,jcs)
c
            if((ijxmin .eq. ijxmax) .and. (ijymin .lt. ijymax) .and.
     &        (ijzmin .lt. ijzmax)) then
c
              jbeg=inx(1,jcs)+1
              Px=basdat(11,jbeg)
              Py=basdat(12,jbeg)
              Pz=basdat(13,jbeg)
c
              jend=inx(5,jcs)
              jcssize=jend-jbeg+1
              jtype=inx(12,jcs)
c
              if(jcs .ne. ics) then
                cji=2.0d0
              else
                cji=1.0d0
              end if
c
              do jsh=jbeg,jend
                eta=basdat(1,jsh)
                if(jtype .eq. 3) then
                  concoef=basdat(3,jsh)
                  ccs=basdat(2,jsh)/concoef
                else
                  concoef=basdat(2,jsh)
                end if
c
c
                call calc_expvaluesv(expvalxv,Px,Lxmin,eta,
     &               cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
                ipyexp=ipyexpstore(jsh)+(ijymin-jyminold)
                ipzexp=ipzexpstore(jsh)+(ijzmin-jzmin)
                npwyexp=npwy-(ijymin-jyminold)
                npwzexp=npwz-(ijzmin-jzmin)
cc
                call make_gridfunctions_f(
     &    jtype,      ccs,        expvalxv,xexpvalxv, xxexpvalxv,
     &    Px,         Py,         Pz,      fockcontr, ijyindexmin,
     &    ijyindexmax,ijzmin,     ijzmax,  Lxmin,     Lymin,
     &    Lzmin,      npwy,       npwz,    griddens,  yexpstore(ipyexp),
     &    npwyexp,    zexpstore(ipzexp),   ro2,       ro4,
     &    npwzexp)
c
                index=1
                do jfunc=jfirstf,jlastf
                  value_fock=fockcontr(index)
                  if(ifunc .ge. jfunc) then
                    ij = ii + jfunc
                  else
                    ij = (jfunc*(jfunc-1))/2 + ifunc
                  endif
                  Fockmx(ij) = Fockmx(ij) + value_Fock
                  index=index+1
                end do
              end do !over the jsh contraction
c
            end if
c
 120        continue
          end do  ! end of the loop over the contracted shell jcs
c
        end do  ! end of the loop over the contracted basis function i
c
 150    continue
      end do  ! end of the loop over the contracted shell ics
c
      return
      end
c***********************************************************************
c
      SUBROUTINE Make_FockL(ncf,scal,Fockl)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Scales FTC Fock matrix ready for DFT
C
C  ARGUMENTS
C
C  ncf     -  dimension of Fock matrix (number of basis functions)
C  scal    -  root scaling factor
C  Fockl   -  on exit lower triangle matrix ready for DFT
C
      Dimension Fockl(*)
c
      scalH = 0.5d0*scal
c
      IJ = 0
      DO 20 I=1,ncf
        DO 10 J=1,I
          IJ = IJ+1
          Fockl(IJ) = scalH*Fockl(IJ)
 10     CONTINUE
        Fockl(IJ) = Fockl(IJ)+Fockl(IJ)
 20   CONTINUE
C
      RETURN
      END
c*************************************************************************
c
c
      SUBROUTINE Make_FockL_tmp(ncf,scal,Fockl,Fockmp,densmx)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Scales FTC Fock matrix ready for DFT
C
C  ARGUMENTS
C
C  ncf     -  dimension of Fock matrix (number of basis functions)
C  scal    -  root scaling factor
C  Fockl   -  on exit lower triangle matrix ready for DFT
C
      Dimension Fockl(*)
      real*8 Fockmp(ncf,ncf)
      real*8 densmx(ncf,ncf)
c
      scalH = 0.5d0*scal
c      scalmp=2.0d-10
      scalmp=2.0d0
c
c      write(6,*)"The Multi-pole, and Density mx:"
c
      IJ = 0
      DO 20 I=1,ncf
        DO 10 J=1,I
          IJ = IJ+1
c      write(6,100) I,J,scalmp*Fockmp(J,I)*1.0d10,densmx(J,I)
          Fockl(IJ) = scalH*Fockl(IJ)+scalmp*Fockmp(J,I)
 10     CONTINUE
        Fockl(IJ) = Fockl(IJ)+Fockl(IJ)
 20   CONTINUE
C
      call f_lush(6)
 100  format(2i3,2f18.8)
      RETURN
      END
c***********************************************************************
c

      Subroutine build_a_gridfunction(
     &    itype,      expvalx,    icomp,      Px,         Py,
     &    Pz,         ro,         iymin,      iymax,      izmin,
     &    izmax,      Lxmin,      Lymin,      Lzmin,      npwy,
     &    npwz,       griddens,   xexpvalx,   xxexpvalx,  yexpstore,
     &    npwyexp,    zexpstore)
c
c This subroutine calculate the values of a given basis function
c on the Fourier grid along the whole range of interest and put them
c into the ro array.
c
c Some additional usefull information:
c
c P: The order is x,y,z and 2*sqrt(eta) is inside the concoef.
c
c D: The order is: (3z**2-r**2),(x**2-y**2),xy,xz,yz.
c    4*eta inside the concoef
c
c D6: The order is: x**2,y**2,z**2,xy,xz,yz. The first three are
c     normalized to 3, the others to 1. 4*eta inside the concoef.
c
c F:  The order is: (5x**2*y-r**2*y),(5x**2*z-r**2*z),(5y**2*x-r**2*x),
c     (5y**2*z-r**2*z),(5z**2*x-r**2*x),(5z**2*y-r**2*y),x*y*z.
c     {(2*eta)^(3/2)}/sqrt(5) inside the concoef.
c
c F10: The order is: x**3,  yx**2, zx**2, xy**2, xyz,
c                    xz**2, y**3,  zy**2, yz**2, z**3

c      The x**3,y**3,z**3 are normalized to 15/24.
c      The yx**2,zx**2,xy**2,xz**2,zy**2,yz**2 are
c                             normalized to 3/24
c      The xyz is normalized to 1/24.
c     {(2*eta)^(3/2)}/sqrt(3) inside the concoef.
c
c Input:
c  itype:        type of the shell
c  eta:          orbital exponent
c  icomp:        component of the shell
c  concoef:      contraction coef.
c  cji:          Mo coeff (if it is needed, set to 1.0d0 in the calling
c                          subroutine otherwise)
c  Px,Py,Px      center of the given basis function
c  npw:          Number of plane-wave in each direction
c  Lbox:         Box length for each direction
c  ro:           empty array with npw*npw*npw dimension
c
c Output:
c  ro:           The result array.
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 ro(npwz,npwy)
      real*8 Lxmin,Lymin,Lzmin
c
      real*8 yexpstore(npwyexp)
      real*8 zexpstore(npwz)
c
      data pi/3.1415926535897932384626433d0/
c
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
      if (itype .eq. 1) then
c
c  s type function
c
c
        valx=expvalx
        do iy=iymin,iymax
c
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 1)) then
c
c Px function
c
        valx=xexpvalx
        do iy=iymin,iymax
c
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 2)) then
c
c Py function
c
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 3)) then
c
c Pz function
c
c
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 1)) then
c
c s part of L
c
        valx=expvalx
        do iy=iymin,iymax
c
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
c
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 2)) then
c
c Px part of L
c
        valx=xexpvalx
        do iy=iymin,iymax
c
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 3)) then
c
c Py part of L
c
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 3) .and. (icomp .eq. 4)) then
c
c Pz part of L
c
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 1)) then
c
c First component of D which is (3*z^2-r^2)
c
        const2=1.0d0/dsqrt(12.0d0)
c
        valx=const2*expvalx
        r2x=x**2
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          r2xy=r2x+y**2
          indexz=1
          z=z0
          do iz=izmin,izmax
c
            a = 2.0d0*z**2 - r2xy
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 2)) then
c
c 2-nd component of D which is (x^2-y^2)
c
        const2=0.5d0
        ax2=x**2
        valx=const2*expvalx
        do iy=iymin,iymax
          a=ax2-y**2
          valxy=a*valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 4) .and. (icomp .eq. 3)) then
c
c 3-th component of D which is x*y
c
        valx=xexpvalx
        do iy=iymin,iymax
          valxy=valx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 4) .and. (icomp .eq. 4)) then
c
c 4-th component of D which is x*z
c
        do iy=iymin,iymax
          valxy=xexpvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
c
      else if ((itype .eq. 4) .and. (icomp .eq. 5)) then
c
c 5-th component of D which is y*z
c
        do iy=iymin,iymax
          valxy=expvalx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 5) .and. (icomp .eq. 1)) then
c
c First component of D6 which is x^2 (normalized to 3 in the PQS)
c
c
        do iy=iymin,iymax
          valxy=xxexpvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
c
      else if ((itype .eq. 5) .and. (icomp .eq. 2)) then
c
c Second component of D6 which is y^2 (normalized to 3 in the PQS)
c
c
        do iy=iymin,iymax
          y2=y**2
          valxy=expvalx*y2*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 5) .and. (icomp .eq. 3)) then
c
c 3-th component of D6 which is z^2 (normalized to 3 in the PQS)
c
        do iy=iymin,iymax
          valxy=expvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
c
            z2=z**2
            ro(iz,iy)=valxy*z2*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 4)) then
c
c 4-th component of D6 which is x*y
c
        do iy=iymin,iymax
          valxy=xexpvalx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 5)) then
c
c 5-th component of D6 which is x*z
c
        do iy=iymin,iymax
          valxy=xexpvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 6)) then
c
c 6-th component of D6 which is y*z
c
        do iy=iymin,iymax
          valxy=expvalx*y*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
c
            ro(iz,iy)=valxy*z*zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 6) .and. (icomp .eq. 1)) then
c
c first component of F which is 5*x^2*y-r^2*y
c
        ax=4.0d0*x**2
        valx=expvalx
        do iy=iymin,iymax
          axy=ax*y-y**3
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy-y*z**2
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 2)) then
c
c Second component of F which is 5*x^2*z-r^2*z
c
        ax=4.0d0*x**2
        valx=expvalx
        do iy=iymin,iymax
          axy=ax-y**2
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy*z-z**3
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 3)) then
c
c 3-th component of F which is 5*y^2*x-r^2*x
c
c
        ax=x**3
        valx=expvalx
        do iy=iymin,iymax
          axy=4.0d0*x*y**2-ax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy-x*z**2
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 6) .and. (icomp .eq. 4)) then
c
c 4-th component of F which is 5*y^2*z-r^2*z
c
c
        ax=x**2
        valx=expvalx
        do iy=iymin,iymax
          axy=4.0d0*y**2-ax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy*z-z**3
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 5)) then
c
c 5-th component of F which is 5*z^2*x-r^2*x
c
        ax=x**3
        valx=expvalx
        do iy=iymin,iymax
          axy=x*y**2+ax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=4.0d0*x*z**2-axy
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 6)) then
c
c 6-th component of F which is 5*z^2*y-r^2*y
c
        ax=x**2
        valx=expvalx
        do iy=iymin,iymax
          axy=ax*y+y**3
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=4.0d0*y*z**2-axy
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 6) .and. (icomp .eq. 7)) then
c
c 7-th component of F which is x*y*z
c
        const2=sqrt(40.0d0)
c
        valx=const2*expvalx
        do iy=iymin,iymax
          axy=x*y
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy*z
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
c
      else if ((itype .eq. 7) .and. (icomp .eq. 1)) then
c
c first component of F10 which is x**3
c
        x3=x**3
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= x3 * valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 2)) then
c
c Second component of F10 which is y*x**2
c
        x2=x**2
        valx=expvalx
        do iy=iymin,iymax
          axy=x2*y
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= axy * valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 3)) then
c
c Third component of F10 which is z*x**2
c
        x2=x**2
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=z*x2
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
c
      else if ((itype .eq. 7) .and. (icomp .eq. 4)) then
c
c 4-th component of F10 which is x*y**2
c
        valx=expvalx
        do iy=iymin,iymax
          axy=x*y**2
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= axy * valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 5)) then
c
c 5-th component of F10 which is x*y*z
c
        valx=expvalx
        do iy=iymin,iymax
          axy=x*y
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=axy*z
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 6)) then
c
c 6-th component of F10 which is x*z**2
c
        ax=x
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=ax*z**2
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          end do
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 7)) then
c
c 7-th component of F10 which is y**3
c
        valx=expvalx
        do iy=iymin,iymax
          ay=y**3
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            ro(iz,iy)= ay * valxy * zexpstore(indexz)
            indexz=indexz+1
c
          end do
          y=y+dx
        end do
c
      else if ((itype.eq.7) .and. (icomp.eq.8)) then
c
c 8-th component of F10 which is z*y**2
c
        valx=expvalx
        do iy=iymin,iymax
          ay=y**2
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=ay*z
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          enddo
          y=y+dx
        enddo
c
      else if ((itype.eq.7) .and. (icomp.eq.9)) then
c
c 9-th component of F10 which is y*z**2
c
        valx=expvalx
        do iy=iymin,iymax
          ay=y
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=ay*z**2
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          enddo
          y=y+dx
        enddo
c
      else if ((itype.eq.7) .and. (icomp.eq.10)) then
c
c 10-th component of F10 which is z**3
c
        valx=expvalx
        do iy=iymin,iymax
          valxy=valx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          z=z0
          do iz=izmin,izmax
            a=z**3
c
            ro(iz,iy)= a * valxy * zexpstore(indexz)
            indexz=indexz+1
            z=z+dx
c
          enddo
        enddo
c
c
      else
        call nerror(1,'build_a_gridfunction',
     $    ' FTC: Can only handle S, P, D and F functions',0,0)
      end if
c
      return
      end
c***************************************************************************
c
      Subroutine make_gridfunctions_dpf(
     &    itype,      ccs,        expvalx,    xexpvalx,   xxexpvalx,
     &    Px,         Py,         Pz,         fockcontr,  iymin,
     &    iymax,      izmin,      izmax,      Lxmin,      Lymin,
     &    Lzmin,      npwy,       npwz,       griddens,   yexpstore,
     &    npwyexp,    zexpstore,  ro5,        ro6,        ro4,
     &    densmx,     npwsz,      npwzdiff,   npwzexp)
c
c  This subroutine forms the given smooth function in 2d and builds
c  up the mixed density contribution + calculates the mixed Coulomb
c  mx. contr. using the smooth potential. Works up to f7 right now.
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      integer itype,iymin,iymax,izmin,izmax,npwy,npwz
      integer npwyexp,npwsz,npwzdiff,npwzexp
      real*8 ccs,xexpvalx,xxexpvalx,Px,Py,Pz
      real*8 Lxmin,Lymin,Lzmin,griddens
      integer*4 ro6(npwsz,npwy) ! sharp grid , input
      real*8 ro4(npwz,npwy)     ! smooth potential, input
c
      real*8 yexpstore(npwyexp)
      real*8 zexpstore(npwzexp)
      real*8 densmx(10)     ! given denity mx elements, input
      real*8 ro5(npwz,npwy) ! mixed density, input and output
      real*8 fockcontr(10)  ! fock mx contribution, output
c
      data pi/3.1415926535897932384626433d0/
c
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
c  S function
c
        d1=densmx(1)
        do iy=iymin,iymax
c
          valxy=expvalx*yexpstore(indexy)
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
c
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            fsm=valxy*zexpstore(indexz)
            indexz=indexz+1
            fctfsm=fc*fsm
            ro5(iz,iy)=ro5(iz,iy)+d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c
          end do
        end do
        fockcontr(1)=focksumm1
c
      else if (itype .eq. 2) then
c
c  P function
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          valxypx=xexpvalx*expy
          valxypy=expvalx*y*expy
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            expz=zexpstore(indexz)
c
            fsm=valxypx*expz
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now px is done
            fsm=valxypy*expz
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now py is done
            fsm=z*valxypz*expz
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3
            indexz=indexz+1
            z=z+dx
c
          end do
          y=y+dx
        end do
        fockcontr(1)=focksumm1
        fockcontr(2)=focksumm2
        fockcontr(3)=focksumm3
c
      else if (itype .eq. 3) then
c
c  SP (L) function
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
          valxypy=expvalx*y*expy
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
          z=z0
c
          do iz=izmin,izmax
c
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            expz=zexpstore(indexz)
c
            fsm=valxys*expz
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now s is done
            fsm=valxypx*expz
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now px is done
            fsm=valxypy*expz
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now py is done
            fsm=z*valxypz*expz
            fctfsm=fc*fsm
            denssumm4=d4*fctfsm
            focksumm4=focksumm4+pot*fctfsm
c now pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+
     &                            denssumm3+denssumm4
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
c
      else if (itype .eq. 4) then
c
c  D5 function
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
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            a = 2.0d0*z**2 - r2xy
c
            fsm1=zexpstore(indexz)
            fsm=a*valxy1*fsm1
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now D1 is done
            fsm=valxy2*fsm1
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now D2 is done
            fsm=valxy3*fsm1
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now D3 is done
            fsm=z*valxy4*fsm1
            fctfsm=fc*fsm
            denssumm4=d4*fctfsm
            focksumm4=focksumm4+pot*fctfsm
c now D4 is done
            fsm=z*valxy5*fsm1
            fctfsm=fc*fsm
            denssumm5=d5*fctfsm
            focksumm5=focksumm5+pot*fctfsm
c now D5 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5
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
c
      else if (itype .eq. 5) then
c
c  D6 function
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
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            fsm1=zexpstore(indexz)
            ztfsm1=z*fsm1
c
            fsm=valxy1*fsm1
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now D61 is done
            fsm=valxy2*fsm1
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now D62 is done
            fsm=valxy3*ztfsm1*z
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now D63 is done
            fsm=valxy4*fsm1
            fctfsm=fc*fsm
            denssumm4=d4*fctfsm
            focksumm4=focksumm4+pot*fctfsm
c now D64 is done
            fsm=valxy5*ztfsm1
            fctfsm=fc*fsm
            denssumm5=d5*fctfsm
            focksumm5=focksumm5+pot*fctfsm
c now D65 is done
            fsm=valxy6*ztfsm1
            fctfsm=fc*fsm
            denssumm6=d6*fctfsm
            focksumm6=focksumm6+pot*fctfsm
c now D66 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5+denssumm6
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
c  F7 function
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
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
            z2=z**2
            z3=z2*z
            a1=axy1-y*z2
            a2=axy2*z-z3
            a3=axy3-x*z2
            a4=axy4*z-z3
c     a5v=fourtx*z2
            a5=fourtx*z2-axy5
            a6=fourty*z2-axy6 !must be 4*y*z^2-x^2*y-y^3
            a7=axy7*z
c
            fsm1=zexpstore(indexz)
            vals=valxy*fsm1
c
            fsm=a1*vals
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now F1 is done
            fsm=a2*vals
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now F2 is done
            fsm=a3*vals
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now F3 is done
            fsm=a4*vals
            fctfsm=fc*fsm
            denssumm4=d4*fctfsm
            focksumm4=focksumm4+pot*fctfsm
c now F4 is done
            fsm=a5*vals
            fctfsm=fc*fsm
            denssumm5=d5*fctfsm
            focksumm5=focksumm5+pot*fctfsm
c now F5 is done
            fsm=a6*vals
            fctfsm=fc*fsm
            denssumm6=d6*fctfsm
            focksumm6=focksumm6+pot*fctfsm
c now F6 is done
            fsm=a7*vals
            fctfsm=fc*fsm
            denssumm7=d7*fctfsm
            focksumm7=focksumm7+pot*fctfsm
c now F7 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                          denssumm4+denssumm5+denssumm6+denssumm7
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
      else if(itype.eq.7) then
c
c  F10 function
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
            fc=ro6(iz-npwzdiff,iy)
            pot=ro4(iz,iy)
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
            a10=z2*z
c
            fsm1=zexpstore(indexz)
            vals=valxy*fsm1
c
            fsm=a1*vals
            fctfsm=fc*fsm
            denssumm1=d1*fctfsm
            focksumm1=focksumm1+pot*fctfsm
c now F1 is done
            fsm=a2*vals
            fctfsm=fc*fsm
            denssumm2=d2*fctfsm
            focksumm2=focksumm2+pot*fctfsm
c now F2 is done
            fsm=a3*vals
            fctfsm=fc*fsm
            denssumm3=d3*fctfsm
            focksumm3=focksumm3+pot*fctfsm
c now F3 is done
            fsm=a4*vals
            fctfsm=fc*fsm
            denssumm4=d4*fctfsm
            focksumm4=focksumm4+pot*fctfsm
c now F4 is done
            fsm=a5*vals
            fctfsm=fc*fsm
            denssumm5=d5*fctfsm
            focksumm5=focksumm5+pot*fctfsm
c now F5 is done
            fsm=a6*vals
            fctfsm=fc*fsm
            denssumm6=d6*fctfsm
            focksumm6=focksumm6+pot*fctfsm
c now F6 is done
            fsm=a7*vals
            fctfsm=fc*fsm
            denssumm7=d7*fctfsm
            focksumm7=focksumm7+pot*fctfsm
c now F7 is done
            fsm=a8*vals
            fctfsm=fc*fsm
            denssumm8=d8*fctfsm
            focksumm8=focksumm8+pot*fctfsm
c now F8 is done
            fsm=a9*vals
            fctfsm=fc*fsm
            denssumm9=d9*fctfsm
            focksumm9=focksumm9+pot*fctfsm
c now F9 is done
            fsm=a10*vals
            fctfsm=fc*fsm
            denssumm10=d10*fctfsm
            focksumm10=focksumm10+pot*fctfsm
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
        call nerror(1,'make_gridfunctions_dpf',
     $    ' FTC: Can only handle S, P, D and F functions',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
      Subroutine fillin_dens(dens,dens2,ncf)
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 dens(ncf,ncf)
      real*8 dens2(ncf*(ncf+1)/2)
c
      index=1
      do icf=1,ncf
        do jcf=1,icf
          dens(icf,jcf)=dens2(index)
          dens(jcf,icf)=dens2(index)
          index=index+1
        end do
      end do
c
      return
      end
c
c

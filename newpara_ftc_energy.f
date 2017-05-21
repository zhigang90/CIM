      Subroutine Para_FTC_ENERGY(
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
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
c
c  This is the main subroutine for the FTC code
c  It computes the FTC contribution to the Fock matrix
c
c  ARGUMENTS
c
c  griddens      grid density (Rydberg  is about 9.8696*griddens**2)
c  Lxmin         minimum x coordinate of box
c  Lymin         minimum y coordinate of box
c  Lzmin         minimum z coordinate of box
c  Lxo           original box length along x
c  Lyo           original box length along y
c  Lzo           original box length along z
c  Lxe           extended box length along x
c  Lye           extended box length along y
c  Lze           extended box length along z
c  PLDmax        Cutoff distance of the Coulomb operator
c  npwx          number of grid points along x in original box
c  npwy          number of grid points along y in original box
c  npwz          number of grid points along z in original box
c  npwxe         number of grid points along x in extended box
c  npwye         number of grid points along y in extended box
c  npwze         number of grid points along z in extended box
c  isharpness    pointer to shell sharpness array
c  iinit         initialization flag
c  ncspl         number of sharp contracted shells
c  igridranges   pointer to shell grid range array
c  igridranges2    ditto
c  icssharps     pointer to array indicating which shells are sharp
c  iicsdiff      pointer to array indicating which shells are diffuse
c  iicspltype    pointer to integer array for the type of every
c                 contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
c  maxovs        maximum number of shells overlapping with the sharp shells	
c  insharpovs    pointer to integer array indicating how many shells
c                 overlap with a given sharp shell
c  isharpovs        ditto
c  denA          alpha/closed-shell density matrix (lower triangle)
c  denB           ditto  beta density matrix
c  fockA         alpha/closed-shell Fock matrix (lower triangle)
c  fockB          ditto  beta Fock matrix
c  thint         integral threshold (not needed in parallel routine)
c  dist2mp       distance cutoff for multipoles
c  imultipolecut pointer to multipole distance cutoff array
c  nomultipole   multipole flag  0 - use multipoles
c                                1 - do not use multipoles
c  conver        SCF converegence flag (logical)
c
c
c -- parallel header files
c     include 'fpvm3.h'
c     include 'txmsg.h'
c
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      real*8 denA(*),denB(*),fockA(*),fockB(*)
      logical conver
cccccc
      common/ftctime/tinit,einit,tsmooth,esmooth,tftc1,eftc1,
     $               tmix,emix,tftc2,eftc2,tfock1,efock1,tmpole,empole
cccccc
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
c
c
      call elapsec(ee1)
      call secund(tt1)
c
      call getival('na  ',na)
      call getival('ndum',ndum)
      natom = na-ndum                   ! number or real atoms
      call getival('ncs ',ncs)          ! number of shells
      call getival('ncf ',ncf)          ! number of basis functions
c
c -- check if the slaves are OK
      call para_check
c
c -- send initialization flag
      If(iinit.EQ.1) Then
        iinit = 0
      Else
        call para_bcast(iinit,TxFTCInit0)
      EndIf
c
      call mmark
c
c    allocate  the array ro4(npwz,npwy,npwx)
c
      npwxyz=npwz*npwy*npwx
      npwyz=npwz*npwy
      call getmem(npwxyz,iro4)
      call zeroit(bl(iro4),npwxyz)
c
      call elapsec(ee2)
      call secund(tt2)
      tinit = tinit + tt2-tt1
      einit = einit + ee2-ee1
c
c ***************************************************************
c   BUILD UP THE SMOOTH DENSITY
c ***************************************************************
c
      call elapsec(ee1)
      call secund(tt1)
c
      nblks=5*nslv
      ntmp=10
  5   CONTINUE
      iblksize=9*npwx/(ntmp*nblks)
      If(iblksize.eq.0) Then
       ntmp = ntmp-1
       if(ntmp.eq.0) call nerror(13,'Para_RFTC_ENERGY',
     $    'Too many Slaves. Resubmit Job with LESS Slaves!',0,0)
       go to 5
      EndIf
      nleft=npwx-nblks*iblksize
c
c -- start by sending one block of data to each slave
      ixstart=1
      do islave=1,nslv
      ixend=ixstart+iblksize-1
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c -- now send blocks on round-robin basis
      do ibl=nslv+1,nblks
      ixend=ixstart+iblksize-1
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(ixstarto,1)
      call para_unpack_int(ixendo,1)
cc      write(6,*) "Master got infos back from  ",islave
cc      write(6,*) "ixstarto,ixendo",ixstarto,ixendo
cc      call f_lush(6)
      iro4p=npwyz*(ixstarto-1)
      call para_unpack_real(bl(iro4+iro4p),
     &                npwy*npwz*(ixendo-ixstarto+1))
c
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c Now send individual indices
      do ibl=1,nleft
      ixend=ixstart
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(ixstarto,1)
      call para_unpack_int(ixendo,1)
cc      write(6,*) "Master got infos back from  ",islave
cc      write(6,*) "ixstarto,ixendo",ixstarto,ixendo
cc      call f_lush(6)
      iro4p=npwyz*(ixstarto-1)
      call para_unpack_real(bl(iro4+iro4p),
     &                npwy*npwz*(ixendo-ixstarto+1))
c
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c Now do the last nslv blocks
c at same time send termination flag
c
      do ibl=1,nslv
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(ixstarto,1)
      call para_unpack_int(ixendo,1)
      iro4p=npwyz*(ixstarto-1)
      call para_unpack_real(bl(iro4+iro4p),
     &                npwy*npwz*(ixendo-ixstarto+1))
c -- slave termination
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCAssign)
cc      write(6,*) ' MASTER:  Just sent termination to slave',islave
cc      call f_lush(6)
      end do
c
      call elapsec(ee2)
      call secund(tt2)
      tsmooth = tsmooth + tt2-tt1
      esmooth = esmooth + ee2-ee1
c
c the smooth density is in ro4
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE SMOOTH DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call elapsec(ee1)
      call secund(tt1)
c
c     allocate the real arrays bfk(npwz,npwy,npwxe),w1(npwxe+15)
c
      call getmem(npwyz*npwxe,ibfk)
      call zeroit(bl(ibfk),npwyz*npwxe)
      call getmem(npwxe+15,iw1)
c
c Now create the Coulomb potential from the smooth density in ro4
c
      call Para_COULOMB_on_FTGRID(
     &    npwx,       npwy,       npwz,       npwxe,
     &    bl(ibfk),   bl(iw1),    bl(iro4))
c
c Now the Coulomb potential of the smooth density is in ro4
c Save the smooth potential if this is the last cycle
      If(conver) CALL SAVE_POT(npwxyz,'potS',bl(iro4))
c
      call retmem(2)
c
      call secund(tt2)
      call elapsec(ee2)
      tftc1 = tftc1 + tt2-tt1
      eftc1 = eftc1 + ee2-ee1
c
c ***************************************************************
c   CALCULATE MIXED (SMOOTH + SHARP) DENSITY
c ***************************************************************
c
      call elapsec(ee1)
      call secund(tt1)
c
c    allocate  the array ro5(npwz,npwy,npwx)
c
      call getmem(npwxyz,iro5)
      call zeroit(bl(iro5),npwxyz)
c
c -- broadcast smooth Coulomb operator to all slaves
      call para_bcast_real(bl(iro4),npwz*npwy*npwx,TxFTCDen)
c
c -- no need to send any more data as slaves already know
c    which shells each slave has
c    simply get back mixed density when ready
c
      call para_reduce(bl(iro5),npwz*npwy*npwx,TxFTCBTM)
c
      call elapsec(ee2)
      call secund(tt2)
      tmix = tmix + tt2-tt1
      emix = emix + ee2-ee1
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE MIXED DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      call elapsec(ee1)
      call secund(tt1)
c
c     allocate the real arrays bfk(npwz,npwy,npwxe),w1(npwxe+15)
c
      call getmem(npwyz*npwxe,ibfk)
      call zeroit(bl(ibfk),npwyz*npwxe)
      call getmem(npwxe+15,iw1)
c
c Now create the Coulomb potential from the mixed density in ro5
c
      call Para_COULOMB_on_FTGRID(
     &    npwx,       npwy,       npwz,       npwxe,    
     &    bl(ibfk),   bl(iw1),    bl(iro5))
c
      call retmem(2)
c
c Now the Coulomb potential of the mixed density is in ro5
c
c Sum smooth and mixed potentials. The result will be in ro5
c Save the potential if this is the last cycle
c
      call daxpy(npwxyz,1.0d0,bl(iro4),1,bl(iro5),1)
      If(conver) CALL SAVE_POT(npwxyz,'potM',bl(iro5))
c
      call elapsec(ee2)
      call secund(tt2)
      tftc2 = tftc2 + tt2-tt1
      eftc2 = eftc2 + ee2-ee1
c
c ***************************************************************
c   CALCULATE FOCK MATRIX FOR SMOOTH + MIXED DENSITY
c ***************************************************************
c
      call elapsec(ee1)
      call secund(tt1)
c
      nblks=5*nslv
      ntmp=10
  6   CONTINUE
      iblksize=9*npwx/(ntmp*nblks)
      If(iblksize.eq.0) Then
       ntmp = ntmp-1
       if(ntmp.eq.0) call nerror(13,'Para_RFTC_ENERGY',
     $    'Too many Slaves. Resubmit Job with LESS Slaves!',0,0)
       go to 6
      EndIf
      nleft=npwx-nblks*iblksize
c
c -- start by sending one block of data to each slave
      ixstart=1
      do islave=1,nslv
      ixend=ixstart+iblksize-1
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
      iro5p=npwyz*(ixstart-1)
      call para_send_real(bl(iro5+iro5p),
     $                   npwy*npwz*(ixend-ixstart+1),islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c -- now send blocks on round-robin basis
      do ibl=nslv+1,nblks
      ixend=ixstart+iblksize-1
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
      iro5p=npwyz*(ixstart-1)
      call para_send_real(bl(iro5+iro5p),
     $                   npwy*npwz*(ixend-ixstart+1),islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c Now send individual indices
      do ibl=1,nleft
      ixend=ixstart
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
      iro5p=npwyz*(ixstart-1)
      call para_send_real(bl(iro5+iro5p),
     $                   npwy*npwz*(ixend-ixstart+1),islave,TxFTCAssign)
cc      write(6,*) "Master will send TxFTCAssign to  ",islave
cc      write(6,*) "ixstart,ixend",ixstart,ixend
cc      call f_lush(6)
      ixstart=ixend+1
      end do
c
c Now send termination flag
c
      do ibl=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCAssign)
cc      write(6,*) ' MASTER:  Just sent termination to slave',islave
cc      call f_lush(6)
      end do
c
      call retmark            ! deallocate ro4,ro5
c
c ***************************************************************
c   MULTIPOLE CONTRIBUTIONS FOR <CD|CD> INTEGRALS
c ***************************************************************
c
      If(nomultipole.EQ.0) Then
c
        call elapsec(em1)
        call secund(tm1)
c
c -- allocate memory for atomicpoles array
        call getmem(35*natom,iatomicpoles)
        call zeroit(bl(iatomicpoles),35*natom)
c
c -- divide ncspl among the slaves
        nsize = ncspl/nslv
        nleft = ncspl - nslv*nsize
cc        write(6,*) ' MASTER:  nsize:',nsize,' nleft:',nleft
cc        call f_lush(6)
c
c -- now send each block of indices to a slave
        istart = 1
        do islave=1,nslv
        iend = istart+nsize-1
        If(nleft.gt.0) Then
          iend = iend+1
          nleft = nleft-1
        EndIf
        call f_lush(6)
        call para_initsend
        call para_pack_int(istart,1)
        call para_pack_int(iend,1)
        call para_send_pack(islave,TxFTCMult)
        istart = iend+1
        enddo
c
c -- get atomicpoles array back from slaves
cc        write(6,*) ' MASTER - About to collect atomicpoles array'
cc        call f_lush(6)
        call para_reduce(bl(iatomicpoles),35*natom,TxFTCBTM)
c
c -- now broadcast completed atomicpoles array to slaves
        call para_bcast_real(bl(iatomicpoles),35*natom,TxFTCDen)
c
c  slaves already have indices
c  simply get completed Fock matrix back from the slaves
c
        ntri = (ncf*(ncf+1))/2
        call para_reduce(fockA,ntri,TxScfFA)
c
        call retmem(1)      ! atomicpoles array
c
        call elapsec(em2)
        call secund(tm2)
        tmult = tm2-tm1
        emult = em2-em1
        tmpole = tmpole + tmult
        empole = empole + emult
c
      Else
c
c -- get FTC Fock matrix back from slaves
        ntri = (ncf*(ncf+1))/2
        call para_reduce(fockA,ntri,TxScfFA)
c
        tmult = 0.0d0
        emult = 0.0d0
c
      EndIf
cc      write(6,*) ' Final Coulomb Fock matrix is:'
cc      do kk=1,ncf
cc      kkk = (kk*(kk-1))/2
cc      write(6,*) (fockA(kkk+jj),jj=1,kk)
cc      enddo
c
c -- determine value of Fock matrix scaling factor
      dv=(Lxo/dfloat(npwx)) * (Lyo/dfloat(npwy)) * (Lzo/dfloat(npwz))
      constfftw=1.0d0/npwxe/npwye/npwze
      call setrval('cfftw',constfftw)
      fskal0 = dv*constfftw
c
c -- scale the Fock matrix
      call Make_FockL(ncf,fskal0,fockA)
c
      call secund(tt2)
      call elapsec(ee2)
      tfock1 = tfock1 + tt2-tt1 - tmult
      efock1 = efock1 + ee2-ee1 - emult
c
      return
      end
c=======================================================================
c
      Subroutine Para_COULOMB_on_FTGRID(
     &    npwx,       npwy,       npwz,       npwxe,  
     &    bfk,        w1,         ftwork)
c
c This subroutine calculates the exact Coulomb matrix elements on the
c Fourier grid using the density and the value of D which is
c the maximum distance in the Coulomb interaction.
c
c ** PARALLEL DRIVER **
c
c Input:
c     npwz,npwy,npwx            - grid dimensions of the original box
c                                 ALL MUST BE EVEN NUMBER !!!
c     npwxe                     - expanded X grid dimension
c                                 ALL MUST BE EVEN NUMBER !!!
c     bfk(npxz,npwy,npwxe)      - work array
c     w1(npwxe+15)              - work array
c     ftwork(npwz,npwy,npwx)    - array with the density
c
c Output:
c     ftwork(npwz,npwy,npwx)    - array with the Coulomb matrix elements
c                                   on the grid points multiplied by
c                                   npwxe*npwxe*npwxe
c

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 bfk(npwz,npwy,npwxe),w1(npwxe+15)
      real*8 ftwork(npwz,npwy,npwx)
c
c
c  All this routine does is divide up the expanded box X-dimension
c  over the number of slaves, send each block to a separate slave,
c  and accumulate the final Coulomb matrix
c
c
c -- parallel header files
c     include 'fpvm3.h'
c     include 'txmsg.h'
c
c -- expand the density onto the extended grid
      ixstart = (npwxe-npwx)/2 + 1
      do ix=1,npwx
      ixe=ixstart+ix-1
      do iy=1,npwy
      do iz=1,npwz
      bfk(iz,iy,ixe)=ftwork(iz,iy,ix)
      end do
      end do
      end do
c
      npwyz = npwz*npwy
      ndex = npwxe/2 + 1       ! total number of expanded X grid points
c
      nblks=5*nslv
      ntmp=10
  5   CONTINUE
      iblksize=9*ndex/(ntmp*nblks)
      If(iblksize.eq.0) Then
       ntmp = ntmp-1
       if(ntmp.eq.0) call nerror(13,'Para_RFTC_ENERGY',
     $    'Too many Slaves. Resubmit Job with LESS Slaves!',0,0)
       go to 5
      EndIf
      nleft=ndex-nblks*iblksize
cc      write(6,*) ' ndex:',ndex,' iblksize:',iblksize,' nleft:',nleft
cc      call f_lush(6)
c
c -- initialize the FFT and forward transform along X
      call vrffti(npwxe,w1)
      call vrfftf(1,npwyz,1,npwxe,bfk,npwyz,w1)
c
c -- start by sending one block of data to each slave
      ixstart = 1
cc      write(6,*) ' MASTER - Started initial blocks'
cc      call f_lush(6)
      do islave=1,nslv
      ixend = ixstart+iblksize-1
      If(ixstart.EQ.1.OR.ixend.EQ.ndex) Then
        ndata = npwyz*(2*(ixend-ixstart)+1)
      Else
        ndata = npwyz*2*(ixend-ixstart+1)
      EndIf
      lstart = MAX(1,2*ixstart-2)
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_pack_int(lstart,1)
      call para_pack_int(ndata,1)
      call para_send_pack(islave,TxFTCFFT)
c -- need to send separately to avoid mem. allocation problems on slave      
      call para_send_real(bfk(1,1,lstart),ndata,islave,TxFTCFFT)
      ixstart=ixend+1
      end do
c
c -- now send blocks on a round-robin basis
cc      write(6,*) ' MASTER - Started round robin'
cc      call f_lush(6)
      do ibl=nslv+1,nblks
      ixend=ixstart+iblksize-1
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(lstart,1)
      call para_unpack_int(ndata,1)
      call para_unpack_real(bfk(1,1,lstart),ndata)
c
      lstart = 2*ixstart-2
      ndata = npwyz*2*(ixend-ixstart+1)
      If(ixend.EQ.ndex) ndata = ndata-npwyz
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_pack_int(lstart,1)
      call para_pack_int(ndata,1)
      call para_send_pack(islave,TxFTCFFT)
c -- need to send separately to avoid mem. allocation problems on slave      
      call para_send_real(bfk(1,1,lstart),ndata,islave,TxFTCFFT)
      ixstart=ixend+1
      end do
c
c -- now send individual indices
cc      write(6,*) ' MASTER - started individual indices'
cc      call f_lush(6)
      do ibl=1,nleft
      ixend=ixstart
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(lstart,1)
      call para_unpack_int(ndata,1)
      call para_unpack_real(bfk(1,1,lstart),ndata)
c
      lstart = 2*ixstart-2
      ndata = npwyz*2
      If(ixend.EQ.ndex) ndata = npwyz
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_pack_int(lstart,1)
      call para_pack_int(ndata,1)
      call para_send_pack(islave,TxFTCFFT)
c -- need to send separately to avoid mem. allocation problems on slave      
      call para_send_real(bfk(1,1,lstart),ndata,islave,TxFTCFFT)
      ixstart=ixend+1
      end do
c
c -- now do the last nslv blocks
c -- at same time send termination flag
      do ibl=1,nslv
      call para_recv_pack(islave,TxFTCBTM)
      call para_unpack_int(lstart,1)
      call para_unpack_int(ndata,1)
      call para_unpack_real(bfk(1,1,lstart),ndata)
c -- slave termination
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCFFT)
cc      write(6,*) ' MASTER:  Just sent FFT termination to slave',islave
cc      call f_lush(6)
      end do
c
c -- back transform along X
      call vrfftb(1,npwyz,1,npwxe,bfk,npwyz,w1)
c
c -- copy back into original array
      ixstart = (npwxe-npwx)/2 + 1
      do ix=1,npwx
      ixe=ixstart+ix-1
      do iy=1,npwy
      do iz=1,npwz
      ftwork(iz,iy,ix) = bfk(iz,iy,ixe)
      end do
      end do
      end do
c
      return
      end

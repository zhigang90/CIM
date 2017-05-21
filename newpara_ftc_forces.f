      Subroutine para_FTC_RHF_FORCES(
     $    na,         ncs,        ncf,        ibas,       ictr,
     $    griddens,   Lxmin,      Lymin,      Lzmin,      Lxo,
     $    Lyo,        Lzo,        Dmax,       npwx,       npwy,
     $    npwz,       isharpgrd,  isharpness, ncspl,   igridranges,
     $  igridranges2, ncfpl,      iicsplsize, icssharps,  iicsdiff,
     $    idens,      iyregions,  npwyregdim, iicspltype, maxovs,
     $    isharpovs,  insharpovs, ii4sharps,  iplbasdat,  icorecut,
     $    iro4,       iro5,    iisharpgridrange1, iisharpgridrange2)

      use memory
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C  This is the main driver for the parallel FTC gradient
C
C  ARGUMENTS
C
C  na                   number of (real) atoms
C  ncs                  number of shells
C  ncf                  number of contracted basis functions
C  ibas                 pointer to BASDAT array
C  ictr                 pointer to INX array (Integer basis set data)
C  griddens             grid density
C  Lxmin                minimum x coordinate of box
C  Lymin                minimum y coordinate of box
C  Lzmin                minimum z coordinate of box
C  Lxo                  original box length along x
C  Lyo                  original box length along y
C  Lzo                  original box length along z
C  DMax                 Cutoff distance of the Coulomb operator
C  npwx                 number of grid points along x in original box
C  npwy                 number of grid points along y in original box
C  npwz                 number of grid points along z in original box
C  iisharpgrd
C  isharpness           pointer to shell sharpness array
C  ncspl                number of sharp contracted shells
C  igridranges          pointer to shell grid range array
C  igridranges2           ditto
C  ncfpl                number of sharp basis functions
C  iicsplsize
C  icssharps            pointer to array indicating which shells are sharp
C  iicsdiff             pointer to array indicating which shells are diffuse
C  idens                pointer to expanded density matrix
C  iyregions            number of y-regions (normally 1)
C  npwyregdim           dimension of y-regions (normally npwy)
C  iicspltype           pointer to integer array for the type of every
C                        contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
C  maxovs               maximum number of shells overlapping with the sharp shells
C  isharpovs            pointer to integer array indicating how many shells
C                        overlap with a given sharp shell
C  insharpovs             ditto
C  ii4sharps
C  iplbasdat
C  icorecut
C  iro4                 pointer to diffuse Coulomb potential
C  iro5                 pointer to diffuse+mixed Coulomb potential
C  iisharpgridrange1
C  iisharpgridrange2
C
C
c -- parallel header files
c     include 'fpvm3.h'
c     include 'txmsg.h'
c
      real*8 Lxmin,Lymin,Lzmin,Lxo,Lyo,Lzo
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' Print out stuff on entry to <para_FTC_RHF_FORCES>'
cc      write(6,*) ' ncf:',ncf,' ncs:',ncs,' ncspl:',ncspl
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' Lxmin:',lxmin,' Lymin:',lymin,' Lzmin:',lzmin
cc      write(6,*) ' DMax:',dmax,' griddens:',griddens,' maxovs:',maxovs
cc      write(6,*) ' icorecut:',icorecut,' na:',na
cc      write(6,*) ' inx array is:'
cc      call prnt_inx(ncs,inx(ictr))
cc      write(6,*) ' BASDAT array is:'
cc      call prntmat(12,13,12,bl(ibas))
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dx=1.0d0/griddens
      dv=(Lxo/dfloat(npwx)) * (Lyo/dfloat(npwy)) * (Lzo/dfloat(npwz))
      call getrval('cfftw',constfftw)
      call getival('nprint',IPRNT)
      ncsdiff = ncs-ncspl
      npwyz = npwy*npwz
c
c -- check if the slaves are OK
      call para_check
c
      call mmark
      call getmem(na,igradx)
      call ZeroIT(bl(igradx),na)
      call getmem(na,igrady)
      call ZeroIT(bl(igrady),na)
      call getmem(na,igradz)
      call ZeroIT(bl(igradz),na)
c
c ***************************************************************
c   GRADIENT FROM THE SMOOTH DENSITY
c ***************************************************************
c
      call elapsec(ee1)
      call secund(tt1)
c
c  we are going to send blocks of the potential to each slave
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
cc      write(6,*) ' MASTER: nblks:',nblks,' iblksize:',iblksize,
cc     $           ' nleft:',nleft
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
      call para_send_real(bl(iro5+iro5p),npwyz*(ixend-ixstart+1),
     $                    islave,TxFTCAssign)
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
      call para_send_real(bl(iro5+iro5p),npwyz*(ixend-ixstart+1),
     $                    islave,TxFTCAssign)
      ixstart=ixend+1
      end do
c
c -- now send individual indices
      do ibl=1,nleft
      ixend=ixstart
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(ixstart,1)
      call para_pack_int(ixend,1)
      call para_send_pack(islave,TxFTCAssign)
      iro5p=npwyz*(ixstart-1)
      call para_send_real(bl(iro5+iro5p),npwyz*(ixend-ixstart+1),
     $                    islave,TxFTCAssign)
      ixstart=ixend+1
      end do
c
c -- send termination flag
      do ibl=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCAssign)
      enddo
c
      call elapsec(ee2)
      call secund(tt2)
      tsmooth = (tt2-tt1)/60.0d0
      esmooth = (ee2-ee1)/60.0d0
c
      write(6,1000) tsmooth,esmooth
c
c ***************************************************************
c   GRADIENT FROM THE MIXED DENSITY
c ***************************************************************
c
c  This is parallelized over the number of sharp shells which is
c  already known to each slave from the initialization step
c
c  Simply broadcast the smooth potential to each slave and get back
c  the final gradient
c
c -- broadcast smooth Coulomb operator to all slaves
      call para_bcast_real(bl(iro4),npwyz*npwx,TxFTCDen)
c
c -- accumulate mixed part of gradient
      call para_reduce(bl(igradx),na,TxFTCGX)
      call para_reduce(bl(igrady),na,TxFTCGY)
      call para_reduce(bl(igradz),na,TxFTCGZ)
c
      call elapsec(ee1)
      call secund(tt1)
      tmix = (tt1-tt2)/60.0d0
      emix = (ee1-ee2)/60.0d0
c
      write(6,1100) tmix,emix
c
      If(IPRNT.GT.2) Then
        write(6,*) ' Total FTC gradient is:'
        C = -2.0d0*constfftw*dv
        do i=0,na-1
        write(6,*) C*bl(igradx+i),C*bl(igrady+i),C*bl(igradz+i)
        enddo
      EndIf
c
c Now add the FTC to the classical-integral forces
c
      if(IPRNT.GT.2)
     $   write(6,*) ' Total 2-e forces (FTC + Classical) are:'
      call getival('lforc2',lforc2)
      C = -2.0d0*constfftw*dv
      do i=0,na-1
      ii = lforc2+3*i
      bl(ii) = bl(ii) + C*bl(igradx+i)
      bl(ii+1) = bl(ii+1) + C*bl(igrady+i)
      bl(ii+2) = bl(ii+2) + C*bl(igradz+i)
      if(IPRNT.GT.2) write(6,*) bl(ii),bl(ii+1),bl(ii+2)
      enddo
c
      call retmark ! all arrays for all gradient parts are destroyed
C
      RETURN
c
 1000 Format('Master CPU time for smooth FTC gradient =',
     $        f8.2,' Elapsed = ',f8.2,' min')
 1100 Format('Master CPU time for mixed FTC gradient = ',
     $        f8.2,' Elapsed = ',f8.2,' min')
c
      END

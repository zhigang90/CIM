      Subroutine PARA_preinit_FTC(
     &       expaccl,   dist2mp,   isharpgrd, isharpness, ncspl,
     &       ncfpl,    iicsplsize, icssharps, iicsdiff,   ilistsd,
     &       iicspltype,iplbasdat, imultipolecut,nomultipole)

      use memory
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
c
c  This subroutine initializes the FTC program.
c  It is called once only at the beginning of the calculation.
c  Output on exit comprises variable values or integer addresses
c  for array pointers
c
c  ARGUMENTS
c
c  expaccl       exponent cutoff parameter: 1-10 bigger number keeps
c                values from 1-10; the larger the value, the more
c                gaussians are retained in the classical space
c                (i.e., the smaller the value of the exponent cut-off)
c  dist2mp       distance cutoff for multipoles
c  nomultipole   multipole flag  0 - multipoles
c                                1 - do not use multipoles
c
c  ON EXIT
c
c  isharpgrd     unit number for PWDFT Coulomb I/O
c  isharpness    pointer to shell sharpness array
c  ncspl         number of sharp contracted shells
c  ncfpl         number of primitive core shells
c  iicsplsize    pointer to size of core shells array
c  icssharps     pointer to array indicating which shells are sharp
c  iicsdiff      pointer to array indicating which shells are diffuse
c  ilistsd       pointer to array listing all shells as either
c                 sharp=1  or  diffuse=0  (needed for classical ints)
c  iicspltype    pointer to integer array for the type of every
c                 contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
c  iplbasdat     modified BASDAT array for plane wave space
c  imultipolecut natoms*natoms distance cutoff array for multipoles
c
      integer expaccl
      character*256 filename,filname1,scrfile
      character*3 ch3
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(5000000)
c
c
      isharpgrd = 51          ! unit no. for FTC I/O
c
      call getival('na',na)   ! GET RID OF DUMMY ATOMS - JB
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('inuc',inuc)
c
c for the multipole stuff
      If(nomultipole.EQ.0) Then
c
c  allocate integer*1 array multipolecut(na,na)
c
        call getint_1(na*na,imultipolecut)
        call setival('imultipolecut',imultipolecut)
        call get_multipolecutmx(na,dist2mp,bl(inuc),
     &                          bl(imultipolecut))
      Else
        imultipolecut = 1
      EndIf
c
c Get the necessary memory for the sharp arrays.
c Some integer arrays are also defined
c
      call getint(ncs,isharpness)
      call IZeroIT(bl(isharpness),ncs)
c
      call getmem(4*2*10,iexpcuts)
      call ZeroIT(bl(iexpcuts),4*2*10)
      call fillin_expcuts(bl(iexpcuts))
c
      call determine_sharpdim(ncs,ncspl,ncfpl,bl(ibas),bl(ictr),
     &                        bl(iexpcuts),bl(isharpness),expaccl)
c
      call retmem(1)     ! return expcuts memory
c
      call getmem(6*ncfpl,iplbasdat)
      call ZeroIT(bl(iplbasdat),6*ncfpl)
      call getint(ncspl,iicsplsize)
      call IZeroIT(bl(iicsplsize),ncspl)
      call getint(ncspl,iicspltype)
      call IZeroIT(bl(iicspltype),ncspl)
      call getint(ncspl,icssharps)
      call IZeroIT(bl(icssharps),ncspl)
      ncsdiff=ncs-ncspl
      call getint(ncsdiff,iicsdiff)
      call IZeroIT(bl(iicsdiff),ncsdiff)
c
      call getint(ncs,icssharps2)
      call IZeroIT(bl(icssharps2),ncs)
c
      call make_sharp_arrays(ncs,ncfpl,ncspl,bl(ibas),bl(ictr),
     &bl(iplbasdat),bl(iicsplsize),bl(iicspltype),bl(isharpness),
     &bl(icssharps),bl(icssharps2),ncsdiff,bl(iicsdiff))
c
      call retmem(1)     ! return cssharp2 memory
c
      call getint(ncs,ilistsd)
      call make_listsd(ncs,ncspl,bl(icssharps),bl(ilistsd))
c
c -- save pointer to listsd array in depository
      call setival('ftc0',ilistsd)
c
c ....................................................................
c -- first bcast the filename and array sizes
      call getchval('scrf',scrfile)
c
      call para_initsend
      call para_pack_string(scrfile,256)
c
      call para_pack_int(ncspl,1)
      call para_pack_int(ncsdiff,1)
      call para_pack_int(ncfpl,1)
      call para_pack_int(nomultipole,1)
      call para_bcast_pack(TxMP2File)
c
c -- now broadcast the arrays
cc      write(6,*) ' MASTER: Broadcasting TxFTCInit:',txftcinit
cc      call f_lush(6)
      call para_initsend
      call para_pack_int(bl(iicsdiff),ncsdiff)
      call para_pack_int(bl(isharpness),ncs)
      call para_pack_int(bl(iicspltype),ncspl)
      call para_pack_int(bl(iicsplsize),ncspl)
      call para_pack_int(bl(icssharps),ncspl)
      call para_pack_int(bl(ilistsd),ncs)
      If(nomultipole.EQ.0) Then
         call para_pack_real(dist2mp,1)
c -- this array is actually integer*1, but it's OK - JB
         call para_pack_byte(bl(imultipolecut),na*na)
      EndIf
      call para_bcast_pack(TxFTCInit)
cc      write(6,*) ' MASTER: Just broadcast TxFTCInit'
cc      call f_lush(6)
c .....................................................................
c
      return
      end
c =====================================================================
c
      Subroutine PARA_FTC_INITIAL(
     &       griddens,  irangeacc, irangeacc2, expaccl,   ncspl,
     &       ncfpl,     isharpgrd, isharpness, iicsplsize,icssharps,
     &       iicspltype,iplbasdat, Lxmin,      Lxmax,     Lymin,
     &       Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
     &       Lzo,       Lxe,       Lye,        Lze,       PLDmax,
     &       npwx,      npwy,      npwz,       npwxe,     npwye,
     &       npwze,     iinit,   igridranges,igridranges2,maxovs,
     &       isharpovs, insharpovs,ii4sharps,  icorecut,
     &    iisharpgridrange1,iisharpgridrange2, griddensf)
c
c  This subroutine initializes the plane wave grid for the FTC program.
c  It is called once at the beginning of the calculation and again
c  after switching to the tight integral threshold (assuming there
c  is such a switch in the SCF)
c  Output on exit comprises variable values or integer addresses
c  for array pointers
c
c  ARGUMENTS
c
c  griddens      grid density (Rydberg  is about 9.8696*griddens**2)
c  irangeacc     range accuracy parameter
c  irangeacc2    range accuracy parameter for the mixed density
c  expaccl       exponent cutoff parameter: 1-10 bigger number keeps
c                values from 1-10; the larger the value, the more
c                gaussians are retained in the classical space
c                (i.e., the smaller the value of the exponent cut-off)
c  ncspl         number of contracted core shells
c  ncfpl         number of primitive core shells
c  isharpgrd     unit number for PWDFT Coulomb I/O
c  isharpness    pointer to shell sharpness array  
c  iicsplsize    pointer to size of core shells array
c  icssharps     pointer to array indicating which shells are sharp
c  iicspltype    pointer to integer array for the type of every
c                 contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
c  iplbasdat     modified BASDAT array for plane wave space
c
c  ON EXIT
c
c  Lxmin         minimum x coordinate of box
c  Lxmax         maximum x coordinate of box
c  Lymin         minimum y coordinate of box
c  Lymax         maximum y coordinate of box
c  Lzmin         minimum z coordinate of box
c  Lzmax         maximum z coordinate of box
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
c  iinit         initialization flag
c  igridranges   pointer to shell grid range array
c  igridranges2    ditto
c  maxovs        maximum number of shells overlapping with the sharp shells	
c  insharpovs    pointer to integer array indicating how many shells
c                 overlap with a given sharp shell
c  isharpovs        ditto
c
c

      use memory
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
c -- parallel header files
c     include 'fpvm3.h'
c     include 'txmsg.h'
c
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      integer expaccl
c
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /intbl/maxsh,inx(100)
c     common /big/bl(5000000)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! P6 polinom coefs
c
c
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
c
c
c Get the necessary memory for the sharp arrays.
c Some integer arrays are also defined
c
c
      call getmem(ncs,iranges)
      call ZeroIT(bl(iranges),ncs)
      call getmem(ncs,iranges2)
      call ZeroIT(bl(iranges2),ncs)
      call getint(6*ncs,igridranges)
      call IZeroIT(bl(igridranges),6*ncs)
      call getint(6*ncs,igridranges2)
      call IZeroIT(bl(igridranges2),6*ncs)
c
      call getmem(4*8,irangeconst)
      call ZeroIT(bl(irangeconst),4*8)
      call fillin_rangeconst(bl(irangeconst))
c
      call calc_ranges(ncs,bl(ibas),bl(ictr),bl(irangeconst),
     &     bl(iranges),
     &     Lxmin,Lxmax,Lymin,Lymax,Lzmin,
     &     Lzmax,bl(igridranges),
     &     irangeacc,griddens,bl(iranges2),bl(igridranges2),
     &     irangeacc2,npwx,npwy,npwz)
c
      call retmem(1)     ! release memory for rangeconst
c
      call calc_Dmax(ncs,bl(ibas),bl(ictr),bl(iranges),PLDmax)
c
      Lxo=Lxmax-Lxmin
      Lyo=Lymax-Lymin
      Lzo=Lzmax-Lzmin
cc
cc  dimension of ii4sharp should be at least the number of
cc  contracted core basis functions times npwx
cc
      call getmem(10*ncspl*npwx,ii4sharps)
      call ZeroIT(bl(ii4sharps),10*ncspl*npwx)
cc
      call determine_expanded_dims(npwx,npwy,npwz,npwxe,npwye,
     &           npwze,PLDmax,griddens,
     &           Lxo, Lyo, Lzo, Lxe, Lye, Lze)
c
      call getint(ncspl,insharpovs)
      call IZeroIT(bl(insharpovs),ncspl)
c
      call calc_overlapdim(ncs,ncspl,bl(ibas),bl(ictr),bl(iranges2),
     &                     bl(icssharps),bl(insharpovs),maxovs)
c
      call getint(ncspl*maxovs,isharpovs)
      call IZeroIT(bl(isharpovs),ncspl*maxovs)
c
      call getint(6*ncspl,iisharpgridrange1)
      call IZeroIT(bl(iisharpgridrange1),6*ncspl)
      call getint(6*ncspl,iisharpgridrange2)
      call IZeroIT(bl(iisharpgridrange2),6*ncspl)
c
      call getint(ncspl,iifilesplit)
      call IZeroIT(bl(iifilesplit),ncspl)
      maxfilesize=240000000  ! redundant - get rid of this   JB
c
c   set the flag iteration, which is used by calc_sharpovs_new
c
      if(griddens.ge.griddensf)then
        iteration=1
      else
        iteration=0
      endif
      call calc_sharpovs_new(ncs,ncspl,bl(ibas),bl(ictr),bl(iranges2),
     &bl(icssharps),maxovs,bl(isharpovs),bl(igridranges2),
     &bl(iifilesplit),bl(iicspltype),maxfilesize,bl(iisharpgridrange1),
     &bl(iisharpgridrange2),icorecut,iteration)
c
      call retmem(1)
c
c
c ***************************************************************
c  CALCULATION OF FTC GRID DATA STARTS HERE
c ***************************************************************
c
      iinit = 1       ! initialization flag
      call para_bcast(iinit,TxFTCInit0)
c
c -- determine value of Fock matrix scaling factor
      dv=(Lxo/dfloat(npwx)) * (Lyo/dfloat(npwy)) * (Lzo/dfloat(npwz))
      constfftw=1.0d0/npwxe/npwye/npwze
      fskal0 = dv*constfftw
c
c -- send FTC grid data
c
c -- send npwx separately - it is needed to set array sizes
c
      call para_bcast(npwx,TxFTCInit1)      
c      
      call para_initsend
      call para_pack_int(icorecut,1)
      call para_pack_int(npwy,1)
      call para_pack_int(npwz,1)
      call para_pack_int(npwxe,1)
      call para_pack_int(npwye,1)
      call para_pack_int(npwze,1)
      call para_pack_int(maxovs,1)
      call para_pack_real(Lxmin,1)
      call para_pack_real(Lymin,1)
      call para_pack_real(Lzmin,1)
      call para_pack_real(Lxe,1)
      call para_pack_real(Lye,1)
      call para_pack_real(Lze,1)
      call para_pack_real(griddens,1)
      call para_pack_real(PLDmax,1)
      call para_pack_real(str,1)
      call para_pack_real(c1r,1)
      call para_pack_real(c2r,1)
      call para_pack_real(c3r,1)
      call para_pack_real(c1l,1)
      call para_pack_real(c2l,1)
      call para_pack_real(c3l,1)
      call para_pack_real(c4l,1)
      call para_pack_real(fskal0,1)
      call para_pack_int(bl(igridranges),6*ncs)
      call para_pack_int(bl(igridranges2),6*ncs)
      call para_pack_int(bl(iisharpgridrange1),6*ncspl)
      call para_pack_int(bl(iisharpgridrange2),6*ncspl)
      call para_pack_int(bl(insharpovs),ncspl)
      call para_pack_int(bl(isharpovs),maxovs*ncspl)
      call para_pack_real(bl(iplbasdat),6*ncfpl)
      call para_bcast_pack(TxFTCInit1)
c
c -- calculate the first primitive shell indices for every
c    sharp contracted shell
c
c
c -- do shells as round robin
c -- send one shell index plus its starting value to each slave
      do ish=1,ncspl
      call calc_pshdx(ncspl,bl(iicsplsize),ish,ipshdx)
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(ish,1)
      call para_pack_int(ipshdx,1)
      call para_send_pack(islave,TxFTCAssign)
      end do
c
c -- Finished. Send termination to each slave
      do ish=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCAssign)
      end do
c
      return
      end
c .........................................................................
c
      subroutine calc_pshdx(ncspl,icsplsize,ish,ipshdx)
      implicit integer(a-z)
      dimension icsplsize(ncspl)
c
      icfpl=0
      if(ish.gt.1)then
        do ishell=1,ish-1
          ncontr=icsplsize(ishell)
          icfpl=icfpl+ncontr
        end do
      endif
      ipshdx=icfpl
c
      return
      end
c =====================================================================
c
      SUBROUTINE para_FTC_INIT_FORCES(
     $       isharpgrd, isharpness, ncspl,    ncfpl,     iicsplsize,
     $       icssharps, iicsdiff,   ilistsd, iicspltype, iplbasdat,
     $       Lxmin,     Lxmax,      Lymin,    Lymax,     Lzmin,
     $       Lzmax,     Lxo,        Lyo,      Lzo,       PLDmax,
     $       npwx,      npwy,      npwz,    igridranges,igridranges2,
     $       maxovs,   isharpovs, insharpovs, ii4sharps, icorecut,
     $     iisharpgridrange1, iisharpgridrange2, griddens)

      use memory
      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine initializes the FTC for the forces.
C  It is based on the scorresponding SCF initialization routines.
C  Output on exit comproses variable values or integer addresses
C  for array pointers
C
C  ARGUMENTS
C
C  isharpgrd     unit number for PWDFT Coulomb I/O
C  isharpness    pointer to shell sharpness array
C  ncspl         number of sharp contracted shells
C  ncfpl         number of primitive core shells
C  iicsplsize    pointer to size of core shells array
C  icssharps     pointer to array indicating which shells are sharp
C  iicsdiff      pointer to array indicating which shells are diffuse
C  ilistsd       pointer to array listing all shells as either
C                 sharp=1  or  diffuse=0  (needed for classical ints)
C  iicspltype    pointer to integer array for the type of every
C                 contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
C  iplbasdat     modified BASDAT array for plane wave space
C  Lxmin         minimum x coordinate of box
C  Lxmax         maximum x coordinate of box
C  Lymin         minimum y coordinate of box
C  Lymax         maximum y coordinate of box
C  Lzmin         minimum z coordinate of box
C  Lzmax         maximum z coordinate of box
C  Lxo           original box length along x
C  Lyo           original box length along y
C  Lzo           original box length along z
C  PLDmax        Cutoff distance of the Coulomb operator
C  npwx          number of grid points along x in original box
C  npwy          number of grid points along y in original box
C  npwz          number of grid points along z in original box
C  igridranges1  pointer to shell grid range array
C  igridranges2    ditto
C  maxovs        maximum number of shells overlapping with the sharp shells
C  insharpovs    pointer to integer array indicating how many shells
C                 overlap with a given sharp shell
C  isharpovs        ditto
C  griddens      final grid density (read from depository)
C
C
c -- parallel header files
c     include 'fpvm3.h'
c     include 'txmsg.h'
c
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo
      character*256 scrfile,filname1,filename
      character*3 ch3
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l ! P6 polinom coefs
c
      isharpgrd = 51          ! unit no. for FTC I/O
c
c -- get values from depository
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
c
c -- get FTC values
      call getival('expl',iexpl)
      call getival('irange',irangeacc)
      irangeacc2 = irangeacc
      call getrval('griddf',griddens)
C
C  Get the memory for the sharp arrays
C
      call getint(ncs,isharpness)
      call IZeroIT(bl(isharpness),ncs)
c
      call getmem(4*2*10,iexpcuts)
      call ZeroIT(bl(iexpcuts),4*2*10)
      call fillin_expcuts(bl(iexpcuts))
c
      call determine_sharpdim(ncs,ncspl,ncfpl,bl(ibas),bl(ictr),
     &                        bl(iexpcuts),bl(isharpness),iexpl)
c
      call retmem(1)     ! return expcuts memory
c
      call getmem(6*ncfpl,iplbasdat)
      call ZeroIT(bl(iplbasdat),6*ncfpl)
      call getint(ncspl,iicsplsize)
      call IZeroIT(bl(iicsplsize),ncspl)
      call getint(ncspl,iicspltype)
      call IZeroIT(bl(iicspltype),ncspl)
      call getint(ncspl,icssharps)
      call IZeroIT(bl(icssharps),ncspl)
      ncsdiff=ncs-ncspl
      call getint(ncsdiff,iicsdiff)
      call IZeroIT(bl(iicsdiff),ncsdiff)
c
      call getint(ncs,icssharps2)
      call IZeroIT(bl(icssharps2),ncs)
c
      call make_sharp_arrays(ncs,ncfpl,ncspl,bl(ibas),bl(ictr),
     &bl(iplbasdat),bl(iicsplsize),bl(iicspltype),bl(isharpness),
     &bl(icssharps),bl(icssharps2),ncsdiff,bl(iicsdiff))
c
      call retmem(1)     ! return cssharp2 memory
c
      call getint(ncs,ilistsd)
      call make_listsd(ncs,ncspl,bl(icssharps),bl(ilistsd))
c
c -- save pointer to listsd array in depository
      call setival('ftc0',ilistsd)
c
      call getmem(ncs,iranges)
      call ZeroIT(bl(iranges),ncs)
      call getmem(ncs,iranges2)
      call ZeroIT(bl(iranges2),ncs)
      call getint(6*ncs,igridranges)
      call IZeroIT(bl(igridranges),6*ncs)
      call getint(6*ncs,igridranges2)
      call IZeroIT(bl(igridranges2),6*ncs)
c
      call getmem(4*8,irangeconst)
      call ZeroIT(bl(irangeconst),4*8)
      call fillin_rangeconst(bl(irangeconst))
c
      call calc_ranges(ncs,bl(ibas),bl(ictr),bl(irangeconst),
     &     bl(iranges),
     &     Lxmin,Lxmax,Lymin,Lymax,Lzmin,
     &     Lzmax,bl(igridranges),
     &     irangeacc,griddens,bl(iranges2),bl(igridranges2),
     &     irangeacc2,npwx,npwy,npwz)
c
      call retmem(1)     ! release memory for rangeconst
c
      call calc_Dmax(ncs,bl(ibas),bl(ictr),bl(iranges),PLDmax)
c
      Lxo=Lxmax-Lxmin
      Lyo=Lymax-Lymin
      Lzo=Lzmax-Lzmin
cc
      call getmem(5*ncspl*npwx,ii4sharps)
      call ZeroIT(bl(ii4sharps),5*ncspl*npwx)
      call getint(ncspl,insharpovs)
      call IZeroIT(bl(insharpovs),ncspl)
c
      call calc_overlapdim(ncs,ncspl,bl(ibas),bl(ictr),bl(iranges2),
     &                     bl(icssharps),bl(insharpovs),maxovs)
c
      call getint(ncspl*maxovs,isharpovs)
      call IZeroIT(bl(isharpovs),ncspl*maxovs)
c
      call getint(6*ncspl,iisharpgridrange1)
      call IZeroIT(bl(iisharpgridrange1),6*ncspl)
      call getint(6*ncspl,iisharpgridrange2)
      call IZeroIT(bl(iisharpgridrange2),6*ncspl)
c
      call getint(ncspl,iifilesplit)
      call IZeroIT(bl(iifilesplit),ncspl)
      maxfilesize=240000000  ! redundant - get rid of this   JB
c
      iteration = 1       ! final grid density
      icorecut = 1
      call calc_sharpovs_new(ncs,ncspl,bl(ibas),bl(ictr),bl(iranges2),
     &bl(icssharps),maxovs,bl(isharpovs),bl(igridranges2),
     &bl(iifilesplit),bl(iicspltype),maxfilesize,bl(iisharpgridrange1),
     &bl(iisharpgridrange2),icorecut,iteration)
c
      call retmem(1)
c
c
c ***************************************************************
c  CALCULATION OF FTC GRID DATA STARTS HERE
c ***************************************************************
c
c
c -- first bcast filenames and array sizes
      call getchval('scrf',scrfile)
      call para_initsend
      call para_pack_string(scrfile,256)
c
      call para_pack_int(ncs,1)
      call para_pack_int(ncspl,1)
      call para_pack_int(ncsdiff,1)
      call para_pack_int(ncfpl,1)
      call para_pack_int(npwx,1)
      call para_pack_int(npwy,1)
      call para_pack_int(npwz,1)
      call para_bcast_pack(TxMP2File)
c -- send FTC initialization data to slaves
      call para_initsend
      call para_pack_int(bl(iicsdiff),ncsdiff)
      call para_pack_int(bl(isharpness),ncs)
      call para_pack_int(bl(iicspltype),ncspl)
      call para_pack_int(bl(iicsplsize),ncspl)
      call para_pack_int(bl(icssharps),ncspl)
      call para_pack_int(bl(ilistsd),ncs)
      call para_pack_int(icorecut,1)
      call para_pack_int(maxovs,1)
      call para_pack_real(Lxmin,1)
      call para_pack_real(Lymin,1)
      call para_pack_real(Lzmin,1)
      call para_pack_real(griddens,1)
      call para_pack_real(PLDmax,1)
      call para_pack_real(str,1)
      call para_pack_real(c1r,1)
      call para_pack_real(c2r,1)
      call para_pack_real(c3r,1)
      call para_pack_real(c1l,1)
      call para_pack_real(c2l,1)
      call para_pack_real(c3l,1)
      call para_pack_real(c4l,1)
      call para_pack_int(bl(igridranges),6*ncs)
      call para_pack_int(bl(igridranges2),6*ncs)
      call para_pack_int(bl(iisharpgridrange1),6*ncspl)
      call para_pack_int(bl(iisharpgridrange2),6*ncspl)
      call para_pack_int(bl(insharpovs),ncspl)
      call para_pack_int(bl(isharpovs),maxovs*ncspl)
      call para_pack_real(bl(iplbasdat),6*ncfpl)
      call para_bcast_pack(TxFTCInit)
c
c -- calculate the first primitive shell indices for every
c    sharp contracted shell
c
c
c -- do shells as round robin
c -- send one shell index plus its starting value to each slave
      do ish=1,ncspl
      call calc_pshdx(ncspl,bl(iicsplsize),ish,ipshdx)
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(ish,1)
      call para_pack_int(ipshdx,1)
      call para_send_pack(islave,TxFTCAssign)
      end do
c
c -- Finished. Send termination to each slave
      do ish=1,nslv
      call para_recv(imsg,islave,TxDftReq)
      call para_initsend
      call para_pack_int(0,1)
      call para_pack_int(0,1)
      call para_send_pack(islave,TxFTCAssign)
      end do
c
      RETURN
      END

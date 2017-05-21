      Subroutine preinit_FTC(
     &       expaccl,   dist2mp,   isharpgrd, isharpness, ncspl,
     &       ncfpl,    iicsplsize, icssharps, iicsdiff,   ilistsd,
     &      iicspltype, iplbasdat, imultipolecut,nomultipole)

      use memory

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
c
      integer expaccl
      character*256 jobname,filename
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(5000000)
c
c
      isharpgrd = 51          ! unit no. for FTC I/O
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('inuc',inuc)
c
c for the multipole stuff
      If(nomultipole.EQ.0) Then
c
c allocate integer*1 array multipolecut(na,na)
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
      return
      end
c =====================================================================
c
      Subroutine FTC_INITIAL(
     &       griddens,  irangeacc, irangeacc2, expaccl,   ncspl,
     &       ncfpl,     isharpgrd, isharpness, iicsplsize,icssharps,
     &       iicspltype,iplbasdat, Lxmin,      Lxmax,     Lymin,
     &       Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
     &       Lzo,       Lxe,       Lye,        Lze,       PLDmax,
     &       npwx,      npwy,      npwz,       npwxe,     npwye,
     &       npwze,     iinit,  igridranges,igridranges2, maxovs,
     &       isharpovs, insharpovs,ii4sharps,  icorecut,
     &       iisharpgridrange1,iisharpgridrange2,griddensf)
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

      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      integer expaccl
      character*256 scrfile,filename
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
c  Get the necessary memory for the sharp arrays.
c  Some integer arrays are also defined
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
c
c  open files
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len1)
      filename = scrfile(1:len1)//'.sharpgrd'
      OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
      filename = scrfile(1:len1)//'.comprcoef'
      OPEN(UNIT=isharpgrd+1,FILE=filename(1:len1+10),
     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
c
      call make_sharpgrids(
     &           ncs,ncspl,bl(icssharps),bl(igridranges2),bl(iplbasdat),
     &           ncfpl,nfsh,bl(iicsplsize),bl(iicspltype),Lxmin,
     &           Lxmax,Lymin,Lymax,Lzmin,Lzmax,
     &           griddens,isharpgrd,bl(ii4sharps),
     &           bl(iisharpgridrange2),icorecut)
c
      CLOSE (UNIT=isharpgrd,STATUS='KEEP')
      CLOSE (UNIT=isharpgrd+1,STATUS='KEEP')
C
cccccccccccccccccccccccccccccc
cc      write(6,*) ' About to leave <FTC_INITIAL>  Things are:'
cc      call print_init(griddens,irangeacc,irangeacc2,expaccl, ncspl,
cc     &     ncfpl,isharpgrd,bl(isharpness),bl(iicsplsize),bl(icssharps),
cc     &     bl(iicspltype),bl(iplbasdat), Lxmin,      Lxmax,     Lymin,
cc     &         Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
cc     &         Lzo,       Lxe,       Lye,        Lze,       PLDmax,
cc     &         npwx,      npwy,      npwz,       npwxe,     npwye,
cc     &         npwze,     iinit,bl(igridranges),bl(igridranges2),maxovs,
cc     &         bl(isharpovs),bl(insharpovs),bl(ii4sharps),icorecut,ncs,
cc     &      bl(iisharpgridrange1),bl(iisharpgridrange2), griddensf)
ccccccccccccccccccccccccccccccc
      return
      end
c=======================================================================
c - TEMPORARY - THIS WHOLE THING WILL HAVE TO GO AS IT SIMPLY DUPLICATES
c - WHAT IS ALREADY DONE IN THE FTC INITIALIZATION STEP
c
      subroutine setup_IFTC_array( ilistsd )

      use memory

      implicit real*8(a-h,o-z)
C
C  This routine sets up the integer array with sharp=1 or diffuse=0
c  for each contracted shell
C
C  ARGUMENTS
C
C  ilistsd       pointer to array listing all shells as either
C                 sharp=1  or  diffuse=0  (needed for classical ints)
C                 contracted core-like shell (s=1,p=2,l=3,d=4 etc...)
C
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
C
c -- get values from depository
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
c
c -- get FTC values
c
      call getival('expl',iexpl)
C
c-----------------------------------------------------
c get memory for final output array 
c     call getmem(ncs/2+1,ilistsd)
c
c save pointer to listsd array in depository
c
c     call setival('ftc0',ilistsd)
c-----------------------------------------------------
      call mmark
c-----------------------------------------------------
C  Get the memory for the sharp arrays
C
      call getint(ncs,isharpness)
      call getmem(4*2*10,iexpcuts)
c
      call IZeroIT(bl(isharpness),ncs/2+1)
      call ZeroIT(bl(iexpcuts),4*2*10)
c
      call fillin_expcuts(bl(iexpcuts))
c
      call determine_sharpdim(ncs,ncspl,ncfpl,bl(ibas),bl(ictr),
     &                        bl(iexpcuts),bl(isharpness),iexpl)
c
      call retmem(1)     ! return expcuts memory
c
C  Get the memory 
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
      call make_listsd(ncs,ncspl,bl(icssharps),bl(ilistsd))
c-----------------------------------------------------
      call retmark
c-----------------------------------------------------
c
      end
c =====================================================================
c
      subroutine print_init(griddens,irangeacc,irangeacc2,expaccl,ncspl,
     &         ncfpl,    isharpgrd, ISHARPNESS, ICSPLSIZE, ICSSHARPS,
     &         ICSPLTYPE, PLBASDAT,  Lxmin,      Lxmax,     Lymin,
     &         Lymax,     Lzmin,     Lzmax,      Lxo,       Lyo,
     &         Lzo,       Lxe,       Lye,        Lze,       PLDmax,
     &         npwx,      npwy,      npwz,       npwxe,     npwye,
     &         npwze,     iinit,  IGRIDRANGES, IGRIDRANGES2,maxovs,
     &         ISHARPOVS, NSHARPOVS, COMPRCOEFS,  icorecut, ncs,
     &      ISHARPGRIDRANGE1,  ISHARPGRIDRANGE2, griddensf)
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo,Lxe,Lye,Lze
      integer expaccl
      INTEGER ISHARPNESS(ncs),ICSPLSIZE(ncspl),ICSSHARPS(ncspl),
     $        ICSPLTYPE(ncspl),IGRIDRANGES(6,ncs),IGRIDRANGES2(6,ncs),
     $        ISHARPOVS(maxovs,ncspl),NSHARPOVS(ncspl)
      INTEGER ISHARPGRIDRANGE1(6,ncspl),ISHARPGRIDRANGE2(6,ncspl)
      REAL*8 PLBASDAT(6,ncfpl),COMPRCOEFS(5*ncspl*npwx)
c
c
      write(6,*) ' griddens:',griddens,' irangeacc:',irangeacc,
     $           ' irangeacc2:',irangeacc2
      write(6,*) ' expaccl:',expaccl,' ncspl:',ncspl,' ncfpl:',ncfpl,
     $           ' ncs:',ncs
      write(6,*) ' Lxmin:',lxmin,' Lxmax:',lxmax
      write(6,*) ' Lymin:',lymin,' Lymax:',lymax
      write(6,*) ' Lzmin:',lzmin,' Lzmax:',lzmax
      write(6,*) ' Lxo:',lxo,' Lyo:',lyo,' Lzo:',lzo
      write(6,*) ' Lxe:',lxe,' Lye:',lye,' Lze:',lze
      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
      write(6,*) ' npwxe:',npwxe,' npwye:',npwye,' npwze:',npwze
      write(6,*) ' PLDmax:',pldmax,' iinit:',iinit,' icorecut:',icorecut
      write(6,*) ' ISHARPNESS array is:'
      write(6,*) (isharpness(i),i=1,ncs)
      write(6,*) ' ICSPLSIZE array is:'
      write(6,*) (icsplsize(i),i=1,ncspl)
      write(8,*) ' ICSSHARPS array is:'
      write(6,*) (icssharps(i),i=1,ncspl)
      write(6,*) ' ICSPLTYPE array is:'
      write(6,*) (icspltype(i),i=1,ncspl)
      write(6,*) ' NSHARPOVS array is:'
      write(6,*) (nsharpovs(i),i=1,ncspl)
      write(6,*) ' IGRIDRANGES array is:'
      do i=1,6
      write(6,*) (igridranges(i,j),j=1,ncs)
      enddo
      write(6,*) ' IGRIDRANGES2 array is:'
      do i=1,6
      write(6,*) (igridranges2(i,j),j=1,ncs)
      enddo
      write(6,*) ' ISHARPOVS array is:'
      do i=1,maxovs
      write(6,*) (isharpovs(i,j),j=1,ncspl)
      enddo
      write(6,*) ' ISHARPGRIDRANGE1 array is:'
      do i=1,6
      write(6,*) (isharpgridrange1(i,j),j=1,ncspl)
      enddo
      write(6,*) ' ISHARPGRIDRANGE2 array is:'
      do i=1,6
      write(6,*) (isharpgridrange2(i,j),j=1,ncspl)
      enddo
      write(6,*) ' PLBASDAT array is:'
      call prntmat(ncfpl,13,ncfpl,plbasdat)
cc      write(6,*) ' COMPRCOEFS array is:'
cc      call prntmat(npwx,5*nscpl,npwx)
c
      return
      end

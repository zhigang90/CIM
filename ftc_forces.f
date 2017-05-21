      SUBROUTINE FTC_INIT_FORCES(
     $       isharpgrd, isharpness, ncspl,    ncfpl,     iicsplsize,
     $       icssharps, iicsdiff,   ilistsd, iicspltype, iplbasdat,
     $       Lxmin,     Lxmax,      Lymin,    Lymax,     Lzmin,
     $       Lzmax,     Lxo,        Lyo,      Lzo,       PLDmax,
     $       npwx,      npwy,      npwz,    igridranges,igridranges2,
     $       maxovs,   isharpovs, insharpovs, ii4sharps, icorecut,
     $     iisharpgridrange1, iisharpgridrange2, griddens)

      use memory

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
      real*8 Lxmin,Lxmax,Lymin,Lymax,Lzmin,Lzmax
      real*8 Lxo,Lyo,Lzo
      character*256 scrfile,filename
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
C
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
C
C  Sharp grids file should already be available from SCF step
C  HOWEVER, need to read back comprcoefs array
C
cc      call getchval('scrf',scrfile)
cc      call rmblan(scrfile,256,len1)
cc      filename = scrfile(1:len1)//'.comprcoef'
cc      OPEN(UNIT=isharpgrd+1,FILE=filename(1:len1+10),
cc     $     FORM='UNFORMATTED',STATUS='OLD')
cc      READ(isharpgrd+1) ifcount
cc      call ReadBinary1(isharpgrd+1,ifcount,bl(ii4sharps))
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
cc      filename = scrfile(1:len1)//'.comprcoef'
cc      OPEN(UNIT=isharpgrd+1,FILE=filename(1:len1+10),
cc     &     FORM='UNFORMATTED',STATUS='UNKNOWN')
c
      call make_sharpgrids(
     &           ncs,ncspl,bl(icssharps),bl(igridranges2),bl(iplbasdat),
     &           ncfpl,nfsh,bl(iicsplsize),bl(iicspltype),Lxmin,
     &           Lxmax,Lymin,Lymax,Lzmin,Lzmax,
     &           griddens,isharpgrd,bl(ii4sharps),
     &           bl(iisharpgridrange2),icorecut)
c
      CLOSE (UNIT=isharpgrd,STATUS='KEEP')
cc      CLOSE (UNIT=isharpgrd+1,STATUS='KEEP')
C
      RETURN
      END
c
c**********************************************************************
c
      SUBROUTINE GenIFTC1

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine is needed because L-Shells cannot be handled
C  properly by the Classical gradient code and need to be
C  segmented into separate S and P shells. There is an array
C  in the FTC which indicates which shells are "core" and
C  which are "diffuse". Once L-shells have been segmented,
C  this array unfortunately needs to be regenerated otherwise
C  it is wrong.
C
C  This array is in fact needed in the gradient module ONLY by
C  the classical integral derivative code. Consequently it is
C  now generated in this routine for both non-L and L-Shell
C  basis sets.
C
c     common /big/bl(3000)
c     common /intbl/maxsh,inx(100)
c
c -- get values from depository
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('expl',iexpl)
c
c -- get memory for the listsd array
      call getint(ncs,ilistsd)
c
c -- get temporary memory
      call getint(ncs,isharpness)
      call IZeroIT(bl(isharpness),ncs)
      call getmem(4*2*10,iexpcuts)
      call ZeroIT(bl(iexpcuts),4*2*10)
c
      call fillin_expcuts(bl(iexpcuts))
      call determine_sharpdim(ncs,ncspl,ncfpl,bl(ibas),bl(ictr),
     &                        bl(iexpcuts),bl(isharpness),iexpl)
c
      call make_listsd1(ncs,bl(isharpness),bl(ilistsd))
c
c  -- remove temporary storage
      call retmem(2)
c
c -- save pointer to listsd array in depository
      call setival('ftc0',ilistsd)
c
      RETURN
      END
c
c**********************************************************************
c
      subroutine make_listsd1(ncs,sharp,listsd)
      IMPLICIT INTEGER(A-Z)
C
C  This subroutine builds up the listsd integer array.
C  It sets 0 for diffuse shells and 1 for compact (sharp).
C
C  ARGUMENTS
C
C  ncs     -  number of shells
C  sharp   -  integer array listing shell sharpness
C  listsd  -  on exit, integer array listing all shells
C
      Dimension sharp(ncs),listsd(ncs)
c
      DO 10 I=1,ncs
      If(sharp(I).EQ.3) Then
        listsd(I) = 1
      Else
        listsd(I) = 0
      EndIf
 10   CONTINUE
C
      RETURN
      END
c
c**********************************************************************
c
      SUBROUTINE FORCE2_FTC(
     $    natom,      ictr,       griddens,   Lxmin,     Lymin,
     $    Lzmin,      Lxo,        Lyo,        Lzo,        Dmax,
     $    npwx,       npwy,       npwz,     isharpgrd,    isharpness,
     $    ncspl,   igridranges, igridranges2, ncfpl,      iicsplsize,
     $    icssharps,  iicsdiff,   iicspltype, maxovs,     isharpovs,
     $    insharpovs, ii4sharps,  iplbasdat,  icorecut,
     $  iisharpgridrange1, iisharpgridrange2)

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Main driving routine for FTC forces
C
C
      real*8 Lxmin,Lymin,Lzmin,Lxo,Lyo,Lzo
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
      Character*256 jobname
c
      Parameter (IUnit=1)
      Common /job/jobname,lenJ
c
c
c -- get FTC and other data from depository
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
c
c -- recover density and expand into full matrix
c    (we are going to read the density matrix back from disk as
c     only god knows what KW has done with it)
      call getmem(ncf*ncf,idens)
      call getival('ldensi',lden)
cc      np4 = 4
cc      call rea(bl(lden),(ncf*(ncf+1))/2,np4,'den0_rhf')
      call expand(ncf,bl(lden),bl(idens))
c
c -- allocate storage for and read in diffuse + mixed potential
      npwxyz=npwz*npwy*npwx
      call getmem(npwxyz,iro4)
      call getmem(npwxyz,iro5)
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.potS',
     $      FORM='UNFORMATTED',STATUS='OLD')
      CALL ReadBinary1(IUnit,npwxyz,bl(iro4))
      CLOSE (UNIT=IUnit,STATUS='keep')
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.potM',
     $      FORM='UNFORMATTED',STATUS='OLD')
      CALL ReadBinary1(IUnit,npwxyz,bl(iro5))
      CLOSE (UNIT=IUnit,STATUS='keep')
c
cc      write(6,*) ' Just read FTC potential files'
c
c -- looks like the y-dimension is divided up into regions
c -- not sure if this is due to excessive memory demands or otherwise
c -- setting values to full range (1 pass)
      iyregions = 1
      npwyregdim = npwy
c
      CALL para_FTC_RHF_FORCES(
     $    natom,      ncs,        ncf,        ibas,       ictr,
     $    griddens,   Lxmin,      Lymin,      Lzmin,      Lxo,
     $    Lyo,        Lzo,        Dmax,       npwx,       npwy,
     $    npwz,       isharpgrd,  isharpness, ncspl,   igridranges,
     $  igridranges2, ncfpl,      iicsplsize, icssharps,  iicsdiff,
     $    idens,      iyregions,  npwyregdim, iicspltype, maxovs,
     $    isharpovs,  insharpovs, ii4sharps,  iplbasdat,  icorecut,
     $    iro4,       iro5,    iisharpgridrange1, iisharpgridrange2)
C
      call retmem(3)
c
c -- delete files once FTC part is complete
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.potS',
     $      FORM='UNFORMATTED',STATUS='OLD')
      CLOSE (UNIT=IUnit,STATUS='DELETE')
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.potM',
     $      FORM='UNFORMATTED',STATUS='OLD')
      CLOSE (UNIT=IUnit,STATUS='DELETE')
c
cc      write(6,*) ' Just deleted FTC potential files'
C 
      return
      end
c
c***********************************************************************
c
      Subroutine FTC_RHF_FORCES(
     $    na,         ncs,        ncf,        ibas,       ictr,
     $    griddens,   Lxmin,      Lymin,      Lzmin,      Lxo,
     $    Lyo,        Lzo,        Dmax,       npwx,       npwy,
     $    npwz,       isharpgrd,  isharpness, ncspl,   igridranges,
     $  igridranges2, ncfpl,      iicsplsize, icssharps,  iicsdiff,
     $    idens,      iyregions,  npwyregdim, iicspltype, maxovs,
     $    isharpovs,  insharpovs, ii4sharps,  iplbasdat,  icorecut,
     $    iro4,       iro5,    iisharpgridrange1, iisharpgridrange2)
c

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      real*8 Lxmin,Lymin,Lzmin,Lxo,Lyo,Lzo
      character*256 scrfile,filename
c
c     common /intbl/maxsh,inx(100)
c     common /big/bl(300)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' Print out stuff on entry to <FTC_RHF_FORCES>'
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
cc      constfftw=1.0d0/npwxe/npwye/npwze
      call getrval('cfftw',constfftw)
      call getival('nprint',IPRNT)
      ncsdiff=ncs-ncspl
c
      call mmark
      call getmem(na,igradx)
      call ZeroIT(bl(igradx),na)
      call getmem(na,igrady)
      call ZeroIT(bl(igrady),na)
      call getmem(na,igradz)
      call ZeroIT(bl(igradz),na)
c
c  Now the smooth part is coming and the 0-th mixed part together
c
      call elapsec(ee1)
      call secund(tt1)
c
      call mmark !arrays for the smooth part
c
      call getmem(3*npwy*npwz,iro1)
      call ZeroIT(bl(iro1),3*npwyregdim*npwz)
      call getmem(3*npwy*npwz,iro2)
      call ZeroIT(bl(iro2),3*npwyregdim*npwz)
      call getmem(npwy*npwz,iro3)
      call ZeroIT(bl(iro3),npwyregdim*npwz)
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges),iystoredim,izstoredim,ipexpdim)
c
      call getmem(iystoredim,iyexpstore)
      call ZeroIT(bl(iyexpstore),iystoredim)
      call getmem(iystoredim,iy2expstore)
      call ZeroIT(bl(iy2expstore),iystoredim)
      call getmem(iystoredim,iy3expstore)
      call ZeroIT(bl(iy3expstore),iystoredim)
c
      call getmem(izstoredim,izexpstore)
      call ZeroIT(bl(izexpstore),izstoredim)
      call getmem(izstoredim,iz2expstore)
      call ZeroIT(bl(iz2expstore),izstoredim)
      call getmem(izstoredim,iz3expstore)
      call ZeroIT(bl(iz3expstore),izstoredim)
c
      call getint(ipexpdim,iipyexpstore)
      call IZeroIT(bl(iipyexpstore),ipexpdim)
      call getint(ipexpdim,iipzexpstore)
      call IZeroIT(bl(iipzexpstore),ipexpdim)
c
      call calc_exp_store_yz(
     &    ncsdiff,        bl(iicsdiff),   ncs,           bl(ictr),
     &    bl(ibas),       bl(igridranges),iystoredim,   bl(iyexpstore),
     &    bl(iy2expstore),bl(iy3expstore),izstoredim,   bl(izexpstore),
     &    bl(iz2expstore),bl(iz3expstore),ipexpdim,    bl(iipyexpstore),
     &   bl(iipzexpstore),Lymin,          Lzmin,         griddens)
c
      ishift=0
      npwregx=1
      do ix=1,npwx
c
        xmin=Lxmin+float(ix-1)*dx
c
        do iyreg=1,iyregions
          iyregmin=(iyreg-1)*npwyregdim + 1
          iyregmax=min(iyregmin+npwyregdim-1,npwy)
          npwregy=iyregmax-iyregmin+1
          ymin=Lymin+float(iyregmin-1)*dx
c
          call calc_smooth_gradients(
     $        ncf,         bl(iro1),       npwregx,        npwregy,
     $        npwz,     bl(igridranges),   bl(iro2),       ncs,
     $     bl(isharpness), ncspl,          bl(iro3),       bl(idens),
     $        bl(ictr),   bl(ibas),       xmin,           ymin,
     $        Lzmin,     bl(iro5+ishift),  ix,             iyregmin,
     &        iyregmax,    griddens,       iystoredim, bl(iyexpstore),
     &     bl(iy2expstore),bl(iy3expstore),izstoredim, bl(izexpstore),
     &     bl(iz2expstore),bl(iz3expstore),ipexpdim,  bl(iipyexpstore),
     &     bl(iipzexpstore), na,           bl(igradx),     bl(igrady),
     &        bl(igradz))
c
          ishift=ishift+npwregy*npwz
c
        end do
c
      end do
c
      call secund(tt2)
      call elapsec(ee2)
      tsmooth = (tt2-tt1)/60.0d0
      esmooth = (ee2-ee1)/60.0d0
c
      write(6,1000) tsmooth,esmooth
c
      If(IPRNT.GT.3) Then
        write(6,*) ' First part of FTC gradient is:'
        C = -2.0d0*constfftw*dv
        do i=0,na-1
        write(6,*) C*bl(igradx+i),C*bl(igrady+i),C*bl(igradz+i)
        enddo
      EndIf
c
      call retmark !arrays for the smooth part are destroyed
c
c  Smooth part of the gradients end.
c  Now one of the mixed part with analytical FT
c  for the sharp function derivatives is coming.
c
      call mmark ! arrays for the first mixed part
c
      call getmem(npwy*npwz,iro3)
      call ZeroIT(bl(iro3),npwy*npwz)
c
c Now figure out some dimensons to pre-calculate the exp functions
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges2),iystoredim,izstoredim,ipexpdim)
c
      call getmem(iystoredim,iyexpstore)
      call ZeroIT(bl(iyexpstore),iystoredim)
      call getmem(iystoredim,iy2expstore)
      call ZeroIT(bl(iy2expstore),iystoredim)
      call getmem(iystoredim,iy3expstore)
      call ZeroIT(bl(iy3expstore),iystoredim)
c
      call getmem(izstoredim,izexpstore)
      call ZeroIT(bl(izexpstore),izstoredim)
      call getmem(izstoredim,iz2expstore)
      call ZeroIT(bl(iz2expstore),izstoredim)
      call getmem(izstoredim,iz3expstore)
      call ZeroIT(bl(iz3expstore),izstoredim)
c
      call getint(ipexpdim,iipyexpstore)
      call IZeroIT(bl(iipyexpstore),ipexpdim)
      call getint(ipexpdim,iipzexpstore)
      call IZeroIT(bl(iipzexpstore),ipexpdim)
c
c Now precalculating the 1D exp functions along y and z for all
c smooth primitive shells
c
      call calc_exp_store_yz(
     &    ncsdiff,       bl(iicsdiff),    ncs,         bl(ictr),
     &    bl(ibas),      bl(igridranges2),iystoredim,  bl(iyexpstore),
     &    bl(iy2expstore),bl(iy3expstore),izstoredim,  bl(izexpstore),
     &    bl(iz2expstore),bl(iz3expstore),ipexpdim,    bl(iipyexpstore),
     &    bl(iipzexpstore),Lymin,         Lzmin,       griddens)
c
      npwregx=1
      xmin0=Lxmin-dx
c
c
      icfpl=0
      do icspl=1,ncspl ! loop over the sharp contracted shells
        call mmark
c
        call AFT_sharpgrid_derivatives_arrange(
     &    ncs,           ncspl,         bl(icssharps), bl(igridranges2),
     &    bl(iicsplsize),bl(iicspltype),Lxmin,         Lymin,
     &    Lzmin,         griddens,      icspl,         ncontr,
     &    rLx,           rLy,           rLz,           xp,
     &    yp,            zp,            npwxs,         npwys,
     &    npwzs,         ics,           itype,         nzeffdim,
     &    iexpfuncx,     iexpfuncy,     iexpfuncz,     icftwork,
     &    icftwork2,     iftwork,       isxmin,        isxmax,
     &    isymin,        isymax,        iszmin,        iszmax,
     &    bl(ictr),     ifirstf,       ilastf,  bl(iisharpgridrange2),
     &    icorecut)
c
        ymin=Lymin+float(isymin-1)*dx
        npwregy=npwys
        ishift1help=3*npwys*npwzs
        ishift2help=npwy*npwz
        valgradx=0.0d0
        valgrady=0.0d0
        valgradz=0.0d0
        ncomp=ilastf-ifirstf+1
        do ifunc=ifirstf,ilastf !loop over the sharp basis functions
                                !inside the shell
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
c
c Now the derivatives (3D) for the given ifunc are in the ftwork array.
c
          if(icomp .ne. ncomp) icfpl=icfpl-ncontr
c
c Now calling the main routine to calculate the
c sum over jfunc desmx(ifunc,jfunc)*bfc(smoot;jfunc) and multiply
c this with the three derivatives of the sharp ifunc and with the
c smooth Coulomb potential.
c
          ishift1=0
          ishift2=(isxmin-1)*npwy*npwz+(isymin-1)*npwz
          xmin=xmin0+(isxmin-1)*dx
c
          do ix=isxmin,isxmax
            xmin=xmin+dx
c
            call calc_mix1_gradients(
     &    ncf,      bl(iftwork+ishift1),npwregx,       npwregy,
     &    npwz,        bl(igridranges2),ncs,           bl(isharpness),
     &    ncspl,         bl(iro3),      bl(idens),     bl(ictr),
     &    bl(ibas),      xmin,          ymin,          Lzmin,
     &    ix,            isymin,        isymax,        icspl,
     &    ifunc,         bl(isharpovs), maxovs,        bl(insharpovs),
     &    npwzs,      bl(iro4+ishift2), griddens,      iystoredim,
     &  bl(iyexpstore),bl(iy2expstore),bl(iy3expstore),izstoredim,
     &  bl(izexpstore),bl(iz2expstore),bl(iz3expstore),ipexpdim,
     &  bl(iipyexpstore),bl(iipzexpstore),iszmin,      iszmax,
     &    valgradx,     valgrady,        valgradz,     icorecut,
     &  bl(iisharpgridrange1), bl(iisharpgridrange2))
c
            ishift1=ishift1+ishift1help
            ishift2=ishift2+ishift2help
c
          end do !over ix
c
        end do !over ifunc
c
        call add_to_gradients(ncs,bl(ictr),ics,na,valgradx,
     &           valgrady,valgradz,bl(igradx),bl(igrady),bl(igradz))
c
        call retmark
      end do ! end loop over the sharp contracted shells
c
      call retmark !tmp arrays for the derivatives of the first
                   !mixed part are destroyed
c
      If(IPRNT.GT.3) Then
        write(6,*) ' Second part of FTC gradient is:'
        C = -2.0d0*constfftw*dv
        do i=0,na-1
        write(6,*) C*bl(igradx+i),C*bl(igrady+i),C*bl(igradz+i)
        enddo
      EndIf
c
c Now the second mixed part with AFT of the cores and derivatives
c of the smooths in the inner cycle is coming.
c
      call mmark ! Arrays for the mix2 gradient part
c
      call getmem(3*npwy*npwz,iro5v)
      call ZeroIT(bl(iro5v),3*npwy*npwz)
c
      call calc_exp_dims(ncsdiff,bl(iicsdiff),ncs,bl(ictr),bl(ibas),
     &           bl(igridranges2),iystoredim,izstoredim,ipexpdim)
c
      call getmem(iystoredim,iyexpstore)
      call ZeroIT(bl(iyexpstore),iystoredim)
      call getmem(iystoredim,iy2expstore)
      call ZeroIT(bl(iy2expstore),iystoredim)
      call getmem(iystoredim,iy3expstore)
      call ZeroIT(bl(iy3expstore),iystoredim)
c
      call getmem(izstoredim,izexpstore)
      call ZeroIT(bl(izexpstore),izstoredim)
      call getmem(izstoredim,iz2expstore)
      call ZeroIT(bl(iz2expstore),izstoredim)
      call getmem(izstoredim,iz3expstore)
      call ZeroIT(bl(iz3expstore),izstoredim)
c
      call getint(ipexpdim,iipyexpstore)
      call IZeroIT(bl(iipyexpstore),ipexpdim)
      call getint(ipexpdim,iipzexpstore)
      call IZeroIT(bl(iipzexpstore),ipexpdim)
c
      call calc_exp_store_yz(
     &    ncsdiff,       bl(iicsdiff),  ncs,           bl(ictr),
     &    bl(ibas),    bl(igridranges2),iystoredim,    bl(iyexpstore),
     &  bl(iy2expstore),bl(iy3expstore),izstoredim,    bl(izexpstore),
     &  bl(iz2expstore),bl(iz3expstore),ipexpdim,      bl(iipyexpstore),
     &  bl(iipzexpstore),Lymin,         Lzmin,         griddens)
c
c open file
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len1)
      filename = scrfile(1:len1)//'.sharpgrd'
      OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &     FORM='UNFORMATTED',STATUS='OLD')
c
      ifunc=0
      npwregx=1
      xmin0=Lxmin-dx
      indxcomprcoef=0
      do icspl=1,ncspl ! Loop over the sharp shells
c
        call mixed_density_help1(
     &    icspl,         bl(iicspltype),ncspl,         bl(icssharps),
     &    ncomp,       bl(igridranges2),ncs,           isxmin,
     &    isxmax,        isymin,        isymax,        iszmin,
     &    iszmax,        bl(ictr),     ifirstf,       ilastf,
     &    icorecut,bl(iisharpgridrange2))
c
        npwsx=isxmax-isxmin+1
        npwsy=isymax-isymin+1
        npwsz=iszmax-iszmin+1
        ymin=Lymin+float(isymin-1)*dx
        npwregy=isymax -isymin+1
        ishifthelp1=(isymin-1)*npwz
        ishifthelp2=npwregy*npwz+(npwy-isymax)*npwz
        ishift0=(isxmin-1)*npwy*npwz
        ishift2help1=npwregy*npwsz
c
        call mmark
        call getint(npwsy*npwsz,iro62d)
c
        do ifunc=ifirstf,ilastf ! Loop for the sharp contracted basis
                                ! functions which are on the disk
          ishift2=0
          ishift=ishift0
          xmin=xmin0+(isxmin-1)*dx
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
     &    npwz,        bl(igridranges2),ncs,           bl(isharpness),
     &    ncspl,         bl(idens),     bl(ictr),     bl(ibas),
     &    xmin,          ymin,          Lzmin,         bl(iro5v),
     &    ix,            isymin,        isymax,        icspl,
     &    ifunc,         bl(isharpovs), maxovs,        bl(insharpovs),
     &    npwsz,         iszmin,      bl(iro4+ishift), griddens,
     &    iystoredim,  bl(iyexpstore),bl(iy2expstore),bl(iy3expstore),
     &    izstoredim,  bl(izexpstore),bl(iz2expstore),bl(iz3expstore),
     &    ipexpdim,  bl(iipyexpstore),bl(iipzexpstore),comprcoef,
     &    na,            bl(igradx),    bl(igrady),     bl(igradz),
     &    icorecut, bl(iisharpgridrange1),bl(iisharpgridrange2))
c
            ishift=ishift+ishifthelp2
            ishift2=ishift2+ishift2help1
c
          end do !end of the loop over ix
c
        end do ! end of the loop over the component
        call retmark ! ro62d is destroyed
      end do ! end of the loop over the sharp shells
c
      CLOSE(UNIT=isharpgrd)
c
      call retmark !arrays for the mix2 part of the gradients
                   !are destroyed
c
      call elapsec(ee1)
      call secund(tt1)
      tmix = (tt1-tt2)/60.0d0
      emix = (ee1-ee2)/60.0d0
c
      write(6,1100) tmix,emix
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
 1000 Format('Master CPU time for smooth FTC gradient =',
     $        f8.2,' Elapsed = ',f8.2,' min')
 1100 Format('Master CPU time for mixed FTC gradient = ',
     $        f8.2,' Elapsed = ',f8.2,' min')
      END
c
c***********************************************************************
c
      Subroutine calc_smooth_gradients(
     &    ncf,           ro1,           npwx,          npwy,
     &    npwz,          gridranges,    ro2,           ncs,
     &    sharpness,     ncspl,         ro3,           dens,
     &    inx,           basdat,        Lxmin,         Lymin,
     &    Lzmin,         ro4,           ixreg,         iyregmin,
     &   iyregmax,       griddens,      iystoredim,    yexpstore,
     &   y2expstore,    y3expstore,     izstoredim,    zexpstore,
     &   z2expstore,    z3expstore,     ipexpdim,      ipyexpstore,
     &   ipzexpstore,    natom,         gradx,         grady,
     &    gradz)
c
c
      IMPLICIT REAL*8(A-H,O-Z)
c
c So npwx=1 in this subroutine and npwy=npwregy in this subroutine
c Actual values of the Lxmin and Lymin are coming from outside
c
      real*8 ro1(npwz,npwy,3) ! help variable for the derivatives
      real*8 ro2(npwz,npwy,3) ! derivatives; 1 -> x; 2 -> y; z -> 3
      real*8 ro3(npwz,npwy)   ! will hold the gridfunctions
      real*8 ro4(npwz,npwy)   ! hold the smooth Coulomb potential
      integer gridranges(6,ncs)
      integer sharpness(ncs)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 dens(ncf,ncf)
      real*8 denspart(10)     ! density matrix contribution
      real*8 Lxmin,Lymin,Lzmin
c
      real*8 yexpstore(iystoredim)
      real*8 y2expstore(iystoredim)
      real*8 y3expstore(iystoredim)
      real*8 zexpstore(izstoredim)
      real*8 z2expstore(izstoredim)
      real*8 z3expstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      real*8 gradx(natom)
      real*8 grady(natom)
      real*8 gradz(natom)
c
c
      istep=1
      ip=0
c
      do ics=1,ncs  ! loop over the contracted shell ics
c
        isharpness=sharpness(ics)
        if(isharpness .eq. 3) goto 150 !the function is sharp
        ixmin=gridranges(1,ics)
        ixmax=gridranges(2,ics)
c
        if((ixmin .gt. ixreg) .or. (ixmax .lt. ixreg)) then
          goto 150                      !the function is zero
        else
          ixminold=ixmin
          ixmin=ixreg
          ixmax=ixreg
        end if
c
        ifirstf=inx(11,ics)+1
        ilastf=inx(10,ics)
        iatom=inx(2,ics)
c
        iymin=gridranges(3,ics)
        iymax=gridranges(4,ics)
        iyminold=iymin
        iymaxold=iymax
c
        call get_y_gridrange(iymin,iymax,ixminold,ixreg,icsrad,
     &                     icsradx2,iycenter)
c
        if((iymin .gt. iyregmax) .or. (iymax .lt. iyregmin)) then
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
        ibeg=inx(1,ics)+1
        Px=basdat(11,ibeg)
        Py=basdat(12,ibeg)
        Pz=basdat(13,ibeg)
        iend=inx(5,ics)
        icssize=iend-ibeg+1
        itype=inx(12,ics)
c
        valgradx=0.0d0
        valgrady=0.0d0
        valgradz=0.0d0
c
        do ifunc=ifirstf,ilastf !loop over contracted basis function i
          icomp=ifunc-ifirstf+1
          cji=1.0d0
c
          if(ibeg .ne. iend) then
            do iy=iyindexmin,iyindexmax
              do iz=izmin,izmax
                ro2(iz,iy,1)=0.0d0
                ro2(iz,iy,2)=0.0d0
                ro2(iz,iy,3)=0.0d0
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
            call calc_expvaluesv(expvalxv,Px,Lxmin,eta,
     &                       cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
              ipyexp=ipyexpstore(ish)+(iymin-iyminold)
              ipzexp=ipzexpstore(ish)
              npwyexp=npwy-(iymin-iyminold)
            if(icssize .eq. 1) then
cc
              call make_gridfunction_derivatives(
     &    itype,         expvalxv,      icomp,         Px,
     &    Py,            Pz,            ro2,           iout,
     &    1,             1,             iyindexmin,    iyindexmax,
     &    izmin,         izmax,         istep,         Lxmin,
     &    Lymin,         Lzmin,         npwx,          npwy,
     &    npwz,          griddens,      xexpvalxv,     xxexpvalxv,
     & yexpstore(ipyexp),y2expstore(ipyexp),y3expstore(ipyexp),npwyexp,
     & zexpstore(ipzexp),z2expstore(ipzexp),z3expstore(ipzexp),eta)
cc
            else
cc
               call make_gridfunction_derivatives(
     &    itype,         expvalxv,      icomp,         Px,
     &    Py,            Pz,            ro1,           iout,
     &    1,             1,             iyindexmin,    iyindexmax,
     &    izmin,         izmax,         istep,         Lxmin,
     &    Lymin,         Lzmin,         npwx,          npwy,
     &    npwz,          griddens,      xexpvalxv,     xxexpvalxv,
     & yexpstore(ipyexp),y2expstore(ipyexp),y3expstore(ipyexp),npwyexp,
     & zexpstore(ipzexp),z2expstore(ipzexp),z3expstore(ipzexp),eta)
cc
            end if
c
            if(icssize .gt. 1) then
              do iy=iyindexmin,iyindexmax
                do iz=izmin,izmax
                  ro2(iz,iy,1)=ro2(iz,iy,1)+ro1(iz,iy,1)
                  ro2(iz,iy,2)=ro2(iz,iy,2)+ro1(iz,iy,2)
                  ro2(iz,iy,3)=ro2(iz,iy,3)+ro1(iz,iy,3)
                end do
              end do
            end if
c
          end do  ! end of the loop over the contractions ish
c
c zero out the ro3 working array
c
          do iy=iyindexmin,iyindexmax
            do iz=izmin,izmax
              ro3(iz,iy)=0.0d0
            end do
          end do
c
c fml
          do jcs=1,ncs   ! loop over the contracted shell jcs
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
            jfirstf=inx(11,jcs)+1
            jlastf=inx(10,jcs)
c
            if((ijxmin .eq. ijxmax) .and. (ijymin .le. ijymax) .and.
     &         (ijzmin .lt. ijzmax)) then
c
c Now ijxmin=ijxmax=ixreg
c
              jbeg=inx(1,jcs)+1
              Pxj=basdat(11,jbeg)
              Pyj=basdat(12,jbeg)
              Pzj=basdat(13,jbeg)
c
              jend=inx(5,jcs)
              jcssize=jend-jbeg+1
              jtype=inx(12,jcs)
c
              cji=1.0d0
              index=1
              do jfunc=jfirstf,jlastf
                denspart(index)=dens(jfunc,ifunc)
                index=index+1
              end do
c
c Here loop for the contractions
c
              do jsh=jbeg,jend  ! loop over the contractions
c
                eta=basdat(1,jsh)
                if(jtype .eq. 3) then
                  concoef=basdat(3,jsh)
                  ccs=basdat(2,jsh)/concoef
                else
                  concoef=basdat(2,jsh)
                end if
c
                call calc_expvaluesv(
     &    expvalxv,      Pxj,           Lxmin,         eta,
     &    cji,           concoef,       griddens,      xexpvalxv,
     &    xxexpvalxv)
c
                ipyexp=ipyexpstore(jsh)+(ijymin-jyminold)
                ipzexp=ipzexpstore(jsh)+(ijzmin-jzmin)
                npwyexp=npwy-(ijymin-jyminold)
                npwzexp=npwz-(ijzmin-jzmin)
c
                call make_gridfunctions_for_derivatives(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Pxj,           Pyj,           Pzj,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,     yexpstore(ipyexp),
     &y2expstore(ipyexp),npwyexp,    zexpstore(ipzexp),ro3,
     &    denspart,      npwzexp)
c
              end do ! end of the loop for the contractions
cc
            end if
c
 120        continue
          end do  ! end of the loop over the contracted shell jcs
c
          do iy=iyindexmin,iyindexmax
            do iz=izmin,izmax
c
              rojmu=ro3(iz,iy)*ro4(iz,iy)
              valgradx=valgradx+rojmu*ro2(iz,iy,1)
              valgrady=valgrady+rojmu*ro2(iz,iy,2)
              valgradz=valgradz+rojmu*ro2(iz,iy,3)
c
            end do
          end do
c
        end do  ! end of the loop over the contracted basis function i
c
        gradx(iatom)=gradx(iatom)+valgradx
        grady(iatom)=grady(iatom)+valgrady
        gradz(iatom)=gradz(iatom)+valgradz
c
 150    continue
      end do  ! end of the loop over the contracted shell ics
c
      return
      end
c
c***********************************************************************
c
      Subroutine AFT_sharpgrid_derivatives_arrange(
     &    ncs,           ncspl,         cssharps,      gridranges,
     &    icsplsize,     icspltype,     Lxmin,         Lymin,
     &    Lzmin,         griddens,      icspl,         ncontr,
     &    Lx,            Ly,            Lz,            xp,
     &    yp,            zp,            npwx,          npwy,
     &    npwz,          ics,           itype,         nzeffdim,
     &    iexpfuncx,     iexpfuncy,     iexpfuncz,     icftwork,
     &    icftwork2,     iftwork,       isxmin,        isxmax,
     &    isymin,        isymax,        iszmin,        iszmax,
     &    inx,           ifirstf,       ilastf,      isharpgridrange2,
     &    icorecut)
c
c This subroutine makes the necessary arrangements for
c make_sharpgrid_derivatives_work subroutine.
c
c Input:
c
c

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      dimension inx(12,ncs)
      integer cssharps(ncspl)
      integer gridranges(6,ncs)
      real*8 Lxmin,Lymin,Lzmin,Lx,Ly,Lz
      integer icsplsize(ncspl)
      integer icspltype(ncspl)
      integer isharpgridrange2(6,ncspl)
c     common /big/bl(1000000)
c
c
      dx=1.0d0/griddens
      cji=1.0d0
c
      ics=cssharps(icspl)
      itype=icspltype(icspl)
      ncontr=icsplsize(icspl)
      ifirstf=inx(11,ics)+1
      ilastf=inx(10,ics)
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
        isxmin=max(ixdmin,ixcmin)
        isymin=max(iydmin,iycmin)
        iszmin=max(izdmin,izcmin)
        isxmax=min(ixdmax,ixcmax)
        isymax=min(iydmax,iycmax)
        iszmax=min(izdmax,izcmax)
c
        npwx=isxmax-isxmin+1
        npwy=isymax-isymin+1
        npwz=iszmax-iszmin+1
        Lx=float(npwx)/griddens
        Ly=float(npwy)/griddens
        Lz=float(npwz)/griddens
c
c determine the coordinates of the origin of the actual box
c
        xp=float(isxmin-1)*dx + Lxmin + Lx/2.0d0
        yp=float(isymin-1)*dx + Lymin + Ly/2.0d0
        zp=float(iszmin-1)*dx + Lzmin + Lz/2.0d0
c
      else
        isxmin=gridranges(1,ics)
        isxmax=gridranges(2,ics)
        isymin=gridranges(3,ics)
        isymax=gridranges(4,ics)
        iszmin=gridranges(5,ics)
        iszmax=gridranges(6,ics)
c
        npwx=isxmax-isxmin+1
        npwy=isymax-isymin+1
        npwz=iszmax-iszmin+1
        Lx=float(npwx)/griddens
        Ly=float(npwy)/griddens
        Lz=float(npwz)/griddens
c
c determine the coordinates of the origin of the actual box
c
        xp=float(gridranges(1,ics)-1)*dx + Lxmin + Lx/2.0d0
        yp=float(gridranges(3,ics)-1)*dx + Lymin + Ly/2.0d0
        zp=float(gridranges(5,ics)-1)*dx + Lzmin + Lz/2.0d0
c
      end if
c
      if(itype .eq. 1) then
        ncomp=1
      else if(itype .eq. 2) then
        ncomp=3
      else if(itype .eq. 3) then
        ncomp=4
      else if(itype .eq. 4) then
        ncomp=5
      else if(itype .eq. 5) then
        ncomp=6
      else if(itype .eq. 6) then
        ncomp=7
      else if(itype .eq. 7) then
        ncomp=10
      else
        call nerror('AFT_sharpgrid_derivatives_arrange',
     $              'Error in itype',0,0)
      end if
c
      call getmem(npwx,iexpfuncx)
      call getmem(npwy,iexpfuncy)
      call getmem(npwz,iexpfuncz)
c -- Now precalculate the exp functions to speed up the analytical FT
      call Precalc_expfuncs(
     &     bl(iexpfuncx),bl(iexpfuncy),bl(iexpfuncz),npwx,npwy,
     &     npwz,Lx,Ly,Lz)
c
      nzeffdim=npwz/2+1 !for complex to real back Fourier trafo
      call getmem(3*2*npwx*npwy*nzeffdim,icftwork) !this is for anal FT
                                                   ! and back FT
      call getmem(2*npwx*npwy*nzeffdim,icftwork2) ! this is a helping
                                                  ! array
      call getmem(3*npwx*npwy*npwz,iftwork) !this will hold the results
c
      return
      end
c
c***********************************************************************
c
      Subroutine calc_mix1_gradients(
     &    ncf,           ro6,           npwx,          npwy,
     &    npwz,          gridranges,    ncs,           sharpness,
     &    ncspl,         ro3,           dens,          inx,
     &    basdat,        Lxmin,         Lymin,         Lzmin,
     &    ixreg,         iyregmin,      iyregmax,      icspl,
     &    ifunc,         sharpovs,      maxovs,        nsharpovs,
     &    npwsz,         ro4,           griddens,      iystoredim,
     &    yexpstore,     y2expstore,    y3expstore,    izstoredim,
     &    zexpstore,     z2expstore,    z3expstore,    ipexpdim,
     &    ipyexpstore,   ipzexpstore,   iszmin,        iszmax,
     &    valgradx,      valgrady,      valgradz,      icorecut,
     &  isharpgridrange1,isharpgridrange2)
c

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c
c So npwx=1 in this subroutine and npwy=npwregy in this subroutine
c Actual values of the Lxmin and Lymin are coming from outside
c
      real*8 ro6(3,npwsz,npwy) ! sharp function derivatives
      real*8 ro3(npwz,npwy)    ! smooth function
      real*8 ro4(npwz,npwy)    ! Coulomb potencial from the smooth density
      real*8 denspart(10)      ! density mx contribution
      integer gridranges(6,ncs)
      integer sharpness(ncs)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 dens(ncf,ncf)
      real*8 Lxmin,Lymin,Lzmin
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
c
      real*8 yexpstore(iystoredim)
      real*8 y2expstore(iystoredim)
      real*8 y3expstore(iystoredim)
      real*8 zexpstore(izstoredim)
      real*8 z2expstore(izstoredim)
      real*8 z3expstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      integer isharpgridrange1(6,ncspl)
      integer isharpgridrange2(6,ncspl)
c     common /big/bl(5000000)
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' Print out stuff on entry to <calc_mix1_gradients>'
cc      write(6,*) ' ncf:',ncf,' ncs:',ncs,' ncspl:',ncspl
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' Lxmin:',lxmin,' Lymin:',lymin,' Lzmin:',lzmin
cc      write(6,*) ' ixreg:',ixreg,' iyregmin:',iyregmin,' iyregmax:',
cc     $             iyregmax
cc      write(6,*) ' icspl:',icspl,' ifunc:',ifunc,' maxovs:',maxovs
cc      write(6,*) ' npwsz:',npwsz,' griddens:',griddens
cc      write(6,*) ' iystoredim:',iystoredim,' izstoredim:',izstoredim,
cc     $           ' ipexpdim:',ipexpdim
cc      write(6,*) ' iszmin:',iszmin,' iszmax:',iszmax,
cc     $           ' icorecut:',icorecut
cc      write(6,*) ' ro6 array is:'
cc      call prntmat(npwsz*npwy,6,npwsz*npwy,ro6)
cc      write(6,*) ' gridranges array is:'
cc      do i=1,6
cc      write(6,*) (gridranges(i,j),j=1,ncs)
cc      enddo
cc      write(6,*) ' sharpness array is:'
cc      write(6,*) (sharpness(i),i=1,ncs)
cc      write(6,*) ' ro3 arrsy is:'
cc      do i=1,npwz
cc      write(6,*) (ro3(i,j),j=1,npwy)
cc      enddo
cc      write(6,*) ' density matrix is:'
cc      call prntmat(ncf,ncf,ncf,dens)
cc      write(6,*) ' inx array is:'
cc      do i=1,12
cc      write(6,*) (inx(i,j),j=1,ncs)
cc      enddo
cc      write(6,*) ' BASDAT array is:'
cc      call prntmat(16,13,16,basdat)
cc      write(6,*) ' sharpovs array is:'
cc      call prntmat(ncspl,maxovs,ncspl,sharpovs)
cc      write(6,*) ' nsharpovs array is:'
cc      write(6,*) (nsharpovs(i),i=1,ncspl)
cc      write(6,*) ' ro4 array is:'
cc      do i=1,npwz
cc      write(6,*) (ro4(i,j),j=1,npwy)
cc      enddo
cc      write(6,*) ' yexpstore, y2expstore and y3expstore arrays:'
cc      do i=1,iystoredim
cc      write(6,*) i,'  ',yexpstore(i),y2expstore(i),y3expstore(i)
cc      enddo
cc      write(6,*) ' zexpstore, z2expstore and z3expstore arrays:'
cc      do i=1,izstoredim
cc      write(6,*) i,'  ',zexpstore(i),z2expstore(i),z3expstore(i)
cc      enddo
cc      write(6,*) ' ipyexpstore and ipzexpstore arrays:'
cc      do i=1,ipexpdim
cc      write(6,*) i,'  ',ipyexpstore(i),ipzexpstore(i)
cc      enddo
cc      write(6,*) ' isharpgridrange1 array is:'
cc      do i=1,6
cc      write(6,*) (isharpgridrange1(i,j),j=1,ncspl)
cc      enddo
cc      write(6,*) ' isharpgridrange2 array is:'
cc      do i=1,6
cc      write(6,*) (isharpgridrange2(i,j),j=1,ncspl)
cc      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      if(icorecut .eq. 1) then
        icorexmin2=isharpgridrange2(1,icspl)
        icorexmax2=isharpgridrange2(2,icspl)
        if(icorexmin2 .gt. ixreg) return !the given x grid does
                                         !not contribute to the core
        if(icorexmax2 .lt. ixreg) return !the given x grid does
                                         !not contribute to the core
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
      npwzdiff=iszmin-1
      cji=1.0d0
c
c zero out the ro3 working array
c
      do iy=1,npwy
        do iz=iszmin,iszmax
          ro3(iz,iy)=0.0d0
        end do
      end do
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
c  new part to calculate spherical reagions
c
        call get_y_gridrange(jymin,jymax,jxmin,ixreg,jcsrad,
     &                     jcsradx2,jycenter)
c
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
        if((ijxmin .eq. ijxmax) .and. (ijymin .le. ijymax) .and.
     &     (ijzmin .le. ijzmax)) then
          ijyindexmin=1+ijymin-iyregmin
          ijyindexmax=npwy-iyregmax+ijymax
          jfirstf=inx(11,jcs)+1
          jlastf=inx(10,jcs)
c
c Now ijxmin=ijxmax=ixreg
c
          jbeg=inx(1,jcs)+1
          Pxj=basdat(11,jbeg)
          Pyj=basdat(12,jbeg)
          Pzj=basdat(13,jbeg)
c
          jend=inx(5,jcs)
          jcssize=jend-jbeg+1
          jtype=inx(12,jcs)
c
          index=1
          do jfunc=jfirstf,jlastf
            denspart(index)=dens(jfunc,ifunc)
            index=index+1
          end do
c
c Here loop for the contractions
c
          do jsh=jbeg,jend  ! loop over the contractions
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
     &    expvalxv,      Pxj,           Lxmin,         eta,
     &    cji,           concoef,       griddens,      xexpvalxv,
     &    xxexpvalxv,    icorexmin1,    icorexmax1,    icorexmin2,
     &    icorexmax2,    ixreg)
c
c    allocate arrays core1y(npwyexp),core2y(npwyexp),core1z(npwzexp)
c
              call mmark
              call getmem(npwyexp,icore1y)
              call getmem(npwyexp,icore2y)
              call getmem(npwzexp,icore1z)
c
              call corecut_1d(
     &    bl(icore1y),   bl(icore2y),   bl(icore1z),   npwyexp,
     & yexpstore(ipyexp),y2expstore(ipyexp),npwzexp,zexpstore(ipzexp),
     &    icoreymin1,    icoreymax1,    icorezmin1,    icorezmax1,
     &    ijzmin,        ijzmax,        ijymin,        ijymax)
c
             call make_gridfunctions_for_derivatives(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Pxj,           Pyj,           Pzj,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,      bl(icore1y),
     &    bl(icore2y),   npwyexp,       bl(icore1z),   ro3,
     &    denspart,      npwzexp)
c
              call retmark  ! deallocate core1y,core2y,core1z
            else
             call calc_expvaluesv(expvalxv,Pxj,Lxmin,eta,
     &         cji,concoef,griddens,xexpvalxv,xxexpvalxv)
c
             call make_gridfunctions_for_derivatives(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Pxj,           Pyj,           Pzj,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,     yexpstore(ipyexp),
     &y2expstore(ipyexp),npwyexp,    zexpstore(ipzexp),ro3,
     &    denspart,      npwzexp)
c
            end if
c
          end do ! end of the loop for the contractions
        end if
c
 120    continue
c
      end do !over jcs
c
      do iy=1,npwy
        do iz=iszmin,iszmax
          izp=iz-npwzdiff
          rojmu=ro3(iz,iy)*ro4(iz,iy)
          valgradx=valgradx+rojmu*ro6(1,izp,iy)
          valgrady=valgrady+rojmu*ro6(2,izp,iy)
          valgradz=valgradz+rojmu*ro6(3,izp,iy)
c
        end do
      end do
c
      return
      end
c
c***********************************************************************
c
      Subroutine add_to_gradients(ncs,inx,ics,natom,valgradx,
     &           valgrady,valgradz,gradx,grady,gradz)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 gradx(natom)
      real*8 grady(natom)
      real*8 gradz(natom)
      dimension inx(12,ncs)
c
      iatom=inx(2,ics)
      gradx(iatom)=gradx(iatom)+valgradx
      grady(iatom)=grady(iatom)+valgrady
      gradz(iatom)=gradz(iatom)+valgradz
c
      return
      end
c
c***********************************************************************
c
      Subroutine calc_mix2_gradients(
     &    ncf,           ro6,           npwx,          npwy,
     &    npwz,          gridranges,    ncs,           sharpness,
     &    ncspl,         dens,          inx,           basdat,
     &    Lxmin,         Lymin,         Lzmin,         ro5,
     &    ixreg,         iyregmin,      iyregmax,      icspl,
     &    ifunc,         sharpovs,      maxovs,        nsharpovs,
     &    npwsz,         iszmin,        ro4,           griddens,
     &    iystoredim,    yexpstore,     y2expstore,    y3expstore,
     &    izstoredim,    zexpstore,     z2expstore,    z3expstore,
     &    ipexpdim,      ipyexpstore,   ipzexpstore,   comprcoef,
     &    natom,         gradx,         grady,         gradz,
     &    icorecut,  isharpgridrange1,isharpgridrange2)
c

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c
c So npwx=1 in this subroutine and npwy=npwregy in this subroutine
c Actual values of the Lxmin and Lymin are coming from outside
c
      integer*4 ro6(npwsz,npwy) ! sharp function
      real*8 ro5(3,npwz,npwy)   ! smooth derivatives
      real*8 ro4(npwz,npwy)     ! Coulomb potencial from the smooth density
      real*8 denspart(10)       ! density mx contribution
      integer gridranges(6,ncs)
      integer sharpness(ncs)
      dimension inx(12,ncs)
      dimension basdat(13,*)
      real*8 dens(ncf,ncf)
      real*8 Lxmin,Lymin,Lzmin
      integer sharpovs(maxovs,ncspl)
      integer nsharpovs(ncspl)
c
      real*8 yexpstore(iystoredim)
      real*8 y2expstore(iystoredim)
      real*8 y3expstore(iystoredim)
      real*8 zexpstore(izstoredim)
      real*8 z2expstore(izstoredim)
      real*8 z3expstore(izstoredim)
      integer ipyexpstore(ipexpdim)
      integer ipzexpstore(ipexpdim)
      real*8 gradx(natom)
      real*8 grady(natom)
      real*8 gradz(natom)
      integer isharpgridrange1(6,ncspl)
      integer isharpgridrange2(6,ncspl)
c     common /big/bl(5000000)
c
      if(icorecut .eq. 1) then
        icorexmin2=isharpgridrange2(1,icspl)
        icorexmax2=isharpgridrange2(2,icspl)
        if(icorexmin2 .gt. ixreg) return !the given x grid does not
                                         ! contribute to the core
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
      npwzdiff=iszmin-1
      cji=1.0d0
      do jcssmooth=1,nsharpovs(icspl)
        jcs=sharpovs(jcssmooth,icspl) ! Loop over the contracted shells
c
        jsharpness=sharpness(jcs)
        if(jsharpness .eq. 3) goto 120 !no sharp is allowed
        jxmin=gridranges(1,jcs)
        if(jxmin .gt. ixreg) goto 120 !the given x grid does not
                                      ! contribute to the function
        jxmax=gridranges(2,jcs)
        if(jxmax .lt. ixreg) goto 120 !the given x grid does not
                                      ! contribute to the function
        jymin=gridranges(3,jcs)
        jyminold=jymin
        jymax=gridranges(4,jcs)
        jzmin=gridranges(5,jcs)
        jzmax=gridranges(6,jcs)
c
c  new part to calculate spherical reagions
c
        call get_y_gridrange(jymin,jymax,jxmin,ixreg,jcsrad,
     &                     jcsradx2,jycenter)
c
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
        if((ijxmin .eq. ijxmax) .and. (ijymin .le. ijymax) .and.
     &     (ijzmin .le. ijzmax)) then
c        write(6,*)"There is common region !"
          ijyindexmin=1+ijymin-iyregmin
          ijyindexmax=npwy-iyregmax+ijymax
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
          index=1
          do jfunc=jfirstf,jlastf
            denspart(index)=dens(jfunc,ifunc)
            index=index+1
          end do
c
c Zero out the ro5 working array
c
          do iy=ijyindexmin,ijyindexmax
            do iz=ijzmin,ijzmax
              ro5(1,iz,iy)=0.0d0
              ro5(2,iz,iy)=0.0d0
              ro5(3,iz,iy)=0.0d0
            end do
          end do
c
          do jsh=jbeg,jend  ! loop over the contractions
c
            eta=basdat(1,jsh)
            if(jtype .eq. 3) then
              concoef=basdat(3,jsh)
              ccs=basdat(2,jsh)/concoef
            else
              concoef=basdat(2,jsh)
            end if
c
c
            if(comprcoef .gt. 0) then
              cji_new=cji/comprcoef    !for the compression to i4
            else
              cji_new=0.0d0
            end if
c
            ipyexp=ipyexpstore(jsh)+(ijymin-jyminold)
            ipzexp=ipzexpstore(jsh)+(ijzmin-jzmin)
            npwyexp=ijyindexmax-ijyindexmin+1
            npwzexp=ijzmax-ijzmin+1
c
            if(icorecut .eq. 1) then
              call calc_expvalx_corecut(
     &    expvalxv,      Px,            Lxmin,         eta,
     &    cji_new,       concoef,       griddens,      xexpvalxv,
     &    xxexpvalxv,    icorexmin1,    icorexmax1,    icorexmin2,
     &    icorexmax2,    ixreg)
c
c    allocate arrays core1y(npwyexp),core2y(npwyexp),core1z(npwzexp)
c
              call mmark
              call getmem(npwyexp,icore1y)
              call getmem(npwyexp,icore2y)
              call getmem(npwzexp,icore1z)
c
              call corecut_1d(
     &    bl(icore1y),   bl(icore2y),   bl(icore1z),   npwyexp,
     & yexpstore(ipyexp),y2expstore(ipyexp),npwzexp,zexpstore(ipzexp),
     &    icoreymin1,    icoreymax1,    icorezmin1,    icorezmax1,
     &    ijzmin,        ijzmax,        ijymin,        ijymax)
c
              call make_gridfunction_derivatives_together(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Px,            Py,            Pz,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,      bl(icore1y),
     &    bl(icore2y),   npwyexp,       bl(icore1z),  npwzexp,
     &    ro5,           denspart,      eta)
c
              call retmark  ! deallocate core1y,core2y,core1z
c
            else
c
              call calc_expvaluesv(expvalxv,Px,Lxmin,eta,
     &         cji_new,concoef,griddens,xexpvalxv,xxexpvalxv)
c
              call make_gridfunction_derivatives_together(
     &    jtype,         ccs,           expvalxv,      xexpvalxv,
     &    xxexpvalxv,    Px,            Py,            Pz,
     &    ijyindexmin,   ijyindexmax,   ijzmin,        ijzmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,     yexpstore(ipyexp),
     &y2expstore(ipyexp),npwyexp,    zexpstore(ipzexp),npwzexp,
     &    ro5,           denspart,      eta)
            endif
c
          end do !over the contractions
c
          valgradx=0.0d0
          valgrady=0.0d0
          valgradz=0.0d0
c
          do iy=ijyindexmin,ijyindexmax
            do iz=ijzmin,ijzmax
              izp=iz-npwzdiff
              rojmu=ro6(izp,iy)*ro4(iz,iy)
              valgradx=valgradx+rojmu*ro5(1,iz,iy)
              valgrady=valgrady+rojmu*ro5(2,iz,iy)
              valgradz=valgradz+rojmu*ro5(3,iz,iy)
            end do
          end do
          iatom=inx(2,jcs)
          gradx(iatom)=gradx(iatom)+valgradx
          grady(iatom)=grady(iatom)+valgrady
          gradz(iatom)=gradz(iatom)+valgradz
c
        end if !for common region
c
 120    continue
c
      end do !over jcssmooth
c
      return
      end
c
c***********************************************************************
c
      Subroutine make_gridfunction_derivatives(
     &    itype,         expvalx,       icomp,         Px,
     &    Py,            Pz,            ro,            iout,
     &    ixmin,         ixmax,         iymin,         iymax,
     &    izmin,         izmax,         istep,         Lxmin,
     &    Lymin,         Lzmin,         npwx,          npwy,
     &    npwz,          griddens,      xexpvalx,      xxexpvalx,
     &    yexpstore,     y2expstore,    y3expstore,    npwyexp,
     &    zexpstore,     z2expstore,    z3expstore,    eta)
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 ro(npwz,npwy,3)
      real*8 Lxmin,Lymin,Lzmin
c
      real*8 yexpstore(npwyexp)
      real*8 y2expstore(npwyexp)
      real*8 y3expstore(npwyexp)
      real*8 zexpstore(npwz)
      real*8 z2expstore(npwz)
      real*8 z3expstore(npwz)
c
      data pi/3.1415926535897932384626433d0/
c
c
      dx=1.0d0/griddens
      xinit=Lxmin - Px - dx
      yinit=Lymin - Py - dx
      zinit=Lzmin - Pz - dx
      ix=1
      x=xinit+dx
      indexy=1
c
      if (itype .eq. 1) then ! s derivatives
c
        valx2=2.0d0*eta*expvalx
        valx1=valx2*x
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 1)) then !Px derivatives
c
        valx1=(2.0d0*eta*x**2-1.0d0)*expvalx
        valx2=2.0d0*eta*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 2)) then !Py derivatives
c
        valx1=2.0d0*eta*xexpvalx
        valx2=expvalx
        valx3=2.0d0*eta*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          yexpy=y*expy
          valxy1=valx1*yexpy
          valxy2=(2.0d0*eta*y**2-1.0d0)*valx2*expy
          valxy3=valx3*yexpy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 2) .and. (icomp .eq. 3)) then !Pz derivatives
c
        etat2=2.0d0*eta
        valx1=etat2*xexpvalx
        valx2=etat2*expvalx
        valx3=expvalx
        do iy=iymin,iymax
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y2expstore(indexy)
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            zexpz=z*expz
            ro(iz,iy,1)=valxy1*zexpz
            ro(iz,iy,2)=valxy2*zexpz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
          indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 1)) then !s of L deriv.
c
        valx2=2.0d0*eta*expvalx
        valx1=valx2*x
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax,istep
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 2)) then !Px of L deriv.
c
        valx1=(2.0d0*eta*x**2-1.0d0)*expvalx
        valx2=2.0d0*eta*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 3)) then !Py of L deriv.
c
        valx1=2.0d0*eta*xexpvalx
        valx2=expvalx
        valx3=2.0d0*eta*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          yexpy=y*expy
          valxy1=valx1*yexpy
          valxy2=(2.0d0*eta*y**2-1.0d0)*valx2*expy
          valxy3=valx3*yexpy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 3) .and. (icomp .eq. 4)) then !Pz of L deriv.
c
        etat2=2.0d0*eta
        valx1=etat2*xexpvalx
        valx2=etat2*expvalx
        valx3=expvalx
        do iy=iymin,iymax
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y2expstore(indexy)
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            zexpz=z*expz
            ro(iz,iy,1)=valxy1*zexpz
            ro(iz,iy,2)=valxy2*zexpz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 1)) then !D1 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=-2.0d0*xexpvalx/sqrt(12.0d0)
        valx2=-2.0d0*expvalx/sqrt(12.0d0)
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          segxy1=etatx2+eta*y**2-1.0d0
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y*expy
          valxy3=valx2*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            segedxyz=segxy1-etat2*z**2
            ro(iz,iy,1)=segedxyz*valxy1*expz
            ro(iz,iy,2)=segedxyz*valxy2*expz
            ro(iz,iy,3)=(segedxyz+3.0d0)*valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 2)) then !D2 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
c      valx1=2.0d0*xexpvalx
c      valx2=2.0d0*expvalx
        valx1=xexpvalx
        valx2=expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          segxy1=etatx2-eta*y**2-1.0d0
          segxy2=segxy1+2.0d0
          segxy3=segxy1+1.0d0
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=segxy1*valxy1*expz
            ro(iz,iy,2)=segxy2*valxy2*expz
            ro(iz,iy,3)=segxy3*valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 3)) then !D3 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=(2.0d0*etatx2-1.0d0)*expvalx
        valx2=xexpvalx
        valx3=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          segxy2=etat2*y**2-1.0d0
          expy=yexpstore(indexy)
          valxy1=valx1*y*expy
          valxy2=segxy2*valx2*expy
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 4)) then !D4 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=(2.0d0*etatx2-1.0d0)*expvalx
        valx3=xexpvalx
        valx2=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx3*expy
          valxy2=valx2*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*z*expz
            ro(iz,iy,2)=valxy2*z*expz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 4) .and. (icomp .eq. 5)) then !D5 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx2=expvalx
        valx3=expvalx
        valx1=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy2=valx2*(etat2*y**2-1.0d0)*expy
          valxy3=valx3*y*expy
          valxy1=valx1*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*z*expz
            ro(iz,iy,2)=valxy2*z*expz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 1)) then !D6,1 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx2=2.0d0*etatx2*expvalx
        valx3=valx2
        valx1=2.0d0*(etatx2-1.0d0)*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy2=valx2*y*expy
          valxy3=valx3*expy
          valxy1=valx1*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=z*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 2)) then !D6,2 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=etat2*xexpvalx
        valx3=etat2*expvalx
        valx2=2.0d0*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          expy=yexpstore(indexy)
          valxy2=(eta*y2-1.0d0)*valx2*y*expy
          valxy3=valx3*y2*expy
          valxy1=valx1*y2*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=z*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 3)) then !D6,3 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=etat2*xexpvalx
        valx3=2.0d0*expvalx
        valx2=etat2*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy2=valx2*y*expy
          valxy3=valx3*expy
          valxy1=valx1*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            expz=zexpstore(indexz)
            expztz2=expz*z2
            ro(iz,iy,1)=valxy1*expztz2
            ro(iz,iy,2)=valxy2*expztz2
            ro(iz,iy,3)=(eta*z2-1.0d0)*valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 4)) then !D6,4 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=(2.0d0*etatx2-1.0d0)*expvalx
        valx2=xexpvalx
        valx3=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          segxy2=etat2*y**2-1.0d0
          expy=yexpstore(indexy)
          valxy1=valx1*y*expy
          valxy2=segxy2*valx2*expy
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 5)) then !D6,5 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx1=(2.0d0*etatx2-1.0d0)*expvalx
        valx3=xexpvalx
        valx2=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx3*expy
          valxy2=valx2*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            expztz=expz*z
            ro(iz,iy,1)=valxy1*expztz
            ro(iz,iy,2)=valxy2*expztz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 5) .and. (icomp .eq. 6)) then !D6,6 deriv.
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        valx2=expvalx
        valx3=expvalx
        valx1=etat2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          expyty=expy*y
          valxy2=valx2*(etat2*y**2-1.0d0)*expy
          valxy3=valx3*expyty
          valxy1=valx1*expyty
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            expztz=expz*z
            ro(iz,iy,1)=valxy1*expztz
            ro(iz,iy,2)=valxy2*expztz
            ro(iz,iy,3)=(etat2*z**2-1.0d0)*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 1)) then !F1 deriv.
c
        segx1=4.0d0*eta*x**2-4.0d0
        valx1=2.0d0*xexpvalx
        fourtx2=4.0d0*x**2
        valx2=expvalx
        segx3=fourtx2*eta+1.0d0
        valx3=2.0d0*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          y4=y2**2
          etaty2=eta*y2
          expy=yexpstore(indexy)
          segxy1=segx1-etaty2
          valxy1=valx1*y*expy
          segxy21=2.0d0*etaty2-1
          segxy22=2.0d0*eta*y4
          segxy23=3.0d0*y2
          valxy2=valx2*expy
          segxy3=segx3-etaty2
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            expz=zexpstore(indexz)
            segxyz1=segxy1-etatz2
            ro(iz,iy,1)=valxy1*segxyz1*expz
            segxyz2=segxy21*(fourtx2-z2)-segxy22+segxy23
            ro(iz,iy,2)=segxyz2*valxy2*expz
            segxyz3=segxy3-etatz2
            ro(iz,iy,3)=z*segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 2)) then !F2 deriv.
c
c 2<->3 change in the helping variables !
c
        etat2=2.0d0*eta
        segx1=4.0d0*eta*x**2-4.0d0
        valx1=2.0d0*xexpvalx
        fourtx2=4.0d0*x**2
        valx2=expvalx
        segx3=fourtx2*eta+1.0d0
        valx3=2.0d0*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          etaty2=eta*y2
          y4=y2**2
          expy=yexpstore(indexy)
          segxy1=segx1-etaty2
          valxy1=valx1*expy
          segxy2=fourtx2-y2
          valxy2=valx2*expy
          segxy3=segx3-etaty2
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            expz=zexpstore(indexz)
            segxyz1=segxy1-etatz2
            ro(iz,iy,1)=valxy1*segxyz1*z*expz
            segxyz2=(etat2*z2-1.0d0)*segxy2-etat2*z2**2+3.0d0*z2
            ro(iz,iy,3)=segxyz2*valxy2*expz
            segxyz3=segxy3-etatz2
            ro(iz,iy,2)=z*segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 3)) then !F3 deriv.
c
c 1<->2 change in the helping variables !
c
        etat4=4.0d0*eta
        segx1=-eta*x**2-4.0d0
        valx1=2.0d0*xexpvalx
        fourtx2=4.0d0*x**2
        segx21=2.0d0*eta*x**2-1.0d0
        segx22=2.0d0*eta*x**4
        segx23=3.0d0*x**2
        valx2=expvalx
        segx3=-eta*x**2+1.0d0
        valx3=2.0d0*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          y4=y2**2
          expy=yexpstore(indexy)
          segxy1=segx1+etat4*y2
          valxy1=valx1*y*expy
          segxy2=4.0d0*y2
          valxy2=valx2*expy
          segxy3=segx3+4.0d0*eta*y2
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            expz=zexpstore(indexz)
            segxyz1=segxy1-etatz2
            ro(iz,iy,2)=valxy1*segxyz1*expz
            segxyz2=(segxy2-z2)*segx21-segx22+segx23
            ro(iz,iy,1)=segxyz2*valxy2*expz
            segxyz3=segxy3-etatz2
            ro(iz,iy,3)=z*segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 4)) then !F4 deriv.
c
c 1<->2 change + 1<->3 change after that in the helping variables  !
c
        etat2=2.0d0*eta
        etat4=4.0d0*eta
        segx1=-eta*x**2-4.0d0
        valx1=2.0d0*expvalx
        fourtx2=4.0d0*x**2
        segx2=-x**2
        valx2=expvalx
        segx3=-eta*x**2+1.0d0
        valx3=2.0d0*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          y4=y2**2
          y2t4=y2*4.0d0
          y2t4teta=eta*y2t4
          expy=yexpstore(indexy)
          segxy1=segx1+y2t4teta
          valxy1=valx1*y*expy
          segxy2=segx2+y2t4
          valxy2=valx2*expy
          segxy3=segx3+y2t4teta
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            expz=zexpstore(indexz)
            segxyz1=segxy1-etatz2
            ro(iz,iy,2)=valxy1*segxyz1*z*expz
            segxyz2=(etat2*z2-1.0d0)*segxy2-etat2*z2**2+3.0d0*z2
            ro(iz,iy,3)=segxyz2*valxy2*expz
            segxyz3=segxy3-etatz2
            ro(iz,iy,1)=z*segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 5)) then !F5 deriv.
c
c 1<->2 change + 2<->3 change after that in the helping variables !
c
        etat4=4.0d0*eta
        segx1=-eta*x**2-4.0d0
        valx1=2.0d0*xexpvalx
        fourtx2=4.0d0*x**2
        segx21=2.0d0*eta*x**2-1.0d0
        segx22=2.0d0*eta*x**4
        segx23=3.0d0*x**2
        valx2=expvalx
        segx3=-eta*x**2+1.0d0
        valx3=2.0d0*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          etaty2=eta*y2
          y4=y2**2
          expy=yexpstore(indexy)
          segxy1=segx1-etaty2
          valxy1=valx1*expy
          segxy2=-y2
          valxy2=valx2*expy
          segxy3=segx3-etaty2
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            z2t4=z2*4.0d0
            z2t4teta=z2t4*eta
            expz=zexpstore(indexz)
            segxyz1=segxy1+z2t4teta
            ro(iz,iy,3)=valxy1*segxyz1*z*expz
            segxyz2=(segxy2+z2t4)*segx21-segx22+segx23
            ro(iz,iy,1)=segxyz2*valxy2*expz
            segxyz3=segxy3+z2t4teta
            ro(iz,iy,2)=segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 6)) then !F6 deriv.
c
c 1<->2 change + 2<->3  + 1<->2 change after that in the helping variables !
c
        etat4=4.0d0*eta
        segx1=-eta*x**2-4.0d0
        valx1=2.0d0*expvalx
        fourtx2=4.0d0*x**2
        segx2=-x**2
        valx2=expvalx
        segx3=-eta*x**2+1.0d0
        valx3=2.0d0*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          etaty2=eta*y2
          expy=yexpstore(indexy)
          segxy1=segx1-etaty2
          valxy1=valx1*y*expy
          segxy21=2.0d0*etaty2-1.0d0
          segxy22=2.0d0*etaty2*y2
          segxy23=3.0d0*y2
          valxy2=valx2*expy
          segxy3=segx3-etaty2
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            z2=z**2
            z2t4=z2*4.0d0
            z2t4teta=z2t4*eta
            expz=zexpstore(indexz)
            segxyz1=segxy1+z2t4teta
            ro(iz,iy,3)=valxy1*segxyz1*z*expz
            segxyz2=(segx2+z2t4)*segxy21-segxy22+segxy23
            ro(iz,iy,2)=segxyz2*valxy2*expz
            segxyz3=segxy3+z2t4teta
            ro(iz,iy,1)=segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 6) .and. (icomp .eq. 7)) then !F7 deriv.
c
        etat2=2.0d0*eta
        segx1=etat2*x**2-1.0d0
        valx1=segx1*expvalx*sqrt(40.0d0)
        valx2=xexpvalx*sqrt(40.0d0)
        valx3=xexpvalx*sqrt(40.0d0)
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          expy=yexpstore(indexy)
          expyty=expy*y
          valxy1=valx1*expyty
          segxy2=etat2*y2-1.0d0
          valxy2=segxy2*valx2*expy
          valxy3=valx3*expyty
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            expztz=expz*z
            ro(iz,iy,1)=valxy1*expztz
            ro(iz,iy,2)=valxy2*expztz
            segxyz3=etat2*z**2-1.0d0
            ro(iz,iy,3)=segxyz3*valxy3*expz
            indexz=indexz+1
          end do
        end do
c
      else if ((itype .eq. 7) .and. (icomp .eq. 1)) then
c
c  1-st component of F10 which is x**3
c
        eta2=2.0d0*eta
        x2=x**2
        valx1=expvalx*(eta2*x2-3.0d0)*x2
        valx2=eta2*x2*xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 2)) then
c
c  2-nd component of F10 which is yx**2
c
        eta2=2.0d0*eta
        x2=x**2
        valx1=(eta2*x2-2.0d0)*xexpvalx
        valx2=xxexpvalx
        valx3=valx2*eta2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*y*expy
          valxy2=valx2*(eta2*y**2-1.0d0)*expy
          valxy3=valx3*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 3)) then
c
c 3-rd component of F10 which is zx**2
c
        eta2=2.0d0*eta
        x2=x**2
        valx1=(eta2*x2-2.0d0)*xexpvalx
        valx3=xxexpvalx
        valx2=valx3*eta2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y*expy
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*z*expz
            ro(iz,iy,2)=valxy2*z*expz
            ro(iz,iy,3)=valxy3*(eta2*z**2-1.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 4)) then
c
c 4-th component of F10 which is xy**2
c
        eta2=2.0d0*eta
        valx1=(eta2*x**2-1.0d0)*expvalx
        valx2=xexpvalx
        valx3=valx2*eta2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          y2=y**2
          valxy1=valx1*y2*expy
          valxy2=valx2*y*(eta2*y2-2.0d0)*expy
          valxy3=valx3*y2*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 5)) then
c
c 5-th component of F10 which is xyz
c
        eta2=2.0d0*eta
        valx1=(eta2*x**2-1.0d0)*expvalx
        valx2=xexpvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*y*expy
          valxy2=valx2*(eta2*y**2-1.0d0)*expy
          valxy3=valx2*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*z*expz
            ro(iz,iy,2)=valxy2*z*expz
            ro(iz,iy,3)=valxy3*(eta2*z**2-1.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 6)) then
c
c 6-th component of F10 which is xz**2
c
        eta2=2.0d0*eta
        valx1=(eta2*x**2-1.0d0)*expvalx
        valx3=xexpvalx
        valx2=valx3*eta2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y*expy
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            z2=z**2
            ro(iz,iy,1)=valxy1*expz*z2
            ro(iz,iy,2)=valxy2*expz*z2
            ro(iz,iy,3)=valxy3*z*(eta2*z2-2.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 7)) then
c
c 7-th component of F10 which is y**3
c
        eta2=2.0d0*eta
        valx1=eta2*xexpvalx
        valx2=expvalx
        valx3=eta2*valx2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          y2=y**2
          y3=y*y2*expy
          valxy1=valx1*y3
          valxy2=valx2*y2*(eta2*y2-3.0d0)*expy
          valxy3=valx3*y3
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*expz
            ro(iz,iy,2)=valxy2*expz
            ro(iz,iy,3)=valxy3*z*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 8)) then
c
c 8-th component of F10 which is zy**2
c
        eta2=2.0d0*eta
        valx1=eta2*xexpvalx
        valx2=expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          y2=y**2
          valxy1=valx1*y2*expy
          valxy2=valx2*y*(eta2*y2-2.0d0)*expy
          valxy3=valx2*y2*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            ro(iz,iy,1)=valxy1*z*expz
            ro(iz,iy,2)=valxy2*z*expz
            ro(iz,iy,3)=valxy3*(eta2*z**2-1.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 9)) then
c
c 9-th component of F10 which is yz**2
c
        eta2=2.0d0*eta
        valx1=eta2*xexpvalx
        valx2=expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*y*expy
          valxy2=valx2*(eta2*y**2-1.0d0)*expy
          valxy3=valx2*y*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            z2=z**2
            ro(iz,iy,1)=valxy1*z2*expz
            ro(iz,iy,2)=valxy2*z2*expz
            ro(iz,iy,3)=valxy3*z*(eta2*z2-2.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
      else if ((itype .eq. 7) .and. (icomp .eq. 10)) then
c
c 10-th component of F10 which is z**3
c
        eta2=2.0d0*eta
        valx1=eta2*xexpvalx
        valx3=expvalx
        valx2=eta2*valx3
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy2=valx2*y*expy
          valxy3=valx3*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
            z2=z**2
            z3=z*z2*expz
            ro(iz,iy,1)=valxy1*z3
            ro(iz,iy,2)=valxy2*z3
            ro(iz,iy,3)=valxy3*z2*(eta2*z2-3.0d0)*expz
            indexz=indexz+1
          enddo
        enddo
c
c
      else
        call nerror(1,'make_gridfunction_derivatives',
     $    ' FTC: Can only handle S, P, D and F functions',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
      Subroutine make_gridfunctions_for_derivatives(
     &    itype,         ccs,           expvalx,       xexpvalx,
     &    xxexpvalx,     Px,            Py,            Pz,
     &    iymin,         iymax,         izmin,         izmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,      yexpstore,
     &    y2expstore,    npwyexp,       zexpstore,     ro5,
     &    densmx,        npwzexp)
c
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 ro5(npwz,npwy)    ! sum over d(mu,nu)*g(mu), input and output
      real*8 Lxmin,Lymin,Lzmin
c
      real*8 yexpstore(npwyexp)
      real*8 y2expstore(npwyexp)
      real*8 zexpstore(npwzexp)
      real*8 densmx(10)        ! density matrix elements, input
c
      data pi/3.1415926535897932384626433d0/
c
c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' On entry to <make_gridfunctions_for_derivatives>'
cc      write(6,*) ' itype:',itype,' ccs:',ccs,' expvalx:',expvalx
cc      write(6,*) ' xexpvalx:',xexpvalx,' xxexpvalx:',xxexpvalx
cc      write(6,*) ' Px:',px,' Py:',py,' Pz:',pz
cc      write(6,*) ' iymin:',iymin,' iymax:',iymax,' izmin:',izmin,
cc     $           ' izmax:',izmax
cc      write(6,*) ' Lxmin:',lxmin,' Lymin:',lymin,' Lzmin:',lzmin
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' griddens:',griddens,' npwyexp:',npwyexp,
cc     $           ' npwzexp:',npwzexp
cc      write(6,*) ' yexpstore and y2expstore arrays are:'
cc      do i=1,npwyexp
cc      write(6,*) i,'  ',yexpstore(i),y2expstore(i)
cc      enddo
cc      write(6,*) ' zexpstore array is:'
cc      do i=1,npwzexp
cc      write(6,*) i,'  ',zexpstore(i)
cc      enddo
cc      write(6,*) ' densmx array is:'
cc      do i=1,10
cc      write(6,*) i,'  ',densmx(i)
cc      enddo
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dx=1.0d0/griddens
      xinit=Lxmin - Px - dx
      yinit=Lymin - Py - dx
      zinit=Lzmin - Pz - dx
      ix=1
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
            fsm=valxy*zexpstore(indexz)
            indexz=indexz+1
            ro5(iz,iy)=ro5(iz,iy)+d1*fsm
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
          valxypy=expvalx*y2expstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
c
            fsm=valxypx*expz
            denssumm1=d1*fsm
c now       px is done
            fsm=valxypy*expz
            denssumm2=d2*fsm
c now       py is done
            fsm=z*valxypz*expz
            denssumm3=d3*fsm
c now       pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3
            indexz=indexz+1
c
          end do
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
          valxypy=expvalx*y2expstore(indexy)
          valxypz=expvalx*expy
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
            expz=zexpstore(indexz)
c
            fsm=valxys*expz
            denssumm1=d1*fsm
c now       s is done
            fsm=valxypx*expz
            denssumm2=d2*fsm
c now       px is done
            fsm=valxypy*expz
            denssumm3=d3*fsm
c now       py is done
            fsm=z*valxypz*expz
            denssumm4=d4*fsm
c now       pz is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4
            indexz=indexz+1
c
          end do
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
        x=xinit+float(ix)*dx
        r2x=x**2
        do iy=iymin,iymax
c
          y=yinit+float(iy)*dx
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
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
            a = 2.0d0*z**2 - r2xy
c
c
            fsm1=zexpstore(indexz)
            fsm=a*valxy1*fsm1
            denssumm1=d1*fsm
c now       D1 is done
            fsm=valxy2*fsm1
            denssumm2=d2*fsm
c now       D2 is done
            fsm=valxy3*fsm1
            denssumm3=d3*fsm
c now       D3 is done
            fsm=z*valxy4*fsm1
            denssumm4=d4*fsm
c now       D4 is done
            fsm=z*valxy5*fsm1
            denssumm5=d5*fsm
c now       D5 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5
            indexz=indexz+1
c
          end do
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
          y=yinit+float(iy)*dx
          valxy1=xxexpvalx*expy
          valxy2=expvalx*expy*y**2
          valxy3=expvalx*expy
          valxy5=xexpvalx*expy
          valxy4=valxy5*y
          valxy6=y*valxy3
c
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
            fsm1=zexpstore(indexz)
            ztfsm1=z*fsm1
c
            fsm=valxy1*fsm1
            denssumm1=d1*fsm
c now       D61 is done
            fsm=valxy2*fsm1
            denssumm2=d2*fsm
c now       D62 is done
            fsm=valxy3*ztfsm1*z
            denssumm3=d3*fsm
c now       D63 is done
            fsm=valxy4*fsm1
            denssumm4=d4*fsm
c now       D64 is done
            fsm=valxy5*ztfsm1
            denssumm5=d5*fsm
c now       D65 is done
            fsm=valxy6*ztfsm1
            denssumm6=d6*fsm
c now       D66 is done
            ro5(iz,iy)=ro5(iz,iy)+denssumm1+denssumm2+denssumm3+
     &                            denssumm4+denssumm5+denssumm6
            indexz=indexz+1
c
          end do
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
          y=yinit+float(iy)*dx
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
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
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
            denssumm1=d1*a1
            denssumm2=d2*a2
            denssumm3=d3*a3
            denssumm4=d4*a4
            denssumm5=d5*a5
            denssumm6=d6*a6
            denssumm7=d7*a7
c
            ro5(iz,iy)=ro5(iz,iy)+(denssumm1+denssumm2+denssumm3+
     &                  denssumm4+denssumm5+denssumm6+denssumm7)*vals
            indexz=indexz+1
c
          end do
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
        x=xinit+float(ix)*dx
        x2=x**2
        x3=x2*x
c
        do iy=iymin,iymax
c
          expy=yexpstore(indexy)
          y=yinit+float(iy)*dx
          y2=y**2
          y3=y2*y
          xy=x*y
          x2y=x*xy
          xy2=y*xy
          valxy=expvalx*expy
c
          indexy=indexy+1
          indexz=1
c
          do iz=izmin,izmax
c
            z=zinit+float(iz)*dx
            z2=z**2
            fsm1=zexpstore(indexz)
            vals=valxy*fsm1
c
            denssumm1=d1*x3
            denssumm2=d2*x2y
            denssumm3=d3*x2*z
            denssumm4=d4*xy2
            denssumm5=d5*xy*z
            denssumm6=d6*x*z2
            denssumm7=d7*y3
            denssumm8=d8*y2*z
            denssumm9=d9*y*z2
            denssumm10=d10*z2*z
c
            ro5(iz,iy)=ro5(iz,iy)+(denssumm1+denssumm2+denssumm3
     &                           + denssumm4+denssumm5+denssumm6
     &                           + denssumm7+denssumm8+denssumm9
     &                           + denssumm10)*vals
            indexz=indexz+1
c
          end do
        end do
c
      else
        call nerror(1,'make_gridfunction_for_derivatives',
     $    ' FTC: Can only handle S, P, D and F functions',0,0)
      end if
c
      return
      end
c
c**********************************************************************
c
      Subroutine make_gridfunction_derivatives_together(
     &    itype,         ccs,           expvalx,       xexpvalx,
     &    xxexpvalx,     Px,            Py,            Pz,
     &    iymin,         iymax,         izmin,         izmax,
     &    Lxmin,         Lymin,         Lzmin,         npwx,
     &    npwy,          npwz,          griddens,      yexpstore,
     &    y2expstore,    npwyexp,       zexpstore,     npwzexp,
     &    ro5,           densmx,        eta)
c
c
      IMPLICIT REAL*8(A-H,O-Z)
      real*8 ro5(3,npwz,npwy) ! sum over d(mu,nu)*dg(mu)/dxi;
                              ! xi=x,y,z; input and output
      real*8 Lxmin,Lymin,Lzmin
c
      real*8 yexpstore(npwyexp)
      real*8 y2expstore(npwyexp)
      real*8 zexpstore(npwzexp)
      real*8 densmx(10)       ! density matrix elements, input
c
      data pi/3.1415926535897932384626433d0/
c
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc      write(6,*) ' On entry to <make_gridfunction_derivatives_together>'
cc      write(6,*) ' itype:',itype,' ccs:',ccs,' expvalx:',expvalx
cc      write(6,*) ' xexpvalx:',xexpvalx,' xxexpvalx:',xxexpvalx
cc      write(6,*) ' Px:',px,' Py:',py,' Pz:',pz
cc      write(6,*) ' iymin:',iymin,' iymax:',iymax,' izmin:',izmin,
cc     $           ' izmax:',izmax
cc      write(6,*) ' Lxmin:',lxmin,' Lymin:',lymin,' Lzmin:',lzmin
cc      write(6,*) ' npwx:',npwx,' npwy:',npwy,' npwz:',npwz
cc      write(6,*) ' griddens:',griddens,' npwyexp:',npwyexp,
cc     $           ' npwzexp:',npwzexp
cc      write(6,*) ' yexpstore and y2expstore arrays are:'
cc      do i=1,npwyexp
cc      write(6,*) i,'  ',yexpstore(i),y2expstore(i)
cc      enddo
cc      write(6,*) ' zexpstore array is:'
cc      do i=1,npwzexp
cc      write(6,*) i,'  ',zexpstore(i)
cc      enddo
cc      write(6,*) ' densmx array is:'
cc      do i=1,10
cc      write(6,*) i,'  ',densmx(i)
cc      enddo
cc      write(6,*) ' eta is:',eta
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      dx=1.0d0/griddens
      xinit=Lxmin - Px - dx
      x=xinit+dx
      yinit=Lymin - Py - dx
      zinit=Lzmin - Pz - dx
      ix=1
      indexy=1
c
      if (itype .eq. 1) then
c
c  s type function
c
        d1=densmx(1)
        valx2=2.0d0*eta*expvalx*d1
        valx1=valx2*x
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          expy=yexpstore(indexy)
          valxy1=valx1*expy
          valxy3=valx2*expy
          valxy2=valxy3*y
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            ro5(1,iz,iy)=ro5(1,iz,iy)+valxy1*expz
            ro5(2,iz,iy)=ro5(2,iz,iy)+valxy2*expz
            ro5(3,iz,iy)=ro5(3,iz,iy)+valxy3*z*expz
            indexz=indexz+1
          end do
        end do
c
      else if (itype .eq. 2) then
c
c P functions
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
c
        etat2=eta*2.0d0
        valxs=etat2*expvalx
        valx1=(etat2*x**2-1.0d0)*expvalx
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          xty=x*y
          expy=yexpstore(indexy)
          expvalxy=expvalx*expy
          d3texpvalxy=d3*expvalxy
          valxys=valxs*expy
          valxy1=valx1*expy
          d1tvalxy1=d1*valxy1
          valxy5=(etat2*y**2-1.0d0)*expy*expvalx
          d2tvalxy5=d2*valxy5
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            valxyzs=valxys*expz
            valxyzstxty=valxyzs*xty
            valxyzstxtz=valxyzs*x*z
            valxyzstytz=valxyzs*y*z
c
            ro5(1,iz,iy)=ro5(1,iz,iy) + d1tvalxy1*expz +
     &                   d2*valxyzstxty + d3*valxyzstxtz
            ro5(2,iz,iy)=ro5(2,iz,iy) + d1*valxyzstxty +
     &                   d2tvalxy5*expz + d3*valxyzstytz
            ro5(3,iz,iy)=ro5(3,iz,iy) + d1*valxyzstxtz +
     &              d2*valxyzstytz +(etat2*z**2-1.0d0)*d3texpvalxy*expz
            indexz=indexz+1
          end do
        end do
c
      else if (itype .eq. 3) then
c
        d1=densmx(1)
        d2=densmx(2)
        d3=densmx(3)
        d4=densmx(4)
c
        etat2=eta*2.0d0
        valxs=etat2*expvalx
        valx1=(etat2*x**2-1.0d0)*expvalx
        d1tx=d1*x
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          xty=x*y
          d1ty=d1*y
          expy=yexpstore(indexy)
          expvalxy=expvalx*expy
          d4texpvalxy=d4*expvalxy
          valxys=valxs*expy
          valxy1=valx1*expy
          d2tvalxy1=d2*valxy1
          valxy5=(etat2*y**2-1.0d0)*expy*expvalx
          d3tvalxy5=d3*valxy5
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            valxyzs=valxys*expz
            valxyzstxty=valxyzs*xty
            valxyzstxtz=valxyzs*x*z
            valxyzstytz=valxyzs*y*z
c
            ro5(1,iz,iy)=ro5(1,iz,iy) + d2tvalxy1*expz +
     &                   d3*valxyzstxty + d4*valxyzstxtz + d1tx*valxyzs
            ro5(2,iz,iy)=ro5(2,iz,iy) + d2*valxyzstxty +
     &                   d3tvalxy5*expz + d4*valxyzstytz + d1ty*valxyzs
            ro5(3,iz,iy)=ro5(3,iz,iy) + d2*valxyzstxtz +
     &            d3*valxyzstytz + (etat2*z**2-1.0d0)*d4texpvalxy*expz+
     &                                  d1*z*valxyzs
            indexz=indexz+1
          end do
        end do
c
      else if (itype .eq. 4) then
c
c D functions
c
        d1=densmx(1)/sqrt(12.0d0)
        d2=densmx(2)/2.0d0
        d3=densmx(3)
        d4=densmx(4)
        d5=densmx(5)
c
        etat2=eta*2.0d0
        etatx2=eta*x**2
        segx1=2.0d0*etatx2-1.0d0
        d3tsegx1=d3*segx1
        d4tsegx1=d4*segx1
        x2=x**2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          etaty2=eta*y2
          segxy5=etatx2+etaty2-1.0d0
          etaty2p1=etaty2+1.0d0
          segxy3=d2*(etatx2-etaty2p1)
          etaty2m1=etaty2-1.0d0
          segxy4=d2*(etatx2-etaty2m1)
          etaty2t2m1=2.0d0*etaty2m1+1.0d0
          d3tetaty2t2m1=d3*etaty2t2m1
          d5tetaty2t2m1=d5*etaty2t2m1
          expy=yexpstore(indexy)
          valxys=expvalx*expy
          valxystx=valxys*x
          valxystxt2=2.0d0*valxystx
          valxysty=valxys*y
          valxystyt2=valxysty*2.0d0
          segxy1=valxystx*y*etat2
          segxy2=etatx2-etaty2
          d2tsegxy2=d2*segxy2
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            etatz2t2=2.0d0*etatz2
            valxyzs=valxys*expz
            valxyzstx=x*valxyzs
            valxyzstxt2=2.0d0*valxyzstx
            valxyzsty=y*valxyzs
            valxyzstyt2=2.0d0*valxyzsty
            valxyzstz=z*valxyzs
            valxyzstzt2=2.0d0*valxyzstz
            seg1=segxy5-etatz2t2
            d1tseg1=d1*seg1
            seg2=seg1+3.0d0
            seg3=segxy1*z*expz
            seg4=etatz2t2-1.0d0
c
            ro5(1,iz,iy)=ro5(1,iz,iy) - d1tseg1*valxyzstxt2 +
     &                                  segxy3*valxyzstxt2 +
     &                        d3tsegx1*valxyzsty + d4tsegx1*valxyzstz +
     &                                  d5*seg3
c
            ro5(2,iz,iy)=ro5(2,iz,iy) - d1tseg1*valxyzstyt2 +
     &                                  segxy4*valxyzstyt2 +
     &                              d3tetaty2t2m1*valxyzstx + d4*seg3 +
     &                                  d5tetaty2t2m1*valxyzstz
c
            ro5(3,iz,iy)=ro5(3,iz,iy) - d1*seg2*valxyzstzt2 +
     &                               d2tsegxy2*valxyzstzt2 + d3*seg3 +
     &                          d4*seg4*valxyzstx + d5*seg4*valxyzsty
c
            indexz=indexz+1
          end do
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
        etat2=eta*2.0d0
        etatx2=eta*x**2
        d1tetatx2=d1*etatx2
        segx1=2.0d0*etatx2-1.0d0
        segx1v=2.0d0*etatx2-2.0d0
        d1tsegx1=d1*segx1v
        d4tsegx1=d4*segx1
        d5tsegx1=d5*segx1
        x2=x**2
        do iy=iymin,iymax
          y=yinit+float(iy)*dx
          y2=y**2
          d2ty2=d2*y2
          etaty2=eta*y2
          d2tetaty2=d2*etaty2
          etaty2p1=etaty2+1.0d0
          etaty2m1=etaty2-1.0d0
          d2tetaty2m1=d2*etaty2m1
          etaty2t2m1=2.0d0*etaty2m1+1.0d0
          d4tetaty2t2m1=d4*etaty2t2m1
          d6tetaty2t2m1=d6*etaty2t2m1
          expy=yexpstore(indexy)
          valxys=expvalx*expy
          segxy1=valxys*x*y*etat2
          segxy2=etatx2-etaty2
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            etatz2t2=2.0d0*etatz2
            seg4=etatz2t2-1.0d0
            seg3=segxy1*z*expz
            valxyzs=valxys*expz
            valxyzstx=x*valxyzs
            valxyzstxt2=2.0d0*valxyzstx
            valxyzsty=y*valxyzs
            valxyzstyt2=2.0d0*valxyzsty
            valxyzstz=z*valxyzs
            valxyzstzt2=2.0d0*valxyzstz
c
            ro5(1,iz,iy)=ro5(1,iz,iy) + d1tsegx1*valxyzstx +
     &                                  (d2ty2+d3*z2)*eta*valxyzstxt2 +
     &                                  d4tsegx1*valxyzsty +
     &                                  d5tsegx1*valxyzstz + d6*seg3
c
            ro5(2,iz,iy)=ro5(2,iz,iy) + d1tetatx2*valxyzstyt2 +
     &                                  d2tetaty2m1*valxyzstyt2 +
     &                                  d3*etatz2*valxyzstyt2 +
     &                            d4tetaty2t2m1*valxyzstx + d5*seg3 +
     &                                  d6tetaty2t2m1*valxyzstz
c
            ro5(3,iz,iy)=ro5(3,iz,iy) + d1tetatx2*valxyzstzt2 +
     &                                  d2tetaty2*valxyzstzt2 +
     &                                  d3*valxyzstzt2*(etatz2-1.0d0) +
     &                                  d4*seg3+
     &                           d5*seg4*valxyzstx + d6*seg4*valxyzsty
c
            indexz=indexz+1
          end do
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
        d7=densmx(7)*sqrt(40.0d0)*0.5d0 !this 0.5 makes
                                        !the computation little faster !
c
        etat2=2.0d0*eta
        etat4=4.0d0*eta
        xt2=2.0d0*x
        x2=x**2
        x2t3=3.0d0*x2
        etatx2=eta*x2
        segx3=2.0d0*etatx2-1.0d0
        segx4=2.0d0*etatx2*x2
        segx5=x2t3-segx4
        segx1=-etatx2-4.0d0
        x2t4=4.0d0*x2
        etatx2t4=eta*x2t4
        segx2=etatx2t4-4.0d0
c
        do iy=iymin,iymax
          expy=yexpstore(indexy)
          y=yinit+float(iy)*dx
          yt2=2.0d0*y
          y2=y**2
          etaty2=eta*y2
          segxy11=segx1+etat4*y2
          segxy12=segx1-etaty2
          segxy2=segx2-etaty2
          x2t4my2=x2t4-y2
          y2t4=4.0d0*y2
          y2t4mx2=y2t4-x2
          etat2ty2=etat2*y2
          segxy3=etat2ty2-1.0d0
          segxy4=etat2ty2*y2
          segxy5=3.0d0*y2-segxy4
          valxys=expvalx*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            z2=z**2
            etatz2=eta*z2
            etatz2t2=2.0d0*etatz2
            etatz2t2m1=etatz2t2-1.0d0
            z2t4=4.0d0*z2
            seg1=segxy2-etatz2
            seg2=seg1+5.0d0
            seg3=segxy11-etatz2
            seg4=seg3+5.0d0
            seg5=segxy12+4.0d0*etatz2
            seg6=seg5+5.0d0
            seg7=etatz2t2-1.0d0
            seg8=etatz2t2*z2
            valxyzs=valxys*expz
            valxyzstxt2=xt2*valxyzs
            valxyzstxtyt2=y*valxyzstxt2
            valxyzstxtzt2=z*valxyzstxt2
            valxyzstytzt2=yt2*z*valxyzs
            seg9=segx3*valxyzs
            seg10=valxyzs*segx5
            seg11=3.0d0*z2-etat2*z2**2
c
            ro5(1,iz,iy) = ro5(1,iz,iy) + d1*seg1*valxyzstxtyt2 +
     &                                    d2*seg1*valxyzstxtzt2 +
     &                                    d3*seg9*(y2t4-z2)+d3*seg10 +
     &                                    d4*seg4*valxyzstxtzt2 +
     &                                    d5*seg9*(z2t4-y2)+d5*seg10 +
     &                                    d6*seg6*valxyzstxtyt2 +
     &                                    d7*segx3*valxyzstytzt2
c
            ro5(2,iz,iy) = ro5(2,iz,iy) +
     &                       d1*valxyzs*(segxy5+segxy3*(x2t4-z2))+
     &                                   d2*seg2*valxyzstytzt2 +
     &                                   d3*seg3*valxyzstxtyt2 +
     &                                   d4*seg3*valxyzstytzt2 +
     &                                   d5*seg6*valxyzstxtyt2 +
     &                            d6*valxyzs*(segxy5+segxy3*(z2t4-x2))+
     &                                   d7*segxy3*valxyzstxtzt2
c
            ro5(3,iz,iy) = ro5(3,iz,iy) +d1*seg2*valxyzstytzt2 +
     &                     d2*valxyzs*(etatz2t2m1*x2t4my2+seg11)+
     &                                   d3*seg4*valxyzstxtzt2 +
     &                     d4*valxyzs*(etatz2t2m1*y2t4mx2+seg11) +
     &                                   d5*seg5*valxyzstxtzt2 +
     &                                   d6*seg5*valxyzstytzt2 +
     &                                   d7*etatz2t2m1*valxyzstxtyt2
c
            indexz=indexz+1
          end do
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
        eta2=2.0d0*eta
        eta2x=eta2*x
        eta2x2=eta2x*x
        x2=x**2
        x3=x2*x
        c1=(eta2x2-1.0d0)
        c2=(eta2x2-2.0d0)*x
        segx1x=(eta2x2-3.0d0)*x2
        segx1y=x3*eta2
        segx1z=x3*eta2
        segx2x=c2
        segx2y=x2
        segx2z=eta2x2
        segx3x=c2
        segx3y=eta2x2
        segx3z=x2
        segx4x=c1
        segx4y=x
        segx4z=eta2x
        segx5x=c1
        segx5y=x
        segx5z=x
        segx6x=c1
        segx6y=eta2x
        segx6z=x
        segx7x=eta2x
        segx7y=1.0d0
        segx7z=eta2
        segx8x=eta2x
        segx8y=1.0d0
        segx8z=1.0d0
        segx9x=eta2x
        segx9y=1.0d0
        segx9z=1.0d0
        segx10x=eta2x
        segx10y=eta2
        segx10z=1.0d0
c
        do iy=iymin,iymax
          expy=yexpstore(indexy)
          y=yinit+float(iy)*dx
          y2=y**2
          y3=y2*y
          c1=(eta2*y2-1.0d0)
          c2=(eta2*y2-2.0d0)*y
          segy1x=segx1x
          segy1y=segx1y*y
          segy1z=segx1z
          segy2x=segx2x*y
          segy2y=segx2y*c1
          segy2z=segx2z*y
          segy3x=segx3x
          segy3y=segx3y*y
          segy3z=segx3z
          segy4x=segx4x*y2
          segy4y=segx4y*c2
          segy4z=segx4z*y2
          segy5x=segx5x*y
          segy5y=segx5y*c1
          segy5z=segx5z*y
          segy6x=segx6x
          segy6y=segx6y*y
          segy6z=segx6z
          segy7x=segx7x*y3
          segy7y=segx7y*(eta2*y2-3.0d0)*y2
          segy7z=segx7z*y3
          segy8x=segx8x*y2
          segy8y=segx8y*c2
          segy8z=segx8z*y2
          segy9x=segx9x*y
          segy9y=segx9y*c1
          segy9z=segx9z*y
          segy10x=segx10x
          segy10y=segx10y*y
          segy10z=segx10z
          valxy=expvalx*expy
          indexy=indexy+1
          indexz=1
          do iz=izmin,izmax
            expz=zexpstore(indexz)
            z=zinit+float(iz)*dx
            z2=z**2
            z3=z2*z
            c1=(eta2*z2-1.0d0)
            c2=(eta2*z2-2.0d0)*z
            segz1x=segy1x
            segz1y=segy1y
            segz1z=segy1z*z
            segz2x=segy2x
            segz2y=segy2y
            segz2z=segy2z*z
            segz3x=segy3x*z
            segz3y=segy3y*z
            segz3z=segy3z*c1
            segz4x=segy4x
            segz4y=segy4y
            segz4z=segy4z*z
            segz5x=segy5x*z
            segz5y=segy5y*z
            segz5z=segy5z*c1
            segz6x=segy6x*z2
            segz6y=segy6y*z2
            segz6z=segy6z*c2
            segz7x=segy7x
            segz7y=segy7y
            segz7z=segy7z*z
            segz8x=segy8x*z
            segz8y=segy8y*z
            segz8z=segy8z*c1
            segz9x=segy9x*z2
            segz9y=segy9y*z2
            segz9z=segy9z*c2
            segz10x=segy10x*z3
            segz10y=segy10y*z3
            segz10z=segy10z*(eta2*z2-3.0d0)*z2
            valxyz=valxy*expz
c
            ro5(1,iz,iy) = ro5(1,iz,iy) +(d1*segz1x + d2*segz2x
     &                      + d3*segz3x + d4*segz4x + d5*segz5x
     &                      + d6*segz6x + d7*segz7x + d8*segz8x
     &                      + d9*segz9x + d10*segz10x)*valxyz
c
            ro5(2,iz,iy) = ro5(2,iz,iy) +(d1*segz1y + d2*segz2y
     &                      + d3*segz3y + d4*segz4y + d5*segz5y
     &                      + d6*segz6y + d7*segz7y + d8*segz8y
     &                      + d9*segz9y + d10*segz10y)*valxyz
c
            ro5(3,iz,iy) = ro5(3,iz,iy) +(d1*segz1z + d2*segz2z
     &                      + d3*segz3z + d4*segz4z + d5*segz5z
     &                      + d6*segz6z + d7*segz7z + d8*segz8z
     &                      + d9*segz9z + d10*segz10z)*valxyz
c
            indexz=indexz+1
          end do
        end do
c
      else
        call nerror(1,'make_gridfunction_derivatives_together',
     $    ' FTC: Can only handle S, P, D and F functions',0,0)
      end if
c
      return
      end
c
c***********************************************************************
c
      Subroutine Precalc_expfuncs(expfuncx,expfuncy,expfuncz,npwx,npwy,
     &                            npwz,Lx,Ly,Lz)
c
      IMPLICIT REAL*8(A-H,O-Z)
c
      real*8 expfuncx(npwx)
      real*8 expfuncy(npwy)
      real*8 expfuncz(npwz)
      real*8 Lx,Ly,Lz
      data pi/3.1415926535897932384626433d0/
c
      ax = 2.0d0 * pi / Lx
      ay = 2.0d0 * pi / Ly
      az = 2.0d0 * pi / Lz
c
      nx=npwx/2+2
      ny=npwy/2+2
      nz=npwz/2+2
c
      do jkx=1,npwx
        kx=(jkx-1)-npwx*int(jkx/nx)
        expfuncx(jkx)=exp(-0.25d0*(ax*kx)**2)
      end do
c
      do jky=1,npwy
        ky=(jky-1)-npwy*int(jky/ny)
        expfuncy(jky)=exp(-0.25d0*(ay*ky)**2)
      end do
c
      do jkz=1,npwz
        kz=(jkz-1)-npwz*int(jkz/nz)
        expfuncz(jkz)=exp(-0.25d0*(az*kz)**2)
      end do
c
      return
      end
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine prnt_inx(ncs,inx)
      implicit integer(a-z)
      dimension inx(12,ncs)
      do i=1,12
      write(6,*) (inx(i,j),j=1,ncs)
      enddo
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

c The routines used in the slaves are here
c slave:         the slave driver (We employ humane treatment)
c                NO LONGER - we are now without mercy    ! JB March 2001
c do_scf:        do scf integrals
c do_post:       do force/nmr integrals
c do_prop:       do nuclear properties
c do_mp2:        do MP2 energy
c do_hess:       do analytical Hessian
c do_cosmo:      do COSMO force
c do_ccsd:       do post-SCF energies
c para_get_info: unpacks the basic data required
c
c=================================================
c
      subroutine slave
c
c=================================================
      
      use memory
      use newpara
      
      integer myjob,ioffset
      integer incore,indisk
      integer*4 license,chklicense
      character*13 statfile
      character*256 pqs_root,scrfile,prodstr,scrdir,jobname
c     common /mmanag/lcore,nreq,maxmem,mark,marks(100),nadr(2)
c
c close stdin (/dev/null)
      close(5)
c
c
c get memory data from master
c
      call para_slave_mem(lcore,incore,indisk,pqs_root,scrfile,
     $                    scrdir,jobname)
      call allocmem(lcore,ioffset)
      call setival('lcore',lcore)
      call setival('incore',incore)
      call setival('indisk',indisk)
      call setival('ioffset',ioffset)
      call setival('nslv',nslv)
c -- set pqs_root and update environment in slaves
      call setchval('pqs_root',pqs_root)
      call set_pqs_env
c -- set jobname and scratch file root
      call setchval('scrdir',scrdir)
      call setchval('jobname',jobname)
c
c append slave id to scratch file name
c
      lens=len_trim(scrfile)
      write(scrfile(lens+1:lens+5),'(''_'',z4.4)')MY_GID
      call setchval('scrf',scrfile)
c
c -- check for valid license
c
      call getchval('prodstr',prodstr)
      license = chklicense(prodstr)
      If(License.lt.1) Then
        call nerror(1,'slave','No valid license found',0,0)
        stop
      EndIf
      If(License.lt.NSLV) Then
        call nerror(2,'slave','Too many processors',NSLV,license)
        stop
      EndIf
c
c -- reserve the integral storage area
      if(incore.ne.0) then
        call getmem(incore,incorex)
        call setival('incorex',incorex)
      endif
c
c  open krzystof timing file (send it to /dev/null)
c
      call setival('tim',0)
      call open_tim
c
c open unit 98 for use in the slaves
c name is unique for the slaves
c
c      statfile(1:5)='para.'
c      write (statfile(6:13),fmt='(I8.8)') MY_GID
c      open(98,file=statfile)
c
c  check the scratch directory
c
      call checkscratch
c
c This is the main thing - get jobtype and do your work
c and do not eat up memory!!!
c
      do
         call mmark
c        call immark
         call para_recv_bcast(myjob,TxJobType)
         if (myjob.eq.TxDoScf) then             ! SCF
           call do_scf
         else if (myjob.eq.TxDoPost) then       ! Gradient/NMR
           call do_post
         else if (myjob.eq.TxDoProp) then       ! nuclear properties
           call do_prop
         else if (myjob.eq.TxDoMP2) then        ! MP2 energies
           call do_mp2
         else if (myjob.eq.TxDoHess) then       ! Analytical Hessian
           call do_hess
         else if (myjob.eq.TxDoCosmoForc) then  ! COSMO gradiemt
           call do_cosmoforc
         else if (myjob.eq.TxDoCCSD) then       ! post-SCF energies
           call do_ccsd
c -- we are ready, shut down
         else if (myjob.eq.TxFinish) then
           exit
         end if
c        call retimark
         call retmark
         call matreset
      end do
      close(91)
      end
c======================================================
c
      subroutine do_scf
c
c======================================================
      
      use memory
      use newpara
      
      implicit real*8 (a-h,o-z)
      character*24 datestr
      character scftype*11, ch3*3
      character*256 conel,scrfile
      Logical rhf
      Logical first
      dimension xintxx(9)
c
c  FTC stuff
c
      logical initialize,lpwave
      character*256 filename,filname1
      real*8 Lxmin,Lymin,Lzmin,Lxe,Lye,Lze
      parameter (isharpgrd=51)            ! unit number for FTC I/O
      parameter (Half=0.5d0)
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l !P6 polinom coefs
      Dimension NFunc(7)
      data NFunc/ 1, 3, 4, 5, 6, 7,10/
      data llsemi/0/
c     common /big/ bl(300)
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
c     common /intbl/ maxsh,inx(100)
      common /datstore/ thres_stored,isto,isecond,ito
      common /lindvec/ lind,idensp
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c
c ..............................................................
      integer*8 xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
c ..............................................................
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     $                ixgd,iwt,iscr,icofm,nscr
c
      COMMON /DFTCoeff/ aXX(18) ! contains 18 density functional coefficients
c ..............................................................
c
      initialize = .False.
      call setival('ftc0',0)     ! JB - NEED TO CHECK THIS
c
c -- set memory mark
      call mmark
c
c get the corresponding host name
      call secund(tstart)
cd     call elapsec(etstart)
cd     call date(datestr)
cd     write(98,*) MY_HOSTNM,' starts at ',datestr
c
c -- unpack the geometry,cpu,symmetry data
      call para_get_info
c
c -- get the remaining initialization data
      call para_recv_bcastpack(TxJobInit)
      call para_unpack_int(iforwhat,1)
      lpwave=(iforwhat.eq.11)                    ! FTC flag
c  ****************************************************************
c  *  This incredible piece of coding is needed because for some  *
c  *  totally inexplicable and monumentally perverse reason       *
c  *  the basis set pointer "ibas" has its value decremented      *
c  *  by 1 somewhere in the bowels of the integral code. This     *
c  *  pointer is communicated via the common block /infob/ and    *
c  *  the "new" (incorrect) value is subsequently picked up in    *
c  *  this routine causing mayhem.                                *
        ibss = ibas     !  save the correct pointer value         *
c  ****************************************************************
c
c -- get number of dummy atoms (if any)
      call getival('ndum',ndum)
c  At this point, na is the total number of atoms (real+dummy)
c  twoint95 further down cahnges na back to the number of real atoms
c  Why this was done is a mystery. In the original code slave_oneint
c  was called with the real atoms, and thus misses the effect of the dummies.
c  natom is used only in the plane wave code for mulipoles
      natom = na-ndum
c     natom=na
c      print *, na,ndum
c
c ***************************************************************
c   INITIALIZATION OF FOURIER-TRANSFORM COULOMB
c ***************************************************************
c
      IF(lpwave) THEN
c
c --   receive from master data for FTC code
c --   first get filename and array sizes
c
        call para_recv_bcastpack(TxMP2File)
        call para_unpack_string(filname1,256)
c make filename unique
        call rmblan(filname1,256,len1)
        write(ch3,'(".",I2.2)') MY_GID
        filname1 = filname1(1:len1)//ch3
        len1 = len1+3

        call para_unpack_int(ncspl,1)
        call para_unpack_int(ncsdiff,1)
        call para_unpack_int(ncfpl,1)
        call para_unpack_int(nomultipole,1)
c
        call setival('nomultipole',nomultipole)
        call getint(ncs,ilistsd)
        call setival('ftc0',ilistsd)
c
c       allocate integer arrays icsdiff(ncsdiff),sharpness(ncs),
c                icspltype(ncspl),icsplsize(ncspl),cssharps(ncspl)
c
        call mmark
        call getint(ncsdiff,iicsdiff)
        call getint(ncs,isharpness)
        call getint(ncspl,iicspltype)
        call getint(ncspl,iicsplsize)
        call getint(ncspl,icssharps)
        If(nomultipole.EQ.0) Then
c
c       allocate integer*1 array multipolecut(natom,natom)
c
          call getint_1(natom*natom,imultipolecut)
          call setival('imultipolecut',imultipolecut)
        EndIf
c
c
c --   now get initialization data
c
        call para_recv_bcastpack(TxFTCInit)
        call para_unpack_int(bl(iicsdiff),ncsdiff)
        call para_unpack_int(bl(isharpness),ncs)
        call para_unpack_int(bl(iicspltype),ncspl)
        call para_unpack_int(bl(iicsplsize),ncspl)
        call para_unpack_int(bl(icssharps),ncspl)
        call para_unpack_int(bl(ilistsd),ncs)
        If(nomultipole.EQ.0) Then
           call para_unpack_real(dist2mp,1)
c -- this array is actually integer*1, but it's OK - JB
           call para_unpack_byte(bl(imultipolecut),natom*natom)
        EndIf
      ENDIF
c----------------------------------------------------------------
c
c -- unpack SCF/DFT data
      call para_recv_bcastpack(TxScfInit)
      call para_unpack_int(nfock,1)
      call para_unpack_int(lsemi,1)
      call para_unpack_int(iforwhat,1)
      lpwave=(iforwhat.eq.11)                    ! FTC flag
      call para_unpack_int(ncachex,1)
      call para_unpack_int(NAlpha,1)
      call para_unpack_int(NBeta,1)
      call para_unpack_int(rhf,1)
      call para_unpack_real(threshx,1)
      call para_unpack_real(threshy,1)
c 
      iccdd=1
      if(lpwave) call para_unpack_int(iccdd,1)
      call setival('ccdd',iccdd)
      call setival('nomultipole',1)      ! disable multipoles - KW broke it
      call para_unpack_int(iroute,1)
      call para_unpack_int(istab,1)
c .......................................................
      call para_unpack_int(idft,1)
      If(idft.gt.0) Then
        call para_unpack_real(aXX,18)
        call para_unpack_int(nrad,1)
        call para_unpack_int(nang,1)
        call para_unpack_int(lrad,1)
        call para_unpack_int(lang,1)
        call para_unpack_int(IradQ,1)
        call para_unpack_int(NBatch,1)
        call setival('nmo ',NAlpha)
c
c -- restore coordinates, atomic charges and get symmetry data
        call getmem(3*natom,ixnc)               ! nuclear coordinates
        call getmem(natom,iqa)                  ! atomic charges/numbers
c
        call getnucdat(natom,bl(ixnc),bl(iqa))
        ax = aXX(1)    ! amount of "exact exchange"
        ityp=1
        If(rhf) ityp=0
        call setup_dft(idft,   ityp,   nfock,  nrad,   nang,
     $                 NBatch, bl(iqa),bl(ixnc))
cc        If(ityp.eq.1) call getmem(nscr,iscr)
        If(ityp.le.1) call getmem(nscr,iscr)
        call getint(natom,iunq)                 ! symmetry-unique atoms
        call getint(natom,islv)                 ! atom/slave array
      EndIf
c ........................................................
c
c get values necessary for initializing the integrals (twoint95)
c most of the things were not sent from the master
c see comments in ptwoint.f
c
      call getival('maxprice',maxpricex)
      iroutex=iroute
      call getival('chec',icheckx)
      call getival('iprn',iprintx)
      call getival('istat',istat)
c
c rest of the things needed for twoint95
c
      call getival('lcore',lcore)
      call getival('nsym',nsym)
c
cd      write(98,9991) MY_HOSTNM,MY_GID
 9991 format('SCF_SLAVE on the ',a60/
     *       'computat.group id=',i10/)
cd     call f_lush_(98)
c
c Initialize the local integral system
c
      call setival('iroute',iroute)
      call setival('stab',istab)
      call setival('iforwhat',iforwhat)
c
      call twoint95(lcore,nsym,maxpricex,ncachex,iroutex,
     *              icheckx,iprintx,threshx,istat,iforwhat,
     *              lsemi,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
cd     write(98,*)'The local integral system has been initialized'
cd     call elapsec(etend)
cd     write(98,'(a,f8.3,a)') 'startup time',etend-etstart,' sec'
cd     etstart=etend
cd     call f_lush_(98)
cd     call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
cd     write(98 ,1100) nreq,nmark,lastadr,mxmem,memtot
cd     call f_lush_(98)
c
c Reserve memory for labels, density and fock matrices
c
      call getmem(maxlabels,labels)
      isize=ncf*(ncf+1)/2
      ntri=isize
      ifsize=isize*nfock
      call getmem(ifsize,ifA)
      call getmem(ifsize,idA)
      call getmem(isize,iden)
c -- unrestricted open shell
      If(.NOT.rhf) Then
        call getmem(isize,ifB)
        call getmem(isize,idB)
      EndIf
c
      call getival('ictr',ictr)
c
c------------------------------------------------------
c  ONE-ELECTRON  INTEGRALS
c------------------------------------------------------
c
c  This was incorrect because na AT THIS POINT (after the call to TWOINT95
c  does not include the dummy atoms which are needed for point charges
c      call slave_oneint(1,     na,  bl(ifA),
c
c      print *, 'before slave_oneint',na,ndum
       call slave_oneint(1, na+ndum, bl(ifA),
     $                bl(ictr),0,      0,   bl(ibss), bl(inuc),
     $                  ncs)
c----------------------------------------------------------------
c  COSMO
c----------------------------------------------------------------
      call getival('cosmo',icosmo)
c
      If(icosmo.NE.0) Then
c
c  initialize COSMO for the slave, and compute surface
c  matrix elements
c
c communication and memory allocation can not be intermixed as
c memory is needed temporarily for MPI buffers              
c
c  allocate memory for cosmo (size known previously)
c
        call mmark
cc        call getmem(3*natom,icoxyz)                ! atomic coordinates
        call getmem(ntri,ivimat)                   ! scratch area for matrix elements
        call getmem(ntri,ih0cos)                   ! COSMO one electron hamiltonian
c        
        call para_recv_bcastpack(TxCosmoDat)              
cc        call para_unpack_int(nratom,1)
        call para_unpack_int(nps,1)
        call para_unpack_int(npspher,1)
        call para_unpack_string(conel,256)
        call para_unpack_real(fepsi,1)
        call para_unpack_real(cthr,1)
cc        call para_unpack_real(bl(icoxyz),3*natom)
cc        call setival('c_nratom',nratom)
        call setival('c_nps',nps)
        call setival('c_npsphe',npspher)
        call setchval('c_onel',conel)
        call setrval('c_fepsi',fepsi)
c
c   need to set this because it is used by the one electron
c   integral routine
c
        call setrval('ithr',cthr)
c
c  allocate memory for cosmo (sizes received now)
c
        call getmem(3*(nps+npspher),icosurf)       ! surface coordinates
        call getmem(max(nps,npspher),iphi)         ! surface potential
        call getmem(nps,iqcos)                     ! surface charges
        call getint(max(nps,npspher),imyseg) ! which segments belong to this slave
        
        call para_recv_bcastreal(bl(icosurf),3*(nps+npspher),TxCosmoDat)
c
c  surface matrix element
c
        call getival('ibas',ibasc)
        if(conel.ne.'direct')then
c
c  generate a unique name for the matrix elements file,
c
          call blankit(conel,256)
          call getchval('scrf',scrfile)
          lens=len_trim(scrfile)
          conel(1:lens)=scrfile(1:lens)
          conel(lens+1:lens+12)='.cosmo_onel.'
          lens=lens+12
          write (conel(lens+1:lens+8),fmt='(z8.8)') MY_GID
          call setchval('c_onel',conel)
c
c  compute and store matrix elements
c
          call slave_cosmo_surfrep(bl(ivimat),bl(icosurf),bl(ictr),
     $                 bl(ibasc),bl(imyseg),
     $                 nps,ncs,ncf,ntri)
        endif
      EndIf
c----------------------------------------------------------------
c------------------------------------------------------
c  TWO-ELECTRON  INTEGRALS
c------------------------------------------------------
      If(lpwave) call update_run_mode(bl,bl(ictr),ncf,lsemi,scftype)
c
c do the integral storage if necessary
c
      IF(scftype.ne.'full-direct') THEN
c
        If(lpwave) call symmoff
c
c open disk files for stored integrals if needed
c
        call getival('indisk',indisk)
        if(indisk.ne.0) then
           call getchval('scrf',scrfile)
           call open4store
        endif
c
c save integral threshold used to calculate stored integrals
c (in common /datstore/ thres_stored   )
c
         thres_stored=threshx
c
c zero out number of stored integrals and quartets
c
         integrals=0
         nqstore=0
         nblstore=0
         xntegrals=0
         integ62=0
         integ72=0
         integx2=0
         integx4=0
         integx8=0
c
c make denspar(ics,jcs)=1.0
c
         call getmem(ncs*ncs,idensp)
         call setup_denschw(ncs,bl(idensp))
c
 110     continue
c
c request work from the master, send number of integrals and quartets stored
c
         call para_initsend
         call para_pack_int(integrals,1)
         call para_pack_int(nqstore,1)
         call para_pack_int(nblstore,1)
c
         if(indisk.gt.0) then
c           rntegrals=xntegrals
c           rnteg62=integ62
c           rnteg72=integ72
c           rntegx2=integx2
c           rntegx4=integx4
c           rntegx8=integx8
            call para_pack_real(xntegrals,1)
            call para_pack_real(integ62,1)
            call para_pack_real(integ72,1)
            call para_pack_real(integx2,1)
            call para_pack_real(integx4,1)
            call para_pack_real(integx8,1)
c        write(98,*)'      xntegrals=',xntegrals
c        write(98,*)'      xntegra72=',integ72
c        write(98,*)'      xntegra62=',integ62
c        write(98,*)'      xntegrag2=',integx2
c        write(98,*)'      xntegrag4=',integx4
c        write(98,*)'      xntegrag8=',integx8
c           call f_lush(98)
         endif
c
         call para_send_pack(0,TxStorReq)
c get work(starting and finishing block, switch if all negatives were stored)
         call para_recv_pack(isource,TxStorAssign)
         call para_unpack_int(ifrom,1)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
ckw
         call para_unpack_int(ipack,1)
         call para_unpack_int(iqpac,1)
         call para_unpack_int(ibpac,1)
         call para_unpack_int(maxsi,1)
ckw
c
         call setival('intpack',ipack)
         call setival('qrtpack',iqpac)
         call setival('blkpack',ibpac)
         call setival('maxsint',maxsi)
c
c      write(98,*)'       : f-t=',ifrom,ito,' sec=',isecond
c      write(98,*)'       : packs=',ipack,iqpac,ibpac,maxsi
c      write(98,*)' before: nblst=',nblstore
c      call f_lush(98)
ckw
c store
         call int_store(ifrom,ito,bl,bl(ictr),bl(idensp),threshx,
     1                  bl(labels),isecond)
c
c stop storing if nothing was sent
         if(ito.ge.ifrom) goto 110
c
c find out if positive superblocks were stored and the last stored
c
         call para_send(MY_GID,0,TxStore)
         call para_recv_pack(isource,TxStore)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
c
c release memory allocated for idensp
c
         call retmem(1)
      ENDIF
c
c ***************************************************************
c   COMPUTATION OF CLASSICAL INTEGRALS
c ***************************************************************
c
c Set initial integral threshold
c
      thresh=threshy
 111  continue
      If(lpwave) call symmoff
      isto=1
c
c get densities
      call para_recv_bcastreal(bl(iden),isize,TxScfDens)
      call para_recv_bcastreal(bl(idA),ifsize,TxScfDens)
      If(.NOT.rhf) call para_recv_bcastreal(bl(idB),isize,TxScfDens)
      call para_recv_bcast(mgo,TxScfDens)
c
c zero the fock matrices
c
      call zeroit(bl(ifA),ifsize)
      If(.NOT.rhf) call zeroit(bl(ifB),isize)
c
c switch integral threshold if asked
c
      if(mgo.eq.2) thresh=threshx
cd     call elapsec(etend)
cd     write(98,'(a,f8.3,a)') 'dens dist time',etend-etstart,' sec'
cd     etstart=etend
c
      call mmark
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0.and.isto.eq.0) exit ! to make sure stored 
                                             ! integrals are used
c
        call int_fock(idft,ax,nblocks,nfock,rhf,ncf,bl(ictr),
     $              thresh,bl(iden),bl(ifA),bl(ifB),bl(idA),bl(idB),
     $              bl(labels),mywork,igran1)

        if(mywork.le.0)exit ! stored integrals have been used
c
cd     call elapsec(etend)
cd     write(98,9992) mywork,mywork+igran1-1,etend-etstart
cd     etstart=etend
 9992 format('blocks from',i10,' to:',i10,f8.3,' sec')
cd     call f_lush_(98)
      end do
      call retmark
c
c ***************************************************************
c   FOURIER-TRANSFORM COULOMB COMPUTATION STARTS HERE
c ***************************************************************
c
      IF(lpwave) THEN
c
c -- get initialization flag from master
c
        call para_recv_bcast(iinit,TxFTCInit0)
c
        IF(iinit.EQ.1) THEN
c
c ***************************************************************
c    SHARP GRIDS INITIALIZATION
c ***************************************************************
c
c -- get initial FTC grid data
c
c -- get npwx separately (to size up memory)
          call para_recv_bcast(npwx,TxFTCInit1)
c
c -- allocate memory before communication!!!
          If(.NOT.initialize) Then
c
c   allocate the integer arrays NSLAVE(ncspl),gridranges(6,ncs),
c        gridrange2(6,ncs),isharpgridrange1(6,ncspl),
c        isharpgridrange2(6,ncspl),nsharpovs(ncspl),
c        sharpovs(ncs,ncspl)
c   and the real array comprcoefs(5*ncspl*npwx)
c
            call getint(ncspl,inslave)
            call getint(6*ncs,igridranges)
            call getint(6*ncs,igridrange2)
            call getint(6*ncspl,iisharpgridrange1)
            call getint(6*ncspl,iisharpgridrange2)
            call getint(ncspl,insharpovs)
            call getint(ncs*ncspl,isharpovs)    ! use ncs to be sure
                                                ! to allocate enough
                                                ! memory
            call getmem(5*ncspl*npwx,icomprcoefs)
          EndIf
c
c   allocate array plbasdat(6,ncfpl)
c
          call mmark
          call getmem(6*ncfpl,iplbasdat)
c
          call para_recv_bcastpack(TxFTCInit1)
          call para_unpack_int(icorecut,1)
          call para_unpack_int(npwy,1)
          call para_unpack_int(npwz,1)
          call para_unpack_int(npwxe,1)
          call para_unpack_int(npwye,1)
          call para_unpack_int(npwze,1)
          npwyz=npwz*npwy
          npwxyz=npwx*npwyz
          call para_unpack_int(maxovs,1)
          call para_unpack_real(Lxmin,1)
          call para_unpack_real(Lymin,1)
          call para_unpack_real(Lzmin,1)
          call para_unpack_real(Lxe,1)
          call para_unpack_real(Lye,1)
          call para_unpack_real(Lze,1)
          call para_unpack_real(griddens,1)
          call para_unpack_real(PLDmax,1)
          call para_unpack_real(str,1)
          call para_unpack_real(c1r,1)
          call para_unpack_real(c2r,1)
          call para_unpack_real(c3r,1)
          call para_unpack_real(c1l,1)
          call para_unpack_real(c2l,1)
          call para_unpack_real(c3l,1)
          call para_unpack_real(c4l,1)
          call para_unpack_real(fskal0,1)
c
          call para_unpack_int(bl(igridranges),6*ncs)
          call para_unpack_int(bl(igridrange2),6*ncs)
          call para_unpack_int(bl(iisharpgridrange1),
     $                                                 6*ncspl)
          call para_unpack_int(bl(iisharpgridrange2),
     $                                                 6*ncspl)
          call para_unpack_int(bl(insharpovs),ncspl)
          call para_unpack_int(bl(isharpovs),maxovs*ncspl)
          call para_unpack_real(bl(iplbasdat),6*ncfpl)
c
          initialize = .True.
c
          call slave_sharpgrid(
     &   griddens,       Lxmin,
     &   Lymin,         Lzmin,          ncfpl,          isharpgrd,
     &   icorecut,      filname1,       len1,           bl(inslave),
     &   ntot,          bl(icssharps),  bl(iicspltype), bl(iicsplsize),
     &   bl(igridrange2),bl(iisharpgridrange2),bl(iplbasdat),
     &                                                  bl(icomprcoefs))
c
          call retmark   ! deallocate plbasdat
c
        ENDIF    ! end of sharp grids initialization
c
c ***************************************************************
c   BUILD UP THE SMOOTH DENSITY
c ***************************************************************
c
c -- scale Fock matrix
        fskal = thresh/fskal0
        call vscal(ntri,fskal,bl(ifA))
c
c
c     allocate arrays ro1(npwz,npwy) and ro2(npwz,npwy)
c
        call mmark
        call getmem(npwyz,iro1)
        call getmem(npwyz,iro2)
        call zeroit(bl(iro1),npwyz)
        call zeroit(bl(iro2),npwyz)
c
        call calc_exp_dims(
     &   ncsdiff,       bl(iicsdiff),  ncs,           bl(ictr),
     &   bl(ibss),     bl(igridranges),iystoredim,    izstoredim,
     &   ipexpdim)
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
     &   ncsdiff,       bl(iicsdiff),  ncs,           bl(ictr),
     &   bl(ibss),     bl(igridranges),iystoredim,    bl(iyexpstore),
     &   izstoredim,    bl(izexpstore),ipexpdim,      bl(iipyexpstore),
     & bl(iipzexpstore),Lymin,         Lzmin,         griddens)
c
        call slave_smoothd(
     &   ncf,           ncs,
     &   ncspl,         npwy,          npwz,          griddens,
     &   Lxmin,         Lymin,         Lzmin,         bl(igridranges),
     &   bl(iro1),      bl(iro2),      bl(isharpness),bl(ida),
     &   bl(ictr),     bl(ibss),      iystoredim,    izstoredim,
     &   ipexpdim,      bl(iyexpstore),bl(izexpstore),bl(iipyexpstore),
     &   bl(iipzexpstore))
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE SMOOTH DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      CALL Slave_COULOMB_on_FTGRID(
     &                 npwx,    npwy,    npwz,
     &                 npwxe,   npwye,   npwze,   PLDmax,  Lxe,
     &                 Lye,     Lze)
c
c ***************************************************************
c   CALCULATE MIXED (SMOOTH + SHARP) DENSITY
c ***************************************************************
c
        call slave_mixed_dpf(
     &   ncf,           ncs,
     &   ncspl,         npwx,          npwy,          npwz,
     &   ntot,          maxovs,        icorecut,      isharpgrd,
     &   filname1,      len1,          griddens,      Lxmin,
     &   Lymin,         Lzmin,         bl(inslave),   bl(iicspltype),
     &   bl(icssharps), bl(igridrange2),bl(ictr),    bl(ibss),
     &bl(iisharpgridrange1),bl(isharpness),bl(isharpovs),bl(insharpovs),
     &bl(iisharpgridrange2),bl(iro1),  bl(iro2),      bl(ida),
     &   bl(ifa),      bl(icomprcoefs),iystoredim,    izstoredim,
     &   ipexpdim,      bl(iyexpstore),bl(izexpstore),bl(iipyexpstore),
     &   bl(iipzexpstore))
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  FFTW COULOMB BUILDING WITH OUR NEW ANALYTICAL OPERATOR TO DECOUPLE
c  EXACTLY THE PERIODIC IMAGES (APPLYING IT FOR THE MIXED DENSITY)
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      CALL Slave_COULOMB_on_FTGRID(
     &                 npwx,    npwy,    npwz,
     &                 npwxe,   npwye,   npwze,   PLDmax,  Lxe,
     &                 Lye,     Lze)
c
c ***************************************************************
c   CALCULATE FOCK MATRIX FOR SMOOTH + MIXED DENSITY
c ***************************************************************
c
c   allocate array ro3(npwz,npwy)
c
        call getmem(npwyz,iro3)
c
        call slave_smooth_f(
     &   ncf,           ncs,           ncspl,
     &   npwx,          npwy,          npwz,          griddens,
     &   Lxmin,         Lymin,         Lzmin,         bl(ictr),
     &   bl(ibss),      bl(iro1),      bl(iro2),      bl(iro3),
     &   bl(ifa),      bl(igridranges),bl(isharpness),iystoredim,
     &   izstoredim,    ipexpdim,      bl(iyexpstore),bl(izexpstore),
     &   bl(iipyexpstore),bl(iipzexpstore))
c
        call retmark  ! deallocate ro1,ro2,ro3,yexpstore,zexpstore,
                      !  ipyexpstore,ipzexpstore
c
c ***************************************************************
c   MULTIPOLE CONTRIBUTIONS FOR <CD|CD> INTEGRALS
c ***************************************************************
c
        If(nomultipole.EQ.0) Then
c
          CALL Slave_Multipoles_Main(
     $         ncs,       bl(ictr),  bl(ibss), 
     $         ncspl, bl(icssharps), natom,     ncf,        bl(ida),
     $     bl(isharpovs), maxovs,bl(insharpovs),dist2mp,bl(isharpness),
     $         bl(ixnc),  bl(imultipolecut),    fskal0,     bl(ifA))
        EndIf

c ***************************************************************
c   FTC END
c ***************************************************************
c
      call symmon
c
      ENDIF
c
c summarize the results
c
      call para_reduce(bl(ifA),ifsize,TxScfFA)
      If(.not.rhf) call para_reduce(bl(ifB),isize,TxScfFB)
c
c .......................................................
c  COSMO
c .......................................................
      If(icosmo.NE.0) Then
c
c  compute the potential on the COSMO surface
c  bl(ih0cos) is used as storage for the density matrix
c
        call slave_cosmo_pot(bl(ih0cos),bl(ictr),bl(ibasc),bl(iphi),
     $             bl(icosurf),bl(ivimat),bl(imyseg),
     $             nps,ncf,ncs,ntri)
c
c  compute COSMO contribution to one electron hamiltonian
c
        call slave_cosmo_h0(bl(ih0cos),bl(ivimat),bl(ictr),bl(ibasc),
     $             bl(icosurf),bl(iqcos),bl(imyseg),fepsi,
     $             nps,ncf,ncs,ntri)
      EndIf
c .......................................................
      IF(idft.GT.0) THEN
cd     call elapsec(etend)
cd     write(98,'(a,f8.3,a)') '2e reduce',(etend-etstart),' sec'
cd     etstart=etend
c
c -- we are going to reuse the SCF memory for DFT
c       bl(ifA) will be used for the alpha exchange-correlation matrix
c       bl(idA) will be used for the alpha density matrix
c       bl(ifB) will be used for the beta exchange-correlation matrix
c       bl(idB) will be used for the beta density matrix
c       bl(iden) will be used for the maximum density element per column
c
        IEntry = 1                      ! set for first entry
        call zeroit(bl(ifA),ifsize)
        exc=zero
        el=zero
        If(.NOT.rhf) call zeroit(bl(ifB),isize)
c
        If(nsym.gt.0) Then
          call getival('nsyo',nsy)
          call getival('ngener',ngen)
          call getival('SymNuPr1',nupr)
        EndIf
c
c -- get the variable DFT data from the master
        call para_recv_bcastpack(TxDftDat)              
        call para_unpack_int(lgrid,1)
        call para_unpack_int(lsemi,1)
        call para_unpack_int(NQ,1)
        call para_unpack_int(bl(iunq),NQ)
        call para_unpack_int(bl(islv),NQ)
        call para_unpack_real(factor,1)
        call para_unpack_real(thrsh,1)
        call para_unpack_real(bl(idA),isize)
        If(.NOT.rhf) call para_unpack_real(bl(idB),isize)
        call para_unpack_real(bl(iden),ncf)
cd     call elapsec(etend)
cd     write(98,'(a,f8.3,a)') 'dft uptime',(etend-etstart),' sec'
cd     etstart=etend
c
        If(lsemi.GT.0) llsemi=1     ! set if semi-direct DFT ever used
c
c -- if requested delete all grid/potential files
        If(lgrid.EQ.2) Then
          Call TidyGRID(MY_GID,NQ,bl(iunq),bl(islv),lsemi)
        EndIf
c
 100    CONTINUE
c
c -- request atom we are currently doing from master
        call para_get_atom(ICntr)
c
        IF(ICntr.gt.0) THEN
          If(rhf) Then
          CALL DFTFockC(idft,   ICntr, natom,   bl(ixnc), bl(iqa),
     $                  nsym,   ngen,  bl(nsy), bl(nupr),  0,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                 bl(idst),bl(iaij),bl(ird),bl(iixx),bl(iiwt),
     $                 bl(ixgd),bl(iwt),thrsh,  ncf,    ncs,
     $                 bl(ibss),bl(ictr),bl(ipre),bl(iexp),bl(idA),
     $                 bl(iden), nscr, bl(iscr),bl(ifA), exc,
     $                  el,     lgrid,  lsemi,  IEntry)
          Else
          CALL DFTFockU(idft,   ICntr, natom,   bl(ixnc), bl(iqa),
     $                  nsym,   ngen,  bl(nsy), bl(nupr),  0,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                 bl(idst),bl(iaij),bl(ird),bl(iixx),bl(iiwt),
     $                 bl(ixgd),bl(iwt),thrsh,  ncf,    ncs,
     $                 bl(ibss),bl(ictr),bl(ipre),bl(iexp),bl(idA),
     $                  bl(idB),bl(iden), nscr, bl(iscr), bl(ifA),
     $                  bl(ifB),exc,    el,    lgrid,   lsemi,  IEntry)
          EndIf
cd     call elapsec(etend)
cd     write(98,9993) ICntr,etend-etstart
cd     etstart=etend
 9993 format('atom',i10,f8.3,' sec')
cd     call f_lush_(98)
          GO TO 100
        ELSE
c
c -- we are done
c -- send results back to master
cc          write(6,*) ' SLAVE:',MY_GID,' About to sum up'
cc          call f_lush(6)
          call para_reduce(bl(ifA),ifsize,TxDftF1)
          call para_reduce(exc,1,TxDftF2)
          call para_reduce(el,1,TxDftF3)
          If(.NOT.rhf) call para_reduce(bl(ifB),isize,TxDftF4)
cc          write(6,*) ' SLAVE:',MY_GID,' Finished summing up'
cc          call f_lush(6)
        ENDIF
cc
cd     call elapsec(etend)
cd     write(98,'(a,f8.3,a)') 'dft reduce',(etend-etstart),' sec'
cd     etstart=etend
      ENDIF
c .......................................................
c
c
c Any more SCF cycles?
c
      call para_recv_bcast(mgo,TxNext)
      If(mgo.eq.1.or.mgo.eq.2) goto 111
c
c ........................................................
c  COSMO
c ........................................................
      If(icosmo.NE.0) Then
c
c  Cosmo outlying charge correction
c
c  compute matrix elements for outer surface
c
        call slave_cosmo_surfrep(bl(ivimat),bl(icosurf+3*nps),
     $               bl(ictr),bl(ibasc),bl(imyseg),
     $               npspher,ncs,ncf,ntri)
c
c  compute the potential on the COSMO outer surface
c  bl(ih0cos) is used as storage for the density matrix
c
        call slave_cosmo_pot(bl(ih0cos),bl(ictr),bl(ibasc),bl(iphi),
     $             bl(icosurf+3*nps),bl(ivimat),bl(imyseg),
     $             npspher,ncf,ncs,ntri)
c
c -- delete cosmo matrix elements file
        call cosmo_del('c_onel  ')
c
c -- release cosmo memory
        call retmark
      EndIf
c ........................................................
c
c -- if DFT, delete all grid/potential files
      If(idft.GT.0) Then
        call TidyGRID(MY_GID,NQ,bl(iunq),bl(islv),llsemi)
c       call retint(2)
      EndIf
c ........................................................
      call getival('indisk',indisk)
      if(indisk.gt.0) then
         call clos4int()
      endif
c ........................................................
c
c send timing back to master then move on
c
c (I have added one more receiving and unpacking here, because
c  now the master executes one more call to para_next, after the
c  cosmo OC correction. This is because para_next expects also to
c  receive the timing data, so in the call at the end of SCF it
c  had to be fooled into not requesting the timing, to allow
c  for the parallel cosmo OC calculation MM 02/13/2004)
c
      call para_recv_bcast(mgo,TxNext)
      call secund(tend)
c
cd     print *,'MY_GID,cputime,hostname',MY_GID,tend-tstart,MY_HOSTNM
c
      call para_initsend
      call para_pack_real(tend-tstart,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
c
      if(lpwave)then
c
c  release FTC memory
c
        call retmark    !  deallocate cssharps,icspltype,icsplsize,
                        !  sharpness,icsdiff,nsharpovs,isharpgridrange2,
                        !  isharpgridrange1,gridrange2,gridranges,
                        !  NSLAVE sharpovs,comprcoefs
c
c -- delete the sharpgrid file
c
        filename = filname1(1:len1)//'.sharpgrd'
        OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),STATUS='OLD')
        CLOSE (UNIT=isharpgrd,STATUS='DELETE')
      endif
c
c -- release all memory
      call retmark
c
c   reset the values of some flags, in order not to screw up
c   subsequent job steps
c
      call setival('nomultipole',1)
cd     write(98,*)'SCF_INTEGRALS : END     '
cd     call f_lush_(98)
c
cd     call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
cd     write(98 ,1100) nreq,nmark,lastadr,mxmem,memtot
cd     call f_lush_(98)
c
 1100 format(' Memory status:'//
     *' request number=',i4,' memory marks=',i3,
     *' last used address=',i9,/,
     *' high water=',i9,' total available memory=',i9/)
c
      end
c=================================================
c
      subroutine para_get_info
c
c=================================================
c
c  receive initialization data from the Master. Updated to separate
c  the sending of single variables and of arrays, because in the MPI
c  version it is better not to allocate memory while packing and
c  unpacking data, as the memory array is used also as buffer for
c  send/receive
c

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      common /cpu/ intsize,iacc,icache,memreal
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /intlim/ limxmem,limblks,limpair
c
c  unpack the data and put in the depository
c  actually, most of the things are used from the commons above
c
c -- initialize variables
      call para_recv_bcastpack(TxInitVar)
c
c  output channel
c
      call para_unpack_int(iout,1)
      call setival('iout',iout)
c
c  geometry and basis set info
c
      call para_unpack_int(na,1)
      call para_unpack_int(ndum,1)
      call para_unpack_int(ncf,1)
      call para_unpack_int(ncs,1)
      call para_unpack_int(nbf,1)
      call para_unpack_int(nsh,1)
      call para_unpack_int(nonredun,1)
c
      call setival('na  ',na)
      call setival('ndum',ndum)
      call setival('ncf ',ncf)
      call setival('ncs ',ncs)
      call setival('nbf',nbf)
      call setival('nsh',nsh)
      call setival('nonredun',nonredun)
c
c  COSMO flag
c
      call para_unpack_int(icosmo,1)
      call setival('cosmo',icosmo)
c
c  CPU related information
c
      call para_unpack_int(intsize,1)
      call para_unpack_int(iacc,1)
      call para_unpack_int(icache,1)
      call para_unpack_int(memreal,1)
c
      call setival('ints',intsize)
      call setival('accu',iacc)
      call setival('cach',icache)
      call setival('memr',memreal)
c
c  the integral slaves use these from common /intlim/
c
      call para_unpack_int(limxmem,1)
      call para_unpack_int(limblks,1)
      call para_unpack_int(limpair,1)
c
c  Symmetry flag
c
      call para_unpack_int(nsym,1)
      call setival('nsym',nsym)
c
c  Reserve memory for the geometry and basis set matrices
c
      call getmem(13*nsh,ibas)
      call getmem(5*na,inuc)
      call getint(12*ncs,ictr)
c
c  receive initialization arrays
c
      call para_recv_bcastpack(TxInitArr)
      call para_unpack_real(bl(ibas),13*nsh)
      call para_unpack_real(bl(inuc),5*na)
      call para_unpack_int(bl(ictr),12*ncs)
      call setival('ibas',ibas)
      call setival('inuc',inuc)
      call setival('ictr',ictr)
c
c  Symmetry data
c
      if (nsym.gt.0) then
        call getint(7,iadr)
        call getint(na*nsym,nupair)
        nupair1 = nupair
        call getint(7*ncf,ifp)
        call getint(7*ncs,ifp1)
        if (ndum.gt.0) then
          natom = na-ndum
          call getint(natom*nsym,nupair1)
        endif
        call para_recv_bcastpack(TxSymData)
        call para_unpack_int(ngener,1)
        call para_unpack_int(bl(iadr),7)
        call para_unpack_int(bl(nupair),na*nsym)
        call para_unpack_int(bl(ifp),7*ncf)
        call para_unpack_int(bl(ifp1),7*ncs)
        if (ndum.gt.0) then
          call para_unpack_int(bl(nupair1),natom*nsym)
        endif
        call setival('ngener',ngener)
        call setival('nsyo',iadr)
        call setival('SymNuPr',nupair)
        call setival('SymNuPr1',nupair1)
        call setival('SymFunPr',ifp)
        call setival('SymShPr',ifp1)
      endif
      end
c=================================================
c
      subroutine do_post
c
c=================================================

      use memory
      use newpara
      
      implicit real*8 (a-h,o-z)
      character scftype*11, ch3*3
      Logical rhf,oneel
c     common /big/ bl(300)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /nmrint/ ichf,ncache,iprint,nogiao,noacce,ngauge,ntrans
c     common /intbl/ maxsh,inx(100)
      common /datstore/ thres_stored,isto,isecond,ito
      common /lindvec/ lind,idensp
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
      integer*8 xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
      common /forcdbl/ thre1,thre2,tchf
      dimension xintxx(9)
c
c  FTC stuff
c
      logical lpwave
      character*256 filename,filname1,scrfile
      real*8 Lxmin,Lymin,Lzmin
      parameter (isharpgrd=51)            ! unit number for FTC I/O
      common /corecut/ str,c1r,c2r,c3r,c1l,c2l,c3l,c4l !P6 polinom coefs
c
c ..............................................................
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     $                ixgd,iwt,iscr,icofm,nscr
c
      COMMON /DFTCoeff/ aXX(18) ! contains 18 density functional coefficients
c
      parameter (Zero=0.0d0)
c ..............................................................
c
c -- set memory mark
      call mmark
c
c get the corresponding host name
      call secund(tstart)
c
c -- unpack the geometry,cpu,symmetry data
      call para_get_info
c
c -- get the rest of the initialization data
      call para_recv_bcastpack(TxJobInit)
      call para_unpack_int(iforwhat,1)
      lpwave=(iforwhat.eq.13)                    ! FTC flag
c
c ***************************************************************
c   INITIALIZATION OF FOURIER-TRANSFORM COULOMB GRADIENT
c ***************************************************************
c
        IF(lpwave) THEN
c
c --   receive from master data for FTC code
c --   first get unique filename
          call para_recv_bcastpack(TxMP2File)
          call para_unpack_string(filname1,256)
c make filename unique
          call rmblan(filname1,256,len1)
          write(ch3,'(".",I2.2)') MY_GID
          filname1 = filname1(1:len1)//ch3
          len1 = len1+3
ckw
          call para_unpack_int(ncs,1)
ckw
          call para_unpack_int(ncspl,1)
          call para_unpack_int(ncsdiff,1)
          call para_unpack_int(ncfpl,1)
          call para_unpack_int(npwx,1)
          call para_unpack_int(npwy,1)
          call para_unpack_int(npwz,1)
          npwyz=npwy*npwz
          npwxyz=npwx*npwyz
c
c --   now get initialization data
c
          call mmark
          call getint(ncs,ilistsd)
          call setival('ftc0',ilistsd)
c
          call getint(ncsdiff,iicsdiff)
          call getint(ncs,isharpness)
          call getint(ncspl,iicspltype)
          call getint(ncspl,iicsplsize)
          call getint(ncspl,icssharps)
c
c   allocate the integer arrays NSLAVE(ncspl),gridranges(6,ncs),
c        gridrange2(6,ncs),isharpgridrange1(6,ncspl),
c        isharpgridrange2(6,ncspl),nsharpovs(ncspl),
c        sharpovs(ncs,ncspl)
c   and the real array comprcoefs(5*ncspl*npwx)
c
          call getint(ncspl,inslave)
          call getint(6*ncs,igridranges)
          call getint(6*ncs,igridrange2)
          call getint(6*ncspl,iisharpgridrange1)
          call getint(6*ncspl,iisharpgridrange2)
          call getint(ncspl,insharpovs)
          call getint(ncs*ncspl,isharpovs)    ! use ncs to be sure
                                              ! to allocate enough
                                              ! memory
          call getmem(5*ncspl*npwx,icomprcoefs)
          call getmem(6*ncfpl,iplbasdat)
c
c
          call para_recv_bcastpack(TxFTCInit)
          call para_unpack_int(bl(iicsdiff),ncsdiff)
          call para_unpack_int(bl(isharpness),ncs)
          call para_unpack_int(bl(iicspltype),ncspl)
          call para_unpack_int(bl(iicsplsize),ncspl)
          call para_unpack_int(bl(icssharps),ncspl)
          call para_unpack_int(bl(ilistsd),ncs)
c
          call para_unpack_int(icorecut,1)
          call para_unpack_int(maxovs,1)
          call para_unpack_real(Lxmin,1)
          call para_unpack_real(Lymin,1)
          call para_unpack_real(Lzmin,1)
          call para_unpack_real(griddens,1)
          call para_unpack_real(PLDmax,1)
          call para_unpack_real(str,1)
          call para_unpack_real(c1r,1)
          call para_unpack_real(c2r,1)
          call para_unpack_real(c3r,1)
          call para_unpack_real(c1l,1)
          call para_unpack_real(c2l,1)
          call para_unpack_real(c3l,1)
          call para_unpack_real(c4l,1)
          call para_unpack_int(bl(igridranges),6*ncs)
          call para_unpack_int(bl(igridrange2),6*ncs)
          call para_unpack_int(bl(iisharpgridrange1),
     $                                                   6*ncspl)
          call para_unpack_int(bl(iisharpgridrange2),
     $                                                   6*ncspl)
          call para_unpack_int(bl(insharpovs),ncspl)
          call para_unpack_int(bl(isharpovs),maxovs*ncspl)
          call para_unpack_real(bl(iplbasdat),6*ncfpl)
c
          call slave_sharpgrid(
     &   griddens,       Lxmin,
     &   Lymin,         Lzmin,          ncfpl,          isharpgrd,
     &   icorecut,      filname1,       len1,           bl(inslave),
     &   ntot,          bl(icssharps),  bl(iicspltype), bl(iicsplsize),
     &   bl(igridrange2),bl(iisharpgridrange2),bl(iplbasdat),
     &                                                  bl(icomprcoefs))
c
        ENDIF    ! end of sharp grids initialization
c----------------------------------------------------------------
c
c get the initial post-SCF data
c
      call para_recv_bcastpack(TxPostInit)
      call para_unpack_int(ncache,1)
      call para_unpack_int(iforwhat,1)
      call para_unpack_int(NAlpha,1)
      call para_unpack_int(NBeta,1)
      call para_unpack_int(rhf,1)
      call para_unpack_real(thre1,1)
      call para_unpack_real(thre2,1)
      call para_unpack_int(iroute,1)
      call para_unpack_int(istab,1)
c2006
      iccdd=1
      if(lpwave) call para_unpack_int(iccdd,1)
      call setival('ccdd',iccdd)
c2006
c
c get the data needed for NMR
c
      if (iforwhat.eq.2) then
         call para_unpack_real(thre3,1)
         call para_unpack_int(ichf,1)
         call para_unpack_int(nogiao,1)
         call para_unpack_int(noacce,1)
         call para_unpack_int(ngauge,1)
         call para_unpack_int(ntrans,1)
      end if
c
c setup the common /ganz/ needed later
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('lcore',lcore)
      call getival('nsym',nsym)
c
c get together the arguments needed for post_scf
c
      call getival('iprn',iprint)
      call getival('ictr',ictr)
      thresh=thre2
c
cd     write(98,9991) MY_HOSTNM,MY_GID
 9991 format('NMR_SLAVE on the ',a60/
     *       'master process id=',i10/
     *       'slave  process id=',i10/
     *       'computat.group id=',i10/)
cd     call f_lush_(98)
c
c prepare the integral system
c
      call setival('iroute',iroute)
      call setival('stab',istab)
c
c
      call post_scf(bl,bl(ictr),ncache,iprint,thresh,
     *              iforwhat,1,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
cd     write(98,*)'The local integral system has been initialized'
cd     call f_lush_(98)
c
c Reserve memory for the labels, the density and the first order
c Fock matrices. Get triple memory for density in NMR (density
c derivatives are put there later)
c
      ntri=ncf*(ncf+1)/2
c
      IF (iforwhat.eq.3.OR.iforwhat.eq.13) THEN
c -- gradient
        ntri3=3*na
        call getmem(ntri,lden)
        If(.not.rhf) call getmem(ntri,mden)
      ELSE IF (iforwhat.eq.2) THEN
c -- NMR
        ntri3=ntri*3
        call getmem(ntri3,lden)
c memory for screening density
        ncs2=2*ncs*ncs
        call getmem(ncs2,idscree)
      ENDIF
      call getmem(ntri3,lfxyz)
      call getmem(maxlabels,labels)
c
c zero out, get density
c
      call zeroit(bl(lfxyz),ntri3)
      call para_recv_bcastpack(TxPostDens)
      call para_unpack_real(bl(lden),ntri)
      if(iforwhat.eq.2) call para_unpack_real(bl(idscree),ncs2)
      If(.not.rhf) call para_unpack_real(bl(mden),ntri)
      call para_unpack_int(idft,1)
      call para_unpack_real(aXX,18)
c
ccccccccccccccccccccc
cc      write(6,*) ' About to start calculation'
cc      call prntbas(ncs,nsh,ncf,nbf,inx(ictr),bl(ibas),bl(lden))
cc      call f_lush(6)
cccccccccccccccccccccc
      If(iforwhat.eq.2) GO TO 99      ! NMR
c
c *************************************************
c  one-electron integrals
c *************************************************
c
c -- first we are going to do the electron-nuclear attraction
c    1-e integrals that were skipped during the serial evaluation
c    These will be parallelized over atoms
c
      If(.not.rhf) Then
c -- add alpha and beta densities together
        call getmem(ntri,llden)
        Call AddVEC(ntri,bl(lden),bl(mden),bl(llden))
      Else
        llden = lden
      EndIf
c
      oneel = .False.
c
 100  CONTINUE
c
c -- request atom we are currently doing from master
      call para_get_atom(IAtm)
c
      If(IAtm.GT.0) Then
c -- Calculate electron-nuclear attraction contributions to the forces
        call intof7(IAtm, bl(ictr), bl(ibas), bl(inuc), ncs,
     $              ncf,  ntri,     bl(llden), bl(lfxyz))
        oneel = .True.
        GO TO 100
      EndIf
c
c -- scale 1-el contribution (as rescaled later in 2-el part)
      If(oneel) CALL VScal(3*na,1.0d0/thresh,bl(lfxyz))
c
      If(.not.rhf) call retmem(1)      ! deallocate llden
c
  99  CONTINUE
c
c *************************************************
c  two-electron integrals
c *************************************************
c
c calculate integrals
c
      ax = aXX(1)    ! amount of "exact exchange"
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0) exit
c
cd     write(98,9992) mywork,mywork+igran-1
 9992 format(' giao_integral blocks from',i10,' to:',i10,' started')
cd     call f_lush_(98)
c
      if (iforwhat.eq.3.or.iforwhat.eq.13) then
        call int_grad(idft,ax,rhf,nblocks,bl,bl(ictr),ntri,thresh,
     *          bl(lden),bl(mden),bl(lfxyz),bl(labels),mywork,igran1)
      elseif (iforwhat.eq.2) then
        call int_giao(idft,ax,nblocks,bl,bl(ictr),ntri,thresh,
     *                bl(idscree),bl(lden),bl(lfxyz),
     *                bl(labels),mywork,igran1)
      end if
      end do
c
c do a global sum of the results
c
      call para_reduce(bl(lfxyz),ntri3,TxPostFock)
c
c ***************************************************************
c   FOURIER-TRANSFORM COULOMB GRADIENT
c ***************************************************************
c
      IF(lpwave) THEN
c
        call mmark
        call getmem(na,igradx)
        call ZeroIT(bl(igradx),na)
        call getmem(na,igrady)
        call ZeroIT(bl(igrady),na)
        call getmem(na,igradz)
        call ZeroIT(bl(igradz),na)
c
c -- expand density into full matrix
        call getmem(ncf*ncf,idens)
        call expand(ncf,bl(lden),bl(idens))
c
c ***************************************************************
c   GRADIENT FROM THE SMOOTH DENSITY
c ***************************************************************
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
        dx = 1.0d0/griddens
        npwregx = 1
        npwregy = npwy
        iyregmin = 1
        iyregmax = npwy
c
        call mmark
        call getmem(3*npwyz,iro1)
        call ZeroIT(bl(iro1),3*npwyz)
        call getmem(3*npwyz,iro2)
        call ZeroIT(bl(iro2),3*npwyz)
        call getmem(npwyz,iro3)
        call ZeroIT(bl(iro3),npwyz)
c
c -- get block of data from master
        call para_recv_pack(islave,TxFTCAssign)
        call para_unpack_int(ixstart,1)
        call para_unpack_int(ixend,1)
c
c -- now start main computation
c
 130    CONTINUE
c
c -- allocate array ro5(npwz,npwy,ixstart:ixend)
        call getmem(npwyz*(ixend-ixstart+1),iro5)
c
        call para_recv_real(bl(iro5),npwyz*(ixend-ixstart+1),
     &                      0,TxFTCAssign)
c
        ishift = 0
        do ix=ixstart,ixend
        xmin = Lxmin+float(ix-1)*dx
c
        call calc_smooth_gradients(
     $        ncf,         bl(iro1),       npwregx,        npwregy,
     $        npwz,     bl(igridranges),   bl(iro2),       ncs,
     $     bl(isharpness), ncspl,          bl(iro3),       bl(idens),
     $        bl(ictr),   bl(ibas),       xmin,           Lymin,
     $        Lzmin,    bl(iro5+ishift),   ix,             iyregmin,
     $        iyregmax,    griddens,       iystoredim, bl(iyexpstore),
     $     bl(iy2expstore),bl(iy3expstore),izstoredim, bl(izexpstore),
     $     bl(iz2expstore),bl(iz3expstore),ipexpdim,  bl(iipyexpstore),
     $     bl(iipzexpstore), na,           bl(igradx),     bl(igrady),
     $        bl(igradz))
c
        ishift = ishift+npwyz
        enddo
c
        call retmem(1)     !  deallocate ro5
c
c -- finished loop; get new data from master
        call para_send(MY_GID,0,TxDftReq)
        call para_recv_pack(islave,TxFTCAssign)
        call para_unpack_int(ixstart,1)
        call para_unpack_int(ixend,1)
c -- test for termination
        If(ixstart.NE.0.AND.ixend.NE.0) GO TO 130
c
        call retmark        ! deallocate ro1, ro2, ro3
c
c ***************************************************************
c   GRADIENT FROM THE MIXED DENSITY
c ***************************************************************
c
c  We know which sharp shells are stored on this slave
c  ntot - number; bl(inslave) - values
c  as obtained from routine <slave_sharpgrid>
c
        call getmem(npwxyz,iro4)
c
c -- first get the full smooth Coulomb potential
        call para_recv_bcastreal(bl(iro4),npwxyz,TxFTCDen)
c
c -- first part of mixed FTC gradient
        call slave_sharp1_grad(
     $       ncs,           ncspl,         icssharps, igridrange2,
     $     bl(iicsplsize),  iicspltype,    Lxmin,     Lymin,
     $       Lzmin,         griddens,      iplbasdat, iisharpgridrange2,
     $       iro4,          isharpness,    idens,     iisharpgridrange1,
     $       maxovs,        insharpovs,    isharpovs, npwx,
     $       npwy,          npwz,          icorecut,  bl(ictr),
     $       iystoredim,    iyexpstore,  iy2expstore, iy3expstore,
     $       izstoredim,    izexpstore,  iz2expstore, iz3expstore,
     $       ipexpdim,    iipyexpstore, iipzexpstore, ibas,
     $       ncf,           ncfpl,         na,        ntot,
     $       bl(inslave),   bl(igradx),    bl(igrady),bl(igradz))
c
c -- second part of mixed FTC gradient
c    open file
        filename = filname1(1:len1)//'.sharpgrd'
        OPEN(UNIT=isharpgrd,FILE=filename(1:len1+9),
     &       FORM='UNFORMATTED',STATUS='OLD')
c
        call slave_sharp2_grad(
     $       ncs,           ncspl,         icssharps, igridrange2,
     $       icomprcoefs,   iicspltype,    Lxmin,     Lymin,
     $       Lzmin,         griddens,      isharpgrd, iisharpgridrange2,
     $       iro4,          isharpness,    idens,     iisharpgridrange1,
     $       maxovs,        insharpovs,    isharpovs, npwx,
     $       npwy,          npwz,          icorecut,  bl(ictr),
     $       iystoredim,    iyexpstore,  iy2expstore, iy3expstore,
     $       izstoredim,    izexpstore,  iz2expstore, iz3expstore,
     $       ipexpdim,    iipyexpstore, iipzexpstore, ibas,
     $       ncf,           na,            ntot,      bl(inslave),
     $       bl(igradx),    bl(igrady),    bl(igradz))
c
        CLOSE(UNIT=isharpgrd)
        call retmem(1)      ! deallocate ro4
c
c -- accumulate total FTC gradient on master
        call para_reduce(bl(igradx),na,TxFTCGX)
        call para_reduce(bl(igrady),na,TxFTCGY)
        call para_reduce(bl(igradz),na,TxFTCGZ)
c
        call retmark      ! deallocate main FTC Gradient storage
        call retmark      ! deallocate FTC initialization storage
c
c -- reset the values of some FTC flags in order not to screw up
c -- subsequent job steps
       call setival('ftc0',0)     ! switch off FTC
       call setival('nomultipole',1)           ! switch off multipoles
c ***************************************************************
c   FTC END
c ***************************************************************
c
      ENDIF
c
c we are ready with coulomb forces
c
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c   PARALLEL DFT SECTION
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
      IF(idft.GT.0.AND.(iforwhat.eq.3.OR.iforwhat.eq.13)) THEN
c
COMMENT
c------
c  There are problems with the basis set data which is corrupted
c  due to expansion of SP shells into separate S and P for coulomb
c  forces - this needs to be looked at.  At the moment this data
c  is sent over from the Master.
c ----------------------------------------------------------------
c
        call mmark
c -- get initial DFT data
        IEntry = 1                   ! set for first entry
        call getival('ibas',ibas)    ! god knows what value this has
ckw
      IF(lpwave) THEN
        call getint(12*ncs,ictr)
        call getmem(13*nsh,ibas)
        call setival('ictr',ictr)
        call setival('ibas',ibas)
      ENDIF
ckw
        call para_recv_bcastpack(TxDftInit)              
        call para_unpack_int(nrad,1)
        call para_unpack_int(nang,1)
        call para_unpack_int(lrad,1)
        call para_unpack_int(lang,1)
        call para_unpack_int(IradQ,1)
        call para_unpack_int(NBatch,1)
        call para_unpack_int(IdWt,1)
c .........................................................
c -- these next lines needed because the BASDAT
c -- array is corrupted in the slaves
        call para_unpack_int(ncs,1)
        call para_unpack_int(nsh,1)
        call para_unpack_int(bl(ictr),12*ncs)
        call para_unpack_real(bl(ibas),13*nsh)
c .........................................................
        call para_unpack_real(factor,1)
        call para_unpack_real(thrsh,1)
c
c -- ensure correct values are in depository
        call setival('nmo ',NAlpha)
        call setival('ncs ',ncs)
        call setival('nsh ',nsh)
c
c -- get number of dummy atoms (if any)
        call getival('ndum',ndum)
        natom = na-ndum
c
c -- restore coordinates and atomic charges
c       call immark
        call getmem(3*natom,ixnc)               !  nuclear coordinates
        call getmem(natom,iqa)                  !  atomic charges/numbers
        call getmem(3*natom,lfxyz)              !  reallocate as deleted
c
        call getnucdat(natom,bl(ixnc),bl(iqa))
c
        ityp = 3
        If(rhf) ityp=2
        call setup_dft(idft,   ityp,   2,      nrad,   nang,
     $                 NBatch, bl(iqa),bl(ixnc))
c
c -- make sure we have all pointers
        If(nsym.gt.0) Then
          call getival('nsyo',isy)
          call getival('ngener',ngen)
          call getival('SymNuPr1',nupr)
        EndIf
c
c -- get potentially quite obscene amounts of memory
c        
c needed to be zeroed out because it could happen in parallel that a
c node does not get an atom - and the code mysteriously and irreducibly
c died. This was happening only for very small systems. (gm)
c        
        ncf2 = ncf**2
        call getint(natom,iatm)              !  no. basis functions per atom
        call getint(7*ncf,ifp)               !  list of basis function pairs
        call getmem(ncf2,igx)                !  partial DFT x derivatives
        call zeroit(bl(igx),ncf2)
        call getmem(ncf2,igy)                !  partial DFT y derivatives
        call zeroit(bl(igy),ncf2)
        call getmem(ncf2,igz)                !  partial DFT z derivatives
        call zeroit(bl(igz),ncf2)
        call getmem(natom*3,igwt)            !  weight derivatives
        call zeroit(bl(igwt),natom*3)
c
        If(.NOT.rhf) Then
c -- yet more memory for open shell
          call getmem(ncf2,igxb)             !  partial beta DFT x derivatives
          call zeroit(bl(igxb),ncf2)
          call getmem(ncf2,igyb)             !  partial beta DFT y derivatives
          call zeroit(bl(igyb),ncf2)
          call getmem(ncf2,igzb)             !  partial beta DFT z derivatives
          call zeroit(bl(igzb),ncf2)
        EndIf

c
c -- unpack second batch of DFT data
ckw
        IF(lpwave) then
           call getmem(ntri,lden)
           If(.NOT.rhf) Then 
             call getmem(ntri,mden)
           endif
        EndIf
ckw
        call para_recv_bcastpack(TxDftDat)              
        call para_unpack_int(bl(iatm),natom)
        call para_unpack_int(bl(ifp),7*ncf)
        call para_unpack_real(bl(lden),ntri)
        If(.NOT.rhf) Then
          call para_unpack_real(bl(mden),ntri)
        EndIf
c
c  zero out dft quantities, because we could actually do no
c  work at all, if there are more slaves that unique atoms,
c  thus we must be sure not to send uninitialized data
c  back to the master during the reduce process
c
        call zeroit(bl(igx),ncf2)
        call zeroit(bl(igy),ncf2)
        call zeroit(bl(igz),ncf2)
        If(.not.rhf) Then
          call zeroit(bl(igxb),ncf2)
          call zeroit(bl(igyb),ncf2)
          call zeroit(bl(igzb),ncf2)
        EndIf
        el=zero
c
 200    CONTINUE
c
c -- request atom we are currently doing from master
        call para_get_atom(ICntr)
c
        IF(ICntr.gt.0) THEN
          If(rhf) Then
          CALL DFTGradC(idft,   ICntr, natom,   bl(ixnc), bl(iqa),
     $                  nsym,   ngen, bl(isy),bl(nupr),  0,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                 bl(idst),bl(iaij),bl(ird),bl(iixx),bl(iiwt),
     $                 bl(ixgd),bl(iwt),thrsh,  ncf,    ncs,
     $                 bl(ibas),bl(ictr),bl(ipre),bl(iexp), IdWt,
     $                  bl(iatm), bl(lden),nscr,  bl(iscr),bl(igx),
     $                  bl(igy),bl(igz),bl(igwt),bl(lfxyz), el,
     $                  IEntry)
          Else
          CALL DFTGradU(idft,   ICntr, natom,   bl(ixnc), bl(iqa),
     $                  nsym,   ngen, bl(isy),bl(nupr),  0,
     $                  lrad,   lang,   IradQ,  factor, NBatch,
     $                 bl(idst),bl(iaij),bl(ird),bl(iixx),bl(iiwt),
     $                 bl(ixgd),bl(iwt),thrsh,  ncf,    ncs,
     $                 bl(ibas),bl(ictr),bl(ipre),bl(iexp), IdWt,
     $                  bl(iatm), bl(lden),bl(mden),nscr,bl(iscr),
     $                 bl(igx), bl(igy), bl(igz),bl(igxb),bl(igyb),
     $                 bl(igzb),bl(igwt),bl(lfxyz), el, IEntry)
          EndIf
          GO TO 200
        ELSE
c
c -- we are done
c -- send results back to master
          call para_reduce(bl(igx),ncf2,TxDftF1)
          call para_reduce(bl(igy),ncf2,TxDftF2)
          call para_reduce(bl(igz),ncf2,TxDftF3)
          call para_reduce(el,1,TxDftF4)
c
          If(.not.rhf) Then
            call para_reduce(bl(igxb),ncf2,TxDftF1b)
            call para_reduce(bl(igyb),ncf2,TxDftF2b)
            call para_reduce(bl(igzb),ncf2,TxDftF3b)
          EndIf

c
c -- remove DFT memory
c         call retimark
          call retmark
        ENDIF
c
c DFT NMR
c
      ELSE IF(idft.GT.0.AND.iforwhat.eq.2) THEN
c
c -- we are going to reuse memory for DFT NMR
c       bl(iden) will be used for the occupied MOs
c
        call zeroit(bl(lfxyz),ntri3)
        IEntry = 1                      ! set for first entry
        call getival('ibas',ibss)       ! god knows what value this has
c
c -- get number of dummy atoms (if any)
        call getival('ndum',ndum)
        natom = na-ndum
c
        If(nsym.gt.0) Then
          call getival('nsyo',isy)
          call getival('ngener',ngen)
          call getival('SymNuPr1',nupr)
        EndIf
c allocate memory for malkin        
        if(Malk.ne.0) then
            call getmem(ncf*(ncf-NAlpha),ivmo)
            call getmem(NAlpha*(ncf-NAlpha),imalk)
            call zeroit(bl(imalk),NAlpha*(ncf-NAlpha))
        else
            imalk=lben     ! just to have some value,
            ivmo=lden      ! these are not used w/o malkin
        endif
c
c -- get the DFT data from the master
        call para_recv_bcastpack(TxDftNDat)              
        call para_unpack_int(nrad,1)
        call para_unpack_int(nang,1)
        call para_unpack_int(lrad,1)
        call para_unpack_int(lang,1)
        call para_unpack_int(IradQ,1)
        call para_unpack_int(NBatch,1)
        call para_unpack_int(ncphf,1)
        call para_unpack_int(Malk,1)
        call para_unpack_real(factor,1)
        call para_unpack_real(thrsh,1)
        call para_unpack_real(bl(lden),ncf*NAlpha)
        if(Malk.ne.0) then
            call para_unpack_real(bl(ivmo),ncf*(ncf-NAlpha))
        endif
        call setival('nmo ',NAlpha)
c
c -- restore coordinates, atomic charges and get symmetry data
        call getmem(3*natom,ixnc)               ! nuclear coordinates
        call getmem(natom,iqa)                  ! atomic charges/numbers
c
        call getnucdat(natom,bl(ixnc),bl(iqa))
        call setup_dft(idft,   4, 1,  nrad,   nang,
     $                 NBatch, bl(iqa),bl(ixnc))
        call getmem(nscr,iscr)
c
c  zero out dft quantities, because we could actually do no
c  work at all, if there are more slaves that unique atoms,
c  thus we must be sure not to send uninitialized data
c  back to the master during the reduce process
c
        call zeroit(bl(lfxyz),ntri3)
        if(Malk.ne.0) then
          call zeroit(bl(imalk),NAlpha*(ncf-NAlpha))
        endif
c
 300    CONTINUE
c
c -- request atom we are currently doing from master
        call para_send(MY_GID,0,TxDftNReq)
        call para_recv(ICntr,isource,TxDftNAssign)
c
        IF(ICntr.gt.0) THEN
           CALL DFTNMRC(idft,   ICntr, natom,   bl(ixnc),bl(iqa),
     *                  nsym,   ngen,  bl(isy), bl(nupr), 0,
     *                  lrad,   lang,   IradQ,  factor, NBatch,
     *                  bl(idst),bl(iaij),bl(ird),bl(iixx),bl(iiwt),
     *                  bl(ixgd),bl(iwt),thrsh, ncf,ncs,
     *                  bl(ibss),bl(ictr),bl(ipre),bl(iexp),NAlpha,
     *                  bl(lden),bl(ivmo),Malk,nscr,bl(iscr),
     *                  bl(lfxyz),bl(imalk),IEntry)
          GO TO 300
        ELSE
c
c -- we are done
c -- send results back to master
cd     write(98,*)'DFTNMR: before reduce '
cd     call f_lush_(98)
          call para_reduce(bl(lfxyz),ntri3,TxDftN)
        if(Malk.ne.0) then
          call para_reduce(bl(imalk),NAlpha*(ncf-NAlpha),TxDftN2)
        endif
        ENDIF
cc
      ENDIF
c :::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c no cphf for gradients, nonhybrid dft, nonhybrid dft with nocphf option
      if (iforwhat.eq.3.or.(idft.gt.0.and.(ax.EQ.Zero.or.ncphf.eq.0)))
     *   goto 115
c----------------------------------------------------------------
c NOW IT'S TIME FOR COUPLED PERTURBED HARTREE-FOCK!!!!!!
c thresh=thre2  ! final integral threshold
c thresh=thre3  ! loose integral threshold
c
cd     write(98,*)'CPHF_INTEGRALS : BEGINNING '
cd     call f_lush_(98)
c
c do the integral storage if necessary
c
      if(scftype.ne.'full-direct') then
c
c open disk files for stored integrals if needed
c
        call getival('indisk',indisk)
        if(indisk.ne.0) then
           call getchval('scrf',scrfile)
           call open4store
        endif
c
c save integral threshold used to calculate stored integrals
c (in common /datstore/ thres_stored   )
c
         thres_stored=thre2
c
c zero out number of stored integrals and quartets
c
         integrals=0
         nqstore=0
         nblstore=0
         xntegrals=0
         integ62=0
         integ72=0
         integx2=0
         integx4=0
         integx8=0
c
c make denspar(ics,jcs)=1.0
c
         call getmem(ncs*ncs,idensp)
         call setup_denschw(ncs,bl(idensp))
c
 110     continue
c
c request work from the master, send number of integrals and quartets stored
c
         call para_initsend
         call para_pack_int(integrals,1)
         call para_pack_int(nqstore,1)
         call para_pack_int(nblstore,1)
c
         if(indisk.gt.0) then
            call para_pack_real(xntegrals,1)
            call para_pack_real(integ62,1)
            call para_pack_real(integ72,1)
            call para_pack_real(integx2,1)
            call para_pack_real(integx4,1)
            call para_pack_real(integx8,1)
         endif
c
         call para_send_pack(0,TxStorReq)
c get work(starting and finishing block, switch if all negatives were stored)
         call para_recv_pack(isource,TxStorAssign)
         call para_unpack_int(ifrom,1)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
ckw
         call para_unpack_int(ipack,1)
         call para_unpack_int(iqpac,1)
         call para_unpack_int(ibpac,1)
         call para_unpack_int(maxsi,1)
ckw
c
         call setival('intpack',ipack)
         call setival('qrtpack',iqpac)
         call setival('blkpack',ibpac)
         call setival('maxsint',maxsi)
c store
         call int_store(ifrom,ito,bl,bl(ictr),bl(idensp),thre2,
     1                  bl(labels),isecond)
c
cd     call date(datestr)
cd     write(98,*) MY_HOSTNM,' starts ',mywork,' block at ',datestr
cd     call f_lush_(98)
c stop storing if nothing was sent
         if(ito.ge.ifrom) goto 110
c
c find out if positive superblocks were stored and the last stored
c
         call para_send(MY_GID,0,TxStore)
         call para_recv_pack(isource,TxStore)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
c
         call retmem(1) ! release memory allocated for idensp
      endif
c
c
c
 113  continue
c
      ncs2=2*ncs*ncs
      call getmem(ncs2,idscree)
c
      call zeroit(bl(lfxyz),ntri3)
c
      call para_recv_bcastpack(TxPostDens)
      call para_unpack_real(bl(lden),ntri3)
      call para_unpack_real(bl(idscree),ncs2 )
      call para_unpack_int(mgo,1)
c
      if(mgo.eq.1) thresh=thre3
      if(mgo.eq.2) thresh=thre2*10.d0
      if(mgo.eq.3) thresh=thre2
      if(mgo.eq.4) thresh=thre2*0.1d0
      if(mgo.eq.5) thresh=thre2*0.01d0
c make sure that stored integrals are used once each cycle
      isto=1
c
c calculate cphf integrals
c
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0.and.isto.eq.0) exit ! make sure that stored 
                                             ! integrals are used
cd     write(98,9993) mywork,mywork+igran-1
 9993 format(' cphf_integral blocks from',i10,' to:',i10,' started')
cd     call f_lush_(98)
c
c
        call int_cphf(nblocks,bl,bl(ictr),ntri,thresh,bl(idscree),
     *              bl(lden),bl(lfxyz),bl(labels),mywork,igran1)

        if(mywork.le.0)exit ! stored integrals have been used

      enddo
c
c do a global sum of the results
c
      call para_reduce(bl(lfxyz),ntri3,TxPostFock)
c
c
      call retmem(1)  ! allocated for idscree
c
c
c Find out whether another iteration is needed
c
 115  continue
      call para_recv_bcast(mgo,TxNext)
c
      if(mgo.ne.0) goto 113
c
c calculate the NMR shieldings atom by atom
c
      if(iforwhat.eq.2) then
         call getmem(ntri,ldn0) !for D00
         call getmem(18,lforc) !for the shielding tensors
         call para_recv_bcastreal(bl(ldn0),ntri,TxNext) !D00
         call para_recv_bcastreal(bl(lden),ntri3,TxNext) !D01
c
c send initial request w/o results
c
         call para_initsend
         call para_pack_int(0,1) ! no results yet
         call para_send_pack(0,TxShReq)
 116     continue
         call zeroit(bl(lforc),18) !zero out shieldings
         call para_recv_pack(isource,TxShAssign)
         call para_unpack_int(nrat,1)!atom no
         if(nrat.eq.0) goto 117  !we are ready
         call para_unpack_int(nra,1) !no among symm unique
c -- WARNING:  somebody modified <shield> to include two new parameters
         call shield(na,bl,bl(ictr),last,0,6,bl(ibas),bl(inuc),ncs,
     *        ncf,ntri,1,bl(ldn0),bl(lden),bl(lforc),bl(lfxyz),1,nrat,
     *        1,1)
         call para_initsend
         call para_pack_int(1,1) ! results
         call para_pack_int(nra,1) ! which among unique atoms
         call para_pack_real(bl(lforc),18)
         call para_send_pack(0,TxShReq)
         goto 116
 117     continue
         call para_recv_bcast(mgo,TxNext)
      endif
c
c send timings back to master then move on
c
      call secund(tend)
      call para_initsend
      call para_pack_real(tend-tstart,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
c
      call getival('indisk',indisk)
      if(indisk.gt.0) then
         call clos4int()
      endif
c
cd     write(98,*)'POSTSCF INTEGRALS : END     '
cd     call f_lush_(98)
c
c release allocated memory
c
      call retmem(3)
c
c -- release all memory
      call retmark
c
cd     call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
cd     write(98 ,1100) nreq,nmark,lastadr,mxmem,memtot
cd     call f_lush_(98)
 1100 format(' Memory status:'//
     *' request number=',i4,' memory marks=',i3,
     *' last used address=',i9,/,
     *' high water=',i9,' total available memory=',i9/)
      end
c=================================================
c
      subroutine do_prop
c
c=================================================

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      Logical rhf,spin,efg
c     common /big/ bl(300)
c     common /intbl/ maxsh,inx(100)
c
c
c -- get the corresponding host name
      call secund(tstart)
c
c -- unpack the initial data
cc      write(6,*) ' On Slave ',MY_GID,'  About to get TxPropInit'
cc      call f_lush(6)
      call para_recv_bcastpack(TxPropInit)
      call para_unpack_int(NAlpha,1)
      call para_unpack_int(NBeta,1)
      call para_unpack_int(rhf,1)
      call para_unpack_int(spin,1)
      call para_unpack_int(efg,1)
      call para_unpack_int(NRad,1)
      call para_unpack_int(NAng,1)
      call para_unpack_int(LMax,1)
      call para_unpack_real(factor,1)
      call para_unpack_real(thrsh,1)
cc      write(6,*) ' On Slave ',MY_GID,'  NAlpha:',nalpha,' NBeta:',nbeta,
cc     $           ' spin: ',spin,' rhf: ',rhf
cc      write(6,*) ' NRad:',nrad,' NAng:',nang,' LMAX:',lmax,
cc     $           ' factor:',factor,' thrsh:',thrsh
cc      call f_lush(6)
c
c -- put dowm a memory mark
      call mmark
c     call immark
c
c -- unpack the geometry,cpu,symmetry data
      call para_get_info
c
c -- get some necessary basic data
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('ictr',ictr)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('nsym',nsym)
c was undefined, but needed      
      call getival('iprn',iprnt)
c
c -- get number of dummy atoms (if any)
      call getival('ndum',ndum)
      natom = na-ndum
c
c -- set maximum number of radial grid points
      numR = MAX(NRad,400)
c
c -- restore coordinates, atomic charges and get symmetry data
      call getmem(3*natom,ixc)                ! nuclear coordinates
      call getmem(natom,iqa)                  ! atomic charges
      call getmem(natom,ian)                  ! atomic numbers
c
      call getnucdat(natom,bl(ixc),bl(ian))
      call getatchrg(natom,bl(iqa))
c
c  Precalculate data for the numerical grid
c  ----------------------------------------
c
      call getmem(ncs,iexp)                   ! minimum exponents per shell
      call getmem(nsh,ipre)                   ! precomputed normalization
      call getmem(natom**2,iaij)              ! precomputed Becke weights
      call getmem(natom,idst)                 ! nearest neighbour distances
      call getmem(natom**2,irst)              ! inverse interatomic distances
      call getmem(numR,irad)                  ! radial grid points
      call getmem(numR,irwt)                  ! radial weights
      call getmem(3*1130,ixxa)                ! angular quadrature grid
      call getmem(1130,iwta)                  ! angular weights
c
c -- calculate inverse atomic distances, Becke aij parameters
c -- and angular grid and weights prior to full grid construction
       call PreGRID(natom,  bl(ixc), bl(ian), rdum,    rdum,
     $             bl(idst),bl(irst),bl(iaij),bl(ixxa),bl(iwta))
c
c -- get array of smallest exponent per shell
       call GetEXP(nsh,ncs,bl(ibas),bl(ictr),bl(iexp))
c
c -- precompute shell normalization factors
       call AOInit(ncs,bl(ibas),bl(ictr),bl(ipre))
c
c -- get memory for MOs
      call getmem(NAlpha*ncf,icmoA)
      If(.NOT.rhf) call getmem(NBeta*ncf,icmoB)
c
c -- get scratch memory
      NScr =  ncf**2 + ncf + 13*ncs + 5*natom
      If(.NOT.rhf) NScr = NScr + ncf + 3*NBeta + 3
      call getmem(NScr,iscr)
c
c -- unpack MO data
cc      write(6,*) ' On Slave ',MY_GID,'  About to get TxPropDat'
cc      call f_lush(6)
      call para_recv_bcastreal(bl(icmoA),NAlpha*ncf,TxPropDat)
      If(NBeta.GT.0)
     $  call para_recv_bcastreal(bl(icmoB),NBeta*ncf,TxPropDat)
cc      write(6,*) ' On Slave ',MY_GID,'  got TxPropDat'
cc      call f_lush(6)
c
      If(nsym.gt.0) Then
        call getival('nsyo',nsy)
        call getival('ngener',ngen)
        call getival('SymNuPr1',nupr)
      Else
        ngen = 0
      EndIf
C
C
C  Charge/Spin Density
C  -------------------
C
      IF(spin) THEN
c
c -- allocate memory
        call mmark
        LL = (LMax+1)*(2*LMax+1)
        call getmem(NAlpha*7,isyA)            ! MO symmetry transformation
        call getmem(LL,iylm)                  ! spherical harmonics
        call getmem(NAlpha*LL,icA)            ! spherically-averaged MOs
        call getmem(2*natom,ichg)             ! charge density
        call zeroit(bl(ichg),2*natom)         ! zero charge density array
        If(.NOT.rhf) Then
          call getmem(NBeta*7,isyB)           ! beta MO symmetry transformation
          call getmem(NBeta*LL,icB)           ! spherically-averaged beta MOs
          call getmem(2*natom,ispn)           ! spin density
          call zeroit(bl(ispn),2*natom)       ! zero spin density array
        EndIf
c
c -- unpack charge/spin data
cc      write(6,*) ' On Slave ',MY_GID,'  About to get TxSpinDat'
cc      call f_lush(6)
        call para_recv_bcastpack(TxSpinDat)
        call para_unpack_int(MSym,1)
        call para_unpack_real(r0f,1)
        call para_unpack_real(bl(isyA),NAlpha*7)
        If(NBeta.GT.0) call para_unpack_real(bl(isyB),NBeta*7)
cc      write(6,*) ' On Slave ',MY_GID,'  got TxSpinDat'
cc      write(6,*) ' On Slave ',MY_GID,'  LMax:',lmax,' MSym:',msym,
cc     $           ' r0f:',r0f
cc      call f_lush(6)
c
c
 100   CONTINUE
c
c -- request atom we are currently doing from master
cc      write(6,*) ' On Slave ',MY_GID,'  About to get atom'
cc      call f_lush(6)
        call para_get_atom(ICntr)
cc      write(6,*) ' On Slave ',MY_GID,'  got atom',icntr
cc      call f_lush(6)
c
        IF(ICntr.gt.0) THEN
          If(rhf) Then
          CALL NucSPINC(ICntr,  natom,  bl(ixc),bl(ian),bl(iqa),
     $                  MSym,   NGen,   bl(nsy),bl(nupr), IPrnt,
     $                  NRad,   NAng,   factor, r0f,   bl(idst),
     $                  bl(iaij),bl(irst),bl(ixxa),bl(iwta),thrsh,
     $                  bl(irad),bl(irwt),ncf,   ncs, bl(ibas),
     $                  bl(ictr),bl(ipre),bl(iexp),LMax, bl(iylm),
     $                  NAlpha, bl(icmoA),bl(icA),bl(isyA), NScr,
     $                  bl(iscr),bl(ichg))
          Else
          CALL NucSPINU(ICntr,  natom,  bl(ixc),bl(ian),bl(iqa),
     $                  MSym,   NGen,   bl(nsy),bl(nupr), IPrnt,
     $                  NRad,   NAng,   factor, r0f,   bl(idst),
     $                  bl(iaij),bl(irst),bl(ixxa),bl(iwta),thrsh,
     $                  bl(irad),bl(irwt),ncf,   ncs, bl(ibas),
     $                  bl(ictr),bl(ipre),bl(iexp),LMax, bl(iylm),
     $                  NAlpha, NBeta, bl(icmoA),bl(icmoB),bl(icA),
     $                  bl(icB),bl(isyA),bl(isyB),NScr, bl(iscr),
     $                  bl(ichg),bl(ispn))
          EndIf
          GO TO 100
        ENDIF
c
c -- we are done
c -- send results back to master
        call para_reduce(bl(ichg),2*natom,TxProp1)
        If(.not.rhf) call para_reduce(bl(ispn),2*natom,TxProp2)
c
c -- return memory
        call retmark
c
      ENDIF
C
C
C  Electric Field Gradient
C  -----------------------
C
      IF(efg) THEN
c
c -- allocate memory
        call mmark
c       call immark
        call getint(natom,iunq)               ! symmetry-unique atoms
        call getmem(9*natom,iefg)             ! electric field gradient
        call zeroit(bl(iefg),9*natom)
c
c -- unpack EFG data
        call para_recv_bcastpack(TxEFGDat)
        call para_unpack_int(IradQ,1)
        call para_unpack_int(NQ,1)
        call para_unpack_int(bl(iunq),natom)
c
c *******************************************
c -- Symmetry not working properly
c -- switch off until fixed
        call TempSYM(natom,nsym,ngen,NQ,bl(iunq))
c *******************************************
c
c
 200   CONTINUE
c
c -- request atom we are currently doing from master
cc      write(6,*) ' On Slave ',MY_GID,'  About to get atom'
        call para_get_atom(ICntr)
cc      write(6,*) ' On Slave ',MY_GID,'  got atom'
c
        If(ICntr.gt.0) Then
          CALL NucEFG(ICntr,  natom, bl(ixc),bl(ian),bl(iqa),
     $                NSym,   NGen,   bl(nsy), bl(nupr), NQ,
     $                bl(iunq),IPrnt,  NRad,   NAng,   IradQ,
     $                factor, bl(idst),bl(iaij),bl(irst),bl(ixxa),
     $                bl(iwta),thrsh,bl(irad),bl(irwt), ncf,
     $                ncs,    bl(ibas),bl(ictr),bl(ipre),bl(iexp),
     $                rhf,    NAlpha, NBeta, bl(icmoA),bl(icmoB),
     $                NScr,   bl(iscr), bl(iefg))
          GO TO 200
        EndIf
c
c -- we are done
c -- send results back to master
        call para_reduce(bl(iefg),9*natom,TxProp3)                
c
c -- return memory
c       call retimark
        call retmark
c
      ENDIF
c
c
      call para_recv_bcast(mgo,TxNext)
c
c -- send timings back to master then move on
      call secund(tend)
      call para_initsend
      call para_pack_real(tend-tstart,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
c
c -- release allocated memory
c     call retimark
      call retmark
c
      RETURN
      END
c=================================================
c
      subroutine do_hess
c
c=================================================

      use memory
      use newpara
      
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 screen
      character*256 scrfile
      Logical rhf
      Logical do_allat
c     common /big/ bl(300)
c     common /intbl/ maxsh,inx(100)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,
     1              ncs,nsy(4),nsym,nganz(35),lopt(30)
      common /datstore/ thres_stored,isto,isecond,ito
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
      integer*8 xntegrals,integ62,integ72,integx2,integx4,integx8
      common /datadisk/xntegrals,integ62,integ72,integx2,integx4,integx8
      common /screen_type/ screen
c ...................................................................
      common /dftpnt/ ird,iaij,idst,iixx,iiwt,ipre,iexp,
     $                ixgd,iwt,iscr,icofm,nscr
c
      COMMON /DFTCoeff/ aXX(18) ! contains 18 density functional coefficients
c ...................................................................
c
      dimension xintxx(9)
      data iprnt/0/, xlvsh/0.0d0/
c
c -- set memory mark
      call mmark
c
c get the corresponding host name
      call secund(tstart)
c
c get the initial data
c
      call para_recv_bcastpack(TxHessInit)
      call para_unpack_int(ncache,1)
      call para_unpack_int(iforwhat,1)
      call para_unpack_int(NAlpha,1)
      call para_unpack_int(NBeta,1)
      call para_unpack_int(rhf,1)
      call para_unpack_real(threg2,1)
      call para_unpack_real(thref1,1)
      call para_unpack_real(thref2,1)
c
c -- unpack the geometry,cpu,symmetry data
      call para_get_info
c
c setup the common /ganz/ needed later
c
      call getival('na  ',na)
      call getival('ndum',ndum)
      na = na-ndum        ! real atoms only for Hessian evaluation
      call setival('na  ',na)
c
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf ',nbf)
      call getival('nsh ',nsh)
      call getival('lcore',lcore)
      call getival('nsym',nsym)
c
      call getival('ictr',ictr)
      call setival('printh',iprnt)
c
c ---------------------------------------------------------------
c -- set up common /negxyz/ symmetry arrays here
      if(nsym.gt.0) then
         call getival('nsyo',nsyo)
         call make_negxyz(nsym,bl(nsyo))
      endif
c ---------------------------------------------------------------
c
c -- allocate memory for density matrices
      nat3 = 3*na
      ntri = (ncf*(ncf+1))/2
      ntri3 = 3*ntri
      call getmem(ntri,lden)
      mden = lden
      If(.not.rhf) call getmem(ntri,mden)
ckw
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,ncenter)
ckw
c
c -- get density
c
      call para_recv_bcastreal(bl(lden),ntri,TxHess1)
      If(.NOT.rhf) call para_recv_bcastreal(bl(mden),ntri,TxHess1)
      call para_recv_bcastpack(TxHess2)
      call para_unpack_real(bl(idensp),ncs*ncs)
      call para_unpack_int(bl(ncenter),ncf)
      call para_unpack_int(idft,1)
      call para_unpack_real(aXX,18)
c
c -----------------------------------------------------------------
c  G(D0,gxy) - DIRECT CONTRIBUTION TO HESSIAN MATRIX
c -----------------------------------------------------------------
c
c -- set memory mark
      call mmark
c
c prepare the integral system
c
      ax = aXX(1)    ! amount of "exact exchange"
      thresh=threg2
      call post_scf(bl,bl(ictr),ncache,iprnt,thresh,
     *              iforwhat,0,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
c -- get memory for integrals
      call getmem(maxlabels,labels)
c
c -- get memory for Hessian and zero out
      call getmem(nat3**2,lhess)
      call zeroit(bl(lhess),nat3**2)
c
c calculate 2nd derivative integrals
c
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0) exit
c
      call int_d0g2(idft,ax,rhf,nblocks,bl,bl(ictr),ntri,threg2,
     *        bl(idensp),bl(ncenter),
     *        bl(lden),bl(mden),bl(lhess),bl(labels),mywork,igran1)
      enddo
c
c do a global sum of the results
c
      call para_reduce(bl(lhess),nat3**2,TxHess1)
c
c we are ready with contribution of integral 2nd derivatives to Hessian
c
c -- release allocated memory
      call retmark
c
c -----------------------------------------------------------------
c  G(D0,gx) - CONTRIBUTION TO DERIVATIVE FOCK MATRICES
c -----------------------------------------------------------------
c
c prepare the integral system
c
      iforwhat=3
      call post_scf(bl,bl(ictr),ncache,iprnt,thref1,
     *              iforwhat,0,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
c -- get memory for integrals
      call getmem(maxlabels,labels)
c -- get memory symmetry info arrays
      call getint(na,listreal)
c
 300  CONTINUE
cc      write(6,*) ' SLAVE:',MY_GID,' AT 300 CONTINUE'
cc      call f_lush(6)
c
      call mmark
c
c -- get data from master
c
      call para_recv_bcastpack(TxHessDat)
      call para_unpack_int(natb,1)
      call para_unpack_int(nate,1)
ckw
c -- are we done?
      If(nate.EQ.-1) GO TO 301
ckw
      call para_unpack_real(bl(idensp),2*ncs*ncs)
      call para_unpack_int(bl(listreal),na)
      call para_unpack_int(do_allat,1)
ckw
cc      write(6,*) ' SLAVE:',MY_GID,' na=',na,' natb,nate=',natb,nate
cc      write(6,*) '         do_allat=',do_allat
cc      call f_lush(6)
c
c -- are we done?
c     If(nate.EQ.-1) GO TO 301
c
      natonce=nate-natb+1
c
c -- get memory for derivative Fock matrices and zero out
c
      call getmem(natonce*3*ntri,lfock1)
      call zeroit(bl(lfock1),natonce*3*ntri)
      lfockB = lfock1
      If(.NOT.rhf) Then
        call getmem(natonce*3*ntri,lfockB)
        call zeroit(bl(lfockB),natonce*3*ntri)
      EndIf
c
c calculate 1st derivative integrals
c
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0) exit
c
      call int_d0g1(idft,ax,rhf,nblocks,bl,bl(ictr),ntri,thref1,
     *              natb,nate,bl(listreal),do_allat,
     *        bl(idensp),bl(ncenter),
     *        bl(lden),bl(mden),bl(lfock1),bl(lfockB),bl(labels),
     *        mywork,igran1)
      enddo
c
c do a global sum of the results
      call para_reduce(bl(lfock1),3*ntri*natonce,TxHess2)
c -- now send back beta derivative Fock matrix
      If(.NOT.rhf) call para_reduce(bl(lfockB),3*ntri*natonce,TxHess3)
c --------------------------------------------------------------------
c
c we are ready with contribution of integral 1st derivatives to
c derivative Fock matrices
c
c -- release allocated memory
      call retmark
c
      go to 300
c
  301 continue
c -- release memory for integrals
      call retmem(1)
c -- release memory symmetry info arrays
      call retmem(1)
c--------------------------------------------------------
c ........................................................
c  DFT PART
c  (contributions to Fock derivative matrices and direct
c   contribution to Hessian matrix)
c
      IF(idft.NE.0) THEN
c
c -- get initial DFT data
        IEntry = 1                   ! set for first entry
        call getival('ibas',ibss)    ! god knows what value this has
        call para_recv_bcastpack(TxDftInit)
        call para_unpack_int(nrad,1)
        call para_unpack_int(nang,1)
        call para_unpack_int(lrad,1)
        call para_unpack_int(lang,1)
        call para_unpack_int(IradQ,1)
        call para_unpack_int(NBatch,1)
        call para_unpack_int(IdWt,1)
c .........................................................
c -- these next lines needed because the BASDAT
c -- array is corrupted in the slaves
        call para_unpack_int(ncs,1)
        call para_unpack_int(nsh,1)
        call para_unpack_int(bl(ictr),12*ncs)
        call para_unpack_real(bl(ibss),13*nsh)
c .........................................................
        call para_unpack_real(factor,1)
        call para_unpack_real(thrsh,1)
c
c -- ensure correct values are in depository
        call setival('nmo ',NAlpha)
        call setival('ncs ',ncs)
        call setival('nsh ',nsh)
        call setival('wghtd',IdWt)
c
c -- get number of dummy atoms (if any)
c    (NO need to do this now, as "na" has correct value)   ! JB  Oct 2006
cc        call getival('ndum',ndum)
cc        natom = na-ndum
        natom = na
c
c -- restore coordinates and atomic charges
        call mmark
c       call immark
        call getmem(3*natom,ixnc)            !  nuclear coordinates
        call getmem(natom,iqa)               !  atomic charges/numbers
c
        call getnucdat(natom,bl(ixnc),bl(iqa))
c
        ityp = 5
        if(.NOT.rhf) ityp = 7
c -- DAMN!! With "na" now being the number of REAL atoms   ! JB  Oct 2006
c -- we have to reset it to all centres for the DFT initialization
        na = na+ndum        ! all atoms only for DFT setup
        call setival('na  ',na)
        call setup_dft(idft,   ityp,   1,      nrad,   nang,
     $                 NBatch, bl(iqa),bl(ixnc))
c -- restore "na"
        na = na-ndum        ! real atoms only
        call setival('na  ',na)
c
c -- make sure we have all pointers
        If(nsym.gt.0) Then
          call getival('nsyo',isy)
          call getival('ngener',ngen)
          call getival('SymNuPr1',nupr)
          call getival('SymFunPr',ifp)
        EndIf
c
        call getint(ncf,iatm)                !  no. basis functions per atom
        call getmem(ncf,idm)                 !  maximum density per row/column
c
c -- unpack second batch of DFT data
        call para_recv_bcastpack(TxDftDat)
        call para_unpack_int(bl(iatm),ncf)
        call para_unpack_real(bl(lden),ntri)
        If(.NOT.rhf) call para_unpack_real(bl(mden),ntri)
        call para_unpack_real(bl(idm),ncf)
c
c -- multiple passes loop.
c
 250    CONTINUE
c
c -- unpack total number of passes, current pass, first and last
c    component to compute
c
        call para_recv_bcastpack(TxDftHmp)
        call para_unpack_int(npass,1)
        call para_unpack_int(ipass,1)
        call para_unpack_int(Nb,1)
        call para_unpack_int(Ne,1)
cc        write(44,*)'npass',npass,' ipass',ipass,' nb',nb,' ne',ne
c
c -- get memory for Hessian and derivative Fock matrices and zero out
        natcurr=Ne-Nb+1
        call mmark
c       call immark
        call getmem(nat3*nat3,lhess)
        call zeroit(bl(lhess),nat3*nat3)
        call getmem(3*natcurr*ntri,lfock1)
        call zeroit(bl(lfock1),3*natcurr*ntri)
        If(.NOT.rhf) Then
          call getmem(3*natcurr*ntri,lfockB)
          call zeroit(bl(lfockB),3*natcurr*ntri)
        EndIf
        El=0.0d0
c
 100    CONTINUE
c
c -- request atom we are currently doing from master
        call para_get_atom(ICntr)
c
        IF(ICntr.gt.0) THEN
          If(rhf) Then
          CALL DFTHessC(idft,   ICntr,   natom,     Nb,       Ne,
     $                bl(ixnc), bl(iqa), nsym,      ngen,     bl(isy),
     $                bl(nupr),0,       lrad,      lang,     IradQ,
     $                factor,   NBatch,  bl(idst),  bl(iaij), bl(ird),
     $                bl(iixx), bl(iiwt),bl(ixgd),  bl(iwt),  thrsh,
     $                ncf,      ncs,     bl(ibss),  bl(ictr),bl(ipre),
     $                bl(iexp),IdWt,     bl(iatm), bl(lden), bl(idm),
     $                nscr,    bl(iscr), bl(lfock1),bl(lhess),el,
     $                IEntry)
          Else
          CALL DFTHessU(idft, ICntr,   natom,    Nb,        Ne,
     $              bl(ixnc), bl(iqa), nsym,     ngen,      bl(isy),
     $              bl(nupr),0,       lrad,     lang,      IradQ,
     $              factor,   NBatch,  bl(idst), bl(iaij),  bl(ird),
     $              bl(iixx), bl(iiwt),bl(ixgd), bl(iwt),   thrsh,
     $              ncf,      ncs,     bl(ibss), bl(ictr), bl(ipre),
     $              bl(iexp), IdWt,    bl(iatm),bl(lden),  bl(mden),
     $              bl(idm),  nscr,    bl(iscr), bl(lfock1),bl(lfockB),
     $              bl(lhess),el,      IEntry)
          EndIf
          GO TO 100
        ENDIF
c
c -- we are done
c -- send results back to master
        call para_reduce(bl(lhess),nat3**2,TxHess4)
        call para_reduce(el,1,TxDftF4)
c --------------------------------------------------------------------
        call para_reduce(bl(lfock1),3*ntri*natcurr,TxHess2)
c -- send back beta Fock matrix derivatives
        If(.NOT.rhf) call para_reduce(bl(lfockB),3*ntri*natcurr,TxHess3)
c
c --  ready for next pass
c
c       call retimark
        call retmark
        if(ipass.ne.npass)goto 250
c --------------------------------------------------------------------
c
c -- remove DFT memory
c       call retimark
        call retmark
        call retmark            ! added as memory leak found    Feb 2009
cc
      ENDIF
c ........................................................
c
c -- release memory for density matrices
      call retmem(1)
      if(.not.rhf) call retmem(1)
c
c
c -----------------------------------------------------------------
c  G(D1,g0) - CONTRIBUTION TO CPHF MATRICES
c -----------------------------------------------------------------
c
c prepare the integral system
c
      iforwhat=1
      call post_scf(bl,bl(ictr),ncache,iprnt,thref1,
     *              iforwhat,1,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
c -- get memory for integrals
      call getmem(maxlabels,labels)
c
c
 200  CONTINUE
cc      write(6,*) ' SLAVE:',MY_GID,' AT 200 CONTINUE'
cc      call f_lush(6)
c
c .....................................................................
c  ** WARNING **  Level shifting has NOT been implemented in parallel
c .....................................................................
c
c -- get number of atoms for CPHF from master
      call para_recv_bcast(natonce,TxHessDat)
c
c -- are we done this pass?
      If(natonce.EQ.-1) GO TO 350
c -- are we done completely?
      If(natonce.EQ.-2) GO TO 999
c
c
c .................................................................
c   CONSTANT PART OF FIRST-ORDER DENSITY & FDS1
c .................................................................
c
c -- get VCD flag
      call para_recv_bcast(ivcd,TxHessDat)
c
c -- allocate memory for constant data
      call getmem(ncf,ind)
      call setup_lind(bl(ind),ncf)
      call getmem(ntri,lden)
      call getmem(ncf**2,lfd)
      call getmem(ncf**2,ivec)
      call getmem(ncf,ival)
      If(.NOT.rhf) Then
        call getmem(ntri,ldenB)
        call getmem(ncf**2,lfdB)
        call getmem(ncf**2,ivecB)
        call getmem(ncf,ivalB)
      EndIf
c
c -- first get constant data
      call para_recv_bcastreal(bl(lden),ntri,TxWdens1)
      call para_recv_bcastreal(bl(lfd),ncf**2,TxWdens1)
      call para_recv_bcastreal(bl(ivec),ncf**2,TxWdens1)
      call para_recv_bcastreal(bl(ival),ncf,TxWdens1)
      If(.NOT.rhf) Then
        call para_recv_bcastreal(bl(ldenB),ntri,TxWdens1)
        call para_recv_bcastreal(bl(lfdB),ncf**2,TxWdens1)
        call para_recv_bcastreal(bl(ivecB),ncf**2,TxWdens1)
        call para_recv_bcastreal(bl(ivalB),ncf,TxWdens1)
      EndIf
c
c -- allocate memory for 1st-order data
      call getmem(ntri3,lfock1)
      call getmem(ntri3,lover1)
      call getmem(ntri3,ir1)
      call getmem(ncf*ncf*3,ifds)
      If(.NOT.rhf) Then
       call getmem(ntri3,lfockB)
       call getmem(ntri3,ir1B)
       call getmem(ncf*ncf*3,ifdsB)
      EndIf
      If(ivcd.NE.0) call getmem(ncf*NAlpha*3,lpve)
c
c send initial request w/o results
c
      call para_send(-1,0,TxWdens2)
c
c -- get first-order matrices and atom
 400  CONTINUE
      call para_recv(IAtm,isource,TxWDens3)
      If(IAtm.EQ.-1) GO TO 495       ! we are done
      call para_recv_real(bl(lfock1),ntri3,0,TxWDens3)
      call para_recv_real(bl(lover1),ntri3,0,TxWDens3)
      If(.NOT.rhf) call para_recv_real(bl(lfockB),ntri3,0,TxWDens3)
c
c .... alpha part
c -- calculate FDS1 and S1DF matrices and F1-(FDS1+S1DF)
c
      call calcF1FDS1(ncf,    ntri,   bl(lfd),bl(lover1),bl(lfock1),
     *                bl(ifds))
c
c -- calculate constant part of 1st-order density matrix
      call d1const_xyz(rhf,    ncf,    ntri,   NAlpha, bl(ind),
     $                bl(ivec),bl(ival),bl(lden),bl(lfock1),bl(lover1),
     $                 bl(ir1),ivcd,   bl(lpve),bl)
c
      If(.NOT.rhf) Then
c .... beta part
c -- calculate FDS1 and S1DF matrices and F1-(FDS1+S1DF)
      call calcF1FDS1(ncf,    ntri,   bl(lfdB),bl(lover1),bl(lfockB),
     $                bl(ifdsB))
c
c -- calculate constant part of 1st-order density matrix
      call d1const_xyz(rhf,    ncf,    ntri,   NBeta,  bl(ind),
     $             bl(ivecB),bl(ivalB),bl(ldenB),bl(lfockB),bl(lover1),
     $                 bl(ir1B),0,     bl,     bl)   ! NO VCD in UHF
      EndIf
c
c -- send results back to master
      call para_send(IAtm,0,TxWdens2)
      call para_send_real(bl(ir1),ntri3,0,TxWDens2)
      call para_send_real(bl(ifds),ncf*ncf*3,0,TxWDens2)
      If(.NOT.rhf) Then
        call para_send_real(bl(ir1B),ntri3,0,TxWDens2)
        call para_send_real(bl(ifdsB),ncf*ncf*3,0,TxWDens2)
      EndIf
      If(ivcd.NE.0)
     *  call para_send_real(bl(lpve),ncf*NAlpha*3,0,TxWDens2)
      GO TO 400
c
 495  CONTINUE
      If(.NOT.rhf) call retmem(7)
      call retmem(9)
      If(ivcd.NE.0) call retmem(1)
c
c .................................................................
c        CPHF PROPER
c .................................................................
c
c
c do the integral storage if necessary
c
      if(scftype.ne.'full-direct') then
c
c open disk files for stored integrals if needed
c
        call getival('indisk',indisk)
        if(indisk.ne.0) then
           call getchval('scrf',scrfile)
           call open4store
        endif
c
c save integral threshold used to calculate stored integrals
c (in common /datstore/ thres_stored   )
c
         thres_stored=thref1
c
c zero out number of stored integrals and quartets
c
         integrals=0
         nqstore=0
         nblstore=0
         xntegrals=0
         integ62=0
         integ72=0
         integx2=0
         integx4=0
         integx8=0
c
c make denspar(ics,jcs)=1.0
c
         call getmem(ncs*ncs,idensp)
         call setup_denschw(ncs,bl(idensp))
c
 110     continue
c
c request work from the master, send number of integrals and quartets stored
c
         call para_initsend
         call para_pack_int(integrals,1)
         call para_pack_int(nqstore,1)
         call para_pack_int(nblstore,1)
c
         if(indisk.gt.0) then
            call para_pack_real(xntegrals,1)
            call para_pack_real(integ62,1)
            call para_pack_real(integ72,1)
            call para_pack_real(integx2,1)
            call para_pack_real(integx4,1)
            call para_pack_real(integx8,1)
         endif
c
         call para_send_pack(0,TxStorReq)
c get work(starting and finishing block, switch if all negatives were stored)
         call para_recv_pack(isource,TxStorAssign)
         call para_unpack_int(ifrom,1)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
ckw
         call para_unpack_int(ipack,1)
         call para_unpack_int(iqpac,1)
         call para_unpack_int(ibpac,1)
         call para_unpack_int(maxsi,1)
ckw
c
         call setival('intpack',ipack)
         call setival('qrtpack',iqpac)
         call setival('blkpack',ibpac)
         call setival('maxsint',maxsi)
c store
         call int_store(ifrom,ito,bl,bl(ictr),bl(idensp),thref1,
     1                  bl(labels),isecond)
c
cd     call date(datestr)
cd     write(98,*) MY_HOSTNM,' starts ',mywork,' block at ',datestr
cd     call f_lush_(98)
c stop storing if nothing was sent
         if(ito.ge.ifrom) goto 110
c
c find out if positive superblocks were stored and the last stored
c
         call para_send(MY_GID,0,TxStore)
         call para_recv_pack(isource,TxStore)
         call para_unpack_int(ito,1)
         call para_unpack_int(isecond,1)
c
         call retmem(1) ! release memory allocated for idensp
      endif
c 
      call getmem(ncf,ind)
      call setup_lind(bl(ind),ncf)
 350  CONTINUE
c
c  CPHF Loop
c
cc      write(6,*) ' SLAVE:',MY_GID,' AT 350 CONTINUE - CPHF loop'
cc      call f_lush(6)
c
c -- set memory mark
      call mmark
c     call immark
c
c -- get data for CPHF from master
      call para_recv_bcast(natonce,TxHessDat)
c
c -- are we done this pass?
      If(natonce.EQ.-1) GO TO 201
c
c -- allocate memory for screening density and 1st-order density matrix
      call getint(natonce,iatm)
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,map_fs)
      call getmem(3*ntri*natonce,lden)
      If(.NOT.rhf) call getmem(3*ntri*natonce,mden)
c
      call para_recv_bcastpack(TxHessDat)
      call para_unpack_int(bl(iatm),natonce)
      call para_unpack_string(screen,4)
      call para_unpack_real(thrsh,1)
      call para_unpack_int(bl(map_fs),ncf)
      call para_unpack_real(bl(idensp),2*ncs*ncs)
      call para_recv_bcastreal(bl(lden),3*ntri*natonce,TxHessDat)
      If(.NOT.rhf) call
     *         para_recv_bcastreal(bl(mden),3*ntri*natonce,TxHessDat)
c
c -- get memory for CPHF Fock and density matrices and zero out
      call getmem(3*ntri*natonce,lfock)
      call zeroit(bl(lfock),3*ntri*natonce)
      If(.NOT.rhf) Then
        call getmem(3*ntri*natonce,lfockB)
        call zeroit(bl(lfockB),3*ntri*natonce)
      EndIf
c make sure that stored integrals are used once each cycle
      isto=1
c
c calculate integrals
c
      do
         call para_get_block(mywork,igran1)
         if (mywork.le.0.and.isto.eq.0) exit ! make sure that stored 
                                             ! integrals are used
c
cc      write(6,*) ' SLAVE:',MY_GID,' calling <int_d1g0nat>  mywork:',
cc     $            mywork,' igran1:',igran1
cc      call f_lush(6)
      If(rhf) Then
      call int_d1g0nat(bl(iatm),natonce,idft, ax,     nblocks,
     *                 bl,     bl(ictr),ntri, thrsh, bl(map_fs),
     *              bl(idensp),bl(lden),bl(lfock),bl(labels),mywork,
     *                 igran1)
      Else
      call int_d1g0nat_uhf(bl(iatm),natonce,idft, ax,     nblocks,
     *                     bl,     bl(ictr),ntri, thrsh, bl(map_fs),
     *              bl(idensp),bl(lden),bl(mden),bl(lfock),bl(lfockB),
     *                  bl(labels),mywork,igran1)
      EndIf

        if(mywork.le.0)exit ! stored integrals have been used
c
      enddo
c
c we are ready with contribution of integrals and 1st order
c densities to CPHF Fock matrices
c summing only at the end
c
c -- rescale fock-derivative matrices from coulomb
        if(rhf) then
           call rescale_nat_fock(bl(lfock),ntri,natonce,bl(ind),
     &                           thrsh*0.5d0)
        else
           call rescale_nat_fock(bl(lfock),ntri,natonce,bl(ind),thrsh)
           call rescale_nat_fock(bl(lfockB),ntri,natonce,bl(ind),thrsh)
        endif
c
c ........................................................
c  DFT PART
c  (contributions to partial Fock derivative matrices in CPHF)
c
      IF(idft.NE.0) THEN
c
c -- restore coordinates and atomic charges
        call getmem(3*natom,ixnc)            !  atomic coordinates
        call getmem(natom,iqa)               !  atomic charges/numbers
        call getmem(3*natom,ixs)             !  symmetry-unique atomic coordinates
c
        call getnucdat(natom,bl(ixnc),bl(iqa))
c
c -- allocate memory for density matrix and maximum 1st-order density elements
        call getmem(ntri,ida)                !  density matrix
        If(.NOT.rhf) call getmem(ntri,idb)   !  beta density matrix
        call getmem(ncf,idm)                 !  maximum density per row/column
c
c -- unpack CPHF DFT data
c -- NOTE: most of the required data is already here from Coulomb term
        IEntry = 1                           ! set for first entry
        call para_recv_bcastpack(TxDftDat)
        call para_unpack_int(lgrid,1)
        call para_unpack_int(natdo,1)
        call para_unpack_real(bl(idm),ncf)
c -- this next needed because the BASDAT array (or its pointer)
c -- is corrupted in the slaves   - SOMEBODY MUST FIX THIS
        call para_unpack_real(bl(ibss),13*nsh)
        call para_unpack_int(nq,1)
        call para_unpack_real(bl(ixs),3*nq)
        call para_unpack_int(bl(iqa),nq)
        call para_unpack_real(bl(ida),ntri)
        If(.NOT.rhf) Then
          call para_unpack_real(bl(idb),ntri)
c no need to receive the perturbed densities, they are the same as in
c the integrals
c       call para_recv_bcastreal(bl(lden),3*ntri*natonce,TxDftDat)
c         call para_recv_bcastreal(bl(mden),3*ntri*natonce,TxDftDat)
        EndIf
c
c -- need to redo preliminary grid stuff using symmetry-unique atoms only
c -- temporarily set natoms to nq
        call setival('na  ',nq+ndum)
        ityp = 6
        If(.NOT.rhf) ityp = 8
        call setup_dft(idft,   ityp,   1,      nrad,   nang,
     $                 NBatch, bl(iqa),bl(ixs))
c -- restore true value
        call setival('na  ',na)
c
 150    CONTINUE
c
c -- request atom we are currently doing from master
        call para_get_atom(ICntr)
c
        IF(ICntr.gt.0) THEN
          If(rhf) Then
          CALL DFTCphfC(idft,   ICntr,  nq,   natonce, natdo,
     $                 bl(ixnc),bl(ixs),bl(iqa), 0,    lrad,
     $                  lang,   IradQ,  factor, NBatch,bl(idst),
     $                 bl(iaij),bl(ird),bl(iixx),bl(iiwt),bl(ixgd),
     $                  bl(iwt),thrsh,  ncf,    ncs,   bl(ibss),
     $                 bl(ictr),bl(ipre),bl(iexp),IdWt, bl(ida),
     $                 bl(lden), bl(idm),nscr, bl(iscr),bl(lfock),
     $                  lgrid,  IEntry)
          Else
          CALL DFTCphfU(idft,   ICntr,  nq,   natonce, natdo,
     $                 bl(ixnc),bl(ixs),bl(iqa), 0,    lrad,
     $                  lang,   IradQ,  factor, NBatch,bl(idst),
     $                 bl(iaij),bl(ird),bl(iixx),bl(iiwt),bl(ixgd),
     $                  bl(iwt),thrsh,  ncf,    ncs,   bl(ibss),
     $                 bl(ictr),bl(ipre),bl(iexp),IdWt, bl(ida),
     $                 bl(idb), bl(lden),bl(mden), bl(idm),nscr,
     $                 bl(iscr),bl(lfock),bl(lfockB),lgrid,IEntry)
          EndIf
          GO TO 150
        ENDIF
cc
      ENDIF
c
c -- we are done
c -- send results back to master
c --------------------------------------------------------------------
      call para_reduce(bl(lfock),3*ntri*natonce,TxHess3)
      If(.NOT.rhf) call para_reduce(bl(lfockB),3*ntri*natonce,TxHess2)
c --------------------------------------------------------------------
c
c -- release allocated memory
c     call retimark
      call retmark
c
      GO TO 350
 201  CONTINUE
      call retmem(1)  ! ind
c
c
c .................................................................
c   CONSTRUCTION OF FINAL FIRST-ORDER WEIGHTED DENSITY
c .................................................................
c
c --  memory has been released so reallocate
      call getmem(ncf,ind)
      call setup_lind(bl(ind),ncf)
      call getmem(ntri,lden)
      call getmem(ncf**2,lfd)
      call getmem(3*ncf**2,lscr)
      call getmem(ntri3,lden1)
      call getmem(ntri3,lfock1)
      call getmem(ntri3,lgmat)
      call getmem(ntri3,ir1)
      If(.NOT.rhf) Then
        call getmem(ntri,ldenB)
        call getmem(ncf**2,lfdB)
        call getmem(ntri3,lden1B)
        call getmem(ntri3,lfock1B)
        call getmem(ntri3,lgmatB)
        call getmem(ntri3,is1)
      EndIf
c
c -- first get constant data
      call para_recv_bcastpack(TxWDens1)
      call para_unpack_real(bl(lden),ntri)
      call para_unpack_real(bl(lfd),ncf**2)
      If(.NOT.rhf) Then
        call para_unpack_real(bl(ldenB),ntri)
        call para_unpack_real(bl(lfdB),ncf**2)
      EndIf
c
c send initial request w/o results
c
      call para_initsend
      call para_pack_int(-1,1)
      call para_pack_int(-1,1)
      call para_send_pack(0,TxWDens2)
c
c -- get first-order matrices and atom
 500  CONTINUE
      call para_recv(IEntry,isource,TxWDens3)
      If(IEntry.EQ.-1) GO TO 595       ! we are done
      call para_recv_pack(isource,TxWDens5)
      call para_unpack_int(IAtm,1)
      call para_unpack_real(bl(lfock1),ntri3)
      call para_unpack_real(bl(lgmat),ntri3)
      call para_unpack_real(bl(lden1),ntri3)
      If(.NOT.rhf) Then
        call para_unpack_real(bl(lfock1B),ntri3)
        call para_unpack_real(bl(lgmatB),ntri3)
        call para_unpack_real(bl(lden1B),ntri3)
      EndIf
c
c -- calculate weighted-density
c -- alpha/closed-shell part
      CALL wdens1(rhf,    ncf,    ntri,   bl(ind),bl(lden),
     $            bl(lfd),bl(lfock1),bl(lgmat),bl(lden1),
     $            bl(lscr),bl(ir1))
c
      If(.NOT.rhf) Then
c -- beta part
        CALL wdens1(rhf,    ncf,    ntri,   bl(ind),bl(ldenB),
     $              bl(lfdB),bl(lfock1B),bl(lgmatB),bl(lden1B),
     $              bl(lscr),bl(is1))
c -- sum up alpha & beta
        CALL AddVEC(ntri3,bl(ir1),bl(is1),bl(ir1))
      EnDIf
c
c -- send result back to master
      call para_initsend
      call para_pack_int(IEntry,1)
      call para_pack_int(IAtm,1)
      call para_send_pack(0,TxWDens2)
      call para_send_real(bl(ir1),ntri3,0,TxWDens4)
      GO TO 500
c
 595  CONTINUE
      If(.NOT.rhf) call retmem(6)
      call retmem(8)
c
      GO TO 200
 999  CONTINUE
c
c -- release memory for integrals
      call retmem(1)
c
c -- release all memory
      call retmark
      call retmark            ! added as memory leak found    Feb 2009
c
c send timings back to master then move on
c
      call para_send(MY_GID,0,TxContinue)
c
      call para_recv_bcast(mgo,TxNext)
      call secund(tend)
      call para_initsend
      call para_pack_real(tend-tstart,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
c
c -- if semidirect, remove integral files
      If(scftype.ne.'full-direct') then
        call getival('indisk',indisk)
        if(indisk.gt.0) call clos4int()
      endif
C
      RETURN
      END

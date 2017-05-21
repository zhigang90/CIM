      Subroutine para_rmp2(nmo)

      use memory
      use newpara
      
      implicit real*8 (a-h,o-z)
      integer stat
      logical LastSymPair,emp2only,scs,igran_calculate
      dimension xintxx(9)
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      common /job/jobname,lenJ
      parameter(nopt=16)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $             filname3,filname4
      character*256 chopv(nopt)
c ............................................................
      integer binlist(nslv)
c ............................................................
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20
      logical nofr,restrt,dualbasis,smal,keep,Test
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension xnmo(2)
      parameter (maxpair=700*701/2,maxgran=maxpair/2)   !FRAGILE!
      integer igran2pair(2,maxgran),ipair2gran(maxpair)
c
      data opnames  /'nofr','orbs','thre','core','rest',
     1               'dual','big ','smal','maxd','keep',
     2               'prin','grad','scs ','gran','type',
     3               'test'/
      data ioptyp / 0, 2,11,11, 0, 0, 0, 0,11, 0, 1, 0,12, 1, 1, 0/
c
c  Explanation of the options:
c  NOFR = no frozen core correlate all orbitals
c  ORBS = correlate orbitals from i to j
c  THRE = integral thresold, n pH form, default is  9, meaning 1.0d-9
c  CORE = the orbital energy separating the cores from the valence, in
c         atomic units, default is -3 Eh; if set in RHF, that value is
c         used
c  REST = restart, meaning that the half-transformed integrals are
c         already calculated
c  DUAL = dual basis MP2
c  BIG  = use the 'big' memory model. This is the default for calculations
c         with more than 30 atoms.
c  SMAL = use the 'small memory model. If BIG and SMAL are both given,
c         then BIG is assumed. Default for calculations with < 30 atoms
c  MAXD = total amount of disk space available per processor (in GB)
c  KEEP = do not delete the half-transformed integral file
c  PRIN = print level, default=2
c  GRAD = calculate the extra quantities needed for an MP2 gradient
c         and write them to disk
c  SCS  = use scaled MP2 (with optional read-in of 2 scale factors)
c  GRAN = granule size for writing bins in bunches
c  TYPE = dual-basis type (1, 2 or 3)
c  TEST = undocumented feature. Will stop the job after the first
c         half-transformation. Used to test the restart.
c
      data IUnit/1/, nfock/1/
      parameter (gigabyte=1.0737d9)
c
      call secund(tt0)
      call elapsec(elaps0)
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c  basis function symmetry pairs are in bl(ifp)
      iprnt=2
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      lastvirt=igetival('nonredun')
      ictr=igetival('ictr')
c-----------------------------------------------------------
c .............................................................
c -- CHECK! Only do MP2 if current wavefunction is RHF
c --   (at the same time get the lowest eigenvalue of the overlap)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,wvfnc)
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(wvfnc(1:3).EQ.'RHF' .or. wvfnc(1:3).EQ.'RMP') then
      else
         call nerror(1,'MP2 module',
     $   'MP2 Requested But Current Wavefunc. is NOT RHF or MP2 !',0,0)
      endif
c .............................................................
c
      call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c
c  'NOFR' no frozen core
      nofr=.false.
      if(ifound(1).gt.0) nofr=.true.
c  ORBS - see below
c  THRESHOLD (integral threshold)
      if (ifound(3).gt.0) then
         thresh=10.0d0**(-ropv(1,3))
      else
         thresh = MIN(1.0d-10,xlow**2)
         If (thresh.LT.1.0d-12) then
            write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final MP2 energy may have greater than mhartree error')
            thresh = 1.0d-12
         EndIf
      end if
c  core is the limit separating core from valence;
c  it is normally set in rhf but it can be set here too
      call tstrval('core',iyes)
      if (iyes.eq.1) then
         core=rgetrval('core')
      else
         core=-3.0d0
      end if
      if(ifound(4).gt.0)  core=ropv(1,4)
c  restart
      restrt=.false.
      if (ifound(5).gt.0) then
         restrt=.true.
      endif
c .............................................................................
c -- File handling --
c  filname1=half-transformed integrals; filname2=sorted bins
c  If gradient: filname3=Tij amplitudes;  filname4=Kov amplitudes
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.htr'
      len1 = len+4
      filname2=scrfile(1:len)//'.bins'
      len2 = len+5
c ............................................................................
      dualbasis=ifound(6).gt.0
      natoms=igetival('na')
      if (natoms.le.30) then
         smal=.true.
      else
cc        smal=.false.
         smal=.true.      ! once again, BIG option disable  ! JB
      end if
      if(ifound(7).gt.0) smal=.false.
      if(ifound(8).gt.0) smal=.true.
      if(ifound(9).gt.0) then
         DiskMax = ropv(1,9)
      else
         DiskMax = 20.0d0
      endif
      keep = ifound(10).gt.0
      if(ifound(11).gt.0) iprnt=iopv(1,11)
c -- SCS scaling  ! SS
      scs=.false.
      iscs=0
      p1=1.2d0
      p2=1.0d0/3.0d0
      If(ifound(13).gt.0) Then
         scs=.true.
         iscs=1
         if(ropv(1,13).gt.0) then
            p11=ropv(1,13)
            p22=ropv(2,13)
            if(abs(p1-p11).gt.1.0d-05) then
               iscs=2
               p1=p11
               p2=p22
               write(6,*)'Scaled MP2 with p1= ',p1,'p2= ',p2
            endif
         else
            write(6,*) 'Conventional SCS-MP2'
         endif
      EndIf
      call setival('iscs',iscs)
      call setrval('p1',p1)
      call setrval('p2',p2)
c-----------------------------------------------------------
c -- determine logical variable emp2only
c    normally FORCE is the next line for and MP2 gradient
c    and we are going to check for this automatically here
      READ(inp,'(A4)',END=95) chopv(12)(1:4)
      call lowercas(chopv(12),4)
      If(chopv(12)(1:4).EQ.'forc') ifound(12)=1
      BACKSPACE inp
c
 95   emp2only=.true.
      mp2only = 1
      If(ifound(12).gt.0) Then
         emp2only=.false.
         mp2only=-1
         filname3=scrfile(1:len)//'.Tij'
         len3 = len+4
         filname4=scrfile(1:len)//'.Kov'
         len4 = len+4
      EndIf
      call setival('mp2only',mp2only)
      if(ifound(14).gt.0) then
         igranulesize=iopv(1,14)
         igran_calculate=.false.
      else
         igran_calculate=.true.
      endif
ckw07
c      dual type (1,2 or 3)
c
      itype=1
      If(ifound(15).gt.0) Then
         itype=iopv(1,15)
         if(itype.gt.3) itype=1
      Endif
c
      Test = ifound(16).gt.0
c-----------------------------------------------------------
  45  format( 72('-'))
  46  format(/72('='))
c-----------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
      write(iout,*) '                            The RMP2 Module  '
      write(iout,*)' '
c-----------------------------------------------------------
c allocate memory for an array showing if a big bs shell
c was present (0 or 1) in a small bs
c
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c
c this is needed for dual basis set but is also used
c in blocking for ordinary MP2 (must have all 1)
c-----------------------------------------------------------
      if(dualbasis) then
c        the current "extended" basis set was already read-in
c        after SCF. Get the smaller basis set info (one used for SCF):
c
         call get_small_basis(bl,ncs,ncf,bl(ictr), ncs_sm,ncf_sm)
c
c-----------------------------------------------------------
         write(iout,45)
         write(iout,*)
     *   '              Dual basis set will be used in MP2 calculations'
         write(iout,*) ' '
         write(iout,47) ncs_sm,ncf_sm
         write(iout,48) ncs   ,ncf
   47    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for occupied MOs')
   48    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for virtual  MOs')
         write(iout,45)
         write(iout,49) itype
         write(iout,45)
   49    format(22x,'dual basis set approach type = ',i2)
c-----------------------------------------------------------
c
      else
         ncs_sm=ncs
         ncf_sm=ncf
         call set_mpres(bl(mpres_in_bg),ncs)
      endif
c-----------------------------------------------------------
      write(iout,*) ' MP2 integral thresh    = ',thresh
      write(iout,*) ' core orbitals with eps < ', core
      write(iout,*) ' '
      call f_lush(iout)
c-----------------------------------------------------------
c  put down a memory marker
      call matreset
      call mmark
      call matmark
c     call immark
c
      np1=1
      np4=4
      call sudat(np4,'nocc_rhf',ni)
      if(ni.gt.0) then
         call rea(xnmo,2,np4,'nocc_rhf')
      else
         call restart(np4,0)
         call restart(np1,0)
         call sudat(np4,'nocc_rhf',ni)
         if (ni.gt.0) then
            call rea(xnmo,2,np4,'nocc_rhf')
         else
            call nerror(2,'MP2 module',
     1     'Program cannot find the number of MOs on <jobname>.14',
     2      ni,ni)
         end if
      end if
      nmo=xnmo(1)
      erhf=xnmo(2)
c-----------------------------------------------------------
c zero out irrelevant dft stuff
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
c-------------------------------------------------------
c get memory for MOs, orbital energies and screening density
c
      np4=4
      call matdef('cano','q',ncf,ncf)
      call matdef('epsi','d',ncf,ncf)
      call matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      idics=mataddr('dsmx')
c-------------------------------------------------------
      if(dualbasis) then
c
c  get eigenvectors & eigenvalues
c
c  the Fock operator is defined in a small basis set
c  but its matrix representation goes over the big one
c  It will be diagonalized WITHOUT allowing occupied
c  and virtuals orbitals to mix.
c
ckw07
c       call get_mix_eigen(bl,nblocks,bl(ictr),labels,thresh,
c    *                     ncf,ncf_sm,ncs,ncs_sm,nmo)
ckw07
        call get_mix_eigen2(bl,nblocks,bl(ictr),labels,thresh,
     *                     ncf,ncf_sm,ncs,ncs_sm,nmo,itype)
c
c       integral threshold can be changed on return if "big"
c       basis set overlap shows linear dependency
      else
        call matread('cano',np4,'evec_rhf')
        call matread('epsi',np4,'eval_rhf')
      endif
c-------------------------------------------------------
c get number of virtuals again (could be changed in get_mix_eigen2)
c
      lastvirt=igetival('nonredun')
c-------------------------------------------------------
c initialize the two-el. integral program
c
      iforwhat=5
      thint = thresh
      call para_jobinit(iforwhat)
      call ptwoint(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             0,      scftype,xintxx, nblocks,maxbuffer,
     *             maxlabels)
c
c  reserve memory for labels ..
c
        call getmem(maxlabels,labels)
c-------------------------------------------------------
      ncore=ncores(ncf,bl(iepsi),core)
      nfirst=ncore+1
      nlast=nmo
      nval=nmo-ncore
cc      nvirt = ncf-nmo
      nvirt = lastvirt-nmo
      if(nofr) then
        nfirst=1
        nlast=nmo
        nval=nmo
      end if
      npairs=nval*(nval+1)/2
      if (npairs.gt.maxpair) call nerror(3,'MP2 module',
     $   'Too many pairs for array allocator size',igranules,maxpair)
c---------------------------------------------------------------------
c check the maximum values of BOTH
c the occupied and virtual eigenvectors
c
      iocco=mataddr('cano')
      iepsi=mataddr('epsi')
      call check_coef(ncs,ncf,ncore,nmo,nval,bl(iocco))
c
c---------------------------------------------------------------------
c Construct "screening density" matrix to be used to estimate the
c importance of integrals. It can be based either on canonical or
c localized orbitals. The latter has the advantage that integral
c prescreening based on localizability is supposed to be more efficient.
c
c For ordinary one-basis set MP2 the "screening density" is made
c out of canonical or localized vectors taken from the disk as
c "evec_rhf" or "loca_rhf" .
c
c For dual-basis set the "screening density" is formed using current
c MOS content. In the input deck, the GUESS=READ statement MUST precede
c the MP2 or LMP2 keyword in order to ensure the "corresponding orbital"
c transformation (projection from small to big basis set) has been
c preformed.
c
c
c---------------------------------------------------------------------
c        check if Localization was performed
         np4=4
         call sudat(np4,'loca_rhf',ni)
c
         if(ni.gt.0) then
           loca_done=1
      write(iout,*)' Localized orbitals are used for integral screening'
          if(.not.dualbasis) then
            if(.not.restrt) 
     1      call DmxMakeL(dualbasis, nval, nmo, ncf, ncs,
     2                    np4,.false.,inx, iprnt, bl(idics))
          else
            call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
     *                    np4, inx, iprnt)
          endif
         else
            loca_done=0
      write(iout,*)' Occupied orbitals are used for integral screening'
          if(.not.dualbasis) then
            if(.not.restrt)
     1      call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
     2                    np4, inx, iprnt)
          else
            call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
     *                    np4, inx, iprnt)
          endif
         endif
c---------------------------------------------------------------------
c
c For dual basis set make screening density using CANONICAL orbitals
c
c     if(.not.dualbasis) then
c
c      check if Localization was performed
c      call sudat(np4,'loca_rhf',ni)
c      if(ni.gt.0) then
c        loca_done=1
c     write(iout,*)' Localized orbitals are used for integral screening'
c         if(.not.restrt)
c    1    call DmxMakeL(dualbasis, nval, nmo, ncf, ncs,
c    2                  np4,.false.,inx, iprnt, bl(idics))
c      else
c         loca_done=0
c         write(iout,*)
c    1 ' Canonical orbitals are used for integral screeing; results',
c    2 ' may be inaccurate for large molecules. Please use the LOCA',
c    3 ' option in SCF (say LOCA=PIPEK) to improve accuracy and',
c    4 ' performance for large systems'
c       call nerror(1,'cmp2',
c    1  'Please use localized orbitals in SCF for screening',0,0)
c         if(.not.restrt)
c    *       call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
c    *                     np4, inx, iprnt)
c      endif
c     else
c      write(iout,*)' Occupied orbitals are used for integral screening'
c        call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
c    *                 np4, inx, iprnt)
c     endif
c---------------------------------------------------------------------
c
c  establish orbital symmetry characters
      if(nsym.gt.0) then
        call getint(nval*nsym,iorbpair)
        call getmem(ncf,itmp)
        call OrbSymPair(nsym, ncf ,   bl(icano), bl(ifp) ,nfirst,
     1                  nlast,iprnt,  bl(itmp),  bl(iorbpair))
        call retmem(1)
      endif
      if(ifound(2).gt.0) then
        nfirst=iopv(1,2)
        nlast=iopv(2,2)
        nval=nlast-nfirst+1
      endif
c  print the data of the calculation early
      write(iout,'("  Number of contracted functions =",i8,/,
     1             "  Number of correlated orbitals  =",i8,/,
     2             "  Number of virtual orbitals     =",i8)')
     3                ncf,nval,lastvirt-nmo
      write(iout,*) ' '
      call f_lush(iout)
c.................................................
c check if there will be a split in integrals calculations
cc      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
      call matsub('occu','cano',nfirst,nlast)
      call matsub('virt','cano',nmo+1,lastvirt)
      ioccu=mataddr('occu')
c.................................................
      if(iprnt.gt.3) then
        call matprint ('occu',6)
        call matprint ('virt',6)
      end if
c .................................................
c -- determine length of direct access file records
c
      ints = 4              ! size of an integer in bytes
      lrec=(nval**2+4)*ints + nval**2
c
      if(iprnt.gt.1) then
        write(iout,2100) lrec
 2100 Format(' Record length of half-transformed DA File is ',I10,
     $       ' bytes')
      endif
c .................................................
c
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c---------------------------------------------------
c     if(dualbasis.and.itype.eq.1) then
c        update content of "basis sets" according to current needs
c        i.e. put small basis set into "basis 2" & "basis 4"
c
c        call update_basis(bl,ncs_sm)
c     endif
c---------------------------------------------------
c
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c
      ENonZero=zero
c
ckw..........
      call init_info('fock')
ckw..........
c
      call elapsec(tstart)
c
c ....................................................................
c
c -- send initial MP2 data to slaves
      call para_initsend
      call para_pack_int(nfirst,1)
      call para_pack_int(nlast,1)
      call para_pack_int(nval,1)
      call para_pack_int(lrec,1)
      call para_pack_int(smal,1)
      call para_pack_int(keep,1)
      call para_pack_int(restrt,1)
      call para_pack_int(Test,1)
      call para_pack_int(mp2only,1)
      call para_pack_int(bl(mpres_in_bg),ncs)
      call para_pack_real(thresh,1)
      call para_pack_real(bl(iepsi),ncf)
      call para_pack_real(bl(icano),ncf*ncf)
      call para_bcast_pack(TxMP2Dat)
c
c -- send scr. den., symm. info and filename separately
      call para_initsend
      call para_pack_real(bl(idics),ncs*ncs)
      If(nsym.gt.0)
     $  call para_pack_int(bl(iorbpair),nval*nsym)

c -- take care with file handling
c -- send each slave the base filenames
c -- they will make it unique based on their gid (1...nslv)
      call para_pack_string(filname1,256)
      call para_pack_int(len1,1)
      call para_pack_string(filname2,256)
      call para_pack_int(len2,1)
      If(.NOT.emp2only) Then
         call para_pack_string(filname3,256)
         call para_pack_int(len3,1)
         call para_pack_string(filname4,256)
         call para_pack_int(len4,1)
      End If
      call para_bcast_pack(TxMP2File)
c .....................................................................
c
      if(restrt) go to 50
c
c .....................................................................
c
c    F I R S T   H A L F - T R A N S F O R M A T I O N
c    -------------------------------------------------
c
      do ics=ncs,1,-1
         do kcs=ics,1,-1
c -- of each (ics,kcs) pair, calculate only the last one
c -- (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
            if(LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c
c -- send ics,kcs pair to slave
            call para_recv(imsg,islave,TxDftReq)
            call para_initsend
            call para_pack_int(ics,1)
            call para_pack_int(kcs,1)
            call para_send_pack(islave,TxMP2Assign)
         end do
      end do
c
c -- all done  stop all slaves
      call elapsec(timo)
      call getrval('imbal',tlost)
      Do ics=1,nslv
         call para_recv(imsg,islave,TxDftReq)
         call para_initsend
         call para_pack_int(0,1)
         call para_pack_int(0,1)
         call para_send_pack(islave,TxMP2Assign)
         call elapsec(timb)
         tlost=tlost+(ics-1)*(timb-timo)
         timo=timb
      EndDo
      call setrval('imbal',tlost)
      call elapsec(tend)
c
      If (iprnt.gt.0) Then
         telap = (tend-tstart)/60.0d0
         write(iout,1000) telap
 1000    Format(' End of First Half-transformation    elapsed time = ',
     $           F9.2,' min')
         call f_lush(iout)
      EndIf
c
c .....................................................
      If (Test) Then
         write(iout,65)
   65    format('  == Calculation Halted at User Request ==')
         STOP
       EndIf
c .....................................................
c
 50   CONTINUE
c
c .....................................................................
c
c    B I N   S O R T  &  S E C O N D   T R A N S F O R M A T I O N
c    -------------------------------------------------------------
c
c    At this point each slave has its own partially transformed integral file,
c    with all w,v (AOs) for each I,J (MOs). This needs to be resorted to give
c    all I,J for each w,v. We use a standard Yoshimine bin sort to reorder the
c    indices; however, as each bin is filled it is NOT written out to a direct
c    access file on that slave, but is instead sent to the slave that will
c    ultimately transform that (w,v) pair. 
c
C    This is done by each slave polling for bins sent to it and writing them
c    to a direct access file. At the end of the bin sort, each slave should
c    only those (w,v) pairs that it will subsequently transform.
c
c
      call elapsec(tstart)
c
c -- check if the slaves are still OK
      call para_check
c ..........................................................................
c -- determine memory available for the bins. They should be as big
c -- as possible but leave space for other allocations
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memvirt=memtot
      memoccu=lastaddr-ioffset
      memavail=memvirt-memoccu
c
c -- memory for the bins is the minimum of the available virtual and
c -- real memory, minus 1MW for miscellaneous
      memtrans=3*nval*(nval+1)/2+nvirt**2+ncf*nvirt
      mem=memavail-memtrans-1000000
c -- lbin is in units of 9-bytes (2 indices           - integer*2
c                                 compressed integral - integer*4
c                                 precision overflow  - integer*1)
      lbin=(8*mem)/(9*npairs)
c     lbin=lbin/2       ! divide as less memory on slaves (TJ ?)
      lbin=lbin*0.9d0   ! divide as less memory on slaves (?)
      If (lbin.gt.ncf*(ncf+1)) lbin = ncf*(ncf+1)
cPP -- limit lbin in case of large basis and large memory (?)
      maxlbin = 524287
      If (lbin.gt.maxlbin) lbin=maxlbin
c
      If (lbin.lt.ncf*(ncf+1).AND.lbin.lt.100) Then
         call nerror(5,'MP2 module',
     1        'memory available for bins leads to bin size < 100',
     2         mem,lbin)
      EndIf
c
c Get amount of available memory on slaves:
      islavemem=2147483646                   ! max positive integer, FRAGILE
      do i=1,nslv
         call para_recv(islavemem1,islave,TxMP2S6)
         if (islavemem1.lt.islavemem) islavemem=islavemem1
      enddo
c
c -- too large bins cause paging in PVM version during bin sort
c -- half actual memory so that bins are smaller
      islavemem = islavemem/2
c     write(iout,'('' Slave memory available:'',I12)'),islavemem
      islavemem=islavemem-nvirt*nvirt ! for fully transformed matrix
c
c -- can we write all the bin sort files to disk at one time?
c -- if not, will need multiple passes
c    [disk space remaining is DiskMax - size of half-transformed
c     file calculated as
c       (total number of indices * size of record)/nslaves]
c
      Usage = dble(ncf*(ncf+1)/2)*lrec/(nslv*gigabyte)  ! current disk usage (GB)
      DiskLeft = DiskMax - Usage
c                9 because of bin size with indices. +2 because ij pars
c      may not be evenly distributed over slaves (overestimation done "on
c       purpose").
      DiskNeed = (9+2)*npairs*dble(ncf*(ncf+1))/(nslv*gigabyte) ! disk required
cc      write(6,*) ' Disk usage: ',usage,' left:',Diskleft,
cc     &           ' need:',diskneed
c
      If (DiskLeft.LT.2.0d0) Then
         call nerror(3,'MP2 module',
     $   'Not Enough Disk Space left for bin sort',0,0)
      EndIf
c
      If (DiskNeed.GT.DiskLeft) Then
         NPass = NINT(DiskNeed/DiskLeft) + 1
      Else
         NPass = 1
      EndIf
c
      if (iprnt.gt.1) then
         write(iout,'(a,a,i4,a)') ' Second Half-transformation will be',
     $                ' done in ',npass,' pass(es)'
         call f_lush(iout)
      endif
CTJ granules of bins will be written
c
      if (igran_calculate) then
         igranulesize=islavemem/(lbin*9/8+1+ncf*ncf)
         do
            igranules=npairs/igranulesize
            if (mod(npairs,igranulesize).ne.0) igranules=igranules+1
            if (igranules.ge.nslv*NPass*2) exit ! 2 granules/slave/pass
            if (igranulesize.lt.1) call nerror(1,'MP2 module',
     $                        'Error in granule determination',0,0)
            if (igranulesize.eq.1) exit
            igranulesize=igranulesize-1
         enddo
      endif
      do
         lrec =(lbin+8*lbin+4)*igranulesize+8  ! record length in bytes
c Limit record size to 2 MB. This pevents overloading the binrdr's
c by large amount of large blocks passed via PVM, which may next cause RAM
c memory overflow or swapping and efficiency loss.
         if (lrec.gt.2097152) then
            if (igranulesize.eq.1) exit
            igranulesize=igranulesize-1
         else
            exit
         endif
         igranules=npairs/igranulesize
         if (mod(npairs,igranulesize).ne.0) igranules=igranules+1
      enddo
c     print *,'igranules:    ',igranules
c     print *,'igranulesize: ',igranulesize
c     call f_lush(iout)
      if (igranules.gt.maxgran) call nerror(3,'MP2 module',
     $   'Too many granules for array allocator size',igranules,maxgran)
      call f_lush(iout)
c build maps: granule -> ipairstart,ipairstop
c             ipair   -> granule
      istart=1
      do i=1,igranules
         istop=istart+igranulesize-1
         if (istop.gt.npairs) istop=npairs
         if (istart.gt.npairs) call nerror(3,'MP2 module',
     $   'Error in granule maps determination ',istart,npairs)
         igran2pair(1,i)=istart
         igran2pair(2,i)=istop
         do j=istart,istop
            ipair2gran(j)=i
         enddo
         istart=istop+1
      enddo
c
      if (iprnt.gt.3) then
         do i=1,igranules
            write(iout,'(''Granule: '',i0,'' start: '',i0,
     &                 '' stop: '',i0)')
     *                 i,igran2pair(1,i),igran2pair(2,i)
         enddo
         call f_lush(iout)
      endif
c -- determine number of bins and bin ownership list
c -- i.e., which pairs will be handled on which slaves
      ics = igranules/nslv          ! initial # of pairs per slave
      left = igranules - ics*nslv   ! # of pairs left
      nbin = ics                    ! maximum amount of gran. per slave
      if (left.gt.0) nbin = nbin+1
c
c -- simply distribute granules evenly to each slave
      do i=1,nslv
         binlist(i) = ics
      enddo
      do i=1,left
         binlist(i) = binlist(i)+1
         if (binlist(i).eq.0) write(iout,*)
     *       'WARNING!!! Some slave will have empty ij list!'
      enddo
cc      write(6,*) ' MASTER: no. of bins per slave:'
cc      write(6,*) (binlist(i),i=1,nslv)
      call f_lush(iout)
c
c -- now convert to actual bins
      do i=2,nslv
         binlist(i) = binlist(i-1) + binlist(i)
      enddo
cc      write(6,*) ' MASTER: final binlist is:'
cc      write(6,*) (binlist(i),i=1,nslv)
cc      write(6,*) ' maximum bin size:',nbin
cc      write(6,*) ' MASTER:  slave IDs and hosts are:'
cc      do i=1,nslv
cc      write(6,*) i,'  ',host(i)(1:20),'  ',slave1(i)
cc      enddo
cc      call f_lush(6)
c
c            overflow  bins    numij               prevrec
      lrec =(lbin     +8*lbin   +4  )*igranulesize +8  ! record length in bytes
c
      if (iprnt.gt.1) then
         write(iout,2200) lrec
         call f_lush(iout)
 2200 Format(' Record length of 2nd-transformed DA File is ',I10,
     $       ' bytes')
      endif
c .........................................................................
c
c
c -- send to the original slaves the number of passes
c    and the initialisation data for second half-transformation
c
      call para_initsend
      call para_pack_int(NPass,1)
      call para_pack_int(lbin,1)
      call para_pack_int(igranulesize,1)
      call para_pack_int(igranules,1)
      call para_pack_int(npairs,1)
      call para_pack_int(igran2pair,2*igranules)
      call para_pack_int(ipair2gran,npairs)
      call para_bcast_pack(TxBlockAssign)
c
c disable direct routing during binsort (applies only to PVM version)
c (some small MP2 jobs show random errors if direct routing
c  is enabled during the parallel binsort in the PVM version)
c This is now done only for small jobs ( less than 200 basis functions)
c because for large jobs the nodirect route might cause the PVM
c deamon to allocate too much memory and the job will either die or
c start to swap.
c
      if(ncf.lt.200)  call para_nodirectroute
c
      If (NPass.GT.1) Then
         do i=1,nslv
            binlist(i) = NINT(dble(binlist(i))/NPass)
            if (i.gt.1) then
               iper_slave=binlist(i)-binlist(i-1)
            else
               iper_slave=binlist(i)
            endif
            if (iper_slave.eq.0) write(iout,*)
     *       'WARNING!!! Some slave will have empty ij list!'
         enddo
         lastbin = binlist(nslv)
      EndIf
      call f_lush(iout)
c
      tbin = zero
c
      DO 75 IPass=1,NPass
         call elapsec(t1)               ! start of bin sort
c
c -- set binlist array
         If (IPass.GT.1) Then
            do i=1,nslv
               binlist(i) = binlist(i) + lastbin
            enddo
         EndIf
         if (nslv.gt.1) then
            if (binlist(nslv-1).gt.igranules) 
     $             call nerror(1,'MP2 master',
     $             'Internal binlist error ',binlist(nslv-1),igranules)
         endif
         If (IPass.EQ.NPass) binlist(nslv) = igranules
cc      write(6,*) ' MASTER:  IPass is ',ipass,' binlist array is:'
cc      do i=1,nslv
cc      write(6,*) i,'  ',binlist(i)
cc      enddo
cc      call f_lush(6)
c
c -- send initialization data and bin ownership to slaves
c -- npairs can be comp-d from binlist      
         call para_initsend
         call para_pack_int(binlist,nslv)
         call para_bcast_pack(TxBinDat)
c
c -- now wait until get message from EACH slave that
c -- the bin sort has been completed
         Do i=1,nslv
            call para_recv_pack(ifrom,TxBinS1)
         EndDO
c
         call elapsec(t2)               ! end of bin sort
         tbin = tbin + (t2-t1)          ! accumulated bin sort time
c
c -- now wait until get message from EACH original slave that
c -- this pass of the second half-transformation has completed
         Do i=1,nslv
            call para_recv(imsg,islave,TxBlockReq)
         EndDo
         if(iprnt.gt.2) write(iout,'(a,i4,a)')' MASTER: Finished pass ',
     $                                     ipass
         call f_lush(iout)
 75   CONTINUE
      call f_lush(iout)
c
c -- check if all the slaves are OK
      call para_check
c
c -- now sum up total MP2 energy and wavefunction norm from each slave
      emp2 = zero
      tnorm = zero
      emps=zero
      empt=zero
      empanti=zero
      emppar=zero
      call para_reduce(emp2,   1,TxMP2S1)
      call para_reduce(tnorm,  1,TxMP2S2)
      call para_reduce(emps,   1,TxMP2S3)
      call para_reduce(empt,   1,TxMP2S4)
      call para_reduce(emppar, 1,TxMP2S5)
      call para_reduce(empanti,1,TxMP2S6)
c
      call elapsec(tend)
      If(iprnt.gt.0) Then
        telap = (tend-tstart)/60.0d0
        tbin = tbin/60.0d0
        write(iout,1200) telap,tbin
 1200   Format(' End of Second Half-Transformation   elapsed time = ',
     $           F9.2,' min',/,
     $         ' time spent in bin sort = ',F9.2,' min')
      EndIf
c
c -- timing stuff
      call para_next(-1)
c
c restore direct routing (applies only to PVM version)
      if(ncf.lt.200) call para_directroute
c
c Print out the final results :
c
      tnorm = tnorm+one
      tnorm=sqrt(tnorm)
c
      esing_exc = zero
      if(dualbasis) call getrval('esing_exc',esing_exc)
      edual = erhf + esing_exc
c
      write(iout,*) ' '
      write(iout,*) '-------------------------------------------------'
      write(iout,*) '                Final SCF/MP2 results '
      write(iout,*) ' '
      write(iout,3001) erhf
      if(dualbasis) then
         write(iout,3002) esing_exc
         write(iout,3022) edual
         write(icond,3022) edual
      write(iout,*) '-------------------------------------------------'
      endif
      write(iout,3003) emp2
      write(iout,*) '-------------------------------------------------'
      write(iout,3004) edual + emp2
      write(iout,3005) tnorm
      write(iout,*) '-------------------------------------------------'
c-------------
c -- Scaled MP2 - see S. Grimme, J. Chem. Phys. 118 (2003) 9095
c
      write(iout,*)'                    Scaled MP2    '
      write(iout,*)'      [See S.Grimme, J.Chem.Phys. 118 (2003) 9095]'
      write(iout,*) ' '
      write(iout,3033) emps,empt
      write(iout,3034) empanti,emppar
      empscal = 1.2d0*empanti+emppar/3.0d0
      if(iscs.eq.2) emp2scal=p1*empanti+p2*emppar
      write(iout,3035) empscal
      write(iout,3036) erhf+empscal
      if(iscs.eq.2) then
        write(iout,3037)emp2scal,p1,p2
        write(iout,3038) erhf+emp2scal
      endif
      write(iout,*) '-------------------------------------------------'
      write(iout,*) ' '
c-------------
c
 3001 format('        Total SCF energy         =',f15.8)
 3002 format('        Single excit. contrib.   =',f15.8)
 3022 format('        Corrected SCF energy     =',f15.8)
 3003 format('        MP2 correlation energy   =',f15.8)
c
 3033 format('        Singlet part:',f15.8,'  Triplet part:',f15.8)
 3034 format('        Antiparallel:',f15.8,'  Parallel:    ',f15.8)
 3035 format('        SCS-MP2 energy           =',f15.8)
 3036 format('        Total SCS-MP2 energy     =',f15.8)
 3037 format('        Scaled MP2 energy        =',f15.8,' p1=',f5.3,
     $        'p2=',f5.3)
 3038 format('        Total Scaled MP2 energy  =',f15.8)
c
 3004 format('        Total MP2 energy         =',f15.8)
 3005 format('        Norm of the wavefunction =',f15.8)
c
      call f_lush(6)
c
c -- set wavefunction type
      wvfnc = 'RMP2'
      If(dualbasis) wvfnc = 'RMP2-dual'
      If(scs) wvfnc = 'RMP2-scs'
      If(scs.AND.dualbasis) wvfnc = 'RMP2-dual-scs'
c
c -- write final MP2 energy to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      If(iscs.eq.0) Then
        Call wrcntrl(IUnit,7,'$energy',2,idum,emp2+edual,wvfnc)
      Else If(iscs.eq.1) Then
        Call wrcntrl(IUnit,7,'$energy',2,idum,erhf+empscal,wvfnc)
      Else
        Call wrcntrl(IUnit,7,'$energy',2,idum,erhf+emp2scal,wvfnc)
      Endif
      If(dualbasis) Call wrcntrl(IUnit,6,'$edual',2,idum,edual,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
c     call retimark
      call matremark
      call retmark
c---------------------------------------------------------------------
      call setival('ncore',ncore)
      call setrval('thresh',thresh)
      RETURN
      END
c .....................................................................
c
      subroutine Para_BinSort( 
     1   nval,      npairs,    ndisk,       ndisk2,    lbin,
     2   nbin,      ijfirst,   ijlast,      binlist,   ibins,
     3   ibin1,     i1pair,    int1,        icounter,  keep,
     4   igran2pair,ipair2gran,igranulesize,igranules, ngran,
     5   imod,      sendb)

      use newpara
      implicit real*8 (a-h,o-z)
c
c
c  This routine sorts the half-transformed integrals
c
c  arguments
c
c  nval    -  number of occupied MOs to be transformed
c             (usually the valence MOs)
c  npairs  -  number of occupied orbital pairs
c             (nval*(nval+1)/2)
c  ndisk   -  unit number for half-transformed integrals
c  ndisk2  -  unit number for the sorted bins
c  lbin    -  length of one bin
c  nbin:      bin location array (which bin granule 
c             on which record in binsort file)
c  ijfirst -  index of first I,J pair to be sorted
c  ijlast  -  index of last  I,J pair to be sorted
c             (the pair number IJ is I*(I-1)/2 + J
c  binlist -  list of which bins granules to send to which slaves
c  ibins   -  integer*4 storage for the bins
c  ibin1   -  integer*1 storage for precision overflow
c             (how many times 2**31 divides the integer 
c              representation of the integral)
c  i1pair  -  integer*4 storage for the exchange matrix to be sorted
c  int1    -  integer*1 storage for precision overflow
c             (how many times 2**31 divides the integer 
c              representation of the integral)
c  icounter   integer*4, shows degree of occupancy of the bins
c  keep    -  logical flag for keeping half-transformed file 
c             after sort
c  igran2pair, ipair2gran - conversion arrays, granule to pair 
c              and pair to granule
c  igranules - total number of granules
c  igranulesize - size of present granule
c  ngran:     number of bin granules in this pass
c  imod      granule number modifier for current slave
c  sendb     buffer for non-blocking send (used only by the MPI version)
c
      integer*1 int1(nval,nval),ibin1(lbin,npairs)
      integer*2 iindex(2),jindex(2)
      integer*4 ibins(2*lbin,npairs),i1pair(nval,nval)
      integer*4 icounter(npairs)
      integer binlist(nslv),nbin(ngran)
      integer igran2pair(2,igranules),ipair2gran(npairs)
      integer*4 longint,longint1
      equivalence(longint,iindex(1)),(longint1,jindex(1))
      logical keep
c
      do i=1,npairs
         icounter(i)=0  ! cannot use izeroit because integer*4
      enddo
      irec=0
      nrec=0            ! record no in binsort file
      islv=nslv         ! number of unfinished slaves
      lbin2 = 2*lbin
c
 50   irec=irec+1
      read(ndisk,rec=irec) mu,lam,int1,i1pair
c
c  internal end of file: mu=lam=-1
c
      if(mu.eq.-1.and.lam.eq.-1) go to 100
c
      iindex(1)=mu
      iindex(2)=lam
      jindex(1)=lam
      jindex(2)=mu
      ij=0

      outer: do i=1,nval
         do j=1,i
            ij=ij+1
            if(ij.lt.ijfirst) cycle
            if(ij.gt.ijlast) exit outer
c
c -- test on integral  (don't store if zero)
c
            If (i1pair(i,j).ne.0.OR.int1(i,j).ne.0) Then
c
c   write bunch of bins even if some of them are not full.
c
               if (icounter(ij).ge.lbin2) then
                  igran=ipair2gran(ij)
                  istart= igran2pair(1,igran)
                  istop = igran2pair(2,igran)
                  isize=istop-istart+1
                                         !  compute destination slave
                  idest=-1
                  do n=1,NSLV
                     if (igran.le.binlist(n)) then
                        idest=n
                        exit
                     endif
                  enddo
                  if (idest.lt.0) call nerror(1,'para_binsort',
     $               'Know not where to send bin ',igran,0)

                          ! check if there are any bins to be received

                  call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                               ndisk2 ,nbin, nrec, islv, .false.)

                  if (MY_GID.eq.idest) then ! destination is this slave
                     call para_write_bin(ibins(1,istart),
     $                                   ibin1(1,istart),
     $                                   icounter(istart),lbin,
     &                                   igranulesize,nbin,ngran,imod,
     $                                   igran,ndisk2,nrec)

                  else               ! destination is another slave
                     call para_sendbin(ibins(1,istart),ibin1(1,istart),
     &                                 lbin,icounter(istart),isize,
     &                                 igran,idest,sendb,igranulesize, 
     &                                 ngran,imod,ndisk2,nbin,nrec,islv)
                  endif
c
                  do iii=istart,istop
                     icounter(iii)=0
                  enddo
               end if
c
               icounter(ij)=icounter(ij)+2
               ibins(icounter(ij)-1,ij)=i1pair(i,j)
               ibins(icounter(ij),ij)=longint
               ibin1(icounter(ij)/2,ij)=int1(i,j)
            EndIf
c
            IF (mu.ne.lam) THEN
c
c -- test on integral  (don't store if zero)
               If (i1pair(j,i).ne.0.OR.int1(j,i).ne.0) Then
c   write bunch of bins even if some of them are not full.
                  if (icounter(ij).ge.lbin2) then
                     igran=ipair2gran(ij)
                     istart= igran2pair(1,igran)
                     istop = igran2pair(2,igran)
                     isize=istop-istart+1
                                         !  compute destination slave
                     idest=-1
                     do n=1,NSLV
                        if (igran.le.binlist(n)) then
                           idest=n
                           exit
                        endif
                     enddo
                     if (idest.lt.0) call nerror(1,'para_binsort',
     $                  'Know not where to send bin ',igran,0)

                          ! check if there are any bins to be received

                     call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                                  ndisk2,nbin,nrec, islv, .false.)

                     if (MY_GID.eq.idest) then ! destination is this slave
                        call para_write_bin(ibins(1,istart),
     &                                      ibin1(1,istart),
     &                                      icounter(istart),lbin,
     $                                      igranulesize,nbin,ngran,
     &                                      imod,igran,ndisk2,nrec)

                     else               ! destination is another slave
                        call para_sendbin(ibins(1,istart),
     &                                    ibin1(1,istart),lbin,
     &                                    icounter(istart),isize,igran,
     &                                    idest,sendb,igranulesize,
     &                                    ngran,imod,ndisk2,nbin,nrec,
     &                                    islv)
                     endif

                     do iii=istart,istop
                        icounter(iii)=0
                     enddo
                  end if
                  icounter(ij)=icounter(ij)+2
                  ibins(icounter(ij)-1,ij)=i1pair(j,i)
                  ibins(icounter(ij),ij)=longint1
                  ibin1(icounter(ij)/2,ij)=int1(j,i)
               EndIf
            ENDIF         !  mu=lam condition
         end do          !  j loop
c
c check for incoming bins
c
         call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                     ndisk2,nbin, nrec, islv, .false.)
      end do outer   !  i loop
c     
      go to 50
  100 continue
c  end-of-file
      If (keep) Then
         CLOSE (Unit=NDisk,STATUS='KEEP')
      Else
         CLOSE (Unit=NDisk,STATUS='DELETE')
      EndIf
c  write out the content of the bins at the end
      do ij=1,npairs
         if (icounter(ij).gt.0) then
c  zero out the rest
c the number of cells occupied is 2*icounter(ij)
            do j=icounter(ij)/2+1,lbin
               ibins(2*j-1,ij)=0
               ibins(2*j,ij)=0
               ibin1(j,ij)=0
            end do
         end if
      end do
c   write bunch of bins even if some of them are not full.
      do ij=1,npairs
         if (icounter(ij).gt.0) then
            igran=ipair2gran(ij)
            istart= igran2pair(1,igran)
            istop = igran2pair(2,igran)
            isize=istop-istart+1
                                         !  compute destination slave
            idest=-1
            do n=1,NSLV
               if (igran.le.binlist(n)) then
                  idest=n
                  exit
               endif
            enddo
            if (idest.lt.0) call nerror(1,'para_binsort',
     $         'Know not where to send bin ',igran,0)

            if (MY_GID.eq.idest) then ! destination is this slave
               call para_write_bin(ibins(1,istart),ibin1(1,istart),
     &                             icounter(istart),lbin,igranulesize,
     &                             nbin,ngran,imod,igran,ndisk2,nrec)
            else               ! destination is another slave
               call para_sendbin(ibins(1,istart),ibin1(1,istart),lbin,
     &                           icounter(istart),isize,igran,idest,
     &                           sendb,igranulesize,ngran,imod,ndisk2,
     $                           nbin,nrec,islv)
            endif

                      ! check if there are any bins to be received

            call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                         ndisk2 ,nbin, nrec, islv, .false.)

            do iii=istart,istop
               icounter(iii)=0
            enddo
         endif
         call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                     ndisk2,nbin, nrec, islv, .false.)
      enddo
c -- tell others that bin sort is complete
      Do i=1,nslv
         if (i.eq.MY_GID) then
            call para_sendbin(ibins,       ibin1,  lbin,  icounter,
     $                        -1,          -1,     0,     sendb,
     $                        igranulesize,ngran,  imod,  ndisk2,
     $                        nbin,        nrec,   islv) ! to master
            islv=islv-1                       ! to self
         else   
            call para_sendbin(ibins,       ibin1,  lbin,  icounter,
     $                        -1,          -1,     i,     sendb,
     $                        igranulesize,ngran,  imod,  ndisk2,
     $                        nbin,        nrec,   islv)
         endif
         call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                  ndisk2 ,nbin, nrec, islv, .false.)
      EndDO
c
      if (islv.gt.0) then
         call para_recv_bin(lbin,  igranulesize,ngran,imod,
     $                      ndisk2,nbin, nrec, islv, .true.)
      endif
c
      return
      end
c .....................................................................
c
      SUBROUTINE GetSymPair(nsym, nfirst, nmo, iorbpair,
     $                      I,    J,      ipair, jpair)
      IMPLICIT INTEGER(A-Z)
      Dimension iorbpair(nsym,nfirst:nmo),ipair(7),jpair(7)
c
      do isym=1,nsym
      ipair(isym) = iorbpair(isym,I)
      jpair(isym) = iorbpair(isym,J)
      enddo
c
      return
      end
c .....................................................................
c
      SUBROUTINE ReadBin(ncf,    lbin,   ndisk,  icount, ngran,
     $                   NBin,   thresh, nsym,   nfirst, NAlpha,
     $                   iorbpair, ifp,    IBins,  IBin1,  xmat,
     $                   igran2pair,ipair2gran,igranule,igranulesize,
     $                   ncore)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Does the second-half transformation and forms the contribution
C  to the total MP2 energy from each (w,v) AO index pair
C
C  ARGUMENTS
C
C  ncf     -  number of basis functions (AOs)
C  lbin    -  number of integrals in each record on bin sort file
C  ndisk   -  unit number of first bin sort file
C  icount  -  location of this record in bin location array
C  ngran   -  number of granules we are dealing with
C  NBin    -  bin location array
C             (location of last record on direct access file for this gran)
C  thresh  -  integral threshold
C  nsym    -  number of (Abelian) symmetry operations
C  nfirst  -  the first correlating orbital
C  NAlpha  -  number of occ orbitals
C  iorbpair-  symm pairs of orbitals
C  ifp     -  symm pairs of bf
C  ifp     -  basis function symmetry pairs
C  IBins   -  integrals & indices (stored as compact integers)
C  IBin1   -  integer*1 storage for precision overflow
C           (how many times 2**31 divides integer representation of integral)
C  xmat    -  on exit, exchange matrix for current IJ index
C  igran2pair, ipair2gran - conversion arrays, granule to pair and pair
C                           to granule
C  igranule - the granule number, which we are completing now
C  igranulesize - size of present granule
C
C
      Dimension xmat(ncf,ncf,igranulesize), NBin(ngran)
      Integer ifp(7,ncf),ipair(7),jpair(7),iorbpair(nsym,*)
      Integer*2 iindex(2)
      Integer*1 IBin1(lbin,igranulesize)
      integer*4 IBins(2,lbin,igranulesize)
      integer*4 icounter(igranulesize)
      integer   igran2pair(2,*),ipair2gran(*)
      integer*4 longint
      Equivalence (longint,iindex(1))
      Parameter (dblmax=2 147 483 648.0d0)
      integer ib1j
C
C
C  read back bins, extract one exchange matrix, and transform it
C
      ijstart=igran2pair(1,igranule)
      ijstop =igran2pair(2,igranule)
      isize=ijstop-ijstart+1
      CALL ZeroIT(xmat,ncf*ncf*isize)
      dblcmp = dblmax*thresh
      KUnit = ndisk
      irec = NBin(icount)         ! location of last record
c     print *,'icount,irec,ngran: ',icount,irec,igranule,ijstart,ijstop
c
 50   CONTINUE
      If(irec.EQ.0) GO TO 45
      READ(UNIT=KUnit,rec=irec) irec,icounter,IBin1,IBins
      ij0=0
      do ij=ijstart,ijstop
      ij0=ij0+1
      numij=icounter(ij0)
      call get_ij_half(ij,ip,jp)
      ip = ip + ncore
      jp = jp + ncore
      If(nsym.gt.0)
     *         call GetSymPair(nsym, nfirst, NAlpha, iorbpair,
     *                         ip,   jp,     ipair,  jpair)
c
      DO 40 J=1,numij/2
      longint = IBins(2,J,ij0)
      mu = iindex(1)
      lam = iindex(2)
      If(mu.EQ.0.OR.lam.EQ.0) cycle     ! do we still need this?   JB
c
c -- decompress the integral -----------------------
cc      call reconstruct(ibins(1,j,ij0),ibin1(j,ij0),xx)
cc      xx = xx*thresh
      If(IBin1(J,ij0).eq.0) Then
        xx = IBins(1,J,ij0)*thresh
      Else If(IBin1(J,ij0).gt.0) Then
        xx = IBins(1,J,ij0)*thresh
        xx = xx + SIGN(IBin1(J,ij0)*dblcmp,xx)
      Else
        xx = IBins(1,J,ij0)*thresh*10.0d0**(-IBin1(J,ij0))
      EndIf
c --------------------------------------------------
      xmat(mu,lam,ij0) = xx
      DO 30 isym=1,nsym
        xxx = xx
        mu1 = ifp(isym,mu)
        mu2 = ABs(mu1)
        If(mu1.LT.0) xxx = -xxx
        lam1 = ifp(isym,lam)
        lam2 = ABs(lam1)
        If(mu2.EQ.mu.AND.lam2.EQ.lam) cycle
        If(lam1.LT.0) xxx = -xxx
        If(ipair(isym).LT.0) xxx = -xxx
        If(jpair(isym).LT.0) xxx = -xxx
        xmat(mu2,lam2,ij0) = xxx
 30   CONTINUE
 40   CONTINUE
      enddo
c
      If(irec.NE.0) GO TO 50       ! get another granule
C
C  last record for bin file
C
 45   CONTINUE
C
      RETURN
      END
c .....................................................................
c
      subroutine PSaveTij(ncf,    nmo,    nval,   nvir,   i,
     1                    j,      ncore,  e,      thresh, ibins,
     2                    int1,   x,      ndisk,  ij)
      implicit real*8 (a-h,o-z)
C
C  This subroutine starts with the internal exchange matrix Kij in
C  (virtual) MO basis, generates the amplitudes Tij by dividing by
C  the usual energy denominator (stored in X(nvirt,nvirt) or 'xmos'
C
C  xmos is scaled by the inverse of the integral threshold
C  and kept in an integer matrix intres(nvir,nvir) which is written to
C  external storage on ndisk, one record nvir*nvir long for each pair.
C
C  Transformation of this matrix to AO basis is now carried out in
C  mp2grad since storage in MO basis takes less space than AO basis
C  and the A1-terms are best calculated using Tij in MO basis.
C
C  Arguments (INPUT):
C  ncf:       # of AO basis functions
C  nmo        # of occupied orbital
C  nvir:      # of virtual orbitals
C  i,j:       occupied orbitals defining the pair(i,j)
C  e:         orbital energies
C  thresh:    integral threshold
C  intres:    integer matrix for residuum Tij in AO basis
C  int1:      integer*1 storage for precision overflow or threshold reduction
C  x:         matrix with Tij in virtual basis
C  ndisk:     unit on which the final matrix is written
C  ij:        record number
C
C  NOTE  x is the same as  'xmos'
C
C         Svein Saebo  Fayetteville Ar, and Starkville, MS May 2002
C
C  Modified from original SS version to save I,J indices to disk
C
      dimension e(*)
      dimension x(nvir,nvir)
      integer*1 int1(nvir,nvir),i1
      integer*4 ibins(nvir,nvir)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c
      end1=e(i)+e(j)
      do ib=1,nvir
      end2=end1-e(nmo+ib)
        do ia=1,nvir
        x(ia,ib)=x(ia,ib)/(end2-e(nmo+ia))
        end do
      end do
C  check magnitude and convert to integers
      call matscal('xmos',one/thresh)
      do ib=1,nvir
      do ia=1,nvir
      val = x(ia,ib)
      If(abs(val).ge.dblmax) Then
        b = val*dblinv
        if(abs(b).ge.d1max) then
          dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
          dfac = LOG10(dfac)
          i1 = -NINT(dfac+0.5d0)
          val = val*10.0d0**i1
          int1(ia,ib) = i1
          ibins(ia,ib) = val
cc          write(6,*) ' Threshold Tij - ia:',ia,' ib:',ib,' i1:',i1
        else
          i1 = abs(b)
          b = val - SIGN(i1*dblmax,val)
          int1(ia,ib) = i1
          ibins(ia,ib) = b
        endif
      Else
        int1(ia,ib) = 0
        ibins(ia,ib) = val
      EndIf
      enddo
      enddo
C  write final matrix to disk unit ndisk
      write(ndisk,rec=ij) i-ncore,j-ncore,int1,ibins
c
      return
      end
c .....................................................................
c
      subroutine PSaveKov(ncf,    nval,   nvir,   x,      ibins,
     1                    int1,   i,      j,      ncore,  thresh,
     2                    length, ndisk,  nrec)
      implicit real*8(a-h,o-z)
C
C  Converts <Kij> to AO basis and saves to local direct access file.
C  Modified from original SS version to save I,J indices to disk
C
      dimension x(*)
      integer*1 int1(length,2),i1
      integer*4 ibins(length,2)
      parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c
      call matpose('xmat')
      call matdef('tss02','r',nvir,ncf)
      call matmmult('tvir','xmat','tss02')
      if(ncore.gt.0) then
        call matmmult('tss02','occa','kijvo')
      else
        call matmmult('tss02','occu','kijvo')
      endif
C  check magnitude and convert to integers
      call matscal('kijvo',one/thresh)
      do imv=1,length
      val = x(imv)
      If(abs(val).ge.dblmax) Then
        b = val*dblinv
        if(abs(b).ge.d1max) then
          dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
          dfac = LOG10(dfac)
          i1 = -NINT(dfac+0.5d0)
          val = val*10.0d0**i1
          int1(imv,2) = i1
          ibins(imv,2) = val
cc          write(6,*) ' Threshold Kij - imv:',imv,' i1:',i1
        else
          i1 = abs(b)
          b = val - SIGN(i1*dblmax,val)
          int1(imv,2) = i1
          ibins(imv,2) = b
        endif
      Else
        int1(imv,2) = 0
        ibins(imv,2) = val
      EndIf
      enddo
cc
      call matpose('xmat')
      call matmmult('tvir','xmat','tss02')
      if(ncore.gt.0) then
        call matmmult('tss02','occa','kijvo')
      else
        call matmmult('tss02','occu','kijvo')
      endif
C  check magnitude and convert to integers
      call matscal('kijvo',one/thresh)
      do imv=1,length
      val = x(imv)
      If(abs(val).ge.dblmax) Then
        b = val*dblinv
        if(abs(b).ge.d1max) then
          dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
          dfac = LOG10(dfac)
          i1 = -NINT(dfac+0.5d0)
          val = val*10.0d0**i1
          int1(imv,1) = i1
          ibins(imv,1) = val
cc          write(6,*) ' Threshold Kij - imv:',imv,' i1:',i1
        else
          i1 = abs(b)
          b = val - SIGN(i1*dblmax,val)
          int1(imv,1) = i1
          ibins(imv,1) = b
        endif
      Else
        int1(imv,1) = 0
        ibins(imv,1) = val
      EndIf
      enddo
      call matrem('tss02')
C  finally write to disk
      write(ndisk,rec=nrec) i-ncore,j-ncore,int1,ibins
      end
c=================================================
c
      subroutine do_mp2
c
c=================================================

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      character*256 filename,filname1,filname2,filname3,filname4,
     $              Restart_File,msgpass
      character ch3*3,scftype*11
      Logical rhf,smal,keep,restrt,keep1,Test,exst,alphabeta
      dimension xintxx(9)
ckw
      dimension ICSpass(2,28), KCSpass(2,28) !28 is 28 comp.of cart. I-function
ckw

c     common /big/ bl(300)
c     common /intbl/ maxsh,inx(100)
      integer*4 SLAVE_ID(nslv),JBIN(nslv),IJSTRT(nslv)
      integer binlist(nslv)
      integer ipair(7),jpair(7)
      integer WaitAndLock, Unlock
      parameter (maxpair=700*701/2,maxgran=maxpair/2)   !FRAGILE!
      integer igran2pair(2,maxgran),ipair2gran(maxpair)
      integer*1 i01
      integer*4 i0
      Parameter (ndisk=50,ndisk2=60)
      Parameter (mdisk1=51,mdisk2=52,mdisk3=53,mdisk4=54)    ! UMP2 unit nos.
      Parameter (zero=0.0d0,one=1.0d0)
      Data iprnt/0/
c
c
c -- put down a memory mark
      call matreset
c
      call mmark
      call matmark
c
c -- get the corresponding host name
      call secund(tstart)
c
c -- get initializing data
c -- we first initialize the integral stuff (which reuses SCF call)
c
c -- unpack the geometry,cpu,symmetry data
      call para_get_info
c
c -- get the remaining iniial data
      call para_recv_bcastpack(TxJobInit)
      call para_unpack_int(iforwhat,1)
c
      call para_recv_bcastpack(TxScfInit)
      call para_unpack_int(nfock,1)
      call para_unpack_int(lsemi,1)
      call para_unpack_int(iforwhat,1)
      call para_unpack_int(ncachex,1)
      call para_unpack_int(NAlpha,1)
      call para_unpack_int(NBeta,1)
      call para_unpack_int(rhf,1)
      call para_unpack_real(threshx,1)
      call para_unpack_real(threshy,1)
      call para_unpack_int(iroute,1)
      call para_unpack_int(istab,1)
      call para_unpack_int(idft,1)
      call setival('nmo ',NAlpha)
c
c -- get some necessary basic data from the depository
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ictr',ictr)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('nsym',nsym)
      If(nsym.gt.0) call getival('SymFunPr',ifp)
      call getival('nonredun',lastvirt)
c
c -- get number of dummy atoms (if any)
      call getival('ndum',ndum)
      natom = na-ndum
c
c -- restore coordinates, atomic charges and get symmetry data
      call getmem(3*natom,ixc)                ! nuclear coordinates
      call getmem(natom,iqa)                  ! atomic charges/numbers
c
      call getnucdat(natom,bl(ixc),bl(iqa))
c
      ncf2 = ncf**2
c .........................................................................
c
c  get values necessary for initializing the integrals (twoint95)
c  most of the things were not sent from the master
c  (see comments in ptwoint.f)
c
      call getival('maxprice',maxpricex)
      iroutex=iroute
      call getival('chec',icheckx)
      call getival('iprn',iprintx)
      call getival('istat',istat)
c
c  rest of the things needed for twoint95
c
      call getival('lcore',lcore)
      call getival('nsym',nsym)
c
c -- initialize the local integral system
c2002
      call setival('iroute',iroute)
      call setival('stab',istab)
c2002
      iforwhat = 5
      call twoint95(lcore,nsym,maxpricex,
     *              ncachex,iroutex,
     *              icheckx,iprintx,threshx,istat,iforwhat,
     *              lsemi,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c ..........................................................................
c
c -- allocate memory for MOs and orbital energies
c    (unfortunately need matrix system as this is used throughout subroutines)
      call matdef('cano','q',ncf,ncf)
      call matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      ieorb=mataddr('epsi')
c
      If(.NOT.rhf) Then
        call matdef('canoB','q',ncf,ncf)
        call matdef('epsiB','d',ncf,ncf)
        icanoB=mataddr('canoB')
        ieorbB=mataddr('epsiB')
      EndIf
c
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c
c -- unpack the initial MP2 data
      call para_recv_bcastpack(TxMP2Dat)
      call para_unpack_int(nfirst,1)
      call para_unpack_int(nlast,1)
      call para_unpack_int(nval,1)
      call para_unpack_int(lrec,1)
      If(.NOT.rhf) Then
        call para_unpack_int(nfirstB,1)
        call para_unpack_int(nlastB,1)
        call para_unpack_int(nvalB,1)
        call para_unpack_int(lrec2,1)
        call para_unpack_int(lrec3,1)
        call para_unpack_byte(alphabeta,1)
      EndIf
      call para_unpack_int(smal,1)
      call para_unpack_int(keep,1)
      call para_unpack_int(restrt,1)
      call para_unpack_int(Test,1)
      call para_unpack_int(mp2only,1)
      call para_unpack_int(bl(mpres_in_bg),ncs)
      call para_unpack_real(thresh,1)
      call para_unpack_real(bl(ieorb),ncf)
      call para_unpack_real(bl(icano),ncf*ncf)
      If(.NOT.rhf) Then
        call para_unpack_real(bl(ieorbB),ncf)
        call para_unpack_real(bl(icanoB),ncf*ncf)
      EndIf
c
c -- allocate integer memory for orbital symmetry pairs
c      
      If(nsym.gt.0) Then
        call getint(nsym*nval,iorbprA)
        if(.NOT.rhf) call getint(nsym*nvalB,iorbprB)
      EndIf
c
c -- set up sub-matrix for actual MOs correlated
      call matsub('occu','cano',nfirst,nlast)
      call matsub('virt','cano',NAlpha+1,lastvirt)
      ioccu=mataddr('occu')
      ivrtA=mataddr('virt')
      If(.NOT.rhf) Then
        call matsub('occuB','canoB',nfirstB,nlastB)
        call matsub('virtB','canoB',NBeta+1,lastvirt)
        ioccuB=mataddr('occuB')
        ivrtB=mataddr('virtB')
      EndIf
c
c -- put down a memory marker
      call mmark
c
c -- allocate memory for screening density
      call matdef('dsmx','q',ncs,ncs)
      idics=mataddr('dsmx')
c
c -- get screening density, symmetry info and file basename
c
      call para_recv_bcastpack(TxMP2File)
      call para_unpack_real(bl(idics),ncs*ncs)
      If(nsym.gt.0) Then
        call para_unpack_int(bl(iorbprA),nsym*nval)
        if(.NOT.rhf) call para_unpack_int(bl(iorbprB),nsym*nvalB)
      EndIf
c
c ................................................................
c -- get base file names from master
      call para_unpack_string(filname1,256)
      call para_unpack_int(len1,1)
      call para_unpack_string(filname2,256)
      call para_unpack_int(len2,1)
c ------------------------------------------------------------
c -- file handling for restart
c    (I don't know if this is going to work or not; the idea
c     is NOT to open the file if another process is using it)
c
      INQUIRE(FILE=filname1(1:len1)//'.rst',EXIST=exst)
      If(.NOT.exst) Then
c -- create restart file
        OPEN (UNIT=39,FILE=filname1(1:len1)//'.rst',FORM='FORMATTED',
     $        ACCESS='APPEND',STATUS='UNKNOWN')
        CLOSE (UNIT=39,STATUS='KEEP')
      EndIf
c
      Restart_File = filname1(1:len1)//'.rst'
      LenR = len1+4
      ILock = WaitAndLock(Restart_File,lenR)
c
c -- open restart file
      If(restrt) Then
      OPEN (UNIT=39,FILE=filname1(1:len1)//'.rst',FORM='FORMATTED',
     $      ACCESS='SEQUENTIAL',STATUS='UNKNOWN')
      Else
      OPEN (UNIT=39,FILE=filname1(1:len1)//'.rst',FORM='FORMATTED',
     $      ACCESS='APPEND',STATUS='UNKNOWN')
      EndIf
c ------------------------------------------------------------
c make filenames unique, based on gid
c
      write(ch3,'(".",I2.2)') MY_GID
c
      IF (rhf) THEN
c -- closed shell RHF
         filname1 = filname1(1:len1)//ch3
         len1 = len1+3
         filname2 = filname2(1:len2)//ch3
         len2 = len2+3
         If (mp2only.EQ.-1) Then
            call para_unpack_string(filname3,256)
            call para_unpack_int(len3,1)
            call para_unpack_string(filname4,256)
            call para_unpack_int(len4,1)
            filname3 = filname3(1:len3)//ch3
            filname4 = filname4(1:len4)//ch3
            len3 = len3+3
            len4 = len4+3
         EndIf
      ELSE
c -- open shell UHF
         filname4 = filname2(1:len2)//ch3
         len4 = len2+3
         filname1 = filname1(1:len1)//'.AA'//ch3
         filname2 = filname1(1:len1)//'.AB'//ch3
         filname3 = filname1(1:len1)//'.BB'//ch3
         len1 = len1+6
      ENDIF
c
c ................................................................
c
      If(restrt) Then
c -- get integral file name from <restart> file
         CALL ReStart_MP2(39,filname1,len1)
         GO TO 50
      Else
c -- write integral file name to <restart> file
         WRITE(39,900) filname1,0
 900     Format(A256,I4)
         CLOSE (UNIT=39,STATUS='KEEP')
      EndIf
c
      ILock = UnLock(Restart_File,lenR)
c
c .....................................................................
c -- open direct access file for half-transformed integrals
c
      If (rhf) Then
         OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec)
      Else
         OPEN (UNIT=mdisk1,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec)
         OPEN (UNIT=mdisk2,FILE=filname2(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec2)
         OPEN (UNIT=mdisk3,FILE=filname3(1:len1),FORM='UNFORMATTED',
     $        ACCESS='DIRECT',RECL=lrec3)
      EndIf
c
      irec = 0         ! current record counter
      nrec = 0         ! total number of records written
c .....................................................................
c
      If(restrt) GO TO 50
c
c .....................................................................
c -- allocate remaining memory for MP2
c
      call getint(ncf,mapf2s)             ! mapping array
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
      call getmem(ncf*ncf,ixadr)             ! AO integral matrix
      call matdef('halftra','q',nval,nval)   ! half-transformed exchange operator
      ihalf=mataddr('halftra')
c
c  bl(icol) and bl(jcol) store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
      call getint(ncf,irow)
      call getint(ncf,icol)
      call getint(ncf,irow1)
      call getint(ncf,icol1)
      call getint(ncf,lzero)
c
CSS new integral storage
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
      kcore=lcore-ioffset
      if(.not.smal) then
        lenindj=ncf*ncf*30
        call getint_2(lenindj,indlj)
        call getint(ncf,if2cr)
        call getint(ncf,if2cc)
        call getmem(1,last_entry)
        mp2mem=kcore-last_entry+ioffset-1000000-5*ncf*ncf
        mp2mem=mp2mem/2
        lmp2_size = mp2mem
        call getmem(lmp2_size,lmp2int)
        lenmax = lenindj/6
      endif
c
      If(.NOT.rhf) Then
c
c  reserve space for the various half-transformed integrals
c  ints     holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
        call getint_4(nval**2,intsAA)
        call getint_1(nval**2,int1AA)
        call getint_4(nval*nvalB,intsAB)
        call getint_1(nval*nvalB,int1AB)
        call getint_4(nvalB*nval,intsBA)
        call getint_1(nvalB*nval,int1BA)
        call getint_4(nvalB**2,intsBB)
        call getint_1(nvalB**2,int1BB)
      EndIf
c .....................................................................
c
c -- turn off symmetry for MP2 integrals
      call symmoff
c
c -- initialize
      mulam=0
      mulamd=0
      ENonZero=zero
c
c
c    F I R S T   H A L F - T R A N S F O R M A T I O N
c    -------------------------------------------------
c
 100  CONTINUE
c
c -- request shell pairs we are currently doing from master
c
      call para_send(MY_GID,0,TxDftReq)
      call para_recv_pack(isource,TxMP2Assign)
      call para_unpack_int(ics,1)
      call para_unpack_int(kcs,1)
c
      IF (ics.gt.0) THEN
cc
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf*ncf*ics_size
         call get_shell_size(bl(ictr),kcs,kcs_size)
c
         lmp2_size=lmp2_siz1*kcs_size
c
c check if a size of the lmp2 integral buffer is not too big
c if it is  then split over kcs ( and possibly ics)
c
c
         call check_size1(lmp2_size,ncf,nval,nval,ntimes)
c
         call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                   ntimes,Icspass,Kcspass,Itimes,Ktimes)
c
         do itime=1,itimes
            icf1=icspass(1,itime)
            icf2=icspass(2,itime)
            iatonce=icf2-icf1+1
            do ktime=1,ktimes
               kcf1=kcspass(1,ktime)
               kcf2=kcspass(2,ktime)
               katonce=kcf2-kcf1+1
c
               lmp2_size=iatonce*katonce*ncf*ncf
c
               If (smal) Then
                  call getmem(lmp2_size,lmp2int)
                  call mmark
                  call int_lmp2(bl,bl(ictr),thresh,ics,icf1,icf2,kcs,
     &                          kcf1,kcf2,bl(mapf2s),bl(idics),iprnt,
     &                          bl(lmp2int),nintotal,nrow,ncol,bl(irow),
     &                          bl(icol),bl(lzero))
                  call retmark
               Else
                  ind=0
                  intstore=0
                  call getmem(lmp2_size,lrestor)
                  call mmark
                  call int_lmp2b(bl, bl(ictr), thresh, ics, kcs,
     1                 bl(mapf2s),bl(idics),iprnt,bl(lmp2int),nintotal,
     2                 bl(icol),bl(irow), bl(if2cc),bl(if2cr),bl(indlj),
     &                 bl(lrestor), nrow, ncol, ind,intstore,
     &                 lmp2_size, lenmax)
                  call retmark
                  indmax=max0(indmax,ind)
                  istmax=max0(istmax,intstore)
                  ttint=ttint+nintotal
                  call retmem(1)               ! lrestor
               EndIf
c
c......................................................................
               IF (rhf) Then
                  call TransOneShell(ncf,ncs,nval,ics,kcs,icf1,icf2,
     &                               kcf1,kcf2,ndisk,bl(ictr),
     &                               bl(lmp2int),bl(ioccu),iprnt,thresh,
     &                               bl(ixadr),bl(ihalf),irec,mulam,
     &                               mulamd,nrow,ncol,bl(irow),bl(icol),
     4                               bl(irow1),bl(icol1),ENonZero,smal)
c......................................................................
                  if (smal) call retmem(1)        ! lmp2int
               ELSE
                  call TransOneShlAB(ncf,ncs,nval,nvalB,ics,kcs,icf1,
     &                               icf2,kcf1,kcf2,bl(ictr),
     &                               bl(lmp2int),bl(ioccu),bl(ioccuB),
     &                               alphabeta,iprnt,thresh,bl(ixadr),
     &                               bl(ihalf), bl(ihalf),irec, 
     &                               bl(intsAA),bl(int1AA),bl(intsAB),
     &                               bl(int1AB),bl(intsBA),bl(int1BA),
     &                               bl(intsBB),bl(int1BB),mulam,mulamd,
     &                               nrow,ncol,bl(irow),bl(icol),
     &                               bl(irow1),bl(icol1),ENonZero,smal,
     &                               ihalf, mdisk1, mdisk2, mdisk3)
c......................................................................
                  call retmem(1)          ! lmp2int
c......................................................................
               ENDIF
            enddo    ! over ktime (selected kcf belonging to kcs shell )
         enddo    ! over itime (selected icf belonging to ics shell )
c......................................................................
c
         GO TO 100
cc
      ELSE
c
c -- half-transformation complete
c
         if(.not.smal) call retmem(5)
c
         If (rhf) Then
c -- write end of file record
            i0=0
            i01=0
            nrec = nrec+irec
            write(ndisk,rec=irec+1) -1,-1,(i01,i=1,nval**2),
     $                                (i0,i=1,nval**2)
            CLOSE (UNIT=ndisk,STATUS='KEEP')
         Else
            close (unit=mdisk1,status='keep')
            close (unit=mdisk2,status='keep')
            close (unit=mdisk3,status='keep')
         EndIf
c
      ENDIF
c
 50   CONTINUE
cc      write(6,*) ' SLAVE:',MY_GID,' ** FINISHED FIRST-HALF TRANS **'
cc      call f_lush(6)
c
c -- release memory used in first half-transformation
      call retmark
c
c .....................................................
      If (Test) STOP
c .....................................................
c
c
c    B I N   S O R T  &  S E C O N D   T R A N S F O R M A T I O N
c    -------------------------------------------------------------
c
c
      call flush(6)
      IF (rhf) THEN
c
c    After the half-transformation we have all I,J (MOs) for each w,v (AOs).
c    We need all w,v for each I,J.  Do a standard Yoshimine bin sort to
c    reorder the indices. We are going to read each w,v bin from the half-
c    transformed file and put each I,J into its own bin. When an I,J bin is
c    full it will NOT be written to a file on the original slave, but will
c    be sent to the slave that will ultimately handle that (I,J) pair.
c
c    Send to master amount of available memory on slave.
         call getmem(0,iaddress)
         call retmem(1)
         call getival('lcore',lcore)
         islavemem=lcore-iaddress
         call para_send(islavemem,0,TxMP2S6)

c -- get number of passes and initialisaion data
c    for second half-transformation
c
         call para_recv_bcastpack(TxBlockAssign)
         call para_unpack_int(NPass,1)
         call para_unpack_int(lbin,1)
         call para_unpack_int(igranulesize,1)
         call para_unpack_int(igranules,1)   ! Total no. of granules
         call para_unpack_int(npairs,1)
         call para_unpack_int(igran2pair,2*igranules)
         call para_unpack_int(ipair2gran,npairs)
c
c disable direct routing during binsort (applies only to PVM version)
c (some small MP2 jobs show random errors if direct routing
c  is enabled during the parallel binsort in the PVM version)
c This is now done only for small jobs ( less than 200 basis functions)
c because for large jobs the nodirect route might cause the PVM
c deamon to allocate too much memory and the job will either die or
c start to swap.
c
         if(ncf.lt.200)  call para_nodirectroute
c
         lrec2 =(lbin     +8*lbin   +4  )*igranulesize +8  ! record length in bytes
c
         nvirt = lastvirt-NAlpha      ! number of virtuals
         ncore = NAlpha-nval          ! number of core orbitals
         emp2 = zero
         tnorm = zero
         emps = zero
         empt = zero
         emppar = zero
         empanti = zero
c
c -- memory
c
         call matdef('xmos','q',nvirt,nvirt)   ! exchange matrix
         ixmo=mataddr('xmos')
c
         lastgp=0         ! last granule of previous pass

         ijfirst = 1      ! first ij pair in this current pass
         DO 95 IPass=1,NPass
c
c -- open the direct access file for writing
c
            OPEN (UNIT=ndisk2,FILE=filname2(1:len2),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec2,STATUS='NEW')
cc      write(6,*) ' SLAVE:',MY_GID,' Opened direct access file'
cc      write(6,*) ' Filename is: ',filname2(1:len2),' unit no:',ndisk2
cc      call flush(6)
c
c -- get binlist from Master (for bin write)
c
            call para_recv_bcastpack(TxBinDat)
            call para_unpack_int(binlist,nslv)

c     ngran=binlist(1)-lastgp  ! no. of granules in current pass BUGGY
            if (MY_GID.eq.1) then
               ngran=binlist(1)-lastgp
            else
               ngran=binlist(MY_GID)-binlist(MY_GID-1)
            endif
c
c   set granule number modifier for the current slave
c
            if (MY_GID.eq.1) then
               imod=lastgp
            else
               imod=binlist(MY_GID-1)
            endif          
c
c -- allocate and zero out the bin locator array
c
            call getint(ngran,ibin)
            call izeroit(bl(ibin),ngran)
c
c -- allocate memory for the bin sort
c
            npairs = nval*(nval+1)/2
            call getint_4(npairs*lbin*2,i1)   ! memory for bins
            call getint_1(npairs*lbin,i2)     !  ditto precision overflow
            call getint_4(nval**2,i3)         ! nval*nval integer matrix
            call getint_1(nval**2,i4)         !  ditto precision overflow
            call getint_4(npairs,i5)          ! bin occupancy counter
c
c  the MPI version needs a buffer for the non-blocking send of bins 
c
            msgpass=''
            call getchval('msgpass',msgpass)
            if (msgpass.eq.'MPI')then
               nbbuf=8+8+4*2*lbin*igranulesize+4*igranulesize+
     &               lbin*igranulesize
            else
               nbbuf=0
            endif
            call getbyte(nbbuf,isendb)
c
c -- reopen the half-transformed integral file
c
            OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $            ACCESS='DIRECT',RECL=lrec,STATUS='OLD')
cc      write(6,*) ' SLAVE:',MY_GID,' Opened half-integral file'
cc      write(6,*) ' Filename is: ',filname1(1:len1),' unit no:',ndisk
cc      call flush(6)
c
c -- now sort the bins
c
            ijlast=igran2pair(2,binlist(nslv))! last ij pair in this pass
            keep1 = (IPass.NE.NPass.OR.keep)
c
            call Para_BinSort(
     1           nval,      npairs,    ndisk,       ndisk2,    lbin,
     2           bl(ibin),  ijfirst,   ijlast,      binlist,   bl(i1),
     3           bl(i2),    bl(i3),    bl(i4),      bl(i5),    keep1,
     4           igran2pair,ipair2gran,igranulesize,igranules, ngran,
     5           imod,      bl(isendb))
c
c -- return bin sort memory
c
            call retmem(6)
c
c check if there are bins left to process for our slave
c and eventually process them
c
            call check_bin_left(bl(ibin),ngran,ibinleft)
            if (ibinleft.ne.0)then
c
c -- reopen the bin sort direct access file
c      
               OPEN (UNIT=ndisk2,FILE=filname2(1:len2),RECL=lrec2,
     &               FORM='UNFORMATTED',ACCESS='DIRECT',STATUS='OLD')
c      
               call getint_4(lbin*igranulesize*2,i1) ! memory for integral bin
               call getint_1(lbin*igranulesize,i2)   !  ditto for overflow
               call getmem(ncf*ncf*igranulesize,ixm) ! for bunch of EXCH matrices
c
c -- form contribution to final MP2 energy
c
               icount = 0
               call secund(t1)
               call elapsec(t1e)
c
c  first ij pair on current slave
c
               if (MY_GID.eq.1)then
                  ijfirsts=ijfirst
               else
                  ijfirsts=igran2pair(2,binlist(MY_GID-1)) + 1
               endif
c
c  last ij pair on current slave
c
               ijlasts=igran2pair(2,binlist(MY_GID))

               istart=ipair2gran(ijfirsts) ! first granule on slave
               istop =ipair2gran(ijlasts)  ! last granule on slave
               DO 70 II=istart,istop        ! over granules
                  icount = icount+1
c
                  Call ReadBin(ncf,lbin,ndisk2,icount,ngran,bl(ibin),
     &                         thresh,nsym,nfirst,NAlpha,bl(iorbprA),
     &                         bl(ifp),bl(i1),bl(i2),bl(ixm),igran2pair,
     &                         ipair2gran,II,igranulesize,ncore)
c
c------------------------------------------------------------------------
c -- Additional stuff for MP2 gradients
c
cc      IF(mp2only.LT.0) THEN
c -- Kij in AO basis now in xmat [i.e., in bl(ixm)]
c -- transform to virtual-occupied block and save for gradients
c
cc        call PSaveKIJ(ncf,    nval,   nvirt, bl(ikij),bl(i1),
cc     $               bl(i2), IJ,     thresh, lengvo, ndisk5,
cc     $               ncore)
ccccccc
c NOTE by JB:  need to define matrix 'tvir' (transpose of 'virt')
c              also the matrix 'kijvo' which is the same as bl(ikij)
c              in the argument list. lengvo is the length of the bins
c              and ndisk5 is the unit number of the bins file, which
c              needs to be opened. Seems to me that the last argument
c              ncore can be omitted - if so, change in serial version
c              need to check dimension of bins files bl(i1),bl(i2)
ccccccc
cc      ENDIF
c------------------------------------------------------------------------
c
                  i_ij_start=igran2pair(1,II)  ! Convert to pair number, first pair
                  i_ij_stop =igran2pair(2,II)  ! Convert to pair number, last  pair
                  do i_ij=i_ij_start,i_ij_stop
                     iaddress=ixm+(i_ij-i_ij_start)*ncf*ncf
                     call matconn('xmat','q',ncf,ncf,iaddress)
                     Call matsimtr('xmat','virt','xmos')
                     call matdisc('xmat')
                     call get_ij_half(i_ij,I,J)     ! get separate indices
                     I = I+ncore
                     J = J+ncore
c
                     Call PairEner(ncf,NAlpha,nvirt,I,J,bl(ixmo),
     &                             bl(ieorb),eij,xnorm,eijs,eijt,eijpar,
     &                             eijanti)
c
                     emp2 = emp2 + eij
                     tnorm = tnorm + xnorm
c -- Singlet and triplet components of the correlation energy
                     emps=emps+eijs
                     empt=empt+eijt
                     emppar=emppar+eijpar
                     empanti=empanti+eijanti
                  enddo
c
c------------------------------------------------------------------------
c -- Additional stuff for MP2 gradients
c
cc      IF(mp2only.LT.0) THEN
c -- for MP2-gradients, backtransform to AO basis after
c -- scaling with energy denominators
cc        call PSaveTIJ(ncf,    NAlpha, nval,   nvirt,  I,
cc     $               J,   bl(ieorb), thresh, bl(i1), bl(i2),
cc     $              bl(ixmo),ndisk3, IJ)
cc      ENDIF
c------------------------------------------------------------------------
c
 70            CONTINUE
               call elapsec(t2e)
               call secund(t2)
               t1 = t2-t1
               t1e = t2e-t1e
c
               call retmem(3)          ! release integral bin memory
c
c -- delete the bin-sort direct access file
c
               CLOSE (UNIT=ndisk2,STATUS='DELETE')
            endif   ! ibinleft.ne.0
c
            call retmem(1)     ! release bin locator memory
c
c -- send to master to start next loop
c
            call para_send(MY_GID,0,TxBlockReq)
            ijfirst = ijlast+1
            lastgp = binlist(nslv)
c
 95      CONTINUE
c
c -- do a global sum of the final results
         call para_reduce(emp2,   1,TxMP2S1)
         call para_reduce(tnorm,  1,TxMP2S2)
         call para_reduce(emps,   1,TxMP2S3)
         call para_reduce(empt,   1,TxMP2S4)
         call para_reduce(emppar, 1,TxMP2S5)
         call para_reduce(empanti,1,TxMP2S6)
c
c -- delete the restart file
c     OPEN (UNIT=39,FILE=Restart_File(1:lenR),STATUS='UNKNOWN')
c     CLOSE (UNIT=39,STATUS='DELETE')
c
c -- send timings back to master then move on
         call secund(tend)
         call para_initsend
         call para_pack_real(tend-tstart,1)
         call para_pack_string(MY_HOSTNM,60)
         call para_send_pack(0,TxTime)
c
c restore direct routing (applies only to PVM version)
         if(ncf.lt.200) call para_directroute
c
c
      ELSE        ! -- UMP2 --
c
c    After the half-transformation we have all I,J (MOs) for each w,v (AOs).
c    We need all w,v for each I,J.  We are going to do a modified bin sort.
c    Each slave will sort into IJ bins all of the integrals that it has.
c    The IJ indices will be divided up roughly equally between the slaves
c    and when a block of bins is full, it will be sent to the slave to which
c    it has been assigned which will do the final transformation, compute
c    the pair energies and accumulate them into a partial UMP2 energy.
c    When all integrals have been transformed, the partial energies will be
c    summed up on the master.
c
c -- get message from master to start
c
      call para_recv_bcastpack(TxBinDat)
      call para_unpack_int(IStart,1)
c
c -- send master the size of the integral file (printout)
      nrec = irec      ! number of w,v records on this slave
      call para_send(nrec,0,TxMP2S4)
c
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
c --  lastvirt is less than ncf if redundant basis functions suppressed
      lastvirt=igetival('nonredun')
      nvirtA = lastvirt-NAlpha
      nvirtB = lastvirt-NBeta
c
c-----------------------------------------------------------------------
c  put down a memory marker
      call mmark
c
c-----------------------------------------------------------------------
c  reserve space for the various half-transformed integrals
c  ints     holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
      call getint_4(nval**2,intsAA)
      call getint_1(nval**2,int1AA)
      call getint_4(nval*nvalB,intsAB)
      call getint_1(nval*nvalB,int1AB)
c
c-----------------------------------------------------------------------
c  reserve space for batch of fully transformed integrals plus intermediate
c  storage
c
      nvirt = MAX(nvirtA,nvirtB)
      call getmem(nvirt**2,ixtrn)
      call getmem(ncf*nvirt,iz)
c-----------------------------------------------------------------------
c  determine memory available for the bins. They should be as big
c  as possible but leave space for other allocations
c
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memoccu=lastaddr-ioffset
      memavail=memtot-memoccu
      membin=0.8*memavail                  ! units are doublewords (8-bytes)
c-----------------------------------------------------------------------
c  send the available memory on this slave to the master 
      call para_send(membin,0,TxDftReq)
c
c  get back from the master the slave ID array
      call para_recv_bcastpack(TxMP2S5)
      call para_unpack_int4(SLAVE_ID,nslv)
c-----------------------------------------------------------------------
c
c  initialize
c
      eAA = zero
      eAB = zero
      eBB = zero
      tnorm= zero
c
      tbin = zero
      tsort = zero
      ttran = zero
c
      itA = ieorb+nfirst-1       ! start of alpha orbital energies
      itB = ieorbB+nfirstB-1     ! start of beta  orbital energies
c
c-----------------------------------------------------------------------
c
c  -- we are going to handle each type separately
c
c -- determine the disk storage used by current integral files
      halfspace1 = float(nrec)*float(lrec1)/8.0d0      ! file size in doublewords
      halfspace2 = float(nrec)*float(lrec2)/8.0d0
      halfspace3 = float(nrec)*float(lrec3)/8.0d0
c
      IF(alphabeta) GO TO 999
c
c
c  (1)  ALPHA-ALPHA
c
c----------------------------------------------------------------------
c  reopen the Alpha-Alpha half-transformed integral file
c
      OPEN (UNIT=mdisk1,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec,STATUS='OLD')
c-----------------------------------------------------------------------
c
      call mmark
c
c -- allocate the Wolinski reordering array
      call getint(nval*nval,listAA)
      call getint(nval*nval,listAA_back)
      call makeListAA(nval,bl(listAA),bl(listAA_back))
c
c -- determine the disk storage currently being used
      halfspace = halfspace1 + halfspace2 + halfspace3
c
      npairs = nval**2             !  total number of IJ pairs
c
c  get the bin sort parameters from the master
c
      call para_recv_bcastpack(TxMP2S4)
      call para_unpack_int(nIJ,1)
      call para_unpack_int(nIntPerBin,1)
      call para_unpack_int(nBinPerRec,1)
c
c
c----------------------------------------------------------------------
c  open the Alpha-Alpha bin sort file
c
      lrec4 = nIntPerBin*(5*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c-----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass=0
c
c  return address for multipasses
c
 2100 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1 
      istrtAA = 1 + (lpass-1)*nIJ
      iendAA = MIN(istrtAA+nIJ-1,npairs)
c
c -- determine which slaves will transform which IJ indices
c    and determine numrec & nrem for THIS slave
      call split_IJ(istrtAA, iendAA,  nslv,  MY_GID, SLAVE_ID,
     $           nBinPerRec, JBIN,   IJSTRT, numrec, nrem,
     $              MyID)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
c
      call elapsec(t1)
      Call Para_BinSendAA(nval,   nrec,   mdisk1, nIntPerBin,nBinPerRec,
     $                   numrec,  nrem,   istrtAA,  iendAA,  bl(intsAA),
     $              bl(ints1AA),bl(imu), bl(ilam),bl(intsIJ),bl(int1IJ),
     $               bl(listAA),  mdisk4, SLAVE_ID, JBIN,    IJSTRT,
     $                    MyID,   mrec,   jstrtAA,  jendAA)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c  upon return from <para_BinSendAA> we should have all the integrals with
c  IJ indices that will be transformed on this slave
c    jstrtAA  -  starting IJ index
c    jendAA   -  ending IJ index
c
cc      If(iprnt.gt.0) write(iout,1237) MY_GID,mrec
 1237 format(1X,' SLAVE:',I6,' Number of records on bin sort file: ',I6)
cc      write(6,*) ' SLAVE:  ',MY_GID,' istrtAA is: ',istrtAA
cc      write(6,*) ' SLAVE:  ',MY_GID,' jstrtAA is: ',jstrtAA
cc      write(6,*) ' SLAVE:  ',MY_GID,' jendAA is: ',jendAA
cc      call f_lush(6)
c
c -- allocate memory for sorted integrals
      call retmem(2)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getmem(ncf2*nBinPerRec,ibAA)
c
      DO iloop=1,numrec
c
      jstrt = jstrtAA + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,jendAA)
      call elapsec(t1)
      Call BinSrt1AA(ncf,      nval,     nsym,     mdisk4, nIntPerBin,
     1             nBinPerRec, numrec,   mrec,     iloop, bl(iorbprA),
     2               bl(ifp),  thresh,   jstrt,    jend,  bl(intsIJ),
     3             bl(int1IJ), bl(imu), bl(ilam), bl(ibAA), bl(listAA),
     4                bl(listAA_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtA,   nval,     jstrt,    jend,
     1                bl(ibAA), bl(ivrtA),bl(itA),  bl(iz),  bl(ixtrn),
     2                eAA,      tnorm,  bl(listAA), iprnt,    ibAA,
     3                ivrtA,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(3)
c
c     go back to make next pass
c
      If(iendAA.LT.npairs) GO TO 2100
c
      call retmark
c
c  -- close integral files
      close (unit=mdisk1,status='delete')
      close (unit=mdisk4,status='delete')
C
C  sum up partial alpha-alpha UMP2 energy contribution on the master
      call para_reduce(eAA,   1,TxMP2S1)
c
c
c  (2) BETA-BETA
c
c----------------------------------------------------------------------
c  reopen the Beta-Beta half-transformed integral file
c
      OPEN (UNIT=mdisk3,FILE=filname3(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec3,STATUS='OLD')
c-----------------------------------------------------------------------
c
      call mmark
c
c -- allocate the Wolinski reordering array
      call getint(nvalB*nvalB,listBB)
      call getint(nvalB*nvalB,listBB_back)
      call makeListAA(nvalB,bl(listBB),bl(listBB_back))
c
c -- determine the disk storage currently being used
      halfspace = halfspace2 + halfspace3
c
      npairs = nvalB**2             !  total number of IJ pairs
c
c  get the bin sort parameters from the master
c
      call para_recv_bcastpack(TxMP2S5)
      call para_unpack_int(nIJ,1)
      call para_unpack_int(nIntPerBin,1)
      call para_unpack_int(nBinPerRec,1)
c
c
c----------------------------------------------------------------------
c  open the Beta-Beta bin sort file
c
      lrec4 = nIntPerBin*(5*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c-----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass=0
c
c  return address for multipasses
c
 2200 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1
      istrtBB = 1 + (lpass-1)*nIJ
      iendBB = MIN(istrtBB+nIJ-1,npairs)
c
c -- determine which slaves will transform which IJ indices
c    and determine numrec & nrem for THIS slave
      call split_IJ(istrtBB, iendBB,  nslv,  MY_GID, SLAVE_ID,
     $           nBinPerRec, JBIN,   IJSTRT, numrec, nrem,
     $              MyID)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
c
      call elapsec(t1)
      Call Para_BinSendAA(nvalB,  nrec,   mdisk3, nIntPerBin,nBinPerRec,
     $                   numrec,  nrem,   istrtBB,  iendBB,  bl(intsAA),
     $              bl(ints1AA),bl(imu), bl(ilam),bl(intsIJ),bl(int1IJ),
     $               bl(listBB),  mdisk4, SLAVE_ID, JBIN,    IJSTRT,
     $                    MyID,   mrec,   jstrtBB,  jendBB)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c  upon return from <para_BinSendAA> we should have all the integrals with
c  IJ indices that will be transformed on this slave
c    jstrtBB  -  starting IJ index
c    jendBB   -  ending IJ index
c
cc      If(iprnt.gt.0) write(iout,1237) MY_GID,mrec
cc      write(6,*) ' SLAVE:  ',MY_GID,' istrtBB is: ',istrtBB
cc      write(6,*) ' SLAVE:  ',MY_GID,' jstrtBB is: ',jstrtBB
cc      write(6,*) ' SLAVE:  ',MY_GID,' jendBB is: ',jendBB
cc      call f_lush(6)
c
c -- allocate memory for sorted integrals
      call retmem(2)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getmem(ncf2*nBinPerRec,ibBB)
c
      DO iloop=1,numrec
c
      jstrt = jstrtBB + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,jendBB)
      call elapsec(t1)
      Call BinSrt1AA(ncf,      nvalB,    nsym,     mdisk4, nIntPerBin,
     1             nBinPerRec, numrec,   mrec,     iloop, bl(iorbprB),
     2               bl(ifp),  thresh,   jstrt,    jend,  bl(intsIJ),
     3             bl(int1IJ), bl(imu), bl(ilam), bl(ibBB), bl(listBB),
     4                bl(listBB_back))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAA(ncf,      nvirtB,   nvalB,    jstrt,    jend,
     1                bl(ibBB), bl(ivrtB),bl(itB),  bl(iz),  bl(ixtrn),
     2                eBB,      tnorm,  bl(listBB), iprnt,    ibBB,
     3                ivrtB,    ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(3)
c
c     go back to make next pass
c
      If(iendBB.LT.npairs) GO TO 2200
c
      call retmark
c
c  -- close integral files
      close (unit=mdisk3,status='delete')
      close (unit=mdisk4,status='delete')
C
C  sum up partial beta-beta UMP2 energy contribution on the master
      call para_reduce(eBB,   1,TxMP2S2)
C
  999 CONTINUE
c
c
c  (3) ALPHA-BETA
c
c----------------------------------------------------------------------
c  reopen the Alpha-Beta half-transformed integral file
c
      OPEN (UNIT=mdisk2,FILE=filname2(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec2,STATUS='OLD')
c-----------------------------------------------------------------------
c
c -- determine the disk storage used by current integral files
      halfspace =  halfspace2
c
      npairs = nval*nvalB            !  total number of IJ pairs
c
c  get the bin sort parameters from the master
c
      call para_recv_bcastpack(TxMP2S6)
      call para_unpack_int(nIJ,1)
      call para_unpack_int(nIntPerBin,1)
      call para_unpack_int(nBinPerRec,1)
c
c
c----------------------------------------------------------------------
c  open the Alpha-Beta bin sort file
c
      lrec4 = nIntPerBin*(10*nBinPerRec + 4) + 4      ! record length in bytes
      OPEN (UNIT=mdisk4,FILE=filname4(1:len4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec4,STATUS='NEW')
c----------------------------------------------------------------------
c
c -- allocate the memory for JMu and JLam
      call getmem(nIntPerBin/4+1,imu)          ! integer*2
      call getmem(nIntPerBin/4+1,ilam)         ! integer*2
c
c  Loop over the number of passes
c
      lpass = 0
c
c     return address for multipasses
c
 2300 CONTINUE
c
c  get the range of pairs this pass
c
      lpass = lpass+1
      istrtAB = 1 + (lpass-1)*nIJ
      iendAB = MIN(istrtAB+nIJ-1,npairs)
c
c -- determine which slaves will transform which IJ indices
c    and determine numrec for THIS slave and nrem
      call split_IJ(istrtAB, iendAB,  nslv,  MY_GID, SLAVE_ID,
     $           nBinPerRec, JBIN,   IJSTRT, numrec, nrem,
     $              MyID)
cc      write(6,*) ' SLAVE: ',MY_GID,' Back from <split_IJ>'
cc      write(6,*) ' SLAVE: ',MY_GID,' numrec is: ',numrec
cc      write(6,*) ' SLAVE: ',MY_GID,' nrem is: ',nrem
cc      call f_lush(6)
c
c -- allocate memory for the bins
      call getint_4(nIntPerBin*(nIJ+nrem),intsIJ)
      call getint_1(nIntPerBin*(nIJ+nrem),int1IJ)
      call getint_4(nIntPerBin*(nIJ+nrem),intsJI)
      call getint_1(nIntPerBin*(nIJ+nrem),int1JI)
c
      call elapsec(t1)
      Call Para_BinSendAB(nval,    nvalB,   nrec,   mdisk2, nIntPerBin,
     $                 nBinPerRec, numrec,  nrem,   istrtAB,  iendAB,
     $             bl(intsAA),bl(ints1AA),bl(intsAB),bl(int1AB),bl(imu),
     $             bl(ilam),bl(intsIJ),bl(int1IJ),bl(intsJI),bl(int1JI),
     $                    mdisk4, SLAVE_ID, JBIN,   IJSTRT,   MyID,
     $                    mrec,   jstrtAB,  jendAB)
      call elapsec(t2)
      tbin = tbin + t2-t1
c
c  upon return from <para_BinSendAB> we should have all the integrals with
c  IJ indices that will be transformed on this slave
c    jstrtAB  -  starting IJ index
c    jendAB   -  ending IJ index
c
cc      If(iprnt.gt.0) write(iout,1237) MY_GID,mrec
cc      write(6,*) ' SLAVE:  ',MY_GID,' istrtAB is: ',istrtAB
cc      write(6,*) ' SLAVE:  ',MY_GID,' jstrtAB is: ',jstrtAB
cc      write(6,*) ' SLAVE:  ',MY_GID,' jendAB is: ',jendAB
cc      call f_lush(6)
c
c -- allocate memory for sorted integrals
      call retmem(4)      !  deallocate memory for integral bins
      call getint_4(nIntPerBin*nBinPerRec,intsIJ)
      call getint_1(nIntPerBin*nBinPerRec,int1IJ)
      call getint_4(nIntPerBin*nBinPerRec,intsJI)
      call getint_1(nIntPerBin*nBinPerRec,int1JI)
      call getmem(ncf2*nBinPerRec,ibAB)
c
      Do iloop=1,numrec
c
      jstrt = jstrtAB + (iloop-1)*nBinPerRec
      jend = MIN(jstrt+nBinPerRec-1,jendAB)
c
      call elapsec(t1)
      Call BinSrt1AB(ncf,      nval,     nvalB,    nsym,    mdisk4,
     1          nIntPerBin, nBinPerRec,  numrec,   mrec,    iloop,
     2          bl(iorbprA),bl(iorbprB),bl(ifp), thresh,    jstrt,
     3               jend, bl(intsIJ),bl(int1IJ),bl(intsJI),bl(int1JI),
     4               bl(imu),  bl(ilam), bl(ibAB))
      call elapsec(t2)
      tsort = tsort + t2-t1
c
      Call TrnsVirtAB(ncf,      nvirtA,   nvirtB,   nval,     nvalB,
     1                jstrt,    jend,     bl(ibAB), bl(ivrtA),bl(ivrtB),
     2                bl(itA),  bl(itB),  bl(iz),   bl(ixtrn),eAB,
     3                tnorm,    iprnt,    ibAB,     ivrtA,    ivrtB,
     4                ixtrn,    iz)
      call elapsec(t1)
      ttran = ttran + t1-t2
c
      EndDo
c
      call retmem(5)
c
c     go back to make next pass
c
      If(iendAB.LT.npairs) GO TO 2300
c
c  -- close integral files
      close (unit=mdisk2,status='delete')
      close (unit=mdisk4,status='delete')
C
C  sum up partial alpha-beta UMP2 energy contribution on the master
      call para_reduce(eAB,   1,TxMP2S3)
C
c-----------------------------------------------------------------------
c
c -- that's it!
c
c -- send detailed timings back to the master
      call para_initsend
      call para_pack_real(tbin,1)
      call para_pack_real(tsort,1)
      call para_pack_real(ttran,1)
      call para_send_pack(0,TxMP2S5)
c
c  sum up final wavefunction norm
      call para_reduce(tnorm, 1,TxMP2S4)
c
c -- send timings back to master then move on
      call secund(tend)
      call para_initsend
      call para_pack_real(tend-tstart,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
c
c
      ENDIF
c
c---------------------------------------------------------------------
      call symmon        ! turn symmetry back on for integrals
c---------------------------------------------------------------------
      call retmark
c---------------------------------------------------------------------
c
c -- release allocated memory
      call matremark
      call retmark
c---------------------------------------------------------------------
c
      RETURN
      END
c ============================================================================
c
      SUBROUTINE para_recv_bin(lbin, igranulesize,ngran,imod,
     $                         ndisk,nbin, nrec, islv, wait)
c
c probes for a bin from any slave, if there is any receives
c and writes it to disk. Will go on receiving until all bins received
c      
c lbin:         length of bin
c igranulesize  size of granule
c ngran:        number of granules in this pass
c imod :        granule number modifier for current slave
c ndisk:        binsort file unit no
c nbin:         bin allocation array (which bin on which record)
c nrec:         present record no
c islv:         number of slaves still sorting      
c wait:         if .true. then receive messages until the last (-1)
c      
      use newpara

      INTEGER NBin(ngran)
      INTEGER*4 IBins(2*lbin*igranulesize) ! dynamic allocation 
      integer*4 icount(igranulesize)       ! dynamic allocation
      INTEGER*1 IBin1(lbin*igranulesize)   ! dynamic allocation
      logical isthere,wait,mywait
c
      mywait=wait
      do
c
c                --  probe for a bin
c
         call para_probe(TxBinS1,isthere)
         if (isthere.or.mywait) then
            call para_recv_pack(islave,TxBinS1)
            call para_unpack_int(igran,1)       ! granule number
            call para_unpack_int(isize,1)       ! granule size
            If (igran.ne.-1) Then
c
c                  -- unpack the bin
c
               call para_unpack_int4(icount,isize) ! number of elements in bin
               call para_unpack_int4(ibins,2*lbin*isize)  ! bin itself
               call para_unpack_byte(ibin1,lbin*isize)    ! precision overflow
c          
c                  --  write bin to the direct access file
c
               call para_write_bin(
     $                     ibins,        ibin1,    icount,    lbin,
     $                     igranulesize, nbin,     ngran,     imod,
     $                     igran,        ndisk,    nrec)
            Else           ! one more slave has finished
               islv=islv-1
               if (islv.eq.0) mywait=.false.  ! this pass is over
            Endif
         else             !no bin arrived and not waiting to finish
            return
         endif
      enddo
c
      end
c
c ============================================================================
c
      SUBROUTINE para_write_bin(
     $                ibins,        ibin1,    icount,    lbin,
     $                igranulesize, nbin,     ngran,     imod,
     $                igran,        ndisk,    nrec)
c
c writes a bin out to binsort file, called when
c a bin is received or when the bin is written on the same slave
c      
c ibins:        bin data
c ibin1:        bin overflow data      
c icount:       nummber of elements in bin
c igranulesize: granularity
c lbin:         length of bin
c nbin:         bin allocation array (which bin on which record)
c ngran         no of granules
c imod:         some data I have to make sense of
c igran:        ditto
c ndisk:        binsort file unit no
c nrec:         present record no
c      
      implicit none

      integer lbin,igranulesize,ngran,imod,igran,ndisk,nrec
      INTEGER NBin(ngran)
      INTEGER*4 IBins(2*lbin*igranulesize)
      INTEGER*1 IBin1(lbin*igranulesize)
      integer*4 icount(igranulesize)
c
c --  write bin to the direct access file
c
      nrec = nrec+1
      WRITE(UNIT=ndisk,rec=nrec) NBin(igran-imod),icount,IBin1,IBins
      NBin(igran-imod) = nrec
cc
      end
c
c ============================================================================
c
      subroutine check_bin_left(nbin,ngran,ibinleft)
c
c  check if there are bin to process for the current slave
c
      implicit none
      integer ngran,ibinleft
      integer nbin(ngran)
      integer i
c
      ibinleft=0
      do i=1,ngran
        if(nbin(i).ne.0)ibinleft=1
      enddo
c
      return
      end
c ============================================================================
c
      function generator()
      integer generator
      generator=3
      end
c ============================================================================
c
      SUBROUTINE ReStart_MP2(IUnit,filname1,len1)
      IMPLICIT NONE
C
C  Finds an appropriate half-transformed integral file to open
C  during an MP2 restart
C
C  ARGUMENTS
C
C  IUnit     -  Unit number of restart file
C               Should already be opened prior to calling this routine
C  filname1  -  on exit filename of file to open
C  len1      -  on exit length of <filname1>
C
      Character*256 FName(96),filname1
      INTEGER IUnit,len1
      INTEGER*2 IOpen(96),Int2,NFile,I
C
C  Read the restart file
C
      NFile = 0
 10   CONTINUE
      READ(IUnit,900,END=95) filname1,Int2
c
      NFile = NFile+1
      FName(NFile) = filname1
      IOpen(NFile) = Int2
      GO TO 10
c
 95   CONTINUE
C
C  hit end of file
C  there are NFile integral files on this node
C
cc      write(6,*) ' There are ',nfile,' files on this node'
cc      do i=1,nfile
cc      write(6,*) ' file ',i,' is ',fname(i)(1:40)
cc      enddo
C
C  find a suitable file to open
C  i.e., one that has not already been opened
C
      Do I=1,NFile
      If(IOpen(I).EQ.0) Then
        IOpen(I) = 1
        GO TO 30
      EndIf
      EndDo
C
C  Should not get here
C  No suitable file to open
C
      call nerror(666,'ReStart_MP2',
     $                 'Cannot Find Suitable File to Open',0,0)
c
 30   CONTINUE
      filname1 = FName(I)
      call rmblan(filname1,256,len1)
C
C  found file name and length
C  write filenames and status back to restart file
C
      REWIND IUnit
      Do I=1,NFile
      WRITE(IUnit,900) FName(I),IOpen(I)
      EndDo
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
C
      RETURN
c
 900  Format(A256,I4)
c
      END

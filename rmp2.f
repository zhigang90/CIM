      Subroutine rmp2(nmo)

      use memory

      implicit real*8 (a-h,o-z)
      integer*1 i01
      integer*4 i0
      character*256 jobname,scrfile,filename,filname1,filname2,
     $              filname3,filname4
      logical LastSymPair,emp2only,scs
      dimension xintxx(9)
c
      Data IUnit/1/                     ! unit number for checkpoint I/O
      common /job/jobname,lenJ
c
      parameter(nopt=16)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 chopv(nopt)
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20,cdum*20
      logical nofr,exst,restrt,dualbasis,smal,Test
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension ifilestat(13)
      integer, parameter :: maxgran = 300*301/2
      integer*4 ipair2gran(  maxgran)
      integer*4 igran2pair(2,maxgran)
ckw
      dimension Ipass(2,28), Kpass(2,28) !28 is 28 comp.of cartisian I-function
c-----------------------------------------------------------------------

      data opnames  /'nofr','orbs','thre','core','rest',
     1               'dual','big ','smal','maxd','keep',
     2               'prin','grad','pmij','scs ','type',
     3               'test'/
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c  11=a single real, 12=two reals, 13=three reals, 21=a character option
      data ioptyp / 0, 2,11,11, 0, 0, 0, 0,11, 0, 1, 0, 2,12, 1, 0/
c-----------------------------------------------------------------------
c  Explanation of the options:
c  NOFR = no frozen core correlate all orbitals
c  ORBS = correlate orbitals from i to j
c  THRE = integral thresold, n pH form, default is  9, meaning 1.0d-9
c  CORE = the orbital energy separating the cores from the valence, in
c         atomic units, default is -3 Eh; if set in RHF, that value is
c         used
c  REST = restart, meaning that the transformed integrals are
c         already calculated
c ---------------------------------------------------------------------
c **WARNING**  The MP2 restart (most likely used for parallel jobs) was
c   found to be unreliable, I think by MM. The problem is to do with
c   how the slaves are assigned on the restart - they have to be the
c   SAME as the originals to guarantee a proper restart and it is far
c   from straightforward to ensure this. Consequently the restart has
c   been eliminated from this version.
c ---------------------------------------------------------------------
c  DUAL = dual basis MP2
c  BIG  = use the 'big' memory model. This is the default for calculations
c         with more than 30 atoms. (**WARNING - Turned Off; does not work)
c  SMAL = use the 'small memory model. If BIG and SMAL are both given,
c         then BIG is assumed. Default for calculations with < 30 atoms
c  MAXD = maximum disk space in gigabytes (default 20 GB)
c  KEEP = do not delete the half-transformed integral file and the
c         sort file after the job
c  PRIN = print level, default=2
c  GRAD = calculate the extra quantities needed for an MP2 gradient
c         and write them to disk
c  PMIJ = print par correlation matrices (beware!!) from i to j
c  SCS  = use scaled MP2 (with optional read-in of 2 scale factors)
c  TYPE = dual-basis type (1, 2 or 3)
c  TEST = undocumented feature. Will stop the job after the first
c         half-transformation. Used to test the restart.
c
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
c  basis function symmetry pairs are in inx(ifp)
      iprnt=1
c .............................................................
c -- CHECK! Only do MP2 if current wavefunction is RHF
c --   (at the same time get the lowest eigenvalue of the overlap
c --    and the SCF energy)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,cdum)
      Call rdcntrl(IUnit,5,'$escf',2,idum,erhf,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(wvfnc(1:3).EQ.'RHF' .or. wvfnc(1:3).EQ.'RMP') then
      else
        call nerror(1,'RMP2 Module',
     $  'MP2 Requested But Current Wavefun. is NOT RHF or RMP2 !',0,0)
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
      if(ifound(3).gt.0) then
        thresh=10.0d0**(-ropv(1,3))
      else
c -- check lowest eigenvalue of overlap (from SCF)
        thresh = MIN(1.0d-10,xlow**2)
ccccc  write(6,*)' xlow=',xlow,' intther=',thresh
        If(thresh.LT.1.0d-12) then
          write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final MP2 energy may have greater than mhartree error')
          thresh = 1.0d-12
        EndIf
      end if
      call setrval('thresh',thresh)
c  core is the limit separating core from valence;
c  it is normally set in rhf but it can be set here too
      call tstrval('core',iyes)
      if(iyes.eq.1) then
        core=rgetrval('core')
      else
        core=-3.0d0
      end if
      if(ifound(4).gt.0)  core=ropv(1,4)
c  restart
      restrt=.false.
      if(ifound(5).gt.0) restrt=.true.
c .............................................................................
c -- File handling --
c  filname1=half-transformed integrals; filname2=sorted bins
c  If gradient: filname3=Tij amplitudes; filname4=Kov amplitudes
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.htr'
      len1 = len+4
      filname2=scrfile(1:len)//'.bins'
      len2 = len+5
c .........................................................................
      dualbasis=ifound(6).gt.0
      natoms=igetival('na')
c  *************************************************
c   WARNING!  I don't think this routine properly
c    takes into account dummy atoms   ! JB Nov 07
c  *************************************************
      if(natoms.le.30) then
        smal=.true.
      else
        smal=.false.
      end if
      smal=.true.       ! always use small option
      if(ifound(7).gt.0) smal=.false.
      if(ifound(8).gt.0) smal=.true.
      if(ifound(9).gt.0) then
c  Disk allotment in gigabytes
        xmaxdisk=ropv(1,9)
        xmaxdisk=125 000 000.0d0*xmaxdisk     ! convert to megawords
      else
c  default disk storage is 20 GB
        xmaxdisk=2 500 000 000.0d0
      end if
      if(ifound(11).gt.0) iprnt=iopv(1,11)
c  'PMIJ' print matrix; first number is 1 (AO Xchange) or 2
c   (both AO & MO) second number: pair to be printed
      ipmij=0
      ipairpr=0
      if(ifound(13).gt.0) then
        ipmij=iopv(1,13)
        ipairpr=iopv(2,13)
      end if
c -- SCS scaling  ! SS
      scs=.false.
      iscs=0
      p1=1.2d0
      p2=1.0d0/3.0d0
      If(ifound(14).gt.0) Then
        scs=.true.
        iscs=1
        if(ropv(1,14).gt.0) then
          p11=ropv(1,14)
          p22=ropv(2,14)
          if(abs(p1-p11).gt.1.0d-05) then
            iscs=2
            p1=p11
            p2=p22
            write(iout,1234) p1,p2
 1234       format(' Scaled MP2 with p1: ',F8.5,' p2: ',F8.5)
          endif
        else
          write(iout,*) ' Conventional SCS-MP2'
        endif
      EndIf
c -- dual type (1, 2 or 3)
      itype=1
      If(ifound(15).gt.0) Then
        itype=iopv(1,15)
        if(itype.gt.3) itype=1
      EndIf
c
      Test = ifound(16).gt.0
c-----------------------------------------------------------
      call setival('iscs',iscs)
      call setrval('p1',p1)
      call setrval('p2',p2)
c-----------------------------------------------------------
c -- determine logical variable emp2only
c    normally FORCE is the next line for an MP2 gradient and
c    we are going to check for this automatically here
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
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
c-----------------------------------------------------------
  45  format(72('-'))
  46  format(/72('='))
      write(iout,46)
      write(iout,*) '                           The RMP2 Module  '
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
c------------------------------------------------------------
      if(dualbasis) then
c        the current "extended" basis set was already read-in
c        after SCF. Get the smaller bais set info (one used for SCF):
c
         call get_small_basis(bl,ncs,ncf,bl(ictr), ncs_sm,ncf_sm)
c
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
      else
         ncs_sm=ncs
         ncf_sm=ncf
         call set_mpres(bl(mpres_in_bg),ncs)
      endif
c-----------------------------------------------------------
      write(iout,*) ' MP2 integral thresh    = ',thresh
      write(iout,*) ' core orbitals with eps < ', core
      write(iout,*) ' max disk storage (MW)  = ',xmaxdisk*1.d-6
      write(iout,*) ' '
      call f_lush(iout)
c-----------------------------------------------------------
c  put down a memory marker
      call matreset
c
      call mmark
      call matmark
c
c-----------------------------------------------------------
c zero out irrelevant dft stuff
      nfock=1
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
c-------------------------------------------------------
c allocate memory for an array mapping contr.func. to contr.shells :
c
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c
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
c get eigenvectors & eigenvalues
c
c the Fock operator is defined in a small basis set
c but its matrix representation goes over the big one
c It will be diagonalized WITHOUT allowing occupied
c and virtuals orbitals to mix.
c
        call get_mix_eigen2(bl,nblocks,bl(ictr),labels,thresh,
     *                     ncf,ncf_sm,ncs,ncs_sm,nmo,itype)
c
c  integral threshold can be changed on return if "big"
c  basis set overlap shows linear dependency
c
      else
c
c -- set up the MOs filename
      call tstchval('mosfname',iyes)
      if(iyes.eq.1) then
        call getchval('mosfname',filename)
      else
        filename = jobname(1:lenJ)
      endif
      call rmblan2(filename,256,lenM)
      if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
      lenM=lenM+4
c
c -- read in the MOs
      itype = 1
      CALL ReadMOS(ncf,bl(icano),bl(iepsi),.True.,lenM,
     $             filename(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(3,'RMP2 Module',
     $   'MOS File Does Not Exist',0,0)
      endif
c-------------------------------------------------------
c initialize the two-el. integral program
c
      thint=thresh
      iforwhat=5
      call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             scftype,xintxx, nblocks,maxbuffer,maxlabels)
c
c-------------------------------------------------------
c
c  determine the orbitals to be correlated
c
      if(ifound(2).gt.0) then
        nfirst=iopv(1,2)
        nlast=iopv(2,2)
        nval=nlast-nfirst+1
        ncore=nmo-nval
      else
        if(nofr) then
          nfirst=1
          nlast=nmo
          nval=nmo
          ncore=0
        else
          ncore=ncores(ncf,bl(iepsi),core)
          nfirst=ncore+1
          nlast=nmo
          nval=nmo-ncore
        endif
      endif
c---------------------------------------------------------------------
c check the maximum values of BOTH
c the occupied and virtual eigenvectors
c
      iocco=mataddr('cano')
      iepsi=mataddr('epsi')
      call check_coef(ncs,ncf,ncore,nmo,nval,bl(iocco) )
c---------------------------------------------------------------------
c Construct "screening density" matrix to be used to estimate the
c importance of integrals. It can be based either on canonical or
c localized orbilatls. The latter has some advantage : integral's
c prescreening based on localizability arguments is supposed to be
c more efficient.
c
c For ordinary one-basis set MP2 the "screening density" is made
c out of canonical or localized vectors taken from the disk as
c "evec_rhf" or "loca_rhf".
c
c For dual-basis set the "screening density is formed using current
c MOS content. In the input deck, the GUESS=READ statement MUST precede
c the MP2 or LMP2 keyword in order to ensure the "corresponding orbital"
c transformation (projection from small to big basis set) has been
c preformed.
c
c---------------------------------------------------------------------
c For dual basis set make screening density using CANONICAL orbitals
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
c      if(.not.dualbasis) then
c        check if Localization was performed
c        np4=4
c        call sudat(np4,'loca_rhf',ni)
c
c        if(ni.gt.0) then
c          loca_done=1
c     write(iout,*)' Localized orbitals are used for integral screening'
c           if(.not.restrt)
c    1      call DmxMakeL(dualbasis, nval, nmo, ncf, ncs,
c    2                    np4,.false.,inx, iprnt, bl(idics))
c        else
c           loca_done=0
c           write(iout,*)
c    1 ' Canonical orbitals are used for integral screeing; results',
c    2 ' may be inaccurate for large molecules. Please use the LOCA',
c    3 ' option in SCF (say LOCA=PIPEK) to improve accuracy and',
c    4 ' performance for large systems'
c           call nerror(4,'RMP2 Module',
c    1  'Please use localized orbitals in SCF for screening',0,0)
c           if(.not.restrt)
c    1    call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
c    2                  np4, inx, iprnt)
c        endif
c
c     else      !  dual-basis
c      write(iout,*)' Occupied orbitals are used for integral screening'
c       call DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
c    *                np4, inx, iprnt)
c     endif
c---------------------------------------------------------------------
cccc    call matprint('dsmx',6)
c---------------------------------------------------------------------
c
c  establish orbital symmetry characters
      if(nsym.gt.0) then
        call getint(nval*nsym,iorbpair)
        call getmem(ncf,itmp)
        call OrbSymPair(nsym, ncf ,   bl(icano), bl(ifp) ,nfirst,
     1                   nlast,iprnt, bl(itmp),  bl(iorbpair),iout)
        call retmem(1)
      endif
c
      lastvirt=igetival('nonredun')
c  print the data of the calculation early
      write(iout,'("  Number of contracted functions =",i8,/,
     1             "  Number of correlated orbitals  =",i8,/,
     2             "  Number of virtual orbitals     =",i8)')
     3                ncf,nval,lastvirt-nmo
      write(iout,*) ' '
      call f_lush(iout)
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
      if(ncore.gt.0) call matsub('occa','cano',1,nmo)
      call  matsub('occu','cano',nfirst,nlast)
      call  matsub('virt','cano',nmo+1,lastvirt)
      ioccu=mataddr('occu')
c.................................................
      if(iprnt.gt.3) then
        call matprint ('occu',6)
        call matprint ('virt',6)
      end if
c .................................................
c  open direct access scratch file to store half-transformed integrals
c
c  the length of a block is given in bytes.
      ints = 4             ! size of an integer in bytes
c     lrec = 4*(nval**2+2) + nval**2
      lrec = 4*(nval**2+4) + nval**2
      nrec1 = ncf*(ncf+1)/2
      write(iout,1235) lrec
 1235 format(' Record length for integral file (bytes): ',I10)
c
      ndisk =41        ! unit number for half-transformed integral files
      ndisk2=42        ! unit number of binsort file
      ndisk3=43        ! unit number for TIJ file (needed for gradient)
      ndisk5=44        ! unit number for KOV file (ditto)
c .................................................
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec)
c .................................................
c
c   reserve space for one AO integral matrix
      call matdef('xmat','q',ncf,ncf)
      ixadr=mataddr('xmat')
c space for one half-transformed  exchange operator
      call matdef('halftra','q',nval,nval)
      ihalftra=mataddr('halftra')
c  put down a memory marker
      call mmark
c
      thint=thresh
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c  nrec is the total number of records written on all files
c  irec is the current counter in the file
      nrec=0
      irec=0
      if(restrt) then
        mulam=ncf*(ncf-1)/2
        mulamd=ncf
        call symmoff
        go to 50
      end if
c---------------------------------------------------
c     if(dualbasis.and.itype.eq.1) then
c        update content of "basis sets" according to current needs
c        i.e. put small basis set into "basis 2" & "basis 4"
c
c        call update_basis(bl,ncs_sm)
c     endif
c---------------------------------------------------
c
      call secund(tt)
      if(iprnt.gt.3) write(iout,*) 'Startup time=',tt-tt0
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c  iktot is the number of ics,kcs pairs really calculated & transformed
c   if an (ics,kcs) pair is skipped because there are no integrals in it
c   then it is NOT incremented
      iktot=0
c  Count the actual number of (mu, lam) pairs, both the diagonal
c   and the non-diagonal ones
      mulam=0
      mulamd=0
      ENonZero=zero
ckw..........
      call secund(txxx1)
      call init_info('fock')
ckw..........
c
c turn off symmetry for mp2 integrals :
      call symmoff
c
c  bl(icol) and bl(jcol) store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
      call getint(ncf,irow)
      call getint(ncf,icol)
      call getint(ncf,irow1)
      call getint(ncf,icol1)
      call getint(ncf,lzero)
C
C     the following for 'BIG' option, SS
C
      if(.not.smal) then
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
      kcore=lcore-ioffset
        call getint(ncf,if2cr)
        call getint(ncf,if2cc)
        call getmem(1,last_entry)
        mem=kcore-last_entry+ioffset-6*ncf*ncf-2000000
        memx=mem/9
        mem=memx*4
        lenindj=memx*4
        call getmem(memx,indlj)
        write(iout,*) ' BIG option used for integrals'
        write(iout,*) ' memory assigned to the job:     ',lcore-ioffset
        write(iout,*) ' memory available for integrals: ',mem, '*2'
        lmp2_size = mem
        call getmem(lmp2_size,lmp2int)
        write(iout,*) ' lmp2int start:',lmp2int-ioffset,mem,' long'
        write(iout,*) ' end integral storage',lmp2int-ioffset+2*mem
        indmax=0
        istmax=0
        lenmax=lenindj/6
      endif
c--------------------------------------------------------------------
      ncf2=ncf*ncf
c--------------------------------------------------------------------
      nskipped=0
      do ics=ncs,1,-1
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf*ncf*ics_size
         do kcs=ics,1,-1
            if(LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c  of each (ics,kcs) pair, calculate only the last one
c  (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
c
            call get_shell_size(bl(ictr),kcs,kcs_size)
c
            IF (smal) THEN
               lmp2_size=lmp2_siz1*kcs_size
c
c check if size of the mp2 integral buffer is not too big
c if it is then split over kcs (and possibly ics)
c
               call check_size1(lmp2_size,ncf,nval,nval,ntimes)
c
               call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                         ntimes,Ipass,Kpass,Itimes,Ktimes)
c
               do itime=1,itimes
                  icf1=ipass(1,itime)
                  icf2=ipass(2,itime)
                  iatonce=icf2-icf1+1
                  do ktime=1,ktimes
c
                     call secund(tt2)
                     call elapsec(telap2)
c
                     kcf1=kpass(1,ktime)
                     kcf2=kpass(2,ktime)
                     katonce=kcf2-kcf1+1
c
                     lmp2_size=iatonce*katonce*ncf2
c
                     call getmem(lmp2_size,lmp2int)
                     call mmark
                     call int_lmp2(bl,bl(ictr),thresh,ics,icf1,icf2,kcs,
     1                             kcf1,kcf2,bl(mapf2s),bl(idics),iprnt,
     2                             bl(lmp2int),nintotal,nrow,ncol,
     3                             bl(irow),bl(icol),bl(lzero))
                     call retmark
                     if (nintotal.eq.0) then
                        nskipped=nskipped+1
                        call retmem(1)
                        cycle
                     else
                        iktot=iktot+1
                     endif
c 
                     call secund(tt3)
                     call elapsec(telap3)
                     tinteg=tinteg+tt3-tt2
                     elapint=elapint+telap3-telap2
c......................................................................
                     call TransOneShell(ncf,ncs,nval,ics,kcs,icf1,icf2,
     &                                  kcf1,kcf2,ndisk,bl(ictr),
     &                                  bl(lmp2int),bl(ioccu),iprnt,
     &                                  thresh,bl(ixadr),bl(ihalftra), 
     &                                  irec,mulam,mulamd,nrow,ncol,
     &                                  bl(irow), bl(icol),bl(irow1),
     &                                  bl(icol1),ENonZero,smal)
c
                     call secund(tt4)
                     call elapsec(telap4)
                     ttrans=ttrans+tt4-tt3
                     elaptrans=elaptrans+telap4-telap3
c......................................................................
                     call retmem(1)          ! lmp2int
c......................................................................
                  enddo    ! over ktime (selected kcf belonging to kcs shell )
               enddo    ! over itime (selected icf belonging to ics shell )
c......................................................................
            ELSE
c
               call secund(tt2)
               call elapsec(telap2)
c
               ind=0
               intstore=0
               call getmem(lmp2_size,lrestor)
               call mmark
               call int_lmp2b(bl,bl(ictr),thresh,ics,kcs,
     1              bl(mapf2s),bl(idics),iprnt,  bl(lmp2int),nintotal,
     2              bl(icol), bl(irow),  bl(if2cc),bl(if2cr),bl(indlj),
     3              bl(lrestor),  nrow,    ncol,     ind,     intstore,
     4              lmp2_size, lenmax)
               call retmark
               ttint=ttint+nintotal
               call retmem(1)          ! lrestor
               if (nintotal.eq.0) then
                  nskipped=nskipped+1
                  cycle
               else
                  iktot=iktot+1
               end if
c
               call secund(tt3)
               call elapsec(telap3)
               tinteg=tinteg+tt3-tt2
               elapint=elapint+telap3-telap2
               call getinfs(bl(ictr),ics,kcs,icf1,icf2,kcf1,kcf2)
c......................................................................
               call TransOneShell(ncf,    ncs,    nval,   ics,    kcs,
     *                     icf1,   icf2,   kcf1,   kcf2,   ndisk,
     1                   bl(ictr),bl(lmp2int),bl(ioccu),iprnt,
     2                   thresh,bl(ixadr),bl(ihalftra), irec, mulam,
     3                   mulamd,   nrow,     ncol,  bl(irow), bl(icol),
     4                   bl(irow1),bl(icol1),ENonZero,smal)
c
               call secund(tt4)
               call elapsec(telap4)
               ttrans=ttrans+tt4-tt3
               elaptrans=elaptrans+telap4-telap3
            ENDIF
         end do     !  over kcs shell
      end do     !  over ics shell
c-------------------------------------------------------------------
      if(nskipped.gt.0) then
        write(iout,*) nskipped,
     1   ' pairs were skipped because no integrals were calculated'
      end if
      call retmem(5)
c......................................................................
c timing here is irrelevant :
      call secund(txxx2)
      txxx0=txxx1
      call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
c
      if(iprnt.gt.2) then
         if(nsym.gt.0) then
            write(iout,*) 'Number of shell pairs omitted by symmetry=',
     1       inegl, ' retained=',iret,'transformed',iktot
         endif
      endif
c......................................................................
      nrec=nrec+irec
      if(iprnt.ge.1) then
         write(iout,'(a,i10)')' Total number of records written=',nrec
      endif
c -- write end-of-file record
      i0=0
      i01=0
      write(ndisk,rec=irec+1) -1,-1,(i01,i=1,nval**2),
     $                                       (i0,i=1,nval**2)
      CLOSE (UNIT=ndisk,STATUS='KEEP')
c
      if(iprnt.ge.1) then
c        Calculate statistics
         PerCent=100.0d0*ENonZero/(dble(mulam+mulamd)*ncf**2)
         write(iout,40) PerCent
      endif
c
  40  format(
     *' Percentage of the AO matrix used in the transformation=',f7.3/)
c......................................................................
c  release memory
  50  continue
      call retmark
      call matrem('halftra')
      call matrem('xmat')
c......................................................................
      write(iout,61)
   61 format('  CPU & Elapsed timings in the MP2 module ')
      write(iout,*) '  '
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format('  Integrals for MP2  =',f8.2,' and ',f8.2,' minutes')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format('  First half transf. =',f8.2,' and ',f8.2,' minutes')
      call f_lush(iout)
ckw---------------------------------------------------------
      if(.not.restrt .and. iprnt.ge.2) then
         call getrval('timeblks',timeblks)
         write(91,123) timeblks/sixty
  123    format('CPU time for new blocking =',f8.2)
      endif
ckw---------------------------------------------------------
      call secund(tt1)
      call elapsec(elaps1)
c
c .....................................................
      If(Test) Then
        write(iout,65)
   65   format('  == Calculation Halted at User Request ==')
        STOP
       EndIf
c .....................................................
c
c..................................................
c  reopen the first half-transformed integral file for the next step
c
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec)
c
c -- determine the file size
      halfspace = float(nrec)*float(lrec)/8.0d0
c..................................................
      TBinSort=zero
c
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
c -- lastvirt is less than ncf if redundant basis functions suppressed
      lastvirt=igetival('nonredun')
      nvirt = lastvirt-nmo
c  determine memory available for the bins. They should be as big
c  as possible but leave space for other allocations
c  virtual memory size
c  available memory
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memvirt=memtot
      memoccu=lastaddr-ioffset
      memavail=memvirt-memoccu
c  memory for the bins is the minimum of the available virtual and real memory
c  size, minus 1 MW for miscellaneous - the space needed for matrix multiplies etc.
      memtrans=3*nval*(nval+1)/2+2*nvirt**2+2*ncf*nvirt
      mem=memavail-memtrans-1 000 000  
      npairs=nval*(nval+1)/2
c  lbin is in units of 9-bytes (2 indices           - integer*2
c                               compressed integral - integer*4
c                               precision overflow  - integer*1)
      lbin=8*mem/(9*npairs)
      if(lbin.lt.100) then
        call nerror(5,'RMP2 Module',
     1  'memory available for bins leads to bin size < 100, mem,lbin',
     2  mem,lbin)
      end if
c  lbin is the length of a bin in 9-byte words
      if(lbin.gt.ncf*(ncf+1)) then
        lbin=ncf*(ncf+1)
      end if
CTJ granules of bins will be written
C*******************************************************************************
      do ! over lbin
        igranulesize=mem/(lbin*9/8+1+ncf*ncf)
        if (igranulesize.gt.npairs) igranulesize=npairs
cc        write(6,'(A,I5)') "Max granule size: ",igranulesize
        do
          lrec =(lbin     +8*lbin       )*igranulesize ! record length in bytes
c Limit record size to 8 MB.
          if (lrec.gt.8 388 608) then
            if (igranulesize.eq.1) exit
            igranulesize=igranulesize-1
          else
            exit
          endif
          igranules=npairs/igranulesize
          if (mod(npairs,igranulesize).ne.0) igranules=igranules+1
        enddo
        call optimal_isize(igranulesize,npairs)
        igranules=npairs/igranulesize
        if (mod(npairs,igranulesize).ne.0) igranules=igranules+1
c
        if ((igranulesize*igranules*lbin*9/8+1).lt.mem) exit
        lbin=lbin-1
cc        write(6,*) "Lbin reduced by one"
        if(lbin.lt.50) then
          call nerror(5,'RMP2 Module',
     1    'memory available for bins leads to bin size < 100, mem,lbin',
     2    mem,lbin)
        end if
      enddo ! lbin
c
      if (igranules.gt.maxgran.or.npairs.gt.maxgran) 
     *   call nerror(6,'RMP2 module',
     $   'Too many granules or pairs for array allocator size',
     $    igranules,npairs)
      call f_lush(iout)
c build maps: granule -> ipairstart,ipairstop
c             ipair   -> granule
      istart=1
      do i=1,igranules
        istop=istart+igranulesize-1
        if (istop.gt.npairs) istop=npairs
        if (istart.gt.npairs) call nerror(7,'RMP2 module',
     $   'Error in granule maps determination ',istart,npairs)
        igran2pair(1,i)=istart
        igran2pair(2,i)=istop
        do j=istart,istop
          ipair2gran(j)=i
        enddo
        istart=istop+1
      enddo
c
      if(iprnt.gt.4) then
        do i=1,igranules
        print *,'Granule: ',i,'   start: ',igran2pair(1,i),'   stop: ',
     *          igran2pair(2,i)
        enddo
        call f_lush(iout)
      endif
C*******************************************************************************
c
c  estimate the total number of bins (must exceed the true value)
c  For each pair, ncf*(ncf-1) non-diagonal & ncf diagonal elements
c  are stored. However, because mu and lambda go over shells,
c  there may be a few more elements
c      NBinsPerPair=((ncf*(ncf+1)+2*ncf-1)/lbin+1)
c  Each non-diagonal (mu, lam) pair is used in each pair twice:
c   as (mu, lam) and as (lam, mu). Diagonal index pairs are used
c   only once
       NBinsPerPair=(2*mulam+mulamd+1)/lbin+1
c  In the symmetrical case, the bins hold ALL (i>=j) pairs
c   but only the symmetry-unique (mu, lambda) indices (by shell,
c   not by basis function, i.e. symmetry-redundant (mu, lambda)
c   pairs may be there if they belong to the same shell
c  The full matrix is restored when the bins are read.
c  Unfortunately, symmetry is not being utilized in the second half
c   transformation.
      NBins=NBinsPerPair*npairs
c  reserve memory for a counter showing which bin belongs to which pair
      call getint(NBins,ibinctr)
c  reserve memory for nvalxnval compressed half-transformed integrals
c     call getmem(nval**2/2+1,i1pair)
      call getint_4(nval**2,i1pair)
c     call getmem(nval**2/8+1,int1)
      call getint_1(nval**2,int1)
c  reserve memory for counter showing occupancy of 1 bin
c     call getmem(npairs/2+1,icounter)
      call getint(npairs,icounter)
c
c
      LDiskPerGran=NBinsPerPair*lbin*igranulesize
ctj
      lrec =(lbin     +8*lbin       )*igranulesize
      open(unit=ndisk2,file=filname2,access='direct',
     1     recl=lrec)
c
c  get a value for nrec if it does not exist
      if(nrec.eq.0) nrec=nrec1
c  define external exchange matrix
      call matdef('xmos','q',nvirt,nvirt)
      ixmos=mataddr('xmos')
c  determine the number of batches for the bin sort (to avoid too much
c  disk space)
c  nbatch is the number of pair batches AFTER THE FIRST (ODD) one
c  halfspace is the space taken up by 1/2 transformed integrals (W)
c
      if(iprnt.gt.1) then
       write(iout,'(a,f15.0)')' Maximum Disk Space requested (bytes):',
     $                  8*xmaxdisk
       write(iout,'(a,f15.0)')' Disk Space used for integrals (bytes):',
     $                   8*halfspace
      endif
c
c  maximum disk available for bin sort
      xmaxd=xmaxdisk-halfspace
c     write(iout,*) 'maxd before limiting',xmaxd
c maxd is limited to 230 MWords (230 000 000 words = 1.8 MB)
c     if(xmaxd.gt.2.3d8) then
c       maxd=230 000 000
c     else
c       maxd=xmaxd
c     end if
c     write(iout,*) 'maxd after limiting',maxd
      if(xmaxd.lt.10 000 000) then
        maxdisk=xmaxdisk
        call nerror(8,'RMP2 Module',
     1 'Less than 10 MW disk storage remains after half transformation',
     2  int(xmaxdisk),int(xmaxd))
      end if
      xldisk=dble(lbin)*dble(nbins)
      if(xldisk.gt.xmaxd) then
        nbatches=xldisk/xmaxd
        NGransInBatch=xmaxd/LDiskPerGran
        NeedDisk=NGransinBatch*LDiskPerGran*8
      else
        NGransInBatch=igranules
        nbatches=0
        NeedDisk=nbins*lbin*9
      end if
c
c  to save space in the sort+ virtual transformation section,
c  only part of the half-transformed integrals are sorted at a time
c  the range of pairs is determined by the available disk space
c  MaxDisk. Subtract the space needed for the half-transformed integrals
c  (or zero for a restart) to give maxd, the disk area (in MegaWords)
c  available for the sorted integrals. Ldisk is the total requirement
c  for sorted integrals (in MW). Recall that the sorted integrals are
c  stored with indices (mu,lambda) in integer (*scale factor) form,
c  and thus each integral takes 4 bytes. The default scale factor is
c  1.0d-9 . Dividing the maximum available disk space into the total
c  required gives the number of batches. NBatch is actually 1 fewer
c  than the number of batches, as there is always at least 1 batch.
c
      if(iprnt.gt.1) then
         write(iout,'(a,i6)') ' Batches used in Yoshimine sort=',
     $                        nbatches+1
         write(iout,*) 'NGransInBatch,LdiskPerGran,NeedDisk',
     1     NGransInBatch,LdiskPerGran,NeedDisk
      endif
C..................................................
      DiskNeed=dble(NeedDisk)/1048576.0d0
      NOdd=igranules-nbatches*NGransInBatch
C..................................................
CSS
      IF(.not.emp2only) THEN
        write(iout,*)' *** Quantitites Stored for MP2 gradient ***'
        nvmo=nval
        if(ncore.gt.0) nvmo=nmo
        llrec=nvirt*nvirt
        krec=5*llrec        ! record length in bytes
        open(unit=ndisk3,file=filname3(1:len3),form='unformatted',
     1       access='direct',recl=krec)
        npars=nval*(nval+1)/2
cc
        If(iprnt.gt.1) then
        Write(iout,*)' record length:',llrec,' integers or ',
     $                 krec,' bytes'
        write(iout,*)' number of pairs:',npars
        endif
c
c -- also open files for virtual occupied blocks of kij
        lxrec=nvmo*nvirt*2
        kxrec=5*lxrec       ! record length in bytes
        open(unit=ndisk5,file=filname4(1:len4),form='unformatted',
     1       access='direct',recl=kxrec)
        npars=nval*(nval+1)/2
cc
        If(iprnt.gt.1) then
        Write(iout,*)' record length:',lxrec,' integers or ',
     $                 kxrec,' bytes'
        write(iout,*)' number of pairs:',npars
        endif
cc
        call matdef('tvir','r',nvirt,ncf)
        call matpose2('virt','tvir','n')
        call f_lush(6)
      ENDIF     !   IF(.not.emp2only) THEN
CSS
c
      if(iprnt.gt.1) then
         write(iout,'("Yoshimine BinSort ",i9," bins",i7," long ",
     1    "Disk space=",f8.2," MBytes")') nbins,lbin,DiskNeed
      endif
C..................................................
      tmp2e=zero
      tmbtr=zero
      ttkij=zero
c  zero out correlation energy, wavefunction norm
      emp2=zero
      emps=zero
      empt=zero
      emppar=zero
      empanti=zero
      tnorm=one

c  write header if orbital energies are to be printed
      if(iprnt.gt.2) then
        write(iout,*) ' '
        write(iout,*) 'Orbitals Pair No. occuPair energy  Pair norm'
      end if
      call secund(tt2)
      call elapsec(elaps2)
      nxrec=0
      mxrec=0
c  do first the odd block
      igrfirst=1
      igrlast=NOdd
      call SortTransf(
     1   ncf,         nval,       npairs,     ndisk,      ndisk2,
     2   nrec,        lbin,       nbins,      igrfirst,   igrlast,
     3   iout,       emp2only,    'virt',     bl(i1pair), bl(int1),
     4   bl(icounter),bl(ibinctr),lrec,       nfirst,     nmo,
     5   nvirt,       thresh,     nsym,       bl(iorbpair),bl(ifp),
     6                iprnt,      ipmij,      ipairpr,    bl(ixmos),
     7   bl(iepsi),   emp2,       tnorm,      emps,       empt,
     8   empanti,     emppar,     TBinSort,   ndisk3,     krec,
     9   tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     $   ttkij,       tmp2e,      ncore,      igran2pair,ipair2gran,
     1   igranulesize)
      igrfirst=NOdd+1
      do ibatches=1,nbatches
c  the first (odd) batch is not included in nbatches
      igrlast=igrfirst-1+NGransInBatch
      call SortTransf(
     1   ncf,         nval,       npairs,     ndisk,      ndisk2,
     2   nrec,        lbin,       nbins,      igrfirst,    igrlast,
     3   iout,       emp2only,    'virt',     bl(i1pair), bl(int1),
     4   bl(icounter),bl(ibinctr),lrec,       nfirst,     nmo,
     5   nvirt,       thresh,     nsym,       bl(iorbpair),bl(ifp),
     6                iprnt,      ipmij,      ipairpr,    bl(ixmos),
     7   bl(iepsi),   emp2,       tnorm,      emps,       empt,
     8   empanti,     emppar,     TBinSort,   ndisk3,     krec,
     9   tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     $   ttkij,       tmp2e,      ncore,      igran2pair,ipair2gran,
     1   igranulesize)
      igrfirst=igrlast+1
      end do
c---------------------------------------------------------------------
      IF(.not.emp2only) THEN
        call matrem('tvir')
        if(iprnt.gt.1) then
          write(iout,1010) tmbtr/sixty,ttkij/sixty,tmp2e/sixty
 1010     format(' Time for saving <Tij> in MO basis:',f8.2,' min',/,
     $           ' Time for saving <Kov> block:      ',f8.2,' min',/,
     $           ' Time for MP2 energy calculation:  ',f8.2,' min')
          write(iout,*) ' <Tij>: ',nxrec,' records written on unit ',
     1                   ndisk3,' record length ',5*llrec
          write(iout,*) ' <Kov>: ',mxrec,' records written on unit ',
     1                   ndisk5,' record length ',5*lxrec
        endif
c -- close all files
        close(unit=ndisk3,status='keep')
        close(unit=ndisk5,status='keep')
      ENDIF
c---------------------------------------------------------------------
c  if ifound(10) is greater than 0 then the KEEP option is on
c ** KEEP OPTION CURRENTLY DISABLED - JB **
      if(ifound(10).le.0) then
        inquire(ndisk,exist=exst)
        if(exst) close(ndisk,status='keep')
        close(unit=ndisk2,status='delete')
      end if
c---------------------------------------------------------------------
c print out timings :
        call secund(tt3)
        call elapsec(elaps3)
c........................................................
      write(91,62)
      write(91,64) tinteg/sixty, elapint/sixty
c
      write(iout,704) TBinSort/sixty,nbins
  704 format('  Binary sorting     =             ',f8.2,' minutes'/
     *       '  (number of bins    =',i8,')' )
      write(iout,705) (tt3-tt2)/sixty,(elaps3-elaps2)/sixty
  705 format('  sort & virt.trans. =',f8.2,' and ',f8.2,' minutes')
      write(iout,*) ' '
c Print out the final results :
c
      tnorm=sqrt(tnorm)
c
      esing_exc=zero
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
      call symmon
c---------------------------------------------------------------------
      call matrem('xmos')

      call retmem(5)
      call matremark
      call retmark
c---------------------------------------------------------------------
      call setival('ncore',ncore)
      call setival('nstrip2',nstrip2)
      call setival('mrcpf',mrcpf)
      call setival('lrcpf',lrcpf)
      call setrval('thrsx',thresh)
c---------------------------------------------------------------------
      call memory_status('end of mp2')
c---------------------------------------------------------------------
      end
c=======================================================================
      subroutine TransOneShell(ncf,  ncs,  nval,    ics,   kcs,
     *                       icf1,icf2, kcf1,kcf2,
     1                       ndisk,  inx,   xint,   CMO,   iprnt,
     2                       thresh, xmat, halftra, nrec,  mulam,
     3                       mulamd, nrow,  ncol,  irow,   icol,
     4                       irow1,  icol1, ENonZero,smal)

      use memory

      implicit real*8 (a-h,o-z)
c  transform a batch of integrals
c  Arguments:
c  INTEN=IN:
c  ncf  number of contracted basis functions
c  ncs  number of contracted shells
c  nval number of orbitals to be transformed (usually the valence ones)
c  ics,kcs: contracted shells for which the transformation is carried out
c  inx: array containing contraction info
c  xint: AO integrals ordered as
c        (l=1:ncf,j=1:ncf,k=kcf1:kcf2,i=icf1:icf2)   where icf1 & icf2
c        are the first and the last contr. funct. in the ICS shell ,
c        similarly for kcs .
c        NOTE THAT THEY ARE SCALED by 1/thresh
c  CMO:  MO coefficients
c  iprnt: print level
c  thresh = integral threshold
c  nrow   = number of non-zero rows of the AO integral matrices
c  ncol   = ditto for columns
c  irow(k),k=1,nrow = the indices of the non-zero rows
c  icol(k),k=1,nco  = ditto for columns
c  STORAGE:
c  xmat: (ncf,ncf) place for 1/4 transformed integrals
c  intent=OUT:
c  halftr: space for half-transformed integrals to be written on disk
c  nrec:  the number of records written  (input/output)
c  mulam: number of non-diagonal (mu, lambda) pairs transformed
c  integrals are {mu,nu|lam,isig);  mu.ne.lam
c  mulamd: number of diagonal pairs (mu,mu) transformed
c  ENonZero = number of the elements in non-zero columns and rows
c
      dimension xint(*),xmat(ncf,ncf),halftra(*),irow(ncf),icol(ncf)
      dimension CMO(ncf,*),irow1(ncf),icol1(ncf)
      dimension inx(12,*)
c     common /big/bl(30000)
      logical smal
c  return if there are no integrals
      if(nrow.eq.0.or.ncol.eq.0) RETURN
c
c  intstore holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
      call getint_4(nval**2,intstore)
      call getint_1(nval**2,int1)
c
c  icol and jcol store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
c  itrunc serves as temporary storage for a compacted coefficient matrix
c  only the rows (or columns in the second quarter transformation) of
c  X which are non-zero are  present in bl(itrunc)
c
      call getmem(ncf*nval,itrunc)
c
c     kcf1=inx(11,kcs)+1
c     kcf2=inx(10,kcs)
c     icf1=inx(11,ics)+1
c     icf2=inx(10,ics)
c
      do icf=icf1,icf2
         do kcf=kcf1,kcf2
            if(kcf.gt.icf) exit
            if (smal) then
               call ExtractOne(ncf,   icf,   kcf,   icf1,  icf2,
     2                         kcf1,  kcf2,  xint,  xmat,  nrow,
     3                         ncol,  irow,  icol,  nrow1, ncol1,
     4                         irow1, icol1)
            else
               call mmark
               call Extract1KI(ncf,   icf,   kcf,  icf1,  icf2,
     2                         kcf1,  kcf2,  xint, xmat,  nrow,
     3                         ncol,  irow,  icol, nrow1, ncol1,
     4                         irow1, icol1)
               call retmark
            endif
c  Determine some statistics:
            ENonZero=ENonZero+dble(nrow*ncol)
c
            if(nrow1.eq.0.or.ncol1.eq.0) cycle
            call IntTransp(ncf, nval,  icf,    kcf,      iprnt,
     1                     ndisk, nrow1, ncol1,  irow1,    icol1,
     2                     thresh,xmat,  CMO, bl(intstore),bl(int1),
     3                     itrunc,halftra,nrec)
            if (icf.ne.kcf) then
               mulam=mulam+1
            else
               mulamd=mulamd+1
            end if
         end do
      end do
      call retmem(3)
      end
c===========================================================================
      subroutine ExtractOne( ncf,   icf,   kcf,   icf1,  icf2,
     1                       kcf1,  kcf2,  xint,  xmat,  nrow,
     2                       ncol,  irow,  icol,  nrow1, ncol1,
     3                       irow1, icol1)
c  This routine extracts a single matrix from the integral array
c  The latter is indexed as xint(l=1:ncf,j=1:ncf,k=kcf1:kcf2,ncf,i=icf1:icf2);
c  elements (sigma,nu,kcf,icf) are put in xmat(nu,sigma)
c  It does not remove the overall scaling. The permutational
c  factors are not put in by the integrals program
c  Arguments:
c  INTENT: IN
c  ncf = number of contracted basis functions
c  icf, kcf = fixed AO indices of the first and third AOs
c  icf1,icf2 = The present shell goes from icf1 to icf2
c  kcf1,kcf2 = Same for kcf
c  xint = integral array
c  nrow = the number of non-zero rows
c  ncol = the number of non-zero columns
c  irow(1:ncf) = the first nrow elements of this array give the
c                original positions of the non-zero rows after
c                compacting
c  icol(1:ncf) = the same as above for the columns
c  INTENT:OUT
c  xmat = the exchange matrix extracted is put in xmat
c  nrow1, ncol1: the actual number of non-zero rows & columns
c  irow1(ncf),icol1(ncf): irow1(k) gives the original row number
c  of row k in the compacted matrix xmat
      implicit real*8 (a-h,o-z)
      dimension xint(ncf,ncf,kcf1:kcf2,icf1:icf2),xmat(ncf,ncf)
      dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
c
      do ll=1,ncol
         l=icol(ll)
         do jj=1,nrow
            j=irow(jj)
            xmat(jj,ll)=xint(l,j,kcf,icf)
         end do
      end do
c  matrix in xmat. Note that xmat is transmitted here as a parameter
c  but it is also a defined matrix
c  now compact xmat again: eliminate zero rows and columns
c  It is easier to do this for the columns
      ncol1=0
      do ll=1,ncol
         l=icol(ll)
         do jj=1,nrow
c  at this point, the integrals are scaled, so quantities below 1.0
c  are negligible
            if (abs(xmat(jj,ll)).gt.two) then
c column ll has at least one non-zero element
               ncol1=ncol1+1
               icol1(ncol1)=l
               go to 100
            endif
         end do
         go to 150
  100    continue
         if (ll.gt.ncol1) then
            do jj=1,nrow
               xmat(jj,ncol1)=xmat(jj,ll)
            end do
         end if
  150    continue
      end do
      nrow1=0
      do jj=1,nrow
         j=irow(jj)
         do ll=1,ncol1
            if (abs(xmat(jj,ll)).gt.two) then
               nrow1=nrow1+1
               irow1(nrow1)=j
               go to 200
            endif
         end do
         go to 250
  200    continue
         if (jj.gt.nrow1) then
            do ll=1,ncol1
               xmat(nrow1,ll)=xmat(jj,ll)
            end do
         end if
  250    continue
      end do
      end
c===========================================================================
      subroutine IntTransp(ncf,    nval,    mu,      lam,     iprnt,
     1                     ndisk,  nrow,    ncol,    irow,    icol,
     2                     thresh, xmat,    coef,    intint,  int1,
     3                     itrunc, halftra, nrec)
c   This routine transforms the integrals for a given mu, lam to MO basis,
c   yielding (mu,j|lam,l) where mu , lam are AO indices and j,l are MO
c   ones. It uses halftra both as a matrix 'halftra' and as an array
c   halftra for storing a temporary transformed matrix
c   The resulting integrals are written to disk as integers
c   ARGUMENTS:
c   INPUT
c  ncf: number of contracted basis functions
c  nval: number of orbitals to be transformed (valence orbitals)
c  mu,lam: fixed basis function indices. The transformed integrals are
c          (mu,i|lam,j) where i,j are MO indices
c  iprnt: print level
c  ndisk: the unit number of disk to write on
c  xmat = xmat(nu.isig)=(mu.nu | lam,isig) originally. However,
c         it is compacted here: xmat(p,q)=(mu,irow(p)|lam,icol(q))
c  note that xmat is also accessed by the matrix 'xmat'
c  coef(ncf,nval): SCF coefficients
c  nrow, ncol: the number of nonzero rows and columns of xmat
c  irow, icol = the indices of nonzero rows + columns
c  thresh = integral threshold
c  STORAGE
c  intint: storage for nval**2 4-byte integers
c  int1:  storage for nval**2 1-byte integers
c  itrunc = address of an (ncf x  nval) array. This holds a matrix
c  of the truncated SCF coefficients
c  OUTPUT
c  halftra: half-transformed integrals
c  nrec: number of records written so far (input/output)
c

      use memory

      implicit real*8 (a-h,o-z)
      integer*4 intint(nval,nval)
      integer*1 int1(nval,nval),i1
c -- masks needed for the modified 5-byte storage
      integer*1 m(7),mm
      dimension halftra(nval,nval),irow(nrow),icol(ncol)
      dimension xmat(ncf,ncf), coef(ncf,nval)
c     common /big/bl(30000)
c  the number below is the largest number to fit in a 4 byte signed integer
      parameter (zero=0.0d0,one=1.0d0,dblmax=2147483648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c  defined in the calling program
ccccccccc
cc      write(6,*) ' input xmat matrix is:'
cc      do i=1,nrow
cc      write(6,1111) (xmat(i,j)*thresh,j=1,ncol)
cc 1111 format(1x,6f12.6)
cc      enddo
cccccccccc
c  xmat has the COMPACTED AO exchange matrix
c    prepare the compacted SCF coefficients for the left side
      call CompactCoef(ncf,nval,coef,nrow,irow,bl(itrunc))
c  perform first matrix multiplication. Note that the matrix
c  system routine matsimtr cannot be used here because the right and
c  left matrices may be different after compaction
c  define temporary matrix product
cc      iout=igetival('iout')
      iout = 6     ! do NOT change this   JB
      call matdef('temp','r',nval,ncol)
c   'ymat' is similar to 'xmat' but is compacted
      call matdef('ymat','r',nrow,ncol)
      call PutTrunc(ncf,nrow,ncol,xmat,bl(mataddr('ymat')))
c  connect the truncated SCF coefficient matrix to the matrix system
      call matconn('trunc','r',nrow,nval,itrunc)
      call matmmul2('trunc','ymat','temp','t','n','n')
      call matdisc('trunc')
      call matrem('ymat')
      call CompactCoef(ncf,nval,coef,ncol,icol,bl(itrunc))
      call matconn('trunc','r',ncol,nval,itrunc)
      call matmmult('temp','trunc','halftra')
      call matdisc('trunc')
      call matrem('temp')
c
      if (iprnt.gt.3) then
         call matscal('halftra',thresh)
         write(iout,1236) mu,lam,nrow,ncol
 1236    format(' mu: ',I6,' lam: ',I6,' nrow: ',I6,' ncol: ',I6)
         call matprint('halftra',iout)
         call matscal('halftra',one/thresh)
      end if
c
c -- initialize the marks
cc      mm = 0
cc      m(1)=ibset(mm,5)
cc      m(2)=ibset(mm,6)
cc      m(3)=ibset(m(2),5)
cc      m(4)=ibset(mm,7)
cc      m(5)=ibset(m(4),5)
cc      m(6)=ibset(m(4),6)
cc      m(7)=ibset(m(6),5)
c
      do i=1,nval
         do j=1,nval
            x=halftra(j,i)
cc        call integerst(m,x,intint(j,i),int1(j,i))
            IF (abs(x).ge.dblmax) THEN
               b = x*dblinv
               if (abs(b).ge.d1max) then
                  dfac = abs(x)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
                  dfac = LOG10(dfac)
                  i1 = -NINT(dfac+0.5d0)
                  x = x*10.0d0**i1
                  int1(j,i) = i1
                  intint(j,i) = x
               else
                  i1 = abs(b)
                  b = x - SIGN(i1*dblmax,x)
                  int1(j,i)  = i1
                  intint(j,i) = b
               endif
            ELSE
               int1(j,i) = 0
               intint(j,i) = x
            ENDIF
         end do
      end do
c
      nrec=nrec+1
      write(unit=ndisk,rec=nrec) mu,lam,int1,intint
c
      return
      end
c==========================================================================
      subroutine TransVirt(ncf,i,j,iprnt,ipmij,ipairpr,virt)

      use memory

      implicit real*8 (a-h,o-z)
c  Arguments:
c  ncf:    # of AOS
c  i,j:    transform to virtual orbitals the pair (i,j)
c  iprnt: print level
c  ipmij:  flag (>0) for printing selected AO exchange matrices (test)
c  ipairpr:pair to be printed
c
      character*(*) virt
c
      if(iprnt.gt.3) then
        write(6,1234) i,j
 1234   format( ' Pair: ',2I6)
        call matprint('xmat',6)
      end if
c
      if(ipmij.ge.1) then
        ij=i*(i-1)/2+j
        if(ij.eq.ipairpr) then
          write(6,1235) i,j
 1235     format(' AO exchange matrix for pair: ',2I6)
          call matprint('xmat',6)
        end if
      end if
c
      call matsimtr('xmat',virt,'xmos')
c
      if(iprnt.gt.3) then
        call matprint('xmos',6)
      end if
      end
c=====================================================================
      subroutine PairEner(ncf,    nmo,    nvirt,  i,      j,
     $                    x,      e,      eij,    tnorm,  eijs,
     $                    eijt,   eijpar, eijanti)
      implicit real*8 (a-h,o-z)
c  Arguments (INPUT):
c  ncf:            # of AO basis functions
c  nvirt:          # of virtual orbitals
c  i,j    :        occupied orbitals defining the pair; e(i,j) is calc.
c  x(nvirt,nvirt): x(a,b)=(a i|j b) MO integrals in matrix form
c  e:              virtual orbital energies
c       OUTPUT:
c  eij:            pair energy
c  eijs:           singlet coupled
c  eijt:           triplet coupled
c  eijpar:         parallel spin contribution for scaled MP2
c  eijanti:        antiparallel spin contribution for scaled MP2
      dimension x(nvirt,nvirt),e(ncf)
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0,three=3.0d0)
      ei=e(i)
      ej=e(j)
      eij=zero
      eijs=zero
      eijt=zero
      eijpar=zero
      eijanti=zero
      tnorm=zero
      do ia=1,nvirt
        ea=e(nmo+ia)
        do ib=1,nvirt
          rdenom=one/(ei+ej-ea-e(nmo+ib))
          tij= x(ia,ib)*rdenom
          tij1=(x(ia,ib)+x(ia,ib)-x(ib,ia))*rdenom
          tnorm=tnorm+tij*tij1
          eij=eij+x(ia,ib)*tij1
c -- separate singlet and triplet contributions
c    singlet part divided by two if i=j
c    triplet part multiplied by three
          eijs=eijs+tij*(x(ia,ib)+x(ib,ia))
          eijt=eijt+tij*(x(ia,ib)-x(ib,ia))
c -- Grimme's parallel and antiparallel components
          eijpar=eijpar+tij*(x(ia,ib)-x(ib,ia))
          eijanti=eijanti+tij*x(ia,ib)
        end do
      end do
      eijs=half*eijs
      eijt=three*eijt
      eijpar=eijpar+eijpar
      if(i.ne.j) then
        eij=eij+eij
        tnorm=tnorm+tnorm
        eijs=eijs+eijs
        eijanti=eijanti+eijanti
      end if
      end
c=====================================================================
      subroutine SaveTij(ncf,  nmo,  nval,  nvir,  i,
     1                   j,    e,    thresh,intres,int1,
     2                   x,   ndisk, ij)
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
      dimension e(*)
      dimension x(nvir,nvir)
      integer*1 int1(nvir,nvir),i1
      integer*4 intres(nvir,nvir)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      parameter (dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c
      end1=e(i)+e(j)
      do ib=1,nvir
         end2=end1-e(nmo+ib)
         do ia=1,nvir
            x(ia,ib)=x(ia,ib)/(end2-e(nmo+ia))
C NZG
C            write(6,"(4I5,F15.8)") i,j,ia,ib,x(ia,ib)
         end do
      end do
C  check magnitude and convert to integers
      call matscal('xmos',one/thresh)
      do ib=1,nvir
         do ia=1,nvir
            val = x(ia,ib)
            If (abs(val).ge.dblmax) Then
               b = val*dblinv
               if (abs(b).ge.d1max) then
                  dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
                  dfac = LOG10(dfac)
                  i1 = -NINT(dfac+0.5d0)
                  val = val*10.0d0**i1
                  int1(ia,ib) = i1
                  intres(ia,ib) = val
cc          write(6,*) ' Threshold Tij - ia:',ia,' ib:',ib,' i1:',i1
               else
                  i1 = abs(b)
                  b = val - SIGN(i1*dblmax,val)
                  int1(ia,ib) = i1
                  intres(ia,ib) = b
               endif
            Else
               int1(ia,ib) = 0
               intres(ia,ib) = val
            EndIf
         enddo
      enddo
C  write final matrix to disk unit ndisk
      write(ndisk,rec=ij) int1,intres
c
      return
      end
C========SaveKij=====================================================
      subroutine saveKij(ncf,nval,nvir,x,intres,int1,
     1                   ij,thresh,length,ndisk,ncore)
      implicit real*8(a-h,o-z)
      dimension x(*)
      integer*1 int1(length,2),i1
      integer*4 intres(length,2)
      parameter (one=1.0d0,dblmax=2147483648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
c
      call matpose('xmat')
      call matdef('tss02','r',nvir,ncf)
      call matmmult('tvir','xmat','tss02')
      if (ncore.gt.0) then
         call matmmult('tss02','occa','kijvo')
      else
         call matmmult('tss02','occu','kijvo')
      endif
C  check magnitude and convert to integers
      call matscal('kijvo',one/thresh)
      do imv=1,length
         val = x(imv)
         If (abs(val).ge.dblmax) Then
            b = val*dblinv
            if(abs(b).ge.d1max) then
               dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
               dfac = LOG10(dfac)
               i1 = -NINT(dfac+0.5d0)
               val = val*10.0d0**i1
               int1(imv,2) = i1
               intres(imv,2) = val
cc          write(6,*) ' Threshold Kij - imv:',imv,' i1:',i1
            else
               i1 = abs(b)
               b = val - SIGN(i1*dblmax,val)
               int1(imv,2) = i1
               intres(imv,2) = b
            endif
         Else
            int1(imv,2) = 0
            intres(imv,2) = val
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
         If (abs(val).ge.dblmax) Then
            b = val*dblinv
            if (abs(b).ge.d1max) then
               dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
               dfac = LOG10(dfac)
               i1 = -NINT(dfac+0.5d0)
               val = val*10.0d0**i1
               int1(imv,1) = i1
               intres(imv,1) = val
cc          write(6,*) ' Threshold Kij - imv:',imv,' i1:',i1
            else
               i1 = abs(b)
               b = val - SIGN(i1*dblmax,val)
               int1(imv,1) = i1
               intres(imv,1) = b
            endif
         Else
            int1(imv,1) = 0
            intres(imv,1) = val
         EndIf
      enddo
      call matrem('tss02')
C  finally write to disk
      write(ndisk,rec=ij) int1,intres
      end
c=====================================================================
      subroutine BinSort(ncf,   nval,   npairs, ndisk,   ndisk2,
     1                   nrec,  lbin,   nbins,  igrfirst, igrlast,
     2                   ibins, ibin1,  i1pair, int1,   icounter,
     3                  ibinctr,lrec,igran2pair,ipair2gran,igrsize)
      implicit real*8 (a-h,o-z)
      integer*1 int1(nval,nval),ibin1(lbin,npairs)
      integer*2 iindex(2),jindex(2)
      integer*4 ibins(2*lbin,npairs),i1pair(nval,nval)
      integer*4 igran2pair(2,*),ipair2gran(*)
      dimension icounter(npairs), ibinctr(nbins)
      integer*4 longint,longint1
      equivalence(longint,iindex(1)),(longint1,jindex(1))
c  This routine sorts the half-transformed integrals
c  Arguments
c  INPUT
c  ncf:    number of AO basis functions
c  nval:   number of occupied MOs to be transformed (usually the
c             valence ones)
c  npairs: number of occupied orbital pairs, nval*(nval+1)/2
c  ndisk:  unit number of file containing half transformed integrals
c  ndisk2: unit number of file containing sorted integrals
c  nrec:   number of records written on the half-transformed integral file
c  lbin:   number of records per bin
c  nbins:  estimate of the number of bins which will be produced
c  ijfirst:first pair to be sorted
c  ijlast: last pair to be sorted
c    the pair number (ij) is i*(i-1)/2+j. Integrals (mu,i|j,nu)
c    where ijfirst<=(ij)<=ijlast will be sorted
c  ibins:  integer*4 storage for the bins
c  ibin1:  integer*1 storage for precision overflow
c    (how many times 2**31 divides integer representation of integral)
c  i1pair: integer*4 storage are for the exchange matrix sorted in
c    the wrong order  (mu,i|j,lam) for all mu>=lam, all i,j
c  int1:   integer*1 storage for precision overflow
c    (how many times 2**31 divides integer representation of integral)
c  icounter: counters showing the degree of occupancy of the bins
c  ibinctr: there is one integer for each bin WRITTEN OUT, showing
c           which pair it contains
c  lrec:     length of each bin on direct access file
c
      call izeroit(icounter,npairs)
      nbin=0                ! no. of bins written so far
      iirec=0
      ijfirst=igran2pair(1,igrfirst)
      ijlast =igran2pair(2,igrlast)
      isize=igrsize*lbin
cc      write(6,'(A,2I5)') "Binsort: ",ijfirst,ijlast
c
      do irec=1,nrec
        iirec=iirec+1
        read(ndisk,rec=iirec) mu,lam,int1,i1pair
c  end of file: mu=lam=-1
        if(mu.eq.-1.and.lam.eq.-1) go to 300
        iindex(1)=mu
        iindex(2)=lam
        jindex(1)=lam
        jindex(2)=mu
        ij=0
        do i=1,nval
        do j=1,i
          ij=ij+1
          if(ij.lt.ijfirst) cycle
          if(ij.gt.ijlast) cycle
c -- test on integral  (don't store if zero)
          If(i1pair(i,j).ne.0.OR.int1(i,j).ne.0) Then
c
            icounter(ij)=icounter(ij)+2
            ibins(icounter(ij)-1,ij)=i1pair(i,j)
            ibins(icounter(ij),ij)=longint
            ibin1(icounter(ij)/2,ij)=int1(i,j)
            if(icounter(ij).ge.2*lbin) then
              igran=ipair2gran(ij)
              ist=igran2pair(1,igran)
              istop=igran2pair(2,igran)
              nbin=nbin+1
              ibinctr(nbin)=igran
              do iii=ist,istop
                do k=icounter(iii)/2+1,lbin
                  ibins(2*k-1,iii)=0
                  ibins(2*k,iii)=0
                  ibin1(k,iii)=0
                end do
                icounter(iii)=0
              enddo
cc              write(6,*) ' About to call <WriteBin>   isize:',isize
              call WriteBin(ndisk2,nbin,isize,ibins(1,ist),ibin1(1,ist))
            end if
          EndIf
c
          IF(mu.ne.lam) THEN
c
c -- test on integral  (don't store if zero)
            If(i1pair(j,i).ne.0.OR.int1(j,i).ne.0) Then
c
              icounter(ij)=icounter(ij)+2
              ibins(icounter(ij)-1,ij)=i1pair(j,i)
              ibins(icounter(ij),ij)=longint1
              ibin1(icounter(ij)/2,ij)=int1(j,i)
              if(icounter(ij).ge.2*lbin) then
                igran=ipair2gran(ij)
                ist=igran2pair(1,igran)
                istop=igran2pair(2,igran)
                nbin=nbin+1
                ibinctr(nbin)=igran
                do iii=ist,istop
                  do k=icounter(iii)/2+1,lbin
                    ibins(2*k-1,iii)=0
                    ibins(2*k,iii)=0
                    ibin1(k,iii)=0
                  end do
                  icounter(iii)=0
                enddo
cc              write(6,*) ' About to call <WriteBin> 2   isize:',isize
              call WriteBin(ndisk2,nbin,isize,ibins(1,ist),ibin1(1,ist))
              end if
            EndIf
          ENDIF   !  mu=lam condition
        end do    !  j loop
        end do    !  i loop
c
      end do      !  irec loop
  300 continue
c  write out the content of the bins at the end
      do ij=1,npairs
        if(icounter(ij).gt.0) then
          igran=ipair2gran(ij)
          ist=igran2pair(1,igran)
          istop=igran2pair(2,igran)
          nbin=nbin+1
          ibinctr(nbin)=igran
          do iii=ist,istop
            do j=icounter(iii)/2+1,lbin
              ibins(2*j-1,iii)=0
              ibins(2*j,iii)=0
              ibin1(j,iii)=0
            end do
            icounter(iii)=0
          enddo
cc        write(6,*) ' About to call <WriteBin> 3   isize:',isize
        call WriteBin(ndisk2,nbin,isize,ibins(1,ist),ibin1(1,ist))
        end if
      end do
      nbins=nbin
      end
c==============================================================
      subroutine WriteBin(ndisk2,nbin,lbin,ibins,ibin1)
c  This routine writes a full bin on disk
c  Arguments:
c  ndisk2: unit number
c  nbin:   position of the bin to be written
c          so far, bins are written as they come. One could write
c          them in a better order
c  lbin:   number of records per bin
c  ibins:  compressed integral + indices
c  ibin1:  precision overflow integer
c
      integer*1 ibin1(lbin)
      integer*4 ibins(2*lbin)
cccccccccc
cc      integer*2 iindex(2)
cc      integer*4 longint
cc      equivalence (longint,iindex(1))
cc      write(6,*) ' In <WriteBin>  lbin:',lbin,'  Integral indices:'
cc      do i=1,lbin
cc      longint = ibins(2*i)
cc      mu = iindex(1)
cc      lam = iindex(2)
cc      write(6,1234) i,mu,lam
cc 1234 format(' integral ',i4,'  indices: ',2I4)
cc      enddo
cccccccccc
      write(ndisk2,rec=nbin) ibin1,ibins
      end
c==============================================================
      subroutine ReadBin1(ncf,   iorbpair,      igran,  ndisk2,
     1                   lbin,   nbins, ibinctr,thresh, ibins,
     2                   ibin1,  nsym,  ifunpair,xmat,  igran2pair,
     3                   ipair2gran,igrsize)
      implicit real*8 (a-h,o-z)
      integer*1 ibin1(lbin,0:igrsize-1)
      integer*2 iindex(2)
      integer*4 ibins(2,lbin,0:igrsize-1)
      integer*4 igran2pair(2,*),ipair2gran(*)
      dimension ipair(7),jpair(7),ifunpair(7,ncf)
      dimension iorbpair(nsym,*)
      dimension xmat(ncf,ncf,0:igrsize-1),ibinctr(nbins)
      integer*4 longint
      equivalence (longint,iindex(1))
      Parameter (dblmax=2 147 483 648.0d0)
c
c  This routine reads all bins belonging to pair ij
c  in this first implementation, all records in ibinctr are
c  scanned for contributions to pair (ij); this should be replaced
c  by a proper sort of the indices ibinctr prior to reading.
c      ARGUMENTS:
c  ncf =    number of AO basis functions
c  ipair(7),jpair(7) symmetry pairs (it should be only + or - i or
c  j, resp.) of the orbitals i and j consituting the pair (ij) under
c  the symm. operations
c  ij =     pair under consideration
c  ndisk2 = unit number of the disk file containing the sorted integrals
c  lbin  =  number of records per bin
c  nbins =  nymber of bins
c  ibinctr(ibin) shows which pair does  bin ibin belong to
c  thresh = integral threshold, needed to remove scaling of the integrals
c  ibins =  compressed integral and indices
c  ibin1 =  integer*1 storage for precision overflow
c    (how many times 2**31 divides integer representation of integral)
c  xmat =   external exchange matrix X(mu,lam) generated
c  nsym = number of symmetry operations
c  ifunpair = an array showing the symmetry pairs of the AOs
c    ifunpair(7,ncf)
c
      call zeroit(xmat,ncf*ncf*igrsize)
      dblcmp = dblmax*thresh
c
      ijstart=igran2pair(1,igran)
      ijstop =igran2pair(2,igran)
      do kk=1,nbins
      if(ibinctr(kk).ne.igran) cycle
        read(ndisk2,rec=kk) ibin1,ibins
          do ij=ijstart,ijstop
            ijind=ij-ijstart
cc            write(6,*) '-- IJ index is: ',ijind,' --'
            call get_ij_half(ij,i,j)
            do isym=1,nsym
              ipair(isym)=iorbpair(isym,i)
            end do
            do isym=1,nsym
              jpair(isym)=iorbpair(isym,j)
            end do
            do jj=1,lbin
            longint=ibins(2,jj,ijind)
            mu=iindex(1)
            lam=iindex(2)
            if(mu.eq.0.or.lam.eq.0) cycle
            if(mu.lt.0.or.mu.gt.ncf.or.lam.lt.0.or.lam.gt.ncf) then
              call nerror(1,'ReadBin1','wrong index', mu,ncf)
            else
c -- decompress the integral -------------------------
cc              call reconstruct(ibins(1,jj,ijind),ibin1(jj,ijind),xx)
cc              xx = xx*thresh
              If(ibin1(jj,ijind).eq.0) Then
                xx = ibins(1,jj,ijind)*thresh
              Else If(ibin1(jj,ijind).gt.0) Then
                xx = ibins(1,jj,ijind)*thresh
                xx = xx + SIGN(ibin1(jj,ijind)*dblcmp,xx)
              Else
                xx = ibins(1,jj,ijind)*thresh*10.0d0**(-ibin1(jj,ijind))
              EndIf
c -----------------------------------------------------
              xmat(mu,lam,ijind)=xx
              do isym=1,nsym
                xxx=xx
                mu1=ifunpair(isym,mu)
                mu2=abs(mu1)
                if(mu1.lt.0) xxx = -xxx
                lam1=ifunpair(isym,lam)
                lam2=abs(lam1)
                if(mu2.eq.mu.and.lam2.eq.lam) cycle
                if(lam1.lt.0) xxx = -xxx
                if(ipair(isym).lt.0) xxx = -xxx
                if(jpair(isym).lt.0) xxx = -xxx
                xmat(mu2,lam2,ijind)=xxx
              end do
cc              write(6,*) '  X(',mu,',',lam,') = ',xmat(mu,lam,ijind)
            end if
          end do  ! jj elements in the bin
        enddo
      end do    ! kk; bins
      end
c===================================================================
      subroutine SortTransf(
     1   ncf,         nval,       npairs,     ndisk,      ndisk2,
     1   nrec,        lbin,       nbins,      igrfirst,   igrlast,
     2   iout,       emp2only,    virt,       i1pair,     int1,
     3   icounter,    ibinctr,    lrec,       nfirst,     nmo,
     4   nvirt,       thresh,     nsym,       iorbpair,   ifunpair,
     5                iprnt,      ipmij,      ipairpr,    xmos,
     6   epsi,        emp2,       tnorm,      emps,       empt,
     7   empanti,     emppar,     TBinSort,   ndisk3,     krec,
     8   tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     9   ttkij,       tmp2e,      ncore, igran2pair,ipair2gran,
     $   igrsize)

      use memory

      implicit real*8 (a-h,o-z)
      integer*4 i1pair(*)
      dimension icounter(*),ibinctr(*)
      dimension epsi(*)
      dimension iorbpair(nsym,nfirst:nmo),ifunpair(7,ncf)
      dimension ipair(7),jpair(7)
      integer*4 igran2pair(2,*),ipair2gran(*)
      integer*1 int1(*)
      logical emp2only
      character*(*) virt
      parameter(zero=0.0d0,one=1.0d0,four=4.0d0)
c
c This routine has a number of arguments because it transmits
c them to BinSort, ReadBin, TransVirt and PairEner
c Parameters 1-20 and 37 belong to BinSort
c Parameters 21-28 belong to ReadBin
c Parameters 29-31 belong to TransVirt
c Parameters 32,33 belong to PairEner
c Parameters 34-36 are output:
c  emp2 = MP2 energy
c  tnorm = total norm of the wavefunction in intermediate normalization
c  TBinSort = elapsed time in BinSort
C
      call mmark()
c  reserve memory for bins
      call getint_4(npairs*lbin*2,ibins)
      call getint_1(npairs*lbin,ibin1)
c
      call elapsec(elap1)
      call BinSort(ncf,      nval,     npairs, ndisk,    ndisk2,
     1             nrec,     lbin,     nbins,  igrfirst,  igrlast,
     2             bl(ibins),bl(ibin1),i1pair, int1,    icounter,
     3            ibinctr,lrec,igran2pair,ipair2gran,igrsize)
      call elapsec(elap2)
      TBinSort=TBinSort+elap2-elap1
      call retmark()
      call mmark()
      call getmem(ncf*ncf*igrsize,ixmat)
      call getint_4(igrsize*lbin*2,ibins)
      call getint_1(igrsize*lbin,ibin1)
c
c -- Additional stuff for RMP2 gradients
      If(.NOT.emp2only) Then
        nval=nmo-nfirst+1
        nvmo=nval
        if(ncore.gt.0) nvmo=nmo
        call matdef('kijvo','r',nvirt,nvmo)
        ikijvo=mataddr('kijvo')
        lengvo=nvirt*nvmo
      EndIf
c
c  read back bins, extract 1 exchange matrix, and transform that
      ij=0
      do igran=igrfirst,igrlast
      call ReadBin1(ncf,       iorbpair,       igran,  ndisk2,
     1              lbin,      nbins, ibinctr, thresh, bl(ibins),
     2              bl(ibin1), nsym,  ifunpair, bl(ixmat),igran2pair,
     3              ipair2gran,igrsize)
c
c-------------------------------------------------------------------
      ijstart=igran2pair(1,igran)
      ijstop =igran2pair(2,igran)
      do ij=ijstart,ijstop
        ijind=ij-ijstart
        call get_ij_half(ij,i,j)
        i=i+nfirst-1
        j=j+nfirst-1
        call matconn('xmat','q',ncf,ncf,ixmat+ijind*ncf*ncf)
c
c -- Additional stuff for MP2 gradients
      If(.NOT.emp2only) Then
c -- Kij in AO basis now in xmat
c -- transform to virtual- occupied block and save for gradients
        call elapsec(tkij1)
        mxrec=mxrec+1
        call SaveKij(ncf,      nvmo,   nvirt,bl(ikijvo),bl(ibins),
     1               bl(ibin1),ij,     thresh, lengvo,  ndisk5,
     2               ncore)
        call elapsec(tkij2)
        ttkij=ttkij+tkij2-tkij1
      EndIf
c--------------------------------------------------------------------
c
      call TransVirt(ncf,i,j,iprnt,ipmij,ipairpr,virt)
      call PairEner(ncf,    nmo,    nvirt,  i,      j,
     $              xmos,   epsi,   eij,    xnorm,  eijs,
     $              eijt,   eijpar, eijanti)
      emp2=emp2+eij
      tnorm=tnorm+xnorm
c -- Singlet and triplet components of the correlation energy
      emps=emps+eijs
      empt=empt+eijt
      emppar=emppar+eijpar
      empanti=empanti+eijanti
c
c-------------------------------------------------------------------
c -- Additional stuff for MP2 gradients
      If(.NOT.emp2only) Then
c -- for MP2-gradients, backtransform to AO basis after
c -- scaling with energy denominators
        call elapsec(tt1)
        tmp2e=tmp2e+tt1-tkij2
        nxrec=nxrec+1
        call SaveTij(ncf,    nmo,    nval,   nvirt,  i,
     1               j,      epsi,   thresh, bl(ibins),bl(ibin1),
     2               xmos,   ndisk3, ij)
        call elapsec(tt2)
        tmbtr=tmbtr+tt2-tt1
      EndIf
c--------------------------------------------------------------------
c
      if(iprnt.gt.2) write(iout,100) i,j,ij,eij,xnorm
  100 format(2i4,i6,f12.7,2x,f12.7)
        call matdisc('xmat')
      end do
      end do
      if(.NOT.emp2only) call matrem('kijvo')
      call retmark()
      end
c=================================================================
      subroutine OrbSymPair(nsym, ncf ,coef, ifp ,ifirst,
     1                   ilast,iprnt, tmp, iorbpair,iout)
c  this routine determines the symmetry pairs of orbitals
c  under the Abelian point group
c
c  INPUT:
c  nsym=number of symmetry operations, max. 7
c  ncf =number of contracted functions
c  coef=orbital coefficient matrix, first index basis fnct, second orbital
c  ifp(7,ncf)= basis function symmetry pairs. abs(ifp(ns,k)) gives
c  the image of basis function k under symmetry operation ns.
c  This number is negative if k changes sign under the symm. operation
c  ifirst,ilast: the first and last orbital to be symmetry-analyzed
c
c  STORAGE:
c  tmp(ncf)
c  OUTPUT
c  iorbpair(nsym,nfirst:nlast) :  iorbpair(isym,i) is the symmetry
c   image of orbital i . It is negative if the symm. image is neg.
c
c  the overlap matrix is not used here, i.e. the AOs are assumed
c  to be orthonormal. This is OK as long as it is not assumed that
c  they are also normalized
      implicit real*8 (a-h,o-z)
      logical breaksymm
      dimension coef(ncf,*),tmp(ncf)
      dimension ifp(7,ncf),iorbpair(nsym,ifirst:ilast)
      parameter(zero=0.0d0,one=1.0d0,eps=2.0d-8)
c
      if(nsym.eq.0) return
      breaksymm=.false.
      do iorb=ifirst,ilast
c       print *,'iorb,ifirst,ilast',iorb,ifirst,ilast
c       Determine the "norm" of the orbital (assuming orthonormal AOs
        xnorm=ddot(ncf,coef(1,iorb),1,coef(1,iorb),1)
        do ns=1,nsym
          iorbpair(ns,iorb)=0
          do k=1,ncf
c generate the symmetry image of orbital k under symmetry op. ns
            ispr=ifp(ns,k)
            isp1=abs(ispr)
            if(ispr.gt.0) then
              tmp(isp1)=coef(k,iorb)
            else if(ispr.lt.0) then
              tmp(isp1)=-coef(k,iorb)
            else
              tmp(isp1)=zero
            end if
          end do     !k
c  Calculate the "overlap" between the symmetry image of the orbital
c   and the other orbitals
          do jorb=ifirst,ilast
            ovrl= ddot(ncf,tmp,1,coef(1,jorb),1)
            if(abs(ovrl-xnorm).lt.eps) then
              iorbpair(ns,iorb)=jorb
              exit
            else if(abs(ovrl+xnorm).lt.eps) then
              iorbpair(ns,iorb)=-jorb
              exit
            end if
          end do   !jorb
          if(iprnt.gt.3) then
      write(iout,*)'ns,iorb,iorbpair(ns,iorb)',ns,iorb,iorbpair(ns,iorb)
          end if
          if(iorbpair(ns,iorb).eq.0) then
            breaksymm=.true.
          end if
        end do     !ns
      end do       !iorb
      if(breaksymm) call nerror(1,'OrbSymPair',
     1  'Orbitals do not conform to symmetry',iorb,jorb)
      end
c=======================================================
      logical function LastSymPair(ics,kcs,nsym,ishpr,inegl,iret)
c  This routine returns a value of true if there is a shell pair
c  symmetry equivalent to (ics,kcs) which is behind (ics,kcs) in
c the canonical order
c  Arguments:
c  INPUT:
c  ics,kcs = contracted shells
c  nsym = number of symmetry operations
c  ishpr(7,ncs) = ishpr(isym,ics) is the symmetry equivalent of shell
c                 ics under symmetry operation isym
c  OUTPUT:
c  LastSymPair = The last symm. equivalent of ics
c  inegl = it is increased by 1 for each pair omitted
c  iret = it is increased by 1 for each pair retained
      dimension ishpr(7,*)
      LastSymPair=.false.
      do isym=1,nsym
        iprime=ishpr(isym,ics)
        kprime=ishpr(isym,kcs)
        if(iprime.gt.ics.or.kprime.gt.ics) then
          LastSymPair=.true.
          go to 100
        else
          ii=max(iprime,kprime)
          kk=min(iprime,kprime)
        end if
        if(ii.gt.ics.or.ii.eq.ics.and.kk.gt.kcs) then
          LastSymPair=.true.
          go to 100
        end if
      end do
      iret=iret+1
c      print *,'pair ',ics,kcs,' retained'
      return
100   continue
      inegl=inegl+1
c      print *,'pair ',ics,kcs,' omitted'
      end
c=======================================================
      subroutine DmxMakeL(dualbasis, nval, nmo, ncf, ncs,
     $                    np4, beta, inx, iprnt,DMax)
c...................................................................
C      Main routine for construction DMax(ncs,ncs) matrix
C      to be used for integral pre-screening for MP2
C      with Canonical orbitals.
C      Dmax(m,l) = max(i,j){abs(L(m,i)*L(l,j)} * rfact(i,j)
C      Where Rfac(i,j) = sqrt(ri*rj)/(rij**3) * 1/(-ei-ej)
C      the factor rfac should be equivalent to eliminating
C      contribution from negligable pairs.
C      Input:
c      dualbasis - logical
C      nval number of valence orbital
C      nmo number of occupied orbitals
C      ncf number of contracted basis functions
C      ncs number of contracted shells
C      np4 unit number , normally 4
C      beta - logical flag for beta spin (otherwise alpha/closed shell)
C      inx array with shell informations
C      Output:
C      DMax - the final dmax matrix
C
C    Svein Saebo  Fayetteville AR, June 1999
c    modified by KW for dual basis set approach
C

      use memory

      implicit real*8(a-h,o-z)
      logical dualbasis,beta
ckw
      character*256 jobname,MOS
      common /job/jobname,lenJ
ckw
      dimension DMax(ncs,ncs),inx(*)
c
      ncore=nmo-nval
c
      call matdef('dicf','q',ncf,ncf)
      call matdef('rfac','q',nval,nval)
c
      idicf=mataddr('dicf')
      irfa =mataddr('rfac')
c
C get orbital centers and radii form disk (written out in LOCA)
c
      call matdef('xmo','d',nmo,nmo)
      call matdef('ymo','d',nmo,nmo)
      call matdef('zmo','d',nmo,nmo)
      call matdef('r2mo','d',nmo,nmo)
      if(.NOT.beta) then
cc        write(6,*) ' In <DmxMakeL>  Reading Alpha localized centroids',
cc     $             ' - nmo is ',nmo
        call matread('xmo',np4,'xcor_loA')
        call matread('ymo',np4,'ycor_loA')
        call matread('zmo',np4,'zcor_loA')
        call matread('r2mo',np4,'rr_loA')
      else
cc        write(6,*) ' In <DmxMakeL>  Reading Beta localized centroids',
cc     $             ' - nmo is ',nmo
        call matread('xmo',np4,'xcor_loB')
        call matread('ymo',np4,'ycor_loB')
        call matread('zmo',np4,'zcor_loB')
        call matread('r2mo',np4,'rr_loB')
      endif
c
      ixmo=mataddr('xmo')+ncore
      iymo=mataddr('ymo')+ncore
      izmo=mataddr('zmo')+ncore
      irra=mataddr('r2mo')+ncore
      call matdef('eval','d',nmo,nmo)
      if(.NOT.beta) then
        call matread('eval',np4,'eloc_rhf')
      else
        call matread('eval',np4,'eloc_uhf')
      endif
      ieps=mataddr('eval')+ncore
c
C     make rfact(i,j):
c
      call rfactor(bl(ixmo),bl(iymo),bl(izmo),bl(irra),nval,bl(irfa),
     *             bl(ieps))
c
C      rfac(i,j) is in bl(irfa)
c
      call matrem('eval')
      call matrem('r2mo')
      call matrem('zmo')
      call matrem('ymo')
      call matrem('xmo')
c
      call matdef('loca','r',ncf,nmo)
c
      if(.not.dualbasis) then
        if(.NOT.beta) then
          call matread('loca',np4,'loca_rhf')
        else
          call matread('loca',np4,'loca_uhf')
        endif
      else
c        get eigenvectors from MOS file :
c        they are supposed to be already transformed
c        from small to big basis set by the guess module
c.....
         call tstchval('mos-file',iyes)
         If(iyes.eq.1) Then
           call getchval('mos-file',MOS)
         Else
           MOS = jobname(1:lenJ)//'.mos'
         EndIf
c.....
c
         iloca=mataddr('loca')
         call rmblan2(MOS,256,lenM)
         itype=1
         call readmos(ncf,bl(iloca),jnk,.False.,lenM,MOS,itype,IErr)
cccc     call matprint ('loca',6)   !  kw test
      endif
c
      call matsub('occo','loca',ncore+1,nmo)
      ioccu=mataddr('occo')
c
      call getmem(ncf*nval,itemp)
      call Dmxij_new(bl(idicf), DMax,  ncf,      ncs,      nval,
     1               inx,   bl(ioccu), bl(irfa), bl(itemp),iprnt)
      call retmem(1)
c-----------------------------------------------------------------
      if(dualbasis) then
         call getival('mpres_in_bg',mpres_in_bg)
         call updateDS(ncs, DMax    , bl(mpres_in_bg))
      endif
c-----------------------------------------------------------------
      call matrem('occo')
      call matrem('loca')
      call matrem('rfac')
      call matrem('dicf')
c
      end
C=====================================================================
      subroutine DmxIJ(D,  DS, ncf,   ncs,   nval,
     1                 inx,C,  Rfact, iprnt)
C      Calculates Max  of C(my,i) * C(lam,j)
C      Used for prescreening of integrals.
C      Output:
C      D(ncf,ncf) (this matrix is removed after return)
C       DS(ncs,ncs) same as D but over shells rather than
C      individual functions
C      Input:
C      ncf - number of basis functions
C      ncs - number of contracted shells
C      nval - number of valence orbitals
C      C(ncf,nval) matrix of MO coefficients (localized orbitals)
C      Svein Saebo, Fayetteville AR June 1999
c      this routine is awfully inefficient. Use its accelerated version
c      dmxij_new

      use memory

      implicit real*8(a-h,o-z)
      dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs)
      dimension rfact(nval,nval)
      dimension inx(*)
      do my=1,ncf
        do lam=1,my
        dmxx=0.0d0
          do imo=1,nval
            do jmo=1,nval
              cprod=abs(C(my,imo)*C(lam,jmo))
              if(imo.ne.jmo) cprod=cprod*rfact(imo,jmo)
              dmxx=max(dmxx,cprod)
            end do
          end do
          D(my,lam)=dmxx
          D(lam,my)=dmxx
        end do
      end do
C      call evaldmax(ncf,D,iprnt)
C      call prntmat(ncf,ncf,ncf,D)
C
C  construct a matrix over contracted shells rather than basis function
C
      ictr=igetival('ictr')
      ipoint=0
      do ics=1,ncs
        call get_shell_size(bl(ictr),ics,isize)
        kpoint=0
        do kcs=1,ncs
          call get_shell_size(bl(ictr),kcs,ksize)
          dmxs=0.0d0
          do my=1,isize
            do lam=1,ksize
            dmx=D(ipoint+my,kpoint+lam)
            dmxs=max(dmxs,dmx)
          end do     !lam in kcs
        end do     ! my in ics
        kpoint=kpoint+ksize
ckw
        dmxs=min(dmxs,1.d0)
ckw
ckw     DS(ics,kcs)=dmxs
        DS(ics,kcs)=dmxs*dmxs
c
      end do     ! kcs
      ipoint=ipoint+isize
      end do     !ics
      end
C=====================================================================
      subroutine DmxIJ_new(D,  DS, ncf,   ncs,   nval,
     1                 inx,C,  Rfact,     dmaxi,  iprnt)
C      Calculates Max  of C(my,i) * C(lam,j)
C      Used for prescreening of integrals.
C      Output:
C      D(ncf,ncf) (this matrix is removed after return)
c        D(mu,lam)=max(i,j)| C(mu,i)*C(lam,j)*rfact(i,j)|
C      DS(ncs,ncs) same as D but over shells rather than
C      individual functions
C      Input:
C      ncf - number of basis functions
C      ncs - number of contracted shells
C      nval - number of valence orbitals
C      C(ncf,nval) matrix of MO coefficients (localized orbitals)
c      Rfact is a factor estimating the interaction between
c        two localized orbitals, i and j (see P. Pulay, S. Saebo,
c        K. Wolinski, Chem. Phys. Lett. 2001)
c      dmaxi = intermediate array, ncf*nval
C      original version Svein Saebo, Fayetteville AR June 1999
c     accelerated version (about 100 times faster for large molecules:
c     P. Pulay, Spring 2002
c
      use memory

      implicit real*8(a-h,o-z)
      dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs)
      dimension rfact(nval,nval),dmaxi(ncf,nval)
      dimension inx(*)
      parameter (zero=0.0d0)
      do lam=1,ncf
        do imo=1,nval
          dmaxj=zero
          do jmo=1,nval
            cj=abs(C(lam,jmo))
            if(imo.ne.jmo) cj=cj*rfact(imo,jmo)
              dmaxj=max(dmaxj,cj)
          end do !j
          dmaxi(lam,imo)=dmaxj
        end do   !i
      end do     !lam
c
      do my=1,ncf
        do lam=1,my
          dmxx=0.0d0
          do imo=1,nval
            cprod=abs(C(my,imo)*dmaxi(lam,imo))
            dmxx=max(dmxx,cprod)
          end do  !imo
          D(my,lam)=dmxx
          D(lam,my)=dmxx
        end do    ! lam
      end do      ! my
CPP
c
C      call evaldmax(ncf,D,iprnt)
C      call prntmat(ncf,ncf,ncf,D)
C
C  construct a matrix over contracted shells rather than basis function
C
      ictr=igetival('ictr')
      ipoint=0
      do ics=1,ncs
        call get_shell_size(bl(ictr),ics,isize)
        kpoint=0
        do kcs=1,ncs
          call get_shell_size(bl(ictr),kcs,ksize)
          dmxs=0.0d0
          do my=1,isize
            do lam=1,ksize
            dmx=D(ipoint+my,kpoint+lam)
            dmxs=max(dmxs,dmx)
          end do     !lam in kcs
        end do     ! my in ics
        kpoint=kpoint+ksize
ckw
        dmxs=min(dmxs,1.d0)
ckw
ckw     DS(ics,kcs)=dmxs
        DS(ics,kcs)=dmxs*dmxs
c
      end do     ! kcs
      ipoint=ipoint+isize
      end do     !ics
      end
c======================================================================
ckw density for pre-screening based on Canonical MOS
c
      subroutine DmxMakeC(dualbasis, nval, nmo, ncf, ncs,
     *                    np4, inx, iprnt)
c-----------------------------------------------------
C      Main routine for construction Dmax(ncs,ncs) matrix
C      to be used for integral pre-screening for MP2
C      with Canonical orbitals.
C      Dmax(m,l) = max(i,j){abs(L(m,i)*L(l,j)}
C      Input:
c   dualbasis - logical : one basis or two are used
C      nval number of valence orbital
C      nmo number of occupied orbitals
C      ncf number of contracted basis functions
C      ncs number of contracted shells
C      np4 unit number , normally 4
C      inx array with shell informations
C      Output:
C      the final dmax is in 'dsmx'
c-----------------------------------------------------
C

      use memory

      implicit real*8(a-h,o-z)
      logical dualbasis
ckw
      character*256 jobname,MOS
      common /job/jobname,lenJ
ckw
      dimension inx(*)
c     common /big/bl(30000)
c----------------------------------------------------
      ncore=nmo-nval
c
cc    call matdef('dsmx','q',ncs,ncs)      ! already defined
      call matdef('dicf','q',ncf,ncf)
c
      idics=mataddr('dsmx')
      idicf=mataddr('dicf')
c-----------------------------------------------------
c for both cases the ordinary and dual-basis
c orbitals are already in 'cano'
c
c check the maximum values of BOTH the occupied and
c virtual eigenvectors
c
      call matsub('occo','cano',ncore+1,nmo)
      ioccu=mataddr('occo')
c
cccc  write(6,*)' Occupied orbitals used for integral screening'
      norbitals=nval
c
      call getmem(ncf*norbitals,itemp)
      call dmxijc_new(bl(idicf),  bl(idics),  ncf,  ncs,  norbitals,
     1                inx,        bl(ioccu),  bl(itemp),  iprnt)
      call retmem(1)
c-----------------------------------------------------------------
      if(dualbasis) then
         call getival('mpres_in_bg',mpres_in_bg)
         call updateDS(ncs,bl(idics), bl(mpres_in_bg))
      endif
c-----------------------------------------------------------------
c     call matprint('dsmx',6)
c-----------------------------------------------------------------
c release memory & remove matrices :
c
      call matrem('occo')
      call matrem('dicf')
c
      return
      end
C=====================================================================
      subroutine updateDS(ncs,ds, mpres_in_bg)
      implicit real*8(a-h,o-z)
      dimension ds(ncs,ncs)
      dimension mpres_in_bg(ncs)
      parameter(one=1.0d0)
c
      do ics=1,ncs
         ip=mpres_in_bg(ics)
         do jcs=1,ncs
            jp=mpres_in_bg(jcs)
            if(ip.eq.0 .or. jp.eq.0) then
                ds(ics,jcs)=one
            endif
         enddo
      enddo
c
      end
C=====================================================================
      subroutine DmxIJc(D,DS,ncf,ncs,nval,inx,C, iprnt)
c------------------------------------------------------------------
C      Calculates Max  of C(my,i) * C(lam,j)
C      Used for prescreening of integrals.
C      Output:
C      D(ncf,ncf) (this matrix is removed after return)
C       DS(ncs,ncs) same as D but over shells rather than
C      individual functions
C      Input:
C      ncf - number of basis functions
C      ncs - number of contracted shells
C      nval - number of valence orbitals
C      C(ncf,nval) matrix of MO coefficients (localized orbitals)
C      Svein Saebo, Fayetteville AR June 1999
c      Use the accelerated version dmxij_new
c------------------------------------------------------------------

      use memory

      implicit real*8(a-h,o-z)
      dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs)
      dimension inx(*)
      do my=1,ncf
        do lam=1,my
        dmxx=0.0d0
          do imo=1,nval
            do jmo=1,nval
              cprod=abs(C(my,imo)*C(lam,jmo))
              dmxx=max(dmxx,cprod)
            end do
          end do
          D(my,lam)=dmxx
          D(lam,my)=dmxx
        end do
      end do
C     call evaldmax(ncf,D,iprnt)
C     call prntmat(ncf,ncf,ncf,D)
C
C  construct a matrix over contracted shells rather than basis function
C
      ictr=igetival('ictr')
      ipoint=0
      do ics=1,ncs
        call get_shell_size(bl(ictr),ics,isize)
        kpoint=0
        do kcs=1,ncs
          call get_shell_size(bl(ictr),kcs,ksize)
          dmxs=0.0d0
          do my=1,isize
            do lam=1,ksize
            dmx=D(ipoint+my,kpoint+lam)
            dmxs=max(dmxs,dmx)
          end do     !lam in kcs
        end do     ! my in ics
        kpoint=kpoint+ksize
ckw
        dmxs=min(dmxs,1.d0)
        DS(ics,kcs)=dmxs*dmxs
c
      end do     ! kcs
      ipoint=ipoint+isize
      end do     !ics
c     call prntmat(ncs,ncs,ncs,DS)
c     call evaldmax(ncs,DS,iprnt)
      end
c======================================================================
      subroutine DmxIJc_new(D,   DS,  ncf,  ncs,  nval,
     1                      inx, C,   dmax, iprnt)
c------------------------------------------------------------------
C      Calculates Max  of C(my,i) * C(lam,j)
C      Used for prescreening of integrals.
C      Output:
C      D(ncf,ncf) (this matrix is removed after return)
C       DS(ncs,ncs) same as D but over shells rather than
C      individual functions
C      Input:
C      ncf - number of basis functions
C      ncf - number of basis functions
C      ncs - number of contracted shells
C      nval - number of valence orbitals
C      C(ncf,nval) matrix of MO coefficients (localized orbitals)
c      dmax  =storage, ncf
c      This is a much faster version of dmxijc
c      P. Pulay, Spring 2002
c------------------------------------------------------------------

      use memory

      implicit real*8(a-h,o-z)
      dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs),dmax(ncf)
      dimension inx(*)
      parameter(zero=0.0d0)
c
      do lam=1,ncf
        dmaxj=zero
        do jmo=1,nval
          cj=abs(C(lam,jmo))
          dmaxj=max(dmaxj,cj)
        end do !j
        dmax(lam)=dmaxj
      end do     !lam
c
      do mu=1,ncf
        do lam=1,mu
          D(mu,lam)=dmax(mu)*dmax(lam)
          D(lam,mu)=D(mu,lam)
        end do ! lam
      end do   ! mu
c
C     call evaldmax(ncf,D,iprnt)
C     call prntmat(ncf,ncf,ncf,D)
C
C  construct a matrix over contracted shells rather than basis function
C
      ictr=igetival('ictr')
      ipoint=0
      do ics=1,ncs
        call get_shell_size(bl(ictr),ics,isize)
        kpoint=0
        do kcs=1,ncs
          call get_shell_size(bl(ictr),kcs,ksize)
          dmxs=0.0d0
          do my=1,isize
            do lam=1,ksize
            dmx=D(ipoint+my,kpoint+lam)
            dmxs=max(dmxs,dmx)
          end do     !lam in kcs
        end do     ! my in ics
        kpoint=kpoint+ksize
ckw
        dmxs=min(dmxs,1.d0)
        DS(ics,kcs)=dmxs*dmxs
c
      end do     ! kcs
      ipoint=ipoint+isize
      end do     !ics
c     call prntmat(ncs,ncs,ncs,DS)
c     call evaldmax(ncs,DS,iprnt)
      end
c======================================================================
      subroutine EvalDmax(ncs,D,iprnt)
c      This subroutine evaluates the magnitudes of the elements of
c      the "Density" matrix used for pre-screening of integrals
c      These matrices are generated in Subroutine DmaxIJ.
c      Called from DmaxIJ
c       the subroutine can of course be used to evaluate any quadratic
c      matrix D of order ncs.
c  Arguments:
c  INTENT: IN
c  ncs= dimension of the matrix
c  D = matrix, D(ncs,ncs)
c  iprnt, print flag, statistics    are printed if it is 3 or greater
c  there are no arguments with INTENT: OUT
C
C      Svein Saebo , Fayetteville, AR  june 1999
C
      implicit real*8(a-h,o-z)
      dimension D(ncs,ncs)
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      do jmo=1,ncs
        do imo=jmo,ncs
          elem=D(imo,jmo)
          if (elem.gt.1.0d-1) then    !  > 0.1
            n1=n1+1
            cycle
          endif
          if (elem.gt.1.0d-2) then    !   >0.01
            n2=n2+1
            cycle
          endif
          if (elem.gt.1.0d-3) then    !   >0.001
            n3=n3+1
            cycle
          endif
          if (elem.gt.1.0d-4) then    !   >0.0001
            n4=n4+1
            cycle
          endif
          n5=n5+1
        end do
      end do
      ntot=ncs*(ncs+1)/2
      if(iprnt.ge.2) then
        iout=igetival('iout')
        write(iout,*)n1,' number of D elements > 0.1'
        write(iout,*)n2,' number of D elements > 0.01 and < 0.1'
        write(iout,*)n3,' number of D elements > 0.001 and < 0.01'
        write(iout,*)n4,' number of D elements > 0.0001 and < 0.001'
        write(iout,*)n5,' number of D elements < 0.0001'
        write(iout,*) 'total number of elements: ', ntot
      end if
      end
C======================================================================
      subroutine rfactor(x,y,z,r2,nval,rfac,epsi)
C.....X,Y,Z,r2 are centers and diameter of localized orbitals
C..... taken from LOCA
      implicit real*8(a-h,o-z)
      dimension x(*),y(*),z(*),r2(*)
      dimension rfac(nval,nval),epsi(*)
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,xmult=1.0d4)
      do imo=1,nval
        do jmo=1,imo
          if(imo.eq.jmo) cycle
C....calculate distance between orbital centers
          x2=(x(imo)-x(jmo))**2
          y2=(y(imo)-y(jmo))**2
          z2=(z(imo)-z(jmo))**2
          dd=(x2+y2+z2)+r2(imo)+r2(jmo)
c          d6=dd**3
          d4=dd**2
          rrI=sqrt(r2(imo)-x(imo)**2-y(imo)**2-z(imo)**2)
          rrJ=sqrt(r2(jmo)-x(jmo)**2-y(jmo)**2-z(jmo)**2)
c          rirj=sqrt(rri*rrj)
c          fact=xmult*rrI*rrJ/(d6*abs(epsi(imo)+epsi(jmo)))
          fact=xmult*rrI*rrJ/(d4*abs(epsi(imo)+epsi(jmo)))
          rfac(imo,jmo)=fact
          rfac(jmo,imo)=fact
        enddo
      enddo
      end
C======================================================================
      subroutine CompactCoef(ncf,nval,coef,nrow,irow,xtrunc)
      implicit real*8(a-h,o-z)
      dimension coef(ncf,nval),xtrunc(nrow,nval),irow(nrow)
c  This routine prepares a submatrix of the SCF coefficient matrix
c  coef(ncf,nval) according to the index array irow.
c  the resulting matrix (nrow x nval) contains only the rows in
c  irow(k)
      do i=1,nval
         do mu1=1,nrow
            mu=irow(mu1)
            xtrunc(mu1,i)=coef(mu,i)
         end do
      end do
      end
c======================================================================
      subroutine PutTrunc(ncf,nrow,ncol,xmat,ymat)
      implicit real*8(a-h,o-z)
      dimension xmat(ncf,ncf),ymat(nrow,ncol)
      do j=1,ncol
        do i=1,nrow
          ymat(i,j)=xmat(i,j)
        end do
      end do
      end
c======================================================================
      subroutine make_mapf2s(mapf2s,ncs,ncf,inx)
      dimension mapf2s(ncf)
      dimension inx(12,*)
c
      do ics=1,ncs
         icfb=inx(11,ics)+1
         icfe=inx(10,ics)
         do icf=icfb,icfe
            mapf2s(icf)=ics
         enddo
      enddo
c
      end
c======================================================================
      subroutine get_small_basis(bl,ncs_bg,ncf_bg,inx_bg,ncs_sm,ncf_sm)

      use memory, block => bl

      implicit real*8(a-h,o-z)
c---------------------------------------------------------------------
      dimension bl(*)
      dimension inx_bg(12,ncs_bg)       ! big basis set
c---------------------------------------------------------------------
c small basis set data :
c
      call getival('ncs_sm',ncs_sm)
      call getival('ncf_sm',ncf_sm)
      call getival('nsh_sm',nsh_sm)
      call getival('nbf_sm',nbf_sm)
c
c big basis set :
c
      call getival('ibas',ibas)       ! big basis set
c---------------------------------------------------------------------
c allocate memory for mapping arrays (from small to big basis set)
c
      call getint(ncs_sm,map_s2b)
      call getint(ncf_sm,maf_s2b)
      call setival('map_s2b',map_s2b)
      call setival('maf_s2b',maf_s2b)
c
c get memory address for an array showing if a big bs shell
c was present (0 or 1) in a small bs
c
      call getival('mpres_in_bg',mpres_in_bg)
c---------------------------------------------------------------------
c allocate integer memory for inx_sm
c
      call getint(12*ncs_sm,inx_sm)
c
c allocate memory for datbas_small & inx_small (to read from a disk)
c
      call matdef('dat_sm','r',13,nsh_sm)
      call matdef('inx_sm','r',12,ncs_sm)
      ibas_sm=mataddr('dat_sm')
      inxs_sm =mataddr('inx_sm')
c
      np4=4
      call matread('dat_sm',np4,'smallbas')
      call matread('inx_sm',np4,'smallinx')
c
c copy back inx_sm from double prec. to integer array :
c
      call get_inx_sm(bl(inxs_sm),12*ncs_sm,bl(inx_sm) )
c
      call matrem('inx_sm')
c---------------------------------------------------------------------
c Now the small basis set data are in 'dat_sm' (datbas_sm) and in inx_sm
c---------------------------------------------------------------------
c the small basis set must be entirely included in the big (current) one
c check if this is the case ;
c
      call check_if_included(bl(ibas),inx_bg,ncs_bg,
     *                      bl(ibas_sm),bl(inx_sm),ncs_sm,
     *                      bl(map_s2b),bl(maf_s2b),
     *                      bl(mpres_in_bg) )
c
c If execution was not stopped in the call above then the whole
c small basis is included in the big one so we may proceed.
c---------------------------------------------------------------------
c release memory :
      call matrem('dat_sm')
      call retmem(1)            ! allocated for integer inx_sm
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine get_inx_sm(dnx_sm,idim,inx_sm )
      implicit real*8(a-h,o-z)
      dimension dnx_sm(idim)
      dimension inx_sm(idim)
c
      do i=1,idim
         inx_sm(i)= int( dnx_sm(i) )
      enddo
c
      end
c======================================================================
      subroutine check_if_included(datbas_bg,inx_bg,ncs_bg,
     *                             datbas_sm,inx_sm,ncs_sm,
     *                             map_s2b,maf_s2b,mpres_in_bg)
      implicit real*8(a-h,o-z)
      logical found
      dimension inx_bg(12,ncs_bg)
      dimension inx_sm(12,ncs_sm)
      dimension datbas_bg(13,*)
      dimension datbas_sm(13,*)
      dimension map_s2b(ncs_sm)
      dimension maf_s2b(*)               !     (ncf_sm)
      dimension mpres_in_bg(ncs_bg)
c
      epsi=abs(epsilon(0.0d0))   ! epsilon  value for comparing real*8
      ncf_sm=0
      do ics_sm=1,ncs_sm
         iat_sm=inx_sm(2,ics_sm)       !   center on
         ism=inx_sm(1,ics_sm)+1        !   starting contr
         exp1_sm=datbas_sm(1,ism)
         csi1_sm=datbas_sm(2,ism)
c contr.func.
         icf1_sm=inx_sm(11,ics_sm)+1
         icf2_sm=inx_sm(10,ics_sm)
c
         ncf_sm=ncf_sm + icf2_sm-icf1_sm + 1
c
         found=.false.
         do ics_bg=1,ncs_bg
            iat_bg=inx_bg(2,ics_bg)       !   center on
            ibg=inx_bg(1,ics_bg)+1        !   starting contr
            exp1_bg=datbas_bg(1,ibg)
            exp1_cmp=abs((exp1_sm-exp1_bg)/max(exp1_sm,exp1_bg))
            csi1_bg=datbas_bg(2,ibg)
            csi1_cmp=abs((csi1_sm-csi1_bg)/max(csi1_sm,csi1_bg))
c
            if(exp1_cmp.le.epsi .and. csi1_cmp.le.epsi .and.
     *         iat_sm.eq.iat_bg) then
               found=.true.
               map_s2b(ics_sm)=ics_bg
c contr.func.
               icf1_bg=inx_bg(11,ics_bg)+1
               icf2_bg=inx_bg(10,ics_bg)
               increm=icf1_bg-icf1_sm
c
               do icf_sm=icf1_sm,icf2_sm
                  icf_bg=icf_sm+increm
                  maf_s2b(icf_sm)=icf_bg
               enddo
c
               exit
            endif
         enddo
c
         if(.not.found) then
            call nerror(1,'check_if_included',
     1 'shell of small basis not found in ext. basis', ics_sm,0)
         endif
      enddo
c----------------------------------------------------------------------
c find out which shell from the big basis set are not present in a small
c i.e. find extension shells
c
      call izeroit(mpres_in_bg,ncs_bg)
c
      do ics_sm=1,ncs_sm
         ics_bg=map_s2b(ics_sm)
         mpres_in_bg(ics_bg)=1
      enddo
c
c     do ics_bg=1,ncs_bg
c     write(6,321) ics_bg,mpres_in_bg(ics_bg)
c321  format('Big bs shell=',i3,' present in Smallbs:',i2)
c     enddo
c
      end
c======================================================================
      subroutine update_basis(bl,ncs_sm)
      implicit real*8(a-h,o-z)
      dimension bl(*)
c
c Place shells from a small basis set into the list of shells
c for the SECOND and the FOURTH basis sets for use in int_lmp2 :
c
      call getival('map_s2b',map_s2b)
      call getival('list_2',list_2)
      call getival('list_4',list_4)
c
      call setival('ncs_2',ncs_sm)
      call setival('ncs_4',ncs_sm)
c
      call update_list(bl(list_2),ncs_sm,bl(map_s2b))
      call update_list(bl(list_4),ncs_sm,bl(map_s2b))
c
      end
c======================================================================
      subroutine update_list(list,ncs,map_s2b)
      dimension list(ncs),map_s2b(ncs)
c
      do i=1,ncs
         list(i)=map_s2b(i)
ccccc    write(6,*)' list : ',i,' list(i)=',list(i)
      enddo
c
      end
c======================================================================
      subroutine get_mix_eigen(bl,nblocks,inx,labels,thresh,
     *                         ncf,ncf_sm,ncs,ncs_sm,nmo)

      use memory, block => bl

      implicit real*8(a-h,o-z)
      dimension bl(*)
      dimension inx(12,*)
c     common /intbl/ifpp,inxx(100)
      logical rhf
      character scftype*11
c----------------------------------------------------------------
c Input :
c bl(*)     - general storage
c nblocks   - total number of blocks (of quartets) for the old
c             integral program used here
c label     - a pointer to the integral's labels array
c inx(12,*) - general basis set info
c ncf, ncs  - dimensions of a big basis set
c ncf_sm,ncs_sm - dimensions of a small basis sets
c----------------------------------------------------------------
c This routine is used only for "dual basis set" MP2 .
c It returns eigenvectors & eigenvalues needed for MP2.
c For occupied orbitals MOs & Mos energies are those from
c ths "small basis set" SCF. For virtuals they are obtained
c by special diagonalization of a modified Fock matrix. "Special"
c means that the mixing between occupied & virtuals is not allowed.
c The modified Fock matrix is defined in the big basis set but the
c Fock operator is the Fock operator of the small basis, i.e.
c the small basis set density is used to build the Fock operator
c
c There are the following steps :
c
c (1) construct "big basis set" density by projection of the small
c     one (& zeros where nedeed)
c     Such density is needed to build the special Fock matrix
c     and it will be used to pre-screen integrals. Integrals
c     needed now have two indecies in a small basis set and
c     two others in a big basis set.
c (2) construct 2-el. part of the special Fock matrix
c     it is done by calling "old" integral program (int_fock
c     renamed here to int_fmp2) running over the big basis set
c     but the corresponding density used in prescreening should
c     eleminate all unneccesary integrals involving shells/functions
c     which were not present in a small basis set.
c     Dkl *{ (ij|kl) - 0.5(il|kj) } where kl belong to small bs only.
c (3) calculate 1-el. part of the full Fock matrix . Call inton
c     to get H0 (over big basis set) and add it to the previously
c     obatained 2-el. part.
c (4) diagonalize the Fock matrix in such a way that occupied and
c     virtuals do not mix. From that we take virtuals only i.e.
c     MOs & energise. The occupied ones are those obtained in ths
c     small basis set SCF.
c (5) calculate "basis set" correction to the energy :
c     E_bsc= - SUM(occ,vir) (Fov)**2/(Ev -Eo)
c
c----------------------------------------------------------------
      call mmark
c----------------------------------------------------------------
         call secund(tgmix0)
         call elapsec(egmix0)
c----------------------------------------------------------------
      np4=4
c
      iout=igetival('iout')
c----------------------------------------------------------------
c Make "big basis set" density :
c transfer density dens(ij) from the smalll basis set to the big basis.
c Put zeros where needed.
c......................................................
c allocate memory for "big"  and small density matrices
c
      ntri=ncf*(ncf+1)/2             ! big basis set
      ntri_sm=ncf_sm*(ncf_sm+1)/2    ! small basis set
c -- get the number of non-redundant basis functions
      lastvirt=igetival('nonredun')
c
      call matdef('densb','s',ncf,ncf)
      idensbig=mataddr('densb')
c
      call matdef('densm','s',ncf_sm,ncf_sm)
      idensmal=mataddr('densm')
c......................................................
c get small basis set density matrix
c
      call rea (bl(idensmal),ntri_sm,np4,'den0_rhf')
c......................................................
c make projection of dens. from small to big basis set
c
      call zeroit(bl(idensbig),ntri)
c
      call getival('maf_s2b',maf_s2b)
      call dens_s2b(bl(idensmal),ncf_sm,bl(idensbig),ncf,bl(maf_s2b))
ctest
c     write(6,*)' Small basis set density matrix'
c     call matprint('densm',6)
c     write(6,*)' Big   basis set density matrix'
c     call matprint('densb',6)
c......................................................
c
      call matrem('densm')
c----------------------------------------------------------------
c DO 1-EL PART which involves overlap S
c
      call getival('ibas',ibas)       ! big basis set
      call getival('inuc',inuc)       ! big basis set
      call getival('na  ',na)
c
c............................................................
c calculate S matrix in a big basis set
c
      call matdef('smat','s',ncf,ncf)
      ismat=mataddr('smat')
      call inton(0,na,bl(ismat),inx,0,0,bl(ibas),bl(inuc),ncs)
c............................................................
c Make  projector = I-(1/2)DS
c where I = unit mtx, D = small basis density
c (projected to the large basis), S = (large basis) overlap
c The matrix "u" is used as scratch for storing the projector
c
c make projector (in u) :
c
c
      call matdef('u','q',ncf,ncf)
c
      call make_proj('densb','smat','u')
c
c write projector on disk
c
      call matwrite('u',np4,0,'projector')
c
c............................................................
c calculate the LU factorization of S
c
      call cholesky('smat','u',ncf,xlow)
c
c write U matrix on disk
c
      call matwrite('u',np4,0,'umatrix')
c
c............................................................
      call matrem('u')
      call matrem('smat')
c----------------------------------------------------------------
c Check if the integral threshold should be sharpened
c check lowest eigenvalue of overlap in big basis set
c
      if(thresh.gt.xlow**2) then
        ithresh=-log10(thresh)
        if(ithresh.eq.10) then     ! default
           thresh = xlow**2
           thresh = max(thresh, 1.0d-12)
           write(iout,2000) thresh
        else
           write(iout,2100)
        endif
      endif
      call f_lush(iout)
 2000 format(/'  Near Linear Dependency in Big Basis Set',/
     *  '  Sharpening integral threshold to ',e14.6)
 2100 format(/'  Near Linear Dependency in Big Basis Set',/
     *  '  integral threshold would be sharpened',/
     *  '  (if not overwritten from the input)')
c----------------------------------------------------------------
c allocate memory for fock
c
      call matdef('fockmp2','s',ncf,ncf)
      ifockmp2=mataddr('fockmp2')
c
      call zeroit(bl(ifockmp2),ntri)
c----------------------------------------------------------------
         call secund(tgmix1)
         call elapsec(egmix1)
c----------------------------------------------------------------
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
c initialize the two-el. integral program
c
c    change for scf integrals iforwhat=1
c
      iforwhat=1
      nfock=1
      thint=thresh
      call getival('nslv',nslv)       ! no. of slaves
      If(nslv.GT.0) call para_JobInit(iforwhat)
c
      call ptwoint(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             0,      scftype,xintxx, nblocks,maxbuffer,
     *             maxlabels)

c----------------------------------------------------------------
c Calculate 1-el. part of the fock matrix : H0
c Construct the full H0 matrix in big basis set
c
      call getival('ibas',ibas)       ! big basis set
      call getival('inuc',inuc)       ! big basis set
      call getival('na  ',na)
c
      call matdef('h0mat','s',ncf,ncf)
      ihmat=mataddr('h0mat')
      call zeroit(bl(ihmat),ntri)
      call para_oneint(1,na,bl(ihmat),inx,0,0,bl(ibas), bl(inuc),ncs)
ctest
cprt  write(6,*)' H0 matrix '
cprt  call matprint('h0mat',6)
c
c add pseudopotential contribution if necessary
c
      call getival('npsp',npsp)
      if(npsp.ne.0)then
        call matdef('h0psp','s',ncf,ncf)
        ihpsp=mataddr('h0psp')
        call zeroit(bl(ihpsp),ntri)
        call getrval('ithr',psptol)
        psptol=psptol*0.01d0       ! psp thres.=0.01*main int. thres
        call psph0(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $             psptol,bl(ihpsp))
        call matadd('h0psp','h0mat')
        call matrem('h0psp')
      endif
c
c
c   COSMO section
c
      call tstival('cosmo',icosmo)    !if COSMO flag is not defined,
      If(icosmo.EQ.0) Then            ! define it now
        call setival('cosmo',0)
      Else
        call getival('cosmo',icosmo)
      EndIf
      if(icosmo.ne.0)then
        call getival('ndum',ndum)       ! no. of dummy atoms
        natom = natoms-ndum             ! no. of real atoms
        call getival('c_nps',nps)       ! no. of surface segments
        call getival('c_npsphe',npspher)! no. of outer surface segments
        call getival('nslv',nslv)       ! no. of slaves
        call getival('ibas',ibas)       ! big basis set
        call mmark
        call getmem(3*natom,icoxyz)     ! atomic ccordinates
        call getmem(3*(nps+npspher),icosurf) ! surface coordinates
        call getmem(nps,iar)            ! surface area
        call getmem(nps,iqcos)          ! surface charges
        call getint(nps,iiatsp)         ! mapping surface --> atom
        call getint(natom+1,icharge)    ! nuclear charges
c
c  read arrays from COSMO data file
c
        call cosmo_rdata(bl(icoxyz),bl(icosurf),bl(iar),bl(iqcos),
     $                   bl(iiatsp),bl(icharge),natom,nps)
        do i=0,3*npspher-1
           bl(icosurf+3*nps+i)=bl(icosurf+i)
        enddo
c
        call getmem(ntri,ivimat) ! electron-surface repulsion
        call getmem(ntri,ih0cos) ! cosmo contribution to H0
        call getmem(max(nps,npspher),iphi) ! surface potential
        call getmem(max(nps,npspher),iphin)! nuclear part of surf. pot.
c
c   compute nuclear part of surface potential
c
        call cosmo_potn(bl(icoxyz),bl(icharge),bl(iphin),
     $                  bl(icosurf),natom,nps)
        if(nslv.eq.0)then
          call cosmo_surfrep(bl(ivimat),bl(icosurf),inx,bl(ibas),
     $                       nps,ncs,ncf,ntri)
        else
          call para_cosmo_data(bl(icosurf),bl(icoxyz),
     $                         natom,nps,npspher)
          call para_cosmo_surfrep(nps)
        endif
      endif
c----------------------------------------------------------------
c setup parameters for int_fock (here int_fmp2) call :
c
      call getmem(maxlabels,labels)
c
      nfock=1
      thres1=thresh
      mywork=0
      igran=0
c
      idft=0
      ax=0.0d0
      rhf=.true.
c     fockB=
c     denB=
      mgo=1
c
      call pfock(idft,ax,nblocks,nfock,rhf,ncf,bl,inx,thres1,
     *           mgo,bl(idensbig),bl(ifockmp2),fockB,
     *           bl(idensbig),DenB,bl(labels),iforwhat)
c
ctj      call int_fmp2(nblocks,nfock,ncf,bl,inx,thres1,
ctj     *              bl(idensbig),bl(ifockmp2),bl(idensbig),
ctj     *              bl(labels),mywork,igran)
c................................................................
c NOTE :
c fock matrix on return is re-scaled by integ.threshold
c for parallel impl. call pfock instead of int_fock (int_fmp2)
c and re-scaling will be done in pfock (para_fock)
ctest
c     write(6,*)' 2-el part of MP2-Fock matrix'
c     call matprint('fockmp2',6)
c................................................................
c
      if(icosmo.ne.0) then
c
c  COSMO potential and one electron contribution
c
        if(nslv.eq.0)then
          call cosmo_pot(bl(idensbig),inx,bl(ibas),bl(iphin),
     $                   bl(iphi),bl(icosurf),bl(ivimat),
     $                   nps,ncf,ncs,ntri)
        else
          call para_cosmo_pot(bl(idensbig),bl(iphin),bl(iphi),
     $                        nps,ntri)
        endif
        call getrval('c_fepsi',fepsi)
        if(nslv.eq.0)then
          call cosmo_h0(bl(ih0cos),bl(ivimat),inx,bl(ibas),
     $                bl(icosurf),bl(iqcos),fepsi,nps,ncf,ncs,ntri)
        else
          call para_cosmo_h0(bl(ih0cos),bl(iqcos),fepsi,nps,ntri)
        endif
        call addvec(ntri,bl(ifockmp2),bl(ih0cos),bl(ifockmp2))
      endif
c
c  this is to fool para_next into not requesting the timing data
c  at this stage
c
      call para_next(5)
c
c  if this is a parallel run with cosmo, the slaves want to
c  compute the OC correction now, so we have to fool them
c
      if(icosmo.ne.0)then
        if(nslv.ne.0)then
          call para_cosmo_surfrep(0)
          call para_cosmo_pot(bl(idensbig),bl(iphin),bl(iphi),
     $                        npspher,ntri)
        endif
        call cosmo_del('c_onel  ')
        call retmark
      endif
c
c   now it is really over, we can get the timing data from the slaves
c
      mgo=0
      call para_next(mgo)
c......................................................
         call secund(tgmix2)
         call elapsec(egmix2)
c----------------------------------------------------------------
c form the full Fock matrix
c
      call matadd('h0mat','fockmp2')
ctest
c     write(6,*)' Final Fock matrix'
c     call matprint('fockmp2',6)
c............................................................
c save 'fockmp2' (needed later for basis set correction):
c
      call matwrite('fockmp2',np4,0,'fockmp2 ')
c----------------------------------------------------------------
c The final step : diagonalize full fock without
c                  mixing between occupied & virtual
c .................................................................
c  allocate memory
c
c use cano & epsi defined already before
c
      call matdef('u','q',ncf,ncf)
      call matdef('uinv','q',ncf,ncf)
c............................................................
c project out the occupied orbital subspace from the basis set
c everything below is in the large basis set
c projector = I-(1/2)DS where I = unit mtx, D = small basis density
c (projected to the large basis), S = (large basis) overlap
c
c get projector (in u) from a disk :
c
      call matread('u',np4,'projector')
c
c make projection :
c
      call matsimtr('fockmp2','u','fockmp2')
c
c............................................................
c get  U matrix from a disk and find its inverse
c
      call matread('u',np4,'umatrix')
      call u_inverse('u','uinv',ncf)
c
c............................................................
c Make the final diagonalization in order to get all virtuals
c
      zero=0.d0
      call geneig('fockmp2','u','uinv','epsi','cano','densb',zero,0,'U')
c
c............................................................
c Define the virtual subspace of this. The occupied orbitals
c give zero eigenvalues. What happens if there are negative
c orbital energy virtuals? one will have to check the epsilons.
c
      nmo=igetival('nmo')
      do i=1,nmo
         call matelem('epsi',i,i,x)
         if(abs(x).gt.1.0d-8) then
           call nerror(1,'get_mix',
     1          'one of the first nmo  eigenvalues is non-zero',i,i)
         endif
      end do
c
      call matrem('uinv')
      call matrem('u')
c
ctest
c     write(6,*)' Big basis set MOs & orbital energies'
c     call matprint('cano',6)
c     call matprint('epsi',6)
c----------------------------------------------------------------
c Form the final eigenvectors & eigenvalues for dual basis set MP2
c
c get small basis set eigens :
c
      call matdef('canosm','q',ncf_sm,ncf_sm)
      call matdef('epsism','d',ncf_sm,ncf_sm)
      call matread('canosm',np4,'evec_rhf')
      call matread('epsism',np4,'eval_rhf')
c
cprt  write(6,*)' Small basis set MOs & orbital energies'
cprt  call matprint('canosm',6)
cprt  call matprint('epsism',6)
c
      icoef_sm=mataddr('canosm')
      iepsi_sm=mataddr('epsism')
c
c big basis set eigens :
c
      icoef_bg=mataddr('cano')
      iepsi_bg=mataddr('epsi')
c
c
      call make_mp2_eigens(nmo,bl(maf_s2b),
     *                         bl(icoef_sm),bl(iepsi_sm),ncf_sm,
     *                         bl(icoef_bg),bl(iepsi_bg),ncf   )
c
c on return 'cano' & 'epsi' contain occupied MOs from small &
c                                   virtuals from big b.s.
c
ctest
cprt  write(6,*)' Final MOs & orbital energies for dual basis set MP2'
cprt  call matprint('cano',6)
cprt  call matprint('epsi',6)
c............................................................
c release memory
c
      call matrem('epsism')
      call matrem('canosm')
c----------------------------------------------------------------
c Calculate "basis set correction" to the SCF energy :
c
c     E_bsc= - SUM(occ,vir) (Fov)**2/(Ev -Eo)
c
c first make <OCC|FOCK|VIRT> block of fockmp2 :
c
c Fov  = Cocc(T) * Fao * Cvir
c
c get back 'fockmp2' :
c
      call matread('fockmp2',np4,'fockmp2 ')
      nvirt=lastvirt-nmo
      call matdef('fockov','r',nmo,nvirt)
      call matsub('occupied','cano',1,nmo)
      call matsub('virtuals','cano',nmo+1,lastvirt)
      call matdef('tempor','r',nmo,ncf)
c
      call matmmul2('occupied','fockmp2','tempor','t','n','n')
      call matmmul2('tempor','virtuals','fockov','n','n','n')
c
      call matrem('tempor')
      call matrem('virtuals')
      call matrem('occupied')
c
ctest
c     write(6,*)' < occ |Fock|vir. block of F'
c     call matprint('fockov',6)
c............................................................
c calculate correction :
c
      ifockov=mataddr('fockov')
      call calc_bsc(nmo,nvirt,bl(ifockov),bl(iepsi_bg),ncf,esing_exc)
      call setrval('esing_exc',esing_exc)
      call matrem('fockov')
c----------------------------------------------------------------
c integral timing :
c
      tgmix=(tgmix2-tgmix1)/60.d0
      egmix=(egmix2-egmix1)/60.d0
c
c total timings :
c
      call secund(tgmix3)
      call elapsec(egmix3)
      tgmit=(tgmix3-tgmix0)/60.d0
      egmit=(egmix3-egmix0)/60.d0
c
      write(iout,*) ' '
      write(iout,101)
  101 format(' CPU & Elapsed timings in dual basis set Fock ')
      write(iout,*) ' '
      write(iout,102) tgmix, egmix
  102 format('  Integral calculations =',f8.2,' and ',f8.2,' minutes')
      write(iout,103) tgmit, egmit
  103 format('  Total for Fock build  =',f8.2,' and ',f8.2,' minutes')
      write(iout,*) ' '
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------
      return
      end
c======================================================================
      subroutine make_mp2_eigens(nmo,maf_s2b,
     *                               coef_sm,epsi_sm,ncf_sm,
     *                               coef_bg,epsi_bg,ncf_bg)
      implicit real*8(a-h,o-z)
      dimension maf_s2b(ncf_sm)
      dimension coef_sm(ncf_sm,ncf_sm),epsi_sm(ncf_sm)
      dimension coef_bg(ncf_bg,ncf_bg),epsi_bg(ncf_bg)
c
      call dcopy(nmo,epsi_sm,1,epsi_bg,1)
c
      do imo=1,nmo
         call zeroit(coef_bg(1,imo),ncf_bg)
         do icf_sm=1,ncf_sm
            icf_bg=maf_s2b(icf_sm)
            coef_bg(icf_bg,imo)=coef_sm(icf_sm,imo)
         enddo
      enddo
c
      end
c======================================================================
      subroutine dens_s2b(densmal,ncf_sm,densbig,ncf_bg,maf_s2b)
      implicit real*8(a-h,o-z)
      dimension densmal(*)         ! input
      dimension densbig(*)         ! output
      dimension maf_s2b(ncf_sm)    ! input
c
      do icf_sm=1,ncf_sm
         icf_bg=maf_s2b(icf_sm)
         iis=icf_sm*(icf_sm-1)/2
         iib=icf_bg*(icf_bg-1)/2
         do jcf_sm=1,icf_sm
            jcf_bg=maf_s2b(jcf_sm)
            ijs=iis+jcf_sm
            ijb=iib+jcf_bg
            densbig(ijb)=densmal(ijs)
         enddo
      enddo
c
      end
c======================================================================
      subroutine int_fmp2(nblocks,nfock,ncf,bl,inx,thres1,
     *                    dens,fock,dn,labels,mywork,igran)
c---------------------------------------------------------------------
c This routine calls two-el. int. block by block
c and constructs the closed-shell Fock matrix .
c
c INPUT:
c  nblocks   -  number of superblocks, not strictly needed
c  nfock     -  number of Fock matrices to contruct
c  bl        -  common bl (free memory)
c  inx       -  common inx (contraction info)
c  thres1    -  threshold for integral*density
c  dens      -  density matrix, used to determine if a given
c               shell quartet is to be calculated
c               for open shell systems, contains alpha+beta densities
c  fock      -  alpha/closed-shell Fock matrix/matrices
c  dn        -  alpha/closed-shell density matrix/matrices
c
c  labels     = memory area for labels
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NOT EXISTING super block
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf,moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,NCX,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      common /cpu/ intsize,iacc,icache,memreal
      dimension inx(12,*)
      dimension bl(*),dens(*),fock(nfock,*),dn(nfock,*)
      dimension labels(*)
c-----------------------------------------------------------------
      call secund(time0)
c----------------------------------------------------------------
      where='    '
c----------------------------------------------------------------
c      nc22=ncf**2/intsize+1
c      call getmem(nc22,iarray)
       call getint(ncf**2,iarray)
       call setiarray(ncf,bl(iarray))
c----------------------------------------------------------------
c transform ordinary density dens(ij) into denspar(ics,jcs)
c
      call getmem(ncs*ncs,idensp)
      call getmem(ncf,map_fs)
      call setup_densp2(inx,ncf,ncs,dens,bl(idensp),
     *                  bl(map_fs))
      call retmem(1)
c----------------------------------------------------------------
      call mmark
c----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
      screen='fock'
c----------------------------------------------------------------
c set integral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c total number of super-blocks  nsupblks=nbl2*(nbl2+1)/2
c
      nsupblks=nblocks
      if(mywork.eq.0) then
         istart=1
         istop=nsupblks
      else
         istart=mywork
         istop=mywork+igran-1
      end if
c----------------------------------------------------------------
c
      DO isupbl=istart,istop
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,bl(idensp),where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c       check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c       integrals arrived in bl(ibuffz); check if they are :
c       nintez is the size of a given quartet.It is set up
c       to zero if there is no integrals
c
        if(nintez.gt.0) then
           call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
           lsh=labels(lab2)*labels(lab2+1)*labels(lab2+2)*labels(lab2+3)
c          ...........................
c           Standard HF
c          ...........................
           If(lsh.gt.1296) Then
             call fock_bldr(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
           Else
            if(labels(lab2+3).gt.1) then
              if(nfock.eq.1) then
                 call fock_bldr1(ncf,bl(ibuffz),fock,dn,
     *                     bl(iarray),nblsiz,ngctoz,nintez,
     *                     labels(lab1),labels(lab2),labels(lab3))
              else
                call fock_bldr2(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              end if
            else
              if(nfock.gt.1) then
                call fock_bldr3(nfock,ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              else
                call fock_bldr4(ncf,bl(ibuffz),fock,dn,
     *                    bl(iarray),nblsiz,ngctoz,nintez,
     *                    labels(lab1),labels(lab2),labels(lab3))
              end if
            end if
           EndIf
        endif
c
c----------------------------------------------------------------
        if(moreint) go to 11
c
 1111 continue
c
      END DO
      call retmark
c----------------------------------------------------------------
c       call dens_fock_prt(dens,fock,'fock -mp2')
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
c  DO NOT rescale the fock matrix here - it is made in para_fock
c
      call resc_nfock(ncf,nfock,fock,thres1)
      call secund(time2)
      call term_info(thres1,time2-time1,time1-time0,where)
c----------------------------------------------------------------
c release memory allocated for idensp
      call retmem(1)
c----------------------------------------------------------------

      end
c=====================================================================
      subroutine make_proj(densb,smat,u)
      implicit real*8 (a-h,o-z)
c-------------------------------------------------------------
c  This routine calculates the projector I-(1/2)DS
c  Arguments:
c  Input:
c    densb (character) = name of the density matrix
c    smat  (character) = overlap matrix
c  Output:
c    u     (character) = I-0.5*DS
c-------------------------------------------------------------
      character*(*) densb,smat,u
      parameter(zero=0.0d0, half=0.5d0, one=1.0d0,two=2.0d0)
c
      call matzero(u)
      call matinfo(u,mshape,idim1,idim2,iaddr,length)
      do i=1,idim1
         call mateset(u,i,i,-two)
      enddo
      call matmmul2(densb,smat,u,'n','n','a')
c  U = U + DS
      call matscal(u,-half)
c
      end
c===================================================================
      subroutine calc_bsc(nmo,nvirt,fockov,epsi,ncf,esing_exc)
      implicit real*8 (a-h,o-z)
c
      dimension fockov(nmo,nvirt), epsi(ncf)
c
      sum=0.d0
      do iocc=1,nmo
         eocc=epsi(iocc)
         do ivirt=1,nvirt
            evir=epsi(nmo+ivirt)
            fov=fockov(iocc,ivirt)
            sum=sum + fov*fov/(evir-eocc)
         enddo
      enddo
c
      esing_exc=-2.d0*sum
c
c factor of 2.0 comes from doubly occupied MO
c
      end
c===================================================================
      subroutine check_size1(lmp2_size,ncf,nval,mval,ntimes)

      use memory

      parameter (limit_size=3465600) !size limit of the mp2 iteg.buffer
c (h2o)20/6-31g*  dd split
c
      call getival('lcore',lcore)
      call getmem(0,last0)
      call retmem(1)
c
c
      nmemory=lcore-last0
      nmemory=nmemory - 1000000               ! leave some mem for integrals
      nmemory=nmemory - (nval*mval)/2+1       ! allocated in transOneShell
      nmemory=nmemory - (nval*mval)/8+1       ! allocated in transOneShell
      nmemory=nmemory - MAX(nval,mval)*ncf    ! allocated in transOneShell
c
c just for tests :
c
c     nmemory=limit_size
c
      if(lmp2_size.gt.nmemory) then
         ntimes=lmp2_size/nmemory
         if(ntimes*nmemory.lt.lmp2_size) ntimes=ntimes+1
      else
         ntimes=1
      endif
c
      end
c=======================================================================
      subroutine check_sizes(inx,ncs,ncf,nval)

      use memory

ctest parameter (limit_size=110976 ) !size limit of the mp2 iteg.buffer
      dimension inx(12,*)
c
      call getival('lcore',lcore)
      call getmem(0,last0)
      call retmem(1)
c
      nmemory=lcore-last0
      nmemory=nmemory - 1000000  !  leave some mem for integrals :
      nmemory=nmemory - nval*nval/2+1 ! allocated in transOneShell
      nmemory=nmemory - nval*nval/8+1 ! allocated in transOneShell
      nmemory=nmemory - nval*ncf      ! allocated in transOneShell
c
c just for test
c
c     nmemory=limit_size
c
c size of (Ics,all|Kcs,all)= Kics,kcs(ncf,ncf) buffers
c
      ncf2=ncf*ncf
      memper1=nmemory/ncf2
c
      isplit=0
      nsplit=0
      do ics=ncs,1,-1
         call get_shell_size(inx,ics,ics_size)
         do kcs=ics,1,-1
            call get_shell_size(inx,kcs,kcs_size)
            ik_size=ics_size*kcs_size
            itimes=ik_size/memper1
            if(itimes*memper1.lt.ik_size) itimes=itimes+1
            if(itimes.gt.1) then
               isplit=isplit+itimes
               nsplit=nsplit+1
            else
               ish=ics
               ksh=kcs
               go to 1234
            endif
         enddo
      enddo
c
 1234 continue
c
      if(isplit.gt.0) then
c        split starts from
         if(ish.eq.ksh) then
            ics=ish+1
            kcs=1
         else
            ics=ish
            kcs=ksh+1
         endif
         call get_shell_size(inx,ics,ics_size)
         call get_shell_size(inx,kcs,kcs_size)
c
         call getival('iout',iout)
         write(iout,10) nmemory,nsplit,ics_size,kcs_size,isplit,nsplit
         call f_lush(iout)
  10  format(
     *'  Because available memory limited to ',i9,' Words',/,
     *'  calculations of the K(ics,kcs)(ncf,ncf) matrices ',/,
     *'  for ',i4,' shell-pairs ( sizes : ',i2,' ,',i2,'   and up )',/,
     *'  will be done in ', i5   ,' passes ( instead of ', i4,' )',/)
      endif
c
      end
c=======================================================================
      subroutine set_passes(inx,ics,kcs,ics_size,kcs_size, ntimes,
     *                      Ipass,Kpass,Itimes,Ktimes)
c
      dimension inx(12,*)
      dimension ipass(2,*),kpass(2,*)
c
      kcf1=inx(11,kcs)+1
      kcf2=inx(10,kcs)
      icf1=inx(11,ics)+1
      icf2=inx(10,ics)
c
      if(ntimes.eq.1) then
         ipass(1,1)=icf1
         ipass(2,1)=icf2
         kpass(1,1)=kcf1
         kpass(2,1)=kcf2
         itimes=1
         ktimes=1
         return
      endif
c
c we need to split ics_size*kcs_size into ntimes pices.
c Try to do it first over kcs only
c
      katonce=kcs_size/ntimes
      if(katonce.gt.0) then
c        redefine ntimes according to maximum katonce:
         ntimes=kcs_size/katonce
         if(katonce*ntimes.LT.kcs_size) ntimes=ntimes+1
         nbeg=kcf1-1
         do itime=1,ntimes-1
           kpass(1,itime)=nbeg+1
           kpass(2,itime)=katonce+nbeg
           nbeg=nbeg+katonce
         enddo
         kpass(1,ntimes)=nbeg+1
         kpass(2,ntimes)=kcf2
c
         ktimes=ntimes
         itimes=1
         ipass(1,1)=icf1
         ipass(2,1)=icf2
c
         return
      endif
c
c if it could not be done over kcs then try to do it over ics only
c
      iatonce=ics_size/ntimes
      if(iatonce.gt.0) then
c        redefine ntimes according to maximum iatonce:
         ntimes=ics_size/iatonce
         if(iatonce*ntimes.LT.ics_size) ntimes=ntimes+1
         nbeg=icf1-1
         do itime=1,ntimes-1
           ipass(1,itime)=nbeg+1
           ipass(2,itime)=iatonce+nbeg
           nbeg=nbeg+iatonce
         enddo
         ipass(1,ntimes)=nbeg+1
         ipass(2,ntimes)=icf2
c
         itimes=ntimes
         ktimes=1
         kpass(1,1)=kcf1
         kpass(2,1)=kcf2
c
         return
      endif
c
      call nerror(1,'set_passes',' spliting over kcs & ics failed',0,0)
         stop
c
      end
c=======================================================================
      subroutine optimal_isize(isize,iwhole_size)
      implicit none
      integer isize,iwhole_size
c
      integer npass,itrysize,ipass
c
      npass=iwhole_size/isize
      if (mod(iwhole_size,isize).ne.0) npass=npass+1
      do itrysize=isize,1,-1
        ipass=iwhole_size/itrysize
        if (mod(iwhole_size,itrysize).ne.0) ipass=ipass+1
        if (ipass.ne.npass) then
          isize=itrysize+1
          exit
        endif
      enddo
      end
c=======================================================================
      subroutine integerst(m,xx,i4,i1)
c  This routine stores xx in a fixed-point representation in i4 and i1
c  Arguments
c  INTENT(IN)
c  xx  = 8-byte floating point value
c  m(7)= 1-byte masks in bits 5-7 of m. The values, if the last 5 bits
c        are omitted, are 1,2,3,4,5,6,7
c  INTENT(out)
c  i4  = 4-byte part of the stored integral
c  i1  = 1-byte part of the sored integral
c
c  This routine encodes a real number in a 4-byte and a 1-byte
c  integer. The 4-byte number (int4) contains the sign and leading bits
c  of the real*8 number. The leftmost 3 bits of int1 contain the exponent
c  in octal form, and the last 5 bits of the real number xx.
c  The number is: 8**ii*(dble(i4) + isign(i2/32,i4) where
c  ii= the leftmost 3 bits of i1 (can be 0 to 7), and i2=the rightmost
c  5 bits of i1 (can be 0 to 31)
c
      real*8 xx,bbig,r64,dblmax,a
      integer*4 i4,ii,mult
      integer*1 i1,m(7)
      parameter(dblmax=2147483648.0d0,r8=0.125d0)
      parameter(r82=r8**2,r83=r8**3,r84=r82**2,r85=r84*r8,r86=r84*r82,
     1          r87=r84*r83)
c
      if(abs(xx).lt.dblmax) then
        i4=xx
        i1=abs((xx-i4)*32)
      else
      bbig=xx
c      write(*,*)  'xx=',xx
      mult=1
      do ii=1,7
        bbig=bbig*r8
        mult=mult*8
        if(abs(bbig).lt.dblmax) exit
      end do
      if(abs(bbig).ge.dblmax) then
        call nerror(6,'MP2 module',
     $      'Integral cannot be encoded - Check your basis set',0,0)
      endif
      i4=bbig
c      write(*,'(3f20.2,f8.1)') abs(xx),abs(dble(i4))*dble(mult),
c     1        abs(dble(i4)),dble(mult)
      i1=abs(bbig-i4)*32
      i1=IOR(i1,m(ii))
c      write(*,*) 'xx,i4,mult,i1=',xx,i4,mult,i1
      end if
      end
c=======================================================================
      subroutine reconstruct(int4,int1,xx)
c  This routine reconstructs a real number from a 4-byte and a 1-byte
c  integer. The 4-byte number (int4) contains the sign and leading bits
c  of the real*8 number. The leftmost 3 bits of int1 contain the exponent
c  in octal form, and the last 5 bits of the real number xx.
c  ARGUMENTS
c  INTENT(IN)
c  int4=4-byte integer
c  int1=1-byte integer
c  INTENT(OUT)
c  xx=resulting real
c
      real*8 xx,bbig,r32,dblmax,yy,zz,r8(7)
      integer*4 int4,ii,mult
      integer*1 int1,m2,int2
      parameter(dblmax=2147483648.0d0,r32=0.03125d0)
      data m2/31/
      data r8/8.0d0,64.0d0,512.0d0,4096.0d0,32768.0d0,262144.0d0,
     1       2097152.0d0/
      ii=0
      ii=ishft(int1,-5)
      int1=IAND(int1,m2)
      xx=dble(int4)+sign(dble(int1)*r32,dble(int4))
      if(ii.eq.0) then
        return
      else
        xx=xx*r8(ii)
      end if
      end
c======================================================================
      subroutine check_coef(ncs,ncf,ncore,nmo,nval,coef)
      implicit real*8 (a-h,o-z)
      dimension coef(ncf,ncf)
c
      ioaomax=0
      iomomax=0
      xomax=0.d0
      do imo=ncore+1,nmo
         call absmax(ncf,coef(1,imo),imax,xmax)
         if(xmax.gt.xomax) then
            xomax=xmax
            ioaomax=imax
            iomomax=imo
         endif
      enddo
c
      ivaomax=0
      ivmomax=0
      xvmax=0.d0
      do imo=nmo+1,ncf
         call absmax(ncf,coef(1,imo),imax,xmax)
         if(xmax.gt.xvmax) then
            xvmax=xmax
            ivaomax=imax
            ivmomax=imo
         endif
      enddo
c
      write(6,66) xomax,iomomax, ioaomax, xvmax,ivmomax,ivaomax
  66  format(/,' Maximum coef ',f12.6,' occupied MO ',i4,' AO: ',i4,/
     *        ,' Maximum coef ',f12.6,' virtual  MO ',i4,' AO: ',i4,/)
c
      end
c=======================================================================
      subroutine set_mpres(mpres_in_bg,ncs)
      dimension mpres_in_bg(ncs)
      do i=1,ncs
         mpres_in_bg(i)=1
      enddo
      end

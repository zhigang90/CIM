#include "maxi_cim.h"
      subroutine para_cimsub

      use memory
      use newpara
      use kinds
      implicit real*8 (a-h,o-z)

      character word*4
      character*256 jobname

c-----------------------------------------------------------------------
c 0=logical option, 1=a single integer, mp2scal=p1*empanti+p2*emppar
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c-----------------------------------------------------------------------
      data ioptyp  / 1,     0,     0,     0,     0,     0,     11   /
      data word    /'prin','rimp','mp2 ','ccd ','ccsd','trip','thre'/
      data IUnit   /1/ ! unit number for checkpoint I/O
      parameter (nopt=7)
      common /job/ jobname,lenJ
      COMMON /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      dimension word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      

 10   format(72('='))
      write(iout,10)  
      write(iout,"(20x,' The Module of CIM Cluster Calculation')")
      write(iout,*) 
      
      iprnt=0
      thresh=1.0D-12
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
      if (ifound(1).gt.0) iprnt=iopv(1,1)
      if (ifound(7).gt.0) thresh=10.0d0**(-ropv(1,7))
C      if (ifound(2).gt.0) call CIM_RIMP2
      if (ifound(3).gt.0) call para_CIM_MP2(iprnt)
      if (ifound(4).gt.0) call para_CIM_CC(9,thresh,iprnt)
      if (ifound(5).gt.0) then
         if (ifound(6).gt.0) then
            call para_CIM_CC(11,thresh,iprnt)
         else
            call para_CIM_CC(10,thresh,iprnt)
         end if 
      end if
      write(iout,*) "End of CIM cluster calculation!"
      end subroutine para_cimsub


C ***************************************************
C * Parallel version of CIM cluster MP2 calculation *
C * Modify the source code of para_rmp2 in PQS      *
C * NZG_7/1/2016 @UARK                              *
C ***************************************************
      subroutine para_CIM_MP2(iprnt)

      use memory
      use newpara
      
      implicit real*8 (a-h,o-z)
      integer stat
      logical LastSymPair,emp2only,scs,igran_calculate
      dimension xintxx(9)
      common /job/ jobname,lenJ
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $             filname3,filname4
c ............................................................
      integer binlist(nslv)
c ............................................................
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20
      logical nofr,restrt,dualbasis,smal,keep,Test
      dimension xnmo(2)
      parameter (maxpair=700*701/2,maxgran=maxpair/2)   !FRAGILE!
      integer igran2pair(2,maxgran),ipair2gran(maxpair)

      integer,allocatable::IAN(:),Z(:),ILST(:,:),INX2(:,:)
      real*8,allocatable::XC(:,:),QA(:),BASDAT(:,:)
      real*8,allocatable::SMO(:,:),SMONEW(:,:),eorb(:),trans(:,:)
      character*8,allocatable::AtSymb(:)
      real*8,allocatable::tmpint(:,:,:),tmpintvir(:,:)
      real*8,allocatable::tranint(:,:,:,:),tmpintpair(:,:,:)
      real*8,allocatable::TR1(:),TR2(:,:),TC1(:),TC2(:,:)
      real*8,allocatable::PCIMA(:,:),PCIMB(:,:),ECIMA(:,:),ECIMB(:,:)
      real*8,allocatable::eqcmo(:,:),elmo(:,:),eorblmo(:)
      real*8,allocatable::x_QCMO(:,:,:),tij1_QCMO(:,:,:)
      real*8,allocatable::xqcmo2(:,:),tqcmo2(:,:)
      real*8,allocatable::x_LMO(:,:,:,:),tij1_LMO(:,:,:,:)
      real*8,allocatable::energypair_CIM(:,:),energyorb(:)

C For local MP2 calculation
      character*4 motype

C For normalization
C      real*8,allocatable::Sij(:,:),SOVER(:,:),tmp(:)
c
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
c-----------------------------------------------------------

C Read the number of orbitals of the cluster
c **Warning** The number of occ MOs cannot be determined by nele/2
C nmo:   number of occupied orbitals
C norb:  number of occupied and virtual orbitals
C trans: transformation matrix that transform QCMO to LMO
      open(unit=iunit,file=jobname(1:lenJ)//'.cim',form='formatted',
     &     status='old')
      call rdcntrl(iunit,5,'MOTYP',3,idum,rdum,motype)
      if (motype=='LMO ') then
         close(unit=iunit,status='keep')
         call para_CIM_CC(1,thresh,iprnt)
         return
      endif
      call rdcntrl(iunit,4,'NOCC',1,nmo,rdum,cdum)
      call rdcntrl(iunit,3,'NMO',1,norb,rdum,cdum)
      call rdcntrl(iunit,4,'NCEN',1,ncen,rdum,cdum)
      allocate(trans(nmo,nmo))
      call rread8(iunit,'$TRMX-A',nmo*nmo,trans(1,1))
      close(unit=iunit,status='keep')

      natoms=igetival('na')
      allocate(XC(3,NAtoms),AtSymb(NAtoms),IAN(NAtoms),QA(NAtoms))
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,jnk) ! In Bohr
      CLOSE(UNIT=IUnit,STATUS='KEEP')

C  get atomic numbers from atomic symbols
      CALL GetAtNo(NAtoms,AtSymb,IAN)
      call GetAtChrg(NAtoms,QA)

      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
      call getival('nsh',nsh)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      allocate(Z(natoms),ILST(4,ncs),BASDAT(13,nsh))
      call rdbasis(IUnit,NAtoms,AtSymb,XC,Z,ILST,BASDAT)
      deallocate(Z,XC,AtSymb,QA)
      CLOSE(UNIT=IUnit,STATUS='KEEP')

      allocate(INX2(12,ncs))
      CALL SortBAS1(NAtoms,ncs,ILST,INX2)
      deallocate(ILST)
      CALL normaliz(ncs,INX2,BASDAT)
      deallocate(BASDAT)

      i1=1
      i2=i1+ncf
      i3=i2+ncs
      IEnd=i3+12*ncs-1
c
      allocate(Z(IEnd))
      Call reorder(ncs,INX2,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinskiorder
      Call SortBAS2(ncs,Z(i3),Z(i2),Z(i1),INX2)        ! per atom

C -- get the overlap matrix for normalizing the orbitals
C      allocate(tmp(127008),Sij(ncf,ncf),SOVER(ncf,ncf))
C      call inton2(0,natoms,Sij,INX2,INX2,0,0,BASDAT,BASDAT,XC,IAN,
C     &            ncs,ncs,ncf,ncf,tmp)
C      call ReorderFock2(ncf,ncf,Z(1),Sij,SOVER)
      deallocate(IAN,INX2)

c -- set up the MOs filename
      call tstchval('mosfname',iyes)
      if (iyes.eq.1) then
         call getchval('mosfname',filename)
      else
         filename=jobname(1:lenJ)
      endif
      call rmblan2(filename,256,lenM)
      if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
      lenM=lenM+4

c -- read in the MOs
      itype=1
      allocate(SMO(ncf,norb),SMONEW(ncf,norb),eorb(norb))
      CALL ReadMOS(ncf,SMO,eorb,.True.,lenM,filename(1:lenM),itype,
     &             IErr)
      If(IErr.NE.0) Call nerror(3,'CIM Cluster RMP2 Module',
     $   'MOS File Does Not Exist',0,0)
       
C -- Reorder the coefficients to KW order
      DO J=1,norb
         DO I=1,ncf
            II=Z(I)
            SMONew(I,J)=SMO(II,J)
         end do
      end do
C      call normorb(ncf,norb,SMONew,SOVER)
      deallocate(SMO,Z)

C -- Write the reordered MOs to XXX.mos file
      filename(lenM+1:lenM+1)='2'
      lenM=lenM+1
      call WriteMOS(ncf,norb,SMONEW,eorb,.true.,lenM,filename,itype)
      deallocate(SMONEW,eorb)
      lastvirt=norb

      thresh=1.0D-10
      call setrval('thresh',thresh)
      restrt=.false.
      keep=.false.
      Test=.false.
c .............................................................................
c -- File handling --
c  filname1=half-transformed integrals; filname2=sorted bins
c  If gradient: filname3=Tij amplitudes;  filname4=Kov amplitudes
C  CIM: filname5=MO integrals
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.htr'
      len1 = len+4
      filname2=scrfile(1:len)//'.bins'
      len2 = len+5
c ............................................................................
      smal=.true.

C NZG Take care of the storage
      DiskMax = 200.0d0

c -- SCS scaling  ! SS
      scs=.false.
      iscs=0
      p1=1.2d0
      p2=1.0d0/3.0d0
      call setival('iscs',iscs)
      call setrval('p1',p1)
      call setrval('p2',p2)
 95   emp2only=.true.
      mp2only = 1
      call setival('mp2only',mp2only)
      igran_calculate=.true.
      itype=1
c-----------------------------------------------------------
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
      ncs_sm=ncs
      ncf_sm=ncf
      call set_mpres(bl(mpres_in_bg),ncs)
c-----------------------------------------------------------
      write(iout,*) ' MP2 integral thresh    = ',thresh
      write(iout,*) ' '
      call f_lush(iout)
c-----------------------------------------------------------
c  put down a memory marker
      call matreset
      call mmark
      call matmark
c
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
      CALL ReadMOS(ncf,bl(icano),bl(iepsi),.True.,lenM,
     $             filename(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(3,'CIM Cluster RMP2 Module',
     $   'MOS File Does Not Exist',0,0)
c
      lastvirt=norb
c-------------------------------------------------------
c initialize the two-el. integral program
c
      iforwhat=8
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
      ncore=0
      nfirst=1
      nlast=nmo
      nval=nmo
      nvirt = lastvirt-nmo
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
      dualbasis=.false.
      call DmxMakeC(dualbasis,nval,nmo,ncf,ncs,np4,inx,iprnt)
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
c -- send initial MP2 data to slaves
      call para_initsend
      call para_pack_int(nfirst,1)
      call para_pack_int(nlast,1)
      call para_pack_int(nval,1)
      call para_pack_int(norb,1)
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
      If(iprnt.gt.0) Then
         telap = (tend-tstart)/60.0d0
         write(iout,1000) telap
 1000    Format(' End of First Half-transformation    elapsed time = ',
     $           F9.2,' min')
         call f_lush(iout)
      EndIf
c
c .....................................................
      If(Test) Then
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
      If(lbin.gt.ncf*(ncf+1)) lbin = ncf*(ncf+1)
cPP -- limit lbin in case of large basis and large memory (?)
      maxlbin = 524287
      If(lbin.gt.maxlbin) lbin=maxlbin
c
      If(lbin.lt.ncf*(ncf+1).AND.lbin.lt.100) Then
        call nerror(5,'MP2 module',
     1       'memory available for bins leads to bin size < 100',
     2        mem,lbin)
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
      If(DiskLeft.LT.2.0d0) Then
        call nerror(3,'MP2 module',
     $   'Not Enough Disk Space left for bin sort',0,0)
      EndIf
c
      If(DiskNeed.GT.DiskLeft) Then
        NPass = NINT(DiskNeed/DiskLeft) + 1
      Else
        NPass = 1
      EndIf
c
      write(iout,'(a,a,i4,a)') ' Second Half-transformation will be',
     $                ' done in ',npass,' pass(es)'
      call f_lush(iout)
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
      lrec =(lbin     +8*lbin   +4  )*igranulesize +8  ! record length in bytes
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
      if(iprnt.gt.3) then
        do i=1,igranules
        write(iout,'(''Granule: '',i0,'' start: '',i0,'' stop: '',i0)')
     *                   i,igran2pair(1,i),igran2pair(2,i)
        enddo
        call f_lush(iout)
      endif
c -- determine number of bins and bin ownership list
c -- i.e., which pairs will be handled on which slaves
      ics = igranules/nslv          ! initial # of pairs per slave
      left = igranules - ics*nslv   ! # of pairs left
      nbin = ics                    ! maximum amount of gran. per slave
      if(left.gt.0) nbin = nbin+1
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
      write(iout,2200) lrec
      call f_lush(iout)
 2200 Format(' Record length of 2nd-transformed DA File is ',I10,
     $       ' bytes')
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
      If(NPass.GT.1) Then
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
c arrays for MO integrals storage
c Move it from outside pass cycle to here
c Zhigang 3/17/2016
c
      allocate(tmpint(nmo*(nmo+1)/2,nvirt,nvirt))
      allocate(tmpintvir(nvirt,nvirt))
      numint=nvirt*nvirt
      allocate(x_QCMO(nmo*(nmo+1)/2,nvirt,nvirt))
      allocate(tij1_QCMO(nmo*(nmo+1)/2,nvirt,nvirt))     

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
     &             call nerror(1,'MP2 master',
     $             'Internal binlist error ',binlist(nslv-1),igranules)
         endif
         If(IPass.EQ.NPass) binlist(nslv) = igranules
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

c -- receive the integrals from slaves for the last step of CIM
c -- Zhigang 3/17/2017 @UARK
         ncount=0
         do 
            call para_recv_pack(ifrom,TxCIMMP2Int)
            call para_unpack_int(i,1)
            if (i==0) then
               ncount=ncount+1
               if (ncount==nslv) then
                  exit           
               else
                  cycle
               endif
            endif
            call para_unpack_int(j,1)
            call para_unpack_real(tmpintvir(1:nvirt,1:nvirt),numint)
            call IntExpr2(ncf,nmo,nvirt,i,j,bl(iepsi),tmpintvir,
     &                    x_QCMO,tij1_QCMO)
         enddo

c -- now wait until get message from EACH original slave that
c -- this pass of the second half-transformation has completed
         Do i=1,nslv
            call para_recv(imsg,islave,TxBlockReq)
         EndDo
         if (iprnt.gt.2) 
     &      write(iout,'(a,i4,a)')' MASTER:  Finished pass ',ipass
         call f_lush(iout)
 75   CONTINUE
      call f_lush(iout)

c -- receive the MO integrals from all the slaves


C      do i=1,nslv
C         call para_recv(imsg,islave,TxCIMMP2Int)
C         write(ch3,'(".",I2.2)') islave
C         filname6 = filname5(1:len5)//ch3
C         len6 = len5+3
C         open(unit=50,file=filname6(1:len6),status='OLD',
C     &        form='UNFORMATTED',access='DIRECT',RECL=8)
C         read(50) tmpintvir(1:nmo,1:nmo,1:nvirt,1:nvirt)
C         close(50,status='DELETE')
C         tmpint=tmpint+tmpintvir
C     enddo

      deallocate(tmpintvir)
c
c -- check if all the slaves are OK
      call para_check
c
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

      call elapsec(tend)
      telap = (tend-tstart)/60.0d0
      tbin = tbin/60.0d0
      write(iout,1200) telap,tbin
1200  Format(' End of Second Half-Transformation   elapsed time = ',
     $       F9.2,' min',/,
     $       ' time spent in bin sort = ',F9.2,' min')
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
      allocate(x_LMO(nmo,nmo,nvirt,nvirt),tij1_LMO(nmo,nmo,nvirt,nvirt))
      x_LMO=0.0D0; tij1_LMO=0.0D0

      do k=1,nvirt
         do l=1,nvirt
            allocate(xqcmo2(nmo,nmo),tqcmo2(nmo,nmo))
            do j=1,nmo
               do ii=1,nmo
                  if (ii>=j) then
                     ij=ii*(ii-1)/2+j
                     xqcmo2(ii,j)=x_QCMO(ij,k,l)
                     tqcmo2(ii,j)=tij1_QCMO(ij,k,l)
                  else
                     ij=j*(j-1)/2+ii
                     xqcmo2(ii,j)=x_QCMO(ij,l,k)
                     tqcmo2(ii,j)=tij1_QCMO(ij,l,k)
                  endif
               enddo
            enddo
            call dgemm('T','N',nmo,nmo,nmo,1.0D0,trans,nmo,
     &                 xqcmo2,nmo,0.0D0,x_LMO(:,:,k,l),nmo)
            call dgemm('T','N',nmo,nmo,nmo,1.0D0,trans,nmo,
     &                 tqcmo2,nmo,0.0D0,tij1_LMO(:,:,k,l),nmo)
            deallocate(xqcmo2,tqcmo2)
         enddo
      enddo

C      do i=1,nmo
C         do j=1,nmo
C            do k=1,nvirt
C               do l=1,nvirt
C                  do ii=1,nmo
C                     x_LMO(i,j,k,l)=x_LMO(i,j,k,l)
C     &                              +x_QCMO(ii,j,k,l)*trans(ii,i)
C                     tij1_LMO(i,j,k,l)=tij1_LMO(i,j,k,l)
C     &                                 +tij1_QCMO(ii,j,k,l)*trans(ii,i)
C                  enddo
C               enddo
C            enddo
C         enddo
C      enddo

      allocate(energypair_CIM(nmo,nmo),energyorb(nmo))
      energypair_CIM=0.0D0; energyorb=0.0D0
      do i=1,nmo
         do j=1,nmo
            do k=1,nvirt
               do l=1,nvirt
                  x=x_LMO(i,j,k,l)
                  tij1=tij1_LMO(i,j,k,l)
                  energypair_CIM(i,j)=energypair_CIM(i,j)+tij1*x
               enddo
            enddo
         enddo
      enddo

      do i=1,nmo
         do j=1,nmo
            energyorb(i)=energyorb(i)+energypair_CIM(i,j)
         enddo
      enddo

      Ecluster=DSUM(nmo,energyorb,1)
      Ecen=PSUM(nmo,energyorb,NCEN)

      write(6,*) "Orbital energies for LMOs:" 
      write(6,'(5F20.9)') energyorb 
      write(6,'(A,F20.9)') ' E(CORR-CENTRAL)=',Ecen
      write(6,'(A,F20.9)') ' E(CORR-SUBSYS) =',Ecluster

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


C *********************************************************
C * Parallel mp2 slave subroutine for CIM-MP2 calculation *
C * Modify the subroutine do_mp2 in PQS                   *
C * NZG_7/1/2016 @UARK                                    *
C *********************************************************
      subroutine do_cimmp2

      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      character*256 filename,filname1,filname2,filname3,filname4,
     $              Restart_File,msgpass
      character ch3*3,scftype*11
      Logical rhf,smal,keep,restrt,keep1,Test,exst,alphabeta
      dimension xintxx(9)
      dimension ICSpass(2,28), KCSpass(2,28) !28 is 28 comp.of cart. I-function

      integer*4 SLAVE_ID(nslv),JBIN(nslv),IJSTRT(nslv)
      integer binlist(nslv)
      integer ipair(7),jpair(7)
      integer WaitAndLock, Unlock
      parameter (maxpair=700*701/2,maxgran=maxpair/2)   !FRAGILE!
      integer igran2pair(2,maxgran),ipair2gran(maxpair)
      integer*1 i01
      integer*4 i0
      real*8,allocatable::tmpint(:,:,:,:)
      
      Parameter (ndisk=50,ndisk2=60,ndisk5=70)
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
C      call getival('nonredun',lastvirt)
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
      call para_unpack_int(norb,1)
      lastvirt=norb
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
      IF(rhf) THEN
c -- closed shell RHF
        filname1 = filname1(1:len1)//ch3
        len1 = len1+3
        filname2 = filname2(1:len2)//ch3
        len2 = len2+3
        If(mp2only.EQ.-1) Then
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
 900   Format(A256,I4)
       CLOSE (UNIT=39,STATUS='KEEP')
      EndIf
c
      ILock = UnLock(Restart_File,lenR)
c
c .....................................................................
c -- open direct access file for half-transformed integrals
c
      If(rhf) Then
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
      IF(ics.gt.0) THEN
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
     *                      ntimes,Icspass,Kcspass,Itimes,Ktimes)
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
        If(smal) Then
          call getmem(lmp2_size,lmp2int)
          call mmark
          call int_lmp2(bl,     bl(ictr),thresh,
     *                  ics,icf1,icf2,kcs,kcf1,kcf2,
     1               bl(mapf2s),bl(idics),iprnt,bl(lmp2int),nintotal,
     2                  nrow,   ncol,    bl(irow),bl(icol),bl(lzero))
          call retmark
        Else
          ind=0
          intstore=0
          call getmem(lmp2_size,lrestor)
          call mmark
          call int_lmp2b(bl,   bl(ictr),thresh,     ics,     kcs,
     1               bl(mapf2s),bl(idics),iprnt, bl(lmp2int),nintotal,
     2               bl(icol),bl(irow), bl(if2cc),bl(if2cr),bl(indlj),
     3               bl(lrestor), nrow,  ncol,     ind,     intstore,
     4               lmp2_size, lenmax)
          call retmark
          indmax=max0(indmax,ind)
          istmax=max0(istmax,intstore)
          ttint=ttint+nintotal
          call retmem(1)               ! lrestor
        EndIf
c
c......................................................................
        IF(rhf) Then
        call TransOneShell(ncf,   ncs,      nval,      ics,   kcs,
     *                            icf1,     icf2,      kcf1,  kcf2,
     1                    ndisk,bl(ictr),bl(lmp2int),bl(ioccu),iprnt,
     2                  thresh,bl(ixadr),bl(ihalf), irec,  mulam,
     3                    mulamd, nrow,  ncol,     bl(irow), bl(icol),
     4                   bl(irow1),bl(icol1),ENonZero,smal)
c......................................................................
        if(smal) call retmem(1)        ! lmp2int
        ELSE
        call TransOneShlAB(ncf,    ncs,    nval,   nvalB,  ics,
     1                     kcs,    icf1,   icf2,   kcf1,   kcf2,
     2             bl(ictr),bl(lmp2int),bl(ioccu),bl(ioccuB),alphabeta,
     3                   iprnt,thresh,bl(ixadr),bl(ihalf), bl(ihalf),
     4             irec, bl(intsAA),bl(int1AA),bl(intsAB),bl(int1AB),
     5                   bl(intsBA),bl(int1BA),bl(intsBB),bl(int1BB),
     6                     mulam,  mulamd, nrow,   ncol, bl(irow),
     7                bl(icol),bl(irow1),bl(icol1),ENonZero,smal, ihalf,
     8                     mdisk1, mdisk2, mdisk3)
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
        If(rhf) Then
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
      If(Test) STOP
c .....................................................
c
c
c    B I N   S O R T  &  S E C O N D   T R A N S F O R M A T I O N
c    -------------------------------------------------------------
c
c
      call flush(6)
      IF(rhf) THEN
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
      
C MO integrals by ij pair index

      DO 95 IPass=1,NPass
c
c -- open the direct access file for writing
c
         OPEN (UNIT=ndisk2,FILE=filname2(1:len2),FORM='UNFORMATTED',
     $         ACCESS='DIRECT',RECL=lrec2,STATUS='NEW')
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
     &                                     lbin*igranulesize
         else
            nbbuf=0
         endif
         call getbyte(nbbuf,isendb)
c
c -- reopen the half-transformed integral file
c
         OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     $         ACCESS='DIRECT',RECL=lrec,STATUS='OLD')
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
     1      nval,      npairs,    ndisk,       ndisk2,    lbin,
     2      bl(ibin),  ijfirst,   ijlast,      binlist,   bl(i1),
     3      bl(i2),    bl(i3),    bl(i4),      bl(i5),    keep1,
     4      igran2pair,ipair2gran,igranulesize,igranules, ngran,
     5      imod,      bl(isendb))
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
            OPEN (UNIT=ndisk2,FILE=filname2(1:len2),FORM='UNFORMATTED',
     $            ACCESS='DIRECT',RECL=lrec2,STATUS='OLD')
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
               Call ReadBin(ncf,    lbin,   ndisk2, icount, ngran,
     $           bl(ibin),thresh, nsym,   nfirst, NAlpha,
     $           bl(iorbprA), bl(ifp),bl(i1), bl(i2), bl(ixm),
     $           igran2pair,    ipair2gran,II,   igranulesize,
     $           ncore)
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
                  call get_ij_half(i_ij,I,J)
                  I=I+ncore
                  J=J+ncore
                  Call PairEner(ncf,  NAlpha,   nvirt,  I,      J,
     $              bl(ixmo),bl(ieorb),   eij,    xnorm,  eijs,
     $              eijt,   eijpar, eijanti)

                  emp2 = emp2 + eij
                  tnorm = tnorm + xnorm
c -- Singlet and triplet components of the correlation energy
                  emps=emps+eijs
                  empt=empt+eijt
                  emppar=emppar+eijpar
                  empanti=empanti+eijanti
                  call para_initsend
                  call para_pack_int(I,1)
                  call para_pack_int(J,1)
                  call para_pack_real(bl(ixmo),nvirt*nvirt)
                  call para_send_pack(0,TxCIMMP2Int)
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
 70         CONTINUE

c tell master that all integrals have been sent
            call para_initsend
            call para_pack_int(0,1)
            call para_send_pack(0,TxCIMMP2Int)

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
 95   CONTINUE
C Send the MO integerals to master

C      OPEN (UNIT=ndisk5,FILE=filname5(1:len5),FORM='UNFORMATTED',
C     $      STATUS='REPLACE',access='DIRECT',recl=8)
C      write(ndisk5) tmpint(1:NAlpha,1:NAlpha,1:nvirt,1:nvirt)
C      call flush(ndisk5)
C      close(ndisk5,status='KEEP')
C      call para_send(MY_GID,0,TxCIMMP2Int)
C      numint=NAlpha*NAlpha*nvirt*nvirt
C      write(*,*) 'numint',numint

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
         Call BinSrt1AA(ncf,nval,nsym,mdisk4,nIntPerBin,nBinPerRec,
     &                  numrec,mrec,iloop,bl(iorbprA),bl(ifp),thresh,
     &                  jstrt,jend,bl(intsIJ),bl(int1IJ),bl(imu),
     &                  bl(ilam),bl(ibAA),bl(listAA),bl(listAA_back))
         call elapsec(t2)
         tsort = tsort + t2-t1
c
         Call TrnsVirtAA(ncf,nvirtA,nval,jstrt,jend,bl(ibAA),bl(ivrtA),
     &                   bl(itA),bl(iz),bl(ixtrn),eAA,tnorm,bl(listAA),
     &                   iprnt,ibAA,ivrtA,ixtrn,iz)
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


C ***************************************************************
C * Parallel version of CIM cluster Coupled-Cluster calculation *
C * Modify the parallel code of CC in PQS                       *
C * NZG_7/2/2016 @UARK                                          *
C ***************************************************************
      subroutine para_CIM_CC(imethod,thresh,iprint)

      use memory
      use newpara
      use kinds
      implicit real*8 (a-h,o-z)
      integer nbf
      logical LastSymPair,emp2only,equal,gauss_seidel,nodisk,fileopen
      logical variational,linear,cc,ccsd,norecalc,byt8,singles
      logical loca,omit,log_diis,mp2last,small,vorb,reset,reduced
      logical cache,calchunks,diis_less_thr,l_triples,qcisd,l_restart
      dimension xintxx(9)
      common /job/jobname,lenJ
      parameter(nopt=20)
      parameter(nmethod=10)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $              filname3,filname4
      character*256 chopv(nopt)
      character*6 methods(nmethod),method
      character opnames*4,scftype*11,wvfnc*20,ch3*3
      logical nofr,exst,restrt,dualbasis,smal,mp3,cep0,cep2,mp2,af
      character*2 int_kind
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension ifilestat(13)
      dimension xnmo(2)
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /extrone_tmp/ af,islaves,ndisktrc,ndisktrx,ndisktre,
     *                     ndisktrtx,ndisktrtc,ndisktrtt,mygid,nbf
      parameter (itim_no=20)
      dimension cctimesij(4),ccsdcpu(itim_no),ccsdela(itim_no)
      integer itids
      logical integerT
      common /GlobalCCSD/ cache,integerT,iprnt
      real*8 elaps(5)
      logical mp4,do_mp4,trans,signum
      integer ispectr(106) ! 0.2 from -5 to +15 (length: 21)
      integer :: info,ibuffer
      integer diis
      integer*4 iseed
      character *80, parameter :: spc='                             '
      character*80 dumpfile
      integer, parameter :: ndisk_dump = 59
      logical dump,last_dump_only,docansym
      integer,allocatable::IAN(:),Z(:),ILST(:,:),INX2(:,:)
      real*8,allocatable::XC(:,:),QA(:),BASDAT(:,:)
      real*8,allocatable::SMO(:,:),SMONEW(:,:),eorb(:),TX(:,:)
      real*8,allocatable::dtmp(:),Fock1(:,:)
      character*8,allocatable::AtSymb(:)
      real*8,allocatable::exc_QCMO(:,:,:),exc_LMO(:,:,:,:)
      real*8,allocatable::coef_QCMO(:,:,:),coef_LMO(:,:,:,:)
      real*8,allocatable::resi_QCMO(:,:,:),resi_LMO(:,:,:,:)
      real*8,allocatable::tmpexc(:,:),tmpcoef(:,:)
      real*8,allocatable::energypair_CIM(:,:),energyorb(:)
      real*8,allocatable::EPairQCMO(:,:),tmpresi(:,:)
      real*8,allocatable::CIMPairEnergy(:,:),CIMOrbEnergy(:)
      real*8,allocatable::eqcmo2(:,:),rqcmo2(:,:),cqcmo2(:,:)
c-----------------------------------------------------------------------
c
c  Explanation of the options:
c
c  CORR - Correlated method requested
c  TRIP - Include perturbative triples (for MP4, QCISD and CCSD)
c  NOFR - Correlate the core orbitals (no frozen core)
c  THRE - Integral threshold in pH format (default 10)
c  CORE - Orbitals with energy lower than this taken as core (default -3.0)
c  ORBS - Range of orbitals for correlation treatment (occupied AND virtual)
c  LVSH - Level shift for CI-based wavefunctions (default 0.0)
c  DIIS - For enhanced convergence, similar to SCF
c          first integer  =  size of DIIS subspace
c          second integer =  iteration to switch DIIS on
c  CONV - Convergence criteria on Energy and residuals, respectively
c         in pH format (default 6.0 in both cases)
c  MEMO - Specifies how much dynamic memory may be used in CORR
c         calculations. Without this option program would use as much
c         dynamic memory as possible which is dangerous - it may slow
c         down all the operating system, even shutdown may be impossible in
c         that case. Default guard-value is 10 MWords which may be too
c         little in most cases, but prevents infinite memory consumption.
c         The value below 2000 - Mword unit, above - word unit, where
c         word = 8 bytes (M=1.0d6, not 2**20).
c         ** Can also use units, MB or GB **
c  PRIN - Print level
c  ITER - maximum number of iterations
c  CHKP - Take checkpoint, i.e., save amplitudes at end of each iteration
c  DUMP - Save amplitudes at end of job
c  REST - Restart option. If specified alone, restart from last checkpoint;
c         if specified with TRIP assume amplitudes converged and compute
c         triples contribution
c
c -- the following are currently not for human consumption
c  CHUN - Manually determine size of matrix blocks
c  LOCA - Use localized orbitals
c  BYT5 - Store integrals as 5-byte quantities
c  AORB - Use AO-based instead of MO-based algorithm
c  CACH - store as much as possible "locally" on slaves disk - minimizes AF
c         usage and speeds up calculation, at the expense of heavy use of
c         disk storage.
c
      data IUnit/1/
c
      info = 0
      ibuffer = 0
      iseed = 0
      call srand(iseed)
      call reset_counters
      strace=0.0d0
      ik=0
      islaves=1
      af=.true.
c
      write(6,*)
      write(6,'(A)') "============================================="//
     *               "==========================="
      write(6,'(A)') "                    CIM Cluster Coupled"//
     *               " Cluster Module"
      write(6,'(A)') "============================================="//
     *               "==========================="
      write(6,*)
c
      call secund(tt0)
      call elapsec(elaps0)
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
      mygid=MY_GID
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c Defaults:
      nodisk=.false.
      gauss_seidel=.false.
      loca=.false.
      omit=.false.
      mp2=.false.
      mp3=.false.
      cep0=.false.
      cep2=.false.
      mp4=.false.
      do_mp4=.false.
      linear=.false.
      smal=.true.
      small=.true.
      cc=.false.
      norecalc=.true.
      byt8=.true.
      log_diis=.false.
      singles=.false.
      qcisd=.false.
      ccsd=.false.
      l_triples=.false.
      vorb=.true.
      cache=.false.
      variational=.true.
      integerT=.false.
      xmaxdisk=2500000000.0d0      ! currently not active
      diis_thres=0.01d0
      diis_less_thr=.false.
      itrstart=0
      nofr=.false.
      log_diis=.true.
      length_diis=6
      ifirst_diis_iter=1
c
c -- determine correlation method
      If(imethod.EQ.4) then
        write(6,'(A)') 'The CID energy will be calculated'
        method='CID'
        lmet=3
      else if(imethod.EQ.2) then
        write(6,'(A)') 'The MP3 energy will be calculated'
        mp3=.true.
        method='MP3'
        lmet=3
      else if(imethod.EQ.6) then
        write(6,'(A)')'The CEPA-0 energy will be calculated (doubles)'
        method='CEPA-0'
        lmet=6
        cep0=.true.
      else if(imethod.EQ.7) then
        write(6,'(A)') 'The CEPA-2 energy will be calculated (doubles)'
        method='CEPA-2'
        lmet=6
        cep2=.true.
      else if(imethod.EQ.5) then
        write(6,'(A)') 'The CISD energy will be calculated'
        method='CISD'
        lmet=4
        singles=.true.
      else if(imethod.EQ.1) then
        write(6,'(A)') 'The MP2 energy will be calculated'
        method='MP2'
        lmet=3
        mp2=.true.
C If one uses this subroutine to calculate MP2 energy, this means local
C MP2 calculation
        loca=.true.
      else if(imethod.EQ.9)  then
        cc=.true.
        write(6,'(A)') 'The CCD energy will be calculated'
        method='CCD'
        lmet=3
      else if(imethod.EQ.8)  then
        write(6,'(A)') 'The QCISD energy will be calculated'
        method='QCISD'
        lmet=5
        singles=.true.
        cc=.true.
        qcisd=.true.
      else if(imethod.EQ.10)  then
        write(6,'(A)')  'The CCSD energy will be calculated'
        method='CCSD'
        lmet=4
        ccsd=.true.
        singles=.true.
        cc=.true.
      else if(imethod.EQ.3) then
        write(6,'(A)')  'The MP4 energy will be calculated'
        method='MP4'
        lmet=3
        mp3=.true.
        mp4=.true.
      else if(imethod.EQ.11) then
        write(6,'(A)')  'The CCSD(T) energy will be calculated'
        l_triples=.true.
      else
        write(6,'(A)')  'No valid correlation method declared'
        STOP 'ERROR! No valid correlation method declared'
      EndIf

c -- set up name of dump file for restart
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      dumpfile=scrfile(1:len)//'.cc_dump'
      lenD = len+8
c
c  basis function symmetry pairs are in bl(ifp)
      natoms=igetival('na')
c
      call matreset
c
c .............................................................
c------------------------------------------------------------------
c -- THRE (integral threshold)
      call setrval('thresh',thresh)
c---------------------------------------------------------------------
c
c -- CORE is the limit separating core from valence
      core=-3.0d0
c
c -- LVSH
      shift=0.0d0
      shiftini=0.0d0      ! initial level shift for mp2 part
c
c -- CONV
      Ethresh=1.0d-6
      Wavethresh=1.0d-6
c
c -- MEMO
      max_dyn_mem=100000000000
 
      call dynamic_init(max_dyn_mem,bl(1))
c
c -- PRINt
      call setival('ccprint',iprint)
c
c -- ITER
      maxiter=40
c
c -- CHKP/DUMP
      dump=.false.
c
c -- RESTart
      l_restart=.false.
c
c -- CHUNk
      calchunks=.true.
c
c -- AORB
      vorb=.true.
      write(6,*) 'Calculation using MO basis'
c
c -- CACHe
      cache=.false.
      write(6,*) 'Calculation using af only'
c-----------------------------------------------------------
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call para_bcast(TxCIMSubCC,TxJobType)
      call fafinit(iresult)
c
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
      call getival('nsh',nsh)
c-----------------------------------------------------------
c -- check basis set (G is highest angular momentum basis function coded)
      call ChkAngMom(ncs,bl(ictr),2)
c-----------------------------------------------------------
      erhf=0.0D0 ! HF energy ! Just set it as zero
c----------------------------------------------
c
C Read the number of orbitals of the cluster
c **Warning** The number of occ MOs cannot be determined by nele/2
C nmo: number of occupied orbitals
C nbf: number of occupied and virtual orbitals
C      i.e. number of nonredundant basis functions (not true in CIM) 
C TX:  transformation matrix that transform QCMO to LMO
      open(unit=iunit,file=jobname(1:lenJ)//'.cim',form='formatted',
     &     status='old')
      call rdcntrl(iunit,4,'NOCC',1,nmo,rdum,cdum)
      call rdcntrl(iunit,3,'NMO',1,norb,rdum,cdum)
      call rdcntrl(iunit,4,'NCEN',1,ncen,rdum,cdum)

C For local calculation, we do not need the transformation matrix
      if (.not.loca) then
         allocate(TX(nmo,nmo))
         call rread8(iunit,'$TRMX-A',nmo*nmo,TX(1,1))
      endif
      LF2=ncf*(ncf+1)/2
      allocate(dtmp(LF2),Fock1(ncf,ncf))
      call rread8(iunit,'$AO-FOCK-A',LF2,dtmp(1))
      call FockFull(ncf,LF2,dtmp,Fock1)
      close(unit=iunit,status='keep')
c
      natoms=igetival('na')
      allocate(XC(3,NAtoms),AtSymb(NAtoms),IAN(NAtoms),QA(NAtoms))
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtoms,AtSymb,XC,-1,jnk) ! In Bohr
      CLOSE(UNIT=IUnit,STATUS='KEEP')

C  get atomic numbers from atomic symbols
      CALL GetAtNo(NAtoms,AtSymb,IAN)
      call GetAtChrg(NAtoms,QA)

      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      allocate(Z(natoms),ILST(4,ncs),BASDAT(13,nsh))
      call rdbasis(IUnit,NAtoms,AtSymb,XC,Z,ILST,BASDAT)
      deallocate(Z)
      CLOSE(UNIT=IUnit,STATUS='KEEP')

      allocate(INX2(12,ncs))
      CALL SortBAS1(NAtoms,ncs,ILST,INX2)
      CALL normaliz(ncs,INX2,BASDAT)

      i1=1
      i2=i1+ncf
      i3=i2+ncs
      IEnd=i3+12*ncs-1
c
      allocate(Z(IEnd))
      Call reorder(ncs,INX2,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinskiorder
      Call SortBAS2(ncs,Z(i3),Z(i2),Z(i1),INX2)        ! per atom

c -- set up the MOs filename
      call tstchval('mosfname',iyes)
      if (iyes.eq.1) then
         call getchval('mosfname',filename)
      else
         filename=jobname(1:lenJ)
      endif
      call rmblan2(filename,256,lenM)
      if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
      lenM=lenM+4

c -- read in the MOs
      itype=1
      allocate(SMO(ncf,norb),SMONEW(ncf,norb),eorb(norb))
      CALL ReadMOS(ncf,SMO,eorb,.True.,lenM,filename(1:lenM),itype,
     &             IErr)
      If(IErr.NE.0) Call nerror(3,'CIM Cluster CC Module',
     $   'MOS File Does Not Exist',0,0)
       
C -- Reorder the coefficients to KW order
      DO J=1,norb
         DO I=1,ncf
            II=Z(I)
            SMONew(I,J)=SMO(II,J)
         end do
      end do
C -- Write the reordered MOs to XXX.mos2 file
      filename(lenM+1:lenM+1)='2'
      lenM=lenM+1
      call WriteMOS(ncf,norb,SMONEW,eorb,.true.,lenM,filename,itype)


      call dynamic_matdef('cano','q',ncf,ncf)
      call dynamic_matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')

      CALL ReadMOS(ncf,bl(icano),bl(iepsi),.True.,lenM,
     $             filename(1:lenM),itype,IErr)
      If (IErr.NE.0) Call nerror(3,'CIM Cluster CC Module',
     $   'MOS File Does Not Exist',0,0)

c !!!!!Caution about nbf NZG
      nbf=norb
c
      call dynamic_lock(bl(icano),i)
      call dynamic_lock(bl(iepsi),i)
c
      write(6,'(A,I5)')'DIIS expansion length: ',length_diis
      write(6,'(A,I5)')'DIIS evaluation start: ',ifirst_diis_iter
      write(6,'(A,F6.1)') 'Level shift set to ',shift
      write(6,*) 'Threshold for energy  convergence:', Ethresh
      write(6,*) 'Threshold for wvfn.   convergence:', Wavethresh
C
C  Determine the orbitals to be correlated
C
      nfirst=1
      nlast=nmo
      nval=nlast
      ncore=0
      nend=norb

      call matsub('genvirt','cano',nmo+1,nbf)
      nvirt=nbf-nmo      ! number of virtual orbitals
c
      idimen=nvirt
c
      call twoelinit
c
      call open_mp3_resid(ndisk_mp3,idimen,.true.,mp4)
c  print the data of the calculation early
      write(iout,'("  Number of contracted functions =",i8,/,
     1             "  Number of core orbitals        =",i8,/,
     1             "  Number of correlated orbitals  =",i8,/,
     2             "  Number of virtual orbitals     =",i8)')
     3                ncf,ncore,nval,nbf-nmo
      write(iout,'(A,E20.5E2)') 'Integral threshold: ',thresh
      write(iout,*) ' '
      call flush(iout)
c
c Check orbital symmetry:
      call dynamic_unlock(bl(icano),i)
      call dynamic_unlock(bl(iepsi),i)
      call dynamic_matdef('overlap','q',ncf,ncf)
      ioverlap=mataddr('overlap')
      call OverlapBuilder(ioverlap)
      norb=nvirt+nval
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      icorr=icano+ncore*ncf
      call symm_memory_allocator(ncf)
      if (.not.loca.and..not.l_restart) then
         call dynamic_getmem(ncf*ncf,iscr1)
         call dynamic_getmem(ncf*ncf,iscr2)
         call dynamic_getmem(ncf*ncf,iscr3)
         call adv_symmetrize(nsym,   ncf,  ncore,   bl(ifp),
     *                       bl(icano),bl(ioverlap),
     *                       bl(iscr1),bl(iscr2),bl(iscr3),bl(iepsi))
         call adv_symmetrize(nsym,ncf,nval,bl(ifp),bl(icorr),
     &                       bl(ioverlap),bl(iscr1),bl(iscr2),bl(iscr3),
     &                       bl(iepsi+ncore))
         call adv_symmetrize(nsym,ncf,nvirt,bl(ifp),bl(icorr+nval*ncf),
     &                       bl(ioverlap),bl(iscr1),bl(iscr2),bl(iscr3),
     &                       bl(iepsi+nmo))
         call dynamic_retmem(3)
      endif
      call orb_pairs(nsym,ncf,nval,nvirt,bl(ifp),bl(ivalpair),
     *               bl(ivirpair),bl(icorr),bl(ioverlap))
cc
      docansym=.false.
      if (nsym.gt.0.and..not.loca) docansym=.true.
c
      call dynamic_getmem(ncf*ncf,ifockAO)
      call ReorderFock2(ncf,ncf,Z(1),Fock1,bl(ifockAO))
      deallocate(dtmp,Fock1,Z)
c
      ichar_count=1
      if (.not.loca) then
         call dynamic_getmem(8*7,ichatacters)
         call dynamic_getmem(8*8,im_table)
         call dynamic_getmem(8*(nval+1),iotable)
         call dynamic_getmem(nval,iorevtable)
         call dynamic_getmem(8*(nvirt+1),ivtable)
         call dynamic_getmem(nvirt,ivrevtable)
         call dynamic_getmem((nval*nval+1)*8,ijtable)
         call elapsec(et0)
         call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                      bl(ichatacters),ichar_count,bl(im_table))
         call elapsec(et1)
         write(6,'(A,I5)') 'Number of irreps: ',ichar_count
c Fock rediagonalization with symmetry:
c
         call dynamic_mmark()
         call dynamic_getmem(2*8,      ibassym)
         call dynamic_getmem(ncf*ncf,  imat1)
         call dynamic_getmem(ncf*ncf*8,imat2)
         call elapsec(et0)
         call symm_adapted(ncf,nsym,ichar_count,bl(ichatacters),bl(ifp),
     *                     bl(imat1),bl(imat2),bl(ibassym))
         call elapsec(et1)
         call dynamic_retmem(1)
         call elapsec(et0)

         if (ncf.eq.nbf) call Fock_Diag(ncf,ifockAO,ioverlap,icano,
     &                                  iepsi,bl(imat1),bl(ibassym),nmo)
         call elapsec(et1)
         call dynamic_retmark()
         call orb_pairs(nsym,ncf,nval,nvirt,bl(ifp),bl(ivalpair),
     *                  bl(ivirpair),bl(icorr),bl(ioverlap))
         call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                      bl(ichatacters),ichar_count,bl(im_table))
c
         call tables_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                         bl(ivalpair),bl(iotable),bl(iorevtable))
         call tables_generator(nsym,nvirt,ichar_count,bl(ichatacters),
     *                         bl(ivirpair),bl(ivtable),bl(ivrevtable))
         call ijtab_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                        bl(ivalpair),
     *                        bl(im_table),bl(iorevtable),bl(ijtable))
      endif
cc
      if (loca) then
         if(ncore.gt.0) call symmetrize_energies(ncf,nsym,ncore,
     &                                           bl(icorpair),bl(iepsi))
         call symmetrize_energies(ncf,nsym,nval, bl(ivalpair),
     &                            bl(iepsi+ncore))
         if(ncore.gt.0) call symmetrizer2(nsym,ncf,ncore,bl(ifp),
     &                                    bl(icorpair),bl(icano))
         call symmetrizer2(nsym,ncf,nval,bl(ifp),bl(ivalpair),bl(icorr))
         call symmetrizer2(nsym,ncf,nvirt,bl(ifp),bl(ivirpair),
     &                     bl(icorr+nval*ncf))
      endif
      call dynamic_lock(bl(iepsi),i)
      call pair_s_initializer(vorb,     nval,     nsym, ifp, ifp1,
     *                           ivalpair, ivirpair, loca,iorevtable,
     *                           ivrevtable,im_table,ichar_count,
     *                           ivtable)
c     call print_symmetries(nvirt,ncf)
      ij=0
      ij_unique=0
c Let us assume that all UNIQUE ij pairs are arranged in row in
c increasing order. Each of them will have a number assigned. Here I
c allocate memory for array which maps ij->ij_unique (from the row
c described above), or ij->0 if ij is a mirror of any pair
      npairs=(nval+1)*nval/2
      nstrong=npairs
      call dynamic_getmem(npairs,isympairs) ! for record numbering, used
c                                             to avoid empty records in files
      call dynamic_getmem(npairs+1,ipairimages) ! how many images given
c                                               pair has, if given pair
c                                               is already an image, = 0
      call izeroit(bl(ipairimages),npairs+1)
      call int_array_write(bl(ipairimages),npairs+1,1)
      do i=1,nval
         do j=1,i
            ij=ij+1
            call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                         ns, trans, signum)
            if (ns.eq.0) then
               ij_unique=ij_unique+1
               call int_array_write(bl(isympairs),ij,ij_unique)
            else
               call int_array_write(bl(isympairs),ij,0)
            endif
            ival=int_array(bl(ipairimages),ijprim)
            ival=ival+1
            call int_array_write(bl(ipairimages),ijprim,ival)
         enddo
      enddo
c 
c Write basic data to slaves:
      call para_initsend
      call para_send_info
c     call para_bcast_pack(TxScfInit)
c
c Distribute symmetry info:
      call para_initsend
      call para_pack_int(nbf,1)
      call para_pack_int(max_dyn_mem,1)
      call para_pack_int(npairs,1)
      call para_pack_int(bl(isympairs),npairs)
      call para_pack_int(bl(ipairimages),npairs+1)
      call para_pack_int(nval,1)
      call para_pack_int(nvirt,1)
      call para_pack_int(bl(ivalpair),7*nval)
      call para_pack_int(bl(ivirpair),7*nvirt)
      call para_pack_int(ij_unique,1)
      call para_pack_int(loca,1)
      call para_pack_int(docansym,1)
      if (.not.loca) then
         call para_pack_int(bl(iotable),8*(nval+1))
         call para_pack_int(bl(ivtable),8*(nvirt+1))
         call para_pack_int(bl(iorevtable),nval)
         call para_pack_int(bl(ivrevtable),nvirt)
         call para_pack_int(bl(im_table),8*8)
         call para_pack_int(bl(ijtable),(nval*nval+1)*8)
         call para_pack_int(ichar_count,1)
      endif
      call para_bcast_pack(TxCCSDInit)
c Write mp2 data:
      call para_initsend
      call para_pack_int(nodisk,1)
      call para_pack_int(gauss_seidel,1)
      call para_pack_int(omit,1)
      call para_pack_int(mp2,1)
      call para_pack_int(mp3,1)
      call para_pack_int(mp4,1)
      call para_pack_int(cep0,1)
      call para_pack_int(cep2,1)
      call para_pack_int(smal,1)
      call para_pack_int(small,1)
      call para_pack_int(cc,1)
      call para_pack_int(norecalc,1)
      call para_pack_int(byt8,1)
      call para_pack_int(log_diis,1)
      call para_pack_int(singles,1)
      call para_pack_int(ccsd,1)
      call para_pack_int(qcisd,1)
      call para_pack_int(vorb,1)
      call para_pack_int(nofr,1)
      call para_pack_int(cache,1)
      call para_pack_int(integerT,1)
      call para_pack_int(l_triples,1)
      call para_pack_int(l_restart,1)
c int:
      call para_pack_int(nmo,1)
      call para_pack_int(nvirt,1)
      call para_pack_int(nfirst,1)
      call para_pack_int(nlast,1)
      call para_pack_int(nval,1)
      call para_pack_int(idimen,1)
      call para_pack_int(npairs,1)
      call para_pack_int(length_diis,1)
      call para_pack_int(ndisk_mp3,1)
      call para_pack_int(ncore,1)
c real
      call para_pack_real(erhf,1)
      call para_pack_real(thresh,1)
      call para_pack_real(xmaxdisk,1)
      call para_pack_real(core,1)
      call para_pack_real(shift,1)
      call para_pack_real(shiftini,1)
      call para_pack_real(bl(icano),ncf*ncf)
      call para_pack_real(bl(iepsi),ncf)
c
      call para_bcast_pack(TxCCSDInit)
      call dynamic_lock(bl(icano),i)
c
      islaves=nslv
c reserve memory for list of strong, weak and distant pairs
      call getint(npairs,ilist)
c All pairs are filled with "1" which means that all are strong at this stage
c of calculation
      call filllist(bl(ilist),npairs)
ctj   First we need to have the coulomb integrals on the disk
ctj   Reserve space for coulomb and exchange integrals indexer
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrc)
      call izeroit(bl(irecadrc),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrx)
      call izeroit(bl(irecadrx),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadre)
      call izeroit(bl(irecadre),npairs*nslv)         !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtx)
      call izeroit(bl(irecadrtx),nval*nval*nslv)     !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtc)
      call izeroit(bl(irecadrtc),nval*nval*nslv)   !zero out the pair records
      call dynamic_getmem((nval*nvirt/intsize+1)*nslv,irecadrtt)
      call izeroit(bl(irecadrtt),nval*nvirt*nslv)   !zero out the pair records
ctj Generate the coulomb integrals
      int_kind='c'
      ipass=0
      call elapsec(ec_start)
      call secund(sc_start)
c     call fafInit(iresult)
      ndiskc=-1
      lrecc=0
      lbinc=0
      ndisktrc=-1
      if (.not.mp2.and..not.(l_restart.and.l_triples)) then
      call master_Gen(ncs,     ncf,           ictr,   nval,  nmo,
     *                nfirst,  nlast,         thresh, core,  xmaxdisk,
     *                ndiskc,  bl(irecadrc),  lrecc,  lbinc, int_kind,
     *                nodisk,  byt8,          small,  vorb,  ndisktrc,
     *                bl(ilist),bl(isympairs),ij_unique,iprnt,nbf)
      endif
c
      call secund(sc_stop)
      call elapsec(ec_stop)
      if (iprint.gt.2)
     *write(6,9) 'Elapsec time for coulomb integral calc. (sec): ',
     * ec_stop-ec_start
  9   format (A,F10.2)
      if (iprint.gt.2)
     *write(6,9) 'CPU time for coulomb integral calc. (sec):     ',
     * sc_stop-sc_start
      ipass=0
ctj Generate the exchange integrals
      int_kind='x'
      call elapsec(ex_start)
      call secund(sx_start)
      call master_Gen(ncs,     ncf,           ictr,   nval,  nmo,
     *                nfirst,  nlast,         thresh, core,  xmaxdisk,
     *                ndiskx,  bl(irecadrx),  lrecx,  lbinx, int_kind,
     *                nodisk,  byt8,          small,  vorb,  ndisktrx,
     *                bl(ilist),bl(isympairs),ij_unique,iprint,nbf)
      call secund(sx_stop)
      call elapsec(ex_stop)
      if (iprint.gt.2)
     *write(6,9) 'Elapsec time for exchange integral calc. (sec): ',
     * ex_stop-ex_start
      call flush(6)
      if (iprint.gt.2)
     *write(6,9) 'CPU time for exchange integral calc. (sec):     ',
     * sx_stop-sx_start
      call flush(6)
c 
      ipass=0
ctj Generate the three-external coulomb integrals
      if (l_triples) then
         int_kind='tt'
         call elapsec(ex_start)
         call secund(sx_start)
         call master_Gen(ncs,     ncf,          ictr,   nval,  nmo,
     *                   nfirst,  nlast,        thresh, core,  xmaxdisk,
     *                   ndisktt, bl(irecadrtt),lrectt, lbintt,int_kind,
     *                   nodisk,  byt8,        small,  vorb,  ndisktrtt,
     *                   bl(ilist),bl(isympairs),ij_unique,iprint,nbf)
         call secund(sx_stop)
         call elapsec(ex_stop)
      endif
      call flush(6)
c 
      call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,   ccsd,
     *                     vorb,  nmo,ivalpair,    ivirpair, ifp,
     *                     nsym,  nval,nvirt,      ij_unique,isympairs,
     *                     l_restart,ndisk_dump,nbf)
      call dynamic_getmem(ncf*ncf,ifockMO)
      call dynamic_matdef('corefock','q',idimen,idimen)
      icorefock=mataddr('corefock')
      call matconn('fockAO','q',idimen,idimen,ifockAO)
c
      if (.not.mp2) then
         call elapsec(et0)
         call secund(t0)
         call Kijk_Vec_Mast(bl(irecadrx),npairs,ndiskx,lbinx,ncf,thresh,
     *                      nfirst,nlast,byt8,nmo,vorb,.true.,kijkndisk,
     *                      nbf)
         call elapsec(et1)
         call secund(t1)
      endif
c
      if (mp2) lbinc=0 ! FRAGILE!
      call elapsec(et0)
      call secund(t0)
      call FockBuilder(bl(irecadrc),bl(irecadrx),nval, npairs, thresh,
     *                 ndiskc,      ndiskx,      lbinc,lbinx,  byt8, 
     *                 nmo,        vorb,     ifockAO, ifockMO,icorefock,
     *                 nbf)
      call sync_barrier
      call elapsec(et1)
      call secund(t1)
      if (iprint.gt.2) then
         write(6,9) 'FockBuilder ELAPS: ',et1-et0
         write(6,9) 'FockBuilder   CPU: ',t1-t0
         call flush(6)
      endif
      call Fock_Vector_Init(ncf,nfirst,nlast,ifockAO,ifockMO,nmo,vorb,
     *                      nbf)
c
      call para_bcast_real(bl(ifockAO),idimen*idimen,TxCCSDInit)
      call para_bcast_real(bl(ifockMO),nbf*nbf,TxCCSDInit)
      call para_bcast_real(bl(ioverlap),ncf*ncf,TxCCSDInit)
      call para_bcast_real(bl(icorefock),idimen*idimen,TxCCSDInit)
      call dynamic_lock(bl(icorefock),iamount)
      call dynamic_lock(bl(ifockAO),iamount)
      call dynamic_lock(bl(ifockMO),iamount)
c
      if (.not.mp2) then
         call CreateKijklDisk(nval,kijklndisk,af)
         call sync_barrier
         call elapsec(et0)
         call secund(t0)
c     call KijklInit(bl(irecadrx),npairs,ndiskx,lbinx,ncf,thresh,
c    *             nfirst,nlast,bl(ifockMO),byt8,nmo,vorb,kijklndisk,af)
         call KijklInit_Master(bl(irecadrx),npairs,ndiskx,lbinx,ncf,
     *                         thresh,nfirst,nlast,bl(ifockMO),byt8,
     *                         nmo,vorb,kijklndisk,af,nslv)
         call sync_barrier
         call elapsec(et1)
         call secund(t1)
         if (iprint.gt.2) then
            write(6,9)'KijklInit ELAPS: ',et1-et0
            write(6,9)'KijklInit   CPU: ',t1-t0
            call flush(6)
         endif
      endif
c
      call initsingles_master(nval,ncf,singles.or.mp4,nmo,vorb,
     *                              l_restart,ndisk_dump,nbf)
      call initS(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call initR(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
c
c  Memory for energy pairs:
      call dynamic_getmem(npairs,iepair)
      call dynamic_getmem(npairs,iepairr)
      call dynamic_getmem(npairs,ixepair)
      call dynamic_getmem(npairs,ixepairr)
c
      call mmark
      call dynamic_mmark
      call matmark
c Space (temporary) for one K_ij matrix
      call dynamic_matdef('resao','q',idimen,idimen)
      iresAO=mataddr('resao')
      call dynamic_matdef('exchao','q',idimen,idimen)
      iexchAO=mataddr('exchao')
c Different instances of HF coeff.
      call dynamic_matdef('invcano','q',ncf,ncf) !inversion of C + transposition
      iinvcano=mataddr('invcano')
      call matcopy('cano','invcano')
c
c     call matinv('invcano')
      call matmmul2('cano','overlap','invcano','t','n','n')

c
c     END of inverse refinement
c
      call dynamic_matdef('rinvcano','q',ncf,ncf) ! inversion of C
      call matcopy('invcano','rinvcano') !store r(eal) inverted cano for singles
      call dynamic_matdef('canotran','q',ncf,ncf)
      call matcopy('cano','canotran')
      call matpose('canotran')
c  This is necessary because of transf: T_MO=U-1 T_AO U-1(T)
c I use later matsimtr, which makes: U(T)A U
      call matpose('invcano')
c checking:
      call dynamic_matdef('unit','q',ncf,ncf)
      junit=mataddr('unit')
      call matmmult('rinvcano','cano','unit')
      do i=0,ncf-1
         bl(junit+i+i*ncf)=bl(junit+i+i*ncf)-1.0d0
      enddo
      summ=0.0d0
      xmax=0.0d0
      do i=0,ncf*ncf-1
         ele=dabs(bl(junit+i))
         if (ele.gt.xmax) xmax=ele
         summ=summ+ele
      enddo
      if (iprint.gt.3)
     * write(6,272) 'Final sum: ',summ,' max. element: ',xmax
 272  FORMAT (2(A,E15.5E2))
      call dynamic_matrem('unit')
c locking:
      i=mataddr('invcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('rinvcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('canotran')
      call dynamic_lock(bl(i),j)
c
      call dynamic_matdef('work1','q',idimen,idimen)
      iwork1=mataddr('work1')
      call dynamic_matdef('work2','q',idimen,idimen)
      iwork2=mataddr('work2')
      call dynamic_matdef('work3','q',idimen,idimen)
      iwork3=mataddr('work3')
      call dynamic_matdef('beta','q',nval,nval)
      ibeta=mataddr('beta')
c
      square=0.0d0
      energy=0.0d0
      energyr=0.0d0
      ij=0
      call matzero('resao')
      call elapsec(e1_start)
      call secund(s1_start)
c Determine if DIIS have to be performed or not
      if (.not.log_diis) then
         call CCDiis(i,       j,          'nodiis',     ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,   
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
      else
         call prepare_diis_files(ndiskdiisr,ndiskdiisc,af,ncf,
     *                           nmo,vorb,nbf)
         i=length_diis
         call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
      endif
      call sync_barrier
      call para_initsend
      call para_pack_int(ndiskdiisr,1)
      call para_pack_int(ndiskdiisc,1)
      call para_bcast_pack(TxCCSDInit)
c END DIIS
c     Singles definitions
c
      call dynamic_matdef('sing_res','r',idimen,1)
      ising_res=mataddr('sing_res')
      call dynamic_matdef('work_v1','r',idimen,1)
      iwork_v1=mataddr('work_v1')
      call dynamic_matdef('singl_am','r',idimen,1)
      isingl_am=mataddr('singl_am')
c
      square=0.0d0
      energy=0.0d0
      senergy=0.0d0
      energyr=0.0d0
      ij=0
      xmaxele=0.0d0
      if (l_restart) goto 126
c Skip if restart:
      do i=1,nval
         do j=1,i
            ij=ij+1
            call para_recv_pack(isgid,TxCCSDMP2Req)
            call para_unpack_int(isltid,1)
            call para_unpack_int(iresult,1)
            if (iresult.eq.1) then
c Result:
               call para_unpack_real(epair,1)
               call para_unpack_real(sepair,1)
               call para_unpack_real(sqij,1)
               call para_unpack_real(xmaxele_sl,1)
               if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
               energy=energy+epair
               senergy=senergy+sepair
               square=square+sqij
            endif
c Work:
            call para_initsend
            iwork=1
            call para_pack_int(iwork,1)
            call para_pack_int(i,1)
            call para_pack_int(j,1)
            call para_pack_int(ij,1)
            call para_send_pack(isgid,TxCCSDMP2Job)
         enddo
      enddo
c End the slaves and store the leftovers of data
      do i=1,nslv
         call para_recv_pack(isgid,TxCCSDMP2Req)
         call para_unpack_int(islid,1)
         call para_unpack_int(iresult,1)
         if (iresult.eq.1) then
c Result:
            call para_unpack_real(epair,1)
            call para_unpack_real(sepair,1)
            call para_unpack_real(sqij,1)
            call para_unpack_real(xmaxele_sl,1)
            if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
            energy=energy+epair
            senergy=senergy+sepair
            square=square+sqij
         endif
c Work:
         call para_initsend
         iwork=0
         call para_send(iwork,isgid,TxCCSDMP2Job)
      enddo
      call secund(s1_stop)
      call elapsec(e1_stop)
      if (iprint.gt.2) then
         write(6,9) 'Elapsec time for first iteration:  ',
     *   e1_stop-e1_start
         write(6,9) 'CPU time for first iteration:      ',
     *   s1_stop-s1_start
      endif
      if (.not.loca) then 
         xmp2=energy
         scsmp2=senergy
         if (iprint.gt.3) write(6,*)scsmp2
         if (mp2) goto 999
      else
         if (iprint.gt.2)
     *   write(6,'(A,F21.17)') 'Nonlocal! - "MP2" energy: ',energy
      endif
      call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,   ccsd,
     *                     vorb,  nmo,ivalpair,    ivirpair, ifp,
     *                     nsym,  nval,nvirt,      ij_unique,isympairs,
     *                     l_restart,ndisk_dump,nbf)
c
      mp2_iter=0
      mp2last=.false.
      write(6,'(A)') "     MP2 module:"
      write(6,'(A)')
      write(6,'(A)')
     * "Iter:  Energy:        Delta E           "//
     * "Max resid:    Err sq:      Elapsed time"
 222  continue
      call f_lush(6)
c
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      energyrold=energyr
      ij=0
      mp2_iter=mp2_iter+1
      energy=0.0d0
      senergy=0.0d0
      energyr=0.0d0
      energy1=0.0d0
      energy2=0.0d0
      sqij=0.0d0
      xmaxele=0.0d0
      if (mp2last) allocate(CIMPairEnergy(nval,nval),CIMOrbEnergy(nval))
      call secund(t0)
      call elapsec(et0)
      idimen1=idimen
      call Open_Scratch_File(ndiskmp2,idimen,idimen1,.true.,'lmp2')
      call sync_barrier
      call para_initsend
      call para_bcast(ndiskmp2,MP2file)
      call Generate_MP2_G_MAS(nval,ncf,ndiskmp2,nmo,nfirst,vorb,.true.,
     *                    bl(ifockMO),nslv,nbf)
      call sync_barrier
      do i=1,nval
         do j=1,i
            ij=ij+1
            kl=0
            call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                         ns, trans, signum)
            if (ns.ne.0) then
               cycle
            endif
            call para_recv_pack(islgid,TxCCSDMP2Req)
            call para_unpack_int(islid,1)
            call para_unpack_int(iresult,1)
            if (iresult.eq.1) then
c Result:
               call para_unpack_int(k,1)
               call para_unpack_int(l,1)
               call para_unpack_int(kl,1)
               call para_unpack_real(epair,1)
               call para_unpack_real(epairr,1)
               call para_unpack_real(xepair,1)
               call para_unpack_real(sepair,1)
               call para_unpack_real(sepairr,1)
               call para_unpack_real(sxepair,1)
               call para_unpack_real(sqij_sl,1)
               call para_unpack_real(xmaxele_sl,1)
               bl(iepairr+kl-1)=epairr
               bl(iepair+kl-1)=epair
               bl(ixepair+kl-1)=xepair
               energyr=energyr+bl(iepairr+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
               energy1=energy1+bl(iepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
               energy2=energy2+bl(ixepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
               senergy=senergy+sxepair
     *                      *dble(int_array(bl(ipairimages),kl))
               sqij=sqij+sqij_sl
     *                      *dble(int_array(bl(ipairimages),kl))
               if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
               if ((mp2_iter.eq.1.or.mp2last).and.iprint.gt.3)
     *            write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',k,l,
     *                                       bl(iepairr+kl-1)
               if (mp2last) then
                  energytmp=bl(iepairr+kl-1)
     &                      *dble(int_array(bl(ipairimages),kl))
                  if (k/=l) then
                     CIMPairEnergy(k,l)=energytmp/2.0D0
                     CIMPairEnergy(l,k)=energytmp/2.0D0
                  else
                     CIMPairEnergy(k,l)=energytmp
                  endif
               endif
            endif
c Work:
            call para_initsend
            iwork=1
            call para_pack_int(iwork,1)
            call para_pack_int(i,1)
            call para_pack_int(j,1)
            call para_pack_int(ij,1)
            call para_send_pack(islgid,TxCCSDMP2Job)
         enddo
      enddo
      do i=1,nslv
         call para_recv_pack(islgid,TxCCSDMP2Req)
         call para_unpack_int(islid,1)
         call para_unpack_int(iresult,1)
         if (iresult.eq.1) then
c Result:
            call para_unpack_int(k,1)
            call para_unpack_int(l,1)
            call para_unpack_int(kl,1)
            call para_unpack_real(epair,1)
            call para_unpack_real(epairr,1)
            call para_unpack_real(xepair,1)
            call para_unpack_real(sepair,1)
            call para_unpack_real(sepairr,1)
            call para_unpack_real(sxepair,1)
            call para_unpack_real(sqij_sl,1)
            call para_unpack_real(xmaxele_sl,1)
            bl(iepairr+kl-1)=epairr
            bl(iepair+kl-1)=epair
            bl(ixepair+kl-1)=xepair
            energyr=energyr+bl(iepairr+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            energy1=energy1+bl(iepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            energy2=energy2+bl(ixepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            senergy=senergy+sxepair
     *                      *dble(int_array(bl(ipairimages),kl))
            sqij=sqij+sqij_sl
     *                      *dble(int_array(bl(ipairimages),kl))
            if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
            if ((mp2_iter.eq.1.or.mp2last).and.kl.ne.0.and.iprint.gt.3)
     *          write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',k,l,
     *                                     bl(iepairr+kl-1)
            if (mp2last) then
               energytmp=bl(iepairr+kl-1)
     &                   *dble(int_array(bl(ipairimages),kl))
               if (k/=l) then
                  CIMPairEnergy(k,l)=energytmp/2.0D0
                  CIMPairEnergy(l,k)=energytmp/2.0D0
               else
                  CIMPairEnergy(k,l)=energytmp
               endif
            endif
         endif
c Work:
         call para_initsend
         iwork=0
         call para_send(iwork,islgid,TxCCSDMP2Job)
      enddo
      call fafClosem(ndiskmp2,0,info)
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
c
      if (diis.eq.1.or.(dabs(xmaxele).lt.5.0d-2.and.
     *                   ifirst_diis_iter.le.mp2_iter)) then
         diis=1
         call sync_barrier
         call para_bcast(diis,TxCCSDInit)
         call CCDiis(i,       j,          'calculate',  ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'mp2',bl(ilist),
     *               bl(ipairimages),nbf)
      else
         iwork=0
         call sync_barrier
         call para_bcast(iwork,TxCCSDInit)
      endif
c
      call sync_barrier
      call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,   ccsd,
     *                     vorb,  nmo,ivalpair,    ivirpair, ifp,
     *                     nsym,  nval,nvirt,      ij_unique,isympairs,
     *                     l_restart,ndisk_dump,nbf)
c
      if (iprint.gt.1) then
      write(6,'(A,I3)')     'MP2 iteration no:              ',mp2_iter
      if (.not.do_mp4) then
      write(6,'(A,F21.17)') 'Iterations - MP2 linear prev:  ',energy1
      write(6,'(A,F21.17)') 'Iterations - MP2 energy+resid: ',energyr
      write(6,'(A,F21.17)') 'Iterations - MP2 linear curr:  ',energy2
      else
      write(6,'(A,F21.17)') 'Iterations - MP4 energy:       ',energy2
      endif
      write(6,*) ' '
      write(6,'(A,F8.2)') 'Total time for iteration CPU:  ',t1-t0
      write(6,'(A,F8.2)') 'Total time for iteration ELA:  ',et1-et0
      write(6,'(A,E20.10E2)') 'Sum of squares of the residuum:',sqij
      write(6,'(A,E20.10E2)') 'The maximum residuum element:  ',xmaxele
      else
      xxenergy=energyr
      if (do_mp4) xxenergy=energy2
      xdelta=energyrold-energyr
      write(6,3232) mp2_iter,xxenergy,xdelta,xmaxele,sqij,(et1-et0)/6d1
      call flush(6)
      endif
 3232 format (I3,2F16.10,2E14.3E2,F15.1)
      call flush(6)
      if (dabs(energyrold-energyr).gt.Ethresh.or.xmaxele.gt.Wavethresh) 
     *  then
         iterate=1
         call para_barrier
         call para_bcast(iterate,TxCCSDIter)
         call para_barrier
         goto 222
      else
         if (.not.mp2last) then
            mp2last=.true.
            iterate=1
            call para_barrier
            call para_bcast(iterate,TxCCSDIter)
            call para_barrier
            goto 222
         endif
      endif
      mp2last=.false.
      iterate=0
      call para_initsend
      call para_barrier
      call para_bcast(iterate,TxCCSDIter)
      call para_barrier
      if (do_mp4) then
         xmp4_d=energy2
      else
         xmp2  =energyr
         scsmp2=senergy
      endif
      if (mp2) then
         ecenmp2=0.0D0; CIMOrbEnergy=0.0D0
         do i=1,nval
            do j=1,nval
               CIMOrbEnergy(i)=CIMOrbEnergy(i)+CIMPairEnergy(i,j)
            enddo
         enddo
         do i=1,ncen
            ecenmp2=ecenmp2+CIMOrbEnergy(i)
         enddo
      endif
      if (mp2.or.do_mp4) then
         goto 999
      endif
c     END loop
c if localization is on, calculate list of strong pairs
      call flush(6)
      call CalcList(bl(ilist),nfirst,nval,nmo,omit,inumber,
     *              weaksepar,distsepar,bl(iepairr),iEnergy)
 126  continue ! restart!
      if (l_restart) then
      call general_read1(ndisk_dump,bl(ilist),npairs*8/intsize)
      call general_read1(ndisk_dump,inumber,8/intsize)
      endif
      if (iprint.gt.1)
     *write(6,*) 'There are ',inumber,' strong pairs.'
      call para_initsend
      call para_pack_int(bl(ilist),npairs)
      call para_pack_int(inumber,1)
      call para_bcast_pack(TxCCSDInit)
c
      call dynamic_matdef('WernerX','q',idimen,idimen)
      iX=mataddr('WernerX')
c
      xiterenergy=energyr
      if (mp3) xiterenergy=0.0d0
      iteration=0
      varenergyold=0.0d0
      diis=0
      call CCDiis(i,       j,          'reset',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      if (iprint.gt.2) write(6,*) "CCDIIS done"
      call flush(6)
c*******************************************************************
      if ((mp4.or.l_restart).and.l_triples) then
c     call mas_trip(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),res,
c    *             ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *             singles.and.cc.and..not.ccsd,af,nslv)
      call elapsec(xt0)
      time_keep=xt0
      if (mygid.ne.0) STOP 'mygid.ne.0 ??!!'
      call sort_amplit(nval,idimen,ndisk_a,af,mygid,
     *                 ichar_count,bl(ivtable),bl(im_table),bl(ijtable),
     *                 bl(ivrevtable),amp_ratio)
      if (iprint.gt.2)
     *write(6,*) "sort_amplit done"
      call flush(6)
      call elapsec(xt1)
      if (iprint.gt.2) then
      write(6,'(A,F10.2)')'Amplitudes sort time:       ',xt1-xt0
      call flush(6)
      endif
      call elapsec(xt0)
ctmp? call fafClosem(ndisktt,0,iresult)
      call master_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                 ndisktt,
     *                 ncf,lbintt,thresh,byt8,nmo,vorb,af,mygid,
     *                 nslv,ext3_ratio)
      call fafClosem(ndisktrtt,0,iresult)
      if (iprint.gt.2)
     *write(6,*)'RATIO: ',ext3_ratio
ccs      call sleep(10)
      call make_3ext_pairs_mas(idimen,nval,af,ndisk_ie,ndisk_ie1,
     *                         nslv)
      call elapsec(xt1)
      if (iprint.gt.2)
     *write(6,'(A,F10.2)')'3ext integrals sort time:   ',xt1-xt0
      call elapsec(xt0)
      call sort_3int(nval,idimen,ndisk_ii,af,mygid)
      call elapsec(xt1)
      if (iprint.gt.2)
     *write(6,'(A,F10.2)')'3int integrals sort time:   ',xt1-xt0
      call dynamic_show_free(mem)
      mem=mem-nval*nval*nval-3*nval*nval
      call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                     nslv,npass,isize)
      call elapsec(xt0)
      call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *               ncf,lbinx,thresh,byt8,nmo,vorb,af,mygid,
     *               isize,npass)
      call elapsec(xt1)
      if (iprint.gt.2) then
      write(6,'(A,F10.2)')'Kext integrals sort time:   ',xt1-xt0
      write(6,'(A,F10.2)')'Total sort time:            ',xt1-time_keep
      endif
      call flush(6)
ccs      call sleep(5)
      call elapsec(xt0)
      call dynamic_getmem(20,itimes)
      call new_tri_mas(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,ndisk_ix,
     *                       res,bl(iepsi+ncore),bl(itimes),20,ccsd,
     *                       qcisd,af,nslv,ext3_ratio,amp_ratio,npass,
     *                       isize,iprint,energys,energyd,itrstart)
      call elapsec(xt1)
      write(6,'(A,F10.1)') 'Triples elapsed time: ',(xt1-xt0)/6d1
      if (iprint.gt.2) then
      write(6,'(A,F20.15,A,F10.2)')
     *  'Triples energy correction: ',res,',   elapsed time: ',xt1-xt0
      write(6,*)'W zero & sort:           ',bl(itimes+0)
      write(6,*)'W build (total):         ',bl(itimes+1)
      write(6,*)'Mult. over virt space:   ',bl(itimes+2)
      write(6,*)'Mult. over occ. space:   ',bl(itimes+3)
      write(6,*)'Reading of ampl. & int:  ',bl(itimes+4)
      write(6,*)'Amplitudes:              ',bl(itimes+5)
      write(6,*)'3int integrals:          ',bl(itimes+6)
      write(6,*)'3ext integrals:          ',bl(itimes+7)
      write(6,*)'Kext integrals:          ',bl(itimes+9)
      write(6,*)'Energy:                  ',bl(itimes+8)
      write(6,*)'Wabc relocate:           ',bl(itimes+10)
      write(6,*)'Inte relocate:           ',bl(itimes+11)
      write(6,*)'3int relocate:           ',bl(itimes+12)
      write(6,*)'Total virt part:         ',bl(itimes+13)
      write(6,*)'Tables reading:          ',bl(itimes+14)
      call flush(6)
      endif
      if (l_restart) goto 999
      endif
c*******************************************************************
      write(6,'(A)')
      write(6,'(A)')
      write(6,'(A)') "     CC module:"
      write(6,'(A)')
      write(6,'(A)')
     * "Iter:  Energy:        Delta E           "//
     * "Max resid:    Err sq:      Elapsed time"

      iterate=1
 666  continue
      call f_lush(6)
      do i=1,itim_no
         ccsdcpu(i)=0.0d0
         ccsdela(i)=0.0d0
      enddo
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      call para_initsend
      if (iterate==1) then
         call para_barrier
         call para_pack_real(xmaxele,1)
         call para_pack_int(iterate,1)
         call para_bcast_pack(TxCCSDIter1)
         call para_barrier
      endif
      iteration=iteration+1
      tot_xmax=xmaxele
      xmaxele=0.0d0
      if (iteration.eq.2) tot_xmax=0.1d0
      call set_thresh(tot_xmax,xthresh_EEO)
      if (cep0) xiterenergy=0.0d0
      call elapsec(eistart)
      call secund(sistart)
      int_kind='e' ! EE operator
      strace=0.0d0
      sconstr=0.0d0
      sort=0.0d0
      extr=0.0d0
      const=0.0d0
      call secund(t0)
      call elapsec(et0)
      if (.not.do_mp4) then
         ndiske=-1
         lrece =-1
         lbine =-1
         call MAS_EEO_INT(ncs,    ncf,     bl(ictr),   nval,  nmo,
     *                    nfirst, nlast,   thresh,
     *                    vorb,   ndisktre, npairs, nslv,iprint,nbf)
      endif
c
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(1)=t1-t0        ! EEO
      ccsdela(1)=et1-et0      ! EEO
c
      if (ccsd) then
         int_kind='tx' ! TEIO operator
         call secund(t0)
         call elapsec(et0)
         call master_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,core,
     &                   xmaxdisk,ndisktx,bl(irecadrtx),lrectx,lbintx,
     &                   int_kind,nodisk,byt8,small,vorb,ndisktrtx,
     *                   bl(ilist),bl(isympairs),ij_unique,iprint,nbf)
         call secund(t1)
         call elapsec(et1)
         ccsdcpu(2)=t1-t0        ! TEIO tx
         ccsdela(2)=et1-et0      ! TEIO tx
c
         int_kind='tc' ! TEIO operator
         call secund(t0)
         call elapsec(et0)
         call master_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,core,
     &                   xmaxdisk,ndisktc,bl(irecadrtc),lrectc,lbintc,
     &                   int_kind,nodisk,byt8,small,vorb,ndisktrtc,
     *                   bl(ilist),bl(isympairs),ij_unique,iprint,nbf)
         call secund(t1)
         call elapsec(et1)
         ccsdcpu(3)=t1-t0        ! TEIO tc
         ccsdela(3)=et1-et0      ! TEIO tc
c
      endif
c
      call secund(t0)
      call elapsec(et0)
      do i=1,5
         elaps(i)=0.0d0
      enddo
      call CreateAlphaDisk(nval,ndiskalpha,af)
      call sync_barrier
      call para_bcast(ndiskalpha,TxCCSDInit)
      if (cc.or.ccsd.or.do_mp4) then
         call CCalphamast(ncf,nval,bl(irecadrx),npairs,ndiskx,lbinx,
     *                    thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                    ndiskalpha,af,nbf)
         call sync_barrier
         call para_reduce(elaps,5,TxCCSDRedu)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(4)=t1-t0        ! CCalphaonce+supplement
      ccsdela(4)=et1-et0      ! CCalphaonce+supplement
      call secund(t0)
      call elapsec(et0)
      call SupplementAlpha_BetaBuild(nval,ndiskalpha,bl(ibeta),ncf,
     *                              nfirst,nmo,vorb,cc,ccsd,bl(ifockMO),
     *                              kijklndisk,af,do_mp4,nbf)
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
      call para_initsend
      call dynamic_unlock(bl(ibeta),ii)
      call para_bcast_real(bl(ibeta),nval*nval,TxCCSDInit)
      if (iprint.gt.2) then
         write(6,9)'SupplementAlpha_BetaBuild, ELAPS: ',et1-et0
         write(6,9)'SupplementAlpha_BetaBuild, CPU:   ',t1-t0
         call flush(6)
      endif
      call dynamic_lock(bl(ibeta),ii)
c
      call secund(t0)
      call elapsec(et0)
      call B41a_G_Master(ncf,    nval,        norecalc, cc,   ccsd,
     *                   nfirst, bl(ifockMO), nmo,      vorb, ndiskG41a,
     *                   .true., ndiskalpha,  bl(ibeta),nbf)
c     close(ndiskalpha,STATUS='delete')
c
c
      call fafClosem(ndiskalpha,0,info)
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(5)=t1-t0        ! Builder41a_G
      ccsdela(5)=et1-et0      ! Builder41a_G
c
      call secund(t0)
      call elapsec(et0)
      if (cc.or.ccsd.or.do_mp4) then
         call CCAOncemast(bl(irecadrx),npairs,nval,    ndiskx,ncf,
     *                    lbinx,       thresh,norecalc,byt8,  ioverlap,
     *                    nmo,         vorb,  nslv,    nbf)
         call sync_barrier
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(6)=t1-t0        ! CCAOnce
      ccsdela(6)=et1-et0      ! CCAOnce
      ij=0
      xlinear_energy=0.0d0
c
c
      call secund(t0)
      call elapsec(et0)
c This is the s and r loop, START
      if (singles.or.ccsd.or.do_mp4) then
         call matcopy('corefock','work1')
c
         if (cc.or.ccsd) then
            call CCA(ncf,iwork2,nmo,vorb,nbf)
            call matpose('work2')
            call matadd1('work2',-1.0d0,'work1')
         endif
c
         call Fock_Vector_Pointer(ifock)
         call pointS(iiS)
         call tfer(bl(ifock),bl(iiS),nval*idimen)
c
         if (.not.do_mp4) then
            call pointersingles(ifock)
            call matconn('sing','r',idimen,nval,ifock)
            call matconn('S_vect','r',idimen,nval,iiS)
            call matmmul2('work1','sing','S_vect','n','n','a')
            call matdisc('S_vect')
            call matdisc('sing')
         endif
c
         call elapsec(esi0)
         call secund(csi0)
         call Tl_gen_mast(ncf,nval,nmo,vorb,ioverlap,iiS,nslv,nbf)
         call elapsec(esi1)
         call secund(csi1)
         ccsdela(15)=esi1-esi0
         ccsdcpu(15)=csi1-csi0
         call elapsec(esi0)
         call secund(csi0)
         call sync_barrier
         call EEO_vector_master(nval,bl(irecadre),npairs,ndiske,ncf,
     *                          lbine,thresh,nfirst,nlast,byt8,
     *                          nmo,vorb,nslv,iiS,nbf)
         call elapsec(esi1)
         call secund(csi1)
         call sync_barrier
         ccsdela(16)=esi1-esi0 ! EEO_vector_extractor
         ccsdcpu(16)=csi1-csi0 ! EEO_vector_extractor
         call para_bcast_real(bl(iiS),idimen*nval,TxCCSDInit)
c  R loop
         call Fock_Vector_Pointer(ifock)
         call pointRR(iiR)
         call tfer(bl(ifock),bl(iiR),nval*idimen)
         if (cc.or.ccsd) then
            call Lt_gen_mast(idimen,  nval,  npairs,  ncf,  vorb,
     *                       iiR,nslv)
            call sync_barrier
         endif
         call para_bcast_real(bl(iiR),idimen*nval,TxCCSDInit)
      endif
c
c This was the s loop, STOP
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(8)=t1-t0        ! singles s, r
      ccsdela(8)=et1-et0      ! singles s, r
c
      call secund(t0)
      call elapsec(et0)
      call buildX(bl(irecadrtx),npairs,ndisktx,lbintx,thresh,byt8,
     *           bl(irecadrtc),ndisktc,lbintc,ioverlap,nval,ncf,cc,ccsd,
     *            ifockAO,iX,nmo,vorb,do_mp4,nbf)
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(7)=t1-t0        ! buildX
      ccsdela(7)=et1-et0      ! buildX
      call sync_barrier
      call para_initsend
      call para_pack_real(bl(iX),idimen*idimen)
      call para_pack_int(ndiskg41a,1)
      call para_pack_real(xiterenergy,1)
      call para_bcast_pack(TxCCSDInit)
c
c Calculate ichunk, jchunk:
c
      call secund(t0)
      call elapsec(et0)
      if (calchunks)
     *      call calculate_ij_chunk(ichunk,jchunk,nslv,idimen,nval)
      if (iprint.gt.2) write(6,*)'ichunk,jchunk: ',ichunk,jchunk
      call flush(6)
 133  FORMAT (4I3,F20.10)
      call prepare_CCYZ_file(ncf,nmo,vorb,ndiskYZ,.true.,nbf)
      call sync_barrier
      if (iprint.gt.3) write(6,*) 'I am master, ndiskYZ: ',ndiskYZ
      call flush(6)
      call para_bcast(ndiskYZ,TxCCSDInit)
      istop=0
      icalc_pairs=0
      reduced=.false.
      call flush(6)
      do i=1,5
         elaps(i)=0.0d0
      enddo
      do
         if (istop.eq.nval) exit
         if (calchunks)  call check_ijchunk_reduce(icalc_pairs,nval,
     &                        ichunk,jchunk,reduced)
        istart=istop+1
        istop =istart+ichunk-1
        if (istop.gt.nval) istop=nval
        if (istart.gt.nval) STOP 'Error with istart!'
        if (istart.gt.istop) STOP 'Error with istart & istop!'
        jstop=0
        do
           if (jstop.eq.nval) exit
           jstart=jstop+1
           jstop =jstart+jchunk-1
           if (jstop.gt.nval) jstop=nval
           if (jstart.gt.nval) STOP 'Error with jstart!'
           if (jstart.gt.jstop) STOP 'Error with jstart & jstop!'
           icalc_pairs=icalc_pairs+(jstop-jstart+1)*(istop-istart+1)
cc          if (iprnt.ge.2) write(91,*)
cc     *         'istart,istop,jstart,jstop: ',istart,istop,jstart,jstop
cc          call flush(6)
           call para_recv_pack(islgid,TxCCSDReq)
           call para_unpack_int(islid,1)
           call para_unpack_int(iresult,1)
           if (iresult.eq.1) then
c Result:
              call para_unpack_int(kstart,1)
              call para_unpack_int(kstop,1)
              call para_unpack_int(lstart,1)
              call para_unpack_int(lstop,1)
           endif
c Work:
           call para_initsend
           iwork=1
           call para_pack_int(iwork,1)
           call para_pack_int(istart,1)
           call para_pack_int(istop,1)
           call para_pack_int(jstart,1)
           call para_pack_int(jstop,1)
           call para_send_pack(islgid,TxCCSDJob)
         enddo
      enddo
      do i=1,nslv
         call para_recv_pack(islgid,TxCCSDReq)
         call para_unpack_int(islid,1)
         call para_unpack_int(iresult,1)
         if (iresult.eq.1) then
c Result:
            call para_unpack_int(kstart,1)
            call para_unpack_int(kstop,1)
            call para_unpack_int(lstart,1)
            call para_unpack_int(lstop,1)
         endif
c Work:
         call para_initsend
         iwork=0
         call para_pack_int(iwork,1)
         call para_send_pack(islgid,TxCCSDJob)
      enddo
      call sync_barrier
      call para_reduce(elaps,5,TxCCSDRedu)
      call para_barrier
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(9)=t1-t0        ! YZ
      ccsdela(9)=et1-et0      ! YZ
c
c
c Q parts:
c
      call secund(t0)
      call elapsec(et0)
      call prepare_Qparts_file(ncf,nmo,vorb,ndiskQ,af,nbf)
      call sync_barrier
      call para_bcast(ndiskQ,TxCCSDInit)
      if (calchunks)
     *       call calculate_ij_chunk(ichunk,jchunk,nslv,idimen,nval)
      istop=0
      icalc_pairs=0
      reduced=.false.
      call flush(6)
      do
         if (istop.eq.nval) exit
         write(91,*) 'icalc_pairs: ',icalc_pairs
         if (calchunks) call check_ijchunk_reduce(icalc_pairs,nval,
     &                       ichunk,jchunk,reduced)
         istart=istop+1
         istop =istart+ichunk-1
         if (istop.gt.nval) istop=nval
         if (istart.gt.nval) STOP 'Error with istart!'
         if (istart.gt.istop) STOP 'Error with istart & istop!'
         jstop=0
         do
            if (jstop.eq.nval) exit
            jstart=jstop+1
            jstop =jstart+jchunk-1
            if (jstop.gt.nval) jstop=nval
            if (jstart.gt.nval) STOP 'Error with jstart!'
            if (jstart.gt.jstop) STOP 'Error with jstart & jstop!'
            icalc_pairs=icalc_pairs+(jstop-jstart+1)*(istop-istart+1)
            write(91,*)
     *           'istart,istop,jstart,jstop: ',istart,istop,jstart,jstop
            call flush(6)
            call para_recv_pack(islgid,TxCCSDReqQ)
            call para_unpack_int(islid,1)
            call para_unpack_int(iresult,1)
            if (iresult.eq.1) then
c Result:
               call para_unpack_int(kstart,1)
               call para_unpack_int(kstop,1)
               call para_unpack_int(lstart,1)
               call para_unpack_int(lstop,1)
            endif
c Work:
            call para_initsend
            iwork=1
            call para_pack_int(iwork,1)
            call para_pack_int(istart,1)
            call para_pack_int(istop,1)
            call para_pack_int(jstart,1)
            call para_pack_int(jstop,1)
            call para_send_pack(islgid,TxCCSDJobQ)
         enddo
      enddo
      do i=1,nslv
         call para_recv_pack(islgid,TxCCSDReqQ)
         call para_unpack_int(islid,1)
         call para_unpack_int(iresult,1)
         if (iresult.eq.1) then
c Result:
            call para_unpack_int(kstart,1)
            call para_unpack_int(kstop,1)
            call para_unpack_int(lstart,1)
            call para_unpack_int(lstop,1)
         endif
c Work:
         iwork=0
         call para_send(iwork,islgid,TxCCSDJobQ)
      enddo
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(10)=t1-t0        ! Q parts
      ccsdela(10)=et1-et0      ! Q parts
c
      call secund(t0)
      call elapsec(et0)

      if (iterate==0) then
         allocate(exc_QCMO(nval*(nval+1)/2,nvirt,nvirt))
         allocate(coef_QCMO(nval*(nval+1)/2,nvirt,nvirt))
         allocate(resi_QCMO(nval*(nval+1)/2,nvirt,nvirt))
         exc_QCMO=0.0D0; coef_QCMO=0.0D0; resi_QCMO=0.0D0
         allocate(tmpexc(nvirt,nvirt))
         allocate(tmpcoef(nvirt,nvirt))
         allocate(tmpresi(nvirt,nvirt))
         tmpexc=0.0D0; tmpcoef=0.0D0; tmpresi=0.0D0
      endif

      energy=0.0d0
      energyr=0.0d0
      square=0.0d0
      totnorm=0.0d0
      ij=0
      do i=1,nval
         do j=1,i
            ij=ij+1
            call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                         ns, trans, signum)
            if (ns.ne.0) cycle
            kl=0
            call para_recv_pack(islgid,TxCCSDReq1)
            call para_unpack_int(islid,1)
            call para_unpack_int(iresult,1)
            if (iresult.eq.1) then
               call para_unpack_int(k,1)
               call para_unpack_int(l,1)
               call para_unpack_int(kl,1)
               call para_unpack_real(epair,1)
               call para_unpack_real(epairr,1)
               call para_unpack_real(xepair,1)
               call para_unpack_real(sqij,1)
               call para_unpack_real(xmaxele_sl,1)
               call para_unpack_real(xnormij,1)
               bl(iepairr+kl-1)=epairr
               bl(iepair+kl-1)=epair
               bl(ixepair+kl-1)=xepair
               xlinear_energy=xlinear_energy+bl(ixepair+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
               energyr=energyr+bl(iepairr+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
               square=square+sqij
     *                   *dble(int_array(bl(ipairimages),kl))
               energy=energy+bl(iepair+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
               totnorm=totnorm+xnormij
     *                   *dble(int_array(bl(ipairimages),kl))
               if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
            endif
c Work:
            call para_initsend
            iwork=1
            call para_pack_int(iwork,1)
            call para_pack_int(i,1)
            call para_pack_int(j,1)
            call para_pack_int(ij,1)
            if (cep2) then
               xiterenergy=bl(iepairr+ij-1)
               call para_pack_real(xiterenergy,1)
            endif
            call para_send_pack(islgid,TxCCSDJob1)
            if (iterate==0) then
               call para_recv_pack(islgid,TxCIMCCSDInt)
               call para_unpack_int(k,1)
               call para_unpack_int(l,1)
               call para_unpack_real(tmpexc,nvirt*nvirt)
               call IntExpr(nval,nvirt,k,l,tmpexc,exc_QCMO)
               call para_unpack_real(tmpcoef,nvirt*nvirt)
               call IntExpr(nval,nvirt,k,l,tmpcoef,coef_QCMO)
               call para_unpack_real(tmpresi,nvirt*nvirt)
               call IntExpr(nval,nvirt,k,l,tmpresi,resi_QCMO)
            endif
         enddo
      enddo
      do i=1,nslv
         kl=0
         call para_recv_pack(islgid,TxCCSDReq1)
         call para_unpack_int(islid,1)
         call para_unpack_int(iresult,1)
         if (iresult.eq.1) then
            call para_unpack_int(k,1)
            call para_unpack_int(l,1)
            call para_unpack_int(kl,1)
            call para_unpack_real(epair,1)
            call para_unpack_real(epairr,1)
            call para_unpack_real(xepair,1)
            call para_unpack_real(sqij,1)
            call para_unpack_real(xmaxele_sl,1)
            call para_unpack_real(xnormij,1)
            bl(iepairr+kl-1)=epairr
            bl(iepair+kl-1)=epair
            bl(ixepair+kl-1)=xepair
            xlinear_energy=xlinear_energy+bl(ixepair+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
            energyr=energyr+bl(iepairr+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
            square=square+sqij
     *                   *dble(int_array(bl(ipairimages),kl))
            energy=energy+bl(iepair+kl-1)
     *                   *dble(int_array(bl(ipairimages),kl))
            totnorm=totnorm+xnormij
     *                   *dble(int_array(bl(ipairimages),kl))
            if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
         endif
c No work:
         iwork=0
         call para_send(iwork,islgid,TxCCSDJob1)
      enddo
      
      if (iterate==0) goto 888
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(11)=t1-t0        ! Residuals
      ccsdela(11)=et1-et0      ! Residuals
c
c calculate variational (I hope) CI energy:
      varenergy=(energyr+totnorm*xiterenergy)/(1.0d0+totnorm)
c
c This is the singles residuum "loop", START
      call secund(t0)
      call elapsec(et0)
      if (do_mp4) then
         call do_mp4_beta(bl(ifockMO),bl(ibeta),nval,ncf,nfirst)
         call dynamic_matdef('spwork1','q',ncf,ncf)
         call dynamic_matdef('spwork2','q',ncf,ncf)
         call dynamic_matdef('spwork3','q',ncf,ncf)
         ispwork1=mataddr('spwork1')
      endif
 665  continue
      sing_max=0.0d0
      if (do_mp4) call matzero('spwork1')
      if (singles.or.ccsd.or.do_mp4) then
         ssq=0.0d0
         iwork=1
         call sync_barrier
         call para_bcast(iwork,TxCCSDInit)
c
         call pointS(iiS)
         call dynamic_matdef('s_resid','r',idimen,nval)
         is_resid=mataddr('s_resid')
         call pointersingles(isingles)
         call pointernewsingles(inewsingles)
         call matconn('full_sing','r',idimen,nval,isingles)
         call tfer(bl(iiS),bl(is_resid),idimen*nval)
         call Tr_gen_mast(idimen,nval,ncf,npairs,vorb,
     *                       ioverlap,nslv,is_resid)
         call beta_t_sumator1(ncf,idimen,nval,ibeta,vorb,ioverlap,
     *                        is_resid)
         if ((.not.cc).and.(.not.ccsd).and.(.not.do_mp4)) then
            call cisd_energy_add(ncf,idimen,nval,xiterenergy,vorb,
     *                           ioverlap,is_resid)
         endif
         if (do_mp4) then
            call matmmul2('fockAO','full_sing','s_resid','n','n','a')
         endif
         call matdisc('full_sing')
         call UpdateSingles1(isingles,is_resid,nmo,nfirst,ncf,npairs,
     &                       nval,bl(iepsi),shift,ssq,vorb,sing_max,
     &                       ndiskdiisr,ndiskdiisc,bl(ilist),
     &                       bl(ipairimages),inewsingles,bl(iorevtable),
     &                       bl(ivrevtable),docansym,t1diagnost,nbf)
c
         if (do_mp4) then
            do i=1,nval
               isin=inewsingles+(i-1)*idimen
               call storesingles(i,bl(isin),bl(ispwork1),ncf,idimen,
     &                           nfirst)
            enddo
         endif
         call dynamic_matrem('s_resid')
         if (do_mp4) then
            call matconn('fockMO','q',ncf,ncf,ifockMO)
            total_singl_en=0.0d0
            call matmmult('spwork1','fockMO','spwork2')
            if (vorb) then
               call matprodtr('spwork2','spwork1',part_singl)
            else
               write(6,*)'AO orbitals are not working yet'
               STOP    'AO orbitals are not working yet'
               call matmmul2('spwork2','spwork1','spwork3','n','t','n')
               call matprodtr('spwork3','overlap',part_singl)
            endif
            total_singl_en=total_singl_en+part_singl
            call matmmult('fockMO','spwork1','spwork2')
            if (vorb) then
               call matprodtr('spwork2','spwork1',part_singl)
            else
               write(6,*)'AO orbitals are not working yet'
               STOP    'AO orbitals are not working yet'
               call matmmul2('spwork2','spwork1','spwork3','t','n','n')
               call matprodtr('spwork3','overlap',part_singl)
            endif
            total_singl_en=total_singl_en-part_singl
            total_singl_en=-2.0d0*total_singl_en
            write(6,*)'Total singles MP4 energy: ', total_singl_en
            call initsingles_local(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
            call matdisc('fockMO')
         endif
         if (iprint.gt.2) write(6,*)'Sum of singles squares: ', ssq
      endif
      if (do_mp4.and.sing_max.gt.Wavethresh) goto 665
      if (singles.or.ccsd.or.do_mp4) then
         iwork=0
         call sync_barrier
         call para_bcast(iwork,TxCCSDInit)
      endif
      xmp4_s=total_singl_en
      if (do_mp4) then
         call dynamic_matrem('spwork3')
         call dynamic_matrem('spwork2')
         call dynamic_matrem('spwork1')
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(13)=t1-t0        ! Singles Residuals
      ccsdela(13)=et1-et0      ! Singles Residuals
      call sync_barrier
c This is the singles residuum loop, STOP
c
  8   FORMAT (A)
 11   FORMAT ((A,F25.15))
 12   FORMAT ((A,F25.15))
      if (iprint.gt.1) then
      write(6,8) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,'(A,I4)') '          The results for iteration:',iteration
      write(6,8) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,11) 'Sum of squares of residuums:         ',square
      write(6,11) 'Maximum residuum element (abs value):',xmaxele
      write(6,12) 'T1 diagnostic: ',t1diagnost
      endif
      if (mp3.or.do_mp4) then 
         if (mp3.and.iprint.gt.1) 
     *      write (*,'(A,F25.16)') 'MP3 Energy is: ', energy
         if (do_mp4) then 
            xmp4_q=energy
         else
            xmp3=energy
         endif
      else
         if (iprint.gt.1) then
      if (cep0) write(6,'(A)') '* * * * * * * * * * * * * * * * CEPA-0'
     *//' Energies: * * * * * * * * * * * * * * *'
      if (cep2) write(6,'(A)') '* * * * * * * * * * * * * * * * CEPA-2'
     *//' Energies: * * * * * * * * * * * * * * *'
      write(6,35) energy
      write(6,36) energyr
      write(6,37) xlinear_energy
      if (iprint.gt.2) then
        write(6,39) varenergy
        write(6,41) varenergy-varenergyold
        write(6,40) totnorm
      endif
      write(6,38) energyr-energyrold
      endif
      endif
 32   FORMAT('Energy:  ',F21.16)
 33   FORMAT('Energy with residuum:  ',F21.16)
 34   FORMAT('MP2 energy:  ',F21.16)
 35   FORMAT('Energy from coefficients from n-1 cycle:  ',F21.16)
 36   FORMAT('Energy from coeff. from n-1 cycle+resid:  ',F21.16)
 37   FORMAT('Energy from coefficients from n   cycle:  ',F21.16)
 38   FORMAT('Energy(n) - energy(n-1) (from quad. form):',F21.16)
 39   FORMAT('Variational energy (only disk ngss!):     ',F21.16)
 40   FORMAT('Wave function norm-1:                     ',F21.16)
 41   FORMAT('Present var. energy-previous var. energy: ',F21.16)
      DeltaE=dabs(energyr-energyrold)
      DeltaE1=    energyr-energyrold
      if (iteration.le.1.and.l_restart)
     *    call general_read1(ndisk_dump,DeltaE,8)
c     DeltaE=dabs(xlinear_energy-energy)
      if (.not.(mp3.or.do_mp4)) then
      energyrold=energyr
      varenergyold=varenergy
      if (variational) then
      xiterenergy=varenergy ! for CID quadratic corrected
      if (iprint.gt.1) then
      write(6,'(A)')'Variational energy will be used for the'
     *          //' next iteration.'
      write(6,'(A)') 'Remember! In fact this energy is not variational'
     *           //' when Gauss-Seidel'
      write(6,'(A)')'                             algorithm was used'
      endif
      else if (linear) then
      xiterenergy=xlinear_energy   ! for CID linear
      if (iprint.gt.1) then
      write(6,'(A)')'Linear energy will be used for the next iteration.'
      endif
      else
      xiterenergy=energyr ! for CID quadratic (default)
      if (iprint.gt.1) then
      write(6,'(A)') 'Quadratic (non-variational!) energy will be used '
     *           //'for the next iteration.'
      endif
      endif
c DIIS!
      call secund(t0)
      call elapsec(et0)
      if (diis.eq.1.or.(dabs(xmaxele).lt.5.0d-2.and.
     *                   ifirst_diis_iter.le.iteration)) then
      if (diis.ne.1.and.iprint.gt.2) 
     *         write(6,*)'Calculation of diis matrix started.'
      diis=1
      call sync_barrier
      call para_bcast(diis,TxCCSDInit)
      if (iprint.gt.2) write(6,*) 'Sent diis info: ',diis
      call CCDiis(i,       j,          'calculate',  ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'ccc',bl(ilist),
     *            bl(ipairimages),nbf)
      else
      iwork=0
      call sync_barrier
      call para_bcast(iwork,TxCCSDInit)
      if (iprint.gt.2) write(6,*) 'Sent diis info: ',iwork
      call flush(6)
      endif
c
c
      if (dabs(xmaxele).lt.diis_thres) diis_less_thr=.true.
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(14)=t1-t0        ! DIIS
      ccsdela(14)=et1-et0      ! DIIS
c DIIS!
c---------------------------------------------------------------------
      endif ! if not mp3.or.do_mp4
      call sync_barrier
      call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,   ccsd,
     *                     vorb,  nmo,ivalpair,    ivirpair, ifp,
     *                     nsym,  nval,nvirt,      ij_unique,isympairs,
     *                     l_restart,ndisk_dump,nbf)
      call initsingles_master(nval,ncf,singles.or.mp4,nmo,vorb,
     *                              l_restart,ndisk_dump,nbf)
c---------------------------------------------------------------------
      call elapsec(eistop)
      call secund(sistop)
      ccsdcpu(12)=sistop-sistart      ! Tot. iteration time
      ccsdela(12)=eistop-eistart      ! Tot. iteration time
c---------------------------------------------------------------------
      slave_cpu=0.0d0
      call sync_barrier
      call para_reduce(slave_cpu,1,TxCCSDRedu)
      call para_barrier
      efficiency=(slave_cpu+ccsdcpu(12))/ccsdela(12)
      if (iprint.gt.2) then
      write(6,'(A,F8.2)')'Efficiency for this iteration: ',efficiency
      write(6,'(A,I5)')  'Number of slaves working:      ',nslv
      endif
c---------------------------------------------------------------------
      call secund(t0)
      call elapsec(et0)
      if (.not.last_dump_only) then
      call dump_data(dump,     ndisk_dump,ncf,    nval,   npairs,
     *               idimen,   icano,     iepsi,  ifockAO,ifockMO,
     *               icorefock,ilist,     inumber,singles,DeltaE)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(17)=t1-t0        ! dump
      ccsdela(17)=et1-et0      ! dump
c
      if (iprint.gt.1) call print_CCSD_results(ccsdela,ccsdcpu,itim_no)
      write(6,3232) iteration,energyr,DeltaE1,xmaxele,square,
     *              ccsdela(12)/6d1
      call flush(6)
c
c     STOP 'xx'
      call sync_barrier
      if ((xmaxele.gt.Wavethresh .or. DeltaE.gt.Ethresh).and..not.mp3
     *     .and..not.do_mp4) then
         if (iteration.ge.maxiter) then
            write(6,*)'            * * * * Maximum number of'//
     *                ' iterations reached * * * *'
            write(6,*)'            * * * *             NO '//
     *                'CONVERGENCE           * * * *'
         else
            goto 666
         endif
      endif

      call para_initsend
      iterate=0
      call para_barrier
      call para_pack_real(xmaxele,1)
      call para_pack_int(iterate,1)
      call para_bcast_pack(TxCCSDIter1)
      call para_barrier
      if (last_dump_only) then
         call secund(t0)
         call elapsec(et0)
         call dump_data(dump,     ndisk_dump,ncf,    nval,   npairs,
     *                  idimen,   icano,     iepsi,  ifockAO,ifockMO,
     *                  icorefock,ilist,     inumber,singles,DeltaE)
      endif
      call secund(t1)
      call elapsec(et1)
      if (iprint.gt.2) then
         write(6,'(A,F10.1)')'Dump CPU time:   ',t1-t0        ! dump
         write(6,'(A,F10.1)')'Dump elaps time: ',et1-et0      ! dump
      endif
      goto 666 !Converged and need to collect integrals
 70   continue
      
888   allocate(exc_LMO(nval,nval,idimen,idimen))
      allocate(coef_LMO(nval,nval,idimen,idimen))
      allocate(resi_LMO(nval,nval,idimen,idimen))
      exc_LMO=0.0D0; coef_LMO=0.0D0; resi_LMO=0.0D0
C      do i=1,nval
C         do j=1,nval
C            do k=1,idimen
C               do l=1,idimen
C                  do ii=1,nval
C                     if (ii>=j) then
C                        ij=ii*(ii-1)/2+j
C                        eqcmo=exc_QCMO(ij,k,l)
C                        cqcmo=coef_QCMO(ij,k,l)
C                        rqcmo=resi_QCMO(ij,k,l)
C                     else
C                        ij=j*(j-1)/2+ii
C                        eqcmo=exc_QCMO(ij,l,k)
C                        cqcmo=coef_QCMO(ij,l,k)
C                        rqcmo=resi_QCMO(ij,l,k)
C                     endif
C                     exc_LMO(i,j,k,l)=exc_LMO(i,j,k,l)+eqcmo*TX(ii,i)
C                     coef_LMO(i,j,k,l)=coef_LMO(i,j,k,l)+cqcmo*TX(ii,i)
C                     resi_LMO(i,j,k,l)=resi_LMO(i,j,k,l)+rqcmo*TX(ii,i)
C                  enddo
C               enddo
C            enddo
C         enddo
C      enddo

      do k=1,idimen
         do l=1,idimen
            allocate(eqcmo2(nval,nval),cqcmo2(nval,nval))
            allocate(rqcmo2(nval,nval))
            do j=1,nval
               do ii=1,nval
                  if (ii>=j) then
                     ij=ii*(ii-1)/2+j
                     eqcmo2(ii,j)=exc_QCMO(ij,k,l)
                     cqcmo2(ii,j)=coef_QCMO(ij,k,l)
                     rqcmo2(ii,j)=resi_QCMO(ij,k,l)
                  else
                     ij=j*(j-1)/2+ii
                     eqcmo2(ii,j)=exc_QCMO(ij,l,k)
                     cqcmo2(ii,j)=coef_QCMO(ij,l,k)
                     rqcmo2(ii,j)=resi_QCMO(ij,l,k)
                  endif
               enddo
            enddo         

            call dgemm('T','N',nval,nval,nval,1.0D0,TX,nval,eqcmo2,
     &                 nval,0.0D0,exc_LMO(:,:,k,l),nval)
            call dgemm('T','N',nval,nval,nval,1.0D0,TX,nval,cqcmo2,
     &                 nval,0.0D0,coef_LMO(:,:,k,l),nval)
            call dgemm('T','N',nval,nval,nval,1.0D0,TX,nval,rqcmo2,
     &                 nval,0.0D0,resi_LMO(:,:,k,l),nval)
            deallocate(eqcmo2,cqcmo2,rqcmo2)
         enddo
      enddo

      allocate(energypair_CIM(nval,nval),energyorb(nval))
      energypair_CIM=0.0D0; energyorb=0.0D0
      do i=1,nval
         do j=1,nval
            do k=1,idimen
               do l=1,idimen
                  xkl=exc_LMO(i,j,k,l)+resi_LMO(i,j,k,l)
                  xlk=exc_LMO(i,j,l,k)+resi_LMO(i,j,l,k)
                  energypair_CIM(i,j)=energypair_CIM(i,j)+
     &                                (2.0D0*xlk-xkl)*coef_LMO(i,j,l,k)
               enddo
            enddo
         enddo
      enddo

      write(6,*) "Pair correlation energies between LMOs:"
      call NJ_prtcol(6,nval,energypair_CIM,1,nval,'f14.10')
      write(6,*)

      do i=1,nval
         do j=1,nval
            energyorb(i)=energyorb(i)+energypair_CIM(i,j)
         enddo
      enddo

      Ecluster=DSUM(nval,energyorb,1)
      Ecen=PSUM(nval,energyorb,NCEN)
      write(6,*) "Orbital energies for LMOs:"
      write(6,'(5F20.9)') energyorb
      write(6,'(A,F20.9)') ' E(CORR-CENTRAL)=',Ecen
      write(6,'(A,F20.9)') ' E(CORR-SUBSYS) =',Ecluster

      if (do_mp4) then 
         call matzero('work1')
         do i=1,nval
            do j=1,i
               call CoefWrite(i,j,iwork1)
            enddo
         enddo
         call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,ccsd,vorb,
     &                        nmo,ivalpair,ivirpair,ifp,nsym,nval,nvirt,
     &                        ij_unique,isympairs,l_restart,ndisk_dump,
     &                        nbf)
         call CCDiis(i,       j,          'reset',      ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
         goto 222
      endif
      if (mp4.and..not.do_mp4) then
        mp3=.false.
        do_mp4=.true.
        goto 666
      endif
      call flush(6)
c
      call dynamic_retmark ! Mark after initR
      call matremark
      call retmark
c
      if ((cc.and.singles).and.l_triples) then
c     call mas_trip(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),res,
c    *             ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *             singles.and.cc.and..not.ccsd,af,nslv)
c     write(6,'(A,F20.15)')'MP4 triples correction: ',res
      call elapsec(xt0)
      time_keep=xt0
      if (mygid.ne.0) STOP 'mygid.ne.0 ??!!'
      call sort_amplit(nval,idimen,ndisk_a,af,mygid,
     *                 ichar_count,bl(ivtable),bl(im_table),bl(ijtable),
     *                 bl(ivrevtable),amp_ratio)
      call elapsec(xt1)
      if (iprnt.gt.2) then
      write(6,'(A,F10.2)')'Amplitudes sort time:       ',xt1-xt0
      endif
      call elapsec(xt0)
ctmp? call fafClosem(ndisktt,0,iresult)
      call master_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                 ndisktt,
     *                 ncf,lbintt,thresh,byt8,nmo,vorb,af,mygid,
     *                 nslv,ext3_ratio)
      call fafClosem(ndisktrtt,0,iresult)
      if (iprnt.gt.2) write(6,*)'RATIO: ',ext3_ratio
ccs      call sleep(10)
      call make_3ext_pairs_mas(idimen,nval,af,ndisk_ie,ndisk_ie1,
     *                         nslv)
      call elapsec(xt1)
      if (iprnt.gt.2)
     *write(6,'(A,F10.2)')'3ext integrals sort time:   ',xt1-xt0
      call elapsec(xt0)
      call sort_3int(nval,idimen,ndisk_ii,af,mygid)
      call elapsec(xt1)
      if (iprnt.gt.2)
     *write(6,'(A,F10.2)')'3int integrals sort time:   ',xt1-xt0
      call dynamic_show_free(mem)
      mem=mem-nval*nval*nval-3*nval*nval
      call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                     nslv,npass,isize)
      call elapsec(xt0)
      call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *               ncf,lbinx,thresh,byt8,nmo,vorb,af,mygid,
     *               isize,npass)
      call elapsec(xt1)
      if (iprnt.gt.2) then
      write(6,'(A,F10.2)')'Kext integrals sort time:   ',xt1-xt0
      write(6,'(A,F10.2)')'Total sort time:            ',xt1-time_keep
      endif
      call flush(6)
      call elapsec(xt0)
      call dynamic_getmem(20,itimes)
      call new_tri_mas(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,ndisk_ix,
     *                       res,bl(iepsi+ncore),bl(itimes),20,ccsd,
     *                       qcisd,af,nslv,ext3_ratio,amp_ratio,npass,
     *                       isize,iprnt,energys,energyd,itrstart)
      call elapsec(xt1)
      write(6,'(A,F10.1)') 'Triples elapsed time: ',(xt1-xt0)/6d1
      if (iprnt.gt.1) then
      write(6,'(A,F20.15,A,F10.2)')
     *  'Triples energy correction: ',res,',   elapsed time: ',xt1-xt0
      write(6,'(A,2F25.15)') "Energies, s & d: ", energys, energyd
      write(6,*)'W zero & sort:           ',bl(itimes+0)
      write(6,*)'W build (total):         ',bl(itimes+1)
      write(6,*)'Mult. over virt space:   ',bl(itimes+2)
      write(6,*)'Mult. over occ. space:   ',bl(itimes+3)
      write(6,*)'Reading of ampl. & int:  ',bl(itimes+4)
      write(6,*)'Amplitudes:              ',bl(itimes+5)
      write(6,*)'3int integrals:          ',bl(itimes+6)
      write(6,*)'3ext integrals:          ',bl(itimes+7)
      write(6,*)'Kext integrals:          ',bl(itimes+9)
      write(6,*)'Energy:                  ',bl(itimes+8)
      write(6,*)'Wabc relocate:           ',bl(itimes+10)
      write(6,*)'Inte relocate:           ',bl(itimes+11)
      write(6,*)'3int relocate:           ',bl(itimes+12)
      write(6,*)'Total virt part:         ',bl(itimes+13)
      write(6,*)'Tables reading:          ',bl(itimes+14)
      call flush(6)
      endif
      endif
c---------------------------------------------------------------------
c
c Finish
 999  continue
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
     *      //   ' * * * * * * * * * * * *'
      write(6,*) '                                Final results:    '
      write(6,*) '                                                     '
        write(6,'(A34,F22.9)') 'SCF energy:                  ',erhf
      if (mp4) then
        write(6,'(A34,F22.9)') 'MP2 correlation energy:      ',xmp2
        write(6,'(A34,F22.9)') 'MP3 correlation energy:      ',xmp3
        write(6,'(A34,F22.9)') 'MP4 singles energy:          ',xmp4_s
        write(6,'(A34,F22.9)') 'MP4 doubles energy:          ',xmp4_d
        wvfnc = 'MP4SDQ'
        if (l_triples) then
          write(6,'(A34,F22.9)') 'MP4 triples energy:          ',res
          wvfnc = 'MP4SDTQ'
        else
          res=0d0
        endif
        write(6,'(A34,F22.9)') 'MP4 quadruples energy:       ',xmp4_q
        write(6,'(A34,F22.9)') 'Total MP4 correlation energy:',
     *     xmp2+xmp3+xmp4_s+xmp4_d+xmp4_q+res
        etot = erhf+xmp2+xmp3+xmp4_s+xmp4_d+xmp4_q+res
      else if (mp3) then
        write(6,'(A34,F22.9)') 'MP2 correlation energy:      ',xmp2
        write(6,'(A34,F22.9)') 'MP3 correlation energy:      ',xmp3
        write(6,'(A34,F22.9)') 'Total MP3 correlation energy:',
     *     xmp2+xmp3
        wvfnc = 'MP3'
        etot = erhf+xmp2+xmp3
      else
       write(6,'(A34,F22.9)')'MP2 correlation energy:      ',xmp2
       write(6,'(A34,F22.9)')'MP2 total energy:            ',xmp2+erhf
       write(6,'(A34,F22.9)')'SCS-MP2 correlation energy:  ',scsmp2
       write(6,'(A34,F22.9)')'SCS-MP2 total energy:        ',scsmp2+erhf
C       write(6,'(A34,F22.9)')'E(CORR-SUBSYS) =             ',xmp2
C       write(6,'(A34,F22.9)')'E(CORR-CENTRAL)=             ',ecenmp2
        if(method(1:3).EQ.'MP2') GO TO 95     ! JB July 2010
       write(6,'(5X,A29,F22.9)')
     *     method(1:lmet)//' correlation energy:'//spc,energyr
        if (l_triples) then
        write(6,'(5X,A29,F22.9)')
     *      method(1:lmet)//' triples correction:'//spc,res
        else
          res=0d0
        endif
        if(ifound(2).eq.0) then
        write(6,'(5X,A29,F22.9)')
     *      'Total '//method(1:lmet)//' energy:'//spc,energyr+erhf+res
        wvfnc = method(1:lmet)
        else
        write(6,'(5X,A29,F22.9)')
     *    'Total '//method(1:lmet)//'(T) energy:'//spc,energyr+erhf+res
        wvfnc = method(1:lmet)//'(T)'
        endif
        etot = energyr+erhf+res
      endif
 95   CONTINUE
      write(6,*) '                            End of CCSD              '
      call f_lush(6)
c 
c     call memory_status('end of CCSD')
      call para_next(-1)
c
c -- write energy to control file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,idum,etot,chopv)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c*******************************************************************
      call fafTerminate(iresult)
      call allclose_and_delete(.true.)
c
      call dynamic_retmark
      call matremark
      call retmark
c---------------------------------------------------------------------
      end subroutine para_CIM_CC


C ********************************************************
C * Subroutine for slaves in parallel CIM-CC calculation *
C * Modify the source code from PQS                      *
C * 7/5/2016_NZG @UARK                                   *
C ********************************************************
      subroutine do_cimcc
      use memory
      use kinds
      use newpara
      implicit none
      integer iMasID,mygid,mytid
c     
      logical nodisk,gauss_seidel,loca,omit,mp2,mp3,cep0,cep2,smal
      logical small,cc,norecalc,byt8,log_diis,singles,ccsd,vorb,nofr
      integer nmo,nvirt,nfirst,nlast,nval,idimen,npairs,icano,iepsi
      integer maxsh,na,ncf,ncs,nbfx,nsh,ibas,inuc,ictr,nsym,ngener
      integer iadr,nupair,nupair1,ifp,ifp1,mataddr,i
      integer ipass,IsWork,lbinc, lrecc, ndiskc, irecadrc
      integer              lbinx, lrecx, ndiskx, irecadrx
      integer              lbine, lrece, ndiske, irecadre
      integer              lbintx,lrectx,ndisktx,irecadrtx
      integer              lbintc,lrectc,ndisktc,irecadrtc
      integer              lbintt,lrectt,ndisktt,irecadrtt
      integer iterate,ndiskYZ,ilist,ising_res,iX,ndiskdiisr
      integer ndiskdiisc,ndiskalpha,length_diis
      logical af,l_triples,qcisd,restart
      real*8  thresh,xmaxdisk,core,shift,shiftini,erhf,tot_xmax
      character*20 chtid
      character*3 chgid,chite
      character*256 hostname
      character*2 int_kind
      integer islaves,ndisktrc,ndisktrx,ndisktre,ndisktrtx,ndisktrtc
      integer ndisktrtt
      integer mygidcp
      common /extrone_tmp/ af,islaves,ndisktrc,ndisktrx,ndisktre,
     *                     ndisktrtx,ndisktrtc,ndisktrtt,mygidcp,nbf
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      logical equal,fileopen
      integer iresAO,iexchAO,iwork1
      integer iwork2,iwork3,iresult,iwork,j,ij,icoeAO,iresult_iter,ibeta
      integer ifockAO,ifockMO,ioverlap,icorefock,icoeffAO,kijklndisk
      integer ndiskg41a,inum,iii,istart,istop,jstart,jstop,ndiskQ
      integer max_dyn_mem,iteration,iad,lcore,memavail,int_array
      real*8 epair,epairr,sqij,xmaxele,xepair,xepairr,xnormij
      real*8 s_ext_matstop,e_ext_matstop,e_ext_mat_tot1,e_ext_matstart
      real*8 eqstart,eqstop,sqstart,sqstop,s_ext_mat_tot1,eq_tot1
      real*8 sq_tot1,egstart,sgstart,egstop,sgstop,eg_tot1,sg_tot1
      real*8 sco_tot1,eco_tot1,ecostop,scostop,scostart,ecostart
      real*8 sen_tot1,een_tot1,t0,t1,eenstop,senstop,CPU_doubles_singles
      real*8 eenstart,senstart,xiterenergy
      real*8 s_ext_matstart,elaps(8),slave_cpu,slave_cpu0,slave_cpu1
      real*8 slave_tot_cpu0,slave_tot_cpu1,slave_tot_cpu
      character*256 home_dir,scrfile
      logical mp4,do_mp4
      integer len,istatus,ndisk_mp3,ndiskmp2,iamount
      real*8 strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
c TEST
      real*8 zeroed,zeroed_max,ext3_ratio
      common /zero_in_EEO/ zeroed,zeroed_max
      integer :: info, ibuffer, nstrong
      integer :: isympairs,ipairimages,ivalpair,ivirpair,ij_unique
      integer iad2,iad3,iinvcano,diis,itmp_mat1,itmp_mat2,icode
      integer ncore,iiR,iiS,iotable,ivtable,ivrevtable,im_table,ijtable
      integer ichar_count,iorevtable,imastgid,itrala,nbf
      real*8 det,sepair,sepairr,sxepair,sxepairr
      logical docansym
      real*8,allocatable::exc_QCMO(:,:,:,:),exc_LMO(:,:,:,:)
      real*8,allocatable::coef_QCMO(:,:,:,:),coef_LMO(:,:,:,:)
      real*8,allocatable::resi_QCMO(:,:,:,:),resi_LMO(:,:,:,:)
      real*8,allocatable::tmpexc(:,:,:,:),tmpcoef(:,:,:,:)
      real*8,allocatable::tmpresi(:,:,:,:)
c
      zeroed=0.0d0
      zeroed_max=0.0d0
      info = 0
      ibuffer = 0
      ichar_count = 1
c
      call fafinit(iresult)
      call reset_counters
      if (iresult.eq.1) return ! I have been an I/O daemon
      if (iresult.ne.0.and.iresult.ne.1)
     * STOP 'fafinit error'
      iMasID=MASTER_ID
      mygid=MY_GID
      mytid=MY_ID
      call secund(slave_tot_cpu0)
      call setival('ccprint',0)
      af=.true.
      do_mp4=.false.
      mygidcp=mygid
c     call getenv('HOME',home_dir)
c     call rmblan(home_dir,80,len) !returns len of string without spaces,
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
      write (chtid,'(I10.10)') mytid
      open(99,FILE=scrfile(1:len)//'.debug_slaves'//chtid)
      write(chgid,'(I3.3)') mygid
c
c Memory markers:
      call matreset
      call mmark
      call matmark
c
      call fafCreates(itrala) ! peer of masters open_mp3_resid
c     Get basic information about geometry etc.
c     call para_recv_bcastpack(TxScfInit)
      call para_get_info
c
c Extract data from para_get_info
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('nbf',nbfx)
      call getival('nsh',nsh)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('ictr',ictr)
      call getival('nsym',nsym)
      if (nsym.gt.0) then
         call getival('ngener',ngener)
         call getival('nsyo',iadr)
         call getival('SymNuPr',nupair)
         call getival('SymNuPr1',nupair1)
         call getival('SymFunPr',ifp)
         call getival('SymShPr',ifp1)
      endif
      call setival('iout',99)
      call setival('icond',99)
c
c Receive symmetry info:
      
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(nbf,1)
      call para_unpack_int(max_dyn_mem,1)
      call dynamic_init(max_dyn_mem,bl(1))
      call para_unpack_int(npairs,1)
      call dynamic_getmem(npairs,isympairs)
      call dynamic_getmem(npairs+1,ipairimages)
      call para_unpack_int(bl(isympairs),npairs)
      call para_unpack_int(bl(ipairimages),npairs+1)
      call para_unpack_int(nval,1)
      call para_unpack_int(nvirt,1)
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      call para_unpack_int(bl(ivalpair),7*nval)
      call para_unpack_int(bl(ivirpair),7*nvirt)
      call para_unpack_int(ij_unique,1)
      call para_unpack_int(loca,1)
      call para_unpack_int(docansym,1)
      if (.not.loca) then
         call dynamic_getmem(8*8,im_table)
         call dynamic_getmem(8*(nval+1),iotable)
         call dynamic_getmem(8*(nvirt+1),ivtable)
         call dynamic_getmem(nval,iorevtable)
         call dynamic_getmem(nvirt,ivrevtable)
         call dynamic_getmem((nval*nval+1)*8,ijtable)
         call para_unpack_int(bl(iotable),8*(nval+1))
         call para_unpack_int(bl(ivtable),8*(nvirt+1))
         call para_unpack_int(bl(iorevtable),nval)
         call para_unpack_int(bl(ivrevtable),nvirt)
         call para_unpack_int(bl(im_table),8*8)
         call para_unpack_int(bl(ijtable),(nval*nval+1)*8)
         call para_unpack_int(ichar_count,1)
      endif
c CCSD data:
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(nodisk,1)
      call para_unpack_int(gauss_seidel,1)
      call para_unpack_int(omit,1)
      call para_unpack_int(mp2,1)
      call para_unpack_int(mp3,1)
      call para_unpack_int(mp4,1)
      call para_unpack_int(cep0,1)
      call para_unpack_int(cep2,1)
      call para_unpack_int(smal,1)
      call para_unpack_int(small,1)
      call para_unpack_int(cc,1)
      call para_unpack_int(norecalc,1)
      call para_unpack_int(byt8,1)
      call para_unpack_int(log_diis,1)
      call para_unpack_int(singles,1)
      call para_unpack_int(ccsd,1)
      call para_unpack_int(qcisd,1)
      call para_unpack_int(vorb,1)
      call para_unpack_int(nofr,1)
      call para_unpack_int(cache,1)
      call para_unpack_int(integerT,1)
      call para_unpack_int(l_triples,1)
      call para_unpack_int(restart,1)
c int:
      call para_unpack_int(nmo,1)
      call para_unpack_int(nvirt,1)
      call para_unpack_int(nfirst,1)
      call para_unpack_int(nlast,1)
      call para_unpack_int(nval,1)
      call para_unpack_int(idimen,1)
      call para_unpack_int(npairs,1)
      call para_unpack_int(length_diis,1)
      call para_unpack_int(ndisk_mp3,1)
      call para_unpack_int(ncore,1)
      call dynamic_mmark
c real
      call para_unpack_real(erhf,1)
      call para_unpack_real(thresh,1)
      call para_unpack_real(xmaxdisk,1)
      call para_unpack_real(core,1)
      call para_unpack_real(shift,1)
      call para_unpack_real(shiftini,1)
      call dynamic_matdef('cano','q',ncf,ncf)
      call dynamic_matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      call para_unpack_real(bl(icano),ncf*ncf)
      call para_unpack_real(bl(iepsi),ncf)
      islaves=nslv
      call setival('nslv',nslv)
      call matsub('genvirt','cano',nmo+1,nbf)
      call dynamic_lock(bl(icano),i)
      call dynamic_lock(bl(iepsi),i)
c
c Initialize symmetry commons
      call symm_memory_allocator(ncf)
      call pair_s_initializer(vorb,     nval,     nsym, ifp, ifp1,
     *                        ivalpair, ivirpair, loca,iorevtable,
     *                        ivrevtable,im_table)
c
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrc)
      call izeroit(bl(irecadrc),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadrx)
      call izeroit(bl(irecadrx),npairs*nslv)    !zero out the pair records
      call dynamic_getmem((npairs/intsize+1)*nslv,irecadre)
      call izeroit(bl(irecadre),npairs*nslv)         !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtx)
      call izeroit(bl(irecadrtx),nval*nval*nslv)     !zero out the pair records
      call dynamic_getmem((nval*nval/intsize+1)*nslv,irecadrtc)
      call izeroit(bl(irecadrtc),nval*nval*nslv)   !zero out the pair records
      call dynamic_getmem((nval*(nbf-nmo)/intsize+1)*nslv,irecadrtt)
      call izeroit(bl(irecadrtt),nval*(nbf-nmo)*nslv)!zero out the pair records
c
      call twoelinit ! initalize callable two-electron program
      call getint(npairs,ilist)
      call filllist(bl(ilist),npairs)
c
      int_kind='c'
      ipass=0
      if (.not.mp2.and..not.(restart.and.l_triples)) then
         call slave_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,core,
     &                  xmaxdisk,ndiskc,bl(irecadrc),lrecc,lbinc,
     &                  int_kind,nodisk,byt8,small,vorb,iMasID,mygid,
     &                  mytid,ndisktrc,tot_xmax,bl(ilist),nstrong,
     &                  bl(isympairs),ij_unique,nbf)
      endif
      call izeroit(bl(irecadrx),npairs*nslv)
      int_kind='x'
      ipass=0
      call slave_Gen(ncs,     ncf,          ictr,    nval,   nmo,
     *               nfirst,  nlast,        thresh,  core,   xmaxdisk,
     *               ndiskx,  bl(irecadrx), lrecx,   lbinx,  int_kind,
     *               nodisk,  byt8,         small,   vorb,   iMasID,
     *               mygid,   mytid,        ndisktrx,tot_xmax,bl(ilist),
     *               nstrong,bl(isympairs),ij_unique,nbf)
      if (l_triples) then
         call izeroit(bl(irecadrtt),nval*(nbf-nmo)*nslv)
         int_kind='tt'
         ipass=0
         call slave_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,core,
     &                  xmaxdisk,ndisktt,bl(irecadrtt),lrectt,lbintt,
     &                  int_kind,nodisk,byt8,small,vorb,iMasID,mygid,
     &                  mytid,ndisktrtt,tot_xmax,bl(ilist),nstrong,
     &                  bl(isympairs),ij_unique,nbf)
      endif
c
c  Initialization:
c
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
      call dynamic_matdef('fockAO','q',idimen,idimen)
      ifockAO=mataddr('fockAO')
      call dynamic_getmem(nbf*nbf,ifockMO)
      call dynamic_matdef('overlap','q',ncf,ncf)
      ioverlap=mataddr('overlap')
      call dynamic_matdef('corefock','q',idimen,idimen)
      icorefock=mataddr('corefock')
c
      if (.not.mp2) then
         call Kijk_Vec_Slave(bl(irecadrx),npairs,ndiskx,lbinx,thresh,
     *                          byt8,kijklndisk,ncf,nmo,nval,
     *                          vorb,af,mygid,nfirst,nlast,
     *                          nslv,iMasID,mytid,nbf)
      endif
      call sync_barrier
c
      call para_recv_bcastreal(bl(ifockAO),idimen*idimen,
     *                               TxCCSDInit)
      call para_recv_bcastreal(bl(ifockMO),nbf*nbf,TxCCSDInit)
      call para_recv_bcastreal(bl(ioverlap),ncf*ncf,TxCCSDInit)
      call para_recv_bcastreal(bl(icorefock),idimen*idimen,
     *                               TxCCSDInit)
      call dynamic_lock(bl(icorefock),iamount)
      call dynamic_lock(bl(ifockAO),iamount)
      call dynamic_lock(bl(ifockMO),iamount)
c
      if (.not.mp2) then
         call fafCreates(itrala)
         call sync_barrier

         call KijklInit_Slave(bl(irecadrx),npairs,ndiskx,lbinx,ncf,
     *                        thresh,nfirst,nlast,bl(ifockMO),byt8,
     *                        nmo,vorb,af,nslv,iMasID,mytid)
         call sync_barrier
      endif
c
      call Fock_Vector_Init(ncf,nfirst,nlast,ifockAO,ifockMO,nmo,vorb,
     *                      nbf)
      call initsingles_slave(nval,ncf,singles.or.mp4,nmo,vorb,iMasID,
     *                       nbf)
      call initS(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call initR(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call mmark
      call dynamic_mmark
      call matmark
c
c Space for one K_ij matrix
      call dynamic_matdef('resao','q',idimen,idimen)
      iresAO=mataddr('resao')
      call dynamic_matdef('exchao','q',idimen,idimen)
      iexchAO=mataddr('exchao')
      call dynamic_matdef('invcano','q',ncf,ncf)
      iinvcano=mataddr('invcano')
      call matzero('invcano')
      call matcopy('cano','invcano')
      call matmmul2('cano','overlap','invcano','t','n','n')
c
      call dynamic_matdef('rinvcano','q',ncf,ncf)
      call matcopy('invcano','rinvcano') !store r(eal) inverted cano for singles
      call dynamic_matdef('canotran','q',ncf,ncf)
      call matcopy('cano','canotran')
      call matpose('canotran')
c  This is necessary because of transf: T_MO=U-1 T_AO U-1(T)
c I use later matsimtr, which makes: U(T)A U
      call matpose('invcano')
c locking:
      i=mataddr('invcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('rinvcano')
      call dynamic_lock(bl(i),j)
      i=mataddr('canotran')
      call dynamic_lock(bl(i),j)
c
      call dynamic_matdef('work1','q',idimen,idimen)
      iwork1=mataddr('work1')
      call dynamic_matdef('work2','q',idimen,idimen)
      iwork2=mataddr('work2')
      call dynamic_matdef('work3','q',idimen,idimen)
      iwork3=mataddr('work3')
      call dynamic_matdef('beta','q',nval,nval)
      ibeta=mataddr('beta')
      if (log_diis) then
         call fafcreates(ndiskdiisr)
         call fafcreates(ndiskdiisc)
      endif
      call sync_barrier
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(ndiskdiisr,1)
      call para_unpack_int(ndiskdiisc,1)
      if (.not.log_diis) then
         call CCDiis(i,       j,          'nodiis',     ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
      else
         i=length_diis
         call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
      endif
c
c  Calculation:
c
      if (restart) goto 126
      iresult=0
      do
         call para_initsend
         call para_pack_int(mytid,1)
         call para_pack_int(iresult,1)
         if (iresult.eq.1) then
            call para_pack_real(epair,1)
            call para_pack_real(sepair,1)
            call para_pack_real(sqij,1)
            call para_pack_real(xmaxele,1)
         endif
         call para_send_pack(0,TxCCSDMP2Req)
         xmaxele=0.0d0
c
         call para_recv_pack(imastgid,TxCCSDMP2Job)
         call para_unpack_int(iwork,1)
         if (iwork.eq.0) exit
         call para_unpack_int(i,1)
         call para_unpack_int(j,1)
         call para_unpack_int(ij,1)
c
c calculation, the kernel from master:
c
         equal=.false.
         if (i.eq.j) equal=.true.
         call ExtrOne(i,    j,     bl(irecadrx),  npairs,     ndiskx,
     *                ncf,  lbinx, thresh,        byt8,       'x',
     *                'mo', nmo,   vorb,          bl(iexchAO))
         call CoefRead('tt',i,j,icoeAO)
         call UpdateCoef(i,       j,     icoeAO,  iexchAO,    bl(iepsi),
     *                   ncf,     nmo,   nfirst,  npairs,     shiftini,
     *                   'nodiis',vorb, singles, ndiskdiisr, ndiskdiisc,
     *                   sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
         call CoefWrite(i,j,icoeAO)
         call PairEnergyMP2(icoeAO,iexchAO,iresAO,epairr,epair,
     *                     ncf,   equal,  nmo,   vorb,  sepairr,
     *                     sepair,nbf)
c       print '(A15,I5,5X,F20.10)', "Pair energy: ",ij,epair
         iresult=1
c
      enddo
      if (mp2.and..not.loca) goto 999
      call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                    vorb,  nmo,iMasID,      mygid, ivalpair,
     *                    ivirpair,ifp,nsym,      nval,  nvirt,
     *                    ij_unique,isympairs,nbf)
c Local MP2:
 222  continue
      do             ! iterative loop
         call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
         iresult=0
         call fafCreates(ndiskmp2)
         call sync_barrier
         call para_recv_bcast(ndiskmp2,MP2file)
         call Generate_MP2_G_SLA(nval,ncf,ndiskmp2,nmo,nfirst,vorb,
     &                           .true.,bl(ifockMO),nslv,iMasID,mytid,
     &                           nbf)
         call sync_barrier
         do           ! ij loop
            call para_initsend
            call para_pack_int(mytid,1)
            call para_pack_int(iresult,1)
            if (iresult.eq.1) then
               call para_pack_int(i,1)
               call para_pack_int(j,1)
               call para_pack_int(ij,1)
               call para_pack_real(epair,1)
               call para_pack_real(epairr,1)
               call para_pack_real(xepair,1)
               call para_pack_real(sepair,1)
               call para_pack_real(sepairr,1)
               call para_pack_real(sxepair,1)
               call para_pack_real(sqij,1)
               call para_pack_real(xmaxele,1)
            endif
            call para_send_pack(0,TxCCSDMP2Req)
            xmaxele=0.0d0
c
            call para_recv_pack(imastgid,TxCCSDMP2Job)
            call para_unpack_int(iwork,1)
            if (iwork.eq.0) exit
            call para_unpack_int(i,1)
            call para_unpack_int(j,1)
            call para_unpack_int(ij,1)
c
c calculation, the kernel from master:
c
            equal=.false.
            if (i.eq.j) equal=.true.
            if (do_mp4) then
               call fafread(ndisk_mp3,bl(iexchAO),8,idimen*idimen,1,ij,
     *                      istatus)
               if (istatus.lt.8*idimen*idimen) then
                  stop "error in slave"
                  call flush(6)
               endif
            else
               call ExtrOne(i,j,bl(irecadrx),npairs,ndiskx,ncf,lbinx,
     &                      thresh,byt8,'x','mo',nmo,vorb,bl(iexchAO))
            endif
            call matcopy('exchao','resao')
            call CoefRead('tt',i,j,icoeffAO)
            call matconn('coefAO','q',idimen,idimen,icoeffAO)
            if (vorb) then
               call matmmul2('fockAO','coefAO','resao','n','n','a')
               call matmmul2('coefAO','fockAO','resao','n','n','a')
            else
               call matmmul2('fockAO','coefAO','work1','n','n','n')
               call matmmul2('work1','overlap','resao','n','n','a')
               call matmmul2('coefAO','fockAO','work1','n','n','n')
               call matmmul2('overlap','work1','resao','n','n','a')
            endif
            call matdisc('coefAO')
            call Read_MP2(i,j,nmo,ncf,nval,ndiskmp2,vorb,.true.,iwork1,
     *                    nbf)
            if (vorb) then
               call matadd('work1','resao')
            else
               call matsimtr('work1','overlap','work2')
               call matadd('work2','resao')
            endif
            call CoefRead('tt',i,j,icoeffAO)
            call PairEnergyMP2(icoeffAO,iexchAO,iresAO,epairr,epair,
     *                         ncf,     equal,  nmo,   vorb,  sepairr,
     *                         sepair,nbf)
            call UpdateCoef(i,j,icoeffAO,iresAO,bl(iepsi),ncf,nmo,
     &                      nfirst,npairs,shift,'diis',vorb,singles,
     &                      ndiskdiisr,ndiskdiisc,sqij,xmaxele,
     &                      bl(ilist),bl(ipairimages),bl(iorevtable),
     &                      bl(ivrevtable),bl(im_table),docansym,nbf)
            call CoefWrite(i,j,icoeffAO)
            call PairEnergyMP2(icoeffAO,iexchAO,iresAO,xepairr,xepair,
     *                         ncf,     equal,  nmo,   vorb,   sxepairr,
     *                         sxepair, nbf)
c
            iresult=1
         enddo        ! ij loop
         call fafCloses
         call sync_barrier
         call sync_barrier
         call para_recv_bcast(diis,TxCCSDInit)
         if (diis.eq.1.and.log_diis) then
            call dynamic_getmem(nvirt*nvirt,itmp_mat1)
            call dynamic_getmem(nvirt*nvirt,itmp_mat2)
            call sla_diis_matrix(af,singles,npairs,bl(ilist),
     &                           bl(ipairimages),ncf,nmo,bl(itmp_mat1),
     &                           bl(itmp_mat2),iMasID,mygid,mytid,nbf)
            call dynamic_retmem(2)
 345        continue
            call sync_barrier
            call para_recv_bcast(icode,TxDIISFlow)
            if (icode.eq.2) goto 346
            if (icode.eq.1) goto 345
            if (icode.ne.0) STOP 'Cannot happen, on slave DIIS'
            call sla_new_coeff(singles,af,idimen,nval,bl(ilist),
     *                        bl(ipairimages),npairs,iMasID,mygid,mytid)

 346        continue
         endif
         call sync_barrier
         call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                       vorb,  nmo,iMasID,      mygid, ivalpair,
     *                       ivirpair,ifp,nsym,      nval,  nvirt,
     *                       ij_unique,isympairs,nbf)
         call para_barrier
         call para_recv_bcast(iterate,TxCCSDIter)
         call para_barrier
         if (iterate.eq.0) exit
      enddo          ! iterative loop
      if (mp2.or.do_mp4) then
         goto 999
      endif
 126  continue
c
      call para_recv_bcastpack(TxCCSDInit)
      call para_unpack_int(bl(ilist),npairs)
      call para_unpack_int(nstrong,1)
      if (info.ne.0) then
         stop "info ne 0"
         call flush(6)
      endif
      call dynamic_matdef('WernerX','q',idimen,idimen)
      iX=mataddr('WernerX')
      call dynamic_matdef('sing_res','r',idimen,1)
      ising_res=mataddr('sing_res')
      if ((mp4.or.restart).and.l_triples) then
         call fafCreates(itrala)
         call fafCreates(itrala)
         call fafCloses
         call slave_3ext(nval,idimen,bl(irecadrtt),npairs,ndisktt,ncf,
     &                   lbintt,thresh,byt8,nmo,vorb,af,iMasID,mygid,
     &                   mytid,ichar_count,bl(ivtable),bl(im_table),
     &                   bl(ivrevtable),bl(iotable),ext3_ratio)
         call fafCloses
         call make_3ext_pairs_sla(idimen,nval,af,ichar_count,
     &                            bl(ivtable),bl(im_table),bl(iotable),
     &                            bl(ivrevtable),iMasID,mygid,mytid)
         call fafCreates(itrala)
         call fafCreates(itrala)
         call fafCloses
         call fafCreates(itrala)
         call fafCreates(itrala)
         call fafCloses
         call new_tri_sla(nval,idimen,bl(iepsi+ncore),ccsd,qcisd,af,
     &                    iMasID,mygid,mytid,bl(iotable),bl(ivtable),
     &                    bl(ivrevtable),bl(im_table),bl(ijtable),
     &                    ichar_count)
         if (restart) goto 999
      endif
c
c The main CCSD loop:
c
      iteration=1
 666  continue   ! for MP4 only
      do
         call secund(slave_cpu0)
         iteration=iteration+1
         call para_barrier
         call para_recv_bcastpack(TxCCSDIter1)
         call para_unpack_real(tot_xmax,1)
         call para_unpack_int(iterate,1)
         call para_barrier

         if (iteration.eq.2) tot_xmax=0.1d0
c
         call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *               npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *               singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *               bl(ipairimages),nbf)
c
         int_kind='e'
         ipass=0
         if (.not.do_mp4) then
            call SLA_EEO_INT(ncs,    ncf,     bl(ictr),   nval,  nmo,
     *                       nfirst, nlast,   thresh,
     *                       vorb,   ndisktre, npairs,
     *                       iMasID, mygid,   mytid,  nslv, nbf)
            ndiske=-1
            lrece =-1
            lbine =-1
         endif
         if (ccsd) then
            int_kind='tx'
            ipass=0
            call slave_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,
     &                     core,xmaxdisk,ndisktx,bl(irecadrtx),lrectx,
     &                     lbintx,int_kind,nodisk,byt8,small,vorb,
     &                     iMasID,mygid,mytid,ndisktrtx,
     &                     tot_xmax,bl(ilist),nstrong,bl(isympairs),
     &                     ij_unique,nbf)
            int_kind='tc'
            ipass=0
            call slave_Gen(ncs,ncf,ictr,nval,nmo,nfirst,nlast,thresh,
     &                     core,xmaxdisk,ndisktc,bl(irecadrtc),lrectc,
     &                     lbintc,int_kind,nodisk,byt8,small,vorb,
     &                     iMasID,mygid,mytid,ndisktrtc,tot_xmax,
     &                     bl(ilist),nstrong,bl(isympairs),
     &                     ij_unique,nbf)
         endif
c
         do i=1,8
            elaps(i)=0.0d0
         enddo
         call fafCreates(ndiskalpha)
         call sync_barrier
         call para_recv_bcast(ndiskalpha,TxCCSDInit)
         if (cc.or.ccsd.or.do_mp4) then
            call CCalphaslave(ncf,nval,bl(irecadrx),npairs,ndiskx,lbinx,
     *                        thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                        ndiskalpha,af,iMasID,mytid,nbf)
            call sync_barrier
            call para_reduce(elaps,5,TxCCSDRedu)
         endif
         call sync_barrier
         call zeroit(bl(ibeta),nval*nval)
         call para_recv_bcastreal(bl(ibeta),nval*nval,TxCCSDInit)
         call B41a_G_Slave(ncf,nval,norecalc,cc,ccsd,nfirst,bl(ifockMO),
     *                     nmo,vorb,ndiskG41a,.true.,ndiskalpha,
     &                     bl(ibeta),iMasID,mytid,nbf)
         call fafCloses
         call sync_barrier
c A matrix
         if (cc.or.ccsd.or.do_mp4) then
            call CCAOnceslav(bl(irecadrx),npairs,nval,ndiskx,ncf,lbinx,
     &                       thresh,norecalc,byt8,ioverlap,nmo,vorb,
     &                       nslv,iMasID,mygid,mytid,nbf)
            call sync_barrier
         endif
c
         if (singles.or.ccsd.or.do_mp4) then
            call Tl_gen_slav(ncf,nval,nmo,vorb,nslv,iMasID,mygid,mytid,
     &                       nbf)
            call sync_barrier
            call EEO_vector_slave(nval,bl(irecadre),ndiske,lbine,
     &                            bl(irecadrx),ndiskx,lbinx,
     &                            bl(irecadrc),ndiskc,lbinc,ncf,npairs,
     &                            thresh,nfirst,nlast,byt8,nmo,vorb,
     &                            nslv,iMasID,mytid,bl(ilist),nbf)
            call sync_barrier
            call unlockS
            call unlockR
            call pointS(iiS)
            call para_recv_bcastreal(bl(iiS),idimen*nval,TxCCSDInit)
            if (cc.or.ccsd) then
               call Lt_gen_slav(idimen,nval,bl(irecadrx),npairs,ndiskx,
     *                          ncf,   lbinx,  thresh,      byt8,  nmo,
     *                          vorb,  nslv,  iMasID,      mygid, mytid)
               call sync_barrier
            endif
            call pointRR(iiR)
            call para_recv_bcastreal(bl(iiR),idimen*nval,TxCCSDInit)
            call lockS
            call lockR
         endif
c  X matrix:
         call sync_barrier
         call para_recv_bcastpack(TxCCSDInit)
         call para_unpack_real(bl(iX),idimen*idimen)
         call para_unpack_int(ndiskg41a,1)
         call para_unpack_real(xiterenergy,1)
c
c  YZ matrices:
c
         do i=1,5
            elaps(i)=0.0d0
         enddo
         iresult=0
         call fafCloses()
         call fafCreates(ndiskYZ)
         call sync_barrier
         call para_recv_bcast(ndiskYZ,TxCCSDInit)
         do
            call para_initsend
            call para_pack_int(mytid,1)
            call para_pack_int(iresult,1)
            if (iresult.eq.1) then
               call para_pack_int(istart,1)
               call para_pack_int(istop,1)
               call para_pack_int(jstart,1)
               call para_pack_int(jstop,1)
            endif
            call para_send_pack(0,TxCCSDReq)
c 
            call para_recv_pack(imastgid,TxCCSDJob)
            call para_unpack_int(iwork,1)
            if (iwork.eq.0) exit
            call para_unpack_int(istart,1)
            call para_unpack_int(istop,1)
            call para_unpack_int(jstart,1)
            call para_unpack_int(jstop,1)
c
c calculation, the kernel from master:
c
            call CCYZ(istart,istop,jstart,jstop,ncf,nval,bl(irecadrx),
     *                bl(irecadrc),npairs,ndiskx,ndiskc,lbinx,lbinc,
     &                bl(irecadrtx),ndisktx,lbintx,bl(irecadrtc),
     &                ndisktc,lbintc,thresh,byt8,ioverlap,cc,ccsd,
     &                bl(ilist),nmo,vorb,ndiskYZ,.true.,elaps,do_mp4,
     *                ichar_count,bl(iorevtable),bl(ivtable),
     *                bl(im_table),nbf)
c
            iresult=1
         enddo
c
c  YZ matrices STOP
         call sync_barrier
         call para_reduce(elaps,5,TxCCSDRedu)
         call para_barrier
c  Q parts:
c
         iresult=0
         call fafCloses()
         call fafCreates(ndiskQ)
         call sync_barrier
         call para_recv_bcast(ndiskQ,TxCCSDInit)
         do
            call para_initsend
            call para_pack_int(mytid,1)
            call para_pack_int(iresult,1)
            if (iresult.eq.1) then
               call para_pack_int(istart,1)
               call para_pack_int(istop,1)
               call para_pack_int(jstart,1)
               call para_pack_int(jstop,1)
            endif
            call para_send_pack(0,TxCCSDReqQ)
c
            call para_recv_pack(imastgid,TxCCSDJobQ)
            call para_unpack_int(iwork,1)
            if (iwork.eq.0) exit
            call para_unpack_int(istart,1)
            call para_unpack_int(istop,1)
            call para_unpack_int(jstart,1)
            call para_unpack_int(jstop,1)
c
c calculation, the kernel from master:
            call Qparts(istart,istop,jstart,jstop,ncf,nval,npairs,byt8,
     *                  cc,ccsd,bl(ilist),nmo,vorb,ndiskYZ,af,ndiskQ,
     *                  ichar_count,bl(iorevtable),bl(ivtable),
     &                  bl(im_table),nbf)
c
            iresult=1
         enddo
         call sync_barrier
c
c  YZ matrices STOP
c
c   The main CCSD doubles loop:
c
         iresult_iter=0
         xmaxele=0.0d0
         do           ! ij loop
            call para_initsend
            call para_pack_int(mytid,1)
            call para_pack_int(iresult_iter,1)
            if (iresult_iter.eq.1) then
               call para_pack_int(i,1)
               call para_pack_int(j,1)
               call para_pack_int(ij,1)
               call para_pack_real(epair,1)
               call para_pack_real(epairr,1)
               call para_pack_real(xepair,1)
               call para_pack_real(sqij,1)
               call para_pack_real(xmaxele,1)
               call para_pack_real(xnormij,1)
            endif
            call para_send_pack(0,TxCCSDReq1)
            if (iresult_iter.eq.1) xmaxele=0.0d0
c
            call para_recv_pack(imastgid,TxCCSDJob1)
            call para_unpack_int(iwork,1)
            if (iwork.eq.0) exit
            call para_unpack_int(i,1)
            call para_unpack_int(j,1)
            call para_unpack_int(ij,1)
            if (cep2) then
               call para_unpack_real(xiterenergy,1)
            endif
c
c  Kernel from master:
            equal=.false.
            if (i.eq.j) equal=.true.
            if (int_array(bl(ilist),ij).eq.1) then
c
c  Read exchange operators
               if (.not.do_mp4) then
                  call ExtrOne(i,j,bl(irecadrx),npairs,ndiskx,ncf,lbinx,
     &                        thresh,byt8,'x','mo',nmo,vorb,bl(iexchAO))
c  Read external echange operators
                  call ExtrOne(i,j,bl(irecadre),npairs,ndiske,ncf,lbine,
     *                         thresh,byt8,'e','mo',nmo,vorb,bl(iwork1))
c  Add it do residuum matrix
                  call matcopy('exchao','resao')
                  call matadd('work1','resao')
               else
                  call matzero('resao')
               endif
c
               call secund(s_ext_matstop)
               call elapsec(e_ext_matstop)
c
              e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
              s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c  Qgen, iwork1=Qij, iwork2=Qji(T), probably OK, order of en. OK
               call elapsec(eqstart)
               call secund(sqstart)
c
               call Q_read_build(i,j,iX,ndiskQ,af,vorb,ncf,nmo,nval,
     *                           iwork1,iwork2,nbf)
c
               call elapsec(eqstop)
               call secund(sqstop)
c
               eq_tot1=eq_tot1+eqstop-eqstart
               sq_tot1=sq_tot1+sqstop-sqstart
c   S*Qij and add to residuum
               call elapsec(e_ext_matstart)
               call secund(s_ext_matstart)
c
               if (vorb) then
                  call matadd('work1','resao')
                  call matadd('work2','resao')
               else
                  call matmmult('overlap','work1','work3')
                  call matadd('work3','resao')
                  call matmmult('work2','overlap','work3')
                  call matadd('work3','resao')
               endif
c   Qji(T)*S and add to residuum
c
               call secund(s_ext_matstop)
               call elapsec(e_ext_matstop)
c
              e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
              s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c   Read Gij and Gji(T) with part of CCSD
               call elapsec(egstart)
               call secund(sgstart)
c
               call ReadCalcG_41a(i,j,norecalc,ncf,nmo,vorb,ndiskG41a,
     *                            iwork1,.true.,nbf)
c
               call elapsec(egstop)
               call secund(sgstop)
c
               eg_tot1=eg_tot1+egstop-egstart
               sg_tot1=sg_tot1+sgstop-sgstart
c
               if (cc.or.ccsd.or.do_mp4) then
                  continue
               else
                  call elapsec(e_ext_matstart)
                  call secund(s_ext_matstart)
c
                  call CoefRead('tt',i,j,icoeffAO)
                  call matconn('coeffAO','q',idimen,idimen,icoeffAO)
                  call matadd1('coeffAO',-xiterenergy,'work1')
                  call matdisc('coeffAO')
               endif
c       S*Sum*S
               if (vorb) then
                  call matadd('work1','resao')
               else
                  call matsimtr('work1','overlap','work2')
                  call matadd('work2','resao')
               endif
c
               call secund(s_ext_matstop)
               call elapsec(e_ext_matstop)
c
              e_ext_mat_tot1=e_ext_mat_tot1+e_ext_matstop-e_ext_matstart
              s_ext_mat_tot1=s_ext_mat_tot1+s_ext_matstop-s_ext_matstart
c
               call elapsec(eenstart)
               call secund(senstart)
c
c Singles part of CISD CCSD start
               call secund(t0)
               if (singles.or.ccsd) then
                  call dynamic_matdef('CISD','q',idimen,idimen)
                  iresult=mataddr('CISD')
                  call singles_CISD(i,j,ncf,nval,ioverlap,bl(irecadre),
     &                              npairs,nfirst,nlast,ndiske,lbine,
     &                              thresh,byt8,ccsd,iresult,nmo,vorb,
     &                              nbf)
                  call matadd('CISD','resao')
                  call dynamic_matrem('CISD')
               endif
               call secund(t1)
               CPU_doubles_singles=CPU_doubles_singles+t1-t0
c Singles part of CISD CCSD stop
            else ! if (int_array(bl(ilist),ij).eq.1) then
               call matzero('resao')
               call ExtrOne(i,j,bl(irecadrx),npairs,ndiskx,ncf,lbinx,
     &                      thresh,byt8,'x','mo',nmo,vorb,bl(iexchAO))
            endif
            call CoefRead('tc',i,j,icoeffAO)
            if (mp3.or.do_mp4) then
               call PairEnergy(icoeffAO,iresAO,iresAO,epairr,epair,ncf,
     &                         equal,nmo,vorb,nbf)
               sqij=0.0d0
               if (mp3.and.mp4) call fafWrite(ndisk_mp3,bl(iresAO),8,
     *                                       idimen*idimen,1,ij,istatus)
            else
               call PairEnergy(icoeffAO,iexchAO,iresAO,epairr,epair,ncf,
     &                         equal,nmo,vorb,nbf)
               if (iterate==0) then
                  call para_initsend
                  call para_pack_int(i,1)
                  call para_pack_int(j,1)
                  call para_pack_real(bl(iexchAO),idimen*idimen)
                  call para_pack_real(bl(icoeffAO),idimen*idimen)
                  call para_pack_real(bl(iresAO),idimen*idimen)
                  call para_send_pack(0,TxCIMCCSDInt)
               endif       
            endif
            if (omit.and.int_array(bl(ilist),ij).ne.1) epairr=epair
c
            call elapsec(eenstop)
            call secund(senstop)
c
            een_tot1=een_tot1+eenstop-eenstart
            sen_tot1=sen_tot1+senstop-senstart
c
            call elapsec(ecostart)
            call secund(scostart)
c
            call CoefRead('tt',i,j,icoeffAO)
            if (.not.(mp3.or.do_mp4)) then
               call CINorm(icoeffAO,xnormij,ncf,equal,nmo,vorb,nbf)
               if (.not.omit.or.int_array(bl(ilist),ij).eq.1) then
                  call UpdateCoef(i,j,icoeffAO,iresAO,bl(iepsi),ncf,nmo,
     &                            nfirst,npairs,shift,'diis',vorb,
     &                            singles,ndiskdiisr,ndiskdiisc,sqij,
     &                            xmaxele,bl(ilist),bl(ipairimages),
     *                            bl(iorevtable),bl(ivrevtable),
     &                            bl(im_table),docansym,nbf)
               else
                  sqij=0.0d0
                  xmaxele=0.0d0
               endif
            endif
c
            call elapsec(ecostop)
            call secund(scostop)
c
            eco_tot1=eco_tot1+ecostop-ecostart
            sco_tot1=sco_tot1+scostop-scostart
c
            call CoefWrite(i,j,icoeffAO) ! Still I write "updated"
c          amplitudes in CoefWrite with "omit" switched on, because 
c          of non Gauss-Seidel algorithm
c
            call elapsec(eenstart)
            call secund(senstart)
c
            call CoefRead('ee',i,j,icoeffAO)
            call PairEnergy(icoeffAO,iexchAO,iresAO,xepairr,xepair,ncf,
     &                      equal,nmo,vorb,nbf)
            if (omit.and.int_array(bl(ilist),ij).ne.1) xepairr=xepair
c
            call elapsec(eenstop)
            call secund(senstop)
c
            een_tot1=een_tot1+eenstop-eenstart
            sen_tot1=sen_tot1+senstop-senstart
c
            iresult_iter=1
         enddo ! of the loop over i,j pairs in main CCSD loop
         if (iterate.eq.0) exit
         call sync_barrier
         if (singles.or.ccsd.or.do_mp4) then
            do
               call sync_barrier
               call para_recv_bcast(iwork,TxCCSDInit)
               if (iwork.eq.0) exit
               call Tr_gen_slav(idimen, nval,  ncf,  vorb, nslv,
     *                          iMasID, mygid, mytid)
            enddo
         endif
         call sync_barrier
c
         if (.not.(mp3.or.do_mp4)) then
            call sync_barrier
            call para_recv_bcast(diis,TxCCSDInit)
            if (diis.eq.1.and.log_diis) then
               call dynamic_getmem(nvirt*nvirt,itmp_mat1)
               call dynamic_getmem(nvirt*nvirt,itmp_mat2)
               call sla_diis_matrix(af,singles,npairs,bl(ilist),
     &                              bl(ipairimages),ncf,nmo,
     &                              bl(itmp_mat1),bl(itmp_mat2),
     *                              iMasID,mygid,mytid,nbf)
               call dynamic_retmem(2)
 745           continue
               call sync_barrier
               call para_recv_bcast(icode,TxDIISFlow)
               if (icode.eq.2) goto 746
               if (icode.eq.1) goto 745
               if (icode.ne.0) STOP 'Cannot happen, on slave DIIS'
               call sla_new_coeff(singles,af,idimen,nval,bl(ilist),
     *                       bl(ipairimages),npairs,iMasID,mygid,mytid)

 746           continue
            endif
         endif
c
         call sync_barrier
         call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                       vorb,  nmo,iMasID,      mygid, ivalpair,
     *                       ivirpair,ifp,nsym,      nval,  nvirt,
     *                       ij_unique,isympairs,nbf)
         call initsingles_slave(nval,ncf,singles.or.mp4,nmo,vorb,iMasID,
     *                          nbf)
         call secund(slave_cpu1)
         slave_cpu=slave_cpu1-slave_cpu0
         call sync_barrier
         call para_reduce(slave_cpu,1,TxCCSDRedu)
         call para_barrier
c
         call sync_barrier
      enddo ! of the main ccsd loop

      if (do_mp4) then 
         call CoefInit_slave(npairs,ncf,gauss_seidel,nodisk,ccsd,
     *                       vorb,  nmo,iMasID,      mygid, ivalpair,
     *                       ivirpair,ifp,nsym,      nval,  nvirt,
     *                       ij_unique,isympairs,nbf)
         goto 222
      endif
      if (mp4.and..not.do_mp4) then
         mp3=.false.
         do_mp4=.true.
         goto 666
      endif
c
      call dynamic_retmark ! Mark after initR
      call matremark
      call retmark
c
      if ((cc.and.singles).and.l_triples) then
c     call sla_trip(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),
c    *          ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *          singles.and.cc.and..not.ccsd,af,nslv,iMasID,mygid,mytid)
c     call sort_3ext_sla(nval,idimen,bl(irecadrtt),npairs,
c    *                   ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,
c    *                   af,iMasID,mygid,mytid,nslv)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
c     call fafCloses
      call slave_3ext(nval,idimen,bl(irecadrtt),npairs,ndisktt,
     *               ncf,lbintt,thresh,byt8,nmo,vorb,af,
     *               iMasID,mygid,mytid,ichar_count,bl(ivtable),
     *               bl(im_table),bl(ivrevtable),bl(iotable),ext3_ratio)
      call fafCloses
      call make_3ext_pairs_sla(idimen,nval,af,ichar_count,bl(ivtable),
     *                         bl(im_table),bl(iotable),bl(ivrevtable),
     *                         iMasID,mygid,mytid)
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call fafCreates(itrala)
      call fafCreates(itrala)
      call fafCloses
      call new_tri_sla(nval,idimen,bl(iepsi+ncore),ccsd,qcisd,af,iMasID,
     *              mygid, mytid,bl(iotable),bl(ivtable),bl(ivrevtable),
     *              bl(im_table),bl(ijtable),ichar_count)
      endif
 999  continue
      call secund(slave_tot_cpu1)
c -- send timings back to master then move on
      slave_tot_cpu=slave_tot_cpu1-slave_tot_cpu0
      call para_initsend
      call para_pack_real(slave_tot_cpu,1)
      call para_pack_string(MY_HOSTNM,60)
      call para_send_pack(0,TxTime)
      call fafTerminate(info)
      call allclose_and_delete(.true.)
c Release memory:
      call matremark
      call retmark
      call dynamic_retmark
      end subroutine do_cimcc



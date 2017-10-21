#include "maxi_cim.h"
      subroutine cimsub
      use memory
      implicit real*8 (a-h,o-z)

      character word*4
      character*256 jobname
      logical emp2only

c-----------------------------------------------------------------------
c 0=logical option, 1=a single integer, mp2scal=p1*empanti+p2*emppar
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c-----------------------------------------------------------------------
      data ioptyp  / 1,     0,     0,     0,     0,     0,     0  /
      data word    /'prin','rimp','mp2 ','ccd ','ccsd','trip','forc'/
      data IUnit   /1/ ! unit number for checkpoint I/O
      parameter (nopt=7)
      common /job/ jobname,lenJ
      COMMON /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /GlobalCCSD/ cache,integerT,iprnt
      dimension word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
c ----------------------------------------------------------------------
c Explanation of the options:
c FORC = calculate the extra quantities needed for an MP2 gradient and
c        write them to disk       !Zhigang 3/23/2017 @UARK
c ----------------------------------------------------------------------

 10   format(72('='))
      write(iout,10)  
      write(iout,"(20x,' The Module of CIM Cluster Calculation')")
      write(iout,*) 
      
      iprnt=0
      emp2only=.true.
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
      if (ifound(1).gt.0) iprnt=iopv(1,1)
      call setival('nprint',iprnt)
      if (ifound(7).gt.0) emp2only=.false.
C      if (ifound(2).gt.0) call CIM_RIMP2
      if (ifound(3).gt.0) call CIM_MP2(emp2only,iprnt)
      if (ifound(4).gt.0) call CIM_CC(9,iprnt)
      if (ifound(5).gt.0) then
         if (ifound(6).gt.0) then
            call CIM_CC(11,iprnt)
         else
            call CIM_CC(10,iprnt)
         end if 
      end if
      write(iout,*) "End of CIM cluster calculation!"
      end subroutine cimsub


C *************************************************
C * Calculate CIM-MP2 energy for the clusters     *
C * Modify the subroutine rmp2 in PQS source code *
C * NZG_5/8/2016 @UARK                            *
C *************************************************
      subroutine CIM_MP2(emp2only,iprnt)

      use memory

      implicit real*8 (a-h,o-z)
 
      integer*1 i01
      integer*4 i0
      character*256 jobname,scrfile,filename,filname1,filname2
      character*256 filname3,filname4
      character*4 motype
      character scftype*11,cdum*20,wvfnc*20
      logical exst,restrt,dualbasis,emp2only,LastSymPair,smal

      integer,allocatable::IAN(:),Z(:),ILST(:,:),INX2(:,:)
      real*8,allocatable::XC(:,:),QA(:),BASDAT(:,:)
      real*8,allocatable::SMO(:,:),SMONEW(:,:),eorb(:),trans(:,:)
      character*8,allocatable::AtSymb(:)

      real*8,allocatable::x_QCMO(:,:,:),tij1_QCMO(:,:,:)
      real*8,allocatable::x_LMO(:,:,:,:),tij1_LMO(:,:,:,:)
      real*8,allocatable::xqcmo2(:,:),tqcmo2(:,:)
      real*8,allocatable::energypair_CIM(:,:),energyorb(:)
      real*8,allocatable::TR1(:),TR2(:,:),TC1(:),TC2(:,:)
      real*8,allocatable::PCIMA(:,:),PCIMB(:,:),ECIMA(:,:),ECIMB(:,:)

C For test
C      real*8,allocatable::Sij(:,:),SOVER(:,:),tmp(:)

      parameter (zero=0.0D0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      integer,parameter::maxgran=300*301/2
      integer*4 ipair2gran(maxgran),igran2pair(2,maxgran)

      character*8 name,names,blank
      character*1 shape,shapes(6),shapec(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat),jtype(nmat)

      common /job/ jobname,lenJ
      data IUnit/1/   ! unit number for checkpoint I/O
      dimension xintxx(9)
      dimension ipass(2,28),kpass(2,28)

      call secund(tt0)
      call elapsec(elaps0)

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
         call CIM_CC(1,iprnt)
         return
      endif

      call rdcntrl(iunit,4,'NOCC',1,nmo,rdum,cdum)
      call rdcntrl(iunit,3,'NMO',1,norb,rdum,cdum)
      call rdcntrl(iunit,4,'NCEN',1,ncen,rdum,cdum)
      call rdcntrl(iunit,5,'NCORE',1,ncore,rdum,cdum)
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
C  Read the basis set data
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')  ! Location of INX array
      call getival('nsh',nsh)
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

C -- get the overlap matrix for normalizing the orbitals
C      allocate(tmp(127008),Sij(ncf,ncf),SOVER(ncf,ncf))
C      call inton2(0,natoms,Sij,INX2,INX2,0,0,BASDAT,BASDAT,XC,IAN,
C     $            ncs,ncs,ncf,ncf,tmp)
C      call ReorderFock2(ncf,ncf,Z(1),Sij,SOVER)

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

C -- Renormalize the orbitals
C      call normorb(ncf,norb,SMONew,SOVER)

C -- Write the reordered MOs to XXX.mos file
      filename(lenM+1:lenM+1)='2'
      lenM=lenM+1
      call WriteMOS(ncf,norb,SMONEW,eorb,.true.,lenM,filename,itype)
      
CC Check if the orbitals are orthogonal to each other
C      
C      write(io,*) "Check the orthonormalization!"
C      KB=nocc
C      allocate(Sij(KB,KB))
C      Sij=0.0D0
C      call MO_over(Sij,SMONew,
C      do j=1,KB
C         do k=1,KB
C            if (j==k) then
C               if (Sij(j,k)-1.0D0>1.0D-6) then
C                  write(io,*) "Fail to normalize!"
C               end if
C            else
C               if (Sij(k,j)>1.0D-6) then
C                  write(io,*) Sij(k,j)
C                  write(io,*) "Fail to orthonormalize!"
C               end if
C            end if
C         end do
C      end do
C      deallocate(Sij)
C
C      stop
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')

c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c  basis function symmetry pairs are in inx(ifp)
c
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
C + For the CIM cluster calculation, most of the parameters will + 
C + use the default values                                       +
C ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c  THRESHOLD (integral threshold)
      thresh=1.0d-10
      call setrval('thresh',thresh)
      restrt=.false.
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
      mp2only=1
      if (emp2only==.false.) then
         mp2only=-1
         filname3=scrfile(1:len)//'.Tij'
         len3 = len+4
         filname4=scrfile(1:len)//'.Kov'
         len4 = len+4
      end if
      call setival('mp2only',mp2only)
c .........................................................................
      smal=.true.       ! always use small option

C Set the disk storage as 200GB
      xmaxdisk=25000000000.0d0

c  'PMIJ' print matrix; first number is 1 (AO Xchange) or 2
c   (both AO & MO) second number: pair to be printed
      ipmij=0
      ipairpr=0

C Set some values which are not needed in energy calculation but
C may be needed in the force calculation -NZG_3/30/2017 @UARK
      iscs=0
      call setival('iscs',iscs)

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
      ncs_sm=ncs
      ncf_sm=ncf
      call set_mpres(bl(mpres_in_bg),ncs)
c-----------------------------------------------------------
      write(iout,"(' MP2 integral thresh   = ',E14.4)") thresh
      write(iout,"(' Max disk storage (MW) = ',F14.4)") xmaxdisk*1.d-6
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
      np4=4
      call matdef('cano','q',ncf,ncf)
      call matdef('epsi','d',ncf,ncf)
      call matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      idics=mataddr('dsmx')
c-------------------------------------------------------
c
      CALL ReadMOS(ncf,bl(icano),bl(iepsi),.True.,lenM,
     $             filename(1:lenM),itype,IErr)
      If(IErr.NE.0) Call nerror(3,'CIM Cluster RMP2 Module',
     $   'MOS File Does Not Exist',0,0)
c initialize the two-el. integral program
      thint=thresh
      iforwhat=5
      call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             scftype,xintxx, nblocks,maxbuffer,maxlabels)

c  In CIM cluster, all the occ MOs are correlated -Fail now
c  If force is to be calculated, we have core orbitals.
c  -NZG_4/28/2017
      nfirst=ncore+1
      nlast=nmo
      nval=nmo-ncore
c
c check the maximum values of BOTH
c the occupied and virtual eigenvectors
c
      iocco=mataddr('cano')
      iepsi=mataddr('epsi')
      call check_coef(ncs,ncf,ncore,nmo,nval,bl(iocco))

c---------------------------------------------------------------------
c Construct "screening density" matrix to be used to estimate the
c importance of integrals. 
c In CIM calculation, we use canonical orbitals for integral screening
c---------------------------------------------------------------------
 
      write(iout,*) 'Occupied orbitals are used for integral screening'
      dualbasis=.false.
      call DmxMakeC(dualbasis, nval, nmo, ncf, ncs, np4, inx2, iprnt)
c
c  establish orbital symmetry characters
      if (nsym.gt.0) then
         call getint(nval*nsym,iorbpair)
         call getmem(ncf,itmp)
         call OrbSymPair(nsym, ncf ,   bl(icano), bl(ifp) ,nfirst,
     1                   nlast,iprnt, bl(itmp),  bl(iorbpair),iout)
         call retmem(1)
      endif
c
      lastvirt=norb
c  print the data of the calculation early
      write(iout,'(" Number of contracted functions =",i8,/,
     1             " Number of correlated orbitals  =",i8,/,
     2             " Number of virtual orbitals     =",i8)')
     3               ncf,nval,lastvirt-nmo
      write(iout,*) ' '
      call f_lush(iout)
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
C      if (ncore.gt.0) call matsub('occa','cano',1,nmo)
C Always define 'occa' matrix -NZG_5/17/2017
      call matsub('occa','cano',1,nmo)
      call matsub('occu','cano',nfirst,nlast)
      call matsub('virt','cano',nmo+1,lastvirt)
      ioccu=mataddr('occu')
c.................................................
      if (iprnt.gt.3) then
         call matprint ('occu',6)
         call matprint ('virt',6)
      end if
c .................................................
c  open direct access scratch file to store half-transformed integrals
c
c  the length of a block is given in bytes.
      ints = 4             ! size of an integer in bytes
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
      if (restrt) then
         mulam=ncf*(ncf-1)/2
         mulamd=ncf
         call symmoff
         goto 50
      end if
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
      ncf2=ncf*ncf
      nskipped=0
      do ics=ncs,1,-1
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf*ncf*ics_size
         do kcs=ics,1,-1
            if (LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c  of each (ics,kcs) pair, calculate only the last one
c  (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
c
            call get_shell_size(bl(ictr),kcs,kcs_size)
c
            lmp2_size=lmp2_siz1*kcs_size
c
c check if size of the mp2 integral buffer is not too big
c if it is then split over kcs (and possibly ics)
c
            call check_size1(lmp2_size,ncf,nval,nval,ntimes)
c
            call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                      ntimes,Ipass,Kpass,Itimes,Ktimes)
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
                  call int_lmp2(bl, bl(ictr), thresh, ics, icf1,
     1                          icf2, kcs, kcf1, kcf2, bl(mapf2s),
     2                    bl(idics),iprnt,bl(lmp2int),nintotal,nrow,
     3                    ncol, bl(irow),bl(icol),bl(lzero))
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
                  call TransOneShell(ncf,   ncs,    nval,   ics,    kcs,
     *                         icf1,   icf2,   kcf1,   kcf2,   ndisk,
     1                    bl(ictr),bl(lmp2int),bl(ioccu),iprnt,
     2                    thresh,bl(ixadr),bl(ihalftra), irec, mulam,
     3                    mulamd,   nrow,     ncol,  bl(irow), bl(icol),
     4                    bl(irow1),bl(icol1),ENonZero,smal)
c
                  call secund(tt4)
                  call elapsec(telap4)
                  ttrans=ttrans+tt4-tt3
                  elaptrans=elaptrans+telap4-telap3
c......................................................................
                  call retmem(1)          ! lmp2int
               enddo    ! over ktime (selected kcf belonging to kcs shell )
            enddo    ! over itime (selected icf belonging to ics shell )
         enddo     !  over kcs shell
      enddo     !  over ics shell
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
C NZG
      write(6,*) 'HERE'
      
c......................................................................
      write(iout,61)
   61 format(' CPU & Elapsed timings in the MP2 module ')
      write(iout,*) '  '
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format(' Integrals for MP2  =',f8.2,' and ',f8.2,' minutes')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format(' First half transf. =',f8.2,' and ',f8.2,' minutes')
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
      lastvirt=norb
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
      mem=memavail-memtrans-1000000  
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
      if (lbin.gt.ncf*(ncf+1)) lbin=ncf*(ncf+1)
CTJ granules of bins will be written
C*******************************************************************************
      do ! over lbin
         igranulesize=mem/(lbin*9/8+1+ncf*ncf)
         if (igranulesize.gt.npairs) igranulesize=npairs
cc        write(6,'(A,I5)') "Max granule size: ",igranulesize
         do
            lrec =(lbin     +8*lbin       )*igranulesize ! record length in bytes
c Limit record size to 8 MB.
            if (lrec.gt.8388608) then
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
         if (lbin.lt.50) then
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
      if (iprnt.gt.4) then
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
      open(unit=ndisk2,file=filname2,access='direct',recl=lrec)
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
      if(xmaxd.lt.10000000) then
        maxdisk=xmaxdisk
        call nerror(8,'RMP2 Module',
     1 'Less than 10 MW disk storage remains after half transformation',
     2  int(xmaxdisk),int(xmaxd))
      end if
      xldisk=dble(lbin)*dble(nbins)
      if (xldisk.gt.xmaxd) then
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
c
      if (.not.emp2only) then
         write(iout,*) 
         write(iout,*) '*** Quantities Stored for MP2 gradient ***'
         write(iout,*) 
         nvmo=nval
         if (ncore.gt.0) nvmo=nmo
         llrec=nvirt*nvirt
         krec=5*llrec
         open(unit=ndisk3,file=filname3(1:len3),form='unformatted',
     &        access='direct',recl=krec)
         npars=nval*(nval+1)/2
c
c -- also open files for virtual occupied blocks of Kij
         lxrec=nvmo*nvirt*2
         kxrec=5*lxrec
         open(unit=ndisk5,file=filname4(1:len4),form='unformatted',
     &        access='direct',recl=kxrec)
         
         call matdef('tvir','r',nvirt,ncf)
         call matpose2('virt','tvir','n')
         call f_lush(6)
      end if

      if(iprnt.gt.1) then
         write(iout,'(" Yoshimine BinSort ",i9," bins",i7," long ",
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
      if (iprnt.gt.2) then
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

      allocate(x_QCMO(nmo*(nmo+1)/2,nvirt,nvirt))
      allocate(tij1_QCMO(nmo*(nmo+1)/2,nvirt,nvirt))
      x_QCMO=0.0D0; tij1_QCMO=0.0D0

      call SortTransf_CIM(
     1   ncf,         nval,       npairs,     ndisk,      ndisk2,
     2   nrec,        lbin,       nbins,      igrfirst,   igrlast,
     3   iout,       emp2only,    'virt',     bl(i1pair), bl(int1),
     4   bl(icounter),bl(ibinctr),lrec,       nfirst,     nmo,
     5   nvirt,       thresh,     nsym,       bl(iorbpair),bl(ifp),
     6   trans,       iprnt,      ipmij,      ipairpr,    bl(ixmos),
     7   bl(iepsi),   emp2,       tnorm,      emps,       empt,
     8   empanti,     emppar,     TBinSort,   ndisk3,     krec,
     9   tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     $   ttkij,       tmp2e,      ncore,      igran2pair,ipair2gran,
     1   igranulesize,ncen,       x_QCMO,     tij1_QCMO)
      igrfirst=NOdd+1
      do ibatches=1,nbatches
c  the first (odd) batch is not included in nbatches
         igrlast=igrfirst-1+NGransInBatch
         call SortTransf_CIM(
     1       ncf,         nval,       npairs,     ndisk,      ndisk2,
     2       nrec,        lbin,       nbins,      igrfirst,    igrlast,
     3       iout,       emp2only,    'virt',     bl(i1pair), bl(int1),
     4       bl(icounter),bl(ibinctr),lrec,       nfirst,     nmo,
     5       nvirt,       thresh,     nsym,       bl(iorbpair),bl(ifp),
     6       trans,       iprnt,      ipmij,      ipairpr,    bl(ixmos),
     7       bl(iepsi),   emp2,       tnorm,      emps,       empt,
     8       empanti,     emppar,     TBinSort,   ndisk3,     krec,
     9       tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     $       ttkij,       tmp2e,      ncore,      igran2pair,ipair2gran,
     1       igranulesize,ncen,       x_QCMO,     tij1_QCMO)
         igrfirst=igrlast+1
      end do

C For the CIM calculation, we only need to transform the integrals and
C amplitudes from QCMOs to the central LMOs.
C So the indices of some matrices have been changed from nmo to ncen.
C If in the future, the LCCE-CIM is implemented, we still need the
C energy of the whole subsystem. Maybe in this case, we do not need any
C transformation, but directly use the QCMO results.
C NZG_5/8/2017 @UARK
      allocate(x_LMO(ncen,nval,nvirt,nvirt))
      allocate(tij1_LMO(ncen,nval,nvirt,nvirt))
      x_LMO=0.0D0; tij1_LMO=0.0D0

      do k=1,nvirt
         do l=1,nvirt
            allocate(xqcmo2(nval,nval),tqcmo2(nval,nval))
            do j=nfirst,nlast
               do ii=nfirst,nlast
                  if (ii>=j) then
                     ij=ii*(ii-1)/2+j
                     xqcmo2(ii-ncore,j-ncore)=x_QCMO(ij,k,l)
                     tqcmo2(ii-ncore,j-ncore)=tij1_QCMO(ij,k,l)
                  else
                     ij=j*(j-1)/2+ii
                     xqcmo2(ii-ncore,j-ncore)=x_QCMO(ij,l,k)
                     tqcmo2(ii-ncore,j-ncore)=tij1_QCMO(ij,l,k)
                  endif
               enddo
            enddo
            call dgemm('T','N',ncen,nval,nval,1.0D0,trans(:,1:ncen),
     &                 nval,xqcmo2,nval,0.0D0,x_LMO(:,:,k,l),ncen)
            call dgemm('T','N',ncen,nval,nval,1.0D0,trans(:,1:ncen),
     &                 nval,tqcmo2,nval,0.0D0,tij1_LMO(:,:,k,l),ncen)
            deallocate(xqcmo2,tqcmo2)
         enddo
      enddo

      allocate(energypair_CIM(ncen,nval),energyorb(ncen))
      energypair_CIM=0.0D0; energyorb=0.0D0
      do i=1,ncen
         do j=1,nval
            do k=1,nvirt
               do l=1,nvirt
                  x=x_LMO(i,j,k,l)
                  tij1=tij1_LMO(i,j,k,l)
                  energypair_CIM(i,j)=energypair_CIM(i,j)+tij1*x
               enddo
            enddo
         enddo
      enddo

      do i=1,ncen
         do j=1,nval
            energyorb(i)=energyorb(i)+energypair_CIM(i,j)
         enddo
      enddo

      Ecen=DSUM(ncen,energyorb,1)
C      Ecen=PSUM(nmo,energyorb,NCEN)


C      allocate(x_LMO(nmo,nmo,nvirt,nvirt),tij1_LMO(nmo,nmo,nvirt,nvirt))
C      x_LMO=0.0D0; tij1_LMO=0.0D0
C      do k=1,nvirt
C         do l=1,nvirt
C            allocate(xqcmo2(nmo,nmo),tqcmo2(nmo,nmo))
C            do j=1,nmo
C               do ii=1,nmo
C                  if (ii>=j) then
C                     ij=ii*(ii-1)/2+j
C                     xqcmo2(ii,j)=x_QCMO(ij,k,l)
C                     tqcmo2(ii,j)=tij1_QCMO(ij,k,l)
C                  else
C                     ij=j*(j-1)/2+ii
C                     xqcmo2(ii,j)=x_QCMO(ij,l,k)
C                     tqcmo2(ii,j)=tij1_QCMO(ij,l,k)
C                  endif
C               enddo
C            enddo
C            call dgemm('T','N',nmo,nmo,nmo,1.0D0,trans,nmo,
C     &                 xqcmo2,nmo,0.0D0,x_LMO(:,:,k,l),nmo)
C            call dgemm('T','N',nmo,nmo,nmo,1.0D0,trans,nmo,
C     &                 tqcmo2,nmo,0.0D0,tij1_LMO(:,:,k,l),nmo)
C            deallocate(xqcmo2,tqcmo2)
C         enddo
C      enddo
C
C      allocate(energypair_CIM(nmo,nmo),energyorb(nmo))
C      energypair_CIM=0.0D0; energyorb=0.0D0
C      do i=1,nmo
C         do j=1,nmo
C            do k=1,nvirt
C               do l=1,nvirt
C                  x=x_LMO(i,j,k,l)
C                  tij1=tij1_LMO(i,j,k,l)
C                  energypair_CIM(i,j)=energypair_CIM(i,j)+tij1*x
C               enddo
C            enddo
C         enddo
C      enddo
C
C      do i=1,nmo
C         do j=1,nmo
C            energyorb(i)=energyorb(i)+energypair_CIM(i,j)
C         enddo
C      enddo
C
C      Ecluster=DSUM(nmo,energyorb,1)
C      Ecen=PSUM(nmo,energyorb,NCEN)

      write(6,*) 
      write(6,*) "Orbital energies for central LMOs:" 
      write(6,'(5F15.8)') energyorb 
      write(6,'(A,F15.9)') ' E(CORR-CENTRAL)=',Ecen
C      write(6,'(A,F15.9)') ' E(CORR-SUBSYS) =',Ecluster
      write(6,*)

      deallocate(x_QCMO,tij1_QCMO,x_LMO,tij1_LMO)
      deallocate(energypair_CIM,energyorb)

      IF (.not.emp2only) THEN
         call matrem('tvir')
         if (iprnt.gt.1) then
            write(iout,1010) tmbtr/sixty,ttkij/sixty,tmp2e/sixty
 1010       format(' Time for saving <Tij> in MO basis:',f8.2,'min',/,
     $             ' Time for saving <Kov> block:      ',f8.2,'min',/,
     $             ' Time for MP2 energy calculation:  ',f8.2,'min')
            write(iout,*) ' <Tij>: ',nxrec,' records written on unit',
     1                    ndisk3,' record length ',5*llrec
            write(iout,*) ' <Kov>: ',mxrec,' records written on unit',
     1                    ndisk5,' record length ',5*lxrec
         endif
c -- close all files
         close(unit=ndisk3,status='keep')
         close(unit=ndisk5,status='keep')
      ENDIF

c---------------------------------------------------------------------
c print out timings :
      call secund(tt3)
      call elapsec(elaps3)
c........................................................
      write(91,62)
      write(91,64) tinteg/sixty, elapint/sixty
c
      write(iout,704) TBinSort/sixty,nbins
  704 format(' Binary sorting     =             ',f8.2,' minutes'/
     *       ' (number of bins    =',i8,')' )
      write(iout,705) (tt3-tt2)/sixty,(elaps3-elaps2)/sixty
  705 format(' sort & virt.trans. =',f8.2,' and ',f8.2,' minutes')
      write(iout,*) ' '
c Print out the final results :
c
      tnorm=sqrt(tnorm)
c
      esing_exc=zero
c
c
      call f_lush(6)
c -- set wavefunction type
      wvfnc = 'RMP2'
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
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
      call symmon
c---------------------------------------------------------------------
      call matrem('xmos')
c---------------------------------------------------------------------
      call setival('ncore',ncore)
      call setival('nstrip2',nstrip2)
      call setival('mrcpf',mrcpf)
      call setival('lrcpf',lrcpf)
      call setrval('thrsx',thresh)
c---------------------------------------------------------------------
C      call memory_status('end of cim_mp2')
c---------------------------------------------------------------------

c -- Entering MP2 gradient calculation module for CIM cluster
c -- Zhigang 3/23/2017 @UARK

      if (.not.emp2only) call force2_MP2_CIM(natoms,nmo,nvirt,ncen,
     &                                       .true.,ictr,iprnt,trans)

      call retmem(5)
      call matremark
      call retmark

      end subroutine CIM_MP2


C ***************************************************************
C * Modify the subroutine SortTransf in PQS for CIM calculation *
C * NZG_5/16/2016 @UARK                                         *
C ***************************************************************
      subroutine SortTransf_CIM(
     1   ncf,         nval,       npairs,     ndisk,      ndisk2,
     1   nrec,        lbin,       nbins,      igrfirst,   igrlast,
     2   iout,       emp2only,    virt,       i1pair,     int1,
     3   icounter,    ibinctr,    lrec,       nfirst,     nmo,
     4   nvirt,       thresh,     nsym,       iorbpair,   ifunpair,
     5   trans,       iprnt,      ipmij,      ipairpr,    xmos,
     6   epsi,        emp2,       tnorm,      emps,       empt,
     7   empanti,     emppar,     TBinSort,   ndisk3,     krec,
     8   tmbtr,       nxrec,      ndisk5,     kxrec,      mxrec,
     9   ttkij,       tmp2e,      ncore, igran2pair,ipair2gran,
     $   igrsize,     ncen,       x_QCMO,     tij1_QCMO)

      use memory

      implicit real*8 (a-h,o-z)
      integer*4 i1pair(*)
      dimension icounter(*),ibinctr(*),epsi(*)
      dimension iorbpair(nsym,nfirst:nmo),ifunpair(7,ncf)
      dimension ipair(7),jpair(7),trans(nmo,nmo)
      integer*4 igran2pair(2,*),ipair2gran(*)
      integer*1 int1(*)
      logical emp2only
      character*(*) virt
      parameter(zero=0.0d0,one=1.0d0,four=4.0d0)
      real*8 x_QCMO(nmo*(nmo+1)/2,nvirt,nvirt)
      real*8 tij1_QCMO(nmo*(nmo+1)/2,nvirt,nvirt)
c
      call mmark()
c  reserve memory for bins
      call getint_4(npairs*lbin*2,ibins)
      call getint_1(npairs*lbin,ibin1)
c
      call elapsec(elap1)
      call BinSort(ncf,      nval,     npairs, ndisk,    ndisk2,
     1             nrec,     lbin,     nbins,  igrfirst,  igrlast,
     2             bl(ibins),bl(ibin1),i1pair, int1,    icounter,
     3             ibinctr,lrec,igran2pair,ipair2gran,igrsize)
      call elapsec(elap2)
      TBinSort=TBinSort+elap2-elap1
      call retmark()
      call mmark()
      call getmem(ncf*ncf*igrsize,ixmat)
      call getint_4(igrsize*lbin*2,ibins)
      call getint_1(igrsize*lbin,ibin1)
c
c -- Additional stuff for RMP2 gradients
      If (.NOT.emp2only) Then
         nval=nmo-nfirst+1
         nvmo=nval
         if (ncore.gt.0) nvmo=nmo
         call matdef('kijvo','r',nvirt,nvmo)
         ikijvo=mataddr('kijvo')
         lengvo=nvirt*nvmo
      EndIf
c
c  read back bins, extract 1 exchange matrix, and transform that
      ij=0
      do igran=igrfirst,igrlast
         call ReadBin1(ncf,       iorbpair,       igran,  ndisk2,
     1                 lbin,      nbins, ibinctr, thresh, bl(ibins),
     2                 bl(ibin1), nsym,  ifunpair, bl(ixmat),igran2pair,
     3                 ipair2gran,igrsize)
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

c -- Additional stuff for MP2 gradients
            If (.NOT.emp2only) Then
c -- Kij in AO basis now in xmat
c -- transform to virtual- occupied block and save for gradients
               call elapsec(tkij1)
               mxrec=mxrec+1
               call SaveKij(ncf,nvmo,nvirt,bl(ikijvo),bl(ibins),
     1                      bl(ibin1),ij,thresh,lengvo,ndisk5,ncore)
               call elapsec(tkij2)
               ttkij=ttkij+tkij2-tkij1
            EndIf

            call TransVirt(ncf,i,j,iprnt,ipmij,ipairpr,virt)
            call IntExpr2(ncf,nmo,nvirt,i,j,epsi,xmos,x_QCMO,tij1_QCMO)

c-------------------------------------------------------------------
c -- Additional stuff for MP2 gradients
            If (.NOT.emp2only) Then
               call elapsec(tt1)
               tmp2e=tmp2e+tt1-tkij2
               nxrec=nxrec+1
               call SaveTij(ncf,    nmo,    nval,   nvirt,  i,
     1                      j,      epsi,   thresh, bl(ibins),bl(ibin1),
     2                      xmos,   ndisk3, ij)
               call elapsec(tt2)
               tmbtr=tmbtr+tt2-tt1
            EndIf
c--------------------------------------------------------------------
            call matdisc('xmat')
         end do
      end do

      if (.not.emp2only) call matrem('kijvo')

      call retmark()
      end subroutine SortTransf_CIM


C **********************************************************
C * Express the two electron integrals to four-index array *
C * It is needed by the next step for transformation       *
C * NZG_5/17/2016 @UARK                                    *
C **********************************************************
      subroutine IntExpr(nmo,nvirt,i,j,xmos,tmpint)
      implicit real*8 (a-h,o-z)
      
      dimension xmos(nvirt,nvirt),tmpint(nmo*(nmo+1)/2,nvirt,nvirt)
      do ib=1,nvirt
         do ia=1,nvirt
            ij=i*(i-1)/2+j
            tmpint(ij,ia,ib)=xmos(ia,ib)
         end do
      end do
      
      end subroutine IntExpr


C *******************************************************************
C * A new subroutine for transforming the two-electron integrals to *
C * four-index array by modifying the subroutine PairEner from PQS  *
C * source code                                                     *
C * NZG_11/18/2016 @NJU                                             *
C *******************************************************************
      subroutine IntExpr2(ncf,nmo,nvirt,i,j,e,x,x_QCMO,tij1_QCMO)
      implicit real*8 (a-h,o-z)
c  Arguments (INPUT):

c  nvirt:          # of virtual orbitals
c  i,j    :        occupied orbitals defining the pair; e(i,j) is calc.
c  x(nvirt,nvirt): x(a,b)=(a i|j b) MO integrals in matrix form
c  e:              orbital energies
      dimension x(nvirt,nvirt),e(ncf)
      dimension x_QCMO(nmo*(nmo+1)/2,nvirt,nvirt)
      dimension tij1_QCMO(nmo*(nmo+1)/2,nvirt,nvirt)
      parameter(one=1.0D0)

      ei=e(i)
      ej=e(j)
      do ia=1,nvirt
         ea=e(nmo+ia)
         do ib=1,nvirt
            rdenom=one/(ei+ej-ea-e(nmo+ib))
            tij1=x(ia,ib)*rdenom
            ij=i*(i-1)/2+j
            x_QCMO(ij,ia,ib)=x(ia,ib)+x(ia,ib)-x(ib,ia)
            tij1_QCMO(ij,ia,ib)=tij1
         enddo
      enddo

C      if (i/=j) then
C         do ia=1,nvirt
C            do ib=1,nvirt
C               x_QCMO(j,i,ia,ib)=x_QCMO(i,j,ib,ia)
C               tij1_QCMO(j,i,ia,ib)=tij1_QCMO(i,j,ib,ia)
C            enddo
C         enddo
C      endif

      end subroutine IntExpr2


      double precision function DSUM(n,x,incr)
      implicit none

      integer n,incr,i,ii
      double precision x(*)
      double precision sum

      sum = 0.0D0
      if(incr.eq.1) then
         do i = 1,n
            sum = sum + x(i)
         end do
      else
         ii = 1
         do i = 1,n
            sum = sum + x(ii)
            ii = ii + incr
         end do
      end if
      DSUM = sum
      end function DSUM


C *********************************************************************
C * function PSUM is different from the original one in GAMESS        *
C * In this version, central MOs are always arranged at the beginning *
C * So we only need the number of central orbitals                    *
C * NZG_5/21/2016 @UARK                                               *
C *********************************************************************
      double precision function PSUM(n,x,ncen)
      implicit none
      integer n,ncen,i
      double precision sum,x(n)
      double precision ZERO
      parameter (ZERO=0.0D+00)

      sum = ZERO
      do i = 1,ncen
         sum = sum + x(i)
      end do
      PSUM = sum
      end function PSUM

      subroutine rwrit(iunit, key, n, A)
      implicit none
      character key*(*)
      integer iunit,n,L,k
      double precision A(n)
C
      L=len(key)
      write(iunit,*) key(1:L), n
      write(iunit,'(1P,5E20.12)') (A(k),k=1,n)
      write(iunit,*) '$END'
      end subroutine rwrit


C ************************************************************
C * Subroutine for CIM-CC calculation                        *
C * Modify the source code from coulomb.F in PQS source code *
C * NZG_5/23/2016 @UARK                                      *
C * NZG_10/19/2016 @NJU  Modify to use this subroutine to do *
C * CIM-LMP2 calculation                                     *
C ************************************************************
      subroutine CIM_CC(imethod,iprint)

      use memory
      use kinds
c testing program for Jij Coulomb integrals 
      implicit real*8 (a-h,o-z)
      logical equal,gauss_seidel,nodisk 
      logical variational,linear,cc,ccsd,norecalc,byt8,singles
      logical loca,omit,log_diis,mp2last,small,vorb,reset
      logical cache,calchunks,l_triples,qcisd
      dimension xintxx(9)
      common /job/jobname,lenJ
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $              filname3,filname4
      character*6 method
      character scftype*11,wvfnc*20,ch3*3
      logical nofr,exst,restrt,dualbasis,smal,mp3,cep0,cep2,mp2,af,mp4
      logical do_mp4
      character*2 int_kind
      dimension ifilestat(13),xnmo(2)
      real*8 elaps(5)
      integer*4 iseed
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /extrone_tmp/ af,islaves,ndisktrc,ndisktrx,ndisktre,
     *                     ndisktrtx,ndisktrtc,ndisktrtt,mygid,nbf
      dimension cctimesij(4),cctimes(4) 
      parameter (itim_no=20)
      dimension ccsdcpu(itim_no),ccsdela(itim_no)
      logical integerT
      common /GlobalCCSD/ cache,integerT,iprnt
      logical trans,signum,docansym,converged
      character*80, parameter :: spc='                             '
      integer,allocatable::IAN(:),Z(:),ILST(:,:),INX2(:,:)
      real*8,allocatable::XC(:,:),QA(:),BASDAT(:,:)
      real*8,allocatable::SMO(:,:),SMONEW(:,:),eorb(:),TX(:,:)
      real*8,allocatable::dtmp(:),Fock1(:,:)
      character*8,allocatable::AtSymb(:)
      real*8,allocatable::exc_QCMO(:,:,:,:),exc_LMO(:,:,:,:)
      real*8,allocatable::coef_QCMO(:,:,:,:),coef_LMO(:,:,:,:)
      real*8,allocatable::resi_QCMO(:,:,:,:),resi_LMO(:,:,:,:)
      real*8,allocatable::energypair_CIM(:,:),energyorb(:)
      real*8,allocatable::EPairQCMO(:,:)
      real*8,allocatable::CIMPairEnergy(:,:),CIMOrbEnergy(:)
      real*8,allocatable::Sij(:,:)
c-----------------------------------------------------------------------
c
c  Explanation of the options:
c
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
      data IUnit/1/
c
      strace=0.0d0
      ik=0
      af=.false.
      islaves=1
      ndiske=-1
      lbine=-1
c
      call signal_default(11)
      call secund(tt0)
      call elapsec(elaps0)
      iseed = 0
      call srand(iseed)
      call reset_counters
      call setival('ccprint',iprnt)
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c Defaults:
      ndisk_mp3=83
      nodisk=.false.
      gauss_seidel=.true.
      loca=.false.
      omit=.false.
      mp2=.false.
      mp3=.false.
      mp4=.false.
      do_mp4=.false.
      cep0=.false.
      cep2=.false.
      linear=.false.
      smal=.true.
      small=.true.
      cc=.false.
      qcisd=.false.
      norecalc=.true.
      byt8=.true.
      singles=.false.
      ccsd=.false.
      l_triples=.false.
      vorb=.true.
      cache=.false.
      variational=.true.
      integerT=.false.
      xmaxdisk=2500000000.0d0      ! currently not active
      maxiter=50
      log_diis=.true.
      length_diis=6
      ifirst_diis_iter=1
      Ethresh=1.0d-6
      Wavethresh=1.0d-6
c
      call matreset
c
c -- determine correlation method
      If (imethod.EQ.4) then
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
C If one use this subroutine to calculate CIM-MP2 energy, that means
C it is CIM-LMP2 calculation.
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
c
c------------------------------------------------------------------
c -- THRE (integral threshold)
      thresh=1.0d-12
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
c -- MEMO  !NZG: Caution about the number, think about it later
      max_dyn_mem=100000000000

      call dynamic_init(max_dyn_mem,bl(1))
      call dynamic_mmark
      call mmark
      call matmark
c
c -- CHUNk
      calchunks=.true.
c
C  Read the basis set data
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')  ! Location of INX array
      call getival('nsh',nsh)
c-----------------------------------------------------------
c -- check basis set (G is highest angular momentum basis function coded)
      call ChkAngMom(ncs,bl(ictr),2)
c-----------------------------------------------------------

      erhf=0.0D0 ! HF energy; 
C                  In CIM cluster calculation, we just set it as zero
c----------------------------------------------
c
C Read the number of orbitals of the cluster
c **Warning** The number of occ MOs cannot be determined by nele/2
C nmo:   number of occupied orbitals
C nbf: number of occupied and virtual orbitals
C      i.e. number of nonredundant basis functions (not true in CIM) 
C TX:  transformation matrix that transform QCMO to LMO
      open(unit=iunit,file=jobname(1:lenJ)//'.cim',form='formatted',
     &     status='old')
      call rdcntrl(iunit,4,'NOCC',1,nmo,rdum,cdum)
      call rdcntrl(iunit,3,'NMO',1,nbf,rdum,cdum)
      call rdcntrl(iunit,4,'NCEN',1,ncen,rdum,cdum)

C If one use this subroutine to calculate CIM-MP2 energy, that means
C it is CIM-LMP2 calculation. So for local calculation, we do not need
C the transformation matrix
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
      norb=nbf
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

      call dynamic_lock(bl(icano),i)
      call dynamic_lock(bl(iepsi),i)
c
      call dynamic_getmem(ncf*ncf,ifockAO)
      call ReorderFock2(ncf,ncf,Z(1),Fock1,bl(ifockAO))
      deallocate(dtmp,Fock1,Z)

      write(6,'(A,I5)') 'DIIS expansion length: ',length_diis
      write(6,'(A,I5)') 'DIIS evaluation start: ',ifirst_diis_iter
c
      write(6,'(A,F6.1)')  'Level shift set to ',shift
      write(6,'(A,E10.3)') 'Threshold for energy convergence:', Ethresh
      write(6,'(A,E10.3)') 
     &      'Threshold for wvfn.  convergence:', Wavethresh
C
C  Determine the orbitals to be correlated
C  In CIM cluster calculation, all the occupied orbitals are correlated
      nfirst=1
      nlast=nmo
      nval=nmo
      ncore=0
      nend=nbf
c
      call matsub('genvirt','cano',nmo+1,nbf)
      nvirt=nbf-nmo      ! number of virtual orbitals
c
      idimen=nvirt
      call twoelinit
c
      call open_mp3_resid(ndisk_mp3,idimen,.false.,mp4)
c  print the data of the calculation early
      write(iout,'("Number of contracted functions =",i8,/,
     1             "Number of correlated orbitals  =",i8,/,
     2             "Number of virtual orbitals     =",i8)')
     3             ncf,nval,nvirt
      write(iout,*) ' '
      call flush(iout)
c
c Check orbital symmetry:
      call dynamic_matdef('overlap','q',ncf,ncf)
      ioverlap=mataddr('overlap')
      call OverlapBuilder(ioverlap)

C      allocate(Sij(ncf,ncf))
C      call MO_over(Sij,bl(icano),bl(ioverlap),ncf,ncf)
C      write(6,*) Sij
C
C      stop
      call dynamic_unlock(bl(icano),i)
      call dynamic_unlock(bl(iepsi),i)
      norb=nvirt+nval
      icorr=icano
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      call symm_memory_allocator(ncf)
      ichar_count=1
      if (.not.loca) then
         call dynamic_getmem(ncf*ncf,iscr1)
         call dynamic_getmem(ncf*ncf,iscr2)
         call dynamic_getmem(ncf*ncf,iscr3)
         call adv_symmetrize(nsym,ncf,ncore,bl(ifp),bl(icano),
     &                       bl(ioverlap),bl(iscr1),bl(iscr2),bl(iscr3),
     &                       bl(iepsi))
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
      docansym=.false.
      if (nsym.gt.0.and..not.loca) docansym=.true.
      if (.not.loca) then
         call dynamic_getmem(8*7,ichatacters)
         call dynamic_getmem(8*8,im_table)
         call dynamic_getmem(8*(nval+1),iotable)
         call dynamic_getmem(nval,iorevtable)
         call dynamic_getmem(8*(nvirt+1),ivtable)
         call dynamic_getmem(nvirt,ivrevtable)
         call dynamic_getmem((nval*nval+1)*8,ijtable)
         call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                      bl(ichatacters),ichar_count,bl(im_table))
         call dynamic_mmark()
         call dynamic_getmem(2*8,      ibassym)
         call dynamic_getmem(ncf*ncf,  imat1)
         call dynamic_getmem(ncf*ncf*8,imat2)
         call symm_adapted(ncf,nsym,ichar_count,bl(ichatacters),bl(ifp),
     *                     bl(imat1),bl(imat2),bl(ibassym))
         call dynamic_retmem(1)
         if (ncf.eq.nbf) call Fock_Diag(ncf,ifockAO,ioverlap,icano,
     &                                  iepsi,bl(imat1),bl(ibassym),nmo)
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
      if (loca) then
         call symmetrize_energies(ncf,nsym,nval, bl(ivalpair),
     *                            bl(iepsi+ncore))
         call symmetrizer2(nsym,  ncf,  nval,  bl(ifp),  bl(ivalpair),
     *                     bl(icorr))
         call symmetrizer2(nsym,  ncf,  nvirt,  bl(ifp),  bl(ivirpair),
     *                     bl(icorr+nval*ncf))
      endif
!      call dynamic_lock(bl(iepsi),i)
      call pair_s_initializer(vorb,     nval,     nsym, ifp, ifp1,
     *                        ivalpair, ivirpair, loca,iorevtable,
     *                        ivrevtable,im_table,ichar_count,
     *                        ivtable)
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
c                                                 pair has, if given pair
c                                                 is already an image, = 0
c                                                 npairs+1 because of DIIS
      call izeroit(bl(ipairimages),npairs+1)
      call int_array_write(bl(ipairimages),npairs+1,1)
c
c     call orb_pairs_check_ao(ncf,bl(icorr),1d-14)
      call dynamic_lock(bl(icano),i)
      do i=1,nval
         do j=1,i
            ij=ij+1
            call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                         ns, trans, signum)
            if (ns.eq.0) then
               if (iprnt.ge.3) 
     &            write(6,'(A,I2,A,I2,I5)') "Pair: ",i,",",j,ij
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

      ij=0 
c 
ctj   First we need to have the coulomb integrals on the disk
ctj   Reserve space for coulomb and exchange integrals indexer
      call getint(npairs,ilist)
      call filllist(bl(ilist),npairs)
      call dynamic_getmem(npairs/intsize+1,irecadrc)
      call izeroit(bl(irecadrc),npairs)         !zero out the pair records
      call dynamic_getmem(npairs/intsize+1,irecadrx)
      call izeroit(bl(irecadrx),npairs)         !zero out the pair records
      call dynamic_getmem(npairs/intsize+1,irecadre)
      call izeroit(bl(irecadre),npairs)         !zero out the pair records
      call dynamic_getmem(nval*nval/intsize+1,irecadrtx)
      call izeroit(bl(irecadrtx),nval*nval)     !zero out the pair records
      call dynamic_getmem(nval*nval/intsize+1,irecadrtc)
      call izeroit(bl(irecadrtc),nval*nval)     !zero out the pair records
      call dynamic_getmem(nval*(nbf-nmo)/intsize+1,irecadrtt)
      call izeroit(bl(irecadrtt),nval*(nbf-nmo))   !zero out the pair records
ctj Generate the coulomb integrals
      int_kind='c'
      ipass=0
      call elapsec(ec_start)
      call secund(sc_start)
      if (.not.mp2) then
         call GenCoulExInt(ncs,    ncf,         ictr,   nval,  nmo,
     *                     nfirst, nlast,       thresh, core,  xmaxdisk,
     *                     ndiskc, bl(irecadrc),lrecc,  lbinc, int_kind,
     *                     nodisk, byt8,        small,  vorb,  ndisktrc,
     *                     bl(ilist),  nstrong,bl(isympairs), ij_unique,
     *                     nbf,iprnt)
      endif
      call secund(sc_stop)
      call elapsec(ec_stop)
      if (iprnt.ge.3) then
         write(6,*) 'Elapsec time for coulomb integral calc. (sec): ',
     *              ec_stop-ec_start
         write(6,*) 'CPU time for coulomb integral calc. (sec):     ',
     *              sc_stop-sc_start
         write(6,*) 'Ipass:  ',ipass
      endif
      ipass=0
ctj Generate the exchange integrals
      int_kind='x'
      call elapsec(ex_start)
      call secund(sx_start)
      call GenCoulExInt(ncs,    ncf,          ictr,   nval,  nmo,
     *                  nfirst, nlast,        thresh, core,  xmaxdisk,
     *                  ndiskx, bl(irecadrx), lrecx,  lbinx, int_kind,
     *                  nodisk, byt8,         small,  vorb,  ndisktrx,
     *                  bl(ilist),   nstrong, bl(isympairs), ij_unique,
     *                  nbf,iprnt)
      call secund(sx_stop)
      call elapsec(ex_stop)
      if (iprnt.ge.3) then
         write(6,*) 'Elapsec time for exchange integral calc. (sec): ',
     *              ex_stop-ex_start
         write(6,*) 'CPU time for exchange integral calc. (sec):     ',
     *              sx_stop-sx_start
         write(6,*) 'Ipass:  ',ipass
      endif
ctj Generate the exchange integrals
      if (l_triples) then
         int_kind='tt'
         call elapsec(ex_start)
         call secund(sx_start)
         call GenCoulExInt(ncs,    ncf,        ictr,   nval,  nmo,
     *                     nfirst, nlast,      thresh, core,  xmaxdisk,
     *                     ndisktt,bl(irecadrtt),lrectt,lbintt,int_kind,
     *                     nodisk, byt8,       small,  vorb,  ndisktrtt,
     *                     bl(ilist), nstrong, bl(isympairs), ij_unique,
     *                     nbf,iprnt)
c ndisktt,irecadrtt,lrectt,lbintt,ndisktrtt
         call secund(sx_stop)
         call elapsec(ex_stop)
         if (iprnt.ge.3) then
            write(6,*) 'Elapsec time for 3ext integral calc. (sec): ',
     *                  ex_stop-ex_start
            write(6,*) 'CPU time for 3ext integral calc. (sec):     ',
     *                  sx_stop-sx_start
            write(6,*) 'Ipass:  ',ipass
         endif
      endif
c 1
 
      call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
     *              vorb,  nmo, ivalpair,     ivirpair,   ifp, 
     *              nsym,  nval,nvirt,        ij_unique,  isympairs,
     *              nbf)
      call matconn('fockAO','q',idimen,idimen,ifockAO)
      call dynamic_getmem(ncf*ncf,ifockMO)
      call dynamic_matdef('corefock','q',idimen,idimen)
      if (mp2) lbinc=0
      icorefock=mataddr('corefock')
      call FockBuilder(bl(irecadrc),bl(irecadrx),nval,npairs,thresh,
     &                 ndiskc,ndiskx,lbinc,lbinx,byt8,nmo,vorb,ifockAO,
     &                 ifockMO,icorefock,nbf)

C For CIM local calculation, we cannot get the diagnal elements of Fock
C matrix from previous calculation. So we need to get it from the Fock
C matrix by hand. --NZG_10/23/2016 @NJU
      if (loca) call getepsi(bl(ifockMO),bl(iepsi),ncf)

      call CreateKijklDisk(nval,kijklndisk,af)
      call KijklInit(bl(irecadrx),npairs,ndiskx,lbinx,ncf,thresh,
     *        nfirst,nlast,bl(ifockMO),byt8,nmo,vorb,kijklndisk,af)
      call Kijk_Vec_Init(bl(irecadrx),npairs,ndiskx,lbinx,ncf,thresh,
     *                 nfirst,nlast,byt8,nmo,vorb,.false.,kijkndisk,nbf)
      call Fock_Vector_Init(ncf,nfirst,nlast,ifockAO,ifockMO,nmo,vorb,
     *                      nbf)
      call initsingles(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call initS(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
      call initR(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
c
c  Memory for energy pairs:
      call dynamic_getmem(npairs,iepair)
      call dynamic_getmem(npairs,iepairr)
      call dynamic_getmem(npairs,ixepair)
      call dynamic_getmem(npairs,ixepairr)
c
c Space for one K_ij matrix
      call dynamic_matdef('resao','q',idimen,idimen)
      iresAO=mataddr('resao')
      call dynamic_matdef('exchao','q',idimen,idimen)
      iexchAO=mataddr('exchao')
      call dynamic_matdef('invcano','q',ncf,ncf)
      iinvcano=mataddr('invcano')
      call matmmul2('cano','overlap','invcano','t','n','n')
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
c
      square=0.0d0
      energy=0.0d0
      scsenergy=0d0
      scsenergy1=0d0
      scsenergy2=0d0
      energyr=0.0d0
      ij=0
      call matzero('resao')
      call elapsec(e1_start)
      call secund(s1_start)

C for DIIS
      call prepare_diis_files(ndiskdiisr,ndiskdiisc,af,ncf,nmo,vorb,nbf)
      i=length_diis
      call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx', bl(ilist),
     *            bl(ipairimages),nbf)
c END DIIS

      xmaxele=0.0d0
      do i=1,nval
         do j=1,i
            ij=ij+1
            equal=.false.
            if (i.eq.j) equal=.true.
            call ExtrOne(i,j,     bl(irecadrx), npairs,       ndiskx,
     *                   ncf,  lbinx, thresh,       byt8,         'x',
     *                   'mo', nmo,   vorb,         bl(iexchAO))
            call CoefRead('tt',i,j,icoeAO)
            call UpdateCoef(i,    j,      icoeAO,  iexchAO,   bl(iepsi),
     *                      ncf,  nmo,    nfirst,  npairs,    shiftini,
     *                      'nodiis',vorb,singles,ndiskdiisr,ndiskdiisc,
     *                      sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                      bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                      docansym,nbf)
            call CoefWrite(i,j,icoeAO)
            call PairEnergyMP2(icoeAO, iexchAO, iresAO, epairr, epair,
     *                         ncf,    equal,   nmo,    vorb,   sepairr,
     *                         sepair,nbf)
            energy=energy+epair
            scsenergy=scsenergy+sepair
            scsenergy1=scsenergy1+sepairr
            scsenergy2=scsenergy2+sepair
            square=square+sqij
         enddo
      enddo
c
c     Singles definitions
c
      call dynamic_matdef('sing_res','r',idimen,1)
      ising_res=mataddr('sing_res')
      call dynamic_matdef('work_v1','r',idimen,1)
      iwork_v1=mataddr('work_v1')
      call dynamic_matdef('singl_am','r',idimen,1)
      isingl_am=mataddr('singl_am')
c
      call secund(s1_stop)
      call elapsec(e1_stop)
      if (iprnt.ge.3) then
         write(6,*) 'Elapsec time for first iteration:  ',
     *   e1_stop-e1_start
         write(6,*) 'CPU time for first iteration:      ',
     *   s1_stop-s1_start
         write(6,*)  'Sum of the squares of the residuums: ',square
         write(6,*)  'Maximum residuum element (abs value) ',xmaxele
         write(6,'(A,F21.17)')  'Nonlocal! - "MP2" energy: ',energy
      endif
      if (.not.loca) then
         xmp2=energy
         scsmp2=scsenergy2
         if (mp2) goto 999
      else
         write(6,'(A,F21.17)')  'For local - MP2 energy: ',energy
      endif
      call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
     *              vorb,  nmo, ivalpair,     ivirpair,   ifp, 
     *              nsym,  nval,nvirt,        ij_unique,  isympairs,
     *              nbf)
c
      mp2_iter=0
      mp2last=.false.
      write(6,'(A)') "     MP2 module:"
      write(6,'(A)')
      write(6,'(A)')
     * "Iter:  Energy:          Delta E         "//
     * "Max resid:    Err sq:      Elapsed time"
 222  continue
c
      ij=0
      mp2_iter=mp2_iter+1
      energy=0.0d0
      if (mp2_iter>1) energyrold=energyr
      energyr=0.0d0
      energy1=0.0d0
      energy2=0.0d0
      square=0.0d0
      xmaxele=0.0d0
      scsenergy1=0d0
      scsenergy2=0d0
      call secund(t0)
      call elapsec(et0)
      ndiskmp2=69
      idimen1=idimen
      call Open_Scratch_File(ndiskmp2,idimen,idimen1,.false.,'lmp2')
      call secund(txx0)
      call Generate_MP2_G(nval,ncf,ndiskmp2,nmo,nfirst,vorb,.false.,
     *                    bl(ifockMO),nbf)
      call secund(txx1)
      reading=0.0d0
      xmultiplication=0.0d0
      update_energy=0.0d0
      update=0.0d0
      copy=0.0d0

C For CIM calculation, we need to collect the energy contribution from
C orbitals to get the energy of a certain subsystem
C NZG_10/24/2016 @NJU
      if (mp2last) allocate(CIMPairEnergy(nval,nval),CIMOrbEnergy(nval))
      do i=1,nval
         do j=1,i
            ij=ij+1
            equal=.false.
            if (i.eq.j) equal=.true.
            if (int_array(bl(isympairs),ij).ne.0) then
               call secund(tx0)
               if (do_mp4) then
                  call reader(ndisk_mp3,idimen,ij,bl(iexchAO))
               else
                  call ExtrOne(i,  j,   bl(irecadrx), npairs,    ndiskx,
     *                         ncf,  lbinx, thresh,   byt8,        'x',
     *                         'mo', nmo,   vorb,         bl(iexchAO))
               endif
               call secund(tx1)
               reading=reading+(tx1-tx0)
               call secund(tx0)
               call matcopy('exchao','resao')
               call secund(tx1)
               copy=copy+(tx1-tx0)
               call secund(tx0)
               call CoefRead('tt',i,j,icoeffAO)
               call secund(tx1)
               reading=reading+(tx1-tx0)
               call matconn('coefAO','q',idimen,idimen,icoeffAO)
               call secund(tx0)
               call matmmul2('fockAO','coefAO','resao','n','n','a')
               call matmmul2('coefAO','fockAO','resao','n','n','a')
               call secund(tx1)
               xmultiplication=xmultiplication+(tx1-tx0)
               call matdisc('coefAO')
               call secund(tx0)
               call Read_MP2(i,j,nmo,ncf,nval,ndiskmp2,vorb,.false.,
     &                       iwork1,nbf)
               call secund(tx1)
               reading=reading+(tx1-tx0)
               call secund(tx0)
               call matadd('work1','resao')
               call secund(tx1)
               copy=copy+(tx1-tx0)
               call secund(tx0)
               call CoefRead('tt',i,j,icoeffAO)
               call secund(tx1)
               reading=reading+(tx1-tx0)
               call secund(tx0)
               call PairEnergyMP2(icoeffAO,iexchAO,iresAO,
     &                            bl(iepairr+ij-1),bl(iepair+ij-1),
     *                            ncf,equal,nmo,vorb,sepairr,sepair,nbf)
               call UpdateCoef(i, j,      icoeffAO,iresAO,    bl(iepsi),
     *                      ncf,    nmo,   nfirst,  npairs,    shift,
     *                      'nodiis',vorb,singles,ndiskdiisr,ndiskdiisc,
     *                      sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                      bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                      docansym,nbf)
               square=square+sqij
               call secund(tx1)
               update_energy=update_energy+(tx1-tx0)
               call secund(tx0)
               call CoefWrite(i,j,icoeffAO)
               call secund(tx1)
               reading=reading+(tx1-tx0)
               call secund(tx0)
               call PairEnergyMP2(icoeffAO,iexchAO,iresAO,
     &                            bl(ixepairr+ij-1),bl(ixepair+ij-1),ncf
     *                            ,equal,nmo,vorb,sxepairr,sxepair,nbf)
               call secund(tx1)
               update_energy=update_energy+(tx1-tx0)
            else !  if (int_array(bl(isympairs),ij).ne.0) then
               bl(iepair  +ij-1)=0.0d0
               bl(iepairr +ij-1)=0.0d0
               bl(ixepair +ij-1)=0.0d0
               bl(ixepairr+ij-1)=0.0d0
            endif            !   if (int_array(bl(isympairs),ij).ne.0)
            if ((mp2_iter.eq.1.or.mp2last).and.iprnt.ge.3) then
               if (do_mp4) then
                  write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',i,j,
     *                                       bl(ixepair+ij-1)
               else
                  write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',i,j,
     *                                       bl(iepairr+ij-1)
               endif
            endif
            energyr=energyr+bl(iepairr+ij-1)
     *               *dble(int_array(bl(ipairimages),ij))
            energy1=energy1+bl(iepair+ij-1)
     *               *dble(int_array(bl(ipairimages),ij))
            energy2=energy2+bl(ixepair+ij-1)
     *               *dble(int_array(bl(ipairimages),ij))
            scsenergy1=scsenergy1+sepairr
     *               *dble(int_array(bl(ipairimages),ij))
            scsenergy2=scsenergy2+sxepair
     *               *dble(int_array(bl(ipairimages),ij))
            if (mp2last) then
               pairtmp=bl(iepairr+ij-1)
     &                    *dble(int_array(bl(ipairimages),ij))
               if (.not.equal) then
                  CIMPairEnergy(i,j)=pairtmp/2.0D0
                  CIMPairEnergy(j,i)=pairtmp/2.0D0
               else
                  CIMPairEnergy(i,j)=pairtmp
               endif
            endif
         enddo
      enddo
      if (iprnt.ge.3) then
         write(6,'(A30,F10.2)') 'Time for Generate_MP2_G: ',txx1-txx0
         write(6,'(A30,F10.2)') 'Reading in MP2: ',reading
         write(6,'(A30,F10.2)') 'Multiplication in MP2:',xmultiplication
         write(6,'(A30,F10.2)') 'Update+energy in MP2: ',update_energy
         write(6,'(A30,F10.2)') 'Copy:                 ',copy
         write(6,'(A30,F10.2)') 'Pure update:          ',update
      endif
      close(ndiskmp2,STATUS='DELETE')
      call secund(t1)
      call elapsec(et1)
c
      call secund(tx0)
      call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
     *              vorb,  nmo, ivalpair,     ivirpair,   ifp, 
     *              nsym,  nval,nvirt,        ij_unique,  isympairs,
     *              nbf)
      call secund(tx1)
      if (iprnt.ge.2) then
         write(6,'(A30,F10.2)') 'CoefInit              ',tx1-tx0
         write(6,'(A,I3)')   'MP2 iteration no:              ', mp2_iter
      endif
      if (iprnt.ge.2) then
         if (.not.do_mp4) then
            write(6,'(A,F21.17)')  'Iterations - MP2 linear prev:  ',
     &                             energy1
            write(6,'(A,F21.17)')  'Iterations - MP2 energy+resid: ',
     &                             energyr
            write(6,'(A,F21.17)')  'Iterations - MP2 linear curr:  ',
     &                             energy2
         else
            write(6,'(A,F21.17)')  'Iterations - MP4 energy:       ',
     &                             energy2
         endif
         write(6,*)  ' '
         write(6,'(A,F6.2)')  'Total time for iteration CPU:  ',t1-t0
         write(6,'(A,F6.2)')  'Total time for iteration ELA:  ',et1-et0
         write(6,'(A,F20.16)')  'Sum of squares of the residuum:',square
         call flush(6)
      endif
      if (do_mp4) then
         xmp4_d = energy2
      else
         xmp2   = energyr
         scsmp2 = scsenergy2
      endif
      xxenergy=energyr
      if (do_mp4) xxenergy=energy2
      xdelta=energyrold-energyr
      write(6,3232) mp2_iter,xxenergy,xdelta,xmaxele,sqij,(et1-et0)/6d1
      call flush(6)
 3232 format (I3,2F16.10,2E14.3E2,F15.1)
      if (dabs(energy2-energy1).gt.Ethresh.or.sqij.gt.Wavethresh) then
         goto 222
      else
         if (.not.mp2last) then
            mp2last=.true.
            goto 222
         endif
      endif
      if (do_mp4) goto 999
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
         goto 999
      endif
c     END loop
c if localization is on, calculate list of strong pairs
      call CalcList(bl(ilist),nfirst,nval,nmo,omit,inumber,
     *              weaksepar,distsepar,bl(iepairr),iEnergy)
      write(6,*)  ''
      write(6,*)  ''
      write(6,*)  'There are ',inumber,' strong pairs.'
c
      call dynamic_matdef('WernerX','q',idimen,idimen)
      iX=mataddr('WernerX')
c
      xiterenergy=energyr
      if (mp3) xiterenergy=0.0d0
      ee_tot=0.0d0
      se_tot=0.0d0
      e_ext_mat_tot=0.0d0
      s_ext_mat_tot=0.0d0
      eq_tot=0.0d0
      sq_tot=0.0d0
      eg_tot=0.0d0
      sg_tot=0.0d0
      een_tot=0.0d0
      sen_tot=0.0d0
      eco_tot=0.0d0
      sco_tot=0.0d0
      iteration=0
      totstrace=0.0d0
      totsconstr=0.0d0
      totsort=0.0d0
      totextr=0.0d0
      totconst=0.0d0
      varenergyold=0.0d0
      diis=0
c Calculate MP4 triples:
      if (mp4.and.l_triples) then
ctmp  call triples(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),res,
ctmp *             ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
ctmp *             singles.and.cc.and..not.ccsd,af,1)
ctmp  write(6,'(A,F20.15)') 'MP4 triples correction: ',res
         call elapsec(xt0)
         call sort_amplit(nval,idimen,ndisk_a,af,0,ichar_count,
     &                    bl(ivtable),bl(im_table),bl(ijtable),
     *                    bl(ivrevtable),c1_ratio)
c     call sort_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,ndisktt,
c    *                     ncf,lbintt,thresh,byt8,nmo,vorb,af,0)
         call simple_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                    ndisktt,ncf,lbintt,thresh,byt8,
     *                    nmo,vorb,af,0,ichar_count,bl(ivtable),
     *                    bl(im_table),bl(ivrevtable),bl(iotable))
         call make_3ext_pairs(idimen,nval,af,ndisk_ie,ndisk_ie1)
         call sort_3int(nval,idimen,ndisk_ii,af,0)
         amp_ratio=1d0
         ext3_ratio=1d0
         nslv=1
         call dynamic_show_free(mem)
         mem=mem-nval*nval*nval-3*nval*nval
         call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                        nslv,npass,isize)
         call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *                  ncf,lbinx,thresh,byt8,nmo,vorb,af,0,isize,npass)
c     call check_kext(nval,idimen,ndisk_ix)
c     call check_kext1(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
c    *                 ncf,lbinx,thresh,byt8,nmo,vorb)
         call elapsec(xt1)
         if (iprnt.ge.4) write(6,*) 'Sort time: ',xt1-xt0
         call flush(6)
         call elapsec(xt0)
         call dynamic_getmem(20,itimes)
c     call dynamic_getmem(nval*nval*idimen*idimen,ixt)
c     call check_kext1(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
c    *                 ncf,lbinx,thresh,byt8,nmo,vorb,bl(ixt))
         call new_triples(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,
     &                    ndisk_ix,res,bl(iepsi+ncore),bl(itimes),20,
     &                    ccsd,qcisd,af,bl(iotable),bl(ivtable),
     &                    bl(ivrevtable),bl(im_table),bl(ijtable),
     &                    ichar_count,npass,isize,energys,energyd,iprnt)
         call elapsec(xt1)
         if (iprnt.ge.4) then
            write(6,*) res,'   Time: ',xt1-xt0
            write(6,*) 'W zero & sort:           ',bl(itimes+0)
            write(6,*) 'W build (total):         ',bl(itimes+1)
            write(6,*) 'Mult. over virt space:   ',bl(itimes+2)
            write(6,*) 'Mult. over occ. space:   ',bl(itimes+3)
            write(6,*) 'Reading of ampl. & int:  ',bl(itimes+4)
            write(6,*) 'Amplitudes:              ',bl(itimes+5)
            write(6,*) '3int integrals:          ',bl(itimes+6)
            write(6,*) '3ext integrals:          ',bl(itimes+7)
            write(6,*) 'Kext integrals:          ',bl(itimes+9)
            write(6,*) 'Energy:                  ',bl(itimes+8)
            write(6,*) 'Wabc relocate:           ',bl(itimes+10)
            write(6,*) 'Inte relocate:           ',bl(itimes+11)
            write(6,*) '3int relocate:           ',bl(itimes+12)
            write(6,*) 'Total virt part:         ',bl(itimes+13)
            write(6,*) 'Tables reading:          ',bl(itimes+14)
            call flush(6)
         endif
      endif
      write(6,'(A)')
      write(6,'(A)')
      write(6,'(A)') "     CC module:"
      write(6,'(A)')
      write(6,'(A)')
     * "Iter:  Energy:        Delta E           "//
     * "Max resid:    Err sq:      Elapsed time"

      converged=.false.
 666  continue  !CC iteration
      call elapsec(eistart)
      call secund(sistart)
      ccsdcpu=0.0d0; ccsdela=0.0d0
c
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
c
      if (iprnt.ge.3) then
         write(icond,"(' The CCSD iteration:',I4)") iteration
      endif
      xmaxele=0.0d0
      iteration=iteration+1
      if (cep0) xiterenergy=0.0d0
      call elapsec(eistart)
      call secund(sistart)
c
      int_kind='e' ! EE operator
c     call izeroit(bl(irecadre),npairs)         !zero out the pair records
      strace=0.0d0
      sconstr=0.0d0
      sort=0.0d0
      extr=0.0d0
      const=0.0d0
      call secund(t0)
      call elapsec(et0)
      if (.not.do_mp4) then
        call NEW_EEO_INT(ncs,ncf,bl(ictr),nval,nmo,nfirst,nlast,thresh,
     *                   vorb,ndisktre,npairs,af,iprnt,nbf)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(1)=t1-t0        ! EEO
      ccsdela(1)=et1-et0      ! EEO
c
      totstrace  = totstrace + strace
      totsconstr = totsconstr+ sconstr
      totsort    = totsort   + sort
      totextr    = totextr   + extr
      totconst   = totconst  + const
      ee_tot=ee_tot+ee_stop-ee_start
      se_tot=se_tot+se_stop-se_start
c
      if (ccsd) then
         call secund(t0)
         call elapsec(et0)
         int_kind='tx' ! TEIO operator
c     call izeroit(bl(irecadrtx),nval*nval)     !zero out the pair records
         call GenCoulExInt(ncs,    ncf,         ictr,  nval,  nmo,
     *                     nfirst, nlast,       thresh,core,  xmaxdisk,
     *                     ndisktx,bl(irecadrtx),lrectx,lbintx,int_kind,
     *                     nodisk, byt8,        small, vorb,  ndisktrtx,
     *                     bl(ilist),  nstrong, bl(isympairs),ij_unique,
     *                     nbf,iprnt)
         call secund(t1)
         call elapsec(et1)
         ccsdcpu(2)=t1-t0        ! TEIO tx
         ccsdela(2)=et1-et0      ! TEIO tx
c
         call secund(t0)
         call elapsec(et0)
         int_kind='tc' ! TEIO operator
c     call izeroit(bl(irecadrtc),nval*nval)     !zero out the pair records
         call GenCoulExInt(ncs,    ncf,          ictr,  nval,  nmo,
     *                     nfirst, nlast,        thresh,core,  xmaxdisk,
     *                     ndisktc,bl(irecadrtc),lrectc,lbintc,int_kind,
     *                     nodisk, byt8,         small, vorb, ndisktrtc,
     *                     bl(ilist),nstrong,   bl(isympairs),ij_unique,
     *                     nbf,iprnt)
         call secund(t1)
         call elapsec(et1)
         ccsdcpu(3)=t1-t0        ! TEIO tc
         ccsdela(3)=et1-et0      ! TEIO tc
      endif
c
      call secund(t0)
      call elapsec(et0)
      elaps(1:5)=0.0d0
      call CreateAlphaDisk(nval,ndiskalpha,af)
      if (cc.or.ccsd.or.do_mp4) then
         call CCalphaonce(ncf,nval,bl(irecadrx),npairs,ndiskx,lbinx,
     *                    thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                    ndiskalpha,af,nbf)
         if (iprnt.ge.3) then
            write(6,'(A)') 'CCalph: ampl:       Kij:     mult:  putres:'
            write(6,'(4X,4F10.2)') elaps(1),elaps(2),elaps(3),elaps(4)
         endif
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(4)=t1-t0        ! CCalphaonce+supplement
      ccsdela(4)=et1-et0      ! CCalphaonce+supplement
c
      call secund(t0)
      call elapsec(et0)
      call secund(t1)
      call elapsec(et1)
      if (iprnt.ge.3) then
         write(6,*) 'SupplementAlpha_BetaBuild, ELAPS: ',et1-et0
         write(6,*) 'SupplementAlpha_BetaBuild, CPU:   ',t1-t0
      endif
c
      call secund(t0)
      call elapsec(et0)
      call Builder41a_G(ncf,nval,norecalc,cc,ccsd,nfirst,bl(ifockMO),
     *               nmo,vorb,ndiskG41a,.false.,ndiskalpha,bl(ibeta),
     *                 nbf)
      close(ndiskalpha,STATUS='delete')
      call secund(t1)
      call elapsec(et1)
c
      ccsdcpu(5)=t1-t0        ! Builder41a_G
      ccsdela(5)=et1-et0      ! Builder41a_G
c 
      call secund(t0)
      call elapsec(et0)
      if (cc.or.ccsd.or.do_mp4) then
         call CCAOnce(bl(irecadrx),npairs,nval,ndiskx,ncf,lbinx,thresh,
     *                norecalc,byt8,ioverlap,nmo,vorb,nbf)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(6)=t1-t0        ! CCAOnce
      ccsdela(6)=et1-et0      ! CCAOnce
c 
      ij=0
      xlinear_energy=0.0d0
c
c This is the s and r (singles part) production code, START
      call secund(t0)
      call elapsec(et0)
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
         call Tl_generator(ncf,nval,nmo,vorb,ioverlap,iiS,nbf)
         call elapsec(esi1)
         call secund(csi1)
         ccsdela(15)=esi1-esi0 ! Tl_generator
         ccsdcpu(15)=csi1-csi0 ! Tl_generator
         call elapsec(esi0)
         call secund(csi0)
         call EEO_vector_extractor(nval,bl(irecadre),npairs,ndiske,ncf,
     *                             lbine,thresh,nfirst,nlast,byt8,
     *                             nmo,vorb,iiS,nbf)
         call elapsec(esi1)
         call secund(csi1)
         ccsdela(16)=esi1-esi0 ! EEO_vector_extractor
         ccsdcpu(16)=csi1-csi0 ! EEO_vector_extractor
c  R:
         call Fock_Vector_Pointer(ifock)
         call pointRR(iiR)
         call tfer(bl(ifock),bl(iiR),nval*idimen)
         if (cc.or.ccsd) then
            call Lt_generator(idimen,nval,   bl(irecadrx),npairs,ndiskx,
     *                        ncf,   lbinx,  thresh,      byt8,  nmo,
     *                        vorb,  iiR)
         endif
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(8)=t1-t0        ! singles s, r
      ccsdela(8)=et1-et0      ! singles s, r
c
c This was the s,r code STOP
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
c
      energy=0.0d0
      energyr=0.0d0
 133  FORMAT (4I3,F20.10)
      square=0.0d0
      totnorm=0.0d0
c
      call secund(t0)
      call elapsec(et0)
      if (calchunks)
     *    call calculate_ij_chunk(ichunk,jchunk,1,idimen,nval)
c
      ipassess=nval/ichunk
      if (mod(nval,ichunk).gt.0) ipassess=ipassess+1
      jpassess=nval/jchunk
      if (mod(nval,jchunk).gt.0) jpassess=jpassess+1
c
      call prepare_CCYZ_file(ncf,nmo,vorb,ndiskYZ,.false.,nbf)
      do ipass=1,ipassess
         istart=(ipass-1)*ichunk+1
         istop = ipass   *ichunk
         if (istop.gt.nval) istop=nval
         do jpass=1,jpassess
            jstart=(jpass-1)*jchunk+1
            jstop = jpass   *jchunk
            if (jstop.gt.nval) jstop=nval
            call CCYZ(istart,istop,jstart,jstop,ncf,nval,bl(irecadrx),
     *                bl(irecadrc),npairs,ndiskx,ndiskc,lbinx,lbinc,
     &                bl(irecadrtx),ndisktx,lbintx,bl(irecadrtc),
     &                ndisktc,lbintc,thresh,byt8,ioverlap,cc,ccsd,
     &                bl(ilist),nmo,vorb,ndiskYZ,.false.,elaps,do_mp4,
     *                ichar_count,bl(iorevtable),bl(ivtable),
     &                bl(im_table),nbf)
         enddo
      enddo
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(9)=t1-t0        ! YZ
      ccsdela(9)=et1-et0      ! YZ
c      
      call secund(t0)
      call elapsec(et0)
      call prepare_Qparts_file(ncf,nmo,vorb,ndiskQ,af,nbf)
      do ipass=1,ipassess
         istart=(ipass-1)*ichunk+1
         istop = ipass   *ichunk
         if (istop.gt.nval) istop=nval
         do jpass=1,jpassess
            jstart=(jpass-1)*jchunk+1
            jstop = jpass   *jchunk
            if (jstop.gt.nval) jstop=nval
            if (iprnt.ge.3)
     *                   write(6,"(' istart,istop,jstart,jstop: ',4I4)")
     &                   istart,istop,jstart,jstop
            call Qparts(istart,istop,jstart,jstop,ncf,nval,npairs,byt8,
     &                  cc,ccsd,bl(ilist),nmo,vorb,ndiskYZ,af,ndiskQ,
     *                  ichar_count,bl(iorevtable),bl(ivtable),
     &                  bl(im_table),nbf)
         enddo
      enddo
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(10)=t1-t0        ! Q parts
      ccsdela(10)=et1-et0      ! Q parts
c
      call secund(t0)
      call elapsec(et0)
c
      if (converged) then
         allocate(exc_QCMO(nval,nval,idimen,idimen))
         allocate(coef_QCMO(nval,nval,idimen,idimen))
         allocate(resi_QCMO(nval,nval,idimen,idimen))
         exc_QCMO=0.0D0; coef_QCMO=0.0D0; resi_QCMO=0.0D0
      endif

      if (.not. allocated(EPairQCMO)) allocate(EPairQCMO(nval,nval))

      do i=1,nval
         do j=1,i
            ij=ij+1
            call flush(6)
            equal=.false.
            if (i.eq.j) equal=.true.
            if (int_array(bl(isympairs),ij).ne.0) then
               if (cep2) xiterenergy=bl(iepairr+ij-1)
c
c  Read exchange operators directly to residuum matrix
               if (.not.do_mp4) then
                  call ExtrOne(i,j,bl(irecadrx),npairs,ndiskx,ncf,lbinx,
     &                        thresh,byt8,'x','mo',nmo,vorb,bl(iexchAO))
c  Read external echange operators
                  call ExtrOne(i,j,bl(irecadre),npairs,ndiske,ncf,lbine,
     &                         thresh,byt8,'e','mo',nmo,vorb,bl(iwork1))
c  Add it do residuum matrix
                  call matcopy('exchao','resao')
                  call matadd('work1','resao')
               else
                  call matzero('resao')
               endif
c
c  Qgen, iwork1=Qij, iwork2=Qji(T), probably OK, order of en. OK
               call Q_read_build(i,j,iX,ndiskQ,af,vorb,ncf,nmo,nval,
     *                           iwork1,iwork2,nbf)
c
c   S*Qij and add to residuum
c
               call matadd('work1','resao')
               call matadd('work2','resao')
c   Qji(T)*S and add to residuum
c
c   Read Gij and Gji(T) with part of CCSD
c
               call ReadCalcG_41a(i,j,norecalc,ncf,nmo,vorb,ndiskG41a,
     *                            iwork1,.false.,nbf)
c
               if (cc.or.ccsd.or.do_mp4) then
                  continue
               else
                  call CoefRead('tt',i,j,icoeffAO)
                  call matconn('coeffAO','q',idimen,idimen,icoeffAO)
                  call matadd1('coeffAO',-xiterenergy,'work1')
                  call matdisc('coeffAO')
               endif
c       S*Sum*S
               call matadd('work1','resao')
c
c Singles part of CISD CCSD start
               if (singles.or.ccsd) then
                  call dynamic_matdef('CISD','q',idimen,idimen)
                  iresult=mataddr('CISD')
                  call singles_CISD(i,j,ncf,nval,ioverlap,bl(irecadre),
     &                         npairs,nfirst,nlast,ndiske,lbine,thresh,
     &                         byt8,ccsd,iresult,nmo,vorb,nbf)
                  call matadd('CISD','resao')
                  call dynamic_matrem('CISD')
               endif
c Singles part of CISD CCSD stop

               call CoefRead('tc',i,j,icoeffAO)
               if (mp3.or.do_mp4) then
                  call PairEnergy(icoeffAO,iresAO,iresAO,
     *                            bl(iepairr+ij-1),bl(iepair+ij-1),ncf,
     &                            equal,nmo,vorb,nbf)
                  if (mp3.and.mp4) call writer(ndisk_mp3,idimen,ij,
     &                                         bl(iresAO))
               else
                  if (converged) then
                     call IntExpress(i,j,bl(iexchAO),bl(icoeffAO),
     &                               bl(iresAO),exc_QCMO,coef_QCMO,
     &                               resi_QCMO,idimen,nval)
                     cycle
                  endif
                  call PairEnergy(icoeffAO,iexchAO,iresAO,
     *                            bl(iepairr+ij-1),bl(iepair+ij-1),ncf,
     &                            equal,nmo,vorb,nbf)
                  if (.not.equal) then
                     EPairQCMO(i,j)=bl(iepairr+ij-1)/2.0D0
                     EPairQCMO(j,i)=bl(iepairr+ij-1)/2.0D0
                  else
                     EPairQCMO(i,j)=bl(iepairr+ij-1)
                  endif
               endif
c
               call CoefRead('tt',i,j,icoeffAO)
               if (.not.(mp3.or.do_mp4)) then
                  call CINorm(icoeffAO,xnormij,ncf,equal,nmo,vorb,nbf)
                  totnorm=totnorm+xnormij
                  call UpdateCoef(i,   j, icoeffAO, iresAO,   bl(iepsi),
     *                   ncf,    nmo,  nfirst,   npairs,     shift,
     *                   'diis', vorb, singles,  ndiskdiisr, ndiskdiisc,
     *                   sqij,   xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
               endif
c
               call CoefWrite(i,j,icoeffAO)
               call CoefRead('ee',i,j,icoeffAO)
               call PairEnergy(icoeffAO,iexchAO,iresAO,
     *                         bl(ixepairr+ij-1),bl(ixepair+ij-1),ncf,
     &                         equal,nmo,vorb,nbf)
c
 777           continue  ! Jump to this point if the pair has to be omitted
            else ! if (int_array(bl(isympairs),ij).ne.0) then
               bl(iepair  +ij-1)=0.0d0
               bl(iepairr +ij-1)=0.0d0
               bl(ixepair +ij-1)=0.0d0
               bl(ixepairr+ij-1)=0.0d0
               sqij=0.0d0
            endif ! if (int_array(bl(isympairs),ij).ne.0) then
            xlinear_energy=xlinear_energy+bl(ixepair+ij-1)
     *                           *dble(int_array(bl(ipairimages),ij))
            energyr=energyr+bl(iepairr+ij-1)
     *                           *dble(int_array(bl(ipairimages),ij))
c  Calculate quadratic energy here:
            energy=energy+bl(iepair+ij-1)
     *                           *dble(int_array(bl(ipairimages),ij))
            square=square+sqij*dble(int_array(bl(ipairimages),ij))
c   Print Quadratic pair energies:
            if (iprnt.ge.3) then
               if (.not.mp3) then
                  write(6,'(A,2I5,F20.12)') 'Pair en:  i,j,energy:',i,j,
     *                               bl(iepairr+ij-1)
               else
                  write(6,'(A,2I5,F20.12)') 'Pair en:  i,j,energy:',i,j,
     *                               bl(iepair+ij-1)
               endif
            endif
            call flush(6)
         enddo
      enddo

      if (converged) then
         if (iprnt>=2) then
            write(6,*) "Pair correlation energies between QCMOs:"
            call NJ_prtcol(6,nval,EPAIRQCMO,1,nval,'f14.10')
            write(6,*)
         endif
         goto 888
      endif

      call secund(t1)
      call elapsec(et1)
      ccsdcpu(11)=t1-t0        ! Residuals
      ccsdela(11)=et1-et0      ! Residuals

c
c calculate variational (I hope) CI energy:
      varenergy=(energyr+totnorm*xiterenergy)/(1.0d0+totnorm)
c
c This is the singles residuum loop, START
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

      sing: if (singles.or.ccsd.or.do_mp4) then
         ssq=0.0d0
         call pointS(iiS)
         call dynamic_matdef('s_resid','r',idimen,nval)
         is_resid=mataddr('s_resid')
         call pointersingles(isingles)
         call pointernewsingles(inewsingles)
         call matconn('full_sing','r',idimen,nval,isingles)
         call tfer(bl(iiS),bl(is_resid),idimen*nval)
         call Tr_generator1(idimen,nval,ncf,vorb,ioverlap,is_resid)
         call beta_t_sumator1(ncf,idimen,nval,ibeta,vorb,ioverlap,
     *                        is_resid)
         if ((.not.cc).and.(.not.ccsd).and.(.not.do_mp4)) then
            call cisd_energy_add(ncf,idimen,nval,xiterenergy,vorb,
     *                       ioverlap,is_resid)
         endif
         if (do_mp4) then
            call matmmul2('fockAO','full_sing','s_resid','n','n','a')
         endif
         call matdisc('full_sing')
         call UpdateSingles1(isingles,is_resid,nmo,nfirst,
     *                       ncf,npairs,nval,bl(iepsi),shift,ssq,vorb,
     *                       sing_max,ndiskdiisr,ndiskdiisc,bl(ilist),
     *                       bl(ipairimages),inewsingles,bl(iorevtable),
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
c
         mp4_if: if (do_mp4) then
            call matconn('fockMO','q',ncf,ncf,ifockMO)
            total_singl_en=0.0d0
            call matmmult('spwork1','fockMO','spwork2')
            call matprodtr('spwork2','spwork1',part_singl)
            total_singl_en=total_singl_en+part_singl
            call matmmult('fockMO','spwork1','spwork2')
            call matprodtr('spwork2','spwork1',part_singl)
            total_singl_en=total_singl_en-part_singl
            total_singl_en=-2.0d0*total_singl_en
            write(6,*) 'Total singles MP4 energy: ', total_singl_en
            call initsingles(nval,ncf,singles.or.mp4,nmo,vorb,nbf)
            call matdisc('fockMO')
            if (iprnt.ge.3) write(6,*) 'Sum of singles squares:   ', ssq
            if (sing_max.gt.Wavethresh) goto 665
            xmp4_s=total_singl_en
         else
            if (iprnt.ge.3) write(6,*) 'Sum of singles squares:   ', ssq
         endif mp4_if
      endif sing
      if (do_mp4) then
         call dynamic_matrem('spwork3')
         call dynamic_matrem('spwork2')
         call dynamic_matrem('spwork1')
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(13)=t1-t0        ! Singles Residuals
      ccsdela(13)=et1-et0      ! Singles Residuals
c This is the singles residuum loop, STOP
c
      if (iprnt.ge.2) then
         write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * *'
         write(6,*) '     The results for iteration: ',iteration
         write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * *'
         write(6,*)  'Sum of the squares of the residuums:',square
         write(6,*)  'Maximum residuum element (abs value)',xmaxele
      endif
      if (mp3.or.do_mp4) then 
         if (mp3.and.iprnt.ge.2) write (6,'(A,F25.16)') 
     &                           'MP3 Energy is: ', energy
         energy_mp3=energy
         if (do_mp4) then 
            xmp4_q=energy
         else
            xmp3=energy
         endif
         goto 70
      endif
      if (iprnt.ge.2) then
         if (cep0) write(6,'(A)') '* * * * * * * * * * * * * * * *'// 
     &             'CEPA-0 Energies: * * * * * * * * * * * * * * *'
         write(6,35) energy
         write(6,36) energyr
         write(6,37) xlinear_energy
         write(6,39) varenergy
         write(6,41) varenergy-varenergyold
         write(6,40) totnorm
         write(6,38) energyr-energyrold
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
 42   FORMAT('SCS energy:  ',F21.16)
      DeltaE=dabs(energyr-energyrold)
c     DeltaE=dabs(xlinear_energy-energy)
      DeltaE1=    energyr-energyrold
      energyrold=energyr
      varenergyold=varenergy
      if (variational) then
         xiterenergy=varenergy ! for CID quadratic corrected
         if (iprnt.ge.3) then
            write(6,'(A)') 'Variational energy will be used for the'
     *          //' next iteration.'
            write(6,'(A)')'In fact this energy is not variational'//
     *                    ' when Gauss-Seidel algorithm was used!'
         endif
      else if (linear) then
         xiterenergy=xlinear_energy   ! for CID linear
         if (iprnt.ge.3) then
            write(6,'(A)')'Linear energy will be used for the'//
     &                    ' next iteration.'
         endif
      else
         xiterenergy=energyr ! for CID quadratic (default)
         if (iprnt.ge.3) then
            write(6,'(A)') 'Quadratic (non-variational!) energy will'//
     &                     ' be used for the next iteration.'
         endif
      endif
c DIIS!
      call secund(t0)
      call elapsec(et0)
      if (diis.eq.1.or.(dabs(xmaxele).lt.5.0d-2.and.
     *                   ifirst_diis_iter.le.iteration)) then
      if (diis.ne.1.and.iprnt.ge.3) 
     *    write(6,*) 'Calculation of diis matrix started.'
      diis=1
      call CCDiis(i,       j,          'calculate',  ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'ccc',bl(ilist),
     *            bl(ipairimages),nbf)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(14)=t1-t0        ! DIIS
      ccsdela(14)=et1-et0      ! DIIS
c DIIS!
c---------------------------------------------------------------------
      call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
     *              vorb,  nmo, ivalpair,     ivirpair,   ifp, 
     *              nsym,  nval,nvirt,        ij_unique,  isympairs,
     *              nbf)
      call initsingles(nval,ncf,singles,nmo,vorb,nbf)
c---------------------------------------------------------------------
      call elapsec(eistop)
      call secund(sistop)
      ccsdcpu(12)=sistop-sistart      ! Tot. time without singles res and DIIS
      ccsdela(12)=eistop-eistart      ! Tot. time without singles res and DIIS
      if (iprnt.ge.3) call print_CCSD_results(ccsdela,ccsdcpu,itim_no)
      write(6,3232) iteration,energyr,DeltaE1,xmaxele,square,
     *              ccsdela(12)/6.0d1

      if (xmaxele.gt.Wavethresh .or. DeltaE.gt.Ethresh) then
         if (iteration.ge.maxiter) then
            write(6,*) '            * * * * Maximum number of'//
     *             ' iterations reached * * * *'
            write(6,*) '            * * * *             NO '//
     *             'CONVERGENCE           * * * *'
         endif
      else
         converged=.true.
      endif
      goto 666

888   allocate(exc_LMO(nval,nval,idimen,idimen))
      allocate(coef_LMO(nval,nval,idimen,idimen))
      allocate(resi_LMO(nval,nval,idimen,idimen))
      exc_LMO=0.0D0; coef_LMO=0.0D0; resi_LMO=0.0D0
      do i=1,nval
         do j=1,nval
            do k=1,idimen
               do l=1,idimen
                  do ii=1,nval
                     exc_LMO(i,j,k,l)=exc_LMO(i,j,k,l)
     &                                +exc_QCMO(ii,j,k,l)*TX(ii,i)
                     coef_LMO(i,j,k,l)=coef_LMO(i,j,k,l)
     &                                +coef_QCMO(ii,j,k,l)*TX(ii,i)
                     resi_LMO(i,j,k,l)=resi_LMO(i,j,k,l)
     &                                +resi_QCMO(ii,j,k,l)*TX(ii,i)
                  enddo
               enddo
            enddo
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
      write(6,'(A,F20.9)') 'E(CORR-SUBSYS) =',Ecluster
      write(6,'(A,F20.9)') 'E(CORR-CENTRAL)=',Ecen

 70   continue
      call flush(6)
      if (do_mp4) goto 222
      if (mp4 .and. .not.do_mp4) then 
         mp3=.false.
         do_mp4=.true.
         goto 666
      endif
      call flush(6)
      if ((cc.and.singles).and.l_triples) then
         call elapsec(xt0)
         call sort_amplit(nval,idimen,ndisk_a,af,0,
     *                 ichar_count,bl(ivtable),bl(im_table),bl(ijtable),
     *                 bl(ivrevtable),c1_ratio)
         call simple_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                 ndisktt,ncf,lbintt,thresh,byt8,
     *                 nmo,vorb,af,0,ichar_count,bl(ivtable),
     *                 bl(im_table),bl(ivrevtable),bl(iotable))
         call make_3ext_pairs(idimen,nval,af,ndisk_ie,ndisk_ie1)
         call sort_3int(nval,idimen,ndisk_ii,af,0)
         amp_ratio=1.0d0
         ext3_ratio=1.0d0
         nslv=1
         call dynamic_show_free(mem)
         mem=mem-nval*nval*nval-3*nval*nval
         call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                        nslv,npass,isize)
         call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *                  ncf,lbinx,thresh,byt8,nmo,vorb,af,0,isize,npass)
         call elapsec(xt1)
         if (iprnt.ge.3) write(6,*) 'Sort time: ',xt1-xt0
         call flush(6)
         call elapsec(xt0)
         call dynamic_getmem(20,itimes)
         call new_triples(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,
     &                    ndisk_ix,res,bl(iepsi+ncore),bl(itimes),20,
     &                    ccsd,qcisd,af,bl(iotable),bl(ivtable),
     &                    bl(ivrevtable),bl(im_table),bl(ijtable),
     &                    ichar_count,npass,isize,energys,energyd,iprnt)
         call elapsec(xt1)
         if (iprnt.ge.3) then
            write(6,*) res,'   Time: ',xt1-xt0
            write(6,'(A,2F25.15)') "Energies, s & d: ", energys, energyd
            write(6,*) 'W zero & sort:           ',bl(itimes+0)
            write(6,*) 'W build (total):         ',bl(itimes+1)
            write(6,*) 'Mult. over virt space:   ',bl(itimes+2)
            write(6,*) 'Mult. over occ. space:   ',bl(itimes+3)
            write(6,*) 'Reading of ampl. & int:  ',bl(itimes+4)
            write(6,*) 'Amplitudes:              ',bl(itimes+5)
            write(6,*) '3int integrals:          ',bl(itimes+6)
            write(6,*) '3ext integrals:          ',bl(itimes+7)
            write(6,*) 'Kext integrals:          ',bl(itimes+9)
            write(6,*) 'Energy:                  ',bl(itimes+8)
            call flush(6)
         endif
      endif
c
 999  continue
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
     *      //   ' * * * * * * * * * * * *'
      write(6,*) '                            Correlation energies:    '
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
            res=0.0d0
         endif
         write(6,'(A34,F22.9)') 'MP4 quadruples energy:       ',xmp4_q
         write(6,'(A34,F22.9)') 'Total MP4 correlation energy:',
     *      xmp2+xmp3+xmp4_s+xmp4_d+xmp4_q+res
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
         write(6,'(A34,F22.9)')'SCS-MP2 total energy:        ',
     &                                                       scsmp2+erhf
         write(6,'(A34,F22.9)')'E(CORR-SUBSYS) =             ',xmp2
         write(6,'(A34,F22.9)')'E(CORR-CENTRAL)=             ',ecenmp2
         if(method(1:3).EQ.'MP2') GO TO 95     ! JB July 2010
         write(6,'(5X,A29,F22.9)')
     *      method(1:lmet)//' correlation energy:'//spc,energyr
         if (l_triples) then
            write(6,'(5X,A29,F22.9)')
     *      method(1:lmet)//' triples correction:'//spc,res
         else
            res=0.0d0
         endif
         etot = energyr+erhf+res
      endif
 95   CONTINUE
c 
      call memory_status('end of MP2 or CCSD')

      if (singles) then
         call pointersingles(ipointer)
         call orb_pairs_check(idimen,bl(ipointer),1d-7)
      endif
c
c -- write energy to control file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,idum,etot,chopv)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      call allclose_and_delete(.false.)
      call dynamic_retmark
      call matremark
      call retmark
      end subroutine CIM_CC


C ****************************************************************
C * Express the integrals by the index of occupied MOs           *
C * It is convenient for transforming integrals from QCMO to LMO *
C * NZG_6/9/2016 @UARK                                           *
C ****************************************************************
      subroutine IntExpress(i,j,exc,coef,resi,exc_QCMO,coef_QCMO,
     &                      resi_QCMO,idimen,nval)
      implicit none
      
      integer i,j,idimen,nval,a,b
      real*8 exc(idimen,idimen),coef(idimen,idimen),resi(idimen,idimen)
      real*8 exc_QCMO(nval,nval,idimen,idimen)
      real*8 coef_QCMO(nval,nval,idimen,idimen)
      real*8 resi_QCMO(nval,nval,idimen,idimen)

      do a=1,idimen
         do b=1,idimen
            exc_QCMO(i,j,a,b)=exc(a,b)
            coef_QCMO(i,j,a,b)=coef(a,b)
            resi_QCMO(i,j,a,b)=resi(a,b)
         enddo
      enddo

      if (i/=j) then
         do a=1,idimen
            do b=1,idimen
               exc_QCMO(j,i,a,b)=exc_QCMO(i,j,b,a)
               coef_QCMO(j,i,a,b)=coef_QCMO(i,j,b,a)
               resi_QCMO(j,i,a,b)=resi_QCMO(i,j,b,a)
            enddo
         enddo
      endif 
 
      end subroutine IntExpress


C ****************************************************************
C * Express Fock matrix from triangle form to symmetric matrix   *
C * This Fock matrix has the atom order, it need to be reordered *
C * NZG_6/13/2016 @UARK                                          *
C ****************************************************************
      subroutine FockFull(ncf,LF2,dtmp,Fock1)
      implicit none
      
      integer ncf,LF2,i,j,k
      real*8 dtmp(LF2),Fock1(ncf,ncf)

      k=0
      do i=1,ncf
         do j=1,i
            k=k+1
            Fock1(j,i)=dtmp(k)
            Fock1(i,j)=dtmp(k)
         enddo
      enddo

      end subroutine FockFull


      SUBROUTINE EECIM(NO,TX,EK,ECIM)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      DIMENSION TX(NO,NO),EK(NO,NO),ECIM(NO)
C
      ECIM=0.0D0  
      DO II=1,NO
         DO I=1,NO
            DO J=1,NO
               ECIM(II)=ECIM(II)+TX(I,II)*TX(J,II)*EK(J,I)
            END DO
         END DO
      END DO
C
      RETURN
      END


      subroutine getepsi(Fock,epsi,ncf)
      implicit none

      integer ncf,i
      real*8 Fock(ncf,ncf),epsi(ncf)
      
      do i=1,ncf
         epsi(i)=Fock(i,i)
      enddo
      
      end subroutine getepsi

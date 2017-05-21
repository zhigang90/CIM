#include "maxi.h"
#define DUMP_AMPLITxUDES
      Subroutine coupled_cluster
c Main parallel CCSD soubroutine. It uses its own memory management.
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
c-----------------------------------------------------------------------
c
      data opnames  /'corr','trip','nofr','thre','core','orbs',
     2               'lvsh','diis','conv','memo','prin','iter',
     3               'chkp','dump','rest','chun','loca','byt5',
     4               'aorb','cach'/
c
      data ioptyp /21, 0, 0,11,11, 2,11, 2,12,21, 1, 1,
     $              0, 0, 0, 2, 0, 0, 0, 1/
c
      data methods/'mp2   ','mp3   ','mp4   ','cid   ','cisd  ',
     $             'cepa-0','cepa-2','qcisd ','ccd   ','ccsd  '/
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
      write(6,'(A)') "                    Closed-Shell Coupled"//
     *               " Cluster module"
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
      vorb=.true.
      cache=.false.
      variational=.true.
      integerT=.false.
      xmaxdisk=2 500 000 000.0d0      ! currently not active
      diis_thres=0.01d0
      diis_less_thr=.false.
c
c -- set up name of dump file for restart
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      dumpfile=scrfile(1:len)//'.cc_dump'
      lenD = len+8
c
c  basis function symmetry pairs are in bl(ifp)
      iprnt=1
      natoms=igetival('na')
c
      call matreset
c
c .............................................................
c -- CHECK! Only do CORR if current wavefunction is RHF or RMP2
c --   (at the same time get the lowest eigenvalue of the overlap)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,wvfnc)
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(wvfnc(1:3).EQ.'RHF' .or. wvfnc(1:3).EQ.'RMP') then
      else
        call nerror(1,'Correlated Energy Module',
     $  'Correlated Energy Requested but System is NOT RHF or RMP2 !',
     $              0,0)
      endif
c .............................................................
c
      call readop1(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c
c -- determine correlation method
      imethod=0
      If(ifound(1).GT.0) Then      ! must be true
        method = chopv(1)(1:6)
        call lowerca2(method,6)
        Do i=1,nmethod
        if(method.eq.methods(i)) then
          imethod = i
          exit
        endif
        EndDo
      EndIf
c
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
      else
        write(6,'(A)')  'No valid correlation method declared'
        STOP 'ERROR! No valid correlation method declared'
      EndIf
c
c -- process other options
c
c -- TRIP
      l_triples=.false.
      itrstart=0
      if(ifound(2).gt.0) l_triples=.true.
cc  NEED TO SORT OUT BELOW
cc      if(ifound(32).gt.0) then
cc        l_triples=.true.
cc        itrstart=iopv(1,32)
cc        write(6,*) "Triples correction, start from abc=",itrstart
cc      else
cc        l_triples=.false.
cc      endif
c
c -- NOFR
      nofr=.false.
      if(ifound(3).gt.0) nofr=.true.
c
c------------------------------------------------------------------
c -- THRE (integral threshold)
      if(ifound(4).gt.0) then
        thresh=10.0d0**(-ropv(1,4))
      else
c -- check lowest eigenvalue of overlap (from SCF)
        thresh = MIN(1.0d-12,xlow**2)
        If(thresh.LT.1.0d-14) then
          write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final energy may have greater than mhartree error')
          thresh = 1.0d-14
        EndIf
      end if
      call setrval('thresh',thresh)
c---------------------------------------------------------------------
c
c -- CORE is the limit separating core from valence
      core=-3.0d0
      if(ifound(5).gt.0)  core=ropv(1,5)
c
c -- ORBS
c     see below
c
c -- LVSH
      if(ifound(7).ne.0) then
        shift=ropv(1,7)
        write(6,'(A,F6.1)')  'Level shift set to ',shift
      else
        shift=0.0d0
        write(6,'(A,F6.1)')  'Default level shift ',shift
      endif
      shiftini=0.0d0      ! initial level shift for mp2 part
c
c -- DIIS
      if(ifound(8).ne.0)  then
        length_diis=iopv(1,8)
        ifirst_diis_iter=iopv(2,8)
        if (ifirst_diis_iter.lt.2) ifirst_diis_iter=1
        if (length_diis.lt.2) then
          log_diis=.false.
        else
          log_diis=.true.
        endif
      else
        log_diis=.true.
        length_diis=6
        ifirst_diis_iter=1
      endif
c
c -- CONV
      if(ifound(9).gt.0) then
        Ethresh=   10.d0**(-ropv(1,9))
        Wavethresh=10.d0**(-ropv(2,9))
      else
        Ethresh=1.0d-6
        Wavethresh=1.0d-6
      endif
c
c -- MEMO
      if(ifound(10).gt.0) then
        max_dyn_mem = memfunc(chopv(10))
        if (max_dyn_mem.lt.2000) max_dyn_mem=max_dyn_mem*1000000
      else
        max_dyn_mem=10000000
      endif
      call dynamic_init(max_dyn_mem,bl(1))
c
c -- PRINt
      if(ifound(11).gt.0) iprnt=iopv(1,11)
      call setival('ccprint',iprnt)
c
c -- ITER
      maxiter=20
      if(ifound(12).gt.0) maxiter=iopv(1,12)
c
c -- CHKP/DUMP
      if(ifound(13).gt.0.or.ifound(14).gt.0) then
        dump=.true.
        if (ifound(14).gt.0) then
          last_dump_only=.true.
        else
          last_dump_only=.false.
        endif
        if(ifound(15).eq.0) then
          OPEN (UNIT=ndisk_dump,FILE=dumpfile(1:lenD),
     1          FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE')
        endif
      else
        dump=.false.
      endif
c
c -- RESTart
      if(ifound(15).gt.0) then
        l_restart=.true.
        OPEN (UNIT=ndisk_dump,FILE=dumpfile(1:lenD),
     #        FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='OLD')
      else
        l_restart=.false.
      endif
c
c -- CHUNk
      if(ifound(16).gt.0) then
        ichunk=iopv(1,16)
        jchunk=iopv(2,16)
        calchunks=.false.
      else
        calchunks=.true.
      endif
c
c -- BYT5
      if (ifound(18).ne.0)  then
        byt8=.false.
        write(6,'(A)') 'Integrals will be stored in 5 byte format'
      else
        write(6,'(A)')
     $          'Integrals will be stored in full double precision'
      endif
c
c -- AORB
      if(ifound(19).gt.0) then
        vorb=.false.
        write(6,*) 'Calculation using AO basis'
      else
        vorb=.true.
        write(6,*) 'Calculation using MO basis'
      endif
c
c -- CACHe
      if(ifound(20).gt.0) then
        cache=.true.
        write(6,*) 'Calculation using local data caching'
      else
        cache=.false.
        write(6,*) 'Calculation using af only'
      endif
c
c-----------------------------------------------------------
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call para_bcast(TxDoCCSD,TxJobType)
      call fafinit(iresult)
c
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
c-----------------------------------------------------------
c -- check basis set (G is highest angular momentum basis function coded)
      call ChkAngMom(ncs,bl(ictr),2)
c-----------------------------------------------------------
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
          call nerror(2,'newtrans',
     1   'Program cannot find the number of MOs on <jobname>.14',
     2    ni,ni)
        end if
      end if
      nmo=xnmo(1)  ! no of occupied orbitals
      erhf=xnmo(2) ! HF energy
c----------------------------------------------
c define some matrices,
      nbf=ncf
      call dynamic_matdef('cano','q',ncf,nbf)
      call dynamic_matdef('epsi','d',nbf,nbf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      if(ifound(17).eq.0) then
        write(6,*)'Canonical orbitals are used'
        if (l_restart) then
          call general_read1(ndisk_dump,bl(icano),ncf*nbf*8)
          call general_read1(ndisk_dump,bl(iepsi),nbf*8)
        else
          call matread('cano',np4,'evec_rhf') ! MOs
          call matread('epsi',np4,'eval_rhf') ! HF eigenvalues
        endif
      else
        write(6,*)'Localized orbitals are used'
        if (l_restart) then
          call general_read1(ndisk_dump,bl(icano),ncf*nbf*8)
          call general_read1(ndisk_dump,bl(iepsi),nbf*8)
        else
          call matread('cano',np4,'loca_rhf') ! MOs
          call matread('epsi',np4,'eval_rhf') ! HF eigenvalues
          call dynamic_matdef('loc-epsi','d',nmo,nmo)
          call matread('loc-epsi',np4,'eloc_rhf') ! HF eigenvalues
          ilocepsi=mataddr('loc-epsi')
          call tfer(bl(ilocepsi),bl(iepsi),nmo)
          call dynamic_matrem('loc-epsi')
          loca=.true.
        endif
      endif
c
c -- test for redundant basis functions
      call tstival('nonredun',iyes)
      if (iyes.ne.0) then 
        nbf=igetival('nonredun')
      else
      do i=1,ncf
        xsum=0d0
        do k=1,ncf
          xx=bl(icano+(k-1)+(i-1)*ncf)
          xsum=xsum+xx*xx
        enddo
        if (xsum.lt.1d-5) then
          nbf=i-1
          exit
        endif
      enddo
      endif
c
      call dynamic_lock(bl(icano),i)
      call dynamic_lock(bl(iepsi),i)
c
      if (log_diis) then 
        write(6,'(A,I5)')'DIIS expansion length: ',length_diis
        write(6,'(A,I5)')'DIIS evaluation start: ',ifirst_diis_iter
      else
        write(6,'(A)') 'Calculation without DIIS'
      endif
      write(6,'(A,F6.1)') 'Level shift set to ',shift
      write(6,*) 'Threshold for energy  convergence:', Ethresh
      write(6,*) 'Threshold for wvfn.   convergence:', Wavethresh
C
C  Determine the orbitals to be correlated
C
      If(ifound(6).gt.0) Then
        nfirst=iopv(1,6)
        if(nfirst.LT.1.or.nfirst.GT.nmo)
     $     call nerror(5,'Correlated Energy Module',
     $      'Lowest Orbital in Correlation Range Incorrect',0,0)
        nlast=nmo
        nval=nlast-nfirst+1
        ncore=nlast-nval
        nend=iopv(2,6)
        if(nend.LE.nmo.or.nend.GT.nbf)
     $     call nerror(6,'Correlated Energy Module',
     $      'Highest Orbital in Correlation Range Incorrect',0,0)
        nbf=nend
      Else
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
      EndIf
c
      call matsub('genvirt','cano',nmo+1,nbf)
      nvirt=nbf-nmo      ! number of virtual orbitals
c
      if (vorb) then
        idimen=nvirt
      else
        idimen=ncf
      endif
      npairs=(nval+1)*nval/2
c
c -- check on the sizes of any requested chunks
      if(ifound(16).gt.0) then
        call dynamic_show_free(memavail)
        if((ichunk*jchunk+ichunk+jchunk)*idimen**2.gt.memavail) then
          call message('**WARNING** from CORR module',
     $                 '  Size of Chunks may cause memory overflow',0,0)
        endif
      endif
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
      if(ncore.gt.0) call dynamic_getmem(7*ncore, icorpair)
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      icorr=icano+ncore*ncf
      call symm_memory_allocator(ncf)
      if (.not.loca.and..not.l_restart) then
      call dynamic_getmem(ncf*ncf,iscr1)
      call dynamic_getmem(ncf*ncf,iscr2)
      call dynamic_getmem(ncf*ncf,iscr3)
      call adv_symmetrize(nsym,   ncf,  ncore,   bl(ifp),
     *                    bl(icano),bl(ioverlap),
     *                    bl(iscr1),bl(iscr2),bl(iscr3),bl(iepsi))
      call adv_symmetrize(nsym,   ncf,  nval,   bl(ifp),
     *                    bl(icorr),bl(ioverlap),
     *                    bl(iscr1),bl(iscr2),bl(iscr3),bl(iepsi+ncore))
      call adv_symmetrize(nsym,   ncf,  nvirt,   bl(ifp),
     *                    bl(icorr+nval*ncf),bl(ioverlap),
     *                    bl(iscr1),bl(iscr2),bl(iscr3),bl(iepsi+nmo))
      call dynamic_retmem(3)
      endif
      if(ncore.gt.0)
     *  call core_pairs(nsym,   ncf,  ncore,   bl(ifp),  bl(icorpair),
     *                  bl(icano),bl(ioverlap))
      call orb_pairs(nsym,ncf,nval,nvirt,bl(ifp),bl(ivalpair),
     *               bl(ivirpair),bl(icorr),bl(ioverlap))
cc
      docansym=.false.
      if (nsym.gt.0.and..not.loca) docansym=.true.
c
      call dynamic_getmem(ncf*ncf,ifockAO)
      if (.not.l_restart) call FockRead(ncf,ifockAO)
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
     *                   bl(ichatacters),ichar_count,bl(im_table))
      call elapsec(et1)
      write(6,'(A,I5)') 'Number of irreps: ',ichar_count
  9   FORMAT (A,F10.2)
      if (iprnt.gt.3)
     *write(6,9)'Time spent in irr_rep_table creation: ',et1-et0
c Fock rediagonalization with symmetry:
c
      if (.not.l_restart) then
      call dynamic_mmark()
      call dynamic_getmem(2*8,      ibassym)
      call dynamic_getmem(ncf*ncf,  imat1)
      call dynamic_getmem(ncf*ncf*8,imat2)
      call elapsec(et0)
      call symm_adapted(ncf,nsym,ichar_count,bl(ichatacters),bl(ifp),
     *                  bl(imat1),bl(imat2),bl(ibassym))
      call elapsec(et1)
      if (iprnt.gt.2)
     *write(6,9)'Time in symm_adapted tranformation finder: ',et1-et0
      call dynamic_retmem(1)
      call elapsec(et0)
      if(ncf.eq.nbf) then
        call Fock_Diag(ncf,ifockAO,ioverlap,icano,iepsi,bl(imat1),
     *                 bl(ibassym),nmo)
      endif
      call elapsec(et1)
      if (iprnt.gt.2)
     *write(6,9)'Time in FockRead+Fock_Diag subroutines:    ',et1-et0
      call dynamic_retmark()
      if(ncore.gt.0)
     *  call core_pairs(nsym,   ncf,  ncore,   bl(ifp),  bl(icorpair),
     *                  bl(icano),bl(ioverlap))
      call orb_pairs(nsym,ncf,nval,nvirt,bl(ifp),bl(ivalpair),
     *               bl(ivirpair),bl(icorr),bl(ioverlap))
      call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                   bl(ichatacters),ichar_count,bl(im_table))
      endif
c
      call tables_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                      bl(ivalpair),bl(iotable),bl(iorevtable))
      call tables_generator(nsym,nvirt,ichar_count,bl(ichatacters),
     *                      bl(ivirpair),bl(ivtable),bl(ivrevtable))
      call ijtab_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                     bl(ivalpair),
     *                     bl(im_table),bl(iorevtable),bl(ijtable))
      endif
cc
      if (loca) then
       if(ncore.gt.0)
     *  call symmetrize_energies(ncf,nsym,ncore,bl(icorpair),bl(iepsi))
        call symmetrize_energies(ncf,nsym,nval, bl(ivalpair),
     *                           bl(iepsi+ncore))
       if(ncore.gt.0)
     *  call symmetrizer2(nsym,  ncf,  ncore, bl(ifp),  bl(icorpair),
     *                    bl(icano))
        call symmetrizer2(nsym,  ncf,  nval,  bl(ifp),  bl(ivalpair),
     *                    bl(icorr))
        call symmetrizer2(nsym,  ncf,  nvirt,  bl(ifp),  bl(ivirpair),
     *                    bl(icorr+nval*ncf))
      endif
      call dynamic_lock(bl(iepsi),i)
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
c                                               pair has, if given pair
c                                               is already an image, = 0
      call izeroit(bl(ipairimages),npairs+1)
      call int_array_write(bl(ipairimages),npairs+1,1)
      do i=1,nval
        do j=1,i
          ij=ij+1
          call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                       ns, trans, signum)
        if (ns.eq.0) then
        if (iprnt.gt.3) write(6,'(A,I2,A,I2,I5)')"Pair: ",i,",",j,ij
        ij_unique=ij_unique+1
        call int_array_write(bl(isympairs),ij,ij_unique)
        else
        call int_array_write(bl(isympairs),ij,0)
        if (iprnt.gt.3)
     *    write(6,'(A,I2,A,I2,I5,A,I2,A,I2,I5,A,I1,L2,L2)')
     *           "Pair: ",i,",",j,ij," image: ",
     *           iprim,",",jprim,ijprim," ns: ",ns,trans,signum
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
      if (iprnt.gt.2)
     *write(6,9) 'Elapsec time for coulomb integral calc. (sec): ',
     * ec_stop-ec_start
      if (iprnt.gt.2)
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
     *                bl(ilist),bl(isympairs),ij_unique,iprnt,nbf)
      call secund(sx_stop)
      call elapsec(ex_stop)
      if (iprnt.gt.2)
     *write(6,9) 'Elapsec time for exchange integral calc. (sec): ',
     * ex_stop-ex_start
      call flush(6)
      if (iprnt.gt.2)
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
      call master_Gen(ncs,     ncf,           ictr,   nval,  nmo,
     *                nfirst,  nlast,         thresh, core,  xmaxdisk,
     *                ndisktt, bl(irecadrtt), lrectt,  lbintt, int_kind,
     *                nodisk,  byt8,          small,  vorb,  ndisktrtt,
     *                bl(ilist),bl(isympairs),ij_unique,iprnt,nbf)
      call secund(sx_stop)
      call elapsec(ex_stop)
      if (iprnt.gt.2)
     *write(6,9) 'Elapsec time for 3ext coulomb integral calc. (sec): ',
     * ex_stop-ex_start
      call flush(6)
      if (iprnt.gt.2)
     *write(6,9) 'CPU time for 3ext coulomb integral calc. (sec):     ',
     * sx_stop-sx_start
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
     *                   nfirst,nlast,byt8,nmo,vorb,.true.,kijkndisk,
     *                   nbf)
      call elapsec(et1)
      call secund(t1)
      if (iprnt.gt.2) then
        write(6,9)'Kijk_Vec_Init ELAPS: ',et1-et0
        write(6,9)'Kijk_Vec_Init   CPU: ',t1-t0
        call flush(6)
      endif
      endif
c
      if (mp2) lbinc=0 ! FRAGILE!
      call elapsec(et0)
      call secund(t0)
      if (l_restart) then
      call general_read1(ndisk_dump,bl(ifockAO),idimen*idimen*8)
      call general_read1(ndisk_dump,bl(ifockMO),ncf*ncf*8)
      call general_read1(ndisk_dump,bl(icorefock),idimen*idimen*8)
      else
      call FockBuilder(bl(irecadrc),bl(irecadrx),nval, npairs, thresh,
     *                 ndiskc,      ndiskx,      lbinc,lbinx,  byt8, 
     *                 nmo,        vorb,     ifockAO, ifockMO,icorefock,
     *                 nbf)
      endif
      call sync_barrier
      call elapsec(et1)
      call secund(t1)
      if (iprnt.gt.2) then
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
     *                      thresh,nfirst,nlast,bl(ifockMO),byt8,
     *                      nmo,vorb,kijklndisk,af,nslv)
      call sync_barrier
      call elapsec(et1)
      call secund(t1)
      if (iprnt.gt.2) then
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
ctmp:
#if 0
      call dynamic_getmem(ncf,iad2)
      call dynamic_getmem(ncf,iad3)
      call osinv(bl(iinvcano),ncf,det,1.0d-18,bl(iad2),bl(iad3))
      call dynamic_retmem(2)
      if (dabs(det).lt.1d-200) then
        write(6,*)'DET: ',det
        STOP 'Determinant too small'
      endif
c
c     inverse refinement:
      call dynamic_matdef('unit','q',ncf,ncf)
      junit=mataddr('unit')
      call dynamic_matdef('delta','q',ncf,ncf)
      do k=1,1
      call matmmult('invcano','cano','unit')
      do i=0,ncf-1
        bl(junit+i+i*ncf)=bl(junit+i+i*ncf)-1.0d0
      enddo
      summ=0.0d0
      xmax=0.0d0
      do i=0,ncf*ncf-1
        ele=dabs(bl(junit+i))
        summ=summ+ele
        if (ele.gt.xmax) xmax=ele
      enddo
      if (iprnt.gt.3)
     *  write(6,272)'Sum_ of abs of epsilon before improvement: '
     *              ,summ,' max: ',xmax
      call matmmult('unit','invcano','delta')
      call matadd1('delta',-1.0d0,'invcano')
      enddo
      call dynamic_matrem('delta')
      call dynamic_matrem('unit')
#endif
 272  FORMAT (2(A,E15.5E2))
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
      if (iprnt.gt.3)
     * write(6,272) 'Final sum: ',summ,' max. element: ',xmax
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
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,   
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      else
      call prepare_diis_files(ndiskdiisr,ndiskdiisc,af,ncf,
     *                        nmo,vorb,nbf)
       i=length_diis
      call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
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
      if (iprnt.gt.2) then
      write(6,9) 'Elapsec time for first iteration:  ',
     * e1_stop-e1_start
      write(6,9) 'CPU time for first iteration:      ',
     * s1_stop-s1_start
      endif
      if (iprnt.gt.3) then
      write(6,*) 'Sum of the squares of the residuums: ',square
      write(6,*) 'Maximum residuum element (abs value) ',xmaxele
      write(6,'(A,F21.17)') 'Nonlocal! - "MP2" energy: ',energy
      endif
      if (.not.loca) then 
         xmp2=energy
         scsmp2=senergy
         if (iprnt.gt.3) write(6,*)scsmp2
         if (mp2) goto 999
      else
      if (iprnt.gt.2)
     *  write(6,'(A,F21.17)') 'Nonlocal! - "MP2" energy: ',energy
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
     *                       ns, trans, signum)
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
            if ((mp2_iter.eq.1.or.mp2last).and.iprnt.gt.3)
     *            write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',k,l,
     *                                       bl(iepairr+kl-1)
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
          if ((mp2_iter.eq.1.or.mp2last).and.kl.ne.0.and.iprnt.gt.3)
     *          write(6,'(A,2I5,F20.12)') 'Pair en: i,j,energy:',k,l,
     *                                     bl(iepairr+kl-1)
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
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'mp2',bl(ilist),
     *            bl(ipairimages),nbf)
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
      if (iprnt.gt.1) then
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
      if (iprnt.gt.1)
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
      if (iprnt.gt.2)
     *write(6,*) "CCDIIS done"
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
      if (iprnt.gt.2)
     *write(6,*) "sort_amplit done"
      call flush(6)
      call elapsec(xt1)
      if (iprnt.gt.2) then
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
      if (iprnt.gt.2)
     *write(6,*)'RATIO: ',ext3_ratio
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
ccs      call sleep(5)
      call elapsec(xt0)
      call dynamic_getmem(20,itimes)
      call new_tri_mas(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,ndisk_ix,
     *                       res,bl(iepsi+ncore),bl(itimes),20,ccsd,
     *                       qcisd,af,nslv,ext3_ratio,amp_ratio,npass,
     *                       isize,iprnt,energys,energyd,itrstart)
      call elapsec(xt1)
      write(6,'(A,F10.1)') 'Triples elapsed time: ',(xt1-xt0)/6d1
      if (iprnt.gt.2) then
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
      iterate=1
      call para_barrier
      call para_pack_real(xmaxele,1)
      call para_pack_int(iterate,1)
      call para_bcast_pack(TxCCSDIter1)
      call para_barrier
      iteration=iteration+1
      if(iprnt.ge.3) then
        write(icond,*) '                 The CCSD iteration: ',iteration
        call f_lush(icond)
      endif
      tot_xmax=xmaxele
      xmaxele=0.0d0
      if (iteration.eq.2) tot_xmax=0.1d0
      call set_thresh(tot_xmax,xthresh_EEO)
      if (cep0) xiterenergy=0.0d0
c     if (iteration.eq.2) then
c       xiterenergy=0.0d0
c       write(6,*)'Iterations starts with energy=0 ("mp3" approximation)'
c     endif
c
      call elapsec(eistart)
      call secund(sistart)
      int_kind='e' ! EE operator
c     call izeroit(bl(irecadre),npairs*nslv)    !zero out the pair records
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
     *                   nfirst, nlast,   thresh,
     *                   vorb,   ndisktre, npairs, nslv,iprnt,nbf)
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
c     call izeroit(bl(irecadrtx),nval*nval)     !zero out the pair records
      call master_Gen(ncs,     ncf,           ictr,   nval,   nmo,
     *                nfirst,  nlast,         thresh, core,   xmaxdisk,
     *                ndisktx, bl(irecadrtx), lrectx, lbintx, int_kind,
     *                nodisk,  byt8,          small,  vorb,   ndisktrtx,
     *                bl(ilist),bl(isympairs),ij_unique,iprnt,nbf)
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(2)=t1-t0        ! TEIO tx
      ccsdela(2)=et1-et0      ! TEIO tx
c
      int_kind='tc' ! TEIO operator
      call secund(t0)
      call elapsec(et0)
c     call izeroit(bl(irecadrtc),nval*nval)     !zero out the pair records
      call master_Gen(ncs,     ncf,           ictr,   nval,   nmo,
     *                nfirst,  nlast,         thresh, core,   xmaxdisk,
     *                ndisktc, bl(irecadrtc), lrectc, lbintc, int_kind,
     *                nodisk,  byt8,          small,  vorb,   ndisktrtc,
     *                bl(ilist),bl(isympairs),ij_unique,iprnt,nbf)
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
     *                 thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                 ndiskalpha,af,nbf)
      call sync_barrier
      call para_reduce(elaps,5,TxCCSDRedu)
      if (iprnt.gt.2) then
      write(6,'(A)')'CCalph:    Read amp: Read Kij:   Transp: Multiply:'
     *            // ' Putresult:'
      write(6,'(10X,5F10.2)')elaps(1),elaps(2),elaps(3),elaps(4),
     *                       elaps(5)
      endif
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(4)=t1-t0        ! CCalphaonce+supplement
      ccsdela(4)=et1-et0      ! CCalphaonce+supplement
      call secund(t0)
      call elapsec(et0)
ccs      call sleep(5)
      call SupplementAlpha_BetaBuild(nval,ndiskalpha,bl(ibeta),ncf,
     *                              nfirst,nmo,vorb,cc,ccsd,bl(ifockMO),
     *                              kijklndisk,af,do_mp4,nbf)
      call sync_barrier
      call secund(t1)
      call elapsec(et1)
      call para_initsend
      call dynamic_unlock(bl(ibeta),ii)
      call para_bcast_real(bl(ibeta),nval*nval,TxCCSDInit)
      if (iprnt.gt.2) then
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
     *                 lbinx,       thresh,norecalc,byt8,  ioverlap,
     *                 nmo,         vorb,  nslv,    nbf)
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
     *                       lbine,thresh,nfirst,nlast,byt8,
     *                       nmo,vorb,nslv,iiS,nbf)
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
      if (iprnt.gt.2) write(6,*)'ichunk,jchunk: ',ichunk,jchunk
      call flush(6)
 133  FORMAT (4I3,F20.10)
      call prepare_CCYZ_file(ncf,nmo,vorb,ndiskYZ,.true.,nbf)
      call sync_barrier
      if (iprnt.gt.3) write(6,*) 'I am master, ndiskYZ: ',ndiskYZ
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
cc        if (iprnt.ge.2) write(91,*) 'icalc_pairs: ',icalc_pairs
cc      call flush(6)
      if (calchunks)
     *  call check_ijchunk_reduce(icalc_pairs,nval,ichunk,jchunk,
     *                            reduced)
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
cc            if (iprnt.ge.2)
cc     *       write(91,'(A,4I6,I10)') 'Block: i,j: ',kstart,kstop,lstart,
cc     *                                         lstop,islid
cc            call flush(6)
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
cc            if (iprnt.ge.2)
cc     *      write(91,'(A,4I6,I10)') 'Block: i,j: ',kstart,kstop,lstart,
cc     *                                         lstop,islid
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
      if (iprnt.gt.3) then
      write(6,*)'YZ matrices:        read,     calc,    write,read kijk'
     * //', read once'
      write(6,'(A,5F10.2)')'               ',
     *                   elaps(1),elaps(2),elaps(3),elaps(4),elaps(5)
      endif
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
      if (calchunks)
     *  call check_ijchunk_reduce(icalc_pairs,nval,ichunk,jchunk,
     *                            reduced)
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
cc            if (iprnt.ge.2)
cc     *      write(91,'(A,4I6,I10)') 'Block: i,j: ',kstart,kstop,lstart,
cc     *                                         lstop,islid
cc            call flush(6)
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
cc            if (iprnt.ge.2)
cc     *      write(91,'(A,4I6,I10)') 'Block: i,j: ',kstart,kstop,lstart,
cc     *                                         lstop,islid
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
      energy=0.0d0
      energyr=0.0d0
      square=0.0d0
      totnorm=0.0d0
      ij=0
c     if (iteration.eq.3) then
c     OPEN(UNIT=99,FILE='paraT3',FORM='UNFORMATTED',
c    *     ACCESS='SEQUENTIAL')
c       do ii=1,nval
c         do jj=1,ii
c           call CoefRead('tt',ii,jj,icoeffAO)
c           call general_write1(99,bl(icoeffAO),8*idimen*idimen)
c         enddo
c       enddo
c     STOP 'test'
c     endif
      do i=1,nval
        do j=1,i
          ij=ij+1
          call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                       ns, trans, signum)
          if (ns.ne.0) then
            cycle
          endif
cc          if (iprnt.gt.2) write(91,*) 'will send ij: ',ij
cc          call flush(6)
          kl=0
          call para_recv_pack(islgid,TxCCSDReq1)
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
            call para_unpack_real(sqij,1)
            call para_unpack_real(xmaxele_sl,1)
            call para_unpack_real(xnormij,1)
            bl(iepairr+kl-1)=epairr
            bl(iepair+kl-1)=epair
            bl(ixepair+kl-1)=xepair
            xlinear_energy=xlinear_energy+bl(ixepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            energyr=energyr+bl(iepairr+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            square=square+sqij
     *                      *dble(int_array(bl(ipairimages),kl))
            energy=energy+bl(iepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
            totnorm=totnorm+xnormij
     *                      *dble(int_array(bl(ipairimages),kl))
cc            write(icond,'(A,2I5,F20.12)') 'Pair en:  i,j,energy:',k,l,
cc     *                               bl(iepairr+kl-1)
cc            call flush(icond)
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
        enddo
      enddo
      do i=1,nslv
        kl=0
        call para_recv_pack(islgid,TxCCSDReq1)
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
          call para_unpack_real(sqij,1)
          call para_unpack_real(xmaxele_sl,1)
          call para_unpack_real(xnormij,1)
          bl(iepairr+kl-1)=epairr
          bl(iepair+kl-1)=epair
          bl(ixepair+kl-1)=xepair
          xlinear_energy=xlinear_energy+bl(ixepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
          energyr=energyr+bl(iepairr+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
          square=square+sqij
     *                      *dble(int_array(bl(ipairimages),kl))
          energy=energy+bl(iepair+kl-1)
     *                      *dble(int_array(bl(ipairimages),kl))
          totnorm=totnorm+xnormij
     *                      *dble(int_array(bl(ipairimages),kl))
cc          write(icond,'(A,2I5,F20.12)') 'Pair en:  i,j,energy:',k,l,
cc     *                             bl(iepairr+kl-1)
cc          call flush(icond)
          if (xmaxele_sl.gt.xmaxele) xmaxele=xmaxele_sl
        endif
c No work:
        iwork=0
        call para_send(iwork,islgid,TxCCSDJob1)
      enddo
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
     *                    is_resid)
      if ((.not.cc).and.(.not.ccsd).and.(.not.do_mp4)) then
        call cisd_energy_add(ncf,idimen,nval,xiterenergy,vorb,
     *                       ioverlap,is_resid)
      endif
      if (do_mp4) then
        call matmmul2('fockAO','full_sing','s_resid','n','n','a')
      endif
      call matdisc('full_sing')
      call UpdateSingles1(isingles,is_resid,nmo,nfirst,
     *                ncf,npairs,nval,bl(iepsi),shift,ssq,vorb,
     *                sing_max,ndiskdiisr,ndiskdiisc,bl(ilist),
     *                bl(ipairimages),inewsingles,
     *                bl(iorevtable),bl(ivrevtable),docansym,t1diagnost,
     *                nbf)
c
      if (do_mp4) then
        do i=1,nval
          isin=inewsingles+(i-1)*idimen
          call storesingles(i,bl(isin),bl(ispwork1),ncf,idimen,nfirst)
        enddo
      endif
      call dynamic_matrem('s_resid')
      if (do_mp4) then
        call matconn('fockMO','q',ncf,ncf,ifockMO)
        total_singl_en=0.0d0
        call matmmult('spwork1','fockMO','spwork2')
c       call compact_print(bl(ifockMO),ncf)
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
      if (iprnt.gt.2) write(6,*)'Sum of singles squares: ', ssq
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
      if (iprnt.gt.1) then
      write(6,8) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,'(A,I4)') '          The results for iteration:',iteration
      write(6,8) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,11) 'Sum of squares of residuums:         ',square
      write(6,11) 'Maximum residuum element (abs value):',xmaxele
      write(6,12) 'T1 diagnostic: ',t1diagnost
      endif
      if (mp3.or.do_mp4) then 
        if (mp3.and.iprnt.gt.1) 
     *    write (*,'(A,F25.16)') 'MP3 Energy is: ', energy
        if (do_mp4) then 
          xmp4_q=energy
        else
          xmp3=energy
        endif
      else
      if (iprnt.gt.1) then
      if (cep0) write(6,'(A)') '* * * * * * * * * * * * * * * * CEPA-0'
     *//' Energies: * * * * * * * * * * * * * * *'
      if (cep2) write(6,'(A)') '* * * * * * * * * * * * * * * * CEPA-2'
     *//' Energies: * * * * * * * * * * * * * * *'
      write(6,35) energy
      write(6,36) energyr
      write(6,37) xlinear_energy
      if (iprnt.gt.2) then
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
      if (iprnt.gt.1) then
      write(6,'(A)')'Variational energy will be used for the'
     *          //' next iteration.'
      write(6,'(A)') 'Remember! In fact this energy is not variational'
     *           //' when Gauss-Seidel'
      write(6,'(A)')'                             algorithm was used'
      endif
      else if (linear) then
      xiterenergy=xlinear_energy   ! for CID linear
      if (iprnt.gt.1) then
      write(6,'(A)')'Linear energy will be used for the next iteration.'
      endif
      else
      xiterenergy=energyr ! for CID quadratic (default)
      if (iprnt.gt.1) then
      write(6,'(A)') 'Quadratic (non-variational!) energy will be used '
     *           //'for the next iteration.'
      endif
      endif
c DIIS!
      call secund(t0)
      call elapsec(et0)
      if (diis.eq.1.or.(dabs(xmaxele).lt.5.0d-2.and.
     *                   ifirst_diis_iter.le.iteration)) then
      if (diis.ne.1.and.iprnt.gt.2) 
     *         write(6,*)'Calculation of diis matrix started.'
      diis=1
      call sync_barrier
      call para_bcast(diis,TxCCSDInit)
      if (iprnt.gt.2) write(6,*) 'Sent diis info: ',diis
      call CCDiis(i,       j,          'calculate',  ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'ccc',bl(ilist),
     *            bl(ipairimages),nbf)
      else
      iwork=0
      call sync_barrier
      call para_bcast(iwork,TxCCSDInit)
      if (iprnt.gt.2) write(6,*) 'Sent diis info: ',iwork
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
      if (iprnt.gt.2) then
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
      if (iprnt.gt.1) call print_CCSD_results(ccsdela,ccsdcpu,itim_no)
      write(6,3232) iteration,energyr,DeltaE1,xmaxele,square,
     *              ccsdela(12)/6d1
      call flush(6)
c
c     STOP 'xx'
      call sync_barrier
      if (.false..and.omit.and.DeltaE.le.1.0d-7) then  ! DISABLED
        call filllist(bl(ilist),npairs)
        write(6,*)'All pairs are switched on!'
        omit=.false.
      call CCDiis(i,       j,          'reset',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
        goto 666
      endif
      if ((xmaxele.gt.Wavethresh .or. DeltaE.gt.Ethresh).and..not.mp3
     *     .and..not.do_mp4) then
        if (iteration.ge.maxiter) then
          write(6,*)'            * * * * Maximum number of'//
     *            ' iterations reached * * * *'
          write(6,*)'            * * * *             NO '//
     *            'CONVERGENCE           * * * *'
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
     *               idimen,   icano,     iepsi,  ifockAO,ifockMO,
     *               icorefock,ilist,     inumber,singles,DeltaE)
      endif
      call secund(t1)
      call elapsec(et1)
      if (iprnt.gt.2) then
      write(6,'(A,F10.1)')'Dump CPU time:   ',t1-t0        ! dump
      write(6,'(A,F10.1)')'Dump elaps time: ',et1-et0      ! dump
      endif
 70   continue
      if (do_mp4) then 
        call matzero('work1')
        do i=1,nval
          do j=1,i
            call CoefWrite(i,j,iwork1)
          enddo
        enddo
      call CoefInit_master(npairs,ncf,gauss_seidel,nodisk,   ccsd,
     *                     vorb,  nmo,ivalpair,    ivirpair, ifp,
     *                     nsym,  nval,nvirt,      ij_unique,isympairs,
     *                     l_restart,ndisk_dump,nbf)
      call CCDiis(i,       j,          'reset',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
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
ccs      call sleep(5)
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
#ifdef DUMP_AMPLITUDES
      irecc=8*(idimen*idimen+3)
c     OPEN(UNIT=99,FILE='amplitudes',FORM='UNFORMATTED',
c    *     ACCESS='DIRECT',RECL=irecc)
      ij=npairs
      call dynamic_getmem(idimen*idimen,icoef)
      do i=1,nval
        do j=1,i
          ij=ij+1
          call CoefRead('tt',i,j,icoef)
          call tmp_writ(i,j,ij,bl(iepairr+ij-npairs),bl(icoef),idimen)
        enddo
      enddo
      call dynamic_retmem(1)
#endif
c*******************************************************************
      call fafTerminate(iresult)
      call allclose_and_delete(.true.)
c
      call dynamic_retmark
      call matremark
      call retmark
c---------------------------------------------------------------------
      end
c=======================================================================
      subroutine tmp_writ(i,j,ij,epairr,matrix,isize)
      implicit none
      integer i,j,ij,isize
      real*8 epairr,matrix(isize,isize)
      write(99,rec=ij) i,j,epairr,matrix
      end
c=======================================================================
      subroutine afwritebin(ndisk, lbin, ibin4, ibin1, bin8, indxbin, 
     1                      irecij,  ij, islvid,byt8,irec_no,irec_last)
c    This routine writes a full bin to an array file
c Parameters:
c   ndisk   - the logical file unit number used here
c   lbin    - bin length
c   ibin4   - the integer*4 one ij bin
c   ibin1   - the integer*1 one ij bin overflow
c   bin8    - the double precision one ij bin
c   indxbin - the bin indices for ij pair
c   irecij  - the number of records written for this ij pair on this slave
c   ij      - the ij pair number
c   islvid  - the slave ID number, this is group number
c   byt8    - the logical variable, if bin8=.true. double precision bin is 
c             used, otherwise the ibin4 + ibin1 pair
c
        use memory
        use newpara
        implicit none
        integer ndisk,lbin,irecij,ij,islvid
        integer*4 ibin4(lbin)
        integer*1 ibin1(lbin)
        real*8    bin8(lbin)
        integer*2 indxbin(2,lbin)
        integer afindex(4),indexx
        integer irec_no,irec_last,istatus,info,mu,nu,icount
        logical byt8
        if (irec_no>irec_last .or. irec_no<0) then
          call para_initsend
          call para_send_pack(0,CCSD_Rec_Req)
          call para_recv_int(irec_no,1,0,CCSD_Rec_Off)
          call para_recv_int(irec_last,1,0,CCSD_Rec_Off)
        endif
        if (byt8) then
         call packer(3,indexx,bl(1),bin8,lbin*8,indxbin,lbin*4,irecij,8)
         irecij=irec_no
         call fafwrite(ndisk,bl(indexx),1,12*lbin+8,1,irecij,istatus)
        else
         STOP 'Not ready yet'
         afindex(4)=0
         call fafwrite(ndisk,ibin4,4,lbin,4,afindex,istatus)
         afindex(4)=1 
         call fafwrite(ndisk,ibin1,1,lbin,4,afindex,istatus)
         afindex(4)=2 
         call fafwrite(ndisk,indxbin,2,lbin*2,4,afindex,istatus)
        endif
        irec_no=irec_no+1
      end
!=======================================================================
c=======================================================================
      subroutine afwritebingr(ndisk, lbin, ibin4, ibin1, bin8, indxbin, 
     1                    irecij,ij,islvid,byt8,isize,irec_no,irec_last)
c    This routine writes a full bin to an array file
c Parameters:
c   ndisk   - the logical file unit number used here
c   lbin    - bin length
c   ibin4   - the integer*4 one ij bunch of bins
c   ibin1   - the integer*1 one ij bunch of bins overflow
c   bin8    - the double precision one ij bunch of bins
c   indxbin - the bin indices for ij pair
c   irecij  - the number of records written for this ij bunch of pairs on 
c             this slave
c   ij      - the granule number
c   islvid  - the slave ID number, this is group number
c   byt8    - the logical variable, if bin8=.true. double precision bin is 
c             used, otherwise the ibin4 + ibin1 pair
c   isize   - the amount of ij pairs written at once
c
        use memory
        use newpara
        implicit none
        integer ndisk,lbin,irecij,ij,islvid,isize,irec_no,irec_last
        integer*4 ibin4(lbin*isize)
        integer*1 ibin1(lbin*isize)
        real*8     bin8(lbin*isize)
        integer*2 indxbin(2,lbin*isize)
        integer afindex(4)
        logical byt8
c
        integer istatus,info,indexx
c
        if (irec_no>irec_last .or. irec_no<0) then
          call para_initsend
          call para_send_pack(0,CCSD_Rec_Req)
          call para_recv_int(irec_no,1,0,CCSD_Rec_Off)
          call para_recv_int(irec_last,1,0,CCSD_Rec_Off)
        endif
        if (byt8) then
         call packer(3,indexx,bl(1),bin8,isize*lbin*8,indxbin,
     1               isize*lbin*4,irecij,8)
         irecij=irec_no
         call fafwrite(ndisk,bl(indexx),1,12*isize*lbin+8,1,irecij,
     1                 istatus)
        else
          STOP 'Not ready yet'
          afindex(4)=0
          call fafwrite(ndisk,ibin4,4,lbin*isize,4,afindex,istatus)
          afindex(4)=1 
          call fafwrite(ndisk,ibin1,1,lbin*isize,4,afindex,istatus)
          afindex(4)=2 
          call fafwrite(ndisk,indxbin,2,lbin*2*isize,4,afindex,istatus)
        endif
        irec_no=irec_no+1
      end
!=======================================================================
      subroutine master_Gen(ncs,     ncf,      ictr,   nval,  nmo,
     *                      nfirst,  nlast,    thresh, core,  xmaxdisk,
     *                      ndisk,   irecadr,  lrec,   lbin,  int_kind,
     *                      nodisk,  byt8,     small,  vorb,  ndisktr,
     *                      list,    isympairs,ij_unique,iprnt,nbf)
c This subroutine calculates all necessary integrals needed for CCSD and
c transforms them. 
c INPUT:
c ncs      - number of contracted shells
c ncf      - number of contracted basis functions
c ictr     - pointer do contraction info etc.
c nval     - number of correlating occupied orbitals.
c nmo      - number of occupied orbitals.
c nfirst   - first correlating orbital.
c nlast    - last  correlating orbital.
c thresh   - integral threshold
c core     - this parameter is used for determination of number of
c            correlating orbitals. Orbitals with energy below core are frozen
c xmaxdisk - not used at present, but it is kept for future use
c irecadr  - irecadr contains the number of records for the ij bin,
c            written on each slave
c            where ij is a pair number. This is in fact two-dimensional
c            table, irecadr(npairs,nslv) where nslv is number of slaves
c int_kind - it is character string indicating the kind of calculated
c            integrals. May be: 'c' coulomb, 'x' exchange, 'e' EEO, 
c            'tc' - three external coulomb, 'tx' - three external exchange
c nodisk   - the parameter indicating if the calculation is performed
c            with all amplitudes kept in RAM memory. It means that
c            symmetrized and antisymmetrized amplitudes will be kept in RAM in
c            addition to usual amplitudes (EEO only)
c byt8     - logical parameter, if .true. all bins are writen to disk in 
c            real*8 precision, otherwise in 5 byte compact format
c small    - how to perform calculations, using symm anti-method or
c            using full matrices of K integrals and amplitudes? (EEO only)
c            usually=.true., false was used for experimantal purposes,
c            may be useful for very lage molecules with sparsity
c vorb     - the logical indicating if intergals will transformed
c            from AO to MO - virt. orbitals (second transformation)
c list     - integer array, for each ij pair the integer number is provided
c            1 = pair is strong, 0 = pair is weak, -1 = pair is distant
c OUTPUT:
c lrec     - the length of bin record
c lbin     - the length of bin
c ndisktr  - contains integrals extracted from bins, if vorb=.true. the
c            integrals are already transformed, otherwise there are
c            rewritten only  from integrals in the bins. But integrals
c            in nodisk still are needed if vorb=.true. because some
c            parts of program requires integrals in AO form.
c            This should be changed, if vorb=.false. the bins should be
c            deleted because are not used at all.
c ndisk    - number of disk unit where bins are stored. This in really
c            output value
      use memory
      use newpara
      use kinds
      use ccounters
      implicit real*8 (a-h,o-z)
      integer nbf
      logical LastSymPair,emp2only,nodisk,keepnodisk,byt8,small,vorb,af
      logical memory_success,igran_calculate,split
      dimension xintxx(9)
      integer, dimension(*) :: list
      common /job/jobname,lenJ
      parameter (max_merg=(MAX_MERGE))  ! fragile!!
      common /indexsaver/ik,ikindex(2,max_merg) ! for EEOGenerator
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*3 ch3
      character scftype*11
      character*256 jobname,scrfile,filename
      logical nofr,exst,restrt,dualbasis,smal 
      character*(*) int_kind
      character*2 int_cp
      dimension xnmo(2)
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
ctj header end
      dimension irecadr(*)
      character*256 filname1
      dimension Ipass(2,28), Jpass(2,28) !28 is 28 comp.of cartisian I-function
      parameter (i_orb_size=28)
      common /zero_in_EEO/ zeroed,zeroed_max
      parameter (ixint_time=20)
      real*8 xint_time(ixint_time)
      equivalence (tott,xint_time(ixint_time))
      equivalence (tot, xint_time(ixint_time-1))
      equivalence (skip,xint_time(ixint_time-2))
      character(len=100) :: err_msg
      integer(kind=i_4) :: inf
      integer isympairs(*),ij_unique
      integer :: icompr = 25
      integer*4 itest_tag
      logical image
c
ctj Int_kind=c coulomb, Int_kind=x exchange, Int_kind=e EEO
ctj -- BEGIN get some common variables
      strace=0.0d0
      sconstr=0.0d0
      sort=0.0d0
      extr=0.0d0
      const=0.0d0
      c1=0.0d0
      c2=0.0d0
      etrace=0.0d0
      econstr=0.0d0
      eort=0.0d0
      eetr=0.0d0
      eonst=0.0d0
      e1=0.0d0
      e2=0.0d0
      call elapsec(total_elaps0)
      call secund(total_cpu0)
      keepnodisk=nodisk
      call getival('iout',iout)
      smal=.TRUE.
ctj -- END   get some common variables
  46  format(/72('='))
      if (iprnt.gt.1) then
      write(iout,46)
      write(91,46)
      if (Int_kind.eq.'c') then
      write(iout,*) '                            The Jij Module  '
      write(91,*)   '                            The Jij Module  '
      else if (Int_kind.eq.'x') then
      write(iout,*) '                            The Kij Module  '
      write(91,*)   '                            The Kij Module  '
      else if (Int_kind.eq.'e') then
      write(iout,*) '                            The EEOij Module  '
      write(91,*)   '                            The EEOij Module  '
      else if (Int_kind.eq.'tx') then
      write(iout,*) '                         The TEIO exch ij Module  '
      write(91,*)   '                         The TEIO exch ij Module  '
      else if (Int_kind.eq.'tc') then
      write(iout,*) '                         The TEIO coul ij Module  '
      write(91,*)   '                         The TEIO coul ij Module  '
      else if (Int_kind.eq.'tt') then
      write(iout,*) '                   Three external int. for triples'
      write(91,*)   '                   Three external int. for triples'
      else
        call nerror(1,'GenCoulExInt',
     1  'Int_kind does not indicate any kind of integrals',0,0)
      endif
      write(iout,*)' '
      endif
c     call getival('nslv',nslv)
c-----------------------------------------------------------
      ncs_sm=ncs
      ncf_sm=ncf
c-----------------------------------------------------------
c This is for use in EEO Generator and for bins construction when we
c decide how much memory we need to reserve
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ncf*(ncf-1)/2
      nvirt=nbf-nmo
c-----------------------------------------------------------
c  put down a memory marker
c
      call mmark
      call matmark
      call dynamic_mmark
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
c SLAVE!
      call dynamic_matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      idics=mataddr('dsmx')
      call dynamic_getmem(ncs*ncs,ischwarz)
#if 0
      call memo1_int(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c
c get memory for screening density
c
      np4=4
c-------------------------------------------------------
c initialize the two-el. integral program
c
      thint=thresh
      iforwhat=5
c     call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
c    *             ax,     nrad,   nang,   lrad,   lang,
c    *             Iradq,  NBatch, .true., nmo,    0,
c    *             scftype,xintxx, nblocks,maxbuffer,maxlabels)
c                     ^       ^       ^       ^          ^
c OUTPUT:          NOT USED PRICES NO. BL. INT. BUF.  LABELS
c                          (ARRAY)           SIZE      SIZE
c-------------------------------------------------------
c
#endif
      call mmark
      call dynamic_mmark
c---------------------------------------------------------------------
c WARNING! Screening is turned off at present!
c.................................................
      if(ncore.gt.0) call matsub('occa','cano',1,nmo)
      call  matsub('occu','cano',nfirst,nlast)
      call  matsub('virt','cano',nmo+1,nbf)
      ioccu=mataddr('occu')
c.................................................
c SLAVE!
#if 0
      call  fill1(ncs,bl(idics))
c
c---------------------------------------------------------------------
c
      icano=mataddr('cano')
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
c   reserve space for one AO integral matrix
#endif
      call dynamic_matdef('xmat','q',ncf,ncf)
      ixadr=mataddr('xmat')
c space for one half-transformed coulomb or exchange operator
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
        call dynamic_matdef('halftra','q',nval,nval)
        ihalftra=mataddr('halftra')
      else if (int_kind.eq.'tt') then
        call dynamic_matdef('halftra','r',nvirt,nval)
        ihalftra=mataddr('halftra')
      else
        call dynamic_matdef('halfmulam','q',nval,nval)
        ihalfmulam=mataddr('halfmulam')
        call dynamic_matdef('halflammu','q',nval,nval)
        ihalflammu=mataddr('halflammu')
      endif
c  put down a memory marker
c#endif
      call mmark
      call dynamic_mmark
c
      write(91,*) 'Mark 1'
      call flush(6)
c SLAVE!
      thint=thresh
      istat=igetival('istat')
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
c  nrec is the total number of records written on all files
c  irec is the current counter in the file
      nrec=0
      irec=0
c---------------------------------------------------
c
      call secund(tt)
      if(iprnt.gt.4) write(iout,*) 'Startup time=',tt-tt0
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
c SLAVE!
      write(91,*) 'Mark 2'
      call flush(6)
#if 0
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
      call matrix_int_zeroit(bl(lzero),ncf)
C
CTJ
#endif
      npairs=(nval**2+nval)/2
#if 1
#if 1
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      call dynamic_getmem(npairs/intsize+1,icountr) !memory for bin counters
      call izeroit(bl(icountr),npairs)         !zero out the bin counters
      else if (int_kind.eq.'tt') then
      call dynamic_getmem(nval*nvirt/intsize,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nvirt)         !zero out the bin counters
      else
      call dynamic_getmem(nval*nval/intsize+1,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nval)         !zero out the bin counters
      endif
#endif
      write(91,*) 'Mark 3'
      call flush(6)
!     Determine how much memory we can allocate to the bins.  They should
!     be as big as possible but leave space for other allocations.
      memreal=igetival('memr')
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memvirt=memtot
      memoccu=lastaddr-ioffset
      memavail=memvirt-memoccu
      call dynamic_show_free(memavail)
      if (iprnt.gt.2) write(6,*)"Memory free: ",memavail
c     temporary changed, dynamic memory:
c  memory for the bins is the minimum of the available virtual and real memory
c  size, minus 1 MW for miscellaneous - the space needed for matrix multiplies etc.
      write(91,*) 'Mark 3a'
      call flush(6)
      maxshell=49 !max contracted shell  squared
cAF memory for first half transform
      if (byt8) then
        lb_multiplier=12
      else
        lb_multiplier=9
      endif
 432  continue
      if (maxshell.lt.10) call nerror(4,'GenCoulExInt',
     1  'Cannot further reduce the maxshell parameter',
     2  memleft,maxshell)
      mem1trans=npairs+nval**2+nval**2/8+3*ncf+ncf**2+
     1          ncf*ncf*maxshell+nval**2
      mem2trans=0
      memleft=memavail-mem1trans-1 500 000
c if amount of avail mem leads to lbon<100:
      if (memleft.lt.lb_multiplier*npairs*100/8) then
        maxshell=maxshell*0.8
        goto 432
      endif
      memavail=memleft
c In the case with EEO it is not desirable to allocate as much memory
c for bins as possible. There will be only one write of result of size
c npairs*ncf*ncf. But, we will have to read many times the amplitudes,
c having the same size, in the icf kcf loop. For small shells merging it
c means that we may read them over 1000 times. The time spent in the bins 
c writing in that case will be negligible.
      write(91,*) 'Mark 3b'
      call flush(6)
      if (Int_kind.eq.'e') then
      if (iprnt.gt.2)
     *  write(6,*) 'Reducing available memory for bins in EEO'//
     *             ' operator procedure'
        memleft =0.3d0*memleft
      endif
c  lbin is in units of 9-bytes (2 indices           - integer*2
c                               compressed integral - integer*4
c                               precision overflow  - integer*1)
      write(91,*) 'Mark 3c'
      call flush(6)
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      lbin=8*memleft/(lb_multiplier*npairs) !Size of bins in 8 byte words
      else if (int_kind.eq.'tt') then
      lbin=8*memleft/(lb_multiplier*nval*nvirt)!Size of bins in 8 byte words
      else
      lbin=8*memleft/(lb_multiplier*nval*nval) !Size of bins in 8 byte words
      endif
      write(91,*) 'Mark 3d'
      call flush(6)
      if(lbin.lt.50) then
        call nerror(4,'GenCoulExInt',
     1  'memory available for bins leads to bin size < 50, mem,lbin',
     2  memleft,lbin)
      end if
c  lbin is the length of a bin in 9-byte words
      if(lbin.gt.ncf*(ncf+1)/4) then ! speed up is achieved from writing
c                                      at different moments, if
c                                      lbin=ncf*(ncf+1) all slaves try
c                                      to write at the same time via
c                                      writeleftovers subroutine and 
c                                      efficiency is much lower, for
c                                      glucose the slowdown is from 7/9 to 4/9
c
        lbin=ncf*(ncf+1)/5
      end if
      imax_rec=2**21    ! record is limited to 2 MB
      if (lb_multiplier*lbin.gt.imax_rec) then
        lbin=imax_rec/lb_multiplier
      endif
c Send slaves bin length:
      write(91,*) 'Mark 4'
      call flush(6)
      call para_bcast(lbin,TxCCSDJob)
      write(91,*) 'Mark 5'
      call flush(6)
c  How much memory left?
#if 1
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      indexmem=npairs*lbin/i4size+1
      if (byt8) then
      memavail=memavail-npairs*lbin
      memavail=memavail-indexmem
      else
      int4mem=npairs*lbin/i4size+1
      int1mem=npairs*lbin/i1size+1
      memavail=memavail-int4mem
      memavail=memavail-int1mem
      memavail=memavail-indexmem
      endif
      if (iprnt.gt.2)
     *write(6,*)'mem,npairs,lbin,nbins',memavail,npairs,lbin,
     *                                  ncf*ncf/lbin
      call flush(6)
      else if (int_kind.eq.'tt') then
      indexmem=nval*nvirt*lbin/i4size+1
      if (byt8) then
      memavail=memavail-nval*nvirt*lbin
      memavail=memavail-indexmem
      else
      int4mem=nval*nvirt*lbin/i4size+1
      int1mem=nval*nvirt*lbin/i1size+1
      memavail=memavail-int4mem
      memavail=memavail-int1mem
      memavail=memavail-indexmem
      endif
      if (iprnt.gt.2)
     *write(6,*)'mem,npairs,lbin,nbins',memavail,npairs,lbin,
     *                                  ncf*ncf/lbin
      call flush(6)
      else !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      indexmem=nval*nval*lbin/i4size+1
      if (byt8) then
      memavail=memavail-nval*nval*lbin
      memavail=memavail-indexmem
      else
      int4mem=nval*nval*lbin/i4size+1
      int1mem=nval*nval*lbin/i1size+1
      memavail=memavail-int4mem
      memavail=memavail-int1mem
      memavail=memavail-indexmem
      endif
      if (iprnt.gt.2)
     *  write(6,*)'mem,npairs,lbin,nbins',memavail,npairs,lbin,
     *                                    ncf*ncf/lbin
      call flush(6)
      endif !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
#endif
c
c
c .................................................
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting length
      if (Int_kind.eq.'c') then
        filname1=scrfile(1:len)//'.coul.bins'
        len1 = len+10
        call fafClosem(ndisk,-1,iresult)
        call fafCreatem(filname1(1:len1)//char(0),ndisk,12*lbin+8)
        call sync_barrier
      else if (Int_kind.eq.'x') then
        filname1=scrfile(1:len)//'.exc.bins'
        len1 = len+9
        call fafClosem(ndisk,-1,iresult)
        call fafCreatem(filname1(1:len1)//char(0),ndisk,12*lbin+8)
        call sync_barrier
      else if (Int_kind.eq.'e') then
c ====granules add
      ioccvirt=npairs ! number of bins
      call dynamic_getmem((ioccvirt/i4size+1)*1,ipair2gran)
      call dynamic_getmem((ioccvirt/i4size+1)*2,igran2pair)
      call granules_calculate(ioccvirt,ncf,byt8,lbin,bl(igran2pair),
     *                        bl(ipair2gran),igranules,igranulesize,
     *                        nslv)
c ====granules stop
c
        ieeocounter=ieeocounter+1
        filname1=scrfile(1:len)//'.eeo.bins'
        len1 = len+9
        if (ieeocounter.gt.1) then 
          call fafClosem(ndisk,0,iresult)
        else
          call fafClosem(ndisk,-1,iresult)
        endif
        call fafCreatem(filname1(1:len1)//char(0),ndisk,12*lbin+8)
        call sync_barrier
      else if (Int_kind.eq.'tx') then
      txcounter=txcounter+1
      filname1=scrfile(1:len)//'.tx.bins'
      len1 = len+8
      if (txcounter.gt.1) then 
              call fafClosem(ndisk,0,iresult)
       else
          call fafClosem(ndisk,-1,iresult)
       endif
      call fafCreatem(filname1(1:len1)//char(0),ndisk,12*lbin+8)
        call sync_barrier
      else if (Int_kind.eq.'tc') then
      tccounter=tccounter+1
      filname1=scrfile(1:len)//'.tc.bins'
      len1 = len+8
      if (tccounter.gt.1) then 
              call fafClosem(ndisk,0,iresult)
       else
          call fafClosem(ndisk,-1,iresult)
       endif
      call fafCreatem(filname1(1:len1)//char(0),ndisk,12*lbin+8)
        call sync_barrier
      else if (Int_kind.eq.'tt') then
c ====granules add
      ioccvirt=nval*nvirt ! number of bins
      call dynamic_getmem((ioccvirt/i4size+1)*1,ipair2gran)
      call dynamic_getmem((ioccvirt/i4size+1)*2,igran2pair)
      call granules_calculate(ioccvirt,ncf,byt8,lbin,bl(igran2pair),
     *                        bl(ipair2gran),igranules,igranulesize,
     *                        nslv)
c ====granules stop
      ttcounter=ttcounter+1
      filname1=scrfile(1:len)//'.tt.bins'
      len1 = len+8
      if (ttcounter.gt.1) then 
              call fafClosem(ndisk,0,iresult)
       else
          call fafClosem(ndisk,-1,iresult)
       endif
      call fafCreatem(filname1(1:len1)//char(0),ndisk,
     *               12*lbin*igranulesize+8)
        call sync_barrier
      endif
c Send slaves ndisk number:
      call para_initsend
      call para_pack_int(ndisk,1)
c and granules:
      if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
      call para_pack_int(igranules,1)
      call para_pack_int(igranulesize,1)
      call para_pack_int(ioccvirt,1)
      call para_pack_int4(bl(igran2pair),2*ioccvirt)
      call para_pack_int4(bl(ipair2gran),1*ioccvirt)
      endif
      call para_bcast_pack(TxCCSDJob)
#endif
      write(91,*) 'Mark 6'
      call flush(6)
c--------------------------------------------------------------------
      icont=0
c TEST!!
c     call dynamic_show_free(memavail)
c--------------------------------------------------------------------
c
      split=.true. ! shells splitting and merging
      if (Int_kind.eq.'e') then
        max_merg_shell_size=keepmax
        if (small) then
        if (.not.nodisk) then
          write(6,*)'Memory available before max_merg_shell_size calc.',
     *             memavail
C allocation for experimental shells merging
c Things which are commented are moved up
          memavail=memavail-30*ncf*ncf ! see above - amplitudes
          max_merg_shell_size=memavail/(ncf*ncf+2*npairs*icompr)
          if (max_merg_shell_size.lt.keepmax*3) then
           write(6,*) 'Splitted max_merg_shell_size would be: ',
     *                            max_merg_shell_size
           split=.false. ! shells merging without splitting
           memavail=memavail+30*ncf*ncf ! Fragile!
           max_merg_shell_size=memavail/(ncf*ncf+2*npairs*icompr)
          write(6,*)'But first try no-splited max_merg_shell_size is: ',
     *                            max_merg_shell_size
          endif
          if (max_merg_shell_size.lt.keepmax*4) then ! if still...
            memavail=memavail+24*ncf*ncf  ! decrease memory for
c                           integral buffer, increases INTEGRAL splitting.
           split=.false. ! shells merging without splitting
           max_merg_shell_size=memavail/(ncf*ncf+2*npairs*icompr)
           write(6,*) 'Second try non-splited max_merg_shell_size is: ',
     *                            max_merg_shell_size
          endif
          call flush(6)
          if (max_merg_shell_size.lt.maxshell) then
            write(6,*)'Estimated number of max_merg_shell_size was:',
     *               max_merg_shell_size
            STOP 'Too less memory for EEO operators'
          endif
          if (max_merg_shell_size.gt.max_merg) then
            max_merg_shell_size=max_merg
          endif
        endif
        write(6,*)'Splitting of shells: ',split
c
        call dynamic_matdef('screen','q',ncf,ncf)
        iscr=mataddr('screen')
        call matzero('screen')
        call GenAmplSymm_Mas(npairs, ncf,       nval,  nodisk,   nfirst,
     *                       nmo,    vorb,     .true., 'master', ndisks,
     *                       ndiska, bl(iscr), list,isympairs,ij_unique,
     *                       nslv)
        call build_idics(bl(iscr),bl(idics),ncf,ncs,bl(ictr))!wrong at present!
c
c Send neccessary data to slaves:
c
        call para_initsend
        call para_pack_int(max_merg_shell_size,1)
c       call para_pack_int(ndisks,1)
c       call para_pack_int(ndiska,1)
        call para_pack_int(split,1)
        call para_pack_real(bl(iscr),ncf*ncf)
        call para_pack_real(bl(idics),ncs*ncs)
        call para_bcast_pack(TxCCSDJob)
        call dynamic_matrem('screen')
c
c       call dynamic_getmem(ntrisymm*max_merg_shell_size,ibasesymm)
c       call dynamic_getmem(ntrianti*max_merg_shell_size,ibaseanti)
        write(6,*)'max_merg_shell_size=  ',max_merg_shell_size
        if (.not.nodisk)
     *write(6,*)'amplitudes will be read',ncf*ncf/2/max_merg_shell_size,
     *         '  times'
        call flush(6)
      else ! if (small) then - no shell merging
        ndiskampl=46
        ndiskresid=47
        filname1=scrfile(1:len)//'.eeo.ampl'
        len1 = len+9
        ireclength=ncf*ncf*8
        if (ieeocounter.eq.1) then 
        OPEN (UNIT=ndiskampl,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL=ireclength)
        filname1=scrfile(1:len)//'.eeo.resi'
        len1 = len+9
        OPEN (UNIT=ndiskresid,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL=ireclength)
        endif
        call AmplitudesGenerator(ncf,nval,ndiskampl,ndiskresid)
        memavail=memavail-ncf*ncf*maxshell ! for integral store
        memavail=memavail-1 500 000 ! for internal reserv in integr.  proc.
        nmatrices=memavail/(ncf*ncf)
            memory_success=.false.
#if 0
c  May be useful in above symm and antisymm algorithm
        if (nmatrices.gt.(nval*nval+max_merg_shell_size)) then
          do
          memory_used=(max_merg_shell_size+nval*nval)*ncf*ncf+
     *                max_merg_shell_size*nval*nval
          if (memory_used.gt.memavail) then
            max_merg_shell_size=max_merg_shell_size-1
            write(6,*) 'Reducing max_merg_shell_size which is: '
     *                                              max_merg_shell_size
            memory_success=.false.
          else
            memory_success=.true.
            exit
          endif
          enddo
            if (max_merg_shell_size.lt.(keepmax-10))
     *                                           memory_success=.false.
            max_amplitudes_size=nval*nval
            if (memory_success) write(6,*) 'All amplit. in memory.'
        endif
#endif
        if (.not.memory_success) then
          max_amplitudes_size=keepmax
          max_merg_shell_size=(nmatrices-max_amplitudes_size)/2 ! I will
c         reserve memory inside EEO subroutine for squeezed copy of
c         integrals, that's why divide by 2
          if (max_merg_shell_size.gt.ncf*(ncf+1)/2)
     *                         max_merg_shell_size=ncf*(ncf+1)/2
          do
          memory_used=(max_merg_shell_size+max_amplitudes_size)*ncf*ncf+
     *                max_merg_shell_size*nval*nval
          if (memory_used.ge.memavail) then
            max_merg_shell_size=max_merg_shell_size-1
c           write(6,*) 'Reducing max_merg_shell_size which is: ',
c    *                                              max_merg_shell_size
            if (max_merg_shell_size.le.1) STOP 'Too little memory EEO'
          else
            exit
          endif
          enddo
        endif
        write(6,*) 'max_merg_shell_size: ',max_merg_shell_size
        write(6,*) 'max_amplitudes_size: ',max_amplitudes_size
        write(6,*) 'Amplitudes will be read: ',
     *                       ncf*ncf/2/max_merg_shell_size,' times.'
        if (max_merg_shell_size.gt.max_merg)
     *                                  max_merg_shell_size=max_merg
        call dynamic_getmem(ncf*ncf*max_merg_shell_size,ishellsaddr)
        call dynamic_getmem(ncf*ncf*max_amplitudes_size,iamplitaddr)
        call dynamic_getmem(max_merg_shell_size*nval*nval,iresultaddr)
      endif ! if (small) then
      endif ! if (Int_kind.eq.'e') then
c
      do i=1,ixint_time
        xint_time(i)=0.0d0
      enddo
C
      write(91,*) 'Mark 7'
      call flush(6)
      ncf2=ncf*ncf
      nskipped=0
      info=0
      istep=max(15,npairs/4)
      istart_rec=1
      istop_rec=istart_rec+istep-1
      do ics=ncs,1,-1
         do jcs=ics,1,-1
         if (image(ics,jcs)) cycle
         do
           call para_recv_anypack(islgid,itest_tag)
           if (itest_tag.eq.TxCCSDReq) then
             call para_unpack_int(iready,1)
             call para_unpack_int(islvid,1)
             if (iready.eq.1) then
               call para_initsend
               iwork=1
               call para_pack_int(iwork,1)
               call para_pack_int(ics,1)
               call para_pack_int(jcs,1)
               call para_send_pack(islgid,TxCCSDJob)
             endif
             exit
           else if (itest_tag.eq.CCSD_Rec_Req) then
             call para_send(istart_rec,islgid,CCSD_Rec_Off)
             call para_send(istop_rec,islgid,CCSD_Rec_Off)
             istart_rec=istop_rec+1
             istop_rec=istart_rec+istep-1
           else
             STOP 'Error, unexpected message received'
           endif
         enddo
        end do     !  over kcs shell
      end do       !  over ics shell
      if (Int_kind.eq.'tx'.or.Int_kind.eq.'tc') then
        irecidx=nval*nval
      else if (Int_kind.eq.'e'.or.Int_kind.eq.'tt') then
        irecidx=igranules
      else
        irecidx=npairs
      endif

      call dynamic_unlock(irecadr,iii)
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
        do iii=1,npairs*nslv
        irecadr(iii)=0
        enddo
      else if (int_kind.eq.'tt') then
        do iii=1,nval*(nbf-nmo)*nslv
        irecadr(iii)=0
        enddo
      else
        do iii=1,nval*nval*nslv
        irecadr(iii)=0
        enddo
      endif


      icompl=0
      iresult=0
      do !tell slaves that this is the end of this part of job, collect
!         results
        call para_recv_anypack(islgid,itest_tag)
        if (itest_tag.eq.TxCCSDReq) then
          call para_unpack_int(iready,1)
          call para_unpack_int(islvid,1)
          if (iready.eq.1) then
            call para_initsend
            iwork=0
            call para_pack_int(iwork,1)
            call para_send_pack(islgid,TxCCSDJob)
            icompl=icompl+1
          endif
        else if (itest_tag.eq.CCSD_Rec_Req) then
          call para_send(istart_rec,islgid,CCSD_Rec_Off)
          call para_send(istop_rec,islgid,CCSD_Rec_Off)
          istart_rec=istop_rec+1
          istop_rec=istart_rec+istep-1
        else if (itest_tag.eq.TxCCSDRes) then
          call para_unpack_int(igid,1)
          if (igid.gt.nslv.or.igid.lt.1) then
            write(6,*) 'igid.gt.nslv.or.igid.lt.1 !!!',igid,nslv
            STOP 'ERROR1, see output.'
          endif
          call para_unpack_int(irecadr(irecidx*(igid-1)+1),irecidx)
          iresult=iresult+1
        else
          STOP 'Error, unexpected message received 1'
        endif
        if (iresult.eq.nslv.and.icompl.eq.nslv) exit
      enddo
      call sync_barrier
      if (Int_kind.eq.'e'.and.(.not.nodisk)) then!files where symm and antysymm 
c                                                 coefficients are written
c       close(78,STATUS='delete')
c       close(79,STATUS='delete')
        call fafClosem(ndisks,0,info)
        call fafClosem(ndiska,0,info)
      endif
      nodisk=keepnodisk
ctj We need to write down leftovers. For CC all coulomb integrals will
c   be used many times during calculation and we have to had them all
c   written on disk
      write(91,*) 'Mark 8'
      call flush(6)

c-------------------------------------------------------------------
c     if(nskipped.gt.0) then
c       write(iout,*) nskipped,
c    1   ' pairs were skipped because no integrals were calculated'
c     end if
c Added transformation here BEGIN
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                 ! 80 is starting lenght
      int_cp=int_kind
      if      (int_kind.eq.'x') then
        icounterx=icounterx+1
        filename=scrfile(1:len)//'.Kij'
        len1=len+4
      else if (int_kind.eq.'c') then
        icounterc=icounterc+1
        filename=scrfile(1:len)//'.Cij'
        len1=len+4
      else if (int_kind.eq.'e') then
        icountere=icountere+1
        filename=scrfile(1:len)//'.Eij'
        len1=len+4
      else if (int_kind.eq.'tx') then
        icountertx=icountertx+1
        filename=scrfile(1:len)//'.TXij'
        len1=len+5
      else if (int_kind.eq.'tc') then
        icountertc=icountertc+1
        filename=scrfile(1:len)//'.TCij'
        len1=len+5
      else if (int_kind.eq.'tt') then
        icountertt=icountertt+1
        filename=scrfile(1:len)//'.TTij'
        len1=len+5
      endif
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
      lrec=idimen*idimen*8
c     fileopen=.false.
      call fafreopen(ndisk)
c     inquire(ndisktr,opened=fileopen)
c     if (fileopen) STOP 'Error in subroutine GenCoul, ndisktr opened.'
cc    close(ndisktr,status='DELETE')
cc    OPEN (UNIT=ndisktr,FILE=filename(1:len1),FORM='UNFORMATTED',
cc   1      ACCESS='DIRECT',RECL=lrec)
      if (int_kind.eq.'x' .and.icounterx .gt.1 .or.
     *    int_kind.eq.'c' .and.icounterc .gt.1 .or.
c    *    int_kind.eq.'e' .and.icountere .gt.1 .or.
     *    int_kind.eq.'tx'.and.icountertx.gt.1 .or.
     *    int_kind.eq.'tc'.and.icountertc.gt.1 .or.
     *    int_kind.eq.'tt'.and.icountertt.gt.1) then
          call fafclosem(ndisktr,0,istatus)
      else
          call fafclosem(ndisktr,-1,istatus)
      endif
      if (int_kind.ne.'e' .or.icountere .eq.1) then
          call fafCreatem(filename(1:len1)//char(0),ndisktr,lrec)
      else
          call fafCreatem(filename(1:len1)//char(0),ndisktr,-1)
      endif
          call sync_barrier
c
      irecskip=irecidx
c     if (int_kind.eq.'tx'.or.int_kind.eq.'tc') irecskip=nval*nval
      islaves=nslv 
      af=.true.
c
c Initialize PVM
      call para_initsend
      call para_pack_int(irecadr,irecskip*islaves)
      call para_pack_int(ndisktr,1)
      call para_bcast_pack(TxCCSDInit)
      call dynamic_lock(irecadr,iii)
c Do job:
      iwork=1
      do ij=1,irecskip
        call para_recv_pack(islgid,TxCCSDReq)
        call para_unpack_int(iready,1)
        if (iready.eq.0) STOP 'Error slave.'
        call para_unpack_int(sltid,1)
        call para_initsend
        call para_pack_int(iwork,1)
        call para_pack_int(ij,1)
        call para_send_pack(islgid,TxCCSDJob)
      enddo
c Tell slaves job is done.
      iwork=0
      do i=1,nslv
      call para_recv_pack(islgid,TxCCSDReq)
      call para_unpack_int(iready,1)
      if (iready.eq.0) STOP 'Error slave.'
      call para_unpack_int(sltid,1)
      call para_send(iwork,islgid,TxCCSDJob)
      enddo
      call fafreopen(ndisktr)
      call sync_barrier
      tot_latency1=0.0d0
      tot_latency2=0.0d0
      total_cpu=0.0d0
      zeroed=0.0d0
      zeroed_max=0.0d0
      call elapsec(total_elaps1)
      call secund(total_cpu1)
      total_elaps=total_elaps1-total_elaps0
      total_cpu_mast=total_cpu1-total_cpu0
      inf=0
      info=0
      call para_reduce(total_cpu,1,TxCCSDRedu1)
      call para_barrier
      call para_reduce(tot_latency1,1,TxCCSDRedu2)
      call para_barrier
      call para_reduce(tot_latency2,1,TxCCSDRedu3)
      call para_barrier
      call para_reduce(zeroed,2,TxCCSDRedu1)
      call para_barrier
      xiare_zeros=0.0d0
      call para_reduce(xiare_zeros,1,TxCCSDRedu2)
      call para_barrier
      call para_reduce(strace,14,TxCCSDRedu3)
      call para_barrier
      call para_reduce(xint_time,ixint_time,TxCCSDRedu1)
      call para_barrier
      call para_bcast(iwork,TxCCSDQuit)
      totint=dble(ncf)*(dble(ncf)+1.0d0)/2.0d0*dble(ncf)*dble(ncf)
      if (iprnt.gt.2) then
      write(6,*)'Max latency for work=0 was:    ',tot_latency1
      write(6,*)'Max latency after work=0 was:  ',tot_latency2
      write(6,*)'Efficiency in this subroutine: ',
     *                          (total_cpu+total_cpu_mast)/total_elaps
      write(6,264)'Number & % of zeroed integrals:      ',zeroed,
     *                          zeroed/dble(totint)*100.0d0
      write(6,264)'Max number & % of zeroed integrals:  ',zeroed_max,
     *                          zeroed_max/dble(totint)*100.0d0
      write(6,264)'Number & % of skipped integrals:     ',
     *           xiare_zeros,xiare_zeros/dble(totint)*100.0d0
      endif
      if (Int_kind.eq.'e') then
      icompr=1/(1d0-zeroed/dble(totint))*2
      if (icompr.gt.max_merg/(max_merg_shell_size))
     *                       icompr=max_merg/(max_merg_shell_size)+1
      endif
c
 264  FORMAT (A,F20.0,F7.2)
c     common /timingstj/ strace,sconstr,sort,extr,const,c1,c2
c    *                   etrace,econstr,eort,eetr,eonst,e1,e2
      if (iprnt.gt.2) then
      write(6,345)
     *     'Matrix multiplication     CPU, ELA: ',strace,etrace,
     *     'Amplitudes reads          CPU, ELA: ',sconstr,econstr,
     *     'Sorting                   CPU, ELA: ',sort,eort,
     *     'Screening                 CPU, ELA: ',extr,eetr,
     *     'Symm-anti integrals build CPU, ELA: ',const,eonst,
     *     'Zero row seek             CPU, ELA: ',c1,e1,
     *     'Zero integrals            CPU, ELA: ',c2,e2
 345  FORMAT (7(A,2F15.2,/))
      write(6,*)'                      ************************'
      write(6,*)'total amount of integrals: ', tott
      write(6,*)'total would be calc: ',tot
      write(6,*)'skipped:             ',skip
      write(6,12)'Time spent in rys,zeroing: ',xint_time(1),
     *         '                  getmem:  ',xint_time(2),
     *         '                  int:     ',xint_time(3),
     *         '                  put res: ',xint_time(4),
     *         '                  irow,icol',xint_time(5),
     *         'Total time inside new int: ',xint_time(6)
 12   FORMAT (6(A,F15.4,/))
      call flush(6)
      write(6,*)'                      ************************'
      endif
cold  call dynamic_retmem(4)
c......................................................................
c timing here is irrelevant :
      call secund(txxx2)
      txxx0=txxx1
      call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
c
c......................................................................
      nrec=nrec+irec
      if(iprnt.gt.2) then
         write(iout,*) ' Total number of records written=',nrec
      endif
c
      if(iprnt.gt.2) then
c        Calculate statistics
         PerCent=100.0d0*ENonZero/(dble(mulam+mulamd)*ncf**2)
         write(iout,40) PerCent
      endif
c
  40  format(
     *'  Percentage of the AO matrix used in the transformation=',f7.3/)
c......................................................................
  50  continue
c  release memory
      call retmark
      call dynamic_retmark ! mmark just before int_coul
c Slave!
#if 0
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e'
     *                   .or.int_kind.eq.'tt') then
      call dynamic_matrem('halftra')
      else
      call dynamic_matrem('halflammu')
      call dynamic_matrem('halfmulam')
      endif
#endif
      call retmark
      call dynamic_retmark ! mmark after irec allocation
c......................................................................
#if 0
      write(iout,61)
   61 format('  CPU & Elapsed timings in the coulomb integrals module ')
      write(iout,*) '  '
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format('  Coulomb Integrals =',f8.2,' and ',f8.2,' minutes')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format('  First half transf. =',f8.2,' and ',f8.2,' minutes')
      call flush(iout)
#endif
c-----------------------------------------------------------
      call secund(tt1)
      call elapsec(elaps1)
      if(iprnt.ge.2) then
        write(iout,*) ' CPU and elapsed time for integrals+first half ',
     1  'transformation ',(tt1-tt0)/sixty,(elaps1-elaps0)/sixty
      end if
c..................................................
c---------------------------------------------------------------------
cslavecall symmon
c---------------------------------------------------------------------
      call dynamic_retmark
      call matremark
      call retmark
      if (int_kind.eq.'tt') call fafclosem(ndisk,0,istatus)
      end
c=======================================================================
      subroutine check_ijchunk_reduce(icalc_pairs,nval,ichunk,jchunk,
     *                                reduced)
c This simple subroutine reduces the sizes of bocks used for
c Y Z and Q computation. This is needed at the end of calculations
c for better parallel balancing. The price for it is that efficiency of
c calculation of such small blocks is MUCH lower. But the choice is either to 
c finish calculations with low efficiency due to small blocks or finish
c with low efficiency because some slaves still finish their job and other 
c slaves do not work at this time.
      implicit none
      integer icalc_pairs,nval,ichunk,jchunk
      logical reduced
      real*8 ratio
      if (reduced) return
      ratio=dble(icalc_pairs)/dble(nval*nval)
      if (ratio.gt.0.75d0) then
        ichunk=ichunk/3
        jchunk=jchunk/3
        if (ichunk.lt.1) ichunk=1
        if (jchunk.lt.1) jchunk=1
        reduced=.true.
      endif
      end
c=======================================================================
      subroutine granules_calculate(g_pairs,ncf,byt8,lbin,igran2pair,
     *                              ipair2gran,igranules,igranulesize,
     *                              nslv)
      implicit none
      integer g_pairs,ncf,lbin,igranules,igranulesize,nslv
      logical byt8
      integer*4 igran2pair(2,g_pairs),ipair2gran(g_pairs)
c
      integer lb_multiplier,islavemem,lrec,istart,istop,i,j
c
      if (byt8) then
        lb_multiplier=12
      else
        lb_multiplier=9
      endif
c
      call dynamic_show_free(islavemem)   ! Assesment of memory after 1st half
      islavemem=islavemem-3*ncf*ncf       ! just in case 3*ncf*ncf
      if (byt8) then
        igranulesize=islavemem/(lbin*3/2+1+ncf*ncf)
      else
        igranulesize=islavemem/(lbin*9/8+1+ncf*ncf)
      endif
      do
        igranules=g_pairs/igranulesize
        if (mod(g_pairs,igranulesize).ne.0) igranules=igranules+1
        if (igranules.ge.nslv*2) exit ! 2 granules/slave
        if (igranulesize.lt.1) call nerror(1,'MP2 module',
     $                      'Error in granule determination',0,0)
        if (igranulesize.eq.1) exit
        igranulesize=igranulesize-1
      enddo
      do
        lrec =lb_multiplier*lbin*igranulesize! record length in bytes
c Limit record size to 2 MB. This pevents overloading AF
c by large amount of large blocks passed via PVM, which may next cause RAM 
c memory overflow or swapping and efficiency loss.
        if (lrec.gt.2 097 152) then 
          if (igranulesize.eq.1) exit
          igranulesize=igranulesize-1
        else
          exit
        endif
        igranules=g_pairs/igranulesize
        if (mod(g_pairs,igranulesize).ne.0) igranules=igranules+1
      enddo
      call flush(6)
      call optimal_isize(igranulesize,g_pairs)
      igranules=g_pairs/igranulesize
      if (mod(g_pairs,igranulesize).ne.0) igranules=igranules+1
c     write(6,*)'igranules:    ',igranules
c     write(6,*)'igranulesize: ',igranulesize
      lrec =lb_multiplier*lbin*igranulesize! record length in bytes
c     write(6,*)'Records size: ',lrec,' bytes'
      call flush(6)
c build maps: granule -> ipairstart,ipairstop
c             ipair   -> granule
      istart=1
      do i=1,igranules
        istop=istart+igranulesize-1
        if (istop.gt.g_pairs) istop=g_pairs
        if (istart.gt.g_pairs) call nerror(3,'MP2 module',
     $   'Error in granule maps determination ',istart,g_pairs)
        igran2pair(1,i)=istart
        igran2pair(2,i)=istop
        do j=istart,istop
          ipair2gran(j)=i
        enddo
        istart=istop+1
      enddo
c
c       do i=1,igranules
c       write(6,*)'Granule: ',i,'   start: ',igran2pair(1,i),'  stop: ',
c    *          igran2pair(2,i)
c       enddo
c       call flush(6)
      end
c===============================================================================
      subroutine dump_data(dump,     ndisk_dump,ncf,    nval,   npairs,
     *                     idimen,   icano,     iepsi,  ifockAO,ifockMO,
     *                     icorefock,ilist,     inumber,singles,DeltaE)
      use memory
      use kinds
      use ccounters
      implicit none
      logical dump,singles
      integer ndisk_dump,ncf,idimen,icano,iepsi,ifockAO,ifockMO
      integer icorefock,ilist,inumber,nval,npairs
      real*8 DeltaE
c
      integer i,j,icoef,isingles
      character*256 :: filename = ' '
      integer :: len1 = 0
      logical open_st
c
      ic_dump_data=ic_dump_data+1
c
      if (dump) then ! * * * * * * * * * * * * DUMP FOR RESTART
      if (ic_dump_data.eq.1) then
        inquire(ndisk_dump,opened=open_st)
        if (.not.open_st) STOP 'File not opened'
        inquire(ndisk_dump,name=filename)
        close(ndisk_dump,status='keep')
        call rmblan(filename,80,len1)
      endif
      call flush(6)
      OPEN (UNIT=ndisk_dump,FILE=filename(1:len1),
     1      FORM='UNFORMATTED',ACCESS='SEQUENTIAL',STATUS='REPLACE')
c
        call general_write1(ndisk_dump,bl(icano),ncf*ncf*8)
        call general_write1(ndisk_dump,bl(iepsi),ncf*8)
        do i=1,nval
          do j=1,i
            call CoefRead('tt',i,j,icoef)
            call general_write1(ndisk_dump,bl(icoef),idimen*idimen*8)
          enddo
        enddo
        call general_write1(ndisk_dump,bl(ifockAO),idimen*idimen*8)
        call general_write1(ndisk_dump,bl(ifockMO),ncf*ncf*8)
        call general_write1(ndisk_dump,bl(icorefock),idimen*idimen*8)
        if (singles) then
        call pointersingles(isingles)
        call general_write1(ndisk_dump,bl(isingles),nval*idimen*8)
c       write(99,'(A)') '==========================================='
c       do i=0,nval*idimen-1
c       write(99,'(F25.16)') bl(isingles+i)
c       enddo
        endif
        call general_write1(ndisk_dump,bl(ilist),npairs*8/intsize)
        call general_write1(ndisk_dump,inumber,8/intsize)
        call general_write1(ndisk_dump,DeltaE,8)
      endif
      close(ndisk_dump,status='keep')
      end
c===============================================================================
      subroutine print_overlap(istart,istop,filename)
      use memory
      implicit none
      integer istart,istop
      character*(*) filename
c
      integer isize,iress,i,j,mataddr
c
      isize=istop-istart+1
      call matdef('ress','q',isize,isize)
      iress=mataddr('ress')
      call matsub('aa','cano',istart,istop)
      call matsimtr('overlap','aa','ress')
      OPEN(97,FILE=filename)
      do i=0,isize-1
        do j=0,isize-1
          write(97,'(2I6,F25.15)') i+1,j+1,bl(iress+i+j*isize)
        enddo
      enddo
      CLOSE(97)
      call matdisc('aa')
      call matrem('ress')
      end
c===============================================================================

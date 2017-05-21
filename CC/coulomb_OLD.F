#include "maxi.h"
      Subroutine coulomb
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
      parameter(nopt=20)
      parameter(nmethod=10)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $              filname3,filname4
      character*256 chopv(nopt)
      character*6 methods(nmethod),method
      character opnames*4,scftype*11,wvfnc*20,ch3*3
      logical nofr,exst,restrt,dualbasis,smal,mp3,cep0,cep2,mp2,af,mp4
      logical do_mp4
      character*2 int_kind
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension ifilestat(13)
      dimension xnmo(2)
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
      logical trans,signum,docansym
      character*80, parameter :: spc='                             '
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
      strace=0.0d0
      ik=0
      af=.false.
      islaves=1
      ndiske=-1
      lbine=-1
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
      call signal_default(11)
      call secund(tt0)
      call elapsec(elaps0)
      iseed = 0
      call srand(iseed)
      call reset_counters
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
      log_diis=.false.
      singles=.false.
      ccsd=.false.
      vorb=.true.
      cache=.false.
      variational=.true.
      integerT=.false.
      xmaxdisk=2 500 000 000.0d0      ! currently not active
c
c  basis function symmetry pairs are in inx(ifp)
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
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,wvfnc) ! xlow - lowest eig. of S
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
      If(wvfnc(1:3).EQ.'RHF' .or. wvfnc(1:3).EQ.'RMP') then
      else
         call nerror(1,'Correlated Energy Module',
     $   'Correlated Energy Requested but System is NOT RHF or MP2 !',
     $    0,0)
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
      if(ifound(2).gt.0) l_triples=.true.
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
      call dynamic_mmark
      call mmark
      call matmark
c
c -- PRINt
      if(ifound(11).gt.0) iprnt=iopv(1,11)
      call setival('ccprint',iprnt)
c
c -- ITER
      maxiter=20
      if(ifound(12).gt.0) maxiter=iopv(1,12)
c
c -- CHKP
c      not implemented in serial mode
c
c -- DUMP
c      not implemented in serial mode
c
c -- RESTart
c      not implemented in serial mode
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
c      not implemented in serial mode
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
          call nerror(2,'Correlated Energy Module',
     1   'Program cannot find the number of MOs on <jobname>.14',
     2    ni,ni)
        end if
      end if
      nmo=xnmo(1)  ! no of occupied orbitals
      erhf=xnmo(2) ! HF energy
c----------------------------------------------
c define some matrices,
      call dynamic_matdef('cano','q',ncf,ncf)
      call dynamic_matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
      if(ifound(17).eq.0) then
        call matread('cano',np4,'evec_rhf') ! MOs
        call matread('epsi',np4,'eval_rhf') ! HF eigenvalues
      else
        call matread('cano',np4,'loca_rhf') ! MOs
        call matread('epsi',np4,'eval_rhf') ! HF eigenvalues
        call dynamic_matdef('loc-epsi','d',nmo,nmo)
        call matread('loc-epsi',np4,'eloc_rhf') ! HF eigenvalues
        ilocepsi=mataddr('loc-epsi')
        call tfer(bl(ilocepsi),bl(iepsi),nmo)
        call dynamic_matrem('loc-epsi')
        loca=.true.
      end if
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
        write(6,'(A,I5)') 'DIIS expansion length: ',length_diis
        write(6,'(A,I5)') 'DIIS evaluation start: ',ifirst_diis_iter
      else
        write(6,'(A)')  'Calculation without DIIS'
      endif
c
      write(6,'(A,F6.1)')  'Level shift set to ',shift
      write(6,*)  'Threshold for energy convergence:', Ethresh
      write(6,*)  'Threshold for wvfn.  convergence:', Wavethresh
C
C  Determine the orbitals to be correlated
C
      If(ifound(6).gt.0) Then
        nfirst=iopv(1,6)
        if(nfirst.LT.1.or.nfirst.GT.nmo)
     $     call nerror(5,'Correlated Energy Module',
     $      'Lowest Orbital in range for Correlation Incorrect',0,0)
        nlast=nmo
        nval=nlast-nfirst+1
        ncore=nlast-nval
        nend=iopv(2,6)
        if(nend.LE.nmo.or.nend.GT.nbf)
     $     call nerror(6,'Correlated Energy Module',
     $      'Highest Orbital in range for Correlation Incorrect',0,0)
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
c
      call twoelinit
c
      call open_mp3_resid(ndisk_mp3,idimen,.false.,mp4)
c  print the data of the calculation early
      write(iout,'("  Number of contracted functions =",i8,/,
     1             "  Number of correlated orbitals  =",i8,/,
     2             "  Number of virtual orbitals     =",i8)')
     3                ncf,nval,nvirt
      write(iout,*) ' '
      call flush(iout)
c
c Check orbital symmetry:
      call dynamic_matdef('overlap','q',ncf,ncf)
      ioverlap=mataddr('overlap')
      call OverlapBuilder(ioverlap)
      call dynamic_unlock(bl(icano),i)
      call dynamic_unlock(bl(iepsi),i)
      norb=nvirt+nval
      icorr=icano+ncore*ncf
      if(ncore.gt.0) call dynamic_getmem(7*ncore,icorpair)
      call dynamic_getmem(7*nval, ivalpair)
      call dynamic_getmem(7*nvirt,ivirpair)
      call symm_memory_allocator(ncf)
      ichar_count=1
      if (.not.loca) then
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
      docansym=.false.
      if (nsym.gt.0.and..not.loca) docansym=.true.
      call dynamic_getmem(ncf*ncf,ifockAO)
      call FockRead(ncf,ifockAO)
      if (.not.loca) then
      call dynamic_getmem(8*7,ichatacters)
      call dynamic_getmem(8*8,im_table)
      call dynamic_getmem(8*(nval+1),iotable)
      call dynamic_getmem(nval,iorevtable)
      call dynamic_getmem(8*(nvirt+1),ivtable)
      call dynamic_getmem(nvirt,ivrevtable)
      call dynamic_getmem((nval*nval+1)*8,ijtable)
      call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                   bl(ichatacters),ichar_count,bl(im_table))
      call dynamic_mmark()
      call dynamic_getmem(2*8,      ibassym)
      call dynamic_getmem(ncf*ncf,  imat1)
      call dynamic_getmem(ncf*ncf*8,imat2)
      call symm_adapted(ncf,nsym,ichar_count,bl(ichatacters),bl(ifp),
     *                  bl(imat1),bl(imat2),bl(ibassym))
      call dynamic_retmem(1)
      if(ncf.eq.nbf) then
        call Fock_Diag(ncf,ifockAO,ioverlap,icano,iepsi,bl(imat1),
     *                 bl(ibassym),nmo)
      endif
      call dynamic_retmark()
      if(ncore.gt.0)
     *  call core_pairs(nsym,   ncf,  ncore,   bl(ifp),  bl(icorpair),
     *                  bl(icano),bl(ioverlap))
      call orb_pairs(nsym,ncf,nval,nvirt,bl(ifp),bl(ivalpair),
     *               bl(ivirpair),bl(icorr),bl(ioverlap))
      call irr_rep_table(bl(ivalpair),bl(ivirpair),nsym,nval,nvirt,
     *                   bl(ichatacters),ichar_count,bl(im_table))
c
      call tables_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                      bl(ivalpair),bl(iotable),bl(iorevtable))
      call tables_generator(nsym,nvirt,ichar_count,bl(ichatacters),
     *                      bl(ivirpair),bl(ivtable),bl(ivrevtable))
      call ijtab_generator(nsym,nval,ichar_count,bl(ichatacters),
     *                     bl(ivalpair),
     *                     bl(im_table),bl(iorevtable),bl(ijtable))
      endif
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
c                                                 pair has, if given pair
c                                                 is already an image, = 0
c                                                 npairs+1 because of DIIS
      call izeroit(bl(ipairimages),npairs+1)
      call int_array_write(bl(ipairimages),npairs+1,1)
c
      if (iprnt.ge.3) then
        do i=1,nval
        call orb_searcher(i, iprim, ns, signum)
        write(6,*)"Orb ",i,", has an image: ",iprim,",  sym: ",ns,signum
        enddo
      endif
c
c     call orb_pairs_check_ao(ncf,bl(icorr),1d-14)
      call dynamic_lock(bl(icano),i)
      do i=1,nval
        do j=1,i
          ij=ij+1
          call pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                       ns, trans, signum)
        if (ns.eq.0) then
        if (iprnt.ge.3) write(6,'(A,I2,A,I2,I5)') "Pair: ",i,",",j,ij
        ij_unique=ij_unique+1
        call int_array_write(bl(isympairs),ij,ij_unique)
        else
        call int_array_write(bl(isympairs),ij,0)
        if (iprnt.ge.3) 
     *         write(6,'(A,I2,A,I2,I5,A,I2,A,I2,I5,A,I1,L2,L2)')
     *         "Pair: ",i,",",j,ij," image: ",
     *         iprim,",",jprim,ijprim," ns: ",ns,trans,signum
        endif
        ival=int_array(bl(ipairimages),ijprim)
        ival=ival+1
        call int_array_write(bl(ipairimages),ijprim,ival)
        enddo
      enddo
      if (iprnt.ge.3) then
      do i=1,npairs+1
      write(6,'(2I5)') i,int_array(bl(ipairimages),i)
      enddo
      endif
      ij=0 
      if (iprnt.ge.3) then
      do i=1,nval
        do j=1,nval
          ij=ij+1
          call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                            ns, signum)
        if (ns.eq.0) then
          write(6,'(A,I2,A,I2,I5)') "Full pair: ",i,",",j,ij
        else
          write(6,'(A,I2,A,I2,I5,A,I2,A,I2,I5,A,I1,L2)') 
     *           "Full pair: ",i,",",j,ij," image: ",
     *           iprim,",",jprim,ijprim," ns: ",ns,signum
        endif
        enddo
      enddo
      endif
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
      call GenCoulExInt(ncs,    ncf,          ictr,   nval,  nmo,
     *                  nfirst, nlast,        thresh, core,  xmaxdisk,
     *                  ndiskc, bl(irecadrc), lrecc,  lbinc, int_kind,
     *                  nodisk, byt8,         small,  vorb,  ndisktrc,
     *                  bl(ilist),   nstrong, bl(isympairs), ij_unique,
     *                  nbf,iprnt)
      endif
      call secund(sc_stop)
      call elapsec(ec_stop)
      if (iprnt.ge.3) then
      write(6,*) 'Elapsec time for coulomb integral calc. (sec): ',
     * ec_stop-ec_start
      write(6,*) 'CPU time for coulomb integral calc. (sec):     ',
     * sc_stop-sc_start
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
     * ex_stop-ex_start
      write(6,*) 'CPU time for exchange integral calc. (sec):     ',
     * sx_stop-sx_start
      write(6,*) 'Ipass:  ',ipass
      endif
ctj Generate the exchange integrals
      if (l_triples) then
      int_kind='tt'
      call elapsec(ex_start)
      call secund(sx_start)
      call GenCoulExInt(ncs,    ncf,          ictr,   nval,  nmo,
     *                  nfirst, nlast,        thresh, core,  xmaxdisk,
     *                  ndisktt, bl(irecadrtt),lrectt, lbintt, int_kind,
     *                  nodisk, byt8,         small,  vorb,  ndisktrtt,
     *                  bl(ilist),   nstrong, bl(isympairs), ij_unique,
     *                  nbf,iprnt)
c ndisktt,irecadrtt,lrectt,lbintt,ndisktrtt
      call secund(sx_stop)
      call elapsec(ex_stop)
      if (iprnt.ge.3) then
      write(6,*) 'Elapsec time for 3ext integral calc. (sec): ',
     * ex_stop-ex_start
      write(6,*) 'CPU time for 3ext integral calc. (sec):     ',
     * sx_stop-sx_start
      write(6,*) 'Ipass:  ',ipass
      endif
      endif
c 
      call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
     *              vorb,  nmo, ivalpair,     ivirpair,   ifp, 
     *              nsym,  nval,nvirt,        ij_unique,  isympairs,
     *              nbf)
      call matconn('fockAO','q',idimen,idimen,ifockAO)
      call dynamic_getmem(ncf*ncf,ifockMO)
      call dynamic_matdef('corefock','q',idimen,idimen)
      if (mp2) lbinc=0
      icorefock=mataddr('corefock')
      call FockBuilder(bl(irecadrc),bl(irecadrx),nval, npairs, thresh,
     *                 ndiskc,      ndiskx,      lbinc,lbinx,  byt8, 
     *                 nmo,        vorb,     ifockAO, ifockMO,icorefock,
     *                 nbf)
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
      call matzero('invcano')
      call matcopy('cano','invcano')
#if 1
c     call matinv('invcano')
      call dynamic_getmem(ncf,iad2)
      call dynamic_getmem(ncf,iad3)
      call osinv(bl(iinvcano),ncf,det,1.0d-18,bl(iad2),bl(iad3))
      call dynamic_retmem(2)
      if (dabs(det).lt.1d-200) then
        write(6,*) 'DET: ',det
        STOP 'Determinant too small'
      endif
#endif
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
c Determine if DIIS have to be performed or not
      if (.not.log_diis) then
      call CCDiis(i,       j,          'nodiis',     ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
      else
      call prepare_diis_files(ndiskdiisr,ndiskdiisc,af,ncf,nmo,vorb,nbf)
      i=length_diis
      call CCDiis(i,       j,          'diis',       ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx', bl(ilist),
     *            bl(ipairimages),nbf)
      endif
c END DIIS
      xmaxele=0.0d0
      do i=1,nval
        do j=1,i
          ij=ij+1
          equal=.false.
          if (i.eq.j) equal=.true.
          call ExtrOne(i,    j,     bl(irecadrx), npairs,       ndiskx,
     *                 ncf,  lbinx, thresh,       byt8,         'x',
     *                 'mo', nmo,   vorb,         bl(iexchAO))
          call CoefRead('tt',i,j,icoeAO)
         call UpdateCoef(i,       j,      icoeAO,  iexchAO,   bl(iepsi),
     *                   ncf,     nmo,    nfirst,  npairs,    shiftini,
     *                   'nodiis',vorb,  singles, ndiskdiisr,ndiskdiisc,
     *                   sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
          call CoefWrite(i,j,icoeAO)
          call PairEnergyMP2(icoeAO, iexchAO, iresAO, epairr, epair,
     *                       ncf,    equal,   nmo,    vorb,   sepairr,
     *                       sepair,nbf)
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
     * e1_stop-e1_start
      write(6,*) 'CPU time for first iteration:      ',
     * s1_stop-s1_start
      write(6,*)  'Sum of the squares of the residuums: ',square
      write(6,*)  'Maximum residuum element (abs value) ',xmaxele
      write(6,'(A,F21.17)')  'Nonlocal! - "MP2" energy: ',energy
      endif
      if (.not.loca) then 
         xmp2=energy
         scsmp2=scsenergy2
         if (mp2) goto 999
      else
      write(6,'(A,F21.17)')  'For nonlocal - MP2 energy: ',energy
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
     * "Iter:  Energy:        Delta E           "//
     * "Max resid:    Err sq:      Elapsed time"
 222  continue
c
      ij=0
      mp2_iter=mp2_iter+1
      energy=0.0d0
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
          call ExtrOne(i,    j,     bl(irecadrx), npairs,      ndiskx,
     *                 ncf,  lbinx, thresh,       byt8,        'x',
     *                 'mo', nmo,   vorb,         bl(iexchAO))
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
          if (vorb) then
          call matmmul2('fockAO','coefAO','resao','n','n','a')
          call matmmul2('coefAO','fockAO','resao','n','n','a')
          else
          call matmmul2('fockAO','coefAO','work1','n','n','n')
          call matmmul2('work1','overlap','resao','n','n','a')
          call matmmul2('coefAO','fockAO','work1','n','n','n')
          call matmmul2('overlap','work1','resao','n','n','a')
          endif
          call secund(tx1)
          xmultiplication=xmultiplication+(tx1-tx0)
          call matdisc('coefAO')
          call secund(tx0)
          call Read_MP2(i,j,nmo,ncf,nval,ndiskmp2,vorb,.false.,iwork1,
     *                  nbf)
          call secund(tx1)
          reading=reading+(tx1-tx0)
          call secund(tx0)
          if (vorb) then
          call matadd('work1','resao')
          else
          call matsimtr('work1','overlap','work2')
          call matadd('work2','resao')
          endif
          call secund(tx1)
          copy=copy+(tx1-tx0)
          call secund(tx0)
          call CoefRead('tt',i,j,icoeffAO)
          call secund(tx1)
          reading=reading+(tx1-tx0)
          call secund(tx0)
          call PairEnergyMP2(icoeffAO,iexchAO,iresAO,bl(iepairr+ij-1),
     *                                                  bl(iepair+ij-1),
     *                    ncf,equal,nmo,vorb,sepairr,
     *                    sepair,nbf)
         call UpdateCoef(i,       j,      icoeffAO,iresAO,    bl(iepsi),
     *                   ncf,     nmo,    nfirst,  npairs,    shift,
     *                   'nodiis',vorb,  singles, ndiskdiisr,ndiskdiisc,
     *                   sqij,    xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
          square=square+sqij
          call secund(tx1)
          update_energy=update_energy+(tx1-tx0)
          call secund(tx0)
          call CoefWrite(i,j,icoeffAO)
          call secund(tx1)
          reading=reading+(tx1-tx0)
          call secund(tx0)
          call PairEnergyMP2(icoeffAO,iexchAO,iresAO,bl(ixepairr+ij-1),
     *                                                 bl(ixepair+ij-1),
     *                    ncf,equal,nmo,vorb,sxepairr,
     *                    sxepair,nbf)
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
        enddo
      enddo
      if (iprnt.ge.3) then
      write(6,'(A30,F10.2)') 'Time for Generate_MP2_G: ',txx1-txx0
      write(6,'(A30,F10.2)') 'Reading in MP2: ',reading
      write(6,'(A30,F10.2)') 'Multiplication in MP2: ',xmultiplication
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
c
      write(6,'(A,I3)')      'MP2 iteration no:              ', mp2_iter
      endif
      if (iprnt.ge.2) then
      if (.not.do_mp4) then
      write(6,'(A,F21.17)')  'Iterations - MP2 linear prev:  ',energy1
      write(6,'(A,F21.17)')  'Iterations - MP2 energy+resid: ',energyr
      write(6,'(A,F21.17)')  'Iterations - MP2 linear curr:  ',energy2
      else
      write(6,'(A,F21.17)')  'Iterations - MP4 energy:       ',energy2
      endif
      write(6,*)  ' '
      write(6,'(A,F6.2)')  'Total time for iteration CPU:  ',t1-t0
      write(6,'(A,F6.2)')  'Total time for iteration ELA:  ',et1-et0
      write(6,'(A,F20.16)')  'Sum of squares of the residuum:',square
      call flush(6)
      endif
      energyrold=energyr
      if (dabs(energy2-energy1).gt.Ethresh.or.sqij.gt.Wavethresh) then
        goto 222
      else
        if (.not.mp2last) then
          mp2last=.true.
          goto 222
        endif
      endif
      mp2last=.false.
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
      if (do_mp4) goto 999
      if (mp2) goto 999
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
      call sort_amplit(nval,idimen,ndisk_a,af,0,
     *                 ichar_count,bl(ivtable),bl(im_table),bl(ijtable),
     *                 bl(ivrevtable),c1_ratio)
c     call sort_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,ndisktt,
c    *                     ncf,lbintt,thresh,byt8,nmo,vorb,af,0)
      call simple_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                 ndisktt,ncf,lbintt,thresh,byt8,
     *                 nmo,vorb,af,0,ichar_count,bl(ivtable),
     *                 bl(im_table),bl(ivrevtable),bl(iotable))
      call make_3ext_pairs(idimen,nval,af,ndisk_ie,ndisk_ie1)
      call sort_3int(nval,idimen,ndisk_ii,af,0)
      amp_ratio=1d0
      ext3_ratio=1d0
      nslv=1
      call dynamic_show_free(mem)
      mem=mem-nval*nval*nval-3*nval*nval
      call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                     nslv,npass,isize)
      call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *               ncf,lbinx,thresh,byt8,nmo,vorb,af,0,isize,npass)
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
      call new_triples(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,ndisk_ix,
     *                 res,bl(iepsi+ncore),bl(itimes),20,ccsd,
     *                 qcisd,af,bl(iotable),bl(ivtable),bl(ivrevtable),
     *                 bl(im_table),bl(ijtable),ichar_count,npass,isize,
     *                 energys,energyd,iprnt)
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
 666  continue
      call elapsec(eistart)
      call secund(sistart)
      do i=1,itim_no       
      ccsdcpu(i)=0.0d0
      ccsdela(i)=0.0d0
      enddo
c
      call CCDiis(i,       j,          'itera',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
c
      if (iprnt.ge.3) then
        write(icond,*) '                 The CCSD iteration: ',iteration
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
        call NEW_EEO_INT(ncs,    ncf,      bl(ictr),   nval,  nmo,
     *                   nfirst, nlast,    thresh,
     *                   vorb,   ndisktre, npairs, af,iprnt,nbf)
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
      call GenCoulExInt(ncs,     ncf,           ictr,  nval,  nmo,
     *                  nfirst,  nlast,         thresh,core,  xmaxdisk,
     *                  ndisktx, bl(irecadrtx), lrectx,lbintx,int_kind,
     *                  nodisk,  byt8,          small, vorb,  ndisktrtx,
     *                  bl(ilist),   nstrong,   bl(isympairs),ij_unique,
     *                  nbf,iprnt)
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(2)=t1-t0        ! TEIO tx
      ccsdela(2)=et1-et0      ! TEIO tx
c
      call secund(t0)
      call elapsec(et0)
      int_kind='tc' ! TEIO operator
c     call izeroit(bl(irecadrtc),nval*nval)     !zero out the pair records
      call GenCoulExInt(ncs,     ncf,           ictr,  nval,  nmo,
     *                  nfirst,  nlast,         thresh,core,  xmaxdisk,
     *                  ndisktc, bl(irecadrtc), lrectc,lbintc,int_kind,
     *                  nodisk,  byt8,          small, vorb,  ndisktrtc,
     *                  bl(ilist),   nstrong,   bl(isympairs),ij_unique,
     *                  nbf,iprnt)
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(3)=t1-t0        ! TEIO tc
      ccsdela(3)=et1-et0      ! TEIO tc
c
c
      endif
c
      call secund(t0)
      call elapsec(et0)
      do i=1,5
      elaps(i)=0.0d0
      enddo
      call CreateAlphaDisk(nval,ndiskalpha,af)
      if (cc.or.ccsd.or.do_mp4) then
      call CCalphaonce(ncf,nval,bl(irecadrx),npairs,ndiskx,lbinx,
     *                 thresh,norecalc,byt8,ccsd,nmo,vorb,elaps,
     *                 ndiskalpha,af,nbf)
      if (iprnt.ge.3) then
      write(6,'(A)') 'CCalph:    ampl:       Kij:     mult:  putres:'
      write(6,'(7X,4F10.2)') elaps(1),elaps(2),elaps(3),elaps(4)
      endif
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(4)=t1-t0        ! CCalphaonce+supplement
      ccsdela(4)=et1-et0      ! CCalphaonce+supplement
c
      call secund(t0)
      call elapsec(et0)
      call SupplementAlpha_BetaBuild(nval,ndiskalpha,bl(ibeta),ncf,
     *                           nfirst,nmo,vorb,cc,ccsd,bl(ifockMO),
     *                           kijklndisk,af,do_mp4,nbf)
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
c
      ccsdcpu(5)=t1-t0        ! Builder41a_G
      ccsdela(5)=et1-et0      ! Builder41a_G
c 
      call secund(t0)
      call elapsec(et0)
      if (cc.or.ccsd.or.do_mp4) then
      call CCAOnce(bl(irecadrx),npairs,nval,ndiskx,ncf,lbinx,thresh,
     *             norecalc,byt8,ioverlap,nmo,vorb,nbf)
      endif
      call secund(t1)
      call elapsec(et1)
      ccsdcpu(6)=t1-t0        ! CCAOnce
      ccsdela(6)=et1-et0      ! CCAOnce
c 
      ij=0
      xlinear_energy=0.0d0
c
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
     *                          lbine,thresh,nfirst,nlast,byt8,
     *                          nmo,vorb,iiS,nbf)
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
     *                    ncf,   lbinx,  thresh,      byt8,  nmo,
     *                    vorb,  iiR)
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
      if (iprnt.ge.3)
     * write(6,*)'istart,istop,jstart,jstop: ',istart,istop,jstart,jstop
        call CCYZ(istart,istop,jstart,jstop,ncf,nval,bl(irecadrx),
     *            bl(irecadrc),npairs,ndiskx,
     *        ndiskc,lbinx,lbinc,bl(irecadrtx),ndisktx,lbintx,
     *        bl(irecadrtc),ndisktc,lbintc,thresh,
     *        byt8,ioverlap,cc,ccsd,bl(ilist),nmo,vorb,ndiskYZ,.false.,
     *        elaps,do_mp4,
     *              ichar_count,bl(iorevtable),bl(ivtable),bl(im_table),
     *        nbf)
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
     * write(6,*)'istart,istop,jstart,jstop: ',istart,istop,jstart,jstop
          call Qparts(istart,istop,jstart,jstop,ncf,nval,npairs,byt8,cc,
     *              ccsd,bl(ilist),nmo,vorb,ndiskYZ,af,ndiskQ,
     *              ichar_count,bl(iorevtable),bl(ivtable),bl(im_table),
     *              nbf)
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
          call ExtrOne(i,    j,     bl(irecadrx), npairs,     ndiskx,
     *                 ncf,  lbinx, thresh,       byt8,       'x',
     *                 'mo', nmo,   vorb,         bl(iexchAO))
c  Read external echange operators
          call ExtrOne(i,    j,     bl(irecadre),   npairs,    ndiske,
     *                 ncf,  lbine, thresh,         byt8,      'e',
     *                 'mo', nmo,   vorb,           bl(iwork1))
c  Add it do residuum matrix
          call matcopy('exchao','resao')
          call matadd('work1','resao')
          else
          call matzero('resao')
          endif
c
c  Qgen, iwork1=Qij, iwork2=Qji(T), probably OK, order of en. OK
          call Q_read_build(i,j,iX,ndiskQ,af,vorb,ncf,nmo,nval,
     *                      iwork1,iwork2,nbf)
c
c   S*Qij and add to residuum
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
c   Read Gij and Gji(T) with part of CCSD
c
          call ReadCalcG_41a(i,j,norecalc,ncf,nmo,vorb,ndiskG41a,
     *                       iwork1,.false.,nbf)
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
          if (vorb) then
          call matadd('work1','resao')
          else
          call matsimtr('work1','overlap','work2')
          call matadd('work2','resao')
          endif
c
c Singles part of CISD CCSD start
          if (singles.or.ccsd) then
          call dynamic_matdef('CISD','q',idimen,idimen)
          iresult=mataddr('CISD')
          call singles_CISD(i,j,ncf,nval,ioverlap,bl(irecadre),npairs,
     *               nfirst,nlast,ndiske,lbine,thresh,byt8,ccsd,iresult,
     *               nmo,vorb,nbf)
          call matadd('CISD','resao')
          call dynamic_matrem('CISD')
          endif
c Singles part of CISD CCSD stop
          call CoefRead('tc',i,j,icoeffAO)
          if (mp3.or.do_mp4) then
          call PairEnergy(icoeffAO,iresAO,iresAO,
     *             bl(iepairr+ij-1),bl(iepair+ij-1),ncf,equal,nmo,vorb,
     *             nbf)
            if (mp3.and.mp4) call writer(ndisk_mp3,idimen,ij,bl(iresAO))
          else
          call PairEnergy(icoeffAO,iexchAO,iresAO,
     *             bl(iepairr+ij-1),bl(iepair+ij-1),ncf,equal,nmo,vorb,
     *             nbf)
          endif
c
          call CoefRead('tt',i,j,icoeffAO)
          if (.not.(mp3.or.do_mp4)) then
          call CINorm(icoeffAO,xnormij,ncf,equal,nmo,vorb,nbf)
          totnorm=totnorm+xnormij
         call UpdateCoef(i,      j,    icoeffAO, iresAO,     bl(iepsi),
     *                   ncf,    nmo,  nfirst,   npairs,     shift,
     *                   'diis', vorb, singles,  ndiskdiisr, ndiskdiisc,
     *                   sqij,   xmaxele,bl(ilist),bl(ipairimages),
     *                   bl(iorevtable),bl(ivrevtable),bl(im_table),
     *                   docansym,nbf)
          endif
c
          call CoefWrite(i,j,icoeffAO)
c
          call CoefRead('ee',i,j,icoeffAO)
          call PairEnergy(icoeffAO,iexchAO,iresAO,
     *            bl(ixepairr+ij-1),bl(ixepair+ij-1),ncf,equal,nmo,vorb,
     *             nbf)
c
 777      continue  ! Jump to this point if the pair has to be omitted
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
c
        enddo
      enddo
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
c
      call pointS(iiS)
      call dynamic_matdef('s_resid','r',idimen,nval)
      is_resid=mataddr('s_resid')
      call pointersingles(isingles)
      call pointernewsingles(inewsingles)
      call matconn('full_sing','r',idimen,nval,isingles)
      call tfer(bl(iiS),bl(is_resid),idimen*nval)
      call Tr_generator1(idimen,nval,ncf,vorb,ioverlap,is_resid)
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
c
      mp4_if: if (do_mp4) then
        call matconn('fockMO','q',ncf,ncf,ifockMO)
        total_singl_en=0.0d0
        call matmmult('spwork1','fockMO','spwork2')
c       call compact_print(bl(ifockMO),ncf)
        if (vorb) then
          call matprodtr('spwork2','spwork1',part_singl)
        else
          write(6,*) 'Not working yet'
          STOP 'Not working yet'
          call matmmul2('spwork2','spwork1','spwork3','n','t','n')
          call matprodtr('spwork3','overlap',part_singl)
        endif
        total_singl_en=total_singl_en+part_singl
        call matmmult('fockMO','spwork1','spwork2')
        if (vorb) then
          call matprodtr('spwork2','spwork1',part_singl)
        else
          write(6,*) 'Not working yet'
          call matmmul2('spwork2','spwork1','spwork3','t','n','n')
          call matprodtr('spwork3','overlap',part_singl)
        endif
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
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,*) '          The results for iteration: ',iteration
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,*)  'Sum of the squares of the residuums:',square
      write(6,*)  'Maximum residuum element (abs value)',xmaxele
      endif
      if (mp3.or.do_mp4) then 
        if (mp3.and.iprnt.ge.2) 
     *     write (*,'(A,F25.16)') 'MP3 Energy is: ', energy
        energy_mp3=energy
        if (do_mp4) then 
          xmp4_q=energy
        else
          xmp3=energy
        endif
        goto 70
      endif
      if (iprnt.ge.2) then
      if (cep0) write(6,'(A)') '* * * * * * * * * * * * * * * * CEPA-0'
     *//' Energies: * * * * * * * * * * * * * * *'
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
      write(6,'(A)')'Remember! In fact this energy is not variational'//
     *             ' when Gauss-Seidel'
      write(6,'(A)') '                             algorithm was used'
      endif
      else if (linear) then
      xiterenergy=xlinear_energy   ! for CID linear
      if (iprnt.ge.3) then
      write(6,'(A)')'Linear energy will be used for the next iteration.'
      endif
      else
      xiterenergy=energyr ! for CID quadratic (default)
      if (iprnt.ge.3) then
      write(6,'(A)') 'Quadratic (non-variational!) energy will be used '
     *           //'for the next iteration.'
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
     *              ccsdela(12)/6d1
      if (.false..and.omit.and.DeltaE.le.1.0d-7) then  ! disabled
        call filllist(bl(ilist),npairs)
        omit=.false.
      call CCDiis(i,       j,          'reset',      ncf,  nmo,
     *            npairs,  bl(iresAO), bl(icoeffAO), nval, vorb,
     *            singles, ndiskdiisr, ndiskdiisc,   'xxx',bl(ilist),
     *            bl(ipairimages),nbf)
        goto 666
      endif
      if (xmaxele.gt.Wavethresh .or. DeltaE.gt.Ethresh) then
        if (iteration.ge.maxiter) then
          write(6,*) '            * * * * Maximum number of'//
     *            ' iterations reached * * * *'
          write(6,*) '            * * * *             NO '//
     *            'CONVERGENCE           * * * *'
        else
          goto 666
        endif
      endif
 70   continue
      call flush(6)
      if (do_mp4) goto 222
      if (mp4.and..not.do_mp4) then 
        mp3=.false.
        do_mp4=.true.
        goto 666
      endif
      call flush(6)
      if ((cc.and.singles).and.l_triples) then
c     call triples(nval,idimen,npairs,ncf,nmo,vorb,bl(iepsi+ncore),res,
c    *             ccsd,bl(irecadrx),ndiskx,lbinx,thresh,byt8,
c    *             singles.and.cc.and..not.ccsd,af,1)
c     write(6,'(A,F20.15)') 'CCSD triples correction: ',res
      call elapsec(xt0)
      call sort_amplit(nval,idimen,ndisk_a,af,0,
     *                 ichar_count,bl(ivtable),bl(im_table),bl(ijtable),
     *                 bl(ivrevtable),c1_ratio)
c     call sort_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,ndisktt,
c    *                     ncf,lbintt,thresh,byt8,nmo,vorb,af,0)
      call simple_3ext(nval,idimen,ndisk_ie,bl(irecadrtt),npairs,
     *                 ndisktt,ncf,lbintt,thresh,byt8,
     *                 nmo,vorb,af,0,ichar_count,bl(ivtable),
     *                 bl(im_table),bl(ivrevtable),bl(iotable))
      call make_3ext_pairs(idimen,nval,af,ndisk_ie,ndisk_ie1)
      call sort_3int(nval,idimen,ndisk_ii,af,0)
      amp_ratio=1d0
      ext3_ratio=1d0
      nslv=1
      call dynamic_show_free(mem)
      mem=mem-nval*nval*nval-3*nval*nval
      call calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                     nslv,npass,isize)
      call sort_Kext(nval,idimen,ndisk_ix,bl(irecadrx),npairs,ndiskx,
     *               ncf,lbinx,thresh,byt8,nmo,vorb,af,0,isize,npass)
      call elapsec(xt1)
      if (iprnt.ge.3) write(6,*) 'Sort time: ',xt1-xt0
      call flush(6)
      call elapsec(xt0)
      call dynamic_getmem(20,itimes)
      call new_triples(nval,idimen,ndisk_a,ndisk_ie1,ndisk_ii,ndisk_ix,
     *                 res,bl(iepsi+ncore),bl(itimes),20,ccsd,
     *                 qcisd,af,bl(iotable),bl(ivtable),bl(ivrevtable),
     *                 bl(im_table),bl(ijtable),ichar_count,npass,isize,
     *                 energys,energyd,iprnt)
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
c---------------------------------------------------------------------
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
     *     method(1:lmet)//' triples correction:'//spc,res
        else
          res=0d0
        endif
        if (ifound(2).eq.0) then
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
c 
      call memory_status('end of CCSD')

#if 0
      do i=1,nval
        do j=1,i
        call CoefRead('tt',i,j,icoeffAO)
          write(6,*) i,j
          do k=1,idimen
            do l=1,idimen-1
              write(6,'(F10.5)',ADVANCE="NO") 
     *                 bl(icoeffAO+(k-1)+(l-1)*idimen)
            enddo
            l=l+1
            write(6,'(F10.5)',ADVANCE="YES") 
     *               bl(icoeffAO+(k-1)+(l-1)*idimen)
            write(6,*) '***********************************************'
          enddo
        enddo
      enddo
#endif
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
c---------------------------------------------------------------------
      end
c=======================================================================
      subroutine matrix_int_zeroit(lzero,ncf)
      implicit real*8 (a-h,o-z)
      dimension lzero(*)
      do i=1,ncf
      lzero(i)=0
      enddo
      end
C==============================================================================
      subroutine ext3_WriteLeftOvers(ndisk, npairs, lbin, icountr, irec,
     *                          indexx, ibin4,   ibin1, bin8,      byt8,
     *                          af, islvid,nval,ncf,nmo,igranules,
     *                          igranulesize,igran2pair,ipair2gran,
     *                          irec_no,irec_last,nbf)
c     Writes leftovers of bins on the disk
c     INPUT:
c     ndisk    - number of disk unit
c     npairs   - number of ij pairs
c     lbin     - length of the bin 
c     icountr  - the array which "remember" how much of the bin is filled
c     indexx   - mu, nu indices for given integral stored in bin
c     ibin4    - bin with integer*4
c     ibin1    - bin with integer*1
c     bin8     - the double precision bin
c     byt8     - logical variable controlling whther we use ibin4 & ibin1 pair 
c                or bin8
c     af       - logical variable, if true this subroutine assumes that 
c                two things happen: job is parallel and we use Array Files
c     islvid   - group ID, used as one of the indices. Master has 0 and slaves 
c                from 1 to # of slaves. Irrelevant in single mode.
c     OUTPUT:
c     irec     - the array with the number of last record for given ij pair
c                (single) or number of records written for given ij (parallel)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      implicit real*8 (a-h,o-z)
      integer*4 ibin4(lbin,nval*(nbf-nmo))
      integer irec(nval*(nbf-nmo)), icountr(nval*(nbf-nmo))
      real*8    bin8(lbin,nval*(nbf-nmo))
      integer*2 indexx(2,lbin,nval*(nbf-nmo))
      integer*1 ibin1(lbin,nval*(nbf-nmo))
      integer islvid
      logical byt8,af
      integer*4 igran2pair(2,*),ipair2gran(*)
c     First put zeros in the non-filled part of bins and index, then
c     write it
      if (af) then
      do igr=1, igranules
      istart=igran2pair(1,igr)
      istop =igran2pair(2,igr)
      isize=istop-istart+1
      do ij=istart,istop
        if (byt8) then
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            bin8(icount,ij)=0.0d0
          enddo
        else
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            ibin4(icount,ij)=0
            ibin1(icount,ij)=0
          enddo
        endif
      enddo
      call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),ibin1(1,istart),
     1                  bin8(1,istart),indexx(1,1,istart),irec(igr),igr,
     2                  islvid,byt8,isize,irec_no,irec_last)
      enddo
      else
      do ij=1, nval*(nbf-nmo)
        if (byt8) then
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            bin8(icount,ij)=0.0d0
          enddo
        else
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            ibin4(icount,ij)=0
            ibin1(icount,ij)=0
          enddo
        endif
        call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij), bin8(1,ij),
     +                  indexx(1,1,ij),irec(ij), ij,nval*(nbf-nmo),byt8)
      enddo
      endif
      end
C==============================================================================
C==============================================================================
      subroutine WriteLeftOvers(ndisk,  npairs,  lbin,  icountr,   irec,
     *                          indexx, ibin4,   ibin1, bin8,      byt8,
     *                          af,     islvid,  irec_no,irec_last)
c     Writes leftovers of bins on the disk
c     INPUT:
c     ndisk    - number of disk unit
c     npairs   - number of ij pairs
c     lbin     - length of the bin 
c     icountr  - the array which "remember" how much of the bin is filled
c     indexx   - mu, nu indices for given integral stored in bin
c     ibin4    - bin with integer*4
c     ibin1    - bin with integer*1
c     bin8     - the double precision bin
c     byt8     - logical variable controlling whther we use ibin4 & ibin1 pair 
c                or bin8
c     af       - logical variable, if true this subroutine assumes that 
c                two things happen: job is parallel and we use Array Files
c     islvid   - group ID, used as one of the indices. Master has 0 and slaves 
c                from 1 to # of slaves. Irrelevant in single mode.
c     OUTPUT:
c     irec     - the array with the number of last record for given ij pair
c                (single) or number of records written for given ij (parallel)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      implicit real*8 (a-h,o-z)
      integer*4 ibin4(lbin,npairs)
      integer irec(npairs), icountr(npairs)
      real*8    bin8(lbin,npairs)
      integer*2 indexx(2,lbin,npairs)
      integer*1 ibin1(lbin,npairs)
      integer islvid
      logical byt8,af
c     First put zeros in the non-filled part of bins and index, then
c     write it
      do ij=1, npairs
        if (byt8) then
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            bin8(icount,ij)=0.0d0
          enddo
        else
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            ibin4(icount,ij)=0
            ibin1(icount,ij)=0
          enddo
        endif
        if (af) then
          call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indexx(1,1,ij),irec(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
        else
        call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij), bin8(1,ij),
     +                  indexx(1,1,ij),irec(ij), ij,npairs,byt8)
        endif
      enddo
      end
C==============================================================================
      subroutine WriteLeftOversGr(ndisk, npairs, lbin,  icountr, irec,
     *                           indexx,ibin4,  ibin1, bin8,    byt8,
     *                           af,    islvid, igranules,igranulesize,
     *                          igran2pair,ipair2gran,irec_no,irec_last)
c     Writes leftovers of bins on the disk, in groups.
c     INPUT:
c     ndisk    - number of disk unit
c     npairs   - number of ij pairs
c     lbin     - length of the bin 
c     icountr  - the array which "remember" how much of the bin is filled
c     indexx   - mu, nu indices for given integral stored in bin
c     ibin4    - bin with integer*4
c     ibin1    - bin with integer*1
c     bin8     - the double precision bin
c     byt8     - logical variable controlling whther we use ibin4 & ibin1 pair 
c                or bin8
c     af       - logical variable, if true this subroutine assumes that 
c                two things happen: job is parallel and we use Array Files
c     islvid   - group ID, used as one of the indices. Master has 0 and slaves 
c                from 1 to # of slaves. Irrelevant in single mode.
c     igranules- number of groups of bins. Might be from 1 to npairs.
c     igranulesize- number of bins in one group of bins (granule)
c     igran2pair - array mapping granule number to two numbers: starting ij 
c                index and closing ij index
c     ipair2gran - as above but in the opposite way, it maps ij pair to 
c                granule number
c     OUTPUT:
c     irec     - the array with the number of last record for given ij pair
c                (single) or number of records written for given ij (parallel)
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
      implicit real*8 (a-h,o-z)
      integer*4 ibin4(lbin,npairs)
      integer irec(npairs),  icountr(npairs)
      real*8    bin8(lbin,npairs)
      integer*2 indexx(2,lbin,npairs)
      integer*1 ibin1(lbin,npairs)
      integer*4 igran2pair(2,igranules)
      integer*2 ipair2gran(npairs)
      integer islvid
      logical byt8,af
c     First put zeros in the non-filled part of bins and index, then
c     write it
      do igr=1, igranules
      istart=igran2pair(1,igr)
      istop =igran2pair(2,igr)
      isize=istop-istart+1
      do ij=istart,istop
        if (byt8) then
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            bin8(icount,ij)=0.0d0
          enddo
        else
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            ibin4(icount,ij)=0
            ibin1(icount,ij)=0
          enddo
        endif
      enddo
        if (af) then
          call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),ibin1(1,istart),
     1                 bin8(1,istart),indexx(1,1,istart),irec(igr),igr,
     2                 islvid,byt8,isize,irec_no,irec_last)
        else
        call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij), bin8(1,ij),
     +                  indexx(1,1,ij),irec(ij), ij,npairs,byt8)
        endif
      enddo
      end
C==============================================================================
      subroutine WriteLeftOversExt(ndisk,npairs,nval,lbin,icountr,
     *                             irec,indexx,ibin4,ibin1,bin8,
     *                             byt8,af,islvid,irec_no,irec_last)
      implicit real*8 (a-h,o-z)
      integer*4 ibin4(lbin,nval*nval)
      integer irec(nval*nval),   icountr(nval*nval)
      real*8    bin8(lbin,nval*nval)
      integer*2 indexx(2,lbin,nval*nval)
      integer*1 ibin1(lbin,nval*nval)
      integer islvid
      logical byt8,af
c     Writes leftovers on the disk
c     See the WriteLeftOvers description. This subroutine writes bin for 
c     TEIO operators (three-external) needed for CCSD. It is different because
c     bin is longer, not npairs dimension, but nval*nval
c
c     First put zeros in the non-filled part of bins and index, then
c     write it
      do ij=1, nval*nval
        if (byt8) then
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            bin8(icount,ij)=0.0d0
          enddo
        else
          do icount=icountr(ij)+1, lbin
            indexx(1,icount,ij)=0
            indexx(2,icount,ij)=0
            ibin4(icount,ij)=0
            ibin1(icount,ij)=0
          enddo
        endif
        if (af) then
          call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indexx(1,1,ij),irec(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
        else
        call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij), bin8(1,ij),
     +                  indexx(1,1,ij),irec(ij), ij,nval*nval,byt8)
        endif
      enddo
      end
C==============================================================================
      subroutine GenCoulExInt(ncs,   ncf,    ictr,   nval,  nmo,
     *                        nfirst,nlast,  thresh, core,  xmaxdisk,
     *                        ndisk, irecadr,lrec,   lbin,  int_kind,
     *                        nodisk,byt8,   small, vorb,    ndisktr,
     *                        list,  nstrong,isympairs,ij_unique,nbf,
     *                        iprnt)
      use memory
      use kinds
      use ccounters
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
c irecadr  - irecadr contains the number of last record of the ij bin,
c            where ij is pair number. Previous record numbers are
c            kept with bins and are read from disk.
c int_kind - it is character string indicating the kind of calculated
c            integrals. May be: 'c' coulomb, 'x' exchange, 'e' EEO, 
c            'tc' - three external coulomb, 'tx' - three external exchange
c            'tt' - three external exchange for triples
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
      implicit real*8 (a-h,o-z)
      logical LastSymPair,emp2only,nodisk,keepnodisk,byt8,small,vorb,af
      logical memory_success,reset
      dimension xintxx(9)
      common /job/jobname,lenJ
      parameter (max_merg=(MAX_MERGE))  ! fragile!!
      common /indexsaver/ik,ikindex(2,max_merg) ! for EEOGenerator
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*3 ch3
      character scftype*11
      character*256 jobname,scrfile,filename
      logical nofr,exst,restrt,dualbasis,smal,fileopen,split
      character*(*) int_kind
      dimension xnmo(2)
ctj header end
      dimension irecadr(*)
      character*256 filname1
      dimension Ipass(2,28), Jpass(2,28) !28 is 28 comp.of cartisian I-function
      parameter (i_orb_size=28)
      parameter (ixint_time=10)
      real*8 xint_time(ixint_time)
      integer list(*),nstrong,isympairs(*),ij_unique
      logical image
c
ctj Int_kind=c coulomb, Int_kind=x exchange, Int_kind=e EEO
ctj -- BEGIN get some common variables
      keepnodisk=nodisk
      dualbasis=.FALSE.
      split=.true.
      call getival('iout',iout)
      smal=.TRUE.
      itott=0
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
ctj -- END   get some common variables
  46  format(/72('='))
c     write(iout,46)
c     write(6,*)  'Int_kind=',Int_kind,'.'
      if (Int_kind.eq.'c') then
      write(iout,*) '                            The Jij Module  '
      else if (Int_kind.eq.'x') then
      write(iout,*) '                            The Kij Module  '
      else if (Int_kind.eq.'e') then
      write(iout,*) '                            The EEOij Module  '
      else if (Int_kind.eq.'tx') then
      if (iprnt.ge.2)
     *write(iout,*) '                         The TEIO exch ij Module  '
      else if (Int_kind.eq.'tc') then
      if (iprnt.ge.2)
     *write(iout,*) '                         The TEIO coul ij Module  '
      else if (Int_kind.eq.'tt') then
      write(iout,*) '                   Three external int. for triples'
      else
        call nerror(1,'GenCoulExInt',
     1  'Int_kind does not indicate any kind of integrals',0,0)
      endif
c     write(iout,*)' '
c-----------------------------------------------------------
      ncs_sm=ncs
      ncf_sm=ncf
c-----------------------------------------------------------
c This is for use in EEO Generator and for bins construction when we
c decide how much memory we need to reserve
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ncf*(ncf-1)/2
c-----------------------------------------------------------
      if (iprnt.ge.3) then
      if (Int_kind.eq.'c') then
      write(iout,*) ' Jij integral thresh    = ',thresh
      else if (Int_kind.eq.'x') then
      write(iout,*) ' Kij integral thresh    = ',thresh
      else if (Int_kind.eq.'e') then
      write(iout,*) ' EEOij integral thresh    = ',thresh
      else if (Int_kind.eq.'tx') then
      write(iout,*) ' TEIO exch ij integral thresh    = ',thresh
      else if (Int_kind.eq.'tc') then
      write(iout,*) ' TEIO coul ij integral thresh    = ',thresh
      else if (Int_kind.eq.'tt') then
      write(iout,*) ' TEI for triples ij integral thresh    = ',thresh
      endif
      write(iout,*) ' core orbitals with eps < ', core
      write(iout,*) ' max disk storage (MW)  = ',xmaxdisk*1.d-6
      write(iout,*) ' '
      call flush(iout)
      endif
c-----------------------------------------------------------
c  put down a memory marker
c
      call mmark
      call dynamic_mmark
      call matmark
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
c Perform schwartz integrals calculation:
      call dynamic_getmem(ncs*ncs,ischwarz)
      call schwarz_calc(bl(ictr),ncs,bl(ischwarz))
c
c get memory for screening density
      np4=4
      call dynamic_matdef('dsmx','q',ncs,ncs) ! filled in <DmxMakeC>/<DmxMakeL>
      idics=mataddr('dsmx')
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
      call mmark
      call dynamic_mmark
c---------------------------------------------------------------------
c
       call sudat(np4,'loca_rhf',ni)
c WARNING! Screening is turned off at present!
      call  fill1(ncs,bl(idics))
      call dynamic_lock(bl(idics),iii)
c
c---------------------------------------------------------------------
c
      icano=mataddr('cano')
c.................................................
c check if there will be a split in integral calculations
c
      call check_sizes(bl(ictr),ncs,ncf,nval)
c.................................................
      if(ncore.gt.0) call matsub('occa','cano',1,nmo)
      call  matsub('occu','cano',nfirst,nlast)
      call  matsub('virt','cano',nmo+1,nbf)
      ioccu=mataddr('occu')
c.................................................
      if(iprnt.gt.3) then
        call matprint ('occu',6)
        call matprint ('virt',6)
      end if
c .................................................
      nvirt=nbf-nmo                        !number of virtual orbitals
c   reserve space for one AO integral matrix
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
      call mmark
      call dynamic_mmark
c
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
      if(iprnt.gt.3) write(iout,*) 'Startup time=',tt-tt0
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c  ijtot is the number of ics,jcs pairs really calculated & transformed
c   if an (ics,jcs) pair is skipped because there are no integrals in it
c   then it is NOT incremented
      ijtot=0
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
      call matrix_int_zeroit(bl(lzero),ncf)
C
CTJ
      npairs=(nval**2+nval)/2
#if 1
c The copy from From Alan Ford file:
CTJ
c
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      call dynamic_getmem(npairs/intsize,icountr) !memory for bin counters
      call izeroit(bl(icountr),npairs)         !zero out the bin counters
      else if (int_kind.eq.'tt') then
      call dynamic_getmem(nval*nvirt/intsize,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nvirt)         !zero out the bin counters
      else
      call dynamic_getmem(nval*nval/intsize,icountr) !memory for bin counters
      call izeroit(bl(icountr),nval*nval)         !zero out the bin counters
      endif
!     Determine how much memory we can allocate to the bins.  They should
!     be as big as possible but leave space for other allocations.
      memreal=igetival('memr')
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      memvirt=memtot
      memoccu=lastaddr-ioffset
      memavail=memvirt-memoccu
      call dynamic_show_free(memavail)
      if (iprnt.ge.3) then
      write(6,*) 'memavail: ',memavail,'  memreal: ',memreal
      endif
c  memory for the bins is the minimum of the available virtual and real memory
c  size, minus 1 MW for miscellaneous - the space needed for matrix multiplies etc.
      maxshell=49 !max contracted shell  squared
cAF memory for first half transform
      mem1trans=npairs+nval**2+nval**2/8+3*ncf+ncf**2+
     1          ncf*ncf*maxshell+nval**2
cAF memory for second half transform
ctj Not needed here.
c     mem2trans=3*npairs+2*nvirt**2+2*ncf*nvirt
      mem2trans=0
      memleft=memavail-mem1trans-mem2trans-1 000 000
cc      memleft=min0(memreal,memavail)-mem1trans-mem2trans-1 000 000
c In the case with EEO it is not desirable to allocate as much memory
c for bins as possible. There will be only one write of result of size
c npairs*ncf*ncf. But, we will have to read many times the amplitudes,
c having the same size, in the icf kcf loop. For small shells merging it
c means that we may read the over 1000. The time spent in bins writing
c in that case will be negligible.
      if (Int_kind.eq.'e') then
        write(6,*) 'Reducing available memory for bins in EEO'//
     *             ' operator procedure'
c       memleft=2000 000
        memleft =0.3d0*memleft
      endif
c  lbin is in units of 9-bytes (2 indices           - integer*2
c                               compressed integral - integer*4
c                               precision overflow  - integer*1)
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      if (byt8) then
      lbin=8*memleft/(12*npairs) !Size of bins in 8 byte words        
      else
      lbin=8*memleft/(9*npairs) !Size of bins in 8 byte words        
      endif
      else if (int_kind.eq.'tt') then
      if (byt8) then
      lbin=8*memleft/(12*nval*nvirt) !Size of bins in 8 byte words        
      else
      lbin=8*memleft/(9*nval*nvirt) !Size of bins in 8 byte words        
      endif
      else
      if (byt8) then
      lbin=8*memleft/(12*nval*nval) !Size of bins in 8 byte words        
      else
      lbin=8*memleft/(9*nval*nval) !Size of bins in 8 byte words        
      endif
      endif
      if(lbin.lt.50) then
        call nerror(4,'GenCoulExInt',
     1  'memory available for bins leads to bin size < 50, mem,lbin',
     2  memleft,lbin)
      end if
c  lbin is the length of a bin in 9-byte (or 12 byte) words
      if(lbin.gt.ncf*(ncf+1)) then
        lbin=ncf*(ncf+1)
      end if
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
      indexmem=npairs*lbin/i4size+1
      if (byt8) then
      call dynamic_getmem(npairs*lbin,ibin8) !memory for bins
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin8),npairs*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=npairs*lbin/i4size+1
      int1mem=npairs*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      if (iprnt.ge.3) then
      write(6,*)'mem,npairs,lbin,nbins',memleft,npairs,lbin,ncf*ncf/lbin
      call flush(6)
      endif
      else if (int_kind.eq.'tt') then
c ====granules add
c     ioccvirt=nval*nvirt ! number of bins
c     nslv=1
c     call dynamic_getmem((ioccvirt/i4size+1)*1,ipair2gran)
c     call dynamic_getmem((ioccvirt/i4size+1)*2,igran2pair)
c     call granules_calculate(ioccvirt,ncf,byt8,lbin,bl(igran2pair),
c    *                        bl(ipair2gran),igranules,igranulesize,
c    *                        nslv)
c ====granules stop
      indexmem=nval*nvirt*lbin/i4size+1
      if (byt8) then
      call dynamic_getmem(nval*nvirt*lbin,ibin8) !memory for bins
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin8),nval*nvirt*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=nval*nvirt*lbin/i4size+1
      int1mem=nval*nvirt*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      if (iprnt.ge.3) then
      write(6,*)'mem,npairs,lbin,nbins',memleft,npairs,lbin,ncf*ncf/lbin
      call flush(6)
      endif
      else !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e','tt') then
      indexmem=nval*nval*lbin/i4size+1
      if (byt8) then
      call dynamic_getmem(nval*nval*lbin,ibin8) !memory for bins
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin8),nval*nval*lbin)
      call izeroit(bl(indexadr),intsize*indexmem)        
      ibin4adr=1
      ibin1adr=1
c
      else
      int4mem=nval*nval*lbin/i4size+1
      int1mem=nval*nval*lbin/i1size+1
      call dynamic_getmem(int4mem,ibin4adr)  !memory for bins
      call dynamic_getmem(int1mem,ibin1adr)  !memory for overflow
      call dynamic_getmem(indexmem,indexadr) !memory for AO indices
      call zeroit(bl(ibin4adr),int4mem)
      call zeroit(bl(ibin1adr),int1mem)
      call izeroit(bl(indexadr),intsize*indexmem)
      ibin8=1
      endif
      if (iprnt.ge.3) then
      write(6,*)'mem,npairs,lbin,nbins',memleft,npairs,lbin,ncf*ncf/lbin
      call flush(6)
      endif
      endif !if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then

c .................................................
      if (byt8) then
      lrec=12*lbin+8/intsize
      else
      lrec=9*lbin+8/intsize
      endif
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      if (Int_kind.eq.'c') then
      filname1=scrfile(1:len)//'.coul.bins'
      len1 = len+10
      ndisk=41
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      else if (Int_kind.eq.'x') then
      ndisk=42
      filname1=scrfile(1:len)//'.exc.bins'
      len1 = len+9
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      else if (Int_kind.eq.'e') then
      ieeocounter=ieeocounter+1
      ndisk=43
      filname1=scrfile(1:len)//'.eeo.bins'
      len1 = len+9
      if (ieeocounter.gt.1) then 
              close(UNIT=ndisk,STATUS='DELETE')
       endif
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      else if (Int_kind.eq.'tx') then
      txcounter=txcounter+1
      ndisk=44
      filname1=scrfile(1:len)//'.tx.bins'
      len1 = len+8
      if (txcounter.gt.1) then 
              close(UNIT=ndisk,STATUS='DELETE')
       endif
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      else if (Int_kind.eq.'tc') then
      tccounter=tccounter+1
      ndisk=45
      filname1=scrfile(1:len)//'.tc.bins'
      len1 = len+8
      if (tccounter.gt.1) then 
              close(UNIT=ndisk,STATUS='DELETE')
       endif
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      else if (Int_kind.eq.'tt') then
      if (byt8) then
      lrectt=12*lbin+8/intsize
      else
      lrectt=9*lbin+8/intsize
      endif
      ttcounter=ttcounter+1
      ndisk=46
      filname1=scrfile(1:len)//'.tt.bins'
      len1 = len+8
      if (ttcounter.gt.1) then 
              close(UNIT=ndisk,STATUS='DELETE')
       endif
      OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      endif
#endif
! From Alan Ford copy END
CTJ
c--------------------------------------------------------------------
      call dynamic_show_free(memavail)
      if (iprnt.ge.3) write(6,*) 'dynamic memory left:  ', memavail
      icont=0
c--------------------------------------------------------------------
c 2*merged_size*npairs = matrix with results in EEO
c (ipass_size=30)*ncf*ncf
c
      if (Int_kind.eq.'e') then
        max_merg_shell_size=keepmax
        if (small) then
        if (.not.nodisk) then
C allocation for experimental shells merging
          memavail=memavail-ncf*ncf*maxshell ! for integral store
          memavail=memavail-1 500 000 ! for internal reserv in integr.  proc.
          memavail=memavail-30*ncf*ncf ! see above - amplitudes
          icompr_fact=3                ! approximate compression factor,
                                       ! it increases real max_merg_shell_size
          max_merg_shell_size=memavail/(ncf*ncf+2*npairs*icompr_fact)
          if (max_merg_shell_size.lt.keepmax) then
            write(6,*) 'Estimated number of max_merg_shell_size was:',
     *               max_merg_shell_size
            STOP 'Too less memory for EEO operators'
          endif
          if ((max_merg_shell_size-keepmax).gt.npairs) then
            nodisk=.true.  ! temporary in this procedure and subproc. only
            max_merg_shell_size=keepmax
            write(6,*)  'All amplitudes will be kept in memory'
          else if (max_merg_shell_size.gt.max_merg) then
            max_merg_shell_size=max_merg
          endif
        endif
        call dynamic_matdef('screen','q',ncf,ncf)
        iscr=mataddr('screen')
        call matzero('screen')
        call GenAmplSymmAnti(npairs,ncf,     nval,   nodisk,  nfirst,
     *                       nmo,   vorb,    .false.,'master',ndisks,
     *                       ndiska,bl(iscr),list,  isympairs,ij_unique)
c       call dynamic_matrem('screen')
        call dynamic_getmem(ntrisymm*max_merg_shell_size,ibasesymm)
        call dynamic_getmem(ntrianti*max_merg_shell_size,ibaseanti)
        write(6,*) 'max_merg_shell_size=  ',max_merg_shell_size
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
        endif
        filname1=scrfile(1:len)//'.eeo.resi'
        len1 = len+9
        OPEN (UNIT=ndiskresid,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL=ireclength)
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
            write(6,*)  'Reducing max_merg_shell_size which is: '
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
            if (memory_success) write(6,*)  'All amplit. in memory.'
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
c           write(6,*)  'Reducing max_merg_shell_size which is: ',
c    *                                              max_merg_shell_size
            if (max_merg_shell_size.le.1) STOP 'Too little memory EEO'
          else
            exit
          endif
          enddo
        endif
        write(6,*)  'max_merg_shell_size: ',max_merg_shell_size
        write(6,*)  'max_amplitudes_size: ',max_amplitudes_size
        write(6,*)  'Amplitudes will be read: ',
     *                       ncf*ncf/2/max_merg_shell_size,' times.'
        if (max_merg_shell_size.gt.max_merg)
     *                                  max_merg_shell_size=max_merg
        call dynamic_getmem(ncf*ncf*max_merg_shell_size,ishellsaddr)
        call dynamic_getmem(ncf*ncf*max_amplitudes_size,iamplitaddr)
        call dynamic_getmem(max_merg_shell_size*nval*nval,iresultaddr)
      endif ! if (small) then
      endif
c
C
c unlock irecadr ! REMEMBER TO ALLOCATE APPRIOPRATE AMOUNT OF MEMORY
c                  OUTSIDE OF THIS SUBROUTINE! FRAGILE!
      call dynamic_unlock(irecadr,iii)
      igranulesize=1
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e') then
        do iii=1,npairs
        irecadr(iii)=0
        enddo
      else if (int_kind.eq.'tt') then
        do iii=1,nval*nvirt
        irecadr(iii)=0
        enddo
      else
        do iii=1,nval*nval
        irecadr(iii)=0
        enddo
      endif
c
      ncf2=ncf*ncf
      nskipped=0
      skip=0.0d0
      tot=0.0d0
      do i=1,ixint_time
        xint_time(i)=0.0d0
      enddo
      do ics=ncs,1,-1
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf*ncf*ics_size
         do jcs=ics,1,-1
c        if (int_kind.ne.'e'.and.image(ics,jcs)) cycle
         if (image(ics,jcs)) cycle
c
           call get_shell_size(bl(ictr),jcs,jcs_size)
c
           call secund(tt2)
           call elapsec(telap2)
c
           lmp2_size=lmp2_siz1*jcs_size
c
c check if a size of the lmp2 integral buffer is not too big
c if it is  then split over kcs ( and possibly ics)
c
c
          call check_size2(lmp2_size,ncf,nval,max_merg_shell_size,
     *                     npairs,int_kind,ntimes,nodisk,split)
c
          call set_passes(bl(ictr),ics,jcs,ics_size,jcs_size,
     *                      ntimes,Ipass,Jpass,Itimes,Jtimes)
c
c         write(6,*)' ics,jcs=',ics,jcs,' passes=',itimes,jtimes
c
          do itime=1,itimes
             icf1=ipass(1,itime)
             icf2=ipass(2,itime)
             iatonce=icf2-icf1+1
          do jtime=1,jtimes
             jcf1=Jpass(1,jtime)
             jcf2=Jpass(2,jtime)
             jatonce=jcf2-jcf1+1
c
             lmp2_size=iatonce*jatonce*ncf2
c
            call dynamic_getmem(lmp2_size,lmp2int)
            if (Int_kind.eq.'e') 
     *      then
              do itj=0,lmp2_size-1
                bl(lmp2int+itj)=0.0d0
              enddo
            endif
            call dynamic_mmark
            call mmark
            if (Int_kind.eq.'c'.or.Int_kind.eq.'tc'.or.Int_kind.eq.'tt')
     *      then
c           call int_coul(bl, bl(ictr), thresh, ics, icf1,
c    1                    icf2, jcs, jcf1, jcf2, bl(mapf2s),
c    2                    bl(idics),iprnt,bl(lmp2int),nintotal,nrow,
c    3                    ncol, bl(irow),bl(icol),bl(lzero))
            call new_coul_int(bl(ictr), thresh,  ics,     icf1,   icf2,
     *                       jcs,  jcf1,   jcf2,  bl(mapf2s), bl(idics),
     *                        iprnt,bl(lmp2int), nintotal, nrow,   ncol,
     *                        bl(irow),  bl(icol), bl(lzero), ncf, ncs,
     *                        bl(ischwarz),skip,tot,xint_time)
            else if (Int_kind.eq.'x' .or.  Int_kind.eq.'e'
     *                                      .or.Int_kind.eq.'tx') then
c           call int_lmp2(bl,bl(ictr),thresh,
c    1                    ics,icf1,icf2,jcs,jcf1,jcf2,
c    2                 bl(mapf2s),bl(idics),iprnt,bl(lmp2int),nintotal,
c    3                    nrow,   ncol,    bl(irow),bl(icol),bl(lzero))
            call new_exch_int(bl(ictr), thresh,  ics,  icf1,   icf2,
     *                        jcs,  jcf1,  jcf2,  bl(mapf2s), bl(idics),
     *                        iprnt, bl(lmp2int), nintotal, nrow,  ncol,
     *                        bl(irow), bl(icol), bl(lzero), ncf, ncs, 
     *                        bl(ischwarz),skip,tot,xint_time)
            endif
            call retmark
            call dynamic_retmark
          if(nintotal.eq.0) then
            nskipped=nskipped+1
            call dynamic_retmem(1)
            cycle
          else
            ijtot=ijtot+1
            itott=itott+nintotal
          endif
      call secund(tt3)
      call elapsec(telap3)
      tinteg=tinteg+tt3-tt2
      elapint=elapint+telap3-telap2
c......................................................................
      if (Int_kind.eq.'e') then
       reset=.false.
       call EEOGenerator (ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),lmp2int,bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), nrow,
     3                    ncol,
     4                             lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .FALSE.,reset,0,         ibasesymm,ibaseanti,
     8               max_merg_shell_size,nodisk,byt8,small,ndiskresid,
     9          max_amplitudes_size,ishellsaddr,iamplitaddr,iresultaddr,
     *          igranules,igranulesize,igran2pair,ipair2gran,split,
     1          bl(iscr),bl(mapf2s),ik_capacity,tot_xmax,list,nstrong,
     2          isympairs,ij_unique)
      else if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'tt')
     *    then
       call TransOneShell_tj(ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),bl(lmp2int),bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), bl(1),mulam, mulamd, 
     *                    nrow,
     3                    ncol, bl(irow), bl(icol), bl(irow1),bl(icol1),
     4                    ENonZero,smal,lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .FALSE.,' ',0,int_kind,byt8,nmo,vorb,igranules,
     8             igranulesize,igran2pair,ipair2gran,irec_no,irec_last,
     *                  nbf)
       else
       call TransOneShell_tj(ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),bl(lmp2int),bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalfmulam), bl(ihalflammu),
     *                    mulam, mulamd,  nrow,
     3                    ncol, bl(irow), bl(icol), bl(irow1),bl(icol1),
     4                    ENonZero,smal,lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .FALSE.,' ',0,int_kind,byt8,nmo,vorb,igranules,
     8             igranulesize,igran2pair,ipair2gran,irec_no,irec_last,
     *                  nbf)
      endif
c
      call secund(tt4)
      call elapsec(telap4)
      ttrans=ttrans+tt4-tt3
      elaptrans=elaptrans+telap4-telap3
c......................................................................
      call dynamic_retmem(1)          ! lmp2int
c......................................................................
       enddo    ! over jtime (selected jcf belonging to jcs shell )
       enddo    ! over itime (selected icf belonging to ics shell )
c......................................................................
        end do     !  over kcs shell
      end do       !  over ics shell
      if (Int_kind.eq.'e') then
       reset=.true.
       call EEOGenerator (ncf,   ncs,      nval,       ics,   jcs,
     *                    icf1,  icf2,     jcf1,       jcf2,  ndisk,
     1                    bl(ictr),lmp2int,bl(ioccu),iprnt, thresh,
     2                    bl(ixadr),bl(ihalftra), nrow,
     3                    ncol,
     4                             lbin, npairs, bl(icountr),
     5               irecadr,bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     6               bl(ibin8),
     7               .FALSE.,reset,0,         ibasesymm,ibaseanti,
     8               max_merg_shell_size,nodisk,byt8,small,ndiskresid,
     9          max_amplitudes_size,ishellsaddr,iamplitaddr,iresultaddr,
     *          igranules,igranulesize,igran2pair,ipair2gran,split,
     1          bl(iscr),bl(mapf2s),ik_capacity,tot_xmax,list,nstrong,
     2          isympairs,ij_unique)
       endif
      if (Int_kind.eq.'e'.and.(.not.nodisk)) then!files where symm and antysymm 
c                                                 coefficients are written
        if (small) then
        close(ndisks,STATUS='delete')
        close(ndiska,STATUS='delete')
        else
        close(ndiskampl,STATUS='delete')
        endif
      endif
      nodisk=keepnodisk
ctj We need to write down leftovers. For CC all coulomb integrals will
c   be used many times during calculation and we have to have them all
c   written to disk
      if (int_kind.eq.'c'.or.int_kind.eq.'x'.or.int_kind.eq.'e') then
      call WriteLeftOvers(ndisk,npairs,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.false.,0)
      else if (int_kind.eq.'tt') then
      call ext3_WriteLeftOvers(ndisk,npairs,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.false.,0,nval,ncf,nmo,
     3                    igranules,igranulesize,igran2pair,ipair2gran,
     *                    irec_no,irec_last,nbf)
      else
      call WriteLeftOversExt(ndisk,npairs,nval,lbin,bl(icountr),irecadr,
     1                    bl(indexadr),bl(ibin4adr),bl(ibin1adr),
     2                    bl(ibin8),byt8,.false.,0)
      endif

c-------------------------------------------------------------------
      if (iprnt.ge.3) then
      if(nskipped.gt.0) then
        write(iout,*) nskipped,
     1   ' pairs were skipped because no integrals were calculated'
      end if
      write(6,*) 'total amount of integrals: ', itott
      write(6,*) 'total would be calc: ',tot
      write(6,*) 'skipped:             ',skip
      write(6,12)'Time spent in rys,zeroing: ',xint_time(1),
     *           '                  getmem:  ',xint_time(2),
     *           '                  int:     ',xint_time(3),
     *           '                  put res: ',xint_time(4),
     *           '                  irow,icol',xint_time(5),
     *           'Total time inside new int: ',xint_time(6)
      endif
 12   FORMAT (6(A,F10.4,/))
c Added transformation here BEGIN
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
      lrec=idimen*idimen*8
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                 ! 80 is starting lenght
      if      (int_kind.eq.'x') then
        ndisktr=98
        filename=scrfile(1:len)//'.Kij'
        len1=len+4
      else if (int_kind.eq.'c') then
        ndisktr=63
        filename=scrfile(1:len)//'.Cij'
        len1=len+4
      else if (int_kind.eq.'e') then
        ndisktr=64
        filename=scrfile(1:len)//'.Eij'
        len1=len+4
        lrec=ncf*ncf*8
      else if (int_kind.eq.'tx') then
        ndisktr=65
        filename=scrfile(1:len)//'.TXij'
        len1=len+5
      else if (int_kind.eq.'tc') then
        ndisktr=66
        filename=scrfile(1:len)//'.TCij'
        len1=len+5
      else if (int_kind.eq.'tt') then
        ndisktr=73
        filename=scrfile(1:len)//'.TTij'
        len1=len+5
      endif
c     fileopen=.false.
c     inquire(ndisktr,opened=fileopen)
c     if (fileopen) STOP 'Error in subroutine GenCoul, ndisktr opened.'
      close(ndisktr,status='DELETE')
      OPEN (UNIT=ndisktr,FILE=filename(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
c
      irecskip=npairs
      if (int_kind.eq.'tx'.or.int_kind.eq.'tc') irecskip=nval*nval
      if (int_kind.eq.'tt') irecskip=nval*(nbf-nmo)
      islaves=1  ! tmp single
      af=.false. ! tmp single
      do ij=1,irecskip
      call transformator(ij,islaves,irecadr,irecskip,lbin,bl(indexadr),
     *         bl(ibin1adr),bl(ibin4adr),bl(ibin8),ndisk,ndisktr,npairs,
     *         ncf,bl(ixadr),thresh,af,byt8,nmo,vorb,int_kind,nbf)
c       if (Int_kind.eq.'tt') then
c      call transformatgr(ij,islaves,irecadr,irecskip,lbin,bl(indexadr),
c    *        bl(ibin1adr),bl(ibin4adr),bl(ibin8),ndisk,ndisktr,npairs,
c    *        ncf,bl(ixadr),thresh,af,byt8,nmo,vorb,int_kind,
c    *   igranules,igranulesize,bl(igran2pair),bl(ipair2gran),list_copy)
c       else
c      call transformator(ij,islaves,irecadr,irecskip,lbin,bl(indexadr),
c    *        bl(ibin1adr),bl(ibin4adr),bl(ibin8),ndisk,ndisktr,npairs,
c    *        ncf,bl(ixadr),thresh,af,byt8,nmo,vorb,int_kind)
c       endif
      enddo
c
c irecadr readonly:
      call dynamic_lock(irecadr,iii)
c
c     call transformator(ij,islaves,irec,irecskip,lbin,indexx,
c    *                   ibin1,ibin4,bin8,ndisk,ndisktr,npairs,
c    *                   ncf,xmat,thresh,af,byt8,nmo,vorb)
c Added transformation here END
c......................................................................
c timing here is irrelevant :
      call secund(txxx2)
      txxx0=txxx1
      call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
c
      if(iprnt.gt.3) then
         if(nsym.gt.0) then
            write(iout,*) 'Number of shell pairs omitted by symmetry=',
     1       inegl, ' retained=',iret,'transformed',iktot
         endif
      endif
c......................................................................
      nrec=nrec+irec
      if(iprnt.ge.2) then
         write(iout,*) ' Total number of records written=',nrec
      endif
c
      if(iprnt.ge.2) then
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
      call retmark ! mmark just before int_coul
      call dynamic_retmark ! mmark just before int_coul
      if (int_kind.eq.'x'.or.int_kind.eq.'c'.or.int_kind.eq.'e'.or.
     *    int_kind.eq.'tt') then
      call dynamic_matrem('halftra')
      else
      call dynamic_matrem('halflammu')
      call dynamic_matrem('halfmulam')
      endif
      call retmark ! mmark after irec allocation
      call dynamic_retmark ! mmark after irec allocation
c......................................................................
      if (iprnt.ge.2) then
      write(iout,61)
   61 format('  CPU & Elapsed timings in the coulomb integrals module ')
      write(iout,*) '  '
c
      write(iout,62) tinteg/sixty, elapint/sixty
   62 format('  Coulomb Integrals =',f8.2,' and ',f8.2,' minutes')
      write(iout,64) ttrans/sixty, elaptrans/sixty
   64 format('  First half transf. =',f8.2,' and ',f8.2,' minutes')
      call flush(iout)
      endif
ckw---------------------------------------------------------
      if(.not.restrt .and. iprnt.ge.5) then
         call getrval('timeblks',timeblks)
         write(91,123) timeblks/sixty
  123    format('CPU time for new blocking =',f8.2)
      endif
ckw---------------------------------------------------------
      call secund(tt1)
      call elapsec(elaps1)
      if(iprnt.ge.2) then
        write(iout,*) ' CPU and elapsed time for integrals+first half ',
     1  'transformation ',(tt1-tt0)/sixty,(elaps1-elaps0)/sixty
      end if
c..................................................
c---------------------------------------------------------------------
      call symmon
c---------------------------------------------------------------------
      call matremark
      call retmark
      call dynamic_retmark
      end
c---------------------------------------------------------------------
c===============================================================================
      subroutine new_coul_int(inx,   thresh,  ics,     icf1,   icf2,
     *                        jcs,  jcf1,    jcf2,     mapf2s, DS,
     *                        iprnt, xmp2int, nintotal, nrow,   ncol,
     *                        irow,  icol,    lzero,    ncf, ncs, 
     *                        schwarz,skip,tot,xint_time)
c This subroutine builds coulomb matrices C(ij), which are put into
c array xmp2int(ncf(l index),ncf(k index),jcf1:jcf2,icf1:icf2)
c INTENT(IN)
c inx             - basis set contraction info
c thresh          - threshold for integrals, used by schwarz inequality and for 
c                   re-scaling integrals 
c ics             - number of contracted shell for i index of coulomb matrix
c icf1, icf2      - if splitting is in action those indices indicate
c                   subset of basis functions for ics shell for which
c                   integrals are calculated. If splitting is switched off,
c                   those parameters provide the same information as you
c                   can get from inx - first and last cf for given ics shell
c jcs,jcf1,jcf2   - as for "i"
c mapf2s          - not used at present, contains mapping mapf2s(icf)->ics
c DS              - screening density, not used
c iprnt           - print level
c ncf,ncs         - n(umber) of c(ontraced) f(unctions) - ncf or s(hells) - ncs
c schwarz         - SQUARE symmetrix array which contains schwarz(i,j)=(ij|ij)
c INTENT(OUT)
c xmp2int         - integrals array
c nintotal        - total number of integrals calculated (including repetitions
c                   from permutations (ij|kl) and (ij|lk) )
c nrow,ncol       - number of nonzero rows & columns for xmp2int
c irow,icol       - irow,icol arrays mapping new idex to old index in xmp2int
c lzero           - not used
c skip            - number of integrals neglected from schw. inequality
c                   (permutationally unique)
c tot             - total number of permutationally unique integrals calculated
c xint_time       - The array which gather different timings from this subr
      use memory
      implicit none
      integer ics,icf1,icf2,jcs,jcf1,jcf2,mapf2s(*),iprnt,nintotal
      integer nrow,ncol,irow(*),icol(*),lzero(*),ncf
      integer inx(12,*),ncs
      real*8 thresh,DS(*),xmp2int(ncf,ncf,jcf1:jcf2,icf1:icf2)
      real*8 xint_time(*)
      real*8 schwarz(ncs,ncs),skip,tot
c
      integer kcs,lcs,integr,ibuf,isize,jsize,ksize,lsize,ind
      integer i,j,k,l,ij,kl,iff,jff,kff,lff,istart,istop,kstart,kstop
      integer jstart,jstop,lstart,lstop,ityp,jtyp,ktyp,ltyp
      integer igetival,ishift,jshift,kshift,iind,jind,kkind
      integer ijsize,ijksize,ijklsize
      real*8 zero,perm,threshinv,eps,xint,t0,t1,perm_start
      logical nonz,zerro
      parameter (zero=0.0d0)
      integer ncfc,jcf1c,jcf2c,icf1c,icf2c
      real*8 threshc,xstart,xstop
      integer idummy
      common /putting_integrals/isize,jsize,ksize,lsize,perm,idummy,
     *                          ncfc,jcf1c,jcf2c,icf1c,icf2c,threshc,
     *                          istart,istop,jstart,jstop,
     *                          kstart,kstop,lstart,lstop
c
      call secund(xstart)
      ncfc =ncf
      jcf1c=jcf1
      jcf2c=jcf2
      icf1c=icf1
      icf2c=icf2
      threshc=thresh
      threshinv=1.0d0/thresh
      ncs=igetival('ncs')
      nintotal=0
      eps=thresh*1.d-1
      if (ics.eq.jcs) then 
        perm_start=2.0d0
      else
        perm_start=1.0d0
      endif
c
      call dynamic_mmark
c
      call secund(t0)
      call izeroit(irow,ncf)
      call izeroit(icol,ncf)
      call izeroit(lzero,ncf)
      call zeroit(xmp2int,ncf*ncf*(jcf2-jcf1+1)*(icf2-icf1+1))
      call secund(t1)
      xint_time(1)=xint_time(1)+t1-t0
      call secund(t0)
      call dynamic_getmem(15**4,ibuf) ! for (GG|GG) cartesian buffer
      call dynamic_getmem(15**4,integr)
      call secund(t1)
      xint_time(2)=xint_time(2)+t1-t0
      isize=inx(3,ics)
      jsize=inx(3,jcs)
      istart=inx(11,ics)+1
      istop =inx(10,ics)
      jstart=inx(11,jcs)+1
      jstop =inx(10,jcs)
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
      ij=ics*(ics-1)/2+jcs ! ASSUME!!! that ics>=jcs
      ijsize=isize*jsize
      if (ityp.ge.9.and.ityp.ne.11.or.jtyp.ge.9.and.jtyp.ne.11) then
        write(6,*) 'H & I functions not supported yet'
        STOP    'H & I functions not supported yet'
      endif
      if (istart.gt.icf1.or.istop.lt.icf2.or.
     *    jstart.gt.jcf1.or.jstop.lt.jcf2) then
        write(6,*) 'Wrong icf1 icf2 jcf1 jcf2 indices'
        STOP    'Wrong icf1 icf2 jcf1 jcf2 indices'
      endif
      do kcs=1,ncs
        ktyp=inx(12,kcs)
        if (ktyp.ge.9.and.ktyp.ne.11) then
          write(6,*) 'H & I functions not supported yet'
          STOP    'H & I functions not supported yet'
        endif
        ksize=inx(3,kcs)
        kstart=inx(11,kcs)+1
        kstop =inx(10,kcs)
        ijksize=ijsize*ksize
        do lcs=1,kcs            ! true for coulomb only!
          ltyp=inx(12,lcs)
          if (ltyp.ge.9.and.ltyp.ne.11) then
            write(6,*) 'H & I functions not supported yet'
            STOP    'H & I functions not supported yet'
          endif
          lsize=inx(3,lcs)
          lstart=inx(11,lcs)+1
          lstop =inx(10,lcs)
          ijklsize=ijksize*lsize
          tot=tot+ijklsize
          if (sqrt(schwarz(ics,jcs)*schwarz(kcs,lcs)).lt.thresh) then 
            skip=skip+ijklsize
            cycle
          endif
          call array_files(22) ! poll AF
          call secund(t0)
          call twoelmain(ics,jcs,kcs,lcs,eps,nonz,bl(integr),bl(ibuf))
          call secund(t1)
          xint_time(3)=xint_time(3)+t1-t0
          if (.not.nonz) then 
            cycle
          endif
          if (kcs.ge.lcs) then 
            kl=kcs*(kcs-1)/2+lcs
          else 
            kl=lcs*(lcs-1)/2+kcs
          endif
c
          perm=perm_start
          if (kcs.eq.lcs) perm=perm*2.0d0
          if (ij.eq.kl) perm=perm*2.0d0
c
          ishift=lsize*ksize*jsize
          jshift=lsize*ksize
          kshift=lsize
c
          call secund(t0)
          call put_int_results(bl(integr),xmp2int,nintotal,irow,icol)
          call secund(t1)
          xint_time(4)=xint_time(4)+t1-t0
c
        enddo
      enddo
c
c Now compact irow and icol
      nrow=0
      ncol=0
      call secund(t0)
      do i=1,ncf
        if(irow(i).gt.0) then
          nrow=nrow+1
          irow(nrow)=i
        end if
        if(icol(i).gt.0) then
          ncol=ncol+1
          icol(ncol)=i
        end if
      end do
      call secund(t1)
      xint_time(5)=xint_time(5)+t1-t0
c
      call dynamic_retmark
c
 23   FORMAT (4I5,3X,F20.15)
      call secund(xstop)
      xint_time(6)=xint_time(6)+xstop-xstart
      end
c===============================================================================
      subroutine put_int_results(xintegr,xmp2int,nintotal,irow,icol)
c This subroutine simply puts the integral array over shells into large
c coulomb matrix, as needed for CCSD or MP2.
c INTENT(IN) common block /putting_integrals/:
c isize,jsize,ksize,lsize 
c              -  sizes of i,j,k,l shells
c perm         - permutation factor defined as (1+delta(i,j))*(1+delta(k,l))*
c                                              (1+delta(ij,kl))
c ncf          - number of contracted functions
c icf1,icf2    - start and stop contr. f. index for coulomb matrices
c jcf1,jcf2    - start and stop contr. f. index for coulomb matrices
c thresh       - integrals threshold, used mainly for re-scaling
c                integrals for backward compatibility with MP2 subroutines
c istart,istop,jstart,jstop,kstart,kstop,lstart,lstop
c              - start and stop contracted function index for shells i,j,k,l
c INTENT(OUT) common block /putting_integrals/:
c xintegr      - the array from callable integral program with integrals
c                for shell quartet (LCS,KCS|JCS,ICS)
c xmp2int      - the array which contains integrals in the order 
c                (l(1:ncf),k(1:ncf)|j(jcf1,jcf2),i(icf1,icf2))=(ij|kl)
c nintotal     - total amount of integrals bigger than thresh which I put 
c                into coul. matrix
c irow,icol    - the arrays indicating the nonzero rows and columns.
c
      implicit none
      integer isize,jsize,ksize,lsize,ncf,jcf1,jcf2,icf1,icf2,nintotal
      integer irow(*),icol(*)
      integer istart,istop,jstart,jstop,kstart,kstop,lstart,lstop
      real*8 xintegr(0:lsize-1,0:ksize-1,0:jsize-1,0:isize-1)
      real*8 thresh,perm,xmp2int(ncf,ncf,jcf1:jcf2,icf1:icf2)
c
      real*8 threshinv,xint
      integer iff,jff,kff,lff,i,j,k,l,iherestop,idummy
      common /putting_integrals/isize,jsize,ksize,lsize,perm,idummy,
     *                          ncf,jcf1,jcf2,icf1,icf2,thresh,
     *                          istart,istop,jstart,jstop,
     *                          kstart,kstop,lstart,lstop
c
      threshinv=1.0d0/thresh
c
          do iff=icf1,icf2
            i=iff-istart
            do jff=jcf1,jcf2
              j=jff-jstart
              do kff=kstart,kstop
                k=kff-kstart
                iherestop=min(kff,lstop)
                do lff=lstart,iherestop
                  l=lff-lstart
                  xint=xintegr(l,k,j,i)*perm*threshinv
                  if (dabs(xint).gt.1.0d0) then 
                    nintotal=nintotal+1
                    irow(kff)=1
                    icol(lff)=1
                    xmp2int(lff,kff,jff,iff)=xint
c                   if (lff.ne.kff) then
                      nintotal=nintotal+1
                      irow(lff)=1
                      icol(kff)=1
                      xmp2int(kff,lff,jff,iff)=xint
c                   endif
c                 else
c                   xmp2int(lff,kff,jff,iff)=0.0d0
c                   if (lff.ne.kff) xmp2int(kff,lff,jff,iff)=0.0d0
                  endif
                enddo
              enddo
            enddo
          enddo
      end
c===============================================================================
      subroutine new_exch_int(inx,   thresh,  ics,     icf1,   icf2,
     *                        kcs,  kcf1,    kcf2,     mapf2s, DS,
     *                        iprnt, xmp2int, nintotal, nrow,   ncol,
     *                        irow,  icol,    lzero,    ncf, ncs, 
     *                        schwarz,skip,tot,xint_time)
c This subroutine builds exchange matrices K(ik), which are put into
c array xmp2int(ncf(l index),ncf(j index),kcf1:kcf2,icf1:icf2)
c INTENT(IN)
c inx             - basis set contraction info
c thresh          - threshold for integrals, used by schwarz inequality and for 
c                   re-scaling integrals 
c ics             - number of contracted shell for i index of coulomb matrix
c icf1, icf2      - if splitting is in action those indices indicate
c                   subset of basis functions for ics shell for which
c                   integrals are calculated. If splitting is switched off,
c                   those parameters provide the same information as you
c                   can get from inx - first and last cf for given ics shell
c kcs,kcf1,kcf2   - as for "i"
c mapf2s          - not used at present, contains mapping mapf2s(icf)->ics
c DS              - screening density, not used
c iprnt           - print level
c ncf,ncs         - n(umber) of c(ontraced) f(unctions) - ncf or s(hells) - ncs
c schwarz         - SQUARE symmetrix array which contains schwarz(i,j)=(ij|ij)
c INTENT(OUT)
c xmp2int         - integrals array
c nintotal        - total number of integrals calculated (including repetitions
c                   from permutations (ij|kl) and (ij|lk) )
c nrow,ncol       - number of nonzero rows & columns for xmp2int
c irow,icol       - irow,icol arrays mapping new idex to old index in xmp2int
c lzero           - not used
c skip            - number of integrals neglected from schw. inequality
c                   (permutationally unique)
c tot             - total number of permutationally unique integrals calculated
c xint_time       - The array which gather different timings from this subr
      use memory
      implicit none
      integer ics,icf1,icf2,kcs,kcf1,kcf2,mapf2s(*),iprnt,nintotal
      integer nrow,ncol,irow(*),icol(*),lzero(*),ncf
      integer inx(12,*),ncs
      real*8 thresh,DS(*),xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
      real*8 xint_time(*)
      real*8 schwarz(ncs,ncs),skip,tot
c
      integer jcs,lcs,integr,ibuf,isize,jsize,ksize,lsize,ind
      integer i,j,k,l,ij,kl,iff,jff,kff,lff,istart,istop,kstart,kstop
      integer jstart,jstop,lstart,lstop,ityp,jtyp,ktyp,ltyp
      integer igetival,ishift,jshift,kshift,iind,jind,kkind
      integer max_shell_size,iksize,ikjsize,ikjlsize
      real*8 zero,perm,threshinv,eps,xint,t0,t1,perm_start,threshc
      integer ncfc,kcf1c,kcf2c,icf1c,icf2c,idummy
      logical nonz,zerro
      parameter (zero=0.0d0,max_shell_size=15)
      common /putting_integrals/isize,jsize,ksize,lsize,perm,idummy,
     *                          ncfc,kcf1c,kcf2c,icf1c,icf2c,threshc,
     *                          istart,istop,jstart,jstop,
     *                          kstart,kstop,lstart,lstop
c
      ncfc =ncf
      kcf1c=kcf1
      kcf2c=kcf2
      icf1c=icf1
      icf2c=icf2
      threshc=thresh
c
      threshinv=1.0d0/thresh
      ncs=igetival('ncs')
      nintotal=0
      eps=thresh*1.d-1
c
      call dynamic_mmark
c
      call secund(t0)
      call izeroit(irow,ncf)
      call izeroit(icol,ncf)
      call izeroit(lzero,ncf)
      call zeroit (xmp2int,ncf*ncf*(kcf2-kcf1+1)*(icf2-icf1+1))
      call secund(t1)
      xint_time(1)=xint_time(1)+t1-t0
      call secund(t0)
      call dynamic_getmem(max_shell_size**4,ibuf)!for (GG|GG) cartesian buffer
      call dynamic_getmem(max_shell_size**4,integr)
      call secund(t1)
      xint_time(2)=xint_time(2)+t1-t0
      isize=inx(3,ics)
      ksize=inx(3,kcs)
      istart=inx(11,ics)+1
      istop =inx(10,ics)
      kstart=inx(11,kcs)+1
      kstop =inx(10,kcs)
      ityp=inx(12,ics)
      ktyp=inx(12,kcs)
      iksize=isize*ksize
      if (ityp.ge.9.and.ityp.ne.11.or.ktyp.ge.9.and.ktyp.ne.11) then
        write(6,*) 'H & I functions not supported yet'
        STOP    'H & I functions not supported yet'
      endif
      do jcs=1,ncs
        jtyp=inx(12,jcs)
        if (jtyp.ge.9.and.jtyp.ne.11) then
          write(6,*) 'H & I functions not supported yet'
          STOP    'H & I functions not supported yet'
        endif
        jsize=inx(3,jcs)
        jstart=inx(11,jcs)+1
        jstop =inx(10,jcs)
        ikjsize=iksize*jsize
        if (ics.ge.jcs) then 
          ij=ics*(ics-1)/2+jcs
        else 
          ij=jcs*(jcs-1)/2+ics
        endif
        do lcs=1,ncs
          ltyp=inx(12,lcs)
          if (ltyp.ge.9.and.ltyp.ne.11) then
            write(6,*) 'H & I functions not supported yet'
            STOP    'H & I functions not supported yet'
          endif
          lsize=inx(3,lcs)
          lstart=inx(11,lcs)+1
          lstop =inx(10,lcs)
          ikjlsize=ikjsize*lsize
          tot=tot+ikjlsize
          if (sqrt(schwarz(ics,jcs)*schwarz(kcs,lcs)).lt.thresh) then 
            skip=skip+ikjlsize
            cycle
          endif
          call secund(t0)
          call twoelmain(ics,jcs,kcs,lcs,eps,nonz,bl(integr),bl(ibuf))
          call secund(t1)
          xint_time(3)=xint_time(3)+t1-t0
          if (.not.nonz) then 
            cycle
          endif
          if (kcs.ge.lcs) then 
            kl=kcs*(kcs-1)/2+lcs
          else 
            kl=lcs*(lcs-1)/2+kcs
          endif
c
          perm=1.0d0
          if (kcs.eq.lcs) perm=perm*2.0d0
          if (ics.eq.jcs) perm=perm*2.0d0
          if (ij.eq.kl) perm=perm*2.0d0
c
c         ishift=lsize*jsize*ksize
c         kshift=lsize*jsize
c         jshift=lsize
c
c
          call secund(t0)
          call put_int_results_ex(bl(integr),xmp2int,nintotal,irow,icol)
          call secund(t1)
          xint_time(4)=xint_time(4)+t1-t0
c
        enddo
      enddo
c
c Now compact irow and icol
      nrow=0
      ncol=0
      call secund(t0)
      do i=1,ncf
        if(irow(i).gt.0) then
          nrow=nrow+1
          irow(nrow)=i
        end if
        if(icol(i).gt.0) then
          ncol=ncol+1
          icol(ncol)=i
        end if
      end do
      call secund(t1)
      xint_time(5)=xint_time(5)+t1-t0
c
      call dynamic_retmark
c
 23   FORMAT (4I5,3X,F20.15)
      end
c===============================================================================
      subroutine put_int_results_ex(xintegr,xmp2int,nintotal,irow,icol)
c This subroutine simply puts the integral array over shells into large
c exchange matrix, as needed for CCSD or MP2.
c INTENT(IN) common block /putting_integrals/:
c isize,jsize,ksize,lsize 
c              -  sizes of i,j,k,l shells
c perm         - permutation factor defined as (1+delta(i,j))*(1+delta(k,l))*
c                                              (1+delta(ij,kl))
c ncf          - number of contracted functions
c icf1,icf2    - start and stop contr. f. index for exchange matrices
c kcf1,kcf2    - start and stop contr. f. index for exchange matrices
c thresh       - integrals threshold, used mainly for re-scaling
c                integrals for backward compatibility with MP2 subroutines
c istart,istop,jstart,jstop,kstart,kstop,lstart,lstop
c              - start and stop contracted function index for shells i,j,k,l
c INTENT(OUT) common block /putting_integrals/:
c xintegr      - the array from callable integral program with integrals
c                for shell quartet (LCS,KCS|JCS,ICS)
c xmp2int      - the array which contains integrals in the order 
c                (l(1:ncf),j(1:ncf)|k(kcf1,kcf2),i(icf1,icf2))=(ij|kl)
c nintotal     - total amount of integrals bigger than thresh which I put 
c                into exch. matrix
c irow,icol    - the arrays indicating the nonzero rows and columns.
c
      implicit none
      integer isize,jsize,ksize,lsize,ncf,kcf1,kcf2,icf1,icf2,nintotal
      integer irow(*),icol(*)
      integer istart,istop,jstart,jstop,kstart,kstop,lstart,lstop
      real*8 xintegr(0:lsize-1,0:ksize-1,0:jsize-1,0:isize-1)
      real*8 thresh,perm,xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c
      real*8 threshinv,xint
      integer iff,jff,kff,lff,i,j,k,l,idummy
      common /putting_integrals/isize,jsize,ksize,lsize,perm,idummy,
     *                          ncf,kcf1,kcf2,icf1,icf2,thresh,
     *                          istart,istop,jstart,jstop,
     *                          kstart,kstop,lstart,lstop
c
      threshinv=1.0d0/thresh
c
          do iff=icf1,icf2
            i=iff-istart
            do kff=kcf1,kcf2
              k=kff-kstart
              do jff=jstart,jstop
                j=jff-jstart
                do lff=lstart,lstop
                  l=lff-lstart
                  xint=xintegr(l,k,j,i)*perm*threshinv
                  if (dabs(xint).gt.1.0d0) then 
                    nintotal=nintotal+1
                    irow(jff)=1
                    icol(lff)=1
                    xmp2int(lff,jff,kff,iff)=xint
                  endif
                enddo
              enddo
            enddo
          enddo
      end

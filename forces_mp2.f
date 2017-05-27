      subroutine force2_mp2(natoms,natom,rhf,ictr)
C    interface for mp2-gradients
C    this subroutine is called from forces if mp2 is the wavefunction.
C    the MP2 correction to the gradient is calculated and added to the
C    HF gradients.
C
C    this subroutine calls mp2_grad which is the main routine for
C    the MP2-gradients.
C
C    Svein Saebo, Fayetteville, AR Summer 2002
C

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      logical rhf
      common /intlim/ limxmem,limblks,limpair
c
      parameter (idft=0)
C
C  currently only closed-shell MP2 gradients are available
      If(.NOT.rhf) Call nerror(1,'MP2 GRADIENT module',
     $   'Sorry - Open-Shell MP2 gradients currently Unavailable',0,0)
C
C  set values appropriate for integral derivatives during MP2
      limxmem=2 000 000
      limblks=300
      limpair=100
C
      np4=4
      np1=1
c
c -- recover values from depository
      ncf=igetival('ncf')
      ncs=igetival('ncs')
      ncore=igetival('ncore')
      call getival('nprint',IPRNT)
      if(IPRNT.GT.0) write(6,*) 'print level = ',iprnt
      call getrval('thresh',thresh)
      lforc2=igetival('lforc2')
C
C  define matrices and reserve memory for SCF-data:
C  epsi - orbital energies, cano - Canonical orbitals,
C  ovla - overlap matrix, fock - Fock matrix
C  ** density matrix defined and read in calling routine **
      call matdef('epsi','d',ncf,ncf)
      call matread('epsi',np4,'eval_rhf')
      call matdef('cano','q',ncf,ncf)
      call matread('cano',np4,'evec_rhf')
      call matdef('ovla','s',ncf,ncf)
      call matread('ovla',np1,'s matrix')
      call matdef('fock','s',ncf,ncf)
      call matread('fock',np4,'fock_rhf')
c
c -- save pointers in depository
      lfock0=mataddr('fock')
      lvec=mataddr('cano')
      lval=mataddr('epsi')
      call setival('lfock0',lfock0)
      call setival('lvec',lvec)
      call setival('lval',lval)
      call getmem(2,locc)
      call rea (bl(locc),1      ,np4,'nocc_rhf')
      nocc=bl(locc)
      nmo=nocc
      call setival('nocc',nocc)
      call retmem(1)
c
c  check if the basis set contains L-shells
c  and if it does then make S,P partitioning
c
      call trans_l_2_sp(rhf,bl,bl(ictr),ictr,lshell,'mp2grad')
c
c ictr is changed on return if l-shells are present
c
      if(ncore.gt.0) call matsub('coro','cano',1,ncore)
      call matsub('occa','cano',1,nmo)
      call matsub('occu','cano',ncore+1,nmo)
      call matsub('virt','cano',nmo+1,ncf)
      nval=nmo-ncore
      nvirt=ncf-nmo
c
      if(IPRNT.GT.0) write(6,*) 'Integral Threshold:',thresh
c-----------------------------------------------------------------
c allocate memory for an array showing if a big bs shell 
c was present (0 or 1) in a small bs
c
      ncs=igetival('ncs')
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c
c this is needed for dual basis set but is also used
c in blocking for ordinary MP2 (must have all 1)
c-----------------------------------------------------------
      call set_mpres(bl(mpres_in_bg),ncs)
c-----------------------------------------------------------------
      call mmark
      call mp2_grad(ncf,    nval,   nvirt,  IPRNT,  thresh,
     2              nmo,    natom,bl(lforc2),ncore)
      call retmark
c-----------------------------------------------------------------
c
c  At this point one- and two-electron contributions to the
c  gradient are known. Put them together and write them to
c  the file <grad>
c
      call force_total(natom,IPRNT,bl)
c
c transfer back original L-shell basis set info
c
      call trans_sp_2_l(bl,bl(ictr),ictr,lshell)
c
c returns the original value of ictr
c
      call matreset
c     call retimark
      call retmark
c ...........................................................
      return
      end
C=======mp2_grad=====================================================
      subroutine mp2_grad(ncf,    nval,   nvir,   iprint, thresh,
     1                    nmo,    natoms, gradv,  ncore)
C
C  Main routine for MP2-gradients.  The residia (Tij) in (virtual)
C  MO basis are on disk, opened via unit ndisk1, stored as 5-byte integers.
C  The virtual-occupied block of the exchange matrices are also on disk,
C  opened via unit ndisk2, also stored as 5-byte integers.
C  The final results for the various contributions are added
C  up in gradv(3,natoms)
C
C     External calls:
C     ATerms    (constructs the A matrix)
C     XWYterms  (constructs the matrices X, W, and Y, except for
C                D1, D2, and D3 contributions)
C     BinSoRev  (A2 contributions Eq.27, Backtransformed amplitudes)
C     D1terms   (matrix Y3, Eq.46, D1-contribution to Y)
C     cphfz      CPHF calculation
C
C     arguments:
C     ncf        number of basis functions
C     nval       number of correlated orbitals
C     nvir       number of virtual irbitals
C     iprint     print level
C     thresh     integral threshold, normally 1.0d-09
C     nmo        number of occupied orbitals
C     natoms     number of real atoms in the system
C     gradv      gradient vector dimension 3*natoms
C     ncore      number of core orbitals
C
C         Svein Saebo, Fayetteville AR and Starkville,MS
C         Summer 2002
C

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      dimension gradv(3,*)
      character*256 scrfile,filname1,filname2,filname3,filname4,filname5
      logical exists
      parameter(sixty=60.0d0,two=2.0d0,onef=0.25d0)
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,four=4.0d0)
c
      call secund(tgr0)
      call elapsec(egr0)
c
c-------------------------------------------------------------------
      call getival('mp2only',mp2only)
      if(mp2only.GT.0) then
        call nerror(13,'MP2 Gradient module',
     $    'Gradient files NOT written in preceeding MP2 step',0,0)
      endif
c-------------------------------------------------------------------
      call mmark
C
      iscs=igetival('iscs')
      if(iscs.eq.1) write(6,*) 'SCS gradient will be calculated'
      if(iscs.eq.2) write(6,*) 'Scaled MP2 gradient will be calculated'
      iout=igetival('iout')
      inuc=igetival('inuc')
      ncs=igetival('ncs')
      ictr=igetival('ictr')
      call setival('printh',2)
      call setival('printh',iprint)
C
      nvrsq=nvir*nvir
      npars=nval*(nval+1)/2
c
c ................................................................
c  determine size of <bins> file for <ATerms>
c
      nindxp=nvir*(nvir+1)/2
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      lcore=igetival('lcore')
      kcore=lcore-ioffset
      memused=lastaddr-ioffset
C  left-  matrices in A1phas -a little extra
      mem=kcore-memused-4*nvrsq
cc
      if(iprint.ge.2) then
        write(iout,*) ' Memory available for sort: ',mem
        write(iout,*) ' Memory assigned for job:   ',lcore-ioffset
      endif
cc
C  Divide available memory into nindxp bins
      lbin=(8*mem)/(10*nindxp) ! length of bin in 10-byte words
      if(lbin.lt.100)
     1   call nerror(2,'MP2 Gradient',
     2                 '<ATerms> Sort bin too small',lbin,npars)
      if(lbin.gt.npars) lbin=npars
C
C  nbins= estimate for the number of bins written
C  (must exceed the true value)
C
      nbinspp=npars/lbin
      nbins=nindxp*(nbinspp+1)
      if(lbin.eq.npars) nbins=nindxp
      nrcpb=nbins/nindxp
cc
      if(iprint.ge.2) then
        write(iout,111) 10*lbin,nbins,10*nindxp*lbin,nindxp
  111 format(' Memory requirements for sort in MO basis:',/,
     2       ' length of one bin:       ',I10,' bytes',/,
     3       ' number of bins:          ',I10,/,
     4       ' memory for bins:         ',I10,' bytes',/,
     5       ' number of index pairs:   ',I10)
      endif
c ................................................................
C
C  get filenames for <Tij> <Kov> and <bins> files
C
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filname1=scrfile(1:len)//'.Tij'
      filname2=scrfile(1:len)//'.Kov'
      filname3=scrfile(1:len)//'.bins'
      filname4=scrfile(1:len)//'.Zai'
      filname5=scrfile(1:len)//'.Gai'
c
      ndisk1 = 41        ! unit number for <Tij> file
      ndisk2 = 42        ! unit number for <Kov> file
      ndisk3 = 43        ! unit number for bins file
c -- ndisk1 reused later for Zai file
c -- ndisk2 reused later for Gai file
C
C  open the <Tij>, <Kov> and <bins> files
      OPEN (UNIT=ndisk1,FILE=filname1(1:len+4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=5*nvir*nvir)
      OPEN (UNIT=ndisk2,FILE=filname2(1:len+4),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=10*nvir*nmo)
      inquire(file=filname3(1:len+5),exist=exists)
      OPEN (UNIT=ndisk3,FILE=filname3(1:len+5),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=10*lbin)
C
C  Calculate Matrix A
      call secund(taik1)
      call elapsec(eaik1)
      call matdef('Aik','s',nval,nval)
      call matzero('Aik')
      call ATerms(nval,   nvir,   ndisk1, ndisk3, iprint,
     $            nvrsq,  nindxp, lbin,   nrcpb,  thresh,
     $            iscs)
      call secund(taik2)
      call elapsec(eaik2)
cc
      if(iprint.ge.2) then
        write(iout,*) ' Construction of Matrix A:'
        write(iout,100) (taik2-taik1)/sixty,(eaik2-eaik1)/sixty
      endif
C
C  the following matrices should be saved to the end
C  matrix W will be contracted with Sx, and X and Aik
C  will be contracted with Fx
C  Y  input for CPHF
C  Zhigang_3/31/2017 @UARK
C  Change the last parameter from ncf to nvir+nmo
C  There was some problem with the dimension of evir when ncf/=nmo+nvir
C      call matsub('evir','epsi',nmo+1,ncf)
      call matsub('evir','epsi',nmo+1,nvir+nmo) !NZG_3/31/2017
      call matsub('eocc','epsi',ncore+1,nmo)
      call matdef('Y','r',ncf,nmo)
      call matsub('Yp','Y',ncore+1,nmo)
      iypadr=mataddr('Yp')
C     call memory_status('Yp')
      call matdef('DDT','s',ncf,ncf)
      call matdef('W','q',ncf,ncf)
      call matdef('X','s',ncf,ncf)
      call matdef('Gmat','s',ncf,ncf)
C  the following matrices are used for temporary storage in
C  XWYterms and can be removed upon exit from this routine
      call matdef('W1','s',ncf,ncf)
      call matdef('W2','s',ncf,ncf)
C
C  part of B-terms but best calculated in XWYterms
      call matdef('B1','r',nvir,nmo)
      call matzero('B1')
C
      call matdef('Ttilda','q',nvir,nvir)
      ittij=mataddr('Ttilda')
      call matdef('Tij','q',nvir,nvir)
      iatij=mataddr('Tij')
      call matdef('Kvo','r',nvir,nmo)
      iovka=mataddr('Kvo')
      call matdef('tvir','r',nvir,ncf)
      call matpose2('virt','tvir','n')
cc
      if(iprint.ge.6) then
        if(ncore.gt.0) then
          call matprint('occa',6)
        else
          call matprint('occu',6)
        endif
        call matprint('virt',6)
      endif
cc
      lbfdim=max0(nvrsq,nmo*nvir*2)
c     call getmem(lbfdim/2+1,ibuf)
      call getint_4(lbfdim,ibuf)
c     call getmem(lbfdim/8+1,i1)
      call getint_1(lbfdim,i1)
      call mmark
C  calculate matrices X, W, A, and Y (except D1-terms)
C  X Eq. 19; A Eq. 21, W Eq. 58 or 60; Y Eq.59 or 61
      call XWYterms(ncf,    nval,   nvir,   ndisk1, ndisk2,
     1              iprint, thresh,bl(iatij),bl(ittij),bl(ibuf),
     2              bl(i1), gradv,  natoms, bl(ibuf),bl(i1),
     3              bl(iovka),ncore,nmo,    iscs)
      call retmark
      call retmem(2)
C  remove matrices for temporary storage
      call matrem('tvir')
      call matrem('Kvo')
      call matrem('Tij')
      call matrem('Ttilda')
      call matrem('B1')
      call matrem('W2')
      call matrem('W1')
C
C  the matrices X, and W, and several contributions to Y
C  are in X, W, Aik, and Y, respectively
C
      if(iprint.ge.6) then
        call matprint('X',6)
        call matprint('W',6)
        call matprint('Aik',6)
      endif
C
      if(iprint.ge.2) then
        call secund(ta1)
        call elapsec(egra1)
        write(iout,*) ' CPU and elapsed time for X, W, and Y-terms:'
        write(iout,100)(ta1-taik2)/sixty,(egra1-eaik2)/sixty
  100   format(1x,1f8.2,' minutes ',1f8.2,' minutes ')
      endif
c
c ..................................................................
c  we have finished with the <Kov> file
      CLOSE (UNIT=ndisk2,STATUS='DELETE')
c  delete current <bins> file
      CLOSE (UNIT=ndisk3,STATUS='DELETE')
c ..................................................................
C
C  construct the two-electron part of the Fock-matrix with DDT
C  as density, result in Gmat.
C  matrix DDT was constructed in XWYterms
      call mmark                                  !gmat
      igadr=mataddr('Gmat')
      idsa=mataddr('DDT')
      call MakeGma2(ncf,nmo,nval,nvir,thresh,
     1              bl,bl(ictr),bl(idsa),bl(igadr))
c
c
      call matscal('Gmat',four)
C
C  multiply with occupied orbitals and add to Y:
      call matmmul2('Gmat','occa','Y','n','n,','a')
cc
      if(iprint.ge.6) then
C       call matprint('Gmat',6)
        write(iout,*) ' Y before D1-terms:'
        call matprint('Y',6)
      endif
cc
      call retmark                                !gmat
      call matrem('Gmat')
      call secund(tgm2)
      call elapsec(egm2)
cc
      if(iprint.ge.2) then
        write(iout,*) ' CPU and elapsed time for construction of Gmat:'
        write(iout,100)(tgm2-ta1)/sixty,(egm2-egra1)/sixty
      endif
C
C  do A2-terms Eq.27
C
C  BinSoRev performs 'reversed' bin sort of the amplitudes
C
C  the fully backtransformed amplitudes Tmylam(ny,sigma) are too
C  numerous to be stored and must be contracted with integral
C  derivatives as soon as they are generated.  The results
C  are added to the gradient vector
C  this means that the two-electron integral derivatives are also
C  calculated here.
C
C  Upon exit from BinSoRev the amplitudes are stored (in bins) as
C  Tmylam(i,j)  there are ncf*(ncf+1)/2 matrices of dimension
C  nval*nval.  These will also be used for the D1-terms (see below)
C
c ................................................................
c  determine size of <bins> file for <BinSoRev>
c
      call mmark
      nindxp=ncf*(ncf+1)/2
      ncfsq=ncf*ncf
      nvrsq=nvir*nvir
      npars=nval*(nval+1)/2
C  check available memory
      call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
      lcore=igetival('lcore')
      kcore=lcore-ioffset
      memused=lastaddr-ioffset
C  left- buffer-matrices in Rphas1 -a little extra
      mem=kcore-memused-3*ncfsq-nvrsq-ncf*nvir
cc
      if(iprint.ge.2) then
        write(iout,*)' Memory available for sort: ',mem
        write(iout,*)' Memory assigned for job:   ',lcore-ioffset
      endif
cc
C  Divide available memory into nindxp bins
      lbin=(8*mem)/(10*nindxp)      ! length of bin in 10-byte words
      if(lbin.lt.100)
     1   call nerror(3,'MP2 Gradient',
     2                 '<BinSoRev> Sort bin too small',lbin,npars)
      if(lbin.gt.npars) lbin=npars
C
C  for ~1000 bf and nval=120: nindxp=500000 and npairs= 14400
C  for each Tij we store 1 integral (no indices needed here)
C  for each indexpair (mylam) we save  npairs integers
C  about 14400 integers.
C  disk-storage: 4*500000*14400 bytes ~29 GB !
C  With 1.6 GB lbin could be ~850 writeout 17 bins per indexpair
C  or 850,000 records, record length 3400 bytes.
C
C  to save disk-space we could sort a limited number of mylam
C  at the time (see binsort for energy calculation)
C  may not be necessary, the above scheme seem to work fine
C
C  nbins= estimate for the number of bins written
C  (must exceed the true value)
C
      nbinspp=npars/lbin
      nbins=nindxp*(nbinspp+1)
      if(lbin.eq.npars) nbins=nindxp
      nrcpb=nbins/nindxp
cc
      if(iprint.ge.2) then
        write(iout,112) 10*lbin,nbins,10*nindxp*lbin,nindxp
  112   format(' Memory requirements for reversed Binsort:',/,
     2  ' length of one bin:       ',I10,' bytes',/,
     3  ' number of bins:          ',I10,/,
     4  ' memory for bins:         ',I10,' bytes',/,
     5  ' number of indexpairs:    ',I10)
      endif
C
C  reopen the <bins> files
      OPEN (UNIT=ndisk3,FILE=filname3(1:len+5),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=10*lbin)
cc
      call BinSoRev(ncf,    nval,   iprint, ndisk1, ndisk3,
     1              thresh, nvir,   natoms, nindxp, lbin,
     2              nrcpb,  ncore,  gradv)
      call retmark
cc
      if(iprint.ge.2) then
        call secund(ta2)
        call elapsec(egra2)
        write(iout,*) ' CPU and elapsed time for binsort, phase 1'
        write(iout,100)(ta2-tgm2)/sixty,(egra2-egm2)/sixty
      endif
c
c ..................................................................
c  we have finished with the <Tij> file
      CLOSE (UNIT=ndisk1,STATUS='DELETE')
c ..................................................................
C
C   calculate D1 terms
C   generates matrix Y3(ny,i) (eq 46)
C   bins still on ndisk3
C
      call mmark                                   !D1terms
c     call getmem(lbin,ibin)
      call getint_4(lbin*2,ibin)
c     call getmem(lbin/4+1,i1)
      call getint_1(lbin*2,i1)
C
      call symmoff
      call D1terms(ndisk3,bl(iypadr),nval,   iprint, thresh,
     1             ncore,  bl(ibin), bl(i1), lbin,   nrcpb,
     2             ncf, nmo,iscs)
      call symmon
      call retmark                                 !D1terms
cc
      if(iprint.ge.2) then
        call secund(td1)
        call elapsec(egrd1)
        write(iout,*) ' CPU and elapsed time  for D1-terms: '
        write(iout,100)(td1-ta2)/sixty,(egrd1-egra2)/sixty
      endif
C
C  calculate Fx and Sx and write to disk
      nfunit=39
      nsunit=40
      call mmark
      call FxSx(natoms,iprint,nfunit,nsunit)
      call retmark
C
C  everything done except CPHF
C
      call secund(tcphf1)
      call elapsec(ecphf1)
C  for frozen core -------------------------------------------------
      if(ncore.gt.0) then
         call matdef('Zic','r',nval,ncore)
         call matdef('Yoo','q',nmo,nmo)
         icic=mataddr('Zic')
         iyoo=mataddr('Yoo')
         if(iprint.ge.6) call matprint('Y',6)
         call matmmul2('occa','Y','Yoo','t','n','n')
         icor=mataddr('epsi')
         call putZic(bl(icic),bl(iyoo),nmo,nval,ncore,bl(icor))
cc
         if(iprint.ge.3) then
            call matprint('Yoo',6)
            call matprint('Zic',6)
         endif
         call matrem('Yoo')
      endif
C  end frozen core stuff -------------------------------------------
C
      call matdef('Zai','r',nvir,nmo)
      call matdef('Yai','r',nvir,nmo)
      izaia=mataddr('Zai')
      iyaia=mataddr('Yai')
      iocca=mataddr('epsi')
      ivira=mataddr('evir')
      call matmmul2('virt','Y','Yai','t','n','n')
C
C  reuse ndisk1 and ndisk2 for Zai and Gz's
C
      open(unit=ndisk1,file=filname4(1:len+4),form='unformatted',
     1     access='direct',recl=8*nvir*nmo)
      open(unit=ndisk2,file=filname5(1:len+4),form='unformatted',
     1     access='direct',recl=8*nvir*nmo)
      call cphfz(ncf,nval,nvir,bl(iocca),bl(ivira),
     1           bl(iyaia),bl(izaia),nmo,thresh,inx,
     2           iprint,ncore,ndisk1,ndisk2)
      close(unit=ndisk1,status='delete')
      close(unit=ndisk2,status='delete')
cc
      call secund(tcphf)
      call elapsec(ecphf)
      write(iout,430) (tcphf-tcphf1)/sixty,(ecphf-ecphf1)/sixty
      call f_lush(iout)
  430 format(/'Master CPU time for   CPHF  solver    = '
     *,f8.2,'  Elapsed = ',f8.2,' min' )
cc
      call matrem('Yai')
      call matrem('Zai')
      if(ncore.gt.0) call matrem('Zic')
C
C  do second phase of reversed binsort here
C
C  this is moved here to avoid recalculating integral derivatives
C  the "A2" contribution will also contain the two-electon part
C  of the HF gradient and the two electron part of the F(x) terms
C
      call matdef('TMO','q',nval,nval)
      itmo=mataddr('TMO')
      call matdef('TAO','q',ncf,ncf)
      itao=mataddr('TAO')
      call matdef('tocc','r',nval,ncf)
      call matpose2('occu','tocc','n')
c
c     call getmem(lbin,ibin)
      call getint_4(lbin*2,ibin)
c     call getmem(lbin/4+1,i1)
      call getint_1(lbin*2,i1)
C
      call symmoff
      call getival('SymFunPr',ifp)
      call Rphas2(ncf,    nval,   ndisk3, nbins,  lbin,
     1            thresh,bl(ibin),bl(i1),bl(itmo),bl(itao),
     2            nsym,  iprint, bl(ictr),nrcpb,bl(ifp),
     3            gradv, iscs,  natoms)
      call symmon
C
      call retmem(2)
      call matrem('tocc')
      call matrem('TAO')
      call matrem('TMO')
C
c ..................................................................
C  Delete the <bins> file
      CLOSE (UNIT=ndisk3,STATUS='delete')
c ..................................................................
cc
      if(iprint.ge.2) then
         call secund(tt2)
         call elapsec(elaps2)
         tt21=tt2-tcphf
         el21=elaps2-ecphf
         t21=tt21/sixty
         e21=el21/sixty
         write(iout,*) ' Binsort Phase2'
         write(iout,100) t21,e21
      endif
cc
      if(iprint.ge.2) then
        Write(iout,*) ' MP2 gradients after A2-terms'
        call torque(NAtoms,0,bl(inuc),gradv )
      endif
C
C  build gradient vector
C
      ntri=ncf*(ncf+1)/2
      call matdef('fxsx','s',ncf,ncf)
      call matdef('fxsy','s',ncf,ncf)
      call matdef('fxsz','s',ncf,ncf)
      ifxsx=mataddr('fxsx')
C
C  do F(x) terms:
C  make it quadratic to simplify trace
      call matdef('XF','q',ncf,ncf)
      call matcopy('DDT','XF')
      ixadr=mataddr('XF')
C
C  add -<X|Fx> to forces:
C  NOTE only one-electron part left of Fx
      call Makegrad(natoms,gradv,bl(ifxsx),nfunit,ntri,
     1             bl(ixadr),ncf)
      call matrem('XF')
cc
      if(iprint.ge.2) then
         Write(iout,*) ' MP2 gradients after X-terms:'
         call torque(NAtoms,0,bl(inuc),gradv )
      endif
C
C  do <SxW> terms
C  subtract 1/4 DYCo to restore orthogonality
C
      if(iprint.ge.3) then
         write(iout,*) ' final W'
         call matprint('W',6)
      endif
cc
      call matdef('DY','r',ncf,nmo)
      call matmmult('den0','Y','DY')
      call matscal('DY',-onef)
      call matmmul2('DY','occa','W','n','t','a')
      call matrem('DY')
      iwad=mataddr('W')
C
C  add <Sx|W> to forces:
C  call matpose('W')
      call Makegrad(natoms,gradv,bl(ifxsx),nsunit,ntri,
     1             bl(iwad),ncf)
C
      call matrem('fxsz')
      call matrem('fxsy')
      call matrem('fxsx')
cc
      if(iprint.ge.3) then
         write(iout,*) ' Final W:'
         call matprint('W',6)
      endif
cc
      if(iprint.ge.2) then
         Write(iout,*) ' MP2 gradients after W-terms:'
         call torque(NAtoms,0,bl(inuc),gradv )
      endif
C
C   before returning calculate the MP2 dipole moments
      call matdef('dip','v',3,3)
      idip=mataddr('dip')
      call mp2dip(bl(idip),ncf,bl(ictr),ncs,natoms)
      call matrem('dip')
c
c
ckw2008kw2008kw2008kw2008kw2008kw2008kw2008kw2008-----
c
c save DMP2 matrix and G(DMP2)
c
      call matdef('GMP2','s',ncf,ncf)
      igmp2=mataddr('GMP2')
      idmp2=mataddr('DDT')
c
c
      call mmark
      call MakeGma2(ncf,nmo,nval,nvir,thresh,
     1              bl,bl(ictr),bl(idmp2),bl(igmp2))
      call retmark
c
c     call saveDmp2(ntri,1)
c
      np4=4
      call matwrite('DDT',np4,0,'dens_mp2')
      call matwrite('GMP2',np4,0,'gock_mp2')
ccc   call matwrite('DDT',np4,0,'den0_rhf')
ccc   call matwrite('GMP2',np4,0,'fock_rhf')
c
      call matrem('GMP2')
c
ckw2008kw2008kw2008kw2008kw2008kw2008kw2008kw2008-----
c
c
C  calculation finished
C  release all memory
C
      call retmark
cc
      if(iprint.ge.2) then
         call secund(tt3)
         call elapsec(et3)
         ttra=(tt3-tt2)/sixty
         etra=(et3-elaps2)/sixty
         write(iout,*)' Time forming gradient vector:'
         write(iout,100) ttra,etra
      endif
C
c ..................................................................
C  delete files with Fx and Sx
      close(unit=nsunit,status='delete')
      close(unit=nfunit,status='delete')
c ..................................................................
C
      return
      end
C======XWYterms=========================================================
      subroutine XWYterms(ncf,    nval,   nvir,   ndisk1, ndisk2,
     1                    iprint, thresh, tij,    ttilda, ibuf,
     2                    i1,     gradv,  natoms, lbuf,   i1bin,
     3                    Kvo,    ncore,  nmo,    iscs)
C
C    calculates the A1 and A3 contributions to the MP2 gradients
C    as well as several contributions to Y
C
C    the A1 contributions ( to be contracted with Fx, eq.15a) are the
C    symmetrical matrices X (eq. 16) and A (eq.18)
C    the A3 contributions ( to be contracted with Sx, eq.15b) are the
C    symmetrical matrices W1 (eq.16)  and W2 (eq 17) as well as some
C    B-type contributons that can be formulated in terms of W, X and
C    Xp
C
C    all matrices are initially generated in (virtual) MO basis, and
C    transformed to AO basis at the end.
C    the Tij in virtual basis are on unit ndisk1
C
C    finally the 'density matrix' used to generate the Fock-matrix for
C    the D2-terms (eq 34) , called 'DDT' is constructed.
C
C    Svein Saebo, Fayetteville, AR summer 2002
C
C   called from mp2grad
C
C    arguments:
C    ncf         number of basis functions
C    nval        number of correlated orbitals
C    nvir        number of virtual orbitals
C    ndisk1      unit for Tij an virtual basis
C    iprint      print level
C    thresh      integral threshold, normally 1.0d-09
C    tij         matrix of dimension nvir*nvir
C    ttilda      matrix of dimension nvir*nvir
C    ibuf        integer matrix of dimension nvir*nvir
C

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30000)
      real*8 Kvo(*)
      dimension tij(*),ttilda(*),gradv(3,natoms)
      integer*4 ibuf(nvir**2),lbuf(nmo*nvir,2)
      integer*1 i1(nvir**2),i1bin(nmo*nvir,2)
      parameter (zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0,four=4.0d0)
      parameter (dblmax=2147483648.0d0)
C
C   all matrices calculated here : X, W, Aik, should be saved
C   to be contracted with Fx and Sx later
C
C    now calculate matrix X (eq.16) W1 (eq 16) W2 (eq 17) Aik (eq.18)
C    and Xp (part of B-terms)
C
C
      dblcmp = dblmax*thresh
c
      call matdef('Tsum1','q',nvir,nvir)
      call matzero('Tsum1')
      call matdef('Tsum2','q',nvir,nvir)
      call matzero('Tsum2')
      call matdef('Tsum4','q',nvir,nvir)
      call matzero('Tsum4')
      call matdef('TT','q',nvir,nvir)
C
      iocca=mataddr('eocc')-1         ! address for orbital energies
C
C    loop over pairs of occupied orbitals: ij
      ttkl=zero
      ij=0
      NTij=0
      NKij=0
c
      do ii=1,nval
         oei=bl(iocca+ii)
         do jj=1,ii
            ij=ij+1
            oej=bl(iocca+jj)
C  read Tij from disk convert to real and put into matrix Tij
            read(ndisk1,rec=ij) i1,ibuf
            do imov=1,nvir**2
c -- decompress the integral ------------------------
               If (i1(imov).eq.0) Then
                  xx = ibuf(imov)*thresh
               Else If(i1(imov).gt.0) Then
                  xx = ibuf(imov)*thresh
                  xx = xx + SIGN(i1(imov)*dblcmp,xx)
               Else
                  xx = ibuf(imov)*thresh*10.0d0**(-i1(imov))
                  write(6,*) ' Threshold Tij-imov:',imov,' i1:',i1(imov)
               EndIf
c ---------------------------------------------------
               Tij(imov)=xx
            enddo
            call matcopy('Tij','Ttilda')
            call atoat2(ttilda,nvir,'y',iscs)
C
C     atoat generates Tbar (not Ttilda) 'y' means transpose
C
C     read virtual-occupied block of Kij for B1-terms
C
            read(ndisk2,rec=ij) i1bin,lbuf
            do imov=1,nvir*nmo
c -- decompress the integral ------------------------
               If (i1bin(imov,1).eq.0) Then
                  xx = lbuf(imov,1)*thresh
               Else If(i1bin(imov,1).gt.0) Then
                  xx = lbuf(imov,1)*thresh
                  xx = xx + SIGN(i1bin(imov,1)*dblcmp,xx)
               Else
                  xx = lbuf(imov,1)*thresh*10.0d0**(-i1bin(imov,1))
cc        write(6,*) ' Threshold Kij - imov:',imov,' i1:',i1bin(imov,1)
               EndIf
c ---------------------------------------------------
               Kvo(imov)=xx
            enddo
C
            call matmmul2('Ttilda','Kvo','B1','n','n','a')
            If (ii.gt.jj) Then
               do imov=1,nvir*nmo
c -- decompress the integral ------------------------
                  If (i1bin(imov,2).eq.0) Then
                     xx = lbuf(imov,2)*thresh
                  Else If(i1bin(imov,2).gt.0) Then
                     xx = lbuf(imov,2)*thresh
                     xx = xx + SIGN(i1bin(imov,2)*dblcmp,xx)
                  Else
                     xx = lbuf(imov,2)*thresh*10.0d0**(-i1bin(imov,2))
cc        write(6,*) ' Threshold Kij - imov:',imov,' i1:',i1bin(imov,2)
                  EndIf
c ---------------------------------------------------
                  Kvo(imov)=xx
               enddo
               call matmmul2('Ttilda','Kvo','B1','t','n','a')
            EndIf
C
C  finished B1 terms  (nvir,nval)
C
            call matmmult('Tij','Ttilda','TT')
            call matadd('TT','Tsum1')  !  Tsum1=sum Tij*Ttilda+
            call matadd1('TT',oei+oej,'Tsum4')
            if (ii.gt.jj) then
               call matmmul2('Tij','Ttilda','TT','t','t','n')
               call matadd('TT','Tsum1')  !  Tsum1=sum Tij*Ttilda+
               call matadd1('TT',oei+oej,'Tsum4')
            endif
C
            call matdef('txx','q',nvir,nvir)
            call matmmult('Tij','evir','txx')
            call matmmult('txx','Ttilda','TT')
            call matadd('TT','Tsum2')
            if (ii.gt.jj) then
               call matmmul2('Tij','evir','txx','t','n','n')
               call matmmul2('txx','Ttilda','TT','n','t','n')
               call matadd('TT','Tsum2')
            endif
            call matrem('txx')
C
C   form matrix Aik:
C     C  SS Sept 2033 Aik now calculated with sorted Ts see Aterms
C     do kk=1,ii
C     ik=iik+kk
C     call gettkj(jj,kk,ndisk1,nvrsq,ibuf,
C    1            mrcpf,bl(itjka),nvir,thresh)
C  note A(ik)=<Tkj Ttildaji>, Ttildaji in ttilda, Tjk bl(itkja)
C  <TkjTji> = ddot(TjkTji)
C
C     bl(iaika+ik)=bl(iaika+ik) +
C    1 ddot(nvrsq,bl(itjka),1,ttilda,1)
C     enddo    !loop over kk
C
         enddo    !loop over jj
      enddo    !loop over ii
C
C  Finally transform to AO basis and save
C
      call matsimtr('Tsum1','tvir','X')
      call matsimtr('Tsum2','tvir','W1')
      call matsimtr('Tsum4','tvir','W2')
cc
      if(iprint.ge.3) then
         call matprint('W1',6)
         call matprint('W2',6)
      endif
C
C  matrices X, W1 and W2 ready in AO basis , remove matrices
C  for temorary storage
C
      call matrem('TT')
      call matrem('Tsum4')
      call matrem('Tsum2')
      call matrem('Tsum1')
C
C  form vector W (ncfxncf)
C
C  W=2W1-2W2-2Cv(B1)Co  {-W1SD+W2SD-XFD} {frozen core}
      call matcopy('W1','W')
      call matscal('W',two)
      call matadd1('W2',-two,'W')
      call matdef('tmw1','q',ncf,ncf)
C  construct B1 contribution to  W and Y
      call matdef('tmss02','r',ncf,nmo)
      call matmmult('virt','B1','tmss02')
      call matscal('tmss02',two)
      call matmmul2('tmss02','occa','tmw1','n','t','n')
C
C  CvB1 saved in tmss02 for Y-terms
C
C  this is the last contribution to W
      if(iprint.ge.6) call matprint('W',6)
C     call matpose('tmw1')
      call matadd1('tmw1',-two,'W')
cc
      if(iprint.ge.3) then
         call matprint('W',6)
         call matprint('X',6)
         call matprint('B1',6)
      endif
C
C  this should complete all contributions to W
C
C  now add Y-contributions
C  note that tmss02 has already been multiplied by 2 above
      call matscal('tmss02',-two)
C  this makes 4, but 2 according to my notes?
      call matmmult('ovla','tmss02','Y')
C  this  is the first contribution to Y :  Y1 = -2(4)SCv(B1)
      if(iprint.ge.6) call matprint('Y',6)
      call matrem('tmss02')
      call matrem('tmw1')
C
C  calculate the first term on eq. 33 and add to Y
      call matdef('CA','r',ncf,nval)
      call matdef('tmxni','r',ncf,nval)
      call matmmult('occu','Aik','CA')
      call matmmult('fock','CA','tmxni')
      call matscal('tmxni',-four)
      call matadd('tmxni','Yp')
      call matrem('tmxni')
c     call matprint('Y',6)
C  calculate the "density matrix' DDT=2[X-CAC] to be used for
C  construction of matrix G(DDT)
      call matmmul2('CA','occu','DDT','n','t','n')
      call matadd1('X',-one,'DDT')
      call matscal('DDT',-two)
      call matrem('CA')
C     call matprint('Y',6)
c
      return
      end
C==============BinSoRev=============================================
      subroutine BinSoRev(ncf,    nval,   iprint, ndisk1, ndisk3,
     1                    thresh, nvir,   natoms, nindxp, lbin,
     2                    nrcpb,  ncore,  gradv)
C
C   reversed binsort used for MP2 gradients.  This is reversed compared to the
C   binsort carried out in the MP2 energy calculation.
C   have Tij(ia,ib) in virtual basis need Tmylam(i,j)
C
C   This is in essence a two phase Yoshimine bin sort, but simpler.
C   Since ALL elements are present, the final Tmylam in (occupied) MO
C   basis is placed in its correct position during Phase 1.  Due to the
C   large number of bins each Tmylam may require several records.
C   the same number of integrals as the normal binsort performed in the
C   MP2 energy calculation, but many more (but smaller) bins in this
C   case
C
C   Arguments:
C   ncf       number of basis function
C   nval      number of correlated orbitals
C   iprint    printlevel
C   ndisk1    unit for <Tij> in virtual MO basis
C   ndisk3    unit number for <bins>
C   thresh    integral threshold normally 1.0d-09
C   nvir      number of virtual orbitals
C   natoms    number of atoms
C   nindxp    number of indices
C   lbin      length of one bin in words
C   nrcpb     number of record for each Tmylam
C   ncore     number of core orbitals
C   gradv     gradient vector (3*natoms)
C
C         Svein Saebo, Starkville, MS and Fayetteville, AR May 2002
C

      use memory

      implicit real*8(a-h,o-z)
c     common /intbl/maxsh,inx(100)
c     common /big/bl(30000)
      dimension gradv(3,natoms)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0,sixty=60.0d0)
C
      call secund(tt0)
      call elapsec(elaps0)
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
      ictr=igetival('ictr')
      iout=igetival('iout')
C
      ncfsq=ncf*ncf
C
C  reserve memory for one Tij
c     call getmem(ncfsq/2+1,ibuf)
      call getint_4(ncfsq,ibuf)
c     call getmem(ncfsq/8+1,i1)
      call getint_1(ncfsq,i1)
C  reserve memory for bins
c     call getmem(nindxp*lbin,ibins)
      call getint_4(nindxp*lbin*2,ibins)
c     call getmem(nindxp*lbin/4+1,i1bin)
      call getint_1(nindxp*lbin*2,i1bin)
C
      call Rphas1(ncf,    nval,   ncore,  nindxp, ndisk1,
     1            ndisk3, lbin,   nbins,bl(ibins),bl(i1bin),
     2           bl(ibuf),bl(i1), bl(ibuf),bl(i1),iprint,
     3           bl(ifp), nsym, bl(ictr),nrcpb, nvir,
     4            thresh)
      call retmem(4)
C
C  Phase 1 finished integrals now on  XXXX.bins, XXXX.bins.01  etc
C
C  Start Phase 2:
C
c  Read back bins, extract 1 amplitude matrix and backtransform
c  to AO  basis
C  will be contracted with integralderivatives on the fly and discarded
C  after it has been used.
C
C  SS Aug 2003
C  Phase2 has been moved to after CPHF to avoid repeated calculation of
C  integral derivatives.
C  The sorted integrals now in bins are used 2 places:
C  for the D1-terms and the A2 terms.
C  the amplitudes in (virtual) MO basis are not need anymore and have
C  been deleted
C
      end
C=========Rphas1===============================================
      subroutine  Rphas1(ncf,    nval,   ncore,  nindxp, ndisk1,
     1                   ndisk3, lbin,   nbin,   ibins,  i1bin,
     2                   ibuf,   int1,   jbuf,   jnt1,   iprint,
     3                  ifunpair,nsym,   inx,    nrcpb,  nvir,
     4                   thresh)

      use memory

      implicit real*8(a-h,o-z)
C
C  Bin-Sort Phase one. This is a 'reversed' binsort used for the
C  MP2 gradients. Based on Subroutine BinSort1 (in CMP2) written
C  by Peter Pulay, and Phase1 in the LMP2 program.
C
C  This subroutine reads the files containing the Residua Tij in
C  virtual basis, as one matrix Tij(ia,ib).  These are on files starting
C  with ndisk and stored as integer*4 in Canonical order without
C  indices approximately 2 GB on each file.
C
C  This matrix is first transformed to AO basis, and the one
C  matrix is sorted at the time
C
C  The memory is divided into ncf*(ncf+1)/2 bins and a given integral
C  Tij(my,lam) is placed in a bin mylam as follows:
C  Since the original Tij are read in Canonical order i.ge.j and
C  ij=i*(i-1)/2+j each bin has 2*npairs position and
C  integral T(my,lam) is put into the (1,ij) position of bin mylam
C  integral T(lam,my) is put into the (2,ij) position of bin mylam
C  the bins can normally only hold -lbin- ij values
C  when the lbin have been sorted, all bins are are full and they
C  are written to direct access file as follows:
C  nrcpb is the number of records need for each Tmylam
C  (~napirs/lbin)
C    record 1..nrcpb              T11
C    record nrcpb+1 ,   2*nrecpb  T21
C    etc.
C  Note when the nrcpb record are read back for a given my,lam
C    Tmylam(i,j) = ibins(1,ij)
C    Tmylam(j,i) = ibins(2,ij)
C  The Tmylam are used several places through out the program
C
C
C  Arguments:
C  ncf        number of basis functions
C  nval       number of valence orbitals
C  ncore      number of core orbitals
C  nindxp     number of index pairs: ncf*(ncf+1)/2
C  ndisk1     unit for <Tij> in virtual MO basis
C  ndisk3     unit number for <bins>
C  lbin       length of bin (each bin contains 2*lbin integers)
C  nbin       total number of bins
C  ibins      integer array for bins ibins(2,lbin,nindxp)
C  i1bin      integer*1 storage for precision overflow in bins
C ** WARNING: ibuf & jbuf and i1 & j1 can share storage **
C  ibuf       integer array nvir*nvir long
C  int1       integer*1 storage for precision overflow in ibuf
C  jbuf       integer array ncf*ncf long
C  jnt1       integer*1 storage for precision overflow in jbuf
C  iprint     printlevel
C  ifunpair   symmetry relations between basis functions
C  nsym       number of symmetry operations
C  inx        basis and contraction information
C  nrcpb      number of records per bin (or per mylam)
C  nvir       number of virtual orbitals
C  thresh     integral threshold, normally 1.0d-09
C
C  Svein Saebo Fayetteville Ar, and Starkville, MS May 2002
C
c     common /big/bl(30000)
      integer*1 i1bin(2,lbin,*),int1(nvir,nvir),jnt1(ncf,ncf)
      integer*4 ibins(2,lbin,*),ibuf(nvir,nvir),jbuf(ncf,ncf)
      dimension ifunpair(7,*)
      dimension inx(12,*)
      parameter(zero=0.0d0,four=4.0d0)
C
C
      tbktr=zero
      ebktr=zero
      ncfsq=ncf*ncf
      npairs=nval*(nval+1)/2
      nbin=0
      iwrb=0
      call matdef('tao','q',ncf,ncf)
      call matdef('tmo','q',nvir,nvir)
      call matdef('tvir','r',nvir,ncf)
      call matpose2('virt','tvir','n')
      imao=mataddr('tao')
      immo=mataddr('tmo')
C
C  loop over valence pairs
      ij=0
      iijj=0
c
      do ii=1,nval
         do jj=1,ii
            ij=ij+1
            iijj=iijj+1
            if (ij.gt.lbin) then
               krec=-nrcpb+1
               do irec=1,nindxp
                  krec=krec+nrcpb
                  mrec=krec+iwrb
                  nbin=nbin+1
                  call WriteBin2(ndisk3,mrec,2*lbin,ibins(1,1,irec),
     $                           i1bin(1,1,irec))
               enddo
               iwrb=iwrb+1
               ij=1
            end if
C
C  read one record and process it
C
C  first read Tij in virtual basis and transform to AO basis
C
            read(ndisk1,rec=iijj) int1,ibuf
            call secund(ttx)
            call elapsec(etx)
C NZG
C            write(6,*) "ii jj",ii,jj
            call transAO(ncf,nvir,bl(imao),bl(immo),ibuf,
     1                   int1,jbuf,jnt1,thresh)
            call secund(tty)
            call elapsec(ety)
            tbktr=tbktr+tty-ttx
            ebktr=ebktr+ety-etx
C  Tij in AO basis returned in jbuf as integer*4 words
            ncs=igetival('ncs')
            do MYS=1,ncs
               my1=inx(11,MYS)+1
               my2=inx(10,MYS)
               do LAS=1,MYS
                  lam1=inx(11,LAS)+1
                  lam2=inx(10,LAS)
                  do my=my1,my2
                     do lam=lam1,lam2
                        if(lam.gt.my) exit
                        mylam=my*(my-1)/2+lam
C
                        ibins(1,ij,mylam)=jbuf(my,lam)
                        ibins(2,ij,mylam)=jbuf(lam,my)
                        i1bin(1,ij,mylam)=jnt1(my,lam)
                        i1bin(2,ij,mylam)=jnt1(lam,my)
C
                     end do      !  lam
                  end do      !  my
               end do
            end do
         end do      !  jj loop
      end do      !  ii loop
c  write out the content of the bins at the end
      krec=-nrcpb+1
      do irec=1,nindxp
         krec=krec+nrcpb
         mrec=krec+iwrb
         nbin=nbin+1
         call WriteBin2(ndisk3,mrec,2*lbin,ibins(1,1,irec),
     $                  i1bin(1,1,irec))
      enddo
cc
      if(iprint.ge.2) then
         write(6,*)' Virtual to AO transformation: '
         write(6,100) tbktr/60.0d0,ebktr/60.0d0
  100    format(1x,1f8.2,' minutes ',1f8.2,' minutes ')
         write(6,*)' Number of bins written: ', nbin,' length ',10*lbin
      endif
cc
      call matrem('tvir')
      call matrem('tmo')
      call matrem('tao')
      end
C========Rphas2=======================================================
      subroutine Rphas2(ncf,    nval,   ndisk3, nbins,  lbin,
     1                  thresh, ibin,   int1,   Tmnmo,  Tmnao,
     2                  nsym,   iprint, inx,    nrcpb, ifunpair,
     3                  gradv,  iscs,   natom)

      use memory

      implicit real*8(a-h,o-z)
C
C      Reversed binsort phase 2. Bins are read and one Tmylam in MO
C      basis (i,j) is extracted.  This matrix is back-transformed to AO
C      basis to form Tmylam(ny,sig).
C      Integrals in bins are scaled by the inverse of the integral
C      threshold and stored as 5-byte integers
C      Note the bin file can be read sequentially, this was made
C      possible by the direct access write in phase 1
C      A scattered write is more efficient than a scattered read
C      since read are buffered
C
C     arguments:
C     ncf         number of basis functions
C     nval        number of correlated orbitals
C     ndisk3      unit number for bins
C     nbins       total number of bins on files
C     lbin        number of ij values per bin
C     thresh      integral threshold, normally 1.0d-09
C     ibin        array for storing one record
C     Tmnmo       Tmylam(i,j) real*8 array nval*nval
C     Tmnao       Tmylam(ny,sig) real*8 array ncf*ncf
C     nsym        number of symmetry operations
C     iprint      print level
C     inx         basis and contraction information
C     nrcpb       number of records per Tmylam(i,j) matrix
C     gigabyte    2.135d9
C     nrcpf       number of records per file (bins)
C     ifunpair    symmetry relations between basis functions
C     itmn        integer matrix nval*nval for temprorary storage
C     natom       number of real atoms
C     gradv       gradient vector 3*natoms
C
C     Svein Saebo, Fayetteville, AR May 2002
C
c     common /big/bl(30000)
      common /lindvec/ lind,idensp
      integer*1 int1(2,lbin)
      integer*4 ibin(2,lbin)
      character*11 scftype
      dimension Tmnao(ncf,ncf),Tmnmo(nval,*)
      dimension ifunpair(7,*)
      dimension gradv(3,*)
      dimension inx(12,*)
      dimension xintxx(9)
      parameter(zero=0.0d0,half=0.5d0,two=2.0d0,four=4.0d0)
      parameter (dblmax=2 147 483 648.0d0)
C
      call secund(trphas2_b)
      call elapsec(erphas2_b)
C
C  initialize two-el programs
C
      ancf=zero
      ncs=igetival('ncs')
      nfock=1
      thint=thresh
      iforwhat=5
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
      call getival('nocc',nmo)
      dblcmp = dblmax*thresh
      call mmark
c
      call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
     *              ax,     nrad,   nang,   lrad,   lang,
     *              Iradq,  NBatch, .true., nmo,    0,
     *              scftype,xintxx, nblocks,maxbuffer,maxlabels)
C
cc      call symmoff
C  reserve memory for labels ..
      call getmem(maxlabels,labels)
      call getival('schwarz',ischwarz)
C
C  allocate memory for an array mapping contr.func. ti contr.shells :
      ictr=igetival('ictr')
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,inx)
c
      nbin=0
      ncfsq=ncf*ncf
      tot=0
      numt=0
      tds=zero
      ttint=zero   ! cpu     time for itegral derivatives * T
      etint=zero   ! elapsed time for itegral derivatives * T
      iocca=mataddr('occu')
C
      call matdef('denf','q',ncf,ncf)
      call matdef('ddtf','q',ncf,ncf)
      call matcopy('den0','denf')
      call matcopy('DDT','ddtf')
      idena=mataddr('denf')
      iddt=mataddr('ddtf')
C
C  allocate memory for screening density
      call matdef('DS','q',ncs,ncs)
      idsad=mataddr('DS')
C
      call init_info('fock')
c----------------------------------------------------------------------
C  TRIAGULAR LOOP OVER SHELLS
c----------------------------------------------------------------------
C
      ttbt=0.0d0
      ttbad=0.0d0
      icta=mataddr('tocc')
      call getint(ncf,ipoint)
      call getint(ncf,jpoint)
c
      do MYS=1,ncs
         call get_shell_size(inx,MYS,MYS_size)
         do LAS=1,MYS
cc      if(LastSymPair(ics,kcs,nsym,ifunpair,inegl,iret)) cycle
            call get_shell_size(inx,LAS,LAS_size)
            lam1=inx(11,LAS)+1
            lam2=inx(10,LAS)
            my1=inx(11,MYS)+1
            my2=inx(10,MYS)
c
            mylsiz=LAS_size*MYS_size*ncfsq
C     if(MYS.eq.LAS) mylsiz=MYS_size*(MYS_size-1)/2
C  reserve memory for Tmylam(ny,sig) for shell pairs MYS,LAS
C  for 2 different D6 and ncf=1000, 6*6*1000*1000 36MW or 288MB
C  should be ok.
            call getmem(mylsiz,iTadr)
            call zeroit(bl(itadr),mylsiz)
C
C  inaddr is a pointer within this large array to a particular T
C  loop over functions within the two shells contruct the T's and store
C  get pointers ipoint and number of functions retained ipr
C
            call TranSetup(bl(icta),ncf,nval,ipr,jpr,
     1                     bl(ipoint),bl(jpoint),bl(ischwarz),ncs,MYS,
     2                     LAS,ancf,bl(mapf2s),thresh)
C
            my3=0
            do my=my1,my2
               my3=my3+1
               lam3=0
               do lam=lam1,lam2
                  lam3=lam3+1
                  if(lam.gt.my) exit
C     istar=istar+nrcpb
                  mylam=my*(my-1)/2+lam
                  istar=(mylam-1)*nrcpb
                  ij=0
                  iwrd=1
                  mrec=istar+1
                  read(ndisk3,rec=mrec) int1,ibin
                  nbin=nbin+1
                  do ii=1,nval
                     do jj=1,ii
                        ij=ij+1
                        if(ij.gt.lbin) then
                           iwrd=iwrd+1
                           mrec=istar+ iwrd
                           nbin=nbin+1
                           read(ndisk3,rec=mrec) int1,ibin
                           ij=1
                        endif
c -- decompress the integral ------------------------
                        If(int1(1,ij).eq.0) Then
                           xx = ibin(1,ij)*thresh
                        Else If(int1(1,ij).gt.0) Then
                           xx = ibin(1,ij)*thresh
                           xx = xx + SIGN(int1(1,ij)*dblcmp,xx)
                        Else
                           xx = ibin(1,ij)*thresh*10.0d0**(-int1(1,ij))
                        EndIf
                        Tmnmo(ii,jj)=xx
                        If(int1(2,ij).eq.0) Then
                           xx = ibin(2,ij)*thresh
                        Else If(int1(2,ij).gt.0) Then
                           xx = ibin(2,ij)*thresh
                           xx = xx + SIGN(int1(2,ij)*dblcmp,xx)
                        Else
                           xx = ibin(2,ij)*thresh*10.0d0**(-int1(2,ij))
                        EndIf
                        Tmnmo(jj,ii)=xx
c ---------------------------------------------------
                     enddo  ! over jj
                  enddo  ! over ii
                  if(my.eq.lam) call matscal('TMO',half)
                  call atoat2(Tmnmo,nval,'y',iscs)
C  backtransform this matrix to AO basis
                  call secund(ttbt1)
                  call BackTrans(Tmnao,ncf,nval,ipr,jpr,
     2                           bl(ipoint),bl(jpoint),bl(iocca))
                  call secund(ttbt2)
                  ttbt=ttbt+ttbt2-ttbt1
CC
CC   SS July 2003
CC   add extra terms to TAO here to avoid repeated integral derivatives
C
C NZG
C                  call matzero('TAO')
C NZG
                  call addtoT1(tmnao,bl(idena),bl(iddt),ncf,my,lam)
                  call moveTsh(Tmnao,bl(iTadr),ncfsq,my3,lam3,
     1                         MYS_size,LAS_size)
C NZG
C                  call zeroit(bl(itadr),mylsize)
C NZF
                  call secund(ttbt3)
                  ttbad=ttbad+ttbt3-ttbt2
                  numt=numt+1
                  tot=tot+ncfsq
               enddo  !over functons in shell LAS
            enddo  !over functons in shell MYS
C  construct screening density, result in 'DS' mataddr=idsad
            call secund(tmkds)
            call getmem(ncfsq,idsax)
            ndens=MYS_size*LAS_size
C  bl(idsad)is final  screening matrix in shells (ncfxncs)
            call mkdscl(bl(idsax),bl(idsax),bl(idsad),bl(iTadr),ncfsq,
     1                  ndens,inx,ncf,ncs,dsmax)
            call setrval('dmx_mp2d',dsmax)
            call retmem(1)
            call secund(tmkds2)
            tds= tds+tmkds2-tmkds
C
C  all TMYLAM for the pair of shells are stored in bl(iTadr)
C
C  Calculate integral derivatives
C_____________________________________________________________________
ckw..........
            call secund(txxx1)
            call elapsec(exxx1)
cccc  call init_info('fock')
ckw..........
C_____________________________________________________________________
            call mmark
            call int_mp2d(bl,inx,thresh,MYS,LAS,
     *                    bl(mapf2s),bl(idsad),bl(itadr),nintotal,gradv)
            call retmark
C-----------------------------------------------------------------
            call secund(txxx2)
            call elapsec(exxx2)
cccc  txxx0=txxx1
ccccc call term_info(thresh,txxx2-txxx1,txxx1-txxx0,'fock')
C-----------------------------------------------------------------
            ttint=ttint+txxx2-txxx1
            etint=etint+exxx2-exxx1
C  density for Screening : bl(idsad)
C  densities to be used  : bl(iTadr)
C  this matrix can be considered a 3 dimensional matrix:
C  T(ncfsq,MYS_size,LAS_size)
C  Current fixed SHells: MYS, LAS  (used to be ics,kcs)
C  similar to MP2, request integral derivaltives
C  (MYS,ny,LAS,si)x  for fixed shells MYS and LAS and all ny,si
C-----------------------------------------------------------------
            call retmem(1) !mylsiz
         enddo  ! shells LAS
      enddo  ! shell MYS
      call retmem(2) !ipoint,jpoint
c----------------------------------------------------------------------
C  END of TRIAGULAR LOOP OVER SHELLS
C-----------------------------------------------------------------
c in gradv we have a very first contributions to the mp2 grad.
c There are contr. from derivative integrals and backtransf
c amplitudes. These contrib. HAVE TO BE rescaled by the integ. thresh.
c
      call dscal(3*natom,4.0d0*thresh,gradv,1)
      call matrem('DS')
      call matrem('ddtf')
      call matrem('denf')
      call retmem(3)
      call retmark
ckw
c------------------------------------------------------------
      write(6,100) ttint/60.0d0,etint/60.0d0
  100 format(/'CPU and elapsed time for deriv.integ. & mult. with T :',
     *       1x,f8.2,' min ',f8.2,' min')
c------------------------------------------------------------
      call secund(trphas2_e)
      call elapsec(erphas2_e)
      trphas2=(trphas2_e-trphas2_b)/60.0d0
      erphas2=(erphas2_e-erphas2_b)/60.0d0
      write(6,110) trphas2 , erphas2
  110 format(/'Total CPU and elapsed time in Rphas2 :',
     *       1x,f8.2,' min ',f8.2,' min')
c------------------------------------------------------------
ckw
cc
      if(iprint.ge.2) then
        ancf=ancf/(ncf*ncs*(ncs+1))
        ancf=2.0d2*ancf
        write(6,*) ' Percent of AOs used for last transformation: ',ancf
        write(6,*) ' time adding extra terms to T: ',ttbad/60.0d0
        write(6,*) ' time second backtransformaton: ',ttbt/60.d0
        write(6,*) ' time constructing DS: ', tds/60.0d0
        write(6,*) numt,' TAO matrices of dimension ',ncfsq,' generated'
        write(6,*) ' total number of TAO(my,lam)s: ',tot,' elements'
      endif
      end
C===========D1terms=======================================
      subroutine D1terms(ndisk3, Y,      nval,   iprint, thresh,
     1                   ncore,  ibins,  int1,   lbin,   nrcpb,
     2                   ncf,    nmo,    iscs)

      use memory

      implicit real*8(a-h,o-z)
C  calculate the D1 contribution to matrix Y, Matrix Y3 Eq. 46
C  this is done by calculating the quarter transformed integrals and
C  contracting them with the sorted amplitude matrix Tmylam(i,j). These
C  were generated for the reversed binsort for the A2 contributions and
C  are still stored in the bins.
C
C  Svein Saebo Fayetteville, AR Summer 2002
C
C  Called once from mp2_grad
C
C  External Calls
C  mp2 integral routines +
C  mkdsqt       constructs screening matrix appropriate for the quarter
C               transformed inegrals.
C  trasigtoj    transforms one index (sigma) from AO-basis to Mo=O basis
C               j.  Same block structure as for mp2 energy calculation
C               When one block (my,ny|lam,j) is calculated the integral is
C               multipled with the (i.j) element of the amplitude Tmylam
C               and accumulated in matrix Y3 according to equation 46
C               see routine trasigtoj and routined called within for
C               further details
C   ARGUMENTS
C   ndisk3   unit for bins containing sorted amplitudes
C   Y        Array(ncf,nval)  Matrix Y Eq 59 or 61, contribution 4*Y3
C            is calculated here
C   nval     number of correlated orbitals
C   iprint   print parameter
C   thresh   integral threshold
C   ncore    number of core-orbitals
C   ibins    integer array(2,lbin) for sorted amplitudes
C   int1     integer*1 storage for precision overflow
C  ****finish later..
c-----------------------------------------------------------------------
ckw
      dimension Ipass(2,28), Kpass(2,28) !28 is 28 comp.of cartisian I-f
c-----------------------------------------------------------------------
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      character*11 scftype
      logical LastSymPair,smal,dualbasis
      integer*1 int1(2,*)
      integer*4 ibins(2,*)
      dimension Y(ncf,*)
      dimension xintxx(9)
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,four=4.0d0,sixty=60.0d0)
C
      call elapsec(elaps0)
      call secund(tt0)
      ncs=igetival('ncs')
      ictr=igetival('ictr')
C
      dualbasis=.false.
C
c  get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c  basis function symmetry pairs are in inx(ifp)
c .............................................................
      natoms=igetival('na')
      if(natoms.le.30) then
         smal=.true.
      else
         smal=.false.
      end if
      smal=.true.       ! always use small option
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
      nfock=1
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
c allocate memory for an array mapping contr.func. to contr.shells :
c
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c
c-------------------------------------------------------
c get memory for MOs, orbital energies and screening density
c
      call matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
      ioccu=mataddr('occu')
      idics=mataddr('dsmx')
C     call matprint('occu',6)
c-------------------------------------------------------
      call getmem(maxlabels,labels)
c---------------------------------------------------------------------
c Construct "screening density" matrix to be used to estimate the
c importance of integrals.
      call getmem(ncf*ncf,idsa)
      call getmem(ncf,ivec)
      call mkdsqt(bl(idsa),bl(ivec),bl(idics),nval,ncf,
     1            bl(ioccu),bl(ictr),ncs,densmax)
      call retmem(2)
c---------------------------------------------------------------------
c  establish orbital symmetry characters
      nfirst= ncore+1
      nlast=nmo
C  this causes problems, jump over this for now.  SS
      if(nsym.gt.1000) then
         call getint(nval*nsym,iorbpair)
         call getmem(ncf,itmp)
         icano=mataddr('cano')
         call OrbSymPair(nsym, ncf ,   bl(icano), bl(ifp) ,nfirst,
     1                   nlast,iprint, bl(itmp),  bl(iorbpair))
         call retmem(1)
      endif
c.................................................
c   reserve space for one AO integral matrix
      call matdef('xmat','q',ncf,ncf)
      ixadr=mataddr('xmat')
c  nrec is the total number of records written on all files
c  irec is the current counter in the file
      nrec=0
      irec=0
c
c---------------------------------------------------
c     if(dualbasis) then
c  update content of "basis sets" according to current needs
c  i.e. put small basis set into "basis 2" & "basis 4"
c
c       call update_basis(bl,ncs_sm)
c     endif
c---------------------------------------------------
c
      elapint=0.0d0
      call secund(tt3)
      call elapsec(telap3)
c  The next 2 counters count retained and omitted contracted shells
      inegl=0
      iret=0
c  iktot is the number of ics,kcs pairs really calculated & transformed
c   if an (ics,kcs) pair is skipped because there are no integrals in it
c   then it is NOT incremented
      iktot=0
cc      call symmoff
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
         call getint(ncf,if2cr)
         call getint(ncf,if2cc)
         call memstat(nreq,nmark,lastaddr,memtot,iceiling,ioffset)
         lcore=igetival('lcore')
         kcore=lcore-ioffset
         memused=lastaddr-ioffset
         mem=kcore-memused-6*ncf*ncf-2000000
         memx=mem/9
         mem=memx*4
         lenindj=memx*4
         call getmem(memx,indlj)
         write(6,*)' BIG option used for integrals'
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
c
c-----------------------------------------------------------------------
      ncf2=ncf*ncf
c-----------------------------------------------------------------------
      ttrans=0.0d0
      elaptrans=0.0d0
      tinteg=0.0d0
      nbin=0
      nskipped=0
      do ics=1,ncs
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf2*ics_size
         do kcs=1,ics
cc      if(LastSymPair(ics,kcs,nsym,inx(ifp1),inegl,iret)) cycle
c  of each (ics,kcs) pair, calculate only the last one
c  (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
c
c
            call get_shell_size(bl(ictr),kcs,kcs_size)
            call secund(tt2)
            call elapsec(telap2)
            ttrans=ttrans+tt2-tt3
            elaptrans=elaptrans+telap2-telap3
c
            if(smal) then
               lmp2_size=lmp2_siz1*kcs_size
c
c check if a size of the lmp2 integral buffer is not too big
c if it is  then split over kcs ( and possibly ics)
c
c
               call check_size1(lmp2_size,ncf,nval,nval,ntimes)
c
               call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                         ntimes,Ipass,Kpass,Itimes,Ktimes)
c
c
               do itime=1,itimes
                  icf1=ipass(1,itime)
                  icf2=ipass(2,itime)
                  iatonce=icf2-icf1+1
                  do ktime=1,ktimes
                     kcf1=kpass(1,ktime)
                     kcf2=kpass(2,ktime)
                     katonce=kcf2-kcf1+1
c
                     lmp2_size=iatonce*katonce*ncf2
c
                     call getmem(lmp2_size,lmp2int)
                     call mmark
                     call int_lmp2(bl,bl(ictr),thresh,ics,icf1,icf2,kcs,
     &                             kcf1,kcf2,bl(mapf2s),bl(idics),
     &                             iprint,bl(lmp2int),nintotal,nrow,
     &                             ncol,bl(irow),bl(icol),bl(lzero))
                     call retmark
                     if (nintotal.eq.0) then
                        call retmem(1)
                        cycle
                     else
                        iktot=iktot+1
                     end if
c
ccc      write(*,*) 'lmp2 int',ics,kcs,(bl(lmp2int+ii-1),ii=1,lmp2_size)
c
                     call secund(tt3)
                     call elapsec(telap3)
                     tinteg=tinteg+tt3-tt2
                     elapint=elapint+telap3-telap2
                     call mmark
                     call matdef('Tmyl','q',nval,nval)
                     itmylam=mataddr('Tmyl')
c......................................................................
                     call Trasigtoj(ncf,ncs,nval,ics,icf1,icf2,kcs,kcf1,
     1                              kcf2,bl(ictr),bl(lmp2int),bl(ioccu),
     2                              iprint,thresh,bl(ixadr),nrow,ncol,
     3                              bl(irow),bl(icol),bl(irow1),
     4                              bl(icol1),smal,Y,bl(itmylam),ibins, 
     5                              int1,lbin,nrcpb,ndisk3,iscs)
                     call matrem('Tmyl')
                     call retmark
c......................................................................
                     call retmem(1)          ! lmp2int
c......................................................................
                  enddo    ! over ktime (selected kcf belonging to kcs shell )
               enddo    ! over itime (selected icf belonging to ics shell )
c........................................................................
            else
               ind=0
               intstore=0
               call mmark
               call getmem(lmp2_size,lrestor)
               call int_lmp2b(bl,bl(ictr),thresh,ics,kcs,bl(mapf2s),
     1                        bl(idics),iprint,bl(lmp2int),nintotal,
     2                        bl(icol),bl(irow),bl(if2cc),bl(if2cr),
     3                        bl(indlj),bl(lrestor),nrow,ncol,ind,
     4                        intstore,lmp2_size,lenmax)
               call retmark
               ttint=ttint+nintotal
               if (nintotal.eq.0) then
                  cycle
               else
                  iktot=iktot+1
               end if
c
ccc      write(*,*) 'lmp2 int',ics,kcs,(bl(lmp2int+ii-1),ii=1,lmp2_size)
c
               call secund(tt3)
               call elapsec(telap3)
               tinteg=tinteg+tt3-tt2
               elapint=elapint+telap3-telap2
               call mmark
               call matdef('Tmyl','q',nval,nval)
               itmylam=mataddr('Tmyl')
C
               call getinfs(bl(ictr),ics,kcs,icf1,icf2,kcf1,kcf2)
c......................................................................
               call Trasigtoj(ncf,ncs,nval,ics,icf1,icf2,kcs,kcf1,kcf2,
     1                        bl(ictr),bl(lmp2int),bl(ioccu),iprint,
     2                        thresh,bl(ixadr),nrow,ncol,bl(irow),
     3                        bl(icol),bl(irow1),bl(icol1),smal,Y,
     4                        bl(itmylam),ibins,int1,lbin,nrcpb,ndisk3,
     5                        iscs)
               call matrem('Tmyl')
               call retmark
c......................................................................
            endif
         end do
      end do
cc
      if(iprint.ge.2) then
         write(iout,61)
         write(iout,*) '  '
c
         write(iout,62) tinteg/sixty, elapint/sixty
         write(iout,64) ttrans/sixty, elaptrans/sixty
      endif
cc
   61 format(' CPU & Elapsed time for D1 contribution')
   62 format(' AO Integrals     =',f8.2,' and ',f8.2,' minutes')
   64 format(' Xnyi matrix      =',f8.2,' and ',f8.2,' minutes')
c
      end
C==============Trasigtoj================================================
      subroutine Trasigtoj(ncf,    ncs,    nval,
     *                     ics,    icf1,   icf2,   kcs,    kcf1,   kcf2,
     1                     inx,    xint,   CMO,    iprint, thresh,
     2                     xmat,   nrow,   ncol,   irow,   icol,
     3                     irow1,  icol1,  smal,   Y,      Tmylam,
     4                     ibins,  int1,   lbin,   nrcpb,  ndisk3,
     5                     iscs)

      use memory

      implicit real*8 (a-h,o-z)
c  Arguments:
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
c  iprint: print level
c  thresh = integral threshold
c  nrow   = number of non-zero rows of the AO integral matrices
c  ncol   = ditto for columns
c  irow(k),k=1,nrow = the indices of the non-zero rows
c  icol(k),k=1,nco  = ditto for columns
c  STORAGE:
c  xmat: (ncf,ncf) place for 1/4 transformed integrals
c  intent=OUT:
c  halftr: space for half-transformed integrals to be written on disk
c  integrals are {mu,nu|lam,isig);  mu.ne.lam
c
      integer*1 int1(2,lbin)
      integer*4 ibins(2,lbin)
      dimension xint(*),xmat(ncf,ncf)
      dimension CMO(ncf,*),irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      dimension Y(ncf,*),Tmylam(nval,*)
      dimension inx(12,*)
c     common /big/bl(30000)
      logical smal
      parameter(half=0.5d0,one=1.0d0,two=2.0d0,four=4.0d0)
      parameter (dblmax=2 147 483 648.0d0)
c
c  return if there are no integrals
      if(nrow.eq.0.or.ncol.eq.0) RETURN
c
       dblcmp = dblmax*thresh
c
c  icol and jcol store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
c  itrunc serves as temporary storage for a compacted coefficient matrix
c  only the rows (or columns in the second quarter transformation) of
c  X which are non-zero are  present in bl(itrunc)
c
      call getmem(ncf*nval,itrunc)
C
ckw   kcf1=inx(11,kcs)+1
ckw   kcf2=inx(10,kcs)
ckw   icf1=inx(11,ics)+1
ckw   icf2=inx(10,ics)
c
      do my=icf1,icf2
         do lam=kcf1,kcf2
            if(lam.gt.my) exit
            if(smal) then
               call ExtractOne(ncf,   my,    lam,   icf1,  icf2,
     2                         kcf1,  kcf2,  xint,  xmat,  nrow,
     3                         ncol,  irow,  icol,  nrow1, ncol1,
     4                         irow1, icol1)
            else
               call mmark
               call ExtractEn(ncf,   my,    lam,  icf1,  icf2,
     2                        kcf1,  kcf2,  xint, xmat,  nrow,
     3                        ncol,  irow,  icol, nrow1, ncol1,
     4                        irow1, icol1)
               call retmark
            endif
c
            if(nrow1.eq.0.or.ncol1.eq.0) cycle
            call matdef('quartra','r',nval,ncol1)
            call matdef('quartr2','r',nrow1,nval)
            iquatr=mataddr('quartra')
            iquat2=mataddr('quartr2')
            call Transbatch(ncf,nval,nrow1,ncol1,irow1,
     1                      icol1,xmat,CMO,itrunc)
            call matscal('quartra',thresh)
            call matscal('quartr2',thresh)
C
C  batch of quartertransformed integrals for current my,lam
C  is returned in quartra (from right) and in quartr2 (from left)
C
C  get Tmylam(i,j) from bins:
C
            mylam=my*(my-1)/2+lam
            istar=(mylam-1)*nrcpb
            mrec=istar+1
            nbin=nbin+1
            read(ndisk3,rec=mrec) int1,ibins
            iwrd=1
            ij=0
            do ii=1,nval
               do jj=1,ii
                  ij=ij+1
                  if (ij.gt.lbin) then
C  end of bin  read next
                     iwrd=iwrd+1
                     mrec=istar+iwrd
                     nbin=nbin+1
                     read(ndisk3,rec=mrec) int1,ibins
                     ij=1
                  endif
c -- decompress the integral ------------------------
                  If (int1(1,ij).eq.0) Then
                     xx = ibins(1,ij)*thresh
                  Else If(int1(1,ij).gt.0) Then
                     xx = ibins(1,ij)*thresh
                     xx = xx + SIGN(int1(1,ij)*dblcmp,xx)
                  Else
                     xx = ibins(1,ij)*thresh*10.0d0**(-int1(1,ij))
                  EndIf
                  Tmylam(ii,jj)=xx
                  If (int1(2,ij).eq.0) Then
                     xx = ibins(2,ij)*thresh
                  Else If(int1(2,ij).gt.0) Then
                     xx = ibins(2,ij)*thresh
                     xx = xx + SIGN(int1(2,ij)*dblcmp,xx)
                  Else
                     xx = ibins(2,ij)*thresh*10.0d0**(-int1(2,ij))
                  EndIf
                  Tmylam(jj,ii)=xx
c ---------------------------------------------------
               enddo
            enddo
C  form T tilda
            call atoat2(Tmylam,nval,'n',iscs)
C
            call matscal('Tmyl',four)
C
C  form matrix Ynyi (eq.46)
C
            call FormXnyi(Y,Tmylam,bl(iquatr),nval,ncol1,icol1,ncf)
            if(my.gt.lam)
     1         call FormXny2(Y,Tmylam,bl(iquat2),nval,nrow1,irow1,ncf)
C
            call matrem('quartr2')
            call matrem('quartra')
         enddo
      enddo
      call retmem(1)
      end
C=====TransBatch===================================================
      subroutine transbatch(ncf,nval,nrow,ncol,irow,
     1                      icol,xmat,coef,itrunc)
C  transform (my,ny,lam,sig) to (my,ny,lam,j)
c   ARGUMENTS:
c  ncf: number of contracted basis functions
c  nval: number of orbitals to be transformed (valence orbitals)
c  mu,lam: fixed basis function indices. The transformed integrals are
c          (mu,ny|lam,j) where j is  MO index
c  iprint: print level
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
      dimension irow(nrow),icol(ncol)
      dimension xmat(ncf,ncf), coef(ncf,nval)
c     common /big/bl(30000)
c  defined in the calling program
c  xmat has the COMPACTED AO exchange matrix
c    prepare the compacted SCF coefficients for the left side
c   'ymat' is similar to 'xmat' but is compacted
      call CompactCoef(ncf,nval,coef,nrow,irow,bl(itrunc))
      call matdef('ymat','r',nrow,ncol)
      call PutTrunc(ncf,nrow,ncol,xmat,bl(mataddr('ymat')))
c  connect the truncated SCF coefficient matrix to the matrix system
      call matconn('trunc','r',nrow,nval,itrunc)
      call matmmul2('trunc','ymat','quartra','t','n','n')
      call matdisc('trunc')
      call CompactCoef(ncf,nval,coef,ncol,icol,bl(itrunc))
c  connect the truncated SCF coefficient matrix to the matrix system
      call matconn('trunc','r',ncol,nval,itrunc)
      call matmmult('ymat','trunc','quartr2')
      call matdisc('trunc')
      call matrem('ymat')
      end
C==========updatz========================================
      subroutine updatz(z,y,zt,ea,ei,
     1                  nvir,nmo,zold,diffm)
C  updates the z-matrix during cphf iterations
C
C  Svein Saebo Fayetteville, AR summer 2002
C
C  called form chfz once for each iteration
C
C     OUTPUT   z   z-matrix dimension(nvir,nval)
C     INPUT    y   y-matrix dimension(nvir,nval)
C     INPUT    zt   1/4 G(DP)  where DP =(Co+)z(Cv)+(Cv+)z+(Co)
C     INPUT    ea,ei  vectors with orbital energies
C     nvir     number of virtual orbitals
C     nval     number of occupied orbitals
C     zold     used for temporary storag eof old z-matrix
C     diffm    large absolute difference between old and new z-mat
C
      implicit real*8(a-h,o-z)
      dimension z(nvir,*),y(nvir,*),zt(nvir,*)
      dimension ea(*),ei(*),zold(*)
      do ii=1,nmo
      do ia=1,nvir
      z(ia,ii)= - (y(ia,ii)+zt(ia,ii))/(ea(ia)-ei(ii))
      enddo
      enddo
      call matscal('Zold',-1.0d0)
      call matadd('Zai','Zold')
      ix=idamax(nvir*nmo,zold,1)
      diffm=abs(zold(ix))
      end
C=======MakeGma2======================================================
      subroutine MakeGma2(ncf,nmo,nval,nvir,thresh,
     1                    bl,inx,dens,gmat)

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
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*)
      logical LastSymPair
      dimension xintxx(9)
c     common /intbl/maxsh,inxx(100)
      logical nofr,exst,restrt,dualbasis,smal
      dimension ifilestat(13)
      dimension dens(*),gmat(*)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
C
      rhf=.true.
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c  get symetry info: number of symm. op. starting address of contr. info
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
      ncs=igetival('ncs')
c  zero out irrelevant dft stuff
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
      iforwhat=1
      thint=thresh
      nfock=1
      call ptwoint (nfock, thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             0,      scftype,xintxx, nblocks,maxbuffer,
     *             maxlabels)
c
c allocate memory for an array mapping contr.func. to contr.shells :
c
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,inx)
c-------------------------------------------------------
c reserve memory for labels ..
c
      call getmem(maxlabels,labels)
      mywork=0
      ax=0.0d0
      igran=0
      nfock=1
      thres1=thresh
      mgo=1
      ntri=ncf*(ncf+1)/2
      call zeroit(gmat,ntri)
c
      call int_fmp2(nblocks,nfock,ncf,bl,inx,
     1              thres1,dens,gmat,dens,bl(labels),
     2              mywork,igran)
c
      nsym=igetival('nsym')
      if(nsym.gt.0) then
        call getival('ngener',ngener)
        call getival('SymFunPr',ifp)
        nfock=1
        call FockSymm_SCF(ngener,ncf,nfock,gmat,bl(ifp))
      endif
C       call matprint('Gmat',6)
      end
C=====FxSx=========================================================
      subroutine FxSx(natoms,iprint,nfocu,novlu)
C   calculates Fock matrix buildt from derivative integrals
C   and overlap derivatives.  The 3*natoms matrices of dimension
C   ncf*(ncf+1)/2 are written to files nfocu and novlu, repectively
C   the parts needed for MP2-gradients are taken from hessana by
C   KW.
C
C   Svein Saebo, Fayetteville, AR, Summer 2003
C
C   called from mp2_grad
C
C   Arguments
C   natoms   number of atoms
C   iprint   print parameter
C   nfocu    unit for Fock-derivatives
C   novlu    unit for overlap-derivatives

      use memory

      implicit real*8 (a-h,o-z)
      character jobname*256,datestr*24,cdum*20
      Logical rhf
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /cpu/ intsize,iacc,icache,memreal
c     common /big/ bl(10000)
c     common /intbl/ifpp,inx(100)
C     data ncall/0/
C     save ncall
c-------------------------------------------
C     For MP2-gradients:
C     calculates Sx and the one-electron part of FX and writes the
C     results to disk (nfocu,novlu)
C     taken from hessana, but removed stuff not needed
C     for mp2-gradients
C     SS, Summer 2002
C
C     now only one-electron contributions are calculated here
C     two electron contributions are now calculated with the
C     A2 terms
C     SS Summer 2003
c-------------------------------------------
      PARAMETER (IUnit=1)               ! unit number for checkpoint I/O
c-----------------------------------------------------------------------
      intsizer=intsize
      iaccr=iacc
      icacher=icache
      memrealr=memreal
C
      rhf=.true.
C       ncall=ncall+1
      call mmark
      ncall=-1
      call readbl(ncall,ictr)
c
      thref1=1.0d-8
      idft=0
      call setrval('thref1',thref1)
c
c get original values of symmetry pair pointers
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
      natoms3=3*natoms
      ntri=ncf*(ncf+1)/2
c---------------------------------------------------------------
c open files for : fock1, overlap1
C   NOTE!  units have been changed by SS for MP2 gradients
c units             39     40
c
      call openFxSx(3*ntri,iprint,nfocu,novlu)
c
      call getmem(natoms3*ntri,lfock1)
      call setival('lfock1',lfock1)
c
      call zeroit(bl(lfock1),natoms3*ntri)
c
c  find symmetry unique atoms :
c
      call find_unqat(bl,natoms)
c
c calculates at once :
c
c       second part of the derivative fock matrices
c       G(D0,gx) stored in lfock1 . This is the very
c       first contrib. there and it is rescaled right
c       away by the integ. thresh.
c       (later the 1-el Hx will be added to it forming
c       the final fock1 (derivative fock) matrix .
c---------------------------------------------------------------
c
c THE ONE-ELECTRON PART OF THE HESSIAN CALCULATIONS:
c
c---------------------------------------------------------------
c
C   calculate one-electron contributions
      call mmark
      call onefock1(rhf,bl,bl(ictr),iprint,nfocu,novlu,natoms)
      call retmark
c
      call retmark
      intsize=intsizer
      iacc=iaccr
      icache=icacher
      memreal=memrealr
      end
C========openfxsx=============================================
      subroutine openfxsx(ndim,iprint,nfile1,nfile2)
      implicit real*8 (a-h,o-z)
      character filename1*256
      character filename2*256
      character scrfile*256
C    this is a simplified version of open4hess,taken from hessana.f
C    SS, June 2002
c----------------------------------------------------
c Input
c ndim - record length
c----------------------------------------------------
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      filename1=scrfile(1:len)//'.fock1'
      filename2=scrfile(1:len)//'.over1'
c----------------------------------------------------
      lrec=ndim*8
c
      open (unit=nfile1,file=filename1,
     *      form='unformatted',access='direct',recl=lrec)
c
      open (unit=nfile2,file=filename2,
     *      form='unformatted',access='direct',recl=lrec)
c----------------------------------------------------
      if(iprint.ge.2)
     1  write(6,*)' files',nfile1,nfile2,' opened for Fx and SX'
      end
C=======onefock1=============================================
C      same as hess1 but removed stuff not needed for MP2
C      gradients.  Calculates  Sx and the one-electron part
C      to Fx only.
C      SS, Summer 2002
C=============================================================
      subroutine onefock1(rhf,bl,inx,iprint,nfocu,novlu,natom)
C   rhf  logical flag, true here
C   bl   common/big/
C   inx  array with contraction information
C   iprint  print parameter
C   nfocu   unit for fock derivatives, note NOT 61 as in hessana
C   novlu   unit for overladerivatives, NOT 62
C   natom   number of real atoms
C
C  call from FxSx once
C
C  extenal calls: several integral routines.

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c external field info :
      common /fieldx/ xfield,xfgrad,elfiel(9)
c Hessian  options :
caug  common /forcdbl/ thre1,thre2
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------------------------
      nprint=igetival('printh')
      call secund(tgrad1)
      call elapsec(elagrad1)
c------------------------------------------------------------------
      ifield=xfield
      ifgrad=xfgrad
c------------------------------------------------------------------
c memory was allocated in hessana : get addresses
c
      call getival('lfock1',lfock1)
c
c lfock1 contains already 2-el. contr.  F(D0,G(D0,g1)} (rescaled)
c
c------------------------------------------------------------------
c check-run only or lack of scf convergence :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c------------------------------------------------------------------
c -- Before calling first derivative integral routines need to mimic
c -- the listreal array - this was added by KW to improve symmetry
c -- handling in the hessian routines, but was overlooked here
c
      call getint(natom,listreal)
      call setup_listreal(natom,bl(listreal))
c
c Calculate T1 (kinetic) first derivative integrals
c for derivative fock matrix
c
      call intonFder(1,natom,inx,ifield,ifgrad,elfiel,
     *            bl(ibas),bl(inuc),ncs,ncf,ntri,
     *            bl(listreal),bl(lfock1) )
c
      if(nprint.ge.6) then
         write(6,*)' Derivative Fock after 1st-order Kinetic '
         call fder_print(bl(lfock1),natom,ntri,1.d0  )
      endif
c------------------------------------------------------------------
c Calculate first-order electron-nuclear attraction contributions
c to the derivative Fock matrix
c
      call elenuc_Fder(natom,1,natom,inx,bl(ibas),bl(inuc),ncs,ncf,ntri,
     *                 bl(listreal),bl(lfock1) )
c
      if(nprint.ge.6) then
         write(6,*)' Derivative Fock after 1st-order Ele-Nuc  '
         call fder_print(bl(lfock1),natom,ntri,1.d0  )
      endif
c------------------------------------------------------------------
c If necessary, calculate first-order pseudopotential contributions
c to the derivative Fock matrix
c
      call getival('npsp',npsp)            ! pseudopotentials
      If(npsp.GT.0) Then
        lfock1B=lfock1                     ! currently closed-shell only
        call pspFder(.true.,natom,natb,nate,npsp,ncf,ncs,inx,bl(inuc),
     $               bl(ibas),bl(listreal),bl(lfock1),bl(lfock1B))
c
        if(nprint.gt.6) then
          write(6,*)' Derivative Fock after pseudopotentials'
          call fder_print(bl(lfock1),natom,ntri,1.d0  )
        endif
      EndIf
c------------------------------------------------------------------
c
c -- First-order derivative Fock matrices are done.
c -- Need to store them on disk in transposed form on unit nfocu
c -- transpose fock1(3,natoms,ntri) into fock1t(ntri,3,natoms)
c
      call getmem(3*ntri,itemp)
      call trsp_disk(nfocu,natom,ntri,bl(itemp),bl(lfock1))
      call retmem(1)
c------------------------------------------------------------------
c Calculate S1 (overlap) first derivative integrals
c Use lfock1 allocation :
c
      call zeroit(bl(lfock1),natom*3*ntri)
      call intonSder(1,natom,inx,bl(ibas),bl(inuc),ncs,ncf,ntri,
     *               bl(listreal),bl(lfock1) )
c
      if(nprint.ge.6) then
         write(6,*)' 1st-order Overlap          '
         call fder_print(bl(lfock1),natom,ntri,1.d0  )
      endif
c
c Transpose over1(3,natoms,ntri) into over1t(ntri, 3,natoms)
c then write them to disk
c
      call getmem(3*ntri,itemp)
      call trsp_disk(novlu,natom,ntri,bl(itemp),bl(lfock1))
      call retmem(1)
c------------------------------------------------------------------
      if(iprint.ge.2) then
      call secund(tgrad2)
      call elapsec(elagrad2)
      total=tgrad2-tgrad1
      elaps=elagrad2-elagrad1
      write(iout ,400) total/60.0d0,elaps/60.0d0
  400 format('Master CPU time for Sx and 1e-part of FX='
     *,f8.2,' Elapsed = ',f8.2,' min')
      endif
      end
C=====Makegrad=============================================
      subroutine MakeGrad(natoms,gradv,FxSx,iunfs,ntri,
     1                    WX,ncf)
C
C  contructs gradient contribution trace(FxSx * WX)
C  the contribution is added to gradv
C  Svein Sabeo, Fayetteville, AR Summer 2002
C  Called from Mp2_grad (twice) once for X- once for W-terms
C  natoms     number of atoms
C  gradv      gradient vector (3,natoms)
C  FxSx       contains either F(x) or S(x), dimension (ntri,3)
C  iunfs      unit for F(x) or S(x)
C  ntri       ncf*(ncf+1)/2
C  WX         ncf*ncf matrix, contains either matrix W or X
C             W is contracted with S(x) X is contracted with F(x)
C  ncf        number of basis functions
C

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30000)
      dimension gradv(3,natoms),FXSX(ntri,3),WX(*)
      call matdef('fxq','q',ncf,ncf)
C fxq matrix for quadratic storage of one derivative matrix for correct
C trace.
      ifxq=mataddr('fxq')
      nqd=ncf*ncf
      do iat=1,natoms
         read(unit=iunfs,rec=iat)FXSX
C     do ixyz=1,3
         call matcopy('fxsx','fxq')
         gradv(1,iat)=gradv(1,iat) -
     1                ddot(nqd,bl(ifxq),1,WX,1)
         call matcopy('fxsy','fxq')
         gradv(2,iat)=gradv(2,iat) -
     1                ddot(nqd,bl(ifxq),1,WX,1)
         call matcopy('fxsz','fxq')
         gradv(3,iat)=gradv(3,iat) -
     1                ddot(nqd,bl(ifxq),1,WX,1)
C     enddo
      enddo
      call matrem('fxq')
      end
C=====Makegra2=============================================
      subroutine MakeGra2(natoms,gradv,FxSx,iunfs,ntri,
     1                    WX,nval,xymo)
C
C  contructs gradient contribution trace(F(x)(MO) * A)
C  the contribution is added to gradv
C  Svein Sabeo, Fayetteville, AR Summer 2002
C  Called from Mp2_grad once for A-terms
C  natoms     number of atoms
C  gradv      gradient vector (3,natoms)
C  FxSx       contains either F(x) , dimension (ntri,3)
C  iunfs      unit for F(x)
C  ntri       ncf*(ncf+1)/2
C  WX         ncf*ncf matrix, contains  matrix A
C  nval       number of valence orbitals
C  xymo       storage for one F(x) in MO basis
C
      implicit real*8(a-h,o-z)
      dimension gradv(3,natoms),FXSX(ntri,3),WX(*),xymo(*)
      nstri=nval*nval
      do iat=1,natoms
      read(unit=iunfs,rec=iat)FXSX
C     do ixyz=1,3
C   need to roll out loop over xyz here
      call matsimtr('fxsx','occu','xyzmo')
      gradv(1,iat)=gradv(1,iat)+ddot(nstri,WX,1,xymo,1)
      call matsimtr('fxsy','occu','xyzmo')
      gradv(2,iat)=gradv(2,iat)+ddot(nstri,WX,1,xymo,1)
      call matsimtr('fxsz','occu','xyzmo')
      gradv(3,iat)=gradv(3,iat)+ddot(nstri,WX,1,xymo,1)
      enddo
      end
C=====transAO==================================================
      subroutine transAO(ncf,nvir,tao,tmo,ibuf,
     1                   int1,jbuf,jnt1,thresh)
C    transforms one Tij originally in virtual MO basis to AO basis
C
C    Svein Saebo, Fayetteville, AR May 2002
C
C    Arguments:
C    ncf   - number of basis functions
C    nvir  - number of virtuals
C    tmo   - array with Tij in MO basis
C    tao   - Tij in AO basis
C    ibuf  - Tij in MO basis as integer*4
C    int1  - integer*1 storage for precision overflow
C    jbuf  - Tij in AO basis as integer*4
C    jnt1  - integer*1 storage for precision overflow
C
      implicit real*8(a-h,o-z)
      integer*1 int1(nvir,*),jnt1(ncf,*)
      integer*4 ibuf(nvir,*),jbuf(ncf,*)
      dimension tmo(nvir,*),tao(ncf,*)
      parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
C
      dblcmp = dblmax*thresh
c
      do ib=1,nvir
         do ia=1,nvir
c -- decompress the integral ------------------------
            If (int1(ia,ib).eq.0) Then
               xx = ibuf(ia,ib)*thresh
            Else If(int1(ia,ib).gt.0) Then
               xx = ibuf(ia,ib)*thresh
               xx = xx + SIGN(int1(ia,ib)*dblcmp,xx)
            Else
               xx = ibuf(ia,ib)*thresh*10.0d0**(-int1(ia,ib))
            EndIf
c ---------------------------------------------------
            tmo(ia,ib)=xx
         enddo
      enddo
c
      call matsimtr('tmo','tvir','tao')
      call matscal('tao',one/thresh)
      do lam=1,ncf
         do my=1,ncf
c -- check magnitude and convert to integers
            val = tao(my,lam)
            If(abs(val).ge.dblmax) Then
               b = val*dblinv
               if(abs(b).ge.d1max) then
                  dfac = abs(val)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
                  dfac = LOG10(dfac)
                  i1 = -NINT(dfac+0.5d0)
                  val = val*10.0d0**i1
                  jnt1(my,lam) = i1
                  jbuf(my,lam) = val
               else
                  i1 = abs(b)
                  b = val - SIGN(i1*dblmax,val)
                  jnt1(my,lam) = i1
                  jbuf(my,lam) = b
               endif
            Else
               jnt1(my,lam) = 0
               jbuf(my,lam) = val
            EndIf
         enddo
      enddo
      end
C=====gettkj====================================================
      subroutine gettkj(k,j,ndisk3,nvrsq,jbuf,
     1                  nrcpf,Tkj,nvir,thresh)
C
C    reads an amplitude matrix Tkj from the direct access file.
C    if j.gt.k  the matrix Tjk is read and then transposed.
C    the matrices are kept as integers on the direct access file
C    in  virtual MO basis (nvir*nvir)
C    the matrix is scaled by the integral threshold and converted
C    to real*8 before returned.
C
C    Svein Saebo Faytteville, AR June 2002
C
C    Arguments:
C    k,j     orbital indices
C    ndisk3  unit number
C    nvrsq   = nvir*nvir
C    jbuf    integer array for temporary storage dimension: nvir*nvir
C    nrcpf   number of records on each direct access file
C    Tkj     result is returned in this matrix.
C    nvir    number of virtual orbitals
C    thresh  integral threshold normally 1.0d-09
C
      implicit real*8(a-h,o-z)
      dimension jbuf(nvrsq),Tkj(*)
C
      if(k.ge.j) then
        kj=k*(k-1)/2+j
        ifil4=(kj-1)/nrcpf
        kjrec=kj-ifil4*nrcpf
        ndisk4=ndisk3+ifil4
        read(ndisk4,rec=kjrec)jbuf
        do imo=1,nvrsq
        Tkj(imo)=thresh*jbuf(imo)
        enddo
C
C    else if j.gt.k
C
      else
        jk=j*(j-1)/2+k
        ifil4=(jk-1)/nrcpf
        jkrec=jk-ifil4*nrcpf
        ndisk4=ndisk3+ifil4
        read(ndisk4,rec=jkrec) jbuf
        do imo=1,nvrsq
        Tkj(imo)=thresh*jbuf(imo)
        enddo
        call transpos(Tkj,nvir)
      endif
      end
C====transpose=============================================
      subroutine transpos(A,n)
C     calculates the transpose of a square matrix dimension n*n
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(n,*)
      do  I=1,N
      do  J=1,I
      HOLD=A(I,J)
      A(I,J)=A(J,I)
      A(J,I)=HOLD
      enddo
      enddo
      END
C====atoat=====================================================
      subroutine atoat(A,ncf,trnp)
C     **** input a matrix A
C     **** output 2A-A+ (Atilda)   returned in A

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      character*1 trnp
c     common/big/bl(30000)
      dimension A(*)
      parameter(two=2.0d0)
C
      nqd=ncf*ncf
      call getmem(nqd,itmst)
      itms=itmst-1
      do icp=1,nqd
      bl(itms+icp)=a(icp)
      end do
C
      call transpos(a,ncf)
      if(trnp.eq.'y'.or.trnp.eq.'Y') then
        do iadd=1,nqd
        a(iadd)=two*a(iadd)-bl(itms+iadd)
        end do
      else
        do iadd=1,nqd
        a(iadd)=two*bl(itms+iadd)-a(iadd)
        end do
      endif
      call retmem(1)
      END
C====atoat2====================================================
      subroutine atoat2(A,ncf,trnp,iscs)
C     **** input a matrix A
C     **** output 2A-A+ (Atilda)   returned in A

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
      character*1 trnp
c     common/big/bl(30000)
      dimension A(*)
C   iscs =0 for conventional MP2, 1 for SCS, 2 for scaled MP2
      parameter(two=2.0d0)
c
      if(iscs.gt.0) then
         call getrval('p1',p11)
         call getrval('p2',p2)
         p1=p11+p2
      endif
C
      nqd=ncf*ncf
      call getmem(nqd,itmst)
      itms=itmst-1
      do icp=1,nqd
         bl(itms+icp)=a(icp)
      end do
C
      if(iscs.eq.0) then
         call transpos(a,ncf)
         if(trnp.eq.'y'.or.trnp.eq.'Y') then
            do iadd=1,nqd
               a(iadd)=two*a(iadd)-bl(itms+iadd)
            end do
         else
            do iadd=1,nqd
               a(iadd)=two*bl(itms+iadd)-a(iadd)
            end do
         endif
      else
         call transpos(a,ncf)
         if(trnp.eq.'y'.or.trnp.eq.'Y') then
            do iadd=1,nqd
               a(iadd)=p1*a(iadd)-p2*bl(itms+iadd)
            end do
         else
            do iadd=1,nqd
               a(iadd)=p1*bl(itms+iadd)-p2*a(iadd)
            end do
         endif
      endif
      call retmem(1)
      END
C====tplustt====================================================
      subroutine tplustt(A,ncf)
C     **** input a matrix A
C     **** output T+Ttrsnp returned in A

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c     common/big/bl(30000)
      dimension A(*)
C
      nqd=ncf*ncf
      call getmem(nqd,itmst)
      itms=itmst-1
      do icp=1,nqd
         bl(itms+icp)=a(icp)
      end do
      call transpos(bl(itmst),ncf)
      do iadd=1,nqd
         a(iadd)=a(iadd)+bl(itms+iadd)
      end do
      call retmem(1)
      END
C=======moveTsh============================================
      subroutine moveTsh(tao,Tshells,nqd,my,lam,mys,las)
      implicit real*8(a-h,o-z)
C
C    Takes one fully backtransformed amplitude matrix Tmylam(ny,isi)
C    and store it in a larger array Tshells(nqs,LAS,MYS)
C
C   Svein Saebo Fayetteville Ar, summer 2002
C
C    tao     INPUT matrix dimension ncf*ncf
C    Tshells OUTPUT matrix dimension(nqd,las,mys)
C    nqd = ncf*ncf
C    my, lam  fixed AO indices within shells MYS and LAS
C    from 1 to MYS, and  to LAS
C    mys,las dimesnions of fixed shells
C
C    called from rphas2 once for each pair of AOS
C
      dimension tao(*),Tshells(nqd,las,mys)
      do imove=1,nqd
         Tshells(imove,lam,my)=tao(imove)
      enddo
      end
C=======mkdscl=================================================
      subroutine mkdscl(DSA,DSX,DS,D,ncfsq,
     1                  ndens,inx,ncf,ncs,densmax)
      implicit real*8(a-h,o-z)
      dimension DSA(*),DSX(ncf,*),D(ncfsq,*),inx(*),DS(ncs,*)
      parameter(zero=0.0d0,one=1.0d0)
      do imo=1,ncfsq
      ix=idamax(ndens,D(imo,1),ncfsq)
      DSA(imo)=abs(D(imo,ix))
      enddo
C construct scaling density in terms of shells
c
      densmax=zero
      ipoint=0
      do ics=1,ncs
        call get_shell_size(inx,ics,isize)
        kpoint=0
        do kcs=1,ncs
          call get_shell_size(inx,kcs,ksize)
          dmxs=zero
          do my=1,isize
          do lam=1,ksize
            dmx=DSX(ipoint+my,kpoint+lam)
            dmxs=max(dmxs,dmx)
          end do     !lam in kcs
          end do     ! my in ics
          kpoint=kpoint+ksize
ckw
          dmxs=min(dmxs,one)
          densmax=max(densmax,dmxs)
          DS(ics,kcs)=dmxs*dmxs
cc          DS(kcs,ics)=DS(ics,kcs)
        end do     ! kcs
        ipoint=ipoint+isize
      end do     !ics
      densmax=densmax*densmax
cccccc
cc      write(6,*) 'In <mkdscl>  DS array is:'
cc      call prntmat(ncs,ncs,ncs,DS)
cccccc
      end
C=======mkdsqt=================================================
      subroutine mkdsqt(DSA,VEC,DS,nval,ncf,
     1                  C,inx,ncs,densmax)
      implicit real*8(a-h,o-z)
      dimension DSA(ncf,*),VEC(*),inx(*),DS(ncs,*),C(ncf,*)
      parameter(zero=0.0d0,one=1.0d0)
      do iao=1,ncf
      ix=idamax(nval,C(iao,1),ncf)
      VEC(iao)=abs(C(iao,ix))
      enddo
C
      do iao=1,ncf
      do jao=1,iao
      DSA(iao,jao)=max(VEC(iao),VEC(jao))
      DSA(jao,iao)=DSA(iao,jao)
      enddo
      enddo
C construct scaling density in terms of shells
c
      ids=1
      densmax=zero
      ipoint=0
      do ics=1,ncs
      call get_shell_size(inx,ics,isize)
      kpoint=0
      do kcs=1,ncs
      call get_shell_size(inx,kcs,ksize)
      dmxs=zero
      do my=1,isize
      do lam=1,ksize
      dmx=DSA(ipoint+my,kpoint+lam)
      dmxs=max(dmxs,dmx)
      end do     !lam in kcs
      end do     ! my in ics
      kpoint=kpoint+ksize
ckw
      dmxs=min(dmxs,one)
      densmax=max(densmax,dmxs)
      DS(ics,kcs)=dmxs*dmxs
cc      DS(kcs,ics)=DS(ics,kcs)
      end do     ! kcs
      ipoint=ipoint+isize
      end do     !ics
cccccc
cc      write(6,*) 'In <mkdsqt>  DS array is:'
cc      call prntmat(ncs,ncs,ncs,DS)
cccccc
      end
C================putZic================================================
      subroutine putZic(zic,yoo,nmo,nval,ncore,ei)
      implicit real*8(a-h,o-z)
      dimension zic(nval,*),yoo(nmo,*),ei(*)
      do ic=1,ncore
         do ii=1,nval
            zic(ii,ic)=(yoo(ic,ncore+ii)-
     1           yoo(ncore+ii,ic))/(ei(ii+ncore)-ei(ic))
         enddo
      enddo
      end
C================makzb=================================
      subroutine  makzb(zb,zai,ei,ncf,nvir,nval)
      implicit real*8(a-h,o-z)
      dimension zb(ncf,ncf),zai(nvir,nval),ei(*)
      do ii=1,nval
      do jj=1,nval
      zb(ii,jj)=0
      if (ii.eq.jj) zb(ii,jj)=2*ei(ii)
      enddo
      enddo
      do ii=1,nval
      do ia=nval+1,nval+nvir
      zb(ii,ia)=zai(ii,ia-nval)*ei(ii)
      zb(ia,ii)=0
      enddo
      enddo
      do ia=nval+1,nval+nvir
      do ib=nval+1,nval+nvir
      zb(ia,ib)=0
      enddo
      enddo
      call matscal('Zb',-1.0d0)
      call matdef('CX','q',ncf,ncf)
      call matmmult('cano','Zb','CX')
      call matmmul2('CX','cano','W','n','t','a')
      call matrem('CX')
      end
C================formxnyi=================================
      subroutine FormXnyi(Y,T,aojm,nval,ncol,icol,ncf)

      use memory

      implicit real*8(a-h,o-z)
c     common/big/bl(30000)
C  accumulates D1 contributions in matrix Y ( Eq. 46)
      dimension Y(ncf,*),T(nval,*),aojm(nval,*),icol(*)
      call matdef('tqtr','r',nval,ncol)
      itqta=mataddr('tqtr')-1
      call matmmul2('Tmyl','quartra', 'tqtr','t','n','n')
      ix=0
      do ia=1,ncol
      ny=icol(ia)
      do ii=1,nval
      ix=ix+1
      Y(ny,ii)=Y(ny,ii)+bl(itqta+ix)
      enddo
      enddo
      call matrem('tqtr')
      end
C================formxny2=================================
      subroutine FormXny2(Y,T,aojm,nval,nrow,irow,ncf)

      use memory

      implicit real*8(a-h,o-z)
c     common/big/bl(30000)
C  accumulates D1 contributions in matrix Y ( Eq. 46)
      dimension Y(ncf,*),T(nval,*),aojm(nrow,*),irow(*)
      call matdef('tqtr','r',nval,nrow)
      itqta=mataddr('tqtr')-1
      call matmmul2('Tmyl','quartr2', 'tqtr','n','t','n')
      ix=0
      do ia=1,nrow
      ny=irow(ia)
      do jj=1,nval
      ix=ix+1
      Y(ny,jj)=Y(ny,jj)+bl(itqta+ix)
      enddo
      enddo
      call matrem('tqtr')
      end
C===============makzcit==================================
      subroutine makzict(Zict,Zic,ei,nval,ncore)
      implicit real*8(a-h,o-z)
C  calculates zict = zic*(ec+ei)/2 used for contribution from
C  frozen core
C
C  SS, April 2003
C
      dimension Zict(nval,*),ei(*),Zic(nval,*)
      parameter (half=0.5d0)
      do imo=1,nval
      do icore=1,ncore
      zict(imo,icore)=Zic(imo,icore)*(ei(ncore+imo)+ei(icore))*half
C     zict(imo,icore)=Zic(imo,icore)*ei(ncore+imo)
      enddo
      enddo
      end
C=====addtoT======================================================
      subroutine addtoT(T,D,DDT,ncf,my,lam)
      implicit real*8(a-h,o-z)
C
C   this subroutine adds extra terms to TAO before contracting with
C   integral derivatives.
C   We add here 2 contributions
C   contributions from the HF gradient:
C   contribution from terms contracted with Fx
C   only 2-electron contributions are handled here.
C   Now integral derivatives only need to be calculated once
C   (instead of 3 times: here + HF gradient + Fx)
C
C   this is quite expensive can it be improved?
C
C   I have made my best effort to move as many operation as far
C   out as possible.
C   this makes it hard to figure out what we are doing.
C   please look in 'addtot_orig'  which does the same as this
C   subroutine.
C
C   Svein Saebo, Fayetteville AR, July 2003
C
C   Arguments:
C   INPUT : T in AO basis (fully backtransformed)
C   OUTPUT: T input matrix plus extra terms to include HF gradient and Fx
C           terms
C   D density matrix 'den0'
C   DDT= 2X - 2C+AC
C
C   called from rphas2
C
      dimension T(ncf,*),D(ncf,*),DDT(ncf,*)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)
      fact=one4
      if(my.eq.lam) fact=one8
      dml=D(my,lam)*half*fact
      ddtml=DDT(my,lam)*one8
      do ny=1,ncf
      dmyny=D(ny,my)*fact
      ddtmyny=DDT(ny,my)*one4
      ddtlany=DDT(lam,ny)*one4
      do isi=1,ncf
      ttr2=ddtml*D(isi,ny)
      add1=
     1 dmyny*D(lam,isi)-dml*D(isi,ny)
      add2=
     2 ddtmyny*D(lam,isi)-ttr2
      T(isi,ny)=T(isi,ny)+add1+add2
      if(my.gt.lam) then
      add3=
     2 ddtlany*D(my,isi)-ttr2
      T(ny,isi)=T(ny,isi) + add3
      endif
      enddo
      enddo
      end
C=====addtoT1=====================================================
      subroutine addtoT1(T,D,DDT,ncf,my,lam)
      implicit real*8(a-h,o-z)
C
C   this subroutine adds extra terms to TAO before contracting with
C   integral derivatives.
C   We add here 2 contributions
C   contributions from the HF gradient:
C   contribution from terms contracted with Fx
C   only 2-electron contributions are handled here.
C   Now integral derivatives only need to be calculated once
C   (instead of 3 times: here + HF gradient + Fx)
C
C   this is quite expensive can it be improved?
C
C   I have made my best effort to move as many operation as far
C   out as possible.
C   this makes it hard to figure out what we are doing.
C   please look in 'addtot_orig'  which does the same as this
C   subroutine.
C
C   Svein Saebo, Fayetteville AR, July 2003
C
C   Arguments:
C   INPUT : T in AO basis (fully backtransformed)
C   OUTPUT: T input matrix plus extra terms to include HF gradient and Fx
C           terms
C   D density matrix 'den0'
C   DDT= 2X - 2C+AC
C
C   called from rphas2
C
C=============
      dimension T(ncf,ncf),D(ncf,*),DDT(ncf,*)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)
C NZG
C      T=0.0D0
C NZG
      fact=one4
      if(my.eq.lam) fact=one8
      dml=D(my,lam)*half*fact+DDT(my,lam)*one8
C     dml=DDT(my,lam)*one8
      ddtml=DDT(my,lam)*one8
c
      do ny=1,ncf
         dmyny=D(ny,my)*fact+DDT(ny,my)*one4
         do isi=1,ncf
            T(isi,ny)=T(isi,ny)+dmyny*D(isi,lam)-dml*D(isi,ny)
         enddo
      enddo
c
      if(my.le.lam) RETURN
      do isi=1,ncf
         disimy=D(isi,my)*one4
         do ny=1,ncf
            T(ny,isi)=T(ny,isi)+disimy*DDT(ny,lam)-ddtml*D(ny,isi)
         enddo
      enddo
      end
C=====addtoTx=====================================================
      subroutine addtoTx(T,D,DDT,ncf,my,
     1                  lam,mylam,lind)
      implicit real*8(a-h,o-z)
C
C    same as addtot but using triagular matrices
C    not as efficient but needs less memory
C
C   this subroutine adds extra terms to TAO before contracting with
C   integral derivatives.
C   We add here 2 contributions
C   contributions from the HF gradient:
C   contribution from terms contracted with Fx
C   only 2-electron contributions are handled here.
C   Now integral derivatives only need to be calculated once
C   (instead of 3 times: here + HF gradient + Fx)
C
C   this is quite expensive can it be improved?
C
C   I have made my best effort to move as many operation as far
C   out as possible.
C   this makes it hard to figure out what we are doing.
C   please look in 'addtot_orig'  which does the same as this
C   subroutine.
C
C   Svein Saebo, Fayetteville AR, July 2003
C
C   Arguments:
C   INPUT : T in AO basis (fully backtransformed)
C   OUTPUT: T input matrix plus extra terms to include HF gradient and Fx
C           terms
C   D density matrix 'den0'
C   DDT= 2X - 2C+AC
C
C   called from rphas2
C
      dimension T(ncf,*),D(*),DDT(*),lind(*)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)
      fact=one4
      if(my.eq.lam) fact=one8
      dml=D(mylam)*half*fact
      ddtml=DDT(mylam)*one8
      myin=lind(my)
      lain=lind(lam)
      do ny=1,ncf
      nyin=lind(ny)
      myny=max0(myin,nyin)+min0(my,ny)
      lamny=max0(lain,nyin)+min0(lam,ny)
      dmyny=D(myny)*fact
      ddtmyny=DDT(myny)*one4
      ddtlany=DDT(lamny)*one4
      do isi=1,ncf
      isin=lind(isi)
      lamsi=max0(lain,isin)+min0(lam,isi)
      nysig=max0(nyin,isin)+min0(ny,isi)
      ttr2=ddtml*D(nysig)
      T(isi,ny)=T(isi,ny)+dmyny*D(lamsi)-dml*D(nysig)+
     1 ddtmyny*D(lamsi)-ttr2
      if(my.gt.lam) then
        mysig=max0(myin,isin)+min0(my,isi)
        T(ny,isi)=T(ny,isi)+ddtlany*D(mysig)-ttr2
      endif
      enddo
      enddo
      end
C==========addtot_orig===================================
      subroutine addtoT_orig(T,D,DDT,ncf,my,lam)
      implicit real*8(a-h,o-z)
      dimension T(ncf,*),D(ncf,*),DDT(ncf,*)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)
      do ny=1,ncf
      do isi=1,ncf
      if(my.gt.lam) then
      T(isi,ny)=T(isi,ny)
     1 +one4*D(my,ny)*D(lam,isi)-one8*D(my,lam)*D(ny,isi)  !  SCF
     2 +one4*DDT(my,ny)*D(lam,isi)-one8*DDT(my,lam)*D(ny,isi)! Fx
      T(ny,isi)=T(ny,isi)
     2 +one4*DDT(lam,ny)*D(my,isi)-one8*DDT(my,lam)*D(ny,isi)! Fx
      else
       T(isi,ny)=T(isi,ny)
     1 +one8*D(my,ny)*D(lam,isi)-on16*D(my,lam)*D(ny,isi)
     2 +one4*DDT(my,ny)*D(lam,isi)-one8*DDT(my,lam)*D(ny,isi)
      endif
      enddo
      enddo
      end
C================BackTrans============================================
      subroutine BackTrans(TAO,ncf,nval,ipr,jpr,
     2                     ipoint,jpoint,CO)

      use memory

      implicit real*8(a-h,o-z)
      dimension TAO(ncf,*),ipoint(*),jpoint(*)
      dimension CO(ncf,*)
c     common /big/bl(30000)
C
C   this subroutine replaces
C     call matsimtr('TMO','tocc','TAO')
C   in rphas2, speeding up second transformation by using
C   compacted matrices
C
C   Svein Saebo, Fayetteville AR, July 2003
C
C  ARGUMENTS
C    OUTPUT:  TAO bactransformed matrix dimension ncf,ncf)
C    ncf  number of contracted basis functions
C    nval number of valence orbitals
C    ipr dimension of compacted coefficient matrix right side
C    jpr dimension of compacted coefficient matrix left side
C    ipoint pointers right side (see comments on Transetup)
C    jpoint pointers left side
C    CO canonical orbitals (ncf,nval)
C
C    note the pointers ipoint and jpoint (and ipr and jpr)
C    are generated in subroutine Transetup
C
C    called from rphas2
C
      call matzero('TAO')
      if(jpr.eq.0.or.ipr.eq.0) return
C
      call matdef('Thalf','r',jpr,nval)
      ittmp=mataddr('Thalf')
      call matdef('ctrunc','r',jpr,nval)
      ictru=mataddr('ctrunc')
      call CompactCoef(ncf,nval,CO,jpr,jpoint,bl(ictru))
      call matmmult('ctrunc','TMO','Thalf')
c     call matmmult('occu','TMO','Thalf')
      call matrem('ctrunc')
C   compress and perform matrix multiplication with compact matrices
      call matdef('ttemp','r',jpr,ipr)
      ittmp=mataddr('ttemp')
      call matdef('ctrunc','r',ipr,nval)
      ictru=mataddr('ctrunc')
      call CompactCoef(ncf,nval,CO,ipr,ipoint,bl(ictru))
      call matmmul2('Thalf','ctrunc','ttemp','n','t','n')
C   Scatter back to full AO dimension
      call Scatmylam(ncf,TAO,bl(ittmp),ipr,jpr,ipoint,jpoint)
      call matrem('ctrunc')
      call matrem('ttemp')
      call matrem('Thalf')
      end
C===============Scatmylam=====================================
      subroutine Scatmylam(ncf,TAO,TTemp,ipr,jpr,ipoint,jpoint)
C
C   backtransformation of Tmylam(i,j) to AO basis is carried out
C   with compacted coeffcient matrices.
C   result is a compacted matrix in AO basis of dimension
C   (jpr,ipr).  This subroutine scatters the compact matrix to
C   TAO if full AO basis.
C   TAO is zeroed out in Backtrans
C
C   called from Backtrans
C
C   SS Aug. 2003
C
C  INPUT: TTemp(jpr,ipr)
C  OUTPUT: TAO(ncf,ncf)
C
C  please refer to subroutine backtrans for remaining arguments
C
      implicit real*8(a-h,o-z)
      dimension TAO(ncf,*),TTemp(jpr,*),ipoint(*),jpoint(*)
      do itrunc=1,ipr
         do jtrunc=1,jpr
            ny=jpoint(jtrunc)
            isi=ipoint(itrunc)
            TAO(ny,isi)=Ttemp(jtrunc,itrunc)
         enddo
      enddo
      end
C==============TranSetup=================================
      subroutine TranSetup(CT,ncf,nval,ipr,jpr,
     1                     ipoint,jpoint,schw,ncs,MYS,
     2                     LAS,ancf,mapf2s,thresh)
      implicit real*8(a-h,o-z)
C   this subroutine construct array ipoint, which in turn is used in
C   subroutine backtrans.
C
C    Svein Saebo Fayetteville, AR July 2003
C
C  ARGUMENTS
C   CT transposed orbital coefficients
C   ncf number of contracted functions
C   nval number of valence orbitals
C   ipr  number of cbfs retained for right side
C   jpr  number of cbfs retained for left side
C   sigma = ipoint(ip), ip=1,ipr  pointers to basis functions
C   ny = jpoint*jp), jp=1,jpr
C   schw(ncs,ncs)  Schwarz integrals
C   ncs  number of contracted shells
C   MYS  fixed shell (my)
C   LAS  fixed shell (lam)
C   ancf variable to calculate % of AO retained
C   mapf2s(ncf) array mapping a function to a shell
C
C   OUTPUT
C   ipr and array ipoint
C
      dimension CT(nval,*),ipoint(*),jpoint(*),schw(ncs,*)
      dimension mapf2s(*)
C
      thre1=thresh
C  use the normal integral threshold for the time being?
      ipr=0
      jpr=0
      do isi=1,ncf
         ishel=mapf2s(isi)
         jmax=idamax(nval,CT(1,isi),1)
         cjsmx=abs(CT(jmax,isi))
         if (cjsmx*schw(MYS,ishel).gt.thre1) then
            ipr=ipr+1
            ipoint(ipr)=isi
         endif
         if (cjsmx*schw(LAS,ishel).gt.thre1) then
            jpr=jpr+1
            jpoint(jpr)=isi
         endif
      enddo
      ancf=ancf+(ipr+jpr)/2
      end
C==================================================
C  a few print routines used for debugging:
C===============Schwp===============================
      subroutine schwp(schw,ncs)
      implicit real*8(a-h,o-z)
      dimension schw(ncs,ncs)
C  prints Schwarz integrals
      do ip=1,ncs
      write(6,*) (schw(ip,iw),iw=1,ncs)
      enddo
      end
c================mapf2s_print==========================================
      subroutine mapf2s_print( mapf2s,ncs,ncf)
      implicit real*8(a-h,o-z)
      dimension mapf2s(ncf)
      do icf=1,ncf
         ics=mapf2s(icf)
         write(6,*)' icf=',icf,' ics=',ics
      enddo
      end
C================printb====================================
      subroutine  printb(Tbig,ncfsq,isize,ksize)
      implicit real*8(a-h,o-z)
      dimension Tbig(ncfsq,ksize,isize)
      do ich=1,isize
      do kch=1,ksize
      write(6,*) ' basispair',ich,kch
      write(6,*) (Tbig(ipm,kch,ich),ipm=1,ncfsq)
      enddo
      enddo
      end
C============cphf=========================================
      subroutine cphfz(ncf,    nval,   nvir,   ei,     ea,
     1                 y,      z,      nmo,    thresh, inx,
     2                 iprint, ncore,  ndisk,  ndisk2)
C
      use memory
      implicit real*8(a-h,o-z)
C
C   main routine for CPHF for MP2-gradients
C
C   SS. April 2003
C   Modified subroutine according to Eqs. 78-81
C   This will eliminate subroutine ztermsn
C
C   Svein Saebo, Fayetteville, AR summer 2002
C   Svein Saebo, Fayetteville, AR summer 2003
C   Svein Saebo, Fayetteville, AR summer 2006
C
C   determines z(a,i) virtual-occupied block,
C   input matrix y in Y
C   output marix z in z
C
C   this subroutine is called once from mp2_grad
C
C   Arguments:
C   ncf           number of contracted basis functions
C   nval          number of valence orbitals
C   nvir          number of virtual orbitals
C   ei(*)         orbital energies occupied prbital
C   ea(*)         orbital energies virtual orbitals
C   y(nvir,nval)  INPUT matrix y.  Note this matrix is generated from the
C                 original Y(ncf,nval):  y=Cv *Y (done in the calling
C                 program)
C   z(nvir,nval)  OUTPUT matrix z
C   nmo           number of occupied mos
C   thresh        integral threshold (normally 1.0e-10)
C   inx           aray with contraction info
C   iprint        printlevel
C   ncore         number of core orbitals
C   ndisk         direct access unit for old iterated (Z's)
C   ndisk2        direct access unit for old G(z)'s
C
C   External Calls:
C   MakeGma2 :    Calculates two-electron part of the Fock matrix
C                 using any symmertical matrix as density
C   updatzk       determines new z's during CPHF-iterations
C                 see subroutine updatzk for formula.
C
c     common /big/bl(30000)
      dimension ei(*),ea(*),y(nvir,*),z(nvir,nmo)
      integer*4 info4
      dimension inx(*)
      logical done,swtr,efit,dz2z
      parameter (half=0.5d0,one=1.0d0,two=2.0d0,four=4.0d0,maxit=40)
C   max CPHF iterations set to 40 here!
      parameter (accur=1.0d-05,accu2=0.000416d0,accu3=0.0000416d0) 
c     accu2 is threshold to switch integral threshold
C
      done=.false.
      swtr=.true.
      diffm=0.d0
      tcphfz=0.0d0
      efit=.false.
      dz2z=.false.
      itsw=7
      if(ncf.gt.600)itsw=6
c
c -- Allocate memory
      nstrip=nvir*nmo
      ictr=igetival('ictr')
      call matdef('Gmat','s',ncf,ncf)
      call matdef('DP','s',ncf,ncf)
      call matdef('GZ','r',nvir,nmo)   ! this is G(z)
      call matdef('dptm','q',ncf,ncf)
      iadmp=mataddr('dptm')
      call matdef('DelZ','r',nvir,nmo)
      idelz=mataddr('DelZ')
      call matdef('Gzold','r',nvir,nmo)
      call matdef('ttxx','r',ncf,nmo)
C  reserve memory for two copies of H and B used in DIIS
      call getmem(maxit*maxit,ih1)
      call getmem(maxit*maxit,ih2)
      call getmem(maxit,iba)
      call matdef('DZ','r',nvir,nmo)
      iadz=mataddr('DZ')
      call matdef('Zaio','r',nvir,nmo)
      izaio=mataddr('Zaio')
C
C     end memory allocation
c----------------------------------------------------------------
      thres1=min(thresh,1.0d-10)
      thres2=thres1*1000.d0           ! loose integral threshold
      thresX=thres2
      if(.not.swtr)thresX=thres1
C  the first iterations with thres2 normally 10**-7, then thres1
C  normally 10**-10
C----------------------------------------------------------------
C  get initial Z's ! current z is z (z is the same as 'Zai')
      do ia=1,nvir
         deno=ea(ia)
         do ii=1,nmo
            z(ia,ii)=-y(ia,ii)/(deno-ei(ii))
         enddo
      enddo
C  initial Delta Z is simply Z
      call matcopy('Zai','DelZ')
C  start CPHF-iterations
      write(6,*) ' Start CPHF'
      write(6,*)
      write(6,*) ' Initial Integral Threshold:  ',thresX
      write(6,*)
      write(6,*)  ' Iter    Max DeltaZ  Elapsed Time   Threshold'
c----------------------------------------------------------------
C    Iterations start here...........
      call writez(ndisk,1,z,nstrip)
      icpdi=0
      icptx=0
      call secund(tcpz1)
      do icpit=1,maxit
         call elapsec(tcpit1)
         icpdi=icpdi+1
         icptx=icptx+1
C  construct  G(z)
C    DP=1/2(Cv * z * Co+ + Co * z+ * Cv+)
C  we are now using DeltaZ instead of Z
         call matmmult('virt','DelZ','ttxx')
         call matmmul2('ttxx','occa','dptm','n','t','n')
         call tplustt(bl(iadmp),ncf)
         call matscal('dptm',half)
         call matcopy('dptm','DP')
         igadr=mataddr('Gmat')
         idsa=mataddr('DP')
C  construct two-electron part of the Fock matrix with DP as density
C  result in 'Gmat' = bl(igadr)
         call mmark
         call MakeGma2(ncf,nmo,nval,nvir,thresX,
     1                 bl,bl(ictr),bl(idsa),bl(igadr))
         call retmark
         call matdef('ztmp','r',nvir,ncf)
         call matmmul2('virt','Gmat','ztmp','t','n','n')
         call matmmult('ztmp','occa','GZ')
         call matrem('ztmp')
         call matscal('GZ',four)
         igz=mataddr('GZ')
         igzo=mataddr('Gzold')

C    COVERGENCE ACCELEFRATION--------------------------------
C   copy previous H into current H with correct dimensions
         if(icpdi.gt.1) call Hcopy(bl(ih1),bl(ih2),icpdi)
C   calculate  new elements of H
         call makdz(z,bl(iadz),nvir,nmo,ea,ei)
C   need Gz(zn)  this can be calculated as
C   Gz(deltaz)+Gz(zn-1)
         if (icpdi.eq.1.or.icptx.eq.1) then
            call matcopy('GZ','Gzold')
         else
            call readz(ndisk2,icpdi-1,bl(igzo),nstrip)
            call matadd('GZ','Gzold')
         endif
         call writez(ndisk2,icpdi,bl(igzo),nstrip)
         call newHB(icpdi,bl(ih1),bl(iba),nstrip,ndisk,
     1              bl(igzo),z,bl(iadz),y)
C     call printH(bl(ih1),bl(iba),icpdi)
C   now ready to solve set of linear equations
C   first make a copy of H
         call dcopy(icpdi*icpdi,bl(ih1),1,bl(ih2),1)
C
         call getint(icpdi,ipiv)
         call DGESV(icpdi,1,bl(ih1),icpdi,bl(ipiv),bl(iba),icpdi,info4)
         info=int(info4)
         call retmem(1)
C  error check
         if (info.ne.0) then
            if (info.lt.0) then
               call nerror(1,'DGESV','wrong argument',info,icpdi)
            else
               call nerror(2,'DGESV','singular matrix',info,icpdi)
            endif
         endif
C  no errors ....continue
C  upon exit from DGESV H is destroyed and B contains the coefficients
C  construct improved z from these coefficients and old iterates
C
         if (iprint.ge.4) then
            write(6,*) ' Coefficients:'
            write(6,12)(bl(iba-1+iw),iw=1,icpdi)
   12       format(1x,5f10.7)
         endif
c
C  construct  unprimed Z's:
         call matscal('Zai',bl(iba-1+icpdi))
         if (icpdi.gt.1) then
            do minit=1,icpdi-1
               call readz(ndisk,minit,bl(iadz),nstrip)
               call matadd1('DZ',bl(iba-1+minit),'Zai')
            enddo
         endif
         call writez(ndisk,icpdi,z,nstrip)
C if not first iteration compare with previous z
         if (icpdi.gt.1) then
            call readz(ndisk,icpdi-1,bl(idelz),nstrip)
            call matscal('DelZ',-1.0d0)
            call matadd('Zai','DelZ')
            ix=idamax(nvir*nmo,bl(idelz),1)
            diffm=abs(bl(idelz+ix-1))
            if (diffm.lt.accur) then
               done=.true.
               goto 200
            endif
C   if done is true we are finished
            if (swtr.and.(diffm.lt.accu2.or.icpit.ge.itsw)) then
               swtr=.false.
               icpdi=0
               call elapsec(tcpit2)
               write(6,11) icpit,diffm,(tcpit2-tcpit1)/60.0d0,thresX
               thresX=thres1
               write(6,*) ' Switching to integral threshold: ',thresX
               write(6,*) ' Restarting DIIS'
               call matcopy('Zai','DelZ')
               call writez(ndisk,1,z,nstrip)
               diffm=0.0d0
               icptx=0
               cycle
            endif
C       if(diffm.lt.accu3.and.icptx.ge.4) then
C       icptx=0
C       dz2z=.true.
C       endif
            if (icpdi.ge.itsw) then
               call elapsec(tcpit2)
               write(6,11) icpit,diffm,(tcpit2-tcpit1)/60.0d0,thresX
               write(6,*) 'Switching to full Z and restarting DIIS'
               call matcopy('Zai','DelZ')
               call writez(ndisk,1,z,nstrip)
               diffm=0.0d0
               icpdi=0
               dz2z=.false.
               icptx=0
               cycle
            endif
         endif
C  the H' and G' must be modified to H and G (fixH and fixGz)
C     call printH(bl(ih2),bl(iba),icpdi)
         call fixH(bl(iba),bl(ih2),bl(ih1),icpdi)
C     call printH(bl(ih1),bl(iba),icpdi)
         call dcopy(icpdi*icpdi,bl(ih1),1,bl(ih2),1)
         call fixGz(bl(iba),ndisk2,icpdi,nvir,nmo)
C  all DIIS stuff updated  goon...
         call readz(ndisk2,icpdi,bl(igzo),nstrip)
         iflg=1
         call updatzk(z,y,bl(igzo),ea,ei,nvir,nmo,iflg)
         call writez(ndisk,icpdi+1,z,nstrip)
         call readz(ndisk,icpdi,bl(idelz),nstrip)
         call matadd1('Zai',-1.0d0,'DelZ')
         call matscal('DelZ',-1.0d0)
         ix=idamax(nvir*nmo,bl(idelz),1)
         diffx=abs(bl(idelz+ix-1))
         if (diffx.lt.accur) then
            write(6,*) 'Convergence criterion met'
            done=.true.
            diffm=diffx
         endif
C
  200    continue
         call elapsec(tcpit2)
         write(6,11) icpit,diffm,(tcpit2-tcpit1)/60.0d0,thresX
   11    format(1x,i3,6x,f10.8,4x,1f6.3,5x,1e10.2)
C
         if (dz2z) then
            write(6,*) 'Switching to full Z'
            call matcopy('Zai','DelZ')
            dz2z=.false.
         endif
         if (done) then
            write(6,*) ' CPHF has converged '
            goto 100
         endif
C
      enddo    ! end of iteration loop
cc
      write(6,*) ' CPHF did not converge: Max diff :',diffm
  100 continue
      call secund(tcpz2)
      tcphfz=tcpz2-tcpz1
      if (efit) then  !  this is not needed currently disabled
         write(6,*)'Do a final full iteration'
C  construct  G(z)
         call matmmult('virt','Zai','ttxx')
         call matmmul2('ttxx','occa','dptm','n','t','n')
         call tplustt(bl(iadmp),ncf)
         call matscal('dptm',half)
         call matcopy('dptm','DP')
         igadr=mataddr('Gmat')
         idsa=mataddr('DP')
C  construct two-electron part of the Fock matrix with DP as density
C  result in 'Gmat' = bl(igadr)
         call mmark
         call secund(tcpz1)
         call MakeGma2(ncf,nmo,nval,nvir,thresX,
     1                 bl,bl(ictr),bl(idsa),bl(igadr))
         call secund(tcpz2)
         tcphfz=tcphfz+tcpz2-tcpz1
         call retmark
         call matdef('ztmp','r',nvir,ncf)
         call matmmul2('virt','Gmat','ztmp','t','n','n')
         call matmmult('ztmp','occa','GZ')
         call matrem('ztmp')
         call matscal('GZ',four)
         iflg=1
         call updatzk(z,y,bl(igz),ea,ei,nvir,nmo,iflg)
      endif
C
      call f_lush(6)
      call matrem('Zaio')
      call matrem('DZ')
      call retmem(3)
      write(6,13)  icpit,tcphfz/60.0d0
   13 format(1x,/,'Elapsed time for ',I2,' CPHF iterations:',f10.2,
     c            ' min.')
C
C  iterations finished
C  Calculate contributions to the gradient, the following replaces
C  subroutine ztermsn
C  Contribution from the first term  (Fx) of Eq. 76
C   the symmetrical matrix Z (DP here) should be added to
C   2X-2CoACo+ before contracting with F(x)  Eq. (77
      call matmmult('virt','Zai','ttxx')
      call matmmul2('ttxx','occa','dptm','n','t','n')
      call matrem('ttxx')
      call tplustt(bl(iadmp),ncf)
      call matscal('dptm',half)
      call matadd('dptm','DDT')
C  Frozen core:
      if (ncore.gt.0) then
         call matdef('zicao','q',ncf,ncf)
         call matdef('zictm','r',ncf,ncore)
         call matmmult('occu','Zic','zictm')
         call matmmul2('zictm','coro','zicao','n','t','n')
         izicao=mataddr('zicao')
         call tplustt(bl(izicao),ncf)
         call matscal('zicao',half)
         call matadd('zicao','DDT')
      endif
C  end of first term frozen core
C
C  also -1/2 DG(Z)D should be added to W before contracting with Sx
C  Eq. 80
C
      call matcopy('dptm','DP')
C    frozen core
      if (ncore.gt.0) call matadd('zicao','DP')
      iZts=mataddr('DP')
      call mmark
      call MakeGma2(ncf,nmo,nval,nvir,thresh,
     1              bl,bl(ictr),bl(iZts),bl(igadr))
      call retmark

      call matmmult('Gmat','den0','dptm')
      call matscal('dptm',-half)
      call matmmul2('den0','dptm','W','n','n','a')
C  end frozen core
C  end -1/2DG(Z)D term
C  form matrix ztilda and add -c*ztildaC+ to W (Eq. 78)
C
      Call matdef('Ztilda','q',ncf,ncf)
      call matdef('ztemp','r',ncf,nmo)
      call matdef('ztil','r',nvir,nmo)
      iztia=mataddr('ztil')
      nstep=-1
      do io=1,nmo
         do ia=1,nvir
            nstep=nstep+1
            bl(iztia+nstep)=z(ia,io)*ei(io)
         enddo
      enddo
      call matmmult('virt','ztil','ztemp')
      call matmmul2('ztemp','occa','Ztilda','n','t','n')
      izta=mataddr('Ztilda')
      call tplustt(bl(izta),ncf)
      call matadd1('Ztilda',-half,'W')
      call matrem('ztil')
      call matrem('ztemp')
      call matrem('Ztilda')
C   now frozen core
      if (ncore.gt.0) then
         call matdef('zict','r',nval,ncore)
         izica=mataddr('zict')
         izico=mataddr('Zic')
         izicao=mataddr('zicao')
         call makzict(bl(izica),bl(izico),ei,nval,ncore)
         call matmmult('occu','zict','zictm')
         call matmmul2('zictm','coro','zicao','n','t','n')
         call tplustt(bl(izicao),ncf)
         call matadd1('zicao',-half,'W')
         call matrem('zict')
         call matrem('zictm')
         call matrem('zicao')
      endif
C    end frozen core
      call matrem('Gzold')
      call matrem('DelZ')
      call matrem('dptm')
      call matrem('GZ')
      call matrem('DP')
      call matrem('Gmat')
      end
C==============Hcopy======================================
      subroutine Hcopy(H,HO,iter)
      implicit real*8(a-h,o-z)
      dimension H(iter,iter),HO(iter-1,iter-1)
      do jcp=1,iter-1
         do icp=1,iter-1
            H(icp,jcp)=HO(icp,jcp)
         enddo
      enddo
      end
C===============newHB==================================
      subroutine newHB(iter,H,B,nstrip,kdisk,GZ,Z,DZ,Y)
      use memory
      implicit real*8(a-h,o-z)
C    calculate new elements of matrix H, and vector B
C    H c =b
C    used in cphf
C
C    ss Aug 2003
C
c     common /big/bl(30000)
      dimension H(iter,*),B(*),Z(*),DZ(*),Y(*),GZ(*)
C  current z in Z and current Dz n DZ
C  save current c
      call getmem(nstrip,icop)
      call dcopy(nstrip,Z,1,bl(icop),1)
      do kiter=iter,1,-1
         if(kiter.lt.iter) call readz(kdisk,kiter,Z,nstrip)
         hik=ddot(nstrip,Z,1,DZ,1)+ddot(nstrip,Z,1,GZ,1)
         H(kiter,iter)=hik
         H(iter,kiter)=hik
         B(kiter)=-ddot(nstrip,z,1,y,1)
      enddo
C  restore Zai
      call dcopy(nstrip,bl(icop),1,Z,1)
      call retmem(1)
      end
C==========updatzk========================================
      subroutine updatzk(z,y,gz,ea,ei,nvir,nmo,iflag)
C  updates the z-matrix during cphf iterations
C
C  Svein Saebo Fayetteville, AR summer 2002
C
C  called from cphfz once for each iteration
C
C     OUTPUT   z   z-matrix dimension(nvir,nval)
C     INPUT    y   y-matrix dimension(nvir,nval)
C     INPUT    gz   1/4 G(DP)  where DP =(Co+)z(Cv)+(Cv+)z+(Co)
C     INPUT    ea,ei  vectors with orbital energies
C     nvir     number of virtual orbitals
C     nval     number of occupied orbitals
C     zold     used for temporary storag eof old z-matrix
C     diffm    large absolute difference between old and new z-mat
C
      implicit real*8(a-h,o-z)
      dimension z(nvir,*),y(nvir,*),gz(nvir,*)
      dimension ea(*),ei(*)
      if (iflag.eq.0) then
         do ia=1,nvir
            deno=ea(ia)
            do ii=1,nmo
               z(ia,ii)= - gz(ia,ii)/(deno-ei(ii))
            enddo
         enddo
      else
         do ia=1,nvir
            deno=ea(ia)
            do ii=1,nmo
               z(ia,ii)= - (y(ia,ii)+gz(ia,ii))/(deno-ei(ii))
            enddo
         enddo
      endif
      end
C=======makdz==============================
      subroutine makdz(z,dz,nvir,nmo,ea,ei)
      implicit real*8(a-h,o-z)
      dimension z(nvir,*),dz(nvir,*),ea(*),ei(*)
      do ia=1,nvir
         deno=ea(ia)
         do ii=1,nmo
            dz(ia,ii)=z(ia,ii)*(deno-ei(ii))
         enddo
      enddo
      end
C=======readz===============================
      subroutine readz(ndisk,nrec,zd,nstrip)
      implicit real*8(a-h,o-z)
      dimension zd(nstrip)
      read(ndisk,rec=nrec) zd
      end
C=======writez==============================
      subroutine writez(ndisk,nrec,zd,nstrip)
      implicit real*8(a-h,o-z)
      dimension zd(nstrip)
      write(ndisk,rec=nrec) zd
      end
C===============================================
      subroutine printH(H,B,idim)
      implicit real*8(a-h,o-z)
      dimension H(idim,idim), B(*)
      write(6,*) 'Matrix H'
      do ip=1,idim
      write(6,10) (H(ip,iw),iw=1,idim)
   10 format(1x,5f10.7)
      enddo
      write(6,*) ' Vector B'
      write(6,10) (B(iw),iw=1,idim)
      end
C===============================================
      subroutine Aterms(nval,   nvir,   ndisk1, ndisk3, iprint,
     $                  nvrsq,  nindxp, lbin,   nrcpb,  thresh,
     $                  iscs)

      use memory

      implicit real*8(a-h,o-z)
C
C   Calculation of the matrix A.
C   A is now calculated by first sorting Tij(ab) to Tab(ij).  This
C   allows calculating A using matrix multiplication. Since DGEMM is
C   significantly more efficient that DDOT for large systems this is
C   worth the extra effort.
C
C   Called from mp2_grad
C
C   Svein Saebo, Starkville, MS September 2003
C
c     common /big/bl(30000)
      parameter(sixty=60.0d0,two=2.0d0,onef=0.25d0)
C
      call secund(taso1)
      call elapsec(easo1)
C
C  reserve memory for one Tij
c     call getmem(nvrsq/2+1,ibuf)
      call getint_4(nvrsq,ibuf)
c     call getmem(nvrsq/8+1,i1)
      call getint_1(nvrsq,i1)
C  reserve memory for bins
c     call getmem(nindxp*lbin,ibins)
      call getint_4(nindxp*lbin*2,ibins)
c     call getmem(nindxp*lbin/4+1,i1bin)
      call getint_1(nindxp*lbin*2,i1bin)
C
      call Aphas1(nval,   nvir,   ndisk1, ndisk3, lbin,
     1            nrcpb,  iprint,bl(ibins),bl(i1bin),bl(ibuf),
     2            bl(i1))
      call retmem(4)
cc
      if(iprint.ge.2) then
        call secund(taso2)
        call elapsec(easo2)
        write(6,*)' CPU and elapsed time for sort for A-terms:'
        write(6,100) (taso2-taso1)/60.0d0,(easo2-easo1)/60.0d0
  100 format(1x,f8.2,' minutes ',f8.2,' minutes ')
      endif
C
C  sort finished - integrals now on <bins> file
C
      call matdef('Tab','q',nval,nval)
      call matdef('Tbar','q',nval,nval)
      itab=mataddr('Tab')
      itbar=mataddr('Tbar')
c     call getmem(lbin,ibin)
      call getint_4(lbin*2,ibin)
c     call getmem(lbin/4+1,i1)
      call getint_1(lbin*2,i1)
      call CalcA(nval,   nvir,   nrcpb,  ndisk3, lbin,
     2           thresh, iprint,bl(ibin),bl(i1), bl(itab),
     1           bl(itbar),iscs)
      call retmem(2)
      call matrem('Tbar')
      call matrem('Tab')
C
      return
      end
C================calca==========================================
      subroutine CalcA(nval,   nvir,   nrcpb,  ndisk3, lbin,
     2                 thresh, iprint, ibins,  i1bin,  Tab,
     2                 Tbar,iscs)
      implicit real*8(a-h,o-z)
      integer*1 i1bin(2,lbin)
      integer*4 ibins(2,lbin)
      dimension Tab(nval,*),Tbar(nval,*)
      parameter (half=0.5d0,dblmax=2147483648.0d0)
c
      dblcmp = dblmax*thresh
c
      iab=0      ! bin record counter
      nbin=0
c
      do ia=1,nvir
         do ib=1,ia
            iab=iab+1
C  read back matrix Tab from bins
            istar=(iab-1)*nrcpb
            iwrd=1
            mrec=istar+1
            read(ndisk3,rec=mrec) i1bin,ibins
            nbin=nbin+1
            ij=0
            do ii=1,nval
               do jj=1,ii
                  ij=ij+1
                  if (ij.gt.lbin) then
                     iwrd=iwrd+1
                     mrec=istar+iwrd
                     nbin=nbin+1
                     read(ndisk3,rec=mrec) i1bin,ibins
                     ij=1
                  endif
c -- decompress the integral ------------------------
                  If (i1bin(1,ij).eq.0) Then
                     xx = ibins(1,ij)*thresh
                  Else If(i1bin(1,ij).gt.0) Then
                     xx = ibins(1,ij)*thresh
                     xx = xx + SIGN(i1bin(1,ij)*dblcmp,xx)
                  Else
                     xx = ibins(1,ij)*thresh*10.0d0**(-i1bin(1,ij))
                  EndIf
                  Tab(ii,jj)=xx
                  If (i1bin(2,ij).eq.0) Then
                     xx = ibins(2,ij)*thresh
                  Else If(i1bin(2,ij).gt.0) Then
                     xx = ibins(2,ij)*thresh
                     xx = xx + SIGN(i1bin(2,ij)*dblcmp,xx)
                  Else
                     xx = ibins(2,ij)*thresh*10.0d0**(-i1bin(2,ij))
                  EndIf
                  Tab(jj,ii)=xx
c ---------------------------------------------------
               enddo  ! over jj
            enddo  ! over ii
            call matcopy('Tab','Tbar')
            call atoat2(Tbar,nval,'n',iscs)
            call matmmul2('Tbar','Tab','Aik','t','n','a')
            if (ia.gt.ib) then
               call matmmul2('Tbar','Tab','Aik','n','t','a')
            endif
         enddo
      enddo
      if(iprint.ge.6) call matprint('Aik',6)
      end
C===============Aphas1====================================================
      subroutine Aphas1(nval,   nvir,   ndisk1, ndisk3, lbin,
     1                  nrcpb,  iprint, ibins,  i1bin,  ibuf,
     2                  i1)

      use memory

      implicit real*8(a-h,o-z)
C
C  The residua are on disk in virtual basis, as one matrix Tij(ia,ib)
C  on ndisk1 stored as 5-byte integers in Canonical order without
C  indices
C
C  The memory is divided into* nvir*(nvir+1)/2 bins and a given integral
C  Tij(a,b) is placed in a bin iab as follows:
C     integral T(a,b) is put into the (1,ij) position of bin IAB
C     integral T(b,a) is put into the (2,ij) position of bin IAB
C  the bins can normally only hold -lbin- ij values
C  when the lbin have been sorted, all bins are are full and they
C  are written to direct access file as follows:
C  nrcpb is the number of records needed for each Tab
C  (~mpairs/lbin; mpairs=nvir*(nvir+1)/2)
C     record 1..nrcpb              T11
C     record nrcpb+1 ,   2*nrecpb  T21    etc...
C  Note when the nrcpb record are read back for a given a,b
C     Tab(i,j) = ibins(1,ij)
C     Tab(j,i) = ibins(2,ij)
C  Tab is only used here and discarded as soon as matrix A has
C  been calculated
C
C     INPUT     Tij in virtual basis on ndisk1
C     OUTPUT    Tab in valence basis on ndisk3
C
C     NOTE the sorted Tab are only needed here and will be discarded,
C     while Tij will be used later and must be kept
C
C  ARGUMENTS
C     nval       number of valence orbitals
C     nvir       number of virtual orbitals
C     ndisk1     unit number for Tij in virtual basis
C     ndisk3     unit number for bins
C     lbin       length of bin (each bin contains 2*lbin integers)
C     nrcpb      number of records per bin
C     iprint     print level
C     ibins      integer array for bins ibins(2,lbin,nindxp)
C     i1bin      integer*1 storage for precision overflow in bins
C     ibuf       integer array nvir*nvir long
C     i1         integer*1 storage for precision overflow in buffer
C
C  Svein Saebo Fayetteville Ar, and Starkville, MS Sept 2003
C
c     common /big/bl(30000)
      integer*1 i1bin(2,lbin,*),i1(nvir,nvir)
      integer*4 ibins(2,lbin,*),ibuf(nvir,nvir)
C
C
      nindxp=nvir*(nvir+1)/2
      nbin=0
      iwrb=0
C
C  loop over valence pairs
C
      ij=0
      iijj=0
      do ii=1,nval
         do jj=1,ii
            ij=ij+1
            iijj=iijj+1
C
C  write out the transformed bin if full
C
            if (ij.gt.lbin) then
               krec=-nrcpb+1
               do irec=1,nindxp
                  krec=krec+nrcpb
                  mrec=krec+iwrb
                  nbin=nbin+1
                  call WriteBin2(ndisk3,mrec,2*lbin,ibins(1,1,irec),
     $                           i1bin(1,1,irec))
               enddo
               iwrb=iwrb+1
               ij=1
            end if
C
C  read Tij in virtual basis and transform to AO basis
C
            read(ndisk1,rec=iijj) i1,ibuf
            iab=0
            do ia=1,nvir
               do ib=1,ia
                  iab=iab+1
                  ibins(1,ij,iab)=ibuf(ia,ib)
                  ibins(2,ij,iab)=ibuf(ib,ia)
                  i1bin(1,ij,iab)=i1(ia,ib)
                  i1bin(2,ij,iab)=i1(ib,ia)
               enddo
            enddo
c
         end do      !  jj loop
      end do      !  ii loop
C
C  write out the content of the bins at the end
C
      krec=-nrcpb+1
      do irec=1,nindxp
         krec=krec+nrcpb
         mrec=krec+iwrb
         nbin=nbin+1
         call WriteBin2(ndisk3,mrec,2*lbin,ibins(1,1,irec),
     $                  i1bin(1,1,irec))
      enddo
cc
      end
C===============WriteBin2=================================================
      subroutine WriteBin2(ndisk,nbin,lbin,ibins,ibin1)
c  This routine writes a full bin on disk
c  Arguments:
c  ndisk:  unit number
c  nbin:   position of the bin to be written
c          so far, bins are written as they come. One could write
c          them in a better order
c  lbin:   number of records per bin
c  ibins:  compressed integral + indices
c  ibin1:  precision overflow integer
c
      integer*1 ibin1(lbin)
      integer*4 ibins(lbin)
      write(ndisk,rec=nbin) ibin1,ibins
      end
C===============setup_listreal========================================
      subroutine setup_listreal(natoms,listreal)
      implicit integer(a-z)
c
c  temporary routine to set up listreal array
c  simply sets each atom to itself as symmetry is currently
c  not used in the MP2 gradient
c
      integer listreal(natoms)
c
      Do 10 I=1,natoms
      listreal(I) = I
 10   Continue
c
      Return
      End
C=====================================================================
C====NOT USED ..FOR TESTING ONLY
      subroutine addto2(T,D,DDT,ncf,my,lam)
      implicit real*8(a-h,o-z)
      dimension T(ncf,*),D(ncf,*),DDT(ncf,*)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)
      parameter(on32=0.03125d0)
      do ny=1,ncf
      do isi=1,ncf
      if(my.gt.lam) then
      T(isi,ny)=T(isi,ny)
     2 +one4*DDT(my,ny)*D(lam,isi)-one8*DDT(my,lam)*D(ny,isi)! Fx
      T(ny,isi)=T(ny,isi)
     2 +one4*DDT(lam,ny)*D(my,isi)-one8*DDT(my,lam)*D(ny,isi)! Fx
      else
       T(isi,ny)=T(isi,ny)
     2 +one4*DDT(my,ny)*D(lam,isi)-one8*DDT(my,lam)*D(ny,isi)
      endif
      enddo
      enddo
      end
      subroutine prinTij(jbuf,TT,nvir,thresh)
      implicit real*8(a-h,o-z)
      integer*4 jbuf(nvir,nvir)
      dimension TT(nvir,*)
      do ii=1,nvir
      do jj=1,nvir
      TT(ii,jj)=thresh*jbuf(ii,jj)
      enddo
      enddo
      call matprint('TT',6)
      end
      subroutine getinfs(inx,ics,kcs,icf1,icf2,kcf1,kcf2)
      implicit integer(a-z)
      dimension inx(12,*)
      kcf1=inx(11,kcs)+1
      kcf2=inx(10,kcs)
      icf1=inx(11,ics)+1
      icf2=inx(10,ics)
      end
c=============ExtractEn===============================================
      subroutine ExtractEn(   ncf,   icf,  kcf,   icf1,  icf2,
     1                        kcf1,  kcf2,  xint,  xmat,  nrow,
     2                        ncol,  irow,  icol,  nrow1, ncol1,
     3                        irow1, icol1)
c  This routine extracts a single matrix from the integral array
c  The latter is indexed as
C  xint(l=1:ncf,j=1:ncf,k=kcf1:kcf2,ncf,i=icf1:icf2);
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
      use memory
      implicit real*8 (a-h,o-z)
      dimension xint(ncol,nrow,kcf1:kcf2,icf1:icf2),xmat(ncf,ncf)
      integer irow(ncf),icol(ncf)
      integer irow1(ncf),icol1(ncf)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
C
      do ll=1,ncol
      do jj=1,nrow
      xmat(jj,ll)=xint(ll,jj,kcf,icf)
      enddo
      enddo
c
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
      if(abs(xmat(jj,ll)).gt.one) then
c column ll has at least one non-zero element
        ncol1=ncol1+1
        icol1(ncol1)=l
        go to 100
      endif
      end do
      go to 150
 100  continue
      if(ll.gt.ncol1) then
        do jj=1,nrow
        xmat(jj,ncol1)=xmat(jj,ll)
        end do
      end if
 150  continue
      end do
      nrow1=0
      do jj=1,nrow
      j=irow(jj)
      do ll=1,ncol1
      if(abs(xmat(jj,ll)).gt.one) then
        nrow1=nrow1+1
        irow1(nrow1)=j
        go to 200
      endif
      end do
      go to 250
 200  continue
      if(jj.gt.nrow1) then
        do ll=1,ncol1
        xmat(nrow1,ll)=xmat(jj,ll)
        end do
      end if
 250  continue
      end do
      end
C==============MP2dip=========================================
      subroutine MP2DIP(dip,ncf,inx,ncs,natom)
      use memory
      implicit real*8(a-h,o-z)
C     common /big/bl(30000)
      dimension dip(3),inx(*)
      PARAMETER (Debye=0.39342658d0)
c  define (temporary) dipole matrix
      call matadd('den0','DDT')
      call matdef('dipp','s',ncf,ncf)
      call getival('inuc',inuc)
      call getival('ibas',ibas)
      call inton(3,natom,bl(mataddr('dipp')),inx,1,0,bl(ibas),
     1           bl(inuc),ncs)
ccc   call matprint('dipp',6)
      call matprodtr('DDT','dipp',dip(1))
C     call matwrite('dipp',np1,0,'aoX     ')
c
      call inton(3,natom,bl(mataddr('dipp')),inx,2,0,bl(ibas),
     1           bl(inuc),ncs)
ccc   call matprint('dipp',6)
      call matprodtr('DDT','dipp',dip(2))
C     call matwrite('dipp',np1,0,'aoY     ')
c
      call inton(3,natom,bl(mataddr('dipp')),inx,3,0,bl(ibas),
     1           bl(inuc),ncs)
ccc   call matprint('dipp',6)
      call matprodtr('DDT','dipp',dip(3))
C     call matwrite('dipp',np1,0,'aoZ     ')
c
      call matrem('dipp')
c
      write(6,*)
      write(6,*) 'MP2 DIPOLE MOMENTS'
c  the electronic part is negative because of the negative charge
      dip(1) = -dip(1)
      dip(2) = -dip(2)
      dip(3) = -dip(3)
c     write(6,100) dip,tot
c  calculate nuclear dipole and add to dip
      call nucdipole(natom,dip,bl(inuc))
c
      tot=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
      write(6,100) dip,tot
  100 format(5x,'dipole/au', 3f10.6,'  total',1f10.6)
      dip(1)=dip(1)/debye
      dip(2)=dip(2)/debye
      dip(3)=dip(3)/debye
      tot=sqrt(dip(1)**2+dip(2)**2+dip(3)**2)
      write(6,200) dip,tot
      write(6,*)
  200 format(5x,'dipole/D ', 3f10.6,'  total',1f10.6)
      end
C==============FixH===========================================
      subroutine FixH(Alfa,H,Hny,idim)
      use memory
      implicit real*8(a-h,o-z)
      dimension alfa(*), H(idim,*),Hny(idim,*)
      if(idim.gt.1)then
      do i=1,idim
      hik=alfa(idim)*H(i,idim)
      do j=1,idim-1
      hik=hik+alfa(j)*H(i,j)
      enddo
      Hny(i,idim)=hik
      Hny(idim,i)=hik
      enddo
      endif
C  do Hkk
      hkk=alfa(idim)*alfa(idim)*H(idim,idim)
      if(idim.gt.1) then
      do i=1,idim-1
      do j=1,idim-1
      hkk=hkk+alfa(i)*alfa(j)*H(i,j)
      enddo
      hkk=hkk+2.0d0*alfa(i)*alfa(idim)*H(i,idim)
      enddo
      endif
      Hny(idim,idim)=hkk
      end
C==============FixGZ==========================================
      subroutine FixGZ(alfa,ndisk,idim,nvir,nmo)
      use memory
      implicit real*8(a-h,o-z)
      dimension alfa(*)
      nstrip=nvir*nmo
      call matdef('Gzk','r',nvir,nmo)
      call matdef('Gzkp','r',nvir,nmo)
      igza=mataddr('Gzk')
      igzap=mataddr('Gzkp')
      call readz(ndisk,idim,bl(igza),nstrip)
      call matscal('Gzk',alfa(idim))
      if(idim.gt.1) then
      do i=1,idim-1
      call readz(ndisk,i,bl(igzap),nstrip)
      call matadd1('Gzkp',alfa(i),'Gzk')
      enddo
      endif
      call writez(ndisk,idim,bl(igza),nstrip)
      call matrem('Gzkp')
      call matrem('Gzk')
      end
c==================================kw=========================

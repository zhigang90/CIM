C ********************************************************
C * subroutine for calculating MP2 force of CIM clusters *
C * Modify the source code of force2_MP2 from PQS        *
C * Zhigang 3/23/2017 @UARK                              *
C ********************************************************

      subroutine force2_MP2_CIM(natom,nmo,nvirt,ncen,
     &                          rhf,ictr,iprnt,trans)
C    interface for cim-mp2-gradients
C    this subroutine is called from cim_mp2 if force is to be calculated.
C    the MP2 correction to the gradient is calculated and added to the
C    HF gradients.
C
C    this subroutine calls mp2_grad which is the main routine for
C    the MP2-gradients.
C
C    Svein Saebo, Fayetteville, AR Summer 2002
C
C    In CIM, we do not have dummy atoms, so only natom is needed here.
C    Zhigang 3/24/2017 @UARK
C
      use memory

      implicit real*8(a-h,o-z)
      logical rhf
      character*256 jobname

      common /intlim/ limxmem,limblks,limpair
      COMMON /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /job/ jobname,lenJ
c
      parameter (idft=0)

      dimension trans(nmo,nmo)
      integer iatom(natom)
      integer,allocatable::IAN(:),Z(:),ILST(:,:),INX2(:,:)
      real*8,allocatable::dtmp(:),Fock1(:,:),XC(:,:),QA(:),BASDAT(:,:)
      real*8,allocatable::tmp(:),SOVER(:,:),Ccen(:,:),Ccen1(:,:)
      character*8,allocatable::AtSymb(:)

C For test
      real*8,allocatable::Sij(:,:),SOVER2(:,:)
      integer same

C Arguments
C ---------

C iatom(natom): Corespondence between atoms in the cluster and in the
C               whole system

      write(iout,*)
 10   format(72('='))
      write(iout,10)
      write(iout,"(15x,' The Module of CIM Cluster Force Calculation')")
      write(iout,*)

C
C  currently only closed-shell MP2 gradients are available
      If(.NOT.rhf) Call nerror(1,'CIM MP2 GRADIENT module',
     $   'Sorry - Open-Shell MP2 gradients currently Unavailable',0,0)
C
C  set values appropriate for integral derivatives during MP2
      limxmem=2000000
      limblks=300
      limpair=100
C
      np4=4
      np1=1
c
c -- recover values from depository
      ncf=igetival('ncf')
      ncs=igetival('ncs')
      call getival('nsh',nsh)
      ncore=igetival('ncore')
      call getrval('thresh',thresh)
      call getmem(natom*3,lforc2)
      call setival('lforc2',lforc2)
C
C  define matrices and reserve memory for SCF-data:
C  epsi - orbital energies, cano - Canonical orbitals,
C  ovla - overlap matrix, fock - Fock matrix
C  ** density matrix defined and read in calling routine **
C  Fock matrix should be read from the XXX.cim file

      iunit=1
      open(unit=iunit,file=jobname(1:lenJ)//'.cim',form='formatted',
     &     status='old')
      LF2=ncf*(ncf+1)/2
      allocate(dtmp(LF2),Fock1(ncf,ncf))
      call rread8(iunit,'$AO-FOCK-A',LF2,dtmp(1))
      call iread8(iunit,'$ATOMS',natom,iatom(1))
      call FockFull(ncf,LF2,dtmp,Fock1)
      close(unit=iunit,status='keep')

      allocate(XC(3,NAtom),AtSymb(NAtom),IAN(NAtom),QA(NAtom))
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtom,AtSymb,XC,-1,jnk) ! In Bohr
      CLOSE(UNIT=IUnit,STATUS='KEEP')

C  get atomic numbers from atomic symbols
      CALL GetAtNo(NAtom,AtSymb,IAN)
      call GetAtChrg(NAtom,QA)

      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      allocate(Z(natom),ILST(4,ncs),BASDAT(13,nsh))
      call rdbasis(IUnit,NAtom,AtSymb,XC,Z,ILST,BASDAT)
      deallocate(Z)
      CLOSE(UNIT=IUnit,STATUS='KEEP')

      allocate(INX2(12,ncs))
      CALL SortBAS1(NAtom,ncs,ILST,INX2)
      CALL normaliz(ncs,INX2,BASDAT)

      i1=1
      i2=i1+ncf
      i3=i2+ncs
      IEnd=i3+12*ncs-1
c
      allocate(Z(IEnd))
      Call reorder(ncs,INX2,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinskiorder
      Call SortBAS2(ncs,Z(i3),Z(i2),Z(i1),INX2)        ! per atom

      call matdef('fock','s',ncf,ncf)
      lfock0=mataddr('fock')

C      allocate(SOVER2(ncf,ncf))
C      call ReorderFock2(ncf,ncf,Z(1),Fock1,SOVER2)

      call ReorderFock3(ncf,Z(1),Fock1,bl(lfock0))

C Read the coefficients of central orbitals (LMO) and reorder to KW
C order. -NZG_5/22/2017 @UARK
      allocate(Ccen1(ncf,ncen),Ccen(ncf,ncen),tmp(ncen))

      itype=1
      lenM=lenJ+4
      call ReadMOS(ncf,Ccen1,tmp,.false.,lenM,jobname(1:lenJ)//'.cen',
     &             itype,IErr)
      If (IErr.NE.0) call nerror(3,'CIM Cluster Force Module',
     &   'MOS File Does Not Exist',0,0)

      do J=1,ncen
         do I=1,ncf
            II=Z(I)
            Ccen(I,J)=Ccen1(II,J)
         enddo
      enddo

      deallocate(dtmp,Fock1,tmp,Ccen1)
           
      call matdef('ovla','s',ncf,ncf)
      lovla=mataddr('ovla')
      allocate(tmp(127008),SOVER(ncf,ncf))
      call inton2(0,natom,SOVER,INX2,INX2,0,0,BASDAT,BASDAT,XC,
     &            IAN,ncs,ncs,ncf,ncf,tmp)
       
      call ReorderFock3(ncf,Z(1),SOVER,bl(lovla))
      deallocate(tmp,Z)

C For debug
C Compare the overlap got by two ways to insure that it is right
C      call dynamic_matdef('overlap','q',ncf,ncf)
C      ioverlap=mataddr('overlap')
C      call overlapbuilder(ioverlap)
C      allocate(SOVER2(ncf,ncf))
C      call ReorderFock(ncf,ncf,Z(1),bl(ioverlap),SOVER2)
C      same=0
C      call checkiden(ncf,ncf,SOVER,SOVER2,same)
C      write(iout,*) "same?",same
c -- save pointers in depository
      lvec=mataddr('cano')
      lval=mataddr('epsi')
      call setival('lfock0',lfock0)
      call setival('lvec',lvec)
      call setival('lval',lval)
      nocc=nmo
      call setival('nocc',nocc)

C calculate density matrix from the MO coefficients
C NZG_5/4/2017 @UARK
      call matdef('den0','s',ncf,ncf)
      lden = mataddr('den0')
      call setival('ldensi',lden)
      
      call NJ_denmat_sym(ncf,ncf,nocc,bl(lvec),bl(lden))

C calculate density matrix from only the central MOs
C NZG_6/1/2017 @UARK
      call matdef('dencen','s',ncf,ncf)
      ldencen=mataddr('dencen')
      call setival('ldencen',ldencen)

      call NJ_denmat_sym(ncf,ncf,ncen,Ccen,bl(ldencen))

c  check if the basis set contains L-shells
c  and if it does then make S,P partitioning
c
      call trans_l_2_sp(rhf,bl,bl(ictr),ictr,lshell,'mp2grad')
c
c ictr is changed on return if l-shells are present
c
      if(ncore.gt.0) call matsub('coro','cano',1,ncore)
C      call matsub('occa','cano',1,nmo)
C      call matsub('occu','cano',ncore+1,nmo)
C      call matsub('virt','cano',nmo+1,ncf)
      nval=nmo-ncore
c
      if(IPRNT.GT.0) write(6,"(' Integral Threshold:',E10.4)") thresh
c-----------------------------------------------------------------
c allocate memory for an array showing if a big bs shell 
c was present (0 or 1) in a small bs
c
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c
c this is needed for dual basis set but is also used
c in blocking for ordinary MP2 (must have all 1)
c-----------------------------------------------------------
      call set_mpres(bl(mpres_in_bg),ncs)
c-----------------------------------------------------------------
      call mmark

C zero-out the two-electron contributions to the force
      call zeroit(bl(lforc2),natom*3)

      call mp2_grad_cim(ncf,nval,nvirt,IPRNT,thresh,nmo,natom,
     &                  bl(lforc2),ncore,trans,Ccen,ncen,iatom)
c-----------------------------------------------------------------
c
c  At this point two-electron contributions to the gradient are known.
C  Print it out and it will be needed to sum up for the whole system.
c  NZG_5/4/2017 @UARK
c
54    format(72('*'))
      write(iout,*) 
      write(iout,54)
      write(iout,*) "Two-electron contributions of the forces:"
      inuc=igetival('inuc')
      call torque_CIM(Natom,0,bl(inuc),bl(lforc2))
      write(iout,54)

      call retmark
c
c transfer back original L-shell basis set info
c
      call trans_sp_2_l(bl,bl(ictr),ictr,lshell)
c
c returns the original value of ictr
c
C      call matreset
c     call retimark
      call retmark
c ...........................................................
      return
      end


C ====================================================================== 
      subroutine mp2_grad_cim(ncf,    nval,   nvir,   iprint, thresh,
     1                        nmo,    natoms, gradv,  ncore,  trans,
     2                        Ccen,   ncen,   iatom)
C
C  Main routine for CIM-MP2-gradients. 
C  In CIM we don't release memory after energy calculation. Some of the
C  matrices do not need to be defined twice. So I modified part of the
C  subroutine. - Zhigang 3/31/2017 @UARK

C  The residia (Tij) in (virtual) MO basis are on disk, opened via unit
C  ndisk1, stored as 5-byte integers.
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
      dimension gradv(3,*),trans(nmo,nmo),Ccen(ncf,ncen),iatom(natoms)
      character*256 scrfile,filname1,filname2,filname3,filname4,filname5
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
  100    format(1x,1f8.2,' minutes ',1f8.2,' minutes ')
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
         call matprint('Gmat',6)
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
  112    format(' Memory requirements for reversed Binsort:',/,
     2   ' length of one bin:       ',I10,' bytes',/,
     3   ' number of bins:          ',I10,/,
     4   ' memory for bins:         ',I10,' bytes',/,
     5   ' number of indexpairs:    ',I10)
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
      call D1terms_CIM(ndisk3,bl(iypadr),nval,   iprint, thresh,
     1                 ncore,  bl(ibin), bl(i1), lbin,   nrcpb,
     2                 ncf, nmo,iscs)
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
      call Rphas2_CIM(ncf,nval,ndisk3,nbins,lbin,thresh,bl(ibin),bl(i1),
     1                bl(itmo),bl(itao),nsym,iprint,bl(ictr),nrcpb,
     2                bl(ifp),gradv,iscs,natoms,trans,ncen,Ccen)
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
C NZG
C      if(iprint.ge.0) then
         Write(iout,*) ' MP2 gradients after A2-terms'
         call torque_CIM(NAtoms,0,bl(inuc),gradv,iatom)
C      endif


      stop
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
      end subroutine mp2_grad_cim


C ======================================================================
      subroutine D1terms_CIM(ndisk3, Y,      nval,   iprint, thresh,
     1                       ncore,  ibins,  int1,   lbin,   nrcpb,
     2                       ncf,    nmo,    iscs)

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
C  Called once from mp2_grad_cim
C 
C  Modified by Zhigang 3/31/2017 @UARK For some dsmx matrix definition
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
C      call matdef('dsmx','q',ncs,ncs)   ! filled in <DmxMakeC>/<DmxMakeL>
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
            if (smal) then
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
     &                              kcf2,bl(ictr),bl(lmp2int),bl(ioccu),
     &                              iprint,thresh,bl(ixadr),nrow,ncol,
     &                              bl(irow),bl(icol),bl(irow1),
     &                              bl(icol1),smal,Y,bl(itmylam),ibins,
     &                              int1,lbin,nrcpb,ndisk3,iscs)
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
               call int_lmp2b(bl,bl(ictr),thresh,     ics,     kcs,
     1             bl(mapf2s),bl(idics),iprint,  bl(lmp2int),nintotal,
     2             bl(icol), bl(irow),  bl(if2cc),bl(if2cr),bl(indlj),
     3             bl(lrestor),  nrow,    ncol,     ind,     intstore,
     4             lmp2_size, lenmax)
               call retmark
               ttint=ttint+nintotal
               if(nintotal.eq.0) then
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
               call Trasigtoj(ncf,    ncs,    nval,
     *                    ics,icf1,icf2,kcs,kcf1,kcf2,
     1              bl(ictr),bl(lmp2int),bl(ioccu),iprint,thresh,
     2              bl(ixadr),nrow,   ncol,   bl(irow),bl(icol),
     3              bl(irow1),bl(icol1),smal, Y,    bl(itmylam),
     4                ibins,  int1,   lbin,   nrcpb,  ndisk3,
     5                iscs)
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


C ======================================================================
      subroutine torque_CIM(natoms,idft,datnuc,forc3,iatomf)

      use memory

      implicit real*8 (a-h,o-z)
c
c  This routine calculates the sum and torque and the max.
c  Cartesian component of the forces & prints them out
c
      common /forcint/ ncache,nprint,maxprice
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /big/bl(10000)
      dimension datnuc(5,*)
      dimension forc3(3,natoms),iatomf(natoms)
      parameter (zero=0.0d0,thrsh=1.0d-6)
c---------------------------------------------------------------------
c INPUT :
c
c natoms    - number of atoms
c idft      - dft flag  0 - no dft
c                      >0 - dft contribution to forces needs to be computed
c forc3     -  space for total forces (will be read in from file)
c iatomf(i) - label of ith atom in the whole system
C
c OUTPUT:
c
c forc3 - total forces on atoms
c
c---------------------------------------------------------------------
C Simply modify the original subroutine in PQS for always printing out
C the energy. -NZG_5/5/2017 @0UARK
      write(iout,151)
  151 format(
     *' Atom  Name FullLabel   force-x      force-y     force-z'/)
c
      sumx=zero
      sumy=zero
      sumz=zero
      torqx=zero
      torqy=zero
      torqz=zero
      do iat=1,natoms
         xat=datnuc(2,iat)
         yat=datnuc(3,iat)
         zat=datnuc(4,iat)
         sumx=sumx+forc3(1,iat)
         sumy=sumy+forc3(2,iat)
         sumz=sumz+forc3(3,iat)
c
c torque= RxF (cross product of the coordinate and the the force
         torqx=torqx+forc3(3,iat)*yat-forc3(2,iat)*zat
         torqy=torqy+forc3(1,iat)*zat-forc3(3,iat)*xat
         torqz=torqz+forc3(2,iat)*xat-forc3(1,iat)*yat
c
         write(iout,200) iat,datnuc(5,iat),iatomf(iat),
     *                   forc3(1,iat),forc3(2,iat),forc3(3,iat)
c
      enddo
c
  200 format(i3,5x,a3,3x,i4,4x,3(f11.7,1x))
c
      formax=abs(max(sumx,sumy,sumz,torqx,torqy,torqz))
      formin=abs(min(sumx,sumy,sumz,torqx,torqy,torqz))
      forabs=max(formin,formax)
c
      if((forabs.gt.thrsh.and.(idft.eq.0.or.nprint.gt.1)).or.
     1   forabs.gt.thrsh*50.0d0) then
         write(iout ,220) sumx,sumy,sumz,torqx,torqy,torqz
  220    format('Sum or total torque of the forces does not vanish',/,
     *   'Sum=',3f10.7,' Torque=',3f10.7)
      endif
       call absmax(3*natoms,forc3,ifmax,fmax)
ckw   ifmax=(ifmax-2)/3+1
ckw   write(iout,230) fmax,ifmax
      iatom=(ifmax-2)/3+1
      iatom_x=(iatom-1)*3+ 1
      iatom_y=iatom_x+1
      iatom_z=iatom_y+1
      if(ifmax.eq.iatom_x.or.ifmax.eq.iatom_y.or.ifmax.eq.iatom_z) then
         iatomx=iatom
      else
         iatomx=iatom+1
      endif
      write(iout,230) fmax,iatomx
  230 format(/,' Maximum Cartesian force component=',f12.7,
     1  ' Eh/a0 on atom ',i4)
c
      end subroutine torque_CIM


C========Rphas2_CIM===================================================
      subroutine Rphas2_CIM(ncf,nval,ndisk3,nbins,lbin,thresh,ibin,int1,
     &                      Tmnmo,Tmnao,nsym,iprint,inx,nrcpb,ifunpair,
     &                      gradv,iscs,natom,trans,ncen,Ccen)

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
c     NZG_5/16/2017 @UARK
c     Modify the origin code for CIM force calculation. In the CIM force
c     calculation, we need to transform the T from QCMO to LMO.

      common /lindvec/ lind,idensp
      integer*1 int1(2,lbin)
      integer*4 ibin(2,lbin)
      character*11 scftype
      dimension Tmnao(ncf,ncf),Tmnmo(nval,*),trans(nval,*),Ccen(ncf,*)
      dimension ifunpair(7,*)
      dimension gradv(3,*)
      dimension inx(12,*)
      dimension xintxx(9)
      parameter(zero=0.0d0,half=0.5d0,two=2.0d0,four=4.0d0)
      parameter (dblmax=2147483648.0d0)
C
      real*8,allocatable::Tlq(:,:)
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

C      call matcopy('den0','denf')
      call matcopy('dencen','denf')
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
C  For large systems, the subroutine may cause large errors. -NZG
C  In CIM calculation, we don't do the prescreening. -NZG 6/1/2017
C            call TranSetup(bl(icta),ncf,nval,ipr,jpr,
C     1                     bl(ipoint),bl(jpoint),bl(ischwarz),ncs,MYS,
C     2                     LAS,ancf,bl(mapf2s),thresh)
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

C  First, transform one of the occupied indices from QCMO to LMO
                  allocate(Tlq(ncen,nval))
                  call TQtoL(Tmnmo,Tlq,trans,ncen,nval)
                  call BackTrans_CIM(Tlq,Tmnao,ncf,nval,ncen,Ccen,
     &                               bl(iocca))
                  deallocate(Tlq)
                  call secund(ttbt2)
                  ttbt=ttbt+ttbt2-ttbt1
CC
CC   SS July 2003
CC   add extra terms to TAO here to avoid repeated integral derivatives
C
C                  call matzero('denf')
                  call matzero('ddtf')

                  call addtoT1_CIM(tmnao,bl(idena),bl(iddt),ncf,my,lam)
                  call moveTsh(Tmnao,bl(iTadr),ncfsq,my3,lam3,
     1                         MYS_size,LAS_size)
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
C ====================================================================== 

     
C ======================================================================
      subroutine TQtoL(Tmnmo,Tlq,trans,ncen,nval)
C this routine transforms one of the indices of Tmnmo index from QCMO to
C LMO.
C NZG_5/16/2017 @UARK

      use memory
      implicit none

      integer ncen,nval
      real*8 Tmnmo(nval,nval),Tlq(ncen,nval),trans(nval,nval)

      Tlq=0.0D0
      call dgemm('T','N',ncen,nval,nval,1.0D0,trans(:,1:ncen),
     &           nval,Tmnmo,nval,0.0D0,Tlq,ncen)

      end subroutine TQtoL     
    

C =====================================================================
      subroutine BackTrans_CIM(Tlq,Tmnao,ncf,nval,ncen,Ccen,Cqcmo)
 
      use memory
      implicit none

      integer ncf,nval,ncen
      real*8 Tlq(ncen,nval),Tmnao(ncf,ncf)
      real*8 Ccen(ncf,ncen),Cqcmo(ncf,nval),CqcmoT(nval,ncf)
      real*8 temp(ncen,ncf)

      CqcmoT=transpose(Cqcmo)
      temp=0.0D0
      call matmul_mkl(Tlq,CqcmoT,temp,ncen,nval,ncf)
      Tmnao=0.0D0
      call matmul_mkl(Ccen,temp,Tmnao,ncf,ncen,ncf)

      end subroutine BackTrans_CIM


      subroutine addtoT1_CIM(T,D,DDT,ncf,my,lam)
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
      dimension T(ncf,*),D(ncf,ncf),DDT(ncf,ncf)
      parameter(half=0.5d0,one4=0.25d0,one8=0.125d0,on16=0.0625d0)

C NZG
      DDT=0.0D0


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
      

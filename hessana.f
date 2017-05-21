c======================================================================
c                THE MAIN ANALYTICAL HESSIAN PROGRAM
c
c Note : because of re-scaling of two-el. integrals & int. deriv.
c        the two-electron part of the calculation is done FIRST
c======================================================================
c PP "Ab Initio Methods in Quantum Chemistry -II, Edited by K.P. Lawley
c 1987 John Wiley & Sons Ltd. page 241
c
c
c Second derivative of energy is given by :
c
c (1) asymmetrical formula :
c
c    ab            ab    ab      ab                        ab
c   E   = 0.5 Tr [h  +  h + G(D,g  )]*D   - 2 Tr (C e C+)*S
c
c                        a       a     b
c        +1.0 Tr [      h + G(D,g  )]*D    !!! missing in PP paper !!!
c
c
c                 a   b            b             b
c          -2 Tr S ( C e C+  +  C e  C+  +  C e C+ )
c                      1st-order weighted density
c
c             ab
c          + Vnuc
c
c   where D=2*CC+   and g_ijkl=(ij|kl) - 0.5 (il|kj)
c
c
c
c The term involving 1st-order weighted density can be expressed in
c terms of the 1st-order ordinary density D1=C1C+ + CC1+ .
c For this one may use
c                       b    b         b           b
c     e= C+F C   and   e  = C+F C + C+F C  +  C+F C
c
c  Then
c
c                  a   b            b             b
c           -2 Tr S ( C e C+  +  C e  C+  +  C e C+ ) =
c
c                  a  b                 a    b
c             -Tr S  D  F D   - 0.5 Tr S  D F  D .
c
c                      a    b           b            b
c or          -0.5 Tr S  [ D  F D +  D F  D  +  D F D  ]
c
c                            b     b        b          b
c  where F = h + G(D,g)  ;  F  =  h   +  G(D,g) + G(D,g )    !!!!
c
c======================================================================
c                    ab
c Final formula for E     derivatives . Terms in order of calculations :
c
c
c
c    ab                 ab
c   E   = 0.5*Tr [ G(D,g  )]*D     term 1  calculated in hess2e.f
c
c             ab
c          + Vnuc                  term 2  calculated in hess1.f
c
c                  ab
c         1.0*Tr [h  ]*D           term 3  calculated in hess1.f
c
c                 ab
c       - 2.0*Tr S *[C e C+]       term 4  calculated in hess1.f
c
c
c  and after solving the CPHF equations
c
c                   a       a     b
c        +1.0*Tr [ h + G(D,g ) ]*D    term 5  calculated in hess_total
c
c                 a    b           b            b
c        -0.5 Tr S  [ D  F D +  D F  D  +  D F D  ]
c
c======================================================================
      subroutine hessana

      use memory

      implicit real*8 (a-h,o-z)
      character jobname*256,cdum*20,GROUP*4
      Logical rhf,Symflag,AtomAx
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,
     1              ncs,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      dimension RM(3,3)
c-------------------------------------------
c for several calls of the hessian program in one run:
      data ncall /0/
      save ncall
c----------------------------------------------------------------------
      PARAMETER (IUnit=1)               ! unit number for checkpoint I/O
c
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
c
      DATA thrsh/1.0d-6/                ! zero threshold
c
c-----------------------------------------------------------------------
      ncall=ncall+1
c-----------------------------------------------------------------------
c get original values of symmetry pair pointers
c
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c-----------------------------------------------------------------------
c Re-store commons big, intbl, and ener
c
      call readbl(ncall,ictr)
c
c output : ictr = new address in inx (texas95)
c-----------------------------------------------------------------------
c  check for COSMO (cannot do analytical Hessian if COSMO on)
      call tstival('cosmo',icosmo)
      If(icosmo.NE.0) Then
        call getival('cosmo',icosmo)
        If(icosmo.NE.0) call nerror(13,'Analytical Hessian Module',
     $   'Cannot do analytical Hessian with COSMO Solvation Model',0,0)
      EndIf
c-----------------------------------------------------------------------
c  read from <control> file
c  total number of atomic centers
c  multiplicity
c  number of alpha electrons
c  number of beta electrons
c  dft flag
c  number of dummy atoms (if any)
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $     FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      CALL RdDFT(IUnit,idft)
c
c
c -- dummy atoms
c
      ndum1=0     ! number of charged dummies
c
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
       NAtom = NAtoms-Ndum1-Ndum2
      Else
       NAtom = NAtoms
      EndIf
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c -- make sure value for real atoms only is in the depository
c
      call setival('ndum1',ndum1)
      call set_Natoms(natom)
c
c -- save number of point charges for use in hess1
c
      call setival('ndum1',ndum1)
c
c -- check there are electrons!
      If(NAlpha.EQ.0) Call nerror(1,'Analytical Hessian Module',
     $                     'There Are No Electrons!',0,0)
c
c -- set logical flag "rhf"
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
c-----------------------------------------------------------------------
c Read in the hessian options
c
      call hess_options(idft,rhf,fact0)
c
c-----------------------------------------------------------------------
c Check if any change in grid factor
      If(idft.GT.0.AND.fact0.NE.0) Then
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call rdcntrl(IUnit,7,'$factor',2,idum,factor,cdum)
c -- write new grid factor on the control file
        call wrcntrl(IUnit,7,'$factor',2,idum,fact0,cdum)
c -- if new grid factor LESS than used in SCF print warning
        if(fact0.LT.factor) call message('HESSIAN',
     $        '**WARNING** DFT grid quality REDUCED from SCF',0,0)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
      EndIf
c---------------------------------------------------------------
c allocate memory for the final hessian matrix hess(3*nat,3*nat)
c
      nat3=3*natom
      call getmem(nat3**2,lhess)
      call setival('lhess',lhess)
c-----------------------------------------------------------------
c zero out the hessian matrix:
c
      call zeroit(bl(lhess),nat3**2)
c---------------------------------------------------------------
c in the lhess location the following will be kept :
c
c 1. + tr D0*G(D0,gxy) = Sum(ijkl) Dij*dkl*gxy,ijkl (2-el)
c 2. nuclear repulsion contributions
c 3. + tr D0*Hxy                                    (1-el)
c 4. + tr W0*Sxy                                    (1-el)
c
c   and after CPHF is solved
c
c 5. + tr Dy*F(D0,gx)  where F(D0,gx)= Hx + G(D0,gx) (Fock derivative)
c 6. + Tr Wy*Sx
c
c These above are ALL terms for the final hessian matrix
c---------------------------------------------------------------
c put down a memory mark
      call mmark
c     call immark
c---------------------------------------------------------------
c read in the density & fock matrix
c as well as eigenvectors & eigenvalues
c make weighted density
c
      ntri=(ncf*(ncf+1))/2
c
      call getmem(ntri,lsmat)
      call setival('lsmat',lsmat)
c
      call getmem(ntri,lden)
      call setival('ldensi',lden)
      if(.not.rhf)then
        call getmem(ntri,mden)
        call setival('ldensB',mden)
      endif
      call getmem(ntri,lfoc)
      call setival('lfock0',lfoc)
      if(.not.rhf)then
        call getmem(ntri,mfoc)
        call setival('lfock0B',mfoc)
      endif
c
      call getmem(ntri,ldew)         ! weighted density
      call setival('ldewsi',ldew)
c
      call getmem(ncf*ncf,lvec)      ! eigenvectors
      call getmem(ncf    ,lval)      ! eigenvalues
      call setival('lvec',lvec)
      call setival('lval',lval)
      if(.not.rhf)then
        call getmem(ncf*ncf,mvec)      ! eigenvectors
        call getmem(ncf    ,mval)      ! eigenvalues
        call setival('lvecB',mvec)
        call setival('lvalB',mval)
      endif
c
      call read_dfeig(rhf,ncf,NAlpha,NBeta)
c
      call setival('nocc',NAlpha)              ! put in depository
      if(.not.rhf)call setival('noccB',NBeta)  ! put in depository
c---------------------------------------------------------------
c check if the basis set contains  L-shells if so, make S,P partitioning
c check if the basis set contains GC-shells and what type
c
      call getival('ncs ',ncs)
      call check4lgc(ncs,bl(ictr),lshell,lgcshell,ltype_gc,ldeep_gc)
c
c
      if(lshell.gt.0) then
         call trans_l_2_sp(rhf,bl,bl(ictr),ictr,lshell,'hessian')
      endif
c
      if(lgcshell.gt.0.and.ltype_gc.le.2.and.ldeep_gc.le.2) then
c
c       if there is a GC basis set but only as deep as 2 and only
c       s-functions are GC then we DO NOT profit from GC, esspecially
c       in the gradient and probably hessian. Thus, segment it
c
         call trans_gc_2_sg(rhf,bl,bl(ictr),ictr,nogcsi,'hessian')
c        write(iout,*)
c        write(iout,*)
c    * '   General Contracted shells have been segmented for 2e-hess'
c              write(iout,*)
c    * '   because there were only s/3 & p/2 gen.contracted orbitals'
c              write(iout,*) ' '
      endif
c
      call f_lush(6)
c
c ictr is changed on return if l-shells or gc-shells are present
c----------------------------------------------------------------------
c open files for: screen                        60
c                 fock1                         61
c                 overlap1                      62
c                 density1                      63
c                 weighted density1             64
c                 constant part of density1     65
c                 part of the hessain (save)    66
c
c  for the uhf case we have also:
c
c                 fock1 beta                    71
c                 density1 beta                 73
c              constant part of density1 beta   75
c
      call open4hess(ncs,ncf,ntri,natom,rhf)
c
c--------------------------------------------------------------
c -- set up common /negxyz/ symmetry arrays here
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('nsyo',nsyo)
         call make_negxyz(nsym,bl(nsyo))
      endif
c---------------------------------------------------------------
c THE TWO-ELECTRON PART OF THE HESSIAN:
c -------------------------------------
c Hartree-Fock or Coulomb DFT
c
c---------------------------------------------------------------
c First do the 2-el. part of the Hessian
c because of re-scaling of 2-el. integ.deriv.
c---------------------------------------------------------------
c (1)   direct contributions to the final hessian
c       Tr { D0*G(D0,gxy)} = Sum ijkl {Dij*Dkl*gxy,ijkl}
c       stored in lhess .
c       This is the very first contribution there and
c       it is rescaled right away by integ. threshold
c
c -- find symmetry unique atoms :
c
      call find_unqat(bl,natom)
c
c----------------------------------------------------------------------
      call getival('resta',irestart)
      if(irestart.gt.0) then
         nfile68=68
         call read1mat(nfile68, 1 ,nat3*nat3,bl(lhess))
         go to 4321
      endif
c----------------------------------------------------------------------
      call mmark
c     call immark
c
      call hess2e2(idft,ax,rhf,bl,bl(ictr))     ! TrD0*G(D0,g2)
      call f_lush(6)
c
      call retmark
c     call retimark
c
      call getival('printh',nprint)
      if(nprint.ge.1) then
         call hess_print(bl(lhess),natom,1.d0,'after hess2e2   ')
      endif
c----------------------------------------------------------------------
c     call memory_status(' after hess2e2')
c---------------------------------------------------------------
c for possible 2-particle density type screening make
c the second screening density as a maximum of
c     Ci*Ca(T)/ei-ea
c over all occupied (i) and virtual (a) orbitals.
c
CCCC  call getival('nocc',nocc)
CCCC  call make_d2hesx(bl,inx(ictr),nocc,ncf,bl(lvec),bl(lval))
      call make_d2hess(rhf,bl,bl(ictr),ncf,ncs)
c---------------------------------------------------------------
c (2)   second part of the derivative fock matrices
c       G(D0,gx) stored in lfock1 . This is the very
c       first contrib. there and it is rescaled right
c       away by the integ. thresh.
c       (later the 1-el Hx will be added to it forming
c       the final fock1 (derivative fock) matrix .
c---------------------------------------------------------------
      call mmark
c     call immark
c
      call hess2e1(idft,ax,rhf,bl,bl(ictr))  ! Fder(D0,g1) for natoms
c
c     call getival('printh',nprint)
c     if(nprint.ge.1) then
c        call hess_print(bl(lhess),natom,1.d0,'after hess2e1   ')
c     endif
c
      call retmark
c     call retimark
c---------------------------------------------------------------
c     call memory_status(' after hess2e1')
c---------------------------------------------------------------
c  DFT CONTRIBUTIONS
c
c  Calculate direct DFT XC contributions to final Hessian matrix
c  and to derivative Fock matrices here
c...............................................................
c
      IF(idft.GT.0) THEN
        call getmem(0,lastx)
        call retmem(1)
        CALL HESSDFT(bl(lhess),lcore-lastx, bl(lastx))
      ENDIF
c---------------------------------------------------------------
c     if(nprint.ge.1) then
c        call hess_print(bl(lhess),natom,1.d0,'after hessDft   ')
c     endif
c----------------------------------------------------------------------
c *************************************
c THE ONE-ELECTRON PART OF THE HESSIAN:
c *************************************
c---------------------------------------------------------------
c calculates :
c
c 1. nuclear repulsion contributions to the final hessian
c    stored in lhess
c    (added to Tr { D0*G(D0,gxy)} which is already there)
c
c 2. direct contributions to the final hessian from
c    Tr{D0*Hxy} and Tr{W0*Sxy}
c    stored in lhess
c    (add to final hess, which already contains
c                     tr(d0*g(d0,g2)) + nuc.rep.)
c
c 3. Matrices Hx and Sx needed for CPHF
c    Hx is a 1-e part of F(D0,gx)=Hx + G(D0,gx)
c    Sx is also needed for CPHF
c    stored in lfock1  and lover1
c
c    Both Hx or precisely F(D0,gx) and Sx contribute
c    later to the final hessian (x,y) as
c
c     +Tr{ Dy*F(D0,gx) } and -Tr{ Wy*Sx }
c
c 4. D0*Mx contribution to the dipole moment derivatives
c    (Mx is the dipole moment integral derivatives)
c---------------------------------------------------------------
      call mmark
c     call immark
c
c
      call hess1(rhf,bl,bl(ictr))
c
      call retmark
c     call retimark
c--------------------------------------------------------
c     call memory_status(' after hess1  ')
c--------------------------------------------------------
c first-order Fock and overlap matrices written on a disk
c in hess1 in trasposed form :
c from fock1(3,natoms,ntri) into fock1t(ntri, 3,natoms)
c--------------------------------------------------------
      call getival('savef',isavefil)
      if(isavefil.gt.0) then
         nfile68=68
         call save1mat(nfile68, 1 ,nat3*nat3,bl(lhess))
      endif
c--------------------------------------------------------
c
 4321 CONTINUE       !   if(irestart.gt.0) go to 4321
c
c--------------------------------------------------------
c THE CPHF PART:
c --------------
c
c
c Use of Symmetry
c ---------------
c  We need only solve CPHF equations for symmetry-unique atoms.
c  Full (non-abelian) point group symmetry can be used here.
c  Abelian symmetry can be used when computing the 2-el integrals.
c  For DFT, the full grid is generated on each symmetry-unique atom.
c
c - read from FULL symmetry file number of symmetry operations
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,7,'$ntrans',1,NTrans,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c -- allocate memory for symmetry operations and equivalent atoms
c
      call getmem(natom,iuq)
      call getmem(NTrans*9,itn)
      call getmem(NTrans*natom,inq)
c
c read all data from <sym> file
c
      Symflag = NTrans.GT.1
      CALL RdSYM(Symflag,natom,  RM,     GROUP,  NTrans,
     $           NDEG,   NQ,     bl(iuq),bl(itn),bl(inq))
c
c---------------------------------------------------------------
c
      call mmark
c     call immark
c
c
c solve the Coupled Perturbed Hartree-Fock equations
c for nuclear displacement perturbations
c
c calculate the first-order density and weighted density matrices
c
      call get_d1w1(idft,  ax,  rhf,    bl,    bl(ictr),
     $              natom, NQ,  bl(iuq))
c
      call retmark
c     call retimark
c
c     call memory_status(' after  cphf  ')
c
c-------------------------------------------
c  FINAL HESSIAN CONSTRUCTION:
c ----------------------------
c
c calculate two remaining parts of the hessian
c    1. Tr Dy*F(D0,gx)
c    2. Tr Wy*Sx
c and add them to the hess matrix
c
      If(NTrans.EQ.1) Then
c
        call hess_total(rhf,natom,ncf,bl)
c
      Else
c
c allocate memory
c
        call mmark
        call getmem(3*ntri,lden1)
        call getmem(3*ntri,lwen1)
        call getmem(3*ntri,lfock1)
        call getmem(3*ntri,lover1)
c
        CALL Hess_Tot(rhf,natom,  NQ,    bl(iuq), ncf,    ntri,
     $               bl(lden1),bl(lwen1),bl(lfock1),bl(lover1),
     $               bl(lhess))
        call retmark
c
c construct and symmetrize final Hessian
c
        call getmem(natom,iss)
        call mmark
        call getmem(nat3**2,lhold)
        call getmem(nat3**2,lpp)
        call getmem(nat3**2,lrr)
c
c -- set up ISYM array
        CALL GetISYM(natom,NQ,bl(iuq),bl(iss))
c
        CALL FulHES(natom,  NTrans, bl(iss),bl(inq),bl(itn),
     $              bl(lpp),bl(lrr),bl(lhess))
c
        CALL CpyVEC(nat3**2,bl(lhess),bl(lhold))
        CALL SymHES(natom,  NTrans, bl(inq),bl(itn),bl(lhold),
     $              bl(lpp),bl(lrr),thrsh,  bl(lhess))
        call retmark
        call tstival('vcd',ivcd)
        if (ivcd.ne.0) then
            call getmem(natom,lidon)  
            call getmem(9*natom,laai)  
            call getmem(9*natom,laao)  
            call zeroit(bl(laao),9*natom)
            call RdAAT(3*natom,bl(laai),AtomAx)
C
C  The atomic axial tensors do NOT transform under symmetry operations in the
C  same way as normal coordinate vectors. Prepare symmetry operators for
C  axial pseudovectors
C
            call getmem(9*NTrans,itnA)
            call AssignSYM(NTrans, bl(itn), bl(itnA))
            call symaat(natom,NTrans,bl(iss),bl(inq),bl(itn),
     $                  bl(itnA),bl(laai),bl(laao),bl(lidon))
C            
C  zero out elements below thrsh
C
            do I=1,9*natom
            If(Abs(bl(laao+I-1)).LT.thrsh) bl(laao+I-1) = 0.0d0
            enddo
          call retmem(5)
        endif
      EndIf
c
c-------------------------------------------
c write the final Hessian matrix to disk
c
      call wrhess(jobname(1:lenJ)//'.hess',lenj+5,0,nat3,bl(lhess))
c
c print out the final hessian matrix if requested
c
      call getival('printh',nprint)
      if(nprint.ge.1) then
      call hess_print(bl(lhess),natom,1.d0,'the final matrix')
      endif
c-------------------------------------------
c write to the <control> file an info about hessian quality
c
c   'hessquality'  : ihessq=-1 (negative) poor, crude obtained from geom.opt
c   'hessquality'  : ihessq=+1 (positive) good, obtained from Hessian calc.
c
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
c
      ihessq=+1
      call  wrcntrl(IUnit,9,'$hessqual',1,ihessq,Rdum,Cdum)
c
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c-------------------------------------------
c Calculate the D1*M contribution to the dipole moment derivatives
c and form the final dipole moment derivatives matrix
c
      Call GetIVal('ibas',ibas)
      Call GetIVal('inuc',inuc)
c
c
c----------------------------------------------------------------------
c
      Call mmark
c
c---- Allocate memory for the 1-st order density matrix (for one atom)
      lDens1=3*ntri
      Call getmem(lDens1,iDens1)
      if(.not.rhf) Call getmem(lDens1,iDens1B)
c---- Allocate memory for dipole moment integrals
      lDipIn=3*ntri
      Call getmem(lDipIn,iDipIn)
c---- Allocate memory for dipole moment derivatives
      lDMDer=3*3*natom
      Call getmem(lDMDer,iDMDer)
c---- Allocate memory for buffer with D0*Mx contribution
      lBuff=3*3*natom
      Call getmem(lBuff,iBuff)
c
c---- Calculate IR intensities
      Call GetIntens(rhf,natom,ncf,ntri,ncs,na,bl(ictr),
     &     bl(iDMDer),bl(iBuff),bl(iDens1),bl(iDens1B),
     $     bl(iDipIn),bl(ibas),bl(inuc))
c
      Call retmark
c
c-------------------------------------------
c end of the story for scf hessian
c-------------------------------------------
c
c close all temporary files
c
      call close4hess(rhf)
c
c----------------------------------------------------------------------
c transfer back original L-shell basis set info
c
      call trans_sp_2_l(bl,bl(ictr),ictr,lshell)
c----------------------------------------------------------------------
c transfer back original GC-shell basis set info
c
      if(lgcshell.gt.0.and.ltype_gc.le.2.and.ldeep_gc.le.2) then
         call trans_sg_2_gc(bl,bl(ictr),ictr,lgcshell)
      endif
c
c returns the original value of ictr
c----------------------------------------------------------------------
c return all memory allocated in the 1 and 2-el hessian
c except the total hessian lhess
c
      call retmark
c     call retimark
c
c restore original values of symmetry pair and atom pointers
c
      call setival('SymFunPr',ifp)
      call setival('SymShPr',ifp1)
      call setival('na',NAtoms)
c
c-------------------------------------------
      call retmem(1)  ! release allocation for lhess
c-------------------------------------------
c memory status :
c
      call memory_status('end of hessian')
c-------------------------------------------
      end
c======================================================================
      subroutine hess_options(idft,rhf,fact0)
c
c  read in hessian options
c
c  ARGUMWNTS
c
c  idft    -  dft flag (input)
c  rhf     -  logical flag for closed-shell system
c  fact0   -  DFT grid factor

      use memory

      implicit real*8 (a-h,o-z)
      logical rhf
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /big/ bl(1)
c---
c options for integral calculations (from readin):
      common /intgop/ ncachx,maxpricx,iiii(2)
c---
      common /forcdbl/ thre1,thre2,tchf
      common /forcint/ ncache,nprint,maxprice
      common /intlim/ limxmem,limblks,limpair
c ...................................................................
      PARAMETER (IUnit=1)        !  unit number for checkpoint I/O
      Character jobname*256,cdum*20
      Common /job/jobname,lenJ
c  ...................................................................
c
      parameter (nopt=16)
      character*4 word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      character*256 chopv(nopt)
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp/1,    11,     13,   11,     1,
     *            1,     0,     1,     3,    11,
     *            0,     0,     0,     0,     1,
     *            11/
c
      data word /'prin','thr1','thr2','thre','iter',
     *           'ncac','noac','rese','limi','lvsh',
     *           'rest','save','nowd','nodd','scre',
     *           'grid'/
c
c---------------------------------------------------------------------
c maxprice (ipay) CANNOT be changed :
      maxprice=maxpricx
c---------------------------------------------------------------------
c thr1 =1.0d-10   ! one-el. 2-nd derivative integral threshold
c
c thr2 =1.0d-9    ! two-el. 2-nd derivative integral threshold
c       1.0d-10   ! two-el. 1-st derivative integral threshold and
c                   final integral threshold in CPHF
c       1.0d-8    ! loose integral threshold in CPHF
c
c thre =1.d0-5    ! threshold for CPHF convergance
c
c iter = 30       ! maximum number of CPHF iterations
c
c ireset = 20     ! reset CPHF if not converged after after 20 iterations
c
c xlvsh=0.0d0     ! level shift in CPHF, added to all energy differences
c---------------------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
      write(iout,*)
     *  '                The Analytical Hessian Module '
      write(iout,*)' '
      if(idft.ne.0)then
        if(rhf)then
          write(iout,*)
     *  '                           RHF/DFT '
        else
          write(iout,*)
     *  '                           UHF/DFT '
        endif
      else if(rhf)then
         write(iout,*)
     *  '                             RHF '
      else
         write(iout,*)
     *  '                             UHF '
      endif
      write(iout,*)' '
c
      call f_lush(iout)
c-----------------------------------------------------------
      call getival('nsym',nsym)
c-----------------------------------------------------------
c
c default values :
c
      nprint=0
      thre1=1.0d-10   !  1-el. second derivative integ.threshold
c
c three 2-el integral thresholds (option 'thr2') :
c
c     threg2=1.0d-9   ! thresh for 2nd-der.2e.int.  G(D,g2)
c     thref1=1.0d-10  ! thresh for 1st-der.2e.int.  G(D,g1)
c                     ! and
c                     ! final thresh for integ. in CPHF   G(D1,g)
c     thref2=1.0d-8   ! loose thresh for integ. in CPHF   G(D1,g)
c
      threg2=1.0d-9
      thref1=1.0d-10
      thref2=1.0d-8
c
c CPHF convergence threshold
      thres=1.0d-5    ! CPHF convergence threshold
c
      maxiter=30
      ncache=ncachx
      noacce=0
      ireset=20                ! reset CPHF after iter=20
      limxmem=2 000 000
c2002 limblks=300
      limblks=0
      limpair=100
      xlvsh=0.0d0
      irestart=0
      isavefil=0
      IdWt = 1             ! weight derivatives in DFT
c2003
      idelt=1              ! use delta density in CPHF
      iscre=1              ! screening type in cphf 1=1-particle-density
c                          ! screening type in cphf 2=2-particle-density
c
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
c
c ..................................................................
      If(ifound(1).gt.0) Then
        nprint=iopv(1,1)
c -- write print flag to <control> file
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call wrcntrl(IUnit,6,'$print',1,nprint,rdum,cdum)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
      EndIf
c ..................................................................
      if(ifound(2).gt.0) thre1=10.d0**(-ropv(1,2))
c
      if(ifound(3).gt.0) threg2=10.d0**(-ropv(1,3))
      if(ifound(3).gt.0) thref1=10.d0**(-ropv(2,3))
      if(ifound(3).gt.0) thref2=10.d0**(-ropv(3,3))
c
      if(ifound(4).gt.0) thres=10.d0**(-ropv(1,4))
      if(ifound(5).gt.0) maxiter =iopv(1,5)
      if(ifound(6).gt.0) ncache  =iopv(1,6)
c
      if(ifound(7).gt.0) noacce  = 1  ! no acceler.
      if(ifound(8).gt.0) ireset  =iopv(1,8)

      if(ifound(9).gt.0) then
         limxmem=iopv(1,9)
         limblks=iopv(2,9)
         limpair=iopv(3,9)
      endif
cNOT  if(ifound(10).gt.0) xlvsh=ropv(1,10)
      if(ifound(10).gt.0) then
         xlvsh=0.0d0
      endif
      if(ifound(11).gt.0) irestart= 1  ! restart calc. from CPHF
      if(ifound(12).gt.0) isavefil= 1  ! save files for next restart
      if(ifound(13).gt.0) IdWt = 0     ! switch off weight derivatives
c
      if(ifound(14).gt.0) idelt=0      ! switch off delta density in cphf
      if(ifound(15).gt.0) iscre=iopv(1,15)
c
c grid quality factor for DFT
c
      fact0=0.0d0                      ! if zero, no change
      if(ifound(16).gt.0) fact0=ropv(1,16)
c
      call setival('printh',nprint)
      call setrval('thr1',thre1)
c
      if(ireset.gt.30) ireset=30
      call setival('reset',ireset)
      call setival('maxit',maxiter)
      call setival('noacc',noacce)
      call setrval('xlvsh',xlvsh)
      call setival('resta',irestart)
      call setival('savef',isavefil)
      call setival('wghtd',IdWt)
c
      call setival('delta',idelt)
      call setival('scree',iscre)
c
c---------------------------------------------------------------------
      lopt(1)=nprint
c---------------------------------------------------------------------
c
c print options :
c
      if(ifound(1).gt.0) write(iout,210) word(1),iopv(1,1)
      if(ifound(2).gt.0) write(iout,220) word(2),ropv(1,2)
      if(ifound(3).gt.0) write(iout,220) word(3),(ropv(ii,3),ii=1,3)
      if(ifound(4).gt.0) write(iout,220) word(4),ropv(1,4)
      if(ifound(5).gt.0) write(iout,210) word(5),iopv(1,5)
      if(ifound(6).gt.0) write(iout,210) word(6),iopv(1,6)
      if(ifound(7).gt.0) write(iout,210) word(7),1
      if(ifound(8).gt.0) write(iout,210) word(8),iopv(1,8)
      if(ifound(9).gt.0) write(iout,210) word(9),(iopv(ii,9),ii=1,3)
      if(ifound(10).gt.0) write(iout,230) word(10),' not active yet'
      if(ifound(11).gt.0) write(iout,210) word(11),1
      if(ifound(12).gt.0) write(iout,210) word(12),1
      if(ifound(13).gt.0) write(iout,230) word(13),
     $                                 ' No Weight Derivatives in DFT'
      if(ifound(14).gt.0) write(iout,230) word(14),
     $                                 ' No Delta Density in CPHF    '
      if(ifound(15).gt.0) then
cc      if(iscre.eq.1) write(iout,230) word(15),
cc   *                                 ' 1-particle density screening'
        if(iscre.ge.2) write(iout,230) word(15),
     *                                 ' 2-particle density screening'
      endif
      if(idft.gt.0) then
        if(ifound(16).gt.0) write(iout,220) word(16),ropv(1,16)
      endif
c
  210 format(' HESS option = ',a4,'   is on ; value = ',3(i10,2x))
  220 format(' HESS option = ',a4,'   is on ; value = ',3(f10.5,1x))
  230 format(' HESS option = ',a4,'   is on ; value =',a35)
c---------------------------------------------------------------------
c check for linear dependencies in the basis set :
c
      if(ifound(3).gt.0) then
         call check4ld3(threg2,thref1,thref2,thres,.false.)
      else
         call check4ld3(threg2,thref1,thref2,thres,.true. )
      endif
c
      call setrval('threg2',threg2)    ! int.thres for G(D0,g2)
      call setrval('thref1',thref1)    ! final int.thresh in CPHF
      call setrval('thref2',thref2)    ! loose int.thresh in CPHF
c
      call setrval('thres',thres)      ! cphf conv. threshold
c---------------------------------------------------------------------
c
      write(iout,*)' '
c---------------------------------------------------------------------
      call f_lush(iout)
c
      end
c======================================================================
c
      subroutine hess_infoprt(iout,icond)
      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /ganz/ lcore,iov,last,lflag(4),lrest(77)
      common /runtype/ scftype,where
      common /forcint/ ncache,nprint,maxprice
      common /fieldx/ xfield,xfgrad,elfiel(9)
c
c electric field and electric field gradient :
c
      ifield=xfield
      ifgrad=xfgrad
c
cc      write (iout,223)
cc      write (iout,280)
      if( ifield .eq.1) write (iout,2820) elfiel(1),elfiel(2),elfiel(3)
      if( ifgrad .eq.1) write (iout,2821) elfiel(4),elfiel(5),elfiel(6),
     *                                    elfiel(7),elfiel(8),elfiel(9)
cc      write (iout,223)
c
cc      write (icond,223)
cc      write (icond,280)
cc      if( ifield .eq.1) write (icond,2820) elfiel(1),elfiel(2),elfiel(3)
cc      if( ifgrad .eq.1) write (icond,2821) elfiel(4),elfiel(5),elfiel(6)
cc     *                                    ,elfiel(7),elfiel(8),elfiel(9)
cc      write (icond,223)
c------------------------------------------------
  223 format(58(1H*))
  280 format
     *('****          CALCULATIONS OF THE HESSIAN             ****')
 2820 format
     *('****           external electric field is             ****'/
     *'****        fx=',f6.4,2x,' fy=',f6.4,2x,' fz=',f6.4,
     *'         ****')
 2821 format
     *('****      external electric field gradient is         ****'/
     *'****        xx=',f6.4,2x,' yx=',f6.4,2x,' yy=',f6.4,
     *'         ****'/
     *'****        zx=',f6.4,2x,' zy=',f6.4,2x,' zz=',f6.4,
     *'         ****')
c
      end
c======================================================================
c
      subroutine hess_memcheck(rhf)

      use memory

      implicit real*8 (a-h,o-z)
      Logical rhf
c---------------------------------------------------------------------
c
c  Total memory needed for gradient run is :
c
c  (0) forces on atoms 3*natoms for one-el., two-el. parts & total
c      3*(3*natoms)
c      ntri - for alpha/closed shell density
c             ditto for beta if open shell
c
c
c  (1) one-electron part ?????
c
c  (2) two-electron part : ncs*ncs  - density for prescreening(int_grad)
c                          limxmem  - two-el.int.deriv
c                          4*nbls*lnijkl*ngcd if ngcd>1 (transpose)
c
c---------------------------------------------------------------------
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /intlim/ limxmem,limblks,limpair
c-----------------------------------------------------------
      call getmem(0,lastx)
      call retmem(1)
c
      ntri=ncf*(ncf+1)/2
c
      mem_c= 3*(3*na) + ntri               ! common f1,f2,f3 + density
      If(.not.rhf) mem_c = mem_c + ntri    ! beta density
      mem_1= 3*(3*ntri)                    ! S1,Kinet1,V1  ????
      mem_2=           ncs*ncs
      mem_2= mem_2 + 10*mem_2/100          ! ??????
c
      needed=lastx + mem_c + mem_1 + mem_2
c
      if(needed.ge.lcore) then
         ioffset=igetival('ioffset')
         write(iout,210) needed-ioffset,lcore-ioffset
         call nerror(3,'Analytical Hessian Module',
     *              'Memory: needed and available ',
     *               needed-ioffset,lcore-ioffset)
      endif
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) then
        write(iout,209) needed,lcore
        return
      endif
c
  209 format (/1x,' Memory in BL for HESSIAN   calculations '/
     *            '  needed ',i10,'  available ',i10/)
  210 format (/1x,'common bl too small for hess  run, required =',i10,
     *3x,' available =',i10/)
c----------------------------
      end
c======================================================================
c
      subroutine read_dfeig(rhf,ncf,NAlpha,NBeta)

      use memory

      implicit real*8 (a-h,o-z)
c
c  reads in density & fock & eigenvec & eigenval.
c
c     common /big/bl(10000)
      Logical rhf
c
      ntri=(ncf*(ncf+1))/2
c----------------------------------------------------
      np1 = 1
      np4 = 4
c
      call getival('lsmat',lsmat)
      call getival('ldensi',lden)
      call getival('lfock0',lfoc)
      call getival('ldewsi',ldew)
      call getival('lvec',lvec)
      call getival('lval',lval)
      if(.not.rhf)then
        call getival('ldensB',mden)
        call getival('lfock0B',mfoc)
        call getival('lvecB',mvec)
        call getival('lvalB',mval)
      endif
c
c -- overlap
      call sudat(np1,'s matrix',ni)
      if(ni.gt.0) then
        call rea(bl(lsmat),ntri,np1,'s matrix')
      else
        call restart(np1,0)
        call sudat(np1,'s matrix',ni)
        if(ni.gt.0) then
           call rea(bl(lsmat),ntri,np1,'s matrix')
        else
          call nerror(1,'In <read_dfeig>',
     1   'Program cannot find the 0th-order overlap on <jobname>.14',
     2    ni,ni)
        endif
      endif
c
c
      if(rhf)then
c..............................................
        call sudat(np4,'den0_rhf',ni)
        if(ni.eq.0) call restart(np4,0)
c..............................................
        call rea(bl(lden),ntri,np4,'den0_rhf')
        call rea(bl(lfoc),ntri,np4,'fock_rhf')
      else
        call rea(bl(lden),ntri,np4,'dena_uhf')
        call rea(bl(lfoc),ntri,np4,'foca_uhf')
        call rea(bl(mden),ntri,np4,'denb_uhf')
        call rea(bl(mfoc),ntri,np4,'focb_uhf')
      endif
c----------------------------------------------------
c read-in the eigenvectors & eigenvalues :
c
      np4=4
      If(rhf) Then
c --     read in closed-shell MOs and orbital energies
         call rea (bl(lvec),ncf*ncf,np4,'evec_rhf')
         call rea (bl(lval),ncf    ,np4,'eval_rhf')
         call make_wdens(ncf,NAlpha,bl(lvec),bl(lval),bl(ldew),rhf)
      Else
c --     read in alpha MOs and orbital energies
         call rea (bl(lvec),ncf*ncf,np4,'evea_uhf')
         call rea (bl(lval),ncf    ,np4,'evaa_uhf')
         call make_wdens(ncf,NAlpha,bl(lvec),bl(lval),bl(ldew),rhf)
c --     now read in beta MOs and orbital energies
c   ** NEED TO CHECK THIS **    ! JB
         call rea (bl(mvec),ncf*ncf,np4,'eveb_uhf')
         call rea (bl(mval),ncf    ,np4,'evab_uhf')
         call make_wdens(ncf,-NBeta,bl(mvec),bl(mval),bl(ldew),rhf)
      EndIf
c----------------------------------------------------
      end
c======================================================================
      subroutine hess_total(rhf,natoms,ncf,bl)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      dimension bl(*)
c---------------------------------------------------------------------
c     natoms  - total number of atoms
c     ncf     - number of basis functions
c     bl      - storage foe everything
c---------------------------------------------------------------------
c  calculates two contributions (last two) to the final hessian
c  tr(d1*f1) - tr(w1*s1)
c
c  makes the full hessian matrix out of its upper triangle
c  writes it on a disk
c---------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
c---------------------------------------------------------------------
      call getival('lhess',lhess)
c
c Now fock1, over1, den1 & wen1 are dimensioned as (ntri,3,natoms)
c They are on a disk
c allocate memory for Fock1, S1, D1 & W1 (weighted density 1)
c for TWO atoms only
c
      call getmem(2*ntri3,lfock1)
      call getmem(2*ntri3,lover1)
      call getmem(2*ntri3,lden1)
      call getmem(2*ntri3,lwen1)
c---------------------------------------------------------------------
c Read in dens1 and wens1 matrices from a disk
c ( Ist-order density & weighted density)
c These were calculated and written on a disk for symmetry unique
c atoms only. Make D1 and W1 for all atoms and write them back on disk
c the lfock1 storage is used to store the beta density in the uhf case
c
      call getival('nsym',nsym)
      if(nsym.gt.0) then
         call getival('listall',listall)
         call getival('listsym',listsym)
         call getival('SymFunPr',ifp)
         call symm_dewe1(rhf,bl(ifp),bl(listall),bl(listsym),natoms,
     *                   bl(lden1),bl(lfock1),bl(lwen1),ncf,ntri)
      endif
c
c  calculates two contributions (last two) to the final hessian
c  tr(d1*f1) - tr(w1*s1)
c
      call calc_d1f1_w1s1(rhf,natoms,ntri,ncf,
     *                    bl(lden1),bl(lfock1),bl(lwen1),bl(lover1),
     *                    bl(lhess) )
c
      call retmem(4) ! lfock1,lover1,lden1, lwen1
c---------------------------------------------------------------------
c make full hessian out of its upper triangle part :
c
      call hess_full_up(bl(lhess),natoms)
c---------------------------------------------------------------------
      return
      end
c======================================================================
      subroutine calc_d1f1_w1s1(rhf,natoms,ntri, ncf,
     *                          den1,fock1, wen1,over1, hess)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c  calculates two contributions (last two) to the final hessian
c  tr(d1*f1) - tr(w1*s1)
c---------------------------------------------------------------------
c Input :
c
c rhf     - flag for rhf/uhf
c natoms  - number of atoms
c ntri    - ncf*(ncf+1)/2
c ncf     - basis set dimension
c  den1() - Ist-order density matrix
c fock1() - Ist-order fock matrix
c  wen1() _ Ist-order weighted density
c over1() - Ist-order overlap matrix
c
c Input/Output  :
c
c hess()  - final hessian
c---------------------------------------------------------------------
      logical rhf
      dimension den1(ntri,3,2),fock1(ntri,3,2)
      dimension wen1(ntri,3,2),over1(ntri,3,2)
      dimension hess(3,natoms,3,natoms)
c-----------------------------------------------------------------
      call getival('printh',nprint)
c-----------------------------------------------------------------
      ntri3=ntri*3
c
      nfile61=61
      nfile71=71
      nfile62=62
      nfile63=63
      nfile73=73
      nfile64=64
c-----------------------------------------------------------------
c Atom=Btom
c
      do iat=1,natoms
         call read1mat(nfile61,iat,ntri3,fock1(1,1,1))
         call read1mat(nfile62,iat,ntri3,over1(1,1,1))
         call read1mat(nfile63,iat,ntri3, den1(1,1,1))
         call read1mat(nfile64,iat,ntri3, wen1(1,1,1))
         do ixyz=1,3
            call spur(den1(1,ixyz, 1 ),fock1(1,ixyz,1),ncf,df1)
            call spur(wen1(1,ixyz, 1 ),over1(1,ixyz,1),ncf,ws1)
            dewe=df1-ws1
            hess(ixyz,iat,ixyz,iat)=hess(ixyz,iat,ixyz,iat)+dewe
            do jxyz=ixyz+1,3
                call spur(den1(1,ixyz, 1 ),fock1(1,jxyz, 1 ),ncf,df1)
                call spur(den1(1,jxyz, 1 ),fock1(1,ixyz, 1 ),ncf,fd1)
                call spur(wen1(1,ixyz, 1 ),over1(1,jxyz, 1 ),ncf,ws1)
                call spur(wen1(1,jxyz, 1 ),over1(1,ixyz, 1 ),ncf,sw1)
                df=df1+fd1
                ws=ws1+sw1
                dewe=df-ws
                dewe= dewe*0.5d0
                hess(ixyz,iat,jxyz,iat)=hess(ixyz,iat,jxyz,iat)+dewe
            enddo
         enddo
      enddo
c
c  add beta component in uhf case (weighted density is already OK)
c
      if(.not.rhf)then
        do iat=1,natoms
           call read1mat(nfile71,iat,ntri3,fock1(1,1,1))
           call read1mat(nfile73,iat,ntri3, den1(1,1,1))
           do ixyz=1,3
              call spur(den1(1,ixyz, 1 ),fock1(1,ixyz,1),ncf,df1)
              dewe=df1
              hess(ixyz,iat,ixyz,iat)=hess(ixyz,iat,ixyz,iat)+dewe
              do jxyz=ixyz+1,3
                  call spur(den1(1,ixyz, 1 ),fock1(1,jxyz, 1 ),ncf,df1)
                  call spur(den1(1,jxyz, 1 ),fock1(1,ixyz, 1 ),ncf,fd1)
                  df=df1+fd1
                  dewe=df
                  dewe= dewe*0.5d0
                  hess(ixyz,iat,jxyz,iat)=hess(ixyz,iat,jxyz,iat)+dewe
              enddo
           enddo
        enddo
      endif
c
c different atoms :
c
      do iat=1,natoms
         call read1mat(nfile61,iat,ntri3,fock1(1,1,1))
         call read1mat(nfile62,iat,ntri3,over1(1,1,1))
         call read1mat(nfile63,iat,ntri3, den1(1,1,1))
         call read1mat(nfile64,iat,ntri3, wen1(1,1,1))
         do jat=iat+1,natoms
            call read1mat(nfile61,jat,ntri3,fock1(1,1,2))
            call read1mat(nfile62,jat,ntri3,over1(1,1,2))
            call read1mat(nfile63,jat,ntri3, den1(1,1,2))
            call read1mat(nfile64,jat,ntri3, wen1(1,1,2))
            do ixyz=1,3
               do jxyz=1,3
                  call spur(den1(1,ixyz, 1 ),fock1(1,jxyz, 2 ),ncf,df1)
                  call spur(den1(1,jxyz, 2 ),fock1(1,ixyz, 1 ),ncf,fd1)
                  call spur(wen1(1,ixyz, 1 ),over1(1,jxyz, 2 ),ncf,ws1)
                  call spur(wen1(1,jxyz, 2 ),over1(1,ixyz, 1 ),ncf,sw1)
                  df=df1+fd1
                  ws=ws1+sw1
                  dewe=df-ws
                  dewe= dewe*0.5d0
                  hess(ixyz,iat,jxyz,jat)=hess(ixyz,iat,jxyz,jat)+dewe
               enddo
            enddo
         enddo
      enddo
c
c  add beta component in uhf case (weighted density is already OK)
c
      if(.not.rhf)then
        do iat=1,natoms
          call read1mat(nfile71,iat,ntri3,fock1(1,1,1))
          call read1mat(nfile73,iat,ntri3, den1(1,1,1))
          do jat=iat+1,natoms
            call read1mat(nfile71,jat,ntri3,fock1(1,1,2))
            call read1mat(nfile73,jat,ntri3, den1(1,1,2))
            do ixyz=1,3
              do jxyz=1,3
                call spur(den1(1,ixyz, 1 ),fock1(1,jxyz, 2 ),ncf,df1)
                call spur(den1(1,jxyz, 2 ),fock1(1,ixyz, 1 ),ncf,fd1)
                df=df1+fd1
                dewe=df
                dewe= dewe*0.5d0
                hess(ixyz,iat,jxyz,jat)=hess(ixyz,iat,jxyz,jat)+dewe
              enddo
            enddo
          enddo
        enddo
      endif
c----------------------------------------------------------------------
      end
c======================================================================
c printing routines :
c======================================================================
      subroutine fder_print(fder,na,ntri,factor)
      implicit real*8 (a-h,o-z)
      dimension fder(3,na,ntri)
c
      write(6,*)'           Derivative Fock matrices '
      write(6,*)'           X          Y           Z '
      do iat=1,na
         write(6,*) 'Atom no : ',iat
         do ij=1,ntri
            write(6,111) ij,fder(1,iat,ij)*factor,
     *                      fder(2,iat,ij)*factor,
     *                      fder(3,iat,ij)*factor
         enddo
      enddo
c
  111 format('ij=',i3,2x,3(f12.6,1x))
c
      end
c=====================================================================
c
      subroutine hess_print(hessian,na,factor,text)
      implicit real*8 (a-h,o-z)
      character*(*)text
      dimension hessian(3,na,3,na)
c
      call getival('iout',iout)
c
      write(iout,100)
  100 format(/61('-'))
c
      write(iout,*)'            The Analytical Hessian Matrix '

      write(iout,*) '                 ',  text
      write(iout,*)'  '
c
c     do iat=1,na
c        do jat=iat,na
      do iat=1,na
         do jat=1  ,na
             write(iout,*) 'Atoms =', iat,jat
c            if(iat.eq.jat) then
c               xx=hessian(1,iat,1,jat)*factor
c               xy=hessian(1,iat,2,jat)*factor
c               xz=hessian(1,iat,3,jat)*factor
c               yy=hessian(2,iat,2,jat)*factor
c               yz=hessian(2,iat,3,jat)*factor
c               zz=hessian(3,iat,3,jat)*factor
c               write(iout,103) 'xx xy xz : ',xx,xy,xz
c               write(iout,102) '   yy yz : ',   yy,yz
c               write(iout,101) '      zz : ',      zz
c            else
                xx=hessian(1,iat,1,jat)*factor
                xy=hessian(1,iat,2,jat)*factor
                xz=hessian(1,iat,3,jat)*factor
                yx=hessian(2,iat,1,jat)*factor
                yy=hessian(2,iat,2,jat)*factor
                yz=hessian(2,iat,3,jat)*factor
                zx=hessian(3,iat,1,jat)*factor
                zy=hessian(3,iat,2,jat)*factor
                zz=hessian(3,iat,3,jat)*factor
                write(iout,103) 'xx xy xz : ',xx,xy,xz
                write(iout,103) 'yx yy yz : ',yx,yy,yz
                write(iout,103) 'zx zy zz : ',zx,zy,zz
c            endif
         enddo
      enddo
c
      write(iout,100)
c
  103 format(a10,3(f15.3,1x) )
c 103 format(a10,3(f15.6,1x) )
c 102 format(a10,13x,2(f15.6,1x) )
c 101 format(a10,26x,1(f15.6,1x) )
c
      end
c==============================================================
c
      subroutine trsp_inplace(bl,array,n1,n2)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension array(*)
c
c a general matrix array(n1,n2) transposed into array(n2,n1)
c
      call getmem(n1*n2,itemp)    !  scratch array
c
      call trspmo(array,n1, bl(itemp),n2)
      call dcopy(n1*n2,bl(itemp),1,array,1)
c
      call retmem(1)    ! release scratch memory
c
      end
c==============================================================
c
      subroutine drumh (d,ncf,ifi,txt)
      implicit real*8 (a-h,o-z)
      character*8 txt
      dimension d(1)
c
      write (ifi,20) txt
c----------------------------------------------------------------------c
c do not print anything if it is too long :
      if(ncf.gt.200) then
          write(ifi,*)
     *  ' it is not printed out for NCF=',ncf,' > 200 '
         return
      endif
c----------------------------------------------------------------------c
      n=ncf
      ja=1
      do 10 i=1,n
         je=ja+i-1
         write (ifi,30) i,(d(j),j=ja,je)
         ja=je+1
   10 continue
      return
c
   20 format (30x,3h***,a8,3h***)
c  30 format (1x,i4,2x,7e12.3,/,(7x,7e12.3))
   30 format (1x,i4,2x,7f12.7,/,(7x,7f12.7))
c
      end
c=======================================================================
      subroutine memory_status(text)

      use memory

      character*(*) text
c-------------------------------------------
      call getival('iout',iout)
c-------------------------------------------
c memory status :
c
      call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
      write(iout,5001)
 5001 format(72('-'))
c
      write(iout,1100) text,nreq,nmark,lastadr-ioffset,mxmem,memtot
 1100 format(' Memory status:',a20,/,
     *' request number        =',i15,/,
     *' memory marks          =',i15,/,
     *' last used address     =',i15,/,
     *' high water            =',i15,/,
     *' total available memory=',i15)
c
      write(iout,5001)
c-------------------------------------------
      call f_lush(iout)
c-------------------------------------------
      end
c=======================================================================
c
      subroutine find_unqat(bl,natoms)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
c
c allocate memory for symmetry unique atoms
c                     symmetry equivalent atoms
c                     symmetry operation relating two atoms
c
      call getint(natoms,listunq)
      call getint(natoms,listall)
      call getint(natoms,listsym)
      call getint(natoms,listrea)
c
      call getival('nsym',nsym)
c
      nupair=1
      if(nsym.gt.0) call getival('SymNuPr1',nupair)
      call make_listunq(natoms,nsym,bl(nupair),natunq,
     *                  bl(listunq),bl(listall),bl(listsym),bl(listrea))
c
c save pointers :
c
      call setival('natunq',natunq)   ! number of unique atoms
      call setival('listunq',listunq)
      call setival('listall',listall)
      call setival('listsym',listsym)
      call setival('listrea',listrea)
c
      end
c============================================================
c
      subroutine make_listunq(natoms,nsym,nupair,
     *                        natunq,lunq, lall,lsym,lrea)
      dimension nupair(natoms,nsym)
      dimension lunq(natoms)   ! symmetry unique atoms
      dimension lall(natoms)   ! symmetry equivalent atoms
      dimension lsym(natoms)   ! related by THIS symmetry operation
      dimension lrea(natoms)   ! mapping from real to unique
c
      do iat=1,natoms
         lunq(iat)=iat
         lall(iat)=iat
         lsym(iat)=0
         lrea(iat)=0
      enddo
      natunq=natoms
      if(nsym.eq.0) then
c        write(6,*)' nsym=',nsym,' natoms=',natoms,' unique at=',natunq
c        do ii=1,natunq
c          write(6,*)' unq at=',ii,'  real atom no=',lunq(ii)
c        enddo
      do iat=1,natoms
         lrea(iat)=iat
      enddo
         return
      endif
c
c finds symmetry unique atoms
c
      iu=1
      do iat=2,natoms
         iax=10000
         do ns=1,nsym
            iat1=nupair(iat,ns)
            if(iat1.lt.iax) iax=iat1
         enddo
         iax=min(iat,iax)
         if(iax.gt.lunq(iu) ) then
            iu=iu+1
            lunq(iu)=iax
         endif
      enddo
      natunq=iu
c
c     write(6,*)' nsym=',nsym,' natoms=',natoms,' unique at=',natunq
c     do ii=1,natunq
c       write(6,*)' unq at=',ii,'  real atom no=',lunq(ii)
c     enddo
c
      do ii=1,natoms
         lall(ii)=0
         lsym(ii)=0
         lrea(ii)=0
      enddo
c
      do iu=1,natunq
         ir=lunq(iu)
         lall(ir)=ir
         lrea(ir)=iu
      enddo
c
c     do ii=1,natoms
c       write(6,*) 'ii=',ii,' listall(ii)=',lall(ii),' symop=',lsym(ii)
c     enddo
c
      do ia=1,natoms
         iat=lall(ia)
         if(iat.eq.0) then
            do ns=1,nsym
               ia1=nupair(ia,ns)
               iat1=lall(ia1)
               if(iat1.gt.0) then
                  lall(ia)=-ia1
                  lsym(ia)=ns
                  go to 123
               endif
            enddo
  123       continue
         endif
      enddo
c
c     do ii=1,natoms
c       write(6,*) 'ii=',ii,' listall(ii)=',lall(ii),' symop=',lsym(ii),
c    *              ' listrea(ii)=',lrea(ii)
c     enddo
c
c     do iat=1,natoms
c        write(6,*)' atom =',iat
c        do ns=1,nsym
c           ia1=nupair(iat,ns)
c           write(6,*)'    symmetry=',ns,' image=',ia1
c        enddo
c     enddo
c
      end
c============================================================
c
      subroutine symm_dewe1(rhf,ifp,lall,lsym,natoms,den1,den1B,
     $                      wen1,ncf,ntri)
      implicit real*8 (a-h,o-z)
      logical rhf
      common /negxyz/ negx(7),negy(7),negz(7)
      dimension ifp(7,ncf)
      dimension lall(natoms)   ! symmetry equivalent atoms
      dimension lsym(natoms)   ! related by THIS symmetry operation
      dimension den1(ntri,3,2),den1B(ntri,3,2),wen1(ntri,3,2)
cccc  dimension den1(ntri,3,natoms),wen1(ntri,3,natoms)
c-------------------------------------------------------------------
c Makes Ist-order density & weighted density matrices for ALL atoms.
c I case of symmetry D1 & W1 have been calculated for UNIQUE atoms only.
c
c INPUT :
c  rhf  flag for rhf/uhf
c  ifp(isym,ifun) - gives symmetry ISYM image of the ifun
c  lall(iat)      - symmetry equivalent atoms
c  lsym(iat)      - by this symmetry operation
c  natoms         - total number of atoms
c
c INPUT/OUTPUT
c  den1           - Ist-order density (closed shell or alpha)
c  den1B          - Ist-order density (uhf beta component)
c  wen1           - Ist-order weighted density
c-------------------------------------------------------------------
      ntri3=ntri*3
c
      do iat=1,natoms
         ia1=lall(iat)
         if(ia1.lt.0) then
            ias=-ia1              ! symmetry unique, equivalent to iat
            isymop=lsym(iat)      ! symmetry oper. making ias->iat
c           read in from a disk d1 & w1 for symm. unique atom ias
            call read1mat(63,ias,ntri3, den1(1,1,1))
            if(.not.rhf)call read1mat(73,ias,ntri3, den1B(1,1,1))
            call read1mat(64,ias,ntri3, wen1(1,1,1))
c
            ij=0
            do icf=1,ncf
               ic1=ifp(isymop,icf)
               fci=1.d0
               if(ic1.lt.0) then
                  ic1=-ic1
                  fci=-1.d0
               endif
               do jcf=1,icf
                  jc1=ifp(isymop,jcf)
                  fcj=1.d0
                  if(jc1.lt.0) then
                     jc1=-jc1
                     fcj=-1.d0
                  endif
c
                  fct=fci*fcj
c
                  ij=ij+1
                  ij1=ic1*(ic1-1)/2+jc1
                  if(jc1.gt.ic1) ij1=jc1*(jc1-1)/2+ic1
c................................................................
c                 den1(ij1,1,iat)=fct*den1(ij,1,ias)*negx(isymop)
c                 den1(ij1,2,iat)=fct*den1(ij,2,ias)*negy(isymop)
c                 den1(ij1,3,iat)=fct*den1(ij,3,ias)*negz(isymop)
c                 wen1(ij1,1,iat)=fct*wen1(ij,1,ias)*negx(isymop)
c                 wen1(ij1,2,iat)=fct*wen1(ij,2,ias)*negy(isymop)
c                 wen1(ij1,3,iat)=fct*wen1(ij,3,ias)*negz(isymop)
c................................................................
c
                  den1(ij1,1, 2 )=fct*den1(ij,1, 1 )*negx(isymop)
                  den1(ij1,2, 2 )=fct*den1(ij,2, 1 )*negy(isymop)
                  den1(ij1,3, 2 )=fct*den1(ij,3, 1 )*negz(isymop)
                  if(.not.rhf)then
                    den1B(ij1,1, 2 )=fct*den1B(ij,1, 1 )*negx(isymop)
                    den1B(ij1,2, 2 )=fct*den1B(ij,2, 1 )*negy(isymop)
                    den1B(ij1,3, 2 )=fct*den1B(ij,3, 1 )*negz(isymop)
                  endif
                  wen1(ij1,1, 2 )=fct*wen1(ij,1, 1 )*negx(isymop)
                  wen1(ij1,2, 2 )=fct*wen1(ij,2, 1 )*negy(isymop)
                  wen1(ij1,3, 2 )=fct*wen1(ij,3, 1 )*negz(isymop)
               enddo    !    jcf=1,icf
            enddo       !    icf=1,ncf
c           write on a disk d1 & w1 for equivalent atom iat
            call save1mat(63,iat,ntri3, den1(1,1,2))
            if(.not.rhf)call save1mat(73,iat,ntri3, den1B(1,1,2))
            call save1mat(64,iat,ntri3, wen1(1,1,2))
         endif
      enddo             ! over atoms
c
      end
c========================================================
      subroutine open4hess(ncs,ncf,ntri,natom,rhf)
      implicit real*8(a-h,o-z)
      logical rhf
      character*256 jobname,scrf,filename
      common /job/jobname,lenJ
c----------------------------------------------------
c open direct access fileis for hessian calculations
c
c Input
c ndim - record length
c rhf  - rhf flag
c----------------------------------------------------
      ncf8=ncf*ncf*8
      ncs8=ncs*ncs*8
c
      ndim=ntri*3
      lrec = ndim*8          ! record length in bytes
      nfile = 60
      nfileu= 70
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
c
      len = len1 + 6
c
c new file for DS,SD and Brillouin matrices, needed for convergence criterion
c
      filename = scrf(1:len1)//'.brill'
      open (unit=nfile+0,file=filename(1:len),
     *      form='unformatted',access='direct',recl=ncf8)
c
      filename = scrf(1:len1)//'.fock1'
      open (unit=nfile+1,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
      if(.not.rhf)then
        filename = scrf(1:len1)//'.fock1B'
        open (unit=nfileu+1,file=filename(1:len+1),
     *        form='unformatted',access='direct',recl=lrec)
      endif
c
      filename = scrf(1:len1)//'.over1'
      open (unit=nfile+2,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.dens1'
      open (unit=nfile+3,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
      if(.not.rhf)then
        filename = scrf(1:len1)//'.dens1B'
        open (unit=nfileu+3,file=filename(1:len+1),
     *        form='unformatted',access='direct',recl=lrec)
      endif
c
      filename = scrf(1:len1)//'.wens1'
      open (unit=nfile+4,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      filename = scrf(1:len1)//'.dcon1'
      open (unit=nfile+5,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
      if(.not.rhf)then
        filename = scrf(1:len1)//'.dcon1B'
        open (unit=nfileu+5,file=filename(1:len+1),
     *        form='unformatted',access='direct',recl=lrec)
      endif
c
      filename = scrf(1:len1)//'.resi1'
      open (unit=nfile+6,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
      if(.not.rhf)then
        filename = scrf(1:len1)//'.resi1B'
        open (unit=nfileu+6,file=filename(1:len+1),
     *        form='unformatted',access='direct',recl=lrec)
      endif
c
      filename = scrf(1:len1)//'.lasd1'
      open (unit=nfile+7,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
      if(.not.rhf)then
        filename = scrf(1:len1)//'.lasd1B'
        open (unit=nfileu+7,file=filename(1:len+1),
     *        form='unformatted',access='direct',recl=lrec)
      endif
c
      ndim1=natom*natom*9
      lrec1= ndim1*8          ! record length in bytes
c
      filename = scrf(1:len1)//'.hess1'
      call getival('resta',irestart)
      if(irestart.eq.0) then
         open (unit=nfile+8,file=filename(1:len),
     *         form='unformatted',access='direct',recl=lrec1)
      else
         open (unit=nfile+8,file=filename(1:len),status='old',
     *         form='unformatted',access='direct',recl=lrec1,err=95)
      endif
c
c new file for screening density
c
      filename = scrf(1:len1)//'.screen'
      open (unit=nfile+9,file=filename(1:len+1),
     *      form='unformatted',access='direct',recl=ncs8)
c
c file for perturbed coefficient matrix in VCD
c
      call tstival('vcd',ivcd)
      If(ivcd.ne.0) then
        ncoefile=79
        call getival('nocc',nocc)
        ncoe3=ncf*nocc*3*8   ! record length for the 3 mx
        filename = scrf(1:len1)//'.pcoef'
        open (unit=ncoefile,file=filename(1:len),
     *        form='unformatted',access='direct',recl=ncoe3)
      EndIf
c
      return
c
   95 continue
      Call nerror(2,'Analytical Hessian module',
     $     'Restart Specified and old <hess1> File Does Not Exist',
     $      0,0)
      end
c======================================================================
      subroutine close4hess(rhf)
      logical rhf
c-----------------------------------------------------------------
c close files open for hess
      call getival('savef',isavefil)
c
      if(isavefil.eq.0) then
         close (60,status='delete')
         close (61,status='delete')
         close (62,status='delete')
         close (63,status='delete')
         close (64,status='delete')
         close (65,status='delete')
         close (66,status='delete')
         close (67,status='delete')
         close (68,status='delete')
         close (69,status='delete')
         if(.not.rhf)then
           close (71,status='delete')
           close (73,status='delete')
           close (75,status='delete')
           close (76,status='delete')
           close (77,status='delete')
         endif
      else
         close (60,status='keep')
         close (61,status='keep')
         close (62,status='keep')
         close (63,status='keep')
         close (64,status='keep')
         close (65,status='keep')
         close (66,status='keep')
         close (67,status='keep')
         close (68,status='keep')
         close (69,status='keep')
         if(.not.rhf)then
           close (71,status='keep')
           close (73,status='keep')
           close (75,status='keep')
           close (76,status='keep')
           close (77,status='keep')
         endif
      endif
c
      call tstival('vcd',ivcd)
      if(ivcd.ne.0) close (79,status='keep')
c---------------------------------------------------------------------
      end
c======================================================================
c
      SUBROUTINE trsp_disk(IUnit,NAtoms,ntri,bl,A)
      IMPLICIT REAL*8(a-h,o-z)
C
C  Transpose either the Fock or overlap derivative matrix (A)
C  and write to direct access file IUnit  PER ATOM
C    A(3,NAtoms,ntri)  --->  A(ntri,3,natoms)
C
C  ARGUMENTS
C
C  IUnit   -  Unit number to write to (must be open)
C  NAtoms  -  total number of atoms
C  ntri    -  number of elements in lower triangle
C  bl      -  scratch storage
C  A       -  matrix to transform
C
C
      DIMENSION A(3,NAtoms,ntri),bl(ntri,3)
C
      DO 20 IAtm=1,NAtoms
      DO 10 I=1,ntri
      bl(I,1) = A(1,IAtm,I)
      bl(I,2) = A(2,IAtm,I)
      bl(I,3) = A(3,IAtm,I)
 10   CONTINUE
C
C  write transposed matrix for this atom to disk
C
      WRITE(Unit=IUnit,rec=IAtm) bl
C
 20   CONTINUE
C
      RETURN
      END
c======================================================================
c
      subroutine read1mat(nfile,iat,ntri3,xmat)
      implicit real*8(a-h,o-z)
      dimension xmat(ntri3)
         read(unit=nfile,rec=iat) xmat
      end
c======================================================================
c
      subroutine save1mat(nfile,iat,ntri3,xmat)
      implicit real*8(a-h,o-z)
      dimension xmat(ntri3)
         write(unit=nfile,rec=iat) xmat
      end
c======================================================================
c
      subroutine hess_full_lo(hess,na)
      implicit real*8 (a-h,o-z)
      dimension hess(3,na,3,na)
c
      do i=1,na
         do ixyz=1,3
            do j=1,i
               jxyz_e=3
               if(j.eq.i) jxyz_e=ixyz
               do jxyz=1,jxyz_e
                  hess(jxyz,j,ixyz,i)=hess(ixyz,i,jxyz,j)
               enddo
            enddo
         enddo
      enddo
c
      end
c======================================================================
c
      subroutine hess_full_up(hess,na)
      implicit real*8 (a-h,o-z)
      dimension hess(3,na,3,na)
c
      do i=1,na
         do ixyz=1,3
            do j=i,na
               if(j.eq.i) then
                  do jxyz=ixyz,3
                     hess(jxyz,j,ixyz,i)=hess(ixyz,i,jxyz,j)
                  enddo
               else
                  do jxyz=1,3
                     hess(jxyz,j,ixyz,i)=hess(ixyz,i,jxyz,j)
                  enddo
               endif
            enddo
         enddo
      enddo
c
      end
c======================================================================
c
      SUBROUTINE Hess_Tot(rhf,NAtoms, NQ,     IUNQ,   ncf,    ntri,
     $                    den1,   wen1,   fock1,  over1,  HESS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Forms the final contribution to the Hessian matrix
C  tr(den1*fock1) - tr(wen1*over1)
C  following solution of CPHF equations
C  ** FORMS HESSIAN COLUMNS FOR SYMMETRY-UNIQUE ATOMS ONLY **
C
C  ARGUMENTS
C
c  rhf     -  flag for rhf/uhf
C  NAtoms  -  total number of atoms
C  NQ      -  number of symmetry unique atoms
C  IUNQ    -  list of symmetry-unique atoms
C  ncf     -  number of basis functions
C  ntri    -  ncf*(ncf+1)/2
C  den1    -  storage for 1st-order density matrix per atom
C  wen1    -    ditto weighted 1st-order density matrix
C  fock1   -    ditto derivative Fock matrix
C  over1   -    ditto derivative overlap matrix
C  HESS    -  Hessian matrix
C
C
      logical rhf
      DIMENSION den1(ntri,3),wen1(ntri,3),fock1(ntri,3),over1(ntri,3)
      DIMENSION IUNQ(NQ),HESS(3,NAtoms,3,NAtoms)
C
      PARAMETER (Zero=0.0d0)
C
C
      ntri3 = 3*ntri
C
C  make full Hessian from upper triangle
C
      call hess_full_up(HESS,NAtoms)
C
C  Loop over symmetry-unique atoms
C
      DO 20 IQ=1,NQ
      IAtm = IUNQ(IQ)
C
C  read in den1/wen1 for this symmetry-unique atom
C
      call read1mat(63,iatm,ntri3,den1(1,1))
      call read1mat(64,iatm,ntri3,wen1(1,1))
C
      DO 10 JAtm=1,NAtoms
C
C  now read in fock1/over1 for each in turn
C
      call read1mat(61,jatm,ntri3,fock1(1,1))
      call read1mat(62,jatm,ntri3,over1(1,1))
C
C  form appropriate trace and add to Hessian
C
      do ixyz=1,3
      do jxyz=1,3
        call spur(den1(1,ixyz),fock1(1,jxyz),ncf,df1)
cc        call spur(den1(1,jxyz),fock1(1,ixyz),ncf,fd1)
        call spur(wen1(1,ixyz),over1(1,jxyz),ncf,ws1)
cc        call spur(wen1(1,jxyz),over1(1,ixyz),ncf,sw1)
cc        df=df1+fd1
cc        ws=ws1+sw1
cc        dewe=df-ws
cc        dewe= dewe*0.5d0
        dewe = df1-ws1
        HESS(ixyz,iatm,jxyz,jatm) = HESS(ixyz,iatm,jxyz,jatm) + dewe
      enddo
      enddo
c
 10   CONTINUE
 20   CONTINUE
c
c  for the uhf case, add the beta component of the density-fock
c  product
c
      if(.not.rhf)then
        do iq=1,nq
          iatm = iunq(iq)
          call read1mat(73,iatm,ntri3,den1(1,1))
          do jatm=1,natoms
            do ixyz=1,3
            do jxyz=1,3
              call read1mat(71,jatm,ntri3,fock1(1,1))
              call spur(den1(1,ixyz),fock1(1,jxyz),ncf,df1)
              dewe = df1
              HESS(ixyz,iatm,jxyz,jatm)=HESS(ixyz,iatm,jxyz,jatm)+dewe
            enddo
            enddo
          enddo
        enddo
      endif
C
C  now copy rows into columns
C
      DO 40 IQ=1,NQ
      IAtm = IUNQ(IQ)
      DO 30 JAtm=1,NAtoms
      do ixyz=1,3
      do jxyz=1,3
      HESS(jxyz,JAtm,ixyz,IAtm) = HESS(ixyz,IAtm,jxyz,JAtm)
      enddo
      enddo
 30   CONTINUE
 40   CONTINUE
C
C  zero out all non-symmetry unique blocks
C
      IQ = 1
      DO 60 IAtm=1,NAtoms
      If(IAtm.EQ.IUNQ(IQ)) Then
        IQ = IQ+1
        GO TO 60
      EndIf
      JQ = 1
      DO 50 JAtm=1,NAtoms
      If(JAtm.EQ.IUNQ(JQ)) Then
        JQ = JQ+1
        GO TO 50
      EndIf
      do ixyz=1,3
      do jxyz=1,3
      HESS(jxyz,JAtm,ixyz,IAtm) = Zero
      enddo
      enddo
 50   CONTINUE
 60   CONTINUE
C
      RETURN
      END
c======================================================================
c
      SUBROUTINE GetISYM(NAtoms,NQ,IUNQ,ISYM)
      IMPLICIT INTEGER(A-Z)
      DIMENSION IUNQ(NQ),ISYM(NAtoms)
C
C  generate ISYM array
C
      CALL IZeroIT(ISYM,NAtoms)
      DO 25 I=1,NQ
      II = IUNQ(I)
      ISYM(II) = 1
 25   CONTINUE
C
      RETURN
      END
c======================================================================
      SubRoutine GetIntens(rhf,nAtoms,nCF,nTri,nCS,na,inx,
     &           DMDer,Buff,Dens1,Dens1B,DipIn,blib,blin)

      use memory

      Implicit None
*---- List of arguments
      logical rhf
      Integer nAtoms,nCF,nTri,na,ncs
      Integer inx(*)
      Real*8 DMDer(3,nAtoms,3),Buff(3,nAtoms,3)
      Real*8 Dens1(nTri,3),Dens1B(nTri,3),DipIn(nTri,3)
      Real*8 blib(13,*),blin(5,*)
*---- Local variables
      Real*8 AtChar,Trace,PlDer(1)
      Integer iAt,iQ,i,ii,iRComp,iDComp
      Logical Tst
      Character*60 Comment
*---- Local variables: symmetry information
      Integer MSymOp
      Parameter (MSymOp=120)
      Character*4 Group
      Integer nAtms,nTrans,nDeg,nQ
      Integer lRM,iRM,lUnq,iUnq,lTrans,iTrans,lEqAtm,iEqAtm
      Integer lScr1,iScr1,lScr2,iScr2,lSym,iSym
      Real*8 Thrsh
*---- Common
      Character*256 JobName,ScrDir
      Integer lenJ,len1
      Common /job/JobName,lenJ
*
c     Real*8 bl
c     Common /big/bl(1000)
*
      Tst=.false.
*
      If (Tst) Then
         Write(6,*)
         Write(6,*) '***** SubRoutine  GetIntens *****'
         Write(6,*) '*** Dipole moment derivatives ***'
         Write(6,*)
         call f_lush(6)
      End If
*
      Call mmark
c
*---- Allocate memory for symmetry information
      lRM=3*3
      Call getmem(lRM,iRM)
      lUnq=nAtoms
      Call getmem(lUnq,iUnq)     ! should be changed to INTBL
      lTrans=3*3*MSymOp
      Call getmem(lTrans,iTrans)
      lEqAtm=nAtoms*MSymOp
      Call getmem(lEqAtm,iEqAtm) ! should be changed to INTBL
      lSym=nAtoms
      Call getmem(lSym,iSym)     ! should be changed to INTBL
*---- Allocate scratch space
      lScr1=9*nAtoms*nAtoms
      Call getmem(lScr1,iScr1)
      lScr2=9*nAtoms*nAtoms
      Call getmem(lScr2,iScr2)
*---- Read the symmetry information
      nAtms=nAtoms
      Call RdSym(.true.,nAtms,bl(iRM),Group,nTrans,nDeg,nQ,
     &           bl(iUnq),bl(iTrans),bl(iEqAtm))
*
      If (Tst) Then
         Write(6,*) 'nAtms = ',nAtms
         Write(6,*) 'Group = ',Group
         Write(6,*) 'nTrans= ',nTrans
         Write(6,*) 'nDeg  = ',nDeg
         Write(6,*) 'nQ    = ',nQ
         call f_lush(6)
      End If
*
*---- Get the dipole moment integrals
      Call IntOn(3,na,DipIn(1,1),inx,1,0,blib,blin,ncs)
      Call IntOn(3,na,DipIn(1,2),inx,2,0,blib,blin,ncs)
      Call IntOn(3,na,DipIn(1,3),inx,3,0,blib,blin,ncs)
*
*---- Open the file with the 1-st order density matrix
      Call GetChVal('scrdir',ScrDir)
      Call RmBlan2(ScrDir,256,len1)
c     Open(unit=63,
c    &     file=ScrDir(1:len1)//JobName(1:lenJ)//'.dens1',
c    &     status='old',
c    &     form='unformatted',
c    &     access='direct',
c    &     recl=3*nTri*8)
c     if(.not.rhf) Open(unit=73,
c    &     file=ScrDir(1:len1)//JobName(1:lenJ)//'.dens1B',
c    &     status='old',
c    &     form='unformatted',
c    &     access='direct',
c    &     recl=3*nTri*8)
*
*---- Set the dipole moment derivatives matrix to zero
      Call DCopy(3*3*nAtoms,0.d0,0,DMDer,1)
*
*---- Get the dipole moment derivatives
*---- Storage: DMDer(1:Rx,Ry,Rz ; 1:nAtoms ; 1:Mx,My,Mz)
      Do iQ=1,nQ
         Call GetAtm(iQ,iAt,bl(iUnq))
         AtChar=blin(1,iAt)
         If (Tst) Then
            Write(6,*)
            Write(6,*) 'For atom ',iAt,' the charge is ',AtChar
         End If
*
*------- Read the 1-st order density matrix for atom iAt
         Call Read1Mat(63,iAt,nTri*3,Dens1)
         if(.not.rhf)then
           Call Read1Mat(73,iAt,nTri*3,Dens1B)
           Call AddVEC(ntri*3,Dens1,Dens1B,Dens1)
         endif
*
         Do iDComp=1,3
            Do iRComp=1,3
*------------- Accumulate the nuclear contribution
               If (iRComp.eq.iDComp) DMDer(iRComp,iAt,iDComp)=
     &         DMDer(iRComp,iAt,iDComp)+AtChar
*------------- Accumulate contribution from the 1st-order density matrix
               Call Spur(Dens1(1,iRComp),DipIn(1,iDComp),nCF,Trace)
               DMDer(iRComp,iAt,iDComp)=DMDer(iRComp,iAt,iDComp)-Trace
            End Do
         End Do
*
*------- Print out the dipole moment derivatives (test purpose)
         If (Tst) Then
            Write(6,*) 'The dipole moment derivatives are'
            Write(6,*)'           X                 Y                 Z'
            Write(6,*)
     &       '        -------------------------------------------------'
            Do iDComp=1,3
               If (iDComp.eq.1) Comment='M(x)'
               If (iDComp.eq.2) Comment='M(y)'
               If (iDComp.eq.3) Comment='M(z)'
               Write(6,'(a4,3(2x,f16.10))') Comment,
     &            (DMDer(iRComp,iAt,iDComp),iRComp=1,3)
            End Do
         End If
*
      End Do
*
      Call GetIS(bl(iSym),bl(iUnq),nAtoms,nQ)
*
      Call FulDip(nAtoms,nTrans,bl(iSym),bl(iEqAtm),bl(iTrans),
     &            bl(iScr1),bl(iScr2),DMDer)
*
      Thrsh=1.d-6
      Call SymDip(nAtoms,nTrans,bl(iEqAtm),bl(iTrans),
     &            bl(iScr1),bl(iScr2),Thrsh,DMDer)
*
*---- Read the D0*Mx centribution and add to the D1*M
      Open(Unit=40,File=jobname(1:lenJ)//'.deriv_temp',
     $     Form='formatted',Status='old')
      Do iAt=1,nAtoms
         Do i=1,3
            Read(40,*) Buff(i,iAt,1),Buff(i,iAt,2),Buff(i,iAt,3)
         End Do
      End Do
      Close(Unit=40,Status='delete')
ckw   Close(Unit=40,Status='keep')
      Call DAXPY(3*nAtoms*3,1.d0,Buff,1,DMDer,1)
*
*---- Write the dipole moment derivatives to disc
      Call WrDeriv(3*nAtoms,DMDer,PlDer,.true.,.false.)
*
*---- Close and delete the 1-st order density matrix file
      Close(63,status='delete')
      if(.not.rhf)Close(73,status='delete')
*
      Call retmark
*
      Return
      End
*-----------------------------------------------------------------------
      SubRoutine GetAtm(iQ,iAt,iUnq)
      Implicit None
*---- List of arguments
      Integer iQ,iAt
      Integer iUnq(*)
*
      iAt=iUnq(iQ)
*
      End
*-----------------------------------------------------------------------
      SubRoutine GetIS(iSym,iUnq,nAtoms,nQ)
      Implicit None
*---- List of arguments
      Integer nAtoms,nQ
      Integer iSym(nAtoms),iUnq(nAtoms)
*---- Local variables
      Integer i,ii
*
      Call IZeroIT(iSym,nAtoms)
      Do i=1,nQ
         ii = iUnq(i)
         iSym(ii) = 1
      End Do
*
      End
*-----------------------------------------------------------------------
      subroutine check4ld3(threg2,thres1,thres2,thrs,doit)
      implicit real*8 (a-h,o-z)
      logical doit
c--------------------------------------------------------------------
c used for analytical hessian only (like check4ld)
c--------------------------------------------------------------------
c inp/out : integral thresholds threg2,thres1,thres2 & cphf threshold thrs
c--------------------------------------------------------------------
c check for linear dependencies in the basis set :
c
      call getival('nsym',nsym)
      call getival('iout',iout)
      call tstrval('xlows',lows)
      if(lows.eq.1) then
         call getrval('xlows',xlow)  ! lowest eigenvalue of S
      else
         xlow=0.1d0                  ! when running without scf
      endif
c--------------------------------------------------------------------
c new thresholds (after testing job from David Anick)
c H41O20/6-311++G-dp (xlow=8.5*10-6)
c
c--------------------------------------------------------------------
      if(0.5d-4.LT.xlow) RETURN  ! everything should be OK
c--------------------------------------------------------------------
c     if(xlow.ge.0.5d-4) then   ! default
c        threg2,thres1,thres2 & cphf threshold thrs
c        10^-9  10^-10 10^-8                   10^-5
c     endif
c--------------------------------------------------------------------
c if the default thresholds are changed via the input then 
c accept what user wants and ignore threshold sharpening
c
c
      if( 0.5d-4.ge.xlow) then
         write(iout,*)'  '
         write(iout,*)
     *   '              Near Linear Dependency in the Basis Set'
       
         if(.not.doit) then
            write(iout,*)
     *   '                integral threshold would be sharper '
            write(iout,*)
     *   '                were it not for input thresholds    '
            return
         endif
      endif
c--------------------------------------------------------------------
c
      if(0.5d-4.ge.xlow .and. xlow.gt.1.0d-5) then 
c def    threg2=1.0d-9
c def    thres1=1.0d-10
         thres2=1.0d-9
c def    thrs  =1.0d-5
         write(iout,153) threg2,thres1,thres2
      endif
c
      if(1.0d-5.ge.xlow .and. xlow.gt.1.0d-6) then 
c def    threg2=1.0d-9
         thres1=1.0d-11
         thres2=1.0d-10
c def    thrs  =1.0d-5
         write(iout,153) threg2,thres1,thres2
      endif
c
      if(1.0d-6.ge.xlow .and. xlow.gt.1.0d-7) then 
         threg2=1.0d-10
         thres1=1.0d-11
         thres2=1.0d-11
c def    thrs  =1.0d-5
         write(iout,153) threg2,thres1,thres2
      endif
c
      if(1.0d-7.ge.xlow                     ) then 
         threg2=1.0d-11
         thres1=1.0d-12
         thres2=1.0d-12
         thrs  =5.0d-5
         write(iout,153) threg2,thres1,thres2
      endif
c
      RETURN
c
c--------------------------------------------------------------------
c following is not used any more
c-------------------------------
c older setup; not good enough in the case of Anick's job
c
c after testing on annulene/6-31++G :
c Lowest eigenvalue of the overlap matrix  0.626716E-07
c
c
      call getival('nsym',nsym)
      rsymm=1.d0
      if(nsym.gt.0) rsymm=dble(nsym)
c
      if(xlow .lt. 5.0d-6 ) then
         write(iout,*)'  '
         write(iout,*)
     *   '              Near Linear Dependency in Basis Set'
       if(.not.doit) then
         write(iout,*)
     *   '              integral threshold would be sharpen '
         write(iout,*)
     *   '              if not overwritten in the input deck'
         return
       endif
c
         accur=xlow*1.0d-6
         accur=accur/rsymm
c
         thres1=min(thres1,accur)
         if(thres1.lt.1.0d-12) thres1=1.0d-12
         thres2=thres1
         threg2=min(threg2,accur)
         if(threg2.lt.1.0d-12) threg2=1.0d-12
c
         write(iout,153) threg2,thres1,thres2
  153 format(' Sharpen integral thresholds = ',3(1pe11.4,2x))
c
         if(xlow.lt.1.0d-7)then
            thrs=thrs*5.d0   ! looser cphf threshold
            write(iout,154) thrs
  154 format('  Loosen cphf conv threshold = ',1pe11.4,2x)
         endif
      endif
c--------------------------------------------------------------------
      end
c======================================================================
      subroutine make_d2hesX(bl,inx,nocc,ncf,vec,val)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension inx(12,*)
      dimension vec(ncf,ncf),val(ncf)
c
      call getival('ncs ',ncs)
c
      call getmem(ncs*ncs,lden2)
      call d2screen(inx,ncs,ncf,nocc,vec,val, bl(lden2),bl)
c
      call save1mat(69, 1 ,ncs*ncs ,bl(lden2))
      call retmem(1)
c
      end
c======================================================================
      subroutine make_d2hess(rhf,bl,inx,ncf,ncs)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      dimension bl(*)
      dimension inx(12,*)
c
      call getival('lvec',lvec)
      call getival('lval',lval)
      call getival('nocc',nocc)
c
      call getmem(ncs*ncs,lden2)
      call d2screen(inx,ncs,ncf,nocc,bl(lvec),bl(lval),bl(lden2),bl)
c
      if(.not.rhf)then
        call getival('lvecB',mvec)
        call getival('lvalB',mval)
        call getival('noccB',mocc)
        call getmem(ncs*ncs,ldem2)
        call d2screen(inx,ncs,ncf,mocc,bl(mvec),bl(mval),bl(ldem2),bl)
        call add2vec(ncs*ncs,bl(lden2),bl(ldem2),bl(lden2))
        call retmem(1)
      endif
c
      call save1mat(69, 1 ,ncs*ncs ,bl(lden2))
      call retmem(1)
c
      end
c======================================================================
      subroutine set_Natoms(natom)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,
     1              ncs,nsy(4),nsym,nganz(35),lopt(30)
c
      call setival('na',natom)
      na=natom
c
      end
c======================================================================

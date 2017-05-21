c======================================================================
c May 2001 KW :
c
c        the main program for analytical dipole polarizabilites
c
c December 2007 KW :
C  add hyperpolarizability and made it working with DFT
c======================================================================
c Formulae :
c
c    Axy = d2E(F)/dFxdFy
c
c    where F=(Fx,Fy,Fz) an external electric field
c
c
c Perturbation theory :
c
c     H=H0 + F*r  ; r=(x,y,z)
c
c               x y       y x
c     Axy = Tr H D  = Tr H D
c
c   x
c  Hij =  (ix|h0|j) + (i|x|j) + (i|h0|jx)
c
c                             T   y       y
c   y            y           Co*[F  - eo*S ]*Cv       T       T
c  D   =  0.5*D*S D + SUMov{ -------------------[ Co*Cv + Cv*Co ] }
c                                   eo - ev
c
c        y    y     y          y
c where F  = H + G(D,g) + G(D,g )
c
c
c Above formulae are general, for field dependent basis sets.
c======================================================================
c
      subroutine polariz

      use memory

      implicit real*8 (a-h,o-z)
      character jobname*256,cdum*20,group*4
      Logical rhf
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icon,itest,npl(9)
c
      character*80 wvfnc
c
c----------------------------------------------------------------------
      PARAMETER (IUnit=1)               ! unit number for checkpoint I/O
      data thrsh/1.0d-6/                ! zero threshold
c
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
c-----------------------------------------------------------------------
c get original values of symmetry pair pointers
c
        call getival('SymFunPr',ifp)
        call getival('SymShPr',ifp1)
c-----------------------------------------------------------------------
c Re-store commons big, intbl, and ener
c
      ncall=1
      call readbl(ncall,ictr)
c
c output : ictr = new address in inx (texas95)
c-----------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c---------------------------------------------------------------
c  read from <control> file
c  total number of atomic centres, including dummies
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
c
      call fdcntrl(iunit,9,'$wavefunc',iend)
      if(iend.eq.0)call rdcntrl(IUnit,9,'$wavefunc',3,idum,dum,wvfnc)
      CALL RdDFT(IUnit,idft)
c
c     dummy atoms
c
      ndum1=0     ! number of charged dummies
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
       NAtom = NAtoms-Ndum1-Ndum2
      Else
       NAtom = NAtoms
      EndIf
c
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------
c save number of real atoms for use in DFT
      call setival('realat',natom)
c---------------------------------------------------------------
c -- check there are electrons!
      If(NAlpha.EQ.0) Call nerror(1,'POLAR Module',
     $                     'There Are No Electrons!',0,0)
c
c -- set logical flag "rhf"
c
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
C
c-----------------------------------------------------------------------
c preper symmetry data for dft 
c
        if(idft.ne.0) call symm4dft(natom)
c-----------------------------------------------------------------------
c Read-in the polarizability options
c
      call pola_options(wvfnc,fact0)
c-----------------------------------------------------------------------
c Check if any change in grid factor
      If(idft.GT.0.AND.fact0.NE.0) Then
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $        FORM='FORMATTED',STATUS='OLD')
        call rdcntrl(IUnit,7,'$factor',2,idum,factor,cdum)
c -- write new grid factor on the control file
        call wrcntrl(IUnit,7,'$factor',2,idum,fact0,cdum)
c -- if new grid factor LESS than used in SCF print warning
        if(fact0.LT.factor) call message('POLAR',
     $        '**WARNING** DFT grid quality REDUCED from SCF',0,0)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
      EndIf
c-----------------------------------------------------------------------
C
C  Currently this code can only handle closed-shell RHF and RDFT
C
      If(wvfnc(1:3).NE.'RHF'.AND.wvfnc(1:4).NE.'RDFT') Then
         call nerror(2,'POLAR module',
     $  'Can only handle RHF or RDFT - use numerical module',
     $   0,0)
      EndIf
c
c-----------------------------------------------------------------------
c open files for Polarizability calculations :
c
      call open4nmr(idft)
c-----------------------------------------------------------------------
c memory checking
c
      call pola_memcheck(rhf)
c---------------------------------------------------------------
c put down a memory mark
      call mmark
c---------------------------------------------------------------
c allocate memory for :
c 1. dipol moment (3)
c 2. polarizability tensor polar(3,3)
c 3. hyperpolarizability tensor hypol(10)
c and zero them  out
c
c
      call getmem(3,ldipol)
      call getmem(9,lpolar)
      call getmem(10,lhypol) ! 10 independent components
C
      call zeroit(bl(ldipol),3)
      call zeroit(bl(lpolar),9)
      call zeroit(bl(lhypol),10)
c---------------------------------------------------------------
c read in the density as well as eigenvectors & eigenvalues
c
      call read_deig(rhf,bl)
c---------------------------------------------------------------
c check if the basis set contains L-shells
c and if it does then make S,P partitioning
c
c??   call trans_l_2_sp(rhf,bl,inx(ictr),ictr,lshell,'polariz')
c
c ictr is changed on return if l-shells are present
c---------------------------------------------------------------
c
c THE ONE-ELECTRON PART OF THE POLARIZABILITY CALCULATIONS:
c
c allocate memory needed in 1-el. part of polariz
c for Hx,Hy and Hz matrices (ntri,3) (dipole integrals)
c These are needed for the final polar. Tr(Hx*Dy) and
c for CPHF as well.
c If field-dependent basis set is used then for CPHF
c we need also (ix|h0|j)+(i|h0|jx) and S1ij=(i1|j)+(i|j1)
c
c So allocate Fxyz for CPHF and Hxyz for the final polariz.tensor
c For field-independent basis set Fxyz=Hxyz + G(Dxyz,g0)
c
      call getmem(3*ntri,lhxyz)
      call getmem(3*ntri,lfxyz)
      call getmem(3*ntri,lsxyz)
c
      call setival('lhxyz',lhxyz)
      call setival('lfxyz',lfxyz)
      call setival('lsxyz',lsxyz)
c
      call zeroit(bl(lhxyz),3*ntri)
      call zeroit(bl(lfxyz),3*ntri)
      call zeroit(bl(lsxyz),3*ntri)
c---------------------------------------------------------------
c
      CALL POLAR1(rhf,bl,bl(ictr),ntri,bl(lsxyz),bl(lhxyz),bl(lfxyz))
c---------------------------------------------------------------
c
c THE TWO-ELECTRON PART OF THE POLARIZABILITY CALCULATIONS:
c
c it supposed to calculate the G(D0,g1) part of the fock1 matrices
c and it is needed only if field-dependent basis set is used
c
c--------------------------------------------------------
c this call is needed even for field-idependent basis set
c because it initializes the 2-el. program used in CPHF
c
      call POLAR2(idft,ax,rhf,bl,bl(ictr))
c
c--------------------------------------------------------
c
c  THE CPHF PART
c
c allocate memory for Ist-order density matrices :
c
      call getmem(ntri*3,ldxyz)
c     call setival('ldxyz',ldxyz)
c
c solve the Coupled Perturbed Hartree-Fock equestion
c for an external electric field
c
c calculate the first-order density matrices
c
      lsemi=0
      call get_d1xyz(idft,ax,rhf,bl,bl(ictr),ncf,ntri,
     *               bl(lsxyz), bl(lfxyz),bl(ldxyz),lsemi)
c
c--------------------------------------------------------
c lfxyz contains only 2el-part, add 1el part :
c
      call add2vec(3*ntri,bl(lhxyz),bl(lfxyz),bl(lfxyz))
c
c--------------------------------------------------------
c At this point we are ready for the final 
c dipole, polarizability and hyperpolarizability :
c
      call getival('ldensi',lden)
      call dipo_moment(bl(lhxyz),bl(lden),ncf,ntri,bl(ldipol),dipol)
c
      call pola_tensor(bl(lhxyz),bl(ldxyz),ncf,ntri,bl(lpolar),Am,Aa)
c
      call hypol_tensor(bl(lfxyz),bl(ldxyz),ncf,ntri,bl(lhypol),HPm,
     *                  wvfnc)
c
c--------------------------------------------------------
c zero out tensor components which are below threshold
c
       call zerobelowt( 3,thrsh,bl(ldipol))
       call zerobelowt( 9,thrsh,bl(lpolar))
       call zerobelowt(10,thrsh,bl(lhypol))
c--------------------------------------------------------
c print out the dipol moment vector
c
       call printDip(iout,icon,bl(ldipol),dipol)
c--------------------------------------------------------
c print out the polarizability and hyperpolarizability 
c
       call printPol(iout,icon,bl(lpolar),am,aa)
c--------------------------------------------------------
c print out the hyperpolarizability tensor
c
       call printHypol(iout,icon,bl(lhypol),hpm)
c--------------------------------------------------------
c-------------------------------------------
c end of the story for polarizablity
c-------------------------------------------
c transfer back original L-shell basis set info
c
c     call trans_sp_2_l(bl,inx(ictr),ictr,lshell)
c
c returns the original value of ictr
c-------------------------------------------
c return all memory allocated in this routine
c
      call retmark
c-------------------------------------------
c restore original values of symmetry pair and atom pointers
c
      call setival('SymFunPr',ifp)
      call setival('SymShPr',ifp1)
c
c memory status :
c
      call memory_status(' end of polarizab.')
c-------------------------------------------
      call clos4nmr()
c
      end
c======================================================================
      subroutine pola_options(wvfnc,fact0)
c
c read in polarizability options
c

      use memory

      implicit real*8 (a-h,o-z)
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
c---------------------------------------------------------------------
      character*80 wvfnc
c---------------------------------------------------------------------
      parameter (nopt=11)
      character*4 word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      character*256 chopv(nopt)
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp/1,    11,     12,   11,     1,
     *            1,     0,     1,     3,    11,
     *            11/
c
      data word /'prin','thr1','thr2','thre','iter',
     *           'ncac','noac','rese','limi','lvsh',
     *           'grid'/
c
c---------------------------------------------------------------------
      maxprice=maxpricx
c---------------------------------------------------------------------
c thre1=1.0d-10   ! one-el. integral threshold
c thre2=1.0d-10   ! two-el. integral threshold and
c                   final two_el. integral threshold in direct-CPHF
c thre3=1.0d-8      loose two-el. integral threshold in direct-CPHF
c
c thre =1.d0-4     threshold for CPHF convergance
c
c iter = 20       maximum number of CPHF iterations
c
c idelt=0 or 1 ; do not use or use delta D10 in CPHF
c
c xlvsh=0.0d0  ! level shift in CPHF, added to all energy differences
c---------------------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
      write(iout,*)
     *  '         The Analytical Module for Electrical Properies'
      write(iout,*)' '
      write(iout,1155) wvfnc(1:10)
 1155 format('                               ',a10)
      write(iout,*)' '
c
      call f_lush(iout)
c-----------------------------------------------------------
c
c default values :
c
      nopola=1            ! no field-dependent basis set
      nprint=0
      thre1=1.0d-10
      thre2=1.0d-10
      thre3=1.0d-8
c
      thres=1.0d-3
      maxiter=30
      ncache=ncachx
      noacce=0
c2002 idelt =1              ! delta density ALWAYS used in new cphf
      ireset=20             ! reset CPHF after iter=20
      limxmem=2 000 000
c2002 limblks=300
      limblks=0
      limpair=100
      xlvsh=0.0d0
c
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
c
      if(ifound(1).gt.0) nprint=iopv(1,1)
      if(ifound(2).gt.0) thre1=10.d0**(-ropv(1,2))
      if(ifound(3).gt.0) thre2=10.d0**(-ropv(1,3))
      if(ifound(3).gt.0) thre3=10.d0**(-ropv(2,3))
      if(ifound(4).gt.0) thres=10.d0**(-ropv(1,4))
      if(ifound(5).gt.0) maxiter =iopv(1,5)
      if(ifound(6).gt.0) ncache  =iopv(1,6)
      if(ifound(7).gt.0) noacce  = 1
      if(ifound(8).gt.0) ireset  =iopv(1,8)
      if(ifound(9).gt.0) then
         limxmem=iopv(1,9)
         limblks=iopv(2,9)
         limpair=iopv(3,9)
      endif
      if(ifound(10).gt.0) xlvsh=ropv(1,10)
      fact0=0.0d0           ! grid factor
      if(ifound(11).gt.0) fact0=ropv(1,11)
c---------------------------------------------------------------------
c check for linear dependencies in basis set :
c
      call tstrval('xlows',lows)
      if(lows.eq.1) then
         call getrval('xlows',xlow)  ! lowest eigenvalue of S
      else
         xlow=0.1d0                  ! when running without scf
      endif
      if(xlow .lt. 5.0d-6 ) then
         if((ifound(3).gt.0)) then 
            write(iout,*)
     *   '              Near Linear Dependency in Basis Set'
            write(iout,*)
     *   '              integral threshold will be sharpened'
            write(iout,*)
     *   '                if not redefined in the input    '
         else
            call check4ld(thre2,thre3,thres)
         endif
      endif
c---------------------------------------------------------------------
      call setival('nopola',nopola)
      call setival('print',nprint)
      call setrval('thr1',thre1)
      call setrval('thr2',thre2)
      call setrval('thr3',thre3)
      call setrval('thres',thres)
      call setival('maxit',maxiter)
      call setival('noacc',noacce)
c
      if(ireset.gt.30) ireset=30
c
      call setival('reset',ireset)
      call setrval('xlvsh',xlvsh)
c---------------------------------------------------------------------
      lopt(1)=nprint
c---------------------------------------------------------------------
c print options :
c
c     write(iout,200)
      do 100 iop=1,nopt
       if(ifound(iop).gt.0) then
          if(iop.eq. 1) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 2) write(iout,220) word(iop),ropv(1,iop)
          if(iop.eq. 3) write(iout,220) word(iop),(ropv(ii,iop),ii=1,2)
          if(iop.eq. 4) write(iout,220) word(iop),ropv(1,iop)
          if(iop.eq. 5) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 6) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 7) write(iout,210) word(iop),1
          if(iop.eq. 8) write(iout,210) word(iop),iopv(1,iop)
          if(iop.eq. 9) write(iout,210) word(iop),(iopv(ii,iop),ii=1,3)
          if(iop.eq.10) write(iout,220) word(iop),ropv(1,iop)
          if(iop.eq.11) write(iout,220) word(iop),ropv(1,iop)
       endif
  100 continue
c
  200 format(72('-'))
  210 format(' POLA option = ',a4,'   is on ; value = ',3(i10,2x))
  220 format(' POLA option = ',a4,'   is on ; value = ',3(f10.5,1x))
c 220 format(' POLA option = ',a4,'   is on ; value = ',3(e12.5,1x))
c---------------------------------------------------------------------
      write(iout,200)
      write(iout,*)' '
c---------------------------------------------------------------------
      call f_lush(iout)
c
      end
c======================================================================
      subroutine pola_memcheck(rhf)

      use memory

      implicit real*8 (a-h,o-z)
      Logical rhf
c---------------------------------------------------------------------
c
c  Total memory needed for polarizability  :
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
         call nerror(3,'POLAR module',
     *              'Memory: needed and available ',
     *               needed-ioffset,lcore-ioffset)
      endif
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) then
        write(iout,209) needed,lcore
        return
      endif
c
  209 format (/1x,' Memory in BL for POLARIZ   calculations '/
     *            '  needed ',i10,'  available ',i10/)
  210 format (/1x,'common bl too small for hess  run, required =',i10,
     *3x,' available =',i10/)
c----------------------------
      end
c======================================================================
c
      subroutine read_deig(rhf,bl)
c
c  reads in density & eigenvec & eigenval.
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension bl(*)
c
      ntri=(ncf*(ncf+1))/2
c
c----------------------------------------------------
c allocate memory for, and read in, density & ovelap matrices
c
      np1 = 1
      np4 = 4
      call getmem(ntri,lden)
      call getmem(ntri,lsmat)
      call getmem(ntri,lfock)
c
      call setival('lfock',lfock)
      call setival('lsmat',lsmat)
      call setival('ldensi',lden)
c
coverlap
      call sudat(np1,'s matrix',ni)
      if(ni.gt.0) then
        call rea(bl(lsmat),ntri,np1,'s matrix')
      else
        call restart(np1,0)
        call sudat(np1,'s matrix',ni)
        if(ni.gt.0) then
           call rea(bl(lsmat),ntri,np1,'s matrix')
        else
          call nerror(4,'POLAR module',
     1   'Cannot find the 0th-order overlap on <jobname>.14',
     2    ni,ni)
        endif
      endif
c
      IF(.not.rhf) Then
         call getmem(ntri,mden)
         call setival('ldensB',mden)
c
         call sudat(np4,'dena_uhf',ni)
         if(ni.gt.0) then
           call rea(bl(lden),ntri,np4,'dena_uhf')
           call rea(bl(mden),ntri,np4,'denb_uhf')
         else
           call restart(np4,0)
           call restart(np1,0)
           call sudat(np4,'dena_uhf',ni)
           if (ni.gt.0) then
              call rea(bl(lden),ntri,np4,'dena_uhf')
              call rea(bl(mden),ntri,np4,'denb_uhf')
           else
             call nerror(5,'POLAR module',
     1      'Cannot find the 0th-order density on <jobname>.14',
     2       ni,ni)
           endif
         endif
      ElSE
         call sudat(np4,'den0_rhf',ni)
         if(ni.gt.0) then
           call rea (bl(lden),ntri,np4,'den0_rhf')
           call rea (bl(lfock),ntri,np4,'fock_rhf')
         else
           call restart(np4,0)
           call restart(np1,0)
           call sudat(np4,'den0_rhf',ni)
           if (ni.gt.0) then
             call rea (bl(lden),ntri,np4,'den0_rhf')
             call rea (bl(lfock),ntri,np4,'fock_rhf')
           else
             call nerror(5,'POLAR module',
     1      'Cannot find the 0th-oder density on <jobname>.14',
     2       ni,ni)
           endif
         endif
      ENDIF
c----------------------------------------------------
c allocate memory for eigen and read in:
c
      call getmem(ncf*ncf,lvec)      ! eigenvectors
      call getmem(ncf    ,lval)      ! eigenvalues
      call getmem(2      ,locc)      ! no of occupied orbitals
c
      call setival('lvec',lvec)
      call setival('lval',lval)
c----------------------------------------------------
c read-in the eigenvectors & eigenvalues :
c
      np4=4
      If(rhf) Then
c --     read in closed-shell MOs and orbital energies
         call rea (bl(lvec),ncf*ncf,np4,'evec_rhf')
         call rea (bl(lval),ncf    ,np4,'eval_rhf')
         call rea (bl(locc),1      ,np4,'nocc_rhf')
         nocc=bl(locc)
         call setival('nocc',nocc)
      Else
c --     read in alpha MOs and orbital energies
         call rea (bl(lvec),ncf*ncf,np4,'evea_uhf')
         call rea (bl(lval),ncf    ,np4,'evaa_uhf')
         call rea (bl(locc),2      ,np4,'nocc_uhf')
         NAlpha = bl(locc)
         NBeta  = bl(locc+1)
c --     now read in beta MOs and orbital energies
         call rea (bl(lvec),ncf*ncf,np4,'eveb_uhf')
         call rea (bl(lval),ncf    ,np4,'evab_uhf')
      EndIf
c----------------------------------------------------
      end
c======================================================================
      subroutine dipo_moment(hxyz,dens,ncf,ntri,dipole,dipolt)
      use memory
      implicit real*8 (a-h,o-z)
      dimension hxyz(ntri,3),dens(ntri),dipole(3)
c
      call spur(hxyz(1,1),dens,ncf,dipole(1))
      call spur(hxyz(1,2),dens,ncf,dipole(2))
      call spur(hxyz(1,3),dens,ncf,dipole(3))
c
      dipole(1)=-dipole(1)
      dipole(2)=-dipole(2)
      dipole(3)=-dipole(3)
c
      call getival('inuc',inuc)
      call getival('na',na)
      call nucdipole(na,dipole,bl(inuc))
c
      dipolt=dipole(1)**2+dipole(2)**2+dipole(3)**2
      dipolt=sqrt(dipolt)
c
      end
c=======================================================================
c
      subroutine pola_tensor(hxyz,dxyz,ncf,ntri,Pol,Amean,Aaniz)
      implicit real*8 (a-h,o-z)
      dimension hxyz(ntri,3),dxyz(ntri,3),Pol(3,3)
c
      call spur(hxyz(1,1),dxyz(1,1),ncf,Pol(1,1))
      call spur(hxyz(1,1),dxyz(1,2),ncf,Pol(1,2))
      call spur(hxyz(1,1),dxyz(1,3),ncf,Pol(1,3))
c
      call spur(hxyz(1,2),dxyz(1,1),ncf,Pol(2,1))
      call spur(hxyz(1,2),dxyz(1,2),ncf,Pol(2,2))
      call spur(hxyz(1,2),dxyz(1,3),ncf,Pol(2,3))
c
      call spur(hxyz(1,3),dxyz(1,1),ncf,Pol(3,1))
      call spur(hxyz(1,3),dxyz(1,2),ncf,Pol(3,2))
      call spur(hxyz(1,3),dxyz(1,3),ncf,Pol(3,3))
c
c what we have now is 2*d2E/dFidFj , polij=-2*d2E/dFidFj
c change the sign
c
      call vscal(9,-1.0d0,Pol)
c
c Experimentally determined polarizability invariants are
c (P. Calaminici, K.Jug, A.M.Koster, J.C.P., 109,(18),7756 (1998)
c
c
c (1) mean polarizability : alpha=1/3(Axx+Ayy+Azz)
c (2) polarizabilty anisotropy: 1/2[(Axx-Ayy)^2+(Axx-Azz)^2+(Ayy-Azz)^2]
c
      amean=0.3333d0*(pol(1,1)+pol(2,2)+pol(3,3))
      aaniz=      (pol(1,1)-pol(2,2))**2
      aaniz=aaniz+(pol(1,1)-pol(3,3))**2
      aaniz=aaniz+(pol(2,2)-pol(3,3))**2
      aaniz=0.5d0*aaniz
c   
c
      return
      end
c======================================================================
c
      subroutine hypol_tensor(fxyz,dxyz,ncf,ntri,hypol,HPmean,wvfnc)
      use memory
      implicit real*8 (a-h,o-z)
      character*80 wvfnc
      dimension fxyz(ntri,3),dxyz(ntri,3),hypol(10)
c----------------------------------------------------------------------
c Input :
c
c fxyz - Fx, Fy and Fz matrices (1st-order Fock matrices)
c dxyz - Dx, Dy and Dz matrices (I-st order density)
c ncf  - basis set dimension
c ntri - ncf*(ncf+1)/2
c
c Output :
c
c hypol  : hyperpolarizability tensor (10 independent components)
c
c           1 xxx
c           2 yxx    3 yyx   4 yyy
c           5 zxx    6 zyx   7 zyy
c           8 zzx    9 zzy  10 zzz
c----------------------------------------------------------------------
c THE (Ist-ORDER) HYPERPOLARIZABILITY in ORTHOGONAL BASIS SET
c
c beta=d3E/dFxdFydFz ; wrt components of an ext. electric field
c
c Hartree-Fock :
c
c Bxyz= -4*Tr[Fx*(YT*Z-Z*YT)] -4*Tr[Fy*(ZT*X-X*ZT)] -4*Tr[Fz*(XT*Y-Y*ZT)]
c                  x y z               y z x                 z x y
c
c where X,Y,Z are projected Dx,Dy,Dz matrices :
c
c X=P1*Dx*P2  ,  Y=P1*Dy*P2  ,  Z=P1*Dz*P2
c
c with projectors P1=1/2*D0   , P2=1-1/2D0  (P2=1-P1)
c----------------------------------------------------------------------
c The above formulas for NON-ORTHOGONAL BASIS SET become 
c
c Bxyz= -4*Tr[Fx*(YT*S*Z-Z*S*YT)] 
c       -4*Tr[Fy*(ZT*S*X-X*S*ZT)] 
c       -4*Tr[Fz*(XT*S*Y-Y*S*ZT)]
c
c where X,Y,Z are projected Dx,Dy,Dz matrices :
c
c X=P1*Dx*P2  ,  Y=P1*Dy*P2  ,  Z=P1*Dz*P2
c
c with projectors P1=1/2*D0*S   , P2=1-1/2D0*S  (P2=1-P1)
c----------------------------------------------------------------------
c
c Independent Beta tensor components :
c
c XXX, YXX, YYX, YYY, ZXX, ZYX, ZYY, ZZX, ZZY, ZZZ  ( 10 of them )
c
c  hp_xxx=-12*x_xx                  = -4 [Fx*XX+Fx*XX+Fx*XX]
c  hp_yxx= -4*y_xx -4*x_xy -4*x_yx  = -4 [Fy*XX+Fx*XY+Fx*YX]
c  hp_yyx= -4*y_yx -4*y_xy -4*x_yy  = -4 [Fy*YX+Fy*XY+Fx*YY]
c  hp_yyy=-12*y_yy
c  hp_zxx= -4*z_xx -4*x_xz -4*x_zx  = -4 [Fz*XX+Fx*XZ+Fx*ZX]
c  hp_zyx= -4*z_yx -4*y_xz -4*x_zy  = -4 [Fz*YX+Fy*XZ+Fx*ZY]
c  hp_zyy= -4*z_yy -4*y_yz -4*y_zy  = -4 [Fz*YY+Fy*YZ+Fy*ZY]
c  hp_zzx= -4*z_zx -4*z_xz -4*x_zz  = -4 [Fz*ZX+Fz*XZ+Fx*ZZ]
c  hp_zzy= -4*z_zy -4*z_yz -4*y_zz  = -4 [Fz*ZY+Fz*YZ+Fy*ZZ]
c  hp_zzz=-12*z_zz
c
c where  "projected density" product combinations are :
c
c  XX=XT*X-X*XT   XY=XT*Y-Y*XT  XZ=XT*Z-Z*XT
c  YX=YT*X-X*YT   YY=YT*Y-Y*YT  YZ=YT*Z-Z*YT
c  ZX=ZT*X-X*ZT   ZY=ZT*Y-Y*ZT  ZZ=ZT*Z-Z*ZT
c
c The products needed :
c
c  XT*X and X*XT           XT*Y=(YT*X)T and Y*XT  XT*Z=(ZT*X)T and Z*XT
c  YT*X and X*YT=(Y*XT)T   YT*Y         and Y*YT  YT*Z=(ZT*Y)T and Z*YT
c  ZT*X and X*ZT=(Z*XT)T   ZT*Y=(YT*Z)T and Y*ZT  ZT*Z         and Z*ZT
c
c
c Express the upper triangle by the lower one using (A*B)T=BT*AT
c 
c
c  XX=XT*X-X*XT    XY=(YT*X)T-(X*YT)T   XZ=(ZT*X)T-(X*ZT)T
c  YX=YT*X-X*YT    YY=YT*Y   -   Y*YT   YZ=(ZT*Y)T-(Y*ZT)T
c  ZX=ZT*X-X*ZT    ZY=ZT*Y   -   Y*ZT   ZZ=ZT*Z   -   Z*ZT
c
c  Using (A+B)T=AT+BT we have e.g. :
c
c  XY=XT*Y-Y*XT=(YT*X)T-(X*YT)T=(YT*X-X*YT)T=(YX)T
c  XZ=                                      =(ZX)T
c  YZ=YT*Z-Z*YT=(ZT*Y)T-(Y*ZT)T=(ZT*Y-Y*ZT)T=(ZY)T
c
c All we need is :
c
c  XX  (YX)T  (ZX)T
c  YX    YY   (ZY)T
c  ZX    ZY     ZZ
c
c  6 product combinations needed 
c
c  XX=XT*X-X*XT
c  YX=YT*X-X*YT    YY=YT*Y-Y*YT
c  ZX=ZT*X-X*ZT    ZY=ZT*Y-Y*ZT   ZZ=ZT*Z-Z*ZT
c
c F1      Prod.Comb.   contributes to Hypol
c------------------------------
c Fx  *      XX          XXX
c Fy  *      XX          YXX
c Fz  *      XX          ZXX
c
c Fx  *      YX+         YXX
c Fx  *      YX          YXX
c Fy  *      YX          YYX
c Fy  *      YX+         YYX
c Fz  *      YX          ZYX
c
c Fy  *      YY          YYY
c Fx  *      YY          YYX
c Fz  *      YY          ZYY
c
c Fx  *      ZX          ZXX
c Fx  *      ZX+         ZXX
c Fy  *      ZX+         ZYX
c Fz  *      ZX          ZZX
c Fz  *      ZX+         ZZX
c
c Fx         ZY          ZYX
c Fy  *      ZY+         ZYY
c Fy         ZY          ZYY
c Fz  *      ZY          ZZY
c Fz  *      ZY+         ZZY
c
c Fz  *      ZZ          ZZZ
c Fx  *      ZZ          ZZX
c Fy  *      ZZ          ZZY
c----------------------------------------------------------------------
c     call getmem(0,last1)
c----------------------------------------------------------------------
c  allocate memory
c
      call matdef('FmaX','s',ncf,ncf)   ! Fx
      call matdef('FmaY','s',ncf,ncf)   ! Fy
      call matdef('FmaZ','s',ncf,ncf)   ! Fz
c
      call matdef('DenX','s',ncf,ncf)
      call matdef('DenY','s',ncf,ncf)
      call matdef('DenZ','s',ncf,ncf)
c
c matrix addresses
c
      ifx=mataddr('FmaX')
      ify=mataddr('FmaY')
      ifz=mataddr('FmaZ')
c
      idx=mataddr('DenX')
      idy=mataddr('DenY')
      idz=mataddr('DenZ')
c
c----------------------------------------------------------------------
c copy H1 and D1 matrices into local storage :
c
      call dcopy(ntri,fxyz(1,1),1,bl(ifx),1)
      call dcopy(ntri,fxyz(1,2),1,bl(ify),1)
      call dcopy(ntri,fxyz(1,3),1,bl(ifz),1)
c
      call dcopy(ntri,dxyz(1,1),1,bl(idx),1)
      call dcopy(ntri,dxyz(1,2),1,bl(idy),1)
      call dcopy(ntri,dxyz(1,3),1,bl(idz),1)
c----------------------------------------------------------------------
c MAKE PROJECTORS  P2=1-1/2D0*S0  and P1=1-P2=1/2D0*S0
c
c allocate memory for projected densities
c 
      call matdef('PenX','q',ncf,ncf)
      call matdef('PenY','q',ncf,ncf)
      call matdef('PenZ','q',ncf,ncf)
c
c...................................................
      call matdef('pro1','q',ncf,ncf) !  P1=  1/2D0*S0
      call matdef('pro2','q',ncf,ncf) !  P2=1-1/2D0*S0
c...................................................
      ip1=mataddr('pro1')
      ip2=mataddr('pro2')
c
      call matdef('sma0','s',ncf,ncf)   ! zero-order overlap
      call matdef('den0','s',ncf,ncf)   ! zero-order density
      np1=1
      np4=4
      call matread('sma0',np1,'s matrix')
      if(wvfnc(1:4).eq.'RMP2') then
         call matread('den0',np4,'dens_mp2')
      else
         call matread('den0',np4,'den0_rhf')
      endif
c
ctest....test.......test.......test.......test.......
c
c idempotency test beginning
c
      call getival('print',nprint)
      if(nprint.ge.3) then
         call matdef('tes1','q',ncf,ncf)
         call matdef('D-1/2DSD','s',ncf,ncf)
         call matzero('tes1')
         call matzero('D-1/2DSD')
         call matmmul2('den0','sma0','tes1','n','n','n')
         call matmmul2('tes1','den0','D-1/2DSD','n','n','n')
c
         call mattrace('D-1/2DSD',trace1)
         call mattrace('den0'    ,trace2)
         write(6,*)' trace of DSD=',trace1
         write(6,*)' trace of 2D =',trace2*2.0d0
         trace2=2.0d0*trace2
         factor=trace1/trace2
         write(6,*)' factor=',factor
c
cccccc   call matscal('den0',factor)

c
         call matadd1('den0',-2.0d0,'D-1/2DSD')
         call matprint('D-1/2DSD',6)
         call mattrace('D-1/2DSD',trace1)
         write(6,*)'trace of D-1/2DSD =',trace1
c
         call matrem('D-1/2DSD')
         call matrem('tes1')
      endif
c
c idempotency test end
c
ctest....test.......test.......test.......test.......
c
c----------------------------------------------------------------------
c make projectors P1 and P2
c
c     call matprint('sma0',6)
c     call matprint('den0',6)
c
      call make_proj('den0','sma0','pro2')   ! from cmp2
c
      call matrem('den0')
      call matrem('sma0')
c
c
c make projector P1=1/2D0*S0   (p1=1-p2)
c
      call make_pro1(bl(ip2),ncf,bl(ip1))
c
c     call matprint('pro1',6)
c     call matprint('pro2',6)
c----------------------------------------------------------------------
c do projection on Dx, Dy, Dz 
c
c projections according to P1*DX*P2+
c
      call matdef('temp','q',ncf,ncf)
c Dx
      call matmmul2('pro1','DenX','temp','n','n','n')
      call matmmul2('temp','pro2','PenX','n','t','n')
c Dy
      call matmmul2('pro1','DenY','temp','n','n','n')
      call matmmul2('temp','pro2','PenY','n','t','n')
c Dz
      call matmmul2('pro1','DenZ','temp','n','n','n')
      call matmmul2('temp','pro2','PenZ','n','t','n')
c
      call matrem('temp')
c
c----------------------------------------------------------------------
c if MP2 then D1 is not =p1 D1 p2 + (p1 D1 p2)T
c but plus  +p1D1p1 +p2D1p2
c
      GO TO 4321
c
      if(wvfnc(1:4).eq.'RMP2') then
         call matdef('temp','q',ncf,ncf)
         call matdef('p1dp1','q',ncf,ncf)
         call matdef('p2dp2','q',ncf,ncf)
c Dx
         call matmmul2('pro1','DenX','temp','n','n','n')
         call matmmul2('temp','pro1','p1dp1','n','t','n')
         call matmmul2('pro2','DenX','temp','n','n','n')
         call matmmul2('temp','pro2','p2dp2','n','t','n')
         call matadd('p1dp1','p2dp2')
         call matscal('p2dp2',0.5d0)
         call matadd('p2dp2','PenX') !P1=p1D1p2+p2D1p1+0.5(p1D1p1+p2D1p2)
c Dy
         call matmmul2('pro1','DenY','temp','n','n','n')
         call matmmul2('temp','pro1','p1dp1','n','t','n')
         call matmmul2('pro2','DenY','temp','n','n','n')
         call matmmul2('temp','pro2','p2dp2','n','t','n')
         call matadd('p1dp1','p2dp2')
         call matscal('p2dp2',0.5d0)
         call matadd('p2dp2','PenY') !P1=p1D1p2+p2D1p1+0.5(p1D1p1+p2D1p2)
c Dz
         call matmmul2('pro1','DenZ','temp','n','n','n')
         call matmmul2('temp','pro1','p1dp1','n','t','n')
         call matmmul2('pro2','DenZ','temp','n','n','n')
         call matmmul2('temp','pro2','p2dp2','n','t','n')
         call matadd('p1dp1','p2dp2')
         call matscal('p2dp2',0.5d0)
         call matadd('p2dp2','PenZ') !P1=p1D1p2+p2D1p1+0.5(p1D1p1+p2D1p2)
c
         call matrem('p2dp2')
         call matrem('p1dp1')
         call matrem('temp')
      endif
 4321 CONTINUE !      GO TO 4321
c----------------------------------------------------------------------
ctest....test.......test.......test.......test.......
c
      if(nprint.ge.3) then
      call matdef('test','q',ncf,ncf)
      call matdef('tess','s',ncf,ncf)
c
      call matprint('DenX',6)
      call matcopy('PenX','test')
      call matpose('test')
      call matadd('PenX','test')
      call matcopy('test','tess')
      call matprint('tess',6)
      call matadd1('DenX',-1.0d0,'tess')
      call matprint('tess',6)
c
      call matprint('DenY',6)
      call matcopy('PenY','test')
      call matpose('test')
      call matadd('PenY','test')
      call matcopy('test','tess')
      call matprint('tess',6)
      call matadd1('DenY',-1.0d0,'tess')
      call matprint('tess',6)
c
      call matprint('DenZ',6)
      call matcopy('PenZ','test')
      call matpose('test')
      call matadd('PenZ','test')
      call matcopy('test','tess')
      call matprint('tess',6)
      call matadd1('DenZ',-1.0d0,'tess')
      call matprint('tess',6)
c
      call matrem('tess')
      call matrem('test')
      endif
c
ctest....test......test........test...........test...
c----------------------------------------------------------------------
c
c----------------------------------------------------------------------
c Now PenX=P1*Dx*P2+ , PenY=P1*Dy*P2+ , PenZ=P1*Dz*P2+
c...................................................
      call matrem('pro2')
      call matrem('pro1')
c...................................................
c----------------------------------------------------------------------
c XXX, YXX, YYX, YYY, ZXX, ZYX, ZYY, ZZX, ZZY, ZZZ  ( 10 of them )
      hp_xxx=0.0d0
      hp_yxx=0.0d0
      hp_yyx=0.0d0
      hp_yyy=0.0d0
      hp_zxx=0.0d0
      hp_zyx=0.0d0
      hp_zyy=0.0d0
      hp_zzx=0.0d0
      hp_zzy=0.0d0
      hp_zzz=0.0d0
c----------------------------------------------------------------------
c calculate "density product" combinations
c
      call matdef('AtSB','q',ncf,ncf)
      call matdef('BSAt','q',ncf,ncf)
      call matdef('sma0','s',ncf,ncf)   ! zero-order overlap
      np1=1
      call matread('sma0',np1,'s matrix')
c
      call matdef('temp','q',ncf,ncf)
c---------------------------------------------------
c XX=XT*S*X - X*S*XT
c
      call matmmul2('PenX','sma0','temp','t','n','n')     ! Xt*S
      call matmmul2('temp','PenX','AtSB','n','n','n')     ! Xt*S*X
c
      call matmmul2('PenX','sma0','temp','n','n','n')     ! X*S
      call matmmul2('temp','PenX','BSAt','n','t','n')     ! X*S*Xt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Xt*S*X-X*S*Xt) in AtSB 
c
c contributions to hyper.pol. 
c Fx  *      XX          XXX
c Fy  *      XX          YXX
c Fz  *      XX          ZXX
c 
c calculate trace of Fx*(Xt*X-X*Xt)
c
      call matprodtr('FmaX','AtSB', hyp) ! trace of Fx*(Xt*S*X-X*S*Xt)
      hp_xxx=hp_xxx+3.0d0*hyp
c
      call matprodtr('FmaY','AtSB', hyp) ! trace of Fy*(Xt*S*X-X*S*Xt)
      hp_yxx=hp_yxx+hyp
c
      call matprodtr('FmaZ','AtSB', hyp) ! trace of Fz*(Xt*S*X-X*S*Xt)
      hp_zxx=hp_zxx+hyp
c
c---------------------------------------------------
c YX=YT*S*X - X*S*YT
c
      call matmmul2('PenY','sma0','temp','t','n','n')   ! Yt*S
      call matmmul2('temp','PenX','AtSB','n','n','n')   ! Yt*S*X
c
      call matmmul2('PenX','sma0','temp','n','n','n')   ! X*S
      call matmmul2('temp','PenY','BSAt','n','t','n')   ! X*S*Yt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Yt*S*X-X*S*Yt) in AtSB 
c
c contributions to hyper.pol. 
c Fx  *      YX+         YXX
c Fx  *      YX          YXX
c Fy  *      YX          YYX
c Fy  *      YX+         YYX
c Fz  *      YX          ZYX
c
ccc   call matprodtr('FmaX','AtSB', hyp) ! trace of FxT*(Yt*S*X-X*S*Yt)
      call matprodtr('AtSB','FmaX', hyp) ! trace of (Yt*S*X-X*S*Yt)t*Fx
      hp_yxx=hp_yxx+hyp
c
      call matprodtr('FmaX','AtSB', hyp) ! trace of FxT*(Yt*S*X-X*S*Yt)
      hp_yxx=hp_yxx+hyp
c
      call matprodtr('FmaY','AtSB', hyp) ! trace of FyT*(Yt*S*X-X*S*Yt)
      hp_yyx=hp_yyx+hyp
c
      call matprodtr('AtSB','FmaY', hyp) ! trace of (Yt*S*X-X*S*Yt)t*Fy
      hp_yyx=hp_yyx+hyp
c
      call matprodtr('FmaZ','AtSB', hyp) ! trace of FzT*(Yt*S*X-X*S*Yt)
      hp_zyx=hp_zyx+hyp
c
c---------------------------------------------------
c
c YY=YT*S*Y-Y*S*YT
c
      call matmmul2('PenY','sma0','temp','t','n','n')   ! Yt*S
      call matmmul2('temp','PenY','AtSB','n','n','n')   ! Yt*S*Y
c
      call matmmul2('PenY','sma0','temp','n','n','n')   ! Y*S
      call matmmul2('temp','PenY','BSAt','n','t','n')   ! Y*S*Yt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Yt*S*Y-Y*S*Yt) in AtSB 
c
c contributions to hyper.pol. 
c Fy  *      YY          YYY
c Fx  *      YY          YYX
c Fz  *      YY          ZYY
c
      call matprodtr('FmaY','AtSB', hyp) ! trace of Fy*(Yt*Y-Y*Yt)
      hp_yyy=hp_yyy+3.0d0*hyp
c
      call matprodtr('FmaX','AtSB', hyp) ! trace of Fx*(Yt*Y-Y*Yt)
      hp_yyx=hp_yyx+hyp
c
      call matprodtr('FmaZ','AtSB', hyp) ! trace of Fz*(Yt*Y-Y*Yt)
      hp_zyy=hp_zyy+hyp
c
c---------------------------------------------------
c
c ZX=ZT*S*X-X*S*ZT
c
      call matmmul2('PenZ','sma0','temp','t','n','n')   ! Zt*S
      call matmmul2('temp','PenX','AtSB','n','n','n')   ! Zt*S*X
c
      call matmmul2('PenX','sma0','temp','n','n','n')   ! X*S
      call matmmul2('temp','PenZ','BSAt','n','t','n')   ! X*S*Zt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Zt*S*X-X*S*Zt) in AtSB 
c
c contributions to hyper.pol. 
c Fx  *      ZX          ZXX
c Fx  *      ZX+         ZXX
c Fy  *      ZX+         ZYX
c Fz  *      ZX          ZZX
c Fz  *      ZX+         ZZX
c
c
      call matprodtr('FmaX','AtSB', hyp) ! trace of FxT*(Zt*X-X*Zt)
      hp_zxx=hp_zxx+hyp
      call matprodtr('AtSB','FmaX', hyp) ! trace of (Zt*X-X*Zt)T * Fx
      hp_zxx=hp_zxx+hyp
      call matprodtr('AtSB','FmaY', hyp) ! trace of (Zt*X-X*Zt)T * Fy
      hp_zyx=hp_zyx+hyp
      call matprodtr('FmaZ','AtSB', hyp) ! trace of FzT*(Zt*X-X*Zt)
      hp_zzx=hp_zzx+hyp
      call matprodtr('AtSB','FmaZ', hyp) ! trace of (Zt*X-X*Zt)T * Fz
      hp_zzx=hp_zzx+hyp
c
c---------------------------------------------------
c
c ZY=ZT*S*Y-Y*S*ZT
c
      call matmmul2('PenZ','sma0','temp','t','n','n')   ! Zt*S
      call matmmul2('temp','PenY','AtSB','n','n','n')   ! Zt*S*Y
c
      call matmmul2('PenY','sma0','temp','n','n','n')   ! Y*S
      call matmmul2('temp','PenZ','BSAt','n','t','n')   ! Y*S*Zt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Zt*S*Y-Y*S*Zt) in AtSB 
c
c contributions to hyper.pol. 
c Fx         ZY          ZYX
c Fy  *      ZY+         ZYY
c Fy         ZY          ZYY
c Fz  *      ZY          ZZY
c Fz  *      ZY+         ZZY
c
      call matprodtr('FmaX','AtSB', hyp) ! trace of FxT*(Zt*Y-Y*Zt)
      hp_zyx=hp_zyx+hyp
      call matprodtr('AtSB','FmaY', hyp) ! trace of (Zt*Y-Y*Zt)T * Fy
      hp_zyy=hp_zyy+hyp
      call matprodtr('FmaY','AtSB', hyp) ! trace of FyT*(Zt*Y-Y*Zt)
      hp_zyy=hp_zyy+hyp
      call matprodtr('FmaZ','AtSB', hyp) ! trace of FzT*(Zt*Y-Y*Zt)
      hp_zzy=hp_zzy+hyp
      call matprodtr('AtSB','FmaZ', hyp) ! trace of (Zt*Y-Y*Zt)T * Fz
      hp_zzy=hp_zzy+hyp
c
c---------------------------------------------------
c
c ZZ=ZT*S*Z-Z*S*ZT
c
      call matmmul2('PenZ','sma0','temp','t','n','n')   ! Zt*S
      call matmmul2('temp','PenZ','AtSB','n','n','n')   ! Zt*S*Z
c
      call matmmul2('PenZ','sma0','temp','n','n','n')   ! Z*S
      call matmmul2('temp','PenZ','BSAt','n','t','n')   ! Z*S*Zt
c
      call matadd1('BSAt',-1.0d0,'AtSB')      ! (Zt*S*Z-Z*S*Zt) in AtSB 
c
c contributions to hyper.pol. 
c Fz  *      ZZ          ZZZ
c Fx  *      ZZ          ZZX
c Fy  *      ZZ          ZZY
c
      call matprodtr('FmaZ','AtSB', hyp) ! trace of FzT*(Zt*Z-Z*Zt)
      hp_zzz=hp_zzz+3.0d0*hyp
      call matprodtr('FmaX','AtSB', hyp) ! trace of FxT*(Zt*Z-Z*Zt)
      hp_zzx=hp_zzx+hyp
      call matprodtr('FmaY','AtSB', hyp) ! trace of FyT*(Zt*Z-Z*Zt)
      hp_zzy=hp_zzy+hyp
c
c----------------------------------------------------------------------
      call matrem('temp')
      call matrem('sma0')
      call matrem('BSAt')
      call matrem('AtSB')
c----------------------------------------------------------------------
      call matrem('PenZ')
      call matrem('PenY')
      call matrem('PenX')
c----------------------------------------------------------------------
c
      call matrem('DenZ')
      call matrem('DenY')
      call matrem('DenX')
c
      call matrem('FmaZ')
      call matrem('FmaY')
      call matrem('FmaX')
c----------------------------------------------------------------------
c     call getmem(0,last2)
c     write(6,*)' last1=',last1,' last2=',last2
c----------------------------------------------------------------------
      hypol(1)= hp_xxx
      hypol(2)= hp_yxx
      hypol(3)= hp_yyx
      hypol(4)= hp_yyy
      hypol(5)= hp_zxx
      hypol(6)= hp_zyx
      hypol(7)= hp_zyy
      hypol(8)= hp_zzx
      hypol(9)= hp_zzy
      hypol(10)= hp_zzz
c----------------------------------------------------------------------
c
c Experimentally determined hyper polarizability invariant is  
c (P. Calaminici, K.Jug, A.M.Koster, J.C.P., 109,(18),7756 (1998)
c
c mean hyperpolarizability HPmean= 3/5* (HP_xxz+HP_yyz+HP_zzz)
c
      HPmean= 0.6d0*(hp_zxx+hp_zyy+hp_zzz)
c----------------------------------------------------------------------
c     write(6,*)
c    *'----------------------------------------------------------------'
c     write(6,*)' '
c     write(6,*) '   ----------- hyperpolarizability -----------'
c     write(6,*)' '
c
c     write(6,61) -hp_xxx
c 61  format('| XXX=',f9.4,' |      ', 9x ,' |      ', 9x ,' |')
c     write(6,62) -hp_yxx,-hp_yyx,-hp_yyy
c 62  format('| YXX=',f9.4,' |  YYX=',f9.4,' |  YYY=',f9.4,' |')
c     write(6,63) -hp_zxx,-hp_zyx,-hp_zyy
c 63  format('| ZXX=',f9.4,' |  ZYX=',f9.4,' |  ZYY=',f9.4,' |')
c     write(6,64) -hp_zzx,-hp_zzy,-hp_zzz
c 64  format('| ZZX=',f9.4,' |  ZZY=',f9.4,' |  ZZZ=',f9.4,' |')
c     write(6,*)' '
c     write(6,*)
c    *'----------------------------------------------------------------'
c----------------------------------------------------------------------
c
      end
c======================================================================
      subroutine make_pro1(p2,ncf,p1)
      implicit real*8 (a-h,o-z)
      dimension p1(ncf,ncf), p2(ncf,ncf)
      parameter (one=1.0d0)
c
      do i=1,ncf
         do j=1,ncf
            p1(i,j)=-p2(i,j)
         enddo
      enddo
c
      do i=1,ncf
         p1(i,i)=one-p2(i,i)
      enddo
c
      end
c======================================================================
      subroutine symm4dft(natoms)
      use memory
      implicit real*8 (a-h,o-z)
      common /symm/nsym,nsy(7)         ! old texas symmetry data
c---------------------------------------------------------------
c Only abelian point group symmety is used 
c---------------------------------------------------------------
c gets number of symmetry unique atoms and its list as needed in DFT
c                                 NQ    and  IUQ
c---------------------------------------------------------------
      call getival('nsym',nsym)
c---------------------------------------------------------------
c  allocate memory for symmetry unique atom list
c
      call getmem(natoms,iuq)
c---------------------------------------------------------------
c  Get abelian point group symmetry data 
c
      call getsymdata_abel(natoms,nsym  ,nq,bl(iuq))
c
c output : nq and list in bl(iuq)
c---------------------------------------------------------------
c save number of symmetry unique atoms and its list (pointer iuq)
c
      call setival('nofuniq',nq)
      call setival('iunqato',iuq)
c---------------------------------------------------------------
c
      end
c======================================================================
      subroutine zerobelowt(ndim,thrsh,array)
      implicit real*8(a-h,o-z)
      dimension array(ndim)
      parameter (zero=0.0d0)
c
      do i=1,ndim
         if(abs(array(i)).lt.thrsh) array(i) = Zero
      enddo
c
      end
c=======================================================================
      subroutine getsymdata_abel(natoms,nsym  ,nq,iuq)
c-----------------------------------------------------------------------
c  only abelian subgroup symmetry data
c-----------------------------------------------------------------------
      use memory
      implicit real*8 (a-h,o-z)
      dimension iuq(*)
c-----------------------------------------------------------------------
c  prepare symmetry data (done even if no symmetry)
c  (be careful as DFT only needs real atoms)
c  only abelian subgroup
c
c iuq  - symmetry-unique atom list
c
      call getmem(3*natoms,ixnc)     ! nuclear coordinates
      call getmem(natoms,ian)        ! atomic charges/numbers
      call getmem(natoms,isc)        ! scratch
c
      call getnucdat(natoms,bl(ixnc),bl(ian))
c
      nupr = 1      ! not set if symmetry not used
      If(nsym.gt.0) call getival('SymNuPr1',nupr)
      call getunq(natoms,nsym,bl(ixnc),bl(nupr),nq,   iuq ,bl(isc)) 
      call retmem(3)
c
c needed : nq and list in iuq
c-----------------------------------------------------------------------
      end
c=======================================================================
      subroutine getsymdata_full(natoms,nq,iuq)
c-----------------------------------------------------------------------
c  Full point group symmetry data
c  output :
c nq - number of symmetry unique atoms
c iuq - list
c-----------------------------------------------------------------------
c
      use memory
      implicit real*8 (a-h,o-z)
      character jobname*256,cdum*20,GROUP*4
      Logical symflag
      dimension rm(3,3)
      dimension iuq(*)
      PARAMETER (IUnit=1)               ! unit number for checkpoint I/O
c-----------------------------------------------------------------------
c
C       Read from the <sym> file
C         number of atoms
C         number of symmetry operations
C
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c
        OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.sym',
     $      FORM='FORMATTED',STATUS='OLD')
cccc    call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
        call rdcntrl(IUnit,7,'$ntrans',1,NTrans,dum,cdum)
        CLOSE(UNIT=IUnit,STATUS='KEEP')
c-----------------------------------------------------------------------
c  Full point group symmetry data (not only abelian subgroup)
c  is availabe from the <sym> file
c
      call getmem(ntrans*9,itn)
      call getmem(ntrans*natoms,inq)
c
c read all full symmetry data from <sym> file
c
      symflag = ntrans.GT.1
      call rdsym(symflag,natoms, rm,   GROUP,  ntrans,
     $           ndeg,   nq,        iuq ,bl(itn),bl(inq))
c
      call retmem(2)
c-----------------------------------------------------------------------
      end
c=======================================================================
      subroutine printPol(iout,icon,polar,am,aa)
      implicit real*8 (a-h,o-z)
      dimension polar(9)
c
c--------------------------------------------------------
c -- print out
c
      write(iout,1200)
      write(iout,1300)
      write(iout,1200)
      write(iout,1400) polar(1),polar(2),polar(5),
     $                 polar(3),polar(6),polar(9)
      write(iout,1200)
      write(iout,1500) Am, Aa
      write(iout,1200)
ctest
c     write(iout,1400) polar(1),polar(4),polar(5),
c    $                 polar(7),polar(8),polar(9)
ctest
c log file 
      write(icon,1200)
      write(icon,1300)
      write(icon,1200)
      write(icon,1400) polar(1),polar(2),polar(5),
     $                 polar(3),polar(6),polar(9)
      write(icon,1200)
      write(icon,1500) Am, Aa
      write(icon,1200)
c--------------------------------------------------------
 1200 format(72('-'))
 1300 format(
     * '--                   DIPOLE POLARIZABILITY',
     * ' TENSOR                     --',/
     * '--                          (atomic',
     * ' units)                            --')
 1400 FORMAT(8X,'XX',10X,'XY',10X,'YY',10X,'XZ',10X,'YZ',10X,'ZZ',/,
     $         6F12.4)
ccc  $         6F12.6)
 1500 FORMAT(6X,'Mean polarizability =',f10.4,'   Anisotropy =',f10.4)
c1500 FORMAT(6X,'Mean polarizability =',f12.6,'   Anisotropy =',f12.6)
c--------------------------------------------------------
c
      end
c=======================================================================
      subroutine printHypol(iout,icon,hypol,hpm)
      implicit real*8 (a-h,o-z)
      dimension hypol(10)
c
c--------------------------------------------------------
c print into log file
      write(iout,1200)
      write(iout,2300)
      write(iout,1200)
      write(iout,2400) 
     * hypol(1),hypol(2),hypol(3),hypol(4),hypol(5)
      write(iout,*)' '
      write(iout,2500) 
     * hypol(6),hypol(7),hypol(8),hypol(9),hypol(10)
      write(iout,1200)
      write(iout,2600) hpm
      write(iout,1200)
c
c print into out file
      write(icon,1200)
      write(icon,2300)
      write(icon,1200)
      write(icon,2400) 
     * hypol(1),hypol(2),hypol(3),hypol(4),hypol(5)
      write(icon,*)' '
      write(icon,2500) 
     * hypol(6),hypol(7),hypol(8),hypol(9),hypol(10)
      write(icon,1200)
      write(icon,2600) hpm
      write(icon,1200)
c
 1200 format(72('-'))
 2300 format(
     * '--                 DIPOLE HYPERPOLARIZABILITY',
     * ' TENSOR                  --',/
     * '--                          (atomic',
     * ' units)                            --')
 2400 FORMAT(10X,'XXX',9X,'YXX',9X,'YYX',9X,'YYY',9x,'ZXX',/,2x,5F12.4)
 2500 FORMAT(10X,'ZYX',9X,'ZYY',9X,'ZZX',9X,'ZZY',9x,'ZZZ',/,2x,5F12.4)
c2600 FORMAT(6X,'Mean hyperpolarizability =',f12.4)
 2600 FORMAT(6X,'Mean hyperpolarizability =',f12.4,
     *   '  [ 3/5(h_xxz+h_yyz+h_zzz) ]' )
c-------------------------------------------
c
      end
c=======================================================================
      subroutine printDip(iout,icon,dipole,dipolt)
      implicit real*8 (a-h,o-z)
      dimension dipole(3)
c--------------------------------------------------------
c -- print out
c
      write(iout,1200)
      write(iout,1100)
      write(iout,1200)
      write(iout,1300)
      write(iout,1200)
      write(iout,1400) dipole(1),dipole(2),dipole(3),dipolt
      write(iout,1200)
c log file 
      write(icon,1200)
      write(icon,1100)
      write(icon,1200)
      write(icon,1300)
      write(icon,1200)
      write(icon,1400) dipole(1),dipole(2),dipole(3),dipolt
      write(icon,1200)
c--------------------------------------------------------
 1100 format(
     * '                  Molecule in Standard Orientation' )
 1200 format(72('-'))
 1300 format(
     * '--                          DIPOLE MOMENT ',
     * '                            --',/
     * '--                          (atomic',
     * ' units)                            --')
 1500 FORMAT(6X,'Total dipole =',f10.4)
 1400 FORMAT(
     * 10X,' X ',9X,'   ',4X,' Y ',9X,'   ',4x,' Z ',14x,'Total',/,
     *           2x,f12.4, 7x,f12.4, 7x,f12.4,7x,f12.4)
c--------------------------------------------------------
      end
c=======================================================================

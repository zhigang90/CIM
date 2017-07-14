C *****************************************************
C *  generate Cluster-in-Molecule cluster input file  *
C *  Modify the code from GAMESS                      *
C *  NZG_4/16/2016 @@UARK                             *
C *****************************************************

* -- LIST OF VARIABLES [integer] ---
* NATOM:    Number of atoms
* NBas:     Number of basis functions
* Nmo:      Number of all MOs (Number of independant functions, Nmo.le.Nbas)
* nocc:     Number of occupied MOs
* icha:     Number of charges
* mult:     Multiplicity
* nel:      Number of electrons
* ncs:      Number of contracted shells
* nsh:      Number of primitive shells
* tbs(nbas):Atomic label for basis
* nfocc:    Number of frozen occ MOs (default: nfocc=0)
* nprint:   Level of printing
*        nprint=0 ~ Print main important information of CIM
*        nprint=1 ~ Print detailed information of CIM
**     
*  nuchar(NATOM): Nuclear charges
*  NV:     The number of vir MOs
*  maxbs:  The number of basis functions of the largest clusters
*  nclu:   The number of clusters
*  numc:   The number of non-hydrogen atoms
*  numh:   The number of hydrogrn atoms
*     SR(NATOM,Nmo):     Sorted Mulliken population of basis in atoms
*     SR1(NATOM,Nmo):    Unsorted Mulliken population of basis in atoms
*  dis(NATOM,NATOM):  Distance matrix of atoms
*  link(NATOM,NATOM): Link matrix of atoms (0 or 1)
*  link2(NATOM,NATOM):Bond distance matrix of atoms (the minimum bonds between two atoms) 
*                        e.g. for C(1)-C(2)-C(3) link2(1,3)=2 (two bonds)
*     SNB(NUW):          The number of atoms for the occupied LMOs
*     SOB(NATOM,NUW):    The labels of atoms for the occupied LMOs
*     ISUB(NUW):         The number of MOs in the central domain
*     ISUB2(NUW):        The number of central MOs in the central domain
*     ISUB3(NUW):        The number of central MOs in the full domain
*     ISUB4(NUW):        The number of MOs in the full domain
*     NB(NATOM):         The first basis function of each atom
*     natc(NATOM):       The labels of all non-hydrogen atoms
*     nath(NATOM):       The labels of all hydrogen atoms
*     NAsub(NUW):        The number of atoms in the full atomic domain
*     NWsub(NUW):        The number of basis function in the full atomic domain
*     INF(NUW,NUW):      The labels of MOs in the central domain
*     INF2(NUW,NUW):     The labels of MOs in the full domain
*     BA(NATOM,NUW):     The labels of atoms in the full domain
*     ZA(NATOM,NUW):     The labels of basis functions in the full domain

* -- LIST OF VARIABLES [real(kind=8)] ---
*  NN:     Nuclear Repulsion Energy',NN,k)
*  ETOT:   Total MP2 correlation energy
*  nuchar(NATOM): Nuclear charges
*  coor(3,NATOM): Cartesian coordinates in Angstrom
*  XC(3,NATOM):   Cartesian coordinates in Bohr
*  eorb(Nmo):     Orbital energies
*  SMO(NBas,Nmo): Alpha MOs coefficients
*  SOVER(NW,NW):  AO overlap matrix
*  HCORE(NW,NW):  AO core-Hamiltonian matrix
*  FK(Nmo,Nmo):   AO Fock matrix
*  FIJ(NUW,NUW):  MO Fock matrix  ! 2008.07.25 Nmo --> NUW
*

      subroutine cimsubgen(nclu)

      use memory
      implicit real*8 (a-h,o-z)

      integer nclu
      PARAMETER (one=1.0D0,UNITS=ONE/0.52917724924D+00)
      COMMON /tape/ inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      integer NA,NB
      real(kind=8) ZAN,C,STAR,FINI
      integer inp,io,igms,itmp,isys,job
      integer i,j,k,l,m,n,i1,i2,i3,i4,j1,j2,j3,j4,k1,k2,k3,k4,ii,jj,kk
      real(kind=8),dimension(:,:),allocatable::SOVER,HCORE,DM,PS
C
CC    ******Added for GEBF idea to calculate energy*****
C      interface
C          subroutine gebf_deri_s_cpp(nprim, nfrag, SS1, coef)
C          !DEC$ ATTRIBUTES C, ALIAS:'gebf_deri_s_cpp'::gebf_deri_s_cpp
C          integer::nprim
C          !DEC$ ATTRIBUTES REFERENCE::nprim
C          integer::nfrag
C          !DEC$ ATTRIBUTES REFERENCE::nfrag
C          integer::SS1
C          !DEC$ ATTRIBUTES REFERENCE::SS1
C          real(kind=8)::coef
C          !DEC$ ATTRIBUTES REFERENCE::coef
C          end subroutine
C
C          subroutine gebf_deri_m_cpp(nprim, nfrag, SS1, coef)
C          !DEC$ ATTRIBUTES C, ALIAS:'gebf_deri_m_cpp'::gebf_deri_m_cpp
C          integer::nprim
C          !DEC$ ATTRIBUTES REFERENCE::nprim
C          integer::nfrag
C          !DEC$ ATTRIBUTES REFERENCE::nfrag
C          integer::SS1
C          !DEC$ ATTRIBUTES REFERENCE::SS1
C          real(kind=8)::coef
C          !DEC$ ATTRIBUTES REFERENCE::coef
C          end subroutine
C
C          subroutine gebf_deri_l_cpp(nprim, nfrag, SS1, coef)
C          !DEC$ ATTRIBUTES C, ALIAS:'gebf_deri_l_cpp'::gebf_deri_l_cpp
C          integer::nprim
C          !DEC$ ATTRIBUTES REFERENCE::nprim
C          integer::nfrag
C          !DEC$ ATTRIBUTES REFERENCE::nfrag
C          integer::SS1
C          !DEC$ ATTRIBUTES REFERENCE::SS1
C          real(kind=8)::coef
C          !DEC$ ATTRIBUTES REFERENCE::coef
C          end subroutine
C      end interface
C      integer,allocatable::SYS_ATOM(:,:),SYS_tmp(:,:),ISY_ATOM(:)
C      integer,allocatable::SYS_ATOMtmp(:,:)
C      integer::inarray1,nprim,maxatom,nprim_dis,ngrp_di
C      integer::COEFOUT,COEFOUT2
C      real(kind=8),allocatable::coef_GEBF(:)
C      integer::ngrp_dis,minc,single_clu,maxgrp,larger,min1,min2
C      integer,allocatable::natm_grp(:),atm_grp(:,:),atm2grp(:)
C      integer,allocatable::natm_grp_tmp(:),atm_grp_tmp(:,:)
C      real(kind=8)::mindis,thre_grp
C      integer,allocatable::ISY_GRPtmp(:),SYS_GRPtmp(:,:)
C      integer,allocatable::ISY_GRP(:),SYS_GRP(:,:)
C      integer,allocatable::min_dis(:,:,:)
C      real(kind=8),allocatable::dis_grp(:,:)  
C      integer,allocatable::norb_atm(:),orb_atm(:,:)
C      integer,allocatable::norb_grp(:),orb_grp(:,:),orb2grp(:)
C      integer,allocatable::SYS_ini(:,:),SYS_orb(:,:)
C      character*4::ENTYP   


      parameter (nopt=14)
      integer nsh,nprint,nbas,NATOM,ncore
      integer nmo,nalpha,nbeta,ierr,nfocc,np,jerr,debug,nocc
      integer ioptyp(nopt)
      character*4 word(nopt),motype
      character*8 submtd,highmtd
      character*7 methods(5)
      character*256 chopv(nopt),jobname,MOS,MOB,cluname,inpname
      character*256 infname,mosname,cenname,act_str,high_str
      character line*100
      dimension iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      integer,allocatable::IAN(:),natc(:),nath(:),link(:,:),link2(:,:)
      integer,allocatable::atmclu(:,:,:),natmclu(:,:),Z(:),INX(:,:)
      integer,allocatable::ILST(:,:),NBAtm(:),batom(:,:),tbs(:)
      integer,allocatable::SR(:,:),SR1(:,:),SOB(:,:)
      integer,allocatable::atmlevl(:),molevl(:),nocc_atm(:),occ_atm(:,:)
      integer,allocatable::orb_clu(:,:),ML(:,:),norb_clu(:),ncenorb(:)
      integer,allocatable::highorlow(:),clumtd(:),atm_clu(:,:),nb_clu(:)
      integer,allocatable::natm_clu(:),bas_clu(:,:),nvir_clu(:)
      integer,allocatable::ZA(:,:),BA(:,:),cenMO(:)
      integer,allocatable::ncen_clu(:),ncore_clu(:)
      real*8,allocatable::XC(:,:),coor(:,:),QA(:),dis(:,:),basdat(:,:)
      real*8,allocatable::SMO(:,:),SMOW(:,:),SMOT(:,:),eorb(:)
      real*8,allocatable::SMOAW(:,:),SMOAT(:,:),SMOA(:,:),norm(:,:)
      real*8,allocatable::SMOBW(:,:),SMOBT(:,:),SMOB(:,:),tmp(:)
      real*8,allocatable::MOCOOR(:,:),modis(:,:),PAO(:,:),MOS1(:,:)
      real*8,allocatable::MO_clu(:,:),MOS2(:,:),FK(:,:),FKW(:,:)
      real*8,allocatable::F2(:,:),FH(:,:),FHH(:),FIA(:,:),TRANS(:,:)
      real*8,allocatable::MAT(:,:),VECT(:,:),VALU(:),VC(:),dtmp(:)
      real*8,allocatable::virtime(:),DIAF(:),S2(:,:),MOS3(:,:)
      real*8,allocatable::density(:,:)
      real*8,external::dtrace2
      character*8,allocatable::AtSymb(:)

c Below are used for define the active space to calculate the
C correlation energy
C NZG_11/20/2016 @@NJU
C NZG_11/29/2016 @@NJU --Combine active space to one cluster
      logical,allocatable::act_atm(:),act_orb(:),high_atm(:)
      logical act,combine
      integer,allocatable::ncenatom(:),cenatom(:,:)
      integer,allocatable::orb_clu_act(:,:),norb_clu_act(:)
      integer,allocatable::ncenorb_act(:)
      integer,allocatable::ncenatom_act(:),cenatom_act(:,:)
      integer lineardep

C Add force option for gradient calculation
      logical calforce,calforce2

C Below are used for a new way of constructing virtual space: construct
C PAOs with a cluster
C NZG_4/7/2017 @@UARK
C Failed finally.... Zhigang 4/19/2017 @@UARK
      logical fpao,cpao

C Add none frozen core option for none frozen core calculation
      logical nofrozen

c-----------------------------------------------------------------------
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c-----------------------------------------------------------------------
      data ioptyp  / 1,     11,    21,    11,    21,    21,    0,
     &               21,    11,    21,    0,     0,     0,     0  /

C                    1      11     21     11     21     21     0
      data word    /'prin','disl','subm','virt','moty','acti','comb',
C                    21     11     21     0      0      0      0
     &              'hmtd','hdis','hatm','fpao','cpao','forc','nofr' /
      data IUnit   /1/ ! unit number for checkpoint I/O
      data methods /'RIMP2  ','MP2    ','CCD    ','CCSD   ','CCSD(T)'/
      common /job/ jobname,lenJ

C Options:
C prin: print level
C disl: distance threshold for occupied orbitals
C subm: correlation level for subsystems
C virt: threshold for projecting PAOs in full pao method
C moty: use QCMO or LMO in the subsystem calculation
C       LMO still has some problem in convergence
C acti: index of atoms in the active region
C comb: combine high level active regions
C hmtd: correlation level for high level subsystems
C hdis: distance threshold for occupied orbitals treated by high level 
C hatm: index of atoms in the high level regions
C fpao: use PAOs of the whole system
C cpao: construct PAOs within a cluster -NZG_4/7/2017 @@UARK
C forc: force will be calculated. In this step, the FORCE option will be
C       written to the input files of the subsystems

 10   format(72('='))
      write(iout,10)  
      write(iout,"(20x,'Generate CIM cluster input files')")
      write(iout,*) 
      dislmo=5.5D0
      dislmo2=4.0D0
      submtd='mp2'
      virt=0.05D0
      nprint=0
      act=.false.
      combine=.false.

C Add LMO version of CIM cluster calculation
C NZG_10/6/2016 @@NJU
      motype='qcmo'

      fpao=.true.
      cpao=.false.

      calforce2=.false.
      nofrozen=.false.

      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
      if (ifound(1).gt.0) nprint=iopv(1,1)
      if (nprint>0) nprint=iout
      if (ifound(2).gt.0) dislmo=ropv(1,2)
      if (ifound(3).gt.0) submtd=chopv(3)
      if (ifound(4).gt.0) virt=ropv(1,4)
      if (ifound(5).gt.0) motype=chopv(5)
      if (ifound(13).gt.0) calforce2=.true.
      if (ifound(14).gt.0) nofrozen=.true.
C !!!!!!!!!!!!!!!!!!!!!!!Caution!!!!!!!!!!!!!!!!!
C Consider the core orbitals later
      calforce=.false.
      call elapsec(tstart)
C
C  Read the Cartesian coordinates
C
      call getival('na',NAtom)
      allocate(act_atm(NAtom))
      act_atm=.false.
      if (ifound(6).gt.0) then
         act_str=chopv(6)
         write(6,*) "Active atoms:"//trim(act_str)
         call get_act_atm(act_str,act_atm,NAtom)
         act=.true.
         if (ifound(7).gt.0) combine=.true.
      endif

      if (ifound(12).gt.0) then
         cpao=.true.
         fpao=.false.
      endif

      allocate(XC(3,NAtom),AtSymb(NAtom),IAN(NAtom),QA(NAtom))
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcoord(IUnit,NAtom,AtSymb,XC,-1,jnk) ! In Bohr
      CLOSE(UNIT=IUnit,STATUS='KEEP')

C  get atomic numbers from atomic symbols
      CALL GetAtNo(NAtom,AtSymb,IAN)
C
C  get actual atomic charges (there may, e.g., be ghost atoms)
      CALL GetAtChrg(NAtom,QA)

C --- natc(numc): non-hydrogen atoms; nath(numh): hydrogen atoms
      allocate(natc(NATOM),nath(NATOM))
      numc=0;numh=0
      do i=1,NATOM
         if (IAN(i).ne.1) then
            numc=numc+1;natc(numc)=i
         else
            numh=numh+1;nath(numh)=i
         endif
      enddo

C --- the distance and linkage between atoms
      allocate(Coor(3,NATOM))
      coor=XC/UNITS    !In Angstrom
      allocate(dis(NATOM,NATOM),link(NATOM,NATOM),link2(NATOM,NATOM))
      call dislink(nprint,NATOM,IAN,Coor,dis,link)
      call NJ_idis(NATOM,link,link2)
      link=link2
      deallocate(link2)
      allocate(atmclu(NATOM,NATOM,7),natmclu(NATOM,7))
      call Atom_clu(NATOM,dis,atmclu,natmclu)

C  Read the basis set data
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.basis',
     $      FORM='FORMATTED',STATUS='OLD')
      call getival('ncs',ncs)
      call getival('nsh',nsh)
      call getival('ncf',nbas)
      allocate(Z(natom),ILST(4,ncs),BASDAT(13,nsh))
      call rdbasis(IUnit,NAtom,AtSymb,XC,Z,ILST,BASDAT)
      deallocate(Z)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
C ....................................................................
C  Be careful with basis order!
C  We need basis ordered PER ATOM for CIM calculation.
C  Wolinski ordering (for integral evaluation and hence for the MOs)
C  is PER SHELL (e.g. all S functions together regardless of atom)
C
C  The following routines are relevant here
C    SortBAS1  -  sorts basis per shell (but NOT full Wolinski)
C    reorder   -  orders basis from SortBAS1 into Wolinski
C    SortBAS2  -  orders basis per atom and supplies integer
C                 ordering array relating atom and Wolinski order
C ....................................................................
C
      allocate(INX(12,ncs),NBAtm(NATOM))
      CALL SortBAS1(NAtom,ncs,ILST,INX)
      CALL normaliz(ncs,INX,BASDAT)
C
C  get number of basis functions per atom
      CALL BasATOM(NAtom,Ncs,ILST,NBAtm)
      ini=1
      allocate(tbs(nbas),batom(nbas,natom))
      do i=1,natom
         ifi=ini+nbatm(i)-1
         tbs(ini:ifi)=i
         batom(1,i)=ini
         if (nbatm(i)>1) then
            do j=2,nbatm(i)
               batom(j,i)=batom(j-1,i)+1
            end do
         end if
         ini=ifi+1
      end do

      if (nprint>0) then
         write(nprint,*) '+++ Number of basis for atoms and labels+++'
         do i=1,NATOM
           write(nprint,'(2i4,2x,20i4)')
     &           i,NBAtm(i),(BATOM(j,i),j=1,NBatm(i))
         enddo
         write(nprint,*)
      endif


C      write(6,*) ' Number of basis functions per atom is:'
C      do i=1,natom
C         write(6,*) i,nbatm(i)
C      enddo

C  need to reorder MOs as need basis functions ordered per atom
C  and MOs have Wolinski special order
C
      i1=1
      i2=i1+NBas
      i3=i2+ncs
      IEnd=i3+12*ncs-1
c
      allocate(Z(IEnd))
      Call reorder(ncs,INX,Z(i1),Z(i2),Z(i3),IAN)     ! Wolinski order
      Call SortBAS2(ncs,Z(i3),Z(i2),Z(i1),INX)        ! per atom
C
C set up the MOs filenames MOS and MOB
C
      call tstchval('mosfname',iyes)
      if (iyes.eq.1) then
         call getchval('mosfname',MOS)
      else
         MOS=jobname(1:lenJ)
      end if
      call rmblan2(MOS,256,lenM)
      if (lenM.le.252) MOS(lenM+1:lenM+4)='.los'
      lenM=lenM+4
      MOB=MOS
      MOB(lenM:lenM)='b'
C
C  NZG_4/16/2016 Get the coefficients of LMOs from XXX.los file
C  Read the old MOs (from the existing binary file)
C  the number of nmo got from getival is not the actual number of MOs
C  but the number of alpha MOs
C  
      call getival('NAlpha',nalpha)  !alpha occ MOs
      call getival('NBeta',nbeta)    !beta  occ MOs
      
      if (nbeta>0) then
         write(iout,*) "We do not treat with UHF."
         stop
      end if     
      nocc=nalpha+nbeta
      call getival('nonredun',nmo)
      lineardep=0
      if (nmo<nbas) lineardep=1
      allocate(SMOAW(nbas,nbas),SMOAT(nmo,nbas),SMOA(nbas,nmo))
C      allocate(SMOBW(nbas,nbeta),SMOBT(nbeta,nbas),SMOB(nbas,nbeta))
      allocate(SMO(nbas,nmo),eorb(nmo))
      itype=1
      CALL ReadMOS(NBas,SMOAW,eorb,.False.,lenM,MOS(1:lenM),itype,ierr)
      If (ierr.NE.0) call nerror(1,'MOS File Does Not Exist',0,0)
C      If (NBeta.GT.0) Then
C         CALL ReadMOS(NBas,SMOBW,eorb(nalpha+1:nocc),
C     &               .False.,lenM,MOB(1:lenM),1,IErr)
C         If(IErr/=0) Call nerror(1,'MOB File Does Not Exist',0,0)
C      EndIf
C
C  Sort MOs to per atom order and form overlap and density matrices
C
      call getival('Multip',imult)
      If (NBeta.EQ.0.AND.IMult.EQ.1) Then
         CALL ReorderMO(NBas,nmo,Z(i1),SMOAW,SMOAT)
         SMO=transpose(SMOAT)
C         CALL FormDEN(NBas,NAlpha,PA,PS)
C         Call Expand(NBas,PS,PA)
C      Else
C         CALL ReorderMO(NBas,NAlpha,Z(i1),SMOAW,SMOAT)
C         CALL ReorderMO(NBas,NBeta,Z(i1),SMOBW,SMOBT)
C         SMOT(1:nalpha,:)=SMOAT(:,:)
C         SMOT(nalpha+1:,:)=SMOBT(:,:)
C         SMO=transpose(SMOT)
C         CALL FormDENU(NBas,NAlpha,NBeta,PA,PB,PS,S)
C         Call Expand(NBas,PS,PA)
C         Call Expand(NBas,S,PB)
      EndIf

C  get the overlap matrix
C
      allocate(tmp(127008),SOVER(nbas,nbas))
      CALL inton2(0,NAtom,SOVER,INX,INX,0,0,BASDAT,BASDAT,XC,
     $            IAN,ncs,ncs,NBas,NBas,tmp)
C
C  get the Fock matrix
C
      allocate(FKW(nbas,nbas),FK(nbas,nbas))
      np4=4
      call matdef('fock','s',nbas,nbas)
      call matread('fock',np4,'fock_rhf')

C Check density matrices from canonical and localized orbitals
C      call matdef('den0','s',nbas,nbas)
C      call matread('den0',np4,'den0_rhf')
C      call matprint('den0',6)

C      allocate(density(nbas,nbas))
C      call NJ_denmat(0,nbas,nmo,nocc,SMOAW,density)
C      write(6,*) density

      iafm=mataddr('fock')
      do i=1,nbas
         do j=1,i
            FKW(j,i)=bl(iafm)
            FKW(i,j)=bl(iafm)
            iafm=iafm+1
         end do
      end do
      CALL ReorderFock(nbas,nbas,Z(i1),FKW,FK)

      if (nofrozen) then
         nfocc=0
      else
         call FRZORB(nprint,NATOM,IAN,nfocc)
      endif

C If CIM-MP2 gradient is to be calculated, we need the information of
C core orbitals. In this case, set nfocc as zero. It doesn't mean the
C number of frozen core orbitals is zero.
C NZG_4/27/2017 @@UARK
C      if (calforce) then
C         nfocc=0
C      else
C         nfocc=ncoreocc
C      endif
      nvalocc=nocc-nfocc

      write(iout,"(' NUMBER OF BASIS FUNCTIONS:',i5)") NBAS
      write(iout,"(' NUMBER OF MOS            :',i5)") nmo
      write(iout,"(' NUMBER OF OCC MOS        :',i5)") NOCC
      write(iout,"(' NUMBER OF FROZEN OCC MOS :',i5)") nfocc
      write(iout,*)

C --- Calculate Mulliken Population and Assign the LMOs to a certain atom---
      allocate(SR(natom,nmo),SR1(natom,nmo),SOB(nmo,natom))
      allocate(atmlevl(natom),molevl(nocc),highorlow(nocc))
C
      call NJ_lower(submtd)
      if (trim(submtd)=='mp2') then
         atmlevl=2
      elseif (trim(submtd)=='ccd') then
         atmlevl=3
      elseif (trim(submtd)=='ccsd') then
         atmlevl=4
      elseif (trim(submtd)=='ccsd(t)') then
         atmlevl=5
      endif

C --- Determine the high level for multi-level calculation ---
C --- NZG_12/20/2016 @@NJU
      allocate(high_atm(natom))
      high_atm=.false.
      if (ifound(8).gt.0) then
         highmtd=chopv(8) 
         if (ifound(9).gt.0) dislmo2=ropv(1,9)
         if (ifound(10).gt.0) then
            high_str=chopv(10)
            write(6,*) "High level atom center:"//trim(high_str)
            call get_act_atm(high_str,high_atm,Natom)
            do i=1,natom
               if (high_atm(i)==.true.) then
                  if (trim(highmtd)=='mp2') then
                     atmlevl(i)=2
                  elseif (trim(highmtd)=='ccd') then
                     atmlevl(i)=3
                  elseif (trim(highmtd)=='ccsd') then
                     atmlevl(i)=4
                  elseif (trim(highmtd)=='ccsd(t)') then
                     atmlevl(i)=5
                  endif
               endif
            enddo
         else
            write(6,*) "Please assign the high-level atoms"
            stop
         endif
      endif

      write(iout,"(' MO DISTANCE THRESH:',f5.2)") dislmo
      write(iout,"(' MO DISTANCE THRESH (high-level):',f5.2)") dislmo2
      write(iout,*)
         
      call CalMPop(nprint,NATOM,nbas,nocc,nmo,SOVER,SMO,nbatm,batom,
     &             SR,SR1,SOB,atmlevl,molevl)

C --- Distance matrix between occ MOs ---
      allocate(MOCOOR(3,nocc),nocc_atm(natom),occ_atm(nocc,natom))
C  nocc_atm(i) : Number of occ MOs centered on atom i
C  occ_atm(j,i): Label of jth occ MO on atom i
C  SOB(j,i)    : Label of atoms that jth occ MO distribute (l to s)
      MOCOOR=0.0D0; nocc_atm=0; occ_atm=0
      do i=1+nfocc,nocc
         j=SOB(i,1)
         MOCOOR(:,i)=coor(:,j)
         nocc_atm(j)=nocc_atm(j)+1
         occ_atm(nocc_atm(j),j)=i
      end do
      
      write(iout,'(1x,a)')'------- DISTRIBUTION OF OCC MOS -------'
      write(iout,'(1x,a13,a9,a10)')'INDEX  ATOM','NUC_CHA','MO_INDEX'
      j=0
      do i=1,NATOM
         if(nocc_atm(i)==0) cycle
         j=j+1
         call NJ_prtlab(line,nocc_atm(i),occ_atm(:nocc_atm(i),i))
         call NJ_trim(line,k1,k2)
         write(iout,'(2x,I4,3x,I4,3x,i3,6x,a)') 
     &                                      j,i,IAN(i),line(k1:k2)
      end do
      write(iout,'(1x,a)')'---------------------------------------'
      write(iout,*)

      allocate(act_orb(NOCC))
      act_orb=.false.
      do i=1+nfocc,nocc
         if (act_atm(SOB(i,1))==.true.) act_orb(i)=.true.
      enddo
     
      highorlow=0
      allocate(modis(nvalocc,nvalocc))
      do i=1+nfocc,nocc
         do j=1+nfocc,nocc
            modis(j-nfocc,i-nfocc)=(MOCOOR(1,i)-MOCOOR(1,j))**2+
     &      (MOCOOR(2,i)-MOCOOR(2,j))**2+(MOCOOR(3,i)-MOCOOR(3,j))**2
            modis(j-nfocc,i-nfocc)=sqrt(modis(j-nfocc,i-nfocc))
         end do
      end do
      allocate(orb_clu(nvalocc,nvalocc),norb_clu(nvalocc))
      allocate(ncenorb(nvalocc))
      call Fragment_Deter(nprint,nocc,nfocc,modis,orb_clu,ncenorb,
     &                    norb_clu,molevl,dislmo,dislmo2,highorlow)
      call Reduce_Frag(nprint,nocc,nfocc,orb_clu,ncenorb,norb_clu,
     &                 nclu,molevl,act_orb)
      write(6,"(' Final number of clusters:',I4)") nclu
      write(6,*)
      allocate(ncenatom(nclu),cenatom(NATOM,nclu))
      ncenatom=0; cenatom=0
      do i=1,nclu
         do j=1,ncenorb(i)
            k=SOB(orb_clu(j,i),1)
            call inarray(k,cenatom(:,i),NATOM,ierr)
            if (ierr==0) then
               ncenatom(i)=ncenatom(i)+1
               cenatom(ncenatom(i),i)=k
            endif
         enddo
      enddo

      write(iout,'(1x,a)')'------ DISTRIBUTION OF CENTRAL ATOMS ------'
      write(iout,*) ' INDEX   CenAtom'
      flush(iout)
      do i=1,nclu
         call NJ_prtlab(line,ncenatom(i),cenatom(:ncenatom(i),i))
         call NJ_trim(line,k1,k2)
         write(iout,'(i6,3x,a)') i,line(k1:k2)
      enddo
      write(iout,'(1x,a)')'-------------------------------------------'
      write(iout,*)

C Now the information of the initial cluster:
C nclu:     number of clusters
C orb_clu:  orbitals in the cluster
C norb_clu: number of orbitals in the cluster
C ncenorb:  number of central orbitals in the cluster
C ncenatom: number of central atoms in the cluster
C cenatom:  labels of central atoms in the cluster

C Remain only the clusters with active central atom
C 11/20/2016_NZG @@NJU
      if (act==.true.) then
         nclu_act=0
         allocate(orb_clu_act(nvalocc,nclu),norb_clu_act(nclu))
         allocate(ncenorb_act(nclu))
         allocate(ncenatom_act(nclu),cenatom_act(NATOM,nclu))
         loop1: do i=1,nclu
            loop2: do j=1,ncenatom(i)
               if (act_atm(cenatom(j,i))==.true.) then
                  nclu_act=nclu_act+1
                  norb_clu_act(nclu_act)=norb_clu(i)
                  orb_clu_act(:,nclu_act)=orb_clu(:,i)
                  ncenorb_act(nclu_act)=ncenorb(i)
                  ncenorb_act(nclu_act)=ncenorb(i)
                  ncenatom_act(nclu_act)=ncenatom(i)
                  cenatom_act(:,nclu_act)=cenatom(:,i)
                  exit loop2
               endif
            enddo loop2
         enddo loop1
         nclu=nclu_act
         norb_clu(1:nclu)=norb_clu_act(1:nclu)
         orb_clu(:,1:nclu)=orb_clu_act(:,1:nclu)
         ncenorb(1:nclu)=ncenorb_act(1:nclu)
         ncenatom(1:nclu)=ncenatom_act(1:nclu)
         cenatom(:,1:nclu)=cenatom_act(:,1:nclu)
         if (combine==.true.) then
            ncenatom(1)=sum(ncenatom(1:nclu))
            ncenorb(1)=sum(ncenorb(1:nclu))
            j=1; k=1; orb_clu=0; norb_clu=0
            do i=1,nclu
               cenatom(j:j+ncenatom_act(i)-1,1)=
     &         cenatom_act(1:ncenatom_act(i),i)
               j=j+ncenatom_act(i)
               orb_clu(k:k+ncenorb_act(i)-1,1)=
     &         orb_clu_act(1:ncenorb_act(i),i)
               k=k+ncenorb_act(i)
            enddo
            if (ncenorb(1)/=k-1) write(6,*) "Wrong combination!"
            norb_clu(1)=ncenorb(1)
            do i=1,nclu
               do ii=ncenorb_act(i)+1,norb_clu_act(i)
                  j=orb_clu_act(ii,i)
                  call inarray(j,orb_clu(:,1),nvalocc,ierr)
                  if (ierr==0) then
                     norb_clu(1)=norb_clu(1)+1
                     orb_clu(norb_clu(1),1)=j
                  endif
               enddo
            enddo
            nclu=1
         endif
      endif

      allocate(clumtd(nclu))
      do i=1,nclu
         clumtd(i)=molevl(orb_clu(1,i))
      end do     

C Construct ML array! Its name is from GAMESS
C ML(i,j): Label of ith occ MO in jth cluster
C It contains all the occ MOs with norb_clu(j) orbitals at the beginning
      allocate(ML(nocc,nclu))
      ML=0
      do i=1,nclu
         ML(1:norb_clu(i),i)=orb_clu(1:norb_clu(i),i)
      end do
      do i=1,nclu
         k=norb_clu(i)+1
         do j=1,nocc
            call inarray(j,ML(:,i),nocc,ierr)
            if (ierr==0) then
               ML(k,i)=j
               k=k+1
            end if
         end do
      end do

C  --- EMPLOY DISTANCE TO DETERMINE AO DOMAIN ----
C  atm_clu(i,j): Label of ith atom in jth cluster
C  natm_clu(j) : Number of atoms in jth cluster

      allocate(atm_clu(natom,nclu),natm_clu(nclu))
      atm_clu=0; natm_clu=0
      do i=1,nclu
         do l=1,ncenorb(i)
            j=SOB(orb_clu(l,i),1)
            do k=1,NATOM
               if (highorlow(orb_clu(l,i))==1) then
                  if (dis(k,j)<=dislmo2) then
                     call inarray(k,atm_clu(:,i),natom,ierr)
                     if (ierr==1) cycle
                     natm_clu(i)=natm_clu(i)+1
                     atm_clu(natm_clu(i),i)=k
                  end if
               else
                  if (dis(k,j)<=dislmo) then
                     call inarray(k,atm_clu(:,i),natom,ierr)
                     if (ierr==1) cycle
                     natm_clu(i)=natm_clu(i)+1
                     atm_clu(natm_clu(i),i)=k
                  end if
               end if
           end do
         end do
      end do
C  --- EMPLOY DISTANCE TO DETERMINE AO DOMAIN OVER ----

      allocate(bas_clu(nbas,nclu),nb_clu(nclu))
      bas_clu=0; nb_clu=0
      do i=1,nclu
         do j=1,natm_clu(i)
            k=atm_clu(j,i)
            ini=nb_clu(i)+1
            nb_clu(i)=nb_clu(i)+nbatm(k)
            bas_clu(ini:nb_clu(i),i)=batom(:nbatm(k),k)
         end do
      end do

C --- Construct ZA matrix. Its name is from GAMESS ---
C ZA(i): Reorder of basis function labels
C basis functions of the cluster at the beginning
      allocate(ZA(nbas,nclu))
      ZA=0
      do i=1,nclu
         ZA(1:nb_clu(i),i)=bas_clu(1:nb_clu(i),i)
      end do
      do i=1,nclu
         k=nb_clu(i)+1
         do j=1,nbas
            call inarray(j,ZA(:,i),nbas,ierr)
            if (ierr==0) then
               ZA(k,i)=j
               k=k+1
            end if
         end do
      end do
      
C --- Construct BA matrix. Its name is from GAMESS ---
C ZA(i): Reorder of atom labels
C atoms of the cluster at the beginning
      allocate(BA(natom,nclu))
      BA=0
      do i=1,nclu
         BA(1:natm_clu(i),i)=atm_clu(1:natm_clu(i),i)
      end do
      do i=1,nclu
         k=natm_clu(i)+1
         do j=1,natom
            call inarray(j,BA(:,i),natom,ierr)
            if (ierr==0) then
               BA(k,i)=j
               k=k+1
            end if
         end do
      end do
      
C --- Print out the information of initial CIM cluster ---
      write(iout,'(11x,26(''-''))')
      write(iout,'(11x,"INITIAL CIM DOMAIN SUMMARY")')
      write(iout,'(11x,26(''-''))')
      write(iout,*) ' Clu_INDX   NOcc   Natom   Nbas   Method   CenAtom'
      flush(iout)
      do i=1,nclu
         call NJ_prtlab(line,ncenatom(i),cenatom(:ncenatom(i),i))
         call NJ_trim(line,k1,k2)
         write(iout,'(i7,1x,3i8,5x,a8,a)') i,norb_clu(i),natm_clu(i),
     &                                     nb_clu(i),methods(clumtd(i)),
     &                                     line(k1:k2)
      enddo
      write(iout,*)

C If using fpao type virtual orbitals, construct the PAOs for the 
C whole system -NZG_4/7/2017 @@UARK
      if (fpao) then
         allocate(PAO(nbas,nbas))
         call dunit(nbas,PAO)
         call pjotorb(nbas,nocc,nbas,SMO(:,1:nocc),PAO,SOVER) 
         call normorb(nbas,nbas,PAO,SOVER)
      endif
   
      call elapsec(tend1)

C --- Construct virtual space for each cluster!
C nvir_clu(i)    : Number of virtual orbitals of cluster i
C MOS1(nbas,nbas): MO coefficients of the cluster. 
C                  Shape: nb_clu(i)*(norb_clu(i)+nvir_clu(i))
C virt           : Threshold for projecting PAOs (default: 0.05)
      allocate(nvir_clu(nclu),virtime(nclu))
      allocate(ncen_clu(nclu),ncore_clu(nclu))
      do i=1,nclu
         if (nprint>0) then
            write(nprint,'(" Construct vir space for clu",5i)') i
         endif
         call elapsec(virt1)
         allocate(MOS1(nbas,nbas),S2(nbas,nbas))
         if (fpao) then
            call COV_FPAO(nprint,i,natom,nmo,batom(1,:),natm_clu(i),
     &                    nb_clu(i),nbas,norb_clu(i),nocc,nvir_clu(i),
     &                    ncen_clu(i),ncenorb(i),ncore_clu(i),ncoreocc,
     &                    ML(:,i),BA(:,i),ZA(:,i),SOVER,S2,SMO,MOS1,dis,
     &                    link,nbatm,batom,atmclu,natmclu,PAO,virt,
     &                    calforce)
         else
            call COV_CPAO(nprint,i,natom,nmo,batom(1,:),natm_clu(i),
     &                    nb_clu(i),nbas,norb_clu(i),ncenorb(i),nocc,
     &                    ncoreocc,nvir_clu(i),ML(:,i),BA(:,i),ZA(:,i),
     &                    SOVER,S2,SMO,MOS1,dis,link,nbatm,batom)
         endif
C These names are from GAMESS
C KB:   Number of occ MOs in the cluster
C KV:   Number of vir MOs in the cluster
C NA:   Number of MOs in the cluster
C JF:   Number of basis functions in the cluster
C MOS2: QCMO coefficients in the cluster
C
         KB=norb_clu(i)
         NCORE=ncore_clu(i)
         NCEN=ncen_clu(i)
         NVAL=KB-NCORE
         KV=nvir_clu(i)
         NA=KB+KV
         JF=nb_clu(i)
         allocate(MO_clu(JF,NA),FHH(NA),F2(JF,JF))
         MO_clu(:,:)=MOS1(1:JF,1:NA)
         FHH=0.0D0; F2=0.0D0
         do k=1,JF
            do j=1,JF
               F2(j,k)=FK(ZA(j,i),ZA(k,i))
            enddo
         enddo
         deallocate(MOS1)

C---------------------------------------------------------------------
C --- Construct QCMOs for each cluster ---
C     Construct molecular Fock matrix of cluster
C---------------------------------------------------------------------

         allocate(FH(NA,NA),TRANS(NVAL,NVAL))
         FH=0.0d0
         
         TRANS=0.0D0
         do j=1,NVAL
            TRANS(j,j)=1.0D0
         enddo
C
C-WL- 2007.11.20 ADD FOR FOCK MATRIX DIAGONALIZATION --- liwei -----

         call MKL_Tfock(nprint,JF,NA,F2,MO_clu,FH)  !-WL- 2007.10.23
C        
         allocate(FIA(KB,KV))
         FIA(:,:)=0.0D0
         do k=1,KB
            do j=1,KV
               FIA(k,1)=FIA(k,1)+FH(k,KB+j)*FH(k,KB+j)
            end do
            FIA(k,1)=FIA(k,1)/real(KV*KV)
         end do
   
         FIA=FH(1:KB,KB+1:KB+KV)

         PP=0.0D0
         do k=1,KV
            do j=1,KB
               PP=PP+FIA(j,k)*FIA(j,k)
            end do
         end do
         PP=PP/real(KB*KV*KB*KV)
         deallocate(FIA)
         if (nprint>0) write(nprint,*)'THE RMSD OF FIA:',PP
C
C Diagonalize the core and valence occupied orbitals separately.
C NZG_5/17/2017 @@UARK
C
C First core occupied MOs
         allocate(MOS2(JF,NA))
         if (NCORE>0) then
            allocate(MAT(NCORE,NCORE),VECT(NCORE,NCORE))
            allocate(VALU(NCORE),VC(NCORE))
            MAT(1:NCORE,1:NCORE)=FH(NVAL+1:KB,NVAL+1:KB)
            call NJ_qr(nprint,MAT,NCORE,VECT,VALU,VC,k,1)
            FHH(1:NCORE)=VALU(1:NCORE)
C
            do j=1,NCORE
               do k=1,JF
                  MOS2(k,j)=0.0d0
                  do L=1,NCORE
                     MOS2(k,j)=MOS2(k,j)+MO_clu(k,L+NVAL)*VECT(L,j)
                  enddo
               enddo
            enddo
            deallocate(MAT,VECT,VALU,VC)
         endif

C Then valence occupied MOs
         allocate(MAT(NVAL,NVAL),VECT(NVAL,NVAL),VALU(NVAL),VC(NVAL))
         MAT(:,:)=FH(:NVAL,:NVAL)
         call NJ_qr(nprint,MAT,NVAL,VECT,VALU,VC,k,1)
         FHH(NCORE+1:KB)=VALU(1:NVAL)
C
         do j=1,NVAL
            do k=1,JF
               MOS2(k,j+NCORE)=0.0d0
               do L=1,NVAL
                  MOS2(k,j+NCORE)=MOS2(k,j+NCORE)+MO_clu(k,L)*VECT(L,j)
               enddo
            enddo
         enddo
 
         if (NCORE>0) then
            if (FHH(NCORE+1)<FHH(NCORE)) then
                write(iout,'(a,2f10.6)')
     &                 'Hcore LARGER THAN Lval:',FHH(NCORE),FHH(NCORE+1)
                write(iout,"('Message from cluster: ',I4)") i
               stop
            endif
         endif
         TRANS=transpose(VECT)
         deallocate(MAT,VECT,VALU,VC)
C

C Then virtual MOs

         allocate(MAT(KV,KV),VECT(KV,KV),VALU(KV),VC(KV))
         do k=1,KV
            do j=1,KV
               MAT(j,k)=FH(j+KB,k+KB)
            enddo
         enddo
         call NJ_qr(nprint,MAT,KV,VECT,VALU,VC,k,1)
         do j=1,KV
            FHH(KB+j)=VALU(j)
         enddo
C
         if (FHH(KB+1)<FHH(KB)) then
            write(iout,'(a,2f10.6)')
     &                 'HOMO LARGER THAN LUMO:',FHH(KB),FHH(KB+1)
            write(iout,"('Message from cluster: ',I4)") i
            stop
         end if
C
         do j=1,KV
            do k=1,JF
               MOS2(k,j+KB)=0.0d0
               do L=1,KV
                  MOS2(k,j+KB)=MOS2(k,j+KB)+MO_clu(k,L+KB)*VECT(L,j)
               enddo
            enddo
         enddo
         deallocate(MAT,VECT,VALU,VC)

C         kk=KB; KB=NA
C         allocate(MAT(KB,KB),VECT(KB,KB),VALU(KB),VC(KB))
C         MAT(:,:)=FH(:KB,:KB)
C         call NJ_qr(nprint,MAT,KB,VECT,VALU,VC,k,1)
C         FHH(1:KB)=VALU(1:KB)
CC
C         allocate(MOS2(JF,NA))
C         do j=1,KB
C            do k=1,JF
C               MOS2(k,j)=0.0d0
C               do L=1,KB
C                  MOS2(k,j)=MOS2(k,j)+MO_clu(k,L)*VECT(L,j)
C               enddo
C            enddo
C         enddo
C         allocate(TRANS(KK,KK))
C         TRANS=0.0D0
C         do j=1,KK
C            TRANS(j,j)=1.0D0
C         end do
C         TRANS=transpose(VECT(1:KK,1:KK))
C
C         KB=KK
C --Normalize the QCMOs
C --NZG_4/11/2017 @@UARK
         call normorb(JF,NA,MOS2,S2(1:JF,1:JF))
         deallocate(S2)
 
         if (nprint>0) write(nprint,"(1x,'Trace(Fock LMO)=',f20.10)")
     &                       dtrace2(NA,FH)
C
C --- Properties of the clusters 
C      nelec=nelec-mult+1 
106      nelec=KB*2  !Since we only deal with ristricted system

         nucl=0
         do k=1,natm_clu(i)
            j=BA(k,i)
            nucl=nucl+IAN(j)
         enddo
         mult=1
C
C   -- Modified By Guoyang --
C      icharg=nucl-nelec
         icharg=nucl
C   -- Modified By Guoyang Over --
         icharg=mod(icharg,2)  !-WL- 27 Aug 2009 4:03PM
C
C ===========================================
C --- Create PQS input file for each cluster
C ===========================================
         write(line,*) i
         call NJ_trim(line,k1,k2)
         cluname=jobname(1:lenJ)//'_Sys-'//line(k1:k2)
         call NJ_trim(cluname,k1,k2)
         inpname=trim(cluname)//'.inp'
         infname=trim(cluname)//'.cim'
         mosname=trim(cluname)//'.mos'
         if (calforce2) cenname=trim(cluname)//'.cen'
         lenM=k2-k1+5
         call NJ_trim(inpname,k1,k2)
         inpclu=123+i
         open(inpclu,file=inpname)

CC --- Estimate the memory needed for the clusters
C      nwords = memgms(scftyp,mplevl,cctyp,JF_sub,KB,NA-KBB,KU)
C      if (cctyp.eq.'CCSD(T) '.or.cctyp.eq.'CR-CCL  ') then
C         k = memgms(scftyp,mplevl,'CCSD    ',JF_sub,KB,NA-KBB,KU)
C         nwords = max(nwords,k)
C      endif
C      mwords = int(nwords/1000d0/1000d0) + 10
C      mbytes = int(nwords/125.d0/1000d0) + 10
C      write(io,190)scftyp,mplevl,cctyp,JF_sub,KB,NA-KBB,KU,nwords,mwords
C 190  format(1x,a8,i4,2x,a8,4i6,i12,' words',i6,' mwords')

C
         mwords=200 !set to 100 temporarily
         call SUBINP(inpclu,NATOM,AtSymb,coor,natm_clu(i),BA(:,i),
     &               clumtd(i),mwords,icharg,mult,lineardep,calforce2)  
         close(inpclu)

C- Print out some important information in XXX.cim file
         open(inpclu,file=infname)
         write(inpclu,*) '$INFO'
         write(inpclu,'(a8,a)') 'SUBMTD  ',methods(clumtd(i))

         if (motype=='qcmo') then
            write(inpclu,'(a8,a8)') 'MOTYP   ','QCMO    '
         else
            write(inpclu,'(a8,a8)') 'MOTYP   ','LMO     '
         endif

         write(inpclu,'(a8,i8)') 'SYS     ', I
         write(inpclu,'(a8,i8)') 'NAT     ', natm_clu(i)
         write(inpclu,'(a8,i8)') 'NELEC   ', nelec
         write(inpclu,'(a8,i8)') 'ICH     ', icharg
         write(inpclu,'(a8,i8)') 'MUL     ', mult
         write(inpclu,'(a8,i8)') 'NBAS    ', JF
         write(inpclu,'(a8,i8)') 'NOCC    ', KB
         write(inpclu,'(a8,i8)') 'NVIR    ', KV
         write(inpclu,'(a8,i8)') 'NMO     ', NA
         write(inpclu,'(a8,i8)') 'NCEN    ', ncen_clu(i)
         write(inpclu,'(a8,i8)') 'NCORE   ', ncore_clu(i)
         write(inpclu,*) '$END' 
C
         call iwrit8(inpclu,'$MO-OCC',KB,ML(1,i))

         if (calforce2) call iwrit8(inpclu,'$ATOMS',natm_clu(i),
     &                              BA(1:natm_clu(i),i))
C         
C         allocate(CenMO(KB))
C         CenMO(1:ncenorb(i))=1       !Central MOs labeled as 1
C         CenMO(ncenorb(i)+1:KB)=0    !Else as 0
         
C         call iwrit8(inpclu,'$MO-CEN',KB,CenMO(1)) !(CenMO(k,I)),k=1,KB)
         call rwrit8(inpclu,'$TRMX-A',KB*KB,TRANS(1,1))
         deallocate(TRANS)

         call rwrit8(inpclu,'$VEC',JF*NA,MOS2(1,1))
         if (motype=='lmo ') then
            allocate(DIAF(NA))
            do j=1,NA
               DIAF(j)=FH(j,j)
            enddo
            call rwrit8(inpclu,'$EIGVAL',NA,DIAF(1))
            deallocate(DIAF)
         else
            call rwrit8(inpclu,'$EIGVAL',NA,FHH(1))
         endif
C
         LF2=JF*(JF+1)/2
         allocate(dtmp(LF2))
         L=0
         do j=1,JF
            do k=1,j
               L=L+1
               dtmp(L)=F2(k,j)
            enddo
         enddo
         call rwrit8(inpclu,'$AO-FOCK-A',LF2,dtmp(1)) !((F2(k,j),k=1,j),j=1,JF)
         deallocate(dtmp,F2,FH)
         close(inpclu)
C
C -- Write MO coefficients and orbital energies to XXX.mos file
         if (motype=='lmo ') MOS2=MO_clu
         itype=1
         call WriteMOS(JF,NA,MOS2,FHH,.true.,lenM,mosname,itype)

C If force is to be calculated, the coefficients of central orbitals are
C needed. -NZG_5/22/2017 @@UARK
C for debug
         if (calforce2) call WriteMOS(JF,ncen,MO_clu(1:JF,1:ncen),FHH,
     &                                .false.,lenM,cenname,itype)
C         if (calforce2) call WriteMOS(JF,KB,MO_clu(1:JF,1:KB),FHH,
C     &                                .false.,lenM,cenname,itype)

         deallocate(MOS2,FHH,MO_clu)
         call elapsec(tvir2)
         virtime(i)=(tend1-tstart+tvir2-tvir1)/60.0D0
      end do 

      write(iout,'(11x,26(''-''))')
      write(iout,'(11x,"FINAL CIM DOMAIN SUMMARY")')
      write(iout,'(11x,26(''-''))')
      write(iout,*) 
     & ' Clu_INDX  Natom   Nbas     Nocc    Nvir   Method  Elap_time'
      flush(iout)
      do i=1,nclu
         write(iout,'(i7,1x,4i8,5x,a8,f8.2)') i,natm_clu(i),nb_clu(i),
     &         norb_clu(i),nvir_clu(i),methods(clumtd(i)),virtime(i)
      enddo
      write(iout,*)

      maxbas=0; maxsys=0
      do i=1,nclu
         if (nb_clu(i)>maxbas) then
            maxbas=nb_clu(i)
            maxsys=i
         endif
      enddo
      write(iout,'(" Max(nbs)=",I6," in sys",I4)') maxbas,maxsys

      end subroutine cimsubgen


C *****************************************************************
C * --Constructing virtual space of the clusters (PAO basis)      *
C *   Updated by NZG_2/26/2016 @@UARK                              *
C * --Change the name from COV to COV_FPAO -NZG_4/7/2017 @@UARK    *
C * --Deal with core orbitals and do Schmidt orthogonalization if *
C *   force is to be calculated. -NZG_4/27/2017 @@UARK             *
C *****************************************************************
      subroutine COV_FPAO(io,I,natom,nmo,nb,J0,JF,JF1,KB,KB1,nvext,
     &                    ncen,ncen1,ncore,ncore1,ML,BA,ZA,S,S2,SMO,
     &                    MOS1,dis,link,nbs_atm,bas_atm,atmclu,natmclu,
     &                    PAO,virthre,calforce)
      implicit none
      integer i,j,k,jzhi,kzhi,l,m,maxdis,nstar,nonzero,ncore,ncore1
      integer io,iout,NATOM,Nmo,J0,JF0,JF,KB,KB1,J1,JF1,nat_init,na
      integer nzero,ncut,nvext,ll,np,ierr,paoout,ncen,ncen1
      integer kb_aao,nat_fao,JF_fao,JF_red
      integer ML(KB1),BA(NATOM),ZA(JF1),nb(NATOM),link(NATOM,NATOM)
      integer atmclu(NATOM,NATOM,7),natmclu(NATOM,7)
      integer nbs_atm(NATOM),bas_atm(JF1,NATOM)
      integer,dimension(:),allocatable::BAA,ZA_AO,num_at_aao,J1_aao
      integer,dimension(:),allocatable::JF_aao,BA_FAO,paored,paored_full
      integer,dimension(:,:),allocatable::BA_AAO
      real(kind=8),parameter::zeps=1.0d-5
      real(kind=8),parameter::occthre=0.01D0
      real(kind=8) thre,P1,P2,eps,epst,cvszhi,virthre
      real(kind=8) S(JF1,JF1),SMO(JF1,nmo),MOS1(JF1,JF1),PAO(JF1,JF1)
      real(kind=8) dis(NATOM,NATOM),S2(JF1,JF1)
      real(kind=8),dimension(:),allocatable::MOSS,wf1,wf2,e,d
      real(kind=8),dimension(:,:),allocatable::co1,S4,tr,cvs,SS,Q,c
      real(kind=8),dimension(:,:),allocatable::smo1,PAO_tmp,MOSSS,Sij
      real(kind=8),dimension(:,:),allocatable::c2,coaao,cvs2,co1_aao
      logical calforce

C KB:            Number of occ MOs in the cluster
C KB1:           Number of occ MOs of the whole molecule
C J0:            Number of atoms in the cluster
C JF:            Number of basis functions in the cluster
C JF1:           Number of basis functions of the whole molecule
C S2:            Overlap of basis in the cluster at the beginning
C --Change it to an input and output to get the overlap within the
C --cluster to do normalization for the QCMOs
C --Zhigang 4/11/2017 @@UARK
C S4:            Overlap of basis in the cluster
C num_at_aao(:): Number of atoms in each AAO domain
C BA_AAO(:,:):   Label of each atom in each AAO domain

      iout=6

C In the CIM-MP2 gradient calculation, we need the core orbitals.
C Get the number of core orbitals in the cluster and rearrange the
C orbitals. The order of the orbitals:
C central valence -> buffer valence -> core -> PAOs
C NZG_4/27/2017 @@UARK

C ncen:   Number of central valence orbitals in the cluster 
C ncen1:  Number of cent orbitals of the cluster including core orbitals
C ncore:  Number of core orbitals in the cluster
C ncore1: Number of core orbitals of the whole molecule
      allocate(co1(JF1,KB))
      if (calforce) then
         ncore=0
         do k=1,KB
            if (ML(k)<=ncore1) ncore=ncore+1
         enddo
         j=0; l=0; ncen=0
         do k=1,KB
            if (ML(k)<=ncore1) then
               co1(:,KB-j)=SMO(:,ML(k))
               j=j+1
            else
               if (k<=ncen1) then
                  ncen=ncen+1
                  co1(:,ncen)=SMO(:,ML(k))
               else
                  co1(:,KB-ncore-l)=SMO(:,ML(k))
                  l=l+1
               endif
            endif
         enddo
      else
         do k=1,KB
            co1(:,k)=SMO(:,ML(k))
         enddo
         ncen=ncen1
         ncore=0
      endif

      allocate(S4(JF,JF))
      MOS1(:,1:KB)=co1(:,1:KB)
      call INDEX_REORDER(io,JF1,JF,KB,S,S2,S4,MOS1(:,:KB),ZA)
      allocate(cvs(JF,KB),wf1(KB),wf2(KB))
      call projorb(io,JF1,JF,KB,KB,MOS1(:,1:KB),cvs,S2,wf1,wf2)
      call locindx_occ(io,JF1,MOS1(:,1:KB),JF,cvs,
     &                                    KB,S4,wf1,occthre,nstar)
      deallocate(wf1,wf2)
      nat_init=J0;thre=0.0D0;JF0=JF
      allocate(BAA(JF1),ZA_AO(JF1))
      BAA=BA
      do while(1)
         BA=BAA
         J0=nat_init    !Zhigang 2016.3.1 
         deallocate(S4)
532      call Add_Atom(io,KB,NATOM,link,dis,BA,J0,J1,thre)
         if (J1==J0 .and. thre>0.5D0) then
            thre=thre+1.0D0
            goto 532
         end if
C         if (J1==NATOM) call nerror(1,'CIM Module',
C     &                  "AO domain contains all the atoms!",0,0)
         call BASIS_REORDER(io,JF1,NATOM,ZA,nb,BA,J1,JF)      
         allocate(S4(JF,JF))
         MOS1(:,:KB)=co1(:,:KB)
         call INDEX_REORDER(io,JF1,JF,KB,S,S2,S4,MOS1(:,1:KB),ZA) 
         deallocate(cvs)
         allocate(cvs(JF,KB),wf1(KB),wf2(KB))
         call projorb(io,JF1,JF,KB,KB,MOS1(:,1:KB),cvs,S2,wf1,wf2)
         call locindx_occ(io,JF1,MOS1(:,1:KB),JF,cvs,KB,S4,wf1,
     &                    occthre,nstar)
         deallocate(wf1,wf2,S4,cvs)
    
         ZA_AO=ZA

         if (io>0) then
1000        write(io,'(" INITIAL NUMBER OF ATOM/BASIS:",I5,I6)') J0,JF0
         endif
         J0=J1
         if (io>0) then
            write(io,'("   FINAL NUMBER OF ATOM/BASIS:",I5,I6)') J0,JF
            write(io,'(" Final",I5," atoms in the AO domain:")') J0
            write(io,'(10I5)') BA(:J0)
         endif

C ***** Constructing the AAO domain of each atom in the AO domain *****
         allocate(J1_aao(J0),JF_aao(J0),BA_AAO(NATOM,J0))
         allocate(co1_aao(JF1,JF1))
         maxdis=7
123      do j=1,J0
            do kzhi=2,maxdis
               do k=1,NATOM
                  BA_AAO(k,j)=atmclu(k,BA(j),kzhi) !construct AAO domain of atm j in the AO domain
               end do
               kb_aao=nbs_atm(BA_AAO(1,j)) ! number of PAOs of this atom
               if (io>0) then
                  write(io,*) "***********************************"
                  write(io,'(" Distance thresh: ",I4," Angstrom")') kzhi
               end if
               J1_aao(j)=natmclu(BA(j),kzhi)
               if (io>0) then
                  write(io,*) "Atomic AO domain",j
                  write(io,*) "Central atom of the AAO:",BA_AAO(1,j)
                  write(io,*) "Num of atoms in the AAO:",J1_aao(j)
               end if
               call BASIS_REORDER(io,JF1,NATOM,ZA,nb,BA_AAO(:,j),
     &                            J1_aao(j),JF_aao(j))
               do k=1,kb_aao
                  co1_aao(:,k)=PAO(:,nb(BA_AAO(1,j))+k-1)
               end do
               allocate(S4(JF_aao(j),JF_aao(j)))
               call INDEX_REORDER(io,JF1,JF_aao(j),kb_aao,
     &                            S,S2,S4,co1_aao(:,1:kb_aao),ZA)
               if (io>0) write(io,*) "Orbitals to be projected:",kb_aao
               allocate(coaao(JF_aao(j),kb_aao),wf1(kb_aao),wf2(kb_aao))
               call projorb(io,JF1,JF_aao(j),kb_aao,kb_aao,
     &                     co1_aao(:,1:kb_aao),coaao,S2,wf1,wf2)
               call locindx_occ(io,JF1,co1_aao(:,1:kb_aao),JF_aao(j),
     &                          coaao,kb_aao,S4,wf1,virthre,nstar)
               
               if (kzhi==maxdis .and. nstar<kb_aao) then
                  write(iout,*) "Warning: AAO domain over",
     &                           maxdis,"Angstroms!"
                  write(iout,*) "Reduce the redundancy:",kb_aao-nstar
                  if (allocated(paored_full)) then
                     allocate(paored(kb_aao))
                  else
                     allocate(paored(kb_aao),paored_full(JF1))
                     paored_full=0
                  end if
                  call locindx_red(io,JF1,co1_aao(:,1:kb_aao),JF_aao(j),
     &                      coaao,kb_aao,S4,wf1,virthre,nstar,paored)
                  paored_full(nb(BA_AAO(1,j)):nb(BA_AAO(1,j))+kb_aao-1)
     &            =paored(:kb_aao)
                  deallocate(paored)
                  if (io>0) then
                     write(io,'(" Final distance:",I2," A")') kzhi
                     write(io,*) "Atoms in the AAO domain:"
                     write(io,"(10I5)") atmclu(:J1_aao(j),BA(j),kzhi)
                     write(io,*) "***********************************"
                  end if
                  write(iout,*)
               end if
               deallocate(wf1,wf2,coaao,S4)
               if (nstar==kb_aao) then
                  if (io>0) then
                     write(io,'(" Final distance:",I2," A")')  kzhi
                     write(io,*) "Atoms in the AAO domain:"
                     write(io,"(10I5)") atmclu(:J1_aao(j),BA(j),kzhi)
                     write(io,*) "***********************************"
                     write(io,*)
                  endif
                  exit
               end if
            end do
         end do
         
         if (allocated(BA_FAO)) deallocate(BA_FAO)
         allocate(BA_FAO(NATOM))
         call Full_AAO(NATOM,J0,J1_aao,BA_AAO,BA_FAO,nat_fao)
         if (io>0) then
            write(io,*) "###### Final AO Domain #####"
            write(io,'(" Num of atoms in the FAO domain:",I4)') nat_fao
            write(io,*) "Atoms in the final AO domain No.",I
            write(io,'(10I5)') BA_FAO(:nat_fao)
         endif
      
         call BASIS_REORDER(io,JF1,NATOM,ZA,nb,BA_FAO,nat_fao,JF_fao)

         if (JF_fao>(JF1*90/100) .and. maxdis>6) then !NZG 2015/10/23 To avoid too large cluster
            maxdis=maxdis-1
            write(iout,*) "Warning: Number of basis functions over 90%."
            write(iout,'(" Message from cluster ",I5)') I
            goto 123
         end if
      
* ==============================================================
*   Project occ LMOs from the AO domain to the full AAO domain
* ==============================================================
         if (io>0) then
            write(io,*) "**** Begin to project occ LMOs ****"
            write(io,*) "KB=",KB
            write(io,*) "KB1=",KB1
         endif
         allocate(S4(JF_fao,JF_fao))
         call INDEX_REORDER(io,JF1,JF_fao,KB,S,S2,S4,co1(:,1:KB),ZA)
         allocate(cvs(JF_fao,KB),wf1(KB),wf2(KB))
         call projorb(io,JF1,JF_fao,KB,KB,co1(:,1:KB),cvs,S2,wf1,wf2)
         call locindx_occ(io,JF1,co1(:,1:KB),JF_fao,cvs,
     &                                 KB,S4,wf1,occthre,nstar)
         if (io>0) then
            write(io,'(" Number of occupied LMOs:",I5)') KB
            write(io,'(" Number of orbitals OK:  ",I5)') nstar
         endif
         deallocate(wf1,wf2,J1_aao,JF_aao,BA_AAO,co1_aao)
         if (KB/=nstar) then
            write(iout,*) "Initial AO domain is not large enough!"
            thre=thre+1.0D0
         else  
            exit
         end if
      end do     ! This end do is from "do while(1)" -*- Zhigang 2015.9.28
      
      call normorb(JF_fao,KB,cvs,S4)
      deallocate(S4,co1)
      MOS1(1:JF_fao,1:KB)=cvs(:,1:KB)
      deallocate(cvs)
      if (io>0) then
         write(io,*) "*************************************"
         write(io,*) "Finished constructing Proj. occ LMOs!"
         write(io,*) "*************************************"
      endif

* =============================================================
*    Project PAOs from the AO domain to the full AAO domain
* =============================================================
      if (allocated(paored_full)) then
         allocate(PAO_tmp(JF1,JF))
         JF_red=0
         do k=1,JF
            if (paored_full(ZA_AO(k))==0) then
               JF_red=JF_red+1
               PAO_tmp(:,JF_red)=PAO(:,ZA_AO(k))
            end if
         end do
         JF=JF_red
         allocate(co1(JF1,JF))
         co1(:,:)=PAO_tmp(:,:JF)
         deallocate(PAO_tmp,paored_full)
      else
         allocate(co1(JF1,JF))
         do k=1,JF
            co1(:,k)=PAO(:,ZA_AO(k))
         end do
      end if
      deallocate(ZA_AO)
      allocate(S4(JF_fao,JF_fao))
      call INDEX_REORDER(io,JF1,JF_fao,JF,
     &                       S,S2,S4,co1(:,1:JF),ZA)
      allocate(cvs(JF_fao,JF),wf1(JF),wf2(JF))
      call projorb(io,JF1,JF_fao,JF,JF,
     &                 co1(:,1:JF),cvs,S2,wf1,wf2)
      call locindx_occ(io,JF1,co1(:,1:JF),JF_fao,cvs,
     &                     JF,S4,wf1,virthre,nstar)
      if (io>0) then
         write(io,'(" Number of PAOs:       ",I5)') JF
         write(io,'(" Number of orbitals OK:",I5)') nstar
         write(io,*) "*************************************"
         write(io,*) "Finished constructing Projected PAOs!"
         write(io,*) "*************************************"
      endif
      deallocate(wf1,wf2)
      
      
      allocate(SS(JF,JF),Q(JF,JF),e(JF),d(JF))
      call MKL_Tfock(0,JF_fao,JF,S4,cvs(:,1:JF),SS)
      call NJ_qr(io,SS,JF,Q,e,d,k,1)
      call numpoint(JF,e,-zeps,zeps,nzero)
      nonzero=JF-nzero; nvext=nonzero; ncut =nzero
C Set MOS1(JF1,JF1) as a zero matrix
      deallocate(SS)
      allocate(cvs2(JF_fao,nvext),tr(JF,nvext))

      cvs2(:,:)=0.0D0
      call NJ_canorth(-1,JF,nvext,e(ncut+1),q(1,ncut+1),tr)
      deallocate(q,e,d)
      P1=1.0D0;P2=0.0D0
      call DGEMM('N','N',JF_fao,nvext,JF,P1,cvs(:,1:JF),
     &                      JF_fao,tr,JF,P2,cvs2,JF_fao)
      call normorb(JF_fao,nvext,cvs2,S4)
      deallocate(tr)
      JF=JF_fao; J0=nat_fao; BA=BA_fao; NA=KB+nvext
      do k=KB+1,NA
         do j=1,JF
            MOS1(j,k)=cvs2(j,k-KB)
         enddo
      enddo

C If doing CIM gradient calculation, we need to orthogonalize the
C orbitals.
C NZG_4/19/2017
C      if (calforce) call Schmidt_orth(JF,NA,MOS1(:JF,1:NA),S4)

      deallocate(cvs2,S4)
      deallocate(cvs)

      return
      end subroutine COV_FPAO


C ************************************************************
C * Constructing virtual space of the clusters (PAO basis)   *
C * In this way, PAOs are constructed within a cluster, by   *
C * projecting out the occupied space of the cluster.        *
C * NZG_4/7/2017 @@UARK                                       *
C ************************************************************
      subroutine COV_CPAO(io,I,natom,nmo,nb,J0,JF,JF1,KB,ncen,KB1,
     &                    ncoreocc,nvext,ML,BA,ZA,S,S2,SMO,MOS1,dis,
     &                    link,nbs_atm,bas_atm)
      implicit none
      integer i,j,k,jzhi,kzhi,l,m,maxdis,nstar,nonzero,ncoreocc,ncore1
      integer io,iout,NATOM,Nmo,J0,JF0,JF,KB,KB1,J1,JF1,npao,na
      integer nzero,ncut,nvext,ll,ini,ifi,np,ierr,paoout,ncen,ncen1
      integer ML(KB1),BA(NATOM),ZA(JF1),nb(NATOM),link(NATOM,NATOM)
      integer nbs_atm(NATOM),bas_atm(JF1,NATOM)
      integer,dimension(:),allocatable::BAA,ZA_AO
      real(kind=8),parameter::zeps=1.0d-5
      real(kind=8),parameter::occthre=0.01D0
      real(kind=8) thre,P1,P2,eps,epst,cvszhi
      real(kind=8) S(JF1,JF1),SMO(JF1,nmo),MOS1(JF1,JF1)
      real(kind=8) dis(NATOM,NATOM),S2(JF1,JF1)
      real(kind=8),dimension(:),allocatable::MOSS,wf1,wf2,e,d
      real(kind=8),dimension(:,:),allocatable::S4,tr,cvs,SS,Q,c
      real(kind=8),dimension(:,:),allocatable::smo1,PAO,MOSSS,Sij
      real(kind=8),dimension(:,:),allocatable::c2,cvs2

C KB:            Number of occ MOs in the cluster
C KB1:           Number of occ MOs of the whole molecule
C J0:            Number of atoms in the cluster
C JF:            Number of basis functions in the cluster
C JF1:           Number of basis functions of the whole molecule
C S2:            Overlap of basis in the cluster at the beginning
C S4:            Overlap of basis in the cluster
C ncoreocc:      Number of core occ LMOs in the whole molecule
C ncore1:        Number of core occ LMOs in the cluster

      iout=6

C Get the number of core orbitals in the cluster and rearrange the core
C orbitals.
C For Schmidt orthogonalization, the order should be
C central valence -> beffer valence -> core
C NZG_4/18/2017 @@UARK
      ncore1=0
      do k=1,KB
         if (ML(k)<=ncoreocc) ncore1=ncore1+1
      enddo
      j=0; l=0; ncen1=0
      do k=1,KB
         if (ML(k)<=ncoreocc) then
            MOS1(:,KB-j)=SMO(:,ML(k))
            j=j+1
         else
            if (k<=ncen) then
               ncen1=ncen1+1
               MOS1(:,ncen1)=SMO(:,ML(k))
            else
               MOS1(:,KB-ncore1-l)=SMO(:,ML(k))
               l=l+1
            endif
         endif
      enddo
      write(iout,*) "nunber of central orbitals:",ncen,ncen1
      ncen=ncen1

      if (io>0) write(iout,'("Number of core orbitals in the cluster:",
     &          6i)') ncore1

      allocate(S4(JF,JF))
      call INDEX_REORDER(io,JF1,JF,KB,S,S2,S4,MOS1(:,:KB),ZA)
      allocate(cvs(JF,KB),wf1(KB),wf2(KB))
      call projorb(io,JF1,JF,KB,KB,MOS1(:,1:KB),cvs,S2,wf1,wf2)
      call locindx_occ(io,JF1,MOS1(:,1:KB),JF,cvs,
     &                                    KB,S4,wf1,occthre,nstar)
      deallocate(wf1,wf2)
      J1=J0; thre=2.0D0; JF0=JF; npao=JF
      allocate(BAA(JF1))
      BAA=BA
      do while(nstar/=KB)
         BA=BAA
C         J0=nat_init    !Zhigang 2016.3.1 
         allocate(MOSS(JF1))
         do k=1,KB
            MOSS(:)=MOS1(:,k)
            do j=1,JF1
               MOS1(ZA(j),k)=MOSS(j)
            end do
         end do
         deallocate(MOSS,S4)
         call Add_Atom(io,KB,NATOM,link,dis,BA,J0,J1,thre)
         call BASIS_REORDER(io,JF1,NATOM,ZA,nb,BA,J1,JF)      
         allocate(S4(JF,JF))
         call INDEX_REORDER(io,JF1,JF,KB,S,S2,S4,MOS1(:,1:KB),ZA) 
         deallocate(cvs)
         allocate(cvs(JF,KB),wf1(KB),wf2(KB))
         call projorb(io,JF1,JF,KB,KB,MOS1(:,1:KB),cvs,S2,wf1,wf2)
         call locindx_occ(io,JF1,MOS1(:,1:KB),JF,cvs,KB,S4,wf1,
     &                    occthre,nstar)
         deallocate(wf1,wf2)
         thre=thre+1.0D0
         if (JF==JF1) then
            write(iout,*) " Warning: The AO domain is the whole system!"
            exit
         endif
      enddo

      if (io>0) then
1000      write(io,'(" INITIAL NUMBER OF ATOM/BASIS:",I5,I6)') J0,JF0
      endif
      J0=J1
      if (io>0) then
         write(io,'("   FINAL NUMBER OF ATOM/BASIS:",I5,I6)') J0,JF
         write(io,'(" Final",I5," atoms in the AO domain:")') J0
         write(io,'(10I5)') BA(:J0)
      endif

C Normalize the orbitals
      call normorb(JF,KB,cvs,S4)

C Check the orthogonality
C      allocate(Sij(KB,KB))
C      call MO_over(Sij,cvs,S4,KB,JF)
C      write(6,*) "KB",KB
C      write(6,*) Sij
C      stop
C The results show that the trancated occ LMOs are weak orthogonal to 
C each other. The lagest overlap is usually 10^-4
C Anyway, orthogonalize the orbitals

C      allocate(SS(KB,KB),Q(KB,KB),e(KB),d(KB))
C      call MKL_Tfock(0,JF,KB,S4,cvs,SS)
C      call NJ_qr(io,SS,KB,Q,e,d,k,1)
C      write(6,*) "after qr"
C      call numpoint(KB,e,-zeps,zeps,nzero)
C      write(6,*) "nzero",nzero
C      nonzero=KB-nzero; nvext=nonzero; ncut =nzero
C      if (ncut/=0) then
C         write(iout,*) " Linear dependance in occ LMOs!"
C         stop
C      endif
C      deallocate(SS)
C      write(6,*) "non-linear"
C      allocate(cvs2(JF,nvext),tr(KB,nvext))
C
C      cvs2(:,:)=0.0D0
C      call NJ_canorth(-1,KB,nvext,e(ncut+1),q(1,ncut+1),tr)
C      deallocate(q,e,d)
C      P1=1.0D0;P2=0.0D0
C      call DGEMM('N','N',JF,nvext,KB,P1,cvs,JF,tr,KB,P2,cvs2,JF)
C      call normorb(JF,nvext,cvs2,S4)
C      deallocate(tr,cvs)

C To keep the central and valence space of the occupied orbitals,
C Schmidt orthogonalization should be used.
      call Schmidt_orth(JF,KB,cvs,S4)
      MOS1(1:JF,1:(KB-ncore1))=cvs(1:JF,1:(KB-ncore1))

C Construct PAOs within this cluster by projecting out the occupied
C space of this cluster.
C -NZG_4/12/2017 @@UARK
      
      allocate(PAO(JF,JF))
      call dunit(JF,PAO)
      call pjotorb(JF,KB,JF,cvs,PAO,S4)
      call normorb(JF,JF,PAO,S4)
      KB=KB-ncore1
      deallocate(cvs)
 
C Orthogonalize the PAOs
      allocate(SS(npao,npao),Q(npao,npao),e(npao),d(npao))
      call MKL_Tfock(0,JF,npao,S4,PAO(:,1:npao),SS)
      call NJ_qr(io,SS,npao,Q,e,d,k,1)
      call numpoint(npao,e,-zeps,zeps,nzero)
      nonzero=npao-nzero; nvext=nonzero; ncut =nzero
      write(iout,*) 'JF,nvext',npao,nvext
      deallocate(SS)
      allocate(cvs2(JF,nvext),tr(npao,nvext))
     
      cvs2(:,:)=0.0D0
      call NJ_canorth(-1,npao,nvext,e(ncut+1),q(1,ncut+1),tr)
      deallocate(q,e,d)
      P1=1.0D0;P2=0.0D0
      call DGEMM('N','N',JF,nvext,npao,P1,PAO(:,1:npao),JF,tr,npao,P2,
     &           cvs2,JF)
      call normorb(JF,nvext,cvs2,S4)
      deallocate(tr)
      
      MOS1(1:JF,KB+1:KB+nvext)=cvs2
      deallocate(cvs2,PAO)

      return
      end subroutine COV_CPAO


C
C      ngrp_dis=0
C      allocate(natm_grp(natom),atm_grp(30,natom))
C      if (numc<=3) then
C         ngrp_dis=1
C         natm_grp(1)=NATOM
C         do i=1,NATOM
C            atm_grp(i,1)=i
C         end do
C      else
C         call divi_mole3(numc,natc,NATOM,dis,ngrp_dis,natm_grp,atm_grp)
C         allocate(natm_grp_tmp(natom),atm_grp_tmp(30,natom))
C         natm_grp_tmp=natm_grp      
C         atm_grp_tmp=atm_grp
C         loopi:do i=1,numh
C            mindis=100.0D0
C            do j=1,numc
C               if (dis(nath(i),natc(j))<mindis) then
C                  mindis=dis(nath(i),natc(j))
C                  minc=natc(j)
C               end if
C            end do
C            do jzhi=1,ngrp_dis
C               do kzhi=1,natm_grp(jzhi)
C                  if (minc==atm_grp(kzhi,jzhi)) then
C                     natm_grp_tmp(jzhi)=natm_grp_tmp(jzhi)+1
C                     atm_grp_tmp(natm_grp_tmp(jzhi),jzhi)=nath(i) 
C                     cycle loopi
C                  end if
C               end do
C            end do
C         end do loopi
C         natm_grp=natm_grp_tmp
C         atm_grp=atm_grp_tmp
C         deallocate(natm_grp_tmp,atm_grp_tmp)
C      end if
C      write(io,*) "------- GROUPS OF THE MOLECULE -------"
C      write(io,*) "INDEX  ATOMS   ATOM_INDX"
C      do i=1,ngrp_dis
C         call NJ_prtlab(line,natm_grp(i),atm_grp(:,i))
C         call NJ_trim(line,k1,k2)
C         write(io,'(2x,I3,3x,I4,4x,a)') i,natm_grp(i),line(k1:k2)
C      end do
C      write(io,*) "--------------------------------------"
C      write(io,*) 
C
CC === Assign each atom to its group ===
C      allocate(atm2grp(NATOM))
C      do i=1,ngrp_dis
C         do j=1,natm_grp(i)
C            atm2grp(atm_grp(j,i))=i
C         end do
C      end do
C
CC === From atoms to orbitals ===
C      allocate(norb_atm(NATOM),orb_atm(NUW-ncore,NATOM))
C      norb_atm=0; orb_atm=0
C      do i=ncore+1,NUW
C         norb_atm(SOB(i,1))=norb_atm(SOB(i,1))+1
C         orb_atm(norb_atm(SOB(i,1)),SOB(i,1))=i
C      end do
C
CC === From groups to orbitals ===
C      allocate(norb_grp(ngrp_dis),orb_grp(NUW-ncore,ngrp_dis))
C      norb_grp=0; orb_grp=0
C      do i=1,ngrp_dis
C         do j=1,natm_grp(i)
C            izhi=norb_grp(i)
C            norb_grp(i)=norb_grp(i)+norb_atm(atm_grp(j,i))
C            orb_grp(izhi+1:norb_grp(i),i)
C     &      =orb_atm(:norb_atm(atm_grp(j,i)),atm_grp(j,i))
C         end do
C      end do
C
CC === From orbitals to groups ===
C      allocate(orb2grp(NUW))
C      orb2grp=0
C      do i=1,ngrp_dis
C         do j=1,norb_grp(i)
C            orb2grp(orb_grp(j,i))=i
C         end do
C      end do
C
CC === Create distance matrix between groups ===
C      allocate(dis_grp(ngrp_dis,ngrp_dis),min_dis(2,ngrp_dis,ngrp_dis))
C      do i=1,ngrp_dis
C         dis_grp(i,i)=0.0D0
C      end do
C      do i=1,ngrp_dis-1
C         do j=i+1,ngrp_dis
C            dis_grp(i,j)=dis(atm_grp(1,i),atm_grp(1,j))
C            dis_grp(j,i)=dis_grp(i,j)
C            min_dis(1,i,j)=atm_grp(1,i)
C            min_dis(1,j,i)=min_dis(1,i,j)
C            min_dis(2,i,j)=atm_grp(1,j)
C            min_dis(2,j,i)=min_dis(2,i,j)
C            do izhi=1,natm_grp(i)
C               do jzhi=1,natm_grp(j)
C                  if (dis(atm_grp(izhi,i),atm_grp(jzhi,j))<dis_grp(i,j))
C     &                                                              then
C                     dis_grp(i,j)=dis(atm_grp(izhi,i),atm_grp(jzhi,j))
C                     dis_grp(j,i)=dis_grp(i,j)
C                     min_dis(1,i,j)=atm_grp(izhi,i)
C                     min_dis(1,j,i)=min_dis(1,i,j)
C                     min_dis(2,i,j)=atm_grp(jzhi,j)
C                     min_dis(2,j,i)=min_dis(2,i,j)
C                  end if
C               end do
C            end do
C         end do
C      end do
C      write(io,*)
CC      write(io,*) "=== DISTANCES BETWEEN GROUPS ==="
CC      write(io,*) "GROUP_PAIR  ATOM_PAIR  DISTANCE"
CC      do i=1,ngrp_dis
CC         do j=1,ngrp_dis
CC            write(io,'(2x,I3,2x,I3,4x,I3,2x,I3,4x,F6.3)') 
CC     &      i,j,min_dis(1,i,j),min_dis(2,j,i),dis_grp(i,j)
CC         end do
CC      end do
CC      write(io,*) "================================"
C      
C      ngroup=NUW-ncore
C      nsys=nfrg
C      allocate(SYS(NUW,10000),ISY(10000),ISY2(10000))
CC
C      ISY=0
C      do i=1,nfrg
C         SYS(1:NUW-ncore,i)=FRG(:,i)
C         ISY2(i)=NF0(i)
C         ISY(i)=NF1(i)
C      end do
C      deallocate(FRG)
C
CC ==== Represent the MO domain by atoms ====
C      allocate(SYS_ATOM(NATOM,10000),SYS_tmp(NUW,nsys))
C      allocate(ISY_ATOM(10000))   ! There are ISY_ATOM(nsys) atoms in the subsystem.
C      SYS_ATOM=0; SYS_tmp=0; ISY_ATOM=0;
C      do i=1,nsys
C         do j=1,ISY(i)
C            SYS_tmp(j,i)=SOB(SYS(j,i),1)
C         end do
C      end do
C      do i=1,nsys
C         SYS_ATOM(1,i)=SYS_tmp(1,i)
C         ISY_ATOM(i)=1
C         do j=2,ISY(i)
C            do k=1,NATOM
C               if (SYS_tmp(j,i)==SYS_ATOM(k,i)) then
C                  inarray1=1
C                  exit
C               else
C                  inarray1=0
C               end if
C            end do
C            if (inarray1==0) then
C               ISY_ATOM(i)=ISY_ATOM(i)+1
C               SYS_ATOM(ISY_ATOM(i),i)=SYS_tmp(j,i)
C            end if
C         end do
C      end do
C
CC      do i=1,nsys
CC         write(io,*) "SYS:",i,ISY_ATOM(i)
CC         write(io,'(15I4)') (SYS_ATOM(j,i),j=1,ISY_ATOM(i))
CC      end do
CC
CC ==== Use subroutine from Jia Junteng to generate derivative subsystems ====      
C      maxatom=ISY_ATOM(1) ! Max number of atoms in the subsystems
C      do i=2,nsys
C         if (ISY_ATOM(i)>maxatom) maxatom=ISY_ATOM(i)
C      end do
C      nprim=nsys
C      ENTYP='CENT'
C      if (ENTYP=='LCCE') then
C      allocate(coef_GEBF(10000),SYS_ATOMtmp(maxatom,10000))
C      SYS_ATOMtmp=0
C      do i=1,nsys
C         SYS_ATOMtmp(:,i)=SYS_ATOM(1:maxatom,i)
C      end do
C      coef_GEBF=0.0D0
C      coef_GEBF(1:nsys)=1.0D0
C      write(io,*) "nsys",nsys
C      write(io,*) "maxatom",maxatom
C      if (maxatom<=10) then
C        call gebf_deri_s_cpp(nsys,maxatom,SYS_ATOMtmp(1,1),coef_GEBF(1))
C      else if (maxatom<=20) then
C        call gebf_deri_m_cpp(nsys,maxatom,SYS_ATOMtmp(1,1),coef_GEBF(1))
C      else 
C        call gebf_deri_l_cpp(nsys,maxatom,SYS_ATOMtmp(1,1),coef_GEBF(1))
C      end if
C      
C      do i=1+nprim,10000
C         if (SYS_ATOMtmp(1,i)/=0) then
C            nsys=nsys+1
C         else
C            exit
C         end if
C      end do
C      write(io,*) "Subsystem before overlap:",nprim
C      write(io,*) "Subsystem after  overlap:",nsys
C      do i=nprim+1,nsys
C         SYS_ATOM(1:maxatom,i)=SYS_ATOMtmp(:,i)
C         do j=1,maxatom
C            if (SYS_ATOMtmp(j,i)/=0) then
C               ISY_ATOM(i)=ISY_ATOM(i)+1
C            else
C               exit
C            end if
C         end do
C      end do
C      single_clu=0
C      do i=1,nsys
C         write(io,*) "Deri_SYS:",i,ISY_ATOM(i)
C         if (ISY_ATOM(i)==1) then
C            write(io,*) "Single atom cluster exists!"
C            single_clu=0
C            exit
C         end if
C         write(io,'(15I4)') (SYS_ATOM(j,i),j=1,ISY_ATOM(i))
C      end do
C      deallocate(SYS_ATOMtmp)
C      
C      if (single_clu==1) then
C        thre_grp=2.5D0
C        allocate(SYS_ini(ngrp_dis,ngrp_dis))
C        deallocate(NF0,NF1)
C        allocate(NF0(ngrp_dis),NF1(ngrp_dis))
C        call Frag_Deter(io,ngrp_dis,dis_grp,norb_grp,SYS_ini,
C     &                  NF0,NF1,thre_grp,ngrp_di)
C        allocate(SYS_orb(ngrp_di,ngrp_di))
C        do i=1,ngrp_di
C           do j=1,ngrp_di
C              SYS_orb(i,j)=SYS_ini(i,j)
C           end do
C        end do
C        deallocate(SYS_ini)
C        call Reduce_Frag(io,ngrp_di,0,SYS_orb,NF0,NF1,nprim_dis)
C        maxgrp=NF1(1)
C        nsys=nprim_dis
C        allocate(SYS_grptmp(maxgrp,10000))
C        SYS_grptmp=0
C        do i=1,nsys
C          SYS_grptmp(:,i)=SYS_orb(:maxgrp,i)
C        end do
C        coef_GEBF=0.0D0
C        coef_GEBF(1:nsys)=1.0D0
C        if (maxgrp<=10) then
C          call gebf_deri_s_cpp(nsys,maxgrp,SYS_grptmp(1,1),coef_GEBF(1))
C        else if (maxgrp<=20) then
C          call gebf_deri_m_cpp(nsys,maxgrp,SYS_grptmp(1,1),coef_GEBF(1))
C        else
C          call gebf_deri_l_cpp(nsys,maxgrp,SYS_grptmp(1,1),coef_GEBF(1))
C        end if
C        
C        do i=1+nprim_dis,10000
C          if (SYS_grptmp(1,i)/=0) then
C            nsys=nsys+1
C          else
C            exit
C          end if
C        end do
C
C        write(io,*) "Subsystem before overlap:",nprim_dis
C        write(io,*) "Subsystem after  overlap:",nsys
C        
C        allocate(ISY_GRP(nsys),SYS_GRP(maxgrp,nsys))
C        ISY_GRP=0; SYS_GRP=0
C        do i=1,nprim_dis
C          ISY_GRP(i)=NF1(i)
C          SYS_GRP(:ISY_GRP(i),i)=SYS_orb(:ISY_GRP(i),i)
C        end do 
C        do i=nprim_dis+1,nsys
C          SYS_GRP(:,i)=SYS_grptmp(:,i)
C          do j=1,maxgrp
C            if (SYS_grptmp(j,i)/=0) then
C               ISY_GRP(i)=ISY_GRP(i)+1
C            else
C               exit
C            end if
C          end do
C        end do
C        
C        do i=1,nsys
C          write(io,*) "Deri_SYS:",i,ISY_GRP(i)
C          write(io,'(15I4)') (SYS_GRP(j,i),j=1,ISY_GRP(i))
C        end do
C
CC === Transform the central groups to central LMOs ===        
C        ISY2=0
C        do i=1,nprim_dis
C          do j=1,NF0(i)
C            do k=ncore+1,NUW
C              if (orb2grp(k)==SYS_GRP(j,i)) then
C                 ISY2(i)=ISY2(i)+1
C              end if
C            end do
C          end do
C        end do
C
CC === Transform the groups to atoms ===
C        SYS_ATOM=0;ISY_ATOM=0
C        do i=1,nsys
C          do j=1,ISY_GRP(i)
C            do k=1,NATOM
C              if (atm2grp(k)==SYS_GRP(j,i)) then
C                 ISY_ATOM(i)=ISY_ATOM(i)+1
C                 SYS_ATOM(ISY_ATOM(i),i)=k
C              end if
C            end do
C          end do
C        end do
C
C      end if
C
CC ==== Transform the atoms to LMOs ====
C      ISY=0;SYS=0
C      do i=1,nsys
C         do j=1,ISY_ATOM(i)
C            do k=ncore+1,NUW
C               if (SOB(k,1)==SYS_ATOM(j,i)) then
C                  ISY(i)=ISY(i)+1
C                  SYS(ISY(i),i)=k
C               end if
C            end do
C         end do
C      end do
C      
C      coefout=308
C      call SEQOPN(COEFOUT,'COEFOUT','UNKNOWN',.FALSE.,'FORMATTED')
C      do i=1,nsys
C         write(COEFOUT,*) coef_GEBF(i)
C      end do
C      deallocate(coef_GEBF)
CC     ISY2(nprim+1:nsys)=1  ! For the derivative subsystems, we assign the first LMO to be central MO.
C      if (single_clu==0) then
C         ISY2(nprim+1:nsys)=ISY(nprim+1:nsys) ! NZG-2015.4.16 Assign all LMOs as central MOs.
C      else
C         ISY2(nprim_dis+1:nsys)=ISY(nprim_dis+1:nsys)   ! NZG-2015.5.10 For the group part, all LMOs are central MOs.
C      end if                                            ! NZG-2015.5.25 For the group part, we have central MOs.
C      
C      end if   ! This if is from if (ENTYP=LCCE) then ... -*- NZG-2015.9.3
C
C      deallocate(FXY)
C      ngroup=nsys
C      ngat = ngroup
C      allocate(FRG(NUW,ngroup))
C      allocate(EF(NUW,ngroup),FF(ngroup,ngroup))
C      allocate(gdis(ngroup,ngroup),glink(ngroup,ngroup))
C      allocate(GSYS(ngroup,ngroup),NG0(ngroup),NG1(ngroup))
C      allocate(GCEN(ngroup,ngroup))
C      allocate(grpcms(3,ngroup),twocms(3,ngroup,ngroup))
C      allocate(labgmo(ngroup))
C
CC
C      allocate(INF2(NUW,KG),BA(NATOM,KG),ZA(NW,KG),nat_sys(KG))
C      allocate(NASUB(KG),NWSUB(KG))
C      INF2=0;BA=0; ZA=0; NASUB=0; NWSUB=0
C
C      call timit(1)
C      write(io,*)
CC
C      if (IWORK(7).eq.1) IWORK(7)=0  !-WL- 2009.10.27 ADD for read occ MO only
CC
C      allocate(ISUB5(KGG),JOBX(KGG))
C      allocate(mbrams(KGG),xtime(KGG))
C      allocate(scftyps(KGG))
C      call VICLR(ISUB5, 1, KGG)
C      call VICLR(JOBX,  1, KGG)
C      call VICLR(mbrams,1, KGG)
C      call VCLR(xtime, 1, KGG)
C      call VCCLR(scftyps, 1, KGG)
CC
CC --- Calculations of each subsystems (1 - KGG) ---
C      write(io,*) '--- GENERATING SUBSYSTEMS FROM CIM DOMAINS ---'
C      maxbs  = 0
CC --- Modified By Guoyang -------
C      Nvi=Nmo-NUW
C      allocate(MOS1(NW,Nvi))
C      MOS1(:,1:Nvi)=SMO(:,NUW+1:Nmo)
C      epsa=1d-6; epst=1d-12; maxcyc=1000000
C      nprint = 0
C      np     = 1
C      nprtcyc= 20
C
C      SMO(:,NUW+1:Nmo)=MOS1(:,:)
CC
C      call rwrit(GEBFINP,'$VEC', NW*Nmo, SMO(1,1))
CC
C      CHUNK=10
C
C      KGG2=0
C      do I=1,KGG
C         if (mtdsys2(I).ge.1) KGG2=KGG2+1
C      enddo
CC
C      write(io,200)
C      NSYGEN = 0
C      memtot = 0
C      timtot = 0.0D+00
C      NERR = 0
C      nomax = 0
C      numax = 0
C      nbsmax = 0
C      imax = 0
C      jmax = 0
C      kmax = 0
C      do I=1,KGG
C         L=mtdsys2(I)  ! 0: NONE; >0: MP/CC; <0 MP/CC with existed input
CC        if (L.le.0) cycle ! Skip NONE and existed subsystems
C         METHOD = METHDS(abs(L))
C         if (L.gt.0) then
C            NSYGEN = NSYGEN + 1
C            JOBX(NSYGEN) = I
C         elseif (L.eq.0) then
C            cycle
C         endif
CC
C         if (scftyps(I).eq.'ROHF    ') then
C            if (METHOD.ne.'CCSD    '.and.METHOD.ne.'CR-CCL') NERR=NERR+1
C         endif
CC
C         memtot=memtot+mbrams(I)
CC
C         nc  = ISY(I)
C         no  = ISY(I)
C         nu  = ISUB5(I)
C         nbs = NWSUB(I)
C         nu4 = nu**4
CC
C         if (no.gt.nomax) then
C            nomax=no
C            imax=i
C         endif
CC
C         if (nu.gt.numax) then
C            numax=nu
C            jmax=i
C         endif
CC
C         if (nbs.gt.nbsmax) then
C            nbsmax=nbs
C            kmax=i
C         endif
CC
C         xo1w4=1.0D-08*dble(nbs*nbs)
C         xo1w4=xo1w4*nbs*nbs*no*1.0D-03
C         xo2u4=1.0D-08*dble(nu4)
C         xo2u4=xo2u4*no*no*1.0D-03
C         xo3u4=xo2u4*no*0.5D+00
CC
C         if (METHOD.eq.'MP2     ' .or. METHOD.eq.'RIMP2   ') then
C            xtime(I) = xo1w4
C         elseif (METHOD.eq.'CCD     ') then
C            xtime(I) = xo2u4*15+xo1w4
C         elseif (METHOD.eq.'CCSD    ') then
C            xtime(I) = xo2u4*30+xo1w4
C         elseif (METHOD.eq.'CCSD(T) ') then
C            xtime(I) = xo3u4*10+xo2u4*30+xo1w4
C         elseif (METHOD.eq.'CR-CCL  ') then
C            xtime(I) = xo3u4*30+xo2u4*60+xo1w4
C         endif
CC
C         if (scftyps(I).eq.'ROHF    ') xtime(I)=xtime(I)*3.0D+00
C         timtot = timtot + xtime(I)
CC
C         if (L.gt.0) then
C            write(io,201) I,scftyps(I),METHOD,nc,no,nu,nat_sys(I),
C     &                    nbs,mbrams(I),xtime(I)
C         else
C            write(io,202) I,METHOD,nc,no,nbs
C         endif
C      enddo
C      write(io,205) memtot,timtot
CC
C      write(io,'('' Max(no) ='',i6,'' in sys'',i5)') nomax, imax
C      write(io,'('' Max(nu) ='',i6,'' in sys'',i5)') numax, jmax
C      write(io,'('' Max(nbs)='',i6,'' in sys'',i5)') nbsmax,kmax
CC
CC --- Get required memory and estimated CPU time for total system
C      yo1w4 = 1.0D-08*dble(NW*NW)
C      yo1w4 = yo1w4*NW*NW*NUTT*1.0D-03
C      yo2u4 = 1.0D-08*dble(NVTT*NVTT*NVTT*NVTT)
C      yo2u4 = yo2u4*NUTT*NUTT*1.0D-03
C      yo3u4 = yo2u4*NUTT*0.5D+00/3.0D+00
CC
C      if (scftyp.eq.'RHF     ') then
C         k1 = memgms(scftyp,2,'NONE    ',NW,NUTT,NVTT,0)
C         k2 = memgms(scftyp,0,'CCD     ',NW,NUTT,NVTT,0)
C         k3 = memgms(scftyp,0,'CCSD    ',NW,NUTT,NVTT,0)
C         k4 = memgms(scftyp,0,'CCSD(T) ',NW,NUTT,NVTT,0)
C         k5 = memgms(scftyp,0,'CR-CCL  ',NW,NUTT,NVTT,0)
C         k4 = max(k3,k4)
C         k5 = max(k3,k5)
C         MBRAM(1) = int(k1/125.d0/1000d0) + 1
C         MBRAM(2) = int(k2/125.d0/1000d0) + 1
C         MBRAM(3) = int(k3/125.d0/1000d0) + 1
C         MBRAM(4) = int(k4/125.d0/1000d0) + 1
C         MBRAM(5) = int(k5/125.d0/1000d0) + 1
C         ytime(1) = yo1w4
C         ytime(2) = yo2u4*15+yo1w4
C         ytime(3) = yo2u4*30+yo1w4
C         ytime(4) = yo3u4*10+yo2u4*30+yo1w4
C         ytime(5) = yo3u4*30+yo2u4*60+yo1w4
C         write(io,212) (scftyp,METHDS(k),NUTT,NVTT,NW,
C     &                 MBRAM(k),ytime(k),k=1,5)
C      elseif (scftyp.eq.'ROHF    ') then
C         k3 = memgms(scftyp,0,'CCSD    ',NW,NUTT,NVTT,0)
C         k5 = memgms(scftyp,0,'CR-CCL  ',NW,NUTT,NVTT,0)
C         k5 = max(k3,k5)
C         MBRAM(3) = int(k3/125.d0/1000d0) + 1
C         MBRAM(5) = int(k5/125.d0/1000d0) + 1
C         ytime(3) = (yo2u4*30+yo1w4)*3.0D+00
C         ytime(5) = (yo3u4*30+yo2u4*60+yo1w4)*3.0D+00
C         write(io,212) scftyp,METHDS(3),NUTT,NVTT,NW,MBRAM(3),ytime(3)
C         write(io,212) scftyp,METHDS(5),NUTT,NVTT,NW,MBRAM(5),ytime(5)
C      endif
C      write(io,215)
C      write(io,*)
CC
CC --- Sorting by the estimated CPU timing ---
C      do idx=1,NSYGEN-1
C         i = JOBX(idx)
C         do jdx=idx+1,NSYGEN
C            j = JOBX(jdx)
CC
C            time1=xtime(i)
C            time2=xtime(j)
C            if (time1.ge.time2) cycle
CC
C            k=JOBX(i)
C            JOBX(i)=JOBX(j)
C            JOBX(j)=k
CC
C            timex=xtime(i)
C            xtime(i)=xtime(j)
C            xtime(j)=timex
C         enddo
C      enddo
CC
C      write(io,*) 'Sorted subsystems by decreasing CPU time:'
C      write(io,224) (JOBX(I),I=1,NSYGEN)
C      write(io,*)
CC
C      if (NERR.NE.0) then
C         write(io,226) NERR
C         call abrt
C      endif
CC
C 200  format(1x,'=== Required memory (MB) and ',
C     &      'estimated CPU time (Relative) ==='/
C     &       1x,'SYS   SCFTYP  METHOD     NC   NO   NU  NATOM NBS  ',
C     &      'MEMORY    CPU TIME')
C 201  format(i4,3x,2a8,5i5,i8,f12.4)
C 202  format(i4,3x,'--',6x,a8,2i5,3x,'--',i5,6x,'--',6x,'SKIPED')
C 205  format(' Tot',3x,8x,8x,5x,15x,i8,f12.4)
C 212  format(' Can.  ',a8,a8,5x,3i5,i8,f12.4)
CC
C 215  format(1x,'Notes: The estimated CPU times are based on the',
C     &       ' following scalings:'/
C     &       1x,'MP2: O(NO*NBS^4)'/
C     &       1x,'CCD or CCSD: O(NO^2*NU^4)'/
C     &       1x,'CCSD(T) or CR-CC(2,3): O(NO^3*NU^4)'/
C     &       1x,'t(ROHF-CC)/t(RHF-CC) = 3'/
C     &       1x,'The CPU times between different methods are',
C     &       ' estimated empirically.'/
C     &       1x,'They may be inaccurate especially for small systems/',
C     &       'subsystems.')
C 224  format(1x,5i5,2x,5i5)
C 226  format(1x,'ERROR(S) IN',i4,' SUBSYSTEMS,',
C     *     ' ROHF CAN ONLY BE USED FOR CIM CCSD AND CR-CCL'/)
CC204  format(1x,'Tot',16x,f7.2,' (100.%)',f8.2,' (100.%)',i8)
CC
C      write(inp,*) '=== SOME INFOMATION FOR DOMAIN ==='
C      write(inp,'(''NUMSYS = '',2I8)') KGG,ngroup2
C      write(inp,'(6i7)') NATOM,NW,Nmo,NCORE,NUW,IWORK(6)
C      write(inp,*) 'KSymm()'
C      write(inp,'(15i6)') (KSymm(i),i=1,NUW)
C      write(inp,'(''NGROUP = '',2I8)') nprim,ngmo
C      write(inp,*)
CC
C      deallocate(KTWO,KTS,CTWO,DTWO)
C      deallocate(SR,SR1,SNB,SOB)
C      deallocate(INF,INF2,ISUB,ISUB2)
C      deallocate(BA,ZA,NB)
C      deallocate(NASUB,NWSUB,natc,nath,ED,dis,link,NCLU)
C      deallocate(ECCA,ECCB,ETCA,ETCB,NCLA,NCLB,KSymm,group,FRG,NF0,FF)
C      deallocate(EF)
C      deallocate(mbrams,ISUB5,JOBX)
C      deallocate(xtime)
C      deallocate(CenMO,molevl,mtdsys,FrgMO)
C      deallocate(FXY,scftyps)
CC
C      deallocate(gdis,glink)
C      deallocate(GSYS,GCEN,NG0,NG1)
C      deallocate(SYS,ISY,ISY2)
C      deallocate(ngrp)
C      deallocate(grpcms,twocms,labgmo)
CC
C      END
CC
CCC
Cc     ##############################################################
Cc     ##  subroutine NJ_prtcol -- print column(c1:c2) of mat(m,n) ##
Cc     ##  2005.12.22 by Wei Li; Update 2005.12.25 by Wei Li       ##
Cc     ##############################################################
Cc
Cc     format: f(x.y) x>7 sugg 12.5,12.7,14.9
CC     based on NJ_prtcol3b but only print m2 lines for each column
Cc
C      subroutine NJ_prtcol3b(io,m,n,mat0,mat,m2,c1,c2,fm)
C      implicit none
C      integer i,j,jj,m,n,io,c1,c2,n5,nf,nc,x,y,k,mat0(n,m),m2
C      real(kind=8) mat(m,n)
C      character fm*(*),ch,fm2*10
C      character*40 fmt1,fmt2,fmt3,fmt4
CC
C      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
C      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
C      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
CC
C      write(fmt1,101) ch,x,y; write(fmt2,102) nf,nf,ch,x,y
C 101  format('(i4,1x,5i4,5',a1,i2,'.',i2,')')
C 102  format('(i4,1x,',i2,'i4,',i2,a1,i2,'.',i2,')')
C      write(fmt3,103) x-7; write(fmt4,104) nf,nf,x-7
C 103  format('(5x,5i4,5(',i2,'x,i7))')
C 104  format('(5x,',i2,'i4,',i2,'(',i2,'x,i7))')
CC
C      do jj=1,n5
C         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1),
C     &                  (j,j=c1+(jj-1)*5,c1+jj*5-1)
C         write(io,fmt1)
C     &  (i,(mat0(j,i),j=c1+(jj-1)*5,c1+jj*5-1),
C     &      (mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=1,m2)
C         if (jj.ne.n5.or.nf.ne.0) write(io,'(1x,74(''-''))')
C      enddo
CC
C      if (nf.ne.0) then
C         write(io,fmt4)(j,j=c1+n5*5,c2),(j,j=c1+n5*5,c2)
C         write(io,fmt2) (i,(mat0(j,i),j=c1+n5*5,c2),
C     &                  (mat(i,j),j=c1+n5*5,c2),i=1,m2)
C      endif
CC
C      end
CC
CC
Cc     ##############################################################
Cc     ##  subroutine NJ_prtcol -- print column(c1:c2) of mat(m,n) ##
Cc     ##  2005.12.22 by Wei Li; Update 2008.08.06 by Wei Li       ##
Cc     ##############################################################
Cc
Cc     format: f(x.y) x>7 sugg 12.5,12.7,14.9
Cc
C      subroutine NJ_prtcol4(io,m,n,mat0,mat,r1,r2,c1,c2,fm)
C      implicit none
C      integer i,j,jj,m,n,io,c1,c2,n5,nf,nc,x,y,k,mat0(m,n),r1,r2
C      real(kind=8) mat(m,n)
C      character fm*(*),ch,fm2*10
C      character*40 fmt1,fmt2,fmt3,fmt4
CC
C      if (r1.lt.1.or.r2.gt.m) then
C         write(*,*) 'Wrong range for r1-r2'
C         write(*,*) 'switch to full row output'
C         r1=1
C         r2=m
C      endif
CC
C      nc=c2-c1+1; n5=nc/5; nf=mod(nc,5)
C      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
C      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
CC
C      write(fmt1,101) ch,x,y; write(fmt2,102) nf,nf,ch,x,y
C 101  format('(i4,1x,5i4,5',a1,i2,'.',i2,')')
C 102  format('(i4,1x,',i2,'i4,',i2,a1,i2,'.',i2,')')
C      write(fmt3,103) x-7; write(fmt4,104) nf,nf,x-7
C 103  format('(5x,5i4,5(',i2,'x,i7))')
C 104  format('(5x,',i2,'i4,',i2,'(',i2,'x,i7))')
CC
C      do jj=1,n5
C         write(io,fmt3) (j,j=c1+(jj-1)*5,c1+jj*5-1),
C     &                  (j,j=c1+(jj-1)*5,c1+jj*5-1)
C         write(io,fmt1)
C     &  (i,(mat0(i,j),j=c1+(jj-1)*5,c1+jj*5-1),
C     &      (mat(i,j),j=c1+(jj-1)*5,c1+jj*5-1),i=r1,r2)
C         if (jj.ne.n5.or.nf.ne.0) write(io,'(1x,74(''-''))')
C      enddo
CC
C      if (nf.ne.0) then
C         write(io,fmt4)(j,j=c1+n5*5,c2),(j,j=c1+n5*5,c2)
C         write(io,fmt2) (i,(mat0(i,j),j=c1+n5*5,c2),
C     &                  (mat(i,j),j=c1+n5*5,c2),i=r1,r2)
C      endif
CC
C      end
CC
Cc     format: f(x.y) x>7 sugg 12.5,12.7,14.9
Cc
CC
Cc     ##############################################################
Cc     ##  subroutine NJ_prtsym  --  print symmetric mat(n,n)      ##
Cc     ##  2005.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
Cc     ##############################################################
Cc
Cc     format: f(x.y) x>7 sugg 12.5,12.7,14.9
Cc
C      subroutine NJ_prtsym2(io,n,m,mat,fm)
C      implicit none
C      integer i,j,jj,m,n,io,n5,nf,nc,x,y,ini,ifi,k
C      real(kind=8) mat(m,m)
C      character fm*(*),ch,fm2*10
C      character*40 fmt1,fmt2,fmt3,fmt4
CC
C      n5=n/5; nf=mod(n,5)
C      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
C      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
CC
C      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
C 101  format('(i7,5',a1,i2,'.',i2,')')
C 102  format('(i7,',i2,a1,i2,'.',i2,')')
C      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
C 103  format('(3x,5(',i2,'x,i7))')
C 104  format('(3x,',i2,'(',i2,'x,i7))')
CC
C      do jj=1,n5
C         ini=1+(jj-1)*5
C         write(io,fmt3) (j,j=ini,jj*5)
C         do k=1+(jj-1)*5,n
C            ifi=min(jj*5,k)
C            write(io,fmt1) k,(mat(k,j),j=ini,ifi)
C         enddo
C      enddo
CC
C      if (nf.ne.0) then
C         ini=n-nf+1
C         write(io,fmt4)(j,j=ini,n)
C         do k=ini,n
C            write(io,fmt2) k,(mat(k,j),j=ini,k)
C         enddo
C      endif
CC
C      end
C
CC
Cc     ##############################################################
Cc     ##  subroutine NJ_prtsym  --  print symmetric mat(n,n)      ##
Cc     ##  2005.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
Cc     ##############################################################
Cc
Cc     format: f(x.y) x>7 sugg 12.5,12.7,14.9
Cc
C      subroutine NJ_prtsym5(io,n,m,mat,idx,fm)
C      implicit none
C      integer i,j,jj,m,n,io,n5,nf,nc,x,y,ini,ifi,k
C      real(kind=8) mat(m,m)
C      character fm*(*),ch,fm2*10
C      character*40 fmt1,fmt2,fmt3,fmt4
C      integer idx(m)
CC
C      n5=n/5; nf=mod(n,5)
C      fm2=fm; ch=fm2(1:1); k=index(fm2,'.')
C      read(fm2(2:k-1),*) x; read(fm2(k+1:10),*) y
CC
C      write(fmt1,101) ch,x,y; write(fmt2,102) nf,ch,x,y
C 101  format('(i7,5',a1,i2,'.',i2,')')
C 102  format('(i7,',i2,a1,i2,'.',i2,')')
C      write(fmt3,103) x-7; write(fmt4,104) nf,x-7
C 103  format('(3x,5(',i2,'x,i7))')
C 104  format('(3x,',i2,'(',i2,'x,i7))')
CC
C      do jj=1,n5
C         ini=1+(jj-1)*5
C         write(io,fmt3) (idx(j),j=ini,jj*5)
C         do k=1+(jj-1)*5,n
C            ifi=min(jj*5,k)
C            write(io,fmt1) idx(k),(mat(k,j),j=ini,ifi)
C         enddo
C      enddo
CC
C      if (nf.ne.0) then
C         ini=n-nf+1
C         write(io,fmt4)(idx(j),j=ini,n)
C         do k=ini,n
C            write(io,fmt2) idx(k),(mat(k,j),j=ini,k)
C         enddo
C      endif
CC
C      end
C
C
C
Cc     ##############################################################
Cc     ##  subroutine NJ_upper  --  make a string uppercase        ##
Cc     ##  2004.04.16 by Wei Li; Update 2005.10.17 by Wei Li       ##
Cc     ##############################################################
Cc
C      subroutine NJ_upper(line)
C      implicit none
C      integer i,k,ich
C      character line*(*),ch
C
C      k=len(line)
C      do i=1,k
C         ich=ichar(line(i:i))
C         if (ich.ge.97.and.ich.le.122)  then
C            line(i:i)=char(ich-32)
C         endif
C      enddo
C
C      end
C
Cc     ##############################################################
Cc     ##  function elemsyl  --  element: order --> symbol         ##
Cc     ##  2004.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
CC     ##  2015.07.05 by Zhigang for ghost atoms (nuclear=0)       ##
Cc     ##############################################################
CC
C      function elemsyl(n)
C      implicit none
C      integer n
C      character elemsyl*2,elem(110)*2
C
C      data elem/'H ','He','Li','Be','B ',  'C ','N ','O ','F ','Ne',
C     &          'Na','Mg','Al','Si','P ',  'S ','Cl','Ar','K ','Ca',
C     &          'Sc','Ti','V ','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',
C     &          'Ga','Ge','As','Se','Br',  'Kr','Rb','Sr','Y ','Zr',
C     &          'Nb','Mo','Tc','Ru','Rh',  'Pd','Ag','Cd','In','Sn',
C     &          'Sb','Te','I ','Xe','Cs',  'Ba','La','Ce','Pr','Nd',
C     &          'Pm','Sm','Eu','Gd','Tb',  'Dy','Ho','Er','Tm','Yb',
C     &          'Lu','Hf','Ta','W ','Re',  'Os','Ir','Pt','Au','Hg',
C     &          'Tl','Pb','Bi','Po','At',  'Rn','Fr','Ra','Ac','Th',
C     &          'Pa','U ','Np','Pu','Am',  'Cm','Bk','Cf','Es','Fm',
C     &          'Md','No','Lr','Rf','Db',  'Sg','Bh','Hs','Mt','Gh'/
C
C      if (n>0.and.n<110) then
C         elemsyl=elem(n); return
C      else if (n==0) then
C         elemsyl=elem(110); return
C      else
C         write(*,*) 'Wrong nuclear charges for',n; stop
C      endif
C
C      end
C
C****** 2004.12.22 element: symbol --> order ******
Cc     ##############################################################
Cc     ##  function elemord  --  element: symbol --> order         ##
Cc     ##  2004.12.22 by Wei Li; Update 2005.10.16 by Wei Li       ##
Cc     ##############################################################
CC
C      function elemord(a)
C      implicit none
C      integer elemord,i,j,k,l,ich
C      character a*(*),aa*2,elem(109)*2
C
C      data elem/'H ','He','Li','Be','B ',  'C ','N ','O ','F ','Ne',
C     &          'Na','Mg','Al','Si','P ',  'S ','Cl','Ar','K ','Ca',
C     &          'Sc','Ti','V ','Cr','Mn',  'Fe','Co','Ni','Cu','Zn',
C     &          'Ga','Ge','As','Se','Br',  'Kr','Rb','Sr','Y ','Zr',
C     &          'Nb','Mo','Tc','Ru','Rh',  'Pd','Ag','Cd','In','Sn',
C     &          'Sb','Te','I ','Xe','Cs',  'Ba','La','Ce','Pr','Nd',
C     &          'Pm','Sm','Eu','Gd','Tb',  'Dy','Ho','Er','Tm','Yb',
C     &          'Lu','Hf','Ta','W ','Re',  'Os','Ir','Pt','Au','Hg',
C     &          'Tl','Pb','Bi','Po','At',  'Rn','Fr','Ra','Ac','Th',
C     &          'Pa','U ','Np','Pu','Am',  'Cm','Bk','Cf','Es','Fm',
C     &          'Md','No','Lr','Rf','Db',  'Sg','Bh','Hs','Mt'/
C
C      elemord=0
C      aa=a
C      ich=ichar(aa(1:1))
C      if (ich.ge.97.and.ich.le.122) aa(1:1)=char(ich-32)
C      ich=ichar(aa(2:2))
C      if (ich.ge.65.and.ich.le.90) aa(2:2)=char(ich+32)
C
C      do i=1,109
C         if (aa==elem(i)) then
C            elemord=i; return
C         endif
C      enddo
C
C      end
C
Cc     ##############################################################
Cc     ##  subroutine NJ_blank -- find the n-th blank, LN=0 if no  ##
Cc     ##  2005.04.12 by Wei Li; Update 2005.10.16 by Wei Li       ##
Cc     ##############################################################
Cc
C****** 2005.4.12 if no LN=0
C      subroutine  NJ_blank(Inp,n,LN)
C      character line*100
C      external Nchar
C
C      rewind(Inp); Nct=0; LN=0
C
C      do 100 k=1,100000
C         read(Inp,'(a)',err=110,end=110) line
C         LN=LN+1
C         i=Nchar(line,1,100)
C         if (i.eq.0) Nct=Nct+1
C         if (Nct.eq.n) return
C 100  enddo
C
C 110  LN=0
C      end
C
CC     2004.04.16 find number of no-blank char in line(ini:ifi)(length<=100) #liwei@@itcc
C      function Nchar(line,ini,ifi)
C      character line*(*)
C      Nchar=0
C      do i=ini,ifi
C         if (line(i:i).ne.' ') Nchar=Nchar+1
C      enddo
C
C      end
C
Cc     ##############################################################
Cc     ##  subroutine NJ_seekkey  --  Seek key line with ch first  ##
Cc     ##  2004.12.24 by Wei Li; Update 2005.10.16 by Wei Li       ##
Cc     ##############################################################
Cc
C      subroutine NJ_seekkey(inp,kch,key,yes)
C      implicit none
C      integer inp,i,j,k,back,yes,ith,kth
C      character key*(*),line*256,kch*(*)
C
C      yes=0
C      back=0
C      ith=len(key)
C      kth=len(kch)
C      call NJ_upper(key)
C      call NJ_upper(kch)
C 200  do
C         read(inp,'(a)',err=800,end=800) line
C         call NJ_upper(line)
C         do j=1,256
C            if(line(j:j).ne.' ') exit
C         enddo
C         if (j>256) cycle
C         k=index(line(j:256),key(1:ith))
C         if (line(j:j+kth-1)==kch(1:kth).and.k.ne.0) then
C            backspace(inp); yes=1; return
C         endif
C      enddo
C
C 800  if (back==0) then
C         rewind(inp); back=1; goto 200
C      endif
C      end
C
CC     
C
CC === Divide a molecule into some groups which contain only two or three atoms ===
CC === Zhigang-2015.5.5 ===
CC === Before calling divi_mole, ngrp should be set as 0 ===
C      recursive subroutine divi_mole(numc,natc,NATOM,dis,ngrp,     
C     &                               natm_grp,atm_grp)
C      implicit none
C      integer::numc,NATOM,ngrp,i,j,max1,max2
C      integer::natc(numc),natm_sys(2),natm_grp(natom),atm_grp(30,natom)
C      real(kind=8)::dis(natom,natom),maxdis
C      integer,allocatable::atm_in_sys(:,:)
C      
C      maxdis=0.0D0
C      if (numc==1 .or. numc==2 .or. numc==3) then
C         ngrp=ngrp+1
C         natm_grp(ngrp)=numc
C         atm_grp(:numc,ngrp)=natc(:numc)
C         return
C      else
C         do i=1,numc-1
C            do j=i+1,numc
C      	       if (dis(natc(i),natc(j))>maxdis) then
C      	          maxdis=dis(natc(i),natc(j))
C      	          max1=natc(i)
C                  max2=natc(j)
C      	       end if
C      	    end do
C         end do
C         natm_sys=0
C         allocate(atm_in_sys(numc,2))
C         atm_in_sys=0
C         do i=1,numc
C            if (dis(natc(i),max1)<dis(natc(i),max2)) then
C      	       natm_sys(1)=natm_sys(1)+1
C      	       atm_in_sys(natm_sys(1),1)=natc(i)
C      	    else
C      	       natm_sys(2)=natm_sys(2)+1
C      	       atm_in_sys(natm_sys(2),2)=natc(i)
C      	    end if
C         end do
C         do i=1,2
C            call divi_mole(natm_sys(i),atm_in_sys(:,i),
C     &                     NATOM,dis,ngrp,natm_grp,atm_grp)
C         end do
C      	 deallocate(atm_in_sys)
C      end if
C      return
C      end subroutine divi_mole
C
CC === For larger smallest group and the smallest group contains 3 atoms - NZG-2015.5.13 ===
CC === Before calling divi_mole, ngrp should be set as 0.
C      recursive subroutine divi_mole3(numc,natc,NATOM,dis,ngrp,    
C     &                               natm_grp,atm_grp)
C      implicit none
C      integer::numc,NATOM,ngrp,i,j,max1,max2
C      integer::natc(numc),natm_sys(2),natm_grp(natom),atm_grp(30,natom)
C      real(kind=8)::dis(natom,natom),maxdis
C      integer,allocatable::atm_in_sys(:,:)
C
C      maxdis=0.0D0
C      if (numc==3 .or. numc==4 .or. numc==5) then
C         ngrp=ngrp+1
C         natm_grp(ngrp)=numc
C         atm_grp(:numc,ngrp)=natc(:numc)
C         return
C      else
C         do i=1,numc-1
C            do j=i+1,numc
C               if (dis(natc(i),natc(j))>maxdis) then
C                  maxdis=dis(natc(i),natc(j))
C                  max1=natc(i)
C                  max2=natc(j)
C               end if
C            end do
C         end do
C         natm_sys=0
C         allocate(atm_in_sys(numc,2))
C         atm_in_sys=0
C         do i=1,numc
C            if (dis(natc(i),max1)<dis(natc(i),max2)) then
C               natm_sys(1)=natm_sys(1)+1
C               atm_in_sys(natm_sys(1),1)=natc(i)
C            else
C               natm_sys(2)=natm_sys(2)+1
C               atm_in_sys(natm_sys(2),2)=natc(i)
C            end if
C         end do
C         if (natm_sys(1)<3 .or. natm_sys(2)<3) then
C            ngrp=ngrp+1
C            natm_grp(ngrp)=numc
C            atm_grp(:numc,ngrp)=natc(:numc)
C            return
C         else
C            do i=1,2
C               call divi_mole3(natm_sys(i),atm_in_sys(:,i),
C     &                     NATOM,dis,ngrp,natm_grp,atm_grp)
C            end do
C         end if
C         deallocate(atm_in_sys)
C      end if
C      return
C      end subroutine divi_mole3
C
C
CC
CC
Cc     ##############################################################
Cc     ##  subroutine NJ_compare -- compare if an array in another ##
Cc     ##  2005.04.04 by Wei Li; Update 2005.10.17 by Wei Li       ##
Cc     ##############################################################
Cc
Cc     if tmp1 in tmp2 ch='<'
Cc     if tmp2 in tmp1 ch='>'
Cc     if tmp1 eq tmp2 ch='='
Cc
C      subroutine NJ_compare(n,tmp1,tmp2,ch)
C      implicit none
C      integer n,tmp1(n),tmp2(n),i,j,k,L,m,n1,n2,ii,jj
C      character ch
C
C      ch='='
C      do i=1,n
C         if (tmp1(i).ne.tmp2(i)) then
C            ch='x'
C            exit
C         endif
C      enddo
CC
C      if (ch.eq.'=') return
CC
C      n1=0; n2=0
C      do i=1,n
C         if (tmp1(i).ne.0) n1=n1+1
C         if (tmp2(i).ne.0) n2=n2+1
C      enddo
CC
CC --- To check if tmp1 included in (<) tmp2
C      if (n1.lt.n2) then
C         k=0
C         do i=1,n1
C            ii=tmp1(i)
C            do j=1,n2
C               jj=tmp2(j)
C               if (ii.eq.jj) then
C                  k=k+1
C                  exit   
C               endif
C            enddo
C         enddo
C         if (k.eq.n1) then
C            ch='<'
C            return
C         endif
CC --- To check if tmp1 included (>) tmp2
C      elseif (n1.gt.n2) then
C         k=0
C         do j=1,n2
C            jj=tmp2(j)
C            do i=1,n1
C               ii=tmp1(i)
C               if (ii.eq.jj) then
C                  k=k+1
C                  exit   
C               endif
C            enddo
C         enddo
C         if (k.eq.n2) then
C            ch='>'
C            return
C         endif
C      else
C         return
C      endif
C
C      end
C
CC
CC
CC
CC --- 2008.02.28 Add for sorting LMOs according to the labels of atoms
C      subroutine SortLMO(io,nat,nbs,no,nob,ncor,nmo,smo,SR,SR1,SNB,SOB,
C     &           nuchar)
C      implicit none
C      integer io,nat,nbs,no,ncor,nmo,i,j,k,L,k1,k2,k3,k4,i1,j1,i2,j2,nob
C      integer SNB(nmo),SOB(nmo,nat)
C      real(kind=8) smo(nbs,nmo),SR(nat,nmo),SR1(nat,nmo),P,nuchar(nat)
C      logical LSort
CC
C      if (io.gt.0) write(io,*)
C     &   'Sort the occupied LMOs according to the labels of atoms'
CC
C      do i=1,ncor-1
C         do j=i+1,ncor
C            i1=SOB(i,1)
C            j1=SOB(j,1)
C            if (i1.gt.j1) then
C               call swapmo(i,j,nat,nbs,no,nmo,smo,SR,SR1,SNB,SOB)
C            endif
C         enddo
C      enddo
CC
C      do i=ncor+1,nob-1
C         do j=i+1,nob
C            LSort=.false.
C            if (SNB(i).eq.1 .and. SNB(j).eq.1) then
C               i1=SOB(i,1)
C               j1=SOB(j,1)
C               if (i1.gt.j1) LSort=.true.
C            elseif (SNB(i).eq.1 .and. SNB(j).gt.1) then
C               i1=SOB(i,1)
C               j1=SOB(j,1)
C               j2=SOB(j,2)
C               if (j1.gt.j2) call swap(j1,j2)
CC
C               k1=nint(nuchar(j1))
C               k2=nint(nuchar(j2))
C               if (k1.eq.1.and.k2.ne.1) call swap(j1,j2)
CC
C               if (i1.gt.j1) LSort=.true.
C            elseif (SNB(i).gt.1 .and. SNB(j).eq.1) then
C               i1=SOB(i,1)
C               i2=SOB(i,2)
C               j1=SOB(j,1)
C               if (i1.gt.i2) call swap(i1,i2)
CC
C               k1=nint(nuchar(i1))
C               k2=nint(nuchar(i2))
C               if (k1.eq.1.and.k2.ne.1) call swap(i1,i2)
CC
C               if (i1.ge.j1) LSort=.true.
C            else
C               i1=SOB(i,1)
C               i2=SOB(i,2)
C               j1=SOB(j,1)
C               j2=SOB(j,2)
C               if (i1.gt.i2) call swap(i1,i2)
C               if (j1.gt.j2) call swap(j1,j2)
CC
C               k1=nint(nuchar(i1))
C               k2=nint(nuchar(i2))
C               if (k1.eq.1.and.k2.ne.1) call swap(i1,i2)
C               k1=nint(nuchar(j1))
C               k2=nint(nuchar(j2))
C               if (k1.eq.1.and.k2.ne.1) call swap(j1,j2)
CC
C               if (i1.gt.j1) LSort=.true.
CC-WL- 2010.02.10 comment the line and replace it with the following lines
CC              if (i1.eq.j1 .and. i2.gt.j2) LSort=.true.
C               if (i1.eq.j1) then
C                  k1=nint(nuchar(i2))
C                  k2=nint(nuchar(j2))
C                  if (k1.eq.1.and.k2.eq.1) then
C                     if (i2.gt.j2) LSort=.true.
C                  elseif (k1.ne.1.and.k2.ne.1) then
C                     if (i2.gt.j2) LSort=.true.
C                  elseif (k1.ne.1.and.k2.eq.1) then
C                     LSort=.true.
C                  endif
C               endif
C            endif
C            if(LSort) call swapmo(i,j,nat,nbs,no,nmo,smo,SR,SR1,SNB,SOB)
C         enddo
C      enddo
CC
C      end
CC
CC --- Swap i and j
C      subroutine swap(i,j)
C      implicit none
C      integer i,j,L
C      L=i
C      i=j
C      j=L
C      end
CC
CC --- Swap MO i and j and other info. ---
C      subroutine swapmo(i,j,nat,nbs,no,nmo,smo,SR,SR1,SNB,SOB)
C      implicit none
C      integer nat,nbs,no,nmo,i,j,k,L
C      integer SNB(nmo),SOB(nmo,nat)
C      real(kind=8) smo(nbs,nmo),SR(nat,nmo),SR1(nat,nmo),P
CC
C      L=SNB(i)
C      SNB(i)=SNB(j)
C      SNB(j)=L
C      do k=1,nat
C         L=SOB(i,k)
C         SOB(i,k)=SOB(j,k)
C         SOB(j,k)=L
C         P=SR(k,i)
C         SR(k,i)=SR(k,j)
C         SR(k,j)=P
C         P=SR1(k,i)
C         SR1(k,i)=SR1(k,j)
C         SR1(k,j)=P
C      enddo
C      do k=1,nbs
C         P=smo(k,i)
C         smo(k,i)=smo(k,j)
C         smo(k,j)=P
C      enddo
CC
C      end
CC
CC---------------------------------------------------------------------
CC
CC --- 2008.03.03 Decide the symmetry of occupied orbitals ---
C      subroutine SymmOrb2(io,nmo,ncor,no,FIJ,IWORK,RWORK,KSymm)
C      implicit none
C      integer io,nmo,ncor,no,i,j,k,L,m,n,LLL,kk,k1,k2,k3,k4,k5,k6
C      integer KSymm(no),IWORK(100)  ! ,IrMO(no)
C      real(kind=8) FIJ(no,no),eps,pii,pjj,pij,p1,p2,RWORK(100)
C      character(len=2000) line,line1,line2,line3
C      integer,allocatable::symo(:,:),nsymo(:),ntmp(:),IrMO(:)
CC
C      allocate(IrMO(no))
CC
C      do i=1,no
C         KSymm(i)=i
C         IrMO(i)=0
C      enddo
CC
C      do i=2,no
C         pii=FIJ(i,i)
C         do j=1,i-1
C            pjj=FIJ(j,j)
C            if (dabs(pjj-pii).lt.RWORK(9)) then
C               KSymm(i)=KSymm(j)
C               exit
C            endif
C         enddo
C      enddo
CC
C      do i=1,no
C         if (KSymm(i).eq.i) KSymm(i)=0
C      enddo
CC
C      LLL=0
C      do i=ncor+1,no
C         if (KSymm(i).eq.0) then
C            LLL=LLL+1
C            IrMO(LLL)=i
C         endif
C      enddo
CC
C      if (LLL.eq.no-ncor) then
CC        write(io,*) 'No symmetric LMOs for this system'
C         return
C      endif
CC
C      allocate(symo(no,LLL),nsymo(LLL))
C      symo=0
C      nsymo=0
C      do i=1,LLL
C         symo(1,i)=IrMO(i)
C         m=IrMO(i)
C         L=1
C         do j=ncor+1,no
C            if (j.eq.m) cycle
C            if (KSymm(j).eq.m) then
C               L=L+1
C               symo(L,i)=j
C            endif
C         enddo
C         nsymo(i)=L
C      enddo
CC
C      write(io,*) '=== Checking the symmetry of occupied LMOs ==='
C      write(io,100) RWORK(9)
C      write(io,*) 'NO LMOs'
C      do i=1,LLL
C         write(line1,*) i
C         call NJ_trim(line1,k1,k2)
C         kk=nsymo(i)
C         call NJ_prtlab(line2,kk,symo(1,i))
C         call NJ_trim(line2,k3,k4)
C         write(io,*) line1(k1:k2)//' ('//line2(k3:k4)//')'
C      enddo
CC
C      if (IWORK(6).ne.0) then
C         IWORK(6)=LLL
C         write(io,*) 'The symmetry of occupied LMOs are used'
C      else
C         do i=1,no
C            KSymm(i)=0
C         enddo
C         write(io,*) 'The symmetry of occupied LMOs are not used'
C      endif
C      write(io,*)
CC
C 100  format(1x,'The threshold of Fock matrix =',d12.4)
C      deallocate(symo,nsymo,IrMO)
CC
C      end      
CC
CC --- Sort integer array ---
C      subroutine isort(n,A)
C      integer n,A(n),i,j,k1,k2,L
CC
C      do i=1,n-1
C         do j=i+1,n
C            if (A(i).gt.A(j)) then
C               L=A(i)
C               A(i)=A(j)
C               A(j)=L
C            endif
C         enddo
C      enddo
CC
C      end
CC
CC --- Sort integer array ---
C      subroutine absisort(n,A)
C      integer n,A(n),i,j,k1,k2,L
CC
C      do i=1,n-1
C         do j=i+1,n
C            if (abs(A(i)).gt.abs(A(j))) then
C               L=A(i)
C               A(i)=A(j)
C               A(j)=L
C            endif
C         enddo
C      enddo
CC
C      end
C
CC  ******************************************************************
CC  Construct the initial AO domain by combining all the central atoms
CC  of each occ LMO of the MO domain.  --Zhigang 2015.04.20
CC  ******************************************************************      
C      subroutine Init_AO(NATOM,NW,nb,SYS_ATOM,ISY_ATOM,J0,BA,JF,ZA)
C      implicit none
C      integer NATOM,NW,J0,JF,ISY_ATOM,J,K1,K2,K
C      integer nb(NATOM),BA(NATOM),ZA(NW),SYS_ATOM(NATOM)
C
C      J0=ISY_ATOM
C      BA(:J0)=SYS_ATOM(:J0)
C      ZA=0;JF=0
C      DO J=1,J0
C         K1=nb(BA(J))
C         if (BA(j).ne.NATOM) then
C            K2=nb(BA(J)+1)-1
C         else
C            K2=NW
C         endif
C         DO K=K1,K2
C            JF=JF+1
C            ZA(JF)=K
C         end do
C      end do
C      
C      end subroutine Init_AO
CC
CC---------------------------------------------------------------------
CC     CONSTRUCT THE ATOMIC DOMAIN OF EACH MO DOMAIN --07.01.26
C!        call COA(NATOM,NW,NMO,NUW,NB,SOB,SNB,INF(1,I),ISUB(I),
C!    &        NAsub(KGG),BA(1,KGG),NWsub(KGG),ZA(1,KGG))
CC---------------------------------------------------------------------
C      subroutine COA(NATOM,NW,NMO,NUW,nb,SOB,SNB,INF,ISUB,J0,BA,JF,ZA)
C      integer SOB(NMO,NATOM),SNB(NMO),INF(NUW),BA(NATOM),ZA(NW)
C      integer nb(NATOM)
CC
C      ii=INF(1)
C      J0=SNB(ii)
CC
C      DO J=1,J0
C         BA(J)=SOB(ii,J)
C      enddo
CC
C      J2=J0
C      DO J=2,ISUB
C         J1=INF(J)
C         DO 3689 K=1,SNB(J1)
C            DO  L=1,J0
C               IF (SOB(J1,K).EQ.BA(L))GO TO 3689
C            enddo
C            J2=J2+1
C            BA(J2)=SOB(J1,K)
C3689     CONTINUE
C         J0=J2
C      enddo
CC
C        ZA=0;JF=0
C        DO 3695 J=1,J0
C           K1=nb(BA(J))
C           if(BA(j).ne.NATOM)then
C            K2=nb(BA(J)+1)-1
C           else
C            K2=NW
C           endif
C           DO 3697 K=K1,K2
C              JF=JF+1
C              ZA(JF)=K
C3697       CONTINUE
C3695    CONTINUE
C      return
C      end
CC
C
CC
C
CC
CC
C
C      subroutine iunit(n,a)
C      integer a(n,n)
C      do i=1,n
C         do j=1,n
C            if (j.eq.i) then
C               a(j,i)=1
C            else
C               a(j,i)=0
C            endif
C         enddo
C      enddo
C      end
C
CC
CC
C
CC********************************************************************
CC CSSTQ method to diagonalize
CC********************************************************************
C        SUBROUTINE CSSTQ_0(N,B,C,Q,EPS,L)
C        DIMENSION B(N),C(N),Q(N,N)
C        DOUBLE PRECISION B,C,Q,D,H,P,R,F,E,S,G,eps
C        C(N)=0.0d0
C        D=0.0d0
C        F=0.0d0
C        EPS=0.0d0
C        DO 50 J=1,N
C          IT=0
C          H=EPS*(ABS(B(J))+ABS(C(J)))
C          IF (H.GT.D) D=H
C          M=J-1
C10        M=M+1
C          IF (M.LE.N) THEN
C            IF (ABS(C(M)).GT.D) GOTO 10
C          END IF
C          IF (M.NE.J) THEN
C15          IF (IT.EQ.60) THEN
C              L=0
C              WRITE(*,18)
C18            FORMAT(1X,'  FAIL')
C              RETURN
C            END IF
C            IT=IT+1
C            G=B(J)
C            P=(B(J+1)-G)/(2.0*C(J))
C            R=SQRT(P*P+1.0)
C            IF (P.GE.0.0) THEN
C              B(J)=C(J)/(P+R)
C            ELSE
C              B(J)=C(J)/(P-R)
C            END IF
C            H=G-B(J)
C            DO 20 I=J+1,N
C20          B(I)=B(I)-H
C            F=F+H
C            P=B(M)
C            E=1.0
C            S=0.0
C            DO 40 I=M-1,J,-1
C              G=E*C(I)
C              H=E*P
C              IF (ABS(P).GE.ABS(C(I))) THEN
C                E=C(I)/P
C                R=SQRT(E*E+1.0)
C                C(I+1)=S*P*R
C                S=E/R
C                E=1.0/R
C              ELSE
C                E=P/C(I)
C                R=SQRT(E*E+1.0)
C                C(I+1)=S*C(I)*R
C                S=1.0/R
C                E=E/R
C              END IF
C              P=E*B(I)-S*G
C              B(I+1)=H+S*(E*G+S*B(I))
C              DO 30 K=1,N
C                H=Q(K,I+1)
C                Q(K,I+1)=S*Q(K,I)+E*H
C                Q(K,I)=E*Q(K,I)-S*H
C30            CONTINUE
C40          CONTINUE
C            C(J)=S*P
C            B(J)=E*P
C            IF (ABS(C(J)).GT.D) GOTO 15
C          END IF
C          B(J)=B(J)+F
C50      CONTINUE
C        DO 80 I=1,N
C          K=I
C          P=B(I)
C          IF (I+1.LE.N) THEN
C            J=I
C60          J=J+1
C            IF (J.LE.N) THEN
C              IF (B(J).LE.P) THEN
C                K=J
C                P=B(J)
C                GOTO 60
C              END IF
C            END IF
C          END IF
C          IF (K.NE.I) THEN
C            B(K)=B(I)
C            B(I)=P
C            DO 70 J=1,N
C              P=Q(J,I)
C              Q(J,I)=Q(J,K)
C              Q(J,K)=P
C70          CONTINUE
C          END IF
C80      CONTINUE
C        L=1
C        RETURN
C        END
C
C        SUBROUTINE CSTRQ_0(A,N,Q,B,C)
C        DIMENSION A(N,N),Q(N,N),B(N),C(N)
C        DOUBLE PRECISION A,Q,B,C,F,H,G,H2
C        DO 10 I=1,N
C        DO 10 J=1,N
C10      Q(I,J)=A(I,J)
C        DO 80 I=N,2,-1
C          H=0.0
C          IF (I.GT.2) THEN
C            DO 20 K=1,I-1
C20          H=H+Q(I,K)*Q(I,K)
C          END IF
C          IF (H+1.0.EQ.1.0) THEN
C            C(I)=0.0d0
C            IF (I.EQ.2) C(I)=Q(I,I-1)
C            B(I)=0.0d0
C          ELSE
C            C(I)=SQRT(H)
C            IF (Q(I,I-1).GT.0.0d0) C(I)=-C(I)
C            H=H-Q(I,I-1)*C(I)
C            Q(I,I-1)=Q(I,I-1)-C(I)
C            F=0.0d0
C            DO 50 J=1,I-1
C              Q(J,I)=Q(I,J)/H
C              G=0.0d0
C              DO 30 K=1,J
C30            G=G+Q(J,K)*Q(I,K)
C              IF (J+1.LE.I-1) THEN
C                DO 40 K=J+1,I-1
C40              G=G+Q(K,J)*Q(I,K)
C              END IF
C              C(J)=G/H
C              F=F+G*Q(J,I)
C50          CONTINUE
C            H2=F/(H+H)
C            DO 70 J=1,I-1
C              F=Q(I,J)
C              G=C(J)-H2*F
C              C(J)=G
C              DO 60 K=1,J
C60            Q(J,K)=Q(J,K)-F*C(K)-G*Q(I,K)
C70          CONTINUE
C            B(I)=H
C          END IF
C80      CONTINUE
C        DO 85 I=1,N-1
C85      C(I)=C(I+1)
C        C(N)=0.0d0
C        B(1)=0.0d0
C        DO 130 I=1,N
C          IF ((B(I).NE.0.0).AND.(I-1.GE.1)) THEN
C            DO 110 J=1,I-1
C              G=0.0d0
C              DO 90 K=1,I-1
C90            G=G+Q(I,K)*Q(K,J)
C              DO 100 K=1,I-1
C100           Q(K,J)=Q(K,J)-G*Q(K,I)
C110         CONTINUE
C          END IF
C          B(I)=Q(I,I)
C          Q(I,I)=1.0d0
C          IF (I-1.GE.1) THEN
C            DO 120 J=1,I-1
C              Q(I,J)=0.0d0
C              Q(J,I)=0.0d0
C120         CONTINUE
C          END IF
C130     CONTINUE
C        RETURN
C        END
C
CC********************************************************************
CC Sort
CC********************************************************************
C        subroutine sort2(n,ra,rb)
C        real(kind=8) ra(n),rra
C        integer rb(n),rrb
C        l=n/2+1;ir=n
C        do
C         if(l.gt.1)then
C          l=l-1;rra=ra(l);rrb=rb(l)
C         else
C          rra=ra(ir);rrb=rb(ir);ra(ir)=ra(1);rb(ir)=rb(1);ir=ir-1
C          if(ir.eq.1)then
C           ra(1)=rra;rb(1)=rrb
C           return
C          endif
C         endif
C         i=l;j=l+l
C         do while(j.le.ir)
C          if(j.lt.ir)then
C           if(ra(j).lt.ra(j+1))j=j+1
C          endif
C          if(rra.lt.ra(j))then
C           ra(i)=ra(j);rb(i)=rb(j);i=j;j=j+j
C          else
C           j=ir+1
C          endif
C         enddo
C         ra(i)=rra;rb(i)=rrb
C        enddo
C        return
C        end
CC
CC
CC     2001.4.22. correct.
CC     Purpose: solve a system of linear equation by using the method of
CC     Gaussian elimination with partial pivoting.
CC
C      SUBROUTINE GAUSS2(io,N,A)
CC
C      real(kind=8) A(N,N+1),AV,T
CC
C      DO 30 K=1,N-1
C         AV=0.0d0
C         DO 10 I=K,N
C            IF (DABS(A(I,K)).LE.DABS(AV)) GOTO 10
C            AV=A(I,K)
C            L=I
C10       CONTINUE
C         IF (DABS(AV).LT.1.0D-8) THEN
C             WRITE(io,*) 'SINGULAR COEFFICIENT MATRIX'
C             STOP
C         ENDIF
C         IF (L.NE.K) THEN
C            DO 15 J=K,N+1
C               T=A(K,J)
C               A(K,J)=A(L,J)
C               A(L,J)=T
C15          CONTINUE
C         ENDIF
C         AV=1.0d0/AV
C         DO 25 J=K+1,N+1
C            A(K,J)=A(K,J)*AV
C            DO 20 I=K+1,N
C               A(I,J)=A(I,J)-A(I,K)*A(K,J)
C20          CONTINUE
C25       CONTINUE
CC
C30    CONTINUE
C      A(N,N+1)=A(N,N+1)/A(N,N)
C      DO 40 K=1,N-1
C         I=N-K
C         AV=0.0d0
C         DO 35 J=I+1,N
C            AV=AV+A(I,J)*A(J,N+1)
C35       CONTINUE
C         A(I,N+1)=A(I,N+1)-AV
C40    CONTINUE
C      RETURN
C      END
CC
CC
CC --- 2005.11.16 --- Transformation matrix for overlap by Symmetric Orthogonalization  ---
C      subroutine NJ_symorth(io,nbs,Evalu,Evect,Tmatr)
C      implicit none
C      integer io,nbs,i,j,k
C      real*8 Evalu(nbs),Evect(nbs,nbs),Tmatr(nbs,nbs),PP
C
C      Tmatr=0.0d0
C      do j=1,nbs
C         do i=1,nbs
C            do k=1,nbs
C               PP=Evect(i,k)*Evect(j,k)/dsqrt(Evalu(k))
C               Tmatr(i,j)=Tmatr(i,j)+PP
C            enddo
C         enddo
C      enddo
C
C      if (io>0) then
C         write(io,*) 'Transformation Matrix by Symmetric',
C     &               ' Orthogonalization'
C         call NJ_prtsym(io,nbs,Tmatr,'d14.6')
C         write(io,*)
C      endif
C
C      end
CC
CC
CC --- 2008.02.13 Add for the orthogonalization of MOs
C      subroutine OrthMO(io,nbs,nmo,mo,s,omo,norm,itype)  ! liweiii
C      implicit none
C      integer io,nbs,nmo,i,j,k,L,m,n,k1,k2,k3,k4,j1,j2,ii,mu,nu,norm
C      integer itype
C      real(kind=8) mo(nbs,nmo),s(nbs,nbs),P1,PP,P2,omo(nbs,nmo)
C      real(kind=8),allocatable::ss(:,:),q(:,:),e(:),vc(:),t(:,:)
CC
C      allocate(ss(nmo,nmo),q(nmo,nmo),e(nmo),vc(nmo),t(nmo,nmo))
CC
CC     write(io,*) 'Occupied LMOs in expanded subsystems'
CC     call NJ_prtcol2(io,nbs,nmo,mo,1,nmo,'f11.6')
CC
C      call NJ_tfock(0,nbs,nmo,s,mo,ss) ! s'=mo^T*s*mo
CC
CC     write(io,*) 'Occupied LMOs in expanded subsystems'
CC     do j1=1,nmo-1
CC        do j2=j1+1,nmo
CC           P1=ss(j2,j1)
CC           if (dabs(P1)>1d-6) write(io,'(2i4,f12.8)') j1,j2,P1
CC        enddo
CC     enddo
CC
C      call NJ_qr(io,ss,nmo,q,e,vc,k,0)
C      write(io,*) 'Eigenvalues:'
C      write(io,'(5f11.6)') e
C      if (itype.eq.1) then
C         write(io,*) 'Symmetric Orthogonalization is used'
C         call NJ_symorth(0,nmo,e,q,t)
C      else if (itype.eq.2) then
C         write(io,*) 'Canonical Orthogonalization is used'
C         call NJ_canorth(0,nmo,nmo,e,q,t)
C      else
C         write(io,*) 'Wrong type for orthogonalization'
C         write(io,*) 'itype in OrthMO() should be 1 or 2'
C         write(io,*) '1. Symmetric Orthogonalization'
C         write(io,*) '2. Canonical Orthogonalization'
C      endif
C      call NJ_matpro(io,nbs,nmo,nmo,mo,t,omo)
C      if (norm.ne.0) call normorb(io,nbs,nmo,omo,s)
CC
CC     write(io,*) 'Orthogonal occupied LMOs in expanded subsystems'
CC     call NJ_prtcol2(io,nbs,nmo,omo,1,nmo,'f11.6')
C      write(io,*) 'C^+SC='
C      call NJ_tfock(0,nbs,nmo,s,omo,ss) ! s'=mo^T*s*mo
CC
C      write(io,*) 'New --> Old     Overlap (> 0.1)'
C      do j=1,nmo
C         do i=1,nmo
C            P1=0d0
C            do mu=1,nbs
C               do nu=1,nbs
C                  P1=P1+mo(mu,i)*omo(nu,j)*s(nu,mu)
C               enddo
C            enddo
C            if (dabs(P1).gt.0.1d0) 
C     &         write(io,'(i4,'' -->'',i4,5x,''S='',f12.8)') j,i,P1
C         enddo
C      enddo
CC
C      deallocate(ss,q,e,vc,t)
C      end
CC
CC
CC --- Project nmo2 of nmo occ. LMOs from total system onto extended subsystems
C      subroutine projorb2(io,nbs,nbs2,nmo,nmo2,mo,mo2,s)
C      integer io,nbs,nbs2,nmo,i,j,k,L
C      real(kind=8) mo(nbs,nmo),mo2(nbs2,nmo2),s(nbs,nbs)
CC
C      if (nmo.eq.nmo2) then
C         write(io,100) nmo,nbs,nbs2
C      else
C         write(io,105) nmo,nmo2,nbs,nbs2
C      endif
CC
C      if (nbs2.gt.nbs) then
C         write(io,110) nbs2,nbs
C         call abrt
C      end if
CC
C      if (nmo2.gt.nmo) then
C         write(io,120) nmo2,nmo
C         call abrt
C      end if
CC
C      if (nbs2.eq.nbs) then
C         do j=1,nmo2
C            do i=1,nbs2
C               mo2(i,j)=mo(i,j)
C            enddo
C         enddo
C      else
C         do j=1,nmo2
C            CALL MINA2(io,nbs2,j,nbs,nmo,mo,s,mo2(1,j))
C         enddo
C      endif
C     
CC
C 100  format(3x,'[Project MOs]    MO=',i4,' AO:',i4,' -->',i4)
C 105  format(3x,'[Project MOs]    MO:',i4,' -->',i4,' AO:',i4,' -->',i4)
C 110  format(1x,'It is impossible to project orbitals on more AOs',2i5)
C 120  format(1x,'It is impossible to project more MOs than olds',2i5)
C      end
CC
CC
CC --- 2008.02.13 Add for the orthogonalization of MOs
C      subroutine orthorb(io,nbs,nmo,mo,s,itype)  ! liweiii
C      implicit none
C      integer io,nbs,nmo,i,j,k,L,m,n,k1,k2,k3,k4,j1,j2,ii,mu,nu,norm
C      integer itype,nzero
C      real(kind=8) mo(nbs,nmo),s(nbs,nbs),P1,PP,P2,zeps
C      real(kind=8),allocatable::ss(:,:),q(:,:),e(:),vc(:),t(:,:)
C      real(kind=8),allocatable::omo(:,:)
CC
C      zeps=1.0d-5
CC
C      allocate(ss(nmo,nmo),q(nmo,nmo),e(nmo),vc(nmo),t(nmo,nmo))
C      allocate(omo(nbs,nmo))
C      call NJ_tfock(0,nbs,nmo,s,mo,ss) ! s'=mo^T*s*mo  overlap over MO
C      call NJ_qr(0,ss,nmo,q,e,vc,k,0) ! QR diagonalization of ss
C      call numpoint(nmo,e,-zeps,zeps,nzero)  ! Count the number of values \in [-zeps,zeps] ("zero" eigenvalues)
C      if (nzero.ne.0) then
C         write(io,200) nzero
C         write(io,210) (e(i),i=1,nmo)
C      end if
CC
C      if (itype.eq.1) then
C         write(io,100) nmo,nbs
C         call NJ_symorth(0,nmo,e,q,t)
C      else if (itype.eq.2) then
C         write(io,110) nmo,nbs
C         call NJ_canorth(0,nmo,nmo,e,q,t)
C      else
C         write(io,*) 'Wrong type for orthogonalization'
C         write(io,*) 'itype in OrthMO() should be 1 or 2'
C         write(io,*) '1. Symmetric Orthogonalization'
C         write(io,*) '2. Canonical Orthogonalization'
C         stop
C      endif
C      call NJ_matpro(io,nbs,nmo,nmo,mo,t,omo)
CC
C      do i=1,nmo
C         do j=1,nbs
C            mo(j,i)=omo(j,i)
C         end do
C      end do
CC
C      deallocate(ss,q,e,vc,t,omo)
CC
C 100  format(3x,'[Symmetric orth] MO=',i4,' AO=',i4)
C 110  format(3x,'[Canonical orth] MO=',i4,' AO=',i4)
C 200  format(1x,'Warning: number of zero eigenvalues=',i4)
C 210  format(1x,5f14.6)
C      end
CC
CC --- Overlap between two set of MOs mo1(nbs,nmo1),mo2(nbs,nmo2)
C      subroutine locindx(io,nbs1,mo1,nbs2,mo2,nmo,s,indx,eta,nstar)
C      implicit none
C      integer io,nbs1,nbs2,nmo,i,j,k,L,i1,i2,j1,j2,k1,k2,nstar
C      real(kind=8) mo1(nbs1,nmo),mo2(nbs2,nmo),s(nbs2,nbs2)
C      real(kind=8) p1,pp,over,diff,pmax,pmin,indx(nmo),eta
C      character(len=1),allocatable::star(:)
C      character(len=1)::ps
C      parameter(diff=0.2d0)
C      integer imax,imin
CC
C      allocate(star(nmo))
CC
C      do i=1,nmo
C         star = ' '
C      enddo
CC
C      nstar=0
C      do i=1,nmo
C         p1=0d0
C         do j2=1,nbs2
C            do j1=1,nbs2
CCC             p1=p1+mo1(j1,i)*mo1(j2,i)*s(j1,j2)
CC              p1=p1+mo1(j1,i)*mo2(j2,i)*s(j1,j2) ! local index
C              p1=p1+mo2(j1,i)*mo2(j2,i)*s(j1,j2) ! overlap
C            end do
C         end do
CC         indx(i)=p1
C         indx(i)=1.0d0-p1
C**       indx(i)=dabs(p1-1.0d0)
C         if (abs(indx(i)).lt.eta) then
Ccc       if (indx(i).le.eta) then
CC        if (indx(i).ge.1d0-eta) then
C            star(i) = '*'
C            nstar = nstar + 1
C         endif
C      end do
CC
C      do i=1,nmo-1
C         do j=i+1,nmo
C            if (abs(indx(i)).le.abs(indx(j))) cycle
CC
C            pp=indx(i)
C            indx(i)=indx(j)
C            indx(j)=pp
CC
C            ps=star(i)
C            star(i)=star(j)
C            star(j)=ps
CC
C            do k=1,nbs2
C               PP=mo2(k,i)
C               mo2(k,i)=mo2(k,j)
C               mo2(k,j)=PP
C            enddo
C
C            do k=1,nbs1
C               PP=mo1(k,i)
C               mo1(k,i)=mo1(k,j)
C               mo1(k,j)=PP
C            enddo
C
C         enddo
C      enddo
CC
C      pmax=indx(1)
C      imax=1
C      pmin=indx(1)
C      imin=1
C      do i=2,nmo
C         if (indx(i).gt.pmax) then
C            pmax=indx(i)
C            imax=i
C         endif
C         if (indx(i).lt.pmin) then
C            pmin=indx(i)
C            imin=i
C         endif
C      enddo
CC
C      if (io.gt.0) then
C         write(io,100) nmo,nbs1,nbs2
C         write(io,110) (i,indx(i),star(i),i=1,nmo)
C         write(io,120) nstar,eta
C         write(io,130) pmax,imax,pmin,imin
C      end if
CC
C      deallocate(star)
C 100  format(1x,'local index of',i5,' MOs with',i5,' ->',i5,' AOs')
C 110  format(5(i5,'-',f7.4,a1))   
C 120  format(1x,'*',i5,' localization indices little than',f7.4)
C 130  format(1x,'Max:',f7.4,' in',i5,' and Min:',f7.4,' in',i5)
C      end
CC
CC
CC --- WL Sort jobs in .jobs file from the most expensive to cheapest one
C      subroutine sortjob(job,njobs,kord)
C      implicit none
C      integer i,j,k,L,k1,k2,k3,k4,job,njobs,kord(njobs)
C      character(len=256),dimension(:),allocatable::cjob
CC
C      allocate(cjob(4*njobs))
C      cjob=' '
CC
C      rewind(job)
C      read(job,*)
C      do j=1,4*njobs
C         read(job,'(a)') cjob(j)
C      enddo
CC
C      rewind(job)
C      read(job,*)
C      do i=1,njobs
C         k=kord(i)
C         k1=4*k-3
C         k2=4*k-2
C         k3=4*k-1
C         k4=4*k
C         write(job,'(a)') trim(cjob(k1))
C         write(job,'(a)') trim(cjob(k2))
C         write(job,'(a)') trim(cjob(k3))
C         write(job,'(a)') trim(cjob(k4))
C      enddo
CC
C      deallocate(cjob)
C      end
CC
CC
CC
C      subroutine NJ_dm(io,nbs,noc,cmo,dm)
C      implicit none
C      integer io,nbs,noc,i,j,k
C      real*8 cmo(nbs,noc),dm(nbs,nbs),PP
C      dm=0d0
C      do i=1,nbs; do j=1,nbs
C         PP=0d0
C         do k=1,noc
C            PP=PP+cmo(i,k)*cmo(j,k)
C         enddo
C         dm(i,j)=2d0*PP
C      enddo; enddo
C
C      if (io>0) then
C         write(io,*) '*** Density matrix ***'
C         call NJ_prtsym(io,nbs,dm,'d14.6')
C         write(io,*)
C      endif
C      return
C      end
CC
CC
CC
CC
C      subroutine ifind2(iunit, key)
C      implicit none
C      character key*(*),line*200
C      integer iunit,n,L,k,iyes,itry
CC
C      iyes=0
C      itry=0
C      L=len(key)
C 50   do
C         read(iunit,'(a)',end=100,err=100) line
C         if (index(line,key(1:L)).ne.0) then
C            iyes=1
C            backspace(iunit)
C            return
C         endif
C      end do
CC
C 100  itry=itry+1
C      if (itry.eq.1) then
C         rewind(iunit)
C         goto 50
C      endif
C
C      end
CC
CC
CC --- Updated on 25 AUG 10 for MP2 memory
CC     a new variable nbf (number of AOs) added
CC --- Calculating required memory in words
CC     RHF MP2,CCSD,CCSD(T),CR-CCL OR ROHF CCSD,CR-CCL
C      function memgms(scftyp,mplevl,cctyp,nbf,no,nu,nc)
C      IMPLICIT DOUBLE PRECISION(A-H,O-Z)
C      integer memgms,ityp,no,nu,nc,nbf
C      character scftyp*8,cctyp*8
C      logical DIRSCF
CC
C      if (nc.eq.0) then
C         ICIM=0
C      else
C         ICIM=2
C      endif
C      no2   = no*no
C      nou   = no*nu
C      nu2   = nu*nu
C      no3   = no2*no
C      no2u  = no2*nu
C      nou2  = no*nu2
C      nu3   = nu2*nu
C      no4   = no3*no
C      no3u  = no3*nu
C      no2u2 = no2*nu2
C      nou3  = no*nu3
C      nu4   = nu3*nu
CC
C      NBF2 = (NBF*NBF+NBF)/2
C      NBF3  = NBF*NBF
C      NOC   = no
C      NVIR  = nu
C      NINTMX= 15000
C      DIRSCF=.FALSE.
C      NINTIC= 0
C      NDIM  = MAX(NBF2,NVIR*NOC)
CC
C      if (scftyp.eq.'RHF     ') then
C         NU3s  = NU*(NU+1)*(NU+2)/6  !-CIM-
C      else if (scftyp.eq.'ROHF    ') then
C         NU3s  = NU*NU*(NU+1)/2  !-CIM-
C      end if
C      NCU3s = NC*NU3s  !-CIM-
C      KMICRO= 6
CC
C      LOADFM = 0
C      if (scftyp.eq.'RHF     ' .and. mplevl.eq.2) then  ! RHF MP2
C         IVEC   = 1    + LOADFM
C         IPTR1  = IVEC + NBF*NBF
C         IPTR2  = IPTR1+ NOC*NBF
C         IENG   = IPTR2+ NVIR*NBF
C         ILAB   = IENG + NBF
C         IIRP   = ILAB + NBF
C         IDEG   = IIRP + NBF
C         LAST   = IDEG + NBF
C         NEEDA  = LAST - LOADFM - 1
C      
C         MINTMX=NINTMX
C         IF(NINTIC.NE.0) MINTMX=0
C         IF(DIRSCF) THEN
C            LGHND  = 1      + LOADFM
C            LXINTS = LGHND  + MAXG
C            LDDIJ  = LXINTS + NSH2
C            LXX    = LDDIJ  + 49*MXG2
C            LIX    = LXX    + NINTMX
C            LAST   = LIX    + NINTMX
C         ELSE
C            LGHND  = 1      + LOADFM
C            LXINTS = LGHND
C            LDDIJ  = LXINTS
C            LXX    = LDDIJ
C            LIX    = LXX    + MINTMX
C            LAST   = LIX    + MINTMX
C            IF(NINTIC.NE.0) THEN
C               LXX  = LBUFPIC
C               LIX  = LIXIC
C            ENDIF
C         END IF
C         NEEDD = LAST - LOADFM - 1
C      
C         NEED  = 2*NBF3 + NVIR*NVIR*NOC*(NOC+1) + NOC*NOC*3
C         if (NEED.GT.800000000) then
C             NGOT = min(int(NEED*1.5),need+1000000000)
C         else
C             NGOT = 1000000000 ! 1000 MWORDS
C         endif
C
C         LEFT  = NGOT - NEED
C         NMIN  = NDIM*NBF
C         MNMEM = NEEDA + NEEDD + NEED + NMIN*1
C         MXMEM = NEEDA + NEEDD + NEED + NMIN*NOC
CC
C         LPASS = MIN(NOC,LEFT/NMIN)
CC        LPASS = MIN(NOC,1000000/NBF2+1)
CC
CC -------------------------
C         LOADFM=0
C         LWRK1 = 1     + LOADFM
C         LWRK2 = LWRK1 + NBF3
C         LPQRJ = LWRK2 + NBF3
C         IF (ICIM.EQ.2) THEN
C            LNOA  = LPQRJ + NDIM*NBF*LPASS  !-CIM- 25 FEB 10
C            LNOB  = LNOA  + NOC             !-CIM- 25 FEB 10
C            LPOA  = LNOB  + NOC             !-CIM-
C            LPOB  = LPOA  + NOC*NOC         !-CIM-
C            LTX   = LPOB  + NOC*NOC         !-CIM- 25 FEB 10
C            LYO2  = LTX   + NOC*NOC         !-CIM- 25 FEB 10
C            LYT2T = LYO2  + NOC             !-CIM- 25 FEB 10
C            LYT2S = LYT2T + NOC             !-CIM- 25 FEB 10
C            LTR1  = LYT2S + NOC             !-CIM- 25 FEB 10
C            LTR2  = LTR1  + NOC
C            LTC1  = LTR2  + NVIR*NVIR*NOC*(NOC+1)/2
C            LTC2  = LTC1  + NOC
C            LAST  = LTC2  + NVIR*NVIR*NOC*(NOC+1)/2
C         ELSE
C            LAST  = LPQRJ + NDIM*NBF*LPASS  !-CIM- 25 FEB 10
C         ENDIF
C         NEEDE = LAST - LOADFM - 1
C         MEM   = NEEDA + NEEDD + NEEDE      ! ~= NBF*NBF2*LPASS
CC
CC -------------------------
CC
C         NEED = MEM
Cc        NEED = 2.5*NVIR*NVIR*NOC*NOC+5*NDIM*NBF   ! NEW DEF. MEMORY
C         NEED = 5*NVIR*NVIR*NOC*NOC+5*NDIM*NBF   ! NEW DEF. MEMORY
C      elseif (scftyp.eq.'RHF     ' .and.
C     &       (cctyp.eq.'CCSD    '.or.cctyp.eq.'CCD     ')) then  ! RHF CCSD/CCD
C         LO1  = LOADFM + 1
C         LT1  = LO1    + NOU
C         LFH  = LT1    + NOU
C         LFPH = LFH    + NO2
C         LFP  = LFPH   + NOU
C         LVHH = LFP    + NU2
C         LVM  = LVHH   + NO4
C         LTI  = LVM    + NO3U
C         LO2  = LTI    + NU3
C         LT2  = LO2    + NO2U2
C         LVL  = LT2    + NO2U2
C         LVR  = LVL    + NO2U2
C         LVPP = LVR    + NO2U2
C         IF (ICIM.EQ.2) THEN  !-CIM-
C            LNO  = LVPP   + NOU3
C            LK   = LNO    + NO   !-CIM-   ECIM(NO)
C            LTX  = LK     + NO2  !-CIM-   EK(NO,NO)
C            LYO2 = LTX    + NO2  !-CIM-   TX(NO,NO)  2009.08.06
C            LYT2 = LYO2   + NO   !-CIM-   YO2(NO)    2009.08.06
C            LAST = LYT2   + NO   !-CIM-   YT2(NO)    2009.08.06
C         ELSE
C            LAST = LVPP   + NOU3
C         END IF
C         NEED = LAST - LOADFM - 1
C      else if (scftyp.eq.'RHF     ' .and. cctyp.eq.'CCSD(T) ') then  ! RHF CCSD(T)
C         I1   = LOADFM + 1
C         I2   = I1     + NOU
C         I3   = I2     + NO2U2
C         I4   = I3     + NO3U
C         I5   = I4     + NOU3
C         I6   = I5     + NU3
C         I7   = I6     + NU3
C         I8   = I7     + NO2U2
C         LAST = I8     + NOU
C         IF (ICIM.EQ.2) THEN  !-CIM-
C            I9   = I8     + NOU
C            I10  = I9     + NO      ! ECIM(NO)  -CIM-
C            I11  = I10    + NO2     ! TX(NO,NO) -CIM-
C            I12  = I11    + NU3*NC  ! XF3       -CIM-
C            I13  = I12    + NU3*NC  ! XT3       -CIM-
C            I14  = I13    + NU*NC   ! XO1       -CIM-
C            LAST = I14    + NU*NC   ! XT1       -CIM-
C         ELSE
C            LAST = I8     + NOU
C         END IF
C         NEED = LAST - LOADFM - 1
C      else if (scftyp.eq.'RHF     ' .and. cctyp.eq.'CR-CCL  ') then  ! RHF CR-CC(2,3)
CC        memgms = (no+nu) + (2*no2+2*nu2+6*nou)+3*nu3+2*nou2+2*no2u
CC        memgms =  memgms + 3*no3u+5*no2u2+3*nou3
CC        if (nc.ne.0) memgms = memgms + no+no2+11*nu3*nc
C         I1   = LOADFM + 1
C         I2   = I1     + NOU
C         I3   = I2     + NO2U2
C         I4   = I3     + NO3U
C         I5   = I4     + NOU3
C         I7   = I5     + NU3
C         I8   = I7     + NOU
C         I9   = I8     + NO2U2
C         I10  = I9     + NOU3
C         I11  = I10    + NO3U
C         I12  = I11    + NO2U2
C         I13  = I12    + NU3
C         I14  = I13    + NO
C         I15  = I14    + NU
C         I16  = I15    + 2*NU2
C         I17  = I16    + 2*NO2
C         I18  = I17    + 2*NOU
C         I19  = I18    + 2*NOU2
C         I20  = I19    + 2*NO2U
C         I21  = I20    + NOU
C         I22  = I21    + NO2U2
C         I23  = I22    + NU3
C         I24  = I23    + NOU3
C         I26  = I24    + NO3U
C         I27  = I26    + NOU
C         IF (ICIM.EQ.2) THEN
C            I29  = I27    + NO2U2    !-CIM-
C            I30  = I29    + NO*2     !-CIM- ECIM  09 MAR 10 NO --> NO*2
C            I31  = I30    + NO2      !-CIM- TX
C            I32  = I31    + NU3*NC   !-CIM- XV3
C            I33  = I32    + NCU3s*20 !-CIM- XBL
C            LAST = I33    + NC       !-CIM- TXs
C         ELSE
C            LAST = I27    + NO2U2
C         END IF
C         NEED = LAST - LOADFM - 1
C      else if (scftyp.eq.'ROHF    ' .and. cctyp.eq.'CCSD    ') then ! ROHF CCSD
CC        memgms = 1*no4+1*no3u+9*no2u2+3*nou3 + 4*nou+2*no2+2*nu2+nu3
C         I1   = LOADFM + 1
C         I2   = I1     + NOU       ! O1AA
C         I3   = I2     + NOU       ! T1
C         I4   = I3     + NOU       ! O1BB
C         I5   = I4     + NO4       ! VHHAA
C         I6   = I5     + NO2U2     ! VHPLAA
C         I7   = I6     + NO2U2     ! VHPLBB
C         I8   = I7     + NO2U2     ! O2AA
C         I9   = I8     + NO2U2     ! O2BB
C         I10  = I9     + NO2U2     ! O2AB
C         I11  = I10    + NOU3      ! VEAB
C         I12  = I11    + NOU3      ! VEBA
C         I13  = I12    + NO2U2     ! VHPRBB
C         I14  = I13    + NO2U2     ! VHPRAA
C         I15  = I14    + NO2U2     ! VHPRAB
C         I16  = I15    + NOU3      ! TI
C         I17  = I16    + NO3U      ! VMAA
C         I18  = I17    + NOU       ! FHP
C         I19  = I18    + NO2       ! FHHAA
C         I20  = I19    + NO2       ! FHHBB
C         I21  = I20    + NU2       ! FPPAA
C         I22  = I21    + NU2       ! FPPBB
C         I23  = I22    + NO2U2     ! T2
C         I24  = I23    + NU3       ! VPP
C         I25  = I24    + (KMICRO+1)**2  ! XMAT
C         I26  = I25    + (KMICRO+1)     ! BVEC
C         IF (ICIM.EQ.2) THEN
C            I27  = I26    + (KMICRO+1)     ! IPVT
C            I28  = I27    + NO       !-CIM- ECIMA
C            I29  = I28    + NO       !-CIM- ECIMB
C            I30  = I29    + NO2      !-CIM- TXA(NOA,NOA)
C            I31  = I30    + NO2      !-CIM- TXB(NOB,NOB)
C            I32  = I31    + NO2      !-CIM- EKA(NOA,NOA)
C            LAST = I32    + NO2      !-CIM- EKB(NOB,NOB)
C         ELSE
C            LAST = I26    + (KMICRO+1)     ! IPVT
C         END IF
C         NEED=LAST-LOADFM-1
C      else if (scftyp.eq.'ROHF    ' .and. cctyp.eq.'CR-CCL  ') then ! ROHF CR-CC(2,3)
CC        memgms = 6*no3u+10*no2u2+6*nou3 
CC        memgms = memgms + 8*nou+3*nu3+5*no2+5*nu2+4*nou2+4*no2u+no4
CC        if (nc.ne.0) memgms = memgms + no4+2*no+2*no2+2*nu3*nc
CC        NEED = memgms
C         I1  = LOADFM + 1
C         I2   = I1     + NOU       ! O1AA
C         I3   = I2     + NOU       ! O1BB
C         I4   = I3     + NO2U2     ! L2AA
C         I5   = I4     + NO2U2     ! L2BB
C         I6   = I5     + NO2U2     ! L2AB
C         I7   = I6     + NO2U2     ! O2
C         I8   = I7     + NO2U2     ! O2AA
C         I9   = I8     + NO2U2     ! O2BB
C         I10  = I9     + NO2U2     ! O2AB
C         I11  = I10    + NOU3      ! VEAA
C         I12  = I11    + NOU3      ! VEBB
C         I13  = I12    + NOU3      ! VEAB
C         I14  = I13    + NOU3      ! VEBA
C         I15  = I14    + NOU3      ! VEAB21
C         I16  = I15    + NOU3      ! VEBA21
C         I17  = I16    + NOU       ! FHPAA
C         I18  = I17    + NOU       ! FHPBB
C         I19  = I18    + NO2U2     ! VHPRAA
C         I20  = I19    + NO2U2     ! VHPRBB
C         I21  = I20    + NU3       ! M3
C         I22  = I21    + NU3       ! L3
C         I23  = I22    + NO2U2     ! VHPRAB
C         I24  = I23    + NU3       ! TI
C         I25  = I24    + NO3U      ! VMAB
C         I26  = I25    + NO3U      ! VMBA
C         I27  = I26    + NO3U      ! VMAA
C         I28  = I27    + NO3U      ! VMBB
C         I29  = I28    + NO3U      ! VMAB21
C         I30  = I29    + NO3U      ! VMBA21
C         I31  = I30    + NO2       ! FHHA
C         I32  = I31    + NO2       ! FHHB
C         I33  = I32    + NU2       ! FPPA
C         I34  = I33    + NU2       ! FPPB
C         I35  = I34    + NO2       ! X1AA
C         I36  = I35    + NO2       ! X1BB
C         I37  = I36    + NO2       ! X1AB
C         I38  = I37    + NU2       ! X2AA
C         I39  = I38    + NU2       ! X2BB
C         I40  = I39    + NU2       ! X2AB
C         I41  = I40    + NOU       ! X3AA
C         I42  = I41    + NOU       ! X3BB
C         I43  = I42    + NOU       ! X3AB
C         I44  = I43    + NOU       ! X3BA
C         I45  = I44    + NOU2      ! X4AAA
C         I46  = I45    + NOU2      ! X4AAB
C         I47  = I46    + NOU2      ! X4BBA
C         I48  = I47    + NOU2      ! X4BBB
C         I49  = I48    + NO2U      ! X5AAA
C         I50  = I49    + NO2U      ! X5AAB
C         I51  = I50    + NO2U      ! X5BBA
C         I52  = I51    + NO2U      ! X5BBB
C         IF (ICIM.EQ.2) THEN
C            I53  = I52    + NO4       !-CIM- VHHHH
C            I54  = I53    + NOA       !-CIM- ECIMA
C            I55  = I54    + NOB       !-CIM- ECIMB
C            I56  = I55    + NOA2      !-CIM- TXA
C            I57  = I56    + NOB2      !-CIM- TXB
C            I58  = I57    + NCU3s     !-CIM- XXM3
C            I59  = I58    + NCU3s     !-CIM- XXL3
C            LAST = I59    + NC        !-CIM- TXs
C         ELSE
C            LAST = I52    + NO4       ! VHHHH
C         END IF
C         NEED=LAST-LOADFM-1
C      end if
CC
C      memgms = NEED
C      return
CC
C      end
CC
CC
Cc
CC
C
CC === Use groups as the building units to construct primitive clusters  ===
CC === Zhigang-2015.5.9 ===
C      subroutine Frag_Deter(io,ngrp,dis,norb_grp,SYS,NS0,NS1,thre,kzhi)
C      implicit none
C      integer::io,ngrp
C      integer::SYS(ngrp,ngrp),NS0(ngrp),NS1(ngrp),norb_grp(ngrp)
C      real(kind=8)::dis(ngrp,ngrp),thre
C      integer::i,j,k,l,tmp(ngrp),KK,k1,k2,kzhi
C      integer,allocatable::Ktmp(:)
C      character::line*1000
C
C      SYS(:,:)=0;NS0(:)=1;NS1(:)=0
C      
C      kzhi=0
C      do i=1,ngrp
C         if (norb_grp(i)==0) cycle
C         kzhi=kzhi+1
C         SYS(1,kzhi)=i
C         NS1(kzhi)=1
C         do j=1,ngrp
C            if (j==i) cycle
C            if (norb_grp(j)==0) cycle
C            if (dis(j,i)<=thre) then
C               NS1(kzhi)=NS1(kzhi)+1
C               SYS(NS1(kzhi),kzhi)=j
C            end if
C         end do
CC
C         tmp(:)=0
C         tmp(2:NS1(kzhi))=SYS(2:NS1(kzhi),i)
CC
CC === For j from 2 to NS1(i), sort from small to large ===
C         do j=2,NS1(kzhi)
C            KK=tmp(2);l=2
C            do k=2,NS1(kzhi)
C               if (tmp(k)>KK) then
C                  KK=tmp(k);l=k
C               end if
C            end do
C            SYS(NS1(kzhi)-j+2,i)=KK;tmp(l)=0
C         end do
C      end do
C
C      if(io>0) then
C        write(io,'(1x,a)')'----- Original SYS by Group -----'
C        write(io,'(1x,a8,a6,a13)')'SYS','NUM','GROUP_INDEX'
C        do i=1,kzhi
C          allocate(Ktmp(NS1(i)))
C          Ktmp(1:NS1(i))=SYS(1:NS1(i),i)
C          call NJ_prtlab(line,NS1(i),Ktmp)
C          deallocate(Ktmp)
C          call NJ_trim(line,k1,k2)
C          write(io,'(1x,I8,I5,4x,a)')i,NS1(i),line(k1:k2)
C        end do
C        write(io,'(1x,a)')'---------------------------------'
C      end if
C
C      flush(io)
C      return
C      end subroutine Frag_Deter
C
C
C
CC ***** Include atoms to the AAO domain 2015.1.8 Zhigang *****
C      subroutine Inc_Atom(io,NATOM,link,dis,BA_AAO,J0,J1,thre)
C      implicit none
C      integer::io,NATOM,i,j,k,l,AA,BB,J0,J1,J00
C      integer::link(NATOM,NATOM),iindex(NATOM),BA_AAO(NATOM)
C      real(kind=8)::dis(NATOM,NATOM),thre
C
C      J00=J0
C      do i=J0+1,NATOM
C          do j=1,J0
C              if (dis(BA_AAO(i),BA_AAO(j))<thre) then
C	          write(io,*) "Add atom ",i
C		  J00=J00+1
C		  iindex(J00)=i
C		  exit
C	      end if
C	  end do
C      end do
C
C      do i=J0+1,J00
C          AA=BA_AAO(i)
C          BA_AAO(i)=BA_AAO(iindex(i))
C          BA_AAO(iindex(i))=AA
C      end do 
C      
C      J1=J00
C      write(io,'(I8,a27,f8.2,a)')
C     &             J00-J0,'ATOM ADDED TO AO DOMAIN IN',thre,' Angstrom'
C      write(io,'(10I8)')(BA_AAO(i),i=J0+1,J00)
C      write(io,'(I8,a27,f8.2,a)')
C     &             J00-J0,'ATOM ADDED TO AO DOMAIN IN',thre,' Angstrom'
C      write(io,*)''
C
C      return
C      end subroutine Inc_Atom
C
C
C
C      subroutine inarraych(element,array,num,in)
C      implicit none
C      character*2::array(num)
C      character*2::element
C      integer::num,i,in
C		      
C      do i=1,num
C	if (element==array(i)) then
C	  in=1
C	  return
C	else
C	  in=0  
C	end if
C      end do
C      return
C      end subroutine inarraych
C      
C
C      subroutine choldc2(a,n,np,da,eps)
C      integer n,np
C      real(kind=8) a(np,np),p(n),da(n,n),eps
C      integer i,j,k
C      real(kind=8),allocatable::sum(:)
CC
C      allocate(sum(n))
C      sum=0.0d0
C      do i=1,n
C         do j=i,n
C            sum(j)=a(i,j)
C         enddo
CC --- notes: if we exchange k and j cycle, the subroutine is much slower
CC C100H202: 6-31G** overlap: 2s --> 43s
C         do k=i-1,1,-1
C            if (dabs(da(i,k)).lt.eps) cycle
C            do j=i,n
C               sum(j)=sum(j)-da(i,k)*da(j,k)
C            enddo
C         enddo
C         if (sum(i).le.0) pause 'choldc failed'
C         da(i,i)=dsqrt(sum(i))
C         do j=i+1,n
C            da(j,i)=sum(j)/da(i,i)
C         enddo
C      enddo
C      deallocate(sum)
C      return
C      end
CC
C      subroutine choldc3(io,a,n,np,da,eps)
C      integer n,np
C      real(kind=8) a(np,np),p(n),da(n,n),PP,eps
C      integer i,j,k
C      real(kind=8),allocatable::sum(:)
C      real*8 Tim0,CPUTim
C      integer system,Wall0,Wall,TIME
C      external CPUTim,system
C      real(kind=8) t0,t1,t2,tt1,tt2
C      integer nmult,ndivd,ncomp
CC
C      nmult=0;ndivd=0
C      allocate(sum(n))
C      do i=1,n
C         sum=0.0D0
C         do k=1,i-1
C            if (dabs(da(i,k)).lt.eps) cycle
C            do j=i,n
C               sum(j)=sum(j)-da(i,k)*da(j,k)
C               nmult=nmult+1
C            enddo
C         enddo
C         sum(i)=sum(i)+a(i,i)
C         if (sum(i).le.0) pause 'choldc failed'
C         da(i,i)=dsqrt(sum(i))
C         if (dabs(da(i,i)).lt.eps) da(i,i)=0d0
C         do j=i+1,n
C            sum(j)=sum(j)+a(i,j)
C            da(j,i)=sum(j)/da(i,i)
C            ndivd=ndivd+1
C            if (dabs(da(j,i)).lt.eps) da(j,i)=0d0
C         enddo
C      enddo
C      deallocate(sum)
C 301  format(1x,'Num(*,/)=',2i12,  3x,'(Tot,Ratio)=', i12,  f7.2)
C 302  format(1x,'Tim(*,/)=',2f12.4,3x,'(Tot,Ratio)=', f12.4,f7.2)
C      return
C      end
C
CC
C 
C      subroutine CalcTwal(io,Twal0,cinit)
C      implicit real(kind=8) (A-H,O-Z)
C      integer TIME,Twal,Twal0
C      character cinit*(*)
C      Twal = TIME()-Twal0
C      NDays= Twal/(3600*24)
C      Twal = Twal-NDays*(3600*24)
C      NHours= Twal/3600
C      Twal = Twal-NHours*3600
C      NMin = Twal/60
C      Twal = Twal-NMin*60
C      write(io,1000) cinit, NDays, NHours, NMin, Twal
C      return
C 1000 format(1X,A,' Wall time:',I3,' days ',I2,' hours ',
C     &        I2,' minutes ',I4,' seconds.')
C      end subroutine CalcTwal
C
C

*************************************************
C *  generate Cluster-in-Molecule cluster input file  *
C *  Modify the code from GAMESS                      *
C *  NZG_4/16/2016 @UARK                              *
C *****************************************************

* -- LIST OF VARIABLES [integer] ---
* NATOM:    Number of atoms
* NBas:     Number of basis functions
* Nmo:      Number of all MOs (Number of independant functions, Nmo.le.Nbas)
* nocc:     Number of occupied MOs
* icha:     Number of charges
* mult:     Multiplicity
* nel:      Number of electrons
* k_alph:   Number of alpha electrons
* k_beta:   Number of beta electrons
* ncs:      Number of contracted shells
* nsh:      Number of primitive shells
* tbs(nbas):Atomic label for basis
* nfocc:    Number of frozen occ MOs (default: nfocc=0)
* kmem:     Max memory for integral transformation (unit: Words, e.g. 8 Bytes)
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
*     ED(NUW):           The correlation energy of each occupied LMO
*     NCLU(NUW):         The calculated times of each occupied LMO
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

      subroutine para_cimsubgen(nclu)

      use memory
      use newpara
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


      parameter (nopt=4)
      integer nsh,nprint,nbas,NATOM
      integer nmo,nalpha,nbeta,ierr,nfocc,np,jerr,debug,nocc,stat
      integer ioptyp(nopt)
      character*4 word(nopt)
      character*8 submtd
      character*7 methods(5)
      character*256 chopv(nopt),jobname,MOS,MOB,cluname,inpname
      character*256 infname,mosname
      character line*100,basline*100
      dimension iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      integer,allocatable::IAN(:),natc(:),nath(:),link(:,:),link2(:,:)
      integer,allocatable::atmclu(:,:,:),natmclu(:,:),Z(:),INX(:,:)
      integer,allocatable::ILST(:,:),NBAtm(:),batom(:,:),tbs(:)
      integer,allocatable::SR(:,:),SR1(:,:),SOB(:,:)
      integer,allocatable::atmlevl(:),molevl(:),nocc_atm(:),occ_atm(:,:)
      integer,allocatable::orb_clu(:,:),ML(:,:),norb_clu(:),ncenorb(:)
      integer,allocatable::highorlow(:),clumtd(:),atm_clu(:,:),nb_clu(:)
      integer,allocatable::natm_clu(:),bas_clu(:,:),nvir_clu(:)
      integer,allocatable::ZA(:,:),BA(:,:),cenMO(:,:)
      integer,allocatable::natm_clu_slv(:),JF_slv(:),KB_slv(:),KV_slv(:)
      real*8,allocatable::XC(:,:),coor(:,:),QA(:),dis(:,:),basdat(:,:)
      real*8,allocatable::SMO(:,:),SMOW(:,:),SMOT(:,:),eorb(:)
      real*8,allocatable::SMOAW(:,:),SMOAT(:,:),SMOA(:,:),norm(:,:)
      real*8,allocatable::SMOBW(:,:),SMOBT(:,:),SMOB(:,:),tmp(:)
      real*8,allocatable::MOCOOR(:,:),modis(:,:),PAO(:,:),MOS1(:,:)
      real*8,allocatable::MO_clu(:,:),MOS2(:,:),FK(:,:),FKW(:,:)
      real*8,allocatable::F2(:,:),FH(:,:),FHH(:),FIA(:,:),TRANS(:,:)
      real*8,allocatable::MAT(:,:),VECT(:,:),VALU(:),VC(:),dtmp(:)
      real*8,allocatable::virtime(:),virtime_slv(:)
      real*8,external::dtrace2
      character*8,allocatable::AtSymb(:)
c-----------------------------------------------------------------------
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c-----------------------------------------------------------------------
      data ioptyp  / 1,     11,    21,    11   /
      data word    /'prin','disl','subm','virt'/
      data IUnit   /1/ ! unit number for checkpoint I/O
      data methods /'RIMP2  ','MP2    ','CCD    ','CCSD   ','CCSD(T)'/
      common /job/ jobname,lenJ

 10   format(72('='))
      write(iout,10)  
      write(iout,"(20x,'Generate CIM cluster input files')")
      write(iout,*) 
      dislmo=5.0D0
      dislmo2=3.5D0
      submtd='mp2'
      virt=0.05D0
      nprint=0

      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
      if (ifound(1).gt.0) nprint=iopv(1,1)
      if (nprint>0) nprint=iout
      if (ifound(2).gt.0) dislmo=ropv(1,2)
      if (ifound(3).gt.0) submtd=chopv(3)
      if (ifound(4).gt.0) virt=ropv(1,4)
      call elapsec(tstart)
C
C  Read the Cartesian coordinates
C
      call getival('na',NAtom)
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
      call setival('nocc_cim',nocc)
      call getival('nonredun',nmo)
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
      iafm=mataddr('fock')
      do i=1,nbas
         do j=1,i
            FKW(j,i)=bl(iafm)
            FKW(i,j)=bl(iafm)
            iafm=iafm+1
         end do
      end do
      CALL ReorderFock(nbas,nbas,Z(i1),FKW,FK)

C      allocate(norm(nmo,nmo))
C      do i=1,nmo
C         do j=1,nmo
C            norm(i,j)=0.0D0
C            do k=1,nbas
C               do l=1,nbas
C                  norm(i,j)=norm(i,j)+SMO(k,i)*FK(k,l)*SMO(l,j)
C               end do
C            enddo
C         end do
C      end do
C      call prntmat(nmo,nmo,nmo,norm)    

      call FRZORB(nprint,NATOM,IAN,nfocc)
      nvalocc=nocc-nfocc  
      write(iout,"(' NUMBER OF BASIS FUNCTIONS:',i5)") NBAS
      write(iout,"(' NUMBER OF MOS            :',i5)") nmo
      write(iout,"(' NUMBER OF OCC MOS        :',i5)") NOCC
      write(iout,"(' NUMBER OF FROZEN OCC MOS :',i5)") nfocc
      write(iout,*)
C --- Calculate Mulliken Population and Assign the LMOs to a certain atom---
      allocate(SR(natom,nmo),SR1(natom,nmo),SOB(nmo,natom))
      allocate(atmlevl(natom),molevl(nocc),highorlow(nocc))

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
     
      write(iout,"(' MO DISTANCE THRESH:',f5.2)") dislmo
      write(iout,"(' MO DISTANCE THRESH (high-level):',f5.2)") dislmo2
      write(iout,*)
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
     &                 nclu,molevl)
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
      write(iout,*) ' Clu_INDX   NOcc   Natom   Nbas   Method'
      flush(iout)
      do i=1,nclu
         write(iout,'(i7,1x,3i8,5x,a8)') i,norb_clu(i),natm_clu(i),
     &                                   nb_clu(i),methods(clumtd(i))
      enddo
      write(iout,*)
      
      allocate(PAO(nbas,nbas))
      call dunit(nbas,PAO)
      call pjotorb(nbas,nocc,nbas,SMO(:,1:nocc),PAO,SOVER) 
      call normorb(nbas,nbas,PAO,SOVER)

      open(unit=iunit,file=jobname(1:lenJ)//'.inp')
      do while(.true.)
         read(unit=iunit,fmt="(A50)",iostat=stat) line
         call NJ_trim(line,L5,L6)
         call NJ_lower(line(L5:L6))
         if (line(L5:L5+4)=='basis') basline=line
         if (line(L5:L5+2)=='cim') exit
      end do
      close(iunit)
      
C --- Construct virtual space for each cluster!
C nvir_clu(i)    : Number of virtual orbitals of cluster i
C MOS1(nbas,nbas): MO coefficients of the cluster. 
C                  Shape: nb_clu(i)*(norb_clu(i)+nvir_clu(i))
C virt           : Threshold for projecting PAOs (default: 0.05)

C Initializing para CIM input generation

      allocate(nvir_clu(nclu),virtime(nclu))
      iforwhat=7
      call getival('nslv',nslv)
      call para_JobInit(iforwhat)

      call para_initsend
      call para_pack_int(nclu,1)
      call para_pack_int(natom,1)
      call para_pack_int(nbas,1)
      call para_pack_int(nmo,1)
      call para_pack_int(nocc,1)
      call para_pack_string(basline,100)
      call para_bcast_pack(TxCIMDat)

      call elapsec(tend1)
C Send data to slaves
      do i=1,nclu
         call para_recv(imsg,islave,TxDftReq)
         call para_initsend
         call para_pack_int(i,1)
         call para_pack_int(batom,nbas*natom)
         call para_pack_int(natm_clu(i),1)
         call para_pack_int(nb_clu(i),1)
         call para_pack_int(norb_clu(i),1)
         call para_pack_int(ML(:,i),nocc)
         call para_pack_int(BA(:,i),natom)
         call para_pack_int(ZA(:,i),nbas)
         call para_pack_real(SOVER,nbas*nbas)
         call para_pack_real(SMO,nbas*nmo)
         call para_pack_real(dis,natom*natom)
         call para_pack_int(link,natom*natom)
         call para_pack_int(nbatm,natom)
         call para_pack_int(atmclu,natom*natom*7)
         call para_pack_int(natmclu,natom*7)
         call para_pack_real(PAO,nbas*nbas)
         call para_pack_real(virt,1)
         call para_pack_real(FK,nbas*nbas)
         call para_pack_int(IAN,natom)
         call para_pack_int(clumtd(i),1)
         call para_pack_int(ncenorb(i),1)
         call para_pack_int(lenJ,1)
         call para_pack_string(jobname,256)
         call para_pack_string(AtSymb,natom*8)
         call para_pack_real(coor,3*natom)
         call para_send_pack(islave,TxCIMVir)
      enddo

      do i=1,nslv
         call para_recv(imsg,islave,TxDftReq)
         call para_initsend
         call para_pack_int(0,1)
         call para_pack_int(batom,nbas*natom)
         call Para_pack_int(natm_clu(i),1)
         call para_pack_int(nb_clu(i),1)
         call para_pack_int(norb_clu(i),1)
         call para_pack_int(ML(:,i),nocc)
         call para_pack_int(BA(:,i),natom)
         call para_pack_int(ZA(:,i),nbas)
         call para_pack_real(SOVER,nbas*nbas)
         call para_pack_real(SMO,nbas*nmo)
         call para_pack_real(dis,natom*natom)
         call para_pack_int(link,natom*natom)
         call para_pack_int(nbatm,natom)
         call para_pack_int(atmclu,natom*natom*7)
         call para_pack_int(natmclu,natom*7)
         call para_pack_real(PAO,nbas*nbas)
         call para_pack_real(virt,1)
         call para_pack_real(FK,nbas*nbas)
         call para_pack_int(IAN,natom)
         call para_pack_int(clumtd(i),1)
         call para_pack_int(ncenorb(i),1)
         call para_pack_int(lenJ,1)
         call para_pack_string(jobname,256)
         call para_pack_string(AtSymb,natom*8)
         call para_pack_real(coor,3*natom)
         call para_send_pack(islave,TxCIMVir)
      enddo

      allocate(natm_clu_slv(nclu),JF_slv(nclu),KB_slv(nclu))
      allocate(KV_slv(nclu),virtime_slv(nclu))

      natm_clu=0; nb_clu=0; norb_clu=0; nvir_clu=0; virtime=0.0D0
      do i=1,nslv
         call para_recv_pack(islave,TxCIMInfo)
         call para_unpack_int(natm_clu_slv,nclu)
         natm_clu=natm_clu+natm_clu_slv
         call para_unpack_int(JF_slv,nclu)
         nb_clu=nb_clu+JF_slv
         call para_unpack_int(KB_slv,nclu)
         norb_clu=norb_clu+KB_slv
         call para_unpack_int(KV_slv,nclu)
         nvir_clu=nvir_clu+KV_slv
         call para_unpack_real(virtime_slv,nclu)
         virtime=virtime+virtime_slv
      enddo

      do j=1,nclu
         virtime(j)=(virtime(j)+tend1-tstart)/60.0D0
      enddo

      write(iout,'(11x,26(''-''))')
      write(iout,'(11x,"FINAL CIM DOMAIN SUMMARY")')
      write(iout,'(11x,26(''-''))')
      write(iout,*) 
     & ' Clu_INDX  Natom   Nbas     Nocc     Nvir   Method  Elap_time'
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
      return
      end subroutine para_cimsubgen


C ****************************************************
C * subroutine for generating CIM cluster input file *
C * NZG_6/28/2016 @UARK                              *
C ****************************************************

      subroutine do_cimgen
      use memory
      use newpara

      implicit real*8 (a-h,o-z)
      
      integer clumtd
      integer,allocatable::batom(:,:),ML(:),BA(:),ZA(:),link(:,:)
      integer,allocatable::nbatm(:),atmclu(:,:,:),natmclu(:,:),IAN(:)
      integer,allocatable::CenMO(:),natm_clu_slv(:),JF_slv(:),KB_slv(:)
      integer,allocatable::KV_slv(:)
      real*8,allocatable::SOVER(:,:),SMO(:,:),MOS1(:,:),dis(:,:)
      real*8,allocatable::PAO(:,:),MO_clu(:,:),FK(:,:),dtmp(:)
      real*8,allocatable::FH(:,:),FHH(:),F2(:,:),TRANS(:,:),MOS2(:,:)
      real*8,allocatable::FIA(:,:),MAT(:,:),VECT(:,:),VALU(:),VC(:)
      real*8,allocatable::coor(:,:),virtime_slv(:)
      character line*100,basline*100
      character*256 cluname,jobname,inpname,infname,mosname
      character*8,allocatable::AtSymb(:)

      call para_recv_bcastpack(TxCIMDat)
      call para_unpack_int(nclu,1)
      call para_unpack_int(natom,1)
      call para_unpack_int(nbas,1)
      call para_unpack_int(nmo,1)
      call para_unpack_int(nocc,1)
      call para_unpack_string(basline,100)

      allocate(natm_clu_slv(nclu),JF_slv(nclu),KB_slv(nclu))
      allocate(KV_slv(nclu),virtime_slv(nclu))

      natm_clu_slv=0; JF_slv=0; KB_slv=0
      KV_slv=0; virtime_slv=0.0D0

      do
         call elapsec(tvir1)
         allocate(batom(nbas,natom),ML(nocc),BA(natom),ZA(nbas))
         allocate(SOVER(nbas,nbas),SMO(nbas,nmo),dis(natom,natom))
         allocate(link(natom,natom),nbatm(natom),atmclu(natom,natom,7))
         allocate(natmclu(natom,7),PAO(nbas,nbas),FK(nbas,nbas))
         allocate(IAN(natom),AtSymb(natom),coor(3,natom))
         call para_send(MY_GID,0,TxDftReq)
         call para_recv_pack(isource,TxCIMVir)
         call para_unpack_int(i,1)
         if (i==0) exit
         call para_unpack_int(batom,nbas*natom)
         call Para_unpack_int(natm_clu,1)
         call para_unpack_int(nb_clu,1)
         call para_unpack_int(norb_clu,1)
         call para_unpack_int(ML,nocc)
         call para_unpack_int(BA,natom)
         call para_unpack_int(ZA,nbas)
         call para_unpack_real(SOVER,nbas*nbas)
         call para_unpack_real(SMO,nbas*nmo)
         call para_unpack_real(dis,natom*natom)
         call para_unpack_int(link,natom*natom)
         call para_unpack_int(nbatm,natom)
         call para_unpack_int(atmclu,natom*natom*7)
         call para_unpack_int(natmclu,natom*7)
         call para_unpack_real(PAO,nbas*nbas)
         call para_unpack_real(virt,1)
         call para_unpack_real(FK,nbas*nbas)
         call para_unpack_int(IAN,natom)
         call para_unpack_int(clumtd,1)
         call para_unpack_int(ncenorb,1)
         call para_unpack_int(lenJ,1)
         call para_unpack_string(jobname,256)
         call para_unpack_string(AtSymb,natom*8)
         call para_unpack_real(coor,3*natom)

         allocate(MOS1(nbas,nbas))
         nprint=0
         call COV(nprint,i,natom,nmo,batom(1,:),natm_clu,nb_clu,nbas,
     &            norb_clu,nocc,nvir_clu,ML(:),BA(:),ZA(:),SOVER,SMO,
     &            MOS1,dis,link,nbatm,batom,atmclu,natmclu,PAO,virt)
         KB=norb_clu; KV=nvir_clu; NA=KB+KV; JF=nb_clu
         allocate(MO_clu(JF,NA))
         MO_clu(:,:)=MOS1(1:JF,1:NA)
         deallocate(MOS1)
C---------------------------------------------------------------------
C     Construct molecular Fock matrix of cluster
C---------------------------------------------------------------------

         allocate(FH(NA,NA),FHH(NA),F2(JF,JF),TRANS(KB,KB))
         FH=0.0d0; FHH=0.0d0; F2=0.0d0 
         do k=1,JF
            do j=1,JF
               F2(j,k)=FK(ZA(j),ZA(k))
            enddo
         enddo

         TRANS=0.0D0
         do j=1,KB
            TRANS(j,j)=1.0D0
         end do
C
C-WL- 2007.11.20 ADD FOR FOCK MATRIX DIAGONALIZATION --- liwei -----

         allocate(MOS2(JF,NA))
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
C First occupied MOs

         allocate(MAT(KB,KB),VECT(KB,KB),VALU(KB),VC(KB))
         MAT(:,:)=FH(:KB,:KB)
         call NJ_qr(nprint,MAT,KB,VECT,VALU,VC,k,1)
         FHH(1:KB)=VALU(1:KB)
C
         do j=1,KB
            do k=1,JF
               MOS2(k,j)=0.0d0
               do L=1,KB
                  MOS2(k,j)=MOS2(k,j)+MO_clu(k,L)*VECT(L,j)
               enddo
            enddo
         enddo
         TRANS=transpose(VECT)
C
         deallocate(MAT,VECT,VALU,VC)

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
            stop
         end if
C
         do j=1,KV
            do k=1,JF
               MOS2(k,j+KB)=0.0d0
               do L=1,KV
                  MOS2(k,j+KB)=MOS2(k,j+KB)+MO_clu(k,L+kB)*VECT(L,j)
               enddo
            enddo
         enddo
         deallocate(MAT,VECT,VALU,VC,MO_clu)
C
         if (nprint>0) write(nprint,"(1x,'Trace(Fock LMO)=',f20.10)")
     &                       dtrace2(NA,FH)
C
C --- Properties of the clusters 
C      nelec=nelec-mult+1 
         nelec=KB*2  !Since we only deal with ristricted system

         nucl=0
         do k=1,natm_clu
            j=BA(k)
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
         mwords=100 !set to 100 temporarily
         call SUBINP_slave(inpclu,NATOM,AtSymb,coor,natm_clu,BA,clumtd,
     &                     mwords,icharg,mult,basline)  
         close(inpclu)

C- Print out some important information in XXX.cim file
         open(inpclu,file=infname)
         write(inpclu,*) '$INFO'
         select case(clumtd)
            case(1)
               write(inpclu,'(a)') 'SUBMTD  RIMP2'
            case(2)
               write(inpclu,'(a)') 'SUBMTD  MP2'
            case(3)
               write(inpclu,'(a)') 'SUBMTD  CCD'
            case(4)
               write(inpclu,'(a)') 'SUBMTD  CCSD'
            case(5)
               write(inpclu,'(a)') 'SUBMTD  CCSD(T)'
         end select
C         if (IWORK(8).eq.0) then
         write(inpclu,'(a8,a8)') 'MOTYP   ','QCMO    ' 
C         else
C            write(inp,'(a8,a,a8)') 'MOTYP   ', '=', 'LMO     ' 
C         endif
         write(inpclu,'(a8,i8)') 'SYS     ', I
         write(inpclu,'(a8,i8)') 'NAT     ', natm_clu
         write(inpclu,'(a8,i8)') 'NELEC   ', nelec
         write(inpclu,'(a8,i8)') 'ICH     ', icharg
         write(inpclu,'(a8,i8)') 'MUL     ', mult
         write(inpclu,'(a8,i8)') 'NBAS    ', JF
         write(inpclu,'(a8,i8)') 'NOCC    ', KB
         write(inpclu,'(a8,i8)') 'NVIR    ', KV
         write(inpclu,'(a8,i8)') 'NMO     ', NA
         write(inpclu,'(a8,i8)') 'NCEN    ', ncenorb
         write(inpclu,*) '$END' 
C
         call iwrit8(inpclu,'$MO-OCC',KB,ML(1))
C         
         allocate(CenMO(KB))
         CenMO(1:ncenorb)=1       !Central MOs labeled as 1
         CenMO(ncenorb+1:KB)=0    !Else as 0
         
         call iwrit8(inpclu,'$MO-CEN',KB,CenMO(1)) !(CenMO(k,I)),k=1,KB)
         call rwrit8(inpclu,'$TRMX-A',KB*KB,TRANS(1,1))
         call rwrit8(inpclu,'$VEC',JF*NA,MOS2(1,1))
         call rwrit8(inpclu,'$EIGVAL',NA,FHH(1))
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
C
         deallocate(dtmp,TRANS)
         close(inpclu)
C
C -- Write QCMO coefficients and orbital energies to XXX.mos file
         itype=1
         call WriteMOS(JF,NA,MOS2,FHH,.true.,lenM,mosname,itype)

         deallocate(F2,MOS2,FH,FHH,CenMO,ML,batom,BA,ZA,SOVER,SMO,dis)
         deallocate(link,nbatm,atmclu,natmclu,PAO,FK,IAN,AtSymb,coor)

         call elapsec(tvir2)
         virtime=tvir2-tvir1
         natm_clu_slv(I)=natm_clu
         JF_slv(I)=JF
         KB_slv(I)=KB
         KV_slv(I)=KV
         virtime_slv(I)=virtime
      end do

      call para_initsend
      call para_pack_int(natm_clu_slv,nclu)
      call para_pack_int(JF_slv,nclu)
      call para_pack_int(KB_slv,nclu)
      call para_pack_int(KV_slv,nclu)
      call para_pack_real(virtime_slv,nclu)
      call para_send_pack(0,TxCIMInfo)

      return
      end subroutine do_cimgen


c  ################################################################
c  ##  Make PQS input of clusters                                ##
c  ##  Originally in GAMESS-2004.12.26; Update 2005.12.26 by WL  ##
C  ##  5/1/2016_NZG @UARK                                        ##
C  ##  7/1/2016_NZG @UARK Update for basis set input             ##
c  ################################################################
      subroutine SUBINP_slave(inp,nat,AtSymb,coor,snat,BA,clumtd,mwords,
     &                  icharg,mult,basline)
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      integer inp,nat,snat,BA(snat),i,j,k,l,m,n,m1,m2,clumtd
      real(kind=8) coor(3,nat)
      integer k1,k2,k3,k4,k5,k6,stat,mplevl,lenJ
      integer memnew,icharg,L1,L2,mult,L3,L4,mwords,L5,L6
      character*256 line,jobname
      character charg*100,mulp*100,basline*100
      character*8 AtSymb(nat)
      data IUnit   /1/ ! unit number for checkpoint I/O
      common /job/ jobname,lenJ

C
      rewind(inp)
C
C --- MEMORY
      write(line,*) mwords
      call NJ_trim(line,k1,k2)
      line='%MEM='//line
      write(inp,100) trim(line)
C
C --- GEOMETRY
      write(charg,*) icharg
      call NJ_trim(charg,L1,L2)
      write(mulp,*) mult
      call NJ_trim(mulp,L3,L4)
      line='GEOM=PQS '//'CHAR='//charg(L1:L2)//' MULT='//mulp(L3:L4)//
     &     ' SYMM=0'
      write(inp,100) trim(line)

      do k=1,snat
         L=BA(k)
         write(inp,120) AtSymb(L),(coor(j,L),j=1,3)
      enddo
C
C --- BASIS
       write(inp,100) trim(basline)
C      write(inp,100) 'BASIS NEXT'
C      open(unit=iunit,file=jobname(1:lenJ)//'.basis')
C      do while(.true.)
C         read(iunit,fmt="(A50)",iostat=stat) line
C         if (stat/=0) exit
C         if (line(1:1)=='$') then
C            cycle
C         else
C            write(inp,100) trim(line)
C         end if
C      end do
C
C --- METHOD
      select case(clumtd)
         case(1)
            line='CIMSUB RIMP2'
         case(2)
            line='CIMSUB MP2'
         case(3)
            line='CIMSUB CCD'
         case(4)
            line='CIMSUB CCSD'
         case(5)
            line='CIMSUB CCSD TRIPLES'
      end select
      write(inp,100) trim(line)

 100  format(a)
 120  format(A8,3f16.8)
C
      return
      end subroutine SUBINP_slave

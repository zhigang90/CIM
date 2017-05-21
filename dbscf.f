      subroutine dualscf
c
      use memory
c
      implicit real*8 (a-h,o-z)
ccc   common /intbl/maxsh,inx(100)
      common /job/jobname,lenJ
      parameter(nopt=3)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0,sixty=60.0d0)
      character*256 jobname,scrfile,filename,filname1,filname2,
     $             filname3,filname4
      character*256 chopv(nopt)
      character*3 ch3
      character opnames*4,scftype*11,wvfnc*20
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension xnmo(2)
c-----------------------------------------------------------------------
      data opnames  /'thre','type','prin'/
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c  11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp / 11,1,1/
c-----------------------------------------------------------------------
c  Explanation of the options:
c
c  THRE = integral thresold, n pH form, default is  9, meaning 1.0d-9
c  TYPE = two types 1 & 2 for different generation of occ & virt spaces
c         and different calculations of the E1 correction
c  PRIN = print level, default=2
c-----------------------------------------------------------------------
      data IUnit/1/
c-----------------------------------------------------------------------
c
      call secund(t_beg)
      call elapsec(e_beg) 
c
c  get input/output file names
c
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
c-----------------------------------------------------------
  45  format(/72('-'))
  46  format(/72('='))
c     write(iout,46)
c     write(iout,*) '                          The Dual SCF Module  '
c     write(iout,*)' '
c-----------------------------------------------------------
c
c  get symetry info: number of symm. op. starting address of contr. info
c
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c  basis function symmetry pairs are in inx(ifp)
c-----------------------------------------------------------------------
c --  read in current wavefunction and 
c --  at the same time get the lowest eigenvalue of the overlap)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,wvfnc)
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c-----------------------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
c
      write(iout,*)
     *  '                The Dual Basis Set SCF Module '
      write(iout,*)' '
      if(wvfnc.eq.'RHF ') then
         write(iout,*)
     *  '                            RHF              '
      endif
      if(wvfnc.eq.'UHF ') then
         write(iout,*)
     *  '                            UHF              '
      endif
      if(wvfnc.eq.'RDFT') then
         write(iout,*)
     *  '                          RHF/DFT            '
      endif
      if(wvfnc.eq.'UDFT') then
         write(iout,*)
     *  '                          UHF/DFT            '
      endif
      if(wvfnc.eq.'RMP2') then
         write(iout,*)
     *  '                          RHF/MP2            '
      endif
      if(wvfnc.eq.'UMP2') then
         write(iout,*)
     *  '                          UHF/MP2            '
      endif
      if(wvfnc.eq.'RMP2-dual') then
         write(iout,*)
     *  '                        RHF/Dual-MP2         '
      endif
      if(wvfnc.eq.'UMP2-dual') then
         write(iout,*)
     *  '                        UHF/Dual-MP2         '
      endif
      write(iout,*)' '
c
      call flush(iout)
c-----------------------------------------------------------------------
c
      call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c
c-----------------------------------------------------------------------
c  THRESHOLD (integral threshold)
      if(ifound(1).gt.0) then
        thresh=10.0d0**(-ropv(1,1))
      else
c --    check lowest eigenvalue of overlap (from SCF)
        thresh = MIN(1.0d-10,xlow**2)
        If(thresh.LT.1.0d-12) then
          write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final energy may have greater than mhartree error')
          thresh = 1.0d-12
        EndIf
      end if
c
      call setrval('thresh',thresh)
c .............................................................................
c   type of correction
c
      itype=1
      if(ifound(2).gt.0) then
         itype=iopv(1,2)
         if(itype.gt.3) itype=3
      endif
c
c .............................................................................
c    print
      iprint=0
      if(ifound(3).gt.0) iprint=iopv(1,3)
c-----------------------------------------------------------------------
       natoms=igetival('na')
c-----------------------------------------------------------
      ncs=igetival('ncs')
      ncf=igetival('ncf')
      ictr=igetival('ictr')
c-----------------------------------------------------------
      call getint(ncs,mpres_in_bg)
      call setival('mpres_in_bg',mpres_in_bg)
c-----------------------------------------------------------
c        the current "extended" basis set was already read-in
c        after SCF. Get the smaller bais set info (one used for SCF):
c
         call get_small_basis(bl,ncs,ncf, bl(ictr), ncs_sm,ncf_sm)
c
c-----------------------------------------------------------
         write(iout,45)
         write(iout,*)
     *   '              Dual basis set will be used in SCF calculations'
         write(iout,*) ' '
         write(iout,47) ncs_sm,ncf_sm
         write(iout,48) ncs   ,ncf
   47    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for occupied MOs')
   48    format(7x,i5,' shells & ',i5,
     *    ' contracted functions for virtual  MOs')
         write(iout,45)
         write(iout,49) itype
         write(iout,45)
   49    format(22x,'dual basis set approach type = ',i2)
c-----------------------------------------------------------
      write(iout,*) ' SCF integral thresh    = ',thresh
      write(iout,*) ' '
      call flush(iout)
c-----------------------------------------------------------
c  put down a memory marker
      call matreset
c
cccc  call immark
      call mmark
      call matmark
c-----------------------------------------------------------
c
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
          call nerror(2,'dualscf ',
     1   'Program cannot find the number of MOs on <jobname>.14',
     2    ni,ni)
        end if
      end if
      nmo=xnmo(1)
      erhf=xnmo(2)
c-----------------------------------------------------------
c get memory for MOs, orbital energies and screening density
c
      np4=4
      call matdef('cano','q',ncf,ncf)
      call matdef('epsi','d',ncf,ncf)
      icano=mataddr('cano')
      iepsi=mataddr('epsi')
c-------------------------------------------------------
c
c get eigenvectors & eigenvalues
c
c       call get_mix_eigen(bl,nblocks,inx(ictr),labels,thresh,
c    *                     ncf,ncf_sm,ncs,ncs_sm,nmo)
c
        call get_mix_eigen2(bl,nblocks, bl(ictr),labels,thresh,
     *                     ncf,ncf_sm,ncs,ncs_sm,nmo,itype)
c-------------------------------------------------------
        nfirst=1
        nlast=nmo
        nval=nmo
c---------------------------------------------------------------------
c  establish orbital symmetry characters
c     if(nsym.gt.0) then
c       call getint(nval*nsym,iorbpair)
c       call getmem(ncf,itmp)
c       call OrbSymPair(nsym, ncf ,   bl(icano), inx(ifp) ,nfirst,
c    1                   nlast,iprint, bl(itmp),  inx(iorbpair))
c       call retmem(1)
c     endif
c---------------------------------------------------------------------
c  print the data of the calculation early
      write(iout,'("  Number of contracted functions =",i8,/,
     1             "  Number of occupied orbitals    =",i8,/,
     2             "  Number of virtual orbitals     =",i8)')
     3                ncf,nval,ncf-nmo
      write(iout,*) ' '
      call flush(iout)
c.................................................
      call  matsub('occu','cano',nfirst,nlast)
      call  matsub('virt','cano',nmo+1,ncf)
      ioccu=mataddr('occu')
c.................................................
      if(iprint.gt.3) then
        call matprint ('occu',6)
        call matprint ('virt',6)
      end if
c--------------------------------------------------------------------
c Print out the final results :
c
      call getrval('esing_exc',esing_exc)
      edual = erhf + esing_exc
c
      write(iout,*) '-------------------------------------------------'
      write(iout,*) '               Final DUAL SCF results '
      write(iout,*) ' '
      write(iout,3001) erhf
      write(iout,3002) esing_exc
      write(iout,*) '-------------------------------------------------'
      write(iout,3003) edual
      write(iout,*) '-------------------------------------------------'
c
 3001 format('        Small basis set SCF energy =',f15.8)
 3002 format('        Single excitations energy  =',f15.8)
 3003 format('        Corrected SCF energy       =',f15.8)
c
      call flush(6)
c---------------------------------------------------------------------
c     xnmo(2)=edual
c     call wri(xnmo,2,np4,-1,'nocc_rhf')
c     call matwrite('cano',np4,0,'evec_rhf')
c     call matwrite('epsi',np4,0,'eval_rhf')
c---------------------------------------------------------------------
c
c -- set wavefunction type
c
      wvfnc = 'RHF-dual'
c
c -- write final daual/scf energy to <control> file
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      Call wrcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
      Call wrcntrl(IUnit,6,'$edual' ,2,idum,edual,wvfnc)
      Call wrcntrl(IUnit,6,'$escfS' ,2,idum,erhf ,wvfnc)
      Call wrcntrl(IUnit,6,'$escf1' ,2,idum,esing_exc ,wvfnc)
      Call wrcntrl(IUnit,6,'$escfB' ,2,idum,edual ,wvfnc)
      Call wrcntrl(IUnit,7,'$energy',2,idum,edual,wvfnc)
      Call wrcntrl(IUnit,5,'$escf',2,idum,edual,wvfnc)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c---------------------------------------------------------------------
      call secund(t_end)
      call elapsec(e_end) 
c
      total=(t_end-t_beg)/60.0d0
      elaps=(e_end-e_beg)/60.0d0
c
      write(iout ,400) total,elaps
  400 format('Master CPU time for dual basis set scf ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c---------------------------------------------------------------------
      call matremark
      call retmark
cccc  call retimark
c---------------------------------------------------------------------
      call memory_status('end of dual scf')
c---------------------------------------------------------------------
      end
c=======================================================================
      subroutine get_mix_eigen2(bl,nblocks,inx,labels,thresh,
     *                         ncf,ncf_sm,ncs,ncs_sm,nmo,idualtype)
c
      use memory, block => bl
c
      implicit real*8(a-h,o-z)
      dimension bl(*)
      dimension xintxx(9)
      dimension inx(12,*)
      dimension field(3)
ccccc common /intbl/ifpp,inxx(100)
      logical rhf
      character scftype*11
      character*256 jobname, filename
      common /job/jobname,lenJ
c----------------------------------------------------------------
c different approachs :
c
c (1)  occupied from small scf i.e. from diagonalization of f(d) in small
c      basis set converged scf, virtuals from diagonalization of PF(d)P
c      where projector P = I-(1/2)dS with d being small density augmented
c      with zeros to match big basis set dimension. Thus,
c     
c      occupied from diag. of f(d), virtuals from diag.of PF(d)P
c      Esing_exc=E_bsc= - SUM(occ,vir) (Fov)**2/(Ev -Eo)
c
c
c (2) occupied from diag. of f(d) (small scf), virtuals from diag. of F(d)
c     Esing_exc=Tr[ F(d) ]*(D-d) where D from diagon. of F(d)
c
c----------------------------------------------------------------
c Input :
c bl(*)     - general storage
c nblocks   - total number of blocks (of quartets) for the old
c             integral program used here
c label     - a pointer to the integral's labels array
c inx(12,*) - general basis set info
c ncf, ncs  - dimensions of a big basis set
c ncf_sm,ncs_sm - dimensions of a small basis sets
c----------------------------------------------------------------
c This routine is used only for "dual basis set" MP2 .
c It returns eigenvectors & eigenvalues needed for MP2.
c For occupied orbitals MOs & Mos energies are those from
c ths "small basis set" SCF. For virtuals they are obtained
c by special diagonalization of a modified Fock matrix. "Special"
c means that the mixing between occupied & virtuals is not allowed.
c The modified Fock matrix is defined in the big basis set but the
c Fock operator is the Fock operator of the small basis, i.e.
c the small basis set density is used to build the Fock operator
c
c There are the following steps :
c
c (1) construct "big basis set" density by projection of the small
c     one (& zeros where nedeed)
c     Such density is needed to build the special Fock matrix
c     and it will be used to pre-screen integrals. Integrals
c     needed now have two indecies in a small basis set and
c     two others in a big basis set.
c (2) construct 2-el. part of the special Fock matrix
c     it is done by calling "old" integral program (int_fock
c     renamed here to int_fmp2) running over the big basis set
c     but the corresponding density used in prescreening should
c     eleminate all unneccesary integrals involving shells/functions
c     which were not present in a small basis set.
c     Dkl *{ (ij|kl) - 0.5(il|kj) } where kl belong to small bs only.
c (3) calculate 1-el. part of the full Fock matrix . Call inton
c     to get H0 (over big basis set) and add it to the previously
c     obatained 2-el. part.
c (4) diagonalize the Fock matrix in such a way that occupied and
c     virtuals do not mix. From that we take virtuals only i.e.
c     MOs & energise. The occupied ones are those obtained in ths
c     small basis set SCF.
c (5) calculate "basis set" correction to the energy :
c     E_bsc= - SUM(occ,vir) (Fov)**2/(Ev -Eo)
c
c----------------------------------------------------------------
      call mmark
c----------------------------------------------------------------
         call secund(tgmix0)
         call elapsec(egmix0)
c----------------------------------------------------------------
      np1=1
      np4=4
c
      iout=igetival('iout')
c----------------------------------------------------------------
c Make "big basis set" density :
c transfer density dens(ij) from the smalll basis set to the big basis.
c Put zeros where nedeed.
c......................................................
c allocate memory for "big"  and small density matrices
c
      ntri=ncf*(ncf+1)/2             ! big basis set
      ntri_sm=ncf_sm*(ncf_sm+1)/2    ! small basis set
c
      call matdef('densb','s',ncf,ncf)
      idensbig=mataddr('densb')
c
      call matdef('densm','s',ncf_sm,ncf_sm)
      idensmal=mataddr('densm')
c......................................................
c get small basis set density matrix
c
      call rea (bl(idensmal),ntri_sm,np4,'den0_rhf')
c
cccc  write(6,*)' ncf_sm=',ncf_sm,' ncf_big=',ncf
cccc  call matprint('densm',6)
c......................................................
c make projection of dens. from small to big basis set
c
      call zeroit(bl(idensbig),ntri)
c
      call getival('maf_s2b',maf_s2b)
      call dens_s2b(bl(idensmal),ncf_sm,bl(idensbig),ncf,bl(maf_s2b))
ctest
c     write(6,*)' Small basis set density matrix'
c     call matprint('densm',6)
c     write(6,*)' Big   basis set density matrix'
c     call matprint('densb',6)
c......................................................
c save big basis set density for NMR use
c
      call matwrite('densb',np4,0,'den0_rhf')
c......................................................
c
      call matrem('densm')
c----------------------------------------------------------------
c DO 1-EL PART which involves overlap S
c
      call getival('ibas',ibas)       ! big basis set
      call getival('inuc',inuc)       ! big basis set
      call getival('na  ',na)
c
c............................................................
c calculate S matrix in a big basis set
c
      call matdef('smat','s',ncf,ncf)
      ismat=mataddr('smat')
      call inton(0,na,bl(ismat),inx,0,0,bl(ibas),bl(inuc),ncs)
c
c save S matrix on a disk
      call matwrite('smat',np1,0,'s matrix')
c............................................................
c Make  projector = I-(1/2)DS 
c where I = unit mtx, D = small basis density
c (projected to the large basis), S = (large basis) overlap
c The matrix "u" is used as scratch for storing the projector
c
c make projector (in u) and write it on a disk:
c
c
      call matdef('u','q',ncf,ncf)
c
      if(idualtype.eq.1) then
         call make_proj('densb','smat','u')
         call matwrite('u',np4,0,'projector')
      endif
c
c............................................................
c calculate the LU factorization of S 
c
      call cholesky('smat','u',ncf,xlow)
c     
c----------------------------------------------------------------
c
      lastvirt=igetival('nonredun')
c
      if(lastvirt.ne.ncf_sm) then
         write(iout,*)'  '
         write(iout,*)
     * 'There was a basis set suppresion in a small basis set SCF'
         write(iout,66) ncf_sm,lastvirt
  66     format(' number of non-redundant basis functions reduced from '
     *          ,i4,' to ',i4)
         if(idualtype.eq.1) then
           write(iout,*)
     *'the small basis set SCF energy correction might be in error'
         write(iout,*)'  '
         endif
      endif
c
      call setival('nonredun',ncf)
c
c----------------------------------------------------------------
c write U matrix on disk
c
      call matwrite('u',np4,0,'umatrix')
c
c............................................................
      call matrem('u')
      call matrem('smat')
c----------------------------------------------------------------
c Check if the integral threshold should be sharpen:
c check lowest eigenvalue of overlap in a big basis set
c
      if(thresh.gt.xlow**2) then
        ithresh=-log10(thresh)
        if(ithresh.eq.10) then     ! default
           thresh = xlow**2
           thresh = max(thresh, 1.0d-12)
           write(iout,2000) thresh
        else
           write(iout,2100)
        endif
      endif
      call flush(iout)
 2000 format(/'  Near Linear Dependency in Big Basis Set',/
     *  '  Sharpening integral threshold to ',e12.4)
 2100 format(/'  Near Linear Dependency in Big Basis Set',/
     *  '  integral threshold would be sharpened',/
     *  '  (if not overwritten from the input)')
c----------------------------------------------------------------
c allocate memory for fock
c
      call matdef('fockmp2','s',ncf,ncf)
      ifockmp2=mataddr('fockmp2')
c
      call zeroit(bl(ifockmp2),ntri)
c----------------------------------------------------------------
         call secund(tgmix1)
         call elapsec(egmix1)
c----------------------------------------------------------------
c zero out irrelevant dft stuff
      idft=0
      ax=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
c----------------------------------------------------------------
c initialize the two-el. integral program
c
c    change for scf integrals iforwhat=1
c
      iforwhat=1
      nfock=1
      thint=thresh
      call getival('nslv',nslv)           ! no. of slaves
      If(nslv.GT.0) call para_JobInit(iforwhat)
c
      call ptwoint(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             0,
     *             scftype,xintxx, nblocks,maxbuffer,maxlabels)
c----------------------------------------------------------------
c Calculate 1-el. part of the fock matrix : H0
c Construct the full H0 matrix in big basis set
c
      call getival('ibas',ibas)       ! big basis set
      call getival('inuc',inuc)       ! big basis set
      call getival('na  ',na)
c
      call matdef('h0mat','s',ncf,ncf)
      ihmat=mataddr('h0mat')
      call zeroit(bl(ihmat),ntri)
      call para_oneint(1,na,bl(ihmat),inx,0,0,bl(ibas), bl(inuc),ncs)
ctest
cprt  write(6,*)' H0 matrix '
cprt  call matprint('h0mat',6)
c
c----------------------------------------------------------------
c add an electric field
c
c -- see if field defined earlier (e.g. in polarizability module)
c
      call getival('field',ifield)
      if(ifield.ne.0) then
        call getrval('fieldX',field(1))
        call getrval('fieldY',field(2))
        call getrval('fieldZ',field(3))
c
c       write(6,*)' External electric field on '
c       write(6,*)'  components : ', field(1),field(2),field(3)
c
        call matdef('dipol','s',ncf,ncf)
        idipo=mataddr('dipol')
c
c       Use 'dipol' as temporary matrix
        do icomp=1,3
          if(abs(field(icomp)).gt.1.0d-6) then
            call inton(3,na,bl(idipo),inx,icomp,0,bl(ibas),
     1                 bl(inuc),ncs)
            call matadd1('dipol',field(icomp),'h0mat')
          end if
        end do
        call matrem('dipol')
      end if
c----------------------------------------------------------------
c add pseudopotential contribution if necessary
c
      call getival('npsp',npsp)
      if(npsp.ne.0)then
        call matdef('h0psp','s',ncf,ncf)
        ihpsp=mataddr('h0psp')
        call zeroit(bl(ihpsp),ntri)
        call getrval('ithr',psptol)
        psptol=psptol*0.01d0       ! psp thres.=0.01*main int. thres
        call psph0(na,npsp,ncf,ncs,inx,bl(inuc),bl(ibas),
     $             psptol,bl(ihpsp))
        call matadd('h0psp','h0mat')
        call matrem('h0psp')
      endif
c
c----------------------------------------------------------------
c save H0 matrix
c
cprt  call matprint('h0mat',6)
      call matwrite('h0mat',np4,0,'h0mat')
      call matrem('h0mat')
c
c----------------------------------------------------------------
c   COSMO section
c
      call tstival('cosmo',icosmo)    !if COSMO flag is not defined,
      If(icosmo.EQ.0) Then            ! define it now
        call setival('cosmo',0)
      Else
        call getival('cosmo',icosmo)
      EndIf
      if(icosmo.ne.0)then
        call getival('c_nratom',nratom) ! no. of real atoms (no ghost)
        call getival('c_nps',nps)       ! no. of surface segments
        call getival('c_npsphe',npspher)! no. of outer surface segments
        call getival('nslv',nslv)       ! no. of slaves
        call getival('ibas',ibas)       ! big basis set
        call mmark
        call getmem(3*nratom,icoxyz)    ! atomic ccordinates
        call getmem(3*(nps+npspher),icosurf) ! surface coordinates
        call getmem(nps,iar)            ! surface area
        call getmem(nps,iqcos)          ! surface charges
        call getmem((nps+1)/2,iiatsp)   ! mapping surface --> atom
        call getmem((nratom+1)/2,icharge) ! nuclear charges
c
c  read arrays from COSMO data file
c
        call cosmo_rdata(bl(icoxyz),bl(icosurf),bl(iar),bl(iqcos),
     $                   bl(iiatsp),bl(icharge),nratom,nps)
        do i=0,3*npspher-1
           bl(icosurf+3*nps+i)=bl(icosusurf+i)
        enddo
c
        call getmem(ntri,ivimat) ! electron-surface repulsion
        call getmem(ntri,ih0cos) ! cosmo contribution to H0
        call getmem(max(nps,npspher),iphi) ! surface potential
        call getmem(max(nps,npspher),iphin)! nuclear part of surf. pot.
c
c   compute nuclear part of surface potential
c
        call cosmo_potn(bl(icoxyz),bl(icharge),bl(iphin),
     $                  bl(icosurf),nratom,nps)
        if(nslv.eq.0)then
          call cosmo_surfrep(bl(ivimat),bl(icosurf),inx,bl(ibas),
     $                       nps,ncs,ncf,ntri)
        else
          call para_cosmo_data(bl(icosurf),bl(icoxyz),
     $                         nslv,nratom,nps,npspher)
          call para_cosmo_surfrep(nslv,nps)
        endif
      endif
c----------------------------------------------------------------
c setup parameters for int_fock (here int_fmp2) call :
c
      call getmem(maxlabels,labels)
c
      nfock=1
      thres1=thresh
      mywork=0
      igran=0
c
      idft=0
      ax=0.0d0
      rhf=.true.
c     fockB=
c     denB=
      mgo=1
c
      call pfock(idft,ax,nblocks,nfock,rhf,ncf,bl,inx,thres1,
     *                 mgo,bl(idensbig),bl(ifockmp2),fockB,
     *                 bl(idensbig),DenB,bl(labels))
c
ctj      call int_fmp2(nblocks,nfock,ncf,bl,inx,thres1,
ctj     *              bl(idensbig),bl(ifockmp2),bl(idensbig),
ctj     *              bl(labels),mywork,igran)
c................................................................
c NOTE : 
c fock matrix on return is re-scaled by integ.threshold
c for parallel impl. call pfock instead of int_fock (int_fmp2)
c and re-scaling will be done in pfock (para_fock)
ctest
c     write(6,*)' 2-el part of MP2-Fock matrix'
c     call matprint('fockmp2',6)
c................................................................
c
      if(icosmo.ne.0)then
c
c  COSMO potential and one electron contribution
c
        if(nslv.eq.0)then
          call cosmo_pot(bl(idensbig),inx,bl(ibas),bl(iphin),
     $                   bl(iphi),bl(icosurf),bl(ivimat),
     $                   nps,ncf,ncs,ntri)
        else
          call para_cosmo_pot(bl(idensbig),bl(iphin),bl(iphi),
     $                        nslv,nps,ntri)
        endif
        call getrval('c_fepsi',fepsi)
        if(nslv.eq.0)then
          call cosmo_h0(bl(ih0cos),bl(ivimat),inx,bl(ibas),
     $                bl(icosurf),bl(iqcos),fepsi,nps,ncf,ncs,ntri)
        else
          call para_cosmo_h0(bl(ih0cos),bl(iqcos),fepsi,nslv,nps,ntri)
        endif
        call addvec(ntri,bl(ifockmp2),bl(ih0cos),bl(ifockmp2))
      endif
c
c  this is to fool para_next into not requesting the timing data
c  at this stage
c
      call para_next(5)
c
c  if this is a parallel run with cosmo, the slaves want to
c  compute the OC correction now, so ve have to fool them
c
      if(icosmo.ne.0)then
        if(nslv.ne.0)then
          call para_cosmo_surfrep(nslv,0)
          call para_cosmo_pot(bl(idensbig),bl(iphin),bl(iphi),
     $                        nslv,npspher,ntri)
        endif
        call cosmo_del('c_onel  ')
        call retmark
      endif
c
c   now it is really over, we can get the timing data from the slaves
c
      mgo=0
      call para_next(mgo)
c......................................................
         call secund(tgmix2)
         call elapsec(egmix2)
c----------------------------------------------------------------
c form the full Fock matrix
c
      call matdef('h0mat','s',ncf,ncf)
      call matread('h0mat',np4,'h0mat')
c
      call matadd('h0mat','fockmp2')
c
      call matrem('h0mat')
c
ctest
c     write(6,*)' Final Fock matrix'
c     call matprint('fockmp2',6)
c............................................................
c save 'fockmp2' (needed later for basis set correction):
c
      call matwrite('fockmp2',np4,0,'fockmp2 ')
c
c save 'fockmp2' under 'fock_rhf' like from SCF (may be needed in NMR)
c
      call matwrite('fockmp2',np4,0,'fock_rhf')
c----------------------------------------------------------------
c The final step : diagonalize full fock without
c                  mixing between occupied & virtual
c .................................................................
c  allocate memory
c
c use cano & epsi defined already before
c
      call matdef('u','q',ncf,ncf)
      call matdef('uinv','q',ncf,ncf)
c............................................................
c project out the occupied orbital subspace from the basis set
c everything below is in the large basis set
c projector = I-(1/2)DS where I = unit mtx, D = small basis density
c (projected to the large basis), S = (large basis) overlap
c
      if(idualtype.eq.1) then
         call matread('u',np4,'projector') ! get projector (in u) 
         call matsimtr('fockmp2','u','fockmp2') ! make projection
      endif
c
c............................................................
c get  U matrix from a disk and find its inverse
c
      call matread('u',np4,'umatrix')
      call u_inverse('u','uinv',ncf)
c
c............................................................
c Make the final diagonalization in order to get all virtuals
c
      zero=0.d0
      call geneig('fockmp2','u','uinv','epsi','cano','densb',zero,0,'U')
c
c above on entry fockmp2 was either projected Fock PF(d)P (itype=1)
c or full Fock F(d) (itype=2 or 3)
c............................................................
c For itype=1 :
c Define the virtual subspace of this. The occupied orbitals
c give zero eigenvalues. What happens if there are negative
c orbital energy virtuals? one will have to check the epsilons.
c
      nmo=igetival('nmo')
      if(idualtype.eq.1) then
         do i=1,nmo
            call matelem('epsi',i,i,x)
            if(abs(x).gt.1.0d-8) then
         write(6,*)'one of the first nmo eigenvalues is non-zero'
         write(6,*)' it is no=',i,' with a value of ',x           
              call nerror(1,'get_mix',
     1          'one of the first nmo  eigenvalues is non-zero',i,i)
            endif
         end do
      endif
c............................................................
c
c Dual basis set energy correction according to Liang &Head-Gordon :
c
c E=Esscf + Tr[ F(d) ](D-d)
c
c............................................................
      if(idualtype.eq.2 .or. idualtype.eq.3) then
c         construct new density from eigenvectors of F(d)
c         and calcualte :
c         Esing_exc= Tr[ F(d) ](D-d)    type=2
c
          call matdef('densnew','s',ncf,ncf)
          call densma('cano','densnew',nmo,.true.)
c         call matprint('densnew',6)
c
c save big basis set density for NMR use
c
          call matwrite('densnew',np4,0,'den0_rhf')
c
          call matadd1('densnew',-1.0d0,'densb')   ! d-D  in 'densb'
          call matrem('densnew')
c         call matprint('densb',6)
c
          call spur(bl(ifockmp2),bl(idensbig),ncf,esing_exc)
c
          esing_exc=-esing_exc   !  Liang&Head-Gordon
c
          call setrval('esing_exc',esing_exc)
c
          write(6,*)' Esing_exc= Tr[F(d)](D-d)=',esing_exc
      endif
c............................................................
c
      call matrem('uinv')
      call matrem('u')
c
ctest
c     write(6,*)' Big basis set MOs & orbital energies'
c     call matprint('cano',6)
c     call matprint('epsi',6)
c----------------------------------------------------------------
c Form the final eigenvectors & eigenvalues for dual basis set SCF
c
c get small basis set eigens :
c
      call matdef('canosm','q',ncf_sm,ncf_sm)
      call matdef('epsism','d',ncf_sm,ncf_sm)
      call matread('canosm',np4,'evec_rhf')
      call matread('epsism',np4,'eval_rhf')
c
c     write(6,*)' Small basis set MOs & orbital energies'
c     call matprint('canosm',6)
c     call matprint('epsism',6)
c
      icoef_sm=mataddr('canosm')
      iepsi_sm=mataddr('epsism')
c
c big basis set eigens :
c
      icoef_bg=mataddr('cano')
      iepsi_bg=mataddr('epsi')
c
      IF(idualtype.eq.1 .or. idualtype.eq.3) THEN
         call make_mp2_eigens(nmo,bl(maf_s2b),
     *                         bl(icoef_sm),bl(iepsi_sm),ncf_sm,
     *                         bl(icoef_bg),bl(iepsi_bg),ncf   )
      ENDIF
c
c     call matprint('cano',6)
c     call matprint('epsi',6)
c............................................................
c save eigenvectors and eigen values on disk for NMR use
c
c     IF(idualtype.eq.2) THEN
         call matwrite('cano',np4,0,'evec_rhf')
         call matwrite('epsi',np4,0,'eval_rhf')
c     ENDIF
c............................................................
c............................................................
c -- set up the MOs filename
c     call tstchval('mosfname',iyes)
c     if(iyes.eq.1) then
c       call getchval('mosfname',filename)
c     else
c       filename = jobname(1:lenJ)
c     endif
c     call rmblan2(filename,256,lenM)
c     if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
c     lenM=lenM+4
c
c -- read in the MOs
c     itype = 1
c     CALL ReadMOS(ncf,bl(icano),bl(iepsi),.True.,lenM,
c    $             filename(1:lenM),itype,IErr)
c............................................................
c............................................................
      filename = jobname(1:lenJ)
      call rmblan2(filename,256,lenM)
      if(lenM.le.252) filename(lenM+1:lenM+4)='.mos'
      lenM=lenM+4
      call WriteMOS(ncf,    ncf ,   bl(icoef_bg),bl(iepsi_bg),.True.,
     $              lenM,filename(1:lenM), 1)
c............................................................
c on return 'cano' & 'epsi' contain :
c
c For itype=1 
c occupied MOs from small & virtuals from diagonalization of PF(d)P
c
c For itype=2 it is irrelevant for scf energy correction but needed for mp2
c occupied MOs  & virtuals from diagonalization of F(d)
c
c For itype=3 
c occupied MOs from small & virtuals from diagonalization of F(d)
c
ctest
c     write(6,*)' Final MOs & orbital energies for dual basis set SCF'
c     call matprint('cano',6)
c     call matprint('epsi',6)
c............................................................
c release memory
c
      call matrem('epsism')
      call matrem('canosm')
c
c----------------------------------------------------------------
c Calculate "basis set correction" to the SCF energy :
c
c     E_bsc= - 2*SUM(occ,vir) (Fov)**2/(Ev -Eo)
c
c first make <OCC|FOCK|VIRT> block of fockmp2 :
c
c Fov  = Cocc(T) * Fao * Cvir
c
c get back 'fockmp2' :
c
      call matread('fockmp2',np4,'fockmp2 ')
c     call matprint('fockmp2',6)
c
      IF(idualtype.eq.1) THEN
c
      nvirt=ncf-nmo
      call matdef('fockov','r',nmo,nvirt)
      call matsub('occupied','cano',1,nmo)
      call matsub('virtuals','cano',nmo+1,ncf)
      call matdef('tempor','r',nmo,ncf)
c
      call matmmul2('occupied','fockmp2','tempor','t','n','n')
      call matmmul2('tempor','virtuals','fockov','n','n','n')
c
      call matrem('tempor')
      call matrem('virtuals')
      call matrem('occupied')
c
ctest
c     write(6,*)' < occ |Fock|vir> block of F'
c     call matprint('fockov',6)
c............................................................
c calculate correction :
c
      ifockov=mataddr('fockov')
      call calc_bsc(nmo,nvirt,bl(ifockov),bl(iepsi_bg),ncf,esing_exc)
      call setrval('esing_exc',esing_exc)
      call matrem('fockov')
c
c     write(6,*)' Esing_exc=E_bsc=-2*SUM(occ,vir)(Fov)**2/(Ev -Eo)=',
c    *    esing_exc
c
      ENDIF
c----------------------------------------------------------------
c integral's timings :
c
      tgmix=(tgmix2-tgmix1)/60.d0
      egmix=(egmix2-egmix1)/60.d0
c
c total timings :
c
      call secund(tgmix3)
      call elapsec(egmix3)
      tgmit=(tgmix3-tgmix0)/60.d0
      egmit=(egmix3-egmix0)/60.d0
c
      write(iout,*) ' '
      write(iout,101)
  101 format(' CPU & Elapsed timings in dual basis set Fock ')
      write(iout,*) ' '
      write(iout,102) tgmix, egmix
  102 format('  Integral calculations =',f8.2,' and ',f8.2,' minutes')
      write(iout,103) tgmit, egmit
  103 format('  Total for Fock build  =',f8.2,' and ',f8.2,' minutes')
      write(iout,*) ' '
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------
      end
c======================================================================

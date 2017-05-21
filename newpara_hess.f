C ---------------------------------------------------------------------
C  Parallel Hessian Routines           Jon Baker    July 2001
C ---------------------------------------------------------------------
C
      SUBROUTINE Para_HessInit(ncache,iforwhat,nblocks)

      use newpara

      IMPLICIT REAL*8(A-H,O-Z)
C
C  initializes the slaves for parallel Hessian calculation
C
      Logical rhf
c
c check if the slaves are OK
c
      call para_check
c
c -- get data from the depository
      call getival('NAlpha',NAlpha)
      call getival('NBeta',NBeta)
      call getival('Multip',IMult)
c
      rhf = NBeta.EQ.0.AND.IMult.EQ.1
c
c tell the slaves that analytical Hessian calculation is next
c
      call para_bcast(TxDoHess,TxJobType)
c
c Obtain and pack slave initialization data
c
      call para_initsend
      call para_pack_int(ncache,1)
      call para_pack_int(iforwhat,1)
      call para_pack_int(NAlpha,1)
      call para_pack_int(NBeta,1)
      call para_pack_int(rhf,1)
c
      call getrval('threg2',threg2)
      call getrval('thref1',thref1)
      call getrval('thref2',thref2)
c
      call para_pack_real(threg2,1)
      call para_pack_real(thref1,1)
      call para_pack_real(thref2,1)
c
      call para_bcast_pack(TxHessInit)
c
c Now send the remaining basic data necessary (geometry,cpu,symmetry)
c
      call para_send_info
c
      end
c ..........................................................................
c
      subroutine para_d0g1(idft,ax,rhf,nblocks,bl,inx,ntri,thref1,
     *                      natb,nate,listreal,do_allat,
     *                     densp, ncenter,
     *                     dens,denB,fder,fderB,labels)

      use newpara

      implicit real*8 (a-h,o-z)
c
c This routine calculates integral 1st derivatives (for derivative Fock)
c block by block and then calculates their contribution to the Hessian.
c Actually, this routine handles the PVM calls and distributes tasks
c among nodes. 
c---------------------------------------------------------------------
c integral first derivatives are calculated in int_d0g1()
c---------------------------------------------------------------------
c INPUT :
c  idft              - dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  rhf               - logical flag for closed-shell
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  ntri              - for dimensioning
c  thref1            - integral threshold
c  natb,nate         - for atoms from NATB to NATE
c  densp             - screening density
c  ncenter           - center of each basis function
c  dens              - closed-shell or alpha density matrix
c  denB              - beta density matrix
c  fder              - on exit contains contribution to alpha derivative Fock
c  fderB             - on exit contains contribution to beta derivative Fock
c  labels            - labels buffer
c---------------------------------------------------------------------
      character*11 scftype
      character*4 where
      Logical rhf,do_allat
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),dens(*),denB(*),fder(*),fderB(*)
      dimension densp(*)
      dimension listreal(*)
      dimension ncenter(*)
      dimension labels(*)
c
      call getival('ncs ',ncs)
      call getival('gran',igran)
      call getival('na  ',na)
      call getival('natunq',natunq) ! number of unique atoms (Abelian)
c     nat3 = 3*na
c
      if(scftype.ne.'full-direct') then
         call nerror(1,'para_d0g1',
     *                  'can only run in full-direct mode',0,0)
      endif
c
ckwkwkwkwkwkwkwkwkw
c
      natonce=nate-natb+1
c
c       write(*,*)' pvm_hess: natb,nate=',natb,nate,' na=',na
c       write(*,*)'                       natunq        =',natunq
c
      call para_initsend
      call para_pack_int(natb,1)
      call para_pack_int(nate,1)
ckw
      call para_pack_real(densp,2*ncs*ncs)
      call para_pack_int(listreal,na)
      call para_pack_int(do_allat,1)
ckw
      call para_bcast_pack(TxHessDat)
c
ckwkwkwkwkwkwkwkwkw
c
c Give work assignments to slaves until all work is done
      call para_distr_blocks(nblocks,igran)
c
c do global sum of results
c
cc        call elapsec(ttt1)
      call para_reduce(fder,3*ntri*natonce,TxHess2)
c -- now get back beta derivative Fock matrix
      If(.NOT.rhf) call para_reduce(fderB,3*ntri*natonce,TxHess3)
c ---------------------------------------------------------------------
cc        call elapsec(ttt2)
cc        ttt1 = (ttt2-ttt1)/60.0d0
cc        write(6,*) ' time spent in <pvmfreduce> is:',ttt1,' min.'
cc        call f_lush(6)
c
c
c  tell slaves we are done :
c
ckw
      nattodo=na
      if(.not.do_allat) nattodo=natunq
ckw
c     write(6,*) ' MASTER:  nate:',nate,' na:',na,' unique=',natunq
c     write(6,*) ' MASTER:  nattodo=',nattodo
c     call f_lush(6)
ckw
ckw   if(nate.eq.na) then
      if(nate.eq.nattodo) then
        natbe=-1
        natee=-1
cc        write(6,*) ' MASTER:  Sending CPHF Termination'
cc        call f_lush(6)
        call para_initsend
        call para_pack_int(natbe,1)
        call para_pack_int(natee,1)
        call para_bcast_pack(TxHessDat)
      endif
c
      end
c=======================================================================
c
      subroutine para_d0g2(idft,ax,rhf,nblocks,bl,inx,ntri,threg2,
     *                     densp, ncenter,
     *                     dens,denB,hess,labels)

      use newpara

      implicit real*8 (a-h,o-z)
c
c This routine calculates integral 2nd derivatives (for analytical Hessian)
c block by block and then calculates their contribution to the Hessian.
c Actually, this routine handles the PVM calls and distributes tasks
c among nodes. 
c---------------------------------------------------------------------
c integral 2nd derivatives are calculated in int_d0g2()
c---------------------------------------------------------------------
c INPUT :
c  idft              - dft flag:  0 - no dft
c  ax                - factor to multiply exchange contribution for ACM
c  rhf               - logical flag for closed-shell
c  nblocks           - number of blocks of contracted shell quartets
c  bl                - scratch for everything
c  inx(12,ncs)       - basis set & nuclear data
c  ntri              - for dimensioning
c  threg2            - integral threshold
c  densp             - screening density
c  ncenter           - center of each basis function
c  dens              - closed-shell or alpha density matrix
c  denB              - beta density matrix
c  hess              - on exit contains Hessian contribution
c  labels            - labels buffer
c---------------------------------------------------------------------
      character*11 scftype
      character*4 where
      Logical rhf
      common /runtype/ scftype,where
      COMMON /DFTCoeff/ aXX(18) ! contains 15 density functional coefficients
      dimension inx(12,*)
      dimension bl(*),dens(*),denB(*),hess(*)
      dimension densp(*)
      dimension ncenter(*)
      dimension labels(*)
c
      call getival('ncf',ncf)
      call getival('ncs',ncs)
c
      call getival('gran',igran)
      call getival('na  ',na)
      nat3=na*3
c
      if(scftype.ne.'full-direct') then
         call nerror(1,'para_d0g2',
     *                  'can only run in full-direct mode',0,0)
      endif
      call para_bcast_real(dens,ntri,TxHess1)
      If(.NOT.rhf) call para_bcast_real(denB,ntri,TxHess1)
      call para_initsend
      call para_pack_real(densp,ncs*ncs)
      call para_pack_int(ncenter,ncf)
      call para_pack_int(idft,1)
      call para_pack_real(aXX,18)
      call para_bcast_pack(TxHess2)
c
c Give work assignments to slaves until all work is done
      call para_distr_blocks(nblocks,igran)
c
c do global sum of results
c
      call para_reduce(hess,nat3**2,TxHess1)
      end
c ..........................................................................
c
      subroutine para_cphf_nat(natoms, natend, natonce,npos,   NQ,
     *                         IUNQ,   idft,   ax,     nblocks,bl,
     *                         inx,    ntri,   thres,  lsemi,  dens,
     *                         fock,   labels)

      use memory, block => bl
      use newpara
      
      implicit real*8 (a-h,o-z)
c
c Three Fock matrices are constructed F(DX,g0),F(DY,g0) and F(DZ,g0)
c for natonce atoms using first-order densities and standard
c (zeroth-order) integrals
c---------------------------------------------------------------------
c integrals are calculated in int_d1g0()
c---------------------------------------------------------------------
c INPUT :
c  natoms         -  total number of atoms
c  natend         -  shows if cphf for each of the atoms treated
c                      is converged (1) or not (0)
c  natonce        -  number of atoms treated at once in cphf
c  npos           -  starting position of current first atom
c  NQ             -  number of symmetry-unique atoms
c  IUNQ           -  list of symmetry-unique atoms
c  idft           -  dft flag:  0 - no dft
c  ax             -  factor to multiply exchange contribution for ACM
c  nblocks        -  number of blocks of contracted shell quartets
c  bl             -  scratch for everything
c  inx(12,ncs)    -  basis set & nuclear data
c  ntri           -  for dimensioning
c  thres          -  integral threshold
c  lsemi          -  flag for reuse of DFT grid
c                    (currently semidirect is NOT an option is parallel)
c  dens           -  first-order closed-shell density matrices
c  fock           -  on exit contains CPHF Fock matrices
c  labels         -  labels buffer
c---------------------------------------------------------------------
      character*11 scftype
      character*4 where,screen
      common /screen_type/ screen
cc      Logical rhf
      common /runtype/ scftype,where
      common /lindvec/ lind,idensp
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      dimension inx(12,*)
      dimension bl(*),dens(ntri,3,natonce),fock(ntri,3,natonce)
      dimension labels(*),natend(natonce)
c
c
c -- prepare for integral recalculation
c----------------------------------------------------------------
c transpose dens from (ntri,3,natonce) into (3,natonce,ntri)
c
      call trsp_inplace(bl,dens,ntri,3*natonce)
c----------------------------------------------------------------
c allocate memory for denspar(ics,jcs) (screening density)
c
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,map_fs)      ! mapping from funct. to shells
c----------------------------------------------------------------
c----------------------------------------------------------------
c------------ screening in CPHF ---------------------------------
c
c
      call getival('delta',idelt)
c
      if(idelt.eq.0) then   ! no delta density in cphf screening
c        in this case d1const is used as a screening density.
c        Screening den. over shells is done in get_d1w1.f
c        and written on unit 69 in the second record
c
         call read1mat(69, 2 ,ncs*ncs,bl(idensp))
c
         if(screen.eq.'2pde') then
c          2-particle type screening :
c          read in the second screening density constructed
c          before in hessana
           call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c
c
c there are 3 density matrices for NATONCE atoms
c select one (max) for screening in calcint2
c
      if(idelt.ne.0) then   ! delta density in cphf screening
         call getmem(ntri,lselect)
         call select_dnat(natend,natonce,dens,ntri,bl(lselect))
         call setup_densp2(inx,ncf,ncs,bl(lselect),
     *                        bl(idensp), bl(map_fs))
         call retmem(1)    ! lselect
         if(screen.EQ.'2pde') then
c           read-in second screening density
            call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c----------------------------------------------------------------
      call getival('gran',igran)
c
      call zeroit(fock,ntri*3*natonce)
c
c -- send CPHF data to slaves
ckw     call elapsec(ttt1)
      call para_bcast(natonce,TxHessDat)
      call para_initsend
      call para_pack_int(natend,natonce)
      call para_pack_string(screen,4)
      call para_pack_real(thres,1)
      call para_pack_int(bl(map_fs),ncf)
      call para_pack_real(bl(idensp),2*ncs*ncs)
      call para_bcast_pack(TxHessDat)
      call para_bcast_real(dens,3*ntri*natonce,TxHessDat)
ckw     call elapsec(ttt2)
ckw     ttt1 = (ttt2-ttt1)/60.0d0
ckw     write(6,*) ' CPHF: Time sending data to slaves:',ttt1,' min.'
c
c Give work assignments to slaves until all work is done
      call para_distr_blocks(nblocks,igran)
c----------------------------------------------------------------
c     call retint(1)    ! map_fs
      call retmem(2)    ! map_fs, idensp
c----------------------------------------------------------------
c
c re-scaling is done on the slaves, results are only summed at
c the end, no more need to transpose      
c
c----------------------------------------------------------------
c     do iat=1,natonce
c        call rescale_fock(fock(1,1,iat),ntri,bl(lind),thres)
c     enddo
c----------------------------------------------------------------
c
c ..............................................................
c  DFT CONTRIBUTIONS
c
c  Calculate DFT XC contributions to partial derivative Fock
c  matrices for CPHF equations
c
      IF(idft.GT.0) THEN
        call getival('lcore',lcore)
        call getmem(0,lastx)
        call retmem(1)
        CALL CPHFDFT(natoms, natonce,natend, npos,    NQ,
     $               IUNQ,   dens,   dens,   fock,    fock,
     $               thres,  lsemi,lcore-lastx,bl(lastx))
      ENDIF
c----------------------------------------------------------------
c
c do global sum of results
c
      call para_reduce(fock,3*ntri*natonce,TxHess3)
c ..............................................................
c transpose back fock1 & d1 from (3,natonce,ntri) into (ntri,3,natonce)
c
      call trsp_inplace(bl,fock,3*natonce,ntri)
      call trsp_inplace(bl,dens,3*natonce,ntri)
c ..............................................................
c
      return
      end
c ..........................................................................
c
      subroutine para_cphf_nat_uhf(natoms, natend, natonce,npos,   NQ,
     *                             IUNQ,   idft,   ax,     nblocks,bl,
     *                             inx,    ntri,   thres,  lsemi,  denA,
     *                             denB,   fockA,  fockB,  labels)

      use memory, block => bl
      use newpara
      
      implicit real*8 (a-h,o-z)
c
c Three Fock matrices are constructed F(DX,g0),F(DY,g0) and F(DZ,g0)
c for natonce atoms using first-order densities and standard
c (zeroth-order) integrals
c---------------------------------------------------------------------
c integrals are calculated in int_d1g0()
c---------------------------------------------------------------------
c INPUT :
c  natoms         -  total number of atoms
c  natend         -  shows if cphf for each of the atoms treated
c                      is converged (1) or not (0)
c  natonce        -  number of atoms treated at once in cphf
c  npos           -  starting position of current first atom
c  NQ             -  number of symmetry-unique atoms
c  IUNQ           -  list of symmetry-unique atoms
c  idft           -  dft flag:  0 - no dft
c  ax             -  factor to multiply exchange contribution for ACM
c  nblocks        -  number of blocks of contracted shell quartets
c  bl             -  scratch for everything
c  inx(12,ncs)    -  basis set & nuclear data
c  ntri           -  for dimensioning
c  thres          -  integral threshold
c  lsemi          -  flag for reuse of DFT grid
c                    (currently semidirect is NOT an option is parallel)
c  denA           -  first-order alpha density matrices
c  denB           -  first-order beta  density matrices
c  fockA          -  on exit contains alpha CPHF Fock matrices
c  fockB          -  on exit contains beta  CPHF Fock matrices
c  labels         -  labels buffer
c---------------------------------------------------------------------
      character*11 scftype
      character*4 where,screen
      common /screen_type/ screen
      common /runtype/ scftype,where
      common /lindvec/ lind,idensp
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      dimension inx(12,*)
      dimension denA(ntri,3,natonce),fockA(ntri,3,natonce)
      dimension denB(ntri,3,natonce),fockB(ntri,3,natonce)
      dimension bl(*),labels(*),natend(natonce)
c
c
c -- prepare for integral recalculation
c----------------------------------------------------------------
c transpose dens from (ntri,3,natonce) into (3,natonce,ntri)
c
      call trsp_inplace(bl,denA,ntri,3*natonce)
      call trsp_inplace(bl,denB,ntri,3*natonce)
c----------------------------------------------------------------
c allocate memory for denspar(ics,jcs) (screening density)
c
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,map_fs)      ! mapping from funct. to shells
c----------------------------------------------------------------
c------------ screening in CPHF ---------------------------------
c
      call getival('delta',idelt)
c
      if(idelt.eq.0) then   ! no delta density in cphf screening
c        in this case d1const is used as a screening density.
c        Screening den. over shells is done in get_d1w1.f
c        and written on unit 69 in the second record
c
         call read1mat(69, 2 ,ncs*ncs,bl(idensp))
c
         if(screen.eq.'2pde') then
c          2-particle type screening :
c          read in the second screening density constructed
c          before in hessana
           call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c
c there are 3 density matrices for NATONCE atoms
c select one (max) for screening in calcint2
c
      if(idelt.ne.0) then   ! delta density in cphf screening
         call getmem(ntri,lselect)
      call select_dnat_uhf(natend,natonce,denA,denB,ntri,bl(lselect))
         call setup_densp2(inx,ncf,ncs,bl(lselect),
     *                        bl(idensp), bl(map_fs))
c                    screening Dens NOT squared
         call retmem(1)    ! lselect
         if(screen.EQ.'2pde') then
c           read-in second screening density
            call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c----------------------------------------------------------------
      call getival('gran',igran)
c
      call zeroit(fockA,ntri*3*natonce)
      call zeroit(fockB,ntri*3*natonce)
c
c -- send CPHF data to slaves
ckw     call elapsec(ttt1)
        call para_bcast(natonce,TxHessDat)
        call para_initsend
        call para_pack_int(natend,natonce)
        call para_pack_string(screen,4)
        call para_pack_real(thres,1)
        call para_pack_int(bl(map_fs),ncf)
        call para_pack_real(bl(idensp),2*ncs*ncs)
        call para_bcast_pack(TxHessDat)
        call para_bcast_real(denA,3*ntri*natonce,TxHessDat)
        call para_bcast_real(denB,3*ntri*natonce,TxHessDat)
ckw     call elapsec(ttt2)
ckw     ttt1 = (ttt2-ttt1)/60.0d0
ckw     write(6,*) ' CPHF: Time sending data to slaves:',ttt1,' min.'
c
c Give work assignments to slaves until all work is done
        call para_distr_blocks(nblocks,igran)
c
c ---------------------------------------------------------------------
c     call retint(1)    ! map_fs
      call retmem(2)    ! map_fs,idensp
c----------------------------------------------------------------
c
c re-scaling is done on the slaves, results are only summed at
c the end, no more need to transpose      
c
c ..............................................................
c  DFT CONTRIBUTIONS
c
c  Calculate DFT XC contributions to partial derivative Fock
c  matrices for CPHF equations
c
      IF(idft.GT.0) THEN
        call getival('lcore',lcore)
        call getmem(0,lastx)
        call retmem(1)
        CALL CPHFDFT(natoms, natonce,natend, npos,    NQ,
     $               IUNQ,   denA,   denB,   fockA,   fockB,
     $               thres,  lsemi,lcore-lastx,bl(lastx))
      ENDIF
      call para_reduce(fockA,3*ntri*natonce,TxHess3)
c -- now send beta Fock derivative matrices
      call para_reduce(fockB,3*ntri*natonce,TxHess2)
c----------------------------------------------------------------
c transpose back fock1 & d1 from (3,natonce,ntri) into (ntri,3,natonce)
c
      call trsp_inplace(bl,fockA,3*natonce,ntri)
      call trsp_inplace(bl,fockB,3*natonce,ntri)
      call trsp_inplace(bl,denA,3*natonce,ntri)
      call trsp_inplace(bl,denB,3*natonce,ntri)
c ..............................................................
c
      return
      end
c ..........................................................................
c
c for polarizability CPHF
c
      subroutine para_cphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres,dens,
     *                         fock,labels)
      implicit real*8 (a-h,o-z)
      call int_d1g0(idft,ax,nblocks,bl,inx,ntri,thres,dens,fock,
     *              labels,0,0)
      end
c ..........................................................................
c
      subroutine para_d1const(natoms, natb,   nate,   listunq,ncf,
     $                        ntri,   nocc,   lind,   dens0,  s1,
     $                        f1,     vec,    val,    bl,     fds1,
     $                        r1)

      use memory, block => bl
      use newpara

      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(ncf),dens0(*),s1(*),r1(*),vec(*),val(*),
     $          bl(*),fds1(ncf,ncf,3),f1(ntri,3,natoms)
      data nfile60/60/, nfile61/61/, nfile62/62/, nfile65/65/
c
c ** This should be parallelized differently - too much communication **
      ncf2 = ncf**2
      ntri3 = 3*ntri
      natodo = nate-natb+1          ! number of atoms to do
c
c -- for VCD reserve memory for constant part of perturbed
c -- wavefunction coefficients
      call tstival('vcd',ivcd)
      If(ivcd.NE.0) Then
        call getmem(ncf*nocc*3,lpve)
        ncoefile = 79    ! unit for the pert. coeff. file
      EndIf
c
      If(nslv.GT.natodo) Then
c -- inefficiency warning
        call message('**WARNING** in Hessian CPHF',
     $  'Not possible to use slaves efficiently: Nslv>Nat:',nslv,natodo)
      EndIf
c
c -- send to each slave the number of atoms to do and the VCD parameter
c   (needed for program flow)
      call para_bcast(natodo,TxHessDat)
      call para_bcast(ivcd,TxHessDat)
c
c -- read in precomputed FD matrix (use r1 temporarily as scratch)
      lrec=1
      call read1mat(nfile60,lrec,ncf2,r1)
c
c -- now broadcast to each slave the constant data
      call para_bcast_real(dens0,ntri,TxWDens1)
      call para_bcast_real(r1,ncf2,TxWDens1)
      call para_bcast_real(vec,ncf2,TxWDens1)
      call para_bcast_real(val,ncf,TxWDens1)
c
      call elapsec(timo)
      call getrval('imbal',tlost)      
      nidl=nslv
c -- now send first-order data for each atom
      Do iat=natb,nate+nslv
         I=iat-natb+1 
         nat = listunq(iat)
c -- get request and results if any         
         call para_recv(iatm,islave,TxWdens2)
         if (iatm.gt.0) then
            call para_recv_real(r1,ntri3,islave,TxWdens2)
            call save1mat(nfile65,iatm,ntri3,r1)
            call para_recv_real(fds1,ncf*ncf*3,islave,TxWdens2)
c -- just got constant part of density and FDS1 - write to disk
            irec=(iatm-1)*3 + lrec
            call save1mat(nfile60,irec+1,ncf2,fds1(1,1,1) )
            call save1mat(nfile60,irec+2,ncf2,fds1(1,1,2) )
            call save1mat(nfile60,irec+3,ncf2,fds1(1,1,3) )
            If(ivcd.NE.0) Then
              call para_recv_real(bl(lpve),ncf*nocc*3,islave,TxWdens2)
              call save1mat(ncoefile,iatm,ncf*nocc*3,bl(lpve))
            EndIf
         else
c -- account for idling slaves that did not get a job yet 
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl-1
            timo=timb
         endif
         if (iat.le.nate) then   !not all done
c -- read in fock1 and s1 matrices from disk
            call read1mat(nfile61,nat,ntri3,f1)
            call read1mat(nfile62,nat,ntri3,s1)
            call para_send(I,islave,TxWDens3)  ! send the atom
            call para_send_real(f1,ntri3,islave,TxWDens3)
            call para_send_real(s1,ntri3,islave,TxWDens3)
         else
            call para_send(-1,islave,TxWDens3)  ! finished
c -- account for slaves that finished early
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl+1
            timo=timb
         endif
      enddo
      call setrval('imbal',tlost)
c
      If(ivcd.NE.0) call retmem(1)
c
      return
      end
c ..........................................................................
c
      subroutine para_d1const_uhf(natoms, natb,   nate,   listunq,ncf,
     $                            ntri,   nalpha, nbeta,  lind,  dens0A,
     $                            dens0B, s1,     f1,     vecA,   vecB,
     $                            valA,   valB,   bl,     fds1A,  fds1B,
     $                            r1)

      use newpara

      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0A(*),dens0B(*),s1(*),f1(*),
     $          vecA(*),vecB(*),valA(*),valB(*),bl(*),
     $          fds1A(ncf,ncf,3),fds1B(ncf,ncf,3),r1(*)
      data nfile60/60/, nfile61/61/, nfile62/62/, nfile65/65/
      data nfile71/71/, nfile75/75/
c
c
      ncf2 = ncf**2
      ntri3 = 3*ntri
      natodo = nate-natb+1          ! number of atoms to do
c
      If(nslv.GT.natodo) Then
c -- inefficiency warning
        call message('**WARNING** in Hessian CPHF',
     $  'Not possible to use slaves efficiently: Nslv>Nat:',nslv,natodo)
      EndIf
c
c -- send to each slave the number of atoms to do
c   (needed for program flow)
      call para_bcast(natodo,TxHessDat)
      ivcd0=0
      call para_bcast(ivcd0,TxHessDat)      ! VCD flag, not in UHF
c
c -- read in precomputed FD matrix (use r1/f1 temporarily as scratch)
      lrec=1
      call read1mat(nfile60,lrec,ncf2,r1)
      lrec=natodo*3+2
      call read1mat(nfile60,lrec,ncf2,f1)
      call dscal(ncf2, 2.d0 ,r1  ,1)
      call dscal(ncf2, 2.d0 ,f1  ,1)
c
c -- now broadcast to each slave the constant data
      call para_bcast_real(dens0A,ntri,TxWDens1)
      call para_bcast_real(r1,ncf2,TxWDens1)
      call para_bcast_real(vecA,ncf2,TxWDens1)
      call para_bcast_real(valA,ncf,TxWDens1)
      call para_bcast_real(dens0B,ntri,TxWDens1)
      call para_bcast_real(f1,ncf2,TxWDens1)
      call para_bcast_real(vecB,ncf2,TxWDens1)
      call para_bcast_real(valB,ncf,TxWDens1)
c
      call elapsec(timo)
      call getrval('imbal',tlost)      
      nidl=nslv
c -- now send first-order data for each atom
      Do iat=natb,nate+nslv
         I=iat-natb+1 
         nat = listunq(iat)
c -- get request and results if any         
         call para_recv(iatm,islave,TxWdens2)
         if (iatm.gt.0) then
c -- alpha data
            call para_recv_real(r1,ntri3,islave,TxWdens2)
            call save1mat(nfile65,iatm,ntri3,r1)
            call para_recv_real(fds1A,ncf*ncf*3,islave,TxWdens2)
c -- just got constant part of density and FDS1 - write to disk
            lrec=1
            irec=(iatm-1)*3 + lrec
            call save1mat(nfile60,irec+1,ncf2,fds1A(1,1,1) )
            call save1mat(nfile60,irec+2,ncf2,fds1A(1,1,2) )
            call save1mat(nfile60,irec+3,ncf2,fds1A(1,1,3) )
c -- beta data
            call para_recv_real(f1,ntri3,islave,TxWdens2)
            call save1mat(nfile75,iatm,ntri3,f1) ! D1const beta
            call para_recv_real(fds1B,ncf*ncf*3,islave,TxWdens2)
c -- just got constant part of density and FDS1 - write to disk
            lrec=natodo*3+2
            irec=(iatm-1)*3 + lrec
            call save1mat(nfile60,irec+1,ncf2,fds1B(1,1,1) )
            call save1mat(nfile60,irec+2,ncf2,fds1B(1,1,2) )
            call save1mat(nfile60,irec+3,ncf2,fds1B(1,1,3) )
         else
c -- account for idling slaves that did not get a job yet 
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl-1
            timo=timb
         endif
         if (iat.le.nate) then   !not all done
c -- read in fock1 and s1 matrices from disk
            call read1mat(nfile61,nat,ntri3,f1)
            call read1mat(nfile62,nat,ntri3,s1)
            call read1mat(nfile71,nat,ntri3,r1)
            call para_send(I,islave,TxWDens3)  ! send the atom
            call para_send_real(f1,ntri3,islave,TxWDens3)
            call para_send_real(s1,ntri3,islave,TxWDens3)
            call para_send_real(r1,ntri3,islave,TxWDens3)
         else
            call para_send(-1,islave,TxWDens3)  ! finished
c -- account for slaves that finished early
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl+1
            timo=timb
         endif
      enddo
      call setrval('imbal',tlost)
C
      return
      end
c ..........................................................................
c
      subroutine para_wdens(natb,   nate,   listunq,ncf,    ntri,
     $                      lind,   dens0,  fd0,    f1,     gmat1,
     $                      d1,     work,   r1)

      use newpara

      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(ncf),dens0(*),fd0(*),work(ncf,ncf,3),
     $          f1(ntri,3),gmat1(ntri,3,*),d1(ntri,3,*),r1(ntri,3)
      data nfile61/61/, nfile63/63/, nfile64/64/
c
      ntri3 = 3*ntri
      natodo = nate-natb+1          ! number of atoms to do
c
      If(nslv.GT.natodo) Then
c -- inefficiency warning
        call message('**WARNING** in Hessian CPHF',
     $  'Not possible to use slaves efficiently: Nslv>Nat:',nslv,natodo)
c -- we are going to temporarily reduce the number of active slaves
      EndIf
c
c -- first broadcast to each slave the constant data
      call para_initsend
      call para_pack_real(dens0,ntri)
      call para_pack_real(fd0,ncf**2)
      call para_bcast_pack(TxWDens1)
c
      call elapsec(timo)
      call getrval('imbal',tlost)      
      nidl=nslv
c -- now send first-order data for each atom
      Do iat=natb,nate+nslv
c -- get request and results if any         
         call para_recv_pack(islave,TxWDens2)
         call para_unpack_int(J,1)       ! get back the actual array entry
         call para_unpack_int(nat,1)     !  and the actual atom
         if (J.gt.0) then
            call para_recv_real(r1,ntri3,islave,TxWDens4)
c -- just got weighted density - write it and d1 to disk
            call save1mat(nfile63,nat,ntri3,d1(1,1,J) )
            call save1mat(nfile64,nat,ntri3,r1)
         else
c -- account for idling slaves that did not get a job yet 
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl-1
            timo=timb
         endif
         I=iat-natb+1 
         nat = listunq(iat)
         if (iat.le.nate) then   !not all done
c -- read in fock1 matrix from disk
            call read1mat(nfile61,nat,ntri3,f1)
            call para_send(I,islave,TxWDens3) ! send the actual array entry
            call para_initsend
            call para_pack_int(nat,1)       !  and the actual atom
            call para_pack_real(f1,ntri3)
            call para_pack_real(gmat1(1,1,I),ntri3)
            call para_pack_real(d1(1,1,I),ntri3)
            call para_send_pack(islave,TxWDens5)
         else
            call para_send(-1,islave,TxWDens3)  ! finished
c -- account for slaves that finished early
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl+1
            timo=timb
         endif
      EndDo
      call setrval('imbal',tlost)
C
      return
      end
c ..........................................................................
c
      subroutine para_wdensu_xyz(natb,   nate,   listunq,ncf,    ntri,
     $                           lind,   dens0A, dens0B, fd0A,   fd0B,
     $                           f1,     gmat1A, gmat1B, d1A,    d1B,
     $                           work,   s1,     r1)

      use newpara

      implicit real*8(a-h,o-z)
      dimension listunq(*),lind(*),dens0A(*),dens0B(*),fd0A(*),
     $          fd0B(*),f1(ntri,3),gmat1A(ntri,3,*),gmat1B(ntri,3,*),
     $          d1A(ntri,3,*),d1B(ntri,3,*),work(ncf,ncf,3),
     $          s1(ntri*3),r1(ntri,3)
      data nfile61/61/, nfile63/63/, nfile64/64/,
     $     nfile71/71/, nfile73/73/
c
      ntri3 = 3*ntri
      natodo = nate-natb+1          ! number of atoms to do
c
      If(nslv.GT.natodo) Then
c -- inefficiency warning
        call message('**WARNING** in Hessian CPHF',
     $  'Not possible to use slaves efficiently: Nslv>Nat:',nslv,natodo)
      EndIf
c
c -- first broadcast to each slave the constant data
      call para_initsend
      call para_pack_real(dens0A,ntri)
      call para_pack_real(fd0A,ncf**2)
      call para_pack_real(dens0B,ntri)
      call para_pack_real(fd0B,ncf**2)
      call para_bcast_pack(TxWDens1)
c
      call elapsec(timo)
      call getrval('imbal',tlost)      
      nidl=nslv
c
c -- now send each slave first-order data for each atom
      Do iat=natb,nate+nslv
c -- get request and results if any         
         call para_recv_pack(islave,TxWdens2)
         call para_unpack_int(J,1)       ! get back the actual array entry
         call para_unpack_int(nat,1)     !  and the actual atom
         if (J.gt.0) then
            call para_recv_real(r1,ntri3,islave,TxWdens4)
c -- just got weighted density - write it and d1 to disk
            call save1mat(nfile64,nat,ntri3,r1)
            call save1mat(nfile63,nat,ntri3,d1A(1,1,J) )
            call save1mat(nfile73,nat,ntri3,d1B(1,1,J) )
         else
c -- account for idling slaves that did not get a job yet 
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl-1
            timo=timb
         endif
         I=iat-natb+1 
         nat = listunq(iat)
         if (iat.le.nate) then   !not all done
c -- read in fock1 matrix from disk
            call read1mat(nfile61,nat,ntri3,f1)
            call read1mat(nfile71,nat,ntri3,s1)
            call para_send(I,islave,TxWDens3) ! send the actual array entry
            call para_initsend
            call para_pack_int(nat,1)       !  and the actual atom
            call para_pack_real(f1,ntri3)
            call para_pack_real(gmat1A(1,1,I),ntri3)
            call para_pack_real(d1A(1,1,I),ntri3)
            call para_pack_real(s1,ntri3)
            call para_pack_real(gmat1B(1,1,I),ntri3)
            call para_pack_real(d1B(1,1,I),ntri3)
            call para_send_pack(islave,TxWDens5)
         else
            call para_send(-1,islave,TxWDens3)  ! finished
c -- account for slaves that finished early
            call elapsec(timb)
            tlost=tlost+nidl*(timb-timo)
            nidl=nidl+1
            timo=timb
         endif
      EndDo
      call setrval('imbal',tlost)
c
      return
      end

      subroutine para_cphf_nat(natoms, natend, natonce,npos,   NQ,
     *                         IUNQ,   idft,   ax,     nblocks,bl,
     *                         inx,    ntri,   thres,  lsemi,  dens,
     *                         fock,   labels)
c
c=====================================================================
c This routine calls two-el. int. block by block and constructs
c closed-shell Fock matrix using 1st-order density.
c---------------------------------------------------------------------
c This will construct THREE G(D1,g0) matrices for NATONCE atoms !!!
c---------------------------------------------------------------------
      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*4 screen
      common /screen_type/ screen
      common /lindvec/ lind,idensp
      common /datstore/ thresh_stored,isto
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      dimension inx(12,*)
      dimension bl(*)
      dimension dens(ntri,3,natonce),fock(ntri,3,natonce)
      dimension natend(natonce)
      dimension iunq(nq)
      dimension labels(*)
c----------------------------------------------------------------
c
      call zeroit(fock,ntri*3*natonce)
c
c setup flag for retrieval of stored integrals
c
      isto=1
c
c prepare for integral recalculation
c----------------------------------------------------------------
c transpose d1 from (ntri,3,natonce) into (3,natonce,ntri)
c
      call trsp_inplace(bl,dens,ntri,3*natonce)
c----------------------------------------------------------------
c allocate memory for denspar(ics,jcs) (screening density)
c
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,map_fs)      ! mapping from funct. to shells
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
c
c       we also need to set upthe the mapping (funcions->shell) array
c
        call setup_mapfs(inx,ncf,ncs,bl(map_fs))
c
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
c                    screening Dens NOT squared
         call retmem(1)    ! lselect
c
         if(screen.EQ.'2pde') then
c           read-in second screening density
            call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c----------------------------------------------------------------
c
      call int_d1g0nat(
     $        natend,      natonce,    idft,       ax,
     $        nblocks,     bl,         inx,        ntri,
     $        thres,       bl(map_fs), bl(idensp), dens,
     $        fock,        labels,     0,          0)
c
c----------------------------------------------------------------
      call retmem(1)    ! map_fs
      call retmem(1)    ! idensp
c----------------------------------------------------------------
c re-scale fock matrices in transposed form:
c
      call rescale_nat_fock(fock,ntri,natonce,bl(lind),thres*0.5d0)
c---------------------------------------------------------------------
c symmetrization of G(D1,g0) is done in the fock_d1g0 builder.
c While integrals g0 keep symmetry the D1 does not . That is why
c we need to do symmetrization explicitly .
c---------------------------------------------------------------------
c
c  DFT CONTRIBUTIONS
c
c  Calculate DFT XC contributions to partial derivative Fock
c  matrices for CPHF equations
c
      IF(idft.GT.0) THEN
        call getival('lcore',lcore)
        call getmem(0,lastx)
        call retmem(1)
        CALL CPHFDFT(natoms, natonce,natend,npos,   NQ,
     $               IUNQ,   dens,   dens,  fock,   fock,
     $               thres,  lsemi,lcore-lastx,bl(lastx))
      ENDIF
c no more transposing while calculating      
      call trsp_inplace(bl,fock,3*natonce,ntri)
      call trsp_inplace(bl,dens,3*natonce,ntri)
c ..............................................................
c
      return
      end
c ========================================================================
c
      subroutine para_cphf_nat_uhf(natoms, natend, natonce,npos,   NQ,
     *                             IUNQ,   idft,   ax,     nblocks,bl,
     *                             inx,    ntri,   thres,  lsemi, densA,
     *                             densB,  fockA,  fockB,  labels)
c
c=====================================================================
c This routine calls two-el. int. block by block and constructs
c open-shell Fock matrices using 1st-order density.
c---------------------------------------------------------------------
c This will construct THREE G(D1,g0) matrices for NATONCE atoms !!!
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*4 screen
      common /screen_type/ screen
      common /lindvec/ lind,idensp
      common /datstore/ thresh_stored,isto
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      dimension inx(12,*)
      dimension bl(*)
      dimension densA(ntri,3,natonce),fockA(ntri,3,natonce)
      dimension densb(ntri,3,natonce),fockB(ntri,3,natonce)
      dimension natend(natonce)
      dimension iunq(nq)
      dimension labels(*)
c----------------------------------------------------------------
c
      call zeroit(fockA,ntri*3*natonce)
      call zeroit(fockB,ntri*3*natonce)
c
c setup flag for retrieval of stored integrals
c
      isto=1
c
c prepare for integral recalculation
c----------------------------------------------------------------
c transpose d1 from (ntri,3,natonce) into (3,natonce,ntri)
c
      call trsp_inplace(bl,densA,ntri,3*natonce)
      call trsp_inplace(bl,densB,ntri,3*natonce)
c----------------------------------------------------------------
c allocate memory for denspar(ics,jcs) (screening density)
c
      call getmem(2*ncs*ncs,idensp)
      call getint(ncf,map_fs)      ! mapping from funct. to shells
c----------------------------------------------------------------
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
c
c       we also need to set upthe the mapping (funcions->shell) array
c
        call setup_mapfs(inx,ncf,ncs,bl(map_fs))
c
      endif
c
c
c there are 3 density matrices for NATONCE atoms
c select one (max) for screening in calcint2
c
      if(idelt.ne.0) then   ! delta density in cphf screening
         call getmem(ntri,lselect)
      call select_dnat_uhf(natend,natonce,densA,densB,ntri,bl(lselect))
         call setup_densp2(inx,ncf,ncs,bl(lselect),
     *                     bl(idensp), bl(map_fs))
c                    screening Dens NOT squared
         call retmem(1)    ! lselect
         if(screen.EQ.'2pde') then
c           read-in second screening density
            call read1mat(69, 1 ,ncs*ncs,bl(ncs*ncs+idensp))
         endif
      endif
c----------------------------------------------------------------
c
      call int_d1g0nat_uhf(
     $        natend,      natonce,    idft,       ax,
     $        nblocks,     bl,         inx,        ntri,
     $        thres,       bl(map_fs), bl(idensp), densA,
     $        densB,       fockA,      fockB,      labels,
     $        0,           0)
c
c----------------------------------------------------------------
      call retmem(1)    ! map_fs
      call retmem(1)    ! idensp
c----------------------------------------------------------------
c re-scale fock matrices in transposed form
c
      call rescale_nat_fock(fockA,ntri,natonce,bl(lind),thres)
      call rescale_nat_fock(fockB,ntri,natonce,bl(lind),thres)
c---------------------------------------------------------------------
c symmetrization of G(D1,g0) is done in the fock_d1g0 builder.
c While integrals g0 keep symmetry the D1 does not . That is why
c we need to do symmetrization explicitly .
c---------------------------------------------------------------------
c
c  DFT CONTRIBUTIONS
c
c  Calculate DFT XC contributions to partial derivative Fock
c  matrices for CPHF equations
c
      IF(idft.GT.0) THEN
        call getival('lcore',lcore)
        call getmem(0,lastx)
        call retmem(1)
        CALL CPHFDFT(natoms, natonce,natend, npos,   NQ,
     $               IUNQ,   densA,  densB,  fockA,  fockB,
     $               thres,  lsemi,lcore-lastx,bl(lastx))
      ENDIF
      call trsp_inplace(bl,fockA,3*natonce,ntri)
      call trsp_inplace(bl,fockB,3*natonce,ntri)
      call trsp_inplace(bl,densA,3*natonce,ntri)
      call trsp_inplace(bl,densB,3*natonce,ntri)
c ..............................................................
c
      return
      end
c=====================================================================
      subroutine setup_mapfs(inx,ncf,ncs,map_fs)
      implicit real*8 (a-h,o-z)
      dimension inx(12,ncs)
      dimension map_fs(ncf)
c--------------------------------------------------------------------
c densmat is assumed to be ncf*(ncf+1)/2 in all cases
c--------------------------------------------------------------------
c input :
c inx(12,ncs)       = basis set info
c ncf, ncs          = contracted functions and shells
c
c output :
c map_fs(icf)=ics     mapping from contr.funct to contr.shells
c--------------------------------------------------------------------
c
      do 10 ics=1,ncs
      icf_b=inx(11,ics)+1
      icf_e=inx(10,ics)
      do 10 icf=icf_b,icf_e
      map_fs(icf)=ics
   10 continue
c
      end
c=====================================================================

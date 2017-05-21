c=====================================================================
      subroutine get_d1xyz(idft,ax,rhf,bl,inx,ncf,ntri,sxyz,fxyz,dxyz,
     *                     lsemi)
c--------------------------------------------------------------------
c calculates Ist-order density for electric filed perturbation (polar)
c--------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
c----------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension sxyz(ntri,3)  ! input Ist-order overlap
      dimension fxyz(ntri,3)  ! input F1=H1 + G(D0,g1)
      dimension dxyz(ntri,3)  ! output Ist-order density
c----------------------------------------------------
      call getival('iout',iout)
c----------------------------------------------------
c     call secund(tchfb)
c     call elapsec(etchfb)
c----------------------------------------------------
c get CPHF parameters :
c
      call getrval('thres',cphf_thresh) ! cphf-threshold
c
c print the info above into the output file :
c
      call cphf_infoxyz(cphf_thresh)
c----------------------------------------------------
      ntri3=ntri*3
c----------------------------------------------------
c check bl() status and memory needed & available :
c
      call cphf_xyzcheck
c----------------------------------------------------
c get zeroth-order density pointer:
c
      call getival('ldensi',lden0)
c
      call getival('lvec',lvec)      ! eigenvectors
      call getival('lval',lval)      ! eigenvalues
      call getival('nocc',nocc)      ! no of occupied orbitals
ctest
c     write(6,*) ' 0th-order density matrix '
c     call drumh(bl(lden0),ncf, 6  ,'density0')
c
c     write(6,*) ' Fxyz marices             '
c     call drumh(fxyz(1,1),ncf, 6  ,'Fx = Hx ')
c     call drumh(fxyz(1,2),ncf, 6  ,'Fy = Hy ')
c     call drumh(fxyz(1,3),ncf, 6  ,'Fz = Hz ')
c
c     write(6,*) ' Sxyz marices             '
c     call drumh(sxyz(1,1),ncf, 6  ,'Sx = 0  ')
c     call drumh(sxyz(1,2),ncf, 6  ,'Sy = 0  ')
c     call drumh(sxyz(1,3),ncf, 6  ,'Sz = 0  ')
c----------------------------------------------------
c put memory mark
c
      call mmark
c----------------------------------------------------
c allocate local memory :
c
      call getmem(ntri3  ,ld1c)      ! constant part of D1
c----------------------------------------------------
c   DO CPHF ;  get D1
c
      call chfsol_xyz1(idft,ax,bl,inx,nocc,ncf,ntri,cphf_thresh,
     *                 bl(lden0),bl(lvec),bl(lval),
     *                 sxyz,fxyz,bl(ld1c), dxyz,lsemi)
c
      call retmark
c----------------------------------------------------------------------
c     call secund(tchfe)
c     call elapsec(etchfe)
c     tchf=(tchfe-tchfb)/60.0d0
c     elaps=(etchfe-etchfb)/60.0d0
c     write(iout,430) tchf,elaps
c 430 format(/'Master CPU time for   CPHF  solver    = '
c    *,f8.2,'  Elapsed = ',f8.2,' min'/)
c----------------------------------------------------------------------
      end
c===================================================================
      subroutine cphf_infoxyz(cphf_thresh)
      implicit real*8 (a-h,o-z)
c
      call getival('iout',iout)
      call getival('print',iprint)
c----------------------------------------------------
      write(iout,100)
  100 format(/72('-'))
      write(iout,*)' '
      write(iout,*)'            The CPHF Solver '
      write(iout,*)' '
      call f_lush(iout)
c----------------------------------------------------
      end
c===================================================================
C    Coupled-perturbed HARTREE-FOCK solver.
c
      subroutine chfsol_xyz1(idft,ax,bl,inx,nocc,ncf,ntri,thrs,
     *                       dens0, vec,val, sxyz,fxyz,d1con,
     *                       dxyz,lsemi)
      implicit real*8 (a-h,o-z)
      common /lindvec/ lind,idensp
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension dens0(ntri)               ! zeroth-order density
      dimension vec(*),val(*)
c------------------------------------------------
      dimension sxyz(ntri,3),fxyz(ntri,3)           !  (ntri,3) overl & Fock
      dimension d1con(ntri,3)
      dimension dxyz(ntri,3)              ! output Ist-order density
c
      parameter (Zero=0.0d0)
c----------------------------------------------------
      call getival('print',iprint)
c------------------------------------------------
      ntri3=3*ntri
c------------------------------------------------
c    Calculate the constant part of the D10
c
c   1. contributions from [ H1 + F1(d0,g1) - ei S1 ]
c   2. contributions from the 0.5 D0*S1*D0
c   The total constant part of the D1 stored in d1con(*)
c
c
c
c for rhf=.true.
c
      call d1const_xyz(.true.,   ncf,    ntri,   nocc,   bl(lind),
     *                 vec,      val,    dens0,  fxyz,   sxyz,
     *                 d1con,    0,      bl,     bl )
c
      irec=1
      call save1mat(65,irec,ntri*3,d1con ) ! d1const
c
ctest if(iprint.ge.2) then
c        write(6,*) '  constant part of D1 matrix '
c        call drumh(d1con(1,1),ncf, 6  ,'dens1-X ')
c        call drumh(d1con(1,2),ncf, 6  ,'dens1-Y ')
c        call drumh(d1con(1,3),ncf, 6  ,'dens1-Z ')
ctest endif
c------------------------------------------------
C07   IF ((idft.eq.0) .or. (ax.NE.Zero)) THEN
c
         cphf_thres=thrs
c
c        Solve the CPHF equation for the D10 matrix
c           starting from D10constant
c
         call solver_xyz1(idft,ax,inx,nocc,cphf_thres,ncf,ntri, bl,
     *                   dxyz,d1con,sxyz,val,vec,fxyz,bl(lind),lsemi)
c.......................................
      if(iprint.ge.2) then
         write(6,*) ' 1st-order density matrix '
         call drumh(dxyz(1,1),ncf, 6  ,'dens1-X ')
         call drumh(dxyz(1,2),ncf, 6  ,'dens1-Y ')
         call drumh(dxyz(1,3),ncf, 6  ,'dens1-Z ')
      endif
c.......................................
C     ELSE
C       do not do cphf for non-hybrid dft : d1=d1cons
C     ENDIF
c------------------------------------------------
      end
c======================================================================
      subroutine cphf_xyzcheck

      use memory

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c
      call getmem(0,last0)
      call retmem(1)
c
      natoms=na
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
c
      nmemory=  ntri3                      ! resulting D1
     *       + 6*ncf                      ! in dens1 for w1,w2(3*ncf)
c
      last1=last0+nmemory
      if(last1 .ge. lcore) then
         call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
         iout=igetival('iout')
         write (iout,2001) last1-ioffset,lcore-ioffset
         call nerror(1,'CPHF ','Memory: needed and available ',
     *               last1-ioffset,lcore-ioffset)
      endif
c
 2001 format(1x,' more memory needed in CPHF : needed=',i5,
     *      '  available=',i5)
c
      end
c======================================================================
      subroutine solver_xyz1(idft,ax,inx,nocc,thrs,ncf,ntri,bl,
     *                       d1 ,d1con, r1 ,val,vec, g1 ,lind,lsemi)
c-------------------------------------------------------------------
c This routine solves the cphf equations .
c-------------------------------------------------------------------
c  1. d1    - resulting first-order density
c  2. d1con - constant part of the first-order density matrix
c  3. r1    - used for residuum r1
c  4. val - eigenvalues of the unperturbed solutions
c  5. vec - eigenvectors of the unperturbed solutions
c  6. g1   - the Fock matrix (first-order)
c  8. lind - the n*(n-1)/2 matrix
c
c thrs - threshold for CPHF (not integrals)
c idft - dft functional type
c ax - factor to multiply exchange contribution in DFT
c-------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical oscylates
      character*4 fext1,fext2      ! name extantion for files
      dimension bl(*)
      dimension inx(12,*)
c-------------------------------------------------------------------
      dimension d1(ntri,3),d1con(ntri,3),r1(ntri,3), g1(ntri,3)
      dimension val(*),vec(*)
      dimension lind(*)
c-------------------------------------------------------------------
      nfile63=63    !  D1 matrix
      nfile64=64    !  R1 matrix
      nfile65=65    ! constant part of D1
c-------------------------------------------------------------------
      residx=1.d+5
      residy=1.d+5
      residz=1.d+5
c-------------------------------------------------------------------
      call getrval('xlvsh',xlvsh)
      if(abs(xlvsh).gt.0.d0) then
         call getival('lsmat',ls0) ! S0
      endif
c-------------------------------------------------------------------
      call mmark
c-------------------------------------------------------------------
      call getival('maxit',mxiter)
c-------------------------------------------------------------------
      call getival('iout',iout)
      call getival('noacc',noacce)
c-------------------------------------------------------------------
      ntr3=3*ntri
c-------------------------------------------------------------------
      factor=2.0d0    ! for dens1_part1 like orbital's occupancy
c-------------------------------------------------------------------
      call secund(tchfb)
      call elapsec(etchfb)
c-------------------------------------------------------------------
      call getmem(ncf**2,lw1)
      nvirt=ncf-nocc
      call getmem(ncf*nvirt,lw2)
c---------------------------------------------------------------
c get input for int_cphf(*) ( fock(d1,g0) builder:
c
      call cphf_input(nblocks,labels)
c
c there is one memory allocation by getmem (for labels)
c---------------------------------------------------------------
c thres1 is the final integral threshold
c thres2 is the integral threshold to begin CPHF with
c
      call getrval('thr2',thres1)
      call getrval('thr3',thres2)
c---------------------------------------------------------------
      mgo=1    ! loose integral threshold
c---------------------------------------------------------------
      if(thres2.eq.thres1) mgo=3    !  final integral threshold
c---------------------------------------------------------------
      write(iout,*) ' '
      write(iout,150) thrs
  150 format('threshold for 1st-order density = ',1pe11.4)
c 150 format('threshold for CPHF convergence  = ',1pe11.4)
c
      thres_i=thres2      ! loose int. thresh. at the begining
      write(iout,151)
  151 format(/'CPHF iter',2x,'res-x',7x,'res-y',7x,'res-z',7x,
     *        'cpu   elapsed  oscillating')
c
      if(mgo.lt.3) then
         write(iout,152) thres_i
      else
         write(iout,154) thres_i
      endif
  152 format('    Loose integral threshold = ',1pe11.4)
  154 format('    Final integral threshold = ',1pe11.4)
c
      call f_lush(iout)
c---------------------------------------------------------------
c last integral threshold used :
c
      thres_last=thres_i
      liter_thre=0     ! counting iterations with the same int. threshold
c---------------------------------------------------------------
c
c           begining of the chf equation solution
c
c---------------------------------------------------------------
      IAT=1
c---------------------------------------------------------------
c The first D1 matrices = D1const
c
      call read1mat(nfile65,iat,ntri*3,r1 ) ! r0=d1const
      call read1mat(nfile65,iat,ntri*3,g1 ) ! r0=d1const
c---------------------------------------------------------------
      call zeroit(d1,ntr3)
c---------------------------------------------------------------
      if(mxiter.eq.0) then
         call read1mat(nfile65,iat,ntri*3,d1 ) ! r0=d1const
cccc     return
         go to 999
      endif
c---------------------------------------------------------------
      call getival('reset',ireset)
      icycle=0
 4321 continue
c---------------------------------------------------------------
      call save1mat(nfile63,iat,ntri*3,r1 )   ! save d0=r0
      call save1mat(nfile64,iat,ntri*3,r1 )   ! save r0=r0
      call save1mat(nfile65,iat,ntri*3,g1)    ! save new r0=r0
c---------------------------------------------------------------
      DO LITER=1,MXITER
         icycle=icycle+1
c
         call get_fext(liter,fext1,fext2)
c                            name.fext1  - files for rl=l(r)
c                            name.fext2  - files for ro=O(rl)
c
c        Calculate the D = C + SUMi [ Li(C)]
c
         call secund(citerb)
         call elapsec(eiterb)
c
c        calculate G(D1,g0) in g1
c
         call f_lush(6)
         call pcphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres_i,mgo,
     *                  r1,g1,bl(labels),lsemi)
c
c next line commented out because the polarizability calculation
c is in fact serial!!   MM 08/18/06
c
c        call para_next(mgo)
c
         if(abs(xlvsh).ne.0.d0) then
c           with level shift calculate SD1S  and add it to the G1
            call getmem(ncf*ncf*3,lsds)
            call sd1s_xyz(bl(lsds),ncf,ntri,lind,bl(ls0),r1,xlvsh,g1)
            call retmem(1)
         endif
c
        do icr=1,3
c         call dens1_1dir1(factor,ncf,ntri,nocc,lind,xlvsh,
c     *                    g1(1,icr),vec,val,r1(1,icr),  !output: r1=L(r0)
c     *                    bl(lw1),bl(lw2) )
         call dens1_1dir1n(factor,ncf,nocc,xlvsh,g1(1,icr),
     1                      vec, val, r1(1,icr),bl(lw1),bl(lw2))
        enddo
c
         call file4cphf_o(ntri*3,iat,fext1,r1,'write') ! save   rl1=L(r1)
c
         IF(NOACCE.EQ.1) THEN
            call calc_d10(nfile65,liter,ntri,iat,d1, r1) ! out: d1,r1
         ELSE
            call make_r_orto(nfile63,liter,ntri,iat,d1con, r1) ! i/o l(r) orthogonal to ..
cc          call chec_r_orto(nfile63,liter,ntri,iat,d1con, r1)
            call file4cphf_o(ntri*3,iat,fext2,r1,'write') ! save  ro1=O(rl1)
c           ......................................
            if(liter.ge.2) then
c                calculate current solution d1
c
                 call getmem(liter*liter,iamat)
                 call getmem(liter*liter,ibmat)
                 call getmem(liter+1    ,iwvec)
                 call getmem(liter      ,ilvec)
                 call getmem(liter      ,imvec)
                 call getmem(liter      ,icoef)
c                use d1con & g1 as scratch arrays
                 call calc_d1r(nfile63,liter,ntri,iat,
     *                        bl(iamat),bl(ibmat),bl(iwvec),
     *                        bl(ilvec),bl(imvec),bl(icoef),
     *                        d1con,r1,  d1)  ! out :current d1  new r1
                 call retmem(6)
c           ......................................
            endif    !  (liter.ge.2) then
            call file4cphf_o(ntri*3,iat,fext2,r1,'read') ! get  r01 for enx
         ENDIF    !   IF(NOACC.EQ.1) THEN
c
         call secund(citere)
         call elapsec(eitere)
         cpuit=citere-citerb
         elait=eitere-eiterb
c
         call cphf_enx(r1,thrs,ntri,mgo,liter,cpuit,elait,lend,
     *                 errmax, residx,residy,residz,oscylates)
c
         IF(LEND.EQ.1) EXIT
c
         if(noacce.eq.1) then
            call file4cphf_o(ntri*3,iat,fext1,r1,'read')  ! get rl1
         else
ccccc       call file4cphf_o(ntri*3,iat,fext2,r1,'read')  ! get ro1=O(rl1)
         endif
c
         call cphf_int_thr(mgo,liter_thre, errmax,oscylates,
     *                     thres_i,thres_last,thres1,thrs)
c
         if(icycle.eq.mxiter) exit
         if(liter.eq.ireset .and. noacce.eq.0) go to 4321
c
      ENDDO    ! end of iterations
c----------------------------------------------
      call file4cphf_c(ntri*3,liter)
c----------------------------------------------
c     end of iterations
c----------------------------------------------
c DO ONE MORE CALCULATION OF D1 USING FULL D1 :
c
      if(liter.eq.icycle) then
         if(liter.eq.1) then  ! 2002 lublin
            call read1mat(nfile65,iat,ntri*3,d1)  ! d1=d1const
         endif
         call pcphf_xyz(idft,ax,nblocks,bl,inx,ntri,thres_i,mgo,
     *                  d1,g1,bl(labels),lsemi)
c
        do icr=1,3
c         call dens1_1dir1(factor,ncf,ntri,nocc, lind,xlvsh,
c     *                    g1(1,icr),vec,val, d1(1,icr),  ! d1= ProjG(D1,g0)
c     *                    bl(lw1),bl(lw2) )
         call dens1_1dir1n(factor,ncf,nocc,xlvsh,g1(1,icr),
     1                      vec, val, d1(1,icr),bl(lw1),bl(lw2))
        enddo
         call read1mat(nfile65,iat,ntri*3,r1 ) ! r1=d1const
         call add2vec(ntri*3,d1,r1,d1 )        ! FINAL D1
      endif
c----------------------------------------------
        if(lend.eq.0) then
           write(iout,440)
           write(iout,441) mxiter,thrs
           write(iout,442)
        endif
c----------------------------------------------------------------------
        call secund(tchfe)
        call elapsec(etchfe)
        tchf=(tchfe-tchfb)/60.0d0
        elaps=(etchfe-etchfb)/60.0d0
        write(iout,430) tchf,elaps
c----------------------------------------------------------------------
  430 format(/'Master CPU time for   CPHF  solver    = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
  440 FORMAT(/
     *'**   YOU DID NOT GET THE CONVERGENCE FOR CHF SOLUTION   **')
  441 FORMAT(
     *'**       mxiter =',i4,'  threshold = ',D12.6,'       **')
  442 FORMAT(
     *'**                        S O R R Y                     **')
c----------------------------------------------------------------------
c     jump here if maxiter=0 i.e. no cphf
 999  continue
c----------------------------------------------------------------------
      call retmark
c
      end
c====================================================================

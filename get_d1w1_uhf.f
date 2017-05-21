c=====================================================================
      subroutine cphf_d1w1_uhf(idft,   ax,     rhf,    bl,     inx,
     $                         ncf,    ntri,   natoms, NQ,     IUNQ)
c--------------------------------------------------------------------
c calculates 1st-order density & weighted density matrices
c for hessian and store them on a disk (units 63,73 & 64)
c ** UHF version
c--------------------------------------------------------------------
c all parameters are input parameters
c
c  idft    -  dft flag
c  ax      -  scaling parameter for exact exchange
c  rhf     -  logical flag for closed-shell
c  bl      -  general storage
c  inx     -  inx(12,ncs)  integer basis set data
c  ncf     -  number of basis functions
c  ntri    -  ncf*(ncf+1)/2
c  natoms  -  total number of atoms
c  NQ      -  number of symmetry-unique atoms
c  IUNQ    -  list of symmetry-unique atoms
c--------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
c----------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
c----------------------------------------------------
      dimension iunq(*)      ! list of symmetry unique atoms
c----------------------------------------------------
      call getival('iout',iout)
c----------------------------------------------------
      call secund(tchfb)
      call elapsec(etchfb)
c----------------------------------------------------
      ntri3=ntri*3
c----------------------------------------------------
c check bl() status and memory needed & available :
c
      call cphf_memcheck
c----------------------------------------------------
c get zeroth-order density & fock pointers:
c
      call getival('ldensi',lden0)   ! alpha density
      call getival('ldensB',mden0)   ! beta density
      call getival('lfock0',lfoc0)   ! alpha Fock matrix
      call getival('lfock0B',mfoc0)  ! beta Fock matrix
c
      call getival('lvec',lvec)      ! alpha eigenvectors
      call getival('lvecB',mvec)     ! beta eigenvectors
      call getival('lval',lval)      ! alpha eigenvalues
      call getival('lvalB',mval)     ! beta eigenvalues
      call getival('nocc',nalpha)    ! no of occupied alpha orbitals
      call getival('noccB',nbeta)    ! no of occupied beta orbitals
c----------------------------------------------------
c put memory mark
c
      call mmark
c----------------------------------------------------
c find out how many atoms can be treated at once in cphf :
c
      call cphf4nat_uhf(NQ,ncf,ntri,ncphf,natonce)
c----------------------------------------------------
c
c allocate local memory :
c
      nmemory=natonce*ntri3
c
      call getmem(nmemory,lgmat1)   ! for alpha G(D1,g0)
      call getmem(nmemory,lgmat1B)  ! for beta G(D1,g0)
      call getmem(ntri3 ,lwork1)   ! for pojected G(D1,g0) in CPHF
c
      call getmem(ntri3 ,ls1   )   ! over1 copy
      call getmem(ntri3 ,lf1   )   ! fock1 copy
c
      call getmem(nmemory,ld1   )   ! working alpha d1 for natonce atoms
      call getmem(nmemory,ld1B  )   ! working beta d1 for natonce atoms
c
      call getint(natonce,natend)
c----------------------------------------------------
c
c     DO CPHF ;  get D1 & W1 for symmetry unique atoms :
c
      call do_ncphf_u(natoms,NQ,natonce,IUNQ,ncphf,
     $                idft,ax,bl,inx,ncf,ntri,nalpha,nbeta,
     $                bl(lden0),bl(mden0),bl(lfoc0),bl(mfoc0),
     $                bl(lvec),bl(mvec),bl(lval),bl(mval),
     $                bl(lgmat1),bl(lgmat1B),bl(lwork1),bl(ls1),bl(lf1),
     $                bl(natend), bl(ld1),bl(ld1B))
c
c----------------------------------------------------------------------
      call secund(tchfe)
      call elapsec(etchfe)
      tchf=(tchfe-tchfb)/60.0d0
      elaps=(etchfe-etchfb)/60.0d0
      write(iout,430) tchf,elaps
  430 format( 'Master CPU time for Dens1 and Wens1   = '
     *,f8.2,'  Elapsed = ',f8.2,' min'/)
c----------------------------------------------------------------------
      call retmark
c
      end
c======================================================================
      subroutine cphf4nat_uhf(natunq,ncf,ntri,ncphf, natonce)

      use memory

c-------------------------------------------------------------
c This routine estimates how many atoms can
c be treated at once in the UHF CPHF procedure
c based on the available memory
c
c Input :
c natunq - number of symmetry unique atoms
c ncf    - number of basis functions
c ntri   - as usually ncf*(ncf+1)/2
c Output :
c ncphf  - how many times the CPHF eq. will be solved
c natonce- how many atoms will be treated at once in CPHF
c-------------------------------------------------------------
      call getival('iout',iout)
c
      call getival('lcore',lcore)
      call getmem(0,last0)
      call retmem(1)
c
      ntri3=ntri*3
      ncf3 =ncf*3
c
c  memory needs for cphf open-shell:
c
c  3*ntri*natonce   alpha first order fock nmatrix
c  3*ntri*natonce   beta first order fock nmatrix
c
c  3*ntri*natonce   alpha first order density nmatrix
c  3*ntri*natonce   beta first order density nmatrix
c
c  3*ntri*natonce   work area for transposing matrices
c                   subroutine trsp_inplace
c
c  7*natonce        residuals of cphf iterations
c
c  3*ntri           projected fock1 in cphf (one atom at a time)
c  3*ntri           first order overlap (one atom at a time)
c  3*ntri           area for fock1 copy (one atom a a time)
c
c  3*ntri           area for w1A in para_wdens
c  3*ntri           area for w1B in para_wdens
c
c  3*ncf*ncf
c-------------------------------------------------------------
c  needed= 15*ntaonce*ntri +11*ntri + 3*ncf*ncf
c         + 13*natonce + 18*natonce*natoncX
c
c    where natoncX=natonce approx. by natoncx=min(natunq,50)
c-------------------------------------------------------------
      memory1=26*ntri + 31 + 3*ncf*ncf
c-------------------------------------------------------------
c We have :
c
      nmemory=lcore-last0
c leave some mem for integrals :
      nmemory=nmemory - 1000000
c-------------------------------------------------------------
      if(memory1.gt.nmemory) then
         write(iout,*)' not enough memory to handle one atom in cphf'
         write(iout,*)' needed =',memory1,' available =',nmemory
         write(iout,*)' increase declared memory by  =',memory1-nmemory
         call para_stop
         stop ' stopped because of memory problem in cphf'
      endif
c-------------------------------------------------------------
      natoncx=min(natunq,50)
      natonce = (nmemory -11*ntri -3*ncf*ncf)/(15*ntri+13+18*natoncx)
c
      if(natonce.gt.natunq) natonce=natunq
c-------------------------------------------------------------
ctest
c     natonce = 2
c     ncphf=1
c     return
c-------------------------------------------------------------
      if(natonce.eq.0) then
         write(iout,*)' not enough memory to handle one atom in cphf'
         write(iout,*)' needed =',memory1,' available =',nmemory
         write(iout,*)' increase declared memory by  =',memory1-nmemory
         call para_stop
         stop ' stopped because of memory problem in cphf'
      else
         ntimes=natunq/natonce
         ncphf =ntimes
         moreat=natunq-ntimes*natonce
         if(moreat.ne.0) ntimes=ntimes+1
c2003 to make more equal pieces :
         natonce=natunq/ntimes
         if(natunq.gt.ntimes*natonce) natonce=natonce+1
c2003
         if(ntimes.eq.1) write(iout,100)
  100    format('cphf will be solved only once ')
         if(ntimes.gt.1) write(iout,101) ntimes
  101    format('cphf will be solved ',i5,' times')
         write(iout,*) '  '
      endif
c-------------------------------------------------------------
c
      end
c===================================================================
      subroutine do_ncphf_u(natoms,natunq,natonce,listunq,ncphf,
     *                    idft,ax,bl,inx,ncf,ntri,nalpha,nbeta,
     *                    dens0A,dens0B,fock0A,fock0B,
     *                    vecA,vecB,valA,valB,
     *                    gmat1A,gmat1B,work1,s1,f1,
     *                    natend, d1A,d1B)
c------------------------------------------------------------------
c This is the open-shell CPHF driver. It loops over batches of
c perturbations (atoms) ncphf times doing CPHF for natonce atoms
c at the time.
c------------------------------------------------------------------
c INPUT :
c
c natoms   - number of atoms
c natunq   - number of symmetry unique atoms (only these are consider)
c natonce  - number of atoms treated at once in cphf
c listunq  - list of symmetry unique atoms
c ncphf    - how many times (batches of atoms) cphf will be performed
c            ncphf or ncphf+1 times
c idft     - dft flag
c  ax      - dft scaling parameter
c  bl      - storage for everything
c inx      - inx(12,ncs)  basis set info
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c nalpha   - number of alpha occupied orbitals
c nbeta    - number of beta occupied orbitals
c dens0A   - unperturbed alpha density matrix
c dens0B   - unperturbed beta density matrix
c fock0A   - unperturbed alpha fock matrix
c fock0B   - unperturbed beta fock matrix
c vecA     - alpha eigenvectors
c vecB     - beta eigenvectors
c valA     - alpha eigenvalues
c valB     - beta eigenvalues
c
c SCRATCH arrays :
c
c gmat1,work1 - for constant part of D1 , G(D1,g0) and scratch
c
c INPUT :
c s1,f1    - Ist-order overlap and Fock matrices (f1=h1+G(d0,g1))
c            for one atom only (reused)
c
c SCRATCH arrays :
c
c natend   - natend(iat)= 0 or 1 showing if cphf for IAT converaged
c d1, r1(ntri,3,natonce) - to store Dens1 & Resid1
c                                      as needed for cphf acceleration
c OUTPUT :
c
c d1A(ntri,3,natonce)  - final alpha 1st-order density
c d1B(ntri,3,natonce)  - final beta 1st-order density
c WRITTEN ON A DISK
c------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /lindvec/ lind,idensp
      dimension bl(*)
      dimension inx(12,*)
      dimension listunq(*)
      dimension dens0A(ntri),dens0B(ntri),fock0A(ntri),fock0B(ntri)
      dimension vecA(ncf,ncf),vecB(ncf,ncf),valA(*),valB(*)
      dimension s1(ntri,3),f1(ntri,3)
      dimension work1(ntri,3)
      dimension gmat1A(ntri,3,NATONCE),gmat1B(ntri,3,NATONCE)
      dimension d1A(ntri,3,NATONCE),d1B(ntri,3,NATONCE)
      dimension natend(natonce)
c-----------------------------------------------------------------------
      call getival('iout',iout)
      call getrval('xlvsh',xlvsh)
      call getival('lsmat',ls0)
c-----------------------------------------------------------------------
      nfile61=61      !   fock1 alpha                to be read
      nfile71=71      !   fock1 beta                 to be read
      nfile62=62      !   overlap1                   to be read
      nfile63=63      !   density1 alpha             to be written
      nfile73=73      !   density1 beta              to be written
      nfile64=64      !   weighted density1          to be written
      nfile65=65      !   constant part of D1 alpha  to be written
      nfile75=75      !   constant part of D1 beta   to be written
c-----------------------------------------------------------------------
c  HIGHLY DUBIOUS - set here; apparently used in main CPHF routine
c  via depository (THIS USAGE SHOULD BE HIGHLY DISCOURAGED   JB)
c
      call getival('atmass',iatmass)
c-----------------------------------------------------------------------
c For Better convergence criterion: tr(h1+g1)d1  and tr[s1(d1fd+dfd1)]
c calculate  FD   ; write on a disk (60)
c
      lrec=1
      call calcFD(lrec,ncf,ntri,fock0A,dens0A,s1) ! s1 as a scratch
c
      lrec=natonce*3+2
      call calcFD(lrec,ncf,ntri,fock0B,dens0B,s1) ! s1 as a scratch
c-----------------------------------------------------------------------
      ncf2=ncf*ncf
c-----------------------------------------------------------------------
c Solve the CPHF equations for NUMBER of atoms at once :
c
      call getrval('thres',thresh)
c
      nat_b=-natonce+1
      nat_e=0
      npos = -natonce
      do icphf=1,ncphf
         nat_b=nat_b+natonce
         nat_e=nat_e+natonce
         npos = npos + natonce
         write(iout,1234) icphf,natonce
 1234    format('  Solving the CPHF equations ',I3,' time for ',I3,
     $          ' atoms')
         call f_lush(iout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         iat1=0
         do iat=nat_b,nat_e
         iat1=iat1+1
         nat=listunq(iat)
c -- get mass of atom nat
         call getmass1(nat,iat1,natoms,bl(iatmass))
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c -- calculate constant part of 1st-order density and compute
c    FDS1 matrices  This can now be done in parallel
c
         call getmem(3*ncf*ncf,ifdsA)
         call getmem(3*ncf*ncf,ifdsB)
         call para_d1const_uhf(natoms, nat_b,  nat_e,  listunq,ncf,
     $                         ntri,   nalpha, nbeta, bl(lind),dens0A,
     $                         dens0B, s1,     f1,    vecA,    vecB,
     $                         valA,   valB,   bl,  bl(ifdsA),bl(ifdsB),
     $                         work1)
         call retmem(2)
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c        ...............................................................
c        SOLVE THE cphf ; GET D1
c
c        use s1(3*ntri) as a scratch for alpha D1constant (1 atom only)
c        use f1(3*ntri) as a scratch for beta  D1constant (1 atom only)
c
         call chfsol_nat_u(natoms,natunq,listunq,natonce,npos,
     $                     idft,ax,bl,inx,nalpha,nbeta,ncf,ntri,
     $                     bl(lind),thresh,vecA,vecB,valA,valB,s1,f1,
     $                     gmat1A,gmat1B,work1, natend,
     $                     d1A,d1B)
c
c        d1A and d1B contain final 1st-order densities for
c        the current batch of natonce atoms. Now compute
c        the corresponding first order weighted density
c        ...............................................................
c
         call para_done(-1)
c
c        read in precomputed FD matrices from disk, use s1/work1 as scratch
         call read1mat(60,  1 ,ncf2,s1)                ! lrec=1
         call read1mat(60, natonce*3+2 ,ncf2,work1)    ! lrec=natonce*3+2
c
c -- now compute 1st-order weighted densities
c    this can now be done in parallel
c
         call getmem(ncf2*3,lww)
         call getmem(ntri*3,lwA)
         call getmem(ntri*3,lwB)
         call para_wdensu_xyz(nat_b,  nat_e,  listunq,ncf,    ntri,
     $                       bl(lind),dens0A, dens0B, s1,     work1,
     $                        f1,     gmat1A, gmat1B, d1A,    d1B,
     $                        bl(lww),bl(lwA),bl(lwB))
         call retmem(3)
c
         write(iout,*)' '
      enddo
c
c
c     DO REMAINING (natunq-nat_e) ATOMS :
c
      natodo=natunq-nat_e
      npos = npos + natodo
      if(natodo.gt.0) then
         natonce=natodo
c-----------------------------------------------------------------------
c calculate  FD   ; write on a disk (60)
c
         lrec=1
         call calcFD(lrec,ncf,ntri,fock0A,dens0A,s1) ! s1 as a scratch
c
         lrec=natonce*3+2
         call calcFD(lrec,ncf,ntri,fock0B,dens0B,s1) ! s1 as a scratch
c-----------------------------------------------------------------------
         write(iout,1234) ncphf+1,natodo
         call f_lush(iout)
c
         nat_b=nat_e+1
         nat_e=natunq
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         iat1=0
         do iat=nat_b,nat_e
         iat1=iat1+1
         nat=listunq(iat)
c -- get mass of atom nat
         call getmass1(nat,iat1,natoms,bl(iatmass))
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c -- calculate constant part of 1st-order density and compute
c    FDS1 matrices  This can now be done in parallel
c
         call getmem(3*ncf*ncf,ifdsA)
         call getmem(3*ncf*ncf,ifdsB)
         call para_d1const_uhf(natoms, nat_b,  nat_e,  listunq,ncf,
     $                         ntri,   nalpha, nbeta, bl(lind),dens0A,
     $                         dens0B, s1,     f1,    vecA,    vecB,
     $                         valA,   valB,   bl,  bl(ifdsA),bl(ifdsB),
     $                         work1)
         call retmem(2)
c
c        ...............................................................
c        SOLVE THE CPHF ; get D1
c
c        use s1(3*ntri) as a scratch for alpha D1constant (1 atom only)
c        use f1(3*ntri) as a scratch for beta  D1constant (1 atom only)
c
         call chfsol_nat_u(natoms,natunq,listunq,natodo,npos,
     $                     idft,ax,bl,inx,nalpha,nbeta,ncf,ntri,
     $                     bl(lind),thresh,vecA,vecB,valA,valB,s1,f1,
     $                     gmat1A,gmat1B,work1, natend,
     $                     d1A,d1B)
c
c        d1A and d1B contain final Ist-order densities for
c        the current batch of natonce atoms. Now compute
c        the corresponding first order weighted density
c        ...............................................................
c
         call para_done(-1)
c
c        read in precomputed FD matrices from disk, use s1/work1 as scratch
         call read1mat(60,  1 ,ncf2,s1)                ! lrec=1
         call read1mat(60, natonce*3+2 ,ncf2,work1)    ! lrec=natonce*3+2
c
c -- now compute 1st-order weighted densities
c    this can now be done in parallel
c
         call getmem(ncf2*3,lww)
         call getmem(ntri*3,lwA)
         call getmem(ntri*3,lwB)
         call para_wdensu_xyz(nat_b,  nat_e,  listunq,ncf,    ntri,
     $                       bl(lind),dens0A, dens0B, s1,     work1,
     $                        f1,     gmat1A, gmat1B, d1A,    d1B,
     $                        bl(lww),bl(lwA),bl(lwB))
         call retmem(3)
c
      endif
c
c -- tell slaves we are completely done
      call para_done(-2)
c
      write(iout,*)' '
c-------------------------------------------------------------------
      end
c===================================================================
      subroutine chfsol_nat_u(natoms,natunq,listunq,natonce,npos,
     $                      idft,ax,bl,inx,nalpha,nbeta,ncf,ntri,
     $                      lind,thrs,vecA,vecB,valA,valB,d1conA,d1conB,
     $                      gmat1A,gmat1B,
     *                      work1, natend, d1A,d1B)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c This routine solves the CPHF equations for NATONCE atoms at once
c At the end it calculates the G(D1,g0) with a final D1 as needed
c ***  Open-shell version
c-----------------------------------------------------------------
c INPUT :
c
c natoms   - number of atoms
c natunq   - number of symmetry-unique atoms (only these are considered)
c listunq  - list of symmetry-unique atoms
c natonce  - number of atoms treated at once in cphf
c npos     - current starting atom
c idft     - dft flag
c  ax      - dft scaling parameter
c  bl      - storage for everything
c inx      - inx(12,ncs)  basis set info
c nalpha   - number of alpha electrons
c nbeta   - number of beta electrons
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c thrs     - threshold for cphf
c vecA     - alpha eigenvectors
c vecB     - beta eigenvectors
c valA     - alpha eigen values
c valB     - beta eigen values
c d1conA   - storage for constant part of alpha D1 for ONE atom only :
c d1conB   - storage for constant part of beta D1 for ONE atom only :
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eoS1)Cv*[CoCv+ + CvCo+]
c gmat1A   - storage for alpha G(D1,g0)
c gmat1B   - storage for alpha G(D1,g0)
c natend   - storage for a vector showing if cphf for a given atom is
c            done or not (natend(iat)=0 or 1 : not converged or converged)
c
c OUTPUT :
c
c d1A(ntri,3,natonce) - final alpha Ist-order density
c d1B(ntri,3,natonce) - final beta Ist-order density
c-----------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension lind(*)
      dimension vecA(*),vecB(*),valA(*),valB(*)
c------------------------------------------------
      dimension d1conA(ntri,3),d1conB(ntri,3), work1(ntri,3)
      dimension gmat1A(ntri,3,NATONCE),gmat1B(ntri,3,NATONCE)
c
      dimension d1A(ntri,3,NATONCE),d1B(ntri,3,NATONCE)
c
      dimension natend(natonce),listunq(natunq)
c------------------------------------------------
      call getival('printh',iprint)
      call getival('iout',iout)
c------------------------------------------------
      call secund(tchfb)
      call elapsec(etchfb)
c------------------------------------------------
      call mmark
c------------------------------------------------
c allocate memory for : residuals in D1 and Hess
c                       resx,resy,resz,
c                       hesx,hesy,hesz,
c                       deltx,delty,deltz,
c                       heltx,helty,heltz
c
      call getmem(natonce,iresx)
      call getmem(natonce,iresy)
      call getmem(natonce,iresz)
c
      call getmem(natonce,ihesx)
      call getmem(natonce,ihesy)
      call getmem(natonce,ihesz)
c
      call getmem(natonce,ideltx)
      call getmem(natonce,idelty)
      call getmem(natonce,ideltz)
c
      call getmem(natonce,iheltx)
      call getmem(natonce,ihelty)
      call getmem(natonce,iheltz)
c------------------------------------------------
      ntri3=3*ntri
c------------------------------------------------
         cphf_thres=thrs
c
c        Solve the CPHF equation for the D10 matrix
c             starting from D10constant
c
         call solver_nat_u(natoms, natonce,npos,   natunq, listunq,
     $                     idft,   ax,     inx,    nalpha,  nbeta,
     $                     ncf,    ntri,   lind,   bl,     cphf_thres,
     $                     vecA,   vecB,   valA,   valB,   d1conA,
     $                     d1conB, gmat1A, gmat1B, natend,
     *                   work1,
     *                   bl(iresx),   bl(iresy),   bl(iresz),
     *                   bl(ihesx),   bl(ihesy),   bl(ihesz),
     *                   bl(ideltx),  bl(idelty),  bl(ideltz),
     *                   bl(iheltx),  bl(ihelty),  bl(iheltz),
     *                     d1A,    d1B)
c
      call secund(tchfe)
      call elapsec(etchfe)
      tchf=(tchfe-tchfb)/60.0d0
      elaps=(etchfe-etchfb)/60.0d0
      write(iout,430) tchf,elaps
  430 format(/'Master CPU time for   CPHF  solver    = '
     *,f8.2,'  Elapsed = ',f8.2,' min' )
c
      call f_lush(iout)
c------------------------------------------------
      call retmark
c------------------------------------------------
c
      end
c======================================================================
      subroutine solver_nat_u(natoms, natonce,npos,   natunq, listunq,
     $                        idft,   ax,     inx,    nalpha, nbeta,
     $                        ncf,    ntri,   lind,   bl,     thrs,
     $                        vecA,   vecB,   valA,   valB,   d1conA,
     $                        d1conB, g1A,    g1B,    natend, work1,
     *                        resx,   resy,   resz,   hesx,   hesy,
     *                        hesz,   deltx,  delty,  deltz,  heltx,
     *                        helty,  heltz,  r1A,    r1B)
c-------------------------------------------------------------------
c This routine solves the chfp equations for NATONCE atoms
c INPUT :
c
c natoms   - total number of atoms
c natonce  - number of atoms treated at once in this cphf batch
c npos     - position of current starting atom
c natunq   - number of symmetry-unique atoms
c listunq  - list of symmetry-unique atoms
c idft     - dft flag
c  ax      - dft scaling parameter
c  bl      - storage for everything
c inx      - inx(12,ncs)  basis set info
c nalpha   - number of alpha electrons
c nbeta    - number of beta electrons
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c thrs     - threshold for cphf
c vecA     - alpha eigenvectors
c vecB     - beta eigenvectors
c valA     - alpha eigenvalues
c valB     - beta eigenvalues
c d1conA   - storage for constant part of alpha D1 for ONE atom only:
c d1conB   - storage for constant part of beta D1 for ONE atom only:
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eoS1)Cv*[CoCv+ + CvCo+]
c g1A      - storage for alpha G(D1,g0)
c g1B      - storage for beta G(D1,g0)
c natend   - storage for a vector showing if cphf for a given atom is
c            done or not (natend(iat)=0 or 1 : not converged or converged)
c resx,resy,resz(iat) - max residals for each atom in cphf
c r1A,r1B(ntri,3,natonce) - to store Dens1 & Resid1
c
c OUTPUT :
c
c R1A(ntri,3,natonce) - final 1st-order alpha density
c R1B(ntri,3,natonce) - final 1st-order beta density
c-------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2 ! name extention for files
c
      character*4 screen
      common /screen_type/ screen
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      logical oscylates
      dimension bl(*)
      dimension inx(12,*)
      dimension lind(*)
      dimension vecA(*),vecB(*),valA(*),valB(*)
c-------------------------------------------------------------------
      dimension natend(natonce),listunq(natunq)
c-------------------------------------------------------------------
      dimension resx(natonce),resy(natonce),resz(natonce)
      dimension hesx(natonce),hesy(natonce),hesz(natonce)
c
      dimension deltx(natonce),delty(natonce),deltz(natonce)
      dimension heltx(natonce),helty(natonce),heltz(natonce)
c-------------------------------------------------------------------
      dimension d1conA(ntri,3),d1conB(ntri,3), work1(ntri,3)
      dimension g1A(ntri,3,natonce),g1B(ntri,3,natonce)
      dimension r1A(ntri,3,natonce),r1B(ntri,3,natonce)
c
c
c output d1:1st-order density
c-------------------------------------------------------------------
      call getival('ncs ',ncs)
c-------------------------------------------------------------------
      call getival('atmass',iatmass)
c-------------------------------------------------------------------
      call getival('lsmat',ls0)
c-------------------------------------------------------------------
c
c  files used during cphf solution:
c
      nfile63=63    ! alpha D1 matrix
      nfile73=73    ! beta D1 matrix
      nfile65=65    ! constant part of alpha D1
      nfile75=75    ! constant part of beta D1
c
      nfile66=66    ! current alpha rO
      nfile76=76    ! current beta rO
      nfile67=67    ! current alpha d1
      nfile77=77    ! current beta d1
c-------------------------------------------------------------------
      do iat=1,natonce
         resx(iat)=1.d+5
         resy(iat)=1.d+5
         resz(iat)=1.d+5
         natend(iat)=0       ! cphf not converged for this atom yet
         hesx(iat)=1.d+5
         hesy(iat)=1.d+5
         hesz(iat)=1.d+5
      enddo
c-------------------------------------------------------------------
      call getival('printh',iprint)
      call getival('iout',iout)
      call getival('noacc',noacce)
c-------------------------------------------------------------------
      call getrval('xlvsh',xlvsh)
c
c     if(xlvsh.NE.0.d0) then
c        call getival('lsmat',ls0)
c        call drumh(bl(ls0     ),ncf, 6  ,'S matrix')
c     endif
c-------------------------------------------------------------------
      call mmark
c-------------------------------------------------------------------
      call getival('maxit',maxiter)
      mxiter=maxiter
c-------------------------------------------------------------------
      call secund(tchfb)
      call elapsec(etchfb)
c-------------------------------------------------------------------
c     begining of the chf equation solution
c-----------------------------------------------------------------------
c for dens1_part1 :
c
      factor=1.0d0           ! for UHF like orbital's occupancy
c-----------------------------------------------------------------------
      call getmem(ncf**2,lw1)      !
      nvirt=ncf-nbeta
c  assume that Nalpha > Nbeta
      call getmem(ncf*nvirt,lw2)      !
c---------------------------------------------------------------
c get input for int_cphf(*) ( fock(d1,g0) builder:
c
      call cphf_input(nblocks,labels)
c
c there is one memory allocation by getmem (for labels)
c---------------------------------------------------------------
c thres1 is the final integral threshold
c thres2 is the integral threshold to start CPHF with
c
      call getrval('thref1',thres1)    ! final int.thresh in CPHF
      call getrval('thref2',thres2)    ! loose int.thresh in CPHF
c
c if the calculation is semidirect, store integrals now
c
      if(scftype.ne.'full-direct') then
        call pstore(nblocks,bl,inx,thres1,bl(labels))
      endif
c----------------------------------------------
      mgo=1    ! loose integral threshold
      if(thres2.eq.thres1) mgo=3
c
c final integral threshold from the begining
c-------------------------------------------------------------------
      liter=0
      icycle=0
      lend=0
      ncit=0
      lsemi=0      ! for reuse of DFT grid
      if(mxiter.eq.0) go to 3000
c-----------------------------------------------------------------------
      write(iout,151)
  151 format(/'CPHF       Res-x       Res-y       Res-z',
     *              '       cpu   elapsed  oscillating'/
     *                    'iter       H(rx)       H(ry)       H(rz) ' )
c-----------------------------------------------------------------------
c     write(iout,151)
c 151 format(/'CPHF iter',2x,'res-x',7x,'res-y',7x,'res-z',7x,
c    *        'cpu   elapsed  oscillating')
c-----------------------------------------------------------------------
c
      thres_i=thres2      ! loose int. thresh. at the begining
c
      if(mgo.lt.3) then
         write(iout,152) thres_i
      else
         write(iout,154) thres_i
      endif
  152 format('    Loose integral threshold = ',1pe11.4)
  154 format('    Final integral threshold = ',1pe11.4)
      call f_lush(iout)
c---------------------------------------------------------------
c last integral threshold used :
c
      thres_last=thres_i
      liter_thre=0     ! counting iterations with the same int. threshold
c---------------------------------------------------------------
c
c     Begining of the solution of the CPHF equations
c
c---------------------------------------------------------------
c screening in cphf :
c
      call getival('delta',idelt)
      call getival('scree',iscre)
c
      screen='fock'
      if(iscre.ne.1) screen='2pde'
c
c---------------------------------------------------------------
c make screening density from D1const :
c
      if(idelt.eq.0) then
c        take current D1const  from unit=65,75
         call make_screenD_uhf(iscre,r1A,r1B,inx,lind,
     *                         work1,natend,natonce,ncs,ncf,ntri)
      endif
c
c---------------------------------------------------------------
c---------------------------------------------------------------
c The first D1 matrices = D1const
c
      do iat=1,natonce
         call read1mat(nfile65,iat,ntri*3,r1A(1,1,iat) ) ! r0A=d1constA
         call read1mat(nfile75,iat,ntri*3,r1B(1,1,iat) ) ! r0B=d1constB
c
         call save1mat(nfile66,iat,ntri*3,r1A(1,1,iat) ) ! save new r0A=r0A
         call save1mat(nfile76,iat,ntri*3,r1B(1,1,iat) ) ! save new r0B=r0B
         call save1mat(nfile67,iat,ntri*3,r1A(1,1,iat) ) ! current d1A
         call save1mat(nfile77,iat,ntri*3,r1B(1,1,iat) ) ! current d1B
      enddo
c---------------------------------------------------------------
      call getival('reset',ireset)
c---------------------------------------------------------------
c     return address for reset
 4321 continue
c---------------------------------------------------------------
      DO LITER=1,MXITER
c
         call setival('liter',liter)
         icycle=icycle+1
         ncit=ncit+1
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
         call para_cphf_nat_uhf(natoms,natend,natonce,npos,natunq,
     $                          listunq, idft, ax, nblocks, bl,
     $                          inx, ntri, thres_i, lsemi, r1A,
     $                          r1B, g1A, g1B, bl(labels))
c
         DO IAT=1,NATONCE
         ireca=(iat-1)*2+1   ! record position for alpha quantities
         irecb=ireca+1       ! record position for beta quantities
         if(natend(iat).eq.0) then
            if(abs(xlvsh).ne.0.d0) then
c              with level shift calculate SD1S  and add it to the G1
               call getmem(ncf*ncf*3,lsds)
               call sd1s_xyz_uhf(bl(lsds),ncf,ntri,lind,bl(ls0),
     $                       r1A(1,1,iat),xlvsh,g1A(1,1,iat))
               call sd1s_xyz_uhf(bl(lsds),ncf,ntri,lind,bl(ls0),
     $                       r1B(1,1,iat),xlvsh,g1B(1,1,iat))
               call retmem(1)
            endif
c
           do icr=1,3
c  alpha
c            call dens1_1dir1(factor,ncf,ntri,nalpha,lind,xlvsh,
c     *                       g1A(1,icr,iat),vecA,valA,
c     *                       r1A(1,icr,iat),bl(lw1),bl(lw2) ) !output: r1a= L(r0a)
             call dens1_1dir1n(factor,ncf,nalpha,xlvsh,g1A(1,icr,iat),
     1                     vecA, valA, r1A(1,icr,iat), bl(lw1),bl(lw2))
c  beta
c            call dens1_1dir1(factor,ncf,ntri,nbeta,lind,xlvsh,
c     *                       g1B(1,icr,iat),vecB,valB,
c     *                       r1B(1,icr,iat),bl(lw1),bl(lw2) ) !output: r1b= L(r0b)
             call dens1_1dir1n(factor,ncf,nbeta,xlvsh,g1B(1,icr,iat),
     1                     vecB, valB, r1B(1,icr,iat), bl(lw1),bl(lw2))
           enddo
c
            call file4cphf_o(ntri*3,ireca,fext1,r1A(1,1,iat),'write') !rl1a=L(r1a)
            call file4cphf_o(ntri*3,irecb,fext1,r1B(1,1,iat),'write') !rl1b=L(r1b)
c.....................................................................
c
c   FROM HEREAFTER kepp current D1 in G1
c

            if(noacce.eq.1) then
c              no acceleration
               call calc_d10_uhf(nfile65,nfile75,liter,ntri,iat,
     $                           ireca,irecb,G1A(1,1,iat),G1B(1,1,iat),
     $                            d1conA     , d1conB     )
c              output :   current d1 in G1A,G1B
c
               call file4cphf_o(ntri*3,ireca,fext1,r1A(1,1,iat),'read') !rl1a=L(r1a)
               call file4cphf_o(ntri*3,irecb,fext1,r1B(1,1,iat),'read') !rl1b=L(r1b)
c
               go to 4000
c
            endif
c
            call make_r_orto_uhf1(nfile66,nfile76,liter,ntri,iat,
     $                           ireca,irecb,d1conA,d1conB,
     $                           r1A(1,1,iat),r1B(1,1,iat))
c
            call file4cphf_o(ntri*3,ireca,fext2,r1A(1,1,iat),'write') !ro1=O(rl1)
            call file4cphf_o(ntri*3,irecb,fext2,r1B(1,1,iat),'write') !ro1=O(rl1)
c           ......................................
            if(liter.ge.2) then
c             calculate current solution d1
c
              call getmem(liter*liter,iamat)
              call getmem(liter*liter,ibmat)
              call getmem(liter+1    ,iwvec)
              call getmem(liter      ,ilvec)
              call getmem(liter      ,imvec)
              call getmem(liter      ,icoef)
c             use d1conA & d1conB as scratch arrays
c
              call calc_d1r_uhf1(nfile66,nfile76,ireca,irecb,liter,
     $                     ntri,iat,bl(iamat),bl(ibmat),bl(iwvec),
     $                     bl(ilvec),bl(imvec),bl(icoef),
     *                     d1conA,d1conB,G1A(1,1,iat),G1B(1,1,iat))
c                              scratch  , out :  current G1
              call retmem(6)
            endif    !  (liter.ge.2) then
c           ......................................
           if(liter.eq.1) then
c              calculate current d1
               call read1mat(nfile67,iat,ntri*3,G1A(1,1,iat) ) ! last d1
               call read1mat(nfile77,iat,ntri*3,G1B(1,1,iat) ) ! last d1
c
               call file4cphf_o(ntri*3,irecA,fext1,r1A(1,1,iat),'read')
               call file4cphf_o(ntri*3,irecB,fext1,r1B(1,1,iat),'read')
c
               call daxpy(ntri*3,1.0d0,r1A(1,1,iat),1,G1A(1,1,iat),1)    !  current d1
               call daxpy(ntri*3,1.0d0,r1B(1,1,iat),1,G1B(1,1,iat),1)    !  current d1
            endif
c           ......................................
c           for the ending test :
            if(idft.eq.0) then
             call file4cphf_o(ntri*3,irecA,fext2,r1A(1,1,iat),'read')
             call file4cphf_o(ntri*3,irecB,fext2,r1B(1,1,iat),'read')
            else
             call read1mat(nfile67,iat,ntri*3,d1conA ) ! last d1
             call read1mat(nfile77,iat,ntri*3,d1conB ) ! last d1
             call sub2vec(ntri*3,G1A(1,1,iat),d1conA,r1A(1,1,iat))
             call sub2vec(ntri*3,G1B(1,1,iat),d1conB,r1B(1,1,iat))
            endif
c          ......................................
 4000       continue   !   if no acceleraion
c          ......................................
c
c           save current D1 :
c
            call save1mat(nfile67,iat,ntri*3,G1A(1,1,iat) ) ! save d1
            call save1mat(nfile77,iat,ntri*3,G1B(1,1,iat) ) ! save d1
c          ......................................
c
         endif ! if(natend(iat).eq.0) then
         ENDDO ! over atoms
c
c--------------------------------------------------------------------
c for convergence test calculate contributions to hessian :
c                FDS1*D1 and S1DF*D1
c
         call getmem(ncf*ncf*3  ,ifds1)
         call getmem(natonce*natonce*9,ihesA)
         call getmem(natonce*natonce*9,ihesB)
c
c   use d1cona, d1conb  as a scratch for f1 :
c
c   alpha :
         lrec=1
         call HessCont2Diag(lrec ,natoms, natend,natonce,ncf,ntri,R1A,
     *                  bl(ifds1),d1conA ,bl(ihesA) )
c   beta :
         lrec=natonce*3+2
         call HessCont2Diag(lrec ,natoms, natend,natonce,ncf,ntri,R1B,
     *                  bl(ifds1),d1conB ,bl(ihesB) )
c
         call daxpy(natonce*natonce*9,1.0d0,bl(ihesB),1,bl(ihesA),1)
c
         call whessContDiag(bl(ihesA),bl(iatmass),natoms,natonce)
c
         do iat=1,natonce
         if(natend(iat).eq.0) then
          call daxpy(ntri*3,1.0d0,r1B(1,1,iat),1,r1A(1,1,iat),1)
         endif
         enddo
c
         call secund(citere)
         call elapsec(eitere)
         cpuit=citere-citerb
         elait=eitere-eiterb
c
        call cphf_enHR(r1A,bl(ihesA),thrs,ntri,ncf,liter,cpuit,elait,
     *                 lend, mgo, natend,errmax,natonce,
     *                 resx,resy,resz,    hesx,hesy,hesz,
     *                 deltx,delty,deltz, heltx,helty,heltz,
     *                 oscylates)
c
        call retmem(3)
c--------------------------------------------------------------------
c        ............................................................
cok       do iat=1,natonce
cok          if(natend(iat).eq.0) then
cok       call add2vec(ntri*3,r1A(1,1,iat),r1B(1,1,iat),r1A(1,1,iat))
cok          endif
cok       enddo
c
cok      call cphf_enw(r1A,thrs,ntri,mgo,liter,cpuit,elait,lend,natend,
cok  *                 errmax,natonce,resx,resy,resz,oscylates)
c        ............................................................
c--------------------------------------------------------------------
c
c
         IF(LEND.EQ.1 .or. ICYCLE.EQ.MXITER) EXIT
c
         if(noacce.eq.0) then
           do iat=1,natonce
             ireca=(iat-1)*2+1   ! record position for alpha quantities
             irecb=ireca+1       ! record position for beta quantities
             if(natend(iat).eq.0) then
               call file4cphf_o(ntri*3,ireca,fext2,r1A(1,1,iat),'read')
               call file4cphf_o(ntri*3,irecb,fext2,r1B(1,1,iat),'read')
c                  r01 for next iter
             endif ! if(natend(iat).eq.0) then
           enddo ! over atoms
         endif
c
         call cphf_int_thr(mgo,liter_thre, errmax,oscylates,
     *                     thres_i,thres_last,thres1,thrs)
c
         IF(ICYCLE.EQ.MXITER) EXIT
c
         if(liter.eq.ireset .and. noacce.eq.0) then
           do iat=1,natonce
             ireca=(iat-1)*2+1   ! record position for alpha quantities
             irecb=ireca+1       ! record position for beta quantities
             if(natend(iat).eq.0) then
               call file4cphf_o(ntri*3,ireca,fext2,r1A(1,1,iat),'read')
               call file4cphf_o(ntri*3,irecb,fext2,r1B(1,1,iat),'read')
c
               call save1mat(nfile66,iat,ntri*3,G1A(1,1,iat) ) ! save d1
               call save1mat(nfile76,iat,ntri*3,G1B(1,1,iat) ) ! save d1
             endif
           enddo ! over atoms
           call file4cphf_c(ntri*3,liter)
           go to 4321
         endif
c
      ENDDO    ! end of iterations
c---------------------------------------------------------------
      call file4cphf_c(ntri*3,liter)
c---------------------------------------------------------------
 3000 continue
c---------------------------------------------------------------
c read in last solutions D1 (R1)
c
         do iat=1,natonce
            call read1mat(nfile67,iat,ntri*3,R1A(1,1,iat) )
            call read1mat(nfile77,iat,ntri*3,R1B(1,1,iat) )
         enddo
c---------------------------------------------------------------
c DO ONE MORE CALCULATION OF D1 USING FULL D1 & full integral thresh.
c
         idelt=1                           ! use delta dens in last iter.
         call setival('delta',idelt)
ccccc    screen='2pd1'
c
         call secund(citerb)
         call elapsec(eiterb)
c---------------------------------------------------------------
         do iat=1,natonce
           natend(iat)=0       ! cphf not converged for this atom yet
         enddo
c
         ncit=ncit+1
c
c d1a and d1b in R1a , R1B
c
         call para_cphf_nat_uhf(natoms,natend,natonce,npos,natunq,
     $                          listunq, idft, ax, nblocks, bl,
     $                          inx, ntri, thres_i, lsemi, R1A,
     $                          R1B, g1A, g1B, bl(labels))
c
       call secund(t1)
       do iat=1,natonce
c.......  if(abs(xlvsh).ne.0.d0) then
c            with level shift calculate SD1S  and add it to the G1
c            call getmem(ncf*ncf*3,lsds)
c            call sd1s_xyz_uhf(bl(lsds),ncf,ntri,lind,bl(ls0),
c    $                     r1A(1,1,iat),xlvsh,g1A(1,1,iat))
c            call sd1s_xyz_uhf(bl(lsds),ncf,ntri,lind,bl(ls0),
c    $                     r1B(1,1,iat),xlvsh,g1B(1,1,iat))
c            call retmem(1)
c.......  endif
c
c
        do icr=1,3
c  alpha
c         call dens1_1dir1(factor,ncf,ntri,nalpha,lind,xlvsh,
c     *                     g1A(1,icr,iat),vecA,valA,
c     *                     R1A(1,icr,iat),bl(lw1),bl(lw2))
             call dens1_1dir1n(factor,ncf,nalpha,xlvsh,g1A(1,icr,iat),
     1                     vecA, valA, r1A(1,icr,iat), bl(lw1),bl(lw2))
c  beta
c         call dens1_1dir1(factor,ncf,ntri,nbeta,lind,xlvsh,
c     *                     g1B(1,icr,iat),vecB,valB,
c     *                     R1B(1,icr,iat),bl(lw1),bl(lw2))
             call dens1_1dir1n(factor,ncf,nbeta,xlvsh,g1B(1,icr,iat),
     1                     vecB, valB, r1B(1,icr,iat), bl(lw1),bl(lw2))
        enddo
         call read1mat(nfile65,iat,ntri*3,d1conA)
         call read1mat(nfile75,iat,ntri*3,d1conB)
         call daxpy(ntri*3,1.0d0,d1conA,1,R1A(1,1,iat),1)
         call daxpy(ntri*3,1.0d0,d1conB,1,R1B(1,1,iat),1)
c
         call save1mat(nfile66,iat,ntri*3,R1A(1,1,iat)) ! FINAL D1
         call save1mat(nfile76,iat,ntri*3,R1B(1,1,iat)) ! FINAL D1
       enddo
c---------------------------------------------------------------
       call secund(citere)
       call elapsec(eitere)
       cpuit=citere-citerb
       elait=eitere-eiterb
       cput=cpuit/60.d0
       elat=elait/60.d0
c
      write(iout,420) liter+1,cput,elat, .false.
  420 format(i3,5x,
     * ' last iteration with full dens1   ',
     * 1x,2(f8.2,1x),l5)
c
      call f_lush(iout)
c---------------------------------------------------------------
      go to 5555
c---------------------------------------------------------------
c calculate final D1-residuum and Hess-residuum
c
c get last D1 in d1con , final D1 is in d1
c
       do iat=1,natonce
          call read1mat(nfile66,iat,ntri*3, r1A(1,1,iat)) ! final d1A
          call read1mat(nfile67,iat,ntri*3,d1conA) ! last  d1A
          call daxpy(ntri*3,-1.0d0,d1conA,1,r1A(1,1,iat),1) ! FINAL R1
c
          call read1mat(nfile76,iat,ntri*3, r1B(1,1,iat)) ! Final d1B
          call read1mat(nfile77,iat,ntri*3,d1conB) ! last  d1B
          call daxpy(ntri*3,-1.0d0,d1conB,1,r1B(1,1,iat),1) ! FINAL R1
       enddo
c
       call getmem(ncf*ncf*3  ,ifds1)
       call getmem(natonce*natonce*9,ihesA)
       call getmem(natonce*natonce*9,ihesB)
c
c      use d1con as a scratch for f1
c
c   alpha :
            lrec=1
            call HessCont2(lrec ,natoms, natend,natonce,ncf,ntri,R1A,
     *                     bl(ifds1),d1conA ,bl(ihesA) )
c   beta :
            lrec=natonce*3+2
            call HessCont2Diag(lrec,natoms, natend,natonce,ncf,ntri,R1B,
     *                     bl(ifds1),d1conB ,bl(ihesB) )
c
          call daxpy(natonce*natonce*9,1.0d0,bl(ihesB),1,bl(ihesA),1)
          do iat=1,natonce
          call daxpy(ntri*3,1.0d0,r1B(1,1,iat),1,r1A(1,1,iat),1)
          enddo
c
            call whessContDiag(bl(ihesA),bl(iatmass),natoms,natonce)
c
c
c     call cphf_enHR(r1,bl(ihes),thrs,ntri,ncf,liter+1,cpuit,elait,lend,
c    *                natend,errmax,natonce,
c    *                resx,resy,resz,    hesx,hesy,hesz,
c    *                deltx,delty,deltz, heltx,helty,heltz,
c    *                oscylates)
c
       call final_HR1(r1A,bl(ihesA),thrs,ntri,ncf,natonce)
       call retmem(3)
c---------------------------------------------------------------
 5555  continue
c---------------------------------------------------------------
       do iat=1,natonce
          call read1mat(nfile66,iat,ntri*3,R1A(1,1,iat) ) ! FINAL d1
          call read1mat(nfile76,iat,ntri*3,R1B(1,1,iat) ) ! FINAL d1
       enddo
c---------------------------------------------------------------
      if(lend.eq.0) then
         write(iout,440)
         write(iout,441) mxiter,thrs
         write(iout,442)
      endif
c----------------------------------------------------------------------
  440 format(/
     *'--   YOU DID NOT GET THE CONVERGENCE FOR CHF SOLUTION   --')
  441 format(
     *'--      maxiter =',i4,'  threshold = ',D12.6,'         --')
  442 format(
     *'--                        S O R R Y                     --')
c----------------------------------------------------------------------
c -- if semidirect, remove integral files
      If(scftype.ne.'full-direct') then
        call getival('indisk',indisk)
        if(indisk.gt.0) call clos4int()
      EndIf
c
      call retmark
      end
c======================================================================
      subroutine make_r_orto_uhf(nfilea,nfileb,liter,ntri,iat,
     $                      ireca,irecb,
     $                       ra,rb,ra_curr,rb_curr)
      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extantion for files
      dimension  ra_curr(ntri,3),rb_curr(ntri,3)   ! current    r
      dimension  ra(ntri,3),rb(ntri,3)             ! previous   r
      data acc,zero,one /1.d-15 , 0.d0 , 1.d0/
c
c       dn+1 = L(dn) - SUM(i=1,n)[ <di|L(dn)>/<di|di> * di ]
c
c e.g.       d2 = L(d1) - <d1|L(d1)>/<d1|d1> * d1
c
c the above is repeated for alpha and beta components
c
      do istep=0,liter-1
         if(istep.eq.0) then
            call read1mat(nfilea,iat,ntri*3,ra) ! very first residuum r0
            call read1mat(nfileb,iat,ntri*3,rb) ! very first residuum r0
         else
            call get_fext(istep,fext1,fext2)
            call file4cphf_o(ntri*3,ireca,fext2,ra,'read') ! previous ro
            call file4cphf_o(ntri*3,irecb,fext2,rb,'read') ! previous ro
         endif
c
c        calculate scalars <r|l(r_curr)> & <r|r>
c
         dxldxa=ddot(ntri,ra(1,1),1,ra_curr(1,1),1)
         dyldya=ddot(ntri,ra(1,2),1,ra_curr(1,2),1)
         dzldza=ddot(ntri,ra(1,3),1,ra_curr(1,3),1)
         dxldxb=ddot(ntri,rb(1,1),1,rb_curr(1,1),1)
         dyldyb=ddot(ntri,rb(1,2),1,rb_curr(1,2),1)
         dzldzb=ddot(ntri,rb(1,3),1,rb_curr(1,3),1)
c
         dx_dxa=ddot(ntri,ra(1,1),1,ra(1,1),1)
         dy_dya=ddot(ntri,ra(1,2),1,ra(1,2),1)
         dz_dza=ddot(ntri,ra(1,3),1,ra(1,3),1)
         dx_dxb=ddot(ntri,rb(1,1),1,rb(1,1),1)
         dy_dyb=ddot(ntri,rb(1,2),1,rb(1,2),1)
         dz_dzb=ddot(ntri,rb(1,3),1,rb(1,3),1)
c
         if(dx_dxa.gt.acc) then
            scxa=dxldxa/dx_dxa
         else
            scxa=zero
         endif
         if(dx_dxb.gt.acc) then
            scxb=dxldxb/dx_dxb
         else
            scxb=zero
         endif
c
         if(dy_dya.gt.acc) then
            scya=dyldya/dy_dya
         else
            scya=zero
         endif
         if(dy_dyb.gt.acc) then
            scyb=dyldyb/dy_dyb
         else
            scyb=zero
         endif
c
         if(dz_dza.gt.acc) then
            scza=dzldza/dz_dza
         else
            scza=zero
         endif
         if(dz_dzb.gt.acc) then
            sczb=dzldzb/dz_dzb
         else
            sczb=zero
         endif
c
c         do i=1,ntri
c            ra_curr(i,1)=ra_curr(i,1)-scxa*ra(i,1)
c            ra_curr(i,2)=ra_curr(i,2)-scya*ra(i,2)
c            ra_curr(i,3)=ra_curr(i,3)-scza*ra(i,3)
c            rb_curr(i,1)=rb_curr(i,1)-scxb*rb(i,1)
c            rb_curr(i,2)=rb_curr(i,2)-scyb*rb(i,2)
c            rb_curr(i,3)=rb_curr(i,3)-sczb*rb(i,3)
c         enddo
         call daxpy(ntri,-scxa,ra(1,1),1,ra_curr(1,1),1)
         call daxpy(ntri,-scya,ra(1,2),1,ra_curr(1,2),1)
         call daxpy(ntri,-scza,ra(1,3),1,ra_curr(1,3),1)
         call daxpy(ntri,-scxb,rb(1,1),1,rb_curr(1,1),1)
         call daxpy(ntri,-scyb,rb(1,2),1,rb_curr(1,2),1)
         call daxpy(ntri,-sczb,rb(1,3),1,rb_curr(1,3),1)
      enddo
c
      end
c======================================================================
      subroutine calc_d10_uhf(nfileA,nfileB,liter,ntri,iat,irecA,irecB,
     $                        dA,dB,rA,rB)
      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extantion for files
      dimension  dA(ntri,3),dB(ntri,3),rA(ntri,3),rB(ntri,3) ! output
c----------------------------------------------------------------------
c read in data :
c
      call read1mat(nfileA,iat,ntri*3,dA)  ! constant part of D1 alpha
      call read1mat(nfileB,iat,ntri*3,dB)  ! constant part of D1 beta
c----------------------------------------------------------------------
        do istep=1,liter
          call get_fext(istep,fext1,fext2)
          call file4cphf_o(ntri*3,irecA,fext1,rA,'read') !  rl alpha
          call file4cphf_o(ntri*3,irecB,fext1,rB,'read') !  rl beta
          call daxpy(ntri*3,1.0d0,rA,1,dA,1) ! d=r0+r1+..
          call daxpy(ntri*3,1.0d0,rB,1,dB,1) ! d=r0+r1+..
        enddo
c----------------------------------------------------------------------
      end
c======================================================================
      subroutine select_dnat_uhf(natend,natonce,densA,densB,ntri,select)
      implicit real*8 (a-h,o-z)
      dimension densA(3,natonce,ntri),densB(3,natonce,ntri),select(ntri)
      dimension natend(natonce)
c
      do i=1,ntri
         dmax=0.d0
         do iat=1,natonce
            if(natend(iat).eq.0) then   ! not converged yet
               x=abs(densA(1,iat,i))+abs(densB(1,iat,i))
               y=abs(densA(2,iat,i))+abs(densB(2,iat,i))
               z=abs(densA(3,iat,i))+abs(densB(3,iat,i))
               dmax=max(x,y,z,dmax)
            endif
         enddo
         select(i)=dmax
      enddo
c
      end
c=====================================================================
      subroutine int_d1g0nat_uhf(
     $        natend,      natonce,    idft,       ax,
     $        nblocks,     bl,         inx,        ntri,
     $        thres1,      map_fs,     denspar,    densA,
     $        densB,       fockA,      fockB,      labels,
     $        mywork,      igran)
c---------------------------------------------------------------------
c for NATONCE atoms at one : cphf for hessian
c   open-shell version
c---------------------------------------------------------------------
c This routine is called from coupled-perturbed HF for all three modes :
c           non -, semi-, and full-direct.
c Three Fock matrices are constructed F(DX,g0),F(DY,g0) and F(DZ,g0)
c using stored integrals first (if any) and then re-calculated integrals
c
c For non - and semi-direct mode the stored integrals are used first
c and corresponding part of Fock matrix is built .Then, for semi-direct
c mode remaining integrals are calculated (with where='fock') and
c building of the Fock matrix is finished.
c
c For recalculated integrals value of where='fock' is set up HERE
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      logical rescale
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
c---------------------------------------------------------------------
c moreint - shows if a given Super-block is
c           finished. If it is true that means
c           that more small blocks are left
c stopnow - if true shows that one tries to get
c           integrals from NOT EXISTING super block
c---------------------------------------------------------------------
c known from int_store   ;
      common /datstore/ thres_stored,isto,isecond,ito
      common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
     *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c---------------------------------------------------------------------
      common /memor1a/ npard,ncost,nsupb, mxsize
c---------------------------------------------------------------------
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*)
      dimension denspar(*),densA(*),densB(*),fockA(*),fockB(*)
      dimension natend(natonce)
      dimension labels(*),map_fs(*)
c
      Parameter (Zero=0.0d0)
c----------------------------------------------------------------
      call secund(time0)
c----------------------------------------------------------------
c symmetry stuff :
c
      call getival('nsym',nsym)
      if(nsym.gt.0) call getival('SymFunPr',ifp)
c----------------------------------------------------------------
c
c  if it is a full-direct run, reset common datstore
c  and skip storage retrieval
c
      if(scftype.eq.'full-direct')then
        isto=0
        ito=0
        isecond=0
      endif
c
c  use stored integrals only once per iteration
c
      if(isto.eq.0) go to 9999
c----------------------------------------------------------------
c First construct a part of fock matrix with stored integrals:
c (those integrals were calculated by calling int_store routine)
c Do it only once each iteration (isto=1 means Use stored!)
c
c       return threshold used for stored integrals:
c
        thres0=thres_stored
c make sure stored integrals are only used once each cycle
        isto=0
c
c get number of stored big-blocks and address of stored quartets array
c
c       first check where integrals are stored (in-core or on-disk):
c
        if(incorex.ne.0) then
c          in-core storage :
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nd1g0_core_nos_uhf(
     $          bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $          bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $          bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $          nblstore,    thres0,     thres1,     bl(lind),
     $          densa,       densb,      focka,      fockb,
     $          ntri,        natonce,    idft,       ax,
     $          denspar,     ncs,        map_fs)
            else
              call nd1g0_core_nosC_uhf(
     $          bl(incorex), bl(ilab),   bl(jlab),   bl(klab),
     $          bl(llab),    bl(isiz),   bl(jsiz),   bl(ksiz),
     $          bl(lsiz),    bl(iqrt),   integrals,  nqstore,
     $          nblstore,    thres0,     thres1,     bl(lind),
     $          densa,       densb,      focka,      fockb,
     $          ntri,       natonce,     denspar,    ncs,
     $          map_fs)
            endif
          else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nd1g0_core_sym_uhf(
     $          nsym,        bl(ifp),    bl(incorex),bl(ilab),
     $          bl(jlab),    bl(klab),   bl(llab),   bl(isiz),
     $          bl(jsiz),    bl(ksiz),   bl(lsiz),   bl(iqrt),
     $          integrals,   nqstore,    nblstore,   thres0,
     $          thres1,      bl(lind),   densa,      densb,
     $          focka,       fockb,      ntri,       natonce,
     $          idft,        ax,         denspar,    ncs,
     $          map_fs)
             else
              call nd1g0_core_symC_uhf(
     $          nsym,        bl(ifp),    bl(incorex),bl(ilab),
     $          bl(jlab),    bl(klab),   bl(llab),   bl(isiz),
     $          bl(jsiz),    bl(ksiz),   bl(lsiz),   bl(iqrt),
     $          integrals,   nqstore,    nblstore,   thres0,
     $          thres1,      bl(lind),   densa,      densb,
     $          focka,       fockb,      ntri,       natonce,
     $          denspar,    ncs,         map_fs)
            endif
          endif
        else
c
c          on-disk storage :
c
           rewind(160)
           rewind(161)
           rewind(162)
           rewind(163)
           rewind(164)
           rewind(165)
           rewind(166)
           rewind(167)
c
           nblcount=0
           rescale=.false.
 9876      continue
c
           call secund(r1)
           call elapsec(e1)
c
           call read_33(160,iblstore,iqstore,nsplit) 
c
           if(iblstore.eq.0) go to 9999
c
           call secund(r2)
           call elapsec(e2)
           readcpu1=readcpu1+(r2-r1)
           readela1=readela1+(e2-e1)
c
           call getint(nsplit*2,ifromto)
           call getint(nsplit*5,nofinteg)
c
           call secund(r1)
           call elapsec(e1)
c
           call read_44(165,bl(ifromto),bl(nofinteg),nsplit)
c
           call secund(r2)
           call elapsec(e2)
           readcpu2=readcpu2+(r2-r1)
           readela2=readela2+(e2-e1)
c
           nblcount=nblcount+iblstore
           if(nblcount.EQ.nblstore) rescale=.true.
c
           call getint_2(4*iblstore , ijklsiz) ! I*2
           call getint_2(  iblstore , nquarts) ! I*2
           call getint_2(4*iqstore  , ijkllab) ! I*2
c
           call getival('intbu72',integbu72)
           call getival('intbu62',integbu62)
           call getival('intbuf2',integbuf2)
           call getival('intbuf4',integbuf4)
           call getival('intbuf8',integbuf8)
c
           call getmem(integbuf8,ixint8)
           call getint_4(integbuf4,ixint4)
           call getint_2(integbuf2,ixint2)
           call getint_2(integbu62,ixin62)
           call getint_2(integbu72,ixin72)
c
c----------------Fock-builders with on-disk integrals------------------
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
                call  nd1g0_disk_nos_uhf(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   densa,      densb,
     $        focka,       fockb,      ntri,       natonce,
     $        idft,        ax,         ncs,        denspar,
     $        map_fs)
            else
                call  nd1g0_disk_nosC_uhf(
     $        ncf,         bl(ixin72), bl(ixin62), bl(ixint2),
     $        bl(ixint4),  bl(ixint8), bl(ijkllab),bl(ijklsiz),
     $        bl(nquarts), iqstore,    iblstore,   rescale,
     $        nsplit,      bl(ifromto),bl(nofinteg),thres0,
     $        thres1,      bl(lind),   densa,      densb,
     $        focka,       fockb,      ntri,       natonce,
     $        ncs,         denspar,    map_fs)
            endif
          else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
                call  nd1g0_disk_sym_uhf(
     $        nsym,        bl(ifp),    ncf,        bl(ixin72),
     $        bl(ixin62),  bl(ixint2), bl(ixint4), bl(ixint8),
     $        bl(ijkllab), bl(ijklsiz),bl(nquarts),iqstore,
     $        iblstore,    rescale,    nsplit,     bl(ifromto),
     $        bl(nofinteg),thres0,     thres1,     bl(lind),
     $        densa,       densb,      focka,      fockb,
     $        ntri,        natonce,    idft,       ax,
     $        ncs,         denspar,    map_fs)
            else
                call  nd1g0_disk_symC_uhf(
     $        nsym,        bl(ifp),    ncf,        bl(ixin72),
     $        bl(ixin62),  bl(ixint2), bl(ixint4), bl(ixint8),
     $        bl(ijkllab), bl(ijklsiz),bl(nquarts),iqstore,
     $        iblstore,    rescale,    nsplit,     bl(ifromto),
     $        bl(nofinteg),thres0,     thres1,     bl(lind),
     $        densa,       densb,      focka,      fockb,
     $        ntri,        natonce,    ncs,        denspar,
     $        map_fs)
            endif
          endif
          call retmem(10)
          if(.not.rescale) go to 9876
        endif
ccccc   call dens_fock_prt(dens,fock,'part_stored')
c
        if(scftype.eq.'non -direct') then
           call secund(time1)
c          call term_info(thres1, 0.d0,time1-time0,where)
C???       return
           go to 8888
        endif
c
 9999 continue
c
c  check if we actually have to compute any integrals,
c  because it might happen that a slave executes this routine
c  even with mywork=-1, in order to retrieve the stored integrals
c  (this is to avoid a race condition MM 08/29/2006 )
c
      if(mywork.lt.0)return
c----------------------------------------------------------------
      call secund(time1)
c----------------------------------------------------------------
c re-calculate integrals
c
      where='fock'
cccc  screen='fock'
c
c set integral timings to zero as well as neglect stuff
c
      call init_info(where)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c total number of super-blocks  nsupblks=nbl2*(nbl2+1)/2
c
      nsupblks=nblocks
      if(mywork.eq.0) then
         istart=1
         istop=nsupblks
      else
         istart=mywork
         istop=mywork+igran-1
      end if
c----------------------------------------------------------------
c
      do isupbl=istart,istop
c
c        get integral's price
c
         call get_price(bl(ncost),isupbl,iprice)
c
c blocks were stored:
c  - if they have negative price and they are before the last stored
c  - if they have negative price and all negative blocks were stored
c  - if all negative blocks were stored and our block is before the last
c
         if((iprice.LE.0.and.isupbl.le.ito).or.
     *      (iprice.le.0.and. isecond.eq.1).or.
     *      (isecond.eq.1.and.isupbl.le.ito)) go to 1111
c
   11 continue
c
        call calcint2(isupbl,bl,inx,thres1,denspar,where,
     *                labels,ibuffz,nblsiz,nintez,ngctoz,
     *                moreint,stopnow)
c
c check if requested super-block was in a list :
c
        if(stopnow) go to 1111
c
c integrals arrived in bl(ibuffz); check if they are :
c
c  nintez is the size of a given quartet.It is set up
c  to zero if there is no integrals
c
        IF(nintez.gt.0) THEN
          call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
          If(nsym.eq.0) Then
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nfock_d1g0_nos_uhf(bl(ibuffz),
     *                            idft,ax,densA,densB,fockA,fockB,
     *                            bl(lind),ntri,natend,natonce,
     *                            denspar,map_fs,ncs,
     *                            nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            Else
c -- "Pure" DFT - no exchange contribution
              call nfock_d1g0_nosC_uhf(bl(ibuffz),densA,densB,fockA,
     $                            fockB,bl(lind),ntri,natend,natonce,
     *                             denspar,map_fs,ncs,
     *                             nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            EndIf
          Else
            If(idft.EQ.0.OR.(idft.GT.0.AND.aX.NE.Zero)) Then
              call nfock_d1g0_sym_uhf(nsym,bl(ifp),bl(ibuffz),
     *                            idft,ax,densA,densB,fockA,fockB,
     *                            bl(lind),ntri,natend,natonce,
     *                            denspar,map_fs,ncs,
     *                            nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            Else
c -- "Pure" DFT - no exchange contribution
              call nfock_d1g0_symC_uhf(nsym,bl(ifp),bl(ibuffz),
     *                             densA,densB,fockA,fockB,
     *                             bl(lind),ntri,natend,natonce,
     *                             denspar,map_fs,ncs,
     *                             nblsiz,ngctoz,nintez,thres1,
     *                          labels(lab1),labels(lab2),labels(lab3))
            EndIf
          EndIf
        ENDIF
c
c----------------------------------------------------------------
c
        if(moreint) go to 11
c
 1111 continue
c
      END DO
c----------------------------------------------------------------
c----------------------------------------------------------------
 8888 continue
c----------------------------------------------------------------
c print integ. timings as well as neglect stuff
c
      call secund(time2)
      call term_info(thres1,time2-time1,time1-time0,where)
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine nfock_d1g0_nos_uhf(buf,idft,ax,densA,densB,fockA,fockB,
     *                          lind,ntri,natend,nat,densp,
     *                          map_fs,ncs,nbls,ngcd,lnijkl,
     *                          thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c NO Symmetry, open-shell
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension densA(3*nat,ntri),fockA(3*nat,ntri)
      dimension densB(3*nat,ntri),fockB(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
      do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
        do 150 iqu=1,ngcq
cccc    do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c          dmax=max(16.d0*dij,16.d0*dkl,dik,dil,djk,djl)
c          dmax=sqrt( dmax )
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
       do 200 iii=1,ilen
          icf=icff+iii
          ii=lind(icf)
       do 200 jjj=1,jlen
          jcf=jcff+jjj
          jj=lind(jcf)
c         output : ijf & density elements xij,yij,zij
       do 200 kkk=1,klen
          kcf=kcff+kkk
          kk=lind(kcf)
       do 200 lll=1,llen
          lcf=lcff+lll
          ll=lind(lcf)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c         magnitude checking again
c
          if(dmax*abs(xint0).gt.0.5d0) then
c
             xint2=xint0*2.d0
             if(idft.gt.0) xint0=xint0*ax
c..........................................................
c            write(6,1234) icf,jcf,kcf,lcf,xint0*thres1
c1234 format(4i3,1x,2(f12.8,1x))
c ....................Coulomb  part ....................
          if(icf.ge.jcf) then
            ijf = ii+jcf
          else
            ijf = jj+icf
          endif
          if(kcf.ge.lcf) then
            klf = kk+lcf
          else
            klf = ll+kcf
          endif
c ....................Exchange part ....................
          if(icf.ge.kcf) then
            ikf = ii+kcf
          else
            ikf = kk+icf
          endif
          if(icf.ge.lcf) then
            ilf = ii+lcf
          else
            ilf = ll+icf
          endif
          if(jcf.ge.kcf) then
            jkf = jj+kcf
          else
            jkf = kk+jcf
          endif
          if(jcf.ge.lcf) then
            jlf = jj+lcf
          else
            jlf = ll+jcf
          endif
c

c..........................................................
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $       (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $       (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $       (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $       (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
       enddo
c ....................Exchange part ....................
c
       do icoordiat=1,nat3
             fockA(icoordiat,ilf)=fockA(icoordiat,ilf) -
     $                            densA(icoordiat,jkf)*xint0
             fockB(icoordiat,ilf)=fockB(icoordiat,ilf) -
     $                            densB(icoordiat,jkf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,jlf)=fockA(icoordiat,jlf) -
     $                            densA(icoordiat,ikf)*xint0
             fockB(icoordiat,jlf)=fockB(icoordiat,jlf) -
     $                            densB(icoordiat,ikf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,ikf)=fockA(icoordiat,ikf) -
     $                            densA(icoordiat,jlf)*xint0
             fockB(icoordiat,ikf)=fockB(icoordiat,ikf) -
     $                            densB(icoordiat,jlf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,jkf)=fockA(icoordiat,jkf) -
     $                            densA(icoordiat,ilf)*xint0
             fockB(icoordiat,jkf)=fockB(icoordiat,jkf) -
     $                            densB(icoordiat,ilf)*xint0
       enddo
c..........................................................
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_sym_uhf(nsym,ifp,buf,idft,ax,
     $                          densA,densB,fockA,fockB,lind,ntri,
     $                          natend,nat,densp,map_fs,ncs,nbls,ngcd,
     $                          lnijkl,thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c SYMMETRY, open-shell
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      dimension buf(nbls,lnijkl,ngcd)
      dimension densA(3*nat,ntri),fockA(3*nat,ntri)
      dimension densB(3*nat,ntri),fockB(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
      do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
        do 150 iqu=1,ngcq
ccc     do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
           dik=densp(ics,kcs)
           dil=densp(ics,lcs)
           djk=densp(jcs,kcs)
           djl=densp(jcs,lcs)
c
           dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c          dmax=max(16.d0*dij,16.d0*dkl,dik,dil,djk,djl)
c          dmax=sqrt( dmax )
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
       do 200 iii=1,ilen
          ic1=icff+iii
          ii1=lind(ic1)
       do 200 jjj=1,jlen
          jc1=jcff+jjj
          jj1=lind(jc1)
       do 200 kkk=1,klen
          kc1=kcff+kkk
          kk1=lind(kc1)
       do 200 lll=1,llen
          lc1=lcff+lll
          ll1=lind(lc1)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c-------------------------------------------------------
c
          IF(dmax*abs(xint0).gt.0.5d0) THEN
c
            xint2=xint0*2.d0
            If(idft.gt.0) xint0=xint0*ax
c-------------------------------------------------------
c ....................Coulomb  part ....................
            if(ic1.ge.jc1) then
              ijf = ii1+jc1
            else
              ijf = jj1+ic1
            endif
            if(kc1.ge.lc1) then
              klf = kk1+lc1
            else
              klf = ll1+kc1
            endif
c ....................Exchange part ....................
            if(ic1.ge.kc1) then
              ikf = ii1+kc1
            else
              ikf = kk1+ic1
            endif
            if(ic1.ge.lc1) then
              ilf = ii1+lc1
            else
              ilf = ll1+ic1
            endif
            if(jc1.ge.kc1) then
              jkf = jj1+kc1
            else
              jkf = kk1+jc1
            endif
            if(jc1.ge.lc1) then
              jlf = jj1+lc1
            else
              jlf = ll1+jc1
            endif
c-------------------------------------------------------
c
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $          (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $          (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $          (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $          (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
       enddo
c ....................Exchange part ....................
c
       do icoordiat=1,nat3
             fockA(icoordiat,ilf)=fockA(icoordiat,ilf) -
     $                            densA(icoordiat,jkf)*xint0
             fockB(icoordiat,ilf)=fockB(icoordiat,ilf) -
     $                            densB(icoordiat,jkf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,jlf)=fockA(icoordiat,jlf) -
     $                            densA(icoordiat,ikf)*xint0
             fockB(icoordiat,jlf)=fockB(icoordiat,jlf) -
     $                            densB(icoordiat,ikf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,ikf)=fockA(icoordiat,ikf) -
     $                            densA(icoordiat,jlf)*xint0
             fockB(icoordiat,ikf)=fockB(icoordiat,ikf) -
     $                            densB(icoordiat,jlf)*xint0
       enddo
c
       do icoordiat=1,nat3
             fockA(icoordiat,jkf)=fockA(icoordiat,jkf) -
     $                            densA(icoordiat,ilf)*xint0
             fockB(icoordiat,jkf)=fockB(icoordiat,jkf) -
     $                            densB(icoordiat,ilf)*xint0
       enddo
c..............................................................
c
       do ns=1,nsym
c
          xcoul=xint2
          xexch=xint0
c
          icf=ifp(ns,ic1)
          jcf=ifp(ns,jc1)
          kcf=ifp(ns,kc1)
          lcf=ifp(ns,lc1)
c
          if(icf.lt.0) then
c            write(6,*)' icf < 0 =',icf
             icf=-icf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(jcf.lt.0) then
c            write(6,*)' jcf < 0 =',jcf
             jcf=-jcf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(kcf.lt.0) then
c            write(6,*)' kcf < 0 =',kcf
             kcf=-kcf
             xcoul=-xcoul
             xexch=-xexch
          endif
          if(lcf.lt.0) then
c            write(6,*)' lcf < 0 =',lcf
             lcf=-lcf
             xcoul=-xcoul
             xexch=-xexch
          endif
c
          ii=lind(icf)
          jj=lind(jcf)
          kk=lind(kcf)
          ll=lind(lcf)
c
c ....................Coulomb  part ....................
          if(icf.ge.jcf) then
            ijf = ii+jcf
          else
            ijf = jj+icf
          endif
          if(kcf.ge.lcf) then
            klf = kk+lcf
          else
            klf = ll+kcf
          endif
c ....................Exchange part ....................
          if(icf.ge.kcf) then
            ikf = ii+kcf
          else
            ikf = kk+icf
          endif
          if(icf.ge.lcf) then
            ilf = ii+lcf
          else
            ilf = ll+icf
          endif
          if(jcf.ge.kcf) then
            jkf = jj+kcf
          else
            jkf = kk+jcf
          endif
          if(jcf.ge.lcf) then
            jlf = jj+lcf
          else
            jlf = ll+jcf
          endif
c
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xcoul
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xcoul
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xcoul
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xcoul
       enddo
c
c ....................Exchange part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ilf)=fockA(icoordiat,ilf) -
     $                            densA(icoordiat,jkf)*xexch
             fockB(icoordiat,ilf)=fockB(icoordiat,ilf) -
     $                            densB(icoordiat,jkf)*xexch
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,jlf)=fockA(icoordiat,jlf) -
     $                            densA(icoordiat,ikf)*xexch
             fockB(icoordiat,jlf)=fockB(icoordiat,jlf) -
     $                            densB(icoordiat,ikf)*xexch
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,ikf)=fockA(icoordiat,ikf) -
     $                            densA(icoordiat,jlf)*xexch
             fockB(icoordiat,ikf)=fockB(icoordiat,ikf) -
     $                            densB(icoordiat,jlf)*xexch
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,jkf)=fockA(icoordiat,jkf) -
     $                            densA(icoordiat,ilf)*xexch
             fockB(icoordiat,jkf)=fockB(icoordiat,jkf) -
     $                            densB(icoordiat,ilf)*xexch
       enddo
c..............................................................
       enddo   !     ns=1,nsym
       ENDIF
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_nosC_uhf(buf,densA,densB,fockA,fockB,
     *                           lind,ntri,natend,nat,densp,
     *                           map_fs,ncs,nbls,ngcd,lnijkl,
     *                           thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c NO Symmetry   ** COULOMB ONLY ** open-shell
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension densA(3*nat,ntri),fockA(3*nat,ntri)
      dimension densB(3*nat,ntri),fockB(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
      do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
        do 150 iqu=1,ngcq
cccc    do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
c
           dmax=4.0d0*max(dij,dkl)
c          dmax=4.0d0*Sqrt(max(dij,dkl))
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
       do 200 iii=1,ilen
          icf=icff+iii
          ii=lind(icf)
       do 200 jjj=1,jlen
          jcf=jcff+jjj
          jj=lind(jcf)
c         output : ijf & density elements xij,yij,zij
       do 200 kkk=1,klen
          kcf=kcff+kkk
          kk=lind(kcf)
       do 200 lll=1,llen
          lcf=lcff+lll
          ll=lind(lcf)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c         magnitude checking again
c
          if(dmax*abs(xint0).gt.0.5d0) then
c
            xint2=xint0*2.d0
c ....................Coulomb  part ....................
          if(icf.ge.jcf) then
            ijf = ii+jcf
          else
            ijf = jj+icf
          endif
          if(kcf.ge.lcf) then
            klf = kk+lcf
          else
            klf = ll+kcf
          endif
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
       enddo
c
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c=====================================================================
      subroutine nfock_d1g0_symC_uhf(nsym,ifp,buf,densA,densB,
     $                           fockA,fockB,lind,ntri,natend,nat,
     $                           densp,map_fs,ncs,nbls,ngcd,lnijkl,
     *                           thres1,labels,length,lgenct)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c SYMMETRY   ** COULOMB ONLY ** open-shell
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      dimension buf(nbls,lnijkl,ngcd)
      dimension densA(3*nat,ntri),fockA(3*nat,ntri)
      dimension densB(3*nat,ntri),fockB(3*nat,ntri)
      dimension natend(nat)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
c----------------------------------------------------------------
      nat3=nat*3
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
      do 100 ijklp=1,nbls
        ngcq=lgenct(ijklp)
        do 150 iqu=1,ngcq
ccc     do 150 iqu=1,ngcd
           icff=labels(1,iqu,ijklp)
           jcff=labels(2,iqu,ijklp)
           kcff=labels(3,iqu,ijklp)
           lcff=labels(4,iqu,ijklp)
c-------------------------------------------------------
c shell :
           ics=map_fs(icff+1)
           jcs=map_fs(jcff+1)
           kcs=map_fs(kcff+1)
           lcs=map_fs(lcff+1)
c
c screening density (squared):
c
           dij=densp(ics,jcs)
           dkl=densp(kcs,lcs)
c
           dmax=4.0d0*max(dij,dkl)
c          dmax=4.0d0*Sqrt(max(dij,dkl))
c-------------------------------------------------------
c  Indices and integrals in the quartet ijkl :
c
          integ=0
       do 200 iii=1,ilen
          ic1=icff+iii
          ii1=lind(ic1)
       do 200 jjj=1,jlen
          jc1=jcff+jjj
          jj1=lind(jc1)
       do 200 kkk=1,klen
          kc1=kcff+kkk
          kk1=lind(kc1)
       do 200 lll=1,llen
          lc1=lcff+lll
          ll1=lind(lc1)
c-------------------------------------------------------
          integ=integ+1
          xint0=buf(ijklp,integ,iqu)
c-------------------------------------------------------
c
          if(dmax*abs(xint0).gt.0.5d0) then
c
             xint2=xint0*2.d0
c ....................Coulomb  part ....................
          if(ic1.ge.jc1) then
            ijf = ii1+jc1
          else
            ijf = jj1+ic1
          endif
          if(kc1.ge.lc1) then
            klf = kk1+lc1
          else
            klf = ll1+kc1
          endif
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xint2
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xint2
       enddo
c..............................................................
c
       do ns=1,nsym
c
          xcoul=xint2
c
          icf=ifp(ns,ic1)
          jcf=ifp(ns,jc1)
          kcf=ifp(ns,kc1)
          lcf=ifp(ns,lc1)
c
          if(icf.lt.0) then
c            write(6,*)' icf < 0 =',icf
             icf=-icf
             xcoul=-xcoul
          endif
          if(jcf.lt.0) then
c            write(6,*)' jcf < 0 =',jcf
             jcf=-jcf
             xcoul=-xcoul
          endif
          if(kcf.lt.0) then
c            write(6,*)' kcf < 0 =',kcf
             kcf=-kcf
             xcoul=-xcoul
          endif
          if(lcf.lt.0) then
c            write(6,*)' lcf < 0 =',lcf
             lcf=-lcf
             xcoul=-xcoul
          endif
c
          ii=lind(icf)
          jj=lind(jcf)
          kk=lind(kcf)
          ll=lind(lcf)
c
c ....................Coulomb  part ....................
          if(icf.ge.jcf) then
            ijf = ii+jcf
          else
            ijf = jj+icf
          endif
          if(kcf.ge.lcf) then
            klf = kk+lcf
          else
            klf = ll+kcf
          endif
c
c ....................Coulomb  part ....................
       do icoordiat=1,nat3
             fockA(icoordiat,ijf)=fockA(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xcoul
             fockB(icoordiat,ijf)=fockB(icoordiat,ijf) +
     $         (densA(icoordiat,klf)+densB(icoordiat,klf))*xcoul
       enddo
       do icoordiat=1,nat3
             fockA(icoordiat,klf)=fockA(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xcoul
             fockB(icoordiat,klf)=fockB(icoordiat,klf) +
     $         (densA(icoordiat,ijf)+densB(icoordiat,ijf))*xcoul
       enddo
c..............................................................
       enddo   !     ns=1,nsym
          endif
c-------------------------------------------------------
  200     continue
  150   continue
  100 continue
c
      end
c======================================================================
      subroutine ds1d_xyz_uhf(work,ncf,ntri,lind,den,ove,hfc)
      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension den(*),ove(*), hfc(*)
      dimension work(3,ncf,ncf)
      data zero,half /0.d0, 0.5d0/
c
      ntr2=ntri*2
c---------------------------------------------
c  calculate  the d0*s1*d0  matrix
c  (this is the uhf version. the corresponding closed-shell version
c   multiplies by an additional factor 0.5)
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               kj=jj+k
               if(k.gt.j) kj=kk+j
               sux=sux+ove(ik)*den(kj)
               suy=suy+ove(ntri+ik)*den(kj)
               suz=suz+ove(ntr2+ik)*den(kj)
            enddo
            work(1,i,j)=sux
            work(2,i,j)=suy
            work(3,i,j)=suz
         enddo
      enddo
c
      ij=0
      do i=1,ncf
         ii=lind(i)
         do j=1,i
            ij=ij+1
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               sux=sux+den(ik)*work(1,k,j)
               suy=suy+den(ik)*work(2,k,j)
               suz=suz+den(ik)*work(3,k,j)
            enddo
c
            hfc(ij)     =hfc(ij)      - sux
            hfc(ntri+ij)=hfc(ntri+ij) - suy
            hfc(ntr2+ij)=hfc(ntr2+ij) - suz
c
         enddo
      enddo
c
      end
c======================================================================
      subroutine f1shift_xyz_uhf(fock1,ncf,ntri,lind,dens0,over0,over1,
     *                       ds1d,bl,xlvsh)
c-----------------------------------------------------------------
c Calculates
c           xlvsh * S0*( D0*S1*D0 )*S0
c needed for F1const with level shift and adds it to Fock1
c
c this is the uhf version. the corresponding rhf routine multiplies
c by an additional factor 0.5
c
c It does it for three directions at once (x,y,z)
c-----------------------------------------------------------------
c INPUT :
c
c fock1    - constant part of F1=H1+G(D0,g1)
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c dens0    - unpreturbed density
c over0    - unpreturbed overlap
c over1    - Ist-order overlap matrix
c ds1d     - scratch for d0*S1*d0 matrix
c bl       - storage for everything
c xlvsh    - level shift
c
c Output :
c
c fock1    - constant part of F1=H1+G(D0,g1) - 0.5*xlvsh*S0*[0.5*D0*S1*D0]*S0
c
c-----------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension bl(*)
      dimension fock1(ntri,3)
      dimension dens0(*)
      dimension over0(*)
      dimension over1(ntri,3)
      dimension ds1d(ntri,3)
c--------------------------------------------------------------
c      return
c--------------------------------------------------------------
c   D0*S1*D0
c
      call getmem(ncf*ncf*3,lwork3)
c
      call zeroit(ds1d,ntri*3)
      call ds1d_xyz_uhf(bl(lwork3),ncf,ntri,lind,dens0, over1,ds1d)
c---------------------------------------------------------------
c        write(6,*)' Fock1 on entry  '
c        call drumh(Fock1(1,1),ncf, 6  ,'Fock1-x ')
c        call drumh(Fock1(1,2),ncf, 6  ,'Fock1-y ')
c        call drumh(Fock1(1,3),ncf, 6  ,'Fock1-z ')
c---------------------------------------------------------------
c we need to add -xlvsh*S( D S1 D ) S  to  Fock1=H1 + G(D0,g1)
c change sign of -DS1D to +
c
      factor=-1.d0
      call dscal(ntri*3,factor,ds1d,1)
c
c - S0*( D0*S1*D0 )*S0* xlvsh  and add it to f1const :
c
      call  sd1s_xyz_uhf(bl(lwork3),ncf,ntri,lind,over0,
     $                   ds1d,xlvsh,fock1)
c---------------------------------------------------------------
c        write(6,*)' Fock1 on exit   '
c        call drumh(Fock1(1,1),ncf, 6  ,'Fock1-x ')
c        call drumh(Fock1(1,2),ncf, 6  ,'Fock1-y ')
c        call drumh(Fock1(1,3),ncf, 6  ,'Fock1-z ')
c---------------------------------------------------------------
c
      call retmem(1)
c--------------------------------------------------------------
      end
c======================================================================
      subroutine sd1s_xyz_uhf(work,ncf,ntri,lind,over0,dens1,xlvsh,hfc)
c
c it calculates S0*D1*S0 as needed for level shifted cphf open shell:
c
      implicit real*8 (a-h,o-z)
      dimension lind(*)
      dimension over0(*),dens1(*), hfc(*)
      dimension work(3,ncf,ncf)
      data zero,half /0.d0, 0.5d0/
c---------------------------------------------
      ntr2=ntri*2
c---------------------------------------------
c  calculate  the  s0*d1*s0  matrix
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               kj=jj+k
               if(k.gt.j) kj=kk+j
               sux=sux+dens1(ik)*over0(kj)
               suy=suy+dens1(ntri+ik)*over0(kj)
               suz=suz+dens1(ntr2+ik)*over0(kj)
            enddo
            work(1,i,j)=sux
            work(2,i,j)=suy
            work(3,i,j)=suz
         enddo
      enddo
c
      ij=0
      do i=1,ncf
         ii=lind(i)
         do j=1,i
            ij=ij+1
            sux=zero
            suy=zero
            suz=zero
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               if(k.gt.i) ik=kk+i
               sux=sux+over0(ik)*work(1,k,j)
               suy=suy+over0(ik)*work(2,k,j)
               suz=suz+over0(ik)*work(3,k,j)
            enddo
c
            hfc(ij)     =hfc(ij)      -xlvsh*sux
            hfc(ntri+ij)=hfc(ntri+ij) -xlvsh*suy
            hfc(ntr2+ij)=hfc(ntr2+ij) -xlvsh*suz
c
         enddo
      enddo
c
      end
c======================================================================
      subroutine calc_d1r_uhf(nfileA,nfileB,irecA,irecB,liter,ntri,iat,
     $                        a,b,w,lv,mv,c,ri,rj,dA,dB,rA,rB)
      implicit real*8 (a-h,o-z)
      character*4 fexi1,fexi2      ! extantion name for files
      character*4 fexj1,fexj2      ! extantion name for files
      dimension ri(ntri,3),rj(ntri,3)
      dimension  dA(ntri,3),dB(ntri,3),rA(ntri,3),rB(ntri,3) ! output
c----------------------------------------------------------------------
      dimension a(liter*liter),b(liter,liter),w(liter+1),c(liter)
      dimension lv(liter), mv(liter)
      dimension liter_dir(3)
c----------------------------------------------------------------------
c input :
c
c liter   - current iteration
c ntri    - ncf*(ncf+1)/2  ; ncf=no of b.f.
c a,b     - arrays (liter x liter) for linear sys. of eq.
c   w       vector (liter)
c lv,mv   - aux. vctors needed for: osinv (b,n,det ,acc,lv,mv)
c
c           a c = w
c
c ri,rj   - storage for residuals from previous iterations
c
c output :
c
c dA,dB,  - resulting 1st-order density matrix (solution at iter=liter)
c rA,rB   - "predicted" residuum , needed ONLY for RESET
c
c
c....................................................................
c Note : the algorithm works like this :
c
c    equation to solve :  D = D0 + L(D) -> Res R= D0 + L(D) -D
c
c  iter
c
c    0      d0=dconst=r0=r0o
c    1      calc: r1=L(r0) ,
c           make r1 orthog. to r0o  r1o=Ort(r1)
c    2      calc: r2=l(r1o),
c           make r2 orthog. to r0o,r1o      r2o=Ort(r2)
c           calc d2=c0*r0o + c1*r1o
c    3      calc: r3=l(r2o),
c           make r3 orthog. to r0o,r1o,r2o  r3o=Ort(r3)
c           calc d3=c0*r0o + c1*r1o + c2*r2o
c    4      calc: r4=l(r3o),
c           make r4 orthog. to r0o,...,r3o,  r4o=Ort(r4)
c           calc d4=c0*r0o + c1*r1o + c2*r2o + c3*r3o
c    5      calc: r5=l(r4o),
c           make r5 orthog. to r0o,...,r4o,  r5o=Ort(r5)
c           calc d5=c0*r0o + c1*r1o + c2*r2o + c3*r3o+c4*r4o
c
c   k+1     calc: r(k+1)=l(rko),
c           make r(k+1) orthog. to r0o,...,rko,  r(k+1)o=Ort(r(k+1))
c           calc d(k+1)=c0*r0o + c1*r1o + ... + c(k)*r(k)o
c....................................................................
c Coefficients { c(0).....c(k) } determined from projections of the
c " predicted " residuum R(K+1) on all previous ORTHOG. resid. {r(k)o}
c
c  <r0o | R(k+1)> =0
c  <r1o | R(k+1)> =0
c  <r2o | R(k+1)> =0
c    .........
c  <rko | R(k+1)> =0
c
c where R(K+1) = r0o + c0*(r1 - r0o)
c                    + c1*(r2 - r1o)
c                    + c2*(r3 - r2o)
c                      ............
c                    + ck*(r(k+1) - rko)
c
c R(K+1) needed to be calculated ONLY because of potential RESET:
c if we want to reset (restart) CPHF after, let say, iter=3
c then using :
c
c    d3=c0*r0o + c1*r1o + c2*r2o and R3=r0o+ c0*(r1 - r0o)
c                                          + c1*(r2 - r1o)
c                                          + c2*(r3 - r2o)
c
c we express :
c
c   d5 = d3  +  c3*r3o+c4*r4o
c
c   R5 = R3  +  c3*(r4-r3o) + c4*(r5-r4o)
c
c  the above procedure is performed in turn for alpha and beta
c  quantities.
c----------------------------------------------------------------------
c
c  alpha component
c
c read in data :
c
      call read1mat(nfileA,iat,ntri*3,dA)  ! get current r0 into d
c----------------------------------------------------------------------
      liter_dir(1)=liter
      liter_dir(2)=liter
      liter_dir(3)=liter
      do istep=1,liter
         call get_fext(istep,fexi1,fexi2)
         call file4cphf_o(ntri*3,irecA,fexi1,ri,'read')   !  rl
         do icr=1,3
            call absmax(ntri,ri(1,icr),ix,rmax)
            if(rmax.le.1.0d-7) liter_dir(icr)=liter_dir(icr)-1
         enddo
      enddo
c----------------------------------------------------------------------
      do icr=1,3
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileA,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecA,fexi2,ri,'read') !  ro
            endif
c
            w(istep)=ddot(ntri,ri(1,icr),1, dA(1,icr),1) !<ro|dcon>
            c(istep)=ddot(ntri,ri(1,icr),1,ri(1,icr),1)
c
            do jstep=1,liter_dir(icr)
               call get_fext(jstep,fexj1,fexj2)
               call file4cphf_o(ntri*3,irecA,fexj1,rj,'read')    !  rl
               b(istep,jstep)=ddot(ntri,ri(1,icr),1,rj(1,icr),1) ! <ro | rl>
            enddo
         enddo
c----------------------------------------------------------
         do ii=1,liter_dir(icr)
            b(ii,ii)=b(ii,ii)-c(ii)
            w(ii)=-w(ii)
         enddo
c
         call make_smallA(b,liter,a,liter_dir(icr))
c
c        ready to solve linear system of equations liter x liter :
c
         call lineqsys(liter_dir(icr),A,w,lv,mv,c)
c
c        calculate density (residuum not needed) :
c
c         d= c1*ro_0 + c2*ro_1 + c3*ro_2 + ...+ c_l*ro_l-1
c
c         r=d0 + c1*(rl_1-ro_0)
c              + c2*(rl_2-ro_1)
c              + c3*(rl_3-ro_2)
c              +...
c              + cl*(rl_l -ro_l-1)
c
c for density
              call zeroit(dA(1,icr),ntri)
c for residuum
c             call read1mat(nfile,iat,ntri*3,ri)     ! r0
c             call dcopy(ntri,ri(1,icr),1,r(1,icr),1)
c
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileA,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecA,fexi2,ri,'read') !  ro
            endif
            call daxpy(ntri,c(istep),ri(1,icr),1,dA(1,icr),1) ! final d
c new residuum :
c           call get_fext(istep,fexj1,fexj2)
c           call file4cphf_o(ntri*3,iat,fexj1,rj,'read')    !  rl
c           call sub2vec(ntri,rj(1,icr),ri(1,icr),rj(1,icr))   ! rl_i - ro_i-1
c           call daxpy(ntri,c(istep), rj(1,icr),1, r(1,icr),1 ) ! final r
c end of new res
c
         enddo
c
      enddo
c----------------------------------------------------------------------
c
c  beta component
c
c read in data :
c
      call read1mat(nfileB,iat,ntri*3,dB)  ! get current r0 into d
c----------------------------------------------------------------------
      liter_dir(1)=liter
      liter_dir(2)=liter
      liter_dir(3)=liter
      do istep=1,liter
         call get_fext(istep,fexi1,fexi2)
         call file4cphf_o(ntri*3,irecB,fexi1,ri,'read')   !  rl
         do icr=1,3
            call absmax(ntri,ri(1,icr),ix,rmax)
            if(rmax.le.1.0d-7) liter_dir(icr)=liter_dir(icr)-1
         enddo
      enddo
c----------------------------------------------------------------------
      do icr=1,3
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileB,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecB,fexi2,ri,'read') !  ro
            endif
c
            w(istep)=ddot(ntri,ri(1,icr),1, dB(1,icr),1) !<ro|dcon>
            c(istep)=ddot(ntri,ri(1,icr),1,ri(1,icr),1)
c
            do jstep=1,liter_dir(icr)
               call get_fext(jstep,fexj1,fexj2)
               call file4cphf_o(ntri*3,irecB,fexj1,rj,'read')    !  rl
               b(istep,jstep)=ddot(ntri,ri(1,icr),1,rj(1,icr),1) ! <ro | rl>
            enddo
         enddo
c----------------------------------------------------------
         do ii=1,liter_dir(icr)
            b(ii,ii)=b(ii,ii)-c(ii)
            w(ii)=-w(ii)
         enddo
c
         call make_smallA(b,liter,a,liter_dir(icr))
c
c        ready to solve linear system of equations liter x liter :
c
         call lineqsys(liter_dir(icr),A,w,lv,mv,c)
c
c        calculate density (residuum not needed) :
c
c         d= c1*ro_0 + c2*ro_1 + c3*ro_2 + ...+ c_l*ro_l-1
c
c         r=d0 + c1*(rl_1-ro_0)
c              + c2*(rl_2-ro_1)
c              + c3*(rl_3-ro_2)
c              +...
c              + cl*(rl_l -ro_l-1)
c
c for density
              call zeroit(dB(1,icr),ntri)
c for residuum
c             call read1mat(nfile,iat,ntri*3,ri)     ! r0
c             call dcopy(ntri,ri(1,icr),1,r(1,icr),1)
c
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileB,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecB,fexi2,ri,'read') !  ro
            endif
            call daxpy(ntri,c(istep),ri(1,icr),1,dB(1,icr),1) ! final d
c new residuum :
c           call get_fext(istep,fexj1,fexj2)
c           call file4cphf_o(ntri*3,iat,fexj1,rj,'read')    !  rl
c           call sub2vec(ntri,rj(1,icr),ri(1,icr),rj(1,icr))   ! rl_i - ro_i-1
c           call daxpy(ntri,c(istep), rj(1,icr),1, r(1,icr),1 ) ! final r
c end of new res
c
         enddo
c
      enddo
      end
c======================================================================
      subroutine calc_d1r_uhf1(nfileA,nfileB,irecA,irecB,liter,ntri,iat,
     $                        a,b,w,lv,mv,c,ri,rj,dA,dB)
      implicit real*8 (a-h,o-z)
      character*4 fexi1,fexi2      ! extantion name for files
      character*4 fexj1,fexj2      ! extantion name for files
      dimension ri(ntri,3),rj(ntri,3)
      dimension  dA(ntri,3),dB(ntri,3)             ! output
c----------------------------------------------------------------------
      dimension a(liter*liter),b(liter,liter),w(liter+1),c(liter)
      dimension lv(liter), mv(liter)
      dimension liter_dir(3)
c----------------------------------------------------------------------
c input :
c
c liter   - current iteration
c ntri    - ncf*(ncf+1)/2  ; ncf=no of b.f.
c a,b     - arrays (liter x liter) for linear sys. of eq.
c   w       vector (liter)
c lv,mv   - aux. vctors needed for: osinv (b,n,det ,acc,lv,mv)
c
c           a c = w
c
c ri,rj   - storage for residuals from previous iterations
c
c output :
c
c dA,dB,  - resulting 1st-order density matrix (solution at iter=liter)
c
c
c....................................................................
c Note : the algorithm works like this :
c
c    equation to solve :  D = D0 + L(D) -> Res R= D0 + L(D) -D
c
c  iter
c
c    0      d0=dconst=r0=r0o
c    1      calc: r1=L(r0) ,
c           make r1 orthog. to r0o  r1o=Ort(r1)
c    2      calc: r2=l(r1o),
c           make r2 orthog. to r0o,r1o      r2o=Ort(r2)
c           calc d2=c0*r0o + c1*r1o
c    3      calc: r3=l(r2o),
c           make r3 orthog. to r0o,r1o,r2o  r3o=Ort(r3)
c           calc d3=c0*r0o + c1*r1o + c2*r2o
c    4      calc: r4=l(r3o),
c           make r4 orthog. to r0o,...,r3o,  r4o=Ort(r4)
c           calc d4=c0*r0o + c1*r1o + c2*r2o + c3*r3o
c    5      calc: r5=l(r4o),
c           make r5 orthog. to r0o,...,r4o,  r5o=Ort(r5)
c           calc d5=c0*r0o + c1*r1o + c2*r2o + c3*r3o+c4*r4o
c
c   k+1     calc: r(k+1)=l(rko),
c           make r(k+1) orthog. to r0o,...,rko,  r(k+1)o=Ort(r(k+1))
c           calc d(k+1)=c0*r0o + c1*r1o + ... + c(k)*r(k)o
c....................................................................
c Coefficients { c(0).....c(k) } determined from projections of the
c " predicted " residuum R(K+1) on all previous ORTHOG. resid. {r(k)o}
c
c  <r0o | R(k+1)> =0
c  <r1o | R(k+1)> =0
c  <r2o | R(k+1)> =0
c    .........
c  <rko | R(k+1)> =0
c
c where R(K+1) = r0o + c0*(r1 - r0o)
c                    + c1*(r2 - r1o)
c                    + c2*(r3 - r2o)
c                      ............
c                    + ck*(r(k+1) - rko)
c
c R(K+1) needed to be calculated ONLY because of potential RESET:
c if we want to reset (restart) CPHF after, let say, iter=3
c then using :
c
c    d3=c0*r0o + c1*r1o + c2*r2o and R3=r0o+ c0*(r1 - r0o)
c                                          + c1*(r2 - r1o)
c                                          + c2*(r3 - r2o)
c
c we express :
c
c   d5 = d3  +  c3*r3o+c4*r4o
c
c   R5 = R3  +  c3*(r4-r3o) + c4*(r5-r4o)
c
c
c read in data :
c
      call read1mat(nfileA,iat,ntri*3,dA)  ! get current r0 into d
      call read1mat(nfileB,iat,ntri*3,dB)  ! get current r0 into d
c----------------------------------------------------------------------
      liter_dir(1)=liter
      liter_dir(2)=liter
      liter_dir(3)=liter
      do istep=1,liter
         call get_fext(istep,fexi1,fexi2)
         call file4cphf_o(ntri*3,irecA,fexi1,ri,'read')   !  rl
         call file4cphf_o(ntri*3,irecB,fexi1,rj,'read')   !  rl
         do icr=1,3
            call absmax(ntri,ri(1,icr),ix,rmaxA)
            call absmax(ntri,rj(1,icr),ix,rmaxB)
            rmax=max(rmaxA,rmaxB)
            if(rmax.le.1.0d-7) liter_dir(icr)=liter_dir(icr)-1
         enddo
      enddo
c----------------------------------------------------------------------
      do icr=1,3
c
c  alpha
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileA,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecA,fexi2,ri,'read') !  ro
            endif
c
            w(istep)=ddot(ntri,ri(1,icr),1, dA(1,icr),1) !<ro|dcon>
            c(istep)=ddot(ntri,ri(1,icr),1,ri(1,icr),1)
c
            do jstep=1,liter_dir(icr)
               call get_fext(jstep,fexj1,fexj2)
               call file4cphf_o(ntri*3,irecA,fexj1,rj,'read')    !  rl
               b(istep,jstep)=ddot(ntri,ri(1,icr),1,rj(1,icr),1) ! <ro | rl>
            enddo
         enddo
c
c  beta
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileB,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecB,fexi2,ri,'read') !  ro
            endif
c
            w(istep)=w(istep)+ddot(ntri,ri(1,icr),1, dB(1,icr),1) !<ro|dcon>
            c(istep)=c(istep)+ddot(ntri,ri(1,icr),1,ri(1,icr),1)
c
            do jstep=1,liter_dir(icr)
               call get_fext(jstep,fexj1,fexj2)
               call file4cphf_o(ntri*3,irecB,fexj1,rj,'read')    !  rl
               b(istep,jstep)=b(istep,jstep)+
     $                   ddot(ntri,ri(1,icr),1,rj(1,icr),1) ! <ro | rl>
            enddo
         enddo
c----------------------------------------------------------
         do ii=1,liter_dir(icr)
            b(ii,ii)=b(ii,ii)-c(ii)
            w(ii)=-w(ii)
         enddo
c
         call make_smallA(b,liter,a,liter_dir(icr))
c
c        ready to solve linear system of equations liter x liter :
c
         call lineqsys(liter_dir(icr),A,w,lv,mv,c)
c
c        calculate density (residuum not needed) :
c
c         d= c1*ro_0 + c2*ro_1 + c3*ro_2 + ...+ c_l*ro_l-1
c
c         r=d0 + c1*(rl_1-ro_0)
c              + c2*(rl_2-ro_1)
c              + c3*(rl_3-ro_2)
c              +...
c              + cl*(rl_l -ro_l-1)
c
c for density
              call zeroit(dA(1,icr),ntri)
              call zeroit(dB(1,icr),ntri)
c for residuum
c             call read1mat(nfile,iat,ntri*3,ri)     ! r0
c             call dcopy(ntri,ri(1,icr),1,r(1,icr),1)
c
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfileA,iat,ntri*3,ri )           ! r0
               call read1mat(nfileB,iat,ntri*3,rj )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,irecA,fexi2,ri,'read') !  ro
               call file4cphf_o(ntri*3,irecB,fexi2,rj,'read') !  ro
            endif
            call daxpy(ntri,c(istep),ri(1,icr),1,dA(1,icr),1) ! final d
            call daxpy(ntri,c(istep),rj(1,icr),1,dB(1,icr),1) ! final d
c new residuum :
c           call get_fext(istep,fexj1,fexj2)
c           call file4cphf_o(ntri*3,iat,fexj1,rj,'read')    !  rl
c           call sub2vec(ntri,rj(1,icr),ri(1,icr),rj(1,icr))   ! rl_i - ro_i-1
c           call daxpy(ntri,c(istep), rj(1,icr),1, r(1,icr),1 ) ! final r
c end of new res
c
         enddo
c
      enddo
      end
c======================================================================
      subroutine make_r_orto_uhf1(nfilea,nfileb,liter,ntri,iat,
     $                      ireca,irecb,
     $                       ra,rb,ra_curr,rb_curr)
      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extantion for files
      dimension  ra_curr(ntri,3),rb_curr(ntri,3)   ! current    r
      dimension  ra(ntri,3),rb(ntri,3)             ! previous   r
      data acc,zero,one /1.d-15 , 0.d0 , 1.d0/
c
c       dn+1 = L(dn) - SUM(i=1,n)[ <di|L(dn)>/<di|di> * di ]
c
c e.g.       d2 = L(d1) - <d1|L(d1)>/<d1|d1> * d1
c
      do istep=0,liter-1
         if(istep.eq.0) then
            call read1mat(nfilea,iat,ntri*3,ra) ! very first residuum r0
            call read1mat(nfileb,iat,ntri*3,rb) ! very first residuum r0
         else
            call get_fext(istep,fext1,fext2)
            call file4cphf_o(ntri*3,ireca,fext2,ra,'read') ! previous ro
            call file4cphf_o(ntri*3,irecb,fext2,rb,'read') ! previous ro
         endif
c
c        calculate scalars <r|l(r_curr)> & <r|r>
c
         dxldx=ddot(ntri,ra(1,1),1,ra_curr(1,1),1)+
     $         ddot(ntri,rb(1,1),1,rb_curr(1,1),1)
         dyldy=ddot(ntri,ra(1,2),1,ra_curr(1,2),1)+
     $         ddot(ntri,rb(1,2),1,rb_curr(1,2),1)
         dzldz=ddot(ntri,ra(1,3),1,ra_curr(1,3),1)+
     $         ddot(ntri,rb(1,3),1,rb_curr(1,3),1)
c
         dx_dx=ddot(ntri,ra(1,1),1,ra(1,1),1)+
     $         ddot(ntri,rb(1,1),1,rb(1,1),1)
         dy_dy=ddot(ntri,ra(1,2),1,ra(1,2),1)+
     $         ddot(ntri,rb(1,2),1,rb(1,2),1)
         dz_dz=ddot(ntri,ra(1,3),1,ra(1,3),1)+
     $         ddot(ntri,rb(1,3),1,rb(1,3),1)
c
         if(dx_dx.gt.acc) then
            scx=dxldx/dx_dx
         else
            scx=zero
         endif
c
         if(dy_dy.gt.acc) then
            scy=dyldy/dy_dy
         else
            scy=zero
         endif
c
         if(dz_dz.gt.acc) then
            scz=dzldz/dz_dz
         else
            scz=zero
         endif
c
         do i=1,ntri
            ra_curr(i,1)=ra_curr(i,1)-scx*ra(i,1)
            ra_curr(i,2)=ra_curr(i,2)-scy*ra(i,2)
            ra_curr(i,3)=ra_curr(i,3)-scz*ra(i,3)
            rb_curr(i,1)=rb_curr(i,1)-scx*rb(i,1)
            rb_curr(i,2)=rb_curr(i,2)-scy*rb(i,2)
            rb_curr(i,3)=rb_curr(i,3)-scz*rb(i,3)
         enddo
      enddo
c
      end
c======================================================================
      subroutine make_screenD_uhf(iscre,d1a, d1b,inx,lind,
     *                            map_fs, natend,natonce,ncs,ncf,ntri)
      implicit real*8 (a-h,o-z)
      dimension d1A(ntri*3,natonce),d1B(ntri*3,natonce)
      dimension natend(natonce)
      dimension inx(12,*), lind(*) , map_fs(*)
c
c     take current D1const  from unit=65,75
c
c     alpha D1 const on iuniA=65
c     beta  D1 const on iuniB=75
c
      iuniA=65
      iuniB=75
c
         do iat=1,natonce
            call read1mat(iuniA,iat,ntri*3,d1A(1,iat))
            call read1mat(iuniB,iat,ntri*3,d1B(1,iat))
            call add2vabs(ntri*3*natonce,d1A,d1B,d1A)
         enddo
c
         call trspmo(d1A,ntri, d1B,3*natonce)
         call dcopy(3*ntri*natonce,d1B, 1 , d1A, 1)
c
         call select_dnat(natend,natonce,d1A,ntri,d1B)
c
         call setup_densp2(inx,ncf,ncs,d1B,d1A,map_fs)
c           screening density NOT squared
c
         call save1mat(69, 2 ,ncs*ncs,d1A)
c
         end
c======================================================================
c======================================================================
      subroutine add2vabs(ndim,A,B,C)
      implicit real*8 (a-h,o-z)
      dimension a(ndim),b(ndim),c(ndim)
      do i=1,ndim
        c(i)=abs(a(i))+abs(b(i))
      enddo
      end
c==============================================================
      subroutine nd1g0_core_nos_uhf(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dena,        denb,       focka,      fockb,
     $        ntri,        nat,        idft,       ax,
     $        densp,       ncs,         map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dena(3*nat,ntri),fockA(3*nat,ntri)
      dimension denB(3*nat,ntri),fockB(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
              ij=lind(max0(icf,jcf))+min0(icf,jcf)
              do kcf=kcff+1,kcff+klen
                ik=lind(max0(icf,kcf))+min0(icf,kcf)
                jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
                do lcf=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint2=xint0+xint0
                    if(idft.gt.0) xint0=xint0*ax

                    il=lind(max0(icf,lcf))+min0(icf,lcf)
                    jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                    kl=lind(max0(kcf,lcf))+min0(kcf,lcf)

c ......................Coulomb  part ....................
                    do iat=1,nat3
                      focka(iat,ij)=focka(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      fockb(iat,ij)=fockb(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      focka(iat,kl)=focka(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                      fockb(iat,kl)=fockb(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                    enddo
c ......................Exchange part ....................
                    call daxpy(nat3,-xint0,dena(1,jl),1,focka(1,ik),1)
                    call daxpy(nat3,-xint0,dena(1,ik),1,focka(1,jl),1)
                    call daxpy(nat3,-xint0,dena(1,il),1,focka(1,jk),1)
                    call daxpy(nat3,-xint0,dena(1,jk),1,focka(1,il),1)
                    call daxpy(nat3,-xint0,denb(1,jl),1,fockb(1,ik),1)
                    call daxpy(nat3,-xint0,denb(1,ik),1,fockb(1,jl),1)
                    call daxpy(nat3,-xint0,denb(1,il),1,fockb(1,jk),1)
                    call daxpy(nat3,-xint0,denb(1,jk),1,fockb(1,il),1)
c..........................................................
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_sym_uhf(
     $        nsym,        ifp,        core,        ilab,
     $        jlab,       klab,        llab,        isiz,
     $        jsiz,       ksiz,        lsiz,        iqrt,
     $        integrals,  nqstore,     nblstore,    thres0,
     $        thres1,     lind,        dena,        denb,
     $        focka,      fockb,       ntri,        nat,
     $        idft,       ax,          densp,       ncs,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do ic1=icff+1,icff+ilen
            do jc1=jcff+1,jcff+jlen
              ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
              do kc1=kcff+1,kcff+klen
                ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
                jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
                do lc1=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint2=xint0+xint0
                    if(idft.gt.0) xint0=xint0*ax

                    il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                    jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                    kl=lind(max0(kc1,lc1))+min0(kc1,lc1)

c ......................Coulomb  part ....................
                    do iat=1,nat3
                      focka(iat,ij)=focka(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      fockb(iat,ij)=fockb(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      focka(iat,kl)=focka(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                      fockb(iat,kl)=fockb(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                    enddo
c ......................Exchange part ....................
                    call daxpy(nat3,-xint0,dena(1,jl),1,focka(1,ik),1)
                    call daxpy(nat3,-xint0,dena(1,ik),1,focka(1,jl),1)
                    call daxpy(nat3,-xint0,dena(1,il),1,focka(1,jk),1)
                    call daxpy(nat3,-xint0,dena(1,jk),1,focka(1,il),1)
                    call daxpy(nat3,-xint0,denb(1,jl),1,fockb(1,ik),1)
                    call daxpy(nat3,-xint0,denb(1,ik),1,fockb(1,jl),1)
                    call daxpy(nat3,-xint0,denb(1,il),1,fockb(1,jk),1)
                    call daxpy(nat3,-xint0,denb(1,jk),1,fockb(1,il),1)
c..........................................................
                    do ns=1,nsym
                      xco=xint2
                      xex=xint0
                      icf=ifp(ns,ic1)
                      jcf=ifp(ns,jc1)
                      kcf=ifp(ns,kc1)
                      lcf=ifp(ns,lc1)
                      if(icf.lt.0) then
                        icf=-icf
                        xco=-xco
                        xex=-xex
                      endif
                      if(jcf.lt.0) then
                        jcf=-jcf
                        xco=-xco
                        xex=-xex
                      endif
                      if(kcf.lt.0) then
                        kcf=-kcf
                        xco=-xco
                        xex=-xex
                      endif
                      if(lcf.lt.0) then
                        lcf=-lcf
                        xco=-xco
                        xex=-xex
                      endif
                      lic=lind(icf)
                      ljc=lind(jcf)
                      lkc=lind(kcf)
                      llc=lind(lcf)
                      ijc=max0(lic,ljc)+min0(icf,jcf)
                      ikc=max0(lic,lkc)+min0(icf,kcf)
                      ilc=max0(lic,llc)+min0(icf,lcf)
                      jkc=max0(ljc,lkc)+min0(jcf,kcf)
                      jlc=max0(ljc,llc)+min0(jcf,lcf)
                      klc=max0(lkc,llc)+min0(kcf,lcf)
c ......................Coulomb  part ....................
                      do iat=1,nat3
                        focka(iat,ijc)=focka(iat,ijc)+
     $                       (dena(iat,klc)+denb(iat,klc))*xco
                        fockb(iat,ijc)=fockb(iat,ijc)+
     $                       (dena(iat,klc)+denb(iat,klc))*xco
                        focka(iat,klc)=focka(iat,klc)+
     $                       (dena(iat,ijc)+denb(iat,ijc))*xco
                        fockb(iat,klc)=fockb(iat,klc)+
     $                       (dena(iat,ijc)+denb(iat,ijc))*xco
                      enddo
c ......................Exchange part ....................
                      call daxpy(nat3,-xex,dena(1,jlc),1,focka(1,ikc),1)
                      call daxpy(nat3,-xex,dena(1,ikc),1,focka(1,jlc),1)
                      call daxpy(nat3,-xex,dena(1,ilc),1,focka(1,jkc),1)
                      call daxpy(nat3,-xex,dena(1,jkc),1,focka(1,ilc),1)
                      call daxpy(nat3,-xex,denb(1,jlc),1,fockb(1,ikc),1)
                      call daxpy(nat3,-xex,denb(1,ikc),1,fockb(1,jlc),1)
                      call daxpy(nat3,-xex,denb(1,ilc),1,fockb(1,jkc),1)
                      call daxpy(nat3,-xex,denb(1,jkc),1,fockb(1,ilc),1)
c..........................................................
                    enddo
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_nosC_uhf(
     $        core,        ilab,       jlab,       klab,
     $        llab,        isiz,       jsiz,       ksiz,
     $        lsiz,        iqrt,       integrals,  nqstore,
     $        nblstore,    thres0,     thres1,     lind,
     $        dena,        denb,       focka,      fockb,
     $        ntri,        nat,        densp,      ncs,
     $      map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, NO Symmetry, Coulomb term only
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do icf=icff+1,icff+ilen
            do jcf=jcff+1,jcff+jlen
              ij=lind(max0(icf,jcf))+min0(icf,jcf)
              do kcf=kcff+1,kcff+klen
                do lcf=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint2=xint0+xint0

                    kl=lind(max0(kcf,lcf))+min0(kcf,lcf)

c ......................Coulomb  part ....................
                    do iat=1,nat3
                      focka(iat,ij)=focka(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      fockb(iat,ij)=fockb(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      focka(iat,kl)=focka(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                      fockb(iat,kl)=fockb(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                    enddo
c..........................................................
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c==============================================================
      subroutine nd1g0_core_symC_uhf(
     $        nsym,        ifp,        core,        ilab,
     $        jlab,       klab,        llab,        isiz,
     $        jsiz,       ksiz,        lsiz,        iqrt,
     $        integrals,  nqstore,     nblstore,    thres0,
     $        thres1,     lind,        dena,        denb,
     $        focka,      fockb,       ntri,        nat,
     $        densp,      ncs,         map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c in core storage, Symmetry, Coulomb term only
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension ifp(7,*)             !    ncf)
      integer*2 ilab(nqstore),jlab(nqstore),klab(nqstore),llab(nqstore)
      integer*2 isiz(*),jsiz(*),ksiz(*),lsiz(*)   ! dim=nblstore
      integer*4 iqrt(*)                           ! dim=nblstore
c
      dimension core(integrals)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
      dimension lind(*)
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c
c----------------------------------------------------------------
c integrals in core have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescaled the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
      intx=0
      ijklp=0
      do isbl=1,nblstore
        ilen=isiz(isbl)
        jlen=jsiz(isbl)
        klen=ksiz(isbl)
        llen=lsiz(isbl)
c
        nqrt1=iqrt(isbl)
c
        do iqrt1=1,nqrt1
          ijklp=ijklp+1
          icff=ilab(ijklp)
          jcff=jlab(ijklp)
          kcff=klab(ijklp)
          lcff=llab(ijklp)

          ics=map_fs(icff+1)
          jcs=map_fs(jcff+1)
          kcs=map_fs(kcff+1)
          lcs=map_fs(lcff+1)
c
c screening density
c
          dij=densp(ics,jcs)
          dkl=densp(kcs,lcs)
          dik=densp(ics,kcs)
          dil=densp(ics,lcs)
          djk=densp(jcs,kcs)
          djl=densp(jcs,lcs)
c
          dmax=max(4.d0*dij,4.d0*dkl,dik,dil,djk,djl)
c------------------------------------------------------
c shell 
          do ic1=icff+1,icff+ilen
            do jc1=jcff+1,jcff+jlen
              ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
              do kc1=kcff+1,kcff+klen
                do lc1=lcff+1,lcff+llen
                  intx=intx+1
                  xint0=core(intx)
                  if(dmax*abs(xint0).gt.0.5d0) then
                    xint2=xint0+xint0

                    kl=lind(max0(kc1,lc1))+min0(kc1,lc1)

c ......................Coulomb  part ....................
                    do iat=1,nat3
                      focka(iat,ij)=focka(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      fockb(iat,ij)=fockb(iat,ij)+
     $                     (dena(iat,kl)+denb(iat,kl))*xint2
                      focka(iat,kl)=focka(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                      fockb(iat,kl)=fockb(iat,kl)+
     $                     (dena(iat,ij)+denb(iat,ij))*xint2
                    enddo
c..........................................................
                    do ns=1,nsym
                      xco=xint2
                      icf=ifp(ns,ic1)
                      jcf=ifp(ns,jc1)
                      kcf=ifp(ns,kc1)
                      lcf=ifp(ns,lc1)
                      if(icf.lt.0) then
                        icf=-icf
                        xco=-xco
                      endif
                      if(jcf.lt.0) then
                        jcf=-jcf
                        xco=-xco
                      endif
                      if(kcf.lt.0) then
                        kcf=-kcf
                        xco=-xco
                      endif
                      if(lcf.lt.0) then
                        lcf=-lcf
                        xco=-xco
                      endif
                      ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                      klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ......................Coulomb  part ....................
                      do iat=1,nat3
                        focka(iat,ijc)=focka(iat,ijc)+
     $                       (dena(iat,klc)+denb(iat,klc))*xco
                        fockb(iat,ijc)=fockb(iat,ijc)+
     $                       (dena(iat,klc)+denb(iat,klc))*xco
                        focka(iat,klc)=focka(iat,klc)+
     $                       (dena(iat,ijc)+denb(iat,ijc))*xco
                        fockb(iat,klc)=fockb(iat,klc)+
     $                       (dena(iat,ijc)+denb(iat,ijc))*xco
                      enddo
c..........................................................
                    enddo
                  endif
                enddo         !   do lll=1,llen
              enddo         !   do kkk=1,klen
            enddo         !   do jjj=1,jlen
          enddo         !   do iii=1,ilen
        enddo         !   do iqrt1=1,nqrt1
      enddo         !   do isbl=1,nblstore
c
c-----------------------------------------------------------
c
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_nos_uhf(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dena,       denb,
     $        focka,       fockb,      ntri,       nat,
     $        idft,        ax,         ncs,        densp,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical rescale
      logical doint72,doint62
c----------------------------------------------------------------
      common /screening/ dens_max_el 
c----------------------------------------------------------------
      integer*2 ijkllab(4,iqstore) !   ilab,jlab,klab,llab
      integer*2 ijklsiz(4,iblstore)!   isiz(isbl),jsiz,ksiz,lsiz  
      integer*2 nquarts(iblstore)     !   no of quartets in each block
c----------------------------------------------------------------
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
c integrals on-disk have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
c units for stored stuff :
c
      iunit0=160 ! for no of stored bl & qr
      iunit1=161 ! sizes
      iunit2=162 ! quarts
      iunit3=163 ! labels
      iunit4=164 ! integrals in  real*8 & i*4 & i*2
      iunit6=166 ! integrals in  i*2
      iunit7=167 ! integrals in  i*2
c----------------------------------------------------------------
c maximum real value of integrals stored in xinteg2() 
c
      xmax7=0.5d-7
      xmax6=0.5d-6
c 
c----------------------------------------------------------------
c kept as xmax7/thres0 ! e.g. 0.5*10**-7/10**-10= 500
c
c if actual threshold thres1 is bigger then or equal to 10**-7
c then we do NOT need these small integrals. They are needed 
c only if thres1 < 10**-7. 
c Integrlas have to be 
c
c         x > 0.5*thres1
c Including maximum density d gives 
c      
c      d* x > 0.5*thres1
c
c Thus we have :
c
c   d* 0.5 *10**-7 > d*x > 0.5*thres1
c
c which finally yields
c
c   d* 10**-7 > thres1
c
c----------------------------------------------------------------
c Thresholds ratio is
c
      t01=thres0/thres1
c----------------------------------------------------------------
      dens_max=4.0d0*sqrt(dens_max_el)
      dens_max=min(dens_max,1.0d0)
c
      doint72=.false.
      if(dens_max*1.0d-7.GT.thres1) doint72=.true.
      doint62=.false.
      if(dens_max*1.0d-6.GT.thres1) doint62=.true.
c----------------------------------------------------------------
c     write(6,666) dens_max,thres1,doint72,doint62
c666  format('   dens_max=',f10.6,' thers1=',1pe8.1' doint=',2(L2,1x))
c----------------------------------------------------------------
      call secund(r1)
      call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
      call secund(r2)
      call elapsec(e2)
      readcpu3=readcpu3+(r2-r1)
      readela3=readela3+(e2-e1)
c----------------------------------------------------------------
      ijklp=0
      do isplit=1,nsplit
        isbl_beg=ifromto(1,isplit)
        isbl_end=ifromto(2,isplit)
        integi72=nofinteg(1,isplit)
        integi62=nofinteg(2,isplit)
        integin2=nofinteg(3,isplit)
        integin4=nofinteg(4,isplit)
        integin8=nofinteg(5,isplit)
c
        call secund(r1)
        call elapsec(e1)
c
        if(doint72) then
        if(integi72.gt.0) then
          call read_int2(iunit6,xinte72,integi72)
        endif
c
        endif
        if(doint62) then
        if(integi62.gt.0) then
          call read_int2(iunit7,xinte62,integi62)
        endif
        endif
c
        if(integin2.gt.0) then
          call read_int2(iunit4,xinteg2,integin2)
        endif
c
        if(integin4.gt.0) then
          call read_int4(iunit4,xinteg4,integin4)
        endif
c
        if(integin8.gt.0) then
          call read_int8(iunit4,xinteg8,integin8)
        endif
c
        call secund(r2)
        call elapsec(e2)
        readcpu4=readcpu4+(r2-r1)
        readela4=readela4+(e2-e1)
c
        intx8=0
        intx4=0
        intx2=0
        int62=0
        int72=0
        do isbl=isbl_beg,isbl_end
          ilen=ijklsiz(1,isbl)
          jlen=ijklsiz(2,isbl)
          klen=ijklsiz(3,isbl)
          llen=ijklsiz(4,isbl)
          nqrt=nquarts(isbl)
          isize=ilen*jlen*klen*llen
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
            icase=5
            if(integin8.gt.0) then
              if(icff.lt.0) then
                icase=1
                icff=-icff
              endif
            endif
c
            if(integin4.gt.0.and.icase.eq.5) then
              if(jcff.lt.0) then
                icase=2
                jcff=-jcff
              endif
            endif
c
            if(integin2.gt.0.and.icase.eq.5) then
              if(kcff.lt.0) then
                icase=3
                kcff=-kcff
              endif
            endif
c
            if(integi62.gt.0.and.icase.eq.5) then
              if(lcff.lt.0) then
                icase=4
                lcff=-lcff
              endif
            endif
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_nos_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        intx8,      xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_nos_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        intx4,      xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_nos_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        intx2,      xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_nos_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       idft,      ax,
     $             dmax,        int62,      xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                  call make_nfock2_nos_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       idft,      ax,
     $             dmax,        int72,      xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_nosC_uhf(
     $        ncf,         xinte72,    xinte62,    xinteg2,
     $        xinteg4,     xinteg8,    ijkllab,    ijklsiz,
     $        nquarts,     iqstore,    iblstore,   rescale,
     $        nsplit,      ifromto,    nofinteg,   thres0,
     $        thres1,      lind,       dena,       denb,
     $        focka,       fockb,      ntri,       nat,
     $        ncs,        densp,       map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry, Coulomb term only
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical rescale
      logical doint72,doint62
c----------------------------------------------------------------
      common /screening/ dens_max_el 
c----------------------------------------------------------------
      integer*2 ijkllab(4,iqstore) !   ilab,jlab,klab,llab
      integer*2 ijklsiz(4,iblstore)!   isiz(isbl),jsiz,ksiz,lsiz  
      integer*2 nquarts(iblstore)     !   no of quartets in each block
c----------------------------------------------------------------
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
c
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
c integrals on-disk have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
c units for stored stuff :
c
      iunit0=160 ! for no of stored bl & qr
      iunit1=161 ! sizes
      iunit2=162 ! quarts
      iunit3=163 ! labels
      iunit4=164 ! integrals in  real*8 & i*4 & i*2
      iunit6=166 ! integrals in  i*2
      iunit7=167 ! integrals in  i*2
c----------------------------------------------------------------
c maximum real value of integrals stored in xinteg2() 
c
      xmax7=0.5d-7
      xmax6=0.5d-6
c 
c----------------------------------------------------------------
c kept as xmax7/thres0 ! e.g. 0.5*10**-7/10**-10= 500
c
c if actual threshold thres1 is bigger then or equal to 10**-7
c then we do NOT need these small integrals. They are needed 
c only if thres1 < 10**-7. 
c Integrlas have to be 
c
c         x > 0.5*thres1
c Including maximum density d gives 
c      
c      d* x > 0.5*thres1
c
c Thus we have :
c
c   d* 0.5 *10**-7 > d*x > 0.5*thres1
c
c which finally yields
c
c   d* 10**-7 > thres1
c
c----------------------------------------------------------------
c Thresholds ratio is
c
      t01=thres0/thres1
c----------------------------------------------------------------
      dens_max=4.0d0*sqrt(dens_max_el)
      dens_max=min(dens_max,1.0d0)
c
      doint72=.false.
      if(dens_max*1.0d-7.GT.thres1) doint72=.true.
      doint62=.false.
      if(dens_max*1.0d-6.GT.thres1) doint62=.true.
c----------------------------------------------------------------
c     write(6,666) dens_max,thres1,doint72,doint62
c666  format('   dens_max=',f10.6,' thers1=',1pe8.1' doint=',2(L2,1x))
c----------------------------------------------------------------
      call secund(r1)
      call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
      call secund(r2)
      call elapsec(e2)
      readcpu3=readcpu3+(r2-r1)
      readela3=readela3+(e2-e1)
c----------------------------------------------------------------
      ijklp=0
      do isplit=1,nsplit
        isbl_beg=ifromto(1,isplit)
        isbl_end=ifromto(2,isplit)
        integi72=nofinteg(1,isplit)
        integi62=nofinteg(2,isplit)
        integin2=nofinteg(3,isplit)
        integin4=nofinteg(4,isplit)
        integin8=nofinteg(5,isplit)
c
        call secund(r1)
        call elapsec(e1)
c
        if(doint72) then
        if(integi72.gt.0) then
          call read_int2(iunit6,xinte72,integi72)
        endif
c
        endif
        if(doint62) then
        if(integi62.gt.0) then
          call read_int2(iunit7,xinte62,integi62)
        endif
        endif
c
        if(integin2.gt.0) then
          call read_int2(iunit4,xinteg2,integin2)
        endif
c
        if(integin4.gt.0) then
          call read_int4(iunit4,xinteg4,integin4)
        endif
c
        if(integin8.gt.0) then
          call read_int8(iunit4,xinteg8,integin8)
        endif
c
        call secund(r2)
        call elapsec(e2)
        readcpu4=readcpu4+(r2-r1)
        readela4=readela4+(e2-e1)
c
        intx8=0
        intx4=0
        intx2=0
        int62=0
        int72=0
        do isbl=isbl_beg,isbl_end
          ilen=ijklsiz(1,isbl)
          jlen=ijklsiz(2,isbl)
          klen=ijklsiz(3,isbl)
          llen=ijklsiz(4,isbl)
          nqrt=nquarts(isbl)
          isize=ilen*jlen*klen*llen
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
            icase=5
            if(integin8.gt.0) then
              if(icff.lt.0) then
                icase=1
                icff=-icff
              endif
            endif
c
            if(integin4.gt.0.and.icase.eq.5) then
              if(jcff.lt.0) then
                icase=2
                jcff=-jcff
              endif
            endif
c
            if(integin2.gt.0.and.icase.eq.5) then
              if(kcff.lt.0) then
                icase=3
                kcff=-kcff
              endif
            endif
c
            if(integi62.gt.0.and.icase.eq.5) then
              if(lcff.lt.0) then
                icase=4
                lcff=-lcff
              endif
            endif
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_nosC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      intx8,
     $         xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_nosC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      intx4,
     $         xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_nosC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      intx2,
     $         xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_nosC_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       dmax,      int62,
     $             xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_nosC_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       dmax,      int72,
     $             xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_sym_uhf(
     $        nsym,        ifp,        ncf,         xinte72,
     $        xinte62,     xinteg2,    xinteg4,     xinteg8,
     $        ijkllab,     ijklsiz,    nquarts,     iqstore,
     $        iblstore,    rescale,    nsplit,      ifromto,
     $        nofinteg,    thres0,     thres1,      lind,
     $        dena,        denb,       focka,       fockb,
     $        ntri,        nat,        idft,        ax,
     $        ncs,         densp,      map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical rescale
      logical doint72,doint62
c----------------------------------------------------------------
      common /screening/ dens_max_el 
c----------------------------------------------------------------
      integer*2 ijkllab(4,iqstore) !   ilab,jlab,klab,llab
      integer*2 ijklsiz(4,iblstore)!   isiz(isbl),jsiz,ksiz,lsiz  
      integer*2 nquarts(iblstore)     !   no of quartets in each block
c----------------------------------------------------------------
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
      dimension ifp(7,*)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
c integrals on-disk have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
c units for stored stuff :
c
      iunit0=160 ! for no of stored bl & qr
      iunit1=161 ! sizes
      iunit2=162 ! quarts
      iunit3=163 ! labels
      iunit4=164 ! integrals in  real*8 & i*4 & i*2
      iunit6=166 ! integrals in  i*2
      iunit7=167 ! integrals in  i*2
c----------------------------------------------------------------
c maximum real value of integrals stored in xinteg2() 
c
      xmax7=0.5d-7
      xmax6=0.5d-6
c 
c----------------------------------------------------------------
c kept as xmax7/thres0 ! e.g. 0.5*10**-7/10**-10= 500
c
c if actual threshold thres1 is bigger then or equal to 10**-7
c then we do NOT need these small integrals. They are needed 
c only if thres1 < 10**-7. 
c Integrlas have to be 
c
c         x > 0.5*thres1
c Including maximum density d gives 
c      
c      d* x > 0.5*thres1
c
c Thus we have :
c
c   d* 0.5 *10**-7 > d*x > 0.5*thres1
c
c which finally yields
c
c   d* 10**-7 > thres1
c
c----------------------------------------------------------------
c Thresholds ratio is
c
      t01=thres0/thres1
c----------------------------------------------------------------
      dens_max=4.0d0*sqrt(dens_max_el)
      dens_max=min(dens_max,1.0d0)
c
      doint72=.false.
      if(dens_max*1.0d-7.GT.thres1) doint72=.true.
      doint62=.false.
      if(dens_max*1.0d-6.GT.thres1) doint62=.true.
c----------------------------------------------------------------
c     write(6,666) dens_max,thres1,doint72,doint62
c666  format('   dens_max=',f10.6,' thers1=',1pe8.1' doint=',2(L2,1x))
c----------------------------------------------------------------
      call secund(r1)
      call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
      call secund(r2)
      call elapsec(e2)
      readcpu3=readcpu3+(r2-r1)
      readela3=readela3+(e2-e1)
c----------------------------------------------------------------
      ijklp=0
      do isplit=1,nsplit
        isbl_beg=ifromto(1,isplit)
        isbl_end=ifromto(2,isplit)
        integi72=nofinteg(1,isplit)
        integi62=nofinteg(2,isplit)
        integin2=nofinteg(3,isplit)
        integin4=nofinteg(4,isplit)
        integin8=nofinteg(5,isplit)
c
        call secund(r1)
        call elapsec(e1)
c
        if(doint72) then
        if(integi72.gt.0) then
          call read_int2(iunit6,xinte72,integi72)
        endif
c
        endif
        if(doint62) then
        if(integi62.gt.0) then
          call read_int2(iunit7,xinte62,integi62)
        endif
        endif
c
        if(integin2.gt.0) then
          call read_int2(iunit4,xinteg2,integin2)
        endif
c
        if(integin4.gt.0) then
          call read_int4(iunit4,xinteg4,integin4)
        endif
c
        if(integin8.gt.0) then
          call read_int8(iunit4,xinteg8,integin8)
        endif
c
        call secund(r2)
        call elapsec(e2)
        readcpu4=readcpu4+(r2-r1)
        readela4=readela4+(e2-e1)
c
        intx8=0
        intx4=0
        intx2=0
        int62=0
        int72=0
        do isbl=isbl_beg,isbl_end
          ilen=ijklsiz(1,isbl)
          jlen=ijklsiz(2,isbl)
          klen=ijklsiz(3,isbl)
          llen=ijklsiz(4,isbl)
          nqrt=nquarts(isbl)
          isize=ilen*jlen*klen*llen
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
            icase=5
            if(integin8.gt.0) then
              if(icff.lt.0) then
                icase=1
                icff=-icff
              endif
            endif
c
            if(integin4.gt.0.and.icase.eq.5) then
              if(jcff.lt.0) then
                icase=2
                jcff=-jcff
              endif
            endif
c
            if(integin2.gt.0.and.icase.eq.5) then
              if(kcff.lt.0) then
                icase=3
                kcff=-kcff
              endif
            endif
c
            if(integi62.gt.0.and.icase.eq.5) then
              if(lcff.lt.0) then
                icase=4
                lcff=-lcff
              endif
            endif
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_sym_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        nsym,       ifp,       intx8,
     $         xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_sym_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        nsym,       ifp,       intx4,
     $         xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_sym_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       idft,      ax,
     $         dmax,        nsym,       ifp,       intx2,
     $         xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_sym_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       idft,      ax,
     $             dmax,        nsym,       ifp,       int62,
     $             xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_sym_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       idft,      ax,
     $             dmax,        nsym,       ifp,       int72,
     $             xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c=====================================================================
      subroutine nd1g0_disk_symC_uhf(
     $        nsym,        ifp,        ncf,         xinte72,
     $        xinte62,     xinteg2,    xinteg4,     xinteg8,
     $        ijkllab,     ijklsiz,    nquarts,     iqstore,
     $        iblstore,    rescale,    nsplit,      ifromto,
     $        nofinteg,    thres0,     thres1,      lind,
     $        dena,        denb,       focka,       fockb,
     $        ntri,        nat,        ncs,         densp,
     $        map_fs)
c----------------------------------------------------------------
c Fock builder ( G(D1,g0) matrix for hessian CPHF )
c
c disk storage, NO Symmetry
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical rescale
      logical doint72,doint62
c----------------------------------------------------------------
      common /screening/ dens_max_el 
c----------------------------------------------------------------
      integer*2 ijkllab(4,iqstore) !   ilab,jlab,klab,llab
      integer*2 ijklsiz(4,iblstore)!   isiz(isbl),jsiz,ksiz,lsiz  
      integer*2 nquarts(iblstore)     !   no of quartets in each block
c----------------------------------------------------------------
      integer*2 xinte72(*)
      integer*2 xinte62(*)
      integer*2 xinteg2(*)
      integer*4 xinteg4(*)
      dimension xinteg8(*)
      dimension ifp(7,*)
      dimension dena(3*nat,ntri),focka(3*nat,ntri)
      dimension denb(3*nat,ntri),fockb(3*nat,ntri)
c
      dimension lind(*)
c----------------------------------------------------------------
      dimension densp(ncs,ncs)
      dimension map_fs(*)
c----------------------------------------------------------------
      dimension nofinteg(5,nsplit),ifromto(2,nsplit)
c----------------------------------------------------------------
      common /readwrite/ readcpu1,readcpu2,readcpu3,readcpu4,
     *                   readela1,readela2,readela3,readela4
c----------------------------------------------------------------
      data half /0.5d0/
c----------------------------------------------------------------
c integrals on-disk have been calculated with threshold thres0
c while for re-calculated integrals thres1 is used :
c we need to rescale the Fock matrix by resc=thres0/thres1
c----------------------------------------------------------------
      nat3=nat*3
c units for stored stuff :
c
      iunit0=160 ! for no of stored bl & qr
      iunit1=161 ! sizes
      iunit2=162 ! quarts
      iunit3=163 ! labels
      iunit4=164 ! integrals in  real*8 & i*4 & i*2
      iunit6=166 ! integrals in  i*2
      iunit7=167 ! integrals in  i*2
c----------------------------------------------------------------
c maximum real value of integrals stored in xinteg2() 
c
      xmax7=0.5d-7
      xmax6=0.5d-6
c 
c----------------------------------------------------------------
c kept as xmax7/thres0 ! e.g. 0.5*10**-7/10**-10= 500
c
c if actual threshold thres1 is bigger then or equal to 10**-7
c then we do NOT need these small integrals. They are needed 
c only if thres1 < 10**-7. 
c Integrlas have to be 
c
c         x > 0.5*thres1
c Including maximum density d gives 
c      
c      d* x > 0.5*thres1
c
c Thus we have :
c
c   d* 0.5 *10**-7 > d*x > 0.5*thres1
c
c which finally yields
c
c   d* 10**-7 > thres1
c
c----------------------------------------------------------------
c Thresholds ratio is
c
      t01=thres0/thres1
c----------------------------------------------------------------
      dens_max=4.0d0*sqrt(dens_max_el)
      dens_max=min(dens_max,1.0d0)
c
      doint72=.false.
      if(dens_max*1.0d-7.GT.thres1) doint72=.true.
      doint62=.false.
      if(dens_max*1.0d-6.GT.thres1) doint62=.true.
c----------------------------------------------------------------
c     write(6,666) dens_max,thres1,doint72,doint62
c666  format('   dens_max=',f10.6,' thers1=',1pe8.1' doint=',2(L2,1x))
c----------------------------------------------------------------
      call secund(r1)
      call elapsec(e1)
c
      call read_int2(iunit1,ijklsiz,4*iblstore)
      call read_int2(iunit2,nquarts,iblstore)
      call read_int2(iunit3,ijkllab,4*iqstore)
c
      call secund(r2)
      call elapsec(e2)
      readcpu3=readcpu3+(r2-r1)
      readela3=readela3+(e2-e1)
c----------------------------------------------------------------
      ijklp=0
      do isplit=1,nsplit
        isbl_beg=ifromto(1,isplit)
        isbl_end=ifromto(2,isplit)
        integi72=nofinteg(1,isplit)
        integi62=nofinteg(2,isplit)
        integin2=nofinteg(3,isplit)
        integin4=nofinteg(4,isplit)
        integin8=nofinteg(5,isplit)
c
        call secund(r1)
        call elapsec(e1)
c
        if(doint72) then
        if(integi72.gt.0) then
          call read_int2(iunit6,xinte72,integi72)
        endif
c
        endif
        if(doint62) then
        if(integi62.gt.0) then
          call read_int2(iunit7,xinte62,integi62)
        endif
        endif
c
        if(integin2.gt.0) then
          call read_int2(iunit4,xinteg2,integin2)
        endif
c
        if(integin4.gt.0) then
          call read_int4(iunit4,xinteg4,integin4)
        endif
c
        if(integin8.gt.0) then
          call read_int8(iunit4,xinteg8,integin8)
        endif
c
        call secund(r2)
        call elapsec(e2)
        readcpu4=readcpu4+(r2-r1)
        readela4=readela4+(e2-e1)
c
        intx8=0
        intx4=0
        intx2=0
        int62=0
        int72=0
        do isbl=isbl_beg,isbl_end
          ilen=ijklsiz(1,isbl)
          jlen=ijklsiz(2,isbl)
          klen=ijklsiz(3,isbl)
          llen=ijklsiz(4,isbl)
          nqrt=nquarts(isbl)
          isize=ilen*jlen*klen*llen
          do iqrt=1,nqrt
            ijklp=ijklp+1
            icff=ijkllab(1,ijklp)
            jcff=ijkllab(2,ijklp)
            kcff=ijkllab(3,ijklp)
            lcff=ijkllab(4,ijklp)
            icase=5
            if(integin8.gt.0) then
              if(icff.lt.0) then
                icase=1
                icff=-icff
              endif
            endif
c
            if(integin4.gt.0.and.icase.eq.5) then
              if(jcff.lt.0) then
                icase=2
                jcff=-jcff
              endif
            endif
c
            if(integin2.gt.0.and.icase.eq.5) then
              if(kcff.lt.0) then
                icase=3
                kcff=-kcff
              endif
            endif
c
            if(integi62.gt.0.and.icase.eq.5) then
              if(lcff.lt.0) then
                icase=4
                lcff=-lcff
              endif
            endif
c------------------------------------------------------
            icff=icff-1
            jcff=jcff-1
            kcff=kcff-1
            lcff=lcff-1
            ics=map_fs(icff+1)
            jcs=map_fs(jcff+1)
            kcs=map_fs(kcff+1)
            lcs=map_fs(lcff+1)
c
c screening density
c
            dij=densp(ics,jcs)
            dkl=densp(kcs,lcs)
            dik=densp(ics,kcs)
            dil=densp(ics,lcs)
            djk=densp(jcs,kcs)
            djl=densp(jcs,lcs)
c
            dmax=max(4.0d0*dij,4.0d0*dkl,dik,dil,djk,djl)
c
            if(integin8.gt.intx8 .and. icase.eq.1) then
              call make_nfock8_symC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      nsym,
     $         ifp,         intx8,      xinteg8)
            endif
c
            if(integin4.gt.intx4 .and. icase.eq.2) then
              call make_nfock4_symC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      nsym,
     $         ifp,         intx4,      xinteg4)
            endif
c
            if(integin2.gt.intx2 .and. icase.eq.3) then
              call make_nfock2_symC_uhf(
     $         ilen,        jlen,       klen,      llen,
     $         icff,        jcff,       kcff,      lcff,
     $         dena,        denb,       focka,     fockb,
     $         lind,        nat3,       dmax,      nsym,
     $         ifp,         intx2,      xinteg2)
            endif
c
            if(doint62) then
              if(integi62.gt.int62 .and. icase.eq.4) then
                if(dmax1*1.0d-6.gt.thres1) then
                  call make_nfock2_symC_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       dmax,      nsym,
     $             ifp,         int62,      xinte62)
                endif
              endif
            endif
c
            if(doint72) then
              if(integi72.gt.int72 .and. icase.eq.5) then
                if(dmax1*1.0d-7.gt.thres1) then
                   call make_nfock2_symC_uhf(
     $             ilen,        jlen,       klen,      llen,
     $             icff,        jcff,       kcff,      lcff,
     $             dena,        denb,       focka,     fockb,
     $             lind,        nat3,       dmax,      nsym,
     $             ifp,         int72,      xinte72)
                endif
              endif
            endif
c
          enddo                  !   do iqrt=1,nqrt
        enddo                 !   do isbl=1,nblstore
      enddo                !   do isplit
c----------------------------------------------------------------
      if(.not.rescale) return
      if(thres0.ne.thres1) then
        resc=thres0/thres1
        call dscal(nat3*ntri,resc,focka,1)
        call dscal(nat3*ntri,resc,fockb,1)
      endif
c
      end
c====================================================================
      subroutine make_nfock8_nos_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        intx8,      xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock8_sym_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        nsym,       ifp,       intx8,
     $        xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dena(1,jlc),1,focka(1,ikc),1)
                  call daxpy(nat3,-xex,dena(1,ikc),1,focka(1,jlc),1)
                  call daxpy(nat3,-xex,dena(1,ilc),1,focka(1,jkc),1)
                  call daxpy(nat3,-xex,dena(1,jkc),1,focka(1,ilc),1)
                  call daxpy(nat3,-xex,denb(1,jlc),1,fockb(1,ikc),1)
                  call daxpy(nat3,-xex,denb(1,ikc),1,fockb(1,jlc),1)
                  call daxpy(nat3,-xex,denb(1,ilc),1,fockb(1,jkc),1)
                  call daxpy(nat3,-xex,denb(1,jkc),1,fockb(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock8_nosC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      intx8,
     $        xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock8_symC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      nsym,
     $        ifp,         intx8,      xinteg8)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      dimension xinteg8(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx8=intx8+1
              xint=xinteg8(intx8)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_nos_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        intx4,      xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_sym_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        nsym,       ifp,       intx4,
     $        xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dena(1,jlc),1,focka(1,ikc),1)
                  call daxpy(nat3,-xex,dena(1,ikc),1,focka(1,jlc),1)
                  call daxpy(nat3,-xex,dena(1,ilc),1,focka(1,jkc),1)
                  call daxpy(nat3,-xex,dena(1,jkc),1,focka(1,ilc),1)
                  call daxpy(nat3,-xex,denb(1,jlc),1,fockb(1,ikc),1)
                  call daxpy(nat3,-xex,denb(1,ikc),1,fockb(1,jlc),1)
                  call daxpy(nat3,-xex,denb(1,ilc),1,fockb(1,jkc),1)
                  call daxpy(nat3,-xex,denb(1,jkc),1,fockb(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock4_nosC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      intx4,
     $        xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock4_symC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      nsym,
     $        ifp,         intx4,      xinteg4)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*4 xinteg4(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx4=intx4+1
              xint=xinteg4(intx4)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_nos_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        intx2,      xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            ik=lind(max0(icf,kcf))+min0(icf,kcf)
            jk=lind(max0(jcf,kcf))+min0(jcf,kcf)
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(icf,lcf))+min0(icf,lcf)
                jl=lind(max0(jcf,lcf))+min0(jcf,lcf)
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_sym_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       idft,      ax,
     $        dmax,        nsym,       ifp,       intx2,
     $        xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            ik=lind(max0(ic1,kc1))+min0(ic1,kc1)
            jk=lind(max0(jc1,kc1))+min0(jc1,kc1)
            do lc1=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                if(idft.gt.0) xint=xint*ax
                il=lind(max0(ic1,lc1))+min0(ic1,lc1)
                jl=lind(max0(jc1,lc1))+min0(jc1,lc1)
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c ..................Exchange part ....................
                call daxpy(nat3,-xint,dena(1,jl),1,focka(1,ik),1)
                call daxpy(nat3,-xint,dena(1,ik),1,focka(1,jl),1)
                call daxpy(nat3,-xint,dena(1,il),1,focka(1,jk),1)
                call daxpy(nat3,-xint,dena(1,jk),1,focka(1,il),1)
                call daxpy(nat3,-xint,denb(1,jl),1,fockb(1,ik),1)
                call daxpy(nat3,-xint,denb(1,ik),1,fockb(1,jl),1)
                call daxpy(nat3,-xint,denb(1,il),1,fockb(1,jk),1)
                call daxpy(nat3,-xint,denb(1,jk),1,fockb(1,il),1)
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  lic=lind(icf)
                  ljc=lind(jcf)
                  lkc=lind(kcf)
                  llc=lind(lcf)
                  ijc=max0(lic,ljc)+min0(icf,jcf)
                  ikc=max0(lic,lkc)+min0(icf,kcf)
                  ilc=max0(lic,llc)+min0(icf,lcf)
                  jkc=max0(ljc,lkc)+min0(jcf,kcf)
                  jlc=max0(ljc,llc)+min0(jcf,lcf)
                  klc=max0(lkc,llc)+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c ..................Exchange part ....................
                  call daxpy(nat3,-xex,dena(1,jlc),1,focka(1,ikc),1)
                  call daxpy(nat3,-xex,dena(1,ikc),1,focka(1,jlc),1)
                  call daxpy(nat3,-xex,dena(1,ilc),1,focka(1,jkc),1)
                  call daxpy(nat3,-xex,dena(1,jkc),1,focka(1,ilc),1)
                  call daxpy(nat3,-xex,denb(1,jlc),1,fockb(1,ikc),1)
                  call daxpy(nat3,-xex,denb(1,ikc),1,fockb(1,jlc),1)
                  call daxpy(nat3,-xex,denb(1,ilc),1,fockb(1,jkc),1)
                  call daxpy(nat3,-xex,denb(1,jkc),1,fockb(1,ilc),1)
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end
c====================================================================
      subroutine make_nfock2_nosC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      intx2,
     $        xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      data half /0.5d0/
c
      do icf=icff+1,icff+ilen
        do jcf=jcff+1,jcff+jlen
          ij=lind(max0(icf,jcf))+min0(icf,jcf)
          do kcf=kcff+1,kcff+klen
            do lcf=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen
c
      end
c====================================================================
      subroutine make_nfock2_symC_uhf(
     $        ilen,        jlen,       klen,      llen,
     $        icff,        jcff,       kcff,      lcff,
     $        dena,        denb,       focka,     fockb,
     $        lind,        nat3,       dmax,      nsym,
     $        ifp,         intx2,      xinteg2)
      implicit real*8 (a-h,o-z)
      dimension dena(nat3,*),focka(nat3,*)
      dimension denb(nat3,*),fockb(nat3,*)
      integer*2 xinteg2(*)
      dimension lind(*)
      dimension ifp(7,*)
      data half /0.5d0/
c
      do ic1=icff+1,icff+ilen
        do jc1=jcff+1,jcff+jlen
          ij=lind(max0(ic1,jc1))+min0(ic1,jc1)
          do kc1=kcff+1,kcff+klen
            do lc1=lcff+1,lcff+llen
              intx2=intx2+1
              xint=xinteg2(intx2)
              if(dmax*abs(xint).gt.half) then
                xin2=xint+xint
                kl=lind(max0(kc1,lc1))+min0(kc1,lc1)
c ..................Coulomb  part ....................
                do iat=1,nat3
                  focka(iat,ij)=focka(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  fockb(iat,ij)=fockb(iat,ij)+
     $                 (dena(iat,kl)+denb(iat,kl))*xin2
                  focka(iat,kl)=focka(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                  fockb(iat,kl)=fockb(iat,kl)+
     $                 (dena(iat,ij)+denb(iat,ij))*xin2
                enddo
c......................................................
                do ns=1,nsym
                  xco=xin2
                  xex=xint
                  icf=ifp(ns,ic1)
                  jcf=ifp(ns,jc1)
                  kcf=ifp(ns,kc1)
                  lcf=ifp(ns,lc1)
                  if(icf.lt.0) then
                    icf=-icf
                    xco=-xco
                    xex=-xex
                  endif
                  if(jcf.lt.0) then
                    jcf=-jcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(kcf.lt.0) then
                    kcf=-kcf
                    xco=-xco
                    xex=-xex
                  endif
                  if(lcf.lt.0) then
                    lcf=-lcf
                    xco=-xco
                    xex=-xex
                  endif
                  ijc=lind(max0(icf,jcf))+min0(icf,jcf)
                  klc=lind(max0(kcf,lcf))+min0(kcf,lcf)
c ..................Coulomb  part ....................
                  do iat=1,nat3
                    focka(iat,ijc)=focka(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    fockb(iat,ijc)=fockb(iat,ijc)+
     $                   (dena(iat,klc)+denb(iat,klc))*xco
                    focka(iat,klc)=focka(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                    fockb(iat,klc)=fockb(iat,klc)+
     $                   (dena(iat,ijc)+denb(iat,ijc))*xco
                  enddo
c......................................................
                enddo
              endif
            enddo         !   over lll=1,llen
          enddo        !   over kkk=1,klen
        enddo       !   over jjj=1,jlen
      enddo      !   over iii=1,ilen

      end

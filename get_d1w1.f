c=====================================================================
      subroutine get_d1w1(idft,   ax,  rhf,  bl,   inx,
     $                    natoms, NQ,  IUNQ)
c--------------------------------------------------------------------
c First it makes blocking of integrals again just for Hessian CPHF
c
c calculates 1st-order density & weighted density matrices
c for hessian and store them on a disk (units 63 & 64)
c--------------------------------------------------------------------
c
c  ARGUMENTS
c
c  idft    -  dft flag
c  ax      -  scaling parameter for exact exchange
c  rhf     -  logical flag for closed-shell
c  bl      -  general storage
c  inx     -  inx(12,ncs)  integer basis set data
c  natoms  -  total number of atoms
c  NQ      -  number of symmetry-unique atoms
c  IUNQ    -  list of symmetry-unique atoms
c--------------------------------------------------------------------
      
      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
c----------------------------------------------------
      character*11 scftype
      character*4 where
c------------------------------------------------
c     common /intbl/ifpp,inxx(100)
c------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /neglect/ eps,eps1,epsr
c------------------------------------------------
      common /runtype/ scftype,where
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /memmax/ ispblx, maxme1,iforwhat
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c------------------------------------------------
c Hessian (Gradient) options :
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------
      dimension bl(*)
      dimension inx(12,*),IUNQ(NQ)
      dimension xintxx(9)
c----------------------------------------------------------------------
c remember the common block /neglect/ eps,eps1,epsr which has been used
c in the INTEG for SCF. It is needed because the integral threshold may
c be different for zero-order and first-order integrals.
c
      reps=eps
      reps1=eps1
      repsr=epsr
      icacher=icache
c----------------------------------------------------------------------
c prepare the integral system for gradient integrals
c----------------------------------------------------------------------
c get integral threshold for CPHF :
c
      call getrval('thref1',thref1)    ! int.thres for G(D1,g0)
      iforwhat=1                       ! integrals for hessian cphf
c
      thresh=thref1
c
      call post_scf(bl,inx,ncache,nprint,thresh,
     *              iforwhat,1,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
c----------------------------------------------------------------------
      call getival('resta',irestart)
      if(irestart.gt.0) then
c -- reset na to total number of centers
        call getival('ndum',ndum)
        na = na+ndum
        call setival('na  ',na)
        call Para_HessInit(ncache,iforwhat,nblocks)
c -- now set it back to number of real atoms
        na = na-ndum
        call setival('na  ',na)
      endif
c----------------------------------------------------------------------
c print cphf info :
c
      call cphf_infoprt(natoms,nq,iunq)
c
c----------------------------------------------------------------------
c get atomic masses :
c
      call getmem(2*natoms,iatmass)
      call setatmass(bl,natoms,bl(iatmass))
      call setival('atmass',iatmass)
c----------------------------------------------------------------------
c
c Calculate
c
c 1st-order density & weighted density matrices
c
      ntri=ncf*(ncf+1)/2
c
c  MM 08/20/2003  I think it is better to branch the uhf
c                 code here
c
      if(rhf)then
        call cphf_d1w1(idft,  ax,   rhf,    bl,  inx,
     $                 ncf,   ntri, natoms, NQ,  IUNQ)
       else
        call cphf_d1w1_uhf(idft,  ax,   rhf,    bl,  inx,
     $                     ncf,   ntri, natoms, NQ,  IUNQ)
       endif
c
c---------------------------------------------------------------------
c  Return the SCF values of eps,eps1,epsr
c
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c---------------------------------------------------------------------
c
c -- If DFT, and single node, delete all DFT grid files
      call getival('nslv',nslv)
      If(idft.GT.0.AND.nslv.EQ.0) Then
c
c -- following chicanery is needed because CPHF uses symmetry-unique
c    atoms differently, and grid files are labelled according to
c    symmetry-unique atom number rather than actual atom number
        call getint(NQ,itmp)
        do i=0,NQ-1
c       inxx(itmp+i) = i+1
        call int_to_bl( bl(itmp), i+1, i+1 )
        enddo
        Call TidyGRID(-1,NQ,bl(itmp),bl,0)
        call retmem(1)
      EndIf
c
      return
      end
c=====================================================================
      subroutine cphf_d1w1(idft,   ax,     rhf,    bl,     inx,
     $                     ncf,    ntri,   natoms, NQ,     IUNQ)
c--------------------------------------------------------------------
c calculates 1st-order density & weighted density matrices
c for hessian and store them on a disk (units 63 & 64)
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
      call getival('ldensi',lden0)
      call getival('lfock0',lfoc0)
c
      call getival('lvec',lvec)      ! eigenvectors
      call getival('lval',lval)      ! eigenvalues
      call getival('nocc',nocc)      ! no of occupied orbitals
c----------------------------------------------------
c put memory mark
c
      call mmark
c----------------------------------------------------
c find out how many atoms can be treated at once in cphf :
c
      call cphf4nat(NQ,ncf,ntri,ncphf,natonce)
c----------------------------------------------------
c
c allocate local memory :
c
      nmemory=natonce*ntri3
c
      call getmem(nmemory,lgmat1)   ! for G(D1,g0)
      call getmem(ntri3 ,lwork1)   ! for pojected G(D1,g0) in CPHF
c
      call getmem(ntri3 ,ls1   )   ! over1 copy
      call getmem(ntri3 ,lf1   )   ! fock1 copy
c
      call getmem(nmemory,ld1   )   ! working d1 for natonce atoms
c
      call getint(natonce,natend)
c----------------------------------------------------
c
c     DO CPHF ;  get D1 & W1 for symmetry unique atoms :
c
      call do_ncphf(natoms,NQ,natonce,IUNQ,ncphf,
     *              idft,ax,bl,inx,ncf,ntri,nocc,
     *              bl(lden0),bl(lfoc0),bl(lvec),bl(lval),
     *              bl(lgmat1),bl(lwork1), bl(ls1),bl(lf1),
     *              bl(natend), bl(ld1) )
c
c----------------------------------------------------------------------
      call secund(tchfe)
      call elapsec(etchfe)
      tchf=(tchfe-tchfb)/60.0d0
      elaps=(etchfe-etchfb)/60.0d0
      write(iout,430) tchf,elaps
  430 format( 'Master CPU time for D1con +  D1 + W1  = '
     *,f8.2,'  Elapsed = ',f8.2,' min'/)
c----------------------------------------------------------------------
      call retmark
c
      end
c===================================================================
      subroutine cphf_infoprt(natoms,natunq,listunq)
      implicit real*8 (a-h,o-z)
      dimension listunq(natoms)
c----------------------------------------------------
c get CPHF parameters :
c
c idelt shows if delta D10 will be used in CPHF or not
c
      call getival('delta',idelt)       ! use delta dens1 (1) or not (0)
c
c new parameter iscre: shows type of screening 1-particle or 2-particle D
c
      call getival('scree',iscre)
c
c new parameter ireset : reset cphf solution after iter=ireset
c
      call getival('reset',ireset)
c
      call getrval('thres',cphf_thresh) ! cphf-threshold
c
c print the info above into the output file :
c
c----------------------------------------------------
c
      call getival('iout',iout)
      call getival('printh',iprint)
c----------------------------------------------------
      write(iout,100)
  100 format(/72('-'))
      write(iout,*)' '
      write(iout,*)'            The CPHF Solver '
      write(iout,*)' '
c
c print the info above into the output file :
c
      write(iout,110) cphf_thresh
  110 format('threshold for 1st-order density = ',1pe11.4)
c
c......................................................
      if(idelt.eq.0) write(iout,120)
      if(idelt.ne.0) write(iout,121)
  120 format('Using delta density in CPHF screening: OFF ')
  121 format('Using delta density in CPHF screening: ON  ')
c......................................................
      if(iscre.eq.1) write(iout,125)
      if(iscre.ne.1) write(iout,126)
  125 format('Using 1-particle density type screening ')
  126 format('Using 2-particle density type screening ')
c......................................................
c
      write(iout,130) natunq
  130 format('number of symmetry unique atoms = ',i3   )
c
c     if(iprint.ge.1) then
c        write(iout,140) (listunq(iat),iat=1,natunq)
c 140    format('here they are : ',  5(i3,1x)/
c    *         5('                ',  5(i3,1x)/))
c     endif
c----------------------------------------------------
      write(iout,*)' '
c
      call f_lush(iout)
c----------------------------------------------------
      end
c===================================================================
      subroutine do_ncphf(natoms,natunq,natonce,listunq,ncphf,
     *                    idft,ax,bl,inx,ncf,ntri,nocc,
     *                    dens0,fock0,vec,val,gmat1,work1, s1,f1,
     *                    natend, d1 )
c------------------------------------------------------------------
c This is the CPHF driver. It loops over batches of perturbations
c (atoms) ncphf times doing CPHF for natonce atoms at the time.
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
c nocc     - number of occupied orbitals
c dens0    - unperturbed density matrix
c fock0    - unperturbed fock    matrix
c vec      - eigenvectors
c val      - eigenvalues
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
c                                      as needed for cphf acceleration
c OUTPUT :
c
c d1(ntri,3,natonce)  - final Ist-order density
c WRITTEN ON A DISK
c------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /lindvec/ lind,idensp
      dimension bl(*)
      dimension inx(12,*)
      dimension listunq(*)
      dimension dens0(ntri), fock0(ntri)
      dimension vec(*),val(*)
      dimension s1(ntri,3),f1(ntri,3)
c
      dimension gmat1(ntri,3,NATONCE)
      dimension work1(ntri,3)
c
      dimension d1(ntri,3,NATONCE)
c
      dimension natend(natonce)
c-----------------------------------------------------------------------
      call getival('iout',iout)
      call getrval('xlvsh',xlvsh)
      call getival('lsmat',ls0)
c-----------------------------------------------------------------------
      nfile61=61      !   fock1                to be read
      nfile62=62      !   overlap1             to be read
      nfile63=63      !   density1             to be written
      nfile64=64      !   weighted density1    to be written
      nfile65=65      !   constant part of D1  to be written
c-----------------------------------------------------------------------
c  HIGHLY DUBIOUS - set here; apparently used in main CPHF routine
c  via depository (THIS USAGE SHOULD BE HIGHLY DISCOURAGED   JB)
c
      call getival('atmass',iatmass)
c-----------------------------------------------------------------------
c For Better convergence criterion :
c tr(h1+g1)d1  and tr[s1(d1fd+dfd1)] or even tr[s1(d1fd+df1d+dfd1)]
c
c calculate  FD   ; write on a disk (60)
c
      lrec=1
      call calcFD(lrec,ncf,ntri,fock0,dens0,s1) ! s1 as a scratch
c
c-----------------------------------------------------------------------
      ncf2=ncf*ncf
c-----------------------------------------------------------------------
c allocate memory for the VCD atomic axial tensors
      call tstival('vcd',ivcd)
      If(ivcd.NE.0) call getmem(9*natoms,laat)
c
c Solve the CPHF equations for NUMBER of atoms at once :
c
      call getrval('thres',thresh)
c
      icphf=0
      do nat_b=1,natunq,natonce
         nat_e=min(natunq,nat_b+natonce-1)   ! last unique atom
         npos=nat_b-1         ! offset for this cycle in the unique list
         icphf=icphf+1        ! cycle no
ckw
         natodo=nat_e-npos
ckw
         write(iout,1234) icphf,nat_e-npos
 1234    format('  Solving the CPHF equations ',I3,' time for ',I3,
     $          ' atoms')
         call f_lush(iout)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
         do iat=nat_b,nat_e
         nat=listunq(iat)
c -- get mass of atom nat
         call getmass1(nat,iat-npos,natoms,bl(iatmass))
         enddo
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c -- calculate constant part of 1st-order density and compute
c    FDS1 matrices  This can now be done in parallel
c
         call getmem(3*ncf*ncf,ifds)
         call para_d1const(natoms, nat_b,  nat_e,  listunq,ncf,
     $                     ntri,   nocc,  bl(lind),dens0,  s1,
     $                     f1,     vec,    val,    bl,    bl(ifds),
     $                     work1)
         call retmem(1)
cccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c        ................................................................
c        SOLVE THE CPHF ; get D1
c
c        use s1(3*ntri) as a scratch (for 1 atom only)
ckw      call chfsol_nat(natoms,natunq,listunq,natonce,npos,
         call chfsol_nat(natoms,natunq,listunq,natodo ,npos,
     *                   idft,ax,bl,inx,nocc,ncf,ntri,bl(lind),thresh,
     *                   vec,val, s1   ,gmat1,work1,
     *                   natend, d1)
c        on return d1 contains final 1st-order density for natonce atoms
c        ................................................................
c
         call para_done(-1)
c
         call getmem(ncf2*3,lww)
c
c -- calculate the atomic axial tensors
         If(ivcd.NE.0)
     *      call get_aat(natoms,nat_b,nat_e, npos,ncf,nocc,ntri,
     *                   listunq,lww,bl(laat))
c
c   read in precomputed FD  matrix from disk, use s1 as scratch
         call read1mat(60,  1 ,ncf2,s1)    ! lrec=1
c
c -- now compute 1st-order weighted densities
c    this can now be done in parallel
c    use work1 as a scratch for wdens1
c
         call para_wdens(nat_b,  nat_e,  listunq,ncf,    ntri,
     $                  bl(lind),dens0,  s1,     f1,     gmat1,
     $                   d1,     bl(lww), work1)
         call retmem(1)
c
         write(iout,*)' '
      enddo
c
      If(ivcd.NE.0) Then
        call WrAAT(natoms,bl(laat))
        call retmem(1)
      EndIf
c
c..............................................................
ctest
cc      write(6,*) ' PRINTOUT IN ROUTINE <DO-NCPHF>'
cc      do iat=1,natoms
cc         call read1mat(nfile63,iat,ntri*3,d1(1,1,iat) )
cc         write(6,*) ' final 1st-order Density matrix  for at=',iat
cc         call drumh(d1(1,1,iat),ncf, 6  ,'dens1-X ')
cc         call drumh(d1(1,2,iat),ncf, 6  ,'dens1-Y ')
cc         call drumh(d1(1,3,iat),ncf, 6  ,'dens1-Z ')
cc         call read1mat(nfile64,iat,ntri*3,gmat1(1,1,iat) )
cc         write(6,*) ' final 1st-order Wensity matrix  for at=',iat
cc         call drumh(gmat1(1,1,iat),ncf, 6  ,'wens1-X ')
cc         call drumh(gmat1(1,2,iat),ncf, 6  ,'wens1-Y ')
cc         call drumh(gmat1(1,3,iat),ncf, 6  ,'wens1-Z ')
cc      enddo
ctest
c..............................................................
c
c -- tell slaves we are completely done
      call para_done(-2)
c
      write(iout,*)' '
c----------------------------------------------------------------------
      end
c===================================================================
      subroutine chfsol_nat(natoms,natunq,listunq,natonce,npos,
     *                      idft,ax,bl,inx,nocc,ncf,ntri,lind,thrs,
     *                      vec,val, d1con,gmat1,work1,
     *                      natend, d1)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c-----------------------------------------------------------------
c This routine solves the CPHF equations for NATONCE atoms at once
c At the end it calculates the G(D1,g0) with a final D1 as needed
c Ist-order weighted density.
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
c nocc     - number of occupied orbitals
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c thrs     - threshold for cphf
c vec      - eigen vectors
c val      - eigen values
c d1con    - storage for constant part of D1 for ONE atom only :
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eoS1)Cv*[CoCv+ + CvCo+]
c gmat1    - storage for G(D1,g0)
c work1    - storage for variable part of D1 resulting from G(D1,g0)
c natend   - storage for a vector showing if cphf for a given atom is
c            done or not (natend(iat)=0 or 1 : not converged or converged)
c d1,r1(ntri,3,natonce) - to store Dens1 & Resid1
c
c OUTPUT :
c
c d1(ntri,3,natonce)  - final Ist-order density
c-----------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension lind(*)
      dimension vec(*),val(*)
c------------------------------------------------
      dimension d1con(ntri,3)
      dimension gmat1(ntri,3,NATONCE)
      dimension work1(ntri,3)
c
      dimension d1(ntri,3,NATONCE)
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
         call solver_nat(natoms, natonce,npos,   natunq, listunq,
     *                   idft,   ax,     inx,    nocc,   ncf,
     *                   ntri,   lind,   bl,  cphf_thres,vec,
     *                   val,    d1con,  gmat1,  natend, work1,
     *                   bl(iresx),   bl(iresy),   bl(iresz),
     *                   bl(ihesx),   bl(ihesy),   bl(ihesz),
     *                   bl(ideltx),  bl(idelty),  bl(ideltz),
     *                   bl(iheltx),  bl(ihelty),  bl(iheltz),
     *                   d1 )
c
c        1st-order density in d1 on return
c------------------------------------------------
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
ctest   if(iprint.ge.3) then
cc        write(6,*) ' PRINTOUT IN SUBROUTINE <CHFSOL_NAT>'
cc        do iat=1,natonce
cc        write(6,*) ' 1st-order density matrix  for at=',iat
cc        call drumh(d1(1,1,iat),ncf, 6  ,'dens1-X ')
cc        call drumh(d1(1,2,iat),ncf, 6  ,'dens1-Y ')
cc        call drumh(d1(1,3,iat),ncf, 6  ,'dens1-Z ')
cc        enddo
ctest   endif
c------------------------------------------------
      call retmark
c------------------------------------------------
      end
c======================================================================
      subroutine cphf_memcheck

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c
      call getival('iout',iout)
c
      call getmem(0,last0)
      call retmem(1)
c
      ntri=ncf*(ncf+1)/2
c     mem_available=lcore - last0
c     write(iout,*)' memory available for cphf is:',mem_available
c     nmat1=mem_available/ntri
c     nmat3= nmat1/3
c     write(iout,*)' nmat1=',nmat1,' nmat3=',nmat3
c
      natoms=na
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
c
      nmemory=  3*ntri3  + 3*ntri3         ! in solver for 3 dens & res.
     *       + 6*ncf                      ! in dens1 for w1,w2(3*ncf)
     *       + ntri3                      ! G(D1,g0) in cphf
c
      last1=last0+nmemory
      if(last1 .ge. lcore) then
        call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
        iout=igetival('iout')
        write (iout,2001) last1-ioffset,lcore-ioffset
        call nerror(1,'CPHF ','Memory: needed and available ',
     *              last1-ioffset,lcore-ioffset)
      endif
c
 2001 format(1x,' more memory needed in CPHF : needed=',i5,
     *      '  available=',i5)
c
      end
c======================================================================
      subroutine add2vec( ndim, avec , bvec , cvec )
c                        dim     a  +   b   =   c
      implicit real*8 (a-h,o-z)
      dimension avec(ndim),bvec(ndim),cvec(ndim)
c
      do i=1,ndim
         cvec(i)=avec(i)+bvec(i)
      enddo
c
      end
c======================================================================
      subroutine sub2vec( ndim, avec , bvec , cvec )
c                        dim     a  -   b   =   c
      implicit real*8 (a-h,o-z)
      dimension avec(ndim),bvec(ndim),cvec(ndim)
c
      do i=1,ndim
         cvec(i)=avec(i)-bvec(i)
      enddo
c
      end
c======================================================================
      subroutine wdens_xyz(rhf,ncf,ntri,lind,
     *                 work1,work2,
     *                 dens0,fock0,
     *                 dens1,fock1,gmat1,
     *                 wens1)
      implicit real*8 (a-h,o-z)
      logical rhf
      parameter (zero=0.d0, half=0.5d0 )
      dimension lind(*)
      dimension work1(ncf,ncf,3), work2(ncf,ncf,3)
      dimension dens0(*), fock0(*)        ! ntri
      dimension dens1(*), fock1(*), gmat1(*) ! ntri,3
      dimension wens1(*)  ! output 1st-order weighted density
c--------------------------------------------------------------
c Input :
c rhf     - logical flag for rhf/uhf
c ncf     - no of contracted basis functions
c ntri    - ncf*(ncf+1)/2
c work1   - scratch array
c work2   - scratch array
c dens0   - 0th-order density
c fock0   - 0th-order fock
c dens1   - 1st-order density
c fock1   - h1 + G(d0,g1)  " derivative fock "
c gmat1   -      G(d1,g0)
c
c Output :
c wens1   - 1st-order weighted density define as
c
c     W1 = D1*F*D + D*F1*D + D*F*D1
c--------------------------------------------------------------
      ntr2=ntri*2
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            dxf=0.d0
            dyf=0.d0
            dzf=0.d0
c
            fdx=0.d0
            fdy=0.d0
            fdz=0.d0
c
            dfx=0.d0
            dfy=0.d0
            dfz=0.d0
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               kj=jj+k
c
               if(k.gt.i) ik=kk+i
               if(k.gt.j) kj=kk+j
c
               dxf=dxf+dens1(ik)     *fock0(kj)
               dyf=dyf+dens1(ik+ntri)*fock0(kj)
               dzf=dzf+dens1(ik+ntr2)*fock0(kj)
c
               fdx=fdx+dens1(kj)     *fock0(ik)
               fdy=fdy+dens1(kj+ntri)*fock0(ik)
               fdz=fdz+dens1(kj+ntr2)*fock0(ik)
c
               dfx=dfx+dens0(ik)*(fock1(kj)     +gmat1(kj     ))
               dfy=dfy+dens0(ik)*(fock1(kj+ntri)+gmat1(kj+ntri))
               dfz=dfz+dens0(ik)*(fock1(kj+ntr2)+gmat1(kj+ntr2))
            enddo
c
            work1(i,j,1)=dxf + dfx
            work1(i,j,2)=dyf + dfy
            work1(i,j,3)=dzf + dfz          ! D1*F0 + D0*F1
c
            work2(i,j,1)=fdx
            work2(i,j,2)=fdy
            work2(i,j,3)=fdz                !  F0*D1
c
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            ij=ii+j
            if(j.gt.i) ij=jj+i
            xdfd=0.d0
            ydfd=0.d0
            zdfd=0.d0
c
            do k=1,ncf
               kk=lind(k)
               kj=jj+k
               if(k.gt.j) kj=kk+j
               ik=ii+k
               if(k.gt.i) ik=kk+i
c
               xdfd=xdfd + work1(i,k,1)*dens0(kj)   ! (D1*F0 + D0*F1)*D0
               ydfd=ydfd + work1(i,k,2)*dens0(kj)
               zdfd=zdfd + work1(i,k,3)*dens0(kj)
c
               xdfd=xdfd + work2(k,j,1)*dens0(ik)   ! D0*( F0*D1 )
               ydfd=ydfd + work2(k,j,2)*dens0(ik)
               zdfd=zdfd + work2(k,j,3)*dens0(ik)
c
            enddo
            wens1(ij     )=xdfd
            wens1(ij+ntri)=ydfd
            wens1(ij+ntr2)=zdfd
c
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c
c    for rhf we need to divide by two
c
      if(rhf)call VScal(ntri*3,0.5d0,wens1)
c
      end
c======================================================================
      subroutine wdens_1dir(rhf,ncf,ntri,lind,
     *                 work1,work2,
     *                 dens0,fock0,
     *                 dens1,fock1,gmat1,
     *                 wens1)

      use memory

      implicit real*8 (a-h,o-z)
      logical rhf
      parameter (zero=0.d0, half=0.5d0 )
      dimension lind(*)
      dimension work1(ncf,ncf), work2(ncf,ncf)
      dimension dens0(ntri), fock0(ntri)
      dimension dens1(ntri), fock1(ntri), gmat1(ntri)
      dimension wens1(ntri)  ! output 1st-order weighted density
c--------------------------------------------------------------
c Input :
c rhf     - logical flag for rhf/uhf
c ncf     - no of contracted basis functions
c ntri    - ncf*(ncf+1)/2
c work1   - scratch array
c work2   - scratch array
c dens0   - 0th-order density
c fock0   - 0th-order fock
c dens1   - 1st-order density
c fock1   - h1 + G(d0,g1)  " derivative fock "
c gmat1   -      G(d1,g0)
c
c Output :
c wens1   - 1st-order weighted density define as
c
c     W1 = D1*F*D + D*F1*D + D*F*D1
c--------------------------------------------------------------
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            dxf=0.d0
            fdx=0.d0
            dfx=0.d0
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               kj=jj+k
               if(k.gt.i) ik=kk+i
               if(k.gt.j) kj=kk+j
               dxf=dxf+dens1(ik)     *fock0(kj)
               fdx=fdx+dens1(kj)     *fock0(ik)
               dfx=dfx+dens0(ik)*(fock1(kj)     +gmat1(kj     ))
            enddo
            work1(i,j)=dxf + dfx
            work2(i,j)=fdx
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            ij=ii+j
            if(j.gt.i) ij=jj+i
            xdfd=0.d0
            do k=1,ncf
               kk=lind(k)
               kj=jj+k
               if(k.gt.j) kj=kk+j
               ik=ii+k
               if(k.gt.i) ik=kk+i
               xdfd=xdfd + work1(i,k)*dens0(kj)   ! (D1*F0 + D0*F1)*D0
               xdfd=xdfd + work2(k,j)*dens0(ik)   ! D0*( F0*D1 )
            enddo
            wens1(ij     )=xdfd
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c
c    for rhf we need to divide by two
c
      if(rhf)call VScal(ntri,0.5d0,wens1)
c
      end
c======================================================================
      subroutine cphf4nat(natunq,ncf,ntri,ncphf, natonce)

      use memory

c-------------------------------------------------------------
c This routine estimates how many atoms can
c be treated at once in the CPHF procedure
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
c-------------------------------------------------------------
c the following allocations will be needed :
c
c     call getmem(ntri*3,ls1   )   ! fock1 copy
c     call getmem(ntri*3,lf1   )   ! over1 copy
c
c     nmemory=ntri*3*natonce
c
c     call getmem(nmemory,lgmat1)   ! for G(d1,g0)
c     call getmem(ntri*3,lwork1)   ! for projected G(d1,g0) in cphf
c
c     call getmem(nmemory,ld1   )   ! working d1 for natonce atoms
c
c     call memo1_int(natonce,natend)
c
c and also
c
c     call getmem(natonce,iresx)   ! maximum residual R for each atom
c     call getmem(natonce,iresy)
c     call getmem(natonce,iresz)
c     call getmem(natonce,ihesx)   ! maximum residual H for each atom
c     call getmem(natonce,ihesy)
c     call getmem(natonce,ihesz)
c
c   and 6 times the same for deltas
c
c   and IMPORTANT ;
c
c    ncf*ncf*3           for  FDS1 matrices
c
c    9*natonce*natonce    for partial Hess
c-------------------------------------------------------------
c To do CPHF for NATONCE atoms at once we need :
c
c First remember that at the begining of the int_d1g0nat routine
c there is a memory alloc. for 3*natonce*ntri when the
c trsp_inplace(bl,fock,ntri,3*natonce) routine is called. Then
c for CPHF we need :
c
c   needed=3*(natonce*ntri*3) + 3*(ntri*3) + 13*natonce
c         + 3*ncf*ncf
c         + 9*natonce*natonce
c
c  hence
c
c   needed=natonce*(9*ntri+13)+natonce*natonce*9
c          + 9*ntri + 3*ncf*ncf
c
c assuming natonce = min(natunq,50)= natoncx
c
c   needed=natonce*(9*ntri+13)+natonce*natoncX*9
c          + 9*ntri + 3*ncf*ncf
c
c   needed=natonce*(9*ntri + 13 + natoncX*9) + 9*ntri + 3*ncf*ncf
c
c hence : natonce = (needed - 9*ntri-3*ncf*ncf)/(9*ntri+13+natoncX*9)
c-------------------------------------------------------------
c To handle one atom in cphf :
c
      memory1=18*ntri+22+3*ncf*ncf
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
      natonce = (nmemory-9*ntri-3*ncf*ncf)/(9*ntri+13+natoncX)
c
      if(natonce.gt.natunq) natonce=natunq
c
ctest
cc      natonce = 3
ctest
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
c======================================================================
      subroutine solver_nat(natoms, natonce,npos,   natunq, listunq,
     *                      idft,   ax,     inx,    nocc,   ncf,
     *                      ntri,   lind,   bl,     thrs,   vec,
     *                      val,    d1con,  g1,     natend, work1,
     *                      resx,   resy,   resz,
     *                      hesx,   hesy,   hesz,
     *                      deltx,  delty,  deltz,
     *                      heltx,  helty,  heltz,
     *                      r1 )
c-------------------------------------------------------------------
c This routine solves the chfp equations for NATONCE atoms
c INPUT :
c
c natoms   - total number of atoms
c natonce  - number of atoms treated at once in cphf
c npos     - position of current starting atom
c natunq   - number of symmetry-unique atoms
c listunq  - list of symmetry-unique atoms
c idft     - dft flag
c  ax      - dft scaling parameter
c inx      - inx(12,ncs)  basis set info
c nocc     - number of occupied orbitals
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c  bl      - storage for everything
c thrs     - threshold for cphf
c vec      - eigen vectors
c val      - eigen values
c d1con    - storage for constant part of D1 for ONE atom only:
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eoS1)Cv*[CoCv+ + CvCo+]
c g1    - storage for G(D1,g0)
c natend   - storage for a vector showing if cphf for a given atom is
c            done or not (natend(iat)=0 or 1 : not converged or converged)
c resx,resy,resz(iat) - max residals for each atom in cphf
c
c OUTPUT :
c
c r1(ntri,3,natonce)  - final 1st-order density
c-------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extention for files
c
      character*4 screen
      common /screen_type/ screen
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      logical oscylates
c
      dimension bl(*)
      dimension inx(12,*)
      dimension lind(*)
      dimension vec(*),val(*)
c-------------------------------------------------------------------
      dimension natend(natonce),listunq(natunq)
c-------------------------------------------------------------------
      dimension resx(natonce),resy(natonce),resz(natonce)
      dimension hesx(natonce),hesy(natonce),hesz(natonce)
c
      dimension deltx(natonce),delty(natonce),deltz(natonce)
      dimension heltx(natonce),helty(natonce),heltz(natonce)
c-------------------------------------------------------------------
      dimension work1(ntri,3)
      dimension d1con(ntri,3)
      dimension g1(ntri,3,natonce)
      dimension r1(ntri,3,natonce)
c
c
c output r1:1st-order density
c-------------------------------------------------------------------
      call getival('ncs ',ncs)
c-------------------------------------------------------------------
      call getival('atmass',iatmass)
c-------------------------------------------------------------------
      call getival('lsmat' ,ls0)    ! S0 matrix
      call getival('lfock0',lf0)    ! F0 matrix
c-------------------------------------------------------------------
      nfile63=63    !  D1 matrix
      nfile64=64    !  R1 matrix
      nfile65=65    ! constant part of D1
c
      nfile66=66    ! current rO
      nfile67=67    ! current d1
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
      factor=2.0d0     ! for RHF , like orbital's occupancy
c-----------------------------------------------------------------------
      nvirt=ncf-nocc
      call getmem(ncf**2,lw1)      !
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
c use r1,g1,work1 as a scratch
c
      if(idelt.eq.0) then
         call make_screenD(iscre, r1,g1,inx,lind,work1,
     *                     natend,natonce,ncs,ncf,ntri)
      endif
c---------------------------------------------------------------
c The first D1 matrices = D1const
c
      do iat=1,natonce
         call read1mat(nfile65,iat,ntri*3,r1(1,1,iat) ) ! r0=d1const
         call save1mat(nfile66,iat,ntri*3,r1(1,1,iat) ) ! save new r0=r0
         call save1mat(nfile67,iat,ntri*3,r1(1,1,iat) ) ! current d1
      enddo
c---------------------------------------------------------------
cc      write(6,*) ' PRINTOUT IN SUBROUTINE <SOLVER_NAT>'
cc      do iat=1,natonce
cc         write(6,*)' first den1=d1const for atom no=',iat
cc         call drumh(r1(1,1,iat),ncf, 6  ,'D1co -x ')
cc         call drumh(r1(1,2,iat),ncf, 6  ,'D1co -y ')
cc         call drumh(r1(1,3,iat),ncf, 6  ,'D1co -z ')
cc      enddo
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
         call para_cphf_nat(natoms,natend,natonce,npos,natunq,listunq,
     *                      idft,ax,nblocks,bl,inx,ntri,thres_i,lsemi,
     *                      r1,g1,bl(labels))
c
         DO IAT=1,NATONCE
         if(natend(iat).eq.0) then
c
            if(abs(xlvsh).ne.0.d0) then
c              with level shift calculate SD1S  and add it to the G1
               call getmem(ncf*ncf*3,lsds)
               call sd1s_xyz(bl(lsds),ncf,ntri,lind,bl(ls0),r1(1,1,iat),
     *                       xlvsh,g1(1,1,iat))
               call retmem(1)
            endif
c
            do icr=1,3
             call dens1_1dir1n(factor,ncf,nocc,xlvsh,g1(1,icr,iat),
     1                        vec, val, r1(1,icr,iat), bl(lw1),bl(lw2))
            enddo

c
            call file4cphf_o(ntri*3,iat,fext1,r1(1,1,iat),'write') !rl1=L(r1)
c
c   FROM HEREAFTER kepp current D1 in G1
c
            if(noacce.eq.1) then
c              no acceleration
               call calc_d10(nfile65,liter,ntri,iat,G1(1,1,iat),work1)
c                                        output :   current d1 in G1
               call file4cphf_o(ntri*3,iat,fext1,r1(1,1,iat),'read') !rl1=L(r1)
               go to 4000
            endif
c
            call make_r_orto(nfile66,liter,ntri,iat,d1con, r1(1,1,iat))
            call file4cphf_o(ntri*3,iat,fext2,r1(1,1,iat),'write') ! ro1=O(rl1)
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
c             use d1con & g1 as scratch arrays
              call calc_d1r(nfile66,liter,ntri,iat,
     *                     bl(iamat),bl(ibmat),bl(iwvec),
     *                     bl(ilvec),bl(imvec),bl(icoef),
     *                     d1con, work1  ,  G1(1,1,iat))
c                              scratch, out :current d1 in G1
              call retmem(6)
            endif    !  (liter.ge.2) then
c           ......................................
c           calculate current d1
c
            if(liter.eq.1) then
c              call read1mat(nfile67,iat,ntri*3, work1) ! last d1
c              call file4cphf_o(ntri*3,iat,fext1,r1(1,1,iat),'read')
c              call add2vec(ntri*3,work1,r1(1,1,iat),G1(1,1,iat)) ! current d1
               call read1mat(nfile67,iat,ntri*3, G1(1,1,iat)) ! last d1
               call file4cphf_o(ntri*3,iat,fext1,r1(1,1,iat),'read')
               call daxpy(ntri*3,1.0d0,r1(1,1,iat),1,G1(1,1,iat),1) ! current d1
            endif
c           ......................................
c           use orthogonal residuum for the cphf end test :
            if(idft.eq.0) then
               call file4cphf_o(ntri*3,iat,fext2,r1(1,1,iat),'read') !ro
            else
c              call read1mat(nfile67,iat,ntri*3, work1    ) ! last d1
c              call sub2vec(ntri*3,G1(1,1,iat),  work1 ,r1(1,1,iat))
               call read1mat(nfile67,iat,ntri*3, r1(1,1,iat)    ) ! last d1
               call daxpy(ntri*3,-1.0d0,G1(1,1,iat),1,r1(1,1,iat),1)
               call dscal(ntri*3,-1.0d0,r1(1,1,iat),1)
            endif
c.....................................................
c   use non-orthogonal residuum for the cphf end test :
c           if(idft.eq.0) then
c              call file4cphf_o(ntri*3,iat,fext1,r1(1,1,iat),'read') !r1
c           else
c              call read1mat(nfile67,iat,ntri*3,work1  ) ! last d1
c              call sub2vec(ntri*3,G1(1,1,iat), work1 ,r1(1,1,iat))
c           endif
c.....................................................
c   use delta density as residuum for the cphf end test :
c       D1(iter)-D1(iter-1)
c              call read1mat(nfile67,iat,ntri*3,work1  ) ! last d1
c              call sub2vec(ntri*3,G1(1,1,iat), work1 ,r1(1,1,iat))
c.....................................................
 4000       continue   !   if no acceleraion
c.....................................................
c
c           save current D1
c
            call save1mat(nfile67,iat,ntri*3,G1(1,1,iat) ) ! save d1
c
c          ......................................
c
         endif ! if(natend(iat).eq.0) then
         ENDDO ! over atoms
c
c        call secund(citere)
c        call elapsec(eitere)
c        cpuit=citere-citerb
c        elait=eitere-eiterb
c        ............................................................
c--------------------------------------------------------------------
c for convergence test calculate contributions to hessian :
c                FDS1*D1 and S1DF*D1
c
            call getmem(ncf*ncf*3  ,ifds1)
            call getmem(natonce*natonce*9,ihes)
c
c with R1   orthogonal
c
c           use work1 as a scratch for f1
c
            lrec=1
            call HessCont2Diag(lrec,natoms, natend,natonce,ncf,ntri,R1,
     *                  bl(ifds1), work1 ,bl(ihes) )
            call whessContDiag(bl(ihes),bl(iatmass),natoms,natonce)
c
         call secund(citere)
         call elapsec(eitere)
         cpuit=citere-citerb
         elait=eitere-eiterb
c        ............................................................
c
        call cphf_enHR(r1,bl(ihes),thrs,ntri,ncf,liter,cpuit,elait,
     *                 lend, mgo, natend,errmax,natonce,
     *                 resx,resy,resz,    hesx,hesy,hesz,
     *                 deltx,delty,deltz, heltx,helty,heltz,
     *                 oscylates)
            call retmem(2)
c--------------------------------------------------------------------
c
         IF(LEND.EQ.1 .or. ICYCLE.EQ.MXITER) EXIT
c
c....................................................................
         if(noacce.eq.0) then
            do iat=1,natonce
               if(natend(iat).eq.0) then
                  call file4cphf_o(ntri*3,iat,fext2,r1(1,1,iat),'read')
c                 r01 for next iter
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
               if(natend(iat).eq.0) then
                  call file4cphf_o(ntri*3,iat,fext2,r1(1,1,iat),'read')
                  call save1mat(nfile66,iat,ntri*3,G1(1,1,iat) ) ! save d1
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
c read in last D1 soluiton (in R1) :
c
         do iat=1,natonce
            call read1mat(nfile67,iat,ntri*3,R1(1,1,iat) )
         enddo
c---------------------------------------------------------------
c DO ONE MORE CALCULATION OF D1 USING FULL D1 & full integral thresh.
c
         idelt=1                           ! use delta dens in last iter.
         call setival('delta',idelt)
ccccc    screen='2pd1'                     ! much faster
c
         call secund(citerb)
         call elapsec(eiterb)
c---------------------------------------------------------------
         do iat=1,natonce
            natend(iat)=0       ! cphf not converged for this atom yet
         enddo
c---------------------------------------------------------------
cccc     thres10=thres1*10.d0
c
         call para_cphf_nat(natoms,natend,natonce,npos,natunq,listunq,
     *                      idft,ax,nblocks,bl,inx,ntri,thres1,lsemi,
cc   *                      idft,ax,nblocks,bl,inx,ntri,thres10,lsemi,
     *                      R1,g1,bl(labels))
c
c---------------------------------------------------------------
c
c  For VCD reserve memory for perturbed wavefunction coefficients
         call tstival('vcd',ivcd)
         If(ivcd.NE.0) Then
           call getmem(ncf*nocc*3,lpve)
           ncoefile = 79    ! unit for the pert. coeff. file
         EndIf
c
         do iat=1,natonce
            if(abs(xlvsh).ne.0.d0) then
c              with level shift calculate SD1S  and add it to the G1
               call getmem(ncf*ncf*3,lsds)
               call sd1s_xyz(bl(lsds),ncf,ntri,lind,bl(ls0),R1(1,1,iat),
     *                       xlvsh,g1(1,1,iat))
               call retmem(1)
            endif
c
c  For VCD read the constant part of the perturbed coefficient matrix
        if(ivcd.NE.0) call read1mat(ncoefile,iat,ncf*nocc*3,bl(lpve))
c
            do icr=1,3
            call dens1_1dir1n(factor,ncf,nocc,xlvsh,g1(1,icr,iat),
     1                        vec, val, r1(1,icr,iat),bl(lw1),bl(lw2))
            If(ivcd.NE.0)
     *         call daxpy(ncf*nocc,1.0d0,bl(lw2),1,
     *                    bl(lpve+(icr-1)*ncf*nocc),1) !full mx
            enddo
        if(ivcd.NE.0) call save1mat(ncoefile,iat,ncf*nocc*3,bl(lpve))
c.........................................................................
           call read1mat(nfile65,iat,ntri*3, work1)  ! d1const
           call daxpy(ntri*3,1.0d0, work1 ,1,R1(1,1,iat),1 ) ! FINAL D1
           call save1mat(nfile66,iat,ntri*3,R1(1,1,iat) )        ! FINAL D1
         enddo    !  do iat=1,natonce
c
         If(ivcd.NE.0) call retmem(1)
c---------------------------------------------------------------
ctiming
       call secund(citere)
       call elapsec(eitere)
       cpuit=citere-citerb
       elait=eitere-eiterb
       cput=cpuit/60.d0
       elat=elait/60.d0
ctiming
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
          call read1mat(nfile66,iat,ntri*3,r1(1,1,iat)) ! FINAL d1
          call read1mat(nfile67,iat,ntri*3,d1con) ! LAST  d1
          call daxpy(ntri*3, -1.0d0 ,d1con ,1,r1(1,1,iat),1 ) ! FINAL R1
       enddo
c
       call getmem(ncf*ncf*3  ,ifds1)
       call getmem(natonce*natonce*9,ihes)
c
c      use d1con as a scratch for f1
c
       lrec=1
       call HessCont2Diag(lrec ,natoms, natend,natonce,ncf,ntri,R1,
     *               bl(ifds1), d1con ,bl(ihes) )
       call whessContDiag(bl(ihes),bl(iatmass),natoms,natonce)
c
       call final_HR1(r1,bl(ihes),thrs,ntri,ncf,natonce)
c
       call retmem(2)
c---------------------------------------------------------------
 5555  continue
c
       do iat=1,natonce
          call read1mat(nfile66,iat,ntri*3,R1(1,1,iat) ) ! FINAL d1
       enddo
c
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
      subroutine file4cphf_o(ndim,iat,fext,xmat,action)
      implicit real*8(a-h,o-z)
      character*256 jobname,scrf,filename
      common /job/jobname,lenJ
      Character*5 action
      Character*4 fext
      dimension xmat(ndim)
c----------------------------------------------------
c ndim - record length
c----------------------------------------------------
      if(action(1:3).eq.'wri' .or. action(1:3).eq.'rea') then
      else
        call nerror(1,'file4cphf_o','wrong action: '//action,0,0)
      endif
c----------------------------------------------------
      lrec = ndim*8          ! record length in bytes
      nfile =97
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
c
      len = len1 + 6
c
      filename = scrf(1:len1)//'.'//fext
      open (unit=nfile,file=filename(1:len),
     *      form='unformatted',access='direct',recl=lrec)
c
      if(action(1:3).eq.'wri') then
         write(unit=nfile,rec=iat) xmat
      endif
      if(action(1:3).eq.'rea') then
         read(unit=nfile,rec=iat) xmat
      endif
c
      close (nfile,status='keep')
c
      end
c======================================================================
      subroutine get_fext(liter,fext1,fext2)
      character*4 fext1,fext2      ! name extention for files
      character*4 name1(30),name2(30),name3(30)      ! name extention for files
      data name1 /'rl1 ','rl2 ','rl3 ','rl4 ','rl5 ',
     *            'rl6 ','rl7 ','rl8 ','rl9 ','rl10',
     *            'rl11','rl12','rl13','rl14','rl15',
     *            'rl16','rl17','rl18','rl19','rl20',
     *            'rl21','rl22','rl23','rl24','rl25',
     *            'rl26','rl27','rl28','rl29','rl30'/
      data name2 /'ro1 ','ro2 ','ro3 ','ro4 ','ro5 ',
     *            'ro6 ','ro7 ','ro8 ','ro9 ','ro10',
     *            'ro11','ro12','ro13','ro14','ro15',
     *            'ro16','ro17','ro18','ro19','ro20',
     *            'ro21','ro22','ro23','ro24','ro25',
     *            'ro26','ro27','ro28','ro29','ro30'/
c
         fext1=name1(liter)
         fext2=name2(liter)
c
      end
c======================================================================
      subroutine lineqsys(n,b,w,lv,mv,c)
      implicit real*8 (a-h,o-z)
      dimension b(n,n),w(n)
      dimension lv(n),mv(n)
      dimension c(n)            ! coefficients output
      data acc,zero,one /1.d-15 , 0.d0 , 1.d0/
c----------------------------------------------------------
c re-normalize rows :
c
      do 110 i=1,n
      bii1=one
      if(abs(b(i,i)).gt.acc) bii1=one/b(i,i)
      w(i)=w(i)*bii1
      do 110 j=1,n
      b(i,j)=b(i,j)*bii1
  110 continue
c----------------------------------------------------------
c     write(6,*)' b-matrix nxn w(n) : n=',n
c     call f_lush(6)
c
c     if(n.eq.2) then
c        do ii=1,n
c           write(6,62)(b(ii,jj),jj=1,n),w(ii)
c        enddo
c     endif
c 62  format(2(f12.6,1x),2x,f12.6)
c     if(n.eq.3) then
c        do ii=1,n
c           write(6,63)(b(ii,jj),jj=1,n),w(ii)
c        enddo
c     endif
c 63  format(3(f12.6,1x),2x,f12.6)
c----------------------------------------------------------
      call osinv (b,n,det ,acc,lv,mv)
c-------------------------------------------------------------
c     write(6,*)'det(n)=',det,' n=',n
c-------------------------------------------------------------
      if( abs(det).gt.acc)then
         do 111 i=1,n
            sx=zero
            do 222 j=1,n
            sx=sx+b(i,j)*w(j)
  222       continue
            c(i)=sx
  111    continue
      else
         do i=1,n
            c(i)=one
         enddo
      endif
c
c     write(6,*)'    n=',n,' det=',det,' coefficients :'
c     write(6,67) n,c
c 67  format('n=',i2,2x,6(f10.6,1x))
c-------------------------------------------------------------
      end
c======================================================================
      subroutine make_smallA(b,liter,a,liter_n0)
      implicit real*8 (a-h,o-z)
      dimension b(liter,liter)
      dimension a(liter_n0,liter_n0)
c
      do j=1,liter_n0
         do i=1,liter_n0
            a(i,j)=b(i,j)
         enddo
      enddo
c
      end
c======================================================================
      subroutine cphf_int_thr(mgo,liter_thre, errmax,oscylates,
     *                        thres_i,thres_last,thres1,thrs)
      implicit real*8 (a-h,o-z)
      logical oscylates
c
      call getival('iout',iout)
c
      if(thres_i.eq.thres_last) then
         liter_thre=liter_thre+1
      else
         liter_thre=1
         thres_last=thres_i
      endif
c.....................................................
c check if loose integ. threshold sould be sharpen :
c.....................................................
      if(mgo.eq.1) then
         if((errmax.le.0.5d-3.and.errmax.gt.thrs*10.d0)
     *                       .or. oscylates
     *                       .or. liter_thre.gt.4) then
            mgo=2
            if(thres_i.gt.thres1*10.d0) then
               thres_i=thres1*10.d0
               write(iout,153) thres_i
            endif
            go to 1234
  153 format('    Sharp integral threshold = ',1pe11.4)
         endif
         if(errmax.le.thrs*10.d0 ) then
c           cases where ALL three directions almost converged
            mgo=2
            if(thres_i.gt.thres1*10.d0) then
               thres_i=thres1*10.d0
               write(iout,153) thres_i
            endif
            go to 1234
         endif
         go to 1234
      endif
c.....................................................
      if(mgo.eq.2) then
         if(errmax.le.thrs*10.d0 .or. oscylates
     *                           .or. liter_thre.gt.3) then
            mgo=3
            thres_i=thres1
            write(iout,154) thres_i
  154 format('    Final integral threshold = ',1pe11.4)
         endif
         go to 1234
      endif
c.....................................................
 1234 continue
c
      end
c======================================================================
      subroutine file4cphf_c(ndim,liter)
      implicit real*8(a-h,o-z)
      character*256 jobname,scrf,filename
      common /job/jobname,lenJ
      character*4 fext1,fext2      ! name extention for files
c----------------------------------------------------
c delete files open for cphf :
c----------------------------------------------------
      lrec = ndim*8          ! record length in bytes
      nfile =97
c
      call getchval('scrf',scrf)
      call rmblan2(scrf,256,len1)
      len = len1 + 6
c
      do iter=1,liter
         call get_fext(iter,fext1,fext2)
         filename = scrf(1:len1)//'.'//fext1
         open (unit=nfile,file=filename(1:len),
     *         form='unformatted',access='direct',recl=lrec)
         close (nfile,status='delete')
c
         filename = scrf(1:len1)//'.'//fext2
         open (unit=nfile,file=filename(1:len),
     *         form='unformatted',access='direct',recl=lrec)
         close (nfile,status='delete')
      enddo
c
      end
c======================================================================
      subroutine chec_r_orto(nfile,liter,ntri,iat,r,r_curr)
      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extention for files
      dimension  r_curr(ntri,3)    ! current    r
      dimension  r(ntri,3)         ! previous   r
c
      do istep=0,liter-1
         if(istep.eq.0) then
            call read1mat(nfile,iat,ntri*3,r)  ! constant part of D1
         else
            call get_fext(istep,fext1,fext2)
            call file4cphf_o(ntri*3,iat,fext2,r ,'read') !  l(r) previous
         endif
c
         dxldx=ddot(ntri,r(1,1),1,r_curr(1,1),1)
         dyldy=ddot(ntri,r(1,2),1,r_curr(1,2),1)
         dzldz=ddot(ntri,r(1,3),1,r_curr(1,3),1)
         write(6,66) iat,istep,liter,dxldx,dyldy,dzldz
  66     format
     *   ('atom=',i2,' <r-', i2  ,' |r-', i2  ,'> x,y,z=',3(f12.8,2x))
      enddo
c
      end
c======================================================================
      subroutine make_r_orto(nfile,liter,ntri,iat,r,r_curr)
      implicit real*8 (a-h,o-z)
      character*4 fext1,fext2      ! name extention for files
      dimension  r_curr(ntri,3)    ! current    r
      dimension  r(ntri,3)         ! previous   r
      data acc,zero,one /1.d-15 , 0.d0 , 1.d0/
c
c       dn+1 = L(dn) - SUM(i=1,n)[ <di|L(dn)>/<di|di> * di ]
c
c e.g.       d2 = L(d1) - <d1|L(d1)>/<d1|d1> * d1
c
c
      do istep=0,liter-1
         if(istep.eq.0) then
            call read1mat(nfile,iat,ntri*3,r)  ! very first residuum r0
         else
            call get_fext(istep,fext1,fext2)
            call file4cphf_o(ntri*3,iat,fext2,r ,'read') ! previous ro
         endif
c
c        calculate scalars <r|l(r_curr)> & <r|r>
c
         dxldx=ddot(ntri,r(1,1),1,r_curr(1,1),1)
         dyldy=ddot(ntri,r(1,2),1,r_curr(1,2),1)
         dzldz=ddot(ntri,r(1,3),1,r_curr(1,3),1)
c
         dx_dx=ddot(ntri,r(1,1),1,r(1,1),1)
         dy_dy=ddot(ntri,r(1,2),1,r(1,2),1)
         dz_dz=ddot(ntri,r(1,3),1,r(1,3),1)
c
         if(dx_dx.gt.acc) then
            scx=dxldx/dx_dx
         else
            scx=zero
         endif
         if(dy_dy.gt.acc) then
            scy=dyldy/dy_dy
         else
            scy=zero
         endif
         if(dz_dz.gt.acc) then
            scz=dzldz/dz_dz
         else
            scz=zero
         endif
c
c        do i=1,ntri
c           r_curr(i,1)=r_curr(i,1)-scx*r(i,1)
c           r_curr(i,2)=r_curr(i,2)-scy*r(i,2)
c           r_curr(i,3)=r_curr(i,3)-scz*r(i,3)
c        enddo
         call daxpy(ntri,-scx,r(1,1),1,r_curr(1,1),1)
         call daxpy(ntri,-scy,r(1,2),1,r_curr(1,2),1)
         call daxpy(ntri,-scz,r(1,3),1,r_curr(1,3),1)
      enddo
c
      end
c======================================================================
      subroutine calc_d10(nfile,liter,ntri,iat,d,r)
      implicit real*8 (a-h,o-z)
c
c  adds up perturbed densities from previous iteration steps
c
c  ARGUMENTS
c
c  nfile   -  unit number of constant part of D1
c  liter   -  no. of iteration steps
c  ntri    -  dimension of vectors
c  iat     -  record number (seems to be always 1)
c
c -- on exit
c  d       -  total perturbed density
c  r       -  last diff. density (also used as scratch)
c
      character*4 fext1,fext2         ! name extention for files
      dimension  d(ntri,3), r(ntri,3) ! output
c----------------------------------------------------------------------
c read in data :
c
      call read1mat(nfile,iat,ntri*3,d)  ! constant part of D1
c----------------------------------------------------------------------
         do istep=1,liter
            call get_fext(istep,fext1,fext2)
            call file4cphf_o(ntri*3,iat,fext1,r,'read')   ! rl
            call daxpy(ntri*3,1.0d0,r(1,1),1,d(1,1),1)    ! use BLAS
cc            call add2vec(ntri,d(1,icr),r(1,icr),d(1,icr)) ! d=r0+r1+..
         enddo
c----------------------------------------------------------------------
      end
c======================================================================
      subroutine calc_d1r(nfile,liter,ntri,iat,a,b,w,lv,mv,c,ri,rj,d)
      implicit real*8 (a-h,o-z)
      character*4 fexi1,fexi2      ! extention name for files
      character*4 fexj1,fexj2      ! extention name for files
      dimension ri(ntri,3),rj(ntri,3)
      dimension  d(ntri,3)             ! output
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
c d       - resulting 1st-order density matrix (solution at iter=liter)
c there is no r among the parameters!!???      
c r       - "predicted" residuum , needed ONLY for RESET
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
c----------------------------------------------------------------------
c read in data :
c
      call read1mat(nfile,iat,ntri*3,d)  ! get current r0 into d
c----------------------------------------------------------------------
      liter_dir(1)=liter
      liter_dir(2)=liter
      liter_dir(3)=liter
      do istep=1,liter
         call get_fext(istep,fexi1,fexi2)
         call file4cphf_o(ntri*3,iat,fexi1,ri,'read')   !  rl
         do icr=1,3
            call absmax(ntri,ri(1,icr),ix,rmax)
            if(rmax.le.1.0d-7) liter_dir(icr)=liter_dir(icr)-1
         enddo
      enddo
c----------------------------------------------------------------------
      do icr=1,3
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfile,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,iat,fexi2,ri,'read') !  ro
            endif
c
            w(istep)=ddot(ntri,ri(1,icr),1, d(1,icr),1) !<ro|dcon>
            c(istep)=ddot(ntri,ri(1,icr),1,ri(1,icr),1)
c
            do jstep=1,liter_dir(icr)
               call get_fext(jstep,fexj1,fexj2)
               call file4cphf_o(ntri*3,iat,fexj1,rj,'read')        !  rl
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
              call zeroit(d(1,icr),ntri)
c
         do istep=1,liter_dir(icr)
            if(istep.eq.1) then
               call read1mat(nfile,iat,ntri*3,ri )           ! r0
            else
               call get_fext(istep-1,fexi1,fexi2)
               call file4cphf_o(ntri*3,iat,fexi2,ri,'read') !  ro
            endif
            call daxpy(ntri,c(istep), ri(1,icr),1, d(1,icr),1 ) ! final d
         enddo
c
      enddo
c----------------------------------------------------------------------
      end
c======================================================================
      subroutine make_screenD(iscre,d1,r1,inx,lind,map_fs,
     *                        natend,natonce,ncs,ncf,ntri)
      implicit real*8 (a-h,o-z)
      dimension d1(ntri*3,natonce),r1(ntri*3,natonce)
      dimension inx(12,*), lind(*) , map_fs(*)
      dimension natend(natonce)
c
      iunit=65    ! D1const
c
      do iat=1,natonce
         call read1mat(iunit,iat,ntri*3,d1(1,iat)) ! D1const
      enddo
c
      call trspmo(d1,ntri, r1,3*natonce )
      call dcopy(3*ntri*natonce, r1 , 1 , d1, 1)
c
      call select_dnat(natend,natonce,d1,ntri,r1)
c
      call setup_densp2(inx,ncf,ncs,r1,d1,map_fs)
c       screening density NOT squared
c
c D1 now contains 1 density matrix over contracted shells
c
c     write screening density D1 on a disk ! lrec=ncs*ncs NOT  3*ntri
c
c     write it on disk as a second record :
c
      call save1mat(69, 2 ,ncs*ncs  ,d1)
c
      end
c======================================================================
      subroutine calcFD(lrec,ncf,ntri,f0,d0,fd)
      implicit real*8 (a-h,o-z)
      dimension f0(ntri),d0(ntri),fd(ncf,ncf)
c
      do i=1,ncf
         ii=i*(i-1)/2
         do j=1,ncf
            jj=j*(j-1)/2
            sum1=0.d0
            do k=1,ncf
               kk=k*(k-1)/2
               if(i.ge.k) then
                  ik=ii+k
               else
                  ik=kk+i
               endif
               if(k.ge.j) then
                  kj=kk+j
               else
                  kj=jj+k
               endif
               sum1=sum1+f0(ik)*d0(kj)
            enddo
            fd(i,j)=sum1
         enddo ! j=1,ncf
      enddo    ! i=1,ncf
c
c   write fd and df on a disk ; unit=60
c
      ncf2=ncf*ncf
      call save1mat(60,lrec,ncf2,fd)
c
c
c       write(6,*)' fd  matrix'
c       do i=1,ncf
c       do j=1,ncf
c          write(6,66) i,j,fd(i,j)
c       enddo
c       enddo
c
  66    format(2(i2,1x),f10.5,1x)
c
      end
c======================================================================
      subroutine calcFDS1(rhf,    ncf,    ntri,   fd,     s1,
     $                    fds1)
      implicit real*8 (a-h,o-z)
      logical rhf
      dimension s1(ntri,3)
      dimension fd(ncf,ncf)
      dimension fds1(ncf,ncf,3)
c------------------------------------------------------------------------
c Input :
c rhf     - logical flag; true for closed-shell
c ncf     - number of basis functions
c ntri    - ncf*(ncf+1)/2
c fd      - Fock*Density matirx (precalculated and store on 60)
c s1      - 1st-order overlap
c
c Output:
c
c FDS1    - Fock*Density*Overlap1 !written on disk , unit 60)
c------------------------------------------------------------------------
c
c in order to save memory calculate first FDS1
c
      do i=1,ncf
         ii=i*(i-1)/2
         do j=1,ncf
            jj=j*(j-1)/2
            sum1=0.d0
            sum2=0.d0
            sum3=0.d0
            do k=1,ncf
               kk=k*(k-1)/2
               if(i.ge.k) then
                  ik=ii+k
               else
                  ik=kk+i
               endif
               if(k.ge.j) then
                  kj=kk+j
               else
                  kj=jj+k
               endif
               sum1=sum1+fd(i,k)*s1(kj,1)
               sum2=sum2+fd(i,k)*s1(kj,2)
               sum3=sum3+fd(i,k)*s1(kj,3)
            enddo !  k=1,ncf
c
            fds1(i,j,1)=sum1
            fds1(i,j,2)=sum2
            fds1(i,j,3)=sum3
c
         enddo !  j=1,ncf
      enddo !  i=1,ncf
c
c  and then S1DF and add it to FDS1
c
      do i=1,ncf
         ii=i*(i-1)/2
         do j=1,ncf
            jj=j*(j-1)/2
            sum1=0.d0
            sum2=0.d0
            sum3=0.d0
            do k=1,ncf
               kk=k*(k-1)/2
               if(i.ge.k) then
                  ik=ii+k
               else
                  ik=kk+i
               endif
               if(k.ge.j) then
                  kj=kk+j
               else
                  kj=jj+k
               endif
               sum1=sum1+fd(i,k)*s1(kj,1)
               sum2=sum2+fd(i,k)*s1(kj,2)
               sum3=sum3+fd(i,k)*s1(kj,3)
            enddo !  k=1,ncf
c
            fds1(j,i,1)=sum1+fds1(j,i,1)
            fds1(j,i,2)=sum2+fds1(j,i,2)
            fds1(j,i,3)=sum3+fds1(j,i,3)
c
         enddo !  j=1,ncf
      enddo !  i=1,ncf
c
c scale it by 2 for UHF
c
      if(.not.rhf) call dscal(3*ncf**2, 2.d0 ,fds1,1)
c
      return
      end
c======================================================================
      subroutine spuXold(xmat1,xmat2,ncf,ntri,trace)
      implicit real*8 (a-h,o-z)
      dimension xmat1(ntri) , xmat2(ncf,ncf)
c
c  calculates trace of X1(ntri)*X2(ncf,ncf)
c
            trace=0.d0
            do i=1,ncf
               ii=i*(i-1)/2
               do k=1,ncf
                  kk=k*(k-1)/2
                  if(k.ge.i) then
                     ki=kk+i
                  else
                     ki=ii+k
                  endif
                  trace=trace+xmat2(i,k)*xmat1(ki)
               enddo
            enddo
c
      end
c======================================================================
      subroutine spuX(xmat1,xmat2,ncf,ntri,trace)
      implicit real*8 (a-h,o-z)
      parameter (zero=0.0d0)
      dimension xmat1(ntri) , xmat2(ncf,ncf)
c
c  calculates trace of X1(ntri)*X2(ncf,ncf)
c  Arguments
c  INTENT(IN)
c  Xmat1    = symmetric matrix stored in the columnwise upper triangle
c  form, Xmat(ij)=X1(i,j) where (ij)=i*(i-1)/2+j, i>=j, X is symmetrical
c  Xmat2    = square matrix
c  ncf      = dimension of the (square) matrices X1 (Xmat1) and Xmat2
c  ntri     = should be ncf*(ncf+1)/2
c  INTENT(OUT)
c  trace    = Trace(X1*Xmat2)

c
      trace=zero
      ii=0
      do i=1,ncf
        do k=1,i
          ik=ii+k
          trace=trace+(xmat2(i,k)+xmat2(k,i))*xmat1(ik)
        enddo
        ii=ii+i
      enddo
c  Subtract the overcounted diagonal
      ii=0
      do i=1,ncf
        ii=ii+i
        trace=trace-xmat2(i,i)*xmat1(ii)
      end do
      end
c======================================================================
      subroutine cphf_enHR(r1,dhess,thrs,ntri,ncf,liter,cpuit,elait,
     *                     lend, mgo, natend,errmax,natonce,
     *                     resx,resy,resz,    hesx,hesy,hesz,
     *                     deltx,delty,deltz, heltx,helty,heltz,
     *                     oscylates)
c----------------------------------------------------------------------
c This routine checks the convergence in the CPHF procedure
c
c INPUT :
c
c natonce            -  number of atoms treated at once in CPHF
c dhess              -  changes in the hessain (delta hessian)
c thrs               -  cphf threshold
c ntri               -  ncf*(ncf+1)/2
c
c INPUT/OUTPUT :
c                       integral threshold is changed upon mgo
c
c INPUT :
c liter              -  number of current CPHF iteration
c cpuit,elait        -  cpu & elapsed time of iteration
c
c OUTPUT
c lend               - shows end of CPHF : 1=converged, 0=not conv.
c
c INPUT/OUTPUT :
c natend             - natend(iatom)=1 or 0 : cophf converged or not for Iatom
c
c OUTPUT
c errmax             - maximum element in delta density
c INPUT/OUTPUT :
c resx,resy,resz     - maximu residuals for each atom
c oscylates          - logical showing smootness of CPHF
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical oscylates
      dimension r1(ntri,3,natonce)     ! last residiuum
      dimension dhess(3,natonce,3,natonce)
c
      dimension resx(natonce),resy(natonce),resz(natonce)
      dimension hesx(natonce),hesy(natonce),hesz(natonce)
c
      dimension heltx(natonce),helty(natonce),heltz(natonce)
      dimension deltx(natonce),delty(natonce),deltz(natonce)
c
      dimension natprt(10000)    ! for local print only
      dimension natend(*)
c
      data zero /0.d0/
      data natconv /0/
      save natconv
c
c..............................................
c natend(iat) = 0 (not done) or 1 (done)
c..............................................
       thrs10= thrs*10.d0
c..............................................
c
      call getival('iout',iout)
      call getival('printh',iprint)
c..............................................
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
           hx_max=zero
           hy_max=zero
           hz_max=zero
ctry       do jat=1,natonce
           do jcr=1,3
              hx= abs( dhess(1,iat, jcr,Iat) )
              hy= abs( dhess(2,iat, jcr,Iat) )
              hz= abs( dhess(3,iat, jcr,Iat) )
              hx_max=max(hx_max,hx)
              hy_max=max(hy_max,hy)
              hz_max=max(hz_max,hz)
           enddo
ctry       enddo
           heltx(iat)=hx_max
           helty(iat)=hy_max
           heltz(iat)=hz_max
        endif
      ENDDO
c..............................................
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
           deltx(iat)=zero
           delty(iat)=zero
           deltz(iat)=zero
           call absmax(ntri,r1(1,1,iat),ix,deltx(iat))
           call absmax(ntri,r1(1,2,iat),iy,delty(iat))
           call absmax(ntri,r1(1,3,iat),iz,deltz(iat))
        endif
      ENDDO
c..............................................
c
      oscylates=.false.
      noscx=0
      noscy=0
      noscz=0
c
      errmax=0.d0
      pelRx=0.d0
      pelRy=0.d0
      pelRz=0.d0
      pelHx=0.d0
      pelHy=0.d0
      pelHz=0.d0
      do iat=1,natonce
        if(natend(iat).eq.0) then
c..............
         if(liter.gt.3) then
           if(resx(iat).gt.thrs10 .and. hesx(iat).gt.thrs10) then
            if(deltx(iat).gt.resx(iat).and.heltx(iat).gt.hesx(iat)) then
               noscx=noscx+1
               oscylates=.true.
            endif
           endif
           if(resy(iat).gt.thrs10 .and. hesy(iat).gt.thrs10) then
            if(delty(iat).gt.resy(iat).and.helty(iat).gt.hesy(iat)) then
               noscy=noscy+1
               oscylates=.true.
            endif
           endif
           if(resz(iat).gt.thrs10 .and. hesz(iat).gt.thrs10) then
            if(deltz(iat).gt.resz(iat).and.heltz(iat).gt.hesz(iat)) then
               noscz=noscz+1
               oscylates=.true.
            endif
           endif
         endif
c..............
           erRiat=max(deltx(iat),delty(iat),deltz(iat))
           erHiat=max(heltx(iat),helty(iat),heltz(iat))
           errorIat=MIN(erRiat,erHiat)
CCCCCCCC   if(errorIat.le.thrs) natend(iat)=1
           if(errorIat.le.thrs .and. mgo.eq.3) natend(iat)=1
           errmax=max(errorIat,errmax)
c
           pelRx=max(     deltx(iat),pelRx)
           pelRy=max(     delty(iat),pelRy)
           pelRz=max(     deltz(iat),pelRz)
c
           pelHx=max(     Heltx(iat),pelHx)
           pelHy=max(     Helty(iat),pelHy)
           pelHz=max(     Heltz(iat),pelHz)
c..............
           resx(iat)=deltx(iat)
           resy(iat)=delty(iat)
           resz(iat)=deltz(iat)
c
           hesx(iat)=heltx(iat)
           hesy(iat)=helty(iat)
           hesz(iat)=heltz(iat)
c..............
        endif      !   (natend(iat).eq.0) then
      enddo
c
      lend=1
      if(errmax.gt.thrs) lend=0
c
      cput=cpuit/60.d0
      elat=elait/60.d0
c
      iatconv=0
      do iat=1,natonce
        if(natend(iat).eq.1) then
           iatconv=iatconv+1
           natprt(iatconv)=iat
        endif
      enddo
c
      if(oscylates) then
         write(iout,421) liter,pelRx,pelRy,pelRz,cput,elat,oscylates,
     *                   noscx,noscy,noscz
         write(iout,4200) pelHx,pelHy,pelHz
      else
         write(iout,420) liter,pelRx,pelRy,pelRz,cput,elat,oscylates
         write(iout,4200) pelHx,pelHy,pelHz
c
c        write(iout,520) liter,pelRx,pelRy,pelRz,
c    *                         pelHx,pelHy,pelHz,
c    *                   cput,elat,oscylates
c 520 format(i3,3x,3(1x,e10.3),1x,3(1x,e10.3),2(f8.2,1x),l5)
      endif
c
      if(iatconv.gt.natconv) then
         if(iprint.ge.1) then
            write(iout,422) iatconv,(natprt(k),k=1,iatconv)
         else
            write(iout,423) iatconv
         endif
      endif
      natconv=iatconv
c
      call f_lush(iout)
c
c
 4200 format(  6x ,3(1x,e11.4) )
  420 format(i3,3x,3(1x,e11.4),1x,2(f8.2,1x),l5)
  421 format(i3,3x,3(1x,e11.4),1x,2(f8.2,1x),l5,
     *        ' xyz=',3(i2,1x))
  422 format(3x,' cphf converged for ',i4,' unique atoms no :',6(i3,1x)/
     *     ,6(3x,'                    ',2x,'                  ',
     *                                                       6(i3,1x)/))
  423 format(3x,' cphf converged for ',i4,' unique atoms ')
c
      end
c======================================================================
      subroutine final_HR1(r1,dhess,thrs,ntri,ncf,natonce)
c----------------------------------------------------------------------
c This routine checks the convergence in the CPHF procedure
c
c INPUT :
c
c natonce            -  number of atoms treated at once in CPHF
c dhess              -  changes in the hessain (delta hessian)
c thrs               -  cphf threshold
c ntri               -  ncf*(ncf+1)/2
c
c INPUT/OUTPUT :
c                       integral threshold is changed upon mgo
c
c OUTPUT
c INPUT/OUTPUT :
c resx,resy,resz     - maximu residuals for each atom
c oscylates          - logical showing smootness of CPHF
c----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*1 icoor
      character*1 icooh,jcooh
      dimension r1(ntri,3,natonce)     ! last residium
      dimension dhess(3,natonce,3,natonce)
c
      data zero /0.d0/
c..............................................
      call getival('iout',iout)
      call getival('printh',iprint)
c..............................................
      nat29=natonce*natonce*9
      call absmax(nat29,dhess,ij_max,hij_max)
c..............................................
      ij=0
      do iat=1,natonce
        do jat=1,natonce
           do icr=1,3
              do jcr=1,3
                 ij=ij+1
                 if(ij.eq.ij_max) then
                    iatH=iat
                    jatH=jat
                    i_h=icr
                    j_h=jcr
                    go to 123
                 endif
              enddo
           enddo
        enddo
      enddo
c
  123 continue
c
c maximum residuum in HESSian is hij_max
c..............................................
      r_max=zero
      DO IAT=1,NATONCE
           call absmax(ntri,r1(1,1,iat),ix,rx)
           call absmax(ntri,r1(1,2,iat),iy,ry)
           call absmax(ntri,r1(1,3,iat),iz,rz)
           if(rx.gt.r_max) then
              r_max=rx
              iatR=iat
              i_r=1
           endif
           if(ry.gt.r_max) then
              r_max=ry
              iatR=iat
              i_r=2
           endif
           if(rz.gt.r_max) then
              r_max=rz
              iatR=iat
              i_r=3
           endif
      ENDDO
c
c..............................................
      if(i_r.eq.1) icoor='x'
      if(i_r.eq.2) icoor='y'
      if(i_r.eq.3) icoor='z'
      if(i_h.eq.1) icooh='x'
      if(i_h.eq.2) icooh='y'
      if(i_h.eq.3) icooh='z'
      if(j_h.eq.1) jcooh='x'
      if(j_h.eq.2) jcooh='y'
      if(j_h.eq.3) jcooh='z'
c..............................................
      write(iout,100)
 100  format(/' Maximum residuals in D1 matrix and Hessian :')
c
      write(iout,101) icoor, r_max, iatR
      write(iout,201) icooh,jcooh, hij_max,iatH,jatH
 101  format('    in D1',a1,2x,1pe11.2,'  on atom   ',i2)
 201  format('    in H', a1,a1,2x,1pe11.2,'  on atoms  ',2(i2,1x))
c
c
      call f_lush(iout)
c
      end
c======================================================================
      subroutine setatmass(bl,natoms,atmass)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      character*256 jobname,scrdir,filename
      common /job/jobname,lenJ
      dimension bl(*)
      dimension atmass(natoms)
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(NAtoms)
c
      iunit=1
c
      call getmem(natoms*3,ixc)
      call getmem(natoms  ,ixq)
c
      jnk=1
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call RdCoordF(IUnit,NAtoms,AtSymb,bl(iXC),-1,jnk,bl(ixq),AtMASS)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
      call retmem(2)
c
      end
c======================================================================
      subroutine getmass1(nat,iat,natoms,atmass)
      implicit real*8 (a-h,o-z)
      dimension atmass(natoms,2)
c
ccc   atmass(iat,2)=atmass(nat,1)
      atmass(iat,2)=1.0d0/sqrt(atmass(nat,1)) ! for weighted hess.
c
      end
c======================================================================
      subroutine whessContDiag(hess,atmass,natoms,natonce)
      implicit real*8 (a-h,o-z)
      dimension hess(3,natonce,3,natonce)
      dimension atmass(natoms,2)
c
      do iat=1,natonce
         xmasi=atmass(iat,2)
         fact=xmasi*xmasi
         do icr=1,3
            do jcr=1,3
               hess(icr,iat,jcr,iat)=hess(icr,iat,jcr,iat)*fact
            enddo
         enddo
      enddo
c
      end
c======================================================================
      subroutine wdens1(rhf,    ncf,    ntri,   lind,   d0,
     $                  fd,     f1,     g1,     d1,     work,
     $                  w1)
      implicit real*8 (a-h,o-z)
c--------------------------------------------------------------
c This routine calculates the  Ist-order weighted density(x,y,z)
c one atom, 3 directions at once
c--------------------------------------------------------------
c At this point we have 1st-order density D1 and G(D1,g0).
c We need to calculate 1-st order "weighted density" W1 which
c contributes to the hessian as -2*Tr W1*S1 . This contribution
c to the final hessian can be expressed in terms of ordinary D1
c and 0th- and 1st-order FULL Fock as follows :
c
c -2 TrSa Wb =  -Tr Sa (Db F D + D Fb  D + D F Db )
c
c where F=h+G(D,g)  and  Fb=hb + G(Db,g)+G(D,gb)
c
c Comparing Tr Sa*Wb with Tr Sa*(Db*F*D+D*Fb*D+D*F*Db)
c one may named the (Db F D + D Fb D + D F Db) matrix
c the 1st-order "weighted density" . Thus we define Wb as :
c
c   Wb = Db*F*D + D*Fb*D + D*F*Db
c
c with full fock :  F=h+G(D,g)  and  Fb=hb + G(Db,g)+G(D,gb)
c
c  G(D1,g0) already done in chfsol_nat
c--------------------------------------------------------------
c INPUT :
c rhf      - logical flag for rhf/uhf
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c lind     - diagonals of i*(i-1)/2
c d0       - unpreturbed density
c fd       - FD matrix
c f1       - part of Ist-order fock matrix f1=h1 + G(D0,g1)
c g1       - G(D1,g0)
c d1       - Ist-order Density
c work     - scratch for D*F1
c
c OUTPUT :
c w1       - resulting Ist-order weighted density
c--------------------------------------------------------------
      logical rhf
      dimension lind(*)
      dimension d0(*)         !   (ntri)
      dimension fd(ncf,ncf)
      dimension work(ncf,ncf,3)
c--------------------------------------------------------------
      dimension f1(ntri,3)
      dimension g1(ntri,3)
c
      dimension d1(ntri,3)
      dimension w1(ntri,3)
c--------------------------------------------------------------
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            dxfd=0.d0
            dyfd=0.d0
            dzfd=0.d0
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               kj=jj+k
               if(k.gt.i) ik=kk+i
               if(k.gt.j) kj=kk+j
               dxfd=dxfd+d1(ik,1)*fd(k,j)
               dyfd=dyfd+d1(ik,2)*fd(k,j)
               dzfd=dzfd+d1(ik,3)*fd(k,j)       !   D1*FD
            enddo
            work(i,j,1)=dxfd
            work(i,j,2)=dyfd
            work(i,j,3)=dzfd               !  D0*F1
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            ij=ii+j
            if(j.gt.i) ij=lind(j)+i
            w1(ij,1)=work(i,j,1)+work(j,i,1)
            w1(ij,2)=work(i,j,2)+work(j,i,2)
            w1(ij,3)=work(i,j,3)+work(j,i,3)
         enddo
      enddo
c--------------------------------------------------------------
      do i=1,ncf
         ii=lind(i)
         do j=1,ncf
            jj=lind(j)
            dfx=0.d0
            dfy=0.d0
            dfz=0.d0
            do k=1,ncf
               kk=lind(k)
               ik=ii+k
               kj=jj+k
               if(k.gt.i) ik=kk+i
               if(k.gt.j) kj=kk+j
               dfx=dfx+d0(ik)*( f1(kj,1)+g1(kj,1) )
               dfy=dfy+d0(ik)*( f1(kj,2)+g1(kj,2) )
               dfz=dfz+d0(ik)*( f1(kj,3)+g1(kj,3) )
            enddo
            work(i,j,1)=dfx
            work(i,j,2)=dfy
            work(i,j,3)=dfz                !  D0*F1
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c--------------------------------------------------------------
      ij=0
      do i=1,ncf
         ii=lind(i)
         do j=1,i
            jj=lind(j)
            ij=ij+1
            dfxd=0.d0
            dfyd=0.d0
            dfzd=0.d0
            do k=1,ncf
               kk=lind(k)
               kj=jj+k
               if(k.gt.j) kj=kk+j
cccccc         ik=ii+k
cccccc         if(k.gt.i) ik=kk+i
               dfxd=dfxd + work(i,k,1)*d0(kj)      ! (D0*F1)*D0
               dfyd=dfyd + work(i,k,2)*d0(kj)
               dfzd=dfzd + work(i,k,3)*d0(kj)
            enddo
            w1(ij,1)=w1(ij,1)+dfxd
            w1(ij,2)=w1(ij,2)+dfyd
            w1(ij,3)=w1(ij,3)+dfzd
         enddo !   j=1,ncf
      enddo    !   i=1,ncf
c--------------------------------------------------------------
c
c     for rhf we need to divide by two
c
      if(rhf)call VScal(ntri*3,0.5d0,w1)
c
c--------------------------------------------------------------
c       Ist-order weighted density in W1 on return
c--------------------------------------------------------------
c     call getival('printh',iprint)
c     if(iprint.ge.3) then
c        write(6,*) ' 1st-order weighted density '
c        call drumh(w1(1,1),ncf, 6  ,'wens1-X ')
c        call drumh(w1(1,2),ncf, 6  ,'wens1-Y ')
c        call drumh(w1(1,3),ncf, 6  ,'wens1-Z ')
c     endif
c--------------------------------------------------------------
      end
c======================================================================
      subroutine HessCont2Diag(lrec ,natoms, natend,natonce,ncf,ntri,d1,
     *                  fds1,f1,hess)
      implicit real*8 (a-h,o-z)
      dimension natend(natonce)
      dimension d1(ntri,3,natonce)
c
c     for 1 atom
      dimension fds1(ncf,ncf,3)
      dimension f1(ntri,3)
      dimension hess(3,natonce,3,natonce)
c---------------------------------------------------------------------
c calculates contributions to the hessian that involve D1
c for atom-diagonal part of the hessian
c---------------------------------------------------------------------
c     fds1 read from file 60 contains : f1 -(fds1+s1df)
c---------------------------------------------------------------------
c     contributions to Hess :
c
c     -0.5*Tr S1*[ D1FD + DF1D  +DFD1]=
c    =-0.5*Tr [D1*FDS1 + S1DF*D1]
c     -0.5*Tr DS1D*F1                 ! this one is not included
c
c                                      (I do not have G(D1,g0) only G(R1,g0)
c     and
c
c    +0.5*Tr [ h1 + G(D,g1) ]*D1     ;  f1=h1+G(d0,g)
c
c---------------------------------------------------------------------
c
      call zeroit(hess,9*natonce*natonce)
c
      ncf2=ncf*ncf
c---------------------------------------------------------------------
c        do iat=1,natonce
c           write(6,*)'iat=',iat,'  natend=',natend(iat)
c        enddo
c---------------------------------------------------------------------
c
c  Iat=Jat :
c
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
c
         irec=(iat-1)*3 + lrec
         call read1mat(60,irec+1,ncf2,fds1(1,1,1) )
         call read1mat(60,irec+2,ncf2,fds1(1,1,2) )
         call read1mat(60,irec+3,ncf2,fds1(1,1,3) )
c
         do ixyz=1,3
            call spuX(d1(1,ixyz,iat),fds1(1,1,ixyz),ncf,ntri,dewe)
            hess(ixyz,iat,ixyz,iat)=hess(ixyz,iat,ixyz,iat)+dewe
            do jxyz=ixyz+1,3
                call spuX(d1(1,ixyz,iat),fds1(1,1,jxyz),ncf,ntri,ws1)
                call spuX(d1(1,jxyz,iat),fds1(1,1,ixyz),ncf,ntri,sw1)
                dewe= ws1+sw1
                dewe= dewe*0.5d0
                hess(ixyz,iat,jxyz,iat)=hess(ixyz,iat,jxyz,iat)+dewe
            enddo
         enddo
        endif    !    (natend(iat).eq.0) then
      ENDDO       !   DO IAT=1,NATONCE
c---------------------------------------------------------------------
c atom-diagonal only
c---------------------------------------------------------------------
c make full hessian out of its upper triangle part :
c
      call hess_full_up(hess,natonce)
c
c---------------------------------------------------------------------
c
c       do iat=1,natonce
c          do jat=1,natonce
c             write(6,*) ' Atoms ',iat,jat
c        write(6,66) ((hess(i,iat,j,jat),i=1,3),j=1,3)
c          enddo
c       enddo
c
c
c 66  format(1x,3(e12.6,1x))
c 66  format(1x,9(f9.5,1x))
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine whessCont(hess,atmass,natoms,natonce)
      implicit real*8 (a-h,o-z)
      dimension hess(3,natonce,3,natonce)
      dimension atmass(natoms,2)
c
      do iat=1,natonce
         xmasi=atmass(iat,2)
ccccc    xmasi=1.d0/sqrt(xmasi)
         do jat=1,natonce
            xmasj=atmass(jat,2)
ccccc       xmasj=1.d0/sqrt(xmasj)
            fact=xmasi*xmasj
            do icr=1,3
               do jcr=1,3
                  hess(icr,iat,jcr,jat)=hess(icr,iat,jcr,jat)*fact
               enddo
            enddo
         enddo
      enddo
c
      end
c======================================================================
      subroutine HessCont2(lrec ,natoms, natend,natonce,ncf,ntri,d1,
     *                  fds1,f1,hess)
      implicit real*8 (a-h,o-z)
      dimension natend(natonce)
      dimension d1(ntri,3,natonce)
c
c     for 1 atom
      dimension fds1(ncf,ncf,3)
      dimension f1(ntri,3)
      dimension hess(3,natonce,3,natonce)
c---------------------------------------------------------------------
c calculates contributions to the hessian that involve D1
c---------------------------------------------------------------------
c     fds1 contains fds1+s1df
c---------------------------------------------------------------------
c     contributions to Hess :
c
c     -0.5*Tr S1*[ D1FD + DF1D  +DFD1]=
c    =-0.5*Tr [D1*FDS1 + S1DF*D1]
c     -0.5*Tr DS1D*F1                 ! this one is not included
c
c                                      (I do not have G(D1,g0) only G(R1,g0)
c     and
c
c    +0.5*Tr [ h1 + G(D,g1) ]*D1     ;  f1=h1+G(d0,g)
c
c---------------------------------------------------------------------
c
      call zeroit(hess,9*natonce*natonce)
c
      ncf2=ncf*ncf
c---------------------------------------------------------------------
c        do iat=1,natonce
c           write(6,*)'iat=',iat,'  natend=',natend(iat)
c        enddo
c---------------------------------------------------------------------
c
c  Iat=Jat :
c
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
c
         irec=(iat-1)*3 + lrec
         call read1mat(60,irec+1,ncf2,fds1(1,1,1) )
         call read1mat(60,irec+2,ncf2,fds1(1,1,2) )
         call read1mat(60,irec+3,ncf2,fds1(1,1,3) )
c
         call read1mat(61,natoms+iat,ntri*3,f1(1,1) )
c
         do ixyz=1,3
            call spur(d1(1,ixyz,iat),  f1(1,ixyz)  ,ncf,d1f1)
            call spuX(d1(1,ixyz,iat),fds1(1,1,ixyz),ncf,ntri,d1fds1)
            dewe=d1f1-d1fds1
            hess(ixyz,iat,ixyz,iat)=hess(ixyz,iat,ixyz,iat)+dewe
            do jxyz=ixyz+1,3
                call spur(d1(1,ixyz,iat),f1(1,jxyz),ncf,d1f1)
                call spur(d1(1,jxyz,iat),f1(1,ixyz),ncf,f1d1)
                call spuX(d1(1,ixyz,iat),fds1(1,1,jxyz),ncf,ntri,ws1)
                call spuX(d1(1,jxyz,iat),fds1(1,1,ixyz),ncf,ntri,sw1)
                df=d1f1+f1d1
                ws=ws1+sw1
                dewe=df-ws
                dewe= dewe*0.5d0
                hess(ixyz,iat,jxyz,iat)=hess(ixyz,iat,jxyz,iat)+dewe
            enddo
         enddo
        endif    !    (natend(iat).eq.0) then
      ENDDO       !   DO IAT=1,NATONCE
c---------------------------------------------------------------------
      if(natonce.eq.1) go to 4321
c---------------------------------------------------------------------
c
c  Iat.NE.Jat :
c
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
         irec=(iat-1)*3 + lrec
         call read1mat(60,irec+1,ncf2,fds1(1,1,1) )
         call read1mat(60,irec+2,ncf2,fds1(1,1,2) )
         call read1mat(60,irec+3,ncf2,fds1(1,1,3) )
         call read1mat(61,natoms+iat,ntri*3,f1(1,1) )
         DO JAT=IAT+1,NATONCE
           if(natend(jat).eq.0) then
            do ixyz=1,3
               do jxyz=1,3
                  call spur(d1(1,jxyz,jat),f1(1,ixyz),ncf,f1d1)
                  call spuX(d1(1,jxyz,jat),fds1(1,1,ixyz),ncf,ntri,sw1)
                  df=f1d1
                  ws=sw1
                  dewe=df-ws
                  dewe= dewe*0.5d0
                  hess(ixyz,iat,jxyz,jat)=hess(ixyz,iat,jxyz,jat)+dewe
               enddo
            enddo
           endif     ! (natend(jat).eq.0) then
         ENDDO       !   DO JAT=,IAT+1,NATONCE
        endif    !    (natend(iat).eq.0) then
      ENDDO       !   DO IAT=1,NATONCE
c
      DO IAT=1,NATONCE
        if(natend(iat).eq.0) then
         DO JAT=IAT+1,NATONCE
           if(natend(jat).eq.0) then
            irec=(jat-1)*3 + lrec
            call read1mat(60,irec+1,ncf2,fds1(1,1,1) )
            call read1mat(60,irec+2,ncf2,fds1(1,1,2) )
            call read1mat(60,irec+3,ncf2,fds1(1,1,3) )
            call read1mat(61,natoms+jat,ntri*3,f1(1,1) )
            do ixyz=1,3
               do jxyz=1,3
                  call spur(d1(1,ixyz,iat),f1(1,jxyz),ncf,d1f1)
                  call spuX(d1(1,ixyz,iat),fds1(1,1,jxyz),ncf,ntri,ws1)
                  df=d1f1
                  ws=ws1
                  dewe=df-ws
                  dewe= dewe*0.5d0
                  hess(ixyz,iat,jxyz,jat)=hess(ixyz,iat,jxyz,jat)+dewe
               enddo
            enddo
           endif     ! (natend(jat).eq.0) then
         ENDDO       !   DO JAT=,IAT+1,NATONCE
        endif    !    (natend(iat).eq.0) then
      ENDDO       !   DO IAT=1,NATONCE
c---------------------------------------------------------------------
 4321 continue
c---------------------------------------------------------------------
c make full hessian out of its upper triangle part :
c
      call hess_full_up(hess,natonce)
c
c---------------------------------------------------------------------
c
c       do iat=1,natonce
c          do jat=1,natonce
c             write(6,*) ' Atoms ',iat,jat
c        write(6,66) ((hess(i,iat,j,jat),i=1,3),j=1,3)
c          enddo
c       enddo
c
c
c 66  format(1x,3(e12.6,1x))
c 66  format(1x,9(f9.5,1x))
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine calcF1FDS1(ncf,    ntri,   fd,     s1,   f1,
     $                      fds1)
      implicit real*8 (a-h,o-z)
      dimension s1(ntri,3),f1(ntri,3)
      dimension fd(ncf,ncf)
      dimension fds1(ncf,ncf,3)
c------------------------------------------------------------------------
c Input :
c ncf     - number of basis functions
c ntri    - ncf*(ncf+1)/2
c fd      - Fock*Density matirx (precalculated and store on 60)
c s1      - 1st-order overlap
c f1      - 1st-order Fock
c
c Output:
c
c FDS1    - Fock1- Fock*Density*Overlap1 !written on disk , unit 60)
c------------------------------------------------------------------------
c
c in order to save memory calculate first FDS1
c
      do i=1,ncf
         ii=i*(i-1)/2
         do j=1,ncf
            jj=j*(j-1)/2
            if(i.ge.j) then
               ij=ii+j
            else
               ij=jj+i
            endif
            sum1=0.d0
            sum2=0.d0
            sum3=0.d0
            do k=1,ncf
               kk=k*(k-1)/2
               if(k.ge.j) then
                  kj=kk+j
               else
                  kj=jj+k
               endif
               sum1=sum1+fd(i,k)*s1(kj,1)
               sum2=sum2+fd(i,k)*s1(kj,2)
               sum3=sum3+fd(i,k)*s1(kj,3)
            enddo !  k=1,ncf
c
            fds1(i,j,1)=f1(ij,1)-sum1
            fds1(i,j,2)=f1(ij,2)-sum2
            fds1(i,j,3)=f1(ij,3)-sum3
c
         enddo !  j=1,ncf
      enddo !  i=1,ncf
c
c  and then S1DF and add it to FDS1
c
      do i=1,ncf
         ii=i*(i-1)/2
         do j=1,ncf
            jj=j*(j-1)/2
            sum1=0.d0
            sum2=0.d0
            sum3=0.d0
            do k=1,ncf
               kk=k*(k-1)/2
               if(k.ge.j) then
                  kj=kk+j
               else
                  kj=jj+k
               endif
               sum1=sum1+fd(i,k)*s1(kj,1)
               sum2=sum2+fd(i,k)*s1(kj,2)
               sum3=sum3+fd(i,k)*s1(kj,3)
            enddo !  k=1,ncf
c
            fds1(j,i,1)=-sum1+fds1(j,i,1)
            fds1(j,i,2)=-sum2+fds1(j,i,2)
            fds1(j,i,3)=-sum3+fds1(j,i,3)
c
         enddo !  j=1,ncf
      enddo !  i=1,ncf
c
      end
c======================================================================
      subroutine rescale_nat_fock(fock,ntri,natoms,lind,factor)
c scale multiple fock matrices during hessian calc in  their
c untransposed form (3,natonce,ntri)
c factor is thresh for uhf, thresh/2 for RHF
c
      implicit real*8 (a-h,o-z)
      dimension fock(3*natoms,ntri)
      dimension lind(*)       
c     
      call getival('ncf',ncf)      
c     
      call dscal(3*natoms*ntri,factor,fock,1)
      do i=1,ncf
         ii=lind(i)+i
         call dscal(3*natoms,2.0d0,fock(1,ii),1)
      enddo
c
      end
c======================================================================
      

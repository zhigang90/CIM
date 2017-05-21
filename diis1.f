      Subroutine Diis(iaction,errvecA,errvecB,fockA,  fockB,
     $                DenA,   DenB,   store,  ndim,   interp,
     $                rhf,    IPRNT,  xlam,   ndii)
      implicit real*8 (a-h,o-z)
c
c  Pulay DIIS procedure
c
c  ARGUMENTS
c
c  iaction -  action flag
c              1 - initialize the DIIS procedure
c              2 - close and remove DIIS files
c              0 - normal DIIS cycle
c             -1 - reduce dimension of DIIS subspace by ndiis
c  errvecA -  error vector (this is the commutator-like quantity [FDS-SDF]
c              for alpha spin (or for the closed shell case)
c  errvecB -  error vector for beta spin
c  fockA   -  alpha/closed-shell Fock matrix
c  fockB   -  beta Fock matrix
c  DenA    -  alpha/closed-shell density matrix
c  DenB    -  beta density matrix
c  store   -  alpha/closed-shell work array (3*ndim)
c  ndim    -  total dimension of matrices
c             (stored as lower triangle arrays)
c  interp  -  logical flag for interpolation
c  rhf     -  logical flag for close/open-shell system
c  IPRNT   -  print flag
c  xlam    -  final squared residue (DIIS error)
c  ndii    -  on input MAY request reduction of size of DIIS subspace
c             on output final dimension of DIIS subspace
c
c  References
c  ----------
c    P. Pulay, Chem.Phys.Letts. 73 (1980) 393
c    P. Pulay, J.Comp.Chem.  3 (1982) 556
c    T. P. Hamilton and P. Pulay,  J.Chem.Phys.  84 (1986) 5728
c
c
c  The first call to <Diis> MUST have iaction=1 to initialize
c  the scratch files
c
c  changes for open-shell implementation            JB  march 1998
c  -------------------------------------
c  minimal change -  UHF Fock and Density matrices written/read
c  on direct access file junit at same time as closed-shell matrices
c  from unit iunit
c
c
      dimension errvecA(ndim),fockA(ndim),DenA(ndim),store(3*ndim)
      dimension errvecB(ndim),fockB(ndim),DenB(ndim)
      logical interp,rhf
c
      parameter (nmaxd=24,nmaxd1=nmaxd+1,nmaxd2=nmaxd1**2)
      parameter (zero=0.0d0,one=1.0d0,tol=1.0d-16)
      parameter (iunit=15,junit=16)
c
      integer counter(nmaxd)
      logical exst
      character*84 dfile
      dimension bmtx(0:nmaxd,0:nmaxd),
     1  dinv(0:nmaxd2),ll(nmaxd1),mm(nmaxd1)
c
      Save icur,ndiis,ifirst,counter,bmtx,ndim1
c
c
      call getival('iout',iout)
      if(iaction.ne.0) go to 100
      if(ndim.ne.ndim1)then
         call nerror(1,'diis',
     1   'diis is called with a different dimension',
     2    ndim,ndim1)
      end if
c
c  counter is a circular counter. Its elements correspond to the
c  nmaxd slots on the external DIIS file on iunit (and junit for UHF)
c  if i is the i-th error + coefficient vector on the disk, then
c  counter(i) gives the temporal sequence number of these vectors
c  among the last (at most nmaxd) iterates
c  ifirst is the oldest element in counter
c  icur is the latest
c  ndiis is the number of the active elements
c  these counters are set initially by a call to diisinit
c
      ndiis=ndiis+1
c -- experimental to prevent overflow   ! JB April 98
      If(ndiis.ge.nmaxd)
     $   call reducediis(ndiis,counter,ifirst,nmaxd,bmtx)
c ---------------------------------------------
      icur=icur+1
      if(icur.gt.nmaxd) then
        icur=icur-nmaxd
      end if
      counter(icur)=ndiis
      if(icur.eq.ifirst.and.ndiis.gt.1) then
        call reducediis(ndiis,counter,ifirst,nmaxd,bmtx)
      end if
      write(iunit,rec=icur) errvecA,fockA,DenA
      if(.not.rhf) write(junit,rec=icur) errvecB,fockB,DenB
      bmtx(ndiis,ndiis)=ddot(ndim,errvecA,1,errvecA,1)
      if(.not.rhf) then
      bmtx(ndiis,ndiis)=bmtx(ndiis,ndiis)+ddot(ndim,errvecB,1,errvecB,1)
      end if
      ii=ifirst-1
      do i=1,ndiis-1
        ii=ii+1
        if(ii.gt.nmaxd) ii=1
        jj=counter(ii)
        read(iunit,rec=ii) store
        bmtx(ndiis,jj)=ddot(ndim,errvecA,1,store,1)
        bmtx(jj,ndiis)=bmtx(ndiis,jj)
      end do
      if(.not.rhf) then
        ii=ifirst-1
        do i=1,ndiis-1
          ii=ii+1
          if(ii.gt.nmaxd) ii=1
          jj=counter(ii)
          read(junit,rec=ii) store
          bmtx(ndiis,jj)=bmtx(ndiis,jj)+ddot(ndim,errvecB,1,store,1)
          bmtx(jj,ndiis)=bmtx(ndiis,jj)
        end do
      end if
c  now carry out the interpolation if instructed so
c  first copy the DIIS matrix & normalize the rows
c  the inverse diis matrix is defined differently because osinv
c  needs a closely packed matrix
      ndii=ndiis
      if(.not.interp) return
c  return address in case it is necessary to reduce the diis matrix
 500  continue
      ndii=ndiis
      if(ndiis.le.1) return
      call diisinv(nmaxd,ndiis,bmtx,dinv,tol,det,ll,mm)
      xlam=dinv(0)
      if(IPRNT.gt.2) then
        write(iout,300) det,(dinv(k),k=1,ndiis)
 300    format('det=',e12.4,' coefs=',/,5(5f12.7,/))
      end if
c  determine the maximum DIIS coefficient;reduce DIIS if it is too large
      call absmax(ndiis,dinv(1),ii,cmax)
      if(cmax.gt.30.0d0.or.abs(det).le.tol) then
        call reducediis(ndiis,counter,ifirst,nmaxd,bmtx)
        go to 500
      end if
c interpolation
      call zeroit(errvecA,ndim)
      call zeroit(fockA,ndim)
      call zeroit(DenA,ndim)
c -- unresticted -----------------
      If(.not.rhf) Then
        call zeroit(errvecB,ndim)
        call zeroit(fockB,ndim)
        call zeroit(DenB,ndim)
      EndIf
c --------------------------------
c  now perform the linear combination using the DIIS coefficients which are
c  in dinv(jj)
      ndim2=2*ndim
      ii=ifirst-1
      do nd=1,ndiis
        ii=ii+1
        if(ii.gt.nmaxd) ii=1
        jj=counter(ii)
        read(iunit,rec=ii) store
        call daxpy(ndim,dinv(jj),store,1,errvecA,1)
        call daxpy(ndim,dinv(jj),store(ndim+1),1,fockA,1)
        call daxpy(ndim,dinv(jj),store(ndim2+1),1,DenA,1)
c -- unrestricted -----------------
        If(.not.rhf) Then
          read(junit,rec=ii) store
          call daxpy(ndim, dinv(jj),store,1,errvecB,1)
          call daxpy(ndim,dinv(jj),store(ndim+1),1,fockB,1)
          call daxpy(ndim,dinv(jj),store(ndim2+1),1,DenB,1)
        EndIf
c ----------------------------------
      end do
      ndii=ndiis
      return
c
c  this is the initialization/closing/reducing part
 100  continue
c     if iaction=  1 then initialize DIIS
c     if iaction > 1 then close DIIS
c     if iaction= -1, reduce the DIIS space by ndii
c
      if(iaction.eq.-1) go to 200
      ndim1=ndim
      call getchval('scrf',dfile)
      call rmblan(dfile,80,len)
      exst=.false.
      inquire(iunit,exist=exst)
      if(exst) close(iunit,status='delete')
c -- unrestricted -----------------
      inquire(junit,exist=exst)
      if(exst) close(junit,status='delete')
c ---------------------------------
c
      if(iaction.eq.1) then
c  the length of a block is given in bytes.
        call getival('double',idoub)  ! size of double-precision word in bytes
c
c -- define and open alpha/closed-shell direct-access file
        leng=3*ndim1*idoub
        open(iunit,file=dfile(1:len)//'.dii',form='unformatted',
     1   access='direct',recl=leng)
c -- if UHF define and open beta direct-access file --------
        if(.not.rhf) then
          leng = 3*ndim1*idoub
          open(junit,file=dfile(1:len)//'.uii',form='unformatted',
     1     access='direct',recl=leng)
        end if
c ----------------------------------------------------------
c
c -- initialize
        ifirst=1
        icur=0
        ndiis=0
c  set the constant part of the DIIS matrix
        bmtx(0,0)=zero
        do i=1,nmaxd
          bmtx(i,0)=-one
          bmtx(0,i)=-one
        end do
      end if
      ndii=ndiis
      return
c
c  this part only reduces the DIIS space
 200  continue
      do ii=1,ndii
        call reducediis(ndiis,counter,ifirst,nmaxd,bmtx)
      end do
      end
c++++++++++++++++++++++++++++++++++++++++
      subroutine diisinv(nmaxd,ndiis,bmtx,dinv,tol,det,ll,mm)
      implicit real*8 (a-h,o-z)
c
c  This routine carries out the actual solution of the DIIS equations
c  The principal reason for separating it from the main routine is that
c  the DIIS matrix could be redefined with the actual dimension ndiis,
c  rather than the fixed dimension
c
c  ARGUMENTS
c
c  nmaxd   -  maximum dimension of DIIS subspace
c  ndiis   -  current dimension of DIIS subspace
c  bmtx    -  DIIS error matrix
c  dinv    -  on exit contains inverse of <bmtx>
c  tol     -  threshold below which DIIS matrix is considered singular
c  det     -  determinant of DIIS matrix
c  ll      -  integer scratch storage for <osinv>
c  mm      -    ditto
c
c  NOTES: The rows of the DIIS matrix are normalized
c         The DIIS coefficients are contained in the first column of
c         the inverse DIIS matrix, from 1..ndiis. dinv(0,0) is the
c         projected residue square after extrapolation
c
      dimension bmtx(0:nmaxd,0:nmaxd),dinv(0:ndiis,0:ndiis),
     1  ll(*),mm(*)
c
      parameter(one=1.0d0)
c
      do j=0,ndiis
        do i=0,ndiis
          dinv(i,j)=bmtx(i,j)
        end do
      end do
      do i=1,ndiis
        xx=one/dinv(i,i)
        do j=0,ndiis
          dinv(i,j)=dinv(i,j)*xx
        end do
      end do
      call osinv(dinv(0,0),ndiis+1,det,tol,ll,mm)
c  the first column of the inverse DIIS matrix gives the estimated
c  residue after interpolation & the interpolation coefficients
      call dscal(ndiis+1,-one,dinv(0,0),1)
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine reducediis(ndiis,counter,ifirst,nmaxd,bmtx)
      implicit real*8 (a-h,o-z)
c
c  This routine reduces the dimension of the DIIS procedure, either
c  because of numerical dependencies or because of the overflow
c  of the DIIS storage
c
c  ARGUMENTS
c
c  ndiis   -  current dimension of DIIS subspace
c             will be reduced by 1
c  counter -  circular counter for DIIS access
c  ifirst  -  position of the first (oldest) DIIS element
c             in the circular buffer
c  nmaxd   -  maximum DIIS dimension (unchanged)
c  bmtx    -  DIIS error matrix,(will be shifted)
c
      integer counter(nmaxd)
      dimension bmtx(0:nmaxd,0:nmaxd)
c
c   it does not make sense to reduce it below 1 - use "call diis"
c   with iaction=1
        if(ndiis.lt.1) return
            do i=ifirst+1,ifirst+ndiis-1
              ii=i
              if(ii.gt.nmaxd) ii=ii-nmaxd
              counter(ii)=counter(ii)-1
            end do
            ifirst=ifirst+1
            if(ifirst.gt.nmaxd) ifirst=1
            do i=1,ndiis-1
              do j=0,ndiis
       bmtx(j,i)=bmtx(j,i+1)
          end do
            end do
            do i=0,ndiis
              do j=1,ndiis-1
       bmtx(j,i)=bmtx(j+1,i)
          end do
            end do
          ndiis=ndiis-1
      end

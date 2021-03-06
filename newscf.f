      subroutine seigvec(smat,n,sthre,upper,uinv,xlow,newn)
c
      use memory
c
      implicit real*8(a-h,o-z)
c This subroutine diagonalizes smat and returns its factors (below)
c     Matrix must be symmetric and positive definite.
C  Arguments:
C  INTENT(IN)
c     smat: name of the  S matrix (character var.)
c     smat must be packed upper triangle stored by columns
c     order (11) (12) (22) (13) (23) (33) (14)...  or, what is the same,
c           (11) (21) (22) (31) (32) (33) ....
c     n is the dimension
c     sthre is the theshold for throwing out eigenvectors of S
C  INTENT(OUT)
c     upper is the name of a matrix; on exit it contains L**(1/2)*V(T)
c       where V are the eigenvectors of S  and (T)=transpose
c     Uinv is the NAME of the inverse of U, Uinv=V*L**(-1/2)
c     xlow is the lowest eigenvalue od Smat that exceeds sthre
c     newn=the number of non-redundant basis functions
      character*(*) smat, upper, uinv
      parameter (zero=0.0d0)
CPP
c      write(6,*) smat, upper,uinv
c      call matprint(smat,6)
      call matinfo(upper,mshape,idim1,idim2,iupperaddr,ilen)
      call matinfo(uinv,mshape1,idim1,idim2,iuinvaddr,ilen1)
        if(mshape.ne.2.or.mshape1.ne.2) then
            call nerror(1,'Cholesky','Matrix U or Uinv is not square
     1',mshape,mshape1)
        end if
        if(ilen.ne.n*n.or.ilen1.ne.n**2) then
          call nerror(2,'Cholesky','Matrix U .or Uinv is not the right
     1    size',ilen, ilen1)
        end if
c      imatrixaddr=mataddr(smat)
c     create a temporary matrix to hold intermediate results
c
c   calculate the eigenvalues of the matrix by diagonalization
      call getmem(n,ifail)
      call getmem(n,ieval)
c      call matdef('temp','v',n*(n+1)/2,1)
      call matdef('temp','s',n,n)
      call getmem(8*n,itempx)
      call getmem(5*n,itempi)
      itempaddr=mataddr('temp')
      call matcopy(smat,'temp')
c      call tfer(bl(imatrixaddr),bl(itempaddr),n*(n+1)/2)
      abstol=2.0d0*dlamch('S')
CPP
c      print *, 'abstol=',abstol
c      call matprint('temp',6)
CPP
      call dspevx('V','A','U',n,bl(itempaddr),VLdum,VUdum,1,1,
     1 abstol,mfound,bl(ieval),bl(iupperaddr),n,bl(itempx),bl(itempi),
     1 bl(ifail), info)
c  Now the eigenvectors of S, SV=VL are in U=V
      call retmem(2)
      call matrem('temp')
CPP
c      call matconn('sevec','d',n,n,ieval)
c      call matprint('sevec',6)
c      call matdisc('sevec')
c      call matprint(upper,6)
CPP
      m=0   ! First eigenvalue exceeding sthre
      do i=1,n
        ss=bl(ieval-1+i)
c        print *, 'loop',i,ss
        if(m.eq.0.and.ss.gt.sthre) then
          m=i
          xlow=ss
        end if
        bl(ifail-1+i)=sqrt(ss)
        if(ss.gt.1.0d-15) then
          bl(ieval-1+i)=1.0d0/ss
        else
          bl(ieval-1+i)=0.0d0
        end if
      end do
c   now bl(ifail) contains lambda**(+1/2). bl(ieval) 1/lambda
c      print *, 'After eigenvalue manipulation'
c      call matprint(upper,6)
      call matconn('sevec','d',n,n,ifail)
CPP
c      call matprint('sevec',6)
c      call matprint(upper,6)
      call matmmult(upper,'sevec',uinv)
C  uinv contains (temporarily) V*L**(+1/2)
c       print *, 'U (not Uinv) full'
c      call matprint(uinv,6)
      call matdisc('sevec')
C  truncate U and Uinv
      newn=n
      if(m.gt.1) then
c   Redundant functions
c  simpler if the eigenvalues would be ordered in descending order
        newn=n-m+1
        nm=n*newn  ! number of elements in the truncated matrix
        iainv=mataddr(uinv)     
        istart=iainv+(m-1)*n
c       call tfer(bl(istart),bl(iainv),nm)
        do i=0,nm-1
          bl(iainv+i)=bl(istart+i)
        enddo
c        print *, 'U transpose after shift'
c        call matprint(uinv,6)
        call matredef(uinv,uinv,'r',n,newn)
        call matredef(upper,upper,'r',newn,n) 
        call matpose2(uinv,upper,'n')
        call matconn('sevec','d',newn,newn,ieval+m-1)
        call matmmult(uinv,'sevec',uinv)
        call matdisc('sevec')
c        call matprint(uinv,6)
      else
        call matpose2(uinv,upper,'n')
        call matconn('sevec','d',n,n,ieval)
        call matmmult(uinv,'sevec',uinv)
        call matdisc('sevec')
C   upper = (V*L**1/2)T;   uinv = V*L(-1/2)
      end if  
CPP
      call retmem(2)
CPP
c      print *, 'Final values'
c      call matprint(upper,6)
c      call matprint(uinv,6)
c      write(6,*) 'xlow=',xlow,'  m=',m
CPP
c
      end
c=======================================================================
      subroutine cholesky1(matrix,n,upper)
c
      use memory
c
      implicit real*8(a-h,o-z)
      integer*4 i4nfo
c This subroutine performs Cholesky factorization on matrix and returns
c     the upper triangular factor upper.  Matrix must be symmetric and
c     positive definite.
C  Arguments:
c     matrix: name of the  S matrix (character var.)
c     n is the dimension of the matrix
c     upper: name of the upper Cholesky factor (character var.)
c     matrix must be packed upper triangle styored by columns
c     order (11) (12) (22) (13) (23) (33) (14)...  or, what is the same,
c           (11) (21) (22) (31) (32) (33) ....
c     upper is a quadratic matrix currently because sminhalf yields a
c      quadratic one
      character*(*) matrix, upper
      parameter (zero=0.0d0)
c First create upper if necessary and get the addresses of the matrices.
      call matsear(upper,jcur)
      if(jcur.eq.0) then
         call matdef(upper,'q',n,n)
         iupperaddr=mataddr(upper)
      else
         call matinfo(upper,mshape,idim1,idim2,iupperaddr,ilen)
         if(mshape.ne.2) then
            call nerror(1,'Cholesky','Target matrix is not square
     1',mshape,0)
         end if
         if(ilen.ne.n*n) then
            call nerror(2,'Cholesky','Target matrix is not the right
     1      size',ilen, n*n)
         end if
      end if
      imatrixaddr=mataddr(matrix)
c     create a temporary matrix to hold intermediate results
      call matdef('temp','v',n*(n+1)/2,1)
      itempaddr=mataddr('temp')
c
c     Transfer matrix into temp (again because dspevx destroys it)
      call tfer(bl(imatrixaddr),bl(itempaddr),n*(n+1)/2)
c     do the Cholesky factorization
      call dpptrf('U',n,bl(itempaddr),i4nfo)
      info=i4nfo
c
      call getival('iout',iout)
c     error checking
      if(info.lt.0) then
         call nerror(3,'Cholesky','Invalid argument in
     1 lapack routine dpptrf : number ',info,0)
      end if
      if(info.gt.0) then
         call nerror(4,'Cholesky','Matrix is not
     1positive definite',info,0)
      end if
c     transfer the result to the matrix upper
      call packed_to_upper(bl(itempaddr),bl(iupperaddr),n)
      call matrem('temp')
      end
c=======================================================================
      subroutine geneig1(fock,u,uinv,diag,coef,dens,xlvsh,nroots,upper)
c
      use memory
c
      implicit real*8 (a-h,o-z)
      character*(*) fock,u,uinv,diag,coef,dens,upper
      character*4 temp,tmp1
      integer*4 info
c  this routine solves the generalized eigenvalue
c  equation FC = SCe   The matrices U and Uinv should satisfy
c  Univ(t)SUinv = I(identity) or S=U(t)U where (t) is transpose;
c  Uinv is the inverse of U.
c  It solves the equation (Uinv(t)FUinv(UC)=(UC)e
c  F =fock, u=U, uinv=Uinv, coef=C, diag=e
c  input data: fock,u,uinv,dens,xlvsh dens is needed for the level shift
c               nroots: number of roots desired. If nroots=0, all
c               roots are calculated, otherwise the lowest nroots
c              'upper' is a character variable, either 'U' or 'L'
c              it is 'u' if the matrix 'u' is an upper triangle,
c              and it is 'L' if it is a lower triangle. this has been
c              added so that it can handle the reverse case which
c              is needed somewhere in the NMR program
c  output: coef,diag
c  note that xlvsh must be half of the intended level shift
c  because it is multiplied by the density at double occupancy
c
c using lapack generalized eigensolver routines would be even faster
c just the levelshift should be applied in AO basis
c
c----------test------------
c      call matprint(fock, 6  )
c      call matprint(u   , 6  )
c      call matprint(uinv, 6  )
c----------test------------
      call getival('ncf',ncf)
      zero=rgetrval('zero')
      one=rgetrval('one')
      temp='temp'
      tmp1='tmp1'
      call matinfo(u,mshape,id1,id2,maddr,length)
      call matdef(temp,'s',id1,id1)
      iadra=mataddr(temp)
      iadrw=mataddr(diag)
      iadrc=mataddr(coef)
      iadru=mataddr(u)
      iadrui=mataddr(uinv)
cc calculate Uinv(t)FUinv 
c Do not assume that Uinv is triangular
      call matsimtr(fock,uinv,temp)
c
c      xlvsh=zero
      if(xlvsh.ne.zero) then
c calculate U*D*U(t)  do not assume U being triangular
        call matdef(tmp1,'r',id1,ncf)
        call matmmult(u,dens,tmp1)
        call matscal(tmp1,-xlvsh)
c        call matinfo(tmp1,mshape,nrow,ncol,maddr,length)
c        print *,'tmp1, shape=',mshape, nrow,ncol,length
c        call matinfo(u,mshape,nrow,ncol,maddr,length)
c        print *,'u, shape=',mshape, nrow,ncol,length
        call matmmul2(tmp1,u,temp,'n','t','a')
        call matrem(tmp1)
c  the above operation subtracts U*D*U(t) from Uinv(t)*F*Uinv
      end if
c
c sdiag2 is a lot slower - use always LAPACK
c     if(nroots.eq.0) then
c       call matdiag(temp,diag,coef)
c     else
      if(nroots.eq.0) then        ! get all eigenvalues
         nroo=id1
      else                        ! get first nroots eigenvalues
         nroo=nroots
      end if
      if(nroo.gt.id1)nroo=id1
      call matdef('coef1','q',id1,id1)
      iadrc1=mataddr('coef1')
      call getmem(8*ncf,iadrtmp)
      call getmem(5*ncf,iwork)
      call getmem(ncf,iadrfail)
      abstol=2.0d0*dlamch('S')
      call dspevx('V','I','U',id1,bl(iadra),zero,zero,1,nroo,
     1            abstol,m,bl(iadrw),bl(iadrc1),id1,bl(iadrtmp),
     2            bl(iwork),bl(iadrfail),info)
c      call dspevx('V','I','U',id1,bl(iadra),zero,zero,1,nroo,
c     1            1.0d-10,m,bl(iadrw),bl(iadrc1),id1,bl(iadrtmp),
c     2            bl(iwork),bl(iadrfail),info)
        call retmem(3)
c calculate Uinv*Coeff taking advantage of Uinv being triangular
c      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iadrui),ncf,
c     *           bl(iadrc1),ncf)
      call matmmult(uinv,'coef1',coef)
      call matrem('coef1')
      call matrem(temp)
      end
c=======================================================================
      subroutine lowesteigv(matrix,n,xlow)
c
      use memory
c
      implicit real*8(a-h,o-z)
      integer*4 i4nfo
c This subroutine determines the lowest eigenvalue of a matrix, 
c in general the overlap matrix
C  Arguments:
c     matrix: name of the  S matrix (character var.)
c     matrix must be packed upper triangle styored by columns
c     order (11) (12) (22) (13) (23) (33) (14)...  or, what is the same,
c           (11) (21) (22) (31) (32) (33) ....
c     n is the dimension
c     xlow is the smallest diagonal element of the matrix
      character*(*) matrix
      parameter (zero=0.0d0)
c
      imatrixaddr=mataddr(matrix)
c     create a temporary matrix to hold intermediate results
      call matdef('temp','v',n*(n+1)/2,1)
      itempaddr=mataddr('temp')
c
c   calculate the lowest eigenvalue of the matrix by diagonalization
c    This is not bad, as the vectors are not needed
      call getmem(8*n,itempx)
      call getmem(5*n,itempi)
      call getmem(n,ifail)
      call getmem(n,ieval)
      call tfer(bl(imatrixaddr),bl(itempaddr),n*(n+1)/2)
      call dspevx('N','I','U',n,bl(itempaddr),VLdum,VUdum,1,1,
     1 1.0d-10,mfound,bl(ieval),Zdummy,n,bl(itempx),bl(itempi),
     1 bl(ifail), i4nfo)
c   xlow is the lowest eigenvalue of the S matrix
      xlow=bl(ieval)
      call retmem(4)
      call matrem('temp')
c
      end
c=======================================================================
      subroutine converged1(u,uinv,density,fock,commut,error,errsq)
c
      use memory
c
      implicit real*8 (a-h,o-z)
c  ARGUMENTS:
c  INTENT(IN)
c  U and Uinv are decompositions of the overlap matrix S such that
c     U(t)U = S Uinv is the inverse of U; U may be S**(1/2),
c     though here it is the Cholesky decomposition of S
c  density is the name of the density matrix D(max. 8 char.)
c  fock is the name of the fock matrix F
c  commut is the commutator; it MUST be defined as antisymmetric,
c  otherwise the results will be wrong
c  INTENT(OUT)
c  error is the norm of the commutator UDFUinv - (UDFUinv)transpose
c
      character*(*) u,uinv,density,fock,commut
      character*4 imed
c  Calculate UDFUinv - its transpose
      call getival('ncf ',ncf)
      ncf2=ncf*(ncf+1)/2
c
      iu=mataddr(u)
      iuinv=mataddr(uinv)
      icommut=mataddr(commut)
      imed='imed'
      call matinfo(u,mshape,id1,id2,maddr,length)
c      print *, 'Dimension of u', id1,id2
      call matdef(imed,'q',ncf,ncf)
      iimed=mataddr(imed)
c
c this is only needed for the intel blas library
c that does not conform to the standard and does not zero out the result
c problems surfaced in the MPI code that used uninitialized memory
      call matzero(imed)
c
c use triangular matrix multiply with U and Uinv
      call matmmult(density,fock,imed)
c      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iu),ncf,
c     *           bl(iimed),ncf)
c      call dtrmm('R','U','N','N',ncf,ncf,1.0d0,bl(iuinv),ncf,
c     *           bl(iimed),ncf)
      call matdef('imed1','r',id1,id2)
      call matmmult(u,imed,'imed1')
      call matmmult('imed1',uinv,commut)
      call matrem('imed1')
c
      call matrem(imed)
      call absmax(ncf2,bl(icommut),iii,error)
      errsq=ddot(ncf2,bl(icommut),1,bl(icommut),1)
      end
c=======================================================================
      subroutine SymORB(nsym,   ncf,  nval,   fpair,  CMO,
     $                  S,      SC,   CSYM,   CSC,    EOrb)
c
c Symmetrize orbitals so that they become one dimensional basis vectors for
c irreps of an Abelian group. They can be mixed in the input "coeff",
c because a higher symmetry is present in the system. This subroutine
c should be able to handle up to five orbitals mixed, although testing
c has not been performed thoroughly for systems with very high symmetry.
c
c Input parameters:
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of orbitals being symmetrized.
c   fpair    - INTEGER(!) array of basis functions symmetry pairs
c   CMO      - full MO coefficient matrix
c   S        - overlap matrix
c   EOrb     - array of orbital energies.
c Scratch space:
c   SC       - SC product (overlap*coeff)
c   CSYM     - C` = Symmetry operation * C (coeff)
c   S        - the product C`SC will be held here
c
      use memory
      implicit none
      integer nsym,ncf,nval,fpair(7,ncf)
      real*8 CMO(ncf,nval),S(ncf,ncf),SC(ncf,nval)
      real*8 CSYM(ncf,nval),CSC(nval,nval),EOrb(nval)
c
      integer indexx,ns,io,jo,j,k,ii,jj,i,itry,ll,isize,kk,ns1
      real*8 xsum
      real*8 overprim(2,nval)
      integer, parameter :: max_block = 5
      real*8 values(max_block), SS(max_block*max_block)
      real*8 vectors(max_block*max_block),newco(ncf,max_block),el
      real*8, parameter :: eps = 1d-12    ! square of max nondiagonal element of S
      real*8, parameter :: deg_thr = 1d-4 ! degeneracy threshold
      real*8 :: eps1 = eps
      logical done
      integer iz,nn,maxc,mm,icf,iover,ires,icoepr,ixov
      real*8 xmax,xval
c
      If (nsym.eq.0) return
c
c calculate the largest non-diagonal element:
c
c Build S*C matrix
      do ns=1,nsym
        call DGEMM('N',    'N',    ncf,    nval,   ncf,
     $              1.0d0,  S,     ncf,    CMO,    ncf,
     $              0.0d0,  SC,    ncf)
c Build C' matrix, which is a symmetry image of C, under operation "ns"
        do io=1,nval
          do j=1,ncf
          indexx=fpair(ns,j)
          CSYM(abs(indexx),io)=CMO(j,io)
          if (indexx.lt.0)
     *        CSYM(abs(indexx),io)=-CSYM(abs(indexx),io)
          enddo
        enddo
c Calculate the overlap C'SC
        call DGEMM('T',    'N',    nval,   nval,   ncf,
     $              1.0d0,  CSYM,  ncf,    SC,     ncf,
     $              0.0d0,  CSC,   nval)
c Look for non-diagonal submatrices to be diagonalized (tricky task)
c They are blocks along diagonal
      ii=0
      DO jj=1,nval
        ii=ii+1
        if (ii.ge.nval) exit
 270    continue
        maxc=ii               ! Max diagonal element connected via
c                               nodiagonal with ii, init to ii itself
        do nn=ii,min(ii+20,nval)     ! n goes over columns
          !investigate column nn, below diag. AND maxc
          do kk=max(nn+1,maxc),min(ii+20,nval)
            el=CSC(kk,nn)*CSC(kk,nn)
            if (el.gt.eps1) then
              maxc=kk
            endif
          enddo
          if (maxc.le.nn) exit
          if (maxc.ge.ii+5) then
            if (eps1.gt.1d-7) then
              write(6,*)ii,maxc
              do mm=ii,min(maxc,nval)
                do kk=ii,min(maxc,nval)
                  write(6,'(E10.2E2)',ADVANCE='no') CSC(mm,kk)
                enddo
                write(6,*)
              enddo
              call nerror(1,'Subroutine SymORB',
     $                    'Too large representation detected',ii,maxc)
            else
              write(6,'(A,2I5,F15.10)')
     *         'Warning, increasing eps due to rep size: ',
     *           ii,maxc,dabs(EOrb(maxc)-EOrb(ii))
              eps1=eps1*2.0d0
              goto 270
            endif
          endif
          if (dabs(EOrb(maxc)-EOrb(ii)).gt.deg_thr) then
c            write(6,'(A,2I5,F15.10)')
c     *         'Warning, increasing eps due to energy difference: ',
c     *           ii,maxc,dabs(EOrb(maxc)-EOrb(ii))
            eps1=eps1*2.0d0
            goto 270
          endif
        enddo
        if (maxc.eq.ii) cycle
        isize=maxc-ii+1
c
        if (isize.eq.1) call nerror(2,'Subroutine SymORB',
     $                    'Rep size cannot be one here',ii,maxc)
        if (isize.gt.5) call nerror(3,'Subroutine SymORB',
     $                    'Rep size cannot be greater than 5',ii,maxc)
c
c Now we know that the block size is at least: isize
c
        do kk=1,isize
          do ll=1,isize
            SS(kk+(ll-1)*isize)=CSC(ii+kk-1,ii+ll-1)
          enddo
        enddo
c Diagonalize the submatrix:
        call sdiag2(isize,isize,SS,values,vectors)
c
c create new linear combination of orbitals:
c
        do i=1,isize
          do j=1,ncf
            newco(j,i)=0.0d0
            do k=1,isize
            newco(j,i)=newco(j,i)+CMO(j,ii+k-1)*vectors(k+(i-1)*isize)
            enddo
          enddo
        enddo
c
c Replace old orbitals by new ones
c
        do k=1,isize
          do j=1,ncf
            CMO(j,ii+k-1)=newco(j,k)
          enddo
        enddo
        ii=ii+isize-1
      EndDO
      enddo
c
      return
      end

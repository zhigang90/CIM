      subroutine geneig1(fock,u,uinv,diag,coef,dens,xlvsh,nroots,upper)
c
      use memory
c
      implicit real*8 (a-h,o-z)
      character*(*) fock,u,uinv,diag,coef,dens,upper
      character*4 temp,tmp1
      integer*4 i4nfo
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
     2            bl(iwork),bl(iadrfail),i4nfo)
c      call dspevx('V','I','U',id1,bl(iadra),zero,zero,1,nroo,
c     1            1.0d-10,m,bl(iadrw),bl(iadrc1),id1,bl(iadrtmp),
c     2            bl(iwork),bl(iadrfail),info)
        call retmem(3)
c calculate Uinv*Coeff taking advantage of Uinv being triangular
c      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iadrui),ncf,
c     *           bl(iadrc1),ncf)
      call matmmult('uinv','coef1','coef')
      call matrem('coef1')
      call matrem(temp)
      end

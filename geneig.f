      subroutine geneig(fock,u,uinv,diag,coef,dens,xlvsh,nroots,upper)

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) fock,u,uinv,diag,coef,dens,upper
      character*4 temp
c     common /big/bl(1000)
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
c introduced triangular matrix multiplies for speed
c using lapack generalized eigensolver routines would be even faster
c just the levelshift should be applied in AO basis
c
c----------test------------
c      write(6,*)' from geneig :'
c      call matprint(fock, 6  )
c      call matprint(u   , 6  )
c      call matprint(uinv, 6  )
c----------test------------
      call getival('ncf',ncf)
      zero=rgetrval('zero')
      one=rgetrval('one')
      temp='temp'
      call matdef(temp,'s',ncf,ncf)
      iadra=mataddr(temp)
      iadrw=mataddr(diag)
      iadrc=mataddr(coef)
      iadru=mataddr(u)
      iadrui=mataddr(uinv)
c calculate Uinv(t)FUinv taking advantage of Uinv being triangular
c  use 'coef' as a temporary matrix
      call matcopy(fock,coef)
      call dtrmm('L',upper,'T','N',ncf,ncf,one,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
      call dtrmm('R',upper,'N','N',ncf,ncf,one,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
      call matcopy(coef,temp)
c     call matsimtr(fock,uinv,temp)
c
      if(xlvsh.ne.zero) then
c calculate U*D*U(t)  taking advantage of U being triangular
c  use 'coef' as a temporary matrix
        call matcopy(dens,coef)
        call dtrmm('L',upper,'N','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matscal(coef,-xlvsh)
        call dtrmm('R',upper,'T','N',ncf,ncf,one,bl(iadru),ncf,
     *             bl(iadrc),ncf)
        call matadd(coef,temp)
c       call matmmult(u,dens,coef)
c       call matscal(coef,-xlvsh)
c       call matmmul2(coef,u,temp,'n','t','a')
c  the above operation subtracts U*D*U(t) from Uinv(t)*F*Uinv
      end if
c
c sdiag2 is a lot slower - use always LAPACK
c     if(nroots.eq.0) then
c       call matdiag(temp,diag,coef)
c     else
      if(nroots.eq.0) then        ! get all eigenvalues
         nroo=ncf
      else                        ! get first nroots eigenvalues
         nroo=nroots
      end if
      call getmem(8*ncf,iadrtmp)
      call getmem(5*ncf,iwork)
      call getmem(ncf,iadrfail)
      call dspevx('V','I','U',ncf,bl(iadra),zero,zero,1,nroo,
     1            1.0d-10,m,bl(iadrw),bl(iadrc),ncf,bl(iadrtmp),
     2            bl(iwork),bl(iadrfail),info)
        info=int(i4nfo)
        call retmem(3)
c calculate Uinv*Coeff taking advantage of Uinv being triangular
      call dtrmm('L','U','N','N',ncf,ncf,1.0d0,bl(iadrui),ncf,
     *           bl(iadrc),ncf)
c     call matmmult(uinv,coef,coef)
      call matrem(temp)
      end

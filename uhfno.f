      SUBROUTINE uhfNO(ncf,    NAlpha, NBeta,  dens,   densB,
     $                 olddA,  u,      uinv,   diag,   coef,
     $                 lprint, sexpv,  xmulti)

      use memory

      implicit real*8(a-h,o-z)
c
c  calculate UHF natural orbitals
c
c  ARGUMENTS
c
c  ncf     -  number of basis functions
c  NAlpha  -  number of occupied alpha MOs
c  NBeta   -  number of occupied Beta MOs
c  dens    -  alpha density matrix
c  densB   -  beta density matrix
c  olddA   -  work storage (for density sum)
c  u       -  S**1/2
c  uinv    -  S**-1/2
c  diag    -  storage for eigenvalues
c  coef    -  storage for UHF natural orbitals
c  lprint  -  print level; if > 1 then the fractionally
c               occupied NOs are printed
c
c  on exit
c
c  sexpv   -  expectation value <S**2>
c  xmulti  -  calculated spin multiplicity
c
c ** WARNING - This subroutine uses TEXAS matrix system **
      character*(*) dens,densB,olddA,u,uinv,diag,coef
c     common /big/bl(30000)
      parameter (Zero=0.0d0,Half=0.5d0,One=1.0d0,Two=2.0d0)
      parameter (np4=4)
c
c  Use the negative of the density matrix so that the largest
c  eigenvalues come out first
      call matcopy(dens,olddA)
      call matadd(densB,olddA)
      call geneig2(olddA,u,uinv,diag,coef,0)
c  restore the sign of diag
      call matscal(diag,-One)
c  write out these quantities
      call matwrite(coef,np4,0,'nato_uhf')
      call matwrite(diag,np4,0,'naoc_uhf')
c    idiag is the starting address of diag
      idiag=mataddr(diag)-1
      eps=0.005d0
      eps1=Two-eps
      nlow=0
      nhi=0
c  expectation value of <S**2>=[(na-nb)/2]**2 + n - Sum[de(i)**2]/2
c  The double of this quantity is calculated here first
      sexpv = Half*(NAlpha-Nbeta)**2 + Two*(NAlpha+Nbeta)
      do i=1,ncf
        if(bl(idiag+i).lt.eps1.and.nlow.eq.0) nlow=i
        if(bl(idiag+i).lt.eps.and.nhi.eq.0)   nhi=i-1
        sexpv=sexpv-bl(idiag+i)**2
      end do
c  The (fractional) multiplicity is sqrt(4*<S**2>+1); note that
c   at this point sexpv is twice its final value
      xmulti=sqrt(Two*sexpv+One)
      sexpv=sexpv*Half
c
      call setrval('<s2>',sexpv)     ! store in depository
c
      iout = igetival('iout')
      write(iout,*)
      write(iout,650) sexpv,xmulti
      If(nlow.gt.0.and.nlow.le.nhi) Then
        write(iout,700) eps,eps1
        do i=nlow,nhi
        write(iout,750) i,bl(idiag+i)
        enddo
        if(lprint.gt.0) then
          call matsub('natorb',coef,nlow,nhi)
          call druwf(bl(idiag+nlow),'natorb',0.0d0,0,0,nhi-nlow+1,
     1     iout)
          call matrem('natorb')
        end if
      EndIf
c
      return
c
 650  format(' Expectation value of S**2 and multiplicity=',2f11.7)
 700  format(' natural occupation numbers between ',f5.3,' and ',f5.3)
 750  format(1x,i4,2x,f10.6)
c
      end

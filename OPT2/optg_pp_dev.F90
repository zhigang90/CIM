 module optpp
 implicit double precision (a-h,o-z)
 include "common/big"
 private
 public optg_pp
 character(16) :: method     !optimization method
 logical       :: use_bmat   !if true, use bmat coordinates
 contains
!---------------------------------------------------------------------------------------------------------------
 subroutine optg_pp(ncen,npt,ityp,iout,atname,x,g,e,h,xnew,enew,iuniq, have_hess,do_hess, iprgrd,iprhes,iprhis)
!---------------------------------------------------------------------------------------------------------------
!
! This routine does a geometry optimization step. The previous geometries and gradients
! are provided in x and g, respectively. Optionally, a hessian is provided in h.
! The new geometry is returned in xnew
! The convengence check is done outside
!
 integer, intent(in)  :: ncen                        !number of centers
 integer, intent(in)  :: npt                         !number of optimization steps
 integer, intent(in)  :: iout                        !print output unit
 integer, intent(in)  :: ityp                        !optimization type: 1=minimum, 2=transition state
 integer, intent(in)  :: iuniq(ncen)                 !points to symmetry unique center for each center
 integer, intent(in)  :: iprhis                      !if positive, print optimization history
 integer, intent(in)  :: iprgrd                      !if positive, print current gradient
 integer, intent(in)  :: iprhes                      !if positive, print current hessian

 character(*), intent(in) :: atname(*)               !center name (for unique atoms)

 double precision, intent(in)  :: x(3,ncen,npt)      !history of coordinates
 double precision, intent(in)  :: g(3,ncen,npt)      !history of gradients
 double precision, intent(in)  :: e(npt)             !history of energies
 double precision, intent(in)  :: h(3*ncen,3*ncen)   !hessian (if have_hess=.true.)
 double precision, intent(out) :: xnew(3,ncen)       !new geometry
 double precision, intent(out) :: enew               !estimated next energy

 logical, intent(in)  :: have_hess                   !if true, hessian is available
 logical, intent(out) :: do_hess                     !if true, compute hessian in next step

 double precision, parameter :: thr=1.d-10           !threshold for zero hessian eigenvalues

 save nq,nprim,inpf3,inpf4,inpf5
 ibase=icorr(0)                                      !remember base address in ibase
 do_hess=.false.
        memxc=icorr(3*ncen)
!pp  In the first round, get options,generate the coordinates, and allocate memory which is needed later
      if(npt.eq.1) then
        ibase1=icorr(0)
        call get_inpi('PRINT', 'OPTG', ipri)
        write(iout,'(1x,A,I4)') "!xasd PRINT LEVEL = ", ipri
        memian=icori(ncen)
        call ConvertCart(ncen,iout,x(1,1,npt),iuniq,atname,  q(memxc),iq(memian))
! In the first round, or (later) restart, generate internal coordinates and write them to inpf3
        call coordgen(ncen,iout,x,iuniq,atname,  q(memxc),iq(memian),nq,nprim)
!  This routine returns two file numbers which are used alternatively to store internal coordinate definitions
        call corlsr(ibase1)
!  These file numbers have been shifted to COMMON /TAPENEU/ as inpf33 and inpf34
!  Define the file for the previous step
        inpf5=35
        call assgn(inpf5,'inthistory',isize,0)
      end if
!  Memory for the old internal coordinates
        memqold=icorr(nq)
!  Memory for the old forces
        memfold=icorr(nq)
!  Memory for the internal Hessian
        memhint=icorr(nq**2)
     if(npt.gt.1) then
          rewind inpf5
! Read in values in the previous step
          read(inpf5,100) (q(memqold+k-1),k=1,nq)
          read(inpf5,100) (q(memfold+k-1),k=1,nq)
          read(inpf5,200) (q(memhint+k-1),k=1,nq**2)
  100   format(f20.10)
  200   format(5f20.10)
!          write(6,*) 'Old internal coordinates upon entering MAIN,npt=',npt
!          write(6,'(5f12.6)') (q(memqold+k-1),k=1,nq)
!          write(6,*) 'Old internal forces upon entering MAIN, npt=',npt
!          write(6,'(5f12.6)') (q(memfold+k-1),k=1,nq)
!          call prntmatn(6,nq,nq,nq,q(memhint),'Input Hessian')
        end if
        memforc=icorr(3*ncen)
        memibcode=icori(6*nprim)
        memibcontr=icori(20*nq)
        memqq=icorr(nq)
        membmat=icorr(54*nq)
        memgmat=icorr(nq**2)
        memproj=icorr(nq**2)
        memphi=icorr(nq)
      call dcopy_x(3*ncen,x(1,1,npt),1,q(memxc),1)
      call dcopy_x(3*ncen,g(1,1,npt),1,q(memforc),1)
      call intforc1(ncen,q(memxc),q(memforc),nq,nprim, &
                    iq(memibcode),iq(memibcontr),q(memqq),q(membmat),q(memgmat),  q(memproj),q(memphi)) 
!pp
!      write(6,*) 'Internal coordinates in step',npt
!      write(6,'(5f12.6)') (q(memqq+k-1),k=1,nq)
!      write(6,*) 'Internal forces in step',npt
!      write(6,'(5f12.6)') (q(memphi+k-1),k=1,nq)
!      if(npt.gt.1) then
!        write(6,*) 'Previous internal coordinates'
!        write(6,'(5f12.6)') (q(memqold+k-1),k=1,nq)
!        write(6,*) 'Previous internal forces'
!        write(6,'(5f12.6)') (q(memfold+k-1),k=1,nq)
!      end if
!  Hessian preparation
        nek=3*ncen
!        write(*,*) 'Have_Hess=',have_hess
        if(have_hess.and.npt.eq.1) then
!Invert it and transform the inverse to redundant internal coordinates
          call TransformHess(nq,nek,h,q(memibcontr),q(membmat),q(memproj),    q(memhint))
        else if(npt.eq.1) then
!   read the force constants generated in genintcor and make initial inverse Hessian 
!         call prntmatn(6,nq,nq,nq,q(memproj),'Projector in the main program')
          call diagfcinv(nq,q(memhint),q(memproj))
        end if 
!  
!  Use BFGS to update the Hessian if npt >1
       if(npt.gt.1) then
!      write(6,*) 'Internal coordinates before entry to BFGS'
!      write(6,'(5f12.6)') (q(memqq+k-1),k=1,nq)
         call bfgs(nq,q(memqq),q(memqold),q(memphi),q(memfold),q(memhint))
       end if
!      write(6,*) 'Internal coordinates after BFGS'
!      write(6,'(5f12.6)') (q(memqq+k-1),k=1,nq)
!
       rewind inpf5
        write(inpf5,100) (q(memqq+k-1),k=1,nq)
        write(inpf5,100) (q(memphi+k-1),k=1,nq)
        write(inpf5,200) (q(memhint+k-1),k=1,nq**2)
!        write(6,*) 'internal coordinates and forces copied to inpf5'
!        write(6,'(5f12.6)') (q(memqq+k-1),k=1,nq)
!        write(6,'(5f12.6)') (q(memphi+k-1),k=1,nq)
!       
       memdelta=icorr(nq)
       memqq1=icorr(nq)
! PRINT is used as print level,value ipri
       coordmax=1.0d-1
       call relaxint(nq,q(memphi),q(memhint),iout,ipri,  &
                     coordmax,q(memdelta),edecr,dismx,icmax)
! Form current internal coordinates
       call add(q(memdelta),q(memqq),nq)
!  Write the current coordinates and forces to the "old" arrays
!       
       memb=icorr(nq*54)
       memnc=icorr(nek)
       pi=3.14159265358979324d0
       call intcvalues(ncen,q(memxc),pi*0.5d0,nqq,q(memqq))
!pp
!  write(*,*) 'In the main program before calling distortnew1, nqq=',nqq      
       call distortnew1(nek,nq,ncen,q(memdelta),q(memqq), q(memqq1),  &
                       q(memxc),iq(memibcode),iq(memibcontr),q(memb),q(memgmat),q(memnc))

        call dcopy_x(nek,q(memnc),1,xnew(1,1),1)
!  Temporary
        enew=e(npt)+edecr
  return
! print history
!
 n=3*ncen
 if(iprhis.gt.0.or.iprgrd.gt.0) then
   i1=npt
   if(iprhis.gt.0) i1=1
   do istep=i1,npt
     write(iout,"(/' Step=',i3,'  Energy=',f15.8/' CENTER',12x,'X',14x,'Y',14x,'Z',16x,'GX',13x,'GY',13x,'GZ')") istep,e(istep)
     do icen=1,ncen
       write(iout,'(1x,a,t10,3f15.8,3x,3f15.8)') atname(iuniq(icen)),(x(k,icen,istep),k=1,3),(g(k,icen,istep),k=1,3)
     end do
   end do
 end if

 if(have_hess) then
!
! Newton-Raphson using hessian
!
   iv=icorr(n*n)                                 !allocate memory for Hessian eigenvectors
   id=icorr(n)                                   !allocate memory for Hessian eigenvalues
   ig=icorr(n)                                   !allocate memory for transformed gradient
   if(iprhes.gt.0) call outsqr(h,n,n,n,'H')      !print hessian
   call fmove(h,q(iv),n*n)                       !copy h to q(iv)
   call diag2(n,n,q(id),q(iv))                   !diagonalize h, eigenvectors on q(iv)
   call mxva(q(iv),n,1,g(1,1,npt),1,q(ig),1,n,n) !transform gradient
   do i=0,n-1
     if(abs(q(id+i)).gt.thr) then
       q(ig+i)=-q(ig+i)/q(id+i)                  !update
     else
       q(ig+i)=0d0
     end if
   end do
   call mxva(q(iv),1,n,q(ig),1,xnew,1,n,n)       !transform update back
   call daxpy_x(n,1.d0,x(1,1,npt),1,xnew,1)      !add update to previous geometry
   call corlsr(iv)                               !release memory

 else
!
! simple steepest descent update, assuming that Hessian is unit matrix
!
   do icen=1,ncen
     do k=1,3
       xnew(k,icen)=x(k,icen,npt)-g(k,icen,npt)
     end do
   end do
 end if
!
! examples for some optimization options. Default values see lib/optg.registry. Also see manual for more options
! New options must be defined in lib/optg.registry.
!
! floating point parameters:
!
 call get_inpf('STEPMAX','OPTG',stepmx)  !max step length in optimization
!
! integer parameters:
!
 call get_inpi('MAXUPD','OPTG',maxupd)   !max number of hessian updates
!
! logical parameters:
!
  call get_inpl('BMAT','OPTG',use_bmat)  !if true, use b-matrix coordinates
!
! string parameters:
!
 call get_inps('METHOD','OPTG',method)   !optimization method (first two letters should be PP to get here).
!
! blas and lapack routines: please append _x to the standard name, e.g.  ! dscal_x, daxpy_x, dgemm_x etc
! (there are interfaces that convert the Molpro integer types to those used in the libraries)

 end subroutine optg_pp

 end module optpp


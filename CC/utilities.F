c This file contains simple subroutines used mainly for testing purposes.
c The subroutines are usually small and self-explaining.
      subroutine RawMatrixWriter(xmatrix,idimx,idimy,ndisk,reset)
      implicit none
      integer idimx,idimy,ndisk
      real*8 xmatrix(idimx,idimy)
      logical reset
      integer irec
      data    irec/0/
      save    irec
      if (reset) then 
      irec=0
      return
      endif
      irec=irec+1
      write(ndisk,rec=irec) xmatrix
      end
c=======================================================================
      subroutine expander(small,ismall,big,ibig)
      implicit none
      integer ismall,ibig
      real*8 small(ismall,ismall),big(ibig,ibig)
c
      integer idif,i,j,ib,jb
c
      idif=ibig-ismall
      if (idif.lt.0) STOP 'Error subroutine expand'
      call zeroit(big,ibig*ibig)
      if (idif.ne.0) then
      do j=1,ismall
        jb=j+idif
        do i=1,ismall
          ib=i+idif
          big(ib,jb)=small(i,j)
        enddo
      enddo
      else
      do j=1,ismall
        do i=1,ismall
          big(i,j)=small(i,j)
        enddo
      enddo
      endif
      end
c=======================================================================
      subroutine spectrum(element,ispectr,istart,istop)
      implicit none
      real*8 element
      integer istart,istop,ispectr(istart:istop)
c
      real*8 xlog
      integer ii
c
      xlog=-log10(element)
      do ii=istart,istop-1
      if (xlog.lt.dble(ii)) then
          ispectr(ii)=ispectr(ii)+1
          return
        endif
      enddo
      ispectr(istop)=ispectr(istop)+1
      end
c=======================================================================
      subroutine spectrum1(element,ispectr,xistart,xistop,xinterv)
      implicit none
      real*8 element,xistart,xistop,xinterv
      integer ispectr(*)
c
      real*8 xlog
      integer ii
      real*8 xepsilon,xi
      parameter (xepsilon=1.d-10)
c
      if (xinterv.le.0.0d0) STOP 'Error spectrum1'
      xlog=-log10(dabs(element))
      xi=xistart
      ii=0
      do
        ii=ii+1
        if (xlog.lt.xi) then
          ispectr(ii)=ispectr(ii)+1
          return
        endif
        xi=xi+xinterv
        if (xi.gt.(xistop+xepsilon)) exit
      enddo
      ispectr(ii)=ispectr(ii)+1
      end
c=======================================================================
      subroutine get_ij(itab,i,j)
      implicit none
      integer itab(2),i,j
      i=itab(1)
      j=itab(2)
      end
      subroutine put_ij(itab,i,j)
      implicit none
      integer itab(2),i,j
      itab(1)=i
      itab(2)=j
      end
c=======================================================================
      subroutine extend_by_unit(coeff,ext,idimc,idime)
      implicit none
      integer idimc,idime
      real*8 coeff(idimc,idimc),ext(idime,idime)
c
      integer i,j,ie,je,iadd
c
      if (idimc.eq.idime) then
        call tfer(coeff,ext,idimc*idimc)
        return
      endif
      if (idimc.gt.idime) STOP 'Error, condition impossible'
      iadd=idime-idimc
      do i=1,idimc
        ie=i+iadd
        do j=1,idimc
          je=j+iadd
          ext(ie,je)=coeff(i,j)
        enddo
      enddo
      do i=1,idime
        do j=1,idime
        if (i.le.iadd.or.j.le.iadd) ext(i,j)=0.0d0
        enddo
      enddo
      do i=1,iadd
        ext(i,i)=1.0d0
      enddo
      end
c=======================================================================
      integer function int_array(array,indexx)
      implicit none
      integer array(*)   ! originally this is array from bl
      integer indexx
      int_array=array(indexx)
      end
c=======================================================================
      subroutine int_array_write(array,indexx,int_value)
      implicit none
      integer array(*)   ! originally this is array from bl
      integer indexx,int_value
      array(indexx)=int_value
      end

      integer function ncores(ncf,ew,core)
      implicit real*8 (a-h,o-z)
      dimension ew(ncf)
      do 700 icore=1,ncf
        if(ew(icore).gt.core) go to 710
 700  continue
 710  ncores=icore-1
c      ncores is the number of core orbitals
      end

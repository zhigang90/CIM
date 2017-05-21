      subroutine precal (ia,ie,ja,je,basdat,pi23,nij,sij,pij)
c  this subroutine calculates the overlap integral of the s factor
c  and the center of the product gaussian px,py,pz for the functions in
c  a contracted shell pair.
c  arguments
c  intent(in)
c ia,ie= beginning and ending  primitive function of the first contracted shell
c ja,je= ditto for the senond shell
c basdat(13,nsh) contains the info on the primitive shells, see comments in 
c    twoel and basin
c pi23 is a constant, (2/pi)**3
c  intent (out)
c  nij= the nymber of primitive shell pairs = (ie-ia+1)*(je-ja+1)
c  sij= the overlap matrix, stored as sij(lj,li) where lj=contraction length
c      for j, li=contraction length for j
c  pij=  center coordinates of the product gaussian, stored as pij(3,lj,li)
c
 
      implicit real*8 (a-h,o-z)
      dimension sij(*), pij(*), ax(3), bx(3), basdat(13,*)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      ij=0
      ij3=-3
	nij=0
      do 30 i=ia,ie
         aa=basdat(1,i)
         do 10 l=1,3
   10    ax(l)=basdat(10+l,i)
      do 30 j=ja,je
         ij=ij+1
	   nij=nij+1
         ij3=ij3+3
         bb=basdat(1,j)
         apb=aa+bb
         r=zero
         do 20 l=1,3
            bx(l)=basdat(10+l,j)
            r=r+(ax(l)-bx(l))**2
            pij(ij3+l)=(aa*ax(l)+bb*bx(l))/apb
   20    continue
         sij(ij)=sqrt(sqrt(aa*bb))**3*exp(-aa*bb*r/apb)*pi23
   30 continue
c
      end

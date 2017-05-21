      subroutine indxyz(indxx,indyy,indzz)
c  This subroutine calculates an indexing array for the Rys integration
c  It also determines the number of roots for Rys integartion
c  Arguments
c  INTENT(IN)
c  It appears that the only input parameters of this function are the Cartesian
c      types of the functions, ityp,jtyp,ktyp,ltyp in COMMON /two1/
c  INTENT(OUT)
c  indxx,indyy,indzz: The index of the x,y,z factors for each function in the shell
c    quartet. Necessary dimension for 4 G functions= 15**4=50625 (dimensioned in 
c    the calling program TWOEL
c
c  NOTE that the results are the same whenever the types are the same, and therefore 
c     this routine needs to be called only once for a block of integrals of the same
c     type; alternatively, one could check the types agains the previous values
c
c     ** this subroutine calculates the indices for the
c     ** factors in the gaussian quadrature. the factors
c     ** are ordered like a 4-dimensional array xfc(5,11,5,5)
c
c  ??? why does the third index go to 9 when the rest go only to 5???
c  the indices (used to simulate this matrix) correspond to 5,11,5,5,
c  NOT 5,9,5,5 as in the original comment (changed now)
c   Note that the indices are ordered backwards: L,K,J,I; 
c
c     ** and the order is xfc(l,k,j,i) where i,j,k,l are the
c     ** four functions. the subroutine takes a list of the
c     ** cartesian components.the ordering of the raw integrals
c     ** is like this:
c     ** for s type s
c     ** for p type x,y,z
c     ** for l type s,x,y,z
c     ** for d type xx,yy,zz,xy,xz,yz
c     ** for f type xxx,xxy,xxz,xyy,xyz,xzz,yyy,yyz,yzz,zzz
c     ** for g type xxxx,xxxy,xxxz,xxyy,xxyz,xxzz,xyyy,xyyz
c     ** xyzz,xzzz,yyyy,yyyz,yyzz,yzzz,zzzz
c
      implicit real*8 (a-h,o-z)
      dimension indxx(*),indyy(*),indzz(*)
      dimension itm275(5),itm55(5),itm5(11)
      dimension nbeg(6),nend(6),nxtyp(35),nytyp(35),nztyp(35),ndege(6)
      common /two1/ iityp,jjtyp,kktyp,lltyp,
     1 ityp,jtyp,ktyp,ltyp,ncount,ncnt,nord,maxl,iwy,ilen1,
     1 jlen1,klen1,len1,ia,ie,ja,je,ka,ke,la,le,lx,eps,eps1,ca3,
     2 ca2,ca1,ca0
      common /rys/ xrys,rysr(10),rysw(10),t,f0,f1,f2,f3,f4,f5,nroots,
     1             nroot2,ie1,ie2,ie3,ie4,ie12,ie34
      data itm275/0,275,550,825,1100/
      data itm55/0,55,110,165,220/
      data itm5/0,5,10,15,20,25,30,35,40,45,50/
      data nbeg/1,2,1,5,11,21/,nend/1,4,4,10,20,35/
      data nxtyp/1, 2,1,1, 3,1,1,2,2,1, 4,3,3,2,2,2,1,1,1,1,
     1           5,4,4,3,3,3,2,2,2,2,1,1,1,1,1/
      data nytyp/1, 1,2,1, 1,3,1,2,1,2, 1,2,1,3,2,1,4,3,2,1,
     1           1,2,1,3,2,1,4,3,2,1,5,4,3,2,1/
      data nztyp/1,1,1,2, 1,1,3,1,2,2, 1,1,2,1,2,3,1,2,3,4,
     1           1,1,2,1,2,3,1,2,3,4,1,2,3,4,5/
      data ndege/1,2,2,3,4,5/
c
      ib1=nbeg(ityp)
      ib2=nbeg(jtyp)
      ib3=nbeg(ktyp)
      ib4=nbeg(ltyp)
      ie1=nend(ityp)
      ie2=nend(jtyp)
      ie3=nend(ktyp)
      ie4=nend(ltyp)
c
c      go to (10,30,50), nord
c  eliminate the integral triplet storage
c
   10 ncount=0
      do 20 ityp1=ib1,ie1
         lax=nxtyp(ityp1)
         lay=nytyp(ityp1)
         laz=nztyp(ityp1)
         indx1=itm275(lax)
         indy1=itm275(lay)
         indz1=itm275(laz)
      do 20 ityp2=ib2,ie2
         lbx=nxtyp(ityp2)
         lby=nytyp(ityp2)
         lbz=nztyp(ityp2)
         indx2=indx1+itm55(lbx)
         indy2=indy1+itm55(lby)
         indz2=indz1+itm55(lbz)
      do 20 ityp3=ib3,ie3
         lcx=nxtyp(ityp3)
         lcy=nytyp(ityp3)
         lcz=nztyp(ityp3)
         indx3=indx2+itm5(lcx)
         indy3=indy2+itm5(lcy)
         indz3=indz2+itm5(lcz)
      do 20 ityp4=ib4,ie4
         ncount=ncount+1
         indxx(ncount)=indx3+nxtyp(ityp4)
         indyy(ncount)=indy3+nytyp(ityp4)
         indzz(ncount)=indz3+nztyp(ityp4)
   20 continue
c     write(32,6666) ib1,ie1,ib2,ie2,ib3,ie3,ib4,ie4,ncount
6666  format(' ib',10i3)
c  Now determine nroots and assoviated quantities
      ie1=ndege(ityp)
      ie2=ndege(jtyp)
      ie3=ndege(ktyp)
      ie4=ndege(ltyp)
      ie12=ie1+ie2-1
      ie34=ie3+ie4-1
      nroot2=(ie1+ie2+ie3+ie4-2)
      nroots=nroot2/2

      return
c

      end

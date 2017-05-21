c=======================================================================
c Oct/Nov 97 KW: The major change in the FIRST_DER routine. Derivatives
c are calculated now in two steps with (1) higher and (2) lower ang.mom.
c For lower ang.mom. contributions loops go over these functions which
c have non-zero pawer of (x-X) ,(y-Y) and (z-Z) . These funcions are found
c by calling FIND_NON0 routine.
c (The previous version of first_der is kept in this file as well under
c  first_der_O name because there is a detailed description of all terms).
c=======================================================================
c            +-------------------------------------+
c            +                                     +
c            + two-electron integral derivatives   +
c            +           for gradient              +
c            +                                     +
c            +         Krzysztof Wolinski          +
c            +             July , 1997             +
c            +-------------------------------------+
c=======================================================================
c===========GRADIENT INTEGRAL DERIVATIVES ROUTINES==================
c
c     It is called from Calcint2 when WHERE='forc'
c
      subroutine force_der(bl,nbls,lnijr,lnklr,npij,ngcd,idx1,ixab)
c

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c
c only for first & second derivatives (for use in amshift):
c
      common /memor4b/ider0,ider1,ider2
c
c dimensions for assembling :
      common /dimasse/ lqijr,lqklr,lqmxr,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
c----------------------------------------------------------
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=lqij
      if(lqkl.gt.lqij) lqmx=lqkl
c----------------------------------------------------------
      ndim=4     ! dimension for buf2(ndim,*) used in first_der
c----------------------------------------------------------
      call getmem(3*nbls,ixabq)
      call conv24r(nbls,npij,bl(idx1),bl(ixab),bl(ixabq))
c----------------------------------------------------------
c find these functions which have non-zero power of x,y,z on A and C

      ijdim=lnij-nfu(nqij)
      kldim=lnkl-nfu(nqkl)
c
      call getint(ijdim,ijvecx)
      call getint(ijdim,ijvecy)
      call getint(ijdim,ijvecz)
c
      call getint(kldim,klvecx)
      call getint(kldim,klvecy)
      call getint(kldim,klvecz)
c
      call find_non0(nqij,lnij,nqkl,lnkl,
     *               nijx,nijy,nijz,nklx,nkly,nklz,
     *               bl(ijvecx),bl(ijvecy),bl(ijvecz),
     *               bl(klvecx),bl(klvecy),bl(klvecz) )
c----------------------------------------------------------
c contracted integrals with higher and lower ang.mom. are in buf now,
c not in buf2. In buf2 derivatives will be placed after call to
c first_der
c
      ibeg =ibuf
      ider1=ibuf2
      incr9=9*ngcd*nbls*lnij*lnkl
      ider0=ider1+incr9
c
      call first_der(ngcd,nbls,bl(ibeg),ndim,
     *               lnijr,lnklr,lnij,lnkl,nqij,nqkl,
     *               bl(ider1),bl(ider0),bl(ixabq),
     *               nijx,nijy,nijz,nklx,nkly,nklz,
     *               bl(ijvecx),bl(ijvecy),bl(ijvecz),
     *               bl(klvecx),bl(klvecy),bl(klvecz) )
c
c-----------------------------------------------------------------
c gradient works ONLY for segmented s,p basis set (NO L-shells)
c     if(lshellt.eq.0) go to 100
c.............
c 100 continue
c-----------------------------------------------------------------
c
      call retmem(6)
      call retmem(1)
c
      end
c=================================================================
      subroutine conv24r(nbls,npij,idx1,xab,xabq)
      implicit real*8 (a-h,o-z)
c
      dimension idx1(nbls)
      dimension xab(npij,3),xabq(nbls,3)
c
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
c     klpar=idx2(ijkl)
        do 150 i=1,3
        xabq(ijkl,i)=xab(ijpar,i)
c       xcdq(ijkl,i)=xcd(klpar,i)
  150   continue
  100 continue
      end
c=================================================================
c this is ok but slow so I keep it with the name first_der_Old
c
      subroutine first_der_O(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
     *                     nqij,nqkl,deriv,der00,xab)
      implicit real*8 (a-h,o-z)
c
c this subroutine is also called when second-derivatives are calculated
c because first-derivatives are needed in the shifting procedure for
c second-der. That's why dimension for buf2(ndim,*,*,*,*) has ndim=4
c for first- and ndim=10 for second-derivatives.
c
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
cccc  dimension buf2(4,nbls,lnijr,lnklr,ngcd) OR buf2(10,etc.)
      dimension buf2(ndim,nbls,lnijr,lnklr,ngcd)
      dimension deriv(9,nbls,lnij,lnkl,ngcd),der00(nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c-------------------------------------------------------------------
c This routine calculates UNFINISHED gradient derivatives
c of the type
c
c d/dAx (a+b,s|c+d,s)=2a*(a+b+1x,s|c+d,s)-n_ab_x*(a+b-1x,s|c+d,s)
c
c d/dBx (a+b,s|c+d,s)=2b*(a+b,s+1x|c+d,s)-  0*(a+b,s-1x|c+d,s)
c                    =2b[(a+b+1x,s|c+d,s)+(Ax-Bx)*(a+b,s|c+d,s)].
c
c d/dCx (a+b,s|c+d,s)=2c*(a+b,s|c+d+1x,s)-n_cd_x*(a+b,s|c+d-1x,s)
c
c-------------------------------------------------------------------
c INPUT buf2(1,nbls,lnijr,lnklr,ngcd) - ordinary 2-el.integ.
c INPUT buf2(2,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_a
c INPUT buf2(3,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_b
c INPUT buf2(4,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_c
c
c OUTPUT :
c
c--->           derAX=deriv(1,nbls,lnij,lnkl,ngcd),
c--->           derBX=deriv(2,nbls,lnij,lnkl,ngcd)
c--->           derCX=deriv(3,nbls,lnij,lnkl,ngcd),
c
c--->           derAY=deriv(4,nbls,lnij,lnkl,ngcd),
c--->           derBY=deriv(5,nbls,lnij,lnkl,ngcd),
c--->           derCY=deriv(6,nbls,lnij,lnkl,ngcd),
c
c--->           derAZ=deriv(7,nbls,lnij,lnkl,ngcd),
c--->           derBZ=deriv(8,nbls,lnij,lnkl,ngcd),
c--->           derCZ=deriv(9,nbls,lnij,lnkl,ngcd),
c
c--->           der00(nbls,lnij,lnkl,ngcd) - ordinary integr.
c               needed in shifting for d/dA,d/dB (1->2) and d/dC (3->4)
c-------------------------------------------------------------------
c
      do 200 kl=nfu(nqkl)+1,lnkl
      n_cd_x=nia(1,kl)
      n_cd_y=nia(2,kl)
      n_cd_z=nia(3,kl)
      klpx=npxyz(1,kl)
      klpy=npxyz(2,kl)
      klpz=npxyz(3,kl)
c
      klmx=0
      klmy=0
      klmz=0
      if(n_cd_x.gt.0) klmx=nmxyz(1,kl)
      if(n_cd_y.gt.0) klmy=nmxyz(2,kl)
      if(n_cd_z.gt.0) klmz=nmxyz(3,kl)
c
      do 200 ij=nfu(nqij)+1,lnij
      n_ab_x=nia(1,ij)
      n_ab_y=nia(2,ij)
      n_ab_z=nia(3,ij)
      ijpx=npxyz(1,ij)
      ijpy=npxyz(2,ij)
      ijpz=npxyz(3,ij)
      ijmx=0
      ijmy=0
      ijmz=0
      if(n_ab_x.gt.0) ijmx=nmxyz(1,ij)
      if(n_ab_y.gt.0) ijmy=nmxyz(2,ij)
      if(n_ab_z.gt.0) ijmz=nmxyz(3,ij)
c
      do 225 iqu=1,ngcd
        do 250 ijkl=1,nbls
c---------------
        two_a_x=buf2(2,ijkl,ijpx,kl,iqu)
        two_a_y=buf2(2,ijkl,ijpy,kl,iqu)
        two_a_z=buf2(2,ijkl,ijpz,kl,iqu)
c
        two_b_0=buf2(3,ijkl,ij,kl,iqu)
        two_b_x=buf2(3,ijkl,ijpx,kl,iqu)
        two_b_y=buf2(3,ijkl,ijpy,kl,iqu)
        two_b_z=buf2(3,ijkl,ijpz,kl,iqu)
c
        two_c_x=buf2(4,ijkl,ij,klpx,iqu)
        two_c_y=buf2(4,ijkl,ij,klpy,iqu)
        two_c_z=buf2(4,ijkl,ij,klpz,iqu)
c
c---------------
        der00(ijkl,ij,kl,iqu)=buf2(1,ijkl,ij,kl,iqu)
c---------------
c--X deriv.
c
        if(n_ab_x.gt.0) then
           x_n_ab=n_ab_x*buf2(1,ijkl,ijmx,kl,iqu)
           deriv(1,ijkl,ij,kl,iqu)=two_a_x - x_n_ab
        else
           deriv(1,ijkl,ij,kl,iqu)=two_a_x
        endif
        deriv(2,ijkl,ij,kl,iqu)=two_b_x + xab(ijkl,1)*two_b_0
        if(n_cd_x.gt.0) then
           x_n_cd=n_cd_x*buf2(1,ijkl,ij,klmx,iqu)
           deriv(3,ijkl,ij,kl,iqu)=two_c_x - x_n_cd
        else
           deriv(3,ijkl,ij,kl,iqu)=two_c_x
        endif
c---------------
c--Y deriv.
c
        if(n_ab_y.gt.0) then
           y_n_ab=n_ab_y*buf2(1,ijkl,ijmy,kl,iqu)
           deriv(4,ijkl,ij,kl,iqu)=two_a_y - y_n_ab
        else
           deriv(4,ijkl,ij,kl,iqu)=two_a_y
        endif
        deriv(5,ijkl,ij,kl,iqu)=two_b_y + xab(ijkl,2)*two_b_0
        if(n_cd_y.gt.0) then
           y_n_cd=n_cd_y*buf2(1,ijkl,ij,klmy,iqu)
           deriv(6,ijkl,ij,kl,iqu)=two_c_y - y_n_cd
        else
           deriv(6,ijkl,ij,kl,iqu)=two_c_y
        endif
c---------------
c--Z deriv.
c
        if(n_ab_z.gt.0) then
           z_n_ab=n_ab_z*buf2(1,ijkl,ijmz,kl,iqu)
           deriv(7,ijkl,ij,kl,iqu)=two_a_z - z_n_ab
        else
           deriv(7,ijkl,ij,kl,iqu)=two_a_z
        endif
        deriv(8,ijkl,ij,kl,iqu)=two_b_z + xab(ijkl,3)*two_b_0
        if(n_cd_z.gt.0) then
           z_n_cd=n_cd_z*buf2(1,ijkl,ij,klmz,iqu)
           deriv(9,ijkl,ij,kl,iqu)=two_c_z - z_n_cd
        else
           deriv(9,ijkl,ij,kl,iqu)=two_c_z
        endif
c---------------
c
  250   continue
  225 continue
  200 continue
c
      end
c=================================================================
      subroutine first_der(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *                     nqij,nqkl,deriv,der00,xab,
     *                     nijx,nijy,nijz,nklx,nkly,nklz,
     *                     ijvecx,ijvecy,ijvecz,klvecx,klvecy,klvecz )
      implicit real*8 (a-h,o-z)
c
c-------------------------------------------------------------------
c
c This routine calculates UNFINISHED gradient derivatives of the type
c
c d/dAx (a+b,s|c+d,s)=2a*(a+b+1x,s|c+d,s)-n_ab_x*(a+b-1x,s|c+d,s)
c
c d/dBx (a+b,s|c+d,s)=2b*(a+b,s+1x|c+d,s)-  0*(a+b,s-1x|c+d,s)
c                    =2b[(a+b+1x,s|c+d,s)+(Ax-Bx)*(a+b,s|c+d,s)].
c
c d/dCx (a+b,s|c+d,s)=2c*(a+b,s|c+d+1x,s)-n_cd_x*(a+b,s|c+d-1x,s)
c
c-------------------------------------------------------------------
c INPUT buf2(1,nbls,lnijr,lnklr,ngcd) - ordinary 2-el.integ.
c INPUT buf2(2,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_a
c INPUT buf2(3,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_b
c INPUT buf2(4,nbls,lnijr,lnklr,ngcd) - 2-el.integ. resc by 2*exp_c
c
c OUTPUT :
c
c               derAX=deriv(1,nbls,lnij,lnkl,ngcd),
c               derBX=deriv(2,nbls,lnij,lnkl,ngcd)
c               derCX=deriv(3,nbls,lnij,lnkl,ngcd),
c
c               derAY=deriv(4,nbls,lnij,lnkl,ngcd),
c               derBY=deriv(5,nbls,lnij,lnkl,ngcd),
c               derCY=deriv(6,nbls,lnij,lnkl,ngcd),
c
c               derAZ=deriv(7,nbls,lnij,lnkl,ngcd),
c               derBZ=deriv(8,nbls,lnij,lnkl,ngcd),
c               derCZ=deriv(9,nbls,lnij,lnkl,ngcd),
c
c               der00(nbls,lnij,lnkl,ngcd) - ordinary integr.
c               needed in shifting for d/dA,d/dB (1->2) and d/dC (3->4)
c-------------------------------------------------------------------
c order of deriv : ax,bx,cx,ay,by,cy,az,bz,cz
c                   1  2  3  4  5  6  7  8  9
c-------------------------------------------------------------------
c this routine is also called when second-derivatives are calculated
c because first-derivatives are needed in the shifting procedure for
c second-der. That's why dimension for buf2(ndim,*,*,*,*) has ndim=4
c for first- and ndim=10 for second-derivatives.
c-------------------------------------------------------------------
c
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
cccc  dimension buf2(4,nbls,lnijr,lnklr,ngcd) OR buf2(10,etc.)
c200  dimension buf2(ndim,nbls,lnijr,lnklr,ngcd)
      dimension buf2(nbls,lnijr,lnklr,ngcd, ndim)
      dimension deriv(9,nbls,lnij,lnkl,ngcd),der00(nbls,lnij,lnkl,ngcd)
      dimension xab(nbls,3)
c
      dimension ijvecx(*),ijvecy(*),ijvecz(*)
      dimension klvecx(*),klvecy(*),klvecz(*)
c-------------------------------------------------------------------
c
      do 200 iqu=1,ngcd
         do 210 kl=nfu(nqkl)+1,lnkl
         klpx=npxyz(1,kl)
         klpy=npxyz(2,kl)
         klpz=npxyz(3,kl)
c
c First calculate contributions with higher (by 1) angular momentum
c
            do 220 ij=nfu(nqij)+1,lnij
            ijpx=npxyz(1,ij)
            ijpy=npxyz(2,ij)
            ijpz=npxyz(3,ij)
               do 230 ijkl=1,nbls
c--------------------------------------------------------------
               two_b_0=buf2(ijkl,ij,kl,iqu,3)
c--------------------------------------------------------------
               der00(ijkl,ij,kl,iqu)=buf2(ijkl,ij,kl,iqu,1)
c
               deriv(1,ijkl,ij,kl,iqu)=buf2(ijkl,ijpx,kl,iqu,2)
               deriv(2,ijkl,ij,kl,iqu)=buf2(ijkl,ijpx,kl,iqu,3)
     *                                    + xab(ijkl,1)*two_b_0
               deriv(3,ijkl,ij,kl,iqu)=buf2(ijkl,ij,klpx,iqu,4)
               deriv(4,ijkl,ij,kl,iqu)=buf2(ijkl,ijpy,kl,iqu,2)
               deriv(5,ijkl,ij,kl,iqu)=buf2(ijkl,ijpy,kl,iqu,3)
     *                                    + xab(ijkl,2)*two_b_0
               deriv(6,ijkl,ij,kl,iqu)=buf2(ijkl,ij,klpy,iqu,4)
               deriv(7,ijkl,ij,kl,iqu)=buf2(ijkl,ijpz,kl,iqu,2)
               deriv(8,ijkl,ij,kl,iqu)=buf2(ijkl,ijpz,kl,iqu,3)
     *                                    + xab(ijkl,3)*two_b_0
               deriv(9,ijkl,ij,kl,iqu)=buf2(ijkl,ij,klpz,iqu,4)
c--------------------------------------------------------------
  230          continue
  220       continue
c--------------------------------------------------------------
c Now calculate contributions with lower (by 1) angular momentum
c
c.......... X deriv. for A center............................
c
            do 221 ijx=1,nijx
            ij=ijvecx(ijx)  ! only these which have n_ab_x >0
c           n_ab_x=nia(1,ij)
            x_n_ab=dble( nia(1,ij) )
            ijmx=nmxyz(1,ij)
               do 231 ijkl=1,nbls
                  deriv(1,ijkl,ij,kl,iqu)=deriv(1,ijkl,ij,kl,iqu)
     *                          - x_n_ab*buf2(ijkl,ijmx,kl,iqu,1)
  231          continue
  221       continue
c
c.......... Y deriv. for A center............................
c
            do 222 ijy=1,nijy
            ij=ijvecy(ijy)  ! only these which have n_ab_y >0
cc          n_ab_y=nia(2,ij)
            y_n_ab=dble( nia(2,ij) )
            ijmy=nmxyz(2,ij)
               do 232 ijkl=1,nbls
                  deriv(4,ijkl,ij,kl,iqu)=deriv(4,ijkl,ij,kl,iqu)
     *                          - y_n_ab*buf2(ijkl,ijmy,kl,iqu,1)
  232          continue
  222       continue
c
c.......... Z deriv. for A center............................
c
            do 223 ijz=1,nijz
            ij=ijvecz(ijz)  ! only these which have n_ab_z >0
c           n_ab_z=nia(3,ij)
            z_n_ab=dble( nia(3,ij) )
            ijmz=nmxyz(3,ij)
               do 233 ijkl=1,nbls
                  deriv(7,ijkl,ij,kl,iqu)=deriv(7,ijkl,ij,kl,iqu)
     *                          - z_n_ab*buf2(ijkl,ijmz,kl,iqu,1)
  233          continue
  223       continue
c
  210    continue     ! over kl
  200 continue        ! over general contraction
c--------------------------------------------------------------
c
c Now calculate contributions with lower (by 1) angular momentum
c on centers CD i.e. with klm>0 .
c
c
      ixyz=1
      call lower_cd(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *              nqij,nklx,klvecx,ixyz,deriv)
      ixyz=2
      call lower_cd(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *              nqij,nkly,klvecy,ixyz,deriv)
      ixyz=3
      call lower_cd(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *              nqij,nklz,klvecz,ixyz,deriv)
c
      end
c===================================================================
      subroutine find_non0(nqij,lnij,nqkl,lnkl,
     *                     nijx,nijy,nijz,nklx,nkly,nklz,
     *                     ijvecx,ijvecy,ijvecz,klvecx,klvecy,klvecz)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
      common /logic11/ npxyz(3,1)
c
      dimension ijvecx(*),ijvecy(*),ijvecz(*)
      dimension klvecx(*),klvecy(*),klvecz(*)
c
      nijx=0
      nijy=0
      nijz=0
      do 100 ij=nfu(nqij)+1,lnij
         n_ab_x=nia(1,ij)
         n_ab_y=nia(2,ij)
         n_ab_z=nia(3,ij)
         if(n_ab_x.gt.0) then
            nijx=nijx+1
            ijvecx(nijx)=ij
         endif
         if(n_ab_y.gt.0) then
            nijy=nijy+1
            ijvecy(nijy)=ij
         endif
         if(n_ab_z.gt.0) then
            nijz=nijz+1
            ijvecz(nijz)=ij
         endif
  100 continue
c
      nklx=0
      nkly=0
      nklz=0
      do 200 kl=nfu(nqkl)+1,lnkl
         n_cd_x=nia(1,kl)
         n_cd_y=nia(2,kl)
         n_cd_z=nia(3,kl)
         if(n_cd_x.gt.0) then
            nklx=nklx+1
            klvecx(nklx)=kl
         endif
         if(n_cd_y.gt.0) then
            nkly=nkly+1
            klvecy(nkly)=kl
         endif
         if(n_cd_z.gt.0) then
            nklz=nklz+1
            klvecz(nklz)=kl
         endif
  200 continue
c
      end
c===================================================================
      subroutine lower_cd(ngcd,nbls,buf2,ndim,lnijr,lnklr,lnij,lnkl,
     *                    nqij,nklxyz,klvecxyz,ixyz, deriv)
      implicit real*8 (a-h,o-z)
      common /logic4/ nfu(1)
      common /logic9/ nia(3,1)
      common /logic10/ nmxyz(3,1)
c
      dimension buf2(nbls,lnijr,lnklr,ngcd,ndim)
      dimension klvecxyz(*)
      dimension deriv(9,nbls,lnij,lnkl,ngcd)
c
c-------------------------------------------------------------------
c      X , Y or Z deriv. for C center
c-------------------------------------------------------------------
      idxyz=3*ixyz  ! shows what derivative (3 for x,6 for y & 9 forz)
c
      do 200 iqu=1,ngcd
         do 210 klxyz=1,nklxyz
         kl=klvecxyz(klxyz)     ! only these with n_cd_x or_y or_z >0
c        n_cd_xyz=nia(ixyz,kl)
         xyz_n_cd= dble( nia(ixyz,kl) )
         klmxyz=nmxyz(ixyz,kl)
            do 220 ij=nfu(nqij)+1,lnij
               do 230 ijkl=1,nbls
               deriv(idxyz,ijkl,ij,kl,iqu)=deriv(idxyz,ijkl,ij,kl,iqu)
     *                           - xyz_n_cd*buf2(ijkl,ij,klmxyz,iqu,1)
  230          continue
  220       continue
  210    continue
  200 continue
c
      end
c===================================================================
      subroutine onecentr(iis,jjs,inx,npij,npkl,npklx,
     *           ibl,kbl,ijbl,nbl2,nijbeg,nijend,nklbeg,nklend,
     *           ipres,ijcent,klcent)
c
      dimension inx(12,*),iis(*),jjs(*),ijbl(nbl2,*)
      dimension ijcent(npij),klcent(npkl)
      dimension ipres(*)
c
c precalculations for the pairs IJ :
c
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl(ibl,ijp)
        ics=iis(ijcs)
        jcs=jjs(ijcs)
        iatom=inx(2,ics)
        jatom=inx(2,jcs)
        ijcent(ijpar)=-12
        if(iatom.eq.jatom) ijcent(ijpar)=iatom
  100 continue
c
c precalculations for the pairs KL :
c
      IF(npklx.ne.0) then
        klpar=0
        do 200 klp=nklbeg,nklend
        klpar=klpar+1
          klcs=ijbl(kbl,klp)
          kcs=iis(klcs)
          lcs=jjs(klcs)
          katom=inx(2,kcs)
          latom=inx(2,lcs)
          klcent(klpar)=-34
          if(katom.eq.latom) klcent(klpar)=katom
  200   continue
      ELSE              !  for a diagonal case :
c                          pairs IJ and KL are the same
        do 300 ijpar=1,npij
           klc=ijcent(ijpar)
           if(klc.eq.-12) klc=-34
           klcent(ijpar)=klc
  300   continue
c
      ENDIF
c
c eliminate contracted quartets of the type (AA|AA) (centers)
c
      ijkl=0
      do 400 ijpar=1,npij
      ijc=ijcent(ijpar)
      npklend=npkl
      if(npklx.eq.0) npklend=ijpar
      do 400 klpar=1,npklend
      klc=klcent(klpar)
      ijkl=ijkl+1
      if(ijc.eq.klc) ipres(ijkl)=0
  400 continue
c
      end
c===================================================================
      subroutine onecentx(inx_1,inx_2,npij,inx_3,inx_4,npkl,npklx,
     *           ibl,kbl,
     *           ijbl_12,nbl2_ijd,
     *           ijbl_34,nbl2_kld,
     *           nijbeg,nijend,nklbeg,nklend,
     *           ipres,ijcent,klcent)
c
      dimension inx_1(12,*),inx_2(12,*),inx_3(12,*),inx_4(12,*)
      dimension ijbl_12(nbl2_ijd,*)
      dimension ijbl_34(nbl2_kld,*)
      dimension ijcent(npij),klcent(npkl)
      dimension ipres(*)
c
c precalculations for the pairs IJ :
c
      ijpar=0
      do 100 ijp=nijbeg,nijend
      ijpar=ijpar+1
        ijcs=ijbl_12(ibl,ijp)
        call get_ij_half(ijcs,ics,jcs)
        iatom=inx_1(2,ics)
        jatom=inx_2(2,jcs)
        ijcent(ijpar)=-12
        if(iatom.eq.jatom) ijcent(ijpar)=iatom
ccc   write(6,*)'ijcs=',ijcs,' ics,jcs=',ics,jcs,' iat,jat=',iatom,jatom
  100 continue
c
c precalculations for the pairs KL :
c
      IF(npklx.ne.0) then
        klpar=0
        do 200 klp=nklbeg,nklend
        klpar=klpar+1
          klcs=ijbl_34(kbl,klp)
          call get_ij_half(klcs,kcs,lcs)
          katom=inx_3(2,kcs)
          latom=inx_4(2,lcs)
          klcent(klpar)=-34
          if(katom.eq.latom) klcent(klpar)=katom
ccc   write(6,*)'klcs=',klcs,' kcs,lcs=',kcs,lcs,' kat,lat=',katom,latom
  200   continue
      ELSE              !  for a diagonal case :
c                          pairs IJ and KL are the same
        do 300 ijpar=1,npij
           klc=ijcent(ijpar)
           if(klc.eq.-12) klc=-34
           klcent(ijpar)=klc
  300   continue
c
      ENDIF
c
c eliminate contracted quartets of the type (AA|AA) (centers)
c
      ijkl=0
      do 400 ijpar=1,npij
      ijc=ijcent(ijpar)
      npklend=npkl
      if(npklx.eq.0) npklend=ijpar
      do 400 klpar=1,npklend
      klc=klcent(klpar)
      ijkl=ijkl+1
      if(ijc.eq.klc) ipres(ijkl)=0
ccc   write(6,*)'ijkl=',ijkl,' ijp,klp=',ijpar,klpar,' ijc,klc=',ijc,klc
  400 continue
c
      end

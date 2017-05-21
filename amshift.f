c tfer subroutine calls were changed to BLAS dcopy calls to eliminate
c overhead
c loops were eliminated where subsequent areas were copied
c still these memory copies take up almost 10% of all 2-el time
c some algorithm change would be welcome here
c
c-----------------------------------------------------------------------
c Nov 97 KW : Coping data in the SHIFT0L and SHIFT0l_DER1 has been
c             reduced by factor of 2. Some data are now directly used
c             in the HORIZ12_AB (SCF+NMRintegrals) and HORIZ12_DER1
c             (gradient integrals).
c
c
c NOTE : similar changes should be made in SHIFTnL routines which handle
c        L-shell cases for SCF & NMR integrals as well as
c        in the SHIFT0l_DER2 for second derivatives.
c-----------------------------------------------------------------------
C
C  THESE ROUTINES SHIFT THE ANGULAR MOMENTUM
C
C         FROM POSITION 1 TO POSITION 2
C
C                   AND
C
C         FROM POSITION 3 TO POSITION 4
c-----------------------------------------------------------------------
c    for re-ordered basis set
c   nqi.ge.nqj  and  nqk.ge.nql
c
c   other cases are not included here !
c-----------------------------------------------------------------------
      subroutine amshift(nbls,l01,l02,npij,npkl,ngcd)

      use memory

      implicit real*8 (a-h,o-z)
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      common /cpu/ intsize,iacc,icache,memreal
cxx
      common /logic4/ nfu(1)
cxx
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
c     common /big/ bl(1)
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
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
C**********************************
c* dimensions for "shifts"
c
      lsmx=max(lnij,lnkl)
      lsjl=max(nfu(nqj+1),nfu(nql+1))
c
      lqij=nfu(nqij+1)
      lqkl=nfu(nqkl+1)
      lqmx=max(lqij,lqkl)
c
c* memory for array INDXX :
c
      mindxx=lnijkl
c     if(intsize.ne.1) mindxx=lnijkl/intsize+1
c
      call getint(mindxx,indxx)
      call make_indxx(bl(indxx),lni,lnj,lnk,lnl)
c------------------------------------------------
c for odinary integrals :
c
          nbuf2=ibuf2
          nbuf=ibuf
          nbuh=ibuf      ! only for hessian
          m=1
c
c for giao integral derivatives :
          if(where.eq.'shif') then
            nmr =ngcd*nbls*lnijkl
            nbuf=ibuf+nmr
            m=6
          endif
c
c for gradient integral derivatives :
          if(where.eq.'forc') then
            m=9
          endif
c
c for hessian  integral derivatives :
          if(where.eq.'hess') then
            nbuh=ibuf + 9*nbls*lnijkl*ngcd
            m=45
          endif
c-
          mnbls=m*nbls
c-
c------------------------------------------------
C   SPECIAL CASES WHERE THE SHIFTING OF ANGULAR
C   MOMENTUM IS NOT NEEDED AT ALL / (DS|SS)..(SS|SD),
C   (XS|YS),(XS|SY),(SX|YS),(SX|SY) /
C
      IF(NQIJ.EQ.NSIJ .AND. NQKL.EQ.NSKL) THEN
c
c         find apropriate matrices with l-shells
c
          ibfijx=ibfij1
          if(lshelij.eq.2) ibfijx=ibfij2
          ibfklx=ibfkl1
          if(lshelkl.eq.2) ibfklx=ibfkl2
          ibf2lx=ibf2l1
          if(lcas2(2).eq.1) ibf2lx=ibf2l2
          if(lcas2(3).eq.1) ibf2lx=ibf2l3
          if(lcas2(4).eq.1) ibf2lx=ibf2l4
c
          if(where.eq.'hess') then
c           input : unfinished integ.derivatves:
c           ider1 - pointer to the 1st-der unfinished integ.
c           ider2 - pointer to the 2ed-der unfinished integ.
c
c           output : ibuf for 1st-deriv
c                    ibuh for 1st-deriv
c                    ibuh=ibuf + 9*nbls*lnijkl*ngcd
c
             incrf=9*nbls*lnijkl
             incrh=45*nbls*lnijkl
c
             incre1=9*nbls*l01*l02   ! increment for input (unfinished)
c                                      1st-deriv integrals
             incre2=45*nbls*l01*l02  ! increment for input (unfinished)
c                                      2ed-deriv integrals
c2002
c            do iqu=1,ngcd
c               jbuf =ibuf +(iqu-1)*incrf     ! output 1st-der.   integ.
c               jbuf2=ider1+(iqu-1)*incre1    ! input 1st-der.unf.integ.
c               call noshift(l01,l02,9*nbls,  bl(jbuf), bl(jbuf2),
c    *          bl(ibfijx),           bl(ibfklx),
c    *          bl(ibf2lx),
c    *          lqij,lqkl,
c    *          bl(indxx),lni*lnj,lnk*lnl)
c            enddo
c2002
             ibuh=ibuf + 9*nbls*lnijkl*ngcd
             do iqu=1,ngcd
                jbuf =ibuh +(iqu-1)*incrh     ! output 2ed-der.   integ.
                jbuf2=ider2+(iqu-1)*incre2    ! input 2ed-der.unf.integ.
                call noshift(l01,l02,45*nbls,  bl(jbuf), bl(jbuf2),
     *          bl(ibfijx),           bl(ibfklx),
     *          bl(ibf2lx),
     *          lqij,lqkl,
     *          bl(indxx),lni*lnj,lnk*lnl)
             enddo
          else           !    everything else but hess :
             incre =mnbls*lnijkl
             incre2=mnbls*l01*l02
             do iqu=1,ngcd
                jbuf =nbuf +(iqu-1)*incre
                jbuf2=nbuf2+(iqu-1)*incre2
                call noshift(l01,l02,mnbls,  bl(jbuf), bl(jbuf2),
     *          bl(ibfijx),           bl(ibfklx),
     *          bl(ibf2lx),
     *          lqij,lqkl,
     *          bl(indxx),lni*lnj,lnk*lnl)
             enddo
          endif
          call retmem(1)
          return
      ENDIF
c------------------------------------------------
c
cccc  if(where.eq.'fock' .or. where.eq.'shif') then
      if(where.ne.'forc') then
         call convr3(bl,m,nbls,npij,npkl,bl(idx1),bl(idx2),
     *               bl(ixab),bl(ixcd),ixabn,ixcdn)
      else
         call getint(3*nbls,ixab_e0_vec)
         call getint(3*nbls,ixcd_e0_vec)
         call getint(3*nbls,ixab_n0_vec)
         call getint(3*nbls,ixcd_n0_vec)
         call getint(3,ixab_e0)
         call getint(3,ixcd_e0)
         call getint(3,ixab_n0)
         call getint(3,ixcd_n0)
c
         call convr3_forc(bl,m,nbls,npij,npkl,bl(idx1),bl(idx2),
     *               bl(ixab),bl(ixcd),ixabn,ixcdn,
     *               bl(ixab_e0),bl(ixab_n0),bl(ixcd_e0),bl(ixcd_n0),
     *               bl(ixab_e0_vec),bl(ixab_n0_vec),
     *               bl(ixcd_e0_vec),bl(ixcd_n0_vec) )
      endif
c
c
       if (lshellt.eq.0) then
c-
        incre =mnbls*lnijkl
        incre2=mnbls*l01*l02
c-
         do 200 iqu=1,ngcd
           jbuf =nbuf +(iqu-1)*incre      ! returning integrals
           jbuf2=nbuf2+(iqu-1)*incre2
c
           if(where.eq.'fock' .or. where.eq.'shif') then
              call shift0l(bl(jbuf),bl(jbuf2),
     *                     l01,l02,bl(iwij),lsmx,lsjl,
     *                     bl(ixij),nfu(nqi+1),lqij,lnkl,mnbls,
     *                     bl(ixabn),bl(ixcdn),
     *                     bl(indxx),lni,lnj,lnk,lnl)
           endif
           if(where.eq.'forc') then
              jbuf0=ider0+(iqu-1)*nbls*l01*l02
              iwij0=iwij +mnbls*lsmx*lsjl
              ixij0=ixij +mnbls*nfu(nqi+1)*lqij*lnkl
              call shift0l_der1(bl(jbuf),bl(jbuf2),bl(jbuf0),
     *                         l01,l02,bl(iwij),bl(iwij0),lsmx,lsjl,
     *                         bl(ixij),bl(ixij0),nfu(nqi+1),lqij,lnkl,
     *                         mnbls,nbls,
     *                         bl(ixabn),bl(ixcdn),
     *                         bl(indxx),lni,lnj,lnk,lnl,
     *               bl(ixab_e0),bl(ixab_n0),bl(ixcd_e0),bl(ixcd_n0),
     *               bl(ixab_e0_vec),bl(ixab_n0_vec),
     *               bl(ixcd_e0_vec),bl(ixcd_n0_vec) )
           endif
           if(where.eq.'hess') then
c
              kbuf =nbuf +(iqu-1)*9*nbls*lnijkl  ! output gradient integ
              kbug =nbuh +(iqu-1)*45*nbls*lnijkl ! output hessian  integ
c
              jbuf0=ider0+(iqu-1)*nbls*l01*l02    ! input 0th-der integ.
              jbuf1=ider1+(iqu-1)*9*nbls*l01*l02  ! input 1st-der integ.
              jbuf2=ider2+(iqu-1)*45*nbls*l01*l02 ! input 2ed-der integ.
c
              iwij1=iwij +mnbls*lsmx*lsjl
              ixij1=ixij +mnbls*nfu(nqi+1)*lqij*lnkl
              iwij0=iwij1 +    9*nbls*lsmx*lsjl
              ixij0=ixij1 +    9*nbls*nfu(nqi+1)*lqij*lnkl
c
c only 2ed    call shift0l_der2(bl(jbuf),bl(jbuf2),bl(jbuf1),bl(jbuf0),
c
c                                  der1     der2     on return
              call shift0l_der2(bl(kbuf),bl(kbug),
     *                         bl(jbuf2),bl(jbuf1),bl(jbuf0),
     *                         l01,l02,bl(iwij),bl(iwij1),bl(iwij0),
     *                         lsmx,lsjl,
     *                         bl(ixij),bl(ixij1),bl(ixij0),
     *                         nfu(nqi+1),lqij,lnkl, mnbls,nbls,
     *                         bl(ixabn),bl(ixcdn),
     *                         bl(indxx),lni,lnj,lnk,lnl)
           endif
  200    continue
         call retmem(3)
         if(where.eq.'forc') call retmem(8)
         return
       endif
c
c- 1 l-shell
       if (lshellt.eq.1) then
c---
          jbuf  = nbuf
          jbuf2 =nbuf2
          jbfijx=ibfij1
          if(lshelij.eq.2) jbfijx=ibfij2
          jbfklx=ibfkl1
          if(lshelkl.eq.2) jbfklx=ibfkl2
c---
          call shift1l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfijx),bl(jbfklx),
     *                 lqij,lqkl,
     *                 bl(iwij),lsmx, bl(ivij),lsjl,
     *                 bl(ixij),nfu(nqi+1),lnkl,bl(iyij),nfu(nqj+1),
     *                 mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c
          call retmem(3)
          return
       endif
ccc
c- 2 l-shell
       if (lshellt.eq.2) then
c-
          jbuf  = nbuf
          jbuf2 =nbuf2
          jbfij1=ibfij1
          jbfij2=ibfij2
          jbfkl1=ibfkl1
          jbfkl2=ibfkl2
c-
          jbfij3=ibfij3
          jbfkl3=ibfkl3
          jbf2l12=ibf2l1
          if(lcas2(2).eq.1) jbf2l12=ibf2l2
          jbf2l34=ibf2l3
          if(lcas2(4).eq.1) jbf2l34=ibf2l4
c---
          call shift2l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l12),bl(jbf2l34),
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *          bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
ccc
c- 3 l-shell
c
       if (lshellt.eq.3) then
c---
          jbuf  = nbuf
          jbuf2 =nbuf2
          jbfij1=ibfij1
          jbfij2=ibfij2
          jbfkl1=ibfkl1
          jbfkl2=ibfkl2
c-
          jbfij3=ibfij3
          jbfkl3=ibfkl3
          jbf2l1=ibf2l1
          jbf2l2=ibf2l2
          jbf2l3=ibf2l3
          jbf2l4=ibf2l4
c-
          jbf3l12=ibf3l1
          ix3l12=ix3l1
          if(lcas3(2).eq.1) then
              jbf3l12=ibf3l2
              ix3l12=ix3l2
          endif
          jbf3l34=ibf3l3
          ix3l34=ix3l3
          if(lcas3(4).eq.1) then
              jbf3l34=ibf3l4
              ix3l34=ix3l4
          endif
c---
          call shift3l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l1),bl(jbf2l2),bl(jbf2l3),bl(jbf2l4),
     *                 bl(jbf3l12),bl(jbf3l34),lqmx,
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *           bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),bl(ix2l2),bl(ix2l3),bl(ix2l4),
     *                 bl(ix3l12),bl(ix3l34),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
cc
c- 4 l-shell
       if (lshellt.eq.4) then
c-
          jbuf  = nbuf
          jbuf2 =nbuf2
          jbfij1=ibfij1
          jbfij2=ibfij2
          jbfkl1=ibfkl1
          jbfkl2=ibfkl2
c-
          jbfij3=ibfij3
          jbfkl3=ibfkl3
          jbf2l1=ibf2l1
          jbf2l2=ibf2l2
          jbf2l3=ibf2l3
          jbf2l4=ibf2l4
c
          jbf3l1=ibf3l1
          jbf3l2=ibf3l2
          jbf3l3=ibf3l3
          jbf3l4=ibf3l4
c-
          jssss =issss
c---
          call shift4l(bl(jbuf),bl(jbuf2),l01,l02,
     *                 bl(jbfij1),bl(jbfij2),bl(jbfkl1),bl(jbfkl2),
     *                 lqij,lqkl,
     *                 bl(jbfij3),bl(jbfkl3),
     *                 bl(jbf2l1),bl(jbf2l2),bl(jbf2l3),bl(jbf2l4),
     *                 bl(jbf3l1),bl(jbf3l2),bl(jbf3l3),bl(jbf3l4),lqmx,
     *                 bl(jssss),
     *                 bl(iwij),lsmx, bl(ivij),bl(iuij),bl(isij),lsjl,
     *           bl(ixij),nfu(nqi+1),lnkl,bl(iyij),bl(izij),nfu(nqj+1),
     *                 bl(ix2l1),bl(ix2l2),bl(ix2l3),bl(ix2l4),
     *                 bl(ix3l1),bl(ix3l2),bl(ix3l3),bl(ix3l4),mnbls,
     *                 bl(ixabn),bl(ixcdn),
     *                 bl(indxx),lni,lnj,lnk,lnl)
c-
          call retmem(3)
          return
       endif
c----------------------
      return
      end
c=====================================================
      subroutine convr3(bl,m,nbls,npij,npkl,idx1,idx2,
     *                   xab,xcd, ixabn,ixcdn)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      dimension bl(*)
      dimension idx1(*),idx2(*)
      dimension xab(npij,3),xcd(npkl,3)
c
      nbls1=nbls
      nbls2=nbls*2
      nbls3=nbls*3
      nbls1=nbls1*m
      nbls2=nbls2*m
      nbls3=nbls3*m
      call getmem(nbls3,ixabn)
      call getmem(nbls3,ixcdn)
c
       ixab1=ixabn-1
       ixcd1=ixcdn-1
c
      ijklnmr=0
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
      klpar=idx2(ijkl)
c
      xab1=xab(ijpar,1)
      xab2=xab(ijpar,2)
      xab3=xab(ijpar,3)
      xcd1=xcd(klpar,1)
      xcd2=xcd(klpar,2)
      xcd3=xcd(klpar,3)
c
        do 100 nmr=1,m
        ijklnmr=ijklnmr+1
        bl(ixab1+ijklnmr)      =xab1
        bl(ixab1+ijklnmr+nbls1)=xab2
        bl(ixab1+ijklnmr+nbls2)=xab3
c
        bl(ixcd1+ijklnmr)      =xcd1
        bl(ixcd1+ijklnmr+nbls1)=xcd2
        bl(ixcd1+ijklnmr+nbls2)=xcd3
c
  100 continue
      return
      end
c=====================================================
      subroutine daxpy3(n,a,z1,z2,z3,y1,y2,y3,x)
c------------------------------------------------
c* performs the vector operations with a stride=1
c*
c*     Z = Y + A*X
c*
c*   for three values of a and three matrices Z and Y and this same X
c*
c------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension z1(n),z2(n),z3(n),y1(n),y2(n),y3(n),x(n),a(n,3)
c
      do 10 i=1,n
      z1(i)=y1(i) + a(i,1)*x(i)
      z2(i)=y2(i) + a(i,2)*x(i)
      z3(i)=y3(i) + a(i,3)*x(i)
   10 continue
c*
      end
c=====================================================
      subroutine noshift(lt1,lt2,mnbls,
     *                   buf,buf2,
     *                   bfijx,      bfklx,
     *                   bf2lx,
     *                   lt3,lt4,
     *                   indxx,ln12,ln34)
c**
c**   special cases where the shifting of angular
c**   momentum is not needed at all / (ds|ss)..(ss|sd),
c**   (xs|ys),(xs|sy),(sx|ys),(sx|sy) /
c**
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
cc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbeg,klbeg
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
cc
      dimension buf2(mnbls,lt1,lt2),
c    * bfij1(mnbls,lt3,lt2),bfij2(mnbls,lt3,lt2),
c    * bfkl1(mnbls,lt1,lt4),bfkl2(mnbls,lt1,lt4),
c    * bf2l1(mnbls,lt3,lt4),bf2l2(mnbls,lt3,lt4),
c    * bf2l3(mnbls,lt3,lt4),bf2l4(mnbls,lt3,lt4)
c---
     * bfijx(mnbls,lt3,lt2),
     * bfklx(mnbls,lt1,lt4),
     * bf2lx(mnbls,lt3,lt4)
       dimension buf(mnbls,*)
       dimension indxx(ln12,ln34)
c-----------------------------------
c
       ijb1=ijbeg-1
       klb1=klbeg-1
c
         ijkl=0
         do 5031 i=1,lni
         ii=(i-1)*lnj
         do 5031 j=1,lnj
         ij=ii+j
         do 5031 k=1,lnk
         kk=(k-1)*lnl
         do 5031 l=1,lnl
         kl=kk+l
         ijkl=ijkl+1
ccccccc  indxx(ij+ijb1,kl+klb1)=ijkl
         indxx(ij     ,kl     )=ijkl
 5031    continue
c
         do 5034 ij=ijbeg,lnij
         do 5034 kl=klbeg,lnkl
         ijkl=indxx(ij-ijb1,kl-klb1)
            do 5034 i=1,mnbls
         buf(i,ijkl)=buf2(i,ij,kl)
 5034    continue
c
           if(lshelkl.eq.1 .or. lshelkl.eq.2) then
              do 5035 ij=ijbeg,lnij
              ijkl=indxx(ij-ijb1,1)
                 do 5035 i=1,mnbls
c--->         buf(i,ijkl)=bfkl1(i,ij,1)
              buf(i,ijkl)=bfklx(i,ij,1)
 5035         continue
           endif
c-------
           if(lshelij.eq.1 .or. lshelij.eq.2) then
              do 6035 kl=klbeg,lnkl
              ijkl=indxx(1,kl-klb1)
                 do 6035 i=1,mnbls
c----->       buf(i,ijkl)=bfij1(i,1,kl)
              buf(i,ijkl)=bfijx(i,1,kl)
 6035         continue
           endif
c---------
           if(lshellt.eq.2) then
                ijkl=indxx(1,1)
                   do 1010 i=1,mnbls
c------>           buf(i,ijkl)=bf2l1(i,1,1)
                   buf(i,ijkl)=bf2lx(i,1,1)
 1010              continue
cc
           endif
      return
      end
c=======================================================================
      subroutine shift0l(buf,buf2,lt1,lt2,
     *                   wij,lt3,lsjl,xij,lt4,lt5,lt6,mnbls,xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1)
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common /logic4/ nfu(1)
c
      dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
      dimension xij(mnbls,lt4,lt5,lt6)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c------------------------------------
c
      nibeg=nfu(nqi)+1
      niend=nfu(nqi+1)
      njbeg=nfu(nqj)+1
      njend=nfu(nqj+1)
c
      nkbeg=nfu(nqk)+1
      nkend=nfu(nqk+1)
      nlbeg=nfu(nql)+1
      nlend=nfu(nql+1)
      ijcount=(niend-nibeg+1)*mnbls
c
      do 100 nkl=nfu(nqkl)+1,nfu(nskl+1)
c
         if(nsij.eq.nqij) then  !  ..the (i,s| case ..; nj=1
c           do 107 nj=njbeg,njend  ! only nj=1
c           do 107 ni=nibeg,niend
c              call tfer(buf2(1,ni,nkl),xij(1,ni,1,nkl),mnbls)
c 107       continue
            call dcopy(ijcount,buf2(1,nibeg,nkl),1,xij(1,nibeg,1,nkl),1)
         else
            call horiz12_ab(buf2(1,1,nkl),
     *                   wij,lt3,lsjl,xab,mnbls,nqi,nqj,nsij1,
     *                   xij(1,1,1,nkl),lt4)
         endif
c
  100 continue
c
c---------------------------------------------------------------------
c this part shifts angular momentum from position 3 to position 4
c---------------------------------------------------------------------
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqk)+1,nfu(nskl+1)
c           call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
            call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  301    continue
c
         if(nskl.gt.nqkl) then
            call horiz12(wij,lt3,lsjl,xcd,mnbls,nqk,nql,nskl1)
         endif
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c
      end
c=======================================================================
      subroutine shift1l(buf,
     *                buf2,lt1,lt2,bijx,bklx,lt3,lt4,
     *                wij,lt5,vij,lt6,
     *                xij,lt7,lt8,yij,lt9,mnbls,xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 1 l-shell is present somewhere
c*          in position 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension buf2(mnbls,lt1,lt2),
     * bijx(mnbls,lt3,lt2),
     * bklx(mnbls,lt1,lt4)
      dimension wij(mnbls,lt5,lt6),vij(mnbls,lt5,lt6)
      dimension xij(mnbls,lt7,lt3,lt8),yij(mnbls,lt7,lt9,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
c
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       icount=(niend-nibeg+1)*mnbls
       nibe2=nibeg
       icoun2=icount
c
       nqix=nqi
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          icount=niend*mnbls
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
ccccc        nqlx=1
          endif
       endif
c
c
c*****
c
      do 10 nkl=nfu(nqkl)+1,nfu(nskl+1)
c
       do 11 ij=nqi,nsij
c      ijbeg=nfu(ij)+1
c      ijend=nfu(ij+1)
c          do 12 nij=ijbeg,ijend
c          call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
c  12      continue
       nij=nfu(ij)+1
       ijcount=(nfu(ij+1)-nfu(ij))*mnbls
       call dcopy(ijcount,buf2(1,nij,nkl),1,wij(1,nij,1),1)
c
   11  continue
c
cccccccc
       call horiz12(wij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
cccccccc
c
       do 16 nj=nfu(nqj)+1,nfu(nqj+1)
c      do 16 ni=nfu(nqi)+1,nfu(nqi+1)
c      call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
       call dcopy(icoun2,wij(1,nibe2,nj),1,xij(1,nibe2,nj,nkl),1)
   16  continue
c
ccc  here lshelij can be eq. 0, 1 or 2 only  ccc
c
      if(lshelij.gt.0) then
          if(lshelij.eq.1) then
              if(jtyp.eq.1) then
c               call tfer(bijx(1,1,nkl),xij(1,1,1,nkl),mnbls)
                call dcopy(mnbls,bijx(1,1,nkl),1,xij(1,1,1,nkl),1)
              else
                call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bijx(1,2,nkl),bijx(1,3,nkl),bijx(1,4,nkl),bijx(1,1,nkl))
              endif
          else
c              do 17 nij=nfu(nqi )+1,nfu(nqi +1)
c              call tfer(bijx(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
c  17          continue
           call dcopy(icoun2,bijx(1,nibe2,nkl),1,xij(1,nibe2,1,nkl),1)
          endif
      endif
c
   10 continue
c******
ccc  here lshelkl can be eq. 0, 1 or 2 only  ccc
      if(lshelkl.gt.0) then
ccc
      do 100 nkl=nfu(nqklx)+1,nfu(nqkl+1)
c
       do 101 ij=nqix,nsij
c      ijbeg=nfu(ij)+1
c      ijend=nfu(ij+1)
c
c            do 103 nij=ijbeg,ijend
c            call tfer(bklx(1,nij,nkl),vij(1,nij,1),mnbls)
c 103        continue
       nij=nfu(ij)+1
       ijcount=(nfu(ij+1)-nfu(ij))*mnbls
       call dcopy(ijcount,bklx(1,nij,nkl),1,vij(1,nij,1),1)
  101  continue
c
cccccc
        call horiz12(vij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
cccccc
            do 1071 nj=njbeg,njend
c           do 1071 ni=nibeg,niend
c           call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
            call dcopy(icount,vij(1,nibeg,nj),1,yij(1,nibeg,nj,nkl),1)
 1071       continue
c
  100 continue
c
      endif
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
c        call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  301    continue
c
        call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
c
      if(lshelkl.gt.0)  then
c
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
c            call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
             call dcopy(mnbls,yij(1,ni,nj,1),1,wij(1,1,1),1)
            else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *    yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
            endif
         else
             do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
c            call tfer(yij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
             call dcopy(mnbls,yij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  312        continue
         endif
      endif
c
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c****
c
      return
      end
c=============================================================
      subroutine shift2l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l12,b2l34,
     *                wij,lt5,vij,uij,sij,lt6,
     *                xij,lt7,lt8,yij,zij,lt9,
     *                x2l,mnbls,xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 2 l-shells are present somewhere
c*          in position 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension x2l(mnbls,lt3,lt4)
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l12(mnbls,lt3,lt4),
     * b2l34(mnbls,lt3,lt4)
       dimension wij(mnbls,lt5,lt6),
     * vij(mnbls,lt5,lt6),uij(mnbls,lt5,lt6),sij(mnbls,lt5,lt6)
      dimension xij(mnbls,lt7,lt3,lt8),yij(mnbls,lt7,lt9,lt4),
     *                             zij(mnbls,lt7,lt9,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
       icount=(niend-nibeg+1)*mnbls
       icoun2=icount
       nibe2=nibeg
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          icount=niend*mnbls
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
ccccc        nqlx=1
          endif
       endif
       icoun3=(nfu(nqij+1)-nfu(nqix))*mnbls
       nibe3=nfu(nqix)+1
c
c
c*****
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
       do 101 ij=nqix,nsij
c      ijbeg=nfu(ij)+1
c      ijend=nfu(ij+1)
c          do 102 nij=ijbeg,ijend
c          call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
c 102      continue
        nij=nfu(ij)+1
        ijcount=(nfu(ij+1)-nfu(ij))*mnbls
        call dcopy(ijcount,buf2(1,nij,nkl),1,wij(1,nij,1),1)
c
        if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
c          do 103 nij=ijbeg,ijend
c          call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
c 103      continue
           call dcopy(ijcount,bkl1(1,nij,nkl),1,vij(1,nij,1),1)
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
c          do 104 nij=ijbeg,ijend
c          call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
c 104      continue
           call dcopy(ijcount,bkl2(1,nij,nkl),1,uij(1,nij,1),1)
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
c          do 105 nij=ijbeg,ijend
c          call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
c 105      continue
           call dcopy(ijcount,bkl3(1,nij),1,sij(1,nij,1),1)
         endif
        endif
c
  101  continue
c
ccccccccccccccccccccccc
c
c
       call horiz12(wij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
c
      if( nkl.le.nfu(nqkl+1)) then
          if(lshelkl.eq.1.or.lshelkl.eq.3) then
            call horiz12(vij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
          if(lshelkl.eq.2.or.lshelkl.eq.3) then
            call horiz12(uij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
      if(nkl.eq.1) then
          if(lshelkl.eq.3) then
            call horiz12(sij,lt5,lt6,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
cccc
c
       do 107 nj=njbeg,njend
c      do 107 ni=nibeg,niend
c      call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
          call dcopy(icount,wij(1,nibeg,nj),1,xij(1,nibeg,nj,nkl),1)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
                do 1071 nj=njbeg,njend
c               do 1071 ni=nibeg,niend
c               call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
            call dcopy(icount,vij(1,nibeg,nj),1,yij(1,nibeg,nj,nkl),1)
 1071           continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
                do 1072 nj=njbeg,njend
c               do 1072 ni=nibeg,niend
c               call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
            call dcopy(icount,uij(1,nibeg,nj),1,zij(1,nibeg,nj,nkl),1)
 1072           continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
c            call tfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
             call dcopy(mnbls,bij1(1,1,nkl),1,xij(1,1,1,nkl),1)
          else
      call daxpy3(mnbls,xab,
     *          xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *          bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
c          do 108 nij=nfu(nqi )+1,nfu(nqi +1)
c          call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
c 108      continue
           call dcopy(icoun2,bij2(1,nibe2,nkl),1,xij(1,nibe2,1,nkl),1)
       endif
       if(lshelij.eq.3) then
c          call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
           call dcopy(mnbls,bij3(1,nkl),1,xij(1,1,1,nkl),1)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
                if(lcas2(1).eq.1 .or. lcas2(2).eq.1) then
c                  call tfer(b2l12(1,1,nkl),x2l(1,1,nkl),mnbls)
                   call dcopy(mnbls,b2l12(1,1,nkl),1,x2l(1,1,nkl),1)
                endif
              else
                if(lcas2(1).eq.1 .or. lcas2(2).eq.1) then
                   call daxpy3(mnbls,xab,
     *             x2l(1,2,nkl),x2l(1,3,nkl),x2l(1,4,nkl),
     *             b2l12(1,2,nkl),b2l12(1,3,nkl),b2l12(1,4,nkl),
     *             b2l12(1,1,nkl))
                endif
              endif
           endif
           if(lshelij.eq.2) then
              if(lcas2(3).eq.1 .or. lcas2(4).eq.1) then
c                 do 109 nij=nfu(nqix)+1,nfu(nqij+1)
c                 call tfer(b2l34(1,nij,nkl),x2l(1,nij,nkl),mnbls)
c 109             continue
              call dcopy(icoun3,b2l34(1,nibe3,nkl),1,x2l(1,nibe3,nkl),1)
              endif
           endif
c
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
c        call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  301    continue
c
cccccccccc
c
       call horiz12(wij,lt5,lt6,xcd,mnbls,nqk,nql,nskl1)
cccccccccc
c
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
c          call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
           call dcopy(mnbls,yij(1,ni,nj,1),1,wij(1,1,1),1)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
c        call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,zij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
c        call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
         call dcopy(mnbls,sij(1,ni,nj),1,wij(1,1,1),1)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2 l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
              if(lcas2(1).eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
c                  call tfer(x2l(1,nj,1),wij(1,1,1),mnbls)
                   call dcopy(mnbls,x2l(1,nj,1),1,wij(1,1,1),1)
               endif
              endif
              if(lcas2(3).eq.1) then
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
c                  call tfer(x2l(1,ni,1),wij(1,1,1),mnbls)
                   call dcopy(mnbls,x2l(1,ni,1),1,wij(1,1,1),1)
               endif
              endif
            else
              if(lcas2(1).eq.1.and.ni.eq.1.and.nj.ge.nfu(nqj)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *            x2l(1,nj,2),x2l(1,nj,3),x2l(1,nj,4),x2l(1,nj,1))
              endif
              if(lcas2(3).eq.1.and.nj.eq.1.and.ni.ge.nfu(nqi)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *            x2l(1,ni,2),x2l(1,ni,3),x2l(1,ni,4),x2l(1,ni,1))
              endif
            endif
         endif
         if(lshelkl.eq.2) then
              do 313 nkl=nfu(nqkx)+1,nfu(nqkl+1)
              if(lcas2(2).eq.1.and.ni.eq.1) then
c                 call tfer(x2l(1,nj,nkl),wij(1,nkl,1),mnbls)
                  call dcopy(mnbls,x2l(1,nj,nkl),1,wij(1,nkl,1),1)
              endif
              if(lcas2(4).eq.1.and.nj.eq.1) then
c                 call tfer(x2l(1,ni,nkl),wij(1,nkl,1),mnbls)
                  call dcopy(mnbls,x2l(1,ni,nkl),1,wij(1,nkl,1),1)
              endif
  313         continue
         endif
       endif
c****
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c
      return
      end
c=============================================================
      subroutine shift3l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l1,b2l2,b2l3,b2l4,
     *                b3l12,b3l34,lt5,
     *                wij,lt6,vij,uij,sij,lt7,
     *                xij,lt8,lt9,yij,zij,lt10,
     *                x2l1,x2l2,x2l3,x2l4,x3l12,x3l34,mnbls,
     *                xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*
c*  when 3 l-shells are present somewhere
c*          in positions 1,2 or 3,4
c***********************************************************
      implicit real*8 (a-h,o-z)
cxxx
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
c
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
c
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
c
cxx
      common /logic4/ nfu(1)
cxx
c
      dimension x2l1(mnbls,lt3,lt4),x2l2(mnbls,lt3,lt4),
     *          x2l3(mnbls,lt3,lt4),x2l4(mnbls,lt3,lt4)
      dimension x3l12(mnbls,lt4),x3l34(mnbls,lt3)
c
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l1(mnbls,lt3,lt4),b2l2(mnbls,lt3,lt4),
     * b2l3(mnbls,lt3,lt4),b2l4(mnbls,lt3,lt4),
     * b3l12(mnbls,lt5),b3l34(mnbls,lt5)
      dimension wij(mnbls,lt6,lt7),
     * vij(mnbls,lt6,lt7),uij(mnbls,lt6,lt7),sij(mnbls,lt6,lt7)
      dimension xij(mnbls,lt8,lt3,lt9),yij(mnbls,lt8,lt10,lt4)
     *                            ,zij(mnbls,lt8,lt10,lt4)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       nibeg=nfu(nqi)+1
       niend=nfu(nqi+1)
       njbeg=nfu(nqj)+1
       njend=nfu(nqj+1)
c
       nkbeg=nfu(nqk)+1
       nkend=nfu(nqk+1)
       nlbeg=nfu(nql)+1
       nlend=nfu(nql+1)
c
       nqix=nqi
c
       icount=(niend-nibeg+1)*mnbls
       icoun2=icount
       nibe2=nibeg
       mnbls3=3*mnbls
ccccc  nqjx=nqj
       if(ityp.eq.3) then
          nibeg=1
          icount=niend*mnbls
          if(jtyp.le.3) nqix=1
       endif
       if(jtyp.eq.3) then
          njbeg=1
ccccc     if(ityp.le.3) nqjx=1
       endif
c
       nqkx=nqk
ccccc  nqlx=nql
       nqklx=nqkl
       if(ktyp.eq.3) then
          nkbeg=1
          if(ltyp.le.3) then
             nqklx=1
             nqkx=1
          endif
       endif
       if(ltyp.eq.3) then
          nlbeg=1
          if(ktyp.le.3) then
             nqklx=1
ccccc        nqlx=1
          endif
       endif
       icoun3=(nfu(nqij+1)-nfu(nqix))*mnbls
       nibe3=nfu(nqix)+1
c
c*****
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
       do 101 ij=nqix,nsij
c      ijbeg=nfu(ij)+1
c      ijend=nfu(ij+1)
           nij=nfu(ij)+1
           ijcount=(nfu(ij+1)-nfu(ij))*mnbls
c          do 102 nij=ijbeg,ijend
c          call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
c 102      continue
           call dcopy(ijcount,buf2(1,nij,nkl),1,wij(1,nij,1),1)
c
       if( nkl.le.nfu(nqkl+1)) then
         if(lshelkl.eq.1.or.lshelkl.eq.3) then
c          do 103 nij=ijbeg,ijend
c          call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
c 103      continue
           call dcopy(ijcount,bkl1(1,nij,nkl),1,vij(1,nij,1),1)
         endif
         if(lshelkl.eq.2.or.lshelkl.eq.3) then
c          do 104 nij=ijbeg,ijend
c          call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
c 104      continue
           call dcopy(ijcount,bkl2(1,nij,nkl),1,uij(1,nij,1),1)
         endif
         if(lshelkl.eq.3.and.nkl.eq.1) then
c          do 105 nij=ijbeg,ijend
c          call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
c 105      continue
           call dcopy(ijcount,bkl3(1,nij),1,sij(1,nij,1),1)
         endif
       endif
c
  101  continue
c
cccccccccccc
c
       call horiz12(wij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
c
      if( nkl.le.nfu(nqkl+1)) then
          if(lshelkl.eq.1.or.lshelkl.eq.3) then
            call horiz12(vij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
          if(lshelkl.eq.2.or.lshelkl.eq.3) then
            call horiz12(uij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
      if(nkl.eq.1) then
          if(lshelkl.eq.3) then
            call horiz12(sij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
          endif
      endif
cccccccccccc
c
       do 107 nj=njbeg,njend
c      do 107 ni=nibeg,niend
c      call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
       call dcopy(icount,wij(1,nibeg,nj),1,xij(1,nibeg,nj,nkl),1)
  107  continue
       if( nkl.le.nfu(nqkl+1)) then
           if(lshelkl.eq.1.or.lshelkl.eq.3) then
             do 1071 nj=njbeg,njend
c            do 1071 ni=nibeg,niend
c            call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
             call dcopy(icount,vij(1,nibeg,nj),1,yij(1,nibeg,nj,nkl),1)
 1071        continue
           endif
           if(lshelkl.eq.2.or.lshelkl.eq.3) then
             do 1072 nj=njbeg,njend
c            do 1072 ni=nibeg,niend
c            call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
             call dcopy(icount,uij(1,nibeg,nj),1,zij(1,nibeg,nj,nkl),1)
 1072        continue
           endif
       endif
c
       if(lshelij.eq.1.or.lshelij.eq.3) then
          if(jtyp.eq.1) then
c            call tfer(bij1(1,1,nkl),xij(1,1,1,nkl),mnbls)
             call dcopy(mnbls,bij1(1,1,nkl),1,xij(1,1,1,nkl),1)
          else
      call daxpy3(mnbls,xab,
     *    xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *    bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
          endif
       endif
       if(lshelij.eq.2.or.lshelij.eq.3) then
c          do 108 nij=nfu(nqi )+1,nfu(nqi +1)
c          call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
c 108      continue
           call dcopy(icoun2,bij2(1,nibe2,nkl),1,xij(1,nibe2,1,nkl),1)
       endif
       if(lshelij.eq.3) then
c          call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
           call dcopy(mnbls,bij3(1,nkl),1,xij(1,1,1,nkl),1)
       endif
c*****
      if( nkl.le.nfu(nqkl+1)) then
c****  2 or 3 l-shells ****
c
           if(lshelij.eq.1) then
              if(jtyp.eq.1) then
c                  call tfer(b2l1(1,1,nkl),x2l1(1,1,nkl),mnbls)
c                  call tfer(b2l2(1,1,nkl),x2l2(1,1,nkl),mnbls)
                   call dcopy(mnbls,b2l1(1,1,nkl),1,x2l1(1,1,nkl),1)
                   call dcopy(mnbls,b2l2(1,1,nkl),1,x2l2(1,1,nkl),1)
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(3).eq.1) then
c                  call tfer(b3l34(1,1),x3l34(1,1),mnbls)
                   call dcopy(mnbls,b3l34(1,1),1,x3l34(1,1),1)
                endif
              else
                call daxpy3(mnbls,xab,
     *          x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *          b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),b2l1(1,1,nkl))
cc
                call daxpy3(mnbls,xab,
     *          x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *          b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),b2l2(1,1,nkl))
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(3).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x3l34(1,2),x3l34(1,3),x3l34(1,4),
     *            b3l34(1,2),b3l34(1,3),b3l34(1,4),b3l34(1,1))
                endif
              endif
           endif
           if(lshelij.eq.2) then
c               do 109 nij=nfu(nqix)+1,nfu(nqij+1)
c                  call tfer(b2l3(1,nij,nkl),x2l3(1,nij,nkl),mnbls)
c                  call tfer(b2l4(1,nij,nkl),x2l4(1,nij,nkl),mnbls)
c 109           continue
              call dcopy(icoun3,b2l3(1,nibe3,nkl),1,x2l3(1,nibe3,nkl),1)
              call dcopy(icoun3,b2l4(1,nibe3,nkl),1,x2l4(1,nibe3,nkl),1)
c-- test ?      if(nkl.eq.1) then
                if(nkl.eq.1 .and. lcas3(4).eq.1) then
c                  do 110 nij=nfu(nqi )+1,nfu(nqi +1)
c                  call tfer(b3l34(1,nij),x3l34(1,nij),mnbls)
c 110              continue
                call dcopy(icoun2,b3l34(1,nibe2),1,x3l34(1,nibe2),1)
                endif
           endif
           if(lshelij.eq.3) then
                if(lcas2(1).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *            b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),
     *            b2l1(1,1,nkl))
cc
c                 call tfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
c                 call tfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
c                 call tfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
                  call dcopy(mnbls3,b2l3(1,2,nkl),1,x2l3(1,2,nkl),1)
                endif
                if(lcas2(2).eq.1) then
                  call daxpy3(mnbls,xab,
     *            x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *            b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),
     *            b2l2(1,1,nkl))
cc
c                 call tfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
c                 call tfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
c                 call tfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
                  call dcopy(mnbls3,b2l4(1,2,nkl),1,x2l4(1,2,nkl),1)
                endif
                if(lcas3(1).eq.1) then
c                 call tfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                  call dcopy(mnbls,b3l12(1,nkl),1,x3l12(1,nkl),1)
                endif
                if(lcas3(2).eq.1) then
c                 call tfer(b3l12(1,nkl),x3l12(1,nkl),mnbls)
                  call dcopy(mnbls,b3l12(1,nkl),1,x3l12(1,nkl),1)
                endif
           endif
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
c        call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  301    continue
c
ccccccc
        call horiz12(wij,lt6,lt7,xcd,mnbls,nqk,nql,nskl1)
ccccccc
cc
       if(lshelkl.eq.1.or.lshelkl.eq.3) then
         if(ltyp.eq.1) then
c          call tfer(yij(1,ni,nj,1),wij(1,1,1),mnbls)
           call dcopy(mnbls,yij(1,ni,nj,1),1,wij(1,1,1),1)
         else
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *   yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
         endif
       endif
       if(lshelkl.eq.2.or.lshelkl.eq.3) then
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
c        call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,zij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  312    continue
       endif
c
       if(lshelkl.eq.3) then
c        call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
         call dcopy(mnbls,sij(1,ni,nj),1,wij(1,1,1),1)
       endif
c
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
c
c****  2,3  l-shells ****
         if(lshelkl.eq.1) then
            if(ltyp.eq.1) then
               if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
c                 call tfer(x2l1(1,nj,1),wij(1,1,1),mnbls)
                  call dcopy(mnbls,x2l1(1,nj,1),1,wij(1,1,1),1)
               endif
               if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
c                 call tfer(x2l3(1,ni,1),wij(1,1,1),mnbls)
                  call dcopy(mnbls,x2l3(1,ni,1),1,wij(1,1,1),1)
               endif
c--test ?      if(ni.eq.1.and.nj.eq.1) then
               if(ni.eq.1.and.nj.eq.1 .and. lcas3(1).eq.1) then
c                 call tfer(x3l12(1,1),wij(1,1,1),mnbls)
                  call dcopy(mnbls,x3l12(1,1),1,wij(1,1,1),1)
               endif
            else
              if(ni.eq.1.and.nj.ge.nfu(nqj)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
              endif
              if(nj.eq.1.and.ni.ge.nfu(nqi)+1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
              endif
c--test ?     if(ni.eq.1.and.nj.eq.1) then
              if(ni.eq.1.and.nj.eq.1 .and. lcas3(1).eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *               x3l12(1,2),x3l12(1,3),x3l12(1,4),x3l12(1,1))
              endif
            endif
c
         endif
         if(lshelkl.eq.2) then
              do 313 nkl=nfu(nqkx)+1,nfu(nqkl+1)
              if(ni.eq.1)  then
c                 call tfer(x2l2(1,nj,nkl),wij(1,nkl,1),mnbls)
                  call dcopy(mnbls,x2l2(1,nj,nkl),1,wij(1,nkl,1),1)
              endif
              if(nj.eq.1)  then
c                 call tfer(x2l4(1,ni,nkl),wij(1,nkl,1),mnbls)
                  call dcopy(mnbls,x2l4(1,ni,nkl),1,wij(1,nkl,1),1)
              endif
c--test ?     if(ni.eq.1.and.nj.eq.1) then
              if(ni.eq.1.and.nj.eq.1 .and. lcas3(2).eq.1) then
c                 call tfer(x3l12(1,nkl),wij(1,nkl,1),mnbls)
                  call dcopy(mnbls,x3l12(1,nkl),1,wij(1,nkl,1),1)
              endif
  313         continue
         endif
c
c****  3 l-shells ****
         if(lshelkl.eq.3) then
              if(lcas2(1).eq.1.and.ni.eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
cc
c             call tfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
c             call tfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
c             call tfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
              call dcopy(mnbls,x2l2(1,nj,2),1,wij(1,2,1),1)
              call dcopy(mnbls,x2l2(1,nj,3),1,wij(1,3,1),1)
              call dcopy(mnbls,x2l2(1,nj,4),1,wij(1,4,1),1)
              endif
              if(lcas2(3).eq.1.and.nj.eq.1) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
c             call tfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
c             call tfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
c             call tfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              call dcopy(mnbls,x2l4(1,ni,2),1,wij(1,2,1),1)
              call dcopy(mnbls,x2l4(1,ni,3),1,wij(1,3,1),1)
              call dcopy(mnbls,x2l4(1,ni,4),1,wij(1,4,1),1)
              endif
c
              if(lcas3(3).eq.1.and.ni.eq.1) then
c                 call tfer(x3l34(1,nj),wij(1,1,1),mnbls)
                  call dcopy(mnbls,x3l34(1,nj),1,wij(1,1,1),1)
              endif
              if(lcas3(4).eq.1.and.nj.eq.1) then
c                 call tfer(x3l34(1,ni),wij(1,1,1),mnbls)
                  call dcopy(mnbls,x3l34(1,ni),1,wij(1,1,1),1)
              endif
         endif
       endif
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c
      return
      end
c=============================================================
      subroutine shift4l(buf,
     *                buf2,lt1,lt2,bij1,bij2,bkl1,bkl2,lt3,lt4,
     *                bij3,bkl3,b2l1,b2l2,b2l3,b2l4,
     *                b3l1,b3l2,b3l3,b3l4,lt5,ssss,
     *                wij,lt6,vij,uij,sij,lt7,
     *                xij,lt8,lt9,yij,zij,lt10,
     *                x2l1,x2l2,x2l3,x2l4,x3l1,x3l2,x3l3,x3l4,mnbls,
     *                xab,xcd,
     *                indxx,lni1,lnj1,lnk1,lnl1)
c***********************************************************
c*        when 4 l-shells are present
c***********************************************************
      implicit real*8 (a-h,o-z)
ctest
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
ctest
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common /logic4/ nfu(1)
c
      dimension x2l1(mnbls,lt3,lt4),x2l2(mnbls,lt3,lt4),
     *          x2l3(mnbls,lt3,lt4),x2l4(mnbls,lt3,lt4)
      dimension x3l1(mnbls,lt4),x3l2(mnbls,lt4),
     *          x3l3(mnbls,lt3),x3l4(mnbls,lt3)
c
      dimension buf2(mnbls,lt1,lt2),
     * bij1(mnbls,lt3,lt2),bij2(mnbls,lt3,lt2),
     * bkl1(mnbls,lt1,lt4),bkl2(mnbls,lt1,lt4),
     * bij3(mnbls,lt2),bkl3(mnbls,lt1),
     * b2l1(mnbls,lt3,lt4),b2l2(mnbls,lt3,lt4),
     * b2l3(mnbls,lt3,lt4),b2l4(mnbls,lt3,lt4),
cc-> * b3l(mnbls,lt5,4)
     * b3l1(mnbls,lt5),b3l2(mnbls,lt5),b3l3(mnbls,lt5),b3l4(mnbls,lt5)
      dimension wij(mnbls,lt6,lt7),
     * vij(mnbls,lt6,lt7),uij(mnbls,lt6,lt7),sij(mnbls,lt6,lt7)
      dimension xij(mnbls,lt8,lt3,lt9),yij(mnbls,lt8,lt10,lt4)
     *                            ,zij(mnbls,lt8,lt10,lt4)
c???? dimension ssss(nbls)
      dimension ssss(mnbls)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c***********************************************************
c
       niend=nfu(nqi+1)
       njend=nfu(nqj+1)
c
       nkend=nfu(nqk+1)
       nlend=nfu(nql+1)
c
       nqix=1
ccccc  nqjx=1
       nibeg=1
       njbeg=1
c
       nqkx=1
ccccc  nqlx=1
       nqklx=1
       nkbeg=1
       nlbeg=1
c
       icount=(niend-nibeg+1)*mnbls
       icoun2=(nfu(nqi+1)-nfu(nqi))*mnbls
       nibe2=nfu(nqi)+1
       icoun3=(nfu(nqij+1)-nfu(nqix))*mnbls
       nibe3=nfu(nqix)+1
       mnbls3=3*mnbls
c*****
c     do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
      do 100 nkl=1,10
c
c      do 101 ij=nqix,nsij
       do 101 ij=1,3
c      ijbeg=nfu(ij)+1
c      ijend=nfu(ij+1)
c          do 102 nij=ijbeg,ijend
c          call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
c 102      continue
       nij=nfu(ij)+1
       ijcount=(nfu(ij+1)-nfu(ij))*mnbls
       call dcopy(ijcount,buf2(1,nij,nkl),1,wij(1,nij,1),1)
c
       if( nkl.le.nfu(nqkl+1)) then
c          do 103 nij=ijbeg,ijend
c          call tfer(bkl1(1,nij,nkl),vij(1,nij,1),mnbls)
c 103      continue
           call dcopy(ijcount,bkl1(1,nij,nkl),1,vij(1,nij,1),1)
c          do 104 nij=ijbeg,ijend
c          call tfer(bkl2(1,nij,nkl),uij(1,nij,1),mnbls)
c 104      continue
           call dcopy(ijcount,bkl2(1,nij,nkl),1,uij(1,nij,1),1)
         if(                 nkl.eq.1) then
c          do 105 nij=ijbeg,ijend
c          call tfer(bkl3(1,nij),sij(1,nij,1),mnbls)
c 105      continue
           call dcopy(ijcount,bkl3(1,nij),1,sij(1,nij,1),1)
         endif
       endif
c
  101  continue
c
ccccccccccc
c
        call horiz12(wij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
c
      if( nkl.le.nfu(nqkl+1)) then
            call horiz12(vij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
            call horiz12(uij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
      endif
      if( nkl.eq.1) then
            call horiz12(sij,lt6,lt7,xab,mnbls,nqi,nqj,nsij1)
      endif
ccccccccccc
c
        do 107 nj=njbeg,njend
c       do 107 ni=nibeg,niend
c       call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
        call dcopy(icount,wij(1,nibeg,nj),1,xij(1,nibeg,nj,nkl),1)
  107   continue
c
      if( nkl.le.nfu(nqkl+1)) then
        do 1071 nj=njbeg,njend
c       do 1071 ni=nibeg,niend
c       call tfer(vij(1,ni,nj),yij(1,ni,nj,nkl),mnbls)
c       call tfer(uij(1,ni,nj),zij(1,ni,nj,nkl),mnbls)
        call dcopy(icount,vij(1,nibeg,nj),1,yij(1,nibeg,nj,nkl),1)
        call dcopy(icount,uij(1,nibeg,nj),1,zij(1,nibeg,nj,nkl),1)
 1071   continue
      endif
c
      call daxpy3(mnbls,xab,
     *   xij(1,1,2,nkl),xij(1,1,3,nkl),xij(1,1,4,nkl),
     *   bij1(1,2,nkl),bij1(1,3,nkl),bij1(1,4,nkl),bij1(1,1,nkl))
c
c          do 108 nij=nfu(nqi )+1,nfu(nqi +1)
c          call tfer(bij2(1,nij,nkl),xij(1,nij,1,nkl),mnbls)
c 108      continue
           call dcopy(icoun2,bij2(1,nibe2,nkl),1,xij(1,nibe2,1,nkl),1)
c
c          call tfer(bij3(1,nkl),xij(1,1,1,nkl),mnbls)
           call dcopy(mnbls,bij3(1,nkl),1,xij(1,1,1,nkl),1)
c
c*****
      if( nkl.le.nfu(nqkl+1)) then
      call daxpy3(mnbls,xab,x2l1(1,2,nkl),x2l1(1,3,nkl),x2l1(1,4,nkl),
     *    b2l1(1,2,nkl),b2l1(1,3,nkl),b2l1(1,4,nkl),b2l1(1,1,nkl))
c               call tfer(b2l3(1,2,nkl),x2l3(1,2,nkl),mnbls)
c               call tfer(b2l3(1,3,nkl),x2l3(1,3,nkl),mnbls)
c               call tfer(b2l3(1,4,nkl),x2l3(1,4,nkl),mnbls)
                call dcopy(mnbls3,b2l3(1,2,nkl),1,x2l3(1,2,nkl),1)
cc
      call daxpy3(mnbls,xab,x2l2(1,2,nkl),x2l2(1,3,nkl),x2l2(1,4,nkl),
     *  b2l2(1,2,nkl),b2l2(1,3,nkl),b2l2(1,4,nkl),b2l2(1,1,nkl))
cc
c               call tfer(b2l4(1,2,nkl),x2l4(1,2,nkl),mnbls)
c               call tfer(b2l4(1,3,nkl),x2l4(1,3,nkl),mnbls)
c               call tfer(b2l4(1,4,nkl),x2l4(1,4,nkl),mnbls)
                call dcopy(mnbls3,b2l4(1,2,nkl),1,x2l4(1,2,nkl),1)
c
c               call tfer(b3l1(1,nkl),x3l1(1,nkl),mnbls)
c               call tfer(b3l2(1,nkl),x3l2(1,nkl),mnbls)
                call dcopy(mnbls,b3l1(1,nkl),1,x3l1(1,nkl),1)
                call dcopy(mnbls,b3l2(1,nkl),1,x3l2(1,nkl),1)
c
      endif
      if(nkl.eq.1) then
c
      call daxpy3(mnbls,xab,x3l3(1,2),x3l3(1,3),x3l3(1,4),
     *  b3l3(1,2),b3l3(1,3),b3l3(1,4),b3l3(1,1))
c
c               do 1101 nij=nfu(nqix)+1,nfu(nqij+1)
c               call tfer(b3l4(1,nij ),x3l4(1,nij),mnbls)
c1101           continue
            call dcopy(icoun3,b3l4(1,nibe3),1,x3l4(1,nibe3),1)
      endif
c
  100 continue
c
c***********************************************************
c*    this part    shifts the angular momentum
c*         from position 3 to position 4
c*                         or
c*         from position 4 to position 3
c***********************************************************
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
c        call tfer(xij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  301    continue
c
ccccccccc
c
           call horiz12(wij,lt6,lt7,xcd,mnbls,nqk,nql,nskl1)
c
cccccc
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * yij(1,ni,nj,2),yij(1,ni,nj,3),yij(1,ni,nj,4),yij(1,ni,nj,1))
cc
         do 312 nkl=nfu(nqkx)+1,nfu(nqkl+1)
c        call tfer(zij(1,ni,nj,nkl),wij(1,nkl,1),mnbls)
         call dcopy(mnbls,zij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
  312    continue
c
c        call tfer(sij(1,ni,nj),wij(1,1,1),mnbls)
         call dcopy(mnbls,sij(1,ni,nj),1,wij(1,1,1),1)
c
c****
       if( ni.le.nfu(nqi +1).and.nj.le.nfu(nqj +1) ) then
              if(ni.eq.1 .and. nj.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x2l1(1,nj,2),x2l1(1,nj,3),x2l1(1,nj,4),x2l1(1,nj,1))
cc
c             call tfer(x2l2(1,nj,2),wij(1,2,1),mnbls)
c             call tfer(x2l2(1,nj,3),wij(1,3,1),mnbls)
c             call tfer(x2l2(1,nj,4),wij(1,4,1),mnbls)
              call dcopy(mnbls,x2l2(1,nj,2),1,wij(1,2,1),1)
              call dcopy(mnbls,x2l2(1,nj,3),1,wij(1,3,1),1)
              call dcopy(mnbls,x2l2(1,nj,4),1,wij(1,4,1),1)
cc
              endif
              if(nj.eq.1 .and. ni.ge.2) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     * x2l3(1,ni,2),x2l3(1,ni,3),x2l3(1,ni,4),x2l3(1,ni,1))
cc
c             call tfer(x2l4(1,ni,2),wij(1,2,1),mnbls)
c             call tfer(x2l4(1,ni,3),wij(1,3,1),mnbls)
c             call tfer(x2l4(1,ni,4),wij(1,4,1),mnbls)
              call dcopy(mnbls,x2l4(1,ni,2),1,wij(1,2,1),1)
              call dcopy(mnbls,x2l4(1,ni,3),1,wij(1,3,1),1)
              call dcopy(mnbls,x2l4(1,ni,4),1,wij(1,4,1),1)
              endif
c
c             if(ni.eq.1) call tfer(x3l3(1,nj),wij(1,1,1),mnbls)
c             if(nj.eq.1) call tfer(x3l4(1,ni),wij(1,1,1),mnbls)
              if(ni.eq.1) call dcopy(mnbls,x3l3(1,nj),1,wij(1,1,1),1)
              if(nj.eq.1) call dcopy(mnbls,x3l4(1,ni),1,wij(1,1,1),1)
       endif
       if( ni.eq.1.and.nj.eq.1 ) then
      call daxpy3(mnbls,xcd,wij(1,1,2),wij(1,1,3),wij(1,1,4),
     *  x3l1(1,2),x3l1(1,3),x3l1(1,4),x3l1(1,1))
cc
c         call tfer(x3l2(1,2),wij(1,2,1),mnbls)
c         call tfer(x3l2(1,3),wij(1,3,1),mnbls)
c         call tfer(x3l2(1,4),wij(1,4,1),mnbls)
          call dcopy(mnbls3,x3l2(1,2),1,wij(1,2,1),1)
cc
c         call tfer(ssss(1),wij(1,1,1),mnbls)
          call dcopy(mnbls,ssss(1),1,wij(1,1,1),1)
       endif
c
ccccccccccccccccccccccccccccccccccc
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c
      return
      end
c=======================================================================
      subroutine horiz12_ab(buf2,wij,lw1,lw2,xab,mnbls,nqi,nqj,nsij1,
     *                      xij,lt4)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij(mnbls,lw1,lw2),xab(mnbls,3)
      dimension buf2(mnbls,*),xij(mnbls,lt4,*)
c-----------------------------------------------------------------------
c for the case with nqj=2
c
      if(nqj.eq.2) then
         j=2
         jbeg=nfu(j)+1
         jend=nfu(j+1)
            do nj=jbeg,jend
            njm=ilast(nj)
            kcr=icool(nj)
               do i=nsij1-j,nqi,-1
               ibeg=nfu(i)+1
               iend=nfu(i+1)
                  do ni=ibeg,iend
                  nij=npxyz(kcr,ni)
                     do n=1,mnbls
                     xij(n,ni,nj)=buf2(n,nij) + xab(n,kcr)*buf2(n,ni)
                     enddo
                  enddo
               enddo
            enddo
         return
      endif
c-----------------------------------------------------------------------
c first only j=2 (p-fun) so njm=1 (s-fun)
c
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
      j=2
      jbeg=nfu(j)+1
      jend=nfu(j+1)
         do nj=jbeg,jend
         njm=ilast(nj)
         kcr=icool(nj)
            do i=nsij1-j,nqi,-1
            ibeg=nfu(i)+1
            iend=nfu(i+1)
               do ni=ibeg,iend
               nij=npxyz(kcr,ni)
                  do n=1,mnbls
                  wij(n,ni,nj)=buf2(n,nij) + xab(n,kcr)*buf2(n,ni)
                  enddo
               enddo
            enddo
         enddo
c-----------------------------------------------------------------------
c now j=3,nqj-1
c
      do j=3,nqj-1
      jbeg=nfu(j)+1
      jend=nfu(j+1)
         do nj=jbeg,jend
         njm=ilast(nj)
         kcr=icool(nj)
            do i=nsij1-j,nqi,-1
            ibeg=nfu(i)+1
            iend=nfu(i+1)
               do ni=ibeg,iend
               nij=npxyz(kcr,ni)
                  do n=1,mnbls
                  wij(n,ni,nj)=wij(n,nij,njm)+xab(n,kcr)*wij(n,ni,njm)
                  enddo
               enddo
            enddo
         enddo
      enddo
c-----------------------------------------------------------------------
c now j=nqj    ; the last step
c
      j=nqj
      jbeg=nfu(j)+1
      jend=nfu(j+1)
         do nj=jbeg,jend
         njm=ilast(nj)
         kcr=icool(nj)
            do i=nsij1-j,nqi+1,-1
            ibeg=nfu(i)+1
            iend=nfu(i+1)
               do ni=ibeg,iend
               nij=npxyz(kcr,ni)
                  do n=1,mnbls
                  wij(n,ni,nj)=wij(n,nij,njm)+xab(n,kcr)*wij(n,ni,njm)
                  enddo
               enddo
            enddo
            i=nqi
            ibeg=nfu(i)+1
            iend=nfu(i+1)
               do ni=ibeg,iend
               nij=npxyz(kcr,ni)
                  do n=1,mnbls
                  xij(n,ni,nj)=wij(n,nij,njm)+xab(n,kcr)*wij(n,ni,njm)
                  enddo
               enddo
         enddo   ! over nj=jbeg,jend
c
      end
c=======================================================================
      subroutine horiz12(wij,lw1,lw2,xab,mnbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij(mnbls,lw1,lw2),xab(mnbls,3)
c---------------------------------------------------
          do 110 j=2,nqj
             jbeg=nfu(j)+1
             jend=nfu(j+1)
             do 115 nj=jbeg,jend
                njm=ilast(nj)
                kcr=icool(nj)
                do 120 i=nsij1-j,nqi,-1
                   ibeg=nfu(i)+1
                   iend=nfu(i+1)
                   do 125 ni=ibeg,iend
                      nij=npxyz(kcr,ni)
                      do 130 n=1,mnbls
        wij(n,ni,nj)=wij(n,nij,njm)+ xab(n,kcr) *wij(n,ni,njm)
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine shift0l_der1(buf,buf2,buf0,lt1,lt2,
     *                   wij,wij0,lt3,lsjl,
     *                   xij,xij0,lt4,lt5,lt6,mnbls,nbls,
     *                   xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1,
     *               ixab_e0,ixab_n0,ixcd_e0,ixcd_n0,
     *               ixab_e0_vec,ixab_n0_vec,
     *               ixcd_e0_vec,ixcd_n0_vec )
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common /logic4/ nfu(1)
c
      dimension buf0(nbls,lt1,lt2),wij0(nbls,lt3,lsjl)
      dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
      dimension xij(mnbls,lt4,lt5,lt6), xij0(nbls,lt4,lt5,lt6)
      dimension buf(mnbls,*)
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
      dimension ixab_e0_vec(3,*),ixab_n0_vec(*)
      dimension ixcd_e0_vec(3,*),ixcd_n0_vec(*)
      dimension ixab_e0(3),ixcd_e0(3)
      dimension ixab_n0(3),ixcd_n0(3)
c------------------------------------
c
      nibeg=nfu(nqi)+1
      niend=nfu(nqi+1)
      icount=(niend-nibeg+1)*mnbls
      icoun2=(niend-nibeg+1)*nbls
      njbeg=nfu(nqj)+1
      njend=nfu(nqj+1)
cc
      nkbeg=nfu(nqk)+1
      nkend=nfu(nqk+1)
      nlbeg=nfu(nql)+1
      nlend=nfu(nql+1)
c
      do 100 nkl=nfu(nqkl)+1,nfu(nskl+1)
c
         if(nsij.eq.nqij) then
            call dcopy(icount,buf2(1,nibeg,nkl),1,xij(1,nibeg,1,nkl),1)
            call dcopy(icoun2,buf0(1,nibeg,nkl),1,xij0(1,nibeg,1,nkl),1)
         else
            call horiz12_der1(buf0(1,1,nkl),buf2(1,1,nkl),
     *                       wij0,wij,lt3,lsjl,xab, nbls,nqi,nqj,nsij1,
     *                       ixab_e0,ixab_n0, ixab_e0_vec,ixab_n0_vec ,
     *                       xij(1,1,1,nkl),xij0(1,1,1,nkl),lt4)
         endif
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
c------------------------------------
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqk)+1,nfu(nskl+1)
c           call tfer( xij(1,ni,nj,nkl), wij(1,nkl,1),mnbls)
c           call tfer(xij0(1,ni,nj,nkl),wij0(1,nkl,1), nbls)
            call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
            call dcopy(nbls,xij0(1,ni,nj,nkl),1,wij0(1,nkl,1),1)
  301    continue
c
         call horiz34_der1(wij0,wij,lt3,lsjl,xcd, nbls,nqk,nql,nskl1,
     *                     ixcd_e0,ixcd_n0, ixcd_e0_vec,ixcd_n0_vec )
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
c     call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
      call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)
  305 continue
  300 continue
c
      end
c======================================================================
      subroutine horiz12_der1(buf0,buf2,
     *                        wij0,wij2,lw1,lw2,xab,nbls,nqi,nqj,nsij1,
     *                        ixab_e0,ixab_n0, ixab_e0_vec,ixab_n0_vec,
     *                        xij2,xij0,lt4)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension buf0(nbls,*),buf2(9,nbls,*)
      dimension xij0(nbls,lt4,*),xij2(9,nbls,lt4,*)
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij2(9,nbls,lw1,lw2),xab(9,nbls,3)
c
      dimension ixab_e0_vec(3,*),ixab_n0_vec(3,*)
      dimension ixab_e0(3)      ,ixab_n0(3)
c
      data zero /0.d0/
c----------------------------------------------------------------------
c order of derivatives :
c
c  wij(1),wij(2),wij(3),wij(4),wij(5),wij(6),wij(7),wij(8),wij(9)
c   Ax     Bx     Cx     Ay     By     Cy     Az     Bx     Cz
c----------------------------------------------------------------------
c the case with j=2 (s-function)
c
       if(nqj.eq.2) then
          j=2
          jbeg=nfu(j)+1
          jend=nfu(j+1)
             do 115 nj=jbeg,jend
             njm=ilast(nj)
             kcr=icool(nj)
             k32=kcr*3-2
             k31=k32+1
                do 120 i=nsij1-j,nqi,-1
                ibeg=nfu(i)+1
                iend=nfu(i+1)
                   do 125 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 1301 ne0=1,ixab_e0(kcr) ! over these with xab(kcr)=0
                         n=ixab_e0_vec(kcr,ne0)
                         xij0(  n,ni,nj)=buf0(  n,nij)
                         xij2(1,n,ni,nj)=buf2(1,n,nij)
                         xij2(2,n,ni,nj)=buf2(2,n,nij)
                         xij2(3,n,ni,nj)=buf2(3,n,nij)
                         xij2(4,n,ni,nj)=buf2(4,n,nij)
                         xij2(5,n,ni,nj)=buf2(5,n,nij)
                         xij2(6,n,ni,nj)=buf2(6,n,nij)
                         xij2(7,n,ni,nj)=buf2(7,n,nij)
                         xij2(8,n,ni,nj)=buf2(8,n,nij)
                         xij2(9,n,ni,nj)=buf2(9,n,nij)
 1301                 continue
                      do 1302 nn0=1,ixab_n0(kcr) ! over these with xab(kcr).ne.0
                         n=ixab_n0_vec(kcr,nn0)
                         xab1=xab(1,n,kcr)
                         xij0(  n,ni,nj)=buf0(  n,nij)+xab1*buf0(  n,ni)
                         xij2(1,n,ni,nj)=buf2(1,n,nij)+xab1*buf2(1,n,ni)
                         xij2(2,n,ni,nj)=buf2(2,n,nij)+xab1*buf2(2,n,ni)
                         xij2(3,n,ni,nj)=buf2(3,n,nij)+xab1*buf2(3,n,ni)
                         xij2(4,n,ni,nj)=buf2(4,n,nij)+xab1*buf2(4,n,ni)
                         xij2(5,n,ni,nj)=buf2(5,n,nij)+xab1*buf2(5,n,ni)
                         xij2(6,n,ni,nj)=buf2(6,n,nij)+xab1*buf2(6,n,ni)
                         xij2(7,n,ni,nj)=buf2(7,n,nij)+xab1*buf2(7,n,ni)
                         xij2(8,n,ni,nj)=buf2(8,n,nij)+xab1*buf2(8,n,ni)
                         xij2(9,n,ni,nj)=buf2(9,n,nij)+xab1*buf2(9,n,ni)
 1302                 continue
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      do 1303 n=1,nbls
                         xij2(k32,n,ni,nj)=xij2(k32,n,ni,nj)+buf0(n,ni)
                         xij2(k31,n,ni,nj)=xij2(k31,n,ni,nj)-buf0(n,ni)
 1303                 continue
  125              continue
  120          continue
  115        continue
          return
       endif
c----------------------------------------------------------------------
c first  part : j=2  ; njm is always s-function :
c
          j=2
          jbeg=nfu(j)+1
          jend=nfu(j+1)
             do 215 nj=jbeg,jend
             njm=ilast(nj)
             kcr=icool(nj)
             k32=kcr*3-2
             k31=k32+1
                do 220 i=nsij1-j,nqi,-1
                ibeg=nfu(i)+1
                iend=nfu(i+1)
                   do 225 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 2301 ne0=1,ixab_e0(kcr) ! over these with xab(kcr)=0
                         n=ixab_e0_vec(kcr,ne0)
                         wij0(  n,ni,nj)=buf0(  n,nij)
                         wij2(1,n,ni,nj)=buf2(1,n,nij)
                         wij2(2,n,ni,nj)=buf2(2,n,nij)
                         wij2(3,n,ni,nj)=buf2(3,n,nij)
                         wij2(4,n,ni,nj)=buf2(4,n,nij)
                         wij2(5,n,ni,nj)=buf2(5,n,nij)
                         wij2(6,n,ni,nj)=buf2(6,n,nij)
                         wij2(7,n,ni,nj)=buf2(7,n,nij)
                         wij2(8,n,ni,nj)=buf2(8,n,nij)
                         wij2(9,n,ni,nj)=buf2(9,n,nij)
 2301                 continue
                      do 2302 nn0=1,ixab_n0(kcr) ! over these with xab(kcr).ne.0
                         n=ixab_n0_vec(kcr,nn0)
                         xab1=xab(1,n,kcr)
                         wij0(  n,ni,nj)=buf0(  n,nij)+xab1*buf0(  n,ni)
                         wij2(1,n,ni,nj)=buf2(1,n,nij)+xab1*buf2(1,n,ni)
                         wij2(2,n,ni,nj)=buf2(2,n,nij)+xab1*buf2(2,n,ni)
                         wij2(3,n,ni,nj)=buf2(3,n,nij)+xab1*buf2(3,n,ni)
                         wij2(4,n,ni,nj)=buf2(4,n,nij)+xab1*buf2(4,n,ni)
                         wij2(5,n,ni,nj)=buf2(5,n,nij)+xab1*buf2(5,n,ni)
                         wij2(6,n,ni,nj)=buf2(6,n,nij)+xab1*buf2(6,n,ni)
                         wij2(7,n,ni,nj)=buf2(7,n,nij)+xab1*buf2(7,n,ni)
                         wij2(8,n,ni,nj)=buf2(8,n,nij)+xab1*buf2(8,n,ni)
                         wij2(9,n,ni,nj)=buf2(9,n,nij)+xab1*buf2(9,n,ni)
 2302                 continue
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      do 2303 n=1,nbls
                         wij2(k32,n,ni,nj)=wij2(k32,n,ni,nj)+buf0(n,ni)
                         wij2(k31,n,ni,nj)=wij2(k31,n,ni,nj)-buf0(n,ni)
 2303                 continue
c
  225              continue
  220          continue
  215        continue
c----------------------------------------------------------------------
c j=3,nqj-1
c
          do 310 j=3,nqj-1
          jbeg=nfu(j)+1
          jend=nfu(j+1)
             do 315 nj=jbeg,jend
             njm=ilast(nj)
             kcr=icool(nj)
             k32=kcr*3-2
             k31=k32+1
                do 320 i=nsij1-j,nqi,-1
                ibeg=nfu(i)+1
                iend=nfu(i+1)
                   do 325 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 3301 ne0=1,ixab_e0(kcr)
                         n=ixab_e0_vec(kcr,ne0)
                         wij0(  n,ni,nj)=wij0(  n,nij,njm)
                         wij2(1,n,ni,nj)=wij2(1,n,nij,njm)
                         wij2(2,n,ni,nj)=wij2(2,n,nij,njm)
                         wij2(3,n,ni,nj)=wij2(3,n,nij,njm)
                         wij2(4,n,ni,nj)=wij2(4,n,nij,njm)
                         wij2(5,n,ni,nj)=wij2(5,n,nij,njm)
                         wij2(6,n,ni,nj)=wij2(6,n,nij,njm)
                         wij2(7,n,ni,nj)=wij2(7,n,nij,njm)
                         wij2(8,n,ni,nj)=wij2(8,n,nij,njm)
                         wij2(9,n,ni,nj)=wij2(9,n,nij,njm)
 3301                 continue
                      do 3302 nn0=1,ixab_n0(kcr)
                         n=ixab_n0_vec(kcr,nn0)
                         xab1=xab(1,n,kcr)
               wij0(  n,ni,nj)=wij0(  n,nij,njm) + xab1*wij0(  n,ni,njm)
               wij2(1,n,ni,nj)=wij2(1,n,nij,njm) + xab1*wij2(1,n,ni,njm)
               wij2(2,n,ni,nj)=wij2(2,n,nij,njm) + xab1*wij2(2,n,ni,njm)
               wij2(3,n,ni,nj)=wij2(3,n,nij,njm) + xab1*wij2(3,n,ni,njm)
               wij2(4,n,ni,nj)=wij2(4,n,nij,njm) + xab1*wij2(4,n,ni,njm)
               wij2(5,n,ni,nj)=wij2(5,n,nij,njm) + xab1*wij2(5,n,ni,njm)
               wij2(6,n,ni,nj)=wij2(6,n,nij,njm) + xab1*wij2(6,n,ni,njm)
               wij2(7,n,ni,nj)=wij2(7,n,nij,njm) + xab1*wij2(7,n,ni,njm)
               wij2(8,n,ni,nj)=wij2(8,n,nij,njm) + xab1*wij2(8,n,ni,njm)
               wij2(9,n,ni,nj)=wij2(9,n,nij,njm) + xab1*wij2(9,n,ni,njm)
 3302                 continue
c
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      do 3303 n=1,nbls
                      wij2(k32,n,ni,nj)=wij2(k32,n,ni,nj)+wij0(n,ni,njm)
                      wij2(k31,n,ni,nj)=wij2(k31,n,ni,nj)-wij0(n,ni,njm)
 3303                 continue
c
  325              continue
  320           continue
  315        continue
  310     continue
c----------------------------------------------------------------------
c j=nqj ; the last step
c
          j=nqj
          jbeg=nfu(j)+1
          jend=nfu(j+1)
             do 415 nj=jbeg,jend
             njm=ilast(nj)
             kcr=icool(nj)
             k32=kcr*3-2
             k31=k32+1
                do 420 i=nsij1-j,nqi+1,-1
                ibeg=nfu(i)+1
                iend=nfu(i+1)
                   do 425 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 4301 ne0=1,ixab_e0(kcr)
                         n=ixab_e0_vec(kcr,ne0)
                         wij0(  n,ni,nj)=wij0(  n,nij,njm)
                         wij2(1,n,ni,nj)=wij2(1,n,nij,njm)
                         wij2(2,n,ni,nj)=wij2(2,n,nij,njm)
                         wij2(3,n,ni,nj)=wij2(3,n,nij,njm)
                         wij2(4,n,ni,nj)=wij2(4,n,nij,njm)
                         wij2(5,n,ni,nj)=wij2(5,n,nij,njm)
                         wij2(6,n,ni,nj)=wij2(6,n,nij,njm)
                         wij2(7,n,ni,nj)=wij2(7,n,nij,njm)
                         wij2(8,n,ni,nj)=wij2(8,n,nij,njm)
                         wij2(9,n,ni,nj)=wij2(9,n,nij,njm)
 4301                 continue
                      do 4302 nn0=1,ixab_n0(kcr)
                         n=ixab_n0_vec(kcr,nn0)
                         xab1=xab(1,n,kcr)
               wij0(  n,ni,nj)=wij0(  n,nij,njm) + xab1*wij0(  n,ni,njm)
               wij2(1,n,ni,nj)=wij2(1,n,nij,njm) + xab1*wij2(1,n,ni,njm)
               wij2(2,n,ni,nj)=wij2(2,n,nij,njm) + xab1*wij2(2,n,ni,njm)
               wij2(3,n,ni,nj)=wij2(3,n,nij,njm) + xab1*wij2(3,n,ni,njm)
               wij2(4,n,ni,nj)=wij2(4,n,nij,njm) + xab1*wij2(4,n,ni,njm)
               wij2(5,n,ni,nj)=wij2(5,n,nij,njm) + xab1*wij2(5,n,ni,njm)
               wij2(6,n,ni,nj)=wij2(6,n,nij,njm) + xab1*wij2(6,n,ni,njm)
               wij2(7,n,ni,nj)=wij2(7,n,nij,njm) + xab1*wij2(7,n,ni,njm)
               wij2(8,n,ni,nj)=wij2(8,n,nij,njm) + xab1*wij2(8,n,ni,njm)
               wij2(9,n,ni,nj)=wij2(9,n,nij,njm) + xab1*wij2(9,n,ni,njm)
 4302                 continue
c
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      do 4303 n=1,nbls
                      wij2(k32,n,ni,nj)=wij2(k32,n,ni,nj)+wij0(n,ni,njm)
                      wij2(k31,n,ni,nj)=wij2(k31,n,ni,nj)-wij0(n,ni,njm)
 4303                 continue
c
  425              continue
  420           continue
c
                i=nqi
                ibeg=nfu(i)+1
                iend=nfu(i+1)
                   do 4251 ni=ibeg,iend
                   nij=npxyz(kcr,ni)
                      do 4304 ne0=1,ixab_e0(kcr)
                         n=ixab_e0_vec(kcr,ne0)
                         xij0(  n,ni,nj)=wij0(  n,nij,njm)
                         xij2(1,n,ni,nj)=wij2(1,n,nij,njm)
                         xij2(2,n,ni,nj)=wij2(2,n,nij,njm)
                         xij2(3,n,ni,nj)=wij2(3,n,nij,njm)
                         xij2(4,n,ni,nj)=wij2(4,n,nij,njm)
                         xij2(5,n,ni,nj)=wij2(5,n,nij,njm)
                         xij2(6,n,ni,nj)=wij2(6,n,nij,njm)
                         xij2(7,n,ni,nj)=wij2(7,n,nij,njm)
                         xij2(8,n,ni,nj)=wij2(8,n,nij,njm)
                         xij2(9,n,ni,nj)=wij2(9,n,nij,njm)
 4304                 continue
                      do 4305 nn0=1,ixab_n0(kcr)
                         n=ixab_n0_vec(kcr,nn0)
                         xab1=xab(1,n,kcr)
               xij0(  n,ni,nj)=wij0(  n,nij,njm) + xab1*wij0(  n,ni,njm)
               xij2(1,n,ni,nj)=wij2(1,n,nij,njm) + xab1*wij2(1,n,ni,njm)
               xij2(2,n,ni,nj)=wij2(2,n,nij,njm) + xab1*wij2(2,n,ni,njm)
               xij2(3,n,ni,nj)=wij2(3,n,nij,njm) + xab1*wij2(3,n,ni,njm)
               xij2(4,n,ni,nj)=wij2(4,n,nij,njm) + xab1*wij2(4,n,ni,njm)
               xij2(5,n,ni,nj)=wij2(5,n,nij,njm) + xab1*wij2(5,n,ni,njm)
               xij2(6,n,ni,nj)=wij2(6,n,nij,njm) + xab1*wij2(6,n,ni,njm)
               xij2(7,n,ni,nj)=wij2(7,n,nij,njm) + xab1*wij2(7,n,ni,njm)
               xij2(8,n,ni,nj)=wij2(8,n,nij,njm) + xab1*wij2(8,n,ni,njm)
               xij2(9,n,ni,nj)=wij2(9,n,nij,njm) + xab1*wij2(9,n,ni,njm)
 4305                 continue
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
                      do 4306 n=1,nbls
                      xij2(k32,n,ni,nj)=xij2(k32,n,ni,nj)+wij0(n,ni,njm)
                      xij2(k31,n,ni,nj)=xij2(k31,n,ni,nj)-wij0(n,ni,njm)
 4306                 continue
 4251              continue
  415        continue
c
c
      end
c======================================================================
      subroutine horiz34_der1(wij0,wij,lw1,lw2,xab,nbls,nqi,nqj,nsij1,
     *                       ixab_e0,ixab_n0, ixab_e0_vec,ixab_n0_vec )
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij0(nbls,lw1,lw2),wij(9,nbls,lw1,lw2),xab(9,nbls,3)
c
      dimension ixab_e0_vec(3,*),ixab_n0_vec(3,*)
      dimension ixab_e0(3)      ,ixab_n0(3)
c
      data zero /0.d0/
c
c----------------------------------------------------------------------
c order of derivatives :
c
c  wij(1),wij(2),wij(3),wij(4),wij(5),wij(6),wij(7),wij(8),wij(9)
c   Ax     Bx     Cx     Ay     By     Cy     Az     Bx     Cz
c----------------------------------------------------------------------
c
      do 210 j=2,nqj
         jbeg=nfu(j)+1
         jend=nfu(j+1)
         do 215 nj=jbeg,jend
            njm=ilast(nj)
            kcr=icool(nj)
            kc3=kcr*3                ! lublin
            do 220 i=nsij1-j,nqi,-1
               ibeg=nfu(i)+1
               iend=nfu(i+1)
               do 225 ni=ibeg,iend
                  nij=npxyz(kcr,ni)
                  do 2301 ne0=1,ixab_e0(kcr) ! over these with xab(kcr)=0
                     n=ixab_e0_vec(kcr,ne0)
                     wij0(n,ni,nj)=wij0(n,nij,njm)
                     wij(1,n,ni,nj)=wij(1,n,nij,njm)
                     wij(2,n,ni,nj)=wij(2,n,nij,njm)
                     wij(3,n,ni,nj)=wij(3,n,nij,njm)
                     wij(4,n,ni,nj)=wij(4,n,nij,njm)
                     wij(5,n,ni,nj)=wij(5,n,nij,njm)
                     wij(6,n,ni,nj)=wij(6,n,nij,njm)
                     wij(7,n,ni,nj)=wij(7,n,nij,njm)
                     wij(8,n,ni,nj)=wij(8,n,nij,njm)
                     wij(9,n,ni,nj)=wij(9,n,nij,njm)
 2301             continue
                  do 2302 nn0=1,ixab_n0(kcr) ! over these with xab(kcr).ne.0
                    n=ixab_n0_vec(kcr,nn0)
                    xab1=xab(1,n,kcr)
                    wij0(n,ni,nj)=wij0(n,nij,njm) + xab1*wij0(n,ni,njm)
                    wij(1,n,ni,nj)=wij(1,n,nij,njm)+xab1*wij(1,n,ni,njm)
                    wij(2,n,ni,nj)=wij(2,n,nij,njm)+xab1*wij(2,n,ni,njm)
                    wij(3,n,ni,nj)=wij(3,n,nij,njm)+xab1*wij(3,n,ni,njm)
                    wij(4,n,ni,nj)=wij(4,n,nij,njm)+xab1*wij(4,n,ni,njm)
                    wij(5,n,ni,nj)=wij(5,n,nij,njm)+xab1*wij(5,n,ni,njm)
                    wij(6,n,ni,nj)=wij(6,n,nij,njm)+xab1*wij(6,n,ni,njm)
                    wij(7,n,ni,nj)=wij(7,n,nij,njm)+xab1*wij(7,n,ni,njm)
                    wij(8,n,ni,nj)=wij(8,n,nij,njm)+xab1*wij(8,n,ni,njm)
                    wij(9,n,ni,nj)=wij(9,n,nij,njm)+xab1*wij(9,n,ni,njm)
 2302             continue
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center C, not over A and B:
c
                  do 2303 n=1, nbls
                    wij(kc3,n,ni,nj)=wij(kc3,n,ni,nj) + wij0(n,ni,njm)
 2303             continue
  225          continue
  220       continue
  215    continue
  210 continue
c
c----------------------------------------------------------------------
c
      end
c======================================================================
c
      subroutine shift0l_der2(buf,buh,
     *                   buf2,buf1,buf0,lt1,lt2,
     *                   wij,wij1,wij0,lt3,lsjl,
     *                   xij,xij1,xij0,lt4,lt5,lt6,mnbls,nbls,
     *                   xab,xcd,
     *                   indxx,lni1,lnj1,lnk1,lnl1)
c------------------------------------
c  when l-shells are not present
c------------------------------------
      implicit real*8 (a-h,o-z)
      common /types/iityp,jjtyp,kktyp,lltyp, ityp,jtyp,ktyp,ltyp
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common/shell/lshellt,lshelij,lshelkl,lhelp,lcas2(4),lcas3(4)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,mmax,
     * nqi,nqj,nqk,nql,nsij,nskl,
     * nqij,nqij1,nsij1,nqkl,nqkl1,nskl1,ijbex,klbex
      common /logic4/ nfu(1)
c
      dimension buf0(nbls,lt1,lt2),wij0(nbls,lt3,lsjl)
      dimension buf1(9*nbls,lt1,lt2),wij1(9*nbls,lt3,lsjl)
      dimension buf2(mnbls,lt1,lt2),wij(mnbls,lt3,lsjl)
      dimension xij(mnbls,lt4,lt5,lt6),xij1(9*nbls,lt4,lt5,lt6),
     *                                 xij0(  nbls,lt4,lt5,lt6)
cccc  dimension buf(mnbls,*)
      dimension buf(9*nbls,*)        ! gradient integrals
      dimension buh(45*nbls,*)       ! hessian  integrals
c
      dimension xab(mnbls,3),xcd(mnbls,3)
      dimension indxx(lni1,lnj1,lnk1,lnl1)
c------------------------------------
c mnbls=45*nbls  ......
c------------------------------------
      nbls9=9*nbls
c
      nibeg=nfu(nqi)+1
      niend=nfu(nqi+1)
      njbeg=nfu(nqj)+1
      njend=nfu(nqj+1)
cc
      nkbeg=nfu(nqk)+1
      nkend=nfu(nqk+1)
      nlbeg=nfu(nql)+1
      nlend=nfu(nql+1)
c
      nqix=nqi
ccccc nqjx=nqj
c
      nqkx=nqk
ccccc nqlx=nql
      nqklx=nqkl
c
      do 100 nkl=nfu(nqklx)+1,nfu(nskl+1)
c
           do 102 nij=nfu(nqix)+1,nfu(nsij+1)
           call tfer(buf2(1,nij,nkl),wij(1,nij,1),mnbls)
           call tfer(buf1(1,nij,nkl),wij1(1,nij,1),nbls9)
           call tfer(buf0(1,nij,nkl),wij0(1,nij,1),nbls)
  102      continue
c
           call horiz12_der2(wij0,wij1,wij,
     *                       lt3,lsjl,xab, nbls,nqi,nqj,nsij1)
c
           do 107 nj=njbeg,njend
           do 107 ni=nibeg,niend
           call tfer(wij(1,ni,nj),xij(1,ni,nj,nkl),mnbls)
           call tfer(wij1(1,ni,nj),xij1(1,ni,nj,nkl),nbls9)
           call tfer(wij0(1,ni,nj),xij0(1,ni,nj,nkl),nbls)
  107      continue
  100 continue
c
c------------------------------------
c this part shifts angular momentum
c   from position 3 to position 4
c------------------------------------
c
      ixyz=0
      do 300 ni=nibeg,niend
      ixyz=ixyz+1
      jxyz=0
      do 300 nj=njbeg,njend
      jxyz=jxyz+1
c
         do 301 nkl=nfu(nqkx)+1,nfu(nskl+1)
c        call tfer( xij(1,ni,nj,nkl), wij(1,nkl,1),mnbls)
c        call tfer(xij1(1,ni,nj,nkl),wij1(1,nkl,1),nbls9)
c        call tfer(xij0(1,ni,nj,nkl),wij0(1,nkl,1), nbls)
         call dcopy(mnbls,xij(1,ni,nj,nkl),1,wij(1,nkl,1),1)
         call dcopy(nbls9,xij1(1,ni,nj,nkl),1,wij1(1,nkl,1),1)
         call dcopy(nbls,xij0(1,ni,nj,nkl),1,wij0(1,nkl,1),1)
  301    continue
c
         call horiz34_der2(wij0,wij1,wij,
     *                     lt3,lsjl,xcd, nbls,nqk,nql,nskl1)
c
      kxyz=0
      do 305 nk=nkbeg,nkend
      kxyz=kxyz+1
      lxyz=0
      do 305 nl=nlbeg,nlend
      lxyz=lxyz+1
      indx=indxx(ixyz,jxyz,kxyz,lxyz)
cold  call tfer(wij(1,nk,nl),buf(1,indx),mnbls)
cold  call dcopy(mnbls,wij(1,nk,nl),1,buf(1,indx),1)  ! only 2ed der.
      call dcopy(nbls9,wij1(1,nk,nl),1,buf(1,indx),1) ! 1st-der
      call dcopy(mnbls, wij(1,nk,nl),1,buh(1,indx),1) ! 2ed-der
  305 continue
  300 continue
c
      end
c=============================================================
      subroutine horiz12_der2(wij0,wij1,wij2,
     *                        lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij1(9,nbls,lw1,lw2)
      dimension wij2(45,nbls,lw1,lw2),xab(45,nbls,3)
c---------------------------------------------------
          do 110 j=2,nqj
             jbeg=nfu(j)+1
             jend=nfu(j+1)
             do 115 nj=jbeg,jend
                njm=ilast(nj)
                kcr=icool(nj)
                do 120 i=nsij1-j,nqi,-1
                   ibeg=nfu(i)+1
                   iend=nfu(i+1)
                   do 125 ni=ibeg,iend
                      nij=npxyz(kcr,ni)
c
                      do 130 n=1, nbls
      wij0(n,ni,nj)=wij0(n,nij,njm) + xab(1,n,kcr)*wij0(n,ni,njm)
                      do m=1,9
      wij1(m,n,ni,nj)=wij1(m,n,nij,njm) + xab(1,n,kcr)*wij1(m,n,ni,njm)
                      enddo
                      do m=1,45
      wij2(m,n,ni,nj)=wij2(m,n,nij,njm) + xab(1,n,kcr)*wij2(m,n,ni,njm)
                      enddo
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center A and B , not C :
c
c  first derivatives:
c
          if(kcr.eq.1) then
             wij1(1,n,ni,nj)=wij1(1,n,ni,nj) + wij0(n,ni,njm)
             wij1(2,n,ni,nj)=wij1(2,n,ni,nj) - wij0(n,ni,njm)
c
          endif
          if(kcr.eq.2) then
             wij1(4,n,ni,nj)=wij1(4,n,ni,nj) + wij0(n,ni,njm)
             wij1(5,n,ni,nj)=wij1(5,n,ni,nj) - wij0(n,ni,njm)
          endif
          if(kcr.eq.3) then
             wij1(7,n,ni,nj)=wij1(7,n,ni,nj) + wij0(n,ni,njm)
             wij1(8,n,ni,nj)=wij1(8,n,ni,nj) - wij0(n,ni,njm)
          endif
c
c second derivatives:
c
          if(kcr.eq.1) then
c block aa:
             wij2(1,n,ni,nj)=wij2(1,n,ni,nj) + wij1(1,n,ni,njm)*2.d0
             wij2(2,n,ni,nj)=wij2(2,n,ni,nj) + wij1(4,n,ni,njm)
             wij2(3,n,ni,nj)=wij2(3,n,ni,nj) + wij1(7,n,ni,njm)
c block ab:
             wij2(7,n,ni,nj) =wij2(7,n,ni,nj) - wij1(1,n,ni,njm)
     *                                        + wij1(2,n,ni,njm)
             wij2(8,n,ni,nj) =wij2(8,n,ni,nj) + wij1(5,n,ni,njm)
             wij2(9,n,ni,nj) =wij2(9,n,ni,nj) + wij1(8,n,ni,njm)
             wij2(10,n,ni,nj)=wij2(10,n,ni,nj)- wij1(4,n,ni,njm)
             wij2(13,n,ni,nj)=wij2(13,n,ni,nj)- wij1(7,n,ni,njm)
c block ac:
             wij2(16,n,ni,nj)=wij2(16,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(17,n,ni,nj)=wij2(17,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(18,n,ni,nj)=wij2(18,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(25,n,ni,nj)=wij2(25,n,ni,nj)- wij1(2,n,ni,njm)*2.d0
             wij2(26,n,ni,nj)=wij2(26,n,ni,nj)- wij1(5,n,ni,njm)
             wij2(27,n,ni,nj)=wij2(27,n,ni,nj)- wij1(8,n,ni,njm)
c block bc:
             wij2(31,n,ni,nj)=wij2(31,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(32,n,ni,nj)=wij2(32,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(33,n,ni,nj)=wij2(33,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
          if(kcr.eq.2) then
c block aa:
             wij2(2,n,ni,nj)=wij2(2,n,ni,nj) + wij1(1,n,ni,njm)
             wij2(4,n,ni,nj)=wij2(4,n,ni,nj) + wij1(4,n,ni,njm)*2.d0
             wij2(5,n,ni,nj)=wij2(5,n,ni,nj) + wij1(7,n,ni,njm)
c block ab:
             wij2(8,n,ni,nj) =wij2(8,n,ni,nj) - wij1(1,n,ni,njm)
             wij2(10,n,ni,nj)=wij2(10,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(11,n,ni,nj)=wij2(11,n,ni,nj)- wij1(4,n,ni,njm)
     *                                        + wij1(5,n,ni,njm)
             wij2(12,n,ni,nj)=wij2(12,n,ni,nj)+ wij1(8,n,ni,njm)
             wij2(14,n,ni,nj)=wij2(14,n,ni,nj)- wij1(7,n,ni,njm)
c block ac:
             wij2(19,n,ni,nj)=wij2(19,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(20,n,ni,nj)=wij2(20,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(21,n,ni,nj)=wij2(21,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(26,n,ni,nj)=wij2(26,n,ni,nj)- wij1(2,n,ni,njm)
             wij2(28,n,ni,nj)=wij2(28,n,ni,nj)- wij1(5,n,ni,njm)*2.d0
             wij2(29,n,ni,nj)=wij2(29,n,ni,nj)- wij1(8,n,ni,njm)
c block bc:
             wij2(34,n,ni,nj)=wij2(34,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(35,n,ni,nj)=wij2(35,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(36,n,ni,nj)=wij2(36,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
          if(kcr.eq.3) then
c block aa:
             wij2(3,n,ni,nj)=wij2(3,n,ni,nj) + wij1(1,n,ni,njm)
             wij2(5,n,ni,nj)=wij2(5,n,ni,nj) + wij1(4,n,ni,njm)
             wij2(6,n,ni,nj)=wij2(6,n,ni,nj) + wij1(7,n,ni,njm)*2.d0
c block ab:
             wij2(9,n,ni,nj) =wij2(9,n,ni,nj) - wij1(1,n,ni,njm)
             wij2(12,n,ni,nj)=wij2(12,n,ni,nj)- wij1(4,n,ni,njm)
             wij2(13,n,ni,nj)=wij2(13,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(14,n,ni,nj)=wij2(14,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(15,n,ni,nj)=wij2(15,n,ni,nj)- wij1(7,n,ni,njm)
     *                                        + wij1(8,n,ni,njm)
c block ac:
             wij2(22,n,ni,nj)=wij2(22,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(23,n,ni,nj)=wij2(23,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(24,n,ni,nj)=wij2(24,n,ni,nj)+ wij1(9,n,ni,njm)
c block bb:
             wij2(27,n,ni,nj)=wij2(27,n,ni,nj)- wij1(2,n,ni,njm)
             wij2(29,n,ni,nj)=wij2(29,n,ni,nj)- wij1(5,n,ni,njm)
             wij2(30,n,ni,nj)=wij2(30,n,ni,nj)- wij1(8,n,ni,njm)*2.d0
c block bc:
             wij2(37,n,ni,nj)=wij2(37,n,ni,nj)- wij1(3,n,ni,njm)
             wij2(38,n,ni,nj)=wij2(38,n,ni,nj)- wij1(6,n,ni,njm)
             wij2(39,n,ni,nj)=wij2(39,n,ni,nj)- wij1(9,n,ni,njm)
c block cc: no contributions of this kind
          endif
  130                 continue
  125              continue
  120          continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine horiz34_der2(wij0,wij1,wij2,
     *                        lw1,lw2,xab,nbls,nqi,nqj,nsij1)
      implicit real*8 (a-h,o-z)
c
      common /logic4/ nfu(1)
      common /logic5/ icoor(1)
      common /logic6/ icool(1)
      common /logic7/ ifrst(1)
      common /logic8/ ilast(1)
      common /logic11/ npxyz(3,1)
c
      dimension wij0(nbls,lw1,lw2)
      dimension wij1(9,nbls,lw1,lw2)
      dimension wij2(45,nbls,lw1,lw2),xab(45,nbls,3)
c---------------------------------------------------
          do 110 j=2,nqj
             jbeg=nfu(j)+1
             jend=nfu(j+1)
             do 115 nj=jbeg,jend
                njm=ilast(nj)
                kcr=icool(nj)
                do 120 i=nsij1-j,nqi,-1
                   ibeg=nfu(i)+1
                   iend=nfu(i+1)
                   do 125 ni=ibeg,iend
                      nij=npxyz(kcr,ni)
                      do 130 n=1, nbls
      wij0(n,ni,nj)=wij0(n,nij,njm) + xab(1,n,kcr)*wij0(n,ni,njm)
                         do m=1,9
      wij1(m,n,ni,nj)=wij1(m,n,nij,njm) + xab(1,n,kcr)*wij1(m,n,ni,njm)
                         enddo
                         do m=1,45
      wij2(m,n,ni,nj)=wij2(m,n,nij,njm) + xab(1,n,kcr)*wij2(m,n,ni,njm)
                         enddo
c
c add an extra term arising from differentiation of the shifting formula
c ONLY for derivatives over center C, not over A and B:
c
c first derivatives:
c
      if(kcr.eq.1) then
         wij1(3,n,ni,nj)=wij1(3,n,ni,nj) + wij0(n,ni,njm)
      endif
      if(kcr.eq.2) then
         wij1(6,n,ni,nj)=wij1(6,n,ni,nj) + wij0(n,ni,njm)
      endif
      if(kcr.eq.3) then
         wij1(9,n,ni,nj)=wij1(9,n,ni,nj) + wij0(n,ni,njm)
      endif
c
c second derivatives:
c
          if(kcr.eq.1) then
c block ac:
             wij2(16,n,ni,nj)=wij2(16,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(19,n,ni,nj)=wij2(19,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(22,n,ni,nj)=wij2(22,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(31,n,ni,nj)=wij2(31,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(34,n,ni,nj)=wij2(34,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(37,n,ni,nj)=wij2(37,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(40,n,ni,nj)=wij2(40,n,ni,nj)+ wij1(3,n,ni,njm)*2.d0
             wij2(41,n,ni,nj)=wij2(41,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(42,n,ni,nj)=wij2(42,n,ni,nj)+ wij1(9,n,ni,njm)
          endif
          if(kcr.eq.2) then
c block ac:
             wij2(17,n,ni,nj)=wij2(17,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(20,n,ni,nj)=wij2(20,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(23,n,ni,nj)=wij2(23,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(32,n,ni,nj)=wij2(32,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(35,n,ni,nj)=wij2(35,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(38,n,ni,nj)=wij2(38,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(41,n,ni,nj)=wij2(41,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(43,n,ni,nj)=wij2(43,n,ni,nj)+ wij1(6,n,ni,njm)*2.d0
             wij2(44,n,ni,nj)=wij2(44,n,ni,nj)+ wij1(9,n,ni,njm)
          endif
          if(kcr.eq.3) then
c block ac:
             wij2(18,n,ni,nj)=wij2(18,n,ni,nj)+ wij1(1,n,ni,njm)
             wij2(21,n,ni,nj)=wij2(21,n,ni,nj)+ wij1(4,n,ni,njm)
             wij2(24,n,ni,nj)=wij2(24,n,ni,nj)+ wij1(7,n,ni,njm)
c block bc:
             wij2(33,n,ni,nj)=wij2(33,n,ni,nj)+ wij1(2,n,ni,njm)
             wij2(36,n,ni,nj)=wij2(36,n,ni,nj)+ wij1(5,n,ni,njm)
             wij2(39,n,ni,nj)=wij2(39,n,ni,nj)+ wij1(8,n,ni,njm)
c block cc:
             wij2(42,n,ni,nj)=wij2(42,n,ni,nj)+ wij1(3,n,ni,njm)
             wij2(44,n,ni,nj)=wij2(44,n,ni,nj)+ wij1(6,n,ni,njm)
             wij2(45,n,ni,nj)=wij2(45,n,ni,nj)+ wij1(9,n,ni,njm)*2.d0
          endif
  130                 continue
  125              continue
  120           continue
  115        continue
  110     continue
c
      end
c=====================================================
      subroutine convr3_forc(bl,m,nbls,npij,npkl,idx1,idx2,
     *                   xab,xcd, ixabn,ixcdn,
     *               ixab_e0,ixab_n0,ixcd_e0,ixcd_n0,
     *               ixab_e0_vec,ixab_n0_vec,
     *               ixcd_e0_vec,ixcd_n0_vec )

      use memory, block => bl

      implicit real*8 (a-h,o-z)
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension idx1(*),idx2(*)
      dimension xab(npij,3),xcd(npkl,3)
      dimension ixab_e0(3),ixcd_e0(3)
      dimension ixab_n0(3),ixcd_n0(3)
      dimension ixab_e0_vec(3,*),ixab_n0_vec(3,*)
      dimension ixcd_e0_vec(3,*),ixcd_n0_vec(3,*)
c
      data zero /0.d0/
c
      nbls1=nbls
      nbls2=nbls*2
      nbls3=nbls*3
      nbls1=nbls1*m
      nbls2=nbls2*m
      nbls3=nbls3*m
      call getmem(nbls3,ixabn)
      call getmem(nbls3,ixcdn)
c
       ixab1=ixabn-1
       ixcd1=ixcdn-1
c
      do i=1,3
         ixab_e0(i)=0
         ixcd_e0(i)=0
         ixab_n0(i)=0
         ixcd_n0(i)=0
      enddo
c
      ijklnmr=0
      do 100 ijkl=1,nbls
      ijpar=idx1(ijkl)
      klpar=idx2(ijkl)
c
      xab1=xab(ijpar,1)
      xab2=xab(ijpar,2)
      xab3=xab(ijpar,3)
      xcd1=xcd(klpar,1)
      xcd2=xcd(klpar,2)
      xcd3=xcd(klpar,3)
c
      if(xab1.eq.zero) then
         ixab_e0(1)=ixab_e0(1)+1
         ixab_e0_vec(1,ixab_e0(1))=ijkl
      else
         ixab_n0(1)=ixab_n0(1)+1
         ixab_n0_vec(1,ixab_n0(1))=ijkl
      endif
      if(xab2.eq.zero) then
         ixab_e0(2)=ixab_e0(2)+1
         ixab_e0_vec(2,ixab_e0(2))=ijkl
      else
         ixab_n0(2)=ixab_n0(2)+1
         ixab_n0_vec(2,ixab_n0(2))=ijkl
      endif
      if(xab3.eq.zero) then
         ixab_e0(3)=ixab_e0(3)+1
         ixab_e0_vec(3,ixab_e0(3))=ijkl
      else
         ixab_n0(3)=ixab_n0(3)+1
         ixab_n0_vec(3,ixab_n0(3))=ijkl
      endif
centers CD :
      if(xcd1.eq.zero) then
         ixcd_e0(1)=ixcd_e0(1)+1
         ixcd_e0_vec(1,ixcd_e0(1))=ijkl
      else
         ixcd_n0(1)=ixcd_n0(1)+1
         ixcd_n0_vec(1,ixcd_n0(1))=ijkl
      endif
      if(xcd2.eq.zero) then
         ixcd_e0(2)=ixcd_e0(2)+1
         ixcd_e0_vec(2,ixcd_e0(2))=ijkl
      else
         ixcd_n0(2)=ixcd_n0(2)+1
         ixcd_n0_vec(2,ixcd_n0(2))=ijkl
      endif
      if(xcd3.eq.zero) then
         ixcd_e0(3)=ixcd_e0(3)+1
         ixcd_e0_vec(3,ixcd_e0(3))=ijkl
      else
         ixcd_n0(3)=ixcd_n0(3)+1
         ixcd_n0_vec(3,ixcd_n0(3))=ijkl
      endif
c
        do 100 nmr=1,m
        ijklnmr=ijklnmr+1
        bl(ixab1+ijklnmr)      =xab1
        bl(ixab1+ijklnmr+nbls1)=xab2
        bl(ixab1+ijklnmr+nbls2)=xab3
c
        bl(ixcd1+ijklnmr)      =xcd1
        bl(ixcd1+ijklnmr+nbls1)=xcd2
        bl(ixcd1+ijklnmr+nbls2)=xcd3
c
  100 continue
c
      end
c=======================================================================
      subroutine make_indxx(indxx,lni,lnj,lnk,lnl)
      dimension indxx(lni,lnj,lnk,lnl)
c
      ncount=0
      do 20 i=1,lni
      do 20 j=1,lnj
      do 20 k=1,lnk
      do 20 l=1,lnl
         ncount=ncount+1
      indxx(i,j,k,l)=ncount
   20 continue
c
      end
c=======================================================================

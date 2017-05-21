c----------------------------------------------------------------
c April 2001 KW The following routins are called for the gradient
c integrals assembling :
c
c for iroute=1 :
c der1_eq1   route=1 nbls=nbls1 (for both firstc= .true. & .false.)
c der1_ne11  route=1 firstc=.true.   nbls>nbls1
c der1_ne12  route=1 firstc=.false.  nbls>nbls1
c
c for iroute=2
c der1_eq21  route=2 nbls=nbls1 firstc=.true.
c der1_eq21  route=2 nbls=nbls1 firstc=.false.
c der1_ne21  route=2 firstc=.true.   nbls>nbls1
c der1_ne22  route=2 firstc=.false.  nbls>nbls1
c----------------------------------------------------------------
c Nov 97, KW : assembling of primitives for gradient integrals has
c              been changed. There used to be three routines to handle
c              contraction of prim. to contracted integrals :
c              der1_NE, der1_EQ and der1_E1.
c              The last one for the case with only one quartet/block
c              has been eliminated.
c
c              The major saving comes from contracting only those
c              primitives which are later needed for derivatives
c              (in first_der). Much fewer primitives (ordinary and
c              rescaled by 2*a, 2*b and 2*c (exponents) ) are
c              contracted now compared to the previous version.
c----------------------------------------------------------------
c All routines of the type name_2 are used when
c          IROUTE=2
c----------------------------------------------------------------
c      ASSEMBLY OF THE 2-EL. INTEGRALS (I+J,0|K+L,0)
c----------------------------------------------------------------
      subroutine assemblx(firstc,nbls,nbls1,l01,l02,
     *                    lci,lcj,lck,lcl,lcij,lckl,npij,npkl)

      use memory

      implicit real*8 (a-h,o-z)
      logical firstc,firstx
c
      common /cpu/intsize,iacc,icache,memreal
      common /route/ iroute
c
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /logic4/ nfu(1)
c
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
c     common /big/ bl(1)
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
c     common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
c2000 needed for iaa,ibb,icc,idd for gradient
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
c
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
c
      common /contr/ ngc1,ngc2,ngc3,ngc4,lc1,lc2,lc3,lc4,lc12,lc34
c
cccc  dimension bl(*)
c----------------------------------------------------------------
c-                  --- for buf2  ---
c
        firstx=firstc
        if(where.eq.'fock' .or. where.eq.'shif') then
           call conbuf2(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                  bl(ibuf2),bl(indx))
        endif
        if(where.eq.'forc') then
           ibut2=ibuf
           if(iroute.eq.1) then
              if(nbls.eq.nbls1) then
                 call der1_eq1(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                         bl(ibut2),bl(indx),
     *                         bl(iaax),bl(ibbx),bl(iccx))
              else
                 nblsnot=nbls-nbls1
                 if(firstx) then
                  call getint(nbls,idxnot)
                  call getint(nbls,idpres)
                  call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
                  call retmem(1)
                  call der1_ne11(nbls,nbls1,bl(iwt0),l01,l02,bl(ibut2),
     *                           bl(idxnot),nblsnot, bl(indx),
     *                           bl(iaax),bl(ibbx),bl(iccx))
                  call retmem(1) ! release idxnot
                  firstx=.false.
                 else
                  call der1_ne12(nbls,nbls1,bl(iwt0),l01,l02,bl(ibut2),
     *                           bl(indx),
     *                           bl(iaax),bl(ibbx),bl(iccx))
                 endif
              endif
           else
c
              if(nbls.eq.nbls1) then
                 if(firstx) then
                   call der1_eq21(nbls,nbls1,bl(iwt0),l01,l02,
     *                            bl(ibut2),
     *                            bl(iaa),lci,bl(ibb),lcj,bl(icc),lck)
                   firstx=.false.
                 else
                   call der1_eq22(nbls,nbls1,bl(iwt0),l01,l02,
     *                                bl(ibut2),
     *                           bl(iaa),lci,bl(ibb),lcj,bl(icc),lck)
                 endif
              else
                 nblsnot=nbls-nbls1
                 if(firstx) then
                  call getint(nbls,idxnot)
                  call getint(nbls,idpres)
                  call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
                  call retmem(1)
                  call der1_ne21(nbls,nbls1,bl(iwt0),l01,l02,bl(ibut2),
     *                           bl(idxnot),nblsnot, bl(indx),
     *                           bl(iaa),lci,bl(ibb),lcj,bl(icc),lck)
                  call retmem(1) ! release idxnot
                  firstx=.false.
                 else
                  call der1_ne22(nbls,nbls1,bl(iwt0),l01,l02,
     *                                bl(ibut2), bl(indx),
     *                           bl(iaa),lci,bl(ibb),lcj,bl(icc),lck)
                 endif
              endif
           endif
        endif
        if(where.eq.'hess') then
           ibut2=ibuf
c2001
           nblsnot=nbls-nbls1
           call getint(nbls,idxnot)
           call getint(nbls,idpres)
           call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
           call retmem(1)
c2001
           call getmem(6*nbls1,iexpo)
           if(iroute.eq.1) then
              call der2_1(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                         bl(ibut2),bl(indx),
     *                         bl(idxnot),nblsnot,
     *                         bl(iaax),bl(ibbx),bl(iccx),bl(iexpo))
           else
              nbl0102=nbls*l01*l02
c             write(6,*)'nbls=',nbls,' l01,l02=',l01,l02,
c    *        ' l01*l02=',l01*l02,
c    *        '  nbls*l01*l02=',nbls*l01*l02,' cache=',icache
              if(icache.gt.10*nbl0102) then
                 call der2_2(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                            bl(ibut2),bl(indx),
     *                            bl(idxnot),nblsnot,
     *                            bl(iaa),lci,bl(ibb),lcj,bl(icc),lck,
     *                            bl(iexpo))
              else
                 call der2_2driv(firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                            bl(ibut2),bl(indx),
     *                            bl(idxnot),nblsnot,
     *                            bl(iaa),lci,bl(ibb),lcj,bl(icc),lck,
     *                            bl(iexpo))
              endif
           endif
           call retmem(1)  ! release iexpo
           call retmem(1)  ! release idxnot
        endif
c
      IF(lshellt.eq.0) go to 100
c
c-----------------------------------------------------------------------
c-                  --- for l-shells ---
c
        firstx=firstc
c
        IF( iroute.eq.1 ) THEN
          call conshel_1(bl,firstx,nbls,nbls1,l01,l02,
     *                   lci,lcj,lck,lcl,lcij,lckl,npij,npkl)
        ELSE
          call conshel_2(bl,firstx,nbls,nbls1,l01,l02,
     *                   lci,lcj,lck,lcl,lcij,lckl)
        ENDIF
c-----------------------------------------------------------------------
c
  100 CONTINUE
c
      firstc=firstx
c
      end
c===============================================================
      subroutine conbuf2(firstc,nbls,nbls1,xt1,lt1,lt2, buf2,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
cnmr
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
cnmr
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension buf2(nbls,lt1,lt2)
C
c-------
      IF(where.eq.'fock') THEN
c
        ijs=nfu(nqij)+1
        IF (FIRSTC) THEN
           DO 501 KL=nfu(nqkl)+1,LNKL
           DO 501 ij=ijs,LNIJ
           do 501 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  501      CONTINUE
C
           FIRSTC=.FALSE.
        ELSE
           DO 601 KL=nfu(nqkl)+1,LNKL
           DO 601 ij=ijs,LNIJ
           do 601 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=BUF2(ijkl,IJ,KL)+XT1(i,IJ,KL)
  601      CONTINUE
        ENDIF
      ENDIF
c
      IF(where.eq.'shif') THEN
c
        ijs=nfu(nqij)+1
        lnijx=nfu(nsij)
        lnklx=nfu(nskl)
        IF (FIRSTC) THEN
           DO 551 KL=nfu(nqkl)+1,nfu(nskl)
           DO 551 ij=ijs,nfu(nsij)
           do 551 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  551      CONTINUE
           DO 552 KL=nfu(nqkl)+1,nfu(nskl)
           DO 552 ij=nfu(nsij)+1,nfu(nsij+1)
           do 552 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  552      CONTINUE
           DO 553 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 553 ij=nfu(nqij)+1,nfu(nsij)
           do 553 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=XT1(i,IJ,KL)
  553      CONTINUE
           FIRSTC=.FALSE.
        ELSE
           DO 651 KL=nfu(nqkl)+1,lnklx
           DO 651 ij=ijs,lnijx
           do 651 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  651      CONTINUE
           DO 652 KL=nfu(nqkl)+1,lnklx
           DO 652 ij=nfu(nsij)+1,nfu(nsij+1)
           do 652 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  652      CONTINUE
           DO 653 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 653 ij=nfu(nqij)+1,nfu(nsij)
           do 653 i=1,nbls1
           ijkl=indx(i)
           BUF2(ijkl,IJ,KL)=buf2(ijkl,ij,kl)+XT1(i,IJ,KL)
  653      CONTINUE
        ENDIF
      ENDIF
c
      end
c=======================================================================
c23456789.123456789.123456789.123456789.123456789.123456789.123456789.12
c=======================================================================
      subroutine der1_ne11(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                          idxnot,nbls2, indx,
     *                          aax,bbx,ccx)
c
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=1 ,  nbls>nbls1
c
c     first=.true. (first cntr. step)
c-----------------------------------------------------------------------
c three regions should be distinguished here :
c 1:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(1,ijkl,ij,kl)=xt1(i,ij,kl)        <-- ordinary integrals
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax(i) <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx(i) <---int.rescaled by 2*b_exp
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx(i) <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 2:
c    do kl=nfu(nskl)+1,nfu(nskl+1)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx(i) <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 3:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nsij)+1,nfu(nsij+1)
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax(i) <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx(i) <---int.rescaled by 2*b_exp
c    enddo
c    enddo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
ccc   logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension idxnot(*)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
cccc  dimension buf2(4,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,4)
      data zero /0.d0/
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
c
ccc   IF (FIRSTC) THEN
c
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,1)=xt1(i,ij,kl)
                  buf2(ijkl,ij,kl,2)=xt1(i,ij,kl)*aax(i)
                  buf2(ijkl,ij,kl,3)=xt1(i,ij,kl)*bbx(i)
                  buf2(ijkl,ij,kl,4)=xt1(i,ij,kl)*ccx(i)
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,1)=zero
                  buf2(ijkl,ij,kl,2)=zero
                  buf2(ijkl,ij,kl,3)=zero
                  buf2(ijkl,ij,kl,4)=zero
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,4)=xt1(i,ij,kl)*ccx(i)
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,4)=zero
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,2)=xt1(i,ij,kl)*aax(i)
                  buf2(ijkl,ij,kl,3)=xt1(i,ij,kl)*bbx(i)
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,2)=zero
                  buf2(ijkl,ij,kl,3)=zero
               enddo
            enddo
         enddo
c
      end
c=======================================================================
      subroutine der1_ne12(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                          indx,
     *                          aax,bbx,ccx)
c
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=1 ,  nbls>nbls1
c
c     first=.false.( not first cntr. step)
c-----------------------------------------------------------------------
c three regions should be distinguished here :
c 1:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(1,ijkl,ij,kl)=xt1(i,ij,kl)        <-- ordinary integrals
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax(i) <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx(i) <---int.rescaled by 2*b_exp
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx(i) <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 2:
c    do kl=nfu(nskl)+1,nfu(nskl+1)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx(i) <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 3:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nsij)+1,nfu(nsij+1)
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax(i) <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx(i) <---int.rescaled by 2*b_exp
c    enddo
c    enddo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
cccc  logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
cccc  dimension buf2(4,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,4)
      data zero /0.d0/
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
c
cccc  IF (FIRSTC) THEN
cccc  ELSE
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,1)=buf2(ijkl,ij,kl,1)+xt1(i,ij,kl)
               buf2(ijkl,ij,kl,2)=buf2(ijkl,ij,kl,2)+xt1(i,ij,kl)*aax(i)
               buf2(ijkl,ij,kl,3)=buf2(ijkl,ij,kl,3)+xt1(i,ij,kl)*bbx(i)
               buf2(ijkl,ij,kl,4)=buf2(ijkl,ij,kl,4)+xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,4)=buf2(ijkl,ij,kl,4)+xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,2)=buf2(ijkl,ij,kl,2)+xt1(i,ij,kl)*aax(i)
               buf2(ijkl,ij,kl,3)=buf2(ijkl,ij,kl,3)+xt1(i,ij,kl)*bbx(i)
               enddo
            enddo
         enddo
cccc  ENDIF
c
      end
c=======================================================================
      subroutine der1_eq1(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *                          aax,bbx,ccx)
c
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=1 , nbls=nbls1
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
cccc  dimension buf2(4,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,4)
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
c
      IF (FIRSTC) THEN
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,1)=xt1(i,ij,kl)
               buf2(i,ij,kl,2)=xt1(i,ij,kl)*aax(i)
               buf2(i,ij,kl,3)=xt1(i,ij,kl)*bbx(i)
               buf2(i,ij,kl,4)=xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,4)=xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               buf2(i,ij,kl,2)=xt1(i,ij,kl)*aax(i)
               buf2(i,ij,kl,3)=xt1(i,ij,kl)*bbx(i)
               enddo
            enddo
         enddo
         FIRSTC=.FALSE.
      ELSE
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,1)=buf2(i,ij,kl,1)+xt1(i,ij,kl)
               buf2(i,ij,kl,2)=buf2(i,ij,kl,2)+xt1(i,ij,kl)*aax(i)
               buf2(i,ij,kl,3)=buf2(i,ij,kl,3)+xt1(i,ij,kl)*bbx(i)
               buf2(i,ij,kl,4)=buf2(i,ij,kl,4)+xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,4)=buf2(i,ij,kl,4)+xt1(i,ij,kl)*ccx(i)
               enddo
            enddo
         enddo
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               buf2(i,ij,kl,2)=buf2(i,ij,kl,2)+xt1(i,ij,kl)*aax(i)
               buf2(i,ij,kl,3)=buf2(i,ij,kl,3)+xt1(i,ij,kl)*bbx(i)
               enddo
            enddo
         enddo
      ENDIF
      end
c=======================================================================
      subroutine der2_1(firstc,nbls,nbls1,xt1,lt1,lt2, buf2,indx,
     *                        idxnot,nblsnot,
     *                        aax,bbx,ccx,expo)
c
cccc  this is called only for where.eq.'hess'
c
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension idxnot(*)
      dimension xt1(nbls1,lt1,lt2)
c
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1),expo(nbls1,6)
C
c2001 dimension buf2(10,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,10)
c
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c end for first derivatives
c               buf2(5,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*b_exp
c               buf2(6,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*c_exp
c               buf2(7,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*c_exp
c               buf2(8,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*a_exp
c               buf2(9,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*b_exp
c               buf2(10nbls,lt1,lt2) - rescaled with 2*c_exp * 2*c_exp
c end for second derivatives
c-------
      data zero /0.d0/
c-------
c multiply exponents :
c
        do i=1,nbls1
          expo(i,1)=aax(i)*bbx(i)
          expo(i,2)=aax(i)*ccx(i)
          expo(i,3)=bbx(i)*ccx(i)
          expo(i,4)=aax(i)*aax(i)
          expo(i,5)=bbx(i)*bbx(i)
          expo(i,6)=ccx(i)*ccx(i)
        enddo
c------------------------------------------------
        ij1=nfu(nqij)
        ij2=nfu(nsij-1)
        ij3=nfu(nsij)
        ij4=nfu(nsij+1)
c
        kl1=nfu(nqkl)
        kl2=nfu(nskl-1)
        kl3=nfu(nskl)
        kl4=nfu(nskl+1)
c------------------------------------------------
        IF (FIRSTC) THEN
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,1)=XT1(i,IJ,KL)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax(i)
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx(i)
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(i,1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(i,2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(i,3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(i,4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(i,5)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(i,6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax(i)
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(i,1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(i,2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(i,3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(i,4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(i,5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(i,1)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(i,4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(i,5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx(i)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(i,2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(i,3)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(i,6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(i,2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(i,3)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(i,6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           FIRSTC=.FALSE.
        ELSE
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,1)=buf2(ijkl,ij,kl,1)+XT1(i,IJ,KL)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax(i)
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx(i)
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(i,1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(i,2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(i,3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(i,4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(i,5)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(i,6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax(i)
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(i,1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(i,2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(i,3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(i,4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(i,5)
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(i,1)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(i,4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(i,5)
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx(i)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(i,2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(i,3)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(i,6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(i,2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(i,3)
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(i,6)
                 enddo
              enddo
           enddo
        ENDIF
c------------------------------------------------
c
      end
c=======================================================================
c
c           Assembling for general contraction.
c
c It uses the transposed buf2 array which is called BUT2 here.
c
c----------------------------------------------------------------
      subroutine assemblg(firstc,nbls,nbls1,l01,l02,
     *                    lci,lcj,lck,lcl,ngcd,ibut2)

      use memory

      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /route/ iroute
      common /runtype/ scftype,where
c
c     common /big/ bl(1)
c
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
c
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
      common /memor5e/ igci,igcj,igck,igcl,indgc,igcoef,
     *                 icfg,jcfg,kcfg,lcfg, igcij,igckl
c new for grad. derivatives:
      common /memor5dd/ iaax,ibbx,iccx
c--------------------------------------------------------
      if(where.eq.'forc') then
         call getmem(ngcd*nbls,igc_ax)
         call getmem(ngcd*nbls,igc_bx)
         call getmem(ngcd*nbls,igc_cx)
         if(iroute.eq.1) then
            call coef_x_exp(nbls,nbls1,bl(indx),ngcd,
     *                      bl(igcoef),bl(iaax),bl(ibbx),bl(iccx),
     *                      bl(igc_ax),bl(igc_bx),bl(igc_cx) )
         else
            call coef_x_exp2(nbls,nbls1,bl(indx),ngcd,
     *                      bl(igcoef),bl(iaa),bl(ibb),bl(icc),
     *                                 lci,lcj,lck,
     *                      bl(igc_ax),bl(igc_bx),bl(igc_cx) )
         endif
c
         if(nbls.gt.nbls1 .and. firstc) then
            nblsnot=nbls-nbls1
            call getint(nbls,idxnot)
            call getint(nbls,idpres)
            call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
            call retmem(1)
           call zerogc(bl(ibut2),ngcd,4,nbls,l01,l02,bl(idxnot),nblsnot)
            call retmem(1)
         endif
         call asselg_der1(firstc,bl(iwt0),l01,l02,nbls, bl(ibut2),
     *                     bl(indx),nbls1, ngcd,bl(igcoef),
     *                     bl(igc_ax),bl(igc_bx),bl(igc_cx) )
c
         call retmem(3)
         return
      endif
c--------------------------------------------------------
      if(where.eq.'hess') then
         call getmem(ngcd*nbls*10,igcexp)
         if(iroute.eq.1) then
            call goef_x_exp1(nbls,nbls1,bl(indx),ngcd,
     *                      bl(igcoef),bl(iaax),bl(ibbx),bl(iccx),
     *                      bl(igcexp) )
         else
            call goef_x_exp2(nbls,nbls1,bl(indx),ngcd,
     *                      bl(igcoef),bl(iaa),bl(ibb),bl(icc),
     *                      lci,lcj,lck, bl(igcexp) )
         endif
c
         if(nbls.gt.nbls1 .and. firstc) then
            nblsnot=nbls-nbls1
            call getint(nbls,idxnot)
            call getint(nbls,idpres)
            call setnot(nbls,nbls1,bl(indx),bl(idxnot),bl(idpres))
            call retmem(1)
          call zerogc(bl(ibut2),ngcd,10,nbls,l01,l02,bl(idxnot),nblsnot)
            call retmem(1)
         endif
c
         call asselg_der2(firstc,bl(iwt0),l01,l02,nbls, bl(ibut2),
     *                     bl(indx),nbls1, ngcd,          bl(igcexp) )
c
         call retmem(1)
         return
      endif
c--------------------------------------------------------
c
      call asselg(firstc,bl(iwt0),l01,l02,nbls, bl(ibut2),
     *            bl(indx),nbls1, ngcd,bl(igcoef) )
c
c--------------------------------------------------------
c
      end
c=======================================================================
      subroutine asselg(firstc,xt1,lt1,lt2,nbls,but2,indx,nbls1,
     *                  ngcd,gcoef)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
cccc  dimension buf2(nbls,lt1,lt2,ngcd)
      dimension but2(ngcd,nbls,lt1,lt2)
c
      dimension gcoef(ngcd,nbls)
c---------------------------------------------------------------------
      ijs=nfu(nqij)+1
      kls=nfu(nqkl)+1
c
      IF (FIRSTC) THEN
        DO 501 KL=kls,LNKL
        DO 501 ij=ijs,LNIJ
        do 501 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
        if(abs(xint).gt.0.d0) then
          do 502 iqu=1,ngcd
          but2(iqu,ijkl,ij,kl)=xint*gcoef(iqu,ijkl)
  502     continue
        else
          do 503 iqu=1,ngcd
          but2(iqu,ijkl,ij,kl)=0.d0
  503     continue
        endif
  501   continue
C
           FIRSTC=.FALSE.
      ELSE
C
        DO 601 KL=kls,LNKL
        DO 601 ij=ijs,LNIJ
        do 601 i=1,nbls1
        ijkl=indx(i)
        xint=xt1(i,ij,kl)
        if(abs(xint).gt.0.d0) then
          do 602 iqu=1,ngcd
          but2(iqu,ijkl,ij,kl)=but2(iqu,ijkl,ij,kl)+xint*gcoef(iqu,ijkl)
  602     continue
        endif
  601   continue
c
      ENDIF
      end
c=======================================================================
      subroutine coef_x_exp(nbls,nbls1,indx,ngcd,
     *                      gcoef,aax,bbx,ccx,
     *                      gc_ax,gc_bx,gc_cx )
      implicit real*8 (a-h,o-z)
      dimension indx(*)
      dimension gcoef(ngcd,nbls)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
      dimension gc_ax(ngcd,nbls),gc_bx(ngcd,nbls),gc_cx(ngcd,nbls)
c
      do i=1,nbls1
         aaxi=aax(i)
         bbxi=bbx(i)
         ccxi=ccx(i)
         ijkl=indx(i)
         do iqu=1,ngcd
            gc_ax(iqu,ijkl)=gcoef(iqu,ijkl)*aaxi
            gc_bx(iqu,ijkl)=gcoef(iqu,ijkl)*bbxi
            gc_cx(iqu,ijkl)=gcoef(iqu,ijkl)*ccxi
         enddo
      enddo
c
      end
c=======================================================================
      subroutine coef_x_exp2(nbls,nbls1,indx,ngcd,
     *                      gcoef,aax,bbx,ccx,lci,lcj,lck,
     *                      gc_ax,gc_bx,gc_cx )
c2000 used only for iroute=2
      implicit real*8 (a-h,o-z)
      dimension indx(*)
      dimension gcoef(ngcd,nbls)
      dimension aax(*),bbx(*),ccx(*)
      dimension gc_ax(ngcd,nbls),gc_bx(ngcd,nbls),gc_cx(ngcd,nbls)
c
         aaxi=aax(lci)
         bbxi=bbx(lcj)
         ccxi=ccx(lck)
c
      do i=1,nbls1
         ijkl=indx(i)
         do iqu=1,ngcd
            gc_ax(iqu,ijkl)=gcoef(iqu,ijkl)*aaxi
            gc_bx(iqu,ijkl)=gcoef(iqu,ijkl)*bbxi
            gc_cx(iqu,ijkl)=gcoef(iqu,ijkl)*ccxi
         enddo
      enddo
c
      end
c=======================================================================
      subroutine asselg_der1(firstc,xt1,lt1,lt2,nbls,but2,indx,nbls1,
     *                       ngcd,gcoef,gc_ax,gc_bx,gc_cx)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
c
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
cccc  dimension buf2(nbls,lt1,lt2,ngcd)
      dimension but2(ngcd,4,nbls,lt1,lt2)
      dimension gcoef(ngcd,nbls)
      dimension gc_ax(ngcd,nbls),gc_bx(ngcd,nbls),gc_cx(ngcd,nbls)
      data zero /0.d0/
c-----------------------------------------------------------------------
c               but2(1,nbls,lt1,lt2) - ordinary contraction
c               but2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               but2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               but2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------------------
c
      ijs=nfu(nqij)+1
      kls=nfu(nqkl)+1
c
      IF(FIRSTC) THEN
        do kl=kls,lnkl
           do ij=ijs,lnij
              do i=1,nbls1
                 ijkl=indx(i)
                 xint=xt1(i,ij,kl)
                 if(abs(xint).gt.zero) then
                    do iqu=1,ngcd
                       but2(iqu,1,ijkl,ij,kl)=xint*gcoef(iqu,ijkl)
                       but2(iqu,2,ijkl,ij,kl)=xint*gc_ax(iqu,ijkl)
                       but2(iqu,3,ijkl,ij,kl)=xint*gc_bx(iqu,ijkl)
                       but2(iqu,4,ijkl,ij,kl)=xint*gc_cx(iqu,ijkl)
                    enddo
                 else
                    do iqu=1,ngcd
                       but2(iqu,1,ijkl,ij,kl)=zero
                       but2(iqu,2,ijkl,ij,kl)=zero
                       but2(iqu,3,ijkl,ij,kl)=zero
                       but2(iqu,4,ijkl,ij,kl)=zero
                    enddo
                 endif
              enddo
           enddo
        enddo
        FIRSTC=.FALSE.
      ELSE
c...... region I
        do kl=kls,nfu(nskl)
           do ij=ijs,nfu(nsij)
              do i=1,nbls1
                 ijkl=indx(i)
                 xint=xt1(i,ij,kl)
                 if(abs(xint).gt.zero) then
                    do iqu=1,ngcd
                       but2(iqu,1,ijkl,ij,kl)=but2(iqu,1,ijkl,ij,kl)
     *                                        + xint*gcoef(iqu,ijkl)
                       but2(iqu,2,ijkl,ij,kl)=but2(iqu,2,ijkl,ij,kl)
     *                                        + xint*gc_ax(iqu,ijkl)
                       but2(iqu,3,ijkl,ij,kl)=but2(iqu,3,ijkl,ij,kl)
     *                                        + xint*gc_bx(iqu,ijkl)
                       but2(iqu,4,ijkl,ij,kl)=but2(iqu,4,ijkl,ij,kl)
     *                                        + xint*gc_cx(iqu,ijkl)
                    enddo
                 endif
              enddo
           enddo
        enddo
c...... region II
        do kl=nfu(nskl)+1,nfu(nskl+1)
           do ij=ijs,nfu(nsij)
              do i=1,nbls1
                 ijkl=indx(i)
                 xint=xt1(i,ij,kl)
                 if(abs(xint).gt.zero) then
                    do iqu=1,ngcd
                       but2(iqu,4,ijkl,ij,kl)=but2(iqu,4,ijkl,ij,kl)
     *                                        + xint*gc_cx(iqu,ijkl)
                    enddo
                 endif
              enddo
           enddo
        enddo
c...... region III
        do kl=kls,nfu(nskl)
           do ij=nfu(nsij)+1,nfu(nsij+1)
              do i=1,nbls1
                 ijkl=indx(i)
                 xint=xt1(i,ij,kl)
                 if(abs(xint).gt.zero) then
                    do iqu=1,ngcd
                       but2(iqu,2,ijkl,ij,kl)=but2(iqu,2,ijkl,ij,kl)
     *                                        + xint*gc_ax(iqu,ijkl)
                       but2(iqu,3,ijkl,ij,kl)=but2(iqu,3,ijkl,ij,kl)
     *                                        + xint*gc_bx(iqu,ijkl)
                    enddo
                 endif
              enddo
           enddo
        enddo
      ENDIF
c
      end
c=======================================================================
c
c======= Duplicated routines for different values of IROUTE ====
c
c----------------------------------------------------------------
c      ASSEMBLY OF THE 2-EL. INTEGRALS (I+J,0|K+L,0)
c               when l shells are present
c              it is called when IROUTE=1
c----------------------------------------------------------------
      subroutine conshel_1(bl,firstc,nbls,nbls1,l01,l02,
     *                     lci,lcj,lck,lcl,lcij,lckl,npij,npkl)
      implicit real*8 (a-h,o-z)
      logical firstc,firstx
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /logic4/ nfu(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
cccc  common /big/ bl(1)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
c dimensions for assembling :
      common /dimasse/ lqij,lqkl,lqmx,lij3,lkl3,l3l,lsss
c
      dimension bl(*)
c----------------------------------------------------------------
        nqijr=nqij
        nqklr=nqkl
      if(where.eq.'shif') then
c-    - for nmr derivatives -
        nqijr=nqij1
        nqklr=nqkl1
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-                     --- for bfij1 s from -> lx/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfij1),lqij,lnkl, bl(icis),npij,lci,
     *           bl(idx1),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-                        --- for bfij2 s from xl/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfij2),lqij,lnkl, bl(icjs),npij,lcj,
     *           bl(idx1),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-                          --- for bfkl1 s from xy/lz ---
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfkl1),lnij,lqkl, bl(icks),npkl,lck,
     *           bl(idx2),bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-                          --- for bfkl2 s from xy/zl ---
c
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *           bl(ibfkl2),lnij,lqkl,  bl(icls),npkl,lcl,
     *           bl(idx2),bl(indx), ijenx,klenx)
      endif
c
      IF(lshellt.eq.1) go to 100
c-----------------------------------------------------------------------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
c-
      if(lshelij.eq.3) then
c-                          --- for bfij3  ss from ll/xy ---
c
                ij3b=1
                kl3b=klbeg
        firstx=firstc
        call conijkl3 (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibfij3),lij3,lnkl, bl(ifij),npij,lcij,
     *                 bl(idx1),bl(indx),ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.3) then
c-                         --- for bfkl3 ss from xy/ll ---
c
                ij3b=ijbeg
                kl3b=1
        firstx=firstc
        call conijkl3 (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibfkl3),lnij,lkl3, bl(ifkl),npkl,lckl,
     *                 bl(idx2),bl(indx), ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(1).eq.1) then
c-                          --- for bf2l1  ss from lx/ly ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l1),lqij,lqkl,
     *               bl(icis),bl(icks),npij,npkl,lci,lck,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(2).eq.1) then
c-                          --- for bf2l2  ss from lx/yl ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l2),lqij,lqkl,
     *               bl(icis),bl(icls),npij,npkl,lci,lcl,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(3).eq.1) then
c-                          --- for bf2l3  ss from xl/ly ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l3),lqij,lqkl,
     *               bl(icjs),bl(icks),npij,npkl,lcj,lck,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(4).eq.1) then
c-                          --- for bf2l4  ss from xl/yl ---
c
        firstx=firstc
        call conb2ln(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(ibf2l4),lqij,lqkl,
     *               bl(icjs),bl(icls),npij,npkl,lcj,lcl,
     *               bl(idx1),bl(idx2),bl(indx))
      endif
c
      IF(lshellt.eq.2) go to 100
c-----------------------------------------------------------------------
c-                      --- for bf3l  ---
c
c     if(lcase.eq.12.or.lcase.eq.16) then
c-
      if(lcas3(1).eq.1) then
c-            --- for bf3l1  sss from ll/lx ---
c
        firstx=firstc
        call conb3la(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icks),bl(ifij),  bl(ibf3l1),l3l,lqmx,
     *               lck,lcij,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.13.or.lcase.eq.16) then
c-
      if(lcas3(2).eq.1) then
c-            --- for bf3l2  sss from ll/xl ---
c
        firstx=firstc
        call conb3la(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icls),bl(ifij),  bl(ibf3l2),l3l,lqmx,
     *               lcl,lcij,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas3(3).eq.1) then
c-            --- for bf3l3  sss from lx/ll ---
c
        firstx=firstc
        call conb3lb(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icis),bl(ifkl),  bl(ibf3l3),lqmx,l3l,
     *               lci,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas3(4).eq.1) then
c-            --- for bf3l4  sss from xl/ll ---
c
        firstx=firstc
        call conb3lb(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *               bl(icjs),bl(ifkl),  bl(ibf3l4),lqmx,l3l,
     *               lcj,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c
      IF(lshellt.eq.3) go to 100
c-----------------------------------------------------------------------
      if(lcase.eq.16) then
c-    -- for ssss(nbls)  ssss from ll/ll --
c
        firstx=firstc
        call conssss (bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                bl(ifij),bl(ifkl), bl(issss),lsss ,
     *                lcij,lckl,npij,npkl,bl(idx1),bl(idx2),bl(indx))
      endif
c
c-----------------------------------------------------------------------
c
  100 CONTINUE
c
      firstc=firstx
c
      end
c===============================================================
      subroutine conijkl1(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij1,lt3,lt4, facti,npij,lci,
     *                    idx1,indx, ijenx,klenx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
      dimension facti(npij,*)
      dimension idx1(*),indx(*)
c
c********
c
      call convr1(bl,nbls1,ifni,facti,lci,npij,idx1,indx)
c
      call assel1(firstc,xt1,lt1,lt2,nbls,
     *            bl(ifni),  bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
c
      call retmem(1)
c
      return
      end
c===============================================================
      subroutine assel1(firstc,xt1,lt1,lt2,nbls,
     *            facti, bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
c
c***
      implicit real*8 (a-h,o-z)
      logical firstc
C***
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
      dimension facti(*)
C
      IF (FIRSTC) THEN
C
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
              do 504 i=1,nbls1
              ijkl=indx(i)
              BFIJ1(ijkl,IJ,KL)=XT1(i,IJ,KL)*FACTI(i)
  504         CONTINUE
C
           FIRSTC=.FALSE.
      ELSE
C
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
              do 604 i=1,nbls1
              ijkl=indx(i)
              BFIJ1(ijkl,IJ,KL)=BFIJ1(ijkl,IJ,KL)+XT1(i,IJ,KL)*FACTI(i)
  604         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conijkl3(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij3,lt3,lt4,  factij,npij,lcij,
     *                    idx1,indx, ij3b,kl3b)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
      dimension factij(npij,*)
      dimension idx1(*),indx(*)
c********
c
      call convr1(bl,nbls1,ifnij,factij,lcij, npij,idx1,indx)
c
      call assel2a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bfij3,lt3,lt4, bl(ifnij), indx, ij3b,kl3b)
c
      call retmem(1)
c
      return
      end
c***************
      subroutine assel2a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij3,lt3,lt4, factij, indx, ij3b,kl3b)
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
      dimension factij(*)
C
      IF (FIRSTC) THEN
c
              do 502 kl=kl3b,lt4
              do 502 ij=ij3b,lt3
              do 502 i=1,nbls1
                 ijkl=indx(i)
                 BFIJ3(ijkl,ij,kl)=XT1(i,ij,KL)*FACTIJ(i)
  502         CONTINUE
c
          FIRSTC=.FALSE.
      ELSE
C
              DO 602 KL=kl3b,lt4
              do 602 ij=ij3b,lt3
              do 602 i=1,nbls1
              ijkl=indx(i)
              BFIJ3(ijkl,ij,kl)=bfij3(ijkl,ij,kl)+XT1(i,ij,KL)*FACTIJ(i)
  602         CONTINUE
C
      ENDIF
      return
      end
c===============================================================
      subroutine conb2ln(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   bf2l1, lt3,lt4,
     *                   facti,factk,npij,npkl,
     *                   lci,lck,idx1,idx2,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
C********************************
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
      dimension facti(npij,*),factk(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c********
c
c     call convr1(bl,nbls1,ifni,facti,lci, npij,idx1,indx)
c     call convr1(bl,nbls1,ifnk,factk,lck, npkl,idx2,indx)
c or
      call convr2(bl,nbls1,ifni ,facti ,lci ,npij,idx1,
     *                  ifnk ,factk ,lck ,npkl,idx2, indx)
c
      call assel2b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bf2l1,lt3,lt4, bl(ifni),bl(ifnk),indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel2b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   bf2l1,lt3,lt4,facti,factk,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
      dimension facti(*),factk(*)
c********
C
             if(where.eq.'shif') then
c            -- for nmr derivatives --
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
                klenx=nfu(nqkl1+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c
      IF (FIRSTC) THEN
c
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
                 do 504 i=1,nbls1
                 ijkl=indx(i)
                 xIJ1=XT1(i,IJ,KL)*FACTI(i)
                 BF2L1(ijkl,ij,kl)=xij1*FACTK(i)
  504         CONTINUE
c*
          FIRSTC=.FALSE.
      ELSE
C
c-----------> DO 604 KL=KLBEG,NFU(NQKL+1)
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
                 do 604 i=1,nbls1
                 ijkl=indx(i)
                 xIJ1=XT1(i,IJ,KL)*FACTI(i)
                 BF2L1(ijkl,ij,kl)=BF2l1(ijkl,ij,kl)+xij1*FACTK(i)
  604         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conb3la(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   factk,factij, bf3l,lt5,lt6,
     *                   lck,lcij,npij,npkl,idx1,idx2,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
c
ccc   common /big/ bl(1)
cc
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension factk(npkl,*), factij(npij,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifnk ,factk ,lck ,npkl,idx2,
     *                  ifnij,factij,lcij,npij,idx1,indx)
c
      call assel3a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bl(ifnk),bl(ifnij), bf3l,lt5,lt6,  indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel3a(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   factk,factij, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension factk(*),factij(*)
c
c********
C
             if(where.eq.'shif') then
c            -- for nmr derivatives --
cccccc          ijenx=nfu(nqij1+1)
c--------nie->  if(nqij1.eq.nsij) ijenx=1
                klenx=nfu(nqkl1+1)
c WHAT ?        if(nqkl1.eq.nskl) klenx=nfu(nqij+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
cccccc          ijenx=nfu(nqij+1)
c--------nie->  if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c--
      IF (FIRSTC) THEN
c
              DO 503 KL=KLBEG,klenx
              do 503 ij=1,lt5
              do 503 i=1,nbls1
                 ijkl=indx(i)
                 XIJ3=XT1(i,ij,KL)*FACTIJ(i)
                 BF3L(ijkl,ij,KL)=xij3*FACTK(i)
  503         CONTINUE
c
          FIRSTC=.FALSE.
      ELSE
c
            DO 603 KL=KLBEG,klenx
            do 603 ij=1,lt5
               do 603 i=1,nbls1
               ijkl=indx(i)
               XIJ3=XT1(i,ij,KL)*FACTIJ(i)
               BF3L(ijkl,ij,KL)=bf3l(ijkl,ij,kl)+xij3*FACTK(i)
  603       CONTINUE
      endif
      return
      end
ccccccc
c===============================================================
      subroutine conb3lb(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   facti,factkl, bf3l,lt5,lt6,
     *                   lci,lckl,npij,npkl,idx1,idx2,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
c
ccc   common /big/ bl(1)
cc
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension facti(npij,*),factkl(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifni ,facti ,lci ,npij,idx1,
     *                  ifnkl,factkl,lckl,npkl,idx2,indx)
c
      call assel3b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bl(ifni),bl(ifnkl), bf3l,lt5,lt6, indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel3b(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   facti,factkl, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2), bf3l(nbls,lt5,lt6)
      dimension facti(*), factkl(*)
c
c********
c--
             if(where.eq.'shif') then
c            -- for nmr derivatives --
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
cccccc          klenx=nfu(nqkl1+1)
c--------nie--> if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
cccccc          klenx=nfu(nqkl+1)
c--------nie--> if(nqkl.eq.nskl) klenx=1
             endif
c--
      IF (FIRSTC) THEN
              DO 511 kl=1,lt6
              DO 511 IJ=IJBEG,ijenx
              do 511 i=1,nbls1
                 ijkl=indx(i)
                 xKL3=XT1(i,IJ,kl)*FACTKL(i)
                 BF3L(ijkl,IJ,kl)=xKL3*FACTI(i)
  511         CONTINUE
c*
          FIRSTC=.FALSE.
      ELSE
C
              DO 611 kl=1,lt6
              DO 611 IJ=IJBEG,ijenx
              do 611 i=1,nbls1
              ijkl=indx(i)
                 xKL3=XT1(i,IJ,kl)*FACTKL(i)
                 BF3L(ijkl,IJ,kl)=BF3L(ijkl,IJ,kl)+xkl3*FACTI(i)
  611         CONTINUE
c
      ENDIF
      return
      end
c===============================================================
      subroutine conssss (bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    factij,factkl, ssss,isdim,
     *                    lcij,lckl,npij,npkl,idx1,idx2,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical firstc
c
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
c--------------------------------------------------------------
c
      dimension factij(npij,*),factkl(npkl,*)
      dimension idx1(*),idx2(*),indx(*)
c
      call convr2(bl,nbls1,ifnij,factij,lcij,npij,idx1,
     *                  ifnkl,factkl,lckl,npkl,idx2,indx)
c
      call assel4(firstc,nbls,nbls1,xt1,lt1,lt2,
     *            bl(ifnij),bl(ifnkl), ssss,isdim, indx)
c
      call retmem(2)
c
      return
      end
c***************
      subroutine assel4(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                  factij,factkl, ssss,isdim,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
C********************************
C
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
      dimension factij(*),factkl(*)
C
      IF (FIRSTC) THEN
c
              do 507 kl=1,isdim
              do 507 ij=1,isdim
              do 507 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)= XT1(i,ij,kl)*FACTIJ(i)*FACTKL(i)
  507         continue
c
           FIRSTC=.FALSE.
      ELSE
C
c------------------------
              do 607 kl=1,isdim
              do 607 ij=1,isdim
              do 607 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)=ssss(ijkl,ij,kl)+
     *                         XT1(i,ij,kl)*FACTIJ(i)*FACTKL(i)
  607         continue
c------------------------
C
      ENDIF
      return
      end
c===============================================================
      subroutine convr1(bl,nbls,ifni,facti,lci,npij,idx1,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension idx1(*),indx(*)
      dimension facti(npij,*)
c
          call getmem(nbls,ifni)
c
       ifni1=ifni-1
       do 100 i=1,nbls
       ijkl=indx(i)
       ijpar=idx1(ijkl)
       bl(ifni1+i)=facti(ijpar,lci)
  100  continue
c
      return
      end
c===============================================================
      subroutine convr2(bl,nbls,ifnk ,factk ,lck ,npkl,idx2,
     *                  ifnij,factij,lcij,npij,idx1,indx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension idx1(*),idx2(*),indx(*)
      dimension factk(npkl,*),factij(npij,*)
c
          call getmem(nbls,ifnk)
          call getmem(nbls,ifnij)
c
       ifnk1=ifnk-1
       ifnij1=ifnij-1
c
       do 100 i=1,nbls
       ijkl=indx(i)
       ijpar=idx1(ijkl)
       klpar=idx2(ijkl)
c
       bl(ifnk1+i)=factk(klpar,lck)
       bl(ifnij1+i)=factij(ijpar,lcij)
  100  continue
c
      return
      end
c===============================================================
c---------------------------------------------------------------
c      ASSEMBLY OF THE 2-EL. INTEGRALS (I+J,0|K+L,0)
c               when l shells are present
c              it is called when IROUTE=2
c---------------------------------------------------------------
      subroutine conshel_2(bl,firstc,nbls,nbls1,l01,l02,
     *                     lci,lcj,lck,lcl,lcij,lckl)
      implicit real*8 (a-h,o-z)
      logical firstc,firstx
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /logic4/ nfu(1)
      COMMON/SHELL/LSHELLT,LSHELIJ,LSHELKL,LHELP,LCAS2(4),LCAS3(4)
      common /lcases/ lcase
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
ccc   common /big/ bl(1)
      common /memor4/ iwt0,iwt1,iwt2,ibuf,ibuf2,
     * ibfij1,ibfij2,ibfkl1,ibfkl2,
     * ibf2l1,ibf2l2,ibf2l3,ibf2l4,ibfij3,ibfkl3,
     * ibf3l,issss,
     * ix2l1,ix2l2,ix2l3,ix2l4,ix3l1,ix3l2,ix3l3,ix3l4,
     * ixij,iyij,izij, iwij,ivij,iuij,isij
c
      common /memor4a/ ibf3l1,ibf3l2,ibf3l3,ibf3l4
      common /memor5a/ iaa,ibb,icc,idd,icis,icjs,icks,icls,
     * ixab,ixp,ixpn,ixpp,iabnia,iapb,i1apb,ifij,icij,isab,
     * ixcd,ixq,ixqn,ixqq,icdnia,icpd,i1cpd,ifkl,ickl,iscd
      common /memor5b/ irppq,
     * irho,irr1,irys,irhoapb,irhocpd,iconst,ixwp,ixwq,ip1234,
     * idx1,idx2,indx
c
c dimensions for assembling :
      common /dimasse/ lqij,lqkl,lqmx,lij3,lkl3,l3l,lsss
      dimension bl(*)
c----------------------------------------------------------------
        nqijr=nqij
        nqklr=nqkl
      if(where.eq.'shif') then
c-    - for nmr derivatives -
        nqijr=nqij1
        nqklr=nqkl1
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 2.or.lcase.eq. 6.or.lcase.eq. 8.or.lcase.eq. 9.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lshelij.eq.1 .or. lshelij.eq.3) then
c-                     --- for bfij1 s from -> lx/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1_2(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *                  bl(ibfij1),lqij,lnkl, bl(icis),lci,
     *                  bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 3.or.lcase.eq. 6.or.lcase.eq.10.or.lcase.eq.11.or.
c    *   lcase.eq.12.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelij.eq.2 .or. lshelij.eq.3) then
c-                        --- for bfij2 s from xl/yz ---
c
                ijenx=nfu(nqijr+1)
                if(nqijr.eq.nsij) then
                   ijenx=1
                   if(where.eq.'shif') ijenx=nfu(nqij+1)
                endif
                klenx=lnkl
        firstx=firstc
        call conijkl1_2(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *                  bl(ibfij2),lqij,lnkl, bl(icjs),lcj,
     *                  bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 4.or.lcase.eq. 7.or.lcase.eq. 8.or.lcase.eq.10.or.
c    *   lcase.eq.12.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.1 .or. lshelkl.eq.3) then
c-                          --- for bfkl1 s from xy/lz ---
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1_2(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *                  bl(ibfkl1),lnij,lqkl, bl(icks),lck,
     *                  bl(indx), ijenx,klenx)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 5.or.lcase.eq. 7.or.lcase.eq. 9.or.lcase.eq.11.or.
c    *   lcase.eq.13.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.2 .or. lshelkl.eq.3) then
c-                          --- for bfkl2 s from xy/zl ---
c
                ijenx=lnij
                klenx=nfu(nqklr+1)
                if(nqklr.eq.nskl) then
                   klenx=1
                   if(where.eq.'shif') klenx=nfu(nqkl+1)
                endif
        firstx=firstc
        call conijkl1_2(bl,firstx,nbls,nbls1,bl(iwt0), l01,l02,
     *                  bl(ibfkl2),lnij,lqkl,  bl(icls),lcl,
     *                  bl(indx), ijenx,klenx)
      endif
c
      IF(lshellt.eq.1) go to 100
c-----------------------------------------------------------------------
c     if(lcase.eq. 6.or.lcase.eq.12.or.lcase.eq.13.or.lcase.eq.16) then
c-
      if(lshelij.eq.3) then
c-                          --- for bfij3  ss from ll/xy ---
c
                ij3b=1
                kl3b=klbeg
        firstx=firstc
        call conijkl3_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                  bl(ibfij3),lij3,lnkl, bl(ifij),lcij,
     *                  bl(indx),ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 7.or.lcase.eq.14.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lshelkl.eq.3) then
c-                         --- for bfkl3 ss from xy/ll ---
c
                ij3b=ijbeg
                kl3b=1
        firstx=firstc
        call conijkl3_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                  bl(ibfkl3),lnij,lkl3, bl(ifkl),lckl,
     *                  bl(indx), ij3b,kl3b)
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 8.or.lcase.eq.12.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(1).eq.1) then
c-                          --- for bf2l1  ss from lx/ly ---
c
        firstx=firstc
        call conb2ln_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibf2l1),lqij,lqkl,
     *                 bl(icis),bl(icks),lci,lck,
     *                 bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq. 9.or.lcase.eq.13.or.lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas2(2).eq.1) then
c-                          --- for bf2l2  ss from lx/yl ---
c
        firstx=firstc
        call conb2ln_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibf2l2),lqij,lqkl,
     *                 bl(icis),bl(icls),lci,lcl,
     *                 bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.10.or.lcase.eq.12.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(3).eq.1) then
c-                          --- for bf2l3  ss from xl/ly ---
c
        firstx=firstc
        call conb2ln_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibf2l3),lqij,lqkl,
     *                 bl(icjs),bl(icks),lcj,lck,
     *                 bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.11.or.lcase.eq.13.or.lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas2(4).eq.1) then
c-                          --- for bf2l4  ss from xl/yl ---
c
        firstx=firstc
        call conb2ln_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ibf2l4),lqij,lqkl,
     *                 bl(icjs),bl(icls),lcj,lcl,
     *                 bl(indx))
      endif
c
      IF(lshellt.eq.2) go to 100
c-----------------------------------------------------------------------
c-                      --- for bf3l  ---
c
c     if(lcase.eq.12.or.lcase.eq.16) then
c-
      if(lcas3(1).eq.1) then
c-            --- for bf3l1  sss from ll/lx ---
c
        firstx=firstc
        call conb3la_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(icks),bl(ifij),  bl(ibf3l1),l3l,lqmx,
     *                 lck,lcij,bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.13.or.lcase.eq.16) then
c-
      if(lcas3(2).eq.1) then
c-            --- for bf3l2  sss from ll/xl ---
c
        firstx=firstc
        call conb3la_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(icls),bl(ifij),  bl(ibf3l2),l3l,lqmx,
     *                 lcl,lcij,bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.14.or.lcase.eq.16) then
c-
      if(lcas3(3).eq.1) then
c-            --- for bf3l3  sss from lx/ll ---
c
        firstx=firstc
        call conb3lb_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(icis),bl(ifkl),  bl(ibf3l3),lqmx,l3l,
     *                 lci,lckl,bl(indx))
      endif
c-----------------------------------------------------------------------
c     if(lcase.eq.15.or.lcase.eq.16) then
c-
      if(lcas3(4).eq.1) then
c-            --- for bf3l4  sss from xl/ll ---
c
        firstx=firstc
        call conb3lb_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(icjs),bl(ifkl),  bl(ibf3l4),lqmx,l3l,
     *                 lcj,lckl,bl(indx))
      endif
c
      IF(lshellt.eq.3) go to 100
c-----------------------------------------------------------------------
      if(lcase.eq.16) then
c-    -- for ssss(nbls)  ssss from ll/ll --
c
        firstx=firstc
        call conssss_2(bl,firstx,nbls,nbls1,bl(iwt0),l01,l02,
     *                 bl(ifij),bl(ifkl), bl(issss),lsss ,
     *                 lcij,lckl,bl(indx))
      endif
c
c-----------------------------------------------------------------------
c
  100 CONTINUE
c
      firstc=firstx
c
      end
c===============================================================
      subroutine conijkl1_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                      bfij1,lt3,lt4, facti,lci,
     *                      indx, ijenx,klenx)
      implicit real*8 (a-h,o-z)
      logical firstc
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
      dimension facti(*)
      dimension indx(*)
c-------------------------------------------
      call assel1_2(firstc,xt1,lt1,lt2,nbls,
     *           facti(lci),bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
c
      end
c===============================================================
      subroutine assel1_2(firstc,xt1,lt1,lt2,nbls,
     *            facti, bfij1,lt3,lt4, indx,nbls1, ijenx,klenx)
      implicit real*8 (a-h,o-z)
      logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij1(nbls,lt3,lt4)
cccc  dimension facti(*)
c--------------------------------------------
C
              facti1=facti
      IF (FIRSTC) THEN
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
              do 504 i=1,nbls1
              ijkl=indx(i)
              bfij1(ijkl,ij,kl)=xt1(i,ij,kl)*facti1
  504         CONTINUE
              FIRSTC=.FALSE.
      ELSE
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
              do 604 i=1,nbls1
              ijkl=indx(i)
              bfij1(ijkl,ij,kl)=bfij1(ijkl,ij,kl)+xt1(i,ij,kl)*facti1
  604         CONTINUE
c
      ENDIF
      end
c===============================================================
      subroutine conijkl3_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                      bfij3,lt3,lt4,  factij,lcij,
     *                      indx, ij3b,kl3b)
      implicit real*8 (a-h,o-z)
      logical firstc
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
      dimension factij(*)
      dimension indx(*)
c----------------------------------------------------------
      call assel2a_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *         bfij3,lt3,lt4,factij(lcij), indx, ij3b,kl3b)
c
      end
c================================================================
      subroutine assel2a_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                    bfij3,lt3,lt4, factij, indx, ij3b,kl3b)
      implicit real*8 (a-h,o-z)
      logical firstc
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bfij3(nbls,lt3,lt4)
ccccc dimension factij(*)
c----------------------------------------------
              factij1=factij
      IF (FIRSTC) THEN
              do 502 kl=kl3b,lt4
              do 502 ij=ij3b,lt3
              do 502 i=1,nbls1
                 ijkl=indx(i)
                 bfij3(ijkl,ij,kl)=xt1(i,ij,kl)*factij1
  502         CONTINUE
              FIRSTC=.FALSE.
      ELSE
              DO 602 KL=kl3b,lt4
              do 602 ij=ij3b,lt3
              do 602 i=1,nbls1
              ijkl=indx(i)
              bfij3(ijkl,ij,kl)=bfij3(ijkl,ij,kl)+xt1(i,ij,kl)*factij1
  602         CONTINUE
      ENDIF
      end
c===============================================================
      subroutine conb2ln_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     bf2l1, lt3,lt4,
     *                     facti,factk,
     *                     lci,lck,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
      dimension facti(*),factk(*)
      dimension indx(*)
c----------------------------------------------
      call assel2b_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             bf2l1,lt3,lt4,facti(lci), factk(lck) ,indx)
c
      return
      end
c================================================================
      subroutine assel2b_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   bf2l1,lt3,lt4,facti,factk,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf2l1(nbls,lt3,lt4)
cccc  dimension facti(*),factk(*)
c--------------------------------------------------------
             if(where.eq.'shif') then
c            -- for nmr derivatives --
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
                klenx=nfu(nqkl1+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c----------------------------------
              factik1=facti*factk
      IF (FIRSTC) THEN
              DO 504 KL=KLBEG,klenx
              DO 504 IJ=IJBEG,ijenx
                 do 504 i=1,nbls1
                 ijkl=indx(i)
                 bf2l1(ijkl,ij,kl)=xt1(i,ij,kl)*factik1
  504         CONTINUE
              FIRSTC=.FALSE.
      ELSE
c-----------> DO 604 KL=KLBEG,NFU(NQKL+1)
              DO 604 KL=KLBEG,klenx
              DO 604 IJ=IJBEG,ijenx
                 do 604 i=1,nbls1
                 ijkl=indx(i)
              bf2l1(ijkl,ij,kl)=bf2l1(ijkl,ij,kl)+xt1(i,ij,kl)*factik1
  604         CONTINUE
      ENDIF
      end
c===============================================================
      subroutine conb3la_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     factk,factij, bf3l,lt5,lt6,
     *                     lck,lcij,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
cccc  common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension factk(*), factij(*)
      dimension indx(*)
c-------------------------------------------------------
      call assel3a_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             factk(lck),factij(lcij), bf3l,lt5,lt6,  indx)
c
      end
c===========================================================
      subroutine assel3a_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   factk,factij, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
cccc  dimension factk(*),factij(*)
c----------------------------------------------------------
             if(where.eq.'shif') then
c            -- for nmr derivatives --
cccccc          ijenx=nfu(nqij1+1)
c--------nie->  if(nqij1.eq.nsij) ijenx=1
                klenx=nfu(nqkl1+1)
c WHAT ??       if(nqkl1.eq.nskl) klenx=nfu(nqij+1)
                if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
cccccc          ijenx=nfu(nqij+1)
c--------nie->  if(nqij.eq.nsij) ijenx=1
                klenx=nfu(nqkl+1)
                if(nqkl.eq.nskl) klenx=1
             endif
c------------------------------------------
             factijk=factij*factk
      IF (FIRSTC) THEN
              DO 503 KL=KLBEG,klenx
              do 503 ij=1,lt5
              do 503 i=1,nbls1
                 ijkl=indx(i)
                 bf3l(ijkl,ij,kl)=xt1(i,ij,kl)*factijk
  503         CONTINUE
              FIRSTC=.FALSE.
      ELSE
            DO 603 KL=KLBEG,klenx
            do 603 ij=1,lt5
               do 603 i=1,nbls1
               ijkl=indx(i)
               bf3l(ijkl,ij,kl)=bf3l(ijkl,ij,kl)+xt1(i,ij,kl)*factijk
  603       CONTINUE
      ENDIF
      end
c===============================================================
      subroutine conb3lb_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     facti,factkl, bf3l,lt5,lt6,
     *                     lci,lckl,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension bf3l(nbls,lt5,lt6)
      dimension facti(*),factkl(*)
      dimension indx(*)
c--------------------------------------------------
      call assel3b_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *             facti(lci),factkl(lckl), bf3l,lt5,lt6, indx)
c
      end
c===============================================================
      subroutine assel3b_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                   facti,factkl, bf3l,lt5,lt6, indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2), bf3l(nbls,lt5,lt6)
cccc  dimension facti(*), factkl(*)
c------------------------------------------------------
             if(where.eq.'shif') then
c            -- for nmr derivatives --
                ijenx=nfu(nqij1+1)
                if(nqij1.eq.nsij) ijenx=nfu(nqij+1)
cccccc          klenx=nfu(nqkl1+1)
c--------nie--> if(nqkl1.eq.nskl) klenx=nfu(nqkl+1)
             else
                ijenx=nfu(nqij+1)
                if(nqij.eq.nsij) ijenx=1
cccccc          klenx=nfu(nqkl+1)
c--------nie--> if(nqkl.eq.nskl) klenx=1
             endif
c------------------------------------------------------
             factkli=factkl*facti
      IF (FIRSTC) THEN
              DO 511 kl=1,lt6
              DO 511 IJ=IJBEG,ijenx
              do 511 i=1,nbls1
                 ijkl=indx(i)
                 bf3l(ijkl,ij,kl)=xt1(i,ij,kl)*factkli
  511         CONTINUE
              FIRSTC=.FALSE.
      ELSE
              DO 611 kl=1,lt6
              DO 611 IJ=IJBEG,ijenx
              do 611 i=1,nbls1
              ijkl=indx(i)
                 bf3l(ijkl,ij,kl)=bf3l(ijkl,ij,kl)+xt1(i,ij,kl)*factkli
  611         CONTINUE
      ENDIF
      end
c===============================================================
      subroutine conssss_2(bl,firstc,nbls,nbls1,xt1,lt1,lt2,
     *                     factij,factkl, ssss,isdim,
     *                     lcij,lckl,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
ccc   common /big/ bl(1)
      dimension bl(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
      dimension factij(*),factkl(*)
      dimension indx(*)
c-----------------------------------------------------------
      call assel4_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *            factij(lcij),factkl(lckl), ssss,isdim, indx)
c
      end
c===============================================================
      subroutine assel4_2(firstc,nbls,nbls1,xt1,lt1,lt2,
     *                  factij,factkl, ssss,isdim,indx)
      implicit real*8 (a-h,o-z)
      logical firstc
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension ssss(nbls,isdim,isdim)
cccc  dimension factij(*),factkl(*)
c------------------------------------------------
              fijkl=factij*factkl
      IF (FIRSTC) THEN
              do 507 kl=1,isdim
              do 507 ij=1,isdim
              do 507 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)= xt1(i,ij,kl)*fijkl
  507         continue
              FIRSTC=.FALSE.
      ELSE
              do 607 kl=1,isdim
              do 607 ij=1,isdim
              do 607 i=1,nbls1
              ijkl=indx(i)
              ssss(ijkl,ij,kl)=ssss(ijkl,ij,kl)+xt1(i,ij,kl)*fijkl
  607         continue
      ENDIF
c
      end
c=======================================================================
      subroutine der1_ne22(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                     indx,aax,lci,bbx,lcj,ccx,lck)
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=2 , nbls>nbls1
c     firstc=.false. ( not first contr. step)
c-----------------------------------------------------------------------
c three regions should be distinguished here :
c 1:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(1,ijkl,ij,kl)=xt1(i,ij,kl)        <-- ordinary integrals
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx <---int.rescaled by 2*b_exp
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 2:
c    do kl=nfu(nskl)+1,nfu(nskl+1)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 3:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nsij)+1,nfu(nsij+1)
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx <---int.rescaled by 2*b_exp
c    enddo
c    enddo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
ccccc logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(*),bbx(*),ccx(*)
      dimension buf2(nbls,lt1,lt2,4)
      data zero /0.d0/
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
      aax1=aax(lci)
      bbx1=bbx(lcj)
      ccx1=ccx(lck)
c
cccc  IF (FIRSTC) THEN
cccc  ELSE
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,1)=buf2(ijkl,ij,kl,1)+xt1(i,ij,kl)
               buf2(ijkl,ij,kl,2)=buf2(ijkl,ij,kl,2)+xt1(i,ij,kl)*aax1
               buf2(ijkl,ij,kl,3)=buf2(ijkl,ij,kl,3)+xt1(i,ij,kl)*bbx1
               buf2(ijkl,ij,kl,4)=buf2(ijkl,ij,kl,4)+xt1(i,ij,kl)*ccx1
               enddo
            enddo
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,2)=buf2(ijkl,ij,kl,2)+xt1(i,ij,kl)*aax1
               buf2(ijkl,ij,kl,3)=buf2(ijkl,ij,kl,3)+xt1(i,ij,kl)*bbx1
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               ijkl=indx(i)
               buf2(ijkl,ij,kl,4)=buf2(ijkl,ij,kl,4)+xt1(i,ij,kl)*ccx1
               enddo
            enddo
         enddo
cccc  ENDIF
c
      end
c=======================================================================
      subroutine der1_ne21(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                           idxnot,nbls2,
     *                           indx,aax,lci,bbx,lcj,ccx,lck)
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=2 , nbls>nbls1
c     firstc=.true (first contr. step)
c-----------------------------------------------------------------------
c three regions should be distinguished here :
c 1:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(1,ijkl,ij,kl)=xt1(i,ij,kl)        <-- ordinary integrals
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx <---int.rescaled by 2*b_exp
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 2:
c    do kl=nfu(nskl)+1,nfu(nskl+1)
c    do ij=nfu(nqij)+1,nfu(nsij)
c       buf2(4,ijkl,ij,kl)=xt1(i,ij,kl)*ccx <---integ.rescaled by 2*c_exp
c    enddo
c    enddo
c 3:
c    do kl=nfu(nqkl)+1,nfu(nskl)
c    do ij=nfu(nsij)+1,nfu(nsij+1)
c       buf2(2,ijkl,ij,kl)=xt1(i,ij,kl)*aax <---int.rescaled by 2*a_exp
c       buf2(3,ijkl,ij,kl)=xt1(i,ij,kl)*bbx <---int.rescaled by 2*b_exp
c    enddo
c    enddo
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
cccc  logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension idxnot(*)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(*),bbx(*),ccx(*)
      dimension buf2(nbls,lt1,lt2,4)
      data zero /0.d0/
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
      aax1=aax(lci)
      bbx1=bbx(lcj)
      ccx1=ccx(lck)
c
cccc  IF (FIRSTC) THEN
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,1)=xt1(i,ij,kl)
                  buf2(ijkl,ij,kl,2)=xt1(i,ij,kl)*aax1
                  buf2(ijkl,ij,kl,3)=xt1(i,ij,kl)*bbx1
                  buf2(ijkl,ij,kl,4)=xt1(i,ij,kl)*ccx1
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,1)=zero
                  buf2(ijkl,ij,kl,2)=zero
                  buf2(ijkl,ij,kl,3)=zero
                  buf2(ijkl,ij,kl,4)=zero
               enddo
            enddo
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,2)=xt1(i,ij,kl)*aax1
                  buf2(ijkl,ij,kl,3)=xt1(i,ij,kl)*bbx1
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,2)=zero
                  buf2(ijkl,ij,kl,3)=zero
               enddo
            enddo
         enddo
c
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
                  ijkl=indx(i)
                  buf2(ijkl,ij,kl,4)=xt1(i,ij,kl)*ccx1
               enddo
               do i=1,nbls2
                  ijkl=idxnot(i)
                  buf2(ijkl,ij,kl,4)=zero
               enddo
            enddo
         enddo
cccc     FIRSTC=.FALSE.
cccc  ELSE
cccc  ENDIF
c
      end
c=======================================================================
      subroutine der2_2(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *                          idxnot,nblsnot,
     *                          aax,lci,bbx,lcj,ccx,lck,expo)
c
cccc  this is called only for where.eq.'hess'
c
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension idxnot(*)
      dimension xt1(nbls1,lt1,lt2)
c
      dimension aax(*),bbx(*),ccx(*),expo(6)
C
c2001 dimension buf2(10,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,10)
c
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c end for first derivatives
c               buf2(5,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*b_exp
c               buf2(6,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*c_exp
c               buf2(7,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*c_exp
c               buf2(8,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*a_exp
c               buf2(9,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*b_exp
c               buf2(10nbls,lt1,lt2) - rescaled with 2*c_exp * 2*c_exp
c end for second derivatives
c-------
      data zero /0.d0/
c-------
c multiply exponents :
c
          aax1=aax(lci)
          bbx1=bbx(lcj)
          ccx1=ccx(lck)
          expo(1)=aax1*bbx1
          expo(2)=aax1*ccx1
          expo(3)=bbx1*ccx1
          expo(4)=aax1*aax1
          expo(5)=bbx1*bbx1
          expo(6)=ccx1*ccx1
c------------------------------------------------
        ij1=nfu(nqij)
        ij2=nfu(nsij-1)
        ij3=nfu(nsij)
        ij4=nfu(nsij+1)
c
        kl1=nfu(nqkl)
        kl2=nfu(nskl-1)
        kl3=nfu(nskl)
        kl4=nfu(nskl+1)
c------------------------------------------------
        IF (FIRSTC) THEN
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,1)=XT1(i,IJ,KL)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax1
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx1
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx1
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax1
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx1
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx1
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           FIRSTC=.FALSE.
        ELSE
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,1)=buf2(ijkl,ij,kl,1)+XT1(i,IJ,KL)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax1
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx1
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx1
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax1
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx1
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx1
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
           enddo
        ENDIF
c------------------------------------------------
      end
c=======================================================================
      subroutine der1_eq21(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                     aax,lci,bbx,lcj,ccx,lck)
c
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=2 , nbls=nbls1
c     firstc=.true.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(*),bbx(*),ccx(*)
      dimension buf2(nbls,lt1,lt2,4)
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
      aax1=aax(lci)
      bbx1=bbx(lcj)
      ccx1=ccx(lck)
c
ctrue IF (FIRSTC) THEN
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,1)=xt1(i,ij,kl)
               buf2(i,ij,kl,2)=xt1(i,ij,kl)*aax1
               buf2(i,ij,kl,3)=xt1(i,ij,kl)*bbx1
               buf2(i,ij,kl,4)=xt1(i,ij,kl)*ccx1
               enddo
            enddo
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               buf2(i,ij,kl,2)=xt1(i,ij,kl)*aax1
               buf2(i,ij,kl,3)=xt1(i,ij,kl)*bbx1
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,4)=xt1(i,ij,kl)*ccx1
               enddo
            enddo
         enddo
ctrue    FIRSTC=.FALSE.
ctrue ELSE
ctrue ENDIF
c
      end
c=======================================================================
      subroutine der1_eq22(nbls,nbls1,xt1,lt1,lt2,buf2,
     *                     aax,lci,bbx,lcj,ccx,lck)
c
c-----------------------------------------------------------------------
c     this is called only for where.eq.'forc'
c     route=2 , nbls=nbls1
c     firstc=.false.
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension xt1(nbls1,lt1,lt2)
      dimension aax(*),bbx(*),ccx(*)
      dimension buf2(nbls,lt1,lt2,4)
      data zero /0.d0/
c-----------------------------------------------------------
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c-----------------------------------------------------------
      aax1=aax(lci)
      bbx1=bbx(lcj)
      ccx1=ccx(lck)
c
         do kl=nfu(nqkl)+1,nfu(nskl)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,1)=buf2(i,ij,kl,1) + xt1(i,ij,kl)
               buf2(i,ij,kl,2)=buf2(i,ij,kl,2) + xt1(i,ij,kl)*aax1
               buf2(i,ij,kl,3)=buf2(i,ij,kl,3) + xt1(i,ij,kl)*bbx1
               buf2(i,ij,kl,4)=buf2(i,ij,kl,4) + xt1(i,ij,kl)*ccx1
               enddo
            enddo
            do ij=nfu(nsij)+1,nfu(nsij+1)
               do i=1,nbls1
               buf2(i,ij,kl,2)=buf2(i,ij,kl,2) + xt1(i,ij,kl)*aax1
               buf2(i,ij,kl,3)=buf2(i,ij,kl,3) + xt1(i,ij,kl)*bbx1
               enddo
            enddo
         enddo
         do kl=nfu(nskl)+1,nfu(nskl+1)
            do ij=nfu(nqij)+1,nfu(nsij)
               do i=1,nbls1
               buf2(i,ij,kl,4)=buf2(i,ij,kl,4) + xt1(i,ij,kl)*ccx1
               enddo
            enddo
         enddo
c
      end
c=======================================================================
c For general contracted hessian : 2002
c=======================================================================
      subroutine goef_x_exp1(nbls,nbls1,indx,ngcd, gcoef,aax,bbx,ccx,
     *                       gcexp )              ! output
      implicit real*8 (a-h,o-z)
      dimension indx(*)
      dimension gcoef(ngcd,nbls)
      dimension aax(nbls1),bbx(nbls1),ccx(nbls1)
      dimension gcexp(ngcd,10,nbls)
c
      do i=1,nbls1
         axi=aax(i)
         bxi=bbx(i)
         cxi=ccx(i)
c
         abxi=axi*bxi
         acxi=axi*cxi
         bcxi=bxi*cxi
c
         aaxi=axi*axi
         bbxi=bxi*bxi
         ccxi=cxi*cxi
c
         ijkl=indx(i)
         do iqu=1,ngcd
            gcexp(iqu, 1,ijkl)=gcoef(iqu,ijkl)
            gcexp(iqu, 2,ijkl)=gcoef(iqu,ijkl)*axi
            gcexp(iqu, 3,ijkl)=gcoef(iqu,ijkl)*bxi
            gcexp(iqu, 4,ijkl)=gcoef(iqu,ijkl)*cxi
            gcexp(iqu, 5,ijkl)=gcoef(iqu,ijkl)*abxi
            gcexp(iqu, 6,ijkl)=gcoef(iqu,ijkl)*acxi
            gcexp(iqu, 7,ijkl)=gcoef(iqu,ijkl)*bcxi
            gcexp(iqu, 8,ijkl)=gcoef(iqu,ijkl)*aaxi
            gcexp(iqu, 9,ijkl)=gcoef(iqu,ijkl)*bbxi
            gcexp(iqu,10,ijkl)=gcoef(iqu,ijkl)*ccxi
         enddo
      enddo
c
      end
c=======================================================================
      subroutine goef_x_exp2(nbls,nbls1,indx,ngcd,
     *                      gcoef,aax,bbx,ccx,lci,lcj,lck,
     *                      gcexp )
c2000 used only for iroute=2
      implicit real*8 (a-h,o-z)
      dimension indx(*)
      dimension gcoef(ngcd,nbls)
      dimension aax(*),bbx(*),ccx(*)
      dimension gcexp(ngcd,10,nbls)
c
         axi=aax(lci)
         bxi=bbx(lcj)
         cxi=ccx(lck)
c
         abxi=axi*bxi
         acxi=axi*cxi
         bcxi=bxi*cxi
c
         aaxi=axi*axi
         bbxi=bxi*bxi
         ccxi=cxi*cxi
c
      do i=1,nbls1
         ijkl=indx(i)
         do iqu=1,ngcd
            gcexp(iqu, 1,ijkl)=gcoef(iqu,ijkl)
            gcexp(iqu, 2,ijkl)=gcoef(iqu,ijkl)*axi
            gcexp(iqu, 3,ijkl)=gcoef(iqu,ijkl)*bxi
            gcexp(iqu, 4,ijkl)=gcoef(iqu,ijkl)*cxi
            gcexp(iqu, 5,ijkl)=gcoef(iqu,ijkl)*abxi
            gcexp(iqu, 6,ijkl)=gcoef(iqu,ijkl)*acxi
            gcexp(iqu, 7,ijkl)=gcoef(iqu,ijkl)*bcxi
            gcexp(iqu, 8,ijkl)=gcoef(iqu,ijkl)*aaxi
            gcexp(iqu, 9,ijkl)=gcoef(iqu,ijkl)*bbxi
            gcexp(iqu,10,ijkl)=gcoef(iqu,ijkl)*ccxi
         enddo
      enddo
c
      end
c=======================================================================
      subroutine asselg_der2(firstc,xt1,lt1,lt2,nbls,but2,indx,nbls1,
     *                       ngcd,gcexp )
      implicit real*8 (a-h,o-z)
      logical firstc
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
      common /logic4/ nfu(1)
      dimension indx(*)
      dimension xt1(nbls1,lt1,lt2)
      dimension but2(ngcd,10,nbls,lt1,lt2)
      dimension gcexp(ngcd,10,nbls)
      data zero /0.d0/
c-----------------------------------------------------------------------
c               but2(1,nbls,lt1,lt2) - ordinary contraction
c               but2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               but2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               but2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c end for first derivatives
c               but2(5,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*b_exp
c               but2(6,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*c_exp
c               but2(7,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*c_exp
c               but2(8,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*a_exp
c               but2(9,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*b_exp
c               but2(10nbls,lt1,lt2) - rescaled with 2*c_exp * 2*c_exp
c end for second derivatives
c-----------------------------------------------------------------------
        IF (firstc) THEN
           DO 551 kl=nfu(nqkl)+1,nfu(nskl)
           DO 551 ij=nfu(nqij)+1,nfu(nsij)
           do 551 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
              do ist=1,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
  551      CONTINUE
c
           DO 552 kl=nfu(nqkl)+1,nfu(nskl)
           DO 552 ij=nfu(nsij)+1,nfu(nsij+1)
           do 552 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
              do ist=2,3
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
              do ist=5,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
c...........  ist=1 & 4 NOT needed
  552      CONTINUE
c
           DO 553 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 553 ij=nfu(nqij)+1,nfu(nsij)
           do 553 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
c...........  ist=1,2,3 NOT needed
              do ist=4,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
  553      CONTINUE
           FIRSTC=.FALSE.
        ELSE
           DO 651 KL=nfu(nqkl)+1,nfu(nskl)
           DO 651 ij=nfu(nqij)+1,nfu(nsij)
           do 651 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
              do ist=1,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=but2(iqu,ist,ijkl,IJ,KL)
     *                                      +xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
  651      CONTINUE
c
           DO 652 KL=nfu(nqkl)+1,nfu(nskl)
           DO 652 ij=nfu(nsij)+1,nfu(nsij+1)
           do 652 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
              do ist=2,3
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=but2(iqu,ist,ijkl,IJ,KL)
     *                                      +xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
              do ist=5,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=but2(iqu,ist,ijkl,IJ,KL)
     *                                      +xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
c...........  ist=1 & 4 NOT needed
  652      CONTINUE
c
           DO 653 KL=nfu(nskl)+1,nfu(nskl+1)
           DO 653 ij=nfu(nqij)+1,nfu(nsij)
           do 653 i=1,nbls1
              ijkl=indx(i)
              xint=XT1(i,IJ,KL)
c...........  ist=1,2,3 NOT needed
              do ist=4,10
                 do iqu=1,ngcd
                    but2(iqu,ist,ijkl,IJ,KL)=but2(iqu,ist,ijkl,IJ,KL)
     *                                      +xint*gcexp(iqu,ist,ijkl)
                 enddo
              enddo
  653      CONTINUE
        ENDIF
c
      end
c=======================================================================
c=======================================================================
      subroutine der2_2driv(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *                          idxnot,nblsnot,
     *                          aax,lci,bbx,lcj,ccx,lck,expo)
c
cccc  this is called only for where.eq.'hess'
c
      implicit real*8 (a-h,o-z)
      logical firstc, firstx
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension idxnot(*)
      dimension xt1(nbls1,lt1,lt2)
c
      dimension aax(*),bbx(*),ccx(*),expo(6)
C
      dimension buf2(nbls,lt1,lt2,10)
c
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c end for first derivatives
c               buf2(5,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*b_exp
c               buf2(6,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*c_exp
c               buf2(7,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*c_exp
c               buf2(8,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*a_exp
c               buf2(9,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*b_exp
c               buf2(10nbls,lt1,lt2) - rescaled with 2*c_exp * 2*c_exp
c end for second derivatives
c
      data zero /0.d0/
c------------------------------------------------
c multiply exponents :
c
          aax1=aax(lci)
          bbx1=bbx(lcj)
          ccx1=ccx(lck)
          expo(1)=aax1*bbx1
          expo(2)=aax1*ccx1
          expo(3)=bbx1*ccx1
          expo(4)=aax1*aax1
          expo(5)=bbx1*bbx1
          expo(6)=ccx1*ccx1
c------------------------------------------------
c first assemble zero- and first-order derivatives :
c
      firstx=firstc
      call der2_21(firstx,nbls,nbls1,xt1,lt1,lt2,buf2(1,1,1,1),indx,
     *             idxnot,nblsnot,aax1,bbx1,ccx1)
c------------------------------------------------
c then assemble second-order derivatives :
c
      firstx=firstc
      call der2_22(firstx,nbls,nbls1,xt1,lt1,lt2,buf2(1,1,1,5),indx,
     *             idxnot,nblsnot,expo)
c------------------------------------------------
      firstc=firstx
c------------------------------------------------
      end
c=======================================================================
      subroutine der2_21(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *                          idxnot,nblsnot,aax1,bbx1,ccx1)
c
cccc  this is called only for where.eq.'hess'
c
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension idxnot(*)
      dimension xt1(nbls1,lt1,lt2)
c
      dimension buf2(nbls,lt1,lt2,4)
c
c               buf2(1,nbls,lt1,lt2) - ordinary contraction
c               buf2(2,nbls,lt1,lt2) - rescaled with 2*a_exp
c               buf2(3,nbls,lt1,lt2) - rescaled with 2*b_exp
c               buf2(4,nbls,lt1,lt2) - rescaled with 2*c_exp
c end for first derivatives
c
      data zero /0.d0/
c------------------------------------------------
        ij1=nfu(nqij)
        ij2=nfu(nsij-1)
        ij3=nfu(nsij)
        ij4=nfu(nsij+1)
c
        kl1=nfu(nqkl)
        kl2=nfu(nskl-1)
        kl3=nfu(nskl)
        kl4=nfu(nskl+1)
c------------------------------------------------
        IF (FIRSTC) THEN
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,1)=XT1(i,IJ,KL)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax1
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx1
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx1
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,2)=XT1(i,IJ,KL)*aax1
                 buf2(ijkl,IJ,KL,3)=XT1(i,IJ,KL)*bbx1
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,4)=XT1(i,IJ,KL)*ccx1
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=1,4
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           FIRSTC=.FALSE.
        ELSE
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,1)=buf2(ijkl,ij,kl,1)+XT1(i,IJ,KL)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax1
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx1
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx1
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,2)=buf2(ijkl,ij,kl,2)+XT1(i,IJ,KL)*aax1
           buf2(ijkl,IJ,KL,3)=buf2(ijkl,ij,kl,3)+XT1(i,IJ,KL)*bbx1
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,4)=buf2(ijkl,ij,kl,4)+XT1(i,IJ,KL)*ccx1
                 enddo
              enddo
           enddo
        ENDIF
c------------------------------------------------
      end
c=======================================================================
      subroutine der2_22(firstc,nbls,nbls1,xt1,lt1,lt2,buf2,indx,
     *                          idxnot,nblsnot,expo)
c
cccc  this is called only for where.eq.'hess'
c
      implicit real*8 (a-h,o-z)
      logical firstc
c
      common/obarai/
     * lni,lnj,lnk,lnl,lnij,lnkl,lnijkl,MMAX,
     * NQI,NQJ,NQK,NQL,NSIJ,NSKL,
     * NQIJ,NQIJ1,NSIJ1,NQKL,NQKL1,NSKL1,ijbeg,klbeg
C
      common /logic4/ nfu(1)
c
      dimension indx(*)
      dimension idxnot(*)
      dimension xt1(nbls1,lt1,lt2)
c
      dimension expo(6)
C
c2001 dimension buf2(10,nbls,lt1,lt2)
      dimension buf2(nbls,lt1,lt2,5:10)
c
c               buf2(5,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*b_exp
c               buf2(6,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*c_exp
c               buf2(7,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*c_exp
c               buf2(8,nbls,lt1,lt2) - rescaled with 2*a_exp * 2*a_exp
c               buf2(9,nbls,lt1,lt2) - rescaled with 2*b_exp * 2*b_exp
c               buf2(10nbls,lt1,lt2) - rescaled with 2*c_exp * 2*c_exp
c end for second derivatives
c
      data zero /0.d0/
c------------------------------------------------
        ij1=nfu(nqij)
        ij2=nfu(nsij-1)
        ij3=nfu(nsij)
        ij4=nfu(nsij+1)
c
        kl1=nfu(nqkl)
        kl2=nfu(nskl-1)
        kl3=nfu(nskl)
        kl4=nfu(nskl+1)
c------------------------------------------------
        IF (FIRSTC) THEN
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,5)=XT1(i,IJ,KL)*expo(1)
                 buf2(ijkl,IJ,KL,8)=XT1(i,IJ,KL)*expo(4)
                 buf2(ijkl,IJ,KL,9)=XT1(i,IJ,KL)*expo(5)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,6)=XT1(i,IJ,KL)*expo(2)
                 buf2(ijkl,IJ,KL,7)=XT1(i,IJ,KL)*expo(3)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
                 buf2(ijkl,IJ,KL,10)=XT1(i,IJ,KL)*expo(6)
                 enddo
                 do i=1,nblsnot
                 ijkl=idxnot(i)
                    do ii=5,10
                    buf2(ijkl,IJ,KL,ii)= zero
                    enddo
                 enddo
              enddo
           enddo
           FIRSTC=.FALSE.
        ELSE
           do kl=kl1+1,kl2
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
                 enddo
              enddo
              do ij=ij3+1,ij4
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,5)=buf2(ijkl,ij,kl,5)+XT1(i,IJ,KL)*expo(1)
           buf2(ijkl,IJ,KL,8)=buf2(ijkl,ij,kl,8)+XT1(i,IJ,KL)*expo(4)
           buf2(ijkl,IJ,KL,9)=buf2(ijkl,ij,kl,9)+XT1(i,IJ,KL)*expo(5)
                 enddo
              enddo
           enddo
           do kl=kl2+1,kl3
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
              do ij=ij2+1,ij3
                 do i=1,nbls1
                 ijkl=indx(i)
           buf2(ijkl,IJ,KL,6)=buf2(ijkl,ij,kl,6)+XT1(i,IJ,KL)*expo(2)
           buf2(ijkl,IJ,KL,7)=buf2(ijkl,ij,kl,7)+XT1(i,IJ,KL)*expo(3)
                 enddo
              enddo
           enddo
           do kl=kl3+1,kl4
              do ij=ij1+1,ij2
                 do i=1,nbls1
                 ijkl=indx(i)
          buf2(ijkl,IJ,KL,10)=buf2(ijkl,ij,kl,10)+XT1(i,IJ,KL)*expo(6)
                 enddo
              enddo
           enddo
        ENDIF
c------------------------------------------------
      end
c=======================================================================

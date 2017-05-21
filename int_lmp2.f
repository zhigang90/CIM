      subroutine int_lmp2(bl,    inx,   thres1,
     *                    icsh,icf1,icf2,  kcsh,kcf1,kcf2,
     *                    mapf2s,DS,    iprnt, xmp2int, nintotal,
     *                    nrow,  ncol,  irow,  icol,      lzero)
c---------------------------------------------------------------------
c This routine calculates the integrals for MP2, (ICSH,jcs|KCSH,lcs) for

c a given  ICSH & KCSH  and all jcs,lcs and puts them into xmp2int(*)

c INPUT:
c  bl          = common bl (free memory)
c  inx         = common inx (contraction info)
c  thres1      = integral threshold
c  icsh,kcsh   = fixed contracted shells
c  mapf2s(*)   = an array mapping from contr. functions to contr.shells
c  DS          = density screening matrix
c  iprnt       = if negative, integrals are printed
c
c OUTPUT :
c xmp2int      = an array holding "lmp2" integrals from (ICS,jcs|KCS,lcs)
c                these integrals are NOT re-scaled to real values
c  nintotal    = total number of integrals generated
c  nrow        = number of non-zero rows
c  ncol        = number of non-zero columns
c  irow        = the first nrow elements hold the indices of non-zero rows
c  icol        =  same as above for columns
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /memmax/ ispblx, maxme1,iforwhat
c---------------------------------------------------------------------
      common /mp2shells/ ics_mp2, kcs_mp2      ! used here and in calcint2.f
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      common /memor1c/ map_ij_bl2
      dimension inx(12,*)
      dimension bl(*),xmp2int(*),DS(*)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
c-----------------------------------------------------------------
c for use in cshneg.f (shows type of screening)
c
      nintotal=0
      screen='mp2 '
c-----------------------------------------------------------------
c zero out irow and icol. These arrays store ultimately the
c indices of rows and columns of xmp2int which contain at least one
c non-zero element
c
      call izeroit(irow,ncf)
      call izeroit(icol,ncf)
      call izeroit(lzero,ncf)
c----------------------------------------------------------------
      ics_mp2=icsh
      kcs_mp2=kcsh
c
      call mmark
c----------------------------------------------------------------
c make blocks of contracted shell quartets for requested Ics,Kcs:
c
      iforwhat=5
      call make_blocks4ik(bl,inx,icsh,kcsh,nbl4,maxlabels,ds)
      call getint(maxlabels,labels)
c----------------------------------------------------------------
c calculate requested MP2 integrals
c
      call do_mp2_int(bl, inx, thres1,xmp2int,
     *               icsh,icf1,icf2,  kcsh,kcf1,kcf2,
     *               bl(labels),ncf,   ncs,    nbl4,
     *               mapf2s,DS,nintotal,irow,icol,lzero,iprnt)
c
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------
c Now compact irow and icol
      nrow=0
      ncol=0
      do i=1,ncf
        if(irow(i).gt.0) then
          nrow=nrow+1
          irow(nrow)=i
        end if
        if(icol(i).gt.0) then
          ncol=ncol+1
          icol(ncol)=i
        end if
      end do
c----------------------------------------------------------------
c
      end
c=====================================================================
      subroutine do_mp2_int(bl,    inx,      thres1,  xmp2int,
     *                      icsh,icf1,icf2,  kcsh,kcf1,kcf2,
     *                      labels,   ncf,     ncs,     nbl4,
     *                      mapf2s,DS,nintotal,irow,icol,lzero,iprnt)
c----------------------------------------------------------------
c  Krzysztof needs to put in a description of the parameters for do_mp2_int
c---------------------------------------------------------------------
c This routine calculates the integrals for MP2, (ICSH,jcs|KCSH,lcs) for
c a given  ICSH & KCSH  and all jcs,lcs and puts them into xmp2int(*)
c INPUT:
c  bl          = common bl (free memory)
c  inx         = common inx (contraction info)
c  thres1      = integral threshold
c  icsh,kcsh   = fixed contracted shells
c  labels      = array for integral's indeces
c  ncf         = number of contracted basis functions
c  ncs         = number of contracted basis shells
c  nbl4        = number of blocks of contracted basis shell quartets
c  mapf2s(*)   = an array mapping from contr. functions to contr.shells
c  DS          = density screening matrix
c  iprnt       = if negative, integrals are printed
c
c OUTPUT :
c xmp2int      = an array holding "lmp2" integrals from (ICS,jcs|KCS,lcs)
c                these integrals are NOT re-scaled to real values
c  nintotal    = total number of integrals generated
c  nrow        = number of non-zero rows
c  ncol        = number of non-zero columns
c  irow        = the first nrow elements hold the indices of non-zero rows
c  icol        =  same as above for columns
c---------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
c----------------------------------------------------------------
      common /screen_type/ screen
      common /memmax/ ispblx, maxme1,iforwhat
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),xmp2int(*)
      dimension mapf2s(*),DS(*)
      dimension labels(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
c----------------------------------------------------------------
c calculate "LMP2" integrals (ICSH,jcs|KCSH,lcs)
c for a given  ICSH & KCSH  and all jcs,lcs
c-----------------------------------------------------------------
cc      idensp=mataddr('dsmx')
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c
c
      DO isupbl=1,nbl4
c
c        find max element of "screening" density for the ISUPBL block
c
         call dmax_find_4b(isupbl,bl,DS,ncs,screen)
c
   22    continue
c
         where='lmp2' ! must be set up here because it is
c                       changed to 'fock' in calcint2
c
         call calcint2_new(isupbl,bl,inx,thres1,DS,where,
     *                     labels,ibuffz,nblsiz,nintez,ngctoz,
     *                     moreint,stopnow)
c
c        check if requested super-block was in a list :
         if(stopnow) go to 2222
c
CKW.....................................................
c        nintez is the size of a given quartet.It is set up
c        to zero if there are no integrals
CKW.....................................................
c
c        integrals arrived in bl(ibuffz); check if they are :
         if(nintez.gt.0) then
              nintotal=nintotal+nintez*ngctoz*nblsiz
              call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c             put "lmp2" integrals into xmp2int array
c
           if(ngctoz.eq.1) then   !     segmented basis set
            if(icsh.eq.kcsh) then
              call lmp2Ebldx(ncf,  bl(ibuffz), nblsiz, ngctoz, nintez,
     *          labels(lab1),labels(lab2),labels(lab3),mapf2s, icsh,
     *                        kcsh, irow, icol, xmp2int,
     *                        icf1,icf2, kcf1,kcf2, lzero)
            else
              call lmp2Nbldx(ncf,  bl(ibuffz), nblsiz, ngctoz, nintez,
     *          labels(lab1),labels(lab2),labels(lab3),mapf2s, icsh,
     *                        kcsh, irow, icol, xmp2int,
     *                        icf1,icf2, kcf1,kcf2, lzero)
            endif
           else                    !    general contracted basis set
              call lmp2_bldx(ncf,  bl(ibuffz), nblsiz, ngctoz, nintez,
     *          labels(lab1),labels(lab2),labels(lab3),mapf2s, icsh,
     *                        kcsh, irow, icol, xmp2int,
     *                        icf1,icf2, kcf1,kcf2, lzero)
           endif
         endif
c----------------------------------------------------------------
         if(moreint) go to 22
 2222    continue
      ENDDO
c
      if(iprnt.lt.0) then
      call lmp2_print(icsh,kcsh,xmp2int,icf1,icf2,kcf1,kcf2,ncf,thres1)
      endif
c
c----------------------------------------------------------------
      end
c==============================================================
      subroutine get_shell_size(inx,ics,ics_size)
c     called only for lmp2 integ.
      dimension inx(12,*)
c
      ics_size=inx(3,ics)*(inx(4,ics)+1)
c
      end
c==============================================================
      subroutine lmp2_print(icsh,kcsh,xmp2int,icf1,icf2,kcf1,kcf2,ncf,
     *                      thres1)
      implicit real*8 (a-h,o-z)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
      data sum_tot /0.d0/
      save sum_tot
c
c     write(6,*)'------------------------------------- '
      write(6,*)' lmp2 integrals for shell =',icsh,kcsh
c     write(6,*)'------------------------------------- '
c
      sum=0.d0
      do icf=icf1,icf2
      do kcf=kcf1,kcf2
          write(6,*)' X mat for icf,kcf=',icf,kcf
      do jcf=1,ncf
      do lcf=1,ncf
         xint=xmp2int(lcf,jcf,kcf,icf)*thres1
c        write(6,55) icf,jcf,kcf,lcf,xint
c55      format(4i3,2x,f12.8)
         write(6,44) jcf,lcf,xint
 44      format(2i3,2x,f12.8)
         sum=sum+xint
      enddo
      enddo
      enddo
      enddo
c
      sum_tot=sum_tot + sum
c
         write(6,66) icsh,kcsh,sum,sum_tot
 66      format
     * (' MP2 shells icsh,kcsh=',2i3,'  sum=',f12.5,'  sum_tot=',f15.5)
c
c     write(6,*)'------------------------------------- '
c
      end
c=====================================================================
      subroutine lmp2_bldr(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contaction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
             icase=0
c------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 icase=1
ccc              write(6,*)' case=1 : Ij Kl'
                 intct=0
                 do 201 icf=icff1,icff2
                 do 201 jcf=jcff1,jcff2
CSS
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
CSS
                 do 201 kcf=kcff1,kcff2
                 do 201 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(jcf)=1
  201            continue
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 icase=2
ccc              write(6,*)' case=2 : Ij lK'
                 intct=0
                 do 202 icf=icff1,icff2
                 do 202 jcf=jcff1,jcff2
CSS
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
CSS
                 do 202 kcf=kcff1,kcff2
                 do 202 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(jcf)=1
  202            continue
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 icase=3
ccc              write(6,*)' case=3 : Ji Kl'
                 intct=0
                 do 203 icf=icff1,icff2
CSS
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
CSS
                 do 203 jcf=jcff1,jcff2
                 do 203 kcf=kcff1,kcff2
                 do 203 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(icf)=1
  203            continue
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 icase=4
ccc              write(6,*)' case=4 : Ji lK'
                 intct=0
                 do 204 icf=icff1,icff2
CSS
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
CSS
                 do 204 jcf=jcff1,jcff2
                 do 204 kcf=kcff1,kcff2
                 do 204 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(icf)=1
  204            continue
             endif
c-------------------------------------------------------
c cases with switched IJbl2 & KLbl2 :
c
c (5) (Kl|Ij) :
             if( (kcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 icase=5
ccc              write(6,*)' case=5 : Kl Ij '
                 intct=0
                 do 205 icf=icff1,icff2
                 do 205 jcf=jcff1,jcff2
                 do 205 kcf=kcff1,kcff2
                 do 205 lcf=lcff1,lcff2
CSS
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
CSS
                    intct=intct+1
                    xmp2int(jcf,lcf,icf,kcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(lcf)=1
  205            continue
             endif
c-------------------------------------------------------
c (6) (lK|Ij) :
             if( (lcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 icase=6
ccc              write(6,*)' case=6 : lK Ij '
                 intct=0
                 do 206 icf=icff1,icff2
                 do 206 jcf=jcff1,jcff2
                 do 206 kcf=kcff1,kcff2
CSS
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
CSS
                 do 206 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(jcf,kcf,icf,lcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(kcf)=1
  206            continue
             endif
c-------------------------------------------------------
c (7) (Kl|jI) :
             if( (kcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 icase=7
ccc              write(6,*)' case=7 : Kl jI '
                 intct=0
                 do 207 icf=icff1,icff2
                 do 207 jcf=jcff1,jcff2
                 do 207 kcf=kcff1,kcff2
                 do 207 lcf=lcff1,lcff2
CSS
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
CSS
                    intct=intct+1
                    xmp2int(icf,lcf,jcf,kcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(lcf)=1
  207            continue
             endif
c-------------------------------------------------------
c (8) (lK|jI) :
             if( (lcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 icase=8
ccc              write(6,*)' case=8 : Kl jI '
                 intct=0
                 do 208 icf=icff1,icff2
                 do 208 jcf=jcff1,jcff2
                 do 208 kcf=kcff1,kcff2
CSS
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
CSS
                 do 208 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(icf,kcf,jcf,lcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(kcf)=1
  208            continue
             endif
c------------------------------------------------------
             if(icase.EQ.0) then
c                write(6,*)' ICASE=0 '
             endif
c------------------------------------------------------
  150      continue
  100   continue
c
      end
c==============================================================
      subroutine make_blocks4ik(bl,inx,ics_r,kcs_r,nbl4,maxlabels,ds)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /memmax/ ispblx, maxme1,iforwhat
      common /infob/ inuc ,ibasx,nax,nbfx,nshx,ncfx,ncsx
      common /intbuf/ maxibuf,maxindx
c
      dimension ds(*)   ! screening "density" for mp2&mp2d integrals
      dimension bl(*),inx(12,*)
      data timetot /0.d0/
      save timetot
c---------------------------------------------------------------
c     call getmem(0,last_0)
c     call retmem(1)
c
      call secund(ttest1)
c---------------------------------------------------------------
c get basis set data : up to 4 different basis sets :
c
      call get_basis_data(inx_1,ibas_1,ncs_1,ncf_1,nsh_1,nbf_1,
     *                    inx_2,ibas_2,ncs_2,ncf_2,nsh_2,nbf_2,
     *                    inx_3,ibas_3,ncs_3,ncf_3,nsh_3,nbf_3,
     *                    inx_4,ibas_4,ncs_4,ncf_4,nsh_4,nbf_4)
c
      call getival('list_1',list_1)
      call getival('list_2',list_2)
      call getival('list_3',list_3)
      call getival('list_4',list_4)
c---------------------------------------------------------------
c
c FOR LEFT side pairs in (ij|kl) :
c
c first set of shells : whole basis set no 1
c placed in basis 1 locations
c
            ibasis=1
c
            call getint(ncs_1*2,nblock1_1)
            call getint(ncs_1  ,nbback1_1)
c
            ncs_b=ics_r
            ncs_e=ics_r
            call blocking1(iforwhat,bl(inuc),ibasis,bl(list_1),
     *                    ncs_1,bl(inx_1),bl(ibas_1), ncs_b,ncs_e,
     *                    bl(nblock1_1),bl(nbback1_1),nbl1_1,maxshell_1)
c
c second set of shells : the whole basis set no 1
c placed in basis 2 locations
c
            ibasis=2
c
            call getint(ncs_2*2,nblock1_2)
            call getint(ncs_2  ,nbback1_2)
c
            ncs_b=1
            ncs_e=ncs_2
            call blocking1(iforwhat,bl(inuc),ibasis,bl(list_2),
     *                    ncs_2,bl(inx_2),bl(ibas_2), ncs_b,ncs_e,
     *                    bl(nblock1_2),bl(nbback1_2),nbl1_2,maxshell_2)
c
c make pair-blocks for left side of (ij|kl):
c
      call blocking2_same(bl,'ijpairs',inx,
     *                                 ds,
     *                                 nbl1_1,bl(nblock1_1),maxshell_1,
     *                                 nbl1_2,bl(nblock1_2),maxshell_2)
c
c    (there are 3 memory allocations in blocking2 )
c
c FOR RIGHT side pairs in (ij|kl) :
c
c first set of shells : whole basis set no 1
c placed in basis 3 locations
c
            ibasis=3
c
            call getint(ncs_3*2,nblock1_3)
            call getint(ncs_3  ,nbback1_3)
c
            ncs_b=kcs_r
            ncs_e=kcs_r
            call blocking1(iforwhat,bl(inuc),ibasis,bl(list_3),
     *                    ncs_3,bl(inx_3),bl(ibas_3), ncs_b,ncs_e,
     *                    bl(nblock1_3),bl(nbback1_3),nbl1_3,maxshell_3)
c
c second set of shells : the whole basis set no 1
c placed in basis 4 locations
c
            ibasis=4
c
            call getint(ncs_4*2,nblock1_4)
            call getint(ncs_4  ,nbback1_4)
c
            ncs_b=1
            ncs_e=ncs_4
            call blocking1(iforwhat,bl(inuc),ibasis,bl(list_4),
     *                    ncs_4,bl(inx_4),bl(ibas_4), ncs_b,ncs_e,
     *                    bl(nblock1_4),bl(nbback1_4),nbl1_4,maxshell_4)
c
c make pair-blocks for right side of (ij|kl):
c
      call blocking2_same(bl,'klpairs',inx,
     *                                 ds,
     *                                 nbl1_3,bl(nblock1_3),maxshell_3,
     *                                 nbl1_4,bl(nblock1_4),maxshell_4)
c
c    (there are 3 memory allocations in blocking2 )
c----------------------------------------------------------
      call secund(ttest2)
c----------------------------------------------------------
c     call getmem(0,last_1)
c     call retmem(1)
c     write(91,*)' Memory checking :last_0 & last_1=',last_0,last_1
c----------------------------------------------------------
c get data concerning pair-blocks :
c
         call getival('nparx',npar_ij)
         call getival('ijblx',ijbl)
c....... call getival('mapijblx',map_ij_bl2_ij)
         call getival('blocksij',nbl2_ij)
         call getival('blpredij',nbl2_ijd)
c
         call getival('npary',npar_kl)
         call getival('ijbly',klbl)
c....... call getival('mapijbly',map_ij_bl2_kl)
         call getival('blockskl',nbl2_kl)
         call getival('blpredkl',nbl2_kld)
c----------------------------------------------------------
         call getival('lcore',lcore)
         call getival('iroute',iroute)
c----------------------------------------------------------
c NO SPLITTING AT ALL
c
c allocate int. memory for an array showing which pair-blocks
c make a given quartet-block :  nblock4(2,nbl4)
c ijbl2=nblock4(1,ibl4)
c klbl2=nblock4(2,ibl4)
c
      call getint(2*nbl2_ij*nbl2_kl, nblock4)
      call setival('nblock4',nblock4)
c
      call blksizer_new(
     * bl(ijbl),nbl2_ijd,nbl2_ij,bl(npar_ij),bl(inx_1),bl(inx_2),
     * bl(klbl),nbl2_kld,nbl2_kl,bl(npar_kl),bl(inx_3),bl(inx_4),
     *                        bl(nblock4), nbl4)
c
c output :
c nbl4 - number of quartet-blocks
c nblock4(2,*) array
c common /intbuf/ maxibuf,maxindx

c----------------------------------------------------------

c size of the current space for labels :

c

       maxlabels=maxindx

c----------------------------------------------------------
      call secund(ttest3)
      timetot=timetot+(ttest3-ttest1)
      call setrval('timeblks',timetot)
c
c     time12=(ttest2-ttest1)/60.d0
c     time124=(ttest3-ttest1)/60.d0
c     timetot=timetot+time124
c     write(91,9988) time12,time124,timetot
c9988 format(' cpu time for new blocking : 12 =',f6.2,' 124 =',f6.2,
c    *       ' total_cum =',f6.2)
c----------------------------------------------------------
      end
c=======================================================
      subroutine zerocol(xx,ncf,kcf1,kcf2,icf1,icf2,j)
C      zeros out a column j in matrix xx
C      replaces call to zeroit for entire xmp2int
      implicit real*8(a-h,o-z)
      dimension xx(ncf,ncf,kcf1:kcf2,icf1:icf2)
      do kcf=kcf1,kcf2
      do icf=icf1,icf2
      call zeroit(xx(1,j,kcf,icf),ncf)
      enddo
      enddo
      end
C=====================================================================
C	NEW int_lmp2 for  large calculations
C	Svein Saebo Nov. 2000
C=====================================================================
      subroutine int_lmp2b(bl, inx, thres1, icsh, kcsh,
     1                     mapf2s, DS, iprnt, Z, nintotal,
     2                     ic2fc, ic2fr, if2cc, if2cr, indlj,
     3                     xmp2int, nrow, ncol, ind, intstore,
     4                     maxstore, indmax)
c---------------------------------------------------------------------
c This routine calculates the integrals for MP2, (ICSH,jcs|KCSH,lcs) for
c a given  ICSH & KCSH  and all jcs,lcs and puts them into xmp2int(*)
c INPUT:
c  bl          = common bl (free memory)
c  inx         = common inx (contraction info)
c  thres1      = integral threshold
c  icsh,kcsh   = fixed contracted shells
c  mapf2s(*)   = an array mapping from contr. functions to contr.shells
c  DS          = density screening matrix
c  iprnt       = if negative, integrals are printed
c
c OUTPUT :
c xmp2int      = an array holding "lmp2" integrals from (ICS,jcs|KCS,lcs)
c                these integrals are NOT re-scaled to real values
c  nintotal    = total number of integrals generated
c  nrow        = number of non-zero rows
c  ncol        = number of non-zero columns
c  irow        = the first nrow elements hold the indices of non-zero rows
c  icol        =  same as above for columns
c---------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
      common /screen_type/ screen
      common /memmax/ ispblx, maxme1,iforwhat
c---------------------------------------------------------------------
      common /mp2shells/ ics_mp2, kcs_mp2 ! used here and in calcint2.f
      common /infob/ inuc,ibas,na,nbf,nsh,ncf,ncs
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      common /memor1c/ map_ij_bl2
      dimension inx(12,*)
      dimension bl(*),xmp2int(*),Z(*)
      dimension mapf2s(*),DS(*)
      dimension if2cr(*),if2cc(*),ic2fr(*),ic2fc(*)
      integer*2 indlj(6,*)
C
      icf1=inx(11,icsh)+1
      icf2=inx(10,icsh)
      kcf1=inx(11,kcsh)+1
      kcf2=inx(10,kcsh)
      ilen=icf2-icf1+1
      klen=kcf2-kcf1+1
c-----------------------------------------------------------------
c for use in cshneg.f (shows type of screening)
c
      nintotal=0
      screen='mp2 '
c-----------------------------------------------------------------
      ics_mp2=icsh
      kcs_mp2=kcsh
      call mmark
c----------------------------------------------------------------
c make blocks of contracted shell quartets for requested Ics,Kcs:
c
      iforwhat=5
      call make_blocks4ik(bl,inx,icsh,kcsh,nbl4,maxlabels,ds)
      call getint(maxlabels,labels)
c----------------------------------------------------------------
c calculate requested MP2 integrals
c
      call do_mp2_intb(bl,inx,thres1,xmp2int,icsh,
     *               kcsh,bl(labels),ncf,ncs,nbl4,
     *               mapf2s,DS,nintotal,indlj,ind,
     3               intstore,iprnt)
c
      if(intstore.ge.maxstore.or.ind.ge.indmax)then
        iout=igetival('iout')
        write(iout,*) ' memory problems :'
        write(iout,*) ' Use SMAL option or increase memory'
        write(iout,*) ind,intstore,indmax,maxstore
        call nerror(1,'int_lmp2b','overflow',intstore,maxstore)
      endif
c----------------------------------------------------------------
      call retmark
c
      if(ind.eq.0) return
c----------------------------------------------------------------
C  determine rows(jcf) and columns(lcf)with nonzero elements
C
      do iz=1,ncf
      if2cr(iz)=0
      if2cc(iz)=0
      enddo
      nrow=0
      ncol=0
      do ic=1,ind
      ny =indlj(1,ic)
      isi=indlj(2,ic)
      if(if2cr(ny).eq.0) then
      nrow=nrow+1
      ic2fr(nrow)=ny
      if2cr(ny)=nrow
      endif
      if(if2cc(isi).eq.0) then
      ncol=ncol+1
      ic2fc(ncol)=isi
      if2cc(isi)=ncol
      endif
      enddo
C
C  now store the non-zero integrals into
C  Z(kcf1:kcf2,icf1:icf2,icol=1,ncol,irow=1,nrow)
C
C  xmp2int(intstore) contains all nonzero intgerals
C  xmp2int generated in subroutie storeint
C
      lentot=ilen*klen*nrow*ncol
      call zeroit(Z,lentot)
C
      call restr(xmp2int,Z,ncf,nrow,ncol,
     1           if2cr,if2cc,ind,intstore,kcf1,
     2           kcf2,icf1,icf2,indlj)
c----------------------------------------------------------------
      end
c=====================================================================
      subroutine do_mp2_intb(bl,inx,thres1,xmp2int,icsh,
     *                      kcsh,labels,ncf,ncs,nbl4,
     *                      mapf2s,DS,nintotal,indlj,ind,
     3                      intstore,iprnt)
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical moreint,stopnow
      character*11 scftype
      character*4 where
      character*4 screen
c----------------------------------------------------------------
      common /screen_type/ screen
      common /memmax/ ispblx, maxme1,iforwhat
      common /lindvec/ lind,idensp
      common /runtype/ scftype,where
      dimension inx(12,*)
      dimension bl(*),xmp2int(*),DS(*)
      dimension mapf2s(*)
      dimension labels(*)
      integer*2 indlj(6,*)
c----------------------------------------------------------------
c calculate "LMP2" integrals (ICSH,jcs|KCSH,lcs)
c for a given  ICSH & KCSH  and all jcs,lcs
c-----------------------------------------------------------------
c denspar(ics,jcs) matrix in a quadratic form
cc      idensp=mataddr('dsmx')
c----------------------------------------------------------------
c get dimension for xmp2int array
c
      icf1=inx(11,icsh)+1
      icf2=inx(10,icsh)
      kcf1=inx(11,kcsh)+1
      kcf2=inx(10,kcsh)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c
      DO isupbl=1,nbl4
c
c        find max element of "screening" density for the ISUPBL block
c
         call dmax_find_4b(isupbl,bl,DS,ncs,screen)
c
   22    continue
c
         where='lmp2' ! must be set up here because it is
c                       changed to 'fock' in calcint2
c

         call calcint2_new(isupbl,bl,inx,thres1,DS,where,
     *                     labels,ibuffz,nblsiz,nintez,ngctoz,
     *                     moreint,stopnow)
c
c        check if requested super-block was in a list :
         if(stopnow) go to 2222
c
CKW.....................................................
c        nintez is the size of a given quartet.It is set up
c        to zero if there are no integrals
CKW.....................................................
c
c        integrals arrived in bl(ibuffz); check if they are :
         if(nintez.gt.0) then
              nintotal=nintotal+nintez*ngctoz*nblsiz
              call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c             put "lmp2" integrals into xmp2int array
c
      call lmp2_bldrb(ncf,bl(ibuffz),nblsiz,ngctoz,nintez,
     1         labels(lab1),labels(lab2),labels(lab3),mapf2s,icsh,
     2         kcsh, xmp2int,icf1, icf2, kcf1,
     3         kcf2,indlj,ind,intstore)
         endif
c----------------------------------------------------------------
         if(moreint) go to 22
 2222    continue
      ENDDO
c----------------------------------------------------------------
      if(iprnt.lt.0) then
      call lmp2_print(icsh,kcsh,xmp2int,icf1,icf2,kcf1,kcf2,ncf,thres1)
      endif
c----------------------------------------------------------------
      end
c==============================================================
      subroutine lmp2_bldrb(ncf,buf,nbls,ngcd,lnijkl,
     1                     labels,length,lgenct, mapf2s,icsh,
     2                     kcsh,xmp2int,icf1,icf2,kcf1,
     3                     kcf2,indlj,ind,intstore)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contaction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
c
      dimension xmp2int(*)
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      integer*2 indlj(6,*)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
c test
c                intct=0
c                do 299 icf=icff1,icff2
c                do 299 jcf=jcff1,jcff2
c                do 299 kcf=kcff1,kcff2
c                do 299 lcf=lcff1,lcff2
c                   intct=intct+1
c                   xint=buf(ijklp,intct,iqu)
c
c                   write(6,59) icf,jcf,kcf,lcf,xint
c59                 format(4i3,2x,f15.1)
c test
c 299            continue
c-------------------------------------------------------
      call storeint(nbls,ijklp,icsh,kcsh,icsc,
     1                 jcsc,kcsc,lcsc,icff1,jcff1,
     2                 kcff1,lcff1,icff2,jcff2,kcff2,
     3                 lcff2,buf(1,1,iqu),xmp2int,ind,intstore,
     4                 indlj)
C----------------------------------
  150      continue
  100   continue
c
      end
C===============restr===========================================
      subroutine restr(xmp2int,Z,ncf,nrow,ncol,
     1                 if2cr,if2cc,ind,intstore,kcf1,
     2                 kcf2,icf1,icf2,indlj)
C
C	put integrals for a batch into an 4-dim array Z
C	Z(nrow,ncol,LAM,MY)
C	MY LAM over 1 fixed shell each
C	nrow,ncol all AOs with non-zero contributions
C
C	Svein Saebo, Fayetteville, AR Nov. 2000
C	fixed for general contraction Jan 2001
C
C 	this will reduce storage from
C	ncf**2 * ilen*klen   to 2*ilen*klen * nrow*ncol
C	where nrow in the number of rows with at least one
C	nonzero elements and ncol is the number of columns
C	with at least 1 non-zero elements.
C	nrow, ncol <<< ncf for large systems!
C
      implicit real*8(a-h,o-z)
      dimension xmp2int(*)
      dimension Z(ncol,nrow,kcf1:kcf2,icf1:icf2)
      dimension if2cr(*),if2cc(*)
      integer*2 indlj(6,*)
C
      intst=0
      do ic=1,ind
      ny=indlj(1,ic)
      isi=indlj(2,ic)
      irow1=if2cr(ny)
      icol1=if2cc(isi)
C
C	NOTE for general contraction kcf1 is not the same as kcff1
C	etc  kcff1,kcff2,icff1,icff2 are saved in storeint
C
      kcff1=indlj(3,ic)
      kcff2=indlj(4,ic)
      icff1=indlj(5,ic)
      icff2=indlj(6,ic)
C
      do lam=kcff1,kcff2
      do my=icff1,icff2
      intst=intst+1
      Z(icol1,irow1,lam,my)=xmp2int(intst)
      enddo
      enddo
C
      enddo
      end
C===================storeint================================
      subroutine storeint(nbls,ijklp,ifirst,ithird,icsc,
     1                        jcsc,kcsc,lcsc,icff1,jcff1,
     2                        kcff1,lcff1,icff2,jcff2,kcff2,
     3                        lcff2,buf,store,ind,intstore,
     4                        indlj)
C	Input Integral buffer as returned from calcint2_new
C	entered here as buf(1,1,iqu)
C	Output store(intstore)
C	since most integrals are zero for large systems
C	intstore is much smaller than ncf*ncf, and the memory
C	requirement for handling mp2-itegrals is reduced
C	from  klen*ilen*ncf*ncf
C	to    klen*ilen*intstore
C	indices are kept in indlj
C	ifirst = icsh
C	ithird = kcsh
C	
C	Svein Saebo Fayteeville AR, Nov. 2000
C	Fix for  general contractions Jan 2001
C
C	called from lmp2_bldrb
C
      implicit real*8(a-h,o-z)
      dimension buf(nbls,lcff1:lcff2,kcff1:kcff2,jcff1:jcff2,
     1              icff1:icff2)
      dimension store(*)
      integer*2 indlj(6,*)
C
C	These are the same 8 cases as in the original lmp2_bldr
C
c (1) (Ij|Kl) :
             if( (icsc.EQ.ifirst).and.(kcsc.EQ.ithird) ) then
C	write(6,*) 'case 1'
C	first i, second j, third k, fourth l
        do ll=lcff1,lcff2
          do jj=jcff1,jcff2
          ind=ind+1
          indlj(1,ind)=jj
          indlj(2,ind)=ll
C
C	for general contractions we'll need to save:
C         kcff1,kcff2,icff1,icff2
C	these are needed in subroutine restr
C
      indlj(3,ind)=kcff1
      indlj(4,ind)=kcff2
      indlj(5,ind)=icff1
      indlj(6,ind)=icff2
C
      do kk=kcff1,kcff2
      do ii=icff1,icff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
      endif
C---------------------------------------------------------
c (2) (Ij|lK) ,
             if( (icsc.EQ.ifirst).and.(lcsc.EQ.ithird) ) then
C	write(6,*) 'case 2'
C	first i, second j, third l, fourth k
      do kk=kcff1,kcff2
      do jj=jcff1,jcff2
      ind=ind+1
      indlj(1,ind)=jj
      indlj(2,ind)=kk
      indlj(3,ind)=lcff1
      indlj(4,ind)=lcff2
      indlj(5,ind)=icff1
      indlj(6,ind)=icff2
      do ll=lcff1,lcff2
      do ii=icff1,icff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
      endif
C------------------------------------------------------------
c (3) (jI|Kl) ,
             if( (jcsc.EQ.ifirst).and.(kcsc.EQ.ithird) ) then
C	write(6,*) 'case 3'
C	first j, second i, third k fourth l
      do ll=lcff1,lcff2
      do ii=icff1,icff2
          ind=ind+1
      indlj(1,ind)=ii
      indlj(2,ind)=ll
      indlj(3,ind)=kcff1
      indlj(4,ind)=kcff2
      indlj(5,ind)=jcff1
      indlj(6,ind)=jcff2
      do kk=kcff1,kcff2
      do jj=jcff1,jcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
          endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.ifirst).and.(lcsc.EQ.ithird) ) then
C	write(6,*) 'case 4'
C	first j, second i, third l, fourth k
      do kk=kcff1,kcff2
      do ii=icff1,icff2
          ind=ind+1
      indlj(1,ind)=ii
      indlj(2,ind)=kk
      indlj(3,ind)=lcff1
      indlj(4,ind)=lcff2
      indlj(5,ind)=jcff1
      indlj(6,ind)=jcff2
      do ll=lcff1,lcff2
      do jj=jcff1,jcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
          endif
c-------------------------------------------------------
c cases with switched IJbl2 & KLbl2 ,
c
c (5) (Kl|Ij) ,
          if( (kcsc.EQ.ifirst).and.(icsc.EQ.ithird) ) then
C	write(6,*) 'case 5'
C	first k, second l, third i, fourth j
      do jj=jcff1,jcff2
      do ll=lcff1,lcff2
          ind=ind+1
      indlj(1,ind)=ll
      indlj(2,ind)=jj
      indlj(3,ind)=icff1
      indlj(4,ind)=icff2
      indlj(5,ind)=kcff1
      indlj(6,ind)=kcff2
      do ii=icff1,icff2
      do kk=kcff1,kcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
          endif
c-------------------------------------------------------
c (6) (lK|Ij) ,
          if( (lcsc.EQ.ifirst).and.(icsc.EQ.ithird) ) then
C	write(6,*) 'case 6'
C	first l, second k, third i, fourth j
      do jj=jcff1,jcff2
      do kk=kcff1,kcff2
          ind=ind+1
      indlj(1,ind)=kk
      indlj(2,ind)=jj
      indlj(3,ind)=icff1
      indlj(4,ind)=icff2
      indlj(5,ind)=lcff1
      indlj(6,ind)=lcff2
      do ii=icff1,icff2
      do ll=lcff1,lcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
      endif
c-------------------------------------------------------
c (7) (Kl|jI) ,
             if( (kcsc  .EQ.ifirst).and.(jcsc  .EQ.ithird) ) then
C	write(6,*) 'case 7'
C	first k,second l, third j fourth i
      do ii=icff1,icff2
      do ll=lcff1,lcff2
          ind=ind+1
      indlj(1,ind)=ll
      indlj(2,ind)=ii
      indlj(3,ind)=jcff1
      indlj(4,ind)=jcff2
      indlj(5,ind)=kcff1
      indlj(6,ind)=kcff2
      do jj=jcff1,jcff2
      do kk=kcff1,kcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
      endif
c-------------------------------------------------------
c (8) (lK|jI) ,
             if( (lcsc  .EQ.ifirst).and.(jcsc  .EQ.ithird) ) then
C	write(6,*) 'case 8'
C	first l, second k, third j, fourth i
      do ii=icff1,icff2
      do kk=kcff1,kcff2
          ind=ind+1
      indlj(1,ind)=kk
      indlj(2,ind)=ii
      indlj(3,ind)=jcff1
      indlj(4,ind)=jcff2
      indlj(5,ind)=lcff1
      indlj(6,ind)=lcff2
      do jj=jcff1,jcff2
      do ll=lcff1,lcff2
      intstore=intstore+1
      store(intstore)=buf(ijklp,ll,kk,jj,ii)
      enddo
      enddo
      enddo
      enddo
      endif
C
      end
c=====================================================================
      subroutine lmp2Ebldr(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Used for Icsh=Kcsh cases :
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contraction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
cprint thres1=1.d-9
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=1
cccc             write(6,*)' case=1 : Ij Kl ; case=5 : Kl Ij '
                 intct=0
                 do 201 icf=icff1,icff2
                 do 201 jcf=jcff1,jcff2
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
                    icol(jcf)=1
                 do 201 kcf=kcff1,kcff2
                 do 201 lcf=lcff1,lcff2
                    if(lzero(lcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
                       lzero(lcf)=1
                    endif
                    intct=intct+1
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
ckw                 irow(jcf)=1
                    xmp2int(jcf,lcf,icf,kcf)=buf(ijklp,intct,iqu)
ckw                 icol(jcf)=1
                    irow(lcf)=1
c                   xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,kcf,jcf,lcf
c 66  format('case=',i2,3x,4(i2,1x),3x,f12.7,2x,4(i2,1x))
  201            continue
                 go to 150
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=2
c                write(6,*)' case=2 : Ij lK'
                 intct=0
                 do 202 icf=icff1,icff2
                 do 202 jcf=jcff1,jcff2
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
                    icol(jcf)=1
                 do 202 kcf=kcff1,kcff2
                    if(lzero(kcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
                       lzero(kcf)=1
                    endif
                    irow(kcf)=1
                    icol(kcf)=1
                 do 202 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
ckw                 icol(kcf)=1
ckw                 irow(jcf)=1
                    xmp2int(jcf,kcf,icf,lcf)=buf(ijklp,intct,iqu)
ckw                 icol(jcf)=1
ckw                 irow(kcf)=1
c                   xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,lcf,jcf,kcf
  202            continue
                 go to 150
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=3
c                write(6,*)' case=3 : Ji Kl'
                 intct=0
                 do 203 icf=icff1,icff2
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
                    icol(icf)=1
                 do 203 jcf=jcff1,jcff2
                 do 203 kcf=kcff1,kcff2
                 do 203 lcf=lcff1,lcff2
                    if(lzero(lcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
                       lzero(lcf)=1
                    endif
                    intct=intct+1
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
ckw                 irow(icf)=1
                    xmp2int(icf,lcf,jcf,kcf)=buf(ijklp,intct,iqu)
ckw                 icol(icf)=1
                    irow(lcf)=1
c                      xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,kcf,icf,lcf
  203            continue
                 go to 150
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=4
c                write(6,*)' case=4 : Ji lK'
                 intct=0
                 do 204 icf=icff1,icff2
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
                    icol(icf)=1
                 do 204 jcf=jcff1,jcff2
                 do 204 kcf=kcff1,kcff2
                    if(lzero(kcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
                       lzero(kcf)=1
                    endif
                    irow(kcf)=1
                    icol(kcf)=1
                 do 204 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
ckw                 icol(kcf)=1
ckw                 irow(icf)=1
                    xmp2int(icf,kcf,jcf,lcf)=buf(ijklp,intct,iqu)
ckw                 icol(icf)=1
ckw                 irow(kcf)=1
c                      xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,lcf,icf,kcf
  204            continue
                 go to 150
             endif
c-------------------------------------------------------
  150      continue
  100   continue
c-------------------------------------------------------
      end
c=====================================================================
      subroutine lmp2Nbldr(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contraction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
cprint thres1=1.d-9
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c test
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
c            icase=0
c------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=1
c                write(6,*)' case=1 : Ij Kl'
                 intct=0
                 do 201 icf=icff1,icff2
                 do 201 jcf=jcff1,jcff2
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
                 do 201 kcf=kcff1,kcff2
                 do 201 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,kcf,jcf,lcf
c 66  format('case=',i2,3x,4(i2,1x),3x,f12.7,2x,4(i2,1x))
  201            continue
                 go to 150
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=2
c                write(6,*)' case=2 : Ij lK'
                 intct=0
                 do 202 icf=icff1,icff2
                 do 202 jcf=jcff1,jcff2
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
                 do 202 kcf=kcff1,kcff2
                    icol(kcf)=1
                 do 202 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,lcf,jcf,kcf
  202            continue
                 go to 150
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=3
c                write(6,*)' case=3 : Ji Kl'
                 intct=0
                 do 203 icf=icff1,icff2
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
                 do 203 jcf=jcff1,jcff2
                 do 203 kcf=kcff1,kcff2
                 do 203 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,kcf,icf,lcf
  203            continue
                 go to 150
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=4
c                write(6,*)' case=4 : Ji lK'
                 intct=0
                 do 204 icf=icff1,icff2
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
                 do 204 jcf=jcff1,jcff2
                 do 204 kcf=kcff1,kcff2
                    icol(kcf)=1
                 do 204 lcf=lcff1,lcff2
                    intct=intct+1
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,lcf,icf,kcf
  204            continue
                 go to 150
             endif
c-------------------------------------------------------
  150      continue
  100   continue
c-------------------------------------------------------
      end
c==============================================================
      subroutine lmp2NbldX(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contraction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
cprint thres1=1.d-9
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c test
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
c            icase=0
c------------------------------------------------------
c (1) (Ij|Kl) :
ckw          if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
      if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh).and.(icsc.ge.kcsc) ) then
c                icase=1
c                write(6,*)' case=1 : Ij Kl'
c
ckwkwkw          do 201 icf=icff1,icff2
                 do 201 icf=icf1 ,icf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 201 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
ckwkwkwkw        do 201 kcf=kcff1,kcff2
                 do 201 kcf=kcf1 ,kcf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 201 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
c
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,kcf,jcf,lcf
c 66  format('case=',i2,3x,4(i2,1x),3x,f12.7,2x,4(i2,1x))
  201            continue
                 go to 150
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
ckw          if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
      if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh).and.(icsc.ge.lcsc) ) then
c                icase=2
c                write(6,*)' case=2 : Ij lK'
c
ckwkwkw          do 202 icf=icff1,icff2
                 do 202 icf=icf1 ,icf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 202 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                    if(lzero(jcf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
                       lzero(jcf)=1
                    endif
                    irow(jcf)=1
                 do 202 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                    icol(kcf)=1
ckwkwkw          do 202 lcf=lcff1,lcff2
                 do 202 lcf=kcf1 ,kcf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
c
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,icf,lcf,jcf,kcf
  202            continue
                 go to 150
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
ckw          if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
      if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh).and.(jcsc.ge.kcsc) ) then
c                icase=3
c                write(6,*)' case=3 : Ji Kl'
c
                 do 203 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
ckwkwkw          do 203 jcf=jcff1,jcff2
                 do 203 jcf=icf1 ,icf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
ckwkwkw          do 203 kcf=kcff1,kcff2
                 do 203 kcf=kcf1 ,kcf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 203 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,kcf,icf,lcf
  203            continue
                 go to 150
             endif
c------------------------------------------------------
c (4) (jI|lK) :
ckw          if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
      if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh).and.(jcsc.ge.lcsc) ) then
c                icase=4
c                write(6,*)' case=4 : Ji lK'
c
                 do 204 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
                    if(lzero(icf).eq.0) then
                       call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
                       lzero(icf)=1
                    endif
                    irow(icf)=1
ckwkwkw          do 204 jcf=jcff1,jcff2
                 do 204 jcf=icf1 ,icf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 204 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                    icol(kcf)=1
ckwkwkw          do 204 lcf=lcff1,lcff2
                 do 204 lcf=kcf1 ,kcf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
c                xa=buf(ijklp,intct,iqu)*1.0d-9
c                write(6,66) icase,icf,jcf,kcf,lcf, xa,jcf,lcf,icf,kcf
  204            continue
                 go to 150
             endif
c-------------------------------------------------------
  150      continue
  100   continue
c-------------------------------------------------------
      end
c==============================================================
      subroutine lmp2Ebldx(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 EQ KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contaction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
             icase=0
c------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 do 201 icf=icf1 ,icf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 201 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
                 do 201 kcf=kcf1 ,kcf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 201 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(jcf)=1
  201            continue
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 intct=0
                 do 202 icf=icf1 ,icf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 202 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
                 do 202 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 202 lcf=kcf1 ,kcf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(jcf)=1
  202            continue
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 do 203 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
                 do 203 jcf=icf1 ,icf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 203 kcf=kcf1 ,kcf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 203 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(icf)=1
  203            continue
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 do 204 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
                 do 204 jcf=icf1 ,icf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 204 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 204 lcf=kcf1 ,kcf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(icf)=1
  204            continue
             endif
c-------------------------------------------------------
c cases with switched IJbl2 & KLbl2 :
c
c (5) (Kl|Ij) :
             if( (kcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 do 205 icf=kcf1 ,kcf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 205 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 205 kcf=icf1 ,icf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 205 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
                    xmp2int(jcf,lcf,icf,kcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(lcf)=1
  205            continue
             endif
c-------------------------------------------------------
c (6) (lK|Ij) :
             if( (lcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 do 206 icf=kcf1 ,kcf2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 206 jcf=jcff1,jcff2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 206 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
                 do 206 lcf=icf1 ,icf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(jcf,kcf,icf,lcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(kcf)=1
  206            continue
             endif
c-------------------------------------------------------
c (7) (Kl|jI) :
             if( (kcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 do 207 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 207 jcf=kcf1 ,kcf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 207 kcf=icf1 ,icf2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
                 do 207 lcf=lcff1,lcff2
                 l=lcf-lcff
                 intct=ijkijk+l
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
                    xmp2int(icf,lcf,jcf,kcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(lcf)=1
  207            continue
             endif
c-------------------------------------------------------
c (8) (lK|jI) :
             if( (lcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 do 208 icf=icff1,icff2
                 i=icf-icff
                 ii=(i-1)*jlen
                 do 208 jcf=kcf1 ,kcf2
                 j=jcf-jcff
                 ij=ii+j
                 ijij=(ij-1)*klen
                 do 208 kcf=kcff1,kcff2
                 k=kcf-kcff
                 ijk=ijij+k
                 ijkijk=(ijk-1)*llen
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
                 do 208 lcf=icf1 ,icf2
                 l=lcf-lcff
                 intct=ijkijk+l
                    xmp2int(icf,kcf,jcf,lcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(kcf)=1
  208            continue
             endif
c------------------------------------------------------
  150      continue
  100   continue
c
      end
c=====================================================================
      subroutine lmp2_bldx(ncf,   buf,   nbls,  ngcd,  lnijkl,
     *                     labels,length,lgenct, mapf2s,icsh,
     *                     kcsh  ,irow,  icol,  xmp2int,
     *                     icf1,icf2,  kcf1,  kcf2,lzero)
c----------------------------------------------------------------
c Makes an array with lmp2 integrals for the ICS_mp2 NE KCS_mp2 :
c          xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integrals buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c  OUTPUT:
c  irow   = if row k has a non-zero element then irow(k)=1
c  icol   = same for the columns
c
c  OUTPUT:
c  xmp2int: an array (ncf,ncf,kcf1:kcf2,icf1,icf2)
c  the contraction KCS goes from kcf1 to kcf2
c  the contaction ICS goes from icf1 to icf2
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      logical doIcf,doKcf
      dimension buf(nbls,lnijkl,ngcd)
      dimension xmp2int(ncf,ncf,kcf1:kcf2,icf1:icf2)     ! output
c  indices are arranged as l=1:ncf, j=1:ncf, k=kcf1:kcf2, i=icf1:icf2
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension irow(ncf),icol(ncf),lzero(ncf)
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' lmp2_bldr : icf1,kcf1=',icf1,kcf1
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
c
           do 150 iqu=1,ngcd
              icff=labels(1,iqu,ijklp)
              jcff=labels(2,iqu,ijklp)
              kcff=labels(3,iqu,ijklp)
              lcff=labels(4,iqu,ijklp)
c
              icff1=icff+1
              jcff1=jcff+1
              kcff1=kcff+1
              lcff1=lcff+1
c
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c-------------------------------------------------------
c Four possible cases without swiching ijblock2 & klblock2
c
c (Ij|Kl) :  if( (icff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (Ij|lK) :  if( (icff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c (jI|Kl) :  if( (jcff1.EQ.icf1).and.(kcff1.EQ.kcf1) )
c (jI|lK) :  if( (jcff1.EQ.icf1).and.(lcff1.EQ.kcf1) )
c
c Four possible cases with    swiching ijblock2 & klblock2
c
c (Kl|Ij) :  if( (kcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (lK|Ij) :  if( (lcff1.EQ.icf1).and.(icff1.EQ.kcf1) )
c (Kl|jI) :  if( (kcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c (lK|jI) :  if( (lcff1.EQ.icf1).and.(jcff1.EQ.kcf1) )
c-------------------------------------------------------
             icase=0
c------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 intct=0
                 do 201 icf=icff1,icff2
               doIcf=.false.
               if(icf.ge.icf1 .and. icf.le.icf2) doicf=.true.
                 do 201 jcf=jcff1,jcff2
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
                 do 201 kcf=kcff1,kcff2
               doKcf=.false.
               if(kcf.ge.kcf1 .and. kcf.le.kcf2) dokcf=.true.
                 do 201 lcf=lcff1,lcff2
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(lcf,jcf,kcf,icf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(jcf)=1
               endif
  201            continue
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 intct=0
                 do 202 icf=icff1,icff2
               doIcf=.false.
               if(icf.ge.icf1 .and. icf.le.icf2) doIcf=.true.
                 do 202 jcf=jcff1,jcff2
      if(lzero(jcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,jcf)
        lzero(jcf)=1
      endif
                 do 202 kcf=kcff1,kcff2
                 do 202 lcf=lcff1,lcff2
               doKcf=.false.
               if(lcf.ge.kcf1 .and. lcf.le.kcf2) doKcf=.true.
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(kcf,jcf,Lcf,icf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(jcf)=1
               endif
  202            continue
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
                 intct=0
                 do 203 icf=icff1,icff2
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
                 do 203 jcf=jcff1,jcff2
               doIcf=.false.
               if(jcf.ge.icf1 .and. jcf.le.icf2) doIcf=.true.
                 do 203 kcf=kcff1,kcff2
               doKcf=.false.
               if(kcf.ge.kcf1 .and. kcf.le.kcf2) doKcf=.true.
                 do 203 lcf=lcff1,lcff2
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(lcf,icf,kcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(lcf)=1
                    irow(icf)=1
               endif
  203            continue
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
                 intct=0
                 do 204 icf=icff1,icff2
      if(lzero(icf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,icf)
        lzero(icf)=1
      endif
                 do 204 jcf=jcff1,jcff2
               doIcf=.false.
               if(jcf.ge.icf1 .and. jcf.le.icf2) doIcf=.true.
                 do 204 kcf=kcff1,kcff2
                 do 204 lcf=lcff1,lcff2
               doKcf=.false.
               if(lcf.ge.kcf1 .and. lcf.le.kcf2) doKcf=.true.
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(kcf,icf,Lcf,Jcf)=buf(ijklp,intct,iqu)
                    icol(kcf)=1
                    irow(icf)=1
               endif
  204            continue
             endif
c-------------------------------------------------------
c cases with switched IJbl2 & KLbl2 :
c
c (5) (Kl|Ij) :
             if( (kcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 intct=0
                 do 205 icf=icff1,icff2
               doKcf=.false.
               if(icf.ge.kcf1 .and. icf.le.kcf2) doKcf=.true.
                 do 205 jcf=jcff1,jcff2
                 do 205 kcf=kcff1,kcff2
               doIcf=.false.
               if(kcf.ge.icf1 .and. kcf.le.icf2) doIcf=.true.
                 do 205 lcf=lcff1,lcff2
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(jcf,lcf,icf,kcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(lcf)=1
               endif
  205            continue
             endif
c-------------------------------------------------------
c (6) (lK|Ij) :
             if( (lcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
                 intct=0
                 do 206 icf=icff1,icff2
               doKcf=.false.
               if(icf.ge.kcf1 .and. icf.le.kcf2) doKcf=.true.
                 do 206 jcf=jcff1,jcff2
                 do 206 kcf=kcff1,kcff2
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
                 do 206 lcf=lcff1,lcff2
               doIcf=.false.
               if(lcf.ge.icf1 .and. lcf.le.icf2) doIcf=.true.
                    intct=intct+1
                    xmp2int(jcf,kcf,icf,lcf)=buf(ijklp,intct,iqu)
                    icol(jcf)=1
                    irow(kcf)=1
  206            continue
             endif
c-------------------------------------------------------
c (7) (Kl|jI) :
             if( (kcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 intct=0
                 do 207 icf=icff1,icff2
                 do 207 jcf=jcff1,jcff2
               doKcf=.false.
               if(jcf.ge.kcf1 .and. jcf.le.kcf2) doKcf=.true.
                 do 207 kcf=kcff1,kcff2
               doIcf=.false.
               if(kcf.ge.icf1 .and. kcf.le.icf2) doIcf=.true.
                 do 207 lcf=lcff1,lcff2
      if(lzero(lcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,lcf)
        lzero(lcf)=1
      endif
                    intct=intct+1
               if(doIcf.and.doKcf) then
                    xmp2int(icf,lcf,jcf,kcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(lcf)=1
               endif
  207            continue
             endif
c-------------------------------------------------------
c (8) (lK|jI) :
             if( (lcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
                 intct=0
                 do 208 icf=icff1,icff2
                 do 208 jcf=jcff1,jcff2
               doKcf=.false.
               if(jcf.ge.kcf1 .and. jcf.le.kcf2) doKcf=.true.
                 do 208 kcf=kcff1,kcff2
      if(lzero(kcf).eq.0) then
        call zerocol(xmp2int,ncf,kcf1,kcf2,icf1,icf2,kcf)
        lzero(kcf)=1
      endif
                 do 208 lcf=lcff1,lcff2
               doIcf=.false.
               if(lcf.ge.icf1 .and. lcf.le.icf2) doIcf=.true.
                    intct=intct+1
                    xmp2int(icf,kcf,jcf,lcf)=buf(ijklp,intct,iqu)
                    icol(icf)=1
                    irow(kcf)=1
  208            continue
             endif
c------------------------------------------------------
  150      continue
  100   continue
c
      end

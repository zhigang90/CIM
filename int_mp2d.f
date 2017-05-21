      subroutine int_mp2d(bl,  inx,   thres1, icsh,     kcsh,
     *                    mapf2s,DS,DT, nintotal,gradv)
c---------------------------------------------------------------------
c This routine calculates the integral DERIVATIVES for MP2 gradient,
c (ICSH,jcs|KCSH,lcs)XYZ
c for a given  ICSH & KCSH  and all jcs,lcs and puts them into xmp2int(*)
c
c INPUT:
c  bl          = common bl (free memory)
c  inx         = common inx (contraction info)
c  thres1      = integral threshold
c  icsh,kcsh   = fixed contracted shells
c  mapf2s(*)   = an array mapping from contr. functions to contr.shells
c  DS          = density screening matrix
c  DT          = density to be used for forces
c
c OUTPUT :
c gradv        = MP2 forces on atoms
c                these integrals are NOT re-scaled to real values
c  nintotal    = total number of integrals generated
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
      dimension bl(*),DS(*),DT(*)
      dimension gradv(3,*)      ! 3,natoms
      dimension mapf2s(*)
c-----------------------------------------------------------------
c for use in cshneg.f (shows type of screening)
c
      nintotal=0
      screen='mp2d'    !  screening for MP2 derivatives (gradient)
c-----------------------------------------------------------------
      ics_mp2=icsh
      kcs_mp2=kcsh
      call mmark
c----------------------------------------------------------------
c make blocks of contracted shell quartets for requested Ics,Kcs:
c
      iforwhat=6
      call make_blocks4ik(bl,inx,icsh,kcsh,nbl4,maxlabels,ds)
      call getint(maxlabels,labels)
cccc  write(91,*)' BLOCKING FOR ICSH,KCSH=',icsh,kcsh,' Nbl4=',nbl4
c----------------------------------------------------------------
c calculate requested MP2 integrals
c
      call do_mp2_der(bl, inx, thres1,icsh,kcsh, bl(labels),
     *                ncf,ncs,nbl4, mapf2s,DS,DT,nintotal,gradv)
c
c----------------------------------------------------------------
      call retmark
c----------------------------------------------------------------

      end
c=====================================================================
      subroutine do_mp2_der(bl,inx,thres1,icsh, kcsh, labels,
     *                     ncf,ncs,  nbl4,mapf2s,DS,DT,nintotal,gradv)
c----------------------------------------------------------------
c This routine calculates the integralderivatives for MP2 gradient,
c (ICSH,jcs|KCSH,lcs)XYZ  for
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
c  DT          = density to be used for forces
c
c OUTPUT :
c                these integrals are NOT re-scaled to real values
c  nintotal    = total number of integrals generated
c---------------------------------------------------------------------

      use memory, block => bl

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
      dimension bl(*),DS(*),DT(*)
      dimension mapf2s(*)
      dimension labels(*)
      dimension gradv(3,*)      ! 3,natoms
c----------------------------------------------------------------
c for a given  ICSH & KCSH  and all jcs,lcs
c-----------------------------------------------------------------
c get dimension for DT density array
c
      icf1=inx(11,icsh)+1
      icf2=inx(10,icsh)
      kcf1=inx(11,kcsh)+1
      kcf2=inx(10,kcsh)
c----------------------------------------------------------------
c     write(6,*)' from int_mp2d : shells requested =',icsh,kcsh,
c    * ' nbl4=',nbl4,' screen=',screen
ccc   call print_ds(ds,ncs)
c----------------------------------------------------------------
c the very first call to calcint2 must be done with :
c
      moreint=.false.
c----------------------------------------------------------------
c allocate memory for mapping array : ncenter(icf)->iatom
c
      call getint(ncf,ncenter)
c
c integer allocation in bl
c----------------------------------------------------------------
c get centers of basis functions :
c
      call get_center(ncs,ncf,inx,bl(ncenter))
c
c----------------------------------------------------------------
c
      DO isupbl=1,nbl4
c
c        find max element of "screening" density for the ISUPBL block
c
         call dmax_find_4b(isupbl,bl,DS,ncs,screen)
cccc     write(6,*)'       block no=',isupbl
c
   22    continue
c
         where='forc' ! must be set up here because it is
c                       changed to 'fock' in calcint2
c
         call calcint2_new(isupbl,bl,inx,thres1,DS,where,
     *                     labels,ibuffz,nblsiz,nintez,ngctoz,
     *                     moreint,stopnow)
c
c        check if requested super-block was in a list :
         if(stopnow) go to 2222
c
c        integrals arrived in bl(ibuffz); check if they are :
         if(nintez.gt.0) then
              nintotal=nintotal+nintez*ngctoz*nblsiz
              call get_lab123(nblsiz,ngctoz,lab1,lab2,lab3)
c
c             put "lmp2" derivative integrals into xmp2int array
c
              call mp2d_bldr(ncf,bl(ibuffz),nblsiz, ngctoz, nintez,
     *                       bl(ncenter),thres1,DT,icf1,icf2, kcf1,kcf2,
     *                       labels(lab1),labels(lab2),labels(lab3),
     *                       mapf2s,icsh, kcsh,gradv)
         endif
c----------------------------------------------------------------
         if(moreint) go to 22
 2222    continue
      ENDDO
c----------------------------------------------------------------
      end
c=====================================================================
      subroutine mp2d_bldr(ncf, buf, nbls, ngcd, lnijkl,
     *                     ncenter,thres1, DT,icf1,icf2, kcf1,kcf2,
     *                     labels,length,lgenct,
     *                     mapf2s, icsh,kcsh,gradv)
c----------------------------------------------------------------
c  input parameters
c  ncf  = number of contracted basis functions
c  buf  = integral derivative buffer
c  nbls = integ. block-size , number of contracted shell quartets
c  ngcd = total length of general contraction
c lnijkl= size of an integral quartet (number of integrals)
c mapf2s()= maps from contracted function to contracted shell
c icsh,kcsh = requested mp2 shells
c
c  OUTPUT:
c----------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      dimension buf(9,nbls,lnijkl,ngcd)
      dimension dt(ncf,ncf,kcf1:kcf2,icf1:icf2)
      dimension gradv(3,*)      ! 3,natoms
c
      dimension labels(4,ngcd,nbls)
      dimension length(4)
      dimension lgenct(nbls)
      dimension mapf2s(*)
      dimension ncenter(*)
c
      dimension xyz(3,4)  ! local
      dimension vec(9)    ! local
c
      parameter(zero=0.0d0)
c----------------------------------------------------------------
      call getival('nprint',nprint)
c----------------------------------------------------------------
ctest call setup_dt(dt,ncf,kcf1,kcf2,icf1,icf2)
c----------------------------------------------------------------
c     write(6,*)
c    *' requested shells : ics=',icsh,'  kcs=',kcsh,' nbls=',nbls
c----------------------------------------------------------------
      ilen=length(1)
      jlen=length(2)
      klen=length(3)
      llen=length(4)
c----------------------------------------------------------------
c     write(6,*)' mp2d_bldr : icf1,kcf1=',icf1,kcf1,' nbls=',nbls
c----------------------------------------------------------------
c  loop over quartets belonging to the block IKBL :
c
        do 100 ijklp=1,nbls
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
              iat=ncenter(icff1)
              jat=ncenter(jcff1)
              kat=ncenter(kcff1)
              lat=ncenter(lcff1)
c-------------------------------------------------------
              icff2=icff+ilen
              jcff2=jcff+jlen
              kcff2=kcff+klen
              lcff2=lcff+llen
c
              icsc=mapf2s(icff1)
              jcsc=mapf2s(jcff1)
              kcsc=mapf2s(kcff1)
              lcsc=mapf2s(lcff1)
c------------------------------------------------------
              permut=1.d0
              if(icsc.eq.jcsc) permut=permut*0.5d0
              if(kcsc.eq.lcsc) permut=permut*0.5d0
              if(icsc.ge.jcsc)then
                 ijcs=icsc*(icsc-1)/2+jcsc
              else
                 ijcs=jcsc*(jcsc-1)/2+icsc
              endif
              if(kcsc.ge.lcsc)then
                 klcs=kcsc*(kcsc-1)/2+lcsc
              else
                 klcs=lcsc*(lcsc-1)/2+kcsc
              endif
              if(ijcs.eq.klcs) permut=permut*0.5d0
c-------------------------------------------------------
ccccccccccc   call zeroit(xyz,12)
              call zeroit(vec,9)
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
c.......................................................
c for print only
c
c
          IF(nprint.ge.5) then
c            put permutation factors in order
c            to compare with scf/grad integrals
              thres1p=thres1*permut
              integ=0
              do 200 icf=icff1,icff2
                 do 250 jcf=jcff1,jcff2
                    do 300 kcf=kcff1,kcff2
                       do 350 lcf=lcff1,lcff2
                          integ=integ+1
c-------------------------------------------------------
                          xa=buf(1,ijklp,integ,iqu)
                          xb=buf(2,ijklp,integ,iqu)
                          xc=buf(3,ijklp,integ,iqu)
c                         xd=-(xa+xb+xc) ! trans. inv.
                          ya=buf(4,ijklp,integ,iqu)
                          yb=buf(5,ijklp,integ,iqu)
                          yc=buf(6,ijklp,integ,iqu)
c                         yd=-(ya+yb+yc) ! trans. inv.
                          za=buf(7,ijklp,integ,iqu)
                          zb=buf(8,ijklp,integ,iqu)
                          zc=buf(9,ijklp,integ,iqu)
c                         zd=-(za+zb+zc) ! trans. inv.
c--------------------------------------------------------
                          xa=thres1p*xa
                          xb=thres1p*xb
                          xc=thres1p*xc
                          xd=-(xa+xb+xc)
c
                          ya=thres1p*ya
                          yb=thres1p*yb
                          yc=thres1p*yc
                          yd=-(ya+yb+yc)
c
                          za=thres1p*za
                          zb=thres1p*zb
                          zc=thres1p*zc
                          zd=-(za+zb+zc)
                          write(6,66) icf,jcf,kcf,lcf, xa,xb,xc,xd
     *                                               , ya,yb,yc,yd
     *                                               , za,zb,zc,zd
  66  format(4(i2,1x),3x,4(f12.7,2x)/15x,4(f12.7,2x)/15x,4(f12.7,2x))
c--------------------------------------------------------
  350                  continue
  300               continue
  250            continue
  200         continue
c
          ENDIF      !      IF(nprint.ge.5) then
c.......................................................
c-------------------------------------------------------
c            icase=0
c-------------------------------------------------------
c (1) (Ij|Kl) :
             if( (icsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=1
cc               write(6,*)' case=1 : Ij Kl'
                 intct=0
                 do 201 icf=icff1,icff2
                 do 201 jcf=jcff1,jcff2
                 do 201 kcf=kcff1,kcff2
                 do 201 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(lcf,jcf,kcf,icf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(lcf,jcf,kcf,icf)=0.d0
  201            continue
             endif
c------------------------------------------------------
c (2) (Ij|lK) :
             if( (icsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=2
cc               write(6,*)' case=2 : Ij lK'
                 intct=0
                 do 202 icf=icff1,icff2
                 do 202 jcf=jcff1,jcff2
                 do 202 kcf=kcff1,kcff2
                 do 202 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(kcf,jcf,Lcf,icf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(kcf,jcf,Lcf,icf)=0.d0
  202            continue
             endif
c------------------------------------------------------
c (3) (jI|Kl) :
             if( (jcsc.EQ.icsh).and.(kcsc.EQ.kcsh) ) then
c                icase=3
cccc             write(6,*)' case=3 : Ji Kl'
                 intct=0
                 do 203 icf=icff1,icff2
                 do 203 jcf=jcff1,jcff2
                 do 203 kcf=kcff1,kcff2
                 do 203 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(lcf,icf,kcf,Jcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(lcf,icf,kcf,Jcf)=0.d0
  203            continue
             endif
c------------------------------------------------------
c (4) (jI|lK) :
             if( (jcsc.EQ.icsh).and.(lcsc.EQ.kcsh) ) then
c                icase=4
cc               write(6,*)' case=4 : Ji lK'
                 intct=0
                 do 204 icf=icff1,icff2
                 do 204 jcf=jcff1,jcff2
                 do 204 kcf=kcff1,kcff2
                 do 204 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(kcf,icf,Lcf,Jcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(kcf,icf,Lcf,Jcf)=0.d0
  204            continue
             endif
c-------------------------------------------------------
c cases with switched IJbl2 & KLbl2 :
c
c (5) (Kl|Ij) :
             if( (kcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
c                icase=5
c                write(6,*)' case=5 : Kl Ij '
                 intct=0
                 do 205 icf=icff1,icff2
                 do 205 jcf=jcff1,jcff2
                 do 205 kcf=kcff1,kcff2
                 do 205 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(jcf,lcf,icf,kcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(jcf,lcf,icf,kcf)=0.d0
  205            continue
             endif
c-------------------------------------------------------
c (6) (lK|Ij) :
             if( (lcsc  .EQ.icsh).and.(icsc  .EQ.kcsh) ) then
c                icase=6
c                write(6,*)' case=6 : lK Ij '
                 intct=0
                 do 206 icf=icff1,icff2
                 do 206 jcf=jcff1,jcff2
                 do 206 kcf=kcff1,kcff2
                 do 206 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(jcf,kcf,icf,lcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(jcf,kcf,icf,lcf)=0.d0
  206            continue
             endif
c-------------------------------------------------------
c (7) (Kl|jI) :
             if( (kcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
c                icase=7
c                write(6,*)' case=7 : Kl jI '
                 intct=0
                 do 207 icf=icff1,icff2
                 do 207 jcf=jcff1,jcff2
                 do 207 kcf=kcff1,kcff2
                 do 207 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(icf,lcf,jcf,kcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(icf,lcf,jcf,kcf)=0.d0
  207            continue
             endif
c-------------------------------------------------------
c (8) (lK|jI) :
             if( (lcsc  .EQ.icsh).and.(jcsc  .EQ.kcsh) ) then
c                icase=8
c                write(6,*)' case=8 : lK jI '
                 intct=0
                 do 208 icf=icff1,icff2
                 do 208 jcf=jcff1,jcff2
                 do 208 kcf=kcff1,kcff2
                 do 208 lcf=lcff1,lcff2
                    intct=intct+1
                    tikjl=dt(icf,kcf,jcf,lcf)
                    call daxpy(9,tikjl,buf(1,ijklp,intct,iqu),1,vec,1)
c
c                   dt(icf,kcf,jcf,lcf)=0.d0
  208            continue
             endif
c------------------------------------------------------
cc    write(6,*)
c    *'   quart=',ijklp,' cent=',iat,jat,kat,lat,' perm=',permut,
c    *'  shells=',icsc,jcsc,kcsc,lcsc
c-------------------------------------------------------
c     write(6,*)
c    *'   quart=',ijklp,' cent=',iat,jat,kat,lat,' perm=',permut
c-------------------------------------------------------
             factor=permut
             call dscal(9,factor,vec,1)
c------------------------------------------------------
             xyz(1,1)=vec(1)
             xyz(1,2)=vec(2)
             xyz(1,3)=vec(3)
             xyz(2,1)=vec(4)
             xyz(2,2)=vec(5)
             xyz(2,3)=vec(6)
             xyz(3,1)=vec(7)
             xyz(3,2)=vec(8)
             xyz(3,3)=vec(9)
c------------------------------------------------------
c            if(icase.EQ.0) then
c                write(6,*)' ICASE=0 '
c            endif
c------------------------------------------------------
c translational invariance :
c
             xyz(1,4)=-xyz(1,1)-xyz(1,2)-xyz(1,3)
             xyz(2,4)=-xyz(2,1)-xyz(2,2)-xyz(2,3)
             xyz(3,4)=-xyz(3,1)-xyz(3,2)-xyz(3,3)
c------------------------------------------------------
             do icr=1,3
                gradv(icr,iat)=gradv(icr,iat)-xyz(icr,1)
                gradv(icr,jat)=gradv(icr,jat)-xyz(icr,2)
                gradv(icr,kat)=gradv(icr,kat)-xyz(icr,3)
                gradv(icr,lat)=gradv(icr,lat)-xyz(icr,4)
             enddo
c--------------------------------------------------------
c            if(nsym.gt.0) then
c               call getival('SymNuPr',nupair)
c               call atforce_symm(na,nsym,inxx(nupair),iat,jat,kat,lat,
c    *                            xyz,gradv)
c            endif
c--------------------------------------------------------
  150      continue
  100   continue
c-------------------------------------------------------
      end
c==============================================================
      subroutine print_ds(ds,ncs)
      implicit real*8 (a-h,o-z)
c
      dimension ds(ncs,ncs)
c
c     write(6,*)' screening density ncs= ',ncs
c     write(6,66) ( (ds(i,j),j=1,ncs),i=1,ncs)
c 66  format(5(f12.6,1x))
c---------------test only---------------------
ccc   do ics=1,ncs
ccc   do jcs=1,ncs
ccc      ds(ics,jcs)=1000000.d0
ccc      ds(ics,jcs)=1000.d0
ccc      ds(ics,jcs)=100.d0
ccc      ds(ics,jcs)=10.d0
ccc      ds(ics,jcs)=1.d0
ccc   enddo
ccc   enddo
c---------------test only---------------------
c
      end
c==============================================================
      subroutine setup_dt(dt,ncf,kcf1,kcf2,icf1,icf2)
      implicit real*8 (a-h,o-z)
      dimension dt(ncf,ncf,kcf1:kcf2,icf1:icf2)
c
      do icf=icf1,icf2
      do kcf=kcf1,kcf2
         do j=1,ncf
         do i=1,ncf
            dt(i,j,kcf,icf)=0.01d0
         enddo
         enddo
      enddo
      enddo
c
      end
c=====================================================================

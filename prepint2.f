c=======================================================================
      subroutine prepint2(bl,eps,inuc,ibas,na,nbf,nsh,ncf,ncs,inx,
     *                    lcore,nsym,lsemi,scftype)

      use memory, block => bl

c---------------------------------------------------------------
c dimensions for logical matrices in commons logicd,logic1-11 :
c
c up to ff,ff :
c     parameter (lpar1=14,lpar2= 455,lpar3= 364,lpar4=5,lpar5=13)
c up to gg,gg :
c     parameter (lpar1=18,lpar2= 969,lpar3= 816,lpar4=6,lpar5=17)
c up to hh,hh :
c     parameter (lpar1=22,lpar2=1771,lpar3=1540,lpar4=7,lpar5=21)
c up to ii,ii :
c     parameter (lpar1=26,lpar2=2925,lpar3=2600,lpar4=8,lpar5=25)
c up to jj,jj :
c     parameter (lpar1=30,lpar2=4495,lpar3=4060,lpar4=9,lpar5=29)
c up to kk,kk (which is ii,ii + 2 )
c
      parameter (lpar1=34,lpar2=6545,lpar3=5984,lpar4=10,lpar5=33)
c
c up to ll,ll :
c     parameter (lpar1=38,lpar2=9139,lpar3=8436,lpar4=11,lpar5=37)
c up to mm,mm :
c     parameter (lpar1=42,lpar2=12341,lpar3=11480,lpar4=12,lpar5=41)
c---------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)
      logical firstd
      logical in_core,on_disk
      character*11 scftype
c     common /intbl/maxsh,inxx(100)
      common /intbuf/ maxibuf,maxindx           ! output from blksizer
      common /memmax/ ispblx, maxme1,iforwhat
      common /route/ iroute
      common /outfile/ ioutput
      common /timing2/ timprep2,timstor2,timreca2,timblok2
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c
      common /logicd/ hnia(3,lpar2)
      common /logic1/ ndege(lpar4)
      common /logic2/ len(lpar4)
      common /logic3/ lensm(lpar5)
c
      common /logic4/ nfu(lpar1)
      common /logic5/ icoor(lpar2)
      common /logic6/ icool(lpar2)
      common /logic7/ ifrst(lpar2)
      common /logic8/ ilast(lpar2)
      common /logic9/ nia(3,lpar2)
      common /logic10/ nmxyz(3,lpar2)
      common /logic11/ npxyz(3,lpar3)
c
      common /lindvec/ lind,idensp
c
      common /cpu/ intsize,iacc,icache,memreal
      common /infor/ icheck,firstd,ndirect,nprint,iblok,nbeg,nend
      common /infob/ inucx,ibasx,nax,nbfx,nshx,ncfx,ncsx
c
      common /memor1/ iisd,jjsd,ijbld
      common /memor1a/ npard,ncost,nsupb, mxsize
      common /memor1b/ nbl2,nbloks
      common /memor2/ nqrtd, nibld,nkbld, nijbd,nijed, nklbd,nkled
c
      common /memors/ nsymx,ijshp,isymm
      common /symmet/ rnsym
c----------------------------------------------------------------
c data for stored integrlas in-core
c
c     common /datacore/ incorex,ilab,jlab,klab,llab,integrals,
c    *                  isiz,jsiz,ksiz,lsiz, nqstore,iqrt,nblstore
c----------------------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
      dimension xintxx(9)
      dimension intnu(15)
      dimension xntnu(15)
      dimension xint_blk(15),cost_blk(15)
c----------------------------------------------------------------
c integrals from blocks more expensive than ...
c xintxx(9) is transfered back
c----------------------------------------------------------------
      call secund(tim1)
c----------------------------------------------------------
c Set up the integral threshold :
c
      call setup_thres(eps)
c----------------------------------------------------------
c zero out timings
c
      timprep2=zero
      timstor2=zero
      timreca2=zero
c----------------------------------------------------------
c
      call fprep
c----------------------------------------------------------
c temporarly : set basis set data :
c
      call set_basis_data(bl,inx,bl(ibas),ncs,ncf,nsh,nbf)
c
c make 4 copies of a basis set for the get_basis_data
c----------------------------------------------------------
c get basis set data : up to 4 different basis sets :
c
      call get_basis_data(inx_1,ibas_1,ncs_1,ncf_1,nsh_1,nbf_1,
     *                    inx_2,ibas_2,ncs_2,ncf_2,nsh_2,nbf_2,
     *                    inx_3,ibas_3,ncs_3,ncf_3,nsh_3,nbf_3,
     *                    inx_4,ibas_4,ncs_4,ncf_4,nsh_4,nbf_4)
c----------------------------------------------------------
c setup the infob common :
c
      inucx=inuc
      ibasx=ibas
      nax  =na
      ndum=igetival('ndum')
      if(ndum.gt.0) nax=nax-ndum
      nbfx =nbf
      nshx =nsh
      ncfx =ncf
      ncsx =ncs
c
      nsymx=nsym
      rnsym=one/(one+float(nsym))
c----------------------------------------------------------
c setup logical matrices : common  /logicd/, logic1-11/
c
       call get_maxtype(bl(inx_1),ncs_1, maxtype_1)
       call get_maxtype(bl(inx_2),ncs_2, maxtype_2)
       call get_maxtype(bl(inx_3),ncs_3, maxtype_3)
       call get_maxtype(bl(inx_4),ncs_4, maxtype_4)
       maxtype=max(maxtype_1,maxtype_2,maxtype_3,maxtype_4)
c
       call datlog(maxtype,lpar1,lpar2,lpar3,lpar4,lpar5)
c
c----------------------------------------------------------
c calculte the lind vector: lind(i)=i*(i-1)/2
c
       ncf_max=max(ncf_1,ncf_2,ncf_3,ncf_4)
       call getmem(ncf_max,lind)
       call setup_lind(bl(lind),ncf_max)
c
c----------------------------------------------------------
      call setival('iforwhat',iforwhat)
c----------------------------------------------------------
       call secund(tim3)
c----------------------------------------------------------
c find out what blocking will be used :
c
      call block_route(iforwhat,nax,bl(inuc),iroute)
c
c output : iroute=1 or 2 (stored in common /route/ iroute
c
      call setival('iroute',iroute)
c----------------------------------------------------------
c
c     Blocking procedure for contracted shells (single shells)
c
      call blockin1(iforwhat,na,ncs,inx,bl,bl(inuc),bl(ibas),nbl1)
c
c     output : nbl1 - number of blocks of single shells
c----------------------------------------------------------
c allocated memory for Schwarz integrals (ij|ij) :
c
      call getmem(ncs*ncs,ischwarz)
      call setival('schwarz',ischwarz)
c----------------------------------------------------------
      call mmark
c----------------------------------------------------------
c     blocking procedure for pairs of contracted shells
c
      call blockin2(lcore,nprint,ncs,inx,bl,nbl1,
     *              nbl2,minpr,maxpr,  xminpr_bl, xmaxpr_bl) ! outputs
c
c nbl2      - number of blocks of shell-pairs
c minpr     - price of the cheapest integral
c maxpr     - price of the most expensive integral
c xminpr_bl - minimum price for the whole big-block
c xmaxpr_bl - maximum price for the whole big-block
c----------------------------------------------------------
       call secund(tim4)
       timblok2=tim4-tim3
c....................................................
       nbloks=nbl2*(nbl2+1)/2     !    number of super-blocks :
c....................................................
c----------------------------------------------------------
       call run_mode(bl,inx,ncf,minpr,maxpr,xminpr_bl,xmaxpr_bl,
     *               lsemi,scftype,1)
c----------------------------------------------------------
c calculate Schwarz integrals (ij|ij) :
cccc  maxbuffer=maxibuf
      maxlabels=maxindx
      call getmem(maxlabels,labels)
         call int_schwarz(bl,inx,bl(ischwarz),bl(labels),'sch1')
      call retmem(1)
c----------------------------------------------------------
c print schwarz integrals
c     call print_schw(ncs,bl(ischwarz))
c----------------------------------------------------------
      call retmark
c----------------------------------------------------------
c Set up the integral threshold as requested for this task:
c (it was overwritten in Schwarz integ. calculations 10**-20
c
      call setup_thres(eps)
c----------------------------------------------------------
c do blocking2 again using schwarz integrals to eliminate some
c contracted shell pairs : re-define pair & quartet blocks :
c
      call blockin2_again(lcore,nprint,ncs,inx,bl,nbl1,
     *              nbl2,minpr,maxpr,  xminpr_bl, xmaxpr_bl) ! outputs
c
      nbloks=nbl2*(nbl2+1)/2     !    number of super-blocks :
c----------------------------------------------------------
c    
      if(iforwhat.lt.11) then    ! for all but ftc
      call run_mode(bl,inx,ncf,minpr,maxpr,xminpr_bl,xmaxpr_bl,
     *              lsemi,scftype,0)
      else
c        for plane waves update_run_mode will be called from rhf5.f
         scftype='full-direct'
      endif
c    
c     call run_mode(bl,inx,ncf,minpr,maxpr,xminpr_bl,xmaxpr_bl,
c    *              lsemi,scftype,0)
c----------------------------------------------------------
      if(iforwhat.eq.5) then
         write(ioutput,*) ' this is going to be a special integral task'
         write(ioutput,*) '       iforwhat=',iforwhat
         write(ioutput,*) '     blocking will be done on a request'
      endif
c----------------------------------------------------------
      call secund(tim2)
      timprep2=tim2-tim1
      write(ioutput,66) timprep2/60.d0
  66  format('  cpu time for preparation of 2-el. step :',f6.2,' min')
c----------------------------------------------------------
      end
c================================================================
      subroutine block_route(iforwhat,natoms,datnuc,iroute)
c---------------------------------------------------------------
c Input
c-------
c iforwhat  - shows type of task (integrals to be calculated)
c natoms    - number of atoms
c datnuc()  - nuclear data
c
c Output
c-------
c iroute    - blocking strategy (1 or 2 i.e. 3 or 5 conditions)
c---------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      character*4 text
      common /outfile/ ioutput
      common /cpu/ intsize,iacc,icache,memreal
      common /intlim/ limxmem,limblks,limpair
      dimension datnuc(5,*)
c---------------------------------------------------------------
c Check out how many atoms of the same type the molecule contains
c and decide which blocking procedure will be used :
c
      if(iroute.eq.0) call whichblk(datnuc,natoms,iroute)
c
c output iroute (=1 or =2) ( stored in the common /route/
c
      if(iforwhat.eq.1) text='SCF '
      if(iforwhat.eq.2) text='GIAO'
      if(iforwhat.eq.3) text='GRAD'
      if(iforwhat.eq.4) text='HESS'
      if(iforwhat.eq.5) text='LMP2'
c
      if(iforwhat.eq.11) text='FTC '
      if(iforwhat.eq.13) text='FTCG'
c
      write(ioutput,900)
      write(ioutput,910) text,iroute
      write(ioutput,900)
c
  900 format('=================================================')
  910 format('Blocking procedure for ',a4,' integrals : route=',i2)
c---------------------------------------------------------------
      write(ioutput,*)'    Maximum sizes of blocks upon restrictions :'
      write(ioutput,*)'         cache memory limit : ',icache
      write(ioutput,*)'         pairs/block2 limit : ',limpair
      if(limblks.GT.0) then
         write(ioutput,*)'         quart/block4 limit : ',limblks
      else
         write(ioutput,*)'         quart/block4 limit : adjusted'
      endif
      write(ioutput,*)'         memory       limit : ',limxmem
      write(ioutput,900)
c---------------------------------------------------------------
      end
c===============================================================
      subroutine whichblk(datnuc,natoms,iroute)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
c---------------------------------------------------------------
c Natoms - total number of atoms
c Ntypes - number of different types of atoms
c
c if one atom contributes to the basis set
c more than 1/ntypes %                       go with iroute=1
c
c if number of atoms of the same type
c is greater then 1/ntypes % of all atoms
c           and
c
c they contribute to the basis set more
c than 1/ntypes %                            go with iroute=2
c---------------------------------------------------------------
c                              charge
c     I-row    H  -  He       1  -   2       2 atoms
c    II-row   Li  -  Ne       3  -  10       8 atoms
c   III-row   Na  -  Ar      11  -  18       8 atoms
c    IV-row    K  -  Kr      19  -  36      18 atoms
c     V-row   Rb  -  Xe      37  -  54      18 atoms
c    VI-row   Cs  -  Rn      55  -  86      32 atoms
c   VII-row   Fr  -  Ha      87  - 105      19 atoms .....
c                                 no of basis
c                                  functions
c---------------------------------------------------------------
      dimension nbas(7)
      dimension nrow(118)
      dimension nats(0:118)     ! number of atoms with a given Q
c
c nbas(irow)-> no of b.f.: number of basis function per atom
c                          with respect to the ROW
c nrow(iq)  -> element iq belongs to row nrow(iq)
c
      data nbas /2,10,18,36,54,86,118/
      data nrow /1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,
     *           4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,4,
     *           5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,5,
     *  6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,6,
     *  7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7,7/
c
      parameter (zero=0.0d0)
c---------------------------------------------------------------
      call izeroit(nats,118)
c---------------------------------------------------------------
      iroute=1
      if(natoms.eq.1) return
c---------------------------------------------------------------
      if(natoms.gt.14) then
        iroute=2
        return
      endif
c---------------------------------------------------------------
c
      iqmax=0
      qtot=zero
      do 10 iat=1,natoms
        iq=nint( datnuc(1,iat) )
        if(iq.gt.0) then
           qtot=qtot+iq
           if(iq.gt.iqmax) iqmax=iq
        endif
   10 continue
c
c
      do 20 iat=1,natoms
        iq=nint( datnuc(1,iat) )
        if(iq.gt.0) then
          nats(iq)=nats(iq)+1
        end if
   20 continue
c
      ntypes=0
      do 30 iq=1,iqmax
      if(nats(iq).gt.0) ntypes=ntypes+1
   30 continue
c---------------------------------------------------------------
      if(ntypes.eq.1) then
         iroute=2
         return
      endif
c---------------------------------------------------------------
c
c calculate maximum number of atoms of the same type
c
      nsame=0
      do iq=1,iqmax
         if(nats(iq).ge.nsame) then
            nsame=nats(iq)
            iqsame=iq
         endif
      enddo
c
c calculate how many functions comes from nsame atoms
c
      nbfsame=nats(iqsame)*nbas( nrow(iqsame) )
c
c calculate number of functions on the heaviest atom
c
      nbfheavy=nbas( nrow(iqmax) )
c
c calculate total number of basis functions and max funct. from the same atoms
c
      nbftotal=0
      nbfq_max=0
      do iq=1,iqmax
        nbfiq   = nats(iq)*nbas( nrow(iq) )
        nbftotal= nbftotal+nbfiq
        if(nbfiq.ge.nbfq_max) then
           nbfq_max=nbfiq
           iq_max=iq
        endif
      enddo
c---------------------------------------------------------------
c     write(6,*) 'total num. of bf=',nbftotal
c
c     write(6,*)' total number of atoms      =',natoms
c     write(6,*)' number of different types  =',ntypes
c     write(6,*)' total number of basis funct=',nbftotal
c     write(6,*)' aver. number of bf on type =',nbftotal/ntypes
c     write(6,*)' basis funct on the heaviest=',nbfheavy
c     write(6,*) nsame,' atoms with Q=',iqsame,' give ',nbfsame,' b.f.'
c     write(6,*)
c    *   nats(iq_max),' atoms with Q=',iq_max,' give ',nbfq_max,' b.f.'
c---------------------------------------------------------------
      if(nbfheavy.GT.nbftotal/ntypes) then
         iroute=1
         RETURN
      endif
c---------------------------------------------------------------
      if(nsame.ge.natoms/ntypes) then
         if(nbfsame.gt.nbftotal/ntypes) then
            iroute=2
            RETURN
         endif
      endif
c---------------------------------------------------------------
      if(nats(iq_max).ge.natoms/ntypes) then
         if(nbfq_max.gt.nbftotal/ntypes) then
            iroute=2
            RETURN
         endif
      endif
c---------------------------------------------------------------
c based no TMS SiC4H12 :
c
      if(nsame.gt.3 .and. nats(iq_max).gt.3) then
c          12H               4C
         if(nsame+nats(iq_max).ge.2*natoms/ntypes) then
c             12H + 4C             2*17/3
             if(nbfsame+nbfq_max.gt.2*nbftotal/ntypes) then
c                 12       40        2* 82/3
                iroute=2
                RETURN
             endif
         endif
      endif
c---------------------------------------------------------------
c at least 3 atoms of the same type :
c     nat_min=1000000
c     do iq=1,iqmax
c        iat=nats(iq)
c        if(iat.gt.0) nat_min=min(iat,nat_min)
c     enddo
c     if(nat_min.ge.3) then
c        iroute=2
c     endif
c---------------------------------------------------------------
c
      end
c===============================================================
      subroutine timepr(tblock)
      implicit real*8 (a-h,o-z)
      common /outfile/ ioutput
c
      common /timex/ tconv1,tconv2,ttrobs
      common /time0/ tprec2
      common /time1/ tpre4n,txwpq ,tassem,tamshf,tdesti
      common /time2/ terint,ttrans
      common /time3/ tzeroi,tspeci
      common /time4/ tderiv
c
      data zero /0.d0/
c
      if(tblock.eq.zero) then
ctime
        tprec2=zero
        tconv1=zero
        tconv2=zero
        ttrobs=zero
        tpre4n=zero
        txwpq =zero
        tassem=zero
        tamshf=zero
        tdesti=zero
        terint=zero
        tzeroi=zero
        tspeci=zero
        ttrans=zero
        tderiv=zero
ctime
      else
c
        write(ioutput,120) tblock,tprec2,tpre4n,tzeroi,txwpq,tconv1,
     *                     tconv2,
     *          ttrobs,tassem,tderiv,tamshf,terint,tspeci,ttrans,tdesti
c
  120   format(
     *  /'***************************',2x,'***************************'/
     *        '   S U B R O U T I N E       T I M I N G  (in sec)     '/
     *   '---------------------------',2x,'---------------------------'/
     *        'time for BLOCKIN=',f10.2,'  time for PRECL2 =',f10.2/
     *   '---------------------------',2x,'---------------------------'/
     *   'time for prec4n =',f10.2/
     *   'time for zeroin =',f10.2/
     *   'time for xwpq   =',f10.2/
     *   'time for conv1x =',f10.2/
     *   'time for conv2x =',f10.2/
     *   'time for trobsa =',f10.2/
     *   'time for assemb =',f10.2/
     *   'time for derivx =',f10.2/
     *   'time for amshif =',f10.2/
     *   '- - - - - - - - - - - - - -'/
     *   'time for ERINTEG=',f10.2,'  time for ERINTSP=',f10.2/
     *   '---------------------------',2x,'---------------------------'/
     *    'time for TRANSFO=',f10.2,'  time for DESTINY=',f10.2/
     *   '***************************',2x,'***************************')
c
      endif
c
      end
c================================================================
      subroutine setup_lind(lind,ncf)
      dimension lind(ncf)
      do 10 i=1,ncf
      lind(i)=i*(i-1)/2
   10 continue
      end
c=============================================================
      subroutine setup_thres(thres)
      implicit real*8 (a-h,o-z)
c Set up the integral threshold :
c----------------------------------------------------------------
c Estimator1=16/(pi)**1/2 * Sab*Scd * 1/(a+b+c+d)**1/2
c
c  where Sab=(a+b)**-1 * (a*b)**+3/4 *exp(-ab/(a+b) * Rab**2 )
c
c
c    Estimator1 > eps    do integrals .
c
c It is done in a quadratic form by  Es1*Es1> eps**2
c
c    Sab**2  *  Scd**2 > pi/256 * eps**2  * (a+b+c+d)
c
c  epsr= pi/256 * eps**2
c
c  epsx=eps -original thresh. for integrals,
c  eps1=1/eps  used to rescale integr. in precalc2a routine
c  epsr=pi/256 *eps**2-used in prec4neg and precspec routines.
c----------------------------------------------------------------
c pi256=pi/256
c
      common /neglect/ eps,eps1,epsr,eps8
      data one,pi, pi256 /1.d0 , 3.1415926535d0 , 0.012271846d0 /
c
      eps =thres
      eps1=one/eps
      eps8=eps1*8.d0
      epsr=pi256*eps*eps
c
c----------------------------------------------------------
      end
c=======================================================================
c NEW routines for new blocking :
c temorary :
      subroutine set_basis_data(bl,inx,datbas,ncs,ncf,nsh,nbf)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c     common /intbl/maxsh,inxx(100)
      dimension bl(*)
      dimension inx(12,ncs)
      dimension datbas(13,nsh)
c----------------------------------
c It makes 4 copies of a given basis set
c----------------------------------
      nsh13=nsh*13
      ncs12=ncs*12
c----------------------------------
c basis set no 1:
c allocate memory:
      call getint(ncs  ,list_1)  ! for list of shells in basis set 1
      call getint(ncs12,inx_1)   ! for inx_1
      call getmem(nsh13,ibas_1)     ! for datbas_1
c save addresses:
      call setival('list_1',list_1)
      call setival('inx_1' ,inx_1 )
      call setival('ibas_1',ibas_1)
      call setival('ncs_1' ,ncs )
      call setival('ncf_1' ,ncf )
      call setival('nsh_1' ,nsh )
      call setival('nbf_1' ,nbf )
c temporary copying data from inx and datbas into :
      call copy_inxs(ncs12,inx, bl(inx_1) )
      call dcopy(nsh13, datbas(1,1),1, bl(ibas_1),1)
c----------------------------------
c basis set no 2:
c allocate memory:
      call getint(ncs  ,list_2)
      call getint(ncs12,inx_2)
      call getmem(nsh13,ibas_2)
c save addresses:
      call setival('list_2',list_2)
      call setival('inx_2' ,inx_2 )
      call setival('ibas_2',ibas_2)
      call setival('ncs_2' ,ncs )
      call setival('ncf_2' ,ncf )
      call setival('nsh_2' ,nsh )
      call setival('nbf_2' ,nbf )
c temporary copying data
      call copy_inxs(ncs12,inx, bl(inx_2) )
      call dcopy(nsh13, datbas(1,1),1, bl(ibas_2),1)
c----------------------------------
c basis set no 3:
c allocate memory:
      call getint(ncs  ,list_3)
      call getint(ncs12,inx_3)
      call getmem(nsh13,ibas_3)
c save addresses:
      call setival('list_3',list_3)
      call setival('inx_3' ,inx_3 )
      call setival('ibas_3',ibas_3)
      call setival('ncs_3' ,ncs )
      call setival('ncf_3' ,ncf )
      call setival('nsh_3' ,nsh )
      call setival('nbf_3' ,nbf )
c temporary copying data
      call copy_inxs(ncs12,inx, bl(inx_3) )
      call dcopy(nsh13, datbas(1,1),1, bl(ibas_3),1)
c----------------------------------
c basis set no 4:
c allocate memory:
      call getint(ncs  ,list_4)
      call getint(ncs12,inx_4)
      call getmem(nsh13,ibas_4)
c save addresses:
      call setival('list_4',list_4)
      call setival('inx_4' ,inx_4 )
      call setival('ibas_4',ibas_4)
      call setival('ncs_4' ,ncs )
      call setival('ncf_4' ,ncf )
      call setival('nsh_4' ,nsh )
      call setival('nbf_4' ,nbf )
c temporary copying data
      call copy_inxs(ncs12,inx, bl(inx_4) )
      call dcopy(nsh13, datbas(1,1),1, bl(ibas_4),1)
c-------------------------------------------------------------------------
c     make list_x ; here we put there number of shells from 1 to ncs
c     it can be overwtitten later by removing some of these shells
c     These lists will be used in blocking1 to make blocks of single-shells
c
      call list_of_shells(bl(list_1),ncs)
      call list_of_shells(bl(list_2),ncs)
      call list_of_shells(bl(list_3),ncs)
      call list_of_shells(bl(list_4),ncs)
c----------------------------------
      end
c----------------------------------
      subroutine list_of_shells(list,ncs)
      dimension list(*)
      do ics=1,ncs
         list(ics)=ics
      enddo
c----------------------------------
      end
c-------------
      subroutine copy_inxs(ncs12,inx, inx_1 )
      dimension inx(ncs12), inx_1(ncs12)
c
      do i=1,ncs12
         inx_1(i)=inx(i)
      enddo
c
      end
c=======================================================================
c NEW routines for new blocking :
c
      subroutine get_basis_data(inx_1,ibas_1,ncs_1,ncf_1,nsh_1,nbf_1,
     *                          inx_2,ibas_2,ncs_2,ncf_2,nsh_2,nbf_2,
     *                          inx_3,ibas_3,ncs_3,ncf_3,nsh_3,nbf_3,
     *                          inx_4,ibas_4,ncs_4,ncf_4,nsh_4,nbf_4)
c---------------------------------------------------------------
c All output parameters :
c inx_X  - address for INX(12,*) arrays
c ibas_X - address for datbas(13,*) arrays
c and
c ncs_x  - number of contracted shells in basis set X
c ncf_x  - number of contracted functions in basis set X
c nsh_x  - number of primitive shells     in basis set X
c nbf_x  - number of primitive  functions in basis set X
c---------------------------------------------------------------
cccc  call getival('iout',iout)
c---------------------------------------------------------------
      call getival('inx_1' ,inx_1 )
      call getival('ibas_1',ibas_1)
      call getival('ncs_1' ,ncs_1 )
      call getival('ncf_1' ,ncf_1 )
      call getival('nsh_1' ,nsh_1 )
      call getival('nbf_1' ,nbf_1 )
c
      call getival('inx_2' ,inx_2 )
      call getival('ibas_2',ibas_2)
      call getival('ncs_2' ,ncs_2 )
      call getival('ncf_2' ,ncf_2 )
      call getival('nsh_2' ,nsh_2 )
      call getival('nbf_2' ,nbf_2 )
c
      call getival('inx_3' ,inx_3 )
      call getival('ibas_3',ibas_3)
      call getival('ncs_3' ,ncs_3 )
      call getival('ncf_3' ,ncf_3 )
      call getival('nsh_3' ,nsh_3 )
      call getival('nbf_3' ,nbf_3 )
c
      call getival('inx_4' ,inx_4 )
      call getival('ibas_4',ibas_4)
      call getival('ncs_4' ,ncs_4 )
      call getival('ncf_4' ,ncf_4 )
      call getival('nsh_4' ,nsh_4 )
      call getival('nbf_4' ,nbf_4 )
c
      end
c=======================================================================
      subroutine get_maxtype(inx, ncs, maxtype)
      dimension inx(12,ncs)
c-----------------------------------------------------------------------
c     correspondance between ityp and ityp1:
c     ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28)
c     ityp1 1    2    3    4    4     5    5      6      7       8
c     ityp  11(g)  12(h)  13(i)
c     ityp1 6       7      8
c-----------------------------------------------------------------------
c     maxtyp1=0
c     do ics=1,ncs
c        itype=inx(12,ics)
c        itype1=itype
c        if(itype.gt.4) itype1=itype-1
c        if(itype1.gt.5) itype1=itype1-1
c        if(itype1.gt.maxtyp1) maxtyp1=itype1
c     enddo
c     maxtype=maxtyp1
c----------------------------------------------
      if(ncs.eq.0) then
         maxtype=0
         return
      endif
c----------------------------------------------
c
c for a reordered basis set the last shell is of the highest type
c
      itype=inx(12,ncs)
      itype1=itype
      if(itype.gt.4) itype1=itype-1
      if(itype1.gt.5) itype1=itype1-1
      if(itype.eq.11) itype1=6
      if(itype.eq.12) itype1=7
      if(itype.eq.13) itype1=8
      maxtype=itype1
c
      end
c=======================================================================
      subroutine print_schw(ncs,schwarz)
      implicit real*8 (a-h,o-z)
      dimension schwarz(ncs,ncs)
      do ics=1,ncs
         do jcs=1,ics
            write(6,66) ics,jcs,schwarz(ics,jcs)
         enddo
      enddo
   66 format(2(i3,2x),5x,f12.6)
      end

c...
c...  MM September-October 2004
c...
c...  thys file contains the driver routines for the calculation
c...  of analytical hessians over ab initio pseudopotentials.
c...
c======================================================================
      subroutine pspFder(rhf,na,natb,nate,npsp,ncf,ncs,inx,
     $                   xnuc,basdat,listreal,fda,fdb)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  driver routine for the psp contribution to the derivative
c...  Fock matrices
c...
c...  rhf     flag for rhf/uhf
c...  na      number of atoms
c...  natb,nate  beginning and end of current pass
c...  npsp    number of pseudopotentials
c...  ncf     number of contracted functions
c...  ncs     number of contracted shells
c...  inx     contraction information
c...  xnuc    nuclear data
c...  basdat  basis set data
c...  fda     alpha (or closed shell) derivative Fock matrix
c...  fdb     beta derivative fock matrix
c...
      logical rhf
      dimension inx(12,*),basdat(13,*),xnuc(5,*),fda(*),fdb(*)
      dimension listreal(*)
c     common /big/bl(1)
      common /forcdbl/ thre1,thre2,tchf
      dimension time(20)
c...
c     call secund(t0)
c     call zeroit(time,20)
c...
c...
c...  get the psp parameters and pointers
c...
      call getival('nradp',nradp)       ! number of radial terms
      call getival('maxlpsp',maxlpsp)   ! maximum l value we can handle
      call getival('iiecpdat',iiecpdat) ! ptr. to psp info
      call getival('inrrad',inrrad)     ! ptr. to r exponents
      call getival('iarad',iarad)       ! ptr. to radial exponents
      call getival('icrad',icrad)       ! ptr. to radial coefficients
c...
c...  parameters for memory setup:
c...
c...  maxlp        maximum angular momentum of psp
c...  maxl         maximum angular momentum of Gaussians
c...  maxd2        maximum number of individual components for
c...               of a contracted shell (number of gaussians times
c...               contraction length or no. of general contractions)
c...
      call pspgetpar(npsp,maxlpsp,ncs,inx,bl(iiecpdat),maxlp,maxl,
     $               maxd2)
c...
c...  increase maxl, because we are dealing with derivatives
c...  of the Gaussians
c...
      maxl=maxl+1
c...
c...  number of atoms for which integral second derivatives
c...  over primitives have to be stored. This is 2 (the centers of
c...  the basis functions) plus the number of pseudopotentials,
c...  as the integrals for a given shell pair are accumulated
c...  over the pseudopotentials
c...
      natder=2+npsp
c...
c...  memory setup
c...
      mlp=max(maxl,maxlp)
      maxpw=max(maxl,2)+1
c...
c...  allocate static memory (this memory is always needed)
c...
      call mmark
      call getmem(maxpw,icax)             ! powers of ca
      call getmem(maxpw,icay)             ! powers of ca
      call getmem(maxpw,icaz)             ! powers of ca
      call getmem(maxpw,icbx)             ! powers of cb
      call getmem(maxpw,icby)             ! powers of cb
      call getmem(maxpw,icbz)             ! powers of cb
      call getmem((maxl+mlp+1)**2,iyra)   ! real sperical harmonics
      call getmem((maxl+maxlp+1)**2,iyrb) ! real sperical harmonics
      call getmem((2*maxl+1)*(maxl+maxlp+1)**2,iq) ! radial integrals
      call getmem((4*mlp+1)**3,ixyzi)     ! precomputed angular terms
      call getmem((maxl+1)**6,iar) ! sums of radial and angular int.
      call getmem(maxlp+1,ilen)           ! array for psp contraction
      call getmem(natder,ilder)           ! mapping of derivatives
c...
c...  now the memory that can be dynamically used
c...
c...  compute free memory
c...
      call getmem(0,last0)
      call retmem(1)
      call getival('lcore',lcore)
      nmemory=lcore-last0
c...
c...  we need two buffers of the same size
c...
      membuf=int(nmemory/2)-1
      npspi=3*natder*(maxd2)**2
c...
c...  npspi is the maximum size of the buffer.
c...
c...  if we have enough memory, we allocate the maximum
c...  size, otherwise we allocate membuf words only, and hope
c...  for the best
c...
      if(npspi.le.membuf)then
        call getmem(npspi,ipspi)            !int. over primitives
        call getmem(npspi,ipspic)           !int. over contracted
      else
        call getmem(membuf,ipspi)            !int. over primitives
        call getmem(membuf,ipspic)           !int. over contracted
      endif
c...
c...  generate a table of double factorial products
c...  for angular integrals
c...
      mijk=4*mlp
      call preomega(mijk,bl(ixyzi))
c...
c...  set thresold to 0.01 times the one electron integral threshold
c...
      tol=0.01d0*thre1
c...
c...  loop over contracted shells
c...
      maxpass=0
      do ics=1,ncs
        do jcs=1,ics
c...
c...  compute integrals
c...
          call psp_Fder(ics,      jcs,       npsp,     nradp,   maxlpsp,
     $                  ncs,      ncf,      maxlp,      maxl,      mijk,
     $         bl(iiecpdat),bl(inrrad), bl(iarad), bl(icrad),  bl(ilen),
     $                 xnuc,       inx,    basdat,  bl(icax),  bl(icay),
     $             bl(icaz),  bl(icbx),  bl(icby),  bl(icbz),  bl(iyra),
     $             bl(iyrb),    bl(iq), bl(ixyzi),   bl(iar), bl(ipspi),
     $           bl(ipspic), bl(ilder),    membuf,       rhf,        na,
     $                 natb,      nate,       tol,  listreal,       fda,
     $                  fdb,      time,   maxpass)
        enddo
      enddo
c...
c...  done. Release memory
c...
      call retmark
c     call secund(t1)
c...
c     call getival('iout',iout)
c     write(iout,*)
c     write(iout,*)'Total time for psp contribution       ',t1-t0
c     write(iout,*)
c     write(iout,*)'Overall integrals over primitives     ',time(8)
c     write(iout,*)'Zeroing integrals arrays              ',time(9)
c     write(iout,*)'Time for radial and angular integrals ',time(1)
c     write(iout,*)'Time for integrals over primitives    ',time(7)
c     write(iout,*)'       Time for summation loop        ',time(2)
c     write(iout,*)'       Time for assembling integrals  ',time(6)
c     write(iout,*)'Time for normalization/contraction    ',time(3)
c     write(iout,*)'Time for cartesian to spherical       ',time(4)
c     write(iout,*)'Time for fillup                       ',time(5)
c     write(iout,*)
c     write(iout,*)'Maximum number of passes              ',maxpass
c     write(iout,*)
c...
      end
c======================================================================
      subroutine psp_Fder(ics,     jcs,    npsp,   nradp,maxlpsp,
     $                    ncs,     ncf,   maxlp,    maxl,   mijk,
     $                iecpdat,   nrrad,   aradp,   cradp,    len,
     $                   xnuc,     inx,  basdat,     cax,    cay,
     $                    caz,     cbx,     cby,     cbz,    yra,
     $                    yrb,       q,    xyzi,      ar,   pspi,
     $                  pspic,    lder,  membuf,     rhf,     na,
     $                   natb,    nate,     tol,listreal,    fda,
     $                    fdb,    time, maxpass)
      implicit real*8 (a-h,o-z)
c...
c...  compute psp contribution to fock derivative matrices for
c..   the contracted shell pair ics,jcs
c...
c...  ics        contracted shell i
c...  jcs        contracted shell j
c...  npsp       number of pseudopotentials
c...  nradp      number of radial terms
c...  maxlpsp    maximum angular momentum of psp we can handle
c...  ncs        number of contracted shells
c...  ncf        number of contracted basis functrions
c...  maxlp      maximum angular momentum of psp
c...  maxl       maximum angular momentum of Gaussians
c...  mijk       dimension for xyzi array
c...  iecpdat    array of psp info
c...  nrrad      array of r exponent of radial terms
c...  aradp      array of radial term exponents
c...  cradp      array of radial term coefficients
c...  len        scratch area for psp contraction length
c...  xnuc       array of nuclear info
c...  inx        array of contraction info
c...  basdat     array of  basis set data
c...  cax        scratch area for powers of cax
c...  cay        scratch area for powers of cay
c...  caz        scratch area for powers of caz
c...  cbx        scratch area for powers of cbx
c...  cby        scratch area for powers of cby
c...  cbz        scratch area for powers of cbz
c...  yra        scratch area for real spherical harmonics
c...  yrb        scratch area for real spherical harmonics
c...  q          scratch area for radial integrals
c...  xyzi       area for precomputed angular terms
c...  ar         scratch area for sums of radial and angular ints.
c...  pspi       scratch area for integrals over Gaussians primitives
c...  pspic      scratch area for integrals over Gaussians contracted
c...  lder       scratch array for mapping derivatives to atoms
c...  membuf     maximum size of the pspi and pspic buffers
c...  rhf        flag for rhf/uhf
c...  na         number of atoms
c...  natb       beginning of this pass
c...  nate       ending of this pass
c...  tol        tolerance threshold
c...  listreal   array of atoms to do (only symmetry unique atoms)
c...  fda        alpha (or closed shell) Fock derivative matrix
c...  fdb        beta derivative Fock matrix
c...
c...  each integral is differentiated with respect the centers
c...  of the two Gaussians involved in the integral. The derivatives
c...  with respect the pseudopotential center are computed by
c...  translational invariance
c...
      dimension iecpdat(maxlpsp+6,npsp),nrrad(nradp),aradp(nradp),
     $          cradp(nradp),len(0:maxlp),lder(*)
      dimension xnuc(5,*),inx(12,*),basdat(13,*)
      dimension cax(0:maxl),cay(0:maxl),caz(0:maxl)
      dimension cbx(0:maxl),cby(0:maxl),cbz(0:maxl)
      dimension pspi(*),pspic(*)
      dimension listreal(*)
      logical rhf
      dimension fda(3,natb:nate,*),fdb(3,natb:nate,*)
      dimension time(*)
c...
      dimension A(3),B(3),C(3),CA(3),CB(3),lima(3),limb(3)
c...
c...  the type of the Gaussian shell is given by ityp:
c...
c...  ityp  1  2  3   4   5   6    7    8   9    10   11  12   13
c...        s  p  l   d5  d6  f7  f10  g15  h21  i28  g9  h11  i13
c...
c...   the integrals are computed over the cartesian components
c...   (d6, f15, g15, h21, i28). The spherical harmonic components,
c...   if needed, will be obtained by linear combination.
c...
c...   data for number of components and angular momenta
c...
      dimension ityp1(13),ng(8),ns(8),lmax(8),mulx(88),muly(88),
     $          mulz(88)
      data ityp1/1,2,3,4,4,5,5,6,7,8,6,7,8/
      data ng/1,3,4,6,10,15,21,28/
      data ns/1,2,5,9,15,25,40,61/
      data lmax/0,1,1,2,3,4,5,6/
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
c...
      logical transa,transb
      logical doiat,dojat,docat
      parameter (zero=0.0d0)
      parameter (argmax=100.0d0)
      tolp=tol
c...
c...  shell i:
c...
      ngci=inx(4,ics)+1          ! no. of general contractions
      ifi=inx(1,ics)             ! first primitive
      ifl=inx(5,ics)             ! last primitive
      ilen=ifl-ifi               ! contraction length
      iatom=inx(2,ics)           ! atomic center
      iau=listreal(iatom)        ! iau=0 if iatom is not unique
      doiat=(iau.ge.natb.and.iau.le.nate)!is atom i in current pass?
      a(1)=xnuc(2,iatom)
      a(2)=xnuc(3,iatom)
      a(3)=xnuc(4,iatom)
      itypa=inx(12,ics)          ! type of shell i
      transa=itypa.eq.4.or.itypa.eq.6.or.itypa.ge.11
      ityp1a=ityp1(itypa)
      lmaxa=lmax(ityp1a)         ! angular momentum
      nga=ng(ityp1a)             ! no. of cartesian components
      nsa=ns(ityp1a)             ! ptr. for angular momenta vector
      nsha=inx(3,ics)            ! shell size
      ifci=inx(11,ics)           ! first contracted (minus one)
c...
c...  shell j:
c...
      ngcj=inx(4,jcs)+1          ! no. of general contractions
      jfi=inx(1,jcs)             ! first primitive
      jfl=inx(5,jcs)             ! last primitive
      jlen=jfl-jfi               ! contraction length
      jatom=inx(2,jcs)           ! atomic center
      jau=listreal(jatom)        ! jau=0 if jatom is not unique
      dojat=(jau.ge.natb.and.jau.le.nate)!is atom j in current pass?
      b(1)=xnuc(2,jatom)
      b(2)=xnuc(3,jatom)
      b(3)=xnuc(4,jatom)
      itypb=inx(12,jcs)          ! type of shell j
      transb=itypb.eq.4.or.itypb.eq.6.or.itypb.ge.11
      ityp1b=ityp1(itypb)
      lmaxb=lmax(ityp1b)         ! angular momentum
      ngb=ng(ityp1b)             ! no. of cartesian components
      nsb=ns(ityp1b)             ! ptr. for angular momenta vector
      nshb=inx(3,jcs)            ! shell size
      ifcj=inx(11,jcs)           ! first contracted (minus one)
      mxpw=max(2,lmaxa+1,lmaxb+1)
c...
c...  derivative to atom mapping
c...
c...  at the end of this section, nder will contain the
c...  number of atomic derivatives that have to be computed
c...  for this shell pair in the current pass
c...
      nder=0
      if(doiat)then
        nder=nder+1
        iat=nder
        lder(nder)=iau
      else
        iat=0
      endif
      if(dojat)then
        if(jau.ne.iau)then
          nder=nder+1
          jat=nder
          lder(nder)=jau
        else
          jat=iat
        endif
      else
        jat=0
      endif
      nderij=nder
      do ipsp=1,npsp
        icatom=iecpdat(3,ipsp)
        icau=listreal(icatom)
        docat=(icau.ge.natb.and.icau.le.nate)
        if(docat)then
          icat=0
          do id=1,nder
            if(icau.eq.lder(id))icat=id
          enddo
          if(icat.eq.0)then
            nder=nder+1
            lder(nder)=icau
          endif
        endif
      enddo
      nderc=nder-nderij
c...
c...  if no integrals have to be computed, we are done
c...
      if(nder.eq.0)return
c...
c...  check if memory buffer is large enough
c...
c...  npass  =1 we can proceed (with multiple passes, if needed)
c...         =0 not ehough memory
c...
c...  nderp     maximum number of atomic derivatives that
c...            can be stored in one pass
c...
      nints=nga*ngb*max(ilen*jlen,ngci*ngcj)
      mder=3*nder*nints
      if(mder.le.membuf)then
        npass=1
        nderp=nder
      else if(nderc.eq.0)then
          npass=0
      else
        nderp=int(membuf/(nints*3))
        if(nderp.ge.nderij+1)then
          npass=1
        else
          npass=0
        endif
      endif
      if(npass.le.0)then
         call getival('iout',iout)
         write(iout,*)
         write(iout,*)'Psp_Fder: not enough memory.'
         write(iout,*)'At least ',9*nints-membuf,' more words needed.'
         write(iout,*)
         call nerror(1,'Psp_Fder','Not enough memory',0,0)
      endif
c...
c...  start passes
c...
      ipass=0
      iniz=1
 1    ipass=ipass+1     ! starting point for  passes
      maxpass=max(maxpass,ipass)
c     call getival('iout',iout)
c     if(ipass.gt.1)write(iout,*)'Psp_Fder is doing multiple passes'
c...
c...  derivative to atom mapping for current pass.
c...
c...  nder    number of atomic derivatives computed by current pass
c...  iniz    first pseudopotential to be included in current pass
c...  ifin     last pseudopotential to be included in current pass
c...
        nder=nderij
        ifin=iniz
        do ipsp=iniz,npsp
          icatom=iecpdat(3,ipsp)
          icau=listreal(icatom)
          docat=(icau.ge.natb.and.icau.le.nate)
          if(docat)then
            icat=0
            do id=1,nder
              if(icau.eq.lder(id))icat=id
            enddo
            if(icat.eq.0)then
              nder=nder+1
              lder(nder)=icau
              if(nder.gt.nderp)then
                nder=nderp   ! memory limit for current pass reached
                goto 2
              endif
            endif
          endif
          ifin=ipsp
        enddo
2       continue
c...
c...  zero out integrals over primitives
c...
        nder3=3*nder
c       call secund(t0)
        call zeroit(pspi,nder3*nints)
        call zeroit(pspic,nder3*nints)
c       call secund(tz)
c       time(9)=time(9)+tz-t0
c...
c...  loop over pseudopotentials
c...
        do ipsp=iniz,ifin
          icatom=iecpdat(3,ipsp)   ! atomic center of psp
          icau=listreal(icatom)    ! icau=0 is icatom is not unique
          docat=(icau.ge.natb.and.icau.le.nate)!is atom c in pass?
          if(docat)then
            do id=1,nder
              if(icau.eq.lder(id))icat=id !pointer to atom in der. array
            enddo
            if(icat.eq.0)
     $        call nerror(2,'Psp_Fder','pointer to icatom is 0 !!!',0,0)
          else
            icat=0
          endif
          lpsp=iecpdat(1,ipsp)     ! maximum l value of psp
          irl=iecpdat(4,ipsp)
          irs=iecpdat(5,ipsp)
c...
c...  do current integrals only if one of the centers is computed
c...  by the current pass
c...
          if(doiat.or.dojat.or.docat)then
c...
c...  compute number of terms of local and non local parts
c...
            len0=irs-irl             ! number of terms in local part
            do l=0,lpsp-1,1
              len(l)=iecpdat(6+l,ipsp)-iecpdat(5+l,ipsp)
            enddo
            len(lpsp)=iecpdat(12,ipsp)-iecpdat(5+lpsp,ipsp)+1
c...
c...   distances CA and CB
c...
            c(1)=xnuc(2,icatom)
            c(2)=xnuc(3,icatom)
            c(3)=xnuc(4,icatom)
            CA(1)=C(1)-A(1)
            CA(2)=C(2)-A(2)
            CA(3)=C(3)-A(3)
            CB(1)=C(1)-B(1)
            CB(2)=C(2)-B(2)
            CB(3)=C(3)-B(3)
c...
c...  compute powers of distances ca and cb and do loops limits
c...
            call zeroit(cax,mxpw+1)
            call zeroit(cay,mxpw+1)
            call zeroit(caz,mxpw+1)
            call zeroit(cbx,mxpw+1)
            call zeroit(cby,mxpw+1)
            call zeroit(cbz,mxpw+1)
            call powcacb(CA(1),CA(2),CA(3),CB(1),CB(2),CB(3),
     $                   mxpw,lima,limb,cax,cay,caz,cbx,cby,cbz)
            ca2=cax(2)+cay(2)+caz(2)
            cb2=cbx(2)+cby(2)+cbz(2)
c...
c...  accumulate integral derivatives over primitives
c...
            call psp_grad12i(
     $             basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $          muly(nsb), mulz(nsb),       nga,      ngb,       lpsp,
     $         nrrad(irl),aradp(irl),cradp(irl),     len0,        len,
     $                cax,       cay,       caz,      cbx,        cby,
     $                cbz,      lima,      limb,     mxpw,    lmaxa+1,
     $            lmaxb+1,      mijk,       yra,      yrb,       xyzi,
     $                  q,        ar,      nder,      iat,        jat,
     $               icat,      tolp,    argmax,     ilen,        ifi,
     $                ca2,      jlen,       jfi,      cb2,       pspi,
     $               time)
          endif
        enddo   ! end loop over pseudopotentials
c       call secund(t1)
c       time(8)=time(8)+t1-t0
c...
c...  contraction and normalization.
c...
c       call secund(t0)
        call psp_contrn(ilen,ifi,ngci,nga,itypa,
     $                  jlen,jfi,ngcj,ngb,itypb,
     $                  basdat,nder3,pspi,pspic)
c       call secund(t1)
c       time(3)=time(3)+t1-t0
c...
c...  eventually transform to spherical components.
c...  pspi is used as scratch space
c...
c       call secund(t0)
        if(transa.or.transb)then
          call cart2sphera(pspic,nder3*ngci*ngcj,itypa,transa,nga,nsha,
     $                             itypb,transb,ngb,nshb,pspi)
        endif
c       call secund(t1)
c       time(4)=time(4)+t1-t0
c...
c...  add contribution to Fock derivatives
c...
c       call secund(t0)
        ici=ifci
        nd3=ngci*ngcj*nga
        nd2=ngci*ngcj
        do igc=1,ngci
          do ica=1,nsha
            ici=ici+1
            ii=ici*(ici-1)/2
            icj=ifcj
            do jgc=1,ngcj
              do icb=1,nshb
                icj=icj+1
                if(ici.ge.icj)then
                   iijj=ii+icj
                  iab=nder3*((icb-1)*nd3+(ica-1)*nd2+(jgc-1)*ngci+igc-1)
                  do id=1,nder
                    idat=lder(id)
                    if(rhf)then
                      iab=iab+1
                      fda(1,idat,iijj)=fda(1,idat,iijj)-pspic(iab)
                      iab=iab+1
                      fda(2,idat,iijj)=fda(2,idat,iijj)-pspic(iab)
                      iab=iab+1
                      fda(3,idat,iijj)=fda(3,idat,iijj)-pspic(iab)
                    else
                      iab=iab+1
                      fda(1,idat,iijj)=fda(1,idat,iijj)-pspic(iab)
                      fdb(1,idat,iijj)=fdb(1,idat,iijj)-pspic(iab)
                      iab=iab+1
                      fda(2,idat,iijj)=fda(2,idat,iijj)-pspic(iab)
                      fdb(2,idat,iijj)=fdb(2,idat,iijj)-pspic(iab)
                      iab=iab+1
                      fda(3,idat,iijj)=fda(3,idat,iijj)-pspic(iab)
                      fdb(3,idat,iijj)=fdb(3,idat,iijj)-pspic(iab)
                    endif
                  enddo
                endif
              enddo
            enddo
          enddo
        enddo
c       call secund(t1)
c       time(5)=time(5)+t1-t0
c...
c...  prepare for next pass
c...
        iniz=ifin+1
      if(iniz.le.npsp)goto 1
c...
      end
c======================================================================
      subroutine psphess(na,npsp,ncf,ncs,inx,xnuc,basdat,dens,hess)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  driver routine for the psp contribution to the hessian
c...
c...  na      number of atoms
c...  npsp    number of pseudopotentials
c...  ncf     number of contracted functions
c...  ncs     number of contracted shells
c...  inx     contraction information
c...  xnuc    nuclear data
c...  basdat  basis set data
c...  dens    density matrix
c...  hess    cumulative hessian
c...
      dimension inx(12,*),basdat(13,*),xnuc(5,*),dens(*),hess(*)
c     common /big/bl(1)
      parameter (zero=0.0d0)
      common /forcdbl/ thre1,thre2,tchf
      dimension time(20)
c...
c     call secund(t0)
c     call zeroit(time,20)
c...
c...
c...  get the psp parameters and pointers
c...
      call getival('nradp',nradp)       ! number of radial terms
      call getival('maxlpsp',maxlpsp)   ! maximum l value we can handle
      call getival('iiecpdat',iiecpdat) ! ptr. to psp info
      call getival('inrrad',inrrad)     ! ptr. to r exponents
      call getival('iarad',iarad)       ! ptr. to radial exponents
      call getival('icrad',icrad)       ! ptr. to radial coefficients
c...
c...  parameters for memory setup:
c...
c...  maxlp        maximum angular momentum of psp
c...  maxl         maximum angular momentum of Gaussians
c...  maxd2        maximum number of individual components for
c...               of a contracted shell (number of gaussians times
c...               contaction length or no. of general contractions)
c...  lcore        total memory available for the calculation
c...
      call pspgetpar(npsp,maxlpsp,ncs,inx,bl(iiecpdat),maxlp,maxl,
     $               maxd2)
c...
c...  increase maxl, because we are dealing with second derivatives
c...  of the Gaussians
c...
      maxl=maxl+2
c...
c...  number of atoms for which integral second derivatives
c...  over primitives have to be stored. This is 2 (the centers of
c...  the basis functions) plus the number of pseudopotentials,
c...  as the integrals for a given shell pair are accumulated
c...  over the pseudopotentials
c...
      natder=2+npsp
c...
c...  memory setup
c...
      mlp=max(maxl,maxlp)
      maxpw=max(maxl,2)+1
c...
c     call getival('iout',iout)
c     write(iout,*)'memory parameters'
c     write(iout,*)
c     write(iout,*)'maxl      =',maxl
c     write(iout,*)'maxlp     =',maxlp
c     write(iout,*)'maxd2     =',maxd2
c     write(iout,*)'natder    =',natder
c     write(iout,*)
c     write(iout,*)'mlp       =',mlp
c     write(iout,*)'maxpw     =',maxpw
c     write(iout,*)'maxd2     =',maxd2
c     write(iout,*)
c     write(iout,*)'memory for powers    =',6*maxpw
c     write(iout,*)'memory for yra       =',(maxl+mlp+1)**2
c     write(iout,*)'memory for yrb       =',(maxl+maxlp+1)**2
c     write(iout,*)'memory for q         =',(2*maxl+1)*(maxl+maxlp+1)**2
c     write(iout,*)'memory for ar        =',(maxl+1)**6
c     write(iout,*)'memory for xyzi      =',(4*mlp+1)**3
c     write(iout,*)'memory for len       =',maxlp+1
c     write(iout,*)
c...
c...  allocate static memory (this memory is always needed)
c...
      call mmark
      call getmem(maxpw,icax)             ! powers of ca
      call getmem(maxpw,icay)             ! powers of ca
      call getmem(maxpw,icaz)             ! powers of ca
      call getmem(maxpw,icbx)             ! powers of cb
      call getmem(maxpw,icby)             ! powers of cb
      call getmem(maxpw,icbz)             ! powers of cb
      call getmem((maxl+mlp+1)**2,iyra)   ! real sperical harmonics
      call getmem((maxl+maxlp+1)**2,iyrb) ! real sperical harmonics
      call getmem((2*maxl+1)*(maxl+maxlp+1)**2,iq) ! radial integrals
      call getmem((4*mlp+1)**3,ixyzi)     ! precomputed angular terms
      call getmem((maxl+1)**6,iar)    ! sums of radial and angular int.
      call getmem(natder,ilder)           ! mapping of derivatives
      call getmem(natder,iatder)          ! derivative processing
      call getmem(maxd2**2,inotc)
      call getmem(maxlp+1,ilen)           ! array for psp contraction
c...
c...  now the memory that can be dynamically used
c...
c...  compute free memory
c...
      call getmem(0,last0)
      call retmem(1)
      call getival('lcore',lcore)
      nmemory=lcore-last0
c...
c...  we need two buffers of the same size
c...
      membuf=int(nmemory/2)-1
      npspi=9*natder**2*(maxd2)**2
c...
c     call getival('iout',iout)
c     write(iout,*)'memory for pspi      =',npspi
c     write(iout,*)'memory for pspic     =',npspi
c     write(iout,*)'memory buffer size   =',membuf
c     write(iout,*)
c...
c...  npspi is the maximum size of the buffer.
c...
c...  if we have enough memory, we allocate the maximum
c...  size, otherwise we allocate membuf words only, and hope
c...  for the best
c...
      if(npspi.le.membuf)then
        call getmem(npspi,ipspi)            !int. over primitives
        call getmem(npspi,ipspic)           !int. over contracted
      else
        call getmem(membuf,ipspi)            !int. over primitives
        call getmem(membuf,ipspic)           !int. over contracted
      endif
c...
c...  generate a table of double factorial products
c...  for angular integrals
c...
      mijk=4*mlp
      call preomega(mijk,bl(ixyzi))
c...
c...  set thresold to 0.01 times the one electron integral threshold
c...
      tol=0.01d0*thre1
c...
c...  loop over contracted shells
c...
      maxpass=0
      do ics=1,ncs
        do jcs=1,ics
c...
c...  compute integrals
c...
          call psp_hess(ics,       jcs,      npsp,     nradp,   maxlpsp,
     $                  ncs,       ncf,     maxlp,      maxl,      mijk,
     $         bl(iiecpdat),bl(inrrad), bl(iarad), bl(icrad),  bl(ilen),
     $                 xnuc,       inx,    basdat,  bl(icax),  bl(icay),
     $             bl(icaz),  bl(icbx),  bl(icby),  bl(icbz),  bl(iyra),
     $             bl(iyrb),    bl(iq), bl(ixyzi),   bl(iar), bl(ipspi),
     $           bl(ipspic), bl(ilder),bl(iatder), bl(inotc),    membuf,
     $                 dens,        na,       tol,      hess,      time,
     $   maxpass)
        enddo
      enddo
c...
c...  done. Release memory
c...
      call retmark
c     call secund(t1)
c...
c     call getival('iout',iout)
c     write(iout,*)
c     write(iout,*)'Total time for psp contribution       ',t1-t0
c     write(iout,*)
c     write(iout,*)'Overall integrals over primitives     ',time(8)
c     write(iout,*)'Zeroing integrals arrays              ',time(9)
c     write(iout,*)'Time for radial and angular integrals ',time(1)
c     write(iout,*)'Time for integrals over primitives    ',time(7)
c     write(iout,*)'       Time for summation loop        ',time(2)
c     write(iout,*)'       Time for assembling integrals  ',time(6)
c     write(iout,*)'Time for normalization/contraction    ',time(3)
c     write(iout,*)'Time for cartesian to spherical       ',time(4)
c     write(iout,*)'Time for fillup                       ',time(5)
c     write(iout,*)
c     write(iout,*)'Maximum number of passes              ',maxpass
c     write(iout,*)
c...
c...
      end
c======================================================================
      subroutine psp_hess(ics,  jcs,   npsp,  nradp,maxlpsp,
     $                    ncs,  ncf,  maxlp,   maxl,   mijk,
     $                iecpdat,nrrad,  aradp,  cradp,    len,
     $                   xnuc,  inx, basdat,    cax,    cay,
     $                    caz,  cbx,    cby,    cbz,    yra,
     $                    yrb,    q,   xyzi,     ar,   pspi,
     $                  pspic, lder,  atder,   notc, membuf,
     $                   dens,   na,    tol,   hess,   time,
     $                maxpass)
      implicit real*8 (a-h,o-z)
c...
c...  compute psp contribution to atomic hessian for the contracted
c...  shell pair ics,jcs
c...
c...  ics        contracted shell i
c...  jcs        contracted shell j
c...  npsp       number of pseudopotentials
c...  nradp      number of radial terms
c...  maxlpsp    maximum angular momentum of psp we can handle
c...  ncs        number of contracted shells
c...  ncf        number of contracted basis functrions
c...  maxlp      maximum angular momentum of psp
c...  maxl       maximum angular momentum of Gaussians
c...  mijk       dimension for xyzi array
c...  iecpdat    array of psp info
c...  nrrad      array of r exponent of radial terms
c...  aradp      array of radial term exponents
c...  cradp      array of radial term coefficients
c...  len        scratch area for psp contraction length
c...  xnuc       array of nuclear info
c...  inx        array of contraction info
c...  basdat     array of  basis set data
c...  cax        scratch area for powers of cax
c...  cay        scratch area for powers of cay
c...  caz        scratch area for powers of caz
c...  cbx        scratch area for powers of cbx
c...  cby        scratch area for powers of cby
c...  cbz        scratch area for powers of cbz
c...  yra        scratch area for real spherical harmonics
c...  yrb        scratch area for real spherical harmonics
c...  q          scratch area for radial integrals
c...  xyzi       area for precomputed angular terms
c...  ar         scratch area for sums of radial and angular ints.
c...  pspi       scratch area for integrals over Gaussians primitives
c...  pspic      scratch area for integrals over Gaussians contracted
c...  lder       scratch array for mapping derivatives to atoms
c...  membuf     maximum size of the pspi and pspic buffers
c...  dens       density matrix
c...  na         number of atoms
c...  tol        tolerance threshold
c...  hess       hessian
c...
c...  each integral is differentiated with respect the centers
c...  of the two Gaussians involved in the integral. The derivatives
c...  with respect the pseudopotential center are computed by
c...  translational invariance
c...
      dimension iecpdat(maxlpsp+6,npsp),nrrad(nradp),aradp(nradp),
     $          cradp(nradp),len(0:maxlp),lder(*)
      logical atder(*)
      dimension xnuc(5,*),inx(12,*),basdat(13,*)
      dimension cax(0:maxl),cay(0:maxl),caz(0:maxl)
      dimension cbx(0:maxl),cby(0:maxl),cbz(0:maxl)
      dimension pspi(*),pspic(*)
      dimension dens(*),hess(*)
      dimension time(*)
c...
      dimension A(3),B(3),C(3),CA(3),CB(3),lima(3),limb(3)
c...
c...  the type of the Gaussian shell is given by ityp:
c...
c...  ityp  1  2  3   4   5   6    7    8   9    10   11  12   13
c...        s  p  l   d5  d6  f7  f10  g15  h21  i28  g9  h11  i13
c...
c...   the integrals are computed over the cartesian components
c...   (d6, f15, g15, h21, i28). The spherical harmonic components,
c...   if needed, will be obtained by linear combination.
c...
c...   data for number of components and angular momenta
c...
      dimension ityp1(13),ng(8),ns(8),lmax(8),mulx(88),muly(88),
     $          mulz(88)
      data ityp1/1,2,3,4,4,5,5,6,7,8,6,7,8/
      data ng/1,3,4,6,10,15,21,28/
      data ns/1,2,5,9,15,25,40,61/
      data lmax/0,1,1,2,3,4,5,6/
      data mulx/0, 1,0,0, 0,1,0,0, 2,0,0,1,1,0, 3,2,2,1,1,1,0,0,0,0,
     1          4,3,3,2,2,2,1,1,1,1,0,0,0,0,0, 5,4,4,3,3,3,2,2,2,2,
     2          1,1,1,1,1,0,0,0,0,0,0,  6,5,5,4,4,4,3,3,3,3,2,2,2,2,2,
     3          1,1,1,1,1,1,0,0,0,0,0,0,0/
      data muly/0, 0,1,0, 0,0,1,0, 0,2,0,1,0,1, 0,1,0,2,1,0,3,2,1,0,
     1          0,1,0,2,1,0,3,2,1,0,4,3,2,1,0, 0,1,0,2,1,0,3,2,1,0,
     2          4,3,2,1,0,5,4,3,2,1,0,  0,1,0,2,1,0,3,2,1,0,4,3,2,1,0,
     3          5,4,3,2,1,0,6,5,4,3,2,1,0/
      data mulz/0, 0,0,1, 0,0,0,1, 0,0,2,0,1,1, 0,0,1,0,1,2,0,1,2,3,
     1          0,0,1,0,1,2,0,1,2,3,0,1,2,3,4, 0,0,1,0,1,2,0,1,2,3,
     2          0,1,2,3,4,0,1,2,3,4,5,  0,0,1,0,1,2,0,1,2,3,0,1,2,3,4,
     3          0,1,2,3,4,5,0,1,2,3,4,5,6/
c...
      logical transa,transb,dorest
      parameter (zero=0.0d0)
      parameter (argmax=100.0d0)
      tolp=tol
c...
c...  shell i:
c...
      ngci=inx(4,ics)+1          ! no. of general contractions
      ifi=inx(1,ics)             ! first primitive
      ifl=inx(5,ics)             ! last primitive
      ilen=ifl-ifi               ! contraction length
      iatom=inx(2,ics)           ! atomic center
      a(1)=xnuc(2,iatom)
      a(2)=xnuc(3,iatom)
      a(3)=xnuc(4,iatom)
      itypa=inx(12,ics)          ! type of shell i
      transa=itypa.eq.4.or.itypa.eq.6.or.itypa.ge.11
      ityp1a=ityp1(itypa)
      lmaxa=lmax(ityp1a)         ! angular momentum
      nga=ng(ityp1a)             ! no. of cartesian components
      nsa=ns(ityp1a)             ! ptr. for angular momenta vector
      nsha=inx(3,ics)            ! shell size
      ifci=inx(11,ics)           ! first contracted (minus one)
c...
c...  shell j:
c...
      ngcj=inx(4,jcs)+1          ! no. of general contractions
      jfi=inx(1,jcs)             ! first primitive
      jfl=inx(5,jcs)             ! last primitive
      jlen=jfl-jfi               ! contraction length
      jatom=inx(2,jcs)           ! atomic center
      b(1)=xnuc(2,jatom)
      b(2)=xnuc(3,jatom)
      b(3)=xnuc(4,jatom)
      itypb=inx(12,jcs)          ! type of shell j
      transb=itypb.eq.4.or.itypb.eq.6.or.itypb.ge.11
      ityp1b=ityp1(itypb)
      lmaxb=lmax(ityp1b)         ! angular momentum
      ngb=ng(ityp1b)             ! no. of cartesian components
      nsb=ns(ityp1b)             ! ptr. for angular momenta vector
      nshb=inx(3,jcs)            ! shell size
      ifcj=inx(11,jcs)           ! first contracted (minus one)
      mxpw=max(2,lmaxa+2,lmaxb+2)
c...
c...  derivative to atom mapping
c...
c...  at the end of this section, nder will contain the
c...  total number of atomic derivatives that have to be computed
c...  for this shell pair in the current pass
c...
      nder=1
      iat=nder
      lder(nder)=iatom
      if(jatom.ne.iatom)then
        nder=nder+1
        jat=nder
        lder(nder)=jatom
      else
        jat=iat
      endif
      nderij=nder
      do ipsp=1,npsp
        icatom=iecpdat(3,ipsp)
        icat=0
        do id=1,nder
          if(icatom.eq.lder(id))icat=id
        enddo
        if(icat.eq.0)then
          nder=nder+1
          lder(nder)=icatom
        endif
      enddo
      nderc=nder-nderij
c...
c...  check if memory buffer is large enough
c...
c...  npass  =1 we can proceed (with multiple passes, if needed)
c...         =0 not ehough memory
c...
c...  nderp     maximum number of atomic derivatives that
c...            can be stored in one pass
c...
      ninti=nga*ngb*ilen*jlen
      nintc=nga*ngb*ngci*ngcj
      nints=max(ninti,nintc)
      mder=9*nder**2*nints
      if(mder.le.membuf)then
        npass=1
        nderp=nder
      else if(nderc.eq.0)then
          npass=0
      else
        nderp=int(sqrt(dfloat(membuf/(nints*9))))
        if(nderp.ge.nderij+1)then
          npass=1
        else
          npass=0
        endif
      endif
      if(npass.le.0)then
         call getival('iout',iout)
         write(iout,*)
         write(iout,*)'Psp_hess: not enough memory.'
         write(iout,*)'At least ',81*nints-membuf,' more words needed.'
         write(iout,*)
         call nerror(1,'Psp_hess','Not enough memory',0,0)
      endif
c...
c...  start passes
c...
      ipass=0
      iniz=1
 1    ipass=ipass+1     ! starting point for  passes
      maxpass=max(maxpass,ipass)
c     call getival('iout',iout)
c     if(ipass.gt.1)write(iout,*)'Psp_hess is doing multiple passes'
c...
c...  derivative to atom mapping for current pass.
c...
c...  nder    number of atomic derivatives computed by current pass
c...  iniz    first pseudopotential to be included in current pass
c...  ifin     last pseudopotential to be included in current pass
c...
        nder=nderij
        ifin=iniz
        do ipsp=iniz,npsp
          icatom=iecpdat(3,ipsp)
          icat=0
          do id=1,nder
            if(icatom.eq.lder(id))icat=id
          enddo
          if(icat.eq.0)then
            nder=nder+1
            lder(nder)=icatom
            if(nder.gt.nderp)then
              nder=nderp   ! memory limit for current pass reached
              goto 2
            endif
          endif
          ifin=ipsp
        enddo
2       continue
c...
c...  zero out integrals
c...
        nder3=3*nder
        na3=3*na
        na2=9*(nder)**2
c       call secund(t0)
c       call zeroit(pspi,na2*ninti)
        do id=1,nder
          atder(id)=.false.
        enddo
c       call secund(tz)
c       time(9)=time(9)+tz-t0
c...
c...  loop over pseudopotentials
c...
        do ipsp=iniz,ifin
          icatom=iecpdat(3,ipsp)   ! atomic center of psp
          icat=0
          do id=1,nder
            if(icatom.eq.lder(id))icat=id !pointer to atom in der. array
          enddo
          if(icat.eq.0)
     $      call nerror(2,'Psp_hess','pointer to icatom is 0 !!!',0,0)
          lpsp=iecpdat(1,ipsp)     ! maximum l value of psp
          irl=iecpdat(4,ipsp)
          irs=iecpdat(5,ipsp)
c...
c...  compute number of terms of local and non local parts
c...
          len0=irs-irl             ! number of terms in local part
          do l=0,lpsp-1,1
            len(l)=iecpdat(6+l,ipsp)-iecpdat(5+l,ipsp)
          enddo
          len(lpsp)=iecpdat(12,ipsp)-iecpdat(5+lpsp,ipsp)+1
c...
c...   distances CA and CB
c...
          c(1)=xnuc(2,icatom)
          c(2)=xnuc(3,icatom)
          c(3)=xnuc(4,icatom)
          CA(1)=C(1)-A(1)
          CA(2)=C(2)-A(2)
          CA(3)=C(3)-A(3)
          CB(1)=C(1)-B(1)
          CB(2)=C(2)-B(2)
          CB(3)=C(3)-B(3)
c         if(debug)then
c           call getival('iout',iout)
c           write(iout,*)'lpsp,len0,len',lpsp,len0,(len(ijz),ijz=0,lpsp)
c           write(iout,*)'nr',(nrrad(irl+ijz),ijz=0,len0-1)
c           write(iout,*)'coeff',(cradp(irl+ijz),ijz=0,len0-1)
c           write(iout,*)'exp',(aradp(irl+ijz),ijz=0,len0-1)
c           write(iout,*)'a',a
c           write(iout,*)'b',b
c           write(iout,*)'c',c
c           write(iout,*)'ca',ca
c           write(iout,*)'cb',cb
c           write(iout,*)'mxpw',mxpw
c         endif
c...
c...  compute powers of distances ca and cb and do loops limits
c...
          call zeroit(cax,mxpw+1)
          call zeroit(cay,mxpw+1)
          call zeroit(caz,mxpw+1)
          call zeroit(cbx,mxpw+1)
          call zeroit(cby,mxpw+1)
          call zeroit(cbz,mxpw+1)
          call powcacb(CA(1),CA(2),CA(3),CB(1),CB(2),CB(3),
     $               mxpw,lima,limb,cax,cay,caz,cbx,cby,cbz)
          ca2=cax(2)+cay(2)+caz(2)
          cb2=cbx(2)+cby(2)+cbz(2)
c...
c...  accumulate integral second derivatives over primitives
c...
           call psp_hess12i(
     $               basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $            muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $           nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $                  cax,       cay,       caz,       cbx,       cby,
     $                  cbz,      lima,      limb,      mxpw,   lmaxa+2,
     $              lmaxb+2,      mijk,       yra,       yrb,      xyzi,
     $                    q,        ar,      nder,       iat,       jat,
     $                 icat,      tolp,    argmax,      ilen,       ifi,
     $                  ca2,      jlen,       jfi,      cb2,       pspi,
     $                atder,      notc,      time)
        enddo   ! end loop over pseudopotentials
c       call secund(t1)
c       time(8)=time(8)+t1-t0
c...
c... do the rest of this pass only if integrals were computed
c...
        dorest=.false.
        do id=1,nder
          if(atder(id))dorest=.true.
        enddo
        if(dorest)then
c...
c...  contraction and normalization.
c...
c         call secund(t0)
          call zeroit(pspic,na2*nintc)
          call psp_contrn2(   ilen,    ifi,   ngci,    nga,  itypa,
     $                        jlen,    jfi,   ngcj,    ngb,  itypb,
     $                      basdat,    na2,   nder,  nder3,   pspi,
     $                       pspic,  atder, nderij)
c         call secund(t1)
c         time(3)=time(3)+t1-t0
c...
c...  eventually transform to spherical components.
c...  pspi is used as scratch space
c...
          if(transa.or.transb)then
            call cart2sphera(pspic,na2*ngci*ngcj,itypa,transa,nga,nsha,
     $                             itypb,transb,ngb,nshb,pspi)
          endif
c...
c...  add contribution to atomic hessian
c...
c         call secund(t0)
          ici=ifci
          nd3=ngci*ngcj*nga
          nd2=ngci*ngcj
          do igc=1,ngci
            do ica=1,nsha
              ici=ici+1
              ii=ici*(ici-1)/2
              icj=ifcj
              do jgc=1,ngcj
                do icb=1,nshb
                  icj=icj+1
                  if(ici.ge.icj)then
                    ind=na2*((icb-1)*nd3+(ica-1)*nd2+(jgc-1)*ngci+igc-1)
                    dij=dens(ii+icj)
                    if(ici.ne.icj)dij=dij+dij
c...
c...  first loop over atoms i and j
c...
                    do id=1,nderij
                      if(atder(id))then
                        idat=lder(id)
                        ida3=3*(idat-1)
                        idx=ida3+1
                        idy=ida3+2
                        idz=ida3+3
                        id3=3*(id-1)
                        icpx=id3+1+ind
                        icpy=id3+2+ind
                        icpz=id3+3+ind
                        do jd=1,nderij
                          if(atder(jd))then
                            jdat=lder(jd)
                            jda3=3*(jdat-1)
                            jdx=na3*jda3
                            jdy=na3*(jda3+1)
                            jdz=na3*(jda3+2)
                            jd3=3*(jd-1)
                            jcpx=nder3*jd3
                            jcpy=nder3*(jd3+1)
                            jcpz=nder3*(jd3+2)
                    hess(idx+jdx)=hess(idx+jdx)+dij*pspic(icpx+jcpx)!xx
                    hess(idy+jdx)=hess(idy+jdx)+dij*pspic(icpy+jcpx)!yx
                    hess(idz+jdx)=hess(idz+jdx)+dij*pspic(icpz+jcpx)!zx
                    hess(idx+jdy)=hess(idx+jdy)+dij*pspic(icpx+jcpy)!xy
                    hess(idy+jdy)=hess(idy+jdy)+dij*pspic(icpy+jcpy)!yy
                    hess(idz+jdy)=hess(idz+jdy)+dij*pspic(icpz+jcpy)!zy
                    hess(idx+jdz)=hess(idx+jdz)+dij*pspic(icpx+jcpz)!xz
                    hess(idy+jdz)=hess(idy+jdz)+dij*pspic(icpy+jcpz)!yz
                    hess(idz+jdz)=hess(idz+jdz)+dij*pspic(icpz+jcpz)!zz
                          endif
                        enddo
                      endif
                    enddo
c...
c...  now loop over the pseudopotential centers
c...
                    do id=nderij+1,nder
                      if(atder(id))then
                        idat=lder(id)
                        ida3=3*(idat-1)
                        idx=ida3+1
                        idy=ida3+2
                        idz=ida3+3
                        id3=3*(id-1)
                        icpx=id3+1+ind
                        icpy=id3+2+ind
                        icpz=id3+3+ind
                        do jd=1,nderij
                          if(atder(jd))then
                            jdat=lder(jd)
                            jda3=3*(jdat-1)
                            jdx=na3*jda3
                            jdy=na3*(jda3+1)
                            jdz=na3*(jda3+2)
                            jd3=3*(jd-1)
                            jcpx=nder3*jd3
                            jcpy=nder3*(jd3+1)
                            jcpz=nder3*(jd3+2)
                    hess(idx+jdx)=hess(idx+jdx)+dij*pspic(icpx+jcpx)!xx
                    hess(idy+jdx)=hess(idy+jdx)+dij*pspic(icpy+jcpx)!yx
                    hess(idz+jdx)=hess(idz+jdx)+dij*pspic(icpz+jcpx)!zx
                    hess(idx+jdy)=hess(idx+jdy)+dij*pspic(icpx+jcpy)!xy
                    hess(idy+jdy)=hess(idy+jdy)+dij*pspic(icpy+jcpy)!yy
                    hess(idz+jdy)=hess(idz+jdy)+dij*pspic(icpz+jcpy)!zy
                    hess(idx+jdz)=hess(idx+jdz)+dij*pspic(icpx+jcpz)!xz
                    hess(idy+jdz)=hess(idy+jdz)+dij*pspic(icpy+jcpz)!yz
                    hess(idz+jdz)=hess(idz+jdz)+dij*pspic(icpz+jcpz)!zz
                          endif
                        enddo
                        jdat=lder(id)
                        jda3=3*(jdat-1)
                        jdx=na3*jda3
                        jdy=na3*(jda3+1)
                        jdz=na3*(jda3+2)
                        jd3=3*(id-1)
                        jcpx=nder3*jd3
                        jcpy=nder3*(jd3+1)
                        jcpz=nder3*(jd3+2)
                    hess(idx+jdx)=hess(idx+jdx)+dij*pspic(icpx+jcpx)!xx
                    hess(idy+jdx)=hess(idy+jdx)+dij*pspic(icpy+jcpx)!yx
                    hess(idz+jdx)=hess(idz+jdx)+dij*pspic(icpz+jcpx)!zx
                    hess(idx+jdy)=hess(idx+jdy)+dij*pspic(icpx+jcpy)!xy
                    hess(idy+jdy)=hess(idy+jdy)+dij*pspic(icpy+jcpy)!yy
                    hess(idz+jdy)=hess(idz+jdy)+dij*pspic(icpz+jcpy)!zy
                    hess(idx+jdz)=hess(idx+jdz)+dij*pspic(icpx+jcpz)!xz
                    hess(idy+jdz)=hess(idy+jdz)+dij*pspic(icpy+jcpz)!yz
                    hess(idz+jdz)=hess(idz+jdz)+dij*pspic(icpz+jcpz)!zz
                        do jd=1,nderij
                          if(atder(jd))then
                            idat=lder(jd)
                            ida3=3*(idat-1)
                            idx=ida3+1
                            idy=ida3+2
                            idz=ida3+3
                            id3=3*(jd-1)
                            icpx=id3+1+ind
                            icpy=id3+2+ind
                            icpz=id3+3+ind
                    hess(idx+jdx)=hess(idx+jdx)+dij*pspic(icpx+jcpx)!xx
                    hess(idy+jdx)=hess(idy+jdx)+dij*pspic(icpy+jcpx)!yx
                    hess(idz+jdx)=hess(idz+jdx)+dij*pspic(icpz+jcpx)!zx
                    hess(idx+jdy)=hess(idx+jdy)+dij*pspic(icpx+jcpy)!xy
                    hess(idy+jdy)=hess(idy+jdy)+dij*pspic(icpy+jcpy)!yy
                    hess(idz+jdy)=hess(idz+jdy)+dij*pspic(icpz+jcpy)!zy
                    hess(idx+jdz)=hess(idx+jdz)+dij*pspic(icpx+jcpz)!xz
                    hess(idy+jdz)=hess(idy+jdz)+dij*pspic(icpy+jcpz)!yz
                    hess(idz+jdz)=hess(idz+jdz)+dij*pspic(icpz+jcpz)!zz
                          endif
                        enddo
                      endif
                    enddo
                  endif
                enddo
              enddo
            enddo
          enddo
c         call secund(t1)
c         time(5)=time(5)+t1-t0
        endif
c...
c...  prepare for next pass
c...
        iniz=ifin+1
      if(iniz.le.npsp)goto 1
c...
      end
c=======================================================================
      subroutine psp_contrn2(   ilen,    ifi,   ngci,    nga,  itypa,
     $                          jlen,    jfi,   ngcj,    ngb,  itypb,
     $                        basdat,   ndim,   nder,  nder3,   pspi,
     $                         pspic,  atder, nderij)
      implicit real*8 (a-h,o-z)
c...
c...  normalizes and contracts a batch of psp integrals over
c...  primitives for both sells I and j
c...
c...  ilen        contraction length of shell i
c...  ifi         initial pointer to primitives of shell i
c...  ngci        number of general contractions on shell i
c...  nga         number of Gaussians of shell i (cartesian)
c...  itypa       type of shell i
c...  jlen        contraction length of shell j
c...  jfi         initial pointer to primitives of shell j
c...  ngcj        number of general contractions on shell j
c...  ngb         number of Gaussians of shell i (cartesian)
c...  itypb       type of shell j
c...  basdat      basis set info
c...  ndim        mumber of components of integrals
c...  nder        number of atomic derivatives handled by current pass
c...  nder3=3*nder
c...  pspi        integrals over primitives
c...  pspic       in exit will contain the integrals over contracted
c...  atder       flag for derivatives to contract
c...  nderij      derivatives with respect basis function center.
c...
c...  a few words on normalization constants:
c...
c...  each Gaussian need to be multiplied by a normalization factor
c...
c...         N1=(2*a/%pi)**(3/4)*(4*a)**(l/2)
c...
c...   where a is the exponent and l is the angular momentum.
c...   (this assuming the integrals are normalized to
c...     (2*n_x-1)!!*(2*n_y-1)!!*(2*n_z-1)!!)
c...
c...  the factor (1/%pi)**(3/4) is already taken into account
c...  during the calculation of the radial integrals, while
c...  the factor (4*a)**(l/2) has already been included
c...  in the normalization of the primitives.
c...  Thus here the integrals need to be multipled only by the factor
c...
c...               NF=(2*a)**(3/4)
c...
c...        note that there are two such factors, one for each
c...        of the two Gaussians involved in the integral.
c...
      dimension basdat(13,*)
      dimension pspi(nder3,nder3,ilen,jlen,nga,ngb)
      dimension pspic(nder3,nder3,ngci,ngcj,nga,ngb)
      logical atder(*)
c...
      do igc=1,ngci
        do jgc=1,ngcj
          do ip=1,ilen
            ipi=ifi+ip
            alpa=basdat(1,ipi)
            if(itypa.eq.3)then
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
              cpa=basdat(3,ipi)*(2.0d0*alpa)**0.75d0
            else
              csa=basdat(igc+1,ipi)*(2.0d0*alpa)**0.75d0
            endif
            do jp=1,jlen
              jpj=jfi+jp
              alpb=basdat(1,jpj)
              if(itypb.eq.3)then
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
                cpb=basdat(3,jpj)*(2.0d0*alpb)**0.75d0
              else
                csb=basdat(jgc+1,jpj)*(2.0d0*alpb)**0.75d0
              endif
              do i=1,nga
                if(itypa.eq.3.and.i.gt.1)then
                  coefi=cpa
                else
                  coefi=csa
                endif
                do j=1,ngb
                  if(itypb.eq.3.and.j.gt.1)then
                    coefj=cpb
                  else
                    coefj=csb
                  endif
                  coefij=coefi*coefj
c...
c...  first loop over components of atoms i and j (centers of
c...  basis functions)
c...
                  do iat=1,nderij
                    if(atder(iat))then
                      ia3=3*(iat-1)
                      do jat=1,nderij
                        if(atder(jat))then
                          ja3=3*(jat-1)
                          pspic(ia3+1,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+1,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+1,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+3,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+3,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+3,ip,jp,i,j)*coefij
                        endif
                      enddo
                    endif
                  enddo
c...
c... now loop over pseudopotential centers
c...
                  do iat=nderij+1,nder
                    if(atder(iat))then
                      ia3=3*(iat-1)
                      pspic(ia3+1,ia3+1,igc,jgc,i,j)=
     $                pspic(ia3+1,ia3+1,igc,jgc,i,j)+
     $                 pspi(ia3+1,ia3+1,ip,jp,i,j)*coefij
                      pspic(ia3+1,ia3+2,igc,jgc,i,j)=
     $                pspic(ia3+1,ia3+2,igc,jgc,i,j)+
     $                 pspi(ia3+1,ia3+2,ip,jp,i,j)*coefij
                      pspic(ia3+1,ia3+3,igc,jgc,i,j)=
     $                pspic(ia3+1,ia3+3,igc,jgc,i,j)+
     $                 pspi(ia3+1,ia3+3,ip,jp,i,j)*coefij
                      pspic(ia3+2,ia3+1,igc,jgc,i,j)=
     $                pspic(ia3+2,ia3+1,igc,jgc,i,j)+
     $                 pspi(ia3+2,ia3+1,ip,jp,i,j)*coefij
                      pspic(ia3+2,ia3+2,igc,jgc,i,j)=
     $                pspic(ia3+2,ia3+2,igc,jgc,i,j)+
     $                 pspi(ia3+2,ia3+2,ip,jp,i,j)*coefij
                      pspic(ia3+2,ia3+3,igc,jgc,i,j)=
     $                pspic(ia3+2,ia3+3,igc,jgc,i,j)+
     $                 pspi(ia3+2,ia3+3,ip,jp,i,j)*coefij
                      pspic(ia3+3,ia3+1,igc,jgc,i,j)=
     $                pspic(ia3+3,ia3+1,igc,jgc,i,j)+
     $                 pspi(ia3+3,ia3+1,ip,jp,i,j)*coefij
                      pspic(ia3+3,ia3+2,igc,jgc,i,j)=
     $                pspic(ia3+3,ia3+2,igc,jgc,i,j)+
     $                 pspi(ia3+3,ia3+2,ip,jp,i,j)*coefij
                      pspic(ia3+3,ia3+3,igc,jgc,i,j)=
     $                pspic(ia3+3,ia3+3,igc,jgc,i,j)+
     $                 pspi(ia3+3,ia3+3,ip,jp,i,j)*coefij
                      do jat=1,nderij
                        if(atder(jat))then
                          ja3=3*(jat-1)
                          pspic(ia3+1,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+1,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+1,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+1,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+1,ja3+3,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+2,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+2,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+2,ja3+3,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+1,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+1,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+1,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+2,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+2,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+2,ip,jp,i,j)*coefij
                          pspic(ia3+3,ja3+3,igc,jgc,i,j)=
     $                    pspic(ia3+3,ja3+3,igc,jgc,i,j)+
     $                     pspi(ia3+3,ja3+3,ip,jp,i,j)*coefij
c...
                          pspic(ja3+1,ia3+1,igc,jgc,i,j)=
     $                    pspic(ja3+1,ia3+1,igc,jgc,i,j)+
     $                     pspi(ja3+1,ia3+1,ip,jp,i,j)*coefij
                          pspic(ja3+1,ia3+2,igc,jgc,i,j)=
     $                    pspic(ja3+1,ia3+2,igc,jgc,i,j)+
     $                     pspi(ja3+1,ia3+2,ip,jp,i,j)*coefij
                          pspic(ja3+1,ia3+3,igc,jgc,i,j)=
     $                    pspic(ja3+1,ia3+3,igc,jgc,i,j)+
     $                     pspi(ja3+1,ia3+3,ip,jp,i,j)*coefij
                          pspic(ja3+2,ia3+1,igc,jgc,i,j)=
     $                    pspic(ja3+2,ia3+1,igc,jgc,i,j)+
     $                     pspi(ja3+2,ia3+1,ip,jp,i,j)*coefij
                          pspic(ja3+2,ia3+2,igc,jgc,i,j)=
     $                    pspic(ja3+2,ia3+2,igc,jgc,i,j)+
     $                     pspi(ja3+2,ia3+2,ip,jp,i,j)*coefij
                          pspic(ja3+2,ia3+3,igc,jgc,i,j)=
     $                    pspic(ja3+2,ia3+3,igc,jgc,i,j)+
     $                     pspi(ja3+2,ia3+3,ip,jp,i,j)*coefij
                          pspic(ja3+3,ia3+1,igc,jgc,i,j)=
     $                    pspic(ja3+3,ia3+1,igc,jgc,i,j)+
     $                     pspi(ja3+3,ia3+1,ip,jp,i,j)*coefij
                          pspic(ja3+3,ia3+2,igc,jgc,i,j)=
     $                    pspic(ja3+3,ia3+2,igc,jgc,i,j)+
     $                     pspi(ja3+3,ia3+2,ip,jp,i,j)*coefij
                          pspic(ja3+3,ia3+3,igc,jgc,i,j)=
     $                    pspic(ja3+3,ia3+3,igc,jgc,i,j)+
     $                     pspi(ja3+3,ia3+3,ip,jp,i,j)*coefij
                        endif
                      enddo
                    endif
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_hess12i(
     $                    basdat,   lla,   mma,   nna,    llb,
     $                       mmb,   nnb,   nga,   ngb,   lpsp,
     $                     nrpsp, alpha,    cc,  len0,    len,
     $                       cax,   cay,   caz,   cbx,    cby,
     $                       cbz,  lima,  limb,  mxpw,  nmaxa,
     $                     nmaxb,  mijk,   yra,   yrb,   xyzi,
     $                         q,    ar,  nder, iatom,  jatom,
     $                    icatom,   tol,argmax,  ilen,    ifi,
     $                       ca2,  jlen,   jfi,   cb2,    psp,
     $                     atder,  notc,  time)
      implicit real*8 (a-h,o-z)
c...
c...   MM august 2004
c...
c...   pseudopotential integral second derivatives over primitives
c...
c...  The integrals are differentiated with respect the centers
c...  of Gaussians A and B. The derivative with respect
c...  the pseudopotential center are obtained by translational
c...  invariance
c...
c...   Input parameters
c...
c...   basdat           basis set data
c...   alpa,alphb       exponents of Gaussians A and B
c...   lla,mma,nna,     arrays of angular momenta for A
c...   llb,mmb,nnb      arrays of angular momenta for B
c...   nga,ngb          number of Gaussian components for A and B
c...                    1=S, 3=P, 4=L, 6=D, 10=F, ...
c...   lpsp             maximum angular momentum of the pseudopotential
c...   nrpsp            array of powers of r of the psudopotential
c...   alpha            array of exponents of the pseudopotential
c...   cc               array of contraction coefficients of psp
c...   len0             contraction length of local part of psp
c...   len              array of contraction length of psp components
c...   cax,cay,caz      powers of CA
c...   cbx,cby,cbz      powers of CB
c...   lima,limb        do loop limits
c...   mxpw             maximum power of CA or CB
c...   nmaxa,nmaxb      maximum values of na and nb
c...   mijk             maximum exponent in angular integrals
c...   yra,yrb          scratch areas for real spherical harmonics
c...   xyzi             table of precomputed values for angular
c...                    integrals
c...   q                scratch area for radial integrals
c...   ar               scratch area for sums of angular and
c...                    radial integrals
c...   nder             number of atomic derivatives to be computed
c...   iatom            center of Gaussian A
c...   jatom            center of Gaussian B
c...   icatom           center of pseudopotential
c...   tol              tolerance threshold
c...   argmax           maximum argument of the exponential term
c...   ilen             contraction length of shell i
c...   ifi              pointer to first primitive of shell i
c...   ca2              distance CA squared
c...   jlen             contraction length of shell j
c...   jfi              pointer to first primitive of shell j
c...   cb2              distance CB squared
c...   notc             scratch area for logical handling of zeroing
c...
c...   Output prameter
c...
c...   atder  (input/output) flag for atomic derivatives
c...   psp              array of pseudopotential integral derivatives
c...
      dimension basdat(13,*)
      dimension notc(ilen,jlen)
      dimension psp(3*nder,3*nder,ilen,jlen,nga,ngb)
      dimension time(*)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      logical atder(*)
      logical compute,cumula,cumulb,cumulc,cumulab,cumulac,cumulbc,
     $        nocumul,csome
c...
c...  compute wheter integrals need to be accumulated or not
c...  (this is to avoid initializing the integrals over primitives)
c...
      csome=.false.
      cumula=atder(iatom)
      if(jatom.eq.iatom)then
        cumulb=.true.
      else
        cumulb=atder(jatom)
      endif
      if(icatom.eq.iatom.or.icatom.eq.jatom)then
        cumulc=.true.
      else
        cumulc=atder(icatom)
      endif
      if(jatom.eq.iatom)then
        cumulab=.true.
      else
        cumulab=cumula.and.cumulb
      endif
      if(icatom.eq.iatom.or.icatom.eq.jatom)then
        cumulac=.true.
      else
        cumulac=cumula.and.cumulc
      endif
      if(icatom.eq.jatom.or.icatom.eq.iatom.or.iatom.eq.jatom)then
        cumulbc=.true.
      else
        cumulbc=cumulb.and.cumulc
      endif
      nocumul=.not.cumula.or..not.cumulb.or..not.cumulc.or.
     $        .not.cumulab.or..not.cumulac.or..not.cumulbc
c...
c...  loop over primitives
c...
      maxn=nmaxa+nmaxb
      do ip=1,ilen
        ipi=ifi+ip
        alpa=basdat(1,ipi)     ! exponent of shell a
        do jp=1,jlen
          jpj=jfi+jp
          alpb=basdat(1,jpj)     ! exponent of shell b
          exparg=alpa*ca2+alpb*cb2
          notc(ip,jp)=0
          if(exparg.le.argmax)then ! test over the exponential
c...
c...  zero out ar array
c...
c           call secund(t0)
            call zeroit(ar,(nmaxa+1)**3*(nmaxb+1)**3)
            compute=.false.
c...
c...  compute radial and angular integrals
c...
            call psp_ar12(
     $                    alpa,  1.0d0,  alpb,  1.0d0,  lpsp,
     $                   nrpsp,  alpha,    cc,   len0,   len,
     $                     cax,    cay,   caz,    cbx,   cby,
     $                     cbz,   mxpw, nmaxa,  nmaxb,  maxn,
     $                    mijk,    yra,   yrb,   xyzi,     q,
     $                  exparg,   tol,     ar,compute,  time)
c           call secund(t1)
c           time(1)=time(1)+t1-t0
            if(compute)then
              notc(ip,jp)=1
              csome=.true.
c...
c...  assemble integrals over primitives
c...
            call psp_lg2(1,nga,1,ngb,alpa,alpb,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,nder,iatom,jatom,icatom,
     $                   cumula,cumulb,cumulc,cumulab,cumulac,cumulbc,
     $                   ip,jp,ilen,jlen,nga,ngb,psp)
            atder(iatom)=.true.
            atder(jatom)=.true.
            atder(icatom)=.true.
            endif
c           call secund(tifi)
c           time(7)=time(7)+tifi-t1
          endif
        enddo
      enddo
c...
c...  some initializing might be needed, after all ...
c...
      if(csome.and.nocumul)then
c           call secund(t1)
      ia3=3*(iatom-1)
      ja3=3*(jatom-1)
      ic3=3*(icatom-1)
      do ip=1,ilen
        do jp=1,jlen
          if(notc(ip,jp).eq.0)then
            do iga=1,nga
            do igb=1,ngb
              if(.not.cumula)then
                psp(ia3+1,ia3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ia3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ia3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ia3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ia3+1,ip,jp,iga,igb)=zero
                psp(ia3+2,ia3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ia3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ia3+2,ip,jp,iga,igb)=zero
                psp(ia3+3,ia3+3,ip,jp,iga,igb)=zero
              endif
              if(.not.cumulb)then
                psp(ja3+1,ja3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ja3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ja3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ja3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ja3+1,ip,jp,iga,igb)=zero
                psp(ja3+2,ja3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ja3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ja3+2,ip,jp,iga,igb)=zero
                psp(ja3+3,ja3+3,ip,jp,iga,igb)=zero
              endif
              if(.not.cumulc)then
                psp(ic3+1,ic3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ic3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ic3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ic3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ic3+1,ip,jp,iga,igb)=zero
                psp(ic3+2,ic3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ic3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ic3+2,ip,jp,iga,igb)=zero
                psp(ic3+3,ic3+3,ip,jp,iga,igb)=zero
              endif
              if(.not.cumulab)then
                psp(ia3+1,ja3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ja3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ja3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ja3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ja3+1,ip,jp,iga,igb)=zero
                psp(ia3+2,ja3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ja3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ja3+2,ip,jp,iga,igb)=zero
                psp(ia3+3,ja3+3,ip,jp,iga,igb)=zero
                psp(ja3+1,ia3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ia3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ia3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ia3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ia3+1,ip,jp,iga,igb)=zero
                psp(ja3+2,ia3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ia3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ia3+2,ip,jp,iga,igb)=zero
                psp(ja3+3,ia3+3,ip,jp,iga,igb)=zero
              endif
              if(.not.cumulac)then
                psp(ia3+1,ic3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ic3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ic3+1,ip,jp,iga,igb)=zero
                psp(ia3+1,ic3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ic3+1,ip,jp,iga,igb)=zero
                psp(ia3+2,ic3+2,ip,jp,iga,igb)=zero
                psp(ia3+2,ic3+3,ip,jp,iga,igb)=zero
                psp(ia3+3,ic3+2,ip,jp,iga,igb)=zero
                psp(ia3+3,ic3+3,ip,jp,iga,igb)=zero
                psp(ic3+1,ia3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ia3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ia3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ia3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ia3+1,ip,jp,iga,igb)=zero
                psp(ic3+2,ia3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ia3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ia3+2,ip,jp,iga,igb)=zero
                psp(ic3+3,ia3+3,ip,jp,iga,igb)=zero
              endif
              if(.not.cumulbc)then
                psp(ja3+1,ic3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ic3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ic3+1,ip,jp,iga,igb)=zero
                psp(ja3+1,ic3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ic3+1,ip,jp,iga,igb)=zero
                psp(ja3+2,ic3+2,ip,jp,iga,igb)=zero
                psp(ja3+2,ic3+3,ip,jp,iga,igb)=zero
                psp(ja3+3,ic3+2,ip,jp,iga,igb)=zero
                psp(ja3+3,ic3+3,ip,jp,iga,igb)=zero
                psp(ic3+1,ja3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ja3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ja3+1,ip,jp,iga,igb)=zero
                psp(ic3+1,ja3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ja3+1,ip,jp,iga,igb)=zero
                psp(ic3+2,ja3+2,ip,jp,iga,igb)=zero
                psp(ic3+2,ja3+3,ip,jp,iga,igb)=zero
                psp(ic3+3,ja3+2,ip,jp,iga,igb)=zero
                psp(ic3+3,ja3+3,ip,jp,iga,igb)=zero
              endif
            enddo
            enddo
          endif
        enddo
      enddo
c           call secund(tifi)
c           time(7)=time(7)+tifi-t1
      endif
      end
c=======================================================================
      subroutine psp_lg2(igas,igaf,igbs,igbf,alpa,alpb,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,nder,iat,jat,icat,
     $                   cumula,cumulb,cumulc,cumulab,cumulac,cumulbc,
     $                   i1,i2,nd1,nd2,nga,ngb,psp)
      implicit real*8 (a-h,o-z)
c...
c...  loops over the Gaussian components and assembling of the
c...  second-order gaussian components
c...
      dimension lla(nga),mma(nga),nna(nga),llb(ngb),mmb(ngb),nnb(ngb)
      dimension lima(3),limb(3)
      dimension psp(3*nder,3*nder,nd1,nd2,nga,ngb)
      logical cumula,cumulb,cumulc,cumulab,cumulac,cumulbc
c...
c...  variables for summation loops
c...
      logical doa,doam1,dob,dobm1,doc,docm1
      logical dod,dodm1,doe,doem1,dof,dofm1
      integer al,af,am1l,am1f,am2l,am2f,ap1l,ap1f,ap2l,ap2f
      integer bl,bf,bm1l,bm1f,bm2l,bm2f,bp1l,bp1f,bp2l,bp2f
      integer cl,cf,cm1l,cm1f,cm2l,cm2f,cp1l,cp1f,cp2l,cp2f
      integer dl,df,dm1l,dm1f,dm2l,dm2f,dp1l,dp1f,dp2l,dp2f
      integer el,ef,em1l,em1f,em2l,em2f,ep1l,ep1f,ep2l,ep2f
      integer fl,ff,fm1l,fm1f,fm2l,fm2f,fp1l,fp1f,fp2l,fp2f
c...
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
c...
      ia3=3*(iat-1)
      ja3=3*(jat-1)
      ic3=3*(icat-1)
      tam=-two*alpa
      tam2=tam*tam
      tbm=-two*alpb
      tbm2=tbm*tbm
      do iga=igas,igaf
        al=lla(iga)
        af=min(al,lima(1))
        fal=dfloat(al)
        doa=al.gt.0
        am1l=al-1
        am1f=min(am1l,lima(1))
        fam1l=dfloat(am1l)
        doam1=am1l.gt.0
        am2l=al-2
        am2f=min(am2l,lima(1))
        ap1l=al+1
        ap1f=min(ap1l,lima(1))
        ap2l=al+2
        ap2f=min(ap2l,lima(1))
c...
        bl=mma(iga)
        bf=min(bl,lima(2))
        fbl=dfloat(bl)
        dob=bl.gt.0
        bm1l=bl-1
        bm1f=min(bm1l,lima(2))
        fbm1l=dfloat(bm1l)
        dobm1=bm1l.gt.0
        bm2l=bl-2
        bm2f=min(bm2l,lima(2))
        bp1l=bl+1
        bp1f=min(bp1l,lima(2))
        bp2l=bl+2
        bp2f=min(bp2l,lima(2))
c...
        cl=nna(iga)
        cf=min(cl,lima(3))
        fcl=dfloat(cl)
        doc=cl.gt.0
        cm1l=cl-1
        cm1f=min(cm1l,lima(3))
        fcm1l=dfloat(cm1l)
        docm1=cm1l.gt.0
        cm2l=cl-2
        cm2f=min(cm2l,lima(3))
        cp1l=cl+1
        cp1f=min(cp1l,lima(3))
        cp2l=cl+2
        cp2f=min(cp2l,lima(3))
        do igb=igbs,igbf
          dl=llb(igb)
          df=min(dl,limb(1))
          fdl=dfloat(dl)
          dod=dl.gt.0
          dm1l=dl-1
          dm1f=min(dm1l,limb(1))
          fdm1l=dfloat(dm1l)
          dodm1=dm1l.gt.0
          dm2l=dl-2
          dm2f=min(dm2l,limb(1))
          dp1l=dl+1
          dp1f=min(dp1l,limb(1))
          dp2l=dl+2
          dp2f=min(dp2l,limb(1))
c...
          el=mmb(igb)
          ef=min(el,limb(2))
          fel=dfloat(el)
          doe=el.gt.0
          em1l=el-1
          em1f=min(em1l,limb(2))
          fem1l=dfloat(em1l)
          doem1=em1l.gt.0
          em2l=el-2
          em2f=min(em2l,limb(2))
          ep1l=el+1
          ep1f=min(ep1l,limb(2))
          ep2l=el+2
          ep2f=min(ep2l,limb(2))
c...
          fl=nnb(igb)
          ff=min(fl,limb(3))
          ffl=dfloat(fl)
          dof=fl.gt.0
          fm1l=fl-1
          fm1f=min(fm1l,limb(3))
          ffm1l=dfloat(fm1l)
          dofm1=fm1l.gt.0
          fm2l=fl-2
          fm2f=min(fm2l,limb(3))
          fp1l=fl+1
          fp1f=min(fp1l,limb(3))
          fp2l=fl+2
          fp2f=min(fp2l,limb(3))
c...
c...  summations
c...
c...  unshifted
c...
          call psp_suml2(al,af,bl,bf,cl,cf,dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppint)
c...
c...  A  xx
c...
          call psp_suml2(ap2l,ap2f,bl,bf,cl,cf,dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppaxp2)
          ppaxm2=zero
          if(doam1)call psp_suml2(am2l,am2f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppaxm2)
c...
c...  A  xy
c...
          call psp_suml2(ap1l,ap1f,bp1l,bp1f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1yp1)
          ppaxp1ym1=zero
          if(dob)call psp_suml2(ap1l,ap1f,bm1l,bm1f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1ym1)
          ppaxm1ym1=zero
          if(doa.and.dob)call psp_suml2(am1l,am1f,bm1l,bm1f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1ym1)
          ppaxm1yp1=zero
          if(doa)call psp_suml2(am1l,am1f,bp1l,bp1f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1yp1)
c...
c...  A  xz
c...
          call psp_suml2(ap1l,ap1f,bl,bf,cp1l,cp1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1zp1)
          ppaxp1zm1=zero
          if(doc)call psp_suml2(ap1l,ap1f,bl,bf,cm1l,cm1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1zm1)
          ppaxm1zm1=zero
          if(doa.and.doc)call psp_suml2(am1l,am1f,bl,bf,cm1l,cm1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1zm1)
          ppaxm1zp1=zero
          if(doa)call psp_suml2(am1l,am1f,bl,bf,cp1l,cp1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1zp1)
c...
c...  A  yy
c...
          call psp_suml2(al,af,bp2l,bp2f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppayp2)
          if(dobm1)then
          call psp_suml2(al,af,bm2l,bm2f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppaym2)
          else
            ppaym2=zero
          endif
c...
c...  A  yz
c...
          call psp_suml2(al,af,bp1l,bp1f,cp1l,cp1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1zp1)
          ppayp1zm1=zero
          if(doc)call psp_suml2(al,af,bp1l,bp1f,cm1l,cm1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1zm1)
          ppaym1zm1=zero
          if(dob.and.doc)call psp_suml2(al,af,bm1l,bm1f,cm1l,cm1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1zm1)
          ppaym1zp1=zero
          if(dob)call psp_suml2(al,af,bm1l,bm1f,cp1l,cp1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1zp1)
c...
c...  A  zz
c...
          call psp_suml2(al,af,bl,bf,cp2l,cp2f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppazp2)
          if(docm1)then
          call psp_suml2(al,af,bl,bf,cm2l,cm2f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppazm2)
          else
            ppazm2=zero
          endif
c...
c...  B  xx
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp2l,dp2f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxp2)
          if(dodm1)then
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dm2l,dm2f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxm2)
          else
            ppbxm2=zero
          endif
c...
c...  B  xy
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp1l,dp1f,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxp1yp1)
          ppbxp1ym1=zero
          if(doe)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp1l,dp1f,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxp1ym1)
          ppbxm1ym1=zero
          if(dod.and.doe)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dm1l,dm1f,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxm1ym1)
          ppbxm1yp1=zero
          if(dod)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dm1l,dm1f,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxm1yp1)
c...
c...  B  xz
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp1l,dp1f,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxp1zp1)
          ppbxp1zm1=zero
          if(dof)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp1l,dp1f,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxp1zm1)
          ppbxm1zm1=zero
          if(dod.and.dof)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dm1l,dm1f,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxm1zm1)
          ppbxm1zp1=zero
          if(dod)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dm1l,dm1f,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbxm1zp1)
c...
c...  B  yy
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,ep2l,ep2f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbyp2)
          if(doem1)then
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,em2l,ep2f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbym2)
          else
            ppbym2=zero
          endif
c...
c...  B  yz
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,ep1l,ep1f,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbyp1zp1)
          ppbyp1zm1=zero
          if(dof)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,ep1l,ep1f,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbyp1zm1)
          ppbym1zm1=zero
          if(doe.and.dof)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,em1l,em1f,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbym1zm1)
          ppbym1zp1=zero
          if(doe)call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,em1l,em1f,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppbym1zp1)
c...
c...  B  zz
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,el,ef,fp2l,fp2f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbzp2)
          if(dofm1)then
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,el,ef,fm2l,fm2f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbzm2)
          else
            ppbzm2=zero
          endif
c...
c...  A x  B x
c...
          call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1bxp1)
          ppaxp1bxm1=zero
          if(dod)call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1bxm1)
          ppaxm1bxm1=zero
          if(doa.and.dod)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1bxm1)
          ppaxm1bxp1=zero
          if(doa)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1bxp1)
c...
c...  A x  B y
c...
          call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1byp1)
          ppaxp1bym1=zero
          if(doe)call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1bym1)
          ppaxm1bym1=zero
          if(doa.and.doe)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1bym1)
          ppaxm1byp1=zero
          if(doa)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1byp1)
c...
c...  A x  B z
c...
          call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1bzp1)
          ppaxp1bzm1=zero
          if(dof)call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxp1bzm1)
          ppaxm1bzm1=zero
          if(doa.and.dof)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1bzm1)
          ppaxm1bzp1=zero
          if(doa)call psp_suml2(am1l,am1f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaxm1bzp1)
c...
c...  A y  B x
c...
          call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1bxp1)
          ppayp1bxm1=zero
          if(dod)call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1bxm1)
          ppaym1bxm1=zero
          if(dob.and.dod)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1bxm1)
          ppaym1bxp1=zero
          if(dob)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1bxp1)
c...
c...  A y  B y
c...
          call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1byp1)
          ppayp1bym1=zero
          if(doe)call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1bym1)
          ppaym1bym1=zero
          if(dob.and.doe)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1bym1)
          ppaym1byp1=zero
          if(dob)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1byp1)
c...
c...  A y  B z
c...
          call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1bzp1)
          ppayp1bzm1=zero
          if(dof)call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppayp1bzm1)
          ppaym1bzm1=zero
          if(dob.and.dof)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1bzm1)
          ppaym1bzp1=zero
          if(dob)call psp_suml2(al,af,bm1l,bm1f,cl,cf,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppaym1bzp1)
c...
c...  A z  B x
c...
          call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1bxp1)
          ppazp1bxm1=zero
          if(dod)call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1bxm1)
          ppazm1bxm1=zero
          if(doc.and.dod)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dm1l,dm1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1bxm1)
          ppazm1bxp1=zero
          if(doc)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1bxp1)
c...
c...  A z  B y
c...
          call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1byp1)
          ppazp1bym1=zero
          if(doe)call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1bym1)
          ppazm1bym1=zero
          if(doc.and.doe)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dl,df,em1l,em1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1bym1)
          ppazm1byp1=zero
          if(doc)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1byp1)
c...
c...  A z  B z
c...
          call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1bzp1)
          ppazp1bzm1=zero
          if(dof)call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazp1bzm1)
          ppazm1bzm1=zero
          if(doc.and.dof)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dl,df,el,ef,fm1l,fm1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1bzm1)
          ppazm1bzp1=zero
          if(doc)call psp_suml2(al,af,bl,bf,cm1l,cm1f,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,
     $                   ppazm1bzp1)
c...
c...  assemble second derivatives
c...
c...
c...  A  A
c...
          axax=tam2*ppaxp2+tam*(two*fal+one)*ppint+fal*fam1l*ppaxm2
          axay=tam2*ppaxp1yp1+tam*(fal*ppaxm1yp1+fbl*ppaxp1ym1)+
     $    fal*fbl*ppaxm1ym1
          axaz=tam2*ppaxp1zp1+tam*(fal*ppaxm1zp1+fcl*ppaxp1zm1)+
     $    fal*fcl*ppaxm1zm1
          ayay=tam2*ppayp2+tam*(two*fbl+one)*ppint+fbl*fbm1l*ppaym2
          ayaz=tam2*ppayp1zp1+tam*(fbl*ppaym1zp1+fcl*ppayp1zm1)+
     $    fbl*fcl*ppaym1zm1
          azaz=tam2*ppazp2+tam*(two*fcl+one)*ppint+fcl*fcm1l*ppazm2
c...
c...  B  B
c...
          bxbx=tbm2*ppbxp2+tbm*(two*fdl+one)*ppint+fdl*fdm1l*ppbxm2
          bxby=tbm2*ppbxp1yp1+tbm*(fdl*ppbxm1yp1+fel*ppbxp1ym1)+
     $    fdl*fel*ppbxm1ym1
          bxbz=tbm2*ppbxp1zp1+tbm*(fdl*ppbxm1zp1+ffl*ppbxp1zm1)+
     $    fdl*ffl*ppbxm1zm1
          byby=tbm2*ppbyp2+tbm*(two*fel+one)*ppint+fel*fem1l*ppbym2
          bybz=tbm2*ppbyp1zp1+tbm*(fel*ppbym1zp1+ffl*ppbyp1zm1)+
     $    fel*ffl*ppbym1zm1
          bzbz=tbm2*ppbzp2+tbm*(two*ffl+one)*ppint+ffl*ffm1l*ppbzm2
c...
c...  A  B
c...
          axbx=tam*tbm*ppaxp1bxp1+fal*fdl*ppaxm1bxm1+
     $         fal*tbm*ppaxm1bxp1+tam*fdl*ppaxp1bxm1
          axby=tam*tbm*ppaxp1byp1+fal*fel*ppaxm1bym1+
     $         fal*tbm*ppaxm1byp1+tam*fel*ppaxp1bym1
          axbz=tam*tbm*ppaxp1bzp1+fal*ffl*ppaxm1bzm1+
     $         fal*tbm*ppaxm1bzp1+tam*ffl*ppaxp1bzm1
          aybx=tam*tbm*ppayp1bxp1+fbl*fdl*ppaym1bxm1+
     $         fbl*tbm*ppaym1bxp1+tam*fdl*ppayp1bxm1
          ayby=tam*tbm*ppayp1byp1+fbl*fel*ppaym1bym1+
     $         fbl*tbm*ppaym1byp1+tam*fel*ppayp1bym1
          aybz=tam*tbm*ppayp1bzp1+fbl*ffl*ppaym1bzm1+
     $         fbl*tbm*ppaym1bzp1+tam*ffl*ppayp1bzm1
          azbx=tam*tbm*ppazp1bxp1+fcl*fdl*ppazm1bxm1+
     $         fcl*tbm*ppazm1bxp1+tam*fdl*ppazp1bxm1
          azby=tam*tbm*ppazp1byp1+fcl*fel*ppazm1bym1+
     $         fcl*tbm*ppazm1byp1+tam*fel*ppazp1bym1
          azbz=tam*tbm*ppazp1bzp1+fcl*ffl*ppazm1bzm1+
     $         fcl*tbm*ppazm1bzp1+tam*ffl*ppazp1bzm1
c...
c...  use translational invariance to get derivatives on C
c...
c...   A  C
c...
          axcx=-axax-axbx
          axcy=-axay-axby
          axcz=-axaz-axbz
          aycx=-axay-aybx
          aycy=-ayay-ayby
          aycz=-ayaz-aybz
          azcx=-axaz-azbx
          azcy=-ayaz-azby
          azcz=-azaz-azbz
c...
c...   B  C
c...
          bxcx=-bxbx-axbx
          bxcy=-bxby-aybx
          bxcz=-bxbz-azbx
          bycx=-bxby-axby
          bycy=-byby-ayby
          bycz=-bybz-azby
          bzcx=-bxbz-axbz
          bzcy=-bybz-aybz
          bzcz=-bzbz-azbz
c...
c...   C  C
c...
          cxcx=-axcx-bxcx
          cxcy=-axcy-bxcy
          cxcz=-axcz-bxcz
          cycy=-aycy-bycy
          cycz=-aycz-bycz
          czcz=-azcz-bzcz
c...
c...  fill the psp1 array
c...
c...
c...   A  A
c...
      if(cumula)then
      psp(ia3+1,ia3+1,i1,i2,iga,igb)=psp(ia3+1,ia3+1,i1,i2,iga,igb)+axax
      psp(ia3+1,ia3+2,i1,i2,iga,igb)=psp(ia3+1,ia3+2,i1,i2,iga,igb)+axay
      psp(ia3+2,ia3+1,i1,i2,iga,igb)=psp(ia3+2,ia3+1,i1,i2,iga,igb)+axay
      psp(ia3+1,ia3+3,i1,i2,iga,igb)=psp(ia3+1,ia3+3,i1,i2,iga,igb)+axaz
      psp(ia3+3,ia3+1,i1,i2,iga,igb)=psp(ia3+3,ia3+1,i1,i2,iga,igb)+axaz
      psp(ia3+2,ia3+2,i1,i2,iga,igb)=psp(ia3+2,ia3+2,i1,i2,iga,igb)+ayay
      psp(ia3+2,ia3+3,i1,i2,iga,igb)=psp(ia3+2,ia3+3,i1,i2,iga,igb)+ayaz
      psp(ia3+3,ia3+2,i1,i2,iga,igb)=psp(ia3+3,ia3+2,i1,i2,iga,igb)+ayaz
      psp(ia3+3,ia3+3,i1,i2,iga,igb)=psp(ia3+3,ia3+3,i1,i2,iga,igb)+azaz
      else
      psp(ia3+1,ia3+1,i1,i2,iga,igb)=axax
      psp(ia3+1,ia3+2,i1,i2,iga,igb)=axay
      psp(ia3+2,ia3+1,i1,i2,iga,igb)=axay
      psp(ia3+1,ia3+3,i1,i2,iga,igb)=axaz
      psp(ia3+3,ia3+1,i1,i2,iga,igb)=axaz
      psp(ia3+2,ia3+2,i1,i2,iga,igb)=ayay
      psp(ia3+2,ia3+3,i1,i2,iga,igb)=ayaz
      psp(ia3+3,ia3+2,i1,i2,iga,igb)=ayaz
      psp(ia3+3,ia3+3,i1,i2,iga,igb)=azaz
      endif
c...
c...   B  B
c...
      if(cumulb)then
      psp(ja3+1,ja3+1,i1,i2,iga,igb)=psp(ja3+1,ja3+1,i1,i2,iga,igb)+bxbx
      psp(ja3+1,ja3+2,i1,i2,iga,igb)=psp(ja3+1,ja3+2,i1,i2,iga,igb)+bxby
      psp(ja3+2,ja3+1,i1,i2,iga,igb)=psp(ja3+2,ja3+1,i1,i2,iga,igb)+bxby
      psp(ja3+1,ja3+3,i1,i2,iga,igb)=psp(ja3+1,ja3+3,i1,i2,iga,igb)+bxbz
      psp(ja3+3,ja3+1,i1,i2,iga,igb)=psp(ja3+3,ja3+1,i1,i2,iga,igb)+bxbz
      psp(ja3+2,ja3+2,i1,i2,iga,igb)=psp(ja3+2,ja3+2,i1,i2,iga,igb)+byby
      psp(ja3+2,ja3+3,i1,i2,iga,igb)=psp(ja3+2,ja3+3,i1,i2,iga,igb)+bybz
      psp(ja3+3,ja3+2,i1,i2,iga,igb)=psp(ja3+3,ja3+2,i1,i2,iga,igb)+bybz
      psp(ja3+3,ja3+3,i1,i2,iga,igb)=psp(ja3+3,ja3+3,i1,i2,iga,igb)+bzbz
      else
      psp(ja3+1,ja3+1,i1,i2,iga,igb)=bxbx
      psp(ja3+1,ja3+2,i1,i2,iga,igb)=bxby
      psp(ja3+2,ja3+1,i1,i2,iga,igb)=bxby
      psp(ja3+1,ja3+3,i1,i2,iga,igb)=bxbz
      psp(ja3+3,ja3+1,i1,i2,iga,igb)=bxbz
      psp(ja3+2,ja3+2,i1,i2,iga,igb)=byby
      psp(ja3+2,ja3+3,i1,i2,iga,igb)=bybz
      psp(ja3+3,ja3+2,i1,i2,iga,igb)=bybz
      psp(ja3+3,ja3+3,i1,i2,iga,igb)=bzbz
      endif
c...
c...   C  C
c...
      if(cumulc)then
      psp(ic3+1,ic3+1,i1,i2,iga,igb)=psp(ic3+1,ic3+1,i1,i2,iga,igb)+cxcx
      psp(ic3+1,ic3+2,i1,i2,iga,igb)=psp(ic3+1,ic3+2,i1,i2,iga,igb)+cxcy
      psp(ic3+2,ic3+1,i1,i2,iga,igb)=psp(ic3+2,ic3+1,i1,i2,iga,igb)+cxcy
      psp(ic3+1,ic3+3,i1,i2,iga,igb)=psp(ic3+1,ic3+3,i1,i2,iga,igb)+cxcz
      psp(ic3+3,ic3+1,i1,i2,iga,igb)=psp(ic3+3,ic3+1,i1,i2,iga,igb)+cxcz
      psp(ic3+2,ic3+2,i1,i2,iga,igb)=psp(ic3+2,ic3+2,i1,i2,iga,igb)+cycy
      psp(ic3+2,ic3+3,i1,i2,iga,igb)=psp(ic3+2,ic3+3,i1,i2,iga,igb)+cycz
      psp(ic3+3,ic3+2,i1,i2,iga,igb)=psp(ic3+3,ic3+2,i1,i2,iga,igb)+cycz
      psp(ic3+3,ic3+3,i1,i2,iga,igb)=psp(ic3+3,ic3+3,i1,i2,iga,igb)+czcz
      else
      psp(ic3+1,ic3+1,i1,i2,iga,igb)=cxcx
      psp(ic3+1,ic3+2,i1,i2,iga,igb)=cxcy
      psp(ic3+2,ic3+1,i1,i2,iga,igb)=cxcy
      psp(ic3+1,ic3+3,i1,i2,iga,igb)=cxcz
      psp(ic3+3,ic3+1,i1,i2,iga,igb)=cxcz
      psp(ic3+2,ic3+2,i1,i2,iga,igb)=cycy
      psp(ic3+2,ic3+3,i1,i2,iga,igb)=cycz
      psp(ic3+3,ic3+2,i1,i2,iga,igb)=cycz
      psp(ic3+3,ic3+3,i1,i2,iga,igb)=czcz
      endif
c...
c...   A  B
c...
      if(cumulab)then
      psp(ia3+1,ja3+1,i1,i2,iga,igb)=psp(ia3+1,ja3+1,i1,i2,iga,igb)+axbx
      psp(ia3+1,ja3+2,i1,i2,iga,igb)=psp(ia3+1,ja3+2,i1,i2,iga,igb)+axby
      psp(ia3+1,ja3+3,i1,i2,iga,igb)=psp(ia3+1,ja3+3,i1,i2,iga,igb)+axbz
      psp(ia3+2,ja3+1,i1,i2,iga,igb)=psp(ia3+2,ja3+1,i1,i2,iga,igb)+aybx
      psp(ia3+2,ja3+2,i1,i2,iga,igb)=psp(ia3+2,ja3+2,i1,i2,iga,igb)+ayby
      psp(ia3+2,ja3+3,i1,i2,iga,igb)=psp(ia3+2,ja3+3,i1,i2,iga,igb)+aybz
      psp(ia3+3,ja3+1,i1,i2,iga,igb)=psp(ia3+3,ja3+1,i1,i2,iga,igb)+azbx
      psp(ia3+3,ja3+2,i1,i2,iga,igb)=psp(ia3+3,ja3+2,i1,i2,iga,igb)+azby
      psp(ia3+3,ja3+3,i1,i2,iga,igb)=psp(ia3+3,ja3+3,i1,i2,iga,igb)+azbz
      psp(ja3+1,ia3+1,i1,i2,iga,igb)=psp(ja3+1,ia3+1,i1,i2,iga,igb)+axbx
      psp(ja3+2,ia3+1,i1,i2,iga,igb)=psp(ja3+2,ia3+1,i1,i2,iga,igb)+axby
      psp(ja3+3,ia3+1,i1,i2,iga,igb)=psp(ja3+3,ia3+1,i1,i2,iga,igb)+axbz
      psp(ja3+1,ia3+2,i1,i2,iga,igb)=psp(ja3+1,ia3+2,i1,i2,iga,igb)+aybx
      psp(ja3+2,ia3+2,i1,i2,iga,igb)=psp(ja3+2,ia3+2,i1,i2,iga,igb)+ayby
      psp(ja3+3,ia3+2,i1,i2,iga,igb)=psp(ja3+3,ia3+2,i1,i2,iga,igb)+aybz
      psp(ja3+1,ia3+3,i1,i2,iga,igb)=psp(ja3+1,ia3+3,i1,i2,iga,igb)+azbx
      psp(ja3+2,ia3+3,i1,i2,iga,igb)=psp(ja3+2,ia3+3,i1,i2,iga,igb)+azby
      psp(ja3+3,ia3+3,i1,i2,iga,igb)=psp(ja3+3,ia3+3,i1,i2,iga,igb)+azbz
      else
      psp(ia3+1,ja3+1,i1,i2,iga,igb)=axbx
      psp(ia3+1,ja3+2,i1,i2,iga,igb)=axby
      psp(ia3+1,ja3+3,i1,i2,iga,igb)=axbz
      psp(ia3+2,ja3+1,i1,i2,iga,igb)=aybx
      psp(ia3+2,ja3+2,i1,i2,iga,igb)=ayby
      psp(ia3+2,ja3+3,i1,i2,iga,igb)=aybz
      psp(ia3+3,ja3+1,i1,i2,iga,igb)=azbx
      psp(ia3+3,ja3+2,i1,i2,iga,igb)=azby
      psp(ia3+3,ja3+3,i1,i2,iga,igb)=azbz
      psp(ja3+1,ia3+1,i1,i2,iga,igb)=axbx
      psp(ja3+2,ia3+1,i1,i2,iga,igb)=axby
      psp(ja3+3,ia3+1,i1,i2,iga,igb)=axbz
      psp(ja3+1,ia3+2,i1,i2,iga,igb)=aybx
      psp(ja3+2,ia3+2,i1,i2,iga,igb)=ayby
      psp(ja3+3,ia3+2,i1,i2,iga,igb)=aybz
      psp(ja3+1,ia3+3,i1,i2,iga,igb)=azbx
      psp(ja3+2,ia3+3,i1,i2,iga,igb)=azby
      psp(ja3+3,ia3+3,i1,i2,iga,igb)=azbz
      endif
c...
c...   A  C
c...
      if(cumulac)then
      psp(ia3+1,ic3+1,i1,i2,iga,igb)=psp(ia3+1,ic3+1,i1,i2,iga,igb)+axcx
      psp(ia3+1,ic3+2,i1,i2,iga,igb)=psp(ia3+1,ic3+2,i1,i2,iga,igb)+axcy
      psp(ia3+1,ic3+3,i1,i2,iga,igb)=psp(ia3+1,ic3+3,i1,i2,iga,igb)+axcz
      psp(ia3+2,ic3+1,i1,i2,iga,igb)=psp(ia3+2,ic3+1,i1,i2,iga,igb)+aycx
      psp(ia3+2,ic3+2,i1,i2,iga,igb)=psp(ia3+2,ic3+2,i1,i2,iga,igb)+aycy
      psp(ia3+2,ic3+3,i1,i2,iga,igb)=psp(ia3+2,ic3+3,i1,i2,iga,igb)+aycz
      psp(ia3+3,ic3+1,i1,i2,iga,igb)=psp(ia3+3,ic3+1,i1,i2,iga,igb)+azcx
      psp(ia3+3,ic3+2,i1,i2,iga,igb)=psp(ia3+3,ic3+2,i1,i2,iga,igb)+azcy
      psp(ia3+3,ic3+3,i1,i2,iga,igb)=psp(ia3+3,ic3+3,i1,i2,iga,igb)+azcz
      psp(ic3+1,ia3+1,i1,i2,iga,igb)=psp(ic3+1,ia3+1,i1,i2,iga,igb)+axcx
      psp(ic3+2,ia3+1,i1,i2,iga,igb)=psp(ic3+2,ia3+1,i1,i2,iga,igb)+axcy
      psp(ic3+3,ia3+1,i1,i2,iga,igb)=psp(ic3+3,ia3+1,i1,i2,iga,igb)+axcz
      psp(ic3+1,ia3+2,i1,i2,iga,igb)=psp(ic3+1,ia3+2,i1,i2,iga,igb)+aycx
      psp(ic3+2,ia3+2,i1,i2,iga,igb)=psp(ic3+2,ia3+2,i1,i2,iga,igb)+aycy
      psp(ic3+3,ia3+2,i1,i2,iga,igb)=psp(ic3+3,ia3+2,i1,i2,iga,igb)+aycz
      psp(ic3+1,ia3+3,i1,i2,iga,igb)=psp(ic3+1,ia3+3,i1,i2,iga,igb)+azcx
      psp(ic3+2,ia3+3,i1,i2,iga,igb)=psp(ic3+2,ia3+3,i1,i2,iga,igb)+azcy
      psp(ic3+3,ia3+3,i1,i2,iga,igb)=psp(ic3+3,ia3+3,i1,i2,iga,igb)+azcz
      else
      psp(ia3+1,ic3+1,i1,i2,iga,igb)=axcx
      psp(ia3+1,ic3+2,i1,i2,iga,igb)=axcy
      psp(ia3+1,ic3+3,i1,i2,iga,igb)=axcz
      psp(ia3+2,ic3+1,i1,i2,iga,igb)=aycx
      psp(ia3+2,ic3+2,i1,i2,iga,igb)=aycy
      psp(ia3+2,ic3+3,i1,i2,iga,igb)=aycz
      psp(ia3+3,ic3+1,i1,i2,iga,igb)=azcx
      psp(ia3+3,ic3+2,i1,i2,iga,igb)=azcy
      psp(ia3+3,ic3+3,i1,i2,iga,igb)=azcz
      psp(ic3+1,ia3+1,i1,i2,iga,igb)=axcx
      psp(ic3+2,ia3+1,i1,i2,iga,igb)=axcy
      psp(ic3+3,ia3+1,i1,i2,iga,igb)=axcz
      psp(ic3+1,ia3+2,i1,i2,iga,igb)=aycx
      psp(ic3+2,ia3+2,i1,i2,iga,igb)=aycy
      psp(ic3+3,ia3+2,i1,i2,iga,igb)=aycz
      psp(ic3+1,ia3+3,i1,i2,iga,igb)=azcx
      psp(ic3+2,ia3+3,i1,i2,iga,igb)=azcy
      psp(ic3+3,ia3+3,i1,i2,iga,igb)=azcz
      endif
c...
c...   B  C
c...
      if(cumulbc)then
      psp(ja3+1,ic3+1,i1,i2,iga,igb)=psp(ja3+1,ic3+1,i1,i2,iga,igb)+bxcx
      psp(ja3+1,ic3+2,i1,i2,iga,igb)=psp(ja3+1,ic3+2,i1,i2,iga,igb)+bxcy
      psp(ja3+1,ic3+3,i1,i2,iga,igb)=psp(ja3+1,ic3+3,i1,i2,iga,igb)+bxcz
      psp(ja3+2,ic3+1,i1,i2,iga,igb)=psp(ja3+2,ic3+1,i1,i2,iga,igb)+bycx
      psp(ja3+2,ic3+2,i1,i2,iga,igb)=psp(ja3+2,ic3+2,i1,i2,iga,igb)+bycy
      psp(ja3+2,ic3+3,i1,i2,iga,igb)=psp(ja3+2,ic3+3,i1,i2,iga,igb)+bycz
      psp(ja3+3,ic3+1,i1,i2,iga,igb)=psp(ja3+3,ic3+1,i1,i2,iga,igb)+bzcx
      psp(ja3+3,ic3+2,i1,i2,iga,igb)=psp(ja3+3,ic3+2,i1,i2,iga,igb)+bzcy
      psp(ja3+3,ic3+3,i1,i2,iga,igb)=psp(ja3+3,ic3+3,i1,i2,iga,igb)+bzcz
      psp(ic3+1,ja3+1,i1,i2,iga,igb)=psp(ic3+1,ja3+1,i1,i2,iga,igb)+bxcx
      psp(ic3+2,ja3+1,i1,i2,iga,igb)=psp(ic3+2,ja3+1,i1,i2,iga,igb)+bxcy
      psp(ic3+3,ja3+1,i1,i2,iga,igb)=psp(ic3+3,ja3+1,i1,i2,iga,igb)+bxcz
      psp(ic3+1,ja3+2,i1,i2,iga,igb)=psp(ic3+1,ja3+2,i1,i2,iga,igb)+bycx
      psp(ic3+2,ja3+2,i1,i2,iga,igb)=psp(ic3+2,ja3+2,i1,i2,iga,igb)+bycy
      psp(ic3+3,ja3+2,i1,i2,iga,igb)=psp(ic3+3,ja3+2,i1,i2,iga,igb)+bycz
      psp(ic3+1,ja3+3,i1,i2,iga,igb)=psp(ic3+1,ja3+3,i1,i2,iga,igb)+bzcx
      psp(ic3+2,ja3+3,i1,i2,iga,igb)=psp(ic3+2,ja3+3,i1,i2,iga,igb)+bzcy
      psp(ic3+3,ja3+3,i1,i2,iga,igb)=psp(ic3+3,ja3+3,i1,i2,iga,igb)+bzcz
      else
      psp(ja3+1,ic3+1,i1,i2,iga,igb)=bxcx
      psp(ja3+1,ic3+2,i1,i2,iga,igb)=bxcy
      psp(ja3+1,ic3+3,i1,i2,iga,igb)=bxcz
      psp(ja3+2,ic3+1,i1,i2,iga,igb)=bycx
      psp(ja3+2,ic3+2,i1,i2,iga,igb)=bycy
      psp(ja3+2,ic3+3,i1,i2,iga,igb)=bycz
      psp(ja3+3,ic3+1,i1,i2,iga,igb)=bzcx
      psp(ja3+3,ic3+2,i1,i2,iga,igb)=bzcy
      psp(ja3+3,ic3+3,i1,i2,iga,igb)=bzcz
      psp(ic3+1,ja3+1,i1,i2,iga,igb)=bxcx
      psp(ic3+2,ja3+1,i1,i2,iga,igb)=bxcy
      psp(ic3+3,ja3+1,i1,i2,iga,igb)=bxcz
      psp(ic3+1,ja3+2,i1,i2,iga,igb)=bycx
      psp(ic3+2,ja3+2,i1,i2,iga,igb)=bycy
      psp(ic3+3,ja3+2,i1,i2,iga,igb)=bycz
      psp(ic3+1,ja3+3,i1,i2,iga,igb)=bzcx
      psp(ic3+2,ja3+3,i1,i2,iga,igb)=bzcy
      psp(ic3+3,ja3+3,i1,i2,iga,igb)=bzcz
      endif
c...
        enddo
      enddo
c...
      end

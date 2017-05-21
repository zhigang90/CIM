c...
c...  MM September-October 2004
c...
c...  this file contains the driver routines for the calculation
c...  of analytical gradients over ab initio pseudopotentials.
c...
c======================================================================
      subroutine pspgrad(na,npsp,ncf,ncs,inx,xnuc,basdat,dens,tol,forc)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  driver routine for the psp contribution to the forces
c...
c...  na      number of atoms
c...  npsp    number of pseudopotentials
c...  ncf     number of contracted functions
c...  ncs     number of contracted shells
c...  inx     contraction information
c...  xnuc    nuclear data
c...  basdat  basis set data
c...  dens    density matrix
c...  tol     tolerance threshold
c...  forc    cumulative forces
c...
      dimension inx(12,*),basdat(13,*),xnuc(5,*),dens(*),forc(3,na)
c     common /big/bl(1)
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
c...  loop over contracted shells
c...
      maxpass=0
      do ics=1,ncs
        do jcs=1,ics
c...
c...  compute integrals
c...
          call psp_forc(   ics,       jcs,     npsp,    nradp,  maxlpsp,
     $                     ncs,       ncf,    maxlp,     maxl,     mijk,
     $            bl(iiecpdat),bl(inrrad),bl(iarad),bl(icrad), bl(ilen),
     $                    xnuc,       inx,   basdat, bl(icax), bl(icay),
     $                bl(icaz),  bl(icbx), bl(icby), bl(icbz), bl(iyra),
     $                bl(iyrb),    bl(iq),bl(ixyzi),  bl(iar),bl(ipspi),
     $              bl(ipspic), bl(ilder),   membuf,     dens,       na,
     $                     tol,      forc,     time,  maxpass)
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
      subroutine psp_forc(ics,  jcs,   npsp,  nradp,maxlpsp,
     $                    ncs,  ncf,  maxlp,   maxl,   mijk,
     $                iecpdat,nrrad,  aradp,  cradp,    len,
     $                   xnuc,  inx, basdat,    cax,    cay,
     $                    caz,  cbx,    cby,    cbz,    yra,
     $                    yrb,    q,   xyzi,     ar,   pspi,
     $                  pspic, lder, membuf,   dens,     na,
     $                    tol, forc,   time,maxpass)
      implicit real*8 (a-h,o-z)
c...
c...  compute psp contribution to atomic forces for the contracted
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
c...  forc       forces
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
      dimension dens(*),forc(3*na)
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
      mxpw=max(2,lmaxa+1,lmaxb+1)
c...
c...  derivative to atom mapping
c...
c...  at the end of this section, nder will contain the
c...  number of atomic derivatives that have to be computed
c...  for this shell pair
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
         write(iout,*)'Psp_Forc: not enough memory.'
         write(iout,*)'At least ',9*nints-membuf,' more words needed.'
         write(iout,*)
         call nerror(1,'Psp_Forc','Not enough memory',0,0)
      endif
c...
c...  start passes
c...
      ipass=0
      iniz=1
 1    ipass=ipass+1     ! starting point for  passes
      maxpass=max(maxpass,ipass)
c     call getival('iout',iout)
c     if(ipass.gt.1)write(iout,*)'Psp_Forc is doing multiple passes'
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
          do id=1,nder
            if(icatom.eq.lder(id))icat=id !pointer to atom in der. array
          enddo
          if(icat.eq.0)
     $      call nerror(2,'Psp_forc','pointer to icatom is 0 !!!',0,0)
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
     $                 mxpw,lima,limb,cax,cay,caz,cbx,cby,cbz)
          ca2=cax(2)+cay(2)+caz(2)
          cb2=cbx(2)+cby(2)+cbz(2)
c...
c...  accumulate integrals over primitives
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
        enddo   ! end loop over pseudopotentials
c       call secund(t1)
c       time(8)=time(8)+t1-t0
c       call secund(t0)
c...
c...  normalization and contraction
c...
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
c...  add contribution to atomic forces
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
                  ind=nder3*((icb-1)*nd3+(ica-1)*nd2+(jgc-1)*ngci+igc-1)
                  dij=dens(ii+icj)
                  if(ici.ne.icj)dij=dij+dij
                  do id=1,nder
                    idat=lder(id)
                    ida3=3*(idat-1)
                    idx=ida3+1
                    idy=ida3+2
                    idz=ida3+3
                    id3=3*(id-1)
                    icpx=id3+1+ind
                    icpy=id3+2+ind
                    icpz=id3+3+ind
                    forc(idx)=forc(idx)+dij*pspic(icpx)
                    forc(idy)=forc(idy)+dij*pspic(icpy)
                    forc(idz)=forc(idz)+dij*pspic(icpz)
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
c=======================================================================
      subroutine psp_grad12i(
     $                    basdat,   lla,   mma,   nna,    llb,
     $                       mmb,   nnb,   nga,   ngb,   lpsp,
     $                     nrpsp, alpha,    cc,  len0,    len,
     $                       cax,   cay,   caz,   cbx,    cby,
     $                       cbz,  lima,  limb,  mxpw,  nmaxa,
     $                     nmaxb,  mijk,   yra,   yrb,   xyzi,
     $                         q,    ar,  nder, iatom,  jatom,
     $                    icatom,   tol,argmax,  ilen,    ifi,
     $                       ca2,  jlen,   jfi,   cb2,    psp,
     $                      time)
      implicit real*8 (a-h,o-z)
c...
c...   MM September 2004
c...
c...   pseudopotential integral derivatives over primitives
c...
c...  The integrals are differentiated with respect the centers
c...  of Gaussians A and B. The derivative with respect
c...  the pseudopotential center are obtained by translational
c...  invariance
c...
c...   Input parameters
c...
c...   basdat           basis set data
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
c...
c...   Output prameter
c...
c...   psp              array of pseudopotential integral derivatives
c...
      dimension basdat(13,*)
      dimension time(*)
c...
      parameter (zero=0.0d0,two=2.0d0)
      logical compute
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
c...
c...  assemble integrals over primitives
c...
c           call secund(tini)
            if(compute)then
            call psp_lg1(1,nga,1,ngb,alpa,alpb,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,nder,iatom,jatom,icatom,
     $                   ip,jp,ilen,jlen,nga,ngb,psp)
            endif
c           call secund(tifi)
c           time(7)=time(7)+tifi-tini
          endif
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_lg1(igas,igaf,igbs,igbf,alpa,alpb,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,nder,iatom,jatom,icatom,
     $                   i1,i2,nd1,nd2,nga,ngb,psp)
      implicit real*8 (a-h,o-z)
c...
c...  loop over the Gaussian component and assembling
c...  of the first-order pseudopotential integrals
c...
      dimension lla(nga),mma(nga),nna(nga),llb(ngb),mmb(ngb),nnb(ngb)
      dimension lima(3),limb(3)
      dimension cax(0:mxpw),cay(0:mxpw),caz(0:mxpw)
      dimension cbx(0:mxpw),cby(0:mxpw),cbz(0:mxpw)
      dimension ar(0:nmaxa,0:nmaxa,0:nmaxa,0:nmaxb,0:nmaxb,0:nmaxb)
      dimension psp(3,nder,nd1,nd2,nga,ngb)
c...
c...  variables for summation loops
c...
      logical doa,dob,doc,dod,doe,dof
      integer a,al,af,aml,amf,apl,apf
      integer b,bl,bf,bml,bmf,bpl,bpf
      integer c,cl,cf,cml,cmf,cpl,cpf
      integer d,dl,df,dml,dmf,dpl,dpf
      integer e,el,ef,eml,emf,epl,epf
      integer f,fl,ff,fml,fmf,fpl,fpf
c...
      parameter (zero=0.0d0,two=2.0d0)
c...
c...   table of binomial coefficients: BC(m,n)=binomial(m,n)
c...
      real*8 BC(0:16,0:16)
      save BC
      data BC/1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0
     $  D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,1.0D0,0.0D0,1.0D0,2.0D0,3
     $  .0D0,4.0D0,5.0D0,6.0D0,7.0D0,8.0D0,9.0D0,1.0D1,1.1D1,1.2D1,1.3D1
     $  ,1.4D1,1.5D1,1.6D1,0.0D0,0.0D0,1.0D0,3.0D0,6.0D0,1.0D1,1.5D1,2.1
     $  D1,2.8D1,3.6D1,4.5D1,5.5D1,6.6D1,7.8D1,9.1D1,1.05D2,1.2D2,0.0D0,
     $  0.0D0,0.0D0,1.0D0,4.0D0,1.0D1,2.0D1,3.5D1,5.6D1,8.4D1,1.2D2,1.65
     $  D2,2.2D2,2.86D2,3.64D2,4.55D2,5.6D2,0.0D0,0.0D0,0.0D0,0.0D0,1.0D
     $  0,5.0D0,1.5D1,3.5D1,7.0D1,1.26D2,2.1D2,3.3D2,4.95D2,7.15D2,1.001
     $  D3,1.365D3,1.82D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,6.0D0,2.1D
     $  1,5.6D1,1.26D2,2.52D2,4.62D2,7.92D2,1.287D3,2.002D3,3.003D3,4.36
     $  8D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,7.0D0,2.8D1,8.4D1,
     $  2.1D2,4.62D2,9.24D2,1.716D3,3.003D3,5.005D3,8.008D3,0.0D0,0.0D0,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,8.0D0,3.6D1,1.2D2,3.3D2,7.92
     $  D2,1.716D3,3.432D3,6.435D3,1.144D4,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
     $  ,0.0D0,0.0D0,0.0D0,1.0D0,9.0D0,4.5D1,1.65D2,4.95D2,1.287D3,3.003
     $  D3,6.435D3,1.287D4,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0
     $  D0,0.0D0,1.0D0,1.0D1,5.5D1,2.2D2,7.15D2,2.002D3,5.005D3,1.144D4,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D
     $  0,1.1D1,6.6D1,2.86D2,1.001D3,3.003D3,8.008D3,0.0D0,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.2D1,7.8D1
     $  ,3.64D2,1.365D3,4.368D3,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D
     $  0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.3D1,9.1D1,4.55D2,1.82D3,
     $  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D
     $  0,0.0D0,0.0D0,1.0D0,1.4D1,1.05D2,5.6D2,0.0D0,0.0D0,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0
     $  ,1.5D1,1.2D2,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0
     $  D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,1.0D0,1.6D1,0.0D0,0.0D0,0
     $  .0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0,0.0D0
     $  ,0.0D0,0.0D0,0.0D0,1.0D0/
c...
c...  loops over the Gaussian components
c...
      talpam=-two*alpa
      talpbm=-two*alpb
      do iga=igas,igaf
        al=lla(iga)
        af=min(al,lima(1))
        fal=dfloat(al)
        doa=al.gt.0
        aml=al-1
        amf=min(aml,lima(1))
        apl=al+1
        apf=min(apl,lima(1))
c...
        bl=mma(iga)
        bf=min(bl,lima(2))
        fbl=dfloat(bl)
        dob=bl.gt.0
        bml=bl-1
        bmf=min(bml,lima(2))
        bpl=bl+1
        bpf=min(bpl,lima(2))
c...
        cl=nna(iga)
        cf=min(cl,lima(3))
        fcl=dfloat(cl)
        doc=cl.gt.0
        cml=cl-1
        cmf=min(cml,lima(3))
        cpl=cl+1
        cpf=min(cpl,lima(3))
        do igb=igbs,igbf
          dl=llb(igb)
          df=min(dl,limb(1))
          fdl=dfloat(dl)
          dod=dl.gt.0
          dml=dl-1
          dmf=min(dml,limb(1))
          dpl=dl+1
          dpf=min(dpl,limb(1))
c...
          el=mmb(igb)
          ef=min(el,limb(2))
          fel=dfloat(el)
          doe=el.gt.0
          eml=el-1
          emf=min(eml,limb(2))
          epl=el+1
          epf=min(epl,limb(2))
c...
          fl=nnb(igb)
          ff=min(fl,limb(3))
          ffl=dfloat(fl)
          dof=fl.gt.0
          fml=fl-1
          fmf=min(fml,limb(3))
          fpl=fl+1
          fpf=min(fpl,limb(3))
c...
c...  x derivatives
c...
          ppintxa=zero
          ppintxb=zero
          do b=bl,bf,-1
            bx=BC(bl,b)*cay(bl-b)
            do c=cl,cf,-1
              bcx=bx*BC(cl,c)*caz(cl-c)
              do e=el,ef,-1
                bcex=bcx*BC(el,e)*cby(el-e)
                do f=fl,ff,-1
                  bcefx=bcex*BC(fl,f)*cbz(fl-f)
c...
c...  Center A
c...
                  do d=dl,df,-1
                    bcdefx=bcefx*BC(dl,d)*cbx(dl-d)
                    if(doa)then
                      do a=aml,amf,-1           ! shifted down
                       ppintxa=ppintxa+bcdefx*fal*BC(aml,a)*cax(aml-a)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do a=apl,apf,-1             ! shifted up
                     ppintxa=ppintxa+bcdefx*talpam*BC(apl,a)*cax(apl-a)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
c...
c...  Center B
c...
                  do a=al,af,-1
                    abcefx=bcefx*BC(al,a)*cax(al-a)
                    if(dod)then
                      do d=dml,dmf,-1           ! shifted down
                       ppintxb=ppintxb+abcefx*fdl*BC(dml,d)*cbx(dml-d)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do d=dpl,dpf,-1             ! shifted up
                     ppintxb=ppintxb+abcefx*talpbm*BC(dpl,d)*cbx(dpl-d)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
c...
c...  y derivatives
c...
          ppintya=zero
          ppintyb=zero
          do a=al,af,-1
            ax=BC(al,a)*cax(al-a)
            do c=cl,cf,-1
              acx=ax*BC(cl,c)*caz(cl-c)
              do d=dl,df,-1
                acdx=acx*BC(dl,d)*cbx(dl-d)
                do f=fl,ff,-1
                  acdfx=acdx*BC(fl,f)*cbz(fl-f)
c...
c...  Center A
c...
                  do e=el,ef,-1
                    acdefx=acdfx*BC(el,e)*cby(el-e)
                    if(dob)then
                      do b=bml,bmf,-1           ! shifted down
                       ppintya=ppintya+acdefx*fbl*BC(bml,b)*cay(bml-b)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do b=bpl,bpf,-1             ! shifted up
                     ppintya=ppintya+acdefx*talpam*BC(bpl,b)*cay(bpl-b)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
c...
c...  Center B
c...
                  do b=bl,bf,-1
                    abcdfx=acdfx*BC(bl,b)*cay(bl-b)
                    if(doe)then
                      do e=eml,emf,-1           ! shifted down
                       ppintyb=ppintyb+abcdfx*fel*BC(eml,e)*cby(eml-e)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do e=epl,epf,-1             ! shifted up
                     ppintyb=ppintyb+abcdfx*talpbm*BC(epl,e)*cby(epl-e)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
c...
c...  z derivatives
c...
          ppintza=zero
          ppintzb=zero
          do a=al,af,-1
            ax=BC(al,a)*cax(al-a)
            do b=bl,bf,-1
              abx=ax*BC(bl,b)*cay(bl-b)
              do d=dl,df,-1
                abdx=abx*BC(dl,d)*cbx(dl-d)
                do e=el,ef,-1
                  abdex=abdx*BC(el,e)*cby(el-e)
c...
c...  Center A
c...
                  do f=fl,ff,-1
                    abdefx=abdex*BC(fl,f)*cbz(fl-f)
                    if(doc)then
                      do c=cml,cmf,-1           ! shifted down
                       ppintza=ppintza+abdefx*fcl*BC(cml,c)*caz(cml-c)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do c=cpl,cpf,-1             ! shifted up
                     ppintza=ppintza+abdefx*talpam*BC(cpl,c)*caz(cpl-c)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
c...
c...  Center B
c...
                  do c=cl,cf,-1
                    abcdex=abdex*BC(cl,c)*caz(cl-c)
                    if(dof)then
                      do f=fml,fmf,-1           ! shifted down
                       ppintzb=ppintzb+abcdex*ffl*BC(fml,f)*cbz(fml-f)*
     $                                        ar(a,b,c,d,e,f)
                      enddo
                    endif
                    do f=fpl,fpf,-1             ! shifted up
                     ppintzb=ppintzb+abcdex*talpbm*BC(fpl,f)*cbz(fpl-f)*
     $                                      ar(a,b,c,d,e,f)
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
c...
c...   assembling of integrals
c...
          if(iatom.ne.0)then
          psp(1,iatom,i1,i2,iga,igb)=psp(1,iatom,i1,i2,iga,igb)+ppintxa
          psp(2,iatom,i1,i2,iga,igb)=psp(2,iatom,i1,i2,iga,igb)+ppintya
          psp(3,iatom,i1,i2,iga,igb)=psp(3,iatom,i1,i2,iga,igb)+ppintza
          endif
          if(jatom.ne.0)then
          psp(1,jatom,i1,i2,iga,igb)=psp(1,jatom,i1,i2,iga,igb)+ppintxb
          psp(2,jatom,i1,i2,iga,igb)=psp(2,jatom,i1,i2,iga,igb)+ppintyb
          psp(3,jatom,i1,i2,iga,igb)=psp(3,jatom,i1,i2,iga,igb)+ppintzb
          endif
          if(icatom.ne.0)then
          psp(1,icatom,i1,i2,iga,igb)=psp(1,icatom,i1,i2,iga,igb)
     $                                -ppintxa-ppintxb
          psp(2,icatom,i1,i2,iga,igb)=psp(2,icatom,i1,i2,iga,igb)
     $                                -ppintya-ppintyb
          psp(3,icatom,i1,i2,iga,igb)=psp(3,icatom,i1,i2,iga,igb)
     $                                -ppintza-ppintzb
          endif
        enddo
      enddo
      end

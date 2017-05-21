c...
c...  MM September-October 2004
c...
c...  thys file contains the driver routines for the pseudopotential
c...  contribution to the First order Fock matrix for the calculation
c...  of NMR chemical shift over a GIAO basis set.
c...
c======================================================================
      subroutine pspFnmr(na,npsp,ncf,ncs,inx,xnuc,basdat,ntri,fd)

      use memory

      implicit real*8 (a-h,o-z)
c...
c...  driver routine for the psp contribution to the derivative
c...  Fock matrix for the nmr chemical shift
c...
c...  na      number of atoms
c...  npsp    number of pseudopotentials
c...  ncf     number of contracted functions
c...  ncs     number of contracted shells
c...  inx     contraction information
c...  xnuc    nuclear data
c...  basdat  basis set data
c...  ntri    dimension of the fd array
c...  fd      derivative Fock matrix
c...
      dimension inx(12,*),basdat(13,*),xnuc(5,*)
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
c...  increase maxl, because we are dealing with shifted
c...  Gaussians
c...
      maxl=maxl+1
c...
c...  memory setup
c...
      mlp=max(maxl,maxlp)
      maxpw=max(maxl,2)+1
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
      call getmem(3*(maxd2)**2,ipspi) !integrals over primitives
      call getmem(3*(maxd2)**2,ipspic)!integrals over contracted
      call getmem(maxlp+1,ilen)             ! array for psp contraction
      call getmem(1,memlim)
      mempsp=memlim-icax
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
      do ics=1,ncs
        do jcs=1,ics
c...
c...  compute integrals
c...
          call psp_Fnmr(ics,       jcs,      npsp,     nradp,   maxlpsp,
     $                  ncs,       ncf,     maxlp,      maxl,      mijk,
     $         bl(iiecpdat),bl(inrrad), bl(iarad), bl(icrad),  bl(ilen),
     $                 xnuc,       inx,    basdat,  bl(icax),  bl(icay),
     $             bl(icaz),  bl(icbx),  bl(icby),  bl(icbz),  bl(iyra),
     $             bl(iyrb),    bl(iq), bl(ixyzi),   bl(iar), bl(ipspi),
     $           bl(ipspic),       tol,      ntri,        fd,      time)
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
c...
      end
c======================================================================
      subroutine psp_Fnmr(ics,  jcs,   npsp,  nradp,maxlpsp,
     $                    ncs,  ncf,  maxlp,   maxl,   mijk,
     $                iecpdat,nrrad,  aradp,  cradp,    len,
     $                   xnuc,  inx, basdat,    cax,    cay,
     $                    caz,  cbx,    cby,    cbz,    yra,
     $                    yrb,    q,   xyzi,     ar,   pspi,
     $                  pspic,  tol,   ntri,     fd,   time)
      implicit real*8 (a-h,o-z)
c...
c...  compute psp contribution to the nmr fock derivative matrix for
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
c...  tol        tolerance threshold
c...  ntri       dimension of fd array
c...  fd         Fock derivative matrix
c...
      dimension iecpdat(maxlpsp+6,npsp),nrrad(nradp),aradp(nradp),
     $          cradp(nradp),len(0:maxlp)
      dimension xnuc(5,*),inx(12,*),basdat(13,*)
      dimension cax(0:maxl),cay(0:maxl),caz(0:maxl)
      dimension cbx(0:maxl),cby(0:maxl),cbz(0:maxl)
      dimension pspi(*),pspic(*)
      dimension fd(ntri,3)
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
c...  we have contribution only if iatom.ne.jatom
c...
      if(iatom.eq.jatom)return
c...
c...  zero out integrals over primitives
c...
c     call secund(t0)
      call zeroit(pspi,3*nga*ngb*ilen*jlen)
c     call secund(tz)
c     time(9)=time(9)+tz-t0
c...
c...  loop over pseudopotentials
c...
      do ipsp=1,npsp
        icatom=iecpdat(3,ipsp)   ! atomic center of psp
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
     $               mxpw,lima,limb,cax,cay,caz,cbx,cby,cbz)
        ca2=cax(2)+cay(2)+caz(2)
        cb2=cbx(2)+cby(2)+cbz(2)
c...
c...  accumulate GIAO derivative integrals over primitives
c...
        call psp_giao12i(
     $            basdat, mulx(nsa), muly(nsa), mulz(nsa), mulx(nsb),
     $         muly(nsb), mulz(nsb),       nga,       ngb,      lpsp,
     $        nrrad(irl),aradp(irl),cradp(irl),      len0,       len,
     $               cax,       cay,       caz,       cbx,       cby,
     $               cbz,      lima,      limb,      mxpw,   lmaxa+1,
     $           lmaxb+1,      mijk,       yra,       yrb,      xyzi,
     $                 q,        ar,         a,         b,      tolp,
     $            argmax,      ilen,       ifi,       ca2,      jlen,
     $               jfi,       cb2,      pspi,      time)
      enddo   ! end loop over pseudopotentials
c     call secund(t1)
c     time(8)=time(8)+t1-t0
c     call secund(t0)
c...
c...  contraction and normalization.
c...
      call zeroit(pspic,3*nga*ngb*ngci*ngcj)
      call psp_contrn(ilen,ifi,ngci,nga,itypa,
     $                jlen,jfi,ngcj,ngb,itypb,
     $                basdat,3,pspi,pspic)
c     call secund(t1)
c     time(3)=time(3)+t1-t0
c...
c...  eventually transform to spherical components.
c...  pspi is used as scratch space
c...
c     call secund(t0)
      if(transa.or.transb)then
        call cart2sphera(pspic,3*ngci*ngcj,itypa,transa,nga,nsha,
     $                           itypb,transb,ngb,nshb,pspi)
      endif
c     call secund(t1)
c     time(4)=time(4)+t1-t0
c...
c...  add contribution to Fock derivatives
c...
c     call secund(t0)
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
              if(ici.gt.icj)then
                ind=3*((icb-1)*nd3+(ica-1)*nd2+(jgc-1)*ngci+igc-1)
                iijj=ii+icj
                fd(iijj,1)=fd(iijj,1)-pspic(1+ind)
                fd(iijj,2)=fd(iijj,2)-pspic(2+ind)
                fd(iijj,3)=fd(iijj,3)-pspic(3+ind)
              endif
            enddo
          enddo
        enddo
      enddo
c     call secund(t1)
c     time(5)=time(5)+t1-t0
c...
      end
c=======================================================================
      subroutine psp_giao12i(
     $                    basdat,   lla,   mma,   nna,   llb,
     $                       mmb,   nnb,   nga,   ngb,  lpsp,
     $                     nrpsp, alpha,    cc,  len0,   len,
     $                       cax,   cay,   caz,   cbx,   cby,
     $                       cbz,  lima,  limb,  mxpw, nmaxa,
     $                     nmaxb,  mijk,   yra,   yrb,  xyzi,
     $                         q,    ar,    ra,    rb,   tol,
     $                    argmax,  ilen,   ifi,   ca2,  jlen,
     $                       jfi,   cb2,   psp,  time)
      implicit real*8 (a-h,o-z)
c...
c...   MM September 2004
c...
c...   pseudopotential integrals over first order GIAO derivatives.
c...
c...   A Gauge Including Atomic Orbital (GIAO) is defined as
c...
c...   chi_a(H)= exp(-%i/(2*c) * (H x R_a) * r) * chi_a(0),
c...
c...   where H is the magnetic field vector, R_a is the position
c...   vector of center a (H x R_a means vector product of H and R_a),
c...   r is the position vector, %i is the imaginary unit, c is
c...   the speed of light, and chi_a(0) is the unperturbed AO
c...   (i.e., a cartesian Gaussian function)
c...
c...   The first derivative of the GIAO with respect to the
c...   magnetic field, computed at the point H=0 is:
c...
c...   chi1_a(0)= -%i/(2*c) * (R_a x r) * chi_a(0)
c...
c...   developing the vector product one obtains (for the x component):
c...
c...   chi1_ax = -%i/(2*c) * (R_ay * z - R_az * y) * chi_a(0),
c...
c...   with analogous expressions for the components y and z.
c...
c...   Now we substitute the explicit form of the cartesian Gaussian:
c...
c...   chi1_ax = -%i/(2*c) * (R_ay * z - R_az * y) * ax^la *
c...             ay^ma * az^na * exp(-alpa * |r_a|^2),
c...
c...   where alpa is the exponent of the Gaussian, and la,ma,na are the
c...   three quantum numbers, r_a is the vector describing the
c...   electron position with respect to center A, and xa, ya and za
c...   are the components of the vector r_a
c...
c...   Given two shells of cartesian Gaussians fuctions A and B,
c...   and a pseudopotential V_psp, this subroutine compute integrals
c...   of the type:
c...
c...   psp_x = <chi1_ax | V_psp | chi0_b> + <chi0_a | V_psp | chi1_bx>
c...
c...   with similar expressions for the components y and z.
c...
c...   To compute the integral, we need to express the factor
c...   (R_ay * z - R_az * y) with respect the centers A or B.
c...   This can be accomplished by the susbtitutions:
c...
c...             x=ax+R_ax    y=ay+R_ay     z=az+R_az
c...    or
c...             x=bx+R_bx    y=by+R_by     z=bz+R_bz
c...
c...   susbtituting with respect to center A, we obtain
c...
c...   psp_xa = (R_by - R_ay) < az*chi_a(0) | V_psp | chi_b(0) > -
c...            (R_bz - R_az) < ay*chi_a(0) | V_psp | chi_b(0) > +
c...   (R_by * R_az - R_bz * R_ay) < chi_a(0) | V_psp | chi_b(0) >
c...
c...  where, for instance:
c...
c...  az*chi_a(0)= ax^la * ay^ma *az^(na+1) * exp(-alpha*|r_a|^2)
c...
c...  is a Gaussian with the quantum number na increased by 1.
c...
c...   In the same way, substituting with respect to center B:
c...
c...   psp_xb = (R_by - R_ay) < chi_a(0) | V_psp | bz*chi_b(0) > -
c...            (R_bz - R_az) < chi_a(0) | V_psp | by*chi_b(0) > +
c...   (R_by * R_az - R_bz * R_ay) < chi_a(0) | V_psp | chi_b(0) >
c...
c...  Due to the form of the pseudopotential term V_psp, the
c...  two expressions psp_xa and psp_xb do not generally give
c...  the same result, thus the integral is computed as the average
c...  of the two above expressions:
c...
c...          psp_x = 1/2 * ( psp_xa + psp_xb )
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
c...   ra,rb            cartesian coordinates of centers A and B
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
c     call secund(t1)
c     time(1)=time(1)+t1-t0
c...
c...  assemble integrals over primitives
c...
            if(compute)then
            call psp_lgiao(1,nga,1,ngb,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,ra,rb,
     $                   ip,jp,ilen,jlen,nga,ngb,psp)
            endif
c           call secund(tifi)
c           time(7)=time(7)+tifi-t1
          endif
        enddo
      enddo
      end
c=======================================================================
      subroutine psp_lgiao(igas,igaf,igbs,igbf,
     $                   lla,mma,nna,llb,mmb,nnb,lima,limb,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxb,ar,ra,rb,
     $                   i1,i2,nd1,nd2,nga,ngb,psp)
      implicit real*8 (a-h,o-z)
c...
c...  loop over the Gaussian component and assembling
c...  of the first-order pseudopotential integrals over GIAO
c...
      dimension lla(nga),mma(nga),nna(nga),llb(ngb),mmb(ngb),nnb(ngb)
      dimension ra(3),rb(3)
      dimension lima(3),limb(3)
      dimension psp(3,nd1,nd2,nga,ngb)
c...
c...  variables for summation loops
c...
      integer al,af,ap1l,ap1f
      integer bl,bf,bp1l,bp1f
      integer cl,cf,cp1l,cp1f
      integer dl,df,dp1l,dp1f
      integer el,ef,ep1l,ep1f
      integer fl,ff,fp1l,fp1f
c...
      parameter (zero=0.0d0,two=2.0d0)
c...
c...  loop over the Gaussian components
c...
      r01=(rb(2)*ra(3)-rb(3)*ra(2))
      r02=(rb(3)*ra(1)-rb(1)*ra(3))
      r03=(rb(1)*ra(2)-rb(2)*ra(1))
      rab1=(rb(1)-ra(1))
      rab2=(rb(2)-ra(2))
      rab3=(rb(3)-ra(3))
      do iga=igas,igaf
        al=lla(iga)
        af=min(al,lima(1))
        ap1l=al+1
        ap1f=min(ap1l,lima(1))
c...
        bl=mma(iga)
        bf=min(bl,lima(2))
        bp1l=bl+1
        bp1f=min(bp1l,lima(2))
c...
        cl=nna(iga)
        cf=min(cl,lima(3))
        fcl=dfloat(cl)
        cp1l=cl+1
        cp1f=min(cp1l,lima(3))
        do igb=igbs,igbf
          dl=llb(igb)
          df=min(dl,limb(1))
          dp1l=dl+1
          dp1f=min(dp1l,limb(1))
c...
          el=mmb(igb)
          ef=min(el,limb(2))
          ep1l=el+1
          ep1f=min(ep1l,limb(2))
c...
          fl=nnb(igb)
          ff=min(fl,limb(3))
          fp1l=fl+1
          fp1f=min(fp1l,limb(3))
c...
c...  summations
c...
c...  unshifted
c...
      call secund(t0)
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppint)
c...
c...  shift on center A
c...
          call psp_suml2(ap1l,ap1f,bl,bf,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppaxp1)
          call psp_suml2(al,af,bp1l,bp1f,cl,cf,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppayp1)
          call psp_suml2(al,af,bl,bf,cp1l,cp1f,
     $                   dl,df,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppazp1)
c...
c...  shift on center B
c...
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dp1l,dp1f,el,ef,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbxp1)
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,ep1l,ep1f,fl,ff,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbyp1)
          call psp_suml2(al,af,bl,bf,cl,cf,
     $                   dl,df,el,ef,fp1l,fp1f,
     $                   cax,cay,caz,cbx,cby,cbz,mxpw,
     $                   nmaxa,nmaxa,nmaxa,nmaxb,nmaxb,nmaxb,ar,ppbzp1)
c...
c...  assemble the integrals
c...
c...  note that we do not multiply by 1/(2*c) at this stage
c...
          x0=r01*ppint
          y0=r02*ppint
          z0=r03*ppint
c...
          xa=rab2*ppazp1-rab3*ppayp1
          ya=rab3*ppaxp1-rab1*ppazp1
          za=rab1*ppayp1-rab2*ppaxp1
c...
          xb=rab2*ppbzp1-rab3*ppbyp1
          yb=rab3*ppbxp1-rab1*ppbzp1
          zb=rab1*ppbyp1-rab2*ppbxp1
c...
c...  fill the psp array, taking the average of the above integrals
c...
          psp(1,i1,i2,iga,igb)=psp(1,i1,i2,iga,igb)+0.5d0*(xa+xb)+x0
          psp(2,i1,i2,iga,igb)=psp(2,i1,i2,iga,igb)+0.5d0*(ya+yb)+y0
          psp(3,i1,i2,iga,igb)=psp(3,i1,i2,iga,igb)+0.5d0*(za+zb)+z0
c...
c         ax=-ra(2)*ppazp1+ra(3)*ppayp1
c         ay=-ra(3)*ppaxp1+ra(1)*ppazp1
c         az=-ra(1)*ppayp1+ra(2)*ppaxp1
c...
c         bx= rb(2)*ppbzp1-rb(3)*ppbyp1
c         by= rb(3)*ppbxp1-rb(1)*ppbzp1
c         bz= rb(1)*ppbyp1-rb(2)*ppbxp1
c...
c         psp(1,iga,igb)=psp(1,iga,igb)+ax+bx
c         psp(2,iga,igb)=psp(2,iga,igb)+ay+by
c         psp(3,iga,igb)=psp(3,iga,igb)+az+bz
        enddo
      enddo
      end

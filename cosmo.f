c
c   this file contains the COSMO code from  Andreas Klamt and
c   Michael Diedenhofen
c
c   MM October-December 2003:
c     - commented out the Fortran 90 statements
c       for dynamic memory allocations. The dynamically
c       allocated arrays are now parameters of the subroutines
c       where they are used, and scratch memory must be allocated
c       for them before the call
c
c     - changed Fortran 90 definition of character variables
c
c     - changed some double precision to real*8
c
c     - changed calls to dzero and izero to zeroit and izeroit
c
c     - added some ierr=-1 before returns any errors
c
c======================================================================
      subroutine ansude(ra,rb,d,rs,ara,arb,aar,abr,arad,arbd)
c   This progrom calculates the area of two intersecting spheres
c   with radii ra and rb at a distance d and a solvent probe
c   radius rs. The two areas are calculated separately (ara,arb).
c   For both areas analytic derivatives with respect to the distance
c   d are calculated (arad,arbd).    (written by Andreas Klamt, 9/9/96)
      implicit real*8 (a-h,o-z)
      pi=2.d0*asin(1.d0)
      qa=ra+rs
      qb=rb+rs
      ca=(qa**2+d**2-qb**2)/(2.d0*qa*d)
      cb=(qb**2+d**2-qa**2)/(2.d0*qb*d)
      sa=sqrt(1.d0-ca*ca)
      sb=sqrt(1.d0-cb*cb)
      ta=pi*sa
      tb=pi*sb
      fza=(1.d0-cos(ta))/2
      fzb=(1.d0-cos(tb))/2
      if (sa.lt.0 .or. sb.lt.0) fza=1.d0
      if (sa.lt.0 .or. sb.lt.0) fzb=1.d0
      xa=fzb**1*rs*(ca+cb)
      xb=fza**1*rs*(ca+cb)
      ya=ra*sa-fzb*rb*sb
      yb=rb*sb-fza*ra*sa
      za=sqrt(xa*xa+ya*ya)
      zb=sqrt(xb*xb+yb*yb)
      ara=pi*ra*(2.d0*(1.d0+ca)*ra+sa*za)
      arb=pi*rb*(2.d0*(1.d0+cb)*rb+sb*zb)
      aar=pi*ra*(sa*za)
      abr=pi*rb*(sb*zb)
c now derivatives
      cad=(qb**2+d**2-qa**2)/(2.d0*qa*d*d)
      cbd=(qa**2+d**2-qb**2)/(2.d0*qb*d*d)
      sad=-ca*cad/sa
      sbd=-cb*cbd/sb
      tad=pi*sad
      tbd=pi*sbd
      fzad=sin(ta)*.5d0
      fzbd=sin(tb)*.5d0
      if (sa.lt.0 .or. sb.lt.0) fzad=0.d0
      if (sa.lt.0 .or. sb.lt.0) fzbd=0.d0
      xad=rs*((ca+cb)*fzbd*tbd+fzb*(cad+cbd))
      xbd=rs*((ca+cb)*fzad*tad+fza*(cad+cbd))
      yad=ra*sad-fzbd*tbd*rb*sb-fzb*rb*sbd
      ybd=rb*sbd-fzad*tad*ra*sa-fza*ra*sad
      zad=(xa*xad+ya*yad)/za
      zbd=(xb*xbd+yb*ybd)/zb
      arad=pi*ra*(sad*za+sa*zad+2.d0*ra*cad)
      arbd=pi*rb*(sbd*zb+sb*zbd+2.d0*rb*cbd)
      end
c =====================================================================
       subroutine consts(coord, nuc, srad, cosurf, iatsp, dirsm,
     $                   dirsmh, dirvec, nar, ar, dirtm, nsetf, nset,
     $                   natoms, maxnps, nspa, nsph, nppa, disex,
     $                   disex2, rsolv, routf, nps, npspher,
     $                   area, volume, phsran, ampran, lcavity, tm,
     $                   ierr, mes, din, xsp, sude, nipa,
     $                   lipa, isude, nn)
c
c ---------------------------------------------------------------------
c     this routine constructs or updates the solvent-accessible
c     surface (sas)
c
c     calls:
c           ansude (no error message, ok)
c           writsude (ok)
c           get_file (ok)
c ---------------------------------------------------------------------
c
      implicit real*8 (a-h,o-z)

      parameter (nsring=20)

      character*100 mes
      dimension xd(3),rvx(3),rvy(3),xja(3),xjb(3),xjc(3),xta(3,3)
      dimension fz(2),yx(2),trp(3),
     &          svx(3),svy(3),tvx(3),tvy(3),cc(3),ss(3),em(3),ee(3,3)
      dimension xx(3),xa(3),xi(3),xj(3),yy(3)
c
c -- this keeps overflowing as systems get bigger   ! MM Aug 2008
c -- let's see if 10 000 can stop the bleeding for a while
c
      integer maxntrp
      parameter (maxntrp=10000)  ! maximum number of triple points
      dimension phiset(2*(maxntrp+1)),iset(2*(maxntrp+1))
c
      dimension    coord(3,natoms), nuc(natoms), srad(natoms)
      dimension    cosurf(3,2*maxnps), iatsp(2*maxnps)
      dimension    dirsm(3,nspa), dirsmh(3,nsph), dirvec(3,nppa)
      dimension    nar(2*maxnps), ar(2*maxnps), dirtm(3,nppa)
      dimension    nsetf(2*maxnps), nset(natoms*nppa)
      dimension    tm(3,3,natoms)

c
c     dynamically allocated arrays. these are now
c     parameters of the subroutine. a23mat actually is not
c     used in this subroutine
c
c     dimensions for max. coordination number 12
c     integer, allocatable ::  nn(:,:), isude(:,:), nipa(:), lipa(:)
c     double precision, allocatable :: sude(:,:)
c
c     double precision, allocatable :: a23mat(:),xsp(:,:)
c     logical*1, allocatable ::    din(:)
c
      dimension xsp(3,2*maxnps),sude(4,6*natoms),nipa(4*natoms),
     $          lipa(4*12*natoms),isude(4,6*natoms),nn(3,natoms)
      logical*1 din(nppa+1)

      data pi/3.14159265358979323844d0/, four/4.d0/
      parameter (TollZero=1.0d-8)
c-----------------------------------------------------------------------
c     dynamic allocation
c     allocate(din(nppa+1), xsp(3,2*maxnps),
c    &         sude(4,6*natoms), nipa(4*natoms), lipa(4*12*natoms),
c    &         isude(4,6*natoms), nn(3,natoms), stat=ierr)
c     if (ierr .ne. 0) then
c        mes='COSMO: failed to allocate memory in subroutine consts.f'
c        return
c     endif
c-----------------------------------------------------------------------
      thrsh=1.0e-15
      fpinppa = four*pi/dble(nppa)
      fdiagr=2.1D0*sqrt(pi)
      sqpi=sqrt(pi)

c    more stable version from A. Klamt
      do 5 i=1,natoms
        do 5 j=1,3
    5   coord(j,i)=coord(j,i)+cos((phsran+5.0d0*i+j)**1.1d0)*ampran
c
      sdis=0.d0
      fdiag=1.05d0*sqrt(dble(nppa))
      inset=1
      ilipa=0
      nps = 0
c     --- initialise total surface area and volume
      area=0.d0
      volume=0.d0
      rmean=0.d0
c     --- loop over all atoms
      do 340 i=1,natoms
          npsoff=nps
          nipi=0
          ds=sqrt(4.d0/nspa)
          if (nuc(i) .eq. 1) ds=2.0d0*ds
          c2ds=cos(2.d0*ds)
          r=srad(i)
          ri=r-rsolv
          fvol=4*pi*ri**2/nppa
          do 20 ix=1,3
   20     xa(ix)=coord(ix,i)
          nps0=nps+1

c search for 3  nearest neighbor atoms
          dist1=1.d20
          dist2=1.d20
          dist3=1.d20
          nn1=0
          nn2=0
          nn3=0
cjwa
          do 70 j=1,natoms
              if (j.eq. i) go to 70
              dist=0.d0
              do 60 ix=1,3
   60         dist=dist+(xa(ix)-coord(ix,j))**2
              if(dist.lt.(r+srad(j))**2) then
                  ilipa=ilipa+1
                  nipi=nipi+1
c                 lipa(ilipa) points on neighbour atom ilipa
                  lipa(ilipa)=j
              end if
              if ((dist+0.05d0-dist3).lt.thrsh) then
                  dist3=dist
                  nn3=j
              end if
              if ((dist3+0.05d0-dist2).lt.thrsh) then
                  dist=dist2
                  dist2=dist3
                  dist3=dist
                  nn3=nn2
                  nn2=j
              end if
              if ((dist2+0.05d0-dist1).lt.thrsh) then
                  dist=dist1
                  dist1=dist2
                  dist2=dist
                  nn2=nn1
                  nn1=j
              end if
   70     continue
c         nipa(i) counts neighbour atoms of atom i
          nipa(i)=nipi
          nn(1,i)=nn1
          nn(2,i)=nn2
          nn(3,i)=nn3
c ---------------------------------------------------------------------
c        build new transformation matrix
c
c        nn1,nn2,nn3 are either atomic indices or 0.
c
c        if they are indices, the direction vectors on tm are built
c        by having the first vector point towards the atom nn1,
c        the second vector lying in the plane of atoms nn1,nn2 and
c        the central atom, and the third vector being the cross
c        product (so nn3 is irrelevant, actually).
c
c        if they are 0, defaults are assumed. note that the default
c        for the second direction is potentially suicidal since if
c        the first direction is 001 (ie nn1.ne.0) but nn2.eq.0, the
c        second direction is ill-defined. enjoy
c ---------------------------------------------------------------------

c         --- the first guy
          if (nn1 .eq. 0) then
              tm(1,1,i)=1.d0
              tm(1,2,i)=0.d0
              tm(1,3,i)=0.d0
          else
              dist1=0.d0
              do 80 ix=1,3
   80         dist1=dist1+(xa(ix)-coord(ix,nn1))**2
              dist=1.0d0/sqrt(dist1)
              tm(1,1,i)=(coord(1,nn1)-xa(1))*dist
              tm(1,2,i)=(coord(2,nn1)-xa(2))*dist
              tm(1,3,i)=(coord(3,nn1)-xa(3))*dist
          end if
c         --- the second guy
   90     if (nn2 .eq. 0) then
              dist=sqrt(9.d0*TM(1,2,I)**2+4.d0*tm(1,1,i)**2+
     &                  tm(1,3,i)**2)
              xx(1)=3.d0*TM(1,2,I)/dist
              xx(2)=TM(1,3,I)/dist
              xx(3)=2.d0*TM(1,1,I)/dist
          else
              dist2=0.d0
              do 100 ix=1,3
  100         dist2=dist2+(xa(ix)-coord(ix,nn2))**2
              dist=1.0d0/sqrt(dist2)
              xx(1)=(coord(1,nn2)-xa(1))*dist
              xx(2)=(coord(2,nn2)-xa(2))*dist
              xx(3)=(coord(3,nn2)-xa(3))*dist
          end if
          sp=xx(1)*tm(1,1,i)+xx(2)*tm(1,2,i)+xx(3)*tm(1,3,i)
          if ((sp*sp- 0.99d0).gt.thrsh) then
              nn2=nn3
              nn3=0
              dist2=dist3
              go to 90
          end if
          sininv=1.d0/sqrt(1.d0-sp*sp)
          tm(2,1,i)=(xx(1)-sp*tm(1,1,i))*sininv
          tm(2,2,i)=(xx(2)-sp*tm(1,2,i))*sininv
          tm(2,3,i)=(xx(3)-sp*tm(1,3,i))*sininv
c         --- the third guy is the cross product of the first two.
          tm(3,1,i)=tm(1,2,i)*tm(2,3,i)-tm(2,2,i)*tm(1,3,i)
          tm(3,2,i)=tm(1,3,i)*tm(2,1,i)-tm(2,3,i)*tm(1,1,i)
          tm(3,3,i)=tm(1,1,i)*tm(2,2,i)-tm(2,1,i)*tm(1,2,i)

c ---------------------------------------------------------------------
c         the atomic transformation matrix is now applied to the
c         vectors of the basic grid.
c ---------------------------------------------------------------------
          do 110 j=1,nppa
              xx(1)=dirvec(1,j)
              xx(2)=dirvec(2,j)
              xx(3)=dirvec(3,j)
              do 110 ix=1,3
                  x=xx(1)*tm(1,ix,i)+xx(2)*tm(2,ix,i)+xx(3)*tm(3,ix,i)
                  dirtm(ix,j)=x
  110     continue

c ---------------------------------------------------------------------
c         find the points of the basic grid on the sas
c         xa = coordinates of current atom
c         r  = radius of sphere around atom
c ---------------------------------------------------------------------
          narea=0
          do 160 j = 1,nppa
              din(j)=.false.
              do 130 ix=1,3
                  xx(ix) = xa(ix) + dirtm(ix,j)* r
  130         continue
c             --- only use points which do not lie inside another atom
c             --- mark those by setting din = .true.
c             --- we need only try those atoms intersecting atom i
              do 150 ik=ilipa-nipa(i)+1,ilipa
                  k = lipa(ik)
                  dist=0.0d0
                  do 140 ix=1,3
                      dist = dist + (xx(ix) - coord(ix,k))**2
  140             continue
                  dist=sqrt(dist)-srad(k)
                  if (dist.lt.0.0d0) go to 160
  150         continue
              narea=narea+1
              volume=volume+fvol*
     &          (dirtm(1,j)*xa(1)+dirtm(2,j)*xa(2)+dirtm(3,j)*xa(3)+ri)
              din(j)=.true.
  160     continue

c         --- narea is the # of basic grid points that survived the
c         --- cruel elimination process.
          if(narea.eq.0) goto 340

c         --- update the area. ri != r (pumpup postprocess)
          area=area+narea*ri*ri

c         --- switch between hydrogen and others.
          if (nuc(i) .eq. 1) then
              do j=1,nsph
                  nps=nps+1
                  if (nps .gt. maxnps) then

                     write(mes,'(2a)') 'COSMO: consts: nps is greater ',
     &               'than maxnps - use bigger maxnps if possible'
                     ierr = -1
                     return
                  endif
                  iatsp(nps)=i
                  xx(1)=dirsmh(1,j)
                  xx(2)=dirsmh(2,j)
                  xx(3)=dirsmh(3,j)
                  cosurf(1,nps)=
     1              xx(1)*tm(1,1,i)+xx(2)*tm(2,1,i)+xx(3)*tm(3,1,i)
                  cosurf(2,nps)=
     1              xx(1)*tm(1,2,i)+xx(2)*tm(2,2,i)+xx(3)*tm(3,2,i)
                  cosurf(3,nps)=
     1              xx(1)*tm(1,3,i)+xx(2)*tm(2,3,i)+xx(3)*tm(3,3,i)
              enddo
          else
              do j=1,nspa
                  nps=nps+1
                  if (nps .gt. maxnps) then
                     write(mes,'(2a)') 'COSMO: consts: nps is greater ',
     &               'than maxnps - use bigger maxnps if possible'
                     ierr = -1
                     return
                  endif
                  iatsp(nps)=i
                  xx(1)=dirsm(1,j)
                  xx(2)=dirsm(2,j)
                  xx(3)=dirsm(3,j)
                  cosurf(1,nps)=
     1              xx(1)*tm(1,1,i)+xx(2)*tm(2,1,i)+xx(3)*tm(3,1,i)
                  cosurf(2,nps)=
     1              xx(1)*tm(1,2,i)+xx(2)*tm(2,2,i)+xx(3)*tm(3,2,i)
                  cosurf(3,nps)=
     1              xx(1)*tm(1,3,i)+xx(2)*tm(2,3,i)+xx(3)*tm(3,3,i)
              enddo
          endif

  200     sdis0=sdis

          do 210 ips=nps0,nps
              nar(ips)=0
              xsp(1,ips)=0.d0
              xsp(2,ips)=0.d0
              xsp(3,ips)=0.d0
  210     continue

c ---------------------------------------------------------------------
c         every segment comprises several basic grid points. the
c         location of the segment is defined as the center of these
c         basic grid points (projected back onto the surface).
c
c         the # of basic grid points for a given segment may change
c         with geometry, resulting in 'jumps' in the segment location.
c         this is a possible problem when comparing analytical against
c         finite difference gradients (cmk).
c
c         the basic grid points enable a more accurate calculation
c         of segment-segment interactions, thereby approximating the
c         surface integrals (see eq.(7) in JCSPT2,799(1993))
c ---------------------------------------------------------------------
          do 250 j=1,nppa
              if (.not. din(j)) go to 250
              spm=-1.0001d0
              x1=dirtm(1,j)
              x2=dirtm(2,j)
              x3=dirtm(3,j)
              do 220 ips=nps0,nps
                  sp=x1*cosurf(1,ips)+x2*cosurf(2,ips)+x3*cosurf(3,ips)
                  if ((sp-spm) .lt.thrsh) go to 220
                  spm=sp
                  ipm=ips
  220         continue
              if ((spm-c2ds).lt.thrsh) then
                  nps=nps+1
                  do 230 ix=1,3
  230             cosurf(ix,nps)=dirtm(ix,j)
                  iatsp(nps)=i
                  go to 200
              end if
              nar(ipm)=nar(ipm)+1
              do 240 ix=1,3
  240         xsp(ix,ipm)=xsp(ix,ipm)+dirtm(ix,j)
  250     continue

          sdis=0.d0
          ips=nps0-1

  260     ips=ips+1
c ---------------------------------------------------------------------
c         if a segment has no basic grid points contributing to it,
c         compress data
c ---------------------------------------------------------------------

  352     if(nar(ips).eq.0)then
              nps=nps-1
              if(nps.lt.ips) goto 200
              do 369 jps=ips,nps
                  nar(jps)=nar(jps+1)
                  xsp(1,jps)=xsp(1,jps+1)
                  xsp(2,jps)=xsp(2,jps+1)
  369             xsp(3,jps)=xsp(3,jps+1)
              goto 352
          endif

          dist=0.d0
          do 280 ix=1,3
              x=xsp(ix,ips)
              dist=dist+x*x
  280     continue
          sdis=sdis+dist
          dist=1.d0/sqrt(dist)

c         --- final result : segment = average of basic grid points.

          do 290 ix=1,3
  290         cosurf(ix,ips)=xsp(ix,ips)*dist

          if(ips.lt.nps) goto 260

          if (abs(sdis-sdis0) .gt. 1.d-5) go to 200

c         --- final assignments of :
c            absolute segment coordinates (xsp)
c            # of basic grid points per segment (nar)
c            offsets of segments wrt nar (nsetf)
c            pointer from all basic grids to generic basic grid (nset)
c            [NOTE : dirvec will be used later on, not dirtm]


          do 310 ips=nps0,nps
              nsetf(ips)=inset
              inset=inset+nar(ips)
              nar(ips)=0
              do 300 ix=1,3
  300         xsp(ix,ips)=xa(ix)+cosurf(ix,ips)*ri
  310     continue

          do 330 j=1,nppa
              if (.not. din(j)) go to 330
              spm=-1.0001d0
              x1=dirtm(1,j)
              x2=dirtm(2,j)
              x3=dirtm(3,j)
              do 320 ips=nps0,nps
                  sp=x1*cosurf(1,ips)+x2*cosurf(2,ips)+x3*cosurf(3,ips)
                  if ((sp-spm) .lt.thrsh) go to 320
                  spm=sp
                  ipm=ips
  320         continue
              if ((spm-c2ds).lt.thrsh) go to 330
              nara=nar(ipm)
              nset(nsetf(ipm)+nara)=j
              nar(ipm)=nara+1
              ar(ipm) = fpinppa*ri*ri*nar(ipm)
  330     continue
          npsat=nps-npsoff
          rmean=rmean+npsat*ri
  340 continue
c     --- end loop over all atoms

c     --- some 4*PI multiplication to get the actual area
      area=area*fpinppa

c     --- calculate disex2
      rmean=rmean/nps
      disex2=4*rmean*rmean*disex*disex
c
c surface closure
      npspher = nps
      if (lcavity.eq.0) goto 3333
      do 2900 i=1,natoms
2900  din(i)=.true.
      do 2901 j=1,nps
2901  din(iatsp(j))=.false.
      nip=0

C generation of segments along the intersection rings
      ilipa=0
      do 3800 ia=1,natoms-1
        if(din(ia)) go to 3800
        ra=srad(ia)
        do 3001 ix=1,3
3001    xta(ix,1)=coord(ix,ia)
        do 3700 iib=ilipa+1,ilipa+nipa(ia)
          ib=lipa(iib)
          if(ib .le.ia) go to 3700
          if(din(ib)) go to 3700
          rb=srad(ib)
          dab=0.d0
          nsab=0
          do 3002 ix=1,3
            xta(ix,2)=coord(ix,ib)
            xx(ix)=xta(ix,2)-xta(ix,1)
3002      dab=dab+xx(ix)**2
          dab=sqrt(dab)
c -- tests added    ! JB Oct 2011
          cosa=(ra**2+dab**2-rb**2)/(2*dab*ra)
       If(Abs(cosa).GE.1.0d0) cosa=SIGN(1.0d0,cosa)
          cosb=(rb**2+dab**2-ra**2)/(2*dab*rb)
       If(Abs(cosb).GE.1.0d0) cosb=SIGN(1.0d0,cosb)
          sina=sqrt(1.d0-cosa**2)
          da=ra*cosa
          hh=ra*sina
          ddd=rsolv*(cosa+cosb)/dab
          fz(1)=(1.d0-cos(hh*pi/ra))/2
          fz(2)=(1.d0-cos(hh*pi/rb))/2
          if(cosa*cosb .lt. 0d0) fz(1)=1.d0
          if(cosa*cosb .lt. 0d0) fz(2)=1.d0
          cosgz=(ra*ra+rb*rb-dab*dab)/(ra*rb)
          dc=rsolv*sqrt(2-cosgz)
          yx(1)=rsolv/srad(ia)
          yx(2)=rsolv/srad(ib)
          do 3005 ix=1,3
3005      xd(ix)=xta(ix,1)+da*xx(ix)/dab
c  create ring vectors
          rvx(1)=xx(2)*dirtm(3,1)-xx(3)*dirtm(2,1)
          rvx(2)=xx(3)*dirtm(1,1)-xx(1)*dirtm(3,1)
          rvx(3)=xx(1)*dirtm(2,1)-xx(2)*dirtm(1,1)
          dist=sqrt(rvx(1)**2+rvx(2)**2+rvx(3)**2)
          do 3006 ix=1,3
3006      rvx(ix)=hh*rvx(ix)/dist
          rvy(1)=(xx(2)*rvx(3)-xx(3)*rvx(2))/dab
          rvy(2)=(xx(3)*rvx(1)-xx(1)*rvx(3))/dab
          rvy(3)=(xx(1)*rvx(2)-xx(2)*rvx(1))/dab
c now all triple points on the ring are searched
          ntrp=0
          ntrp2=0
          do 3390 iic=ilipa+1,ilipa+nipa(ia)
            ic=lipa(iic)
            if (ic.eq.ib) go to 3390
            rc=srad(ic)
            dabc=0.d0
            sp=0.d0
            do 3302 ix=1,3
              xxx=coord(ix,ic)-xd(ix)
              xta(ix,3)=coord(ix,ic)
              sp=sp+xxx*xx(ix)
3302        dabc=dabc+xxx**2
            dabc=sqrt(dabc)
            cosa=sp/dab/dabc
c -- tests added    ! JB Oct 2011
          If(Abs(cosa).GE.1.0d0) cosa=SIGN(1.0d0,cosa)
            sina=sqrt(1.d0-cosa*cosa)
          If(sina.LT.TollZero) GO TO 3390
            cj=(dabc*dabc+hh*hh-rc*rc)/(2*dabc*hh*sina)
            if (cj .lt. 1.d0) ntrp2=ntrp2+2
            if (cj .gt. 1.d0 .or. cj .lt. -1.d0) go to 3390
            sj=sqrt(1.d0-cj*cj)
c  create ring vectors
            do 3305 ix=1,3
3305        tvx(ix)=(xta(ix,3)-xd(ix))-cosa*dabc*xx(ix)/dab
            dist=sqrt(tvx(1)**2+tvx(2)**2+tvx(3)**2)
            do 3306 ix=1,3
3306        tvx(ix)=hh*tvx(ix)/dist
            tvy(1)=(xx(2)*tvx(3)-xx(3)*tvx(2))/dab
            tvy(2)=(xx(3)*tvx(1)-xx(1)*tvx(3))/dab
            tvy(3)=(xx(1)*tvx(2)-xx(2)*tvx(1))/dab
            do 3380 l=-1,1,2
              il=ntrp+1
              do 3310 ix=1,3
3310          trp(ix)=xd(ix)+cj*tvx(ix)+sj*tvy(ix)*l
              do 3320 ik=ilipa+1,ilipa+nipa(ia)
                k=lipa(ik)
                if(k.eq.ib .or. k.eq.ic) go to 3320
                dabck=0.d0
                do 3315 ix=1,3
3315            dabck=dabck+(trp(ix)-coord(ix,k))**2
                dabck=sqrt(dabck)
                if (dabck .lt. srad(k)) go to 3380
3320          continue
              ntrp=ntrp+1
              spx=0.d0
              spy=0.d0
              do 3322 ix=1,3
                spx=spx+rvx(ix)*(trp(ix)-xd(ix))
                spy=spy+rvy(ix)*(trp(ix)-xd(ix))
3322          continue
              phi=acos(spx/hh**2)
              if(spy.lt.0.d0) phi=-phi
              phiset(il)=phi
              sp=0.d0
              do 3324 ix=1,3
3324          sp=sp+(-spy*rvx(ix)+spx*rvy(ix))
     &               *(trp(ix)-xta(ix,3))
              iset(ntrp)=1
              if(sp.lt.0.d0) iset(ntrp)=-1
c if the triple point is new the corresponding surface patches are added
c
c midi use only triple points which fulfil the following checks
              dbc=0.d0
              dac=0.d0
              do  ix=1,3
               dbc=dbc+(coord(ix,ic)-coord(ix,ib))**2
               dac=dac+(coord(ix,ic)-coord(ix,ia))**2
              enddo
              dbc = sqrt(dbc)
              dac = sqrt(dac)
              if(srad(ia)+srad(ib)-rsolv .lt. dab) then
               go to 3380
              endif
              if(srad(ic)+srad(ib)-rsolv .lt. dbc) then
               go to 3380
              endif
              if(srad(ia)+srad(ic)-rsolv .lt. dac) then
               go to 3380
              endif
c end midi
              if(ic.le.ib) go to 3380
              do 3325 ix=1,3
                svx(ix)=xta(ix,2)-xta(ix,1)
                svy(ix)=xta(ix,3)-xta(ix,1)
3325          em(ix)=0.d0
              do 3330 iiia=1,3
                dist=0.d0
                do 3326 ix=1,3
                  xxx=xta(ix,iiia)-trp(ix)
                  dist=dist+xxx*xxx
3326            ee(ix,iiia)=xxx
                dist=sqrt(dist)
                do 3328 ix=1,3
                  xxx=ee(ix,iiia)/dist
                  em(ix)=em(ix)+xxx
3328            ee(ix,iiia)=xxx
3330          continue
              dist=0.d0
              spn=0.d0
              do 3334 iiia=1,3
                dist=dist+em(iiia)**2
                iiib=mod(iiia,3)+1
                iiic=6-iiia-iiib
                xxx=svx(iiib)*svy(iiic)-svx(iiic)*svy(iiib)
                spn=spn+xxx*(xta(iiia,1)-trp(iiia))
                yy(iiia)=xxx
                sp=0.d0
                do 3332 ix=1,3
3332            sp=sp+ee(ix,iiib)*ee(ix,iiic)
                cc(iiia)=sp
3334          ss(iiia)=sqrt(1.d0-sp*sp)
              sar=-pi
              dist=sqrt(dist)
              do 3338 iiia=1,3
                em(iiia)=em(iiia)/dist
                iiib=mod(iiia,3)+1
                iiic=6-iiia-iiib
3338          sar=sar+
     &        acos((cc(iiia)-cc(iiib)*cc(iiic))/ss(iiib)/ss(iiic))
              sar=sar*rsolv**2/3
c if segment is less than 1% of a normal segment, drop it
              if (sar.lt. .5d0/nspa) go to 3380
              do 3370 iiia=1,3
                dist=0.d0
                spa=0.d0
                do 3345 ix=1,3
                  xxx=1.4d0*em(ix)+ee(ix,iiia)
                  dist=dist+xxx*xxx
                  spa=spa+yy(ix)*xxx
3345            cc(ix)=xxx
                dist=1/sqrt(dist)
                fact=rsolv*dist
c if the segment center lies below the plane of the three atoms drop it
                if(spa*fact/spn .gt. 1) then
                  go to 3370
                end if
                nps=nps+1
                iii=ia
                if(iiia.eq.2) iii=ib
                if(iiia.eq.3) iii=ic
                iatsp(nps)=iii
                ar(nps)=sar
                spnn=0.d0
                do 3348 ix=1,3
                  cosurf(ix,nps)=-cc(ix)*dist
                  xsp(ix,nps)=trp(ix)+fact*cc(ix)
                  spnn=spnn+xsp(ix,nps)*cosurf(ix,nps)
3348            continue
                area=area+sar
                volume=volume+sar*spnn
3370          continue
3380        continue
3390      continue
c sort the set of triple points on the ring
          if(mod(ntrp,2).ne.0) then
             mes='COSMO problem: odd ntrp in consts.f'
             ierr=-1
             return
          endif
cc          if(ntrp.gt.18) then
          if(ntrp.gt.maxntrp) then   ! MM Aug 2008  now a parameter
             mes='COSMO problem: ntrp  too big in consts.f'
             ierr=-1
             return
          endif
          if(ntrp+ntrp2.eq.0) then
            phiset(1)=0
            phiset(2)=2*pi
            iset(1)=1
            iset(2)=-1
            ntrp=2
          end if

3400      ic=0
          do 3410 l=2,ntrp
          if (phiset(l).lt.phiset(l-1)) then
            phi=phiset(l)
            iii=iset(l)
            phiset(l)=phiset(l-1)
            iset(l)=iset(l-1)
            phiset(l-1)=phi
            iset(l-1)=iii
            ic=ic+1
          end if
3410      continue
          if(ic.gt.0) go to 3400
          if(iset(1) .eq. -1) then
            phiset(1)=phiset(1)+2*pi
            go to 3400
          end if
c now for each continuous section of the ring triangles are created
          sumphi=0.d0
          ips0=nps
          do 3600 l=2,ntrp,2
            k=l-1
            phiu=phiset(k)
            phio=phiset(l)
            nsa=(phio-phiu)/2/pi*nsring
            nsa=max(nsa+1,2)
            sumphi=sumphi+phio-phiu
            dp=(phio-phiu)/(nsa-1)
            do 3550 ich=1,2
            do 3550 ja=ich,nsa,2
              jb=max(ja-1,1)
              jc=min(ja+1,nsa)
              phi=phiu+(ja-1)*dp
              cphi=cos(phi)
              sphi=sin(phi)
              do 3545 ix=1,3
                ca=xd(ix)+(cphi*rvx(ix)+sphi*rvy(ix))*fz(ich)
                ca=ca+(xta(ix,ich)-ca)*yx(ich)
3545          xja(ix)=ca+(xta(ix,3-ich)-xta(ix,ich))*ddd*
     &                (1.d0-fz(ich))
              phi=phiu+(jb-1)*dp
              cphi=cos(phi)
              sphi=sin(phi)
              do 3546 ix=1,3
                ca=xd(ix)+cphi*rvx(ix)+sphi*rvy(ix)
3546          xjb(ix)=ca+(xta(ix,3-ich)-ca)*yx(3-ich)
              phi=phiu+(jc-1)*dp
              cphi=cos(phi)
              sphi=sin(phi)
              do 3547 ix=1,3
                ca=xd(ix)+cphi*rvx(ix)+sphi*rvy(ix)
3547          xjc(ix)=ca+(xta(ix,3-ich)-ca)*yx(3-ich)
              sp=0.d0
              d1=0.d0
              d2=0.d0
              do ix=1,3
                xx1=xjc(ix)-xjb(ix)
                xx2=xja(ix)-xjb(ix)
                d1=d1+xx1*xx1
                d2=d2+xx2*xx2
                sp=sp+xx1*xx2
              enddo
              dl=sqrt(d2-sp*sp/d1)
              sar=(jc-jb)*dp*hh*(1.d0-yx(3-ich))*dl/2
c if segment is less than 1% of a normal segment, drop it
              if (sar.lt. .5d0/nspa) go to 3550
              nps=nps+1
              nsab=nsab+1
              spn=0.d0
              dist=0.d0
              iatsp(nps)=ib
              if(ich.eq.2) iatsp(nps)=ia
              iat=iatsp(nps)
              do 3548 ix=1,3
                xsp(ix,nps)=(xja(ix)*0.5d0+xjb(ix)+xjc(ix))/2.5d0
                i2=mod(ix,3)+1
                i3=mod(i2,3)+1
                cosurf(ix,nps)=(xjc(i2)-xjb(i2))*(xja(i3)-xjb(i3))
     &                        -(xja(i2)-xjb(i2))*(xjc(i3)-xjb(i3))
                dist=dist+cosurf(ix,nps)**2
                spn=spn+cosurf(ix,nps)*(xsp(ix,nps)-coord(ix,iat))
3548          continue
              dist=1.d0/sqrt(dist)
              if(spn .lt. 0.d0) dist=-dist
              spn=0.d0
              do 3549 ix=1,3
                cosurf(ix,nps)=cosurf(ix,nps)*dist
3549          spn=spn+cosurf(ix,nps)*xsp(ix,nps)
              ar(nps)=sar
              volume=volume+ar(nps)*spn
              area=area+ar(nps)
3550        continue
3600      continue
          if (sumphi.gt.1.d-10) then
c           nip counts intersections between atoms
c           ia, ib denote intersecting atoms
c           nsab is the number of gridpoints on the intersection ring
c             of atoms ia, ib
c           ips0 is the last gridpoint at atom ia before the
c             intersection points
            nip=nip+1
            isude(1,nip)=ia
            isude(2,nip)=ib
            isude(3,nip)=nsab
            isude(4,nip)=ips0
            call ansude(ra-rsolv,rb-rsolv,dab,rsolv,
     &                aa,ab,aar,abr,aad,abd)
            sumphi=sumphi/(2*pi)
            sude(1,nip)=aar*sumphi
            sude(2,nip)=abr*sumphi
            sude(3,nip)=aad*sumphi
            sude(4,nip)=abd*sumphi
          endif
3700    continue
3800  ilipa=ilipa+nipa(ia)
      volume=volume/3

cds   save sude-stuff on file sude
      call writsude(nip,sude,isude,ierr,mes)
      if(ierr .ne.0 ) return

c size for outersphere is npspher
c size for innershpere is nps
c
3333  continue
c replace nps by npspher
      do 449 ips=1,npspher
          jps=nps+ips
          i=iatsp(ips)
          iatsp(jps)=i
          ri=srad(i)+(routf-1)*rsolv
          nar(jps)=nar(ips)
          nsetf(jps)=nsetf(ips)
          do 448 ix=1,3
c             cosurf(ix,jps)=cosurf(ix,ips)
              xsp(ix,jps)=coord(ix,i)+
     $                    (xsp(ix,ips)-coord(ix,i))*ri/(srad(i)-rsolv)
448       continue
          ar(jps)=ar(ips)*(ri/(srad(i)-rsolv))**2
449   continue

      if (nps .gt. maxnps) then
        write (mes,'(a,2i5,a)') ' nps is greater than maxnps (',
     &    nps,maxnps,') --choose bigger maxnps if possible'
        ierr=-1
        return
      endif
      npsd= nps + npspher
c
c.... end of surface construction
c
c.... save xsp to cosurf
      do l=1,npsd
          cosurf(1,l) = xsp(1,l)
          cosurf(2,l) = xsp(2,l)
          cosurf(3,l) = xsp(3,l)
      enddo
c---------------------------------------------------------------------c
c... write surface for debug purposes
c     call get_file(ifcx,mes)
c     if(ifcx .eq. -1) then
c       ierr=-1
c       return
c     endif
c     open(ifcx,file='cosurf',status='unknown',err=900)
c     write(ifcx,'(a)') '# Information about the cavity'
c     write(ifcx,'(a)') '# DISTORTED!!! atomic coordinates'
c     do i=1,natoms
c         write(ifcx,1000) (coord(j,i),j=1,3),nuc(i)
c     enddo
c
c    write segment information
c
c     write(ifcx,'(a)') '# Segment information'
c     write(ifcx,'(a)') '# Number of segments (nps)'
c     write(ifcx,*) nps
c     write(ifcx,'(a)') '# he="normal seg."; x=triple p. or ring'
c     write(ifcx,'(a)') '#  n   atom                position
c    &(X, Y, Z)[bohr]      area [bohr**2]'
c     do i=1,nps
c        if(i .le. npspher) then
c
c           points on spheres
c
c           write(ifcx,1001) i,iatsp(i),(xsp(j,i),j=1,3),
c    &      ar(i),'he'
c         else
c
c           ring and triple points
c
c           write(ifcx,1001) i,iatsp(i),(xsp(j,i),j=1,3),
c    &        ar(i),'x'
c         endif
c     enddo
c1000 format(3(f20.14,2x),4x,i5)
c1001 format(i5,i5,4f15.9,a5)
c     close(ifcx)
c---------------------------------------------------------------------c
c     deallocate(din, xsp, sude, nipa, lipa, isude, nn, stat=ierr)
c     if (ierr /= 0) then
c        mes='COSMO consts: cannot deallocate mem.'
c        ierr=-1
c        return
c     endif

      return
  900 mes = 'COSMO consts: cannot open file cosurf'
      ierr=-1
      return

      end
c =========================================================================
       subroutine coschol1(a,n,info)
c written by A. Klamt, 12/97
cas modifed version !

c this routine performs a Cholesky factorization
c input:   a =  packed lower triangle of a
c               symmetric positive definite n*n matrix
c output:  a =  lower triangle of Cholesky matrix ( inverse pivot elements)

       implicit none

       double precision a(*), summe
       integer n,info,i,j,k,indk,indi,is,kk

       indk=0
       info = 0
       do k = 1, n
         indi=indk
         kk=k+indk
         do i = k, n
            summe = 0.d0
            do j = 1, k-1
               summe = summe + a(j+indi)*a(j+indk)
            enddo
            summe = a(k+indi) - summe
            if (i.eq.k) then
               if (summe.lt.0.0d0 ) then
                  info = -1
                  summe = a(kk)
               endif
               a(kk) = 1.d0/sqrt(summe)
            else
               a(k+indi) = summe*a(kk)
            endif
            indi=indi+i
         enddo
         indk=indk+k
       enddo

       end
       subroutine coschol2(a,x,y,n,scal)
c written by A. Klamt, 12/97
cas modifed version !

c this routine solves the linear system Cx = y based on
c Cholesky factorization
c input:   a =  upper triangle of Cholesky matrix
c input:   y =  vector of length n
c output:  x =  vector of length n

cas    modified for inclusion of a scaling factor scal, i.e.  a*x = scal*y

       implicit none

       double precision a(*), x(*), y(*), summe, scal
       integer n,i,k,indk,indi

       indk=0
       do k = 1, n
         summe=scal*y(k)
         do i = k-1, 1,-1
           summe=summe-a(i+indk)*x(i)
         enddo
         x(k)=summe*a(k+indk)
         indk=indk+k
       enddo

       indk=n*(n-1)/2
       do k = n,1,-1
         summe=x(k)
         indi=indk
         do i = k+1, n
           indi=indi+i-1
           summe=summe-a(k+indi)*x(i)
         enddo
         x(k)=summe*a(k+indk)
         indk=indk-k+1
       enddo

       end
      subroutine cosmo_oc(ediel,phi,qcos,a1mat,natoms,phio,
     &                   fepsi,nps,npspher,qsum,qsumo,de,
     &                   qcosc,phic,ierr,mes,qcoso,a2mat,a3mat,
     &                   jobname,lenj)
c----------------------------------------------------------------------c
c     outlying charge correction
c     calls:
c            get_file
c            coschol1
c            coschol2
c----------------------------------------------------------------------c
c
      implicit real*8 (a-h,o-z)
c

      character*256 jobname
      character*100 mes
      dimension phi(nps), qcos(nps), phio(npspher)
      dimension a1mat(nps*(nps+1)/2)
      dimension qcosc(nps), phic(nps)
c     double precision, allocatable :: qcoso(:)
c     double precision, allocatable :: a2mat(:,:),a3mat(:)
      dimension qcoso(npspher),a2mat(npspher,nps),
     $          a3mat(npspher*(npspher+1)/2)
c

      parameter(zero=0.d0, half=0.5d0)

c----------------------------------------------------------------------c
c     allocate memory
c     allocate(qcoso(npspher),
c    &      a2mat(npspher,nps),a3mat(npspher*(npspher+1)/2), stat=ierr)
c     if(ierr .ne. 0) then
c        ierr=-1
c        mes='COSMO cosmo_oc: cannot alloc. memory'
c        return
c     endif
c----------------------------------------------------------------------c
      call zeroit(qcoso,npspher)
      call zeroit(a2mat,npspher*nps)
      call zeroit(a3mat,npspher*(npspher+1)/2)
      call zeroit(a1mat,nps*(nps+1)/2)
      call get_file(ifca,mes)
      if(ifca.eq.1) return

c     read lower triangle of A1 matrix from file
      open(ifca,file=jobname(1:lenj)//'.amat',form='unformatted',
     &     status='old',err=9010)
c     read(ifca) (a1mat(i),i=1,nps*(nps+1)/2)
      read(ifca) a1mat
      close(ifca,status='delete')

c     read a2mat
      open(ifca,file=jobname(1:lenj)//'.a2mat',form='unformatted',
     &     status='old',err=9020)
      do i=1,nps
          read(ifca) (a2mat(j,i),j=1,npspher)
      enddo
      close(ifca,status='delete')

c     read a3mat
      open(ifca,file=jobname(1:lenj)//'.a3mat',form='unformatted',
     &     status='old',err=9030)
      do i=1,npspher
        read(ifca) (a3mat(i+j*(j-1)/2),j=i,npspher)
      enddo
      close(ifca,status='delete')

c----------------------------------------------------------------------c
c... Add the inner screening charges
      do i=1,nps
          do j=1,npspher
              phio(j)=phio(j)+a2mat(j,i)*qcos(i)
          enddo
      enddo

c... calc. outer screening charges qcoso = - A3**(-1) * phio

      call coschol1(a3mat,npspher,info)
      if (info.lt.0)  then
         mes='COSMO cosmo_oc: a3mat not positive definite'
         ierr=-1
         return
      endif
      scal=-1.d0
      call coschol2(a3mat,qcoso,phio,npspher,scal)

c----------------------------------------------------------------------c
c     calculate corrected charges (qcosc)
      do i=1, npspher
          qcosc(i)=qcos(i)+qcoso(i)
      enddo
      do i=npspher+1, nps
          qcosc(i)=qcos(i)
      enddo

c----------------------------------------------------------------------c
c     calculate  total screening charge
      qsum  = zero
      qsumo = zero
      do i = 1, nps
          qsum  = qsum  + qcos(i)
      enddo
      do i = 1, npspher
          qsumo = qsumo + qcoso(i)
      enddo
c----------------------------------------------------------------------c
      de = zero
      call zeroit(phic,nps)

      ia1=0
      do i=1,nps
          do j=1,i-1
              ia1=ia1+1
              phic(i) = phic(i) - qcosc(j) * a1mat(ia1)
              phic(j) = phic(j) - qcosc(i) * a1mat(ia1)
              de = de - qcosc(i) * qcosc(j) * a1mat(ia1)
          enddo
          ia1=ia1+1
          phic(i) = phic(i) - qcosc(i) * a1mat(ia1)
          de = de - half * qcosc(i)**2 * a1mat(ia1)
      enddo
      de  =  fepsi * de - ediel
c---------------------------------------------------------------------c
c     deallocate memory
c     deallocate(qcoso,a2mat,a3mat, stat=ierr)
c     if(ierr .ne. 0) then
c        ierr=-1
c        mes='COSMO cosmo_oc: cannot dealloc. memory'
c        return
c     endif
c---------------------------------------------------------------------c

      return

 9000 mes='COSMO cosmo_oc: cannot open cosmo output file '
      ierr=-1
      return
 9010 mes='COSMO cosmo_oc: cannot open cosmo file '
      ierr=-1
      return
 9020 mes='COSMO cosmo_oc: cannot open cosmo file a2mat.tmp'
      ierr=-1
      return
 9030 mes='COSMO cosmo_oc: cannot open cosmo file a3mat.tmp'
      ierr=-1
      return
      end

      subroutine dvfill(nppa,dirvec,ierr,mes)
      implicit real*8 (a-h,o-z)

c ---------------------------------------------------------------------
c     construction of a mesh on a unit sphere
c ---------------------------------------------------------------------
c
      dimension dirvec(3,nppa)
      integer fset(3,20), kset(2,30), ierr
      character*100 mes

      data kset/ 1, 2, 1, 3, 1, 4, 1, 5, 1, 6,
     1            12,11,12,10,12, 9,12, 8,12, 7,
     2             2, 3, 3, 4, 4, 5, 5, 6, 6, 2,
     3             7, 8, 8, 9, 9,10,10,11,11, 7,
     4             2,7,7,3,3,8,8,4,4,9,9,5,5,10,10,6,6,11,11,2/
      data fset/ 1, 2, 3, 1, 3, 4, 1, 4, 5, 1, 5, 6, 1, 6, 2,
     1            12,11,10,12,10, 9,12, 9, 8,12, 8, 7,12, 7,11,
     2             2, 3, 7, 3, 4, 8, 4, 5, 9, 5, 6,10, 6, 2,11,
     3             7, 8, 3, 8, 9, 4, 9,10, 5,10,11, 6,11, 7, 2/

      ierr = 0
      dirvec (1,1) =  -1.d0
      dirvec (2,1) =   0.d0
      dirvec (3,1) =   0.d0
      nd=1
      r=sqrt(0.8d0)
      h=sqrt(0.2d0)
      do 10 i= -1,1,2
         do 10 j= 1,5
            nd=nd+1
            beta=1.d0+ j*1.25663706d0 + (i+1)*0.3141593d0
            dirvec(2,nd)=r*cos(beta)
            dirvec(3,nd)=r*sin(beta)
            dirvec(1,nd)=i*h
   10 continue
      dirvec (2,12) =  0.d0
      dirvec (3,12) =  0.d0
      dirvec (1,12) =  1.d0
      nd=12
c  nppa=10*3**k*4**l+2
      m=(nppa-2)/10
      do 20 k=0,10
         if ((m/3)*3 .ne. m) go to 30
   20 m=m/3
   30 do 40 l=0,10
         if ((m/4)*4 .ne. m) go to 50
   40 m=m/4
   50 if (10*3**k*4**l+2 .ne. nppa) then
       mes='COSMO: value of nspa, nsph, or nppa not allowed:
     & it must be 10*3**k*4**l+2'
      ierr=-1
      return
      endif
      kh=k/2
      m=2**l*3**kh
c create on each edge 2**l*3**kh-1 new points
      do 70 i=1,30
         na=kset(1,i)
         nb=kset(2,i)
         do 70 j=1,m-1
            nd=nd+1
            do 60 ix=1,3
   60       dirvec(ix,nd)=dirvec(ix,na)*(m-j)+dirvec(ix,nb)*j
   70 continue
c create points within each triangle
      do 90 i=1,20
         na=fset(1,i)
         nb=fset(2,i)
         nc=fset(3,i)
         do 90 j1=1,m-1
            do 90 j2=1,m-j1-1
               nd=nd+1
               do 80 ix=1,3
   80          dirvec(ix,nd)=dirvec(ix,na)*(m-j1-j2)
     1                     +dirvec(ix,nb)*j1+dirvec(ix,nc)*j2
   90 continue
      if (k .eq. 2*kh) go to 140
c create to additional subgrids
      t=1.0d0/3.0d0
      do 110 i=1,20
         na=fset(1,i)
         nb=fset(2,i)
         nc=fset(3,i)
         do 110 j1=0,m-1
            do 110 j2=0,m-j1-1
               nd=nd+1
               do 100 ix=1,3
  100          dirvec(ix,nd)=dirvec(ix,na)*(m-j1-j2-2.0d0*t)
     1                 +dirvec(ix,nb)*(j1+t)+dirvec(ix,nc)*(j2+t)
  110 continue
      t=2.0d0/3.0d0
      do 130 i=1,20
         na=fset(1,i)
         nb=fset(2,i)
         nc=fset(3,i)
         do 130 j1=0,m-2
            do 130 j2=0,m-j1-2
               nd=nd+1
               do 120 ix=1,3
  120          dirvec(ix,nd)=dirvec(ix,na)*(m-j1-j2-2*t)
     1                  +dirvec(ix,nb)*(j1+t)+dirvec(ix,nc)*(j2+t)
  130 continue

c     --- normalize all vectors
  140 do 170 i=1,nppa
         dist=0.d0
         do 150 ix=1,3
  150    dist=dist+dirvec(ix,i)**2
         dist=1./sqrt(dist)
         do 160 ix=1,3
  160    dirvec(ix,i)=dirvec(ix,i)*dist
c      print *,'di ',i,dirvec(1,i),dirvec(2,i),dirvec(3,i)
  170 continue

      return
      end
c
c----------------------------------------------------------------------
      subroutine fill_grd1082(basgrd)
c----------------------------------------------------------------------
c     fills array with 1082 basis grid points
      double precision basgrd
      dimension    basgrd(3*1082)
      basgrd(1)=  -0.954912781715393d0
      basgrd(2)=   0.296886593103409d0
      basgrd(3)=   0.000084500003140d0
      basgrd(4)=  -0.593343913555145d0
      basgrd(5)=  -0.409820109605789d0
      basgrd(6)=   0.692813515663147d0
      basgrd(7)=  -0.672271013259888d0
      basgrd(8)=  -0.664492487907410d0
      basgrd(9)=  -0.326345413923264d0
      basgrd(10)=  -0.410722494125366d0
      basgrd(11)=   0.181889593601227d0
      basgrd(12)=  -0.893433392047882d0
      basgrd(13)=  -0.170661404728890d0
      basgrd(14)=   0.958825409412384d0
      basgrd(15)=  -0.226998805999756d0
      basgrd(16)=  -0.283331006765366d0
      basgrd(17)=   0.593901574611664d0
      basgrd(18)=   0.752996981143951d0
      basgrd(19)=   0.170589506626129d0
      basgrd(20)=  -0.958851218223572d0
      basgrd(21)=   0.226944103837013d0
      basgrd(22)=   0.283320993185043d0
      basgrd(23)=  -0.593910276889801d0
      basgrd(24)=  -0.752993881702423d0
      basgrd(25)=   0.593272507190704d0
      basgrd(26)=   0.409606099128723d0
      basgrd(27)=  -0.693001210689545d0
      basgrd(28)=   0.672566473484039d0
      basgrd(29)=   0.664503097534180d0
      basgrd(30)=   0.325714498758316d0
      basgrd(31)=   0.410804688930512d0
      basgrd(32)=  -0.181782394647598d0
      basgrd(33)=   0.893417418003082d0
      basgrd(34)=   0.954448103904724d0
      basgrd(35)=  -0.298376888036728d0
      basgrd(36)=  -0.000033500000427d0
      basgrd(37)=  -0.973186671733856d0
      basgrd(38)=   0.183344498276711d0
      basgrd(39)=   0.138897001743317d0
      basgrd(40)=  -0.959407925605774d0
      basgrd(41)=   0.060520999133587d0
      basgrd(42)=   0.275451689958572d0
      basgrd(43)=  -0.911280691623688d0
      basgrd(44)=  -0.067142099142075d0
      basgrd(45)=   0.406275212764740d0
      basgrd(46)=  -0.830864727497101d0
      basgrd(47)=  -0.191389992833138d0
      basgrd(48)=   0.522526204586029d0
      basgrd(49)=  -0.723877489566803d0
      basgrd(50)=  -0.306913793087006d0
      basgrd(51)=   0.617904007434845d0
      basgrd(52)=  -0.988995671272278d0
      basgrd(53)=   0.132551506161690d0
      basgrd(54)=  -0.065708801150322d0
      basgrd(55)=  -0.990818977355957d0
      basgrd(56)=  -0.040664300322533d0
      basgrd(57)=  -0.128934293985367d0
      basgrd(58)=  -0.957531392574310d0
      basgrd(59)=  -0.216375693678856d0
      basgrd(60)=  -0.190565198659897d0
      basgrd(61)=  -0.890303909778595d0
      basgrd(62)=  -0.383900105953217d0
      basgrd(63)=  -0.244907498359680d0
      basgrd(64)=  -0.793691396713257d0
      basgrd(65)=  -0.534194409847260d0
      basgrd(66)=  -0.291016012430191d0
      basgrd(67)=  -0.936085820198059d0
      basgrd(68)=   0.302436888217926d0
      basgrd(69)=  -0.179653406143189d0
      basgrd(70)=  -0.886284470558167d0
      basgrd(71)=   0.295561999082565d0
      basgrd(72)=  -0.356571108102799d0
      basgrd(73)=  -0.803864002227783d0
      basgrd(74)=   0.279391497373581d0
      basgrd(75)=  -0.525112390518189d0
      basgrd(76)=  -0.693184196949005d0
      basgrd(77)=   0.254003196954727d0
      basgrd(78)=  -0.674520611763001d0
      basgrd(79)=  -0.559542298316956d0
      basgrd(80)=   0.222004100680351d0
      basgrd(81)=  -0.798515319824219d0
      basgrd(82)=  -0.887385785579681d0
      basgrd(83)=   0.458791106939316d0
      basgrd(84)=  -0.045356098562479d0
      basgrd(85)=  -0.789565980434418d0
      basgrd(86)=   0.606883406639099d0
      basgrd(87)=  -0.090984500944614d0
      basgrd(88)=  -0.661834180355072d0
      basgrd(89)=   0.737595498561859d0
      basgrd(90)=  -0.133896499872208d0
      basgrd(91)=  -0.511188387870789d0
      basgrd(92)=   0.842123806476593d0
      basgrd(93)=  -0.171796098351479d0
      basgrd(94)=  -0.344993203878403d0
      basgrd(95)=   0.916490077972412d0
      basgrd(96)=  -0.202547803521156d0
      basgrd(97)=  -0.910504400730133d0
      basgrd(98)=   0.384692400693893d0
      basgrd(99)=   0.151636406779289d0
      basgrd(100)=  -0.834995508193970d0
      basgrd(101)=   0.460652798414230d0
      basgrd(102)=   0.300967603921890d0
      basgrd(103)=  -0.728512585163117d0
      basgrd(104)=   0.522332370281220d0
      basgrd(105)=   0.443213492631912d0
      basgrd(106)=  -0.596495211124420d0
      basgrd(107)=   0.566093325614929d0
      basgrd(108)=   0.568974316120148d0
      basgrd(109)=  -0.446372091770172d0
      basgrd(110)=   0.589471399784088d0
      basgrd(111)=   0.673257291316986d0
      basgrd(112)=   0.936357080936432d0
      basgrd(113)=  -0.301573485136032d0
      basgrd(114)=   0.179690897464752d0
      basgrd(115)=   0.886357188224793d0
      basgrd(116)=  -0.295602709054947d0
      basgrd(117)=   0.356356590986252d0
      basgrd(118)=   0.803808987140656d0
      basgrd(119)=  -0.279597193002701d0
      basgrd(120)=   0.525087177753449d0
      basgrd(121)=   0.692752897739410d0
      basgrd(122)=  -0.254478007555008d0
      basgrd(123)=   0.674784719944000d0
      basgrd(124)=   0.559929907321930d0
      basgrd(125)=  -0.221578702330589d0
      basgrd(126)=   0.798361718654633d0
      basgrd(127)=   0.989067614078522d0
      basgrd(128)=  -0.132420599460602d0
      basgrd(129)=   0.064884603023529d0
      basgrd(130)=   0.990777611732483d0
      basgrd(131)=   0.040625199675560d0
      basgrd(132)=   0.129264593124390d0
      basgrd(133)=   0.957545399665833d0
      basgrd(134)=   0.216209098696709d0
      basgrd(135)=   0.190684005618095d0
      basgrd(136)=   0.890233695507050d0
      basgrd(137)=   0.383737802505493d0
      basgrd(138)=   0.245416507124901d0
      basgrd(139)=   0.793765425682068d0
      basgrd(140)=   0.534071922302246d0
      basgrd(141)=   0.291039109230042d0
      basgrd(142)=   0.973104476928711d0
      basgrd(143)=  -0.183540403842926d0
      basgrd(144)=  -0.139213800430298d0
      basgrd(145)=   0.959453403949738d0
      basgrd(146)=  -0.060241401195526d0
      basgrd(147)=  -0.275354504585266d0
      basgrd(148)=   0.911367118358612d0
      basgrd(149)=   0.067354299128056d0
      basgrd(150)=  -0.406046211719513d0
      basgrd(151)=   0.831071913242340d0
      basgrd(152)=   0.191786095499992d0
      basgrd(153)=  -0.522051274776459d0
      basgrd(154)=   0.723503291606903d0
      basgrd(155)=   0.306832194328308d0
      basgrd(156)=  -0.618382573127747d0
      basgrd(157)=   0.910289704799652d0
      basgrd(158)=  -0.385277211666107d0
      basgrd(159)=  -0.151440292596817d0
      basgrd(160)=   0.835057675838471d0
      basgrd(161)=  -0.460336208343506d0
      basgrd(162)=  -0.301279306411743d0
      basgrd(163)=   0.728647470474243d0
      basgrd(164)=  -0.522065818309784d0
      basgrd(165)=  -0.443305909633637d0
      basgrd(166)=   0.596837103366852d0
      basgrd(167)=  -0.565697312355042d0
      basgrd(168)=  -0.569009780883789d0
      basgrd(169)=   0.446017503738403d0
      basgrd(170)=  -0.589910209178925d0
      basgrd(171)=  -0.673107981681824d0
      basgrd(172)=   0.887485086917877d0
      basgrd(173)=  -0.458539605140686d0
      basgrd(174)=   0.045953001827002d0
      basgrd(175)=   0.789576709270477d0
      basgrd(176)=  -0.606891810894013d0
      basgrd(177)=   0.090834997594357d0
      basgrd(178)=   0.661794602870941d0
      basgrd(179)=  -0.737673819065094d0
      basgrd(180)=   0.133660793304443d0
      basgrd(181)=   0.511098921298981d0
      basgrd(182)=  -0.842231273651123d0
      basgrd(183)=   0.171535298228264d0
      basgrd(184)=   0.345062196254730d0
      basgrd(185)=  -0.916366696357727d0
      basgrd(186)=   0.202988103032112d0
      basgrd(187)=  -0.665464222431183d0
      basgrd(188)=  -0.499639898538590d0
      basgrd(189)=   0.554542481899262d0
      basgrd(190)=  -0.717494308948517d0
      basgrd(191)=  -0.574116408824921d0
      basgrd(192)=   0.394452005624771d0
      basgrd(193)=  -0.746370017528534d0
      basgrd(194)=  -0.629360973834992d0
      basgrd(195)=   0.216417506337166d0
      basgrd(196)=  -0.747793495655060d0
      basgrd(197)=  -0.663207828998566d0
      basgrd(198)=   0.030988000333309d0
      basgrd(199)=  -0.721587419509888d0
      basgrd(200)=  -0.675301313400269d0
      basgrd(201)=  -0.152577102184296d0
      basgrd(202)=  -0.683451771736145d0
      basgrd(203)=  -0.558005511760712d0
      basgrd(204)=  -0.470662713050842d0
      basgrd(205)=  -0.673277795314789d0
      basgrd(206)=  -0.429748207330704d0
      basgrd(207)=  -0.601675629615784d0
      basgrd(208)=  -0.639551877975464d0
      basgrd(209)=  -0.283360004425049d0
      basgrd(210)=  -0.714619219303131d0
      basgrd(211)=  -0.581798911094666d0
      basgrd(212)=  -0.127481803297997d0
      basgrd(213)=  -0.803279817104340d0
      basgrd(214)=  -0.503554403781891d0
      basgrd(215)=   0.029554000124335d0
      basgrd(216)=  -0.863457918167114d0
      basgrd(217)=  -0.401748090982437d0
      basgrd(218)=   0.354923099279404d0
      basgrd(219)=  -0.844173014163971d0
      basgrd(220)=  -0.379645586013794d0
      basgrd(221)=   0.519598901271820d0
      basgrd(222)=  -0.765432000160217d0
      basgrd(223)=  -0.344034403562546d0
      basgrd(224)=   0.669683873653412d0
      basgrd(225)=  -0.658151805400848d0
      basgrd(226)=  -0.295150607824326d0
      basgrd(227)=   0.796346127986908d0
      basgrd(228)=  -0.527938425540924d0
      basgrd(229)=  -0.236370101571083d0
      basgrd(230)=   0.893391907215118d0
      basgrd(231)=  -0.382073402404785d0
      basgrd(232)=  -0.210144594311714d0
      basgrd(233)=   0.976298391819000d0
      basgrd(234)=  -0.051775500178337d0
      basgrd(235)=  -0.243006497621536d0
      basgrd(236)=   0.961418390274048d0
      basgrd(237)=   0.128927901387215d0
      basgrd(238)=  -0.268602401018143d0
      basgrd(239)=   0.912409901618958d0
      basgrd(240)=   0.308805495500565d0
      basgrd(241)=  -0.284378707408905d0
      basgrd(242)=   0.831273615360260d0
      basgrd(243)=   0.477611690759659d0
      basgrd(244)=  -0.289085388183594d0
      basgrd(245)=   0.723142921924591d0
      basgrd(246)=   0.627290904521942d0
      basgrd(247)=  -0.372728109359741d0
      basgrd(248)=   0.448902994394302d0
      basgrd(249)=   0.812132894992828d0
      basgrd(250)=  -0.451533406972885d0
      basgrd(251)=   0.285832405090332d0
      basgrd(252)=   0.845232188701630d0
      basgrd(253)=  -0.517620921134949d0
      basgrd(254)=   0.109313897788525d0
      basgrd(255)=   0.848598301410675d0
      basgrd(256)=  -0.564411103725433d0
      basgrd(257)=  -0.070982798933983d0
      basgrd(258)=   0.822436273097992d0
      basgrd(259)=  -0.589245676994324d0
      basgrd(260)=  -0.246764004230499d0
      basgrd(261)=   0.769348502159119d0
      basgrd(262)=   0.210013896226883d0
      basgrd(263)=  -0.976330220699310d0
      basgrd(264)=   0.051706798374653d0
      basgrd(265)=   0.243146494030953d0
      basgrd(266)=  -0.961376428604126d0
      basgrd(267)=  -0.128977805376053d0
      basgrd(268)=   0.268418997526169d0
      basgrd(269)=  -0.912491202354431d0
      basgrd(270)=  -0.308724999427795d0
      basgrd(271)=   0.284179806709290d0
      basgrd(272)=  -0.831273913383484d0
      basgrd(273)=  -0.477729588747025d0
      basgrd(274)=   0.289351493120194d0
      basgrd(275)=  -0.723097681999207d0
      basgrd(276)=  -0.627220511436462d0
      basgrd(277)=   0.372696489095688d0
      basgrd(278)=  -0.448787391185761d0
      basgrd(279)=  -0.812211275100708d0
      basgrd(280)=   0.451811999082565d0
      basgrd(281)=  -0.285763502120972d0
      basgrd(282)=  -0.845106601715088d0
      basgrd(283)=   0.517483890056610d0
      basgrd(284)=  -0.109470598399639d0
      basgrd(285)=  -0.848661601543427d0
      basgrd(286)=   0.564076602458954d0
      basgrd(287)=   0.071169398725033d0
      basgrd(288)=  -0.822649717330933d0
      basgrd(289)=   0.589564383029938d0
      basgrd(290)=   0.246520698070526d0
      basgrd(291)=  -0.769182205200195d0
      basgrd(292)=   0.665301084518433d0
      basgrd(293)=   0.499959886074066d0
      basgrd(294)=  -0.554449677467346d0
      basgrd(295)=   0.717881023883820d0
      basgrd(296)=   0.573740422725678d0
      basgrd(297)=  -0.394295394420624d0
      basgrd(298)=   0.746526777744293d0
      basgrd(299)=   0.629135191440582d0
      basgrd(300)=  -0.216533303260803d0
      basgrd(301)=   0.747379124164581d0
      basgrd(302)=   0.663678526878357d0
      basgrd(303)=  -0.030907399952412d0
      basgrd(304)=   0.721801519393921d0
      basgrd(305)=   0.675050497055054d0
      basgrd(306)=   0.152674302458763d0
      basgrd(307)=   0.683133423328400d0
      basgrd(308)=   0.558319985866547d0
      basgrd(309)=   0.470752090215683d0
      basgrd(310)=   0.673569321632385d0
      basgrd(311)=   0.429404705762863d0
      basgrd(312)=   0.601594626903534d0
      basgrd(313)=   0.639734387397766d0
      basgrd(314)=   0.283428102731705d0
      basgrd(315)=   0.714428722858429d0
      basgrd(316)=   0.581489503383637d0
      basgrd(317)=   0.127562701702118d0
      basgrd(318)=   0.803490996360779d0
      basgrd(319)=   0.503574371337891d0
      basgrd(320)=  -0.029718199744821d0
      basgrd(321)=   0.863440573215485d0
      basgrd(322)=   0.401416391134262d0
      basgrd(323)=  -0.354910999536514d0
      basgrd(324)=   0.844335913658142d0
      basgrd(325)=   0.379967510700226d0
      basgrd(326)=  -0.519627511501312d0
      basgrd(327)=   0.765252888202667d0
      basgrd(328)=   0.344023495912552d0
      basgrd(329)=  -0.669562578201294d0
      basgrd(330)=   0.658280909061432d0
      basgrd(331)=   0.294941395521164d0
      basgrd(332)=  -0.796432375907898d0
      basgrd(333)=   0.527925193309784d0
      basgrd(334)=   0.236476600170136d0
      basgrd(335)=  -0.893398106098175d0
      basgrd(336)=   0.381993085145950d0
      basgrd(337)=  -0.495673507452011d0
      basgrd(338)=  -0.558793604373932d0
      basgrd(339)=   0.664874076843262d0
      basgrd(340)=  -0.379798889160156d0
      basgrd(341)=  -0.692314922809601d0
      basgrd(342)=   0.613557577133179d0
      basgrd(343)=  -0.248228400945664d0
      basgrd(344)=  -0.804310977458954d0
      basgrd(345)=   0.539876282215118d0
      basgrd(346)=  -0.107589401304722d0
      basgrd(347)=  -0.887595474720001d0
      basgrd(348)=   0.447882503271103d0
      basgrd(349)=   0.033585701137781d0
      basgrd(350)=  -0.939008593559265d0
      basgrd(351)=   0.342249602079392d0
      basgrd(352)=   0.017524400725961d0
      basgrd(353)=  -0.990307927131653d0
      basgrd(354)=   0.137779399752617d0
      basgrd(355)=  -0.139831602573395d0
      basgrd(356)=  -0.989291489124298d0
      basgrd(357)=   0.041827298700809d0
      basgrd(358)=  -0.295368403196335d0
      basgrd(359)=  -0.953669428825378d0
      basgrd(360)=  -0.057203300297260d0
      basgrd(361)=  -0.440783113241196d0
      basgrd(362)=  -0.884351611137390d0
      basgrd(363)=  -0.153728693723679d0
      basgrd(364)=  -0.567758023738861d0
      basgrd(365)=  -0.786244094371796d0
      basgrd(366)=  -0.243866994976997d0
      basgrd(367)=  -0.543661594390869d0
      basgrd(368)=  -0.713239014148712d0
      basgrd(369)=  -0.442404985427856d0
      basgrd(370)=  -0.394345492124558d0
      basgrd(371)=  -0.739198982715607d0
      basgrd(372)=  -0.545963823795319d0
      basgrd(373)=  -0.228054404258728d0
      basgrd(374)=  -0.739633202552795d0
      basgrd(375)=  -0.633193373680115d0
      basgrd(376)=  -0.053753599524498d0
      basgrd(377)=  -0.714144706726074d0
      basgrd(378)=  -0.697931110858917d0
      basgrd(379)=   0.118531301617622d0
      basgrd(380)=  -0.664576113224030d0
      basgrd(381)=  -0.737759411334992d0
      basgrd(382)=   0.170413896441460d0
      basgrd(383)=  -0.494423687458038d0
      basgrd(384)=  -0.852352201938629d0
      basgrd(385)=   0.049297798424959d0
      basgrd(386)=  -0.375985592603684d0
      basgrd(387)=  -0.925313174724579d0
      basgrd(388)=  -0.075535997748375d0
      basgrd(389)=  -0.242223098874092d0
      basgrd(390)=  -0.967275679111481d0
      basgrd(391)=  -0.198252394795418d0
      basgrd(392)=  -0.100528903305531d0
      basgrd(393)=  -0.974982023239136d0
      basgrd(394)=  -0.311620593070984d0
      basgrd(395)=   0.042380798608065d0
      basgrd(396)=  -0.949261009693146d0
      basgrd(397)=  -0.247370406985283d0
      basgrd(398)=   0.244400799274445d0
      basgrd(399)=  -0.937590599060059d0
      basgrd(400)=  -0.073338001966476d0
      basgrd(401)=   0.299904495477676d0
      basgrd(402)=  -0.951146125793457d0
      basgrd(403)=   0.107764601707459d0
      basgrd(404)=   0.346581608057022d0
      basgrd(405)=  -0.931809008121491d0
      basgrd(406)=   0.284937202930450d0
      basgrd(407)=   0.380385994911194d0
      basgrd(408)=  -0.879839301109314d0
      basgrd(409)=   0.448488712310791d0
      basgrd(410)=   0.401179701089859d0
      basgrd(411)=  -0.798694372177124d0
      basgrd(412)=   0.495634704828262d0
      basgrd(413)=   0.558957517147064d0
      basgrd(414)=  -0.664765179157257d0
      basgrd(415)=   0.379453599452972d0
      basgrd(416)=   0.692484617233276d0
      basgrd(417)=  -0.613579690456390d0
      basgrd(418)=   0.248121395707130d0
      basgrd(419)=   0.804369628429413d0
      basgrd(420)=  -0.539838194847107d0
      basgrd(421)=   0.107968002557755d0
      basgrd(422)=   0.887540698051453d0
      basgrd(423)=  -0.447900086641312d0
      basgrd(424)=  -0.033889401704073d0
      basgrd(425)=   0.938962101936340d0
      basgrd(426)=  -0.342347413301468d0
      basgrd(427)=  -0.017252299934626d0
      basgrd(428)=   0.990292072296143d0
      basgrd(429)=  -0.137927204370499d0
      basgrd(430)=   0.139438897371292d0
      basgrd(431)=   0.989346921443939d0
      basgrd(432)=  -0.041826099157333d0
      basgrd(433)=   0.295437812805176d0
      basgrd(434)=   0.953644812107086d0
      basgrd(435)=   0.057254400104284d0
      basgrd(436)=   0.441056698560715d0
      basgrd(437)=   0.884218871593475d0
      basgrd(438)=   0.153707295656204d0
      basgrd(439)=   0.568026602268219d0
      basgrd(440)=   0.785977423191071d0
      basgrd(441)=   0.244100898504257d0
      basgrd(442)=   0.543357670307159d0
      basgrd(443)=   0.713251709938049d0
      basgrd(444)=   0.442757785320282d0
      basgrd(445)=   0.394012898206711d0
      basgrd(446)=   0.739372313022614d0
      basgrd(447)=   0.545969307422638d0
      basgrd(448)=   0.227988794445992d0
      basgrd(449)=   0.739639282226563d0
      basgrd(450)=   0.633210003376007d0
      basgrd(451)=   0.054225899279118d0
      basgrd(452)=   0.714066624641419d0
      basgrd(453)=   0.697974503040314d0
      basgrd(454)=  -0.118776001036167d0
      basgrd(455)=   0.664754688739777d0
      basgrd(456)=   0.737559199333191d0
      basgrd(457)=  -0.170168802142143d0
      basgrd(458)=   0.494537889957428d0
      basgrd(459)=   0.852334976196289d0
      basgrd(460)=  -0.049561198800802d0
      basgrd(461)=   0.376092612743378d0
      basgrd(462)=   0.925255715847015d0
      basgrd(463)=   0.075599297881126d0
      basgrd(464)=   0.242237105965614d0
      basgrd(465)=   0.967267215251923d0
      basgrd(466)=   0.198477596044540d0
      basgrd(467)=   0.100399896502495d0
      basgrd(468)=   0.974949479103088d0
      basgrd(469)=   0.311431288719177d0
      basgrd(470)=  -0.042502298951149d0
      basgrd(471)=   0.949317693710327d0
      basgrd(472)=   0.247509807348251d0
      basgrd(473)=  -0.244682192802429d0
      basgrd(474)=   0.937480390071869d0
      basgrd(475)=   0.072844602167606d0
      basgrd(476)=  -0.299759298563004d0
      basgrd(477)=   0.951229691505432d0
      basgrd(478)=  -0.107761196792126d0
      basgrd(479)=  -0.346572488546372d0
      basgrd(480)=   0.931812822818756d0
      basgrd(481)=  -0.284454613924027d0
      basgrd(482)=  -0.380631893873215d0
      basgrd(483)=   0.879889190196991d0
      basgrd(484)=  -0.448512911796570d0
      basgrd(485)=  -0.400829493999481d0
      basgrd(486)=   0.798856616020203d0
      basgrd(487)=  -0.996499896049500d0
      basgrd(488)=  -0.017182800918818d0
      basgrd(489)=   0.081809103488922d0
      basgrd(490)=  -0.982108294963837d0
      basgrd(491)=  -0.187839493155479d0
      basgrd(492)=   0.013403000310063d0
      basgrd(493)=  -0.931238591670990d0
      basgrd(494)=  -0.361656397581101d0
      basgrd(495)=  -0.044714599847794d0
      basgrd(496)=  -0.848621785640717d0
      basgrd(497)=  -0.521448493003845d0
      basgrd(498)=  -0.089064799249172d0
      basgrd(499)=  -0.966082215309143d0
      basgrd(500)=  -0.137317895889282d0
      basgrd(501)=   0.218698292970657d0
      basgrd(502)=  -0.932948827743530d0
      basgrd(503)=  -0.325651288032532d0
      basgrd(504)=   0.153485700488091d0
      basgrd(505)=  -0.858389079570770d0
      basgrd(506)=  -0.504185199737549d0
      basgrd(507)=   0.094686202704906d0
      basgrd(508)=  -0.901239991188049d0
      basgrd(509)=  -0.265011191368103d0
      basgrd(510)=   0.342834502458572d0
      basgrd(511)=  -0.843641877174377d0
      basgrd(512)=  -0.458153098821640d0
      basgrd(513)=   0.279935985803604d0
      basgrd(514)=  -0.807126522064209d0
      basgrd(515)=  -0.388227313756943d0
      basgrd(516)=   0.444776713848114d0
      basgrd(517)=  -0.955810427665710d0
      basgrd(518)=   0.114013500511646d0
      basgrd(519)=  -0.270974814891815d0
      basgrd(520)=  -0.890308320522308d0
      basgrd(521)=   0.108410701155663d0
      basgrd(522)=  -0.442264914512634d0
      basgrd(523)=  -0.791631817817688d0
      basgrd(524)=   0.088033303618431d0
      basgrd(525)=  -0.604623317718506d0
      basgrd(526)=  -0.671117782592773d0
      basgrd(527)=   0.053371198475361d0
      basgrd(528)=  -0.739426970481873d0
      basgrd(529)=  -0.942444026470184d0
      basgrd(530)=  -0.061608098447323d0
      basgrd(531)=  -0.328639388084412d0
      basgrd(532)=  -0.856761813163757d0
      basgrd(533)=  -0.078578196465969d0
      basgrd(534)=  -0.509690701961517d0
      basgrd(535)=  -0.734996080398560d0
      basgrd(536)=  -0.105050101876259d0
      basgrd(537)=  -0.669884502887726d0
      basgrd(538)=  -0.891810297966003d0
      basgrd(539)=  -0.235004097223282d0
      basgrd(540)=  -0.386584401130676d0
      basgrd(541)=  -0.782294809818268d0
      basgrd(542)=  -0.259133189916611d0
      basgrd(543)=  -0.566449284553528d0
      basgrd(544)=  -0.807900309562683d0
      basgrd(545)=  -0.390404999256134d0
      basgrd(546)=  -0.441453307867050d0
      basgrd(547)=  -0.843478322029114d0
      basgrd(548)=   0.475457787513733d0
      basgrd(549)=  -0.249968603253365d0
      basgrd(550)=  -0.732556402683258d0
      basgrd(551)=   0.616495370864868d0
      basgrd(552)=  -0.288608014583588d0
      basgrd(553)=  -0.590696394443512d0
      basgrd(554)=   0.736588001251221d0
      basgrd(555)=  -0.329417586326599d0
      basgrd(556)=  -0.433169901371002d0
      basgrd(557)=   0.822478473186493d0
      basgrd(558)=  -0.368636608123779d0
      basgrd(559)=  -0.780652403831482d0
      basgrd(560)=   0.459831714630127d0
      basgrd(561)=  -0.423245400190353d0
      basgrd(562)=  -0.646050691604614d0
      basgrd(563)=   0.601377904415131d0
      basgrd(564)=  -0.470067113637924d0
      basgrd(565)=  -0.481487393379211d0
      basgrd(566)=   0.713278591632843d0
      basgrd(567)=  -0.509316802024841d0
      basgrd(568)=  -0.683246076107025d0
      basgrd(569)=   0.438351899385452d0
      basgrd(570)=  -0.583971202373505d0
      basgrd(571)=  -0.524859428405762d0
      basgrd(572)=   0.571962296962738d0
      basgrd(573)=  -0.630382180213928d0
      basgrd(574)=  -0.558949470520020d0
      basgrd(575)=   0.414406597614288d0
      basgrd(576)=  -0.718221783638001d0
      basgrd(577)=  -0.814821898937225d0
      basgrd(578)=   0.567725419998169d0
      basgrd(579)=   0.117273896932602d0
      basgrd(580)=  -0.727530479431152d0
      basgrd(581)=   0.632435500621796d0
      basgrd(582)=   0.265941292047501d0
      basgrd(583)=  -0.606588125228882d0
      basgrd(584)=   0.685797274112701d0
      basgrd(585)=   0.402160197496414d0
      basgrd(586)=  -0.464116811752319d0
      basgrd(587)=   0.722723722457886d0
      basgrd(588)=   0.512119114398956d0
      basgrd(589)=  -0.703997492790222d0
      basgrd(590)=   0.706884980201721d0
      basgrd(591)=   0.068565301597118d0
      basgrd(592)=  -0.592185199260712d0
      basgrd(593)=   0.775003314018250d0
      basgrd(594)=   0.220650404691696d0
      basgrd(595)=  -0.448390603065491d0
      basgrd(596)=   0.819903314113617d0
      basgrd(597)=   0.355955809354782d0
      basgrd(598)=  -0.563202679157257d0
      basgrd(599)=   0.825883209705353d0
      basgrd(600)=   0.026828600093722d0
      basgrd(601)=  -0.427419006824493d0
      basgrd(602)=   0.886298775672913d0
      basgrd(603)=   0.178290203213692d0
      basgrd(604)=  -0.404313594102860d0
      basgrd(605)=   0.914618492126465d0
      basgrd(606)=  -0.001894599990919d0
      basgrd(607)=  -0.909571111202240d0
      basgrd(608)=   0.262643486261368d0
      basgrd(609)=   0.322022914886475d0
      basgrd(610)=  -0.882299304008484d0
      basgrd(611)=   0.133743494749069d0
      basgrd(612)=   0.451287806034088d0
      basgrd(613)=  -0.817346215248108d0
      basgrd(614)=   0.006322099827230d0
      basgrd(615)=   0.576112091541290d0
      basgrd(616)=  -0.721229076385498d0
      basgrd(617)=  -0.108218401670456d0
      basgrd(618)=   0.684190988540649d0
      basgrd(619)=  -0.819292724132538d0
      basgrd(620)=   0.335824310779572d0
      basgrd(621)=   0.464738190174103d0
      basgrd(622)=  -0.770075500011444d0
      basgrd(623)=   0.200730696320534d0
      basgrd(624)=   0.605550110340118d0
      basgrd(625)=  -0.681682407855988d0
      basgrd(626)=   0.066955797374249d0
      basgrd(627)=   0.728578090667725d0
      basgrd(628)=  -0.698426425457001d0
      basgrd(629)=   0.389460802078247d0
      basgrd(630)=   0.600433886051178d0
      basgrd(631)=  -0.624847292900085d0
      basgrd(632)=   0.248850598931313d0
      basgrd(633)=   0.740026473999023d0
      basgrd(634)=  -0.557699322700501d0
      basgrd(635)=   0.417982906103134d0
      basgrd(636)=   0.717120528221130d0
      basgrd(637)=   0.955843806266785d0
      basgrd(638)=  -0.114249803125858d0
      basgrd(639)=   0.270757406949997d0
      basgrd(640)=   0.942751884460449d0
      basgrd(641)=   0.061342500150204d0
      basgrd(642)=   0.327804893255234d0
      basgrd(643)=   0.891754209995270d0
      basgrd(644)=   0.234839394688606d0
      basgrd(645)=   0.386813789606094d0
      basgrd(646)=   0.807701170444489d0
      basgrd(647)=   0.390947788953781d0
      basgrd(648)=   0.441337287425995d0
      basgrd(649)=   0.889850914478302d0
      basgrd(650)=  -0.108384199440479d0
      basgrd(651)=   0.443190902471542d0
      basgrd(652)=   0.856887578964233d0
      basgrd(653)=   0.078702099621296d0
      basgrd(654)=   0.509460091590881d0
      basgrd(655)=   0.782245218753815d0
      basgrd(656)=   0.259107112884522d0
      basgrd(657)=   0.566529691219330d0
      basgrd(658)=   0.792142927646637d0
      basgrd(659)=  -0.087629102170467d0
      basgrd(660)=   0.604012310504913d0
      basgrd(661)=   0.734840512275696d0
      basgrd(662)=   0.104854799807072d0
      basgrd(663)=   0.670085728168488d0
      basgrd(664)=   0.671307921409607d0
      basgrd(665)=  -0.053294401615858d0
      basgrd(666)=   0.739260017871857d0
      basgrd(667)=   0.996546506881714d0
      basgrd(668)=   0.016390899196267d0
      basgrd(669)=  -0.081402502954006d0
      basgrd(670)=   0.966223120689392d0
      basgrd(671)=   0.137170195579529d0
      basgrd(672)=  -0.218167901039124d0
      basgrd(673)=   0.900963604450226d0
      basgrd(674)=   0.265118986368179d0
      basgrd(675)=  -0.343477100133896d0
      basgrd(676)=   0.806736171245575d0
      basgrd(677)=   0.388470113277435d0
      basgrd(678)=  -0.445272505283356d0
      basgrd(679)=   0.981971085071564d0
      basgrd(680)=   0.188529804348946d0
      basgrd(681)=  -0.013760999776423d0
      basgrd(682)=   0.933140575885773d0
      basgrd(683)=   0.324862003326416d0
      basgrd(684)=  -0.153991296887398d0
      basgrd(685)=   0.843617916107178d0
      basgrd(686)=   0.458396911621094d0
      basgrd(687)=  -0.279608994722366d0
      basgrd(688)=   0.931219220161438d0
      basgrd(689)=   0.361704289913178d0
      basgrd(690)=   0.044730801135302d0
      basgrd(691)=   0.858283519744873d0
      basgrd(692)=   0.504396796226502d0
      basgrd(693)=  -0.094516001641750d0
      basgrd(694)=   0.849210917949677d0
      basgrd(695)=   0.520598173141480d0
      basgrd(696)=   0.088422000408173d0
      basgrd(697)=   0.909677684307098d0
      basgrd(698)=  -0.262331008911133d0
      basgrd(699)=  -0.321976512670517d0
      basgrd(700)=   0.819799184799194d0
      basgrd(701)=  -0.335445702075958d0
      basgrd(702)=  -0.464117914438248d0
      basgrd(703)=   0.697874188423157d0
      basgrd(704)=  -0.389874786138535d0
      basgrd(705)=  -0.600807070732117d0
      basgrd(706)=   0.557675600051880d0
      basgrd(707)=  -0.418078988790512d0
      basgrd(708)=  -0.717082917690277d0
      basgrd(709)=   0.881988823413849d0
      basgrd(710)=  -0.133651703596115d0
      basgrd(711)=  -0.451921403408051d0
      basgrd(712)=   0.770100176334381d0
      basgrd(713)=  -0.200819104909897d0
      basgrd(714)=  -0.605489373207092d0
      basgrd(715)=   0.624965727329254d0
      basgrd(716)=  -0.248694002628326d0
      basgrd(717)=  -0.739979088306427d0
      basgrd(718)=   0.817378818988800d0
      basgrd(719)=  -0.006242999807000d0
      basgrd(720)=  -0.576066792011261d0
      basgrd(721)=   0.681771397590637d0
      basgrd(722)=  -0.066777601838112d0
      basgrd(723)=  -0.728511095046997d0
      basgrd(724)=   0.721353590488434d0
      basgrd(725)=   0.107727102935314d0
      basgrd(726)=  -0.684137284755707d0
      basgrd(727)=   0.814422786235809d0
      basgrd(728)=  -0.568202316761017d0
      basgrd(729)=  -0.117735400795937d0
      basgrd(730)=   0.704616487026215d0
      basgrd(731)=  -0.706302583217621d0
      basgrd(732)=  -0.068206898868084d0
      basgrd(733)=   0.562912106513977d0
      basgrd(734)=  -0.826085627079010d0
      basgrd(735)=  -0.026694100350142d0
      basgrd(736)=   0.404308885335922d0
      basgrd(737)=  -0.914620399475098d0
      basgrd(738)=   0.001960299909115d0
      basgrd(739)=   0.727052628993988d0
      basgrd(740)=  -0.632898509502411d0
      basgrd(741)=  -0.266146510839462d0
      basgrd(742)=   0.591887474060059d0
      basgrd(743)=  -0.775162279605866d0
      basgrd(744)=  -0.220890298485756d0
      basgrd(745)=   0.427700310945511d0
      basgrd(746)=  -0.886152982711792d0
      basgrd(747)=  -0.178340598940849d0
      basgrd(748)=   0.606748998165131d0
      basgrd(749)=  -0.685755908489227d0
      basgrd(750)=  -0.401988089084625d0
      basgrd(751)=   0.448320508003235d0
      basgrd(752)=  -0.819860577583313d0
      basgrd(753)=  -0.356142401695252d0
      basgrd(754)=   0.464112699031830d0
      basgrd(755)=  -0.722764670848846d0
      basgrd(756)=  -0.512065112590790d0
      basgrd(757)=   0.843236684799194d0
      basgrd(758)=  -0.475735902786255d0
      basgrd(759)=   0.250254303216934d0
      basgrd(760)=   0.781166493892670d0
      basgrd(761)=  -0.459371507167816d0
      basgrd(762)=   0.422796308994293d0
      basgrd(763)=   0.683324813842773d0
      basgrd(764)=  -0.438166499137878d0
      basgrd(765)=   0.584018170833588d0
      basgrd(766)=   0.558947324752808d0
      basgrd(767)=  -0.414234101772308d0
      basgrd(768)=   0.718323111534119d0
      basgrd(769)=   0.732087492942810d0
      basgrd(770)=  -0.616853296756744d0
      basgrd(771)=   0.289032608270645d0
      basgrd(772)=   0.645748794078827d0
      basgrd(773)=  -0.601580977439880d0
      basgrd(774)=   0.470221996307373d0
      basgrd(775)=   0.524907886981964d0
      basgrd(776)=  -0.572089076042175d0
      basgrd(777)=   0.630226790904999d0
      basgrd(778)=   0.591070890426636d0
      basgrd(779)=  -0.736291706562042d0
      basgrd(780)=   0.329408288002014d0
      basgrd(781)=   0.481404602527618d0
      basgrd(782)=  -0.713346004486084d0
      basgrd(783)=   0.509300470352173d0
      basgrd(784)=   0.433061510324478d0
      basgrd(785)=  -0.822513103485107d0
      basgrd(786)=   0.368686795234680d0
      basgrd(787)=  -0.555007815361023d0
      basgrd(788)=  -0.668024301528931d0
      basgrd(789)=   0.495691299438477d0
      basgrd(790)=  -0.429905414581299d0
      basgrd(791)=  -0.789535105228424d0
      basgrd(792)=   0.437967687845230d0
      basgrd(793)=  -0.291748404502869d0
      basgrd(794)=  -0.888424694538117d0
      basgrd(795)=   0.354378908872604d0
      basgrd(796)=  -0.155646994709969d0
      basgrd(797)=  -0.955260694026947d0
      basgrd(798)=   0.251497507095337d0
      basgrd(799)=  -0.600628972053528d0
      basgrd(800)=  -0.729225099086762d0
      basgrd(801)=   0.327834606170654d0
      basgrd(802)=  -0.461556494235992d0
      basgrd(803)=  -0.851542592048645d0
      basgrd(804)=   0.248678103089333d0
      basgrd(805)=  -0.307511508464813d0
      basgrd(806)=  -0.939728319644928d0
      basgrd(807)=   0.149490103125572d0
      basgrd(808)=  -0.616222202777863d0
      basgrd(809)=  -0.774403512477875d0
      basgrd(810)=   0.143420398235321d0
      basgrd(811)=  -0.461458414793015d0
      basgrd(812)=  -0.885725080966950d0
      basgrd(813)=   0.050468798726797d0
      basgrd(814)=  -0.598678588867188d0
      basgrd(815)=  -0.800115585327148d0
      basgrd(816)=  -0.037402201443911d0
      basgrd(817)=  -0.529941499233246d0
      basgrd(818)=  -0.588316917419434d0
      basgrd(819)=  -0.610774278640747d0
      basgrd(820)=  -0.372514009475708d0
      basgrd(821)=  -0.605952799320221d0
      basgrd(822)=  -0.702890098094940d0
      basgrd(823)=  -0.200168296694756d0
      basgrd(824)=  -0.593625128269196d0
      basgrd(825)=  -0.779449820518494d0
      basgrd(826)=  -0.030707199126482d0
      basgrd(827)=  -0.551283597946167d0
      basgrd(828)=  -0.833752691745758d0
      basgrd(829)=  -0.513753890991211d0
      basgrd(830)=  -0.449691891670227d0
      basgrd(831)=  -0.730639517307282d0
      basgrd(832)=  -0.336925506591797d0
      basgrd(833)=  -0.450469493865967d0
      basgrd(834)=  -0.826775908470154d0
      basgrd(835)=  -0.147565200924873d0
      basgrd(836)=  -0.422862291336060d0
      basgrd(837)=  -0.894098401069641d0
      basgrd(838)=  -0.467891693115234d0
      basgrd(839)=  -0.296358585357666d0
      basgrd(840)=  -0.832615792751312d0
      basgrd(841)=  -0.274746209383011d0
      basgrd(842)=  -0.282504886388779d0
      basgrd(843)=  -0.919078588485718d0
      basgrd(844)=  -0.395839691162109d0
      basgrd(845)=  -0.144624501466751d0
      basgrd(846)=  -0.906859815120697d0
      basgrd(847)=  -0.212325707077980d0
      basgrd(848)=   0.440924793481827d0
      basgrd(849)=  -0.872068285942078d0
      basgrd(850)=  -0.034668900072575d0
      basgrd(851)=   0.487085014581680d0
      basgrd(852)=  -0.872666180133820d0
      basgrd(853)=   0.147117599844933d0
      basgrd(854)=   0.527473092079163d0
      basgrd(855)=  -0.836736798286438d0
      basgrd(856)=   0.313184589147568d0
      basgrd(857)=   0.559683382511139d0
      basgrd(858)=  -0.767248272895813d0
      basgrd(859)=  -0.189149007201195d0
      basgrd(860)=   0.597847402095795d0
      basgrd(861)=  -0.778974413871765d0
      basgrd(862)=   0.004314499907196d0
      basgrd(863)=   0.650188207626343d0
      basgrd(864)=  -0.759760916233063d0
      basgrd(865)=   0.195087596774101d0
      basgrd(866)=   0.684752583503723d0
      basgrd(867)=  -0.702178478240967d0
      basgrd(868)=  -0.146802902221680d0
      basgrd(869)=   0.738281011581421d0
      basgrd(870)=  -0.658323705196381d0
      basgrd(871)=   0.055920299142599d0
      basgrd(872)=   0.784260571002960d0
      basgrd(873)=  -0.617906391620636d0
      basgrd(874)=  -0.088265202939510d0
      basgrd(875)=   0.847823917865753d0
      basgrd(876)=  -0.522880375385284d0
      basgrd(877)=  -0.041256900876760d0
      basgrd(878)=   0.996635079383850d0
      basgrd(879)=   0.070825800299644d0
      basgrd(880)=   0.117095701396465d0
      basgrd(881)=   0.979720473289490d0
      basgrd(882)=   0.162592306733131d0
      basgrd(883)=   0.270790308713913d0
      basgrd(884)=   0.926381826400757d0
      basgrd(885)=   0.261704802513123d0
      basgrd(886)=   0.401551395654678d0
      basgrd(887)=   0.842503488063812d0
      basgrd(888)=   0.359088093042374d0
      basgrd(889)=  -0.075381897389889d0
      basgrd(890)=   0.965573191642761d0
      basgrd(891)=   0.248969793319702d0
      basgrd(892)=   0.091270402073860d0
      basgrd(893)=   0.929783403873444d0
      basgrd(894)=   0.356612592935562d0
      basgrd(895)=   0.247175306081772d0
      basgrd(896)=   0.853237092494965d0
      basgrd(897)=   0.459228396415710d0
      basgrd(898)=  -0.096757203340530d0
      basgrd(899)=   0.899540603160858d0
      basgrd(900)=   0.425986886024475d0
      basgrd(901)=   0.073409296572208d0
      basgrd(902)=   0.840517580509186d0
      basgrd(903)=   0.536787927150726d0
      basgrd(904)=  -0.101280301809311d0
      basgrd(905)=   0.805675089359283d0
      basgrd(906)=   0.583635091781616d0
      basgrd(907)=  -0.252299785614014d0
      basgrd(908)=   0.312042295932770d0
      basgrd(909)=   0.915955483913422d0
      basgrd(910)=  -0.126796200871468d0
      basgrd(911)=   0.190663799643517d0
      basgrd(912)=   0.973432123661041d0
      basgrd(913)=  -0.000521100009792d0
      basgrd(914)=   0.051185000687838d0
      basgrd(915)=   0.998689115047455d0
      basgrd(916)=   0.111481100320816d0
      basgrd(917)=  -0.093628197908401d0
      basgrd(918)=   0.989346086978912d0
      basgrd(919)=  -0.329458296298981d0
      basgrd(920)=   0.145694002509117d0
      basgrd(921)=   0.932861506938934d0
      basgrd(922)=  -0.196815103292465d0
      basgrd(923)=   0.001666700001806d0
      basgrd(924)=   0.980439186096191d0
      basgrd(925)=  -0.063606202602386d0
      basgrd(926)=  -0.150688305497170d0
      basgrd(927)=   0.986532986164093d0
      basgrd(928)=  -0.387266814708710d0
      basgrd(929)=  -0.035602398216724d0
      basgrd(930)=   0.921280026435852d0
      basgrd(931)=  -0.246388405561447d0
      basgrd(932)=  -0.191062405705452d0
      basgrd(933)=   0.950151503086090d0
      basgrd(934)=  -0.416676104068756d0
      basgrd(935)=  -0.212497398257256d0
      basgrd(936)=   0.883869886398315d0
      basgrd(937)=   0.041219998151064d0
      basgrd(938)=  -0.996641576290131d0
      basgrd(939)=  -0.070756800472736d0
      basgrd(940)=  -0.116947896778584d0
      basgrd(941)=  -0.979738593101502d0
      basgrd(942)=  -0.162589594721794d0
      basgrd(943)=  -0.270508885383606d0
      basgrd(944)=  -0.926447391510010d0
      basgrd(945)=  -0.261763513088226d0
      basgrd(946)=  -0.401443511247635d0
      basgrd(947)=  -0.842541277408600d0
      basgrd(948)=  -0.359120100736618d0
      basgrd(949)=   0.075765602290630d0
      basgrd(950)=  -0.965560376644135d0
      basgrd(951)=  -0.248903095722199d0
      basgrd(952)=  -0.091307103633881d0
      basgrd(953)=  -0.929808616638184d0
      basgrd(954)=  -0.356537491083145d0
      basgrd(955)=  -0.247341394424439d0
      basgrd(956)=  -0.853136897087097d0
      basgrd(957)=  -0.459325313568115d0
      basgrd(958)=   0.096485301852226d0
      basgrd(959)=  -0.899550080299377d0
      basgrd(960)=  -0.426028400659561d0
      basgrd(961)=  -0.073645897209644d0
      basgrd(962)=  -0.840444386005402d0
      basgrd(963)=  -0.536870121955872d0
      basgrd(964)=   0.101245798170567d0
      basgrd(965)=  -0.805717289447784d0
      basgrd(966)=  -0.583582818508148d0
      basgrd(967)=   0.252238988876343d0
      basgrd(968)=  -0.312179595232010d0
      basgrd(969)=  -0.915925383567810d0
      basgrd(970)=   0.127005904912949d0
      basgrd(971)=  -0.190736800432205d0
      basgrd(972)=  -0.973390400409699d0
      basgrd(973)=   0.000519799999893d0
      basgrd(974)=  -0.051125101745129d0
      basgrd(975)=  -0.998692095279694d0
      basgrd(976)=  -0.111413799226284d0
      basgrd(977)=   0.093655601143837d0
      basgrd(978)=  -0.989351093769074d0
      basgrd(979)=   0.329600393772125d0
      basgrd(980)=  -0.145705401897430d0
      basgrd(981)=  -0.932809472084045d0
      basgrd(982)=   0.196763902902603d0
      basgrd(983)=  -0.001796000055037d0
      basgrd(984)=  -0.980449318885803d0
      basgrd(985)=   0.063660003244877d0
      basgrd(986)=   0.150646492838860d0
      basgrd(987)=  -0.986535906791687d0
      basgrd(988)=   0.387062400579453d0
      basgrd(989)=   0.035640601068735d0
      basgrd(990)=  -0.921364426612854d0
      basgrd(991)=   0.246054396033287d0
      basgrd(992)=   0.191412493586540d0
      basgrd(993)=  -0.950167596340179d0
      basgrd(994)=   0.416655391454697d0
      basgrd(995)=   0.212424993515015d0
      basgrd(996)=  -0.883897006511688d0
      basgrd(997)=   0.554838776588440d0
      basgrd(998)=   0.668085396289825d0
      basgrd(999)=  -0.495798110961914d0
      basgrd(1000)=   0.430376797914505d0
      basgrd(1001)=   0.789387404918671d0
      basgrd(1002)=  -0.437770903110504d0
      basgrd(1003)=   0.291598290205002d0
      basgrd(1004)=   0.888499021530151d0
      basgrd(1005)=  -0.354316204786301d0
      basgrd(1006)=   0.155681505799294d0
      basgrd(1007)=   0.955271899700165d0
      basgrd(1008)=  -0.251433491706848d0
      basgrd(1009)=   0.600600004196167d0
      basgrd(1010)=   0.729302883148193d0
      basgrd(1011)=  -0.327714711427689d0
      basgrd(1012)=   0.461643397808075d0
      basgrd(1013)=   0.851453900337219d0
      basgrd(1014)=  -0.248820498585701d0
      basgrd(1015)=   0.307629793882370d0
      basgrd(1016)=   0.939690113067627d0
      basgrd(1017)=  -0.149487093091011d0
      basgrd(1018)=   0.616129279136658d0
      basgrd(1019)=   0.774490773677826d0
      basgrd(1020)=  -0.143348097801209d0
      basgrd(1021)=   0.461036205291748d0
      basgrd(1022)=   0.885959088802338d0
      basgrd(1023)=  -0.050221301615238d0
      basgrd(1024)=   0.598827302455902d0
      basgrd(1025)=   0.800007998943329d0
      basgrd(1026)=   0.037325199693441d0
      basgrd(1027)=   0.529993891716003d0
      basgrd(1028)=   0.588458120822907d0
      basgrd(1029)=   0.610592722892761d0
      basgrd(1030)=   0.372833490371704d0
      basgrd(1031)=   0.605608820915222d0
      basgrd(1032)=   0.703017115592957d0
      basgrd(1033)=   0.200061693787575d0
      basgrd(1034)=   0.593664824962616d0
      basgrd(1035)=   0.779446899890900d0
      basgrd(1036)=   0.030695199966431d0
      basgrd(1037)=   0.551244378089905d0
      basgrd(1038)=   0.833778977394104d0
      basgrd(1039)=   0.513854503631592d0
      basgrd(1040)=   0.449585795402527d0
      basgrd(1041)=   0.730634093284607d0
      basgrd(1042)=   0.337098389863968d0
      basgrd(1043)=   0.450643599033356d0
      basgrd(1044)=   0.826610624790192d0
      basgrd(1045)=   0.147569298744202d0
      basgrd(1046)=   0.422796607017517d0
      basgrd(1047)=   0.894128799438477d0
      basgrd(1048)=   0.467734307050705d0
      basgrd(1049)=   0.296364486217499d0
      basgrd(1050)=   0.832702100276947d0
      basgrd(1051)=   0.274422198534012d0
      basgrd(1052)=   0.282555907964706d0
      basgrd(1053)=   0.919159710407257d0
      basgrd(1054)=   0.395973801612854d0
      basgrd(1055)=   0.144747793674469d0
      basgrd(1056)=   0.906781613826752d0
      basgrd(1057)=   0.212354198098183d0
      basgrd(1058)=  -0.440843105316162d0
      basgrd(1059)=   0.872102677822113d0
      basgrd(1060)=   0.035041298717260d0
      basgrd(1061)=  -0.487280398607254d0
      basgrd(1062)=   0.872542202472687d0
      basgrd(1063)=  -0.147180899977684d0
      basgrd(1064)=  -0.527492523193359d0
      basgrd(1065)=   0.836713492870331d0
      basgrd(1066)=  -0.313206702470779d0
      basgrd(1067)=  -0.559703171253204d0
      basgrd(1068)=   0.767224788665772d0
      basgrd(1069)=   0.189358502626419d0
      basgrd(1070)=  -0.597856879234314d0
      basgrd(1071)=   0.778916180133820d0
      basgrd(1072)=  -0.004297200124711d0
      basgrd(1073)=  -0.650083899497986d0
      basgrd(1074)=   0.759850323200226d0
      basgrd(1075)=  -0.195179194211960d0
      basgrd(1076)=  -0.684763312339783d0
      basgrd(1077)=   0.702142596244812d0
      basgrd(1078)=   0.146599605679512d0
      basgrd(1079)=  -0.738271892070770d0
      basgrd(1080)=   0.658379077911377d0
      basgrd(1081)=  -0.056144498288631d0
      basgrd(1082)=  -0.784226119518280d0
      basgrd(1083)=   0.617929697036743d0
      basgrd(1084)=   0.088267602026463d0
      basgrd(1085)=  -0.847751379013062d0
      basgrd(1086)=   0.522997498512268d0
      basgrd(1087)=  -0.981259107589722d0
      basgrd(1088)=   0.190601006150246d0
      basgrd(1089)=   0.028315400704741d0
      basgrd(1090)=  -0.998927712440491d0
      basgrd(1091)=   0.025483399629593d0
      basgrd(1092)=  -0.038653701543808d0
      basgrd(1093)=  -0.983459115028381d0
      basgrd(1094)=  -0.148843899369240d0
      basgrd(1095)=  -0.103217199444771d0
      basgrd(1096)=  -0.932625293731690d0
      basgrd(1097)=  -0.322989106178284d0
      basgrd(1098)=  -0.160897493362427d0
      basgrd(1099)=  -0.850076913833618d0
      basgrd(1100)=  -0.482813596725464d0
      basgrd(1101)=  -0.210381299257278d0
      basgrd(1102)=  -0.742531001567841d0
      basgrd(1103)=  -0.622283518314362d0
      basgrd(1104)=  -0.247812300920487d0
      basgrd(1105)=  -0.982966125011444d0
      basgrd(1106)=   0.076667398214340d0
      basgrd(1107)=   0.167032092809677d0
      basgrd(1108)=  -0.986813187599182d0
      basgrd(1109)=  -0.121288396418095d0
      basgrd(1110)=   0.107185602188110d0
      basgrd(1111)=  -0.955336391925812d0
      basgrd(1112)=  -0.293090194463730d0
      basgrd(1113)=   0.037822101265192d0
      basgrd(1114)=  -0.887370109558106d0
      basgrd(1115)=  -0.460985094308853d0
      basgrd(1116)=  -0.008189699612558d0
      basgrd(1117)=  -0.777231514453888d0
      basgrd(1118)=  -0.625240921974182d0
      basgrd(1119)=  -0.070604600012302d0
      basgrd(1120)=  -0.952023923397064d0
      basgrd(1121)=  -0.048403598368168d0
      basgrd(1122)=   0.302171498537064d0
      basgrd(1123)=  -0.939646124839783d0
      basgrd(1124)=  -0.242044493556023d0
      basgrd(1125)=   0.241825506091118d0
      basgrd(1126)=  -0.882611513137817d0
      basgrd(1127)=  -0.435110509395599d0
      basgrd(1128)=   0.177977204322815d0
      basgrd(1129)=  -0.789470791816711d0
      basgrd(1130)=  -0.602904081344605d0
      basgrd(1131)=   0.115076102316380d0
      basgrd(1132)=  -0.886988699436188d0
      basgrd(1133)=  -0.175924196839333d0
      basgrd(1134)=   0.426968097686768d0
      basgrd(1135)=  -0.859096825122833d0
      basgrd(1136)=  -0.371376693248749d0
      basgrd(1137)=   0.352181911468506d0
      basgrd(1138)=  -0.774177312850952d0
      basgrd(1139)=  -0.557899713516235d0
      basgrd(1140)=   0.298994213342667d0
      basgrd(1141)=  -0.793080806732178d0
      basgrd(1142)=  -0.296995699405670d0
      basgrd(1143)=   0.531804919242859d0
      basgrd(1144)=  -0.734284222126007d0
      basgrd(1145)=  -0.491394788026810d0
      basgrd(1146)=   0.468356490135193d0
      basgrd(1147)=  -0.675649225711823d0
      basgrd(1148)=  -0.407306909561157d0
      basgrd(1149)=   0.614490985870361d0
      basgrd(1150)=  -0.966679811477661d0
      basgrd(1151)=   0.237567096948624d0
      basgrd(1152)=  -0.095352202653885d0
      basgrd(1153)=  -0.932256698608398d0
      basgrd(1154)=   0.239794105291367d0
      basgrd(1155)=  -0.270917505025864d0
      basgrd(1156)=  -0.865890383720398d0
      basgrd(1157)=   0.229634106159210d0
      basgrd(1158)=  -0.444411903619766d0
      basgrd(1159)=  -0.768491923809052d0
      basgrd(1160)=   0.208023399114609d0
      basgrd(1161)=  -0.605100393295288d0
      basgrd(1162)=  -0.645207405090332d0
      basgrd(1163)=   0.178497701883316d0
      basgrd(1164)=  -0.742863416671753d0
      basgrd(1165)=  -0.507140696048737d0
      basgrd(1166)=   0.139699503779411d0
      basgrd(1167)=  -0.850466012954712d0
      basgrd(1168)=  -0.985367596149445d0
      basgrd(1169)=   0.068667098879814d0
      basgrd(1170)=  -0.155998393893242d0
      basgrd(1171)=  -0.934020221233368d0
      basgrd(1172)=   0.049693200737238d0
      basgrd(1173)=  -0.353746801614761d0
      basgrd(1174)=  -0.851090073585510d0
      basgrd(1175)=   0.042921699583530d0
      basgrd(1176)=  -0.523262202739716d0
      basgrd(1177)=  -0.741771876811981d0
      basgrd(1178)=   0.009559200145304d0
      basgrd(1179)=  -0.670584082603455d0
      basgrd(1180)=  -0.588742911815643d0
      basgrd(1181)=  -0.015250800177455d0
      basgrd(1182)=  -0.808176517486572d0
      basgrd(1183)=  -0.970162212848663d0
      basgrd(1184)=  -0.106755897402763d0
      basgrd(1185)=  -0.217688798904419d0
      basgrd(1186)=  -0.903855919837952d0
      basgrd(1187)=  -0.127434194087982d0
      basgrd(1188)=  -0.408417791128159d0
      basgrd(1189)=  -0.794965624809265d0
      basgrd(1190)=  -0.150444105267525d0
      basgrd(1191)=  -0.587704181671143d0
      basgrd(1192)=  -0.656565427780151d0
      basgrd(1193)=  -0.173520699143410d0
      basgrd(1194)=  -0.734038472175598d0
      basgrd(1195)=  -0.919372677803040d0
      basgrd(1196)=  -0.280559599399567d0
      basgrd(1197)=  -0.275753706693649d0
      basgrd(1198)=  -0.834168374538422d0
      basgrd(1199)=  -0.290042787790299d0
      basgrd(1200)=  -0.469082295894623d0
      basgrd(1201)=  -0.702866971492767d0
      basgrd(1202)=  -0.326373994350433d0
      basgrd(1203)=  -0.632026910781860d0
      basgrd(1204)=  -0.836579620838165d0
      basgrd(1205)=  -0.439239799976349d0
      basgrd(1206)=  -0.327418595552445d0
      basgrd(1207)=  -0.725766181945801d0
      basgrd(1208)=  -0.463410794734955d0
      basgrd(1209)=  -0.508442580699921d0
      basgrd(1210)=  -0.728081703186035d0
      basgrd(1211)=  -0.576516389846802d0
      basgrd(1212)=  -0.370844811201096d0
      basgrd(1213)=  -0.927329599857330d0
      basgrd(1214)=   0.364003390073776d0
      basgrd(1215)=  -0.086956001818180d0
      basgrd(1216)=  -0.845884382724762d0
      basgrd(1217)=   0.517498373985291d0
      basgrd(1218)=  -0.129132106900215d0
      basgrd(1219)=  -0.732806622982025d0
      basgrd(1220)=   0.658319115638733d0
      basgrd(1221)=  -0.172076895833015d0
      basgrd(1222)=  -0.592174828052521d0
      basgrd(1223)=   0.777040302753449d0
      basgrd(1224)=  -0.213394895195961d0
      basgrd(1225)=  -0.432546913623810d0
      basgrd(1226)=   0.866684317588806d0
      basgrd(1227)=  -0.248518496751785d0
      basgrd(1228)=  -0.264987885951996d0
      basgrd(1229)=   0.923009693622589d0
      basgrd(1230)=  -0.278988599777222d0
      basgrd(1231)=  -0.894683122634888d0
      basgrd(1232)=   0.359988093376160d0
      basgrd(1233)=  -0.264482110738754d0
      basgrd(1234)=  -0.787156224250794d0
      basgrd(1235)=   0.522795975208283d0
      basgrd(1236)=  -0.327214688062668d0
      basgrd(1237)=  -0.660561680793762d0
      basgrd(1238)=   0.657574474811554d0
      basgrd(1239)=  -0.362289786338806d0
      basgrd(1240)=  -0.509920597076416d0
      basgrd(1241)=   0.757714271545410d0
      basgrd(1242)=  -0.407246887683868d0
      basgrd(1243)=  -0.323142707347870d0
      basgrd(1244)=   0.843288779258728d0
      basgrd(1245)=  -0.429468095302582d0
      basgrd(1246)=  -0.829277217388153d0
      basgrd(1247)=   0.347157388925552d0
      basgrd(1248)=  -0.437928199768066d0
      basgrd(1249)=  -0.709179699420929d0
      basgrd(1250)=   0.500505387783051d0
      basgrd(1251)=  -0.496546506881714d0
      basgrd(1252)=  -0.552132129669190d0
      basgrd(1253)=   0.633674919605255d0
      basgrd(1254)=  -0.541854381561279d0
      basgrd(1255)=  -0.375991404056549d0
      basgrd(1256)=   0.731177687644959d0
      basgrd(1257)=  -0.569218397140503d0
      basgrd(1258)=  -0.731624305248261d0
      basgrd(1259)=   0.325923502445221d0
      basgrd(1260)=  -0.598748505115509d0
      basgrd(1261)=  -0.595075309276581d0
      basgrd(1262)=   0.481893301010132d0
      basgrd(1263)=  -0.643167316913605d0
      basgrd(1264)=  -0.418605685234070d0
      basgrd(1265)=   0.591195404529572d0
      basgrd(1266)=  -0.689388990402222d0
      basgrd(1267)=  -0.608789622783661d0
      basgrd(1268)=   0.298093795776367d0
      basgrd(1269)=  -0.735197484493256d0
      basgrd(1270)=  -0.449058592319489d0
      basgrd(1271)=   0.432573288679123d0
      basgrd(1272)=  -0.781809926033020d0
      basgrd(1273)=  -0.467279613018036d0
      basgrd(1274)=   0.265834897756577d0
      basgrd(1275)=  -0.843197286128998d0
      basgrd(1276)=  -0.917525708675385d0
      basgrd(1277)=   0.395533293485642d0
      basgrd(1278)=   0.041231200098991d0
      basgrd(1279)=  -0.859421312808991d0
      basgrd(1280)=   0.473967105150223d0
      basgrd(1281)=   0.191703394055367d0
      basgrd(1282)=  -0.768962681293488d0
      basgrd(1283)=   0.541904270648956d0
      basgrd(1284)=   0.339169710874558d0
      basgrd(1285)=  -0.648313820362091d0
      basgrd(1286)=   0.595432996749878d0
      basgrd(1287)=   0.474498510360718d0
      basgrd(1288)=  -0.505972981452942d0
      basgrd(1289)=   0.629720926284790d0
      basgrd(1290)=   0.589442908763886d0
      basgrd(1291)=  -0.350865900516510d0
      basgrd(1292)=   0.645922482013702d0
      basgrd(1293)=   0.677995026111603d0
      basgrd(1294)=  -0.835944175720215d0
      basgrd(1295)=   0.548770785331726d0
      basgrd(1296)=  -0.006925799883902d0
      basgrd(1297)=  -0.749559223651886d0
      basgrd(1298)=   0.643889904022217d0
      basgrd(1299)=   0.153514593839645d0
      basgrd(1300)=  -0.647323310375214d0
      basgrd(1301)=   0.700391113758087d0
      basgrd(1302)=   0.300707101821899d0
      basgrd(1303)=  -0.512743294239044d0
      basgrd(1304)=   0.748870611190796d0
      basgrd(1305)=   0.419865608215332d0
      basgrd(1306)=  -0.347534000873566d0
      basgrd(1307)=   0.764541387557983d0
      basgrd(1308)=   0.542859613895416d0
      basgrd(1309)=  -0.723388075828552d0
      basgrd(1310)=   0.688437879085541d0
      basgrd(1311)=  -0.052563499659300d0
      basgrd(1312)=  -0.624182283878326d0
      basgrd(1313)=   0.774451196193695d0
      basgrd(1314)=   0.103062197566032d0
      basgrd(1315)=  -0.490021109580994d0
      basgrd(1316)=   0.833849370479584d0
      basgrd(1317)=   0.254115194082260d0
      basgrd(1318)=  -0.335918009281158d0
      basgrd(1319)=   0.860402703285217d0
      basgrd(1320)=   0.383231490850449d0
      basgrd(1321)=  -0.582644999027252d0
      basgrd(1322)=   0.807323396205902d0
      basgrd(1323)=  -0.093561403453350d0
      basgrd(1324)=  -0.472317188978195d0
      basgrd(1325)=   0.878419101238251d0
      basgrd(1326)=   0.072776600718498d0
      basgrd(1327)=  -0.314614593982697d0
      basgrd(1328)=   0.926390171051025d0
      basgrd(1329)=   0.206926897168160d0
      basgrd(1330)=  -0.423543214797974d0
      basgrd(1331)=   0.896930813789368d0
      basgrd(1332)=  -0.126989498734474d0
      basgrd(1333)=  -0.287032812833786d0
      basgrd(1334)=   0.957584619522095d0
      basgrd(1335)=   0.025376100093126d0
      basgrd(1336)=  -0.254495203495026d0
      basgrd(1337)=   0.955222010612488d0
      basgrd(1338)=  -0.150940701365471d0
      basgrd(1339)=  -0.950729012489319d0
      basgrd(1340)=   0.288791000843048d0
      basgrd(1341)=   0.112757302820683d0
      basgrd(1342)=  -0.954073071479797d0
      basgrd(1343)=   0.169550105929375d0
      basgrd(1344)=   0.246976196765900d0
      basgrd(1345)=  -0.924109220504761d0
      basgrd(1346)=   0.042504198849201d0
      basgrd(1347)=   0.379757195711136d0
      basgrd(1348)=  -0.858936429023743d0
      basgrd(1349)=  -0.084622502326965d0
      basgrd(1350)=   0.505041897296906d0
      basgrd(1351)=  -0.764034986495972d0
      basgrd(1352)=  -0.204322293400765d0
      basgrd(1353)=   0.611966371536255d0
      basgrd(1354)=  -0.646422326564789d0
      basgrd(1355)=  -0.309472411870956d0
      basgrd(1356)=   0.697398722171783d0
      basgrd(1357)=  -0.890818417072296d0
      basgrd(1358)=   0.372620314359665d0
      basgrd(1359)=   0.259993493556976d0
      basgrd(1360)=  -0.873500704765320d0
      basgrd(1361)=   0.244410097599030d0
      basgrd(1362)=   0.421022891998291d0
      basgrd(1363)=  -0.829890370368958d0
      basgrd(1364)=   0.111665599048138d0
      basgrd(1365)=   0.546637594699860d0
      basgrd(1366)=  -0.746396005153656d0
      basgrd(1367)=  -0.005028599873185d0
      basgrd(1368)=   0.665483117103577d0
      basgrd(1369)=  -0.628467619419098d0
      basgrd(1370)=  -0.143092006444931d0
      basgrd(1371)=   0.764560818672180d0
      basgrd(1372)=  -0.799592971801758d0
      basgrd(1373)=   0.442698299884796d0
      basgrd(1374)=   0.405794590711594d0
      basgrd(1375)=  -0.767344117164612d0
      basgrd(1376)=   0.313330203294754d0
      basgrd(1377)=   0.559470415115356d0
      basgrd(1378)=  -0.694815695285797d0
      basgrd(1379)=   0.172552004456520d0
      basgrd(1380)=   0.698181211948395d0
      basgrd(1381)=  -0.591836929321289d0
      basgrd(1382)=   0.035454701632261d0
      basgrd(1383)=   0.805277585983276d0
      basgrd(1384)=  -0.679021477699280d0
      basgrd(1385)=   0.495864093303680d0
      basgrd(1386)=   0.541339695453644d0
      basgrd(1387)=  -0.635753810405731d0
      basgrd(1388)=   0.350063204765320d0
      basgrd(1389)=   0.687948286533356d0
      basgrd(1390)=  -0.534490525722504d0
      basgrd(1391)=   0.216015607118607d0
      basgrd(1392)=   0.817102909088135d0
      basgrd(1393)=  -0.537528693675995d0
      basgrd(1394)=   0.528457880020142d0
      basgrd(1395)=   0.657111227512360d0
      basgrd(1396)=  -0.463100492954254d0
      basgrd(1397)=   0.386735498905182d0
      basgrd(1398)=   0.797479510307312d0
      basgrd(1399)=  -0.382881492376328d0
      basgrd(1400)=   0.539724826812744d0
      basgrd(1401)=   0.749732494354248d0
      basgrd(1402)=   0.966864824295044d0
      basgrd(1403)=  -0.237155795097351d0
      basgrd(1404)=   0.094496399164200d0
      basgrd(1405)=   0.985437929630280d0
      basgrd(1406)=  -0.068713203072548d0
      basgrd(1407)=   0.155533492565155d0
      basgrd(1408)=   0.970363676548004d0
      basgrd(1409)=   0.106627300381661d0
      basgrd(1410)=   0.216852098703384d0
      basgrd(1411)=   0.919273495674133d0
      basgrd(1412)=   0.280133605003357d0
      basgrd(1413)=   0.276516407728195d0
      basgrd(1414)=   0.836413323879242d0
      basgrd(1415)=   0.439554214477539d0
      basgrd(1416)=   0.327421694993973d0
      basgrd(1417)=   0.727891385555267d0
      basgrd(1418)=   0.576797723770142d0
      basgrd(1419)=   0.370781004428864d0
      basgrd(1420)=   0.932314515113831d0
      basgrd(1421)=  -0.239295303821564d0
      basgrd(1422)=   0.271159499883652d0
      basgrd(1423)=   0.933808684349060d0
      basgrd(1424)=  -0.050440099090338d0
      basgrd(1425)=   0.354199290275574d0
      basgrd(1426)=   0.903921425342560d0
      basgrd(1427)=   0.127630695700645d0
      basgrd(1428)=   0.408211290836334d0
      basgrd(1429)=   0.834238529205322d0
      basgrd(1430)=   0.290072113275528d0
      basgrd(1431)=   0.468939602375031d0
      basgrd(1432)=   0.725421428680420d0
      basgrd(1433)=   0.463674992322922d0
      basgrd(1434)=   0.508693695068359d0
      basgrd(1435)=   0.865888118743897d0
      basgrd(1436)=  -0.229391306638718d0
      basgrd(1437)=   0.444541811943054d0
      basgrd(1438)=   0.851478815078735d0
      basgrd(1439)=  -0.042284000664949d0
      basgrd(1440)=   0.522681474685669d0
      basgrd(1441)=   0.794843673706055d0
      basgrd(1442)=   0.150282695889473d0
      basgrd(1443)=   0.587910413742065d0
      basgrd(1444)=   0.703099429607391d0
      basgrd(1445)=   0.326316714286804d0
      basgrd(1446)=   0.631798028945923d0
      basgrd(1447)=   0.768040180206299d0
      basgrd(1448)=  -0.208480894565582d0
      basgrd(1449)=   0.605516195297241d0
      basgrd(1450)=   0.741778492927551d0
      basgrd(1451)=  -0.009546199813485d0
      basgrd(1452)=   0.670576989650726d0
      basgrd(1453)=   0.656409084796906d0
      basgrd(1454)=   0.173556596040726d0
      basgrd(1455)=   0.734169781208038d0
      basgrd(1456)=   0.645877778530121d0
      basgrd(1457)=  -0.178160294890404d0
      basgrd(1458)=   0.742361605167389d0
      basgrd(1459)=   0.588729500770569d0
      basgrd(1460)=   0.015079099684954d0
      basgrd(1461)=   0.808189511299133d0
      basgrd(1462)=   0.507141292095184d0
      basgrd(1463)=  -0.139618098735809d0
      basgrd(1464)=   0.850479006767273d0
      basgrd(1465)=   0.980976104736328d0
      basgrd(1466)=  -0.191969498991966d0
      basgrd(1467)=  -0.028872299939394d0
      basgrd(1468)=   0.982953071594238d0
      basgrd(1469)=  -0.076644100248814d0
      basgrd(1470)=  -0.167119398713112d0
      basgrd(1471)=   0.952080070972443d0
      basgrd(1472)=   0.047909699380398d0
      basgrd(1473)=  -0.302073091268539d0
      basgrd(1474)=   0.887257516384125d0
      basgrd(1475)=   0.176348701119423d0
      basgrd(1476)=  -0.426233887672424d0
      basgrd(1477)=   0.792412698268890d0
      basgrd(1478)=   0.297270596027374d0
      basgrd(1479)=  -0.532646477222443d0
      basgrd(1480)=   0.675429821014404d0
      basgrd(1481)=   0.407395005226135d0
      basgrd(1482)=  -0.614673793315888d0
      basgrd(1483)=   0.998940289020538d0
      basgrd(1484)=  -0.025280199944973d0
      basgrd(1485)=   0.038461100310087d0
      basgrd(1486)=   0.987010896205902d0
      basgrd(1487)=   0.120705701410770d0
      basgrd(1488)=  -0.106017403304577d0
      basgrd(1489)=   0.939278721809387d0
      basgrd(1490)=   0.242372199892998d0
      basgrd(1491)=  -0.242922306060791d0
      basgrd(1492)=   0.859110891819000d0
      basgrd(1493)=   0.371154397726059d0
      basgrd(1494)=  -0.352381795644760d0
      basgrd(1495)=   0.734156906604767d0
      basgrd(1496)=   0.491739004850388d0
      basgrd(1497)=  -0.468194812536240d0
      basgrd(1498)=   0.983419775962830d0
      basgrd(1499)=   0.149559602141380d0
      basgrd(1500)=   0.102554000914097d0
      basgrd(1501)=   0.955367386341095d0
      basgrd(1502)=   0.292918205261231d0
      basgrd(1503)=  -0.038368199020624d0
      basgrd(1504)=   0.882893800735474d0
      basgrd(1505)=   0.434825211763382d0
      basgrd(1506)=  -0.177272498607636d0
      basgrd(1507)=   0.774103522300720d0
      basgrd(1508)=   0.557971477508545d0
      basgrd(1509)=  -0.299051314592362d0
      basgrd(1510)=   0.932579278945923d0
      basgrd(1511)=   0.322670400142670d0
      basgrd(1512)=   0.161801099777222d0
      basgrd(1513)=   0.887360811233521d0
      basgrd(1514)=   0.461014002561569d0
      basgrd(1515)=   0.007538499776274d0
      basgrd(1516)=   0.789314985275269d0
      basgrd(1517)=   0.603085696697235d0
      basgrd(1518)=  -0.115193299949169d0
      basgrd(1519)=   0.850335478782654d0
      basgrd(1520)=   0.482490807771683d0
      basgrd(1521)=   0.210076406598091d0
      basgrd(1522)=   0.777557492256165d0
      basgrd(1523)=   0.624804079532623d0
      basgrd(1524)=   0.070881798863411d0
      basgrd(1525)=   0.743068397045136d0
      basgrd(1526)=   0.621808588504791d0
      basgrd(1527)=   0.247393295168877d0
      basgrd(1528)=   0.950501322746277d0
      basgrd(1529)=  -0.289512902498245d0
      basgrd(1530)=  -0.112825803458691d0
      basgrd(1531)=   0.890996396541596d0
      basgrd(1532)=  -0.372450202703476d0
      basgrd(1533)=  -0.259627103805542d0
      basgrd(1534)=   0.799667418003082d0
      basgrd(1535)=  -0.442628502845764d0
      basgrd(1536)=  -0.405724108219147d0
      basgrd(1537)=   0.679468989372253d0
      basgrd(1538)=  -0.495324403047562d0
      basgrd(1539)=  -0.541272103786469d0
      basgrd(1540)=   0.536864817142487d0
      basgrd(1541)=  -0.528910517692566d0
      basgrd(1542)=  -0.657289803028107d0
      basgrd(1543)=   0.383011400699616d0
      basgrd(1544)=  -0.539755821228027d0
      basgrd(1545)=  -0.749643921852112d0
      basgrd(1546)=   0.954032897949219d0
      basgrd(1547)=  -0.169432997703552d0
      basgrd(1548)=  -0.247211694717407d0
      basgrd(1549)=   0.873804807662964d0
      basgrd(1550)=  -0.243584200739861d0
      basgrd(1551)=  -0.420870393514633d0
      basgrd(1552)=   0.766968429088593d0
      basgrd(1553)=  -0.313959389925003d0
      basgrd(1554)=  -0.559632897377014d0
      basgrd(1555)=   0.635765910148621d0
      basgrd(1556)=  -0.349997997283936d0
      basgrd(1557)=  -0.687970280647278d0
      basgrd(1558)=   0.463171303272247d0
      basgrd(1559)=  -0.386646598577499d0
      basgrd(1560)=  -0.797481477260590d0
      basgrd(1561)=   0.923867106437683d0
      basgrd(1562)=  -0.042534600943327d0
      basgrd(1563)=  -0.380342513322830d0
      basgrd(1564)=   0.829886376857758d0
      basgrd(1565)=  -0.111879602074623d0
      basgrd(1566)=  -0.546599984169006d0
      basgrd(1567)=   0.694963276386261d0
      basgrd(1568)=  -0.172340303659439d0
      basgrd(1569)=  -0.698086619377136d0
      basgrd(1570)=   0.534706711769104d0
      basgrd(1571)=  -0.215957403182983d0
      basgrd(1572)=  -0.816976785659790d0
      basgrd(1573)=   0.859103977680206d0
      basgrd(1574)=   0.085299097001553d0
      basgrd(1575)=  -0.504642784595490d0
      basgrd(1576)=   0.746298909187317d0
      basgrd(1577)=   0.005013999994844d0
      basgrd(1578)=  -0.665592074394226d0
      basgrd(1579)=   0.591585099697113d0
      basgrd(1580)=  -0.035555001348257d0
      basgrd(1581)=  -0.805458188056946d0
      basgrd(1582)=   0.764360010623932d0
      basgrd(1583)=   0.203857198357582d0
      basgrd(1584)=  -0.611715674400330d0
      basgrd(1585)=   0.628719925880432d0
      basgrd(1586)=   0.143068596720696d0
      basgrd(1587)=  -0.764357686042786d0
      basgrd(1588)=   0.646212697029114d0
      basgrd(1589)=   0.309264093637466d0
      basgrd(1590)=  -0.697685420513153d0
      basgrd(1591)=   0.917074501514435d0
      basgrd(1592)=  -0.396605193614960d0
      basgrd(1593)=  -0.040973000228405d0
      basgrd(1594)=   0.836342871189117d0
      basgrd(1595)=  -0.548160910606384d0
      basgrd(1596)=   0.007080500014126d0
      basgrd(1597)=   0.723484098911285d0
      basgrd(1598)=  -0.688323020935059d0
      basgrd(1599)=   0.052746899425983d0
      basgrd(1600)=   0.582683682441711d0
      basgrd(1601)=  -0.807337403297424d0
      basgrd(1602)=   0.093198999762535d0
      basgrd(1603)=   0.423074305057526d0
      basgrd(1604)=  -0.897111773490906d0
      basgrd(1605)=   0.127273604273796d0
      basgrd(1606)=   0.254787504673004d0
      basgrd(1607)=  -0.955131173133850d0
      basgrd(1608)=   0.151022106409073d0
      basgrd(1609)=   0.859105408191681d0
      basgrd(1610)=  -0.474321603775024d0
      basgrd(1611)=  -0.192242205142975d0
      basgrd(1612)=   0.749211311340332d0
      basgrd(1613)=  -0.644172012805939d0
      basgrd(1614)=  -0.154029101133347d0
      basgrd(1615)=   0.624307274818420d0
      basgrd(1616)=  -0.774409770965576d0
      basgrd(1617)=  -0.102615296840668d0
      basgrd(1618)=   0.472334414720535d0
      basgrd(1619)=  -0.878400802612305d0
      basgrd(1620)=  -0.072884596884251d0
      basgrd(1621)=   0.287012308835983d0
      basgrd(1622)=  -0.957587897777557d0
      basgrd(1623)=  -0.025481799617410d0
      basgrd(1624)=   0.768856704235077d0
      basgrd(1625)=  -0.542035818099976d0
      basgrd(1626)=  -0.339199811220169d0
      basgrd(1627)=   0.647287905216217d0
      basgrd(1628)=  -0.700507700443268d0
      basgrd(1629)=  -0.300511807203293d0
      basgrd(1630)=   0.489802598953247d0
      basgrd(1631)=  -0.833924591541290d0
      basgrd(1632)=  -0.254289686679840d0
      basgrd(1633)=   0.315028101205826d0
      basgrd(1634)=  -0.926266789436340d0
      basgrd(1635)=  -0.206850498914719d0
      basgrd(1636)=   0.648393511772156d0
      basgrd(1637)=  -0.595152378082275d0
      basgrd(1638)=  -0.474741488695145d0
      basgrd(1639)=   0.512678980827332d0
      basgrd(1640)=  -0.748812198638916d0
      basgrd(1641)=  -0.420048207044601d0
      basgrd(1642)=   0.335636794567108d0
      basgrd(1643)=  -0.860490083694458d0
      basgrd(1644)=  -0.383281499147415d0
      basgrd(1645)=   0.506496191024780d0
      basgrd(1646)=  -0.629669308662415d0
      basgrd(1647)=  -0.589048504829407d0
      basgrd(1648)=   0.347637087106705d0
      basgrd(1649)=  -0.764443278312683d0
      basgrd(1650)=  -0.542931795120239d0
      basgrd(1651)=   0.350451111793518d0
      basgrd(1652)=  -0.646095514297485d0
      basgrd(1653)=  -0.678044676780701d0
      basgrd(1654)=   0.927337884902954d0
      basgrd(1655)=  -0.363718807697296d0
      basgrd(1656)=   0.088051296770573d0
      basgrd(1657)=   0.894970297813416d0
      basgrd(1658)=  -0.359646290540695d0
      basgrd(1659)=   0.263974696397781d0
      basgrd(1660)=   0.829584300518036d0
      basgrd(1661)=  -0.346789509057999d0
      basgrd(1662)=   0.437637895345688d0
      basgrd(1663)=   0.731584012508392d0
      basgrd(1664)=  -0.326210111379623d0
      basgrd(1665)=   0.598641693592072d0
      basgrd(1666)=   0.608275115489960d0
      basgrd(1667)=  -0.297959387302399d0
      basgrd(1668)=   0.735677719116211d0
      basgrd(1669)=   0.467624694108963d0
      basgrd(1670)=  -0.265675097703934d0
      basgrd(1671)=   0.843056321144104d0
      basgrd(1672)=   0.845634222030640d0
      basgrd(1673)=  -0.517813503742218d0
      basgrd(1674)=   0.129506394267082d0
      basgrd(1675)=   0.786805391311646d0
      basgrd(1676)=  -0.523264110088348d0
      basgrd(1677)=   0.327310293912888d0
      basgrd(1678)=   0.709291279315949d0
      basgrd(1679)=  -0.500275671482086d0
      basgrd(1680)=   0.496618688106537d0
      basgrd(1681)=   0.595237672328949d0
      basgrd(1682)=  -0.481920689344406d0
      basgrd(1683)=   0.642996609210968d0
      basgrd(1684)=   0.448900699615479d0
      basgrd(1685)=  -0.432699590921402d0
      basgrd(1686)=   0.781830728054047d0
      basgrd(1687)=   0.732755780220032d0
      basgrd(1688)=  -0.658330380916596d0
      basgrd(1689)=   0.172249898314476d0
      basgrd(1690)=   0.660704910755158d0
      basgrd(1691)=  -0.657271027565002d0
      basgrd(1692)=   0.362579494714737d0
      basgrd(1693)=   0.551962077617645d0
      basgrd(1694)=  -0.633774876594544d0
      basgrd(1695)=   0.541910707950592d0
      basgrd(1696)=   0.418909907341003d0
      basgrd(1697)=  -0.591060817241669d0
      basgrd(1698)=   0.689319670200348d0
      basgrd(1699)=   0.591973423957825d0
      basgrd(1700)=  -0.777281701564789d0
      basgrd(1701)=   0.213074401021004d0
      basgrd(1702)=   0.509807109832764d0
      basgrd(1703)=  -0.757852792739868d0
      basgrd(1704)=   0.407131314277649d0
      basgrd(1705)=   0.375738888978958d0
      basgrd(1706)=  -0.731258690357208d0
      basgrd(1707)=   0.569281101226807d0
      basgrd(1708)=   0.433040112257004d0
      basgrd(1709)=  -0.866388618946075d0
      basgrd(1710)=   0.248690500855446d0
      basgrd(1711)=   0.323122292757034d0
      basgrd(1712)=  -0.843321800231934d0
      basgrd(1713)=   0.429418712854385d0
      basgrd(1714)=   0.264656305313110d0
      basgrd(1715)=  -0.923077106475830d0
      basgrd(1716)=   0.279080003499985d0
      basgrd(1717)=  -0.586927413940430d0
      basgrd(1718)=  -0.505001902580261d0
      basgrd(1719)=   0.632842183113098d0
      basgrd(1720)=  -0.479517698287964d0
      basgrd(1721)=  -0.644283294677734d0
      basgrd(1722)=   0.595786690711975d0
      basgrd(1723)=  -0.354927688837051d0
      basgrd(1724)=  -0.767434000968933d0
      basgrd(1725)=   0.533920705318451d0
      basgrd(1726)=  -0.217332303524017d0
      basgrd(1727)=  -0.865919888019562d0
      basgrd(1728)=   0.450499087572098d0
      basgrd(1729)=  -0.076692499220371d0
      basgrd(1730)=  -0.933489978313446d0
      basgrd(1731)=   0.350306510925293d0
      basgrd(1732)=   0.058503400534391d0
      basgrd(1733)=  -0.969352126121521d0
      basgrd(1734)=   0.238608002662659d0
      basgrd(1735)=  -0.650773108005524d0
      basgrd(1736)=  -0.584334790706635d0
      basgrd(1737)=   0.484816700220108d0
      basgrd(1738)=  -0.530069410800934d0
      basgrd(1739)=  -0.737646400928497d0
      basgrd(1740)=   0.418215602636337d0
      basgrd(1741)=  -0.394962996244431d0
      basgrd(1742)=  -0.848669886589050d0
      basgrd(1743)=   0.351800709962845d0
      basgrd(1744)=  -0.260236293077469d0
      basgrd(1745)=  -0.931690096855164d0
      basgrd(1746)=   0.253437608480454d0
      basgrd(1747)=  -0.092709496617317d0
      basgrd(1748)=  -0.985206186771393d0
      basgrd(1749)=   0.144130706787109d0
      basgrd(1750)=  -0.692767381668091d0
      basgrd(1751)=  -0.648679375648499d0
      basgrd(1752)=   0.315100491046906d0
      basgrd(1753)=  -0.566524803638458d0
      basgrd(1754)=  -0.788040578365326d0
      basgrd(1755)=   0.240918293595314d0
      basgrd(1756)=  -0.411148488521576d0
      basgrd(1757)=  -0.899503529071808d0
      basgrd(1758)=   0.147818699479103d0
      basgrd(1759)=  -0.249493896961212d0
      basgrd(1760)=  -0.967345297336578d0
      basgrd(1761)=   0.044675499200821d0
      basgrd(1762)=  -0.708103477954865d0
      basgrd(1763)=  -0.693778812885284d0
      basgrd(1764)=   0.131379097700119d0
      basgrd(1765)=  -0.560059428215027d0
      basgrd(1766)=  -0.826380729675293d0
      basgrd(1767)=   0.058552298694849d0
      basgrd(1768)=  -0.401598602533341d0
      basgrd(1769)=  -0.914239287376404d0
      basgrd(1770)=  -0.053712699562311d0
      basgrd(1771)=  -0.694803774356842d0
      basgrd(1772)=  -0.717216908931732d0
      basgrd(1773)=  -0.053362600505352d0
      basgrd(1774)=  -0.539078772068024d0
      basgrd(1775)=  -0.829476416110992d0
      basgrd(1776)=  -0.146160900592804d0
      basgrd(1777)=  -0.654999077320099d0
      basgrd(1778)=  -0.720001399517059d0
      basgrd(1779)=  -0.229290604591370d0
      basgrd(1780)=  -0.630659997463226d0
      basgrd(1781)=  -0.646287679672241d0
      basgrd(1782)=  -0.429627805948257d0
      basgrd(1783)=  -0.492116808891296d0
      basgrd(1784)=  -0.685429692268372d0
      basgrd(1785)=  -0.536662995815277d0
      basgrd(1786)=  -0.333649098873138d0
      basgrd(1787)=  -0.699857294559479d0
      basgrd(1788)=  -0.631567895412445d0
      basgrd(1789)=  -0.161868393421173d0
      basgrd(1790)=  -0.687184393405914d0
      basgrd(1791)=  -0.708220422267914d0
      basgrd(1792)=   0.011657600291073d0
      basgrd(1793)=  -0.648432612419128d0
      basgrd(1794)=  -0.761182785034180d0
      basgrd(1795)=   0.176732495427132d0
      basgrd(1796)=  -0.586371302604675d0
      basgrd(1797)=  -0.790527820587158d0
      basgrd(1798)=  -0.633521616458893d0
      basgrd(1799)=  -0.528598010540009d0
      basgrd(1800)=  -0.565008401870728d0
      basgrd(1801)=  -0.471201092004776d0
      basgrd(1802)=  -0.548658490180969d0
      basgrd(1803)=  -0.690610885620117d0
      basgrd(1804)=  -0.303806811571121d0
      basgrd(1805)=  -0.556329727172852d0
      basgrd(1806)=  -0.773433029651642d0
      basgrd(1807)=  -0.133980706334114d0
      basgrd(1808)=  -0.523911416530609d0
      basgrd(1809)=  -0.841169416904450d0
      basgrd(1810)=   0.064118698239326d0
      basgrd(1811)=  -0.476831287145615d0
      basgrd(1812)=  -0.876653075218201d0
      basgrd(1813)=  -0.612833380699158d0
      basgrd(1814)=  -0.390434205532074d0
      basgrd(1815)=  -0.687019884586334d0
      basgrd(1816)=  -0.445788890123367d0
      basgrd(1817)=  -0.399533212184906d0
      basgrd(1818)=  -0.801027715206146d0
      basgrd(1819)=  -0.252103894948959d0
      basgrd(1820)=  -0.386068314313889d0
      basgrd(1821)=  -0.887352705001831d0
      basgrd(1822)=  -0.058410998433828d0
      basgrd(1823)=  -0.348986893892288d0
      basgrd(1824)=  -0.935305416584015d0
      basgrd(1825)=  -0.566798210144043d0
      basgrd(1826)=  -0.237520307302475d0
      basgrd(1827)=  -0.788874983787537d0
      basgrd(1828)=  -0.381151497364044d0
      basgrd(1829)=  -0.249143004417419d0
      basgrd(1830)=  -0.890309691429138d0
      basgrd(1831)=  -0.183915302157402d0
      basgrd(1832)=  -0.210025101900101d0
      basgrd(1833)=  -0.960241973400116d0
      basgrd(1834)=  -0.498021095991135d0
      basgrd(1835)=  -0.081127703189850d0
      basgrd(1836)=  -0.863361597061157d0
      basgrd(1837)=  -0.303429692983627d0
      basgrd(1838)=  -0.067775696516037d0
      basgrd(1839)=  -0.950440406799316d0
      basgrd(1840)=  -0.410804599523544d0
      basgrd(1841)=   0.069998398423195d0
      basgrd(1842)=  -0.909032285213471d0
      basgrd(1843)=  -0.345934987068176d0
      basgrd(1844)=   0.274638891220093d0
      basgrd(1845)=  -0.897163629531860d0
      basgrd(1846)=  -0.178584501147270d0
      basgrd(1847)=   0.330264985561371d0
      basgrd(1848)=  -0.926840126514435d0
      basgrd(1849)=  -0.000102300000435d0
      basgrd(1850)=   0.380370706319809d0
      basgrd(1851)=  -0.924834072589874d0
      basgrd(1852)=   0.181090295314789d0
      basgrd(1853)=   0.420994311571121d0
      basgrd(1854)=  -0.888802587985992d0
      basgrd(1855)=   0.351623713970184d0
      basgrd(1856)=   0.449212610721588d0
      basgrd(1857)=  -0.821321427822113d0
      basgrd(1858)=   0.502149283885956d0
      basgrd(1859)=   0.467336088418961d0
      basgrd(1860)=  -0.727628290653229d0
      basgrd(1861)=  -0.333996385335922d0
      basgrd(1862)=   0.441805392503738d0
      basgrd(1863)=  -0.832619011402130d0
      basgrd(1864)=  -0.141809999942780d0
      basgrd(1865)=   0.516518115997315d0
      basgrd(1866)=  -0.844451904296875d0
      basgrd(1867)=   0.040631398558617d0
      basgrd(1868)=   0.555853307247162d0
      basgrd(1869)=  -0.830286800861359d0
      basgrd(1870)=   0.212933003902435d0
      basgrd(1871)=   0.596564710140228d0
      basgrd(1872)=  -0.773802280426025d0
      basgrd(1873)=   0.399469614028931d0
      basgrd(1874)=   0.608036875724793d0
      basgrd(1875)=  -0.686086773872376d0
      basgrd(1876)=  -0.305997014045715d0
      basgrd(1877)=   0.599825680255890d0
      basgrd(1878)=  -0.739307105541229d0
      basgrd(1879)=  -0.115020997822285d0
      basgrd(1880)=   0.666970014572144d0
      basgrd(1881)=  -0.736153006553650d0
      basgrd(1882)=   0.088313102722168d0
      basgrd(1883)=   0.712629616260529d0
      basgrd(1884)=  -0.695959627628326d0
      basgrd(1885)=   0.275941789150238d0
      basgrd(1886)=   0.732307791709900d0
      basgrd(1887)=  -0.622560381889343d0
      basgrd(1888)=  -0.263556391000748d0
      basgrd(1889)=   0.739698588848114d0
      basgrd(1890)=  -0.619180083274841d0
      basgrd(1891)=  -0.058363799005747d0
      basgrd(1892)=   0.791237711906433d0
      basgrd(1893)=  -0.608717083930969d0
      basgrd(1894)=   0.138180896639824d0
      basgrd(1895)=   0.830963075160980d0
      basgrd(1896)=  -0.538893699645996d0
      basgrd(1897)=  -0.208962306380272d0
      basgrd(1898)=   0.851629078388214d0
      basgrd(1899)=  -0.480689793825150d0
      basgrd(1900)=  -0.004034100100398d0
      basgrd(1901)=   0.897840201854706d0
      basgrd(1902)=  -0.440303087234497d0
      basgrd(1903)=  -0.144275605678558d0
      basgrd(1904)=   0.931676626205444d0
      basgrd(1905)=  -0.333411693572998d0
      basgrd(1906)=  -0.127263292670250d0
      basgrd(1907)=   0.983877420425415d0
      basgrd(1908)=  -0.125655800104141d0
      basgrd(1909)=   0.027428399771452d0
      basgrd(1910)=   0.998931407928467d0
      basgrd(1911)=  -0.037199299782515d0
      basgrd(1912)=   0.185099795460701d0
      basgrd(1913)=   0.980886816978455d0
      basgrd(1914)=   0.059992101043463d0
      basgrd(1915)=   0.337954908609390d0
      basgrd(1916)=   0.927699029445648d0
      basgrd(1917)=   0.158622100949287d0
      basgrd(1918)=   0.474094212055206d0
      basgrd(1919)=   0.843328595161438d0
      basgrd(1920)=   0.253044605255127d0
      basgrd(1921)=   0.585440695285797d0
      basgrd(1922)=   0.735126078128815d0
      basgrd(1923)=   0.341831505298615d0
      basgrd(1924)=  -0.166605100035667d0
      basgrd(1925)=   0.984770476818085d0
      basgrd(1926)=   0.049698099493980d0
      basgrd(1927)=   0.002889599883929d0
      basgrd(1928)=   0.985790193080902d0
      basgrd(1929)=   0.167956098914146d0
      basgrd(1930)=   0.163101702928543d0
      basgrd(1931)=   0.951984405517578d0
      basgrd(1932)=   0.259082108736038d0
      basgrd(1933)=   0.301893591880798d0
      basgrd(1934)=   0.881811320781708d0
      basgrd(1935)=   0.362310796976090d0
      basgrd(1936)=   0.449788898229599d0
      basgrd(1937)=   0.769920885562897d0
      basgrd(1938)=   0.452671796083450d0
      basgrd(1939)=  -0.196598693728447d0
      basgrd(1940)=   0.952963709831238d0
      basgrd(1941)=   0.230671092867851d0
      basgrd(1942)=  -0.030939500778914d0
      basgrd(1943)=   0.937760710716248d0
      basgrd(1944)=   0.345901310443878d0
      basgrd(1945)=   0.139792501926422d0
      basgrd(1946)=   0.878625273704529d0
      basgrd(1947)=   0.456591308116913d0
      basgrd(1948)=   0.291471391916275d0
      basgrd(1949)=   0.782800495624542d0
      basgrd(1950)=   0.549788892269135d0
      basgrd(1951)=  -0.217758700251579d0
      basgrd(1952)=   0.887192785739899d0
      basgrd(1953)=   0.406780213117600d0
      basgrd(1954)=  -0.037931300699711d0
      basgrd(1955)=   0.856921792030335d0
      basgrd(1956)=   0.514048874378204d0
      basgrd(1957)=   0.119260303676128d0
      basgrd(1958)=   0.769971191883087d0
      basgrd(1959)=   0.626834392547607d0
      basgrd(1960)=  -0.227417603135109d0
      basgrd(1961)=   0.792037606239319d0
      basgrd(1962)=   0.566531181335449d0
      basgrd(1963)=  -0.055118698626757d0
      basgrd(1964)=   0.732814371585846d0
      basgrd(1965)=   0.678192377090454d0
      basgrd(1966)=  -0.223621204495430d0
      basgrd(1967)=   0.675322592258453d0
      basgrd(1968)=   0.702803611755371d0
      basgrd(1969)=  -0.275744587182999d0
      basgrd(1970)=   0.502744495868683d0
      basgrd(1971)=   0.819275915622711d0
      basgrd(1972)=  -0.157824501395226d0
      basgrd(1973)=   0.397286087274551d0
      basgrd(1974)=   0.904021680355072d0
      basgrd(1975)=  -0.034025400876999d0
      basgrd(1976)=   0.271151989698410d0
      basgrd(1977)=   0.961934983730316d0
      basgrd(1978)=   0.091820001602173d0
      basgrd(1979)=   0.132203102111816d0
      basgrd(1980)=   0.986960709095001d0
      basgrd(1981)=   0.209251493215561d0
      basgrd(1982)=  -0.011194899678230d0
      basgrd(1983)=   0.977797806262970d0
      basgrd(1984)=   0.310767889022827d0
      basgrd(1985)=  -0.153308898210526d0
      basgrd(1986)=   0.938040316104889d0
      basgrd(1987)=  -0.361738592386246d0
      basgrd(1988)=   0.351538389921188d0
      basgrd(1989)=   0.863461613655090d0
      basgrd(1990)=  -0.236465096473694d0
      basgrd(1991)=   0.211114495992661d0
      basgrd(1992)=   0.948427617549896d0
      basgrd(1993)=  -0.105990499258041d0
      basgrd(1994)=   0.083933003246784d0
      basgrd(1995)=   0.990818500518799d0
      basgrd(1996)=   0.009114200249314d0
      basgrd(1997)=  -0.062731198966503d0
      basgrd(1998)=   0.997988820075989d0
      basgrd(1999)=   0.145382106304169d0
      basgrd(2000)=  -0.214628100395203d0
      basgrd(2001)=   0.965815126895905d0
      basgrd(2002)=  -0.435617089271545d0
      basgrd(2003)=   0.181357502937317d0
      basgrd(2004)=   0.881672978401184d0
      basgrd(2005)=  -0.310338914394379d0
      basgrd(2006)=   0.038345001637936d0
      basgrd(2007)=   0.949852287769318d0
      basgrd(2008)=  -0.168954193592072d0
      basgrd(2009)=  -0.117453500628471d0
      basgrd(2010)=   0.978600621223450d0
      basgrd(2011)=  -0.033441301435232d0
      basgrd(2012)=  -0.267472296953201d0
      basgrd(2013)=   0.962985098361969d0
      basgrd(2014)=  -0.492747992277145d0
      basgrd(2015)=   0.001143699977547d0
      basgrd(2016)=   0.870171308517456d0
      basgrd(2017)=  -0.348167985677719d0
      basgrd(2018)=  -0.142612397670746d0
      basgrd(2019)=   0.926520824432373d0
      basgrd(2020)=  -0.214163303375244d0
      basgrd(2021)=  -0.308155894279480d0
      basgrd(2022)=   0.926916420459747d0
      basgrd(2023)=  -0.527965784072876d0
      basgrd(2024)=  -0.177779898047447d0
      basgrd(2025)=   0.830449521541596d0
      basgrd(2026)=  -0.385218590497971d0
      basgrd(2027)=  -0.334098607301712d0
      basgrd(2028)=   0.860223710536957d0
      basgrd(2029)=  -0.539444386959076d0
      basgrd(2030)=  -0.345495492219925d0
      basgrd(2031)=   0.767875313758850d0
      basgrd(2032)=   0.127580299973488d0
      basgrd(2033)=  -0.983846724033356d0
      basgrd(2034)=   0.125574097037315d0
      basgrd(2035)=  -0.027820400893688d0
      basgrd(2036)=  -0.998920083045960d0
      basgrd(2037)=   0.037211801856756d0
      basgrd(2038)=  -0.185110807418823d0
      basgrd(2039)=  -0.980887889862061d0
      basgrd(2040)=  -0.059940099716187d0
      basgrd(2041)=  -0.337707608938217d0
      basgrd(2042)=  -0.927787125110626d0
      basgrd(2043)=  -0.158633798360825d0
      basgrd(2044)=  -0.473829090595245d0
      basgrd(2045)=  -0.843479394912720d0
      basgrd(2046)=  -0.253038406372070d0
      basgrd(2047)=  -0.585387587547302d0
      basgrd(2048)=  -0.735334277153015d0
      basgrd(2049)=  -0.341474413871765d0
      basgrd(2050)=   0.166604503989220d0
      basgrd(2051)=  -0.984772920608521d0
      basgrd(2052)=  -0.049651999026537d0
      basgrd(2053)=  -0.002672699978575d0
      basgrd(2054)=  -0.985792696475983d0
      basgrd(2055)=  -0.167945504188538d0
      basgrd(2056)=  -0.162859201431274d0
      basgrd(2057)=  -0.952022016048431d0
      basgrd(2058)=  -0.259096503257752d0
      basgrd(2059)=  -0.301810592412949d0
      basgrd(2060)=  -0.881809175014496d0
      basgrd(2061)=  -0.362384885549545d0
      basgrd(2062)=  -0.449979007244110d0
      basgrd(2063)=  -0.769826114177704d0
      basgrd(2064)=  -0.452643990516663d0
      basgrd(2065)=   0.196876093745232d0
      basgrd(2066)=  -0.952928602695465d0
      basgrd(2067)=  -0.230579406023026d0
      basgrd(2068)=   0.030882000923157d0
      basgrd(2069)=  -0.937787413597107d0
      basgrd(2070)=  -0.345833897590637d0
      basgrd(2071)=  -0.139927893877029d0
      basgrd(2072)=  -0.878561079502106d0
      basgrd(2073)=  -0.456673502922058d0
      basgrd(2074)=  -0.291802108287811d0
      basgrd(2075)=  -0.782654523849487d0
      basgrd(2076)=  -0.549821376800537d0
      basgrd(2077)=   0.217589899897575d0
      basgrd(2078)=  -0.887206315994263d0
      basgrd(2079)=  -0.406841099262238d0
      basgrd(2080)=   0.037573799490929d0
      basgrd(2081)=  -0.856916427612305d0
      basgrd(2082)=  -0.514084219932556d0
      basgrd(2083)=  -0.119242697954178d0
      basgrd(2084)=  -0.769980072975159d0
      basgrd(2085)=  -0.626826822757721d0
      basgrd(2086)=   0.227428093552589d0
      basgrd(2087)=  -0.792083680629730d0
      basgrd(2088)=  -0.566462695598602d0
      basgrd(2089)=   0.055615700781345d0
      basgrd(2090)=  -0.732818901538849d0
      basgrd(2091)=  -0.678147017955780d0
      basgrd(2092)=   0.223239198327065d0
      basgrd(2093)=  -0.675354480743408d0
      basgrd(2094)=  -0.702894508838654d0
      basgrd(2095)=   0.276006400585175d0
      basgrd(2096)=  -0.502714216709137d0
      basgrd(2097)=  -0.819206297397614d0
      basgrd(2098)=   0.157499894499779d0
      basgrd(2099)=  -0.397252798080444d0
      basgrd(2100)=  -0.904092907905579d0
      basgrd(2101)=   0.033979199826717d0
      basgrd(2102)=  -0.271220088005066d0
      basgrd(2103)=  -0.961917400360107d0
      basgrd(2104)=  -0.091511599719524d0
      basgrd(2105)=  -0.132233396172524d0
      basgrd(2106)=  -0.986985325813294d0
      basgrd(2107)=  -0.209296107292175d0
      basgrd(2108)=   0.011220600455999d0
      basgrd(2109)=  -0.977787911891937d0
      basgrd(2110)=  -0.310852795839310d0
      basgrd(2111)=   0.153104007244110d0
      basgrd(2112)=  -0.938045680522919d0
      basgrd(2113)=   0.361840397119522d0
      basgrd(2114)=  -0.351519405841827d0
      basgrd(2115)=  -0.863426685333252d0
      basgrd(2116)=   0.236789301037788d0
      basgrd(2117)=  -0.211201503872871d0
      basgrd(2118)=  -0.948327302932739d0
      basgrd(2119)=   0.105925798416138d0
      basgrd(2120)=  -0.083908401429653d0
      basgrd(2121)=  -0.990827500820160d0
      basgrd(2122)=  -0.009128700010478d0
      basgrd(2123)=   0.062666498124599d0
      basgrd(2124)=  -0.997992813587189d0
      basgrd(2125)=  -0.145236402750015d0
      basgrd(2126)=   0.214528605341911d0
      basgrd(2127)=  -0.965859115123749d0
      basgrd(2128)=   0.435572892427445d0
      basgrd(2129)=  -0.181442394852638d0
      basgrd(2130)=  -0.881677329540253d0
      basgrd(2131)=   0.310259014368057d0
      basgrd(2132)=  -0.038472600281239d0
      basgrd(2133)=  -0.949873328208923d0
      basgrd(2134)=   0.169092103838921d0
      basgrd(2135)=   0.117533296346664d0
      basgrd(2136)=  -0.978567183017731d0
      basgrd(2137)=   0.032988101243973d0
      basgrd(2138)=   0.267689198255539d0
      basgrd(2139)=  -0.962940394878388d0
      basgrd(2140)=   0.492677599191666d0
      basgrd(2141)=  -0.001045700046234d0
      basgrd(2142)=  -0.870211303234100d0
      basgrd(2143)=   0.347601801156998d0
      basgrd(2144)=   0.142801195383072d0
      basgrd(2145)=  -0.926704287528992d0
      basgrd(2146)=   0.214173600077629d0
      basgrd(2147)=   0.308104485273361d0
      basgrd(2148)=  -0.926931083202362d0
      basgrd(2149)=   0.528055191040039d0
      basgrd(2150)=   0.177592605352402d0
      basgrd(2151)=  -0.830432713031769d0
      basgrd(2152)=   0.385856300592423d0
      basgrd(2153)=   0.334045410156250d0
      basgrd(2154)=  -0.859958410263062d0
      basgrd(2155)=   0.539162278175354d0
      basgrd(2156)=   0.345473706722260d0
      basgrd(2157)=  -0.768083274364471d0
      basgrd(2158)=   0.587004601955414d0
      basgrd(2159)=   0.504969000816345d0
      basgrd(2160)=  -0.632797002792358d0
      basgrd(2161)=   0.479131013154984d0
      basgrd(2162)=   0.644578874111176d0
      basgrd(2163)=  -0.595778107643127d0
      basgrd(2164)=   0.354791104793549d0
      basgrd(2165)=   0.767505705356598d0
      basgrd(2166)=  -0.533908605575562d0
      basgrd(2167)=   0.217791900038719d0
      basgrd(2168)=   0.865857720375061d0
      basgrd(2169)=  -0.450396597385407d0
      basgrd(2170)=   0.076495200395584d0
      basgrd(2171)=   0.933516204357147d0
      basgrd(2172)=  -0.350279986858368d0
      basgrd(2173)=  -0.058539699763060d0
      basgrd(2174)=   0.969324111938477d0
      basgrd(2175)=  -0.238713204860687d0
      basgrd(2176)=   0.650557219982147d0
      basgrd(2177)=   0.584556996822357d0
      basgrd(2178)=  -0.484838485717773d0
      basgrd(2179)=   0.530629813671112d0
      basgrd(2180)=   0.737296521663666d0
      basgrd(2181)=  -0.418121904134750d0
      basgrd(2182)=   0.394735097885132d0
      basgrd(2183)=   0.848786771297455d0
      basgrd(2184)=  -0.351774305105209d0
      basgrd(2185)=   0.260277509689331d0
      basgrd(2186)=   0.931666195392609d0
      basgrd(2187)=  -0.253482997417450d0
      basgrd(2188)=   0.092900402843952d0
      basgrd(2189)=   0.985182523727417d0
      basgrd(2190)=  -0.144169703125954d0
      basgrd(2191)=   0.692610800266266d0
      basgrd(2192)=   0.648839890956879d0
      basgrd(2193)=  -0.315114408731461d0
      basgrd(2194)=   0.566652119159699d0
      basgrd(2195)=   0.787948906421661d0
      basgrd(2196)=  -0.240919202566147d0
      basgrd(2197)=   0.411309689283371d0
      basgrd(2198)=   0.899446010589600d0
      basgrd(2199)=  -0.147719904780388d0
      basgrd(2200)=   0.249073401093483d0
      basgrd(2201)=   0.967458724975586d0
      basgrd(2202)=  -0.044566400349140d0
      basgrd(2203)=   0.708077490329742d0
      basgrd(2204)=   0.693815410137177d0
      basgrd(2205)=  -0.131325706839562d0
      basgrd(2206)=   0.559486925601959d0
      basgrd(2207)=   0.826767027378082d0
      basgrd(2208)=  -0.058573000133038d0
      basgrd(2209)=   0.401702105998993d0
      basgrd(2210)=   0.914192378520966d0
      basgrd(2211)=   0.053737100213766d0
      basgrd(2212)=   0.695081770420075d0
      basgrd(2213)=   0.716940820217133d0
      basgrd(2214)=   0.053452998399735d0
      basgrd(2215)=   0.539405584335327d0
      basgrd(2216)=   0.829227089881897d0
      basgrd(2217)=   0.146369695663452d0
      basgrd(2218)=   0.655115306377411d0
      basgrd(2219)=   0.719930529594421d0
      basgrd(2220)=   0.229181200265884d0
      basgrd(2221)=   0.630646586418152d0
      basgrd(2222)=   0.646527409553528d0
      basgrd(2223)=   0.429286897182465d0
      basgrd(2224)=   0.491661399602890d0
      basgrd(2225)=   0.685511827468872d0
      basgrd(2226)=   0.536975383758545d0
      basgrd(2227)=   0.333570092916489d0
      basgrd(2228)=   0.699911475181580d0
      basgrd(2229)=   0.631549596786499d0
      basgrd(2230)=   0.162301704287529d0
      basgrd(2231)=   0.686992824077606d0
      basgrd(2232)=   0.708307206630707d0
      basgrd(2233)=  -0.011848499998450d0
      basgrd(2234)=   0.648522794246674d0
      basgrd(2235)=   0.761102974414825d0
      basgrd(2236)=  -0.176642000675201d0
      basgrd(2237)=   0.586501479148865d0
      basgrd(2238)=   0.790451526641846d0
      basgrd(2239)=   0.633333086967468d0
      basgrd(2240)=   0.528723776340485d0
      basgrd(2241)=   0.565102100372315d0
      basgrd(2242)=   0.471828311681747d0
      basgrd(2243)=   0.548452794551849d0
      basgrd(2244)=   0.690346002578735d0
      basgrd(2245)=   0.303616493940353d0
      basgrd(2246)=   0.556289911270142d0
      basgrd(2247)=   0.773536384105682d0
      basgrd(2248)=   0.134037494659424d0
      basgrd(2249)=   0.523956775665283d0
      basgrd(2250)=   0.841132104396820d0
      basgrd(2251)=  -0.064035296440125d0
      basgrd(2252)=   0.476818203926086d0
      basgrd(2253)=   0.876666307449341d0
      basgrd(2254)=   0.612858712673187d0
      basgrd(2255)=   0.390303999185562d0
      basgrd(2256)=   0.687071323394775d0
      basgrd(2257)=   0.445841103792191d0
      basgrd(2258)=   0.399573802947998d0
      basgrd(2259)=   0.800978481769562d0
      basgrd(2260)=   0.252257108688355d0
      basgrd(2261)=   0.386078298091888d0
      basgrd(2262)=   0.887304782867432d0
      basgrd(2263)=   0.058033999055624d0
      basgrd(2264)=   0.348991006612778d0
      basgrd(2265)=   0.935327410697937d0
      basgrd(2266)=   0.566856086254120d0
      basgrd(2267)=   0.237595707178116d0
      basgrd(2268)=   0.788810789585114d0
      basgrd(2269)=   0.380725204944611d0
      basgrd(2270)=   0.249311402440071d0
      basgrd(2271)=   0.890444874763489d0
      basgrd(2272)=   0.183980792760849d0
      basgrd(2273)=   0.210037499666214d0
      basgrd(2274)=   0.960226714611054d0
      basgrd(2275)=   0.498030602931976d0
      basgrd(2276)=   0.081083498895168d0
      basgrd(2277)=   0.863360285758972d0
      basgrd(2278)=   0.303744286298752d0
      basgrd(2279)=   0.067669302225113d0
      basgrd(2280)=   0.950347423553467d0
      basgrd(2281)=   0.410575896501541d0
      basgrd(2282)=  -0.069955699145794d0
      basgrd(2283)=   0.909138977527618d0
      basgrd(2284)=   0.346289485692978d0
      basgrd(2285)=  -0.274604201316834d0
      basgrd(2286)=   0.897037386894226d0
      basgrd(2287)=   0.178020507097244d0
      basgrd(2288)=  -0.330252110958099d0
      basgrd(2289)=   0.926953196525574d0
      basgrd(2290)=   0.000057599998399d0
      basgrd(2291)=  -0.380336195230484d0
      basgrd(2292)=   0.924848318099976d0
      basgrd(2293)=  -0.180613100528717d0
      basgrd(2294)=  -0.421217203140259d0
      basgrd(2295)=   0.888794124126434d0
      basgrd(2296)=  -0.351754099130631d0
      basgrd(2297)=  -0.449144095182419d0
      basgrd(2298)=   0.821303009986877d0
      basgrd(2299)=  -0.501990377902985d0
      basgrd(2300)=  -0.467135608196259d0
      basgrd(2301)=   0.727866828441620d0
      basgrd(2302)=   0.333940386772156d0
      basgrd(2303)=  -0.441711604595184d0
      basgrd(2304)=   0.832691192626953d0
      basgrd(2305)=   0.142343893647194d0
      basgrd(2306)=  -0.516568005084992d0
      basgrd(2307)=   0.844331502914429d0
      basgrd(2308)=  -0.040705099701881d0
      basgrd(2309)=  -0.555892288684845d0
      basgrd(2310)=   0.830257117748261d0
      basgrd(2311)=  -0.212945893406868d0
      basgrd(2312)=  -0.596570372581482d0
      basgrd(2313)=   0.773794412612915d0
      basgrd(2314)=  -0.399550497531891d0
      basgrd(2315)=  -0.607988595962524d0
      basgrd(2316)=   0.686082482337952d0
      basgrd(2317)=   0.306128799915314d0
      basgrd(2318)=  -0.599840879440308d0
      basgrd(2319)=   0.739240229129791d0
      basgrd(2320)=   0.115001700818539d0
      basgrd(2321)=  -0.666880786418915d0
      basgrd(2322)=   0.736236810684204d0
      basgrd(2323)=  -0.088348202407360d0
      basgrd(2324)=  -0.712625086307526d0
      basgrd(2325)=   0.695959806442261d0
      basgrd(2326)=  -0.276276499032974d0
      basgrd(2327)=  -0.732205629348755d0
      basgrd(2328)=   0.622532129287720d0
      basgrd(2329)=   0.263500094413757d0
      basgrd(2330)=  -0.739666819572449d0
      basgrd(2331)=   0.619242072105408d0
      basgrd(2332)=   0.058069601655006d0
      basgrd(2333)=  -0.791203618049622d0
      basgrd(2334)=   0.608789622783661d0
      basgrd(2335)=  -0.138212993741036d0
      basgrd(2336)=  -0.830910623073578d0
      basgrd(2337)=   0.538966298103333d0
      basgrd(2338)=   0.208887696266174d0
      basgrd(2339)=  -0.851623892784119d0
      basgrd(2340)=   0.480731397867203d0
      basgrd(2341)=   0.004429399967194d0
      basgrd(2342)=  -0.897843599319458d0
      basgrd(2343)=   0.440292209386826d0
      basgrd(2344)=   0.143959894776344d0
      basgrd(2345)=  -0.931712508201599d0
      basgrd(2346)=   0.333447694778442d0
      basgrd(2347)=  -0.994666099548340d0
      basgrd(2348)=   0.086979798972607d0
      basgrd(2349)=   0.055441800504923d0
      basgrd(2350)=  -0.996585428714752d0
      basgrd(2351)=  -0.081714503467083d0
      basgrd(2352)=  -0.011847499758005d0
      basgrd(2353)=  -0.963661193847656d0
      basgrd(2354)=  -0.255876392126083d0
      basgrd(2355)=  -0.076709598302841d0
      basgrd(2356)=  -0.895890772342682d0
      basgrd(2357)=  -0.425445914268494d0
      basgrd(2358)=  -0.127967000007629d0
      basgrd(2359)=  -0.800434529781342d0
      basgrd(2360)=  -0.575109720230103d0
      basgrd(2361)=  -0.168977692723274d0
      basgrd(2362)=  -0.980560123920441d0
      basgrd(2363)=  -0.030510300770402d0
      basgrd(2364)=   0.193832397460938d0
      basgrd(2365)=  -0.965697526931763d0
      basgrd(2366)=  -0.224637195467949d0
      basgrd(2367)=   0.130255401134491d0
      basgrd(2368)=  -0.915138483047485d0
      basgrd(2369)=  -0.396341711282730d0
      basgrd(2370)=   0.073720701038837d0
      basgrd(2371)=  -0.823475182056427d0
      basgrd(2372)=  -0.567232728004456d0
      basgrd(2373)=   0.011644899845123d0
      basgrd(2374)=  -0.932180881500244d0
      basgrd(2375)=  -0.155878603458405d0
      basgrd(2376)=   0.326711803674698d0
      basgrd(2377)=  -0.901177227497101d0
      basgrd(2378)=  -0.351090997457504d0
      basgrd(2379)=   0.254194289445877d0
      basgrd(2380)=  -0.820685088634491d0
      basgrd(2381)=  -0.535592675209045d0
      basgrd(2382)=   0.199038803577423d0
      basgrd(2383)=  -0.852145075798035d0
      basgrd(2384)=  -0.282854914665222d0
      basgrd(2385)=   0.440274685621262d0
      basgrd(2386)=  -0.793860673904419d0
      basgrd(2387)=  -0.477530390024185d0
      basgrd(2388)=   0.376496791839600d0
      basgrd(2389)=  -0.745989024639130d0
      basgrd(2390)=  -0.399931192398071d0
      basgrd(2391)=   0.532499313354492d0
      basgrd(2392)=  -0.966841995716095d0
      basgrd(2393)=   0.176573798060417d0
      basgrd(2394)=  -0.184494704008102d0
      basgrd(2395)=  -0.916858792304993d0
      basgrd(2396)=   0.174883499741554d0
      basgrd(2397)=  -0.358867406845093d0
      basgrd(2398)=  -0.834321916103363d0
      basgrd(2399)=   0.161830097436905d0
      basgrd(2400)=  -0.526989519596100d0
      basgrd(2401)=  -0.722634673118591d0
      basgrd(2402)=   0.133844107389450d0
      basgrd(2403)=  -0.678148090839386d0
      basgrd(2404)=  -0.592822670936585d0
      basgrd(2405)=   0.096927501261234d0
      basgrd(2406)=  -0.799478828907013d0
      basgrd(2407)=  -0.969990491867065d0
      basgrd(2408)=   0.003480199957266d0
      basgrd(2409)=  -0.243117898702622d0
      basgrd(2410)=  -0.900394797325134d0
      basgrd(2411)=  -0.014010200276971d0
      basgrd(2412)=  -0.434848189353943d0
      basgrd(2413)=  -0.803669989109039d0
      basgrd(2414)=  -0.035562299191952d0
      basgrd(2415)=  -0.594011723995209d0
      basgrd(2416)=  -0.666292071342468d0
      basgrd(2417)=  -0.060038000345230d0
      basgrd(2418)=  -0.743269979953766d0
      basgrd(2419)=  -0.937600612640381d0
      basgrd(2420)=  -0.173102304339409d0
      basgrd(2421)=  -0.301563888788223d0
      basgrd(2422)=  -0.850028574466705d0
      basgrd(2423)=  -0.185779392719269d0
      basgrd(2424)=  -0.492886811494827d0
      basgrd(2425)=  -0.723021507263184d0
      basgrd(2426)=  -0.218342900276184d0
      basgrd(2427)=  -0.655413091182709d0
      basgrd(2428)=  -0.869305491447449d0
      basgrd(2429)=  -0.339502394199371d0
      basgrd(2430)=  -0.359229892492294d0
      basgrd(2431)=  -0.759054303169251d0
      basgrd(2432)=  -0.363346904516220d0
      basgrd(2433)=  -0.540199577808380d0
      basgrd(2434)=  -0.772585391998291d0
      basgrd(2435)=  -0.485988110303879d0
      basgrd(2436)=  -0.408567488193512d0
      basgrd(2437)=  -0.890455782413483d0
      basgrd(2438)=   0.422193288803101d0
      basgrd(2439)=  -0.169827103614807d0
      basgrd(2440)=  -0.793963909149170d0
      basgrd(2441)=   0.570446074008942d0
      basgrd(2442)=  -0.210267797112465d0
      basgrd(2443)=  -0.666621923446655d0
      basgrd(2444)=   0.702109992504120d0
      basgrd(2445)=  -0.250313401222229d0
      basgrd(2446)=  -0.514521598815918d0
      basgrd(2447)=   0.806630611419678d0
      basgrd(2448)=  -0.290885895490646d0
      basgrd(2449)=  -0.351422190666199d0
      basgrd(2450)=   0.877690792083740d0
      basgrd(2451)=  -0.325823903083801d0
      basgrd(2452)=  -0.842958211898804d0
      basgrd(2453)=   0.411995887756348d0
      basgrd(2454)=  -0.345949202775955d0
      basgrd(2455)=  -0.720687687397003d0
      basgrd(2456)=   0.565805375576019d0
      basgrd(2457)=  -0.400591313838959d0
      basgrd(2458)=  -0.580860793590546d0
      basgrd(2459)=   0.683639883995056d0
      basgrd(2460)=  -0.441856592893601d0
      basgrd(2461)=  -0.405061900615692d0
      basgrd(2462)=   0.783099412918091d0
      basgrd(2463)=  -0.471890002489090d0
      basgrd(2464)=  -0.762207925319672d0
      basgrd(2465)=   0.392832309007645d0
      basgrd(2466)=  -0.514511287212372d0
      basgrd(2467)=  -0.623532474040985d0
      basgrd(2468)=   0.545085191726685d0
      basgrd(2469)=  -0.560436725616455d0
      basgrd(2470)=  -0.451927989721298d0
      basgrd(2471)=   0.656245291233063d0
      basgrd(2472)=  -0.604237675666809d0
      basgrd(2473)=  -0.649629712104797d0
      basgrd(2474)=   0.370478004217148d0
      basgrd(2475)=  -0.663872897624970d0
      basgrd(2476)=  -0.490397810935974d0
      basgrd(2477)=   0.505023777484894d0
      basgrd(2478)=  -0.710254192352295d0
      basgrd(2479)=  -0.516724109649658d0
      basgrd(2480)=   0.342201292514801d0
      basgrd(2481)=  -0.784789383411408d0
      basgrd(2482)=  -0.870988070964813d0
      basgrd(2483)=   0.484789192676544d0
      basgrd(2484)=   0.079745002090931d0
      basgrd(2485)=  -0.798328280448914d0
      basgrd(2486)=   0.556500375270844d0
      basgrd(2487)=   0.230172097682953d0
      basgrd(2488)=  -0.693305075168610d0
      basgrd(2489)=   0.615923881530762d0
      basgrd(2490)=   0.374119907617569d0
      basgrd(2491)=  -0.559368014335632d0
      basgrd(2492)=   0.662029683589935d0
      basgrd(2493)=   0.498822808265686d0
      basgrd(2494)=  -0.410015493631363d0
      basgrd(2495)=   0.688486278057098d0
      basgrd(2496)=   0.598225712776184d0
      basgrd(2497)=  -0.774665474891663d0
      basgrd(2498)=   0.631620585918427d0
      basgrd(2499)=   0.030803799629211d0
      basgrd(2500)=  -0.674897789955139d0
      basgrd(2501)=   0.713388681411743d0
      basgrd(2502)=   0.188651904463768d0
      basgrd(2503)=  -0.555116176605225d0
      basgrd(2504)=   0.766797304153442d0
      basgrd(2505)=   0.322285592556000d0
      basgrd(2506)=  -0.401009112596512d0
      basgrd(2507)=   0.796635508537293d0
      basgrd(2508)=   0.452286988496780d0
      basgrd(2509)=  -0.648038923740387d0
      basgrd(2510)=   0.761452913284302d0
      basgrd(2511)=  -0.015333199873567d0
      basgrd(2512)=  -0.534516811370850d0
      basgrd(2513)=   0.832108080387116d0
      basgrd(2514)=   0.147945895791054d0
      basgrd(2515)=  -0.382510006427765d0
      basgrd(2516)=   0.879498720169067d0
      basgrd(2517)=   0.283139914274216d0
      basgrd(2518)=  -0.496056407690048d0
      basgrd(2519)=   0.866828083992004d0
      basgrd(2520)=  -0.050369501113892d0
      basgrd(2521)=  -0.359728902578354d0
      basgrd(2522)=   0.927440702915192d0
      basgrd(2523)=   0.102219998836517d0
      basgrd(2524)=  -0.331955105066299d0
      basgrd(2525)=   0.940168321132660d0
      basgrd(2526)=  -0.076742097735405d0
      basgrd(2527)=  -0.935442328453064d0
      basgrd(2528)=   0.277412593364716d0
      basgrd(2529)=   0.219066098332405d0
      basgrd(2530)=  -0.923786282539368d0
      basgrd(2531)=   0.152717202901840d0
      basgrd(2532)=   0.351135909557343d0
      basgrd(2533)=  -0.877286911010742d0
      basgrd(2534)=   0.022936899214983d0
      basgrd(2535)=   0.479418009519577d0
      basgrd(2536)=  -0.795153498649597d0
      basgrd(2537)=  -0.099442802369595d0
      basgrd(2538)=   0.598199009895325d0
      basgrd(2539)=  -0.687673687934876d0
      basgrd(2540)=  -0.209589198231697d0
      basgrd(2541)=   0.695109486579895d0
      basgrd(2542)=  -0.860243976116180d0
      basgrd(2543)=   0.356409996747971d0
      basgrd(2544)=   0.364625900983810d0
      basgrd(2545)=  -0.826919317245483d0
      basgrd(2546)=   0.223306298255920d0
      basgrd(2547)=   0.516080200672150d0
      basgrd(2548)=  -0.762157976627350d0
      basgrd(2549)=   0.098359100520611d0
      basgrd(2550)=   0.639875471591950d0
      basgrd(2551)=  -0.660059690475464d0
      basgrd(2552)=  -0.038458198308945d0
      basgrd(2553)=   0.750227987766266d0
      basgrd(2554)=  -0.753583610057831d0
      basgrd(2555)=   0.420744299888611d0
      basgrd(2556)=   0.505060374736786d0
      basgrd(2557)=  -0.706661522388458d0
      basgrd(2558)=   0.276728004217148d0
      basgrd(2559)=   0.651192128658295d0
      basgrd(2560)=  -0.610786616802216d0
      basgrd(2561)=   0.142595499753952d0
      basgrd(2562)=   0.778849303722382d0
      basgrd(2563)=  -0.621441781520844d0
      basgrd(2564)=   0.461835205554962d0
      basgrd(2565)=   0.632865190505981d0
      basgrd(2566)=  -0.547401905059815d0
      basgrd(2567)=   0.319857388734818d0
      basgrd(2568)=   0.773331999778748d0
      basgrd(2569)=  -0.473477303981781d0
      basgrd(2570)=   0.481519699096680d0
      basgrd(2571)=   0.737535119056702d0
      basgrd(2572)=   0.966974079608917d0
      basgrd(2573)=  -0.176420897245407d0
      basgrd(2574)=   0.183947503566742d0
      basgrd(2575)=   0.969919681549072d0
      basgrd(2576)=  -0.004127399995923d0
      basgrd(2577)=   0.243389904499054d0
      basgrd(2578)=   0.937752902507782d0
      basgrd(2579)=   0.172574996948242d0
      basgrd(2580)=   0.301392287015915d0
      basgrd(2581)=   0.869104385375977d0
      basgrd(2582)=   0.339699387550354d0
      basgrd(2583)=   0.359530001878738d0
      basgrd(2584)=   0.772198915481567d0
      basgrd(2585)=   0.486423194408417d0
      basgrd(2586)=   0.408780395984650d0
      basgrd(2587)=   0.916698694229126d0
      basgrd(2588)=  -0.174796804785729d0
      basgrd(2589)=   0.359318196773529d0
      basgrd(2590)=   0.900838792324066d0
      basgrd(2591)=   0.014590299688280d0
      basgrd(2592)=   0.433908492326737d0
      basgrd(2593)=   0.849835395812988d0
      basgrd(2594)=   0.185787096619606d0
      basgrd(2595)=   0.493216902017593d0
      basgrd(2596)=   0.759419083595276d0
      basgrd(2597)=   0.363279610872269d0
      basgrd(2598)=   0.539731919765472d0
      basgrd(2599)=   0.833809494972229d0
      basgrd(2600)=  -0.162178695201874d0
      basgrd(2601)=   0.527692973613739d0
      basgrd(2602)=   0.803708612918854d0
      basgrd(2603)=   0.035674199461937d0
      basgrd(2604)=   0.593952715396881d0
      basgrd(2605)=   0.722991406917572d0
      basgrd(2606)=   0.218149796128273d0
      basgrd(2607)=   0.655510604381561d0
      basgrd(2608)=   0.723245203495026d0
      basgrd(2609)=  -0.133661493659020d0
      basgrd(2610)=   0.677532970905304d0
      basgrd(2611)=   0.666010797023773d0
      basgrd(2612)=   0.060163900256157d0
      basgrd(2613)=   0.743511915206909d0
      basgrd(2614)=   0.593019783496857d0
      basgrd(2615)=  -0.097018897533417d0
      basgrd(2616)=   0.799321472644806d0
      basgrd(2617)=   0.994576513767242d0
      basgrd(2618)=  -0.087965399026871d0
      basgrd(2619)=  -0.055494401603937d0
      basgrd(2620)=   0.980652213096619d0
      basgrd(2621)=   0.030827200040221d0
      basgrd(2622)=  -0.193315505981445d0
      basgrd(2623)=   0.932510316371918d0
      basgrd(2624)=   0.155420094728470d0
      basgrd(2625)=  -0.325989395380020d0
      basgrd(2626)=   0.851574122905731d0
      basgrd(2627)=   0.283373802900314d0
      basgrd(2628)=  -0.441045194864273d0
      basgrd(2629)=   0.745563089847565d0
      basgrd(2630)=   0.400366485118866d0
      basgrd(2631)=  -0.532768487930298d0
      basgrd(2632)=   0.996565222740173d0
      basgrd(2633)=   0.081868700683117d0
      basgrd(2634)=   0.012464299798012d0
      basgrd(2635)=   0.965731620788574d0
      basgrd(2636)=   0.223911300301552d0
      basgrd(2637)=  -0.131248906254768d0
      basgrd(2638)=   0.900889098644257d0
      basgrd(2639)=   0.351737409830093d0
      basgrd(2640)=  -0.254321813583374d0
      basgrd(2641)=   0.794312477111816d0
      basgrd(2642)=   0.476961404085159d0
      basgrd(2643)=  -0.376265197992325d0
      basgrd(2644)=   0.963449299335480d0
      basgrd(2645)=   0.256639689207077d0
      basgrd(2646)=   0.076821200549603d0
      basgrd(2647)=   0.915381491184235d0
      basgrd(2648)=   0.395770490169525d0
      basgrd(2649)=  -0.073773503303528d0
      basgrd(2650)=   0.820885121822357d0
      basgrd(2651)=   0.535341203212738d0
      basgrd(2652)=  -0.198890507221222d0
      basgrd(2653)=   0.896209776401520d0
      basgrd(2654)=   0.424769699573517d0
      basgrd(2655)=   0.127979397773743d0
      basgrd(2656)=   0.823020279407501d0
      basgrd(2657)=   0.567890226840973d0
      basgrd(2658)=  -0.011760899797082d0
      basgrd(2659)=   0.801004528999329d0
      basgrd(2660)=   0.574423193931580d0
      basgrd(2661)=   0.168611094355583d0
      basgrd(2662)=   0.935489177703857d0
      basgrd(2663)=  -0.277362197637558d0
      basgrd(2664)=  -0.218929395079613d0
      basgrd(2665)=   0.860543191432953d0
      basgrd(2666)=  -0.356018602848053d0
      basgrd(2667)=  -0.364302307367325d0
      basgrd(2668)=   0.754030227661133d0
      basgrd(2669)=  -0.420134902000427d0
      basgrd(2670)=  -0.504901111125946d0
      basgrd(2671)=   0.620874106884003d0
      basgrd(2672)=  -0.462294399738312d0
      basgrd(2673)=  -0.633087098598480d0
      basgrd(2674)=   0.473358809947968d0
      basgrd(2675)=  -0.481462687253952d0
      basgrd(2676)=  -0.737648427486420d0
      basgrd(2677)=   0.923859477043152d0
      basgrd(2678)=  -0.152065306901932d0
      basgrd(2679)=  -0.351226210594177d0
      basgrd(2680)=   0.826567173004150d0
      basgrd(2681)=  -0.224111095070839d0
      basgrd(2682)=  -0.516295373439789d0
      basgrd(2683)=   0.706711888313294d0
      basgrd(2684)=  -0.276650607585907d0
      basgrd(2685)=  -0.651170313358307d0
      basgrd(2686)=   0.547578275203705d0
      basgrd(2687)=  -0.319883197546005d0
      basgrd(2688)=  -0.773196399211884d0
      basgrd(2689)=   0.877125680446625d0
      basgrd(2690)=  -0.022385800257325d0
      basgrd(2691)=  -0.479738891124725d0
      basgrd(2692)=   0.762327373027802d0
      basgrd(2693)=  -0.098271697759628d0
      basgrd(2694)=  -0.639687120914459d0
      basgrd(2695)=   0.610936522483826d0
      basgrd(2696)=  -0.142412006855011d0
      basgrd(2697)=  -0.778765320777893d0
      basgrd(2698)=   0.795427381992340d0
      basgrd(2699)=   0.099295802414417d0
      basgrd(2700)=  -0.597859203815460d0
      basgrd(2701)=   0.659589886665344d0
      basgrd(2702)=   0.038412500172853d0
      basgrd(2703)=  -0.750643491744995d0
      basgrd(2704)=   0.688159823417664d0
      basgrd(2705)=   0.209255293011665d0
      basgrd(2706)=  -0.694728910923004d0
      basgrd(2707)=   0.870572984218597d0
      basgrd(2708)=  -0.485449790954590d0
      basgrd(2709)=  -0.080257102847099d0
      basgrd(2710)=   0.775140702724457d0
      basgrd(2711)=  -0.631040215492249d0
      basgrd(2712)=  -0.030744200572371d0
      basgrd(2713)=   0.648166894912720d0
      basgrd(2714)=  -0.761349678039551d0
      basgrd(2715)=   0.015042600221932d0
      basgrd(2716)=   0.495720297098160d0
      basgrd(2717)=  -0.867007970809937d0
      basgrd(2718)=   0.050581701099873d0
      basgrd(2719)=   0.331745803356171d0
      basgrd(2720)=  -0.940250277519226d0
      basgrd(2721)=   0.076642997562885d0
      basgrd(2722)=   0.798034608364105d0
      basgrd(2723)=  -0.556702315807343d0
      basgrd(2724)=  -0.230701997876167d0
      basgrd(2725)=   0.674437284469605d0
      basgrd(2726)=  -0.713971376419067d0
      basgrd(2727)=  -0.188093706965447d0
      basgrd(2728)=   0.534844994544983d0
      basgrd(2729)=  -0.831898212432861d0
      basgrd(2730)=  -0.147939696907997d0
      basgrd(2731)=   0.359872907400131d0
      basgrd(2732)=  -0.927390694618225d0
      basgrd(2733)=  -0.102166503667831d0
      basgrd(2734)=   0.693240523338318d0
      basgrd(2735)=  -0.615791678428650d0
      basgrd(2736)=  -0.374457299709320d0
      basgrd(2737)=   0.554978728294373d0
      basgrd(2738)=  -0.766792714595795d0
      basgrd(2739)=  -0.322533011436462d0
      basgrd(2740)=   0.382354199886322d0
      basgrd(2741)=  -0.879494905471802d0
      basgrd(2742)=  -0.283362001180649d0
      basgrd(2743)=   0.559640884399414d0
      basgrd(2744)=  -0.661980807781220d0
      basgrd(2745)=  -0.498581409454346d0
      basgrd(2746)=   0.400750011205673d0
      basgrd(2747)=  -0.796752691268921d0
      basgrd(2748)=  -0.452310293912888d0
      basgrd(2749)=   0.410415202379227d0
      basgrd(2750)=  -0.688285827636719d0
      basgrd(2751)=  -0.598182201385498d0
      basgrd(2752)=   0.890244483947754d0
      basgrd(2753)=  -0.422470897436142d0
      basgrd(2754)=   0.170244395732880d0
      basgrd(2755)=   0.843210518360138d0
      basgrd(2756)=  -0.411951392889023d0
      basgrd(2757)=   0.345387011766434d0
      basgrd(2758)=   0.762408614158630d0
      basgrd(2759)=  -0.392958194017410d0
      basgrd(2760)=   0.514117717742920d0
      basgrd(2761)=   0.649405717849731d0
      basgrd(2762)=  -0.370465189218521d0
      basgrd(2763)=   0.664099216461182d0
      basgrd(2764)=   0.516222000122070d0
      basgrd(2765)=  -0.342312902212143d0
      basgrd(2766)=   0.785071074962616d0
      basgrd(2767)=   0.793639123439789d0
      basgrd(2768)=  -0.570809423923492d0
      basgrd(2769)=   0.210508003830910d0
      basgrd(2770)=   0.720682322978973d0
      basgrd(2771)=  -0.565303385257721d0
      basgrd(2772)=   0.401309192180634d0
      basgrd(2773)=   0.623670816421509d0
      basgrd(2774)=  -0.545060396194458d0
      basgrd(2775)=   0.560306906700134d0
      basgrd(2776)=   0.490850687026978d0
      basgrd(2777)=  -0.504851579666138d0
      basgrd(2778)=   0.710063695907593d0
      basgrd(2779)=   0.666391670703888d0
      basgrd(2780)=  -0.702421307563782d0
      basgrd(2781)=   0.250052690505981d0
      basgrd(2782)=   0.580631792545319d0
      basgrd(2783)=  -0.683867216110230d0
      basgrd(2784)=   0.441805809736252d0
      basgrd(2785)=   0.451695293188095d0
      basgrd(2786)=  -0.656417489051819d0
      basgrd(2787)=   0.604224681854248d0
      basgrd(2788)=   0.514898717403412d0
      basgrd(2789)=  -0.806341171264648d0
      basgrd(2790)=   0.291020989418030d0
      basgrd(2791)=   0.404949486255646d0
      basgrd(2792)=  -0.783091127872467d0
      basgrd(2793)=   0.472000211477280d0
      basgrd(2794)=   0.351482599973679d0
      basgrd(2795)=  -0.877726912498474d0
      basgrd(2796)=   0.325661689043045d0
      basgrd(2797)=  -0.574198603630066d0
      basgrd(2798)=  -0.590250492095947d0
      basgrd(2799)=   0.567362606525421d0
      basgrd(2800)=  -0.457240492105484d0
      basgrd(2801)=  -0.721351325511932d0
      basgrd(2802)=   0.520176291465759d0
      basgrd(2803)=  -0.324115395545960d0
      basgrd(2804)=  -0.832714617252350d0
      basgrd(2805)=   0.448927313089371d0
      basgrd(2806)=  -0.185243204236031d0
      basgrd(2807)=  -0.916581511497498d0
      basgrd(2808)=   0.354349195957184d0
      basgrd(2809)=  -0.049209598451853d0
      basgrd(2810)=  -0.967843472957611d0
      basgrd(2811)=   0.246692895889282d0
      basgrd(2812)=  -0.629358410835266d0
      basgrd(2813)=  -0.660773813724518d0
      basgrd(2814)=   0.409006088972092d0
      basgrd(2815)=  -0.498923689126968d0
      basgrd(2816)=  -0.799269974231720d0
      basgrd(2817)=   0.335026293992996d0
      basgrd(2818)=  -0.362767785787582d0
      basgrd(2819)=  -0.897077381610870d0
      basgrd(2820)=   0.252292811870575d0
      basgrd(2821)=  -0.201245903968811d0
      basgrd(2822)=  -0.968324601650238d0
      basgrd(2823)=   0.147809401154518d0
      basgrd(2824)=  -0.661052227020264d0
      basgrd(2825)=  -0.713988482952118d0
      basgrd(2826)=   0.230717003345490d0
      basgrd(2827)=  -0.513431727886200d0
      basgrd(2828)=  -0.844123184680939d0
      basgrd(2829)=   0.154414802789688d0
      basgrd(2830)=  -0.356402188539505d0
      basgrd(2831)=  -0.933194696903229d0
      basgrd(2832)=   0.046098999679089d0
      basgrd(2833)=  -0.659991085529327d0
      basgrd(2834)=  -0.749914407730103d0
      basgrd(2835)=   0.045169100165367d0
      basgrd(2836)=  -0.503252208232880d0
      basgrd(2837)=  -0.862800180912018d0
      basgrd(2838)=  -0.048094399273396d0
      basgrd(2839)=  -0.630755305290222d0
      basgrd(2840)=  -0.764443397521973d0
      basgrd(2841)=  -0.133319199085236d0
      basgrd(2842)=  -0.583454608917236d0
      basgrd(2843)=  -0.620762407779694d0
      basgrd(2844)=  -0.523674309253693d0
      basgrd(2845)=  -0.434691786766052d0
      basgrd(2846)=  -0.649775922298431d0
      basgrd(2847)=  -0.623565793037415d0
      basgrd(2848)=  -0.268127799034119d0
      basgrd(2849)=  -0.652683794498444d0
      basgrd(2850)=  -0.708598077297211d0
      basgrd(2851)=  -0.094741098582745d0
      basgrd(2852)=  -0.624797701835632d0
      basgrd(2853)=  -0.775017380714417d0
      basgrd(2854)=   0.073038898408413d0
      basgrd(2855)=  -0.572193205356598d0
      basgrd(2856)=  -0.816860020160675d0
      basgrd(2857)=  -0.577180325984955d0
      basgrd(2858)=  -0.492267608642578d0
      basgrd(2859)=  -0.651563882827759d0
      basgrd(2860)=  -0.406243413686752d0
      basgrd(2861)=  -0.502484202384949d0
      basgrd(2862)=  -0.763201117515564d0
      basgrd(2863)=  -0.236425399780273d0
      basgrd(2864)=  -0.489807188510895d0
      basgrd(2865)=  -0.839161515235901d0
      basgrd(2866)=  -0.042014699429274d0
      basgrd(2867)=  -0.452830493450165d0
      basgrd(2868)=  -0.890606105327606d0
      basgrd(2869)=  -0.546594023704529d0
      basgrd(2870)=  -0.344131797552109d0
      basgrd(2871)=  -0.763418793678284d0
      basgrd(2872)=  -0.360694497823715d0
      basgrd(2873)=  -0.351880192756653d0
      basgrd(2874)=  -0.863759100437164d0
      basgrd(2875)=  -0.165704295039177d0
      basgrd(2876)=  -0.316727906465530d0
      basgrd(2877)=  -0.933930218219757d0
      basgrd(2878)=  -0.486211985349655d0
      basgrd(2879)=  -0.189593702554703d0
      basgrd(2880)=  -0.853025317192078d0
      basgrd(2881)=  -0.291031688451767d0
      basgrd(2882)=  -0.176190003752708d0
      basgrd(2883)=  -0.940349817276001d0
      basgrd(2884)=  -0.406025707721710d0
      basgrd(2885)=  -0.038248799741268d0
      basgrd(2886)=  -0.913060903549194d0
      basgrd(2887)=  -0.280765205621719d0
      basgrd(2888)=   0.360294699668884d0
      basgrd(2889)=  -0.889583408832550d0
      basgrd(2890)=  -0.107327498495579d0
      basgrd(2891)=   0.411116003990173d0
      basgrd(2892)=  -0.905242800712585d0
      basgrd(2893)=   0.074388600885868d0
      basgrd(2894)=   0.454668313264847d0
      basgrd(2895)=  -0.887548923492432d0
      basgrd(2896)=   0.251096993684769d0
      basgrd(2897)=   0.491504490375519d0
      basgrd(2898)=  -0.833890616893768d0
      basgrd(2899)=   0.409645587205887d0
      basgrd(2900)=   0.516426026821137d0
      basgrd(2901)=  -0.751993775367737d0
      basgrd(2902)=  -0.263300687074661d0
      basgrd(2903)=   0.522838175296783d0
      basgrd(2904)=  -0.810748398303986d0
      basgrd(2905)=  -0.068727202713490d0
      basgrd(2906)=   0.586754322052002d0
      basgrd(2907)=  -0.806843221187592d0
      basgrd(2908)=   0.109495803713799d0
      basgrd(2909)=   0.627098917961121d0
      basgrd(2910)=  -0.771205306053162d0
      basgrd(2911)=   0.299069494009018d0
      basgrd(2912)=   0.650068700313568d0
      basgrd(2913)=  -0.698547184467316d0
      basgrd(2914)=  -0.230985701084137d0
      basgrd(2915)=   0.673225700855255d0
      basgrd(2916)=  -0.702433407306671d0
      basgrd(2917)=  -0.026773400604725d0
      basgrd(2918)=   0.724989295005798d0
      basgrd(2919)=  -0.688239574432373d0
      basgrd(2920)=   0.168925702571869d0
      basgrd(2921)=   0.763961374759674d0
      basgrd(2922)=  -0.622757673263550d0
      basgrd(2923)=  -0.179126098752022d0
      basgrd(2924)=   0.799882709980011d0
      basgrd(2925)=  -0.572801411151886d0
      basgrd(2926)=   0.025600800290704d0
      basgrd(2927)=   0.846191823482513d0
      basgrd(2928)=  -0.532263100147247d0
      basgrd(2929)=  -0.116815298795700d0
      basgrd(2930)=   0.894683480262756d0
      basgrd(2931)=  -0.431156098842621d0
      basgrd(2932)=  -0.084826499223709d0
      basgrd(2933)=   0.996023297309876d0
      basgrd(2934)=  -0.027241699397564d0
      basgrd(2935)=   0.072470501065254d0
      basgrd(2936)=   0.995377600193024d0
      basgrd(2937)=   0.063018903136253d0
      basgrd(2938)=   0.230611503124237d0
      basgrd(2939)=   0.959879517555237d0
      basgrd(2940)=   0.159529805183411d0
      basgrd(2941)=   0.374848186969757d0
      basgrd(2942)=   0.890116393566132d0
      basgrd(2943)=   0.259193986654282d0
      basgrd(2944)=   0.496235489845276d0
      basgrd(2945)=   0.793621420860291d0
      basgrd(2946)=   0.352016091346741d0
      basgrd(2947)=  -0.121894098818302d0
      basgrd(2948)=   0.981123626232147d0
      basgrd(2949)=   0.150127902626991d0
      basgrd(2950)=   0.047682300209999d0
      basgrd(2951)=   0.963326275348663d0
      basgrd(2952)=   0.264062196016312d0
      basgrd(2953)=   0.197808593511581d0
      basgrd(2954)=   0.911192119121552d0
      basgrd(2955)=   0.361387193202972d0
      basgrd(2956)=   0.350805908441544d0
      basgrd(2957)=   0.816609025001526d0
      basgrd(2958)=   0.458350300788879d0
      basgrd(2959)=  -0.150530397891998d0
      basgrd(2960)=   0.932001292705536d0
      basgrd(2961)=   0.329718202352524d0
      basgrd(2962)=   0.027084600180388d0
      basgrd(2963)=   0.898568928241730d0
      basgrd(2964)=   0.437995791435242d0
      basgrd(2965)=   0.184851095080376d0
      basgrd(2966)=   0.815399289131165d0
      basgrd(2967)=   0.548592805862427d0
      basgrd(2968)=  -0.163372099399567d0
      basgrd(2969)=   0.850857079029083d0
      basgrd(2970)=   0.499351412057877d0
      basgrd(2971)=   0.008724099956453d0
      basgrd(2972)=   0.791486620903015d0
      basgrd(2973)=   0.611124217510223d0
      basgrd(2974)=  -0.163062706589699d0
      basgrd(2975)=   0.745155990123749d0
      basgrd(2976)=   0.646647572517395d0
      basgrd(2977)=  -0.265450000762940d0
      basgrd(2978)=   0.409535109996796d0
      basgrd(2979)=   0.872821509838104d0
      basgrd(2980)=  -0.143057703971863d0
      basgrd(2981)=   0.295968502759934d0
      basgrd(2982)=   0.944424211978912d0
      basgrd(2983)=  -0.015468999743462d0
      basgrd(2984)=   0.164072602987289d0
      basgrd(2985)=   0.986326992511749d0
      basgrd(2986)=   0.105013199150562d0
      basgrd(2987)=   0.019827900454402d0
      basgrd(2988)=   0.994273126125336d0
      basgrd(2989)=   0.212311506271362d0
      basgrd(2990)=  -0.123601302504540d0
      basgrd(2991)=   0.969353675842285d0
      basgrd(2992)=  -0.347634613513947d0
      basgrd(2993)=   0.250227689743042d0
      basgrd(2994)=   0.903623998165131d0
      basgrd(2995)=  -0.218203797936440d0
      basgrd(2996)=   0.106527402997017d0
      basgrd(2997)=   0.970071613788605d0
      basgrd(2998)=  -0.094216100871563d0
      basgrd(2999)=  -0.030725199729204d0
      basgrd(3000)=   0.995077490806580d0
      basgrd(3001)=   0.041333999484777d0
      basgrd(3002)=  -0.183548793196678d0
      basgrd(3003)=   0.982141196727753d0
      basgrd(3004)=  -0.416522413492203d0
      basgrd(3005)=   0.074688300490379d0
      basgrd(3006)=   0.906052291393280d0
      basgrd(3007)=  -0.273850411176682d0
      basgrd(3008)=  -0.071015596389771d0
      basgrd(3009)=   0.959146916866303d0
      basgrd(3010)=  -0.139922693371773d0
      basgrd(3011)=  -0.233180001378059d0
      basgrd(3012)=   0.962314188480377d0
      basgrd(3013)=  -0.460972994565964d0
      basgrd(3014)=  -0.107237502932549d0
      basgrd(3015)=   0.880910873413086d0
      basgrd(3016)=  -0.317870497703552d0
      basgrd(3017)=  -0.264119803905487d0
      basgrd(3018)=   0.910603702068329d0
      basgrd(3019)=  -0.480463594198227d0
      basgrd(3020)=  -0.280173301696777d0
      basgrd(3021)=   0.831058084964752d0
      basgrd(3022)=   0.084543600678444d0
      basgrd(3023)=  -0.996047914028168d0
      basgrd(3024)=   0.027221400290728d0
      basgrd(3025)=  -0.072493098676205d0
      basgrd(3026)=  -0.995373606681824d0
      basgrd(3027)=  -0.063057199120522d0
      basgrd(3028)=  -0.230424597859383d0
      basgrd(3029)=  -0.959924280643463d0
      basgrd(3030)=  -0.159530192613602d0
      basgrd(3031)=  -0.374559789896011d0
      basgrd(3032)=  -0.890224516391754d0
      basgrd(3033)=  -0.259239792823792d0
      basgrd(3034)=  -0.496119290590286d0
      basgrd(3035)=  -0.793715596199036d0
      basgrd(3036)=  -0.351967602968216d0
      basgrd(3037)=   0.122028797864914d0
      basgrd(3038)=  -0.981126725673676d0
      basgrd(3039)=  -0.149997904896736d0
      basgrd(3040)=  -0.047386601567268d0
      basgrd(3041)=  -0.963364779949188d0
      basgrd(3042)=  -0.263974994421005d0
      basgrd(3043)=  -0.197766304016113d0
      basgrd(3044)=  -0.911159574985504d0
      basgrd(3045)=  -0.361492395401001d0
      basgrd(3046)=  -0.350990086793900d0
      basgrd(3047)=  -0.816504597663879d0
      basgrd(3048)=  -0.458395212888718d0
      basgrd(3049)=   0.150445401668549d0
      basgrd(3050)=  -0.931985616683960d0
      basgrd(3051)=  -0.329801589250565d0
      basgrd(3052)=  -0.027312800288200d0
      basgrd(3053)=  -0.898584783077240d0
      basgrd(3054)=  -0.437949001789093d0
      basgrd(3055)=  -0.185130193829536d0
      basgrd(3056)=  -0.815301477909088d0
      basgrd(3057)=  -0.548644006252289d0
      basgrd(3058)=   0.163075402379036d0
      basgrd(3059)=  -0.850921928882599d0
      basgrd(3060)=  -0.499337911605835d0
      basgrd(3061)=  -0.008667900227010d0
      basgrd(3062)=  -0.791437387466431d0
      basgrd(3063)=  -0.611188709735870d0
      basgrd(3064)=   0.163418292999268d0
      basgrd(3065)=  -0.745187819004059d0
      basgrd(3066)=  -0.646521091461182d0
      basgrd(3067)=   0.265319794416428d0
      basgrd(3068)=  -0.409477591514587d0
      basgrd(3069)=  -0.872888028621674d0
      basgrd(3070)=   0.143090292811394d0
      basgrd(3071)=  -0.295900613069534d0
      basgrd(3072)=  -0.944440603256226d0
      basgrd(3073)=   0.015742199495435d0
      basgrd(3074)=  -0.164145201444626d0
      basgrd(3075)=  -0.986310601234436d0
      basgrd(3076)=  -0.105124600231648d0
      basgrd(3077)=  -0.019837500527501d0
      basgrd(3078)=  -0.994261205196381d0
      basgrd(3079)=  -0.212336897850037d0
      basgrd(3080)=   0.123466201126576d0
      basgrd(3081)=  -0.969365298748016d0
      basgrd(3082)=   0.347746193408966d0
      basgrd(3083)=  -0.250418514013290d0
      basgrd(3084)=  -0.903528213500977d0
      basgrd(3085)=   0.218155995011330d0
      basgrd(3086)=  -0.106680497527123d0
      basgrd(3087)=  -0.970065593719482d0
      basgrd(3088)=   0.094297297298908d0
      basgrd(3089)=   0.030824799090624d0
      basgrd(3090)=  -0.995066821575165d0
      basgrd(3091)=  -0.041170600801706d0
      basgrd(3092)=   0.183561593294144d0
      basgrd(3093)=  -0.982145726680756d0
      basgrd(3094)=   0.416561812162399d0
      basgrd(3095)=  -0.074573896825314d0
      basgrd(3096)=  -0.906043589115143d0
      basgrd(3097)=   0.273732393980026d0
      basgrd(3098)=   0.070935599505901d0
      basgrd(3099)=  -0.959186494350433d0
      basgrd(3100)=   0.139430493116379d0
      basgrd(3101)=   0.233312293887138d0
      basgrd(3102)=  -0.962353587150574d0
      basgrd(3103)=   0.460455000400543d0
      basgrd(3104)=   0.107335403561592d0
      basgrd(3105)=  -0.881169795989990d0
      basgrd(3106)=   0.318075001239777d0
      basgrd(3107)=   0.264266103506088d0
      basgrd(3108)=  -0.910489797592163d0
      basgrd(3109)=   0.480827689170837d0
      basgrd(3110)=   0.279862493276596d0
      basgrd(3111)=  -0.830952286720276d0
      basgrd(3112)=   0.573881328105927d0
      basgrd(3113)=   0.590416193008423d0
      basgrd(3114)=  -0.567511320114136d0
      basgrd(3115)=   0.457060813903809d0
      basgrd(3116)=   0.721543192863464d0
      basgrd(3117)=  -0.520068109035492d0
      basgrd(3118)=   0.324556797742844d0
      basgrd(3119)=   0.832555472850800d0
      basgrd(3120)=  -0.448903411626816d0
      basgrd(3121)=   0.185008004307747d0
      basgrd(3122)=   0.916637718677521d0
      basgrd(3123)=  -0.354326695203781d0
      basgrd(3124)=   0.049186598509550d0
      basgrd(3125)=   0.967816472053528d0
      basgrd(3126)=  -0.246803507208824d0
      basgrd(3127)=   0.629750072956085d0
      basgrd(3128)=   0.660424470901489d0
      basgrd(3129)=  -0.408967405557632d0
      basgrd(3130)=   0.498823106288910d0
      basgrd(3131)=   0.799293577671051d0
      basgrd(3132)=  -0.335119694471359d0
      basgrd(3133)=   0.362880885601044d0
      basgrd(3134)=   0.897072374820709d0
      basgrd(3135)=  -0.252148002386093d0
      basgrd(3136)=   0.201453804969788d0
      basgrd(3137)=   0.968284785747528d0
      basgrd(3138)=  -0.147786498069763d0
      basgrd(3139)=   0.661239683628082d0
      basgrd(3140)=   0.713880121707916d0
      basgrd(3141)=  -0.230515196919441d0
      basgrd(3142)=   0.513390004634857d0
      basgrd(3143)=   0.844129681587219d0
      basgrd(3144)=  -0.154517501592636d0
      basgrd(3145)=   0.355951100587845d0
      basgrd(3146)=   0.933367311954498d0
      basgrd(3147)=  -0.046090599149466d0
      basgrd(3148)=   0.659547328948975d0
      basgrd(3149)=   0.750302910804749d0
      basgrd(3150)=  -0.045198898762465d0
      basgrd(3151)=   0.503486573696137d0
      basgrd(3152)=   0.862651407718658d0
      basgrd(3153)=   0.048309899866581d0
      basgrd(3154)=   0.630965173244476d0
      basgrd(3155)=   0.764294981956482d0
      basgrd(3156)=   0.133176699280739d0
      basgrd(3157)=   0.583324909210205d0
      basgrd(3158)=   0.621047675609589d0
      basgrd(3159)=   0.523480474948883d0
      basgrd(3160)=   0.434445410966873d0
      basgrd(3161)=   0.649718582630158d0
      basgrd(3162)=   0.623797178268433d0
      basgrd(3163)=   0.268585890531540d0
      basgrd(3164)=   0.652529895305634d0
      basgrd(3165)=   0.708566427230835d0
      basgrd(3166)=   0.094523601233959d0
      basgrd(3167)=   0.624814391136169d0
      basgrd(3168)=   0.775030493736267d0
      basgrd(3169)=  -0.072967201471329d0
      basgrd(3170)=   0.572320103645325d0
      basgrd(3171)=   0.816777527332306d0
      basgrd(3172)=   0.577592790126801d0
      basgrd(3173)=   0.491998493671417d0
      basgrd(3174)=   0.651401579380035d0
      basgrd(3175)=   0.406408011913300d0
      basgrd(3176)=   0.502621114253998d0
      basgrd(3177)=   0.763023316860199d0
      basgrd(3178)=   0.236395999789238d0
      basgrd(3179)=   0.489679813385010d0
      basgrd(3180)=   0.839244127273560d0
      basgrd(3181)=   0.042149901390076d0
      basgrd(3182)=   0.452822685241699d0
      basgrd(3183)=   0.890603721141815d0
      basgrd(3184)=   0.546575188636780d0
      basgrd(3185)=   0.343986690044403d0
      basgrd(3186)=   0.763497591018677d0
      basgrd(3187)=   0.360686004161835d0
      basgrd(3188)=   0.352005213499069d0
      basgrd(3189)=   0.863711714744568d0
      basgrd(3190)=   0.165439903736115d0
      basgrd(3191)=   0.316892594099045d0
      basgrd(3192)=   0.933921217918396d0
      basgrd(3193)=   0.486045807600021d0
      basgrd(3194)=   0.189849093556404d0
      basgrd(3195)=   0.853063225746155d0
      basgrd(3196)=   0.291057795286179d0
      basgrd(3197)=   0.176050901412964d0
      basgrd(3198)=   0.940367698669434d0
      basgrd(3199)=   0.406175196170807d0
      basgrd(3200)=   0.038210701197386d0
      basgrd(3201)=   0.912995994091034d0
      basgrd(3202)=   0.280388712882996d0
      basgrd(3203)=  -0.360132902860642d0
      basgrd(3204)=   0.889767587184906d0
      basgrd(3205)=   0.107192799448967d0
      basgrd(3206)=  -0.411204308271408d0
      basgrd(3207)=   0.905218601226807d0
      basgrd(3208)=  -0.073934599757195d0
      basgrd(3209)=  -0.454811185598373d0
      basgrd(3210)=   0.887513577938080d0
      basgrd(3211)=  -0.251210808753967d0
      basgrd(3212)=  -0.491509795188904d0
      basgrd(3213)=   0.833853185176849d0
      basgrd(3214)=  -0.409631103277206d0
      basgrd(3215)=  -0.516359388828278d0
      basgrd(3216)=   0.752047419548035d0
      basgrd(3217)=   0.263756006956101d0
      basgrd(3218)=  -0.522789120674133d0
      basgrd(3219)=   0.810631990432739d0
      basgrd(3220)=   0.068852201104164d0
      basgrd(3221)=  -0.586667299270630d0
      basgrd(3222)=   0.806895792484283d0
      basgrd(3223)=  -0.109541699290276d0
      basgrd(3224)=  -0.627157628536224d0
      basgrd(3225)=   0.771151006221771d0
      basgrd(3226)=  -0.299137204885483d0
      basgrd(3227)=  -0.650056481361389d0
      basgrd(3228)=   0.698529481887817d0
      basgrd(3229)=   0.230894193053246d0
      basgrd(3230)=  -0.673290491104126d0
      basgrd(3231)=   0.702401399612427d0
      basgrd(3232)=   0.026683200150728d0
      basgrd(3233)=  -0.724901497364044d0
      basgrd(3234)=   0.688335478305817d0
      basgrd(3235)=  -0.169200897216797d0
      basgrd(3236)=  -0.763878107070923d0
      basgrd(3237)=   0.622785091400147d0
      basgrd(3238)=   0.178987994790077d0
      basgrd(3239)=  -0.799801588058472d0
      basgrd(3240)=   0.572957873344421d0
      basgrd(3241)=  -0.025643700733781d0
      basgrd(3242)=  -0.846206605434418d0
      basgrd(3243)=   0.532237589359283d0
      basgrd(3244)=   0.117036201059818d0
      basgrd(3245)=  -0.894675195217133d0
      basgrd(3246)=   0.431113392114639d0
      return
      end
c.... search and return a free file unit

      subroutine get_file(ifc,mes)

c     mes : error mesage
c     ifc : file unit or -1

      character (len=100) :: mes
      logical           isda
      ifc=0
c
      do i=10,99
          inquire(i,opened=isda)
          if (.not.isda) then
              ifc = i
              goto 10
          endif
      enddo
      ifc = -1
      mes = 'COSMO: get_file: cannot find free filehandle'
 10   return
      end
c----------------------------------------------------------------------c
      subroutine setamat(coord,natoms,iatsp,dirvec,
     &              nar, ar, a1mat, nsetf, nset, xsp, disex2, rsolv,
     &              routf, nppa, nps, npspher, srad, tm,
     &              ierr, mes, npsd, a23mat,jobname,lenj)
c----------------------------------------------------------------------c
c     set up the A matrices and perform first step of choleski for a1mat
c
c     xsp:    cosurf from consts
c     coord:  must be the distorted coordinates from consts
c     calls:
c           get_file (ok)
c           coschol1 (ok)
c----------------------------------------------------------------------c
      implicit real*8 (a-h,o-z)

      logical  newvers
      character*100 mes
      character*256 jobname

      dimension    coord(3,natoms)
      dimension    iatsp(2*nps)
      dimension    dirvec(3,nppa)
      dimension    nar(nps+npspher), ar(nps+npspher)
      dimension    xsp(3,nps+npspher), nsetf(nps+npspher)
      dimension    nset(natoms*nppa), a1mat(nps*(nps+1)/2)
      dimension    xx(3),xa(3),xi(3),xj(3)
      dimension    tm(3,3,natoms), srad(natoms)

c     double precision, allocatable :: a23mat(:)
c
c     a23mat is now a parameter
c
      dimension a23mat(npsd)

      parameter(pi=3.14159265358979323844d0)
      fdiagr=2.1D0*sqrt(pi)
      sqpi=sqrt(pi)

c...  switch on the newest amat version here
c      newvers=.false.
      newvers=.true.
c...  total number of seg. (outer+inner)
c     npsd=nps+npspher
c---------------------------------------------------------------------c
c     allocate(a23mat(npsd),stat=ierr)
c     if (ierr .ne. 0) then
c         mes=' COSMO: Cannot alloc. a23mat in setamat'
c         ierr = -1
c         return
c     endif

c----------------------------------------------------------------------
c set A-1

      fdiag=1.05d0*sqrt(nppa+0.d0)
c
c ---------------------------------------------------------------------
c     the A matrix is now ready to be set
c ---------------------------------------------------------------------

      call get_file(ifca2,mes)
      if(ifca2.eq.-1) then
        ierr=-1
        return
      endif
      open(ifca2,file=jobname(1:lenj)//'.a2mat',form='unformatted',
     &     status='unknown',err=910)

      call get_file(ifca3,mes)
      if(ifca3.eq.-1) then
        ierr=-1
        return
      endif
      open(ifca3,file=jobname(1:lenj)//'.a3mat',form='unformatted',
     &     status='unknown',err=911)

      do 450 ips=1,npsd
c skip torus
        if(ips.gt.npspher.and.ips.le.nps) then
          aa=fdiagr/sqrt(ar(ips))
        else

          i=iatsp(ips)
          ri=srad(i)-rsolv

c         if (ips .gt. nps) ri=ri+rsolv-0.05d0*rsolv
          if (ips .gt. nps) ri=ri+routf*rsolv
          nari=nar(ips)
c         write(6,*) ' segment ',ips,' # points ',nari
          nsetfi=nsetf(ips)
          aa=0.d0
          do 350 k=nsetfi,nsetfi+nari-1
            j1=nset(k)
            aa=aa+fdiag
            x1=dirvec(1,j1)
            x2=dirvec(2,j1)
            x3=dirvec(3,j1)
            do 350 l=nsetfi,k-1
              j2=nset(l)
              aa=aa+2.d0/sqrt((x1-dirvec(1,j2))**2+
     1           (x2-dirvec(2,j2))**2+(x3-dirvec(3,j2))**2)
  350     continue
          aa=aa/ri/nari/nari
        endif

        if (ips .le. nps) then
          ia1=ips*(ips+1)/2
          a1mat(ia1)=aa
        endif
        a23mat(ips)=aa

        do 360 ix=1,3
          xi(ix)=coord(ix,i)
  360   xa(ix)=xsp(ix,ips)

        do 440 jps=ips+1,npsd

          narj=nar(jps)
          nsetfj=nsetf(jps)
          j=iatsp(jps)
          dist=0.d0
          do 370 ix=1,3
            xj(ix)=coord(ix,j)-xi(ix)
  370       dist=dist+(xsp(ix,jps)-xa(ix))**2
          if (dist .lt. disex2 .and.
     &      (ips.le.npspher.or.ips.gt.nps) .and.
     &      (jps.le.npspher.or.jps.gt.nps)) then
            rj=srad(j)-rsolv
c           if(jps .gt. nps) rj=rj+rsolv-0.05d0*rsolv
            if(jps .gt. nps) rj=rj+routf*rsolv

            aij=0.d0
            do 430 k=nsetfi,nsetfi+nari-1
              j1=nset(k)
              do 380 ix=1,3
  380           xx(ix)=dirvec(ix,j1)*ri
              if (i .ne. j) then
                x1=xx(1)*tm(1,1,i)+xx(2)*tm(2,1,i)+xx(3)*tm(3,1,i)-xj(1)
                x2=xx(1)*tm(1,2,i)+xx(2)*tm(2,2,i)+xx(3)*tm(3,2,i)-xj(2)
                x3=xx(1)*tm(1,3,i)+xx(2)*tm(2,3,i)+xx(3)*tm(3,3,i)-xj(3)
                do 400 l=nsetfj,nsetfj+narj-1
                  j2=nset(l)
                  do 390 ix=1,3
  390               xx(ix)=dirvec(ix,j2)*rj
                  y1=xx(1)*tm(1,1,j)+xx(2)*tm(2,1,j)+xx(3)*tm(3,1,j)-x1
                  y2=xx(1)*tm(1,2,j)+xx(2)*tm(2,2,j)+xx(3)*tm(3,2,j)-x2
                  y3=xx(1)*tm(1,3,j)+xx(2)*tm(2,3,j)+xx(3)*tm(3,3,j)-x3
                  aij=aij+1.d0/sqrt(y1*y1+y2*y2+y3*y3)
  400           continue
              else
                do 420 l=nsetfj,nsetfj+narj-1
                  j2=nset(l)
                  aa=(dirvec(1,j2)*rj-xx(1))**2+(dirvec(2,j2)*rj
     &               -xx(2))**2+(dirvec(3,j2)*rj-xx(3))**2
                  aij = aij + 1.0d0/sqrt(aa)
c                 aij=aij+((dirvec(1,j2)*rj-xx(1))**2+(dirvec(2,j2)*rj
c    &                -xx(2))**2+(dirvec(3,j2)*rj-xx(3))**2)**(-.5d0)
  420           continue
              end if

  430       continue
            aij=aij/nari/narj

c... start of new part
c
c ips is located on the inner spheres  amd jps is a trp or a ring segment
c
          else if (dist .lt. disex2 .and.
     &      (ips.le.npspher) .and. newvers .and.
     &      (jps.gt.npspher.and.jps.le.nps)) then
            aij=0.d0

c calc. radius of basis segment rk sqrt(area/pi)
c area of basis segmnet: ri**2*pi*4/nppa

            rk = sqrt(4.d0*ri**2/nppa)
            do k=nsetfi,nsetfi+nari-1
              j1=nset(k)
              do ix=1,3
                xx(ix)=dirvec(ix,j1)*ri
              end do
              x1=xx(1)*tm(1,1,i)+xx(2)*tm(2,1,i)+xx(3)*tm(3,1,i)+xi(1)
              x2=xx(1)*tm(1,2,i)+xx(2)*tm(2,2,i)+xx(3)*tm(3,2,i)+xi(2)
              x3=xx(1)*tm(1,3,i)+xx(2)*tm(2,3,i)+xx(3)*tm(3,3,i)+xi(3)
              y1=xsp(1,jps)-x1
              y2=xsp(2,jps)-x2
              y3=xsp(3,jps)-x3
c use the same damping function as below to handle to close trp/ring
c basis grid point distances
              bdist = sqrt(y1*y1+y2*y2+y3*y3)
              rskj=(rk+sqrt(ar(jps)/pi))/2.d0
              if (bdist.gt.rskj) then
                 aij=aij+1.d0/bdist
              else
                 ak00=fdiagr/sqpi/rskj
                 aij=aij+ak00+(bdist/rskj)*(1.d0/rskj-ak00)
              endif
            end do
            aij=aij/nari
c...    end of new part

          else

c ak midi  aug. 2001
c now the center appr. is used, with damping for close contacts
c For all pairs with a distance < 0.5 sum of segment radius a damping
c function is used.
c (0.5 seems to be a good approximation)

            if (dist .ge. disex2) then
               aij=1.d0/sqrt(dist)
            else
               rsij=(sqrt(ar(ips)/pi)+sqrt(ar(jps)/pi))/2.d0
               dist=sqrt(dist)
               if (dist.gt.rsij) then
                 aij=1.d0/dist
               else
                 a00=fdiagr/sqpi/rsij
                 aij=a00+(dist/rsij)*(1.d0/rsij-a00)
               end if
             end if
          end if

          if (ips .le. nps .and. jps .le. nps) then
            ia1=ia1+jps-1
            a1mat(ia1)=aij
          end if
          a23mat(jps)=aij

  440   continue
        if (ips.le.nps) then
c         write row of full A2 matrix
          write (ifca2) (a23mat(jps),jps=nps+1,npsd)
        else
c         write row of triangular A3 matrix
          write (ifca3) (a23mat(jps),jps=ips,npsd)
        endif

  450 continue
      close(ifca2)
      close(ifca3)

c     write lower triangle of A1 matrix
      call get_file(ifca,mes)
      if(ifca.eq.-1) then
        ierr=-1
        return
      endif
c     lamat = index(amat,' ')-1

      open(ifca,file=jobname(1:lenj)//'.amat',form='unformatted',
     &     status='unknown',err=900)

c     write (ifca) (a1mat(i),i=1,nps*(nps+1)/2)
      write (ifca) a1mat
      close(ifca)
cas   Cholesky decomposition of A1
      call coschol1(a1mat,nps,info)
      if (info.lt.0) then
         ierr=-1
         mes='COSMO problem: a1mat not positive definite in consts.f'
         return
      endif
c
c-----------------------------------------------------------
c     deallocate(a23mat, STAT=ierr)
c     if (ierr /=0 ) then
c        ierr=-1
c        mes='COSMO setamat: cannot deallocate mem.'
c        return
c     endif
      return
c
c-----------------------------------------------------------
  900 mes ='COSMO setamat: could not open amat file'
      ierr = -1
      return
  910 mes='subroutine setamat: could not open file a2mat.tmp'
      ierr = -1
      return
  911 mes='subroutine setamat: could not open file a3mat.tmp'
      ierr = -1
      return
      end
      subroutine writsude(nip,sude,isude,ierr,mes)
c
c    MM Dec 2003  now name of sude file is taken from depository,
c                 and the number of intersection pairs (nip) is
c                 also stored into depository
c                 I have also eliminated the implicit loops

      implicit real*8 (a-h,o-z)

      character*100 mes
      character*256 csude
      dimension sude(4*nip),isude(4*nip)
      ierr=0
c
c -- write file only if there is something to store   ! MM Aug 2008
      If(nip.GT.0) Then
        call getchval('c_sude',csude)
        call get_file(ifc,mes)
        if(ifc.eq.-1) then
          ierr=-1
          return
        endif
        open(ifc,file=csude(1:len_trim(csude)),status='unknown',err=900)
        write(ifc,*) nip
        call setival('c_nip',nip)
        write(ifc,*) sude
        write(ifc,*) isude
        close(ifc)
      EndIf
      return

  900 mes='COSMO: writsude  cannot open file sude'
      ierr=-1
      return
      end
      subroutine reasude(nip,sude,isude,ierr,mes)
c
c    MM Dec 2003
c

      implicit real*8 (a-h,o-z)

      character*100 mes
      character*256 csude
      dimension sude(4*nip),isude(4*nip)
      ierr=0
c
c -- read file only if there is something to store   ! MM Aug 2008
      If(nip.GT.0) Then
        call getchval('c_sude',csude)
        call get_file(ifc,mes)
        if(ifc.eq.-1) then
          ierr=-1
          return
        endif
        open(ifc,file=csude(1:len_trim(csude)),status='old',err=900)
        read(ifc,*)   ! skip this record, as nip is read from depository
        read(ifc,*) sude
        read(ifc,*) isude
        close(ifc)
      EndIf
      return

  900 mes='COSMO: reasude  cannot open file sude'
      ierr=-1
      return
      end
      subroutine  cosmo_write(bl, icharge, iatsp ,etot, xyz,
     $                        ierr, cerm, natoms)
      implicit real*8 (a-h,o-z)
c
c   print final output to file
c
c   this is an adaptation of the Fortran 90 module
c
c   bl  real memory segment
c   icharge : atomic changes
c   iatsp: array assigning surface segments to atoms
c   etot : scf energy (corrected by ediel)
c   xyz : atomic coordinates without dummies
c   ierr : /= 0 if error
c   cerm : error message
c   natoms : number of atoms
c
      character*100 cerm
      character*256 outfile,info
      dimension bl(*),icharge(natoms),iatsp(*)
      CHARACTER*2 symbol_lc, symbol_uc
      CHARACTER*5 ctmp
      CHARACTER*6 cavity
      character*8 atsy
      parameter (idim=200)
      dimension icnt(idim)
      INTEGER  at
      DIMENSION xyz(3,natoms)
c
      character*26 alpbet,betalp
      data alpbet / 'abcdefghijklmnopqrstuvwxyz' /
      data betalp / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /

      !
      call getrval('angs',angs)
      toang=1.0d0/angs
      ierr = 0
      !
      !... open cosmo file
      !
      call getchval('c_outfil',outfile)
      CALL get_file(ifco,cerm)
      IF(ifco .eq. -1) then
        ierr=-1
        RETURN
      endif
      OPEN(ifco,FILE=outfile(1:len_trim(outfile)),
     $     FORM='formatted',STATUS='unknown',ERR=1000)
      !
      !... print input
      !
      call getchval('c_info',info)
      write(ifco,'(a)') '$info'
      write(ifco,'(a)') info
      write(ifco,'(a)') '$cosmo'

      call getrval('c_fepsi',fepsi)
      call getrval('c_eps',eps)
      IF (fepsi .eq. 1.0d0) THEN
         write(ifco,'(a)') '  epsilon=infinity'
      ELSE
         write(ifco,'(a,f5.2)') '  epsilon=', eps
      END IF
      call getival('c_nppa',nppa)
      write(ifco,'(a,i5)') '  nppa=', nppa
      call getival('c_nspa',nspa)
      write(ifco,'(a,i5)') '  nspa=', nspa
      call getrval('c_disex',disex)
      write(ifco,'(a,g12.6)') '  disex=', disex
      call getrval('c_rsolv',rsolv)
      write(ifco,'(a,f5.2)') '  rsolv=', rsolv*toang
      call getrval('c_routf',routf)
      write(ifco,'(a,f5.2)') '  routf=', routf
      call getival('c_lcavit',lcavity)
      IF (lcavity .eq. 0) THEN
          cavity = 'open'
      ELSE
          cavity = 'closed'
      END IF
      write(ifco,'(a,a6)') '  cavity ', cavity
      call getrval('c_phsran',phsran)
      write(ifco,'(a,g9.2)') '  phsran=',phsran
      call getrval('c_ampran',ampran)
      write(ifco,'(a,g9.2)') '  ampran=',ampran
      !
      !... write calculated parameters and variables $cosmo_data
      !
      write(ifco,'(a)')'$cosmo_data'
      write(ifco,'(a,f9.6)') '  fepsi=',fepsi
      call getrval('c_disex2',disex2)
      write(ifco,'(a,g12.6)')'  disex2=',disex2
      call getival('c_maxnps',maxnps)
      write(ifco,'(a,i5)')   '  maxnps=',maxnps
      call getival('c_nsph',nsph)
      write(ifco,'(a,i5)')   '  nsph=',nsph
      call getival('c_nps',nps)
      write(ifco,'(a,i5)')   '  nps=',nps
      call getival('c_npsphe',npspher)
      write(ifco,'(a,i5)')   '  npsd=',nps+npspher
      write(ifco,'(a,i5)')   '  npspher=',npspher
      call getrval('c_area',area)
      write(ifco,'(a,f8.2)') '  area=',area
      call getrval('c_volume',volume)
      write(ifco,'(a,f8.2)') '  volume=',volume
      !
      !... print final geometry (coordinates and radii)
      !
      call getival('c_israd',israd)
      write(ifco,'(a)') '$coord_rad'
      write(ifco,'(a,a)') '#atom   x                  y          ',
     $                   '        z             element  radius [A]'
      DO i=1,natoms
         sradi=bl(israd+i-1)
         call getatsym1(i,natoms,atsy)
         write(ifco,10) i,(xyz(j,i),j=1,3),atsy,
     $               sradi*toang
      END DO
 10   FORMAT(1x,i3,3(1x,f18.14),2x,a2,2x,f10.5)
      !... for use in cosmors print coordinates in .car format (archive 3)
      write(ifco,'(a,/,a,/,a,/,a,/,a)')'$coord_car','!BIOSYM archive 3',
     $ 'PBC=OFF','coordinates from COSMO calculation','!DATE '

      call izeroit(icnt,idim)
      DO i=1,natoms
         call getatsym1(i,natoms,atsy)
         symbol_lc=atsy(1:2)
         symbol_uc=symbol_lc
         ii=index(alpbet,symbol_lc(2:2))
         if(ii.gt.0)symbol_uc(2:2)=betalp(ii:ii)
         ii=index(alpbet,symbol_lc(1:1))
         if(ii.gt.0)symbol_uc(1:1)=betalp(ii:ii)
         j=icharge(i)
         IF (j .gt. idim) THEN
            cerm='cosmo_write: unreasonable nuc. charge'
            ierr = -1
            RETURN
         END IF
         icnt(j)=icnt(j)+1
         IF(symbol_lc(2:2) == ' ') THEN
            write(ctmp(1:5),'(a,i1)') symbol_uc(1:1),icnt(j)
         ELSE
            write(ctmp(1:5),'(a,i1)') symbol_uc(1:2),icnt(j)
         END IF
         write(ifco,20) ctmp,xyz(1,i)*toang,xyz(2,i)*toang,
     $                  xyz(3,i)*toang,symbol_lc,symbol_uc,0.d0
      END DO
 20   FORMAT(a5,3f15.9,' COSM 1      ',a,6x,a,f7.3)
      write(ifco,'(a/a)') 'end','end'
      !
      !... write charges
      !
        call getrval('c_de',de)
        call getrval('c_dq',dq)
        call getrval('c_qsum',qsum)
        call getrval('c_ediel',ediel)
      write(ifco,'(a)')'$screening_charge'
      write(ifco,'(a,f10.6)') '  cosmo      = ',qsum
      write(ifco,'(a,f10.6)') '  correction = ',dq
      write(ifco,'(a,f10.6)') '  total      = ',qsum+dq
      !
      !... write energies
      !
      write(ifco,30) etot,etot+de,etot+0.5d0*de,
     $               ediel,ediel+de
 30   FORMAT('$cosmo_energy',/,
     $  '  Total energy [a.u.]            =   ', f17.10,/,
     $  '  Total energy + OC corr. [a.u.] =   ', f17.10,/,
     $  '  Total energy corrected [a.u.]  =   ', f17.10,
     $  ' Note: incorrect value',
     $  ' contained for downward compatibility',/,
     $  '  Dielectric energy [a.u.]       =   ', f17.10,/,
     $  '  Diel. energy + OC corr. [a.u.] =   ', f17.10)
      !
      !... write segment information
      !
      write(ifco,40)
 40   FORMAT('$segment_information',/,
     $'# n             - segment number',/,
     $'# atom          - atom associated with segment n',/,
     $'# position      - segment coordinates [a.u.]',/,
     $'# charge        - segment charge (corrected)',/,
     $'# area          - segment area [A**2]',/,
     $'# potential     - solute potential on segment (A length scale)',
     $/,
     $'#',/,
     $'#  n   atom',14x,'position (X, Y, Z)',19x,
     $'charge         area        charge/area     potential',
     $/,'#',/,'#')
      call getival('c_nps',nps)
      call getival('c_iqcosc',iqcosc)
      call getival('c_iphic',iphic)
      call getival('c_iar',iar)
      call getival('c_icosur',icosurf)
      DO i=1,nps
        at = iatsp(i)
        q = bl(iqcosc+i-1)
        p = bl(iphic+i-1)/toang
        a = bl(iar+i-1)*toang*toang
        cx=bl(icosurf+(i-1)*3)
        cy=bl(icosurf+(i-1)*3+1)
        cz=bl(icosurf+(i-1)*3+2)
        write(ifco,'(i5,i5,7f15.9)') i,at,cx,cy,cz,q,a,q/a,p
      END DO

      CLOSE(ifco)

      RETURN
 1000    cerm='Cannot open '//outfile(1:len_trim(outfile))
         ierr = -1
      RETURN
      END
      subroutine cavnucgrd(natoms,xyz,nps,cosurf,iatsp,qcos,ar,
     &                  nip,sude,isude,lcavity,gcos,icharge)

      implicit real*8 ( a-h, o-z )

c
c  MM Dec 2003: changed the charge array from real to integer,
c               to be consistent with the other routines.
c               Use subroutine vscal to change sign
c.........................................................................
c  calculation of:
c               -  gradient contribution of the A matrix
c               -  the surface derivative term
c               -  controbution from the nuclear potential
c  natoms:    number of atoms
c  xyz:       atomic coordinates
c  nps:       number of segments (only inner surface)
c  cosurf:    segment coordinates
c  iatsp:     mapping segment -> atom
c  qcos:      screening charges
c  ar:        area of segments
c  nips:      dimension for sude and isude
c  sude, isude: array for surface derivative term
c  lcavity:     open or closed cavity/surface
c  gcos:        gradient array
c  icharge:     nuclear charges
c.........................................................................

      parameter(pi=3.141592653589793238d0)
      parameter(dmin=1.d-5)

      dimension   gcos(3,natoms)
      dimension   cosurf(3,nps), qcos(nps), iatsp(nps), ar(nps)
      dimension   xk(3),xl(3)
      dimension   xyz(3,natoms), icharge(natoms)
      dimension   sude(4,nip), isude(4,nip)

      fdiagr=2.1d0*sqrt(pi)
c
c....  contribution of A matrix (q*grad(A)*q)
c
      do k=1,nps
         iak=iatsp(k)
         do ix=1,3
            xk(ix)=cosurf(ix,k)
         enddo
         qsk=qcos(k)
         do l=k+1,nps
            ial=iatsp(l)
            if (ial.ne.iak) then
               dist2=0.d0
               do ix=1,3
                  xxx=cosurf(ix,l)-xk(ix)
                  xl(ix)=xxx
                  dist2=dist2+xxx*xxx
               enddo
               ff=qsk*qcos(l)*dist2**(-1.5d0)
               do ix=1,3
                  gcos(ix,iak)=gcos(ix,iak)-xl(ix)*ff
                  gcos(ix,ial)=gcos(ix,ial)+xl(ix)*ff
               enddo
            endif
         enddo
      enddo
c
c...  change of surface with atom movement
c
      if (lcavity .eq. 1) then
         do i=1,nip
            ips0 = isude(4,i)
            ia=isude(1,i)
            ib=isude(2,i)
            suma=0.d0
            sumb=0.d0
            nsab=isude(3,i)
            do ips=ips0+1,ips0+nsab
               iat=iatsp(ips)
               if(iat .eq. ia) then
                  suma=suma+qcos(ips)**2/sqrt(ar(ips))
               else
                  sumb=sumb+qcos(ips)**2/sqrt(ar(ips))
               endif
            enddo
            deab=-0.5d0*fdiagr*(suma*sude(3,i)/sude(1,i)
     &                     +sumb*sude(4,i)/sude(2,i) )**2
            xk(1)=xyz(1,ib)-xyz(1,ia)
            xk(2)=xyz(2,ib)-xyz(2,ia)
            xk(3)=xyz(3,ib)-xyz(3,ia)
            deab=deab/sqrt(xk(1)**2+xk(2)**2+xk(3)**2)
            do ix=1,3
               gcos(ix,ia)=gcos(ix,ia)-xk(ix)*deab
               gcos(ix,ib)=gcos(ix,ib)+xk(ix)*deab
            enddo
         enddo
      endif
c    the upper contribution is -1/2 q*grad(A)*q
      call vscal(3*natoms,-1.d0,gcos)
c
c...  contribution from nuclear potential
c
      do k=1,nps
         kat=iatsp(k)
         do n=1,natoms
            if (n.ne.kat) then
               dx = cosurf(1,k) - xyz(1,n)
               dy = cosurf(2,k) - xyz(2,n)
               dz = cosurf(3,k) - xyz(3,n)
               dist = dx*dx + dy*dy + dz*dz
               dd = sqrt(dist)
               if (dd.gt.dmin) then
                  dd = qcos(k)/(dd*dist)
                  dx = dx*dd
                  dy = dy*dd
                  dz = dz*dd
                  gcos(1,n) = gcos(1,n) + float(icharge(n))*dx
                  gcos(2,n) = gcos(2,n) + float(icharge(n))*dy
                  gcos(3,n) = gcos(3,n) + float(icharge(n))*dz
                  gcos(1,kat) = gcos(1,kat) - float(icharge(n))*dx
                  gcos(2,kat) = gcos(2,kat) - float(icharge(n))*dy
                  gcos(3,kat) = gcos(3,kat) - float(icharge(n))*dz
               endif
            endif
         enddo
      enddo

      return
      end

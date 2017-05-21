c ==================================================================
c  SYMMETRY ROUTINES            JB   August 1997
c ==================================================================
c
      subroutine getsym(xyz,    xyzsh,  rtoc,   rtocsh, trans,
     $                  univec, labsh,  npersh, istype, invt,
     $                  ipath,  jpath,  n,      nt,     logish,
     $                  mtype,  mtypsh, sflies, ntrans, eps,
     $                  adjst,  ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     assign Schoenflies symbol to given set of coordinates
c     adjust coordinates if necessary to standard orientation
c =====================================================================
      dimension xyz(3,n),xyzsh(3,n),rtoc(n),rtocsh(n),
     1          labsh(n),npersh(n)
      dimension trans(9,ntrans),univec(3,ntrans),gen(9)
      dimension istype(ntrans),invt(ntrans)
      dimension ipath(ntrans),jpath(ntrans)
      dimension e1(3),e2(3),e3(3)
      character*8 mtype(n),mtypsh(n)
      character*4 sflies
      logical logish(n),adjst
c
      dimension aline(3),aplane(3)   ! least-squares line and plane
      parameter (zero=0.0d0,one=1.0d0)
c
      sflies='?'//'?'//'?'//'?'      ! to avoid error in IBM preprocessor
c ---------------------------------------------------------------------
c     prepare data for group construction
c ---------------------------------------------------------------------
      ngen=1
      thrsym=min(0.1d+00,eps)
      ntold=1
      nt=1
c -- all point groups have identity
      call SetDiagMat(3,one,trans(1,1))
c ---------------------------------------------------------------------
c     identify linear/planar molecules
c ---------------------------------------------------------------------
      call lipl(xyz,n,aline,aplane,dline,dplane)
      if(dline.lt.eps) then
        e3(1)=aline(1)
        e3(2)=aline(2)
        e3(3)=aline(3)
        call e1e2e3(e1,e2,e3,ierr)
        if(ierr.ne.0) return
        call traxyz(xyz,e1,e2,e3,n)
        luck=1
        do 110 i=1,n
          if ( luck.eq.1.and.xyz(3,i).gt.eps ) then
             lentil=0
             do 120 j=1,n
                if ( lentil.eq.0.and.xyz(3,j).lt.-eps) then
                   if ( abs(xyz(3,j)+xyz(3,i)).lt.eps.and.
     1                  mtype(i).eq.mtype(j) ) lentil=1
                endif
  120           continue
             luck=lentil
          endif
  110     continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        if ( luck.eq.1 ) then      !  this section changed (jb)
           sflies='d*h '
           if(n.eq.1) sflies='ih  '
c -- set up as if d2h
           nt = 8
           call SetDiagMat(3,one,trans(1,2))
           call SetDiagMat(3,one,trans(1,3))
           call SetDiagMat(3,one,trans(1,4))
           trans(1,2) = -one
           trans(5,3) = -one
           trans(9,4) = -one
           call SetDiagMat(3,-one,trans(1,5))
           call SetDiagMat(3,-one,trans(1,6))
           call SetDiagMat(3,-one,trans(1,7))
           call SetDiagMat(3,-one,trans(1,8))
           trans(1,5) = one
           trans(5,6) = one
           trans(9,7) = one
        else
           sflies='c*v '
c -- set up as if c2v
           nt = 4
           call SetDiagMat(3,one,trans(1,2))
           call SetDiagMat(3,one,trans(1,3))
           call SetDiagMat(3,-one,trans(1,4))
           trans(1,2) = -one
           trans(5,3) = -one
           trans(9,4) =  one
        endif
c -----------------------------------------------------------------
c       make sure linear molecule lies ONLY along Z-axis
c -----------------------------------------------------------------
        do 150 i=1,n
        xyz(1,i) = zero
        xyz(2,i) = zero
  150   continue
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ierr=0
        return
c
      elseif(dplane.lt.eps) then
c ---------------------------------------------------------------------
c       planar molecule must lie in xy plane
c ---------------------------------------------------------------------
        e3(1)=aplane(1)
        e3(2)=aplane(2)
        e3(3)=aplane(3)
        call gencs(e3,gen)
        call grpgrp(gen,ngen,trans,ntold,nt,ntrans,thrsym,
     1              ipath,jpath,ierr)
        if(ierr.ne.0) return
        if(nt.gt.ntold) then
          ngen=ngen+1
          ntold=nt
        endif
      endif
c ---------------------------------------------------------------------
c     distance from center of mass -> rtoc(1...n)
c ---------------------------------------------------------------------
      call getdis(xyz,rtoc,n,ierr)
      if(ierr.ne.0) return
c ---------------------------------------------------------------------
c     build shells of atoms (within each shell, distance between
c     atom and center of mass is equal with an absolute accuracy
c     of <eps>)
c ---------------------------------------------------------------------
      epsshl=0.5d0*eps
      call rshell(xyz,xyzsh,mtype,mtypsh,rtoc,rtocsh,labsh,
     1            npersh,nlab,n,epsshl)
      icatom=ionc(rtoc,n)
c ---------------------------------------------------------------------
c      search center of symmetry
c ---------------------------------------------------------------------
      luck=1
      do 200 i=1,n
         x = xyz(1,i)
         y = xyz(2,i)
         z = xyz(3,i)
         do 210 j=1,n
            if ( mtype(j).eq.mtype(i) .and.
     1           abs(x+xyz(1,j)).lt.eps .and.
     2           abs(y+xyz(2,j)).lt.eps .and.
     3           abs(z+xyz(3,j)).lt.eps ) goto 200
            if ( j.eq.n ) then
               luck = 0
               goto 220
            endif
  210       continue
  200    continue
  220 continue
      if ( luck.eq.1 ) then
        gen(1) = -one
        gen(2) =  zero
        gen(3) =  zero
        gen(4) =  zero
        gen(5) = -one
        gen(6) =  zero
        gen(7) =  zero
        gen(8) =  zero
        gen(9) = -one
        call grpgrp(gen,ngen,trans,ntold,nt,ntrans,thrsym,
     1              ipath,jpath,ierr)
        if ( ierr.ne.0 ) return
        if ( nt.gt.ntold ) ntold=nt
      endif
c ---------------------------------------------------------------------
c      search c(n)/s(n) axes to build up group
c ---------------------------------------------------------------------
      call cnsrch(xyzsh,rtocsh,trans,univec,istype,ipath,jpath,labsh,
     1            logish,npersh,n,nlab,ntrans,ncnout,icatom,eps,
     2            ntold,nt,ngen,thrsym,ierr)
      if(ierr.ne.0) return
c ---------------------------------------------------------------------
c      search mirror plane(s) to build up group
c ---------------------------------------------------------------------
      call cssrch(xyzsh,trans,univec,istype,ipath,jpath,
     1            labsh,logish,npersh,n,nlab,ntrans,ncsout,
     2            eps,ntold,nt,ngen,thrsym,ierr)
      if(ierr.ne.0) return
c ----- ngen = number of generators encountered (possibly not minimal)
      ngen=ngen-1
c ---------------------------------------------------------------------
c     assign index array for inverse operations
c ---------------------------------------------------------------------
      call syminv(trans,nt,invt,thrsym,ierr)
      if(ierr.ne.0) return
c ---------------------------------------------------------------------
c     classify symmetry operations
c ---------------------------------------------------------------------
      do 9990 it=1,nt
 9990   istype(it)=0
      ncnmax=0
      icnmax=0
      nsnmax=0
      isnmax=0
      call analys(trans,univec,thrsym,invt,istype,nt,
     1            icnmax,ncnmax,isnmax,nsnmax,0,ierr)
      if(ierr.ne.0) return
c
      iz=0
      nz=0
      iy=0
      ix=0
      ic31=0
      ic32=0
c ---------------------------------------------------------------------
c     identify group and assign 'canonical' axes
c ---------------------------------------------------------------------
      call ggroup(univec,istype,nt,ncnmax,icnmax,nsnmax,isnmax,
     1            iz,nz,iy,ix,ic31,ic32,thrsym,sflies,ierr)
      if(ierr.ne.0) return
c
      do 900 i=1,n
        xyzsh(1,i)=xyz(1,i)
        xyzsh(2,i)=xyz(2,i)
  900   xyzsh(3,i)=xyz(3,i)
c ---------------------------------------------------------------------
c     adjust geometry to meet the requirement of canonical axes and/or
c     mirror planes; modified coordinates on <xyz>
c ---------------------------------------------------------------------
      if(adjst.and.nt.gt.1)
     $ call adjust1(trans,univec,xyz,nt,n,iz,iy,ix,ic31,ic32,ierr)
c
      return
      end
c ======================================================================
c
      subroutine adjust1(t,univec,xyz,nt,n,iz,iy,ix,ic31,ic32,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      adjust geometry
c      iz != 0 => symmetry element iz must match z axis
c      iy != 0 => symmetry element iy must match y axis
c                 (e.g. xz mirror plane)
c      ix != 0 => symmetry element ix must match x axis
c      ic31 != 0 => s.e. ic31 must match c3(1,1,1) for t,o groups
c      ic32 != 0 => s.e. ic32 must match c3(strange) for i groups
c =====================================================================
      dimension t(3,3,nt),univec(3,nt),xyz(3,n)
      dimension e1(3),e2(3),e3(3)
cc      parameter (adjtol=0.001d0)
      parameter (adjtol=0.0001d0)      ! JB  Jan 2008
c
      ierr = -1
c
      if(ic31.gt.0.and.ic32.gt.0) return
      ic3=ic31+ic32
      itra=1
c
      if(iz.gt.0.and.iy.eq.0.and.ix.eq.0.and.ic3.eq.0) then
c ---------------------------------------------------------------------
c       Cn / Sn must match z-axis
c ---------------------------------------------------------------------
        e3(1)=univec(1,iz)
        e3(2)=univec(2,iz)
        e3(3)=univec(3,iz)
        call e1e2e3(e1,e2,e3,ierr)
        if(ierr.ne.0) return
c
      elseif(iz.gt.0.and.iy.gt.0.and.ix.eq.0.and.ic3.eq.0) then
c ---------------------------------------------------------------------
c       Cn/Sn must match z-axis, sigma_v plane must match xz-plane
c ---------------------------------------------------------------------
        e2(1)=univec(1,iy)
        e2(2)=univec(2,iy)
        e2(3)=univec(3,iy)
        e3(1)=univec(1,iz)
        e3(2)=univec(2,iz)
        e3(3)=univec(3,iz)
        call cross(e2,e3,e1)
c
      elseif(iz.gt.0.and.ix.gt.0) then
c ---------------------------------------------------------------------
c       Cn/Sn must match z-axis, C2 axis must match x-axis
c ---------------------------------------------------------------------
        e1(1)=univec(1,ix)
        e1(2)=univec(2,ix)
        e1(3)=univec(3,ix)
        e3(1)=univec(1,iz)
        e3(2)=univec(2,iz)
        e3(3)=univec(3,iz)
        call cross(e3,e1,e2)
c
      elseif(iz.gt.0.and.ic31.gt.0) then
c ---------------------------------------------------------------------
c       C3(1,1,1) : build i and j unit vectors from projecting
c       k out of C3 and subsequent rotation (by plus/minus 45
c       degres) of the resultant vector
c ---------------------------------------------------------------------
        e3(1)=univec(1,iz)
        e3(2)=univec(2,iz)
        e3(3)=univec(3,iz)
        scalp=univec(1,iz)*univec(1,ic31)+
     1        univec(2,iz)*univec(2,ic31)+
     2        univec(3,iz)*univec(3,ic31)
        e1(1)=univec(1,ic31)-scalp*univec(1,iz)
        e1(2)=univec(2,ic31)-scalp*univec(2,iz)
        e1(3)=univec(3,ic31)-scalp*univec(3,iz)
        call anorm(e1,re1,ierr)
        if(ierr.eq.-1) return
        call cross(e3,e1,e2)
        rt2inv=1.d0/sqrt(2.d0)
        emi1=e1(1)-e2(1)
        epl1=e1(1)+e2(1)
        e1(1)=rt2inv*emi1
        e2(1)=rt2inv*epl1
        emi2=e1(2)-e2(2)
        epl2=e1(2)+e2(2)
        e1(2)=rt2inv*emi2
        e2(2)=rt2inv*epl2
        emi3=e1(3)-e2(3)
        epl3=e1(3)+e2(3)
        e1(3)=rt2inv*emi3
        e2(3)=rt2inv*epl3
c
      elseif(iz.gt.0.and.ic32.gt.0) then
c ---------------------------------------------------------------------
c     C3 for I,Ih is part of the xz-plane => get i unit vector from
c     projecting k out of C3 ; j = vector product k times i
c ---------------------------------------------------------------------
        e3(1)=univec(1,iz)
        e3(2)=univec(2,iz)
        e3(3)=univec(3,iz)
        scalp=univec(1,iz)*univec(1,ic32)+
     1        univec(2,iz)*univec(2,ic32)+
     2        univec(3,iz)*univec(3,ic32)
        e1(1)=univec(1,ic32)-scalp*univec(1,iz)
        e1(2)=univec(2,ic32)-scalp*univec(2,iz)
        e1(3)=univec(3,ic32)-scalp*univec(3,iz)
        call anorm(e1,re1,ierr)
        if(ierr.eq.-1) return
        call cross(e3,e1,e2)
      else
        itra=0
      endif
c
c ---------------------------------------------------------------------
c     check whether the transformation e1,e2,e3 changes <t>
c ---------------------------------------------------------------------
      ktra=0
      do 100 i = 1,nt
         urxx = e1(1)*t(1,1,i) + e1(2)*t(2,1,i) + e1(3)*t(3,1,i)
         urxy = e2(1)*t(1,1,i) + e2(2)*t(2,1,i) + e2(3)*t(3,1,i)
         urxz = e3(1)*t(1,1,i) + e3(2)*t(2,1,i) + e3(3)*t(3,1,i)
         uryx = e1(1)*t(1,2,i) + e1(2)*t(2,2,i) + e1(3)*t(3,2,i)
         uryy = e2(1)*t(1,2,i) + e2(2)*t(2,2,i) + e2(3)*t(3,2,i)
         uryz = e3(1)*t(1,2,i) + e3(2)*t(2,2,i) + e3(3)*t(3,2,i)
         urzx = e1(1)*t(1,3,i) + e1(2)*t(2,3,i) + e1(3)*t(3,3,i)
         urzy = e2(1)*t(1,3,i) + e2(2)*t(2,3,i) + e2(3)*t(3,3,i)
         urzz = e3(1)*t(1,3,i) + e3(2)*t(2,3,i) + e3(3)*t(3,3,i)
         urutxx = urxx*e1(1) + uryx*e1(2) + urzx*e1(3)
         urutxy = urxy*e1(1) + uryy*e1(2) + urzy*e1(3)
         urutxz = urxz*e1(1) + uryz*e1(2) + urzz*e1(3)
         urutyx = urxx*e2(1) + uryx*e2(2) + urzx*e2(3)
         urutyy = urxy*e2(1) + uryy*e2(2) + urzy*e2(3)
         urutyz = urxz*e2(1) + uryz*e2(2) + urzz*e2(3)
         urutzx = urxx*e3(1) + uryx*e3(2) + urzx*e3(3)
         urutzy = urxy*e3(1) + uryy*e3(2) + urzy*e3(3)
         urutzz = urxz*e3(1) + uryz*e3(2) + urzz*e3(3)
         do 200 j = 1,nt
            if ( abs(urutxx - t(1,1,j)).lt.adjtol.and.
     1           abs(urutxy - t(1,2,j)).lt.adjtol.and.
     2           abs(urutxz - t(1,3,j)).lt.adjtol.and.
     3           abs(urutyx - t(2,1,j)).lt.adjtol.and.
     4           abs(urutyy - t(2,2,j)).lt.adjtol.and.
     5           abs(urutyz - t(2,3,j)).lt.adjtol.and.
     6           abs(urutzx - t(3,1,j)).lt.adjtol.and.
     7           abs(urutzy - t(3,2,j)).lt.adjtol.and.
     8           abs(urutzz - t(3,3,j)).lt.adjtol ) go to 100
  200       continue
         ktra=1
         GO TO 190
  100    continue
  190 continue
c
cc      if ( itra.eq.1.and.ktra.eq.0 ) write(6,1000)
cc 1000    format(/,' adjustment skipped since molecule oriented ',/)
c
      if ( itra.eq.1.and.ktra.eq.1 ) call traxyz(xyz,e1,e2,e3,n)
c
      ierr = 0
c
      return
      end
c ======================================================================
c
      subroutine analys(trans,univec,thrsym,invt,istype,nt,
     1                  icnmax,ncnmax,isnmax,nsnmax,iprint,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     starting from the representation matrices <trans(9,nt)>,
c     classify symmetry operations -> <istype>
c         1  -> unit operation
c        -1  -> reflection
c         2,... -> rotation
c        -2,... -> improper rotation
c     nt = number of symmetry operations
c     icnmax,ncnmax = no. and <n> of C<n> axis with maximum <n>
c     isnmax,nsnmax = no. and <n> of S<n> axis with maximum <n>
c     univec(3,nt) = unit vectors representing either rotation axis
c                    or axis orthogonal to mirror plane
c =====================================================================
      dimension trans(9,nt),univec(3,nt),invt(nt),istype(nt),
     1          gen(6),eigvek(9)
c
      eppes = thrsym
c
      do 9991 it=1,nt
        itinv=invt(it)
        if(itinv.gt.it) goto 9991
c ---------------------------------------------------------------------
c        build R + R**(-1) (symmetric!)
c ---------------------------------------------------------------------
        gen(1)=trans(1,it)
        gen(2)=0.5d0*(trans(2,it)+trans(4,it))
        gen(3)=trans(5,it)
        gen(4)=0.5d0*(trans(3,it)+trans(7,it))
        gen(5)=0.5d0*(trans(6,it)+trans(8,it))
        gen(6)=trans(9,it)
        call SetDiagMat(3,1.0d0,eigvek)
        erdeps=1.d-14
        call erduw(gen,eigvek,3,erdeps,ierr)
        if(ierr.ne.0) return
        eig1=gen(1)
        eig2=gen(3)
        eig3=gen(6)
        party=deter(trans(1,it))
        trace=eig1+eig2+eig3
c ---------------------------------------------------------------------
        if(party.gt.0.d+00) then
          iparty=1
        else
          iparty=-1
        endif
c ---------------------------------------------------------------------
        if(abs(trace-3.d0).lt.eppes) then
          istype(it)=1
        elseif(abs(trace+3.d0).lt.eppes) then
          istype(it)=-2
        elseif(abs(trace-1.d0).lt.eppes) then
          if(iparty.eq.1) then
            istype(it)=4
          else
            istype(it)=-1
          endif
        elseif(abs(trace+1.d0).lt.eppes) then
          if(iparty.eq.1) then
            istype(it)=2
          else
            istype(it)=-4
          endif
        else
c ---------------------------------------------------------------------
c         the classification of S operations is extremely difficult.
c         if n is odd, S^(2n) = 1.
c         if n is even, S^(n) = 1.
c         therefore, an Sn and an S(2n) operation for n odd cannot
c         be distinguished based on for which x S^x yields 1.
c
c         assume it is found that Sy^x = 1 for x = n+n with n odd.
c         then, there are two possibilities :
c
c         (a) y = 2n, n odd.
c
c          => S(PI/n) * i = C(PI/n + PI). since n is odd,
c             n repeats of C will yield PI + n*PI = 0 mod 2*PI.
c
c         (b) y = n, n odd.
c
c          => S(2PI/n) * i = C(2PI/n + PI). n repeats of C
c             will yield 2*PI + n*PI = PI mod 2*PI since n is odd.
c
c         in conclusion, if y is even, then n repeats of C
c         suffice. if y is odd, then one needs n+n repeats.
c
c         This means that an Sn operation with n odd is classified
c         as such IFF it takes n+n cycles for Sn to yield 1 and
c         it takes n+n cycles for C=Sn*i to yield 1.
c
c         if it takes only n cycles for C=Sn*i to yield 1 then it
c         is an S2n operation (with n odd).
c
c         An S2n operation with n even needs 2n cycles to yield 1.
c         S2n*i needs 2n cycles as well, since (PI/n+PI)*n for n
c         even yields PI + n*PI = PI mod 2*PI.
c ---------------------------------------------------------------------
          if ( iparty.ne.1 ) then
            nrot=ncycle(trans(1,it),eppes)
            do 210 i=1,9
  210          trans(i,it)=-trans(i,it)
            nrotc=ncycle(trans(1,it),eppes)
            do 211 i=1,9
  211          trans(i,it)=-trans(i,it)
            if ( nrotc.eq.nrot.and.mod(nrot/2,2).eq.1 ) nrot = nrot/2
          else
            nrot=ncycle(trans(1,it),eppes)
          endif
c
          if(nrot.eq.100) then
c           write(6,701) it
c 701       format(/,' operation no. ',i4,
c    1      ' will generate cyclic group of order 100 or larger ',/,
c    2      ' try to use a sloppier threshold ... ',/)
            nt = 1
            return
          endif
c
          if(iparty.eq.1) then
            istype(it)=nrot
          else
            istype(it)=-nrot
          endif
c
        endif
c ---------------------------------------------------------------------
        if(iparty.eq.1) then
          if(abs(eig1-1.d0).lt.thrsym) then
            k=0
          elseif(abs(eig2-1.d0).lt.thrsym) then
            k=3
          elseif(abs(eig3-1.d0).lt.thrsym) then
            k=6
          endif
          if(istype(it).gt.ncnmax) then
            icnmax=it
            ncnmax=istype(it)
          endif
        else
          if(abs(eig1+1.d0).lt.thrsym) then
            k=0
          elseif(abs(eig2+1.d0).lt.thrsym) then
            k=3
          elseif(abs(eig3+1.d0).lt.thrsym) then
            k=6
          endif
c         ---  -2 operation has no unique axis assigned to it.
          if(istype(it).lt.-1.and.istype(it).lt.nsnmax) then
            isnmax=it
            nsnmax=istype(it)
          endif
        endif
        univec(1,it)=eigvek(k+1)
        univec(2,it)=eigvek(k+2)
        univec(3,it)=eigvek(k+3)
c       --- inverse operation ...
        if(itinv.ne.it) then
          istype(itinv)=istype(it)
          univec(1,itinv)=eigvek(k+1)
          univec(2,itinv)=eigvek(k+2)
          univec(3,itinv)=eigvek(k+3)
        endif
c ------ debug --------------------------------------------------------
c       if(iprint.ne.0)
c    1    write(6,'(2i5,3f20.14)') it,istype(it),(univec(k,it),k=1,3)
c ---------------------------------------------------------------------
 9991   continue
      nsnmax=-nsnmax
c
      if ( nsnmax.eq.2 ) then
         do 9992 it=1,nt
            if ( istype(it).eq.2 ) isnmax = it
 9992       continue
      endif
c
      return
      end
c ======================================================================
c
      subroutine anorm(a,ra,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      normalize vector a
c      eps = normalization threshold ( = 1.d-3 )
c      ierr = -1 if normalization failed and = 0 otherwise
c =====================================================================
      dimension a(3)
      parameter (eps=1.0d-3)
c
      ierr=-1
c
      ra = sqrt(a(1)*a(1)+a(2)*a(2)+a(3)*a(3))
      if (ra.lt.eps) return
c
      a(1)=a(1)/ra
      a(2)=a(2)/ra
      a(3)=a(3)/ra
      ierr=0
c
      return
      end
c ======================================================================
c
      subroutine ccheck(xyzsh,e1,e2,e3,trans,istype,
     1                  eps,rotang,alpha2,labsh,logish,
     2                  ncn,n,ntrans,iatom,ifound,ifold,rotor)
      implicit real*8 (a-h,o-z)
c
      dimension xyzsh(3,n),e1(3),e2(3),e3(3),trans(3,ntrans),
     1          istype(ntrans)
      dimension labsh(n)
      logical logish(n),parall
c
      dimension xyznew(3)
      external rotor
c
c     note that "rotang" has to be initialized by the calling routine !
c
      ifound=0
c
      do 100 irot=1,ncn-1
        rotang=rotang+alpha2
        do 110 k=1,n
          logish(k)=.false.
  110     continue
        logish(iatom)=.true.
        k=0
  120   k=k+1
        if(k.gt.n) goto 100
        if(logish(k)) goto 120
        call rotor(e1,e2,e3,xyzsh(1,k),xyznew,rotang)
        call hitit(xyzsh,labsh,xyznew,labsh(k),eps,n,imatch)
        if(imatch.eq.0) return
        logish(k)=.true.
        logish(imatch)=.true.
        goto 120
  100   continue
c
      l=0
  140 l=l+1
      if(l.eq.ntrans) then
        trans(1,l)=e3(1)
        trans(2,l)=e3(2)
        trans(3,l)=e3(3)
        ifound=1
c ----- debug ---------------------------------------------------------
c       write(6,9101) ifold,e3
c9101   format(/,' found ',i4,'-fold axis : ',3f10.6,/)
c ---------------------------------------------------------------------
        return
      endif
      if(istype(l).ne.ifold.or.
     1   .not.parall(e3,trans(1,l),eps)) goto 140
c
      return
      end
c ======================================================================
c
      subroutine checks(xyzsh,e3,trans,eps,labsh,logish,
     1                  n,ntrans,iatom,ifound)
      implicit real*8 (a-h,o-z)
c
      dimension xyzsh(3,n),e3(3),trans(3,ntrans)
      dimension labsh(n)
      logical logish(n),parall
c
      dimension xyznew(3),center(3)
c
      center(1) = 0.0d0
      center(2) = 0.0d0
      center(3) = 0.0d0
      ifound=0
c
      do 110 k=1,n
        logish(k)=.false.
  110   continue
      logish(iatom)=.true.
      k=0
  120 k=k+1
      if(k.gt.n) goto 130
      if(logish(k)) goto 120
      call rotcs(center,e3,xyzsh(1,k),xyznew)
      call hitit(xyzsh,labsh,xyznew,labsh(k),eps,n,imatch)
      if(imatch.eq.0) goto 130
      logish(k)=.true.
      logish(imatch)=.true.
      goto 120
  130 continue
      if(k.le.n) return
c
      l=0
  140 l=l+1
      if(l.eq.ntrans) then
        trans(1,l)=e3(1)
        trans(2,l)=e3(2)
        trans(3,l)=e3(3)
        ifound=1
c ----- debug ---------------------------------------------------------
c       write(6,9101) e3
c9101   format(/,' found plane of symmetry orthogonal to ',3f10.6,/)
c ---------------------------------------------------------------------
        return
      endif
      if(.not.parall(e3,trans(1,l),eps)) goto 140
c
      return
      end
c ======================================================================
c
      subroutine cnsrch(xyzsh,rtocsh,trans,univec,istype,ipath,jpath,
     1                  labsh,logish,npersh,n,nlab,ntrans,ncnout,
     2                  icatom,eps,ntold,nt,ngen,thrsym,ierr)
      implicit real*8 (a-h,o-z)
c ---------------------------------------------------------------------
c      try to identify C<n> axes
c ---------------------------------------------------------------------
c     --- max. n for a c(n) axis
      parameter ( nofcn = 10 )
      parameter (epsab=1.0d-6, pi=3.14159265358979323844d0)
c
      dimension xyzsh(3,n),rtocsh(n),npersh(nlab),labsh(n)
      dimension trans(9,ntrans),univec(3,ntrans)
      dimension istype(ntrans),ipath(ntrans),jpath(ntrans)
      logical logish(n),parall
      dimension e1(3),e2(3),e3(3)
      dimension bmina(3),ctoab(3),crocus(3)
      dimension gen(9),igen(3)
c
      external rotc,rots
c
      iatom=0
      itrans=0
      igener=0
      monax=0
c ---------------------------------------------------------------------
c     first step : analyze bins containing one atom only
c     if this results in an axis definition, monax differs from zero
c ---------------------------------------------------------------------
      do 100 ilab=1,nlab
        nofat=npersh(ilab)
        if(nofat.gt.1) then
          iatom=iatom+nofat
          goto 100
        endif
c  --- single atom is either at center of mass (skip that)
c      or defines a possible symmetry axis ---
        iatom=iatom+1
        if(iatom.eq.icatom) goto 100
c  --- build axis
        e3(1)=xyzsh(1,iatom)
        e3(2)=xyzsh(2,iatom)
        e3(3)=xyzsh(3,iatom)
        call anorm(e3,e3norm,ierr)
        if(ierr.ne.0) goto 100
        if(itrans+1.gt.ntrans) return
        if(monax.eq.0) goto 120
c  --- compatible symmetry axes => next bin
        if(parall(e3,univec(1,monax),eps)) goto 100
c  --- incompatible symmetry axes => there is none at all !
cc        write(6,9100) ilab
cc 9100   format(/,' symmetry axes are incompatible with bin ',i4,/)
        ncnout=0
        return
  120   continue
c ------- debug ---------------------------------------------------------
cc        write(6,9110) e3
cc 9110   format(/,' found x-fold axis : ',3f10.6,/)
c ---------------------------------------------------------------------
        itrans=itrans+1
        monax=itrans
        univec(1,monax)=e3(1)
        univec(2,monax)=e3(2)
        univec(3,monax)=e3(3)
        istype(monax)=0
        call e1e2e3(e1,e2,e3,ierr)
        if(ierr.ne.0) return
  100   continue
c ---------------------------------------------------------------------
c     continue with those bins containing more than one atom
c ---------------------------------------------------------------------
      ierr = -1
      iatom=0
c
      do 110 ilab=1,nlab
        nofat=npersh(ilab)
        if(nofat.le.0) return
        if(nofat.eq.1) goto 110

c  --- group of atoms
        rmean=0.d0
        do 190 i=iatom+1,iatom+nofat
  190     rmean=rmean+rtocsh(i)
        rmean=rmean/nofat
c ---------------------------------------------------------------------
c    Z = (A+B)/2  M = center of mass  MC = rotation axis looked for
c ---------------------------------------------------------------------
        do 200 i=iatom+1,iatom+nofat-1
          xa=xyzsh(1,i)
          ya=xyzsh(2,i)
          za=xyzsh(3,i)
          ra=rtocsh(i)
          katom=i
          do 300 j=i+1,iatom+nofat
            xb=xyzsh(1,j)
            yb=xyzsh(2,j)
            zb=xyzsh(3,j)
            rb=rtocsh(j)
c  --- |MA|+|MB| = rmin !=! 2*R (ra = rb with accuracy eps !)
            rmin=ra+rb
c  --- ZM -> "ctoab"
            ctoab(1)=-0.5d0*(xa+xb)
            ctoab(2)=-0.5d0*(ya+yb)
            ctoab(3)=-0.5d0*(za+zb)
c  --- |ZM| = abc
            call anorm(ctoab,abc,ispecl)
c ------------- debug -------------------------------------------------
c           if(ispecl.ne.0) then
c             write(6,9300)
c    1          i,j,(xyzsh(k,i),k=1,3),(xyzsh(k,j),k=1,3),
c    2          (ctoab(k),k=1,3)
c9300         format(/,' center of atom pair ',i4,'-',i4,
c    1               /,' with coordinates ',/,3f10.6,/,3f10.6,
c    2               /,' is located at ',3f10.6,' but ',
c    1               /,' coincides with center of mass ',/)
c             goto 300
c           endif
c ---------------------------------------------------------------------
c  --- BA -> "bmina"
            bmina(1)=xb-xa
            bmina(2)=yb-ya
            bmina(3)=zb-za
c  --- |BA| = abdist
            call anorm(bmina,abdist,ierr)
            if(ierr.ne.0) return
c  --- e(MC') = e(ZM) x e(BA) -> "crocus"
            call cross(ctoab,bmina,crocus)
            call anorm(crocus,cmk,icmk)
c  --- center of A and B = M
            if ( icmk.ne.0 ) goto 300
c ------------- debug -------------------------------------------------
c           if(icmk.ne.0) then
c             write(6,9310) i,j
c9310         format(/,' center of mass and atoms ',i4,',',i4,
c    1               /,' share the same line ',/)
c             goto 300
c           endif
c ---------------------------------------------------------------------
c  --- n=2 => e(MC') is candidate
            e3(1)=ctoab(1)
            e3(2)=ctoab(2)
            e3(3)=ctoab(3)
c  --- get the three (orthogonal) rotation axes
            call e1e2e3(e1,e2,e3,ierr)
            if(ierr.ne.0) return
c ---------------------------------------------------------------------
c  if monax.ne.0, check this axis on consistency
c ---------------------------------------------------------------------
            if(monax.ne.0.and.
     1         .not.parall(e3,univec(1,monax),eps)) then
              goto 130
            endif
            if(itrans+1.gt.ntrans) return
            alpha2=pi
            alphag=alpha2
            ncn=2
            rotang=0.d0
            call ccheck(xyzsh,e1,e2,e3,univec,istype,eps,
     1           rotang,alpha2,labsh,logish,ncn,n,itrans+1,katom,
     2           ifound,2,rotc)
c
            if(ifound.eq.1) then
              istype(itrans+1)=2
              ifacs=1
c --------------------------------------------------------------------
c    check for s(4) if c(2) has been found and
c    overwrite istype if s(4) has been found, too
c ---------------------------------------------------------------------
              alpha2=pi
c  --- ncn=3 since ncn-1 = 2 = 4/2 ! (see ccheck)
              ncn=3
              rotang=-0.5d0*pi
              agensn=-rotang
              call ccheck(xyzsh,e1,e2,e3,univec,istype,eps,
     1             rotang,alpha2,labsh,logish,ncn,n,itrans+1,katom,
     2             ifound,-4,rots)
c
              if(ifound.eq.1) then
                istype(itrans+1)=-4
                ifacs=-1
                alphag=agensn
              endif
c
              itrans=itrans+1
c ---------------------------------------------------------------------
c    try to use this symmetry element as a generator
c ---------------------------------------------------------------------
              call gencn(e1,e2,e3,alphag,ifacs,gen)
              call grpgrp(gen,ngen,trans,ntold,nt,ntrans,thrsym,
     1                    ipath,jpath,ierr)
              if(ierr.ne.0) then
                ierr=0
                nt=ntold
                itrans=itrans-1
                goto 130
              endif
              if(nt.gt.ntold) then
                ngen=ngen+1
                igener=igener+1
                igen(igener)=istype(itrans)
                ntold=nt
                if(igener.eq.2) then
                  i1=abs(igen(1))
                  i2=abs(igen(2))
                  if(i1.gt.i2) then
                    igen(1)=i1
                    igen(2)=i2
                  else
                    igen(1)=i2
                    igen(2)=i1
                  endif
                  if((igen(1).eq.4.and.igen(2).eq.3).or.
     1               (igen(1).eq.5.and.igen(2).eq.3).or.
     2               (igen(1).eq.6.and.igen(2).eq.4).or.
     3               (igen(1).eq.10.and.igen(2).eq.6)) return
                endif
                if(igener.eq.3) return
              endif
              if ( nt.eq.120 ) return
            endif
  130       continue
c ---------------------------------------------------------------------
c  --- skip if there are no more than two atoms
            if(nofat.le.2) goto 300
c ---------------------------------------------------------------------
c  --- ncnmax = pi / arctan(|BA|/(|MA|+|MB|))
            ncnmax=nint(pi/atan(abdist/rmin))
            if(ncnmax.gt.nofcn) ncnmax=nofcn
c ----- debug ---------------------------------------------------------
c           write(6,661) i,j,ilab,ncnmax
c 661       format(' atoms ',2i4,' out of bin ',i4,' -> ncn = ',i4)
c ---------------------------------------------------------------------
c  --- n = 3,...,ncnmax (there MUST be a better guess)
            do 400 icn=3,ncnmax
c  --- pi/n => tan(pi/n)
              alpha=pi/icn
              talpha=tan(alpha)
              xiofn=abdist/(talpha+talpha)
c  --- xi(n) = [ |BA|/(2*tan(pi/n)) ] ** 2
              xiofn=xiofn*xiofn
c  --- a,b
              a=xiofn/abc
              b=xiofn-a*a
              if(b.lt.-epsab) goto 400
              if(b.lt.0.d0) b=0.d0
              b=sqrt(b)
              do 450 ipm=-1,1,2
c  --- MC (ipm is -1 or 1)
              e3(1)=(a-abc)*ctoab(1)+ipm*b*crocus(1)
              e3(2)=(a-abc)*ctoab(2)+ipm*b*crocus(2)
              e3(3)=(a-abc)*ctoab(3)+ipm*b*crocus(3)
c  --- axis of rotation (if M=C, we use crocus)
              call anorm(e3,e3norm,ierr)
              if(ierr.ne.0) then
                e3(1)=crocus(1)
                e3(2)=crocus(2)
                e3(3)=crocus(3)
              endif
c  --- get the rotation dreibein
              call e1e2e3(e1,e2,e3,ierr)
              if(ierr.ne.0) return
              if(monax.ne.0.and.
     1           .not.parall(e3,univec(1,monax),eps)) then
                goto 450
              endif
c
              if(itrans+1.gt.ntrans) return
c
c  --- rotate (by 2*pi/n)
              alpha2=alpha+alpha
              alphag=alpha2
              ncn=icn
              rotang=0.d0
              call ccheck(xyzsh,e1,e2,e3,univec,istype,eps,
     1             rotang,alpha2,labsh,logish,ncn,n,itrans+1,katom,
     2             ifound,icn,rotc)
c
              if(ifound.eq.1) then
                istype(itrans+1)=icn
                ifacs=1
c ---------------------------------------------------------------------
c    check for s(2n) if c(n) has been found and
c    overwrite istype if s(2n) has been found, too
c ---------------------------------------------------------------------
                alpha2=alpha+alpha
                ncn=icn+1
                rotang=-alpha
                agensn=alpha
                call ccheck(xyzsh,e1,e2,e3,univec,istype,eps,
     1               rotang,alpha2,labsh,logish,ncn,n,itrans+1,katom,
     2               ifound,-icn*2,rots)
c
                if(ifound.eq.1) then
                  istype(itrans+1)=-icn*2
                  ifacs=-1
                  alphag=agensn
                endif
c
                itrans=itrans+1
c ---------------------------------------------------------------------
c    try to use this symmetry element as a generator
c ---------------------------------------------------------------------
                call gencn(e1,e2,e3,alphag,ifacs,gen)
                call grpgrp(gen,ngen,trans,ntold,nt,ntrans,thrsym,
     1                      ipath,jpath,ierr)
                if(ierr.ne.0) then
                  ierr=0
                  nt=ntold
                  itrans=itrans-1
                  goto 400
                endif
                if(nt.gt.ntold) then
                  ngen=ngen+1
                  igener=igener+1
                  igen(igener)=istype(itrans)
                  ntold=nt
                  if(igener.eq.2) then
                    i1=abs(igen(1))
                    i2=abs(igen(2))
                    if(i1.gt.i2) then
                      igen(1)=i1
                      igen(2)=i2
                    else
                      igen(1)=i2
                      igen(2)=i1
                    endif
                    if((igen(1).eq.4.and.igen(2).eq.3).or.
     1                 (igen(1).eq.5.and.igen(2).eq.3).or.
     2                 (igen(1).eq.6.and.igen(2).eq.4).or.
     3                 (igen(1).eq.10.and.igen(2).eq.6)) return
                  endif
                  if(igener.eq.3) return
                endif
                if ( nt.eq.120 ) return
              endif
  450         continue
  400         continue
  300       continue
  200     continue
  110   iatom=iatom+nofat
      ncnout=itrans
cc      if(monax.ne.0.and.ncnout.eq.1) write(6,901)
cc  901 format(/,' concludingly, the molecule is either linear or ',
cc     1       /,' possesses no symmetry axis at all ',/)
c
      return
      end
c ======================================================================
c
      subroutine cssrch(xyzsh,trans,univec,istype,ipath,jpath,
     1                  labsh,logish,npersh,n,nlab,ntrans,ncsout,eps,
     2                  ntold,nt,ngen,thrsym,ierr)
      implicit real*8 (a-h,o-z)
c ---------------------------------------------------------------------
c      try to identify cs planes
c ---------------------------------------------------------------------
      dimension xyzsh(3,n),trans(9,ntrans),univec(3,ntrans),
     1          istype(ntrans),ipath(ntrans),jpath(ntrans),npersh(nlab)
      dimension labsh(n)
      logical logish(n)
      dimension bmina(3),gen(9)
c
      ierr = 0
      iatom=0
      itrans=0
c
      do 100 ilab=1,nlab
        nofat=npersh(ilab)
        if(nofat.eq.1) then
c  --- single atom is useless for definition of symmetry planes
c      (but it must be shared by all of them!)
          iatom=iatom+1
        elseif(nofat.gt.1) then
c  --- group of atoms
c ---------------------------------------------------------------------
c    Z = (A+B)/2 is candidate for point on plane, e(B-A) => vector
c    orthogonal to plane
c ---------------------------------------------------------------------
          do 200 i=iatom+1,iatom+nofat-1
            xa=xyzsh(1,i)
            ya=xyzsh(2,i)
            za=xyzsh(3,i)
            katom=i
            do 300 j=i+1,iatom+nofat
              xb=xyzsh(1,j)
              yb=xyzsh(2,j)
              zb=xyzsh(3,j)
c  --- BA -> "bmina"
              bmina(1)=xb-xa
              bmina(2)=yb-ya
              bmina(3)=zb-za
c  --- |BA| = abdist
              call anorm(bmina,abdist,iab)
              if(iab.ne.0) goto 300
  400         continue
              if(itrans+1.gt.ntrans) return
              ndim=itrans+1
              call checks(xyzsh,bmina,univec,
     1             eps,labsh,logish,n,ndim,katom,ifound)
              if(ifound.eq.1) then
                itrans=itrans+1
                istype(itrans)=-1
                call gencs(univec(1,itrans),gen)
                call grpgrp(gen,ngen,trans,ntold,nt,ntrans,thrsym,
     1                      ipath,jpath,ierr)
                if(ierr.ne.0) then
                  ierr=0
                  nt=ntold
                  itrans=itrans-1
                  goto 300
                endif
                if(nt.gt.ntold) then
                  ngen=ngen+1
cc                  write(6,6660) (univec(k,itrans),k=1,3)
cc 6660             format(/,' new generator : cs ',3f10.6)
                  ntold=nt
                  if (nt.eq.ntrans) return
                endif
                if(itrans.eq.3) return
                if(itrans.eq.2) then
                  call cross(univec(1,1),univec(1,2),bmina)
                  call anorm(bmina,abdist,iab)
                  if(iab.ne.0) then
                    itrans=1
                    goto 300
                  endif
                  goto 400
                endif
              endif
  300         continue
  200       continue
          iatom=iatom+nofat
        else
          call message('**WARNING** in SYMMETRY module',
     $                 '  Problems Identifying Cs Planes',0,0)
          ierr = -1
          return
        endif
  100   continue
c
      ncsout=itrans
c
      return
      end
c ======================================================================
c
      double precision function deter(a)
      implicit real*8 (a-h,o-z)
      dimension a(9)
c
      deter=a(1)*(a(5)*a(9)-a(6)*a(8))-
     1      a(2)*(a(4)*a(9)-a(6)*a(7))+
     2      a(3)*(a(4)*a(8)-a(5)*a(7))
c
      return
      end
c ======================================================================
c
      double precision function dif(a,b,n)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
c
      s = 0.0d0
      do 10 k=1,n
  10  s = s + (b(k)-a(k))*(b(k)-a(k))
      dif = sqrt(s)
c
      return
      end
c ======================================================================
c
      subroutine e1e2e3(e1,e2,e3,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     for a given unit vector e3, provide a complete set of
c     orthonormal vectors e1,e2,e3
c
c     return code ierr = -1 indicates failure
c =====================================================================
      dimension e1(3),e2(3),e3(3)
      dimension unit(3,3)
c
      parameter (eps=1.0d-6)
c
      call SetDiagMat(3,1.0d0,unit)
      ierr=-1
      i=0
c
  100 i=i+1
      if(i.gt.3) return
      x=unit(1,i)
      y=unit(2,i)
      z=unit(3,i)
      sp=x*e3(1)+y*e3(2)+z*e3(3)
      x=x-sp*e3(1)
      y=y-sp*e3(2)
      z=z-sp*e3(3)
      r=sqrt(x*x+y*y+z*z)
      if(r.lt.eps) goto 100
      rinv=1.d0/r
      e1(1)=x*rinv
      e1(2)=y*rinv
      e1(3)=z*rinv
      call cross(e3,e1,e2)
c
      ierr=0
c
      return
      end
c ======================================================================
c
      subroutine erduw(a,b,na,epslon,ierr)
      implicit real*8(a-h,o-z)
      dimension a(na*(na+1)/2), b(na*na)
c ------------------------------------------------------------------
c   matrix diagonalization by the jacobi method.
c     a = real symmetric matrix to be diagonalized. it is stored
c       by columns with all sub-diagonal elements omitted, so a(i,j) is
c       stored as a((j*(j-1))/2+i).
c     b = matrix to be multiplied by the matrix of eigenvectors.
c     na = dimension of the matrices.
c     epslon is the convergence criterion for off-diagonal elements.
c     program keeps track of non-diagonal norm and adjusts thresh
c     dynamically. the eventually quadratic convergence of jacobi is
c     exploited to reduce thresh faster at the end.
c ------------------------------------------------------------------
      parameter (zero=0.0d0,half=0.5d0)
c
      ierr=0
      loopc=64
c -- calculate of all ref. nondiag. matrix elements squared
      sumnd=zero
      sum=zero
      ij=1
      do 24 i=1,na
        do 16 j=1,i
          term=a(ij)*a(ij)
          sum=sum+term
            if(i.ne.j) then
            sumnd=sumnd+term
          endif
          ij=ij+1
   16   continue
   24 continue
c ---------------------------------------------------------------------
c |   avoid absolute precision measures - look for THE zero matrix    |
c ---------------------------------------------------------------------
      if (sum.le.zero) then
        call message('**WARNING** in SYMMETRY module',
     $               '  Attempting to Diagonalize Zero Matrix',0,0)
        return
      endif
c
      if(na.eq.1) return
c
      thrshg=sqrt((sum+sumnd)/na)*epslon
      small=sumnd*epslon*na*na
c
c ***** The cycling begins here - will it ever come to an end ? ******
c
   32 continue
      if (sumnd.ge.small) then
        thresh=sqrt(sumnd+sumnd)/na
      else
        thresh=thresh*0.001d0
      endif
      thresh=dmax1(thresh,thrshg)
      n=0
      ij=2
      jj=1
      do 112 j=2,na
      jj=jj+j
      jm1=j-1
      ii=0
      do 104 i=1,jm1
      ii=ii+i
      if(abs(a(ij)).lt.thresh) go to 104
      n=n+1
      sumnd=sumnd-a(ij)*a(ij)
      sum=half*(a(jj)+a(ii))
      term=half*(a(jj)-a(ii))
      amax=sign(sqrt(term*term+a(ij)*a(ij)),term)
      c=sqrt((amax+term)/(amax+amax))
      s=a(ij)/(c*(amax+amax))
      a(ii)=sum-amax
      a(jj)=sum+amax
      a(ij)=zero
      im1=i-1
      if(im1 .gt. 0) then
        ki=ii-i
        kj=jj-j
        do 48 k=1,im1
          ki=ki+1
          kj=kj+1
          term=c*a(ki)-s*a(kj)
          a(kj)=s*a(ki)+c*a(kj)
          a(ki)=term
   48   continue
      endif
   56 if(jm1.ne.i) then
        ip1=i+1
        ik=ii+i
        kj=ij
        do 64 k=ip1,jm1
          kj=kj+1
          term=c*a(ik)-s*a(kj)
          a(kj)=s*a(ik)+c*a(kj)
          a(ik)=term
          ik=ik+k
   64   continue
      endif
   72 if(j.eq.na) go to 88
      jp1=j+1
      ik=jj+i
      jk=jj+j
      do 80 k=jp1,na
      term=c*a(ik)-s*a(jk)
      a(jk)=s*a(ik)+c*a(jk)
      a(ik)=term
      ik=ik+k
   80 jk=jk+k
   88 ki=im1*na
      kj=jm1*na
      do 96 k=1,na
      ki=ki+1
      kj=kj+1
      term=c*b(ki)-s*b(kj)
      b(kj)=s*b(ki)+c*b(kj)
   96 b(ki)=term
  104 ij=ij+1
  112 ij=ij+1
      loopc=loopc-1
cc     write(6,*) 'loopc=',loopc,'  thresh=',thresh
      if(loopc.le.0) then
        ierr=1
        call message('**WARNING** in SYMMETRY module',
     $       '  Problems attaining Diagonalization threshold',0,0)
        return
      endif
  120 if(thresh.gt.thrshg.or.n.gt.0) goto 32
c
      return
      end
c ======================================================================
c
      subroutine gencn(e1,e2,e3,alpha,ifacs,gen)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     3 x 3 representation of proper (ifacs=1) or improper (ifacs=-1)
c     rotation by alpha around axis specified by unit vector e3
c =====================================================================
      dimension e1(3),e2(3),e3(3),gen(3,3)
c
      salpha=sin(alpha)
      calpha=cos(alpha)
c
      gen11 = calpha
      gen21 = salpha
      gen31 = 0.d0
      gen12 = -salpha
      gen22 = calpha
      gen32 = 0.d0
      gen13 = 0.d0
      gen23 = 0.d0
      gen33 = dble(ifacs)
c
      tmp11 = gen11*e1(1)+gen12*e2(1)+gen13*e3(1)
      tmp21 = gen21*e1(1)+gen22*e2(1)+gen23*e3(1)
      tmp31 = gen31*e1(1)+gen32*e2(1)+gen33*e3(1)
      tmp12 = gen11*e1(2)+gen12*e2(2)+gen13*e3(2)
      tmp22 = gen21*e1(2)+gen22*e2(2)+gen23*e3(2)
      tmp32 = gen31*e1(2)+gen32*e2(2)+gen33*e3(2)
      tmp13 = gen11*e1(3)+gen12*e2(3)+gen13*e3(3)
      tmp23 = gen21*e1(3)+gen22*e2(3)+gen23*e3(3)
      tmp33 = gen31*e1(3)+gen32*e2(3)+gen33*e3(3)
c
      gen(1,1) = e1(1)*tmp11+e2(1)*tmp21+e3(1)*tmp31
      gen(2,1) = e1(2)*tmp11+e2(2)*tmp21+e3(2)*tmp31
      gen(3,1) = e1(3)*tmp11+e2(3)*tmp21+e3(3)*tmp31
      gen(1,2) = e1(1)*tmp12+e2(1)*tmp22+e3(1)*tmp32
      gen(2,2) = e1(2)*tmp12+e2(2)*tmp22+e3(2)*tmp32
      gen(3,2) = e1(3)*tmp12+e2(3)*tmp22+e3(3)*tmp32
      gen(1,3) = e1(1)*tmp13+e2(1)*tmp23+e3(1)*tmp33
      gen(2,3) = e1(2)*tmp13+e2(2)*tmp23+e3(2)*tmp33
      gen(3,3) = e1(3)*tmp13+e2(3)*tmp23+e3(3)*tmp33
c
      return
      end
c ======================================================================
c
      subroutine gencs(e3,gen)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     3 x 3 representation of arbitrary reflection at a plane
c     orthogonal to unit vector e3
c =====================================================================
      dimension e3(3),gen(3,3)
      parameter (one=1.0d0,two=2.0d0)
c
      gen(1,1) = one - two*e3(1)*e3(1)
      gen(2,1) = -two*e3(2)*e3(1)
      gen(3,1) = -two*e3(3)*e3(1)
      gen(2,2) = one - two*e3(2)*e3(2)
      gen(3,2) = -two*e3(3)*e3(2)
      gen(3,3) = one - two*e3(3)*e3(3)
      gen(1,2) = gen(2,1)
      gen(1,3) = gen(3,1)
      gen(2,3) = gen(3,2)
c
      return
      end
c ======================================================================
c
      subroutine getdis(xyz,rtoc,n,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      calculate distance from origin for all coordinates <xyz>
c      return distances on <rtoc>
c
c      ierr = -1 if n.le.0
c =====================================================================
      dimension xyz(3,n),rtoc(n)
c
      ierr=-1
      if(n.le.0) return
c
      do 100 i=1,n
        distx=xyz(1,i)
        disty=xyz(2,i)
        distz=xyz(3,i)
        rtoc(i)=sqrt(distx*distx+disty*disty+distz*distz)
  100   continue
c
      ierr=0
c
      return
      end
c ======================================================================
c
      subroutine ggroup(univec,istype,nt,ncnmax,icnmax,nsnmax,isnmax,
     1                  iz,nz,iy,ix,ic31,ic32,thrsym,sflies,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     identify group by its generators and assign Schoenflies symbol
c
c     group elements are represented by an axis (rotation axis or axis
c     orthogonal to mirror plane) and an integer which is 1 for the
c     unit operation, -1 for a reflection, +n for a Cn axis and -n
c     for a Sm axis where (Sm)**n = unit (thus, we have -6 for S3 and
c     S6).
c =====================================================================
      dimension univec(3,nt),istype(nt)
      character*4 sflies
c
      dimension genih(3)
c
      ierr = -1
c ---------------------------------------------------------------------
c     the trivial case : one symmetry element => c1
c ---------------------------------------------------------------------
      if(nt.eq.1) then
        sflies='c1  '
c ---------------------------------------------------------------------
c     two elements : cs (mirror), c2 (rotation), ci (inversion)
c ---------------------------------------------------------------------
      elseif(nt.eq.2) then
        if(istype(nt).eq.-1) then
          sflies='cs  '
          iz=2
        elseif(istype(nt).eq.2) then
          sflies='c2  '
          iz=2
        elseif(istype(nt).eq.-2) then
          sflies='ci  '
        endif
c ---------------------------------------------------------------------
c     cyclic groups Cn, Sn ( n > 2 ) ( having n elements )
c ---------------------------------------------------------------------
      elseif(ncnmax.gt.0.and.nt.eq.ncnmax) then
        sflies='c'
        if(ncnmax.lt.10) then
          write(sflies(2:),2101) ncnmax
 2101     format(i1)
        elseif(ncnmax.lt.100) then
          write(sflies(2:),2102) ncnmax
 2102     format(i2)
        else
          sflies(2:2)='n'
        endif
        iz=icnmax
        nz=ncnmax
      elseif(nsnmax.gt.0.and.nt.eq.nsnmax) then
        sflies='s'
        if(nsnmax.lt.10) then
          write(sflies(2:),2101) nsnmax
        elseif(nsnmax.lt.100) then
          write(sflies(2:),2102) nsnmax
        else
          sflies(2:2)='n'
        endif
        iz=isnmax
        nz=-nsnmax
      else
c ---------------------------------------------------------------------
c     all other groups
c       Cnv, Cnh
c       Dn, Dnd, Dnh
c       T, Th, Td
c       O, Oh
c       I, Ih
c     need some more inspection due to those C2 and C3 axes
c ---------------------------------------------------------------------
        ic2=0
        is2=0
        ich=0
        icv=0
        ic3=0
c
        call symel(univec,istype,nt,icnmax,isnmax,
     1             ic2,is2,ich,icv,ic3,thrsym)

c  --- Dnd
        if(nsnmax.gt.ncnmax.and.is2.gt.0.and.
     1     ic3.eq.0.and.ich.eq.0.and.icv.gt.0) then
          sflies='d'
          if(ncnmax.lt.10) then
            write(sflies(2:),2111) ncnmax,'d'
 2111       format(i1,a1)
          elseif(ncnmax.lt.100) then
            write(sflies(2:),2112) ncnmax,'d'
 2112       format(i2,a1)
          endif
          iz=isnmax
          nz=-nsnmax
          ix=is2
c  --- Dn
        elseif(ncnmax.gt.nsnmax.and.ic2.gt.0.and.
     1         ic3.eq.0.and.ich.eq.0.and.icv.eq.0) then
          sflies='d'
          if(ncnmax.lt.10) then
            write(sflies(2:),2101) ncnmax
          elseif(ncnmax.lt.100) then
            write(sflies(2:),2102) ncnmax
          endif
          iz=icnmax
          nz=ncnmax
          ix=ic2
c  --- Cnh
        elseif(ncnmax.ge.nsnmax.and.ic2.eq.0.and.
     1         ic3.eq.0.and.ich.gt.0.and.icv.eq.0) then
          sflies='c'
          if(ncnmax.lt.10) then
            write(sflies(2:),2111) ncnmax,'h'
          elseif(ncnmax.lt.100) then
            write(sflies(2:),2112) ncnmax,'h'
          endif
          iz=icnmax
          nz=ncnmax
c  --- Cnv
        elseif(ncnmax.gt.nsnmax.and.ic2.eq.0.and.
     1         ic3.eq.0.and.ich.eq.0.and.icv.gt.0) then
          sflies='c'
          if(ncnmax.lt.10) then
            write(sflies(2:),2111) ncnmax,'v'
          elseif(ncnmax.lt.100) then
            write(sflies(2:),2112) ncnmax,'v'
          endif
          iz=icnmax
          nz=ncnmax
          iy=icv
c  --- Ih (probably, nsnmax=5 and nt=120 suffices)
        elseif(ncnmax.eq.5.and.ic3.gt.0.and.nt.eq.120) then
          sflies='ih  '
          iz=icnmax
          nz=ncnmax
          call ihgen(genih,ierr)
          if(ierr.ne.0) return
          ic32=0
          do 912 it=1,nt
            if(ic32.eq.0.and.istype(it).eq.3) then
              scalp=univec(1,iz)*univec(1,it)+
     1              univec(2,iz)*univec(2,it)+
     2              univec(3,iz)*univec(3,it)
              if(abs(scalp-genih(3)).lt.thrsym) ic32=it
            endif
  912       continue
          if(ic32.eq.0) then
           call message('**WARNING** in SYMMETRY module',
     $                  '  C3-axis for group Ih is missing',0,0)
           return
          endif
c  --- I (probably, ncnmax=5 and nt=60 suffices)
        elseif(ncnmax.eq.5.and.ic3.gt.0.and.nt.eq.60) then
          sflies='i   '
          iz=icnmax
          nz=ncnmax
          call ihgen(genih,ierr)
          if(ierr.ne.0) return
          ic32=0
          do 911 it=1,nt
            if(ic32.eq.0.and.istype(it).eq.3) then
              scalp=univec(1,iz)*univec(1,it)+
     1              univec(2,iz)*univec(2,it)+
     2              univec(3,iz)*univec(3,it)
              if(abs(scalp-genih(3)).lt.thrsym) ic32=it
            endif
  911       continue
          if(ic32.eq.0) then
           call message('**WARNING** in SYMMETRY module',
     $                  '  C3-axis for group I is missing',0,0)
           return
         endif
c  --- Td (probably, nsnmax=4 and nt=24 suffices)
        elseif(nsnmax.eq.4.and.ic3.gt.0.and.is2.gt.0) then
          sflies='td  '
          iz=isnmax
          nz=-nsnmax
          ix=is2
          ic31=ic3
c  --- Th (probably, ncnmax=3 and nt=24 suffices)
        elseif(ncnmax.eq.3.and.ic3.gt.0.and.nt.eq.24) then
          sflies='th  '
          iz=itgrep(istype,nt,2)
          if(iz.eq.0) then
           call message('**WARNING** in SYMMETRY module',
     $                  '  Cannot find C2-axis for group Th',0,0)
           return
          endif
          nz=2
          ic31=icnmax
c  --- T
        elseif(ncnmax.eq.3.and.ic3.gt.0.and.nt.eq.12) then
          sflies='t   '
          iz=itgrep(istype,nt,2)
          if(iz.eq.0) then
           call message('**WARNING** in SYMMETRY module',
     $                  '  Cannot find C2-axis for group T',0,0)
           return
          endif
          nz=2
          ic31=icnmax
c  --- Oh (probably, nsnmax=6 and nt.eq.48 suffices)
        elseif(nsnmax.eq.6.and.ncnmax.eq.4.and.ich.gt.0) then
          sflies='oh  '
          iz=icnmax
          nz=ncnmax
          ic31=isnmax
c  --- O (probably, ncnmax=4 and nt=24 suffices)
        elseif(ncnmax.eq.4.and.ic3.gt.0.and.nt.eq.24) then
          sflies='o   '
          iz=icnmax
          nz=ncnmax
          ic31=ic3
c  --- Dnh
        elseif((nsnmax.eq.ncnmax.or.nsnmax.eq.ncnmax*2).and.
     1         ic3.eq.0.and.is2.gt.0.and.ich.gt.0.and.icv.gt.0) then
          sflies='d'
          if(ncnmax.lt.10) then
            write(sflies(2:),2111) ncnmax,'h'
          elseif(ncnmax.lt.100) then
            write(sflies(2:),2112) ncnmax,'h'
          endif
          iz=isnmax
          nz=-nsnmax
          ix=is2
        else
          call message('**WARNING** in SYMMETRY module',
     $     '  Symmetry Operations found correspond to no known group',
     $        0,0)
          return
        endif
      endif
c
      ierr = 0
c
      return
      end
c ======================================================================
c
      subroutine gram(ikont,natoms,x,c)
      implicit real*8 (a-h,o-z)
c ------------------------------------------------------------------
c     build up a right-handed orthonormal system (on array <c>)
c     which is uniquely tied to the orientation of the molecule
c ------------------------------------------------------------------
      dimension x(3,natoms),c(3,3),xhelp(3,3)
      parameter (eps=1.0d-5)
c
      ikont = 0
c
      if(natoms.eq.0) then
        ikont=4
        return
      endif
c
      do 10 i=1,3
        do 20 j=1,3
   20     xhelp(j,i)=0.0d0
   10   xhelp(i,i)=1.0d0
c
      i=0
      idim=0
  100 i=i+1
      if(i.le.natoms) then
        xnext=x(1,i)
        ynext=x(2,i)
        znext=x(3,i)
      elseif(i.le.natoms+3) then
        xnext=xhelp(1,i-natoms)
        ynext=xhelp(2,i-natoms)
        znext=xhelp(3,i-natoms)
      elseif(i.gt.natoms+3) then
        ikont=4
        return
      endif
      do 110 k=1,idim
        a = xnext*c(1,k)+ynext*c(2,k)+znext*c(3,k)
        xnext = xnext - a*c(1,k)
        ynext = ynext - a*c(2,k)
        znext = znext - a*c(3,k)
  110 continue
      rrr=sqrt(xnext*xnext+ynext*ynext+znext*znext)
      if(rrr.gt.eps) then
c  --- a new vector
        idim=idim+1
        rinv = 1.0d0 / rrr
        c(1,idim)=xnext*rinv
        c(2,idim)=ynext*rinv
        c(3,idim)=znext*rinv
c  --- gotcha !
        if(idim.eq.3) return
      endif
      if(i.eq.natoms) ikont=3-idim
      goto 100
c
      end
c ======================================================================
c
      subroutine grpgrp(gen,ngen,u,ntold,nt,ntrans,thrsym,ipath,jpath,
     1                  ierr)
      implicit real*8(a-h,o-z)
c =====================================================================
c     add a generator candidate gen(9) to actual group u(9,1...nt)
c     hopefully enlarging it
c =====================================================================
      dimension u(9,ntrans),gen(9),scr(9)
      dimension ipath(ntrans),jpath(ntrans)
c
      ierr=0
c -- set neutral element
c -- (every point group has the neutral element as first operation)
      nt=ntold
      if(ngen.eq.1) then
        ipath(nt)=0
        jpath(nt)=0
      endif
c
      ncoset = 1
c
      do 220 i=1,9
 220    scr(i) = gen(i)
c ---------------------------------------------------------------------
c     check if generator is already contained in (sub)group
c     of order <ntold>
c ---------------------------------------------------------------------
      inew=1
      do 230 j=1,ntold
        if(dif(u(1,j),scr(1),9).lt.thrsym) inew=0
 230    continue
      if(inew.eq.0) return
c
      ipath(nt+1)=-ngen
      jpath(nt+1)=0
c ---------------------------------------------------------------------
c     form coset for new group element
c ---------------------------------------------------------------------
 245  do 250 it=1,ntold
        nt = nt + 1
        if (nt.gt.ntrans) goto 900
        if (it.gt.1) then
          ipath(nt)=ntold*ncoset+1
          jpath(nt)=it
        endif
        CALL DGemm('N','N', 3, 3, 3, 1.0d0, scr, 3, u(1,it),
     $              3, 0.0d0, u(1,ntold*ncoset+it), 3)
 250  continue
      ncoset = ncoset + 1
c ---------------------------------------------------------------------
c     check if subgroup plus cosets form a group
c ---------------------------------------------------------------------
      do 330 it = 1,nt
        do 320 jcoset = 1,ncoset-1
          CALL DGemm('N','N', 3, 3, 3, 1.0d0, u(1,it), 3,
     $                u(1,1+ntold*jcoset), 3, 0.0d0, scr, 3)
          do 310 kt = 1,nt
            if (dif(u(1,kt),scr(1),9).lt.thrsym) goto 320
 310        continue
c  --- new symmetry operation found !
          ipath(nt+1)=it
          jpath(nt+1)=1+ntold*jcoset
          goto 245
 320      continue
 330    continue
cc      write(6,'(/,i5,a)') nt,' symmetry operations found :'
      return
c
 900  continue
      call message('**WARNING** in SYMMETRY module',
     $   '  Number of operations found greater than maximum',ntrans,0)
      ierr=1
c
      return
      end
c ======================================================================
c
      subroutine hitit(xyzsh,labsh,xyznew,label,eps,nofat,imatch)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      check if coordinate triple <xyznam> is "identical" (with
c      absolute accuracy "eps") to one out of <xyzsh(i),i=1,nofat>
c      no match => imatch = 0
c      otherwise match with atom <imatch>
c =====================================================================
      dimension xyzsh(3,nofat),xyznew(3)
      dimension labsh(nofat)
c
      imatch=0
      k=0
  100 k=k+1
      if(k.gt.nofat) return
      if(labsh(k).ne.label) goto 100
      distx=xyzsh(1,k)-xyznew(1)
      disty=xyzsh(2,k)-xyznew(2)
      distz=xyzsh(3,k)-xyznew(3)
      rdist=sqrt(distx*distx+disty*disty+distz*distz)
      if(rdist.gt.eps) goto 100
      imatch=k
c
      return
      end
c ======================================================================
c
      subroutine ihgen(e3,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     c3 axis for i,ih -> e3
c =====================================================================
      dimension e3(3),gen(6),eigvek(9)
      parameter (pi=3.14159265358979323844d0)
c
      a1=1.d+00
      a2=2.d+00
      a3=3.d+00
      a5=5.d+00
c
      sqr3 = sqrt(a3)
      cosd = sin(a2*pi/a5)/((a1-cos(a2*pi/a5))*sqr3)
      sind = sqrt(a1-cosd**2)
c
      gen1 = a1 - a3*cosd**2/a2
      gen2 = sqr3*cosd/a2
      gen3 = a3*sind*cosd/a2
cc      gen4 = -gen2
      gen5 = -a1/a2
      gen6 = sqr3*sind/a2
cc      gen7 = gen3
cc      gen8 = -gen6
      gen9 = a1-a3*sind**2/a2
c
      gen(1)=gen1
      gen(2)=0.d+00
      gen(3)=gen5
      gen(4)=gen3
      gen(5)=0.d+00
      gen(6)=gen9
c
      call SetDiagMat(3,1.0d0,eigvek)
      erdeps=1.d-14
      call erduw(gen,eigvek,3,erdeps,ierr)
      if(ierr.ne.0) return
c
      thrsym=1.d-9
      if(abs(gen(1)-1.d0).lt.thrsym) then
        k=0
      elseif(abs(gen(3)-1.d0).lt.thrsym) then
        k=3
      elseif(abs(gen(6)-1.d0).lt.thrsym) then
        k=6
      endif
      e3(1)=eigvek(k+1)
      e3(2)=eigvek(k+2)
      e3(3)=eigvek(k+3)
c
      return
      end
c ======================================================================
c
      function ionc(rtoc,n)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      check if any atom coincides with the origin
c      if so, return its index ( absolute accuracy = 1.d-3 )
c =====================================================================
      integer ionc
      dimension rtoc(n)
      parameter (eps=1.0d-3)
c
      ionc=0
      i=0
  100 i=i+1
      if(i.gt.n) return
      if(abs(rtoc(i)).lt.eps) then
        ionc=i
        return
      endif
      goto 100
c
      end
c ======================================================================
c
      function itgrep(istype,nt,itwant)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      itgrep is the first index which fulfills istype(itgrep)=itwant
c =====================================================================
      integer itgrep
      dimension istype(nt)
c
      it=0
  100 it=it+1
      if(it.gt.nt) then
        itgrep=0
        return
      endif
      if(istype(it).ne.itwant) goto 100
      itgrep=it
c
      return
      end
c ======================================================================
c
      subroutine lipl(xyz,n,aline,aplane,dline,dplane)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      try to fit set of coordinates on <xyz> by a straight line or
c      by a flat plane
c      the resulting axes are on aline/aplane;
c      the root mean deviations are on dline/dplane
c =====================================================================
      dimension xyz(3,n),aline(3),aplane(3)
      dimension audick(6),auweia(3,3)
      dimension devlin(3),devpln(3)
      parameter (zero=0.0d0,Two=2.0d0)
c
      call SetDiagMat(3,1.0d0,auweia)
      xx=zero
      xy=zero
      xz=zero
      yy=zero
      yz=zero
      zz=zero
c
      do 400 i=1,n
        x=xyz(1,i)
        y=xyz(2,i)
        z=xyz(3,i)
        xx=xx+x*x
        xy=xy+x*y
        xz=xz+x*z
        yy=yy+y*y
        yz=yz+y*z
        zz=zz+z*z
  400   continue
      rr=xx+yy+zz
c
      audick(1)=xx
      audick(2)=xy
      audick(3)=yy
      audick(4)=xz
      audick(5)=yz
      audick(6)=zz
c
      epsln=1.0d-10
      call erduw(audick,auweia,3,epsln,ierr)
c
      a=auweia(1,1)
      b=auweia(2,1)
      c=auweia(3,1)
      devpln(1)=a*a*xx+b*b*yy+c*c*zz+Two*(a*b*xy+a*c*xz+b*c*yz)
      devlin(1)=rr-devpln(1)
      a=auweia(1,2)
      b=auweia(2,2)
      c=auweia(3,2)
      devpln(2)=a*a*xx+b*b*yy+c*c*zz+Two*(a*b*xy+a*c*xz+b*c*yz)
      devlin(2)=rr-devpln(2)
      a=auweia(1,3)
      b=auweia(2,3)
      c=auweia(3,3)
      devpln(3)=a*a*xx+b*b*yy+c*c*zz+Two*(a*b*xy+a*c*xz+b*c*yz)
      devlin(3)=rr-devpln(3)
c
      iline=1
      if(devlin(2).lt.devlin(iline)) iline=2
      if(devlin(3).lt.devlin(iline)) iline=3
c
      iplane=1
      if(devpln(2).lt.devpln(iplane)) iplane=2
      if(devpln(3).lt.devpln(iplane)) iplane=3
c
      aline(1)=auweia(1,iline)
      aline(2)=auweia(2,iline)
      aline(3)=auweia(3,iline)
      dline=devlin(iline)
c
      aplane(1)=auweia(1,iplane)
      aplane(2)=auweia(2,iplane)
      aplane(3)=auweia(3,iplane)
      dplane=devpln(iplane)
c
      return
      end
c ======================================================================
c
      function ncycle(gen,thrsym)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      find ncycle such that gen**(ncyle) = unit matrix
c      (3,3)-matrices assumed
c =====================================================================
      dimension gen(3,3),scr(3,3),unit(3,3)
      logical doof
c
      call SetDiagMat(3,1.0d0,unit)
      call cpyvec(9,gen,scr)
c
      n=0
   20 n=n+1
      if(n.eq.100.or.doof(scr,unit,9,thrsym)) goto 90
      tmp11=gen(1,1)*scr(1,1)+gen(1,2)*scr(2,1)+gen(1,3)*scr(3,1)
      tmp12=gen(1,1)*scr(1,2)+gen(1,2)*scr(2,2)+gen(1,3)*scr(3,2)
      tmp13=gen(1,1)*scr(1,3)+gen(1,2)*scr(2,3)+gen(1,3)*scr(3,3)
      tmp21=gen(2,1)*scr(1,1)+gen(2,2)*scr(2,1)+gen(2,3)*scr(3,1)
      tmp22=gen(2,1)*scr(1,2)+gen(2,2)*scr(2,2)+gen(2,3)*scr(3,2)
      tmp23=gen(2,1)*scr(1,3)+gen(2,2)*scr(2,3)+gen(2,3)*scr(3,3)
      tmp31=gen(3,1)*scr(1,1)+gen(3,2)*scr(2,1)+gen(3,3)*scr(3,1)
      tmp32=gen(3,1)*scr(1,2)+gen(3,2)*scr(2,2)+gen(3,3)*scr(3,2)
      tmp33=gen(3,1)*scr(1,3)+gen(3,2)*scr(2,3)+gen(3,3)*scr(3,3)
      scr(1,1)=tmp11
      scr(2,1)=tmp21
      scr(3,1)=tmp31
      scr(1,2)=tmp12
      scr(2,2)=tmp22
      scr(3,2)=tmp32
      scr(1,3)=tmp13
      scr(2,3)=tmp23
      scr(3,3)=tmp33
      goto 20
   90 ncycle=n
      if (n.eq.100) call message('**WARNING** in SYMMETRY module',
     $    '  Potential Problems with C3 axes',0,0)
c
      return
      end
c ---------------------------------------------------------------------
      function doof(a,b,n,thr)
      implicit real*8 (a-h,o-z)
      dimension a(n),b(n)
      logical doof
c
      do 100 i = 1,n
         if ( abs(a(i)-b(i)).gt.thr ) then
            doof = .false.
            return
         endif
  100    continue
c
      doof = .true.
c
      return
      end
c ======================================================================
c
      function orthog(a,b,thr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     orthog is .true. if unit vector a is orthogonal to unit vector b
c =====================================================================
      logical orthog
      dimension a(3),b(3)
c
      tmp=abs(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))
      orthog=tmp.lt.thr
c
      return
      end
c ======================================================================
c
      function parall(a,b,thr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     parall is .true. if unit vector a is parallel to unit vector b
c =====================================================================
      logical parall
      dimension a(3),b(3)
c
      tmp=abs(a(1)*b(1)+a(2)*b(2)+a(3)*b(3))
      tmp=abs(tmp-1.0d+00)
      parall=tmp.lt.thr
c
      return
      end
c ======================================================================
c
      subroutine rotc(e1,e2,e3,rin,rout,alpha)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      rout = rin rotated by alpha degrees with respect to axis e3
c      (Cn operation)
c =====================================================================
      dimension e1(3),e2(3),e3(3)
      dimension rin(3),rout(3)
c
      salpha=sin(alpha)
      calpha=cos(alpha)
c  --- scalar products of rin with e1,e2,e3
      rte1 = rin(1)*e1(1)+rin(2)*e1(2)+rin(3)*e1(3)
      rte2 = rin(1)*e2(1)+rin(2)*e2(2)+rin(3)*e2(3)
      rte3 = rin(1)*e3(1)+rin(2)*e3(2)+rin(3)*e3(3)
c  --- rotate in (e1,e2) plane
      rnewe1 = rte1 * calpha - rte2 * salpha
      rnewe2 = rte1 * salpha + rte2 * calpha
c  --- build rout
      rout(1) = rnewe1*e1(1) + rnewe2*e2(1) + rte3*e3(1)
      rout(2) = rnewe1*e1(2) + rnewe2*e2(2) + rte3*e3(2)
      rout(3) = rnewe1*e1(3) + rnewe2*e2(3) + rte3*e3(3)
c ---------------------------------------------------------------------
c     matrix representation
c ---------------------------------------------------------------------
      return
      end
c ======================================================================
c
      subroutine rotcs(cplane,e3,rin,rout)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      rout = rin reflected at the plane specified by the point
c      (cplane) and the unit vector (e3) orthogonal to the plane
c =====================================================================
      dimension cplane(3),e3(3)
      dimension rin(3),rout(3)
c
c  --- scalar product of (rin-cplane) with e3
      scalp = (rin(1)-cplane(1))*e3(1)+
     1        (rin(2)-cplane(2))*e3(2)+
     2        (rin(3)-cplane(3))*e3(3)
      scalp = -(scalp+scalp)
c  --- build rout
      rout(1) = rin(1) + scalp*e3(1)
      rout(2) = rin(2) + scalp*e3(2)
      rout(3) = rin(3) + scalp*e3(3)
c
      return
      end
c ======================================================================
c
      subroutine rots(e1,e2,e3,rin,rout,alpha)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      rout = rin rotated by alpha degrees with respect to e3
c      and reflected by the (e1,e2)-plane (s operation)
c =====================================================================
      dimension e1(3),e2(3),e3(3)
      dimension rin(3),rout(3)
c
      salpha=sin(alpha)
      calpha=cos(alpha)
c  --- scalar products of rin with e1,e2,e3
      rte1 = rin(1)*e1(1)+rin(2)*e1(2)+rin(3)*e1(3)
      rte2 = rin(1)*e2(1)+rin(2)*e2(2)+rin(3)*e2(3)
      rte3 = rin(1)*e3(1)+rin(2)*e3(2)+rin(3)*e3(3)
c  --- rotate in (e1,e2) plane
      rnewe1 = rte1 * calpha - rte2 * salpha
      rnewe2 = rte1 * salpha + rte2 * calpha
c  --- build rout
      rout(1) = rnewe1*e1(1) + rnewe2*e2(1) - rte3*e3(1)
      rout(2) = rnewe1*e1(2) + rnewe2*e2(2) - rte3*e3(2)
      rout(3) = rnewe1*e1(3) + rnewe2*e2(3) - rte3*e3(3)
c
      return
      end
c ======================================================================
c
      subroutine rshell(xyz,xyzsh,mtype,mtypsh,
     1                  rtoc,rtocsh,labsh,npersh,nlab,n,eps)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      put coordinate triples into shells ( accuracy = eps )
c      and return shell indices on labsh
c      nlab = number of resulting shells
c      npersh = number of atoms per shell
c      xyzsh = coordinates sorted according to shells
c =====================================================================
      dimension xyz(3,n),xyzsh(3,n)
      character*8 mtype(n),mtypsh(n),ossi
      dimension rtoc(n),rtocsh(n),labsh(n),npersh(n)
c
      do 200 i=1,n
        xyzsh(1,i)=xyz(1,i)
        xyzsh(2,i)=xyz(2,i)
        xyzsh(3,i)=xyz(3,i)
        rtocsh(i)=rtoc(i)
        mtypsh(i)=mtype(i)
        labsh(i)=i
  200   continue
c
      do 210 i=1,n-1
        do 220 j=2,n-i+1
          if ( ( mtype(labsh(j-1)).eq.mtype(labsh(j)).and.
     1           rtoc(j-1).gt.rtoc(j)+eps ) .or.
c -- note: "+eps" added by MM    ! June 2008
     1         lgt(mtype(labsh(j-1)),mtype(labsh(j))) ) then
             tmp = xyzsh(1,j-1)
             xyzsh(1,j-1) = xyzsh(1,j)
             xyzsh(1,j) = tmp
             tmp = xyzsh(2,j-1)
             xyzsh(2,j-1) = xyzsh(2,j)
             xyzsh(2,j) = tmp
             tmp = xyzsh(3,j-1)
             xyzsh(3,j-1) = xyzsh(3,j)
             xyzsh(3,j) = tmp
             tmp = rtoc(j-1)
             rtoc(j-1) = rtoc(j)
             rtoc(j) = tmp
             k = labsh(j-1)
             labsh(j-1) = labsh(j)
             labsh(j) = k
             ossi = mtypsh(j-1)
             mtypsh(j-1) = mtypsh(j)
             mtypsh(j) = ossi
          endif
  220     continue
  210   continue
c
      iatom = 1
      ilab = 1
      do 240 i=2,n
         if ( mtypsh(i).eq.mtypsh(i-1).and.
     1        rtoc(i).lt.rtoc(i-1)+eps ) then
            iatom = iatom + 1
         else
            npersh(ilab) = iatom
            ilab = ilab + 1
            iatom = 1
         endif
  240    continue
      npersh(ilab) = iatom
c
      nlab=ilab
c
      iatom = 0
      do 290 ilab = 1,nlab
        do 291 j = 1,npersh(ilab)
          iatom = iatom + 1
  291     labsh(iatom) = ilab
  290   continue
c ------- debug -------------------------------------------------------
cc      do 400 i=1,n
cc        write(6,602) i,rtoc(i),labsh(i)
cc  602   format(' atom ',i4,f10.6,' ---> bin ',i4)
cc  400   continue
cc      write(6,*)
c ---------------------------------------------------------------------
c
      return
      end
c ======================================================================
c
      subroutine symel(univec,istype,nt,icnmax,isnmax,
     1                 ic2,is2,ich,icv,ic3,thrsym)
      implicit real*8 (a-h,o-z)
c =====================================================================
c      try to assign
c        ic2 = index of first c2 axis orthogonal to cn main axis
c        is2 = index of first c2 axis orthogonal to sn main axis
c        ich = index of first sigma_h plane
c        icv = index of first sigma_v plane
c        ic3 = index of first c3 axis not parallel to main axis
c =====================================================================
      dimension univec(3,nt),istype(nt)
      logical orthog,parall
c
      do 100 it=1,nt
        istyp=istype(it)
c ---------------------------------------------------------------------
c       c2 axes orthogonal to main axis
c       (including the case that a c3 axis is neither parallel nor
c        orthogonal to a c2 axis)
c ---------------------------------------------------------------------
        if(istyp.eq.2) then
          if(ic2.eq.0.and.icnmax.gt.0.and.
     1       orthog(univec(1,it),univec(1,icnmax),thrsym)) ic2=it
          if(is2.eq.0.and.isnmax.gt.0.and.
     1       orthog(univec(1,it),univec(1,isnmax),thrsym)) is2=it
          if(ic3.eq.0.and.icnmax.eq.3.and.
     1       .not.parall(univec(1,it),univec(1,icnmax),thrsym).and.
     2       .not.orthog(univec(1,it),univec(1,icnmax),thrsym)) ic3=it
c ---------------------------------------------------------------------
c       sigma_v/sigma_h mirror planes
c ---------------------------------------------------------------------
        elseif(istyp.eq.-1) then
          if(icv.eq.0) then
            if(icnmax.gt.0.and.
     1         orthog(univec(1,it),univec(1,icnmax),thrsym)) icv=it
            if(isnmax.gt.0.and.
     1         orthog(univec(1,it),univec(1,isnmax),thrsym)) icv=it
          endif
          if(ich.eq.0) then
            if(icnmax.gt.0.and.
     1         parall(univec(1,it),univec(1,icnmax),thrsym)) ich=it
            if(isnmax.gt.0.and.
     1         parall(univec(1,it),univec(1,isnmax),thrsym)) ich=it
          endif
c ---------------------------------------------------------------------
c       c3 axes not matching main axis
c ---------------------------------------------------------------------
        elseif(istyp.eq.3.and.ic3.eq.0) then
          if(icnmax.gt.0.and.
     1       .not.parall(univec(1,it),univec(1,icnmax),thrsym).and.
     2       .not.orthog(univec(1,it),univec(1,icnmax),thrsym)) ic3=it
          if(isnmax.gt.0.and.
     1       .not.parall(univec(1,it),univec(1,isnmax),thrsym).and.
     2       .not.orthog(univec(1,it),univec(1,isnmax),thrsym)) ic3=it
        endif
  100   continue
cc      write(6,691) ic2,is2,ich,icv,ic3
cc  691 format(/,' c2 axis orthogonal to cn axis         : ',i4,
cc     1       /,' c2 axis orthogonal to sn axis         : ',i4,
cc     2       /,' mirror plane orthogonal to cn/sn axis : ',i4,
cc     3       /,' mirror plane parallel to cn/sn axis   : ',i4,
cc     4       /,' c3 axis different from cn/sn axis     : ',i4,/)
c
      return
      end
c ======================================================================
c
      subroutine syminv(u,nt,invt,thrsym,ierr)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     invt(it) = inverse operation of it
c =====================================================================
      dimension u(9,nt),scr(9),invt(nt)

c
      ierr = -1
c
      do 100 it=1,nt
  100   invt(it)=0
c
      do 200 it=1,nt
        scr(1)=u(1,it)
        scr(2)=u(4,it)
        scr(3)=u(7,it)
        scr(4)=u(2,it)
        scr(5)=u(5,it)
        scr(6)=u(8,it)
        scr(7)=u(3,it)
        scr(8)=u(6,it)
        scr(9)=u(9,it)
        jt=0
  300   jt=jt+1
        if(jt.gt.nt) then
         call message('**WARNING** in SYMMETRY module',
     $                '  Unable to find Inverse Operation',0,0)
         return
        endif
        if(dif(scr,u(1,jt),9).gt.thrsym) goto 300
  200   invt(it)=jt
c
      ierr = 0
c
      return
      end
c ======================================================================
c
      subroutine traxyz(xyz,e1,e2,e3,n)
      implicit real*8 (a-h,o-z)
c =====================================================================
c     transform coordinates according to U = (e1,e2,e3) and write them
c     to file <filnam>, descriptor iunit
c =====================================================================
      dimension xyz(3,n),e1(3),e2(3),e3(3)
c
      do 100 i=1,n
        x=e1(1)*xyz(1,i)+e1(2)*xyz(2,i)+e1(3)*xyz(3,i)
        y=e2(1)*xyz(1,i)+e2(2)*xyz(2,i)+e2(3)*xyz(3,i)
        z=e3(1)*xyz(1,i)+e3(2)*xyz(2,i)+e3(3)*xyz(3,i)
        xyz(1,i)=x
        xyz(2,i)=y
        xyz(3,i)=z
  100   continue
c
      return
      end

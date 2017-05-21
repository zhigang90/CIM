c======================================================================
c
c  MM 2003
c
c  this file contains several numerical quadrature subroutines
c  for computing the dft contributions to the constant part of
c  derivative fock matrices and hessian (without weight derivatives).
c
c  there are two sets of subroutines one for computing
c  both Fock derivatives and hessian at the same time,
c  and one for computing fock derivatives only.
c
c  Set 1: All components of the Hessian and some or all
c         components of the derivative Fock matrices are
c         computed in the same run (this cover Cases 1 and 2
c         of the multiple passes implementation
c
c          RHF case -->  mkmpfh
c          UHF case -->  mkmpfhu
c
c  Set 2: Only components of the derivative Fock matrices
c         are computed (Case 3 of multiple passes)
c
c          RHF case --> mkmpf
c          UHF case --> mkmpfu
c
c =====================================================================
      subroutine mkmpfh(dft,    npp,    nbas,   nbf,    natoms,
     $                  nb,     ne,     nbatm,  thrsh,  da,
     $                  dm,     gden,   wght,   pra,    prara,
     $                  prarb,  pga,    pgc,    praga,  pragb,
     $                  pragc,  pgaga,  pgagb,  pgagc,  pgcgc,
     $                  vao,    vaox,   vaoxx,  vaoxxx, inb,
     $                  vm,     indx,   denx,   denxx,  gdx,
     $                  gdxx,   gx,     gxx,    sv,     sw,
     $                  icntr,  fda,    hess,   td1,    tg1,
     $                  tsw,    tqf,    tqh)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into derivative Fock matrices and Hessian
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C            1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   -  number of contributing (non-zero) grid points this batch
C  NBas  -  total number of basis functions
C  nbf   -  indexing array to location of "non-zero" basis functions
C  NAtoms-  number of atoms
C  Nb,Ne -  First and last component of Fock derivatives to compute
C  NBAtm -  basis functions --> atom index
C  thrsh -  threshold for neglect of contribution
C  DA    -  closed-shell density matrix (lower triangle)
C  DM    -  maximum density matrix element per column
C  GDEN  -  density gradient at grid points (non local only, dft > 3)
C  WGHT  -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density (dft > 3)
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
C  VAO   -  "non-zero" basis function values at grid points
C  VAOX  -  basis function 1st derivatives at grid points
C  VAOXX -  basis function 2nd derivatives at grid points
C  VAOXXX-  basis function 3rd derivatives at grid points (dft > 3)
C  INB   -  indexing array for non-zero entries to VAO
C  VM    -  array containing maximum magnitude AO per grid point
C  INDX  -  index into contributing columns of VAO
C  DenX  -  scratch storage for in situ atomic gradient of the
C           density
C  DenXX -  ditto for atomic Hessian of the density
C  GDX   -  ditto for atomic gradient of density gradient (dft > 3)
C  GDXX  -  ditto for atomic Hessian of density gradient (dft > 3)
C  GX    -  ditto for atomic gradient of gradient invariant (dft > 3)
C  GXX   -  ditto for atomic hessian of gradient invariant (dft > 3)
C  SV    -  storage for coefficient of vao(i)*vao(j) in quadrature
C        -  formula (dft > 3)
C  SW    -  storage for coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  ICntr -  current atomic center
C
C  on exit
C
C  FDA     -  contribution to derivative Fock matrices
C  HESS    -  contribution to Hessian matrix
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),VM(*),
     $          DA(NBas*(NBas+1)/2),DM(NBas)
      DIMENSION DenX(3,NAtoms),DenXX(3,NAtoms,3,NAtoms)
      DIMENSION VAOXXX(10,*)
      DIMENSION GDEN(3,*)
      DIMENSION GDX(3,3,NAtoms),GDXX(3,3,Natoms,3,Natoms),
     $          GX(3,NAtoms),GXX(3,NAtoms,3,NAtoms)
      DIMENSION SV(3,NAtoms),SW(3,3,NAtoms)
      INTEGER   dft,nbf(*),INB(*),INDX(NPP),NBAtm(NBas)
      DIMENSION FDA(3,Nb:Ne,*),HESS(3,NAtoms,3,NAtoms)
      dimension pra(npp),pga(npp),pgc(npp),prara(npp)
      dimension prarb(npp),praga(npp),pragb(npp),pragc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0)
      PARAMETER (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,dogxx,doval
C
C
      NAt3 = 3*NAtoms
c
      IF(dft.LE.3) THEN
C
C  local only
C
      DO 50 IP=1,NPP
      IPP = INDX(IP)
      WG=WGHT(IPP)
      ra = pra(IP)*WG
      rara = prara(IP)*WG
      VMx = VM(IPP)
      CALL ZeroIT(DenX,NAt3)
      CALL ZeroIT(DenXX,NAt3**2)
c
c -- form gradient and Hessian density at current grid point
c
c    NOTE:
c    for the closed shell case, the array DA contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the gradient and hessian
c    of the total closed shell density
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrara=abs(rara)
      dodenx=abrara.gt.epsi
      thrx=unbelpo
      if(dodenx)thrx=thrsh/abrara
      thrx1=thrx/vmx2
      abra=abs(ra)
      dodenxx=abra.gt.epsi
      thrxx=unbelpo
      if(dodenxx)thrxx=thrsh/abra
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        DMx = DM(II)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        ValX = VAOX(1,I)
        ValY = VAOX(2,I)
        ValZ = VAOX(3,I)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx*vmx
        dox=(dodenx.and.(valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        DO 18 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          JAtm = NBAtm(JJ)
          DAIJ = DA(IJ)*VAO(J)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
            DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
            DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
            DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
          endif
          if(doxx.and.abdaijm.gt.thrxx)then
c
c -- atomic Hessian of density
c -- (a) IAtm with IAtm
            DenXX(1,IAtm,1,IAtm)=DenXX(1,IAtm,1,IAtm)+DAIJ*VAOXX(1,I)
            xyc=DAIJ*VAOXX(2,I)
            DenXX(1,IAtm,2,IAtm)=DenXX(1,IAtm,2,IAtm)+xyc
            DenXX(2,IAtm,1,IAtm)=DenXX(2,IAtm,1,IAtm)+xyc
            DenXX(2,IAtm,2,IAtm)=DenXX(2,IAtm,2,IAtm)+DAIJ*VAOXX(3,I)
            xzc=DAIJ*VAOXX(4,I)
            DenXX(1,IAtm,3,IAtm)=DenXX(1,IAtm,3,IAtm)+xzc
            DenXX(3,IAtm,1,IAtm)=DenXX(3,IAtm,1,IAtm)+xzc
            yzc=DAIJ*VAOXX(5,I)
            DenXX(2,IAtm,3,IAtm)=DenXX(2,IAtm,3,IAtm)+yzc
            DenXX(3,IAtm,2,IAtm)=DenXX(3,IAtm,2,IAtm)+yzc
            DenXX(3,IAtm,3,IAtm)=DenXX(3,IAtm,3,IAtm)+DAIJ*VAOXX(6,I)
          endif
          If(JAtm.LT.IAtm) GO TO 18
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
c -- (b) IAtm with JAtm
            DAIJ = DA(IJ)*VAOX(1,J)
            DenXX(1,IAtm,1,JAtm)=DenXX(1,IAtm,1,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,1,JAtm)=DenXX(2,IAtm,1,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,1,JAtm)=DenXX(3,IAtm,1,JAtm)+DAIJ*ValZ
            DAIJ = DA(IJ)*VAOX(2,J)
            DenXX(1,IAtm,2,JAtm)=DenXX(1,IAtm,2,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,2,JAtm)=DenXX(2,IAtm,2,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,2,JAtm)=DenXX(3,IAtm,2,JAtm)+DAIJ*ValZ
            DAIJ = DA(IJ)*VAOX(3,J)
            DenXX(1,IAtm,3,JAtm)=DenXX(1,IAtm,3,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,3,JAtm)=DenXX(2,IAtm,3,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,3,JAtm)=DenXX(3,IAtm,3,JAtm)+DAIJ*ValZ
          endif
 18     CONTINUE
        DO 19 J=I+1,nbf(IPP+1)
          JJ = INB(J)
          IJ = (JJ*(JJ-1))/2 + II
          JAtm = NBAtm(JJ)
          DAIJ = DA(IJ)*VAO(J)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
          if(dox)then
c
c -- atomic gradient of density
            DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
            DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
            DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
          endif
          if(doxx.and.abdaijm.gt.thrxx)then
c
c -- atomic Hessian of density
c -- (a) IAtm with IAtm
            DenXX(1,IAtm,1,IAtm)=DenXX(1,IAtm,1,IAtm)+DAIJ*VAOXX(1,I)
            xyc=DAIJ*VAOXX(2,I)
            DenXX(1,IAtm,2,IAtm)=DenXX(1,IAtm,2,IAtm)+xyc
            DenXX(2,IAtm,1,IAtm)=DenXX(2,IAtm,1,IAtm)+xyc
            DenXX(2,IAtm,2,IAtm)=DenXX(2,IAtm,2,IAtm)+DAIJ*VAOXX(3,I)
            xzc=DAIJ*VAOXX(4,I)
            DenXX(1,IAtm,3,IAtm)=DenXX(1,IAtm,3,IAtm)+xzc
            DenXX(3,IAtm,1,IAtm)=DenXX(3,IAtm,1,IAtm)+xzc
            yzc=DAIJ*VAOXX(5,I)
            DenXX(2,IAtm,3,IAtm)=DenXX(2,IAtm,3,IAtm)+yzc
            DenXX(3,IAtm,2,IAtm)=DenXX(3,IAtm,2,IAtm)+yzc
            DenXX(3,IAtm,3,IAtm)=DenXX(3,IAtm,3,IAtm)+DAIJ*VAOXX(6,I)
          endif
c -- (b) IAtm with JAtm
          If(JAtm.LT.IAtm) GO TO 19
          abdaijm=abs(da(ij))*valm*vmx
          if(doxx.and.abdaijm.gt.thrxx)then
            DAIJ = DA(IJ)*VAOX(1,J)
            DenXX(1,IAtm,1,JAtm)=DenXX(1,IAtm,1,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,1,JAtm)=DenXX(2,IAtm,1,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,1,JAtm)=DenXX(3,IAtm,1,JAtm)+DAIJ*ValZ
            DAIJ = DA(IJ)*VAOX(2,J)
            DenXX(1,IAtm,2,JAtm)=DenXX(1,IAtm,2,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,2,JAtm)=DenXX(2,IAtm,2,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,2,JAtm)=DenXX(3,IAtm,2,JAtm)+DAIJ*ValZ
            DAIJ = DA(IJ)*VAOX(3,J)
            DenXX(1,IAtm,3,JAtm)=DenXX(1,IAtm,3,JAtm)+DAIJ*ValX
            DenXX(2,IAtm,3,JAtm)=DenXX(2,IAtm,3,JAtm)+DAIJ*ValY
            DenXX(3,IAtm,3,JAtm)=DenXX(3,IAtm,3,JAtm)+DAIJ*ValZ
          endif
 19     CONTINUE
 20   CONTINUE
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c    Numerical quadrature for the derivative Fock matrix.
c
c -- get maximum absolute value of 1st-order density at this grid point
      Call absmax(NAt3,DenX,iixx,DMaxyz)
c
c -- global threshold testing
c
      VMax = Max(Abra*VMx2,Abrara*DMaxyz*VMx2)
      If(VMax.LT.thrsh) GO TO 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        Val = VAO(I)*ra
        VDXC = VAO(I)*rara
        ValX = VAOX(1,I)*ra
        ValY = VAOX(2,I)*ra
        ValZ = VAOX(3,I)*ra
        abval= abs(val)
        doval= abval.gt.thrsh1
        valm = max(abs(valx),abs(valy),abs(valz))
        Valt = MAX(abval,valm,Abs(vdxc)*dmaxyz)
        IF(Valt.GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            ValJ = VAO(J)
            if(abs(valj)*valm.gt.thrsh)then
              if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - ValX*ValJ
                FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - ValY*ValJ
                FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - ValZ*ValJ
              endif
            endif
            if(doval)then
              if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - Val*VAOX(1,J)
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - Val*VAOX(2,J)
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - Val*VAOX(3,J)
              endif
            endif
c
c -- contribution of functional derivative
c
            VDJX = VDXC*ValJ
            If(Abs(VDJX)*DMaxyz.GT.thrsh) Then
              Do KAtm=Nb,Ne
                FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ) + VDJX*DENX(1,KAtm)
                FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ) + VDJX*DENX(2,KAtm)
                FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ) + VDJX*DENX(3,KAtm)
              EndDo
            EndIf
 30       CONTINUE
        ENDIF
 40   CONTINUE
c
c    numerical quadrature for the Hessian matrix
c
 45   CONTINUE
        call secund(t3)
        tqf=tqf+t3-t2
c
c -- direct contribution to Hessian matrix.
c    a factor two is applied
c
      rara2=two*rara
      ra2=two*ra
      Do IAtm=1,NAtoms
        HValx=rara2*DenX(1,IAtm)
        HValy=rara2*DenX(2,IAtm)
        HValz=rara2*DenX(3,IAtm)
        Do JAtm=IAtm,NAtoms
          HESS(1,IAtm,1,JAtm) = HESS(1,IAtm,1,JAtm) +
     $     DenX(1,Jatm)*HValx + ra2*DenXX(1,IAtm,1,JAtm)
          HESS(2,IAtm,1,JAtm) = HESS(2,IAtm,1,JAtm) +
     $     DenX(1,Jatm)*HValy + ra2*DenXX(2,IAtm,1,JAtm)
          HESS(3,IAtm,1,JAtm) = HESS(3,IAtm,1,JAtm) +
     $     Denx(1,Jatm)*Hvalz + ra2*DenXX(3,IAtm,1,JAtm)
          HESS(1,IAtm,2,JAtm) = HESS(1,IAtm,2,JAtm) +
     $     DenX(2,JAtm)*Hvalx + ra2*DenXX(1,IAtm,2,JAtm)
          HESS(2,IAtm,2,JAtm) = HESS(2,IAtm,2,JAtm) +
     $     DenX(2,JAtm)*Hvaly + ra2*DenXX(2,IAtm,2,JAtm)
          HESS(3,IAtm,2,JAtm) = HESS(3,IAtm,2,JAtm) +
     $     DenX(2,JAtm)*Hvalz + ra2*DenXX(3,IAtm,2,JAtm)
          HESS(1,IAtm,3,JAtm) = HESS(1,IAtm,3,JAtm) +
     $     DenX(3,JAtm)*Hvalx + ra2*DenXX(1,IAtm,3,JAtm)
          HESS(2,IAtm,3,JAtm) = HESS(2,IAtm,3,JAtm) +
     $     DenX(3,JAtm)*Hvaly + ra2*DenXX(2,IAtm,3,JAtm)
          HESS(3,IAtm,3,JAtm) = HESS(3,IAtm,3,JAtm) +
     $     DenX(3,JAtm)*Hvalz + ra2*DenXX(3,IAtm,3,JAtm)
        EndDo
      EndDo
      call secund(t4)
      tqh=tqh+t4-t3
c
 50   CONTINUE
cc
      ELSE
C
C   non-local dft
C
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        ra = pra(ip)*wg
        ga = pga(ip)*wg
        gc = pgc(ip)*wg
        rara = (prara(ip)+prarb(ip))*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
        prg=raga+ragb+ragc
        pgg=gaga+gagb+Two*gagc+Half*gcgc
        pg=ga+Half*gc
c
c  density gradient at current point
c
        DX=GDEN(1,IPP)
        DY=GDEN(2,IPP)
        DZ=GDEN(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenX,NAt3)
        CALL ZeroIT(DenXX,NAt3**2)
        CALL ZeroIT(GDX,3*NAt3)
        CALL ZeroIT(GDXX,3*NAt3**2)
c
c    form the atomic gradient and Hessian  of density
c    and density gradient at current grid point
c
c    NOTE:
c    for the closed shell case, the array DA contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient and hessian densities is omitted here, so
c    denx and denxx will contain half the atomic gradient and hessian
c    of the total closed shell density. Likewise, gdx and gdxx will
c    contain half the atomic gradient and Hessian of the total closed
c    shell density gradient invariant
c
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          DMx = DM(II)       ! max. element of first order density
          if(dmx.gt.epsi)then
             thtest=thrsh/(vmx*dmx)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJT=DA(IJ)
            abijt=abs(daijt)
            DAIJ = DAIJT*VALJ
            abij=abs(daij)
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
              DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
              DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              GDX(1,1,IAtm)=GDX(1,1,IAtm)-DAIJT*(VALJ*VALXX+VALX*VALJX)
              GDX(1,2,IAtm)=GDX(1,2,IAtm)-DAIJT*(VALJ*VALXY+VALY*VALJX)
              GDX(1,3,IAtm)=GDX(1,3,IAtm)-DAIJT*(VALJ*VALXZ+VALZ*VALJX)
              GDX(2,1,IAtm)=GDX(2,1,IAtm)-DAIJT*(VALJ*VALXY+VALX*VALJY)
              GDX(2,2,IAtm)=GDX(2,2,IAtm)-DAIJT*(VALJ*VALYY+VALY*VALJY)
              GDX(2,3,IAtm)=GDX(2,3,IAtm)-DAIJT*(VALJ*VALYZ+VALZ*VALJY)
              GDX(3,1,IAtm)=GDX(3,1,IAtm)-DAIJT*(VALJ*VALXZ+VALX*VALJZ)
              GDX(3,2,IAtm)=GDX(3,2,IAtm)-DAIJT*(VALJ*VALYZ+VALY*VALJZ)
              GDX(3,3,IAtm)=GDX(3,3,IAtm)-DAIJT*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c          If(IAtm.NE.icntr) then
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXX(1,IAtm,1,IAtm)=DenXX(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXX(1,IAtm,2,IAtm)=DenXX(1,IAtm,2,IAtm)+xyc
              DenXX(2,IAtm,1,IAtm)=DenXX(2,IAtm,1,IAtm)+xyc
              DenXX(2,IAtm,2,IAtm)=DenXX(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXX(1,IAtm,3,IAtm)=DenXX(1,IAtm,3,IAtm)+xzc
              DenXX(3,IAtm,1,IAtm)=DenXX(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXX(2,IAtm,3,IAtm)=DenXX(2,IAtm,3,IAtm)+yzc
              DenXX(3,IAtm,2,IAtm)=DenXX(3,IAtm,2,IAtm)+yzc
              DenXX(3,IAtm,3,IAtm)=DenXX(3,IAtm,3,IAtm)+DAIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
              GDXX(1,1,IAtm,1,IAtm)=GDXX(1,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXX(1,1,IAtm,2,IAtm)=GDXX(1,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXX(1,2,IAtm,1,IAtm)=GDXX(1,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXX(1,2,IAtm,2,IAtm)=GDXX(1,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXX(1,1,IAtm,3,IAtm)=GDXX(1,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXX(1,3,IAtm,1,IAtm)=GDXX(1,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXX(1,2,IAtm,3,IAtm)=GDXX(1,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXX(1,3,IAtm,2,IAtm)=GDXX(1,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXX(1,3,IAtm,3,IAtm)=GDXX(1,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXX(2,1,IAtm,1,IAtm)=GDXX(2,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXX(2,1,IAtm,2,IAtm)=GDXX(2,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXX(2,2,IAtm,1,IAtm)=GDXX(2,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXX(2,2,IAtm,2,IAtm)=GDXX(2,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXX(2,1,IAtm,3,IAtm)=GDXX(2,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXX(2,3,IAtm,1,IAtm)=GDXX(2,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXX(2,2,IAtm,3,IAtm)=GDXX(2,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXX(2,3,IAtm,2,IAtm)=GDXX(2,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXX(2,3,IAtm,3,IAtm)=GDXX(2,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXX(3,1,IAtm,1,IAtm)=GDXX(3,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXX(3,1,IAtm,2,IAtm)=GDXX(3,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXX(3,2,IAtm,1,IAtm)=GDXX(3,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXX(3,2,IAtm,2,IAtm)=GDXX(3,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXX(3,1,IAtm,3,IAtm)=GDXX(3,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXX(3,3,IAtm,1,IAtm)=GDXX(3,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXX(3,2,IAtm,3,IAtm)=GDXX(3,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXX(3,3,IAtm,2,IAtm)=GDXX(3,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXX(3,3,IAtm,3,IAtm)=GDXX(3,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
c
            If(JAtm.LT.IAtm) GO TO 218
c
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJT*VALJX
              DenXX(1,IAtm,1,JAtm)=DenXX(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,1,JAtm)=DenXX(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,1,JAtm)=DenXX(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJT*VALJY
              DenXX(1,IAtm,2,JAtm)=DenXX(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,2,JAtm)=DenXX(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,2,JAtm)=DenXX(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJT*VALJZ
              DenXX(1,IAtm,3,JAtm)=DenXX(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,3,JAtm)=DenXX(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,3,JAtm)=DenXX(3,IAtm,3,JAtm)+DAIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
              GDXX(1,1,IAtm,1,JAtm)=GDXX(1,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXX(1,1,IAtm,2,JAtm)=GDXX(1,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXX(1,2,IAtm,1,JAtm)=GDXX(1,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXX(1,2,IAtm,2,JAtm)=GDXX(1,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXX(1,1,IAtm,3,JAtm)=GDXX(1,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXX(1,3,IAtm,1,JAtm)=GDXX(1,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXX(1,2,IAtm,3,JAtm)=GDXX(1,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXX(1,3,IAtm,2,JAtm)=GDXX(1,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXX(1,3,IAtm,3,JAtm)=GDXX(1,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXX(2,1,IAtm,1,JAtm)=GDXX(2,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXX(2,1,IAtm,2,JAtm)=GDXX(2,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXX(2,2,IAtm,1,JAtm)=GDXX(2,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXX(2,2,IAtm,2,JAtm)=GDXX(2,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXX(2,1,IAtm,3,JAtm)=GDXX(2,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXX(2,3,IAtm,1,JAtm)=GDXX(2,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXX(2,2,IAtm,3,JAtm)=GDXX(2,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXX(2,3,IAtm,2,JAtm)=GDXX(2,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXX(2,3,IAtm,3,JAtm)=GDXX(2,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXX(3,1,IAtm,1,JAtm)=GDXX(3,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXX(3,1,IAtm,2,JAtm)=GDXX(3,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXX(3,2,IAtm,1,JAtm)=GDXX(3,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXX(3,2,IAtm,2,JAtm)=GDXX(3,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXX(3,1,IAtm,3,JAtm)=GDXX(3,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXX(3,3,IAtm,1,JAtm)=GDXX(3,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXX(3,2,IAtm,3,JAtm)=GDXX(3,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXX(3,3,IAtm,2,JAtm)=GDXX(3,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXX(3,3,IAtm,3,JAtm)=GDXX(3,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJT= DA(IJ)
            abijt=abs(daijt)
            DAIJ = DAIJT*VALJ
            abij=abs(daij)
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
              DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
              DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              GDX(1,1,IAtm)=GDX(1,1,IAtm)-DAIJT*(VALJ*VALXX+VALX*VALJX)
              GDX(1,2,IAtm)=GDX(1,2,IAtm)-DAIJT*(VALJ*VALXY+VALY*VALJX)
              GDX(1,3,IAtm)=GDX(1,3,IAtm)-DAIJT*(VALJ*VALXZ+VALZ*VALJX)
              GDX(2,1,IAtm)=GDX(2,1,IAtm)-DAIJT*(VALJ*VALXY+VALX*VALJY)
              GDX(2,2,IAtm)=GDX(2,2,IAtm)-DAIJT*(VALJ*VALYY+VALY*VALJY)
              GDX(2,3,IAtm)=GDX(2,3,IAtm)-DAIJT*(VALJ*VALYZ+VALZ*VALJY)
              GDX(3,1,IAtm)=GDX(3,1,IAtm)-DAIJT*(VALJ*VALXZ+VALX*VALJZ)
              GDX(3,2,IAtm)=GDX(3,2,IAtm)-DAIJT*(VALJ*VALYZ+VALY*VALJZ)
              GDX(3,3,IAtm)=GDX(3,3,IAtm)-DAIJT*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c          If(IAtm.NE.icntr) then
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXX(1,IAtm,1,IAtm)=DenXX(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXX(1,IAtm,2,IAtm)=DenXX(1,IAtm,2,IAtm)+xyc
              DenXX(2,IAtm,1,IAtm)=DenXX(2,IAtm,1,IAtm)+xyc
              DenXX(2,IAtm,2,IAtm)=DenXX(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXX(1,IAtm,3,IAtm)=DenXX(1,IAtm,3,IAtm)+xzc
              DenXX(3,IAtm,1,IAtm)=DenXX(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXX(2,IAtm,3,IAtm)=DenXX(2,IAtm,3,IAtm)+yzc
              DenXX(3,IAtm,2,IAtm)=DenXX(3,IAtm,2,IAtm)+yzc
              DenXX(3,IAtm,3,IAtm)=DenXX(3,IAtm,3,IAtm)+DAIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
              GDXX(1,1,IAtm,1,IAtm)=GDXX(1,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXX(1,1,IAtm,2,IAtm)=GDXX(1,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXX(1,2,IAtm,1,IAtm)=GDXX(1,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXX(1,2,IAtm,2,IAtm)=GDXX(1,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXX(1,1,IAtm,3,IAtm)=GDXX(1,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXX(1,3,IAtm,1,IAtm)=GDXX(1,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXX(1,2,IAtm,3,IAtm)=GDXX(1,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXX(1,3,IAtm,2,IAtm)=GDXX(1,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXX(1,3,IAtm,3,IAtm)=GDXX(1,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXX(2,1,IAtm,1,IAtm)=GDXX(2,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXX(2,1,IAtm,2,IAtm)=GDXX(2,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXX(2,2,IAtm,1,IAtm)=GDXX(2,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXX(2,2,IAtm,2,IAtm)=GDXX(2,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXX(2,1,IAtm,3,IAtm)=GDXX(2,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXX(2,3,IAtm,1,IAtm)=GDXX(2,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXX(2,2,IAtm,3,IAtm)=GDXX(2,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXX(2,3,IAtm,2,IAtm)=GDXX(2,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXX(2,3,IAtm,3,IAtm)=GDXX(2,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXX(3,1,IAtm,1,IAtm)=GDXX(3,1,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXX(3,1,IAtm,2,IAtm)=GDXX(3,1,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXX(3,2,IAtm,1,IAtm)=GDXX(3,2,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXX(3,2,IAtm,2,IAtm)=GDXX(3,2,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXX(3,1,IAtm,3,IAtm)=GDXX(3,1,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXX(3,3,IAtm,1,IAtm)=GDXX(3,3,IAtm,1,IAtm)+
     $                 DAIJT * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXX(3,2,IAtm,3,IAtm)=GDXX(3,2,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXX(3,3,IAtm,2,IAtm)=GDXX(3,3,IAtm,2,IAtm)+
     $                 DAIJT * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXX(3,3,IAtm,3,IAtm)=GDXX(3,3,IAtm,3,IAtm)+
     $                 DAIJT * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
            If(JAtm.LT.IAtm) GO TO 219
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJT*VALJX
              DenXX(1,IAtm,1,JAtm)=DenXX(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,1,JAtm)=DenXX(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,1,JAtm)=DenXX(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJT*VALJY
              DenXX(1,IAtm,2,JAtm)=DenXX(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,2,JAtm)=DenXX(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,2,JAtm)=DenXX(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJT*VALJZ
              DenXX(1,IAtm,3,JAtm)=DenXX(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXX(2,IAtm,3,JAtm)=DenXX(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXX(3,IAtm,3,JAtm)=DenXX(3,IAtm,3,JAtm)+DAIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
              GDXX(1,1,IAtm,1,JAtm)=GDXX(1,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXX(1,1,IAtm,2,JAtm)=GDXX(1,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXX(1,2,IAtm,1,JAtm)=GDXX(1,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXX(1,2,IAtm,2,JAtm)=GDXX(1,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXX(1,1,IAtm,3,JAtm)=GDXX(1,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXX(1,3,IAtm,1,JAtm)=GDXX(1,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXX(1,2,IAtm,3,JAtm)=GDXX(1,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXX(1,3,IAtm,2,JAtm)=GDXX(1,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXX(1,3,IAtm,3,JAtm)=GDXX(1,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXX(2,1,IAtm,1,JAtm)=GDXX(2,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXX(2,1,IAtm,2,JAtm)=GDXX(2,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXX(2,2,IAtm,1,JAtm)=GDXX(2,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXX(2,2,IAtm,2,JAtm)=GDXX(2,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXX(2,1,IAtm,3,JAtm)=GDXX(2,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXX(2,3,IAtm,1,JAtm)=GDXX(2,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXX(2,2,IAtm,3,JAtm)=GDXX(2,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXX(2,3,IAtm,2,JAtm)=GDXX(2,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXX(2,3,IAtm,3,JAtm)=GDXX(2,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXX(3,1,IAtm,1,JAtm)=GDXX(3,1,IAtm,1,JAtm)+
     $                    DAIJT * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXX(3,1,IAtm,2,JAtm)=GDXX(3,1,IAtm,2,JAtm)+
     $                    DAIJT * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXX(3,2,IAtm,1,JAtm)=GDXX(3,2,IAtm,1,JAtm)+
     $                    DAIJT * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXX(3,2,IAtm,2,JAtm)=GDXX(3,2,IAtm,2,JAtm)+
     $                    DAIJT * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXX(3,1,IAtm,3,JAtm)=GDXX(3,1,IAtm,3,JAtm)+
     $                    DAIJT * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXX(3,3,IAtm,1,JAtm)=GDXX(3,3,IAtm,1,JAtm)+
     $                    DAIJT * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXX(3,2,IAtm,3,JAtm)=GDXX(3,2,IAtm,3,JAtm)+
     $                    DAIJT * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXX(3,3,IAtm,2,JAtm)=GDXX(3,3,IAtm,2,JAtm)+
     $                    DAIJT * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXX(3,3,IAtm,3,JAtm)=GDXX(3,3,IAtm,3,JAtm)+
     $                    DAIJT * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and Hessian of density gradient invariant
c
c
        DO IAtm=1,NAtoms
          GX(1,IAtm)=two*(DX*GDX(1,1,IAtm)+
     $                   DY*GDX(2,1,IAtm)+DZ*GDX(3,1,IAtm))
          GX(2,IAtm)=two*(DX*GDX(1,2,IAtm)+
     $                   DY*GDX(2,2,IAtm)+DZ*GDX(3,2,IAtm))
          GX(3,IAtm)=two*(DX*GDX(1,3,IAtm)+
     $                   DY*GDX(2,3,IAtm)+DZ*GDX(3,3,IAtm))
        EndDO
c
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
            GXX(1,IAtm,1,JAtm)=Two*(
     $        DX*GDXX(1,1,IAtm,1,JAtm)+GDX(1,1,IAtm)*GDX(1,1,JAtm)+
     $        DY*GDXX(2,1,IAtm,1,JAtm)+GDX(2,1,IAtm)*GDX(2,1,JAtm)+
     $        DZ*GDXX(3,1,IAtm,1,JAtm)+GDX(3,1,IAtm)*GDX(3,1,JAtm))
            GXX(1,IAtm,2,JAtm)=Two*(
     $        DX*GDXX(1,1,IAtm,2,JAtm)+GDX(1,1,IAtm)*GDX(1,2,JAtm)+
     $        DY*GDXX(2,1,IAtm,2,JAtm)+GDX(2,1,IAtm)*GDX(2,2,JAtm)+
     $        DZ*GDXX(3,1,IAtm,2,JAtm)+GDX(3,1,IAtm)*GDX(3,2,JAtm))
            GXX(2,IAtm,1,JAtm)=Two*(
     $        DX*GDXX(1,2,IAtm,1,JAtm)+GDX(1,2,IAtm)*GDX(1,1,JAtm)+
     $        DY*GDXX(2,2,IAtm,1,JAtm)+GDX(2,2,IAtm)*GDX(2,1,JAtm)+
     $        DZ*GDXX(3,2,IAtm,1,JAtm)+GDX(3,2,IAtm)*GDX(3,1,JAtm))
            GXX(2,IAtm,2,JAtm)=Two*(
     $        DX*GDXX(1,2,IAtm,2,JAtm)+GDX(1,2,IAtm)*GDX(1,2,JAtm)+
     $        DY*GDXX(2,2,IAtm,2,JAtm)+GDX(2,2,IAtm)*GDX(2,2,JAtm)+
     $        DZ*GDXX(3,2,IAtm,2,JAtm)+GDX(3,2,IAtm)*GDX(3,2,JAtm))
            GXX(1,IAtm,3,JAtm)=Two*(
     $        DX*GDXX(1,1,IAtm,3,JAtm)+GDX(1,1,IAtm)*GDX(1,3,JAtm)+
     $        DY*GDXX(2,1,IAtm,3,JAtm)+GDX(2,1,IAtm)*GDX(2,3,JAtm)+
     $        DZ*GDXX(3,1,IAtm,3,JAtm)+GDX(3,1,IAtm)*GDX(3,3,JAtm))
            GXX(3,IAtm,1,JAtm)=Two*(
     $        DX*GDXX(1,3,IAtm,1,JAtm)+GDX(1,3,IAtm)*GDX(1,1,JAtm)+
     $        DY*GDXX(2,3,IAtm,1,JAtm)+GDX(2,3,IAtm)*GDX(2,1,JAtm)+
     $        DZ*GDXX(3,3,IAtm,1,JAtm)+GDX(3,3,IAtm)*GDX(3,1,JAtm))
            GXX(2,IAtm,3,JAtm)=Two*(
     $        DX*GDXX(1,2,IAtm,3,JAtm)+GDX(1,2,IAtm)*GDX(1,3,JAtm)+
     $        DY*GDXX(2,2,IAtm,3,JAtm)+GDX(2,2,IAtm)*GDX(2,3,JAtm)+
     $        DZ*GDXX(3,2,IAtm,3,JAtm)+GDX(3,2,IAtm)*GDX(3,3,JAtm))
            GXX(3,IAtm,2,JAtm)=Two*(
     $        DX*GDXX(1,3,IAtm,2,JAtm)+GDX(1,3,IAtm)*GDX(1,2,JAtm)+
     $        DY*GDXX(2,3,IAtm,2,JAtm)+GDX(2,3,IAtm)*GDX(2,2,JAtm)+
     $        DZ*GDXX(3,3,IAtm,2,JAtm)+GDX(3,3,IAtm)*GDX(3,2,JAtm))
            GXX(3,IAtm,3,JAtm)=Two*(
     $        DX*GDXX(1,3,IAtm,3,JAtm)+GDX(1,3,IAtm)*GDX(1,3,JAtm)+
     $        DY*GDXX(2,3,IAtm,3,JAtm)+GDX(2,3,IAtm)*GDX(2,3,JAtm)+
     $        DZ*GDXX(3,3,IAtm,3,JAtm)+GDX(3,3,IAtm)*GDX(3,3,JAtm))
          EndDO
        EndDO
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
c
        prg2=two*prg
        pgg2=two*pgg
        pg2=two*pg
        SXX=pg2*DX
        SXY=pg2*DY
        SXZ=pg2*DZ
        DO IAtm=1,NAtoms
          vd1=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          vd2=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          vd3=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          sv(1,iatm)=rara*denx(1,iatm)+prg*gx(1,iatm)
          sv(2,iatm)=rara*denx(2,iatm)+prg*gx(2,iatm)
          sv(3,iatm)=rara*denx(3,iatm)+prg*gx(3,iatm)
          SW(1,1,IAtm)=pg2*GDX(1,1,IAtm)+DX*VD1
          SW(2,1,IAtm)=pg2*GDX(2,1,IAtm)+DY*VD1
          SW(3,1,IAtm)=pg2*GDX(3,1,IAtm)+DZ*VD1
          SW(1,2,IAtm)=pg2*GDX(1,2,IAtm)+DX*VD2
          SW(2,2,IAtm)=pg2*GDX(2,2,IAtm)+DY*VD2
          SW(3,2,IAtm)=pg2*GDX(3,2,IAtm)+DZ*VD2
          SW(1,3,IAtm)=pg2*GDX(1,3,IAtm)+DX*VD3
          SW(2,3,IAtm)=pg2*GDX(2,3,IAtm)+DY*VD3
          SW(3,3,IAtm)=pg2*GDX(3,3,IAtm)+DZ*VD3
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    Fock Matrix.
c
c     get the maximum absolute value of coefficients V and W
c
        Call absmax(NAt3,SV,isv,svmax)
        Call absmax(3*NAt3,SW,isw,swmax)
        tswmax=Three*swmax
c -- global threshold testing
        abra=abs(ra)
        vmx2=vmx*vmx
        SMax = Abs(SXX+SXY+SXZ)
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abra+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVI=SXX*ValIX+SXY*ValIY+SXZ*ValIZ+ra*ValI
          SXXI=SXX*ValI
          SXYI=SXY*ValI
          SXZI=SXZ*ValI
          XXVX=SXX*VAOXX(1,I)+SXY*VAOXX(2,I)+SXZ*VAOXX(4,I)+ra*VALIX
          XXVY=SXX*VAOXX(2,I)+SXY*VAOXX(3,I)+SXZ*VAOXX(5,I)+ra*VALIY
          XXVZ=SXX*VAOXX(4,I)+SXY*VAOXX(5,I)+SXZ*VAOXX(6,I)+ra*VALIZ
          valm=abs(vali)
          valm1=abs(xvi)
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1=max(abs(xxvx),abs(xxvy),abs(xxvz))
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJ=SXX*ValJX+SXY*ValJY+SXZ*ValJZ
              if(valmx1*abs(valj)+valmx*abs(xvj).gt.thrsh)then
               if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                 FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - XXVX*ValJ - ValIX*XVJ
                 FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - XXVY*ValJ - ValIY*XVJ
                 FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - XXVZ*ValJ - ValIZ*XVJ
               endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
               if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - XVI*ValJX -
     $          VAOXX(1,J)*SXXI - VAOXX(2,J)*SXYI - VAOXX(4,J)*SXZI
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - XVI*ValJY -
     $          VAOXX(2,J)*SXXI - VAOXX(3,J)*SXYI - VAOXX(5,J)*SXZI
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - XVI*ValJZ -
     $          VAOXX(4,J)*SXXI - VAOXX(5,J)*SXYI - VAOXX(6,J)*SXZI
               endif
              endif
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VALJX*VALI
              VIJY=VALIY*VALJ+VALJY*VALI
              VIJZ=VALIZ*VALJ+VALJZ*VALI
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                Do KAtm=Nb,Ne
                  FDA(1,KAtm,IJ) = FDA(1,KAtm,IJ) + (SV(1,KAtm)*VIJ+
     $            SW(1,1,KAtm)*VIJX+SW(2,1,KAtm)*VIJY+SW(3,1,KAtM)*VIJZ)
                  FDA(2,KAtm,IJ) = FDA(2,KAtm,IJ) + (SV(2,KAtm)*VIJ+
     $            SW(1,2,KAtm)*VIJX+SW(2,2,KAtm)*VIJY+SW(3,2,KAtM)*VIJZ)
                  FDA(3,KAtm,IJ) = FDA(3,KAtm,IJ) + (SV(3,KAtm)*VIJ+
     $            SW(1,3,KAtm)*VIJX+SW(2,3,KAtm)*VIJY+SW(3,3,KatM)*VIJZ)
                EndDo
              EndIf
 230        CONTINUE
          ENDIF
 240    CONTINUE
c
c   Numerical quadrature for direct contribution to Hessian matrix
c
 245    CONTINUE
        call secund(t5)
        tqf=tqf+t5-t4
c
c -- direct contribution to Hessian matrix.
c    a factor two is applied
c
        ra2=two*ra
        rara2=two*rara
        Do IAtm=1,NAtoms
          hdx=rara2*denx(1,iatm)+prg2*gx(1,iatm)
          hdy=rara2*denx(2,iatm)+prg2*gx(2,iatm)
          hdz=rara2*denx(3,iatm)+prg2*gx(3,iatm)
          hgx=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          hgy=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          hgz=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          Do JAtm=IAtm,NAtoms
       HESS(1,IAtm,1,JAtm)=HESS(1,IAtm,1,JAtm)+Hdx*DenX(1,Jatm)+
     $ Hgx*GX(1,JAtm)+ra2*DenXX(1,IAtm,1,JAtm)+pg2*gxx(1,iatm,1,jatm)
       HESS(2,IAtm,1,JAtm)=HESS(2,IAtm,1,JAtm)+Hdy*DenX(1,Jatm)+
     $ Hgy*GX(1,JAtm)+ra2*DenXX(2,IAtm,1,JAtm)+pg2*gxx(2,iatm,1,jatm)
       HESS(3,IAtm,1,JAtm)=HESS(3,IAtm,1,JAtm)+Hdz*DenX(1,Jatm)+
     $ Hgz*GX(1,JAtm)+ra2*DenXX(3,IAtm,1,JAtm)+pg2*gxx(3,iatm,1,jatm)
       HESS(1,IAtm,2,JAtm)=HESS(1,IAtm,2,JAtm)+Hdx*DenX(2,Jatm)+
     $ Hgx*GX(2,JAtm)+ra2*DenXX(1,IAtm,2,JAtm)+pg2*gxx(1,iatm,2,jatm)
       HESS(2,IAtm,2,JAtm)=HESS(2,IAtm,2,JAtm)+Hdy*DenX(2,Jatm)+
     $ Hgy*GX(2,JAtm)+ra2*DenXX(2,IAtm,2,JAtm)+pg2*gxx(2,iatm,2,jatm)
       HESS(3,IAtm,2,JAtm)=HESS(3,IAtm,2,JAtm)+Hdz*DenX(2,Jatm)+
     $ Hgz*GX(2,JAtm)+ra2*DenXX(3,IAtm,2,JAtm)+pg2*gxx(3,iatm,2,jatm)
       HESS(1,IAtm,3,JAtm)=HESS(1,IAtm,3,JAtm)+Hdx*DenX(3,Jatm)+
     $ Hgx*GX(3,JAtm)+ra2*DenXX(1,IAtm,3,JAtm)+pg2*gxx(1,iatm,3,jatm)
       HESS(2,IAtm,3,JAtm)=HESS(2,IAtm,3,JAtm)+Hdy*DenX(3,Jatm)+
     $ Hgy*GX(3,JAtm)+ra2*DenXX(2,IAtm,3,JAtm)+pg2*gxx(2,iatm,3,jatm)
       HESS(3,IAtm,3,JAtm)=HESS(3,IAtm,3,JAtm)+Hdz*DenX(3,Jatm)+
     $ Hgz*GX(3,JAtm)+ra2*DenXX(3,IAtm,3,JAtm)+pg2*gxx(3,iatm,3,jatm)
          EndDo
        EndDo
        call secund(t6)
        tqh=tqh+t6-t5
c
 200  CONTINUE
      ENDIF
C
      RETURN
      END
c =====================================================================
      subroutine mkmpfhu(dft,    npp,    nbas,   nbf,    natoms,
     $                   nb,     ne,     nbatm,  thrsh,  da,
     $                   db,     dm,     gdena,  gdenb,  wght,
     $                   pra,    prb,    prara,  prbrb,  prarb,
     $                   pga,    pgb,    pgc,    praga,  pragb,
     $                   pragc,  prbga,  prbgb,  prbgc,  pgaga,
     $                   pgagb,  pgagc,  pgbgb,  pgbgc,  pgcgc,
     $                   vao,    vaox,   vaoxx,  vaoxxx, inb,
     $                   vm,     indx,   denxa,  denxb,  denxxa,
     $                   denxxb, gdxa,   gdxb,   gdxxa,  gdxxb,
     $                   gxaa,   gxbb,   gxab,   gxxaa,  gxxbb,
     $                   gxxab,  sva,    svb,    swa,    swb,
     $                   sga,    sgb,    sgc,    icntr,  fda,
     $                   fdb,    hess,   td1,    tg1,    tsw,
     $                   tqf,    tqh)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into derivative Fock matrices and Hessian
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C            1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   -  number of contributing (non-zero) grid points this batch
C  NBas  -  total number of basis functions
C  nbf   -  indexing array to location of "non-zero" basis functions
C  NAtoms-  number of atoms
C  Nb,Ne -  First and last component of Fock derivatives to compute
C  NBAtm -  basis functions --> atom index
C  thrsh -  threshold for neglect of contribution
C  DA    -  alpha density matrix (lower triangle)
C  DB    -  beta density matrix (lower triangle)
C  DM    -  maximum density matrix element per column
C  GDENA -  alpha density gradient at grid points (non local, dft > 3)
C  GDENB -  beta density gradient at grid points (non local, dft > 3)
C  WGHT  -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prb   -  Functional derivative w.r.t. beta density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgb   -  Funct. deriv. w.r.t. beta gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  prbga -  Funct. 2nd. deriv. w.r.t. beta dens. and alpha grad.
C           (dft > 3)
C  prbgb -  Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (dft > 3)
C  prbgc -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgbgb -  Funct. 2nd. deriv. w.r.t. beta grad. (dft > 3)
C  pgbgc -  Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
C  VAO   -  "non-zero" basis function values at grid points
C  VAOX  -  basis function 1st derivatives at grid points
C  VAOXX -  basis function 2nd derivatives at grid points
C  VAOXXX-  basis function 3rd derivatives at grid points (dft > 3)
C  INB   -  indexing array for non-zero entries to VAO
C  VM    -  array containing maximum magnitude AO per grid point
C  INDX  -  index into contributing columns of VAO
C  DenXA -  scratch storage for in situ atomic gradient of the
C           alpha density
C  DenXB -  ditto for atomic gradient of beta density
C  DenXXA-  ditto for atomic Hessian of alpha density
C  DenXXB-  ditto for atomic Hessian of beta density
C  GDXA  -  ditto for atomic gradient of alpha dens. gradient (dft > 3)
C  GDXB  -  ditto for atomic gradient of beta dens. gradient (dft > 3)
C  GDXXA -  ditto for atomic Hessian of alpha density gradient (dft > 3)
C  GDXXB -  ditto for atomic Hessian of beta density gradient (dft > 3)
C  GXAA  -  ditto for atomic gradient of alpha grad. invariant (dft > 3)
C  GXBB  -  ditto for atomic gradient of beta grad. invariant (dft > 3)
C  GXAB  -  ditto for atomic gradient of alpha beta grad. invariant
C           (dft > 3)
C  GXXAA -  ditto for atomic hessian of alpha grad. invariant (dft > 3)
C  GXXBB -  ditto for atomic hessian of beta grad. invariant (dft > 3)
C  GXXAB -  ditto for atomic hessian of alpha beta grad. invariant
C           (dft > 3)
C  SVA   -  storage for alpha coefficient of vao(i)*vao(j) in quadrature
C  SVB   -  storage for beta coefficient of vao(i)*vao(j) in quadrature
C  SWA   -  storage for alpha coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  SWB   -  storage for beta coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  SGA   -  storage for partial sum for hessian quadrature (dft > 3)
C  SGB   -  storage for partial sum for hessian quadrature (dft > 3)
C  SGC   -  storage for partial sum for hessian quadrature (dft > 3)
C  ICntr -  current atomic center
C
C  on exit
C
C  FDA     -  contribution to alpha derivative Fock matrices
C  FDB     -  contribution to beta derivative Fock matrices
C  HESS    -  contribution to Hessian matrix
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),VM(*),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas)
      DIMENSION DenXA(3,NAtoms),DenXXA(3,NAtoms,3,NAtoms)
      DIMENSION DenXB(3,NAtoms),DenXXB(3,NAtoms,3,NAtoms)
      DIMENSION VAOXXX(10,*)
      DIMENSION GDENA(3,*),GDENB(3,*)
      DIMENSION GDXA(3,3,NAtoms),GDXXA(3,3,Natoms,3,Natoms)
      DIMENSION GDXB(3,3,NAtoms),GDXXB(3,3,Natoms,3,Natoms)
      DIMENSION GXAA(3,NAtoms),GXXAA(3,NAtoms,3,NAtoms)
      DIMENSION GXBB(3,NAtoms),GXXBB(3,NAtoms,3,NAtoms)
      DIMENSION GXAB(3,NAtoms),GXXAB(3,NAtoms,3,NAtoms)
      DIMENSION SVA(3,NAtoms),SWA(3,3,NAtoms)
      DIMENSION SVB(3,NAtoms),SWB(3,3,NAtoms)
      DIMENSION sga(3,NAtoms),sgb(3,NAtoms),sgc(3,NAtoms)
      INTEGER   dft,nbf(*),INB(*),INDX(NPP),NBAtm(NBas)
      DIMENSION FDA(3,Nb:Ne,*),FDB(3,Nb:Ne,*),HESS(3,NAtoms,3,NAtoms)
      dimension pra(npp),prb(npp),pga(npp),pgb(npp),pgc(npp)
      dimension prara(npp),prbrb(npp),prarb(npp)
      dimension praga(npp),pragb(npp),pragc(npp)
      dimension prbga(npp),prbgb(npp),prbgc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgbgb(npp),pgbgc(npp),
     $          pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0)
      PARAMETER (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,dogxx,doval
C
C
      NAt3 = 3*NAtoms
c
      IF(dft.LE.3) THEN
C
C  local only
C
      DO 50 IP=1,NPP
      IPP = INDX(IP)
      WG=WGHT(IPP)
      ra = pra(IP)*WG
      rb = prb(IP)*WG
      rara = prara(IP)*WG
      rbrb = prbrb(IP)*WG
      rarb = prarb(IP)*WG
      VMx = VM(IPP)
      CALL ZeroIT(DenXA,NAt3)
      CALL ZeroIT(DenXB,NAt3)
      CALL ZeroIT(DenXXA,NAt3**2)
      CALL ZeroIT(DenXXB,NAt3**2)
c
c -- form gradient and Hessian density at current grid point
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrr=max(abs(rara),abs(rbrb))+abs(rarb)
      dodenx=abrr.gt.epsi
      thrx=unbelpo
      if(dodenx)thrx=thrsh/abrr
      thrx1=thrx/vmx2
      abr=max(abs(ra),abs(rb))
      dodenxx=abr.gt.epsi
      thrxx=unbelpo
      if(dodenxx)thrxx=thrsh/abr
      if(.not.(dodenx.or.dodenxx))goto 25 ! should never happen!
c
c   now compute density derivatives
c
      DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        DMx2 = DM(II)*two
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        ValX = VAOX(1,I)
        ValY = VAOX(2,I)
        ValZ = VAOX(3,I)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx2*vmx
        dox=(dodenx.and.(two*valt.gt.thrx1.or.valt**2.gt.thrx))
        doxx=(dodenxx.and.dmx2*vmx2.gt.thrxx)
        if(.not.(dox.or.doxx)) goto 20
        DO 18 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          JAtm = NBAtm(JJ)
          DAIJ2=two*DA(IJ)
          DBIJ2=two*DB(IJ)
          DAIJ = DAIJ2*VAO(J)
          DBIJ = DBIJ2*VAO(J)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
            DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
            DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
            DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
            DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
            DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
            DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
          endif
          if(doxx.and.abdijm.gt.thrxx)then
c
c -- atomic Hessian of density
c -- (a) IAtm with IAtm
            DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VAOXX(1,I)
            xyc=DAIJ*VAOXX(2,I)
            DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
            DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
            DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VAOXX(3,I)
            xzc=DAIJ*VAOXX(4,I)
            DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
            DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
            yzc=DAIJ*VAOXX(5,I)
            DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
            DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
            DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VAOXX(6,I)
c
            DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VAOXX(1,I)
            xyc=DBIJ*VAOXX(2,I)
            DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
            DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
            DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VAOXX(3,I)
            xzc=DBIJ*VAOXX(4,I)
            DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
            DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
            yzc=DBIJ*VAOXX(5,I)
            DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
            DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
            DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VAOXX(6,I)
          endif
          If(JAtm.LT.IAtm) GO TO 18
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
c -- (b) IAtm with JAtm
            DAIJ = DAIJ2*VAOX(1,J)
            DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
            DAIJ = DAIJ2*VAOX(2,J)
            DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
            DAIJ = DAIJ2*VAOX(3,J)
            DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
            DBIJ = DBIJ2*VAOX(1,J)
            DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
            DBIJ = DBIJ2*VAOX(2,J)
            DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
            DBIJ = DBIJ2*VAOX(3,J)
            DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
          endif
 18     CONTINUE
        DO 19 J=I+1,nbf(IPP+1)
          JJ = INB(J)
          IJ = (JJ*(JJ-1))/2 + II
          JAtm = NBAtm(JJ)
          DAIJ2 = two*DA(IJ)
          DBIJ2 = two*DB(IJ)
          DAIJ = DAIJ2*VAO(J)
          DBIJ = DBIJ2*VAO(J)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
          if(dox)then
c
c -- atomic gradient of density
            DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
            DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
            DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
            DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
            DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
            DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
          endif
          if(doxx.and.abdijm.gt.thrxx)then
c
c -- atomic Hessian of density
c -- (a) IAtm with IAtm
            DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VAOXX(1,I)
            xyc=DAIJ*VAOXX(2,I)
            DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
            DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
            DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VAOXX(3,I)
            xzc=DAIJ*VAOXX(4,I)
            DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
            DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
            yzc=DAIJ*VAOXX(5,I)
            DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
            DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
            DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VAOXX(6,I)
c
            DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VAOXX(1,I)
            xyc=DBIJ*VAOXX(2,I)
            DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
            DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
            DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VAOXX(3,I)
            xzc=DBIJ*VAOXX(4,I)
            DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
            DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
            yzc=DBIJ*VAOXX(5,I)
            DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
            DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
            DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VAOXX(6,I)
          endif
c -- (b) IAtm with JAtm
          If(JAtm.LT.IAtm) GO TO 19
          abdijm=max(abs(da(ij)),abs(db(ij)))*valm*vmx
          if(doxx.and.abdijm.gt.thrxx)then
            DAIJ = DAIJ2*VAOX(1,J)
            DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
            DAIJ = DAIJ2*VAOX(2,J)
            DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
            DAIJ = DAIJ2*VAOX(3,J)
            DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
            DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
            DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
            DBIJ = DBIJ2*VAOX(1,J)
            DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
            DBIJ = DBIJ2*VAOX(2,J)
            DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
            DBIJ = DBIJ2*VAOX(3,J)
            DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
            DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
            DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
          endif
 19     CONTINUE
 20   CONTINUE
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7436 of Johnson and Frisch).
c
        DO IAtm=1,NAtoms
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)
          svb(1,iatm)=rbrb*denxb(1,iatm)+rarb*denxa(1,iatm)
          svb(2,iatm)=rbrb*denxb(2,iatm)+rarb*denxa(2,iatm)
          svb(3,iatm)=rbrb*denxb(3,iatm)+rarb*denxa(3,iatm)
        EndDO
c
c    Numerical quadrature for the derivative Fock matrix.
c
c -- get maximum absolute value of sva and svb
      Call absmax(NAt3,sva,iixx,svamax)
      Call absmax(NAt3,svb,iixx,svbmax)
      svmax=max(svamax,svbmax)
c
c -- global threshold testing
c
      VMax = Max(Abr*VMx2,Abrr*svmax*VMx2)
      If(VMax.LT.thrsh) GO TO 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        ValI = VAO(I)
        Vala = ValI*ra
        Valb = ValI*rb
        ValXa = VAOX(1,I)*ra
        ValYa = VAOX(2,I)*ra
        ValZa = VAOX(3,I)*ra
        ValXb = VAOX(1,I)*rb
        ValYb = VAOX(2,I)*rb
        ValZb = VAOX(3,I)*rb
        abval= max(abs(vala),abs(valb))
        doval= abval.gt.thrsh1
        valm = max(abs(vaox(1,I)),abs(vaox(2,I)),abs(vaox(3,I)))*abr
        Valt = MAX(abval,valm,Abs(vali)*svmax)
        IF(Valt.GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            ValJ = VAO(J)
            if(abs(valj)*valm.gt.thrsh)then
              if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - ValXa*ValJ
                FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - ValYa*ValJ
                FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - ValZa*ValJ
                FDB(1,IAtm,IJ) = FDB(1,IAtm,IJ) - ValXb*ValJ
                FDB(2,IAtm,IJ) = FDB(2,IAtm,IJ) - ValYb*ValJ
                FDB(3,IAtm,IJ) = FDB(3,IAtm,IJ) - ValZb*ValJ
              endif
            endif
            if(doval)then
              if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - Vala*VAOX(1,J)
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - Vala*VAOX(2,J)
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - Vala*VAOX(3,J)
                FDB(1,JAtm,IJ) = FDB(1,JAtm,IJ) - Valb*VAOX(1,J)
                FDB(2,JAtm,IJ) = FDB(2,JAtm,IJ) - Valb*VAOX(2,J)
                FDB(3,JAtm,IJ) = FDB(3,JAtm,IJ) - Valb*VAOX(3,J)
              endif
            endif
c
c -- contribution of functional derivative
c
            ValIJ = ValI*ValJ
            If(Abs(ValIJ)*svmax.GT.thrsh) Then
              Do KAtm=Nb,Ne
                FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ) + ValIJ*sva(1,KAtm)
                FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ) + ValIJ*sva(2,KAtm)
                FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ) + ValIJ*sva(3,KAtm)
                FDB(1,KAtm,IJ)=FDB(1,KAtm,IJ) + ValIJ*svb(1,KAtm)
                FDB(2,KAtm,IJ)=FDB(2,KAtm,IJ) + ValIJ*svb(2,KAtm)
                FDB(3,KAtm,IJ)=FDB(3,KAtm,IJ) + ValIJ*svb(3,KAtm)
              EndDo
            EndIf
 30       CONTINUE
        ENDIF
 40   CONTINUE
c
c    numerical quadrature for the Hessian matrix
c
 45   CONTINUE
        call secund(t3)
        tqf=tqf+t3-t2
c
c -- direct contribution to Hessian matrix.
c
      Do IAtm=1,NAtoms
        Do JAtm=IAtm,NAtoms
          HESS(1,IAtm,1,JAtm) = HESS(1,IAtm,1,JAtm) +
     $     sva(1,IAtm)*DenXA(1,Jatm) + svb(1,Iatm)*DenXB(1,Jatm)+
     $     ra*DenXXA(1,IAtm,1,JAtm) + rb*DenXXB(1,IAtm,1,Jatm)
          HESS(2,IAtm,1,JAtm) = HESS(2,IAtm,1,JAtm) +
     $     sva(2,IAtm)*DenXA(1,Jatm) + svb(2,Iatm)*DenXB(1,Jatm)+
     $     ra*DenXXA(2,IAtm,1,JAtm) + rb*DenXXB(2,IAtm,1,Jatm)
          HESS(3,IAtm,1,JAtm) = HESS(3,IAtm,1,JAtm) +
     $     sva(3,IAtm)*DenXA(1,Jatm) + svb(3,Iatm)*DenXB(1,Jatm)+
     $     ra*DenXXA(3,IAtm,1,JAtm) + rb*DenXXB(3,IAtm,1,Jatm)
          HESS(1,IAtm,2,JAtm) = HESS(1,IAtm,2,JAtm) +
     $     sva(1,IAtm)*DenXA(2,Jatm) + svb(1,Iatm)*DenXB(2,Jatm)+
     $     ra*DenXXA(1,IAtm,2,JAtm) + rb*DenXXB(1,IAtm,2,Jatm)
          HESS(2,IAtm,2,JAtm) = HESS(2,IAtm,2,JAtm) +
     $     sva(2,IAtm)*DenXA(2,Jatm) + svb(2,Iatm)*DenXB(2,Jatm)+
     $     ra*DenXXA(2,IAtm,2,JAtm) + rb*DenXXB(2,IAtm,2,Jatm)
          HESS(3,IAtm,2,JAtm) = HESS(3,IAtm,2,JAtm) +
     $     sva(3,IAtm)*DenXA(2,Jatm) + svb(3,Iatm)*DenXB(2,Jatm)+
     $     ra*DenXXA(3,IAtm,2,JAtm) + rb*DenXXB(3,IAtm,2,Jatm)
          HESS(1,IAtm,3,JAtm) = HESS(1,IAtm,3,JAtm) +
     $     sva(1,IAtm)*DenXA(3,Jatm) + svb(1,Iatm)*DenXB(3,Jatm)+
     $     ra*DenXXA(1,IAtm,3,JAtm) + rb*DenXXB(1,IAtm,3,Jatm)
          HESS(2,IAtm,3,JAtm) = HESS(2,IAtm,3,JAtm) +
     $     sva(2,IAtm)*DenXA(3,Jatm) + svb(2,Iatm)*DenXB(3,Jatm)+
     $     ra*DenXXA(2,IAtm,3,JAtm) + rb*DenXXB(2,IAtm,3,Jatm)
          HESS(3,IAtm,3,JAtm) = HESS(3,IAtm,3,JAtm) +
     $     sva(3,IAtm)*DenXA(3,Jatm) + svb(3,Iatm)*DenXB(3,Jatm)+
     $     ra*DenXXA(3,IAtm,3,JAtm) + rb*DenXXB(3,IAtm,3,Jatm)
        EndDo
      EndDo
      call secund(t4)
      tqh=tqh+t4-t3
c
 50   CONTINUE
cc
      ELSE
C
C   non-local dft
C
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        ra = pra(ip)*wg
        rb = prb(ip)*wg
        ga = pga(ip)*wg
        gb = pgb(ip)*wg
        gc = pgc(ip)*wg
        rara = prara(ip)*wg
        rbrb = prbrb(ip)*wg
        rarb = prarb(ip)*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        rbga = prbga(ip)*wg
        rbgb = prbgb(ip)*wg
        rbgc = prbgc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gbgb = pgbgb(ip)*wg
        gbgc = pgbgc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
        DAX=GDENA(1,IPP)
        DAY=GDENA(2,IPP)
        DAZ=GDENA(3,IPP)
        DBX=GDENB(1,IPP)
        DBY=GDENB(2,IPP)
        DBZ=GDENB(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenXA,NAt3)
        CALL ZeroIT(DenXXA,NAt3**2)
        CALL ZeroIT(GDXA,3*NAt3)
        CALL ZeroIT(GDXXA,3*NAt3**2)
        CALL ZeroIT(DenXB,NAt3)
        CALL ZeroIT(DenXXB,NAt3**2)
        CALL ZeroIT(GDXB,3*NAt3)
        CALL ZeroIT(GDXXB,3*NAt3**2)
c
c    form the atomic gradient and Hessian  of density
c    and density gradient at current grid point
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we check whether the second derivatives of the
c    density gradient should be computed at all.
c
        abgmax=max(abs(ga),abs(gb),abs(gc))
        dogradxx=abgmax.gt.epsi
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          DMx2 = DM(II)*two   ! max. element of first order density
          if(dmx2.gt.epsi)then
             thtest=thrsh/(vmx*dmx2)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2=two*DA(IJ)
            DBIJ2=two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c  alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c  beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
c
            If(JAtm.LT.IAtm) GO TO 218
c
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c    beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            JAtm = NBAtm(JJ)
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2= two*DA(IJ)
            DBIJ2= two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
c
c -- (a) one center terms: IAtm with IAtm
c -- atomic Hessian of density
c
            if(abij*valmm.gt.thrsh)then
              DenXXA(1,IAtm,1,IAtm)=DenXXA(1,IAtm,1,IAtm)+DAIJ*VALXX
              xyc=DAIJ*VALXY
              DenXXA(1,IAtm,2,IAtm)=DenXXA(1,IAtm,2,IAtm)+xyc
              DenXXA(2,IAtm,1,IAtm)=DenXXA(2,IAtm,1,IAtm)+xyc
              DenXXA(2,IAtm,2,IAtm)=DenXXA(2,IAtm,2,IAtm)+DAIJ*VALYY
              xzc=DAIJ*VALXZ
              DenXXA(1,IAtm,3,IAtm)=DenXXA(1,IAtm,3,IAtm)+xzc
              DenXXA(3,IAtm,1,IAtm)=DenXXA(3,IAtm,1,IAtm)+xzc
              yzc=DAIJ*VALYZ
              DenXXA(2,IAtm,3,IAtm)=DenXXA(2,IAtm,3,IAtm)+yzc
              DenXXA(3,IAtm,2,IAtm)=DenXXA(3,IAtm,2,IAtm)+yzc
              DenXXA(3,IAtm,3,IAtm)=DenXXA(3,IAtm,3,IAtm)+DAIJ*VALZZ
c
              DenXXB(1,IAtm,1,IAtm)=DenXXB(1,IAtm,1,IAtm)+DBIJ*VALXX
              xyc=DBIJ*VALXY
              DenXXB(1,IAtm,2,IAtm)=DenXXB(1,IAtm,2,IAtm)+xyc
              DenXXB(2,IAtm,1,IAtm)=DenXXB(2,IAtm,1,IAtm)+xyc
              DenXXB(2,IAtm,2,IAtm)=DenXXB(2,IAtm,2,IAtm)+DBIJ*VALYY
              xzc=DBIJ*VALXZ
              DenXXB(1,IAtm,3,IAtm)=DenXXB(1,IAtm,3,IAtm)+xzc
              DenXXB(3,IAtm,1,IAtm)=DenXXB(3,IAtm,1,IAtm)+xzc
              yzc=DBIJ*VALYZ
              DenXXB(2,IAtm,3,IAtm)=DenXXB(2,IAtm,3,IAtm)+yzc
              DenXXB(3,IAtm,2,IAtm)=DenXXB(3,IAtm,2,IAtm)+yzc
              DenXXB(3,IAtm,3,IAtm)=DenXXB(3,IAtm,3,IAtm)+DBIJ*VALZZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvj*vmx+abvvj*valmm).gt.thrsh)then
c    alpha
              GDXXA(1,1,IAtm,1,IAtm)=GDXXA(1,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXA(1,1,IAtm,2,IAtm)=GDXXA(1,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,1,IAtm)=GDXXA(1,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXA(1,2,IAtm,2,IAtm)=GDXXA(1,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXA(1,1,IAtm,3,IAtm)=GDXXA(1,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,3,IAtm,1,IAtm)=GDXXA(1,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXA(1,2,IAtm,3,IAtm)=GDXXA(1,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,2,IAtm)=GDXXA(1,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXA(1,3,IAtm,3,IAtm)=GDXXA(1,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXA(2,1,IAtm,1,IAtm)=GDXXA(2,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXA(2,1,IAtm,2,IAtm)=GDXXA(2,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,1,IAtm)=GDXXA(2,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXA(2,2,IAtm,2,IAtm)=GDXXA(2,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXA(2,1,IAtm,3,IAtm)=GDXXA(2,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,3,IAtm,1,IAtm)=GDXXA(2,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXA(2,2,IAtm,3,IAtm)=GDXXA(2,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,2,IAtm)=GDXXA(2,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXA(2,3,IAtm,3,IAtm)=GDXXA(2,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXA(3,1,IAtm,1,IAtm)=GDXXA(3,1,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXA(3,1,IAtm,2,IAtm)=GDXXA(3,1,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,1,IAtm)=GDXXA(3,2,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXA(3,2,IAtm,2,IAtm)=GDXXA(3,2,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXA(3,1,IAtm,3,IAtm)=GDXXA(3,1,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,3,IAtm,1,IAtm)=GDXXA(3,3,IAtm,1,IAtm)+
     $                 DAIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXA(3,2,IAtm,3,IAtm)=GDXXA(3,2,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,2,IAtm)=GDXXA(3,3,IAtm,2,IAtm)+
     $                 DAIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXA(3,3,IAtm,3,IAtm)=GDXXA(3,3,IAtm,3,IAtm)+
     $                 DAIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
c    beta
              GDXXB(1,1,IAtm,1,IAtm)=GDXXB(1,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(1,I)*VALJ+VALXX*VALJX)
              GDXXB(1,1,IAtm,2,IAtm)=GDXXB(1,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,1,IAtm)=GDXXB(1,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXY*VALJX)
              GDXXB(1,2,IAtm,2,IAtm)=GDXXB(1,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALYY*VALJX)
              GDXXB(1,1,IAtm,3,IAtm)=GDXXB(1,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,3,IAtm,1,IAtm)=GDXXB(1,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXZ*VALJX)
              GDXXB(1,2,IAtm,3,IAtm)=GDXXB(1,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,2,IAtm)=GDXXB(1,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALYZ*VALJX)
              GDXXB(1,3,IAtm,3,IAtm)=GDXXB(1,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALZZ*VALJX)
c
              GDXXB(2,1,IAtm,1,IAtm)=GDXXB(2,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(2,I)*VALJ+VALXX*VALJY)
              GDXXB(2,1,IAtm,2,IAtm)=GDXXB(2,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,1,IAtm)=GDXXB(2,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(3,I)*VALJ+VALXY*VALJY)
              GDXXB(2,2,IAtm,2,IAtm)=GDXXB(2,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(4,I)*VALJ+VALYY*VALJY)
              GDXXB(2,1,IAtm,3,IAtm)=GDXXB(2,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,3,IAtm,1,IAtm)=GDXXB(2,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXZ*VALJY)
              GDXXB(2,2,IAtm,3,IAtm)=GDXXB(2,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,2,IAtm)=GDXXB(2,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYZ*VALJY)
              GDXXB(2,3,IAtm,3,IAtm)=GDXXB(2,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALZZ*VALJY)
c
              GDXXB(3,1,IAtm,1,IAtm)=GDXXB(3,1,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(5,I)*VALJ+VALXX*VALJZ)
              GDXXB(3,1,IAtm,2,IAtm)=GDXXB(3,1,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,1,IAtm)=GDXXB(3,2,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(6,I)*VALJ+VALXY*VALJZ)
              GDXXB(3,2,IAtm,2,IAtm)=GDXXB(3,2,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(7,I)*VALJ+VALYY*VALJZ)
              GDXXB(3,1,IAtm,3,IAtm)=GDXXB(3,1,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,3,IAtm,1,IAtm)=GDXXB(3,3,IAtm,1,IAtm)+
     $                 DBIJ2 * (VAOXXX(8,I)*VALJ+VALXZ*VALJZ)
              GDXXB(3,2,IAtm,3,IAtm)=GDXXB(3,2,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,2,IAtm)=GDXXB(3,3,IAtm,2,IAtm)+
     $                 DBIJ2 * (VAOXXX(9,I)*VALJ+VALYZ*VALJZ)
              GDXXB(3,3,IAtm,3,IAtm)=GDXXB(3,3,IAtm,3,IAtm)+
     $                 DBIJ2 * (VAOXXX(10,I)*VALJ+VALZZ*VALJZ)
            endif
c
c -- (b) Two center terms: IAtm with JAtm
            If(JAtm.LT.IAtm) GO TO 219
c -- atomic Hessian of density
c
            if(abijt*abvvj*valm.gt.thrsh)then
              DAIJ = DAIJ2*VALJX
              DenXXA(1,IAtm,1,JAtm)=DenXXA(1,IAtm,1,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,1,JAtm)=DenXXA(2,IAtm,1,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,1,JAtm)=DenXXA(3,IAtm,1,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJY
              DenXXA(1,IAtm,2,JAtm)=DenXXA(1,IAtm,2,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,2,JAtm)=DenXXA(2,IAtm,2,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,2,JAtm)=DenXXA(3,IAtm,2,JAtm)+DAIJ*ValZ
              DAIJ = DAIJ2*VALJZ
              DenXXA(1,IAtm,3,JAtm)=DenXXA(1,IAtm,3,JAtm)+DAIJ*ValX
              DenXXA(2,IAtm,3,JAtm)=DenXXA(2,IAtm,3,JAtm)+DAIJ*ValY
              DenXXA(3,IAtm,3,JAtm)=DenXXA(3,IAtm,3,JAtm)+DAIJ*ValZ
c
              DBIJ = DBIJ2*VALJX
              DenXXB(1,IAtm,1,JAtm)=DenXXB(1,IAtm,1,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,1,JAtm)=DenXXB(2,IAtm,1,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,1,JAtm)=DenXXB(3,IAtm,1,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJY
              DenXXB(1,IAtm,2,JAtm)=DenXXB(1,IAtm,2,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,2,JAtm)=DenXXB(2,IAtm,2,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,2,JAtm)=DenXXB(3,IAtm,2,JAtm)+DBIJ*ValZ
              DBIJ = DBIJ2*VALJZ
              DenXXB(1,IAtm,3,JAtm)=DenXXB(1,IAtm,3,JAtm)+DBIJ*ValX
              DenXXB(2,IAtm,3,JAtm)=DenXXB(2,IAtm,3,JAtm)+DBIJ*ValY
              DenXXB(3,IAtm,3,JAtm)=DenXXB(3,IAtm,3,JAtm)+DBIJ*ValZ
            endif
c
c -- atomic Hessian of density gradient
c
            if(dogradxx.and.abijt*(abvvj*valmm+vmx*valm).gt.thrsh)then
c     alpha
              GDXXA(1,1,IAtm,1,JAtm)=GDXXA(1,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXA(1,1,IAtm,2,JAtm)=GDXXA(1,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXA(1,2,IAtm,1,JAtm)=GDXXA(1,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXA(1,2,IAtm,2,JAtm)=GDXXA(1,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXA(1,1,IAtm,3,JAtm)=GDXXA(1,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXA(1,3,IAtm,1,JAtm)=GDXXA(1,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXA(1,2,IAtm,3,JAtm)=GDXXA(1,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXA(1,3,IAtm,2,JAtm)=GDXXA(1,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXA(1,3,IAtm,3,JAtm)=GDXXA(1,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXA(2,1,IAtm,1,JAtm)=GDXXA(2,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXA(2,1,IAtm,2,JAtm)=GDXXA(2,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXA(2,2,IAtm,1,JAtm)=GDXXA(2,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXA(2,2,IAtm,2,JAtm)=GDXXA(2,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXA(2,1,IAtm,3,JAtm)=GDXXA(2,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXA(2,3,IAtm,1,JAtm)=GDXXA(2,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXA(2,2,IAtm,3,JAtm)=GDXXA(2,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXA(2,3,IAtm,2,JAtm)=GDXXA(2,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXA(2,3,IAtm,3,JAtm)=GDXXA(2,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXA(3,1,IAtm,1,JAtm)=GDXXA(3,1,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXA(3,1,IAtm,2,JAtm)=GDXXA(3,1,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXA(3,2,IAtm,1,JAtm)=GDXXA(3,2,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXA(3,2,IAtm,2,JAtm)=GDXXA(3,2,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXA(3,1,IAtm,3,JAtm)=GDXXA(3,1,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXA(3,3,IAtm,1,JAtm)=GDXXA(3,3,IAtm,1,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXA(3,2,IAtm,3,JAtm)=GDXXA(3,2,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXA(3,3,IAtm,2,JAtm)=GDXXA(3,3,IAtm,2,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXA(3,3,IAtm,3,JAtm)=GDXXA(3,3,IAtm,3,JAtm)+
     $                    DAIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
c     beta
              GDXXB(1,1,IAtm,1,JAtm)=GDXXB(1,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXX*VALJX+VAlX*VAOXX(1,J))
              GDXXB(1,1,IAtm,2,JAtm)=GDXXB(1,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXX*VALJY+VAlX*VAOXX(2,J))
              GDXXB(1,2,IAtm,1,JAtm)=GDXXB(1,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlY*VAOXX(1,J))
              GDXXB(1,2,IAtm,2,JAtm)=GDXXB(1,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlY*VAOXX(2,J))
              GDXXB(1,1,IAtm,3,JAtm)=GDXXB(1,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXX*VALJZ+VAlX*VAOXX(4,J))
              GDXXB(1,3,IAtm,1,JAtm)=GDXXB(1,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlZ*VAOXX(1,J))
              GDXXB(1,2,IAtm,3,JAtm)=GDXXB(1,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlY*VAOXX(4,J))
              GDXXB(1,3,IAtm,2,JAtm)=GDXXB(1,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlZ*VAOXX(2,J))
              GDXXB(1,3,IAtm,3,JAtm)=GDXXB(1,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlZ*VAOXX(4,J))
c
              GDXXB(2,1,IAtm,1,JAtm)=GDXXB(2,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXY*VALJX+VAlX*VAOXX(2,J))
              GDXXB(2,1,IAtm,2,JAtm)=GDXXB(2,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXY*VALJY+VAlX*VAOXX(3,J))
              GDXXB(2,2,IAtm,1,JAtm)=GDXXB(2,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYY*VALJX+VAlY*VAOXX(2,J))
              GDXXB(2,2,IAtm,2,JAtm)=GDXXB(2,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYY*VALJY+VAlY*VAOXX(3,J))
              GDXXB(2,1,IAtm,3,JAtm)=GDXXB(2,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXY*VALJZ+VAlX*VAOXX(5,J))
              GDXXB(2,3,IAtm,1,JAtm)=GDXXB(2,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlZ*VAOXX(2,J))
              GDXXB(2,2,IAtm,3,JAtm)=GDXXB(2,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYY*VALJZ+VAlY*VAOXX(5,J))
              GDXXB(2,3,IAtm,2,JAtm)=GDXXB(2,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlZ*VAOXX(3,J))
              GDXXB(2,3,IAtm,3,JAtm)=GDXXB(2,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlZ*VAOXX(5,J))
c
              GDXXB(3,1,IAtm,1,JAtm)=GDXXB(3,1,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJX+VAlX*VAOXX(4,J))
              GDXXB(3,1,IAtm,2,JAtm)=GDXXB(3,1,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJY+VAlX*VAOXX(5,J))
              GDXXB(3,2,IAtm,1,JAtm)=GDXXB(3,2,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJX+VAlY*VAOXX(4,J))
              GDXXB(3,2,IAtm,2,JAtm)=GDXXB(3,2,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJY+VAlY*VAOXX(5,J))
              GDXXB(3,1,IAtm,3,JAtm)=GDXXB(3,1,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALXZ*VALJZ+VAlX*VAOXX(6,J))
              GDXXB(3,3,IAtm,1,JAtm)=GDXXB(3,3,IAtm,1,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJX+VAlZ*VAOXX(4,J))
              GDXXB(3,2,IAtm,3,JAtm)=GDXXB(3,2,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALYZ*VALJZ+VAlY*VAOXX(6,J))
              GDXXB(3,3,IAtm,2,JAtm)=GDXXB(3,3,IAtm,2,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJY+VAlZ*VAOXX(5,J))
              GDXXB(3,3,IAtm,3,JAtm)=GDXXB(3,3,IAtm,3,JAtm)+
     $                    DBIJ2 * (VALZZ*VALJZ+VAlZ*VAOXX(6,J))
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient and Hessian of density gradient invariant
c
c
        DO IAtm=1,NAtoms
c   alpha alpha
          GXAA(1,IAtm)=two*(DAX*GDXA(1,1,IAtm)+
     $                   DAY*GDXA(2,1,IAtm)+DAZ*GDXA(3,1,IAtm))
          GXAA(2,IAtm)=two*(DAX*GDXA(1,2,IAtm)+
     $                   DAY*GDXA(2,2,IAtm)+DAZ*GDXA(3,2,IAtm))
          GXAA(3,IAtm)=two*(DAX*GDXA(1,3,IAtm)+
     $                   DAY*GDXA(2,3,IAtm)+DAZ*GDXA(3,3,IAtm))
c   beta beta
          GXBB(1,IAtm)=two*(DBX*GDXB(1,1,IAtm)+
     $                   DBY*GDXB(2,1,IAtm)+DBZ*GDXB(3,1,IAtm))
          GXBB(2,IAtm)=two*(DBX*GDXB(1,2,IAtm)+
     $                   DBY*GDXB(2,2,IAtm)+DBZ*GDXB(3,2,IAtm))
          GXBB(3,IAtm)=two*(DBX*GDXB(1,3,IAtm)+
     $                   DBY*GDXB(2,3,IAtm)+DBZ*GDXB(3,3,IAtm))
c   alpha beta
          GXAB(1,IAtm)=DAX*GDXB(1,1,IAtm)+DBX*GDXA(1,1,IAtm)+
     $                 DAY*GDXB(2,1,IAtm)+DBY*GDXA(2,1,IAtm)+
     $                 DAZ*GDXB(3,1,IAtm)+DBZ*GDXA(3,1,IAtm)
          GXAB(2,IAtm)=DAX*GDXB(1,2,IAtm)+DBX*GDXA(1,2,IAtm)+
     $                 DAY*GDXB(2,2,IAtm)+DBY*GDXA(2,2,IAtm)+
     $                 DAZ*GDXB(3,2,IAtm)+DBZ*GDXA(3,2,IAtm)
          GXAB(3,IAtm)=DAX*GDXB(1,3,IAtm)+DBX*GDXA(1,3,IAtm)+
     $                 DAY*GDXB(2,3,IAtm)+DBY*GDXA(2,3,IAtm)+
     $                 DAZ*GDXB(3,3,IAtm)+DBZ*GDXA(3,3,IAtm)
        EndDO
c
c  compute second derivatives of gradient invariants
c  only of needed
c
        if(dogradxx)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
c     alpha alpha
            GXXAA(1,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXA(3,1,JAtm))
            GXXAA(1,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXA(3,2,JAtm))
            GXXAA(2,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXA(3,2,JAtm))
            GXXAA(1,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,1,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXA(3,1,JAtm))
            GXXAA(2,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXA(3,3,JAtm))
            GXXAA(3,IAtm,2,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXA(3,2,JAtm))
            GXXAA(3,IAtm,3,JAtm)=Two*(
     $        DAX*GDXXA(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXA(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXA(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXA(3,3,JAtm))
c     beta beta
            GXXBB(1,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXB(3,1,JAtm))
            GXXBB(1,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXB(3,2,JAtm))
            GXXBB(2,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXB(3,2,JAtm))
            GXXBB(1,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,1,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXB(3,1,JAtm))
            GXXBB(2,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXB(3,3,JAtm))
            GXXBB(3,IAtm,2,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXB(3,2,JAtm))
            GXXBB(3,IAtm,3,JAtm)=Two*(
     $        DBX*GDXXB(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBY*GDXXB(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBZ*GDXXB(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXB(3,3,JAtm))
          EndDO
        EndDO
c
c     alpha beta is done separately, as only few functionals
c     actually need this one
c
        if(abs(gc).gt.epsi)then
        DO IAtm=1,NAtoms
          DO JAtm=IAtm,NAtoms
            GXXAB(1,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,1,JAtm)+GDXA(1,1,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,1,JAtm)+GDXB(1,1,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,1,JAtm)+GDXA(2,1,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,1,JAtm)+GDXB(2,1,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,1,JAtm)+GDXA(3,1,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,1,JAtm)+GDXB(3,1,IAtm)*GDXA(3,1,JAtm)
            GXXAB(2,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,1,JAtm)+GDXA(1,2,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,1,JAtm)+GDXB(1,2,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,1,JAtm)+GDXA(2,2,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,1,JAtm)+GDXB(2,2,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,1,JAtm)+GDXA(3,2,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,1,JAtm)+GDXB(3,2,IAtm)*GDXA(3,1,JAtm)
            GXXAB(3,IAtm,1,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,1,JAtm)+GDXA(1,3,IAtm)*GDXB(1,1,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,1,JAtm)+GDXB(1,3,IAtm)*GDXA(1,1,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,1,JAtm)+GDXA(2,3,IAtm)*GDXB(2,1,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,1,JAtm)+GDXB(2,3,IAtm)*GDXA(2,1,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,1,JAtm)+GDXA(3,3,IAtm)*GDXB(3,1,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,1,JAtm)+GDXB(3,3,IAtm)*GDXA(3,1,JAtm)
            GXXAB(1,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,2,JAtm)+GDXA(1,1,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,2,JAtm)+GDXB(1,1,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,2,JAtm)+GDXA(2,1,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,2,JAtm)+GDXB(2,1,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,2,JAtm)+GDXA(3,1,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,2,JAtm)+GDXB(3,1,IAtm)*GDXA(3,2,JAtm)
            GXXAB(2,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,2,JAtm)+GDXA(1,2,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,2,JAtm)+GDXB(1,2,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,2,JAtm)+GDXA(2,2,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,2,JAtm)+GDXB(2,2,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,2,JAtm)+GDXA(3,2,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,2,JAtm)+GDXB(3,2,IAtm)*GDXA(3,2,JAtm)
            GXXAB(3,IAtm,2,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,2,JAtm)+GDXA(1,3,IAtm)*GDXB(1,2,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,2,JAtm)+GDXB(1,3,IAtm)*GDXA(1,2,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,2,JAtm)+GDXA(2,3,IAtm)*GDXB(2,2,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,2,JAtm)+GDXB(2,3,IAtm)*GDXA(2,2,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,2,JAtm)+GDXA(3,3,IAtm)*GDXB(3,2,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,2,JAtm)+GDXB(3,3,IAtm)*GDXA(3,2,JAtm)
            GXXAB(1,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,1,IAtm,3,JAtm)+GDXA(1,1,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,1,IAtm,3,JAtm)+GDXB(1,1,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,1,IAtm,3,JAtm)+GDXA(2,1,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,1,IAtm,3,JAtm)+GDXB(2,1,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,1,IAtm,3,JAtm)+GDXA(3,1,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,1,IAtm,3,JAtm)+GDXB(3,1,IAtm)*GDXA(3,3,JAtm)
            GXXAB(2,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,2,IAtm,3,JAtm)+GDXA(1,2,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,2,IAtm,3,JAtm)+GDXB(1,2,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,2,IAtm,3,JAtm)+GDXA(2,2,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,2,IAtm,3,JAtm)+GDXB(2,2,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,2,IAtm,3,JAtm)+GDXA(3,2,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,2,IAtm,3,JAtm)+GDXB(3,2,IAtm)*GDXA(3,3,JAtm)
            GXXAB(3,IAtm,3,JAtm)=
     $        DAX*GDXXB(1,3,IAtm,3,JAtm)+GDXA(1,3,IAtm)*GDXB(1,3,JAtm)+
     $        DBX*GDXXA(1,3,IAtm,3,JAtm)+GDXB(1,3,IAtm)*GDXA(1,3,JAtm)+
     $        DAY*GDXXB(2,3,IAtm,3,JAtm)+GDXA(2,3,IAtm)*GDXB(2,3,JAtm)+
     $        DBY*GDXXA(2,3,IAtm,3,JAtm)+GDXB(2,3,IAtm)*GDXA(2,3,JAtm)+
     $        DAZ*GDXXB(3,3,IAtm,3,JAtm)+GDXA(3,3,IAtm)*GDXB(3,3,JAtm)+
     $        DBZ*GDXXA(3,3,IAtm,3,JAtm)+GDXB(3,3,IAtm)*GDXA(3,3,JAtm)
          EndDO
        EndDO
        else
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        else
          CALL ZeroIT(GXXAA,NAt3**2)
          CALL ZeroIT(GXXBB,NAt3**2)
          CALL ZeroIT(GXXAB,NAt3**2)
        endif
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
C    and also compute partial sums used in Hessian quadrature
c
        ga2= ga+ga
        gb2= gb+gb
        DAX2=DAX+DAX
        DAY2=DAY+DAY
        DAZ2=DAZ+DAZ
        DBX2=DBX+DBX
        DBY2=DBY+DBY
        DBZ2=DBZ+DBZ
c
c     coeff. X
c
        SXXA=ga2*DAX+gc*DBX
        SXYA=ga2*DAY+gc*DBY
        SXZA=ga2*DAZ+gc*DBZ
        SXXB=gb2*DBX+gc*DAX
        SXYB=gb2*DBY+gc*DAY
        SXZB=gb2*DBZ+gc*DAZ
c
        DO IAtm=1,NAtoms
c
c     coeff. V
c
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)+
     $        raga*gxaa(1,iatm)+ragb*gxbb(1,iatm)+ragc*gxab(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)+
     $        raga*gxaa(2,iatm)+ragb*gxbb(2,iatm)+ragc*gxab(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)+
     $        raga*gxaa(3,iatm)+ragb*gxbb(3,iatm)+ragc*gxab(3,iatm)
          svb(1,iatm)=rarb*denxa(1,iatm)+rbrb*denxb(1,iatm)+
     $        rbga*gxaa(1,iatm)+rbgb*gxbb(1,iatm)+rbgc*gxab(1,iatm)
          svb(2,iatm)=rarb*denxa(2,iatm)+rbrb*denxb(2,iatm)+
     $        rbga*gxaa(2,iatm)+rbgb*gxbb(2,iatm)+rbgc*gxab(2,iatm)
          svb(3,iatm)=rarb*denxa(3,iatm)+rbrb*denxb(3,iatm)+
     $        rbga*gxaa(3,iatm)+rbgb*gxbb(3,iatm)+rbgc*gxab(3,iatm)
c
c  partial sums for Hessian
c
          sga(1,iatm)=raga*denxa(1,iatm)+rbga*denxb(1,iatm)+
     $        gaga*gxaa(1,iatm)+gagb*gxbb(1,iatm)+gagc*gxab(1,iatm)
          sga(2,iatm)=raga*denxa(2,iatm)+rbga*denxb(2,iatm)+
     $        gaga*gxaa(2,iatm)+gagb*gxbb(2,iatm)+gagc*gxab(2,iatm)
          sga(3,iatm)=raga*denxa(3,iatm)+rbga*denxb(3,iatm)+
     $        gaga*gxaa(3,iatm)+gagb*gxbb(3,iatm)+gagc*gxab(3,iatm)
c
          sgb(1,iatm)=ragb*denxa(1,iatm)+rbgb*denxb(1,iatm)+
     $        gagb*gxaa(1,iatm)+gbgb*gxbb(1,iatm)+gbgc*gxab(1,iatm)
          sgb(2,iatm)=ragb*denxa(2,iatm)+rbgb*denxb(2,iatm)+
     $        gagb*gxaa(2,iatm)+gbgb*gxbb(2,iatm)+gbgc*gxab(2,iatm)
          sgb(3,iatm)=ragb*denxa(3,iatm)+rbgb*denxb(3,iatm)+
     $        gagb*gxaa(3,iatm)+gbgb*gxbb(3,iatm)+gbgc*gxab(3,iatm)
c
          sgc(1,iatm)=ragc*denxa(1,iatm)+rbgc*denxb(1,iatm)+
     $        gagc*gxaa(1,iatm)+gbgc*gxbb(1,iatm)+gcgc*gxab(1,iatm)
          sgc(2,iatm)=ragc*denxa(2,iatm)+rbgc*denxb(2,iatm)+
     $        gagc*gxaa(2,iatm)+gbgc*gxbb(2,iatm)+gcgc*gxab(2,iatm)
          sgc(3,iatm)=ragc*denxa(3,iatm)+rbgc*denxb(3,iatm)+
     $        gagc*gxaa(3,iatm)+gbgc*gxbb(3,iatm)+gcgc*gxab(3,iatm)
c
c     coeff. W
c
          swa(1,1,iatm)=ga2*gdxa(1,1,iatm)+gc*gdxb(1,1,Iatm)+
     $        DAX2*sga(1,iatm)+DBX*sgc(1,iatm)
          swa(1,2,iatm)=ga2*gdxa(1,2,iatm)+gc*gdxb(1,2,Iatm)+
     $        DAX2*sga(2,iatm)+DBX*sgc(2,iatm)
          swa(1,3,iatm)=ga2*gdxa(1,3,iatm)+gc*gdxb(1,3,Iatm)+
     $        DAX2*sga(3,iatm)+DBX*sgc(3,iatm)
          swa(2,1,iatm)=ga2*gdxa(2,1,iatm)+gc*gdxb(2,1,Iatm)+
     $        DAY2*sga(1,iatm)+DBY*sgc(1,iatm)
          swa(2,2,iatm)=ga2*gdxa(2,2,iatm)+gc*gdxb(2,2,Iatm)+
     $        DAY2*sga(2,iatm)+DBY*sgc(2,iatm)
          swa(2,3,iatm)=ga2*gdxa(2,3,iatm)+gc*gdxb(2,3,Iatm)+
     $        DAY2*sga(3,iatm)+DBY*sgc(3,iatm)
          swa(3,1,iatm)=ga2*gdxa(3,1,iatm)+gc*gdxb(3,1,Iatm)+
     $        DAZ2*sga(1,iatm)+DBZ*sgc(1,iatm)
          swa(3,2,iatm)=ga2*gdxa(3,2,iatm)+gc*gdxb(3,2,Iatm)+
     $        DAZ2*sga(2,iatm)+DBZ*sgc(2,iatm)
          swa(3,3,iatm)=ga2*gdxa(3,3,iatm)+gc*gdxb(3,3,Iatm)+
     $        DAZ2*sga(3,iatm)+DBZ*sgc(3,iatm)
c
          swb(1,1,iatm)=gb2*gdxb(1,1,iatm)+gc*gdxa(1,1,Iatm)+
     $        DBX2*sgb(1,iatm)+DAX*sgc(1,iatm)
          swb(1,2,iatm)=gb2*gdxb(1,2,iatm)+gc*gdxa(1,2,Iatm)+
     $        DBX2*sgb(2,iatm)+DAX*sgc(2,iatm)
          swb(1,3,iatm)=gb2*gdxb(1,3,iatm)+gc*gdxa(1,3,Iatm)+
     $        DBX2*sgb(3,iatm)+DAX*sgc(3,iatm)
          swb(2,1,iatm)=gb2*gdxb(2,1,iatm)+gc*gdxa(2,1,Iatm)+
     $        DBY2*sgb(1,iatm)+DAY*sgc(1,iatm)
          swb(2,2,iatm)=gb2*gdxb(2,2,iatm)+gc*gdxa(2,2,Iatm)+
     $        DBY2*sgb(2,iatm)+DAY*sgc(2,iatm)
          swb(2,3,iatm)=gb2*gdxb(2,3,iatm)+gc*gdxa(2,3,Iatm)+
     $        DBY2*sgb(3,iatm)+DAY*sgc(3,iatm)
          swb(3,1,iatm)=gb2*gdxb(3,1,iatm)+gc*gdxa(3,1,Iatm)+
     $        DBZ2*sgb(1,iatm)+DAZ*sgc(1,iatm)
          swb(3,2,iatm)=gb2*gdxb(3,2,iatm)+gc*gdxa(3,2,Iatm)+
     $        DBZ2*sgb(2,iatm)+DAZ*sgc(2,iatm)
          swb(3,3,iatm)=gb2*gdxb(3,3,iatm)+gc*gdxa(3,3,Iatm)+
     $        DBZ2*sgb(3,iatm)+DAZ*sgc(3,iatm)
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    Fock Matrix.
c
c     get the maximum absolute value of coefficients V and W
c
        Call absmax(NAt3,sva,iixx,svamax)
        Call absmax(NAt3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        Call absmax(3*NAt3,swa,iixx,swamax)
        Call absmax(3*NAt3,swb,iixx,swbmax)
        swmax=max(swamax,swbmax)
        tswmax=Three*swmax
c
c -- global threshold testing
c
        abr=max(abs(ra),abs(rb))
        vmx2=vmx*vmx
        SMax = max(Abs(SXXA+SXYA+SXZA),abs(SXXB+SXYB+SXZB))
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abr+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVIA=SXXA*ValIX+SXYA*ValIY+SXZA*ValIZ+ra*ValI
          XVIB=SXXB*ValIX+SXYB*ValIY+SXZB*ValIZ+rb*ValI
          SXXIA=SXXA*ValI
          SXYIA=SXYA*ValI
          SXZIA=SXZA*ValI
          SXXIB=SXXB*ValI
          SXYIB=SXYB*ValI
          SXZIB=SXZB*ValI
          XXVXA=SXXA*VAOXX(1,I)+SXYA*VAOXX(2,I)+SXZA*VAOXX(4,I)+ra*VALIX
          XXVYA=SXXA*VAOXX(2,I)+SXYA*VAOXX(3,I)+SXZA*VAOXX(5,I)+ra*VALIY
          XXVZA=SXXA*VAOXX(4,I)+SXYA*VAOXX(5,I)+SXZA*VAOXX(6,I)+ra*VALIZ
          XXVXB=SXXB*VAOXX(1,I)+SXYB*VAOXX(2,I)+SXZB*VAOXX(4,I)+rb*VALIX
          XXVYB=SXXB*VAOXX(2,I)+SXYB*VAOXX(3,I)+SXZB*VAOXX(5,I)+rb*VALIY
          XXVZB=SXXB*VAOXX(4,I)+SXYB*VAOXX(5,I)+SXZB*VAOXX(6,I)+rb*VALIZ
          valm=abs(vali)
          valm1=max(abs(xvia),abs(xvib))
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1a=max(abs(xxvxa),abs(xxvya),abs(xxvza))
          valmx1b=max(abs(xxvxb),abs(xxvyb),abs(xxvzb))
          valmx1=max(valmx1a,valmx1b)
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJA=SXXA*ValJX+SXYA*ValJY+SXZA*ValJZ
              XVJB=SXXB*ValJX+SXYB*ValJY+SXZB*ValJZ
              abxvj=max(abs(xvja),abs(xvjb))
              if(valmx1*abs(valj)+valmx*abxvj.gt.thrsh)then
                if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                  FDA(1,IAtm,IJ)=FDA(1,IAtm,IJ)-XXVXA*ValJ-ValIX*XVJA
                  FDA(2,IAtm,IJ)=FDA(2,IAtm,IJ)-XXVYA*ValJ-ValIY*XVJA
                  FDA(3,IAtm,IJ)=FDA(3,IAtm,IJ)-XXVZA*ValJ-ValIZ*XVJA
                  FDB(1,IAtm,IJ)=FDB(1,IAtm,IJ)-XXVXB*ValJ-ValIX*XVJB
                  FDB(2,IAtm,IJ)=FDB(2,IAtm,IJ)-XXVYB*ValJ-ValIY*XVJB
                  FDB(3,IAtm,IJ)=FDB(3,IAtm,IJ)-XXVZB*ValJ-ValIZ*XVJB
                endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
                if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                  FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - XVIA*ValJX -
     $            VAOXX(1,J)*SXXIA - VAOXX(2,J)*SXYIA - VAOXX(4,J)*SXZIA
                  FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - XVIA*ValJY -
     $            VAOXX(2,J)*SXXIA - VAOXX(3,J)*SXYIA - VAOXX(5,J)*SXZIA
                  FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - XVIA*ValJZ -
     $            VAOXX(4,J)*SXXIA - VAOXX(5,J)*SXYIA - VAOXX(6,J)*SXZIA
c
                  FDB(1,JAtm,IJ) = FDB(1,JAtm,IJ) - XVIB*ValJX -
     $            VAOXX(1,J)*SXXIB - VAOXX(2,J)*SXYIB - VAOXX(4,J)*SXZIB
                  FDB(2,JAtm,IJ) = FDB(2,JAtm,IJ) - XVIB*ValJY -
     $            VAOXX(2,J)*SXXIB - VAOXX(3,J)*SXYIB - VAOXX(5,J)*SXZIB
                  FDB(3,JAtm,IJ) = FDB(3,JAtm,IJ) - XVIB*ValJZ -
     $            VAOXX(4,J)*SXXIB - VAOXX(5,J)*SXYIB - VAOXX(6,J)*SXZIB
                endif
              endif
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VALJX*VALI
              VIJY=VALIY*VALJ+VALJY*VALI
              VIJZ=VALIZ*VALJ+VALJZ*VALI
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                Do KAtm=Nb,Ne
                  FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ)+(SVA(1,KAtm)*VIJ+
     $                           SWA(1,1,KAtm)*VIJX+SWA(2,1,KAtm)*VIJY+
     $                           SWA(3,1,KAtM)*VIJZ)
                  FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ)+(SVA(2,KAtm)*VIJ+
     $                           SWA(1,2,KAtm)*VIJX+SWA(2,2,KAtm)*VIJY+
     $                           SWA(3,2,KAtM)*VIJZ)
                  FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ)+(SVA(3,KAtm)*VIJ+
     $                           SWA(1,3,KAtm)*VIJX+SWA(2,3,KAtm)*VIJY+
     $                           SWA(3,3,KatM)*VIJZ)
c
                  FDB(1,KAtm,IJ)=FDB(1,KAtm,IJ)+(SVB(1,KAtm)*VIJ+
     $                           SWB(1,1,KAtm)*VIJX+SWB(2,1,KAtm)*VIJY+
     $                           SWB(3,1,KAtM)*VIJZ)
                  FDB(2,KAtm,IJ)=FDB(2,KAtm,IJ)+(SVB(2,KAtm)*VIJ+
     $                           SWB(1,2,KAtm)*VIJX+SWB(2,2,KAtm)*VIJY+
     $                           SWB(3,2,KAtM)*VIJZ)
                  FDB(3,KAtm,IJ)=FDB(3,KAtm,IJ)+(SVB(3,KAtm)*VIJ+
     $                           SWB(1,3,KAtm)*VIJX+SWB(2,3,KAtm)*VIJY+
     $                           SWB(3,3,KatM)*VIJZ)
                EndDo
              EndIf
 230        CONTINUE
          ENDIF
 240    CONTINUE
c
c   Numerical quadrature for direct contribution to Hessian matrix
c
 245    CONTINUE
        call secund(t5)
        tqf=tqf+t5-t4
        Do IAtm=1,NAtoms
          Do JAtm=IAtm,NAtoms
            HESS(1,IAtm,1,JAtm) = HESS(1,IAtm,1,JAtm) +
     $       sva(1,IAtm)*DenXA(1,JAtm) + svb(1,IAtm)*DenXB(1,JAtm) +
     $       sga(1,IAtm)*gxaa(1,JAtm)  + sgb(1,IAtm)*gxbb(1,JAtm) +
     $       sgc(1,IAtm)*gxab(1,JAtm)  + ra*DenXXA(1,IAtm,1,JAtm) +
     $       rb*DenXXB(1,IAtm,1,JAtm)  + ga*gxxaa(1,IAtm,1,JAtm) +
     $       gb*gxxbb(1,IAtm,1,JAtm)   + gc*gxxab(1,IAtm,1,JAtm)
            HESS(2,IAtm,1,JAtm) = HESS(2,IAtm,1,JAtm) +
     $       sva(2,IAtm)*DenXA(1,JAtm) + svb(2,IAtm)*DenXB(1,JAtm) +
     $       sga(2,IAtm)*gxaa(1,JAtm)  + sgb(2,IAtm)*gxbb(1,JAtm) +
     $       sgc(2,IAtm)*gxab(1,JAtm)  + ra*DenXXA(2,IAtm,1,JAtm) +
     $       rb*DenXXB(2,IAtm,1,JAtm)  + ga*gxxaa(2,IAtm,1,JAtm) +
     $       gb*gxxbb(2,IAtm,1,JAtm)   + gc*gxxab(2,IAtm,1,JAtm)
            HESS(3,IAtm,1,JAtm) = HESS(3,IAtm,1,JAtm) +
     $       sva(3,IAtm)*DenXA(1,JAtm) + svb(3,IAtm)*DenXB(1,JAtm) +
     $       sga(3,IAtm)*gxaa(1,JAtm)  + sgb(3,IAtm)*gxbb(1,JAtm) +
     $       sgc(3,IAtm)*gxab(1,JAtm)  + ra*DenXXA(3,IAtm,1,JAtm) +
     $       rb*DenXXB(3,IAtm,1,JAtm)  + ga*gxxaa(3,IAtm,1,JAtm) +
     $       gb*gxxbb(3,IAtm,1,JAtm)   + gc*gxxab(3,IAtm,1,JAtm)
            HESS(1,IAtm,2,JAtm) = HESS(1,IAtm,2,JAtm) +
     $       sva(1,IAtm)*DenXA(2,JAtm) + svb(1,IAtm)*DenXB(2,JAtm) +
     $       sga(1,IAtm)*gxaa(2,JAtm)  + sgb(1,IAtm)*gxbb(2,JAtm) +
     $       sgc(1,IAtm)*gxab(2,JAtm)  + ra*DenXXA(1,IAtm,2,JAtm) +
     $       rb*DenXXB(1,IAtm,2,JAtm)  + ga*gxxaa(1,IAtm,2,JAtm) +
     $       gb*gxxbb(1,IAtm,2,JAtm)   + gc*gxxab(1,IAtm,2,JAtm)
            HESS(2,IAtm,2,JAtm) = HESS(2,IAtm,2,JAtm) +
     $       sva(2,IAtm)*DenXA(2,JAtm) + svb(2,IAtm)*DenXB(2,JAtm) +
     $       sga(2,IAtm)*gxaa(2,JAtm)  + sgb(2,IAtm)*gxbb(2,JAtm) +
     $       sgc(2,IAtm)*gxab(2,JAtm)  + ra*DenXXA(2,IAtm,2,JAtm) +
     $       rb*DenXXB(2,IAtm,2,JAtm)  + ga*gxxaa(2,IAtm,2,JAtm) +
     $       gb*gxxbb(2,IAtm,2,JAtm)   + gc*gxxab(2,IAtm,2,JAtm)
            HESS(3,IAtm,2,JAtm) = HESS(3,IAtm,2,JAtm) +
     $       sva(3,IAtm)*DenXA(2,JAtm) + svb(3,IAtm)*DenXB(2,JAtm) +
     $       sga(3,IAtm)*gxaa(2,JAtm)  + sgb(3,IAtm)*gxbb(2,JAtm) +
     $       sgc(3,IAtm)*gxab(2,JAtm)  + ra*DenXXA(3,IAtm,2,JAtm) +
     $       rb*DenXXB(3,IAtm,2,JAtm)  + ga*gxxaa(3,IAtm,2,JAtm) +
     $       gb*gxxbb(3,IAtm,2,JAtm)   + gc*gxxab(3,IAtm,2,JAtm)
            HESS(1,IAtm,3,JAtm) = HESS(1,IAtm,3,JAtm) +
     $       sva(1,IAtm)*DenXA(3,JAtm) + svb(1,IAtm)*DenXB(3,JAtm) +
     $       sga(1,IAtm)*gxaa(3,JAtm)  + sgb(1,IAtm)*gxbb(3,JAtm) +
     $       sgc(1,IAtm)*gxab(3,JAtm)  + ra*DenXXA(1,IAtm,3,JAtm) +
     $       rb*DenXXB(1,IAtm,3,JAtm)  + ga*gxxaa(1,IAtm,3,JAtm) +
     $       gb*gxxbb(1,IAtm,3,JAtm)   + gc*gxxab(1,IAtm,3,JAtm)
            HESS(2,IAtm,3,JAtm) = HESS(2,IAtm,3,JAtm) +
     $       sva(2,IAtm)*DenXA(3,JAtm) + svb(2,IAtm)*DenXB(3,JAtm) +
     $       sga(2,IAtm)*gxaa(3,JAtm)  + sgb(2,IAtm)*gxbb(3,JAtm) +
     $       sgc(2,IAtm)*gxab(3,JAtm)  + ra*DenXXA(2,IAtm,3,JAtm) +
     $       rb*DenXXB(2,IAtm,3,JAtm)  + ga*gxxaa(2,IAtm,3,JAtm) +
     $       gb*gxxbb(2,IAtm,3,JAtm)   + gc*gxxab(2,IAtm,3,JAtm)
            HESS(3,IAtm,3,JAtm) = HESS(3,IAtm,3,JAtm) +
     $       sva(3,IAtm)*DenXA(3,JAtm) + svb(3,IAtm)*DenXB(3,JAtm) +
     $       sga(3,IAtm)*gxaa(3,JAtm)  + sgb(3,IAtm)*gxbb(3,JAtm) +
     $       sgc(3,IAtm)*gxab(3,JAtm)  + ra*DenXXA(3,IAtm,3,JAtm) +
     $       rb*DenXXB(3,IAtm,3,JAtm)  + ga*gxxaa(3,IAtm,3,JAtm) +
     $       gb*gxxbb(3,IAtm,3,JAtm)   + gc*gxxab(3,IAtm,3,JAtm)
          EndDo
        EndDo
        call secund(t6)
        tqh=tqh+t6-t5
c
 200  CONTINUE
      ENDIF
C
      RETURN
      END
c======================================================================
      subroutine mkmpf(dft,    npp,    nbas,   nbf,    natoms,
     $                 nb,     ne,     nbatm,  thrsh,  da,
     $                 dm,     gden,   wght,   pra,    prara,
     $                 prarb,  pga,    pgc,    praga,  pragb,
     $                 pragc,  pgaga,  pgagb,  pgagc,  pgcgc,
     $                 vao,    vaox,   vaoxx,  inb,    vm,
     $                 indx,   denx,   gdx,    gx,     sv,
     $                 sw,     icntr,  fda,    td1,    tg1,
     $                 tsw,    tqf)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into derivative Fock matrices
C  ** CLOSED SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C            1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   -  number of contributing (non-zero) grid points this batch
C  NBas  -  total number of basis functions
C  nbf   -  indexing array to location of "non-zero" basis functions
C  NAtoms-  number of atoms
C  Nb,Ne -  First and last component of Fock derivatives to compute
C  NBAtm -  basis functions --> atom index
C  thrsh -  threshold for neglect of contribution
C  DA    -  closed-shell density matrix (lower triangle)
C  DM    -  maximum density matrix element per column
C  GDEN  -  density gradient at grid points (non local only, dft > 3)
C  WGHT  -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density (dft > 3)
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
C  VAO   -  "non-zero" basis function values at grid points
C  VAOX  -  basis function 1st derivatives at grid points
C  VAOXX -  basis function 2nd derivatives at grid points
C  INB   -  indexing array for non-zero entries to VAO
C  VM    -  array containing maximum magnitude AO per grid point
C  INDX  -  index into contributing columns of VAO
C  DenX  -  scratch storage for in situ atomic gradient of the
C           density
C  GDX   -  ditto for atomic gradient of density gradient (dft > 3)
C  GX    -  ditto for atomic gradient of gradient invariant (dft > 3)
C  SV    -  storage for coefficient of vao(i)*vao(j) in quadrature
C        -  formula (dft > 3)
C  SW    -  storage for coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  ICntr -  current atomic center
C
C  on exit
C
C  FDA     -  contribution to derivative Fock matrices
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),VM(*),
     $          DA(NBas*(NBas+1)/2),DM(NBas)
      DIMENSION DenX(3,NAtoms)
      DIMENSION GDEN(3,*)
      DIMENSION GDX(3,3,NAtoms),GX(3,NAtoms)
      DIMENSION SV(3,NAtoms),SW(3,3,NAtoms)
      INTEGER   dft,nbf(*),INB(*),INDX(NPP),NBAtm(NBas)
      DIMENSION FDA(3,Nb:Ne,*)
      dimension pra(npp),pga(npp),pgc(npp),prara(npp)
      dimension prarb(npp),praga(npp),pragb(npp),pragc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0)
      PARAMETER (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,dogxx,doval
C
C
      NAt3 = 3*NAtoms
c
      IF(dft.LE.3) THEN
C
C  local only
C
      DO 50 IP=1,NPP
      IPP = INDX(IP)
      WG=WGHT(IPP)
      ra = pra(IP)*WG
      rara = prara(IP)*WG
      VMx = VM(IPP)
      CALL ZeroIT(DenX,NAt3)
c
c -- form gradient density at current grid point
c
c    NOTE:
c    for the closed shell case, the array DA contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient densities is omitted here, so
c    denx will contain half the gradient of the total closed
c    shell density
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrara=abs(rara)
      dodenx=abrara.gt.epsi
      thrx=unbelpo
      if(dodenx)thrx=thrsh/abrara
      thrx1=thrx/vmx2
      abra=abs(ra)
      if(.not.dodenx)goto 25 ! should never happen!
c
c   now compute density derivatives
c
      DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        DMx = DM(II)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        if(IAtm.lt.Nb.or.Iatm.gt.Ne)goto 20
        ValX = VAOX(1,I)
        ValY = VAOX(2,I)
        ValZ = VAOX(3,I)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx*vmx
        dox=(dodenx.and.(valt.gt.thrx1.or.valt**2.gt.thrx))
        if(.not.dox) goto 20
        DO 18 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          DAIJ = DA(IJ)*VAO(J)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
c
c -- atomic gradient of density
          DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
          DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
          DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
 18     CONTINUE
        DO 19 J=I+1,nbf(IPP+1)
          JJ = INB(J)
          IJ = (JJ*(JJ-1))/2 + II
          DAIJ = DA(IJ)*VAO(J)
          abdaij=abs(daij)
          abdaijm=abdaij*vmx
c
c -- atomic gradient of density
          DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
          DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
          DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
 19     CONTINUE
 20   CONTINUE
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c    Numerical quadrature for the derivative Fock matrix.
c
c -- get maximum absolute value of 1st-order density at this grid point
      Call absmax(NAt3,DenX,iixx,DMaxyz)
c
c -- global threshold testing
c
      VMax = Max(Abra*VMx2,Abrara*DMaxyz*VMx2)
      If(VMax.LT.thrsh) GO TO 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        Val = VAO(I)*ra
        VDXC = VAO(I)*rara
        ValX = VAOX(1,I)*ra
        ValY = VAOX(2,I)*ra
        ValZ = VAOX(3,I)*ra
        abval= abs(val)
        doval= abval.gt.thrsh1
        valm = max(abs(valx),abs(valy),abs(valz))
        Valt = MAX(abval,valm,Abs(vdxc)*dmaxyz)
        IF(Valt.GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            ValJ = VAO(J)
            if(abs(valj)*valm.gt.thrsh)then
              if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - ValX*ValJ
                FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - ValY*ValJ
                FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - ValZ*ValJ
              endif
            endif
            if(doval)then
              if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - Val*VAOX(1,J)
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - Val*VAOX(2,J)
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - Val*VAOX(3,J)
              endif
            endif
c
c -- contribution of functional derivative
c
            VDJX = VDXC*ValJ
            If(Abs(VDJX)*DMaxyz.GT.thrsh) Then
              Do KAtm=Nb,Ne
                FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ) + VDJX*DENX(1,KAtm)
                FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ) + VDJX*DENX(2,KAtm)
                FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ) + VDJX*DENX(3,KAtm)
              EndDo
            EndIf
 30       CONTINUE
        ENDIF
 40   CONTINUE
 45   CONTINUE
      call secund(t3)
      tqf=tqf+t3-t2
 50   CONTINUE
cc
      ELSE
C
C   non-local dft
C
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        ra = pra(ip)*wg
        ga = pga(ip)*wg
        gc = pgc(ip)*wg
        rara = (prara(ip)+prarb(ip))*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  some sums of the potentials that will be used later
c
        prg=raga+ragb+ragc
        pgg=gaga+gagb+Two*gagc+Half*gcgc
        pg=ga+Half*gc
c
c  density gradient at current point
c
        DX=GDEN(1,IPP)
        DY=GDEN(2,IPP)
        DZ=GDEN(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenX,NAt3)
        CALL ZeroIT(GDX,3*NAt3)
c
c    form the atomic gradient of density
c    and density gradient at current grid point
c
c    NOTE:
c    for the closed shell case, the array DA contains the
c    total density matrix, thus the factor two in the definition
c    of the gradient densities is omitted here, so
c    denx will contain half the atomic gradient
c    of the total closed shell density. Likewise, gdx will
c    contain half the atomic gradient of the total closed
c    shell density gradient invariant
c
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          DMx = DM(II)       ! max. element of first order density
          if(dmx.gt.epsi)then
             thtest=thrsh/(vmx*dmx)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          if(IAtm.lt.Nb.or.Iatm.gt.Ne)goto 220
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJT=DA(IJ)
            abijt=abs(daijt)
            DAIJ = DAIJT*VALJ
            abij=abs(daij)
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
              DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
              DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              GDX(1,1,IAtm)=GDX(1,1,IAtm)-DAIJT*(VALJ*VALXX+VALX*VALJX)
              GDX(1,2,IAtm)=GDX(1,2,IAtm)-DAIJT*(VALJ*VALXY+VALY*VALJX)
              GDX(1,3,IAtm)=GDX(1,3,IAtm)-DAIJT*(VALJ*VALXZ+VALZ*VALJX)
              GDX(2,1,IAtm)=GDX(2,1,IAtm)-DAIJT*(VALJ*VALXY+VALX*VALJY)
              GDX(2,2,IAtm)=GDX(2,2,IAtm)-DAIJT*(VALJ*VALYY+VALY*VALJY)
              GDX(2,3,IAtm)=GDX(2,3,IAtm)-DAIJT*(VALJ*VALYZ+VALZ*VALJY)
              GDX(3,1,IAtm)=GDX(3,1,IAtm)-DAIJT*(VALJ*VALXZ+VALX*VALJZ)
              GDX(3,2,IAtm)=GDX(3,2,IAtm)-DAIJT*(VALJ*VALYZ+VALY*VALJZ)
              GDX(3,3,IAtm)=GDX(3,3,IAtm)-DAIJT*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJT= DA(IJ)
            abijt=abs(daijt)
            DAIJ = DAIJT*VALJ
            abij=abs(daij)
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenX(1,IAtm) = DenX(1,IAtm) - DAIJ*ValX
              DenX(2,IAtm) = DenX(2,IAtm) - DAIJ*ValY
              DenX(3,IAtm) = DenX(3,IAtm) - DAIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
              GDX(1,1,IAtm)=GDX(1,1,IAtm)-DAIJT*(VALJ*VALXX+VALX*VALJX)
              GDX(1,2,IAtm)=GDX(1,2,IAtm)-DAIJT*(VALJ*VALXY+VALY*VALJX)
              GDX(1,3,IAtm)=GDX(1,3,IAtm)-DAIJT*(VALJ*VALXZ+VALZ*VALJX)
              GDX(2,1,IAtm)=GDX(2,1,IAtm)-DAIJT*(VALJ*VALXY+VALX*VALJY)
              GDX(2,2,IAtm)=GDX(2,2,IAtm)-DAIJT*(VALJ*VALYY+VALY*VALJY)
              GDX(2,3,IAtm)=GDX(2,3,IAtm)-DAIJT*(VALJ*VALYZ+VALZ*VALJY)
              GDX(3,1,IAtm)=GDX(3,1,IAtm)-DAIJT*(VALJ*VALXZ+VALX*VALJZ)
              GDX(3,2,IAtm)=GDX(3,2,IAtm)-DAIJT*(VALJ*VALYZ+VALY*VALJZ)
              GDX(3,3,IAtm)=GDX(3,3,IAtm)-DAIJT*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient of density gradient invariant
c
c
        DO IAtm=1,NAtoms
          GX(1,IAtm)=two*(DX*GDX(1,1,IAtm)+
     $                   DY*GDX(2,1,IAtm)+DZ*GDX(3,1,IAtm))
          GX(2,IAtm)=two*(DX*GDX(1,2,IAtm)+
     $                   DY*GDX(2,2,IAtm)+DZ*GDX(3,2,IAtm))
          GX(3,IAtm)=two*(DX*GDX(1,3,IAtm)+
     $                   DY*GDX(2,3,IAtm)+DZ*GDX(3,3,IAtm))
        EndDO
c
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
c
        prg2=two*prg
        pgg2=two*pgg
        pg2=two*pg
        SXX=pg2*DX
        SXY=pg2*DY
        SXZ=pg2*DZ
        DO IAtm=1,NAtoms
          vd1=prg2*denx(1,iatm)+pgg2*gx(1,iatm)
          vd2=prg2*denx(2,iatm)+pgg2*gx(2,iatm)
          vd3=prg2*denx(3,iatm)+pgg2*gx(3,iatm)
          sv(1,iatm)=rara*denx(1,iatm)+prg*gx(1,iatm)
          sv(2,iatm)=rara*denx(2,iatm)+prg*gx(2,iatm)
          sv(3,iatm)=rara*denx(3,iatm)+prg*gx(3,iatm)
          SW(1,1,IAtm)=pg2*GDX(1,1,IAtm)+DX*VD1
          SW(2,1,IAtm)=pg2*GDX(2,1,IAtm)+DY*VD1
          SW(3,1,IAtm)=pg2*GDX(3,1,IAtm)+DZ*VD1
          SW(1,2,IAtm)=pg2*GDX(1,2,IAtm)+DX*VD2
          SW(2,2,IAtm)=pg2*GDX(2,2,IAtm)+DY*VD2
          SW(3,2,IAtm)=pg2*GDX(3,2,IAtm)+DZ*VD2
          SW(1,3,IAtm)=pg2*GDX(1,3,IAtm)+DX*VD3
          SW(2,3,IAtm)=pg2*GDX(2,3,IAtm)+DY*VD3
          SW(3,3,IAtm)=pg2*GDX(3,3,IAtm)+DZ*VD3
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    Fock Matrix.
c
c     get the maximum absolute value of coefficients V and W
c
        Call absmax(NAt3,SV,isv,svmax)
        Call absmax(3*NAt3,SW,isw,swmax)
        tswmax=Three*swmax
c -- global threshold testing
        abra=abs(ra)
        vmx2=vmx*vmx
        SMax = Abs(SXX+SXY+SXZ)
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abra+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVI=SXX*ValIX+SXY*ValIY+SXZ*ValIZ+ra*ValI
          SXXI=SXX*ValI
          SXYI=SXY*ValI
          SXZI=SXZ*ValI
          XXVX=SXX*VAOXX(1,I)+SXY*VAOXX(2,I)+SXZ*VAOXX(4,I)+ra*VALIX
          XXVY=SXX*VAOXX(2,I)+SXY*VAOXX(3,I)+SXZ*VAOXX(5,I)+ra*VALIY
          XXVZ=SXX*VAOXX(4,I)+SXY*VAOXX(5,I)+SXZ*VAOXX(6,I)+ra*VALIZ
          valm=abs(vali)
          valm1=abs(xvi)
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1=max(abs(xxvx),abs(xxvy),abs(xxvz))
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJ=SXX*ValJX+SXY*ValJY+SXZ*ValJZ
              if(valmx1*abs(valj)+valmx*abs(xvj).gt.thrsh)then
               if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                 FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - XXVX*ValJ - ValIX*XVJ
                 FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - XXVY*ValJ - ValIY*XVJ
                 FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - XXVZ*ValJ - ValIZ*XVJ
               endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
               if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - XVI*ValJX -
     $          VAOXX(1,J)*SXXI - VAOXX(2,J)*SXYI - VAOXX(4,J)*SXZI
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - XVI*ValJY -
     $          VAOXX(2,J)*SXXI - VAOXX(3,J)*SXYI - VAOXX(5,J)*SXZI
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - XVI*ValJZ -
     $          VAOXX(4,J)*SXXI - VAOXX(5,J)*SXYI - VAOXX(6,J)*SXZI
               endif
              endif
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VALJX*VALI
              VIJY=VALIY*VALJ+VALJY*VALI
              VIJZ=VALIZ*VALJ+VALJZ*VALI
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                Do KAtm=Nb,Ne
                  FDA(1,KAtm,IJ) = FDA(1,KAtm,IJ) + (SV(1,KAtm)*VIJ+
     $            SW(1,1,KAtm)*VIJX+SW(2,1,KAtm)*VIJY+SW(3,1,KAtM)*VIJZ)
                  FDA(2,KAtm,IJ) = FDA(2,KAtm,IJ) + (SV(2,KAtm)*VIJ+
     $            SW(1,2,KAtm)*VIJX+SW(2,2,KAtm)*VIJY+SW(3,2,KAtM)*VIJZ)
                  FDA(3,KAtm,IJ) = FDA(3,KAtm,IJ) + (SV(3,KAtm)*VIJ+
     $            SW(1,3,KAtm)*VIJX+SW(2,3,KAtm)*VIJY+SW(3,3,KatM)*VIJZ)
                EndDo
              EndIf
 230        CONTINUE
          ENDIF
 240    CONTINUE
 245    CONTINUE
        call secund(t5)
        tqf=tqf+t5-t4
c
 200  CONTINUE
      ENDIF
C
      RETURN
      END
c =====================================================================
      subroutine mkmpfu(dft,    npp,    nbas,   nbf,    natoms,
     $                   nb,     ne,     nbatm,  thrsh,  da,
     $                   db,     dm,     gdena,  gdenb,  wght,
     $                   pra,    prb,    prara,  prbrb,  prarb,
     $                   pga,    pgb,    pgc,    praga,  pragb,
     $                   pragc,  prbga,  prbgb,  prbgc,  pgaga,
     $                   pgagb,  pgagc,  pgbgb,  pgbgc,  pgcgc,
     $                   vao,    vaox,   vaoxx,  inb,    vm,
     $                   indx,   denxa,  denxb,  gdxa,   gdxb,
     $                   gxaa,   gxbb,   gxab,   sva,    svb,
     $                   swa,    swb,    sga,    sgb,    sgc,
     $                   icntr,  fda,    fdb,    td1,    tg1,
     $                   tsw,    tqf)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Carries out numerical integration and accumulates contribution
C  from current grid point into derivative Fock matrices and Hessian
C  ** OPEN SHELL **
C
C  ARGUMENTS
C
C  dft   -  method flag (NOTE: All methods include Slater exchange)
C            1 - 3 - local correlation
C           >3 - nonlocal
C  NPP   -  number of contributing (non-zero) grid points this batch
C  NBas  -  total number of basis functions
C  nbf   -  indexing array to location of "non-zero" basis functions
C  NAtoms-  number of atoms
C  Nb,Ne -  First and last component of Fock derivatives to compute
C  NBAtm -  basis functions --> atom index
C  thrsh -  threshold for neglect of contribution
C  DA    -  alpha density matrix (lower triangle)
C  DB    -  beta density matrix (lower triangle)
C  DM    -  maximum density matrix element per column
C  GDENA -  alpha density gradient at grid points (non local, dft > 3)
C  GDENB -  beta density gradient at grid points (non local, dft > 3)
C  WGHT  -  grid quadrature weights
C  pra   -  Functional derivative w.r.t. alpha density at grid points
C  prb   -  Functional derivative w.r.t. beta density at grid points
C  prara -  Functional 2nd deriv. w.r.t. alpha density at grid points
C  prbrb -  Functional 2nd deriv. w.r.t. beta density at grid points
C  prarb -  Funct. 2nd deriv. w.r.t. alpha and beta density
C  pga   -  Funct. deriv. w.r.t. alpha gradient (dft > 3)
C  pgb   -  Funct. deriv. w.r.t. beta gradient (dft > 3)
C  pgc   -  Funct. deriv. w.r.t. alpha beta gradient (dft > 3)
C  praga -  Funct. 2nd. deriv. w.r.t. alpha dens. and  grad. (dft > 3)
C  pragb -  Funct. 2nd. deriv. w.r.t. alpha dens. and  beta grad.
C           (dft > 3)
C  pragc -  Funct. 2nd. deriv. w.r.t. alpha dens. and  alpha beta grad.
C           (dft > 3)
C  prbga -  Funct. 2nd. deriv. w.r.t. beta dens. and alpha grad.
C           (dft > 3)
C  prbgb -  Funct. 2nd. deriv. w.r.t. beta dens. and  grad. (dft > 3)
C  prbgc -  Funct. 2nd. deriv. w.r.t. beta dens. and  alpha beta grad.
C           (dft > 3)
C  pgaga -  Funct. 2nd. deriv. w.r.t. alpha grad. (dft > 3)
C  pgagb -  Funct. 2nd. deriv. w.r.t. alpha and beta grad. (dft > 3)
C  pgagc -  Funct. 2nd. deriv. w.r.t. alpha and alpha beta grad.
C           (dft > 3)
C  pgbgb -  Funct. 2nd. deriv. w.r.t. beta grad. (dft > 3)
C  pgbgc -  Funct. 2nd. deriv. w.r.t. beta and alpha beta grad.
C           (dft > 3)
C  pgcgc -  Funct. 2nd. deriv. w.r.t. alpha beta grad. (dft > 3)
C  VAO   -  "non-zero" basis function values at grid points
C  VAOX  -  basis function 1st derivatives at grid points
C  VAOXX -  basis function 2nd derivatives at grid points
C  INB   -  indexing array for non-zero entries to VAO
C  VM    -  array containing maximum magnitude AO per grid point
C  INDX  -  index into contributing columns of VAO
C  DenXA -  scratch storage for in situ atomic gradient of the
C           alpha density
C  DenXB -  ditto for atomic gradient of beta density
C  GDXA  -  ditto for atomic gradient of alpha dens. gradient (dft > 3)
C  GDXB  -  ditto for atomic gradient of beta dens. gradient (dft > 3)
C  GXAA  -  ditto for atomic gradient of alpha grad. invariant (dft > 3)
C  GXBB  -  ditto for atomic gradient of beta grad. invariant (dft > 3)
C  GXAB  -  ditto for atomic gradient of alpha beta grad. invariant
C           (dft > 3)
C  SVA   -  storage for alpha coefficient of vao(i)*vao(j) in quadrature
C  SVB   -  storage for beta coefficient of vao(i)*vao(j) in quadrature
C  SWA   -  storage for alpha coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  SWB   -  storage for beta coefficient of
C             (vaox(i)*vao(j)+vao(i)*vaox(j))  (dft > 3)
C  SGA   -  storage for partial sum for hessian quadrature (dft > 3)
C  SGB   -  storage for partial sum for hessian quadrature (dft > 3)
C  SGC   -  storage for partial sum for hessian quadrature (dft > 3)
C  ICntr -  current atomic center
C
C  on exit
C
C  FDA     -  contribution to alpha derivative Fock matrices
C  FDB     -  contribution to beta derivative Fock matrices
C
C
      DIMENSION WGHT(*),VAO(*),VAOX(3,*),VAOXX(6,*),VM(*),
     $          DA(NBas*(NBas+1)/2),DB(NBas*(NBas+1)/2),DM(NBas)
      DIMENSION DenXA(3,NAtoms)
      DIMENSION DenXB(3,NAtoms)
      DIMENSION GDENA(3,*),GDENB(3,*)
      DIMENSION GDXA(3,3,NAtoms)
      DIMENSION GDXB(3,3,NAtoms)
      DIMENSION GXAA(3,NAtoms)
      DIMENSION GXBB(3,NAtoms)
      DIMENSION GXAB(3,NAtoms)
      DIMENSION SVA(3,NAtoms),SWA(3,3,NAtoms)
      DIMENSION SVB(3,NAtoms),SWB(3,3,NAtoms)
      DIMENSION sga(3,NAtoms),sgb(3,NAtoms),sgc(3,NAtoms)
      INTEGER   dft,nbf(*),INB(*),INDX(NPP),NBAtm(NBas)
      DIMENSION FDA(3,Nb:Ne,*),FDB(3,Nb:Ne,*)
      dimension pra(npp),prb(npp),pga(npp),pgb(npp),pgc(npp)
      dimension prara(npp),prbrb(npp),prarb(npp)
      dimension praga(npp),pragb(npp),pragc(npp)
      dimension prbga(npp),prbgb(npp),prbgc(npp)
      dimension pgaga(npp),pgagb(npp),pgagc(npp),pgbgb(npp),pgbgc(npp),
     $          pgcgc(npp)
C
      PARAMETER (zero=0.0d0,Two=2.0d0,Three=3.0d0,Four=4.0d0,six=6.0d0)
      parameter (twelve=12.0d0)
      PARAMETER (Half=0.5d0)
      parameter (epsi=1.0d-15,unbelpo=1.0d17)
      logical  dodenx,dodenxx,dogradxx,dox,doxx,dogxx,doval
C
C
      NAt3 = 3*NAtoms
c
      IF(dft.LE.3) THEN
C
C  local only
C
      DO 50 IP=1,NPP
      IPP = INDX(IP)
      WG=WGHT(IPP)
      ra = pra(IP)*WG
      rb = prb(IP)*WG
      rara = prara(IP)*WG
      rbrb = prbrb(IP)*WG
      rarb = prarb(IP)*WG
      VMx = VM(IPP)
      CALL ZeroIT(DenXA,NAt3)
      CALL ZeroIT(DenXB,NAt3)
c
c -- form gradient density at current grid point
c
      call secund(t1)
c
c   compute cutoffs for density derivatives
c
      vmx2=vmx*vmx
      abrr=max(abs(rara),abs(rbrb))+abs(rarb)
      dodenx=abrr.gt.epsi
      thrx=unbelpo
      if(dodenx)thrx=thrsh/abrr
      thrx1=thrx/vmx2
      abr=max(abs(ra),abs(rb))
      if(.not.dodenx)goto 25 ! should never happen!
c
c   now compute density derivatives
c
      DO 20 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        DMx2 = DM(II)*two
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        if(IAtm.lt.Nb.or.Iatm.gt.Ne)goto 20
        ValX = VAOX(1,I)
        ValY = VAOX(2,I)
        ValZ = VAOX(3,I)
        valm=max(abs(valx),abs(valy),abs(valz))
        valt=valm*dmx2*vmx
        dox=(dodenx.and.(two*valt.gt.thrx1.or.valt**2.gt.thrx))
        if(.not.dox) goto 20
        DO 18 J=nbf(IPP)+1,I
          JJ = INB(J)
          IJ = IT + JJ
          DAIJ2=two*DA(IJ)
          DBIJ2=two*DB(IJ)
          DAIJ = DAIJ2*VAO(J)
          DBIJ = DBIJ2*VAO(J)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
c
c -- atomic gradient of density
          DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
          DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
          DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
          DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
          DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
          DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
 18     CONTINUE
        DO 19 J=I+1,nbf(IPP+1)
          JJ = INB(J)
          IJ = (JJ*(JJ-1))/2 + II
          DAIJ2 = two*DA(IJ)
          DBIJ2 = two*DB(IJ)
          DAIJ = DAIJ2*VAO(J)
          DBIJ = DBIJ2*VAO(J)
          abdij=max(abs(daij),abs(dbij))
          abdijm=abdij*vmx
c
c -- atomic gradient of density
          DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
          DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
          DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
          DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
          DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
          DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
 19     CONTINUE
 20   CONTINUE
 25   continue
      call secund(t2)
      td1=td1+t2-t1
c
c -- now form the coefficients that multiply the basis functions
c    in the expression for the Fock Matrix
c    (i.e., quantities V at page 7436 of Johnson and Frisch).
c
        DO IAtm=1,NAtoms
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)
          svb(1,iatm)=rbrb*denxb(1,iatm)+rarb*denxa(1,iatm)
          svb(2,iatm)=rbrb*denxb(2,iatm)+rarb*denxa(2,iatm)
          svb(3,iatm)=rbrb*denxb(3,iatm)+rarb*denxa(3,iatm)
        EndDO
c
c    Numerical quadrature for the derivative Fock matrix.
c
c -- get maximum absolute value of sva and svb
      Call absmax(NAt3,sva,iixx,svamax)
      Call absmax(NAt3,svb,iixx,svbmax)
      svmax=max(svamax,svbmax)
c
c -- global threshold testing
c
      VMax = Max(Abr*VMx2,Abrr*svmax*VMx2)
      If(VMax.LT.thrsh) GO TO 45
c
c --  now perform quadrature
c
      thrsh1=thrsh/vmx
      DO 40 I=nbf(IPP)+1,nbf(IPP+1)
        II = INB(I)
        IT = (II*(II-1))/2
        IAtm = NBAtm(II)
        ValI = VAO(I)
        Vala = ValI*ra
        Valb = ValI*rb
        ValXa = VAOX(1,I)*ra
        ValYa = VAOX(2,I)*ra
        ValZa = VAOX(3,I)*ra
        ValXb = VAOX(1,I)*rb
        ValYb = VAOX(2,I)*rb
        ValZb = VAOX(3,I)*rb
        abval= max(abs(vala),abs(valb))
        doval= abval.gt.thrsh1
        valm = max(abs(vaox(1,I)),abs(vaox(2,I)),abs(vaox(3,I)))*abr
        Valt = MAX(abval,valm,Abs(vali)*svmax)
        IF(Valt.GT.thrsh1) THEN
          DO 30 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            JAtm = NBAtm(JJ)
            ValJ = VAO(J)
            if(abs(valj)*valm.gt.thrsh)then
              if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                FDA(1,IAtm,IJ) = FDA(1,IAtm,IJ) - ValXa*ValJ
                FDA(2,IAtm,IJ) = FDA(2,IAtm,IJ) - ValYa*ValJ
                FDA(3,IAtm,IJ) = FDA(3,IAtm,IJ) - ValZa*ValJ
                FDB(1,IAtm,IJ) = FDB(1,IAtm,IJ) - ValXb*ValJ
                FDB(2,IAtm,IJ) = FDB(2,IAtm,IJ) - ValYb*ValJ
                FDB(3,IAtm,IJ) = FDB(3,IAtm,IJ) - ValZb*ValJ
              endif
            endif
            if(doval)then
              if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - Vala*VAOX(1,J)
                FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - Vala*VAOX(2,J)
                FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - Vala*VAOX(3,J)
                FDB(1,JAtm,IJ) = FDB(1,JAtm,IJ) - Valb*VAOX(1,J)
                FDB(2,JAtm,IJ) = FDB(2,JAtm,IJ) - Valb*VAOX(2,J)
                FDB(3,JAtm,IJ) = FDB(3,JAtm,IJ) - Valb*VAOX(3,J)
              endif
            endif
c
c -- contribution of functional derivative
c
            ValIJ = ValI*ValJ
            If(Abs(ValIJ)*svmax.GT.thrsh) Then
              Do KAtm=Nb,Ne
                FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ) + ValIJ*sva(1,KAtm)
                FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ) + ValIJ*sva(2,KAtm)
                FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ) + ValIJ*sva(3,KAtm)
                FDB(1,KAtm,IJ)=FDB(1,KAtm,IJ) + ValIJ*svb(1,KAtm)
                FDB(2,KAtm,IJ)=FDB(2,KAtm,IJ) + ValIJ*svb(2,KAtm)
                FDB(3,KAtm,IJ)=FDB(3,KAtm,IJ) + ValIJ*svb(3,KAtm)
              EndDo
            EndIf
 30       CONTINUE
        ENDIF
 40   CONTINUE
 45   CONTINUE
      call secund(t3)
      tqf=tqf+t3-t2
c
 50   CONTINUE
cc
      ELSE
C
C   non-local dft
C
      DO 200 IP=1,NPP
        IPP = INDX(IP)
        WG = WGHT(IPP)
        VMx = VM(IPP)
        ra = pra(ip)*wg
        rb = prb(ip)*wg
        ga = pga(ip)*wg
        gb = pgb(ip)*wg
        gc = pgc(ip)*wg
        rara = prara(ip)*wg
        rbrb = prbrb(ip)*wg
        rarb = prarb(ip)*wg
        raga = praga(ip)*wg
        ragb = pragb(ip)*wg
        ragc = pragc(ip)*wg
        rbga = prbga(ip)*wg
        rbgb = prbgb(ip)*wg
        rbgc = prbgc(ip)*wg
        gaga = pgaga(ip)*wg
        gagb = pgagb(ip)*wg
        gagc = pgagc(ip)*wg
        gbgb = pgbgb(ip)*wg
        gbgc = pgbgc(ip)*wg
        gcgc = pgcgc(ip)*wg
c
c  density gradient at current point
c
        DAX=GDENA(1,IPP)
        DAY=GDENA(2,IPP)
        DAZ=GDENA(3,IPP)
        DBX=GDENB(1,IPP)
        DBY=GDENB(2,IPP)
        DBZ=GDENB(3,IPP)
c
c  zero out derivatives of densities and gradients
c
        CALL ZeroIT(DenXA,NAt3)
        CALL ZeroIT(GDXA,3*NAt3)
        CALL ZeroIT(DenXB,NAt3)
        CALL ZeroIT(GDXB,3*NAt3)
c
c    form the atomic gradient of density
c    and density gradient at current grid point
c
c    here it is too complex to compute cutoffs ad hoc for
c    the contribution to the various  densities, thus we
c    will only check the contributions against the global threshold.
c    in addition, we check whether the second derivatives of the
c    density gradient should be computed at all.
c
        abgmax=max(abs(ga),abs(gb),abs(gc))
        call secund(t1)
        DO 220 I=nbf(IPP)+1,nbf(IPP+1)
          II = INB(I)
          DMx2 = DM(II)*two   ! max. element of first order density
          if(dmx2.gt.epsi)then
             thtest=thrsh/(vmx*dmx2)
          else
            goto 220
          endif
          IT = (II*(II-1))/2
          IAtm = NBAtm(II)
          if(IAtm.lt.Nb.or.Iatm.gt.Ne)goto 220
          ValX = VAOX(1,I)
          ValY = VAOX(2,I)
          ValZ = VAOX(3,I)
          valm = max(abs(valx),abs(valy),abs(valz))
          ValXX = VAOXX(1,I)
          ValXY = VAOXX(2,I)
          ValYY = VAOXX(3,I)
          ValXZ = VAOXX(4,I)
          ValYZ = VAOXX(5,I)
          ValZZ = VAOXX(6,I)
          valmm1 =max(abs(valxx),abs(valxy),abs(valyy))
          valmm2 =max(abs(valxz),abs(valyz),abs(valzz))
          valmm=max(valmm1,valmm2)
          if(vmx+valmm.lt.thtest)goto 220
          DO 218 J=nbf(IPP)+1,I
            JJ = INB(J)
            IJ = IT + JJ
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2=two*DA(IJ)
            DBIJ2=two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
c
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 218      CONTINUE
          DO 219 J=I+1,nbf(IPP+1)
            JJ = INB(J)
            IJ = (JJ*(JJ-1))/2 + II
            VALJ=VAO(J)
            abvj=abs(valj)
            VALJX=VAOX(1,J)
            VALJY=VAOX(2,J)
            VALJZ=VAOX(3,J)
            abvvj=max(abs(valjx),abs(valjy),abs(valjz))
            DAIJ2= two*DA(IJ)
            DBIJ2= two*DB(IJ)
            abijt=max(abs(daij2),abs(dbij2))
            DAIJ = DAIJ2*VALJ
            DBIJ = DBIJ2*VALJ
            abij=max(abs(daij),abs(dbij))
c
c -- atomic gradient of density
            if(abij*valm.gt.thrsh)then
              DenXA(1,IAtm) = DenXA(1,IAtm) - DAIJ*ValX
              DenXA(2,IAtm) = DenXA(2,IAtm) - DAIJ*ValY
              DenXA(3,IAtm) = DenXA(3,IAtm) - DAIJ*ValZ
c
              DenXB(1,IAtm) = DenXB(1,IAtm) - DBIJ*ValX
              DenXB(2,IAtm) = DenXB(2,IAtm) - DBIJ*ValY
              DenXB(3,IAtm) = DenXB(3,IAtm) - DBIJ*ValZ
            endif
c
c -- atomic gradient of density gradient
c
            if(abijt*(abvj*valmm+abvvj*valm).gt.thrsh)then
             GDXA(1,1,IAtm)=GDXA(1,1,IAtm)-DAIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXA(1,2,IAtm)=GDXA(1,2,IAtm)-DAIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXA(1,3,IAtm)=GDXA(1,3,IAtm)-DAIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXA(2,1,IAtm)=GDXA(2,1,IAtm)-DAIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXA(2,2,IAtm)=GDXA(2,2,IAtm)-DAIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXA(2,3,IAtm)=GDXA(2,3,IAtm)-DAIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXA(3,1,IAtm)=GDXA(3,1,IAtm)-DAIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXA(3,2,IAtm)=GDXA(3,2,IAtm)-DAIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXA(3,3,IAtm)=GDXA(3,3,IAtm)-DAIJ2*(VALJ*VALZZ+VALZ*VALJZ)
c
             GDXB(1,1,IAtm)=GDXB(1,1,IAtm)-DBIJ2*(VALJ*VALXX+VALX*VALJX)
             GDXB(1,2,IAtm)=GDXB(1,2,IAtm)-DBIJ2*(VALJ*VALXY+VALY*VALJX)
             GDXB(1,3,IAtm)=GDXB(1,3,IAtm)-DBIJ2*(VALJ*VALXZ+VALZ*VALJX)
             GDXB(2,1,IAtm)=GDXB(2,1,IAtm)-DBIJ2*(VALJ*VALXY+VALX*VALJY)
             GDXB(2,2,IAtm)=GDXB(2,2,IAtm)-DBIJ2*(VALJ*VALYY+VALY*VALJY)
             GDXB(2,3,IAtm)=GDXB(2,3,IAtm)-DBIJ2*(VALJ*VALYZ+VALZ*VALJY)
             GDXB(3,1,IAtm)=GDXB(3,1,IAtm)-DBIJ2*(VALJ*VALXZ+VALX*VALJZ)
             GDXB(3,2,IAtm)=GDXB(3,2,IAtm)-DBIJ2*(VALJ*VALYZ+VALY*VALJZ)
             GDXB(3,3,IAtm)=GDXB(3,3,IAtm)-DBIJ2*(VALJ*VALZZ+VALZ*VALJZ)
            endif
 219      CONTINUE
 220    CONTINUE
      call secund(t2)
      td1=td1+t2-t1
c
c  form atomic gradient  of density gradient invariant
c
c
        DO IAtm=1,NAtoms
c   alpha alpha
          GXAA(1,IAtm)=two*(DAX*GDXA(1,1,IAtm)+
     $                   DAY*GDXA(2,1,IAtm)+DAZ*GDXA(3,1,IAtm))
          GXAA(2,IAtm)=two*(DAX*GDXA(1,2,IAtm)+
     $                   DAY*GDXA(2,2,IAtm)+DAZ*GDXA(3,2,IAtm))
          GXAA(3,IAtm)=two*(DAX*GDXA(1,3,IAtm)+
     $                   DAY*GDXA(2,3,IAtm)+DAZ*GDXA(3,3,IAtm))
c   beta beta
          GXBB(1,IAtm)=two*(DBX*GDXB(1,1,IAtm)+
     $                   DBY*GDXB(2,1,IAtm)+DBZ*GDXB(3,1,IAtm))
          GXBB(2,IAtm)=two*(DBX*GDXB(1,2,IAtm)+
     $                   DBY*GDXB(2,2,IAtm)+DBZ*GDXB(3,2,IAtm))
          GXBB(3,IAtm)=two*(DBX*GDXB(1,3,IAtm)+
     $                   DBY*GDXB(2,3,IAtm)+DBZ*GDXB(3,3,IAtm))
c   alpha beta
          GXAB(1,IAtm)=DAX*GDXB(1,1,IAtm)+DBX*GDXA(1,1,IAtm)+
     $                 DAY*GDXB(2,1,IAtm)+DBY*GDXA(2,1,IAtm)+
     $                 DAZ*GDXB(3,1,IAtm)+DBZ*GDXA(3,1,IAtm)
          GXAB(2,IAtm)=DAX*GDXB(1,2,IAtm)+DBX*GDXA(1,2,IAtm)+
     $                 DAY*GDXB(2,2,IAtm)+DBY*GDXA(2,2,IAtm)+
     $                 DAZ*GDXB(3,2,IAtm)+DBZ*GDXA(3,2,IAtm)
          GXAB(3,IAtm)=DAX*GDXB(1,3,IAtm)+DBX*GDXA(1,3,IAtm)+
     $                 DAY*GDXB(2,3,IAtm)+DBY*GDXA(2,3,IAtm)+
     $                 DAZ*GDXB(3,3,IAtm)+DBZ*GDXA(3,3,IAtm)
        EndDO
        call secund(t3)
        tg1=tg1+t3-t2
c
c -- now form the coefficients that multiply the basis functions and
c    their derivatives in the expression for the Fock Matrix
c    (i.e., quantities V, W and X at page 7436 of Johnson and Frisch).
C    and also compute partial sums used in Hessian quadrature
c
        ga2= ga+ga
        gb2= gb+gb
        DAX2=DAX+DAX
        DAY2=DAY+DAY
        DAZ2=DAZ+DAZ
        DBX2=DBX+DBX
        DBY2=DBY+DBY
        DBZ2=DBZ+DBZ
c
c     coeff. X
c
        SXXA=ga2*DAX+gc*DBX
        SXYA=ga2*DAY+gc*DBY
        SXZA=ga2*DAZ+gc*DBZ
        SXXB=gb2*DBX+gc*DAX
        SXYB=gb2*DBY+gc*DAY
        SXZB=gb2*DBZ+gc*DAZ
c
        DO IAtm=1,NAtoms
c
c     coeff. V
c
          sva(1,iatm)=rara*denxa(1,iatm)+rarb*denxb(1,iatm)+
     $        raga*gxaa(1,iatm)+ragb*gxbb(1,iatm)+ragc*gxab(1,iatm)
          sva(2,iatm)=rara*denxa(2,iatm)+rarb*denxb(2,iatm)+
     $        raga*gxaa(2,iatm)+ragb*gxbb(2,iatm)+ragc*gxab(2,iatm)
          sva(3,iatm)=rara*denxa(3,iatm)+rarb*denxb(3,iatm)+
     $        raga*gxaa(3,iatm)+ragb*gxbb(3,iatm)+ragc*gxab(3,iatm)
          svb(1,iatm)=rarb*denxa(1,iatm)+rbrb*denxb(1,iatm)+
     $        rbga*gxaa(1,iatm)+rbgb*gxbb(1,iatm)+rbgc*gxab(1,iatm)
          svb(2,iatm)=rarb*denxa(2,iatm)+rbrb*denxb(2,iatm)+
     $        rbga*gxaa(2,iatm)+rbgb*gxbb(2,iatm)+rbgc*gxab(2,iatm)
          svb(3,iatm)=rarb*denxa(3,iatm)+rbrb*denxb(3,iatm)+
     $        rbga*gxaa(3,iatm)+rbgb*gxbb(3,iatm)+rbgc*gxab(3,iatm)
c
c  partial sums for Hessian
c
          sga(1,iatm)=raga*denxa(1,iatm)+rbga*denxb(1,iatm)+
     $        gaga*gxaa(1,iatm)+gagb*gxbb(1,iatm)+gagc*gxab(1,iatm)
          sga(2,iatm)=raga*denxa(2,iatm)+rbga*denxb(2,iatm)+
     $        gaga*gxaa(2,iatm)+gagb*gxbb(2,iatm)+gagc*gxab(2,iatm)
          sga(3,iatm)=raga*denxa(3,iatm)+rbga*denxb(3,iatm)+
     $        gaga*gxaa(3,iatm)+gagb*gxbb(3,iatm)+gagc*gxab(3,iatm)
c
          sgb(1,iatm)=ragb*denxa(1,iatm)+rbgb*denxb(1,iatm)+
     $        gagb*gxaa(1,iatm)+gbgb*gxbb(1,iatm)+gbgc*gxab(1,iatm)
          sgb(2,iatm)=ragb*denxa(2,iatm)+rbgb*denxb(2,iatm)+
     $        gagb*gxaa(2,iatm)+gbgb*gxbb(2,iatm)+gbgc*gxab(2,iatm)
          sgb(3,iatm)=ragb*denxa(3,iatm)+rbgb*denxb(3,iatm)+
     $        gagb*gxaa(3,iatm)+gbgb*gxbb(3,iatm)+gbgc*gxab(3,iatm)
c
          sgc(1,iatm)=ragc*denxa(1,iatm)+rbgc*denxb(1,iatm)+
     $        gagc*gxaa(1,iatm)+gbgc*gxbb(1,iatm)+gcgc*gxab(1,iatm)
          sgc(2,iatm)=ragc*denxa(2,iatm)+rbgc*denxb(2,iatm)+
     $        gagc*gxaa(2,iatm)+gbgc*gxbb(2,iatm)+gcgc*gxab(2,iatm)
          sgc(3,iatm)=ragc*denxa(3,iatm)+rbgc*denxb(3,iatm)+
     $        gagc*gxaa(3,iatm)+gbgc*gxbb(3,iatm)+gcgc*gxab(3,iatm)
c
c     coeff. W
c
          swa(1,1,iatm)=ga2*gdxa(1,1,iatm)+gc*gdxb(1,1,Iatm)+
     $        DAX2*sga(1,iatm)+DBX*sgc(1,iatm)
          swa(1,2,iatm)=ga2*gdxa(1,2,iatm)+gc*gdxb(1,2,Iatm)+
     $        DAX2*sga(2,iatm)+DBX*sgc(2,iatm)
          swa(1,3,iatm)=ga2*gdxa(1,3,iatm)+gc*gdxb(1,3,Iatm)+
     $        DAX2*sga(3,iatm)+DBX*sgc(3,iatm)
          swa(2,1,iatm)=ga2*gdxa(2,1,iatm)+gc*gdxb(2,1,Iatm)+
     $        DAY2*sga(1,iatm)+DBY*sgc(1,iatm)
          swa(2,2,iatm)=ga2*gdxa(2,2,iatm)+gc*gdxb(2,2,Iatm)+
     $        DAY2*sga(2,iatm)+DBY*sgc(2,iatm)
          swa(2,3,iatm)=ga2*gdxa(2,3,iatm)+gc*gdxb(2,3,Iatm)+
     $        DAY2*sga(3,iatm)+DBY*sgc(3,iatm)
          swa(3,1,iatm)=ga2*gdxa(3,1,iatm)+gc*gdxb(3,1,Iatm)+
     $        DAZ2*sga(1,iatm)+DBZ*sgc(1,iatm)
          swa(3,2,iatm)=ga2*gdxa(3,2,iatm)+gc*gdxb(3,2,Iatm)+
     $        DAZ2*sga(2,iatm)+DBZ*sgc(2,iatm)
          swa(3,3,iatm)=ga2*gdxa(3,3,iatm)+gc*gdxb(3,3,Iatm)+
     $        DAZ2*sga(3,iatm)+DBZ*sgc(3,iatm)
c
          swb(1,1,iatm)=gb2*gdxb(1,1,iatm)+gc*gdxa(1,1,Iatm)+
     $        DBX2*sgb(1,iatm)+DAX*sgc(1,iatm)
          swb(1,2,iatm)=gb2*gdxb(1,2,iatm)+gc*gdxa(1,2,Iatm)+
     $        DBX2*sgb(2,iatm)+DAX*sgc(2,iatm)
          swb(1,3,iatm)=gb2*gdxb(1,3,iatm)+gc*gdxa(1,3,Iatm)+
     $        DBX2*sgb(3,iatm)+DAX*sgc(3,iatm)
          swb(2,1,iatm)=gb2*gdxb(2,1,iatm)+gc*gdxa(2,1,Iatm)+
     $        DBY2*sgb(1,iatm)+DAY*sgc(1,iatm)
          swb(2,2,iatm)=gb2*gdxb(2,2,iatm)+gc*gdxa(2,2,Iatm)+
     $        DBY2*sgb(2,iatm)+DAY*sgc(2,iatm)
          swb(2,3,iatm)=gb2*gdxb(2,3,iatm)+gc*gdxa(2,3,Iatm)+
     $        DBY2*sgb(3,iatm)+DAY*sgc(3,iatm)
          swb(3,1,iatm)=gb2*gdxb(3,1,iatm)+gc*gdxa(3,1,Iatm)+
     $        DBZ2*sgb(1,iatm)+DAZ*sgc(1,iatm)
          swb(3,2,iatm)=gb2*gdxb(3,2,iatm)+gc*gdxa(3,2,Iatm)+
     $        DBZ2*sgb(2,iatm)+DAZ*sgc(2,iatm)
          swb(3,3,iatm)=gb2*gdxb(3,3,iatm)+gc*gdxa(3,3,Iatm)+
     $        DBZ2*sgb(3,iatm)+DAZ*sgc(3,iatm)
        EndDO
        call secund(t4)
        tsw=tsw+t4-t3
c
c    we are ready to perform the quadrature for the derivative
c    Fock Matrix.
c
c     get the maximum absolute value of coefficients V and W
c
        Call absmax(NAt3,sva,iixx,svamax)
        Call absmax(NAt3,svb,iixx,svbmax)
        svmax=max(svamax,svbmax)
        Call absmax(3*NAt3,swa,iixx,swamax)
        Call absmax(3*NAt3,swb,iixx,swbmax)
        swmax=max(swamax,swbmax)
        tswmax=Three*swmax
c
c -- global threshold testing
c
        abr=max(abs(ra),abs(rb))
        vmx2=vmx*vmx
        SMax = max(Abs(SXXA+SXYA+SXZA),abs(SXXB+SXYB+SXZB))
        VMax3 = (svmax+tswmax+tswmax)*VMx2
        VMax1 = (Abr+two*smax)*VMx2
        vmax=max(vmax1,vmax3)
        If(VMax.LT.thrsh)  GO TO 245
c
c -- numerical quadrature  for derivative Fock Matrix
c
        DO 240 I=nbf(IPP)+1,nbf(IPP+1)
          ValI = VAO(I)
          ValIX = VAOX(1,I)
          ValIY = VAOX(2,I)
          ValIZ = VAOX(3,I)
          XVIA=SXXA*ValIX+SXYA*ValIY+SXZA*ValIZ+ra*ValI
          XVIB=SXXB*ValIX+SXYB*ValIY+SXZB*ValIZ+rb*ValI
          SXXIA=SXXA*ValI
          SXYIA=SXYA*ValI
          SXZIA=SXZA*ValI
          SXXIB=SXXB*ValI
          SXYIB=SXYB*ValI
          SXZIB=SXZB*ValI
          XXVXA=SXXA*VAOXX(1,I)+SXYA*VAOXX(2,I)+SXZA*VAOXX(4,I)+ra*VALIX
          XXVYA=SXXA*VAOXX(2,I)+SXYA*VAOXX(3,I)+SXZA*VAOXX(5,I)+ra*VALIY
          XXVZA=SXXA*VAOXX(4,I)+SXYA*VAOXX(5,I)+SXZA*VAOXX(6,I)+ra*VALIZ
          XXVXB=SXXB*VAOXX(1,I)+SXYB*VAOXX(2,I)+SXZB*VAOXX(4,I)+rb*VALIX
          XXVYB=SXXB*VAOXX(2,I)+SXYB*VAOXX(3,I)+SXZB*VAOXX(5,I)+rb*VALIY
          XXVZB=SXXB*VAOXX(4,I)+SXYB*VAOXX(5,I)+SXZB*VAOXX(6,I)+rb*VALIZ
          valm=abs(vali)
          valm1=max(abs(xvia),abs(xvib))
          valmx=max(abs(valix),abs(valiy),abs(valiz))
          valmx1a=max(abs(xxvxa),abs(xxvya),abs(xxvza))
          valmx1b=max(abs(xxvxb),abs(xxvyb),abs(xxvzb))
          valmx1=max(valmx1a,valmx1b)
          val1=(valmx1+smax*valmx)*vmx
          smaxmx=smax*valm*vmx
          val2=valm1*vmx+smaxmx
          val3=(svmax*valm+tswmax*(valmx+valm))*vmx
          valtest=max(val1,val2,val3)
          if(valtest.GT.thrsh)then
            II = INB(I)
            IT = (II*(II-1))/2
            IAtm = NBAtm(II)
            DO 230 J=nbf(IPP)+1,I
              JJ = INB(J)
              IJ = IT + JJ
              JAtm = NBAtm(JJ)
              ValJ = VAO(J)
              ValJX = VAOX(1,J)
              ValJY = VAOX(2,J)
              ValJZ = VAOX(3,J)
              valmxj=max(abs(valjx),abs(valjy),abs(valjz))
              XVJA=SXXA*ValJX+SXYA*ValJY+SXZA*ValJZ
              XVJB=SXXB*ValJX+SXYB*ValJY+SXZB*ValJZ
              abxvj=max(abs(xvja),abs(xvjb))
              if(valmx1*abs(valj)+valmx*abxvj.gt.thrsh)then
                if(IAtm.ge.Nb.and.IAtm.le.Ne)then
                  FDA(1,IAtm,IJ)=FDA(1,IAtm,IJ)-XXVXA*ValJ-ValIX*XVJA
                  FDA(2,IAtm,IJ)=FDA(2,IAtm,IJ)-XXVYA*ValJ-ValIY*XVJA
                  FDA(3,IAtm,IJ)=FDA(3,IAtm,IJ)-XXVZA*ValJ-ValIZ*XVJA
                  FDB(1,IAtm,IJ)=FDB(1,IAtm,IJ)-XXVXB*ValJ-ValIX*XVJB
                  FDB(2,IAtm,IJ)=FDB(2,IAtm,IJ)-XXVYB*ValJ-ValIY*XVJB
                  FDB(3,IAtm,IJ)=FDB(3,IAtm,IJ)-XXVZB*ValJ-ValIZ*XVJB
                endif
              endif
              if(valm1*valmxj+smaxmx.gt.thrsh)then
                if(JAtm.ge.Nb.and.JAtm.le.Ne)then
                  FDA(1,JAtm,IJ) = FDA(1,JAtm,IJ) - XVIA*ValJX -
     $            VAOXX(1,J)*SXXIA - VAOXX(2,J)*SXYIA - VAOXX(4,J)*SXZIA
                  FDA(2,JAtm,IJ) = FDA(2,JAtm,IJ) - XVIA*ValJY -
     $            VAOXX(2,J)*SXXIA - VAOXX(3,J)*SXYIA - VAOXX(5,J)*SXZIA
                  FDA(3,JAtm,IJ) = FDA(3,JAtm,IJ) - XVIA*ValJZ -
     $            VAOXX(4,J)*SXXIA - VAOXX(5,J)*SXYIA - VAOXX(6,J)*SXZIA
c
                  FDB(1,JAtm,IJ) = FDB(1,JAtm,IJ) - XVIB*ValJX -
     $            VAOXX(1,J)*SXXIB - VAOXX(2,J)*SXYIB - VAOXX(4,J)*SXZIB
                  FDB(2,JAtm,IJ) = FDB(2,JAtm,IJ) - XVIB*ValJY -
     $            VAOXX(2,J)*SXXIB - VAOXX(3,J)*SXYIB - VAOXX(5,J)*SXZIB
                  FDB(3,JAtm,IJ) = FDB(3,JAtm,IJ) - XVIB*ValJZ -
     $            VAOXX(4,J)*SXXIB - VAOXX(5,J)*SXYIB - VAOXX(6,J)*SXZIB
                endif
              endif
              VIJ=VALI*VALJ
              VIJX=VALIX*VALJ+VALJX*VALI
              VIJY=VALIY*VALJ+VALJY*VALI
              VIJZ=VALIZ*VALJ+VALJZ*VALI
              if(svmax*abs(vij)+swmax*abs(vijx+vijy+vijz).GT.thrsh)then
                Do KAtm=Nb,Ne
                  FDA(1,KAtm,IJ)=FDA(1,KAtm,IJ)+(SVA(1,KAtm)*VIJ+
     $                           SWA(1,1,KAtm)*VIJX+SWA(2,1,KAtm)*VIJY+
     $                           SWA(3,1,KAtM)*VIJZ)
                  FDA(2,KAtm,IJ)=FDA(2,KAtm,IJ)+(SVA(2,KAtm)*VIJ+
     $                           SWA(1,2,KAtm)*VIJX+SWA(2,2,KAtm)*VIJY+
     $                           SWA(3,2,KAtM)*VIJZ)
                  FDA(3,KAtm,IJ)=FDA(3,KAtm,IJ)+(SVA(3,KAtm)*VIJ+
     $                           SWA(1,3,KAtm)*VIJX+SWA(2,3,KAtm)*VIJY+
     $                           SWA(3,3,KatM)*VIJZ)
c
                  FDB(1,KAtm,IJ)=FDB(1,KAtm,IJ)+(SVB(1,KAtm)*VIJ+
     $                           SWB(1,1,KAtm)*VIJX+SWB(2,1,KAtm)*VIJY+
     $                           SWB(3,1,KAtM)*VIJZ)
                  FDB(2,KAtm,IJ)=FDB(2,KAtm,IJ)+(SVB(2,KAtm)*VIJ+
     $                           SWB(1,2,KAtm)*VIJX+SWB(2,2,KAtm)*VIJY+
     $                           SWB(3,2,KAtM)*VIJZ)
                  FDB(3,KAtm,IJ)=FDB(3,KAtm,IJ)+(SVB(3,KAtm)*VIJ+
     $                           SWB(1,3,KAtm)*VIJX+SWB(2,3,KAtm)*VIJY+
     $                           SWB(3,3,KatM)*VIJZ)
                EndDo
              EndIf
 230        CONTINUE
          ENDIF
 240    CONTINUE
 245    CONTINUE
        call secund(t5)
        tqf=tqf+t5-t4
c
 200  CONTINUE
      ENDIF
C
      RETURN
      END

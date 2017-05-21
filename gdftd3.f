      SUBROUTINE GDFTD3(IDft,   IDisp,  IPrnt,  NAtoms, IAN,
     $                  XC,     VALUES, MaxEL,  MaxC,   R0AB,
     $                  C6AB,   CN,     MXC,    ITMP,   noabc,
     $                  tz,     DispE,  GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Implementation of Grimme's dispersion correction to DFT energies
C  Depends on atomic number and geometry only and has been parametrized
C  for several DFT functionals
C
C  ** This routine computes the dispersion energy + gradient **
C
C  ARGUMENTS
C
C  IDft    -  DFT flag
C              0 - no dft contribution (i.e., Hartree-Fock only)
C              1 - Slater exchange only
C              2 - VWN (Slater exchange + local RPA correlation)
C              3 - VWN5 (ditto but RECOMMENDED Ceperley-Alder fit)
C              4 - Becke 88 exchange
C              5 - Becke 88 + VWN (local correlation)
C              6 - Becke 88 + VWN5 (local correlation)
C              7 - Becke 88 + Perdew 86
C              8 - Becke 88 + Perdew-Wang 91
C              9 - Becke 88 + LYP (correlation)
C             10 - Becke 88 + VWN5 (local) + Perdew 86 (nonlocal)
C             11 - Handy/Cohen exchange
C             12 - Handy/Cohen + VWN
C             13 - Handy/Cohen + VWN5
C             14 - Handy/Cohen + Perdew 86
C             15 - Handy/Cohen + Perdew-Wang 91
C             16 - Handy/Cohen + LYP
C             17 - PW91 exchange + correlation
C             18 - PBE exchange + correlation
C             19 - O3LYP   (VWN5 + OPTX     + LYP  + HF exchange)
C             20 - B3LYP   (VWN  + Becke 88 + LYP  + HF exchange)
C             21 - B3PW91  (VWN5 + Becke 88 + PW91 + HF exchange)
C             22 - WAH     (VWN5 + Becke 88 + LYP  + HF exchange;
C                           modified B3LYP specifically for NMR)
C             23 - user-defined functional
C -- The following functionals are self-contained and were downloaded from the web
C             24 - B97    (Becke 97 exchange-correlation functional)
C             25 - B97-1  ( ditto, but reparametrized)
C             26 - B97-2  ( ditto, but reparametrized again)
C             27 - HCTH   (Hamprecht-Cohen-Tozer-Handy)
C                  (this is similar in form to B97, but has no HF exchange)
C    ** NOTE:  NOT ALL OF THE ABOVE FUNCTIONALS HAVE BEEN PARAMETRIZED **
C  IDisp   -  Type of dispersion correction
C             -1 - No dispersion  (This routine should NOT be called)
C              0 - default dispersion  (=4)
C              1 - ????
C              2 - DFT-D2
C              3 - DFT-D3
C              4 - DFT-D3 with Becke-Johnson finite-damping   (default)
C              50 - User-defined parametrization using default dispersion
C              52 - User-defined parametrization using DFT-D2
C              53 - User-defined parametrization using DFT-D3
C              54 - User-defined parametrization using DFT-D3 + BJ damping
C  IPrnt   -  print flag
C  NAtoms  -  number of atoms
C  IAN     -  atomic numbers
C  XC      -  Cartesian coordinates
C  VALUES  -  user-defined values of the five main dispersion parameters
C              1 - rs6   2 - s18   3 - rs18   4 - s6   5 - alp
C  MaxEL   -  maximum number of elements for which method is parametrized (94)
C  MaxC    -  maximum coordination number
C  R0AB    -  storage for cutoff radii for all element pairs
C  C6AB    -  storage for c6 parameters for all element pairs
C  CN      -  storage for coordination number for all elements
C  MXC     -  maximum coordination number for each element
C  ITMP    -  temporary integer storage
C  noabc   -  logical flag for use of third-order term
C  tz      -  logical flag for TZ-type basis parametrization
C  DispE   -  previously computed dispersion energy
C  GC      -  on exit contains dispersion gradient
C
C  ----------------------------------------------------------------------------------
C  references
C
C    "A consistent and accurate ab initio parameterization of density functional
C     dispersion correction (DFT-D) for the 94 elements H-Pu"
C     S. Grimme, J. Antony, S. Ehrlich and H. Krieg
C     J. Chem. Phys, 132 (2010) 154104
C
C    "Effect of the damping function in dispersion corrected density functional theory"
C     S. Grimme, S. Ehrlich and L. Goerigk    J. Comput. Chem, 32 (2011) 1456
C ------------------------------------------------------------------------------------
C
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c  coversion factors
      parameter (autokcal=627.509541d0)
      parameter (c6conv=1.0d-3/2625.4999d0/(0.052917726d0)**6)   ! J/mol nm**6 --> au
c
      REAL*8 XC(3,NAtoms),GC(3,NAtoms),VALUES(5)
      REAL*8 R0AB(MaxEL,MaxEL),C6AB(MaxEL,MaxEL,MaxC,MaxC,3),CN(NAtoms)
      INTEGER IAN(NAtoms),MXC(MaxEL),ITMP(MaxEL)
      REAL*8 RCOV(94),DCOV(94)              ! covalent radii
      REAL*8 R2R4(94),D2D4(94)              ! atomic (r2/r4) values
      REAL*8 k1,k2,k3
      Logical Error,noabc,numgrad,tz
      real*8 dum6(86)
c
      PARAMETER (k1=16.0d0,k2=4.0d0/3.0d0,k3=-4.0d0)
c
c  PBE0/def2-QZVP atomic values 
      data R2R4/
     $  8.0589,  3.4698, 29.0974, 14.8517, 11.8799,  7.8715,  5.5588,
     $  4.7566,  3.8025,  3.1036, 26.1552, 17.2304, 17.7210, 12.7442,
     $  9.5361,  8.1652,  6.7463,  5.6004, 29.2012, 22.3934, 19.0598,
     $ 16.8590, 15.4023, 12.5589, 13.4788, 12.2309, 11.2809, 10.5569,
     $ 10.1428,  9.4907, 13.4606, 10.8544,  8.9386,  8.1350,  7.1251,
     $  6.1971, 30.0162, 24.4103, 20.3537, 17.4780, 13.5528, 11.8451,
     $ 11.0355, 10.1997,  9.5414,  9.0061,  8.6417,  8.9975, 14.0834,
     $ 11.8333, 10.0179,  9.3844,  8.4110,  7.5152, 32.7622, 27.5708,
     $ 23.1671, 21.6003, 20.9615, 20.4562, 20.1010, 19.7475, 19.4828,
     $ 15.6013, 19.2362, 17.4717, 17.8321, 17.4237, 17.1954, 17.1631,
     $ 14.5716, 15.8758, 13.8989, 12.4834, 11.4421, 10.2671,  8.3549,
     $  7.8496,  7.3278,  7.4820, 13.5124, 11.6554, 10.0959,  9.7340,
     $  8.8584,  8.0125, 29.8135, 26.3157, 19.1885, 15.8542, 16.1305,
     $ 15.6161, 15.1226, 16.1576 /                                       
c
c  covalent radii (taken from Pyykko and Atsumi, Chem. Eur. J. 15, 2009, 188-197)
c  values for metals decreased by 10 %
      data RCOV/
     $  0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
     $  1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,
     $  1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
     $  1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,
     $  1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,
     $  1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,
     $  1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
     $  1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
     $  1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,
     $  1.52, 1.53, 1.54, 1.55 /
c
c
c  scale and convert to au
cc       Call VScal(MaxEL,k2*ANTOAU,RCOV)
      DSkal = k2*ANTOAU
      Do I=1,MaxEL
      DCOV(I) = DSkal*RCOV(I)
      EndDO
c
c  Cutoff r^2 thresholds for the gradient in bohr^2.
c  rthr influences the O(N^2) part of the gradient.
c  rthr2 influences the O(N^3) part of the gradient. When using
c  dftd3 in combination with semi-empirical methods or FFs, and large
c  (>1000 atoms) systems, rthr2 is crucial for speed:
c  Recommended values are 20^2 to 25^2 bohr.
c
      rthr=20000.0d0
      rthr2=1600.0d0

c -- set radii
      call setr0ab(MaxEL,1.0d0/ANTOAU,R0AB)
cc      call prntmat(maxel,maxel,maxel,r0ab)
c
c -- load C6AB values
      call copyc6(MaxC,MaxEL,C6AB,MXC)         
c
c  the analytical E(3) grad is not available yet
      if(.not.noabc) numgrad=.true.
c     
c  normalize rthr2 to r0HH to achieve more reasonable cutoffs
       rthr2=rthr2/R0AB(1,1)
c
c -- set parameters for functionals
        rs6 = VALUES(1)
        s18 = VALUES(2)
        rs18 = VALUES(3)
        s6 = VALUES(4)
        alp = VALUES(5)
cccccccccc
      write(6,*) ' Dispersion parameter values are:'
      write(6,*) '  IDisp is:  ',IDisp
      write(6,*) '  rs6:  ',rs6
      write(6,*) '  s18:  ',s18
      write(6,*) '  rs18: ',rs18
      write(6,*) '  s6:   ',s6
      write(6,*) '  alp:  ',alp
cccccccccc

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C all calculations start here
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

c  CNs for output
      call ncoord(NAtoms,DCOV,IAN,XC,CN)
c
      If(IDisp.eq.2) Then
        If(IPrnt.GT.2) WRITE(6,*) '  Using DFT-D2 Parameters'
        call loadoldpar(1.0d0/ANTOAU,MaxEL,MaxC,C6AB,R0AB,dum6)
c  number of CNs for each element
        Do I=1,MaxEL
        MXC(I) = 1
        EndDo
c -- convert to au
        CALL VScal(MaxEL*MaxEL*MaxC*MaxC*3,c6conv,C6AB)
      EndIf
c
c  scale r4/r2 values of the atoms by sqrt(Z) 
c  sqrt is also globally close to optimum
c  together with the factor 1/2 this yield reasonable
c  c8 for he, ne and ar. for larger Z, C8 becomes too large
c  which effectively mimics higher R^n terms neglected due
c  to stability reasons
c
      Do I=1,MaxEL
      dum = 0.5d0*R2R4(I)*SQRT(DFloat(I))
c -- store it as sqrt because the geom. av. is taken
      D2D4(I)=SQRT(dum)                         
      EndDo
c
c -- set fixed or dependent parameters
      rs8  = rs18       
      rs10 = rs18
      alp6 = alp
      alp8 = alp+2.0d0
      alp10= alp8+2.0d0 
c  note: if version=4 (Becke-Johnson), a1=rs6 and a2=rs18
c        and alp* have no meaning
c
c  check if all parameters have been loaded and are reasonable
      Error = .False.
      DO I=1,NAtoms-1
      IAT = IAN(I)
      DO J=I+1,NAtoms
      JAT = IAN(J)
      If(R0AB(JAT,IAT).LT.0.1) Then
        WRITE(6,*) ' Radius Missing for atom pair: ',I,J
        Error = .True.
      EndIf
      If(IDisp.EQ.2) Then
        c6 = C6AB(JAT,IAT,1,1,1)
      Else
        call getc6(MaxC,MaxEL,C6AB,MXC,IAT,JAT,CN(I),CN(J),c6)
      EndIf
      If(c6.lt.1.d-6) Then
        WRITE(6,*) ' c6 Missing for atom pair: ',I,J
        Error = .True.
      EndIf
      EndDo
      EndDo
c
      If(Error) call nerror(13,'FORCE Module',
     $           'Missing Parameters in DFT Dispersion',0,0)

ccccccccccccccccccccccccccc
c  gradient
ccccccccccccccccccccccccccc

      Call ZeroIT(GC,3*NAtoms)
      call elapsec(t1)
      call gdisp(MaxEL,  MaxC,   NAtoms, XC,     IAN,
     $           C6AB,   MXC,    D2D4,   R0AB,   DCOV,
     $           s6,     s18,    rs6,    rs8,    rs10,
     $           alp6,   alp8,   alp10,  noabc,  rthr,
     $          numgrad, IDisp,  IPrnt,  GC,     dispg,
     $           gnorm,  rthr2)
      call elapsec(t2)
      If(IPrnt.GT.3)
     $   WRITE(6,*) ' time for dispersion gradient: ',t2-t1
c -- check if gdisp yields same energy as edisp
      write(6,*) ' Back from <gdisp>  dispg is: ',dispg
      If(ABS((DispE-dispg)/DispE).GT.1.0d-5) then
        write(6,*) DispE,dispg
        call nerror(14,'SCF Module',
     $    'Internal error in DFT Dispersion ',0,0)
      endif
C
      RETURN
      END

cccccccccccccccccccccccccc
c  compute gradient
cccccccccccccccccccccccccc

      subroutine gdisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .                 s6,s18,rs6,rs8,rs10,alp6,alp8,alp10,noabc,rthr,
     .                 num,version,IPrnt,g,disp,gnorm,rthr2)
      implicit none  
      integer n,iz(*),max_elem,maxc,version,IPrnt,mxc(max_elem)
      real*8 xyz(3,*),r0ab(max_elem,max_elem),r2r4(*)
      real*8 c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8 g(3,*),s6,s18,rcov(max_elem)
      real*8 rs6,rs8,rs10,alp10,alp8,alp6,a1,a2,r2ik        
      logical noabc,num
 
      integer iat,jat,i,j,kat
      real*8 R0,C6,alp,R42,disp,x1,y1,z1,x2,y2,z2,rr,e6abc  
      real*8 dx,dy,dz,r2,r,r4,r6,r8,r10,r12,t6,t8,t10,damp1
      real*8 damp6,damp8,damp10,e6,e8,e10,e12,gnorm,tmp1
      real*8 s10,s8,gC6(3),term,step,dispr,displ,r235,tmp2
      real*8 cn(n),gx1,gy1,gz1,gx2,gy2,gz2,rthr,c8,rthr2
      real*8 dcn2(3,n),dcn3(3,n,n),rthr3

c this is the crucial threshold to reduce the N^3 to an
c effective N^2. 

      rthr3=rthr2
c      write(*,*)'rthr=',rthr,'rthr2=',rthr2,'rthr3=',rthr3

      if(IPrnt.GT.3) write(6,*) 
      if(version.eq.2)then
      if(IPrnt.GT.3)write(6,*) 'doing analytical gradient O(N^2) ...'
      disp=0
      do iat=1,n-1
         do jat=iat+1,n
            R0=r0ab(iz(jat),iz(iat))*rs6
            dx=(xyz(1,iat)-xyz(1,jat))
            dy=(xyz(2,iat)-xyz(2,jat))
            dz=(xyz(3,iat)-xyz(3,jat))
            r2  =dx*dx+dy*dy+dz*dz             
c           if(r2.gt.rthr) cycle
            r235=r2**3.5                       
            r   =sqrt(r2)
            damp6=exp(-alp6*(r/R0-1.0d0))
            damp1=1.+damp6           
            c6=c6ab(iz(jat),iz(iat),1,1,1)*s6
            tmp1=damp6/(damp1*damp1*r235*R0)
            tmp2=6./(damp1*r*r235)
            gx1=alp6* dx*tmp1-tmp2*dx
            gx2=alp6*(-dx)*tmp1+tmp2*dx
            gy1=alp6* dy*tmp1-tmp2*dy
            gy2=alp6*(-dy)*tmp1+tmp2*dy
            gz1=alp6* dz*tmp1-tmp2*dz
            gz2=alp6*(-dz)*tmp1+tmp2*dz
            g(1,iat)=g(1,iat)-gx1*c6
            g(2,iat)=g(2,iat)-gy1*c6
            g(3,iat)=g(3,iat)-gz1*c6
            g(1,jat)=g(1,jat)-gx2*c6  
            g(2,jat)=g(2,jat)-gy2*c6      
            g(3,jat)=g(3,jat)-gz2*c6      
            disp=disp+c6*(1./damp1)/r2**3
         enddo
      enddo
      disp=-disp
      goto 999
      endif

cNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNNN
      if(num) then
      if(IPrnt.GT.3)write(6,*) 'doing numerical gradient O(N^3) ...'

      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,
     .     e6,e8,e10,e12,e6abc)
      disp=-s6*e6-s18*e8-s6*e6abc

      step=2.d-5     

      do i=1,n
      do j=1,3
      xyz(j,i)=xyz(j,i)+step        
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,
     .     e6,e8,e10,e12,e6abc)
      dispr=-s6*e6-s18*e8-s6*e6abc
      xyz(j,i)=xyz(j,i)-2*step      
      call edisp(max_elem,maxc,n,xyz,iz,c6ab,mxc,r2r4,r0ab,rcov,
     .     rs6,rs8,rs10,alp6,alp8,alp10,version,noabc,rthr,
     .     e6,e8,e10,e12,e6abc)
      displ=-s6*e6-s18*e8-s6*e6abc
      g(j,i)=0.5*(dispr-displ)/step  
      xyz(j,i)=xyz(j,i)+step        
      enddo
      enddo

      else

      if(IPrnt.GT.3)write(6,*) 'doing analytical gradient O(N^3) ...'
c precompute for analytical part
      call ncoorda(n,rcov,iz,xyz,cn,dcn2,dcn3)

c 333333333333333333333333333333333333333333333333333333333333333333333333333
c standard correction
      if (version.eq.3) then
      s8 =s18
      s10=s18

      disp=0

      do iat=1,n
         x1=xyz(1,iat)
         y1=xyz(2,iat)
         z1=xyz(3,iat)
         do jat=1,n
            if(iat.eq.jat) cycle
            x2=xyz(1,jat)
            y2=xyz(2,jat)
            z2=xyz(3,jat)
            dx = (x1-x2)**2
            dy = (y1-y2)**2
            dz = (z1-z2)**2
            r2 = dx+dy+dz
CTHR
            if(r2.gt.rthr) cycle

            R0=r0ab(iz(jat),iz(iat))
c stored as sqrt
            r42=r2r4(iz(iat))*r2r4(iz(jat))
            call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                    cn(iat),cn(jat),C6)
           
c analytically
            call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              iat,jat,iat,gC6)

      r = dsqrt(r2)
      t6 = (r/(rs6*R0))**(-alp6)
      damp6 =1.d0/( 1.d0+6.d0*t6 )
      t8 = (r/(rs8*R0))**(-alp8)
      damp8 =1.d0/( 1.d0+6.d0*t8 )
 
      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
      r10 = r4*r6
c     r12 = r6**2

      dx = 2.D0*x1-2.D0*x2

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dx
     &       -1.D0*damp6*s6*gC6(1)/r6
     &       +3.D0*damp6*s6*C6/r8*dx
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dx
     &       -3.D0*damp8*s8*gC6(1)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dx
c    &       -11.025D0*gC6(1)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dx
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dx
      g(1,iat)=g(1,iat)+term

      dy = 2.D0*y1-2.D0*y2       

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dy
     &       -1.D0*damp6*s6*gC6(2)/r6
     &       +3.D0*damp6*s6*C6/r8*dy
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dy
     &       -3.D0*damp8*s8*gC6(2)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dy
c    &       -11.025D0*gC6(2)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dy
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dy
      g(2,iat)=g(2,iat)+term

      dz = 2.D0*z1-2.D0*z2

      term = -3.D0*damp6*damp6*s6*C6/r8*t6*alp6*dz
     &       -1.D0*damp6*s6*gC6(3)/r6
     &       +3.D0*damp6*s6*C6/r8*dz
     &       -9.D0*damp8*damp8*s8*C6*R42/r10*t8*alp8*dz
     &       -3.D0*damp8*s8*gC6(3)*R42/r8
     &       +12.D0*damp8*s8*C6*R42/r10*dz
c    &       -11.025D0*gC6(3)*R42**2*damp10*s10/r10
c    &       -33.075D0*C6*R42**2*damp10*damp10*s10/r12*t10*alp10*dz
c    &       +55.125D0*C6*R42**2*damp10*s10/r12*dz
      g(3,iat)=g(3,iat)+term

      term = -1.D0/(1.D0+6.D0*t6)*s6*C6/r6
     &       -3.D0/(1.D0+6.D0*t8)*s8*C6*R42/r8
c    &       -11.025D0*C6*R42**2/(1.D0+6.D0*t10)*s10/r10
      if(iat.lt.jat)then
         disp=disp+term
      endif

         enddo
       
         do jat=2,n
         if(iat.eq.jat) cycle
         x1=xyz(1,jat)
         y1=xyz(2,jat)
         z1=xyz(3,jat)
CTHR IJ MOST IMPORTANT
         if( (xyz(1,iat)-x1)**2+
     .       (xyz(2,iat)-y1)**2+
     .       (xyz(3,iat)-z1)**2 .gt. rthr2*
     .       r0ab(iz(iat),iz(jat))) cycle


            do kat=1,jat-1
            if(iat.eq.kat) cycle
            x2=xyz(1,kat)
            y2=xyz(2,kat)
            z2=xyz(3,kat)
            dx = (x1-x2)**2
            dy = (y1-y2)**2
            dz = (z1-z2)**2
            r2 = dx+dy+dz

CTHR IK
            if( (xyz(1,kat)-xyz(1,iat))**2+
     .          (xyz(2,kat)-xyz(2,iat))**2+
     .          (xyz(3,kat)-xyz(3,iat))**2 
     .           .gt.rthr2*r0ab(iz(iat),iz(kat))) cycle
CTHR JK 
            if(r2.gt.rthr3*r0ab(iz(jat),iz(kat))) cycle


            R0=r0ab(iz(kat),iz(jat))
            R42=r2r4(iz(jat))*r2r4(iz(kat))
c analytically
            call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              jat,kat,iat,gC6)

      r = dsqrt(r2)
      rr=R0/r
      t6 = (rr*rs6)**alp6
      damp6 =1.d0/( 1.d0+6.d0*t6 )
      t8 = (rr*rs8)**alp8
      damp8 =1.d0/( 1.d0+6.d0*t8 )
c     t10 = (rr*rs10)**alp10
c     damp10=1.d0/( 1.d0+6.d0*t10 )

      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
c     r10 = r4*r6
c     r12 = r6**2

      term = -1.D0*damp6*s6*gC6(1)/r6
     &       -3.D0*damp8*s8*gC6(1)*R42/r8
c    &       -11.025D0*gC6(1)*R42**2*damp10*s10/r10
        g(1,iat)=g(1,iat)+term

      term = -1.D0*damp6*s6*gC6(2)/r6
     &       -3.D0*damp8*s8*gC6(2)*R42/r8
c    &       -11.025D0*gC6(2)*R42**2*damp10*s10/r10
        g(2,iat)=g(2,iat)+term

      term = -1.D0*damp6*s6*gC6(3)/r6
     &       -3.D0*damp8*s8*gC6(3)*R42/r8
c    &       -11.025D0*gC6(3)*R42**2*damp10*s10/r10
        g(3,iat)=g(3,iat)+term

            enddo
         enddo

      enddo
      endif

c BJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJBJ 
c Becke-Johnson finite damping 
      if (version.eq.4) then 
      a1 =rs6
      a2 =rs8
      s8 =s18

      disp=0

      do iat=1,n
         x1=xyz(1,iat)
         y1=xyz(2,iat)
         z1=xyz(3,iat)
         do jat=1,n
            if(iat.eq.jat) cycle
            x2=xyz(1,jat)
            y2=xyz(2,jat)
            z2=xyz(3,jat)
            dx = (x1-x2)**2
            dy = (y1-y2)**2
            dz = (z1-z2)**2
            r2 = dx+dy+dz
CTHR
            if(r2.gt.rthr) cycle
            call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),
     .                                    cn(iat),cn(jat),C6)
c stored as sqrt
            r42=r2r4(iz(iat))*r2r4(iz(jat))
c use BJ radius
            R0=a1*sqrt(3.0d0*r42)+a2
c analytically
            call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              iat,jat,iat,gC6)

      r = dsqrt(r2)
 
      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
      r10 = r4*r6
c     r12 = r6**2
      t6=(r6+R0**6)
      t8=(r8+R0**8)

      dx = 2.D0*x1-2.D0*x2

      term = 3.D0*s6*C6*r4*dx/t6**2
     &      -1.D0*s6*gC6(1)/t6
     &      +12.D0*C6*R42*s8*r6*dx/t8**2
     &      -3.D0*gC6(1)*R42*s8/t8

      g(1,iat)=g(1,iat)+term

      dy = 2.D0*y1-2.D0*y2       

      term = 3.D0*s6*C6*r4*dy/t6**2
     &      -1.D0*s6*gC6(2)/t6
     &      +12.D0*C6*R42*s8*r6*dy/t8**2
     &      -3.D0*gC6(2)*R42*s8/t8

      g(2,iat)=g(2,iat)+term

      dz = 2.D0*z1-2.D0*z2

      term = 3.D0*s6*C6*r4*dz/t6**2
     &      -1.D0*s6*gC6(3)/t6
     &      +12.D0*C6*R42*s8*r6*dz/t8**2
     &      -3.D0*gC6(3)*R42*s8/t8
     
      g(3,iat)=g(3,iat)+term

      term = -1.D0*s6*C6/t6
     &       -3.D0*s8*C6*R42/t8
      if(iat.lt.jat)then
         disp=disp+term
      endif
         enddo

         do jat=2,n
         if(iat.eq.jat) cycle
         x1=xyz(1,jat)
         y1=xyz(2,jat)
         z1=xyz(3,jat)
CTHR IJ MOST IMPORTANT
         if( (xyz(1,iat)-x1)**2+
     .       (xyz(2,iat)-y1)**2+
     .       (xyz(3,iat)-z1)**2 .gt. rthr2*
     .       r0ab(iz(iat),iz(jat))) cycle


            do kat=1,jat-1
            if(iat.eq.kat) cycle
            x2=xyz(1,kat)
            y2=xyz(2,kat)
            z2=xyz(3,kat)
            dx = (x1-x2)**2
            dy = (y1-y2)**2
            dz = (z1-z2)**2
            r2 = dx+dy+dz
CTHR JK 
            if(r2.gt.rthr3*r0ab(iz(iat),iz(jat))) cycle
CTHR IK
            if( (xyz(1,kat)-xyz(1,iat))**2+
     .          (xyz(2,kat)-xyz(2,iat))**2+
     .          (xyz(3,kat)-xyz(3,iat))**2 
     .           .gt.rthr2*r0ab(iz(iat),iz(kat))) cycle


            R42=r2r4(iz(jat))*r2r4(iz(kat))
c use BJ radius
            R0=a1*sqrt(3.0d0*R42)+a2
c analytically
            call anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,iz,c6ab,mxc,
     .              jat,kat,iat,gC6)
c           if(sum(abs(gC6)).lt.thr3)cycle

      r = dsqrt(r2)

      r4 = r2**2
      r6 = r2*r4
      r8 = r4**2
c     r10 = r4*r6
c     r12 = r6**2

      t6=(r6+R0**6)
      t8=(r8+R0**8)
      s8=s18
      term =-1.D0*s6*gC6(1)/t6
     &      -3.D0*gC6(1)*R42*s8/t8

      g(1,iat)=g(1,iat)+term

       term =-1.D0*s6*gC6(2)/t6
     &       -3.D0*gC6(2)*R42*s8/t8

      g(2,iat)=g(2,iat)+term

      term =-1.D0*s6*gC6(3)/t6
     &      -3.D0*gC6(3)*R42*s8/t8
     
      g(3,iat)=g(3,iat)+term
      
            enddo
         enddo

      enddo
      endif


      endif

 999  continue
      gnorm=sum(abs(g(1:3,1:n)))
      if(IPrnt.GT.3)then
      write(6,*)
      write(6,*)'|G|=',gnorm
      endif

      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c calculate dC6/dr analytically
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine anagrdc6(max_elem,maxc,n,cn,dcn2,dcn3,
     .           iz,c6ab,mxc,iat,jat,kat,anag)
      implicit   none
      integer    n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
      real*8     cn(*),dcn2(3,n),anag(3)
      real*8     c6ab(max_elem,max_elem,maxc,maxc,3)
      real*8     term1,term2,term3,term4
      real*8     dterm2(3),dterm3(3),dcn3(3,n,n)
      real*8     zaehler,nenner,dzaehler(3),dnenner(3)
      integer    i,j,k
      real*8 k1,k2,k3
      PARAMETER (k1=16.0d0,k2=4.0d0/3.0d0,k3=-4.0d0)

      if (iat.eq.kat) then
        dterm2=dcn2(:,iat)
        dterm3=dcn3(:,iat,jat)
      else
        dterm2=dcn3(:,kat,iat)
        dterm3=dcn3(:,kat,jat)
      endif
      zaehler=0.0d0
      nenner=0.0d0
      dzaehler=0.0d0
      dnenner=0.0d0
      do i=1,mxc(iz(iat))
        do j=1,mxc(iz(jat))
          term3=c6ab(iz(iat),iz(jat),i,j,3)-cn(jat)
          term2=c6ab(iz(iat),iz(jat),i,j,2)-cn(iat)
          term1=exp(k3*(term2*term2+term3*term3))
          zaehler=zaehler+c6ab(iz(iat),iz(jat),i,j,1)*term1
          nenner=nenner+term1
          term4=term1*k3*2.0d0
          do k=1,3
            dzaehler(k)=dzaehler(k)+c6ab(iz(iat),iz(jat),i,j,1)*term4
     .                  *(term2*dterm2(k)+term3*dterm3(k))
             dnenner(k)=dnenner(k)+term4    
     .                  *(term2*dterm2(k)+term3*dterm3(k))
          enddo
        enddo
      enddo
      if (nenner.gt.1.0d-99) then
        term4=1.0d0/(nenner*nenner)
        do k=1,3
           anag(k)=(dzaehler(k)*nenner-dnenner(k)*zaehler)*term4
        enddo
      else
        anag=0.0d0
      endif
      end


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c ncoord derivative
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine ncoorda(natoms,rcov,iz,xyz,cn,dcn2,dcn3)
      implicit none
      integer iz(*),natoms,i,max_elem
      real*8 xyz(3,*),cn(*),dcn2(3,natoms),rcov(94)
      real*8 dcn3(3,natoms,natoms)
      integer iat
      real*8 dx,dy,dz,r,damp,xn,rr,rrr,rco,tmp1,tmp2,tmp3
      real*8 k1,k2,k3
      PARAMETER (k1=16.0d0,k2=4.0d0/3.0d0,k3=-4.0d0)

      dcn2=0.0d0
      dcn3=0.0d0
      do i=1,natoms
      xn=0.0d0
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            rco=rcov(iz(i))+rcov(iz(iat))
            rr=rco/r
            rrr=1.0d0/(r*r*r)

            tmp1=exp(-k1*(rr-1.0d0))
            tmp2=1.0d0/(tmp1+1.0d0)
            tmp3=tmp1*tmp2*tmp2*k1*rco*rrr

            xn=xn+tmp2
            dcn3(1,iat,i)=tmp3*dx
            dcn3(2,iat,i)=tmp3*dy
            dcn3(3,iat,i)=tmp3*dz
            dcn2(1,i)=dcn2(1,i)+tmp3*dx
            dcn2(2,i)=dcn2(2,i)+tmp3*dy
            dcn2(3,i)=dcn2(3,i)+tmp3*dz
         endif
      enddo
      cn(i)=xn
      enddo
      dcn2=-1.0d0*dcn2
      end

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C gradient of C6(iat,jat) wrt to xyz of kat or iat                          
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

      subroutine grdc6iji(max_elem,maxc,n,dum,cn,
     .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
      implicit none  
      integer n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
      real*8  cn(*),dum(3,*),g(3),rcov(max_elem)
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)

      real*8 xi,xj,c6r,c6l,st,xxx(2),yi,yj,x0i,x0j
      integer j,jjj(2)

      jjj(1)=iat
      jjj(2)=jat

      st=1.d-5       
      do j=1,3
        dum(j,kat)=dum(j,kat)+st

        call ncoord12(n,rcov,iz,dum,jjj,xxx)
        xi=xxx(1)
        xj=xxx(2)
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6r)

        dum(j,kat)=dum(j,kat)-st*2.0d0

        call ncoord12(n,rcov,iz,dum,jjj,xxx)
        xi=xxx(1)
        xj=xxx(2)
        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6l)
        g(j)=0.5d0*(c6r-c6l)/st

        dum(j,kat)=dum(j,kat)+st
      enddo
      
      end
      

      subroutine grdc6ijk(max_elem,maxc,n,dum,cn,
     .           rcov,iz,c6ab,mxc,iat,jat,kat,g)
      implicit none  
      integer n,iz(*),max_elem,maxc,iat,jat,kat,mxc(max_elem)
      real*8  cn(*),dum(3,*),g(3),rcov(max_elem)
      real*8  c6ab(max_elem,max_elem,maxc,maxc,3)

      real*8 xi,xj,c6r,c6l,st,xxx(2),yi,yj,x0i,x0j
      integer j

c x0 is the contribution of kat to the CN of iat and jat
      call ncoord11(iat,kat,rcov,iz,dum,x0i)
      call ncoord11(jat,kat,rcov,iz,dum,x0j)

      st=1.d-5       
      do j=1,3
        dum(j,kat)=dum(j,kat)+st

        call ncoord11(iat,kat,rcov,iz,dum,yi)
        call ncoord11(jat,kat,rcov,iz,dum,yj)
        xi=cn(iat)-x0i+yi
        xj=cn(jat)-x0j+yj

        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6r)

        dum(j,kat)=dum(j,kat)-st*2.0d0

        call ncoord11(iat,kat,rcov,iz,dum,yi)
        call ncoord11(jat,kat,rcov,iz,dum,yj)
        xi=cn(iat)-x0i+yi
        xj=cn(jat)-x0j+yj

        call getc6(maxc,max_elem,c6ab,mxc,iz(iat),iz(jat),xi,xj,c6l)
        g(j)=0.5d0*(c6r-c6l)/st

        dum(j,kat)=dum(j,kat)+st
      enddo
      
      end

c special grad routine
      subroutine ncoord12(natoms,rcov,iz,xyz,ija,cn)
      implicit none  
      integer iz(*),natoms,i,ija(2)
      real*8  xyz(3,*),cn(2)
      real*8  rcov(94)

      integer iat,jjj
      real*8 dx,dy,dz,r,damp,rr,rco

      real*8 k1,k2,k3
      PARAMETER (k1=16.0d0,k2=4.0d0/3.0d0,k3=-4.0d0)

      cn=0.0d0
      do jjj=1,2 
      i=ija(jjj)
      do iat=1,natoms
         if(iat.ne.i)then
            dx=xyz(1,iat)-xyz(1,i)
            dy=xyz(2,iat)-xyz(2,i)
            dz=xyz(3,iat)-xyz(3,i)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            rr=(rcov(iz(i))+rcov(iz(iat)))/r
            damp=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))
            cn(jjj)=cn(jjj)+damp
         endif
      enddo
      enddo

      end

c special grad routine
      subroutine ncoord11(i,kat,rcov,iz,xyz,cn)
      implicit none  
      integer iz(*),i,kat
      real*8  xyz(3,*),cn
      real*8  rcov(94)

      real*8 dx,dy,dz,r,damp,rr

      real*8 k1,k2,k3
      PARAMETER (k1=16.0d0,k2=4.0d0/3.0d0,k3=-4.0d0)

c contribution of kat to CN of i
            dx=xyz(1,kat)-xyz(1,i)
            dy=xyz(2,kat)-xyz(2,i)
            dz=xyz(3,kat)-xyz(3,i)
            r=sqrt(dx*dx+dy*dy+dz*dz)
            rr=(rcov(iz(i))+rcov(iz(kat)))/r
            cn=1.d0/(1.d0+exp(-k1*(rr-1.0d0)))

      end

cc ................................................................
cc  Added routines to read/write/delete grid for each atom to disk
cc  JB 18 Oct. 2001
cc  Added new routine for symmetrizing grid  <SymGRID>
cc  JB  7 July 1997
cc ................................................................
      SUBROUTINE GRID(NAtoms, XA,     IAN,    IPRNT,  NRad,
     $                NAng,   factor, IradQ,  Dist1,  PP,
     $                IND,    DISTN,  RDIST,  AIJ,    IAtom,
     $                XXA,    WTA,    NPoint, XGRID,  WGRID)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates DFT Numerical Grid and weights for a given atom
C
C  See especially
C  "A multicenter numerical integration scheme for polyatomic molecules"
C   A.D.Becke   J.Chem.Phys.  88 (1988) 2547
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XA      -  nuclear coordinates
C  IAN     -  atomic numbers
C  IPRNT   -  print flag
C  NRad    -  number of radial points
C             (if zero, will be determined in this routine)
C             MAXIMUM VALUE: 200
C  NAng    -  angular quadrature order for Lebedev quadrature
C             (if zero, will be determined in this routine for
C              each radial shell)
C             MAXIMUM VALUE: 7  (order 35 - 434 grid points)
C  factor  -  determines grid quality (normally set to unity)
C             increase/decrease if you want proportionally more/less
C             grid points
C  IradQ   -  radial quadrature
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  Dist1   -  scratch array (size NAtoms) used to generate grid
C  PP      -   ditto
C  IND     -   ditto (indexing array)
C  DISTN   -  distance to nearest neighbour for each atom
C  RDIST   -  inverse internuclear distance array
C  AIJ     -  precomputed array of Becke weights
C  IAtom   -  atom for which grid is to be generated
C  XXA     -  precomputed generic angular quadrature grid
C  WTA     -  precomputed angular quadrature weights
C
C  on exit
C
C  NPoint  -  total number of grid points for this atom
C  XGRID   -  X,Y,Z coordinates of grid
C  WGRID   -  quadrature weights
C
      parameter (maxang=7,maxrad=400,kbecke=3)
C
      dimension XA(3,NAtoms),IAN(NAtoms)
      dimension Dist1(NAtoms),PP(NAtoms),RDIST(NAtoms,NAtoms),
     $          DISTN(NAtoms),AIJ(NAtoms,NAtoms),IND(NAtoms)
      dimension XGRID(3,*),WGRID(*)
      dimension radii(maxrad),radweight(maxrad)
      DIMENSION INDX(8),XXA(3,1130),WTA(1130)
c
      data INDX / 1, 15, 41, 91, 201, 395, 697, 1131/
C
C
      NPoint = 0            ! number of grid points
c
      ino = IAN(IAtom)
      RX = XA(1,IAtom)
      RY = XA(2,IAtom)
      RZ = XA(3,IAtom)
      If(NAng.GT.maxang) NAng = maxang
C
C  -----------------------------------
C    RADIAL GRID
C  -----------------------------------
C  determine number of radial grid points
C
      thenum = factor*50.0d0
c
      IF(NRad.LE.0) THEN
        IF(ino.LE.2) THEN
         num = thenum*0.6d0                     ! 30 points
        ELSE IF(3.LE.ino.AND.ino.LE.10) THEN
         num = thenum                           ! 50 points
        ELSE IF(11.LE.ino.AND.ino.LE.18) THEN
         num = thenum*1.4d0                     ! 70 points
        ELSE IF(19.LE.ino.AND.ino.LE.36) THEN
         num = thenum*1.8d0                     ! 90 points
        ELSE
         num = thenum*2.4d0                     ! 120 points
        ENDIF
      ELSE
        num = NRad
      ENDIF
c
      If(num.GT.maxrad) num=maxrad
      If (num.LT.0) then
        num = 100
        write(6,1000) IAtom,num
      EndIf
C
C  get Bragg-Slater atomic radius
C
      CALL braggslaterradius(ino,rbs)            ! rbs in au
      RMax = 20.0d0*rbs
      rm = rbs
C
C  now get the radial grid
C
      CALL RadialGRID(num,IradQ,ino,rm,radii,radweight)
cc       write(6,*) ' In <GRID>  no. radial points: ',num
cc       write(6,*)' factor: ',factor,' rbs: ',rbs
cc       write(6,*)' Rmax: ',Rmax,' NAng: ',NAng
cc       do ll=1,num
cc       write(6,*) ll,'  radius is ',radii(ll),' weight is ',
cc     $                  radweight(ll)
cc       enddo
C
C  "Prune" the grid by dividing into regions with different
C   numbers of angular grid points per radial shell
C
      If(NAng.EQ.0)
     $   CALL PruneGRID(ino,    rbs,    factor, R1,     R2,
     $                  R3,     lang)
C
C  ------------------------------------
C    ANGULAR GRID
C  ------------------------------------
C  standard Lebedev quadrature
C  loop over all radial points
C
      DO 20 k=1,num
c
        rrr = radii(k)
        If(rrr.GT.RMax) exit           ! ignore if radius too large
        rwght = radweight(k)
C
C  determine angular grid order
C  (i.e., number of grid points in this shell)
C
        IF(NAng.EQ.0) THEN
          If(rrr.LT.R1) Then
            iang = lang+2
          Else If(rrr.LT.R2.OR.rrr.GT.R3) Then
            iang = lang+4
          Else
            iang = lang+6
          EndIf
        ELSE
          iang = Nang
        ENDIF
C
C  determine actual grid points
C
        IA1 = INDX(iang)
        IA2 = INDX(iang+1)-1
c
        DO 10 l=IA1,IA2
        NPoint = NPoint+1
        XGRID(1,NPoint) = rrr*XXA(1,l) + RX
        XGRID(2,NPoint) = rrr*XXA(2,l) + RY
        XGRID(3,NPoint) = rrr*XXA(3,l) + RZ
        CALL WBeckeG(IAtom,XGRID(1,NPoint),XGRID(2,NPoint),
     $               XGRID(3,NPoint),NAtoms,XA,RDIST,AIJ,
     $               Dist1,PP,IND,weight)
cc        CALL WStratG(IAtom,XGRID(1,NPoint),XGRID(2,NPoint),
cc     $               XGRID(3,NPoint),rrr,NAtoms,XA,RDIST,AIJ,
cc     $               DISTN,Dist1,PP,IND,weight)
c
        WGRID(NPoint) = WTA(l)*rwght*weight
c
 10     CONTINUE
cc      write(6,*) ' Grid for Atomic Number: ',ino
cc      write(6,1234) k,rrr,ia2-ia1+1,NPoint
cc 1234 format(1X,' point: ',I4,' radius: ',F10.6,' no. ang. points: ',I4,
cc     $          ' grid pts: ',I6)
 20   CONTINUE
cc      write(6,*) ' Number of radial points: ',k-1
C
      RETURN
c
 1000 FORMAT('**WARNING** Atom ',I4,' has a Zero or Negative',
     $        ' Nuclear Charge',/,
     $        ' Taking 100 radial points for DFT Numerical Grid')
c
      END
*********************************************************
      SUBROUTINE RadialGRID(num,IradQ,ino,rm,radii,rweight)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Calculates radial grid and weights for a given atom
C
C  ARGUMENTS
C
C  num     -  number of radial points
C  IradQ   -  radial quadrature
C              0 - Euler-Maclaurin (Handy)          default
C              1 - Chebyshev second type (Becke mapping)
C  ino     -  nuclear charge (atomic number)
C  rm      -  Bragg-Slater radius of current atom
C  radii   -  on exit contains radii
C  rweight -  on exit contains radial weights
C
C
C  ----------------------------------------------------------------
C  A leading reference is
C  "Handbook of Mathematical Functions with Formulas, Graphs and
C   Mathematical Tables"  Eds. M.Abramowitz and I.A.Stegun
C  (National Bureau of Standards, Applied Mathematics Series 55)
C
C  See Numerical Analysis section, integration formulas of Gaussian
C  type, starting on page 887
C
C  Here we use ideas pioneered in this context by Axel Becke
C  "A multicenter numerical integration scheme for polyatomic molecules"
C   A.D.Becke  J.Chem.Phys. 88 (1988) 2547
C
C  See also
C  "Quadrature schemes for integrals of density functional theory"
C   C.W.Murray, N.C.Handy and G.J.Laming  Mol.Phys.  78 (1993) 997
C
C  "Improved radial grids for quadrature in molecular DFT calculations"
C   M.E.Mura and P.J.Knowles  J.Chem.Phys.  104 (1996) 9848
C  -------------------------------------------------------------------
C
      DIMENSION radii(*),rweight(*)
C
      PARAMETER (Half=0.5d0,One=1.0d0,Two=2.0d0,Three=3.0d0)
      PARAMETER (PI=3.14159 26535 89793d0)
      PARAMETER (FourPI=4.0d0*PI)
C
C
      IF(IradQ.EQ.0) THEN
C
C  Euler-Maclaurin (Handy)
C  -----------------------
C  Handy mapping from x =0,1 to r = 0,infinity is
C      r = rm*(x/(1-x))**2
C  where rm is the Bragg-Slater atomic radius
C
      rnum1 = Float(num+1)
      DO 10 I=1,num
        xi = Float(I)
        zz = rnum1-xi
        ri = rm*xi*xi/(zz*zz)
        wi = rm*Two*xi*rnum1/(zz**3)
c
        radii(I) = ri
        rweight(I) = wi*FourPI*ri*ri
 10   CONTINUE
cc
      ELSE IF(IradQ.EQ.1) THEN
C
C  Standard Becke quadrature
C  -------------------------
C  Becke mapping from x =-1,1 to r = 0,infinity is
C      r = rm*(1+x)/(1-x)
C  where rm is half the Bragg-Slater atomic radius (except H)
C  Combined with Chebyshev quadrature of the second kind
C  (see Abramowitz p.889 Type 25.4.40)
C
      If(ino.GT.2) rm = Half*rm              ! see Becke  Eq.25
c
      rnum1 = Float(num+1)
      DO 20 I=1,num
        arg = (rnum1-Float(I))*PI/rnum1
        xi = cos(arg)
        wi = sin(arg)
        wi = wi*wi*PI/rnum1
c
        ri = rm*(One+xi)/(One-xi)
        radii(I) = ri
C
C  calculate radial weight
C  ----------------------------------------------------------------
C  This includes the following terms:
C    1. standard Chebyshev weight wi calculated above
C    2. derivative w.r.t.x of the mapping function
C    3. inverse of the weighting function in the integral
C       (which in this case is sqrt(1-x*x))
C  ** NOTE ** The standard r*r radial factor AND the angular
C             factor of 4*PI are also included here
C  ----------------------------------------------------------------
C
        wi = wi*rm*Two*(ri/(One-xi))**2
        rweight(I) = wi*(One/SQRT(One-xi*xi))*FourPI
 20   CONTINUE
cc
      ELSE IF(IRadQ.EQ.2) THEN
C
C  the Log3 grid of Mura and Knowles
C  ** WARNING **  I CLEARLY DON'T KNOW WHAT I'M DOING HERE
C
      alpha = 5.0d0
cc      If(ino.EQ.3.OR.ino.EQ.4.OR.ino.EQ.11.OR.ino.EQ.12.OR.
cc     $   ino.EQ.19.OR.ino.EQ.20) alpha = 7.0d0
c
      rnum1 = Float(num+1)
      DO 25 I=1,num
        xi = Float(I)/rnum1
        arg = (One-xi*xi*xi)
        RLog3 = LOG(arg)
        ri = -alpha*RLog3
        wi = alpha*Three*xi*xi/arg
c
        radii(I) = ri
        rweight(I) = wi*FourPI*ri*ri
 25   CONTINUE
cc
      ENDIF
C
      RETURN
      END
*************************************************************
      subroutine braggslaterradius(ino,r)
      implicit real*8 (a-h,o-z)
c
c  sets Bragg-Slater atomic radii for atomic number ino
c  (see J.C.Slater,  J.Chem.Phys. 41 (1964) 3199)
c  ** WARNING **  Some values are estimates
c
      dimension radius(96)
      parameter (ANTOAU=1.88972599d0)     ! not needed accurately
cc
      data (radius(i),i=1,36)
c              H       He      Li      Be      B       C
     $     / 0.25d0, 0.20d0, 1.45d0, 1.05d0, 0.85d0, 0.70d0,
c              N       O       F       Ne      Na      Mg
     $       0.65d0, 0.60d0, 0.50d0, 0.45d0, 1.80d0, 1.50d0,
c              Al      Si      P       S       Cl      Ar
     $       1.25d0, 1.10d0, 1.00d0, 1.00d0, 1.00d0, 0.95d0,
c              K       Ca      Sc      Ti      V       Cr
     $       2.20d0, 1.80d0, 1.60d0, 1.40d0, 1.35d0, 1.40d0,
c              Mn      Fe      Co      Ni      Cu      Zn
     $       1.40d0, 1.40d0, 1.35d0, 1.35d0, 1.35d0, 1.35d0,
c              Ga      Ge      As      Se      Br      Kr
     $       1.30d0, 1.25d0, 1.15d0, 1.15d0, 1.15d0, 1.10d0/
      data (radius(i),i=37,72)
c              Rb      Sr      Y       Zr      Nb      Mo
     $     / 2.35d0, 2.00d0, 1.80d0, 1.55d0, 1.45d0, 1.45d0,
c              Tc      Ru      Rh      Pd      Ag      Cd
     $       1.35d0, 1.30d0, 1.35d0, 1.40d0, 1.60d0, 1.55d0,
c              In      Sn      Sb      Te      I       Xe
     $       1.55d0, 1.45d0, 1.45d0, 1.40d0, 1.40d0, 1.35d0,
c              Cs      Ba      La      Ce      Pr      Nd
     $       2.60d0, 2.15d0, 1.95d0, 1.85d0, 1.85d0, 1.85d0,
c              Pm      Sm      Eu      Gd      Tb      Dy
     $       1.85d0, 1.85d0, 1.85d0, 1.80d0, 1.75d0, 1.75d0,
c              Ho      Er      Tm      Yb      Lu      Hf
     $       1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.55d0/
      data (radius(i),i=73,96)
c              Ta      W       Re      Os      Ir      Pt
     $     / 1.45d0, 1.35d0, 1.35d0, 1.30d0, 1.35d0, 1.35d0,
c              Au      Hg      Tl      Pb      Bi      Po
     $       1.35d0, 1.50d0, 1.90d0, 1.60d0, 1.60d0, 1.90d0,
c              At      Rn      Fr      Ra      Ac      Th
     $       1.90d0, 1.90d0, 2.50d0, 2.15d0, 1.95d0, 1.80d0,
c              Pa      U       Np      Pu      Am      Cm
     $       1.80d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0, 1.75d0/
C
C
      If(ino.GT.96) Then
       r = 1.75d0*ANTOAU
      Else
       r = radius(ino)*ANTOAU
      EndIf
C
      RETURN
      END
***********************************************************
      SUBROUTINE PruneGRID(ino,    rm,     factor, R1,     R2,
     $                     R3,     iang)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  divides the radial grid (r) into "regions" for which different
C  numbers of angular grid points needed for accurate quadrature
C
C  ARGUMENTS
C
C  ino     -  nuclear charge (atomic number)
C  rm      -  Bragg-Slater atomic radius
C  factor  -  determines grid quality
C             (degree of grid pruning)
C
C  regions are (on exit)
C
C              r < R1             low order quadrature
C             R1 < r < R2         intermediate order quadrature
C             R2 < r < R3         high order quadrature
C              r > R3             intermediate order quadrature
C
C  iang    -  modification to "default" order of angular quadrature
C
      PARAMETER (Standard=1.0d0)
c
      IF(factor.GT.Standard) THEN
C
C  ** BEST **
C   50(26)   points for shells near nucleus
C  194(110)  points for shells at intermediate distance
C  434(302)  points all other shells
C  (NOTE:  Always 1 order lower for Hydrogen)
C
        If(ino.LT.3) Then
          R1 = 0.25d0*rm
          R2 = 0.75d0*rm
          R3 = 100.0d0*rm
          iang = 0
        Else
          R1 = 0.175d0*rm
          R2 = 0.625d0*rm
          R3 = 100.0d0*rm
          iang = +1
        EndIf
cc
      ELSE IF(factor.LT.Standard) THEN
C
C  super pruning
C  for early iterations
C  (a) first and second row
C   26(14)   points for shells near nucleus
C   110(50)  points for all other shells
C  (NOTE:  Always 1 order lower for Hydrogen)
C  (b) third row and higher
C     14     points for shells near nucleus
C     50     points for shells at intermediate distance
C    194     points all other shells
C
        If(ino.LT.3) Then
          R1 = 0.7d0*rm
          R2 = 100.0d0*rm
          R3 = 0.7d0*rm
          iang = -1
        Else If(ino.LT.19) Then
          R1 = 0.35d0*rm
          R2 = 100.0d0*rm
          R3 = 0.35d0*rm
          iang = 0
        Else
          R1 = 0.25d0*rm
          R2 = 0.75d0*rm
          R3 = 5.0d0*rm
          iang = -1
        EndIf
cc
      ELSE
C
C  "standard" pruning
C   26(14)  points for shells near nucleus
C  110(50)  points for shells at intermediate distance and very far
C  302(194) points all other shells
C  (NOTE:  Always 1 order lower for Hydrogen)
C
        If(ino.LT.3) Then
          R1 = 0.5d0*rm
          R2 = rm
          R3 = 7.5d0*rm
          iang = -1
        Else
          R1 = 0.25d0*rm
          R2 = rm
          R3 = 7.5d0*rm
          iang = 0
        EndIf
cc
      ENDIF
C
      RETURN
      END
**************************************************************
      SUBROUTINE AngGRID(WTA,XXA)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Angular grid points and weights on the surface of a sphere
C  for various Lebedev quadratures.
C  Currently all possible grid quadratures are calculated
C  These are
C     order  5     14 grid points
C     order  7     26 grid points
C     order 11     50 grid points
C     order 17    110 grid points
C     order 23    194 grid points
C     order 29    302 grid points
C     order 35    434 grid points
C
C  ARGUMENTS (on exit only)
C
C  WTA     -  angular quadrature weight
C  XXA     -  angular points on spherical surface
C
C  ............................................................
C  references
C
C    A.H.Stroud,  Approximate Calculation of Multiple Integrals
C     (Prentice Hall, 1971)    for low order grids (to order 11)
C
C    V.I.Lebedev,  Zh.Vychisl.Mat.Fiz.  15 (1975) 48
C       note corrections in             16 (1976) 293
C
C    V.I.Lebedev,  Sibirsk.Mat.Zh.      18 (1977) 132
C  ............................................................
C
C
      DIMENSION WTA(1130),XXA(3,1130)
      dimension nk17(7),nk23(7),nk29(7),nk35(7),
     $          pa17(10),pa23(16),pa29(24),pa35(33)
c
c     110 points, order 17, refined Delley feb 89
      data nk17,pa17/  17,  1,  1,  0,  3,  1,  0,
     |  3.8282704949371616D-03,  9.7937375124875125D-03,
     |  1.8511563534473617D-01,  8.2117372831911110D-03,
     |  3.9568947305594191D-01,  9.5954713360709628D-03,
     |  6.9042104838229218D-01,  9.9428148911781033D-03,
     |  4.7836902881215020D-01,  9.6949963616630283D-03/
c
c     194 points, order 23, refined Delley feb 89
      data nk23,pa23/  23,  1,  1,  1,  4,  1,  1,
     |  1.7823404472446112D-03,  5.5733831788487380D-03,
     |  5.7169059499771019D-03,  4.4469331787174373D-01,
     |  5.5187714672736137D-03,  2.8924656275754386D-01,
     |  5.1582377118053831D-03,  6.7129734426952263D-01,
     |  5.6087040825879968D-03,  1.2993354476500669D-01,
     |  4.1067770281693941D-03,  3.4577021976112827D-01,
     |  5.0518460646148085D-03,  1.5904171053835295D-01,
     |  5.2511857244364202D-01,  5.5302489162330937D-03/
c
c     302 points, order 29, refined Delley feb 89
      data nk29,pa29/  29,  1,  1,  0,  6,  2,  2,
     |  8.5459117251281481D-04,  3.5991192850255715D-03,
     |  7.0117664160895449D-01,  3.6500458076772554D-03,
     |  6.5663294102196118D-01,  3.6048226014198817D-03,
     |  4.7290541325810046D-01,  3.5767296617433671D-03,
     |  3.5156403455701051D-01,  3.4497884243058833D-03,
     |  2.2196452362941784D-01,  3.1089531224136753D-03,
     |  9.6183085226147838D-02,  2.3521014136891644D-03,
     |  5.7189558918789607D-01,  3.6008209322164603D-03,
     |  2.6441528870606625D-01,  2.9823449631718039D-03,
     |  2.5100347517704651D-01,  5.4486773725807738D-01,
     |  3.5715405542733871D-03,  1.2335485325833274D-01,
     |  4.1277240831685310D-01,  3.3923122050061702D-03/
c
c     434 points, order 35, Delley feb 89
      data nk35,pa35 /  35,  1,  1,  1,  7,  2,  4,
     |  5.2658979682244362D-04,  2.5123174189273072D-03,
     |  2.5482199720026072D-03,  6.9093463075091106D-01,
     |  2.5304038011863550D-03,  6.4566647074242561D-01,
     |  2.5132671745975644D-03,  4.9143426377847465D-01,
     |  2.5017251684029361D-03,  3.9272597633680022D-01,
     |  2.4453734373129800D-03,  2.8612890103076384D-01,
     |  2.3026947822274158D-03,  1.7748360546091578D-01,
     |  2.0142790209185282D-03,  7.5680843671780184D-02,
     |  1.4624956215946138D-03,  2.1027252285730696D-01,
     |  1.9109512821795323D-03,  4.7159869115131592D-01,
     |  2.4174423756389808D-03,  9.9217696364292373D-02,
     |  3.3443631453434549D-01,  2.2366077604378487D-03,
     |  2.0548236964030437D-01,  4.5023303825826254D-01,
     |  2.4169300443247753D-03,  3.1042840351665415D-01,
     |  5.5501523610768072D-01,  2.4966440545530860D-03,
     |  1.0680182607580483D-01,  5.9051570489252711D-01,
     |  2.5122368545634951D-03/
C
C
      I = 0
C
C  first set simple low-order grids
C
      CALL SetF3(I,WTA,XXA)
C
C  now set the Lebedev grids
C
      CALL SetLebedev(I,nk17,pa17,WTA,XXA)
      CALL SetLebedev(I,nk23,pa23,WTA,XXA)
      CALL SetLebedev(I,nk29,pa29,WTA,XXA)
      CALL SetLebedev(I,nk35,pa35,WTA,XXA)
C
      RETURN
      END
********************************************************************
      SUBROUTINE SetF3(I,WT,RR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  This routine sets the grid points and the integration weights for
C  low-order radial grids: 14 points (integrates through L=5)
C  26 points (good through L=7) and 50 points (L=11)
C  The first loop sets the grid for the 6 face centers of the cube
C  (ix=1,3, iy=1,-1,2)
C  The next loop (ix=1,-1,,-1, iy=1,-1,-1, iz=1,-1,-1) sets the
C  8 corners of the cube. The weights are 1/15 and 3/40,
C  the ratio is 8/9)
C
      DIMENSION WT(*),RR(3,*)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0,Three=3.0d0)
      PARAMETER (d840=840.0d0,dws=725760.0d0)
      PARAMETER (W1=One/15.0d0,W2=One/21.0d0,W3=9216.0d0/dws,
     $           W4=Three/40.0d0,W5=27.0d0/d840,W6=15309.0d0/dws,
     $           W7=32.0d0/d840,W8=16384.0d0/dws,W9=14641.0d0/dws)
C
C
C  Cubic grid (14 points)
C
      DO 20 IX=1,3
      DO 10 IY=1,-1,-2
      I = I+1
      WT(I) = W1
      WT(I+14) = W2
      WT(I+40) = W3
      RR(1,I) = Zero
      RR(2,I) = Zero
      RR(3,I) = Zero
      RR(IX,I) = Float(IY)
 10   CONTINUE
 20   CONTINUE

      Fc = One/DSQRT(Three)
      DO 50 IX=1,-1,-2
      DO 40 IY=1,-1,-2
      DO 30 IZ=1,-1,-2
      I = I+1
      WT(I) = W4
      WT(I+14) = W5
      WT(I+40) = W6
      RR(1,I) = IX*Fc
      RR(2,I) = IY*Fc
      RR(3,I) = IZ*Fc
 30   CONTINUE
 40   CONTINUE
 50   CONTINUE
c
      IX = I-13
      DO 60 K=IX,I
      KK = K+14
      RR(1,KK) = RR(1,K)
      RR(2,KK) = RR(2,K)
      RR(3,KK) = RR(3,K)
 60   CONTINUE
c
      I = I+14
      Fc = One/DSQRT(Two)
      DO 90 IX=1,-1,-2
      DO 80 IY=1,-1,-2
      DO 70 IZ=1,3
      I = I+1
      WT(I) = W7
      WT(I+26) = W8
      RR(IZ,I) = IX*Fc
      J = MOD(IZ,3)+1
      RR(J,I) = IY*Fc
      J = 6-IZ-J
      RR(J,I) = Zero
 70   CONTINUE
 80   CONTINUE
 90   CONTINUE
c
      IX = I-25
      DO 100 K=IX,I
      KK = K+26
      RR(1,KK) = RR(1,K)
      RR(2,KK) = RR(2,K)
      RR(3,KK) = RR(3,K)
 100  CONTINUE
c
      I = I+26
      Fc = DSQRT(One/11.0d0)
      Fv = DSQRT(9.0d0/11.0d0)
      DO 140 IX=1,-1,-2
      DO 130 IY=1,-1,-2
      DO 120 IZ=1,-1,-2
      DO 110 J=1,3
      I = I+1
      WT(I) = W9
      RR(1,I) = Fc
      RR(2,I) = Fc
      RR(3,I) = Fc
      RR(J,I) = Fv
      RR(1,I) = RR(1,I)*IX
      RR(2,I) = RR(2,I)*IY
      RR(3,I) = RR(3,I)*IZ
 110  CONTINUE
 120  CONTINUE
 130  CONTINUE
 140  CONTINUE
C
      RETURN
      END
**************************************************************
      SUBROUTINE SetLebedev(I,NK,P,WT,RR)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Sets Lebedev angular grids
C
      DIMENSION NK(*),P(*),WT(*),RR(3,*)
C
      PARAMETER (Zero=0.0d0,One=1.0d0,Two=2.0d0,Three=3.0d0)

      ip = 0
      if(nk(2).gt.0) then
        ip = ip + 1
        do 600 ix=1,3
          do 601 iy=1,-1,-2
            i = i + 1
            WT(i) = P(ip)
            RR(1,i) = Zero
            RR(2,i) = Zero
            RR(3,i) = Zero
            RR(ix,i) = iy
 601      continue
 600    continue
      end if

      if(nk(3).gt.0) then
        ip = ip + 1
        c = One/DSQRT(Three)
        do 81 ix=1,-1,-2
          do 80 iy=1,-1,-2
            do 8 iz=1,-1,-2
              i = i + 1
              WT(i) = P(ip)
              RR(1,i) = ix*c
              RR(2,i) = iy*c
              RR(3,i) = iz*c
  8         continue
 80       continue
 81     continue
      end if

      if(nk(4).gt.0) then
        ip = ip + 1
        c = One/DSQRT(Two)
        do 121 ix=1,-1,-2
          do 120 iy=1,-1,-2
            do 12 iz=1,3
              i = i + 1
              WT(i) = P(ip)
              RR(iz,i) = ix*c
              j = mod(iz,3) + 1
              RR(j,i) = iy*c
              j = 6 - iz - j
              RR(j,i) = Zero
  12        continue
 120      continue
 121    continue
      end if

      n1 = nk(5)

      do 254 jj=1,n1
        ip = ip + 1
        uu = P(ip)
        ip = ip + 1
        vv = One - Two*uu*uu
        vv = DSQRT(vv)
        do 2442 ix=1,-1,-2
          do 2441 iy=1,-1,-2
            do 2440 iz=1,-1,-2
              do 244 j=1,3
                i = i + 1
                WT(i) = P(ip)
                RR(1,i) = uu
                RR(2,i) = uu
                RR(3,i) = uu
                RR(j,i) = vv
                RR(1,i) = RR(1,i)*ix
                RR(2,i) = RR(2,i)*iy
                RR(3,i) = RR(3,i)*iz
  244         continue
 2440       continue
 2441     continue
 2442   continue
  254 continue

      n1 = nk(6)

      do 255 jj=1,n1
        ip = ip + 1
        pp = P(ip)
        ip = ip + 1
        qq = One - pp*pp
        qq = DSQRT(qq)
        do 2452 ix=1,-1,-2
          do 2451 iy=1,-1,-2
            do 2450 ii=0,1
              do 245 j=1,3
                i = i + 1
                WT(i) = P(ip)
                j1 = mod(j+ii,3)+1
                RR(j1,i) = pp*ix
                j1 = mod(j+1-ii,3) + 1
                RR(j1,i) = qq*iy
                RR(j,i) = Zero
 245          continue
2450        continue
2451      continue
2452    continue
 255  continue

      n1 = nk(7)

      do 4810 jj=1,n1
        ip = ip + 1
        rrr = P(ip)
        ip = ip + 1
        ss = P(ip)
        ip = ip + 1
        tt = DSQRT(One - rrr*rrr - ss*ss)
        do 483 ix=1,-1,-2
          do 482 iy=1,-1,-2
            do 481 iz=1,-1,-2
              do 480 j=1,3
                do 48 ii=0,1
                  i = i + 1
                  WT(i) = P(ip)
                  RR(j,i) = rrr*ix
                  j1 = mod(j+ii,3) + 1
                  RR(j1,i) = ss*iy
                  j1 = mod(j+1-ii,3) + 1
                  RR(j1,i) = tt*iz
  48            continue
 480          continue
 481        continue
 482      continue
 483    continue
4810  continue
C
      RETURN
      END
********************************************************************
      SUBROUTINE WBeckeG(IAtom,  rX,     rY,     rZ,     NAtoms,
     $                   XNuc,   RDIST,  AIJ,    Dist1,  PP,
     $                   IND,    wbecke)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates Becke weights at the current grid point
C  see  A.D.Becke, J.Chem.Phys.  88 (1988) 2547
C
C  ARGUMENTS
C
C  IAtom   -  Atom associated with current grid point
C  rX
C  rY      -  coordinates of current grid point
C  rX
C  NAtoms  -  total number of atoms
C  XNuc    -  nuclear coordinates
C  RDIST   -  inverse interatomic distance matrix
C  AIJ     -  Becke atomic size adjustments
C             (see appendix in Becke article)
C  Dist1   -  scratch array for point-atom distances
C  PP      -  scratch array for cell functions
C  IND     -  indexing array for contributing atoms
C  wbecke  -  on exit contains Becke weight
C
C
      DIMENSION XNuc(3,NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms)
      dimension Dist1(NAtoms),PP(NAtoms),IND(NAtoms)
C
      PARAMETER (Zero=0.0d0,Quarter=0.25d0,Half=0.5d0,One=1.0d0,
     $           Three=3.0d0)
      parameter (ANTOAU=1.88972599d0)     ! not needed accurately
      PARAMETER (CutOff=One/(20.0d0*ANTOAU))
C
C
C  If only one atom, weight is one
C
      If(NAtoms.EQ.1) Then
        wbecke = One
        RETURN
      EndIf
C
C  get list of atoms within Cutoff of current atom
C  at same time calculate distance from current grid point
C
      NA = 0
      DO 10 I=1,NAtoms
      If(RDIST(I,IAtom).GT.CutOff.OR.I.EQ.IAtom) Then
       NA = NA+1
       IND(NA) = I
c
       Dist1(NA) = SQRT( (rX - XNuc(1,I))**2 + (rY - XNuc(2,I))**2
     $                       + (rZ - XNuc(3,I))**2 )
       PP(NA) = One
       If(I.EQ.IAtom) JAtom=NA       ! new value for IAtom in list
c
      EndIf
 10   CONTINUE
c
      DO 30 I=2,NA
      DistI = Dist1(I)
      II = IND(I)
      DO 20 J=1,I-1
      DistJ = Dist1(J)
      JJ = IND(J)
      fmu = (DistI-DistJ)*RDIST(II,JJ)
      fmu = fmu + AIJ(II,JJ)*(One - fmu*fmu)  ! atomic size adjustments
c
cc      fmu = Half*fmu*(Three - fmu**2)
cc      fmu = Half*fmu*(Three - fmu**2)
cc      fmu = Quarter*fmu*(Three - fmu**2)
      fmu = fmu*(Three - fmu**2)    ! this is the double of mu
      fmu = fmu*(12.0d0 - fmu**2)   ! this is 16 times of it
      fmu = 6.103515625d-5*fmu*(768.0d0 - fmu**2)
c
      PP(I) = PP(I)*(Half-fmu)
      PP(J) = PP(J)*(Half+fmu)
 20   CONTINUE
 30   CONTINUE
C
      wbecke = Zero
      DO 40 I=1,NA
      wbecke = wbecke + PP(I)
 40   CONTINUE
c
      wbecke = PP(JAtom)/wbecke
c
      RETURN
      END
*******************************************************************
      SUBROUTINE WStratG(IAtom,  rX,     rY,     rZ,     rrr,
     $                   NAtoms, XNuc,   RDIST,  AIJ,    DISTN,
     $                   Dist1,  PP,     IND,    weight)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  calculates Stratmann weights at the current grid point
C  see  R.E.Stratmann, G.E.Scuseria and M.J.Frisch
C       Chem.Phys.Letts.  257 (1996) 213
C
C  ARGUMENTS
C
C  IAtom   -  Atom associated with current grid point
C  rX
C  rY      -  coordinates of current grid point
C  rX
C  rrr     -  radius of current shell
C  NAtoms  -  total number of atoms
C  XNuc    -  nuclear coordinates
C  RDIST   -  inverse interatomic distance matrix
C  AIJ     -  Becke atomic size adjustments
C             (see appendix in Becke article)
C  DISTN   -  distance to nearest neighbour for each atom
C  Dist1   -  scratch array for point-atom distances
C  PP      -  scratch array for cell functions
C  IND     -  indexing array for contributing atoms
C  weight  -  on exit contains Stratmann weight
C
C
      DIMENSION XNuc(3,NAtoms),RDIST(NAtoms,NAtoms),AIJ(NAtoms,NAtoms)
      dimension DISTN(NAtoms),Dist1(NAtoms),PP(NAtoms),IND(NAtoms)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0,One=1.0d0)
      PARAMETER (Quarter=0.25d0,Three=3.0d0)
cc  ...........................................................
cc  A is the parameter which governs the cutoff of grid points
cc  A=1 corresponds to the full Becke scheme
cc  ...........................................................
      PARAMETER (A=0.8d0,C1=35.0d0/32.0d0,C2=21.0d0/32.0d0,
     $           C3=5.0d0/32.0d0)
      PARAMETER (ATest=Half*(One-A),A1=One/A)
      parameter (ANTOAU=1.88972599d0)     ! not needed accurately
      PARAMETER (CutOff=One/(30.0d0*ANTOAU))
C
C
C  first determine if the weight is equal to one
C  check on distance to nearest neighbour
C
      IF(rrr.LT.ATest*DISTN(IAtom)) THEN
        weight = One
        RETURN
      ENDIF
C
C  get list of atoms within Cutoff of current atom
C  at same time calculate distance from current grid point
C
      NA = 0
      DO 10 I=1,NAtoms
      If(RDIST(I,IAtom).GT.CutOff.OR.I.EQ.IAtom) Then
       NA = NA+1
       IND(NA) = I
c
       Dist1(NA) = SQRT( (rX - XNuc(1,I))**2 + (rY - XNuc(2,I))**2
     $                       + (rZ - XNuc(3,I))**2 )
       PP(NA) = One
c
      EndIf
 10   CONTINUE
c
      DO 30 I=1,NA
      DistI = Dist1(I)
      II = IND(I)
      DO 20 J=1,NA
      DistJ = Dist1(J)
      JJ = IND(J)
      fmu = (DistI-DistJ)*RDIST(II,JJ)
c
      fmu = fmu + AIJ(II,JJ)*(One - fmu*fmu)  ! atomic size adjustments
c
      If(fmu.GE.A) Then
       PP(I) = Zero
       GO TO 30
      Else If(fmu.LE.-A) Then
c -- do nothing, as equivalent to multiplication by one
      Else
c
       fmu = fmu*A1
cc       fmu = Half*fmu*(Three - fmu**2)
cc       fmu = Half*fmu*(Three - fmu**2)
cc       fmu = Quarter*fmu*(Three - fmu**2)
       fmu = fmu*(Three - fmu**2)    ! this is the double of mu
       fmu = fmu*(12.0d0 - fmu**2)   ! this is 16 times of it
       fmu = 6.103515625d-5*fmu*(768.0d0 - fmu**2)
cc  ..................................................
cc  this is the original Stratmann scheme
cc       f2 = fmu*fmu
cc       fmu = fmu*(C1-f2*(C1-f2*(C2-f2*C3)))
cc  ..................................................
c
       PP(I) = PP(I)*(Half-fmu)
      EndIf
 20   CONTINUE
 30   CONTINUE
C
      weight = Zero
      DO 40 I=1,NA
      weight = weight + PP(I)
 40   CONTINUE
c
      weight = PP(IAtom)/weight
c
      RETURN
      END
*******************************************************************
      SUBROUTINE SymGRID(NAtoms, IAtom,  X,      Y,      Z,
     $                   NSym,   NGen,   IGEN,   NEqATM, NPoint,
     $                   XGRID,  WGHT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  This routine prunes the grid on the current symmetry-unique
C  atom, eliminating grid points which are equivalent by symmetry
C  and adjusting the quadrature weights for those points exactly
C  aligned with symmetry elements
C
C  ARGUMENTS
C
C  IAtom   -  current atom
C  X       -  x-coordinate of atom
C  Y       -  y-coordinate of atom
C  Z       -  z-coordinate of atom
C  NSym    -  total number of symmetry operations
C  NGen    -  number of generators of the group
C  IGEN    -  list of group generators
C             these are highest non-abelian subgroup and are:
C               1 - reflection in YZ plane
C               2 - reflection in XZ plane
C               3 - reflection in XY plane
C               4 - C2 rotation about Z-axis
C               5 - C2 rotation about Y-axis
C               6 - C2 rotation about X-axis
C               7 - inversion through origin
C             (contains ALL symmetry operations; first NGen are generators)
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  NPoint  -  number of grid points (may be modified on exit)
C  XGRID   -  X,Y,Z grid points
C  WGHT    -  grid quadrature weights
C
C
      DIMENSION IGEN(NSym),NEqATM(NAtoms,*),XGRID(3,NPoint),
     $          WGHT(NPoint)
C
      PARAMETER (Zero=0.0d0,Half=0.5d0)
C
C
C  Loop over number of generators
C  ** WARNING - ASSUMES MOLECULE IS CORRECTLY ORIENTED **
C
      IPt = 0
c
      DO 90 IOp=1,NGen
C
C  If the symmetry operation takes atom IAtom into a new atom
C  then do nothing
C
      If(NEqATM(IAtom,IOp).NE.IAtom) GO TO 90
      ISym = IGEN(IOp)
C
C  If molecule has been oriented correctly then we can prune
C  the grid on every atom which has a zero coordinate
C
      IF(ISYM.EQ.1.AND.X.EQ.Zero) THEN
C
C  atom is in YZ plane
C  eliminate all coordinates with negative X
C
        IPt = 0
        DO 10 IPoint=1,NPoint
        XX = XGRID(1,IPoint)
        If(XX.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XX
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XX.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 10     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISYM.EQ.2.AND.Y.EQ.Zero) THEN
C
C  atom is in XZ plane
C  eliminate all coordinates with negative Y
C
        IPt = 0
        DO 20 IPoint=1,NPoint
        YY = XGRID(2,IPoint)
        If(YY.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = YY
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(YY.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 20     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.3.AND.Z.EQ.Zero) THEN
C
C  atom is in XY plane
C  eliminate all coordinates with negative Z
C
        IPt = 0
        DO 30 IPoint=1,NPoint
        ZZ = XGRID(3,IPoint)
        If(ZZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = ZZ
          WGHT(IPt) = WGHT(IPoint)
          If(ZZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 30     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.4.AND.X.EQ.Zero.AND.Y.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        IPt = 0
        DO 40 IPoint=1,NPoint
        XX = XGRID(1,IPoint)
        If(XX.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XX.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 40     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.5.AND.X.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Y-axis
C
        IPt = 0
        DO 50 IPoint=1,NPoint
        ZZ = XGRID(3,IPoint)
        If(ZZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(ZZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 50     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.6.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on X-axis
C
        IPt = 0
        DO 60 IPoint=1,NPoint
        YY = XGRID(2,IPoint)
        If(YY.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(YY.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 60     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.7.AND.X.EQ.Zero.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        IPt = 0
        DO 70 IPoint=1,NPoint
        XYZ = XGRID(1,IPoint)*XGRID(2,IPoint)*XGRID(3,IPoint)
        If(XYZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XYZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 70     CONTINUE
        NPoint = IPt
cc
      ENDIF
 90   CONTINUE
C
      If(IPt.GT.0) RETURN
C
C ------------------------------------------------------------
C  WARNING - EXPERIMENTAL CODE
C
C  No pruning has been done
C  HOWEVER, there MAY be a symmetry operation of the (Abelian) group
C  that takes an atom into itself, even though the generators do not.
C  In this case we need to prune
C
      DO 91 IOp=NGen+1,NSym
      If(NEqATM(IAtom,IOp).NE.IAtom) GO TO 91
      ISym = IGEN(IOp)
C
C  If molecule has been oriented correctly then we can prune
C  the grid on every atom which has a zero coordinate
C
      IF(ISYM.EQ.1.AND.X.EQ.Zero) THEN
C
C  atom is in YZ plane
C  eliminate all coordinates with negative X
C
        IPt = 0
        DO 11 IPoint=1,NPoint
        XX = XGRID(1,IPoint)
        If(XX.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XX
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XX.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 11     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISYM.EQ.2.AND.Y.EQ.Zero) THEN
C
C  atom is in XZ plane
C  eliminate all coordinates with negative Y
C
        IPt = 0
        DO 21 IPoint=1,NPoint
        YY = XGRID(2,IPoint)
        If(YY.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = YY
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(YY.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 21     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.3.AND.Z.EQ.Zero) THEN
C
C  atom is in XY plane
C  eliminate all coordinates with negative Z
C
        IPt = 0
        DO 31 IPoint=1,NPoint
        ZZ = XGRID(3,IPoint)
        If(ZZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = ZZ
          WGHT(IPt) = WGHT(IPoint)
          If(ZZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 31     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.4.AND.X.EQ.Zero.AND.Y.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        IPt = 0
        DO 41 IPoint=1,NPoint
        XX = XGRID(1,IPoint)
        If(XX.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XX.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 41     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.5.AND.X.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Y-axis
C
        IPt = 0
        DO 51 IPoint=1,NPoint
        ZZ = XGRID(3,IPoint)
        If(ZZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(ZZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 51     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.6.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on X-axis
C
        IPt = 0
        DO 61 IPoint=1,NPoint
        YY = XGRID(2,IPoint)
        If(YY.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(YY.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 61     CONTINUE
        NPoint = IPt
cc
      ELSE IF(ISym.EQ.7.AND.X.EQ.Zero.AND.Y.EQ.Zero.AND.Z.EQ.Zero) THEN
C
C  atom is on Z-axis
C
        IPt = 0
        DO 71 IPoint=1,NPoint
        XYZ = XGRID(1,IPoint)*XGRID(2,IPoint)*XGRID(3,IPoint)
        If(XYZ.LE.Zero) Then
          IPt = IPt+1
          XGRID(1,IPt) = XGRID(1,IPoint)
          XGRID(2,IPt) = XGRID(2,IPoint)
          XGRID(3,IPt) = XGRID(3,IPoint)
          WGHT(IPt) = WGHT(IPoint)
          If(XYZ.EQ.Zero) WGHT(IPt) = WGHT(IPt)*Half
        EndIf
 71     CONTINUE
        NPoint = IPt
cc
      ENDIF
 91   CONTINUE
C
C  ----------------------------------------------------------------
C
      RETURN
      END
*******************************************************************
      SUBROUTINE RdGRID(Filename,len,NPoint,XGRID,WGHT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads grid for current atom from disk
C
C  ARGUMENTS
C
C  Filename  -  full name of grid file
C  len       -  number of characters in Filename
C
C  on exit
C
C  NPoint    -  number of grid points
C  XGRID     -  X,Y,Z grid points
C  WGHT      -  grid quadrature weights
C
C
      DIMENSION XGRID(*),WGHT(*)
      CHARACTER*256 Filename,CHAR
C
      OPEN (UNIT=40,FILE=Filename(1:len),FORM='UNFORMATTED',
     $      STATUS='OLD',ERR=95)
c
      READ(40,ERR=95) NPoint
      CALL ReadBinaryG(40,NPoint,XGRID,WGHT)
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 95   CONTINUE
      Char = 'Error Reading Grid file '//Filename(1:len)
      Call nerror(13,'File IO routine <RdGRID>',Char,0,0)
C
      END
*******************************************************************
      subroutine ReadBinaryG(IUnit,NPoint,XGRID,WGHT)
      real*8 XGRID(3,NPoint),WGHT(NPoint)
c
      READ(IUnit) XGRID,WGHT
c
      RETURN
      END
*******************************************************************
      SUBROUTINE WrGRID(Filename,len,NPoint,XGRID,WGHT)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Writes grid for current atom to disk
C
C  ARGUMENTS
C
C  Filename  -  full name of grid file
C  len       -  number of characters in Filename
C
C  on exit
C
C  NPoint    -  number of grid points
C  XGRID     -  X,Y,Z grid points
C  WGHT      -  grid quadrature weights
C
C
      DIMENSION XGRID(3,NPoint),WGHT(NPoint)
      CHARACTER*256 Filename
C
      OPEN (UNIT=40,FILE=Filename(1:len),FORM='UNFORMATTED',
     $      STATUS='UNKNOWN',ERR=95)
c
      WRITE(40,ERR=96) NPoint
      WRITE(40,ERR=96) XGRID,WGHT
c
      CLOSE (UNIT=40,STATUS='KEEP')
      RETURN
c
 95   CONTINUE
      Call message('**WARNING from routine <WrGRID>',
     $ '  Unable to open grid file  unit =',40,0)
      RETURN
c
 96   Call message('**WARNING from routine <WrGRID>',
     $ '  Unable to write to grid file  unit =',40,0)
      CLOSE (UNIT=40,STATUS='KEEP')
C
      RETURN
      END
*******************************************************************
      SUBROUTINE TidyGRID(islave,NQ,IUNQ,NSLAVE,lsemi)
      IMPLICIT INTEGER(A-Z)
C
C  Deletes all DFT grid files for all atoms
C
C  ARGUMENTS
C
C  islave  -  ID number of current slave
C              -1 for serial mode
C  NQ      -  number of symmetry-unique atoms
C  IUNQ    -  list of symmetry-unique atoms
C  NSLAVE  -  atom/slave list for parallel use
C  lsemi   -  logical flag for semidirect DFT
C
C
      integer*4 islave
      DIMENSION IUNQ(NQ),NSLAVE(NQ)
      CHARACTER Filename*256,scrfile*256,ch3*3
      Logical there
c
c -- get root scratch file from depository
      Call getchval('scrf',scrfile)
      Call rmblan(scrfile,256,len)
      len1 = len+9
C
C  loop over all atoms
C
      DO 10 IAtm=1,NQ
      IF(islave.EQ.-1.OR.islave.EQ.NSLAVE(IAtm)) THEN
        ICntr = IUNQ(IAtm)
        Write(ch3,'(I3)') ICntr
        If(ICntr.LT.100) ch3(1:1) = '0'
        If(ICntr.LT.10)  ch3(2:2) = '0'
c
c -- grid file
        Filename = scrfile(1:len)//'.grid.'//ch3
        INQUIRE(FILE=Filename(1:len1),EXIST=there)
cc      write(6,*) ' islave:',islave,' IAtm:',iatm,' ICntr:',icntr
cc      write(6,*) ' Filename: ',Filename(1:len1),' EXIST: ',there
c
        If(there) Then
          OPEN (UNIT=40,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $          STATUS='OLD')
          CLOSE (UNIT=40,STATUS='DELETE')
cc          write(6,*) ' Just closed and deleted file'
        EndIf
c
c -- potential and density file
        If(lsemi.GT.0) Then
          Filename = scrfile(1:len)//'.pots.'//ch3
          INQUIRE(FILE=Filename(1:len1),EXIST=there)
c
          If(there) Then
            OPEN (UNIT=40,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $            STATUS='OLD')
            CLOSE (UNIT=40,STATUS='DELETE')
          EndIf
c
          Filename = scrfile(1:len)//'.dens.'//ch3
          INQUIRE(FILE=Filename(1:len1),EXIST=there)
c
          If(there) Then
            OPEN (UNIT=40,FILE=Filename(1:len1),FORM='UNFORMATTED',
     $            STATUS='OLD')
            CLOSE (UNIT=40,STATUS='DELETE')
          EndIf
        EndIf
      ENDIF
 10   CONTINUE
C
      RETURN
      END

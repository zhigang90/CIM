      subroutine nuclrea(inp,gcard)

      use memory

c  this routine reads the parameters for the nuclear cards
      implicit real*8 (a-h,o-z)
      parameter (nopt=16)
      character*4 options(nopt),style
      character*256 chopval,filnam
      logical adjst,cms,efld,gcard
      dimension ioptyp(nopt)
      dimension iopval(3,nopt),ropval(3,nopt),chopval(nopt),
     $          ifound(nopt)
c     common /big/bl(100)
      common /symm/nsym,nsy(7)
      common /fieldx/ xfield,xfgrad,field(9)
      data options/'styl', 'bohr', 'symm', 'file', 'axes', 'geop',
     $             'char', 'mult', 'prin', 'd2hs', 'noor', 'nocm',
     $             'dist', 'geom', 'nucl', 'fiel'/
c
c  style is the style(e.g. pqs,tx90,gaus,pdb,mopac,gamess)
c
c  bohr shows that the Cartesians & bond lengths are given in bohrs
c
c  symm is the threshold for symmetrization. If symm=0.0d0,
c  no nuclear symmetry will be used. Otherwise, approximate
c  symmetry will be found & the molecule exactly symmetrized
c
c  file is an external filename on which the nuclear info resides
c  if this parameter is given, the data are read from a file
c
c  if axes is given, the molecule will be transformed to principal
c  axes before symmetry is checked
c
c  if geop is specified then the bond lengths etc. are calculated
c   and printed on the long output file only
c
c  char and mult are charge and multiplicity, given as real + integer
c   e.g.  CHARGE=0.0 MULT=1
c
c  the "D2hs" option calls the old (PP's) symmetrizer first
c
c  noor(ient) determines the symmetry but does not reorient the axes
c
c  nocm does not switch to centre-of-mass orientation
c   (actually centre of geometry as all atomic masses taken as unity)
c   Usually a debug keyword only
c
c  if dist is given, all interatomic distances less than the specified
c  value will be printed on the long output
c
c  field is an external electric field
c
      data ioptyp/21,0,11,21,0,0,11,1,1,0,0,0,11,21,21,13/
c
c -- default symmetry threshold = 10**-5
      sthre=0.00001d0
c -- default new geometry read from input file
      inpf=inp
c -- default print flag
      iprnt=2
c -- default charge/multiplicity
      Charg = 0.0d0
      IMult = 1
c -- default number of molecules
      NMol = 1
c -- default maximum number of atomic centres
      ncentre = 2000
c
      call getival('iout',iout)
c
      If(gcard) Then
        call readop1(inp,    nopt,   options, ioptyp, iopval,
     $               ropval, chopval,ifound)
      Else
        call izeroit(ifound,nopt)
        ifound(1) = 1
        chopval(1) = 'read'
      EndIf
c
c -- input style for geometry
c -- default is to read old geometry from existing <control> file
      style=chopval(14)
      if(ifound(15).ne.0) style=chopval(15)
      if(style.eq.'    ') style='read'
      if(ifound(1).ne.0) style=chopval(1)
      call lowerca2(style,4)
      if (style.eq.'read') then
        istyl=0
      else if (style.eq.'pqs'.or.style.eq.'tx92') then
        istyl=1
      else if (style.eq.'tx90') then
        istyl=2
      else if (style.eq.'game') then
        istyl=3
      else if (style.eq.'xray') then
        istyl=4
      else if (style.eq.'pdb') then
        istyl=5
      else if (style.eq.'turb') then
        istyl=6
      else if (style.eq.'cadp') then
        istyl=7
      else if (style.eq.'gaus') then
        istyl=8
      else if (style.eq.'mopa') then
        istyl=9
      else if (style.eq.'car') then
        istyl=10
      else if (style.eq.'zmat') then
        istyl=11
      else if (style.eq.'hin') then
        istyl=12
      else if (style.eq.'mol') then
        istyl=13
      else if (style.eq.'pqb') then
        istyl=14
      else
        call nerror(1,'GEOMETRY module',
     $        'Unknown input style for geometry in <nuclrea>',0,0)
      end if
c
c -- units (angstroms or atomic units - irrelevant if istyl=0)
      if (ifound(2).eq.1) then
       ibohr=1
       call setival('bohr',1)
      else
       ibohr=0
      end if
c
c -- ** WARNING **  Z-matrix input MUST be in angstrom
      if(ibohr.eq.1.and.istyl.eq.11)
     $   call nerror(1,'GEOMETRY module',
     $    'Z-Matrix input MUST be in angstroms and degrees',0,0)
c
c -- symmetry threshold
      if (ifound(3).eq.1) then
       sthre=ropval(1,3)
      end if
c
c -- global print flag
      if(ifound(9).eq.1) iprnt = iopval(1,9)
c
c -- geometry to be read from external file (irrelevant if istyl=0)
      if (ifound(4).eq.1) then
       filnam=chopval(4)
       nucfile=34
       call rmblan(filnam,256,len)
       open (unit=nucfile,file=filnam(1:len),status='old',err=96)
       inpf=nucfile
      else
c
c -- parse input file to predetermine how many atomic centres there are
       call parsegeom(inp,istyl,ncentre)
      end if
c
c -- reserve (usually way too much) memory for the nuclei
       call getmem(5*ncentre,iadr)    ! note1
       call getmem(ncentre,iadr3)     ! note2
       call getmem(ncentre,imol)      ! note3
c
c -- get the geometry
      IF(istyl.NE.0) THEN
c -- read in the nuclear info from file inpf, by style istyl, int bl(iadr)
c -- ibohr=1 for bohr units, na is the returned number of nuclei
       call nuclin(inpf,   istyl,  bl(iadr),ibohr, iprnt,
     $             na,   bl(iadr3),NMol, bl(imol))
c -- move any dummy atom to END of coordinate list
       call getmem(5*na,iadr4)
       call getmem(na,iadr2)
       call getmem(na,iadr5)
       call SortGEOM(na,bl(iadr),bl(iadr3),bl(iadr4),bl(iadr2),
     $               bl(iadr5), ndum1, ndum2)
       call retmem(3)
      ELSE
c -- read in geometry,charge & multiplicity from old <control> file
c    (also read in atomic masses)
       NMol = 0
       call rdgeom(na,     ndum1,  ndum2, bl(iadr),bl(iadr3),
     $             NMol,  bl(imol),Charg, IMult)
      ENDIF
c
       ndum = ndum1+ndum2
c
c -- now reserve the correct amount of memory
      call retmem(3)
c
      if(na.le.0) call nerror(2,'GEOMETRY module',
     $       'Zero nuclei read in <nuclrea>',0,0)
c
c -- na = number of atoms
c -- inuc = starting address of the nuclear info
c ****************************************************************
c   Normally the GEOMETRY module, and this subroutine, will be
c   invoked only ONCE per job. For a geometry optimization in
c   Z-matrix coordinates the GEOMETRY module will be called EACH
c   optimization cycle. After the first call the memory pointers
c   for the nuclear data are already set and are simply reused.
c *****************************************************************
      call tstival('inuc',iexist)
      if(iexist.eq.1) then
       call getival('inuc',inuc)
       call tfer(bl(iadr),bl(inuc),5*na)
       call getival('mass',iadr2)
      else
       call setival('na',na)
       call getmem(5*na,iadr)         ! note1 (reuse)
       call setival('inuc',iadr)
       call getmem(na,iadr2)          ! note2
       call setival('mass',iadr2)
       If(iadr3.ne.iadr2) call tfer(bl(iadr3),bl(iadr2),na)
      end if
c
c -- total molecular charge
      if(ifound(7).eq.1) Charg = ropval(1,7)
c
c -- multiplicity
      if(ifound(8).eq.1) IMult = iopval(1,8)
c
c -- determine empirical formula
      call EmpForm(na-ndum,bl(iadr),iprnt)
c
      if(ifound(5).eq.1) then
c   calculate the principal axes of inertia
       call princax(na,bl(iadr),bl(iadr2))
      end if
c
      if(ifound(10).eq.1) then
        call symtrize(na,bl(iadr),iprnt,sthre,nsym,nsy)
      end if
c
c -- geometry print options
      if(ifound(6).eq.1) then
cc       call geop(na,bl(iadr))
        IPRNT = 3
      end if
c
      rthre = 0.0d0
      if(ifound(13).eq.1) then
        rthre = ropval(1,13)
        If(rthre.le.0.0d0) rthre = 4.0d0   ! default
      end if
c
c -- no axis reorientation
      adjst = .true.
      if(ifound(11).eq.1) adjst=.false.
c
c -- do not put system in "centre-of-mass"
      cms = .true.
      if(ifound(12).eq.1) cms=.false.
c
c -- external field?
      if(ifound(16).gt.0) then
        efld = .true.
        field(1)=ropval(1,16)
        field(2)=ropval(2,16)
        field(3)=ropval(3,16)
      else
        efld = .false.
      endif
c
c  determine the total charge, nuclear energy and an approximation to
c  to the total energy, as well as the molecular weight and sum formula
      call nucparam(na-ndum,bl(iadr),bl(iadr2))
      call calcnuc(na-ndum,bl(iadr),enuc)
      if(iprnt.gt.0) write(iout,1000) enuc
 1000 Format(' nuclear repulsion energy is ',f16.9,' au')
c
c  determine the number of alpha/closed-shell and beta occupied MOs
      call GetNAB(na-ndum,bl(iadr),Charg,IMult,NAlpha,NBeta)
c
c ====================================================================
c  write nuclear data to <control> file
      CALL wrgeom(na,     ndum1,  ndum2, bl(iadr),NMol,
     $         bl(imol),bl(iadr2),Charg, IMult,   NAlpha,
     $            NBeta,  IPRNT,  rthre, sthre,   adjst,
     $            cms,    efld,   field)
c ====================================================================
c
c -- also put data in Texas depository
      call setival('NAlpha',NAlpha)
      call setival('NBeta',NBeta)
      call setrval('charge',Charg)
      call setival('Multip',IMult)
      call setival('ndum1',ndum1)
      call setival('ndum2',ndum2)
      call setival('ndum',ndum)
      call setival('NMol',NMol)
c
c -- make sure field set elsewhere is NOT overwritten
c  **WARNING** There is also a common block!!!
      call tstival('field',iyes)
      If(iyes.eq.0) call setival('field',ifound(16))
      If(efld) Then
        call setrval('fieldX',field(1))
        call setrval('fieldY',field(2))
        call setrval('fieldZ',field(3))
        xfield = 1.1d0     ! what on earth is this!!??
        xfgrad = 0.0d0
      EndIf
c
      return
c
c -- error handling
 96   continue
      Write(iout,*) ' Geometry input file ',filnam(1:len),
     $              ' does not exist!'
      call nerror(20,'nuclrea','file not found',0,0)
c
      end
c ====================================================================
c
      subroutine nucparam(na,xnuc,xmass)
c  this routine calculates the nuclear energy, molecular weight, total
c  charge and an approximation of the total energy. it should be called
c  after symmetrization
      implicit real*8 (a-h,o-z)
      dimension xnuc(5,na),xmass(na)
      logical error
c
      call getrval('zero',zero)
      call getrval('one',one)
      call getival('iout',iout)
      call getrval('10-8',tenm8)
      error = .False.
c  total charge
      qtot=zero
c  molecular weight
      weig=zero
c  a crude estimate of the total energy
      eapr=zero
      do 100 i=1,na
      weig=weig+xmass(i)
      q=xnuc(1,i)
      qtot=qtot+q
      eapr=eapr-q*(q-0.6d0)*(one+0.038d0*q)
      x=xnuc(2,i)
      y=xnuc(3,i)
      z=xnuc(4,i)
      do 200 j=i+1,na
        r=sqrt((xnuc(2,j)-x)**2+(xnuc(3,j)-y)**2+(xnuc(4,j)-z)**2)
        if(r.lt.tenm8) then
         write(iout,250) i,j,r
 250     format(' nuclei i and j are too close ',2i4,f15.10)
         error = .True.
        end if
 200  continue
 100  continue
      If(error) call nerror(20,'nucparam','Error',0,0)
      call setrval('qtot',qtot)
      call setrval('weig',weig)
      call setrval('eapr',eapr)
      end
c ====================================================================
c
      subroutine calcnuc(na,xnuc,enuc)
      implicit real*8(a-h,o-z)
      parameter (eps=1.0d-10)
      dimension xnuc(5,na)
c
c  calculate nuclear repulsion energy
c
      enuc = 0.0d0
      do 100 i=1,na
      q=xnuc(1,i)
      x=xnuc(2,i)
      y=xnuc(3,i)
      z=xnuc(4,i)
      do 200 j=i+1,na
        qj=xnuc(1,j)
        r=sqrt((xnuc(2,j)-x)**2+(xnuc(3,j)-y)**2+(xnuc(4,j)-z)**2)
        if(abs(q).ge.eps.and.abs(qj).ge.eps.and.r.ge.eps) then
          enuc=enuc+q*qj/r
        end if
 200  continue
 100  continue
c
      return
      end
c ====================================================================
c
      subroutine printnuc(na,xnuc)
      implicit real*8 (a-h,o-z)
      dimension xnuc(5,*)
c  printing of the nuclei
      call getival('iout',iout)
c      call getival('icon',icond)
      call getrval('angs',xang)
      write(iout,1300)
c      write(icond,1300)
 1300 format(' nuclear coordinates in Angstrom')
      call getrval('angs',angs)
      do 1400 i=1,na
      write(iout,1500) i,xnuc(5,i),(xnuc(k,i)/angs,k=2,4)
c      write(icond,1500) i,xnuc(5,i),(xnuc(k,i)/angs,k=2,4)
 1500   format(1x,i4,5x,a8,2x,3f10.6)
 1400 continue
      end
c ====================================================================
c
      subroutine nuclin(inpf,   istyl,  xnuc,   ibohr,  iprnt,
     $                  na,     xmass,  nmol,   imol)
      implicit real*8(a-h,o-z)
c
c   this subroutine reads in nuclear info
c   parameters: inpf i the unit number of the input file
c   istyl is the style: the following are available:
c   Cartesian inputs
c   1 =  pqs (default)
c   2 =  tx90  3 = gamess  4 = xray  5 = pdb   6 = turbo  7 = cadp
c   8 = Gaussian  9 = Mopac  10 = Biosym car file  11 = Z-matrix
c   xnuc is the resulting nuclear info xnuc(5,i)
c   xnuc(1,i) contains the nuclear charge, (2-4) the x,y,z coordinates,
c   5 (as a floating-point number) the atom name
c   ibohr (input)  = 1 if the data are in  bohr units
c   iprnt (input) - print flag
c   na (output) is the number of atoms
c   xmass (output) is an array, the atomic masses
c   nmol (output) is the number of molecules
c   imol (output) is an array listing pointers to the start/end
c                 of each molecule
c
      dimension xnuc(5,*),xmass(*),imol(*)
c ...........................................................
c -- likely max # atoms  (temporary character array)
      parameter (MaxAT=500)
c ...........................................................
c
      call getrval('angs',xmult)
      if (ibohr.eq.1.or.istyl.eq.11) xmult=1.0d0
c
      IF (istyl.eq.1) THEN
        call readpqs(inpf,iprnt,xnuc,na,xmass,nmol,imol)
      ELSE IF (istyl.eq.2) THEN
        call readtx90(inpf,xnuc,na)
      ELSE IF (istyl.eq.5) THEN
        call readpdb(inpf,xnuc,na)
      ELSE IF (istyl.eq.9) THEN
        call readmopac(inpf,xnuc,na)
      ELSE IF (istyl.eq.10) THEN
        call readcarfile(inpf,xnuc,na,nmol,imol)
      ELSE IF (istyl.eq.11) THEN
        call extract_zmat(inpf,iprnt,xnuc,MaxAT,na)
      ELSE IF (istyl.eq.12) THEN
        call readhyper(inpf,xnuc,na)
      ELSE IF (istyl.eq.13) THEN
        call readmol(inpf,xnuc,na)
      ELSE IF (istyl.eq.14) THEN
        call readpqb(inpf,xnuc,na)
      ENDIF
c
      if(istyl.ne.1) call setmass(na,xnuc,xmass)
      do 700 i=1,na
      do 700 k=2,4
       xnuc(k,i)=xnuc(k,i)*xmult
 700  continue
c
      end
c ====================================================================
c
      subroutine setmass(na,xnuc,xmass)
      implicit real*8 (a-h,o-z)
      dimension xnuc(5,*),xmass(na)
      dimension amass(104)
c -- isotope-average masses      ! revised & fixed    JB Nov. 2005
      data amass/
c                       H            He           Li           Be
     $   0.0d0,       1.00794d0,   4.00260d0,   6.941d0,     9.012182d0,
c          B            C            N            O            F
     $  10.811d0,    12.01115d0,  14.00674d0,  15.9994d0,   18.998403d0,
c          Ne           Na           Mg           Al           Si
     $  20.1797d0,   22.989768d0, 24.305d0,    26.981539d0, 28.0855d0,
c          P            S            Cl           Ar           K
     $  30.973762d0, 32.066d0,    35.4527d0,   39.948d0,    39.098d0,
c          Ca           Sc           Ti           V            Cr
     $  40.078d0,    44.955910d0, 47.867d0,    50.9415d0,   51.9961d0,
c          Mn           Fe           Co           Ni           Cu
     $  54.938047d0, 55.847d0,    58.933198d0, 58.693d0,    63.546d0,
c          Zn           Ga           Ge           As           Se
     $  65.39d0,     69.723d0,    72.61d0,     74.921594d0, 78.96d0,
c          Br           Kr           Rb           Sr           Y
     $  79.904d0,    83.80d0,     85.4678d0,   87.62d0,     88.905d0,
c          Zr           Nb           Mo           Tc           Ru
     $  91.224d0,    92.906377d0, 95.94d0,     97.907215d0,101.07d0,
c          Rh           Pd           Ag           Cd           In
     $ 102.905500d0,106.42d0,    107.8682d0,  112.411d0,   114.818d0,
c          Sn           Sb           Te           I            Xe
     $ 118.71d0,    121.757d0,   127.60d0,    126.904473d0,131.29d0,
c          Cs           Ba           La           Ce           Pr
     $ 132.905d0,   137.327d0,   138.906d0,   140.116d0,   140.908d0,
c          Nd           Pm           Sm           Eu           Gd
     $ 144.24d0,    145.0d0,     150.36d0,    151.964d0,   157.25d0,
c          Tb           Dy           Ho           Er           Tm
     $ 158.925d0,   162.50d0,    164.930d0,   167.26d0,    168.934d0,
c          Yb           Lu           Hf           Ta           W
     $ 173.04d0,    174.967d0,   178.49d0,    180.947d0,   183.84d0,
c          Re           Os           Ir           Pt           Au
     $ 186.207d0,   190.23d0,    192.217d0,   195.078d0,   196.967d0,
c          Hg           Tl           Pb           Bi           Po
     $ 200.59d0,    204.383d0,   207.2d0,     208.980d0,   209.0d0,
c          At           Rn           Fr           Ra           Ac
     $ 210.0d0,     222.0d0,     223.0d0,     226.0d0,     227.0d0,
c          Th           Pa           U            Np           Pu
     $ 232.038d0,   231.036d0,   238.029d0,   237.0d0,     244.0d0,
c          Am           Cm           Bk           Cf           Es
     $ 243.0d0,     247.0d0,     247.0d0,     251.0d0,     252.0d0,
c          Fm           Md           No           Lr
     $ 257.0d0,     258.0d0,     259.0d0,     262.0d0 /
c
      do 100 i=1,na
       iz=xnuc(1,i)
       xmass(i)=amass(iz+1)
 100  continue
c
      end
c ====================================================================
c
      SUBROUTINE ReadPQS(inp,    IPRNT,  XNuc,   NAtoms,  XMass,
     $                   NMol,   IMOL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads input geometry in general default format
C    atomic symbol     X     Y     Z      atomic charge  At. mass
C
C  Any "$molecule" line found within the geometry input denotes
C  the start/end of a new molecule for cluster/surface optimizations
C
C  ARGUMENTS
C
C  inp     -  unit number of input file
C  IPRNT   -  print flag
C  XNuc    -  nuclear coordinates (Texas format)
C              XNuc(1,i) -  charge (used for point charges on dummy atoms)
C              XNuc(2,i) -  X coordinate
C              XNuc(3,i) -  Y coordinate
C              XNuc(4,i) -  Z coordinate
C              XNuc(5,i) -  atomic symbol
C  NAtoms  -  number of atoms
C  XMass   -  atomic mass
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C  IMOL    -  pointers to start/end of molecules in XNuc array
C
C
      DIMENSION XNuc(5,*),XMass(*),IMOL(*)
      Character*80 char
      Character*8 name
      Parameter (ncmds=40)
      Character*4 commands(ncmds),nextcard
      character elmt(0:99)*2,digit(10)*1
      Equivalence (xname,name)
C
      Parameter (Zero=0.0d0)
c
      data commands /'titl','file','cpu ','text','calc',
     2               'geom','basi','inte','gues','scf ',
     3               'forc','intc','freq','nbo ','pop ',
     4               'pop=','semi','opti','mass','nmr ',
     5               'lmp2','numh','rest','nucl','mp2 ',
     6               'mem=','%mem','jump','clea','stop',
     7               'mtst','dyna','anfc','corr','ffld',
     8               'scan','    ','    ','    ','    '/
c
      data digit/'0','1','2','3','4','5','6','7','8','9'/
      data elmt /'x ','h ','he','li','be','b ','c ','n ',
     1           'o ','f ','ne','na','mg','al','si',
     2           'p ','s ','cl','ar','k ',
     1           'ca','sc','ti','v ','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','y ','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','i ','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','w ','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','u ','np','pu',
     8           'am','cm','bk','cf','es'/
c
      NAtoms = 0
      NMol = 0
      IMOL(1) = 0
c
 10   CONTINUE
      read(inp,900,end=30) char
 900  format(a80)
      call leadblan2(char,80,len1)
      char(80:80)=' '
c
c -- if char is one of the input keywords (indicating the end of the
c -- geometry input) or if end of file is reached, then geometry input is
c -- assumed to be complete
c
      nextcard=char(1:4)
      call lowercas(nextcard,4)
      Do icmd=1,ncmds
      If(nextcard.eq.commands(icmd).or.
     1   nextcard.eq.'    ') go to 30
      EndDo
c
c -- trap any '$coord' or '$end'
c -- (start and end of geometry input)
c
      If(char(1:6).eq.'$coord'.OR.char(1:4).eq.'$end') GO TO 10
c
c -- if the current line begins with '$' then it is assumed to
c -- designate the beginning/end of a molecule
c
      If(char(1:1).eq.'$') Then
       NMol = NMol+1
       IMOL(NMol+1) = NAtoms
       GO TO 10
      EndIf
c
c -- assume we have a geometry line
      NAtoms = NAtoms+1
      Charge = Zero
c
c -- further assume that the atomic symbol starts in column 1
c
      Do i=1,9
c -- atomic symbol should be at most 8 characters
      If(char(i:i).eq.' ') go to 20
      EndDo
c
c -- should never get here
      Call nerror(8,'GEOMETRY module',
     $  'Atomic Symbol Too Long for input atom',NAtoms,0)
c
 20   CONTINUE
      name = char(1:i-1)
      call lowercas(name,i-1)
c
c -- now determine how many fields are in the rest of the line
      Call NumFIELD(char(i:i),81-i,NumF)
c
      IF(NumF.EQ.3) THEN
c -- assume line contains geometry only
        read(char(i:80),*) X,Y,Z
c
c -- check if the name is a number. If yes, read this number, and replace
c    it with the corresponding atomic symbol
        do i=1,10
        if(name(1:1).eq.digit(i)) then
          read(name(1:2),*,err=25) INum
          go to 27
  25      call nerror(10,'ReadPQS',
     1        'Numerical atomic symbol is incomprehensible, atom=',
     2         NAtoms,NAtoms)
  27      if(INum.lt.0) then
               call nerror(11,'ReadPQS',
     1      'Negative atomic number, atom=',NAtoms,INum)
          end if
          name(1:2)=elmt(INum)
          go to 28
        end if
        end do
c -- get charge from atomic symbol and set default mass
        Call GetAtNo(1,name(1:2),INum)
  28    Continue
        Charge = Float(INum)
        Call DefMassn(INum,AtMass)
c
      ELSE IF(NumF.EQ.4) THEN
c -- line contains user-defined charge
        read(char(i:80),*) X,Y,Z,Charge
        AtMass = Zero
c
      ELSE IF(NumF.EQ.5) THEN
c -- line contains user-defined charge and mass
        read(char(i:80),*) X,Y,Z,Charge,AtMass
c
      ELSE
c -- Error in input line
      Call nerror(9,'GEOMETRY module',
     $  'Problems reading geometry - check your input',0,0)
c
      ENDIF
c
c -- now set values
      XNuc(1,NAtoms) = Charge
      XNuc(2,NAtoms) = X
      XNuc(3,NAtoms) = Y
      XNuc(4,NAtoms) = Z
      XNuc(5,NAtoms) = xname
      XMass(NAtoms)  = AtMass
c
c -- read next line
      GO TO 10
C
 30   CONTINUE
      NMol = NMol+1
      IMOL(NMol+1) = NAtoms
      backspace inp
C
      END
c ====================================================================
c
      SUBROUTINE ReadPQB(ifil,   XNuc,   NAtoms)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads input geometry from PQB file (PQSMol build file)
C     X     Y     Z      atomic symbol
C
C  ARGUMENTS
C
C  ifil    -  unit number of MOL file
C  XNuc    -  nuclear coordinates (Texas format)
C              XNuc(1,i) -  charge (used for point charges on dummy atoms)
C              XNuc(2,i) -  X coordinate
C              XNuc(3,i) -  Y coordinate
C              XNuc(4,i) -  Z coordinate
C              XNuc(5,i) -  atomic symbol
C  NAtoms  -  number of atoms
C
C
      character*2 symb,elmt(96)
      character*80 char
      character*8 blank,symb1
      dimension xnuc(5,*)
      equivalence (xsymb,symb1)
      dimension nu(96)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1           's ','k ','v ','y ','i ','w ','u ',
     2           'he','li','be','ne','na','mg','al','si','cl','ar',
     1           'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu /0, 1, 1, 5, 6, 7, 8, 9,15,16,19,23,39,53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
      data blank/'        '/
      I=1
      factor=1.0
      char=''
c
 10   CONTINUE
      read(ifil,900,end=95) char
 900  format(a80)
      char(80:80)=' '
c
      If(char(1:6). eq.'$coord') then
         READ(char(8:80),'(f15.10)', end=95) factor
         go to 10
      end if
      If(char(1:13).eq.'$connectivity') go to 96
      If(char(1:4) .eq.'$end')          go to 96
      If(char(1:1) .eq.'$'.or.char(1:8).eq.blank) go to 10

      READ(char(1:80),1000,end=95) x,y,z,symb
      call lowerca2(symb,2)
      symb1=symb
      XNuc(5,I)=xsymb
      do j=1,96
       if(symb.eq.elmt(j)) exit
      end do
      if(j.gt.96) then
        call nerror(1,'readpqb',
     1        'unidentified atomic symbol '//symb, i,i)
      end if
      XNuc(1,I)=nu(j)
      XNuc(2,I)=x*factor
      XNuc(3,I)=y*factor
      XNuc(4,I)=z*factor
      I=I+1
C
      GO TO 10

 95   CONTINUE
      call nerror(2,'readpqb',
     $              'End of file when reading PQB file',0,0)
C     all done.
 96   CONTINUE
      NAtoms=I-1
      RETURN
c
 1000 format(f15.10,1x,f15.10,1x,f15.10,6x,a2)
      end
c ====================================================================
c
      subroutine readtx90(ifil,xnuc,na)
      implicit real*8 (a-h,o-z)
c  reads the  tx90 style nuclear input
c  ifil is the input file
c  xnuc is the nuclear info, see above, na is the number of
c  nuclei (both outputs)
      character*2 nx,nn
      character*8 name
      dimension xnuc(5,*)
      equivalence (xname,name)
      data nn/'n='/
c   return address
      i=0
 100  continue
      i=i+1
      read(ifil,200,end=1000) nx
 200  format(a2)
      call lowerca2(nx,2)
      backspace ifil
      if (nx.eq.nn) then
      read(ifil,300,end=1000) name,q,x,y,z
      call lowerca2(name,8)
 300    format(2x,a8,4f10.6)
      xnuc(5,i)=xname
      xnuc(1,i)=q
      xnuc(2,i)=x
      xnuc(3,i)=y
      xnuc(4,i)=z
      go to 100
      else
      na=i-1
      return
      end if
 1000 na=i-1
      end
c ====================================================================
c
      SUBROUTINE ReadCarFile(ifil,XNuc,NAtoms,NMol,IMOL)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads a Biosym CAR file
C
C  ARGUMENTS
C
C  ifil    -  unit number of CAR file
C  XNuc    -  nuclear coordinates (Texas format)
C              XNuc(1,i) -  charge (used for point charges on dummy atoms)
C              XNuc(2,i) -  X coordinate
C              XNuc(3,i) -  Y coordinate
C              XNuc(4,i) -  Z coordinate
C              XNuc(5,i) -  atomic symbol
C  NAtoms  -  number of atoms
C  NMol    -  number of molecules (for, e.g., cluster optimization)
C  IMOL    -  pointers to start/end of molecules in XNuc array
C
      character*2 symb,elmt(96)
      character*8 name,symb1
      dimension xnuc(5,*),IMOL(*)
      equivalence (xname,name),(xsymb,symb1)
      dimension nu(96)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1           's ','k ','v ','y ','i ','w ','u ',
     2           'he','li','be','ne','na','mg','al','si','cl','ar',
     1           'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu /0, 1, 1, 5, 6, 7, 8, 9,15,16,19,23,39,53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
c
c -- read 4 cards at the beginning
      idummy=4
      inp = igetival('inp')
      if(ifil.eq.inp) idummy=2    ! two lines removed by input processor
      do i=1,idummy
        read(ifil,'(a1)')
      end do
      symb1='        '
c
c -- start reading file
      NAtoms = 0
      NMol = 0
      IMOL(1) = 0
c
 100  continue
      read(ifil,2000,end=910) name,x,y,z,symb
 110  call lowerca2(symb,2)
      call lowerca2(name,8)
      If(name(1:3).eq.'end') go to 900
      NAtoms = NAtoms+1
      symb1=symb
      XNuc(5,NAtoms)=xsymb
      do j=1,96
       if(symb.eq.elmt(j)) exit
      end do
       XNuc(1,NAtoms)=nu(j)
       XNuc(2,NAtoms)=x
       XNuc(3,NAtoms)=y
       XNuc(4,NAtoms)=z
c
c -- read next line
       go to 100
c
 900  CONTINUE
      NMol = NMol+1
      IMOL(NMol+1) = NAtoms
      read(ifil,2000,end=910) name,x,y,z,symb
      If(name(1:3).eq.'end') RETURN
      go to 110
 910  RETURN
c
 2000 format(a5,3f15.8,21x,a2)
c
      end
c ====================================================================
c
      subroutine princax(na,xnuc,xmass)
c     calculates the principal moments of inertia and
c     transforms the molecule to this axis
c  na=number of nuclei
c  xnuc(5,*)= nuclear info: q,x,y,z,name for each atom
c  xmass(na)= masses
      implicit real*8 (a-h,o-z)
      dimension xnuc(5,na),xmass(na),bb(6),u(3,3),dii(3),cm(3)
      data zero,one/0.0d0,1.0d0/
      call getival('iout',iout)
      call zeroit(bb,6)
      call zeroit(cm,3)
      xmm=zero
      do 100 i=1,na
      xm=xmass(i)
      xmm=xmm+xm
      do 110 k=1,3
        cm(k)=cm(k)+xm*xnuc(k+1,i)
 110    continue
 100  continue
      do 120 k=1,3
      cm(k)=cm(k)/xmm
      do 130 i=1,na
        xnuc(k+1,i)=xnuc(k+1,i)-cm(k)
 130    continue
 120  continue
      do 300 i=1,na
      xm=xmass(i)
      kl=0
      do 200 k=1,3
            xk=xnuc(k+1,i)
       do 200 l=1,k
         xl=xnuc(l+1,i)
         kl=kl+1
         if(k.eq.l) then
           bb(kl)=bb(kl)+xm*(xnuc(2,i)**2+xnuc(3,i)**2+xnuc(4,i)**2)
         end if
        bb(kl)=bb(kl)-xm*xk*xl
 200     continue
 300   continue
       call matdef('inertia','s',3,3)
       call matfrom('inertia',bb,6)
c      call matprint('inertia',iout)
       call matdef ('princmom','d',3,3)
       call matdef('frame','q',3,3)
       call matdiag('inertia','princmom','frame')
       call mattobl('princmom',dii,3)
c  calculate constant (hbar*1E-2)/(4c*1E-20*amu)
       call getrval('hbar',hbar)
       call getrval('c   ',c)
       call getrval('pi  ',pi)
       call getrval('angs',angs)
       call getrval('amu ',amu)
       call getrval('four',four)
       const=hbar*1.0d14*angs**2/(four*pi*amu)
c   rot const=const/I
       do 350 k=1,3
      if(dii(k).gt.1.0d-8) dii(k)=const/dii(k)
 350   continue
       write(iout,360) dii
 360   format(' rotational constants in MHz  ',3f15.2)
       const=1.0d4/c
       call mult(dii,const,3)
       write(iout,370) dii
 370   format(' rotational constants in cm-1 ',3f15.7)
c  rearrange u so that the orientation is similar to the original one
       call mattobl('frame',u,9)
       do 400 k=1,3
       l1=idamax(4-k,u(k,k),3)+k-1
      if(l1.ne.k) then
        call dswap(3,u(1,k),1,u(1,l1),1)
      end if
      if(u(k,k).lt.zero) then
        call mult(u(1,k),-one,3)
      end if
 400    continue
       call matfrom('frame',u,9)
c      call matprint('princmom',iout)
c      call matprint('frame',iout)
c  define a nucdata matrix
       call matdef('nucdata','r',5,na)
c  define its transpose
       call matdef('nuctr','r',na,5)
c  put xnuc into nucdata
       call matfrom('nucdata',xnuc,5*na)
c  transpose nucdata into nuctr
       call matpose2('nucdata','nuctr','n')
c  define the colums 2-4 in nuctr as the Cartesians
       call matsub('cartes','nuctr',2,4)
       call matmmult('cartes','frame','cartes')
       call matpose2('nuctr','nucdata','n')
       call mattobl('nucdata',xnuc,5*na)
       call matrem('cartes')
       call matrem('nuctr')
       call matrem('nucdata')
       call matrem('frame')
       call matrem('princmom')
       call matrem('inertia')
       end
c ====================================================================
c
      subroutine geop(na,xnuc)

      use memory

      implicit real*8 (a-h,o-z)
c     common /big/bl(100)
c  this routine is an interface to geopar
c  xnuc nucl. data, na= number of nuclei
      dimension xnuc(5,na)
      parameter (maxatom=500)    ! for dynamic memory allocation
      call getmem(3*na,iadr)
      call getmem(na,iadr1)
c  geopar needs the nuclear info transposed
      call getrval('angs',ang)
      ang1=1.0d0/ang
      ia1=iadr-1
      ia2=ia1+na
      ia3=ia2+na
      do 100 i=1,na
      bl(ia1+i)=xnuc(2,i)*ang1
      bl(ia2+i)=xnuc(3,i)*ang1
      bl(ia3+i)=xnuc(4,i)*ang1
      bl(iadr1+i-1)=xnuc(1,i)
 100  continue
      call getival('iout',iout)
      call geopar(bl(iadr),bl(iadr1),na,iout,na,ang,maxatom)
      call retmem(2)
      end
c ====================================================================
c
      subroutine geopar (x,q,ind,ik,idim,ang,maxatom)
      implicit real*8 (a-h,o-z)
ckwolinski:
      dimension x(idim,3), q(idim)
c
c
c     .... idim is the corresponding dimension statement in the calling
c     .... it may also read as x(idim),y(idim),z(idim) consecutively
c     .... maximum 100 atoms... easy to change
c
c ..................................................................
c -- replaced by dynamic allocation using F90        JB nov 97
      dimension ian(maxatom)
      dimension li(maxatom*(maxatom-1)/2)
      dimension lj(maxatom*(maxatom-1)/2)
      dimension rr(maxatom*(maxatom-1)/2)
      dimension icheck(maxatom*(maxatom-1)/2)
c
c assuming no more than 10k pairs of atoms of a given type.
      dimension ncc(20*maxatom),nnn(20*maxatom),noo(20*maxatom)
      dimension nch(20*maxatom),nnh(20*maxatom),noh(20*maxatom)
      dimension nnc(20*maxatom),noc(20*maxatom),non(20*maxatom)
c ....................................................................
c
c
      dimension radius(19)
      dimension u(3), v(3), w(3), z(3)
      data radius(1),radius(2),radius(3),radius(4),radius(5),radius(6),r
     1adius(7),radius(8),radius(9),radius(10)/0.4,0.3,1.25,0.97,0.82,0.7
     25,0.69,0.66,0.62,0.4/
      data radius(11),radius(12),radius(13),radius(14),radius(15),radius
     1(16),radius(17),radius(18),radius(19)/1.7,1.6,1.2,1.05,1.0,1.0,1.0
     2,1.0,1.3/
c----------------------------------------------------------------------
c
c     .... x contains the cartesian coordinates in angstrom
c     .... ind=number of atoms, ian(i)=atomic numbers, ik=output file
c     .... uses the routines normal, arcos
c
c----------------------------------------------------------------------
      if (ind.lt.2) return
      if (ind.le.maxatom) go to 10
      write (ik,200) ind
      return
c-------------------------------------------------
   10 continue
      do 20 i=1,ind
   20 ian(i)=q(i)+0.5
c-------------------------------------------------
      write (ik,210)
      nek=ind*3
      nq=0
      nq1=0
      do 30 i=1,ind
      do 30 j=1,i
       if (i.eq.j) go to 30
       nq1=nq1+1
       icheck(nq1)=0
       ii=ian(i)
       jj=ian(j)
       if (ii.eq.0) ii=2
       if (jj.eq.0) jj=2
       if (ii.gt.19) ii=19
       if (jj.gt.19) jj=19
       r=sqrt((x(i,1)-x(j,1))**2+(x(i,2)-x(j,2))**2+(x(i,3)-x(j,3))**2
     1   )
       if (r.gt.(1.25d0*(radius(ii)+radius(jj)))) go to 30
       ni=nq*nek-ind
       nq=nq+1
       li(nq)=i
       lj(nq)=j
       rr(nq)=r
       icheck(nq)=1
   30 continue
c-------------------------------------------------
c
c     calculation of geometry parameters
c
c-------------------------------------------------
c calculate bond lengths:
c
      write (ik,220)
c-------------------------------------------------
c
c select pairs of atoms of the same kind (C-C,C-N etc.)
c
      icc=0
      inn=0
      ioo=0
c
      ich=0
      inh=0
      ioh=0
c
      inc=0
      ioc=0
      ion=0
c
      do 410 m=1,nq
      iat=li(m)
      jat=lj(m)
      iat_q=ian(iat)
      jat_q=ian(jat)
c
c     Pairs of the same atoms :
      if(iat_q.eq.jat_q) then
       if(iat_q.eq.6) then
          icc=icc+1
          if(icc.gt.10000) return
          ncc(icc)=m
       endif
       if(iat_q.eq.7) then
          inn=inn+1
          if(inn.gt.10000) return
          nnn(inn)=m
       endif
       if(iat_q.eq.8) then
          ioo=ioo+1
          if(ioo.gt.10000) return
          noo(ioo)=m
       endif
      endif
c
c     Pairs of different atoms :
c
      ijat_q=iat_q*jat_q
      ijat_d=max(iat_q,jat_q)-min(iat_q,jat_q)
c
c          C and H
      if(ijat_q.eq.6 .and. ijat_d.eq.5) then
       ich=ich+1
       if(ich.gt.10000) return
       nch(ich)=m
      endif
c          N and H
      if(ijat_q.eq.7 .and. ijat_d.eq.6) then
       inh=inh+1
       if(inh.gt.10000) return
       nnh(inh)=m
      endif
c          O and H
      if(ijat_q.eq.8 .and. ijat_d.eq.7) then
       ioh=ioh+1
       if(ioh.gt.10000) return
       noh(ioh)=m
      endif
c          N and C
      if(ijat_q.eq.42 .and. ijat_d.eq.1) then
       inc=inc+1
       if(inc.gt.10000) return
       nnc(inc)=m
      endif
c          O and C
      if(ijat_q.eq.48 .and. ijat_d.eq.2) then
       ioc=ioc+1
       if(ioc.gt.10000) return
       noc(ioc)=m
      endif
c          O and N
      if(ijat_q.eq.56 .and. ijat_d.eq.1) then
       ion=ion+1
       if(ion.gt.10000) return
       non(ion)=m
      endif
  410 continue
c
c print it out :
c
      if(icc.gt.0) then
       write(ik,*)'========= Carbon-Carbon bonds ========='
       rmax=0.d0
       rmin=1000.d0
       do 466 i6=1,icc
       m=ncc(i6)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  466    continue
       write(ik,4000) rmin,rmax
 4000 format(' bond length (Ang): min.=',f7.4,' max.=',f7.4/)
      endif
      if(inn.gt.0) then
       write(ik,*)'======= Nitrogen-Nitrogen bonds ======='
       rmax=0.d0
       rmin=1000.d0
       do 477 i7=1,inn
       m=ncc(i7)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  477    continue
       write(ik,4000) rmin,rmax
      endif
      if(ioo.gt.0) then
       write(ik,*)'========= Oxygen-Oxygen bonds ========='
       rmax=0.d0
       rmin=1000.d0
       do 488 i8=1,ioo
       m=ncc(i8)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  488    continue
       write(ik,4000) rmin,rmax
      endif
      if(ich.gt.0) then
       write(ik,*)'======== Carbon-Hydrogen bonds ========'
       rmax=0.d0
       rmin=1000.d0
       do 461 i6=1,ich
       m=nch(i6)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  461    continue
       write(ik,4000) rmin,rmax
      endif
      if(inh.gt.0) then
       write(ik,*)'======= Nitrogen-Hydrogen bonds ======='
       rmax=0.d0
       rmin=1000.d0
       do 471 i7=1,inh
       m=nnh(i7)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  471    continue
       write(ik,4000) rmin,rmax
      endif
      if(ioh.gt.0) then
       write(ik,*)'======== Oxygen-Hydrogen bonds ========'
       rmax=0.d0
       rmin=1000.d0
       do 481 i8=1,ioh
       m=noh(i8)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  481    continue
       write(ik,4000) rmin,rmax
      endif
      if(inc.gt.0) then
       write(ik,*)'======== Nitrogen-Carbon bonds ========'
       rmax=0.d0
       rmin=1000.d0
       do 476 i7=1,inc
       m=nnc(i7)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  476    continue
       write(ik,4000) rmin,rmax
      endif
      if(ioc.gt.0) then
       write(ik,*)'========= Oxygen-Carbon bonds ========='
       rmax=0.d0
       rmin=1000.d0
       do 486 i8=1,ioc
       m=noc(i8)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  486    continue
       write(ik,4000) rmin,rmax
      endif
      if(ion.gt.0) then
       write(ik,*)'======== Oxygen-Nitrogen bonds ========'
       rmax=0.d0
       rmin=1000.d0
       do 487 i8=1,ion
       m=non(i8)
       rr1=rr(m)*ang
       icheck(m)=0
       write (ik,230) li(m),lj(m),rr(m),rr1
       if(rr(m).gt.rmax) rmax=rr(m)
       if(rr(m).lt.rmin) rmin=rr(m)
  487    continue
       write(ik,4000) rmin,rmax
      endif
c-------------------------------------------------
c Remaining bond lengths:
c
      write(ik,*)'=========== Remaining bonds ==========='
      nq1=0
      do 40 m=1,nq
       if(icheck(m).gt.0) then
         nq1=nq1+1
         rr1=rr(m)*ang
         write (ik,230) li(m),lj(m),rr(m),rr1
       endif
   40 continue
      if(nq1.eq.0)then
      write(ik,*)'====== there are no other  bonds ======'
      write(ik,*)'                                       '
      endif
c-------------------------------------------------
c
c calculate bond angles :
c
      if(ind.gt.100) return
c
      write (ik,240)
      do 110 m=1,nq
      do 110 k=1,m
       if (m.eq.k) go to 110
       im=li(m)
       jm=lj(m)
       ki=li(k)
       jk=lj(k)
       j=0
       if (im.eq.ki) j=1
       if (im.eq.jk) j=2
       if (jm.eq.ki) j=3
       if (jm.eq.jk) j=4
       if (j.eq.0) go to 110
       j1=im
       j2=jm
       j3=ki
       go to (50,60,70,80), j
   50    j3=jk
   60    go to 90
   70    j3=jk
   80    j1=jm
       j2=im
   90    s=0.0d0
       do 100 l=1,3
          s=s+(x(j2,l)-x(j1,l))*(x(j3,l)-x(j1,l))
  100    continue
       s=s/(rr(m)*rr(k))
       s1=arco(s)
       s=s1/57.29577951d0
       write (ik,250) j2,j1,j3,s1,s
  110 continue
      do 190 m=1,nq
      do 190 n=1,nq
      do 190 k=1,n
       if (m.eq.n.or.m.eq.k.or.n.eq.k) go to 190
       j1=li(m)
       j2=lj(m)
       j3=0
       j5=0
       j6=0
       j4=0
       if (li(n).eq.j1) j3=lj(n)
       if (lj(n).eq.j1) j3=li(n)
       if (li(n).eq.j2) j4=lj(n)
       if (lj(n).eq.j2) j4=li(n)
       if (j3.ne.0.and.li(k).eq.j1) j5=lj(k)
       if (j3.ne.0.and.lj(k).eq.j1) j5=li(k)
c
c     .... j2 out of  j3-j5-j1 (centre) plane
c
       if (j4.ne.0.and.li(k).eq.j2) j6=lj(k)
       if (j4.ne.0.and.lj(k).eq.j2) j6=li(k)
       if (j5.ne.0.or.j6.ne.0) go to 150
c
c     .... j1 out of j4-j6-j2 (centre) planr
c
       if (li(k).eq.j1) j3=lj(k)
       if (lj(k).eq.j1) j3=li(k)
       if (li(k).eq.j2) j4=lj(k)
       if (lj(k).eq.j2) j4=li(k)
       if (j3.eq.0.or.j4.eq.0) go to 150
       if (j3.eq.j4) go to 150
c
c      j3-j1-j2-j4 torsion
c
       do 120 l=1,3
          u(l)=x(j3,l)-x(j1,l)
          v(l)=x(j2,l)-x(j1,l)
          w(l)=x(j2,l)-x(j4,l)
  120    continue
       call norb (u,v,z)
       call norb (v,w,u)
       s=0.0d0
       do 130 l=1,3
  130    s=s+z(l)*u(l)
       s1=arco(s)
       call norb (u,v,w)
       s2=0.0d0
       do 140 l=1,3
  140    s2=s2+w(l)*z(l)
       if (s2.lt.0.0d0) s1=-s1
       s=s1/57.29577951d0
       write (ik,260) j3,j1,j2,j4,s1,s
  150    if (j6.eq.0) go to 160
       j5=j6
       j3=j4
       jj=j2
       j2=j1
       j1=jj
  160    if (j5.eq.0) go to 190
       do 170 l=1,3
          u(l)=x(j5,l)-x(j1,l)
          v(l)=x(j2,l)-x(j1,l)
          w(l)=x(j3,l)-x(j1,l)
  170    continue
       call norb (u,w,z)
       s=0.0d0
       do 180 l=1,3
  180    s=v(l)*z(l)+s
       s=s/rr(m)
       s1=arco(s)
       s1=90.0d0-s1
       s=s1/57.29577951d0
       write (ik,270) j2,j3,j5,j1,s1,s
  190 continue
      return
c
  200 format (/1x,56htoo many nuclei, geometry parameters are not calcul
     1ated ,i5,/)
  210 format (/,20x,25h** GEOMETRY PARAMETERS **,/)
  220 format (/,1x,35hbond lengths in angstroms and bohrs,/)
  230 format (1x,2i5,2f14.7)
  240 format (/,1x,35hbond angles in degrees and radians ,/)
  250 format (1x,3i5,f12.4,2x,f12.6)
  260 format (5x,7htorsion,2x,4i5,f12.4,2x,f12.6)
  270 format (1x,i5,1x,6hout of,3i5,6h plane,f12.4,2x,f12.6)
c
      end
c ====================================================================
c
      subroutine norb (u,v,w)
      implicit real*8 (a-h,o-z)
      dimension u(3), v(3), w(3)
c
c     w is perpendicular to the plane uv
c
      w(1)=u(2)*v(3)-u(3)*v(2)
      w(2)=u(3)*v(1)-u(1)*v(3)
      w(3)=u(1)*v(2)-u(2)*v(1)
      r=sqrt(w(1)**2+w(2)**2+w(3)**2)
      if (r.gt.1.0d-7) r=1.0d0/r
      do 10 l=1,3
   10 w(l)=w(l)*r
      return
c
      end
c ====================================================================
c
      function arco (s)
c  arc cos of s
      implicit real*8 (a-h,o-z)
      if (abs(s).lt.1.0d-6) go to 10
      s1=1.0d0-s**2
      if (s1.lt.0.0d0) s1=0.0d0
      s1=sqrt(s1)/s
      s1=atan(s1)*57.29577951d0
      if (s.lt.0.0d0) s1=s1+180.0d0
      arco=s1
      return
   10 arco=90.0d0
      return
c
      end
c ====================================================================
c
      subroutine EmpForm(na,xnuc,iprnt)
c  this routine calculates and prints the empirical chemical
c  formula
c  Parameters:
c  Input:
c       na=number of nuclei
c       xnuc(5,na): nuclear data
c       xnuc(1,i) is the charge of the i-th nucleus
c       xnuc(2-4,i) are the x,y,z coordinates in atomic units
c       xnuc(5,i) is the name, as a Hollerith (character) type data
      implicit real*8 (a-h,o-z)
      character*1 form(70)
      character*3 nn
      character*2 elmt(100)
      dimension xnuc(5,na),iform(0:100)
      data elmt/' C',' H',' N',' O',' F','He','Li','Be',' B','Ne',
     1          'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca',
     2          'Sc','Ti',' V','Cr','Mn','Fe','Co','Ni','Cu','Zn',
     3          'Ga','Ge','As','Se','Br','Kr','Rb','Sr',' Y','Zr',
     4          'Nb','Mo','Tc','Ru','Rh','Pd','Ag','Cd','In','Sn',
     5          'Sb','Te',' I','Xe','Cs','Ba','La','Ce','Pr','Nd',
     6          'Pm','Sm','Eu','Gd','Tb','Dy','Ho','Er','Tm','Yb',
     7          'Lu','Hf','Ta',' W','Re','Os','Ir','Pt','Au','Hg',
     8          'Tl','Pb','Bi','Po','At','Rn','Fr','Ra','Ac','Th',
     9          'Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es','Fm'/
c
c -----------------------------------------------------------
c ** WARNING **
c    This routine may not work properly if the atomic charge
c    has been modified
c -----------------------------------------------------------
c
      call getival('icond',icond)
      call getival('iout',iout)
c  iform(i) is the number of the atoms with atomic number i
c  e.g. iform(6)=number of carbons
      do i=0,100
        iform(i)=0
      end do
      do i=1,na
        iq=nint(xnuc(1,i))
        if(iq.gt.0.and.iq.le.100) iform(iq)=iform(iq)+1
      end do
c rearrange iform in the order C,H,N,O,F,He,Li,Be,B,Ne
      ii=iform(1)
      iform(1)=iform(6)
      iform(6)=ii
      do i=2,5
        j=i+4
        ii=iform(i)
        iform(i)=iform(j)
        iform(j)=ii
      end do
      ipos=0
      call blankit(form,70)
      do i=1,100
        if(iform(i).gt.0) then
          write(nn,'(i3)') iform(i)
          ipos=ipos+1
          if(elmt(i)(1:1).ne.' ') then
         form(ipos)=elmt(i)(1:1)
            ipos=ipos+1
          end if
          form(ipos)=elmt(i)(2:2)
          ipos=ipos+1
          if(iform(i).gt.99) then
         form(ipos)=nn(1:1)
            ipos=ipos+1
          end if
          if(iform(i).gt.9) then
         form(ipos)=nn(2:2)
            ipos=ipos+1
          end if
          if(iform(i).gt.1) then
         form(ipos)=nn(3:3)
          else
            ipos=ipos-1
          end if
        end if
      end do
      if(iprnt.gt.0) write(iout,100) (form(k),k=1,ipos)
      if(iprnt.gt.1) write(icond,100)(form(k),k=1,ipos)
  100 format(/,' Empirical Formula:',1x,50a1)
      end
c ====================================================================
c
      subroutine extract_zmat(inp,iprnt,xnuc,ndim,na)
      implicit real*8(a-h,o-z)
c
c  reads the z-matrix from the input file, writes it to the <zmat>
c  file and gets the initial Cartesian coordinates in Z-matrix orientation
c ..................................................
c -- automatic allocation of arrays in F90
c -- **WARNING**  ndim set to maximum no. of atoms
      Character*8 ZZ(34*ndim)
      Real*8 Z(26*ndim)
c ..................................................
      dimension xnuc(5,*)
      character*256 jobname
      character*80 char
      Parameter (ncmds=40)
      Character*4 commands(ncmds)
      character*1 ch,blank
      Data commands /'titl','file','cpu ','text','calc',
     2               'geom','basi','inte','gues','scf ',
     3               'forc','intc','freq','nbo ','pop ',
     4               'pop=','semi','opti','mass','nmr ',
     5               'lmp2','numh','rest','nucl','mp2 ',
     6               'mem=','%mem','jump','clea','stop',
     7               'mtst','dyna','anfc','corr','ffld',
     8               'scan','    ','    ','    ','    '/
c
      data blank/' '/
      logical found
      parameter (IUnit=1)
c
      Common /job/jobname,lenJ
c
c
c -- does the <zmat> file already exist?
c
c      inquire(file=jobname(1:lenJ)//'.zmat',exist=found)
c
c -- read z-matrix from input file and write to <zmat> file
c
c count the number of cards read
       ncards=0
       open (unit=IUnit,file=jobname(1:lenJ)//'.zmat',
     $       form='formatted',status='unknown')
c
c
 10    CONTINUE
       read(inp,900,end=20) char
 900   Format(A80)
c
c -- if char is one of the input keywords (indicating the end of the
c -- z-matrix) or if end of file is reached, then z-matrix input is
c -- assumed to be complete
c
       call lowerca2(char,80)
       call leadblan2(char,80,len)
c  if len is zero, the line is empty - stop
       if(len.eq.0) go to 20
       do k=1,ncmds
        if(char(1:4).eq.commands(k)) go to 20
       end do
c
       ncards=ncards+1
       if(ncards.eq.1) then
         write(IUnit,'(a)') '$zmatrix'
       end if
c -- assume we have a z-matrix line
c
       If(char(1:3).eq.'var') Then
         write(IUnit,'(a)') '  '    ! blank line assumed in <zmat> file
       Else
         call rmblan(char,80,len)
c -- before writing to <zmat> file, replace all commas, equal signs
c -- and colons by blanks
        Do k=1,len
          ch = char(k:k)
          If(ch.EQ.','.OR.ch.EQ.'='.OR.ch.EQ.':') char(k:k) = blank
        EndDo
        write(IUnit,'(a)') char(1:len)
      Endif
      go to 10
c
 20    CONTINUE
      if(ncards.ne.0) then
         write(IUnit,'(a)') '$end'
      end if
      close (unit=IUnit,status='keep')
c
      backspace inp
c
c -- how many z-matrix lines?
c
      call scanzmat(NZ,IErr)
      If(IErr.NE.0) Call nerror(6,'GEOMETRY module',
     $  'Z-Matrix input style but no Z-Matrix found!',0,0)
c
c -- now allocate storage and get the Cartesian coordinates
c -- from the Z-matrix
c
c -- set up character pointers
c
      IZS = 1
      IZM = IZS + NZ
      IZC = IZM + 7*NZ
      IZV = IZC + 6*NZ
c
c -- set up real pointers
c
      IGE = 1
      IGO = IGE + 3*NZ
      IG  = IGO + 4*NZ
      IXC = IG  + 3*NZ
c
c -- get the full z-matrix and Cartesian coordinates
c
      call getzmat(NZ,     ZZ(IZS), ZZ(IZM), ZZ(IZC), Z(IGE),
     $             Z(IGO), Z(IG),   ZZ(IZV), iprnt,   NVar,
     $             na,     Z(IXC))
c
c -- put coordinates in xnuc
c
      call getatno(na,ZZ(IZS),Z(IGE))
      call putnucdat(na,ZZ(IZS),Z(IGE),Z(IXC),xnuc)
c
      return
c
c
      end
c ====================================================================
c
      SUBROUTINE GetZMAT(NZ,     ZSymb,  ZMTINP, CINTNL, GEO,
     $                   IGEO,   IG,     VARNAM, IPRNT,  NVar,
     $                   NAtoms, XZ)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads Z-matrix from <zmat> file and returns Cartesian coordinates
C  in Z-matrix orientation and number of (real) atoms
C
C  ARGUMENTS
C
C  NZ      -  number of atomic centres (including dummy atoms)
C             in Z-matrix
C  ZSymb   -  Z-matrix symbols (real & dummy atoms)
C  ZMTINP  -  intermediate scratch storage for the up to 7 possible
C             separate data types per Z-matrix line
C  CINTNL  -  character scratch array for input parameters
C  GEO     -  Z-matrix parameters (bond lengths, angles & dihedrals)
C  IGEO    -  Z-matrix connectivity
C  IG      -  array determining what to do with Z-matrix parameter
C                 0 - optimize it
C                 J - assign same value as previous (Jth) variable
C                -J - assign same value, opposite sign
C              1000 - fixed
C  VARNAM  -  names of all "variables" (including fixed)
C  IPRNT   -  print flag
C  NVar    -  number of variables
C  NAtoms  -  number of real atoms
C  XZ      -  Cartesian coordinates of real atoms
C
C
      DIMENSION GEO(NZ,3),IGEO(NZ,4),IG(3*NZ),XZ(3,NZ)
      CHARACTER*8 ZSymb(NZ),ZMTINP(NZ,7),CINTNL(3*NZ,2),VARNAM(3*NZ)
C
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree
c
      ToRAD = PI/180.0d0
      iout=igetival('iout')
C
C
C  first read the z-matrix (XZ used as scratch)
C
      CALL RdZMAT(NZ,   ZSymb,  ZMTINP, CINTNL, XZ,
     $            GEO,  IGEO,   IG,     VARNAM, NVar)
c
      If(IPRNT.GT.2) CALL PrntZMAT(iout,NZ,ZSymb,GEO,IGEO)
C
C  convert angles and dihedrals to radians and
C  distances to bohr
C
      DO 10 I=1,NZ
      GEO(I,1) = ANTOAU*GEO(I,1)
      GEO(I,2) = ToRAD*GEO(I,2)
      GEO(I,3) = ToRAD*GEO(I,3)
 10   CONTINUE
C
C  Get Cartesian coordinates
C
      CALL GMETRY(NZ,GEO,IGEO,XZ)
C
C  get atomic symbols from Z-matrix symbols
C  and assign dummy atoms
C
      CALL GetAtSym(NZ,ZSymb,ZMTINP)
c
      If(IPRNT.GT.2) Then
        WRITE(iout,1000)
        CALL PrntCAR(iout,0,NZ,ZMTINP,XZ)
      EndIf
C
C  remove dummy atoms
C
      CALL RmDummy(NZ,ZMTINP,XZ,NAtoms,ZSymb,XZ)
C
      RETURN
c
 1000 FORMAT(/,4X,' Cartesian Coordinates in Z-Matrix Orientation')
c
      END
c ====================================================================
c
      subroutine putnucdat(na,AtSymb,ian,xc,xnuc)
      implicit real*8(a-h,o-z)
c
c  interface with old texas for geometry, atomic numbers & symbols
c
      dimension ian(na),xc(3,na),xnuc(5,*)
      character*8 AtSymb(na),chrnam
      equivalence (chrnam,xname)
c
      do i=1,na
      xnuc(1,i) = dfloat(ian(i))
      xnuc(2,i) = xc(1,i)
      xnuc(3,i) = xc(2,i)
      xnuc(4,i) = xc(3,i)
c -- this ridiculous chicanery needed because of PC
      chrnam = AtSymb(i)
      xnuc(5,i) = xname
      enddo
c
      return
      end
c
c ====================================================================
c
      subroutine readpdb(ifil,xnuc,na)
      implicit real*8 (a-h,o-z)
c  reads the pdb (Protein Database) style nuclear input
c  ifil is the input file
c  xnuc is the nuclear info, see above, na is the number of
c  nuclei (both outputs)
      character*2 symb,elmt(96)
      character*8 name,symb1
      dimension xnuc(5,*)
      equivalence (xname,name),(xsymb,symb1)
      dimension nu(96)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1's ','k ','v ','y ','i ','w ','u ',
     2            'he','li','be','ne','na','mg','al','si','cl','ar',
     1            'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2            'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3            'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4            'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5            'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6            'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7            'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu /0,1,1,5, 6, 7,  8,  9, 15, 16, 19, 23, 39,  53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
c   return address
      symb1='        '
  50  continue
      read(ifil,200) name
      call lowerca2(name,8)
c      if(name.eq.'header'.or.name.eq.'remark') go to 50
c      if(name.eq.'compnd'.or.name.eq.'author') go to 50
      if(name(1:4).ne.'atom'.and.name(1:4).ne.'heta') go to 50
      backspace ifil
      i=0
 100  continue
      read(ifil,200,end=1000) name,ii,symb,x,y,z
 200  format(a6,i5,1x,a2,16x,3f8.3)
      call lowerca2(symb,2)
      call lowerca2(name,8)
      if(symb(1:1).eq.' ') then
        symb(1:1)=symb(2:2)
        symb(2:2)=' '
      end if
      if(name(1:3).eq.'end') go to 1000
      if(name(1:6).ne.'hetatm'.and.name(1:4).ne.'atom') go to 100
c       if(name(1:4).eq.'atom') symb(2:2)=' '
      i=i+1
      symb1=symb
      xnuc(5,i)=xsymb
      do j=1,96
        if(symb.eq.elmt(j)) exit
      end do
       xnuc(1,i)=nu(j)
       xnuc(2,i)=x
       xnuc(3,i)=y
       xnuc(4,i)=z
       go to 100
 1000 na=i
      end
c ====================================================================
      subroutine readhyper(ifil,xnuc,na)
      implicit real*8 (a-h,o-z)
c  reads the HYPERCHEM (TM) *.hin style nuclear input
c  ifil is the input file
c  xnuc is the nuclear info, see above, na is the number of
c  nuclei (both outputs)
      character*2 symb,elmt(96)
      character*80 name
      character*8 symb1
      dimension xnuc(5,*)
      equivalence (xsymb,symb1)
      dimension nu(96)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1           's ','k ','v ','y ','i ','w ','u ',
     2           'he','li','be','ne','na','mg','al','si','cl','ar',
     1           'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu / 0, 1, 1, 5, 6, 7, 8, 9,15,16,19,23,39,53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
c   return address
      i=0
 100  continue
      read(ifil,200,end=110) name
 200  format(a80)
 110  if(name(1:6).eq.'endmol') go to 1000
      if(name(1:5).ne.'atom') go to 100
      i=i+1
c locate fields
      ifield=1
      symb='  '
      do k=2,80
      if(name(k:k).eq.' '.and.name(k-1:k-1).ne.' ') then
        ifield=ifield+1
c   this is the counter for writing out a condensed version
        if(ifield.eq.4.or.ifield.eq.8) k1=0
      end if
      if(ifield.eq.4.and.name(k:k).ne.' '.and.k1.lt.2) then
        k1=k1+1
        symb(k1:k1)=name(k:k)
      end if
      if(ifield.eq.8.or.ifield.eq.9.or.ifield.eq.10) then
        k1=k1+1
        name(k1:k1)=name(k:k)
      end if
      end do
      k1=k1+1
      name(k1:k1)=' '
      read(name,*) x,y,z
      call lowerca2(symb,2)
      symb1=symb//'      '
      kk=k
      xnuc(5,i)=xsymb
      do j=1,96
      if(symb.eq.elmt(j)) exit
      end do
      if(j.gt.96) then
        call nerror(1,'readhyper',
     1        'unidentified atomic symbol '//symb, i,i)
      end if
      xnuc(1,i)=nu(j)
      xnuc(2,i)=x
      xnuc(3,i)=y
      xnuc(4,i)=z
      go to 100
 1000 continue
      read(ifil,'(a3)',end=1100) name(1:3)
      go to 100
 1100 na=i
      end
c ====================================================================
      SUBROUTINE GetNAB(NAtoms,XNuc,Chrge,IMult,NAlpha,NBeta)
      IMPLICIT INTEGER(A-Z)
C
C
C  From the charge and multiplicity, determine the number of
C  singly and doubly occupied MOs and subsequently the number
C  of alpha/beta occupied spin orbitals
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  XNuc    -  Cartesian coordinates      (old Texas format)
C             (also includes atomic symbols and numbers)
C  Chrge   -  charge
C  IMult   -  multiplicity
C
C  on exit
C
C  NAlpha  -  number of occupied closed-shell/alpha spin orbitals
C  NBeta   -  number of occupied beta spin orbitals
C
C
      REAL*8 XNuc(5,NAtoms),xname,Chrge,Totel,CThrsh
      Character*8 name
      Equivalence (xname,name)
c
      PARAMETER(CThrsh=1.0d-6)
C
C
C  get the total number of electrons
C
      Totel = 0.0d0
      DO 10 I=1,NAtoms
      xname = XNuc(5,I)
      If(name(1:1).NE.'x'.OR.(name(1:1).EQ.'x'.AND.name(2:2).EQ.'e'))
     $                        Totel = Totel + XNuc(1,I)
 10   CONTINUE
      Totel = Totel - Chrge
      NTotel = NINT(Totel)
c
c -- check for integer total charge
      If(Abs(Totel-DBLE(NTotel)).GT.CThrsh) Then
        Call nerror(13,'GEOMETRY module',
     $    'Requested Charge Results in non-integer number of electrons',
     $     0,0)
      EndIf
c
      NUn = IMult-1
      NOcc = (NTotel-NUn)/2
C
C  check
C
      If( (NUn.LT.0). OR. (2*NOcc+NUn.NE.NTotel) ) Then
        Call nerror(4,'GEOMETRY module',
     $     'Requested Charge/Multiplicity Impossible in this Molecule',
     $      0,0)
      EndIf
C
C  determine NAlpha/NBeta
C
      NAlpha = NOcc + NUn       ! number occupied alpha MOs
      If(NUn.GT.0) Then
       NBeta = NOcc             ! number occupied beta MOs
      Else
       NBeta = 0
      EndIf
C
      RETURN
      END
c ====================================================================
c
      SUBROUTINE ReadMOL(ifil,   XNuc,   NAtoms)
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  Reads input geometry from MOL file
C  Not sure of exact format - appears to be one (or more?)
C  blank lines, followed by a line on which the first number
C  is the number of atoms, followed by lines of the form
C     X     Y     Z      atomic symbol
C
C  ARGUMENTS
C
C  ifil    -  unit number of MOL file
C  XNuc    -  nuclear coordinates (Texas format)
C              XNuc(1,i) -  charge (used for point charges on dummy atoms)
C              XNuc(2,i) -  X coordinate
C              XNuc(3,i) -  Y coordinate
C              XNuc(4,i) -  Z coordinate
C              XNuc(5,i) -  atomic symbol
C  NAtoms  -  number of atoms
C
C
      character*2 symb,elmt(96)
      character*8 Char,blank,symb1
      dimension xnuc(5,*)
      equivalence (xsymb,symb1)
      dimension nu(96)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1           's ','k ','v ','y ','i ','w ','u ',
     2           'he','li','be','ne','na','mg','al','si','cl','ar',
     1           'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu /0, 1, 1, 5, 6, 7, 8, 9,15,16,19,23,39,53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
      data blank/'        '/
C
C
C  read blank lines, if any
C
 5    CONTINUE
      Char = blank
      READ(ifil,'(a8)',end=6) Char
 6    If(Char.EQ.blank) GO TO 5
      READ(Char,*) NAtoms
C
C  now read the geometry
C
      DO 10 I=1,NAtoms
      READ(ifil,1000,end=95) x,y,z,symb
      call lowerca2(symb,2)
      symb1=symb
      XNuc(5,I)=xsymb
      do j=1,96
       if(symb.eq.elmt(j)) exit
      end do
      if(j.gt.96) then
        call nerror(1,'readmol',
     1        'unidentified atomic symbol '//symb, i,i)
      end if
      XNuc(1,I)=nu(j)
      XNuc(2,I)=x
      XNuc(3,I)=y
      XNuc(4,I)=z
 10   CONTINUE
C
      RETURN
C
 95   CONTINUE
      call nerror(2,'readmol',
     $              'End of file when reading MOL file',0,0)
c
 1000 format(3f10.4,1x,a2)
c
      END
c ====================================================================
c
      subroutine readmopac(ifile,xnuc,na)
c  this routine reads an ascii MOPAC input file and extracts
c  the atomic charges, symbols and Cartesian coordinates from it
c  Parameters:
c  Input: ifile (file number)
c  Output	xnuc(5,na),na
c   xnuc must be declared as xnuc(5,nmax) in the calling program
c  where nmax>=na. na is the number of atoms in the MOPAC file
      implicit real*8 (a-h,o-z)
      dimension xnuc(5,*),a(3),b(3),c(3),w(3)
      character*80 line
      character*2 symb,elmt
      character*8 symb1
      equivalence (symb1,xsymb)
      parameter(mendeleev=96)
      dimension elmt(mendeleev),nu(mendeleev)
      data elmt /'x ','h ','d ','b ','c ','n ','o ','f ','p ',
     1           's ','k ','v ','y ','i ','w ','u ',
     2           'he','li','be','ne','na','mg','al','si','cl','ar',
     1           'ca','sc','ti','cr','mn','fe','co','ni','cu','zn',
     2           'ga','ge','as','se','br','kr','rb','sr','zr','nb',
     3           'mo','tc','ru','rh','pd','ag','cd','in','sn','sb',
     4           'te','xe','cs','ba','la','ce','pr','nd','pm','sm',
     5           'eu','gd','tb','dy','ho','er','tm','yb','lu','hf',
     6           'ta','re','os','ir','pt','au','hg','tl','pb','bi',
     7           'po','at','rn','fr','ra','ac','th','pa','np','pu'/
      data nu / 0, 1, 1, 5, 6, 7, 8, 9,15,16,19,23,39,53,74,92,
     &          2, 3, 4,10,11,12,13,14,17,18,
     1         20,21,22,24,25,26,27,28,29,30,
     2         31,32,33,34,35,36,37,38,40,41,
     3         42,43,44,45,46,47,48,49,50,51,
     4         52,54,55,56,57,58,59,60,61,62,
     5         63,64,65,66,67,68,69,70,71,72,
     6         73,75,76,77,78,79,80,81,82,83,
     7         84,85,86,87,88,89,90,91,93,94/
      data zero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/
c
      symb1='        '
      pi4=atan(one)
      degree=pi4/45.0d0
c  The first 3 cards are dismissed
      do i=1,3
        read(ifile,100)
        end do
      na=0
   50 continue
      read(ifile,100,end=1000) line
  100 format(a80)
c  remove leading blanks
      i1=0
      if(line(1:1).eq.' ') then
        do i=1,80
          if(line(i:i).ne.' '.or.i1.ge.1) then
            i1=i1+1
          end if
          line(i1:i1)=line(i:i)
          line(i:i)=' '
        end do
      end if
  200 symb=line(1:2)
      call lowercas(symb,2)
c  CHANGE
      do i=1,mendeleev
        if(symb.eq.elmt(i)) then
          q=nu(i)
          na=na+1
          goto 300
        end if
      end do
      goto 1000
  300 continue
      xnuc(1,na)=q
      symb1(1:2)=symb
      xnuc(5,na)=xsymb
      do i=3,80
        line(i-2:i-2)=line(i:i)
      end do
      if(na.eq.1) then
        xnuc(2,1)=zero
        xnuc(3,1)=zero
        xnuc(4,1)=zero
        go to 50
      else if(na.eq.2) then
        read(line,*) r
        xnuc(2,2)=zero
        xnuc(3,2)=zero
        xnuc(4,2)=r
        goto 50
      else if(na.eq.3) then
        read(line,*,end=400) r,ii,theta,jj,tau,kk,i,j,k
        go to 500
  400   continue
        read(line,*,err=1000) r,ii,theta,jj,i,j,k
  500   continue
        theta=theta*degree
c note that 2=x, 3=y,4=z!!
        xnuc(2,3)=r*sin(theta)
        xnuc(3,3)=zero
        if(i.eq.2.and.j.eq.1) then
          xnuc(4,3)=xnuc(4,2)-r*cos(theta)
        else if (i.eq.1.and.j.eq.2) then
          xnuc(4,3)=r*cos(theta)
        else
          call nerror(1,'readmopac',
     1       'incorrect connecting atoms for atom #3',i,j)
        end if
        go to 50
        else
          read(line,*) r,ii,theta,jj,tau,kk,i,j,k
          if(i.ge.na.or.j.ge.na.or.k.ge.na.or.i.eq.j.or.i.eq.k
     1     .or.j.eq.k) then
             call nerror(2,'readmopac',
     1        'error in connectng atom definition at ',na,k)
        end if
        do l=1,3
          a(l)=xnuc(l+1,k)
          b(l)=xnuc(l+1,j)
          c(l)=xnuc(l+1,i)
        end do
        theta=theta*degree
        tau=tau*degree
        call nextatom(a,b,c,r,theta,tau,w)
        do l=1,3
          xnuc(l+1,na)=w(l)
        end do
        go to 50
      end if
 1000 continue
      end
c ====================================================================
c
      subroutine nextatom(a,b,c,r,theta,tau,w)
c  this routine determines the position of atom w
c  located at distance r from atom c, the angle b-c=w is theta
c  and the torsional angle a-b-c-w is tau. a,b and c should not
c  be collinear
c  parameters: input a(3),b(3),c(3),r,theta,tau
c              output: w(3)
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),c(3),u(3),v(3),w(3),axis(3)
      data zero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/
      pi4=atan(1.0d0)
      call bondvec(b,a,u)
      call bondvec(b,c,v)
c  determine the rotation needed to rotate v to the +z axis
      vv=sqrt(v(1)**2+v(2)**2)
      if(vv.lt.1.0d-8) then
        if(v(3).gt.zero) then
c  v is oriented along +z already
          call azimuth(u,r,theta,tau,w)
          go to 1000
        else
          axis(1)=one
          axis(2)=zero
          axis(3)=zero
          cgamma=-one
          sgamma=zero
          go to 500
        end if
      end if
      axis(1)=v(2)/vv
      axis(2)=-v(1)/vv
      axis(3)=zero
      cgamma=v(3)
      sgamma=sqrt(one-cgamma**2)
c rotating v counterclockwise around axis moves it to +z
c now rotate u
  500 continue
      call rotax(u,axis,cgamma,sgamma)
c determine the azimuthal angle of u after rotation, and set up
c  the w vector which includes theta angle with +z and tau torsion
c  with u
      call azimuth(u,r,theta,tau,w)
c  now rotate w back
      call rotax(w,axis,cgamma,-sgamma)
 1000 continue
      w(1)=w(1)+c(1)
      w(2)=w(2)+c(2)
      w(3)=w(3)+c(3)
      end
c ====================================================================
c
      subroutine azimuth(u,r,theta,tau,w)
c  this routine determines the azimuthal (phi) angle the vector
c  u makes with the positive z axis. It then constructs the
c  vector w as a vector of length r, polar angle theta and making
c  a torsional angle tau with u
c parameters: input u(3),r,theta,tau
c             output w(3)
      implicit real*8 (a-h,o-z)
      dimension u(3),w(3)
      data zero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/
      uu=sqrt(u(1)**2+u(2)**2)
      if(uu.lt.1.0d-4) then
          call nerror(1,'azimuth',
     1    'vectors in azimuth are almost collinear',0,0)
      end if
      pi2=two*atan(one)
      pi=two*pi2
      if(abs(u(1)).gt.1.0d-6) then
        phi=atan(u(2)/u(1))
        if(u(1).lt.zero) phi=phi+pi
      else
        if(u(2).gt.zero) then
          phi=pi2
        else
          phi=-pi2
        end if
      end if
      w(1)=r*sin(theta)*cos(phi+tau)
      w(2)=r*sin(theta)*sin(phi+tau)
      w(3)=-r*cos(theta)
      end
c ====================================================================
c
      subroutine rotax(x,axis,cgamma,sgamma)
c  this routine rotates the vector x around the  vector axis
c  sgamma and cgamma are the sine and cosine of the rotational angle
c  parameters: input: x(3) vector to be rotated, replaced by the result
c                     axis(3): unit vector of rotational axis
c  output: x
      implicit real*8 (a-h,o-z)
      dimension x(3),axis(3),y(3),z(3)
      proj=scalar(axis,x)
      do l=1,3
        z(l)=proj*axis(l)
      end do
      call crosspr(axis,x,y)
      do l=1,3
        x(l)=z(l)+cgamma*(x(l)-z(l))+sgamma*y(l)
      end do
      end
c ====================================================================
c
      subroutine crosspr(a,b,c)
c  c(3) is the cross (vector) product of a(3) and b(3)
c  input a,b; output c
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),c(3)
      c(1)=a(2)*b(3)-a(3)*b(2)
      c(2)=a(3)*b(1)-a(1)*b(3)
      c(3)=a(1)*b(2)-a(2)*b(1)
      end
c ====================================================================
c
      subroutine bondvec(a,b,v)
c  builds the unit vector directed from a to b in v
c  parameters
c  input: a(3),b(3)
c  output: v(3)
c  a and b should not coincide
      implicit real*8 (a-h,o-z)
      dimension a(3),b(3),v(3)
      data zero,half,one,two/0.0d0,0.5d0,1.0d0,2.0d0/
      s=zero
      do l=1,3
        v(l)=b(l)-a(l)
        s=s+v(l)**2
      end do
      s=sqrt(s)
      if(s.gt.1.0d-6) then
        do l=1,3
          v(l)=v(l)/s
        end do
      else
            call nerror(1,'bondvec',
     1      'bondvec called with coinciding atoms',0,0)
      end if
      end
c ====================================================================
c
      SUBROUTINE DefMASSN(IAN,AtMASS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Assign default atomic mass
C
      parameter (nelem=92)
      dimension amass(nelem)
c
c --- define atomic weights for the first <nelem> elements ---
c --- last update:  12.12.97    JB   ** CHECKED UP TO Xe (element 54) **
c --- data from "Quantities, Units and Symbols in Physical Chemistry"
c --             IUPAC. second edition, 1993
c
c                         H - Zr
      data (amass(i), i = 1 , 40) /
     1  1.00794d+00,  4.0026d+00,  6.941d+00, 9.012182d+00, 10.811d+00,
     2 12.01115d+00, 14.00674d+00,15.9994d+00,18.998403d+00,20.1797d+00,
     3 22.989768d+00,24.305d+00,26.981539d+00,28.0855d+00,30.973762d+00,
     4 32.066d+00,   35.4527d+00, 39.948d+00,  39.0983d+00, 40.078d+00,
     5 44.95591d+00, 47.88d+00,  50.9415d+00, 51.9961d+00,54.938047d+00,
     6 55.847d+00,  58.933198d+00, 58.34d+00,  63.546d+00,  65.39d+00,
     7 69.723d+00,   72.61d+00,   74.921594d+00, 78.96d+00, 79.904d+00,
     8 83.80d+00,    85.4678d+00, 87.62d+00, 88.905849d+00, 91.224d+00 /
c                         Nb - Hg
      data (amass(i), i = 41 , 80) /
     1  92.906377d+00, 95.94d+00,97.907215d+00,101.07d+00, 102.9055d+00,
     2 106.42d+00,  107.8682d+00,112.411d+00,  114.818d+00, 118.71d+00,
     3 121.757d+00, 127.60d+00, 126.904473d+00, 131.29d+00, 132.905d+00,
     4 137.33d+00,  138.91d+00,  140.115d+00,  140.908d+00 ,144.24d+00,
     5 146.92d+00,  150.36d+00,  151.965d+00,  157.25d+00,  158.925d+00,
     6 162.50d+00,  164.93d+00,  167.26d+00,   168.93d+00,  173.04d+00,
     7 174.97d+00,  178.49d+00,  180.95d+00,   183.85d+00,  186.21d+00,
     8 190.2d+00,   192.22d+00,  195.08d+00,   196.07d+00,  200.59d+00 /
c                         Tl - U
      data (amass(i), i = 81 , 92) /
     1 204.38d+00,  207.2d+00,   208.980d+00,  208.98d+00,  209.99d+00,
     2 222.02d+00,  223.02d+00,  226.03d+00,   227.03d+00,  232.04d+00,
     3 231.04d+00,  238.03d+00  /
c
      if(ian.lt.1.or.ian.gt.nelem) then
        atmass=0.0d0
      else
        atmass=amass(ian)
      end if
c
      END
c ====================================================================
c
      SUBROUTINE ParseGEOM(IUnit,IStyl,NCentre)
      IMPLICIT INTEGER(A-Z)
C
C  This routine parses the input file (already open on unit IUnit)
C  to get a very good estimate of how many atomic centres there are
C
C  ARGUMENTS
C
C  IUnit   -  unit number for input file
C  IStyl   -  input style for geometry
C  NCentre -  on exit number of atomic centres OR default maximum
C
C  --------------------------------------------------------------
C  ** WARNING **  Currently only PQS and TEXAS input is checked
C                 All other input formats default
C  --------------------------------------------------------------
C
      PARAMETER (NumC=9)
      CHARACTER*80 Char
      CHARACTER*4 Keyword(NumC)
c
      Data Keyword/'basi','semi','ffld','numh','opti',
     $             'scan','path','scf ','freq'/
C
C
C  set up default maximum number of atomic centres
C
      NCentre = 2000
c
      If(IStyl.NE.1.AND.IStyl.NE.2) RETURN
c
      NCentre = 0
      nlines=0
c
 10   CONTINUE
      nlines=nlines+1
      READ(IUnit,900,End=95) Char
      call lowerca2(Char,4)
C
C  if <Char> starts with one of the keywords, then there is no
C  more geometry input
C
      DO I=1,NumC
      If(Char(1:4).EQ.Keyword(I)) GO TO 100
      EndDO
c
      NCentre = NCentre+1
      GO TO 10
C
 100  CONTINUE
C
C  we've finished
C  need to backspace input file to beginning of geometry
C
      DO I=1,nlines+1
      BACKSPACE IUnit
      EndDO
c
 20   CONTINUE
      READ(IUnit,900,End=95) Char
      call lowerca2(Char,4)
      If(Char(1:4).EQ.'geom'.OR.Char(1:4).EQ.'nucl') RETURN
      GO TO 20
C
C  ...........................................
C  Error handling
C
 95   CONTINUE
      Call nerror(6,'GEOMETRY module',
     $  'Problems determining number of atoms - check your input',0,0)
C
  900 Format(A80)
c
      END

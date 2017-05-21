      subroutine setdata
c  sets the constants
      implicit real*8 (a-h,o-z)
c  this is retained for historical reasons but it is better to use the
c  options
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c-----------------------------------------------------------------------
c commons for integral package:
c
      common /intlim/ limxmem,limblks,limpair
      common /intgop/ ncache,maxprice,iprint,iblock
      common /outfile/ ioutput
      common /route/ iroute
c-----------------------------------------------------------------------
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      five=5.0d0
      ten=10.0d0
      ten6=1.0d6
      tenm8=1.0d-8
      pi=3.14159 26535 89793d0
      call setrval('zero',zero)
      call setrval('half',half)
      call setrval('one',one)
      call setrval('two',two)
c  note spelling
      call setrval('tree',three)
      call setrval('four',four)
      call setrval('five',five)
      call setrval('ten',ten)
      call setrval('10+6',ten6)
      call setrval('10-8',tenm8)
      call setrval('pi',pi)
      call getival('accu',ii)
      acc=10.0d0**(-ii)
      call setrval('acc',acc)
c  default time limit is 10 M seconds=115 days
      call setrval('time',1.0d7)
c
c  Plack's constant, speed of light,electron charge and mass,
c  atomic mass unit, vacuum permittivity
c  all in SI units
      hbar=1.054 572 66d-34                   ! J*s  Plancks' const.
      c=2.99792458d8                          ! m/s speed of light
      enul=1.602 177 33d-19                   ! C (Coulomb) e
      xme=0.910 938 97d-30                    ! kg (electron mass)
      amu=1.660 540 2d-27                     ! kg (atomic mass unit)
      eps0=8.854 187 816d-12                  ! vacuum permittivity
      avog=6.022 136 7d23                     ! Avogadro-Loschmidt no.
      bolt=1.380 658d-23                      ! Boltzmann const
      bohr=4.0d0*pi*eps0*hbar**2/(xme*enul**2)! m  Bohr radius
      angs=1.0d-10/bohr                       ! Angstrom/bohr, dim.less
      hartree=enul**2/(bohr*4.0d0*pi*eps0)    ! J atomic energy unit
      ajoule=1.0d-18/hartree                  ! hartree/aJ, dim.less
      evolt=enul/hartree                      ! eV/hartree=1/27.21
c kJ/mol
      xkjmol=1.0d3/(avog*hartree)             ! (kJ/mol)/hartree=1/2550
      xkcal=xkjmol*4.184d0                    ! (kcal/mol)/hartree=1/627
      dkel=bolt/hartree                       ! k*1K/(Hartree) k= Boltz
      cmm1=200.0d0*pi*hbar*c/hartree          ! 1 (hc/cm)/hartree=1/219k
      hertz=200*pi*hbar/hartree               ! 1 (h/s)/hartree
c  Angstrom, Debye, Cb*m,
      call setrval('hbar',hbar)
      call setrval('c   ',c)
      call setrval('enul',enul)
      call setrval('me  ',xme)
      call setrval('amu ',amu)
      call setrval('eps0',eps0)
      call setrval('avog',avog)
      call setrval('bolt',bolt)
      call setrval('angs',angs)
      call setrval('deby',0.39342658d0)
      call setrval('cbm',0.117946332d30)
      call setrval('ajou',ajoule)
      call setrval('evol',0.036749309d0)
      call setrval('kcal/mol',xkcal)
      call setrval('kJ/mol',xkjmol)
      call setrval('dkel',dkel)
      call setrval('cmm1',cmm1)
      call setrval('Hertz',hertz)
c
c
c-----------------------------------------------------------------------
c set the default parameters for two-el. inetgrals
c
      threshold=1.0d-10
      call setup_thres(threshold)   ! setup common /neglect/ eps,eps1,epsr
c
      maxprice=1
      ncache=1
      iroute=0
      icheck=0
      iprint=0
      istat= 91         ! istat= statistics file, temporarily 91
      ioutput=istat
      igran= 20         ! granularity for parallel mode
c
c     blocking parameters
c
c     limxmem=500 000
      limxmem=300 000
c2002 limblks=300
      limblks=0
      limpair=150
c
      call setival('maxprice',maxprice)
      call setival('ncache',ncache)
      call setival('iroute',iroute)
      call setival('chec',icheck)
      call setival('iprn',iprint)
      call setival('istat',istat)
      call setival('gran',igran)
      call setival('limxmem',limxmem)
      call setival('limblks',limblks)
      call setival('limpair',limpair)
c2002
c  set default for instability
      call setival('stab',+1)       ! assume everything is STABLE
c
      end
c===================================================================
      subroutine setconst
      implicit real*8(a-h,o-z)
c
c -- sets various constants (in atomic units)
c -- see IUPAC handbook (second edition, 1993)
c -- communicates via common block /CONSTANTS/
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
c
c
      PI = 3.14159 26535 89793d0
      hbar = 1.054 572 66d-34                ! Planck constant/2PI Js
      c = 2.99792458d+10                     ! speed of light (cm/s)
      enul = 1.602 177 33d-19                ! elementary charge C
      xme = 9.10 938 97d-31                  ! electron rest mass kg
      amu = 1.660 540 2d-27                  ! atomic mass unit kg
      eps0 = 8.854 187 816d-12               ! vacuum permittivity
      avogad = 6.022 136 7d+23               ! Avogadro's number
      boltz = 1.380 658d-23                  ! Boltzmann constant J/(K*molecule)
      bohr = 4.0d0*PI*eps0*hbar**2/(xme*enul**2)
      hartree = enul**2/(bohr*4.0d0*pi*eps0) ! hartree (energy unit in au) J
      caljou = 4.184d+03                     ! kcalories --> joules
      R = 8.3144126d0                        ! gas constant J/(K*mol)
c
c -- conversion factor angstroms --> au
      ANTOAU = 1.0d-10/bohr
c
      end
c===================================================================
      subroutine setdatan
c  sets the constants
      implicit real*8 (a-h,o-z)
c  this is retained for historical reasons but it is better to use the
c  options
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
c-----------------------------------------------------------------------
c commons for integral package:
c
      common /intlim/ limxmem,limblks,limpair
      common /intgop/ ncache,maxprice,iprint,iblock
      common /outfile/ ioutput
      common /route/ iroute
c-----------------------------------------------------------------------
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      five=5.0d0
      ten=10.0d0
      ten6=1.0d6
      tenm8=1.0d-8
      pi=3.14159 26535 89793d0
      call setrval('zero',zero)
      call setrval('half',half)
      call setrval('one',one)
      call setrval('two',two)
c  note spelling
      call setrval('tree',three)
      call setrval('four',four)
      call setrval('five',five)
      call setrval('ten',ten)
      call setrval('10+6',ten6)
      call setrval('10-8',tenm8)
      call setrval('pi',pi)
      call getival('accu',ii)
      acc=10.0d0**(-ii)
      call setrval('acc',acc)
c  default time limit is 10 M seconds=115 days
      call setrval('time',1.0d7)
c
c  Plack's constant, speed of light,electron charge and mass,
c  atomic mass unit, vacuum permittivity
c  all in SI units
      hbar = 1.054 571 628d-34                ! Planck constant/2PI Js
      c=2.99792458d8                          ! m/s speed of light
      enul = 1.602 176 487d-19                ! elementary charge C
      xme = 9.109 382 15d-31                  ! electron rest mass kg
      amu = 1.660 538 782d-27                  ! atomic mass unit kg
      eps0 = 8.854 187 817 620d-12               ! vacuum permittivity
      avog = 6.022 141 79d+23               ! Avogadro's number
      bolt = 1.380 6504d-23                ! Boltzmann constant J/(K*molecule)
      bohr=4.0d0*pi*eps0*hbar**2/(xme*enul**2)! m  Bohr radius
      angs=1.0d-10/bohr                       ! Angstrom/bohr, dim.less
      hartree=enul**2/(bohr*4.0d0*pi*eps0)    ! J atomic energy unit
      ajoule=1.0d-18/hartree                  ! hartree/aJ, dim.less
      evolt=enul/hartree                      ! eV/hartree=1/27.21
c kJ/mol
      xkjmol=1.0d3/(avog*hartree)             ! (kJ/mol)/hartree=1/2550
      xkcal=xkjmol*4.184d0                    ! (kcal/mol)/hartree=1/627
      dkel=bolt/hartree                       ! k*1K/(Hartree) k= Boltz
      cmm1=200.0d0*pi*hbar*c/hartree          ! 1 (hc/cm)/hartree=1/219k
      hertz=200.0d0*pi*hbar/hartree               ! 1 (h/s)/hartree
c  Angstrom, Debye, Cb*m,
      call setrval('hbar',hbar)
      call setrval('c   ',c)
      call setrval('enul',enul)
      call setrval('me  ',xme)
      call setrval('amu ',amu)
      call setrval('eps0',eps0)
      call setrval('avog',avog)
      call setrval('bolt',bolt)
      call setrval('angs',angs)
      call setrval('deby',0.39342658d0)
      call setrval('cbm',0.117946332d30)
      call setrval('ajou',ajoule)
      call setrval('evol',0.036749309d0)
      call setrval('kcal/mol',xkcal)
      call setrval('kJ/mol',xkjmol)
      call setrval('dkel',dkel)
      call setrval('cmm1',cmm1)
      call setrval('Hertz',hertz)
c
c
c-----------------------------------------------------------------------
c set the default parameters for two-el. inetgrals
c
      threshold=1.0d-10
      call setup_thres(threshold)   ! setup common /neglect/ eps,eps1,epsr
c
      maxprice=1
      ncache=1
      iroute=0
      icheck=0
      iprint=0
      istat= 91         ! istat= statistics file, temporarily 91
      ioutput=istat
      igran= 20         ! granularity for parallel mode
c
c     blocking parameters
c
c     limxmem=500 000
      limxmem=300 000
c2002 limblks=300
      limblks=0
      limpair=150
c
      call setival('maxprice',maxprice)
      call setival('ncache',ncache)
      call setival('iroute',iroute)
      call setival('chec',icheck)
      call setival('iprn',iprint)
      call setival('istat',istat)
      call setival('gran',igran)
      call setival('limxmem',limxmem)
      call setival('limblks',limblks)
      call setival('limpair',limpair)
c2002
c  set default for instability
      call setival('stab',+1)       ! assume everything is STABLE
c
      end
c===================================================================
      subroutine setconstn
      implicit real*8(a-h,o-z)
c
c -- sets various constants (in atomic units)
c -- see IUPAC handbook (second edition, 1993)
c -- communicates via common block /CONSTANTS/
c
      COMMON /CONSTANTS/ PI,ANTOAU,hbar,c,enul,xme,amu,eps0,
     $                   avogad,boltz,hartree,caljou,R
c
c
      PI = 3.14159 26535 89793d0
      hbar = 1.054 571 628d-34                ! Planck constant/2PI Js
      c = 2.99792458d+10                     ! speed of light (cm/s)
      enul = 1.602 176 487d-19                ! elementary charge C
      xme = 9.109 382 15d-31                  ! electron rest mass kg
      amu = 1.660 538 782d-27                  ! atomic mass unit kg
      eps0 = 8.854 187 817 620d-12               ! vacuum permittivity
      avogad = 6.022 141 79d+23               ! Avogadro's number
      boltz = 1.380 6504d-23                  ! Boltzmann constant J/(K*molecule)
      bohr = 4.0d0*PI*eps0*hbar**2/(xme*enul**2)
      hartree = enul**2/(bohr*4.0d0*pi*eps0) ! hartree (energy unit in au) J
      caljou = 4.184d+03                     ! kcalories --> joules
      R = 8.314 472d0                        ! gas constant J/(K*mol)
c
c -- conversion factor angstroms --> au
      ANTOAU = 1.0d-10/bohr
c
      end
c===================================================================

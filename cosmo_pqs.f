c======================================================================
c
c  MM October2003-February 2004
c
c  this file contains subroutines that are specific to
c  the PQS implementation of COSMO.
c
c  Note that the current version of the COSMO implementation
c  does not support ghost atoms, although some steps have
c  already been taken to use them in future versions.
c
c  The parameters of the COSMO calculation are shared
c  using the depository system implemented in the PQS code.
c
c  Here follows a description of the parameters
c
c     Name         type                comments
c
c     cosmo       integer       Global cosmo flag. It is set to
c                               a value different from zero in
c                               case of a cosmo run. For a
c                               normal (no cosmo) run, it is either
c                               set to zero, or not set at all.
c
c    c_outfil     Char*256      Name of the cosmo output file
c    c_onel       Char*256      Name of the cosmo matrix elements file.
c                               If it is set to 'direct', integrals
c                               will not be saved, but they will
c                               be recomputed as needed.
c    c_data       Char*256      Name of the cosmo data file
c    c_sude       Char*256      Name of the file for surface derivatives
c    c_radi       Char*256      Name of the file for user defined radii
c    c_cradt      Char*256      type of radii to be used
c    c_csol       Char*256      Solvent name
c    c_info       Char*256      Info string for cosmo output
c
c    c_nratom     integer       Number of real atoms (no ghost atoms)
c                                Note: at the moment ghost atoms
c                                are not supported, thus if nratom
c                                not equal natom, the program stops
c    c_nppa       integer       Total number of basis points
c    c_nspa       integer       Number of segment for heavy atoms
c    c_nsph       integer       Number of segment for hydrogen atoms
c    c_lcavit     integer       Type of cavity 1=closed, 0=open
c    c_maxnps     integer       Maximum number of surface segments
c    c_nps        integer       Number of surface segments
c    c_nip        integer       Number of intersection pairs
c    c_npsphe     integer       Number of segments of outer surface
c    c_npsd       integer       = nps + npspher
c    c_gran       integer       Granularity for parallel execution
c
c    c_eps         real         Dielectric constant (it should be
c                               infinity for a COSMO-RS calculation)
c    c_fepsi       real         =(eps-1.0)/(eps+0.5) (=1.0 if eps=inf)
c    c_disex       real         cutoff distance for A matrix calculation
c    c_disex2      real         =(disex*mean_atomic_diameter)**2
c    c_routf       real         factor for outer sphere construction
c    c_phsran      real         phase offset for coordinate randomizat.
c    c_ampran      real         amplitude factor for coord. randomizat.
c    c_rsolv       real         solvent radius for surface construction
c    c_area        real         surface area
c    c_volume      real         cavity volume
c    c_qsum        real         cosmo charge
c    c_dq          real         outlying charge (OC)
c    c_ediel       real         dielectric energy
c    c_de          real         OC correction to enery
c    cpucosmo      real         cumulative slave CPU time for a parallel
c                               cosmo step. Currently used to add the
c                               time for the parallel cosmo gradient
c                               to the total slave time for the gradient
c                               step. Might be undefined, so test it
c                               before using it.
c
c   The following quantities are pointers to arrays used by COSMO.
c   The arrays are built into the PQS array system.
c   NOTE: at the moment these pointers are meaningfull only
c         within the SCF module.
c
c    name    type of array   dimensions             comments
c
c    c_israd      real       natom          COSMO radii
c    c_inuc       integer    natom          nuclear charge
c    c_icosur     real       3,2*maxnps     coord. of surface segment
c    c_iar        real       2*maxnps       area of surface segment
c    c_iiatsp     integer    2*maxnps       atom associated to segment
c    c_ia1mat     real       nps*(nps+1)/2  A matrix
c    c_iqcos      real       nps            surface charges
c    c_iphi       real       nps            surface potentials
c    c_iphin      real       nps            nuclear part of surface
c                                           potentials
c    c_iphio      real       npspher        outer surface potentials
c    c_iphic      real       nps            OC correction to potential
c    c_qcosc      real       nps            OC correction to charge
c
c
c   a sequential file is used to store and communicate data between
c   the program modules.
c   The name of the cosmo data file is stored in c_data, and the
c   structure of the file is the following:
c
c     write(ich)xyz       nuclear coordinates (real atoms only),
c                         real array (3,nratom)
c     write(ich)icharge   nuclear charges (real atoms only),
c                         integer array (nratom)
c     write(ich)cosurf    coordinates of surface segments,
c                         real array (3,nps)
c     write(ich)ar        area of surface segments, real array (nps)
c     write(ich)iatsp     mapping surface segment -> atoms
c                         integer array (nps)
c     write(ich)qcos      surface charges (unscaled and uncorrected)
c                         real array (nps)
c
c
c  Description of the input of the COSMO module:
c
c  COSMO command
c
c  The COSMo command requests a calculation using the Conductor-like
c  Screening Model (COSMO) [1] to model the effect of a solvent on the
c  system under study. In this model, the solute forms a cavity in a
c  dielectric continuum representing the solvent. The size of the cavity
c  is defined by the solvent accessible surface (SAS), which is
c  constructed on the basis of the molecular geometry.
c
c  Currently the COSMO model is available for SCF and FORCE calculations
c  using HF and DFT methods (second derivatives are accessible via
c  the NUMHESS command). The presence of a COSMo keyword in the
c  input sequence will affect all the SCF and FORCE commands
c  following it. Note that the COSMO model performs at his best for
c  polar solvents (large permittivity values), while results
c  for weak dielectrics might be less reliable.
c
c  The default settings of the COSMO module have been tailored
c  for the application of the COSMO-RS (COSMO for Real Solvents)
c  theory [2]. In particular, the COSMO module will produce an
c  output file (<jobname>.cosmo) that can be used as input for
c  the COSMOtherm suite of programs for the calculation of solvation
c  mixture thermodynamics. The COSMOterm software is distributed
c  by CSOMOlogic GmbH & Co.KG (www.cosmologic.de).
c  The COSMO-RS parameters have been optimized for DFT calculations
c  using the BP86 functional and the svp_ahlrichs and tzvp_ahlrichs
c  basis sets, thus the recommended settings for running a COSMO-RS
c  calculation are DFTP=BP86 in the SCF card, and BASIS=svp_ahlrichs
c  or BASIS=tzvp_ahlrichs as basis set choiche.
c
c  Options:
c
c  SOLV=<string>: The name of the solvent to be used. The default
c                 corresponds to the settings for a COSMO-RS
c                 calculation (in other words, for a proper COSMO-RS
c                 calculation, the option SOLV must not be present).
c                 Specifying a solvent will set the values of
c                 the options EPSI and RSOL (see below) to the
c                 predefined  value for the chosen solvent,
c                 according to the following table. Note that
c                 the values of EPSI and RSOL that will actually be
c                 used in the calculation can be changed using
c                 the appropriate options (see below).
c
c                   name           alias(es)      rsol (A)      epsi
c
c           water                 h2o              1.385       78.39
c           methanol              ch3oh            1.855       32.63
c           ethanol               c2h5oh ch3ch2oh  2.180       24.55
c           acetone               ch3coch3         2.38        20.7
c           ether                 ch3ch2och2ch3    2.785        4.335
c           acetonitrile          ch3cn            2.155       36.64
c           nitromethane          ch3no2           2.155       38.20
c           carbontetrachloride   ccl4             2.685        2.228
c           chloroform            chcl3            2.48         4.90
c           dichloromethane       ch2cl2           2.27         8.93
c           dichloroethane        ch2clch2cl       2.505       10.36
c           dimethylsulfoxide     dmso             2.455       46.7
c           aniline               c6h5nh2          2.80         6.89
c           benzene               c6h6             2.63         2.247
c           toluene               c6h5ch3          2.82         2.379
c           chlorobenzene         c6h5cl           2.805        5.621
c           tetrahydrofuran       thf              2.56         7.58
c           cyclohexane           c6h12            2.815        2.023
c           n-heptane             heptane   c7h16  3.125        1.92
c
c                Solvents that are not listed in the table can be
c                simulated by explicitly entering the corresponding
c                values of EPSI and RSOL as indicated below
c
c  EPSI=<real>: the electric permittivity of the dielectric continuum.
c               The default is infinity if no solvent is specified,
c               otherwise it is the value for the chosen solvent.
c  RSOL=<real>: Additional radius (A) for SAS construction. The
c               default value is 1.3 A if no solvent is specified,
c               oterwise it is the value for the chosen solvent.
c  RADI=<string>: Type of atomic radii to be used for SAS construction
c        COSMO   (default) use the COSMO optimized radius if available,
c                otherwise use 1.17*Bondi radius
c        BONDI   Use van der Waals radii from A.Bondi, J.Phys.Chem. 68,
c                441-451 (1964).
c                Note that the internally defined atomic radii that
c                can be accessed with the <COSMO> and <BONDI> options
c                do not cover all the periodic table (although all the
c                commonly used elements are defined). If the
c                program encounters an atom for which an atomic
c                radius cannot be defined, it will stop with an error.
c                To proceed in such a case, the corresponding radius
c                must be entered as user defined value with one
c                of the three methods described below.
c        USER    (requires additional input) user-defined radii. The
c                user defined radii are read from one or more additional
c                input lines immediately following the COSMo keyword.
c                These input lines must begin with the sequence of
c                characters $RADI and are formed by one or more fields
c                of the type <symbol>=<radius>, where <symbol> is a
c                string identifying an atomic center (the use of numbers
c                or special characters to match the symbols used in
c                the geometry section is permitted) and <radius> is
c                the corresponding atomic radius (A).
c                Note that the user defined radii do not need to cover
c                all the atomic centers of the system under study.
c                Every atomic radius still undefined after parsing
c                the user defined input, will be assigned the default
c                value.
c        USERB   (requires additional input). The same as the previous
c                case, except that every atom not defined by input will
c                be given the Bondi radius.
c     <filename> If the string following the RADI option does not match
c                any of the four previous cases, it will be assumed to
c                indicate the name of a file containing the user
c                specified atomic radii. The format of this user
c                supplied file is the same as for the additional input
c                lines required by the USER case, except that the
c                initial $RADI sequence may be omitted (if present,
c                it will be ignored). The existence of the user supplied
c                file will be tested, and an error generated if it is
c                not found.
c
c            Examples of user specified atomic radii:
c
c            --  GEOM=PQS
c                H       0.533943571    0.000000000   -0.763962625
c                O      -0.067278458    0.000000000    0.000001331
c                H       0.533960196    0.000000000    0.763941500
c                COSMO  radi=USER
c                $radi h=1.5 o=1.8
c
c                will assign a radius of 1.5 A to hydrogen, and 1.8 A
c                to oxigen
c
c            --  GEOM=PQS
c                H       0.533943571    0.000000000   -0.763962625
c                O      -0.067278458    0.000000000    0.000001331
c                H       0.533960196    0.000000000    0.763941500
c                COSMO  radi=USER
c                $radi h=1.5
c
c                will assign a radius of 1.5 A to hydrogen, oxigen
c                will have the default value (1.72 A)
c
c            --  GEOM=PQS
c                H       0.533943571    0.000000000   -0.763962625
c                O      -0.067278458    0.000000000    0.000001331
c                H1      0.533960196    0.000000000    0.763941500
c                COSMO  radi=USER
c                $radi h=1.5 o=1.8
c                $radi h1=1.3
c
c                first hydrogen will have 1.5 A, oxigen 1.8 A,
c                second hydrogen (h1) will have 1.3 A
c
c            --  GEOM=PQS
c                H       0.533943571    0.000000000   -0.763962625
c                O      -0.067278458    0.000000000    0.000001331
c                H1      0.533960196    0.000000000    0.763941500
c                COSMO  radi=USER
c                $radi h1=1.3
c
c                first hydrogen  and oxygen will have default value,
c                second hydrogen (h1) will have 1.3 A
c
c            --  GEOM=PQS
c                H       0.533943571    0.000000000   -0.763962625
c                O      -0.067278458    0.000000000    0.000001331
c                H1      0.533960196    0.000000000    0.763941500
c                COSMO  radi=USER
c                $radi h=1.3
c
c                both hydrogens will have 1.3 A, oxygen will have
c                default
c
c  The following options are for fine tuning, and are meant for
c  advanced users:
c
c  ROUT=<real>: Factor for outer sphere construction (default 0.85).
c               The outer sphere is used for the outlying charge
c               correction.
c  DISE=<real>: Cutoff factor for use of basis grid point
c               during SAS construction (default 10).
c  LCAV=<integer>: Type of cavity. 1=closed cavity, 0=open cavity.
c                  default: 1 (closed cavity)
c  NPPA=<integer>: Total number of basis grid points (default 1082).
c  NSPA=<integer>: Number of segments for non hydrogen atoms
c                  (default 92).
c  AMPR=<real>: Amplitude factor for coordinate randomization during
c               SAS construction. Needs to be 0.00001 (default) or
c               smaller.
c  PHSR=<real>: Phase offset for coordinate randomization (default 0.0)
c  ONEL       : This flag will switch on the direct calculation
c               of the COSMO one electron matrix elements. If this
c               option is present, the COSMO one-electron integrals
c               will not be saved on file, but they will be recomputed
c               every time they are needed (that is, two times for
c               each SCF iteration). The default is to save the
c               integrals on file <jobname>.cosmo_onel. This file
c               will require nps*ncf*(ncf+1)*4 bytes of disk storage,
c               where ncf is the number of basis functions and nps
c               is the number of surface segments. The current
c               implementation of the direct calculation of the COSMO
c               matrix elements is rather time consuming, thus it is
c               recommended not to use the ONEL option, unless there is
c               not enough disk space available for the default
c               settings.
c  GRAN=<integer>: Granularity for parallel execution of the COSMO
c                  algorithm (default 5). Needs to be 1 or bigger. Very
c                  large jobs might benefit from a smaller granularity.
c
c  References
c
c  [1] a) A. Klamt and G. Schuurmann, J. Chem. Soc. Perkin Trans. II
c         799 (1993).
c      b) J. Andzelm, C. Kolmel and A. Klamt,
c         J. Chem. Phys. 103, 9312 (1995).
c      c) A. Klamt and V. Jonas, J. Chem. Phys. 105, 9972 (1996).
c      d) K. Baldridge and A. Klamt, J. Chem. Phys. 106, 66622 (1997).
c      e) A. Schafer, A. Klamt, D. Sattel, J.C.W. Lohrenz and
c         F. Eckert, Phys. Chem. Chem. Phys. 2 2187 (2003).
c
c  [2] a) A. Klamt, J. Phys. Chem. 99, 2224 (1995).
c      b) A. Klamt, V. Jonas, T. Burger, and J.C.W. Lohrenz,
c         J. Phys. Chem. A 102, 5074 (1998).
c
c======================================================================
      subroutine cosmo_init(jobname,lenj)
      implicit real*8 (a-h,o-z)
c
c   this subroutine reads the options of the COSMO calculation
c   and initialises the cosmo parameters.
c   The COSMO radii are not assigned here. This will be done
c   before the calculation is carried out.
c
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      parameter (nopt=14)
      character*4 word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      character*256 jobname,scrf
      character*256 chopv(nopt)
      character*256 outfile,info,conel,cdata,csude,cradi,cradts,csol
      character*128 inpline
      parameter (lwrd=128)
      character*5 cradt,cradr
c
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp/11,11,11,11,1,1,1,11,11,21,21,0,1,0/
c
      data word/'epsi','rsol','rout','dise','lcav','nppa','nspa','ampr',
     $          'phsr','radi','solv','onel','gran','off '/
      logical isthere,busy
c
c   get some physical constants
c
      call getrval('angs',angs)
c
c   defaults
c
      eps=-1.0d0        ! a value of -1.0 means eps=infinity
      fepsi=1.0d0       ! fepsi=(eps-1.0)/(eps+0.5) default is eps=infinity
      nppa=1082         ! total number of basis grid points
      nspa=92           ! number of segments on non hydrogen atoms
      disex=10.0d0      ! for distances smaller than disex*mean atomic
                        ! diameter, the A-matrix elements will be computed
                        ! using the basis grid points of the segment
      routf=0.85d0      ! factor for outer sphere construction
      lcavity=1         ! 1=closed cavity, 0=open cavity
      phsran=0.0d0      ! phase offset for coordinate randomization
      ampran=0.00001d0  ! amplitude factor for coordinate randomization
      rsolv=1.30d0*angs ! solvent radius
      csol='unknown'    ! solvent name
      cradt='cosmo'     ! type of radii
      outfile=jobname(1:lenj)//'.cosmo'
      call getchval('scrf',scrf)
      conel=scrf(1:len_trim(scrf))//'.cosmo_onel'
      cdata=scrf(1:len_trim(scrf))//'.cosmo_data'
      csude=scrf(1:len_trim(scrf))//'.cosmo_sude'
      cradi=scrf(1:len_trim(scrf))//'.cosmo_radi'
      info='PQS COSMO module - Implementation of 2 Feb. 2004'
      icgran=5          ! initial granularity for parallel execution
c
      call getival('iout',iout)
      call getival('icond',icon)
      write(iout,'(61(''=''),/,20x,a,/)')'COSMO initialization'
c
c  read input options
c
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
c
c -- first check if we are switching off COSMO
      If(ifound(14).gt.0) Then
        write(iout,*) ' Switching Off COSMO'
        write(icon,*) ' Switching off COSMO'
        call setival('cosmo',0)
        return
      EndIf
c
      if(ifound(3).gt.0)then
        routf=ropv(1,3)
        write(iout,'(6x,a,f10.4)')
     $      'Option '//word(3)//' is on, value=',routf
      endif
      if(ifound(4).gt.0)then
        disex=ropv(1,4)
        write(iout,'(6x,a,f10.4)')
     $      'Option '//word(4)//' is on, value=',disex
      endif
      if(ifound(5).gt.0)then
        lcavity=iopv(1,5)
        write(iout,'(6x,a,i6)')
     $      'Option '//word(5)//' is on, value=',lcavity
      endif
      if(ifound(6).gt.0)then
        nppa=iopv(1,6)
        write(iout,'(6x,a,i6)')
     $      'Option '//word(6)//' is on, value=',nppa
      endif
      if(ifound(7).gt.0)then
        nspa=iopv(1,7)
        write(iout,'(6x,a,i6)')
     $      'Option '//word(7)//' is on, value=',nspa
      endif
      if(ifound(8).gt.0)then
        ampran=ropv(1,8)
        write(iout,'(6x,a,e10.4)')
     $      'Option '//word(8)//' is on, value=',ampran
      endif
      if(ifound(9).gt.0)then
        phsran=ropv(1,9)
        write(iout,'(6x,a,e10.4)')
     $      'Option '//word(9)//' is on, value=',phsran
      endif
      if(ifound(10).gt.0)then
        cradr=chopv(10)(1:5)
        call lowerca2(cradr,5)
        if(cradr.eq.'cosmo'.or.cradr.eq.'bondi'.or.cradr.eq.'user '.or.
     $     cradr.eq.'userb')then
          cradt=cradr
          write(iout,'(6x,a,a)')
     $        'Option '//word(10)//' is on, value=',cradt
        else
          cradt='file '
          cradi=chopv(10)
          call leadblan2(cradi,len_trim(cradi),nlen)
          write(iout,'(6x,a,a)')
     $    'Option '//word(10)//' is on, value=',cradi(1:len_trim(cradi))
          inquire(file=cradi(1:len_trim(cradi)),exist=isthere)
          if(.not.isthere)then
            call nerror(1,'cosmo_init',
     $      'COSMO: file '//cradi(1:len_trim(cradi))//' not found',0,0)
          endif
        endif
        if(cradt.eq.'user '.or.cradt.eq.'userb')then
          do i=20,99
            inquire(i,opened=busy)
            if(.not.busy)then
               ich=i
               goto 10
             endif
          enddo
          call nerror(2,'cosmo_init','cannot find a free channel',0,0)
 10       continue
          open(unit=ich,file=cradi(1:len_trim(cradi)),status='unknown',
     $         form='formatted')
 20       call blankit(inpline,lwrd)
          read(inp,'(a)')inpline
          if(inpline(1:5).eq.'$radi')then
            write(ich,'(a)')inpline(6:lwrd)
            goto 20
          else
            backspace(inp)
          endif
          close(ich,status='keep')
        endif
      endif
      if(ifound(11).gt.0)then
        call blankit(csol,256)
        csol=chopv(11)
        call leadblan2(csol,len_trim(csol),nlen)
        call lowerca2(csol,len_trim(csol))
        write(iout,'(6x,a,a)')
     $  'Option '//word(11)//' is on, value=',csol(1:len_trim(csol))
        call cosmo_sol(csol,eps,rsolv)
      endif
      if(ifound(1).gt.0)then
        eps=ropv(1,1)
        write(iout,'(6x,a,f10.4)')
     $      'Option '//word(1)//' is on, value=',eps
      endif
      if(ifound(2).gt.0)then
        rsolv=ropv(1,2)*angs
        write(iout,'(6x,a,f10.4)')
     $      'Option '//word(2)//' is on, value=',rsolv/angs
      endif
      if(ifound(12).gt.0)then
         call blankit(conel,256)
         conel='direct'
        write(iout,'(6x,a)')
     $      'Option '//word(12)//' is on'
      endif
      if(ifound(13).gt.0)then
        icgran=iopv(1,13)
        write(iout,'(6x,a,i6)')
     $      'Option '//word(13)//' is on, value=',icgran
      endif
      write(iout,*)
c
c  compute fepsi
c
      if(eps.le.0.0d0)then
        fepsi=1.0d0
      else
        fepsi=(eps-1.0d0)/(eps+0.5d0)
      endif
c
c   recompute nspa and nsph
c
      x0=log(float(nspa)*0.1d0-0.199999d0)
      z3=log(3.0d0)
      z4=log(4.0d0)
      i4=int(x0/z4)
c
      nspb=0
      do i=0,i4
        x=x0-i*z4
        n=3**int(x/z3)*4**i
        if(n.gt.nspb)nspb=n
      enddo
      nsph=nspb/3
      if(mod(nspb,3).ne.0)nsph=nspb/4
      nspa=10*nspb+2
      nsph=max(12,nsph*10+2)
c
c  print COSMO parameters
c
      write(iout,'(6x,a)')'COSMO parameters:'
      if(csol.ne.'unknown')write(iout,'(9x,a)')
     $  'solvent                       = '//csol(1:len_trim(csol))
      if(eps.le.0.0d0)then
        write(iout,'(9x,a)')  'electric permittivity         = infinity'
      else
        write(iout,'(9x,a,f11.4)')'electric permittivity         =',eps
      endif
      write(iout,'(9x,a,f11.4)')  'scaling factor f(epsi)        =',
     $                             fepsi
      if(lcavity.eq.0)then
        write(iout,'(9x,a)')    'type of cavity                = open'
      else
        write(iout,'(9x,a)')    'type of cavity                = closed'
      endif
      write(iout,'(9x,a)')  'type of atomic radii          = '//cradt
      if(cradt.eq.'file ')write(iout,'(9x,a)')
     $  'User defined radii file       = '//cradi(1:len_trim(cradi))
      write(iout,'(9x,a,f11.4)')  'distance of solvent sphere (A)=',
     $                             rsolv/angs
      write(iout,'(9x,a,f11.4)')  'factor for outer sphere       =',
     $                             routf
      write(iout,'(9x,a,f11.4)')  'cutoff factor                 =',
     $                             disex
      write(iout,'(9x,a,e11.4)')  'phase offset for randomization=',
     $                             phsran
      write(iout,'(9x,a,e11.4)')  'amplitude factor for randomiz.=',
     $                             ampran
      if(conel.eq.'direct')then
        write(iout,'(9x,a)')'one electron integrals        = direct'
      else
        write(iout,'(9x,a)')'one electron integrals        = save'
      endif
      if(ampran.gt.1.0d-5)then
        call nerror(3,'cosmo_init','ampran must be <= 1.0d-5',0,0)
      endif
      write(iout,'(9x,a,i6)')     'no. of points per atom        =',nppa
      write(iout,'(9x,a,i6)')     'no. of segments per atom      =',nspa
      write(iout,'(9x,a,i6)')     'No. of segments for hydrogen  =',nsph
      call getival('nslv',nslv)
      if(nslv.ne.0)then
        write(iout,'(9x,a,i6)')'granularity                   =',icgran
        if(icgran.lt.1)then
           call nerror(4,'cosmo_init','granularity is less than 1',0,0)
        endif
      endif
      write(iout,*)
c
c  store COSMO parameters into depository
c
      call setrval('c_eps',eps)
      call setrval('c_fepsi',fepsi)
      call setrval('c_disex',disex)
      call setrval('c_routf',routf)
      call setrval('c_phsran',phsran)
      call setrval('c_ampran',ampran)
      call setrval('c_rsolv',rsolv)
c
      call setival('c_nppa',nppa)
      call setival('c_nspa',nspa)
      call setival('c_nsph',nsph)
      call setival('c_lcavit',lcavity)
      call setival('c_gran',icgran)
c
      call setchval('c_outfil',outfile)
      call setchval('c_onel',conel)
      call setchval('c_data',cdata)
      call setchval('c_sude',csude)
      call setchval('c_info',info)
      call blankit(cradts,256)
      cradts(1:5)=cradt
      call setchval('c_cradt',cradts)
      call setchval('c_cradi',cradi)
      call setchval('c_csol',csol)
c
c   switch on the global COSMO flag
c
      call setival('cosmo',1)
      write(icon,*)
      write(icon,*) ' Using COSMO Solvation Model'
c
      write(iout,'(61(''='')),')
      return
      end
c =================================================================

      subroutine cosmo_sol(csol,eps,rsolv)
      implicit real*8 (a-h,o-z)
c
c  data for various solvents
c
c  input
c
c   csol    solvent name
c
c  output
c
c   eps     permittivity
c   rsolv   solvent radius for SAS construction
c
      character*256 csol
c
c  conversion factor
c
      call getrval('angs',angs)
c
c  water
c
      if(csol.eq.'water'.or.csol.eq.'h2o')then
        rsolv=1.385d0*angs
        eps=78.39d0
c
c  methanol
c
      else if(csol.eq.'methanol'.or.csol.eq.'ch3oh')then
        rsolv=1.855d0*angs
        eps=32.63d0
c
c  ethanol
c
      else if(csol.eq.'ethanol'.or.csol.eq.'c2h5oh'
     $    .or.csol.eq.'ch3ch2oh')then
        rsolv=2.180d0*angs
        eps=24.55d0
c
c  acetone
c
      else if(csol.eq.'acetone'.or.csol.eq.'ch3coch3')then
        rsolv=2.38d0*angs
        eps=20.7d0
c
c  ether
c
      else if(csol.eq.'ether'.or.csol.eq.'ch3ch2och2ch3')then
        rsolv=2.785d0*angs
        eps=4.335d0
c
c  acetonitrile
c
      else if(csol.eq.'acetonitrile'.or.csol.eq.'ch3cn')then
        rsolv=2.155d0*angs
        eps=36.64d0
c
c  nitromethane
c
      else if(csol.eq.'nitromethane'.or.csol.eq.'ch3no2')then
        rsolv=2.155d0*angs
        eps=38.20d0
c
c  carbontetrachloride
c
      else if(csol.eq.'carbontetrachloride'.or.csol.eq.'ccl4')then
        rsolv=2.685d0*angs
        eps=2.2280d0
c
c  chloroform
c
      else if(csol.eq.'chloroform'.or.csol.eq.'chcl3')then
        rsolv=2.48d0*angs
        eps=4.90d0
c
c  dichloromethane
c
      else if(csol.eq.'dichloromethane'.or.csol.eq.'ch2cl2')then
        rsolv=2.27d0*angs
        eps=8.93d0
c
c  dichloroethane
c
      else if(csol.eq.'dichloroethane'.or.csol.eq.'ch2clch2cl')then
        rsolv=2.505d0*angs
        eps=10.36d0
c
c  dimethylsulfoxide
c
      else if(csol.eq.'dimethylsulfoxide'.or.csol.eq.'dmso')then
        rsolv=2.455d0*angs
        eps=46.7d0
c
c  aniline
c
      else if(csol.eq.'aniline'.or.csol.eq.'c6h5nh2')then
        rsolv=2.80d0*angs
        eps=6.89d0
c
c  benzene
c
      else if(csol.eq.'benzene'.or.csol.eq.'c6h6')then
        rsolv=2.63d0*angs
        eps=2.247d0
c
c  toluene
c
      else if(csol.eq.'toluene'.or.csol.eq.'c6h5ch3')then
        rsolv=2.82d0*angs
        eps=2.379d0
c
c  chlorobenzene
c
      else if(csol.eq.'chlorobenzene'.or.csol.eq.'c6h5cl')then
        rsolv=2.805d0*angs
        eps=5.621d0
c
c  tetrahydrofuran
c
      else if(csol.eq.'tetrahydrofuran'.or.csol.eq.'thf')then
        rsolv=2.56d0*angs
        eps=7.58d0
c
c  cyclohexane
c
      else if(csol.eq.'cyclohexane'.or.csol.eq.'c6h12')then
        rsolv=2.815d0*angs
        eps=2.023d0
c
c  eptane
c
      else if(csol.eq.'n-heptane'.or.csol.eq.'heptane'.or.
     $        csol.eq.'c7h16') then
        rsolv=3.125d0*angs
        eps=1.92d0
      else
        write(6,*)
        call nerror(1,'cosmo_sol',
     $  'solvent '//csol(1:len_trim(csol))//' not found',0,0)
      endif
      end
c =================================================================

      subroutine cosmo_parseradii(cradt,cradi,srad,natom)
      implicit real*8 (a-h,o-z)
c
c  this subroutine parses the user defined radii file to get
c  the atomic radii to be used in the calculation
c
c  cradt  type of the radii file (char*5), can be 'user ',or 'file '
c  cradi  radii file (char*256)
c  srad   in exit will contain the radii that have been assigned
c         by the user
c  natom number of atoms
c
c  format of the radii file:
c  a sequence of fields of the type xxxxxxxx=rad, where xxxxxxxx is an
c  extended atomic symbol (up to eight characters),
c  and rad is the corresponding radius (A).
c  there can be several fields per line (maximum length of the
c  line is 120 characters).
c  the first 5 characters of each line may contain the string
c  "$radi'.
c
      character*5 cradt
      character*256 cradi
      character*8 atsy,atsy1,xx
      dimension srad(natom)
      parameter (lwrd=120)
      character*120 line
      logical busy
      character*25 special
      data special/'~!@#$%^&*+=<>?"0123456789'/
c
c  look for a free channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(1,'cosmo_parseradii','cannot find a free channel',0,0)
 10   continue
c
c  open radii file
c
      open(unit=ich,file=cradi(1:len_trim(cradi)),status='old',
     $     form='formatted')
c
c  read a line
c
 100  call blankit(line,lwrd)
      read(ich,'(a)',end=999)line
c
c   lower case
c
      call lowerca2(line,lwrd)
      ist=1
c
c skip first 5 characters if $radi
c
      if(line(1:5).eq.'$radi')then
        line(1:5)='     '
        ist=6
      endif
c
c  look for '=' character
c
 200  if(line(ist:ist).eq.'=')then
c
c  get atomic symbol in xx
c
        ieq=ist
        ixe=ieq-1
        ixi=ixe
        do ic=ixe,1,-1
          if(line(ic:ic).eq.' ')then
            ixi=ic+1
            goto 210
          endif
        enddo
        ixi=1
 210    continue
        if(ixi.gt.ixe.or.ixe-ixi+1.gt.8)then
          write(6,*)'cosmo_parseradii, atomic symbol not valid in line'
          write(6,*)line
          write(6,*)ixi,ixe
          call nerror(1,'cosmo_parseradii','invalid atomic symbol',0,0)
        endif
        call blankit(xx,8)
        ixl=ixe-ixi+1
        xx(1:ixl)=line(ixi:ixe)
c
c  get corresponding radius in radxx
c
        iri=ieq+1
        ire=iri
        do ic=iri,lwrd
          if(line(ic:ic).eq.' ')then
            ire=ic-1
            goto 220
          endif
        enddo
        ire=lwrd
 220    continue
        if(iri.gt.ire)then
          write(6,*)'cosmo_parseradii, radius not valid in line'
          write(6,*)line
          write(6,*)iri,ire
          call nerror(2,'cosmo_parseradii','invalid radius',0,0)
        endif
        read(line(iri:ire),*)radxx
c
c   first we check the extended symbols for a match with xx,
c   and we aasign the corresponding radius
c
        do i=1,natom
          call getatsymex(i,natom,atsy1)
          if(atsy1.eq.xx)srad(i)=radxx
        enddo
c
c  if the lenght of xx is less or equal to 2, and
c  xx does not contain special characters, we test it also
c  against the standard symbols. We assign the radius only if
c  it has not already been assigned, so that standard symbols
c  will not overwrite extended symbols.
c
        if(len_trim(xx).le.2.and.index(special,xx).le.0)then
          do i=1,natom
            call getatsym1(i,natom,atsy)
            if(atsy.eq.xx.and.srad(i).eq.0.0d0)srad(i)=radxx
          enddo
        endif
c
c  ready for next field
c
        ist=ire
        if(ist.le.lwrd)goto 200
      else
        ist=ist+1
        if(ist.le.lwrd)goto 200
      endif
      goto 100
 999  continue
      close(unit=ich,status='keep')
      end
c =================================================================
c
      subroutine cosmo_radii(natom,IAN,srad,icharge,maxnps)

      use memory

      implicit real*8 (a-h,o-z)
C
C   This routine computes the maximum number of surface elements
C   and the COSMO radii for the atoms.
C
C  ARGUMENTS
C
C  natom   -  number of real atoms
C  IAN     -  atomic numbers
C
C  on exit
C
C  srad    -  COSMO radius of each atom
C  icharge -  array of integer atomic charges
C  maxnps  -  maximum number of surface elements
C
C  NOTES
C  The icharge array should be changed from integer to real
C  in which case we can use the existing Atomic charges array
C  and eliminate this altogether.  Currently we simply copy
C  the atomic numbers into icharge.
C
C  Dummy atoms are NOT used in COSMO
c
      dimension IAN(natom),srad(natom),icharge(natom)
      dimension rduser(110)
      character*256 cradts,cradi
      character*8 atsy
      character*5 cradt
c
c   Some values of COSMO radii (in Angstroms) taken from Andreas
c   Klamt notes. This table will have to be completed, I guess.
c
      dimension rad(110)
      save rad
      data rad/
     $ 1.3d0,0.0d0,                                              ! H -He
     $ 0.0d0,0.0d0,2.0475d0, 2.0d0, 1.83d0, 1.72d0, 1.72d0,0.0d0,! Li-Ne
     $ 0.0d0,0.0d0,0.0d0, 2.457d0, 2.106d0, 2.16d0, 2.05d0,0.0d0,! Na-Ar
     $                   14*0.0d0, 2.223d0, 2.223d0,2.16d0,0.0d0,! K -Kr
     $                                    16*0.0d0, 2.32d0,0.0d0,! Rb-Xe
     $                                                  38*0.0d0,! Cs-U
     $                                                  18*0.0d0/
c
c  van der Waals radi for atoms
c  most of the values are from A.Bondi, J.Phys.Chem. 68, 441-451 (1964).
c
      dimension Bondirad(110)
      save Bondirad
      data Bondirad/
     $ 1.20d0,1.40d0,                                          ! H-He
     $ 1.82d0,1.45d0,1.80d0,1.70d0,1.55d0,1.52d0,1.47d0,1.54d0,! Li-Ne
     $ 2.27d0,1.73d0,2.30d0,2.10d0,1.80d0,1.80d0,1.75d0,1.88d0,! Na-Ar
     $ 2.75d0,0.00d0,                                          ! K-Ca
     $ 0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,                     ! Sc-Mn
     $ 0.00d0,0.00d0,1.63d0,1.40d0,1.39d0,                     ! Fe-Zn
     $               1.87d0,2.19d0,1.85d0,1.90d0,1.85d0,2.02d0,! Ga-Kr
     $ 0.00d0,0.00d0,                                          ! Rb-Sr
     $ 0.00d0,0.00d0,0.00d0,0.00d0,0.00d0,                     ! Y-Tc
     $ 0.00d0,0.00d0,1.63d0,1.72d0,1.58d0,                     ! Ru-Cd
     $               1.93d0,2.27d0,0.00d0,2.06d0,1.98d0,2.16d0,! In-Xe
     $ 0.00d0,0.00d0,                                          ! Cs-Ba
     $ 0.00d0,                                                 ! La
     $ 14*0.00d0,                                              ! Ce-Lu
     $        0.00d0,0.00d0,0.00d0,0.00d0,                     ! Hf-Re
     $ 0.00d0,0.00d0,1.75d0,1.66d0,1.55d0,                     ! Os-Hg
     $               1.96d0,2.02d0,0.00d0,0.00d0,0.00d0,0.00d0,! Tl-Rn
     $ 0.00d0,0.00d0,                                          ! Fr-Ra
     $ 0.00d0,                                                 ! Ac
     $ 0.00d0,0.00d0,1.86d0,                                   ! Th-U
     $                                               18*0.00d0/!
c
c    get conversion factor
c
      call getrval('angs',angs)
c
c    set radii
c
c   the atomic radii are chosen as follow:
c   1) if the user has defined an input value for the atom,
c      the used defined radius is used. Every atom still
c      not defined will be assigned the 'cosmo' radius
c
c   2) if no input values is present, the value indicated
c      by cradt ('cosmo' or 'bondi') is used.
c
c      'cosmo' means either the value stored in the rad array,
c       or  1.17*bondirad
c      'bondi' means the value stored in bondirad
c
c   3) if after 1) and  2)  the radius is still undefined, an
c      error message is generated
c
      call getchval('c_cradt',cradts)
      cradt=cradts(1:5)
      call zeroit(srad,natom)
c
c  in case of user defined radii, we have to parse the radii file
c
      if(cradt.eq.'user '.or.cradt.eq.'userb'.or.cradt.eq.'file ')then
        call getchval('c_cradi',cradi)
        call cosmo_parseradii(cradt,cradi,srad,natom)
      endif
c
      ih = 0      ! number of hydrogen atoms
      do i=1,natom
        iq = IAN(i)
        icharge(i) = iq
        if(iq.EQ.1) ih=ih+1
        ri=srad(i)
        if(ri.eq.0.0d0)then
          if(cradt.eq.'bondi'.or.cradt.eq.'userb')then
            ri=bondirad(iq)
          else
            ri=rad(iq)
            if(ri.eq.0.0d0)ri=1.17d0*bondirad(iq)
          endif
        endif
        if(ri.eq.0.0d0)then
           write(6,*)
           call nerror(1,'cosmo_radii',
     $       'radius for atom '//atsy//' not defined',0,0)
        endif
        srad(i)=ri*angs
      enddo
c
c   compute maximum number of surface elements (ih= number of hydrogens)
c
      call getival('c_nspa',nspa)
      call getival('c_nsph',nsph)
      maxnps=(natom-ih)*nspa+ih*nsph
      if(maxnps.lt.500)maxnps=int(1.2d0*float(maxnps))
      call setival('c_maxnps',maxnps)
c
      return
      end
c======================================================================
      subroutine cosmo_ediel(qcos,phi,fepsi,nps,ediel)
      implicit real*8 (a-h,o-z)
c
c   this routine computes the COSMO interaction energy
c
c    qcos    unscaled screening charges (nps)
c    phi     surface potential (nps)
c    fepsi   scaling factor
c    nps     number of surface segments
c    ediel   in output will contain the interaction energy
c
c    ediel  = 0.5 * fepsi * scalar(qcos,phi)
c
      dimension qcos(nps),phi(nps)
      parameter (zero=0.0d0,half=0.5d0)
c
      scalp=zero
      do i=1,nps
         scalp=scalp+qcos(i)*phi(i)
      enddo
      ediel=half*fepsi*scalp
c
      return
      end
c======================================================================
      subroutine cosmo_print(iout,etot,cosmotime,cosmoelap)
      implicit real*8 (a-h,o-z)
c
c     prints the COSMO results on channel iout
c
      call getrval('c_de',de)
      call getrval('c_dq',dq)
      call getrval('c_qsum',qsum)
      call getrval('c_ediel',ediel)
c
      write(iout,100)qsum,dq,qsum+dq,etot,etot+de,ediel,ediel+de
c
      write(iout,101)cosmotime/60.0d0,cosmoelap/60.0d0
      return
c
 100  format(/,
     $'====================== COSMO RESULTS =========================',
     $     /,/,1x,'  Screening charge:',/,
     $         1x,'       Cosmo      =',f12.6,/,
     $         1x,'       Correction =',f12.6,/,
     $         1x,'       Total      =',f12.6,/,/,
     $         1x,'  Total energy                 =',f20.9,' Eh',/,
     $         1x,'  Total energy + OC corr.      =',f20.9,/,
     $         1x,'  Dielectric energy            =',f20.9,/,
     $         1x,'  Dielectric energy + OC corr. =',f20.9)
 101  format(/,'Master CPU time for COSMO = ',f7.2,' Elapsed = ',f7.2,
     $' min')
c
      end
c======================================================================
      subroutine cosmo_store(xyz,cosurf,ar,qcos,iatsp,icharge,
     $                       natom,nps)
      implicit real*8 (a-h,o-z)
c
c  store the COSMO surface data on a file
c
      dimension xyz(3,natom),cosurf(3,nps),ar(nps)
      dimension qcos(nps),iatsp(nps),icharge(natom)
      character *256 cdata
      logical busy
c
c  get name of cosmo data file
c
      call getchval('c_data',cdata)
c
c   look for a free I/O channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(1,'cosmo_store','cannot find a free channel',0,0)
 10   continue
c
c   open the cosmo data file
c
      open(unit=ich,file=cdata(1:len_trim(cdata)),status='unknown',
     $     form='unformatted')
c
c   write data
c
      rewind(ich)
      write(ich)xyz      ! nuclear coordinates (real atoms only)
      write(ich)icharge  ! nuclear charges (real atoms only)
      write(ich)cosurf   ! coordinates of surface segments
      write(ich)ar       ! area of surface segments
      write(ich)iatsp    ! mapping surface segment -> atoms
      write(ich)qcos     ! surface charges (unscaled and uncorrected)
c
c    close file
c
      close(ich)
      return
      end
c======================================================================
      subroutine cosmo_rdata(xyz,cosurf,ar,qcos,iatsp,icharge,
     $                       natom,nps)
      implicit real*8 (a-h,o-z)
c
c  retrieve the COSMO surface data from a file
c
      dimension xyz(3,natom),cosurf(3,nps),ar(nps)
      dimension qcos(nps),iatsp(nps),icharge(natom)
      character *256 cdata
      logical busy
c
c  get name of cosmo data file
c
      call getchval('c_data',cdata)
c
c   look for a free I/O channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(1,'cosmo_rdata','cannot find a free channel',0,0)
 10   continue
c
c   open the cosmo data file
c
      open(unit=ich,file=cdata(1:len_trim(cdata)),status='old',
     $     form='unformatted')
c
c   read data
c
      rewind(ich)
      read(ich)xyz      ! nuclear coordinates (real atoms only)
      read(ich)icharge  ! nuclear charges (real atoms only)
      read(ich)cosurf   ! coordinates of surface segments
      read(ich)ar       ! area of surface segments
      read(ich)iatsp    ! mapping surface segment -> atoms
      read(ich)qcos     ! surface charges (unscaled and uncorrected)
c
c    close file
c
      close(ich)
      return
      end
c======================================================================
      subroutine cosmo_surf(natom,maxnps,nspa,nsph,nppa,xyz,nuc,
     $                      srad,cosurf,ar,iatsp,dirsm,dirsmh,dirvec,
     $                      dirtm,xyzpert,radtmp,tm,nar,nsetf,nset)

      use memory

      implicit real*8 (a-h,o-z)
c
c   this routine drives the construction of the COSMO surface
c
      dimension xyz(3,natom),xyzpert(3,natom),srad(natom),nuc(natom),
     $          cosurf(3,2*maxnps),ar(2*maxnps),iatsp(2*maxnps),
     $          dirsm(3*nspa),dirsmh(3*nsph),dirvec(3*nppa),
     $          dirtm(3*nppa),radtmp(natom),tm(3,3,natom),nar(2*maxnps),
     $          nsetf(2*maxnps),nset(natom,nppa)
      character*100 amat,cerm
c
c     common /big/bl(30000)
c
c
c  fill the direction vectors
c
      ierr=0
      call dvfill(nspa,dirsm,ierr,cerm)
      if(ierr.ne.0)then
        call nerror(1,'cosmo_surf',cerm,0,0)
      endif
      call dvfill(nsph,dirsmh,ierr,cerm)
      if(ierr.ne.0)then
        call nerror(2,'cosmo_surf',cerm,0,0)
      endif
c
c  fill the basis grid vectors
c
      if(nppa.eq.1082)then
        call fill_grd1082(dirvec)
      else
        call dvfill(nppa,dirvec,ierr,cerm)
        if(ierr.ne.0)then
          call nerror(3,'cosmo_surf',cerm,0,0)
        endif
      endif
c
c   copy coordinates to a temporary array (the COSMO surface
c   construction uses distorted coordinates)
c
      do i=1,natom
        xyzpert(1,i)=xyz(1,i)
        xyzpert(2,i)=xyz(2,i)
        xyzpert(3,i)=xyz(3,i)
      enddo
c
c  setup temporary array for radii
c
      call getrval('c_rsolv',rsolv)
      do i=1,natom
        radtmp(i)=srad(i)+rsolv
      enddo
c
c   surface construction
c
      call mmark
      call getmem(3*2*maxnps,ixsp)
      call getmem(4*6*natom,isude)
      call getint(nppa+1,idin)
      call getint(4*natom,inipa)
      call getint(4*12*natom,ilipa)
      call getint(4*6*natom,iisude)
      call getint(3*natom,inn)
c
      call getrval('c_disex',disex)
      call getrval('c_routf',routf)
      call getrval('c_phsran',phsran)
      call getrval('c_ampran',ampran)
      call getival('c_lcavit',lcavity)
c
      ierr=0
      call consts(xyzpert,nuc,radtmp,cosurf,iatsp,dirsm,
     $            dirsmh,dirvec,nar,ar,dirtm,nsetf,nset,
     $            natom,maxnps,nspa,nsph,nppa,disex,
     $            disex2,rsolv,routf,nps,npspher,
     $            area,volume,phsran,ampran,lcavity,tm,
     $            ierr,cerm,bl(idin),bl(ixsp),bl(isude),bl(inipa),
     $            bl(ilipa),bl(iisude),bl(inn))
      if(ierr.ne.0)then
        call nerror(4,'cosmo_surf',cerm,0,0)
      endif
      call retmark
c
c  store some values into the depository
c
      call setival('c_nps',nps)  ! number of surface segments
      call setival('c_npsphe',npspher) ! outer surface segments
      npsd=nps+npspher
      call setival('c_npsd',npsd)
      call setrval('c_disex2',disex2)
      call setrval('c_area',area) ! surface area
      call setrval('c_volume',volume) ! surface volume
c
c  print out surface information
c
      call getival('iout',iout)
      write(iout,'(6x,a,i10)')   'No. of surface segments =',nps
      write(iout,'(6x,a,f15.4)') 'Surface area (au)       =',area
      write(iout,'(6x,a,f15.4)') 'Cavity volume (au)      =',volume
      write(iout,*)
c
      return
      end
c======================================================================
      subroutine cosmo_dump(filn,l,a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      character*120 filn
      logical busy
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(2,'cosmo_dump','cannot find a free channel',0,0)
 10   continue
      open(unit=ich,file=filn(1:l),status='unknown',form='unformatted')
      write(ich)a
      close(unit=ich,status='keep')
      return
      end
c======================================================================
      subroutine cosmo_rest(filn,l,a,n)
      implicit real*8 (a-h,o-z)
      dimension a(n)
      character*120 filn
      logical busy
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(2,'cosmo_rest','cannot find a free channel',0,0)
 10   continue
      open(unit=ich,file=filn(1:l),status='old',form='unformatted')
      read(ich)a
      close(unit=ich,status='delete')
      return
      end
c======================================================================
      subroutine cosmo_del(name)
      implicit real*8 (a-h,o-z)
      character*8 name
      character*256 filn
      logical busy,isthere
      call getchval(name,filn)
      inquire(file=filn(1:len_trim(filn)),exist=isthere)
      if(.not.isthere)return
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(2,'cosmo_del','cannot find a free channel',0,0)
 10   continue
      open(unit=ich,file=filn(1:len_trim(filn)),status='old')
      close(unit=ich,status='delete')
      return
      end
c======================================================================
      subroutine cosmo_h0(hmat,vi,inx,basdat,cosurf,qcos,
     $                    fepsi,nps,ncf,ncs,ntri)
      implicit real*8 (a-h,o-z)
c
c   this subroutine computes the COSMO contribution to the one electron
c   Hamiltonian
c
c   the electron-surface matrix elements are either read from file
c   or directly computed, according to the value of the COSMO parameter
c   c_onel
c
c   hmat            one-electron Hamiltonian (ntri)
c   vi              area for electron-surface matrix elements (ntri)
c   qcos            unscaled screening charges (nps)
c   fepsi           scaling factor
c   nps             number of surface segments
c   ntri            dimension of triangular arrays
c
      dimension hmat(ntri),vi(ntri)
      dimension qcos(nps),cosurf(3,nps)
      dimension inx(12,*),basdat(13,*)
      character*256 conel
      logical busy,conven
      parameter (zero=0.0d0)
c
      call getchval('c_onel',conel)
      conven=conel.ne.'direct'
c
c   look for a free I/O channel
c
      if(conven)then
        do i=20,99
          inquire(i,opened=busy)
          if(.not.busy)then
             ich=i
             goto 10
           endif
        enddo
        call nerror(2,'cosmo_h0','cannot find a free channel',0,0)
 10     continue
c
c   open the direct access file containing the surface matix elements
c
        open(unit=ich,file=conel(1:len_trim(conel)),status='old',
     $       form='unformatted',access='direct',recl=ntri*8)
      endif
c
c   loop over surface segment
c
      call zeroit(hmat,ntri)
      do ip=1,nps
c
c   read or compute matrix elements for surface segment ip
c
        if(conven)then
          read(unit=ich,rec=ip)vi
        else
          call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
        endif
c
c  add contribution of surface segment to one-electron Hamiltonian
c
c       write(6,*)'segment',ip,' fepsi ',fepsi,' qcos ',qcos(ip)
c       write(6,*)'Vi matrix'
c       write(6,*)vi
        fq=-fepsi*qcos(ip)
        do i=1,ntri
          hmat(i)=hmat(i)+fq*vi(i)
        enddo
c       call daxpy(ntri,fq,vi,1,hmat,1)
c       write(6,*)'new H0 matrix'
c       write(6,*)hmat
      enddo
c
c  close the cosmo matrix elements file
c
      if(conven)close(ich,status='keep')
c
      return
      end
c======================================================================
      subroutine cosmo_pot(den,inx,basdat,phin,phi,cosurf,vi,
     $                     nps,ncf,ncs,ntri)
      implicit real*8 (a-h,o-z)
c
c   this subroutine computes the COSMO surface potential
c
c   for each segment of the surface, this subroutine reads or
c   computes the corresponding matrix elements, takes the trace
c   with the current density matrix, adds the nuclear term and
c   store the result in the array phi:
c
c   phi_i = tr[D*v_i] + phin_i
c
c   where:
c          phi_i     is the potential at surface segment i
c          phin_i    is the nuclear part of potential at
c                    surface segment i
c          D         is the density matrix
c          v_i       are the electron-surface repulsion matrix
c                    elements of surface segment i
c
c  the matrix elements are treated according to the content of
c  the COSMO parameter c_onel, which is gathered from the depository
c  system. if c_onel='direct' the matrix elements are to be computed,
c  otherwise c_onel will contain the name of the file that stores
c  the one electron integrals
c
c   den           density matrix (ntri)
c   inx           contraction info
c   basdat        basis set data
c   phin          nuclear contribution to surface potentials (nps)
c   phi           in output will contain surface potentials (nps)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   vi            area for matrix elements for a surface segment (ntri)
c   nps           number of surface segments
c   ncf           number of basis function
c   ncs           number of contracted shells
c   ntri          ncf*(ncf+1)/2
c
      dimension den(ntri),vi(ntri),cosurf(3,nps)
      dimension phin(nps),phi(nps)
      dimension inx(12,*),basdat(13,*)
      character*256 conel
      logical busy,conven
      parameter (zero=0.0d0)
c
      call getchval('c_onel',conel)
c
c   set the type of calculation (direct or conventional) according
c   to the value of the variable conel.
c
      conven=conel.ne.'direct'
c
c   for conventional calculation, open the matrix element file
c
      if(conven)then
c
c   look for a free I/O channel
c
        do i=20,99
          inquire(i,opened=busy)
          if(.not.busy)then
             ich=i
             goto 10
           endif
        enddo
        call nerror(2,'cosmo_pot','cannot find a free channel',0,0)
 10     continue
c
c   open the direct access file containing the surface matix elements
c
        open(unit=ich,file=conel(1:len_trim(conel)),status='old',
     $       form='unformatted',access='direct',recl=ntri*8)
      endif
c
c  loop over surface segments
c
      do ip=1,nps
c
c   read or compute matrix elements for surface segment ip
c
        if(conven)then
          read(unit=ich,rec=ip)vi
        else
          call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
        endif
c
c   take trace with density
c
        call spur(den,vi,ncf,trace)
c
c   add nuclear term and fill phi array
c
        phi(ip)=-trace+phin(ip)
      enddo
c
c  close the cosmo matrix elements file
c
      if(conven)close(ich,status='keep')
c
      return
      end
c======================================================================
      subroutine cosmo_potn(xyz,iz,phin,cosurf,natom,nps)
      implicit real*8 (a-h,o-z)
c
c   this subroutine computes the nuclear part of the COSMO
c   surface potential
c
c   phin_i =  Sum_j Z_j / |R_j-s_i|
c
c   where:
c          phin_i    is the nuclear part of potential at
c                    surface segment i
c          j         runs over the atoms (real atoms only)
c          Z_j       atomic charge
c          R_j       atomic coordinates
c          s_i       coordinates of surface segment i
c
c   xyz           atomic coordinates (real atoms only)
c   iz            atomic charges (real atoms only)
c   phin          in output will contain surface potentials (nps)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   natom         number of real atoms
c   nps           number of surface segments
c
      dimension cosurf(3,nps)
      dimension xyz(3,natom),iz(natom),phin(nps)
      parameter (zero=0.0d0)
c
c  loop over surface segments
c
      do ip=1,nps
c
c   loop over atoms
c
        sumj=zero
        do j=1,natom
c
c   compute distance between current atom and current surface segment
c
          dx=xyz(1,j)-cosurf(1,ip)
          dy=xyz(2,j)-cosurf(2,ip)
          dz=xyz(3,j)-cosurf(3,ip)
          d=sqrt(dx*dx+dy*dy+dz*dz)
c
c   add nuclear term
c
          sumj=sumj+float(iz(j))/d
        enddo
c
c   fill phin array
c
        phin(ip)=sumj
      enddo
c
      return
      end
c======================================================================
      subroutine cosmo_eint(den,xyz,iz,cosurf,h0cos,qcos,fepsi,eint,
     $                      natom,nps,ncf,ntri)
      implicit real*8 (a-h,o-z)
c
c   this subroutine computes the cosmo interaction energy.
c   this is just for testing. The actual energy correction is better
c   done using the dielectric energy computed by cosmo_ediel
c
c
c   e_int = tr[D*h0_cos] + fepsi Sum_i q_i * [Sum_j Z_j / |R_j-s_i|]
c
c   where:
c          D         is the density matrix
c          h0_cos    is the cosmo contribution to the one electron
c                    hamiltonian
c          fepsi     is the function of the dielectric costant
c          i         runs over surface segments
c          q_i       is the charge at surface segment i
c          j         runs over the atoms (real atoms only)
c          Z_j       atomic charge
c          R_j       atomic coordinates
c          s_i       coordinates of surface segment i
c
c   den           density matrix (ntri)
c   xyz           atomic coordinates (real atoms only)
c   iz            atomic charges (real atoms only)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   h0cos         COSMO contribution to one electron Hamiltonian (ntri)
c   qcos          surface charges (nps)
c   fepsi         fepsi
c   eint          in exit will contain the interaction energy
c   natom         number of atoms
c   nps           number of surface segments
c   ncf           number of basis function
c   ntri          ncf*(ncf+1)/2
c
      dimension den(ntri),h0cos(ntri),cosurf(3,nps)
      dimension xyz(3,natom),iz(natom),qcos(nps)
      parameter (zero=0.0d0)
c
c  trace with density matrix
c
      call spur(den,h0cos,ncf,eint)
c
c  loop over surface segments
c
      sumi=zero
      do ip=1,nps
c
c   loop over atoms
c
        sumj=zero
        do j=1,natom
c
c   compute distance between current atom and current surface segment
c
          dx=xyz(1,j)-cosurf(1,ip)
          dy=xyz(2,j)-cosurf(2,ip)
          dz=xyz(3,j)-cosurf(3,ip)
          d=sqrt(dx*dx+dy*dy+dz*dz)
c
c   add nuclear term
c
          sumj=sumj+float(iz(j))/d
        enddo
c
c   multiply by surface charge
c

        sumi=sumi+qcos(ip)*sumj
      enddo
c
c   scale by fepsi, and accumulate into eint
c
      eint=eint+fepsi*sumi
c
      return
      end
c======================================================================
      subroutine cosmo_surfrep(vi,cosurf,inx,basdat,
     $                         nps,ncs,ncf,ntri)
      implicit real*8 (a-h,o-z)
c
c   this subroutine drives the calculation of the matrix elements
c   of electron-surface repulsion.
c
c   the matrix elements for each segment are computed and
c   stored in a file to be used later.
c   the calculation is based on the one-electron integral
c   code (subroutines inton, onel, incons and nucpot). The
c   integrals are ultimately computed as nuclear potential
c   integrals, calling subroutine nucpot, except that the
c   coordinates of the surface segments are used in place of the
c   nclear coordinates
c
c   for the moment the external loop is over surface segments
c   I think a faster solution might be to have the loop over
c   surface segments as the inner loop, so the loops over
c   the contraction are executed only once. Although at the
c   moment it is not clear to me if it is possible to do this
c   without using far too much memory, i.e., without storing
c   up to 63504*nps integrals over primitives (nps is the
c   number of surface segments, easily a few hundreds)
c
c   vi            area for matrix elements for a surface segment(ntri)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   inx           contraction info
c   basdat        basis set data
c   natom         number of atoms
c   nps           number of surface segments
c   ncs           number of contracted shells
c   ncf           number of basis function
c   ntri          ncf*(ncf+1)/2
c
      dimension vi(ntri),cosurf(3,nps)
      dimension inx(12,*),basdat(13,*)
      character*256 conel
      logical busy
c
      call getchval('c_onel',conel)
c
c  if conel='direct' the integrals will be recomputed each time
c  they are needed, thus nothing has to be done here
c
      if(conel.eq.'direct')return
c
c
c   look for a free I/O channel
c
      do i=20,99
        inquire(i,opened=busy)
        if(.not.busy)then
           ich=i
           goto 10
         endif
      enddo
      call nerror(2,'cosmo_surfrep','cannot find a free channel',0,0)
 10   continue
c
c   open a direct access file
c
      open(unit=ich,file=conel(1:len_trim(conel)),status='unknown',
     $     form='unformatted',access='direct',recl=ntri*8)
c
c  loop over surface segments
c
      do ip=1,nps
c
c   compute matrix elements for surface segment ip
c
        call cosmo_inton(ip,vi,inx,basdat,cosurf,ncs)
c       write(6,*)'Segment ',ip,' coordinates',cosurf(1,ip),
c    $  cosurf(2,ip),cosurf(3,ip)
c       write(6,*)'vi matrix'
c       write(6,*)vi
c       write(6,*)
c       call f_lush(6)
c
c   save the current batch of matrix elements
c
        write(unit=ich,rec=ip)vi
      enddo
c
c  close the cosmo matrix elements file
c
      close(ich,status='keep')
c
      return
      end
c======================================================================
      subroutine cosmo_pot_forc(den,dcos,qcos,cosurf,inx,basdat,iatsp,
     $                         nps,ncs,ncf,ntri)
      implicit real*8 (a-h,o-z)
c
c   this subroutine drives the calculation of the derivatives
c   of the cosmo potential. Only the electron-surface contribution
c   to the gradient is computed here; the nuclear-surface term is
c   computed by subroutine cavnucgrd.
c
c   the calculation is based on the one-electron gradient
c   code (subroutines intof7, onef7, inconf7 and nucpof). The
c   integrals are ultimately computed as nuclear potential
c   integrals derivative, calling onef7, except that the
c   coordinates and charge of the surface segments are used in
c   place of nuclear coordinates and charges
c
c   den           density matrix in triangular form (ntri)
c   dcos          gradient array (3,*)
c   qcos          surface charges (nps)
c   cosurf        cartesian coordinates of surface segments (3*nps)
c   inx           contraction info
c   basdat        basis set data
c   iatsp         mapping surface segment -> atoms (nps)
c                  NOTE: COSMO does not use ghost atoms for surface
c                        construction, thus the mapping of iatsp
c                        works only if the centers of the calculation
c                        are ordered with the real atoms first,
c                        followed by the ghost atoms.
c   nps           number of surface segments
c   ncs           number of contracted shells
c   ncf           number of basis function
c   ntri          ncf*(ncf+1)/2
c
      dimension den(ntri),cosurf(3,nps),qcos(nps),iatsp(nps)
      dimension inx(12,*),basdat(13,*),dcos(3,*)
c
c  loop over surface segments
c
      do ip=1,nps
c
c   compute matrix elements for surface segment ip
c
        call cosmo_intof7(ip,inx,basdat,cosurf(1,ip),iatsp(ip),
     $                    ncs,ncf,ntri,den,dcos,qcos(ip))
      enddo
c
      return
      end
c=======================================================================
c
      subroutine cosmo_inton (ipt,oneint,inx,basdat,cosurf,ncs)

      use memory

      implicit real*8 (a-h,o-z)
c
c     ipt is the surface point for which the matrix elements are
c     to be computed
c     oneint (out) is the one-electron integral matrix, stored as the
c     upper triangular in Fortran, i.e. ((h(i,j),i=1,j),j=1,n)
c     inx (in) holds the contraction info
c     basdat hold the basis set data
c     cosurf holds the crtesian coordinates of the surface points
c     ncs is the total number of contracted shells
c
      dimension inx(12,*)
      dimension oneint(*),basdat(13,*),cosurf(3,*)
c     common /big/bl(1)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /onethr/onethre
c     reserve memory for the array used to store the general contr.
c     integrals. This is quite big because we may have up to i functions
c     (28 components), i.e. 28**2=784, and a maximum of 9 general
c     contractions for each function, i.e.  81*784=63504
      call getmem(63504,is)
      iprint=igetival('iprn')
c     iout=igetival('iout')
c     if (iprint.eq.1) write (iout,80)'el-surf',ipt
   80 format (//,1x,a7,i6,3x,2a1/)
c     total number of 1-el. integrals
      iza=0
c     ifu=0
c  setup threshold for electron-surface potential
c  it is 1/100th of the main int. threshold
      onethre=0.01d0*rgetrval('ithr')
c
c     cycle through the contracted shells and calculate the integrals
      do 60 ics=1,ncs
c       number of gen. contractions
        ngci=inx(4,ics)
ckwol atoms=centers of contracted shells
        iatom=inx(2,ics)
c
         jfu=0
         len1=inx(3,ics)
         do 50 jcs=1,ics
           ngcj=inx(4,jcs)
           jatom=inx(2,jcs)
            len2=inx(3,jcs)
            len=len1*len2*(ngci+1)*(ngcj+1)
            call cosmo_onel(ics,jcs,len2,len1,ngcj,ngci,basdat,
     $                      cosurf,bl(is),inx,ipt)
          iij=-1
          ifu=inx(11,ics)-len1
          do 45 igc=0,ngci
            ifu=ifu+len1
            jfu=inx(11,jcs)-len2
            do 45 jgc=0,ngcj
              jfu=jfu+len2
              iff=ifu
              do 40 i1=1,len1
                iff=iff+1
                jff=jfu
           ii=iff*(iff-1)/2
           do 40 j1=1,len2
             jff=jff+1
                  iij=iij+1
                  if (jff.le.iff) then
                    ij=ii+jff
                    oneint(ij)=bl(is+iij)
                    iza=iza+1
                  end if
 40           continue
 45       continue
c.....................
 50   continue
 60   continue
      n1=1
      if(iza.gt.0) then
c        call wri (oneint,iza,n1,0,name)
        if (iprint.eq.1) then
          iout = igetival('iout')
          write(iout,*)' Matrix from cosmo_inton'
          ncf=igetival('ncf')
          ii=0
          do i=1,ncf
            write(iout,90) (oneint(k),k=ii+1,ii+i)
            ii=ii+i
   90       format (2x,5f14.9)
            write(iout,90)
          end do
        end if
      end if
c     return memory
      call retmem(1)
c
c
      end
c
c=======================================================================
c
      subroutine cosmo_onel (ics,jcs,jlen,ilen,ngcj,ngci,basdat,
     1  cosurf,ss,inx,nn)
      implicit real*8 (a-h,o-z)
      dimension basdat(13,*),cosurf(3,*)
      dimension xa(3), xb(3), inx(12,*), s(784),xint(784), p(3)
      dimension ss(jlen,ilen,0:ngcj,0:ngci)
      common /tape/ inp,inp2,iout
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /onethr/onethre
      data sqtw/0.2886 7513 4594 8128 8d0/
c
c this is 1/sqrt(12)
c
c     data twopi/0.6366197723675d0/
c
c this subroutine calculates the COSMO one-electron integrals for a
c contracted pair of shells,including all general contractions
c ics,jcs are the indices of contracted shells
c jlen and ilen are the shell sizes (3 for p...) for the second
c   and first contractions - they are deliberately interchanged
c ngcj and ngci are the general contraction lengths, starting at 0
c basdat and inx contain the shell info,
c cosurf contains the cartesian coordinates of the surface segments
c ss holds the integrals, for all general contractions over the shells
c format is ss(jlen,ilen,0:ngci,0:ngcj)
c ics and jcs, in the order s(j,i,ngcj,ngci)
c inx is the contraction info
c nn  is the surface segment under consideration
c
      twopi=two/pi
      ityp=inx(12,ics)
      jtyp=inx(12,jcs)
      len=ilen*jlen
      ilen1=ilen
      jlen1=jlen
      if (ityp.eq.4) ilen1=6
      if (jtyp.eq.4) jlen1=6
c  f functions - 10 long before transformation
      if(ityp.eq.6) ilen1=10
      if(jtyp.eq.6) jlen1=10
c  g functions - 15 long before transformation
      if(ityp.eq.11) ilen1=15
      if(jtyp.eq.11) jlen1=15
c  h functions - 21 long before transformation
      if(ityp.eq.12) ilen1=21
      if(jtyp.eq.12) jlen1=21
c  i functions - 15 long before transformation
      if(ityp.eq.13) ilen1=28
      if(jtyp.eq.13) jlen1=28
c
c      len1=ilen1*jlen1
c
c for d functions, 6 spaces are needed prior to the transformation
c
c     zero out the buffer
      call zeroit(ss,len*(ngci+1)*(ngcj+1))
      ia=inx(1,ics)+1
      ie=inx(5,ics)
      ja=inx(1,jcs)+1
      je=inx(5,jcs)
      do 60 i=ia,ie
         a=basdat(1,i)
         sqa=sqrt(a)
         do 20 l=1,3
           xa(l)=basdat(l+10,i)
   20    continue
         do 50 j=ja,je
            b=basdat(1,j)
            do 24 l=1,3
               xb(l)=basdat(l+10,j)
  24        continue
            sqb=sqrt(b)
            apb=a+b
            r=zero
            e=a*b/apb
            do 30 l=1,3
c              xb(l)=bl(jj+l+4)
               p(l)=(a*xa(l)+b*xb(l))/apb
               r=r+(xa(l)-xb(l))**2
   30       continue
            s0=sqrt(twopi*e/(sqa*sqb))**3*exp(-e*r)
c
c     write(iout,1000) ics,jcs,ityp,jtyp,i,j,ilen1,jlen1,csa,cpa,csb,cpb
c
      ityp1=ityp
      jtyp1=jtyp
c correspondance between ityp and ityp1:
c ityp  1(s) 2(p) 3(l) 4(d) 5(d6) 6(f) 7(f10) 8(g15) 9(h21) 10(i28) 11(g) 12(h) 13(i)
c ityp1 1    2    3    4    4     5    5      6      7       8        6     7    8
      if(ityp.ge.5) ityp1=ityp-1
      if(jtyp.ge.5) jtyp1=jtyp-1
      if(ityp.ge.7) ityp1=ityp-2
      if(jtyp.ge.7) jtyp1=jtyp-2
        if(ityp.eq.11) ityp1=6
        if(jtyp.eq.11) jtyp1=6
        if(ityp.eq.12) ityp1=7
        if(jtyp.eq.12) jtyp1=7
        if(ityp.eq.13) ityp1=8
        if(jtyp.eq.13) jtyp1=8
c
c     *** no special type is needed for d6 at this point
c
c Provide an estimate for the integrals
          estim=s0*max(sqa,sqb)
c         estim is the inverse square root of the larger exponent
c         times the primitive overlap s0
          if(estim.lt.onethre) go to 50
c  jump out of the loop (CYCLE)
            call cosmo_incons(ityp1,jtyp1,nn,a,b,s0,xa,xb,
     1      s,cosurf)
c
c     transform the first subscript if needed
c   d
        if (ityp.eq.4) then
           call dtran1a(xint,s,jlen1)
        end if
c  f
        if (ityp.eq.6) then
          call ftran1a(xint,s,jlen1)
        end if
c transform the 15 Cartesian g functions to 9 spherical harmonics components
        if(ityp.eq.11) then
            call gtran1a(xint,s,jlen1)
          end if
c transform the 21 Cartesian h functions to 11 spherical harmonics components
        if(ityp.eq.12) then
            call htran1a(xint,s,jlen1)
          end if
c transform the 28 Cartesian i functions to 13 spherical harmonics components
        if(ityp.eq.13) then
c	    call itran1a(xint,s,jlen1)
          end if
c
c  transform the second subscript
c  d
        if(jtyp.eq.4) then
           call dtran2a(xint,s,ilen)
        end if
c  f
        if (jtyp.eq.6) then
          call ftran2a (xint,s,ilen)
        end if
c g
          if(jtyp.eq.11) then
            call gtran2a(xint,s,ilen)
          end if
c  h
          if(jtyp.eq.12) then
            call htran2a(xint,s,ilen)
          end if
c
          if(jtyp.eq.13) then
c	    call itran2a(xint,s,ilen)
          end if
c
c  *** add the primitive integrals to the contracted ones
        do 45 igc=0,ngci
          csa=basdat(igc+2,i)
c         cpa is needed only for the l type
          cpa=basdat(3,i)
          do 45 jgc=0,ngcj
              csb=basdat(jgc+2,j)
              cpb=basdat(3,j)
              ij=0
              do 40 i1=1,ilen
                coefi=csa
                if (ityp.eq.3.and.i1.gt.1) coefi=cpa
                do 40 j1=1,jlen
                  ij=ij+1
                  coefj=csb
                  if (jtyp.eq.3.and.j1.gt.1) coefj=cpb
                  ss(j1,i1,jgc,igc)=ss(j1,i1,jgc,igc)+s(ij)*coefi*coefj
   40         continue
   45     continue
   50   continue
   60 continue
c
      end
c
c=======================================================================
c
      subroutine cosmo_incons (ityp,jtyp,nna,a,b,s0,xa,xb,xint,cosurf)
c    this subroutine constructs a set of COSMO one-electron integrals
c    over two primitive shells
c
c    the integrals are computed as the nuclear-electron repulsion
c    integrals (by subroutine nucpot), except that the coordinates
c    of a surface segment are used instead of the atomic coordinates
c
c   ARGUMENTS:
c
c  ityp,jtyp: the types of the shells, s=1,p=2,l=3,d=4,f=5,g=6,h=7,i=8
c    note that only 6-component d and 10-component f are used here
c    transformation to d5 and f7 is later
c  a and b are the Gaussian exponents
c  s0 is the overlap integral of the spherical parts of the Gaussians
c  xa(3) and xb(3) are the center coordinates of the Gaussians
c  cosurf stores the coordinates of the surface segments
c  nna is the current surface segment
c
c  xint(1:lenj,1:leni) gives the integrals; leni and lenj are the
c   sizes (lengths) of the shells, e.g. 3 for P, 4 for L, 6 for D6
c

      use memory

      implicit real*8 (a-h,o-z)
      dimension cosurf(3,*)
      dimension xint(784),xa(3), xb(3), yint(784), len(8), la(3), lb(3),
     1ll(3)
      dimension mulx(88), muly(88), mulz(88), nfu(9), xc(3)
c     common /big/bl(1)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /tape/ inp,inp2,iout
      data len/1,3,4,6,10,15,21,28/
      data nfu/0,1,4,8,14,24,39,60,88/
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
c this threshold was not used at all fetching it took a lot of time
c  get the main integral threshold
c     thre=rgetrval('ithr')
c  the nuclear potential threshold is 10 times less
c     thre=0.1d0*thre
      ln=len(ityp)*len(jtyp)
      do 10 i=1,ln
         yint(i)=zero
   10 continue
      do 20 i=1,3
         la(i)=0
         lb(i)=0
         ll(i)=0
   20 continue
c
      xc(1)=cosurf(1,nna)
      xc(2)=cosurf(2,nna)
      xc(3)=cosurf(3,nna)
        call nucpot (ityp,jtyp,a,b,s0,xa(1),xa(2),xa(3),
     1  xb(1),xb(2),xb(3),xc(1),xc(2),xc(3),la,lb,0,xint)
c
c     write(iout,1100) key,k1,k2,a,b
c
c
c     write(iout,1200) (xint(i),i=1,ln)
c
      return
c
      end
c======================================================================
      subroutine cosmo_intof7(iss,inx,datbas,xss,nra,
     $                        ncs,ncf,ntri,dn,atfor,qss)
      implicit real*8 (a-h,o-z)
c---------------------------------------------------------------------
c
c This routine calculates Ist derivatives of the COSMO
c surface potential (electronic contribution only)
c the integrals are computed as electron-nuclear attaction
c integral derivatives, using subroutine onef7,
c excep that coordinates and charges of the surface segment
c are used instead of nuclear coordinates and charges.
c
c---------------------------------------------------------------------
c INPUT :
c iss             - current surface segment
c inx(12,ncs)     - basis set data
c datbas(13,nsh)  - basis set data
c xss(3)          - cartesian coordinates of current segment
c nra             - atom to which current segment belongs to
c ncs,ncf,ntri    - contr. shell, contr. functions & ntri=ncf*(ncf+1)/2
c dn(ntri)        - density matrix
c atfor(3,natoms) - forces on atoms - INPUT/OUTPUT
c qss             - charge of current surface segment
c---------------------------------------------------------------------
      dimension datbas(13,*)
      dimension inx(12,*)
      dimension dn(ntri)            ! density
      dimension atfor(3,*)          ! atfor(3,natoms)
clocal :
      dimension sa(3*784),sb(3*784) ! 3*(28x28) (i|i)
      dimension xss(3)
      data zero /0.d0/
c---------------------------------------------------------------------
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
            jfu=0
            do 50 j=1,i
               len2=inx(3,j)
               jat=inx(2,j)
               ngcj=inx(4,j)
               ngcjx=ngcj
               if(j.eq.i) ngcjx=igc
               len=len1*len2
               len3=len*3
               do 45 jgc=0,ngcjx
c no need to zero out - zeroed out in onef7
c                  do l=1,len3
c                     sa(l)=zero
c                     sb(l)=zero
c                  enddo
c
c compute integrals
c
                     call onef7(i,j,igc,jgc,datbas,sa,sb,inx,xss,qss)
c
                  iij=0
                  iff=ifu
                  do 40 i1=1,len1
                     iff=iff+1
                     jff=jfu
                     ii=iff*(iff-1)/2
                  do 40 j1=1,len2
                     jff=jff+1
                     iij=iij+1
                     if (jff.gt.iff) go to 40
                     ij=ii+jff
                     dij=dn(ij)
                     if(iff.ne.jff) dij=2.d0*dij
               do icr=1,3
                  iij_icr=iij+(icr-1)*len
                  atfor(icr,nra)=atfor(icr,nra)
     *                           -(sa(iij_icr)+sb(iij_icr))*dij*qss
                  atfor(icr,iat)=atfor(icr,iat)+sa(iij_icr)*dij*qss
                  atfor(icr,jat)=atfor(icr,jat)+sb(iij_icr)*dij*qss
               enddo
   40             continue
                  jfu=jfu+len2
c
   45          continue
   50       continue
            ifu=ifu+len1
   55    continue
   60 continue
c
c---------------------------------------------------------------------
c
      end
c =================================================================
      subroutine getatsym1(n,na,atsymb)

      use memory

      implicit real*8 (a-h,o-z)
c
c   this routine returns the atomic symbols of atom n in the
c   character*8 variable atsymb. Any number or special character
c   is blanked, and the output is converted to lower case, but
c   the routine does not check if the resulting symbol is a
c   valid  atomic symbol
c
c   n                 atom for which symbol has to be returned
c   na                number of atoms
c   atsymb  (char*8)  in exit will contain the requested symbol
c
      character *8 atsymb
c     common /big/bl(30000)
      character*25 special
      data special/'~!@#$%^&*+=<>?"0123456789'/
      character*26 alpbet,betalp
      data alpbet / 'abcdefghijklmnopqrstuvwxyz' /
      data betalp / 'ABCDEFGHIJKLMNOPQRSTUVWXYZ' /
c
      if(n.le.0.or.n.gt.na)then
         call nerror(1,'getatsym1','invalid atom number',n,na)
      endif
c
c  compute pointer for atom n
c
      call getival('inuc',inn)
      inn=inn+5*(n-1)
c   convert symbol from real to character, keep only
c   first two characters
      call blankit(atsymb,8)
      write(atsymb,'(a2,''      '')')bl(inn+4)
c   if second character is a number or a special symbol, blank it
      ii=index(special,atsymb(2:2))
      if(ii.gt.0)atsymb(2:2)=' '
c    convert to lower case
      ii=index(betalp,atsymb(2:2))
      if(ii.gt.0)atsymb(2:2)=alpbet(ii:ii)
      ii=index(betalp,atsymb(1:1))
      if(ii.gt.0)atsymb(1:1)=alpbet(ii:ii)
      return
      end
c =================================================================
      subroutine getatsymex(n,na,atsymb)

      use memory

      implicit real*8 (a-h,o-z)
c
c   this routine returns the extended atomic symbol of atom n
c   in the character*8 variable  atsymb.
c   the output is converted to lower case, but
c   the routine does not check if the resulting symbol is a
c   valid  atomic symbol
c
c   n                 atom for which symbol has to be returned
c   na                number of atoms
c   atsymb  (char*8)  in exit will contain the requested symbol
c
      character*8 atsymb,dummy
c     common /big/bl(30000)
c
      if(n.le.0.or.n.gt.na)then
         call nerror(1,'getatsymex','invalid atom number',n,na)
      endif
c
c  compute pointer for atom n
c
      call getival('inuc',inn)
      inn=inn+5*(n-1)
c
c  get symbol and lower case
c
      call blankit(atsymb,8)
      write(atsymb,'(a8)')bl(inn+4)
      call lowerca2(atsymb,8)
      return
      end
c =================================================================
      subroutine rem_num(xx,n)
      implicit real*8 (a-h,o-z)
      character*(*) xx
      character*10 special
      data special/'0123456789'/
c
      do i=1,n
 10     ii=index(special,xx(i:i))
        if(ii.gt.0)then
          do j=i+1,n
            jj=j-1
            xx(jj:jj)=xx(j:j)
          enddo
          xx(n:n)=' '
          goto 10
        endif
      enddo
      end

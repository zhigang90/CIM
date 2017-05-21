c======================================================================
c July 7,97 KW : the main gradient(forces) program
c August 2003  : revised for use with MP2 forces  (JB)
c April 2005   : revised for use with (parallel) FTC  (JB)
c February 2009: revised for simulated Atomic Force Microscopy  (JB)
c November 2011: revised for dispersion term in DFT  (JB)
c======================================================================
      subroutine forces

      use memory

      implicit real*8 (a-h,o-z)
      character jobname*256,datestr*24,cdum*20,wvfnc*20
      real*8 values(5)
      Logical rhf,noabc,tz,disp
      parameter (MaxEL=94,MaxC=5)  ! DFT dispersion
c     common /big/ bl(10000)
c     common /intbl/ifpp,inx(100)
c----------------------------------------------------------
c for several calls of the forces program in one run:
      data ncall /0/
      save ncall
c----------------------------------------------------------
      PARAMETER (IUnit=1)          ! unit number for checkpoint I/O
c----------------------------------------------------------
      iout=igetival('iout')
      inp=igetival('inp')
c----------------------------------------------------------
      ncall=ncall+1
c----------------------------------------------------------
c  read from <control> file
c  wavefunction type - if semiempirical, then exit
c  total number of atomic centers
c  multiplicity
c  number of alpha electrons
c  number of beta electrons
c  lowest eigenvalue of overlap matrix
c  dft flag
c  dispersion parameters
c  number of dummy atoms (if any)
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $     FORM='FORMATTED',STATUS='OLD')
      Call rdcntrl(IUnit,9,'$wavefunc',3,idum,rdum,wvfnc)
c
c -- check for semiempirical (FORCE command included erroneously)
      If(wvfnc(1:4).eq.'Semi'.or.wvfnc(1:4).eq.'semi') Then
        CLOSE (UNIT=IUnit,STATUS='KEEP')
c -- Skip over the FORCE card (It has been backspaced)
        read(inp,*)
        RETURN
      EndIf
c
      call rdcntrl(IUnit,7,'$natoms',1,NAtoms,dum,cdum)
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      Call rdcntrl(IUnit,5,'$xlow',2,idum,xlow,cdum)
      CALL RdDFT(IUnit,idft)
c
      call fdcntrl(IUnit,11,'$dispersion',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,910) IDisp
       READ(IUnit,911) edisp
       Call RdDISP(IUnit,noabc,tz,values,disp)
       If(.NOT.disp) call nerror(5,'FORCES Module',
     $        'Problems Reading Dispersion Parameters',0,0)
  910  Format(12X,I6)
  911  Format(7X,F20.10)
      EndIf
c
c -- dummy atoms
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) Ndum1,Ndum2
  900 Format(9X,I4,2X,I4)
       NAtom = NAtoms-Ndum1-Ndum2
      Else
       NAtom = NAtoms
      EndIf
c
      NAt3 = 3*NAtoms
c
c -- check there are electrons!
      If(NAlpha.EQ.0) Call nerror(1,'FORCE Module',
     $                     'There Are No Electrons!',0,0)
c
c -- set logical flag "rhf"
      rhf = (NBeta.EQ.0.AND.IMult.EQ.1)
c
c----------------------------------------------------------
c get original values of symmetry pair pointers
c
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
c----------------------------------------------------------
c Re-store commons big, intbl, and ener
c
      call readbl(ncall,ictr)
c
c output : ictr = new address in inx (texas95)
c----------------------------------------------------------
c  Read in the FORCE options
c
      call force_options(wvfnc,xlow,IdWt)
c
c -- write print and DFT weight derivative flag to <control> file
      call getival('nprint',IPRNT)
      if(idft.gt.0) call wrcntrl(IUnit,7,'$weight',1,IdWt,rdum,cdum)
      call wrcntrl(IUnit,6,'$print',1,IPRNT,rdum,cdum)
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
c----------------------------------------------------------
c Print info
c
      if(ncall.eq.1) call force_infoprt
c----------------------------------------------------------
c
c allocate memory for the total forces on atoms
c
      call getmem(NAt3,lforc3)
c
c put down a memory mark
      call mmark
c     call immark
c
c allocate memory for
c (2) one-electron forces on atoms
c (3) two-electron forces on atoms
c (4) any applied external forces
      call getmem(NAt3,lforc1)
      call getmem(NAt3,lforc2)
      call getmem(NAt3,iafm)
c
c save these addresses as parameters
      call setival('lforc1',lforc1)
      call setival('lforc2',lforc2)
      call setival('lforc3',lforc3)
c
c -- allocate memory for, and read in, density matrices
c ** FOR MP2 Texas matrix system is used **
c   (Note: Memory markers are cleared in <forces2_MP2>)
c
      np4 = 4
      call getival('ncf ',ncf)
      IF(wvfnc(2:3).EQ.'MP') THEN
        call matreset
        call matmark
        call matdef('den0','s',ncf,ncf)
        call matread('den0',np4,'den0_rhf')
        lden = mataddr('den0')
        call setival('ldensi',lden)
      ELSE
        ntri=(ncf*(ncf+1))/2
        call getmem(ntri,lden)
        call setival('ldensi',lden)
        If(.not.rhf) Then
          call getmem(ntri,mden)
          call setival('ldensB',mden)
          call rea(bl(lden),ntri,np4,'dena_uhf')
          call rea(bl(mden),ntri,np4,'denb_uhf')
        Else
          call rea(bl(lden),ntri,np4,'den0_rhf')
        EndIf
      ENDIF
c
c zero-out forces
      call zeroit(bl(lforc1),NAt3)
      call zeroit(bl(lforc2),NAt3)
      call zeroit(bl(lforc3),NAt3)
c
c----------------------------------------------------------
c one-electron part of the gradient
c
      call force1(rhf,bl,bl(ictr),IPRNT)
c----------------------------------------------------------
c
c----------------------------------------------------------
c two-electron part of the gradient
c
      IF(wvfnc(2:3).EQ.'MP') THEN
c -- MP2 forces
        call force2_MP2(natoms,natom,rhf,ictr)
      ELSE
c -- HF/DFT forces
        call force2_SCF(idft,natoms,natom,rhf,ictr)
      ENDIF
c----------------------------------------------------------
c dispersion
c
      If(IDisp.GE.0) Then
        call mmark
        call getmem(3*natom,ixnc)
        call getmem(natom,ian)
        call getnucdat(natom,bl(ixnc),bl(ian))
        call getmem(MaxEL**2,ir0ab)
        call getmem(MaxEL*MaxEL*MaxC*MaxC*3,ic6ab)
        call getmem(natom,icn)
        call getmem(MaxEL,imxc)
        call getmem(MaxEL,itmp)
        write(6,*) ' About to call <gdftd3>  edisp is: ',edisp
        call gdftd3(idft,     IDisp,   iprnt,    natom,    bl(ian),
     $              bl(ixnc), values,  MaxEL,    MaxC,     bl(ir0ab),
     $              bl(ic6ab),bl(icn), bl(imxc), bl(itmp), noabc,
     $              tz,       edisp,   bl(iafm))
        write(6,*) ' Dispersion gradient:'
        call prntmat(natom,3,natom,bl(iafm))
        call retmark
c -- now add the dispersion gradient and write to file
        call RdGRAD(natom,bl(lforc3),'save')      ! recover the gradient
        write(6,*) ' Current gradient:'
        call prntmat(natom,3,natom,bl(lforc3))
        call AddVEC(NAt3,bl(lforc3),bl(iafm),bl(lforc3))
        call WrGRAD(NAtoms,bl(lforc3))
      EndIf
c----------------------------------------------------------
c
      call getival('inuc',inuc)      ! geometry
c
c  Are there external forces? (simulated atomic force microscopy)
      call afm_external(inp,NAtom,bl(inuc),bl(iafm),bl(lforc3))
c
c  Now calculate the sum and torque and print forces
      call force_prn(idft,NAtom,bl)
c----------------------------------------------------------
c
c -- return memory for the total forces
      call retmem(1)
c
c memory status :
c
      call memory_status('end of forces')
c     call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
c
c     write(iout,1100) nreq,nmark,lastadr-ioffset,mxmem,memtot
c1100 format(/' Memory status:'/
c    *' request number=',i4,' memory marks=',i3,
c    *' last used address=',i15,/,
c    *' high water=',i15,' total available memory=',i15)
c----------------------------------------------------------
c restore original values of symmetry pair and atom pointers
c
      call setival('SymFunPr',ifp)
      call setival('SymShPr',ifp1)
c----------------------------------------------------------
      end
c======================================================================
      subroutine force2_SCF(idft,natoms,natom,rhf,ictr)

      use memory

      implicit real*8 (a-h,o-z)
      character jobname*256,datestr*24,wvfnc*20
      Logical rhf
c------------------------------------------------
      common /cpu/ intsize,iacc,icache,memreal
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,
     1              ncs,nsy(4),nsym,nganz(35),lopt(30)
      common /neglect/ eps,eps1,epsr
c------------------------------------------------
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
      common /memmax/ ispblx, maxme1,iforwhat
c------------------------------------------------
c Gradient options :
      common /forcdbl/ thre1,thre2,tchf
      common /forcint/ ncache,nprint,maxprice
c------------------------------------------------
      dimension xintxx(9)
c     common /big/ bl(10000)
c     common /intbl/ifpp,inx(100)
c----------------------------------------------------------------------
      COMMON /DFTCoeff/ aX, XS, XVWN, XVWN5, XB88, XP86L, XP86NL,
     $                  XPW91X, XPW91L, XPW91NL, XLYP, XOptX,
     $                  XPBEX,  XPBEC, EFS, VFS, EFSC, VFSC
c----------------------------------------------------------------------
      call getival('iout',iout)
      call getival('nprint',nprint)
c----------------------------------------------------------------------
c -- check for FTC
      call tstival('ftc0',iyes)
      if(iyes.eq.0) call setival('ftc0',0)
      call getival('ftc0',iftc1)
      if(iftc1.ne.0) write(6,1111)
 1111 Format(' Calculating non-FTC Classical integral derivatives')
c
c -- save some values; not sure why
      reps=eps
      reps1=eps1
      repsr=epsr
      icacher=icache
c----------------------------------------------------------------------
c prepare the integral system for gradient integrals (iforwhat=3) :
c (like twoint95 for ordinary integrals)
c Full list of parameters for twoint95 (13) is reduced here to only
c FOUR . The remaning parameters are taken from commons set up for
c the SCF . Hopefully, it will be always the case for post-scf calc.
c I want to use this interface for all post-scf applications (nmr,grad
c and second derivatives). Four input parameters will be taken from
c appropriate commons set up when a given application is called. Here
c four input parameters taken from FORC option commons :
c
      thresh=thre2
      iforwhat=3                   !  gradient derivative integrals
      if(iftc1.ne.0) iforwhat=13   !  gradient derivative non-ftc integral
c
c value of 13 needed for blocking and (gt.10) for filter in cshneg
c
      call getival('ncs ',ncs)
      call check4lgc(ncs,bl(ictr),lshells,lgcshell,ltype_gc,ldeep_gc)
c
c if l-shells are present then make S,P partitioning
c
      if(lshells.gt.0) then
        call trans_l_2_sp(rhf,bl,bl(ictr),ictr,lshell,'forces')
      endif
c
c if there is a GC basis set but only as deep as 2 and only
c s-functions are GC then we DO NOT profit from GC, especially
c in the gradient and probably hessian. Thus, segment it
c
      if(lgcshell.gt.0.and.ltype_gc.eq.1.and.ldeep_gc.eq.1) then
        call trans_gc_2_sg(rhf,bl,bl(ictr),ictr,nogcsi,'forces')
        write(iout,*)
        write(iout,*)
     * '   General Contracted shells have been segmented for forces'
        write(iout,*)
     * '   because there were only s-doubly gen.contracted orbitals'
        write(iout,*) ' '
      else if(lgcshell.gt.0.and.iftc1.ne.0) then
        call nerror(8,'FORCE module',
     $        'Basis MUST be segmented for FTC forces',0,0)
      endif
c
c ictr is changed on return if l-shells or gc-shells are present
c ..............................................................
c
c ::::::::::::::::::::::::::::::::::::::::::
c -- parallel?
      call getival('nslv',nslv)
      If(nslv.GT.0) call para_JobInit(iforwhat)
c ::::::::::::::::::::::::::::::::::::::::::
c
      If(iftc1.NE.0) call para_FTC_INIT_FORCES(
     $       isharpgrd, isharpness, ncspl,    ncfpl,     iicsplsize,
     $       icssharps, iicsdiff,   ilistsd,  iicspltype, iplbasdat,
     $       rLxmin,    rLxmax,     rLymin,   rLymax,    rLzmin,
     $       rLzmax,    rLx,        rLy,      rLz,       PLDmax,
     $       npwx,      npwy,       npwz,   igridranges,igridranges2,
     $       maxovs,   isharpovs, insharpovs, ii4sharps, icorecut,
     $     iisharpgridrange1, iisharpgridrange2, griddens)
c
      call post_scf(bl,bl(ictr),ncache,nprint,thresh,
     *              iforwhat,0,nblocks,maxbuffer,maxlabels,
     *              scftype,xintxx)
c
      if(scftype.ne.'full-direct')then
        call getival('incore',incore)
        if(incore.gt.0)then
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with core storage'
        else
          write(6,'(3a)')'Integral calculation is ',scftype,
     $              ' with disk storage'
        endif
      endif
c
c The previous call sets up the local integral system
c the following initializes the slaves if we are running
c in parallel. Many input parameters are taken from commons
c or not sent (see comments in ptwoint.f)
c
      call para_post_scf(ncache,iforwhat,nblocks)
c
c *******************************************************
c   contribution to gradient from classical integrals
c *******************************************************
c
c -- reserve memory for labels array
      call getmem(maxlabels,labels)
c
c get addresses for density & two-electron forces on atoms
c (memory allocated in forces.f)
      call getival('ldensi',lden  )
      call getival('lforc2',lforc2)
      If(.not.rhf) call getival('ldensB',mden)
c
c -- calculate two-electron integral derivatives for gradient
c
      call secund(tgrad1)
      call elapsec(egrad1)
c
      ntri = ncf*(ncf+1)/2
      call mmark
c     call immark
c
c -- switch off symmetry if FTC
      If(iftc1.NE.0) call symmoff
      call para_grad(idft,aX,rhf,nblocks,bl,bl(ictr),ntri,thresh,
     *               bl(lden),bl(mden),bl(lforc2),bl(labels))
      If(iftc1.NE.0) call symmon
c
c     call retimark
      call retmark
c
      call secund(tgrad2)
      call elapsec(egrad2)
      total=tgrad2-tgrad1
      elaps=egrad2-egrad1
      write(iout,400) total/60.0d0,elaps/60.0d0
      call f_lush(iout)
  400 format('Master CPU time for 2e part of gradient ='
     *,f8.2,' Elapsed = ',f8.2,' min')
c----------------------------------------------------------------------
c re-scale these parts of fock matrices :
c
      lforc21=lforc2-1
      do 250 ii=1,3*na
      bl(lforc21+ii)=bl(lforc21+ii)*thresh
  250 continue
c------------------------------------------------
c
c -- release memory for labels
      call retmem(1)
c
c -- restore the original SCF values
      eps=reps
      eps1=reps1
      epsr=repsr
      icache=icacher
c
      If(nprint.GT.2) Then
        write(6,*) ' Classical integral contribution to forces is:'
        do i=0,natom-1
        ii = lforc2+3*i
        write(6,*) bl(ii),bl(ii+1),bl(ii+2)
        enddo
      EndIf
c
c -- FTC part done here
c
      If(iftc1.NE.0) call FORCE2_FTC(
     $    natom,      ictr,       griddens,   rLxmin,     rLymin,
     $    rLzmin,     rLx,        rLy,        rLz,        PLDmax,
     $    npwx,       npwy,       npwz,     isharpgrd,    isharpness,
     $    ncspl,   igridranges, igridranges2, ncfpl,      iicsplsize,
     $    icssharps,  iicsdiff,   iicspltype, maxovs,     isharpovs,
     $    insharpovs, ii4sharps,  iplbasdat,  icorecut,
     $  iisharpgridrange1, iisharpgridrange2)
c
c At this point one- and two-electron contributions to the
c gradient are known. Put them together and write them to
c the file <grad>
c
      call force_total(Natom,nprint,bl)
c
c ....................................................................
c transfer back original L-shell or GC-shell basis set info
c
      if(lshells.gt.0) then
        call trans_sp_2_l(bl,bl(ictr),ictr,lshell)
      endif
      if(lgcshell.gt.0.and.ltype_gc.eq.1.and.ldeep_gc.eq.1) then
        call trans_sg_2_gc(bl,bl(ictr),ictr,lgcshell)
      endif
c
c returns the original value of ictr
c ....................................................................
c
c----------------------------------------------------------------------
c   now call the DFT forces
c ..........................................
      IF(idft.GT.0) THEN
c -- DFT part of gradient
        call getmem(0,lastx)
        call retmem(1)
        lcore=igetival('lcore')
        call GRADDFT(lcore-lastx,bl(lastx))
c  DFT forces (if any) calculated & written on <grad>
      ENDIF
c----------------------------------------------------------------------
c return all memory allocated in the 1 and 2-el forces except
c the total force forc3
      call retmark
c     call retimark
c ..........................................
c
c -- reset the values of some FTC flags in order not to screw up
c -- subsequent job steps
      if(iftc1.ne.0) then
        call setival('ftc0',0)       ! switch off FTC
        call setival('nomultipole',1)! switch off multipoles
      endif
      return
      end
c======================================================================
      subroutine force_options(wvfnc,xlow,IdWt)
      implicit real*8 (a-h,o-z)
C
C  read in force options
C
C  ARGUMENTS
C
C  wvfnc   -  wavefunction type
C  xlow    -  lowest eigenvalue of overlap matrix
C             (used to set 2-el integral threshold)
C  IdWt    -  DFT weight derivative flag
C
      character wvfnc*20
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c---------------------------------------------------------------------
c options for integral calculations (from readin):
      common /intgop/ ncachx,maxpricx,iiii(2)
c---------------------------------------------------------------------
      common /forcdbl/ thre1,thre2,tchf
      common /forcint/ ncache,nprint,maxprice
      common /intlim/ limxmem,limblks,limpair
c---------------------------------------------------------------------
      parameter (nopt=6)
      character*4 word(nopt)
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      character*256 chopv(nopt)
      character*20 cdum
c---------------------------------------------------------------------
c 0=logical option, 1=a single integer, 2=two integers, 3=three integers
c 11=a single real, 12=two reals, 13=three reals, 21=a character option
c
      data ioptyp/1,     11,    11,     1,     3,     0/
      data word /'prin','thr1','thr2','ncac','limi','wght'/
c---------------------------------------------------------------------
c maxprice (ipay) CAN NOT be changed :
c
      maxprice=maxpricx
c
c---------------------------------------------------------------------
c thre1=1.0d-11   ! one-el. derivative integral threshold
c thre2=1.0d-10   ! two-el. derivative integral threshold
c thre3=1.0d-9    ! MP2 two-el integral threshold
c---------------------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
c
      write(iout,*)
     *  '                The Analytical Forces Module '
      write(iout,*)' '
      if(wvfnc.eq.'RHF ') then
         write(iout,*)
     *  '                            RHF              '
      endif
      if(wvfnc.eq.'UHF ') then
         write(iout,*)
     *  '                            UHF              '
      endif
      if(wvfnc.eq.'RDFT') then
         write(iout,*)
     *  '                          RHF/DFT            '
      endif
      if(wvfnc.eq.'UDFT') then
         write(iout,*)
     *  '                          UHF/DFT            '
      endif
      if(wvfnc.eq.'RMP2') then
         write(iout,*)
     *  '                          RHF/MP2            '
      endif
      if(wvfnc.eq.'UMP2') then
         write(iout,*)
     *  '                          UHF/MP2            '
      endif
      if(wvfnc.eq.'RMP2-dual') then
         write(iout,*)
     *  '                        RHF/Dual-MP2         '
      endif
      if(wvfnc.eq.'UMP2-dual') then
         write(iout,*)
     *  '                        UHF/Dual-MP2         '
      endif
      write(iout,*)' '
c
      call f_lush(iout)
c---------------------------------------------------------------------
c default values :
c
      nprint=1
      thre1=1.0d-11
      thre2=1.0d-10
      thre2=MIN(thre2,xlow**2)
      If(thre2.LT.1.0d-11) then
        write(iout,2000) xlow
 2000 Format('**WARNING** Smallest eigenvalue of overlap is ',d12.4,/,
     $       '  Final gradient may be numerically inaccurate')
        thre2 = 1.0d-11
      EndIf
      ncache=ncachx
      IdWt = 0           ! no weight derivatives in DFT gradient
c
      limxmem=800000
      limblks=0
      limpair=100
c
      call readopt(inp,nopt,word,ioptyp,iopv,ropv,chopv,ifound)
c
      if(ifound(1).gt.0) nprint=iopv(1,1)
      if(ifound(2).gt.0) thre1=10.d0**(-ropv(1,2))
      if(ifound(3).gt.0) thre2=10.d0**(-ropv(1,3))
      if(ifound(4).gt.0) ncache=iopv(1,4)
      if(ifound(5).gt.0) then
        limxmem=iopv(1,5)
        limblks=iopv(2,5)
        limpair=iopv(3,5)
      endif
      if(ifound(6).gt.0) IdWt=1
c---------------------------------------------------------------------
      call setival('nprint',nprint)
      call setrval('thr1',thre1)
      call setrval('thr2',thre2)
      call setrval('thr3',thre2)  ! not used in grad., needed in post_scf
c---------------------------------------------------------------------
c print options :
c
      ii = 0
      do iop=1,nopt
      ii = ii + ifound(iop)
      enddo
c
      IF(ii.gt.0) THEN
        write(iout,200)
        do 100 iop=1,nopt
         if(ifound(iop).gt.0) then
           if(iop.eq.1) write(iout,210) word(iop),iopv(1,iop)
           if(iop.eq.2) write(iout,220) word(iop),thre1
           if(iop.eq.3) write(iout,220) word(iop),thre2
           if(iop.eq.4) write(iout,220) word(iop),iopv(1,iop)
           if(iop.eq.5) write(iout,230) word(iop),iopv(1,iop),
     *                                           iopv(2,iop),
     *                                           iopv(3,iop)
           if(iop.eq.6) write(iout,240)
         endif
  100   continue
        write(iout,200)
      ENDIF
c---------------------------------------------------------------------
      lopt(1)=nprint
c---------------------------------------------------------------------
  200 format(/58('-'))
  210 format ('Force option = ',a4,'  is on ; value = ',i10  )
  220 format ('Force option = ',a4,'  is on ; value = ',e12.5)
  230 format ('Force option = ',a4,'  is on ; value = ',3(i10,2x))
  240 format (' Weight Derivatives will be used in DFT gradient')
c---------------------------------------------------------------------
      end
c======================================================================
      subroutine force_infoprt
      implicit real*8 (a-h,o-z)
      common /fieldx/ xfield,xfgrad,elfiel(9)
c----------------------------------------------------------------------
      call getival('iout',iout)
      call getival('icond',icond)
c----------------------------------------------------------------------
c
c electric field and electric field gradient :
c
      ifield=xfield
      ifgrad=xfgrad
c
      if( ifield .eq.1) write (iout,2820) elfiel(1),elfiel(2),elfiel(3)
      if( ifgrad .eq.1) write (iout,2821) elfiel(4),elfiel(5),elfiel(6),
     *                                    elfiel(7),elfiel(8),elfiel(9)
c------------------------------------------------
  223 format(58(1H*))
  280 format
     *('              CALCULATIONS OF THE GRADIENT                ')
 2820 format
     *('               external electric field is                 '/
     *'            fx=',f6.4,2x,' fy=',f6.4,2x,' fz=',f6.4,
     *'             '/)
 2821 format
     *('          external electric field gradient is             '/
     *'            xx=',f6.4,2x,' yx=',f6.4,2x,' yy=',f6.4,
     *'             '/
     *'            zx=',f6.4,2x,' zy=',f6.4,2x,' zz=',f6.4,
     *'             '/)
c
      end
c======================================================================
      subroutine trans_l_2_sp(rhf,bl,inx,ictr_sp,lshell,called4)
c
c     transform basis set containing L-shell into S,P ones
c     and the density matrix from L- to SP-shell basis form
c
      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      character*(*) called4
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c     common /intbl/ifpp,inxx(100)
      dimension bl(*)
      dimension inx(12,*)
c-----------------------------------------------------------
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ictr',ictr)
c
c-----------------------------------------------------------
c
c first get ncs_sp & nsh_sp  from  ncs_l & nsh_l
c
      lshell=0
      ics_sp=0
      ish_sp=0
      do ics=1,ncs
         ics_sp=ics_sp+1
         lcontr=inx(5,ics)-inx(1,ics)
         ish_sp=ish_sp+lcontr
         itype=inx(12,ics)
         if(itype.eq.3) then
            lshell=lshell+1
            ics_sp=ics_sp+1
            ish_sp=ish_sp+lcontr
         endif
      enddo
c
      if(lshell.eq.0) return
c
      ncs_sp=ics_sp
      nsh_sp=ish_sp
c-----------------------------------------------------------
c allocate memory for sp-basis set :
c
      call getmem(nsh_sp*13,ibas_sp)
      call getint(ncs_sp*12,ictr_sp)
c
c these two allocations are released in trans_sp_2_L
c-----------------------------------------------------------
c zero out the new basis set data (to avoid weird behavior
c with some compilers)      ! MM  Sep. 2006
c
      call zeroit(bl(ibas_sp),nsh_sp*13)
      call izeroit(bl(ictr_sp),ncs_sp*12)
c
c make segmented basis set out of l-shell type one
c
      call make_sp_bas(ncs,nsh,inx,bl(ibas),
     *                 ncs_sp,nsh_sp,bl(ictr_sp),bl(ibas_sp))
c-----------------------------------------------------------
c reorder sp-basis set
c
      call getint(ncf,ireor)
      call getint(ncs_sp,iswap)
      call getint(12*ncs_sp,ncso)
      call getmem(3*na,ix)
      call getmem(na,ian)
      call getnucdat(na,bl(ix),bl(ian))
      call reorder(ncs_sp,bl(ictr_sp),bl(ireor),
     *             bl(iswap),bl(ncso),bl(ian))
      call retmem(4) ! ian, ix, ncso, iswap
c-----------------------------------------------------------
c get the L-basis set info
c
      call getival('ncs ',ncs_lsh)
      call getival('ibas',ibas_lsh)
      call getival('nsh',nsh_lsh)
      call getival('ictr',ictr_lsh)
c
c save the L-basis set info
c
      call setival('ncs_lsh ',ncs_lsh)
      call setival('ibas_lsh',ibas_lsh)
      call setival('nsh_lsh',nsh_lsh)
      call setival('ictr_lsh',ictr_lsh)
c
c replace L-shell by SP-shell basis set info :
c
      call setival('ncs ',ncs_sp)
      call setival('ibas',ibas_sp)
      call setival('nsh',nsh_sp)
      call setival('ictr',ictr_sp)
c
c make the same replacement in the common /ganz/
c
      ibas=ibas_sp
      ncs =ncs_sp
      nsh =nsh_sp
c
c-----------------------------------------------------------
      if(called4.eq.'basiss') then
c        call retint(1)     ! ireor
         RETURN
      endif
c-----------------------------------------------------------
c-----------------------------------------------------------
c setup symmetry info
c
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
      call setival('sfpr_sp',ifp)
      call setival('sspr_sp',ifp1)
      call symfunc(.false.)
c-----------------------------------------------------------
c transfer the original "L-shell" density into sp-basis set and
c put it in the old (l-shell) location (orig.dens. destroyed)
c
      ntri=ncf*(ncf+1)/2
c-----------------------------------------------------------
c called for NMR Chemical Shift 
c
      if(called4.eq.'nmrshift') then
         call getival('lden',lden)  ! from chshift
         call getival('lvec',lvec)  ! from chshift
c
         call getmem(ntri, lden_sp)
         call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
         call retmem(1)
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         call retmem(1)
         return
      endif
c
c-----------------------------------------------------------
c...........................................................
CSS   the following is added for mp2_gradients
c
      if(called4.eq.'mp2grad') then
         lden=mataddr('den0')
         call getmem(ntri, lden_sp)
         call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
         call retmem(1)
         lfoc=mataddr('fock')
         call getmem(ntri, lfoc_sp)
         call make_sp_den(bl(ireor),ncf,bl(lfoc),bl(lfoc_sp))
         call retmem(1)
         lovl=mataddr('ovla')
         call getmem(ntri, lovl_sp)
         call make_sp_den(bl(ireor),ncf,bl(lovl),bl(lovl_sp))
         call retmem(1)
         lvec=mataddr('cano')
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         call retmem(1)
         return
      endif  !  if(called4.eq.'mp2grad') then
c
c-----------------------------------------------------------
c...........................................................
c     called4.eq.'forces' and   called4.eq.'hessian'
c
c density :
         call getival('ldensi',lden)
         call getmem(ntri, lden_sp)
         call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
         If(.not.rhf) Then
            call getival('ldensB',mden)
            call make_sp_den(bl(ireor),ncf,bl(mden),bl(lden_sp))
         EndIf
         call retmem(1)     ! lden_sp
c
c
      if(called4.eq.'hessian') then
c
c overlap ( needed only for level shifted CPHF):
         call tstival('lsmat',is0)
         if(is0.eq.1) then
            call getival('lsmat',ls0)
            call getmem(ntri, ls0_sp)
            call make_sp_den(bl(ireor),ncf,bl(ls0),bl(ls0_sp))
            call retmem(1)     ! ls0_sp
         endif
c
         call getival('lfock0',lfoc)
         call getmem(ntri, lfoc_sp)
         call make_sp_den(bl(ireor),ncf,bl(lfoc),bl(lfoc_sp))
         if(.not.rhf)then
           call getival('lfock0B',mfoc)
           call make_sp_den(bl(ireor),ncf,bl(mfoc),bl(lfoc_sp))
         endif
         call retmem(1)     ! lfoc_sp
c
c        weighted density :
         call getival('ldewsi',ldew)
         call getmem(ntri, ldew_sp)
         call make_sp_den(bl(ireor),ncf,bl(ldew),bl(ldew_sp))
         call retmem(1)     ! ldew_sp
c
         call getival('lvec',lvec)
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         if(.not.rhf) then
           call getival('lvecB',mvec)
           call make_sp_eig(bl(ireor),ncf,bl(mvec),bl(lvec_sp))
         endif
         call retmem(1)
c -- for VCD rearrange perturbed magnetic wavefunction on disk
         call tstival('vcd',ivcd)
c        If(ivcd.NE.0) Then
c           call getival('nocc',nocc)
c           call getmem(ncf*nocc*3,lpve)
c           call getmem(ncf*nocc,lpve_sp)
c           call rea(bl(lpve),3*ncf*nocc,2,'c1_magn')
c           call make_sp_occ(bl(ireor),nocc,ncf,bl(lpve),bl(lpve_sp))
c           call wri(bl(lpve),ncf*nocc*3,2,1,'c1_magn')
c           call retmem(2)
c        EndIf
      endif
c-----------------------------------------------------------
c
      end
c======================================================================
      subroutine make_sp_occ(ireor,nocc,ncf,vec, vec_sp)
c reorders three perturbed orb. coeff. matrices
      implicit real*8 (a-h,o-z)
      dimension vec(ncf,nocc,3)
      dimension vec_sp(ncf,nocc)
      dimension ireor(ncf)
c
      do m=1,3
      do jor=1,nocc
         do j=1,ncf
            jnew=ireor(j)
            vec_sp(jnew,jor)=vec(j,jor,m)
         enddo
      enddo
c
      call dcopy(ncf*nocc,vec_sp(1,1),1,vec(1,1,m),1)
      enddo
c
      end
c======================================================================
      subroutine make_sp_bas(ncs,nsh,inx,datbas,
     *                       ncs_sp,nsh_sp,inx_sp,datbas_sp)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*), inx_sp(12,*)
      dimension datbas(13,*),datbas_sp(13,*)
      data zero/0.0d0/
c
      ips_pp=0   ! counter of p's from l-shells
      ics_sp=0
      ips_sp=0
      icf_sp=0
      lshell=0
      do ics=1,ncs
         ics_sp=ics_sp+1
         itype=inx(12,ics)
         ish_b=inx(1,ics)+1
         ish_e=inx(5,ics)
         nprim=inx(5,ics)-inx(1,ics)
         if(itype.eq.3) then
            lshell=lshell+1
c contracted functions from (11,*)+1 to (10,*)
            inx_sp(11,ics_sp)=inx(11,ics)     ! s-part
            inx_sp(10,ics_sp)=inx(11,ics)+1   ! s-part
            inx_sp(11,ics_sp+1)=inx(11,ics)+1 ! p-part
            inx_sp(10,ics_sp+1)=inx(10,ics)   ! p-part
c primitive shells  from (1,*)+1 to (5,*)
c leave s-primitives where l's were and put p-primitives
c at the end starting from nsh:
            inx_sp(1,ics_sp)=inx(1,ics)                 ! s-part
            inx_sp(5,ics_sp)=inx(5,ics)                 ! s-part
            inx_sp(1,ics_sp+1)=nsh + ips_pp             ! p-part
            inx_sp(5,ics_sp+1)=inx_sp(1,ics_sp+1)+nprim ! p-part
c
c make datbas:
            ii1=0
            do isl=ish_b,ish_e
               ii1=ii1+1
               iss=inx_sp(1,ics_sp  ) + ii1
               datbas_sp(1,iss)=datbas(1,isl)   ! exponent
               datbas_sp(2,iss)=datbas(2,isl)   ! s-coeff
c
               isp=inx_sp(1,ics_sp+1) + ii1
               datbas_sp(1,isp)=datbas(1,isl)   ! exponent
               datbas_sp(2,isp)=datbas(3,isl)   ! p-coeff
                  do izero=3,10
                     datbas_sp(izero,iss)=zero
                     datbas_sp(izero,isp)=zero
                  enddo
c s-part
               datbas_sp(11,iss)=datbas(11,isl)
               datbas_sp(12,iss)=datbas(12,isl)
               datbas_sp(13,iss)=datbas(13,isl)
c p-part
               datbas_sp(11,isp)=datbas(11,isl)
               datbas_sp(12,isp)=datbas(12,isl)
               datbas_sp(13,isp)=datbas(13,isl)
            enddo
c
c type of shells
            inx_sp(12,ics_sp)=1
            inx_sp(12,ics_sp+1)=2
c length of shells
            inx_sp(3,ics_sp)=1
            inx_sp(3,ics_sp+1)=3
c general contr.
            inx_sp(4,ics_sp)=0
            inx_sp(4,ics_sp+1)=0
c center of shells
            inx_sp(2,ics_sp)=inx(2,ics)
            inx_sp(2,ics_sp+1)=inx(2,ics)
c
            ics_sp=ics_sp+1
            ips_sp=ips_sp+nprim
            ips_pp=ips_pp+nprim
         else
c contracted functions from (11,*)+1 to (10,*)
            inx_sp(11,ics_sp)=inx(11,ics)
            inx_sp(10,ics_sp)=inx(10,ics)
c primitive shells  from (1,*)+1 to (5,*)
            inx_sp(1,ics_sp)=inx(1,ics)
            inx_sp(5,ics_sp)=inx(5,ics)
c
c make datbas:
            ii1=0
            do isl=ish_b,ish_e
               ii1=ii1+1
               iss=inx_sp(1,ics_sp  ) + ii1
               datbas_sp(1,iss)=datbas(1,isl)   ! exponent
               datbas_sp(2,iss)=datbas(2,isl)   ! coeff
c
                  do izero=3,10
ccc                  datbas_sp(izero,iss)=zero
                     datbas_sp(izero,iss)=datbas(izero,isl)
                  enddo
               datbas_sp(11,iss)=datbas(11,isl)
               datbas_sp(12,iss)=datbas(12,isl)
               datbas_sp(13,iss)=datbas(13,isl)
            enddo
c
c type of shells
            inx_sp(12,ics_sp)=inx(12,ics)
c length of shells
            inx_sp(3,ics_sp)=inx(3,ics)
c general contr.
            inx_sp(4,ics_sp)=inx(4,ics)
c center of shells
            inx_sp(2,ics_sp)=inx(2,ics)
         endif
         ips_sp=ips_sp+nprim
      enddo
c----------------------------------------------------------------
c
      end
c======================================================================
      subroutine make_sp_den(ireor,ncf,den_lsh,den_sp)
      implicit real*8 (a-h,o-z)
      dimension den_lsh(*)
      dimension den_sp(*)
      dimension ireor(ncf)
c
      ijcf=0
      do icf=1,ncf
         inew=ireor(icf)
         iine=inew*(inew-1)/2
         do jcf=1,icf
            jnew=ireor(jcf)
            ijcf=ijcf+1
            ijne=iine+jnew
            if(jnew.gt.inew) ijne=jnew*(jnew-1)/2+inew
            den_sp(ijne)=den_lsh(ijcf)
         enddo
      enddo
c
      do ij=1,ncf*(ncf+1)/2
         den_lsh(ij)=den_sp(ij)
      enddo
c
      end
c======================================================================
      subroutine force_total(NAtom,iprnt,bl)
      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension bl(*)
c
      call getival('lforc1',lforc1)
      call getival('lforc2',lforc2)
      call getival('lforc3',lforc3)
c
      call grad_tot(NAtom,ncs,iprnt,bl(inuc),bl(lforc1),
     $              bl(lforc2),bl(lforc3))
c
      end
c======================================================================
      subroutine grad_tot(natoms, ncs,    iprnt,  datnuc,  forc1,
     $                    forc2,  forc3)

      use memory

      implicit real*8 (a-h,o-z)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /big/bl(10000)
      dimension datnuc(5,*)
      dimension forc1(3,natoms),forc2(3,natoms),forc3(3,natoms)
      parameter (zero=0.0d0,thrsh=1.0d-6)
c---------------------------------------------------------------------
c INPUT :
c
c natoms - number of atoms
c ncs    - number of shells
c forc1  - one-electron part of the forces
c forc2  - two-electron part of the forces
c forc3  - nuclear repulsion forces
c
c OUTPUT:
c
c forc3 - total forces on atoms
c
c---------------------------------------------------------------------
c
      if(iprnt.gt.1) write(iout,152)
  152 format(
     *'Atom  Name       force-x    force-y    force-z     Coordinates'/)
c
      do iat=1,natoms
c
         fx1=forc1(1,iat)
         fy1=forc1(2,iat)
         fz1=forc1(3,iat)
c
         fx2=forc2(1,iat)
         fy2=forc2(2,iat)
         fz2=forc2(3,iat)
c
         fx3=forc3(1,iat)
         fy3=forc3(2,iat)
         fz3=forc3(3,iat)
c
         forc3(1,iat)=fx1+fx2+fx3
         forc3(2,iat)=fy1+fy2+fy3
         forc3(3,iat)=fz1+fz2+fz3
c
        if(iprnt.gt.1) then
          xat=datnuc(2,iat)
          yat=datnuc(3,iat)
          zat=datnuc(4,iat)
          write(iout,210) iat,fx1,fy1,fz1,xat,yat,zat,
     *                        fx2,fy2,fz2,
     *                        fx3,fy3,fz3,
     *                forc3(1,iat),forc3(2,iat),forc3(3,iat)
        endif
c
      end do
c
c -- forces summed up from 3 contributions (overlap, 1-el, 2-el)
c -- write them to the <control> file
      call wrgrad(natoms,forc3)
c
      return
c
  210 format(/i3,'  (1-el.)  ',3(f10.5,1x),1x,3(f9.5,1x)/
     *        3x,'  (2-el.)  ',3(f10.5,1x)/
     *        3x,'  (nucl.)  ',3(f10.5,1x)/
     *        3x,' -----------------------------------------',/
     *        3x,'  sum      ',3(f10.5,1x))
c
      end
c======================================================================
      subroutine force_prn(idft,NAtoms,bl)
c  This is the main routine for torque
      implicit real*8 (a-h,o-z)
      dimension bl(*)
c
      call getival('lforc3',lforc3)
      call getival('inuc',inuc)
c
c   Read back the gradients from the <grad> file
      call rdgrad(NAtoms,bl(lforc3),'save')
c
      call torque(NAtoms,idft,bl(inuc),bl(lforc3))
c
      end
c======================================================================
      subroutine torque(natoms,idft,datnuc,forc3)

      use memory

      implicit real*8 (a-h,o-z)
c
c  This routine calculates the sum and torque and the max.
c  Cartesian component of the forces & prints them out
c
      common /forcint/ ncache,nprint,maxprice
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
c     common /big/bl(10000)
      dimension datnuc(5,*)
      dimension forc3(3,natoms)
      parameter (zero=0.0d0,thrsh=1.0d-6)
c---------------------------------------------------------------------
c INPUT :
c
c natoms - number of atoms
c idft   - dft flag  0 - no dft
c                   >0 - dft contribution to forces needs to be computed
c forc3  -  space for total forces (will be read in from file)
c
c OUTPUT:
c
c forc3 - total forces on atoms
c
c---------------------------------------------------------------------
       if(nprint.ge.1) write(iout,151)
  151 format(
     *'Atom  Name       force-x    force-y    force-z'/)
c
      sumx=zero
      sumy=zero
      sumz=zero
      torqx=zero
      torqy=zero
      torqz=zero
      do iat=1,natoms
         xat=datnuc(2,iat)
         yat=datnuc(3,iat)
         zat=datnuc(4,iat)
         sumx=sumx+forc3(1,iat)
         sumy=sumy+forc3(2,iat)
         sumz=sumz+forc3(3,iat)
c
c torque= RxF (cross product of the coordinate and the the force
         torqx=torqx+forc3(3,iat)*yat-forc3(2,iat)*zat
         torqy=torqy+forc3(1,iat)*zat-forc3(3,iat)*xat
         torqz=torqz+forc3(2,iat)*xat-forc3(1,iat)*yat
c
         if(nprint.ge.1) then
            write(iout,200) iat,datnuc(5,iat),
     *             forc3(1,iat),forc3(2,iat),forc3(3,iat)
         endif
c
      enddo
c
  200 format(i3,5x,a3,3x,3(f10.7,1x))
c
      formax=abs(max(sumx,sumy,sumz,torqx,torqy,torqz))
      formin=abs(min(sumx,sumy,sumz,torqx,torqy,torqz))
      forabs=max(formin,formax)
c
      if((forabs.gt.thrsh.and.(idft.eq.0.or.nprint.gt.1)).or.
     1   forabs.gt.thrsh*50.0d0) then
         write(iout ,220) sumx,sumy,sumz,torqx,torqy,torqz
  220    format('Sum or total torque of the forces does not vanish',/,
     *   'Sum=',3f10.7,' Torque=',3f10.7)
      endif
       call absmax(3*natoms,forc3,ifmax,fmax)
ckw   ifmax=(ifmax-2)/3+1
ckw   write(iout,230) fmax,ifmax
      iatom=(ifmax-2)/3+1
      iatom_x=(iatom-1)*3+ 1
      iatom_y=iatom_x+1
      iatom_z=iatom_y+1
      if(ifmax.eq.iatom_x.or.ifmax.eq.iatom_y.or.ifmax.eq.iatom_z) then
         iatomx=iatom
      else
         iatomx=iatom+1
      endif
      write(iout,230) fmax,iatomx
  230 format(/,' Maximum Cartesian force component=',f12.7,
     1  ' Eh/a0 on atom ',i4)
c
      end
c======================================================================
      subroutine trans_sp_2_l(bl,inx,ictr,lshell)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension bl(*)
      dimension inx(12,*)
c-----------------------------------------------------------
      if(lshell.eq.0) return
c-----------------------------------------------------------
c get back the L-basis set info seved in trans_l_2_sp
c
      call getival('ncs_lsh ',ncs)
      call getival('ibas_lsh',ibas)
      call getival('nsh_lsh',nsh)
      call getival('ictr_lsh',ictr)
c
c replace SP-shell by L-shell basis set info :
c
      call setival('ncs ',ncs)
      call setival('ibas',ibas)
      call setival('nsh',nsh)
      call setival('ictr',ictr)
c
c-----------------------------------------------------------
c get back symmetry info
c
      call getival('sfpr_sp',ifp)
      call getival('sspr_sp',ifp1)
      call setival('SymFunPr',ifp)
      call setival('SymShPr',ifp1)
cc      call symfunc(.false.)
c-----------------------------------------------------------
c release memory alloctated in trans_l_2_sp
c
c     call retint(1)
      call retmem(1)
c-----------------------------------------------------------
c
      end
c======================================================================
      subroutine make_sp_eig(ireor,ncf,vec, vec_sp)
      implicit real*8 (a-h,o-z)
      dimension vec(*)            ! inp: lsh  / out: sp shell
      dimension vec_sp(*)
      dimension ireor(ncf)
c
      do jor=1,ncf
         je=jor*ncf
         ja=je-ncf
         do j=1,ncf
            jnew=ireor(j)
            vec_sp(ja+jnew)=vec(ja+j)
         enddo
      enddo
c
      call dcopy(ncf*ncf,vec_sp(1),1,vec(1),1)
c
      end
c======================================================================
      subroutine trans_gc_2_sg(rhf,bl,inx,ictr_sg,nogcs,called4)
c
c  transforms basis set containing Gen.Contr.basis set into segmented one
c
      use memory, block => bl

      implicit real*8 (a-h,o-z)
      Logical rhf
      character*(*) called4
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c     common /intbl/ifpp,inxx(100)
      dimension bl(*)
      dimension inx(12,*)
c-----------------------------------------------------------
clublin
c     call getival('lcore',lcore)
c     call getival('nsym',nsym)
clublin
c
      call getival('na  ',na)
      call getival('ncf ',ncf)
      call getival('ncs ',ncs)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('nbf',nbf)
      call getival('nsh',nsh)
      call getival('ictr',ictr)
c
c-----------------------------------------------------------
c first get ncs_sg & nsh_sg  from  ncs_gc & nsh_gc
c
      igcshell=0
      ics_sg=0
      ish_sg=0
      do ics=1,ncs
         lcontr=inx(5,ics)-inx(1,ics)
         lgenc1=inx(4,ics)+1                 ! 1 for segmented
         ics_sg=ics_sg+lgenc1                ! segmented contracted shells
         ish_sg=ish_sg+lgenc1*lcontr         ! segmented primitive shells
         if(lgenc1.gt.1) igcshell=igcshell+1 ! general contracted shells
      enddo
c
      nogcs=igcshell
      if(igcshell.eq.0) return
c
      ncs_sg=ics_sg
      nsh_sg=ish_sg
c-----------------------------------------------------------
c     write(6,*)' no of gc shells=', nogcs
c     write(6,*)'0 ncs,nsh=',ncs,nsh,' ncs_sg,nsh_sg=',ncs_sg,nsh_sg
c-----------------------------------------------------------
c allocate memory for segmented sg-basis set :
c
      call getmem(nsh_sg*13,ibas_sg)
      call getint(ncs_sg*12,ictr_sg)
c
c these two allocations are released in trans_sg_2_gc
c-----------------------------------------------------------
      call getint(ncf,ireor)
c-----------------------------------------------------------
      call getint(ncf,ireor1)
      call getint(ncf,ireor2)
c-----------------------------------------------------------
c zero out the new basis data, in order to avoid some weird
c behavior with some compilers
c
      call zeroit(bl(ibas_sg),nsh_sg*13)
      call izeroit(bl(ictr_sg),ncs_sg*12)
c
c make segmented basis set out of gc-shell type one
c
      call make_sg_bas(ncs,nsh,inx,bl(ibas),bl(ireor1),
     *                 ncs_sg,nsh_sg,bl(ictr_sg),bl(ibas_sg))
c-----------------------------------------------------------
c reorder sg-basis set
c
      call getint(ncs_sg,iswap)
      call getint(12*ncs_sg,ncso)
      call getmem(3*na,ix)
      call getmem(na,ian)
      call getnucdat(na,bl(ix),bl(ian))
      call reorder(ncs_sg,bl(ictr_sg),bl(ireor2),
     *             bl(iswap),bl(ncso),bl(ian))
      call retmem(4)   ! ix,ian, iswap,ncso
c-----------------------------------------------------------
      call mapf_gc_sg_ord(ncf,bl(ireor1),bl(ireor2),bl(ireor) )
c
c     inp: ireor1 - from GC to NOT ordered SG
c     inp: ireor2 - from NOT ordered SG to ordered SG
c     out: ireor  - from GC to ordered SG
c
c     call retint(2)   ! ireor2,ireor1
c-----------------------------------------------------------
c get the GC-basis set info
c
      call getival('ncs ',ncs_gcsh)
      call getival('ibas',ibas_gcsh)
      call getival('nsh',nsh_gcsh)
      call getival('ictr',ictr_gcsh)
c
c save the GC-basis set info
c
      call setival('ncs_gcsh ',ncs_gcsh)
      call setival('ibas_gcsh',ibas_gcsh)
      call setival('nsh_gcsh',nsh_gcsh)
      call setival('ictr_gcsh',ictr_gcsh)
c
c replace GC-shell by SG-shell basis set info :
c
      call setival('ncs ',ncs_sg)
      call setival('ibas',ibas_sg)
      call setival('nsh',nsh_sg)
      call setival('ictr',ictr_sg)
c
c make the same replacement in the common /ganz/
c
      ibas=ibas_sg
      ncs =ncs_sg
      nsh =nsh_sg
c
c-----------------------------------------------------------
c setup symmetry info
c
      call symfunc(.false.)
c-----------------------------------------------------------
      if(called4.eq.'basiss') then
c        call retint(1)     ! ireor
         RETURN
      endif
c-----------------------------------------------------------
c transfer the original "GC-shell" density into SG-basis set and
c put it in the old (gc-shell) location (orig.dens. destroyed)
c
      ntri=ncf*(ncf+1)/2
c-----------------------------------------------------------
c called for NMR Chemical Shift 
c
      if(called4.eq.'nmrshift') then
         call getival('lden',lden)  ! from chshift
         call getival('lvec',lvec)  ! from chshift
c
         call getmem(ntri, lden_sp)
         call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
         call retmem(1)
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         call retmem(1)
         return
      endif
c
c-----------------------------------------------------------
c...............................................................
CSS   the following is added for mp2_gradients
      if(called4.eq.'mp2grad') then
         lden=mataddr('den0')
         call getmem(ntri, lden_sp)
         call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
         call retmem(1)
         lfoc=mataddr('fock')
         call getmem(ntri, lfoc_sp)
         call make_sp_den(bl(ireor),ncf,bl(lfoc),bl(lfoc_sp))
         call retmem(1)
         lovl=mataddr('ovla')
         call getmem(ntri, lovl_sp)
         call make_sp_den(bl(ireor),ncf,bl(lovl),bl(lovl_sp))
         call retmem(1)
         lvec=mataddr('cano')
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         call retmem(1)
c        weighted density :
         call getival('ldewsi',ldew)
         call getmem(ntri, ldew_sp)
         call make_sp_den(bl(ireor),ncf,bl(ldew),bl(ldew_sp))
         call retmem(1)     ! ldew_sp
c        call retint(1)     ! ireor
         return
      endif  !  if(called4.eq.'mp2grad') then
c
c     called4.eq.'forces' and   called4.eq.'hessian'
c
      call getival('ldensi',lden)
      call getmem(ntri, lden_sp)
      call make_sp_den(bl(ireor),ncf,bl(lden),bl(lden_sp))
      If(.not.rhf) Then
         call getival('ldensB',mden)
         call make_sp_den(bl(ireor),ncf,bl(mden),bl(lden_sp))
      EndIf
      call retmem(1)     ! lden_sp
c
      if(called4.eq.'hessian') then
         call getival('lfock0',lfoc)
         call getmem(ntri, lfoc_sp)
         call make_sp_den(bl(ireor),ncf,bl(lfoc),bl(lfoc_sp))
         if(.not.rhf)then
           call getival('lfock0B',mfoc)
           call make_sp_den(bl(ireor),ncf,bl(mfoc),bl(lfoc_sp))
         endif
         call retmem(1)     ! lfoc_sp
c
c        weighted density :
         call getival('ldewsi',ldew)
         call getmem(ntri, ldew_sp)
         call make_sp_den(bl(ireor),ncf,bl(ldew),bl(ldew_sp))
         call retmem(1)     ! ldew_sp
c
         call getival('lvec',lvec)
         call getmem(ncf*ncf,lvec_sp)
         call make_sp_eig(bl(ireor),ncf,bl(lvec),bl(lvec_sp))
         if(.not.rhf)then
           call getival('lvecB',mvec)
           call make_sp_eig(bl(ireor),ncf,bl(mvec),bl(lvec_sp))
         endif
         call retmem(1)
c -- for VCD rearrange perturbed magnetic wavefunction on disk
c        call tstival('vcd',ivcd)
c        If(ivcd.NE.0) Then
c           call getival('nocc',nocc)
c           call getmem(ncf*nocc*3,lpve)
c           call getmem(ncf*nocc,lpve_sp)
c           call rea(bl(lpve),3*ncf*nocc,2,'c1_magn')
c           call make_sp_occ(bl(ireor),nocc,ncf,bl(lpve),bl(lpve_sp))
c           call wri(bl(lpve),ncf*nocc*3,2,1,'c1_magn')
c           call retmem(2)
c        EndIf
      endif
c-----------------------------------------------------------
c release memory
c
c     call retint(1)     ! ireor
c-----------------------------------------------------------
c
      end
c======================================================================
      subroutine make_sg_bas(ncs,nsh,inx,datbas,ireor,
     *                       ncs_sg,nsh_sg,inx_sg,datbas_sg)
      implicit real*8 (a-h,o-z)
      dimension inx(12,*), inx_sg(12,*),ireor(*)
      dimension datbas(13,*),datbas_sg(13,*)
      data zero/0.0d0/
c
      ics_sg=0                     ! seg. contr. shells
      ips_sg=0                     ! seg. primi. shells
      icf_sg=0                     !
      do ics=1,ncs
         isize =inx(3,ics)
         igenc =inx(4,ics)
         igenc1=inx(4,ics)+1
         ish_b=inx(1,ics)+1
         ish_e=inx(5,ics)
         lcontr=inx(5,ics)-inx(1,ics)
         do igc=0,0
            ics_sg=ics_sg+1
c contracted functions from (11,*)+1 to (10,*)
            inx_sg(11,ics_sg)=icf_sg
            inx_sg(10,ics_sg)=inx_sg(11,ics_sg)+isize
c
            incr_gc=inx(11,ics)+igc*isize
            do icf=1,isize
               icf_gc=icf+incr_gc
               ireor(icf_gc)=icf_sg+icf
            enddo
c
c primitive shells  from (1,*)+1 to (5,*)
            inx_sg(1,ics_sg)=ips_sg
            inx_sg(5,ics_sg)=inx_sg(1,ics_sg)+lcontr
c           ....................
            ips_sg=ips_sg+lcontr
            icf_sg=icf_sg+isize
c           ....................
c type of shells
            inx_sg(12,ics_sg)=inx(12,ics)
c           write(6,*)'ics=',ics,' type=',inx(12,ics),' ics_sg=',ics_sg
c length of shells
            inx_sg(3,ics_sg)=inx(3,ics)
c general contr.
            inx_sg(4,ics_sg)=0
c center of shells
            inx_sg(2,ics_sg)=inx(2,ics)
c make datbas:
            ii=0
            do ipf=ish_b,ish_e
               ii=ii+1
               ips=inx_sg(1,ics_sg) + ii
               datbas_sg(1,ips)=datbas(1,ipf)       ! exponent
               datbas_sg(2,ips)=datbas(igc+2,ipf)   ! coeff
c
               datbas_sg(11,ips)=datbas(11,ipf)
               datbas_sg(12,ips)=datbas(12,ipf)
               datbas_sg(13,ips)=datbas(13,ipf)
            enddo
         enddo
      enddo
      do ics=1,ncs
         isize =inx(3,ics)
         igenc =inx(4,ics)
         igenc1=inx(4,ics)+1
         ish_b=inx(1,ics)+1
         ish_e=inx(5,ics)
         lcontr=inx(5,ics)-inx(1,ics)
         do igc=1,igenc
            ics_sg=ics_sg+1
            inx_sg(11,ics_sg)=icf_sg
            inx_sg(10,ics_sg)=inx_sg(11,ics_sg)+isize
c
            incr_gc=inx(11,ics)+igc*isize
            do icf=1,isize
               icf_gc=icf+incr_gc
               ireor(icf_gc)=icf_sg+icf
            enddo
c
            inx_sg(1,ics_sg)=ips_sg
            inx_sg(5,ics_sg)=inx_sg(1,ics_sg)+lcontr
c           ....................
            ips_sg=ips_sg+lcontr
            icf_sg=icf_sg+isize
c           ....................
            inx_sg(12,ics_sg)=inx(12,ics)
            inx_sg(3,ics_sg)=inx(3,ics)
            inx_sg(4,ics_sg)=0
            inx_sg(2,ics_sg)=inx(2,ics)
            ii=0
            do ipf=ish_b,ish_e
               ii=ii+1
               ips=inx_sg(1,ics_sg) + ii
               datbas_sg(1,ips)=datbas(1,ipf)       ! exponent
               datbas_sg(2,ips)=datbas(igc+2,ipf)   ! coeff
c
               datbas_sg(11,ips)=datbas(11,ipf)
               datbas_sg(12,ips)=datbas(12,ipf)
               datbas_sg(13,ips)=datbas(13,ipf)
            enddo
         enddo
      enddo
c     ncs_sg=ics_sg                       ! seg. contr. shells
c     nsh_sg=ips_sg                       ! seg. primi. shells
c     ncf_sg=icf_sg
c----------------------------------------------------------------
c     write(6,*) 'from GC to not ordered SG functions'
c     ncf=inx(10,ncs)
c     write(6,*) (ireor(i),i=1,ncf)
c
      end
c======================================================================
      subroutine trans_sg_2_gc(bl,inx,ictr,lshell)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      dimension bl(*)
      dimension inx(12,*)
c-----------------------------------------------------------
      if(lshell.eq.0) return
c-----------------------------------------------------------
c get back the GC-basis set info seved in trans_gc_2_sg
c
      call getival('ncs_gcsh ',ncs)
      call getival('ibas_gcsh',ibas)
      call getival('nsh_gcsh',nsh)
      call getival('ictr_gcsh',ictr)
c
c replace SG-shell by GC-shell basis set info :
c
      call setival('ncs ',ncs)
      call setival('ibas',ibas)
      call setival('nsh',nsh)
      call setival('ictr',ictr)
c
c
c-----------------------------------------------------------
c setup symmetry info
c
      call symfunc(.false.)
c-----------------------------------------------------------
c release memory alloctated in trans_gc_2_sg
c
c     call retint(1)
      call retmem(1)
c-----------------------------------------------------------
c
      end
c======================================================================
      subroutine  mapf_gc_sg_ord(ncf,ireor1,ireor2,ireor )
c     inp: ireor1 - from GC to NOT ordered SG
c     inp: ireor2 - from NOT ordered SG to ordered SG
c     out: ireor  - from GC to ordered SG
      dimension ireor1(ncf),ireor2(ncf),ireor(ncf)
c
      do icf_gc=1,ncf
         icf_sgn=ireor1(icf_gc)
         icf_sgr=ireor2(icf_sgn)
         ireor(icf_gc)=icf_sgr
      enddo
c
c     write(6,*) 'from GC to ordered SG functions'
c     write(6,*) (ireor(i),i=1,ncf)
      end
c======================================================================
      SUBROUTINE AFM_External(inp,NAtoms,XNuc,GAFM,GC)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Reads in and adds an external force to the gradient vector.
C  This is principally to simulate atomic force microscopy.
C  Applied external forces are directed along a line joining
C  any two specified atoms, either as "push" (to push the two
C  atoms together) or "pull" (to pull them apart).
C
C  ARGUMENTS
C
C  inp     -  unit number of input file
C  NAtoms  -  number of real atoms
C  XNuc    -  Cartesian coordinates (Old Texas Format)
C  GAFM    -  array for external forces
C  GC      -  array for existing gradient vector
C             (will be updated on exit)
C
      DIMENSION XNuc(5,NAtoms),GAFM(3,NAtoms),GC(3,NAtoms)
      REAL*8 V(3)
      Character*80 Char,Direction
      Logical jubel
      Parameter (Zero=0.0d0)
C
C
C  Are there any external forces?
C  These should be embedded in the input file following the FORCe command
C
C  $force
C  atom1   atom2   force(au)   PUSH/PULL
C  atom1   atom2   force(au)   PUSH/PULL
C   "       "       "           "
C  $endforce
C
      READ(inp,900,END=96) Char
      If(Char(1:1).NE.'$') Then
        BACKSPACE inp
  5     RETURN
      EndIf
C
C  We have something to read
C
      call RdGRAD(NAtoms,GC,'save')      ! recover the gradient
      CALL ZeroIT(GAFM,3*NAtoms)
c
 10   CONTINUE
      READ(inp,900,ERR=97) Char
      If(Char(1:4).EQ.'$end') GO TO 95
      READ(Char,*,ERR=97) IAt1,IAt2,force,Direction
C
C  remove ALL blanks from Character string
C
      call leadblan2(Direction,80,len1)
      call rmblan(Direction,len1,len)
      call uppercase(Direction,len)
c
      V(1) = XNuc(2,IAt2) - XNuc(2,IAt1)
      V(2) = XNuc(3,IAt2) - XNuc(3,IAt1)
      V(3) = XNuc(4,IAt2) - XNuc(4,IAt1)
c
      call nom(V,3,jubel,1.0d-10)      ! normalize force direction
c
      If(Direction(1:4).EQ.'PUSH') Then
        GAFM(1,IAt1) =  V(1)*force
        GAFM(2,IAt1) =  V(2)*force
        GAFM(3,IAt1) =  V(3)*force
        GAFM(1,IAt2) = -V(1)*force
        GAFM(2,IAt2) = -V(2)*force
        GAFM(3,IAt2) = -V(3)*force
      Else                               ! assumed PULL
        GAFM(1,IAt1) = -V(1)*force
        GAFM(2,IAt1) = -V(2)*force
        GAFM(3,IAt1) = -V(3)*force
        GAFM(1,IAt2) =  V(1)*force
        GAFM(2,IAt2) =  V(2)*force
        GAFM(3,IAt2) =  V(3)*force
      EndIf
c
      GO TO 10
 95   CONTINUE
C
C  All external force input completed
C  Print out external force, at same time add to gradient vector
C
      IOut = igetival('iout')
      WRITE(IOut,1000)
      DO I=1,NAtoms
      If(GAFM(1,I).NE.Zero.OR.GAFM(2,I).NE.Zero
     $    .OR.GAFM(3,I).NE.Zero) Then
         WRITE(IOut,1100) I,GAFM(1,I),GAFM(2,I),GAFM(3,I)
         GC(1,I) = GC(1,I)+GAFM(1,I)
         GC(2,I) = GC(2,I)+GAFM(2,I)
         GC(3,I) = GC(3,I)+GAFM(3,I)
      EndIf
      EndDO
      WRITE(IOut,*)
C
C  write gradient to <grad> file
C
      CALL WrGRAD(NAtoms,GC)
 96   CONTINUE
      RETURN
C  ..............................................
C    ERROR SECTION
C
 97   CONTINUE          ! missing data of line
      Call nerror(39,'FORCE module',
     $  'Problems Reading External Force - Check your input!',
     $   0,0)
c
  900 Format(A80)
 1000 Format(/,' === The Following External Forces Will Be Applied ===')
 1100 Format('  Atom: ',I4,2X,3F15.8)
c
      END

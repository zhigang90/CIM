C================================================================
      subroutine newtransloc
C       Local MP2  program.
C
C       This program performs the integral transformation for
C       Local MP2, at the end of loctrans MP2iter is called to
C         carry out       the Local MP2 iterations.
C
C       External calls (excluding utility routines)
C       For details see comments in the appropriate subroutine.
C     SORTCOEFF - makes table istab(ncf,nval), this table is used
C       of decreasing magnitude.
C     EVALC - evaluates and prints out the magnitudes of the
C       MO coefficients.  Should probably be removed in the final
C       version.
C     PLIST - generates the pair lists, and classifies the pairs
C       as strong weak and negligible. See comments in Plist for
C       description of the various pair-lists.
C     PRINTLIST - prints the pairlists
C     DMAXIJ - constructs matrix 'dsmx' . This is used for pre-sceening
C       of integrals.
C     MKTAB - makes table istab1(jmax,nval,ncf).  For details
C       see subroutine mktab
C       The following six subroutines are integral and related
C       routines written  by Krys. Wolinski
C     PTWOINT  - preparation for integral calculation (comments below)
C     INT_SCHWARZ - Schwarz integrals
C     DMAX_FING_BB_MP2 - setup for prescreening of integrals
C     INT_LMP2 - calculates a batch of integrals needed for MP2
C       (M,my,L,sig) M, and L Shells , all my and sig
C     INIT_INFO and TERM_INFO - statistics results in file 'timings'
C     TRA1SHEL - transforms one block of iintegrals, removes
C       all negligibel transformed integrals , converts to integers
C       and writes the trasnformed integrals out to scratch disk
C       xxx.htr.  This is where most of the work is done.
C     PHASE1 binsort phase 1 for localized orbitals
C     PHASE2 binsort phase 2 for localized orbitals       
C       in this subroutine the final internal exchange matrices are
C       projected and reduced to local dimensions and then written out
C       on external storage.
C     MP2ITER - performs MP2 iterations and calculates the MP2 energy.
C              
C       Svein Saebo Fayetteville AR Summer 1999
C

      use memory

      implicit real*8 (a-h,o-z)
      dimension xintxx(9)
c     common /big/bl(30000)
c     common /intbl/maxsh,inx(100)
      parameter(nopt=20)
      parameter (zero=0.0d0)
      character*256 chopv(nopt),scrfile,filnam,filename(5),filbin(5)
      character*256 filname1,filname0
      character opnames*4,scftype*11
      character*3 ch3
      logical nofr,exst,full,relb,tcou,lalb,restar,redu,weak,coro
      logical LastSymPair,remf,smal
      logical dualbasis
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      dimension opnames(nopt)
      dimension xnmo(2)
C   added by ckw for smal option
C  see cmp2.f
      dimension Ipass(2,28),Kpass(2,28)
      data opnames  /'prin','nofr','pmij','orbs','thre','core','epsi',
     *               'dual','full','relb','tcou','lalb','clim','rest',
     *               'redu','nowe','coro','nrem','smal','eqor'/
C       set ioptyp according the the type of input option:
      data ioptyp /1,0,2,2,13,11,12,0,0,2,0,0,11,1,0,0,0,0,0,11/
C      0=logical option, 1=a single integer, 2=two integers,
C      3=three integers
C      11=a single real, 12=two reals, 13=three reals,
C      21=a character option
C----------------------------------------------------------------------
      parameter(gigabyte=2.147d9,giga=2.683d8)
C   2gb in bytes and in real*8 words
      call secund(tt0)
      call elapsec(elaps0)
c
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
      write(iout,*)'=================================================='
      write(iout,*)'            The Local MP2 Module  '
      write(iout,*)'=================================================='
      write(iout,*)' '
      iprint=1
      write(iout,*) '+++++++++++++ JOB PARAMETERS ++++++++++++++++++++'
C       Process input options:
C--------------------------------------------------------------------
C      ideally no input options should be neccessary for LMP2
C      most options below are for testing and debugging purposes
C      see below for description of the various options
C--------------------------------------------------------------------
C
      call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
C     'PRIN' sets print level, (3 or higher produces a lot of output)
      if(ifound(1).gt.0) then
        iprint=iopv(1,1)
      end if
C        'NOFR' no frozen core
      write(6,*) 'printlevel =',iprint
      nofr=.false.
      if(ifound(2).gt.0) nofr=.true.
C        'PMIJ' print matrix; first number is 1 (AO Xchange) or 2
C        (both AO& MO) second number: pair to be printed
      ipmij=0
      ipairpr=0
      if(ifound(3).gt.0) then
        ipmij=iopv(1,3)
        ipairpr=iopv(2,3)
      end if
C     ORBS - see below
C       THRESHOLDS
C       first :
C         integral threshold and exit threshold in
C         integral transformation
C       second:
C         threshold for neglecting unprojected kij elements
C         third:
C       threshold for neglecting elements (and reducing local
C       domians) applied after projection
C       give thresholds as PH or as a small number with decimal point
      if(ifound(5).gt.0) then
       if(ropv(1,5).gt.1)then
         thresh=10.0d0**(-ropv(1,5))
       else
       thresh=ropv(1,5)
       endif
       thres1=thresh
         if(ropv(2,5).gt.0) then
       if(ropv(2,5).gt.1)then
         thres2=10.0d0**(-ropv(2,5))
       else
       thres2=ropv(2,5)
       endif
         else
            thres2=1.0d-7
         endif
         if(ropv(3,5).gt.0) then
       if(ropv(3,5).gt.1)then
         thres3=10.0d0**(-ropv(3,5))
       else
       thres3=ropv(3,5)
       endif
         else
            thres3=1.0d-5
         endif
      else
         thresh=5.0d-9
         thres1=thresh
         thres2=1.0d-7       
         thres3=1.0d-5
      end if
C       EPSI  thresholds for negelect and weak pairs
C       give as PH: argument 1 epsnp argument 2 epswp
C       eij estimated pair energy based on the volumes and the
C       distance between orbitals i and j
C       if eij < epsnp  the pair is dropped from the calculation
C       if epsnp < eij < epswp  the pair energy is calculated using
C       a multipole expansion (not ready yet June 99)
      if(ifound(7).gt.0) then
       if(ropv(1,7).gt.1)then
         epsnp=10.0d0**(-ropv(1,7))
       else
       epsnp=ropv(1,7)
       endif
         if(ropv(2,7).gt.0) then
       if(ropv(2,7).gt.1)then
         epswp=10.0d0**(-ropv(2,7))
       else
       epswp=ropv(2,7)
       endif
         else
            epswp=1.0d-5
         endif
      else
         epsnp=1.0d-6
         epswp=1.0d-5
      end if
C       these thresholds seem to work well when weak
C       pairs are calculated with multi-pole expansion
C       Rij>7 Bohrs
ckw
      dualbasis=.false.
      if(ifound(8).gt.0) then
         call nerror(1,'newtransloc',
     *   'The Local MP2 module is NOT ready for dual basis set yet',0,0)
      endif
      full=.false.
      if(ifound(9).gt.0) then
       write(6,*) ' "FULL" calculation local basis not used'
      full=.true.
c       mainly for testing purposes
      endif
       relb=.false.
       if(ifound(10).gt.0) then
       relb=.true.
C       give local domain for a MO as input, mainly for testing
C       purposes Input 2 integers; first valence orbital number,
C       second number of AO in domain
       ilbmo=iopv(1,10)
       idimo=iopv(2,10)
       write(6,*) 'Local domain for orbital: ',ilbmo
       write(6,*) 'given as input, dimension :', idimo
       call getmem(2+idimo,lbsin)   !  MMOK
       bl(lbsin)=ilbmo
       bl(lbsin+1)=idimo
       read(inp,*) (bl(lbsin+1+ile),ile=1,idimo)
       write(6,*) (bl(lbsin+1+ile),ile=1,idimo)
       else
       call getmem(1,lbsin)   !  MMOK
       endif
       tcou=.false.
       if(ifound(11).gt.0) then
       tcou=.true.
      write(6,*)'Multiplication with S inside the loop'
       else
      write(6,*)'Multiplication with S outside the loop'
       endif
       lalb=.false.
       if(ifound(12).gt.0) then
       lalb=.true.
      write(6,*)'large local basis will be used'
       endif
       clim=3.0d-6
       if(ifound(13).gt.0)then
       clim=ropv(1,13)
       endif
       write(6,*) 'convergence criterion in MP2 iterations: ',clim
C       restart
       restar=.false.
       if(ifound(14).gt.0) then
       restar=.true.
       irst1=iopv(1,14)
       if(irst1.le.1) then
       irst1=1
       write(6,*) 'Restart from Internal exchange matrices'
       write(6,*) 'File xxxx.kij required'
       endif
       if(irst1.eq.2) then
       write(6,*) 'Restart from bins file xxx.bins required'
       endif
       if(irst1.ge.3) then
       write(6,*) 'Restart from halftransformed integrals'
       write(6,*) 'File xxxx.hftr required'
       endif
       endif
       redu=.false.
       if(ifound(15).gt.0) then
       write(6,*) 'unprojected formalism will be used'
       redu=.true.
       endif
       weak=.true.
       if(ifound(16).gt.0) then
      write(6,*) 'weak pairs will not be included in this calculation'
       weak=.false.
       else
       write(6,*) 'weak pairs will  be included'
       endif
       if(ifound(17).gt.0) then
       write(6,*) 'amplitudes stored in core'
       coro=.true.
       else
       coro=.false.
       write(6,*) 'amplitudes stored on disk'
       endif
       remf=.true.
       if(ifound(18).gt.0) then
       write(6,*) ' Redundant functions will not be removed'
       remf=.false.
       endif
        natoms=igetival('na')
        if(natoms.gt.30) then
       smal=.false.
       else
       smal=.true.
       endif
       smal=.false.
       if(ifound(19).gt.0) then
       write(6,*) ' SMAL option for integral storage forced'
       smal=.true.
       endif
          syeqv=1.0d-8
          if(ifound(20).gt.0) then
          if(ropv(1,20).gt.1) then
      syeqv=10.d0**(-ropv(1,20))
      else
      syeqv=ropv(1,20)
      endif
      endif
      write(6,*) 'Threshold for symmetry equivalent orbitals:', syeqv
C
       write(6,*) 'MP2 integral threshold: ', thresh
       write(6,*) 'MP2 exit threshold: ', thres1
       write(6,*) 'MP2 truncation threshold(unprojected): ', thres2
       write(6,*) 'MP2 truncation threshold(projected): ', thres3
        write(6,*) 'Pair neglect threshold: ',epsnp
        write(6,*) 'Weak pair threshold: ',epswp
C
      call tstrval('core',iyes)
      if(iyes.eq.1) then
        core=rgetrval('core')
      else
        core=-3.0d0
      end if
      if(ifound(6).gt.0)  core=ropv(1,6)
C-------------- input options processed -----------------
C        preset some variables to zero
      call matreset
      idft=0
      ax=zero
      trmz=zero
      nrad=0
      nang=0
      lrad=0
      lang=0
      Iradq=0
      NBatch=0
      nwo=0
      xn1=zero
      xn2=zero
      xn3=zero
      xn4=zero
      xn5=zero
C       zeroout counters
      xmult=zero
      xnumint=zero
      actint=zero
      transin2=zero
      iktot=0
C
      tss21=zero
C        Get values for  some common variables
      ncs=igetival('ncs')
      ncf=igetival('ncf')
C     nmo=igetival('nmo')
       np4=4
       call sudat(np4,'nocc_rhf',ni)
       if(ni.gt.0) then
         call rea(xnmo,2,np4,'nocc_rhf')
       else
         call restart(np4,0)
         call rea(xnmo,2,np4,'nocc_rhf')
       end if
       nmo=xnmo(1)
       erhf=xnmo(2)
C     SYMMETRY related stuff
      nsym=igetival('nsym')
      call getival('SymFunPr',ifp)
      call getival('SymShPr',ifp1)
C    READ Localized  orbitals and orbital energies
C       we'll only need localized valence orbitals
C       and orbital energies.  The orbital energies are only used to
C       determine the number of core orbitals and will be
C       immediately removed.
C       Coulson energies for the localized orbitals will be needed
C       later in the MP2 iterations , and will be read in then
       call matdef('loca','r',ncf,nmo)
       call matread('loca',np4,'loca_rhf')
       call matdef('epsi','d',nmo,nmo)
       call matread('epsi',np4,'eval_rhf')
       iepsi=mataddr('epsi')
C       determine the number of core orbitals, these will
C       not be correlated.
C       default: orbitals with orbital energies smaller than
C       -3 are core.
      ecore=-3.0d0
      ncore=ncores(ncf,bl(iepsi),ecore)
      nfirst=ncore+1
      nlast=nmo
      nval=nmo-ncore
      if(nofr) then
        nfirst=1
        nlast=nmo
        nval=nmo
      end if
      if(ifound(4).gt.0) then
        nfirst=iopv(1,4)
        nlast=iopv(2,4)
        nval=nlast-nfirst+1
      end if
C       remove orbital energies
       write(6,*) 'Number of contracted basis functions: ',ncf
       write(6,*) 'Number of correlated orbitals:        ',nval
      call matrem('epsi')
      call matsub('occu','loca',nfirst,nlast)
C
C      iloca=mataddr('loca')
C       call memo1_int(ncore,lcortb)
C          call procore(bl(iloca),ncore,ncf,bl(lcortb),iprint)
C          call retmem(1)
C
      jloca=mataddr('occu')
C        Reserve memory for tables, these should be kept for the
C       duration of the calculation.
C       isort will contain sort matrix (see subroutine Sortcoeff)
C       basis function symmetry pairs are in inx(ifp)
C       read in overlap matrix in triangular form
       call secund(tx1)
       if(nsym.gt.0) then
       call getint(nval*nsym,iorrel)
       call matdef('ovl','q',ncf,ncf)
       issad=mataddr('ovl')
       call matdef('smat','s',ncf,ncf)
       np1=1
       call sudat(np1,'s matrix',ni)
       if(ni.gt.0) then
       call matread('smat',np1,'s matrix')
       else
         call restart(np1,0)
       call matread('smat',np1,'s matrix')
       end if
       call matcopy('smat','ovl')
       call matrem('smat')
C       generated tables for symmetry handling
C       construct table isorb(nval,nsym)
C       with symmetryrelations between valence orbitals
C       VSC used in symeq which is used by orbpr
       call matdef('VSC','d',ncf,ncf)
       call orbpr(ncf,nval,nsym,bl(jloca),bl(issad),bl(ifp),
     1 bl(iorrel),iprint,syeqv)
       call matrem('VSC')
       call matrem('ovl')
       else
       call getmem(1,iorrel)   !  MMOK
       endif
C       end symmetry setup
        call secund(tsy)
        ttsy=(tsy-tx1)/60.0d0
        write(6,*) 'symmetry setup',ttsy
       call getint(nval*ncf,isort)
C
      call matdef('tloca','r',nval,ncf)
      iloca=mataddr('tloca')
      call matpose2('occu','tloca','n')
c  Sort the localized orbital coefficients by magnitude
C     call matprint('tloca',6)
      call SortCoeff(ncf,nval,bl(iloca),bl(isort))
      call matpose2('occu','tloca','n')
C     check and print the magnitude of the coeffcients
       if (iprint.ge.2) then
       call EvalC(bl(jloca),nval,ncf)
       endif
       ixma=idamax(nval*ncf,bl(jloca),1)
       xma=abs(bl(jloca-1+ixma))
C        find the largest MO coefficient
      if(iprint.ge.2) write(6,*)'max MO coefficient: ',xma
C       reserve memory for pair-lists (Keep)
        npairs=nval*(nval+1)/2
       nvasq=nval**2
       call getint(npairs,iplis)
       call getint(npairs,ilis)
       call getint(npairs,jlis)
       call getint(nvasq,invli)
c      call getmem(nvasq/8+1,nstrong)
       call getbyte(nvasq,nstrong)
C
C.....Read centers and radii of localized orbitals to
C.....contruct pair list
        call matdef('xmo','d',nmo,nmo)
        call matdef('ymo','d',nmo,nmo)
        call matdef('zmo','d',nmo,nmo)
        call matdef('rro','d',nmo,nmo)
       call matread('xmo',np4,'xcor_loc')
       call matread('ymo',np4,'ycor_loc')
       call matread('zmo',np4,'zcor_loc')
       call matread('rro',np4,'rr_loc')
       ixmo=mataddr('xmo')+ncore
       iymo=mataddr('ymo')+ncore
       izmo=mataddr('zmo')+ncore
       irra=mataddr('rro')+ncore
        call secund(tss1)
        twq=(tss1-tsy)/60.0d0
        write(6,*) 'sorting', twq
      call Plist(bl(ixmo),bl(iymo),bl(izmo),bl(irra),bl(iplis),
     1          bl(ilis),bl(jlis),bl(invli),bl(nstrong),nval,
     2          jcmax,mstrong,epsnp,epswp,nud,
     3          wken,weak,nnegl,nweak,iprint)
        call secund(tss2)
        tss=(tss2-tss1)/60.0d0
        write(6,*) 'Plist: ',tss
C     Print pair-lists:
       mprin=mstrong
       if(weak) mprin=mprin+nweak
      call printlist(bl(iplis),bl(ilis),bl(jlis),bl(invli),nval,iprint,
     1 mprin)
      write(iout,*)'++++++++++++++++++++++++++++++++++++++++++++++++++'
       if(weak.and.(nweak.eq.0)) then
       weak=.false.
       write(6,*)' There are no weak pairs in this system'
       write(6,*)' weak-option turnd off'
       endif
C       remove stuff not needed
       call matrem('rro')
       call matrem('zmo')
       call matrem('ymo')
       call matrem('xmo')
       mpairs=mstrong
       if(weak) mpairs=mstrong+nweak
       call getint(mpairs,iunik)
       call getint(2*mpairs,iekvi)
C       set up a table over weak pairs
       if(weak) then
c      call getmem((nval**2)/4+1,iwkli)
       call getint_2(nval**2,iwkli)
       call mkwkli(nval,bl(iplis),bl(iwkli))
       endif
C
       call equip(bl(iekvi),bl(iunik),nval,bl(iorrel),nsym,
     1 bl(nstrong),bl(invli),bl(iplis),munik,iprint,bl(iwkli),weak,
     2 mstrong,munk,munndia)
       muweak=munk-munik
       if(iprint.ge.2) then
          write(6,*) 'number of symmetry unique pairs       : ',munk
          write(6,*) 'number of symmetry-unique strong pairs: ',munik
       write(6,*) 'number of symmetry unique n.d.   pairs: ',munndia
       write(6,*) 'number of symmetry-unique weak pairs  : ',muweak
       endif
       if(nsym.eq.0) then
       frac=1.0d0
       else
       frac=float(munndia)/float(nud)
       endif
       write(6,*) 'frac=',frac
       call getint(2*nval,iunor)
       call unmos(bl(iunor),nval,nsym,bl(iorrel),iprint)
C
C       call chectab(bl(iekvi),bl(iunik),bl(iunor),
C    1    bl(nstrong),bl(iplis),bl(invli),nval,nsym,iprint)
C
C       construct matrix 'dsmx'
C       this will be used for prescreening of integrals
       call matdef('dsmx','q',ncs,ncs)
       call getmem(ncf**2,idmxx)  ! MMOK
       idmss=mataddr('dsmx')
C         call dmaxIJ(bl(idmxx),bl(idmss),ncf,ncs,nval,bl(nstrong),
C    1                inx,bl(jloca),densmax,iprint)
        call getmem(ncf*nval,itemp) ! MMOK
          call dmaxIJ_new(bl(idmxx),bl(idmss),ncf,ncs,nval,bl(nstrong),
     1                inx,bl(jloca),densmax,iprint,bl(itemp))
          call retmem(1)
       call retmem(1)
c
c         densmax is the maximum "screening" density element
c         in the integral program only those int. are calculated
c         which   Dmax*X >= Thres
c         Here only those integrals are used which pass
c                 Cmox*Cmox*X>=Thres1
c         Thus, threshold for MP2 inregrals (thres) should be adjusted
c         to thres1 , Dmax and Cmox :
c
c         thres/Dmax = thres1/Cmox2   (Cmox2=Cmox*Cmox)
c         thres = thres1* Dmax/Cmox2
ckw................................................................
       if(iprint.ge.2) then
       write(6,*)'max screening dens: ',densmax
       write(6,*)'max MO coefficient: ',xma
       write(6,*)'re-def.integ.thres=thres1*(maxdens/(MOmaxcoef)2)'
       write(6,*)'   original integ.thres.=',thresh
       endif
          xma2=xma*xma
          thresx=thres1*densmax/xma2
       thresh=thres1
       if(iprint.ge.2) Then
       write(6,*)'  redefined integ.thres.=',thresx
       write(6,*) 'thresh=thres1= ',thresh
       endif
ckw................................................................
C       generate  table istab1 (see comments in subroutine mktab
C       for details)
        jjmax=jcmax+1
c     call getmem((ncf*nval*jjmax)/4+1,itab1)
      call getint_2(ncf*nval*jjmax,itab1)
      call secund(tmkt)
      tsy2=(tmkt-tss1)/60.0d0
      write(6,*) 'SYM2',tsy2
      call mktab(bl(isort),bl(itab1),ncf,nval,jjmax,bl(nstrong)
     1       ,iprint)
       call secund(tx2)
        tmktab=(tx2-tmkt)/60.0d0
        write(6,*) 'mktab', tmktab
C       if(iprint.ge.2) then
        ttab=(tx2-tx1)/60.0d0
       write(6,*) 'time for generating tables: ', ttab,' min.'
C       endif
C   reserve emory for local domains
C    The following two arrays will contain final local domais +
C   dimensions and must be saved
c     call getmem((ncf*nval)/4+1,lbasp)
      call getint_2(ncf*nval,lbasp)
c     call getmem(nval/4+1,idimp)
      call getint_2(nval,idimp)
c     call getmem(nval/4+1,jdimp)
      call getint_2(nval,jdimp)
c     call getmem((ncf*nud)/4+1,lbasij)
      call getint_2(ncf*nud,lbasij)
c     call getmem(nud/4+1,idiij)
      call getint_2(nud,idiij)
c     call getmem(nud/4+1,jdiij)
      call getint_2(nud,jdiij)
C       reserve memory for table with pointers to T's
C       must be kept for the duration of the calculation
       call getint(mstrong+nweak+1,itsts)
c  open scratch file to store half-transformed integrals
C       matrix definitions for matrices used in mp2iter, moved here
C       so temporary storage used for binsort can be removed
       if(ncore.gt.0) then
       lencor=ncore/2+1
       call matdef('tab','d',lencor,lencor)
       else
       call matdef('tab','d',1,1)
       endif
       call matdef('ovl','q',ncf,ncf)
       call matdef('fck','q',ncf,ncf)
       call matdef ('proj','q',ncf,ncf)
        call matdef('floc','q',nval,nval)
       if(redu) call matdef('sds','q',ncf,ncf)
C       end matrix definitions
C
c      call getmem(npairs/4+1,invpd)
       call getint_2(npairs,invpd)
       call invndp(bl(nstrong),nval,bl(invpd))
C
C       put down a memory mark
C       all tables, orbitals, orbital energies before this mark
       call mmark
C
      ndisk=17
      iunit=ndisk
      call getchval('scrf',scrfile)
      call rmblan(scrfile,256,len)
      exst=.false.
      inquire(iunit,exist=exst)
      if(exst) close(iunit,status='delete')
c  the length of a block is given in bytes.
      lint=igetival('ints')  ! # of integer word in a double-precision
      idoub=igetival('double')
      ints=idoub/lint   ! size of an integer in bytes
C
C       in this implementation individual transformed integrals
C       are written out with a total of 4 indices (i,j,my,lam)
C       For each integral a total of 3 integer*4 words are written:
C       1)my and lam packed into one integer
C       2)i and j packed into one integer
C       3)Integral as an integer
C       see subroutine Remzero and Skriv for further details
C       we use fixed record length and direct access?
C        (the file is read sequentially though)
C
       lrecw=4096
          lrec=lrecw*ints
C       fixed record length
C
      filname0=scrfile(1:len)//'.htr'
          len1=len+4
      open(unit=ndisk,file=filname0(1:len1),form='unformatted'
     1 ,access='direct',recl=lrec)
          nrcpf=gigabyte/dble(lrec)
      Write(6,*)'record length: ',lrec, 'number of records per file: '
     1, nrcpf
          write(6,*) 'filelength:', nrcpf*dble(lrec)
C
C open 4 more files for now, I don't know how many I'll need
C
      Nstripe=4
      Do ihtr=1,NStripe
      write(ch3,'(I3)') ihtr
      ch3(1:1) = '.'
      If(ihtr.lt.10) ch3(2:2) = '0'
      KUnit = ndisk+ihtr
      filname1 = filname0(1:len1)//ch3
      len = len1+3
      OPEN (UNIT=KUnit,FILE=filname1(1:len),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec)
      EndDo
C
c  space for ao integrals (temp)
       call matdef('xmat','q',ncf,ncf)
       ixadr=mataddr('xmat')
C  space for buffer with Kmylam matrices (temp)
      call getint(lrecw,ibufa)
c space for one half-transformed  exchange operator (temp)
      call matdef('halftra','q',nval,nval)
      ihalftra=mataddr('halftra')
C  prepare for integral calculation ************************
      thint=thresh
      istat=igetival('istat')
      iforwhat=5
c..............................
      call ptwoint1(nfock,  thresh, thint, iforwhat,idft,
     *             ax,     nrad,   nang,   lrad,   lang,
     *             Iradq,  NBatch, .true., nmo,    0,
     *             scftype,xintxx, nblocks,maxbuffer,maxlabels)
c
c---------------------------------------------------
c reserve memory for labels ..
c
           call getmem(maxlabels,labels) ! MMOK
c
c---------------------------------------------------
c Output from the call above :
c nblocks   - total number of blocks of c.s.q.
c maxbuffer - maximum size of the integral buffer
c maxlabels - maximum size of 3 integer arrays for integral labels :
c             lsize(4), lindex(4,ngcd,nbls),lgenc(nbls)
c             lsize(1)=ics_len, lsize(2)=jcs_len etc.
c             lindex(1,iqu,ijkl)=icf_start,
c             lindex(2,iqu,ijkl)=jcf_start etc.
c scftype   - scf mode to be run
c xintxx(9)  contains number of integrals more expensive than :
c 0.0%, 0.1%, 0.5%, 1%, 5%, 10%, 25%, 50%, 75% of maximum price
c calculate Schwarz integrals (ij|ij) :
c
C     call getmem(ncs*ncs,ischwarz)
C     call setival('schwarz',ischwarz)
      ictr=igetival('ictr')
C     call int_schwarz(bl,inx(ictr),bl(ischwarz),bl(labels),'sch1')
c---------------------------------------------------
c allocate memory for an array mapping contr.func. ti contr.shells :
      call getint(ncf,mapf2s)
      call make_mapf2s(bl(mapf2s),ncs,ncf,bl(ictr))
c---------------------------------------------------
      call secund(tt)
      if(iprint.ge.3) write(iout,*) 'Startup time=',tt-tt0
      tinteg=zero
      ttrans=zero
      elapint=zero
      elaptrans=zero
      ncfncf=ncf*ncf
      call secund(tt3)
      call elapsec(telap3)
c
      call init_info('fock')
C                                       ***************************
c  zero out record counter
      nrec=0
        call symmoff
       intsize=igetival('ints')
C
       call getint(nval,icol2)
       call getint(nval*jjmax,jfun3)
c      call getmem(nval/4+1,izero)
       call getint_2(nval,izero)
       call mkjfun3(nval,jjmax,bl(nstrong),bl(jfun3))
       call getint(ncf,lzero)
C       the following array  for indices (1,ncf) integer*2
C       should be sufficient
       call getint(ncf,irow)
       call getint(ncf,icol)
       call getint(ncf,irow1)
       call getint(ncf,icol1)
C
      lcore=igetival('lcore')
      ioffset=igetival('ioffset')
      kcore=lcore-ioffset
C
C     the following for 'BIG' option, SS
C
      if(.not.smal) then
        call getint(ncf,if2cr)
        call getint(ncf,if2cc)
        call getmem(1,last_entry)   ! MMOK
        mem=kcore-last_entry+ioffset-6*ncf*ncf-2000000
        memx=mem/9
        mem=memx*4
        lenindj=memx*4
        call getmem(memx,indlj) ! used as int*2, ?? MM
        write(6,*)' BIG option used for integrals'
        write(iout,*) ' memory assigned to the job:     ',lcore-ioffset
        write(iout,*) ' memory available for integrals: ',mem, '*2'
        lmp2_size = mem
        call getmem(lmp2_size,lmp2int) ! MMOK
        write(iout,*) ' lmp2int start:',lmp2int-ioffset,mem,' long'
        write(iout,*) ' end integral storage',lmp2int-ioffset+2*mem
       indmax=0
       istmax=0
        lenmax=lenindj/6
      endif
      ncf2=ncf*ncf
      call secund(ttsta)
      tstart=(ttsta-tt0)/60.0d0
      write(6,*) 'Startuptime:',tstart
C
C    Transformation Starts HERE...........
C       Triagular loop over shells
      do ics=1,ncs
         call get_shell_size(bl(ictr),ics,ics_size)
         lmp2_siz1=ncf2*ics_size
         do kcs=1,ics
           if(LastSymPair(ics,kcs,nsym,bl(ifp1),inegl,iret)) cycle
c of each (ics,kcs) pair, calculate only the last one
c (ics',kcs')>(ics,kcs) if ics'>ics or ics'=ics and kcs'>kcs
           call secund(tt2)
           call elapsec(telap2)
           ttrans=ttrans+tt2-tt3
           elaptrans=elaptrans+telap2-telap3
C
           if(smal) then
         call get_shell_size(bl(ictr),kcs,kcs_size)
           lmp2_size=lmp2_siz1*kcs_size
c
c check if a size of the lmp2 integral buffer is not too big
c if it is  then split over kcs ( and possibly ics)
c
c
          call check_size1(lmp2_size,ncf,nval,nval,ntimes)
c
          call set_passes(bl(ictr),ics,kcs,ics_size,kcs_size,
     *                      ntimes,Ipass,Kpass,Itimes,Ktimes)
c
          do itime=1,itimes
             icf1=ipass(1,itime)
             icf2=ipass(2,itime)
             iatonce=icf2-icf1+1
          do ktime=1,ktimes
             kcf1=kpass(1,ktime)
             kcf2=kpass(2,ktime)
             katonce=kcf2-kcf1+1
c
             lmp2_size=iatonce*katonce*ncf2
c
            call getmem(lmp2_size,lmp2int) ! MMOK
            call mmark
            call int_lmp2(bl, bl(ictr), thresh, ics, icf1,
     1                    icf2, kcs, kcf1, kcf2, bl(mapf2s),
     2                    bl(idmss),iprnt,bl(lmp2int),nintotal,nrow,
     3                    ncol, bl(irow),bl(icol),bl(lzero))
            call retmark
          if(nintotal.eq.0) then
            nskipped=nskipped+1
           call retmem(1)
           cycle
          else
            iktot=iktot+1
            endif
           call secund(tt3)
           call elapsec(telap3)
           tinteg=tinteg+tt3-tt2
           elapint=elapint+telap3-telap2
c......................................................................
C       evaluate the magnitudes of the integrals
C       this takes considerable time,use only for testing
C       call evalI(bl(lmp2int),lmp2_size,thresh,xn1,xn2,xn3,xn4,xn5)
c......................................................................
C      Find max value of integral in this block
Cs           iimax=idamax(lmp2_size,bl(lmp2int),1)
Cs          aimax=bl(lmp2int-1+iimax)*thresh
Cs       aimax=abs(aimax)
C           max integral in block * max mo coefficient = xmoxint
Cs          xmoxint=aimax*xma
C           write(6,*) 'max integral in block :',aimax
C           if (max integral) * (max coeff)**2 < thre1 cycle
Cs          if((xmoxint*xma).lt.thres1)then
Cs             call retmem(1)     !lmp2int
Cs              cycle
Cs           endif
c................................................................
C      transform one block, remove negligible contributions
C       convert to integers and write out on xxxx.htr
C................................................................
      call Tra1Shel(xma,ncf,ncs,nval,ics,
     *             icf1,icf2,kcf1,kcf2,
     1              kcs,ndisk,bl(ictr),bl(lmp2int),iprint,
     2              bl(ixadr),bl(ixadr),bl(ihalftra),nrec,bl(isort),
     3              bl(itab1),jjmax,ncore,bl(iplis),bl(invli),
     4              thresh,thres1,thres2,nwo,lrecw,
     5              ibufa,mult,bl(nstrong),nsym,bl(ifp),
     6              bl(iorrel),bl(iloca),trmz,transin2,nrow,
     7              ncol,bl(irow),bl(icol),bl(irow1),bl(icol1),
     8              bl(icol2),bl(izero),bl(jfun3),smal,nrcpf)
C
C       xmult is a counter of actual multiplications
C       in the inner loop and the number of time this loop is
C       entered. Should be removed in the final version.
C       mult and loop are incremented in tra1shel.
C     xmult=xmult+mult
C
      call retmem(1)              ! lmp2int
C
      enddo  ! loop over ktimes
      enddo  ! loop over itimes
C
        else   !if(smal)
C
       ind=0
       intstore=0
       call getmem(lmp2_size,lrestor) ! MMOK
      call int_lmp2b(bl,bl(ictr),thresh,ics,kcs,
     1              bl(mapf2s),bl(idmss),iprint,bl(lmp2int),nintotal,
     2              bl(icol),bl(irow),bl(if2cc),bl(if2cr),
     3              bl(indlj),bl(lrestor),nrow,ncol,ind,
     4              intstore,lmp2_size,lenmax)
       indmax=max0(indmax,ind)
       istmax=max0(istmax,intstore)
       ttint=ttint+nintotal
       call retmem(1)              !lrestor
C       if(ind.eq.0) write(6,*) 'ind',ind,intstore
C
           call secund(tt3)
           call elapsec(telap3)
           tinteg=tinteg+tt3-tt2
           elapint=elapint+telap3-telap2
       if(nintotal.eq.0.or.nrow.eq.0.or.ncol.eq.0) then
C       print *,'shell pair ',ics,kcs,' skipped'
       if(smal)call retmem(1)
       cycle
       else
       iktot=iktot+1
       actint=actint+nintotal
       endif
c......................................................................
C       evaluate the magnitudes of the integrals
C       this takes considerable time,use only for testing
C       call evalI(bl(lmp2int),lmp2_size,thresh,xn1,xn2,xn3,xn4,xn5)
c......................................................................
C
c................................................................
C      transform one block, remove negligible contributions
C       convert to integers and write out on xxxx.htr
C................................................................
      call getinfs(bl(ictr),ics,kcs,icf1,icf2,kcf1,kcf2)
      call Tra1Shel(xma,ncf,ncs,nval,ics,
     *             icf1,icf2,kcf1,kcf2,
     1              kcs,ndisk,bl(ictr),bl(lmp2int),iprint,
     2              bl(ixadr),bl(ixadr),bl(ihalftra),nrec,bl(isort),
     3              bl(itab1),jjmax,ncore,bl(iplis),bl(invli),
     4              thresh,thres1,thres2,nwo,lrecw,
     5              ibufa,mult,bl(nstrong),nsym,bl(ifp),
     6              bl(iorrel),bl(iloca),trmz,transin2,nrow,
     7              ncol,bl(irow),bl(icol),bl(irow1),bl(icol1),
     8              bl(icol2),bl(izero),bl(jfun3),smal,nrcpf)
C
C       xmult is a counter of actual multiplications
C       in the inner loop and the number of time this loop is
C       entered. Should be removed in the final version.
C       mult and loop are incremented in tra1shel.
C     xmult=xmult+mult
C
C
      endif    !if(smal)
         enddo                  ! over shells kcs
      enddo                      ! over shells ics
       if(.not.smal) call retmem(1)              !lmp2int
C     statistics results in file 'timings'
      call term_info(thresh,0.0        ,0.0        ,'fock')
c
ckw   write(6,*) 'number of integrals > 0.1: ',xn1
ckw   write(6,*) 'number of integrals > 0.01 and <0.1: ',xn2
ckw   write(6,*) 'number of integrals > 0.001 and <0.01: ',xn3
ckw   write(6,*) 'number of integrals > 0.0001 and <0.001: ',xn4
ckw   write(6,*) 'number of integrals < 0.0001: ',xn5
C  release memory
        call retmem(3)              !labels,ischwartz,mapf2s
C   write out stuff left in buffer
       if(nwo.gt.0) then
       nrec=nrec+1
          nfil=(nrec-1)/nrcpf
          nnrc=nrec-nfil*nrcpf
          nndsk=ndisk+nfil
       call Skriv(nndsk,bl(ibufa),lrecw,nnrc,nwo)
       endif
          write(6,*) nrec, 'records written on ', nfil+1, ' files'
Ckw...................................................................
      if(iprint.ge.2) then
       write(6,99)trmz/60.0d0
   99 format(1x,'Time Remzero: ',1f9.2, 'min.')
        write(iout,*) 'Disk space for half-transformed integrals= '
     1  ,dble(nrec*lrec)/1048576.0d0,'MB'
       write(Iout,100) tinteg/60.0d0,ttrans/60.0d0
  100 format('CPU times     Integrals:',f10.2,
     1 ' 1st half transf.: ',f10.2,' min')
       write(Iout,200) elapint/60.0d0,elaptrans/60.0d0
  200 format('Elapsed times Integrals:',f10.2,
     1 ' 1st half transf.: ',f10.2,' min')
       endif
Ckw...................................................................
                      write(91,100) tinteg/60.0d0,ttrans/60.0d0
                      write(91,200) elapint/60.0d0,elaptrans/60.0d0
Ckw...................................................................
              if(iprint.ge.2) then
       write(6,*) 'max ind and instore values: ',indmax,istmax,ncf**2
       write(6,*) '*********************************************'
       write(6,*) 'Statistics from integral transformation'
       write(6,*) '---------------------------------------------'
      write(6,*) 'shell pairs retained: ',iktot,'out of', ncs*(ncs+1)/2
       write(6,*) 'Total number of integrals:',xnumint
       write(6,*) 'Actual number of integrals:',actint
       write(6,*) 'Number of transformed integrals:',transin2
       write(6,*) '*********************************************'
       endif
C****************************** BINSORT starts here    ********
      call secund(tt1)
      call elapsec(elaps1)
       ttot=(tt1-tt0)/60.0d0
       etot=(elaps1-elaps0)/60.0d0
      write(6,*)'Total CPU and elaped time through halftransformation'
      write(iout,110)ttot,etot
  110 format(1x,'CPU time :',f10.2,'min. Elapsed time:',f10.2,'min.')
          call f_lush(6)
C
       call getmem(1,mmsub) ! MMOK
      mmsub=mmsub-ioffset
       mem=kcore-mmsub-(nval**2)-3*munik
       if(iprint.ge.2) Then
      write(6,*) 'Memory available for sort: ',mem
      write(6,*) 'Memory assigned for job: ',kcore
       endif
       call retmem(1)
C change npairs to mstrong for sort
C Munik = mstrong if no symmetry
      write(6,*)nrec,' records of ',lrecw,' integers written,file htr'
       nestim=(nrec*lrecw+2)/3
      write(6,*)'estimated number of nonzero htrs.', nestim
       ninest=nestim*4
C       this is estimate in integer*4 words
       ninesr=ninest/2+1
      lbina=mem/munik
      lbinb=ninesr/nval
C       most nonzero integrals contribute to the diagonal Kijs
       write(6,*)'lbinab',lbina,lbinb
       lbin=min0(lbina,lbinb)
      if(lbin.lt.100) call nerror(1,'lmp2trans','bin too small',
     1lbina,lbinb)
C   nbins= estimate for the number of bins written
C   (must exceed the true value)
  333    mleft=mem-munik*lbin
       if(mleft.lt.0)then
       msub1=-1
       msub2=mleft/munik
       msub=min0(msub1,msub2)
          lbin=lbin+msub
          goto 333
          endif
       lbin2=lbin*2
C       lbin2 is length of 1 bin in integer*4 words
       nbins=Ninest/lbin2+1+munik
C       this is a generous estimate nbins > true value
C       ndisk2 bin file currently 17+4=21
      ndisk2=ndisk+nstripe+1
C
       if(iprint.ge.2) then
      write(6,111)mem,Ninest,lbin2,nbins,lbin*munik,munik
  111 format(' Storage estimates  for Binsort in words: ',/,
     1'available memory:   ', 1I20,/,
     1'length of bin file: ', 1I20,/,
     2'length of one bin : ', 1I20,/,
     3'number of bins:     ',1I20,/,
     4'memory for bins:    ',1I20,/,
     5'number of pairs:    ',1I20)
       endif
C
      bbubyte=dble(lbin2)*ints
      bunbyte=bbubyte*dble(nbins)
      nbinf=bunbyte/gigabyte+1
      call rmblan(scrfile,256,len)
      filname0=scrfile(1:len)//'.bins'
          len1=len+5
      write(6,*) 'binfls: bytes: ',nbins,lbin2,ints,bunbyte,nbinf
      open  (unit=ndisk2,file=filname0(1:len1),form='unformatted',
     1       access='direct',recl=ints*lbin2)
          if(nbinf.gt.1) then
          nbextr=nbinf-1
      Do ibinf=1,nbextr
      write(ch3,'(I3)') ibinf
      ch3(1:1) = '.'
      If(ibinf.lt.10) ch3(2:2) = '0'
      LUnit = ndisk2+ibinf
      filname1 = filname0(1:len1)//ch3
      len = len1+3
      OPEN (UNIT=LUnit,FILE=filname1(1:len),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=ints*lbin2)
      EndDo
          endif
             nrcpb=gigabyte/bbubyte
      write(6,*) 'A total of ', nbinf,' files  opened for bins unit='
     1,ndisk2
      write(6,*) 'Number of records per file: ', nrcpb
Css Aug 01
       if(iprint.ge.2) then
      write(6,112) nbins,lbin2,bunbyte/1048576.0d0
  112 format(1x,'Yoshimine BinSort=',1I20,' bins',1I20,' long
     1 Disk space needed=',1f10.2,'MB')
       end if
c  reserve memory for a counter showing which bin belongs to which pair
       call getint(nbins,ibinctr)
c  reserve memory for counter showing occupancy of 1 bin
       call getint(munik,icounter)
c  reserve memory for bins
       call getint_4(munik*lbin2,iabins)
C       iabins(lbin2,munik)
C
      call Phase1(ncf,nval,ncore,munik,ndisk,
     1            ndisk2,nrec,bl(iplis),lbin2,nbins,
     2            bl(iabins),bl(icounter),bl(ibinctr),bl(ibufa),lrecw,
     3            bl(idimp),bl(lbasp),bl(nstrong),bl(invli),iprint,
     4            bl(ifp),bl(iorrel),bl(iunik),nsym,bl(iunor),
     5            full,bl(ictr),lalb,bl(lbasij),bl(idiij),
     6            bl(invpd),nud,frac,nfiles,nrcpf,
     7            nrcpb,gigabyte)
C
      call secund(tt2)
      call elapsec(elaps2)
      if(iprint.ge.2) then
        write(iout,*) 'CPU and elapsed time for bin sort phase 1: ',
     1  (tt2-tt1)/60.00d0,(elaps2-elaps1)/60.00d0,' min.'
      write(6,113) nbins,lbin2,(bbubyte*dble(nbins))/1048576.0d0
  113 format(1x,'Yoshimine BinSort=',1I20,' bins',1I20,' long
     1 Disk space used=',1f10.2,'MB')
      end if
      call retmem(2)          ! iabins, icounter
C   remove the half-transformed integral file   (XXXXX.htr)
      close(unit=ndisk,status='delete')
      do iclo=1,nstripe
      close(unit=ndisk+iclo,status='delete')
      enddo
C
C***** end binsort phase 1 integrals now on :  XXXXX.bins
C
C**** statistics:
C       avemult=xmult/xnumint
C       write(6,*)'mult:',xmult,'numint:',xnumint
C       xmultav=avemult/(nval**2)
C       write(6,*) 'fraction of integrals transformed:', xmultav
       call secund(tt4)
       call elapsec(elaps4)
C
C       Bin sort Phase 2:
c  read back bins, extract 1 exchange matrix, and transform that
c  to projected basis, and write out on file XXXXX.kij
C
C         construct projection matrix and projected overlap and fock
C
       np1=1
        ifloc=mataddr('floc')
       call matdef ('dnmt','s',ncf,ncf)
       call matdef('ovls','s',ncf,ncf)
C       call matread('ovls',np4,'ovlap')
       call sudat(np1,'s matrix',ni)
       if(ni.gt.0) then
       call matread('ovls',np1,'s matrix')
       else
         call restart(np1,0)
       call matread('ovls',np1,'s matrix')
       end if
       call matdef('fcks','s',ncf,ncf)
       call matread('fcks',np4,'fock_rhf')
       ifoca=mataddr('fcks')
        call matsimtr('fcks','occu','floc')
       call matread('dnmt',np4,'den0_rhf')
       if(redu) then
C       construct SDS (only for redundant formulation)
C
       call matdef('tttt','q',ncf,ncf)
       call matmmult('dnmt','ovls','tttt')
       call matmmult('ovls','tttt','sds')
       call matrem('tttt')
C       call matprint('sds',6)
       else
C       construct projection matrix I-0.5*SD
C
       call matdef('unit','d',ncf,ncf)
       do iset =1,ncf
       call mateset('unit',iset,iset,1.0d0)
       end do
       call matmmult('dnmt','ovls','proj')
        call matscal('proj',-0.5d0)
       call matadd('unit','proj')
       call matrem('unit')
C
C       projection matrix in  'proj'
C
C
       call matsimtr('ovls','proj','ovls')
       call matsimtr('fcks','proj','fcks')
       endif
C       
       call matcopy('ovls','ovl')
       call matcopy('fcks','fck')
       call matrem('fcks')
       call matrem('ovls')
       call matrem('dnmt')
C
C       projected S and F matrices in 'ovl', and 'fck'
C       
      ndisk4=ndisk2+nbinf
      call rmblan(scrfile,256,len)
       filname0=scrfile(1:len)//'.kij'
          len1=len+4
      open(unit=ndisk4,file=filname0(1:len1),form='unformatted')
          nfexr=nfiles-1
          if(nfexr.gt.0) then
      Do ikij=1,nfexr
      write(ch3,'(I3)') ikij
      ch3(1:1) = '.'
      If(ikij.lt.10) ch3(2:2) = '0'
      MUnit = ndisk4+ikij
      filname1 = filname0(1:len1)//ch3
      len = len1+3
      OPEN (UNIT=MUnit,FILE=filname1(1:len),FORM='UNFORMATTED')
      EndDo
          endif
C
       write(6,*) nfiles,ndisk4
          ndisk5=ndisk4+nfiles
C       do phase 2 of bin sort + project the final Kij + write out
C       in local dimension.
C
c  reserve mempory for 1 bin
       call getint_4(lbin*2,i1bin)
C
       ncda=0
       ttda=zero
       ifil=1
       wrdw=0.0d0
          ndskij=ndisk4-1
      call Phase2(ncf,nval,ndisk2,ndskij,lbin2,
     1            nbins,bl(nstrong),bl(ibinctr),lbin,thresh,
     2            ncore,bl(i1bin),bl(ixadr),bl(idimp),bl(lbasp),
     3            bl(iunik),bl(iorrel),bl(iunor),bl(ifp),munik,
     4            nsym,iprint,thres3,bl(ictr),bl(lbasij),
     5            bl(idiij),bl(iekvi),bl(invpd),redu,nud,
     6            ncda,ttda,remf,wrdw,ifil,
     7            nrcpb,giga)
      if(nfiles.gt.ifil) then
      do iclo=ifil+1,nfiles
      close(unit=ndskij+iclo,status='delete')
      enddo
      nfiles=ifil
      endif
      write(6,*)'Number of files used for Kij-operators: ',nfiles
       ttda=ttda/60.0
       write(6,*) ncda,' diagonalizations for ',ttda,' min.'
C
              call retmem(1)       ! i1bin
C
       close(unit=ndisk2,status='delete')
       if(nbinf.gt.1) then
       do iclose=2,nbinf
       close(unit=ndisk2+iclose-1,status='delete')
       enddo
       endif
       call secund(tt7)
       call elapsec(elaps7)
       if(iprint.ge.2) then
       write(6,*)' Projection of Internal Exchange Operators: '
       tt74=tt7-tt4
       el74=elaps7-elaps4
       t74=tt74/60.0d0
       e74=el74/60.0d0
          write(6,*) 'Removal of redundant basis functions:'
       write(6,300) t74,e74
       endif
       t70=tt7-tt0
       e70=elaps7-elaps0
      write(6,*)'Total CPU and elapsed time through sort and projection'
       write(iout,110) t70/60.d0,e70/60.d0
C       release all memory up to the last memory mark
       call retmark
C       at this point you should have:
C
C       projected overlap in quadrtatic form: 'ovl'
C       projected fock matrix in quadratic form  : 'fck'
C       projection matrix in 'proj'
C       fock matrix in localized basis quadratic form: 'floc'
C       All amplitudes stored as a single vector: 'tful'
C
C       The internal exchange matrices in projected local
C       basis are on disk: unit ndisk4,nkunit
C open file for temporaray MO basis used for iterations
C        noru  .orb file should be unit 30
        noru=ndisk5
      call rmblan(scrfile,256,len)
       filname0=scrfile(1:len)//'.orb'
          len1=len+4
      open(unit=noru,file=filname0,form='unformatted')
C       noru+1 xxxxx.orb1 file currently unit 30
          if(nfiles.gt.1) then
          noext=nfiles-1
      Do ioinf=1,noext
      write(ch3,'(I3)') ioinf
      ch3(1:1) = '.'
      If(ioinf.lt.10) ch3(2:2) = '0'
      nUnit = ndisk5+ioinf
      filname1 = filname0(1:len1)//ch3
      len = len1+3
      OPEN (UNIT=nUnit,FILE=filname1(1:len),FORM='UNFORMATTED')
      EndDo
      endif
          write(6,*) 'orbitals', nfiles,ndisk5
C
       if(weak) then
C
       call secund(ttweak1)
       nmo=nval+ncore
C       construct local domains for weak pairs using
C       Pulay Boughton Scheme
c      call getmem((ncf*nval)/4+1,lbasw)
       call getint_2(ncf*nval,lbasw)
c      call getmem(nval/4+1,idmw)
       call getint_2(nval,idmw)
      call pcksetup(nmo,bl(ictr),ncore,bl(lbasw),bl(idmw),
     1   iprint)
C
       epsx=1.0d-07
      call lbnwea(ncf,nval,iprint,bl(iunik),bl(ifp),
     1            nsym,bl(lbasw),bl(idmw),epsx,bl(iekvi),
     2            ncda,ttda)
C
      call Weakpair(ncf,nval,bl(lbasw),bl(idmw),ncore,
     1              iprint,nweak,ndskij,bl(iplis),mstrong,
     2              bl(iunik),ifil,wrdw,giga)
       call secund(ttweak2)
       tweak=ttweak2-ttweak1
       tweak=tweak/60.0d0
       write(6,*) 'time weak pairs: ',tweak
       endif
C
C       read in distance matrix (used for output)
       call matdef('dsijx','q',nval,nval)
       call matdef('dsijs','s',nval,nval)
       call matread('dsijs',np4,'Rijinv')
       call matcopy('dsijs','dsijx')
       ixija=mataddr('dsijx')
       call matrem('dsijs')
C
C       perform MP2-uterations:
C
       call getint(munik,irec)
       call getint(munik,kfil)
      call mp2iter(ncore,ncf,nval,bl(nstrong),bl(ilis),
     1             bl(jlis),bl(invli),bl(idimp),bl(lbasp),bl(itsts),
     2             bl(iplis),ndisk4-1,mstrong,iprint,bl(iekvi),
     3             bl(iunik),bl(ifp),bl(iorrel),nsym,noru-1,
     4             clim,erhf,weak,bl(lbasij),bl(idiij),
     5             nud,bl(invpd),wken,tcou,redu,
     6             nnegl,nweak,bl(iwkli),bl(idmw),bl(lbasw),
     7             bl(ixija),coro,remf,bl(jdimp),bl(jdiij),
     8             bl(kfil),bl(irec),munik,nfiles,giga,
     9             gigabyte)
C
      call secund(tt5)
      call elapsec(elaps5)
       t54=(tt5-tt7)/60.0d0
       e54=(elaps5-elaps7)/60.0d0
       write(6,*)' MP2-iterations:'
       write(6,300) t54,e54
  300 format(1x,'CPU-time  ',1f10.2,' minutes elapsed time',
     1 1f10.2,' minutes')
       do iifl=1,nfiles
       ndisk=ndisk4+iifl-1
       nor=noru+iifl-1
      close(unit=ndisk,status='delete')
      close(unit=nor,status='delete')
       enddo
      end
c==========tra1shel==================================================
      subroutine Tra1Shel(xmo,ncf,ncs,nval,ics,
     *                    icf1,icf2,kcf1,kcf2,
     1                    kcs,ndisk,inx,xint,iprint,
     2                    xmat,xxmat,halftra,nrec,istab,
     3                    istab1,jcmax,ncore,listp,invli,
     4                    thresh,epsi,thres2,nwo,lrec,
     5                    ibufa,mult,nstro,nsym,ifunpair,
     6                    isorb,CT,trmz,transin2,nrow,
     7                    ncol,irow,icol,irow1,icol1,
     8                    icol2,izero,jfun3,smal,nrcpf)
C  Transform two indices of a batch of integrals.
c  Arguments:
C  xmo max absolute mo coefficient
c  ncf  number of contracted basis functions
c  ncs  number of contracted shells
c  nval number of orbitals to be transformed (usually the valence ones)
c  ics,kcs: current contracted shells  loop in transloc
c  inx: array containing contraction info
c  xint: AO integrals ordered like (icf1:icf2,ncf,kcf1:kcf2,ncf)
c  where icf1 is the first contr function in the shell ics, and icf2 is
c  the last; similarly for kcs . NOTE THAT THEY ARE SCALED by 1/thresh
c  iprint: print level
C  xmat - xmat(ncf,ncf) for storage of one matrix of AO integrals, this
C       matrix is extracted from the bloc of integrals in Extract1IK
C  xxmat - same as xmat but one dimensional
C  halftra - array for haldtransformed integrals halftra(nval,nval)
C       for a given my,lam halftra(i,j) = (my,i,lam,j)
c  nrec:  the number of records written  (input/output)
C  istab(ncf,nval)  C(ny,istab(ny,imo)) will give the C(ny,imo) in
C  istab1(jcmax,nval,ncf) 3-dim. sort array similar to istab. see mktab
C  for details Note istab1(1,ii,my) gives the number of j's to follow
C  order of decreasing magnitude
C  jcmax  maximum number of j's that form strong pairs with any i +1
C  listp - pair list triangular
C       negligible pairs: -1
C         weak pairs:     :  0
C       strong pairs    :  a positive number
C  invli - invli(nval,nval)   inverse pair list invli(i,j) and
C          invli(j,i) will give the pairnumber for pair i,j
C              = (j-1)*nval+ii
C  thresh,epsi,thres2 threshold
C       thresh intgreal threshold
C    epsi  threshold for negligible contributions to the transformed
C             integrals (used in exit and cycle tests below. (these
C             are essential for this algorithm which is formally ON6)
C       thres2 threshold for removing negligible contributions,used
C             in remzero
C   nwo  word counter in buffer of haltransformed integrals
C   lrec  size of buffer in 4-byte words
C   ibufa start address for buffer
C    mult,  counters for testpurposes will be removed in the final
C               version
C   nstro   logical*1 array nstro(i,j) is false  if ij is a strong pair
c  OUTPUT:
c  halftr: space for half-transformed integrals to be written on disk
C
C       called from transloc
C
C       External calls:
C   Extract1IK  extract one matrix xmat(ncf,ncf) of integrals from
C               the block
C   REMZERO     from  the matrxi fo halftransformed integrals
C               halftra(nval,nval) all elements smaller than thres2
C               are removed.  The surviving integrals are written out
C               on xxx.htra with 4 indices i,j and my,ny
C

      use memory

      implicit real*8 (a-h,o-z)
c     common /big/bl(30000)
      dimension  xint(*),xmat(ncf,ncf),halftra(jcmax-1,*),xxmat(*)
      dimension inx(12,*),CT(nval,*)
       logical*1 nstro(nval,nval)
       dimension ifunpair(7,*),isorb(nval,*)
       dimension istab(nval,ncf)
       integer*2 istab1(jcmax,nval,ncf)
       dimension listp(*),invli(nval,nval)
       dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
       dimension icol2(nval),jfun3(jcmax,nval)
       integer*2 izero(nval),null
       logical smal
       parameter (zero=0.0d0,one=1.0d0)
       parameter (null=0)
C
C       if(nrow.eq.0.or.ncol.eq.0) RETURN
       jcm=jcmax-1
      eps2=epsi/thresh
      eps3=thres2/thresh
       nvasq=nval**2
C     kcf1=inx(11,kcs)+1
C     kcf2=inx(10,kcs)
C     icf1=inx(11,ics)+1
C     icf2=inx(10,ics)
      mult=0
C     do kcf=kcf1,kcf2
C      do icf=icf1,icf2
C     if(icf.lt.kcf) cycle
        do icf=icf1,icf2
      do kcf=kcf1,kcf2
      if(kcf.gt.icf) exit
C
       if(smal) then
       call Extract1IK(ncf,icf,kcf,icf1,icf2,
     2      kcf1,kcf2,xint,xmat,nrow,
     3      ncol,irow,icol,nrow1,ncol1,irow1,icol1)
       else
       call Extract1KI(ncf,icf,kcf,icf1,icf2,
     2      kcf1,kcf2,xint,xmat,nrow,
     3      ncol,irow,icol,nrow1,ncol1,irow1,icol1)
       endif
C
       if(nrow1.eq.0.or.ncol1.eq.0) cycle
C
       rintm=zero
       do icx=1,ncol1
       ixmxx=idamax(nrow1,xmat(1,icx),1)
       rintm=max(rintm,abs(xmat(ixmxx,icx)))
       enddo
C
C       ixmxx=idamax(nrow1*ncol1,xxmat,1)
C       rintm=abs(xxmat(ixmxx))
         xmoxint=xmo*rintm
          if(xmo*xmoxint.lt.eps2)cycle
C
C  Now transform integrals: generate Kmyla(i,j)
C     call zeroit(halftra,nvasq)
C
       do iz=1,nval
       izero(iz)=null
          enddo
       ncol2=0
C
C         do ny=1,ncf
       do irr=1,nrow1
       ny=irow1(irr)
            if(abs(xmoxint*CT(istab(1,ny),ny)).lt.eps2) cycle
C       do isig=1,ncf
       do icc=1,ncol1
       isig=icol1(icc)
           zint=xmat(irr,icc)
C        integrals at this stage are scaled
       if(abs(zint).lt.one) cycle
       if(abs(zint*CT(istab(1,ny),ny)*CT(istab(1,isig),isig)).lt.eps2)
     1 cycle
         do imo=1,nval
         ii=istab(imo,ny)
          cnyix=zint*CT(ii,ny)
C       mult=mult+1
        if(abs(cnyix*CT(istab(1,isig),isig)).lt.eps2) exit
        if(abs(cnyix*CT(istab1(2,ii,isig),isig)).lt.eps2) cycle
C   *****for  most (at least stretched)
C   ***** molecules there is only about 6-8 js left and
C   ***** it does not make sense to go through gjmax2 to determine
C   ***** jmax we use istab1(1,ii,isig) instead of jmax
C       call gjmax2(CT(1,isig),cnyix,istab1,isig,jmax, cmax,
C    1  ncf,nval,eps2,ii)
C       do jmo=1,jmax
C
       jmax=jfun3(1,ii)-1
       if(izero(ii).eq.null) then
       ncol2=ncol2+1
       icol2(ncol2)=ii
          izero(ii)=ncol2
       call zeroit(halftra(1,ncol2),jmax)
       endif
       ncol3=izero(ii)
C      do jmo=1,istab1(1,ii,isig)
C       jj=istab1(jmo+1,ii,isig)
C       mult=mult+1
       do jmo=1,jmax
       jj=jfun3(jmo+1,ii)
C**********  transformation carried out here *********************
C     halftra(jj,ii)=halftra(jj,ii)+cnyix*CT(jj,isig)
      halftra(jmo,ncol3)=halftra(jmo,ncol3)+cnyix*CT(jj,isig)
C*****************************************************************
C
       end do !do jmo
        end do !do imo
       end do !do isig
      end do !do ny
C   Convert to integers and write out on disk
       call secund(trz1)
       if(ncol2.eq.0)cycle
      call remZero(halftra,bl(ibufa),eps3,nval,icf,
     1             kcf,jcm,iprint,lrec,nrec,
     2             nwo,ndisk,ncore,icol2,jfun3,
     3             ncol2,transin2,nrcpf)
       call secund(trz2)
       trmz=trmz+trz2-trz1
      end do
      end do
      end
c=============Extract1IK===============================================
      subroutine Extract1IK(   ncf,   icf,  kcf,   icf1,  icf2,
     1                        kcf1,  kcf2,  xint,  xmat,  nrow,
     2                        ncol,  irow,  icol,  nrow1, ncol1,
     3                        irow1, icol1)
c  This routine extracts a single matrix from the integral array
c  The latter is indexed as
C  xint(l=1:ncf,j=1:ncf,k=kcf1:kcf2,ncf,i=icf1:icf2);
c  elements (sigma,nu,kcf,icf) are put in xmat(nu,sigma)
c  It does not remove the overall scaling. The permutational
c  factors are not put in by the integrals program
c  Arguments:
c  INTENT: IN
c  ncf = number of contracted basis functions
c  icf, kcf = fixed AO indices of the first and third AOs
c  icf1,icf2 = The present shell goes from icf1 to icf2
c  kcf1,kcf2 = Same for kcf
c  xint = integral array
c  nrow = the number of non-zero rows
c  ncol = the number of non-zero columns
c  irow(1:ncf) = the first nrow elements of this array give the
c                original positions of the non-zero rows after
c                compacting
c  icol(1:ncf) = the same as above for the columns
c  INTENT:OUT
c  xmat = the exchange matrix extracted is put in xmat
c  nrow1, ncol1: the actual number of non-zero rows & columns
c  irow1(ncf),icol1(ncf): irow1(k) gives the original row number
c  of row k in the compacted matrix xmat
      implicit real*8 (a-h,o-z)
      dimension xint(ncf,ncf,kcf1:kcf2,icf1:icf2),xmat(ncf,ncf)
      dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
      do ll=1,ncol
        l=icol(ll)
        do jj=1,nrow
          j=irow(jj)
          xmat(jj,ll)=xint(l,j,kcf,icf)
        end do
      end do
c  matrix in xmat. Note that xmat is transmitted here as a parameter
c  but it is also a defined matrix
c  now compact xmat again: eliminate zero rows and columns
c  It is easier to do this for the columns
      ncol1=0
      do ll=1,ncol
        l=icol(ll)
        do jj=1,nrow
c  at this point, the integrals are scaled, so quantities below 1.0
c  are negligible
          if(abs(xmat(jj,ll)).gt.one) then
c column ll has at least one non-zero element
             ncol1=ncol1+1
             icol1(ncol1)=l
             go to 100
          endif
        end do
        go to 150
  100   continue
        if(ll.gt.ncol1) then
          do jj=1,nrow
            xmat(jj,ncol1)=xmat(jj,ll)
          end do
        end if
  150   continue
      end do
      nrow1=0
      do jj=1,nrow
        j=irow(jj)
        do ll=1,ncol1
          if(abs(xmat(jj,ll)).gt.one) then
            nrow1=nrow1+1
            irow1(nrow1)=j
            go to 200
          endif
        end do
        go to 250
  200   continue
        if(jj.gt.nrow1) then
          do ll=1,ncol1
            xmat(nrow1,ll)=xmat(jj,ll)
          end do
        end if
  250   continue
      end do
      end
c=============Extract1KI===============================================
      subroutine Extract1KI(   ncf,   icf,  kcf,   icf1,  icf2,
     1                        kcf1,  kcf2,  xint,  xmat,  nrow,
     2                        ncol,  irow,  icol,  nrow1, ncol1,
     3                        irow1, icol1)
c  This routine extracts a single matrix from the integral array
c  The latter is indexed as
C  xint(l=1:ncf,j=1:ncf,k=kcf1:kcf2,ncf,i=icf1:icf2);
c  elements (sigma,nu,kcf,icf) are put in xmat(nu,sigma)
c  It does not remove the overall scaling. The permutational
c  factors are not put in by the integrals program
c  Arguments:
c  INTENT: IN
c  ncf = number of contracted basis functions
c  icf, kcf = fixed AO indices of the first and third AOs
c  icf1,icf2 = The present shell goes from icf1 to icf2
c  kcf1,kcf2 = Same for kcf
c  xint = integral array
c  nrow = the number of non-zero rows
c  ncol = the number of non-zero columns
c  irow(1:ncf) = the first nrow elements of this array give the
c                original positions of the non-zero rows after
c                compacting
c  icol(1:ncf) = the same as above for the columns
c  INTENT:OUT
c  xmat = the exchange matrix extracted is put in xmat
c  nrow1, ncol1: the actual number of non-zero rows & columns
c  irow1(ncf),icol1(ncf): irow1(k) gives the original row number
c  of row k in the compacted matrix xmat
      implicit real*8 (a-h,o-z)
      dimension xint(ncol,nrow,kcf1:kcf2,icf1:icf2),xmat(ncf,ncf)
      dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
C
       do ll=1,ncol
       do jj=1,nrow
       xmat(jj,ll)=xint(ll,jj,kcf,icf)
       enddo
       enddo
c  matrix in xmat. Note that xmat is transmitted here as a parameter
c  but it is also a defined matrix
c  now compact xmat again: eliminate zero rows and columns
c  It is easier to do this for the columns
      ncol1=0
      do ll=1,ncol
        l=icol(ll)
        do jj=1,nrow
c  at this point, the integrals are scaled, so quantities below 1.0
c  are negligible
          if(abs(xmat(jj,ll)).gt.one) then
c column ll has at least one non-zero element
             ncol1=ncol1+1
             icol1(ncol1)=l
             go to 100
          endif
        end do
        go to 150
  100   continue
        if(ll.gt.ncol1) then
          do jj=1,nrow
            xmat(jj,ncol1)=xmat(jj,ll)
          end do
        end if
  150   continue
      end do
      nrow1=0
      do jj=1,nrow
        j=irow(jj)
        do ll=1,ncol1
          if(abs(xmat(jj,ll)).gt.one) then
            nrow1=nrow1+1
            irow1(nrow1)=j
            go to 200
          endif
        end do
        go to 250
  200   continue
        if(jj.gt.nrow1) then
          do ll=1,ncol1
            xmat(nrow1,ll)=xmat(jj,ll)
          end do
        end if
  250   continue
      end do
      end
C==========writebin1============================================
      subroutine WriteBin1(ndisk2,nbin,lbin2,ibin)
      integer*4 ibin(lbin2)
      write(ndisk2,rec=nbin) ibin
      end
C===========ReadBin_loc============================================
      subroutine ReadBin_loc(ncf,nval,ij,ndisk2,lbin2,
     1                       nbins,ibinctr,lbin,thresh,i1bin,
     2                       xmat,ymat,ndiag,nrcpb)
      implicit real*8 (a-h,o-z)
      integer*2 iindex(2)
      integer*4 i1bin(2,lbin)
      dimension xmat(ncf,ncf),ibinctr(*)
       dimension ymat(*)
       logical ndiag
      equivalence (longint,iindex(1))
c  in this first implementation, all records are scanned for
c  contributions to pair (ij); this should be replaced by a
c  proper sort
      call zeroit(ymat,ncf**2)
C       if(ndiag) then
      do k=1,nbins
        if(ibinctr(k).ne.ij) cycle
        nnbfil=(k-1)/nrcpb
        nndisk=ndisk2+nnbfil
        newrec=k-nnbfil*nrcpb
C         read(ndisk2,rec=k) i1bin
          read(nndisk,rec=newrec) i1bin
          do j=1,lbin
            longint=i1bin(2,j)
            if(longint.le.0) exit
            mu=iindex(1)
            lam=iindex(2)
C           if(mu.eq.0.or.lam.eq.0) cycle
            if(mu.le.0.or.mu.gt.ncf.or.lam.le.0.or.lam.gt.ncf) then
              call nerror(1,'readbin',
     1        'wrong index', mu,ncf)
            else
               xmat(mu,lam)=i1bin(1,j)*thresh
C       write(6,*)'X(',mu,',',lam,')=',xmat(mu,lam)
            end if
          end do  ! j elements in the bin
        end do    ! k; bins
C       else
c     do k=1,nbins
c       if(ibinctr(k).ne.ij) cycle
c         read(ndisk2,rec=k) i1bin
c         do j=1,lbin
c           longint=i1bin(2,j)
c           if(longint.le.0) exit
c           mu=iindex(1)
c           lam=iindex(2)
c       if(lam.gt.mu) cycle
C           if(mu.eq.0.or.lam.eq.0) cycle
c          if(mu.le.0.or.mu.gt.ncf.or.lam.le.0.or.lam.gt.ncf) then
c             call nerror(1,'readbin',
c    1        'wrong index', mu,ncf)
c           else
c              xmat(mu,lam)=i1bin(1,j)*thresh
c              xmat(lam,mu)=i1bin(1,j)*thresh
C       write(6,*)'X(',mu,',',lam,')=',xmat(mu,lam)
c           end if
c         end do  ! j elements in the bin
c       end do    ! k; bins
c       endif
      end
C============sortcoeff===============================================
      subroutine SortCoeff(ncf,nval,C,isort)
      implicit real*8 (a-h,o-z)
C       this subroutine creates the sort-matrix isort(nval,ncf)
C       used throughout the integrals transformation
C       for a given AO iao isort(1,iao) = imo
C       and C(iao,imo) will be the largest coefficient in the row
C       iao.   C(iao,isort(2,iao)) will be the second  largest etc.
C       This allows us to process the C's in order of
C       decreasing magnitude.
C       note that the input matrix C is the transpose of the
C       origninal C.  This matrix is also sorted in the process
C       but the sorted C + will not be used later.
C       ncf - number of contraced basis functions
C       nval - number of valence orbitals
      dimension C(nval,ncf),isort(nval,ncf)
C
      do nu=1,ncf
        do j=1,nval
          isort(j,nu)=j
        end do
 100    itrp=0
        do j=1,nval-1
          if(abs(C(j,nu)).lt.abs(C(j+1,nu))) then
            x=C(j+1,nu)
            C(j+1,nu)=C(j,nu)
            C(j,nu)=x
            ir=isort(j+1,nu)
            isort(j+1,nu)=isort(j,nu)
            isort(j,nu)=ir
            itrp=1
          end if
        end do !j
        if(itrp.eq.1) go to 100
      end do   !nu
      end
C==========testso===================================================
      subroutine testso(isort,ncf,nval)
C       not called anymore, simply prints the sort table
       integer isort(nval,ncf)
       write(6,*) 'sortindices:'
       do icol=1,ncf
       write(6,*) (isort(iro,icol),iro=1,nval)
       enddo
       end
c==========evalso====================================================
       subroutine evalso(C,isort,ipar,ncf,nval)
C       evaluates the effectivenss of the sort decribed in 'SortCoeff'
C       Calculates the number of pairs to which a given AO
C             will contribute
C       not used any more
       implicit real*8(a-h,o-z)
       dimension C(nval,ncf),isort(nval,ncf),ipar(*)
       iptot=0
       do iba=1,ncf
       ij=0
       ip=0
       do  iva=1,nval
       do jva=1,iva
       ij=ij+1
       fact=c(iva,iba)*c(jva,iba)
       if(abs(fact).gt.1.0d-05) then
       ip=ip+1
       ipar(ip)=ij
       endif
       end do
       end do
       write(6,*)'AO no: ',iba,' contributes to: ',ip,' pairs'
       write(6,*) (ipar(iii),iii=1,ip)
       iptot=iptot+ip
      enddo
      rpav=float(iptot)/ncf
      nptot=nval*(nval+1)/2
      write(6,*) 'Each AO contributes to an average of ',rpav,' pairs'
      write(6,*) 'Total number of pairs = ', nptot
             end
C===========plist======================================================
       subroutine plist(x,y,z,r2,listp,
     1                     ilist,jlist,invl,nstro,nval,
     2                     jcmax,nstrong,thre1,thre2,nud,
     3                     wken,weak,nnegl,nweak,iprint)
C.....construct pair-lists, classifies pairs as strong, weak,
C.....and negligible.
C.....X,Y,Z,r2 are centers and diameter of localized orbitals
C..... taken from LOCA
C
C       at the end the pair list is changed in the following
C       manner:
C.....negligible pairs flagged with : -1
C.....weak pairs flagged with       : 0
C.....for a strong pair Listp(lpar) = a positive number
C.....listp(lpar) = mstrong
C     lpar is strong pair # mstrong
C
C    there are 5 lists
C    ilist (ipar) = imo
C    For valence pair ipar (i,j) i=imo
C    jlist(ipar) = jmo
C    Invli(nval,nval)  Inverse pair list. NOTE this is
C    quadratic to allow for j<i
C    Invli(i,j) = ipar
C       nstro(nval,nval)   character*1 array
C       nstro(ii,jj) = .true.    if the pair ii,jj is not
C       a strong pair
C    nval number of valence orbitals : occupied- core
C
C    Svein Saebo  Fayetteville AR, May 1999
C

      use memory

      implicit real*8(a-h,o-z)
      dimension X(*),Y(*),Z(*),R2(*)
       dimension listp(*),ilist(*),jlist(*),invl(nval,nval)
       logical*1 nstro(nval,nval)
c     common /big/bl(30000)
       logical weak
       parameter(zero=0.0d0,one=1.0d0)
C
C        thre1=2.0d-6
C       thre2=3.0d-4
C    pairs with estimated pair energy < thre1 are negligible
C    pairs with estimated pair energy < thre2 (but > thre1 are weak
C       
C       save distance matrix for multipole expansion of week pairs
C       we are actually saving the inverse distance matrix
C       we'll also need the Rij vector components
C
       wken=zero
       trwke=thre2
       if(weak) trwke=thre1
       call matdef('r2x','d',nval,nval)
       call matdef('r2y','d',nval,nval)
       call matdef('r2z','d',nval,nval)
       ir2x=mataddr('r2x')-1
       ir2y=mataddr('r2y')-1
       ir2z=mataddr('r2z')-1
       call matdef('Rijx','s',nval,nval)
       call matdef('Rijy','s',nval,nval)
       call matdef('Rijz','s',nval,nval)
       call matdef('dmo','s',nval,nval)
       idmoa=mataddr('dmo')-1
       irijx=mataddr('Rijx')-1
       irijy=mataddr('Rijy')-1
       irijz=mataddr('Rijz')-1
C
       npar=0
       do imo=1,nval
       bl(ir2x+imo)=x(imo)
       bl(ir2y+imo)=y(imo)
       bl(ir2z+imo)=z(imo)
       do jmo=1,imo
       npar=npar+1
       bl(idmoa+npar)=zero
       bl(irijx+npar)=zero
       bl(irijy+npar)=zero
       bl(irijz+npar)=zero
       invl(imo,jmo)=npar
       invl(jmo,imo)=npar
C....calculate distance between orbital centers
      if(imo.eq.jmo) then
      fact=1000.d0
      goto 100
      endif
       rx2=x(jmo)-x(imo)
       ry2=y(jmo)-y(imo)
       rz2=z(jmo)-z(imo)
       x2=rx2**2
       y2=ry2**2
       z2=rz2**2
       d2=x2+y2+z2
       R=sqrt(d2)
       bl(idmoa+npar)=one/R
       bl(irijx+npar)=rx2/R
       bl(irijy+npar)=ry2/R
       bl(irijz+npar)=rz2/R
       d6=d2**3
       rrI=sqrt(r2(imo)-x(imo)**2-y(imo)**2-z(imo)**2)
       rrJ=sqrt(r2(jmo)-x(jmo)**2-y(jmo)**2-z(jmo)**2)
       alpI=rrI**3
       alpJ=rrJ**3
       Fact=0.625d0*alpI*alpJ/d6
  100 listp(npar)=2
       if(Fact.lt.thre2)then
       listp(npar)=1
       endif
       if(Fact.lt.thre1)then
       listp(npar)=0
       endif
       if(fact.lt.trwke) wken=wken-fact
       enddo
       enddo
C
C       write out inverse distance matrix
       np4=4
       call matwrite('r2x',np4,0,'r2x')
       call matwrite('r2y',np4,0,'r2y')
       call matwrite('r2z',np4,0,'r2z')
       call matwrite('Rijx',np4,0,'Rijx')
       call matwrite('Rijy',np4,0,'Rijy')
       call matwrite('Rijz',np4,0,'Rijz')
       call matwrite('dmo',np4,0,'Rijinv')
       call matrem('dmo')
       call matrem('Rijz')
       call matrem('Rijy')
       call matrem('Rijx')
       call matrem('r2z')
       call matrem('r2y')
       call matrem('r2x')
C
       npairs=npar
C......
       nweak=0
       nstrong=0
       nnegl=0
C set up logical array of pairs to be neglected during transformation
       jcmax=0
       do iii=1,nval
       jcount=0
       do jjj=1,nval
       nstro(iii,jjj)=.true.
       nstro(jjj,iii)=.true.
       npar=invl(iii,jjj)
       if(listp(npar).eq.2) then
       jcount=jcount+1
       nstro(iii,jjj)=.false.
       nstro(jjj,iii)=.false.
       endif
       end do
       jcmax=max0(jcmax,jcount)
       end do
       if(iprint.ge.2)
     1 write(6,*) 'max no of j belonging to a strong pair: ',jcmax
       do ipa=1,npairs
       goto(1,2,3) listp(ipa)+1
    1  nnegl=nnegl+1
          cycle
    2 nweak=nweak+1
          cycle
    3 nstrong=nstrong+1
         enddo
C   change pair list to give pairnumbers for strong pairs
       nstr=0
       do ipa=1,npairs
       if(listp(ipa).eq.0)then
              listp(ipa)=-1
       cycle
       end if
       if(listp(ipa).eq.1)then
              listp(ipa)=0
       cycle
       endif
       nstr=nstr+1
       listp(ipa)=nstr
       end do
       lpar=0
       do ii=1,nval
       do jj=1,ii
       if(nstro(ii,jj)) cycle
       lpar=lpar+1
       llpa=invl(ii,jj)
       llpar=listp(llpa)
       if(lpar.ne.llpar) then
       write(6,*) ii,jj,lpar,llpa,llpar
       call nerror(1,'plist','error',0,0)
       endif
       ilist(lpar)=ii
       jlist(lpar)=jj
       end do
       end do
       if(weak)then
       ij=0
       do ii=1,nval
       do jj=1,ii
       ij=ij+1
       if(listp(ij).ne.0) cycle
       lpar=lpar+1
       ilist(lpar)=ii
       jlist(lpar)=jj
       enddo
       enddo
       endif
       write(6,*)'Number of neglected pairs: ',nnegl
       write(6,*)'Number of weak pairs     : ',nweak
       write(6,*)'Number of strong pairs   : ',nstrong
       write(6,*)'Total number of pairs    : ',npairs
       nud=nstrong-nval
       write(6,*) 'Number of non-diagonal pairs included: ', nud
         end
C========printlist===============================================
        subroutine printlist(listp,ilist,jlist,invli,nval,iprint,mstron)
C    Prints the pair lists.  For testing purposes.
C    there are 4 lists
C    listp listp (ipar) should be -1,0 or a positive number
C    -1- neglected pair
C    0 - weak pair
C    pos  - strong pair number
C    ilist (ipar) = imo
C    For valence pair ipar (i,j) i=imo
C    jlist(ipar) = jmo
C    Invli(nval,nval)  Inverse pair list. NOTE this is
C    quadratic to allow for j<i
C    Invli(i,j) = ipar
C    nval number of valence orbitals : occupied- core
C
C    Svein Saebo  Fayetteville AR, May 1999
C
       implicit real*8(a-h,o-z)
       dimension listp(*),ilist(*),jlist(*),invli(nval,nval)
       if(iprint.ge.2) then
         write(6,*) 'Pair List:'
         ijp=0
         do imo=1,nval
         write(6,*) Imo,':',(listp(ijp+jmo),jmo=1,imo)
         ijp=ijp+imo
       end do
       endif
       npairs=nval*(nval+1)/2
       if(iprint.ge.3) then
       write(6,*) 'ILIST:'
       write(6,*) (ilist(ii),ii=1,mstron)
       write(6,*) 'JLIST:'
       write(6,*) (jlist(ii),ii=1,mstron)
       endif
       if(iprint.ge.3) then
       write(6,*) 'Inverse Pair List:'
       do imo=1,nval
       write(6,*) (invli(imo,jmo),jmo=1,nval)
       end do
       end if
       end
C=========Skriv==============================================
       subroutine Skriv(nunit,a,ldim,nrec,nwo)
C        zeros out unused spaces in a buffer and writes it out
C
C       Svein Saebo Fayetteville Ar May 1999
C
       implicit integer(a-z)
       dimension a(ldim)
C
       left=ldim-nwo
       if(left.le.0) then
       call nerror(1,'skriv','wrong wordcount',ldim,nwo)
       else
       a(nwo+1)=-1
       nwo=0
C       if(left.ge.2) then
C       do ip=1,left-1
C       a(nwo+ip)=0
C        end do
C       end if
       end if
       write(unit=nunit,rec=nrec) a
       end
C=======gjmax==================================================
       subroutine gjmax(C,cnyix,istab,isig,jmax,ncf,nval,eps2)
C
       implicit real*8(a-h,o-z)
       dimension C(ncf,*)
       integer istab(nval,*)
C       estimates the maximum jmo for the innermost loop in
C       the integral transformation.
C       ARGUMENTS:
C       Input:
C       C(ncf,nval)       MO coefficients
C       cnyix              C(ny,i)*x(ny,sigma)
C       istab(nval,ncf) table with sortindices
C       isig              AO index sigma
C       ncf              no. of contracted basis functions
C       nval              no. of valence MOs
C       eps2              epsi/thresh (epsi for unscaled integrals)
C       Output:
C       jmax              max value of jmo for inner loop
C       called from tra1shel inside
C       do ny
C       do imo
C       do isig
C       Svein Saebo, Fayetteville, AR June 1999
       j0=1
       j1old=nval
       j1=nval/2
  100 continue
        if(abs(C(isig,istab(j1,isig))*cnyix).lt.eps2) then
C                     if yes:
              j1old=j1
              j1=j1-(j1-j0)/2
              else
C                     if no:
              j0=j1
              j1=j1old
              endif
              if((j1-j0).gt.1)goto 100
C              jmax found within 1
         jmax=j1old
       end
C=======mktab=========================================================
       subroutine mktab(istab,istab1,ncf,nval,jcmax,nstro,iprint)
C       creates table istab1(j,i,my)
C       C(my,istab1(j,i,my)) is  like a C(my,j) but the elements come in
C       in decreasing order for the
C       j's that form strong pairs with i
C       see subroutine sortcoeff that creates the related table
C       istab(i,my)
C
C       Input:
C       istab(nval,ncf)  table created in sortcoeff
C       ncf number of contracted basis functions
C       nval number of basis functions
C       jcmax max number of j's that form strong pairs with any i
C       nstro lohical arrays nstro(i,j) is false if ij is a strong pair
C       iprint       printlevel
C       Output
C       istab1(jmax,nval,ncf)
C
C       Svein Saebo Fayetteville AR Summer 1999
C
C   mktab is called from loctrans and the table is used in htra1shel
       implicit real*8(a-h,o-z)
       logical*1 nstro(nval,nval)
       integer*2 istab1(jcmax,nval,ncf)
       integer istab(nval,ncf)
       if(iprint.ge.4) then
       write(6,*) 'istab'
       do jjjj=1,ncf
       write(6,*) (istab(iw,jjjj),iw=1,nval)
       end do
       endif
       icnt=0
       do imo=1,nval
       do isig=1,ncf
       jcount=1
       do jmo=1,nval
       jj=istab(jmo,isig)
       if(nstro(jj,imo))cycle
       jcount=jcount+1
       istab1(jcount,imo,isig)=jj
       end do ! jmo
       if(jcount.gt.jcmax) call nerror(1,'mktab','toobig',jcount,
     1 jcmax)
       istab1(1,imo,isig)=jcount-1
       icnt=icnt+jcount-1
       end do  ! isig
       end do  !  imo
       avej=float(icnt)/float(ncf*nval)
       if(iprint.ge.2)
     1 write(6,*) 'Average jmax-value in inner loop: ', avej
C   print table if desired (warning 3-d but jcmax is normally small~10
       if(iprint.ge.6) then
       do imo=1,nval
       do isig=1,ncf
       write(6,*)  'istab1(',imo,',',isig,')','=',
     1 (istab1(jj,isig,imo),jj=1,jcmax)
       end do
       end do
       endif
       end
C========gjmax2================================================
       subroutine gjmax2(C,cnyix,istab1,isig,jmax,jcmax,ncf,nval
     1  ,eps2,imo)
C
       implicit real*8(a-h,o-z)
       dimension C(*)
       integer*2 istab1(jcmax,nval,ncf)
C       estimates the maximum jmo for the innermost loop in
C       the integral transformation.
C       ARGUMENTS:
C       Input:
C       C(ncf,nval)       MO coefficients
C       cnyix              C(ny,i)*x(ny,sigma)
C       isig              AO index sigma
C       ncf              no. of contracted basis functions
C       nval              no. of valence MOs
C       eps2              epsi/thresh (epsi for unscaled integrals)
C       Output:
C       jmax              max value of jmo for inner loop
C       called from tra1shel inside
C       do ny
C       do imo
C       do isig
C             Svein Saebo, Fayetteville, AR June 1999
C
       j0=1
       j1old=istab1(1,imo,isig)
       j1=j1old/2
  100 continue
        if(abs(C(istab1(j1+1,imo,isig))*cnyix).lt.eps2) then
C                     if yes:
              j1old=j1
              j1=j1-(j1-j0)/2
              else
C                     if no:
              j0=j1
              j1=j1old
              endif
              if((j1-j0).gt.1)goto 100
C              jmax found within 1
         jmax=j1old
       end
C======evalc==========================================================
       Subroutine EvalC(C,nval,ncf)
       implicit real*8(a-h,o-z)
       dimension C(ncf,nval)
       n1=0
       n2=0
       n3=0
       n4=0
       n5=0
       do imo=1,nval
       do iao=1,ncf
       elem=abs(C(iao,imo))
       if (elem.gt.1.0d-1) then    !  > 0.1
       n1=n1+1
       cycle
       endif
       if (elem.gt.1.0d-2) then    !   >0.01
       n2=n2+1
       cycle
       endif
       if (elem.gt.1.0d-3) then    !   >0.001
       n3=n3+1
       cycle
       endif
       if (elem.gt.1.0d-4) then    !   >0.0001
       n4=n4+1
       cycle
       endif
       n5=n5+1
       end do
       end do
       tot=ncf*nval
       p1=(float(n1)/tot)*100
       p2=(float(n2)/tot)*100
       p3=(float(n3)/tot)*100
       p4=(float(n4)/tot)*100
       p5=(float(n5)/tot)*100
       write(6,*)p1,' % of coeffcients > 0.1'
       write(6,*)p2,' % of coeffcients > 0.01 and <0.1'
       write(6,*)p3,' % of coeffcients > 0.001 and <0.01'
       write(6,*)p4,' % of coeffcients > 0.0001 and <0.001'
       write(6,*)p5,' % of coeffcients < 0.0001'
       end
C=======evali=====================================================
       subroutine EvalI(X,nsize,thresh,xn1,xn2,xn3,xn4,xn5)
       implicit real*8(a-h,o-z)
       dimension X(nsize)
       n1=0
       n2=0
       n3=0
       n4=0
       n5=0
       do ic=1,nsize
       elem=abs(x(ic)*thresh)
       if(elem.ge.1.0d-1) then
       n1=n1+1
       cycle
       endif
       if(elem.ge.1.0d-2) then
       n2=n2+1
       cycle
       endif
       if(elem.ge.1.0d-3) then
       n3=n3+1
       cycle
       endif
       if(elem.ge.1.0d-4) then
       n4=n4+1
       cycle
       endif
       n5=n5+1
       end do
       xn1=xn1+n1
       xn2=xn2+n2
       xn3=xn3+n3
       xn4=xn4+n4
       xn5=xn5+n5
      end
C=======check_Int============================================
       subroutine Check_Int(ncs,schwarz,ics,jcs,val)
       implicit real*8(a-h,o-z)
       dimension schwarz(ncs,ncs)
       val=schwarz(ics,jcs)
       end
C========dmaxij===============================================
      subroutine DmaxIJ(D,DS,ncf,ncs,nval,nstron,inx,C,densmax,iprint)
C       Calculates Max (over ij strong) of C(my,i) * C(lam,j)
C       as described in paper by Schutz, Hetzer, and Werner
C       Used for prescreening of integrals.
C       Output:
C       D(ncf,ncf) (this matrix is removed after return)
C         DS(ncs,ncs) same as D but over shells rather than
c         densmax - maximum element of "screening" density
c
C       individual functions
C       Input:
C       ncf - number of basis functions
C       ncs - number of contracted shells
C       nval - number of valence orbitals
C       nstron - logical array nstron(i,j) is true if i,j not a
C       strong pair
C       C(ncf,nval) matrix of MO coefficients (localized orbitals)
C
C       Svein Saebo, Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
       dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs)
       dimension inx(*)
       logical*1 nstron(nval,nval)
       ids=0
C
       do my=1,ncf
       do lam=1,ncf
       dmxx=0.0d0
       do imo=1,nval
       do jmo=1,nval
          if(nstron(imo,jmo)) cycle
       cprod=abs(C(my,imo)*C(lam,jmo))
       dmxx=max(dmxx,cprod)
       end do
       end do
       D(my,lam)=dmxx
       end do
       end do
c
C       print statstics
       if(iprint.ge.2)call evaldmax_loc(ncf,D,ids)
c
C  construct a matrix over contracted shells rather than basis function
c
          ictr=igetival('ictr')
c
       ids=1
          densmax=0.d0
          ipoint=0
          do ics=1,ncs
             call get_shell_size(bl(ictr),ics,isize)
             kpoint=0
             do kcs=1,ncs
                call get_shell_size(bl(ictr),kcs,ksize)
                dmxs=0.0d0
                do my=1,isize
                   do lam=1,ksize
                      dmx=D(ipoint+my,kpoint+lam)
                      dmxs=max(dmxs,dmx)
                   end do     !lam in kcs
                end do     ! my in ics
                kpoint=kpoint+ksize
ckw
                dmxs=min(dmxs,1.d0)
c
                densmax=max(densmax,dmxs)
c
ckw             DS(ics,kcs)=dmxs
                DS(ics,kcs)=dmxs*dmxs
             end do     ! kcs
             ipoint=ipoint+isize
          end do     !ics
c
C       print statistics
      if(iprint.ge.2) call evaldmax_loc(ncs,DS,ids)
c
       end
C========dmaxij_new===============================================
      subroutine DmaxIJ_new(D,DS,ncf,ncs,nval,
     1                      nstron,inx,C,densmax,iprint,
     2                      dmaxi)
C       Calculates Max (over ij strong) of C(my,i) * C(lam,j)
C       as described in paper by Schutz, Hetzer, and Werner
C       Used for prescreening of integrals.
C       Output:
C       D(ncf,ncf) (this matrix is removed after return)
C         DS(ncs,ncs) same as D but over shells rather than
c         densmax - maximum element of "screening" density
c
C       individual functions
C       Input:
C       ncf - number of basis functions
C       ncs - number of contracted shells
C       nval - number of valence orbitals
C       nstron - logical array nstron(i,j) is true if i,j not a
C       strong pair
C       C(ncf,nval) matrix of MO coefficients (localized orbitals)
C
C       Svein Saebo, Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
       dimension D(ncf,ncf),C(ncf,nval),DS(ncs,ncs),dmaxi(ncf,*)
       dimension inx(*)
       logical*1 nstron(nval,nval)
      do lam=1,ncf
        do imo=1,nval
          dmaxj=0.0d0
          do jmo=1,nval
          if(nstron(imo,jmo)) cycle
            cj=abs(C(lam,jmo))
C           if(imo.ne.jmo) cj=cj*rfact(imo,jmo)
              dmaxj=max(dmaxj,cj)
          end do !j
          dmaxi(lam,imo)=dmaxj
        end do   !i
      end do     !lam
c
      do my=1,ncf
        do lam=1,my
          dmxx=0.0d0
          do imo=1,nval
            cprod=abs(C(my,imo)*dmaxi(lam,imo))
            dmxx=max(dmxx,cprod)
          end do  !imo
          D(my,lam)=dmxx
          D(lam,my)=dmxx
        end do    ! lam
      end do      ! my
CPP
       ids=0
C
C       print statstics
       if(iprint.ge.2)call evaldmax_loc(ncf,D,ids)
c
C  construct a matrix over contracted shells rather than basis function
c
          ictr=igetival('ictr')
c
       ids=1
          densmax=0.d0
          ipoint=0
          do ics=1,ncs
             call get_shell_size(bl(ictr),ics,isize)
             kpoint=0
             do kcs=1,ncs
                call get_shell_size(bl(ictr),kcs,ksize)
                dmxs=0.0d0
                do my=1,isize
                   do lam=1,ksize
                      dmx=D(ipoint+my,kpoint+lam)
                      dmxs=max(dmxs,dmx)
                   end do     !lam in kcs
                end do     ! my in ics
                kpoint=kpoint+ksize
ckw
                dmxs=min(dmxs,1.d0)
c
                densmax=max(densmax,dmxs)
c
ckw             DS(ics,kcs)=dmxs
                DS(ics,kcs)=dmxs*dmxs
             end do     ! kcs
             ipoint=ipoint+isize
          end do     !ics
c
C       print statistics
      if(iprint.ge.2) call evaldmax_loc(ncs,DS,ids)
c
       end
C===========evaldmax_loc=============================================
       subroutine evaldmax_loc(ncs,D,ids)
C       This subroutine evaluates the magnitudes of the elements of
C       the "Density" matrix used for pre-screening of integrals
C       These matrices are generated in Subroutine DmaxIJ.
C       Called from DmaxIJ
C        the subroutine can of course be used to evaluate any quadratic
C       matrix D of order ncs.
C
C       Svein Saebo , Fayetteville, AR  june 1999
C
       implicit real*8(a-h,o-z)
       dimension D(ncs,ncs)
       n1=0
       n2=0
       n3=0
       n4=0
       n5=0
       trx1=1.0d-1
       trx2=1.0d-2
       trx3=1.0d-3
       trx4=1.0d-4
       if(ids.eq.1) then
       trx1=trx1*trx1
       trx2=trx2*trx2
       trx3=trx3*trx3
       trx4=trx4*trx4
       endif
C       DS is squared !
       do imo=1,ncs
       do jmo=1,imo
       elem=D(imo,jmo)
       if (elem.gt.trx1) then    !  > 0.1
       n1=n1+1
       cycle
       endif
       if (elem.gt.trx2) then    !   >0.01
       n2=n2+1
       cycle
       endif
       if (elem.gt.trx3) then    !   >0.001
       n3=n3+1
       cycle
       endif
       if (elem.gt.trx4) then    !   >0.0001
       n4=n4+1
       cycle
       endif
       n5=n5+1
       end do
       end do
       ntot=ncs*(ncs+1)/2
       write(6,*)n1,' number of D elements > 0.1'
       write(6,*)n2,' number of D elements > 0.01 and < 0.1'
       write(6,*)n3,' number of D elements > 0.001 and < 0.01'
       write(6,*)n4,' number of D elements > 0.0001 and < 0.001'
       write(6,*)n5,' number of D elements < 0.0001'
       write(6,*) 'total number of elements: ', ntot
       end
C=======compab==lci utility=========================================
      SUBROUTINE COMPAB(A,B,NF,ID,JD,IFF,JFF)
C
C     **** THIS SUBROUTINE COMPRESSES A INTO B
C     **** A FULL MATRIX OF ORDER NCF
C     **** B ARRAY OF DIMENSION ID X JD
C     **** THE REDUCED BASIS SET ARE IN IFF(ID) AND JFF(JD)
C     Svein Sabeo Fayetteville, AR may 1999
C     From original lci program, slightly modified
C
      IMPLICIT REAL*8(A-H,O-Z)
       integer*2 iff(*),jff(*)
      DIMENSION A(NF,NF),B(ID,JD)
      DO  JJ=1,JD
      DO  II=1,ID
      B(II,JJ)=A(iff(ii),jff(jj))
       enddo
       enddo
      END
C=======prtst===========================================
       subroutine prtst(itst,npairs)
       implicit real*8(a-h,o-z)
       dimension itst(*)
       write(6,*) (itst(ixy),ixy=1,npairs)
       end
C==========remzero===============================
      subroutine remZero(rkmyla,icmp,epsi,nval,my,
     1                   lam,jcm,iprint,nrel,nrec,
     2                   nwo,ndisk,ncore,icol2,jfun3,
     3                   ncol2,transin2,nrcpf)
C
C      finds non-zero elements in matrix, and put this element in buffer
C      followed by two words: MY,nY and imo,jmo
C      RKmyla - matrix to be evaluated: dimension nval*nval
C      Thresh integral threshold to be used for unscaling
C      of the matrix elements.
C      This subroutine is called from tra1Shel once for each
C     my,lam (all)
C        thresh mp2 integrals threshold
C       thres2 threshold for neglecting a contribution
C
C      Svein Saebo Fayetteville, AR May 1999
C
       implicit real*8(a-h,o-z)
       dimension rkmyla(jcm,*)
       integer*2 iindex(2),itest(2)
       dimension jfun3(jcm+1,*),icol2(*)
       dimension icmp(nrel)
       equivalence (longint,iindex(1))
        equivalence (jtest,itest(1))
C
       data dblemax /2 147 483 648.0d0/
       icount=0
C
C       do jmo=1,nval       
C       do imo=1,nval
       do jj=1,ncol2
       jmo=icol2(jj)
       icount=icount+jfun3(1,jmo)-1
       do ii=1,jfun3(1,jmo)-1
       imo=jfun3(ii+1,jmo)
        if(jmo.gt.imo.and.my.eq.lam) cycle
       elem=rkmyla(ii,jj)
       if(abs(elem).gt.epsi) then
C
       if(abs(elem).gt.dblemax) then
       write(6,*) 'elem=',elem
       call nerror(1,'remzero','integral too large',imo,jmo)
       end if
C
       if(nwo+3.ge.nrel) then
C       buffer full write it out
       icmp(nwo+1)=-1
       nrec=nrec+1
          nfil=(nrec-1)/nrcpf
          nnrc=nrec-nfil*nrcpf
          nndsk=ndisk+nfil
       write(unit=nndsk,rec=nnrc) icmp
       nwo=0
       endif
       iindex(1)=my
       iindex(2)=lam
       nwo=nwo+1
       icmp(nwo)=longint
       itest(1)=imo
       itest(2)=jmo
       nwo=nwo+1
       icmp(nwo)=jtest
       nwo=nwo+1
       icmp(nwo)=elem
       endif
       end do
       end do
       transin2=transin2+icount
       end
C=========Phase1===============================================
      subroutine Phase1(ncf,nval,ncore,npairs,ndisk,
     1                  ndisk2,nrec, listp,lbin2, nbins,
     2                  ibins,icounter,ibinctr,ibuf,lrec,
     3                  idm,lbas,nstron,invl,iprint,
     4                  ifunpair,isorb,uniq,nsym,unorb,
     5                  full,inx,lalb,lbasij,idmij,
     6                  inndp,nud,frac,nfiles,nrcpf,
     7                  nrcpb,gigabyte)
C
C       Bin-Sort Phase one.  Based on Subroutine BinSort1 written
C       by Peter Pulay.
C       This subroutine reads the file containing the halftransformed
C       integrals (stored as compact KMYLAM(i,j))
C       An integral KMYLAM(I,J) is put into bin
C       IJ with indices My and Lam.
C       Arguments:
C       ncf - number of basis functions
C       nval - number of valence orbitals
C       npairs - number of pairs included in the transformation
C                = number of symmetry unique strong pairs
C       ndisk, ndisk2  unit numbers
C       nrec - number of records on the half-trasnformed file
C         nbin - the number of bins written so far
C       listp - pair list
C       ibuf - buffer with halftransformed intgrals
C       lrec - record length in words for file with halftransformed
C       integrals
C       nstron - quadratic logical*1 array of orer nval. nstron(i,j)
C                  is true if ij is not a strong pair.
C       invl - inverse pair list
C
C       Svein Saebo Fayetteville Ar, June 1999
C
      implicit real*8 (a-h,o-z)
      logical full,lalb
      logical*1 nstron(nval,*)
      integer*2 idm(*),lbas(ncf,*),en,null
      integer*2 lbasij(ncf,*),idmij(*),inndp(*)
      integer*2 iindex(2),jindex(2),itest(2)
      integer*4 ibins(lbin2,*)
      dimension ibuf(lrec),invl(nval,*),listp(*)
      dimension icounter(*)
      dimension ibinctr(*)
      dimension ifunpair(7,*),isorb(nval,*)
      integer uniq(*),unorb(2,*)
      dimension inx(12,*)
      equivalence(longint,iindex(1)),(longint1,jindex(1))
      equivalence(jtest,itest(1))
      parameter(null=0,en=1)
      parameter(zero=0.0d0)
      write(6,*) 'Phase1'
      nvir=ncf-nval-ncore
      write(6,*) ncf,nval,nvir,lrec,lbin2
C
C       Zero out array with local domains
       do ix=1,nval
       do  iy=1,ncf
       lbas(iy,ix)=null
       end do
       end do
       do iy=1,nud
       do iz=1,ncf
       lbasij(iz,iy)=null
       enddo
       enddo
C       zero out word-counter for bins
      call izeroit(icounter,npairs)
      nbin=0
      do irec=1,nrec
c read one record and process it, it should contain
c several pairs
      ifil=(irec-1)/nrcpf
      krec=irec-ifil*nrcpf
      nndsk=ndisk+ifil
      read(nndsk,rec=krec)ibuf
      nwo=1
C     return here for more pairs in the buffer
  100 continue
         if(ibuf(nwo).eq.-1) cycle
C       end of information in buffer should be -1
C       end of record
       if(ibuf(nwo).eq.0) call nerror(96,'phase1','error',0,0)
C       this should ot happen..
C       Start taking information from buffer
       longint=ibuf(nwo)
       my=iindex(1)
       lam=iindex(2)
       jindex(1)=iindex(2)
       jindex(2)=iindex(1)
       nwo=nwo+1
       jtest=ibuf(nwo)
       ii=itest(2)
       jj=itest(1)
C   IJ will now only go over symmetry unique strong pairs
C*****  Set initial local Domains************************
       if(lalb) then
       if(ii.ge.jj) then
       lbas(my,ii)=en
       lbas(lam,jj)=en
       endif
       if(jj.ge.ii) then
       lbas(my,jj)=en
       lbas(lam,ii)=en
       endif
       else
C******** this is the default initial local domain
       if(ii.eq.jj) then
       lbas(my,ii)=en
       lbas(lam,ii)=en
       endif
       endif
C******* Local Domains finished*************************
          nwo=nwo+1
         iijj=invl(ii,jj)
         ijs=listp(iijj)
       if(ii.ne.jj) then
       lpa=inndp(ijs)
       lbasij(my,lpa)=en
       lbasij(lam,lpa)=en
       endif
       ij=uniq(ijs)
       if(ij.le.0) goto 1000
C
C       Start putting integrals into bins:
C       ijs       strong pair list
C       ij       symmetry unique pair list       
C              
C       always check if bin is full first
C       
C
           if(icounter(ij).ge.lbin2) then
              nbin=nbin+1
           nbnfil=(nbin-1)/nrcpb
           nnbn=nbin-nbnfil*nrcpb
              ibinctr(nbin)=ij
              call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
            end if              !bin full
C
            icounter(ij)=icounter(ij)+2
       if(icounter(ij).gt.lbin2)call nerror(1,'binsloc','overflow'
     1 ,ij,icounter(ij))
            ibins(icounter(ij)-1,ij)=ibuf(nwo)
       if(my.ne.lam) then
        if(ii.gt.jj) then
            ibins(icounter(ij),ij)=longint
            endif              !iigtjj
        if(jj.gt.ii) then
            ibins(icounter(ij),ij)=longint1
            endif                    !jjgtii
               if(ii.eq.jj) then
            ibins(icounter(ij),ij)=longint
C       check if bin is full
           if(icounter(ij).ge.lbin2) then
              nbin=nbin+1
           nbnfil=(nbin-1)/nrcpb
           nnbn=nbin-nbnfil*nrcpb
              ibinctr(nbin)=ij
              call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
            end if              !bin full
C
           icounter(ij)=icounter(ij)+2
      if(icounter(ij).gt.lbin2)call nerror(2,'binsloc','overflow'
     1 ,ij,icounter(ij))
            ibins(icounter(ij)-1,ij)=ibuf(nwo)
            ibins(icounter(ij),ij)=longint1
          endif            !iieqjj
          else                  !my.ne.lam
          ibins(icounter(ij),ij)=longint
           endif                  !my.ne.lam
 1000 continue
C      symmetry handelling
      if(nsym.gt.0) then
      do isym=1,nsym
      iip=isorb(ii,isym)
      jjp=isorb(jj,isym)
      my1=ifunpair(isym,my)
      lam1=ifunpair(isym,lam)
C
C      ichs=iip*jjp*my1*lam1
C      jsign=isign(1,ichs)
      ichs1=iip*my1
      jsig1=isign(1,ichs1)
      ichs2=jjp*lam1
      jsig2=isign(1,ichs2)
      jsign=jsig1*jsig2
C
      iip2=iabs(iip)
      jjp2=iabs(jjp)
      my2=iabs(my1)
      lam2=iabs(lam1)
      iindex(1)=my2
      iindex(2)=lam2
      jindex(1)=lam2
      jindex(2)=my2
C
      if(nstron(iip2,jjp2)) then
      write(6,*) 'weak pair:',iip2,jjp2,isym,ii,jj,nwo
      call nerror(96,'phase1','error',0,0)
      endif
C      local domains *********************************
      if(lalb) then
      if(iip2.ge.jjp2) then
      lbas(my2,iip2)=en
      lbas(lam2,jjp2)=en
      endif
      if(jjp2.ge.iip2) then
      lbas(my2,jjp2)=en
      lbas(lam2,iip2)=en
      endif
          else
      if(iip2.eq.jjp2) then
      lbas(my2,iip2)=en
      lbas(lam2,iip2)=en
      endif
      endif
C  *********  end local domains*****************
        iijj=invl(iip2,jjp2)
        ijs=listp(iijj)
      if(ijs.le.0.or.ijs.ge.nval*nval) then
      write(6,*) iip2,jjp2,ijs
        call nerror(1,'phase1','error',0,0)
                    endif
      ij=uniq(ijs)
      if(iip2.ne.jjp2) then
      lpa=inndp(ijs)
      lbasij(my2,lpa)=en
      lbasij(lam2,lpa)=en
      endif
      if(ij.le.0) cycle
C      check if bin is full
           if(icounter(ij).ge.lbin2) then
              nbin=nbin+1
          nbnfil=(nbin-1)/nrcpb
          nnbn=nbin-nbnfil*nrcpb
              ibinctr(nbin)=ij
              call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
            end if            !bin full
C
            icounter(ij)=icounter(ij)+2
            ibins(icounter(ij)-1,ij)=jsign*ibuf(nwo)
      if(my2.ne.lam2) then
       if(iip2.gt.jjp2) then
            ibins(icounter(ij),ij)=longint
           endif
       if(jjp2.gt.iip2) then
            ibins(icounter(ij),ij)=longint1
           endif
              if(iip2.eq.jjp2) then
            ibins(icounter(ij),ij)=longint
C      check if bin is full
           if(icounter(ij).ge.lbin2) then
              nbin=nbin+1
          nbnfil=(nbin-1)/nrcpb
          nnbn=nbin-nbnfil*nrcpb
              ibinctr(nbin)=ij
              call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
            end if
C
           icounter(ij)=icounter(ij)+2
            ibins(icounter(ij)-1,ij)=jsign*ibuf(nwo)
            ibins(icounter(ij),ij)=longint1
          endif
          else
          ibins(icounter(ij),ij)=longint
            endif
          end do
          end if
          nwo=nwo+1
C   go back for more integrals in buffer
      if(nwo.lt.lrec)      goto 100
      end do      !  irec loop
c  write out the content of the bins at the end
      do ij=1,npairs
           if(icounter(ij).eq.lbin2) then
              nbin=nbin+1
          nbnfil=(nbin-1)/nrcpb
          nnbn=nbin-nbnfil*nrcpb
              ibinctr(nbin)=ij
              call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
            end if
c  zero out the rest
        if(icounter(ij).gt.0) then
          nbin=nbin+1
          nbnfil=(nbin-1)/nrcpb
          nnbn=nbin-nbnfil*nrcpb
          ibinctr(nbin)=ij
      icounter(ij)=icounter(ij)+2
            ibins(icounter(ij)-1,ij)=0
            ibins(icounter(ij),ij)=0
C      put a 0 for labels to indicate end of information
C         do j=icounter(ij)+1,lbin2
c           ibins(j,ij)=0
C           ibins(j,ij)=0
C         end do
          call WriteBin1(ndisk2+nbnfil,nnbn,lbin2,ibins(1,ij))
              icounter(ij)=0
        end if
      end do
      if(nbin.gt.nbins) then
      write(6,111)nbins,nbin
  111 format(1x,'******Number of bins greater than estimate',/,
     1'estimate= ',1I10,'actual= ',1I10, 'stopping........ ')
      call nerror(96,'phase1','error',0,0)
      endif
      nbins=nbin
C      determine domain dimensions
C      and final ibas
      if(full) then
      do imo=1,nval
      do iao=ncore+1,ncf
      lbas(iao,imo)=en
      enddo
      enddo
      endif
      if(iprint.ge.3) then
          write(6,*) '***** Local Domains before checking*******'
      endif
      dimtot=zero
      do imo=1,nval
      if(unorb(1,imo).gt.0) then
      call  lbshel(lbas(1,imo),inx)
      iao=0
      idd=0
      do iy=1,ncf
      if(lbas(iy,imo).eq.en) then
      iao=iao+1
      idd=idd+1
      lbas(iao,imo)=iy
      endif
      enddo
      idm(imo)=idd
          if(iprint.ge.3) then
             write(6,*) 'Local domain for Valence Orbital no.: ',imo
             write(6,*) 'Dimension: ',idd
             write(6,*) 'AOs: ',(lbas(ibo,imo),ibo=1,idd)
          end if
          idx=min0(idd,nvir)
      dimtot=dimtot+idx*idx
      endif
      enddo
C      initial basis for non-diagonal pairs
      do ipp=1,nud
      iao=0
      iddm=0
      do iy=1,ncf
      if(lbasij(iy,ipp).eq.en) then
      iao=iao+1
      iddm=iddm+1
      lbasij(iao,ipp)=iy
      endif
      enddo
      idmij(ipp)=iddm
          iddx=min0(iddm,nvir)
      dimnd=iddx*iddx
      dimtot=dimtot+dimnd*frac
      enddo
      if(nsym.gt.0)
     1  call symchec(unorb,ifunpair,lbas,idm,nval,ncf,iprint)
      write(6,*) 'Initial Total local dimension: ',dimtot
      dimbt=dimtot*8
      nfiles=dimbt/gigabyte+1
          write(6,*) 'number of files needed for kij: ',nfiles,dimbt
          if(nfiles.gt.9) then
                  Write(6,*)'warning no of files=', nfiles
C         call nerror(1,'phase1','too many files',nfiles,dimbt)
      endif
      end
C===========dtest=======================================
        subroutine dtest(dmxx,ixs,kxs,DS,ncs)
      implicit real*8(a-h,o-z)
      dimension DS(ncs,ncs)
      dmxx=DS(ixs,kxs)
      end
C=============phase2=============================================
      Subroutine Phase2(ncf,nval,ndisk2,ndiskx,lbin2,
     1                  nbins,nstron,ibinctr,lbin,thresh,
     2                  ncore,i1bin,xmat,idm,lbas,
     2                  uniq,isorb,unorb,ifunpair,npars,
     4                  nsym,iprint,thres3,inx,lbasij,
     5                  idmij,equi,inndp,redu,nud,
     6                  ncda,ttda,remf,wrdw,ifil,
     7                  nrcpb,giga)
C
C      binsort phase 2
C
C      Svein Saebo, Fayetteville, AR summer 1999
C

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30000)
      dimension xmat(*),ibinctr(*),inx(12,*)
      integer*2 idm(*),lbas(*)
      integer*2 lbasij(*),idmij(*),inndp(*)
      integer*4 i1bin(2,*)
      dimension isorb(nval,*),ifunpair(7,*)
      logical*1 nstron(nval,*)
      logical ndiag,redu,remf,red2
      integer uniq(*),unorb(2,*),equi(2,*)
C
C      calculate average local dimension
C
C      threshold for removing redundant basis functions
      epsix=1.0d-07
      tpro=0.0d0
      red2=.not.redu.and.remf
C
      if(iprint.ge.3) then
      write(6,*)'******Local Domains after checking******'
      write(6,*) 'Thresholds = ',thres3,epsix
      endif
      ndiag=.false.
      lpt=0
      ijs=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      ijs=ijs+1
C      diagonal pairs first
      if(ii.ne.jj) cycle
      ij=uniq(ijs)
      if(ij.le.0) cycle
      call ReadBin_loc(ncf,nval,ij,ndisk2,lbin2,
     1                 nbins,ibinctr,lbin,thresh,i1bin,
     2                 xmat,xmat,ndiag,nrcpb)
      if(redu) then
      idd=idm(ii)
      ifu=ncf*(ii-1)+1
      call matdef('xlox','q',idd,idd)
      ixlox=mataddr('xlox')
      call compab(xmat,bl(ixlox),ncf,idd,idd,
     1                lbas(ifu),lbas(ifu))
      wrdw=wrdw+idd*idd+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=idd*idd+idd
      endif
      ndisk4=ndiskx+ifil
      call wriarr(bl(ixlox),idd*idd,ndisk4)
      call matrem('xlox')
      else
      call writeK3(ncf,idm,lbas,ii,ndiskx,
     1             nval,iprint,thres3,inx,unorb,
     2             epsix,tpro,remf,ncda,ttda,
     3             ifil,wrdw,giga)
      endif
        end do
      end do
      if(nsym.gt.0) then
      call symchec(unorb,ifunpair,lbas,idm,nval,
     1             ncf,iprint)
      endif
      write(6,*)'Diagonal pairs:'
      write(6,*) ncda,ttda/60.0
C
C      construct local bases  for non-diagonal pairs
C
C      save initial local basis for non-diagonal pairs
c     call getmem((nud*ncf)/4+1,lbijt)
      call getint_2(nud*ncf,lbijt)
c     call getmem(nud/4+1,idijt)
      call getint_2(nud,idijt)
      call bascop(lbasij,idmij,bl(lbijt),bl(idijt),nud,ncf)
c     call getmem(ncf/4+1,lbmap)
      call getint_2(ncf,lbmap)
      call mklbmap(lbas,idm,bl(lbmap),ncf,nval,
     1             lbasij,idmij,nstron)
      call retmem(1)
C      remove redundancies
       if(red2) then
       call lbndia(ncf,nval,nstron,iprint,uniq,
     1                ifunpair,nsym,lbasij,idmij,epsix,
     2                equi,inndp,ncda,ttda)
      else
      if(nsym.gt.0)
     1  call symche2(equi,ifunpair,lbasij,idmij,nval,
     2               ncf,iprint,nstron,inndp)
      endif
      ijs=0
      ndiag=.true.
      lpt=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      ijs=ijs+1
C      then nondiagonal pairs
      if(ii.eq.jj) cycle
      lpt=lpt+1
      ij=uniq(ijs)
      if(ij.le.0) cycle
      call ReadBin_loc(ncf,nval,ij,ndisk2,lbin2,
     1                 nbins,ibinctr,lbin,thresh,i1bin,
     2                 xmat,xmat,ndiag,nrcpb)
      if(redu) then
      idd=idmij(lpt)
      ifu=ncf*(lpt-1)+1
      call matdef('xlox','q',idd,idd)
      ixlox=mataddr('xlox')
      call compab(xmat,bl(ixlox),ncf,idd,idd,
     1                lbasij(ifu),lbasij(ifu))
      wrdw=wrdw+idd*idd+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=idd*idd+idd
      endif
      ndisk4=ndiskx+ifil
      call wriarr(bl(ixlox),idd*idd,ndisk4)
      call matrem('xlox')
      else
      call writeK2(ncf,ndiskx,idmij,lbasij,lpt,
     1             tpro,bl(lbijt),bl(idijt),wrdw,ifil,
     2             giga)
      endif
        end do
      end do
      call retmem(2)
      write(6,*) 'Projection: ',tpro/60.0d0
      end
c=======================================================
      subroutine normv(v,nvir,ncf)
      implicit real*8(a-h,o-z)
      dimension v(ncf,nvir)
      do iv=1,nvir
      sum=0.0d0
      do ia=1,ncf
      sum = sum +v(ia,iv)**2
      end do
      s =sqrt(sum)
      do ia=1,ncf
      v(ia,iv)=v(ia,iv)/s
      end do
      end do
      end
C===============================================
      subroutine checore(ics,tab,ncore,icycle)
      implicit real*8(a-h,o-z)
      dimension tab(*)
      icycle=0
      do lx=1,ncore
      itbb=tab(lx)
      if(ics.eq.itbb) then
      icycle=1
      exit
      end if
      end do
      end
C========================
      logical function lastaop(nsym,iao,jao,ifunpair)
      dimension ifunpair(7,*)
      lastaop=.false.
      if(nsym.eq.0) return
      do isym=1,nsym
      iaop=ifunpair(isym,iao)
      jaop=ifunpair(isym,jao)
      iaop1=iabs(iaop)
      jaop1=iabs(jaop)
      if(iaop1.gt.iao.or.jaop1.gt.jao) then
      lastaop=.true.
      goto 100
      else
      iii=max(iaop1,jaop1)
      jjj=min(iaop1,jaop1)
      endif
      if(iii.gt.iao.or.iii.eq.iao.and.jjj.gt.jao) then
      lastaop=.true.
      goto 100
      endif
      enddo
C      last pair
      write(6,*) 'pair retained', iao,jao
      return
  100 continue
C      pair neglected
      write(6,*) 'pair neglected',iao,jao
      end
      logical function lastao(nsym,iao,ifunpair)
      dimension ifunpair(7,*)
      lastao=.false.
      if(nsym.eq.0) return
      do isym=1,nsym
      iaop=ifunpair(isym,iao)
      iaop1=iabs(iaop)
      if(iaop1.gt.iao) then
      lastao=.true.
      return
      endif
      end do
C      last ao
      end
c=========lbshel===================================================
      subroutine lbshel(lbas,inx)
      implicit integer(a-z)
      dimension inx(12,*)
      integer*2 lbas(*),en
C
      en=1
      ncs=igetival('ncs')
      idd=0
      do ics=1,ncs
      ishel=0
      ifirst=inx(11,ics)+1
      ilast=inx(10,ics)
      do icf=ifirst,ilast
      if(lbas(icf).eq.en) then
      ishel=1
C      shell ics should be included
      exit
      endif
      enddo
C      Make sure shell ics is included
      if(ishel.gt.0) then
      do icf =ifirst,ilast
      lbas(icf)=en
      idd=idd+1
      enddo
      endif
      enddo
C      write(6,*) ' new dimension: ',idd
      end
C==========ilash=========================================
      integer function Ilash(ics,inx)
      integer inx(12,*)
      ilash=inx(10,ics)
      return
      end
C======writeK2===================================================
      subroutine writeK2(ncf,ndiskx,idmij,lbasij,lpt,
     1                   tpro,lbsijt,idmijt,wrdw,ifil,
     2                   giga)
C      
C      projects and
C      convert a full Kij to local basis and writes it out on disk
C      unprojected (sparce but full) Kij is in 'xmat'
C      Since unsymmetrical local bases are used K-0.5K(tr) is
C      calculated before contraction to local basis.
C            
C      Svein Saebo  Starkville, MS Dec 1999
C            

      use memory

      implicit real*8(a-h,o-z)
      integer*2 idmij(*),lbasij(*),lbsijt(*),idmijt(*)
c     common/big/bl(30000)
C
      ixadr=mataddr('xmat')
      iprad=mataddr('proj')
C
      idd=idmij(lpt)
      ifu=ncf*(lpt-1)+1
      ijdm=idd*idd
      iddi=idmijt(lpt)
      call matdef('rkc','q',idd,idd)
      ikija=mataddr('rkc')
C
      call secund(ttpro1)
      call matdef('xloc','q',iddi,iddi)
      ixloa=mataddr('xloc')
      call compab(bl(ixadr),bl(ixloa),ncf,iddi,iddi,
     1 lbsijt(ifu),lbsijt(ifu))
C
      call matdef('pki','r',iddi,idd)
      call matdef('pfj','r',idd,iddi)
      call matdef('pfi','r',iddi,idd)
      ipfia=mataddr('pfi')
      call compab(bl(iprad),bl(ipfia),ncf,iddi,idd,
     1 lbsijt(ifu),lbasij(ifu))
      call matmmult('xloc','pfi','pki')      
      call matpose2('pfi','pfj','n')
      call matrem('pfi')
      call matmmult('pfj','pki','rkc')
      call matrem('pfj')
      call matrem('pki')
      call secund(ttpro2)
      tpro=tpro+ttpro2-ttpro1
      wrdw=wrdw+ijdm+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=ijdm+idd
      endif
      ndisk4=ndiskx+ifil
      call wriarr(bl(ikija),ijdm,ndisk4)
      call matrem('xloc')
      call matrem('rkc')
      end
C======writeK3===================================================
      subroutine writeK3(ncf,idm,lbas,ii,ndiskx,
     1             nval,iprint,thres3,inx,unorb,
     2             epsix,tpro,remf,ncda,ttda,
     3             ifil,wrdw,giga)
C      
C      projects and
C      convert a full Kij to local basis and writes it out on disk
C      unprojected (sparce but full) Kij is in 'xmat'
C      Since unsymmetrical local bases are used K-0.5K(tr) is
C      calculated before contraction to local basis.
C      file xxx.kij will thus contain K-0.5K(tr) in projected local
C      basis.
C            
C      Svein Saebo  Starkville, MS Dec 1999
C            

      use memory

      implicit real*8(a-h,o-z)
      integer*2 idm(*),lbas(*)
      dimension inx(12,*)
      logical remf
      integer unorb(*)
c     common/big/bl(30000)
C
      ixadr=mataddr('xmat')
      iprad=mataddr('proj')
C
      idd=idm(ii)
      ifu=ncf*(ii-1)+1
      call secund(ttpro1)
      call matdef('rkc','q',idd,idd)
      ikija=mataddr('rkc')
      call matdef('xloc','q',idd,idd)
      ixloa=mataddr('xloc')
      call compab(bl(ixadr),bl(ixloa),ncf,idd,idd,
     1 lbas(ifu),lbas(ifu))
C
      call matdef('pfi','q',idd,idd)
      call matdef('pki','q',idd,idd)
      ipfia=mataddr('pfi')
      call compab(bl(iprad),bl(ipfia),ncf,idd,idd,
     1 lbas(ifu),lbas(ifu))
      call matmmult('xloc','pfi','pki')      
      call matpose('pfi')
      call matmmult('pfi','pki','rkc')
      call matrem('pki')
      call matrem('pfi')
      call secund(ttpro2)
      tpro=tpro+ttpro2-ttpro1
C      projected Kij in rkc for this pair, check if dimensions
C      can be reduced.
c     call getmem(ncf/4+1,lbaso)
      call getint_2(ncf,lbaso)
c     call getmem(idd/4+1,lbcom)
      call getint_2(idd,lbcom)
      call getint(idd,jrem)
C      call remzerk(bl(ikija),bl(ikija),ncf,idd,idm,
C    1                 lbas,ii,nval,bl(jrem),iprint,
C    2                 thres3,inx)
       call RemZerK2(bl(ikija),idd,lbas,nval,bl(jrem),
     1                 iprint,thres3,inx,bl(lbcom),bl(lbaso),
     2                 idm,ii,ifu,ncf)
      call retmem(3)
      ijdm=idd*idd
      if(remf) then
C
C      remove redundencies NOW!
c     call getmem(ncf/4+1,lbaso)
      call getint_2(ncf,lbaso)
      call lbdiag(ncf,nval,lbas(ifu),idd,iprint,
     1                unorb,ifunpair,nsym,bl(lbaso),idmo,
     2                ii,epsix,ncda,ttda)
      idm(ii)=idd
c     call getmem(idd/4+1,lbcom)
      call getint_2(idd,lbcom)
      call regenK2(idd,lbas(ifu),idmo,bl(lbaso),bl(ikija),
     1                ndiskx,bl(lbcom),ifil,wrdw,giga)
      call retmem(1)
C      call regenK(ncf,nval,idm,lbas(ifu),idmo,
C    1                bl(lbaso),nstron,ii,bl(ikija),ndisk4,
      call retmem(1)
      else
      wrdw=wrdw+ijdm+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=ijdm+idd
      endif
      ndisk4=ndiskx+ifil
      call wriarr(bl(ikija),ijdm,ndisk4)
      endif
      call matrem('xloc')
C
      call matrem('rkc')
      end
C=======remzerk=========================================
      subroutine RemZerK(RK,rk2,ncf,iddd,idm,
     1                       lbas,ixx,nval,jrem,iprint,
     2                       epsix,inx)
C
C      remove additional negligible elements after projection
C
C      Svein Saebo, Faytteville, AR summer 1999
C
      use memory

      implicit real*8(a-h,o-z)
      dimension RK(iddd,iddd),rk2(*),jrem(*),inx(12,*)
      integer*2 idm(*),lbas(ncf,nval)
c     common/big/bl(30000)
C
      idd=iddd
      ix=ixx
C
      ic=0
      call matdef('ktfu','q',ncf,ncf)
      ktfua=mataddr('ktfu')
      call matzero('ktfu')
      call scatad(bl(ktfua),RK,ncf,idd,idd,lbas(1,ix),lbas(1,ix))
      do jj=1,iddd
      izero=0
      do ii=1,iddd
      if(abs(RK(ii,jj)).gt.epsix) then
      izero=1
      exit
      endif
      end do !over ii
C
      if(izero.eq.0)then
C      column zero
C      remove function no jj from lbas
      ic=ic+1
      jrem(ic)=jj
      endif
      end do   !over jj
C
C      remove the ic functions from lbas
      do krem=ic,1,-1
      jst=jrem(krem)
      if(jst.eq.idd) then
      idd=idd-1
      else
      idd=idd-1
      do irem=jst,idd
      lbas(irem,ix)=lbas(irem+1,ix)
      end do
      endif
      enddo   !over krem
C
      idm(ix)=idd
c     call getmem(ncf/4+1,lbasn)
      call getint_2(ncf,lbasn)
      call fixshel(lbas(1,ix),bl(lbasn),idm,ix,ncf,nval,inx)
      idd=idm(ix)
      call retmem(1)
      call compab(bl(ktfua),rk2,ncf,idd,idd,lbas(1,ix),lbas(1,ix))
      if(iddd.ne.idd.and.iprint.ge.3) then
      write(6,*) 'new local dimension for :',ix,' - ',ic,idm(ix)
      endif
       call matrem('ktfu')
      iddd=idd
      end
C=======remzerk2========================================
      subroutine RemZerK2(RK,iddd,lbas,nval,jrem,
     1                       iprint,epsix,inx,lbcom,lbaso,
     2                       idm,ix,ifu,ncf)
C
C      remove additional negligible elements after projection
C
C      Svein Saebo, Faytteville, AR summer 1999
C      modified by SS fall 2000
C

      use memory

      implicit real*8(a-h,o-z)
      dimension RK(iddd,iddd),jrem(*),inx(12,*)
      integer*2 idm(*),lbas(*),lbcom(*),lbaso(*)
c     common/big/bl(30000)
C
      idmo=iddd
      idd=iddd
      ic=0
      do jj=1,iddd
      izero=0
      do ii=1,iddd
      if(abs(RK(ii,jj)).gt.epsix) then
      izero=1
      exit
      endif
      end do !over ii
C
      if(izero.eq.0)then
C      column zero
C      remove function no jj from lbas
      ic=ic+1
      jrem(ic)=jj
      endif
      end do   !over jj
      if(ic.eq.0) return
C
      do imv=1,iddd
      lbaso(imv)=lbas(ifu-1+imv)
      enddo
C
C      remove the ic functions from lbas
      do krem=ic,1,-1
      jst=jrem(krem)
      if(jst.eq.idd) then
      idd=idd-1
      else
      idd=idd-1
      do irem=jst,idd
      lbas(ifu-1+irem)=lbas(ifu+irem)
      end do
      endif
      enddo   !over krem
      if(idd.le.0)call nerror(0,'writek3','dim zero',idd,ii)
C
      idm(ix)=idd
c     call getmem(ncf/4+1,lbasn)
      call getint_2(ncf,lbasn)
      call fixshel(lbas(ifu),bl(lbasn),idm,ix,ncf,nval,inx)
      idd=idm(ix)
      call retmem(1)
      if(idd.eq.idmo) return
      call regenK3(idd,lbas(ifu),idmo,lbaso,RK,
     1             lbcom)
      if(iddd.ne.idd.and.iprint.ge.3) then
      write(6,*) 'new local dimension for :',ix,' - ',ic,idm(ix)
      endif
      iddd=idd
      end
C=========regenK============================================
      subroutine regenK(ncf,nval,idm,lbas,ido,
     1                  lbaso,nstron,ii,rkij,nkunix,
     2                  ifil,wrdw,giga)

      use memory

      implicit real*8(a-h,o-z)
      dimension rkij(*)
      integer*2 idm(*),lbas(*),lbaso(*)
      logical*1 nstron(nval,nval)
c     common /big/bl(30000)
      call matdef('kfull','q',ncf,ncf)
      kfadr=mataddr('kfull')
      ijo=ido*ido
      idd=idm(ii)
      ijn=idd*idd
      call matdef('kij','q',idd,idd)
      call matzero('kfull')
      kladr=mataddr('kij')
      call scatad(bl(kfadr),rkij,ncf,ido,ido,
     1 lbaso,lbaso)
        call compab(bl(kfadr),bl(kladr),ncf,idd,idd,lbas,
     1 lbas)      
      wrdw=wrdw+ijn+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=ijn+idd
      endif
      nkuni2=nkunix+ifil
      call wriarr(bl(kladr),ijn,nkuni2)
      call matrem('kij')
      call matrem('kfull')
      end
C=========regenK2===========================================
      subroutine regenK2(idd,lbas,ido,lbaso,rkij,
     1                  nkunix,lbdif,ifil,wrdw,giga)

      use memory

      implicit real*8(a-h,o-z)
      dimension rkij(*)
      integer*2 lbas(*),lbaso(*),lbdif(*)
c     common /big/bl(30000)
      icom=0
      iny=1
      do ixx=1,ido
      iao=lbaso(ixx)
      nao=lbas(iny)
      if(iao.eq.nao) then
      icom=icom+1
      iny=iny+1
      lbdif(icom)=ixx
      endif
      if(icom.eq.idd) exit
      enddo
      call matdef('kij','q',idd,idd)
      kladr=mataddr('kij')
        call compab(rkij,bl(kladr),ido,idd,idd,lbdif,
     1 lbdif)
      ijn=idd*idd
      wrdw=wrdw+ijn+idd
      if(wrdw.gt.giga) then
      ifil=ifil+1
      wrdw=ijn+idd
      endif
      nkuni2=nkunix+ifil
      call wriarr(bl(kladr),ijn,nkuni2)
      call matrem('kij')
      end
C=========regenK3===========================================
      subroutine regenK3(idd,lbas,ido,lbaso,rkij,
     1                  lbdif)

      use memory

      implicit real*8(a-h,o-z)
      dimension rkij(*)
      integer*2 lbas(*),lbaso(*),lbdif(*)
c     common /big/bl(30000)
      icom=0
      iny=1
      do ixx=1,ido
      iao=lbaso(ixx)
      nao=lbas(iny)
      if(iao.eq.nao) then
      icom=icom+1
      iny=iny+1
      lbdif(icom)=ixx
      endif
      if(icom.eq.idd) exit
      enddo
      call matdef('kij','q',idd,idd)
      kladr=mataddr('kij')
        call compab(rkij,bl(kladr),ido,idd,idd,lbdif,
     1 lbdif)
      iddsq=idd*idd
      call move(bl(kladr),rkij,iddsq)
      call matrem('kij')
      end
C==========fixshel==========================================
      subroutine fixshel(lbas,lbasn,idm,imo,ncf,nval,inx)
      implicit integer(a-z)
      integer*2 lbas(*),lbasn(*),idm(*),null,en
      parameter(en=1,null=0)
      idd=idm(imo)
      do ize=1,ncf
      lbasn(ize)=null
      enddo
      do iset=1,idd
      jfu=lbas(iset)
      lbasn(jfu)=en
      enddo
      call lbshel(lbasn,inx)
C      set new lbas
      idny=0
      do iche=1,ncf
      if(lbasn(iche).eq.null) cycle
      idny=idny+1
      lbas(idny)=iche
      enddo
      idm(imo)=idny
      end
C========lbdiag================================================
      subroutine lbdiag(ncf,nval,lbasi,idd,iprint,
     1                      unorb,ifunpair,nsym,lbaso,idmo,
     2                      ii,epsx,ncda,ttda)
C
C      constrcuts projection matrix and remove redundant functions
C      for diagonal pairs
C
C      Svein Saebo, Fayetteville, AR Summer 2000
C

      use memory

        implicit real*8(a-h,o-z)
      dimension ifunpair(7,*)
      logical*1 nstron(nval,nval)
      integer*2 lbasi(*),lbaso(*)
      integer unorb(2,*)
c     common/big/bl(30000)
      parameter(zero=0.0d0)
C
      do iix=1,idd
      lbaso(iix)=lbasi(iix)
      enddo
      idmo=idd
C      Remove redundant basis functions
      if(unorb(1,ii).gt.0) then
C         call getmem(ncf/2+1,jrem)
      call checLB(idd,lbasi,ncf,epsx,ii,
     1            iprint,ncda,ttda)
C                 bl(jrem))
C         call retmem(1)
      endif
      end
C=========lbndia===========================================
      subroutine lbndia(ncf,nval,nstron,iprint,uniq,
     1                      ifunpair,nsym,lbasij,idmij,epsx,
     2                      equi,inndp,ncda,ttda)
C
C      remove redundant functions
C      for non-diagonal pairs
C
C      Svein Saebo, Fayetteville, AR Summer 2000
C

      use memory

        implicit real*8(a-h,o-z)
      dimension ifunpair(7,*)
c     common/big/bl(30000)
      logical*1 nstron(nval,nval)
      integer*2 lbasij(ncf,*),idmij(*),inndp(*)
      integer uniq(*),equi(2,*)
C
C      check local basis for non-diagonal pairs for redundencies
C
      lpnd=0
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      if(ii.eq.jj) cycle
      lpnd=lpnd+1
C      Remove redundant basis functions
      iuniq=uniq(lpar)
      if(iuniq.gt.0) then
      idd=idmij(lpnd)
C      call getmem(ncf/2+1,jrem)
      call checLB(idd,lbasij(1,lpnd),ncf,epsx,-lpnd,
     1            iprint,ncda,ttda)
C    2            bl(jrem))
C      call retmem(1)
      idmij(lpnd)=idd
      endif
      end do
      enddo
C      call retmem(1)
C      final check of local domains if symmetry is used
      if(nsym.gt.0)
     1  call symche2(equi,ifunpair,lbasij,idmij,nval,
     2               ncf,iprint,nstron,inndp)
      end
C=========mklbmap========================================
      subroutine mklbmap(lbas,idm,lbmap,ncf,nval,
     1                   lbasij,idmij,nstron)
C
C      constructs  a template of local domains
C      and the union if the domains of ij
C
C      Svein Saebo Fayetteville, AR summer 2000
C
      implicit integer (a-z)
      integer*2 lbas(ncf,*),idm(*),lbmap(ncf,*),null,isum
      integer*2 lbasij(ncf,*),idmij(*),en
      logical*1 nstron(nval,*)
      parameter(null=0,en=1)
C
      do ii=1,nval
      do iao=1,ncf
      lbmap(iao,ii) = null
      enddo
      enddo
      do ii=1,nval
      idd=idm(ii)
      do ild=1,idd
      iao=lbas(ild,ii)
      lbmap(iao,ii)=  en
      enddo
      enddo
C      now construct lbasij for strong nondiagonal pairs
      ij=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      if(ii.eq.jj) cycle
      ij=ij+1
      idij=0
      do ik=1,ncf
      isum=lbmap(ik,ii)+lbmap(ik,jj)
      if(isum.gt.null) then
      idij=idij+1
      lbasij(idij,ij)=ik
      endif
      enddo
      idmij(ij)=idij
      enddo
      enddo
      end
C=========checlb================================================
      subroutine ChecLB(idd,lbas,ncf,epsi,ii,
     1                      iprint,ncda,ttda)
C      this subroutine checks a local domain for redundant functions
C      idd  dimension of local domain
C      lbas local domain
C      ncf number of basis functions
C      epsi  threshold for removing redundant functions
C
C      Svein Saebo Fayetteville, AR June 1999
C      

      use memory

      implicit real*8(a-h,o-z)
      integer*2 lbas(*)
c     common/big/bl(30000)
C
      Parameter(zero=0.0d0)
      ilop=0
C
      iovl=mataddr('ovl')
  100 continue
       call matdef('Oiis','s',idd,idd)
      call matdef('ttms','d',idd*8,idd*8)
      indim=(idd*5)/2+1
      call matdef('iwork','d',indim,indim)
       ifad=idd/2+1
       call matdef('ifail','d',ifad,ifad)
       call matdef('Oii','q',idd,idd)
       isiia=mataddr('Oii')
       call compab(bl(iovl),bl(isiia),ncf,idd,idd,lbas(1),lbas(1))
       call matdef('orbi','q',idd,idd)
       iorbi=mataddr('orbi')-1
       call matdef('eigi','d',idd,idd)
       ieigi=mataddr('eigi')
       call matcopy('Oii','Oiis')
       iossa=mataddr('Oiis')
      ittms=mataddr('ttms')
        iiwork=mataddr('iwork')
       ifai=mataddr('ifail')
      abstol=2*Dlamch('S')
c      call sdiag2(idd,idd,bl(isiia),bl(ieigi),bl(iorbi+1))
      call secund(tttx)
      ncda=ncda+1
       call dspeVX('V','I','U',idd,bl(iossa),zero,zero,1,1,
     1 abstol,mrem,bl(ieigi),bl(iorbi+1),idd,bl(ittms),
     2 bl(iiwork),bl(ifai),info)
      call secund(ttty)
       ttda=ttda+ttty-tttx
       if(info.ne.0) call nerror(1,'dspeVX','diagonalization failed',
     1 info,ii)
      if(bl(ieigi).lt.epsi)then
      vmax=zero
      do lts=1,idd
      val=abs(bl(iorbi+lts))
        if(val.lt.vmax) cycle
      imx=lts
      vmax=val
      end do
      ilop=ilop+1
      if(imx.eq.idd) then
      call matrem('eigi')
      call matrem('orbi')
      call matrem('Oii')
       call matrem('ifail')
      call matrem('iwork')
      call matrem('ttms')
       call matrem('Oiis')
C      go back and check if you need to remove another function
      idd=idd-1
      goto 100
      endif
        do irm=imx,idd-1
      lbas(irm)=lbas(irm+1)
      end do
      call matrem('eigi')
      call matrem('orbi')
      call matrem('Oii')
       call matrem('ifail')
      call matrem('iwork')
      call matrem('ttms')
       call matrem('Oiis')
      idd=idd-1
C      go back and check if you need to remove another function
      goto 100
      end if
      call matrem('eigi')
      call matrem('orbi')
      call matrem('Oii')
       call matrem('ifail')
      call matrem('iwork')
      call matrem('ttms')
       call matrem('Oiis')
      if(iprint.ge.3) then
      write(6,*) ilop,' functions removed from domain for ',ii
      write(6,*) 'Dimension: ',idd
      write(6,*) 'AOs: ', (lbas(iwr),iwr=1,idd)
      endif
      end
C==========mkjfun3======================================
      subroutine mkjfun3(nval,jcmax,nstron,jfun3)
      implicit integer(a-z)
      logical*1 nstron(nval,*)
      dimension jfun3(jcmax,*)
      do ii=1,nval
      jdim=1
      do jj=1,nval
      if(nstron(ii,jj)) cycle
      jdim=jdim+1
      jfun3(jdim,ii)=jj
      end do
      jfun3(1,ii)=jdim
      end do
      end
C=========bascop=======================================:1
      subroutine bascop(lbasij,idmij,lbasijt,idmijt,nud,ncf)
      implicit integer(a-z)
      integer*2 lbasij(ncf,*),idmij(*),lbasijt(ncf,*),idmijt(*)
      do iii=1,nud
      idd=idmij(iii)
      idmijt(iii)=idd
      do imv=1,idd
      lbasijt(imv,iii)=lbasij(imv,iii)
      enddo
      enddo
      end
      subroutine mklw2(lisw2,listp,npairs)
      implicit integer(a-h)
      dimension listp(*)
      integer*2 lisw2(*)
      nwe=0
      do lpar=1,npairs
      ij=listp(lpar)
      if(ij.gt.0) cycle
      nwe=nwe+1
      lisw2(nwe)=lpar
      write(6,*)nwe,lisw2(nwe)
      enddo      
      end
C=========wriarr========================================
      subroutine wriarr(A,len,iunit)
      implicit real*8(a-h,o-z)
      dimension A(len)
      write(unit=iunit) A
      end
      subroutine reaarr(A,len,iunit)
      implicit real*8(a-h,o-z)
      dimension A(len)
      read(unit=iunit) A
      end
      subroutine move(a,b,n)
      implicit real*8(a-h,o-z)
      dimension a(*),b(*)
      do i=1,n
      b(i)=a(i)
      enddo
      end
C=========checlbn===============================================
      subroutine ChecLBn(idd,lbas,ncf,epsi,ii,
     1                      iprint,ncda,ttda,jrem)
C      this subroutine checks a local domain for redundant functions
C      idd  dimension of local domain
C      lbas local domain
C      ncf number of basis functions
C      epsi  threshold for removing redundant functions
C
C      Svein Saebo Fayetteville, AR June 1999
C      

      use memory

      implicit real*8(a-h,o-z)
      integer*2 lbas(*),jrem(*)
c     common/big/bl(30000)
C
      Parameter(zero=0.0d0,vl=-1.0d-06,vu=1.0d-06)
C
      call matdef('xqq','v',idd,idd)
      call matdef('ovlu','q',idd,idd)
      iovlu=mataddr('ovlu')
      call matdef('ovluf','q',ncf,ncf)
      iovluf=mataddr('ovluf')
      call matdef('ovls','s',ncf,ncf)
      call matread('ovls',1,'s matrix')
      call matcopy('ovls','ovluf')
      call matrem('ovls')
      call compab(bl(iovluf),bl(iovlu),ncf,idd,
     1                idd,lbas,lbas)
      call matrem('ovluf')
C
      iovl=mataddr('ovl')
       call matdef('orbi','q',idd,idd)
       iorbi=mataddr('orbi')
       call matdef('eigi','d',idd,idd)
       ieigi=mataddr('eigi')
       call matdef('Oii','q',idd,idd)
       isiia=mataddr('Oii')
      call matdef('ttms','d',idd*8,idd*8)
      indim=(idd*5)/2+1
      call matdef('iwork','d',indim,indim)
       ifad=idd/2+1
       call matdef('ifail','d',ifad,ifad)
       call compab(bl(iovl),bl(isiia),ncf,idd,idd,lbas(1),lbas(1))
       call matdef('Oiis','s',idd,idd)
       call matcopy('Oii','Oiis')
       iossa=mataddr('Oiis')
      ittms=mataddr('ttms')
        iiwork=mataddr('iwork')
       ifai=mataddr('ifail')
      abstol=2*Dlamch('S')
C      call sdiag2(idd,idd,bl(isiia),bl(ieigi),bl(iorbi))
       call dspeVX('V','V','U',idd,bl(iossa),vl,vu,1,1,
     1 abstol,mrem,bl(ieigi),bl(iorbi),idd,bl(ittms),
     2 bl(iiwork),bl(ifai),info)
       if(info.ne.0) call nerror(1,'dspeVX','diagonalization failed',
     1 info,ii)
      call matrem('Oiis')
      call matrem('ifail')
      call matrem('iwork')
      call matrem('ttms')
C      mrem=0
C      do ix=1,idd
C      if(bl(ieigi-1+ix).lt.1.0d-06) then
C      mrem=mrem+1
C      else
C      exit
C      endif
C      enddo
       write(6,*) mrem,' functions must be removed'
      if(mrem.eq.0) then
      call matrem('Oii')
      call matrem('eigi')
      call matrem('orbi')
      call matrem('ovlu')
      call matrem('xqq')
      return
      endif
      call matscal('eigi',1.0d08)
      mkeep=idd-mrem
      ixqqa=mataddr('xqq')
       call matsub('u','orbi',1,mrem)
          call matdef('ut','r',mrem,idd)
       call matpose2('u','ut','n')
       call matdef('uut','q',idd,idd)
      iuuts=mataddr('uut')
      iust=mataddr('u')
       call matmmult('u','ut','uut')
C
      idia=-idd-1
      do ix=1,idd
      idia=idia+idd+1
       bl(ixqqa-1+ix) = bl(iuuts+idia)/bl(isiia+idia)
      enddo
      call matprint('xqq',6)
      call matrem('uut')
      call matrem('ut')
       call matrem('u')
      call matrem('Oii')
      call matrem('eigi')
      call matrem('orbi')
      call matrem('ovlu')
      do irem=1,mrem
      krem=idamax(idd,bl(ixqqa),1)
      bl(ixqqa-1+krem)=zero
      jrem(irem)=krem
      enddo
C      remove the mrem functions from lbas
      call lbsort(jrem,mrem)
C      write(6,*)'jrem',(jrem(iw),iw=1,mrem)
      do krem=mrem,1,-1
      jst=jrem(krem)
      if(jst.eq.idd) then
      idd=idd-1
      else
      idd=idd-1
      do irem=jst,idd
      lbas(irem)=lbas(irem+1)
      end do
      endif
      enddo   !over krem
C
      write(6,*) mrem,' functions removed from domain for ',ii
      write(6,*) 'Dimension: ',idd
      write(6,*) 'AOs: ', (lbas(iwr),iwr=1,idd)
      call matrem('xqq')
      end
C=========checlb2===============================================
      subroutine ChecLB2(idd,lbas,ncf,epsi,ii,
     1                      iprint,ncda,ttda,jrem)
C      this subroutine checks a local domain for redundant functions
C      idd  dimension of local domain
C      lbas local domain
C      ncf number of basis functions
C      epsi  threshold for removing redundant functions
C
C      Svein Saebo Fayetteville, AR June 1999
C      

      use memory

      implicit real*8(a-h,o-z)
      integer*2 lbas(*),jrem(*)
c     common/big/bl(30000)
C
      Parameter(zero=0.0d0,vl=-1.0d-06,vu=1.0d-06)
C
      call matdef('ovlu','q',idd,idd)
      iovlu=mataddr('ovlu')
      call matdef('ovluf','q',ncf,ncf)
      iovluf=mataddr('ovluf')
      call matdef('ovls','s',ncf,ncf)
      call matread('ovls',1,'s matrix')
      call matcopy('ovls','ovluf')
      call matrem('ovls')
      call compab(bl(iovluf),bl(iovlu),ncf,idd,
     1                idd,lbas,lbas)
      call matrem('ovluf')
C
      iovl=mataddr('ovl')
       call matdef('orbi','q',idd,idd)
       iorbi=mataddr('orbi')
       call matdef('eigi','d',idd,idd)
       ieigi=mataddr('eigi')
       call matdef('Oii','q',idd,idd)
       isiia=mataddr('Oii')
      call matdef('ttms','d',idd*8,idd*8)
      indim=(idd*5)/2+1
      call matdef('iwork','d',indim,indim)
       ifad=idd/2+1
       call matdef('ifail','d',ifad,ifad)
       call compab(bl(iovl),bl(isiia),ncf,idd,idd,lbas(1),lbas(1))
       call matdef('Oiis','s',idd,idd)
       call matcopy('Oii','Oiis')
       iossa=mataddr('Oiis')
      ittms=mataddr('ttms')
        iiwork=mataddr('iwork')
       ifai=mataddr('ifail')
      abstol=2*Dlamch('S')
C      call sdiag2(idd,idd,bl(isiia),bl(ieigi),bl(iorbi))
       call dspeVX('V','V','U',idd,bl(iossa),vl,vu,1,1,
     1 abstol,mrem,bl(ieigi),bl(iorbi),idd,bl(ittms),
     2 bl(iiwork),bl(ifai),info)
       if(info.ne.0) call nerror(1,'dspeVX','diagonalization failed',
     1 info,ii)
      call matrem('Oiis')
      call matrem('ifail')
      call matrem('iwork')
      call matrem('ttms')
C      mrem=0
C      do ix=1,idd
C      if(bl(ieigi-1+ix).lt.1.0d-06) then
C      mrem=mrem+1
C      else
C      exit
C      endif
C      enddo
       write(6,*) mrem,' functions must be removed'
      if(mrem.eq.0) then
      call matrem('Oii')
      call matrem('eigi')
      call matrem('orbi')
      call matrem('ovlu')
      return
      endif
      call getmem(idd,ivec) ! MMOK
      call findf(bl(iorbi),jrem,bl(ivec),idd,mrem)
      call retmem(1)
C
      call matrem('Oii')
      call matrem('eigi')
      call matrem('orbi')
      call matrem('ovlu')
C      remove the mrem functions from lbas
      call lbsort(jrem,mrem)
       write(6,*)'jrem',(jrem(iw),iw=1,mrem)
      do krem=mrem,1,-1
      jst=jrem(krem)
      if(jst.eq.idd) then
      idd=idd-1
      else
      idd=idd-1
      do irem=jst,idd
      lbas(irem)=lbas(irem+1)
      end do
      endif
      enddo   !over krem
C
      write(6,*) mrem,' functions removed from domain for ',ii
      write(6,*) 'Dimension: ',idd
      write(6,*) 'AOs: ', (lbas(iwr),iwr=1,idd)
      end
      subroutine findf(u,jrem,vec,idd,mrem)
      implicit real*8(a-h,o-z)
      dimension U(idd,mrem),vec(*)
      integer*2 jrem(*)
      do iq=1,idd
      sum=0.0
      do irem=1,mrem
      sum=sum+U(iq,irem)**2
      enddo
      vec(iq)=sum
      enddo
      write(6,*)'vec',(vec(iw),iw=1,idd)
      do irem=1,mrem
      krem=idamax(idd,vec,1)
      jrem(irem)=krem
      vec(krem)=0.0d0
      enddo
      write(6,*)'jrem',(jrem(iw),iw=1,mrem)
      end
      subroutine ortho(col1,col2,idd)
      implicit real*8(a-h,o-z)
      dimension col1(*),col2(*)
      sum=0.0d0
      do ii=1,idd
      sum=sum+col1(ii)*col1(ii)
      enddo
      do ii=1,idd
      col1(ii)=col1(ii)/sqrt(sum)
      enddo
      sub=0.0d0
      do ii=1,idd
      sub=sub+col2(ii)*col1(ii)
      enddo
      do ii=1,idd
      col2(ii) =sub*col1(ii)-col2(ii)
      enddo
      prod=0.0d0
      do ii=1,idd
      prod=prod+col1(ii)*col2(ii)
      enddo
      write(6,*) 'prod',prod
      end

      subroutine nmr_tensor(bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
c     common /intbl/ifpp,inxx(100)
      common /aux/ nreq,nucnr(2000)
c
      common /cux/ lden,lfoc,lforc,love,lh01, lhfc,lval,lvec,lrem
c
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      common /nmrint/ ichf,ncache,nprint,nogiao,noacce,ngauge,ntrans
      common /gauge/ gauger(3)
      dimension bl(*)
      dimension inx(12,*)
      character cdum*20
      character*256 jobname
      Logical Symflag
c
      parameter (IUnit=1)
c----------------------------------------------------------------------
c memory allocations :
c
      call memres2(bl,ncf,na,lden,lh01,lforc,ldn1)
c
      lhna=lh01
c----------------------------------------------------------------------
c update value of last  :
      call getmem(1,last)
      call retmem(1)
c----------------------------------------------------------------------
c
c  get jobname header
c
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c----------------------------------------------------------------------
c set up the gauge-origin for both the CPHF and GIAO
c
      call rgauge(bl(inuc),ngauge,gauger)
c----------------------------------------------------------------------
c check-run only or lack of scf convergance :
c
      if (lflag(2).eq.1.or.lflag(3).eq.1) return
c----------------------------------------------------------------------
      ntri=ncf*(ncf+1)/2
c----------------------------------------------------------------------
c select symmetry unique atoms for nmr shielding
c and make the atom/operation equivalence list
c
      open (unit=IUnit,file=jobname(1:lenJ)//'.sym',
     $      form='formatted',status='old')
      call rdcntrl(IUnit,7,'$ntrans',1,ntrans,dum,cdum)
      close (unit=IUnit,status='keep')
c
c are there any type 2 dummy atoms?
c (these are additional centres to calculate chemical shift)
c
      open (unit=IUnit,file=jobname(1:lenJ)//'.control',
     $      form='formatted',status='old')
      call fdcntrl(IUnit,7,'$ndummy',idum)
      If(idum.EQ.0) Then
       backspace IUnit
       READ(IUnit,900) ndum1,ndum2
  900 Format(9X,I4,2X,I4)
      Else
       ndum1 = 0
       ndum2 = 0
      EndIf
      close (unit=IUnit,status='keep')
c----------------------------------------------------------------------
c allocate memory for rotation matrix and symmetry-unique atoms list
c
      call mmark
c
      call getmem(9,irm)
      call getmem(9,ip)
      call getint(na,iunq)
      call getint(3*na,ilst)
c
c----------------------------------------------------------------------
c check for symmetry usage
c
      If(ntrans.eq.1) Then
       Symflag = .False.
       itrn = 1
       ineq = 1
       itrnA = 1
      Else
c
c allocate memory for symmetry operations etc..
c
       Symflag = .True.
       call getmem(9*ntrans,itrn)
       call getmem(na*ntrans,ineq)
      EndIf
c
c now read in FULL symmetry information
c
      natom = na-ndum1-ndum2
      Call RdSYM(Symflag,natom,  bl(irm), cdum,   ntrans,
     $           ndeg,   nq,  bl(iunq),bl(itrn),bl(ineq))
c
c do not calculate shift for charged dummies
c
      If(ndum1.gt.0) Call ModNucnr(natom,ndum1,nreq,bl(ilst),nucnr)
c
c make atom/operation equivalence list
c
      Call MakeList(na,     ndum2,  natom,  nq,   bl(iunq),
     $              ntrans,bl(ineq),bl(ilst))
c----------------------------------------------------------------------
c print out symmetry unique atoms
c
      if(ntrans.gt.1) call printsymat(bl(iunq),nq)
c----------------------------------------------------------------------
      call para_shift_start(bl,lden,ldn1,ntri)
c----------------------------------------------------------------------
c
      call secund(t2b)
      call elapsec(et2b)
c
c loop over symmetry-unique nuclei
c
      do 100 nra=1,nq
         call int_from_bl( bl(iunq), nra , nrat )
c        nrat=inxx(iunq+nra)
c
c  check if this atom has to be calculated
c
         do ii=1,nreq
           iam=nucnr(ii)
           if(iam.eq.nrat) go to 90
         end do
         go to 100
  90     continue
c
c
         call para_shi(na,bl,inx,last,nprint,iout,bl(ibas),bl(inuc),
     *                 ncs,ncf,ntri,nreq,bl(lden),bl(ldn1),bl(lforc),
     *                 bl(lhna),nra,nrat,ldn1,lhna)
c
  100 continue
c
      call para_shift_end(bl(lforc),nreq)
      call para_next(0)
c----------------------------------------------------------------------
      call secund(t2e)
      call elapsec(et2e)
      total=(t2e-t2b)/60.0d0
      elaps=(et2e-et2b)/60.0d0
      write(iout,400) total,elaps
c----------------------------------------------------------------------
c
      call CalcShift(na,  bl(lforc),  bl,     bl(irm),bl(ip),
     $               ntrans,bl(itrn),bl(ilst),nreq, nucnr,
     2               iout,   icond, nprint)
cc      call symmsh(bl,iout,icond,nreq,nucnr,ntri,lforc)
c----------------------------------------------------------------------
c NICS points
c
c allocate memory for NICS calculations
      call getival('nics',nics)
      if(nics.gt.0) then
         call getmem(3*nics  ,inicx)
         call getmem(9*nics*2,inics)
         call zeroit(bl(inics),nics*18)
c
c get atoms that define NICS plane
         call getival('iatom1',iat1)
         call getival('iatom2',iat2)
         call getival('iatom3',iat3)
         call getival('ngrid',ngrid)
         call getrval('disp',disp)
         call getrval('step',step)
c
       call find_nicsgrid(iat1,iat2,iat3,bl(inuc),na,disp,ngrid,
     *                    step,bl(inicx),nics,igeoc)
c----------------------------------------------------------------------
c
         do ipo=1,nics
            call shnics(na,bl,inx,last,nprint,iout,bl(ibas),bl(inicx),
     *                  ncs,ncf,ntri,nics,bl(lden),bl(ldn1),bl(inics),
     *                  bl(lhna),ipo)
         enddo
c
         call calcsnics(nics,igeoc,bl(inics),bl(inicx),bl,iout,icond,
     *                  nprint)
c
      endif
c
c----------------------------------------------------------------------
c return all memory allocated in this routine
      call retmark
c----------------------------------------------------------------------
  400 format('Master CPU time for GIAO 1e integrals = '
     *,f8.2,'  Elapsed = ',f8.2,' min'/)
c
      end
c===================================================================
      subroutine memres2(bl,ncf,na,lden,lh01,lforc,ldn1)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /lindvec/ lind,idensp
      dimension bl(*)
c-----------------------------------------------------------------------
c LFORC is the starting address for the nuclear shielding tensor
c
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
      na18=na*18
c-----------------------------------------------------------------------
      call getival('lden',lden) ! from chshift
c
      call getmem(na18 ,lforc)
      call getmem(ntri3, lh01) ! only for one nucleus at the time
      call getmem(ntri3,ldn1 )
c
c-----------------------------------------------------------------------
      call getmem(ncf  ,lind )
      call setup_lind(bl(lind),ncf)
c------------------------------------------------
c zero out : Forc, H01
c
      call zeroit(bl(lh01),ntri3)
      call zeroit(bl(lforc),na18)
c------------------------------------------------
c read-in the density matrix :
c
      np4=4
      call rea(bl(ldn1),ntri3,np4,'den1_rhf')
c------------------------------------------------
c update value of last  :
c     call getmem(1,last)
c     call retmem(1)
c-----------------------------------------------------------------------
      end
c=======================================================================
      subroutine shield(nna,bl,inx,last,iprint,iout,datbas,datnuc,ncs,
     *                ncf,ntri,nreq,den,dn1,shie,hna,nra,nrat,ldn1,lhna)
c-----------------------------------------------------------------------
c        this subroutine calculates 1-el integral derivatives which
c                  are connected with 1/r1n operator.
c
c             there are the following integral derivatives:
c
c 1.   (i00/h11n/j00) + (i10/h01n/j00) + (i00/h01n/j10)     matrix
c
c 2.   (i00 / h01n/j00)                                     matrix
c
c These calculations are performed for ONE nucleus (NRAT) at the time
c-----------------------------------------------------------------------
c nna 	- number of atoms
c bl	- real storage area
c inx	- integer storage
c last	- pointer to top of free mem.
c iprint- print level
c iout	- unit of stdout
c datbas- basis set data struct.
c datnuc- nuclear data struct.
c ncs,ncf - no of contr shells and func
c nreq	- number of atoms shielding need to be calc. for
c den	-  D00
c dn1	- pointer to D10
c shie	- shieldings (9*nreq diam.+9*nreq para
c hna	- H01
c nra 	- atom worked on (position in list of symm equ. atoms)
c nrat	- atom worked on (number)
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /lindvec/ lind,idensp
      dimension bl(1)
      dimension den(ntri),dn1(ntri,3),hna(ntri,3),shie(3,3,*)
      dimension datbas(13,*),datnuc(5,*)
      dimension inx(12,1)
c needed 12*I*I : 12*28*28
      dimension s2(9408),xra(3)
c
c units for shielding tensor :
c  1 au = 26.6256691d+00 ppm
c
      data f12c2p / 26.6256691d+00 /
c
c remember value of last
      lenter=last
c
c zero out h01,set up pointers to dia and paramagn.
c
      call zeroit(hna,ntri*3)
c
c nrat is nucleus of interest (nrat>0) or
c dummy atom at (0,0,0)
c
      if(nrat.gt.0) then
         xra(1)=datnuc(2,nrat)
         xra(2)=datnuc(3,nrat)
         xra(3)=datnuc(4,nrat)
      else
         xra(1)=zero
         xra(2)=zero
         xra(3)=zero
      endif
c
      kk=0
      k2=0
c
      key=4
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
            len=len1*len2
            ls2=12*len
            ld=12
c
        do 45 jgc=0,ngcjx
c
      call onesh(i,j,igc,jgc,datbas,bl,s2,inx,key,nna,kk,k2,xra ,ls2,ld)
c
            iijxx=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
c
               iijxx=iijxx+1
               iijxy=iijxx+len
               iijxz=iijxy+len
c
               iijyx=iijxz+len
               iijyy=iijyx+len
               iijyz=iijyy+len
c
               iijzx=iijyz+len
               iijzy=iijzx+len
               iijzz=iijzy+len
         if(jff.gt.iff) go to 40
               ij=ii+jff
               d00=den(ij)
               if(iff.ne.jff) d00=two*d00
c
            hna(ij,1)=hna(ij,1) + s2(len*9+iijxx)
            hna(ij,2)=hna(ij,2) + s2(len*9+iijxy)
            hna(ij,3)=hna(ij,3) + s2(len*9+iijxz)
c
            shie(1,1,nra)=shie(1,1,nra) + s2(iijxx)*d00
            shie(2,1,nra)=shie(2,1,nra) + s2(iijxy)*d00
            shie(3,1,nra)=shie(3,1,nra) + s2(iijxz)*d00

            shie(1,2,nra)=shie(1,2,nra) + s2(iijyx)*d00
            shie(2,2,nra)=shie(2,2,nra) + s2(iijyy)*d00
            shie(3,2,nra)=shie(3,2,nra) + s2(iijyz)*d00
c
            shie(1,3,nra)=shie(1,3,nra) + s2(iijzx)*d00
            shie(2,3,nra)=shie(2,3,nra) + s2(iijzy)*d00
            shie(3,3,nra)=shie(3,3,nra) + s2(iijzz)*d00
c
   40       continue
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c
c return value of last
      last=lenter
c
      if(iprint.gt.2) then
        write(iout,444)
        write(iout,445) nrat
        call druma(hna(1,1),ncf,iout,'h01mat x ')
        call druma(hna(1,2),ncf,iout,'h01mat y')
        call druma(hna(1,3),ncf,iout,'h01mat z')
      endif
c
 444  format('*******************************')
 445  format(/'***  total matrix h01 for atom number  ',i2,'     ***')
c
c Calculate the paramagnetic part of the shielding tensor
c
c the trace of the product of two antisymmetric matrices is really
c simple  Tr(AB)= - Aik Bik   (with the Einstein convention)
c in the present case we also have a factor of -1
c as both matrices are imaginary
c if both square matrices are stored in the triangular form
c the whole thing is a simple blas dot product (the diagonal is zero!)
      do 100 l=1,3
      do 100 m=1,3
         shie(m,l,nra+nreq)=ddot(ntri,dn1(1,l),1,hna(1,m),1)*two
 100  continue
      If(iprint.ge.2) call anal_pshie(shie,nra,nreq,ncf,ldn1,lhna)
c
c Convert diamagnetic and paramagnetic shielding
c from au to ppm
c
      call dscal(9,f12c2p,shie(1,1,nra),1)
      call dscal(9,f12c2p,shie(1,1,nra+nreq),1)
      end
c=======================================================================
      subroutine printsh(bl, nrat, shw,shn, iout,shiel)
      implicit real*8 (a-h,o-z)
      character*5 shiel
      dimension bl(*)
      dimension shw(3,3),shn(3,3)
      data three /3.d0/
c
      if (shiel.eq.'total') then
         aivst=( shn(1,1)+shn(2,2)+shn(3,3) )/three
c
           write(iout,100)
           write(iout,250) nrat
           write(iout,100)
           write(iout,350) shn(1,1),shn(2,1),shn(3,1)
           write(iout,350) shn(1,2),shn(2,2),shn(3,2)
           write(iout,350) shn(1,3),shn(2,3),shn(3,3)
           write(iout,100)
           write(iout,450) aivst
           write(iout,100)
c
      else if(shiel.eq.'short') then
           aivst=( shn(1,1)+shn(2,2)+shn(3,3) )/three
           write(iout,'( "Isotropic shielding of atom ",i4,
     1       " = ", f10.4)') nrat,aivst
      else
c
         aivsp=( shw(1,1)+shw(2,2)+shw(3,3) )/three
         aivsd=( shn(1,1)+shn(2,2)+shn(3,3) )/three
c
           write(iout,200) nrat
           write(iout,100)
           write(iout,300) shn(1,1),shn(2,1),shn(3,1),
     *                     shw(1,1),shw(2,1),shw(3,1)
           write(iout,300) shn(1,2),shn(2,2),shn(3,2),
     *                     shw(1,2),shw(2,2),shw(3,2)
           write(iout,300) shn(1,3),shn(2,3),shn(3,3),
     *                     shw(1,3),shw(2,3),shw(3,3)
           write(iout,100)
           write(iout,400) aivsd,aivsp
c
      endif
c------------------------------------
  100 format(72(1h-))
  200 format('--',10X,12HDiamagnetic ,7X,5HATOM ,I2,13X,12HParamagnetic,
     *6X,' --')
  300 format(F10.4,2(2X,F10.4),1X,1H*,1X,F10.4,2(2X,F10.4))
  400 format('- AVERAGE =',F12.5,25X,F12.5, 9X,'  -' )
c------------------------------------
  250 format('--',18X,' Total shielding for atom no ',i2,18X,' --')
  350 format(18x,F10.4,2(2X,F10.4) )
  450 format('-',' AVERAGE =',18x,F12.5,28X,'  -')
c------------------------------------
      end
c=================================================================
      subroutine druma3(bl,lhfc,ncf,ntri,iout,text)
      implicit real*8 (a-h,o-z)
      character*1 text
      dimension bl(*)
c
         write(iout,444)
         if(text.eq.'1' .or. text.eq.'2') then
           write(iout,5550)
           write(iout,5551)
           if(text.eq.'1') write(iout,5552)
           if(text.eq.'2') write(iout,5553)
         endif
         if(text.eq.'3') write(iout,5554)
c
      call druma(bl(lhfc),ncf,iout,'DENS1-X ')
      call druma(bl(lhfc+ntri),ncf,iout,'DENS1-Y ')
      call druma(bl(lhfc+ntri*2),ncf,iout,'DENS1-Z ')
c
         write(iout,444)
c
  444 format(72(1H*))
 5550 format(/'***  CONSTANT PART OF I-ORDER DENSITY MATRIX  *** ')
 5551 format ('***            contributions from             *** ')
 5552 format ('***         D10[ F10(D0,G1) - e*S10 ]         *** '/)
 5553 format ('***  D10[ F10(D0,G1)-e*S10 ] + 0.5*D0*S10*D0  *** '/)
 5554 format(/'***    THE FINAL FIRST-ORDER DENSITY MARIX    *** '/)
c
      end
c====================================================================
      SUBROUTINE CalcShift(NAtoms, shie,  bl,     shs,    p,
     $                     NTrans, TRANS, ILST,   nreq,   nucnr,
     2                     IOut,  ICon,   iprint)
      IMPLICIT REAL*8(A-H,O-Z)
c--------------------------------------------------------------------
C  Calculates shielding in principal coordinate system for atoms
c  requested
C  ** WARNING **
C     direction cosines may be have degenerate components in
C     symmetrical molecules; in this case the degenerate
c     components are ill-defined
C
C  ARGUMENTS
C
C  NAtoms  -  number of atoms
C  lforc   -  pointer into bl array for shielding tensor
C  bl      -  main storage array
C  shs     -  3x3 matrix for shielding tensor
C  p       -  3x3 work array
C  NTrans  -  number of symmetry operations
C  TRANS   -  symmetry operations as 3x3 transformation matrices
C  ILST    -  atom/operation equivalence list
C              (I,1) - which symmetry-unique atom will generate atom I
C              (I,2) - which symmetry operation is involved
C                      (0 means atom is unique)
C              (I,3) - position in list of symmetry-unique atoms
C              nreq - requested atoms (if the FOR option is given)
C              nucnr(1:nreq): atoms to be calculated
C  IOut    -  unit for standard output
C  ICon    -  unit for summary output
C  iprint  -  print level
C
C   Note: the FOR option may be incompatible with symmetry
C   It is best to use only one of them
c--------------------------------------------------------------------
c -- automatic allocation of an array
      CHARACTER*8 AtSymb(NAtoms)
c--------------------------------------------------------------------
      dimension bl(*),shie(3,3,*)
      dimension TRANS(3,3,NTrans),ILST(NAtoms,3)
      dimension p(3,3),prs(3),shs(3,3)
      dimension nucnr(*)
c
      Parameter (Zero=0.0d0,half=0.5d0,three=3.0d0,one=1.0d0)
c
c--------------------------------------------------------------------
c get atomic symbols
c
      call getatsymb(bl,atsymb,natoms)
c--------------------------------------------------------------------
c -- initial print headers
      write(iout,444)
      write(iout,611)
      write(iout,612)
      write(iout,613)
      write(iout,444)
      write(icon,444)
      write(icon,'( 30x,"NMR SHIELDINGS")' )
      call tstrval('xlvsh',iyes)
      if(iyes.eq.1) then
        xlvsh=rgetrval('xlvsh')
        if(abs(xlvsh).gt.1.0D-8) then
          write(iout,*) 'Level shift=',xlvsh
          write(icon,*) 'Level shift=',xlvsh
        end if
      end if
      write(icon,444)
c
c **************************************************************
c  WARNING:  Following code is not ideal
c  Loops over all atoms and (unnecessarily) rediagonalizes
c  shielding tensor. Can be rearranged to loop over unique
c  atoms only, but either need extra storage to retain all
c  quantities prior to printing OR print out in jumbled
c  atom order
c **************************************************************
c
      do 600 IAtm=1,NAtoms
c  check if this atom has to be calculated
        do ii=1,nreq
          iam=nucnr(ii)
          if(iam.eq.iatm) go to 90
        end do
        go to 600
  90    continue
c
         nra = ILST(IAtm,3)       ! position in list of symmetry-unique atoms
         nshw = nra + nreq
c
         do 100 j=1,3
         do 100 i=1,3
c -- symmetric part of shielding tensor
           shs(i,j)=(shie(i,j,nra)+shie(j,i,nra)+shie(i,j,nshw)+
     *              shie(j,i,nshw))*half
 100    continue
c
c -- diagonalize shielding tensor
         nd=3
         call sdiag2(nd,nd,shs,prs,shs)
c
         averag=(prs(1)+prs(2)+prs(3))/three
         anisot= prs(3)-(prs(2)+prs(1))*half
c
c -- if this atom is unique, nothing more to do
c -- otherwise need to transform shielding tensor according to
c -- symmetry operation
c
         If(ILST(IAtm,2).NE.0) Then
           iop = ILST(IAtm,2)
c
           do 5 i=1,3
           do 5 j=1,3
           p(i,j) = Zero
           do 5 k=1,3
           p(i,j) = p(i,j) + TRANS(i,k,iop)*shs(k,j)
 5         continue
c
           do 6 i=1,3
           do 6 j=1,3
           shs(i,j) = Zero
           do 6 k=1,3
           shs(i,j) = shs(i,j) + p(i,k)*TRANS(j,k,iop)
 6         continue
         EndIf
c
c -- final printout
c
c Print out diamagnetic, paramagnetic and total
c shielding tensor (only for symm. unique)
c
         If(ILST(IAtm,2).EQ.0) Then
            if(iprint.ge.1) then
               call printsh(bl,IAtm,shie(1,1,nshw),
     *                      shie(1,1,nra), iout,'blank')
            endif
            call daxpy(9,one,shie(1,1,nshw),1,shie(1,1,nra),1)
            call zeroit(shie(1,1,nshw),9)
            if(iprint.ge.1) then
               call printsh(bl,IAtm,shie(1,1,nshw),
     *                      shie(1,1,nra), iout,'total')
cc             call printsh(bl,IAtm,shie(1,1,nshw),
cc   *                      shie(1,1,nra), icon,'short')
            endif
         EndIf
cccc     write(iout,614) IAtm
         write(iout,614) atsymb(iatm)(1:2),IAtm
         write(iout,615)
         write(iout,610) prs(1),(shs(i,1),i=1,3)
         write(iout,610) prs(2),(shs(i,2),i=1,3)
         write(iout,610) prs(3),(shs(i,3),i=1,3)
         write(iout,616) averag,anisot
         write(iout,444)
         write(icon,'(a2,2x, "Atom=",i4," Isotropic shielding=",f12.5,
     1    " Anisotropy=",f12.5)' )atsymb(iatm)(1:2),IAtm,averag,anisot
c
 600     continue
         write(icon,444)
c***********************************************************************
 444  format(72(1H-))
 610  format( '--  ',F13.5,12X,3(F12.7,1X),'  --')
 611  format('--',16X,'     The Total Shielding Tensor',21X,'--')
 612  format('--',9X,'in principal axes and direction cosines between',
     *12X,'--')
 613  format('--',10X,' molecular and laboratory coordinate systems',14X
     1,'--')
C614  format('--   Shielding',9X,'Atom ',I4,15X,'Direction Cosines',
 614  format('--   Shielding',6X,a2,' Atom ',I4,15X,'Direction Cosines',
     *6X,'--')
 615  format(72(1H-))
 616  format(/'-- Average =',F12.5,' , Anisotropy=',F12.5,20X,'--')
c***********************************************************************
      end
c =====================================================================
c
      SUBROUTINE ModNucnr(NAtom,  Ndum,   NReq,   IVec,   Nucnr)
      IMPLICIT INTEGER(A-Z)
C
C  Modify list of which atomic centres to calculate shift by
C  specifically excluding charged dummy atoms
C
C  ARGUMENTS
C
C  NAtom   -  number of real atoms
C  Ndum    -  number of charged dummy atoms
C  NReq    -  number of centres to calculate
C             (may be modified on exit)
C  IVec    -  scratch vector
C  Nucnr   -  list of which centres to calculate
C             (may be modified on exit)
C
C
      DIMENSION IVec(NReq),Nucnr(NReq)
C
C
C  initialize
C
      NS = 0
      CALL IZeroIT(IVec,NReq)
C
C  loop over all dummy atoms
C
      DO 20 IAtm=NAtom+1,NAtom+Ndum
C
C  loop over list of atoms
C
      DO 10 I=1,NReq
      If(Nucnr(I).EQ.IAtm) Then
        NS = NS+1
        IVec(I) = 1
      EndIf
 10   CONTINUE
 20   CONTINUE
C
C  If there were any dummies in the list, remove them
C
      IF(NS.GT.0) THEN
        II = 0
        DO 30 I=1,NReq
        If(IVec(I).EQ.0) Then
          II = II+1
          Nucnr(II) = Nucnr(I)
        EndIf
 30     CONTINUE
        NReq = NReq-NS
      ENDIF
C
      RETURN
      END
c=======================================================================
      SUBROUTINE MakeList(NAtoms, Ndum,   NAtom,  NQ,     IUNQ,
     $                    NTrans, NEqATM, ILST)
      IMPLICIT INTEGER(A-Z)
C
C  Make atom/operation equivalence list
C  i.e., for each atom, which operation acting on which symmetry
C  unique atom will generate that atom
C
C  ARGUMENTS
C
C  NAtoms  -  total number of atomic centres (including dummy atoms)
C  Ndum    -  number of type 2 (i.e., no charge) dummy atoms
C  NAtom   -  number of real atoms
C  NQ      -  number of symmetry-unique real atoms
C  IUNQ    -  list of symmetry-unique atoms
C  NTrans  -  number of symmetry operations
C  NEqATM  -  list of atomic equivalences under symmetry operations
C  ILST    -  atom/operation equivalence list
C              (I,1) - which symmetry-unique atom will generate atom I
C              (I,2) - which symmetry operation is involved
C                      (0 means atom is unique)
C              (I,3) - position in list of symmetry-unique atoms
C
C
      Dimension IUNQ(NQ),NEqATM(NAtom,NTrans),ILST(NAtoms,3)
C
C
C  Clear the ILST array
C
      Call IZeroIT(ILST,3*NAtoms)
C
C  If NQ=NAtom, ILST array is trivial
C  (ALL atoms must be calculated plus any dummies)
C
      If(NQ.EQ.NAtom) Then
       DO 10 I=1,NAtoms
       ILST(I,1) = I
       ILST(I,3) = I
 10    CONTINUE
       If(Ndum.GT.0) Then
        NS = NAtoms-Ndum
        DO 15 J=1,Ndum
        JJ = NS+J
        NQ = NQ+1
        ILST(JJ,1) = NQ
        ILST(JJ,3) = NQ
        IUNQ(NQ) = JJ
 15     CONTINUE
       EndIf
       RETURN
      EndIf
C
C  Put symmetry-unique atoms into ILST
C
      DO 20 I=1,NQ
      II = IUNQ(I)
      ILST(II,1) = II
      ILST(II,3) = I
 20   CONTINUE
C
C  Loop over symmetry-unique atoms and find which non-unique
C  atoms can be generated by it
C
      ITot = NQ        ! number of atoms accounted for
c
      DO 50 I=1,NQ
      II = IUNQ(I)
C
C  loop over symmetry operations
C
      DO 30 IOP=1,NTrans
C
C  if operation generates an atom other than itself
C  we have a candidate for inclusion in ILST
C
      JJ = NEqATM(II,IOP)
      If(JJ.NE.II.AND.ILST(JJ,1).EQ.0) Then
       ITot = ITot+1
       ILST(JJ,1) = II        ! equivalent symmetry-unique atom to atom JJ
       ILST(JJ,2) = IOP       ! symmetry operation that generates JJ
       ILST(JJ,3) = I         ! position in list of symmetry-unique atoms
      EndIf
 30   CONTINUE
C
C  have we finished?
C
      If(ITot.EQ.NAtom) Then
C
C  finished - add dummy atoms if any
C
       If(Ndum.GT.0) Then
        NS = NAtoms-Ndum
        DO 40 J=1,Ndum
        JJ = NS+J
        NQ = NQ+1
        ILST(JJ,1) = NQ
        ILST(JJ,3) = NQ
        IUNQ(NQ) = JJ
 40     CONTINUE
       EndIf
c
       RETURN
      EndIf
C
 50   CONTINUE
C
C  should not get here!!
C
      Call nerror(1,'NMR Chemical Shifts',
     $              'Problems making atom/operation equivalence list',
     $               0,0)
C
      END
c =====================================================================
      subroutine getatsymb(bl,atsymb,natoms)

      use memory, block => bl

      implicit real*8(a-h,o-z)
      character*8 atsymb(natoms)
      character*256 jobname
      common /job/jobname,lenJ
      dimension bl(*)
c
      iunit=1
c
      call getmem(natoms*3,ixc)
      call getmem(natoms  ,ixq)
      call getmem(natoms  ,ima)
c
      jnk=1
C
C  Read the Cartesian coordinates
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.coord',
     $      FORM='FORMATTED',STATUS='OLD')
      call RdCoordF(IUnit,NAtoms,AtSymb,bl(iXC),-1,jnk,bl(ixq),bl(ima))
      CLOSE(UNIT=IUnit,STATUS='KEEP')
c
      call retmem(3)
c
      do iat=1,natoms
c        write(6,66) iat,atsymb(iat)(1:2)
         call uppercase(atsymb(iat),8)
c        write(6,66) iat,atsymb(iat)(1:2)
      enddo

  66  format(' atom ',i2,4x,a2)
C
      end
c======================================================================
      subroutine printsymat(iunq,nq)
      common /tape/inp,inp2,iout,ipun,iarc,icond,itest,npl(9)
      dimension IUNQ(NQ)
c
      write(iout,*) ' '
      write(iout,*)' Symmetry unique atoms :'
      write(iout,*) ' '
c
      write(iout,666) (iunq(iuat),iuat=1,nq)
c
  666 format(15(i3,1x))
c
      end
c======================================================================
      subroutine anal_pshie(shie,nra,nreq,ncf,ldn1,lhna)
c This routine partitions the paramagnetic shielding between individual
c excitations between canonical orbitals. The top contributions are
c printed out.
c The paramagnetic is almost always the interesting part wrt the
c structure. Formula: Tr(D1*H01), both antisymmetric matrices.
c The trace is the elementwise sum, individual elements are the
c contributions themselves.      
c As the matrices are in AO basis, they need to be transformed to
c the MO basis.      
c If there are localized orbitals, it will use them.
c Does not work in parallel yet, most data available only on master now.
c
      use memory
      
      implicit real*8(a-h,o-z)
      character*1 xyz
      dimension shie(3,3,*),con(10),norb(2,10),xyz(3)
      data xyz/'x','y','z'/
      
      ntri=ncf*(ncf+1)/2
      np1=1
      np4=4
      write(6,'("Partitioning of the paramagnetic shielding of atom: ",
     #          i2)') nra
      call mmark
      call matmark
      call getmem(ntri,liso)
      call zeroit(bl(liso),ntri)
      call matdef('smat','s',ncf,ncf)
      call matread('smat',1,'s matrix')
      call matdef('coef','q',ncf,ncf)
      call sudat(np4,'loca_rhf',ni)
      if(ni.gt.0) then
        call matread('coef',4,'loca_rhf')
        write(6,*)
     1  ' Contributions from localized orbitals'
      else
        call matread('coef',4,'evec_rhf')
        write(6,*)
     1  ' Contributions from canonical orbitals'
      endif
      call matdef('traf','q',ncf,ncf)
      call matmmult('smat','coef','traf')
      do 500 l=1,3
      do 500 m=1,3
      if(abs(shie(m,l,nra+nreq)).gt.0.02d0)then
      shiel=100/shie(m,l,nra+nreq)
      call mmark
      call matmark
      call matconn('h01','a',ncf,ncf,lhna+(m-1)*ntri)
      call matconn('d1','a',ncf,ncf,ldn1+(l-1)*ntri)
      call matdef('d1mo','a',ncf,ncf)
      call matdef('h01mo','a',ncf,ncf)
      call matsimtr('d1','traf','d1mo')
      call matsimtr('h01','coef','h01mo')
      ld1mo=mataddr('d1mo')
      lh1mo=mataddr('h01mo')
      ii=0
      call zeroit(con,10)
      lcon=1
      indx=lcon
      do 400 i=1,ncf
      do 400 j=1,i
         ii=ii+1
         contr=bl(ld1mo+ii-1)*bl(lh1mo+ii-1)*2.0d0
         if(m.eq.l) bl(liso+ii-1)=bl(liso+ii-1)+contr*0.333333d0
         if(abs(contr*shiel).gt.5.0d0)then
            do 200 ind=lcon,1,-1
               if (abs(contr).lt.abs(con(ind))) then
                  indx=ind+1
                  goto 250
               end if
 200        end do
            indx=1
 250        if(indx.lt.10) then
                do 300 ltra=9,indx,-1
                   con(ltra+1)=con(ltra)
                   norb(1,ltra+1)=norb(1,ltra)
                   norb(2,ltra+1)=norb(2,ltra)
 300            enddo
                con(indx)=contr
                norb(1,indx)=j
                norb(2,indx)=i
            endif
            if(lcon.lt.10) lcon=lcon+1
         endif
 400  enddo
      write(6,'("  Tensor element ",2a,": ",f8.2,a4)') xyz(m),xyz(l), 
     #       shie(m,l,nra+nreq)*26.6256691d+00,' ppm'
      write(6,'("  Top contributions:")')
      psum=0
      do 100 ltra=1,lcon-1
         write(6,'(4x,f8.2," ppm",i4,"-->",i4)') con(ltra)*26.6256691d0,
     #         norb(1,ltra),norb(2,ltra)
         psum=psum+con(ltra)
 100  enddo
      write(6,'(f8.2," % from these excitations")') psum*shiel
      write(6,*)
      call matremark
      call retmark
      endif
 500  continue
      shiel=(shie(1,1,nra+nreq)+shie(2,2,nra+nreq)+shie(3,3,nra+nreq))/3
      shiso=100/shiel
      ii=0
      call zeroit(con,10)
      lcon=1
      indx=lcon
      do  i=1,ncf
      do  j=1,i
         ii=ii+1
         contr=bl(liso+ii-1)
         if(abs(contr*shiso).gt.2.0d0)then
            do ind=lcon,1,-1
               if (abs(contr).lt.abs(con(ind))) then
                  indx=ind+1
                  goto 550
               end if
            end do
            indx=1
 550           if(indx.lt.10) then
                do ltra=9,indx,-1
                   con(ltra+1)=con(ltra)
                   norb(1,ltra+1)=norb(1,ltra)
                   norb(2,ltra+1)=norb(2,ltra)
                enddo
                con(indx)=contr
                norb(1,indx)=j
                norb(2,indx)=i
            endif
            if(lcon.lt.10) lcon=lcon+1
         endif
      enddo
      enddo
      write(6,'("  Isotropic shielding: ",f8.2)') 
     #       shiel*26.6256691d+00
      write(6,'("  Top contributions:")')
      psum=0
      do ltra=1,lcon-1
         write(6,'(4x,f8.2," ppm",i4,"-->",i4)') con(ltra)*26.6256691d0,
     #         norb(1,ltra),norb(2,ltra)
         psum=psum+con(ltra)
      enddo
      write(6,'(4x,f8.2," % from these excitations")') psum*shiso
      write(6,*)
      call matremark
      call retmark
c
      end
c=======================================================================
      subroutine shnics(nna,bl,inx,last,iprint,iout,datbas,datnics,ncs,
     *                  ncf,ntri,nics,den,dn1,shie,hna,nra)
c-----------------------------------------------------------------------
c        this subroutine calculates 1-el integral derivatives which
c                  are connected with 1/r1n operator.
c
c             there are the following integral derivatives:
c
c 1.   (i00/h11n/j00) + (i10/h01n/j00) + (i00/h01n/j10)     matrix
c
c 2.   (i00 / h01n/j00)                                     matrix
c
c These calculations are performed for ONE nucleus (NRAT) at the time
c-----------------------------------------------------------------------
c nna 	- number of atoms
c bl	- real storage area
c inx	- integer storage
c last	- pointer to top of free mem.
c iprint- print level
c iout	- unit of stdout
c datbas- basis set data struct.
c datnics- nuclear data struct.
c ncs,ncf - no of contr shells and func
c nics	- number of atoms shielding need to be calc. for
c den	-  D00
c dn1	- pointer to D10
c shie	- shieldings (9*nics diam.+9*nics para
c hna	- H01
c nra	- atom worked on (number)
c
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /lindvec/ lind,idensp
      dimension bl(1)
      dimension den(ntri),dn1(ntri,3),hna(ntri,3),shie(3,3,*)
      dimension datbas(13,*),datnics(3,nics)
      dimension inx(12,1)
c needed 12*I*I : 12*28*28
      dimension s2(9408),xra(3)
c
c units for shielding tensor :
c  1 au = 26.6256691d+00 ppm
c
      data f12c2p / 26.6256691d+00 /
c
c remember value of last
      lenter=last
c
c zero out h01,set up pointers to dia and paramagn.
c
      call zeroit(hna,ntri*3)
c
c nra is NICS point of interest 
c
         xra(1)=datnics(1,nra)
         xra(2)=datnics(2,nra)
         xra(3)=datnics(3,nra)
c
      kk=0
      k2=0
c
      key=4
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
            len=len1*len2
            ls2=12*len
            ld=12
c
        do 45 jgc=0,ngcjx
c
      call onesh(i,j,igc,jgc,datbas,bl,s2,inx,key,nna,kk,k2,xra ,ls2,ld)
c
            iijxx=0
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
c
               iijxx=iijxx+1
               iijxy=iijxx+len
               iijxz=iijxy+len
c
               iijyx=iijxz+len
               iijyy=iijyx+len
               iijyz=iijyy+len
c
               iijzx=iijyz+len
               iijzy=iijzx+len
               iijzz=iijzy+len
         if(jff.gt.iff) go to 40
               ij=ii+jff
               d00=den(ij)
               if(iff.ne.jff) d00=two*d00
c
            hna(ij,1)=hna(ij,1) + s2(len*9+iijxx)
            hna(ij,2)=hna(ij,2) + s2(len*9+iijxy)
            hna(ij,3)=hna(ij,3) + s2(len*9+iijxz)
c
            shie(1,1,nra)=shie(1,1,nra) + s2(iijxx)*d00
            shie(2,1,nra)=shie(2,1,nra) + s2(iijxy)*d00
            shie(3,1,nra)=shie(3,1,nra) + s2(iijxz)*d00

            shie(1,2,nra)=shie(1,2,nra) + s2(iijyx)*d00
            shie(2,2,nra)=shie(2,2,nra) + s2(iijyy)*d00
            shie(3,2,nra)=shie(3,2,nra) + s2(iijyz)*d00
c
            shie(1,3,nra)=shie(1,3,nra) + s2(iijzx)*d00
            shie(2,3,nra)=shie(2,3,nra) + s2(iijzy)*d00
            shie(3,3,nra)=shie(3,3,nra) + s2(iijzz)*d00
c
   40       continue
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c
c return value of last
      last=lenter
c
      if(iprint.gt.2) then
        write(iout,444)
        write(iout,445) nra
        call druma(hna(1,1),ncf,iout,'h01mat x ')
        call druma(hna(1,2),ncf,iout,'h01mat y')
        call druma(hna(1,3),ncf,iout,'h01mat z')
      endif
c
 444  format('*******************************')
 445  format(/'***  total matrix h01 for atom number  ',i2,'     ***')
c
c Calculate the paramagnetic part of the shielding tensor
c
c the trace of the product of two antisymmetric matrices is really
c simple  Tr(AB)= - Aik Bik   (with the Einstein convention)
c in the present case we also have a factor of -1
c as both matrices are imaginary
c if both square matrices are stored in the triangular form
c the whole thing is a simple blas dot product (the diagonal is zero!)
c
      do 100 l=1,3
      do 100 m=1,3
         shie(m,l,nra+nics)=ddot(ntri,dn1(1,l),1,hna(1,m),1)*two
 100  continue
c
c Convert diamagnetic and paramagnetic shielding
c from au to ppm
c
      call dscal(9,f12c2p,shie(1,1,nra),1)
      call dscal(9,f12c2p,shie(1,1,nra+nics),1)
      end
c=======================================================================
      subroutine calcsnics(NICS,igeoc,shie,datnix,bl,IOut,ICon,iprint)
      implicit real*8(a-h,o-z)
c--------------------------------------------------------------------
C  Calculates shielding in principal coordinate system for atoms
c  requested
C  ** WARNING **
C     direction cosines may be have degenerate components is
C     symmetrical molecules; in this case the degenerate
c     components are ill-defined
C
C  ARGUMENTS
C
c  nics    -  number of NICS points
c  shie    -  shielding tensor 
c  datnix  -  (3,nics) coordinates of nics points
C  bl      -  main storage array
C  IOut    -  unit for standard output
C  ICon    -  unit for summary output
C  iprint  -  print level
C
c--------------------------------------------------------------------
      dimension datnix(3,nics)
      dimension bl(*),shie(3,3,*)
      dimension prs(3),shs(3,3)
c
      Parameter (Zero=0.0d0,half=0.5d0,three=3.0d0,one=1.0d0)
c--------------------------------------------------------------------
c
c -- initial print headers
      write(iout,444)
      write(iout,611)
      write(iout,612)
      write(iout,613)
      write(iout,444)
      write(icon,444)
ccc   write(icon,'( 30x,"NMR SHIELDINGS")' )
CCC   write(icon,'( 30x,"NICS (SHIELDINGS)")' )
      write(icon,608) nics,igeoc
 608  format(24x,'NICS (SHIELDINGS) at ',i5,' points'/
     *       24x,'grid center at point ',i5)
      write(icon,444)
      write(icon,609) 
 609  format(35x,'coordinates',15x,' shielding')
c
      call tstrval('xlvsh',iyes)
      if(iyes.eq.1) then
        xlvsh=rgetrval('xlvsh')
        if(abs(xlvsh).gt.1.0D-8) then
          write(iout,*) 'Level shift=',xlvsh
          write(icon,*) 'Level shift=',xlvsh
        end if
      end if
      write(icon,444)
c
c
c     sumsh=0.d0
      do 600 IAtm=1,NICS
c
         nra = iatm               ! position in list of symmetry-unique atoms
         nshw = nra + nics
c
         do 100 j=1,3
         do 100 i=1,3
c -- symmetric part of shielding tensor
           shs(i,j)=(shie(i,j,nra)+shie(j,i,nra)+shie(i,j,nshw)+
     *              shie(j,i,nshw))*half
 100    continue
c
c -- diagonalize shielding tensor
c
         nd=3
         call sdiag2(nd,nd,shs,prs,shs)
c
         averag=(prs(1)+prs(2)+prs(3))/three
         anisot= prs(3)-(prs(2)+prs(1))*half
c
c -- final printout
c
c Print out diamagnetic, paramagnetic and total
c shielding tensor 
c
            if(iprint.ge.1) then
               call printsh(bl,IAtm,shie(1,1,nshw),
     *                      shie(1,1,nra), iout,'blank')
            endif
            call daxpy(9,one,shie(1,1,nshw),1,shie(1,1,nra),1)
            call zeroit(shie(1,1,nshw),9)
            if(iprint.ge.1) then
               call printsh(bl,IAtm,shie(1,1,nshw),
     *                      shie(1,1,nra), iout,'total')
cc             call printsh(bl,IAtm,shie(1,1,nshw),
cc   *                      shie(1,1,nra), icon,'short')
            endif
      write(iout,614) IAtm,datnix(1,iatm),datnix(2,iatm),datnix(3,iatm)
      write(icon,618) IAtm,datnix(1,iatm),datnix(2,iatm),
     $                datnix(3,iatm),averag
         write(iout,615)
         write(iout,610) prs(1),(shs(i,1),i=1,3)
         write(iout,610) prs(2),(shs(i,2),i=1,3)
         write(iout,610) prs(3),(shs(i,3),i=1,3)
         write(iout,616) averag,anisot
         write(iout,444)
c
c        sumsh=sumh+averag
c
 600     continue
c
c----------------------------------------------------------------------
c        write(icon,*) 'integ. sh=',sumsh
c----------------------------------------------------------------------
         write(icon,444)
 444  format(72(1H-))
 610  format( '--  ',F13.5,12X,3(F12.7,1X),'  --')
 611  format('--',16X,'     The Total Shielding Tensor',21X,'--')
 612  format('--',9X,'in principal axes and direction cosines between',
     *12X,'--')
 613  format('--',10X,' molecular and laboratory coordinate systems',14X
     1,'--')
 614  format
     *('--  Shielding',' at NICS',' point ',I5,
     *' (',3(f10.5,1x),') --')
 615  format(72(1H-))
 616  format(/'-- Average =',F12.5,' , Anisotropy=',F12.5,20X,'--')
 618  format
     *('-- NICS at point ',I5,
     *' (',3(f10.5,1x),')',1x,f10.3,' --')
c-----------------------------------------------------------------------
      end
c======================================================================
      subroutine find_nicsgrid(iat1,iat2,iat3,datnuc,nat,disp,ngrid,
     *                         step,datnic,nics,igeoc)
      implicit real*8 (a-h,o-z)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      logical linear
      dimension datnuc(5,*)
      dimension datnic(3,*)
      dimension p1(3),p2(3),p3(3), vnor(3),plane(4)
      dimension itonplane(nat)
      dimension itonplan1(nat)
      data toler /1.0d-6/
c-----------------------------------------------------------------------
c
c iat1,iat2,iat3 - atoms that define plane P for NICS 
c-----------------------------------------------------------------------
c if NON of iat1,iat2,iat3 is 0 then grid over a plane (square)
c if ONE of iat1,iat2,iat3 is 0 then grid along line    
c-----------------------------------------------------------------------
c disp           - displace P by disp making plane P' parallel to P
c ngrid          - gives number of grid NICS points on a given plane 
c                  number of points for NICS is ngrid*ngrid (+1)
c                  so if ngrid is 0 then only one point (geometrical
c                  center of P plane) will be used for NICS
c step          - step size for a grid
c nat           - number of atoms in a molecule
c datnic        - coordinates of NICS points
c nics          - number of NICS points
c
c igeoc         - central point of a grid
c-----------------------------------------------------------------------
      call getival('iout',iout)
c-----------------------------------------------------------------------
      natmin=min(iat1,iat2,iat3)
      natomx=max(iat1,iat2,iat3)
c-----------------------------------------------------------------------
      IF(natomx.gt.nat.and.natmin.gt.0) then
         write(iout,*)' specified atom no ',natomx,
     *                ' does not belong to molecule'
         write(iout,*)' atoms for NICS plane redefined'
         if(nat.eq.1) then
            iat1=1
            iat2=1
            iat3=1
         endif
         if(nat.eq.2) then
            iat1=1
            iat2=2
            iat3=2
         endif
         if(nat.ge.3) then
            iat1=1
            iat2=2
            iat3=3
         endif
      endif
c-----------------------------------------------------------------------
      if(natmin.eq.0) then
         if(iat1.gt.0.and.iat2.gt.0.and.iat1.ne.iat2) then
            p1(1)=datnuc(2,iat1)
            p1(2)=datnuc(3,iat1)
            p1(3)=datnuc(4,iat1)
            p2(1)=datnuc(2,iat2)
            p2(2)=datnuc(3,iat2)
            p2(3)=datnuc(4,iat2)
            xg=0.5d0*(datnuc(2,iat1)+datnuc(2,iat2))
            yg=0.5d0*(datnuc(3,iat1)+datnuc(3,iat2))
            zg=0.5d0*(datnuc(4,iat1)+datnuc(4,iat2))
            go to 300
         endif
         if(iat1.gt.0.and.iat3.gt.0.and.iat1.ne.iat3) then
            p1(1)=datnuc(2,iat1)
            p1(2)=datnuc(3,iat1)
            p1(3)=datnuc(4,iat1)
            p2(1)=datnuc(2,iat3)
            p2(2)=datnuc(3,iat3)
            p2(3)=datnuc(4,iat3)
            xg=0.5d0*(datnuc(2,iat1)+datnuc(2,iat3))
            yg=0.5d0*(datnuc(3,iat1)+datnuc(3,iat3))
            zg=0.5d0*(datnuc(4,iat1)+datnuc(4,iat3))
            go to 300
         endif
         if(iat2.gt.0.and.iat3.gt.0.and.iat2.ne.iat3) then
            p1(1)=datnuc(2,iat2)
            p1(2)=datnuc(3,iat2)
            p1(3)=datnuc(4,iat2)
            p2(1)=datnuc(2,iat3)
            p2(2)=datnuc(3,iat3)
            p2(3)=datnuc(4,iat3)
            xg=0.5d0*(datnuc(2,iat2)+datnuc(2,iat3))
            yg=0.5d0*(datnuc(3,iat2)+datnuc(3,iat3))
            zg=0.5d0*(datnuc(4,iat2)+datnuc(4,iat3))
            go to 300
         endif
c
c        if it gets here it means that two "nonzero" atoms are the same
c
         if(nat.eq.1) then
            p1(1)=datnuc(2,iat1)
            p1(2)=datnuc(3,iat1)
            p1(3)=datnuc(4,iat1)
c
            xg=p1(1)
            yg=p1(2)
            zg=p1(3)
c
            pnorm=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
c 
            if(pnorm.le.tol) then
c            point p1 at (0,0,0)
               p2(1)=0.0d0
               p2(2)=0.0d0
               p2(3)=1.0d0
            else
               p2(1)=0.0d0
               p2(2)=0.0d0
               p2(3)=0.0d0
            endif
            go to 300
         endif
         if(nat.ge.2) then
            p1(1)=datnuc(2,1)
            p1(2)=datnuc(3,1)
            p1(3)=datnuc(4,1)
            p2(1)=datnuc(2,2)
            p2(2)=datnuc(3,2)
            p2(3)=datnuc(4,2)
            xg=0.5d0*(datnuc(2,2)+datnuc(2,1))
            yg=0.5d0*(datnuc(3,2)+datnuc(3,1))
            zg=0.5d0*(datnuc(4,2)+datnuc(4,1))
            go to 300
         endif
c
      endif
c-----------------------------------------------------------------------
c check if all atoms are on a line 
c
      call check_linearity(nat,datnuc,p1,p2,p3,linear)
c
cc     write(6,*)' linearity=',linear
c-----------------------------------------------------------------------
      if(linear  ) then
c
         call add_points(nat,datnuc,p1,p2,p3)
c
c        write(6,*)' P1=',p1
c        write(6,*)' P2=',p2
c        write(6,*)' P3=',p3
c
c remember 3 points 
c
            x1=p1(1)
            y1=p1(2)
            z1=p1(3)
c
            x2=p2(1)
            y2=p2(2)
            z2=p2(3)
c
            x3=p3(1)
            y3=p3(2)
            z3=p3(3)
c
c find plane
            p2(1)=p2(1)-p1(1)
            p2(2)=p2(2)-p1(2)
            p2(3)=p2(3)-p1(3)
c
            p3(1)=p3(1)-p1(1)
            p3(2)=p3(2)-p1(2)
            p3(3)=p3(3)-p1(3)
c
            call cross(p2,p3,vnor)
            d=-vnor(1)*p1(1)-vnor(2)*p1(2)-vnor(3)*p1(3)
            plane(1)=vnor(1)            
            plane(2)=vnor(2)            
            plane(3)=vnor(3)            
            plane(4)=d
c central point 
            sumx=0.d0
            sumy=0.d0
            sumz=0.d0
            do iat=1,nat   
              itonplane(iat)=iat
              x=datnuc(2,iat)
              y=datnuc(3,iat)
              z=datnuc(4,iat)
              sumx=sumx+x
              sumy=sumy+y
              sumz=sumz+z
            enddo
            deno=1.0d0/dble(nat)
            xg=sumx*deno
            yg=sumy*deno
            zg=sumz*deno
c first point on a plane
            if(nat.eq.1) then
               p1(1)=x2
               p1(2)=y2
               p1(3)=z2
            endif
            if(nat.ge.3) then
c              check if G and P1 are not the same
               x1g=abs(p1(1)-xg)
               y1g=abs(p1(2)-yg)
               z1g=abs(p1(3)-zg)
               if(x1g.le.tole.and.y1g.le.toler.and.z1g.le.toler) then
                  p1(1)=x2
                  p1(2)=y2
                  p1(3)=z2
               endif
            endif
c
            maxat=nat
         go to 200
      endif
c-----------------------------------------------------------------------
      idiff=3
      IF(iat1.EQ.iat2) idiff=2
      IF(iat1.EQ.iat3) idiff=2
      IF(iat2.EQ.iat3) idiff=2
      IF(iat1.EQ.iat2 .AND. iat2.EQ.iat3) idiff=1
c
c
      IF(idiff.eq.3) THEN
c           find plane a*x + b*y + c*z + d = 0 that contains these 3 atoms
c
            p1(1)=datnuc(2,iat1)
            p1(2)=datnuc(3,iat1)
            p1(3)=datnuc(4,iat1)
c
            p2(1)=datnuc(2,iat2)-p1(1)
            p2(2)=datnuc(3,iat2)-p1(2)
            p2(3)=datnuc(4,iat2)-p1(3)
c
            p3(1)=datnuc(2,iat3)-p1(1)
            p3(2)=datnuc(3,iat3)-p1(2)
            p3(3)=datnuc(4,iat3)-p1(3)
c
            call cross(p2,p3,vnor)
c
            if(abs(vnor(1)).lt.toler.and.
     *         abs(vnor(2)).lt.toler.and.
     *         abs(vnor(3)).lt.toler) then
                  write(iout,*) 'specified atoms lay on a line'
                  write(iout,*) '  program will find a plane'
                  write(iout,*) '    with most atoms on it'
                  go to 10
            endif
c
            d=-vnor(1)*p1(1)-vnor(2)*p1(2)-vnor(3)*p1(3)
c
            ile=0
            do iat=1,nat
               dat=-vnor(1)*datnuc(2,iat)
     *             -vnor(2)*datnuc(3,iat)
     *             -vnor(3)*datnuc(4,iat)
               if(abs(dat-d).le.toler) then
                  ile=ile+1
                  itonplane(ile)=iat
                  maxat=ile
               endif
cccc        write(6,*)'iat=',iat,' ile=',ile,' maxat=',maxat
            enddo
c
            plane(1)=vnor(1)            
            plane(2)=vnor(2)            
            plane(3)=vnor(3)            
            plane(4)=d
c
            go to 100
      ENDIF
c-----------------------------------------------------------------------
   10 continue
c-----------------------------------------------------------
c
c     find plane a*x + b*y + c*z + d = 0 that contains most atoms
c
      i3max=1
      maxat=0
      i3=0
      do i=1,nat
         p1(1)=datnuc(2,i)
         p1(2)=datnuc(3,i)
         p1(3)=datnuc(4,i)
ccccc    do j=1,i-1
         do j=i+1,nat
            p2(1)=datnuc(2,j)-p1(1)
            p2(2)=datnuc(3,j)-p1(2)
            p2(3)=datnuc(4,j)-p1(3)
ccccc       do k=1,j-1
            do k=j+1,nat
               p3(1)=datnuc(2,k)-p1(1)
               p3(2)=datnuc(3,k)-p1(2)
               p3(3)=datnuc(4,k)-p1(3)
c
               call cross(p2,p3,vnor)
c
               if(abs(vnor(1)).lt.toler.and.
     *            abs(vnor(2)).lt.toler.and.
     *            abs(vnor(3)).lt.toler) go to 20
c
               d=-vnor(1)*p1(1)-vnor(2)*p1(2)-vnor(3)*p1(3)
c
               ile=0
               do iat=1,nat
                  dat=-vnor(1)*datnuc(2,iat)
     *                -vnor(2)*datnuc(3,iat)
     *                -vnor(3)*datnuc(4,iat)
                  if(abs(dat-d).le.toler) then
                     ile=ile+1
                     itonplan1(ile)=iat
                  endif
               enddo
c
               i3=i3+1
               if(ile.gt.maxat) then
                  maxat=ile
                  i3max=i3
                  plane(1)=vnor(1)            
                  plane(2)=vnor(2)            
                  plane(3)=vnor(3)            
                  plane(4)=d
                  do iat=1,ile
                     itonplane(iat)=itonplan1(iat)
                  enddo
               endif
c
   20          continue
c
            enddo
         enddo
      enddo
c-----------------------------------------------------------
c setup the first point on a plane
c
         iatp1=itonplane(1)
         p1(1)=datnuc(2,iatp1)
         p1(2)=datnuc(3,iatp1)
         p1(3)=datnuc(4,iatp1)
c-----------------------------------------------------------
  100 continue
c-----------------------------------------------------------
c
c     we have a plane P:  a*x + b*y + c*z + d = 0 for NICS points
c     this palne contains maxat atoms
c
c norm of the normal vector
c
c     xnorm=plane(1)**2+plane(2)**2+plane(3)**2
c     xnorm=sqrt(xnorm)
c
c     t=disp/xnorm 
c
c    ! point on P' parallel to P distant by disp
c
c       xp'=xp+t*a
c       yp'=yp+t*b
c       zp'=zp+t*c
c
c-----------------------------------------------------------
c find geometrical center of MAXAT atoms on the P plane
c
      sumx=0.d0
      sumy=0.d0
      sumz=0.d0
      do iatp=1,maxat
        iat=itonplane(iatp)
        x=datnuc(2,iat)
        y=datnuc(3,iat)
        z=datnuc(4,iat)
        sumx=sumx+x
        sumy=sumy+y
        sumz=sumz+z
      enddo
      deno=1.0d0/dble(maxat)
      xg=sumx*deno
      yg=sumy*deno
      zg=sumz*deno
c
c-----------------------------------------------------------
  200 continue
c-----------------------------------------------------------
c we have a plane P:  a*x + b*y + c*z + d = 0 
c norm of the normal vector
c
      xnorm=plane(1)**2+plane(2)**2+plane(3)**2
      xnorm=sqrt(xnorm)
c
      t=disp/xnorm 
c
c    ! point on P' parallel to P distant by disp
c
c       xp'=xp+t*a
c       yp'=yp+t*b
c       zp'=zp+t*c
c-----------------------------------------------------------
c plane P' parallel to P and distant from P by disp
c
      A=plane(1)
      B=plane(2)
      C=plane(3)
      D0=plane(4)
c
c geometrical center on P
      xg0=xg
      yg0=yg
      zg0=zg
c
c geometrical center on P' (image of x0)
      xg=xg0+t*a
      yg=yg0+t*b
      zg=zg0+t*c
c
      d=-(a*xg+b*yg+c*zg)
      plane(4)=d
c-----------------------------------------------------------
c setup atom no 1 (P1) for grid on plane
c image of atom1 on P'
c
      p1(1)=p1(1)+t*a
      p1(2)=p1(2)+t*b
      p1(3)=p1(3)+t*c
c-----------------------------------------------------------
c     write(6,*)' point G=',xg,yg,zg
c     write(6,*)' point P=',p1(1),p1(2),p1(3)
c-----------------------------------------------------------
      write(iout,5001)
 5001 format(/72('='))
      write(iout,301)
  301 format(
     * '----     Calculations of the Nuclear',
     * ' Independent Chemical Shift     ----')
      write(iout,302)
  302 format(
     * '----                         N I C S',
     * '                                ----')
      write(iout,303)
  303 format(
     * '----          Magnetic Shielding at the points',
     * ' on the plane         ----')
      write(iout,444)
c-----------------------------------------------------------
      write(iout,72) a,b,c,d    
  72  format(/'NICS palne P  : ',
     *f9.5,'*X + ',f9.5,'*Y + ',f9.5,'*Z + ',f9.5,' = 0')
      write(iout,73) a,b,c,d0, maxat,(itonplane(i),i=1,maxat)
  73  format(/'parallel to R : ',
     *f9.5,'*X + ',f9.5,'*Y + ',f9.5,'*Z + ',f9.5,' = 0'//,
     *'containing ',i4,' atoms : ',10(i3,1x))
      write(iout,74) xg0,yg0,zg0
  74  format(/'geoemtrical center of atoms on R       at:',3(f9.5,1x))
      write(iout,75) xg,yg,zg
  75  format(/'geoemtrical center of atom images on P at:',3(f9.5,1x))
      write(iout,76) disp
  76  format(/'NICS plane P distant from molecular plane R by ',f9.5)
c-----------------------------------------------------------
c The NICS plane
c 
      A=plane(1)
      B=plane(2)
      C=plane(3)
      D=plane(4)
c
c-----------------------------------------------------------
  300 continue   ! jump here if grid along line
c-----------------------------------------------------------
c setup points for NICS calculations
c
      nics=1
      datnic(1,1)=xg
      datnic(2,1)=yg
      datnic(3,1)=zg
c
      h= step 
c
c-----------------------------------------------------------
      if(ngrid.gt.1) then
         if(natmin.eq.0) then
            call setup_nic1(datnic,p1,p2,p3,vnor,h,ngrid,
     *                      disp,nics,igeoc)
         else
            call setup_nic2(datnic,p1,p2,p3,vnor,h,ngrid,
     *                      a,b,c,d,nics,igeoc)
         endif
      else
         igeoc=1
         nics=1
      endif
c-----------------------------------------------------------
      write(iout,*)' '
      write(iout,444)
 444  format(72(1H-))
c
      if(natmin.eq.0) then
         write(iout,77) ngrid,step, nics,igeoc
      else
         write(iout,78) ngrid,step, nics,igeoc
      endif
  77  format('          GRID for NICS along line (2*N+1)  with N=',i5,/
     *       '                   grid step size is     =',f5.3,' bohr'/
     *       '                   Number of NICS points =',i5,/
     *       '                   Central grid point no =',i5/)
  78  format('          GRID for NICS over square (N x N) with N=',i5,/
     *       '                   grid step size is     =',f5.3,' bohr'/
     *       '                   Number of NICS points =',i5,/
     *       '                   Central grid point no =',i5/)
c-----------------------------------------------------------
c
      end
c======================================================================
      subroutine setup_nic1(datnic,p1,p2,p3,vnor,h,ngrid,dis,nics,igeoc)
      implicit real*8 (a-h,o-z)
      dimension datnic(3,*) 
      dimension p1(3),p2(3)    ! points for a line
      dimension p3(3), vnor(3) ! working vectors
      data tol /1.0d-6/
c-----------------------------------------------------------
c  datnic(3,i)  coordinates on i-nics point
c  p1(3),p2(3)  coordinates of two points defining a line
c  p3,vnor    - working vectors
c  h          - step size for nics grid
c  ngrid      - number of points on a square , total N*N
c  dis        - displace a line by dis
c  nics         number of NICS grid points
c-----------------------------------------------------------
c
      call getival('iout',iout)
c
c geometrical center of the nics grid on a line
c
      xg=datnic(1,1)
      yg=datnic(2,1)
      zg=datnic(3,1)
c
      a=p1(1)-p2(1)
      b=p1(2)-p2(2)
      c=p1(3)-p2(3)
c
      abc=sqrt(a*a+b*b+c*c)
c
      write(6,*)'a,b,c=',a,b,c
c
c-----------------------------------------------------------
c make requested displacement of a line by dis
c
        wx= b
        wy=-a
        wz= 0.0d0
        wxyz=sqrt(wx*wx+wy*wy+wz*wz)
        if(wxyz.le.tol) then
           wx= c
           wy= 0.0d0
           wz= -a
           wxyz=sqrt(wx*wx+wy*wy+wz*wz)
        endif
c
        disp=dis/wxyz
c
        wx=wx*disp
        wy=wy*disp
        wz=wz*disp
c
        xg=xg+wx
        yg=yg+wy
        zg=zg+wz
c-----------------------------------------------------------
c

      t=h/abc
c
c  points online p1-p2
c
      dx=t*a
      dy=t*b
      dz=t*c
c
c  begining 
c
      x0=xg-ngrid*dx
      y0=yg-ngrid*dy
      z0=zg-ngrid*dz
c
      do i=1,2*ngrid+1
         dxi=(i-1)*dx
         dyi=(i-1)*dy
         dzi=(i-1)*dz
         x1=x0+dxi
         y1=y0+dyi
         z1=z0+dzi
         datnic(1,i)=x1
         datnic(2,i)=y1
         datnic(3,i)=z1
      enddo
      igeoc=ngrid+1
      nics=2*ngrid+1
      datnic(1,ngrid+1)=xg
      datnic(2,ngrid+1)=yg
      datnic(3,ngrid+1)=zg
c
      end
c======================================================================
      subroutine setup_nic2(datnic,p1,p2,p3,vnor,h,ngrid,a,b,c,d,nics,
     *                      igeoc)
      implicit real*8 (a-h,o-z)
      dimension datnic(3,*) 
      dimension p1(3)  ! atom1 on a plane  
      dimension p2(3), p3(3), vnor(3) ! working vectors
c-----------------------------------------------------------
c  datnic(3,i)  coordinates on i-nics point
c  p1(3)        coordinates of the first point def. NICS palne
c  p2,p3,vnor - working vectors
c  h          - step size for nics grid
c  ngrid      - number of points on a square , total N*N
c  nics         number of NICS grid points
c  a,b,c,d    - plane
c
c-----------------------------------------------------------
      step=h
c-----------------------------------------------------------
c geometrical center of atoms on main plane Ax+By+Cz+D=0 
c
      xg=datnic(1,1)
      yg=datnic(2,1)
      zg=datnic(3,1)
c-----------------------------------------------------------
c GP1(xp1-xg,yp1-yg,zp1-zg) gives direction of gp1 line
c find point P on gp1 distant from G by RR=ngrid*h/sqrt(2)
c
      rr=dble(ngrid-1)*h/sqrt(2.d0)
c
      p2(1)=p1(1)-xg
      p2(2)=p1(2)-yg
      p2(3)=p1(3)-zg
c
      gP1norm= sqrt(p2(1)*p2(1)+p2(2)*p2(2)+p2(3)*p2(3))
      t=rr/gP1norm
c
c point P is
c
      xp=xg+t*p2(1)
      yp=yg+t*p2(2)
      zp=zg+t*p2(3)
c-----------------------------------------------------------
c cross product of GP1(xp-xg,yp-yg,zp-zg) and normal (A,B,C)
c gives direction of a line GR1 perpendicualr to GP1
c find point R on GR1 distant from G by RR=sqrt(2)ngrid*h
c
c
      p3(1)=a
      p3(2)=b
      p3(3)=c
c
      call cross(p2,p3,vnor)
c
      gR1norm= sqrt(vnor(1)*vnor(1)+vnor(2)*vnor(2)+vnor(3)*vnor(3))
      t=rr/gR1norm
c
c point R is
c
      xr=xg+t*vnor(1)
      yr=yg+t*vnor(2)
      zr=zg+t*vnor(3)
c-----------------------------------------------------------
c  P and R distant by R=ngrid*step
c-----------------------------------------------------------
c find direction of a line PQ perpendicular to PR
c
c
      p2(1)=xp-xr
      p2(2)=yp-yr
      p2(3)=zp-zr
c
      call cross(p2,p3,vnor)
c
      PQnorm= sqrt(vnor(1)*vnor(1)+vnor(2)*vnor(2)+vnor(3)*vnor(3))
      t=sqrt(2.d0)*rr/PQnorm
c
c point Q is
c
      xq=xp-t*vnor(1)
      yq=yp-t*vnor(2)
      zq=zp-t*vnor(3)
c-----------------------------------------------------------
      p3(1)=xp-xq
      p3(2)=yp-yq
      p3(3)=zp-zq
c-----------------------------------------------------------
c check distance PR
c     dpg=sqrt((xg-xp)**2+(yg-yp)**2+(zg-zp)**2)
c     drg=sqrt((xg-xr)**2+(yg-yr)**2+(zg-zr)**2)
c     dqg=sqrt((xg-xq)**2+(yg-yq)**2+(zg-zq)**2)
c     write(6,*)' Pg,Rg,Qg=',dpg,drg,dqg
c     dpr=sqrt((xr-xp)**2+(yr-yp)**2+(zr-zp)**2)
c     write(6,*)' PR=',dpr,'  RR=',sqrt(2.d0)*rr
c     dpq=sqrt((xq-xp)**2+(yq-yp)**2+(zq-zp)**2)
c     write(6,*)' PQ=',dpq,'  RR=',sqrt(2.d0)*rr
c-----------------------------------------------------------
c     write(6,99)' P :',xp,yp,zp
c     write(6,99)' R :',xr,yr,zr
c     write(6,99)' Q :',xq,yq,zq
c 99  format(a5,3(f10.5,1x))
c-----------------------------------------------------------
      t=1.d0/dble(ngrid-1)
c
c  points online PR 
c
      xpr=p2(1)
      ypr=p2(2)
      zpr=p2(3)
c
      xpq=p3(1)
      ypq=p3(2)
      zpq=p3(3)
c
      dxpr=-t*xpr
      dypr=-t*ypr
      dzpr=-t*zpr
c
      dxpq=-t*xpq
      dypq=-t*ypq
      dzpq=-t*zpq
c
c-----------------------------------------------------------
      do i=1,ngrid
         ii=(i-1)*ngrid+1
         x1=xp+(i-1)*dxpr 
         y1=yp+(i-1)*dypr 
         z1=zp+(i-1)*dzpr 
         datnic(1,ii)=x1
         datnic(2,ii)=y1
         datnic(3,ii)=z1
         do j=1,ngrid
            x2=x1+(j-1)*dxpq
            y2=y1+(j-1)*dypq
            z2=z1+(j-1)*dzpq
            nics=nics+1
            ij=ii-1+j
            datnic(1,ij  )=x2
            datnic(2,ij  )=y2
            datnic(3,ij  )=z2
         enddo
      enddo
c
c  if ngrid is odd then the geometry ceneter is one of the
c  already geneterted grid points 
c
      n2=ngrid/2
      if(2*n2.NE.ngrid) then
        nics=nics-1
        ij=n2+1
        igeoc=(ij-1)*ngrid+ij
      else
         datnic(1,nics)=xg
         datnic(2,nics)=yg
         datnic(3,nics)=zg
         igeoc=nics
      endif
c-----------------------------------------------------------
      end
c======================================================================
      subroutine add_points(nat,datnuc,p1,p2,p3)
      implicit real*8 (a-h,o-z)
      dimension datnuc(5,*)
      dimension p1(3),p2(3),p3(3)
      data toler /1.0d-6/
c
c
c this routine is called only if a system was found linear
c
c
      call getival('iout',iout)
c
      p1(1)=datnuc(2,1)
      p1(2)=datnuc(3,1)
      p1(3)=datnuc(4,1)
c
      if(nat.eq.1) then
         if(abs(p1(1)).le.toler.and.
     *      abs(p1(2)).le.toler.and.
     *      abs(p1(3)).le.toler) then
            p2(1)= 1.0d0
            p2(2)= 0.0d0
            p2(3)= 0.0d0
c
            p3(1)= 0.0d0
            p3(2)= 1.0d0
            p3(3)= 0.0d0
c
cccc            write(6,*)'x,y,z=0'
c          
            RETURN
         endif
c
c
         if(abs(p1(1)).gt.toler.and.
     *      abs(p1(2)).gt.toler.and.
     *      abs(p1(3)).gt.toler) then
            p2(1)= 1.0d0
            p2(2)= 0.0d0
            p2(3)= 0.0d0
c
            p3(1)= 0.0d0
            p3(2)= 1.0d0
            p3(3)= 0.0d0
c
cccc            write(6,*)'x,y,z>0'
c          
            RETURN
         endif
c
c
         if(abs(p1(1)).le.toler) then
            p2(1)= 0.0d0
            p2(2)=p1(2)+1.0d0
            p2(3)=p1(3)
c     
            p3(1)= 0.0d0
            p3(2)=p1(2)
            p3(3)=p1(3)+1.0d0
cccc            write(6,*)'x    =0'
c          
            RETURN
         endif
c
         if(abs(p1(2)).le.toler) then
            p2(2)= 0.0d0
            p2(1)=p1(1)+1.0d0
            p2(3)=p1(3)
c     
            p3(2)= 0.0d0
            p3(1)=p1(1)
            p3(3)=p1(3)+1.0d0
cccc            write(6,*)'  y  =0'
c          
            RETURN
         endif
c
         if(abs(p1(3)).le.toler) then
            p2(3)= 0.0d0
            p2(1)=p1(1)+1.0d0
            p2(2)=p1(2)
c     
            p3(3)= 0.0d0
            p3(1)=p1(1)
            p3(2)=p1(2)+1.0d0
cccc            write(6,*)'    z=0'
c          
            RETURN
         endif
      endif        !(nat.eq.1) then
c
c two atoms
c
      p2(1)=datnuc(2,2)
      p2(2)=datnuc(3,2)
      p2(3)=datnuc(4,2)
c
      x12=p1(1)-p2(1)
      y12=p1(2)-p2(2)
      z12=p1(3)-p2(3)
      if(nat.eq.2) then
         if(abs(x12).le.toler) then
            p3(1)=p1(1)+1.0d0
            p3(2)=p1(2)
            p3(3)=p1(3)
            RETURN
         endif
         if(abs(y12).le.toler) then
            p3(1)=p1(1)
            p3(2)=p1(2)+1.0d0
            p3(3)=p1(3)
            RETURN
         endif
         if(abs(z12).le.toler) then
            p3(1)=p1(1)
            p3(2)=p1(2)
            p3(3)=p1(3)+1.d0
            RETURN
         endif
         if(abs(x12).gt.toler.and.
     *      abs(y12).gt.toler.and.
     *      abs(z12).gt.toler)     then
            p3(1)=p1(1)+1.0d0
            p3(2)=p1(2)+1.0d0
            p3(3)=p1(3)+1.d0
            RETURN
         endif
      endif        !(nat.eq.2) then
c
c three or more atoms
c
            x1=p1(1)
            y1=p1(2)
            z1=p1(3)
c
            x2=p2(1)
            y2=p2(2)
            z2=p2(3)
c
            x3=datnuc(2,3)
            y3=datnuc(3,3)
            z3=datnuc(4,3)
c
            p3(1)=datnuc(2,3)
            p3(2)=datnuc(3,3)
            p3(3)=datnuc(4,3)
c
c           write(6,*)' point 1=',p1
c           write(6,*)' point 2=',p2
c           write(6,*)' point 3=',p3

c
            p2(1)=p2(1)-p1(1)
            p2(2)=p2(2)-p1(2)
            p2(3)=p2(3)-p1(3)
c
            p3(1)=p3(1)-p1(1)
            p3(2)=p3(2)-p1(2)
            p3(3)=p3(3)-p1(3)
c
cccc         write(6,*)' vector 12=',p2
cccc         write(6,*)' vector 13=',p3
c
            call cross(p2,p3, p1 )
c
cccc      write(6,*)' normal p1=',p1
c
            if(abs(p1(1)).le.toler.and.
     *         abs(p1(2)).le.toler.and.
     *         abs(p1(3)).le.toler) then
c                 write(iout,*) 'specified atoms lay on a line'
c                 write(iout,*) 'program will add extra point'
c
               p1(1)=x1
               p1(2)=y1
               p1(3)=z1
c
               p2(1)=x2
               p2(2)=y2
               p2(3)=z2
c this is extra point
               p3(1)=p1(1)+1.0d0
               p3(2)=p1(2)+1.0d0
               p3(3)=p1(3)
cccc           p3(3)=p1(3)+1.0d0
            else
c
               p1(1)=x1
               p1(2)=y1
               p1(3)=z1
c
               p2(1)=x2
               p2(2)=y2
               p2(3)=z2
c
               p3(1)=x3
               p3(2)=y3
               p3(3)=z3
c
            endif
c
      end
c======================================================================
      subroutine check_linearity(nat,datnuc,p1,p2,p3,linear)
      implicit real*8 (a-h,o-z)
      logical linear
      dimension datnuc(5,*)
      dimension p1(3),p2(3),p3(3)
      data tol /1.0d-6/
c
      if(nat.le.2) then
         linear=.true.
         return
      endif
c
      call getival('iout',iout)
c
      p1(1)=datnuc(2,1)
      p1(2)=datnuc(3,1)
      p1(3)=datnuc(4,1)
c
      p2(1)=datnuc(2,2)
      p2(2)=datnuc(3,2)
      p2(3)=datnuc(4,2)
c
      a=p1(1)-p2(1)
      b=p1(2)-p2(2)
      c=p1(3)-p2(3)
c
c       write(6,*)'a,b,c=',a,b,c
c
      x1=p1(1)
      y1=p1(2)
      z1=p1(3)
c
      linear=.true.
      if(abs(a).gt.tol.and.abs(b).gt.tol.and.abs(c).gt.tol) then
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
            yp1=y-y1  
            zp1=z-z1
c
            tx=xp1/a
            ty=yp1/b
            tz=zp1/c
c
            txy=abs( tx/ty-1.d0 )
            txz=abs( tx/tz-1.d0 )
            tyz=abs( ty/tz-1.d0 )
c
            if(.not.(txy.le.tol.and. txz.le.tol .and.tyz.le.tol)) then
               linear=.false.
               go to 100
            endif
         enddo
      endif
c
      if(abs(b).gt.tol.and.abs(c).gt.tol) then
c               yz
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
            yp1=y-y1  
            zp1=z-z1
c
            if( abs(xp1).gt.tol) then
               linear=.false.
               go to 100
            else
               ty=yp1/b
               tz=zp1/c
               tyz=abs( ty/tz-1.d0 )
               if( tyz.gt.tol) then
                  linear=.false.
                  go to 100
               endif
            endif
         enddo
      endif
c
      if(abs(a).gt.tol.and.abs(c).gt.tol) then
c               xz
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
            yp1=y-y1  
            zp1=z-z1
c
            if( abs(yp1).gt.tol) then
               linear=.false.
               go to 100
            else
               tx=xp1/a
               tz=zp1/c
               txz=abs( tx/tz-1.d0 )
               if( txz.gt.tol) then
                  linear=.false.
                  go to 100
               endif
            endif
         enddo
      endif
c
      if(abs(a).gt.tol.and.abs(b).gt.tol) then
c               xy
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
            yp1=y-y1  
            zp1=z-z1
c
            if( abs(zp1).gt.tol) then
               linear=.false.
               go to 100
            else
               tx=xp1/a
               ty=zp1/b
               txy=abs( tx/ty-1.d0 )
               if( txy.gt.tol) then
                  linear=.false.
                  go to 100
               endif
            endif
         enddo
      endif
c
      if(abs(a).gt.tol) then
c            b and c = zero
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
cccc        xp1=x-x1
            yp1=y-y1  
            zp1=z-z1
            if( abs(yp1).gt.tol.or.abs(zp1).gt.tol) then
               linear=.false.
               go to 100
            endif
         enddo
      endif
      if(abs(b).gt.tol) then
c            a and c = zero
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
ccccc       yp1=y-y1  
            zp1=z-z1
            if( abs(xp1).gt.tol.or.abs(zp1).gt.tol) then
               linear=.false.
               go to 100
            endif
         enddo
      endif
      if(abs(c).gt.tol) then
c            a and b = zero
         do iat=3,nat
            x=datnuc(2,iat)
            y=datnuc(3,iat)
            z=datnuc(4,iat)
            xp1=x-x1
            yp1=y-y1  
ccccc       zp1=z-z1
            if( abs(xp1).gt.tol.or.abs(yp1).gt.tol) then
               linear=.false.
               go to 100
            endif
         enddo
      endif
c
  100 continue
c
      end

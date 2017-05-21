C===========mp2iter============================================
      subroutine mp2iter(ncore,ncf,nval,nstron,ilist,
     1                   jlist,invl,idm,lbas,itst,
     2                   listp,nkunix,npairs,iprint,equi,
     3                   uniq,ifunpair,isorb,nsym,norx,
     4                   clim,erhf,weak,lbasij,idmij,
     5                   npu,inndp,enegl,tcou,redu,
     6                   nnegl,nweak,listw,idmw,lbasw,
     7                   dsijx,coro,remf,jdm,jdmij,
     8                   ifil,krec,munik,nfiles,giga,
     9                   gigabyte)
C              
C       main program for MP2-iterations
C       
C       Svein Saebo, Fayetteville AR Summer 1999
C       and summer and fall 2000
C
C       throughout this program the diagonal pairs are treated first
C       then the non-diagonal pairs and finally the weak pairs.
C       this improves convergence of the dynamical updating
C
C       when you enter this subroutine you should have:
C
C       projected overlap in quadratic form in 'ovl'
C       projected fock matrix in quadratic form in  'fck'
C       projection matrix in 'proj'
C       fock matrix in localized basis quadratic form in'floc'
C       All amplitudes will be stored as a single vector in 'tful'
C
C       The internal exchange matrices in projected local
C       basis are on disk: unit ndisk4,nkunit  (xxxx.kij)
C

      use memory

        implicit real*8(a-h,o-z)
c     common /intbl/maxsh,inx(100)
c     common/big/bl(30000)
       logical*1 nstron(nval,nval)
       logical ndiag,nconv,maxit,econv,econv2,weak,tcou,redu
       logical coro,nrem,remf,ewkc
       dimension isorb(nval,*),ifunpair(7,*),dsijx(nval,*)
       dimension ilist(*),jlist(*),invl(nval,nval),listp(*)
       integer*2 idmw(*),lbasw(*),listw(nval,*),jdm(*),jdmij(*)
       integer equi(2,*)
       integer uniq(*)
       dimension itst(*),ifil(*),krec(*)
       integer*2 lbas(*),idm(*),lbasij(*),idmij(*),inndp(*)
      character*3 ch3
      character*256 scrfile,filname0,filname1
        parameter(zero=0.0d0)
C
       nrem=.not.remf
C         get input/output file names
      inp=igetival('inp')
      iout=igetival('iout')
      icond=igetival('icond')
      idoub=igetival('double')
       tno2=zero
       trij=zero
       tupd=zero
C
       ntunit=norx+nfiles+1
      write(6,*) 'kij: ', nkunix+1
      write(6,*) 'orb: ', norx+1
      write(6,*) 'Tij: ', ntunit
      Write(6,*) 'nfiles=',nfiles
C
       nvir=ncf-nval-ncore
       nmo=nval+ncore
       idmx=0
       do ii=1,nval
       idd=idm(ii)
       idmx=max0(idd,idmx)
       enddo
C
       write(6,*) 'max local domain: ',idmx
       do llp=1,npu
       idd=idmij(llp)
       idmx=max0(idd,idmx)
       enddo
       write(6,*) 'max local domain: ',idmx
       if(.not.coro) then
       lrec=idmx*idmx
       lrec=lrec*idoub
C
C       determine how many files are needed
C       our current linux fortran only allows files up to 2GB
C       limited to 11  files for now, can easlily be extended
C       numbers too big here for integer arithmetic
C
       tost=float(lrec)*float(munik)
       write(6,*) 'munik=',munik,lrec,idoub,tost
        nnfil=tost/gigabyte +1
       write(6,*) nnfil,' files needed for amplitudes',ntunit
C         open nnfil files for amplitudes
       call getchval('scrf',scrfile)
          call rmblan(scrfile,256,len)
       filname0=scrfile(1:len)//'.Tij'
          len1=len+4
      open(unit=ntunit,file=filname0(1:len1),form='unformatted'
     1 ,access='direct',recl=lrec)
          if(nnfil.gt.1) then
          nstripe=nnfil-1
      Do itij=1,NStripe
      write(ch3,'(I3)') itij
      ch3(1:1) = '.'
      If(itij.lt.10) ch3(2:2) = '0'
      KUnit =ntunit+itij
      filname1 = filname0(1:len1)//ch3
      len = len1+3
      OPEN (UNIT=KUnit,FILE=filname1(1:len),FORM='UNFORMATTED',
     $      ACCESS='DIRECT',RECL=lrec)
      EndDo
          endif
C
       if(nnfil.gt.11.or.nnfil.le.0) then
       write(6,*) nnfil,' files needed but only 11 allowed'
       call nerror(1,'mp2iter','too many files',nnfil,itost)
       endif
C
       nopef=munik/nnfil+1
       kfil=ntunit
       nopf=0
       do ipppx=1,munik
       if(ipppx.gt.nopef*(kfil-ntunit+1)) then
       kfil=kfil+1
       nopf=0
       endif
       nopf=nopf+1
       ifil(ipppx)=kfil
       krec(ipppx)=nopf
       enddo
C       write(6,*) (ifil(iw),iw=1,munik)
C       write(6,*) (krec(iw),iw=1,munik)
C
       call matdef('Tmax','q',idmx,idmx)
       itmxa=mataddr('Tmax')
       idmxs=idmx*idmx
       endif
C
       if(nsym.gt.0) then
c      call getmem((nsym*idmx*(nval+npu+nval))/4+1,lbtb)
       call getint_2(nsym*idmx*(nval+npu+nval),lbtb)
      call mkltab(ifunpair,lbas,idm,nval,ncf,
     1            bl(lbtb),idmx,nsym,idmij,lbasij,
     2            lbasw,idmw,weak,nstron,listw,
     3            npairs)
       else
       call getmem(1,lbtb)
       endif
C       determine the total dimension of the amplitudes
       call tsetup(ijfum,nval,nstron,idm,itst,
     1                iprint,uniq,idmij,idmw,nweak,
     2                listw,weak,coro,iextr)
C
      if(iprint.ge.2)write(6,*)'tful',ijfum
C       zero out space for amplitudes
       if(coro) then
       call matdef('tful','v',ijfum,ijfum)
       call matzero('tful')
       else
       if(weak) then
       call matdef('tful','v',iextr+1,iextr+1)
       call matzero('tful')
       endif
       endif
C
       if(weak) then
C        construct  orbitals for weak pairs and leave in memory
C       for the duration of the calculation
       call getint(nval+1,iorst)
       call getint(nval,ipps)
      call estimem(memor,memei,nval,idmw,bl(iorst),
     1             bl(ipps),ncf,ncore)
       call matdef('morbs','d',memor,memor)
       call matdef('meigs','d',memei,memei)
       iora=mataddr('morbs')
       ieia=mataddr('meigs')
      call neworb3(ncf,nval,idmw,lbasw,bl(ipps),
     1             bl(iorst),bl(iora),bl(ieia))
       endif
C
       np4=4
       call matdef('eval','d',nmo,nmo)
       call matread('eval',np4,'eloc_rhf')
       call matdef('x1','q',ncf,ncf)
C
C       calculate average local dimension
       meme=0
       do ixyz=1,nval
       meme=meme+idm(ixyz)
       enddo
       do ixxx=1,npu
       meme=meme+idmij(ixxx)
       enddo
       aver=float(meme)/(float(nval)+float(npu))
       prct=100*(aver/(ncf-nval-ncore))
       write(6,*) 'Average dimension of local domains: ',aver
       write(6,*) prct,' % of full'
                    call f_lush(6)
C
C       rewind files
C
       do iffilx=1,nfiles
       noru=norx+iffilx
       nkunit=nkunix+iffilx
       rewind noru
       rewind nkunit
       enddo
C
       irmt=0
       niter=0
       nconv=.false.
       maxit=.false.
       econv=.false.
       econv2=.false.
       ewkc=.false.
       nitmax=20
C       convergence should be reached in 20 iterations unless
C       something is wrong
       ifocl=mataddr('floc')
       emp2=zero
       eold=zero
       eweako=zero
       dtmax=zero
       tmax=zero
       tcopl=zero
       ttmat=zero
       ieoa=mataddr('eval')-1
C
C       ********end preliminary stuff *******
C
C       internal exchange matrices on :     nkunit
C       temporary MO basis on :             noru
C       (after pass through n'rewind',eworb2)
C
C       generate initial Tijs and leave in memory for the duration
C       of the calculation
C
       if(redu) call matdef('smaq','q',ncf,ncf)
C
C       diagonal pairs first
C
C
       lfil=1
       wrdw=zero
       lpt=0
       lpar=0
       ndiag=.false.
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj))cycle
       lpar=lpar+1
       if(ii.ne.jj) cycle
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       idd=idm(ii)
       iff=ncf*(ii-1)+1
       ijdm=idd*idd
       call matdef('Rij','q',idd,idd)
       ikadr=mataddr('Rij')
C
       if(nrem) then
        call neworb2n(ncf,lbas(iff),idd,jdd,noru,
     1                  tno2,irmt)
       jdm(ii)=jdd
       else
       call neworb2(ncf,lbas(iff),idd,norx,tno2,
     1    redu,lfil,wrdw,giga)
       endif
C
C       read internal exchange operator from xx.kij (unit nkunut)
C       note this file is opened in transloc.
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C       transform to MO basis and get initial T's
       call secund(ttup1)
C
       if(nrem) then
      call updaten(ncore,ncf,nval,nvir,idd,
     1           jdd,ii,jj,dtmax,niter)
       call matrem('orbj')
       else
      call update(ncore,ncf,nval,nvir,idd,
     1            ii,jj,dtmax,niter)
       endif
C
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       if(coro) then
       call moveT(lpar,idd,itst)
C
       else
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call wrida(bl(ikadr),bl(itmxa),ijdm,idmxs,irec,nunit)
       endif
C
       call matrem('Rij')
       end do
       end do
       if(iprint.ge.2)write(6,*) 'dtmax= ', dtmax
C
C       non-diagonal pairs
C
       ndiag=.true.
       lpar=0
       lpt=0
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj))cycle
       lpar=lpar+1
       if(ii.eq.jj) cycle
       lpt=lpt+1
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       idd=idmij(lpt)
       iff=ncf*(lpt-1)+1
       ijdm=idd*idd
       call matdef('Rij','q',idd,idd)
       ikadr=mataddr('Rij')
C
       if(nrem) then
        call neworb2n(ncf,lbasij(iff),idd,jdd,noru,
     1                  tno2,irmt)
       jdmij(lpt)=jdd
       else
       call neworb2(ncf,lbasij(iff),idd,norx,tno2,
     1    redu,lfil,wrdw,giga)
       endif
C
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C       transform to MO basis and get initial T's
       call secund(ttup1)
C
       if(nrem) then
      call updaten(ncore,ncf,nval,nvir,idd,
     1           jdd,ii,jj,dtmax,niter)
       call matrem('orbj')
       else
      call update(ncore,ncf,nval,nvir,idd,
     1            ii,jj,dtmax,niter)
       endif
C
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       if(coro) then
       call moveT(lpar,idd,itst)
       else
       nunit=ifil(iuniq)
       irec=krec(iuniq)
       
C
      call wrida(bl(ikadr),bl(itmxa),ijdm,idmxs,irec,nunit)
       endif
C
       call matrem('Rij')
       end do
       end do
       if(iprint.ge.2) then
       write(6,*) 'dtmax= ', dtmax
      write(6,*) 'A total of ',irmt,'redundant basis functions ignored'
       endif
C
C       weak  pairs
C
       if(weak) then
C
       dtmaxw=zero
       do ii=1,nval
       do jj=1,ii
       nwea=listw(ii,jj)
       if(nwea.le.0) cycle
       lpar=npairs+nwea
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       if(coro) then
       lppr=lpar
       else
       lppr=nwea
       endif
       idd=idmw(ii)
       jdd=idmw(jj)
       ijdm=idd*jdd
       call matdef('Rij','r',idd,jdd)
       ikadr=mataddr('Rij')
C
       wrdw=wrdw+idd*jdd
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=idd*jdd
       endif
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C
C       transform to MO basis and get initial T's
       call secund(ttup1)
       call updatw(ncore,ncf,lbasw,ii,jj,
     1                idd,jdd,dtmaxw,bl(iorst),bl(ipps))
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       call moveTw(lppr,ijdm,itst)
       call matrem('Rij')
       end do
       end do
       if(iprint.ge.2) then
       write(6,*) 'dtmaxw= ', dtmaxw,lpar
       endif
       endif
       if(redu)call matrem('smaq')
       ttno2=tno2/60.0d0
       write(6,*) 'Construction of MO-basis: ',ttno2,' min.'
C       Start MP2 iterations...
       write(6,*) 'Start MP2 - iterations '
       if(weak) then
      write(6,*) ' ITER         ESTRONG            EWEAK       DELTAT
     1     DELTATW'
       else
      write(6,*) ' ITER          EMP2         DELTAT '
       endif
      call f_lush(6)
       do niter=1,nitmax
       if(niter.eq.(nitmax)) maxit=.true.
       emp2=zero
C       first calculate correlation energy form old C's
       if(nconv.or.maxit.or.econv) then
       write(6,*) '      I    J    PAIR-ENERGY      MULTIPLICITY'
       endif
C
C       rewind files
C
       do iffilx=1,nfiles
       noru=norx+iffilx
       nkunit=nkunix+iffilx
       rewind noru
       rewind nkunit
       enddo
C
C
C       diagonal pairs
C
       lpar=0
       ldia=0
       lpt=0
       disij=zero
       ndiag=.false.
       lfil=1
       wrdw=zero
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
C       diagonal pairs
       if(ii.ne.jj) cycle
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       multi=equi(1,lpar)
       idd=idm(ii)
       ijdm=idd*idd
       call matdef('Rij','q',idd,idd)
C       read in Kij and put into Rij
       ikadr=mataddr('Rij')
C
       wrdw=wrdw+ijdm+idd
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=ijdm+idd
       endif
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C
C       calculate pair energy for this pair
       if(.not.coro) then
       nunit=ifil(iuniq)
       irec=krec(iuniq)
       endif
      call pener(lpar,ii,jj,eij,idd,
     1           itst,irec,coro,nunit,idmxs)
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       write(6,500) ii,jj,eij,multi,disij
  500 format(1x,2I6,1f13.8,6x,'         ',1I3,2x,1f7.4)
       endif
       ldia=ldia+multi
       emp2=emp2+eij*multi
       call matrem('Rij')
       end do
       end do
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       write(6,*)ldia,' diagonal pairs'
       write(6,501)zero,emp2
       endif
                    call f_lush(6)
C
C       Non-diagonal pairs
C
       dsav=zero
       ndiag=.true.
       emps=zero
       lpar=0
       lpt=0
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
C       nondiagonal pairs
       if(ii.eq.jj) cycle
      lpt=lpt+1
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       multi=equi(1,lpar)
       idd=idmij(lpt)
       ijdm=idd*idd
       call matdef('Rij','q',idd,idd)
C       read in Kij and put into Rij
       ikadr=mataddr('Rij')
       wrdw=wrdw+ijdm+idd
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=ijdm+idd
       endif
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C       calculate pair energy for this pair
       if(.not.coro) then
       nunit=ifil(iuniq)
       irec=krec(iuniq)
       endif
      call pener(lpar,ii,jj,eij,idd,
     1           itst,irec,coro,nunit,idmxs)
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       disij=0.529167d0/dsijx(ii,jj)
       dsav=dsav+disij*multi
       write(6,500) ii,jj,eij,multi,disij
       endif
       emps=emps+eij*multi
       call matrem('Rij')
       end do
       end do
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       dsav=dsav/float(npu)
       write(6,*) npu, ' nondiagonal strong pairs'
       write(6,501) dsav,emps
  501 format(1x,'Average distance: ',1f7.4, ' energy: ',1f12.7)
       endif
       emp2=emp2+emps
                    call f_lush(6)
C
C       weak pairs
C
       dsav=zero
       eweak=zero
       if(weak) then
       do ii=1,nval
       do jj=1,ii
       nwea=listw(ii,jj)
       if(nwea.le.0) cycle
       lpar=npairs+nwea
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       if(coro) then
       llpr=lpar
       else
       llpr=nwea
       endif
       multi=equi(1,lpar)
       idd=idmw(ii)
       jdd=idmw(jj)
       ijdm=idd*jdd
C       read in Kij and put into Rij
       call matdef('Rij','r',idd,jdd)
       ikadr=mataddr('Rij')
       wrdw=wrdw+ijdm
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=ijdm
       endif
       nkunit=nkunix+lfil
       call reaarr(bl(ikadr),ijdm,nkunit)
C       calculate pair energy for this pair
      call penw(llpr,ii,jj,eij,idd,
     1          jdd,itst)
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       disij=0.529167d0/dsijx(ii,jj)
       dsav=dsav+disij*multi
       write(6,500) ii,jj,eij,multi,disij
       endif
       eweak=eweak+eij*multi
       call matrem('Rij')
       end do
       end do
       weakd=abs(eweak-eweako)
       if(weakd.le.5.0d-06) ewkc=.true.
        if(econv.or.nconv.or.maxit.or.iprint.ge.3)then
       dsav=dsav/float(nweak)
       write(6,*) nweak,' weak pairs'
       write(6,501) dsav,eweak
       endif
       emp2=emp2+eweak
       endif
       if(econv.and.abs(eold-emp2).lt.1.0d-6)econv2=.true.
       if(abs(eold-emp2).lt.1.0d-6)econv=.true.
       eold=emp2
       emp2=emp2-eweak
       if(weak) then
       write(6,50) niter,emp2,eweak,dtmax,dtmaxw
       else
       write(6,51) niter,emp2,dtmax
       endif
   50   format(1x,i4,5x,1f12.7,5x,1f12.7,5x,1f9.7,5x,1f9.7)
   51   format(1x,i4,5x,1f12.7,5x,1f9.7)
          call f_lush(6)
       if(econv2.or.nconv.or.maxit) then
C       print final results:
C     erhf=rgetrval('erhf')
      write(iout,*) '-------------------------------------------------'
       if(weak) then
      write(iout,*) '                Final SCF/MP2 results '
       else
      write(iout,*) '           SCF/MP2 results for strong pairs'
       endif
      write(iout,*) ' '
      write(iout,3001) erhf
      write(iout,*) '-------------------------------------------------'
       write(6,*) npairs, ' Strong Pairs:'
      write(iout,3003) emp2
      write(icond,3003) emp2
       if(weak)  then
       write(6,*) nweak,' Weak Pairs:'
       write(iout,3005) eweak
       write(icond,3005) eweak
       write(6,*) nnegl,' Neglected Pairs:'
       write(iout,3006) enegl
       write(icond,3006) enegl
       else
       write(6,*) nweak+nnegl,' Neglected Pairs:'
       write(iout,3006) enegl
       write(icond,3006) enegl
       endif
      write(iout,*) '-------------------------------------------------'
       emp2=emp2+enegl+eweak
      write(iout,3003) emp2
      write(icond,3003) emp2
      write(iout,*) '-------------------------------------------------'
      write(iout,3004) erhf + emp2
      write(icond,3004) erhf + emp2
      write(iout,*) '-------------------------------------------------'
 3005 format('        Weak pairs:              =',f15.8)
 3006 format('        Neglected Pairs (est.):  =',f15.8)
 3001 format('        Total SCF Energy         =',f15.8)
 3003 format('        MP2 correlation Energy   =',f15.8)
 3002 format('        Weak pairs (estimated)   =',f15.8)
 3004 format('        Total MP2 Energy         =',f15.8)
       if(econv2.and..not.nconv) then
      write(6,*) 'WARNING: The correlation energy seem to be converged'
      write(6,*) 'but not the amplitudes'
       endif
       exit
       endif
C       loop over pairs, construct second order residuum + update
C       the coefficients for one pair at the time.
C
C       rewind files
C
          call f_lush(6)
       do iffilx=1,nfiles
       nkunit=nkunix+iffilx
       rewind nkunit
       noru=norx+iffilx
       rewind noru
       enddo
C
       lpar=0
       lpt=0
       ndiag=.false.
       dtmax=zero
       dtmaxw=zero
       lfil=1
       wrdw=zero
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
C       first diagonal pairs
       if(ii.ne.jj) cycle
       iuniq=uniq(lpar)
       if(iuniq.le.0)cycle
       idd=idm(ii)
       iff=ncf*(ii-1)+1
      call RijMP2(lpar,ii,jj,idm,lbas,
     1            ncf,nval,bl(ifocl),nstron,nconv,
     2            invl,listp,nkunix,itst,ifunpair,
     3            isorb,uniq,equi,nsym,ilist,
     4            jlist,trij,bl(lbtb),idmx,idmij,
     5            lbasij,lpt,ndiag,inndp,idd,
     6            iff,tcopl,tcou,ttmat,redu,
     7            idmw,lbasw,weak,listw,npairs,
     8            coro,krec,ifil,niter,lfil,
     9            wrdw,giga)
C       update amplitudes for this pair
       call secund(ttup1)
C
       if(nrem) then
       jdd=jdm(ii)
      call getorbn(idd,jdd,noru)
      call updaten(ncore,ncf,nval,nvir,idd,
     1           jdd,ii,jj,dtmax,niter)
       else
       call getorb(idd,norx,lfil)
      call update(ncore,ncf,nval,nvir,idd,
     1            ii,jj,dtmax,niter)
       endif
C
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       if(coro) then
       call moveT(lpar,idd,itst)
       else
C
       call matdef('Tij','q',idd,idd)
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),idd*idd,idmxs,irec,nunit)
       call matadd('Rij','Tij')
      call wrida(bl(itij),bl(itmxa),idd*idd,idmxs,irec,nunit)
       call matrem('Tij')
       endif
C
       call matrem('Rij')
       enddo
       enddo
       lpar=0
       lpt=0
       ndiag=.true.
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
       if(ii.eq.jj) cycle
       lpt=lpt+1
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       idd=idmij(lpt)
       iff=ncf*(lpt-1)+1
      call RijMP2(lpar,ii,jj,idm,lbas,
     1            ncf,nval,bl(ifocl),nstron,nconv,
     2            invl,listp,nkunix,itst,ifunpair,
     3            isorb,uniq,equi,nsym,ilist,
     4            jlist,trij,bl(lbtb),idmx,idmij,
     5            lbasij,lpt,ndiag,inndp,idd,
     6            iff,tcopl,tcou,ttmat,redu,
     7            idmw,lbasw,weak,listw,npairs,
     8            coro,krec,ifil,niter,lfil,
     9            wrdw,giga)
C       update amplitudes for this pair
       call secund(ttup1)
C
       if(nrem) then
       jdd=jdmij(lpt)
      call getorbn(idd,jdd,noru)
      call updaten(ncore,ncf,nval,nvir,idd,
     1           jdd,ii,jj,dtmax,niter)
       else
       call getorb(idd,norx,lfil)
      call update(ncore,ncf,nval,nvir,idd,
     1            ii,jj,dtmax,niter)
       endif
C
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       if(coro) then
       call moveT(lpar,idd,itst)
       else
C
       call matdef('Tij','q',idd,idd)
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),idd*idd,idmxs,irec,nunit)
       call matadd('Rij','Tij')
      call wrida(bl(itij),bl(itmxa),idd*idd,idmxs,irec,nunit)
       call matrem('Tij')
       endif
C
       call matrem('Rij')
       end do     !  loop  over  ii
       end do
C
C       weak pairs
C              
       if(weak.and..not.ewkc) then
       do ii=1,nval
       do jj=1,ii
       nwea=listw(ii,jj)
       if(nwea.le.0) cycle
       lpar=npairs+nwea
       iuniq=uniq(lpar)
       if(iuniq.le.0) cycle
       if(coro) then
       lppw=lpar
       else
       lppw=nwea
       endif
       idd=idmw(ii)
       jdd=idmw(jj)
       iff=ncf*(ii-1)+1
       jff=ncf*(jj-1)+1
      call RijMP2w(lppw,ii,jj,idm,lbas,
     1             ncf,nval,bl(ifocl),nstron,nconv,
     2             invl,listp,nkunix,itst,ifunpair,
     3             isorb,uniq,equi,nsym,ilist,
     4             jlist,trij,bl(lbtb),idmx,idmij,
     5             lbasij,lpt,ndiag,inndp,idd,
     6             jdd,iff,jff,tcopl,ttmat,
     7             redu,idmw,lbasw,listw,npairs,
     8             coro,krec,ifil,lfil,wrdw,
     9             giga)
C       update amplitudes for this pair
       call secund(ttup1)
       call updatw(ncore,ncf,lbasw,ii,jj,
     1                idd,jdd,dtmaxw,bl(iorst),bl(ipps))
       call secund(ttup2)
       tupd=tupd+ttup2-ttup1
       call moveTw(lppw,idd*jdd,itst)
C
       ikadr=mataddr('Rij')
C
       call matrem('Rij')
       end do     !  loop  over  ii
       end do
       endif
C
       if(dtmax.lt.clim.and.dtmaxw.lt.clim) nconv=.true.
       eweako=eweak
       end do     !  mp2 terations
       if(nconv.or.econv2) then
        write(6,*) 'Convergence after ', niter, 'iterations'
       else
       write(6,*) 'No Convergence EMP2=',emp2
       endif
C
C       close and delete scratch files
       do iiunit=1,nnfil
       nunit=ntunit+iiunit-1
       close(unit=nunit,status='delete')
       enddo
C
       if(iprint.ge.2) then
       ttno2=tno2/60.00d0
       write(6,*) 'time generating orbitals: ',ttno2,' min.'
       ttcop=tcopl/60.0d0
       write(6,*) 'calculation of coupling terms:',ttcop,' min.'
       tmat=ttmat/60.0d0
       write(6,*) 'matrix multiplications', tmat
       ttrij=trij/60.00d0
       write(6,*) 'Generation of MP2-residua: ',ttrij,' min.'
       ttupd=tupd/60.00d0
       write(6,*) 'Updating of amplitudes: ', ttupd,' min.'
       endif
C       This completes the calculation
C       Note all scratch files are closed in newtransloc
       end
C================rijmp2============================================
      subroutine RijMP2(lpar,ii,jj,idm,lbas,
     1                  ncf,nval,FF,nstron,nconv,
     2                  invl,listp,nkunix,itst,ifunpair,
     3                  isorb,uniq,equi,nsym,ilist,
     4                  jlist,trij,lbtab,idmx,idmij,
     5                  lbasij,lpt,ndiag,inndp,idd,
     6                  iff,tcopl,tcou,ttmat,redu,
     7                  idmw,lbasw,weak,listw,npairs,
     8                  coro,krec,ifil,niter,lfil,
     9                  wrdw,giga)
C
C       Calculates MP2 residuum using generator state formalism
C
C       Calculates Second order residuum for one pair.
C       lpar pair number
C       ii, jj orbital indices
C       idm(*) dimension of local domains
C       lbas(ncf,*)  contains the local bases
C       ncf  number of contraced basis functions
C       nval number of valence orbitals
C       FF FOck matrix in localized orbitals
C       ilist and jlist give  i and j for a given pair
C       'X1' used for temporary storage ncf*ncf defined in the
C       calling program
C
C       The final residuum will be in 'Rij' in local dimension
C
C       Svein Saebo, Fayetteville, AR Summer 1999
C
C       called from mp2iter
C

      use memory

       implicit real*8(a-h,o-z)
       dimension ifunpair(7,*),isorb(nval,*)
       integer equi(2,*)
       integer uniq(*)
       integer*2 idmw(*),lbasw(*)
       dimension ilist(*),jlist(*)
       dimension krec(*),ifil(*)
       logical coro
       logical*1 nstron(nval,nval)
       logical nconv,ndiag,tcou,redu,weak
       dimension invl(nval,*),listp(*)
       dimension FF(nval,nval)
       integer*2 lbas(*),idm(*),idmij(*),lbasij(*)
       integer*2 lbtab(idmx,nsym,*),listw(nval,*)
       integer*2 inndp(*)
       integer itst(*)
c      common/big/bl(30000)
C
       itmxa=mataddr('Tmax')
       idmxs=idmx*idmx
C
       call secund(t1rij)
       ifoca=mataddr('fck')
       iovad=mataddr('ovl')
       ijdim=idd*idd
       call matdef('Rij','q',idd,idd)
       irad=mataddr('Rij')
       call matdef('Tij','q',idd,idd)
C
       wrdw=wrdw+idd*idd+idd
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=idd*idd+idd
       endif
       nkunit=nkunix+lfil
       call reaarr(bl(irad),idd*idd,nkunit)
C
       itij=mataddr('Tij')
       if(coro) then
       call getT(lpar,idd,itst)
       else
C
       nunit=ifil(uniq(lpar))
       irec=krec(uniq(lpar))
      call reada(bl(itij),bl(itmxa),ijdim,idmxs,irec,nunit)
       endif
       if(redu) then
       call penalty(lpar,ncf,idd,iff,lbas,
     1                 itst,lpt,lbasij,ndiag)
       if(niter.le.5)then
       callmatrem('Tij')
       return
       endif
       endif
C
C       first calculate FCS and SCF terms
       call matdef('FCS','q',idd,idd)
       ifcsa=mataddr('FCS')
C       FCijS-term
C    calculate F*Cij --> Y1
       call matdef('FCij','q',idd,idd)
       call matdef('Fii','q',idd,idd)
       ifiia=mataddr('Fii')
       if(ndiag) then
       call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbasij(iff),lbasij(iff))
       else
       call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbas(iff),lbas(iff))
       endif
       call matmmult('Fii','Tij','FCij')
       call matrem('Fii')
C    multiply with S from the right and add to Rij
       call matdef('Sjj','q',idd,idd)
       isjja=mataddr('Sjj')
       if(ndiag) then
       call compab(bl(iovad),bl(isjja),ncf,idd,idd,
     1             lbasij(iff),lbasij(iff))
       else
       call compab(bl(iovad),bl(isjja),ncf,idd,idd,
     1             lbas(iff),lbas(iff))
       endif
       call matmmult('FCij','Sjj','FCS')
       call matrem('Sjj')
       call matadd('FCS','Rij')
C       SCijF-terms
C       calculate SCij --> FCij
       call matdef('Sii','q',idd,idd)
       isiia=mataddr('Sii')
       if(ndiag) then
       call compab(bl(iovad),bl(isiia),ncf,idd,idd,
     1             lbasij(iff),lbasij(iff))
       else
       call compab(bl(iovad),bl(isiia),ncf,idd,idd,
     1             lbas(iff),lbas(iff))
       endif
       call matmmult('Sii','Tij','FCij')
       call matrem('Sii')
C    multiply with F from the right and add to Rij
       call matdef('Fjj','q',idd,idd)
       ifjja=mataddr('Fjj')
       if(ndiag) then
       call compab(bl(ifoca),bl(ifjja),ncf,idd,idd,
     1             lbasij(iff),lbasij(iff))
       else
       call compab(bl(ifoca),bl(ifjja),ncf,idd,idd,
     1             lbas(iff),lbas(iff))
       endif
       call matmmult('FCij','Fjj','FCS')
       call matrem('Fjj')
       call matrem('FCij')
       call matadd('FCS','Rij')
       call matrem('FCS')
       call matrem('Tij')
C       now calculate coupling terms
C       4 subroutines:
C       TcouplT  multiplications with S outside sum over k(default)
C       Tcoupl2  multiplications with S inside sum over k
C
C       tcoupls is the same as tcouplt if no symmetry (default)
C       tcoupl1 is the same as tcoupl2 if no symmetry
C
C       tcoupls and tcoupl1 can probabbly be eliminated,
C       introduced mainly
C       for debugging purposes but some logic and table
C       lookup eliminated
C
       call secund(ttco1)
       call matdef('FCS','q',idd,idd)
       if(nsym.le.0) then
       if(tcou) then
C        no symmetry multiply by S inside loop over k
      call Tcoupl1(ii,jj,nval,nstron,idm,
     1             invl,listp,FF,lbas,ncf,
     2             itst,ilist,jlist,idmij,lbasij,
     3             lpt,ndiag,inndp,idd,iff,
     4             ttmat,lbasw,weak,idmw,listw,
     5             npairs,coro,idmxs,ifil,krec)
       else
c      call getmem(ncf/4+1,lbun)
       call getint_2(ncf,lbun)
C        no symmetry multiply by S ouside loop over k
      call Tcoupls(ii,jj,nval,nstron,idm,
     1             invl,listp,FF,lbas,ncf,
     2             itst,ilist,jlist,idmij,lbasij,
     3             lpt,ndiag,inndp,idd,iff,
     4             ttmat,bl(lbun),lbasw,weak,
     5             idmw,listw,npairs,coro,idmxs,
     6             ifil,krec)
       call retmem(1)
       endif
       else
       if(tcou) then
C       symmetry  used, multiply by S inside loop over k
      call Tcoupl2(ii,jj,nval,nstron,idm,
     1             invl,listp,FF,lbas,ncf,
     2             itst,ifunpair,isorb,uniq,equi,
     3             nsym,ilist,jlist,lbtab,idmx,
     4             idmij,lbasij,lpt,ndiag,inndp,
     5                idd,iff,ttmat,lbasw,weak,
     6             idmw,listw,npairs,coro,ifil,
     7             krec)
       else
c      call getmem(ncf/4+1,lbun)
       call getint_2(ncf,lbun)
C       symmetry  used, multiply by S ouside loop over k
      call Tcouplt(ii,jj,nval,nstron,idm,
     1             invl,listp,FF,lbas,ncf,
     2             itst,ifunpair,isorb,uniq,equi,
     3             nsym,ilist,jlist,lbtab,idmx,
     4             idmij,lbasij,lpt,ndiag,inndp,
     5                idd,iff,ttmat,bl(lbun),lbasw,
     6             weak,idmw,listw,npairs,coro,
     7             ifil,krec)
       call retmem(1)
       endif
       endif
       call secund(ttco2)
       tcopl=tcopl+ttco2-ttco1
       call matrem('FCS')
C
C
       call secund(t2rij)
       trij=trij+t2rij-t1rij
        end
C==============update=========================================
      subroutine update(ncore,ncf,nval,nvir,idd,
     1                  ii,jj,dtmax,niter)
C
C       this subroutine transforms the mp2 residuums to a temporary
C       MO basis and updates the amplitudes.
C       default:
C       The MO basis was  constructed in either in neworb2
C       and kept on disk noru or
C
C       Svein Saebo Fayetteville AR June 1999
C
C       and modified in Starkville, MS 1999,and 2000
C
C       called from mp2iter
C       this subroutine was called neworb in older verions
C

      use memory

       implicit real*8(a-h,o-z)
c      common/big/bl(30000)
C
C       R for this pair is in 'Rij' T for this pair in 'tful'
C       position itst(lpar)
C
C       orbitals in orbi
C       and the eigenvalues in 'eigi'
C       normally read from noru in subroutine neworb2
       iega=mataddr('eigi')
       call matdef('Rmo','q',idd,idd)
C       
C       transform residuum to temporary MO basis
C
       call traMO(idd)
C
C       the Rij in MO basis is now in  'Rmo'
       irad=mataddr('Rmo')
C       calculate deltat in MO basis
       ieoa=mataddr('eval')-1
       EO=-bl(ieoa+ii+ncore)-bl(ieoa+jj+ncore)
C       level shift removed not needed
       call DeltaT(bl(irad),bl(irad),idd,bl(iega),EO,
     1                tmax)
       dtmax=max(dtmax,tmax)
C       transform deltaT back to AO basis
       call traao(idd)
       call matrem('Rmo')
       call matrem('eigi')
       call matrem('orbi')
       end
C==============updaten========================================
      subroutine updaten(ncore,ncf,nval,nvir,idd,
     1                 jdd,ii,jj,dtmax,niter)
C
C       this subroutine transforms the mp2 residuums to a temporary
C       MO basis and updates the amplitudes.
C       default:
C       The MO basis was  constructed in either in neworb2
C       and kept on disk noru or
C
C       Svein Saebo Fayetteville AR June 1999
C
C       and modified in Starkville, MS 1999,and 2000
C
C       called from mp2iter
C       this subroutine was called neworb in older verions
C

      use memory

       implicit real*8(a-h,o-z)
c      common/big/bl(30000)
C
C       R for this pair is in 'Rij' T for this pair in 'tful'
C       position itst(lpar)
C
C       orbitals in orbi
C       and the eigenvalues in 'eigi'
C       normally read from noru in subroutine neworb2
       iega=mataddr('eigi')
       call matdef('Rmo','q',jdd,jdd)
C       
C       transform residuum to temporary MO basis
C
       call traMOn(idd,jdd)
C
C       the Rij in MO basis is now in  'Rmo'
       irad=mataddr('Rmo')
C       calculate deltat in MO basis
       ieoa=mataddr('eval')-1
       EO=-bl(ieoa+ii+ncore)-bl(ieoa+jj+ncore)
C       level shift removed not needed
       call DeltaT(bl(irad),bl(irad),jdd,bl(iega),EO,
     1                tmax)
       dtmax=max(dtmax,tmax)
C       transform deltaT back to AO basis
       call traaon(idd,jdd)
       call matrem('Rmo')
       call matrem('eigi')
       call matrem('orbi')
       end
C================neworb2===========================================
      subroutine neworb2(ncf,lbas,idd,norx,tno2,
     1                   redu,lfil,wrdw,giga)
C
C       this subroutine determines the temporary MO basis for
C       one pair
C       lbas is the array containg the local bases
C       the orbitals will be in  'orbi' and 'orbj' and the
C       eigenvalues in 'eigi' and 'eigj'
C
C       Svein Saebo Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
       integer*2 lbas(*)
       logical redu
c      common/big/bl(30000)
C
       call secund(ttno2)
       iovad=mataddr('ovl')
       ifoca=mataddr('fck')
       call matdef('orbi','q',idd,idd)
       call matdef('eigi','v',idd,idd)
       iorba=mataddr('orbi')
       ieiga=mataddr('eigi')
C
       if(redu) then
       call matdef('sma','s',ncf,ncf)                     
       call matread('sma',1,'s matrix')
       call matcopy('sma','smaq')
       call matrem('sma')
       iovad=mataddr('smaq')
       endif
C
       call matdef('Sii','q',idd,idd)
       isiia=mataddr('Sii')
       call matdef('Fii','q',idd,idd)
       ifiia=mataddr('Fii')
       call matdef('SFii','q',idd,idd)
C       construct S**-1/2
       call compab(bl(iovad),bl(isiia),ncf,idd,idd,
     1             lbas,lbas)
       call matdef('Siis','s',idd,idd)
       call matcopy('Sii','Siis')
       isiis=mataddr('Siis')
      call invsq(bl(isiis),bl(isiia),bl(iorba),bl(ieiga),idd)
       call matrem('Siis')
C       form S**-1/2 * F
       call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbas,lbas)
       call matmmult('Sii','Fii','SFii')
       call matmmult('SFii','Sii','Fii')
       call matrem('SFii')
C......diagonolize Fii to find i part of orbitals and eigenvalues
C********  replace sdiag2 with dspevx **********
       call matdef('fsym','s',idd,idd)
       call matcopy('Fii','fsym')
       ifsym=mataddr('fsym')
       call eig1(bl(ifsym),bl(ifiia),bl(ieiga),idd,idd)
C       call sdiag2(idd,idd,bl(ifiia),bl(ieiga),bl(ifiia))
       call matrem('fsym')
C********  replace sdiag2 with dspevx **********
       call matmmult('Sii','Fii','orbi')
       call matrem('Fii')
       call matrem('Sii')
       idsq=idd**2
       wrdw=wrdw+idsq+idd
       if(wrdw.gt.giga) then
       lfil=lfil+1
       wrdw=idsq+idd
       endif
       noru=norx+lfil
       call wriarr(bl(iorba),idsq,noru)
       call wriarr(bl(ieiga),idd,noru)
       call secund(t2no2)
       tno2=tno2+t2no2-ttno2
       end
C=========deltaT==============================================
       subroutine DeltaT(R,R2,idd,egi,EO,tmax)
C       calculates DeltaR in MO basis
C       R in R (R2 the same as R)
C       EO= -eps(i)-eps(j)
C       egi orbital energies
C
C       Svein Saebo Fayetteville AR summer 1999
C
       implicit real*8(a-h,o-z)
       dimension R(idd,idd),egi(*)
       dimension R2(*)
       parameter(zero=0.0d0)
       tmax=zero
       do jj=1,idd
       del1=eo+egi(jj)
       do ii=1,idd
        del=del1+egi(ii)
       R(ii,jj)=-R(ii,jj)/del
       end do
       end do
C       find the absolute largest element in deltat
       imaxa=idamax(idd*idd,R2,1)
       tmax=abs(R2(imaxa))
       end
C================neworb2n==========================================
      subroutine neworb2n(ncf,lbas,idd,jdd,noru,
     1                    tno2,irmt)
C
C       this subroutine determines the temporary MO basis for
C       one pair
C       lbas is the array containg the local bases
C       the orbitals will be in  'orbi' and 'orbj' and the
C       eigenvalues in 'eigi' and 'eigj'
C
C       Svein Saebo Fayetteville AR June 1999
C
C       this subroutine checks a local domain for redundant functions
C       idd  dimension of local domain
C       lbas local domain
C       ncf number of basis functions
C       epsi  threshold for removing redundant functions
C       The eigenvectors to S are now renormalized as described in
C       J. Chem. Phys. 104 (1996) 6286
C       by Hempel and Werner, we eliminate repeated matrix
C       digonolizations
C
C       Svein Saebo Fayetteville, AR November 2000
C       

      use memory

       implicit real*8(a-h,o-z)
       integer*2 lbas(*)
c      common/big/bl(30000)
C
       Parameter(zero=0.0d0,one=1.0d0)
C
       call secund(ttno2)
       epsi=1.0d-06
C
       iovl=mataddr('ovl')
        call matdef('orbj','q',idd,idd)
        iorbi=mataddr('orbj')
        call matdef('Oiis','s',idd,idd)
        call matdef('Oii','q',idd,idd)
        isiia=mataddr('Oii')
        call compab(bl(iovl),bl(isiia),ncf,idd,idd,lbas,lbas)
        call matdef('eigi','d',idd,idd)
        ieigi=mataddr('eigi')
        call matcopy('Oii','Oiis')
        iossa=mataddr('Oiis')
c       call sdiag2(idd,idd,bl(isiia),bl(ieigi),bl(iorbi+1))
       call eig1(bl(iossa),bl(iorbi),bl(ieigi),idd,idd)
       irem=0
       icol=iorbi-1-idd
       do ibbb=1,idd
       icol=icol+idd
       if(bl(ieigi-1+ibbb).lt.epsi)then
       irem=irem+1
       cycle
       endif
       ffakt=one/sqrt(bl(ieigi-1+ibbb))
       do iaaa=1,idd
      bl(icol+iaaa)=bl(icol+iaaa)*ffakt
       enddo
       enddo
       jdd=idd-irem
       call matrem('eigi')
       call matrem('Oii')
        call matrem('Oiis')
       call matsub('orbi','orbj',irem+1,idd)
       irmt=irmt+irem
C
       iovad=mataddr('ovl')
       ifoca=mataddr('fck')
       call matdef('eigi','v',jdd,jdd)
       ieiga=mataddr('eigi')
C
       call matdef('XFX','q',jdd,jdd)
       call matdef('XF','r',idd,jdd)
       call matdef('Fii','q',idd,idd)
       ifiia=mataddr('Fii')
       call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbas,lbas)
       call matmmult('Fii','orbi','XF')
       call matrem('Fii')
       call matdef('orbx','r',jdd,idd)
       call matpose2('orbi','orbx','n')
       call matmmult('orbx','XF','XFX')
       call matrem('orbx')
       call matrem('XF')
C......diagonolize Fii to find i part of orbitals and eigenvalues
C********  replace sdiag2 with dspevx **********
       call matdef('fsym','s',jdd,jdd)
       call matcopy('XFX','fsym')
       call matdef('Fjj','q',jdd,jdd)
       ifjja=mataddr('Fjj')
       ifsym=mataddr('fsym')
       call eig1(bl(ifsym),bl(ifjja),bl(ieiga),jdd,jdd)
C       call sdiag2(idd,idd,bl(ifiia),bl(ieiga),bl(ifiia))
C********  replace sdiag2 with dspevx **********
       call matdef('orbx','r',idd,jdd)
       call matmmult('orbi','Fjj','orbx')
       call matcopy('orbx','orbi')
       call matrem('orbx')
       call matrem('Fjj')
       call matrem('fsym')
       call matrem('XFX')
       iorba=mataddr('orbi')
       jdsq=idd*jdd
       call wriarr(bl(iorba),jdsq,noru)
       call wriarr(bl(ieiga),jdd,noru)
       call secund(t2no2)
       tno2=tno2+t2no2-ttno2
       end
C=========traMO================================================
       subroutine tramo(idd)
C       transforms a residuum matrix to temporary MO basis
C       the residuum must be in 'Rij'
C
C       Svein Saebo Fayetteville AR summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
c       common/big/bl(30000)
       call matpose('orbi')
       call matdef('tmp1','q',idd,idd)
       call matmmult('orbi','Rij','tmp1')
       call matpose('orbi')
       call matmmult('tmp1','orbi','Rmo')
C
        call matrem('tmp1')
       end
C=========traMOn================================================
       subroutine tramon(idd,jdd)
C       transforms a residuum matrix to temporary MO basis
C       the residuum must be in 'Rij'
C
C       Svein Saebo Fayetteville AR summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
c       common/big/bl(30000)
       call matdef('orbt','r',jdd,idd)
       call matpose2('orbi','orbt','n')
       call matdef('tmp1','r',jdd,idd)
       call matmmult('orbt','Rij','tmp1')
       call matmmult('tmp1','orbi','Rmo')
C
        call matrem('tmp1')
       call matrem('orbt')
       end
C=========traAO===============================================
       subroutine traAO(idd)
C       transforms R back to Ao basis
C       R in MO basis in 'Rmo'
C       R in AO basis put in 'Rij'
C
C       Svein Saebo, Fayetteville Ar summer 1999
C
       implicit real*8(a-h,o-z)
       call matdef('tmp1','q',idd,idd)
       call matmmult('orbi','Rmo','tmp1')
       call matpose('orbi')
       call matmmult('tmp1','orbi','Rij')
       call matrem('tmp1')
       end
C=========traAOn===============================================
       subroutine traAOn(idd,jdd)
C       transforms R back to Ao basis
C       R in MO basis in 'Rmo'
C       R in AO basis put in 'Rij'
C
C       Svein Saebo, Fayetteville Ar summer 1999
C
       implicit real*8(a-h,o-z)
       call matdef('tmp1','r',idd,jdd)
       call matmmult('orbi','Rmo','tmp1')
       call matdef('orbt','r',jdd,idd)
       call matpose2('orbi','orbt','n')
       call matmmult('tmp1','orbt','Rij')
       call matrem('orbt')
       call matrem('tmp1')
       end
C============================================================
      SUBROUTINE SCATAD(A,B,NCF,IDIM,JDIM,
     1                  IFUN,JFUN)
C     **** A IS A FULL MATRIX OF ORDER NCF. B IS A REDUCED
C     **** MATRIX OF DIMENSION IDIM X JDIM.
C     **** THIS SUBROUTINE EXPANDS B TO A FULL MATRIX AND
C     **** ADDS IT TO WHAT IS ALLREADY IN A.
C     **** THE REDUCED BASIS SET ARE IN IFUN(IDIM) AND JFUN(JDIM)
C
C       Svein Saebo Fayetteville Ar, 1999
C       from the original lci program slightly modified for pqs
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NCF,NCF),B(IDIM,JDIM)
      integer*2 IFUN(*),JFUN(*)
      DO  JJ=1,JDIM
       jff=jfun(jj)
      DO  II=1,IDIM
      A(IFUN(ii),jff)=A(IFUN(ii),jff)+B(II,JJ)
      end do
      end do
      end
C==============================================================
      SUBROUTINE INVSQ(as,A,B,D,N)
      IMPLICIT REAL*8(A-H,O-Z)
C       from the texas program, slightly modified by ss 1999
      DIMENSION A(n,n),B(N,N),D(*),as(*)
       parameter(zero=0.0d0,one=1.0d0)
C     **** REPLACES A BY A TO -1/2
C     **** A SHOULD BE A POSITIVE DEFINITE SYMMETRIC MATRIX
C       as same a a but triangular
C     CALL SDIAG2(N,N,A,D,B)
       call eig1(as,b,d,n,n)
      IF(D(1).LT.1.0D-07) WRITE(6,50) D(1)
   50 FORMAT(//,'  ********* WARNING   EIGENVALUE OF OVERLAP:',1F15.10)
      DO I=1,N
      D(I)=ONE/SQRT(D(I))
       end do
      DO I=1,N
      DO J=1,I
      S=ZERO
      DO K=1,N
      S=S+B(I,K)*B(J,K)*D(K)
       end do
         A(I,J)=S
         A(J,I)=S
        end do
        end do
         end
C=====tsetup==================================================
       subroutine tsetup(ijdim,nval,nstron,idm,itst,
     1                      iprint,uniq,idmij,idmw,nweak,
     2                      listw,weak,coro,iextr)
C
C       generates pointer to amplitudes
C       itst(lpar) points to start address for the amplitude for
C       pair lpar.
C       If the pair is a symmetry dependent one the number is
C       meaningless (same as the next)
C       this subroutine also returns ijdim which is the total
C       space needed to store the ampliudes.
C       This is used in mp2iter to reserve memory for the amplitudes
C       
C       called from mp2iter (in the beginning)
C
C       Svein Saebo, Fayetteville, AR summer 1999
C
       implicit integer(a-z)
       logical*1 nstron(nval,*)
       integer*2 idm(*),idmij(*),idmw(*),listw(nval,*)
       logical coro
       dimension itst(*),uniq(*)
       logical weak
C
       ijdim=0
       ijdfm=0
       lpar=0
       itst(1)=0
       lpt=0
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
       iuniq=uniq(lpar)
       if(ii.eq.jj) then
       idd=idm(ii)
       else
       lpt=lpt+1
       idd=idmij(lpt)
          endif
       ijdfm=ijdfm+idd*idd
       if(iuniq.le.0)then
       itst(lpar)=ijdim
       cycle
       endif
       itst(lpar)=ijdim
       ijdim=ijdim+idd*idd
       end do
       end do
C
C       weak
       if(weak) then
       if(.not.coro) then
       mstrong=lpar
          ijdim=0
       lpar=0
       endif
       iextr=0
       do ii=1,nval
       do jj=1,ii
       nwea=listw(ii,jj)
       if(nwea.eq.0) cycle
       lpar=lpar+1
       idd=idmw(ii)
       jdd=idmw(jj)
       idsq=idd*jdd
       ijdfm=ijdfm+idsq
C
       if(coro) then
       iuniq=uniq(lpar)
       else
       iuniq=uniq(lpar+mstrong)
       endif
       if(iuniq.le.0) then
       itst(lpar)=ijdim
       cycle
       endif
C
       itst(lpar)=ijdim
       ijdim=ijdim+idsq
       iextr=iextr+idsq
       enddo
       enddo
       endif
       ijdim=ijdim+1
       if(iprint.ge.3) then
       write(6,*) 'itst(lpar),uniq(lpar),lpar....'
       write(6,666)(itst(iw),uniq(iw),iw,iw=1,lpar)
  666 format(1x,1I10,2I5,1I10,2I5,1I10,2I5,1I10,2I5)
       endif
       if(weak.and.iprint.ge.2) then
       write(6,*)'Additional storage for weak pairs: ',iextr
       write(6,*)'Dimension of all amplitudes: ',ijdim,itst(lpar),lpar
       write(6,*)'Formal dimension of all amplitudes: ',ijdfm
       endif
       end
C======getbas======================================================
       subroutine getbas(ii,jj,idd,jdd,ifu,jfu,idm,ncf)
C       currently not used!!
       implicit integer(a-z)
        integer*2 idm(*)
C       ncf=igetival('ncf')
       idd=idm(ii)
       jdd=idm(jj)
       ifu=ncf*(ii-1)+1
       jfu=ncf*(jj-1)+1
       end
C=======pener======================================================
      subroutine pener(lpar,ii,jj,eij,idd,
     1                itst,irec,coro,nunit,idmxs)
C       calculates the pair energy for pair:ii,jj
C       result in eij
C
C       Svein Saebo, Fayetteville Ar , summer 1999
C


      use memory

        implicit real*8(a-h,o-z)
       logical coro
       dimension itst(*)
c      common/big/bl(30000)
C
         itmxa=mataddr('Tmax')
C
       call matdef('Tij','q',idd,idd)
       itij=mataddr('Tij')
       ijdim=idd*idd
       if(coro) then
       call gett(lpar,idd,itst)
       else
C
      call reada(bl(itij),bl(itmxa),ijdim,idmxs,irec,nunit)
       endif
C
        irij=mataddr('Rij')
       call atoap(bl(irij),idd)
       call ptrac2(bl(irij),bl(itij),ijdim,eij)
       eij=eij+eij
       if(ii.ne.jj) eij=eij+eij
C
C       trap for error condition:
C
        if(eij.lt.-0.5d0.or.eij.gt.1.0d-7) Then
        write(6,*) 'Pair energy absurd:', eij
        write(6,*)ii,jj,idd,lpar,itst(lpar)
        imxxt=idamax(ijdim,bl(itij),1)
        imxxu=idamax(ijdim,bl(irij),1)
        write(6,*) 'T:',imxxt,bl(itij-1+imxxt)
        write(6,*) 'R:',imxxu,bl(irij-1+imxxu)
        call nerror(1,'pener','energy',ii,jj)
        endif
C
       call matrem('Tij')
       end
C========ptrac2==============================================
       subroutine ptrac2(a,b,mn,tr)
C       matrix trace of a the product A+*B
       implicit real*8 (a-h,o-z)
       dimension a(*),b(*)
       parameter(zero=0.0d0)
       sum=zero
       do ij=1,mn
       sum=sum+a(ij)*b(ij)
       end do
       tr=sum
       end
C=======moveT=================================================
       subroutine moveT(lpar,idd,itst)
C       adds deltaT to current T in memory
C
C       Svein Saebo Fayetteville AR summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
       dimension itst(*)
c      common/big/bl(30000)
C
       itfus=mataddr('tful')+itst(lpar)-1
       idtst=mataddr('Rij')-1
       idij=idd*idd
        do imv=1,idij
       bl(itfus+imv)=bl(itfus+imv)+bl(idtst+imv)
       end do
       end
C=======GetT===================================================
       subroutine getT(lpar,idd,itst)
C       fetches the amplitudes for pair lpar from memory
C       result in 'Tij'
C
C       Svein Saebo, Faytteville Ar summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
       dimension itst(*)
c      common/big/bl(30000)
C
       itfus=mataddr('tful')+itst(lpar)-1
       itijs=mataddr('Tij')-1
       idij=idd*idd
       do imv=1,idij
       bl(itijs+imv)=bl(itfus+imv)
       end do
       end
C=========tcoupT====================================================
      subroutine TcouplT(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ifunpair,isorb,uniq,equi,
     3                   nsym,ilist,jlist,lbtab,idmx,
     4                   idmij,lbasij,lpt,ndiag,inndp,
     5                   idd,iff,ttmat,lbasun,lbasw,
     6                   weak,idmw,listw,npairs,coro,
     7                   ifil,krec)
C       calculates coupling terms for the MP2 residuum.
C       with symmetry
C       with this subroutine the multiplications with S are outside
C       the sum over k.
C       Input:
C       ii,jj orbital indices
C       nval number of valence orbitals
C       nstron logical array telling if a pair is strong or not
C       idm  dimensions of the local domains
C       invl inverse pair list
C       listp pair list
C       FF  Fock matrix in localized orbitals
C       lbas local domains
C       ncf number of contracted basis functions
C       itst pointer to the start of T for a given pair
C       ifunpair symmtry relations between AOs
C       isorb  symmetry relations between MOs
C       uniq gives symmetry uniq pair no for unique pairs
C       (see equip in lmp2_sym.f)
C
C       Svein Saebo, Fayetteville AR June 1999
C       Modified by SS August 200
C              

      use memory

       implicit real*8(a-h,o-z)
       logical*1 nstron(nval,*)
       logical coro
       dimension invl(nval,*),listp(*),ilist(*),jlist(*)
       dimension ifunpair(7,*),isorb(nval,*)
       integer uniq(*),equi(2,*)
       dimension ifil(*),krec(*)
       dimension FF(nval,nval)
       integer*2 lbas(*),idm(*),lbtab(idmx,nsym,*),lbasw(*)
       integer*2 idmij(*),lbasij(*),idmw(*),listw(nval,*)
       integer*2 inndp(*),lbasun(*)
       integer itst(*)
       logical ndiag,weak,weak2
       logical ltrnsp
c      common/big/bl(30000)
C
       idmxs=idmx*idmx
       itmxa=mataddr('Tmax')
C
       ixadr=mataddr('x1')
       ioadr=mataddr('ovl')
      call mkunios(jj,ncf,nval,lbas,idm,
     1             lbasij,idmij,invl,listp,nstron,
     2             inndp,lbasun,idun,equi,lbtab,
     3             idmx,nsym,ilist,jlist,weak,
     4             idmw,listw,lbasw,npairs)
       call zerox1(bl(ixadr),ncf,idun,lbasun)
C       loop over valence orbitals kk
       do kk=1,nval
        if(nstron(kk,jj)) goto 100
       weak2=.false.
        kjpa=invl(kk,jj)
       kjpar=listp(kjpa)
       iuniq=uniq(kjpar)
       if(kk.eq.jj)then
       lpx=0
       kdkj=idm(kk)
       kfkj=ncf*(kk-1)+1
       else
       lpx=inndp(kjpar)
       kdkj=idmij(lpx)
       kfkj=ncf*(lpx-1)+1
       endif
       kjdim=kdkj*kdkj
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','q',kdkj,kdkj)
       call matdef('Tij','q',kdkj,kdkj)
       if(iuniq.gt.0)then
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
C
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       
C
       else
C       getts does the same as gett (puts T into 'Tij') when kjpar
C       is not one of the symmetry related parents stored
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
       endif
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','q',kdkj,kdkj)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
C
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
C
       else
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
       endif
       call matredef('Tij','Tkj','q',kdkj,kdkj)
      endif
C       all this just to get Tkj which now is in 'Tkj'
C
       goto 200
C
  100 continue
C       weak pairs
       if(weak) then
       numw=listw(kk,jj)
       if(numw.eq.0) cycle
       weak2=.true.
       kjpar=npairs+numw
       iuniq=uniq(kjpar)
       kdkj=idmw(kk)
       jdkj=idmw(jj)
       kjdim=kdkj*jdkj
       kfu=ncf*(kk-1)+1
       jfu=ncf*(jj-1)+1
       if(coro) then
       kjpp=kjpar
       else
       kjpp=numw
       endif
      if(jj.gt.kk) then              ! case 1
       ltrnsp=.true.
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','r',kdkj,jdkj)
       if(iuniq.gt.0) then
       call matdef('Tij','r',jdkj,kdkj)
       call gettw(kjpp,kjdim,itst)
       else
      call getTs3(kjpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else                     !jjgtkk
C       Case 2              simply get Tkj
       ltrnsp=.false.
       if(iuniq.gt.0)then
       call matdef('Tij','r',kdkj,jdkj)
       call gettw(kjpp,kjdim,itst)
       else
      call getTs3(kjpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matredef('Tij','Tkj','r',kdkj,jdkj)
C       all this just to get Tkj which now is in 'Tkj'
       endif
          call matscal('Tkj',2.0d0)
       else
       cycle
       endif
  200 continue
C       end weak
        fik=-FF(ii,kk)
       call matscal('Tkj',fik)
       itadr=mataddr('Tkj')
C
       if(weak2)then
       if(iuniq.gt.0) then
      call scatad(bl(ixadr),bl(itadr),ncf,kdkj,jdkj,lbasw(kfu),
     1 lbasw(jfu))
       else
       ipt=npairs+iix
       jpt=npairs+jjx
       if(ltrnsp) then
        call scatad(bl(ixadr),bl(itadr),ncf,kdkj,jdkj,
     1 lbtab(1,isym,jpt),lbtab(1,isym,ipt))
       else
        call scatad(bl(ixadr),bl(itadr),ncf,kdkj,jdkj,
     1 lbtab(1,isym,ipt),lbtab(1,isym,jpt))
       endif
       endif
       else
C
        if(iuniq.gt.0) then
       if(kk.eq.jj) then
       call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbas(kfkj),
     1 lbas(kfkj))
        else
       call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbasij(kfkj),
     1 lbasij(kfkj))
       endif
       else
C
       if(lpx.gt.0) ipt=nval+lpx
C
        call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,
     1 lbtab(1,isym,ipt),lbtab(1,isym,ipt))
        endif
       endif
C
       call matrem('Tkj')
       end do
C
C       now multiply by S from both sides in combined local dimension
C
       call matdef('xun','q',idun,idun)
       ixun=mataddr('xun')
       call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
       call matdef('oxm','r',idd,idun)
       call matdef('ovlfj','r',idun,idd)
       call matdef('ovlif','r',idd,idun)
       iovifa=mataddr('ovlif')
       if(ndiag) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasij(iff),lbasun)
       else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbas(iff),lbasun)
       endif
       call secund(ttmat1)
       call matmmult('ovlif','xun','oxm')
       call matpose2('ovlif','ovlfj','n')
       call matrem('ovlif')
       call matmmult('oxm','ovlfj','FCS')
       call secund(ttmat2)
       ttmat=ttmat+ttmat2-ttmat1
       call matrem('ovlfj')
       call matrem('oxm')
       call matrem('xun')
        call matadd('FCS','Rij')
C
       if(ii.eq.jj) then
       call matpose('FCS')
       call matadd('FCS','Rij')
       return
       endif
C       if(ii.eq.jj.and.(.not.weak2)) return
      call mkunios(ii,ncf,nval,lbas,idm,
     1             lbasij,idmij,invl,listp,nstron,
     2             inndp,lbasun,idun,equi,lbtab,
     3             idmx,nsym,ilist,jlist,weak,
     4             idmw,listw,lbasw,npairs)
       call zerox1(bl(ixadr),ncf,idun,lbasun)
C left half done  do right part
       do kk=1,nval
       if(nstron(ii,kk)) goto 300
       weak2=.false.
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
       iuniq=uniq(ikpar)
       if(kk.eq.ii) then
       lpx=0
       idik=idm(ii)
       ifik=ncf*(ii-1)+1
       else
       lpx=inndp(ikpar)
       idik=idmij(lpx)
       ifik=ncf*(lpx-1)+1
       endif
       ikdim=idik*idik
      if(kk.gt.ii) then              !Case 1
       call matdef('Tik','q',idik,idik)
       call matdef('Tij','q',idik,idik)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(ikpar,idik,itst)
       else
C
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif

       else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else                            !Case 2
       call matdef('Tij','q',idik,idik)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
C
       else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
       endif
       call matredef('Tij','Tik','q',idik,idik)
      endif
       goto 400
  300 continue
C       weak pairs
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       weak2=.true.
       ikpar=npairs+numw
       iuniq=uniq(ikpar)
       idik=idmw(ii)
       kdik=idmw(kk)
       ifu=ncf*(ii-1)+1
       kfu=ncf*(kk-1)+1
       ikdim=idik*kdik
       if(coro) then
       ikpp=ikpar
       else
       ikpp=numw
       endif
      if(kk.gt.ii) then              ! case 1
       ltrnsp=.true.
       call matdef('Tik','r',idik,kdik)
       if(iuniq.gt.0) then
       call matdef('Tij','r',kdik,idik)
       call gettw(ikpp,ikdim,itst)
       else
      call getTs3(ikpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else
       ltrnsp=.false.
       if(iuniq.gt.0) then
       call matdef('Tij','r',idik,kdik)
       call gettw(ikpp,ikdim,itst)
       else
      call getTs3(ikpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
        call matredef('Tij','Tik','r',idik,kdik)
       endif
          call matscal('Tik',2.0d0)
       else
       cycle
       endif
  400 continue
C       end weak
        fkj=-FF(kk,jj)       
       call matscal('Tik',fkj)
       itadr=mataddr('Tik')
C
       if(weak2)then
       if(iuniq.gt.0) then
      call scatad(bl(ixadr),bl(itadr),ncf,idik,kdik,lbasw(ifu),
     1 lbasw(kfu))
       else
       ipt=npairs+iix
       jpt=npairs+jjx
       if(ltrnsp) then
        call scatad(bl(ixadr),bl(itadr),ncf,idik,kdik,
     1 lbtab(1,isym,jpt),lbtab(1,isym,ipt))
       else
        call scatad(bl(ixadr),bl(itadr),ncf,idik,kdik,
     1 lbtab(1,isym,ipt),lbtab(1,isym,jpt))
       endif
       endif
       else
C
        if(iuniq.gt.0)then
C
       if(ii.eq.kk) then
       call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbas(ifik),
     1 lbas(ifik))
       else
       call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbasij(ifik),
     1 lbasij(ifik))
       endif
C
        else
C
       if(lpx.gt.0) ipt=nval+lpx
C
        call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,
     1 lbtab(1,isym,ipt),lbtab(1,isym,ipt))
        endif
       endif
       call matrem('Tik')
       end do                            ! end loop over kk
C
C       now multiply by S from both sides in combined local dimension
C
       call matdef('xun','q',idun,idun)
       ixun=mataddr('xun')
       call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
       call matdef('oxm','r',idd,idun)
       call matdef('ovlfj','r',idun,idd)
       call matdef('ovlif','r',idd,idun)
       iovifa=mataddr('ovlif')
       if(ndiag) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasij(iff),lbasun)
       else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbas(iff),lbasun)
       endif
       call secund(ttmat1)
       call matmmult('ovlif','xun','oxm')
       call matpose2('ovlif','ovlfj','n')
       call matrem('ovlif')
       call matmmult('oxm','ovlfj','FCS')
       call secund(ttmat2)
       ttmat=ttmat+ttmat2-ttmat1
       call matrem('ovlfj')
       call matrem('oxm')
       call matrem('xun')
        call matadd('FCS','Rij')
       end
C=======removf========================================
       subroutine removf(irem,orbit,idd,rem,rmv)
       implicit real*8(a-h,o-z)
       dimension orbit(idd,irem)
       dimension rem(*),rmv(*)
       do iy=1,idd
       sum= 0.0d0
       do ix=1,irem
       sum=sum+orbit(iy,ix)**2
       end do
       rem(iy)=sum
       end do
C       search for the irem largest elements in rem
        do jj=1,idd
       rmv(jj)=jj
       end do
 100    itrp=0
        do j=1,idd-1
          if(rem(j).lt.rem(j+1)) then
            x=rem(j+1)
            rem(j+1)=rem(j)
              rem(j)=x
              ii=rmv(j+1)
              rmv(j+1)=rmv(j)
              rmv(j)=ii
            itrp=1
          end if
        end do !j
        if(itrp.eq.1) go to 100
  200 itrp=0
       do ird=1,irem-1
       if(rmv(ird).lt.rmv(ird+1)) then
       ix=rmv(ird+1)
       rmv(ird+1)=rmv(ird)
       rmv(ird)=ix
       itrp=1
       end if
       end do
       if(itrp.eq.1) goto 200
       write(6,*) 'rmv', (rmv(ii),ii=1,irem)
      end
C==========procore===============================
       subroutine procore(C,ncore,ncf,ltab,iprint)
       implicit real*8(a-h,o-z)
       dimension C(ncf,ncore),ltab(*)
       parameter(zero=0.0d0)
C
       do icore=1,ncore
       imax=idamax(ncf,C(1,icore),1)
       if(iprint.ge.3)write(6,*) 'core:',icore,imax
       do jj=1,ncore
       C(imax,jj)=zero
       end do
       ltab(icore)=imax
       end do
C************the procedure above fails for some systems************
C       sort ltab decending order
  100 iret=0
         do isor=1,ncore-1
              if(ltab(isor).lt.ltab(isor+1)) then
              ix=ltab(isor+1)
              ltab(isor+1)=ltab(isor)
              ltab(isor)=ix
              iret=1
              endif
       end do
       if(iret.eq.1) goto 100
C       temporary fix for testing purposes
C************the procedure above fails for some systems************
C       do iz=1,ncore
C       ltab(ncore-iz+1)=iz
C       end do
C       this can not be assumed to be correct in general! remove
C******************************************************************
       if(iprint.ge.2) then
       write(6,*)'CORE AOS:', (ltab(ii),ii=1,ncore)
       endif
       end
C========atoap================================================
      subroutine atoap(A,ncf)
C     **** input a matrix A in full dimension
C     **** output A-0.5A+  returned in A

      use memory

      IMPLICIT REAL*8(A-H,O-Z)
c     common/big/bl(30000)
       dimension A(*)
      parameter(half=0.5d0)
C
      nqd=ncf**2
      call matdef('txp67x','q',ncf,ncf)
       itmst=mataddr('txp67x')-1
      do icp=1,nqd
      bl(itmst+icp)=a(icp)
      end do
       call matpose('txp67x')
       call matscal('txp67x',-half)
      do iadd=1,nqd
       a(iadd)=a(iadd)+bl(itmst+iadd)
      end do
       call matrem('txp67x')
      END
C========getts2==================================================
      subroutine getTs2(lpar,ncf,nval,itst,idm,
     1                  lbas,equi,ifunpair,isorb,ilist,
     2                  jlist,iix,isym,idmij,lbasij,
     3                  lpx,uniq,coro,idmxs,ifil,
     4                  krec)
C
C       only used when nsym gt 0
C       this subroutine is called when pair is not unique: iuniq<0
C       it returns pair lpar in local dimension in Tij
C         if not a unique pair, get its parent and
C         the symmetry-operation yielding the pair form its parent
C       make the symmetry image of the parent to get the correct T
C
C       Svein Saebo, Fayetteville, AR summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
        integer equi(2,*),uniq(*)
       logical ktrnsp
       dimension ifunpair(7,*),isorb(nval,*)
       dimension ifil(*),krec(*)
       dimension ilist(*),jlist(*)
       integer*2 lbas(ncf,*),idm(*)
       integer*2 lbasij(*),idmij(*)
       integer itst(*)
       logical coro
c      common/big/bl(30000)
C
C       get the parent pair
C
       itmxa=mataddr('Tmax')
C
       ktrnsp=.false.
        iprnt=-equi(1,lpar)
       isym=equi(2,lpar)
       ii1=ilist(iprnt)
       ii2=isorb(ii1,isym)
       jj1=jlist(iprnt)
       jj2=isorb(jj1,isym)
       if(lpx.eq.0) then
       iix=ii1
       idd=idm(ii1)
       else
       idd=idmij(lpx)
       ii2a=iabs(ii2)
       ii3=ilist(lpar)
          if(ii2a.ne.ii3) ktrnsp=.true.
       endif
C
C       ktrnsp=.true.  when we take iprnt and perform operation isym
C       the result is the transpose of lpar not lpar itself
C
C
       call matrem('Tij')
C...........................................................
C       get pair: iprnt and perform operation isym on this
       call matdef('Tij','q',idd,idd)
        itt1s=mataddr('Tij')
       if(coro) then
       call gett(iprnt,idd,itst)
       else
       nunit=ifil(uniq(iprnt))
       irec=krec(uniq(iprnt))
      call reada(bl(itt1s),bl(itmxa),idd*idd,idmxs,irec,nunit)
       endif
C
       jsig=isign(1,(ii2*jj2))
       call matdef('tt1','q',idd,idd)
       itt1=mataddr('tt1')
C
C       get symmetry-image of iprnt under operation isym
       if(lpx.eq.0)then
      call symtl(bl(itt1s),bl(itt1),ncf,isym,ii1,
     1           idd,lbas,ifunpair,jsig,ktrnsp)
       else
      call symtl(bl(itt1s),bl(itt1),ncf,isym,lpx,
     1           idd,lbasij,ifunpair,jsig,ktrnsp)
       endif
       call matrem('tt1')
       end
C========getts3==================================================
      subroutine getTs3(lpar,ncf,nval,itst,equi,
     1                  ilist,jlist,isym,idmw,lbasw,
     2                  ifunpair,isorb,iix,jjx,mstrong,
     3                  coro)
C
C       only used when nsym gt 0
C       this subroutine is called when pair is not unique: iuniq<0
C       it returns pair lpar in local dimension in Tij
C         if not a unique pair, get its parent and
C         the symmetry-operation yielding the pair form its parent
C       make the symmetry image of the parent to get the correct T
C
C       Svein Saebo, Fayetteville, AR summer 1999
C

      use memory

       implicit real*8(a-h,o-z)
        integer equi(2,*)
       logical ktrnsp,coro
       dimension ifunpair(7,*),isorb(nval,*)
       dimension ilist(*),jlist(*)
       integer*2 lbasw(ncf,*),idmw(*)
       integer itst(*)
c      common/big/bl(30000)
C
C       get the parent pair
C
       ktrnsp=.false.
        iprnt=-equi(1,lpar)
       isym=equi(2,lpar)
       ii1=ilist(iprnt)
       ii2=isorb(ii1,isym)
       jj1=jlist(iprnt)
       jj2=isorb(jj1,isym)
       iix=ii1
       jjx=jj1
       idd=idmw(ii1)
       jdd=idmw(jj1)
       ii2a=iabs(ii2)
       ii3=ilist(lpar)
       jj3=jlist(lpar)
       id3=idmw(ii3)
       jd3=idmw(jj3)
          if(ii2a.ne.ii3)then
          ktrnsp=.true.
        iix=jj1
        jjx=ii1
       endif
       call matdef('tt1','r',id3,jd3)
       itt1=mataddr('tt1')
C
C       ktrnsp=.true.  when we take iprnt and perform operation isym
C       the result is the transpose of lpar not lpar itself
C
C       get pair: iprnt and perform operation isym on this
       call matdef('Tij','r',idd,jdd)
        itt1s=mataddr('Tij')
       if(.not.coro) iprnt=iprnt-mstrong
       call gettw(iprnt,idd*jdd,itst)
       jsig=isign(1,(ii2*jj2))
       if(ktrnsp) then
       call matdef('tt2','r',jd3,id3)
       itt3=mataddr('tt2')
       else
       itt3=itt1
       endif

C
C       get symmetry-image of iprnt under operation isym
      call symtlw(bl(itt1s),bl(itt1),bl(itt3),ncf,isym,
     1           ii1,jj1,idd,jdd,id3,
     2           jd3,lbasw,ifunpair,jsig,ktrnsp)
       if(ktrnsp) then
       call matpose2('tt2','tt1','n')
       call matrem('tt2')
       endif
       call matrem('Tij')
       call matredef('tt1','Tij','r',id3,jd3)
       end
C========invndp======================================
       subroutine invndp(nstron,nval,indp)
C       generates an iverse list over non diagonal pairs
C       indp(ii,jj) = lpx
C       ii,jj orbital indices, lpx  non-diagonal pair # lpx
C       used to find the correct local basis for this pair
C
C       ss June 2000
C
       implicit integer(a-z)
       logical*1 nstron(nval,*)
       integer*2 indp(*),null
       null=0
       lpn=0
       lpar=0
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
       indp(lpar)=null
       if(ii.eq.jj) cycle
       lpn=lpn+1
       indp(lpar)=lpn
       enddo
       enddo
       end
C============mkltab================================
      subroutine mkltab(ifunpair,lbas,idm,nval,ncf,
     1                  lbtab,idmx,nsym,idmij,lbasij,
     2                  lbasw,idmw,weak,nstron,listw,
     3                  npairs)
C
C       Svein Saebo Fayetteville AR Summer 2000
C
      implicit real*8(a-h,o-z)
       logical*1 nstron(nval,*)
       dimension ifunpair(7,*)
      integer*2 lbas(ncf,*),idm(*),lbtab(idmx,nsym,*)
       integer*2 lbasij(ncf,*),idmij(*)
       integer*2 lbasw(ncf,*),idmw(*),listw(nval,*)
       logical weak
C
      do ii=1,nval
       idd=idm(ii)
       do icp=1,idd
        iao=lbas(icp,ii)
       do irop=1,nsym
        jao=iabs(ifunpair(irop,iao))
        lbtab(icp,irop,ii)=jao
       enddo
       enddo
       enddo
C       nondiagonal pairs
       lpn=0
       do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       if(ii.eq.jj) cycle
       lpn=lpn+1
       idd=idmij(lpn)
       do icp=1,idd
       iao=lbasij(icp,lpn)
       do irop=1,nsym
       jao=iabs(ifunpair(irop,iao))
       lbtab(icp,irop,nval+lpn)=jao
       enddo
       enddo
       enddo
       enddo
C       weak pairs
       if(weak) then
       do ii=1,nval
       idd=idmw(ii)
       do icp=1,idd
       iao=lbasw(icp,ii)
       do irop=1,nsym
       jao=iabs(ifunpair(irop,iao))
       lbtab(icp,irop,npairs+ii)=jao
       enddo
       enddo
       enddo
       endif
       end
C============symche2================================
      subroutine symche2(equi,ifunpair,lbas,idm,nval,
     1                   ncf,iprint,nstron,inndp)
C       if symmetry is used, make the local domain the
C       symmetry image of the local domain for its symmetry
C         unique parent.
C
C       Svein Saebo Starkville MS. fall 1999
C
C       called from mp2init after cheklb
C
      implicit real*8(a-h,o-z)
      integer equi(2,*)
       dimension ifunpair(7,*)
      integer*2 lbas(ncf,*),idm(*),inndp(*)
       logical*1 nstron(nval,*)
      parameter(zero=0.0d0)
C
       lpn=0
       lpar=0
      do ii=1,nval
       do jj=1,ii
       if(nstron(ii,jj)) cycle
       lpar=lpar+1
       if(ii.eq.jj) cycle
       lpn=lpn+1
        if(equi(1,lpar).lt.0) then
       iprnt=-equi(1,lpar)
       irop=equi(2,lpar)
C    ii is the symmetry image of iprnt under operation irop
       ndpar=inndp(iprnt)
       idd=idm(ndpar)
          idm(lpn)=idd
       do icp=1,idd
        iao=lbas(icp,ndpar)
        jao=iabs(ifunpair(irop,iao))
        lbas(icp,lpn)=iao
       end do
C       call lbsort(lbas(1,ndpar),idd)
       endif
C
      end do
       enddo
       end
C=========tcoupl2===================================================
      subroutine Tcoupl2(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ifunpair,isorb,uniq,equi,
     3                   nsym,ilist,jlist,lbtab,idmx,
     4                   idmij,lbasij,lpt,ndiag,inndp,
     5                   idd,iff,ttmat,lbasw,weak,
     6                   idmw,listw,npairs,coro,ifil,
     7                   krec)
C       calculates coupling terms for the MP2 residuum.
C       with this subroutine the multiplications with S are inside
C       the sum over k.
C       Input:
C       ii,jj orbital indices
C       nval number of valence orbitals
C       nstron logical array telling if a pair is strong or not
C       idm  dimensions of the local domains
C       invl inverse pair list
C       listp pair list
C       FF  Fock matrix in localized orbitals
C       lbas local domains
C       ncf number of contracted basis functions
C       itst pointer to the start of T for a given pair
C
C       Svein Saebo, Fayetteville AR August 2000
C              

      use memory

       implicit real*8(a-h,o-z)
       logical*1 nstron(nval,*)
       dimension invl(nval,*),listp(*),ilist(*),jlist(*)
       dimension ifunpair(7,*),isorb(nval,*)
       integer uniq(*),equi(2,*)
       dimension FF(nval,nval)
       integer*2 lbas(*),idm(*),lbtab(idmx,nsym,*),lbasw(ncf,*)
       integer*2 idmij(*),lbasij(*),listw(nval,*)
       integer*2 inndp(*),idmw(*)
       dimension ifil(*),krec(*)
       integer itst(*)
       logical ndiag,weak,weak2,ltrnsp,coro
c      common/big/bl(30000)
C
       idmxs=idmx*idmx
       itmxa=mataddr('Tmax')
C
       ioadr=mataddr('ovl')
C       loop over valence orbitals kk
       do kk=1,nval
       if(nstron(kk,jj)) goto 100
       weak2=.false.
        kjpa=invl(kk,jj)
       kjpar=listp(kjpa)
       iuniq=uniq(kjpar)
       if(kk.eq.jj)then
       kdkj=idm(kk)
       kfkj=ncf*(kk-1)+1
       lpx=0
       else
       lpx=inndp(kjpar)
       kdkj=idmij(lpx)
       kfkj=ncf*(lpx-1)+1
       endif
       kjdim=kdkj*kdkj
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','q',kdkj,kdkj)
       call matdef('Tij','q',kdkj,kdkj)
       if(iuniq.gt.0)then
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       else
C       getts does the same as gett (puts T into 'Tij') when kjpar
C       is not one of the symmetry related parents stored
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
       endif
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','q',kdkj,kdkj)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       else
      call getTs2(kjpar,ncf,nval,itst,idm,
     1           lbas,equi,ifunpair,isorb,ilist,
     2           jlist,ipt,isym,idmij,lbasij,
     3           lpx,uniq,coro,idmxs,ifil,
     4           krec)
       endif
       call matredef('Tij','Tkj','q',kdkj,kdkj)
      endif
C       all this just to get Tkj which now is in 'Tkj'
       goto 200
  100 continue
C       weak pairs
C
       if(weak) then
       numw=listw(kk,jj)
       if(numw.eq.0) cycle
       weak2=.true.
       kjpar=npairs+numw
       iuniq=uniq(kjpar)
       kdkj=idmw(kk)
       jdkj=idmw(jj)
       kjdim=kdkj*jdkj
       if(coro) then
       kjpp=kjpar
       else
       kjpp=numw
       endif
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','r',kdkj,jdkj)
       if(iuniq.gt.0) then
       call matdef('Tij','r',jdkj,kdkj)
       call gettw(kjpp,kjdim,itst)
       else
       ltrnsp=.true.
      call getTs3(kjpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else                     !jjgtkk
C       Case 2              simply get Tkj
       if(iuniq.gt.0)then
       call matdef('Tij','r',kdkj,jdkj)
       call gettw(kjpp,kjdim,itst)
       else
       ltrnsp=.false.
      call getTs3(kjpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matredef('Tij','Tkj','r',kdkj,jdkj)
       endif
C       all this just to get Tkj which now is in 'Tkj'
        call matscal('Tkj',2.0d0)
       else
       cycle
       endif
  200 continue
C       end weak
        fik=-FF(ii,kk)       
       call matscal('Tkj',fik)
C
C       multiply by S in local dimension inside the loop
C
       if(weak2) then
       call matdef('oxm','r',idd,jdkj)
       else
       call matdef('oxm','r',idd,kdkj)
       endif
       call matdef('ovlif','r',idd,kdkj)
       iovifa=mataddr('ovlif')
C       Weak Pairs............
       if(weak2) then
       if(iuniq.gt.0) then
C
       if(ndiag) then
        call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbasw(1,kk))
       else
        call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbasw(1,kk))
       endif !ndiag
C
       else  !iuniq
       ipt=npairs+iix
       if(ltrnsp) ipt=npairs+jjx
C
       if(ndiag) then
        call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbtab(1,isym,ipt))
       else
        call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbtab(1,isym,ipt))
       endif   !ndiag
C
       endif   !iuniq
       else !     weak2
       if(ndiag) then
        if(iuniq.gt.0) then
       if(kk.eq.jj) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbas(kfkj))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbasij(kfkj))
       endif     !  kkeqjj
       else
C
       if(lpx.gt.0) ipt=nval+lpx
C
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbtab(1,isym,ipt))
        endif      !  iuniq
       else        ! ndiag
        if(iuniq.gt.0) then
       if(kk.eq.jj) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbas(kfkj))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbasij(kfkj))
       endif
       else
C
       if(lpx.gt.0) ipt=nval+lpx
C
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbtab(1,isym,ipt))
        endif
       endif
       endif
       call secund(ttmat1)
       call matmmult('ovlif','Tkj','oxm')
       if(weak2) then
       call matdef('ovlit','r',jdkj,idd)
       iovjfa=mataddr('ovlit')
       if(iuniq.gt.0) then
C
       if(ndiag) then
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbasw(1,jj),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbasw(1,jj),lbas(iff))
       endif
C
       else
       ipt=npairs+jjx
       if(ltrnsp) ipt=npairs+iix
       if(ndiag) then
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbtab(1,isym,ipt),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbtab(1,isym,ipt),lbas(iff))
       endif
       endif
       else
       call matdef('ovlit','r',kdkj,idd)
       call matpose2('ovlif','ovlit','n')
       endif
       call matmmult('oxm','ovlit','FCS')
       call secund(ttmat2)
       call matrem('ovlit')
       call matrem('ovlif')
       call matrem('oxm')
       call matadd('FCS','Rij')
       if(ii.eq.jj) then
       call matpose('FCS')
       call matadd('FCS','Rij')
       endif
       ttmat=ttmat+ttmat2-ttmat1
C
       call matrem('Tkj')
       enddo   !loop over kk
       if(ii.eq.jj) return
C left half done  do right part
C       loop over valence orbitals kk
       do kk=1,nval
       if(nstron(ii,kk)) goto 300
       weak2=.false.
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
       iuniq=uniq(ikpar)
       if(kk.eq.ii) then
       lpx=0
       idik=idm(ii)
       ifik=ncf*(ii-1)+1
       else
       lpx=inndp(ikpar)
       idik=idmij(lpx)
       ifik=ncf*(lpx-1)+1
       endif
       ikdim=idik*idik
      if(kk.gt.ii) then              !Case 1
       call matdef('Tik','q',idik,idik)
       call matdef('Tij','q',idik,idik)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else                            !Case 2
       call matdef('Tij','q',idik,idik)
       if(iuniq.gt.0) then
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(iuniq)
       irec=krec(iuniq)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       else
      call getTs2(ikpar,ncf,nval,itst,idm,
     1            lbas,equi,ifunpair,isorb,ilist,
     2            jlist,ipt,isym,idmij,lbasij,
     3            lpx,uniq,coro,idmxs,ifil,
     4            krec)
       endif
       call matredef('Tij','Tik','q',idik,idik)
      endif
       goto 400
  300 continue
C       weak pairs
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       weak2=.true.
       ikpar=npairs+numw
       iuniq=uniq(ikpar)
       idik=idmw(ii)
       kdik=idmw(kk)
       ikdim=idik*kdik
       if(coro) then
       ikpp=ikpar
       else
       ikpp=numw
       endif
      if(kk.gt.ii) then              ! case 1
       call matdef('Tik','r',idik,kdik)
       if(iuniq.gt.0) then
       call matdef('Tij','r',kdik,idik)
       call gettw(ikpp,ikdim,itst)
       else
       ltrnsp=.true.
      call getTs3(ikpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else                     !jjgtkk
C       Case 2              
       if(iuniq.gt.0)then
       call matdef('Tij','r',idik,kdik)
       call gettw(ikpp,ikdim,itst)
       else
       ltrnsp=.false.
      call getTs3(ikpar,ncf,nval,itst,equi,
     1            ilist,jlist,isym,idmw,lbasw,
     2            ifunpair,isorb,iix,jjx,npairs,
     3            coro)
       endif
       call matredef('Tij','Tik','r',idik,kdik)
       endif
        call matscal('Tik',2.0d0)
       else
       cycle
       endif
  400 continue
C       end weak
        fkj=-FF(kk,jj)       
       call matscal('Tik',fkj)
       itadr=mataddr('Tik')
       if(weak2) then
       call matdef('oxm','r',idd,kdik)
       else
       call matdef('oxm','r',idd,idik)
       endif
       call matdef('ovlif','r',idd,idik)
       iovifa=mataddr('ovlif')
       if(weak2) then
       if(iuniq.gt.0) then
       if(ndiag) then
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbasw(1,ii))
       else
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbasw(1,ii))
       endif
       else
       ipt=npairs+iix
       if(ltrnsp) ipt=npairs+jjx
       if(ndiag) then
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbtab(1,isym,ipt))
       else
        call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbtab(1,isym,ipt))
       endif
       endif
       else
       if(ndiag) then
        if(iuniq.gt.0) then
       if(ii.eq.kk) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbas(ifik))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbasij(ifik))
       endif
       else
C
       if(lpx.gt.0) ipt=nval+lpx
C
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbtab(1,isym,ipt))
        endif
       else
        if(iuniq.gt.0) then
       if(ii.eq.kk) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbas(ifik))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbasij(ifik))
       endif
       else
C
       if(lpx.gt.0) ipt=nval+lpx
C
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbtab(1,isym,ipt))
        endif
       endif
       endif
       call secund(ttmat1)
       call matmmult('ovlif','Tik','oxm')
       if(weak2) then
       call matdef('ovlit','r',kdik,idd)
       iovfka=mataddr('ovlit')
       if(iuniq.gt.0) then
       if(ndiag) then
       call compab(bl(ioadr),bl(iovfka),ncf,kdik,idd,
     1                lbasw(1,kk),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovfka),ncf,kdik,idd,
     1                lbasw(1,kk),lbas(iff))
       endif
       else
       ipt=npairs+jjx
       if(ltrnsp) ipt=npairs+iix
       if(ndiag) then
       call compab(bl(ioadr),bl(iovfka),ncf,kdik,idd,
     1                lbtab(1,isym,ipt),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovfka),ncf,kdik,idd,
     1                lbtab(1,isym,ipt),lbas(iff))
       endif
       endif
       else
       call matdef('ovlit','r',idik,idd)
       call matpose2('ovlif','ovlit','n')
       endif
       call matmmult('oxm','ovlit','FCS')
       call secund(ttmat2)
       ttmat=ttmat+ttmat2-ttmat1
       call matrem('ovlit')
       call matrem('ovlif')
       call matrem('oxm')
       call matadd('FCS','Rij')
C
       call matrem('Tik')
       end do                            ! end loop over kk
       end
C================getorb===========================================
      subroutine getorb(idd,norx,lfil)
C
C       this subroutine determines the temporary MO basis for
C       the pair ii,jj
C       lbas is the array containg the local bases
C       the orbitals will be in  'orbi' and 'orbj' and the
C       eigenvalues in 'eigi' and 'eigj'
C
C       Svein Saebo Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
c      common/big/bl(30000)
C
       call matdef('orbi','q',idd,idd)
       call matdef('eigi','v',idd,idd)
       iorba=mataddr('orbi')
       ieiga=mataddr('eigi')
C
       idsq=idd**2
       noru=norx+lfil
       call reaarr(bl(iorba),idsq,noru)
       call reaarr(bl(ieiga),idd,noru)
       end
C================getorbn==========================================
      subroutine getorbn(idd,jdd,noru)
C
C       this subroutine determines the temporary MO basis for
C       the pair ii,jj
C       lbas is the array containg the local bases
C       the orbitals will be in  'orbi' and 'orbj' and the
C       eigenvalues in 'eigi' and 'eigj'
C
C       Svein Saebo Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
c      common/big/bl(30000)
C
       call matdef('orbi','r',idd,jdd)
       call matdef('eigi','v',jdd,jdd)
       iorba=mataddr('orbi')
       ieiga=mataddr('eigi')
C
       idsq=idd*jdd
       call reaarr(bl(iorba),idsq,noru)
       call reaarr(bl(ieiga),jdd,noru)
       end
C=========tcoupl1===================================================
      subroutine Tcoupl1(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ilist,jlist,idmij,lbasij,
     3                   lpt,ndiag,inndp,idd,iff,
     4                   ttmat,lbasw,weak,idijw,listw,
     5                   npairs,coro,idmxs,ifil,krec)
C       calculates coupling terms for the MP2 residuum.
C       with this subroutine the multiplications with S are inside
C       the sum over k.
C       simplified when no symmetry
C       Input:
C       ii,jj orbital indices
C       nval number of valence orbitals
C       nstron logical array telling if a pair is strong or not
C       idm  dimensions of the local domains
C       invl inverse pair list
C       listp pair list
C       FF  Fock matrix in localized orbitals
C       lbas local domains
C       ncf number of contracted basis functions
C       itst pointer to the start of T for a given pair
C
C       Svein Saebo, Fayetteville AR August 2000
C              

      use memory

       implicit real*8(a-h,o-z)
       logical*1 nstron(nval,*)
       dimension invl(nval,*),listp(*),ilist(*),jlist(*)
       dimension FF(nval,nval)
       integer*2 lbas(*),idm(*)
       integer*2 idmij(*),lbasij(*),lbasw(ncf,*),idijw(*)
       integer*2 inndp(*),listw(nval,*)
       integer itst(*)
       logical coro
       dimension ifil(*),krec(*)
       logical ndiag,weak,weak2
c      common/big/bl(30000)
C
       itmxa=mataddr('Tmax')
C
       ioadr=mataddr('ovl')
C       loop over valence orbitals kk
       do kk=1,nval
       if(nstron(kk,jj)) goto 100
       weak2=.false.
        kjpa=invl(kk,jj)
       kjpar=listp(kjpa)
       if(kk.eq.jj)then
       kdkj=idm(kk)
       kfkj=ncf*(kk-1)+1
       lpx=0
       else
       lpx=inndp(kjpar)
       kdkj=idmij(lpx)
       kfkj=ncf*(lpx-1)+1
       endif
       kjdim=kdkj*kdkj
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','q',kdkj,kdkj)
       call matdef('Tij','q',kdkj,kdkj)
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(kjpar)
       irec=krec(kjpar)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','q',kdkj,kdkj)
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(kjpar)
       irec=krec(kjpar)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       call matredef('Tij','Tkj','q',kdkj,kdkj)
      endif
C       all this just to get Tkj which now is in 'Tkj'
       goto 200
  100 continue
C       weak pairs
C
       if(weak) then
       numw=listw(kk,jj)
       if(numw.eq.0) cycle
       weak2=.true.
       kjpar=npairs+numw
       kdkj=idijw(kk)
       jdkj=idijw(jj)
       kjdim=kdkj*jdkj
       if(coro) then
       kjpp=kjpar
       else
       kjpp=numw
       endif
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','r',kdkj,jdkj)
       call matdef('Tij','r',jdkj,kdkj)
       call gettw(kjpp,kjdim,itst)
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','r',kdkj,jdkj)
       call gettw(kjpp,kjdim,itst)
       call matredef('Tij','Tkj','r',kdkj,jdkj)
C       all this just to get Tkj which now is in 'Tkj'
       endif
        call matscal('Tkj',2.0d0)
       else
       cycle
       endif
  200 continue
C       end weak
        fik=-FF(ii,kk)       
       call matscal('Tkj',fik)
C
C       multiply by S in local dimension inside the loop
C
       if(weak2) then
       call matdef('oxm','r',idd,jdkj)
       else
       call matdef('oxm','r',idd,kdkj)
       endif
       call matdef('ovlif','r',idd,kdkj)
       iovifa=mataddr('ovlif')
C
       if(weak2) then
       if(ndiag) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,
     1                lbasij(iff),lbasw(1,kk))
       else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,
     1                lbas(iff),lbasw(1,kk))
       endif
       else    !weak2
C
       if(ndiag) then
       if(kk.eq.jj) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbas(kfkj))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbasij(iff),
     1 lbasij(kfkj))
       endif     !  kkeqjj
       else        ! ndiag
       if(kk.eq.jj) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbas(kfkj))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,kdkj,lbas(iff),
     1 lbasij(kfkj))
       endif       !kkeqjj
       endif
       endif               !weak2
       call secund(ttmat1)
       call matmmult('ovlif','Tkj','oxm')
       if(weak2) then
       call matdef('ovlit','r',jdkj,idd)
       iovjfa=mataddr('ovlit')
       if(ndiag) then
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbasw(1,jj),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovjfa),ncf,jdkj,idd,
     1                lbasw(1,jj),lbas(iff))
       endif
       else
       call matdef('ovlit','r',kdkj,idd)
       call matpose2('ovlif','ovlit','n')
       endif
       call matmmult('oxm','ovlit','FCS')
       call secund(ttmat2)
       call matrem('ovlit')
       call matrem('ovlif')
       call matrem('oxm')
       call matadd('FCS','Rij')
       if(ii.eq.jj) then
       call matpose('FCS')
       call matadd('FCS','Rij')
       endif
       ttmat=ttmat+ttmat2-ttmat1
C
       call matrem('Tkj')
       enddo   !loop over kk
       if(ii.eq.jj) return
C left half done  do right part
C       loop over valence orbitals kk
       do kk=1,nval
       if(nstron(ii,kk)) goto 300
       weak2=.false.
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
       if(kk.eq.ii) then
       lpx=0
       idik=idm(ii)
       ifik=ncf*(ii-1)+1
       else
       lpx=inndp(ikpar)
       idik=idmij(lpx)
       ifik=ncf*(lpx-1)+1
       endif
       ikdim=idik*idik
      if(kk.gt.ii) then              !Case 1
       call matdef('Tik','q',idik,idik)
       call matdef('Tij','q',idik,idik)
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(ikpar)
       irec=krec(ikpar)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else                            !Case 2
       call matdef('Tij','q',idik,idik)
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(ikpar)
       irec=krec(ikpar)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       call matredef('Tij','Tik','q',idik,idik)
      endif
       goto 400
  300 continue
C       weak pairs
C
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       weak2=.true.
       ikpar=npairs+numw
       idik=idijw(ii)
       kdik=idijw(kk)
       ikdim=idik*kdik
       if(coro) then
       ikpp=ikpar
       else
       ikpp=numw
       endif
      if(kk.gt.ii) then              ! case 1
       call matdef('Tik','r',idik,kdik)
       call matdef('Tij','r',kdik,idik)
       call gettw(ikpp,ikdim,itst)
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else
       call matdef('Tij','r',idik,kdik)
       call gettw(ikpp,ikdim,itst)
       call matredef('Tij','Tik','r',idik,kdik)
       endif
        call matscal('Tik',2.0d0)
       else
       cycle
       endif
  400 continue
C       end weak
        fkj=-FF(kk,jj)       
       call matscal('Tik',fkj)
       itadr=mataddr('Tik')
       if(weak2) then
       call matdef('oxm','r',idd,kdik)
       else
       call matdef('oxm','r',idd,idik)
       endif
       call matdef('ovlif','r',idd,idik)
       iovifa=mataddr('ovlif')
C
       if(weak2) then
       if(ndiag) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,
     1                lbasij(iff),lbasw(1,ii))
       else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,
     1                lbas(iff),lbasw(1,ii))
       endif
       else    !weak2
C
       if(ndiag) then
       if(ii.eq.kk) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbas(ifik))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbasij(iff),
     1 lbasij(ifik))
       endif
       else
       if(ii.eq.kk) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbas(ifik))
        else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idik,lbas(iff),
     1 lbasij(ifik))
       endif
       endif
       endif
       call secund(ttmat1)
       call matmmult('ovlif','Tik','oxm')
       if(weak2) then
       call matdef('ovlit','r',kdik,idd)
       iovjfa=mataddr('ovlit')
       if(ndiag) then
       call compab(bl(ioadr),bl(iovjfa),ncf,kdik,idd,
     1                lbasw(1,kk),lbasij(iff))
       else
       call compab(bl(ioadr),bl(iovjfa),ncf,kdik,idd,
     1                lbasw(1,kk),lbas(iff))
       endif
       else
       call matdef('ovlit','r',idik,idd)
       call matpose2('ovlif','ovlit','n')
       endif
       call matmmult('oxm','ovlit','FCS')
       call secund(ttmat2)
       ttmat=ttmat+ttmat2-ttmat1
       call matrem('ovlit')
       call matrem('ovlif')
       call matrem('oxm')
       call matadd('FCS','Rij')
C
       call matrem('Tik')
       end do                            ! end loop over kk
       end
C=========tcoupls===================================================
      subroutine Tcoupls(ii,jj,nval,nstron,idm,
     1                   invl,listp,FF,lbas,ncf,
     2                   itst,ilist,jlist,idmij,lbasij,
     3                   lpt,ndiag,inndp,idd,iff,
     4                   ttmat,lbasun,lbasw,weak,
     5                   idmw,listw,npairs,coro,idmxs,
     6                   ifil,krec)
C       calculates coupling terms for the MP2 residuum.
C       with this subroutine the multiplications with S are outside
C       the sum over k.
C       no symmetry
C       Input:
C       ii,jj orbital indices
C       nval number of valence orbitals
C       nstron logical array telling if a pair is strong or not
C       idm  dimensions of the local domains
C       invl inverse pair list
C       listp pair list
C       FF  Fock matrix in localized orbitals
C       lbas local domains
C       ncf number of contracted basis functions
C       itst pointer to the start of T for a given pair
C
C       Svein Saebo, Fayetteville AR August 2000
C              

      use memory

       implicit real*8(a-h,o-z)
       logical*1 nstron(nval,*)
       dimension invl(nval,*),listp(*),ilist(*),jlist(*)
       dimension FF(nval,nval)
       integer*2 lbas(*),idm(*)
       integer*2 idmij(*),lbasij(*),lbasw(ncf,*)
       integer*2 inndp(*),lbasun(*),idmw(*),listw(nval,*)
       integer itst(*)
       dimension ifil(*),krec(*)
       logical coro
       logical ndiag,weak,weak2
c      common/big/bl(30000)
C
       ioadr=mataddr('ovl')
       ixadr=mataddr('x1')
       itmxa=mataddr('Tmax')
C
C       first determine the combined local local basis
C
      call mkunion(jj,ncf,nval,lbas,idm,
     1            lbasij,idmij,invl,listp,nstron,
     2            inndp,lbasun,idun,weak,listw,
     3            idmw,lbasw)
C
       call zerox1(bl(ixadr),ncf,idun,lbasun)
C       loop over valence orbitals kk
       do kk=1,nval
       if(nstron(kk,jj)) goto 100
       weak2=.false.
        kjpa=invl(kk,jj)
       kjpar=listp(kjpa)
       if(kk.eq.jj)then
       lpx=0
       kdkj=idm(kk)
       kfkj=ncf*(kk-1)+1
       else
       lpx=inndp(kjpar)
       kdkj=idmij(lpx)
       kfkj=ncf*(lpx-1)+1
       endif
       kjdim=kdkj*kdkj
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','q',kdkj,kdkj)
       call matdef('Tij','q',kdkj,kdkj)
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(kjpar)
       irec=krec(kjpar)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','q',kdkj,kdkj)
       if(coro) then
       call gett(kjpar,kdkj,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(kjpar)
       irec=krec(kjpar)
      call reada(bl(itij),bl(itmxa),kjdim,idmxs,irec,nunit)
       endif
       call matredef('Tij','Tkj','q',kdkj,kdkj)
      endif
       goto 200
  100 continue
C       weak pairs
C
       if(weak) then
       numw=listw(kk,jj)
       if(numw.eq.0) cycle
       weak2=.true.
       kjpar=npairs+numw
       kdkj=idmw(kk)
       jdkj=idmw(jj)
       kjdim=kdkj*jdkj
       if(coro) then
       kjpp=kjpar
       else
       kjpp=numw
       endif
      if(jj.gt.kk) then              ! case 1
C       Case 1:  need pair Tkj but jj>kk so we'll have to
C       get Tjk and transpose it
       call matdef('Tkj','r',kdkj,jdkj)
       call matdef('Tij','r',jdkj,kdkj)
       call gettw(kjpp,kjdim,itst)
       call matpose2('Tij','Tkj','n')
       call matrem('Tij')
      else
C       Case 2              simply get Tkj
       call matdef('Tij','r',kdkj,jdkj)
       call gettw(kjpp,kjdim,itst)
       call matredef('Tij','Tkj','r',kdkj,jdkj)
C       all this just to get Tkj which now is in 'Tkj'
       endif
        call matscal('Tkj',2.0d0)
       else
       cycle
       endif
  200 continue
C       end weak
        fik=-FF(ii,kk)       
       call matscal('Tkj',fik)
       itadr=mataddr('Tkj')
C
       if(weak2) then
      call scatad(bl(ixadr),bl(itadr),ncf,kdkj,jdkj,lbasw(1,kk),
     1 lbasw(1,jj))
        else
C
        if(kk.eq.jj) then
      call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbas(kfkj),
     1 lbas(kfkj))
        else
      call scatad(bl(ixadr),bl(itadr),ncf,kdkj,kdkj,lbasij(kfkj),
     1 lbasij(kfkj))
        endif
        endif
       call matrem('Tkj')
        end do
C
C       now multiply by S from both sides in combined local dimension
C
        call matdef('xun','q',idun,idun)
        ixun=mataddr('xun')
      call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
        call matdef('oxm','r',idd,idun)
        call matdef('ovlfj','r',idun,idd)
        call matdef('ovlif','r',idd,idun)
        iovifa=mataddr('ovlif')
        if(ndiag) then
      call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasij(iff),lbasun)
        else
      call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbas(iff),lbasun)
        endif
        call secund(ttmat1)
        call matmmult('ovlif','xun','oxm')
        call matpose2('ovlif','ovlfj','n')
        call matrem('ovlif')
        call matmmult('oxm','ovlfj','FCS')
        call secund(ttmat2)
        ttmat=ttmat+ttmat2-ttmat1
        call matrem('ovlfj')
        call matrem('oxm')
        call matrem('xun')
        call matadd('FCS','Rij')

        if(ii.eq.jj) then
       call matpose('FCS')
       call matadd('FCS','Rij')
       return
       endif
C left half done  do right part
C
      call mkunion(ii,ncf,nval,lbas,idm,
     1            lbasij,idmij,invl,listp,nstron,
     2            inndp,lbasun,idun,weak,listw,
     3            idmw,lbasw)
       call zerox1(bl(ixadr),ncf,idun,lbasun)
       do kk=1,nval
       if(nstron(ii,kk)) goto 300
       weak2=.false.
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
       if(kk.eq.ii) then
       lpx=0
       idik=idm(ii)
       ifik=ncf*(ii-1)+1
       else
       lpx=inndp(ikpar)
       idik=idmij(lpx)
       ifik=ncf*(lpx-1)+1
       endif
       ikdim=idik*idik
      if(kk.gt.ii) then              !Case 1
       call matdef('Tik','q',idik,idik)
       call matdef('Tij','q',idik,idik)
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(ikpar)
       irec=krec(ikpar)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else                            !Case 2
       call matdef('Tij','q',idik,idik)
       if(coro) then
       call gett(ikpar,idik,itst)
       else
       itij=mataddr('Tij')
       nunit=ifil(ikpar)
       irec=krec(ikpar)
      call reada(bl(itij),bl(itmxa),ikdim,idmxs,irec,nunit)
       endif
       call matredef('Tij','Tik','q',idik,idik)
      endif
       goto 400
  300 continue
C       weak pairs
C
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       weak2=.true.
       ikpar=npairs+numw
       idik=idmw(ii)
       kdik=idmw(kk)
       ikdim=idik*kdik
       if(coro) then
       ikpp=ikpar
       else
       ikpp=numw
       endif
      if(kk.gt.ii) then              ! case 1
       call matdef('Tik','r',idik,kdik)
       call matdef('Tij','r',kdik,idik)
       call gettw(ikpp,ikdim,itst)
       call matpose2('Tij','Tik','n')
       call matrem('Tij')
      else
       call matdef('Tij','r',idik,kdik)
       call gettw(ikpp,ikdim,itst)
       call matredef('Tij','Tik','r',idik,kdik)
       endif
        call matscal('Tik',2.0d0)
       else
       cycle
       endif
  400 continue
C       end weak
        fkj=-FF(kk,jj)       
       call matscal('Tik',fkj)
       itadr=mataddr('Tik')
C
       if(weak2) then
      call scatad(bl(ixadr),bl(itadr),ncf,idik,kdik,lbasw(1,ii),
     1 lbasw(1,kk))
       else
       if(ii.eq.kk) then
       call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbas(ifik),
     1 lbas(ifik))
       else
       call scatad(bl(ixadr),bl(itadr),ncf,idik,idik,lbasij(ifik),
     1 lbasij(ifik))
       endif
       endif
C
       call matrem('Tik')
       end do                            ! end loop over kk
C
C       now multiply by S from both sides in combined local dimension
       call matdef('xun','q',idun,idun)
       ixun=mataddr('xun')
       call compab(bl(ixadr),bl(ixun),ncf,idun,idun,
     1                lbasun,lbasun)
       call matdef('oxm','r',idd,idun)
       call matdef('ovlfj','r',idun,idd)
       call matdef('ovlif','r',idd,idun)
       iovifa=mataddr('ovlif')
       if(ndiag) then
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbasij(iff),lbasun)
       else
       call compab(bl(ioadr),bl(iovifa),ncf,idd,idun,
     1               lbas(iff),lbasun)
       endif
       call secund(ttmat1)
       call matmmult('ovlif','xun','oxm')
       call matpose2('ovlif','ovlfj','n')
       call matrem('ovlif')
       call matmmult('oxm','ovlfj','FCS')
       call secund(ttmat2)
       ttmat=ttmat+ttmat2-ttmat1
       call matrem('ovlfj')
       call matrem('oxm')
       call matrem('xun')
        call matadd('FCS','Rij')
       end
C=======compro====================================================
      SUBROUTINE COMPRO(A,B,NCF,IDIM,IFUN)
C
C     **** A FULL MATRIX OF ORDER NCF
C     **** B MATRIX REDUCED IN THE ROWS,DIMENSION IDIM X NCF
C     **** IFUN(IDIM) THE REDUCED BASIS SET
C     **** THE SUBROUTINE REDUCES A IN THE ROWS,RESULT IN B
C
C       utility routine from original lci program
C       slightly modified for PQS , SS 1999
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NCF,NCF),B(IDIM,NCF)
       integer*2 ifun(*)
      DO JJ=1,NCF
      DO II=1,IDIM
C     IFUNC=IFUN(II)
      B(II,JJ)=A(IFUN(ii),JJ)
       end do
       end do
      END
C==============compp===============================================
      SUBROUTINE COMPP(A,B,NCF,MORB,IDIM,IFUN)
C
C     **** SAME AS COMPRO (SEE ABOVE) BUT A IS RECTANGULAR A(NCF,MORB)
C     **** B REDUCES IN ROWS B (IDIM,MORB)
C     **** IFUN(IDIM) THE REDUCED BASIS SET
C     **** THE SUBROUTINE REDUCES A IN THE ROWS,RESULT IN B
C
C       utility routine from original lci program
C       slightly modified for PQS , SS 1999
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NCF,MORB),B(IDIM,MORB)
       integer*2 ifun(*)
      DO  JJ=1,MORB
      DO  II=1,IDIM
      IFUNC=IFUN(II)
      B(II,JJ)=A(IFUNC,JJ)
       end do
       end do
      END
C==============compco=============================================
      SUBROUTINE COMPCO(A,B,NCF,JDIM,JFUN)
C
C     **** A FULL MATRIX OF ORDER NCF
C     **** B MATRIX REDUCED IN THE COLUMNS, DIMENION JDIM X NCF
C     **** THE REDUCED BASIS IS IN JFUN(JDIM)
C     **** THE SUBROUTINE REDUCES A IN THE COLUMNS, RESULT IN B
C
C       utility routine from original lci program
C       slightly modified for PQS , SS 1999
C
C
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NCF,NCF),B(NCF,JDIM)
       integer*2 jfun(*)
      DO  JJ=1,JDIM
      JFUNC=JFUN(JJ)
      DO  II=1,NCF
      B(II,JJ)=A(II,JFUNC)
       enddo
       enddo
      END
C========mkunion========================================
      subroutine mkunion(ii,ncf,nval,lbas,idm,
     1                   lbasij,idmij,invl,listp,nstron,
     2                   inndp,lbmap,idun,weak,listw,
     3                   idmw,lbasw)
C
C       constructs  a combined  local basis
C       for a given MO ii and all MOS kk
C       same as mkunios when no symmetry
C
C         output lbmap(ncf), and idun
C       Svein Saebo Fayetteville, AR Fall 2000
C
       implicit integer (a-z)
       integer*2 lbas(ncf,*),idm(*),lbmap(*),null
       integer*2 lbasij(ncf,*),idmij(*),inndp(*),en
       dimension invl(nval,*),listp(*)
       logical*1 nstron(nval,*)
       integer*2 lbasw(ncf,*),idmw(*),listw(nval,*)
       logical weak
       parameter(null=0,en=1)
C
       do iao=1,ncf
       lbmap(iao) = null
       enddo
C
       do kk=1,nval
       if(nstron(ii,kk)) goto 100
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
C
       if(kk.eq.ii) then
       idik=idm(ii)
       do ild=1,idik
       iao=lbas(ild,ii)
       lbmap(iao)=  en
       enddo
       else
       lpx=inndp(ikpar)
       idik=idmij(lpx)
C
       do ild=1,idik
       iao=lbasij(ild,lpx)
       lbmap(iao)=  en
       enddo
C
       endif
  100 continue
C
C       add weak pairs
C
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       idik=idmw(ii)
       kdik=idmw(kk)
       do ild=1,idik
       iao=lbasw(ild,ii)
       lbmap(iao)= en
       enddo
       do ild=1,kdik
       iao=lbasw(ild,kk)
       lbmap(iao)=en
       enddo
       endif
       enddo       !loop over kk
C
       idun=0
       do iao=1,ncf
       if(lbmap(iao).eq.en) then
       idun=idun+1
       lbmap(idun)=iao
       endif
       enddo
       end
C========mkunios========================================
      subroutine mkunios(ii,ncf,nval,lbas,idm,
     1                   lbasij,idmij,invl,listp,nstron,
     2                   inndp,lbmap,idun,equi,lbtab,
     3                   idmx,nsym,ilist,jlist,weak,
     4                   idmw,listw,lbasw,mstrong)
C
C       constructs  a combined  local basis
C       for a given MO ii and all MOS kk
C
C         output lbmap(ncf), and idun
C
C       Svein Saebo Fayetteville, AR fall 2000
C
       implicit integer (a-z)
       integer*2 lbas(ncf,*),idm(*),lbmap(*),null,idmw(*)
       integer*2 lbasij(ncf,*),idmij(*),inndp(*),en,listw(nval,*)
       integer*2 lbasw(ncf,*)
       dimension invl(nval,*),listp(*),equi(2,*),ilist(*),jlist(*)
       logical*1 nstron(nval,*)
       logical weak
       integer*2 lbtab(idmx,nsym,*)
       parameter(null=0,en=1)
C
       do iao=1,ncf
       lbmap(iao) = null
       enddo
C
       do kk=1,nval
       if(nstron(ii,kk)) goto 100
       ikpa=invl(ii,kk)
       ikpar=listp(ikpa)
       iprnt=-equi(1,ikpar)
       if(iprnt.gt.0) then
       isym=equi(2,ikpar)
       ipt=ilist(iprnt)
       endif
C
       if(kk.eq.ii) then
       idik=idm(ii)
C
       if(iprnt.gt.0) then
       do ild=1,idik
       iao=lbtab(ild,isym,ipt)
       lbmap(iao)= en
       enddo
C
       else
C
       do ild=1,idik
       iao=lbas(ild,ii)
       lbmap(iao)=  en
       enddo
       endif
C
       else
C
       lpx=inndp(ikpar)
       idik=idmij(lpx)
       if(iprnt.gt.0) then
       ipt=lpx+nval
       do ild=1,idik
       iao=lbtab(ild,isym,ipt)
       lbmap(iao)= en
       enddo
C
       else
C
       do ild=1,idik
       iao=lbasij(ild,lpx)
       lbmap(iao)=  en
       enddo
C
       endif
       endif
       cycle
  100 continue
C
C       add weak pairs
C
       if(weak) then
       numw=listw(ii,kk)
       if(numw.eq.0) cycle
       ikpar=mstrong+numw
       iprnt=-equi(1,ikpar)
       if(iprnt.gt.0) then
       isym=equi(2,ikpar)
       ipt=ilist(iprnt)
       jpt=jlist(iprnt)
       idik=idmw(ipt)
       kdik=idmw(jpt)
       ipt=ipt+mstrong
       jpt=jpt+mstrong
       else
       idik=idmw(ii)
       kdik=idmw(kk)
       endif
C
       if(iprnt.gt.0) then
       do ild=1,idik
       iao=lbtab(ild,isym,ipt)
       lbmap(iao)= en
       enddo
       else
       do ild=1,idik
       iao=lbasw(ild,ii)
       lbmap(iao)=  en
       enddo
       endif
       if(iprnt.gt.0) then
       do ild=1,kdik
       iao=lbtab(ild,isym,jpt)
       lbmap(iao)= en
       enddo
       else
       do ild=1,kdik
       iao=lbasw(ild,kk)
       lbmap(iao)=en
       enddo
       endif
       endif
       enddo       !loop over kk
C
       idun=0
       do iao=1,ncf
       if(lbmap(iao).eq.en) then
       idun=idun+1
       lbmap(idun)=iao
       endif
       enddo
       end
C========zerox1=====================================================
       subroutine zerox1(x1,ncf,idun,lbasc)
C       zeroes out the part of the full matrix used to
C       store couplings terms
C       called from tcouplt and tcoupls
C       
C       Svein Saebo, Fall 2000
C
       implicit real*8(a-h,o-z)
       dimension x1(ncf,*)
       integer*2 lbasc(*)
       parameter(zero=0.0d0)
       do jj=1,idun
       icol=lbasc(jj)
C       call zeroit(x1(1,icol),ncf)
       do kk=1,idun
       irow=lbasc(kk)
       x1(irow,icol)=zero
       enddo
       enddo
       end
C========penalty========================================
       subroutine penalty(lpar,ncf,idd,ifu,lbas,
     1                       itst,lpt,lbasij,ndiag)
C       calculates penalty function when redundant local basis
C       is used
C

      use memory

       implicit real*8(a-h,o-z)
       integer*2 lbas(*),lbasij(*)
       logical ndiag
       dimension itst(*)
c      common/big/bl(30000)
       parameter(fact=1.0d0)
C
       isdsa=mataddr('sds')
C
C       matrix sds is defined in full dimension early in lmp2trans
C       matrix is costructed in mp2init
C
       call matprint('Rij',6)
       call matdef('TSDS','q',idd,idd)
       call matdef('lsds','q',idd,idd)
       ilsds=mataddr('lsds')
       if(ndiag) then
         call compab(bl(isdsa),bl(ilsds),ncf,idd,idd,
     1  lbasij(ifu),lbasij(ifu))       
       else
         call compab(bl(isdsa),bl(ilsds),ncf,idd,idd,
     1  lbas(ifu),lbas(ifu))       
       endif
       call matprint('lsds',6)
       call matprint('Tij',6)
       call matmmult('lsds','Tij','TSDS')
          call matscal('TSDS',fact)
       call matprint('TSDS',6)
       call matadd('TSDS','Rij')       
       call matrem('lsds')
       call matrem('TSDS')
       end
       subroutine mkwkli(nval,listp,listw)
       implicit integer(a-z)
       dimension listp(*)
       integer*2 listw(nval,*)
       ij=0
       nwea=0
       do ii=1,nval
       do jj=1,ii
       ij=ij+1
       if(listp(ij).eq.0) then
       nwea =nwea +1
       listw(ii,jj)=nwea
       listw(jj,ii)=nwea
       else
       listw(ii,jj)=0
       listw(jj,ii)=0
       endif
       enddo
       enddo
       end
C=======estimem===================================================
      subroutine estimem(memor,meme,nval,idm,iorbst,
     2                   ipps,ncf,ncore)
C       generates pointers for MO basis used for the weak pairs
C       will be stored in ipps(nval),iorbs(nval)
C
C       Svein Saebo Fayetteville AR fall 2000
C
       implicit real*8(a-h,o-z)
       integer*2 idm(*)
       dimension iorbst(*),ipps(*)
       memor=0
       meme=0
       do ii=1,nval
       ipps(ii)=meme
       iorbst(ii)=memor
       idd=idm(ii)
       memor=memor+idd**2
       meme=meme+idd
       end do
       iorbst(nval+1)=memor
       avdim=float(meme)/float(nval)
       prcnt=100.d0*avdim/(ncf-nval-ncore)
       write(6,*) 'Average domains-size for weak pairs:',avdim,
     1 prcnt,'% of full'
       end
C================neworb3===========================================
       subroutine neworb3(ncf,nval,idm,lbas,ipps,
     1                      iorbst,orbs,eigs)
C       this subroutine determines the temporary MO basis for
C       the correlated orbitals
C
C       Svein Saebo Fayetteville AR June 1999
C

      use memory

       implicit real*8(a-h,o-z)
       integer*2 idm(*),lbas(*)
       dimension iorbst(*),ipps(*)
       dimension orbs(*),eigs(*)
c      common/big/bl(30000)
C
       iovad=mataddr('ovl')
       ifoca=mataddr('fck')
C
C       call matprint('fck',6)
       do ii=1,nval
       idd=idm(ii)
       ifu=ncf*(ii-1)+1
       ipoint=iorbst(ii)
       ipp=ipps(ii)
       call matdef('orbi','q',idd,idd)
       call matdef('eigi','d',idd,idd)
       iorba=mataddr('orbi')
       ieiga=mataddr('eigi')
       call matdef('Sii','q',idd,idd)
       isiia=mataddr('Sii')
       call matdef('Fii','q',idd,idd)
       ifiia=mataddr('Fii')
       call matdef('SFii','q',idd,idd)
C       construct S**-1/2
       call matdef('Siis','s',idd,idd)
       isiis=mataddr('Siis')
      call compab(bl(iovad),bl(isiia),ncf,idd,idd,
     1            lbas(ifu),lbas(ifu))
       call matcopy('Sii','Siis')
      call invsq(bl(isiis),bl(isiia),bl(iorba),bl(ieiga),idd)
       call matrem('Siis')
C       form S**-1/2 * F
      call compab(bl(ifoca),bl(ifiia),ncf,idd,idd,
     1             lbas(ifu),lbas(ifu))
C       call matprint('Fii',6)
       call matmmult('Sii','Fii','SFii')
       call matmmult('SFii','Sii','Fii')
       call matrem('SFii')
C......diagonolize Fii to find i part of orbitals and eigenvalues
C       call matprint('Fii',6)
C********  replace sdiag2 with dspevx **********
       call matdef('fsym','s',idd,idd)
       call matcopy('Fii','fsym')
       ifsym=mataddr('fsym')
       call eig1(bl(ifsym),bl(ifiia),bl(ieiga),idd,idd)
C       call sdiag2(idd,idd,bl(ifiia),bl(ieiga),bl(ifiia))
       call matrem('fsym')
C********  replace sdiag2 with dspevx **********
C       call matprint('Fii',6)
       call matmmult('Sii','Fii','orbi')
C       call matprint('orbi',6)
C       call matprint('eigi',6)
       call matrem('Fii')
       call matrem('Sii')
       do jm=1,idd*idd
       orbs(ipoint+jm)=bl(iorba-1+jm)
       end do
       do jm=1,idd
       eigs(ipp+jm)=bl(ieiga-1+jm)
       end do
       call matrem('eigi')
       call matrem('orbi')
       end do
       end
       subroutine wrida(a,b,ilen,mlen,nrec,nunit)
       implicit real*8(a-h,o-z)
       dimension a(ilen),b(mlen)
       call move(a,b,ilen)
       write(unit=nunit,rec=nrec) b
       end
       subroutine reada(a,b,ilen,mlen,nrec,nunit)
       implicit real*8(a-h,o-z)
       dimension a(ilen),b(mlen)
       read(unit=nunit,rec=nrec) b
       call move(b,a,ilen)
       end
C==============norma===========================
       subroutine norma(a,idd,jdd)
       implicit real*8(a-h,o-z)
       dimension A(idd,jdd)
       do jj=1,jdd
       sum=0.0d0
       do ii=1,idd
       sum=sum+A(ii,jj)**2
       enddo
       sum=1.0d0/sqrt(sum)
       do ii=1,idd
       A(ii,jj)=A(ii,jj)*sum
       enddo
       enddo
       end

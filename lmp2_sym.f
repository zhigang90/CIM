C========orbpr===============================================
      subroutine orbpr(ncf,nval,nsym,C,S,
     1     ifunpair,isorb,iprint,thrsy)
      implicit real*8(a-h,o-z)
      dimension C(ncf,nval),ifunpair(7,*),isorb(nval,*)
      dimension S(ncf,ncf)
          logical halt
          halt=.false.
      do jj=1,nsym
      do ii=1,nval
      isorb(ii,jj)=0
      end do
      end do
      if (nsym.le.0) return
      if(iprint.ge.3) write(6,*)'Symmetry related orbitals:'
C      loop over symmetry operations
      do isym=1,nsym
C      loop over valence orbitals
      do ii=1,nval
      do jj=1,ii
      call symeq(C(1,ii),C(1,jj),S,ncf,ifunpair,
     1 isym,nsym,neq,thrsy)
        if (neq.eq.1) then
      isorb(ii,isym)=jj
      isorb(jj,isym)=ii
      cycle
      endif
      if(neq.eq.-1) then
      isorb(ii,isym)=-jj
      isorb(jj,isym)=-ii
      cycle
      endif
      end do
      end do
      end do
      do isym=1,nsym
      do ii=1,nval
          if(isorb(ii,isym).eq.0) then
          halt=.true.
      write(6,*)'orbital no ',ii,' operation ',isym,' has no relative'
          endif
          if(iprint.ge.3) write(6,100) isym,ii,isorb(ii,isym)
  100  format(1x,'Symmetry Op #.:',I4,'  VO ',I4,'  -> ',I4)
      enddo
      enddo
          if(halt) then
      write(6,*)'There are orbitals without symmetry equivalent '
      Write(6,*)'relatives.  The program will stop!'
      write(6,*)'This may happen if you are running the program'
      write(6,*)'with a non-Abelian symmerty group, if so distort'
      write(6,*)'the geometry to make the point group Abelian.'
      write(6,*)
      Write(6,*)'This may also be  caused by poorly converged '
      write(6,*)'localization. Each localized orbital must transform'
      write(6,*)'into itself or one of the other localized orbitals'
      write(6,*)
      write(6,*)'Try rerunning with EQOR=N where N <8 eg 6 or 5'
      write(6,*) 'Rather than stopping we will turn symmetry off'
      write(6,*) '********WARNING   SYMMETRY TURNED OFF***********'
C     write(6,*)'stopping......SORRY'
C     stop 'orbpr'
      nsym=0
      endif
      end
C=======symeq===========================================
      subroutine symeq(C1,C2,S,ncf,ifunpair,isym,nsym,neq,dev)

      use memory

      implicit real*8(a-h,o-z)
      dimension C1(*),C2(*),ifunpair(7,*),S(ncf,ncf)
c     common/big/bl(30000)
      parameter(zero=0.0d0,one=1.0d0)
C         parameter(dev=1.0d-7,zero=0.0d0,one=1.0d0)
      neq=0
C      maximum deviation from + or - 1 to be accepted as symmetry
C      equivalent
      ivcst=mataddr('VSC')-1
      do my=1,ncf
      sum=zero
      do ny=1,ncf
      sum=sum+S(ny,my)*C2(ny)
      end do
      bl(ivcst+my)=sum
      end do
      sum=zero
      do my=1,ncf
      ny=ifunpair(isym,my)
        if(ny.gt.0) sum=sum+C1(ny)*bl(ivcst+my)
      if(ny.lt.0) sum=sum-C1(iabs(ny))*bl(ivcst+my)      
      end do
      if(abs(sum-one).lt.dev) neq=1
      if(abs(sum+one).lt.dev) neq=-1
      end
C=======equip=========================================
      subroutine equip(equi,uniq,nval,isorb,nsym,nstron,invl,
     1 listp,munik,iprint,listw,weak,npairs,munk,munndia)
C      this subroutine generates a list of unique pairs.
C      if a pair is not the last of a set of symmetry related
C      pairs equi(lpar) is negative - parent.operation
C      see below in output section
C
C      Input:
C      nval number of valence orbitals
C      isorb array with symmetry relations bewteen MOS
C      nsym      number of symmetry operations
C      nstron logical array, false if pair is strong
C      invl  inverse pair list
C      listp      pair list
C      see comments in plist for details
C      Output:
C      equi
C      uniq
C      munik      number of symmetry unique pairs
C
C      Svein Saebo, Fayetteville AR Summer 1999
C      September 2000 added weak pairs
C            
C      Originally written by SS for Disco, modified for PQS
      implicit integer(a-z)
      dimension equi(2,*)
      dimension uniq(*)
      logical weak
      dimension isorb(nval,*),invl(nval,nval),listp(*)
      integer*2 listw(nval,*)
      logical*1 nstron(nval,nval)
C
      lpar=0
      munik=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      equi(1,lpar)=0
      equi(2,lpar)=0
      uniq(lpar)=0
      enddo
      enddo
C                  Weak
      if(weak) then
      do ii=1,nval
      do jj=1,ii
      ij=invl(ii,jj)
      ijw=listp(ij)
      if(ijw.ne.0)cycle
      lpar=lpar+1
      equi(1,lpar)=0
      equi(2,lpar)=0
      uniq(lpar)=0
      enddo
      enddo
      endif
      write(6,*) 's+w',lpar
C                  Weak
      if(nsym.gt.0) then
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      if(ii.ne.jj) cycle
      do isym=1,nsym
      iny=isorb(ii,isym)
      jny=isorb(jj,isym)
      in1=iabs(iny)
      jn1=iabs(jny)
      ijny=invl(in1,jn1)
      kpar=listp(ijny)
      if(kpar.le.0) then
      write(6,*) 'something is wrong: ',kpar,ii,jj,iny,jny,isym
      write(6,*) 'This error condition is somtimes encoutered'
      write(6,*) 'when the SCF calculation was done in a previous'
      write(6,*) 'run and the files on /scr/user are not compatible'
      write(6,*) 'with the current MP2 calculation'
      Write(6,*) 'Redo the SCF calcuation and try again'
      call nerror(1,'equip','not strong',0,0)
      endif
      if(kpar.gt.lpar) goto 1000
      enddo
      mult=1
       do jsym=nsym,1,-1
       iny=isorb(ii,jsym)
       jny=isorb(jj,jsym)
       in1=iabs(iny)
       jn1=iabs(jny)
       ijny=invl(in1,jn1)
       kpar=listp(ijny)
      if(kpar.lt.lpar.and.equi(1,kpar).eq.0) then
      mult=mult+1
C      equi(kpar)=-float(lpar)-0.01d0*jsym
      equi(1,kpar)=-lpar
      equi(2,kpar)=jsym
      endif
       enddo
      equi(1,lpar)=mult
      munik=munik+1
      uniq(lpar)=munik
 1000 continue
        enddo
      enddo
      mundiag=munik
C      nondiagonal pairs
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      if(ii.eq.jj) cycle
      do isym=1,nsym
      iny=isorb(ii,isym)
      jny=isorb(jj,isym)
      in1=iabs(iny)
      jn1=iabs(jny)
      ijny=invl(in1,jn1)
      kpar=listp(ijny)
      if(kpar.le.0) then
      write(6,*) 'something is wrong: ',kpar,ii,jj,iny,jny,isym
      write(6,*) 'This error condition is somtimes encoutered'
      write(6,*) 'when the SCF calculation was done in a previous'
      write(6,*) 'run and the files on /scr/user are not compatible'
      write(6,*) 'with the current MP2 calculation'
      Write(6,*) 'Redo the SCF calcuation and try again'
      call nerror(1,'equip','not strong',0,0)
      endif
      if(kpar.gt.lpar) goto 2000
      enddo
      mult=1
       do jsym=nsym,1,-1
       iny=isorb(ii,jsym)
       jny=isorb(jj,jsym)
       in1=iabs(iny)
       jn1=iabs(jny)
       ijny=invl(in1,jn1)
       kpar=listp(ijny)
      if(kpar.lt.lpar.and.equi(1,kpar).eq.0) then
      mult=mult+1
C      equi(kpar)=-float(lpar)-0.01d0*jsym
      equi(1,kpar)=-lpar
      equi(2,kpar)=jsym
      endif
       enddo
      equi(1,lpar)=mult
      munik=munik+1
      uniq(lpar)=munik
 2000 continue
        enddo
      enddo
      munk=munik
      munndia=munk-mundiag
C                  Weak
      if(weak) then
      do ii=1,nval
      do jj=1,ii
      numw=listw(ii,jj)
      if(numw.eq.0)cycle
      lpar=npairs+numw
      do isym=1,nsym
      iny=isorb(ii,isym)
      jny=isorb(jj,isym)
      in1=iabs(iny)
      jn1=iabs(jny)
      nwe2=listw(in1,jn1)
      if(nwe2.le.0) then
      write(6,*) 'something is wrong: ',nwe2,ii,jj,iny,jny,isym
      write(6,*) 'This error condition is somtimes encoutered'
      write(6,*) 'when the SCF calculation was done in a previous'
      write(6,*) 'run and the files on /scr/user are not compatible'
      write(6,*) 'the current MP2 calculation'
      Write(6,*) 'Redo the SCF calcuation and try again'
      call nerror(1,'equip','not weak',0,0)
      endif
      if(nwe2.gt.numw) goto 3000
      enddo
      mult=1
       do jsym=nsym,1,-1
       iny=isorb(ii,jsym)
       jny=isorb(jj,jsym)
       in1=iabs(iny)
       jn1=iabs(jny)
      nwe2=listw(in1,jn1)
      kpar=npairs+nwe2
      if(kpar.lt.lpar.and.equi(1,kpar).eq.0) then
      mult=mult+1
      equi(1,kpar)=-lpar
      equi(2,kpar)=jsym
      endif
       enddo
      equi(1,lpar)=mult
      munk=munk+1
      uniq(lpar)=munk
 3000 continue
        enddo
      enddo
      endif
C            Weak end
C      Output Section
      if(iprint.ge.3) then
      write(6,*) 'Symmetry related Pairs:'
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      iuniq=uniq(lpar)
      if(iuniq.gt.0) then
      write(6,*) 'Pair no ',lpar,' is unique pair no ',iuniq
      write(6,*) 'Multiplicity = ',equi(1,lpar)
      else
      iprnt=-equi(1,lpar)
      irop=equi(2,lpar)
C      irop=-(equi(lpar)+iprnt-0.001d0)*100
      write(6,*) 'Pair no ',lpar,' is related to ',iprnt
      write(6,*) 'symmetry operation no: ',irop
      endif
      enddo
      enddo
C            Weak
      if(weak) then
      do ii=1,nval
      do jj=1,ii
      ij=listw(ii,jj)
      if(ij.eq.0) cycle
      lpar=lpar+1
      iuniq=uniq(lpar)
      if(iuniq.gt.0) then
      write(6,*) 'Pair no ',lpar,' is unique pair no ',iuniq
      write(6,*) 'Multiplicity = ',equi(1,lpar)
      else
      iprnt=-equi(1,lpar)
      irop=equi(2,lpar)
C      irop=-(equi(lpar)+iprnt-0.001d0)*100
      write(6,*) 'Pair no ',lpar,' is related to ',iprnt
      write(6,*) 'symmetry operation no: ',irop
      endif
      enddo
      enddo
      endif
      endif
C                              weak
       else
C      below when no symmetry:
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj))  cycle
      lpar=lpar+1
      equi(1,lpar)=1
      uniq(lpar)=lpar
      enddo
      enddo
      munik=lpar
C            Weak
      if(weak)then
      do ii=1,nval
      do jj=1,ii
      ij=listw(ii,jj)
      if(ij.eq.0) cycle
      lpar=lpar+1
      equi(1,lpar)=1
      uniq(lpar)=lpar
      enddo
      enddo
      munk=lpar
      endif
C            Weak
      endif
      end
C=======symt===================================================
      subroutine symt(t1,t2,ncf,isym,ifunpair,jsign)
C      this subroutine generates the symmetry image of t1
C      onto t2 under operation isym.
C      jsign is sign(isorb(i,isym)*isorb(j,isym))
C      ifunpair gives symmetry relations bewteen AOs
C      the subroutine can be used for both residia and
C      amplitudes.
C      from the oroginal lci program.
C      modified for PQS
C
C      Svein Saebo, Fayetteville AR, summer 1999
C
      implicit real*8(a-h,o-z)
      dimension t1(ncf,ncf),t2(ncf,ncf)
      dimension ifunpair(7,*)
      do iao=1,ncf
      iny=ifunpair(isym,iao)
      iny1=iabs(iny)
      do jao=1,ncf
      jny=ifunpair(isym,jao)
      jny1=iabs(jny)
      tt=t1(iao,jao)
      ksig=jsign*iny*jny
      jsig=isign(1,ksig)
      if(jsig.eq.-1) tt=-tt
      t2(iny1,jny1)=tt
      enddo
      enddo
      end
C========unmos================================================
      subroutine unmos(unorb,nval,nsym,isorb,iprint)
      implicit integer(a-z)
      dimension unorb(2,*),isorb(nval,*)
C
      if(nsym.gt.0) then
      do ii=1,nval
      unorb(1,ii)=0
      unorb(2,ii)=0
      enddo
      nounik=0
      do ii=1,nval
      do isym=1,nsym
      iny=isorb(ii,isym)
      iny1=iabs(iny)
      if(iny1.gt.ii) goto 100
      enddo
      do jsym=nsym,1,-1
      iny=isorb(ii,jsym)
      iny1=iabs(iny)
      if(iny1.lt.ii.and.unorb(1,iny1).eq.0) then
      unorb(1,iny1)=-ii
      unorb(2,iny1)=jsym
C      unorb(iny1)=-float(ii)-0.01*jsym
      endif
      enddo
      nounik=nounik+1
      unorb(1,ii)=nounik
  100 continue
      enddo
C      output section
      if(iprint.ge.3) then
      write(6,*) 'mo list:'
      do ii=1,nval
      iuniq=unorb(1,ii)
      if(iuniq.gt.0) then
      write(6,*)'orbital no ',ii,' is unique orbital # ',iuniq
      else
      iprnt=-unorb(1,ii)
      irop=unorb(2,ii)
C      irop=-(unorb(ii)+iprnt-0.001)*100
      write(6,*) 'orbital no ',ii,' is related to ',iprnt,
     1' symmetry operation ',irop
           endif
      enddo
      endif
      else
      do ii=1,nval
      unorb(1,ii)=ii
      enddo
      endif
      end
C============symchec================================
      subroutine symchec(unorb,ifunpair,lbas,idm,nval,ncf,iprint)
C      if symmetry is used, make the local domain the
C      symmetry image of the local domain for its symmetry
C         unique parent.
C
C      Svein Saebo Starkville MS. fall 1999
C
C      called from mp2init after cheklb
C
      implicit real*8(a-h,o-z)
      integer unorb(2,*)
      dimension ifunpair(7,*)
      integer*2 lbas(ncf,*),idm(*)
      parameter(zero=0.0d0)
C
      do ii=1,nval
        if(unorb(1,ii).lt.0) then
      iprnt=-unorb(1,ii)
      irop=unorb(2,ii)
C    ii is the symmetry image of iprnt under operation irop
      idd=idm(iprnt)
          idm(ii)=idd
       do icp=1,idd
        iao=lbas(icp,iprnt)
        jao=iabs(ifunpair(irop,iao))
        lbas(icp,ii)=jao
       end do
C      call lbsort(lbas(1,ii),idd)
      if(iprint.ge.3) then
      write(6,*) 'local domain for ',ii,' constructed from ',iprnt
          write(6,*) 'symmetry operation ', irop
      WRITE(6,*) 'Dimension: ', idd
          write(6,*) (lbas(iwr,ii),iwr=1,idd)
      endif
      endif
      end do
       end
C============================================================
      subroutine lbsort(lbas,idd)
      implicit integer(a-z)
      integer*2 lbas(*),lbte
  100 continue
      iret=0
      do iao=1,idd-1
      if(lbas(iao).gt.lbas(iao+1)) then
      lbte=lbas(iao)
      lbas(iao)=lbas(iao+1)
      lbas(iao+1)=lbte
      iret=1
      endif      
      enddo
      if(iret.eq.1) goto 100
      end
C===============chectab===========================================
      subroutine chectab(equi,uniq,unorb,nstron,listp,invli,
     1nval,nsym,iprint)
C     check tables for debugging only
      implicit integer(a-z)
      dimension uniq(*),equi(2,*), unorb(2,*)
      dimension listp(*),invli(nval,*)
      logical*1 nstron(nval,nval)
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
      ijx=invli(ii,jj)
      ijs=listp(ijx)
      if(ijs.ne.lpar) then
      write(6,*) 'ERROR', ii,jj,lpar,ijs
      write(6,*) 'invli'
      do kk=1,nval
      write(6,*) (invli(kk,ll),ll=1,kk)
      enddo
      do kk=1,nval
      write(6,*) (nstron(kk,ll),ll=1,kk)
      enddo
      stop 67
      endif
C
      if(ii.ne.jj) cycle
      ij=uniq(ijs)
      iiun=unorb(1,ii)
      if(ij.gt.0.and.iiun.le.0) then
      write(6,*)'WARNING',ii,jj,ijs,ij,iiun,unorb(2,ii)
      endif
      if(ij.le.0.and.iiun.gt.0) then
      write(6,*)'WARNING',ii,jj,ijs,ij,iiun,unorb(2,ii)
      endif
      enddo
      enddo
      end
C=======symtl===================================================
      subroutine symtl(t1,t2,ncf,isym,ii,idd,
     1                     lbas,ifunpair,jsign,ktrnsp)
C      this subroutine generates the symmetry image of t1
C      under operation isym.
C      result in t1
C      t1 amplitude in local dimension
C      it is assumed that the local basis of a symmetry
C      related pair is the exact symmetry image of the
C      local basis of the parent
C      jsign is sign(isorb(i,isym)*isorb(j,isym))
C      ifunpair gives symmetry relations bewteen AOs
C      the subroutine can be used for both residia and
C      amplitudes.
C
C      Svein Saebo, Fayetteville AR, summer 2000
C
      implicit real*8(a-h,o-z)
      logical ktrnsp
      dimension t1(idd,*),t2(idd,*)
      integer*2 lbas(ncf,*)
      dimension ifunpair(7,*)
C
      do jx=1,idd
      jao=lbas(jx,ii)
      jny=isign(1,ifunpair(isym,jao))
      do ix=1,idd
      iao=lbas(ix,ii)
      iny=isign(1,ifunpair(isym,iao))
      ksig=jsign*iny*jny
      if(ktrnsp) then
      t2(ix,jx)=ksig*t1(ix,jx)
      else
      t1(ix,jx)=ksig*t1(ix,jx)
      endif
      enddo
      enddo
      if(ktrnsp) call matpose2('tt1','Tij','n')
      end
C=======symtlw==================================================
      subroutine symtlw(t1,t2,t3,ncf,isym,
     1                     ii,jj,idd,jdd,id3,
     2                     jd3,lbas,ifunpair,jsign,ktrnsp)
C      this subroutine generates the symmetry image of t1
C      under operation isym.
C      result in t1
C      t1 amplitude in local dimension
C      it is assumed that the local basis of a symmetry
C      related pair is the exact symmetry image of the
C      local basis of the parent
C      jsign is sign(isorb(i,isym)*isorb(j,isym))
C      ifunpair gives symmetry relations bewteen AOs
C      the subroutine can be used for both residia and
C      amplitudes.
C
C      Svein Saebo, Fayetteville AR, summer 2000
C
      implicit real*8(a-h,o-z)
      logical ktrnsp
      dimension t1(idd,jdd),t2(id3,jd3),t3(jd3,id3)
      integer*2 lbas(ncf,*)
      dimension ifunpair(7,*)
C
      do jx=1,jdd
      jao=lbas(jx,jj)
      jny=isign(1,ifunpair(isym,jao))
      do ix=1,idd
      iao=lbas(ix,ii)
      iny=isign(1,ifunpair(isym,iao))
      ksig=jsign*iny*jny
      if(ktrnsp) then
      t3(ix,jx)=ksig*t1(ix,jx)
      else
      t2(ix,jx)=ksig*t1(ix,jx)
      endif
      enddo
      enddo
      end
C============mklbtab================================
      subroutine mklbtab(ifunpair,lbas,idm,nval,ncf,lbtab,idmx
     1 ,nsym)
C
C      Svein Saebo Fayetteville AR Summer 2000
C
      implicit real*8(a-h,o-z)
      dimension ifunpair(7,*)
      integer*2 lbas(ncf,*),idm(*),lbtab(idmx,nsym,nval)
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
C      do ii=1,nval
C      do krop=1,nsym
C      write(6,*) 'orb= ',ii,' oper= ',krop
C      idd=idm(ii)
C      write(6,*) (lbtab(jcp,krop,ii),jcp=1,idd)
C      enddo
C      enddo
       end
C=======equiw=========================================
      subroutine equipw(equiw,nval,isorb,nsym,nstron,
     1                    invl,listp,munik,iprint)
C      this subroutine generates a list of unique weak pairs.
C      if a pair is not the last of a set of symmetry related
C      pairs equi(lpar) is negative - parent.operation
C      see below in output section
C
C      Input:
C      nval number of valence orbitals
C      isorb array with symmetry relations bewteen MOS
C      nsym      number of symmetry operations
C      nstron logical array, false if pair is strong
C      invl  inverse pair list
C      listp      pair list
C      see comments in plist for details
C      Output:
C
C      Svein Saebo, Fayetteville AR Fall 2000
C            
      implicit integer(a-z)
      dimension equiw(*)
      dimension isorb(nval,*),invl(nval,nval),listp(*)
      logical*1 nstron(nval,nval)
C
      lpar=0
      munik=0
      nweak=0
      do ii=1,nval
      do jj=1,ii
      lpar=lpar+1
      if(listp(lpar).ne.0) cycle
      nweak=nweak+1
      equiw(nweak)=0
      enddo
      enddo
      if(nsym.gt.0) then
      lpar=0
      nweak=0
      do ii=1,nval
      do jj=1,ii
      lpar=lpar+1
      if(listp(lpar).ne.0) cycle
      nweak=nweak+1
      do isym=1,nsym
      iny=isorb(ii,isym)
      jny=isorb(jj,isym)
      in1=iabs(iny)
      jn1=iabs(jny)
      ijny=invl(in1,jn1)
      kpar=listp(ijny)
      if(kpar.ne.0) then
      write(6,*) ii,jj,iny,jny,isym,kpar
      call nerror(1,'equiw','not weak',nweak,lpar)
      endif
      if(ijny.gt.lpar) goto 1000
      enddo
      mult=1
       do jsym=nsym,1,-1
       iny=isorb(ii,jsym)
       jny=isorb(jj,jsym)
       in1=iabs(iny)
       jn1=iabs(jny)
       ijny=invl(in1,jn1)
      if(ijny.lt.lpar.and.equiw(kpar).eq.0) then
      mult=mult+1
      equiw(kpar)=-lpar
      equiw(kpar)=jsym
      endif
       enddo
      equiw(lpar)=mult
      munik=munik+1
C      uniq(lpar)=munik
 1000 continue
        enddo
      enddo
C      Output Section
      if(iprint.ge.3) then
      write(6,*) 'Symmetry related Pairs:'
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj)) cycle
      lpar=lpar+1
C      iuniq=uniq(lpar)
      if(iuniq.gt.0) then
      write(6,*) 'Pair no ',lpar,' is unique pair no ',iuniq
      write(6,*) 'Multiplicity = ',equiw(lpar)
      else
      iprnt=-equiw(lpar)
      irop=equiw(lpar)
      write(6,*) 'Pair no ',lpar,' is related to ',iprnt
      write(6,*) 'symmetry operation no: ',irop
      endif
      enddo
      enddo
      endif
       else
      lpar=0
      do ii=1,nval
      do jj=1,ii
      if(nstron(ii,jj))  cycle
      lpar=lpar+1
      equiw(lpar)=1
C      uniq(lpar)=lpar
      enddo
      enddo
      munik=lpar
      endif
      end

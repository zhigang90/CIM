c-----------------------------------------------------------------------
      subroutine bmat(iter,iconv)
c-----------------------------------------------------------------------
c     Optimization in natural internal (BMAT) coordinates.
c     F.Eckert Jan 96 - Jun 96
c
      implicit double precision (a-h,o-z)
      include "common/tapes"
      include "common/big"
      include "common/cmpp"
      include "common/cgeom"
      include "common/czmat"
      include "common/coptg"
      common /cgadr_opt/ ix,ig,ie,ih
      common /tapeneu/ inpx,inp2,ioutx,ipunx,ixx,icond,inpf3,inpf4
      common /cunit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     >pi,acc
      common /opti/ itopti,mxopti,chmax,thract,iangs
c
c      dimension cart(3,zncen),xyz(3,zncen),cgrad(3,zncen)
      dimension ian(200)
c      dimension symbol(zncen)
      character*2 sym, symbol(200),unit*8
      character*32 tempf
      logical hesinput
      hesinput=modelhes .or. idohes.ge.0
c
      if(iout.ne.6.and.firstiter) then
        if(opttyp(1:3).eq.'QSD') then
          if (opttyp(1:4).ne.'QSDP' .and. iroot.eq.1) then
            write (6,'('' Quadratic Steepest Descent -'',
     >               '' Minimum Search'')')
          else if (opttyp(1:4).ne.'QSDP' .and. iroot.gt.1) then
            write (6,'(/'' Quadratic Steepest Descent - '',
     >               '' Transition State Search'')')
          else
            if (idohes.eq.1 .or. hdetect) then
              write (6,'(/'' Quadratic Steepest Descent - '',
     >           '' Reaction Path Following using exact Hessian'')')
            else
              write (6,'(/'' Quadratic Steepest Descent - '',
     >           '' Reaction Path Following using updated Hessian'')')
            end if
          end if
        else if(opttyp(1:2).eq.'RF') then
          write (6,'(/'' Rational Function Geometry Optimization'')')
        end if
      end if
c
      ibase=icorr(0)
      inpx=inp
      inp2=31
      ioutx=iout
      ipunx=ipun
      ixx=inp
      icond=32
      inpf3=33
      inpf4=34
      if(iter.eq.1) then
c     open files for bmat
c     31=inp2, 32=icond, 33=inpf3
        call tmpnm(inp2,tempf)
        call assgn(inp2,tempf,isize,0)
        if(fileb.ne.' ') then
          call assgn(icond,fileb,isize,1)
        else
          call tmpnm(icond,tempf)
          call assgn(icond,tempf,isize,0)
        end if
        call tmpnm(inpf3,tempf)
        call assgn(inpf3,tempf,isize,0)
        call tmpnm(inpf4,tempf)
        call assgn(inpf4,tempf,isize,0)
      end if
      if(iprocs.ne.0) icond=0
c
      zero=0.0d0
      half=0.5d0
      one=1.0d0
      two=2.0d0
      three=3.0d0
      four=4.0d0
      five=5.0d0
      ten=10.0d0
      ten6=1.0d+6
      tenm8=1.0d-8
      acc=1.0d-14
c
      call getvar('PI',pi,unit,ity,nv,1,1)
      call getvar('ANG',ang,unit,ity,nv,1,1)
      ajoule=0.22936757d0
c.....following constants are not used here!
      debye=0.39342658d0
      cbm=0.117946332d30
      evolt=0.036749026d0
      ckalmo=1.5936018d-3
      dkel=3.1667911d-6
      cmm1=4.5563352d-6
      hertz=1.51982983d-16
c
      itopti=iter
      mxopti=100
      chmax=1.0d0
c
c     Optimization History
c
      n3=3*zncen
c...
      mxpt=zngeo(1)
      ix = icorr(n*mxpt)
      ig = icorr(n*mxpt)
      ie = icorr(mxpt)
      ixnew=icorr(3*zncen)
      ixstr=icorr(3*zncen)
      ignew=icorr(n)
      icart=icorr(3*zncen)
      igrad=icorr(3*zncen)
c.....scale initial hessian
      call sclhes(q(ih),1)
c.....Restart Option
      if(irest.ne.0.and.recrs.ne.recsv) then
        call gethis(q(ix),q(ig),q(ie),1)
        call sclgrd(q(ix),dum,q(ig),1)
        call gethes(q(ix),q(ig),q(ie),q(ih))
      end if
c...  Get History
      call gethis(q(ix),q(ig),q(ie),2)
      call sclgrd(q(ix),dum,q(ig),1)
c...  model hessian
      if (modelhes) then
       call model(q(ih),n3)
       write (iout,'(/'' Making model hessian for BMAT coordinates'')')
       if(iout.ne.6.and.firstiter)
     >    write (6,'(/'' Making model hessian for BMAT coordinates'')')
      end if
c
      idone=0
      enext=0d0
      conv=0d0
      convl=0d0
c
      etot=zenergy
c
c     stop if all gradients are very small
c
      call fmove(q(ig+n*(npt-1)),q(ignew),n)
      call absmax(n,q(ignew),iimax,gradmax)
      if (gradmax.lt.1.0d-6) then
        write(iout,12) gradmax
        if(icond.ne.0) write(icond,12) gradmax
 12     format(/1x,'Convergence reached: max. gradient',f10.7)
        call fmove(q(ix+n*(npt-1)),q(ixnew),n)
        goto 99
      end if
c...  check for torques
      if (n.eq.n3) then
       itorque=0
       call torque(n,npt,q(ix),q(ig),itorque,iprhis,iout)
      end if
c
c     Prepare Molpro Input for BMAT: It works internally in Angstrom
c     and writes it's own history file fort.33!
c     New: do NOT sort out Dummy-Atoms from Z-Matrix-Input !!
c
      na=0
      do 10 i=1,zncen
      if (zcharge(i).eq.0d0) modelhes=.false.
       na=na+1
       nek=3*na
c      Coordinates have to be in Angstrom
       q(icart+nek-3)=zcoord(1,i)/ang
       q(icart+nek-2)=zcoord(2,i)/ang
       q(icart+nek-1)=zcoord(3,i)/ang
c      Forces (negative Gradients) have to be in [aJoule/Angstrom]=mdyne
       q(igrad+nek-3)=-zgrad(1,i)*ang/ajoule
       q(igrad+nek-2)=-zgrad(2,i)*ang/ajoule
       q(igrad+nek-1)=-zgrad(3,i)*ang/ajoule
       ian(na)=int(zcharge(i))
       sym=ztag(i)(1:min(len(ztag(i)),len(sym)))
       if (sym(2:2).ge.'0'.and.sym(2:2).le.'9') sym (2:2)=' '
       symbol(na)=sym
c      write(6,'(1x,a,i3,3a,3f12.8,a,3f12.8)')
c    >    'na=',na,'  sym=',symbol(na),
c    >    '  coord=',q(icart+nek-3),q(icart+nek-2),q(icart+nek-1),
c    >    '  grad=', q(igrad+nek-3),q(igrad+nek-2),q(igrad+nek-1)
c       q(isymb+na-1)=sym
 10   continue
      call fmove(q(ix),q(ixstr),nek)
      call dscal_X(nek,one/ang,q(ixstr),1)
c
c     Call Pulays program
c
cpp
      write(*,'("initial Cart.",/,(3f12.6))') (q(icart+k-1),k=1,nek)
c
      call ovly50(q(ixstr),q(icart),ian,symbol,q(igrad),q(ixnew),ickeep,
     >            etot,na,enext,opttyp,q(ih),hesinput)
c
      write(*,'("final Cart.",/,(3f12.6))') (q(ixnew+k-1),k=1,nek)
c          
c
c     scale new Coordinates to [a.u]
c
      nek=3*zncen
      do i=1,nek
       q(ixnew+i-1)=q(ixnew+i-1)*ang
       do j=1,nek
       q(ih+(i-1)*nek+j-1)=q(ih+(i-1)*nek+j-1)*ajoule/(ang*ang)
       end do
      end do
c...  back scale positions, gradient, and hessian
 99   call sclgrd(q(ix),q(ixnew),q(ig),-1)
      call sclhes(q(ih),-1)
c
c.....check convergence and set variables to new values
      convl=0d0
      iconv=0
c.....symmetry clean coordinates. This assumes that coordinates are in z-matrix order
      call clean_zmat(q(ixnew))
      call optcnv(q(ie),q(ix),q(ixnew),q(ig),q(ih),enext,conv,convl,
     >            iconv)
      if(idone.ne.0) then
        conv=0d0
        iconv=1
      end if
c...  determine ihflag for numerical hessian
      nsteps=nsteps+1
      if (nsteps.lt.idohes) then
        ihflag=0
      else
        nsteps=0
        ihflag=1
      end if
      if (idohes.eq.0) ihflag=0
      if (opttyp(1:3).eq.'QSD'.and.itrap.ne.0.and.nohess.eq.0) then
        nsteps=0
        ihflag=1
      end if
      call corlsr(ibase)
      return
      end
c-----------------------------------------------------------------------
      subroutine ovly50(x,xa,ia,symb,f,xy,ickeep,etot,na,enext,opttyp,
     > hes,modelhes)
c-----------------------------------------------------------------------
c     THIS IS A TOTALLY NEW VERSION OF BMAT, MADE MODULAR
c     IT INCLUDES AUTOMATIC GEOMETRY OPTIMIZATION AND AUTOMATIC
c     INTERNAL COORDINATE GENERATION.
c     This was originally written by P.Pulay and G.Fogarasi.
c     Implementation in Molpro and some error fixes by F.Eckert.
c     The code in lowercase is the new one.
c........................................................................
      implicit real*8 (a-h,o-z)
      real*8 l,m
      logical writ,lgdi,lfdi,lmur,qonly,modelhes
      character*16 opttyp,line*128
      include "common/big"
      include "common/cmpp"
      common /tapeneu/ inp,inp2,iout,ipun,ix,icond,inpf3,inpf4
      common /cunit/ ang,debye,cbm,ajoule,evolt,ckalmo,dkel,cmm1,hertz
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,
     >pi,acc
      parameter (ma=200,lq=3*ma,mq=lq-6,md=100,maxcard=3*mq)
c     MA IS THE MAX. NUMBER OF ATOMS, LQ IS THE MAX. NUMBER OF
c     CARTESISANS, MQ IS THE MAX. NUMBER OF INTERNAL COORDINATES,
c     MD IS THE MAX. NUMBER OF GDIIS GEOMETRY STEPS (OR BFGS)
      dimension xa(3,ma),f(lq),fi(mq),ibcode(6,maxcard),ibcontr(20,mq)
      dimension xy(3,ma),x(3,ma)
      dimension qq(mq),qq1(mq),fc(mq),cq(mq)
      dimension b(54,mq)
      dimension qd(lq,md),fd(lq,md),a(md,md),energy(md)
c      common bl(10000)
      common /opti/ itopti,mxopti,chmax,thract,iangs
c     THESE ARE THE PARAMETERS FOR THE GLOBAL GEOMETRY OPTIMIZATION
c      DIMENSION L(LQ), M(LQ),IA(MA),SYMB(MA)
      dimension l(lq), m(lq),ia(ma)
      dimension hes(3*na,3*na)
      character*2 symb(ma)
      dimension ifix(md)
      dimension lopt(30)
c     THIS IS THE NUMBER OF GEOMETRIES OMITTED FROM GDIIS:
      save mexclu
      data mexclu/0/
c
c     SET DEFAULTS
c
      ibase=icorr(0)
      et=0.d0
c      pi=dacos(-1.0d0)
      ix=inp2
c     IX=INPF3
      nfix=0
      ngeo=0
      ickeep=0
      if (iangs.eq.0) then
         unitl=ang
      else
         unitl=one
      end if
c     SET THE UNITS CORRECTLY - ANGSTROMS OR BOHR
c     THIS PART WORKS INTERNALLY IN ANGSTROMS, UNLIKE THE OTHER PARTS
      lgdi=.false.
      lfdi=.false.
      lmur=.false.
      writ=.false.
      qonly=.false.
c       NO WRITING OF THE B MATRIX
c
c     SET DEFAULT SHIFT POINT FOR MACHB
c
      shftpt=0.5d0*pi
c     SET THE DEFAULT MAXIMUM COORDINATE CHANGE
      coormx=0.1d0
c
c       DIMENSION OF THE B MATRIX MUST CONFORM TO B(NCMAX,NQMAX)
c
         deltae=zero
c        DELTAE IS THE ENERGY ASSOCIATED WITH THE OFFSET FORCES
c        THE ZERO IS ARBITRARY BUT IT IS TAKEN AS sum Q(I)*F(I)
c
      do 8 i=1,30
       lopt(i)=0
 8    continue
      call readopx (lopt,na,etot,nfix,ifix,ngeo,lgdi,lfdi,shftpt,
     1  lmur,writ,coormx,inp,inp2,iout,ix)
c     if(lopt(18).eq.0.and.itopti.eq.1) then
c       close(33,status='delete')
c       open(33,file='fort.33')
c     end if
      if (opttyp(1:2).eq.'DI') lgdi=.true.
      if (lopt(19).eq.1) then
        ickeep=1
      end if
      etot1=etot
      if (lopt(7).eq.1) then
        iangs=1
        unitl=one
      end if
      if (na.gt.0.and.na.le.ma) go to 90
      write (iout,1030) na
      if(icond.ne.0) write (icond,1030) na
 1030 format (/1x,'TOO MANY OR TOO FEW NUCLEI' ,i5)
c      CALL ERROR
       call fehler
   90 continue
      ipu=0
      if (lopt(5).eq.1) ipu=ipun
      if (lopt(16).eq.0) qonly=.true.
        nek=3*na
        if (lopt(16).eq.1) call intcoor(x,ia,na)
c
  200   continue
c
c      call bread(na,inp,iout,nq,nprim,ibcode,ibcontr,writ)
      call bread(na,ix,iout,nq,nprim,ibcode,maxcard,ibcontr,mq,writ)
      nq1=na*3-6
      if (nq.ne.nq1.and.na.gt.2) then
         write (iout,210) nq,nq1
         if(icond.ne.0) write (icond,210) nq,nq1
 210  format (/' WARNING: THERE ARE ',i4, ' INT.COOR., 3N-6= ',i4)
      end if
      call machbnew(na,xa,nq,qonly,shftpt,ibcode,ibcontr,b,qq)
c
      memgmat=icorr(nq**2)
      memproj=icorr(nq**2)
      memginv=icorr(nq**2)
c  Form the G matrix to get its generalized inverse and the nonred. projector
      call formgmat(nq,nek,b,ibcontr,q(memgmat))
c  Form generalized inverse
      write(*,*) 'G matrix' 
c      call prntmat(nq,nq,nq,q(memgmat))
c
      tol=1.0d-6
      call dcopy_x(nq**2,q(memgmat),1,q(memginv),1)
      memvec=icorr(nq**2)
      memx=icorr(nq)
      memy=icorr(nq)
      call geninv(nq,q(memginv),det,tol,q(memvec),
     &            q(memx),q(memy),ninv,q(memproj))
      nred=nq-ninv
      write(*,*) 'There are ',nred,' redundancies and ',ninv,
     &           ' nonredundant coordinates',' 3N-6=',nek-6 
c     write(*,*) 'generalized inverse of G in bmat'
c      call prntmat(nq,nq,nq,q(memginv))
c     write(*,*) 'projector of G in bmat'
c      call prntmat(nq,nq,nq,q(memproj))
      call sdiag2(nq,nq,q(memginv),q(memx),q(memvec))
      call outvec(q(memx),nq,'Eigenvalues of G- in bmat')
c      call outsqr(q(memginv),nq,nq,nq,'Eigenvectors of G- in bmat')
c      call outsqr(q(memproj),nq,nq,nq,
c     &           'Projector to the nonredundant subspace of G')
       write(*,*) 'projector to the nonredundant space of G in bmat'
c      call prntmat(nq,nq,nq,q(memproj))
c
c  The following subroutine returns the (possibly redundant) internal forces
c  in phi, the generalized inverse of the G matrix in q(igmat), and the projector GG-
c  in q(iginv), as well as the maximum internal force in fmax1, and its coordinate icmy
c
      call intforc(nq,nek,b,ibcontr,f, fi,qq,iout,icond,q(memgmat),
     1             q(memginv),icmy,fmax1)
cpp
c      write(*,*)  'Internal coordinates and forces'
c      do kk=1,nq
c        write(*,95) kk,qq(kk),fi(kk)
  95  format(i5,2f14.6)
c      end do
      call outvec(fi,nq,'internal forces')
      do 330 i=1,nq
         if (lopt(5).eq.1) write (ipun,1170) qq(i),fi(i),i
 1170    format (2f10.6,i4,6x,f10.6)
 330    continue
c
 335   continue
      if (lopt(2).lt.1) then
c       NO FMAT ,NO GEOMETRY RELAXATION
         call corlsr(ibase)
         return
      end if
c     READ IN THE FORCE CONSTANTS AND INVERT THEM. ALSO DEAL WITH FIXC
c...    READ INTERNAL COORDINATES AND FORCES FOR DIIS OR BFGS
        call readdiag(fc,nq,nfix,ifix,ix,iout,icond)
        if(iprocs.eq.0) rewind inpf3
        ihq=icorr(nq*nq)
        if(modelhes) then
         ibb=icorr(nq*nek)
         call hes_bmat(hes,q(ihq),q(ibb),b,ibcontr,nq,nek)
cpp
c       write(*,*) 'hes_bmat in ovly50'
c         call prntmat(nq,nq,nq,hes)
c
        end if
c       First optimization step
        if (itopti.lt.2) then
c        nur Notloesung ! lopt(18): 'CONT'=Restart-Option
         if (lopt(18).eq.0) then
          if(iprocs.eq.0) then
           write(inpf3,345) itopti,etot
           do 789 j=1,nq
             write(inpf3,1170) qq(j),fi(j),j,qq(j)
 789       continue
           rewind inpf3
          end if
         else
c         Normalerweise liest bmat die History immer aus INP ?!!
          do 790 i=1,ngeo
           if(iprocs.eq.0) read(ix,'(a)') line
           if(nprocs.gt.1) call global_brdcst_str(line)
           read(line,795) irnd,energy(i),i1,fi1,i2,fi2,i3,fi3
 795       format(6x,i3,7x,f15.8,7x,i3,f10.6,7x,i3,f10.6,7x,i3,f10.6)
           if(iprocs.eq.0) write(inpf3,345) irnd,energy(i),i1,fi1,i2,
     >                                      fi2,i3,fi3
c          HEADER CONTAINING ROUND NO., TOTAL ENERGY, MAX. FORCE,
c          MAX. COORD. CHANGE, MAX. ESTIMATED ENERGY CHANGE
          do 790 j=1,nq
            if(iprocs.eq.0) read(ix,'(a)') line
            if(nprocs.gt.1) call global_brdcst_str(line)
            read(line,1170) qd(j,i),fd(j,i),ii,qd(j,i+1)
            if(iprocs.eq.0) write(inpf3,1170) qd(j,i),
     >                           fd(j,i),j,qd(j,i+1)
c           WRITE A COPY ON INPF3 SO THAT IT COULD BE CONTINUED
  790     continue
          if(iprocs.eq.0) rewind inpf3
          end if
         end if
c         AUTOMATIC GEOMETRY OPTIMIZATION - READ FROM INPF3 (FILE 33)
c         AUTOMATIC OPTIMIZATION PAST THE FIRST STEP
          i=1
 680      if(iprocs.eq.0) then
            read(inpf3,'(a)',end=681) line
            goto 682
 681        line='---'
          end if
 682      continue
          if (nprocs.gt.1) call global_brdcst_str(line)
          if(line(1:3).eq.'---') goto 720
          read(line,795) iii,energy(i)
          do 700 j=1,nq
           if(iprocs.eq.0) then
             read(inpf3,'(a)',end=701) line
             goto 702
 701         line(1:3)='---'
           end if
 702       continue
           if(nprocs.gt.1) call global_brdcst_str(line)
           if(line(1:3).eq.'---') goto 720
           read(line,1170) qd(j,i),fd(j,i),ii,qd(j,i+1)
 700      continue
          i=i+1
          go to 680
 720      continue
cstart cray
c;          if(iprocs.eq.0) backspace inpf3
cend
c         old expression:NGEO=I-1
c         new expression: NGEO=I-2
          ngeo=i-2
      if (ngeo.gt.0) then
c       TRY TO CORRECT THE INTERNAL COORDINATES IN THE CASE THAT
c       A CHANGE BY 2*PI OCCURRED (FOR A TORSION)
        if(lopt(1).eq.0) then
          dxmax=0d0
          do 730 j=1,nq
            call fixtor(qd(j,ngeo+2),qq(j))
            dxmax=max(dxmax,abs(qd(j,ngeo+2)-qq(j)))
            qd(j,ngeo+2)=qq(j)
            fd(j,ngeo+2)=fi(j)
c            CALL FIXTOR(QD(J,NGEO+1),QQ(J))
c            QD(J,NGEO+1)=QQ(J)
c            FD(J,NGEO+1)=FI(J)
            energy(ngeo+2)=etot1
 730      continue
        end if
      end if
      etot1=etot+deltae
c     Optimization step
      ngeo=ngeo+1
cpp
c      write(*,*) 'Line 504 in OVLY50, ngeo=',ngeo,' ninv=',ninv
c     First step is always Quasi Newton step
      if (ngeo.lt.2) then
        if (ninv.lt.nq) then
c   If there is redundancy, transfer the diagonal matrix to q(memginv)
          call transfdfc(nq,fc,q(memginv))
          call relaxqred(nq,fi,q(memginv),q(memproj),tol,
     &                    cq,edecr,dismax,icmx,one)
        else
          write(iout,'(/'' Quasi Newton step'')')
cpp
       write(*,*) 'Calling relaxd'
          call relaxd(nq,fi,fc,cq,edecr,dismax,icmx,one)
        end if
cpp
c       write(*,*) (cq(kk),kk=1,nq)
       edec=edecr*ajoule
       et=etot1+edec
c       write (iout,490) edecr,edec,et
       if(icond.ne.0) write(icond,490) edecr,edec,et
       goto 555
      end if
      if (opttyp(1:2).eq.'AH'.or.opttyp(1:2).eq.'AK'.or.
     >    opttyp(1:2).eq.'RF'.or.opttyp(1:2).eq.'SR'.or.
     >    opttyp(1:2).eq.'DI'.or.opttyp(1:3).eq.'QSD') then
        iqd=icorr(nq*ngeo)
        ifd=icorr(nq*ngeo)
        ienerg=icorr(ngeo)
        do jj=1,ngeo
         q(ienerg+jj-1)=energy(jj+1)
         do kk=1,nq
          q(iqd+nq*(jj-1)+kk-1)=qd(kk,jj+1)
          q(ifd+nq*(jj-1)+kk-1)=fd(kk,jj+1)
         end do
        end do
c        call outsqr(q(iQD),nq,nq,ngeo,'coordinates')
c        call outsqr(q(iFD),nq,nq,ngeo,'forces')
        call opt_mol(nq,ngeo,q(iqd),q(ifd),q(ienerg),fc,cq,
     >   et,dismax,icmx,q(ihq))
      else
       write(iout,'(/'' Quasi Newton step'')')
cpp
       write(*,'(/'' Quasi Newton step'')')
       if (modelhes) then
         if(ninv.lt.nq) then
           call relaxqred(nq,fi,q(ihq),cq,q(memproj),
     &                    tol,edecr,dismax,icmx,one)
         else
           call relaxq(nq,fi,q(ihq),cq,edecr,dismax,icmx,one)
         end if
       else
         if(ninv.lt.nq) then
           call transfdfc(nq,fc,q(memginv))
           call relaxqred(nq,fi,q(memginv),q(memproj),tol,
     &                     cq,edecr,dismax,icmx,one)
         end if
         call relaxd(nq,fi,fc,cq,edecr,dismax,icmx,one)
       end if
cpp
       write(*,*) 'In bmat, after relaxation (line 563)'
       edec=edecr*ajoule
       et=etot1+edec
c       write (iout,490) edecr,edec,et
       if(icond.ne.0) write(icond,490) edecr,edec,et
cpp
        write(*,490) edecr,edec,et
       
  490 format (/1x,'Projected energy lowering=',f13.8,' ajoule =',f13.8,
     1' hartree.',3x,' Total energy=',f13.8,' hartree')
      end if
 555     if (dismax.gt.coormx) then
          scaleback=coormx/dismax
          call mult(cq,scaleback,nq)
c         scale back the whole vector
          dismax=coormx
        end if
        call tfer (qq,qq1,nq)
        call add (cq,qq1,nq)
c        call twocol(iout,nq,cq,qq1,
c     1   'GEOMETRY CHANGE AND NEW INTERNAL COORDINATES')
        call twocol(icond,nq,cq,qq1,
     1    'GEOMETRY CHANGE AND NEW INTERNAL COORDINATES')
c       ZERO OUT THE FORCE IN THE DIRECTION OF THE FIXED COORDINATES
        do  740 i=1,nfix
          ii=ifix(i)
          fi(ii)=zero
 740    continue
        call enmax (nq,fi,cq,thract,imax)
c       CALCULATE THRACT = ESTIMATED ENERGY CHANGE PER COORDINATE
c       OR RATHER THE MAXIMUM THEREOF
c        write(iout,340) itopti,etot,etot1,icmy,fmax1,icmx,cq(icmx),
c     1   imax,thract
        if(icond.ne.0) write(icond,340)itopti,etot,etot1,icmy,fmax1,
     1      icmx,cq(icmx),imax,thract
  340   format(/' ROUND',i4,2f15.8,' MAX. FORCE=',i4,f11.7,/,
     1  ' MAX. COORD. CHANGE= ', i3,f11.7,' MAX. ENERGY CHANGE',
     2  i4,e12.5)
        if (iprocs.eq.0.and.lopt(2).gt.0) then
c       IF THERE IS GEOMETRY OPT., WRITE OUT THE EXPECTED GEOM. TOO
#ifdef MOLPRO_NOBACKSPACE
#else
           backspace inpf3
#endif
           write(inpf3,345) itopti,etot1,icmy,fmax1,icmx,cq(icmx),imax,
     1     thract
c         ROUND No., TOTAL E, MAX. FORCE, MAX. CHANGE, MAX. EN. CHANGE
 345       format('ROUND=',i3,' ETOT1=',f15.8,' FORCE=',i3,f10.6,
     1     ' DELTA=',i3,f10.6,' ENER=',i3,e12.5)
           do 350 i=1,nq
             write(inpf3,1170) qq(i),fi(i),i,qq(i)+cq(i)
 350       continue
        end if
cpp
c        write(*,*) 'Line 606 in OVLY50'
CPP
C  Check if the distortion vector is orthogonal to thr null space of G
c        call dgemv_x('n',nq,nq,1.0d0,q(memproj),nq,cq,1,0.0d0,q(memy),1)
c        call add1(cq,-1.0d0,q(memy),nq)
c        write(*,*) '(I-P)*Delta(q)'
c        write(*,*) (q(memy+kk-1),kk=1,nq)
cpp

        
        call distortnew(nek,nq,na,cq,qq,  
     &               qq1,xa,ibcode,ibcontr,b,
     &               q(memginv),shftpt,xy)
cpp
c        write(*,*) 'After calling distort in OVLY50'
        i0=0
          if (itopti.gt.1.and.dismax.gt.chmax.and.itopti.lt.mxopti) then
            icnd=0
          else
            icnd=icond
          end if
c     DO NOT PRINT THE COORDINATES IF AUTOMATIC GEOM.OPT. IS ON
c     EXCEPT IN THE LAST ROUND
      call dispri (i0,iout,icnd,inp2,ipu,na,ia,symb,xy)
      enext=et
      call corlsr(ibase)
      return
      end
c
c-----------------------------------------------------------------------
      subroutine readopx (lopt,na,etot,nfix,ifix,ngeo,lgdi,lfdi,shftpt,
     1  lmur,writ,coormx,inp,inp2,iout,ix)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (nopt=19)
      character*4 wk(nopt),wo
      logical writ,lgdi,lfdi,lmur
      dimension lopt(*),ifix(*)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      data wk/'CARD','FMAT','DISP','OLDF','PUNC',
     1        'BPRI','ANGS','FIXC','FINT','EXFI',
     2        'SHFT','STEP','GDII','FDIA','BFGS',
     3        'AUTO','ROOT','CONT','COND'/
 10   call input(ncol)
      call gets(1,wo,1)
      do 20 i=1,nopt
        if (wk(i).eq.wo) goto 30
 20   continue
      call backspace_inp
      goto 80
 30   call gets (1,wo,1)
      call getf (2,x,1)
      lopt(i)=1
      write (iout,1020) wo,x
 1020 format (1x,a4,2x,'Option is on ',f10.5)
      go to (32,10,34,10,10,
     1       35,10,36,37,39,
     2       42,44,40,60,70,
     3       38,41,10,10),i
c       CARD OPTION  (READ NUCLEAR COORDINATES FROM CARDS)
 32     na=int(x+0.01d0)
        etot=zero
      go to 10
c
c       DISP OPTION  (GENERATE DISPLACEMENTS)
 34     lopt(3)=int(x+0.01d0)
        if (lopt(3).eq.0) lopt(3)=1
      go to 10
c
c       PRINT OPTION - PRINT B MATRIX
 35     writ=.true.
      go  to 10
c
c       FIXC OPTION   (FIX SOME INTERNAL COORDINATES)
 36     nfix=nfix+1
        ifix(nfix)=int(x+0.01d0)
      go to 10
c
c   FINT OPTION - ADD INTERNAL FORCES (OFFSET VALUES)
 37    lopt(9)=int(x+0.01d0)
       if (lopt(9).eq.0) lopt(9)=-1
c      MINUS 1  SIGNALS THAT THERE ARE AS MANY INTERNAL FORCES
c      AS COORDINATES (DEFAULT)
       go to 10
c
c      AUTO OPTION - INCLUDES THE EXFI OPTION
 38    lopt(10)=1
c
c
c       EXFI OPTION (READS THE INTERNAL COORD. FROM AN EXTERNAL FILE)
 39     ix=inp2
      go to 10
c
c      ROOT OPTION -SADDLE POINTS
 41   if(x.eq.2d0) then
       lopt(17)=1
      else
       lopt(17)=0
      end if
      go to 10
c       SHIF OPTION (CHANGES THE SHIFT POINT FOR TORSIONS)
 42     shftpt=x
      go to 10
c
c     OPTI (OPTIMIZE) OPTION. THE PARAMETER IS THE MAXIMUM CHANGE
c     IN INTERNAL COORDINATES, CHMAX
 44     coormx=x
c      MAXIMUM CHANGE PERMITTED IN THE COORDINATES. DEFAULT IS 0.3
      go to 10
c
c     GDIIS OPTION (GEOMETRY DIIS)
 40      ngeo=int(x+1.d-5)
        lgdi=.true.
      go to 10
   60 continue
      lfdi=.true.
      lopt(2)=1
      go to 10
   70 continue
      lmur=.true.
      go to 10
 80   continue
      lfdi=.true.
      lopt(2)=1
      lopt(10)=1
      lopt(16)=1
      end
c-----------------------------------------------------------------------
      subroutine twocol(ifile,n,q,f,title)
c-----------------------------------------------------------------------
      implicit real*8 (a-h,o-z)
c     THIS SUBROUTINE PRINTS INTERNAL COORDINATES AND FORCES, OR
c     SIMILAR QUANTITIES, IN TWO COLUMNS
c     PARAMETERS:
c     IFILE: OUTPUT FILE
c     N:(INPUT) NUMBER OF COORDINATES
c     Q(1:N),F(1:N) (INPUT) COORDINATES & FORCES
c     TITLE (CHARACTER STRING), (INPUT) TITLE
      character*(*) title
      include "common/cmpp"
      dimension q(*),f(*)
      if(ifile.eq.0.or.iprocs.ne.0) return
      write(ifile,100) title
 100  format(/1x,a/)
      nn=(n+1)/2
      do 200 i=1,nn
      ii=i+nn
      if (ii.le.n) then
        write (ifile,300) i,q(i),f(i),ii,q(ii),f(ii)
      else
        write(ifile,300) i,q(i),f(i)
      end if
 200  continue
 300  format(2(1x,i3,2x,f12.8,1x,f12.8,4x))
      end
c----------------------------------------------------------------------
        subroutine enmax(nq,fi,cc,thract,imax)
c----------------------------------------------------------------------
        implicit real*8 (a-h,o-z)
c       ESTIMATE THE LARGEST ENERGY CHANGE
        dimension fi(nq),cc(nq)
        emax=0.0d0
        imax=0
c      IF THE DISPLACEMENT IS VETRY SMALL BUT THERE IS A SIZABLE FORCE
c      THEN DO NOT STOP YET. THE MINIMUM DISPLACEMENT HERE IS
c      0.0015 A. WITH THE LIMIT 1.0E-6, this limits the forces to
c      6.7E-4
        do 100 i=1,nq
          ccc=cc(i)
          if(abs(ccc).lt.0.0015d0) ccc=0.0015d0
          ee=abs(fi(i)*ccc)
          if (ee.gt.emax) then
            emax=ee
            imax=i
          end if
 100    continue
        thract=0.5d0*emax
        end
c----------------------------------------------------------------------
      subroutine clean_zmat(zcrd)
c----------------------------------------------------------------------
      implicit double precision (a-h,o-z)
      logical do_linesearch
      include "common/tapes"
      include "common/big"
      include "common/cgeom"
      include "common/czmat"
      include "common/coptg"
      dimension zcrd(3,*)
      data tol/1.d-5/
c
      icen=0
      do i=1,znuniq
        irow=zuniq(i)
        x=abs(zcrd(1,irow))
        y=abs(zcrd(2,irow))
        z=abs(zcrd(3,irow))
        do j=2,zuniqd(i)
          irow=zequiv(irow)
          x=x+abs(zcrd(1,irow))
          y=y+abs(zcrd(2,irow))
          z=z+abs(zcrd(3,irow))
        end do
        irow=zuniq(i)
        x=x/dble(max(1,zuniqd(i)))
        y=y/dble(max(1,zuniqd(i)))
        z=z/dble(max(1,zuniqd(i)))
        if(x.lt.tol) x=0d0
        if(y.lt.tol) y=0d0
        if(z.lt.tol) z=0d0
        if((abs(x)-abs(zcrd(1,irow))).gt.tol) goto 100
        if((abs(y)-abs(zcrd(2,irow))).gt.tol) goto 100
        if((abs(z)-abs(zcrd(3,irow))).gt.tol) goto 100
        zcrd(1,irow)=sign(x,zcrd(1,irow))
        zcrd(2,irow)=sign(y,zcrd(2,irow))
        zcrd(3,irow)=sign(z,zcrd(3,irow))
        icen=icen+1
        rr(1,icen)=zcrd(1,irow)
        rr(2,icen)=zcrd(2,irow)
        rr(3,icen)=zcrd(3,irow)
        do j=2,zuniqd(i)
          irow=zequiv(irow)
          if((abs(x)-abs(zcrd(1,irow))).gt.tol) goto 100
          if((abs(y)-abs(zcrd(2,irow))).gt.tol) goto 100
          if((abs(z)-abs(zcrd(3,irow))).gt.tol) goto 100
          zcrd(1,irow)=sign(x,zcrd(1,irow))
          zcrd(2,irow)=sign(y,zcrd(2,irow))
          zcrd(3,irow)=sign(z,zcrd(3,irow))
          icen=icen+1
          rr(1,icen)=zcrd(1,irow)
          rr(2,icen)=zcrd(2,irow)
          rr(3,icen)=zcrd(3,irow)
        end do
      end do
      return
100   write(6,*) 'Symmetry broken for Z-matrix row ',irow
      write(6,*) 'x=',zcrd(1,irow),
     >         '  y=',zcrd(2,irow),
     >         '  z=',zcrd(3,irow)
      call error('Broken symmetry','clean_zmat')
      end

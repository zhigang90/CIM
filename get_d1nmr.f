      subroutine get_d1nmr(idft,ax,bl,inx)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /aux/ natom,nucnr(200)
      common /cux/ lden,lfoc,lforc,love,lh01,lhfc,lval,lvec,lrem
c---
c nmr options : only tchf needed here
c
      common /forcdbl/ thre1,thre2,tchf
c---
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
c----------------------------------------------------
      dimension bl(*)
      dimension inx(12,*)
c----------------------------------------------------
      if(lflag(2).eq.1 .or. lflag(3).eq.1) return
c----------------------------------------------------
      call getival('iout',iout)
c----------------------------------------------------
c check bl() status :
c
      call getmem(0,last)
      call retmem(1)
c----------------------------------------------------
      call mmark
c----------------------------------------------------
c memory allocation  :
c
      ntri=ncf*(ncf+1)/2
      ntri3=ntri*3
c
c2006 call getmem(ncf*ncf,lvec)
c2006 call getmem(ncf    ,lval)
      call getmem(ntri3  ,lhfc)
      call getmem(ntri3  ,ldn1)
c2006 call getmem(ntri3  ,lden)
      call getmem(ntri3  ,lfoc)
      call getmem(ntri3  ,love)
      call getmem(2     ,iocc)
c----------------------------------------------------
      call zeroit(bl(lhfc),ntri3)
c----------------------------------------------------
c Read-in :  no of occ orb.,
c            1st-order fock
c            1st-order overlap
c
      np4=4         ! not used in the tempror
      call rea (bl(iocc),1 ,np4,'nocc_rhf')
      nocc=bl(iocc)
c
      irec=1
      call read1mat(61,irec,ntri3,bl(lfoc))
      call read1mat(62,irec,ntri3,bl(love))
c----------------------------------------------------
c     nvirt=ncf-nocc
c     call druwf(ncf,bl(lval+1),bl(lvec+1),nvirt,nocc,iout,0)
c----------------------------------------------------
c for vcd reserve memory for perturbed wf coefficients
c ivcd - vcd flag, lpve - memory pointer
      call tstival('vcd',ivcd)
      If(ivcd.NE.0) call getmem(ncf*nocc*3,lpve)
c
c check memory :
c----------------------------------------------------
      call getmem(0,last)
      call retmem(1)
c----------------------------------------------------
      lrem=last
      last1=last + 6*ncf            ! in dens1 for w1,w2(3*ncf)
c
      if(last1.ge.lcore) then
         call memstat(nreq,nmark,lastadr,memtot,mxmem,ioffset)
         write (iout,2001) last1-ioffset,lcore-ioffset
         call nerror(1,'CPHF ','Memory: needed and available ',
     *               last1-ioffset,lcore-ioffset)
      endif
c
 2001 format(1x,' more memory needed in CPHF : needed=',i5,
     *      '  available=',i5)
c----------------------------------------------------
c solve the  CPHF system of equations :
c
      thresh=tchf
      call setrval('thres',thresh)
      call getival('ncphf',ncphf)
c
      call chfsol(idft,  ax,    ncphf, bl,    inx,
     $            nocc,  thresh,ldn1,  ivcd,  lpve)
c
c----------------------------------------------------
c save D1 (on disk)
c
      call wri(bl(ldn1),ntri3,np4,0,'den1_rhf')
c----------------------------------------------------
c release memory :
c
      call retmark
      call getmem(0,last)
      call retmem(1)
c     last=lrem
c----------------------------------------------------
      end
c==============================================================
C    Coupled-perturbed HARTREE-FOCK solver.
c
      subroutine chfsol(idft,  ax,    ncphf, bl,    inx,
     $                  nocc,  thrs,  ldn1,  ivcd,  lpve)

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      common /aux/ natom,nucnr(200)
      common /cux/ lden,lfoc,lforc,love,lh01,lhfc,lval,lvec,lrem
      common /ganz/ lcore,iov,last,lflag(4),inuc,ibas,na,nbf,nsh,ncf,ncs
     1,nsy(4),nsym,nganz(35),lopt(30)
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      common /nmrint/ ichf,ncache,iprint,nogiao,noacce,ngauge,ntrans
c translation
      common /translx/ trans(3)
c
      dimension bl(*)
      dimension inx(12,*)
c------------------------------------------------
      call getival('iout',iout)
c------------------------------------------------
      ncf3=ncf*3
      ntri=ncf*(ncf+1)/2
      ntr3=3*ntri
c------------------------------------------------
      call mmark
c------------------------------------------------
      call secund(td1cb)
      call elapsec(ed1cb)
c------------------------------------------------
c    Calculate the constant part of the D10
c
c   1. contributions from [ F1(d0,g1) - ei S1 ]
c   2. contributions from the 0.5 D0*S1*D0
c   The total constant part of the D1 stored in bl(lhfc+1)
c
      call d1const(.true.,ncf,ntri,nocc,
     *             bl(lvec),bl(lval),bl(lden),
     *             bl(lfoc),bl(love),bl(lhfc),
     *             bl(lpve),ivcd,bl)
c
c I disapprove of this, disk writing with hard wired unit nos
c passing information from one part of the code to the other
c when same info is also sent via parameters?
c and it messes up the supposedly optimized CPHF start         GM
c
      irec=1
      call save1mat(65,irec,ntr3,bl(lhfc))
c
c     love location not needed anymore
c------------------------------------------------------------
      call secund(td1ce)
      call elapsec(ed1ce)
      td1c=(td1ce-td1cb)/60.d0
      ed1c=(ed1ce-ed1cb)/60.d0
  400 format('Master CPU time for D1 constant part  = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
      write(iout,400) td1c,ed1c
c------------------------------------------------------------
      IF((idft.eq.0) .or. (ax.NE.Zero .and. ncphf.eq.1)) THEN
c        ----------------------------------------
c        For non-zero start in CPHF. Use memory alloc. for love
c
         call zeroit(bl(love),ntr3)
c        ----------------------------------------
         if(ntrans.eq.3) then
c           Calculate the translational part of the D10
c                 [ (Ri - Rj)xT ]*D0
c           to be used as a start for CPHF
c
            call transld1(bl(inuc),trans,inx,bl(lden),
     *                    bl(love),iprint,ncs,ncf,ntri)
c was a bug here, pointer was passed instead of data
c it did not crash - though it was writing to arbitrary memory since
c no one is using TRSL - does not seem to give a good result anyhow
         endif
c        ----------------------------------------
         write(iout,100)
  100    format(/72('-'))
         write(iout,*)' '
         write(iout,*)'            The CPHF Solver '
c        ----------------------------------------
c        Solve the CPHF equation for the D10 matrix
c        starting from D1const or D1const+D1start
c
         cphf_thres=thrs      !  Accuracy for D10
         call solver(idft,ax,inx,nocc,cphf_thres,ncf,ntri,
     *               bl,bl(ldn1),bl(lhfc),bl(love),bl(lval),
     *               bl(lvec),bl(lfoc),ivcd,bl(lpve))
c        ----------------------------------------
      ELSE
c        do not do cphf for non-hybrid dft: D1=D1const
         call dcopy(ntr3,bl(lhfc),1,bl(ldn1),1)
      ENDIF
c
c for the VCD data I will use the old rea/wri system
c and unit 12, nothing uses that unit apparently now
c and considerable data needs to be transferred
c save perturbed magnetic densities to unit 12
      if(ivcd.NE.0) call wri(bl(lpve),ncf*nocc*3,2,0,'c1_magn')
c------------------------------------------------------------
      call retmark
      call getmem(0,last)
      call retmem(1)
c------------------------------------------------
      END
c=================================================================
      subroutine transld1(datnuc,trans,inx,den,dn1,iprint,ncs,ncf,ntri)
c
      implicit real*8 (a-h,o-z)
      character*8 matrix
      common /number/ zero,half,one,two,three,four,five,ten,ten6,tenm8,p
     1i,acc
      dimension datnuc(5,*)
      dimension den(ntri),dn1(ntri,3)
      dimension inx(12,*)
      dimension trans(3)
c--------------------------------
c Calculate the [ (Ri - Rj)xT ]*D0
c--------------------------------
c
      ld1x=ldn1
      ld1y=ld1x+ntri
      ld1z=ld1y+ntri
c
      ifu=0
      do 60 i=1,ncs
         len1=inx(3,i)
         iat=inx(2,i)
         xi=datnuc(2,iat)
         yi=datnuc(3,iat)
         zi=datnuc(4,iat)
         ngci=inx(4,i)
         do 55 igc=0,ngci
c
         jfu=0
         do 50 j=1,i
            len2=inx(3,j)
            jat=inx(2,j)
            xj=datnuc(2,jat)
            yj=datnuc(3,jat)
            zj=datnuc(4,jat)
            ngcj=inx(4,j)
            ngcjx=ngcj
            if(j.eq.i) ngcjx=igc
c
c-------------------------
c calculate (Ri-Rj)xT *D0
c
             xij=xi-xj
             yij=yi-yj
             zij=zi-zj
c
             xfactr=yij*trans(3) - zij*trans(2)
             yfactr=zij*trans(1) - xij*trans(3)
             zfactr=xij*trans(2) - yij*trans(1)
c
             len=len1*len2
c
c----
        do 45 jgc=0,ngcjx
            iff=ifu
            do 40 i1=1,len1
               iff=iff+1
               jff=jfu
               ii=iff*(iff-1)/2
            do 40 j1=1,len2
               jff=jff+1
               if (jff.gt.iff) go to 40
               ij=ii+jff
c
            dn1(ij,1)= dn1(ij,1)-den(ij)*xfactr
            dn1(ij,2)= dn1(ij,2)-den(ij)*yfactr
            dn1(ij,3)= dn1(ij,3)-den(ij)*zfactr
c
   40       continue
            jfu=jfu+len2
   45    continue
   50    continue
         ifu=ifu+len1
   55 continue
   60 continue
c
c--------------------------------
c
      if(iprint.gt.2) then
c
      iout=6
ccccc call druma(bl(lden+1),ncf,iout,'D0-matrx')
c
      write(iout,400)
      call druma(dn1(1,1),ncf,iout,'DX-start')
      call druma(dn1(1,2),ncf,iout,'DY-start')
      call druma(dn1(1,3),ncf,iout,'DZ-start')
c
      endif
c
 400  format(/'*** Translational part of D1 ',8x,'   ***'/)
c
      end
c======================================================================
      subroutine cphf_input(nblocks,labels)
c
c     all output parameters: needed for int_cphf(*)
c     which is a new version of fock10(*).
c

      use memory

      implicit real*8 (a-h,o-z)
c known from post_scf(*):
      common /intbuf/ maxibuf,maxindx
      common /memor1b/ nbl2,nbloks
c---------------------------------------
      maxlabels=maxindx
      nblocks=nbl2*(nbl2+1)/2
      call getmem(maxlabels,labels)
c
      end
c======================================================================
      subroutine check4ld(thres1,thres2,thrs)
      implicit real*8 (a-h,o-z)
c--------------------------------------------------------------------
c inp/out : integral thresholds thres1,thres2 & cphf threshold thrs
c--------------------------------------------------------------------
c check for linear dependencies in the basis set :
c
      call getival('nsym',nsym)
      call getival('iout',iout)
      call tstrval('xlows',lows)
      if(lows.eq.1) then
         call getrval('xlows',xlow)  ! lowest eigenvalue of S
      else
         xlow=0.1d0                  ! when running without scf
      endif
c--------------------------------------------------------------------
      call getival('nsym',nsym)
      rsymm=1.d0
      if(nsym.gt.0) rsymm=dble(nsym)
c--------------------------------------------------------------------
c after testing on annulene/6-31++G :
c Lowest eigenvalue of the overlap matrix  0.626716E-07
c
      if(xlow .lt. 5.0d-6 ) then
         write(iout,*)'  '
         write(iout,*)
     *   '              Near Linear Dependency in Basis Set'
         write(iout,*)
     *   '                sharpen threshold for integrals  '
c
         accur=xlow*1.0d-6
         accur=accur/rsymm
c
         thres1=min(thres1,accur)
c2005    if(thres1.lt.1.0d-14) thres1=1.0d-14
         if(thres1.lt.1.0d-12) thres1=1.0d-12
         thres2=thres1
         if(xlow.lt.1.0d-7)then
            thrs=thrs*5.d0   ! looser cphf threshold
            write(iout,*)
     *      '                loosen  threshold for cphf conv '
         endif
      endif
c--------------------------------------------------------------------
      end
c====================================================================
      subroutine solver(idft,  ax,    inx,   nocc,  thrs,
     *                  ncf,   ntri,  bl,    d1,    d1con,
     *                  r1,    val,   vec,   g1,    ivcd, pve)
c-------------------------------------------------------------------
c This routine solves the cphf equations .
c-------------------------------------------------------------------
c  1. d1    - resulting first-order density
c  2. d1con - constant part of the first-order density matrix
c  3. r1   - starting part of ------''------''-------''-----
c            (zero or our new approach [ (Ri-Rj)xT ]*D0 )
c            later used for residuum r1
c contrary to claims here, gets overwritten with d1con
c this code is really confusing, and comments are missing or misleading
c i spent several days trying to make sense of it MG
c
c  4. val - eigenvalues of the unperturbed solutions
c  5. vec - eigenvectors of the unperturbed solutions
c  6. g1   - the Fock matrix (first-order)
c
c thrs - threshold for CPHF (not integrals)
c idft - dft functional type
c ax - factor to multiply exchange contribution in DFT
c-------------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical oscylates
      character*4 fext1,fext2      ! name extantion for files
      common /nmrint/ ichf,ncache,iprint,nogiao,noacce,ngauge,ntrans
      character*11 scftype
      character*4 where
      common /runtype/ scftype,where
c
      dimension bl(*)
      dimension inx(12,*)
c-------------------------------------------------------------------
      dimension d1(ntri,3),d1con(ntri,3),r1(ntri,3), g1(ntri,3)
      dimension val(*),vec(*)
      dimension pve(ncf,nocc,3)
c-------------------------------------------------------------------
      nfile63=63    !  D1 matrix
      nfile64=64    !  R1 matrix    nothing seems to use this unit??
      nfile65=65    ! constant part of D1
c-------------------------------------------------------------------
      fact=2.0d0     ! factor for dens1_1nmr call (occupancy)
      one=1.0d0
c-------------------------------------------------------------------
      residx=1.d+5
      residy=1.d+5
      residz=1.d+5
c-------------------------------------------------------------------
      call mmark
c-------------------------------------------------------------------
      mxiter=ichf
c-------------------------------------------------------------------
      call getival('iout',iout)
c-------------------------------------------------------------------
      ntr3=3*ntri
      ncf3=3*ncf
c-------------------------
      call getrval('xlvsh',xlvsh)
      call getival('malk',malk)
      if (malk.ne.0)then
          call getival('imalk',imalk)
      else
          imalk=1  ! will not be used
      endif
c-------------------------------------------------------------------
      call secund(tchfb)
      call elapsec(etchfb)
c-------------------------------------------------------------------
c     call getmem(ncf3,lw1)
c     call getmem(ncf3,lw2)
c
      nvirt=ncf-nocc
      call getmem(ncf*ncf,lw1)
      call getmem(ncf*nvirt,lw2)
c---------------------------------------------------------------
c get input for int_cphf(*) ( fock(d1,g0) builder:
c
      call cphf_input(nblocks,labels)
c
c there is one memory allocation by getmem (for labels)
c---------------------------------------------------------------
c thres1 is the final integral threshold
c thres2 is the integral threshold to begin CPHF with
c
      call getrval('thr2',thres1)
      call getrval('thr3',thres2)
c
c  if the calculation is semi-direct, store integrals now
c
      if(scftype.ne.'full-direct') then
        call pstore(nblocks,bl,inx,thres1,bl(labels))
      endif

c---------------------------------------------------------------
      if(ntrans.eq.3) then
c        calculate contributions to D10 from starting-D1
c
         mgo=3    ! full integral threshold
         call pcphf(idft,ax,nblocks,bl,inx,ntri,thres1,mgo,
     *              r1,g1,bl(labels))
c
         call para_next(mgo)
c
c
       do icr=1,3
         call dens1_1nmr1n(fact ,ncf  ,nocc, xlvsh, g1(1,icr),
     1                        vec ,val ,r1(1,icr),bl(lw1),bl(lw2),
     2                        malk,bl(imalk))
       enddo
c
         call daxpy(ntr3,one,r1,1,d1con,1)
      endif
c---------------------------------------------------------------
      mgo=1    ! loose integral threshold
c---------------------------------------------------------------
      if(thres2.eq.thres1) mgo=3    !  final integral threshold
c---------------------------------------------------------------
      write(iout,*) ' '
      write(iout,150) thrs
  150 format('threshold for 1st-order density = ',1pe11.4)
c 150 format('threshold for CPHF convergence  = ',1pe11.4)
c
      thres_i=thres2      ! loose int. thresh. at the begining
      write(iout,151)
  151 format(/'CPHF iter',2x,'res-x',7x,'res-y',7x,'res-z',7x,
     *        'cpu   elapsed  oscillating')
c
      if(mgo.lt.3) then
         write(iout,152) thres_i
      else
         write(iout,154) thres_i
      endif
  152 format('    Loose integral threshold = ',1pe11.4)
  154 format('    Final integral threshold = ',1pe11.4)
c
      call f_lush(iout)
c---------------------------------------------------------------
c last integral threshold used :
c
      thres_last=thres_i
      liter_thre=0     ! counting iterations with the same int. threshold
c---------------------------------------------------------------
c
c           begining of the solution ofthe CPHF equation
c
c---------------------------------------------------------------
      IAT=1
c---------------------------------------------------------------
c The first D1 matrices = D1const
c
      call read1mat(nfile65,iat,ntr3,r1 ) ! r0=d1const
c      call read1mat(nfile65,iat,ntri*3,g1 ) ! r0=d1const
      call dcopy(ntr3,r1,1,g1,1)
c---------------------------------------------------------------
      call zeroit(d1,ntr3)
c---------------------------------------------------------------
      call getival('reset',ireset)
      icycle=0
 4321 continue
c---------------------------------------------------------------
c if iteration is reset, this part breaks as r1 and g1 are
c residuals and fock matrices stored in place of const part of d1
c (gm)
      call save1mat(nfile63,iat,ntr3,r1 )   ! save d0=r0
c commented out next line as nfile64 is not used (gm)
c      call save1mat(nfile64,iat,ntr3,r1 )   ! save r0=r0
      call save1mat(nfile65,iat,ntr3,g1)    ! save new r0=r0
c---------------------------------------------------------------
      DO LITER=1,MXITER
         icycle=icycle+1
c
         call get_fext(liter,fext1,fext2)
c                            name.fext1  - files for rl=l(r)
c                            name.fext2  - files for ro=O(rl)
c
c        Calculate the D = C + SUMi [ Li(C)]
c
         call secund(citerb)
         call elapsec(eiterb)
c
c        calculate G(D1,g0) in g1
c
         call pcphf(idft,ax,nblocks,bl,inx,ntri,thres_i,mgo,
     *               r1,g1,bl(labels))
c
c
c
       do icr=1,3
         call dens1_1nmr1n(fact ,ncf  ,nocc, xlvsh, g1(1,icr),
     1                        vec ,val ,r1(1,icr),bl(lw1),bl(lw2),
     2                        malk,bl(imalk))
       enddo
c
c
            call file4cphf_o(ntri*3,iat,fext1,r1,'write') ! save   rl1=L(r1)
c
            IF(NOACCE.EQ.1) THEN
               call calc_d10(nfile65,liter,ntri,iat,d1, r1) ! out: d1,r1
            ELSE
               call make_r_orto(nfile63,liter,ntri,iat,d1con, r1) ! i/o l(r) ort
ccccc          call chec_r_orto(nfile63,liter,ntri,iat,d1con, r1)
               call file4cphf_o(ntri*3,iat,fext2,r1,'write') ! save  ro1=O(rl1)
c              ......................................
               if(liter.ge.2) then
c                calculate current solution d1
c
                 call getmem(liter*liter,iamat)
                 call getmem(liter*liter,ibmat)
                 call getmem(liter+1    ,iwvec)
                 call getmem(liter      ,ilvec)
                 call getmem(liter      ,imvec)
                 call getmem(liter      ,icoef)
c                use d1con & g1 as scratch arrays
c strange comment as g1 (fock mx derivative) is not used here (gm)
                 call calc_d1r(nfile63,liter,ntri,iat,
     *                        bl(iamat),bl(ibmat),bl(iwvec),
     *                        bl(ilvec),bl(imvec),bl(icoef),
     *                        d1con,r1,  d1)  ! out :current d1  new r1
                 call retmem(6)
               endif    !  (liter.ge.2) then
c              ......................................
               call file4cphf_o(ntri*3,iat,fext2,r1,'read') ! get ro1 for enx
            ENDIF    !   IF(NOACC.EQ.1) THEN
c
         call secund(citere)
         call elapsec(eitere)
         cpuit=citere-citerb
         elait=eitere-eiterb
c
         call cphf_enx(r1,thrs,ntri,mgo,liter,cpuit,elait,lend,
     *                 errmax, residx,residy,residz,oscylates)
c
         IF(LEND.EQ.1) EXIT
c
c the following lines read back what was just written out (mg)
c        if(noacce.eq.1) then
c           call file4cphf_o(ntri*3,iat,fext1,r1,'read')  ! get rl1
c        else
ccccc       call file4cphf_o(ntri*3,iat,fext2,r1,'read') ! get  ro1=O(rl1)
c        endif
c
         call cphf_int_thr(mgo,liter_thre, errmax,oscylates,
     *                         thres_i,thres_last,thres1,thrs)
c
         if(icycle.eq.mxiter) exit
c
         call para_next(mgo)
         if(liter.eq.ireset .and. noacce.eq.0) go to 4321
c
      ENDDO    ! end of iterations
c----------------------------------------------
      call file4cphf_c(ntri*3,liter)
c----------------------------------------------
c     end of iterations
c----------------------------------------------
c DO ONE MORE CALCULATION OF D1 USING FULL D1 :
c
      if(liter.eq.icycle) then
         call para_next(mgo)
         if(liter.eq.1) then  ! 2002 lublin
            call read1mat(nfile65,iat,ntri*3,d1)  ! d1=d1const
         endif
         call pcphf(idft,ax,nblocks,bl,inx,ntri,thres_i,mgo,
     *              d1,g1,bl(labels))
c
c
       do icr=1,3
         call dens1_1nmr1n(fact ,ncf  ,nocc, xlvsh, g1(1,icr),
     1                        vec ,val ,d1(1,icr),bl(lw1),bl(lw2),
     2                        malk,bl(imalk))
c -- add contribution to perturbed coefficient matrix for VCD
         if(ivcd.ne.0)
     $     call daxpy(ncf*nocc,one,bl(lw2),1,pve(1,1,icr),1)
       enddo
c
c
         call read1mat(nfile65,iat,ntri*3,r1 ) ! r1=d1const
         call daxpy(ntri*3,one,r1,1,d1,1)        ! BLAS!!!
      endif
c----------------------------------------------
      if(iprint.gt.2) then
        call druma(d1(1,1),ncf,iout,'DX-final')
        call druma(d1(1,2),ncf,iout,'Dy-final')
        call druma(d1(1,3),ncf,iout,'Dz-final')
      endif
c
      if(lend.eq.0) then
         write(iout,440)
         write(iout,441) mxiter,thrs
         write(iout,442)
      endif
c
      call getival('indisk',indisk)
      if(indisk.gt.0) then
         call clos4int()
      endif
c----------------------------------------------------------------------
      call secund(tchfe)
      call elapsec(etchfe)
      tchf=(tchfe-tchfb)/60.0d0
      elaps=(etchfe-etchfb)/60.0d0
      write(iout,430) tchf,elaps
      call f_lush(iout)
c----------------------------------------------------------------------
  430 format(/'Master CPU time for   CPHF  solver    = '
     *,f8.2,'  Elapsed = ',f8.2,' min')
  440 FORMAT(/
     *'**   YOU DID NOT GET THE CONVERGENCE FOR CHF SOLUTION   **')
  441 FORMAT(
     *'**       mxiter =',i4,'  threshold = ',D12.6,'       **')
  442 FORMAT(
     *'**                        S O R R Y                     **')
c----------------------------------------------------------------------
      call retmark
c
      end
c====================================================================
      subroutine dens1_1nmr1n(fact ,ncf  ,nocc, xlv,  fock,
     1                        cvec ,eval ,d1  , w1 ,  w2,
     2                        malk,xmalk)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates a part of the first-order density matrix :
c
c D1 = SUM(jocc,avir){Cj(T)*F(D1,g0)*Ca/(ej-ea-xlv)*[Cj*Ca(T)-Ca*Cj(T)]}
c
c  (T) = transpose
c-----------------------------------------------------------------------
c  INTENT(IN)
c  fact    = occupancy : 2 for RHF, 1 for UHF
c  ncf     = number of contracted basis functions
c  nocc    = number of occupied orbitals
c  xlv     = level shift  (see below)
c  fock    = first-order Fock matrix ,F(D1,g0) above, in triangular form
c  cvec    = MO coefficients, C, the first nocc columns are occupied
c  eval    = orbital energies e above
c  malk    = Malkin's correction flag
c  xmalk   = Malkin's corrections
c  INTENT(OUT)
c  d1      = the above part of the first-order density
c  STORAGE
c  w1      = a working array ncf**2 long
c  w2      = a working array ncf*nvirt long  (nvirt=ncf-nocc)
c---------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.0d0,two=2.d0)
      dimension fock(*),cvec(ncf,ncf),eval(*),d1(*)
      dimension xmalk(ncf-nocc,nocc)
      dimension w1(*),w2(*)
c
      nvirt=ncf-nocc
c  expand the Fock matrix to quadratic
      call quad(fock,w1,-one,ncf)
c  W2=Cocc(T)*F
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w2, nocc)
c  W1=Cocc(T)*F*Cvirt=W2*Cvirt  nocc x nvirt matrix
      call dgemm('n','n', nocc, nvirt,ncf,
     1            one, w2, nocc, cvec(1,nocc+1), ncf,
     2            zero,w1, nocc)
c
c  Now scale W1=Cocc(T)FCvirt with 1.0/(e(i)-xlv-e(a)-xmalk(a,i))
c
      if(malk.eq.0) then
         call scalebydenom(nocc,nvirt,eval,eval(nocc+1),xlv,w1)
      else
         call scalebydeno1(nocc,nvirt,eval,eval(nocc+1),xlv,xmalk,w1)
      endif
c
C Calculate Cvirt*W1(T).
c Result, W2(ncf,nocc) is the perturbed coefficient matrix!!
      call dgemm('n','t', ncf, nocc, nvirt,
     1            one, cvec(1,nocc+1), ncf, w1, nocc,
     2            zero,w2, ncf)
c Calculate W2*Cocc(T), factor of -1 since this is the transpose
c Multiply with the factor (occupancy) as well
      call dgemm('n','t',ncf,ncf,nocc,
     1            -fact,w2,ncf,cvec,ncf,
     2            zero, w1,ncf)
c  Result, in quadratic form, is in W1.
c  Add transpose to W1 and transfer it to the triangular array d1
c
      call symtria1(ncf,w1,d1)
c
c     call printd1('new',ncf,d1)
      end
c======================================================================
      subroutine scalebydeno1(nocc,nvirt,eocc,evirt,xlv,xmalk,f)
c  This routine divides element (i,a) of the matrix F by
c  (eocc(i)-eocc(a)-xlv); F(i,a)=F(i,a)/(eocc(i)-eocc(a)-xlv)
c
c  Arguments:
c  INTENT(IN)
c  nocc     = number of occupied orbitals, number of rows of F
c  nvirt    = number of virtual orbitals, number of columns of F
c  eocc     = occupied orbital energies
c  evirt    = virtual orbital energies
c  INTENT(INOUT)
c  f        = Fock matrix (occupied x virtual part in MO basis)
      implicit real*8 (a-h,o-z)
      integer a
      dimension f(nocc,nvirt),eocc(nocc),evirt(nvirt),xmalk(nvirt,nocc)
      do a=1,nvirt
        xx=evirt(a)+xlv
        do i=1,nocc
          yy=eocc(i) - xx - xmalk(a,i)
          f(i,a)=f(i,a)/yy
        end do
      end do
      end
c======================================================================
      subroutine symtria1(n,a,b)
c  This routine adds A-A(T) to the symmetrical matrix B stored as
c  the upper triangle row-wise.
c  Arguments:
c  INTENT(IN)
c  n - dimension of the square matrix A
c  A(n,n) - square matrix
c  INTENT(OUT)
c  B(1:n*(n+1)/2): triangular matrix, it i,j element is A(i,j)+A(j,i)
      implicit real*8 (a-h,o-z)
      dimension a(n,n),b(*)
      ij=0
      do i=1,n
        do j=1,i
          ij=ij+1
          b(ij)=a(i,j)-a(j,i)
        end do
      end do
      end
c======================================================================
      subroutine d1const(rhf,   ncf,   ntri,  nocc,  vec,
     *                   val,   d0,    f1,    s1,    d1c,
     *                   pve,   ivcd,  bl)
c-----------------------------------------------------------------
c Calculates the constant part of the D1 matrix for NMR :
c-----------------------------------------------------------------
c INPUT :
c
c rhf      - rhf or uhf flag
c ncf      - number of basis function
c ntri     -  ncf*(ncf+1)/2
c nocc     - number of occupied orbitals
c vec      - eigen vectors
c val      - eigen values
c d0       - unperturbed density
c f1       - Ist-order fock matrix
c s1       - Ist-order overlap matrix
c ivcd     - is VCD calculated?
c bl       - storage for everything
c OUTPUT :
c
c d1c      - constant part of D1 :
c            D1const=-1/2*DS1D+
c            2*Sum(o,v)(eo-ev)^-1*Co+(h1+G(d0,g1)-eo*S1)Cv*[CoCv+ - CvCo+]
c pve      - perturbed coefficients (only for VCD)
c-----------------------------------------------------------------

      use memory, block => bl

      implicit real*8 (a-h,o-z)
      logical rhf
      dimension bl(*)
      dimension d0(*), vec(*), val(*)
      dimension s1(ntri,3), f1(ntri,3)
      dimension d1c(ntri,3)                      ! output
      dimension pve(ncf,nocc,3)
c--------------------------------------------------------------
c    Calculate the constant part of D10
c   ----------------------------------------
c 1. contributions from [ F1(d0,g1) - ei S1 ]
c 2. contributions from the 0.5 D0*S1*D0
c--------------------------------------------------------------
      call getrval('xlvsh',xlvsh)
      call getival('malk',malk)
      if (malk.ne.0) call getival('imalk',imalk)
c--------------------------------------------------------------
      fact1=2.0d0     ! for dens1_part1 like occupancy number
      fact2=-0.5d0     ! for dS1d
      if(.not.rhf) then
        fact1=1.0d0
        fact2=-1.0d0
      endif
c--------------------------------------------------------------
c 1. calculate contributions from [ F1(d0,g1) - ei S1 ]
c
c one direction at the time :
c
      call getmem(ncf**2,lw1)
      nvirt=ncf-nocc
      memvirt = ncf*MAX(nvirt,nocc)
      call getmem(memvirt,lw2)
      call getmem(ncf*nocc,lw3)
c
      do icr=1,3
         call dens1_1nmr2n(fact1,ncf, nocc, xlvsh, f1(1,icr),
     *                    s1(1,icr),vec, val, d1c(1,icr), bl(lw1),
     *                    bl(lw2),bl(lw3), malk,bl(imalk))
         if(ivcd.NE.0)
c part of the der. coeff. mx related to [ F1(d0,g1) - ei S1 ]
     $     call dcopy(ncf*nocc,bl(lw2),1,pve(1,1,icr),1)
c      call printd1('new',ncf,d1c(1,icr))
      enddo
c
      call retmem(3)
c
c 2.  Calculate  0.5*D0*S1*D0 and add the result
c     to the constant part of the D10.
c     This is the total constsnt part of D10.
      call getmem(ncf*ncf,lw0)
      call getmem(ncf*ncf,lw1)
      call getmem(ncf*ncf,lw2)
c
      do icr=1,3
       call ds1d_1d(bl(lw0),bl(lw1),bl(lw2),ncf,d0,s1(1,icr),d1c(1,icr))
       if(ivcd.NE.0)
c 1/8 D*S1*Co part of the pert. coeff. mx.
     $   call dgemm('n','n',ncf,nocc,ncf,0.25d0,bl(lw2),ncf,
     $              vec,ncf,1.0d0,pve(1,1,icr),ncf)
      enddo
c
      call retmem(3)
c---------------------------------------------------------------
      end
c======================================================================
      subroutine dens1_1nmr2n(fact ,ncf,nocc, xlv,  fock,
     1                       smat, cvec ,eval ,d1, w1,
     2                       w2,   w3,   malk, xmalk)
c-----------------------------------------------------------------------
c one direction at the time
c-----------------------------------------------------------------------
c This routine calculates the constant part of the first-order density matrix :
c
c D1 = SUM(jocc,avir){Cj(T)*[F(D0,g1)-ej*S]*Ca/(ej-ea-xlv)*[Cj*Ca(T)-Ca*Cj(T)]}
c
c  It is calculated here as
c
cD1(p,q)=SUM(j,a)[C(T)*F-Eocc*C(T)S]*Cvirt](j,a)/(ej-ea-xlv)*[Cpj*Cqa-Cpa*Cqj]
c
c  (T) = transpose
c-----------------------------------------------------------------------
c  INTENT(IN)
c  fact    = occupancy : 2 for RHF, 1 for UHF
c  ncf     = number of contracted basis functions
c  nocc    = number of occupied orbitals
c  xlv     = level shift  (should not be used here, set it to zero)
c  fock    = first-order Fock matrix ,F(D0,g1) above, in triangular form
c  smat    = first-order AO overlap matrix in triangular form
c  cvec    = MO coefficients, C, the first nocc columns are occupied
c  eval    = orbital energies e above
c  INTENT(OUT)
c  d1      = the above part of the first-order density
c  STORAGE
c  w1      = a working array ncf**2 long
c  w2      = a working array ncf*nvirt long
c  w3      = a working array ncf*nocc long
c---------------------------------------------------
      implicit real*8 (a-h,o-z)
      parameter (zero=0.d0,one=1.0d0,two=2.d0)
      dimension fock(*),cvec(ncf,ncf),eval(*),d1(*)
      dimension w1(ncf,ncf),w2(*),w3(nocc,ncf)
      dimension xmalk(ncf-nocc,nocc)
c
      nvirt=ncf-nocc
c      ivirtst=nocc*ncf+1
c  expand the Fock matrix to quadratic
      call quad(fock,w1,-one,ncf)
c  W2=Cocc(T)*F
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w2, nocc)
c  expand the overlap matrix to quadratic
      call quad(smat,w1,-one,ncf)
c  W3=Cocc(T)*S
      call dgemm('t','n',nocc,ncf,ncf,
     1           one, cvec, ncf, w1, ncf,
     2           zero, w3, nocc)
c  Multiply the rows of W3 by the occupied orbital energies and subtract
c  them from W2
      call RowMultiply(nocc,ncf,eval,w3,w2)
c  Build W1=[Cocc(T)*F-Eocc*Cocc(T)*S]Cvirt=W2*Cvirt
      call dgemm('n','n', nocc, nvirt,ncf,
     1            one, w2, nocc, cvec(1,nocc+1), ncf,
     2            zero,w1, nocc)
c  Now scale W1=Cocc(T)FCvirt with 1.0/(e(i)-e(a)-xlv)
c
      if(malk.eq.0) then
         call scalebydenom(nocc,nvirt,eval,eval(nocc+1),xlv,w1)
      else
         call scalebydeno1(nocc,nvirt,eval,eval(nocc+1),xlv,xmalk,w1)
      endif
C Calculate Cvirt*W1(T). Result, W2 is part of the perturbed C mx
      call dgemm('n','t', ncf, nocc, nvirt,
     1            one, cvec(1,nocc+1), ncf, w1, nocc,
     2            zero,w2, ncf)
c Multiply with the factor (occupancy)
c Calculate W2*Cocc(T), factor of -1 !
      call dgemm('n','t',ncf,ncf,nocc,
     1            -fact,w2,ncf,cvec,ncf,
     2            zero, w1,ncf)
c  Result, in quadratic form, is in W1.
c  Add transpose to  W1 and transfer it to the triangular array d1
      call symtria1(ncf,w1,d1)
c     call printd1('new',ncf,d1)
      end
c======================================================================
      subroutine ds1d_1d(w0,w1,w2,ncf,d0,s1,d1c)
      implicit real*8 (a-h,o-z)
c
      dimension d0(*),s1(*)
      dimension w0(ncf,ncf),w1(ncf,ncf),w2(ncf,ncf)
      dimension d1c(*)         ! inp/out
c
      data zero,half,one /0.d0, 0.5d0, 1.0d0/
c
c  expand the dens and overlap1 matrices to quadratic
c
      call quad(d0,w0, one,ncf)
      call quad(s1,w1,-one,ncf)
c  calculate W2=D0*S1
      call dgemm('n','n',ncf,ncf,ncf,
     1           one,w0 , ncf, w1, ncf,
     2           zero, w2, ncf)
c  calculate w1=0.5*W2*D0=0.5*D0*S1*D0
      call dgemm('n','n',ncf,ncf,ncf,
     1           half,w2,ncf,w0,ncf,
     2           zero, w1, ncf)
c
      ij=0
      do i=1,ncf
         do j=1,i
            ij=ij+1
            d1c(ij)=d1c(ij)+w1(j,i)
         enddo
c        ii=i*(i+1)/2
c        d1c(ii)=zero
      enddo
c
      end

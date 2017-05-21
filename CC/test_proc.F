#include "maxi.h"
C==============================================================================
      function imap(ij)
c This function returns the i index of triangular matrix, i>=j.
c Not very efficient but rarely used
      if (ij.le.0) stop 'Error in imap function'
      do i=1, (CORR_ORB)
      j=ij-i*(i-1)/2
      if (j.gt.0 .and. j.le.i) then
        imap=i
        return
      endif
      enddo
      stop 'Error in imap function'
      end
C==============================================================================
      function jmap(ij)
c This function returns the j index of triangular matrix, i>=j.
c Not very efficient but rarely used
      if (ij.le.0) stop 'Error in jmap function'
      do i=1, (CORR_ORB)
      j=ij-i*(i-1)/2
      if (j.gt.0 .and. j.le.i) then
        jmap=j
        return
      endif
      enddo
      stop 'Error in jmap function'
      end
C==============================================================================
      subroutine fill1(ncf,xmatrix)
c It fills given xmatrix of dimension ncf*ncf with fixed large numer.
      implicit real*8 (a-h,o-z)
      dimension xmatrix(ncf,ncf)
      do i=1,ncf
      do j=1,ncf
      xmatrix(j,i)=1000.d0
      enddo
      enddo
      end
c==============================================================================
      subroutine FockBuilder(irecadrc,irecadrx,nval,   npairs,  thresh,
     *                       ndiskc,  ndiskx,  lbinc,  lbinx,   byt8,
     *                       nmo,     vorb,    ifockAO,ifockMO, icoreAO,
     *                       nbf)
      use memory
c This subroutine produces Fock matrix in AO and MO form and calculates
c core fock matrix by subtracting all two-electron components from
c the correlating orbitals from full Fock matrix
c The fock matrix is read from the disk, always in AO form.
c INPUT:
c irecadrc - integer array containing (single) the last record number for given
c             ij for coulomb integrals, the rest of record numbers will be read 
c             from disk, irecadrc(npairs) in this case.
c            in parallel case is irecadrc(npairs*nslv) it contains the
c             number of records written in each slave, the reords are
c             numbered there from 1 to irecadr(ij,islave) where islave is
c             the slave number.
c irecadrx - see irecadrc, exchange integrals.
c nval     - number of correlating orbitals
c npairs   - nval*(nval+1)/2
c thresh   - integral threshold
c ndiskc   - the disk number with coulomb  integrals
c ndiskx   - the disk number with exchange integrals
c lbinc    - the bin length for coulomb  integrals
c lbinx    - the bin length for exchange integrals
c byt8     - logical variable indicating how integrals are stored: 8
c            byte of 5 byte format (bins only are influenced by this
c            parameter
c nmo      - number of occupied orbitals
c vorb     - logical variable indicating how integrals are calculated: 
c            in MO or in AO (vorb=.true. or vorb=.false. Influences the 
c            output - the form of ifockAO and icoreAO, for vorb=.true.
c            ifockAO becomes identical to virt x virt subblock of ifockMO
c INPUT/OUTPUT:
c ifockAO  - The bl pointer to AO Fock Matrix or virt x virt transformed
c            subblock of it - see above discussion
c OUTPUT:
c ifockMO  - full Fock matrix in MO form
c icoreAO  - The core-fock matrix in AO or MO (virt x virt) form
      implicit real*8 (a-h,o-z)
      dimension irecadrc(npairs),irecadrx(npairs)
      logical byt8,vorb
c
      call getival('ncf',ncf)
c
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call getival('SymFunPr',ifp)
      nsym=igetival('nsym')
      call check_operator_symmetry(nsym,bl(ifp),ncf,bl(ifockAO))
c
      call dynamic_matdef('workAO','q',ncf,ncf)
      call matconn('fock','q',ncf,ncf,ifockAO)
      call matcopy('fock','workAO')
      call matdisc('fock')
      call matconn('fock','q',idimen,idimen,ifockAO)
      call matconn('fockMO','q',ncf,ncf,ifockMO)
c
      if (vorb) then
        call matsimtr('workAO','genvirt','fock')
      endif
      call matsimtr('workAO','cano','fockMO')
      if (lbinc.ne.0) then
      call matconn('corecore','q',idimen,idimen,icoreAO)
      call dynamic_matdef('AOcore','q',ncf,ncf)
      call matcopy('workAO','AOcore')
      call dynamic_matdef('xmat','q',ncf,ncf)
      ixmat=mataddr('xmat')
      do i=1,nval   ! This is subtraction of correlated orbitals!
      j=i
      call ExtrOne(i,    j,     irecadrc,   npairs,    ndiskc,
     *             ncf,  lbinc, thresh,     byt8,      'c',
     *             'ao', nmo,   vorb,       bl(ixmat))
      call matadd1('xmat',-2.0d0,'AOcore')
      call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *             ncf,  lbinx, thresh,     byt8,      'x',
     *             'ao', nmo,   vorb,       bl(ixmat))
      call matadd1('xmat',+1.0d0,'AOcore')
      enddo
      if (vorb) then
        call matsimtr('AOcore','genvirt','corecore')
      else
        call matcopy('AOcore','corecore')
      endif
      endif
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine ExtrOne(iorb,  jorb,  irec,    npairs,  ndisk,
     *                   ncf,   lbin,  thresh,  byt8,    typ,   
     *                   AOMO,  nmo,   vorb,    xmat)
      use memory
c One of the oldest and most frequently used procedures in CCSD. It
c reads the integrals from the disk and puts them in xmat. Depending on 
c options it reads integrals from the bin file or integral file, where
c integrals are stored in matrix form. This subroutine is responsible
c for all integral caching too. Formerly the transformation AO->MO was
c performed here, but this was moved to GenCoulExInt during parallelization
c INPUT:
c iorb     - each matrix is produced for i,j pair. This is i index
c jorb     -                                               j index
c irec     - the array with records data, see GenCoulExInt FockBuilder descr.
c npairs   - number of correlating orbitals unique pairs.
c ndisk    - the disk no with integrals stored in bins form.
c ncf      - number of basis function
c lbin     - bin legth
c thresh   - integral threshold
c byt8     - .true. one integral (in bin only) occupy 5 bytes, .false. - 
c            occupy 8 bytes, the real*8 number in the latter case.
c typ      - character string indicating the kind of integral, may be:
c            'x'  - exchange operator integrals
c            'c'  - coulomb  operator integrals
c            'e'  - external exchange operator integrals
c            'tx' - three external indices exchange operator integals
c            'tx' - three external indices coulomb  operator integals
c AOMO     - which form is requested, AO from or MO form? Identical in
c            case vorb=.false.
c            'ao' - get me AO form of integrals (ncf x ncf)
c            'mo' - get me MO form of integrals (virt x virt)
c nmo      - number of occupied HF orbitals
c vorb     - logical variable, .true. means that program runs in the MO
c            mode, .false. in the AO mode
c COMMON BLOCK INPUT:
c /extrone_tmp/ common:
c af       - logical indicating if we use Array Files here
c islaves  - number of slaves working. Should be "1" in single mode - 
c            it is always used!
c ndisktr* - there are disk unit numbers with transformed integrals.
c mygid    - the GID number, used in parallel when cache is on
c /GlobalCCSD/ common:
c cache    - logical, if .true. do local disk caching in parallel.
c            Should be set to false in single mode.
c OUTPUT:
c xmat     - AO or MO form of integrals requested. Note the matrix has
c            different size depending on vorb mode and AOMO requested
      implicit none
      integer iorb,jorb,npairs,ndisk,ncf,lbin,nmo
      integer irec(*),ibin4(lbin)
      integer*1 ibin1(lbin)
      real*8    bin8(lbin), thresh, xmat(ncf,ncf)
      integer*2 indexx(2,lbin)
      character*(*) typ ! may be x,c,e,tx,tc (exchange, coulomb,EEO)
      character*(*) AOMO ! may be mo,ao, if vorb=.false. - ignored
      logical byt8,vorb
c local:
      real*8 xint,one,dblinv,d1max,dblmax,rnpairs
c
      integer ij,irecord,jrecord,i,j,icount,nval
      integer idimen,idimen2,ncf2,iafindex(4),islaves,istatus
      character*120 message
      parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
      parameter (dblinv=one/dblmax,d1max=128.0d0)
      integer irecskip,itmp,imap,jmap,mygid
      integer ndisktr,ndisktrc,ndisktrx,ndisktre,ndisktrtx,ndisktrtc
      integer ndisktrtt,nbf
      logical af
      common /extrone_tmp/ af,islaves,ndisktrc,ndisktrx,ndisktre,
     *                     ndisktrtx,ndisktrtc,ndisktrtt,mygid,nbf
      integer len,len1,lrec
      integer itmpEEO1,itmpEEO2,mataddr
c Fragile:
      character*256 scrfile,filename
      character*2 ch_mygid
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      integer ns,iscr,iprim,jprim,ijprim,ii,jj,mu,nu,iii,list(3,10)
      integer irecij,iadd
      logical trans1,signum,ex3
      real*8 diff
c
c
#define SYMMETRY
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
      idimen2=idimen*idimen
      ncf2=ncf*ncf
      rnpairs=npairs
      nval=(-1.0d0+sqrt(8.0d0*npairs+1))/2.d0+0.001d0
      irecskip=npairs
      ex3=.false.
      if (typ.eq.'tx'.or.typ.eq.'tc') irecskip=nval*nval
c
      if      (typ.eq.'x') then
        ndisktr=ndisktrx
      else if (typ.eq.'c') then
        ndisktr=ndisktrc
      else if (typ.eq.'e') then
        ndisktr=ndisktre
      else if (typ.eq.'tx') then
        ndisktr=ndisktrtx
      else if (typ.eq.'tc') then
        ndisktr=ndisktrtc
      else if (typ.eq.'tt') then
        ndisktr=ndisktrtt
        ex3=.true.
      endif
c
      if (typ.eq.'x'.or.typ.eq.'e'.or.typ.eq.'c') then
c     ij - lower triangle i>=j
      if (iorb.ge.jorb) then 
        ij=iorb*(iorb-1)/2+jorb
      else
        ij=jorb*(jorb-1)/2+iorb
      endif
      if (ij.gt.npairs) STOP 'ERROR 1 in ExtrOne subroutine'
      else
        ij=(iorb-1)*nval+jorb
      if (typ.ne.'tt'.and.ij.gt.nval*nval) then 
              print *,'ij,nval,i,j:   ',ij,nval,iorb,jorb
              STOP 'ERROR 1a in ExtrOne subroutine'
       else if (typ.eq.'tt'.and.ij.gt.idimen*nval) then 
              print *,'ij,nval,i,j:   ',ij,nval,iorb,jorb
              STOP 'ERROR 2a in ExtrOne subroutine'
       endif 
      endif
c
c
#ifdef SYMMETRY
      call pair_searcher(iorb, jorb,   iprim, jprim, ijprim,
     *                   ns,   trans1, signum)
#else
      ijprim=ij
#endif
      if (.not.vorb.or.AOMO.eq.'mo') then
        if (.not.vorb.and.typ.eq.'tt') STOP '3ext works with MO only'
        if (ns.ne.0.and.typ.eq.'tt') STOP '3ext works without symmetry'
        if (af) then
          if      (typ.eq.'x'.or.typ.eq.'c') then
            call fafread(ndisktr,xmat,8,idimen*idimen,1,ijprim,istatus)
          if (istatus.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
#ifdef SYMMETRY
            call symm_transform(idimen,ns,trans1,signum,xmat)
#endif
          else if (typ.eq.'e') then ! EEO needs special treatment:
            if (vorb) then
              call dynamic_matdef('tmpEEO1','q',ncf,ncf)
              itmpEEO1=mataddr('tmpEEO1')
              call dynamic_matdef('tmpEEO2','q',idimen,idimen)
              itmpEEO2=mataddr('tmpEEO2')
           call fafread(ndisktr,bl(itmpEEO1),8,ncf*ncf,1,ijprim,istatus)
          if (istatus.lt.8*ncf*ncf) then
              ERROR
              call flush(6)
          endif
#ifdef SYMMETRY
              call symm_transform_ao(ncf,ns,trans1,signum,bl(itmpEEO1))
#endif
c
              call matsimtr('tmpEEO1','genvirt','tmpEEO2')
              call tfer(bl(itmpEEO2),xmat,idimen*idimen)
              call dynamic_matrem('tmpEEO2')
              call dynamic_matrem('tmpEEO1')
            else
              call fafread(ndisktr,xmat,8,ncf*ncf,1,ijprim,istatus)
          if (istatus.lt.8*ncf*ncf) then
              ERROR
              call flush(6)
          endif
#ifdef SYMMETRY
              call symm_transform_ao(ncf,ns,trans1,signum,xmat)
#endif
c
            endif
          else ! (typ.eq.'x'.or.typ.eq.'c') then (== 'tx' or 'tc' or 'tt')
#ifdef SYMMETRY
          if (typ.ne.'tt') then
            call full_pair_searcher(iorb,  jorb,  iprim, jprim, ijprim,
     *                              ns, signum)
          else
            ijprim=ij
          endif
#else
      ijprim=ij
#endif
            call fafread(ndisktr,xmat,8,idimen*idimen,1,ijprim,istatus)
            if (istatus.lt.idimen*idimen*8) then
                ERROR
                print *,'idis,ij,idime,info: ',ndisktr,ij,idimen,istatus
                call flush(6)
            endif
#ifdef SYMMETRY
            if (ns.ne.0.and.typ.ne.'tt') then
              call full_symm_transform(idimen,ns,signum,xmat)
            endif
#endif
          endif
        else !  if (af) then
          if (typ.eq.'e') then ! EEO needs special treatment:
            if (vorb) then
              call dynamic_matdef('tmpEEO1','q',ncf,ncf)
              itmpEEO1=mataddr('tmpEEO1')
              call dynamic_matdef('tmpEEO2','q',idimen,idimen)
              itmpEEO2=mataddr('tmpEEO2')
              call reader(ndisktr,ncf,ijprim,bl(itmpEEO1))
#ifdef SYMMETRY
              call symm_transform_ao(ncf,ns,trans1,signum,bl(itmpEEO1))
#endif
              call matsimtr('tmpEEO1','genvirt','tmpEEO2')
              call tfer(bl(itmpEEO2),xmat,idimen*idimen)
              call dynamic_matrem('tmpEEO2')
              call dynamic_matrem('tmpEEO1')
            else
              call reader(ndisktr,ncf,ijprim,xmat)
#ifdef SYMMETRY
              call symm_transform_ao(ncf,ns,trans1,signum,xmat)
#endif
c
            endif
          else if (typ.eq.'x'.or.typ.eq.'c') then
            call reader(ndisktr,idimen,ijprim,xmat)
#ifdef SYMMETRY
            call symm_transform(idimen,ns,trans1,signum,xmat)
#endif
          else    ! means that typ.eq.'tx'.or.typ.eq.'tc' .or.  typ.eq.'tt'
#ifdef SYMMETRY
          if (typ.ne.'tt') then
            call full_pair_searcher(iorb,  jorb,  iprim, jprim, ijprim,
     *                              ns, signum)
          else
            ijprim=ij
          endif
#else
            ijprim=ij
#endif
            call reader(ndisktr,idimen,ijprim,xmat)
#ifdef SYMMETRY
            if (ns.ne.0.and.typ.ne.'tt') then
              call full_symm_transform(idimen,ns,signum,xmat)
            endif
#endif
          endif
        endif
#ifndef SYMMETRY
        if (typ.eq.'x'.or.typ.eq.'e') then
          if (iorb.lt.jorb) call trans(xmat,idimen)
        endif
#endif
        return
      endif
c********************************************
      if (typ.eq.'e') then ! 'ao' & EEO need special treatment:
        if (af) then
          call fafread(ndisktr,xmat,8,ncf*ncf,1,ijprim,istatus)
          if (istatus.lt.8*ncf*ncf) then
              ERROR
              call flush(6)
          endif
        else
          call reader(ndisktr,ncf,ijprim,xmat)
        endif
#ifdef SYMMETRY
c     stop 'should not be here'
      call symm_transform_ao(ncf,ns,trans1,signum,xmat)
#endif
c
#ifndef SYMMETRY
        if (iorb.lt.jorb) call trans(xmat,ncf)
#endif
        return
      endif
c********************************************
      if (typ.eq.'tt') STOP '3ext cannot be read from bins in extrone'
c
c
      if (typ.eq.'e') STOP 'Error with EEO reading'
      call zeroit(xmat,ncf*ncf)
      if (af) then
      iafindex(2)=ij
      if (byt8) then
      do j=1,islaves
        irecij=irec(ij+irecskip*(j-1))
        do
          call bufpointer(iadd,bl(1),lbin*12+8)
          irecord=irecij
          if (irecord.eq.0) exit
          call fafread(ndisk,bl(iadd),1,lbin*12+8,1,irecord,istatus)
          if (istatus.lt.lbin*12+8) then
              ERROR
              call flush(6)
          endif
          call unpacker(3,bin8,lbin*8,indexx,lbin*4,irecij,8)
          do icount=1, lbin
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            if (mu.eq.0.or.nu.eq.0) cycle
            xmat(mu,nu)=bin8(icount)*thresh
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=bin8(icount)*thresh*list(3,iii)
            enddo
          enddo
        enddo
      enddo
      else
      iafindex(2)=ij
      do j=1,islaves
        iafindex(1)=j
        irecord=irec(ij+irecskip*(j-1))
        do i=1,irecord
          iafindex(3)=i
          iafindex(4)=0
          call fafread(ndisk,ibin4,4,lbin,4,iafindex,istatus)
          if (istatus.lt.4*lbin) then
              ERROR
              call flush(6)
          endif
          iafindex(4)=1
          call fafread(ndisk,ibin1,1,lbin,4,iafindex,istatus)
          if (istatus.lt.lbin) then
              ERROR
              call flush(6)
          endif
          iafindex(4)=2
          call fafread(ndisk,indexx,2,2*lbin,4,iafindex,istatus)
          if (istatus.lt.4*lbin) then
              ERROR
              call flush(6)
          endif
          do icount=1, lbin
            if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) 
     +          cycle
            if (ibin1(icount).eq.0) then 
              xint=ibin4(icount)
            else if (ibin1(icount).gt.0) then
              xint=ibin4(icount)
              xint=xint+SIGN(ibin1(icount)*dblmax,xint)
            else if (ibin1(icount).lt.0) then
              xint=ibin4(icount)
              xint=xint*10.0d0**(-ibin1(icount))
            endif
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=xint*thresh
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=xint*thresh*list(3,iii)
            enddo
          enddo
        enddo
      enddo
      endif
      else   ! if (af) then
      irecord=irec(ij)
      if (byt8) then
      do
        read(ndisk,rec=irecord) bin8,indexx,jrecord
        do icount=1, lbin
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            if (mu.eq.0.or.nu.eq.0) cycle
            xmat(mu,nu)=bin8(icount)*thresh
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=bin8(icount)*thresh*list(3,iii)
            enddo
        enddo
        if (jrecord.eq.0) then
          exit
        else
         irecord=jrecord
        endif
      enddo
      else
      do
        read(ndisk,rec=irecord) ibin4,ibin1,indexx,jrecord
        do icount=1, lbin
          if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) 
     +        cycle
          if (ibin1(icount).eq.0) then 
            xint=ibin4(icount)
          else if (ibin1(icount).gt.0) then
            xint=ibin4(icount)
            xint=xint+SIGN(ibin1(icount)*dblmax,xint)
          else if (ibin1(icount).lt.0) then
            xint=ibin4(icount)
            xint=xint*10.0d0**(-ibin1(icount))
          endif
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=xint*thresh
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=xint*thresh*list(3,iii)
            enddo
        enddo
        if (jrecord.eq.0) then
          exit
        else
         irecord=jrecord
        endif
      enddo
      endif
      endif ! if (af) then
c
      if (iorb.lt.jorb.and.(typ.eq.'x'.or.typ.eq.'e')) then
        call trans(xmat,ncf)
      endif
      end
c==============================================================================
      subroutine AddLeftToRight(xleft,xright,ncf)
      implicit none
      integer ncf
c
      integer i,j
      real*8 xleft(ncf,ncf),xright(ncf,ncf)
      do j=1,ncf
        do i=1,ncf
        xright(i,j)=xright(i,j)+xleft(i,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine desymm(xsym,xdesymm,idimen)
c "Desymmetrize" symmetrix matrix xsym written in triangular form and copy it
c  into square matrix xdesymm.
      implicit real*8 (a-h,o-z)
      dimension xdesymm(idimen,idimen), xsym(idimen*(idimen+1)/2)
      ntri=idimen*(idimen+1)/2
      ij=0
      do j=1,idimen
      do i=1,j
      ij=ij+1
      xdesymm(i,j)=xsym(ij)
      xdesymm(j,i)=xsym(ij)
      enddo
      enddo
      end
c==============================================================================
c The reverse of above procedure is:
      subroutine symme(xnosym,xsym,idimen)
      implicit real*8 (a-h,o-z)
      dimension xnosym(idimen,idimen),xsym(idimen*(idimen+1)/2)
      ntri=idimen*(idimen+1)/2
      ij=0
      do j=1,idimen
      do i=1,j
      ij=ij+1
      xsym(ij)=xnosym(i,j)
      if (dabs(xnosym(i,j)-xnosym(j,i)).ge.1.0d-8)
     * call nerror(1,'symm procedure',
     * 'Trying to symmetrize non-symmetric matrix',0,0)
      enddo
      enddo
      end
c==============================================================================
      subroutine UpdateCoef(iorb,   jorb,   icoef,   iresiduals, epsi,
     *                      ncf,    nmo,    nfirst,  npairs,     shift,
     *                      dodiis, vorb,   singles, ndiskr,     ndiskc,
     *                      sqij,   xmaxele,list,ipairimages,
     *                      iorev,ivrev,mtab,docansym,nbf)
      use memory
c INPUT: 
c  iorb,iorb  - orbital numbers which identify given pair ij
c  icoef      - address in bl to matrix with T_ij coefficients (AO or MO),
c       !!!!!!  This is input and OUTPUT - it is changed in place. !!!!!
c  iresiduals - as above for residuals R_ij (AO or MO)
c  epsi       - diagonal matrix with HF eigenvalues
c  ncf        - number of contracted basis functions
c  nmo        - number of occupied HF orbitals
c  nfirst     - first correlating occupied orbital
c  npairs     - number of unique ij pairs of correlating orbitals
c  shift      - level shift used for new coefficient estimation
c  dodiis     - character string determining whether to perform DIIS,
c               must be equal to 'diis' if you want diis to have performed
c vorb        - logical variable determining whether program works in MO
c               or AO mode: .true. = MO mode. Matrices have different
c               dimensions depending on vorb value!
c singles     - argument passed to diis subroutine
c ndiskr      -    as above
c ndiskc      -    as above
c OUTPUT:
c  sqij       - sum of squares of MO form of residual matrix
c  xmaxele    - max residual element (MO form)
      implicit none
      integer iorb,jorb,icoef,iresiduals,ncf,nmo,nfirst,npairs
      integer ndiskr,ndiskc,list(*),ipairimages(*),nbf
      integer iorev(*),ivrev(*),mtab(8,8)
      real*8 epsi(ncf),sqij,shift,xmaxele
      character*(*) dodiis
      logical vorb,singles,docansym
c
      real*8 enerocc,elemc,elemr,enervirt,residuum,zero
      integer icoeffMO,mataddr,iresidMO,icoeffxx
      integer i,j,nval,idimen,ii,jj,iadd,jshift,iabrep,ijrep,is,js
      data zero /0.0d0/
c
c     if (iorb.ge.jorb) then 
c       ij=iorb*(iorb-1)/2+jorb
c     else
c       ij=jorb*(jorb-1)/2+iorb
c     endif
c
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      nval=nmo-nfirst+1
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('coeffMO','q',idimen,idimen)
      icoeffMO=mataddr('coeffMO')
      call dynamic_matdef('residMO','q',idimen,idimen)
      iresidMO=mataddr('residMO')
      call matconn('coeff','q',idimen,idimen,icoef)
      call matconn('resid','q',idimen,idimen,iresiduals)
      if (docansym) then
        ijrep=mtab(iorev(iorb),iorev(jorb))
      endif
c Convert coeff to MO:
      if (vorb) then
        call matcopy('coeff','coeffMO')
        call matcopy('resid','residMO')
      else
        call matsimtr('coeff','invcano','coeffMO')
        call matsimtr('resid','cano','residMO')
      endif
      enerocc=epsi(iorb-1+nfirst)+epsi(jorb-1+nfirst)
      if (vorb) then
      call squaresum(bl(iresidMO),0,idimen,sqij)
      else
      call squaresum(bl(iresidMO),nmo,idimen,sqij)
      endif
c
      if (vorb) then
        iadd=nmo
      else
        iadd=0
      endif
c
      do j=1,idimen
        jj=iadd+j
        jshift=(j-1)*idimen-1
        js=jj-nmo
        do i=1,idimen
          ii=iadd+i
          if (docansym) then
            is=ii-nmo
            iabrep=mtab(ivrev(js),ivrev(is))
            if (iabrep.ne.ijrep) then
              bl(icoeffMO+i+jshift)=0d0
              cycle
            endif
          endif
          if ((i.le.nmo .or. j.le.nmo).and.(.not.vorb)) then
            cycle
          endif
          elemc=bl(icoeffMO+i+jshift)
          elemr=bl(iresidMO+i+jshift)
          enervirt=epsi(ii)+epsi(jj)
          residuum=elemr/(enerocc-enervirt-shift)
          elemc=elemc+residuum
          bl(icoeffMO+i+jshift)=elemc
          if (dabs(residuum).gt.xmaxele) xmaxele=dabs(residuum)
        enddo
      enddo
      if (vorb) then
         call matcopy('coeffMO','coeff')
      else
        call matsimtr('coeffMO','canotran','coeff')
      endif
      if (dodiis.eq.'diis') then
      call CCDiis(iorb,    jorb,         'store',   ncf,  nmo,
     *            npairs,  bl(iresidMO), bl(icoef), nval, vorb,
     *            singles, ndiskr,       ndiskc,    'xxx',list,
     *            ipairimages,nbf)
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
c     call CoefInit(npairs,ncf, gauss_seidel, nodisk,     ccsd,
c    *              vorb,  nmo, ivalpair,     ivirpair,   nval,
c    *              nvirt)
      subroutine CoefInit(npairs, ncf, gauss_seidel,nodisk,   ccsd_inp,
     *                    vorbpar,nmo, ivalpair,    ivirpair, ifp, 
     *                    nsym,   nval,nvirt,       ij_unique,isympairs,
     *                    nbf)
      use memory
      use ccounters
c This is also, like ExtrOne, one of the oldest procedures, present here
c since first version of CID was build. It works with companion of
c CoefRead and CoefWrite, which handle all doubles amplitudes input and
c output. This subroutine prepares environment for those reading and
c writing subroutines, reserves memory, opens files, prepares common block
c which is used by CoefRead and CoefWrite. It determines how CoefRead and
c CoefWrite will behave: write to disk or keep amplitudes in RAM memory.
c Use Gaus Seidel Algorithm or not?
c 
c ALL PARAMETERS ARE INPUT:
c  npairs       - number of unique orbital pairs
c  ncf          - number of contracted basis functions
c  gauss_seidel - this flag indicates how amplitudes will be written to
c                 disk, if .true. all new amplitudes will immediately
c                 overwrite the old ones (one file is present), if .false. 
c                 there will exist two separate files, one is with old
c                 amplitudes, other is filled with new amplitudes. If
c                 nodisk=.true. this option is irrelevant.
c nodisk        - this option determines the location of storage for
c                 amplitudes, if .true. amplitudes will be kept in RAM,
c                 if false on the disk
c                 this options is directly passed via common block to 
c                 CoefRead and CoefWrite
c ccsd_inp      - ccsd flag determines behavior of CoefRead subroutine in ccsd,
c                 passed directly to them via common.
c vorbpar       - this options determines if program is working in MO
c                 mode or AO mode. It changes sizes of matrices!
c nmo           - number of all occupied HF orbitals.
c OUTPUT is passed via common /ampldisk/
c ndiska        - unit number for disk with old amplitudes
c ndiskb        - unit number for disk with new amplitudes
c ncfa          = ncf
c internal_addr - this is address in bl where all amplitudes are kept if
c                 nodisk=.true., otherwise irrelevant and not used.
c nodiskk       = nodisk
c ireturn       - the addres in memory, reserved here for returning the
c                 amplitudes. This is strange legacy: CoefRead and
c                 CoefWrite does not require external memory reservation,
c                 they return result in its own internal location.
c                 Should be changed in future because it increases memory
c                 copying.
c icallno_propagate - the number of this procedure calls. Ised for files
c                     naming
c npairs_propagate = npairs
c ccsd         - ccsd_inp
c ndiskr       - it was disk unit number with amplitudes differences, 
c                ampl from present iter. - ampl from previous iter., not
c                used
c ivirt        - number of virtual orbitals
c vorb         = vorbpair
c af           - is af used?
c ndisk_cache  - not used in single processor version
c local_rec    - array not used in single proc. version FRAGILE, should
c                be removed from common block, maybe put in separate
c                common frof parallel use only?
c
      implicit real*8 (a-h,o-z)
      character*256 filenamea,filenameb,scrfile,iternumbera,iternumberb
      character*256 filenamer
      dimension zerocoeff(ncf,ncf),zeroMO(nbf-nmo,nbf-nmo) !FRAGILE
      logical gauss_seidel,nodisk,nodiskk,ccsd,ccsd_inp,fileopen
      logical vorb,vorbpar,af
      common /ampldisk/ndiska,ndiskb,ncfa,internal_addr,nodiskk,ireturn,
     *                 icallno_propagate,npairs_propagate,ccsd,ndiskr,
     *                 ivirt,vorb,af
      common /amplsymm/ ivalpair1,ivirpair1,ifp1,nsym1,nval1,nvirt1,
     *                  ij_unique1,isympairs1
      save /ampldisk/
      af=.false.
      ivirt=nbf-nmo
      vorb=vorbpar
      ivalpair1 = ivalpair
      ivirpair1 = ivirpair
      ifp1      = ifp
      nsym1     = nsym
      nval1     = nval
      nvirt1    = nvirt
      ij_unique1= ij_unique
      isympairs1=isympairs
      if (vorb) then
        idimen=ivirt
      else
        idimen=ncf
      endif
      ic_CoefInit=ic_CoefInit+1
      icallno_propagate=ic_CoefInit
      npairs_propagate=npairs
      ccsd=ccsd_inp
c     print *,'Coupled: ',ccsd
      if (ic_CoefInit.eq.1) call dynamic_getmem(idimen*idimen,ireturn)
      if (nodisk) then
        if (ic_CoefInit.eq.1) then 
          nodiskk=nodisk
          ncfa=ncf
c         call dynamic_getmem(npairs*idimen*idimen,ii)
          call dynamic_getmem(ij_unique*idimen*idimen,ii)
          internal_addr=ii
          do i=0,ij_unique*idimen*idimen-1
          bl(internal_addr+i)=0.0d0
          enddo
        endif
        return
      endif
      internal_addr=0
      if (gauss_seidel) then
      if (ic_CoefInit.eq.1) then
      ncfa=ncf
      ndiska=60
      ndiskb=60
      irecord=idimen*idimen*8
      nodiskk=nodisk
      write(iternumbera,'(I3.3)') ic_CoefInit
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filenamea=scrfile(1:len)//'.ampl.'//iternumbera(1:3)
      len1=len+6+3
      inquire(ndiska,opened=fileopen)
      if (fileopen) STOP 'CoefRead file opened 1'
      OPEN(UNIT=ndiska,FILE=filenamea(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=irecord)
      endif
      else
c=================
      irecord=idimen*idimen*8
      write(iternumbera,'(I3.3)') ic_CoefInit
      write(iternumberb,'(I3.3)') ic_CoefInit+1
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filenamea=scrfile(1:len)//'.ampl.'//iternumbera(1:3)
      filenameb=scrfile(1:len)//'.ampl.'//iternumberb(1:3)
      filenamer=scrfile(1:len)//'.differen'
      len1=len+6+3
      if (ic_CoefInit.eq.1) then
        ncfa=ncf
        ndiska=60
        ndiskb=61
        ndiskr=62
        nodiskk=nodisk
      endif
c
      if (ic_CoefInit.gt.1) then
        close(ndiska,STATUS='DELETE')
        close(ndiskb,STATUS='KEEP')
      endif
      inquire(ndiska,opened=fileopen)
      if (fileopen) STOP 'CoefRead file opened 2'
      inquire(ndiskb,opened=fileopen)
      if (fileopen) STOP 'CoefRead file opened 3'
      OPEN(UNIT=ndiska,FILE=filenamea(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=irecord)
      OPEN(UNIT=ndiskb,FILE=filenameb(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=irecord)
c=================
      endif
      if (ic_CoefInit.eq.1) then
        do j=1,idimen
          do i=1,idimen
          if (vorb) then 
            zeroMO(i,j)=0.0d0
          else
            zerocoeff(i,j)=0.0d0
          endif
          enddo
        enddo
        do ij=1,npairs
        if (vorb) then 
        write(ndiska,REC=ij) zeroMO
        else
        write(ndiska,REC=ij) zerocoeff
        endif
        enddo
      endif
      end
c==============================================================================
      subroutine CoefRead(ckind,iorb,jorb,iaddress)
      use memory
c This subroutine reads amplitudes. Before first use it must be
c initialized by CoefInit subroutine, which gives common /ampldisk/
c setup. After each iteration the CoefInit must be issued, which updates
c common /ampldisk/ and performs some task, especially for non
c Gauss-Seidel algorithm.
c This subroutine uses RAM or disk storage, depending on CoefInit setup
c INTENT(IN)
c ckind      - string variable 
c              refer to J. Chem. Phys. Vol. 140, No. 16, page 6286:
c                 ee - new coeficients for energy calc. purposes 
c                 tb - c bar
c                 tc - c tilda
c                 tt - usual non-modified double excitation coeficients
c                 te - new usual non-modified double excitation coeficients
c iorb,jorb  - i & j occupied orbital index
c common block /ampldisk/ explained in CoefInit
c INTENT(OUT)
c iaddress   - here is an address where coeff. matrix is returned, it is
c              allocated inside coefinit!
      implicit real*8 (a-h,o-z)
      logical nodisk,ccsd,vorb,af
      character*(*) ckind
      common /ampldisk/ idiska,idiskb,ncf,iaddr,nodisk,ireturn,icallno,
     *                  npairs,ccsd,ndiskr,ivirt,vorb,af
      common /amplsymm/ ivalpair,ivirpair,ifp,nsym,nval,nvirt,
     *                  ij_unique,isympairs
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      logical trans1,signum
      iaddress=ireturn
      if (vorb) then
        idimen=ivirt
      else
        idimen=ncf
      endif
      if (iorb.ge.jorb) then 
        ij=iorb*(iorb-1)/2+jorb
      else
        ij=jorb*(jorb-1)/2+iorb
      endif
      call pair_searcher(iorb, jorb,   iprim, jprim, ijprim,
     *                   ns,   trans1, signum)
c     ijprim=ij
      ijrec=int_array(bl(isympairs),ijprim)
      if (nodisk) then
        iref=idimen*idimen*(ijrec-1)+iaddr
        do i=0,idimen*idimen-1
          bl(iaddress+i)=bl(iref+i)
        enddo
      else
        if (ckind.eq.'tb'.or.ckind.eq.'tc'.or.ckind.eq.'tt') then
          idis=idiska
          if (af) then
           call fafRead(idis,bl(iaddress),8,idimen*idimen,1,ijrec,info)
            if (info.lt.idimen*idimen*8) then 
              print *, 'idis,idimen,ij,info: ',idis,idimen,ijrec,info
              print *,'      ABORT!'
              ERROR
              call abort()
c             write(*,'(F10.x)') bl(iaddress)
            endif
          else
            call reader(idis,idimen,ijrec,bl(iaddress))
          endif
        else ! 'ee' or 'te'
          idis=idiskb
          if (af) then
          call fafRead(idis,bl(iaddress),8,idimen*idimen,1,ijrec,info)
          if (info.lt.idimen*idimen*8) then
              ERROR
              print *,'idis,ij,idimen,info: ',idis,ijrec,idimen,info
              call flush(6)
          endif
          else
            call reader(idis,idimen,ijrec,bl(iaddress))
          endif
        endif
      endif
c
c     if (iorb.lt.jorb) then
c       call trans(bl(iaddress),idimen)
c     endif
      call symm_transform(idimen,ns,trans1,signum,bl(iaddress))
      if (ckind.eq.'tt') goto 555
      if (ckind.eq.'te') goto 555
      call matconn('internal','q',idimen,idimen,iaddress)
      if (ckind.eq.'tb') then
        call matscal('internal',0.5d0)
      endif
c     print *,'Coupled: ',ccsd
      if (ccsd) then
        call dynamic_matdef('vecti','r',idimen,1)
        ivecti=mataddr('vecti')
        call dynamic_matdef('vectj','r',idimen,1)
        ivectj=mataddr('vectj')
        call readsingles(iorb,bl(ivecti))
        call readsingles(jorb,bl(ivectj))
        call matmmul2('vecti','vectj','internal','n','t','a')
        call dynamic_matrem('vectj')
        call dynamic_matrem('vecti')
      endif
      call matdisc('internal')
 555  continue
      end
c==============================================================================
      subroutine CoefWrite(iorb,jorb,iaddress)
      use memory
c It writes PURE doubles coefficients into file or memory location.
c see CoefRead and CoefInit explanations, also for parameters and common
      implicit real*8 (a-h,o-z)
      logical nodisk,ccsd,vorb,af
      common /ampldisk/ idiska,idiskb,ncf,iaddr,nodisk,ireturn,icallno,
     *                  npairs,ccsd,ndiskr,ivirt,vorb,af
      common /amplsymm/ ivalpair,ivirpair,ifp,nsym,nval,nvirt,
     *                  ij_unique,isympairs
      if (vorb) then
        idimen=ivirt
      else
        idimen=ncf
      endif
      if (iorb.ge.jorb) then 
        ij=iorb*(iorb-1)/2+jorb
      else
        STOP 'CoefWrite: Does it make any sense to write for i<j?'
      endif
      ijrec=int_array(bl(isympairs),ij)
      if (ijrec.eq.0) return ! do not write symmetry images
      if (nodisk) then
        iref=idimen*idimen*(ijrec-1)+iaddr
        if (iaddress.eq.iref) return
        do i=0,idimen*idimen-1
        bl(iref+i)=bl(iaddress+i)
        enddo
      else
        idis=idiskb
        if (af) then
          call fafWrite(idis,bl(iaddress),8,idimen*idimen,1,ijrec,info)
        else
          call writer(idis,idimen,ijrec,bl(iaddress))
        endif
      endif
      end
c==============================================================================
      subroutine reader(idisk,idimen,ij,xmatrix)
      implicit real*8 (a-h,o-z)
      logical nodisk,ccsd,vorb,af
      common /ampldisk/ idiska,idiskb,ncf,iaddr,nodisk,ireturn,icallno,
     *                  npairs,ccsd,ndiskr,ivirt,vorb,af
      dimension xmatrix(idimen,idimen)
      read(idisk,REC=ij) xmatrix
      end
c==============================================================================
      subroutine writer(idisk,idimen,ij,xmatrix)
      implicit real*8 (a-h,o-z)
      logical nodisk,ccsd,vorb,af
      common /ampldisk/ idiska,idiskb,ncf,iaddr,nodisk,ireturn,icallno,
     *                  npairs,ccsd,ndiskr,ivirt,vorb,af
      dimension xmatrix(idimen,idimen)
      write(idisk,REC=ij) xmatrix
      end
c==============================================================================
c not used now
      subroutine ReaderCoefLin(idisk,  ij,    xmatrix,  idimen,   ncf,
     *                         af,     typ)
c 
c This subroutine is used for reading transformed amplitudes needed for EEO.
c "Transformed amplitudes" means that doubles are added up to singles and 
c symmetrized and antisymmetrized matrices are used. Caching is also included, 
c but not used at present. Also reading matrices transformed to integers is 
c possible
c Input parameters:
c idisk     - disk unit number where amplitudes are stored
c ij        - ij pair (=record number) amplitudes are read for
c idimen    - dimension of xmatrix
c ncf       - number of contracted functions
c af        - logical variable indicating whether we run single job with local
c             I/O (.false.) or parallel job using Array Files.
c typ       - if typ='s' it means that we are reading symmetrized amplitudes,
c             typ='a' means that antisymmetrized, needed only with cache.
c Output parameters:
c xmatrix   - the read matrix is deposited here
c
c COMMONS: /SymmAntiCache/, /GlobalCCSD/
c reset_s_a - this variable if .true. means that local cache should be deleted
c             and rebuild, it is relevant only with cache=.true.
c cache     - if cache=.true. local caching of amplitudes is enabled, and 
c             reset_s_a becomes relevant
c integerT  - if integerT=.true. it means that integrals are stored in integer
c             form
      use memory
      use ccounters
      implicit none
      integer idisk,ij,idimen,ncf
      real*8 xmatrix(idimen)
      logical af
      character*1 typ
c fragile!
      integer*1 irecords(PAR_MAX),irecorda(PAR_MAX)
c
      character*256 scrfile,filename
      character*2  ch_mygid
      integer ndisk_caches,ndisk_cachea,len,len1,info
      integer ntrisymm,ntrianti,i,irecord,mygid,num
      integer intmatrix
      logical reset_s_a,cache,integerT
      common /SymmAntiCache/reset_s_a
      common /GlobalCCSD/ cache,integerT
      parameter (ndisk_caches=78,ndisk_cachea=79)
c     assume abs(composed "D" matrix element)<xlargest
      real*8 xlargest
      parameter (xlargest=20.0d0)
      save irecords,irecorda
c
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (integerT) then
       num=4
      else
       num=8
      endif
      call dynamic_getmem(idimen/2+1,intmatrix)
c Local caching, when af=.true.
      if (af.and.reset_s_a.and.cache) then
        reset_s_a=.false.
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len)
        call getival('mgid',mygid)
        write(ch_mygid,'(I2.2)') mygid
        filename=scrfile(1:len)//'.symm--'//ch_mygid
        len1=len+9
        if (ic_ReaderCoefLin.gt.1) close(ndisk_caches,STATUS='delete')
        irecord=ntrisymm*num
       OPEN(UNIT=ndisk_caches,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
        filename=scrfile(1:len)//'.anti--'//ch_mygid
        len1=len+9
        if (ic_ReaderCoefLin.gt.1) close(ndisk_cachea,STATUS='delete')
        irecord=ntrianti*num
       OPEN(UNIT=ndisk_cachea,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
       do i=1, (PAR_MAX)
         irecords(i)=0
       enddo
       do i=1, (PAR_MAX)
         irecorda(i)=0
       enddo
      endif
c
c
      if (af) then
        if ((.not.cache).or.(irecords(ij).eq.0.and.typ.eq.'s'.or.
     *      irecorda(ij).eq.0.and.typ.eq.'a')) then
          if (integerT) then
            call fafread(idisk,bl(intmatrix),4,idimen,1,ij,info)
          else
            call fafread(idisk,xmatrix,8,idimen,1,ij,info)
          endif
          if (info.lt.num*idimen) then
            print *,'idisk,idimen,ij,info: ',idisk,idimen,ij,info
            ERROR
          endif
          if (cache.and.typ.eq.'s') then 
            irecords(ij)=1
            if (integerT) then
              call general_write(ndisk_caches,ij,bl(intmatrix),idimen*4)
            else
              write(ndisk_caches,REC=ij) xmatrix
            endif
          else if (cache.and.typ.eq.'a') then 
            irecorda(ij)=1
            if (integerT) then
              call general_write(ndisk_cachea,ij,bl(intmatrix),idimen*4)
            else
              write(ndisk_cachea,REC=ij) xmatrix
            endif
          endif
        else ! if (irecords(ij).eq.0.and.typ.eq.'s'.or...
          if (typ.eq.'s') then 
            if (integerT) then
              call general_read(ndisk_caches,ij,bl(intmatrix),idimen*4)
            else
              read(ndisk_caches,REC=ij) xmatrix
            endif
          else if (typ.eq.'a') then 
            if (integerT) then
              call general_read(ndisk_cachea,ij,bl(intmatrix),idimen*4)
            else
              read(ndisk_cachea,REC=ij) xmatrix
            endif
          endif
        endif ! if (irecords(ij).eq.0.and.typ.eq.'s'.or...
      else ! if (af) then
        if (integerT) then
          call general_read(idisk,ij,bl(intmatrix),idimen*4)
        else
          read(idisk,REC=ij) xmatrix
        endif
      endif ! if (af) then
      if (integerT) then
        call i2r(bl(intmatrix),xmatrix,idimen,xlargest)
      endif
      call dynamic_retmem(1)
      end
c==============================================================================
      subroutine WriterCoefLin(idisk,ij,xmatrix,idimen,af,xmax)
c This subroutine writes transformed EEO amplitudes (sum of doubles amplitudes
c and singles amplitudes, symmetrized and antisymmetrized). It also is capable 
c of transforming these amplitudes to integer form.
c Input parameters:
c idisk     - disk unit number where amplitudes are stored
c ij        - ij pair (=record number) amplitudes are read for
c idimen    - dimension of xmatrix
c af        - logical variable indicating whether we run single job with local
c             I/O (.false.) or parallel job using Array Files.
c Output parameters:
c xmatrix   - the read matrix is deposited here
c xmax      - the maximum amplitude value
c
c COMMONS: /SymmAntiCache/, /GlobalCCSD/
c Input parameters:
c reset_s_a - this variable if .true. means that local cache should be deleted
c             and rebuild, it is relevant only with cache=.true.
c cache     - if cache=.true. local caching of amplitudes is enabled, and 
c             reset_s_a becomes relevant
c integerT  - if integerT=.true. it means that integrals are stored in integer
c             form
c
      use memory
      use kinds
      implicit none
      integer idisk,ij,idimen
      logical af
      real*8 xmatrix(idimen),xmax
c
      integer intmatrix,i,info
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
c     assume abs(composed "D" matrix element)<xlargest
      real*8 xlargest
      parameter (xlargest=20.0d0)
      call dynamic_getmem(idimen/i4size+1,intmatrix)
      if (.not.integerT) then
        do i=1,idimen
          if (dabs(xmatrix(i)).ge.xmax) xmax=dabs(xmatrix(i))
        enddo
      endif
      if (integerT) then
        do i=1,idimen
          if (dabs(xmatrix(i)).ge.xmax) xmax=dabs(xmatrix(i))
          if (dabs(xmatrix(i)).ge.xlargest) then
         print *, 'WriterCoefLin, too large amplitude:   ',
     *   xmatrix(i)
        call nerror(1,'WriterCoefLin',
     1  'Too large amplitude in WriterCoefLin',1,1)
          endif
c         bl(intmatrix+i-1)=nint(xmatrix(i)*xmult)
        enddo
        call r2i(bl(intmatrix),xmatrix,idimen,xlargest)
      endif
      if (af) then
        if (integerT) then
          call fafwrite(idisk,bl(intmatrix),4,idimen,1,ij,info)
        else
          call fafwrite(idisk,xmatrix,8,idimen,1,ij,info)
        endif
      else
        if (integerT) then
c         write(idisk,REC=ij) intmatrix
          call general_write(idisk,ij,bl(intmatrix),4*idimen)
        else
          write(idisk,REC=ij) xmatrix
        endif
      endif
      call dynamic_retmem(1)
      end
c==============================================================================
      subroutine PairEnergy(icoeffAO, iexchAO, iresiAO, epairr, epair,
     *                      ncf,      equal,   nmo,     vorb,   nbf)
c This subroutine calculates pair energy for given ij pair. It outputs 
c both linear and quadratic energy value
c Input parameters:
c icoeffAO  - bl index to coefficients. In spite of the name they might be in 
c             AO or MO form depending on a program option
c iexchAO   - bl index to exchange integrals matrix for this ij. It also might
c             be in AO or MO form.
c iresiAO   - residuum matrix (AO or MO)
c ncf       - number of contracted functions
c equal     - logical variable equal=.true. if i=j for given ij
c nmo       - number of molecular orbitals
c vorb      - logical variable, if vorb=.true. all matrices are in MO form
c             otherwise in AO form
c Output parameters:
c epairr    - quadratic pair energy
c epair     - linear pair energy
c
      use memory
      implicit none
c
      integer icoeffAO,iexchAO,iresiAO,ncf,nmo,nbf
      real*8 epairr,epair
      logical equal,vorb
c
      integer idimen
c
      call mmark
      call dynamic_mmark
      call matmark
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
#if 0
      call matconn('coeAO','q',idimen,idimen,icoeffAO)
      call matconn('exchAO','q',idimen,idimen,iexchAO)
      call dynamic_matdef('exchwrk','q',idimen,idimen)
      call matcopy('exchAO','exchwrk')
      call matconn('resiAO','q',idimen,idimen,iresiAO)
      call dynamic_matdef('work11','q',idimen,idimen)
      iwork11=mataddr('work11')
c
      call matadd1('exchwrk',4.0d0,'work11')
      call matpose('exchwrk')
      call matadd1('exchwrk',-2.0d0,'work11')
      call matprodtr('coeAO','work11',epair)
#endif
      call energy_contr(bl(iexchAO),bl(icoeffAO),idimen,epair)
      if (.not.equal) epair=epair*2.0d0
c
#if 0
      call matzero('work11')
      call matpose('exchwrk')
      call matadd1('resiAO',+1.0d0,'exchwrk')
      call matadd1('exchwrk',4.0d0,'work11')
      call matpose('exchwrk')
      call matadd1('exchwrk',-2.0d0,'work11')
      call matprodtr('coeAO','work11',epairr)
#endif
      call energy_contr_r(bl(iexchAO),bl(icoeffAO),bl(iresiAO),
     *                    idimen,epairr)
      if (.not.equal) epairr=epairr*2.0d0
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine PairEnergyMP2(icoeffAO,iexchAO,iresiAO,epairr,epair,
     *                         ncf,     equal,  nmo,    vorb,  sepairr,
     *                         sepair,  nbf)
c This subroutine calculates pair energy for given ij pair. It outputs 
c both linear and quadratic energy value
c Input parameters:
c icoeffAO  - bl index to coefficients. In spite of the name they might be in 
c             AO or MO form depending on a program option
c iexchAO   - bl index to exchange integrals matrix for this ij. It also might
c             be in AO or MO form.
c iresiAO   - residuum matrix (AO or MO)
c ncf       - number of contracted functions
c equal     - logical variable equal=.true. if i=j for given ij
c nmo       - number of molecular orbitals
c vorb      - logical variable, if vorb=.true. all matrices are in MO form
c             otherwise in AO form
c Output parameters:
c epairr    - quadratic pair energy
c epair     - linear pair energy
c
      use memory
      implicit real*8 (a-h,o-z)
      data zero /0.0d0/
      logical equal,vorb
      call mmark
      call dynamic_mmark
      call matmark
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call matconn('coeAO','q',idimen,idimen,icoeffAO)
      call matconn('exchAO','q',idimen,idimen,iexchAO)
      call dynamic_matdef('exchwrk','q',idimen,idimen)
      call matcopy('exchAO','exchwrk')
      call matconn('resiAO','q',idimen,idimen,iresiAO)
      call dynamic_matdef('work11','q',idimen,idimen)
      call dynamic_matdef('work||','q',idimen,idimen)
      call dynamic_matdef('work|-','q',idimen,idimen)
      call matzero('work11')
      call matzero('work||')
      call matzero('work|-')
c
      call matadd1('exchwrk',4.0d0,'work11')
      call matadd1('exchwrk',2.0d0,'work|-')
      call matadd1('exchwrk',2.0d0,'work||')
      call matpose('exchwrk')
      call matadd1('exchwrk',-2.0d0,'work11')
      call matadd1('exchwrk',-2.0d0,'work||')
      call matprodtr('coeAO','work11',epair)
      call matprodtr('coeAO','work||',scs_par)
      call matprodtr('coeAO','work|-',scs_ant)
      sepair=6d0/5d0*scs_ant+1d0/3d0*scs_par
      if (equal) epair=epair/2.0d0
      if (equal) sepair=sepair/2.0d0
c
      call matzero('work11')
      call matzero('work||')
      call matzero('work|-')
      call matpose('exchwrk')
      call matadd1('resiAO',+1.0d0,'exchwrk')
      call matadd1('exchwrk',4.0d0,'work11')
      call matadd1('exchwrk',2.0d0,'work|-')
      call matadd1('exchwrk',2.0d0,'work||')
      call matpose('exchwrk')
      call matadd1('exchwrk',-2.0d0,'work11')
      call matadd1('exchwrk',-2.0d0,'work||')
      call matprodtr('coeAO','work11',epairr)
      call matprodtr('coeAO','work||',scs_par)
      call matprodtr('coeAO','work|-',scs_ant)
      sepairr=6d0/5d0*scs_ant+1d0/3d0*scs_par
      if (equal) epairr=epairr/2.0d0
      if (equal) sepairr=sepairr/2.0d0
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine decompact(xmat,decomp,nrow,ncol,irow,icol,ncf)
      implicit real*8 (a-h,o-z)
      dimension irow(nrow),icol(ncol), xmat(ncf,ncf),
     *          decomp(ncf,ncf)
      data zero /0.0d0/
      do j=1,ncf
        do i=1,ncf
          decomp(i,j)=zero
        enddo
      enddo
      do j=1,ncol
        jj=icol(j)
        do i=1,nrow
          ii=irow(i)
          decomp(ii,jj)=xmat(i,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine OverlapBuilder(ioverlap)
c This subroutine simply reads overlap matrix to ioverlap bl location
      use memory
      implicit real*8 (a-h,o-z)
c
      call getival('na  ',na)
      call getival('ictr',ictr)
      call getival('ibas',ibas)
      call getival('inuc',inuc)
      call getival('ncs',ncs)
      call getival('ncf',ncf)
      call getival('nsh',nsh)
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('core','s',ncf,ncf)
c     call dynamic_matdef('vect','q',ncf,ncf)
c     call dynamic_matdef('val','d',ncf,ncf)
      call matzero('core')
      icore=mataddr('core')
c     call general_write1(99,bl(ictr),12*ncs*8)
c     call general_write1(99,bl(ibas),13*nsh*8)
c     call general_write1(99,bl(inuc),5*na*8)
c
      call inton(0,na,bl(icore),bl(ictr),0,0,bl(ibas),
     1           bl(inuc),ncs)
c     call general_write1(99,bl(icore),ncf*(ncf+1)*4)
      call desymm(bl(icore),bl(ioverlap),ncf)
#if 0
      call matdiag('core','val','vect')
      write(*,20) 
      write(*,*) '                           overlap                 ' 
      do i=1,ncf
        call matelem('val',i,i,element)
        write (*,'(F20.10)') element
      enddo
      write(*,20) 
 20   FORMAT ('================================================'//
     *        '===============================')
#endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine KijklInit(irecadrx,  npairs,  ndiskx,  lbinx,   ncf,
     *                     thresh,    ifirst,  last,    F,       byt8,
     *                     nmo,       vorb,    ndisk,   af)
c This subroutine generates integrals, fully transformed into occupied MO 
c orbital basis (ij|kl) and stores them on disk. These integrals are needed
c for alpha paramters construction.
c Input parameters:
c irecadrx  - array containing the record info for K(ij) integrals
c npairs    - number of pairs
c ndiskx    - the disk unit number with Kij integrals. It might happen that
c             it is not used at all, because integrals are stored in compact 
c             form in different file. It depends on the vorb value
c lbinx     - bin length for integrals stored in ndiskx
c ncf       - number of contracted functions
c thresh    - integrals threshold
c ifirst    - first correlating MO orbital
c last      - last orrelating MO orbital, usually HOMO
c F         - Fock matrix, actually not needed
c byt8      - if byt8=.true. double precision bins are used intead of 
c             integer*4 + integer*1
c nmo       - total number of occupied orbitlas
c vorb      - if vorb=.true. then program is working in MO mode - everything 
c             is transformed to MO orbitals. Still it is possible to use AO 
c             matrices if explicitly stated
c ndisk     - logical disk unit for (ij|kl) correlating MO integrals.
c af        - logical variable indicating whether we run single job with local
c             I/O (.false.) or parallel job using Array Files.
c
      use memory
      implicit none
      integer npairs,ndiskx,lbinx,ncf,ifirst,last,nmo,ndisk
      integer irecadrx(npairs)
      real*8 thresh,F(ncf,ncf)  ! F-ock in MO basis
      logical byt8,vorb,af
c
      integer nval,imatrix,mataddr,ikij,ij,i,j,istatus
c
      nval=last-ifirst+1
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('matrix','q',ncf,ncf)
      imatrix=mataddr('matrix')
      call dynamic_matdef('kij','q',nval,nval)
      ikij=mataddr('kij')
      call matsub('occu','cano',ifirst,last)
      ij=0
      do i=1,nval
        do j=1,i
          ij=ij+1
          call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'ao', nmo,   vorb,       bl(imatrix))
          call matsimtr('matrix','occu','kij')
          call matpose('kij')
          if (af) then
            call fafwrite(ndisk,bl(ikij),8,nval*nval,1,ij,istatus)
          else
            call writer(ndisk,nval,ij,bl(ikij))
          endif
        enddo
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      function delta(i,j)
      REAL*8 delta
      delta=0.0d0
      if (i.eq.j) delta=1.0d0
      end
c==============================================================================
      subroutine squaresum(xmatrix,istart,idimen,square)
      implicit real*8 (a-h,o-z)
      dimension xmatrix(idimen,idimen)
      square=0.0d0
      do j=istart+1,idimen
        do i=istart+1,idimen
          square=square+xmatrix(i,j)*xmatrix(i,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine tjmatprodtr(iaddres1,iaddres2,trace,idimen)
      use memory
      implicit real*8 (a-h,o-z)
      idimen2=idimen*idimen
      trace=0.0d0
      do i=0,idimen2-1
      trace=trace+bl(iaddres1+i)*bl(iaddres2+i)
      enddo
      end
c==============================================================================
      subroutine tjmatprodtr1(xmatrix1,xmatrix2,trace,idimen)
      use memory
      implicit none
      integer idimen,i,j
      real*8 xmatrix1(idimen,idimen),xmatrix2(idimen,idimen),trace
      trace=0.0d0
      do j=1,idimen
        do i=1,idimen
          trace=trace+xmatrix1(i,j)*xmatrix2(i,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine PointerCoef(ij,coefkind,iaddress,nodisk,ipairs)
      use memory
c There is some probably unnecessary complication: when this procedure 
c is called with nodisk=.true. it returns the pointer to the memory, because it
c keeps amplitudes in memory instead of disk. The reason is to avoid memory
c copying.
c Otherwise it NEEDS reservation of memory outside and it puts data into
c space indicated by variable iaddress. Variable ipair is meaningful
c only when nodisk=.false. It gives the amount of pairs which have to be
c put in iaddress
c Input parameters:
c ij       - ij pair, this is a starting ij pair because the caller expectes
c            to get "ipairs" pairs starting from ij
c coefkind - if 's' it means we need     symmetric part of amplitudes
c            if 'a' it means we need antisymmetric part of amplitudes
c iaddress - input or output bl address, see discussion above
c nodisk   - if nodisk=.true. amplitudes are kept in memory, common block
c            /SymmAnti/ stores starting pointers, otherwise amplitudes are
c            read from disk
c ipairs   - How many amplitudes matrices should I return, sharting from ij
c Output parameters:
c iaddress - input or output bl address, see discussion above
c Common block /SymmAnti/:
c Input parameters:
c istartsymm - starting memory location of all symmetrized amplitudes, 
c              assigned and filled in GenAmplSymmAnti
c istartanti - as above for antisymmetrized
c ncf        - number of contracted functions
c ndisks     - disk unit number where symm ampl. are stored
c ndiska     - as above for anti
c af        - logical variable indicating whether we run single job with local
c             I/O (.false.) or parallel job using Array Files.
c
      implicit real*8 (a-h,o-z)
      common/SymmAnti/istartsymm,istartanti,ncf,ndisks,ndiska,af,npairs
      logical nodisk,af
      character*1 coefkind
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (nodisk) then
        if (af) STOP 'nodisk not supported'
        if (coefkind.eq.'s') then
          iaddress=istartsymm+(ij-1)*ntrisymm
        else if (coefkind.eq.'a') then
          iaddress=istartanti+(ij-1)*ntrianti
        else
          STOP 'Error in PointerCoef'
        endif
      else
        if (coefkind.eq.'s') then
          do i=0,ipairs-1
            irefsymm=iaddress+i*ntrisymm
            ijrec=ij+i
            call ReaderCoefLin(ndisks,ijrec,bl(irefsymm),ntrisymm,
     *                         ncf,af,'s')
          enddo
        else if (coefkind.eq.'a') then
          do i=0,ipairs-1
            irefanti=iaddress+i*ntrianti
            ijrec=ij+i
            call ReaderCoefLin(ndiska,ijrec,bl(irefanti),ntrianti,
     *                         ncf,af,'a')
          enddo
        else
          STOP 'Error in PointerCoef'
        endif
      endif
      end
c==============================================================================
c 2006_01_20 not used, but might be needed, for nonblocking reading
      subroutine GetCoef(ij,coefkind,iaddress,nodisk,ipairs)
      use memory
      implicit none
      integer ij,iaddress,ipairs
      logical nodisk
      character*1 coefkind
c
      logical cache,integerT,af
      integer istartsymm,istartanti,ncf,ndisks,ndiska,ireqIDs,ireqIDa
      integer npairs
      common/SymmAnti/istartsymm,istartanti,ncf,ndisks,ndiska,af,npairs
      common /GlobalCCSD/ cache,integerT
      common /requestID/ ireqIDs(100),ireqIDa(100) ! FRAGILE?
      integer ntrisymm,ntrianti,ijrec,i,irefsymm,irefanti,info
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (integerT) STOP 'tmp no support for integer'
        if (coefkind.eq.'s') then
          do i=0,ipairs-1
            irefsymm=iaddress+i*ntrisymm
            ijrec=ij+i
            if (ireqIDs(i+1).le.0) STOP 'Error in getcoef'
            call fafnbread(ireqIDs(i+1),bl(irefsymm),8,ntrisymm,info)
            if (info.le.0) STOP 'Error xre'
            ireqIDs(i+1)=-1
          enddo
        else if (coefkind.eq.'a') then
          do i=0,ipairs-1
            irefanti=iaddress+i*ntrianti
            ijrec=ij+i
            if (ireqIDa(i+1).le.0) STOP 'Error in getcoef'
            call fafnbread(ireqIDa(i+1),bl(irefanti),8,ntrianti,info)
            if (info.le.0) STOP 'Error xre'
            ireqIDa(i+1)=-1
          enddo
        else
          STOP 'Error in PointerCoef'
        endif
      end
c==============================================================================
c 2006_01_20 not used, but might be needed, for nonblocking reading
      subroutine OrderCoef(ij,coefkind,iaddress,nodisk,ipairs)
      use memory
      implicit none
      integer ij,iaddress,ipairs
      logical nodisk
      character*1 coefkind
c
      logical cache,integerT,af
      integer istartsymm,istartanti,ncf,ndisks,ndiska,ireqIDs,ireqIDa
      integer npairs
      common/SymmAnti/istartsymm,istartanti,ncf,ndisks,ndiska,af,npairs
      common /GlobalCCSD/ cache,integerT
      common /requestID/ ireqIDs(100),ireqIDa(100)
      data ireqIDs/100*-1/,ireqIDa/100*-1/
      save /requestID/
      integer ntrisymm,ntrianti,ijrec,i,irefsymm,irefanti,info
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (integerT) STOP 'tmp no support for integer'
        if (coefkind.eq.'s') then
          do i=0,ipairs-1
            irefsymm=iaddress+i*ntrisymm
            ijrec=ij+i
            if (ireqIDs(i+1).gt.0) STOP 'Error'
            call fafrequest(ndisks,8,ntrisymm,1,ijrec,ireqIDs(i+1))
            if (ireqIDs(i+1).le.0) STOP 'Error'
          enddo
        else if (coefkind.eq.'a') then
          do i=0,ipairs-1
            irefanti=iaddress+i*ntrianti
            ijrec=ij+i
            if (ireqIDa(i+1).gt.0) STOP 'Error'
            call fafrequest(ndiska,8,ntrianti,1,ijrec,ireqIDa(i+1))
            if (ireqIDa(i+1).le.0) STOP 'Error'
          enddo
        else
          STOP 'Error in PointerCoef'
        endif
      end
c==============================================================================
c not used now
      subroutine GenAmplSymmAnti(npairs, ncf,  nval, nodisk,    nfirst,
     *                           nmo,    vorb, af,   mastslave, ndisks,
     *                           ndiska, scr,  list,isympairs,ij_unique)
      use memory
c This subroutine constructs amplitudes needed for EEO generation. It takes
c doubles amplitudes and combines them with singles amplitudes, then it builds
c AO amplitudes and resolves each matrix into symmetric and antisymmetric part.
c If called on slave it sets up common block only.
c
c Input parameters:
c npairs     - number of correlated pairs
c ncf        - number of contracted functions, passed to common block below
c nval       - number of correlated orbitals
c nodisk     - if nodisk=.true. transformed amplitudes are kept in memory,
c              otherwise on disk
c nfirst     - first correlated orbital (orbitals are ordered according to 
c              SCF energy)
c nmo        - number of occupied SCF orbitals
c vorb       - if vorb=.true. all matrices are in transformed form. This has 
c              the following implication for this procedure that amplitudes 
c              need to be transformed from MO to AO basis
c af         - logical variable indicating whether we run single job with loca
c              I/O (.false.) or parallel job using Array Files.
c mastslave  - this subroutine is called on master and on slave. It performs
c              different tasks on slave than on master. This character variable
c              passes information whether we are on slave or on master.
c ndisks     - the disk unit number where the symmetrized amplitudes are written
c ndiska     - the disk unit number where antisymm. amplitudes are written
c Output parameters:
c scr        - here the screening matrix for EEO integarls is contructed. scr
c              contains, for given mu,nu max of T(ij)(mu,nu) for all ij.
c              To be more specific, first ncf*(ncf+1)/2 elements contains
c              max values for symmetrized T and the rest of scr contains max
c              values of antisymmetrized T.
c Common block /SymmAnti/:
c Output parameters:
c istartsymm - starting memory location of all symmetrized amplitudes, 
c              assigned and filled here if nodisk=.true.
c istartanti - as above for antisymmetrized
c incf       - number of contracted functions
c indisks    - disk unit number where symm ampl. are stored
c indiska    - as above for anti
c aff        - logical variable indicating whether we run single job with local
c              I/O (.false.) or parallel job using Array Files.
c Common block /SymmAnti/:
c Input parameters:
c cache      - if .true. it instructs transformed matrices I/O subsystem to 
c              cache amplitudes taken from AF on local disk.
c integerT   - if .true. it instructs this subroutine to convert double 
c              precision amplitudes into integer numbers (for space saving)
c              It turned out not to work, especially for diffuse basis sets 
c              because of large amount of big amplitudes. There is still
c              possible to store few large amplitudes in real*8 form, but 
c              experience shows that for good convergence we need integrals 
c              accurate to 12 decimal digits, the same applies PROBABLY to
c              amplitudes which have only 9 digits of accuracy in integer form
c
      implicit real*8 (a-h,o-z)
      common/SymmAnti/istartsymm,istartanti,incf,indisks,indiska,aff,
     *                npairsa
      integer, dimension(*) :: list
      logical nodisk,fileopen,vorb,af,aff
      character*256 scrfile,filenamea,filenames
      integer nmo,isingleMO
      character*(*) mastslave
      logical cache,integerT
      real*8 scr(ncf,ncf)
      common /GlobalCCSD/ cache,integerT
      integer isympairs(*),ij_unique
      integer :: icountere = 0
      aff=af
      xmax=0.0d0
      icountere=icountere+1
c
      STOP 'Switched off'
      npairsa=npairs
      if (vorb) then
        idimen=ncf-nmo
      else
        idimen=ncf
      endif
c
      fileopen=.false.
      incf=ncf
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      if (mastslave.eq.'slave') then
        indisks=ndisks
        indiska=ndiska
        return
      endif
      if (nodisk) then
        if (af) STOP 'nodisk not supported in symmanti'
        call dynamic_getmem(npairs*ntrisymm,istartsymm)
        call dynamic_getmem(npairs*ntrianti,istartanti)
        if (vorb) then 
          call dynamic_matdef('worktran','r',ncf,idimen)
          call dynamic_matdef('singleMO','r',idimen,1)
          isingleMO=mataddr('singleMO')
        endif
        call dynamic_matdef('singlvec','r',ncf,1)
        isinglvec=mataddr('singlvec')
        call dynamic_matdef('workhere','q',ncf,ncf)
        iworkhere=mataddr('workhere')
        ij=0
        ijprim=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            if (isympairs(ij).eq.0) cycle
            ijprim=ijprim+1
            irefsymm=istartsymm+(ijprim-1)*ntrisymm
            irefanti=istartanti+(ijprim-1)*ntrianti
            call CoefRead('tc',i,j,icoeff)
            if (vorb) then
            call matconn('coeff','q',idimen,idimen,icoeff)
            call matmmult('genvirt','coeff','worktran')
            call matmmul2('worktran','genvirt','workhere','n','t','n')
            call matdisc('coeff')
            else
            call matconn('coeff','q',ncf,ncf,icoeff)
            call matcopy('coeff','workhere')
            call matdisc('coeff')
            endif
c Singles part of EEO:
#if 1
            if (vorb) then
            call readsingles(j,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(j,bl(isinglvec))
            endif
            call matsub('singlocc','occu',i,i)
            call matmmul2('singlocc','singlvec','workhere','n','t','a')
            call matdisc('singlocc')
c
            if (vorb) then
            call readsingles(i,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(i,bl(isinglvec))
            endif
            call matsub('singlocc','occu',j,j)
            call matmmul2('singlvec','singlocc','workhere','n','t','a')
            call matdisc('singlocc')
#endif
c Singles part of EEO end
            call SymAntiBuild(bl(iworkhere),bl(irefsymm),
     *                        bl(irefanti),ncf,scr)
          enddo
        enddo
      else   !!! * * * * * * * * * DISK * * * * * * * * *
        if (vorb) then 
          call dynamic_matdef('worktran','r',ncf,idimen)
          call dynamic_matdef('singleMO','r',idimen,1)
          isingleMO=mataddr('singleMO')
        endif
        call dynamic_matdef('singlvec','r',ncf,1)
        isinglvec=mataddr('singlvec')
        call dynamic_matdef('workhere','q',ncf,ncf)
        iworkhere=mataddr('workhere')
        inquire(78,opened=fileopen)
        if (fileopen) STOP 'File already opened!'
        inquire(79,opened=fileopen)
        if (fileopen) STOP 'File already opened!'
        ndisks=78
        ndiska=79
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                 ! 80 is starting lenght
        filenames=scrfile(1:len)//'.ampl.'//'symm.'
        filenamea=scrfile(1:len)//'.ampl.'//'anti.'
        len1=len+6+5
        if (integerT) then
         iword=4
        else
         iword=8
        endif
        length=ntrisymm*iword
        if (af) then
          call fafCreatem(filenames(1:len1)//char(0),ndisks)
          call fafCreatem(filenamea(1:len1)//char(0),ndiska)
        else
        OPEN(UNIT=ndisks,FILE=filenames(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=length)
        length=ntrianti*iword
        OPEN(UNIT=ndiska,FILE=filenamea(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=length)
        endif
c
        call dynamic_getmem(ntrisymm,irefsymm)
        call dynamic_getmem(ntrianti,irefanti)
        ij=0
        ij_rec=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            if ((list(ij).eq.1.or.icountere.eq.1)
     *                        .and.isympairs(ij).ne.0) then
            call CoefRead('tc',i,j,icoeff)
            else
            cycle
            call CoefRead('tt',i,j,icoeff)
            endif
            ij_rec=ij_rec+1
            if (vorb) then
            call matconn('coeff','q',idimen,idimen,icoeff)
            call matmmult('genvirt','coeff','worktran')
            call matmmul2('worktran','genvirt','workhere','n','t','n')
            call matdisc('coeff')
            else
            call matconn('coeff','q',ncf,ncf,icoeff)
            call matcopy('coeff','workhere')
            call matdisc('coeff')
            endif
c Singles part of EEO:
            if (list(ij).eq.1) then
            if (vorb) then
            call readsingles(j,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(j,bl(isinglvec))
            endif
            call matsub('singlocc','occu',i,i)
            call matmmul2('singlocc','singlvec','workhere','n','t','a')
            call matdisc('singlocc')
c
            if (vorb) then
            call readsingles(i,bl(isingleMO))
            call matmmult('genvirt','singleMO','singlvec')
            else
            call readsingles(i,bl(isinglvec))
            endif
            call matsub('singlocc','occu',j,j)
            call matmmul2('singlvec','singlocc','workhere','n','t','a')
            call matdisc('singlocc')
            endif
c Singles part of EEO end
            call SymAntiBuild(bl(iworkhere),bl(irefsymm),
     *                        bl(irefanti),ncf,scr)
            call WriterCoefLin(ndisks,ij_rec,bl(irefsymm),ntrisymm,af,
     *                         xmax)
            call WriterCoefLin(ndiska,ij_rec,bl(irefanti),ntrianti,af,
     *                         xmax)
          enddo
        enddo
        call dynamic_retmem(2)
      endif
      call dynamic_matrem('workhere')
      call dynamic_matrem('singlvec')
      if (vorb) call dynamic_matrem('singleMO')
      if (vorb) call dynamic_matrem('worktran')
      indisks=ndisks
      indiska=ndiska
      end
c==============================================================================
c not used, but kept.
      subroutine EEOGenerator(ncf,   ncs,     nval,    ics,    kcs,
     1                        icf1,  icf2,    kcf1,    kcf2,   ndisk,
     2                        inx,   lmp2int, CMO,     iprnt,  thresh,
     3                        xmat,  halftra, nrow,   ncol,   lbin,   
     4                        npairs,icounter,irecord,indxbin,ibin4,
     5                        ibin1, bin8,    afflag, reset,  islvid,
     6                        ibasesymm,ibaseanti,max_merg_shell_size,
     7                                                      nodisk,byt8,
     8                        small, ndiskampl,max_amplitudes_size, 
     9                                          ishellsaddr,iamplitaddr,
     *                        iresultaddr,igranules,igranulesize,
     1                                            igran2pair,ipair2gran,
     2                        split,screen,mapf2s,ik_capacity,tot_xmax,
     3                        list,nstrong, isympairs,ij_unique)
      use memory
      implicit real*8 (a-h,o-z)
c Generate and transform EEO integrals 
c 
c Input parameters:
c ncf        - number of contracted basis functions
c ncs        - number of contracted shells
c nval       - number of orbitals to be transformed (usually the valence ones)
c ics,kcs:   - contracted shells for which the transformation is carried out
c icf1,icf2;kcf1,kcf2
c            - the first and the last contracted function for ics and kcs
c              shell, respectively; these parameters might describe the
c              border of subshell if shell splitting is in action
c ndisk      - the disk file unit number, where transformed int are written
c 
c inx        - array containing contraction info
c lmp2int    - bl index to AO integrals ordered as
c              (l=1:ncf,j=1:ncf,k=kcf1:kcf2,i=icf1:icf2)   where icf1 & icf2
c              are the first and the last contr. funct. in the ICS shell ,
c              similarly for kcs . NOTE THAT THEY ARE SCALED by 1/thresh
c CMO        - correlated occupied MO coefficients
c iprnt      - print level
c thresh     - integral threshold
c xmat       - ncf x ncf array used for temporary storage
c halftra    - nval x nval array used for temporary storage
c nrow       - number of non-zero rows of the AO integral matrices
c              actually not used here
c ncol       - ditto for columns
c              actually not used here
c lbin       - bin length
c npairs     - number of pairs
c icounter   - the array which ij element indicate the occupation of given 
c              bin(ij)
c irecord    - the array showing the disk record info for bin, different
c              information is kept for single mode and parallel because
c              of mulidimensioning which is possible in AF
c indxbin    - the array containing mu lam indexes of given integral
c ibin4      - this integer*4 bin array containing integrals is used if 
c              byt8=.false., otherwise not used and no memory is reserved here
c ibin1      - if byt8=.false. it contains ibin4 overflow
c bin8       - the real*8 array containing integrals is used if byt8=.true., 
c              otherwise not used and no memory is reserved here
c afflag     - the flag which informs whether AF is used or simple FORTAN IO
c              it informs also whether this is single or parallel mode.
c reset      - this flag if true informs EEOGen to transform any
c              leftover integrals and empty the integral buffer
c islvid     - ID of this slave, not used in single mode (if afflag=.false.)
c ibasesymm  - the base bl address of storage for symmetrized integrals,
c              stored here
c ibaseanti  - the base bl address of storage for antisymmetrized integrals
c max_merg_shell_size 
c            - the size of ibasesymm & ibaseanti storages in matrix
c              units (different for ibasesymm & ibaseanti)
c nodisk     - this variable indicate whether symmetrized and
c              antisymmetrized amplitudes are already kept in memory
c byt8       - if .true. the real*8 bin8 is used, otherwise bin4 and
c              bin1
c small      - this parameter if set to .true. uses symmetrized &
c              antisymmetrized matrices, otherwise full matrices are
c              multiplied, which is 2 times slower. It might be useful 
c              in case of big molecules because of better sparsity, but I 
c              don't think so (TJ). 2005-01-20: may not work
c ndiskampl  - if small=.false. it is a disk unit where full amplitude
c              matrices are kept.
c max_amplitudes_size 
c            - the same as max_merge_shell_size for "small" mode
c ishellsaddr- the pointer to storage area for integrals if small=.false.
c iamplitaddr- the pointer to storage area for amplitudes if small=.false.
c iresultaddr- the sorage for result in "big" case (small=.false.)
c igranules  - number of big chunks if ij pairs. This parameter is used
c              by sorting routine, all ij pairs bins which belong to a given
c              granule are written at once.
c igranulesize - size of chunk of ij's
c igran2pair - the two-dimensional array which maps granule to its starting
c              and ending pair addresses, e.g. granule starts at ij pair
c              no 4 and ends at ij pair 7. In this case igranulesize
c              will be 7-4+1=4
c ipair2gran - the one-dimensional array which indicates which
c              granule given pair belongs to
c split      - this variable if set to .true. allows triggering of calculation
c              of EEOs even if not all integrals from given shell are
c              (anti)symmetrized. After EEOs calculation the integral storage 
c              is flushed and collecting of integrals is continued
c screen     - this ncf**2 one-dimensional array contains max over all ij
c              of all amplitudes for given mulam composite index, used
c              for screening, see GenAmplSymmAnti for more details
c mapf2s     - the array which, for given contraced function (index) 
c              returns the shell number this function belongs to.
c ik_capacity - the number which shows how much memory for
c               (anti)symmetrized integrals is still available, in units 
c               which stands for amount of matrices which fit this free area.
c               It is used mainly when split=.false. for memory
c               assessment outside this subroutine
c tot_xmax   -  the largest MO residuum, used for variable screening threshold
c               assessment.
c               The variable threshold feature is switched off at present.
c Commons:
c /timingstj/ contains timings variables
c /indexsaver/:
c ik         - the counter which counts the amount of integral AO matrices
c              stored for EEO generation.
c ikindex(2,*)- Here for given ik the contraced function i,k indices are kept
c               ikindex(1,ik) = icf, ikindex(2,ik) = kcf
c               used later when results are passed to sort subroutine
c /zero_in_EEO/:
c zeroed      - counts the number of integrals omitted by screening
c zeroed_max  - the number of integrals which could be omitted if screening
c               would be applied to each icf,kcf matrix separately, not
c               in shell blocks
      dimension xmat(ncf,ncf),halftra(nval,nval),irow(ncf),
     *          icol(ncf),screen(*)
      dimension CMO(ncf,*),irow1(ncf),icol1(ncf),icounter(*)
      dimension inx(12,*),mapf2s(*)
c ALL ELAPS TIMES:
c                        mult   reads   sort zeros matrix construction
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /indexsaver/ik,ikindex(2,50) ! really is more, up to 10 000
      integer*2 indxbin(2,npairs,lbin)
      integer*4 ibin4(npairs,lbin)
      integer*1 ibin1(npairs,lbin)
      integer*4 igran2pair(2,*),ipair2gran(*)
      integer irecord(*),list(*),nstrong
      real*8    bin8(npairs,lbin),zeroed
      logical nodisk,byt8,small,afflag,reset,split,str_reset
      common /zero_in_EEO/ zeroed,zeroed_max
      integer   ishell_map(2,5),ishell_address_s(3,(MAX_MERGE)), ! FRAGILE
     *                          ishell_address_a(3,(MAX_MERGE))  ! FRAGILE
      integer*2 ishell_index(2,(MAX_MERGE))                      ! FRAGILE
      integer   igroup_no,imax_stripsize,iplus,iminus,icallno
      integer isympairs(*),ij_unique
      save   ishell_map,ishell_address_s,ishell_address_a,ishell_index
      save   igroup_no,imax_stripsize,iplus,iminus
      data ishell_map/10*0/
      data ics_prev/0/, kcs_prev/0/, iwhich_shell/0/,icallno/0/
      data ireal_shell_size_left/0/
      save ireal_shell_size_left
      save ics_prev,    kcs_prev,    iwhich_shell,   icallno
      real*8 tot_xmax
c  return if there are no integrals
      icallno=icallno+1
      if(nrow.eq.0.or.ncol.eq.0) RETURN
      if (ik.eq.0) then
        ics_prev=0
        kcs_prev=0
        iplus =ibasesymm
        iminus=ibaseanti
        if (icallno.gt.1) then
        if (.not.nodisk.and.iwhich_shell.ne.0.or.
     *           nodisk.and.iwhich_shell.ne.1) then 
          print *, 'Error 5yfx23',iwhich_shell
          STOP 'Error 5yfx23'
        endif
        endif
      endif
c
      ncf2=ncf*ncf
c
      call mmark
      call dynamic_mmark
      call matmark
c
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ncf*(ncf-1)/2
c  Determine the product of the shell sizes (primitive method):
      isubshell_size=0
      do icf=icf1,icf2  ! lower triangle!
        do kcf=kcf1,kcf2
          if(kcf.gt.icf) exit
          isubshell_size=isubshell_size+1
        enddo
      enddo
      if (reset) then ! Flush all ik accumulated
      if (small) then
        if (ik.ge.1) then
            print *, 'The last FLUSHED ik will be: ',ik
ctmp        call dynamic_lock(bl(ibasesymm),isize)
ctmp        call dynamic_lock(bl(ibaseanti),isize)
c SCREENING SUPPLEMENT:
            ishell_map(2,iwhich_shell)=ik ! stop
            str_reset=.true.
            ics_prev=0
            kcs_prev=0
            ik_capacity=max_merg_shell_size
            if (.not.nodisk) then
            call strip_builder(str_reset,iwhich_shell,ishell_map,
     *                         ibasesymm,ibaseanti,ncf,
     *                         ishell_address_s,ishell_address_a,
     *                         ishell_index,igroup_no,imax_stripsize,
     *                         screen,iplus,iminus,tot_xmax)
            endif
            iplus =ibasesymm
            iminus=ibaseanti
c SCREENING END.
            call MultiplySort(ntrisymm,ntrianti,ik,ibasesymm,ibaseanti,
     *                        npairs,ikindex, halftra,nval,lbin,
     *                        ndisk,thresh,icounter,irecord,indxbin,
     *                        ibin4,ibin1,bin8,afflag,islvid,
     *                        nodisk,byt8,igranules,igranulesize,       !->
     *                                                       igran2pair,
     *                        ipair2gran,ishell_address_s,              !->
     *                          ishell_address_a,ishell_index,igroup_no,
     *                        imax_stripsize,list,nstrong,isympairs,
     *                        ij_unique)
ctmp        call dynamic_unlock(bl(ibasesymm),isize)
ctmp        call dynamic_unlock(bl(ibaseanti),isize)
            print *, 'The last FLUSHED ik was: ',ik
            call flush(6)
        endif
            ik=0
            call matremark
            call retmark
            call dynamic_retmark
            return
      else ! if (small) then
        if (ik.ge.1) then
            call BigMuliplyEEO(nval,ncf,npairs,lbin,ndisk,halftra,
     1                         thresh,icounter,irecord,indxbin,ibin4,
     2                         ibin1,bin8,afflag,reset,islvid,byt8,
     3                         lmp2int,ndiskampl,ik,
     4                         max_amplitudes_size,ishellsaddr,
     5                         iamplitaddr,iresultaddr,ikindex) 
            print *, 'The last FLUSHED ik was: ',ik
            call flush(6)
        endif
            call matremark
            call retmark
            call dynamic_retmark
            ik=0
            return
      endif! if (small) then
      endif ! reset
c
      if (small) then !shell merging with symm-anti trick
c Determine current shell size (shell, not this part size)
c  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *  *
c SCREENING:
      if ((ics.ne.ics_prev.or.kcs.ne.kcs_prev)) then
        if (ics_prev.ne.0.and.kcs_prev.ne.0.and.
     *       ireal_shell_size_left.ne.0) then
          STOP 'Error y4823h48fh3'
        endif
        if (ireal_shell_size_left.eq.0) then
        ishell_start=inx(11,ics)+1
        ishell_stop =inx(10,ics)
        kshell_start=inx(11,kcs)+1
        kshell_stop =inx(10,kcs)
        do i=ishell_start,ishell_stop
          if (i.lt.icf1) then
          print *, '58ehfs: I do not expect this condition to be true'
          print *, i,icf1
          call flush(6)
          endif
          do k=kshell_start,kshell_stop
            if (k.lt.kcf1) then
            print *, '638dhw: I do not expect this condition to be true'
            print *, k,kcf1
            call flush(6)
            endif
            if (k.gt.i) exit
            ireal_shell_size_left=ireal_shell_size_left+1
          enddo
        enddo
        endif
        ics_prev=ics
        kcs_prev=kcs
        if (.not.nodisk) then 
          iwhich_shell=iwhich_shell+1
        else 
          iwhich_shell=1
        endif
        ishell_map(1,iwhich_shell)=ik+1                  ! start
c Check capacity:
        imemlefts=ibasesymm+(max_merg_shell_size)*ntrisymm-iplus
        imemlefta=ibaseanti+(max_merg_shell_size)*ntrianti-iminus
        ik_scr_capacity=min(imemlefts/ntrisymm,imemlefta/ntrianti)
        if (ik_scr_capacity.lt.1) then
          STOP 'Error, too little space for KExtrConstr 1'
        endif
        if (ik_scr_capacity.ge.ireal_shell_size_left) then
          ishell_map(2,iwhich_shell)=ik+ireal_shell_size_left ! stop
        else                         ! shell size does not fit available memory
          ishell_map(2,iwhich_shell)=ik+ik_scr_capacity  ! tmp for
c                                                          cutting assessment
        endif
        if (iwhich_shell.eq.2.and..not.nodisk) then
          str_reset=.false.
          call strip_builder(str_reset,iwhich_shell,ishell_map,
     *                       ibasesymm,ibaseanti,ncf,
     *                     ishell_address_s,ishell_address_a,
     *                     ishell_index,igroup_no,imax_stripsize,
     *                     screen,iplus,iminus,tot_xmax)
        endif
c Check capacity again, because iplus,iminus are updated
        imemlefts=ibasesymm+(max_merg_shell_size)*ntrisymm-iplus
        imemlefta=ibaseanti+(max_merg_shell_size)*ntrianti-iminus
        ik_scr_capacity=min(imemlefts/ntrisymm,imemlefta/ntrianti)
        if (ik_scr_capacity.lt.1) then
          STOP 'Error, too little space for KExtrConstr 2'
        endif
        if (ik_scr_capacity.ge.ireal_shell_size_left) then
          ishell_map(2,iwhich_shell)=ik+ireal_shell_size_left
        else
          ishell_map(2,iwhich_shell)=ik+ik_scr_capacity
        endif
      endif
c SCR END
      do icf=icf1,icf2  ! lower triangle!
        do kcf=kcf1,kcf2
          if(kcf.gt.icf) exit
          ik=ik+1
          ireal_shell_size_left=ireal_shell_size_left-1
          if (ireal_shell_size_left.lt.0) STOP 'Error 6h4gd'
          ikindex(1,ik)=icf
          ikindex(2,ik)=kcf
c Check capacity:
          imemlefts=ibasesymm+(max_merg_shell_size)*ntrisymm-iplus
          imemlefta=ibaseanti+(max_merg_shell_size)*ntrianti-iminus
          ik_capacity=min(imemlefts/ntrisymm,imemlefta/ntrianti)
          if (ik_capacity.lt.1) then
            STOP 'Error, too little space for KExtrConstr 3'
          endif
c
          call secund(co0)
          call elapsec(eco0)
          call KExtrConstr(ncf,  icf,  kcf,  bl(lmp2int), icf1, 
     *                     icf2, kcf1, kcf2, bl(iplus),   bl(iminus))
          iplus =iplus +ntrisymm
          iminus=iminus+ntrianti
          call secund(co1)
          call elapsec(eco1)
          const=const+co1-co0
          eonst=eonst+eco1-eco0
          if (icf1.eq.1.and.kcf1.eq.1) then 
            print *, 'The ik in the last shells was: ',ik
          endif
c         if (ik.gt.max_merg_shell_size) then
c           STOP 'Error in EEO generator'
c         endif
c  The difference between top=ibasesymm+(max_merg_shell_size)*ntri????
c  and present location of iplusminus should allow putting one matrix - 
c ik_capacity should be >= 1
c
          imemlefts=ibasesymm+(max_merg_shell_size)*ntrisymm-iplus
          imemlefta=ibaseanti+(max_merg_shell_size)*ntrianti-iminus
          ik_capacity=min(imemlefts/ntrisymm,imemlefta/ntrianti)
          if ((ik_capacity.lt.1.or.ik.ge.(MAX_MERGE)).and.split) then
C SCR:
c           print *, '***********************  CUT!! '//
c    *               '************************************'
            if (ishell_map(2,iwhich_shell).ne.ik) STOP 'Error 93bhg6'
            str_reset=.true.
            ik_capacity=max_merg_shell_size
            if (.not.nodisk) then
            call strip_builder(str_reset,iwhich_shell,ishell_map,
     *                         ibasesymm,ibaseanti,ncf,
     *                         ishell_address_s,ishell_address_a,
     *                         ishell_index,igroup_no,imax_stripsize,
     *                         screen,iplus,iminus,tot_xmax)
            endif
            iplus =ibasesymm
            iminus=ibaseanti
c if there are no ik pairs to do in this icf,kcf loop:
        if (icf.eq.icf2.and.(kcf.eq.kcf2.or.kcf.eq.icf2)) then
          ics_prev=0
          kcs_prev=0
          ishell_stop =inx(10,ics)
          kshell_stop =inx(10,kcs)
          if (icf.eq.ishell_stop.and.(kcf.eq.ishell_stop.or.
     *        kcf.eq.kshell_stop).and.ireal_shell_size_left.ne.0) then
            print *, "Error: ", ireal_shell_size_left
            call flush(6)
c           write(*,'(F20.x)') tot_xmax
            call abort()
            STOP 'Error 53h5jdfe00'
          endif
        else
c After this we have to determine starting point for ishell_map and
c predicted stop point, because we are inside ik loop.
          if (.not.nodisk) then 
            iwhich_shell=iwhich_shell+1
          else 
            iwhich_shell=1
          endif
c Set capacity:
          ik_scr_capacity=max_merg_shell_size
c Set the ishell_map array which has to be set here, because no entrance
c point will be executed
          ishell_map(1,iwhich_shell)=1
          if (ik_scr_capacity.ge.ireal_shell_size_left) then
            ishell_map(2,iwhich_shell)=ireal_shell_size_left ! ik=0
          else
            ishell_map(2,iwhich_shell)=ik_scr_capacity       ! ik=0
          endif
        endif ! if (icf.eq.icf2.and.(kcf.eq.kcf2.or.kcf.eq.icf2))
ctmp        call dynamic_lock(bl(ibasesymm),isize)
ctmp        call dynamic_lock(bl(ibaseanti),isize)
            call MultiplySort(ntrisymm,ntrianti,ik,ibasesymm,ibaseanti,
     *                        npairs,ikindex,halftra,nval,lbin,
     *                        ndisk,thresh,icounter,irecord,indxbin,
     *                        ibin4,ibin1,bin8,afflag,islvid,
     *                        nodisk,byt8,igranules,igranulesize,       !->
     *                                                       igran2pair,
     *                        ipair2gran,ishell_address_s,              !->
     *                          ishell_address_a,ishell_index,igroup_no,
     *                        imax_stripsize,list,nstrong,isympairs,
     *                        ij_unique)
ctmp        call dynamic_unlock(bl(ibasesymm),isize)
ctmp        call dynamic_unlock(bl(ibaseanti),isize)
            ik=0
          endif
        end do
      end do
      else !if (small) then - no shell merging
c The bl(lmp2int) matrix is set of TRANSPOSED K() matrices
      do icf=icf1,icf2  ! lower triangle!
        do kcf=kcf1,kcf2
          if(kcf.gt.icf) exit
          ik=ik+1
          ikindex(1,ik)=icf
          ikindex(2,ik)=kcf
          iKone =ishellsaddr+(ik-1)*ncf*ncf
          call Kextractor(ncf, icf, kcf, bl(lmp2int), nrow,
     1                    ncol, irow, icol, icf1, icf2, kcf1, kcf2,
     2                    bl(iKone))
          if (icf1.eq.1.and.kcf1.eq.1) then 
            print *, 'The ik in the last shells was: ',ik
          endif
          if (ik.gt.max_merg_shell_size) then
            STOP 'Error in EEO generator'
          endif
          if (ik.eq.max_merg_shell_size.and.split) then
            call BigMuliplyEEO(nval,ncf,npairs,lbin,ndisk,halftra,
     1                         thresh,icounter,irecord,indxbin,ibin4,
     2                         ibin1,bin8,afflag,reset,islvid,byt8,
     3                         lmp2int,ndiskampl,max_merg_shell_size,
     4                         max_amplitudes_size,ishellsaddr,
     5                         iamplitaddr,iresultaddr,ikindex) 
            ik=0
          endif
        end do
      end do
c
      endif !if (small) then
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine KExtrConstr(ncf,  icf,  kcf,  xint,  icf1,
     1                       icf2, kcf1, kcf2, xplus, xminus)
c This subroutine resolves given matrix into its symmetric and anisymmetric
c parts
c Input parameters:
c ncf       - number of contracted functions
c icf,kcf   - the indices of extracted matrix
c xint      - the integrals array taken from integral subroutine
c icf1,icf2 - parameters which determine borders of xint, for given xint
c             icf indices go from icf1 to icf2.
c kcf1,kcf2 - as above for kcf indices
c Output parameters:
c xplus     - the symmetric part of K(icf,kcf)
c xminus    - the antisymmetric part of K(icf,kcf)
      implicit none
      integer ncf,icf,kcf,icf1,icf2,kcf1,kcf2
      real*8    xint(ncf,ncf,kcf1:kcf2,icf1:icf2)
      real*8    xplus(*),xminus(*)
c
      real*8 zero,one,two,xmult,a,b,symmetric,antisymmetric
      integer ij,i,j
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
c
      ij=0
      xmult=sqrt(2.0d0)/2.0d0
      do j=2,ncf
        do i=1,j-1
          ij=ij+1
          a=xint(j,i,kcf,icf)
          b=xint(i,j,kcf,icf)
          symmetric=a+b
          antisymmetric=a-b
c         xplus(ij) = (xint(j,i,kcf,icf)+xint(i,j,kcf,icf))*xmult
c         xminus(ij)= (xint(j,i,kcf,icf)-xint(i,j,kcf,icf))*xmult
          xplus(ij) = symmetric*xmult
          xminus(ij)= antisymmetric*xmult
        enddo
      enddo
      do i=1,ncf
        ij=ij+1
        symmetric=xint(i,i,kcf,icf)
        xplus(ij)=symmetric
      enddo
      end
c==============================================================================
c==============================================================================
      subroutine Kextractor( ncf,   icf,   kcf,    xint,  nrow,
     1                       ncol,  irow,  icol,   icf1,  icf2, 
     2                       kcf1,  kcf2,  xmatrix)
c This subroutine is some old piece of legacy code. It should not be used
c in production program, because it performs superflous memory copying. Instead
c of using this subroutine I should read integrals directly from xint array.
c The only reason why this subroutine might be useful is for extraction one
c matrix from compressed xint. For EEO I do not use zero rows and columns
c of xint because of different array arrangement for EEO, xint is zeroed
c in integral subroutine and "zero" columns contain real zeros, not garbage.
c This subroutine is called by "big" model of EEO, which is not used currently.
c Input parameters:
c ncf       - number of contracted functions
c icf,kcf   - the indices of extracted matrix
c xint      - the integrals array taken from integral subroutine
c nrow      - number of nonzero rows in xint (NOT USED)
c ncol      - number of nonzero columns in xint (NOT USED)
c irow      - the array resolving compressed row number (index) into real 
c             row number 
c icol      - as above for columns
c icf1,icf2 - parameters which determine borders of xint, for given xint
c             icf indices go from icf1 to icf2.
c kcf1,kcf2 - as above for kcf indices
c Output parameters:
c xmatrix   - to this location the matrix with given i,k index is written
c
      implicit real*8 (a-h,o-z)
      dimension xint(ncf,ncf,kcf1:kcf2,icf1:icf2)
      dimension irow(ncf),icol(ncf),irow1(ncf),icol1(ncf)
      dimension xmatrix(ncf,ncf)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
c
      do j=1,ncf
        do i=1,ncf
          xmatrix(i,j)=xint(j,i,kcf,icf)
        enddo
      enddo
      end
c==============================================================================
      subroutine TraditMult(A,B,C,kk,ll,mm)
      implicit real*8 (a-h,o-z)
      dimension A(ll,kk),B(ll,mm),C(kk,mm)
      do k=1,kk
        do m=1,mm
          C(k,m)=0.0d0
          do l=1,ll
          C(k,m)=C(k,m)+A(l,k)*B(l,m)
          enddo
        enddo
      enddo
      end
c==============================================================================
      subroutine MultiplySort(ntrisymm,ntrianti,merged_size,ibasesymm, !->
     *                                                        ibaseanti,
     *                        npairs,ikindex, halftra, nval,   lbin,
     *                        ndisk, thresh,  icounter,irecord,indxbin,
     *                        ibin4, ibin1,   bin8,    afflag, islvid,
     *                        nodisk,byt8,    igranules,igranulesize,  !->
     *                                                       igran2pair,
     *                        ipair2gran,ishell_address_s,      !->
     *                          ishell_address_a,ishell_index,igroup_no,
     *                        imax_stripsize,list,nstrong,isympairs,
     *                        ij_unique)
c This subroutine performs multiplication of integrals by amplitudes, and
c passes results to sorting subroutine, which finally writes them to disk.
c MultiplySort works in two modes: different code is used when amplitudes are
c kept in memory and different when are read from disk. In the latter case 
c integrals are assumed to be kept in compressed form (screening)
c Parameters:
cI  ntrisymm       - = ncf*(ncf+1)/2, where ncf - number of contracted fun.
cI  ntrianti       - = ncf*(ncf-1)/2
cI  merged_size    - how many of integrals exchange matrices are merged when
c                    this subroutine is called.
cI  ibasesymm      - the bl index to memory location where transformed 
c                    symmetrized integrals are stored.
cI  ibaseanti      - ditto for antisymmetrized int.
cI  npairs         - number of correlated pairs.
cI  ikindex        - this array contains, for given ik, the icf and kcf indices
c                    of matrix stored as "ik" matrix. ik goes from 1 to 
c                    merged_size, and after each MultiplySort is set to 0. We 
c                    need to keep track of real exchange matrices indices.
cIO halftra        - this matrix is used temporarily for collection of results
cI  nval           - number of correlated orbitals
cI  lbin           - length of bin
cI  ndisk          - disk where sort subroutine writes transformed EEO 
c                    integrals.
cI  thresh         - integral threshold
cIO icounter       - counter array indicating bins occupation
cIO irecord        - the array where sort subroutine stores records info
cIO indxbin        - the indices of integral written to bin
cIO ibin4          - the integer*4 bin storge, used if byt8=.false.
cIO ibin1          - the overflow of integral stored in ibin4
cIO bin8           - double precision bin, used if byt8=.true.
cI  afflag         - logical variable indicating whether we run single job 
c                    with loca I/O (.false.) or parallel job using Array Files.
cI  islvid         - the slave group number used for AF indexing
cI  nodisk         - .true.: all amplitudes are kept in memory and different 
c                    is code used here
cI  byt8           - if .true. bin8 are used otherwise ibin4, ibin1 pair
cI  igranules      - number of bin granules, used for writing bunches of bins
cI  igranulesize   - the size of bin granule (bunch of bins) in no of bins units
cI  igran2pair     - the two-dimensional array which maps granule to its 
c                    starting and closing pair addresses, e.g. granule starts 
c                    at ij pair no 4 and ends at ij pair 7. In this case 
c                    igranulesize will be 7-4+1=4
cI  ipair2gran     - the one-dimensional array which indicates which
c                    granule given pair belongs to
cI  ishell_address_s-see definitions below
cI  ishell_address_a-see definitions below
cI  ishell_index   - this array for given strip (bunch of ik matrices) returns
c                    its starting and closing address.
cI  igroup_no      - number of stripes of ik
cI  imax_stripsize - max of all sizes of all stripes.
c
      use memory
      implicit real*8 (a-h,o-z)
      dimension icounter(*)
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      integer*2 indxbin(2,npairs,lbin)
      integer*4 ibin4(npairs,lbin)
      integer*4 igran2pair(2,*),ipair2gran(*)
      integer*1 ibin1(npairs,lbin)
      real*8    bin8(npairs,lbin)
      logical smal,nodisk,byt8,afflag,non_block,ish_opened
      dimension ikindex(2,merged_size), halftra(nval,nval)
      character*2 afname
      integer isympairs(*),ij_unique
      integer*2 ishell_index(2,*)
      integer   ishell_address_s(3,*)
      integer   ishell_address_a(3,*)! defined as for ishell_address_s
      integer, dimension(*) :: list
c 1 matrix addr.,                               <-----|
c 2 nrow =-1 if not compressed, 0 - zero matrix       |
c 3 compression info addr relative to element 1 of this array
c compresion info: if ishell_address(3,*)=-1, matrix is not compressed,
c 0 - zero matrix
c calculate number of symmetry unique strong pairs:
      ij=0
      nstrong_unique=0
      if (nstrong.eq.npairs) then
        do i=1,nval
          do j=1,i
            ij=ij+1
            if (isympairs(ij).eq.0) cycle
            nstrong_unique=nstrong_unique+1
          enddo
        enddo
      else
        do i=1,nval
          do j=1,i
            ij=ij+1
            if (list(ij).ne.1.or.isympairs(ij).eq.0) cycle
            nstrong_unique=nstrong_unique+1
          enddo
        enddo
      endif
c
      if (afflag) then
      non_block=.true.
      else
      non_block=.false.
      endif
      non_block=.false.
c
      call mmark
      call dynamic_mmark
      call matmark
c 
      if (nodisk) then
      call matconn('Kplus' ,'r',ntrisymm,merged_size,ibasesymm)
      call matconn('Kminus','r',ntrianti,merged_size,ibaseanti)
      call elapsec(exread0)
      call secund(xread0)
      call PointerCoef(1,'s',icoefsymm,nodisk,ipairs)
      call PointerCoef(1,'a',icoefanti,nodisk,ipairs)
      call elapsec(exread1)
      call secund(xread1)
      sconstr=sconstr+xread1-xread0
      econstr=econstr+exread1-exread0
      call matconn('Cplus' ,'r',ntrisymm,ij_unique,icoefsymm)
      call matconn('Cminus','r',ntrianti,ij_unique,icoefanti)
      call dynamic_matdef('trans1','r',merged_size,ij_unique)
      call dynamic_matdef('trans2','r',merged_size,ij_unique)
      it1=mataddr('trans1')
      it2=mataddr('trans2')
      call secund(xmult0)
      call elapsec(exmult0)
      call matmmul2('Kplus','Cplus','trans1','t','n','n')
      call matmmul2('Kminus','Cminus','trans2','t','n','n')
      call secund(xmult1)
      call elapsec(exmult1)
      strace=strace+xmult1-xmult0
      etrace=etrace+exmult1-exmult0
!  Sort the half-transformed integrals
      call elapsec(exsort0)
      call secund(xsort0)
      do ik=1,merged_size
          icf=ikindex(1,ik)
          kcf=ikindex(2,ik)
          ij=0
          ijprim=0
          do i=1,nval
            do j=1,i
              ij=ij+1
              if (isympairs(ij).eq.0) then 
                halftra(i,j)=0.0d0
                halftra(j,i)=0.0d0
                cycle
              endif
              ijprim=ijprim+1
              call matelem('trans1',ik,ijprim,elem1)
              call matelem('trans2',ik,ijprim,elem2)
              if (i.eq.j) then 
                halftra(i,j)=elem1
                cycle
              endif
              if (icf.eq.kcf) then 
                halftra(i,j)=elem1
                halftra(j,i)=elem1
                cycle
              endif
              halftra(i,j)=elem1+elem2
              halftra(j,i)=elem1-elem2
            enddo
          enddo
      if (igranulesize.le.1) then ! Have to be true for single processor!
          call BinSort_tj( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1,  bin8, afflag, afname,
     3                  islvid,byt8)
      else
          call BinSortGr( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1,  bin8, afflag, afname,
     3                  islvid,byt8,
     *          igranules,igranulesize,igran2pair,ipair2gran)
      endif
      enddo
      call elapsec(exsort1)
      call secund(xsort1)
      eort=eort+exsort1-exsort0
      sort=sort+xsort1-xsort0
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c *  WITH DISK! * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
      else
      if (.not.afflag) then
      nstrong=npairs
      endif
ctmp  call matconn('Kplus' ,'r',ntrisymm,merged_size,ibasesymm)
ctmp  call matconn('Kminus','r',ntrianti,merged_size,ibaseanti)
c
      call dynamic_getmem(merged_size*npairs,ires1)
      call dynamic_getmem(merged_size*npairs,ires2)
c FRAGILE!
      print *,'imax_stripsize: ',imax_stripsize
      call dynamic_getmem(imax_stripsize*25,ipartres1)
      call dynamic_getmem(imax_stripsize*25,ipartres2)!25 max amplitudes stripe
c
      call dynamic_show_free(memavail)
      memavail=memavail-200 000 ! FRAGILE - I take fragmentation into accout
      ipas_siz=memavail/(ntrisymm+ntrianti)-1 !Pointercoef allocates memory (-1)
      if (ipas_siz.gt.50) ipas_siz=50
c     if (ipas_siz.gt.npairs) ipas_siz=npairs
      if (ipas_siz.gt.nstrong_unique) ipas_siz=nstrong_unique
      if (ipas_siz.ne.1) ipas_siz=ipas_siz/2 ! for compression
      if (ipas_siz.gt.2*nstrong_unique) ipas_siz=nstrong_unique 
c                                          !? Why ipas_siz.gt.2*nstrong?
      print *, 'Amplitudes matrix size: ',ipas_siz, memavail,ntrisymm,
     *                                    ntrianti
      if (ipas_siz.le.4) then 
      call dynamic_dump()
      endif
      if (ipas_siz.lt.1) then 
      STOP 'EEO error in ipas_siz determination'
      endif
c
cAFil if (non_block) then
cAFil call dynamic_getmem(ipas_siz*ntrisymm,icoefsymm1)
cAFil call dynamic_getmem(ipas_siz*ntrianti,icoefanti1)
cAFil call dynamic_getmem(ipas_siz*ntrisymm,icoefsymm2)
cAFil call dynamic_getmem(ipas_siz*ntrianti,icoefanti2)
cAFil else
      call dynamic_getmem(ipas_siz*ntrisymm,icoefsymm)
      call dynamic_getmem(ipas_siz*ntrianti,icoefanti)
c
c compressed:
      call dynamic_getmem(ipas_siz*ntrisymm,icoefsymm_c)
      call dynamic_getmem(ipas_siz*ntrianti,icoefanti_c)
cAFil endif
c
c
      iremainder=  mod(nstrong_unique,ipas_siz)
      npass     =  nstrong_unique/ipas_siz
      if (iremainder.gt.0) npass=npass+1
c
c
      if (non_block) then
        icoefsymm_predict=icoefsymm
        icoefanti_predict=icoefanti
        ij_start_predict= 1
        ij_stop_predict = ipas_siz
        if (ij_stop_predict.gt.npairs) ij_stop_predict=npairs
        ipasssize_predict=ij_stop_predict-ij_start_predict+1
        iswitcher=1
        call OrderCoef(ij_start_predict,'s',icoefsymm,nodisk,
     *                 ipasssize_predict)
        call OrderCoef(ij_start_predict,'a',icoefanti,nodisk,
     *                 ipasssize_predict)
      endif
c
c   LOOP:
c
      pass: do ipaser=1,npass
      ijstart=(ipaser-1)*ipas_siz+1
      if (ijstart.gt.nstrong_unique) STOP 'Straaange errror'
      if (ipaser.ne.npass) then
      ijstop =ipaser*ipas_siz
      else
      ijstop =nstrong_unique
      endif
      iprev_shell=1
      ipasssize=ijstop-ijstart+1
      iresstart1=ires1+(ijstart-1)*merged_size
      iresstart2=ires2+(ijstart-1)*merged_size
c
      if (non_block) then
        ij_start_predict=(ipaser  )*ipas_siz+1
        ij_stop_predict= (ipaser+1)*ipas_siz
        if (ij_start_predict.gt.nstrong_unique) then
          ij_start_predict=-1
          ij_stop_predict =-1
        else if (ij_stop_predict.gt.nstrong_unique) then
          ij_stop_predict=nstrong_unique
        endif
        ipasssize_predict=ij_stop_predict-ij_start_predict+1
cAFil   if (iswitcher.eq.1) then
cAFil     icoefsymm=icoefsymm1
cAFil     icoefanti=icoefanti1
cAFil     icoefsymm_predict=icoefsymm2
cAFil     icoefanti_predict=icoefanti2
cAFil     iswitcher=2
cAFil   else
cAFil     icoefsymm=icoefsymm2
cAFil     icoefanti=icoefanti2
cAFil     icoefsymm_predict=icoefsymm1
cAFil     icoefanti_predict=icoefanti1
cAFil     iswitcher=1
cAFil   endif
      endif !  if (non_block) then
c
      call matconn('Cplus' ,'r',ntrisymm,ipasssize,icoefsymm)
      call matconn('Cminus','r',ntrianti,ipasssize,icoefanti)
      call elapsec(exread0)
      call secund(xread0)
      if (non_block) then
        call GetCoef(ijstart,'s',icoefsymm,nodisk,ipasssize)
        call GetCoef(ijstart,'a',icoefanti,nodisk,ipasssize)
        if (ij_start_predict.ne.-1) then
          call OrderCoef(ij_start_predict,'s',icoefsymm_predict,nodisk,
     *                   ipasssize_predict)
          call OrderCoef(ij_start_predict,'a',icoefanti_predict,nodisk,
     *                   ipasssize_predict)
        endif
      else
        call PointerCoef(ijstart,'s',icoefsymm,nodisk,ipasssize)
        call PointerCoef(ijstart,'a',icoefanti,nodisk,ipasssize)
      endif
      call elapsec(exread1)
      call secund(xread1)
      econstr=econstr+exread1-exread0
      sconstr=sconstr+xread1-xread0
c     call matconn('trans1','r',merged_size,ipasssize,iresstart1)
c     call matconn('trans2','r',merged_size,ipasssize,iresstart2)
c COMPRESSION * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
c 1 matrix addr., 
c 2 nrow =-1 if not compressed
c 3 compression info addr relative to 1, 
      do ishells=1,igroup_no
        istr_start=ishell_index(1,ishells)
        istr_stop =ishell_index(2,ishells)
        istripsize=istr_stop-istr_start+1
        iads=ishell_address_s(1,ishells)
        iada=ishell_address_a(1,ishells)
        iksymm_size=ishell_address_s(2,ishells)
        iksymm_compr_info=iads+ishell_address_s(3,ishells)
        ikanti_size=ishell_address_a(2,ishells)
        ikanti_compr_info=iada+ishell_address_a(3,ishells)
        call matconn('partial1','r',istripsize,ipasssize,ipartres1)
        call matconn('partial2','r',istripsize,ipasssize,ipartres2)
        call elapsec(exzeros0)
        call secund(xzeros0)
        if (iksymm_size.eq.-1) then
          call matconn('Kplus', 'r',ntrisymm,istripsize,iads)
        else
          if (iksymm_size.ne.0) then
          call matconn('Kplus', 'r',iksymm_size,istripsize,iads)
          call C_compress(bl(icoefsymm),   ntrisymm,    ipasssize,
     *                    bl(icoefsymm_c), iksymm_size, ipasssize,
     *                    bl(iksymm_compr_info))
          call matconn('Ccompr+','r',iksymm_size,ipasssize,icoefsymm_c)
          endif
        endif
        if (ikanti_size.eq.-1) then
          call matconn('Kminus','r',ntrianti,istripsize,iada)
        else
          if (ikanti_size.ne.0) then
          call matconn('Kminus', 'r',ikanti_size,istripsize,iada)
          call C_compress(bl(icoefanti),   ntrianti,    ipasssize,
     *                    bl(icoefanti_c), ikanti_size, ipasssize,
     *                    bl(ikanti_compr_info))
          call matconn('Ccompr-','r',ikanti_size,ipasssize,icoefanti_c)
          endif
        endif
        call elapsec(exzeros1)
        call secund(xzeros1)
        eetr=eetr+exzeros1-exzeros0
        extr=extr+xzeros1-xzeros0
        call elapsec(exmult0)
        call secund(xmult0)
        if (iksymm_size.eq.-1) then
          call matmmul2('Kplus','Cplus','partial1','t','n','n')
        else
          if (iksymm_size.eq.0) then
          call matzero('partial1')
          else
          call matmmul2('Kplus','Ccompr+','partial1','t','n','n')
          endif
        endif
        if (ikanti_size.eq.-1) then
          call matmmul2('Kminus','Cminus','partial2','t','n','n')
        else
          if (ikanti_size.eq.0) then
          call matzero('partial2')
          else
          call matmmul2('Kminus','Ccompr-','partial2','t','n','n')
          endif
        endif
        call elapsec(exmult1)
        call secund(xmult1)
        etrace=etrace+exmult1-exmult0
        strace=strace+xmult1-xmult0
        if (ikanti_size.gt.0) call matdisc('Ccompr-')
        if (ikanti_size.ne.0) call matdisc('Kminus')
        if (iksymm_size.gt.0) call matdisc('Ccompr+')
        if (iksymm_size.ne.0) call matdisc('Kplus')
        call matdisc('partial2')
        call matdisc('partial1')
        call elapsec(exzeros0)
        call secund(xzeros0)
c                          source,        dimensions
        call alloc_results(bl(ipartres1), istripsize,ipasssize,
c                          target,        dimensions
     *                     bl(iresstart1),merged_size,ipasssize,
     *                     istr_start) ! shift
c
c                          source,        dimensions
        call alloc_results(bl(ipartres2), istripsize,ipasssize,
c                          target,        dimensions
     *                     bl(iresstart2),merged_size,ipasssize,
     *                     istr_start) ! shift
        call elapsec(exzeros1)
        call secund(xzeros1)
        eetr=eetr+exzeros1-exzeros0
        extr=extr+xzeros1-xzeros0
      enddo
c END COMPRESSION * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

c     call matdisc('trans2')
c     call matdisc('trans1')
      call matdisc('Cminus')
      call matdisc('Cplus')
      enddo pass ! do ipaser=1,npass
c
      call dynamic_retmem(6)
c
c  Sort the half-transformed integrals
      call matconn('trans1','r',merged_size,nstrong_unique,ires1)
      call matconn('trans2','r',merged_size,nstrong_unique,ires2)
      call elapsec(exsort0)
      call secund(xsort0)
      do ik=1,merged_size
          icf=ikindex(1,ik)
          kcf=ikindex(2,ik)
          ij=0
          ij_compr=0
          do i=1,nval
            do j=1,i
              ij=ij+1
              if (list(ij).eq.1.and.isympairs(ij).ne.0) then 
              ij_compr=ij_compr+1
              call matelem('trans1',ik,ij_compr,elem1)
              call matelem('trans2',ik,ij_compr,elem2)
              if (i.eq.j) then 
                halftra(i,j)=elem1
                cycle
              endif
              if (icf.eq.kcf) then 
                halftra(i,j)=elem1
                halftra(j,i)=elem1
                cycle
              endif
              halftra(i,j)=elem1+elem2
              halftra(j,i)=elem1-elem2
              else  !  if (list(ij).eq.1)
              halftra(i,j)=0.0d0
              halftra(j,i)=0.0d0
              endif !  if (list(ij).eq.1)
            enddo
          enddo
c         do i=1,nval
c         do j=1,nval
c         write(78,'(F70.16)') halftra(i,j)
c         enddo
c         enddo
      if (igranulesize.le.1) then
          call BinSort_tj( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1,  bin8, afflag, afname,
     3                  islvid,byt8)
      else
          call BinSortGr( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1, bin8,  afflag, afname,
     3                  islvid,byt8,
     *          igranules,igranulesize,igran2pair,ipair2gran)
      endif
      enddo
      call elapsec(exsort1)
      call secund(xsort1)
      eort=eort+exsort1-exsort0
      sort=sort+xsort1-xsort0
      call dynamic_retmem(2)
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine SymAntiBuild(source,symmetric,antisymmetric,ncf,scr)
c This subroutine resolves given matrix 'source' into its symmetric and 
c antisymmetric components. Note xmult, which ensures correct results when
c such matrix is then mutipled by its counterpars which also should have the 
c same xmult supplied (it is built in KExtrConstr).
c The EEO screening matrix is also constructed here.
c Parameters:
cI  source     - the source matrix
cO  symmetric  - symmetric part of source, lower triangle
cO  antisymmetric- ditto for anti
cI  ncf        - number of contracted functions
cO  scr        - screening matrix, first ncf*(ncf+1)/2 contains max elements
c                from symm matrix, the rest of it contains max elements from
c                anti matrix, total size is ncf*ncf
      implicit real*8 (a-h,o-z)
      dimension source(ncf,ncf),symmetric(*),antisymmetric(*),scr(*)
      parameter (zero=0.0d0,one=1.0d0,two=2.0d0)
c
      ntrisymm=ncf*(ncf+1)/2
      ij=0
      ij_anti=ntrisymm
      xmult=sqrt(2.0d0)/2.0d0
      do j=2,ncf
        do i=1,j-1
          ij=ij+1
          ij_anti=ij_anti+1
          symmetric(ij)     = (source(i,j)+source(j,i))*xmult
          if (dabs(symmetric(ij)).gt.scr(ij)) 
     *                             scr(ij)=dabs(symmetric(ij))
          antisymmetric(ij) = (source(i,j)-source(j,i))*xmult
          if (dabs(antisymmetric(ij)).gt.scr(ij_anti)) 
     *                             scr(ij_anti)=dabs(antisymmetric(ij))
        enddo
      enddo
c
      do i=1,ncf
        ij=ij+1
        symmetric(ij)       = source(i,i)
        if (dabs(symmetric(ij)).gt.scr(ij)) 
     *                           scr(ij)=dabs(symmetric(ij))
      enddo
      if (ij_anti.gt.ncf*ncf) STOP 'Error rifhvh55'
      end
c==============================================================================
      subroutine store(ijstart,ijstop,ncf,nval,ibaseaddress)
c This subroutine stores bunch of amplitudes into given address.
c
c Parameters:
cI  ijstart       - starting ij pair identifying first T matrix
cI  ijstop        - closing ij pair identifying last T matrix
cI  ncf           - number of contracted functions
cI  nval          - number of correlated pairs
cO  ibaseaddress  - the bl location where all amplitudes are deposited.
      use memory
      implicit none
      integer ijstart,ijstop,ncf,ibaseaddress,nval
      integer i,j,ij,iaddr,icoeff
      ij=0
      outer: do i=1,nval
        do j=1,i
          ij=ij+1
          if (ij.lt.ijstart) cycle
          if (ij.gt.ijstop) exit outer
          iaddr=ibaseaddress+(ij-ijstart)*ncf*ncf
          call CoefRead('tc',i,j,icoeff)
          call tfer(bl(icoeff),bl(iaddr),ncf*ncf)
        enddo
      enddo outer
      end
c==============================================================================
      subroutine storeK(irecadrx, ndiskx,lbinx, thresh,       ijstart,
     *                  ijstop,   ncf,   nval,  ibaseaddress, byt8,
     *                  nmo,      vorb,  nbf)
c This subroutine stores bunch of exchange integrals into given memory 
c location.
c Parameters:
cI  irecadrx      - record info for Kext
cI  ndiskx        - logical unit of disk number where Kext are stored
cI  lbinx         - length of bin in file ndiskx
cI  thresh        - integral threshold
cI  ijstart       - starting ij pair identifying first T matrix
cI  ijstop        - closing ij pair identifying last T matrix
cI  ncf           - number of contracted functions
cI  nval          - number of correlated pairs
cO  ibaseaddress  - the bl location where all amplitudes are deposited.
cI  byt8          - logical viariable if .true. bins are stored in real*8 form
cI  nmo           - number of occupied HF orbitals
cI  vorb          - if true all matrices are kept in MO form as default, 
c                   which impacts also their dimension
      use memory
      implicit none
      integer irecadrx(*),ndiskx,lbinx,ijstart,ijstop,ncf,ibaseaddress 
      integer nval,nmo,nbf
      real*8 thresh
      integer i,j,ij,iaddr,npairs,idimen
      logical byt8,vorb
      ij=0
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
      npairs=nval*(nval+1)/2
      outer: do i=1,nval
        do j=1,i
          ij=ij+1
          if (ij.lt.ijstart) cycle
          if (ij.gt.ijstop) exit outer
          iaddr=ibaseaddress+(ij-ijstart)*idimen*idimen
          call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(iaddr))
        enddo
      enddo outer
      end
c==============================================================================
      subroutine CalcList(list,   nfirst,   nval,     nmo,     loca,
     *                    inumber,weaksepar,distsepar,energies,iEnergy)
c This subroutine builds array 'list' which idicates which orbital pairs
c can be classified as strong (+1), weak (0) or distant (-1). This 
c classification is based either on energy or space separation criterion,
c depending on 'iEnergy' value. The values of these criteria are given in
c 'weaksepar' and 'distsepar' variables.
c Parameters:
cO  list         - the integer array where for each ij element the number is 
c                  assigned. 
c                  if list(ij)==+1 the pair is strong, 
c                  if list(ij)== 0 the pair is weak, 
c                  if list(ij)==-1 the pair is distant.
cI  nfirst       - the number of first correlated orbital, would be 1 in the
c                  'nofrozen' mode
cI  nval         - number of correlated orbitals
cI  nmo          - number of occupied HF orbitals.
cI  loca         - logical variable indicating whether pair ommision will be 
c                  performed. If .false. no pair selection is enabled and 
c                  all pairs are considered to be strong.
cO  inumber      - number of strong pairs
cI  weaksepar    - this criterion separates strong and weak pairs
cI  distsepar    - this criterion separates weak and distant pairs
cI  energies     - pair energies used when energy criterion is chosen
cI  iEnergy      - iEnergy==0  => Distance pair selection
c                  otherwise   => Energy   pair selection
      use memory
      implicit none
      integer list(*)
      integer nfirst,nval,nmo,inumber,iEnergy
      logical loca
      real*8 weaksepar,distsepar,energies(*)
      integer mataddr
      integer np4,ij,i,j,ix,iy,iz,i2,iabsol,jabsol
      real*8  xi,yi,zi,xj,yj,zj,r2i,r2j,radiusi,radiusj,distanceij
c     parameter (weaksepar=1.4d0,distsepar=4.0d0) ! Bohr!
      call mmark
      call dynamic_mmark
      call matmark
      if (loca) then
      np4=4
c
      call dynamic_matdef('xmo','d',nmo,nmo)  ! average (i|x|i)   for all occupied i
      call dynamic_matdef('ymo','d',nmo,nmo)  ! average (i|y|i)
      call dynamic_matdef('zmo','d',nmo,nmo)  ! average (i|z|i)
      call dynamic_matdef('r2mo','d',nmo,nmo) ! average (i|r**2|i) 
c
      call matread('xmo',np4,'xcor_loc')
      call matread('ymo',np4,'ycor_loc')
      call matread('zmo',np4,'zcor_loc')
      call matread('r2mo',np4,'rr_loc')
c
      ix=mataddr('xmo')
      iy=mataddr('ymo')
      iz=mataddr('zmo')
      i2=mataddr('r2mo')
c
      ij=0
      inumber=0
      if (iEnergy.eq.0) then
              print *, 'Distance pair selection.'
      else
              print *, 'Energy   pair selection.'
      endif
c
      do i=1,nval
        iabsol=nfirst+i-1
        xi=bl(ix+iabsol-1)
        yi=bl(iy+iabsol-1)
        zi=bl(iz+iabsol-1)
        r2i=bl(i2+iabsol-1)
        radiusi=sqrt(r2i-xi**2-yi**2-zi**2)
c
        do j=1,i
          ij=ij+1
          jabsol=nfirst+j-1
          xj=bl(ix+jabsol-1)
          yj=bl(iy+jabsol-1)
          zj=bl(iz+jabsol-1)
          r2j=bl(i2+jabsol-1)
          radiusj=sqrt(r2j-xj**2-yj**2-zj**2)
          distanceij=sqrt((zj-zi)**2+(yj-yi)**2+(xj-xi)**2)
c
          if (iEnergy.eq.0) then
            if (distanceij.gt.(radiusi+radiusj+distsepar)) then
              list(ij)=-1
              print '(A,2I2,A)','Pair ',i,j,' distant'
            else if (distanceij.gt.(radiusi+radiusj+weaksepar)) then
              list(ij)=0
              print '(A,2I2,A)','Pair ',i,j,' weak'
            else
              list(ij)=1
              inumber=inumber+1
              print '(A,2I2,A)','Pair ',i,j,' strong'
            endif
          else
            if (dabs(energies(ij)).lt.distsepar) then
              list(ij)=-1
              print '(A,2I2,A)','Pair ',i,j,' distant'
            else if (dabs(energies(ij)).lt.weaksepar) then
              list(ij)=0
              print '(A,2I2,A)','Pair ',i,j,' weak'
            else
              list(ij)=1
              inumber=inumber+1
              print '(A,2I2,A)','Pair ',i,j,' strong'
            endif
          endif
c
        enddo
c
      enddo
c
      else
c
      ij=0
      inumber=0
      do i=1,nval
        do j=1,i
        ij=ij+1
        list(ij)=1
        inumber=inumber+1
        enddo
      enddo
c
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine filllist(list,npairs)
      implicit none
      integer list(*),npairs
      integer i
      do i=1,npairs
      list(i)=1
      enddo
      end
c==============================================================================
      subroutine CINorm(icoeffAO,xnormij,ncf,equal,nmo,vorb,nbf)
c This subroutine calculates the contribution to CI norm from all T(ij)
c amplitudes. The result is roughly equivalent to \sum_{i}(T_i*T_i), 
c where i stands for mu-lam composite index.
c Parameters:
cI  icoeffAO     - bl location of coefficients
cO  xnormij      - the result (pair squared norm)
cI  ncf          - number of cotracted functions
cI  equal        - this logical variable is true if i==j
cI  nmo          - number of occupied HF orbitals
cI  vorb         - if equal to .true. it means that all matrices are representedc                  in MO instead AO basis, also their dimension is smaller
c
      use memory
      implicit none
      integer icoeffAO,ncf,nmo,nbf
      REAL*8 xnormij
      REAL*8 zero
      data zero /0.0d0/
      logical equal,vorb
      integer idimen
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('MOtr','q',idimen,idimen)
      call dynamic_matdef('work11','q',idimen,idimen)
      if (vorb) then
        call matconn('MO','q',idimen,idimen,icoeffAO)
      else
        call matconn('AO','q',ncf,ncf,icoeffAO)
        call dynamic_matdef('MO','q',ncf,ncf)
        call matsimtr('AO','invcano','MO')
      endif
      call matcopy('MO','MOtr')
      call matpose('MOtr')
c     do j=1,ncf
c       do i=1,ncf
c         if (i.le.nmo .or. j.le.nmo) then
c           call mateset('exchMO',i,j,zero)
c           call mateset('coeffMO',i,j,zero)
c           call mateset('resiMO',i,j,zero)
c         endif
c       enddo
c     enddo
c
      call matzero('work11')
      call matadd1('MO',4.0d0,'work11')
      call matadd1('MOtr',-2.0d0,'work11')
      call matprodtr('MO','work11',xnormij)
      if (equal) xnormij=xnormij/2.0d0
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine fill22(xmatrix,idimen)
      implicit real*8 (a-h,o-z)
      dimension xmatrix(idimen)
      do i=1,idimen
      xmatrix(i)=1.0d20
      enddo
      end
c==============================================================================
      subroutine ijtriangle(i,j,ij)
      implicit none
      integer i,j,ij
      if (i.ge.j) then
      ij=i*(i-1)/2+j
      else
      ij=j*(j-1)/2+i
      endif
      end
c==============================================================================
      subroutine differ_generator(ndiskold,ndisknew,ndiskdiffer,ncf,
     *                            npairs,ccsd)
c Temporarily not used
      use memory
      implicit none
c
      integer ndiskold,ndisknew,ndiskdiffer,ncf,npairs
      logical ccsd
c
      integer iold,inew,idif,i,j,ii,iorb,jorb,imap,jmap,ivecti,ivectj
      integer mataddr
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('old','q',ncf,ncf)
      call dynamic_matdef('new','q',ncf,ncf)
      call dynamic_matdef('dif','q',ncf,ncf)
      iold=mataddr('old')
      inew=mataddr('new')
      idif=mataddr('dif')
      call dynamic_matdef('vecti','r',ncf,1)
      ivecti=mataddr('vecti')
      call dynamic_matdef('vectj','r',ncf,1)
      ivectj=mataddr('vectj')
c
      do ii=1,npairs
        call reader(ndiskold,ncf,ii,bl(iold))
        call reader(ndisknew,ncf,ii,bl(inew))
        do j=0,ncf*ncf-1
          bl(idif+j)=0.5d0*(bl(inew+j)-bl(iold+j))
        enddo
c
c     print *,'Coupled: ',ccsd
        if (ccsd) then
          iorb=imap(ii)
          jorb=jmap(ii)
          call readsingles(iorb,bl(ivecti))
          call readsingles(jorb,bl(ivectj))
          call matmmul2('vecti','vectj','old','n','t','n')
          call readsinglesnew(iorb,bl(ivecti))
          call readsinglesnew(jorb,bl(ivectj))
          call matmmul2('vecti','vectj','new','n','t','n')
          do j=0,ncf*ncf-1
            bl(idif+j)=bl(idif+j)+bl(inew+j)-bl(iold+j)
          enddo
        endif
c
        call writer(ndiskdiffer,ncf,ii,bl(idif))
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine calculate_omits(ncf,npairs,ndiskr,ndiska,threshold,
     *                           percent,percentrel)
c Temporarily not used
      use memory
      implicit none
      integer ncf,npairs,ndiskr,ndiska
      real*8 percent,percentrel,threshold
c
      integer i,imatrix,iabsol,nrow,ncol,nrowrel,ncolrel
c
      percent=0.0d0
      percentrel=0.0d0
      call dynamic_getmem(ncf*ncf,imatrix)
      call dynamic_getmem(ncf*ncf,iabsol)
      do i=1,npairs
c       call reader(ndiskr,ncf,i,bl(imatrix))
        call reader(ndiska,ncf,i,bl(iabsol))
        call ZeroRowCol(bl(imatrix),bl(iabsol),threshold,ncf,nrow,ncol,
     *                  nrowrel,ncolrel)
        percent=percent+1.0d0*nrow*ncol/(1.0d0*ncf*ncf)
c       percentrel=percentrel+1.0d0*nrowrel*ncolrel/(1.0d0*ncf*ncf)
      enddo
      percent=percent/npairs
      percentrel=percentrel/npairs
      call dynamic_retmem(2)
      end
c==============================================================================
      subroutine calculate_zero_rows(xmatrix,idim1,idim2,screen,nrow)
c Temporarily not used
      implicit none
      integer idim1,idim2,nrow
      real*8 xmatrix(idim1,idim2),screen(idim1)
c
      integer i,j
c
      nrow=0
      do i=1,idim1
        do j=1,idim2
          if (dabs(xmatrix(i,j)*screen(i)).ge.1.d0) then
            nrow=nrow+1
            exit
          endif
        enddo
      enddo
      end
c==============================================================================
      subroutine ZeroRowCol(xmatrix,xrel,thresh,ncf,nrow,ncol,
     *                      nrowrel,ncolrel)
c Temporarily not used
      implicit none
      integer ncf,nrow,ncol,nrowrel,ncolrel
      real*8 xmatrix(ncf,ncf),xrel(ncf,ncf),thresh
c
      logical zerocol,zerocolrel
      integer i,j
c
      nrow=0
      ncol=0
      nrowrel=0
      ncolrel=0
      do i=1,ncf
        zerocol=.true.
        zerocolrel=.true.
        do j=1,ncf
c       if (dabs(xmatrix(i,j)).gt.thresh) zerocol=.false.
c       if (dabs(xmatrix(i,j)/xrel(i,j)).gt.thresh) zerocolrel=.false.
        if (dabs(xrel(i,j)).gt.thresh) zerocol=.false.
        enddo
        if (.not.zerocol)    nrow=nrow+1
        if (.not.zerocolrel) nrowrel=nrowrel+1
      enddo
      do i=1,ncf
        zerocol=.true.
        zerocolrel=.true.
        do j=1,ncf
        if (dabs(xmatrix(j,i)).gt.thresh) zerocol=.false.
        if (dabs(xmatrix(j,i)/xrel(j,i)).gt.thresh) zerocolrel=.false.
        enddo
        if (.not.zerocol) ncol=ncol+1
        if (.not.zerocolrel) ncolrel=ncolrel+1
      enddo
c     print *, 'thresh,nrow,ncol: ',thresh,nrow,ncol
      end
c==============================================================================
      subroutine compactor(xmatrix,thresh,idim1,idim2,irow,icol,
     *                     nrow,ncol)
c Temporarily not used
      implicit none
      integer idim1,idim2,icol(*),irow(*),ncol,nrow
      real*8 xmatrix(idim1,idim2),thresh
c
      integer i,j,iold,jold
c
c
c  compact columns:
      ncol=0
      do j=1,idim2
c  check if there is nonzero element:
        do i=1,idim1
        if (dabs(xmatrix(i,j)).gt.thresh) then
          ncol=ncol+1
          icol(ncol)=j
          exit
        endif
        enddo
      enddo
c  compact rows:
      nrow=0
      do i=1,idim1
c  check if there is nonzero element:
        do j=1,idim2
        if (dabs(xmatrix(i,j)).gt.thresh) then
          nrow=nrow+1
          irow(nrow)=i
          exit
        endif
        enddo
      enddo
      if (nrow.ge.1.and.ncol.ge.1) then
c The last compactation, squeeze the free spaces in columns, remember,
c the final matrix cannot be addressed as a xmatrix(ncf,ncf)!
ctest   call Squeeze(xmatrix,ncf,irow,nrow,ncol)
        call Squeeze('both',xmatrix,idim1,idim2,irow,icol,nrow,ncol)
      endif
      end
c==============================================================================
      subroutine EEORowCompactor(xmatrix, idim1,     idim2,screen,irow,
     *                           nrow,    icompraddr,itop, tot_xmax)
c This subroutine compacts matrix xmatrix (rows only) and stores compression
c info in the space released by compressed matrx. If this space is not large
c enough, the matrix is not compressed. It returns two addresses: one is
c the index where the compression info is stored, second is the address of 
c free space, where next matrix can be stored by caller. Matrix 'xmatrix' is
c compressed in place.
cIO xmatrix      - matrix to be compressed
cI  idim1, idim2 - dimensions of the matrix
cI  screen       - screening matrix, the dimension of it should be equal 
c                  to idim1 (only rows are compacted)
cO  irow         - compression info array
cO  nrow         - number of rows after compression, -1 = no compression
cO  icompraddr   - the address relative to xmatrix start, where compression 
c                  will be copied (from temporary irow array), if 0 the whole 
c                  matrix is null, if -1: no compression at all
cO  itop         - gives the address of top of the (matrix)+(compression info), 
c                  this address is then used as starting address for storing
c                  next matrix, itop may be 0 which means that the whole matrix
c                  is zeroed (null matrix)
cI  tot_xmax     - this is the maximum residuum element from previous iteration
      use kinds
      implicit none
      integer idim1,idim2,nrow,icompraddr,itop
      integer*4 irow(*)
      real*8 xmatrix(idim1,idim2),screen(*),tot_xmax
c
      integer i,j,ii,istart,inewcolshift,ioldcolshift,imat_size,iold
      integer icompr_size,icompr_info_size,idim2_copy,icol(10),modifier
      real*8 scrj,x0,x1,ex0,ex1,percent_factor,thresh
      logical zero,first
      real*8             strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
c     integer     intsize,iacc,icache,memreal
c     common /cpu/intsize,iacc,icache,memreal
      real*8 zeroed,zeroed_max
      common /zero_in_EEO/ zeroed,zeroed_max
c
c
      nrow=0
      call set_thresh(tot_xmax,thresh)
      call secund(x0)
      call elapsec(ex0)
      do i=1,idim1
c  check if there is nonzero element:
c It could be done more efficiently, but I check now whether zeroing ALL
c integrals below certain threshold is better that zeroing only zero
c rows. My preliminary experience shows that zeroing all integrals could
c give better, more stable & consistent results.
        scrj=screen(i)
        first=.true.
        do j=1,idim2
        if (dabs(xmatrix(i,j)*scrj).gt.thresh) then
          if (first) then
            nrow=nrow+1
            irow(nrow)=i
            first=.false.
          endif
        else
          xmatrix(i,j)=0.0d0
          zeroed_max=zeroed_max+1.0d0
        endif
        enddo
      enddo
      call secund(x1)
      call elapsec(ex1)
      c1=c1+x1-x0
      e1=e1+ex1-ex0
      imat_size  =idim1*idim2*8
      icompr_size=nrow *idim2*8
      icompr_info_size=nrow*4
c     percent_factor=dble(idim2)/(dble(idim2)+12.0d0)
c     percent_factor=dble(idim2)/(dble(idim2)+10.0d0)
c     modifier, see subr. definition
      modifier=1
      call calculate_percent(modifier,idim2,percent_factor)
c  if compression is below percent_factor
      if (dble(nrow)/dble(idim1).gt.percent_factor) then 
        nrow=-1
        itop=idim1*idim2
        icompraddr=-1
        return
      else if ((icompr_info_size+icompr_size+20).gt.imat_size) then
        nrow=-1
        itop=idim1*idim2
        icompraddr=-1
        return
      else !matrix is compressed enough and compression info fits free space
        call secund(x0)
        call elapsec(ex0)
        if (nrow.gt.0) then
        idim2_copy=idim2
        call Squeeze('rows',xmatrix,idim1,idim2,irow,icol,
     *                   nrow,idim2_copy)
c       copy compression info into the free space fried by old matrix
        icompraddr=nrow*idim2
        call copy_compression_info(xmatrix(icompraddr+1,1),irow,nrow)
        itop=icompraddr+nrow/i4size+1
        else
          icompraddr=0
          itop=icompraddr
        endif
        call secund(x1)
        call elapsec(ex1)
        c2=c2+x1-x0
        e2=e2+ex1-ex0
      endif
c
      return ! the rest of the code is not executed now
c
c zero integrals: test only, should be squeezed
      istart=1
      outer: do i=1,idim1
        zero=.true.
        do ii=istart,nrow
        if (irow(ii).eq.i) then
          istart=ii+1
          zero=.false.
          exit
        else if (irow(ii).gt.i) then
          istart=ii
          exit
        endif
        enddo
        if (zero) then
          do j=1,idim2
            xmatrix(i,j)=0.0d0
          enddo
        endif
      enddo outer
      end
c==============================================================================
      subroutine RowCompactor(xmatrix,thresh,idim1,idim2,irow,icol,
     *                     nrow,ncol)
c Temporarily not used
      implicit none
      integer idim1,idim2,icol(*),irow(*),ncol,nrow
      real*8 xmatrix(idim1,idim2),thresh
c
      integer i,j,iold,jold
c
c
c  compact columns:
      ncol=idim2
c  compact rows:
      nrow=0
      do i=1,idim1
c  check if there is nonzero element:
        do j=1,idim2
        if (dabs(xmatrix(i,j)).gt.thresh) then
          nrow=nrow+1
          irow(nrow)=i
          exit
        endif
        enddo
      enddo
      if (nrow.ge.1.and.ncol.ge.1) then
c The last compactation, squeeze the free spaces in columns, remember,
c the final matrix cannot be addressed as a xmatrix(ncf,ncf)!
ctest   call Squeeze(xmatrix,ncf,irow,nrow,ncol)
        call Squeeze('rows',xmatrix,idim1,idim2,irow,icol,nrow,ncol)
      else if (nrow.eq.0.or.ncol.eq.0) then
        print *, 'WARNING!!!   nrow equal 0!'
      endif
      end
c==============================================================================
      subroutine Squeeze(WhatToDo, xmatrix, idim1,  idim2,  irow,
     *                   icol,     nrow,    ncol)
c This subroutine squeezes the matrix according to the info provided in
c irow,icol & nrow,ncol variables. It might squeeze only rows, only columns
c or both, depending on the 'WhatToDo' variable value. Squeezing is performed 
c "in place". If you need to squeeze both rows and columns it is better to 
c do it using 'both' option, not by use of consecutive 'colu' and then 'rows'. 
c The reason is that compression in 'both' case is more efficient - rows
c and columns are compressed at once in one loop.
c Parameters:
cI  WhatToDo     - Character variable, might be 'colu', 'rows', 'both', which 
c                  means: compress columns, rows or both rows & columns.
cIO xmatrix      - matrix to be compressed
cI  idim1        - row dimension of input matrix
cI  idim2        - columns dimension of input matrix
cI  irow         - this array contains mapping of new row numbers of
c                  compressed matrix into old (original) row number of source
c                  xmatrix
cI  nrow         - number of rows of compressed matrix
cI  
      implicit none
      integer idim1,idim2,nrow,ncol
      integer*4 irow(*),icol(*)
      real*8 xmatrix(*)
      integer i,j,iold,jold,inewcolshift,ioldcolshift
      character*(*) WhatToDo   ! can be: 'colu', 'rows', 'both'
      if (ncol.eq.0.or.nrow.eq.0) return
      if (WhatToDo.eq.'both') then
        if (ncol.eq.idim2.and.nrow.eq.idim1) return
        do j=1,ncol
          jold=icol(j)
          inewcolshift=nrow*(j-1)
          ioldcolshift=idim1*(jold-1)
          do i=1,nrow
            iold=irow(i)
            xmatrix(i+inewcolshift)=xmatrix(iold+ioldcolshift)
          enddo
        enddo
      else if (WhatToDo.eq.'colu') then
        if (ncol.eq.idim2) return
        do j=1,ncol
          jold=icol(j)
          inewcolshift=nrow*(j-1)
          ioldcolshift=idim1*(jold-1)
          do i=1,nrow
            xmatrix(i+inewcolshift)=xmatrix(i+ioldcolshift)
          enddo
        enddo
      else if (WhatToDo.eq.'rows') then
        if (nrow.eq.idim1) return
        do j=1,ncol
          inewcolshift=nrow*(j-1)
          ioldcolshift=idim1*(j-1)
          do i=1,nrow
            iold=irow(i)
            xmatrix(i+inewcolshift)=xmatrix(iold+ioldcolshift)
          enddo
        enddo
      else
        STOP 'Should be: colu, rows or both'
      endif
ctest do j=1,ncol
ctest   do i=1,nrow
ctest   xmatrix((j-1)*nrow+i)=xmatrix((j-1)*ncf+i)
ctest   enddo
ctest enddo
      end
c==============================================================================
      subroutine SecondIsRightCopyCom(xmatrixs,xmatrixd,idim1,idim2,
     *                                irow,nrow)
c Not used, probably more general 'squeeze' should be used in the future.
      implicit none
      integer idim1,idim2,irow(*),nrow
      real*8 xmatrixs(idim1,idim2),xmatrixd(nrow,idim2)
      integer i,j,iold
c
      if (nrow.eq.idim1) return
      do j=1,idim2
        do i=1,nrow
          iold=irow(i)
          xmatrixd(i,j)=xmatrixs(iold,j)
        enddo
      enddo
      end
c==============================================================================
c not used, but kept
      subroutine BigMuliplyEEO(nval,    ncf,   npairs,  lbin,   ndisk,
     1                         halftra, thresh,icounter,irecord,indxbin,
     2                         ibin4,   ibin1, bin8,    afflag, reset,
     3                         islvid,  byt8,  lmp2int, ndiskampl, ! ->
     4                                              max_merg_shell_size,
     5                         max_amplitudes_size,ishellsaddr,    ! ->
     6                                 iamplitaddr,iresultaddr,ikindex)
c Not used at present, it was used once for testing purposes (EEO big option)
      use memory
      implicit none
      integer nval,ncf,npairs,lbin,ndisk
      integer islvid,lmp2int,ishellsiz,ndiskampl
      integer icounter(*),irecord(*),ikindex(2,*)
      integer max_merg_shell_size,max_amplitudes_size,ishellsaddr
      integer iamplitaddr,iresultaddr
      real*8 thresh,bin8(npairs,lbin),halftra(nval,nval)
      integer*2 indxbin(2,npairs,lbin)
      integer*4 ibin4(npairs,lbin)
      integer*1 ibin1(npairs,lbin)
      logical byt8,afflag,reset
      character*2 afname
      real*8 strace,sconstr,sort,extr,const,t0,t1,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
c
      real*8 element,threshint,xneglect,resneg(20),threshold
      integer nval2,ij,ik,i,j,icf,kcf,ibase,ll,icoeff,iTbig
      integer mataddr,isinglvec,ncf2,k,l,mm
      integer irow(ncf*ncf),icol(nval*nval)
      integer ijpair,nrow,ncol,itmp,nrowexp,iKbigSq
      logical ishere,calculate_omits,iamzero
c
      integer iTpass,iTnoPass,iTstart,iTstop,iTsize,iresult
c
      calculate_omits=.false.
      call mmark
      call dynamic_mmark
      call matmark
      nval2=nval*nval
      ncf2=ncf*ncf
      xneglect=1.d-8
      call matconn('Kbig','r',ncf2,max_merg_shell_size,ishellsaddr)
c     Calculate number of T passess
      iTpass=nval2/max_amplitudes_size
      if (max_amplitudes_size*iTpass.lt.nval2) iTpass=iTpass+1
c     if (dabs(strace).le.1.d-5) then Statistical purposes only
c       calculate_omits=.true.
c         do itmp=1,20
c           resneg(itmp)=0.0d0
c         enddo
c     endif
c
      do iTnoPass=1,iTpass
        iamzero=.false.
        iTstart=(iTnoPass-1)*max_amplitudes_size+1
        iTstop=iTnoPass*max_amplitudes_size
        if (iTstop.gt.nval2) iTstop=nval2
        iTsize=iTstop-iTstart+1
        call secund(t0)
        call storeEEOampl(ndiskampl,iTstart,iTstop,nval,ncf,
     *                                             bl(iamplitaddr))
c
c       if (calculate_omits) then ! Statistical purposes only
c         do itmp=10,18
c         threshold=1.d1**(-1.d0*itmp/2.d0)
c         call calculate_zero_rows(ncf,npairs,bl(iamplitaddr),ncf2,
c    *                             iTsize,threshold,nrowexp)
c         resneg(itmp)=resneg(itmp)+nrowexp*1.d0/(ncf2*1.d0)
c         enddo
c       endif
c
        call secund(t1)
        extr=extr+t1-t0
        iresult=iresultaddr+(iTstart-1)*max_merg_shell_size
        call matconn('result','r',max_merg_shell_size,iTsize,iresult)
csqu    call RowCompactor(bl(iamplitaddr),xneglect,ncf2,iTsize,
csqu *                    irow,icol,nrow,ncol)
csqu    if (nrow.eq.0) then
csqu      iamzero=.true.
csqu      goto 444
csqu    endif
        call matconn('Tbig','r',ncf2,iTsize,iamplitaddr)
csqu    if (nrow.ne.ncf2) then
csqu    call dynamic_matdef('KbigSq','r',nrow,max_merg_shell_size)
csqu    iKbigSq=mataddr('KbigSq')
csqu    call SecondIsRightCopyCom(bl(ishellsaddr),bl(iKbigSq),ncf2,
csqu *                            max_merg_shell_size,irow,nrow)
csqu    endif
csqu    call secund(t0)
csqu    if (nrow.ne.ncf2) then
csqu    call matmmul2('KbigSq','Tbig','result','t','n','n')
csqu    else
csqu    call matmmul2('Kbig','Tbig','result','t','n','n')
csqu    endif
        call matmmul2('Kbig','Tbig','result','t','n','n')
        call secund(t1)
        strace=strace+t1-t0
csqu    if (nrow.ne.ncf2) call dynamic_matrem('KbigSq')
        call matdisc('Tbig')
c444    continue
csqu    if (iamzero) call matzero('result')
        call matdisc('result')
      enddo
c     if (calculate_omits) then
c       do itmp=10,18
c       threshold=1.d1**(-1.d0*itmp/2.d0)
c       resneg(itmp)=resneg(itmp)/(iTpass*1.d0)
c       print '(A,E10.2E2,A,F10.3)', 'Threshold: ',threshold,
c    *                               ' neglect: ',resneg(itmp)
c       enddo
c     endif
      call matconn('results','r',max_merg_shell_size,nval2,iresultaddr)

      do ik=1,max_merg_shell_size
        icf=ikindex(1,ik)
        kcf=ikindex(2,ik)
c  Dense added:
c         ishere=.false.
c         do mm=1,ncol
c           if (icol(mm).eq.ik) then
c           ishere=.true.
c           exit
c           endif
c         enddo
c         if (.not.ishere) cycle
c  Dense added END
          ij=0
          do i=1,nval
            do j=1,i
              ij=ij+1
              call matelem('results',ik,ij,element) ! dense change!
              halftra(i,j)=element
              if (i.ne.j) then
                ij=ij+1
                call matelem('results',ik,ij,element)  ! dense change!
                halftra(j,i)=element
              endif
            enddo
          enddo
          call BinSort_tj( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1, bin8,  afflag, afname,
     3                  islvid,byt8)
      enddo
      call matdisc('results')
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
c not used
      subroutine AmplitudesGenerator(ncf,nval,ndiskampl,ndiskres)
c Not used at present, it was used once for testing purposes (EEO big option)
      use memory
      implicit none
      integer ncf,nval,ndiskampl,ndiskres
c
      integer ij,i,j,icoeff,isinglvec,mataddr,iolds,idiffe
      integer icallno
      data icallno/0/
      save icallno
c
      icallno=icallno+1
c
      call dynamic_matdef('singlvec','r',ncf,1)
      isinglvec=mataddr('singlvec')
      call dynamic_matdef('olds','q',ncf,ncf)
      iolds=mataddr('olds')
      call dynamic_matdef('diffe','q',ncf,ncf)
      idiffe=mataddr('diffe')
c
      ij=0
      do i=1,nval
        do j=1,i
          ij=ij+1
          call CoefRead('tc',i,j,icoeff)
c Singles part of EEO:
          call matconn('coeff','q',ncf,ncf,icoeff)
          call readsingles(j,bl(isinglvec))
          call matsub('singlocc','occu',i,i)
          call matmmul2('singlocc','singlvec','coeff','n','t','a')
          call matdisc('singlocc')
c
          call readsingles(i,bl(isinglvec))
          call matsub('singlocc','occu',j,j)
          call matmmul2('singlvec','singlocc','coeff','n','t','a')
          call matdisc('singlocc')
          if (icallno.gt.1) call reader(ndiskampl,ncf,ij,bl(iolds))
          call writer(ndiskampl,ncf,ij,bl(icoeff))
          call matcopy('coeff','diffe')
          if (icallno.gt.1) call matadd1('olds',-1.d0,'diffe')
          call ZeroRelElements(bl(icoeff),bl(idiffe),ncf)
          call writer(ndiskres,ncf,ij,bl(idiffe))
          call matdisc('coeff')
        enddo
      enddo
      call dynamic_matrem('diffe')
      call dynamic_matrem('olds')
      call dynamic_matrem('singlvec')
c for non-differential builder:
      i=ndiskampl
      ndiskampl=ndiskres
      ndiskres=i
      end
c==============================================================================
      subroutine ZeroRelElements(xabs,xrel,idim1)
c Not used at present, it was used once for testing purposes (EEO big option)
      implicit none
      integer idim1
      real*8 xabs(idim1,idim1),xrel(idim1,idim1)
      integer i,j
      do j=1,idim1
        do i=1,idim1
          if (dabs(xrel(i,j)/xabs(i,j)).lt.1.d-2
     *           .and.dabs(xabs(i,j)).lt.1.d-4) xrel(i,j)=0.d0
        enddo
      enddo
      end
c==============================================================================
      subroutine storeEEOampl(ndiskampl, iTstart, iTstop,   nval,   ncf,
     *                        xmatrix)
c Not used at present, it was used once for testing purposes (EEO big option)
      implicit none
      integer ndiskampl,iTstart,iTstop,nval,ncf
      real*8 xmatrix(*)
c
      integer ncf2,ij,ijtri,i,j,iaddress
      ncf2=ncf*ncf
      ij=0
      ijtri=0
      outer: do i=1,nval
        do j=1,i
        ijtri=ijtri+1
        ij=ij+1
        if (ij.gt.iTstop) exit outer
        if (ij.ge.iTstart) then
          iaddress=(ij-iTstart)*ncf2+1
          call reader(ndiskampl,ncf,ijtri,xmatrix(iaddress))
        endif
        if (i.ne.j) then
          ij=ij+1
          if (ij.gt.iTstop) exit outer
          if (ij.ge.iTstart) then
            iaddress=(ij-iTstart)*ncf2+1
            call reader(ndiskampl,ncf,ijtri,xmatrix(iaddress))
            call trans(xmatrix(iaddress),ncf)
          endif
        endif
        enddo
      enddo outer
      end
c==============================================================================
      subroutine transformator(ij,      islaves, irec,  irecskip,lbin,
     *                         indexx,  ibin1,   ibin4, bin8,    ndisk,
     *                         ndisktr, npairs,  ncf,   xmat,    thresh,
     *                         af,      byt8,    nmo,   vorb,    typ,
     *                         nbf)
c This subroutine is used at the end of an integral generation program. It
c performs two kinds of tasks, if vorb=.true. it reads bins, builds the
c integral matrix and transforms it to MO basis. If vorb=.false. it also
c reads bins, builds integral matrix and then writes the matrix onto
c disk again. This improves efficiency of reading such matrix throughout the
c rest of program because each record is equivalent to one matrix and no 
c matrix build from many small bin records is neccessary.
c EEO is an exception - in this case we do not transform matrices at all, the
c matrices are only built and stored. The transformation is performed on the
c fly when EEO are needed, because each transformed EEO is needed only
c once per ij index, and we need non-transformed EEO many times per ij index,
c even of vorb=.true.
c Parameters:
cI  ij         - the pair number
cI  islaves    - number of slaves
cI  irec       - the record info array, maps ij -> record number or number of
c                records written on given slave (the record info is different
c                in single mode and in parallel mode.
cI  irecskip   - this variable is equal to size of record table it had on each
c                slave. If ve move in the record table by this irecskip
c                value we get record info for the same ij pair on the next 
c                slave. This variable is equal to npairs in case of exchange,
c                coulomb and EEO integrals, and nval*nval in case on three
c                external integrals which don't have ij <-> ji symmetry.
cI  lbin       - length of bin
cIO indexx     - this temporary storage is used for reading. The indices of
c                integrals are stored here.
cIO ibin1      - ditto, for overflow of integrals converted to int
cIO ibin4      - ditto for the integrals converrted to int
cIO bin8       - ditto for integrals which are not converted.
cI  ndisk      - the disk logical unit number with not transformed integrals
cI  ndisktr    - the disk unit number where transformed integrals are written
cI  npairs     - number of ij pairs
cI  ncf        - number of contracted functions
cIO xmat       - the memory location for integral matrix construction
cI  thresh     - integral threshold
cI  af         - .true.:  job is parallel and AF are used as IO device.
c                .false.: job is single processor and local IO is used. Also
c                the irec size and structure is different
cI  byt8       - if true we use double precision bins, otherwise bin
c                contains integrals in int4+int1 format
cI  nmo        - number of occupied HF orbitals
cI  vorb       - if .true. the program uses MO virtual space, and this 
c                procedure performs integral 'second-half' transformation,
c                otherwise the virtual space is AO and copying is performed 
c                here
cI  typ        - the type of integrals: exchange, coulomb, EEO, coulomb TEIO,
c                exchange TEIO (Three-External-Integral-Operator). In fact if 
c                'typ' is specified 'irecskip' is not neccessary because it is 
c                possible to evaluate it here, but 'typ' was added lately for 
c                differentiation of EEO, so I kept 'irecskip' for convenience.
      use memory
      use kinds
      implicit none
      integer ij,islaves,irec(*),irecskip,lbin,ndisk,ndisktr,npairs,ncf
      integer nbf
      character*(*) typ
      integer*2 indexx(2,lbin)
      integer*1 ibin1(lbin)
      integer*4 ibin4(lbin)
      real*8    bin8(lbin)
      integer nmo
      real*8 xmat(ncf,ncf) ! also matrix defined outside. (xmat)
      real*8 thresh
      logical af,byt8,vorb
c
      integer iafindex(4),i,j,irecord,jrecord,istatus,icount
      integer idimen,iresult,mataddr,info,itmp,itid,igid,imap,jmap
      real*8 xint,one,dblmax
      parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
      integer list(3,10),iorb,jorb,iii,mu,nu,nval,iadd,irecij
      logical ex3
c
      ex3=.false.
      if (typ.eq.'c'.or.typ.eq.'x'.or.typ.eq.'e') then
        call get_ij_half(ij,iorb,jorb)
      else if (typ.eq.'tx'.or.typ.eq.'tc') then
        nval=(-1.0d0+sqrt(8.0d0*npairs+1))/2.d0+0.001d0
        call get_ij_full(ij,nval, iorb, jorb)
      else
        nval=(-1.0d0+sqrt(8.0d0*npairs+1))/2.d0+0.001d0
        call get_ij_full(ij,nval, iorb, jorb)
        ex3=.true.
      endif
      call zeroit(xmat,ncf*ncf)
      if (af) then ! change for three-external in parallel!
      if (byt8) then
      do j=1,islaves
        irecij=irec(ij+irecskip*(j-1))
        do ! infinite loop
          if (irecij.eq.0) exit
          call bufpointer(iadd,bl(1),lbin*12+8)
          call fafread(ndisk,bl(iadd),1,lbin*12+8,1,irecij,istatus)
          call unpacker(3,bin8,lbin*8,indexx,lbin*4,irecij,8)
          do icount=1, lbin
            if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) 
     +          cycle
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=bin8(icount)*thresh
c           if (typ.eq.'e'.or.typ.eq.'x') then
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=bin8(icount)*thresh*list(3,iii)
            enddo
c           endif
          enddo
        enddo
      enddo
      else
      STOP 'not working now'
      iafindex(2)=ij
      do j=1,islaves
        iafindex(1)=j
        irecord=irec(ij+irecskip*(j-1))
        do i=1,irecord
          iafindex(3)=i
          iafindex(4)=0
          call fafread(ndisk,ibin4,4,lbin,4,iafindex,istatus)
          iafindex(4)=1
          call fafread(ndisk,ibin1,1,lbin,4,iafindex,istatus)
          iafindex(4)=2
          call fafread(ndisk,indexx,2,2*lbin,4,iafindex,istatus)
          do icount=1, lbin
            if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) 
     +          cycle
            if (ibin1(icount).eq.0) then 
              xint=ibin4(icount)
            else if (ibin1(icount).gt.0) then
              xint=ibin4(icount)
              xint=xint+SIGN(ibin1(icount)*dblmax,xint)
            else if (ibin1(icount).lt.0) then
              xint=ibin4(icount)
              xint=xint*10.0d0**(-ibin1(icount))
            endif
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=xint*thresh
c           if (typ.eq.'e'.or.typ.eq.'x') then
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=xint*thresh*list(3,iii)
            enddo
c           endif
          enddo
        enddo
      enddo
      endif
      else   ! if (af) then
      irecord=irec(ij)
      if (byt8) then
      do
        read(ndisk,rec=irecord) bin8,indexx,jrecord
        do icount=1, lbin
          if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) then
              cycle
           endif
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=bin8(icount)*thresh
c           if (typ.eq.'e'.or.typ.eq.'x') then
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=bin8(icount)*thresh*list(3,iii)
            enddo
c           endif
        enddo
        if (jrecord.eq.0) then
          exit
        else
         irecord=jrecord
        endif
      enddo
      else
      do
        read(ndisk,rec=irecord) ibin4,ibin1,indexx,jrecord
        do icount=1, lbin
          if (indexx(1,icount).eq.0 .or. indexx(2,icount).eq.0) 
     +        cycle
          if (ibin1(icount).eq.0) then 
            xint=ibin4(icount)
          else if (ibin1(icount).gt.0) then
            xint=ibin4(icount)
            xint=xint+SIGN(ibin1(icount)*dblmax,xint)
          else if (ibin1(icount).lt.0) then
            xint=ibin4(icount)
            xint=xint*10.0d0**(-ibin1(icount))
          endif
            mu=indexx(1,icount)
            nu=indexx(2,icount)
            xmat(mu,nu)=xint*thresh
c           if (typ.eq.'e'.or.typ.eq.'x') then
            call images_list(iorb,jorb,mu,nu,list,ex3)
            do iii=1,10
              if (list(3,iii).eq.0) exit
              mu=list(1,iii)
              nu=list(2,iii)
              xmat(mu,nu)=xint*thresh*list(3,iii)
            enddo
c           endif
        enddo
        if (jrecord.eq.0) then
          exit
        else
         irecord=jrecord
        endif
      enddo
      endif
      endif ! if (af) then
c Now we have in xint the proper integrals, transform them:
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
c
      call dynamic_matdef('result','q',idimen,idimen)
      iresult=mataddr('result')
c
      if (typ.ne.'e') then
        if (vorb) then
          call matsimtr('xmat','genvirt','result')
        else
          call tfer(xmat,bl(iresult),idimen*idimen)
        endif
      endif
c
      if (typ.ne.'e') then
        if (af) then
          call fafwrite(ndisktr,bl(iresult),8,idimen*idimen,1,ij,info)
        else
          call writer(ndisktr,idimen,ij,bl(iresult))
        endif
      else
        if (af) then
          call fafwrite(ndisktr,xmat,8,ncf*ncf,1,ij,info)
        else
          call writer(ndisktr,ncf,ij,xmat)
        endif
      endif
      call dynamic_matrem('result')
      end
C===============================================================================
      subroutine transformatgr(ij,     islaves,irec, irecskip,lbin,
     *                         indexx, ibin1,  ibin4,bin8,    ndisk,
     *                         ndisktr,npairs, ncf,  xmat,    thresh,
     *                         af,     byt8,   nmo,  vorb,    typ,
     *                         igranules,igranulesize,igran2pair,  ! ->
     *                                                  ipair2gran,list,
     *                         nbf)
c This subroutine differs slightly from the previous transformator only in
c one detail: it uses files where instead of single bin written to one record 
c bunch of bins is stored in one record.
c Additional parameters:
cI  igranules     - how many granules (groups of bins) do we have?
cI  igranulesize  - the size of single granule in ij units
cI  igran2pair    - the array which, for each granule, stores two values:
c                   the first ij pair constituting this granule and the
c                   last ij pair which belongs to this granule.
cI  ipair2gran    - the reverse mapping, this array returns the granule number 
c                   given pair belongs to
      use memory
      use kinds
      implicit none
      integer ij,islaves,irec(*),irecskip,lbin,ndisk,ndisktr,npairs,ncf
      integer nbf
      integer igranules,igranulesize
      character*(*) typ
      integer*2 indexx(2,lbin,igranulesize)
      integer*4 igran2pair(2,*),ipair2gran(*)
      integer list(*)
      integer*1 ibin1(lbin,igranulesize)
      integer   ibin4(lbin,igranulesize)
      real*8    bin8(lbin,igranulesize)
      integer nmo
      real*8 xmat(ncf,ncf) ! also matrix defined outside. (xmat)
      real*8 thresh
      logical af,byt8,vorb
c
      integer iafindex(4),i,j,irecord,jrecord,istatus,icount,ixma
      integer idimen,iresult,mataddr,info,itmp,itid,igid,imap,jmap
      real*8 xint,one,dblmax
      parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
      integer kk,ikk,istart,istop,isize,ishift,iii
      logical eee,l_print,ex3
      integer list1(3,10),iorb,jorb,mu,nu,nval,iadd,irecij
c
      istatus=0
      eee=.false.
      if (typ.eq.'e') eee=.true.
      ex3=.false.
      if (typ.eq.'tt') ex3=.true.
      nval=(-1.0d0+sqrt(8.0d0*npairs+1))/2.d0+0.001d0
      call dynamic_getmem(ncf*ncf*igranulesize,ixma)
      call zeroit(bl(ixma),ncf*ncf*igranulesize)
      istart=igran2pair(1,ij)
      istop =igran2pair(2,ij)
      isize =istop-istart+1
      if (af) then
      iafindex(2)=ij
      if (byt8) then
      do j=1,islaves
        irecij=irec(ij+irecskip*(j-1))
        do
          if (irecij.eq.0) exit
          call bufpointer(iadd,bl(1),lbin*isize*12+8)
         call fafread(ndisk,bl(iadd),1,lbin*isize*12+8,1,irecij,istatus)
          if (istatus.lt.lbin*isize*12+8) then
          ERROR
          print *,istatus,lbin*isize*12+8
          endif
         call unpacker(3,bin8,lbin*isize*8,indexx,lbin*isize*4,irecij,8)
          if (irecij.lt.0) STOP 'Error v6tfkwsmmrmfjrtr'
          do kk=istart,istop
          if (typ.eq.'e') then
            call get_ij_half(kk,iorb,jorb)
          else
            call get_ij_full(kk,nval, iorb, jorb)
          endif
          if (eee) then 
            if (list(kk).ne.1) cycle
          endif
          ishift=(kk-istart)*ncf*ncf
          ikk=kk-istart+1
          do icount=1, lbin
            if (indexx(1,icount,ikk).eq.0.or.indexx(2,icount,ikk).eq.0) 
     +          cycle
            mu=indexx(1,icount,ikk)
            nu=indexx(2,icount,ikk)
            bl(ixma+(mu-1)+(nu-1)*ncf+ishift)=bin8(icount,ikk)*thresh
c           if (typ.eq.'e') then
            call images_list(iorb,jorb,mu,nu,list1,ex3)
            do iii=1,10
              if (list1(3,iii).eq.0) exit
              mu=list1(1,iii)
              nu=list1(2,iii)
              bl(ixma+(mu-1)+(nu-1)*ncf+ishift)=bin8(icount,ikk)*thresh
     *              *list1(3,iii)
            enddo
c           endif
          enddo ! icount=1, lbin
          enddo ! kk=istart,istop
        enddo ! i=1,irecord
      enddo ! j=1,islaves
      else
      iafindex(2)=ij
      do j=1,islaves
        iafindex(1)=j
        irecord=irec(ij+irecskip*(j-1))
        do i=1,irecord
          iafindex(3)=i
          iafindex(4)=0
          call fafread(ndisk,ibin4,4,lbin*isize,4,iafindex,istatus)
          iafindex(4)=1
          call fafread(ndisk,ibin1,1,lbin*isize,4,iafindex,istatus)
          iafindex(4)=2
          call fafread(ndisk,indexx,2,2*lbin*isize,4,iafindex,istatus)
          do kk=istart,istop
          if (typ.eq.'e') then
            call get_ij_half(kk,iorb,jorb)
          else
            call get_ij_full(kk,nval, iorb, jorb)
          endif
          if (eee) then 
            if (list(kk).ne.1) cycle
          endif
          ishift=(kk-istart)*ncf*ncf
          ikk=kk-istart+1
          do icount=1, lbin
            if (indexx(1,icount,ikk).eq.0.or.indexx(2,icount,ikk).eq.0) 
     +          cycle
            if (ibin1(icount,ikk).eq.0) then 
              xint=ibin4(icount,ikk)
            else if (ibin1(icount,ikk).gt.0) then
              xint=ibin4(icount,ikk)
              xint=xint+SIGN(ibin1(icount,ikk)*dblmax,xint)
            else if (ibin1(icount,ikk).lt.0) then
              xint=ibin4(icount,ikk)
              xint=xint*10.0d0**(-ibin1(icount,ikk))
            endif
            mu=indexx(1,icount,ikk)
            nu=indexx(2,icount,ikk)
            bl(ixma+(mu-1)+(nu-1)*ncf+ishift)=xint*thresh
c           if (typ.eq.'e') then
            call images_list(iorb,jorb,mu,nu,list1,ex3)
            do iii=1,10
              if (list1(3,iii).eq.0) exit
              mu=list1(1,iii)
              nu=list1(2,iii)
              bl(ixma+(mu-1)+(nu-1)*ncf+ishift)=xint*thresh*list1(3,iii)
            enddo
c           endif
          enddo ! icount=1, lbin
          enddo ! kk=istart,istop
        enddo ! i=1,irecord
      enddo ! j=1,islaves
      endif
      else   ! if (af) then
      STOP 'Error! Error! Error! Error!'
      endif ! if (af) then
c Now we have in xint the proper integrals, transform them:
      if (vorb) then
      idimen=nbf-nmo
      else
      idimen=ncf
      endif
c
      call dynamic_matdef('result','q',idimen,idimen)
      iresult=mataddr('result')
c
      do kk=istart,istop
      if (eee) then 
        if (list(kk).ne.1) cycle
      endif
      ishift=(kk-istart)*ncf*ncf
      ikk=kk-istart+1
      call matconn('xma','q',ncf,ncf,ixma+ishift)
      if (typ.ne.'e') then
        if (vorb) then
          call matsimtr('xma','genvirt','result')
        else
          call tfer(bl(ixma+ishift),bl(iresult),idimen*idimen)
        endif
      endif
c
      if (typ.ne.'e') then
        if (af) then
          call fafwrite(ndisktr,bl(iresult),8,idimen*idimen,1,kk,info)
        else
          call writer(ndisktr,idimen,kk,bl(iresult))
        endif
      else
        if (af) then
          call fafwrite(ndisktr,bl(ixma+ishift),8,ncf*ncf,1,kk,info)
        else
          call writer(ndisktr,ncf,kk,bl(ixma+ishift))
        endif
      endif
      call matdisc('xma')
      enddo ! do kk=istart,istop
      call dynamic_matrem('result')
      call dynamic_retmem(1)
      end
c===============================================================================
      subroutine Qparts(istart, istop,  jstart, jstop,   ncf,
     *                  nval,   npairs, byt8,   cc,      ccsd,
     *                  list,   nmo,    vorb,   ndiskYZ, af,
     *                  ndiskQ, nrep,   otable, vtable, mtable,
     *                  nbf)
c This subroutine calculates part of Q expression for CCSD
c                    (see: J. Chem. Phys. 104, 6286, (1996), G in this article)
c The calculation is performed in blocks, the block borders are determined by
c istart,istop,jstart,jstop variables. The result is written on disk. The 
c result is only part of full Q formula, which need to be constructed from this
c 'substrats' and other singles contributions.
c
c Parameters:
cI  istart   - the first i for which Q is calculated
cI  istop    - the last i
cI  jstart   - the first j for which Q is calculated
cI  istop    - the last j
cI  ncf      - number of contracted functions
cI  nval     - number of corralated orbitals
cI  npairs   - the number of correlated pairs
cI  byt8     - if .true. we use double presision bins, otherwise int4+int1
c              bins
cI  cc       - if .true. the coupled cluster function is calculated, otherwise
c              we calculate CID, CISD or CEPA-0, or CEPA-2
cI  ccsd     - if .true. the cc=.true. and CCSD singles are included.
cI  list     - the list of strong, weak & distant pairs, for future use.
cI  nmo      - number of occupied HF orbitals
cI  vorb     - if .true. the program uses MO virtual space
c              otherwise the virtual space is AO
cI  ndiskYZ  - the logical disk unit number where Y and Z matrices are kept
c              see: J. Chem. Phys. 104, 6286, (1996)
cI  af       - .true.:  job is parallel and AF are used as IO device.
c              .false.: job is single processor and local IO is used.
cI  ndiskQ   - the logical disk unit number where results are written
c
      use memory
      implicit none
      integer istart,istop,jstart,jstop,nbf
      integer ncf,nval,npairs,list(*),nmo,ndiskYZ,ndiskQ
      integer nrep,otable(*),vtable(*),mtable(8,8)
      logical byt8,cc,ccsd,vorb,af
c
      integer idimen,icoeffi,isize,jsize,i12memory,ik_mem
      integer kj_mem,i12tot,ik_base,kj_base,ik_addr,kj_addr,i1,i2
      integer irec,irecpart1,irecpart2,k,istatus,i,j,idimen2,ij
      logical reset
      integer iprim,jprim,ijprim,ns,nvirt,ikrep,jkrep,irep1,irep2
      logical signum
c
      nvirt=nbf-nmo
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      idimen2=idimen*idimen
c
      reset=.false.
c
      call mmark
      call dynamic_mmark
      call matmark
c
      isize=istop-istart+1
      jsize=jstop-jstart+1
      i12memory=isize*jsize*(idimen*idimen)
      ik_mem=isize*(idimen*idimen)
      kj_mem=jsize*(idimen*idimen)
c     call dynamic_getmem(i12memory,i12tot)
      call dynamic_getmem(i12memory,i12tot)
      call zeroit(bl(i12tot),i12memory)
c     call dynamic_getmem(ik_mem,ik_base)
      call dynamic_getmem(ik_mem,ik_base)
c     call dynamic_getmem(kj_mem,kj_base)
      call dynamic_getmem(kj_mem,kj_base)
c
      do k=1,nval
        do i=istart,istop
          ik_addr=ik_base+(i-istart)*idimen2
          call CoefRead('tt',i,k,icoeffi)
          call CCLoperator(bl(icoeffi),idimen,bl(ik_addr))
        enddo
c
        do j=jstart,jstop
          kj_addr=kj_base+(j-jstart)*idimen2
          call CCYZ_read(k,j,ncf,nval,kj_addr,nmo,vorb,ndiskYZ,af,
     *                   reset,'Y',nbf)
        enddo
        if (nrep.gt.1) irep2=otable(k)
        do i=istart,istop
          ik_addr=ik_base+(i-istart)*idimen2
          call matconn('left','q',idimen,idimen,ik_addr)
          if (nrep.gt.1) then
            irep1=otable(i)
            ikrep=mtable(irep1,irep2)
            if (ikrep.lt.1.or.ikrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
          endif
          do j=jstart,jstop
            if (nrep.gt.1) then
              irep1=otable(j)
              jkrep=mtable(irep1,irep2)
              if (jkrep.lt.1.or.jkrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
            endif
            call array_files(25)
            call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
            call ijtriangle(i,j,ij)
            if (list(ij).ne.1) cycle
            kj_addr=kj_base+(j-jstart)*idimen2
            call matconn('right','q',idimen,idimen,kj_addr)
            i1=i12tot+((j-jstart)+(i-istart)*jsize)*idimen2
            call matconn('1matrix','q',idimen,idimen,i1)
c     subroutine canonical_symm_mult(xmat1, xmat2,  xmat3, id1, id2,
c    *                               id3,   ir1,    ir2,   nvirt,nrep,
c    *                               vtable,m_table,op)
          call canonical_symm_mult(bl(ik_addr),bl(kj_addr),bl(i1),
     *                             idimen, idimen,idimen,ikrep,jkrep,
     *                             nvirt,nrep,vtable, mtable, vorb, 'a')
c           call matmmul2('left','right','1matrix','n','n','a')
            call matdisc('1matrix')
            call matdisc('right')
          enddo
          call matdisc('left')
        enddo
      enddo
c
      do i=istart,istop
        do j=jstart,jstop
            call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
          irec=((i-1)*nval+j)*2
          irecpart1=irec-1
          i1=i12tot+((j-jstart)+(i-istart)*jsize)*idimen2
          if (af) then
            call fafwrite(ndiskQ,bl(i1),8,idimen*idimen,1,
     *                    irecpart1,istatus)
          else
            call writer(ndiskQ,idimen,irecpart1,bl(i1))
          endif
        enddo
      enddo
c
      call zeroit(bl(i12tot),i12memory)
c
      do k=1,nval
        do i=istart,istop
          ik_addr=ik_base+(i-istart)*idimen2
          call CoefRead('tt',k,i,icoeffi)
          call tfer (bl(icoeffi),bl(ik_addr),idimen*idimen)
        enddo
c
        do j=jstart,jstop
          kj_addr=kj_base+(j-jstart)*idimen2
          call CCYZ_read(k,j,ncf,nval,kj_addr,nmo,vorb,ndiskYZ,af,
     *                   reset,'Z',nbf)
        enddo
        if (nrep.gt.1) irep2=otable(k)
        do i=istart,istop
          if (nrep.gt.1) then
            irep1=otable(i)
            ikrep=mtable(irep1,irep2)
            if (ikrep.lt.1.or.ikrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
          endif
          ik_addr=ik_base+(i-istart)*idimen2
          call matconn('left','q',idimen,idimen,ik_addr)
          do j=jstart,jstop
            if (nrep.gt.1) then
              irep1=otable(j)
              jkrep=mtable(irep1,irep2)
              if (jkrep.lt.1.or.jkrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
            endif
            call array_files(26)
            call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
            call ijtriangle(i,j,ij)
            if (list(ij).ne.1) cycle
            kj_addr=kj_base+(j-jstart)*idimen2
            call matconn('right','q',idimen,idimen,kj_addr)
            i2=i12tot+((j-jstart)+(i-istart)*jsize)*idimen2
            call matconn('2matrix','q',idimen,idimen,i2)
          call canonical_symm_mult(bl(ik_addr),bl(kj_addr),bl(i2),
     *                             idimen, idimen,idimen,ikrep,jkrep,
     *                             nvirt,nrep,vtable, mtable, vorb, 'a')
c           call matmmul2('left','right','2matrix','n','n','a')
            call matdisc('2matrix')
            call matdisc('right')
          enddo
          call matdisc('left')
        enddo
      enddo
c
      do i=istart,istop
        do j=jstart,jstop
            call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
          irec=((i-1)*nval+j)*2
          irecpart2=irec
          i2=i12tot+((j-jstart)+(i-istart)*jsize)*idimen2
          if (af) then
            call fafwrite(ndiskQ,bl(i2),8,idimen*idimen,1,
     *                    irecpart2,istatus)
          else
            call writer(ndiskQ,idimen,irecpart2,bl(i2))
          endif
        enddo
      enddo
c
      call dynamic_retmem(3)
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c===============================================================================
      subroutine prepare_Qparts_file(ncf,nmo,vorb,ndiskQ,af,nbf)
c Opens file for Qparts subroutine use.
c Parameters:
cI  ncf      - number of contracted functions
cI  nmo      - number of occupied HF orbitals
cI  vorb     - if .true. the program uses MO virtual space
c              otherwise the virtual space is AO
cO  ndiskQ   - the logical disk unit number
cI  af       - .true.:  job is parallel and AF are used as IO device.
c              .false.: job is single processor and local IO is used.
c
      use ccounters
      implicit none
      integer ncf,nmo,ndiskQ,nbf
      logical vorb,af
c
      character*256 filename,scrfile
      integer len,len1,lenrec,idimen,irecord
      integer info
      logical fileopen
      save filename,len1
c
      ic_prepare_Qparts_file=ic_prepare_Qparts_file+1
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      if (.not.af) ndiskQ=72
      irecord=idimen*idimen*8
      if (ic_prepare_Qparts_file.eq.1) then
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filename=scrfile(1:len)//'.Qparts'
      len1=len+7
      endif
      if (af) then
        if (ic_prepare_Qparts_file.gt.1) then
          call fafClosem(ndiskQ,0,info)
        else
          call fafClosem(ndiskQ,-1,info)
        endif
      else
        if (ic_prepare_Qparts_file.eq.1) then
          inquire(ndiskQ,opened=fileopen)
          if (fileopen) STOP 'Unit 72 was used by somebody'
        else
          CLOSE(ndiskQ,STATUS='DELETE')
        endif
      endif
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndiskQ,irecord)
      else
      OPEN(UNIT=ndiskQ,FILE=filename(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=irecord)
      endif
      end
C==============================================================================
      subroutine Q_read_build(i,    j,   iX,  ndiskQ, af,
     *                        vorb, ncf, nmo, nval,   iqij,
     *                        iqji, nbf)
c This subroutine reads the matrices produced by Qparts and constructs
c most of the Q defined by J. Chem. Phys. 104, 6286, (1996), G in this article
c Parameters:
cI  i       - the i correlared orbital index
cI  j       - the j correlared orbital index
cI  iX      - the X matrix (see the same article)
cI  ndiskQ  - the logical disk unit number where we keep parts of Q expression
cI  af      - .true.:  job is parallel and AF are used as IO device.
c             .false.: job is single processor and local IO is used.
cI  vorb    - if .true. the program uses MO virtual space
c             otherwise the virtual space is AO
cI  ncf     - number of contracted functions
cI  nmo     - number of molecular orbitals
cI  nval    - number of correlated orbitals
cO  iqij    - the first result Qij
cO  iqji    - the second result Qji
c
      use memory
      implicit none
      integer i,j,iX,ndiskQ,ncf,nmo,nval,iqij,iqji,nbf
      logical af,vorb
c
      integer idimen,i1ij,i2ij,i1ji,i2ji,iaddress,mataddr
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('1ij','q',idimen,idimen)
      call dynamic_matdef('2ij','q',idimen,idimen)
      call dynamic_matdef('1ji','q',idimen,idimen)
      call dynamic_matdef('2ji','q',idimen,idimen)
      i1ij=mataddr('1ij')
      i2ij=mataddr('2ij')
      i1ji=mataddr('1ji')
      i2ji=mataddr('2ji')
      call matconn('qij','q',idimen,idimen,iqij)
      call matconn('qji','q',idimen,idimen,iqji)
      call matconn('X','q',idimen,idimen,iX)
c
      call read_Q_parts(i,j,nval,'1',af,ndiskQ,ncf,nmo,vorb,i1ij,nbf)
      call read_Q_parts(i,j,nval,'2',af,ndiskQ,ncf,nmo,vorb,i2ij,nbf)
      call read_Q_parts(j,i,nval,'1',af,ndiskQ,ncf,nmo,vorb,i1ji,nbf)
      call read_Q_parts(j,i,nval,'2',af,ndiskQ,ncf,nmo,vorb,i2ji,nbf)
      call matzero('qij')
      call matadd('1ij','qij')
      call matadd1('2ij',-0.5d0,'qij')
      call matadd1('2ji',-1.0d0,'qij')
      call CoefRead('tt',i,j,iaddress)
      call matconn('coeffQ','q',idimen,idimen,iaddress)
      call matmmul2('coeffQ','X','qij','n','n','a')
      call matdisc('coeffQ')
c
      call matzero('qji')
      call matadd('1ji','qji')
      call matadd1('2ji',-0.5d0,'qji')
      call matadd1('2ij',-1.0d0,'qji')
      call CoefRead('tt',j,i,iaddress)
      call matconn('coeffQ','q',idimen,idimen,iaddress)
      call matmmul2('coeffQ','X','qji','n','n','a')
      call matdisc('coeffQ')
      call matpose('qji')
      call matremark
      call retmark
      call dynamic_retmark
c
      end
C==============================================================================
      subroutine read_Q_parts(i,      j,   nval, typ,  af,
     *                        ndiskQ, ncf, nmo,  vorb, istorage,
     *                        nbf)
c This subroutine reads from disk output of Qparts for given i,j. The output
c is deposited into istorage bl index.
c Parameters:
c
cI  i        - the i correlared orbital index
cI  j        - the j correlared orbital index
cI  nval     - number of correlated orbitals
cI  typ      - this variable tells read_Q_parts which part callers needs.
cI  af       - .true.:  job is parallel and AF are used as IO device.
c              .false.: job is single processor and local IO is used.
cI  ndiskQ   - the logical disk unit number where we keep parts of Q expression
cI  ncf      - number of contracted functions
cI  nmo      - number of molecular orbitals
cI  vorb     - if .true. the program uses MO virtual space
c              otherwise the virtual space is AO
cO  istorage - the result is written to this bl index.
c 
      use memory
      implicit none
c
      integer i,j,nval,ndiskQ,ncf,nmo,istorage,nbf
      logical af,vorb
      character*1 typ
c
c
      integer idimen,irec,irecpart1,irecpart2,irealrec,istatus
      integer iprim,jprim,ijprim,ns,iscr
      logical signum
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
c     irec=((i-1)*nval+j)*2
      call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                        ns, signum)
      irec=ijprim*2
      irecpart1=irec-1
      irecpart2=irec
      if      (typ.eq.'1') then
        irealrec=irecpart1
      else if (typ.eq.'2') then
        irealrec=irecpart2
      endif
      if (af) then
        call fafread(ndiskQ,bl(istorage),8,idimen*idimen,1,
     *                irealrec,istatus)
          if (istatus.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
      else
      call reader(ndiskQ,idimen,irealrec,bl(istorage))
      endif
      if (ns.ne.0) then
        call full_symm_transform(idimen,ns,signum,bl(istorage))
      endif
      end
C==============================================================================
      subroutine print_CCSD_results(ccsdela,ccsdcpu,itim_no)
      implicit none
      integer itim_no
      real*8 ccsdela(itim_no),ccsdcpu(itim_no)
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,*) '                    The timings:     '
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      write(6,'(A,F10.2)') 'Elaps. EEO operators:  ',ccsdela(1)
      write(6,'(A,F10.2)') 'CPU    EEO operators:  ',ccsdcpu(1)
      write(6,'(A,F10.2)') 'Elaps. tx operators:   ',ccsdela(2)
      write(6,'(A,F10.2)') 'CPU    tx operators:   ',ccsdcpu(2)
      write(6,'(A,F10.2)') 'Elaps. tc operators:   ',ccsdela(3)
      write(6,'(A,F10.2)') 'CPU  . tc operators:   ',ccsdcpu(3)
      write(6,'(A,F10.2)') 'Elaps. CCalphaonce:    ',ccsdela(4)
      write(6,'(A,F10.2)') 'CPU    CCalphaonce:    ',ccsdcpu(4)
      write(6,'(A,F10.2)') 'Elaps. Builder41a_G:   ',ccsdela(5)
      write(6,'(A,F10.2)') 'CPU    Builder41a_G:   ',ccsdcpu(5)
      write(6,'(A,F10.2)') 'Elaps. CCAOnce:        ',ccsdela(6)
      write(6,'(A,F10.2)') 'CPU    CCAOnce:        ',ccsdcpu(6)
      write(6,'(A,F10.2)') 'Elaps. buildX:         ',ccsdela(7)
      write(6,'(A,F10.2)') 'CPU    buildX:         ',ccsdcpu(7)
      write(6,'(A,F10.2)') 'Elaps. singles s, r:   ',ccsdela(8)
      write(6,'(A,F10.2)') 'CPU    singles s, r:   ',ccsdcpu(8)
      write(6,'(A,F10.2)') 'Elaps. YZ generation:  ',ccsdela(9)
      write(6,'(A,F10.2)') 'CPU    YZ generation:  ',ccsdcpu(9)
      write(6,'(A,F10.2)') 'Elaps. Q parts:        ',ccsdela(10)
      write(6,'(A,F10.2)') 'CPU    Q parts:        ',ccsdcpu(10)
      write(6,'(A,F10.2)') 'Elaps. Rij resid:      ',ccsdela(11)
      write(6,'(A,F10.2)') 'CPU    Rij resid:      ',ccsdcpu(11)
      write(6,'(A,F10.2)') 'Elaps. the whole iter: ',ccsdela(12)
      write(6,'(A,F10.2)') 'CPU    the whole iter: ',ccsdcpu(12)
      write(6,'(A,F10.2)') 'Elaps. singles resid:  ',ccsdela(13)
      write(6,'(A,F10.2)') 'CPU    singles resid:  ',ccsdcpu(13)
      write(6,'(A,F10.2)') 'Elaps. DIIS:           ',ccsdela(14)
      write(6,'(A,F10.2)') 'CPU    DIIS:           ',ccsdcpu(14)
      write(6,'(A,F10.2)') 'Elaps. Tl_generator:   ',ccsdela(15)
      write(6,'(A,F10.2)') 'CPU    Tl_generator:   ',ccsdcpu(15)
      write(6,'(A,F10.2)') 'Elaps. EEO_vector_ex:  ',ccsdela(16)
      write(6,'(A,F10.2)') 'CPU    EEO_vector_ex:  ',ccsdcpu(16)
      write(6,'(A,F10.2)') 'Elaps. dump checkpoint:',ccsdela(17)
      write(6,'(A,F10.2)') 'CPU    dump checkpoint:',ccsdcpu(17)
      write(6,*) '* * * * * * * * * * * * * * * * * * * * * * * * * * *'
      end
C==============================================================================
      subroutine open_mp3_resid(ndisk,idimen,af,lopen)
c Creates and opens file for MP3 residuals needed for MP4.
c Parameters:
cO  ndisk  - the disk unit number
cO  idimen - tis variable is used for record size calculation
cI  af     - .true.:  job is parallel and AF are used as IO device.
c            .false.: job is single processor and local IO is used.
c
      implicit none
      character*256 scrfile,filname1
      integer len,len1,ndisk,lrec,idimen
      logical af,lopen
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
      filname1=scrfile(1:len)//'.mp3_res'
      len1=len+8
      lrec=idimen*idimen*8
      if (af) then
        if (lopen) then
        call fafCreatem(filname1(1:len1)//char(0),ndisk,lrec)
        else
        call fafCreatem(filname1(1:len1)//char(0),ndisk,-1)
        endif
      else
        if (lopen)
     *  OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL=lrec)
      endif
      end
c==============================================================================
      subroutine do_mp4_beta(fock,beta,nval,ncf,nfirst)
c This subroutine fills 'beta' array with MO fock matrix elements, in fact
c the only task performed here is picking up apprioprate elements of fock
c and copying them into 'beta'.
c Parameters:
cI  fock   - The full MO Fock matrix
cO  beta   - The resulting beta matrix
cI  nval   - number of correlated orbitals
cI  ncf    - number of contracted functions
cI  nfirst - first correlated HF orbital (this with lowest orbital energy 
c            amongst all correlated orbitals).
c
      implicit none
      integer nval,ncf,nfirst
      real*8 fock(ncf,ncf),beta(nval,nval)
      integer i,j,incf,jncf
      call dynamic_unlock(beta,i)
      do j=1,nval
        jncf=nfirst+j-1
        do i=1,nval
          incf=nfirst+i-1
          beta(i,j)=fock(incf,jncf)
        enddo
      enddo
      call dynamic_lock(beta,i)
      end
c==============================================================================
      subroutine storesingles(i,vector,xmatrix,ncf,idimen,nfirst)
c This subroutine, used in MP4 singles contribution calculation stores
c singles 'i' vector in apprioprate column of square matrix 'xmatrix'
cI  i       - the number of singles vector (the number of correlated pair)
cI  vector  - the source singles vector
cO  xmatrix - the destination matrix
cI  ncf     - number of contracted basis functions
cI  idimen  - dimension of vector
cI  nfirst  - first correlated HF orbital (this with lowest orbital energy 
c             amongst all correlated orbitals).
c
      implicit none
      integer i,ncf,idimen,nfirst
      real*8 vector(*),xmatrix(ncf,ncf)
      integer ii,incf,ishift
      incf=nfirst+i-1
      ishift=ncf-idimen
      do ii=1,idimen
        xmatrix(incf,ii+ishift)=vector(ii)
      enddo
      end
c==============================================================================
      function log2int(integ)
      implicit none
      integer integ,log2int
      log2int=abs(integ)
      end
c==============================================================================
      subroutine Generate_MP2_G(nval, ncf, ndisk, nmo,  nfirst,
     *                          vorb, af,  F,     nbf)
c This subroutine will probably fail for some big calculation in the future, 
c because it needs fixed amount of memory (2*idimen*idimen*nval).
c It calculates "G" part in MP2. The "G" part is defined as
c \sum_k \mathbf{T}^{kj}*F_{ik}+\mathbf{T}^{ik}*F_{kj}
c The calculation is performed in blocs and matrix summatrion is replaced
c by matrix multiplication
cI  nval   - the number of correlated orbitals
cI  ncf    - the number of basis functions
cI  ndisk  - the logical disk unit number where the results are written
cI  nmo    - the number of occupied HF orbitals
cI  nfirst - The first correlated oribtal, i.e. the correlated orbital with
c            lowest HF energy
cI  vorb   - if .true. the program uses MO virtual space
c            otherwise the virtual space is AO
cI  af     - .true.:  job is parallel and AF are used as IO device.
c            .false.: job is single processor and local IO is used.
cI  F      - The Fock matrix in MO basis (diagonal in canonical case)
      use memory
      implicit none
      integer nval,ncf,ndisk,nmo,nfirst,nbf
      logical vorb,af
      real*8 F(ncf,ncf)
c local:
c
      integer i,j,k,iresult,iamplit,ka,ir,ij,idimen,ifp,jfp,itmp,ifock
      integer istatus,iprim,jprim,ijprim,ns
      logical signum,jump
      integer*1 iunique(nval,nval)
c
      do i=1,nval
        do j=1,nval
          call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                            ns, signum)
          iunique(i,j)=ns
        enddo
      enddo
      call dynamic_mmark()
      call matmark()
      call dynamic_getmem(nval*nval,ifock)
      call matconn('fock_occ','q',nval,nval,ifock)
      do j=1,nval
        jfp=j+nfirst-1
        do i=1,nval
          ifp=i+nfirst-1
          bl(ifock+(i-1)+(j-1)*nval)=F(ifp,jfp)
        enddo
      enddo
      call matscal('fock_occ',-1.0d0)
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call dynamic_getmem(idimen*idimen*nval,iresult)
      call dynamic_getmem(idimen*idimen*nval,iamplit)
      do j=1,nval   ! result column
        jump=.true.
        do i=1,nval
          if (iunique(i,j).eq.0) jump=.false.
        enddo
        if (jump) cycle
        do k=1,nval
          ka=iamplit+(k-1)*idimen*idimen
          call CoefRead('tt',k,j,itmp)
          call tfer(bl(itmp),bl(ka),idimen*idimen)
        enddo
        call matconn('result','r',idimen*idimen,nval,iresult)
        call matconn('amplit','r',idimen*idimen,nval,iamplit)
        call matmmult('amplit','fock_occ','result')
        call matdisc('amplit')
        call matdisc('result')
        do i=1,nval ! result row
          ir=iresult+(i-1)*idimen*idimen
          ij=(j-1)*nval+i
          if (iunique(i,j).ne.0) cycle
          if (af) then
            call fafwrite(ndisk,bl(ir),8,idimen*idimen,1,ij,istatus)
          else
          call writer(ndisk,idimen,ij,bl(ir))
          endif
        enddo
      enddo
      call matremark()
      call dynamic_retmark()
      end
c==============================================================================
      subroutine Open_Scratch_File(ndisk,idim1,idim2,af,suffix)
c This is a general procedure which should be used for opening all scratch files
c using both AF and local files. The ndisk is output in parallel (af=.true.)
c but input in single mode. I should use ONLY this subroutine in the future for
c all files opening, because I could have control over file numbers in 
c this subroutine in single processor case. Right now files numbering is 
c a huge mess in single mode.
cIO ndisk  - the disk unit number, the value must be provided for local disk
c            opening, but it is returned in if AF is used. I should change this
c            in order to make files numbering unifirm in both cases.
cI  idim1  - the files are indended for direct access, where one record is
c            one matrix. This is the first dimension of this matrix. Also 
c            I should simplify this and directly specify record length
cI  idim2  - ditto for second dimension
cI  af     - .true.:  job is parallel and AF are used as IO device.
c            .false.: job is single processor and local IO is used.
cI  suffix - this character suffix is added to the default file name
c
      implicit none
      integer ndisk,idim1,idim2
      logical af
      character*(*) suffix
      character*256 scrfile,filname1
      integer len0,len1,lrec
      logical fileopen
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len0) !returns len0 of string without spaces,
      filname1=scrfile(1:len0)//'.'//suffix
      len1=len0+1+len(suffix)
      lrec=idim1*idim2*8
      if (af) then
        call fafCreatem(filname1(1:len1)//char(0),ndisk,lrec)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) then 
        print *,fileopen,ndisk
        call nerror(1,'Open_Scratch_File','File already opened',1,1)
        endif
        OPEN (UNIT=ndisk,FILE=filname1(1:len1),FORM='UNFORMATTED',
     1        ACCESS='DIRECT',RECL=lrec)
      endif
      end
c==============================================================================
      subroutine Read_MP2(i,     j,    nmo,   ncf,     nval,
     *                    ndisk, vorb, af,    iresult, nbf)
c This subroutine is used for reading the matrices calculated by Generate_MP2_G
c and building the real MP2 "G" part from them. OPTIMIZE: when I finish adding
c this comments I should remove one superflous memory copying here.
cI  i,j     - the G i,j indices
cI  nmo     - the number of occupied HF orbitals
cI  ncf     - the number of basis functions
cI  nval    - the number of correlated orbitals
cI  ndisk   - the logical disk unit number where the results are written
cI  vorb    - if .true. the program uses MO virtual space
c             otherwise the virtual space is AO
cI  af      - .true.:  job is parallel and AF are used as IO device.
c             .false.: job is single processor and local IO is used.
cO  iresult - the resulting matrix is deposited here
c
c
      use memory
      implicit none
      integer i,j,nmo,ncf,nval,ndisk,iresult,nbf
      logical vorb,af
c local
      integer idimen,ij,ji,iwork11,mataddr,istatus,ns,iprim,jprim,ijprim
      integer iwork22,k,l
      logical signum
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call matconn('ress','q',idimen,idimen,iresult)
      call dynamic_matdef('work11','q',idimen,idimen)
      iwork11=mataddr('work11')
      ij=(j-1)*nval+i
      ji=(i-1)*nval+j
      call full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                        ns, signum)
      ijprim=(jprim-1)*nval+iprim
      if (af) then
        call fafread(ndisk,bl(iwork11),8,idimen*idimen,1,ijprim,istatus)
          if (istatus.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
      else
        call reader(ndisk,idimen,ijprim,bl(iwork11))
      endif
      call full_symm_trans_cop(idimen,ns,signum,bl(iwork11),bl(iresult))
      call full_pair_searcher(j,  i,     jprim, iprim, ijprim,
     *                        ns, signum)
      ijprim=(iprim-1)*nval+jprim
      if (af) then
        call fafread(ndisk,bl(iwork11),8,idimen*idimen,1,ijprim,istatus)
          if (istatus.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
      else
        call reader(ndisk,idimen,ijprim,bl(iwork11))
      endif
      call full_symm_tr_add_tr(idimen,ns,signum,bl(iwork11),bl(iresult))
      call dynamic_matrem('work11')
      call matdisc('ress')
      end
c============================================================================
      subroutine general_read(ndisk,irec,xmat,ibyte_dim)
      implicit none
      integer ndisk,irec,ibyte_dim
      integer*1 xmat(ibyte_dim)
      read(ndisk,REC=irec) xmat
      end
c============================================================================
      subroutine general_write(ndisk,irec,xmat,ibyte_dim)
      implicit none
      integer ndisk,irec,ibyte_dim
      integer*1 xmat(ibyte_dim)
      write(ndisk,REC=irec) xmat
      end
c============================================================================
      subroutine general_read1(ndisk,xmat,ibyte_dim)
      implicit none
      integer ndisk,ibyte_dim
      integer*1 xmat(ibyte_dim)
      read(ndisk) xmat
      end
c============================================================================
      subroutine general_write1(ndisk,xmat,ibyte_dim)
      implicit none
      integer ndisk,ibyte_dim
      integer*1 xmat(ibyte_dim)
      write(ndisk) xmat
      end
c============================================================================
      subroutine i2r(intmat,realmat,idimen,xmax)
c the opposite of r2i described below.
      implicit none
      integer idimen
      integer*4 intmat(idimen)
      real*8 realmat(idimen),xmax
c
      integer i
      real*8 xmult,xmultinv
c
      xmult=2147483648.d0/xmax
      xmultinv=1.d0/xmult
c
      do i=1,idimen
      realmat(i)=dble(intmat(i))*xmultinv
      enddo
      end
c============================================================================
      subroutine r2i(intmat,realmat,idimen,xmax)
c This subroutine converts the double precision array into integer array,
c with the assumption that the largest element of real*8 array is 'xmax'
cO  intmat  - the resulting integer matrix
cI  realmat - the source double precision matrix
cI  idimen  - total length of the matrix being converted, because conversion
c             does not need know anything about shape of the matrix, it is
c             sufficient to assume that matrix is linear.
cI  xmax    - the assumed maximum element in source matrix
      implicit none
      integer idimen
      integer*4 intmat(idimen)
      real*8 realmat(idimen),xmax
c
      integer i
      real*8 xmult,xmultinv
c
      xmult=2147483648.d0/xmax
      xmultinv=1.d0/xmult
c
      do i=1,idimen
      intmat(i)=nint(realmat(i)*xmult)
      enddo
      end
c============================================================================
      subroutine BinSortGr(nval,   lbin,  npairs,ndisk,   halftra,
     1                     mu,     lam,   thresh,icounter,irecord,
     2                     indxbin,ibin4, ibin1, bin8,    afflag,  
     3                     afname, islvid,byt8,  igranules,igranulesize,
     4                     igran2pair,ipair2gran)
c See BinSort_tj description. This subroutine differs from the BinSort_tj in
c few details. It writes a group of bins at once instead of one bin as 
c traditional BinSort does. The groups info is the same through program
c and described here by additional parameters descibed below. For other
c parameters see BinSort_tj.
c Parameters:
cI  igranules     - how many granules (groups of bins) do we have?
cI  igranulesize  - the size of single granule in ij units
cI  igran2pair    - the array which, for each granule, stores two values:
c                   the first ij pair constituting this granule and the
c                   last ij pair which belongs to this granule.
cI  ipair2gran    - the reverse mapping, this array returns the granule number 
c                   given pair belongs to
c
        implicit real*8 (a-h,o-z)
        real*8 halftra(nval,nval)
        integer*2 indxbin(2,lbin,npairs)
        integer*4 ibin4(lbin,npairs)
        integer*1 ibin1(lbin,npairs)
        real*8    bin8(lbin,npairs)
        integer irecord(npairs),icounter(npairs)
        logical afflag,byt8
        parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
        parameter (dblinv=one/dblmax,d1max=128.0d0)
        character*(*) afname
        integer*4 igran2pair(2,igranules),ipair2gran(npairs)
        ij=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            abhlftra=abs(halftra(i,j))
            if(abhlftra.gt.thresh) then
              if (icounter(ij).eq.lbin) then
                igranule=ipair2gran(ij)
                istart=igran2pair(1,igranule)
                istop =igran2pair(2,igranule)
                isize=istop-istart+1
                if(afflag) then 
                  call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),
     1                            ibin1(1,istart),bin8(1,istart),
     2                            indxbin(1,1,istart),irecord(igranule),
     3                            igranule,islvid,byt8,isize)
                else ! FRAGILE not finished for single version!
                  ERROR
                  STOP 'ERROR BinSortGr'
                  call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1              bin8(1,ij), indxbin(1,1,ij),irecord(ij),ij,npairs,
     2              byt8)
                endif
                do kk=istart,istop
                  icounter(kk)=0
                enddo
              end if
              ii=icounter(ij)+1
              icounter(ij)=ii
              indxbin(1,ii,ij)=mu
              indxbin(2,ii,ij)=lam
              x=halftra(i,j)
              byte8: if (byt8) then
              bin8(ii,ij)=halftra(i,j)
              else
              IF(abs(x).ge.dblmax) THEN
                b = x*dblinv
                if(abs(b).ge.d1max) then
                  dfac = abs(x)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
                  dfac = LOG10(dfac)
                  i1 = -NINT(dfac+0.5d0)
                  x = x*10.0d0**i1
                  ibin1(ii,ij)=i1
                  ibin4(ii,ij)=x
                else
                  i1 = abs(b)
                  b = x - SIGN(i1*dblmax,x)
                  ibin1(ii,ij)  = i1
                  ibin4(ii,ij) = b
                endif
              ELSE
                ibin1(ii,ij) = 0
                ibin4(ii,ij) = x
              ENDIF
              endif byte8
            end if
          end do
        end do
        if (mu.ne.lam) then
          ij=0
          do i=1,nval
            do j=1,i
              ij=ij+1
              abhlftra=abs(halftra(j,i))
              if(abhlftra.gt.thresh) then
                if (icounter(ij).eq.lbin) then
                  igranule=ipair2gran(ij)
                  istart=igran2pair(1,igranule)
                  istop =igran2pair(2,igranule)
                  isize=istop-istart+1
                  if(afflag) then 
                  call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),
     1                            ibin1(1,istart),bin8(1,istart),
     2                            indxbin(1,1,istart),irecord(igranule),
     3                            igranule,islvid,byt8,isize)
                  else ! FRAGGILE not finished for single version!
                    call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1                bin8(1,ij),indxbin(1,1,ij),irecord(ij), ij,npairs,
     2                byt8)
                  endif
                  do kk=istart,istop
                    icounter(kk)=0
                  enddo
                endif
                ii=icounter(ij)+1
                icounter(ij)=ii
                indxbin(1,ii,ij)=lam
                indxbin(2,ii,ij)=mu
                x=halftra(j,i)
                byte81: if (byt8) then
                bin8(ii,ij)=halftra(j,i)
                else
                IF(abs(x).ge.dblmax) THEN
                  b = x*dblinv
                  if(abs(b).ge.d1max) then
                    dfac = abs(x)/dblmax
c -- we are going to reduce the threshold by powers of ten for this
c    integral until it can be stored in 4-bytes
                    dfac = LOG10(dfac)
                    i1 = -NINT(dfac+0.5d0)
                    x = x*10.0d0**i1
                    ibin1(ii,ij)=i1
                    ibin4(ii,ij)=x
                  else
                    i1 = abs(b)
                    b = x - SIGN(i1*dblmax,x)
                    ibin1(ii,ij)  = i1
                    ibin4(ii,ij) = b
                  endif
                ELSE
                  ibin1(ii,ij) = 0
                  ibin4(ii,ij) = x
                ENDIF
                endif byte81
              end if
            end do
          end do
        end if
      end
c==============================================================================
      subroutine check_size2(lmp2_size, ncf,    nval,   merg,   npairs,
     *                       typ,       ntimes, nodisk, split)
c This subroutine measures how much memory do we have and whether
c splitting of integrals calculation into passes will be necessary.
c Parameters:
cI  lmp2_size - The calculated size of resulting matrix for the whole shell
cI  ncf       - the number of contracted basis functions
cI  nval      - the number of correlated orbitals
cI  merg      - not used
cI  npairs    - number of ij pairs = nval*(nval+1)/2
cI  typ       - the kind of integrals: exchange 'x', coulomb 'c', EEO 'e',
c               coulomb TEIO 'tc', exchange TEIO 'tx'.
cO  ntimes    - the division factor, how many parts do we need to divide
c               lmp2_size into?
cI  nodisk    - this logical variable bears information about transformed
c               amplitudes location: if .true. all transformed for EEO
c               amplitudes are kept in memory, otherwise on disk
cI  split     - if true merged shells are splited, otherwise not.
c
      character*(*) typ
      logical nodisk,split
c
      call dynamic_show_free(memory)
c
c
cBEGIN I have to think whether this still is needed:
      memory=memory - 1500000  !  leave some mem for integrals :
cEND
      if (typ.eq.'e'.and..not.nodisk.and.split) then
        memory=memory-30*ncf*ncf    ! amplitudes in EEO fragile!
      else if (typ.eq.'e'.and..not.nodisk.and..not.split) then
        memory=memory-ncf*ncf/2     ! for screening irow(nrow)
      endif
c
c
      if(lmp2_size.gt.memory) then
         ntimes=lmp2_size/memory 
         if(ntimes*memory.lt.lmp2_size) ntimes=ntimes+1
         if (ntimes.lt.1) then
         print *,'check_size2 error, ntimes<1: ',ntimes,memory,lmp2_size
         STOP    'check_size2 error, ntimes<1  '
         endif
         write(*,*) 'Splitting shell... ',ntimes,' times.'
         write(*,*) 'Memory was, needed: ',memory,lmp2_size
      else
         ntimes=1
      endif
c
      end
c==============================================================================
      subroutine i1zeroit(matrix,length)
      implicit none
      integer*1 matrix(*)
      integer length,i
      do i=1,length
        matrix(i)=0
      enddo
      end
c==============================================================================
      subroutine i2zeroit(matrix,length)
      implicit none
      integer*2 matrix(*)
      integer length,i
      do i=1,length
        matrix(i)=0
      enddo
      end
c==============================================================================
      subroutine i4zeroit(matrix,length)
      implicit none
      integer*4 matrix(*)
      integer length,i
      do i=1,length
        matrix(i)=0
      enddo
      end
c==============================================================================
      subroutine build_idics(screen,shellscreen,ncf,ncs,inx)
c Assuming that screening matrix consists of screening elements for each
c contracted function icf,kcf for all jcf and lcf in (ijkl) integral, this
c subroutine costructs screening matrix for ics and kcs shells by picking up
c maximum element in screen for all icf and kcf which belong to given ics and
c kcs pair
c Parameters:
cI  screen      - the source ncf*ncf screen matrix
cO  shellscreen - the resulting ncs*ncs screen matrix
cI  ncf         - the number of contracted functions
cI  ncs         - the number of contracted shells
cI  inx         - the contraction info
c
      implicit none
      integer ncf,ncs
      integer inx(12,*)
      real*8 screen(ncf,ncf),shellscreen(ncs,ncs)
c
      real*8 xmax
      integer ics,jcs,icfb,icfe,jcfb,jcfe,icf,jcf
c
      do ics=1,ncs
        icfb=inx(11,ics)+1
        icfe=inx(10,ics)
        do jcs=1,ncs
          jcfb=inx(11,jcs)+1
          jcfe=inx(10,jcs)
          xmax=0.0d0
          do icf=icfb,icfe
            do jcf=jcfb,jcfe
            if (screen(icf,jcf).gt.xmax) xmax=screen(icf,jcf)
            enddo
          enddo
          shellscreen(ics,jcs)=xmax
        enddo
      enddo
      end
c
c======================================================================
      subroutine alloc_results(xsource,isource1,isource2,
     *                         xtarget,itarget1,itarget2,
     *                         irowshift)
c This subroutine is used by EEO with screening. It copies results located 
c in smaller matrix xsource into matrix xtarget. Because those matrices differ
c in the row dimension, the copying cannot be simple memory move.
c Parameters:
cI  xsource   - the source matrix
cI  isource1  - the first dimension of xsource
cI  isource2  - ditto the secnd dimension
cO  xtarget   - the target (larger) matrix
cI  itarget1  - the first dimension of xtarget
cI  itarget2  - ditto the secnd dimension
cI  irowshift - this variable indicates the shift in bigger matrix in row
c               index. All values from smaller matrix will be put into
c               rows starting from irowshift.
c
      implicit none
      integer isource1,isource2,itarget1,itarget2,irowshift
      real*8 xsource(isource1,isource2),xtarget(itarget1,itarget2)
c
      integer i,j,itarg,imyshift
      imyshift=irowshift-1
      do j=1,isource2
        do i=1,isource1
          itarg=i+imyshift
          xtarget(itarg,j)=xsource(i,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine copy_compression_info(itarget,irow,nrow)
c Simple subroutine which allows to write integers into real*8 array.
      implicit none
      integer*4 itarget(*),irow(*)
      integer nrow
      integer i
      do i=1,nrow
      itarget(i)=irow(i)
      enddo
      end
c==============================================================================
      subroutine C_compress(source,   idims1,    idims2,
     *                      dest,     idimd1,    idimd2,
     *                      irow)
c This simple subroutine compresses source matrix having dimensions
c idims1,idims2 into dest matrix with targest dimensions idimd1,idimd2
c according to compression info in irow integer array, which maps compressed
c matrix indices into original (old)m atrix indices. Compression affects rows
c only.
      implicit none
      integer idims1,idims2,idimd1,idimd2
      integer*4 irow(*)
      real*8 source(idims1,idims2),dest(idimd1,idimd2)
c
      integer i,j,inewcolshift,ioldcolshift,iold
c
      if (idims1.eq.idimd1) return
      if (idims2.ne.idimd2) STOP 'Error in C_compress'
      do j=1,idims2
        do i=1,idimd1
          iold=irow(i)
          dest(i,j)=source(iold,j)
        enddo
      enddo
      end
c==============================================================================
      subroutine what_to_do(icur_size,isize,ik,ik_prev,ihow_many,lstop)
c This subroutine returns logical variable lstop. This variable, when set to
c .true. informs caller that it should cut the bunch of integrals here for 
c separate strip. The strip cannot consist of more that 'max_how_many' shells
c with size bigger than min_shell, it cannot be bigger than max_merg etc. lot 
c of conditions is here and I hope they are self-explaining.
cI  icur_size - current shell size (in fact this is ics*kcs size)
cI  isize     - next shell size (this shell is not transformed yet)
cI  ik        - the starting ik of present shell which is considered to be
c               appended to old strip or not.
cI  ik_prev   - the starting ik of current strip, not closed yet.
cI  ihow_many - how many bigger shell can be merged
cO  lstop     - if .true. it means that we need to cut here.
c
      implicit none
      integer icur_size,isize,ik,ik_prev,ihow_many
      logical lstop
      integer min_shell,max_merg,max_how_many
ctmp  parameter (min_shell=15,max_merg=20)
      parameter (min_shell=5,max_merg=12,max_how_many=1)
c store up to max_merg functions from shell smaller than min_shell
      if (icur_size.lt.min_shell.and.isize.lt.min_shell.and.
     *        (ik-ik_prev).ge.max_merg) then
         lstop=.true.
         return
c if we are storing large shells and next is smaller than min_shell - cut
      else if (icur_size.ge.min_shell.and.isize.lt.min_shell) then
         lstop=.true.
         return
c if we are merging large shells and next is large, but we have merged 2
c shells or the merged size is greater than max_merg - cut
      else if (icur_size.ge.min_shell.and.isize.ge.min_shell.and.
     *     (ihow_many.gt.max_how_many.or.(ik-ik_prev).ge.max_merg)) then
         lstop=.true.
         return
c if current is large but next is bigger than max_merg - cut
      else if (icur_size.ge.min_shell.and.isize.ge.max_merg) then
         lstop=.true.
         return
c if we are storing small shells (up to max_merg) but next is large - cut
c we do not merge small shells with large shells
      else if (icur_size.lt.min_shell.and.isize.ge.min_shell) then
         lstop=.true.
         return
      else
c Default:
        lstop=.false.
c     print *,'                                     ',
c    *  icur_size,isize,ik-ik_prev,ihow_many,lstop
      endif
      end
c==============================================================================
c not used
      subroutine strip_builder(str_reset,iwhich_shell,ishell_map,
     *                         ibasesymm,ibaseanti,ncf,
     *                         ishell_address_s,ishell_address_a,
     *                         ishell_index,igroup_no,imax_stripsize,
     *                         scr,iwrite_scpy,iwrite_acpy,tot_xmax)
c
c The AO exchange matrices are provided to EEO program as a continuous 
c stream of consecutive matrices with increasing ik index, where ik is the
c pair number for i and k index of matrices K(ik)=(ij|kl) for all j and l.
c For screening purposes this set of matrices is divided internally into 
c smaller pieces - stripes - in order to achieve better sparsity. Ideally 
c would be, for the best sparsity, if no shell merging were performed and all
c matrices were multiplied in shell blocks. This, unfortunately, extremely
c increases the amount of IO, because for each such shell full set of
c amplitudes would have to be read. Also matrix multiplication efficiency
c deteriorates significantly. The solution is to store integrals, as much as
c possible, and then, internally, to divide them into smaller stripes and to
c compress each such stripe independetly. Small 's' shells are merged into
c bigger stripes, but larger shells constitute single stripes. All stipes are
c kept in memory so that each read portion of T amplitudes is multiplied by all
c stripes. It requires, unfortunalety, that the same portion of T's must be
c compressed according to different pattern, because each stripe has different
c sparsity pattern.
c This subroutine constructs such stripes and all additional arrays which 
c contain different information about those stripes: their memory locations,
c locations of compression info for each stripe and the first and last ik
c index for each stripe. Things are more complicated here, because ik is not
c exactly the ik index as described above, but it is zeroed after each
c matrix multiplication, after which all integrals are deleted. So, it shows
c the matrix number in current stream after integral buffer flush and it also
c has its own array mapping this ik into real 'i' and 'k' index. The former
c array is used by caller and MultiplySort subroutine, but not here. We use
c 'current buffer' ik index here.
c I am not fully satisfied with the clearness of this subroutine and screening
c design, this is probably not the example of good programing. But somehow
c it works, although I still have doubts: why?
c There is some amount of static (saved) variables.
c The compression of integral matrices is performed here too.
c Parameters:
cI  str_reset        - this logical variable if true informs strip_builder
c                      that it should compress all matrices described by
c                      ishell_map array. This variable is set to true directly
c                      before MultiplySort call.
cI  iwhich_shell     - the number of shells "stored" by caller. This variable 
c                      could be 2 or 
c                      1. If it is 2 & str_reset=.false. it means that shell 
c                      1 is calculated and shell 2 isn't, but its size is 
c                      already known and may be used by this subroutine for 
c                      stripe border calculation.  If str_reset=.true.
c                      shell 2 is stored too, so we need to calculate stripe 
c                      border and compress everything.
cIO ishell_map       - this array contains the ik borders of shells. Remember:
c                      ik is zeroed after each MultiplySort call, so we do not
c                      really have real ik pair number.
cI  ibasesymm        - the base address of the storage area for symm K matrices
cI  ibaseanti        - ditto for antisimmetrized K matrices
cI  ncf              - number of contracted basis functions
cO  ishell_address_s - the array, shape 3 x number of stripes. For each 
c                      symmetrized stripe it contains:
c                      1: the bl index for starting address of this (compressed)
c                         stripe
c                      2: the number of rows of (compressed) stripe (up to 
c                         ntrisymm if there is no compression), if 0 it means
c                         that all matrix is compressed to zero size, if -1
c                         it means that no compression is performed
c                      3: the address, relative to address stored in 1 of the
c                         compression info array (integer*4 array stored in the
c                         memory space released by compressing the integrals),
c                         if 0 empty matrix nrow=0, if -1 no compression
c                      
cO  ishell_address_a - ditto for antisymmetric matrices
cO  ishell_index     - This array has dimension 2 x number of stripes. The 
c                      elements for each stripe contain:
c                      1: first ik for this stripe
c                      2: last ik for this stripe
cO  igroup_no        - the number of stripes constructed so far.
cO  imax_stripsize   - maximum size of stripe for all stripes constructed so far
cI  scr              - the screening matrix
cO  iwrite_scpy      - the bl index address where free memory space starts
c                      after compression (for symmetrized matrices)
cO  iwrite_acpy      - ditto for antisymmetrized matrices
cI  tot_xmax         - the maximum residuum element, used for 
      use memory
      use kinds
      implicit none
      logical str_reset
      integer iwhich_shell,ishell_map(2,*),ibasesymm,ibaseanti,ncf
      integer ishell_address_s(3,*),ishell_address_a(3,*),imax_stripsize
      integer*2 ishell_index(2,*)
      integer igroup_no,iwrite_scpy,iwrite_acpy
      real*8 scr(*),tot_xmax
c
      integer ihow_many_shells,   ifor_comp_s,   ifor_comp_a
      data    ihow_many_shells/0/,ifor_comp_s/0/,ifor_comp_a/0/
      save    ihow_many_shells,   ifor_comp_s,   ifor_comp_a
      integer icur_size,isize,ik,ik_old,ii,ntrisymm,ntrianti,irow
      integer iaddrs,iaddra,istart,istop,istripsize,itops,itopa,nrow
      logical cond1,i_am_resetted
      data          i_am_resetted/.true./
      save          i_am_resetted
      real*8 zeroed,zeroed_max,exzeros0,exzeros1,xzeros0,xzeros1
      common /zero_in_EEO/ zeroed,zeroed_max
      real*8 strace,sconstr,sort,extr,const,t0,t1,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
      common /timingstj/ strace,sconstr,sort,extr,const,c1,c2,
     *                   etrace,econstr,eort,eetr,eonst,e1,e2
c ishell_address_s,ishell_address_a description:
c 1 matrix addr.,                               <-----|
c 2 nrow =-1 if not compressed, 0 - zero matrix       |
c 3 compression info addr relative to element 1 of this array
c compresion info: if ishell_address(3,*)=-1, matrix is not compressed,
c 0 - zero matrix
 333  continue
      ntrisymm=ncf*(ncf+1)/2
      ntrianti=ntrisymm-ncf
      icur_size=ishell_map(2,1)-ishell_map(1,1)+1
      if (icur_size.le.0) then 
        print *, 'ishell_maps: ',ishell_map(2,1),ishell_map(1,1)
        call flush(6)
c       print '(F10.x)', icur_size
        STOP 'Do not expect icur_size .le. 0'
      endif
      if (iwhich_shell.ge.2) then 
        isize=ishell_map(2,2)-ishell_map(1,2)+1
        if (isize.le.0) then 
          print *, 'ishell_maps: ',ishell_map(2,2),ishell_map(1,2)
          call flush(6)
c         print '(F10.x)', icur_size
          STOP 'Do not expect isize .le. 0'
        endif
      else
        isize=0
      endif
      ik=ishell_map(2,1)
      if (isize.gt.0.and.ishell_map(1,2).ne.(ik+1)) 
     *                            STOP 'Error 5383hdjs84'
      if (i_am_resetted.and.ishell_map(1,1).ne.1) STOP 'Error 42gr7ct3'
      if (ishell_map(1,1).eq.1) then !     initialize
        i_am_resetted=.false.
        ifor_comp_s=ibasesymm
        ifor_comp_a=ibaseanti
        ihow_many_shells=0
        igroup_no=1
        ishell_index(1,igroup_no)=1        ! initial start
        ishell_index(2,igroup_no)=1        ! initial stop
c       print *, '* * * * * * start ik: ',1
        ishell_address_s(1,igroup_no)=ibasesymm
        ishell_address_s(2,igroup_no)=-1
        ishell_address_a(1,igroup_no)=ibaseanti
        ishell_address_a(2,igroup_no)=-1
        imax_stripsize=0
      endif
      ihow_many_shells=ihow_many_shells+1
c update external table:
      iwhich_shell=iwhich_shell-1
      do ii=1,iwhich_shell
        ishell_map(1,ii)=ishell_map(1,ii+1)
        ishell_map(2,ii)=ishell_map(2,ii+1)
      enddo
c
c     print *, 'current size, next size: ', icur_size,isize
      ik_old=ishell_index(1,igroup_no)
      call what_to_do(icur_size,isize,ik+1,ik_old,
     *                ihow_many_shells,cond1)
      if (cond1.or.isize.eq.0) then
        ishell_index(2,igroup_no)=ik ! stop
c here I compress: * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
        call elapsec(exzeros0)
        call secund(xzeros0)
c
        call dynamic_getmem(ntrisymm/i4size+1,irow)
        iaddrs=ishell_address_s(1,igroup_no)
        iaddra=ishell_address_a(1,igroup_no)
        istart=ishell_index(1,igroup_no)
        istop =ishell_index(2,igroup_no)
        istripsize=istop-istart+1
        if (istripsize.le.0) STOP 'Error 6h37fhew'
        if (istripsize.ge.2000) STOP 'Error 5h37g63'
        if (istripsize.gt.imax_stripsize) imax_stripsize=istripsize
c 1 matrix addr., 
c 2 nrow =-1 if not compressed, 0 all matrix is eq. 0.0d0
c 3 compression info addr relative to 1, 
        if (ifor_comp_s.ne.iaddrs) STOP 'Error ifor_comp_s'
        if (ifor_comp_a.ne.iaddra) STOP 'Error ifor_comp_a'
        call EEORowCompactor(bl(iaddrs),ntrisymm,istripsize,
     *                     scr,bl(irow),ishell_address_s(2,igroup_no),
     *                     ishell_address_s(3,igroup_no),itops,tot_xmax)
        nrow=ishell_address_s(2,igroup_no)
        if (itops.gt.ntrisymm*istripsize) STOP 'squeeze error'
        if (nrow.ne.-1) zeroed=zeroed+(ntrisymm-nrow)*istripsize
        call EEORowCompactor(bl(iaddra),ntrianti,istripsize,
     *           scr(ntrisymm+1),bl(irow),ishell_address_a(2,igroup_no),
     *           ishell_address_a(3,igroup_no),itopa,tot_xmax)
        if (itopa.gt.ntrianti*istripsize) STOP 'squeeze error'
        nrow=ishell_address_a(2,igroup_no)
        if (nrow.ne.-1) zeroed=zeroed+(ntrianti-nrow)*istripsize
        call dynamic_retmem(1)
        call elapsec(exzeros1)
        call secund(xzeros1)
        extr=extr+xzeros1-xzeros0
        eetr=eetr+exzeros1-exzeros0
c END OF compress: * * * * * * * * * * * * * * * * * * * * * * * * * * * * 
c       print *, '* * * * * * stop ik: ',ik
        if (isize.ne.0) then 
c         print *, '* * * * * * start ik: ',ik+1
          igroup_no=igroup_no+1
          ishell_index(1,igroup_no)=ik+1   ! start
          ishell_address_s(1,igroup_no)=ifor_comp_s+itops
          ishell_address_a(1,igroup_no)=ifor_comp_a+itopa
c         call elapsec(exzeros0)
c         call secund(xzeros0)
c         if ((ifor_comp_s+ntrisymm*istripsize).ne.(ifor_comp_s+itops))
c    *    then
c        call move(bl(ifor_comp_s+ntrisymm*istripsize),
c    *             bl(ifor_comp_s+itops),ntrisymm*isize)
c         endif
c         if ((ifor_comp_a+ntrianti*istripsize).ne.(ifor_comp_a+itopa))
c    *    then
c        call move(bl(ifor_comp_a+ntrianti*istripsize),
c    *             bl(ifor_comp_a+itopa),ntrianti*isize)
c         endif
c         call elapsec(exzeros1)
c         call secund(xzeros1)
c         extr=extr+xzeros1-xzeros0
c         eetr=eetr+exzeros1-exzeros0
          ifor_comp_s=ifor_comp_s+itops
          ifor_comp_a=ifor_comp_a+itopa
          iwrite_scpy=ifor_comp_s
          iwrite_acpy=ifor_comp_a
          ishell_address_s(2,igroup_no)=-1
          ishell_address_a(2,igroup_no)=-1
        else                                   !  if (isize.ne.0)
          iwrite_scpy=-1
          iwrite_acpy=-1
        endif                                  !  if (isize.ne.0)
        ihow_many_shells=0
      endif
      if (str_reset.and.iwhich_shell.gt.0) then
        goto 333
      endif
      if (str_reset.and.iwhich_shell.eq.0) then
        i_am_resetted=.true.
      endif
      end
c==============================================================================
      subroutine set_thresh(tot_xmax,thresh)
      implicit none
      real*8 tot_xmax,thresh
c tmp switched off
      thresh=1.d0
      return
c END of tmp switched off
c
c
      thresh=dabs(tot_xmax*1.d4)
c
      if (thresh.gt.1.d3) then
        thresh=1.d3
        return
      endif
      if (thresh.lt.1.d0) then
        thresh=1.d0
        return
      endif
c
      return
c OLD CODE:
      if      (tot_xmax.gt.1.d0 ) then! for integral thresh=10**-10 it gives
        thresh=1.d6                   ! 10**-4
      else if (tot_xmax.gt.1.d-1) then! 
        thresh=1.d5                   ! 10**-5
      else if (tot_xmax.gt.1.d-2) then
        thresh=1.d4                   ! 10**-6
      else if (tot_xmax.gt.1.d-3) then
        thresh=1.d3                   ! 10**-7
      else if (tot_xmax.gt.1.d-4) then
        thresh=1.d2                   ! 10**-8
      else if (tot_xmax.gt.1.d-5) then
        thresh=1.d1                   ! 10**-9
      else if (tot_xmax.gt.1.d-6) then
        thresh=1.d0                   ! 10**-10
      else if (tot_xmax.gt.1.d-7) then
        thresh=1.d0                   ! 10**-10
      else if (tot_xmax.gt.1.d-8) then
        thresh=1.d0                   ! 10**-10
      else
        thresh=1.d0                   ! 10**-10
      endif
      end
c==============================================================================
      subroutine calculate_percent(modifier,idimen,percent_factor)
c This subroutine tries to calculate the fraction (range 0 to 1) which
c might be used for assessment what fraction of integral matrix must be non-zero
c to retain the same speed of multiplication as without compression. If
c smaller fraction of integral matrix is nonzero, we get speed-up with
c multiplication. If larger fraction of integral matrix is nonzero, the
c compression + multiplication would take more time than multiplication of 
c non-compressed matrices.
c This fraction is a function of idimen and parameters described below
c INTENT(IN)
c modifier - if 1, we return percent_factor according to the above
c            descr., if 2, we artificially increase speed of reading by 2, 
c            if 3 we increase
c            speed of reading by 3 and so on, which increses percent_factor. It
c            might be useful for small memory jobs, where compresion due to disk
c            I/O is more important than mult. efficiency.
c            if modifier >= 1 000 000 return percent_factor == 1.0d0
c idimen    -  matrix dimension 
c INTENT(OUT)
c percent_factor - see descr. above, fraction in (0.0d0,1.0d0) range
      implicit none
      real*8 percent_factor
      integer modifier,idimen
      real*8 a,b,d,s,speed
      parameter (a=3486.97d0,b=7.24354d0,d=-122.248d0)
      if (modifier .ge. 1 000 000) then
        percent_factor=1.0d0
        return
      endif
c Speed of multiplication is function of integral matrix dimension and
c assuming that amplitudes matrix dimension is constant (25), it is
c independent of other dimensions. It is fitted by function:
c speed (flops/s)= (a*idim/(b+idim)+d)*10**6
c
      s=85 000 000d0*dble(modifier) ! speed of read. real*8(doublewords)/sec
      speed=(a*idimen/(b+idimen)+d)*10**6
      percent_factor= 2*s*dble(idimen)/(2*s*dble(idimen)+speed)
      if (percent_factor.lt.0.0d0.or.percent_factor.gt.1.0d0) then
        print *, 'Error 5490348rhfdk38',idimen,a,b,d,s,percent_factor
        STOP     'Error 5490348rhfdk38'
      endif
      end
c==============================================================================
      subroutine schwarz_calc(inx,ncs,schwarz)
c This subroutine calculates schwarz matrix, used for screening. The matrix
c element i,j contains maximum value of all (ii,jj|ii,jj) integrals where
c ii is basis function which belongs to i contracted shell and jj is function 
c in j shell.
cI  inx     - the contraction info
cI  ncs     - number of contracted shells
cO  schwarz - the result
      use memory
      implicit none
      integer inx(12,*),ncs
      real*8 schwarz(ncs,ncs)
c
      real*8 eps,xmax,perm,xint
      integer ibuf,ics,isize,jcs,jsize,integr,istart,istop,jstart,jstop
      integer ishift,jshift,kshift,icf,jcf,indexx,i,j,iii
      logical nonz
c
      call zeroit(schwarz,ncs*ncs)
c
      call dynamic_mmark
c
      eps=0.0d0
      call dynamic_getmem(15**4,ibuf)
      call dynamic_getmem(15**4,integr)
      do ics=1,ncs
        isize=inx(3,ics)
        do jcs=1,ics
          jsize=inx(3,jcs)
          call twoelmain(ics,jcs,ics,jcs,eps,nonz,bl(integr),bl(ibuf))
          if (.not.nonz) then
            cycle
          endif
          istart=inx(11,ics)+1
          istop =inx(10,ics)
          jstart=inx(11,jcs)+1
          jstop =inx(10,jcs)
          ishift=isize*isize*jsize
          jshift=isize*jsize
          kshift=jsize
          xmax=0.0d0
c
          perm=2.0d0
          if (ics.eq.jcs) perm=8.0d0
c
          do iii=0,isize*isize*jsize*jsize-1
              xint=bl(integr+iii)*perm
c             if (ics.eq.20.and.jcs.eq.2) write(*,'(E20.8E2)') xint/perm
              if (dabs(xint).gt.xmax) xmax=dabs(xint)
          enddo
c
          schwarz(ics,jcs)=xmax
          if (ics.ne.jcs) schwarz(jcs,ics)=xmax
c
        enddo
      enddo
      call dynamic_retmark
      end
c==============================================================================
      subroutine energy_contr(exc,coef,idimen,epair)
      implicit none
      integer idimen
      real*8 exc(idimen,idimen),coef(idimen,idimen),epair
      integer i,j
      epair=0d0
      do i=1,idimen
        do j=1,idimen
          epair=epair+(2d0*exc(j,i)-exc(i,j))*coef(j,i)
        enddo
      enddo
      end
c==============================================================================
      subroutine energy_contr_r(exc,coef,resi,idimen,epair)
      implicit none
      integer idimen
      real*8 exc(idimen,idimen),coef(idimen,idimen),epair
      real*8 resi(idimen,idimen),xij,xji
      integer i,j
      epair=0d0
      do i=1,idimen
        do j=1,idimen
          xij=exc(i,j)+resi(i,j)
          xji=exc(j,i)+resi(j,i)
          epair=epair+(2d0*xji-xij)*coef(j,i)
        enddo
      enddo
      end
c==============================================================================
      subroutine bzeroit(xmatrix,isize)
      implicit none
      integer*1 xmatrix(*)
      integer isize
c
      integer i
c
      do i=1,isize
        xmatrix(i)=0
      enddo
      end
c===============================================================================
      subroutine ijresolver(nval,ij,i,j)
      implicit none
      integer ij,i,j,nval
      i=mod(ij,nval)
      if (ij.lt.0.or.nval.lt.0) STOP 'ij.lt.0.or.nval.lt.0'
      if (i.eq.0) i=nval
      j=(ij-i)/nval+1
      if (i.gt.nval.or.j.gt.nval) then
        print *,i,j,ij
        STOP 'i.gt.nval.or.j.gt.nval'
      endif
      end
c===============================================================================
      subroutine allclose_and_delete(af)
      implicit none
      integer i
      logical op,af
      character*256 na
        ! single proc:
        if (af) then
        close(76,STATUS='delete')  ! pqs_scr_propane.kijk-???
        close(99,STATUS='delete')  ! pqs_scr_propane.debug_slaves
        else
        close(41,STATUS='delete')  ! pqs_scr_propane.coul.bins
        close(42,STATUS='delete')  ! pqs_scr_propane.exc.bins
        close(44,STATUS='delete')  ! pqs_scr_propane.tx.bins
        close(45,STATUS='delete')  ! pqs_scr_propane.tc.bins
        close(46,STATUS='delete')  ! pqs_scr_propane.tt.bins
        close(51,STATUS='delete')  ! pqs_scr_propane.Ta
        close(53,STATUS='delete')  ! pqs_scr_propane.Jc
        close(54,STATUS='delete')  ! pqs_scr_propane.Kext_ab
        close(55,STATUS='delete')  ! pqs_scr_propane.Kbc1
        close(60,STATUS='delete')  ! pqs_scr_propane.ampl.001
        close(63,STATUS='delete')  ! pqs_scr_propane.Cij
        close(64,STATUS='delete')  ! pqs_scr_propane.EEOc
        close(65,STATUS='delete')  ! pqs_scr_propane.TXij
        close(66,STATUS='delete')  ! pqs_scr_propane.TCij
        close(68,STATUS='delete')  ! pqs_scr_propane.xkijkl
        close(70,STATUS='delete')  ! pqs_scr_propane.YZ
        close(71,STATUS='delete')  ! pqs_scr_propane.41a
        close(72,STATUS='delete')  ! pqs_scr_propane.Qparts
        close(73,STATUS='delete')  ! pqs_scr_propane.TTij
        close(74,STATUS='delete')  ! pqs_scr_propane.resid
        close(75,STATUS='delete')  ! pqs_scr_propane.diis
        close(76,STATUS='delete')  ! pqs_scr_propane.kijk
        close(98,STATUS='delete')  ! pqs_scr_propane.Kij
        endif
        return
c   Not used now, but may be needed: it is an automatic detector of files
c   that are still opened. Its use is risky in general, because it may
c   delete files that are still needed (e.g. from other PQS modules). As a
c   programmer's tool it is very useful.
        do i=10,99
          inquire (unit=i,OPENED=op)
          if (op) then
            inquire (unit=i,NAMED=op)
            if (op) then
              inquire (unit=i,NAME=na)
              write(6,*) i,', deleting ',na(1:len_trim(na))
              na=''
            else
              write(6,*) i,', deleting unnamed file'
            endif
            call flush(6)
            close(i,STATUS='delete',ERR=222)
            goto 333
 222        continue
            write(6,*) '                            Closing failed'
 333        continue
          endif
        enddo
      end

#include "maxi.h"
      function alpha2(k,      l,          i,          j,    iRstart,
     *                iRstop, ialphabase, ndiskalpha, nval, reset,
     *                af)
c This function evolved from something which had different purposes previously.
c The purpose of this function is to return single double precision value of
c alpha parameter. Because alpha parameter is stored on disk in the form that
c for given PAIR ij all k and l are stored as matrix m(k,l) in a single
c record, it must somehow avoid reading the whole records for obtaining
c particular k,l element from it.  This function knows, that callers has
c external loop over ij from iRstart to iRstop, and all alphas from this ij
c range will be needed before any other alphas from different ij range. So,
c it reads all alphas which belong to iRstart,iRstop into memory location
c ialphabase which should be capable to keep (iRstop-iRstart+1)*nval*nval
c numbers. This memory is allocated elsewhere. As long as the iRstart & iRstop 
c are the same numbers, it simply returns alpha number from its cache. If
c iRstart changes, it reads new alphas into cache, and then it returns them
c until iRstart changes and so on. The reset variable is set to .true. at the
c beginning of alpha calculations in order to inform this subroutine that 
c new iteration started and the alphas must be read first time. This additional
c variable is needed because, with large memory and small system, it is very 
c likely that iRstart will be always equal to 1.
cI  k,l,i,j    - alpha indices
cI  iRstart    - the external subroutine has external loop over ij in this range
c                starting from iRstart
cI  iRstop     - ditto, last element of ij range
cI  ialphabase - the alphas cache
cI  ndiskalpha - the file unit number, where alphas are located
cI  nval       - number of correlated orbitals
cIO reset      - if .true. it means that new iteration just started. This 
c                subroutine immediately sets this variable to .false. If next
c                call brings reset=.true. again, it means that new iteration
c                started.
cI  af         - if .true. the job is parallel & Array Files is used as IO
c                device, otherwise job is single processor and local IO is used.
c
      use memory
      use ccounters
      implicit none
      real*8 alpha2
      integer k,l,i,j,iRstart,iRstop,ialphabase,ndiskalpha,nval
      logical reset,af
c
      integer nval2,ij,ijrec,icuraddr,idata,istatus
      logical reading
      nval2=nval*nval
      if (iRstart.ne.iRstartold.or.reset) then
        reading=.true.
        reset=.false.
        iRstartold=iRstart
      else
        reading=.false.
      endif
      if (i.ge.j) then
        ij=i*(i-1)/2+j
      else
        ij=j*(j-1)/2+i
      endif
      if (ij.lt.iRstart.or.ij.gt.iRstop) STOP 'alpha2 error'
c read into cache:
      if (reading) then
        do ijrec=iRstart,iRstop
          icuraddr=ialphabase+(ijrec-iRstart)*nval2
          if (af) then
            call fafread(ndiskalpha,bl(icuraddr),8,nval*nval,1,ijrec,
     *                   istatus)
          if (istatus.lt.8*nval*nval) then
              ERROR
              call flush(6)
          endif
          else
            call reader(ndiskalpha,nval,ijrec,bl(icuraddr))
          endif
        enddo
      endif
c
      idata=ialphabase+(ij-iRstart)*nval2+(k-1)*nval+(l-1)
      alpha2=bl(idata)
      end
C===============================================================================
      subroutine CCA(ncf,iA,nmo,vorb,nbf)
c This subroutine calculated A matrix (J. Chem. Phys. 81, 1901 (1984)) from 
c scrach each time when amplitudes were updated. Right now it only transfers
c this matrix, calculated elsewehere, only once per iteration.
c Parameters:
cI  ncf  - number of cotracted basis functions
cO  iA   - the result is stranferres here
cI  nmo  - number of occupied HF orbitals
cI  vorb - .true.: the dimension of virtual space is equal to no of virt.
c          orbitals, and MO form of integral matrices is used, otherwise,
c          the virtual space has ncf dimension and all calculations are 
c          performed for matrices in AO form.
      use memory
      implicit none
      integer ncf,iA,nmo,nbf
      real*8 half
      integer iexchAO,icoeffAO,k,l,mataddr,iaddress,idimen
      logical doit,byt8,vorb
      parameter (half=0.5d0)
      common /CCAOnceCom/ iaddress,doit
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
        call tfer(bl(iaddress),bl(iA),idimen*idimen)
      end
C==============================================================================
      subroutine CCLoperator(xmatrix,idimen,xloperator)
c Simple subroutine which builds the matrix defined as:
c L(i,k)=2*X(i,k)-X(k,i)
      implicit none
      integer idimen
      REAL*8 xmatrix(idimen,idimen),xloperator(idimen,idimen)
      integer i,j
      do j=1,idimen
        do i=1,idimen
          xloperator(i,j)=2.0d0*xmatrix(i,j)-xmatrix(j,i)
        enddo
      enddo
      end
C==============================================================================
      subroutine CCalphaonce(ncf,   nval,   irecadrx, npairs,    ndiskx,
     *                       lbinx, thresh, doit,     byt8,      ccsd,
     *                       nmo,   vorb,   elaps,    ndiskalpha,af,
     *                       nbf)
c This subroutine calculates part of alpha parameters. It calculates the
c most expensive part which involves calculation of <T*K> traces. The rest
c of alphas is added elsewhere to the results of this subroutine.
c The caclulations are performed in blocks, the matrices are logically 
c reshaped and all matrix traces are turned into matrix multiplications.
c Parameters:
c   ncf        - number of contracted basis functions
c   nval       - number of correlated orbitals
c   irecadrx   - the record info for file containing Kext bins. Note, that 
c                this array has different size and different information is
c                contained in it in single processor mode and in parallel mode.
c   npairs     - number of correlated orbital pairs
c   ndiskx     - the file unit number cotaining Kext integrals
c   lbinx      - the bin length for Kext
c   thresh     - integral threshold
c   doit       - have to be always .true. !!! It was used earlier, but since
c                alphas are not recalculated after each ij amplitude matrix
c                update, it is not expected to be .false.
c   byt8       - .true.: double precision bins are used. .false.: int4+int1
c                bins are used
c   ccsd       - .true.: the specific for CCSD singles contributions are 
c                calculated. This variable also imply that all terms needed for
c                CCSD are calculated: CI singles, CCD, QCISD singles etc.
c   nmo        - number of occupied HF orbitals
cI  vorb       - .true.: the dimension of virtual space is equal to no of virt.
c                orbitals, and MO form of integral matrices is used, otherwise,
c                the virtual space has ncf dimension and all calculations are 
c                performed for matrices in AO form.
c   elaps      - this array collects different timings from this subroutine
c   ndiskalpha - the file unit number where alphas are deposited. The alphas
c                are stored as records, for each i,j PAIR all k and l are
c                written as one records. The order of k,l element follow 
c                FORTRAN rules of memory occupation.
cI  af         - if .true. the job is parallel & Array Files is used as IO
c                device, otherwise job is single processor and local IO is used.
c
      use memory
      implicit none
      integer ncf,nval,irecadrx(*),npairs,ndiskx,lbinx,nmo,ndiskalpha
      REAL*8 thresh,alpha,elaps(*)
      integer iexchAO,icoeffAO,mataddr,i,j,k,l,iaddress,nbf
      integer inval,icurraddr,nval2,nval3,initial_chunk
      logical doit,idoit,byt8,ccsd,vorb,af
      parameter (initial_chunk=25)
      REAL*8 elem,elaps0,elaps1
      common /CCAlphaCom/iaddress,inval,idoit
c new
      integer memavail,lcore,iad,ncf2,number_matrices,ichunk_T,ichunk_K
      integer iKstore,iTstore,iresult,iKpass,iTpass,iT,iK,iTstart,iTstop
      integer iTsize,iKstart,iKstop,iKsize,ij,kl,number_matrices1,idimen
      integer idimen2,ialphabuffer,itransposition
c
      save /CCAlphaCom/
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      idimen2=idimen*idimen
c
      inval=nval
      idoit=doit
      ncf2=ncf*ncf
      nval2=nval*nval
      nval3=nval*nval*nval
      if (doit) then
      call mmark
      call dynamic_mmark
      call matmark
c
c  Calculate space available for K and T blocks.
c
      call dynamic_show_free(memavail)
c
      number_matrices=memavail/(idimen2) ! +1 for result matrix
c If all T matrices and initial_chunk K matrices will fit the memory
      looper: do
        if ((npairs+initial_chunk).le.number_matrices) then
         ichunk_T=npairs
         number_matrices1=number_matrices-ichunk_T
         if (number_matrices1.gt.npairs) then 
           ichunk_K=npairs
         else
           ichunk_K=number_matrices1
         endif
        else
          ichunk_T=number_matrices-initial_chunk
          ichunk_K=initial_chunk
        endif
        if (ichunk_T.lt.initial_chunk.or.ichunk_K.lt.initial_chunk) then
          if (ichunk_T.ne.npairs.and.ichunk_K.ne.npairs) then
            print *, 'ichunk_T: ',ichunk_T
            print *, 'ichunk_K: ',ichunk_K
            print *, 'idimen: ',idimen
            print *, 'Too little memory in alpha builder'
            STOP 'Too little memory in alpha builder'
          endif
        endif
        if ((idimen2*(ichunk_T+ichunk_K)+ichunk_T*ichunk_K+
     *                             ichunk_T*nval*nval).le.memavail)
     *    then
          exit
        else
         number_matrices=number_matrices-1
        endif
      enddo looper
c
      call dynamic_getmem(idimen2*ichunk_K,iKstore)
      call dynamic_getmem(idimen2*ichunk_T,iTstore)
      call dynamic_getmem(ichunk_T*ichunk_K,iresult)
      call dynamic_getmem(ichunk_T*nval*nval,ialphabuffer)
c
      iKpass=npairs/ichunk_K
      if (mod(npairs,ichunk_K).ne.0) iKpass=iKpass+1
      iTpass=npairs/ichunk_T !!
      if (mod(npairs,ichunk_T).ne.0) iTpass=iTpass+1
c Loops over passes:
      do iT=1,iTpass
        iTstart=ichunk_T*(iT-1)+1
        iTstop =ichunk_T*iT
        if (iTstop.gt.npairs) iTstop=npairs
        iTsize=iTstop-iTstart+1
        call elapsec(elaps0)
        call store(iTstart,iTstop,idimen,nval,iTstore)
        call elapsec(elaps1)
        elaps(1)=elaps(1)+elaps1-elaps0 ! amplitudes storing
        call matconn('BigT','r',idimen2,iTsize,iTstore)
        do iK=1,iKpass
          iKstart=ichunk_K*(iK-1)+1
          iKstop =ichunk_K*iK
          if (iKstop.gt.npairs) iKstop=npairs
          iKsize=iKstop-iKstart+1
          call elapsec(elaps0)
          call storeK(irecadrx,ndiskx,lbinx,thresh,iKstart,iKstop,
     *                ncf,nval,iKstore,byt8,nmo,vorb,nbf)
          call elapsec(elaps1)
          elaps(2)=elaps(2)+elaps1-elaps0 ! Kij storing
          call matconn('BigK','r',idimen2,iKsize,iKstore)
          do itransposition=1,2
          call elapsec(elaps0)
          if (itransposition.eq.2) then
            call BigTransposition(iKstart,iKstop,idimen,nval,iKstore)
          endif
          call matconn('result','r',iTsize,iKsize,iresult)
          call matmmul2('BigT','BigK','result','t','n','n')
          call elapsec(elaps1)
          elaps(3)=elaps(3)+elaps1-elaps0 ! multiplication
          call elapsec(elaps0)
          call putresult(idimen,nval,bl(iresult),iTstart,iTstop,iKstart,
     *                   iKstop,ccsd,ialphabuffer,ndiskalpha,af,
     *                   itransposition)
          call elapsec(elaps1)
          elaps(4)=elaps(4)+elaps1-elaps0 ! putresult
          call matdisc('result')
          enddo
          call matdisc('BigK')
        enddo
        call matdisc('BigT')
      enddo
c
      call matremark
      call dynamic_retmark
      call retmark
      endif
      end
C==============================================================================
      subroutine CCAOnce(irecadrx, npairs, nval, ndiskx, ncf,
     *                   lbinx,    thresh, doit, byt8,   ioverlap,
     *                   nmo,      vorb,   nbf)
c This subroutine calculates A matrix, as it is defined in 
c J. Chem. Phys. 81, 1901 (1984) - formula no. (38)
c Parameters:
c   irecadrx - the record info for file containing Kext bins. Note, that 
c              this array has different size and different information is
c              contained in it in single processor mode and in parallel mode.
c   npairs   - number of correlated pairs
c   nval     - number of correlated orbitals
c   ndiskx   - file unit number where Kext bins are kept
c   ncf      - number of contracted functions
c   lbinx    - Kext bin length
c   thresh   - integral threshold
c   doit     - it is passed to common /CCAOnceCom/, but not used anywhere.
c   byt8     - .true.: double precision bins are used. .false.: int4+int1
c              bins are used
c   ioverlap - the overlap matrix bl index
c   nmo      - number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      use ccounters
      implicit none
      integer irecadrx(*),npairs,nval,ndiskx,ncf,lbinx,iA,ioverlap,nmo
      integer nbf
      real*8 thresh
      real*8 half
      integer iexchAO,icoeffAO,k,l,mataddr,idimen
      logical doit,idoit,byt8,vorb
      parameter (half=0.5d0)
      common /CCAOnceCom/ iA,idoit
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ic_CCAOnce=ic_CCAOnce+1
      idoit=doit
      if (ic_CCAOnce.eq.1) call dynamic_getmem(idimen*idimen,iA)
      call mmark
      call dynamic_mmark
      call matmark
      call matconn('A','q',idimen,idimen,iA)
      call matzero('A')
      call dynamic_matdef('exchAO','q',idimen,idimen)
      call dynamic_matdef('coetilda','q',idimen,idimen)
      call dynamic_matdef('work33','q',idimen,idimen)
      iexchAO=mataddr('exchAO')
c A is calculated according to the recipe taken from Pualy article. It is at
c the end transposed and multiplied by overlap in order to be identical
c with A from Werner article.
      do k=1,nval
        do l=1,k
          call ExtrOne(k,    l,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(iexchAO))
          call CoefRead('tt',k,l,icoeffAO)
          call matconn('coeAO','q',idimen,idimen,icoeffAO)
          call matzero('coetilda')
          call matadd1('coeAO',+4.0d0,'coetilda')
          call matpose('coeAO')
          call matadd1('coeAO',-2.0d0,'coetilda')
          call matpose('coeAO')
          call matdisc('coeAO')
          if (k.eq.l) call matscal('coetilda',half)
          call matmmul2('coetilda','exchAO','work33','t','n','n')
          call matadd1('work33',half,'A')
          call matmmul2('coetilda','exchAO','work33','n','t','n')
          call matadd1('work33',half,'A')
        enddo
      enddo
! A is changed in order to be the same as in Werner article
      call matpose('A')
      if (vorb) then
        continue
c       call matcopy('work33','A')
      else
        call matcopy('A','work33')
        call matconn('overl','q',ncf,ncf,ioverlap)
        call matmmult('work33','overl','A')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
C==============================================================================
      subroutine putresult(ncf,     nval,   resul, iTstart,  iTstop,
     *                     iKstart, iKstop, ccsd,  iaddress, ndisk,
     *                     af,      itransp)
c This subroutine caches data obtained in CCalphaonce (the alpha parameters).
c When the whole records are collected, it flushes cache into file. The record
c is defined as a matrix composed from all alphas with common i,j index 
c (all k,l). The alphas are defined by J. Chem. Phys. 81, 1901 (1984).
c Parameters:
cI  ncf       - number of contracted functions
cI  nval      - number of correlated orbitals
cI  resul     - the partial result passed by CCalphaonce, this result has
c               ij pair range iTstart,iTstop & kl pair range iKstart,iKstop, 
c               but if itransp=2 the k,l pair is considered as l,k pair.
cI  iTstart   - see explanation for resul
cI  iTstop    - see explanation for resul
cI  iKstart   - see explanation for resul
cI  iKstop    - see explanation for resul
cI  ccsd      - .true.: the specific for CCSD singles contributions are 
c               calculated. This variable also imply that all terms needed for
c               CCSD are calculated: CI singles, CCD, QCISD singles etc.
cIO iaddress  - the bl index address for cache. It has to have size
c               at least (iTstop-iTstart+1)*nval*nval
cI  ndisk     - the file unit number where alphas are stored
cI  af        - if .true. the job is parallel & Array Files is used as IO
c               device, otherwise job is single processor and local IO is used.
cI  itransp   - if itransp.eq.1 kl means kl pair number for k.ge.l
c               if itransp.eq.2 kl it means that k and l are reversed for the 
c               same pair number, i.e. kl is pair number for l.ge.k
      use memory
      implicit none
      integer ncf,nval,iTstart,iTstop,iKstart,iKstop,iaddress,ndisk
      integer itransp
      logical ccsd,af
      real*8 resul(iTstop-iTstart+1,*)
      real*8 element
      integer iB,jB,iE,jE,lB,kB,lE,kE,ij,kl,ijmatrix,klmatrix
      integer i,j,k,l,imap,jmap,icurraddr,nval2,nval3,icopy
      integer icoefvi,ikopv,mataddr,ijaddress,istatus,kl_orig,npairs
      nval3=nval*nval*nval
      nval2=nval*nval
      npairs=nval*(nval+1)/2
c
      iB=imap(iTstart)
      iE=imap(iTstop)
c
      kB=imap(iKstart)
      kE=imap(iKstop)
c
cold  lB=mod(iKstart,nval)
cold  if (lB.eq.0) lB=nval
cold  kB=(iKstart-lB)/nval+1 ! K is in the outer loop, but must be
c                            calculated from l
cold  lE=mod(iKstop,nval)
cold  if (lE.eq.0) lE=nval
cold  kE=(iKstop-lE)/nval+1
c
      ij=(iB-1)*iB/2
      ijmatrix=0
      outermost: do i=iB,iE
        do j=1,i
          ij=ij+1
          if (ij.lt.iTstart) cycle
          if (ij.gt.iTstop) exit outermost
          ijmatrix=ijmatrix+1
          klmatrix=0
cold      kl=nval*(kB-1)
          kl=(kB-1)*kB/2
          ijaddress=iaddress+(ijmatrix-1)*nval*nval
          outer: do k=kB,kE
            do l=1,k
             kl=kl+1
             if (kl.lt.iKstart) cycle
             if (kl.gt.iKstop) exit outer
             if (itransp.eq.1) then 
               kl_orig=(k-1)*nval+l
             else
               kl_orig=(l-1)*nval+k
             endif
             klmatrix=klmatrix+1
             icurraddr=ijaddress+(kl_orig-1)
             bl(icurraddr)=resul(ijmatrix,klmatrix)
            enddo
          enddo outer
          if (iKstop.eq.npairs.and.itransp.eq.2) then  !    flush
          if (af) then
            call fafwrite(ndisk,bl(ijaddress),8,nval*nval,1,ij,istatus)
          else
            call writer(ndisk,nval,ij,bl(ijaddress))
          endif
          endif
        enddo
      enddo outermost
      end
C==============================================================================
      subroutine Builder41a_G(ncf,    nval,       norecalc, cc,   ccsd,
     *                        nfirst, F,          nmo,      vorb, ndisk,
     *                        af,     ndiskalpha, betam,    nbf)
c This subroutine calculates  the sum over kl in formula (39)
c J. Chem. Phys. 104, 6286 (1996)
c It does it as matrix multiplication instead of matrix summation,
c calculation goes in blocs of matrices, which are logically reashaped
c Parameters:
c   ncf        - number of contracted functions
c   nval       - number of correlated orbitals
c   norecalc   - not used any more
c   cc         - still passed to storeAlpha, but nod needed here (legacy param)
c   ccsd       - ditto
c   nfirst     - ditto
c   F          - ditto
c   nmo        - the number of occupied HF orbitals
cI  vorb       - .true.: the dimension of virtual space is equal to no of virt.
c                orbitals, and MO form of integral matrices is used, otherwise,
c                the virtual space has ncf dimension and all calculations are 
c                performed for matrices in AO form.
c   ndisk      - the result is written here
cI  af         - if .true. the job is parallel & Array Files is used as IO
c                device, otherwise job is single processor and local IO is used.
c   ndiskalpha - alphas are kept here
c   betam      - the beta matrix, as defined by (41) the article referred above
      use memory
      implicit none
      integer ncf,nval,nfirst,nmo,ndisk,ndiskalpha,nbf
      logical norecalc,cc,ccsd,vorb,af
C
      integer iad,lcore,memavail,initial_chunk,number_matrices,ncf2
      integer npairs,ichunk_T,number_matrices1,ichunk_R,iTstore,iRstore
      integer ialphas,iT,iTpass,iTstart,iTstop,iTsize,iR,iRpass
      integer iRstart,iRstop,iRsize,itwo,idimen,idimen2,ialphabase
      logical transposed,reset
      real*8 F(ncf,ncf),betam(nval,nval)
      parameter (initial_chunk=25)
C
c     if (.not.norecalc) return
c
      reset=.true. ! for storeAlpha
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      ncf2=ncf*ncf
      idimen2=idimen*idimen
      npairs=nval*(nval+1)/2
      call dynamic_show_free(memavail)
c We need memory for alpha matrix and C matrix and result matrix.
c How much we can get?
      number_matrices=memavail/(idimen2)
c If all T matrices and initial_chunk Res matrices will fit the memory
      looper: do
        if ((npairs+initial_chunk).le.number_matrices) then
         ichunk_R=npairs
         number_matrices1=number_matrices-ichunk_R
         if (number_matrices1.gt.npairs) then 
           ichunk_T=npairs
         else
           ichunk_T=number_matrices1
         endif
        else
          ichunk_R=number_matrices-initial_chunk
          ichunk_T=initial_chunk
        endif
        if (ichunk_R.lt.initial_chunk.or.ichunk_T.lt.initial_chunk) then
          if (ichunk_R.ne.npairs.and.ichunk_T.ne.npairs) then
            print *, 'ichunk_R: ',ichunk_R
            print *, 'ichunk_T: ',ichunk_T
            print *, 'idimen: ',idimen
            print *, 'Too little memory in 41a builder'
            STOP 'Too little memory in 41a builder'
          endif
        endif
        if ((idimen2*(ichunk_R+ichunk_T)+ichunk_R*ichunk_T+
     *        ichunk_R*nval*nval).le.memavail)
     *    then
          exit
        else
         number_matrices=number_matrices-1
        endif
      enddo looper
c
      call dynamic_getmem(idimen2*ichunk_T,iTstore)
      call dynamic_getmem(idimen2*ichunk_R,iRstore)
      call dynamic_getmem(ichunk_T*ichunk_R,ialphas)
      call dynamic_getmem(nval*nval*ichunk_R,ialphabase)
c Changed above check!
c
c  Initialize file for 41a matrices
      call Fileinit41a(idimen,ndisk,af)
c calculate passes:
      iRpass=npairs/ichunk_R
      if (mod(npairs,ichunk_R).ne.0) iRpass=iRpass+1
      iTpass=npairs/ichunk_T !!
      if (mod(npairs,ichunk_T).ne.0) iTpass=iTpass+1
c Loops over passes:
      do iR=1,iRpass
        iRstart=ichunk_R*(iR-1)+1
        iRstop =ichunk_R*iR
        if (iRstop.gt.npairs) iRstop=npairs
        iRsize=iRstop-iRstart+1
        call matconn('BigPrR','r',idimen2,iRsize,iRstore)
        call matzero('BigPrR')
c
        do iT=1,iTpass
          iTstart=ichunk_T*(iT-1)+1
          iTstop =ichunk_T*iT
          if (iTstop.gt.npairs) iTstop=npairs
          iTsize=iTstop-iTstart+1
          call store(iTstart,iTstop,idimen,nval,iTstore)
          call matconn('BigT','r',idimen2,iTsize,iTstore)
c
          transposition: do itwo=1,2
          if (itwo.eq.1) then 
c Without transposition:
            transposed=.false.
          else
c With transposition:
          transposed=.true.
          call BigTransposition(iTstart,iTstop,idimen,nval,iTstore)
          endif
c
          call storeAlpha(ncf,nval,iTstart,iTstop,iRstart,
     *          iRstop,bl(ialphas),transposed,cc,ccsd,nfirst,F,nmo,vorb,
     *          ialphabase,ndiskalpha,betam,reset,af,nbf)
          call matconn('alpha','r',iTsize,iRsize,ialphas)
          call matmmul2('BigT','alpha','BigPrR','n','n','A')
c===========================
          call matdisc('alpha')
        enddo transposition
        call matdisc('BigT')
        enddo
c
        call writeresult(ndisk,iRstart,iRstop,idimen,nval,iRstore,af)
        call matdisc('BigPrR')
      enddo
      call matremark
      call dynamic_retmark
      call retmark
      end
C==============================================================================
      subroutine Fileinit41a(ncf,ndisk,af)
      use ccounters
      implicit none
      integer ncf,ndisk
c
      integer lrec,len,len1,info
      character*256 scrfile,filename
      logical fileopen,af
      save filename,len1
      ic_Fileinit41a=ic_Fileinit41a+1
      lrec=ncf*ncf*8
      if (.not.af) ndisk=71
      if (ic_Fileinit41a.eq.1) then
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filename=scrfile(1:len)//'.41a'
      len1=len+4
      endif
c
      if (af) then
        if (ic_Fileinit41a.gt.1) then
          call fafClosem(ndisk,0,info)
        else
          call fafClosem(ndisk,-1,info) ! dummy close
        endif
        call fafCreatem(filename(1:len1)//char(0),ndisk,lrec)
      else
      if (ic_Fileinit41a.eq.1) then
        fileopen=.false.
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 71 was used by somebody'
      else
        CLOSE(ndisk,STATUS='DELETE')
      endif
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=lrec)
      endif
      end
C==============================================================================
      subroutine 
     *  storeAlpha(ncf,        nval,  iTstart,    iTstop, iRstart,
     *             iRstop,     resul, transposed, cc,     ccsd,
     *             nfirst,     F,     nmo,        vorb,   ialphabase,
     *             ndiskalpha, betam, reset,      af,     nbf)
c
c This simple subroutine stores alphas from kl range [iTstart,iTstop] and 
c ij range [iRstart,iRstop] into memory location resul. The kl and ij are
c considered as orbital pairs, but if transposed=.true. kl is reversed, i.e.
c k.le.l instead of normal ordering k.ge.l or i.ge.j
c Parameters:
c   ncf        - number of contracted basis functions
c   nval       - number of correlated orbitals
c   iTstart    - see descr. above
c   iTstop     - ditto
c   iRstart    - ditto
c   iRstop     - ditto
c   resul      - ditto
c   transposed - ditto
c   cc         - not used any more
c   ccsd       - ditto
c   nfirst     - ditto
c   F          - ditto
c   nmo        - number of molecular orbitals
cI  vorb       - .true.: the dimension of virtual space is equal to no of virt.
c                orbitals, and MO form of integral matrices is used, otherwise,
c                the virtual space has ncf dimension and all calculations are 
c                performed for matrices in AO form.
c   ialphabase - the location of memory cache for alphas, see descr. of alpha2
c                function, where it is used.
c   ndiskalpha - the file unit number where alphas are stored
c   betam      - the beta matrix
c   reset      - passed to alpha2, see secr. of alpha2 function
cI  af         - if .true. the job is parallel & Array Files is used as IO
c                device, otherwise job is single processor and local IO is used.
      use memory
      implicit none
      integer ncf,nval,iTstart,iTstop,iRstart,iRstop,nfirst,nmo,nbf
      integer ialphabase,ndiskalpha
      real*8 resul(iTstop-iTstart+1,*),F(ncf,ncf)
      real*8 delta,beta,alpha1,betam(nval,nval),alpha2
      integer iB,jB,iE,jE,lB,kB,lE,kE,ij,kl,ijmatrix,klmatrix
      integer i,j,k,l,imap,jmap,icurraddr,nval2,nval3,icopy,iaddress
      integer i1,i2,idimen
      logical transposed,cc,ccsd,vorb,logic,reset,af
      common /CCAlphaCom/iaddress,i1,logic
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      nval3=nval*nval*nval
      nval2=nval*nval
      iB=imap(iRstart)
      iE=imap(iRstop)
      kB=imap(iTstart)
      kE=imap(iTstop)
      ij=(iB-1)*iB/2
      ijmatrix=0
      outermost: do i=iB,iE
        do j=1,i
          ij=ij+1
          if (ij.lt.iRstart) cycle
          if (ij.gt.iRstop) exit outermost
          ijmatrix=ijmatrix+1
          klmatrix=0
          kl=(kB-1)*kB/2
          outer: do k=kB,kE
            do l=1,k
             kl=kl+1
             if (kl.lt.iTstart) cycle
             if (kl.gt.iTstop) exit outer
             klmatrix=klmatrix+1
c     print *,'Coupled: ',ccsd
             if (transposed) then
             resul(klmatrix,ijmatrix)=alpha2(l,k,i,j,iRstart,iRstop,
     *                              ialphabase,ndiskalpha,nval,reset,af)
             if (j.eq.k) resul(klmatrix,ijmatrix)=
     *                   resul(klmatrix,ijmatrix)-betam(l,i)
             if (l.eq.i) resul(klmatrix,ijmatrix)=
     *                   resul(klmatrix,ijmatrix)-betam(k,j)
             else
             resul(klmatrix,ijmatrix)=alpha2(k,l,i,j,iRstart,iRstop,
     *                              ialphabase,ndiskalpha,nval,reset,af)
             if (j.eq.l) resul(klmatrix,ijmatrix)=
     *                   resul(klmatrix,ijmatrix)-betam(k,i)
             if (k.eq.i) resul(klmatrix,ijmatrix)=
     *                   resul(klmatrix,ijmatrix)-betam(l,j)
             endif
             if (k.eq.l) 
     *          resul(klmatrix,ijmatrix)=resul(klmatrix,ijmatrix)*0.5d0
            enddo
          enddo outer
        enddo
      enddo outermost
      end
C==============================================================================
      subroutine writeresult(ndisk,        ijstart, ijstop, ncf, nval,
     *                       ibaseaddress, af)
c This general subroutine writes matrices numbered by pair index 'ij',
c ranging from 'ijstart' to 'ijstop', located in bl(ibaseaddress) into file 
c unit number 'ndisk'. Each matrix constitutes separate record and is written 
c into record number 'ij'.
c Parameters:
c   ndisk        - see descr.
c   ijstart      - see descr.
c   ijstop       - see descr.
c   ncf          - number of contracted basis functions
c   nval         - number of correlated pairs
c   ibaseaddress - see descr.
cI  af          - if .true. the job is parallel & Array Files is used as IO
c                 device, otherwise job is single processor and local IO is
c                 used.
      use memory
      implicit none
      integer ijstart,ijstop,ncf,ibaseaddress,nval
      logical af
      integer i,j,ij,iaddr,icoeff,ndisk,info
      ij=0
      outer: do i=1,nval
        do j=1,i
          ij=ij+1
          if (ij.lt.ijstart) cycle
          if (ij.gt.ijstop) exit outer
          iaddr=ibaseaddress+(ij-ijstart)*ncf*ncf
          if (af) then
c           print *, 'ndisk,ncf,ij', ndisk,ncf,ij
c           call flush()
            call fafWrite(ndisk,bl(iaddr),8,ncf*ncf,1,ij,info)
          else
            call writer(ndisk,ncf,ij,bl(iaddr))
          endif
        enddo
      enddo outer
      end
C==============================================================================
      subroutine 
     *   BigTransposition(ijstart,ijstop,idimen,nval,ibaseaddress)
c This subroutine transposes all matrices located linearly in bl(ibaseaddress)
c ranging from ijstart to ijstop.
c Parameters:
c   ijstart      - see descr.
c   ijstop       - see descr.
c   idimen       - the dimension of square matrices
c   nval         - number of correlated orbitals
c   ibaseaddress - see descr.
      use memory
      implicit none
      integer ijstart,ijstop,idimen,ibaseaddress,nval
      integer iamount,iaddr,i
      iamount=ijstop-ijstart+1
      do i=1,iamount
        iaddr=ibaseaddress+(i-1)*idimen*idimen
        call matconn('Coetrans','q',idimen,idimen,iaddr)
        call matpose('Coetrans')
        call matdisc('Coetrans')
      enddo
      end
C==============================================================================
      subroutine ReadCalcG_41a(i,    j,     norecalc, ncf, nmo,
     *                         vorb, ndisk, imat,     af,  nbf)
c This subroutine reads results deposited on disk by Builder41a_G. It is 
c simply added to residuum matrix
c Parameters:
c   i        - first index of matrix to be read
c   j        - second one
c   norecalc - legacy, obsolete, should be removed.
c   ncf      - number od contracted functions
c   nmo      - number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
c   ndisk    - file logical unit number
c   imat     - the bl location where result is deposited
cI  af       - if .true. the job is parallel & Array Files is used as IO
c              device, otherwise job is single processor and local IO is
c              used.
      use memory
      implicit none
      integer i,j,ncf,nmo,ndisk,imat,nbf
      logical norecalc,vorb,af
c
      integer ij,idimen,info
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      if (norecalc) then
c Here are G and 41a
        ij=i*(i-1)/2+j
        if (af) then
           call fafRead(ndisk,bl(imat),8,idimen*idimen,1,ij,info)
           if (info.lt.8*idimen*idimen) then
           ERROR 
           endif
        else
          call reader(ndisk,idimen,ij,bl(imat))
        endif
      else
              STOP 'vbnreutlskdjfjslsldjkfhds'
      endif
      end
c==============================================================================
      subroutine CCYZ(kstart, kstop,     jstart,   jstop,     ncf,
     *                nval,   irecadrx,  irecadrc, npairs,    ndiskx,
     *                ndiskc, lbinx,     lbinc,    irecadrtx, ndisktx,
     *                lbintx, irecadrtc, ndisktc,  lbintc,    thresh,
     *                byt8,   ioverlap,  cc,       ccsd,      list,
     *                nmo,    vorb,      ndiskYZ,  af,        elaps,
     *                mp4,    nrep,      otable,   vtable,    mtable,
     *                nbf)
c This subroutine calculates the Y(k,l) and Z(k,l) matrices as defined by
c formulas (47) and (48) J. Chem. Phys. 104, 6286 (1996). The calculation is
c performed in blocks, i.e. all matrices where k is in interval [kstart,kstop]
c and j in [jstart,jstop] are calculated simultaneuosly. It results in
c huge IO savings, because each read matrix can be used many times.
c   kstart    - see descr.
c   kstop     - see descr.
c   jstart    - see descr.
c   jstop     - see descr.
c   ncf       - number of contracted basis functions
c   nval      - number of correlated orbitals
c   irecadrx  - the record info for file containing Kext bins. Note, that 
c               this array has different size and different information is
c               contained in it in single processor mode and in parallel mode.
c   irecadrc  - the record info for file containing coulomb bins. See descr for
c               irecadrx
c   npairs    - number of correlated pairs
c   ndiskx    - file unit number, where Kext bins are kept
c   ndiskc    - ditto for coulomb
c   lbinx     - bin length for Kext
c   lbinc     - ditto for coulomb
c   irecadrtx - the record info for exchange TEIO, see irecadrx also
c   ndisktx   - file unit number, where exchange TEIO are kept
c   lbintx    - bin length for exchange TEIO
c   irecadrtc - ditto for coulomb TEIO
c   ndisktc   - ditto for coulomb TEIO
c   lbintc    - ditto for coulomb TEIO
c   thresh    - the integral thrshold
c   byt8      - .true.: 
c   ioverlap  - the overlap matrix
c   cc        - if .true. all CCD contributions are calculated, if singles
c               is .true. also QCISD singles are calculated
c   ccsd      - if .true. full CCSD is calculated. It amplies also the cc=.true.
c   list      - the list of strong, weak and distand pairs.
c   nmo       - total number of occupied HF orbitals
cI  vorb      - .true.: the dimension of virtual space is equal to no of virt.
c               orbitals, and MO form of integral matrices is used, otherwise,
c               the virtual space has ncf dimension and all calculations are 
c               performed for matrices in AO form.
c   ndiskYZ   - the file logican number where Y & Z matrices will be kept.
cI  af        - if .true. the job is parallel & Array Files is used as IO
c               device, otherwise job is single processor and local IO is
c               used.
c   elaps     - the sotage for different timings
c   mp4       - .true. - MP4 energy is calculated.
      use memory
      implicit none
      integer kstart,kstop,jstart,jstop,ncf,nval,irecadrx(*),npairs
      integer ndiskx,lbinx,nbf
      integer              irecadrc(*),       ndiskc,lbinc
      integer             irecadrtx(*),      ndisktx,lbintx
      integer             irecadrtc(*),      ndisktc,lbintc
      integer ioverlap,list(*),nmo,nrep,otable(*),vtable(*),mtable(8,8)
      real*8 thresh,elaps(5) ! read, calculate, write
      logical byt8,cc,ccsd,vorb,af,mp4
c
      logical fileopen
      integer ndiskYZ,irecord,mataddr,l
      integer ireqrecord,iwork33,icoeffAO
      integer ilvec,icvec,kl,kj,lenrec
      integer iGresult,info
c fragile!
      integer nrow,ncol,ncf1,idimen
      integer idimen1,iii,iY,iZ
c New blocking data:
      integer ksize,jsize,iYZmemory,kl_mem,lj_mem,idimen2,iYZtot
      integer kl_base,lj_base,k,j,kl_addr,lj_addr
      integer kprim,jprim,ns,kjprim,nvirt,irep1,irep2,klrep,ljrep
      logical signum
      logical calc_diff
      real*8 threshcoef,elaps0,elaps1,elaps00,elaps11
c
c Elaps: 1 reading for blocking
c        2 calculations, includes 4
c        3 writing
c        4 reading singles integrals kijk
c        5 reading out of blocking
c
      ncf1=ncf
      calc_diff=.false.
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      idimen1=idimen
      idimen2=idimen*idimen
      nvirt=nbf-nmo
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('work33','q',idimen,idimen)
      iwork33=mataddr('work33')
      call matconn('YZ_overl','q',ncf,ncf,ioverlap)
      if (ccsd) then
        call dynamic_matdef('lvec','r',idimen,1)
        ilvec=mataddr('lvec')
        call dynamic_matdef('cvec','r',idimen,1)
        icvec=mataddr('cvec')
      endif
c
c Reservation of memory for Y matrices:
c The amount: (kstop-kstart+1)*(jstop-jstart+1)*(idimen*idimen)
      ksize=kstop-kstart+1
      jsize=jstop-jstart+1
      iYZmemory=ksize*jsize*(idimen*idimen)
      kl_mem=ksize*(idimen*idimen)
      lj_mem=jsize*(idimen*idimen)
c     call dynamic_getmem(iYZmemory,iYZtot)
      call dynamic_getmem(iYZmemory,iYZtot)
      call zeroit(bl(iYZtot),iYZmemory)
c     call dynamic_getmem(kl_mem,kl_base)
      call dynamic_getmem(kl_mem,kl_base)
c     call dynamic_getmem(lj_mem,lj_base)
      call dynamic_getmem(lj_mem,lj_base)
c
      if (cc.or.ccsd.or.mp4) then
      do l=1,nval
c Read matrices needed for k and j:
        call elapsec(elaps0)
        do k=kstart,kstop
          kl_addr=kl_base+(k-kstart)*idimen2
          call ExtrOne(k,    l,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(iwork33))
          call CCLoperator(bl(iwork33),idimen,bl(kl_addr))
        enddo
c
        do j=jstart,jstop
          lj_addr=lj_base+(j-jstart)*idimen2
          call CoefRead('tb',l,j,icoeffAO)
          call CCLoperator(bl(icoeffAO),idimen,bl(lj_addr))
        enddo
        call elapsec(elaps1)
        elaps(1)=elaps(1)+elaps1-elaps0
c
        call elapsec(elaps0)
        if (nrep.gt.1) irep2=otable(l)
        do k=kstart,kstop
          if (nrep.gt.1) then
            irep1=otable(k)
            klrep=mtable(irep1,irep2)
            if (klrep.lt.1.or.klrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
          endif
          kl_addr=kl_base+(k-kstart)*idimen2
          call matconn('left','q',idimen,idimen,kl_addr)
          do j=jstart,jstop
            if (nrep.gt.1) then
              irep1=otable(j)
              ljrep=mtable(irep1,irep2)
              if (ljrep.lt.1.or.ljrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
            endif
            call array_files(24)
            call full_pair_searcher(k,  j,     kprim, jprim, kjprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
            lj_addr=lj_base+(j-jstart)*idimen2
            call matconn('right','q',idimen,idimen,lj_addr)
            iY=iYZtot+((j-jstart)+(k-kstart)*jsize)*idimen2
            call matconn('Ymatrix','q',idimen,idimen,iY)
c           call matmmul2('left','right','Ymatrix','n','n','a')
          call canonical_symm_mult(bl(kl_addr),bl(lj_addr),bl(iY),
     *                             idimen, idimen,idimen,klrep,ljrep,
     *                             nvirt,nrep,vtable, mtable, vorb, 'a')
            if (ccsd) then
              call elapsec(elaps00)
              call xlijk(k,l,j,bl(ilvec))
              call elapsec(elaps11)
              elaps(4)=elaps(4)+elaps11-elaps00
              call readsingles(l,bl(icvec))
              call matscal('lvec',-1.0d0)
              call matmmul2('lvec','cvec','Ymatrix','n','t','a')
            endif
            call matdisc('Ymatrix')
            call matdisc('right')
          enddo ! do j=jstart,jstop
          call matdisc('left')
        enddo   ! do k=kstart,kstop
        call elapsec(elaps1)
        elaps(2)=elaps(2)+elaps1-elaps0
      enddo     ! do l=1,nval
      endif     ! if (cc.or.ccsd.or.mp4)
c
      do k=kstart,kstop
        do j=jstart,jstop
          call full_pair_searcher(k,  j,     kprim, jprim, kjprim,
     *                            ns, signum)
          if (ns.ne.0) cycle
        call elapsec(elaps0)
        iY=iYZtot+((j-jstart)+(k-kstart)*jsize)*idimen2
        call matconn('Ymatrix','q',idimen,idimen,iY)
        if (cc.or.ccsd.or.mp4) then
          call matscal('Ymatrix',0.5d0)
          if (ccsd) then
            call Fock_Vector_Read(k,bl(ilvec))
            call readsingles(j,bl(icvec))
            call matmmul2('lvec','cvec','Ymatrix','n','t','a')
          endif
          if (.not.vorb) then
          call matmmult('Ymatrix','YZ_overl','work33')
          call matcopy('work33','Ymatrix')
          endif
        endif
        call elapsec(elaps1)
        elaps(2)=elaps(2)+elaps1-elaps0
c
        call elapsec(elaps0)
        if (.not.mp4) then
        call ExtrOne(k,    j,     irecadrx,   npairs,    ndiskx,
     *               ncf,  lbinx, thresh,     byt8,      'x',
     *               'mo', nmo,   vorb,       bl(iwork33))
        call matadd('work33','Ymatrix')
        call ExtrOne(k,    j,     irecadrc,   npairs,    ndiskc,
     *               ncf,  lbinc, thresh,     byt8,      'c',
     *               'mo', nmo,   vorb,       bl(iwork33))
        call matadd1('work33',-0.5d0,'Ymatrix')
        if (ccsd) then
          call ExtrOne(k,    j,     irecadrtx,  npairs,    ndisktx, 
     *                 ncf,  lbintx,thresh,     byt8,      'tx',
     *                 'mo', nmo,   vorb,       bl(iwork33))
          call matadd('work33','Ymatrix')
          call ExtrOne(k,    j,     irecadrtc,  npairs,    ndisktc,
     *                 ncf,  lbintc,thresh,     byt8,      'tc',
     *                 'mo', nmo,   vorb,       bl(iwork33))
c        call CoulombBuilder(k,j,idimen,npairs,ncf,nmo,vorb,bl(iwork33))
          call matadd1('work33',-0.5d0,'Ymatrix')
        endif
      endif ! if (.not.mp4)
        call elapsec(elaps1)
        elaps(5)=elaps(5)+elaps1-elaps0
        ireqrecord=((k-1)*nval+j)*2
        call elapsec(elaps0)
        if (af) then
          call fafwrite(ndiskYZ,bl(iY),8,idimen*idimen,1,ireqrecord-1,
     *                  info)
        else
          call writer(ndiskYZ,idimen,ireqrecord-1,bl(iY))
        endif
        call elapsec(elaps1)
        elaps(3)=elaps(3)+elaps1-elaps0
        call matdisc('Ymatrix')
        enddo
      enddo
c
      call zeroit(bl(iYZtot),iYZmemory)
c  Z matrix
c
      if (cc.or.ccsd.or.mp4) then
c
      do l=1,nval
c Read matrices needed for k and j:
        call elapsec(elaps0)
        do k=kstart,kstop
          kl_addr=kl_base+(k-kstart)*idimen2
          call ExtrOne(l,    k,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'mo', nmo,   vorb,       bl(kl_addr))
        enddo
c
        do j=jstart,jstop
          lj_addr=lj_base+(j-jstart)*idimen2
          call CoefRead('tb',j,l,icoeffAO)
          call tfer(bl(icoeffAO),bl(lj_addr),idimen2)
        enddo
        call elapsec(elaps1)
        elaps(1)=elaps(1)+elaps1-elaps0
c
        call elapsec(elaps0)
        if (nrep.gt.1) irep2=otable(l)
        do k=kstart,kstop
          if (nrep.gt.1) then
            irep1=otable(k)
            klrep=mtable(irep1,irep2)
            if (klrep.lt.1.or.klrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
          endif
          kl_addr=kl_base+(k-kstart)*idimen2
          call matconn('left','q',idimen,idimen,kl_addr)
          do j=jstart,jstop
            if (nrep.gt.1) then
              irep1=otable(j)
              ljrep=mtable(irep1,irep2)
              if (ljrep.lt.1.or.ljrep.gt.8) STOP 'sdf9874l3th3t32l5iuth'
            endif
            call array_files(25)
            call full_pair_searcher(k,  j,     kprim, jprim, kjprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
            lj_addr=lj_base+(j-jstart)*idimen2
            call matconn('right','q',idimen,idimen,lj_addr)
            iZ=iYZtot+((j-jstart)+(k-kstart)*jsize)*idimen2
            call matconn('Zmatrix','q',idimen,idimen,iZ)
c
c           call matmmul2('left','right','Zmatrix','n','n','a')
          call canonical_symm_mult(bl(kl_addr),bl(lj_addr),bl(iZ),
     *                             idimen, idimen,idimen,klrep,ljrep,
     *                             nvirt,nrep,vtable, mtable, vorb, 'a')
            if (ccsd) then
              call elapsec(elaps00)
              call xkijk(l,k,j,bl(ilvec))
              call elapsec(elaps11)
              elaps(4)=elaps(4)+elaps11-elaps00
              call readsingles(l,bl(icvec))
              call matmmul2('lvec','cvec','Zmatrix','n','t','a')
            endif
            call matdisc('Zmatrix')
            call matdisc('right')
          enddo
          call matdisc('left')
        enddo
        call elapsec(elaps1)
        elaps(2)=elaps(2)+elaps1-elaps0
      enddo
      endif ! if (cc.or.ccsd)
c
      do k=kstart,kstop
        do j=jstart,jstop
            call full_pair_searcher(k,  j,     kprim, jprim, kjprim,
     *                              ns, signum)
            if (ns.ne.0) cycle
          call elapsec(elaps0)
          iZ=iYZtot+((j-jstart)+(k-kstart)*jsize)*idimen2
          call matconn('Zmatrix','q',idimen,idimen,iZ)
          if (ccsd.or.cc.or.mp4) then
            if (.not.vorb) then
              call matmmult('Zmatrix','YZ_overl','work33')
              call matcopy('work33','Zmatrix')
            endif
            call matscal('Zmatrix',-1.0d0)
          endif
          call elapsec(elaps1)
          elaps(2)=elaps(2)+elaps1-elaps0
          call elapsec(elaps0)
          if (ccsd) then
            call ExtrOne(k,    j,      irecadrtc,  npairs,    ndisktc,
     *                   ncf,  lbintc, thresh,     byt8,      'tc',
     *                   'mo', nmo,    vorb,       bl(iwork33))
            call matadd('work33','Zmatrix')
          endif
          if (.not.mp4) then
          call ExtrOne(k,    j,     irecadrc,   npairs,    ndiskc,
     *                 ncf,  lbinc, thresh,     byt8,      'c',
     *                 'mo', nmo,   vorb,       bl(iwork33))
          call matadd('work33','Zmatrix')
          endif ! if (.not.mp4)
          call elapsec(elaps1)
          elaps(5)=elaps(5)+elaps1-elaps0
          ireqrecord=((k-1)*nval+j)*2
          call elapsec(elaps0)
          if (af) then
            call fafwrite(ndiskYZ,bl(iZ),8,idimen*idimen,1,ireqrecord,
     *                    info)
          else
            call writer(ndiskYZ,idimen,ireqrecord,bl(iZ))
          endif
          call elapsec(elaps1)
          elaps(3)=elaps(3)+elaps1-elaps0
          call matdisc('Zmatrix')
        enddo
      enddo
c
      call dynamic_retmem(3)
c
      call matremark
      call retmark
      call dynamic_retmark
      end
C==============================================================================
      subroutine buildX(irecadrtx, npairs,    ndisktx, lbintx, thresh,
     *                  byt8,      irecadrtc, ndisktc, lbintc, ioverl,
     *                  nval,      ncf,       cc,      ccsd,   ifockAO,
     *                  iX,        nmo,       vorb,    mp4,    nbf)
c Thissubroutine builds X matrix, as defined by formula no. (46) 
c J. Chem. Phys. 104, 6286 (1996)
c Parameters:
c   irecadrx  - the record info for file containing TEIO exchange bins. Note,
c               that this array has different size and different information
c               is contained in it in single processor mode and in parallel
c               mode.
c   npairs    - number of correlated pairs
c   ndisktx   - the file logical unit number where TEIO exchange integrals are
c               kept
c   lbintx    - TEIO exchange integrals bin length
c   thresh    - integral threashold
c   byt8      - .true.: double precision bin are used, otherwise int4+int1
c   irecadrtc - like irecadrx for coulomb TEIO
c   ndisktc   - like ndisktx for coulomb TEIO
c   lbintc    - like lbintx for coulomb TEIO
c   ioverl    - the overlap matrix
c   nval      - number of correlated orbitals
c   ncf       - number of contracted basis functions
c   cc        - .true. = CCD contributions are calculated
c   ccsd      - .true. = full CCSD is calculated, also cc has to be true
c   ifockAO   - the Fock matrix, either n AO or in MO form, MO form 
c               contains only virtual part.
c   iX        - the resulting matrix
c   nmo       - total number of occupied HF orbitals
cI  vorb      - .true.: the dimension of virtual space is equal to no of virt.
c               orbitals, and MO form of integral matrices is used, otherwise,
c               the virtual space has ncf dimension and all calculations are 
c               performed for matrices in AO form.
c   mp4       - .true. - MP4 energy is calculated.
      use memory
      implicit none
      integer npairs,ndisktx,lbintx,ndisktc,lbintc,ioverl,nval,ncf
      integer ifockAO,iX,nmo,nbf
      integer irecadrtx(*),irecadrtc(*)
      logical byt8,cc,ccsd,vorb,mp4
      real*8 thresh
c
c
      integer iA,mataddr,iwork11,iwork22,k,iRvec,isinglv,idimen,kcp
c
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
      call dynamic_matdef('A','q',idimen,idimen)
      iA=mataddr('A')
      call matconn('FockAO','q',idimen,idimen,ifockAO)
      call matconn('output','q',idimen,idimen,iX)
      if (mp4) then
        call matzero('output')
      else
        call matcopy('FockAO','output')
      endif
c     print *,'Coupled: ',ccsd
      if (ccsd) then
        call dynamic_matdef('work11','q',idimen,idimen)
        iwork11=mataddr('work11')
        call dynamic_matdef('work22','q',idimen,idimen)
        iwork22=mataddr('work22')
        call dynamic_matdef('Rvec','r',idimen,1)
        iRvec=mataddr('Rvec')
        call dynamic_matdef('singlv','r',idimen,1)
        isinglv=mataddr('singlv')
        call matconn('overl','q',idimen,idimen,ioverl)
      endif
c     print *,'Coupled: ',ccsd
      if (cc.or.ccsd.or.mp4) then
      call CCA(ncf,iA,nmo,vorb,nbf) !OK?
      call matadd1('A',-1.0d0,'output')
      endif
c     print *,'Coupled: ',ccsd
      if (ccsd) then
        do k=1,nval
           kcp=k
          call ExtrOne(k,    kcp,   irecadrtx,  npairs,    ndisktx,
     *                 ncf,  lbintx,thresh,     byt8,      'tx',
     *                 'mo', nmo,   vorb,       bl(iwork11))
          call matadd1('work11',-1.0d0,'output')
          call ExtrOne(k,    kcp,   irecadrtc,  npairs,    ndisktc,
     *                 ncf,  lbintc,thresh,     byt8,      'tc',
     *                 'mo', nmo,   vorb,       bl(iwork11))
          call matadd1('work11',2.0d0,'output')
        enddo
        call matzero('work11')
        do k=1,nval
          call readR(k,bl(iRvec))
          call readsingles(k,bl(isinglv))
          call matmmul2('Rvec','singlv','work11','n','t','a')
        enddo
c       call general_write1(99,1,8)
c       call general_write1(99,bl(iwork11),idimen*idimen*8)
c       call flush(99)
        if (vorb) then
          call matadd1('work11',-1.0d0,'output')
        else
          call matmmult('work11','overl','work22')
          call matadd1('work22',-1.0d0,'output')
        endif
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
C==============================================================================
      subroutine CCYZ_read(k,j,ncf,nval,imat,nmo,vorb,ndiskYZ,af,reset,
     *                     typ,nbf)
c This simple subroutine reads Y and Z matrices constructed by CCYZ, which are
c used in Qparts.
c Parameters:
c   k       - the k index (first index) of Y or Z matrix
c   j       - the i index (second) of Y or Z matrix
c   ncf     - number of contracted basis functions
c   nval    - number of correlated orbitals
c   imat    - the read matrix is deposited here
c   nmo     - total number of occupied HF orbitals
cI  vorb    - .true.: the dimension of virtual space is equal to no of virt.
c             orbitals, and MO form of integral matrices is used, otherwise,
c             the virtual space has ncf dimension and all calculations are 
c             performed for matrices in AO form.
c   ndiskYZ - the file unit number, where Y & Z matrices are stored.
cI  af      - if .true. the job is parallel & Array Files is used as IO
c             device, otherwise job is single processor and local IO is
c             used.
c   reset   - this variable is set to .true. at the beginning of iteration, it
c             informs this subroutine that old cache files should be deleted,
c             and cache data resetted. Useful only with 'cache' option provided
c             via common block /GlobalCCSD/
c   typ     - may be 'Y' or 'Z', it determined which matrix is needed.
      use memory
      use ccounters
      implicit none
      integer k,j,ncf,nval,imat,nmo,ndiskYZ,nbf
      logical vorb,af,reset
c
      integer idimen,ireqrecord,ireqrec,info,ndiskYZ_cache,len
      integer len1,lrec,ii,mygid,irealrec
      integer kprim,jprim,kjprim,ns,iscr
      logical signum
c FRAGILE!
      integer*1 irec_tab(PAR_MAX)
      character*256 scrfile,filename
      character*2 ch_mygid
      character*1 typ ! may be 'Y' or 'Z'
      save irec_tab,ndiskYZ_cache
c
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
c
      ic_CCYZ_read=ic_CCYZ_read+1
c
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      if (cache.and.af.and.ic_CCYZ_read.eq.1) then
        call getival('mgid',mygid)
        write(ch_mygid,'(I2.2)') mygid
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len)
        lrec=idimen*idimen*8
        filename=scrfile(1:len)//'.XYij--'//ch_mygid
        len1=len+9
        ndiskYZ_cache=70
      OPEN (UNIT=ndiskYZ_cache,FILE=filename(1:len1),FORM='UNFORMATTED',
     1      ACCESS='DIRECT',RECL=lrec)
      endif
c
      if (cache.and.reset) then
        do ii=1, (PAR_MAX)
         irec_tab(ii)=0
        enddo
        reset=.false.
      endif
c
      call full_pair_searcher(k,  j,     kprim, jprim, kjprim,
     *                              ns, signum)
      ireqrec=kjprim
      ireqrecord=ireqrec*2
      if (typ.eq.'Y') irealrec=ireqrecord-1
      if (typ.eq.'Z') irealrec=ireqrecord
c
      if (af) then
      if ((.not.cache).or.irec_tab(irealrec).eq.0) then
      call fafRead(ndiskYZ,bl(imat),8,idimen*idimen,1,irealrec,info)
          if (info.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
      if (cache) call writer(ndiskYZ_cache,idimen,irealrec,bl(imat))
      if (cache) irec_tab(irealrec)=1
      else
      call reader(ndiskYZ_cache,idimen,irealrec,bl(imat))
      endif
      else
        call reader(ndiskYZ,idimen,irealrec,bl(imat))
      endif
      if (ns.ne.0) then
        call full_symm_transform(idimen, ns, signum,bl(imat))
      endif
      end
C==============================================================================
      subroutine prepare_CCYZ_file(ncf,nmo,vorb,ndiskYZ,af,nbf)
c This subroutine figures out what is common "core" name of scratch files
c and creates file for Y & Z matrices.
c Parameters:
c   ncf     - number of contracted basis functions
c   nmo     - total number of occupied HF orbitals
cI  vorb    - .true.: the dimension of virtual space is equal to no of virt.
c             orbitals, and MO form of integral matrices is used, otherwise,
c             the virtual space has ncf dimension and all calculations are 
c             performed for matrices in AO form.
c   ndiskYZ - the file unit number of the created file
cI  af      - if .true. the job is parallel & Array Files is used as IO
c             device, otherwise job is single processor and local IO is
c             used.
      use ccounters
      implicit none
      integer ncf,nmo,ndiskYZ,nbf
      logical vorb,af
c
      character*256 filename,scrfile
      integer ndisk_Y,ndisk_Z,len,len1,lenrec,idimen,irecord
      integer info
      logical calc_diff,fileopen
      save filename,len1
c
      ic_CCYZ_file=ic_CCYZ_file+1
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      calc_diff=.false.
c WARNING !!!! Not changed for AF!
      if (ic_CCYZ_file.eq.1.and.calc_diff) then
        ndisk_Y=81
        ndisk_Z=82
        inquire(ndisk_Y,OPENED=fileopen)
        if (fileopen) STOP 'File 81 opened in CCYZ'
        inquire(ndisk_Z,OPENED=fileopen)
        if (fileopen) STOP 'File 82 opened in CCYZ'
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                 ! 80 is starting length
        len1=len+8
        lenrec=idimen*idimen*8
        filename=scrfile(1:len)//'.CCY_res'
        OPEN(UNIT=ndisk_Y,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lenrec)
        filename=scrfile(1:len)//'.CCZ_res'
        OPEN(UNIT=ndisk_Z,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lenrec)
      endif
c
      if (.not.af) ndiskYZ=70
      irecord=idimen*idimen*8
      if (ic_CCYZ_file.eq.1) then
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filename=scrfile(1:len)//'.YZ'
      len1=len+3
      endif
      if (af) then
        if (ic_CCYZ_file.gt.1) then
          call fafClosem(ndiskYZ,0,info)
        else
          call fafClosem(ndiskYZ,-1,info)
        endif
      else
        if (ic_CCYZ_file.eq.1) then
          inquire(ndiskYZ,opened=fileopen)
          if (fileopen) STOP 'Unit 70 was used by somebody'
        else
          CLOSE(ndiskYZ,STATUS='DELETE')
        endif
      endif
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndiskYZ,irecord)
      else
      OPEN(UNIT=ndiskYZ,FILE=filename(1:len1),FORM='UNFORMATTED',
     *     ACCESS='DIRECT',RECL=irecord)
      endif
      end
C==============================================================================
      subroutine prepare_diis_files(ndiskr,ndiskc,af,ncf,nmo,vorb,nbf)
c The two files unsed by DIIS subroutine are prepared here: one stores
c residua and second coefficients. It might be optimized in the future
c because I think it is possible not to store coefficients with the same
c result.
c Parameters:
c   ndiskr - the file with residua
c   ndiskc - the file with coefficients (might be removed in the future)
cI  af     - if .true. the job is parallel & Array Files is used as IO
c            device, otherwise job is single processor and local IO is
c            used.
c   ncf    - number of contracted functions
c   nmo    - total number of occupied HF orbitals
cI  vorb   - .true.: the dimension of virtual space is equal to no of virt.
c            orbitals, and MO form of integral matrices is used, otherwise,
c            the virtual space has ncf dimension and all calculations are 
c            performed for matrices in AO form.
      implicit none
      integer ndiskr,ndiskc,ncf,nmo,nbf
      logical af,vorb
c
      integer irecord,ivirt,len1,len,idimen
      character*256 scrfile,filename
      logical fileopen
c
      ivirt=nbf-nmo
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      if (.not.af) then
        ndiskr=74
        ndiskc=75
        irecord=ivirt*ivirt*8
        filename=scrfile(1:len)//'.resid'
        len1=len+6
        inquire(ndiskr,opened=fileopen)
        if (fileopen) STOP 'Unit 74 was used by somebody'
        OPEN(UNIT=ndiskr,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
c
        filename=scrfile(1:len)//'.diis'
        len1=len+5
        irecord=idimen*idimen*8
        inquire(ndiskc,opened=fileopen)
        if (fileopen) STOP 'Unit 75 was used by somebody'
        OPEN(UNIT=ndiskc,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      else
        filename=scrfile(1:len)//'.resid'
        len1=len+6
        call fafCreatem(filename(1:len1)//char(0),ndiskr,ivirt*ivirt*8)
        filename=scrfile(1:len)//'.diis'
        len1=len+5
       call fafCreatem(filename(1:len1)//char(0),ndiskc,idimen*idimen*8)
      endif
      end
C==============================================================================
      subroutine calculate_ij_chunk(ichunk,jchunk,nslv,idimen,nval)
c This subroutine calculates the block size for Qparts and  CCYZ sybroutines.
c the vlock size has two dimensions ichunk,jchunk.
c Parameters:
c   ichunk - the "X" dimension of block size
c   jchunk - the "Y" dimension of block size
c   nslv   - number of slaves
c   idimen - the dimension of matrices used for calculation, depends on
c            the basis set size and the "mode" of calculation, it might be 
c            AO mode or MO mode.
c   nval   - number of correlated orbitals
      implicit none
c
      integer ichunk,jchunk,nslv,idimen,nval
c
      integer iad,lcore,memavail,max_matrix,max_matrix1,max_ichunk
      integer max_jchunk
c
c
      call dynamic_show_free(memavail)
c     max cache:
      max_matrix=memavail/(idimen*idimen)
      ichunk=int(sqrt(dble(max_matrix)))
      jchunk=int(sqrt(dble(max_matrix)))
      do
        if ((ichunk*jchunk+ichunk+jchunk).le.max_matrix) exit
        ichunk=ichunk-1
        jchunk=jchunk-1
      enddo
c I am going to find such distribution, that I will have as small amount
c of blocks as possible. But, I shouldn't have one greater than nslv*n for
c example. The number of iddle slaves should not be greater than say,
c 20% (1/5). If such division will be possible, no reducing of sizes
c will be performed.
      if (nslv.gt.1) then
        max_ichunk=int( sqrt( dble(nval*nval/(5*nslv)) ) )
        max_jchunk=int( sqrt( dble(nval*nval/(5*nslv)) ) )
      else
        max_ichunk=nval
        max_jchunk=nval
      endif
      if (ichunk.gt.max_ichunk) ichunk=max_ichunk
      if (jchunk.gt.max_jchunk) jchunk=max_jchunk
      if (ichunk.lt.1) ichunk=1
      if (jchunk.lt.1) jchunk=1
c     if (ichunk.le.2.or.jchunk.le.2)
c    *  print *,'WARNING! chunk small, less or equal to 2!'
      end
!=======================================================================
      subroutine CreateAlphaDisk(nval,ndisk,af)
c the file for alphas is created here
c Parameters:
c   nval  - number of correlated orbitals
c   ndisk - the returned file unit number
cI  af    - if .true. the job is parallel & Array Files is used as IO
c           device, otherwise job is single processor and local IO is
c           used.
      implicit none
      integer nval,ndisk
      logical af
c
      integer len,len1,irecord
      logical fileopen
      character*256 scrfile,filename
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                  80 is starting length
      filename=scrfile(1:len)//'.alphas'
      len1=len+7
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk,nval*nval*8)
      else
        ndisk=67
        irecord=nval*nval*8
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 67 was used by somebody'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
      end
!=======================================================================
      subroutine 
     * SupplementAlpha_BetaBuild(nval,      ndiskalpha,beta,ncf, nfirst,
     *                           nmo,       vorb,      cc,  ccsd,F,
     *                           kijklndisk,af,        mp4, nbf)
c This subroutine reads all "pure" CCD alphas calculated before, where alpha is 
c defined as in J. Chem. Phys. 81, 1901 (1984), and singles contributions
c are added as needed for CCSD, also the K(ijkl) integrals over occupied 
c correlated MO orbitals are added here. The beta scalar coefficients, 
c as defined by J. Chem. Phys. 104, 6286 (1996) are also built here in matrix
c form. Betas are constructed partly from alphas and Fock matrix elements
c are added, also some contribution to them is calculated independetly.
c This procedure is not parallelized, because the total time of following
c operations is negligible.
c Parameters:
c   nval       - number of correlated orbitals
c   ndiskalpha - the file unit nuber, where the "pre-computed" alphas are stored
c                and where full alphas will be deposited
c   beta       - the direct memory location where beta will be stored
c   ncf        - number of contracted functions
c   nfirst     - the first correlated orbital
c   nmo        - total number of occupied HF orbitals
cI  vorb       - .true.: the dimension of virtual space is equal to no of virt.
c                orbitals, and MO form of integral matrices is used, otherwise,
c                the virtual space has ncf dimension and all calculations are 
c                performed for matrices in AO form.
c   cc         - .true.: the CCD is calculated or QCISD
c   ccsd       - .true.: CCSD
c   F          - Fock matrix in MO
c   kijklndisk - the file unit number where K(ijkl) integrals are kept
cI  af         - if .true. the job is parallel & Array Files is used as IO
c                device, otherwise job is single processor and local IO is
c                used.
c   mp4        - guess what does it mean ;)
      use memory
      implicit none
      integer nval,ndiskalpha,ncf,nfirst,nmo,kijklndisk,nbf
      real*8 beta(nval,nval),F(ncf,ncf)
      logical vorb,cc,ccsd,af,mp4
c
      integer nval2,nval3,nval4,ialphas,ialphamat,mataddr,icoefvi,ikopv
      integer iaddress,i,j,k,l,ib,kb,ij,ijrec,idimen,ifockv,icoefvv
      integer iresuu,iresu,ijaddress,istatus,ikaddr,iladdr,ijaddr
      integer iallsingl,invalmat,ivari_res,iintegral,iindex,iii
      real*8 element,beta1,xkijkl
c
      nval4=nval*nval*nval*nval
      nval3=nval*nval*nval
      nval2=nval*nval
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_getmem(nval4,ialphas)
      call dynamic_matdef('alphamat','q',nval,nval)
      ialphamat=mataddr('alphamat')
      call dynamic_matdef('allsingl','r',idimen,nval)
      iallsingl=mataddr('allsingl')
      call dynamic_matdef('nvalmat','r',nval,nval)
      invalmat=mataddr('nvalmat')
      call dynamic_matdef('integral','r',idimen,nval)
      iintegral=mataddr('integral')
      do i=1,nval
        call readsingles(i,bl(iallsingl+(i-1)*idimen))
      enddo
c Store alphas:
      ijrec=0
      if (cc.or.ccsd.or.mp4) then
      do i=1,nval
        do j=1,i
          ijrec=ijrec+1
          if (af) then
            call fafread(ndiskalpha,bl(ialphamat),8,nval*nval,1,ijrec,
     *                    istatus)
          if (istatus.lt.8*nval*nval) then
              ERROR
              call flush(6)
          endif
          else
            call reader(ndiskalpha,nval,ijrec,bl(ialphamat))
          endif
c
          ij=(i-1)*nval+j
          iaddress=ialphas+(ij-1)*nval2
          call tfer(bl(ialphamat),bl(iaddress),nval2)
c
          if (i.ne.j) then
            ij=(j-1)*nval+i
            iaddress=ialphas+(ij-1)*nval2
            call matpose('alphamat')
            call tfer(bl(ialphamat),bl(iaddress),nval2)
          endif
        enddo
      enddo
      else
        call zeroit(bl(ialphas),nval4)
      endif
c Add the stuff needed for beta
      if (ccsd) then
      call dynamic_matdef('coefvi','r',idimen,1)
      icoefvi=mataddr('coefvi')
      call dynamic_matdef('kopv','r',idimen,1)
      ikopv=mataddr('kopv')
      call dynamic_matdef('resu','q',1,1)
      do l=1,nval
        do k=1,nval
          do i=1,nval
            call xkijk(l,k,i,bl(iintegral+(i-1)*idimen))
          enddo
          call matmmul2('allsingl','integral','nvalmat','t','n','n')
          do i=1,nval
            iindex=(i-1)*nval
            do j=1,nval
               iaddress=ialphas+(i-1)*nval3+(j-1)*nval2+(k-1)*nval+(l-1)
c              call readsingles(j,bl(icoefvi))
c              call matmmul2('coefvi','kopv','resu','t','n','n')
c              call matelem('resu',1,1,element)
               bl(iaddress)=bl(iaddress)+bl(invalmat+iindex+(j-1))
            enddo
          enddo
        enddo
      enddo
      call dynamic_matrem('resu')
      call dynamic_matrem('kopv')
      call dynamic_matrem('coefvi')
      endif
c Build beta:
      if (ccsd) then
      call dynamic_matdef('fockv','r',idimen,1)
      ifockv=mataddr('fockv')
      call dynamic_matdef('coefvv','r',idimen,1)
      icoefvv=mataddr('coefvv')
      call dynamic_matdef('resuu','q',1,1)
      iresuu=mataddr('resuu')
      endif
      call dynamic_unlock(beta,iii)
      do kb=1,nval
        do ib=1,nval
        if (mp4) then
          beta1=0.0d0
        else
          beta1=F(nfirst+kb-1,nfirst+ib-1)
        endif
        if (cc.or.ccsd.or.mp4) then
          do l=1,nval
      beta1=beta1+2.0d0*bl(ialphas+nval3*(ib-1)+nval2*(l-1)+nval*(kb-1)+
     *     (l-1))-bl(ialphas+nval3*(ib-1)+nval2*(l-1)+nval*(l-1)+(kb-1))
          enddo
        endif
          if (ccsd) then
          call Fock_Vector_Read(kb,bl(ifockv))
          call readsingles(ib,bl(icoefvv))
          call matmmul2('fockv','coefvv','resuu','t','n','n')
          call matelem('resuu',1,1,element)
          beta1=beta1+element
          endif
          beta(kb,ib)=beta1
        enddo
      enddo
      call dynamic_lock(beta,iii)
      if (ccsd) then
      call dynamic_matrem('resuu')
      call dynamic_matrem('coefvv')
      call dynamic_matrem('fockv')
      endif
c Finish alphas:
      if (ccsd) then
      call dynamic_matdef('coefvi','r',idimen,1)
      icoefvi=mataddr('coefvi')
      call dynamic_matdef('kopv','r',idimen,1)
      ikopv=mataddr('kopv')
      call dynamic_matdef('resu','q',1,1)
      iresu=mataddr('resu')
      call matzero('resu')
      endif
      ijrec=0
      if (.not.mp4) then
      do i=1,nval
        do j=1,i
          ijrec=ijrec+1
          ijaddress=ialphas+(i-1)*nval3+(j-1)*nval2
          if (af) then
            call fafread(kijklndisk,bl(ialphamat),8,nval*nval,1,ijrec,
     *                   istatus)
          if (istatus.lt.8*nval*nval) then
              ERROR
              call flush(6)
          endif
          else
            call reader(kijklndisk,nval,ijrec,bl(ialphamat))
          endif
          call matconn('purealph','q',nval,nval,ijaddress)
          call matadd('alphamat','purealph')
          call matdisc('purealph')
        enddo
      enddo
      endif ! if (.not.mp4)
      if (ccsd) then
      do k=1,nval
        ikaddr=ialphas+(k-1)*nval
        do l=1,nval
          iladdr=ikaddr+(l-1)
          do j=1,nval
            call xkijk(k,l,j,bl(iintegral+(j-1)*idimen))
          enddo
            call matmmul2('allsingl','integral','nvalmat','t','n','n')
          do j=1,nval
            ijaddr=iladdr+(j-1)*nval2
            iindex=(j-1)*nval
            do i=1,nval
              iaddress=ijaddr+(i-1)*nval3
c             call readsingles(i,bl(icoefvi))
c             call matmmul2('coefvi','kopv','resu','t','n','n')
c             call matelem('resu',1,1,element)
              bl(iaddress)=bl(iaddress)+bl(invalmat+iindex+(i-1))
            enddo
          enddo
        enddo
      enddo
      endif
c
      ijrec=0
      do i=1,nval
        do j=1,i
          ijrec=ijrec+1
          ijaddress=ialphas+(i-1)*nval3+(j-1)*nval2
          if (af) then
            call fafwrite(ndiskalpha,bl(ijaddress),8,nval*nval,1,ijrec,
     *                    istatus)
          else
            call writer(ndiskalpha,nval,ijrec,bl(ijaddress))
          endif
        enddo
      enddo
      call matremark
      call dynamic_retmark
      call retmark
      end
c=======================================================================
      subroutine CreateKijklDisk(nval,ndisk,af)
c Parameters:
c   nval  - number of correlated orbitals
c   ndisk - the result : file unit number
cI  af    - if .true. the job is parallel & Array Files is used as IO
c           device, otherwise job is single processor and local IO is
c           used.
      implicit none
      integer nval,ndisk
      logical af
c
      integer len,len1,irecord
      logical fileopen
      character*256 scrfile,filename
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
c                                  80 is starting length
      filename=scrfile(1:len)//'.xkijkl'
      len1=len+7
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk,nval*nval*8)
      else
        ndisk=68
        irecord=nval*nval*8
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 68 was used by somebody'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
      end
!=======================================================================
      subroutine CCDiis(iorb,    jorb,    WhatToDo, ncf,  nmo,
     *                  npairs,  residMO, coefAO,   nval, vorba,
     *                  singles, ndiskr,  ndiskc,   mp2_cc, list,
     *                  ipairimages,nbf)
c This subroutine should be improved or, better, written from scratch.
c Right now it is too weird.
c It performs different tasks depending on WhatToDo string value. It might
c be initialization, where memory is reserved once and kept in static
c (C language term) variable, storing, where amplitudes & residua are stored
c on disk by call from updatecoef subroutine, and calculation of DIIS matrices
c at the end of iteration, associated with calculation of new, improved
c coefficient. Also it might be resetting, which means that DIIS matrix
c will be built anew.
c Parameters:
c   iorb     - in storing mode the first index of residuum matrix & ampl.
c              in 'diis' mode it means the maximum size of DIIS matrix
c   jorb     - ditto, second index
c   WhatToDo - might have the following values:
c               'reset':       - sets the local logic. variab. to specific value
c                                which causes DIIS to start building its matrix
c                                from scratch
c               'nodiis'       - set local logical variable to specific value,
c                                each next call will return immediately without
c                                any action taken
c               'diis'         - the opposite of previous, it set the size of
c                                DIIS matrix from iorb value and allocates 
c                                memory
c               'itera'        - this is message for DIIS subroutine that next
c                                iteration just started
c               'store'        - store mode, is stores residua and amplitudes
c                                on disk
c               'singlesstore' - ditto for singles ampl. and residua
c               'calculate'    - performs actual DIIS calculations, builds
c                                matrix and calculates DIIS coefficients, then
c                                new improved amplitudes
c   ncf      - number of contracted basis functions
c   nmo      - the total number of HF occupied orbitals
c   npairs   - number of correlated orbital pairs
c   residMO  - residua are passed here
c   coefAO   - coefficient matrix
c   nval     - number of correlated orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
c   singles  - if true singles are also calculated by program
c   ndiskr   - the file unit number for residua
c   ndiskc   - ditto for amplitudes
c   mp2_cc   - if 'mp2' the mp2 diis is performed (no singles), of 'ccc' then
c              coupled cluster part is performed, in 'store' more the value
c              of this variable does not matter, and this variable might be
c              set to 'xxx'
      use memory
      use ccounters
      implicit none
      external check_for_delete
      character*(*) WhatToDo,mp2_cc
      integer iorb,jorb,ncf,nmo,npairs,nval,nbf
      real*8 residMO(ncf,ncf),coefAO(ncf,ncf)
      integer ndiskr,irecord,len,len1,i,j,ivirt,ij,incf,jncf
      integer ndiskc,ivirt1,list(*),ipairimages(*)
      character*256 scrfile,filename
      logical fileopen,nodisk,singles
      logical vorba,vorb,ccsd,af
      integer iresidVOcurr,iresidVOprev
c fragile!
      real*8 diiscoeff(DIIS_LENGTH)
      integer ndiska,ndiskb,ncfa,internal_addr,ireturn,ii
      integer icallno_propagate,npairs_propagate,ndiskrco,ndisk_cache
      common /ampldisk/ndiska,ndiskb,ncfa,internal_addr,nodisk,ireturn,
     *                 icallno_propagate,npairs_propagate,ccsd,ndiskrco,
     *                 ivirt1,vorb,af
      integer icoeff,ifinalcoeff,mataddr
      integer iad2,iad3,idiisinv
      integer isiz,ical_diffe,idimen,ishift,icond,igetival,info
      integer k
      integer isize,inew,keep,icp
      logical check_for_delete
      real*8 trace,element,det,sm_element,elem
c
      icond=91
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      if (mp2_cc.ne.'mp2'.and.mp2_cc.ne.'ccc'.and.mp2_cc.ne.'xxx') 
     *    STOP 'DIIS 234 error'
c
      ivirt=nbf-nmo
      if (WhatToDo.eq.'reset') then 
      diis_restart=.true.
      return
      endif
      if (WhatToDo.eq.'nodiis') then 
       do_diis=.false.
      else if (WhatToDo.eq.'diis') then
      call dynamic_getmem(ivirt*ivirt,isiglresstore)
      call dynamic_getmem(idimen*idimen,isiglcoestore)
c Here iorb means 'max length of the DIIS expansion'.
      max_rec=iorb
      if (max_rec.gt.abs_max_rec) then
        print *, 'DIIS error, too large max_rec requested.'
        STOP 'DIIS error, too large max_rec requested'
      endif
      endif
      if (.not.do_diis) return
      if (WhatToDo.eq.'itera') then 
        ishift_record=ishift_record+1
        keep=irecord_tab(1) !           rotate START
        do i=2,max_rec
          irecord_tab(i-1)=irecord_tab(i)
        enddo
        irecord_tab(max_rec)=keep !     rotate STOP
        if (indexer.lt.max_rec) then 
          indexer=indexer+1
          irecord_tab(max_rec)=indexer
        endif
        return
      endif
c Storing part: store residuals and predicted amplitues, I think it
c should be called in UpdateCoef. Amplitudes are stored in AO format,
c residuals in MO, only ncf-nlast part (virtual part)
      if (WhatToDo.eq.'store') then
       call dynamic_getmem(ncf*ncf,iresidVOcurr)
       call dynamic_getmem(ncf*ncf,iresidVOprev)
        ishift_record=irecord_tab(max_rec)
        if (iorb.ge.jorb) then 
          ij=iorb*(iorb-1)/2+jorb
        else
          ij=jorb*(jorb-1)/2+iorb
        endif
        if (.not.vorb) then
          do j=1,ivirt
            jncf=j+nmo
            do i=1,ivirt
              incf=i+nmo
              bl(iresidVOcurr+(i-1)+(j-1)*ivirt)=residMO(incf,jncf)
            enddo
          enddo
        endif
        if (singles) then 
          ij=(npairs+1)*(ishift_record-1)+ij
        else
          ij=npairs*(ishift_record-1)+ij
        endif
        if (vorb) then
          if (af) then
            call fafWrite(ndiskr,residMO,8,ivirt*ivirt,1,ij,info)
          else
            call writer(ndiskr,ivirt,ij,residMO)
          endif
        else
          if (af) then
          call fafWrite(ndiskr,bl(iresidVOcurr),8,ivirt*ivirt,1,ij,info)
          else
            call writer(ndiskr,ivirt,ij,bl(iresidVOcurr))
          endif
        endif
        if (af) then
          call fafWrite(ndiskc,coefAO,8,idimen*idimen,1,ij,info)
        else
          call writer(ndiskc,idimen,ij,coefAO)
        endif
        call dynamic_retmem(2)
      endif
      if (WhatToDo.eq.'singlesstore') then
        ishift_record=irecord_tab(max_rec)
        if (vorb) then
          ishift=0
        else
          ishift=nmo
        endif
        do i=1,ivirt
          incf=i+ishift
          bl(isiglresstore+(jorb-1)*ivirt+i-1)=residMO(incf,1)
        enddo
        do i=1,idimen
          bl(isiglcoestore+(jorb-1)*idimen+i-1)=coefAO(i,1)
        enddo
        if (iorb.eq.nval.and.jorb.eq.nval) then
          do i=1,ivirt
            do j=nval+1,ivirt
              bl(isiglresstore+ivirt*(j-1)+i-1)=0.0d0
            enddo
          enddo
          do i=1,idimen
            do j=nval+1,idimen
              bl(isiglcoestore+(j-1)*idimen+i-1)=0.0d0
          enddo
        enddo
c ij=(npairs+1)*(ishift_record-1)+(npairs+1)=(npairs+1)*ishift_record
        ij=(npairs+1)*ishift_record
        if (af) then
          call fafWrite(ndiskr,bl(isiglresstore),8,ivirt*ivirt,1,
     *                  ij,info)
          call fafWrite(ndiskc,bl(isiglcoestore),8,idimen*idimen,1,
     *                  ij,info)
        else
          call writer(ndiskr,ivirt,ij,bl(isiglresstore))
          call writer(ndiskc,idimen,ij,bl(isiglcoestore))
        endif
        endif
      endif
c
      if (WhatToDo.eq.'calculate') then
        isize=(max_rec+2)*(max_rec+1)/2
        inew=1
        if (icalculate.eq.max_rec) then
c                rearange the matrix with Bij coef.
          do i=1,isize
            bl(idiis+inew-1)=bl(idiis+i-1)
            bl(idiis+i-1)=0.0d0
          if (.not.check_for_delete(i)) inew=inew+1
          enddo
        endif
        if (icalculate.lt.max_rec) icalculate=icalculate+1
        if (icalculate.gt.max_rec) STOP 'Too big DIIS matrix'
c Reserve memory for DIIS matrix, triangular, isize max=4*5/2=
        if (icalculate.eq.1.or.diis_restart) then
          if (icalculate.eq.1) call dynamic_getmem(max(isize,66),idiis)
          icalculate=1
          diis_restart=.false.
          do i=0,isize-1
          bl(idiis+i)=0.0d0
          enddo
          bl(idiis)=0.0d0
          idiis_last=idiis+1
          ifirst_rec=ishift_record-1
        endif
 444    continue
        bl(idiis_last)=-1.0d0
c Calculate DIIS matrix:
        if (af) then ! af .equiv. job is parallel
ccs        call sleep(5)
        call par_diis_matrix(singles,mp2_cc,npairs,
     *                       irecord_tab,max_rec,ndiskr,
     *                       idiis_last,icalculate,nbf)
        else
        call dynamic_getmem(ncf*ncf,iresidVOcurr)
        call dynamic_getmem(ncf*ncf,iresidVOprev)
        call diis_matrix(af,singles,mp2_cc,npairs,list,
     *                   ipairimages,irecord_tab,max_rec,ndiskr,ncf,
     *                   nmo,idiis_last,icalculate,bl(iresidVOcurr),
     *                   bl(iresidVOprev),nbf)
        call dynamic_retmem(2)
        endif
c
        if (icalculate.lt.max_rec) idiis_last=idiis_last+icalculate+1
        write (icond,*) 'icalculate',icalculate
        do i=1,(icalculate+1)*(icalculate+2)/2
          write (icond,*), 'Matrix: ' ,bl(idiis+i-1)
        enddo
        call flush(icond)
#if 1
555     continue
        if (icalculate.gt.1) then
          call matconn('diismat','s',icalculate+1,icalculate+1,idiis)
          call dynamic_matdef('diisinv','q',icalculate+1,icalculate+1)
          call matcopy('diismat','diisinv')
          idiisinv=mataddr('diisinv')
          call dynamic_getmem(1000,iad2)
          call dynamic_getmem(1000,iad3)
          call osinv(bl(idiisinv),icalculate+1,det,1.0d-20,
     *               bl(iad2),bl(iad3))
          call dynamic_retmem(2)
          write (icond,*) 'DET: ',det
          call flush(icond)
          if (dabs(det).lt.1.0d-200) then
            if (icalculate.eq.2) then
              write (icond,*) 'Determinant too small, '//
     *                        'DIIS fully restarted.'
              call flush(icond)
              call dynamic_matrem('diisinv')
              call matdisc('diismat')
              diis_restart=.true.
              if (af) then
              call diis_slave_next(2)
              endif
              return
            endif
            write (icond,*) 'Determinant too small, '//
     *                      'DIIS partly restarted.'
            call flush(icond)
            call dynamic_matrem('diisinv')
            call matdisc('diismat')
            isiz=icalculate+1
            bl(idiis)=0.0d0
            bl(idiis+1)=-1.0d0
            bl(idiis+2)=bl(idiis+((isiz-1)*(isiz-2)/2+isiz-2))
            bl(idiis+3)=-1.0d0
            bl(idiis+4)=bl(idiis+(isiz*(isiz-1)/2+isiz-2))
            bl(idiis+5)=bl(idiis+(isiz*(isiz-1)/2+isiz-1))
            do i=6,isize
            bl(idiis+i)=0.0d0
            enddo
            idiis_last=idiis+6
            ical_diffe=icalculate-2
            ifirst_rec=ifirst_rec+ical_diffe
            icalculate=2
            if (af) then
            call diis_slave_next(1)
            endif
            goto 555
          endif
c         call matinv('diisinv')
          do i=1,icalculate
          call matelem('diisinv',i+1,1,element)
          diiscoeff(i)=-element
          write (icond,*) 'Element ',diiscoeff(i)
          enddo
          call flush(icond)
          call dynamic_matrem('diisinv')
          call matdisc('diismat')
c
        endif !  if (icalculate.gt.1) then
c
         if (af) then
           call diis_slave_next(0)
ccs           call sleep(5)
           call par_new_coeff(singles,af,mp2_cc,idimen,nval,list,
     *                        ipairimages,icalculate,max_rec,
     *                        irecord_tab,npairs,ndiskc,diiscoeff)
ccs           call sleep(5)
         else
           call calc_new_coeff(singles,af,mp2_cc,idimen,nval,list,
     *                         ipairimages,icalculate,max_rec,
     *                         irecord_tab,npairs,ndiskc,diiscoeff)
         endif
#endif
      endif
      end
C==============================================================================
      function check_for_delete(ij)
c This function assignes a logical value to the number provided. This value 
c is used by DIIS subroutine for determining whether this element of triangular
c matrix should be deleted, see picture below, Y are element which should
c be kept, X are element for deletion
c O
c X X
c O X O
c O X O O
c O X O O O
c O X O O O O
c O X O O O O O
c O X O O O O O O    etc, we delete second row and column of triangular matrix
c if the ij number is from series: 2,3,5,8,12,17,23,30,38,47,57,68,80,93,107...
c the .true. value is returned
      implicit none
      logical check_for_delete
      integer ij,i,ii
      check_for_delete=.false. ! Default
      if (ij.eq.  2) then 
        check_for_delete=.true.
        return
      endif
      ii=0
      i=2
      do ! check if ij=2,3,5,8,12,17,23,30,38,47,57,68,80,93,107...
c                   we start from 3, 2 was checked before
        ii=ii+1
        i=i+ii
        if (ij.eq.i) then
          check_for_delete=.true.
          return
        else if (ij.lt.i) then
          check_for_delete=.false.
          return
        endif
      enddo
      print *,'Error in check_for_delete function'
      STOP    'Error in check_for_delete function'
      end
C===============================================================================
      subroutine diis_matrix(af,singles,mp2_cc,npairs,list,
     *                       ipairimages,irecord_tab,max_rec,ndiskr,ncf,
     *                       nmo,idiis_last,icalculate,residVOcurr,
     *                       residVOprev,nbf)
      use memory
      implicit none
      logical af,singles
      character*(*) mp2_cc
      integer npairs,list(*),ipairimages(*),irecord_tab(*),max_rec,nbf
      integer ndiskr,ncf,nmo,idiis_last,icalculate
      real*8 residVOcurr(nbf-nmo,nbf-nmo)
      real*8 residVOprev(nbf-nmo,nbf-nmo)
c
      integer ivirt,ij,irecord,info,i,k
      real*8 trace
c
        ivirt=nbf-nmo
c
        if (singles) then
        do ij=1,npairs+1
          if (ij.eq.(npairs+1).and. mp2_cc.eq.'mp2') exit
          if (ij.le.npairs) then
            if (list(ij).ne.1) cycle
            if (ipairimages(ij).eq.0) cycle
          endif
          irecord=(irecord_tab(max_rec)-1)*(npairs+1)+ij
          if (af) then
            call fafRead(ndiskr,residVOcurr,8,ivirt*ivirt,
     *                    1,irecord,info)
          if (info.lt.8*ivirt*ivirt) then
              ERROR
              call flush(6)
          endif
          else
            call reader(ndiskr,ivirt,irecord,residVOcurr)
          endif
          do i=1,icalculate ! OPTIMIZE! LAST RECORD = residVOcurr
            k=i+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*(npairs+1)+ij
            if (af) then
              call fafRead(ndiskr,residVOprev,8,ivirt*ivirt,
     *                      1,irecord,info)
          if (info.lt.8*ivirt*ivirt) then
              ERROR
              call flush(6)
          endif
            else
              call reader(ndiskr,ivirt,irecord,residVOprev)
            endif
              call tjmatprodtr1(residVOcurr,residVOprev,trace,ivirt)
            bl(idiis_last+i)=bl(idiis_last+i)+trace*ipairimages(ij)
          enddo
        enddo
        else !  if (singles) then
        do ij=1,npairs
          irecord=(irecord_tab(max_rec)-1)*npairs+ij
          if (ij.le.npairs) then
            if (list(ij).ne.1) cycle
            if (ipairimages(ij).eq.0) cycle
          endif
          if (af) then
            call fafRead(ndiskr,residVOcurr,8,ivirt*ivirt,
     *                    1,irecord,info)
          if (info.lt.8*ivirt*ivirt) then
              ERROR
              call flush(6)
          endif
          else
            call reader(ndiskr,ivirt,irecord,residVOcurr)
          endif
          do i=1,icalculate
            k=i+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*npairs+ij
            if (af) then
              call fafRead(ndiskr,residVOprev,8,ivirt*ivirt,
     *                      1,irecord,info)
          if (info.lt.8*ivirt*ivirt) then
              ERROR
              call flush(6)
          endif
            else
              call reader(ndiskr,ivirt,irecord,residVOprev)
            endif
            call tjmatprodtr1(residVOcurr,residVOprev,trace,ivirt)
            bl(idiis_last+i)=bl(idiis_last+i)+trace*ipairimages(ij)
          enddo
        enddo
        endif
      end
C===============================================================================
      subroutine calc_new_coeff(singles,af,mp2_cc,idimen,nval,list,
     *                          ipairimages,icalculate,max_rec,
     *                          irecord_tab,npairs,ndiskc,diiscoeff)
      use memory
      implicit none
      logical singles,af
      character mp2_cc*(*)
      integer idimen,nval,list(*),ipairimages(*),icalculate,max_rec
      integer irecord_tab(*),npairs,ndiskc
      real*8 diiscoeff(*)
c
      integer icoeff,ifinalcoeff,mataddr,ij,i,j,ii,k,irecord,info
c
c
      if (icalculate.le.1) return
      call dynamic_matdef('finalco','q',idimen,idimen)
      call dynamic_matdef('diiscoef','q',idimen,idimen)
      icoeff=mataddr('diiscoef')
      ifinalcoeff=mataddr('finalco')
      ij=0
      if (singles) then
      do i=1,nval
        do j=1,i
          ij=ij+1
          if (list(ij).ne.1) cycle
          if (ipairimages(ij).eq.0) cycle
          call matzero('finalco')
          do ii=1,icalculate
            k=ii+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*(npairs+1)+ij
            if (af) then
              call fafRead(ndiskc,bl(icoeff),8,idimen*idimen,
     *                        1,irecord,info)
          if (info.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
            else
              call reader(ndiskc,idimen,irecord,bl(icoeff))
            endif
            call matadd1('diiscoef',diiscoeff(ii),'finalco')
          enddo
          call CoefWrite(i,j,ifinalcoeff)
        enddo
      enddo
      if (mp2_cc.eq.'ccc') then
          ij=npairs+1
          call matzero('finalco')
          do ii=1,icalculate
c             irecord=(ifirst_rec+ii-1)*(npairs+1)+ij
            k=ii+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*(npairs+1)+ij
            if (af) then
              call fafRead(ndiskc,bl(icoeff),8,idimen*idimen,
     *                        1,irecord,info)
          if (info.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
            else
              call reader(ndiskc,idimen,irecord,bl(icoeff))
            endif
            call matadd1('diiscoef',diiscoeff(ii),'finalco')
          enddo
          do i=1,nval
c           call CoefWrite(i,j,ifinalcoeff)
          call writesingles(i,bl(ifinalcoeff+(i-1)*idimen))
          enddo
      endif
      else
      do i=1,nval
        do j=1,i
          ij=ij+1
          if (list(ij).ne.1) cycle
          if (ipairimages(ij).eq.0) cycle
          call matzero('finalco')
          do ii=1,icalculate
c             irecord=(ifirst_rec+ii-1)*npairs+ij
            k=ii+max_rec-icalculate
            irecord=(irecord_tab(k)-1)*npairs+ij
            if (af) then
              call fafRead(ndiskc,bl(icoeff),8,idimen*idimen,
     *                        1,irecord,info)
          if (info.lt.8*idimen*idimen) then
              ERROR
              call flush(6)
          endif
            else
              call reader(ndiskc,idimen,irecord,bl(icoeff))
            endif
            call matadd1('diiscoef',diiscoeff(ii),'finalco')
          enddo
          call CoefWrite(i,j,ifinalcoeff)
        enddo
      enddo
      endif
      call dynamic_matrem('diiscoef')
      call dynamic_matrem('finalco')
c
      end
C===============================================================================
      subroutine CoulombBuilder(i,j,idimen,npairs,ncf,nmo,vorb,xmatrix)
      use memory
      implicit none
c
      integer i,j,idimen,npairs,ncf,nmo
      logical vorb
      real*8  xmatrix(*)
c
      integer ivect,icoul,kk,ii,irecadrtc,ndisktc,lbintc,iii,kkk
      logical byt8
      real*8 thresh,xmult,x,xmax
c
      call dynamic_getmem(idimen,ivect)
      call dynamic_getmem(idimen*idimen,icoul)
      call readsingles(j,bl(ivect))
      thresh=1d-12
      byt8=.true.
      do kk=1,idimen*idimen
        xmatrix(kk)=0.0d0
      enddo
      do ii=1,idimen
        xmult=bl(ivect+ii-1)
        call ExtrOne(ii,    i,     irecadrtc,  npairs,    ndisktc,
     *               ncf,  lbintc,thresh,     byt8,      'tt',
     *               'mo', nmo,   vorb,       bl(icoul))
c check if it is coulomb:
        xmax=0d0
        do iii=1,idimen
          do kkk=1,idimen
          x=bl(icoul+(iii-1)*idimen+kkk-1)
          if (dabs(x).gt.xmax) xmax=dabs(x)
         x=bl(icoul+(iii-1)*idimen+kkk-1)-bl(icoul+(kkk-1)*idimen+iii-1)
            if (dabs(x).gt.1d-9) print *,'asym: ',ii,j,iii,kkk
          enddo
        enddo
        if (xmax.lt.1e-7) print *,'very small: ',ii,j
c end check
c
        do kk=1,idimen*idimen
          xmatrix(kk)=xmatrix(kk)+bl(icoul+kk-1)*xmult
        enddo
      enddo
      call dynamic_retmem(2)
      end

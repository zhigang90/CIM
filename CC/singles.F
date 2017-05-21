#include "maxi.h"
C===============================================================================
      subroutine initsingles(nval,ncf,singles,nmo,vorba,nbf)
c This subroutine initializes singles amplitudes IO subsystem. Singles 
c amplitudes are kept in memory, because they occupy little space.
c Initialization encompass memory allocation and commpn block setup. This
c common block is used by readsingles and writesingles subroutines.
c Parameters:
cI  nval     - the number of correlated orbitals
cI  ncf      - the number of contracted basis functions
cI  singles  - logical, if true sigles are calculated, otherwise this subroutine
c              is not needed.
cI  nmo      - the number of occupied HF orbitals
cI  vorba    - logical variable indicating how integrals and amplitudes 
c              are calculated in MO virtual space or in AO space.
c common block /singleskeep/, all parameters are output:
c iaddress0 - the bl memory address used for stoting amplitudes to be read
c ncfa      - the ncf is passed here
c singlesa  - the singles is passed here
c iaddress1 - the bl memory address used for stoting amplitudes to be written
c idimen    - the virtual space dimension, ncf or no of virt orbitals depending
c             on vorb value
c vorb      - vorba is passed here
      use memory
      use ccounters
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba
c
      logical singlesa
      integer iaddress0,iaddress1,ncfa,i,idimen,icheck0,icheck1
      integer iii
      logical vorb
      common /singleskeep/ iaddress0,ncfa,singlesa,iaddress1,idimen,
     *                     vorb
      save /singleskeep/
      logical signum
      integer ns,iprim,ishift0,ishiftprim1,iorb
c
      ic_initsingles=ic_initsingles+1
c
      vorb=vorba
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ncfa=ncf
      singlesa=singles
      if (singles.and.ic_initsingles.eq.1) then
        call dynamic_getmem(nval*idimen,iaddress0)
        call dynamic_getmem(nval*idimen,iaddress1)
        do i=0,nval*idimen-1
        bl(iaddress0+i)=0.0d0
c       bl(iaddress0+i)=1.d-10*dble(irand())
        enddo
      else if (singles.and.ic_initsingles.gt.1) then
        call dynamic_show_alloc(bl(iaddress0),icheck0)
        if (icheck0.ne.nval*idimen) 
     *    call nerror(1,'initsingles','Memory allocation error',
     *                ic_initsingles,0)
        call dynamic_show_alloc(bl(iaddress1),icheck1)
        if (icheck1.ne.nval*idimen) 
     *    call nerror(1,'initsingles','Memory allocation error',
     *                ic_initsingles,0)
#if 0
c Symmetrization
c       call symmetrize_the_same_coeff(ncf,bl(iaddress1))
        do iorb=1,nval
          ishift0=iaddress0+(iorb-1)*idimen
          call orb_searcher(iorb, iprim, ns, signum)
          ishiftprim1=iaddress1+(iprim-1)*idimen
          do i=0,idimen-1
            bl(ishift0+i)=bl(ishiftprim1+i)
          enddo
c         call symm_trans_vec(idimen, ns, signum, bl(ishift0))
          call transform_orbital(idimen,signum,ns,bl(ishift0))
        enddo
#else
c No symmetrization
        do i=0,nval*idimen-1
          bl(iaddress0+i)=bl(iaddress1+i)
        enddo
        call symm2_s_ampl(idimen, bl(iaddress0))
#endif
c     call ampl_symm(ncf, bl(iaddress0))
      endif
      end
C===============================================================================
C===============================================================================
      subroutine initsingles_local(nval,ncf,singles,nmo,vorba,nbf)
c This subroutine does almost exactly the same things as inisingles, but it is
c used by MP4. It was simplified for some reason, but I do not remember what
c reason. Maybe it should be simplre removed and replaced by inisingles.
      use memory
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba
c
      logical singlesa
      integer iaddress0,iaddress1,ncfa,i,idimen
      logical vorb
      common /singleskeep/ iaddress0,ncfa,singlesa,iaddress1,idimen,
     *                     vorb
c
c
      vorb=vorba
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ncfa=ncf
      singlesa=singles
c     if (singles.and.icounter.eq.1) then
c       call dynamic_getmem(nval*idimen,iaddress0)
c       call dynamic_getmem(nval*idimen,iaddress1)
c       do i=0,nval*idimen-1
c       bl(iaddress0+i)=0.0d0
c       bl(iaddress0+i)=1.d-10*dble(irand())
c       enddo
c     endif
c      if (singles.and.icounter.gt.0) then
        do i=0,nval*idimen-1
        bl(iaddress0+i)=bl(iaddress1+i)
        enddo
        call symm2_s_ampl(idimen, bl(iaddress0))
c      endif
      end
C===============================================================================
      subroutine readsingles(iorbital,vector)
c The singles reading. For common's description see initsingles
c Parameters:
cI  iorbital - the number of correlated orbital we need singles vector for
cO  vector   - the output singles vector, it is copied into this vector matrix
c
      use memory
      implicit none
      integer iorbital
      integer iaddress0,ncf,ibase,i,iaddress1,idimen
      logical singles,vorb
      common /singleskeep/ iaddress0,ncf,singles,iaddress1,idimen,vorb
      real*8 vector(idimen)
      ibase=iaddress0+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      vector(i)=bl(ibase+i)
      enddo
      else
      do i=1,idimen
      vector(i)=0.0d0
      enddo
      endif
c     call symm_trans_vec(idimen, ns, signum, vector)
      end
C===============================================================================
      subroutine readsinglesnew(iorbital,vector)
c The singles reading. It reads the singles amplitudes from CURRENT not 
c previous iteration. For common's description see initsingles
c Parameters:
cI  iorbital - the number of correlated orbital we need singles vector for
cO  vector   - the output singles vector, it is copied into this vector matrix
c
      use memory
      implicit none
      integer iorbital
      integer iaddress0,ncf,ibase,i,iaddress1,idimen
      logical singles,vorb
      common /singleskeep/ iaddress0,ncf,singles,iaddress1,idimen,vorb
      real*8 vector(idimen)
      ibase=iaddress1+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      vector(i)=bl(ibase+i)
      enddo
      else
      do i=1,idimen
      vector(i)=0.0d0
      enddo
      endif
      end
C===============================================================================
      subroutine writesingles(iorbital,vector)
c The singles writing. For common's description see initsingles
c Parameters:
cI  iorbital - the number of correlated orbital we need singles vector for
cI  vector   - the input singles vector, it is copied into internal memory 
c              location
c
      use memory
      implicit none
      integer iorbital
      integer iaddress0,ncf,ibase,i,iaddress1,idimen
      logical singles,vorb
      common /singleskeep/ iaddress0,ncf,singles,iaddress1,idimen,vorb
      real*8 vector(idimen)
      ibase=iaddress1+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      bl(ibase+i)=vector(i)
      enddo
      endif
      end
C===============================================================================
      subroutine pointersingles(ipointer)
c This subroutine returns the location of internal memory storage of singles 
c amplitudes.
      implicit none
      integer ipointer
      integer iaddress,ncf,iaddres1,idimen
      logical singles,vorb
      common /singleskeep/ iaddress,ncf,singles,iaddres1,idimen,vorb
      ipointer=iaddress
      end
C===============================================================================
      subroutine pointernewsingles(ipointer)
c This subroutine returns the location of internal memory storage of singles 
c amplitudes.
      implicit none
      integer ipointer
      integer iaddress,ncf,iaddres1,idimen
      logical singles,vorb
      common /singleskeep/ iaddress,ncf,singles,iaddres1,idimen,vorb
      ipointer=iaddres1
      end
C===============================================================================
c==============================================================================
      subroutine Kijk_Vec_Init(irecadrx, npairs, ndiskx, lbinx, ncf,
     *                         thresh,   ifirst, last,   byt8,  nmo,
     *                         vorb,     af,     ndiska, nbf)
c This subroutine calculates k(ijk) vectors, as
c defined by J. Chem. Phys. 104, 6286, (1996), and stores them on disk.
c It sets up common block /kijk/ which is used later by xkijk & xlijk
c reading subroutines.
c Parameters:
cI  irecadrx - the reocrd info for Kext integrals, note that both size of this
c              array and the contents differ in single and parallel mode
cI  npairs   - the number of correlated pairs
cI  ndiskx   - the file unit number where bins of Kext are kept
cI  lbinx    - the length of bin for Kext
cI  ncf      - the number of contracted basis functions
cI  thresh   - the integral threshold
cI  ifirst   - the first correlated orbital
cI  last     - the last correlated orbital usually HOMO
cI  byt8     - if .true. it means that double precision bins are used, otherwise
c              int4+int1 pair
cI  nmo      - the number of occupied orbitals
cI  vorb     - logical variable indicating how integrals and amplitudes 
c              are calculated: in MO virtual space or in AO space.
cI  af       - if .true. the job is parallel, Array Files is used as IO
c              device, otherwise job is single processor and local IO is used.
cO  ndiska   - the file unit number where k(ijk) are stored.
c   The common block /kijk/, all parameters are output:
c ndisk       - the same as ndiska above
c ncfa        - the same af ncf
c nval        - the number of corelated orbitals
c idimen      - the dimension od virtual space, might be ncf or ncf-nmo, 
c               depending on vorb value
c vorba       - the same as vorb
c afa         - af
c ndisk_cache - this file is used for local caching.
c
      use memory
      implicit none
      integer irecadrx,npairs,ndiskx,lbinx,ncf,ifirst,last,nmo,ndiska
      real*8 thresh
      logical byt8,vorb,vorba,af,afa
      dimension irecadrx(npairs)
c
      integer ndisk,nval,nval2,irecord,len,len1,mataddr,imatrix,iresult
      integer nrec,i,j,k,ncfa,idimen,iresultMO,info,ndisk_cache,iaddress
      integer nbf,nbfa
      common /kijk/ ndisk,ncfa,nval,idimen,vorba,afa,ndisk_cache,nbfa
      logical fileopen
      character*256 scrfile,filename
      character*2 ch_mygid
      integer mygid
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
      save /kijk/
c
      vorba=vorb
      afa=af
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ndisk_cache=76
c//   if (cache.and.af) then
      if (af) then
        call getival('mgid',mygid)
        write(ch_mygid,'(I2.2)') mygid
        call getchval('scrf',scrfile)
        call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
        filename=scrfile(1:len)//'.kijk--'//ch_mygid
        len1=len+9
        irecord=idimen*8
        OPEN(UNIT=ndisk_cache,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
c
      nval=last-ifirst+1
      nval2=nval*nval
      ncfa=ncf
      nbfa=nbf
      ndisk=76
      irecord=idimen*8
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                  ! 80 is starting lenght
      filename=scrfile(1:len)//'.kijk'
      len1=len+5
c
      if (af) then
         STOP 'This subroutine cannot be started in parallel mode'
         call fafCreatem(filename(1:len1)//char(0),ndisk,irecord)
         print *,'Kijk_Vec_Init fafcreate: ',ndisk
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 76 was used by somebody'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=irecord)
      endif
c
      ndiska=ndisk
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('matrix','q',ncf,ncf)
      imatrix=mataddr('matrix')
      call dynamic_matdef('result','r',ncf,nval)
      iresult=mataddr('result')
      call matsub('occu','cano',ifirst,last)
      if (vorb) then
        call dynamic_matdef('resultMO','r',idimen,nval)
        iresultMO=mataddr('resultMO')
      endif
      do i=1,nval
        do j=1,nval
          call ExtrOne(i,    j,     irecadrx,   npairs,    ndiskx,
     *                 ncf,  lbinx, thresh,     byt8,      'x',
     *                 'ao', nmo,   vorb,       bl(imatrix))
c         call matsub('occuK','occu',k,k)
          call matmmult('matrix','occu','result')
          if (vorb) then
            call matmmul2('genvirt','result','resultMO','t','n','n')
          endif
c         call matdisc('occuK')
          do k=1,nval
            nrec=(i-1)*nval2+(j-1)*nval+k
            if (vorb) then
c             call matmmul2('genvirt','result','resultMO','t','n','n')
              iaddress=iresultMO+idimen*(k-1)
            else
              iaddress=iresult+idimen*(k-1)
            endif
            if (af) then
            call fafWrite(ndisk,bl(iaddress),8,idimen,1,nrec,info)
            call vectorwrite(ndisk_cache,nrec,
     *                                    bl(iaddress),idimen)
            else
            call vectorwrite(ndisk,nrec,bl(iaddress),idimen)
            endif
          enddo
        enddo
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine vectorwrite(ndisk,nrec,vector,idimen)
      implicit none
      integer ndisk,nrec,idimen
      real*8 vector(idimen)
      write(ndisk,REC=nrec) vector
      end
c==============================================================================
      subroutine vectorread(ndisk,nrec,vector,idimen)
      implicit none
      integer ndisk,nrec,idimen
      real*8 vector(idimen)
      read(ndisk,REC=nrec) vector
      end
c==============================================================================
      subroutine xkijk(i,j,k,vector)
c This subroutine simply reads K(ijk) integrals from local file and deposits
c result in vector. See descr. of commons in Kijk_Vec_Init
      implicit none
      integer i,j,k
      integer ndisk,ncf,nval,idimen,ndisk_cache,nbf
      logical vorb,af
      common /kijk/ ndisk,ncf,nval,idimen,vorb,af,ndisk_cache,nbf
      real*8 vector(idimen)
      integer nrec,info
      logical cache,integerT
      common /GlobalCCSD/ cache,integerT
c
c
      nrec=(i-1)*nval*nval+(j-1)*nval+k
      if (af) then
c//     if (cache) then
          call vectorread(ndisk_cache,nrec,vector,idimen)
c//     else
c//       call fafRead(ndisk,vector,8,idimen,1,nrec,info)
c//     endif
      else
        call vectorread(ndisk,nrec,vector,idimen)
      endif
      end
c==============================================================================
      subroutine xlijk(i,j,k,vector)
c This subroutine simply reads K(ijk) integrals from local file,
c builds l(ijk) and deposits result in vector. See descr. of commons
c in Kijk_Vec_Init
      implicit none
      integer i,j,k
      integer ndisk,ncf,nval,idimen,ndisk_cache,nbf
      logical vorb,af
      common /kijk/ ndisk,ncf,nval,idimen,vorb,af,ndisk_cache,nbf
      real*8 vector(idimen), work1(idimen),work2(idimen)
      integer ii
      call xkijk(i,j,k,work1)
      call xkijk(j,i,k,work2)
      do ii=1,idimen
        vector(ii)=2.0d0*work1(ii)-work2(ii)
      enddo
      end
c==============================================================================
      subroutine Fock_Vector_Init(ncf,ifirst,last,ifockAO,ifockMO,nmo,
     *                            vorba,nbf)
c This subroutine prepares CCSD program for some special usage of Fock matrix. 
c In few places
c of CCSD formulas there are needed vectors of Fock matrix. i-th vector
c is defined as composed of F(mu,i) elements, where i is fixed MO index and 
c mu is the vector AO index.
c In MO space it is F(a,i), where 'a' stands for virtual orbital. 
c If Fock is the solution of Fock equations (both localized or not), it is
c equal to zero.
c Parameters:
c ncf     - the number of contracted functions
c ifirst  - the first contracted orbital
c last    - the last contracted orbital
c ifockAO - the bl index of Fock matrix. See comments below.
c nmo     - the number of occupied HF orbitals
cI  vorba - logical variable indicating how integrals and amplitudes 
c           are calculated: in MO virtual space or in AO space.
c Commons: /fockHalf/, all parameters are output:
c ifockHalf - the memory location where the bunch of Fock vectors is
c             stored
c ncfa      - = ncf
c nval      - the number of correlated orbitals
c idimen    - the dimension of virtual space.
c vorb      - = vorba
c
      use memory
      implicit none
      integer ncf,ifirst,last,ifockAO,nmo,ifockMO,nbf
      logical vorb,vorba
      integer ifockHalf,ncfa,nval,idimen,iii,ii,jj
      common /fockHalf/ ifockHalf,ncfa,nval,idimen,vorb
      save /fockHalf/
c
      vorb=vorba
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      nval=last-ifirst+1
      ncfa=ncf
      call dynamic_getmem(idimen*nval,ifockHalf)
c
      call matconn('fockAO1','q',ncf,ncf,ifockAO) ! invalid for vorb
      call dynamic_matdef('halfhalf','r',ncf,nval)
      call matconn('fockHalf','r',idimen,nval,ifockHalf)
      call matsub('occu','cano',ifirst,last)
      if (vorb) then
        do ii=1,nval
          do jj=1,idimen
            bl(ifockHalf+(jj-1)+(ii-1)*idimen)=
     *      bl(ifockMO+  (last+jj-1)+(ifirst+ii-2)*ncf)
          enddo
        enddo
      else
        call matmmult('fockAO1','occu','halfhalf')  !check the col order
        call matcopy('halfhalf','fockHalf')
      endif
      call matdisc('occu')
      call matdisc('fockHalf')
      call dynamic_matrem('halfhalf')
      call matdisc('fockAO1')
      call dynamic_lock(bl(ifockHalf),iii)
      end
c==============================================================================
      subroutine Fock_Vector_Read(i,vector)
      use memory
      implicit none
      integer i
      real*8 vector(*)
      logical vorb
      integer ifockHalf,ncf,nval,iaddress,ii,idimen,isize
      common /fockHalf/ ifockHalf,ncf,nval,idimen,vorb
      call dynamic_show_alloc(bl(ifockHalf),isize)
      if (isize.ne.idimen*nval) call nerror(2,'Fock_Vector_Read',
     *                              'Memory allocation error',0,0)
      iaddress=ifockHalf+(i-1)*idimen-1
      do ii=1,idimen
      vector(ii)=bl(iaddress+ii)
      enddo
      end
c==============================================================================
      subroutine Fock_Vector_Pointer(iaddr)
      use memory
      implicit none
      integer iaddr
      logical vorb
      integer ifockHalf,ncf,nval,iaddress,ii,idimen,isize
      common /fockHalf/ ifockHalf,ncf,nval,idimen,vorb
      iaddr=ifockHalf
      end
c==============================================================================
      subroutine singles_CISD(i,        j,      ncf,    nval,  ioverlap,
     *                        irecadre, npairs, nfirst, nlast, ndiske,
     *                        lbine,    thresh, byt8,   ccsd,  iresult,
     *                        nmo,      vorb,   nbf)
c This subroutine calculates the singles contributions to G term as defined in
c J. Chem. Phys. 104, 6286, (1996), where by singles contributions I mean the
c terms which contribute directly to G term, but not those which contribute via
c auxillary matrices, as e.g. Y Z or X.
c Parameters:
cI  i        - the number of correlatet i orbital, first index of G matrix
cI  j        - ditto, second index
cI  ncf      - the number of contracted functions
cI  nval     - the number of correlated orbitals
cI  ioverlap - the bl index to overlap matrix
cI  irecadre - the records info for EEO operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - the number of correlated pairs
cI  nfirst   - the first correlated pair
cI  nlast    - the last correlated pair
cI  ndiske   - the file unit number, where EEO are stored in form of bins
cI  lbine    - the EEO's bin length
cI  thresh   - the integral threshold
cI  byt8     - if .true. bins are double precision, otherwise int4+int1
cI  ccsd     - if .true. the singles contribution of CCSD are calculated, this
c              variable imples also in general that CCD is calculated
cO  iresult  - the result is returned here, into this bl location
cI  nmo      - the number of occupied HF orbitals
cI  vorb     - logical variable indicating how integrals and amplitudes 
c              are calculated: in MO virtual space or in AO space.
c
      use memory
      implicit none
      integer i,j,ncf,nval,ioverlap,irecadre,npairs,nfirst,nlast,nmo
      integer ndiske,lbine,iresult,nbf
      real*8 thresh
      logical byt8,ccsd,vorb
      integer isinglevec,ifock,mataddr,k,iEEOread,idimen,iEEOrao,ii
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
      call matconn('S','q',ncf,ncf,ioverlap)
      call matconn('finalij','q',idimen,idimen,iresult)
      call dynamic_matdef('work22','q',idimen,idimen)
      call dynamic_matdef('singlvec','r',idimen,1)
      call dynamic_matdef('fockvec','r',idimen,1)
      isinglevec=mataddr('singlvec')
      ifock     =mataddr('fockvec')
      if (ccsd) then
        call dynamic_matdef('EEOrao','q',ncf,ncf)
        iEEOrao=mataddr('EEOrao')
        call dynamic_matdef('EEOread','r',idimen,ncf)
        iEEOread=mataddr('EEOread')
        call matsub('occu','cano',nfirst,nlast)
      endif
c
      call readsingles(i,bl(isinglevec))
      if (ccsd) then
        call readS(j,bl(ifock))
      else
        call Fock_Vector_Read(j,bl(ifock))
      endif
      call matmmul2('singlvec','fockvec','work22','n','t','n')
      if (ccsd) then
        call ExtrOne(j,    i,     irecadre,   npairs,    ndiske,
     *               ncf,  lbine, thresh,     byt8,      'e',
     *               'ao', nmo,   vorb,       bl(iEEOrao))
        if (vorb) then
          call matmmul2('genvirt','EEOrao','EEOread','t','n','n')
        else
          call matcopy('EEOrao','EEOread')
        endif
      endif
      do k=1,nval
      call readsingles(k,bl(isinglevec))
      call xkijk(j,i,k,bl(ifock))
      if (ccsd) then
        call matsub('onecol','occu',k,k)
        call matmmul2('EEOread','onecol','fockvec','n','n','a')
        call matdisc('onecol')
      endif
      call matscal('singlvec',-1.0d0)
      call matmmul2('singlvec','fockvec','work22','n','t','a')
      enddo
      if (vorb) then
        call matcopy('work22','finalij')
      else
        call matmmul2('S','work22','finalij','n','n','n')
      endif
c
      call readsingles(j,bl(isinglevec))
      if (ccsd) then
        call readS(i,bl(ifock))
      else
      call Fock_Vector_Read(i,bl(ifock))
      endif
      call matmmul2('fockvec','singlvec','work22','n','t','n')
      if (ccsd) then
        call ExtrOne(i,    j,     irecadre,   npairs,    ndiske,
     *               ncf,  lbine, thresh,     byt8,      'e',
     *               'ao', nmo,   vorb,       bl(iEEOrao))
        if (vorb) then
          call matmmul2('genvirt','EEOrao','EEOread','t','n','n')
        else
          call matcopy('EEOrao','EEOread')
        endif
      endif
      do k=1,nval
      call readsingles(k,bl(isinglevec))
      call xkijk(i,j,k,bl(ifock))
      if (ccsd) then
        call matsub('onecol','occu',k,k)
        call matmmul2('EEOread','onecol','fockvec','n','n','a')
        call matdisc('onecol')
      endif
      call matscal('singlvec',-1.0d0)
      call matmmul2('fockvec','singlvec','work22','n','t','a')
      enddo
      if (vorb) then
      call matadd('work22','finalij')
      else
      call matmmul2('work22','S','finalij','n','n','a')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine Tl_generator(ncf,nval,nmo,vorb,ioverl,imatrix,nbf)
c This suroutine calculates the T x l term, according to 
c J. Chem. Phys. 104, 6286, (1996) it is part of formula no 44. 
c The output is written using writeS subroutine to storage which stores
c later the whole (44) term.
c Parameters:
cI  ncf    - the number of contracted functions
cI  nval   - the number of correlated orbitals
cI  nmo    - the number of occupied molecular orbitals
cI  vorb   - .true.: the dimension of virtual space is equal to no of virt.
c            orbitals, and MO form of integral matrices is used, otherwise,
c            the virtual space has ncf dimension and all calculations are 
c            performed for matrices in AO form.
      use memory
      implicit none
      integer ncf,nval,nmo,imatrix,ioverl,nbf
      logical vorb
      integer k,l,i,iresult,idoubles,iwork11,iwork22,mataddr,idimen
      integer isingles,iwork
c
      call mmark
      call dynamic_mmark
      call matmark
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
c     call dynamic_getmem(idimen*nval,isingles)
      call matconn('fin_res','r',idimen,nval,imatrix)
      call matconn('xoverlap','q',ncf,ncf,ioverl)
      call dynamic_matdef('result','r',idimen,nval)
      isingles=mataddr('result')
      call zeroit(bl(isingles),idimen*nval)
c
      call dynamic_matdef('work11','r',idimen,nval)
      iwork11=mataddr('work11')
      do k=1,nval
        do l=1,k
          call CoefRead('tt',l,k,idoubles)
          call matconn('doubl_am','q',idimen,idimen,idoubles)
          do i=1,nval
            iwork=iwork11+(i-1)*idimen
            call xlijk(k,l,i,bl(iwork))
          enddo
          call matmmul2('doubl_am','work11','result','n','n','a')
          if (k.ne.l) then
            do i=1,nval
              iwork=iwork11+(i-1)*idimen
              call xlijk(l,k,i,bl(iwork))
            enddo
            call matmmul2('doubl_am','work11','result','t','n','a')
          endif
          call matdisc('doubl_am')
        enddo
      enddo
      if (vorb) then
        call matadd1('result',-1d0,'fin_res')
      else
        call matscal('result',-1.0d0)
        call matmmul2('xoverlap','result','fin_res','n','n','a')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine 
     *   EEO_vector_extractor(nval,  irecadre, npairs,  ndiske, ncf,
     *                        lbine, thresh,   nfirst,  nlast,  byt8,
     *                        nmo,   vorb,     iresult, nbf)
c The subroutine calculates those contributions to 's' vector which involves
c EEO matrices - it construcrs EEO vectors, according to 
c J. Chem. Phys. 104, 6286, (1996)
cI  nval     - the number of correlated orbitals
cI  irecadre - the records info for EEO operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - the number of correlated pairs
cI  ndiske   - the file unit number, where EEO bins are kept.
cI  ncf      - the number of contracted basis functions
cI  lbine    - the length of EEO bins.
cI  thresh   - the integral threshold
cI  nfirst   - the fist correlated orbital
cI  nlast    - the last correlated orbital
cI  byt8     - .true.: integrals bins are REAL*8, otherwise ibin4+ibin1
cI  nmo      - the number of occupied HF orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
cO  iresult  - the result is deposited here, into this bl array index location
c
      use memory
      implicit none
      integer nval,npairs,ndiske,nbf
      integer ncf,lbine,nfirst,nlast,iresult,nmo
      integer irecadre(npairs)
      real*8 thresh
      logical byt8,vorb
c
      integer i,k,iEEOread,iEEOconv,mataddr,idimen,intermid,ires
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
      call matsub('occu','cano',nfirst,nlast)
c
      call dynamic_matdef('EEOread','q',ncf,ncf)
      iEEOread=mataddr('EEOread')
c
      call dynamic_matdef('EEOconv','q',ncf,ncf)
      iEEOconv=mataddr('EEOconv')
c
      call dynamic_matdef('intermid','r',ncf,nval)
      call matzero('intermid')
      intermid=mataddr('intermid')
c
      call matconn('t_result','r',idimen,nval,iresult)
c
      do k=1,nval
        do i=1,k
          call ExtrOne(i,    k,     irecadre,   npairs,    ndiske,
     *                 ncf,  lbine, thresh,     byt8,      'e',
     *                 'ao', nmo,   vorb,       bl(iEEOread))
          call CCLoperator(bl(iEEOread),ncf,bl(iEEOconv))
          call matsub('onecol','occu',k,k)
          ires=intermid+(i-1)*ncf
          call matconn('result','r',ncf,1,ires)
          call matmmul2('EEOconv','onecol','result','n','n','a')
          call matdisc('result')
          call matdisc('onecol')
c
        if (i.ne.k) then
          call matsub('onecol','occu',i,i)
          ires=intermid+(k-1)*ncf
          call matconn('result','r',ncf,1,ires)
          call matmmul2('EEOconv','onecol','result','t','n','a')
          call matdisc('result')
          call matdisc('onecol')
        endif
        enddo
      enddo
      if (vorb) then
        call matmmul2('genvirt','intermid','t_result','t','n','a')
      else
        call matadd('intermid','t_result')
      endif
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine Tr_generator(i,       nval,  ncf,    irecadrx, npairs,
     *                        ndiskx,  lbinx, thresh, byt8,     cc,
     *                        iresult, nmo,   vorb)
c This subroutine calculates the second term of (38) formula taken from
c J. Chem. Phys. 104, 6286, (1996)
c Parameters:
cI  i        - the requested vector number
cI  nval     - number of correlated orbitals
cI  ncf      - number of contracted basis functions
cI  irecadrx - the records info for Kext operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - number of correlated pairs
cI  ndiskx   - the file unit number where Kext bins are kept.
cI  lbinx    - the length of Kext bin
cI  thresh   - the integral threshold
cI  byt8     - .true.: double precision bins, otherwise int4+int1
cI  cc       - .true.: the CCD contributions are calculated, otherwise CID or 
c              CISD
cO  iresult  - the result is deposited here
cI  nmo      - the number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      implicit none
      integer i,nval,ncf,npairs,ndiskx,lbinx,iresult,nmo
      integer irecadrx(npairs)
      real*8 thresh
      logical byt8,cc,vorb
c
      integer k,l,idoubles,ifock_k,mataddr,icoef_dif,iwork11,idimen
      integer icoeffic,iL,iexchint
c
      call mmark
      call dynamic_mmark
      call matmark
c
      STOP 'Switched off'
      if (vorb) then
        idimen=ncf-nmo
      else
        idimen=ncf
      endif
c
      call dynamic_matdef('coef_dif','q',idimen,idimen)
      icoef_dif=mataddr('coef_dif')
      call matconn('result','r',idimen,1,iresult)
      call matzero('result')
      call dynamic_matdef('fock_k','r',idimen,1)
      ifock_k=mataddr('fock_k')
      do k=1,nval
        call matzero('coef_dif')
        call CoefRead('tt',i,k,idoubles)
        call matconn('coeffi','q',idimen,idimen,idoubles)
        call matadd1('coeffi',2.0d0,'coef_dif')
        call matpose('coeffi')
        call matadd1('coeffi',-1.0d0,'coef_dif')
        call matpose('coeffi')
        call matdisc('coeffi')
        call readR(k,bl(ifock_k))
        call matmmul2('coef_dif','fock_k','result','n','n','a')
      enddo
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine Tr_generator1(idimen,nval,ncf,vorb,ioverlap,iresult)
c This subroutine calculates the second term of (38) formula taken from
c J. Chem. Phys. 104, 6286, (1996)
c Parameters:
cI  i        - the requested vector number
cI  nval     - number of correlated orbitals
cI  ncf      - number of contracted basis functions
cI  irecadrx - the records info for Kext operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - number of correlated pairs
cI  ndiskx   - the file unit number where Kext bins are kept.
cI  lbinx    - the length of Kext bin
cI  thresh   - the integral threshold
cI  byt8     - .true.: double precision bins, otherwise int4+int1
cI  cc       - .true.: the CCD contributions are calculated, otherwise CID or 
c              CISD
cO  iresult  - the result is deposited here
cI  nmo      - the number of molecular orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
      use memory
      implicit none
      integer idimen,nval,ncf,iresult,ioverlap
      logical vorb
c
      integer iiR,icoef_dif,ires_tmp,i,k,idoubles,ir_vec,ires_vec
      integer mataddr
c
      call mmark
      call dynamic_mmark
      call matmark
c
      call matconn('fin_res','r',idimen,nval,iresult)
      call matconn('overl','q',ncf,ncf,ioverlap)
      call dynamic_matdef('coef_dif','q',idimen,idimen)
      icoef_dif=mataddr('coef_dif')
      call pointRR(iiR)
      call dynamic_matdef('res_tmp','r',idimen,nval)
      call matzero('res_tmp')
      ires_tmp=mataddr('res_tmp')
#if 1
      do i=1,nval
        do k=1,i
          call CoefRead('tt',i,k,idoubles)
          call CCLoperator(bl(idoubles),idimen,bl(icoef_dif))
          ir_vec=iiR+(k-1)*idimen
          ires_vec=ires_tmp+(i-1)*idimen
          call matconn('RR','r',idimen,1,ir_vec)
          call matconn('result','r',idimen,1,ires_vec)
          call matmmul2('coef_dif','RR','result','n','n','a')
          call matdisc('result')
          call matdisc('RR')
          if (i.ne.k) then
            ir_vec=iiR+(i-1)*idimen
            ires_vec=ires_tmp+(k-1)*idimen
            call matconn('RR','r',idimen,1,ir_vec)
            call matconn('result','r',idimen,1,ires_vec)
            call matmmul2('coef_dif','RR','result','t','n','a')
            call matdisc('result')
            call matdisc('RR')
          endif
        enddo
      enddo
#endif
      if (vorb) then
        call matadd('res_tmp','fin_res')
      else
        call matmmul2('overl','res_tmp','fin_res','n','n','a')
      endif
c
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine beta_t_sumator(i,     ncf,   irecadrx,npairs, ndiskx,
     *                          lbinx, thresh,byt8,    nval,   FMO,
     *                          nfirst,cc,    ccsd,    iresult,nmo,
     *                          vorb,  betam)
c This subroutine calculates the third term of (38) formula taken from
c J. Chem. Phys. 104, 6286, (1996)
c Parameters:
c     Most of them are legacy things and are not used rigth now.
c  I keep them because it doesn't cost anything and for possible future use.
c
cI  i        - the vector number
cI  ncf      - number of contracted basis functions
cI  irecadrx - the records info for Kext operator, note: different size and
c              different information is contained in this array for single mode 
c              and parallel mode.
cI  npairs   - number of correlated pairs
cI  ndiskx   - the file unit number where Kext bins are kept.
cI  lbinx    - the length of Kext bin
cI  thresh   - the integral threshold
cI  byt8     - .true.: double precision bins, otherwise int4+int1
cI  nval     - number of correlated orbitals
cI  FMO      - The Fock matrix in MO basis
cI  nfirst   - the first correlated orbital
cI  cc       - .true.: the CCD contributions are calculated, otherwise CID or 
cI  ccsd     - .true.: the CC has to be true and CCSD singles contributions are
c               added
cO  iresult  - the reult
cI  nmo      - the number of occupied HF orbitals
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
cI  betam    - beta matrix (defined by (41))
      use memory
      implicit none
      integer i,ncf,npairs,ndiskx,lbinx,nval,nfirst,iresult,nmo
      integer irecadrx(npairs)
      real*8 FMO(ncf,ncf),thresh,betam(nval,nval)
      logical byt8,cc,ccsd,vorb
c
      integer k,ip,kp,mataddr,icoeffic,idimen
      real*8 beta
c
      stop 'Switched off'
      if (vorb) then
        idimen=ncf-nmo
      else
        idimen=ncf
      endif
c     call print_pointer(betam)
c
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('coeffic','r',idimen,1)
      icoeffic=mataddr('coeffic')
      call matconn('result','r',idimen,1,iresult)
      call matzero('result')
      do k=1,nval
        call readsingles(k,bl(icoeffic))
        call matadd1('coeffic',-betam(k,i),'result')
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine UpdateSingles(i,       icoefAO, iresidAO,nmo,   nfirst,
     *                         ncf,     npairs,  nval,    epsi,  shift,
     *                         squares, vorb,    xmax,    ndiskr,ndiskc,
     *                         list,ipairimages,nbf)
c This subroutine performs singles update using residua and HF orbital
c energies.
c Parameters:
cI  i        - the updated vector index
cIO icoefAO  - the bl location of amplitudes, in spite of name it can be 
c              AO or MO, depending on vorb value
cI  iresidAO - the bl location of residua
cI  nmo      - the number of occupied HF orbitals
cI  nfirst   - the first correlated orbital
cI  ncf      - number of contracted basis functions
cI  npairs   - number of correlated pairs
cI  nval     - number of correlated orbitals
cI  epsi     - the linear vector containing all HF eigenvalues
cI  shift    - the level shift used for better separation of occ and virt. 
c              orbitals
cO  squares  - the sum of singles residua squares
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
cO  xmax     - the largest residuum element
cI  ndiskr   - the disk unit numbers passed to DIIS subroutine.
cI  ndiskc   - ditto
c
      use memory
      implicit none
      integer i,icoefAO,iresidAO,nmo,nfirst,ncf,npairs,nval,nbf
      integer ndiskr,ndiskc
      real*8 epsi(ncf),shift,squares,xmax
      logical vorb,singles
c
      integer ip,ia,mataddr,iresMO,idimen,iadd,list(*),ipairimages(*)
      real*8 res_elem,coe_elem,residuum
c
c
      singles=.true.
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      call mmark
      call dynamic_mmark
      call matmark
      call matconn('veccoef','r',idimen,1,icoefAO)
      call matconn('vecresi','r',idimen,1,iresidAO)
      call dynamic_matdef('resMO','r',idimen,1)
      call dynamic_matdef('coeMO','r',idimen,1)
      iresMO=mataddr('resMO')
      if (vorb) then
        call matcopy('vecresi','resMO')
        call matcopy('veccoef','coeMO')
      else
        call matmmult('canotran','vecresi','resMO')
        call matmmult('rinvcano','veccoef','coeMO')
      endif
      ip=i+nfirst-1
c
      if (vorb) then
        iadd=nmo
      else
        iadd=0
      endif
c
      do ia=1,idimen
        if (ia.le.nmo.and.(.not.vorb)) then
          call mateset('coeMO',ia,1,0.0d0)
          cycle
        endif
        call matelem('resMO',ia,1,res_elem)
        call matelem('coeMO',ia,1,coe_elem)
        residuum=-res_elem/(epsi(ia+iadd)-epsi(ip)+shift)
        coe_elem=coe_elem+residuum
        call mateset('coeMO',ia,1,coe_elem)
ctmp    call mateset('resMO',ia,1,residuum)
        squares=squares+res_elem*res_elem
        if (dabs(res_elem).gt.xmax) xmax=res_elem
      enddo
      if (vorb) then
      call matcopy('coeMO','veccoef')
      else
      call matmmult('cano','coeMO','veccoef')
      endif
      call CCDiis(i,       i,          'singlesstore', ncf,  nmo,
     *            npairs,  bl(iresMO), bl(icoefAO),    nval, vorb,
     *            singles, ndiskr,     ndiskc,         'xxx',list,
     *            ipairimages,nbf)
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine IntTEIO   (ncf,  nval,      mu,        lam,  iprnt,
     1                      nrow, ncol,      irow,      icol, thresh, 
     2                      xmat, halfmulam, halflammu, nmo,  vorb,nbf)
c This subroutine transforms AO exchange or coulomb integrals into
c Three-External-Integral-Operators (TEIO), by multiplication by 
c correlated HF orbitals and singles amplitudes. This subroutine is not
c optimized for speed, as it cost has rank size**5 not size**6 as the most
c expensive parts of CCSD. Most of the time goes into integral calculation
c so far.
c Here only one exchange or coulomb matrix is transformed, and two output
c matrices are generated. Those matrices are then sorted, but sorting
c subroutine is subsequent to this subroutine call.
c Parameters:
c   ncf       - number of contracted basis functions
c   nval      - number of correlated orbitals
c   mu        - the AO first index of xmat matrix
c   lam       - the AO second index of xmat matrix
c   iprnt     - the print level
c   nrow      - the number of rows of compressed xmat
c   ncol      - the number of columns of compressed xmat
c   irow      - the contraction info for rows
c   icol      - ditto for columns
c   thresh    - the integral threshold
c   xmat      - the matrix to be transformed (mu, lam matrix)
c   halfmulam - the transformed matrix (all ij, mu lam index)
c   halflammu - ditto lam mu index, unlike MP2 integrals there is no symmetry
c               that K(ij)(t)=K(ji), so we need both matrices.
c   nmo       - the number of occupied molecular orbitals
cI  vorb      - .true.: the dimension of virtual space is equal to no of virt.
c               orbitals, and MO form of integral matrices is used, otherwise,
c               the virtual space has ncf dimension and all calculations are 
c               performed for matrices in AO form.
c
      use memory
      implicit none
      integer ncf,nval,mu,lam,iprnt,nrow,ncol,irow,icol,nmo,nbf
      real*8 thresh,xmat(ncf,ncf),halfmulam(nval,nval) 
      real*8 halflammu(nval,nval)
      logical vorb
c
      integer iexchAO,mataddr,iaddress,iaddressMO,i,j
      call mmark
      call dynamic_mmark
      call matmark
      call dynamic_matdef('exchAO','q',ncf,ncf)
      iexchAO=mataddr('exchAO')
      call dynamic_matdef('work11','r',nval,ncf)
      call decompact(xmat,bl(iexchAO),nrow,ncol,irow,icol,ncf) ! temporary ?
      if (vorb) then
        call pointersingles(iaddressMO)
        call matconn('singlMO','r',nbf-nmo,nval,iaddressMO)
        call dynamic_getmem(ncf*nval,iaddress)
        call matconn('singles','r',ncf,nval,iaddress)
        call matmmult('genvirt','singlMO','singles')
        call matdisc('singles')
        call matdisc('singlMO')
      else
        call pointersingles(iaddress)
      endif
      call matconn('singles','r',ncf,nval,iaddress)
      call matmmul2('occu','exchAO','work11','t','n','n')
      call matmmult('work11','singles','halfmulam')
c     call matprint('halfmulam',99)
      call matpose('exchAO')
      call matmmul2('occu','exchAO','work11','t','n','n')
      call matmmult('work11','singles','halflammu')
c     call matprint('halflammu',99)
c     STOP
 133  FORMAT (4I3,F20.10)
      call matremark
      call retmark
      call dynamic_retmark
      end
c==============================================================================
      subroutine BinSortExt(nval, lbin,npairs, ndisk, halfmulam,
     *                   halflammu,
     1                   mu,      lam,   thresh, icounter,irecord,
     2                   indxbin, ibin4, ibin1, bin8,  afflag,  afname,
     3                   islvid, byt8,irec_no,irec_last)
c See BinSort desciption, this subroutine is slightly modified in order to
c be able to handle IntTEIO subroutine output. 
        implicit real*8 (a-h,o-z)
        real*8 halfmulam(nval,nval),halflammu(nval,nval)
        integer*2 indxbin(2,lbin,nval*nval)
        integer*4 ibin4(lbin,nval*nval)
        integer*1 ibin1(lbin,nval*nval)
        real*8    bin8(lbin,nval*nval)
        integer irecord(nval*nval),icounter(nval*nval),irec_no,irec_last
        logical afflag,byt8
        parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
        parameter (dblinv=one/dblmax,d1max=128.0d0)
c
        ij=0
        do i=1,nval
          do j=1,nval
            ij=ij+1
            abhlftra=abs(halfmulam(i,j))
            if(abhlftra.gt.thresh) then
              ii=icounter(ij)+1
              icounter(ij)=ii
              indxbin(1,ii,ij)=mu
              indxbin(2,ii,ij)=lam
              x=halfmulam(i,j)
              byte8: if (byt8) then
              bin8(ii,ij)=halfmulam(i,j)
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
              if (ii.eq.lbin) then
                if(afflag) then 
                  call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indxbin(1,1,ij),irecord(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
                else
                  call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij), indxbin(1,1,ij),irecord(ij),ij,nval*nval,
     2              byt8)
                endif
                icounter(ij)=0
              end if
            end if
          end do
        end do
        if (mu.ne.lam) then
          ij=0
          do i=1,nval
            do j=1,nval
              ij=ij+1
              abhlftra=abs(halflammu(i,j))
              if(abhlftra.gt.thresh) then
                ii=icounter(ij)+1
                icounter(ij)=ii
                indxbin(1,ii,ij)=lam
                indxbin(2,ii,ij)=mu
                x=halflammu(i,j)
                byte81: if (byt8) then
                bin8(ii,ij)=halflammu(i,j)
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
                if (ii.eq.lbin) then
                  if(afflag) then 
                  call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indxbin(1,1,ij),irecord(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
                  else
                    call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indxbin(1,1,ij),irecord(ij), ij,nval*nval,
     2                byt8)
                  endif
                  icounter(ij)=0
                endif
              end if
            end do
          end do
        end if
      end
C===============================================================================
      subroutine lockS
      use memory
      implicit none
      integer iaddress,ncf,idimen,nval,iii
      logical singles,vorb
      common /skeep/ iaddress,ncf,singles,idimen,vorb,nval
      call dynamic_lock(bl(iaddress),iii)
      end
C===============================================================================
      subroutine unlockS
      use memory
      implicit none
      integer iaddress,ncf,idimen,nval,iii
      logical singles,vorb
      common /skeep/ iaddress,ncf,singles,idimen,vorb,nval
      call dynamic_unlock(bl(iaddress),iii)
      end
C===============================================================================
      subroutine lockR
      use memory
      implicit none
      integer iaddress,ncf,idimen,nval,iii
      logical singles,vorb
      common /rkeep/ iaddress,ncf,singles,idimen,vorb,nval
      call dynamic_lock(bl(iaddress),iii)
      end
C===============================================================================
      subroutine unlockR
      use memory
      implicit none
      integer iaddress,ncf,idimen,nval,iii
      logical singles,vorb
      common /rkeep/ iaddress,ncf,singles,idimen,vorb,nval
      call dynamic_unlock(bl(iaddress),iii)
      end
C===============================================================================
      subroutine initS(nval,ncf,singles,nmo,vorba,nbf)
c See initsingles subroutine, this one does the same thing but for s
c verctors J. Chem. Phys. 104, 6286, (1996), formula no. 44.
      use memory
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba
c
      logical singlesa,vorb
      integer iaddress,ncfa,i,idimen,nvala
      common /skeep/ iaddress,ncfa,singlesa,idimen,vorb,nvala
      save /skeep/
c
c
      vorb=vorba
      nvala=nval
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
      ncfa=ncf
      singlesa=singles
      if (singles) then
        call dynamic_getmem(nval*idimen,iaddress)
        do i=0,nval*idimen-1
        bl(iaddress+i)=0.0d0
        enddo
      endif
      end
C===============================================================================
      subroutine pointS(iadd)
c See initS, IO subroutine for s vectors
      use memory
      implicit none
      integer iadd
      integer iaddress,ncf,idimen,nval
      logical singles,vorb
      common /skeep/ iaddress,ncf,singles,idimen,vorb,nval
      iadd=iaddress
      end
C===============================================================================
      subroutine readS(iorbital,vector)
c See initS, IO subroutine for s vectors
      use memory
      implicit none
      integer iorbital
      integer iaddress,ncf,ibase,i,idimen,iamount,nval
      logical singles,vorb
      common /skeep/ iaddress,ncf,singles,idimen,vorb,nval
      real*8 vector(idimen)
      call dynamic_show_alloc(bl(iaddress),iamount)
      if (iamount.ne.idimen*nval) call nerror(2,'readS',
     *                              'Memory allocation error',0,0)
      ibase=iaddress+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      vector(i)=bl(ibase+i)
      enddo
      else
      do i=1,idimen
      vector(i)=0.0d0
      enddo
      endif
      end
C===============================================================================
      subroutine writeS(iorbital,vector)
c See initS, IO subroutine for s vectors
      use memory
      implicit none
      integer iorbital
      integer iaddress,ncf,ibase,i,idimen,iamount,nval
      logical singles,vorb
      common /skeep/ iaddress,ncf,singles,idimen,vorb,nval
      real*8 vector(idimen)
      call dynamic_show_alloc(bl(iaddress),iamount)
      if (iamount.ne.idimen*nval) call nerror(2,'writeS',
     *                              'Memory allocation error',0,0)
      ibase=iaddress+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      bl(ibase+i)=vector(i)
      enddo
      endif
      end
C===============================================================================
C===============================================================================
      subroutine initR(nval,ncf,singles,nmo,vorba,nbf)
c See initsingles subroutine, this one does the same thing but for r
c verctors J. Chem. Phys. 104, 6286, (1996), formula no. 43.
      use memory
      implicit none
      integer nval,ncf,nmo,nbf
      logical singles,vorba
c
      logical singlesa,vorb
      integer iaddress,ncfa,i,idimen,nvala
      common /rkeep/ iaddress,ncfa,singlesa,idimen,vorb,nvala
      save /rkeep/
c
      vorb=vorba
      nvala=nval
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
c
c
      ncfa=ncf
      singlesa=singles
      if (singles) then
        call dynamic_getmem(nval*idimen,iaddress)
        do i=0,nval*idimen-1
        bl(iaddress+i)=0.0d0
        enddo
      endif
      end
C===============================================================================
      subroutine readR(iorbital,vector)
c See initR, IO subroutine for r vectors
      use memory
      implicit none
      integer iorbital
      integer iaddress,ncf,ibase,i,idimen,iamount,nval
      logical singles,vorb
      common /rkeep/ iaddress,ncf,singles,idimen,vorb,nval
      real*8 vector(idimen)
      call dynamic_show_alloc(bl(iaddress),iamount)
      if (iamount.ne.idimen*nval) call nerror(2,'readR',
     *                              'Memory allocation error',0,0)
      ibase=iaddress+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      vector(i)=bl(ibase+i)
      enddo
      else
      do i=1,idimen
      vector(i)=0.0d0
      enddo
      endif
      end
C===============================================================================
      subroutine writeRR(iorbital,vector)
c See initR, IO subroutine for r vectors
      use memory
      implicit none
      integer iorbital
      integer iaddress,ncf,ibase,i,idimen,iamount,nval
      logical singles,vorb
      common /rkeep/ iaddress,ncf,singles,idimen,vorb,nval
      real*8 vector(idimen)
      call dynamic_show_alloc(bl(iaddress),iamount)
      if (iamount.ne.idimen*nval) call nerror(2,'writeRR',
     *                              'Memory allocation error',0,0)
      ibase=iaddress+(iorbital-1)*idimen
      ibase=ibase-1      ! because of the method of counting in the loop below
      if (singles) then
      do i=1,idimen
      bl(ibase+i)=vector(i)
      enddo
      endif
      end
C===============================================================================
      subroutine pointRR(iaddr)
c See initR, IO subroutine for r vectors
      use memory
      implicit none
      integer iaddr
      integer iaddress,ncf,ibase,i,idimen,iamount,nval
      logical singles,vorb
      common /rkeep/ iaddress,ncf,singles,idimen,vorb,nval
      iaddr=iaddress
      end
C==============================================================================
      subroutine Lt_generator(idimen,nval,   irecadrx,npairs,ndiskx,
     *                        ncf,   lbinx,  thresh,  byt8,  nmo,
     *                        vorb,  iresult)
      use memory
      implicit none
c
      integer idimen,nval,irecadrx(*),npairs,ndiskx,ncf,lbinx
      real*8 thresh
      logical byt8,vorb
      integer nmo,iresult
c
      integer isingles,iL,k,l,iwork11,ising,ires,mataddr
c
      call dynamic_mmark
      call matmark
c
      call dynamic_matdef('L','q',idimen,idimen)
      iL=mataddr('L')
      call dynamic_matdef('work11','q',idimen,idimen)
      iwork11=mataddr('work11')
      call pointersingles(isingles)
#if 1
      do k=1,nval
        do l=1,k
          call ExtrOne(k,    l,     irecadrx,     npairs,     ndiskx,
     *                 ncf,  lbinx, thresh,       byt8,       'x',
     *                 'mo', nmo,   vorb,         bl(iwork11))
          call CCLoperator(bl(iwork11),idimen,bl(iL))
c
          ising=isingles+(l-1)*idimen
          ires =iresult +(k-1)*idimen
          call matconn('singl1','r',idimen,1,ising)
          call matconn('res1','r',idimen,1,ires)
          call matmmul2('L','singl1','res1','n','n','a')
          call matdisc('res1')
          call matdisc('singl1')
          if (k.ne.l) then
            ising=isingles+(k-1)*idimen
            ires =iresult +(l-1)*idimen
            call matconn('singl1','r',idimen,1,ising)
            call matconn('res1','r',idimen,1,ires)
            call matmmul2('L','singl1','res1','t','n','a')
            call matdisc('res1')
            call matdisc('singl1')
          endif
        enddo
      enddo
#endif
      call matremark
      call dynamic_retmark
      end
C==============================================================================
      subroutine beta_t_sumator1(ncf,idimen,nval,ibeta,vorb,ioverlap,
     *                          iresult)
      use memory
      implicit none
c
      integer ncf,idimen,nval,ibeta,ioverlap,iresult
      logical vorb
c
      integer iloc_res,isingles,mataddr
c
      call dynamic_mmark
      call matmark
c
      call matconn('loc_beta','q',nval,nval,ibeta)
      call matconn('overl','q',ncf,ncf,ioverlap)
      call matconn('ext_res','r',idimen,nval,iresult)
      call dynamic_matdef('loc_res','r',idimen,nval)
      iloc_res=mataddr('loc_res')
      call pointersingles(isingles)
      call matconn('loc_sing','r',idimen,nval,isingles)
      call matmmul2('loc_sing','loc_beta','loc_res','n','n','n')
      if (vorb) then
        call matadd1('loc_res',-1d0,'ext_res')
      else
        call matscal('loc_res',-1d0)
        call matmmul2('overl','loc_res','ext_res','n','n','a')
      endif
c
      call matremark
      call dynamic_retmark
      end
C==============================================================================
      subroutine cisd_energy_add(ncf,idimen,nval,xiterenergy,vorb,
     *                           ioverlap,iresult)
      use memory
      implicit none
c
      integer ncf,idimen,nval,ioverlap,iresult
      logical vorb
      real*8 xiterenergy
c
      integer isingles
c
c
      call dynamic_mmark
      call matmark
c
      call pointersingles(isingles)
      call matconn('loc_sin','r',idimen,nval,isingles)
      call matconn('loc_ove','q',ncf,ncf,ioverlap)
      call matconn('loc_res','r',idimen,nval,iresult)
      call dynamic_matdef('loc_work','r',idimen,nval)
      if (vorb) then
        call matadd1('loc_sin',-xiterenergy,'loc_res')
      else
        call matmmul2('loc_ove','loc_sin','loc_work','n','n','n')
        call matadd1('loc_work',-xiterenergy,'loc_res')
      endif
c
      call matremark
      call dynamic_retmark
      end
c==============================================================================
      subroutine UpdateSingles1(icoefAO, iresidAO,nmo,   nfirst,
     *                         ncf,     npairs,  nval,    epsi,  shift,
     *                         squares, vorb,    xmax,    ndiskr,ndiskc,
     *                         list,ipairimages,icoefAO_final,
     *                         iorev,ivrev,docansym,t1diagnostic,nbf)
c This subroutine performs singles update using residua and HF orbital
c energies.
c Parameters:
cI  i        - the updated vector index
cIO icoefAO  - the bl location of amplitudes, in spite of name it can be 
c              AO or MO, depending on vorb value
cI  iresidAO - the bl location of residua
cI  nmo      - the number of occupied HF orbitals
cI  nfirst   - the first correlated orbital
cI  ncf      - number of contracted basis functions
cI  npairs   - number of correlated pairs
cI  nval     - number of correlated orbitals
cI  epsi     - the linear vector containing all HF eigenvalues
cI  shift    - the level shift used for better separation of occ and virt. 
c              orbitals
cO  squares  - the sum of singles residua squares
cI  vorb     - .true.: the dimension of virtual space is equal to no of virt.
c              orbitals, and MO form of integral matrices is used, otherwise,
c              the virtual space has ncf dimension and all calculations are 
c              performed for matrices in AO form.
cO  xmax     - the largest residuum element
cI  ndiskr   - the disk unit numbers passed to DIIS subroutine.
cI  ndiskc   - ditto
c
      use memory
      implicit none
      integer icoefAO,iresidAO,nmo,nfirst,ncf,npairs,nval,nbf
      integer ndiskr,ndiskc,icoefAO_final
      real*8 epsi(ncf),shift,squares,xmax,t1diagnostic
      integer iorev(*),ivrev(*)
      logical vorb,singles,docansym
c
      integer ip,ia,mataddr,iresMO,idimen,iadd,list(*),ipairimages(*)
      integer i,ire,ico,irep,iap,iarep
      real*8 res_elem,coe_elem,residuum
c
c
      singles=.true.
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      t1diagnostic=0d0
c
      call mmark
      call dynamic_mmark
      call matmark
      call matconn('veccoef','r',idimen,nval,icoefAO)
      call matconn('final_c','r',idimen,nval,icoefAO_final)
      call matconn('vecresi','r',idimen,nval,iresidAO)
      call dynamic_matdef('resMO','r',idimen,nval)
      call dynamic_matdef('coeMO','r',idimen,nval)
      iresMO=mataddr('resMO')
      if (vorb) then
        call matcopy('vecresi','resMO')
        call matcopy('veccoef','coeMO')
      else
        call matmmult('canotran','vecresi','resMO')
        call matmmult('rinvcano','veccoef','coeMO')
      endif
c
      if (vorb) then
        iadd=nmo
      else
        iadd=0
      endif
c
      do i=1,nval
      ip=i+nfirst-1
      if (docansym) irep=iorev(i)
      do ia=1,idimen
        iap=ia+iadd
        if (ia.le.nmo.and.(.not.vorb)) then
          call mateset('coeMO',ia,i,0.0d0)
          cycle
        endif
        if (docansym) then
          iarep=ivrev(iap-nmo)
          if (irep.ne.iarep) then 
            call mateset('coeMO',ia,i,0.0d0)
            cycle
          endif
        endif
        call matelem('resMO',ia,i,res_elem)
        call matelem('coeMO',ia,i,coe_elem)
        residuum=-res_elem/(epsi(iap)-epsi(ip)+shift)
        coe_elem=coe_elem+residuum
        call mateset('coeMO',ia,i,coe_elem)
ctmp    call mateset('resMO',ia,1,residuum)
        squares=squares+res_elem*res_elem
        t1diagnostic=t1diagnostic+coe_elem*coe_elem
        if (dabs(res_elem).gt.xmax) xmax=res_elem
      enddo
      enddo
      if (vorb) then
      call matcopy('coeMO','final_c')
      else
      call matmmult('cano','coeMO','final_c')
      endif
      do i=1,nval
      ire=iresMO       +(i-1)*idimen
      ico=icoefAO_final+(i-1)*idimen
      call CCDiis(i,       i,          'singlesstore', ncf,  nmo,
     *            npairs,  bl(ire), bl(ico),    nval, vorb,
     *            singles, ndiskr,     ndiskc,         'xxx',list,
     *            ipairimages,nbf)
      enddo
      call matremark
      call retmark
      call dynamic_retmark
      t1diagnostic=dsqrt(t1diagnostic/2d0/dble(nmo-nfirst+1))
      end
c==============================================================================

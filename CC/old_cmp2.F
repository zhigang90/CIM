c This is an old MP2 stuff, slightly modified for CCSD purposes.
c The MP2 procedures were modifed by Alan Ford for Array Files.
c I extracted them in order to be more independednt of possible MP2
c modifications. This way CCSD might be merged with PQS more easily.
c It might be removed in the future, when CCSD will become "builtin" part 
c of PQS. Then any MP2 modifications will have to take into account the fact,
c that CCSD uses those subroutines.
c But they will stay here if MP2 won't use Array Files in the future and old MP2
c implementation will be maintained.
c==============================================================================
c 2007-01-08    print *,'Ext3 Binsort' removed
      subroutine TransOneShell_tj(ncf,   ncs,     nval,   ics,   kcs,
     *                        icf1,  icf2,     kcf1,   kcf2,  ndisk,
     1                        inx,   xint,     CMO,    iprnt, thresh,
     2                        xmat,  halftra,  halftra1, mulam,  mulamd,
     *                        nrow,
     3                        ncol,  irow,     icol,   irow1,  icol1,
     4                    ENonZero,  smal,     lbin,   npairs, icounter,
     5                     irecord,  indxbin,  ibin4,  ibin1,  bin8, 
     6                      afflag,  afname,  islvid,  int_kind,byt8,
     7                      nmo, vorb,igranules,igranulesize,igran2pair,
     8                      ipair2gran,irec_no,irec_last,nbf)
      use memory
      implicit real*8 (a-h,o-z)
c  transform a batch of integrals
c  Arguments:
c  INTEN=IN:
c  ncf  number of contracted basis functions
c  ncs  number of contracted shells
c  nval number of orbitals to be transformed (usually the valence ones)
c  ics,kcs: contracted shells for which the transformation is carried out
c  inx: array containing contraction info
c  xint: AO integrals ordered as
c        (l=1:ncf,j=1:ncf,k=kcf1:kcf2,i=icf1:icf2)   where icf1 & icf2
c        are the first and the last contr. funct. in the ICS shell ,
c        similarly for kcs .
c        NOTE THAT THEY ARE SCALED by 1/thresh
c  CMO:  MO coefficients
c  iprnt: print level
c  thresh = integral threshold
c  nrow   = number of non-zero rows of the AO integral matrices
c  ncol   = ditto for columns
c  irow(k),k=1,nrow = the indices of the non-zero rows
c  icol(k),k=1,nco  = ditto for columns
c  STORAGE:
c  xmat: (ncf,ncf) place for 1/4 transformed integrals
c  intent=OUT:
c  halftr: space for half-transformed integrals to be written on disk
c  nrec:  the number of records written  (input/output)
c  mulam: number of non-diagonal (mu, lambda) pairs transformed
c  integrals are {mu,nu|lam,isig);  mu.ne.lam
c  mulamd: number of diagonal pairs (mu,mu) transformed
c  ENonZero = number of the elements in non-zero columns and rows
c
      dimension xint(*),xmat(ncf,ncf),halftra(*),halftra1(*),irow(ncf)
      dimension icol(ncf)
      dimension CMO(ncf,*),irow1(ncf),icol1(ncf),icounter(*)
      dimension inx(12,*)
c     common /big/bl(30000)
      integer*2 indxbin(*)
      integer*4 ibin4(*)
      real*8    bin8(*)
      integer*1 ibin1(*)
      logical smal,byt8,vorb,afflag
      character*(*) int_kind,afname
      integer irec_no,irec_last
c  return if there are no integrals
      if(nrow.eq.0.or.ncol.eq.0) RETURN
c
c  intstore holds the half-transformed integrals in integer form
c  int1     holds a 1-byte integer
c          (mimic 5-byte integers to allow lower integral threshold)
c
      call getint(nval**2,intstore)
      call getmem(nval**2/8+1,int1)
c
c  icol and jcol store integer arrays which  hold the indices
c  of the non-zero elements of the AO exchange matrix X=(mu,nu|lam,isig)
c  itrunc serves as temporary storage for a compacted coefficient matrix
c  only the rows (or columns in the second quarter transformation) of
c  X which are non-zero are  present in bl(itrunc)
c
      ivirt=mataddr('virt')
      ioccu=mataddr('occu')
      if (int_kind.eq.'tt') then
        call getmem(ncf*ncf,itrunc)
      else
        call getmem(ncf*nval,itrunc)
      endif
c
      do icf=icf1,icf2
      do kcf=kcf1,kcf2
      if(kcf.gt.icf) exit
      if(smal) then
        call ExtractOne_tj(ncf,   icf,   kcf,  icf1,  icf2,
     2                  kcf1,  kcf2,  xint,  xmat,  nrow,
     3                  ncol,  irow,  icol,  nrow1, ncol1,
     4                  irow1, icol1)
      endif
c  Determine some statistics:
          ENonZero=ENonZero+dble(nrow*ncol)
c
          if(nrow1.eq.0.or.ncol1.eq.0) cycle
          if (int_kind.eq.'c'.or.int_kind.eq.'x') then
          call IntTransp_tj(ncf, nval,  icf,    kcf,      iprnt,
     1                 ndisk, nrow1, ncol1,  irow1,    icol1,
     2                 thresh,xmat,  CMO, bl(intstore),bl(int1),
     3                 itrunc,halftra,nrec)
!  Sort the half-transformed integrals
          call BinSort_tj( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1, bin8, afflag, afname,
     3                  islvid,byt8,irec_no,irec_last)
          else if (int_kind.eq.'tt') then
          call ext3_IntTransp(ncf, nval,  icf,    kcf,      iprnt,
     1                 ndisk, nrow1, ncol1,  irow1,    icol1,
     2                 thresh,xmat,  bl(ioccu), bl(ivirt), bl(intstore),
     3                 bl(int1), itrunc,halftra,nrec,nmo,nbf)
!  Sort the half-transformed integrals
          call ext3_BinSort( nval,     lbin, npairs, ndisk, halftra,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1, bin8, afflag, afname,
     3                  islvid,byt8,ncf,nmo,igranules,igranulesize,
     4                  igran2pair,ipair2gran,irec_no,irec_last,nbf)
          else
      call IntTEIO   (ncf,    nval,    mu,      lam,     iprnt,
     1                      nrow1,    ncol1,    irow1,    icol1,
     2                      thresh,  xmat, halftra,  halftra1,nmo,vorb,
     3                      nbf)
!  Sort the half-transformed integrals
          call BinSortExt(nval,lbin, npairs, ndisk, halftra,halftra1,
     1                  icf,      kcf,  thresh,icounter, irecord,
     2                  indxbin,  ibin4,ibin1, bin8, afflag, afname,
     3                  islvid,byt8,irec_no,irec_last)
          endif
          if(icf.ne.kcf) then
            mulam=mulam+1
          else
            mulamd=mulamd+1
          end if
        end do
      end do
      call retmem(3)
      end
c===========================================================================
      subroutine ExtractOne_tj( ncf,   icf,   kcf,   icf1,  icf2,
     1                       kcf1,  kcf2,  xint,  xmat,  nrow,
     2                       ncol,  irow,  icol,  nrow1, ncol1,
     3                       irow1, icol1)
c  This routine extracts a single matrix from the integral array
c  The latter is indexed as xint(l=1:ncf,j=1:ncf,k=kcf1:kcf2,ncf,i=icf1:icf2);
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
c
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
          if(abs(xmat(jj,ll)).gt.zero) then
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
          if(abs(xmat(jj,ll)).gt.zero) then
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
c===========================================================================
      subroutine ext3_BinSort(nval,    lbin,  npairs, ndisk,   halftra,
     1                   mu,      lam,   thresh, icounter,irecord,
     2                   indxbin, ibin4, ibin1, bin8,  afflag,  afname,
     3                   islvid, byt8,ncf,nmo,igranules,igranulesize,
     4                   igran2pair,ipair2gran,irec_no,irec_last,nbf)
        implicit real*8 (a-h,o-z)
        real*8 halftra(nbf-nmo,nval)
        integer*2 indxbin(2,lbin,nval*(nbf-nmo))
        integer*4 ibin4(lbin,nval*(nbf-nmo))
        integer*1 ibin1(lbin,nval*(nbf-nmo))
        real*8    bin8(lbin,nval*(nbf-nmo))
        integer irecord(nval*(nbf-nmo)),icounter(nval*(nbf-nmo))
        logical afflag,byt8
        parameter (one=1.0d0,dblmax=2 147 483 648.0d0)
        parameter (dblinv=one/dblmax,d1max=128.0d0)
        character*(*) afname
        integer*4 igran2pair(2,*),ipair2gran(*)
        integer irec_no,irec_last
!  This routine sorts the half-transformed integrals
!  Arguments
!  INPUT
!  nval:   number of occupied MOs to be transformed (usually the
!             valence ones)
!  lbin:  number of records per bin
!  irecord: stores the number of the last written record for an ij
!  icounter: the counter for the number of items written to and ij bin 
!  npairs: number of occupied orbital pairs, nval*(nval+1)/2
!  ndisk: unit number or the disk containing sorted integrals
!  halftra:  the half transformed integrals
!  mu:    an atomic orbital descriptor
!  lam:   an atomic orbital descriptor
!  thresh:  integral threshold
!  indxbin:  stores the sorted indices to the AO's (i.e. mu, lam)
!  ibin4:  stores the sorted integrals
!  ibin1:  stores overflow
!  afflag:  determines if array files should be used.
!  afname:  name of array file
!  sid:  a unique index assigned to each slave used with array files
        ij=0
        itotb=nval*(nbf-nmo)
        do i=1,nbf-nmo
          do j=1,nval
            ij=ij+1
            if (icounter(ij).eq.lbin) then
              if(afflag) then 
                igranule=ipair2gran(ij)
                istart=igran2pair(1,igranule)
                istop =igran2pair(2,igranule)
                isize=istop-istart+1
                call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),
     1                          ibin1(1,istart),bin8(1,istart),
     2                          indxbin(1,1,istart),irecord(igranule),
     3                          igranule,islvid,byt8,isize,irec_no,
     4                          irec_last)
                do kk=istart,istop
                  icounter(kk)=0
                enddo
              else
                call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij), indxbin(1,1,ij),irecord(ij),ij,itotb,
     2            byt8)
                icounter(ij)=0
              endif
            end if
            abhlftra=abs(halftra(i,j))
            if(abhlftra.gt.thresh) then
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
          do i=1,nbf-nmo
            do j=1,nval
              ij=ij+1
              if (icounter(ij).eq.lbin) then
                if(afflag) then 
                  igranule=ipair2gran(ij)
                  istart=igran2pair(1,igranule)
                  istop =igran2pair(2,igranule)
                  isize=istop-istart+1
                  call AfWriteBinGr(ndisk,lbin,ibin4(1,istart),
     1                            ibin1(1,istart),bin8(1,istart),
     2                            indxbin(1,1,istart),irecord(igranule),
     3                            igranule,islvid,byt8,isize,irec_no,
     4                            irec_last)
                  do kk=istart,istop
                    icounter(kk)=0
                  enddo
                else
                  call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1              bin8(1,ij), indxbin(1,1,ij),irecord(ij),ij,itotb,
     2              byt8)
                  icounter(ij)=0
                endif
              end if
              abhlftra=abs(halftra(i,j)) ! because they are coulomb type
              if(abhlftra.gt.thresh) then
                ii=icounter(ij)+1
                icounter(ij)=ii
                indxbin(1,ii,ij)=lam
                indxbin(2,ii,ij)=mu
                x=halftra(i,j)
                byte81: if (byt8) then
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
                endif byte81
              end if
            end do
          end do
        end if
      end
!====================================================================
      subroutine IntTransp_tj(ncf,    nval,    mu,      lam,     iprnt,
     1                     ndisk,  nrow,    ncol,    irow,    icol,
     2                     thresh, xmat,    coef,    intint,  int1,
     3                     itrunc, halftra,nrec)
c   This routine transforms the integrals for a given mu, lam to MO basis,
c   yielding (mu,j|lam,l) where mu , lam are AO indices and j,l are MO
c   ones. It uses halftra both as a matrix 'halftra' and as an array
c   halftra for storing a temporary transformed matrix
c   The resulting integrals are written to disk as integers
c   ARGUMENTS:
c   INPUT
c  ncf: number of contracted basis functions
c  nval: number of orbitals to be transformed (valence orbitals)
c  mu,lam: fixed basis function indices. The transformed integrals are
c          (mu,i|lam,j) where i,j are MO indices
c  iprnt: print level
c  ndisk: the unit number of disk to write on
c  xmat = xmat(nu.isig)=(mu.nu | lam,isig) originally. However,
c         it is compacted here: xmat(p,q)=(mu,irow(p)|lam,icol(q))
c  note that xmat is also accessed by the matrix 'xmat'
c  coef(ncf,nval): SCF coefficients
c  nrow, ncol: the number of nonzero rows and columns of xmat
c  irow, icol = the indices of nonzero rows + columns
c  thresh = integral threshold
c  STORAGE
c  intint: storage for nval**2 4-byte integers
c  int1:  storage for nval**2 1-byte integers
c  itrunc = address of an (ncf x  nval) array. This holds a matrix
c  of the truncated SCF coefficients
c  OUTPUT
c  halftra: half-transformed integrals
c  nrec: number of records written so far (input/output)
c
      use memory
      implicit real*8 (a-h,o-z)
      integer*4 intint(nval,nval)
      integer*1 int1(nval,nval),i1
      dimension halftra(nval,nval),irow(nrow),icol(ncol)
      dimension xmat(ncf,ncf), coef(ncf,nval)
c     common /big/bl(30000)
      parameter (one=1.0d0)
c  defined in the calling program
c  xmat has the COMPACTED AO exchange matrix
CPP  print nrow,ncol,irow,icol
c     print *, 'mu,lam=',mu,lam,' nrow,ncol=',nrow,ncol
c     print *, 'irow',irow
c     print *, 'icol',icol
c    prepare the compacted SCF coefficients for the left side
      call CompactCoef_tj(ncf,nval,coef,nrow,irow,bl(itrunc))
c  perform first matrix multiplication. Note that the matrix
c  system routine matsimtr cannot be used here because the right and
c  left matrices may be different after compaction
c  define temporary matrix product
cc      iout=igetival('iout')
      iout = 6     ! do NOT change this   JB
      call matdef('temp','r',nval,ncol)
c   'ymat' is similar to 'xmat' but is compacted
      call matdef('ymat','r',nrow,ncol)
      call PutTrunc_tj(ncf,nrow,ncol,xmat,bl(mataddr('ymat')))
c  connect the truncated SCF coefficient matrix to the matrix system
      call matconn('trunc','r',nrow,nval,itrunc)
      call matmmul2('trunc','ymat','temp','t','n','n')
      call matdisc('trunc')
      call matrem('ymat')
      call CompactCoef_tj(ncf,nval,coef,ncol,icol,bl(itrunc))
      call matconn('trunc','r',ncol,nval,itrunc)
      call matmmult('temp','trunc','halftra')
      call matdisc('trunc')
      call matrem('temp')
c
      if(iprnt.gt.4) then
        call matscal('halftra',thresh)
        print *, mu,lam,nrow,ncol
        call matprint('halftra',iout)
        call matscal('halftra',one/thresh)
      end if
      return
      end
c==========================================================================
      subroutine BinSort_tj(nval,    lbin,  npairs, ndisk,   halftra,
     1                   mu,      lam,   thresh, icounter,irecord,
     2                   indxbin, ibin4, ibin1, bin8,  afflag,  afname,
     3                   islvid, byt8,irec_no,irec_last)
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
        integer irec_no,irec_last
!  This routine sorts the half-transformed integrals
!  Arguments
!  INPUT
!  nval:   number of occupied MOs to be transformed (usually the
!             valence ones)
!  lbin:  number of records per bin
!  irecord: stores the number of the last written record for an ij
!  icounter: the counter for the number of items written to and ij bin 
!  npairs: number of occupied orbital pairs, nval*(nval+1)/2
!  ndisk: unit number or the disk containing sorted integrals
!  halftra:  the half transformed integrals
!  mu:    an atomic orbital descriptor
!  lam:   an atomic orbital descriptor
!  thresh:  integral threshold
!  indxbin:  stores the sorted indices to the AO's (i.e. mu, lam)
!  ibin4:  stores the sorted integrals
!  ibin1:  stores overflow
!  afflag:  determines if array files should be used.
!  afname:  name of array file
!  sid:  a unique index assigned to each slave used with array files
        ij=0
        do i=1,nval
          do j=1,i
            ij=ij+1
            abhlftra=abs(halftra(i,j))
            if(abhlftra.gt.thresh) then
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
              if (ii.eq.lbin) then
                if(afflag) then 
                  call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indxbin(1,1,ij),irecord(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
                else
                  call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1              bin8(1,ij), indxbin(1,1,ij),irecord(ij),ij,npairs,
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
            do j=1,i
              ij=ij+1
              abhlftra=abs(halftra(j,i))
              if(abhlftra.gt.thresh) then
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
                if (ii.eq.lbin) then
                  if(afflag) then 
                  call AfWriteBin(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1            bin8(1,ij),indxbin(1,1,ij),irecord(ij),ij,islvid,byt8,
     2            irec_no,irec_last)
                  else
                    call WriteBin_tj(ndisk,lbin,ibin4(1,ij),ibin1(1,ij),
     1                bin8(1,ij),indxbin(1,1,ij),irecord(ij), ij,npairs,
     2                byt8)
                  endif
                  icounter(ij)=0
                endif
              end if
            end do
          end do
        end if
      end
!====================================================================
      SUBROUTINE WRITEBIN_TJ(NDISK,  LBIN, IBIN4, IBIN1, BIN8, INDXBIN,
     1                    IRECIJ, IJ,   NPAIRS, BYT8)
!     This routine writes a full bin on disk
!     Arguments:
!     INPUT:
!     ndisk: unit number
!     lbin: number of records per bin
!     ibin4:  compressed integrals
!     ibin1: integral overflow
!     indxbin:  sorted AO indices
!     OUTPUT:
!     iirecij:  the record written for the ij pair
        IMPLICIT INTEGER(I-N)
        INTEGER*4 IBIN4(LBIN)
        INTEGER*1 IBIN1(LBIN)
        REAL*8    BIN8(LBIN)
        INTEGER*2 INDXBIN(2,LBIN)
        INTEGER*4 IPOS
        LOGICAL BYT8
        PARAMETER(NSEQ=10) !number of sequential records for a pair
       
        IF(IRECIJ.EQ.0) THEN
          IPOS=(IJ-1)*NSEQ+1
        ELSE
          LEVL=(IRECIJ-1)/(NPAIRS*NSEQ)
          LVLSTART=LEVL*NPAIRS*NSEQ
          IRELREC=IRECIJ-LVLSTART
          IF(IRELREC.EQ.IJ*NSEQ) THEN
            IPOS=LVLSTART+(NPAIRS+IJ-1)*NSEQ+1
          ELSE
            IPOS=IRECIJ+1
          END IF
        END IF
! GET THE POSITION IN THE FILE
        if (byt8) then
        WRITE(NDISK,REC=IPOS) BIN8,INDXBIN,IRECIJ
        else
        WRITE(NDISK,REC=IPOS) IBIN4,IBIN1,INDXBIN,IRECIJ
        endif
        IRECIJ=IPOS
      END
!====================================================================
      subroutine CompactCoef_vir(ncf,nval,coef,nrow,irow,xtrunc,nmo,nbf)
      implicit real*8(a-h,o-z)
      dimension coef(ncf,*),xtrunc(nrow,nval),irow(nrow)
c  This routine prepares a submatrix of the SCF coefficient matrix
c  coef(ncf,nval) according to the index array irow.
c  the resulting matrix (nrow x nval) contains only the rows in
c  irow(k)
      do i=1,nbf-nmo
        do mu1=1,nrow
          mu=irow(mu1)
          xtrunc(mu1,i)=coef(mu,i)
        end do
      end do
      end
c======================================================================
!====================================================================
      subroutine CompactCoef_tj(ncf,nval,coef,nrow,irow,xtrunc)
      implicit real*8(a-h,o-z)
      dimension coef(ncf,nval),xtrunc(nrow,nval),irow(nrow)
c  This routine prepares a submatrix of the SCF coefficient matrix
c  coef(ncf,nval) according to the index array irow.
c  the resulting matrix (nrow x nval) contains only the rows in
c  irow(k)
      do i=1,nval
        do mu1=1,nrow
          mu=irow(mu1)
          xtrunc(mu1,i)=coef(mu,i)
        end do
      end do
      end
c======================================================================
      subroutine PutTrunc_tj(ncf,nrow,ncol,xmat,ymat)
      implicit real*8(a-h,o-z)
      dimension xmat(ncf,ncf),ymat(nrow,ncol)
      do j=1,ncol
        do i=1,nrow
          ymat(i,j)=xmat(i,j)
        end do
      end do
      end
c======================================================================
      subroutine ext3_IntTransp(ncf,   nval,   mu,     lam,    iprnt,
     1                     ndisk,  nrow,    ncol,    irow,    icol,
     2                     thresh, xmat, coefoccu, coefvirt, intint,
     3                     int1, itrunc, halftra,nrec,nmo,nbf)
c   This routine transforms the integrals for a given mu, lam to MO basis,
c   yielding (mu,j|lam,l) where mu , lam are AO indices and j,l are MO
c   ones. It uses halftra both as a matrix 'halftra' and as an array
c   halftra for storing a temporary transformed matrix
c   The resulting integrals are written to disk as integers
c   ARGUMENTS:
c   INPUT
c  ncf: number of contracted basis functions
c  nval: number of orbitals to be transformed (valence orbitals)
c  mu,lam: fixed basis function indices. The transformed integrals are
c          (mu,i|lam,j) where i,j are MO indices
c  iprnt: print level
c  ndisk: the unit number of disk to write on
c  xmat = xmat(nu.isig)=(mu.nu | lam,isig) originally. However,
c         it is compacted here: xmat(p,q)=(mu,irow(p)|lam,icol(q))
c  note that xmat is also accessed by the matrix 'xmat'
c  coef(ncf,nval): SCF coefficients
c  nrow, ncol: the number of nonzero rows and columns of xmat
c  irow, icol = the indices of nonzero rows + columns
c  thresh = integral threshold
c  STORAGE
c  intint: storage for nval**2 4-byte integers
c  int1:  storage for nval**2 1-byte integers
c  itrunc = address of an (ncf x  nval) array. This holds a matrix
c  of the truncated SCF coefficients
c  OUTPUT
c  halftra: half-transformed integrals
c  nrec: number of records written so far (input/output)
c
      use memory
      implicit real*8 (a-h,o-z)
      integer*4 intint(nval,nval)
      integer*1 int1(nval,nval),i1
      dimension halftra(nval,nval),irow(nrow),icol(ncol)
      dimension xmat(ncf,ncf),coefoccu(ncf,nval),coefvirt(ncf,*)
c     common /big/bl(30000)
      parameter (one=1.0d0)
c  defined in the calling program
c  xmat has the COMPACTED AO coulomb matrix
c    prepare the compacted SCF coefficients for the left side
      call CompactCoef_vir(ncf,nval,coefvirt,nrow,irow,bl(itrunc),nmo,
     *                     nbf)
c  perform first matrix multiplication. Note that the matrix
c  system routine matsimtr cannot be used here because the right and
c  left matrices may be different after compaction
c  define temporary matrix product
      iout = 6     ! do NOT change this   JB
      call matdef('temp','r',nbf-nmo,ncol)
c   'ymat' is similar to 'xmat' but is compacted
      call matdef('ymat','r',nrow,ncol)
      call PutTrunc_tj(ncf,nrow,ncol,xmat,bl(mataddr('ymat')))
c  connect the truncated SCF coefficient matrix to the matrix system
      call matconn('trunc','r',nrow,nbf-nmo,itrunc)
      call matmmul2('trunc','ymat','temp','t','n','n')
      call matdisc('trunc')
      call matrem('ymat')
      call CompactCoef_tj(ncf,nval,coefoccu,ncol,icol,bl(itrunc))
      call matconn('trunc','r',ncol,nval,itrunc)
      call matmmult('temp','trunc','halftra')
      call matdisc('trunc')
      call matrem('temp')
c
      if(iprnt.gt.4) then
        call matscal('halftra',thresh)
        print *, mu,lam,nrow,ncol
        call matprint('halftra',iout)
        call matscal('halftra',one/thresh)
      end if
      return
      end
c==========================================================================

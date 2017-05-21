#include "maxi.h"
      subroutine sort_amplit(nval,idimen,ndisk,af,igid,
     *                       nrep,vtable,m_table,ijtable,vrevtable,
     *                       c_ratio)
      use memory
      use kinds
      implicit none
c
      integer nval,idimen,ndisk,igid
      integer nrep,vtable(*),m_table(*),ijtable(*),vrevtable(*)
      logical af
      real*8 c_ratio
c
c
      integer memory1,irecord,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,info,nrec
      character*256 scrfile,filename
      logical fileopen
c
      call matmark()
      call dynamic_mmark()
c
      call dynamic_getmem(idimen,irecord)
      call izeroit(bl(irecord),idimen*intsize)
      call dynamic_getmem(idimen,ifill_level)
      call izeroit(bl(ifill_level),idimen*intsize)
      call dynamic_show_free(memory1)
c
c I need memory for bins, and indexes (3 indexes)
c For given virtual "a", all the bins will occuppy the amount of memory:
c idimen*lbin*(8+3*2)
c
      memory1=memory1-300000 ! possible fragmentation
      lbin=8*memory1/(idimen*14)
      if (lbin.gt.nval*nval*idimen) lbin=nval*nval*idimen
      no_bins=nval*nval*idimen/lbin
      if (mod(nval*nval*idimen,lbin).ne.0) no_bins=no_bins+1
      call dynamic_mmark()
      call dynamic_getmem(idimen*lbin,ibin)
      call dynamic_getmem(6*idimen*lbin/8+1,iindex)
      call izeroit(bl(iindex),(6*idimen*lbin/8+1)*intsize)
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
      filename=scrfile(1:len)//'.sorttmp'
      len1=len+8
      ndisk_sort=50
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk_sort,lbin*14+8)
      else
        inquire(ndisk_sort,opened=fileopen)
        if (fileopen) STOP 'Unit 50 opened'
        OPEN(UNIT=ndisk_sort,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lbin*14+8/intsize)
      endif
c
      nrec=0
      call amplit_sort(bl(ibin),bl(iindex),nval,bl(ifill_level),idimen,
     *                 bl(irecord),ndisk_sort,lbin,af,igid,nrec)
      call leftovers(idimen,   bl(ifill_level),lbin,    bl(iindex),nrec,
     *             bl(irecord),ndisk_sort,     bl(ibin),af,        igid,
     *             idimen,3)
ccs      call sleep(5)
c
      filename=scrfile(1:len)//'.Ta'
      len1=len+3
      ndisk=51
      if (af) then
      call fafCreatem(filename(1:len1)//char(0),ndisk,
     *                nval*nval*idimen*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 51 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*nval*idimen*8)
      endif
      call dynamic_retmark()
c
      call dynamic_mmark()
      call dynamic_getmem(lbin,ibin)
      call dynamic_getmem(6*lbin/8+1,iindex)
      call dynamic_getmem(nval*nval*idimen,ixmat)
      call collect_bins(lbin,bl(ibin),bl(iindex),ndisk_sort,
     *                  ndisk,bl(irecord),idimen,bl(ixmat),nval,af,
     *                  nrep,vtable,m_table,ijtable,vrevtable,c_ratio)
      if (af) then
        call fafClosem(ndisk_sort,0,info)
      else
        close(ndisk_sort,status='delete')
      endif
      call dynamic_retmark()
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine amplit_sort(bin,iindex,nval,ifill_level,idimen,irecord,
     *                       ndisk_bin,lbin,af,igid,nrec)
      use memory
      implicit none
      integer lbin,nval,idimen,igid,nrec
      real*8 bin(lbin,idimen)
      integer*2 iindex(3,lbin,idimen)
      integer ifill_level(idimen),irecord(idimen),ndisk_bin
      logical af
c
      integer i,j,ia,ib,ix,ic_ind,itmp
c
      do i=1,nval
        do j=1,nval
          call CoefRead('tt',i,j,itmp)
          ic_ind=-1
          do ib=1,idimen
            do ia=1,idimen
             ic_ind=ic_ind+1
             ix=ifill_level(ia)+1
             ifill_level(ia)=ix
c
             bin(ix,ia)=bl(itmp+ic_ind)
             iindex(1,ix,ia)=i
             iindex(2,ix,ia)=j
             iindex(3,ix,ia)=ib
c
             if (ix.ge.lbin) then
             if (ix.gt.lbin) STOP 'Error ddfdfdfdf'
               ifill_level(ia)=0
               call amp_bin_wrt(nrec,ia,irecord(ia),ndisk_bin,
     *                          lbin,bin(1,ia),iindex(1,1,ia),af,igid,3)
             endif
            enddo
          enddo
        enddo
      enddo
      end
c===============================================================================
      subroutine leftovers(idimen,  ifill_level, lbin, iindex, nrec,
     *                     irecord, ndisk_bin,   bin,  af,     igid,
     *                     no_bins,no_ind)
      implicit none
      integer idimen,ifill_level(*),lbin,nrec,irecord(*),ndisk_bin,igid
      integer no_bins,no_ind
      integer*2 iindex(no_ind,lbin,idimen)
      real*8 bin(lbin,idimen)
      logical af
c
      integer ia,i
c
      do ia=1,no_bins
        do i=ifill_level(ia)+1,lbin
          iindex(1,i,ia)=0
        enddo
        call amp_bin_wrt(nrec,ia,irecord(ia),ndisk_bin,
     *                   lbin,bin(1,ia),iindex(1,1,ia),af,igid,no_ind)
      enddo
      end
c===============================================================================
      subroutine leftoversgr(idimen,  ifill_level, lbin, iindex, nrec,
     *                       irecord, ndisk_bin,   bin,  af,     igid,
     *                       no_bins,no_ind,igranules,igranulesize,
     *                       igran2pair,ipair2gran)
      implicit none
      integer idimen,ifill_level(*),lbin,nrec,irecord(*),ndisk_bin,igid
      integer no_bins,no_ind,igranules,igranulesize
      integer*2 iindex(no_ind,lbin,idimen)
      integer*4 igran2pair(2,*),ipair2gran(*)
      real*8 bin(lbin,idimen)
      logical af
c
      integer ia,i,istart,istop,igr
c
      do igr=1,igranules
        istart=igran2pair(1,igr)
        istop =igran2pair(2,igr)
        do ia=istart,istop
          do i=ifill_level(ia)+1,lbin
            iindex(1,i,ia)=0
          enddo
        enddo
        call amp_bin_wrt(nrec,igr,irecord(igr),ndisk_bin,
     *lbin*igranulesize,bin(1,istart),iindex(1,1,istart),af,igid,no_ind)
      enddo
      end
c===============================================================================
      subroutine amp_bin_wrt(nrec, ibin,   irecord, ndisk, lbin,
     *                       bin,  iindex, af,      igid,no_ind)
      use memory
      implicit none
      integer nrec,irecord,lbin,ndisk,igid,ibin,no_ind
      real*8 bin(lbin)
      integer*2 iindex(no_ind,lbin)
      logical af
      integer info,indexx
c
      nrec=nrec+1
      if (af) then
        call packer(3,indexx,bl(1),bin,lbin*8,iindex,2*lbin*no_ind,
     1              irecord,8)
        call fafwrite(ndisk,bl(indexx),1,(8+no_ind*2)*lbin+8,1,nrec,
     1                info)
        irecord=nrec
      else
        write(ndisk,rec=nrec) irecord,bin,iindex
        irecord=nrec
      endif
      end
c===============================================================================
      subroutine collect_bins(lbin,bin,iindex,ndisk_s,ndisk_t,irecord,
     *                        idimen,xmat,nval,af,nrep,vtable,m_table,
     *                        ijtable,vrevtable,c_ratio)
      use memory
      implicit none
      integer lbin,ndisk_s,ndisk_t,idimen,nval,nrep,vtable(idimen+1,*)
      integer m_table(8,8),ijtable(nval*nval+1,8),vrevtable(*)
      real*8 bin(lbin),xmat(nval,nval,idimen),c_ratio
      integer*2 iindex(3,lbin)
      integer irecord(idimen)
      logical af
c
      integer ia,nrec,nrec1,ii,i,j,ib,ind(4),igid,info,iscr
      integer id_rep,id_size,ijrep,ijsize,iaddr,id,ids,ij,ijs,icount
      real*8 local_ratio
      integer isize,irec,iindexx
c
      c_ratio=0d0
      call dynamic_mmark()
      call dynamic_getmem(nval*nval*idimen,iscr)
      if (af) then
      do ia=1,idimen
        call zeroit(xmat,nval*nval*idimen)
          irec=irecord(ia)
          do
            if (irec.eq.0) exit
            call bufpointer(iindexx,bl(1),lbin*14+8)
            call fafread(ndisk_s,bl(iindexx),1,lbin*14+8,1,irec,info)
            if (info.lt.lbin*14+8) then
              ERROR
              call flush(6)
            endif
            call unpacker(3,bin,lbin*8,iindex,6*lbin,irec,8)
            if (irec.lt.0) STOP 'Error fdjh5e39w3o9thfrh'
            do ii=1,lbin
              i= iindex(1,ii)
              j= iindex(2,ii)
              ib=iindex(3,ii)
              if (i.eq.0) exit
              if (i.lt.0) STOP 'Error 33k3nfn'
              xmat(i,j,ib)=bin(ii)
            enddo
          enddo ! do nrec=1
      if (nrep.gt.1) then ! Symmetry compression
        call sym_comp_a(xmat,bl(iscr),nval,idimen,nrep,ia,m_table,
     *                      ijtable,vtable,vrevtable,isize)
        local_ratio=dble(isize)/dble(idimen*nval*nval)
        if (local_ratio.gt.c_ratio) c_ratio=local_ratio
        call fafwrite(ndisk_t,bl(iscr),8,isize,1,ia,info)
      else
        c_ratio=1d0
        call fafwrite(ndisk_t,xmat,    8,idimen*nval*nval,1,ia,info)
      endif
      enddo
      else
      do ia=1,idimen
        call zeroit(xmat,nval*nval*idimen)
        nrec=irecord(ia)
        do
          read(ndisk_s,rec=nrec) nrec1,bin,iindex
          do ii=1,lbin
            i= iindex(1,ii)
            j= iindex(2,ii)
            ib=iindex(3,ii)
            if (i.eq.0) exit
            if (i.lt.0) STOP 'Error 33k3nfn'
            xmat(i,j,ib)=bin(ii)
          enddo
          nrec=nrec1
          if (nrec.eq.0) exit
        enddo
      if (nrep.gt.1) then
        call sym_comp_a(xmat,bl(iscr),nval,idimen,nrep,ia,m_table,
     *                      ijtable,vtable,vrevtable,isize)
        call general_write(ndisk_t,ia,bl(iscr),idimen*nval*nval*8)
      else
        write(ndisk_t,rec=ia) xmat
      endif
        c_ratio=1d0
      enddo
      endif
      call dynamic_retmark()
      end
c===============================================================================
c Integrals - 3ext J(d,i)(a,b) type -> K(b,d)(a,i) type
c===============================================================================
      subroutine sort_3ext(nval,idimen,ndisk,irecadrtt,npairs,ndisktt,
     *                     ncf,lbintt,thresh,byt8,nmo,vorb,af,igid,
     *                      nrep,vtable,m_table,vrevtable,otable)
      use memory
      use kinds
      implicit none
c
      integer nval,idimen,ndisk,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nrep,vtable(*),m_table(*),vrevtable(*),otable(*)
      integer nmo,igid
      logical byt8,vorb
      real*8 thresh
      logical af
c
c
      integer memory1,irecord,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,info,nrec
      integer ipasssize,npass,ipass,istart,istop,nslv
      character*256 scrfile,filename
      logical fileopen
c
      nslv=0
      call matmark()
      call dynamic_mmark()
c
      call dynamic_getmem(idimen*idimen,irecord)
      call izeroit(bl(irecord),idimen*idimen*intsize)
      call dynamic_getmem(idimen*idimen,ifill_level)
      call izeroit(bl(ifill_level),idimen*idimen*intsize)
      call dynamic_show_free(memory1)
c
c I need memory for bins, and indexes (2 indexes)
c All the bins will occuppy the amount of memory:
c idimen*idimen*lbin*(8+2*2)
c
      memory1=memory1-300000-idimen*idimen ! possible fragmentation
      memory1=memory1-idimen*idimen ! one memory reservation in ext3_sort
      lbin=8*memory1/(idimen*idimen*12)
      if (lbin.gt.nval*idimen) lbin=nval*idimen
      no_bins=nval*idimen/lbin
      if (mod(nval*idimen,lbin).ne.0) no_bins=no_bins+1
      call dynamic_mmark()
      call dynamic_getmem(idimen*idimen*lbin,ibin)
      call dynamic_getmem(idimen*idimen*lbin/2+1,iindex)
      call izeroit(bl(iindex),(idimen*idimen*lbin/2+1)*intsize)
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
      filename=scrfile(1:len)//'.sorttmp'
      len1=len+8
      ndisk_sort=50
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk_sort)
      else
        inquire(ndisk_sort,opened=fileopen)
        if (fileopen) STOP 'Unit 50 opened'
        OPEN(UNIT=ndisk_sort,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lbin*12+8/intsize)
      endif
c
      nrec=0
      ipasssize=12 ! The framework for parallelization
      npass=idimen*nval/ipasssize
      if (mod(idimen*nval,ipasssize).ne.0) npass=npass+1
      do ipass=1,npass
        istart=(ipass-1)*ipasssize+1
        istop= (ipass  )*ipasssize
        if (istart.gt.idimen*nval) STOP 'Error 392jfallaer'
        if (istop .gt.idimen*nval) istop=idimen*nval
      call ext3_sort(bl(ibin),bl(iindex),nval,bl(ifill_level),idimen,
     *               bl(irecord),ndisk_sort,lbin,irecadrtt,npairs,
     *               ndisktt,ncf,lbintt,thresh,byt8,nmo,vorb,af,igid,
     *               nrec,istart,istop)
      enddo
      call leftovers(idimen,   bl(ifill_level),lbin,    bl(iindex),nrec,
     *             bl(irecord),ndisk_sort,     bl(ibin),af,        igid,
     *             idimen*idimen,2)
ccs      call sleep(5)
c
      filename=scrfile(1:len)//'.Kbc'
      len1=len+4
      ndisk=52
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 52 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*8)
      endif
      call dynamic_retmark()
c
      call dynamic_mmark()
      call dynamic_getmem(lbin,ibin)
      call dynamic_getmem(lbin/2+1,iindex)
      call dynamic_getmem(nval*idimen,ixmat)
c
      ipasssize=12 ! The framework for parallelization
      npass=idimen*idimen/ipasssize
      if (mod(idimen*idimen,ipasssize).ne.0) npass=npass+1
      do ipass=1,npass
        istart=(ipass-1)*ipasssize+1
        istop= (ipass  )*ipasssize
        if (istart.gt.idimen*idimen) STOP 'Error 392jfallaer'
        if (istop .gt.idimen*idimen) istop=idimen*idimen
        call coll_ext3_bins(lbin,bl(ibin),bl(iindex),ndisk_sort,
     *                      ndisk,bl(irecord),idimen,bl(ixmat),
     *                      nval,af,istart,istop,nslv,
     *                      nrep,vtable,m_table,vrevtable,otable)
      enddo
      if (af) then
        call fafClosem(ndisk_sort,0,info)
      else
        close(ndisk_sort,status='delete')
      endif
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine ext3_sort(bin,iindex,nval,ifill_level,idimen,irecord,
     *                     ndisk_bin,lbin,irecadrtt,npairs,ndisktt,ncf,
     *                     lbintt,thresh,byt8,nmo,vorb,af,igid,nrec,
     *                     istart,istop)
      use memory
      implicit none
      integer lbin,nval,idimen,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,igid,nrec,istart,istop
      logical byt8,vorb,af
      real*8 bin(lbin,idimen*idimen),thresh
      integer*2 iindex(2,lbin,idimen*idimen)
      integer ifill_level(idimen*idimen),irecord(idimen*idimen)
      integer ndisk_bin
c
      integer i,j,ia,ib,ix,ic_ind,itmp,ibin,id,iid
c
      call dynamic_getmem(idimen*idimen,itmp)
c
      iid=0
      do i=1,nval
        do id=1,idimen
          iid=iid+1
          if (iid.lt.istart) cycle
          if (iid.gt.istop) exit
          call ExtrOne(id,    i,     irecadrtt,  npairs,    ndisktt,
     *                 ncf,   lbintt,thresh,     byt8,      'tt',
     *                 'mo',  nmo,   vorb,       bl(itmp))
          ic_ind=-1
          do ib=1,idimen
            do ia=1,idimen
             ibin=(ib-1)*idimen+id
             ic_ind=ic_ind+1
             ix=ifill_level(ibin)+1
             ifill_level(ibin)=ix
c
             bin(ix,ibin)=bl(itmp+ic_ind)
             iindex(1,ix,ibin)=ia
             iindex(2,ix,ibin)=i
c
             if (ix.ge.lbin) then
             if (ix.gt.lbin) STOP 'Error ddfdfdfdf4433'
               ifill_level(ibin)=0
               call amp_bin_wrt(nrec,ibin,irecord(ibin),ndisk_bin,
     *                      lbin,bin(1,ibin),iindex(1,1,ibin),af,igid,2)
             endif
            enddo
          enddo
        enddo
      enddo
      call dynamic_retmem(1)
      end
c===============================================================================
      subroutine coll_ext3_bins(lbin,bin,iindex,ndisk_s,ndisk_t,irecord,
     *                        idimen,xmat,nval,af,istart,istop,nslv,
     *                        nrep,vtable,m_table,vrevtable,otable)
      use memory
      implicit none
      integer lbin,ndisk_s,ndisk_t,idimen,nval,istart,istop
      integer nrep,vtable(*),m_table(*),vrevtable(*),otable(*)
      real*8 bin(lbin),xmat(idimen,nval)
      integer*2 iindex(2,lbin)
      integer irecord(idimen*idimen,*),nslv
      logical af
c
      integer ia,nrec,nrec1,ii,i,j,ib,iab,id,ind(4),igid,info,iscr
c
      call dynamic_mmark()
      call dynamic_getmem(nval*idimen,iscr)
      if (af) then
      iab=0
      do ia=1,idimen
        do ib=1,idimen
        iab=iab+1
        if (iab.lt.istart) cycle
        if (iab.gt.istop) exit
        call zeroit(xmat,nval*idimen)
        do igid=1,nslv
          do nrec=1,irecord(iab,igid) ! irecord should refer to igid too
            ind(1)=igid
            ind(2)=iab          ! bin no
            ind(3)=nrec
            ind(4)=1
            call fafread(ndisk_s,bin,   8,lbin,  4,ind,info)
          if (info.lt.8*lbin) then
              ERROR
              call flush(6)
          endif
            ind(4)=2
            call fafread(ndisk_s,iindex,2,lbin*2,4,ind,info)
          if (info.lt.2*lbin*2) then
              ERROR
              call flush(6)
          endif
            do ii=1,lbin
              id= iindex(1,ii)
              i = iindex(2,ii)
              if (id.le.0) then
                if (id.lt.0) STOP 'Error 33k3nfn5rh'
                exit
              endif
              xmat(id,i)=bin(ii)
            enddo
          enddo
        enddo
c
        if (nrep.gt.1) then
          call symm_compress(xmat,bl(iscr),nrep,idimen,nval,vtable,
     *                       m_table,otable,vrevtable,ia,ib)
          call fafwrite(ndisk_t,bl(iscr),8,idimen*nval,1,iab,info)
        else
          call fafwrite(ndisk_t,xmat,8,idimen*nval,1,iab,info)
        endif
c
        enddo
      enddo
      else
      iab=0
      do ia=1,idimen
        do ib=1,idimen
        iab=iab+1
        if (iab.lt.istart) cycle
        if (iab.gt.istop) exit
        call zeroit(xmat,nval*idimen)
        nrec=irecord(iab,1)
        do
          read(ndisk_s,rec=nrec) nrec1,bin,iindex
          do ii=1,lbin
            id= iindex(1,ii)
            i = iindex(2,ii)
            if (id.le.0) then
              if (id.lt.0) STOP 'Error 33k3nfn5rh'
              exit
            endif
            xmat(id,i)=bin(ii)
          enddo
          nrec=nrec1
          if (nrec.eq.0) exit
        enddo
        if (nrep.gt.1) then
          call symm_compress(xmat,bl(iscr),nrep,idimen,nval,vtable,
     *                       m_table,otable,vrevtable,ia,ib)
          call general_write(ndisk_t,iab,bl(iscr),idimen*nval*8)
        else
          write(ndisk_t,rec=iab) xmat
        endif
        enddo
      enddo
      endif
      call dynamic_retmark()
      end
c===============================================================================
c Integrals - 3internal  K(i,j,k)(c) type -> J(c)(j,k,i) type
c                          k i j                  i j k  
c===============================================================================
      subroutine sort_3int(nval,idimen,ndisk,af,igid)
      use memory
      use kinds
      implicit none
c
      integer nval,idimen,ndisk,igid
      logical af
c
c
      integer memory1,irecord,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,info,nrec
      character*256 scrfile,filename
      logical fileopen
c
      call matmark()
      call dynamic_mmark()
c
      call dynamic_getmem(idimen,irecord)
      call izeroit(bl(irecord),idimen*intsize)
      call dynamic_getmem(idimen,ifill_level)
      call izeroit(bl(ifill_level),idimen*intsize)
      call dynamic_show_free(memory1)
c
c I need memory for bins, and indexes (3 indexes)
c For given virtual "a", all the bins will occuppy the amount of memory:
c idimen*lbin*(8+3*2)
c
      memory1=memory1-300000 ! possible fragmentation
      lbin=8*memory1/(idimen*14)
      if (lbin.gt.nval*nval*nval) lbin=nval*nval*nval
      no_bins=nval*nval*nval/lbin
      if (mod(nval*nval*nval,lbin).ne.0) no_bins=no_bins+1
      call dynamic_mmark()
      call dynamic_getmem(idimen*lbin,ibin)
      call dynamic_getmem(6*idimen*lbin/8+1,iindex)
      call izeroit(bl(iindex),(6*idimen*lbin/8+1)*intsize)
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
      filename=scrfile(1:len)//'.sorttmp'
      len1=len+8
      ndisk_sort=50
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk_sort,lbin*14+8)
      else
        inquire(ndisk_sort,opened=fileopen)
        if (fileopen) STOP 'Unit 50 opened'
        OPEN(UNIT=ndisk_sort,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lbin*14+8/intsize)
      endif
c
      nrec=0
      call int3_sort(bl(ibin),bl(iindex),nval,bl(ifill_level),idimen,
     *               bl(irecord),ndisk_sort,lbin,af,igid,nrec)
      call leftovers(idimen,   bl(ifill_level),lbin,    bl(iindex),nrec,
     *             bl(irecord),ndisk_sort,     bl(ibin),af,        igid,
     *             idimen,3)
ccs      call sleep(5)
c
      filename=scrfile(1:len)//'.Jc'
      len1=len+3
      ndisk=53
      if (af) then
       call fafCreatem(filename(1:len1)//char(0),ndisk,nval*nval*nval*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 53 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*nval*nval*8)
      endif
      call dynamic_retmark()
c
      call dynamic_mmark()
      call dynamic_getmem(lbin,ibin)
      call dynamic_getmem(6*lbin/8+1,iindex)
      call dynamic_getmem(nval*nval*nval,ixmat)
      call coll_int3_bins(lbin,bl(ibin),bl(iindex),ndisk_sort,
     *                    ndisk,bl(irecord),idimen,bl(ixmat),nval,af)
      if (af) then
        call fafClosem(ndisk_sort,0,info)
      else
        close(ndisk_sort,status='delete')
      endif
      call dynamic_retmark()
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine int3_sort(bin,iindex,nval,ifill_level,idimen,irecord,
     *                       ndisk_bin,lbin,af,igid,nrec)
      use memory
      implicit none
      integer lbin,nval,idimen,nrec
      real*8 bin(lbin,idimen)
      integer*2 iindex(3,lbin,idimen)
      integer ifill_level(idimen),irecord(idimen),ndisk_bin,igid
      logical af
c
      integer i,j,k,ia,ix,itmp
c
      call dynamic_getmem(idimen,itmp)
      do i=1,nval
        do j=1,nval
          do k=1,nval
          call xkijk(i,j,k,bl(itmp))
          do ia=1,idimen
             ix=ifill_level(ia)+1
             ifill_level(ia)=ix
c
             bin(ix,ia)=bl(itmp+ia-1)
             iindex(1,ix,ia)=j
             iindex(2,ix,ia)=k
             iindex(3,ix,ia)=i
c
             if (ix.ge.lbin) then
             if (ix.gt.lbin) STOP 'Error ddfdfdfdf'
               ifill_level(ia)=0
               call amp_bin_wrt(nrec,ia,irecord(ia),ndisk_bin,lbin,
     *                          bin(1,ia),iindex(1,1,ia),af,igid,3)
             endif
            enddo
          enddo
        enddo
      enddo
      call dynamic_retmem(1)
      end
c===============================================================================
      subroutine coll_int3_bins(lbin,bin,iindex,ndisk_s,ndisk_t,irecord,
     *                          idimen,xmat,nval,af)
      use memory
      implicit none
      integer lbin,ndisk_s,ndisk_t,idimen,nval
      real*8 bin(lbin),xmat(nval,nval,nval)
      integer*2 iindex(3,lbin)
      integer irecord(idimen),iindexx,irec
      logical af
c
      integer ia,nrec,nrec1,ii,i,j,k,ind(4),igid,info
c
      if (af) then
      do ia=1,idimen
        call zeroit(xmat,nval*nval*nval)
          irec=irecord(ia)
          do
            if (irec.eq.0) exit
            call bufpointer(iindexx,bl(1),lbin*14+8)
            call fafread(ndisk_s,bl(iindexx),1,lbin*14+8,1,irec,info)
            if (info.lt.lbin*14+8) then
              ERROR
              call flush(6)
            endif
            call unpacker(3,bin,lbin*8,iindex,6*lbin,irec,8)
            if (irec.lt.0) STOP 'Error fdjh5e6533k2thfrh'
            do ii=1,lbin
              i=iindex(1,ii)
              j=iindex(2,ii)
              k=iindex(3,ii)
              if (i.le.0) then
                if (i.lt.0) STOP 'Error 33k3nfn4gh'
                exit
              endif
              xmat(i,j,k)=bin(ii)
            enddo
          enddo
        call fafwrite(ndisk_t,xmat,8,nval*nval*nval,1,ia,info)
      enddo
      else
      do ia=1,idimen
        call zeroit(xmat,nval*nval*nval)
        nrec=irecord(ia)
        do
          read(ndisk_s,rec=nrec) nrec1,bin,iindex
          do ii=1,lbin
            i=iindex(1,ii)
            j=iindex(2,ii)
            k=iindex(3,ii)
            if (i.le.0) then
              if (i.lt.0) STOP 'Error 33k3nfn4gh'
              exit
            endif
            xmat(i,j,k)=bin(ii)
          enddo
          nrec=nrec1
          if (nrec.eq.0) exit
        enddo
        write(ndisk_t,rec=ia) xmat
      enddo
      endif
      end
c===============================================================================
      subroutine new_triples(nval,idimen,ndisk_a,ndisk_ie,ndisk_ii,
     *                       ndisk_ix,energy,eps,etimes,n,ccsd,
     *                       qcisd,af,otable,vtable,vrevtable,m_table,
     *                       ijtable,nrep,npass,isize,energys,energyd,
     *                       iprnt)
      use memory
      implicit none
c
      integer nval,idimen,ndisk_a,ndisk_ie,ndisk_ii,ndisk_ix,n,nrep
      integer otable(*),vtable(*),vrevtable(*),m_table(*),ijtable(*)
      integer npass,isize,iprnt
      real*8 energy,eps(*),etimes(n),energys,energyd,enes,ened
      logical ccsd,qcisd,af
c
      integer iWabc,ia,ib,ic,i,irec_amp,n_ampl,iamplit,ip,ip_dim
      integer irec_int,n_int,i3int_int,irec_ext,n_ext,i3ext_int,isingles
      integer iap,ibp,icp,ia1,ia2,ib1,ib2,ic1,ic2,hb,hc
      integer iKab,iKac,iKbc,irec,info,iK,irec_K,n_K,nvirpairs
      real*8 en,t0,t1,x0,x1,amp_ratio,ext3_ratio
      logical reset1,reset2,reset3,reset4
      integer ip_array,fragmentation,iampl_store,i3int_store,i3ext_store
      integer iK_store,ialloc_amount,mem,ireal_size,ir_ampl,nslv
c
      call matmark()
      call dynamic_mmark()
c
c etimes: 1 the time spend in W zeroing and sorting
c         2 W building, total time
c         3 in W: multiplication - over virt
c         4 in W: multiplication - over occ
c         5 reading
c         6 reading locate amplitudes
c         7 reading locate 3int integrals
c         8 reading locate 3ext integrals
      energy=0d0
      energys=0d0
      energyd=0d0
      do i=1,n
        etimes(i)=0d0
      enddo
      nvirpairs=idimen*(idimen+1)/2
c
      call dynamic_getmem(nval*nval*nval,iWabc)
      if (qcisd.or.ccsd) then
        call dynamic_getmem(nval*nval,iKab)
        call dynamic_getmem(nval*nval,iKac)
        call dynamic_getmem(nval*nval,iKbc)
      endif
c MEMORY MANAGEMENT:
      if (iprnt.ge.3) then
      write(6,*) 'calculated chunk size, no of passes: ',isize,npass
      endif
c END OF MEMORY MANAGEMENT
c
      ip_dim=6*isize*isize  ! is it woth to reduce to integer*4?
      call dynamic_getmem(ip_dim,ip)
c locate_ampl setup:
      call dynamic_getmem(idimen,irec_amp)
      reset1=.true.
      n_ampl=isize*3
      ir_ampl=nval*nval*idimen
      call dynamic_getmem(nval*nval*idimen*n_ampl,iamplit)
c locate_3int setup:
      call dynamic_getmem(idimen,irec_int)
      reset2=.true.
      n_int=isize*3
      call dynamic_getmem(nval*nval*nval*n_int,i3int_int)
c locate_3ext setup:
      call dynamic_getmem(nvirpairs,irec_ext)
      reset3=.true.
      n_ext=3*isize*isize
      ireal_size=nval*idimen
      call dynamic_getmem(2*ireal_size*n_ext,i3ext_int)
c locate_K setup:
      call dynamic_getmem(nvirpairs,irec_K)
      reset4=.true.
      n_K=6*isize*isize
      call dynamic_getmem(nval*nval*n_K,iK)
c
      call pointersingles(isingles)
c
      do iap=1,npass
        ia1=(iap-1)*isize+1
        ia2=(iap  )*isize
        if (ia1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
        if (ia2.gt.idimen) ia2=idimen
        do ibp=1,iap
          ib1=(ibp-1)*isize+1
          ib2=(ibp  )*isize
          if (ib1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
          if (ib2.gt.idimen) ib2=idimen
          do icp=1,ibp
            ic1=(icp-1)*isize+1
            ic2=(icp  )*isize
            if (ic1.gt.idimen) STOP 'Error 054h32jvhsjdw;a'
            if (ic2.gt.idimen) ic2=idimen
            call elapsec(x0)
            call elapsec(t0)
            call locate_ampl(ia1,   ia2,    ib1,         ib2,   ic1,
     *                       ic2,   iamplit,bl(irec_amp),reset1,n_ampl,
     *                       idimen,nval,   ndisk_a,     af,    m_table,
     *                       ijtable,vtable,   vrevtable,nrep,ir_ampl)
            call elapsec(x1)
            etimes(6)=etimes(6)+x1-x0
            call elapsec(x0)
            call locate_3int(ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   i3int_int,bl(irec_int),reset2,n_int,
     *                       idimen,nval,     ndisk_ii,    af)
            call elapsec(x1)
            etimes(7)=etimes(7)+x1-x0
            call elapsec(x0)
            call locate_3ext(ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   i3ext_int,bl(irec_ext),reset3,n_ext,
     *                       idimen,nval,     ndisk_ie,    bl(ip),af,
     *                       nrep,vtable,m_table,otable,vrevtable,
     *                       ireal_size)
            call elapsec(x1)
            etimes(8)=etimes(8)+x1-x0
            call elapsec(x0)
            if (qcisd.or.ccsd) then
            call locate_K   (ia1,   ia2,      ib1,         ib2,   ic1,
     *                       ic2,   iK,       bl(irec_K),  reset4,n_K,
     *                       idimen,nval,     ndisk_ix,    bl(ip),af,
     *                       isize, npass,    iap,         ibp,   icp)
            endif
            call elapsec(x1)
            etimes(10)=etimes(10)+x1-x0
            call elapsec(t1)
            etimes(5)=etimes(5)+t1-t0
c
c
      do ia=ia1,ia2
        hb=min(ia,ib2)
        do ib=ib1,hb
          hc=min(ib,ic2)
          do ic=ic1,hc
            if (ia.eq.ib.and.ib.eq.ic) cycle
c           print *,ia,ib,ic
            if (qcisd.or.ccsd) then
c             irec=ia*(ia-1)/2+ib
c             if (af) then
c             call fafread(ndisk_ix,bl(iKab),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKab),nval*nval*8)
c             endif
c             irec=ia*(ia-1)/2+ic
c             if (af) then
c             call fafread(ndisk_ix,bl(iKac),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKac),nval*nval*8)
c             endif
c             irec=ib*(ib-1)/2+ic
c             if (af) then
c             call fafread(ndisk_ix,bl(iKbc),8,nval*nval,1,irec,info)
c             else
c             call general_read(ndisk_ix,irec,bl(iKbc),nval*nval*8)
c             endif
              call return_addresses(ia,ib,ic,bl(irec_K),iKab,iKac,iKbc)
            endif
            call elapsec(t0)
            call zeroit(bl(iWabc),nval*nval*nval)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ia,      ib,       ic,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ia,      ic,       ib,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ic,      ia,       ib,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ic,      ib,       ia,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ib,      ic,       ia,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,2,3)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Wabc_new(ib,      ia,       ic,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,bl(irec_amp),bl(irec_int),bl(irec_ext),
     *                    otable,vtable,vrevtable,m_table,ijtable,nrep,
     *                    af)
            call elapsec(t1)
            etimes(2)=etimes(2)+t1-t0
            call elapsec(t0)
            call sorter(bl(iWabc),nval,1,2)
            call elapsec(t1)
            etimes(1)=etimes(1)+t1-t0
            call elapsec(t0)
            call Energy_fromW(ia,  ib, ic,          nval,     bl(iWabc),
     *                        eps, en, bl(isingles),ccsd,     qcisd,
     *                        bl(iKab),bl(iKac),bl(iKbc),idimen,enes,
     *                        ened)
            call elapsec(t1)
            etimes(9)=etimes(9)+t1-t0
            energy=energy+en
            energys=energys+enes
            energyd=energyd+ened
          enddo
        enddo
      enddo
          enddo
        enddo
      enddo
c
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine Wabc_new(ia,      ib,       ic,nval,  idimen,
     *                    ndisk_a, ndisk_ie, ndisk_ii,  iWabc,
     *                    etimes,irec_amp,irec_int,irec_ext,
     *                    otable,vtable,vrevtable,m_table,ijtable,
     *                    nrep,af)
      use memory
      implicit none
      integer ia,ib,ic,nval,idimen,ndisk_a,ndisk_ie,ndisk_ii,iWabc
      integer irec_amp(idimen),irec_int(idimen),irec_ext(idimen*idimen)
      logical af
      real*8 etimes(*)
      integer otable(nval+1,*),vtable(idimen+1,*),vrevtable(*)
      integer m_table(8,8),ijtable(nval*nval+1,*),nrep
c
      integer icoef_base,iinte,mataddr,irec,k,inte_base,iWij,ii,kk
      real*8 t0,t1,xt0,xt1
      integer id_rep,id_size,irep,ik_size,isinte,id,ids,ik,iks,iscoef,ij
      integer iWABCs,ijsize,ijs,ijrep,iks_1_nval2,ik_1_ijsize
      integer ico,il_rep,l_size,i_size,iTsm,l,i,ls,is,iIsm,jk,isize
      integer jk_rep,jks,jk_size,iacount,iicount,icoe,iaddr,iWW,nval2
      logical done
c
      call matmark()
      call dynamic_mmark()
      call mmark()
      nval2=nval*nval
c
      call matconn('WABC','r',nval*nval,nval,iWabc)
      icoef_base=irec_amp(ia)
      if (icoef_base.eq.1) STOP 'Error 45jerkfhr'
      call matconn('coeff','r',nval*nval,idimen,icoef_base)
      if (ib.ge.ic) then
        irec=(ib-1)*ib/2+ic
        iinte=irec_ext(irec)
      else
        irec=(ic-1)*ic/2+ib
        iinte=irec_ext(irec)
        if (af) then
          call ext3_sizer(nrep,idimen,nval,vtable,m_table,
     *                    otable,vrevtable,ib,ic,isize)
          iinte=iinte+isize
        else
          iinte=iinte+nval*idimen
        endif
      endif
      if (iinte.eq.1) STOP 'Error 4k509dfhj4096'
c
c  Begin virt
c
#if 1
      call elapsec(xt0)
      if (nrep.gt.1) then ! symmetry
      iacount=0
      iicount=0
      do id_rep=1,nrep  ! loop over virtual representations
      call elapsec(t0)
        id_size=vtable(1,id_rep)
        irep=m_table(id_rep,vrevtable(ib))
        irep=m_table(irep,vrevtable(ic))
        ijrep=m_table(id_rep,vrevtable(ia))
        ijsize=ijtable(1,ijrep)
        ik_size =otable(1,irep)
      call elapsec(t1)
      etimes(15)=etimes(15)+t1-t0
        if (id_size.eq.0.or.ijsize.eq.0.or.ik_size.eq.0) then 
          iicount=iicount+id_size*ik_size
          iacount=iacount+id_size*ijsize
          cycle
        endif
#if 0
        call matdef('sinte','r',id_size,ik_size)
        isinte=mataddr('sinte')
      call elapsec(t0)
        do ik=1,ik_size
          do id=1,id_size
            ids=vtable(id+1,id_rep)
            iks=otable(ik+1,irep)
         bl(isinte+id-1+(ik-1)*id_size)=bl(iinte+(ids-1)+(iks-1)*idimen)
          enddo
        enddo
      call elapsec(t1)
      etimes(12)=etimes(12)+t1-t0
#else
      call matconn('sinte','r',id_size,ik_size,iinte+iicount)
      iicount=iicount+id_size*ik_size
#endif
#if 0
        call matdef('scoef','r',ijsize,id_size)
        iscoef=mataddr('scoef')
        do id=1,id_size
          ids=vtable(id+1,id_rep)
          do ij=1,ijsize
            ijs=ijtable(ij+1,ijrep)
            bl(iscoef+(ij-1)+(id-1)*ijsize)=
     *                      bl(icoef_base+(ijs-1)+(ids-1)*nval*nval)
          enddo
        enddo
#else
        call matconn('scoef','r',ijsize,id_size,icoef_base+iacount)
        iacount=iacount+id_size*ijsize
#endif
      call elapsec(t0)
      call matdef('WABCs','r',ijsize,ik_size)
      iWABCs=mataddr('WABCs')
      call matmmul2('scoef','sinte','WABCs','n','n','n')
      call elapsec(t1)
      etimes(3)=etimes(3)+t1-t0
c relocate result to WABC
      call elapsec(t0)
      do ik=1,ik_size
        iks=otable(ik+1,irep)
        iks_1_nval2=(iks-1)*nval2
        ik_1_ijsize=(ik-1)*ijsize
        do ij=1,ijsize
          ijs=ijtable(ij+1,ijrep)
          bl(iWabc+(ijs-1)+iks_1_nval2)=
     *    bl(iWabc+(ijs-1)+iks_1_nval2)+
     *    bl(iWABCs+(ij-1)+ik_1_ijsize)
        enddo
      enddo
      call elapsec(t1)
      etimes(11)=etimes(11)+t1-t0
      call matrem('WABCs')
      call matdisc('scoef')
      call matdisc('sinte')
      enddo ! over d representations
      else ! no symmetry:
      call matconn('inte','r',idimen,nval,iinte)
      call elapsec(t0)
      call matmmul2('coeff','inte','WABC','n','n','a')
      call elapsec(t1)
      etimes(3)=etimes(3)+t1-t0
      call matdisc('inte')
      call matdisc('coeff')
      call matdisc('WABC')
      endif
      call elapsec(xt1)
      etimes(14)=etimes(14)+xt1-xt0
#endif 
c
c End of virt part
c
      inte_base=irec_int(ic)
      if (inte_base.eq.1) STOP 'erj59009fgjew'
      ico=icoef_base+(ib-1)*nval*nval
# if 0
      if (nrep.gt.1) then ! symmetry
      do il_rep=1,nrep  ! loop over virtual representations
        l_size=otable(1,il_rep)
        irep=m_table(il_rep,vrevtable(ia))
        irep=m_table(irep,vrevtable(ib))
        jk_rep=m_table(il_rep,vrevtable(ic))
        i_size=otable(1,irep)
        jk_size=ijtable(1,jk_rep)
        if (l_size.eq.0.or.i_size.eq.0.or.jk_size.eq.0) cycle
        call matdef('Tsm','r',i_size,l_size)
        iTsm=mataddr('Tsm')
        do l=1,l_size
          ls=otable(l+1,il_rep)
          do i=1,i_size
            is=otable(i+1,irep)
            bl(iTsm+(i-1)+(l-1)*i_size)=bl(ico+(is-1)+(ls-1)*nval)
          enddo
        enddo
        call matdef('Ism','r',l_size,jk_size)
        iIsm=mataddr('Ism')
        do jk=1,jk_size
          jks=ijtable(jk+1,jk_rep)
          do l=1,l_size
            ls=otable(l+1,il_rep)
            bl(iIsm+(l-1)+(jk-1)*l_size)=
     *                            bl(inte_base+(ls-1)+(jks-1)*nval)
          enddo
        enddo
        call matdef('WABCs','r',i_size,jk_size)
        iWABCs=mataddr('WABCs')
        call elapsec(t0)
        call matmmul2('Tsm','Ism','WABCs','n','n','n')
        call elapsec(t1)
        etimes(4)=etimes(4)+t1-t0
        do jk=1,jk_size
          jks=ijtable(jk+1,jk_rep)
          do i=1,i_size
            is=otable(i+1,irep)
            bl(iWabc+(is-1)+(jks-1)*nval)=
     *      bl(iWabc+(is-1)+(jks-1)*nval)-
     *      bl(iWABCs+(i-1)+(jk-1)*i_size)
          enddo
        enddo
        call matrem('WABCs')
        call matrem('Ism')
        call matrem('Tsm')
      enddo
      else
#endif
      call elapsec(t0)
      if (nrep.gt.1) then ! symmetry
        call matdef('coef1','q',nval,nval)
        call matzero('coef1')
        icoe=mataddr('coef1')
        iicount=0
        done=.false.
        do id_rep=1,nrep  ! loop over virtual representations
          id_size=vtable(1,id_rep)
          ijrep=m_table(id_rep,vrevtable(ia))
          ijsize=ijtable(1,ijrep)
          if (id_size.eq.0.or.ijsize.eq.0) then
            cycle
          endif
          if (vrevtable(ib).eq.id_rep) then
            iaddr=icoef_base+iicount
            do id=1,id_size
              if (vtable(id+1,id_rep).eq.ib) then
                do ij=1,ijsize
                  ijs=ijtable(ij+1,ijrep)
                  bl(icoe+ijs-1)=bl(iaddr+ij-1)
                enddo
                done=.true.
              else
                iaddr=iaddr+ijsize
              endif
            enddo
            exit
          endif
          iicount=iicount+id_size*ijsize
        enddo
cc        if (.not.done) STOP '49823klsdfhwkh4jj'
      else
        call matconn('coef1','q',nval,nval,ico)
      endif
      call elapsec(t1)
      etimes(13)=etimes(13)+t1-t0
      call elapsec(t0)
      call matconn('inte','r',nval,nval*nval,inte_base)
      call matconn('Wij','r',nval,nval*nval,iWabc)
c
      call matdef('WW','r',nval,nval*nval)
      iWW=mataddr('WW')
      call matmmul2('coef1','inte','WW','n','n','n')
c       write(99) ia,ib
c       do ii=0,nval*nval*nval-1
c         write(99) bl(iWW+ii)
c       enddo
      call matadd1('WW',-1d0,'Wij')
c     call matmmul2('coef1','inte','Wij','n','n','s')
      call matrem('WW')
c
      call matdisc('Wij')
      call matdisc('inte')
      call elapsec(t1)
      etimes(4)=etimes(4)+t1-t0
c
      call matremark()
      call dynamic_retmark()
      call retmark()
      end
c===============================================================================
      subroutine Energy_fromW(ia,   ib,     ic,       nval,     W,
     *                        eps,  energy, xsin,     ccsd,     qcisd,
     *                        Kab,  Kac,    Kbc,      idimen,   energys,
     *                        energyd)
      implicit none
      integer nval,ia,ib,ic,ndisk_ix,idimen
      real*8 W(nval,nval,nval),energy,eps(*),Kab(nval,nval)
      real*8 Kac(nval,nval),Kbc(nval,nval),xsin(idimen,nval)
      real*8 energys,energyd,doubles_numerator
      logical ccsd,qcisd
c
      integer i,j,k
      real*8 den_virt,en_i,en_j,en_k,denominator,delta,s
      logical singles
c
      energy=0d0
      energys=0d0
      energyd=0d0
      den_virt=eps(nval+ia)+eps(nval+ib)+eps(nval+ic)
      singles=.false.
      s=0d0
      if (ccsd.or.qcisd) singles=.true.
      do i=1,nval
        en_i=eps(i)
        do j=1,nval
          en_j=eps(j)
          do k=1,nval
            if (singles) then
              s=xsin(ia,i)*Kbc(j,k)+xsin(ib,j)*Kac(i,k)+
     *          xsin(ic,k)*Kab(i,j)
              if (qcisd) s=s*2d0
c             s=xsin(ia,i)*xt(j,ib,k,ic)+
c    *          xsin(ib,j)*xt(k,ic,i,ia)+
c    *          xsin(ic,k)*xt(i,ia,j,ib)
            endif
            en_k=eps(k)
            denominator=en_i+en_j+en_k-den_virt
            doubles_numerator= (
     *      +4d0*W(i,j,k)
     *      -2d0*W(i,k,j)
     *      -2d0*W(k,j,i)
     *      -2d0*W(j,i,k)
     *      +    W(j,k,i)
     *      +    W(k,i,j)          )
            energys=energys+(s)*doubles_numerator/denominator
            energyd=energyd+(W(i,j,k))*doubles_numerator/denominator
          enddo
        enddo
      enddo
      energys=energys*(2d0-delta(ia,ib)-delta(ib,ic))
      energyd=energyd*(2d0-delta(ia,ib)-delta(ib,ic))
      energy=energys+energyd
      end
c===============================================================================
      subroutine locate_ampl(a1,     a2,       b1,       b2,    c1,
     *                       c2,     i3ext_int,irecords, reset, n,  
     *                       idimen, nval,     ndisk,    af,    m_table,
     *                      ijtable,vtable,   vrevtable,nrep,ireal_size)
      use memory
      implicit none
      integer a1,a2,b1,b2,c1,c2,i3ext_int,irecords(*),n,idimen,nval
      integer ndisk,m_table(*),ijtable(*),vtable(*),vrevtable(*)
      integer ireal_size,nrep
      logical reset,af
c
      integer iaddr,iaddr1,irep,irep1
      integer addresses(n),iad_max,ia,ll,ii,id1,id2,id3
      integer ip(3*idimen),i,ip_len,info,isize
      logical shall_I_exit,cycle1
c
      i=0
      do ia=a1,a2
        i=i+1
        ip(i)=ia
      enddo
      do ia=b1,b2
        i=i+1
        ip(i)=ia
      enddo
      do ia=c1,c2
        i=i+1
        ip(i)=ia
      enddo
      ip_len=i
      if (reset) then
        reset=.false.
        do ia=1,n
          iaddr        =i3ext_int+(ia-1)*ireal_size
          addresses(ia)=iaddr
          if (iaddr.eq.1) 
     *  STOP 'Cannot be! I use assumption that the address = 1 is null!'
        enddo
        iad_max=n
        do ia=1,idimen
          irecords(ia)=1
        enddo
c I assume that a2>=a1, b2>=b1, c2>=c1
        if (a2.gt.n.or.b2.gt.n.or.c2.gt.n) then ! load a,b,c first
          do irep=1,ip_len
          if (irecords(ip(irep)).ne.1) cycle
          if (iad_max.le.0) STOP 'Severe error 548gh33'
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ip(irep))=iaddr
          if (af) then
            call ampl_sizer(nval,idimen,nrep,ip(irep),m_table,
     *                      ijtable,vtable,vrevtable,isize)
           if (isize.gt.ireal_size) STOP '50a20ektjtre'
        call fafread(ndisk,bl(iaddr),8,isize,1,ip(irep),info)
          if (info.lt.8*isize) then
              ERROR
              call flush(6)
          endif
          else
          call general_read(ndisk,ip(irep),bl(iaddr),nval*nval*idimen*8)
          endif
          enddo
        endif
        do ia=1,n
          if (iad_max.le.0) exit
          if (irecords(ia).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ia)=iaddr
          if (af) then
            call ampl_sizer(nval,idimen,nrep,ia,m_table,
     *                      ijtable,vtable,vrevtable,isize)
           if (isize.gt.ireal_size) STOP '5020euuktjtre'
            call fafread(ndisk,bl(iaddr),8,isize,1,ia,info)
          if (info.lt.8*isize) then
              ERROR
              call flush(6)
          endif
          else
            call general_read(ndisk,ia,bl(iaddr),nval*nval*idimen*8)
          endif
        enddo
      endif
c
      shall_I_exit=.true.
      z: do irep=1,ip_len
           if (irecords(ip(irep)).eq.1) then
             shall_I_exit=.false.
             exit z
           endif
         enddo z
      if (shall_I_exit) return
c
c
      do irep=1,ip_len
      if (irecords(ip(irep)).eq.1) then
      g: do ia=idimen,1,-1 ! search for largest record which can be replaced
          if (irecords(ia).ne.1) then
            if (ia.eq.ip(irep)) STOP 'Cannot be: ia.eq.i'
cccccccc    if (ia.eq.ind2.or.ia.eq.ind3) cycle
            do ii=1,ip_len
              if (ip(ii).eq.ia) cycle g
            enddo
            iaddr=irecords(ia)
            irecords(ia)=1
            if (irecords(ip(irep)).ne.1) STOP 'irecords(i).ne.1'
            irecords(ip(irep))=iaddr
            if (af) then
            call ampl_sizer(nval,idimen,nrep,ip(irep),m_table,
     *                      ijtable,vtable,vrevtable,isize)
           if (isize.gt.ireal_size) STOP '502xxx0ektjtre'
        call fafread(ndisk,bl(iaddr),8,isize,1,ip(irep),info)
          if (info.lt.8*isize) then
              ERROR
              call flush(6)
          endif
            else
          call general_read(ndisk,ip(irep),bl(iaddr),nval*nval*idimen*8)
            endif
            exit
          endif
        enddo g
        if (irecords(ip(irep)).eq.1) call abort()
      endif
c
      enddo
c
cDEB  do ia=1,idimen
cDEB    if (irecords(ia).ne.1) write(*,'(I3)',ADVANCE='NO') ia
cDEB  enddo
cDEB  write(*,*)
c
      shall_I_exit=.true.
      m: do irep=1,ip_len
           if (irecords(ip(irep)).eq.1) then
             shall_I_exit=.false.
             exit m
           endif
         enddo m
      if (shall_I_exit) return
      call abort()
      STOP 'Cannot leave this procedure (int) with such condition .tr'
      end
c===============================================================================
      subroutine locate_3int(a1,    a2,       b1,      b2,   c1,
     *                       c2,    i3ext_int,irecords,reset,n,
     *                       idimen,nval,     ndisk,   af)
      use memory
      implicit none
      integer a1,b1,c1,a2,b2,c2,i3ext_int,irecords(*),n,idimen,nval
      integer ndisk
      logical reset,af
c
      integer iaddr,iaddr1,irep
      integer addresses(n),iad_max,ia,ip(3*idimen)
      logical shall_I_exit
      integer i,ip_len,ii,info
c
      i=0
      do ia=a1,a2
        i=i+1
        ip(i)=ia
      enddo
      do ia=b1,b2
        i=i+1
        ip(i)=ia
      enddo
      do ia=c1,c2
        i=i+1
        ip(i)=ia
      enddo
      ip_len=i
c
      if (reset) then
        reset=.false.
        do ia=1,n
          iaddr        =i3ext_int+(ia-1)*nval*nval*nval
          addresses(ia)=iaddr
          if (iaddr.eq.1) 
     *  STOP 'Cannot be! I use assumption that the address = 1 is null!'
        enddo
        iad_max=n
        do ia=1,idimen
          irecords(ia)=1
        enddo
        if (a2.gt.n.or.b2.gt.n.or.c2.gt.n) then ! load a,b,c first
          do irep=1,ip_len
          if (irecords(ip(irep)).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ip(irep))=iaddr
          if (af) then
        call fafread(ndisk,bl(iaddr),8,nval*nval*nval,1,ip(irep),info)
          if (info.lt.8*nval*nval*nval) then
              ERROR
              call flush(6)
          endif
          else
          call general_read(ndisk,ip(irep),bl(iaddr),nval*nval*nval*8)
          endif
          enddo
        endif
        do ia=1,n
          if (iad_max.le.0) exit
          if (irecords(ia).ne.1) cycle
          iaddr=addresses(iad_max)
          iad_max=iad_max-1
          irecords(ia)=iaddr
          if (af) then
            call fafread(ndisk,bl(iaddr),8,nval*nval*nval,1,ia,info)
          if (info.lt.8*nval*nval*nval) then
              ERROR
              call flush(6)
          endif
          else
            call general_read(ndisk,ia,bl(iaddr),nval*nval*nval*8)
          endif
        enddo
      endif
c
      shall_I_exit=.true.
      z: do irep=1,ip_len
           if (irecords(ip(irep)).eq.1) then
             shall_I_exit=.false.
             exit z
           endif
         enddo z
      if (shall_I_exit) return
c
      do irep=1,ip_len
      if (irecords(ip(irep)).eq.1) then
      g:do ia=idimen,1,-1 ! search for largest record which can be replaced
          if (irecords(ia).ne.1) then
            if (ia.eq.ip(irep)) STOP 'Cannot be: ia.eq.i'
ccccc       if (ia.eq.ind2.or.ia.eq.ind3) cycle
            do ii=1,ip_len
              if (ip(ii).eq.ia) cycle g
            enddo
            iaddr=irecords(ia)
            irecords(ia)=1
            if (irecords(ip(irep)).ne.1) STOP 'irecords(i).ne.1'
            irecords(ip(irep))=iaddr
            if (af) then
          call fafread(ndisk,bl(iaddr),8,nval*nval*nval,1,ip(irep),info)
          if (info.lt.8*nval*nval*nval) then
              ERROR
              call flush(6)
          endif
            else
            call general_read(ndisk,ip(irep),bl(iaddr),nval*nval*nval*8)
            endif
            exit
          endif
        enddo g
      endif
c
      enddo
c
cDEB  do ia=1,idimen
cDEB    if (irecords(ia).ne.1) write(*,'(I3)',ADVANCE='NO') ia
cDEB  enddo
cDEB  write(*,*)
      shall_I_exit=.true.
      m: do irep=1,ip_len
           if (irecords(ip(irep)).eq.1) then
             shall_I_exit=.false.
             exit m
           endif
         enddo m
      if (shall_I_exit) return
c
      call abort()
      STOP 'Cannot leave this procedure (int) with such condition .tr'
      end
c===============================================================================
      subroutine locate_3ext(a1,     a2,        b1,       b2,    c1,
     *                       c2,     i3ext_int, irecords, reset, n,
     *                       idimen, nval,      ndisk,    ip,    af,
     *                       nrep,vtable,m_table,otable,vrevtable,
     *                       ireal_size)
      use memory
      implicit none
      integer a1,b1,c1,a2,b2,c2,i3ext_int,irecords(*),n,idimen,nval
      integer ndisk,ip(*),ireal_size
      integer nrep,vtable(*),m_table(*),otable(*),vrevtable(*)
      logical reset,af
c
      integer, save :: pairs(6,3)
      integer igranpair(5,3),nvirpairs,last_add,i,j,k,irem,ia,ib
      integer ip_len,aa1,aa2,bb1,bb2,ii,ibase,irec,iaddr,isize,info
      logical i_present,jremove,will_I_return
c
      call fill_ip1(a1,a2,b1,b2,c1,c2,idimen,ip,ip_len)
      nvirpairs=idimen*(idimen+1)/2
      last_add=i3ext_int+(n-1)*ireal_size*2
      if (reset) then
        reset=.false.
        pairs(1,1)=i3ext_int
        pairs(2,1)=0
        pairs(1,2)=i3ext_int+n/3*ireal_size*2
        pairs(2,2)=0
        pairs(1,3)=i3ext_int+n/3*2*ireal_size*2
        pairs(2,3)=0
c
        pairs(3,1)=0
        pairs(4,1)=0
        pairs(5,1)=0
        pairs(6,1)=0
c
        pairs(3,2)=0
        pairs(4,2)=0
        pairs(5,2)=0
        pairs(6,2)=0
c
        pairs(3,3)=0
        pairs(4,3)=0
        pairs(5,3)=0
        pairs(6,3)=0
c
        do ia=1,nvirpairs
          irecords(ia)=1
        enddo
      endif
c
      igranpair(1,1)=a1*(a1-1)/2+b1
        igranpair(2,1)=a1
        igranpair(3,1)=a2
        igranpair(4,1)=b1
        igranpair(5,1)=b2
      igranpair(1,2)=a1*(a1-1)/2+c1
        igranpair(2,2)=a1
        igranpair(3,2)=a2
        igranpair(4,2)=c1
        igranpair(5,2)=c2
      igranpair(1,3)=b1*(b1-1)/2+c1
        igranpair(2,3)=b1
        igranpair(3,3)=b2
        igranpair(4,3)=c1
        igranpair(5,3)=c2
c
      do i=1,3
c       Do we have the igranpair(i) already loaded?
        i_present=.false.
        do j=1,3
          if (igranpair(1,i).eq.pairs(2,j)) i_present=.true.
        enddo
        if (i_present) cycle ! if yes, do not do anything
c
c       Find pair which can be removed:
        irem=0
        do j=1,3
          jremove=.true.
          do k=1,3
            if (pairs(2,j).eq.igranpair(1,k)) jremove=.false.
          enddo
          if (jremove) then
            irem=j
            exit
          endif
        enddo
        if (irem.eq.0) STOP '4523kjsdfhj396ksjv7tjr'
c put "1" in old record locations:
        aa1=pairs(3,irem)
        aa2=pairs(4,irem)
        bb1=pairs(5,irem)
        bb2=pairs(6,irem)
        do ia=aa1,aa2
          if (ia.eq.0) exit
          do ib=bb1,bb2
            if (ib.eq.0) exit
            if (ib.gt.ia) exit
            irec=(ia-1)*ia/2+ib
            irecords(irec)=1
          enddo
        enddo
c
c  Now we replace pair irem by i
        pairs(2,irem)=igranpair(1,i)
        pairs(3,irem)=igranpair(2,i)
        pairs(4,irem)=igranpair(3,i)
        pairs(5,irem)=igranpair(4,i)
        pairs(6,irem)=igranpair(5,i)
c     
        aa1=igranpair(2,i)
        aa2=igranpair(3,i)
        bb1=igranpair(4,i)
        bb2=igranpair(5,i)
        ibase=pairs(1,irem)
        ii=-1
        do ia=aa1,aa2
          do ib=bb1,bb2
            if (ib.gt.ia) exit
            ii=ii+1
            irec=(ia-1)*ia/2+ib
            iaddr=ibase+ii*ireal_size*2
            if (iaddr.gt.last_add) STOP '4o2lklxdgr4368efeus87s'
            if (iaddr.lt.i3ext_int) STOP '4o2lklccxdgr4368efeus87s'
            irecords(irec)=iaddr
c
            if (af) then
                call ext3_sizer(nrep,idimen,nval,vtable,m_table,
     *                          otable,vrevtable,ia,ib,isize)
                call fafread(ndisk,bl(iaddr),8,isize*2,1,irec,info)
              if (info.lt.8*isize*2) then
                ERROR
                call flush(6)
              endif
              if (isize.gt.ireal_size) STOP '5020ektjtre'
            else
              call general_read(ndisk,irec,bl(iaddr),nval*idimen*16)
              if (nval*idimen.gt.ireal_size) STOP '5020ektjtre'
            endif
c
          enddo
        enddo
c
      enddo
c
      will_I_return=.true.
      do i=1,ip_len
        if (irecords(ip(i)).eq.1) then
          will_I_return=.false.
        exit
        endif
      enddo
      if (will_I_return) return
      print *, 'Cannot leave this procedure (int) with condition .tr'
      call flush(6)
      call abort()
      STOP 'Cannot leave this procedure (int) with such condition .tr'
      end
C===============================================================================
      subroutine locate_K   (a1,     a2,        b1,       b2,    c1,
     *                       c2,     i3ext_int, irecords, reset, n,
     *                       idimen, nval,      ndisk,    ip,    af,
     *                       isize, npass,    iap,         ibp,   icp)
      use memory
      implicit none
      integer a1,b1,c1,a2,b2,c2,i3ext_int,irecords(*),n,idimen,nval
      integer ndisk,ip(*),isize,npass,iap,ibp,icp
      logical reset,af
c
      integer, save :: pairs(6,3)
      integer igranpair(6,3),nvirpairs,last_add,i,j,k,irem,ia,ib
      integer ip_len,aa1,aa2,bb1,bb2,ii,ibase,irec,iaddr,info,ipass
      logical i_present,jremove,will_I_return
      real*8 xtmp(nval*nval,isize)
c
c
      call fill_ip1(a1,a2,b1,b2,c1,c2,idimen,ip,ip_len)
      nvirpairs=idimen*(idimen+1)/2
      last_add=i3ext_int+(n-1)*nval*nval
      if (reset) then
        reset=.false.
        pairs(1,1)=i3ext_int
        pairs(2,1)=0
        pairs(1,2)=i3ext_int+n/3*nval*nval
        pairs(2,2)=0
        pairs(1,3)=i3ext_int+n/3*2*nval*nval
        pairs(2,3)=0
c
        pairs(3,1)=0
        pairs(4,1)=0
        pairs(5,1)=0
        pairs(6,1)=0
c
        pairs(3,2)=0
        pairs(4,2)=0
        pairs(5,2)=0
        pairs(6,2)=0
c
        pairs(3,3)=0
        pairs(4,3)=0
        pairs(5,3)=0
        pairs(6,3)=0
c
        do ia=1,nvirpairs
          irecords(ia)=1
        enddo
      endif
c
      igranpair(1,1)=a1*(a1-1)/2+b1
        igranpair(2,1)=a1
        igranpair(3,1)=a2
        igranpair(4,1)=b1
        igranpair(5,1)=b2
        igranpair(6,1)=ibp
      igranpair(1,2)=a1*(a1-1)/2+c1
        igranpair(2,2)=a1
        igranpair(3,2)=a2
        igranpair(4,2)=c1
        igranpair(5,2)=c2
        igranpair(6,2)=icp
      igranpair(1,3)=b1*(b1-1)/2+c1
        igranpair(2,3)=b1
        igranpair(3,3)=b2
        igranpair(4,3)=c1
        igranpair(5,3)=c2
        igranpair(6,3)=icp
c
      do i=1,3
c       Do we have the igranpair(i) already loaded?
        i_present=.false.
        do j=1,3
          if (igranpair(1,i).eq.pairs(2,j)) i_present=.true.
        enddo
        if (i_present) cycle ! if yes, do not do anything
c
c       Find pair which can be removed:
        irem=0
        do j=1,3
          jremove=.true.
          do k=1,3
            if (pairs(2,j).eq.igranpair(1,k)) jremove=.false.
          enddo
          if (jremove) then
            irem=j
            exit
          endif
        enddo
        if (irem.eq.0) STOP '4523kjsdfhj396ksjv7tjr'
c put "1" in old record locations:
        aa1=pairs(3,irem)
        aa2=pairs(4,irem)
        bb1=pairs(5,irem)
        bb2=pairs(6,irem)
        do ia=aa1,aa2
          if (ia.eq.0) exit
          do ib=bb1,bb2
            if (ib.eq.0) exit
            if (ib.gt.ia) exit
            irec=(ia-1)*ia/2+ib
            irecords(irec)=1
          enddo
        enddo
c
c  Now we replace pair irem by i
        pairs(2,irem)=igranpair(1,i)
        pairs(3,irem)=igranpair(2,i)
        pairs(4,irem)=igranpair(3,i)
        pairs(5,irem)=igranpair(4,i)
        pairs(6,irem)=igranpair(5,i)
c     
        aa1=igranpair(2,i)
        aa2=igranpair(3,i)
        bb1=igranpair(4,i)
        bb2=igranpair(5,i)
        ipass=igranpair(6,i)
        ibase=pairs(1,irem)
        ii=-1
        do ia=aa1,aa2
          call irec_loc(ia,ipass,isize,npass,idimen,irec)
c
            if (af) then
                call fafread(ndisk,xtmp,8,nval*nval*isize,1,irec,info)
              if (info.lt.8*nval*nval*isize) then
                ERROR
                call flush(6)
              endif
            else
              call general_read(ndisk,irec,xtmp,nval*nval*isize*8)
            endif
c
          do ib=bb1,bb2
            if (ib.gt.ia) exit
            ii=ii+1
            iaddr=ibase+ii*nval*nval
            irec=(ia-1)*ia/2+ib
            if (iaddr.gt.last_add) STOP '4o2lklxdgr4368efeus87s'
            if (iaddr.lt.i3ext_int) STOP '4o2lklccxdgr4368efeus87s'
            irecords(irec)=iaddr
            call tfer(xtmp(1,ib-bb1+1),bl(iaddr),nval*nval)
          enddo
        enddo
c
      enddo
c
      will_I_return=.true.
      do i=1,ip_len
        if (irecords(ip(i)).eq.1) then
          will_I_return=.false.
        exit
        endif
      enddo
      if (will_I_return) return
      print *, 'Cannot leave this procedure (int) with condition .tr'
      call flush(6)
      call abort()
      STOP 'Cannot leave this procedure (int) with such condition .tr'
      end
C===============================================================================
      subroutine irec_loc(iaa,ip,isize,npass,idimen,irec)
c For given iaa virtual orbital and ip pass (passes of isize size and
c number of them = npass) it returns the record number where the bunch
c of K(iaa,b) matrices are stored, for all b from the range defined by
c ip, but b<=iaa. This is a
c little bit tricky, as all the K(a,b) are defined for a>=b only, so 
c passes must take this triangular form into account
      implicit none
      integer iaa,ip,isize,npass,irec,idimen
c
      integer ia,ipass,istart
      logical ia_log
c
      irec=0
      ia_log=.false.
      x:do ia=1,idimen
        if (ia.eq.iaa) ia_log=.true.
        do ipass=1,npass
          istart=(ipass-1)*isize+1
          if (istart.gt.idimen) STOP '34k5hfwskwir4ut'
          if (istart.gt.ia) exit
          irec=irec+1
          if (ia_log.and.ipass.eq.ip) exit x
        enddo
      enddo x
      if (.not.ia_log) STOP '4945ygfssi4euryth'
      end
C===============================================================================
      subroutine fill_ip(a1,a2,b1,b2,c1,c2,idimen,ip)
      implicit none
      integer a1,a2,b1,b2,c1,c2,ip(*),idimen
c
      integer ia,ib,ic,i
c
      i=0
      do ia=a1,a2
        do ib=b1,b2
          i=i+1
          ip(i)=(ia-1)*idimen+ib
          i=i+1
          ip(i)=(ib-1)*idimen+ia
        enddo
      enddo
      do ia=a1,a2
        do ic=c1,c2
          i=i+1
          ip(i)=(ia-1)*idimen+ic
          i=i+1
          ip(i)=(ic-1)*idimen+ia
        enddo
      enddo
      do ib=b1,b2
        do ic=c1,c2
          i=i+1
          ip(i)=(ib-1)*idimen+ic
          i=i+1
          ip(i)=(ic-1)*idimen+ib
        enddo
      enddo
      end
c===============================================================================
      subroutine fill_ip1(a1,a2,b1,b2,c1,c2,idimen,ip,ip_len)
      implicit none
      integer a1,a2,b1,b2,c1,c2,ip(*),idimen,ip_len
c
      integer ia,ib,ic,i
c
      i=0
      do ia=a1,a2
        do ib=b1,b2
          if (ib.gt.ia) exit
          i=i+1
          ip(i)=(ia-1)*ia/2+ib
        enddo
      enddo
      do ia=a1,a2
        do ic=c1,c2
          if (ic.gt.ia) exit
          i=i+1
          ip(i)=(ia-1)*ia/2+ic
        enddo
      enddo
      do ib=b1,b2
        do ic=c1,c2
          if (ic.gt.ib) exit
          i=i+1
          ip(i)=(ib-1)*ib/2+ic
        enddo
      enddo
      ip_len=i
      end
c===============================================================================
      subroutine sort_Kext(nval,idimen,ndisk,irecadr,npairs,ndiskx,
     *                     ncf,lbinx,thresh,byt8,nmo,vorb,af,igid,
     *                     isize,npass)
      use memory
      use kinds
      implicit none
c
      integer nval,idimen,ndisk,irecadr(*),npairs,ndiskx,ncf,lbinx
      integer nmo,igid,isize,npass
      logical byt8,vorb,af
      real*8 thresh
c
c
      integer memory1,irecord,ifill_level,lbin,ibin,iindex,len,len1
      integer ndisk_sort,ixmat,no_bins,nvirpairs,info,nrec,istore,itable
      integer no_gr,igranules,igranulesize,igran2pair,ipair2gran,maxrec
      integer matsize,memavail,max_itable_dim
      character*256 scrfile,filename
      logical fileopen
c
      call matmark()
      call dynamic_mmark()
c
      nvirpairs=idimen*(idimen+1)/2
c
      call dynamic_getmem(nvirpairs,irecord)
      call izeroit(bl(irecord),nvirpairs*intsize)
      call dynamic_getmem(nvirpairs,ifill_level)
      call izeroit(bl(ifill_level),nvirpairs*intsize)
      call dynamic_show_free(memory1)
c
c I need memory for bins, and indexes (2 indexes)
c All the bins will occuppy the amount of memory:
c idimen*idimen*lbin*(8+2*2)
c
      memory1=memory1-300000-idimen*idimen ! possible fragmentation
      memory1=memory1-idimen*idimen ! one memory reservation in ext3_sort
      lbin=8*memory1/((nvirpairs+nvirpairs/8)*12) ! because of granules
      if (lbin.gt.nval*nval) lbin=nval*nval
      no_bins=nval*nval/lbin
      if (mod(nval*nval,lbin).ne.0) no_bins=no_bins+1
c
      maxrec=4194304 ! 2**22
      matsize=nval*nval
      call dynamic_getmem(nvirpairs,igran2pair)
      call dynamic_getmem(nvirpairs/i4size+1,ipair2gran)
      call dynamic_show_free(memavail)
      no_gr=idimen/isize
      if (mod(idimen,isize).ne.0) no_gr=no_gr+1
      memavail=memavail-nval*nval*isize
      memavail=memavail-idimen*no_gr/i4size-1
      call granules_calculate1(nvirpairs,idimen,lbin,bl(igran2pair),
     *                         bl(ipair2gran),igranules,igranulesize,
     *                         memavail,maxrec,matsize)
c
      call flush(6)
      call dynamic_mmark()
      call dynamic_getmem(igranules*igranulesize*lbin,ibin)
      call dynamic_getmem(igranules*igranulesize*lbin/2+1,iindex)
      call izeroit(bl(iindex),(nvirpairs*lbin/2+1)*intsize)
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len) !returns len of string without spaces,
                                    ! 80 is starting lenght
      filename=scrfile(1:len)//'.sorttmp'
      len1=len+8
      ndisk_sort=50
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk_sort,
     *                 lbin*12*igranulesize+8)
      else
        inquire(ndisk_sort,opened=fileopen)
        if (fileopen) STOP 'Unit 50 opened'
        OPEN(UNIT=ndisk_sort,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=lbin*12*igranulesize+8/intsize)
      endif
c
      nrec=0
      call K_sort(bl(ibin),bl(iindex),nval,bl(ifill_level),idimen,
     *            bl(irecord),ndisk_sort,lbin,irecadr,npairs,
     *            ndiskx,ncf,lbinx,thresh,byt8,nmo,vorb,af,igid,nrec,
     *            igranules,igranulesize,bl(igran2pair),bl(ipair2gran))
      call leftoversgr(idimen, bl(ifill_level),lbin,    bl(iindex),nrec,
     *             bl(irecord),ndisk_sort,     bl(ibin),af,        igid,
     *             nvirpairs,2,igranules,igranulesize,bl(igran2pair),
     *             bl(ipair2gran))
ccs      call sleep(5)
c
      filename=scrfile(1:len)//'.Kext_ab'
      len1=len+8
      ndisk=54
      if (af) then
      call fafCreatem(filename(1:len1)//char(0),ndisk,isize*nval*nval*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 54 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=isize*nval*nval*8)
      endif
      call dynamic_retmark()
c
      call dynamic_mmark()
      call dynamic_getmem(igranulesize*lbin,ibin)
      call dynamic_getmem(igranulesize*lbin/2+1,iindex)
      call dynamic_getmem(nval*nval*igranulesize,ixmat)
      call dynamic_getmem(nval*nval*isize,istore)
      max_itable_dim=(1+npass)*npass*isize/2
      call dynamic_getmem(max_itable_dim,itable) ! change! BUG!
      call coll_Kext_bins(lbin,bl(ibin),bl(iindex),ndisk_sort,
     *                  ndisk,bl(irecord),idimen,bl(ixmat),nval,af,
     *                  isize,npass,bl(istore),bl(itable),igranules,
     *                  igranulesize,bl(igran2pair),bl(ipair2gran))
      if (af) then
        call fafClosem(ndisk_sort,0,info)
      else
        close(ndisk_sort,status='delete')
      endif
      call dynamic_retmark()
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine K_sort(bin,iindex,nval,ifill_level,idimen,irecord,
     *                  ndisk_bin,lbin,irecadrx,npairs,ndiskx,ncf,
     *                  lbinx,thresh,byt8,nmo,vorb,af,igid,nrec,
     *                  igranules,igranulesize,igran2pair,ipair2gran)
      use memory
      implicit none
      integer lbin,nval,idimen,irecadrx(*),npairs,ndiskx,ncf,lbinx
      integer nmo,igid,nrec,igranules,igranulesize
      integer*4 igran2pair(2,*),ipair2gran(*)
      logical byt8,vorb,af
      real*8 bin(lbin,*),thresh
      integer*2 iindex(2,lbin,*)
      integer ifill_level(*),irecord(*)
      integer ndisk_bin,istart,istop
c
      integer i,j,ia,ib,ix,ic_ind,itmp,ibin,nvirpairs,ii,iab,isiz,igran
c
      call dynamic_getmem(idimen*idimen,itmp)
c
      nvirpairs=idimen*(idimen+1)/2
c
      isiz=igranulesize*lbin
      do i=1,nval
        do j=1,i
          call ExtrOne(i,    j,     irecadrx,  npairs,    ndiskx,
     *                 ncf,   lbinx,thresh,     byt8,      'x',
     *                 'mo',  nmo,   vorb,       bl(itmp))
          ic_ind=-1
          do ib=1,idimen
            do ia=1,idimen
             if (ia.ge.ib) then
               ibin=ia*(ia-1)/2+ib
             else
               ibin=ib*(ib-1)/2+ia
             endif
c
             ix=ifill_level(ibin)
             if (ix.ge.lbin) then
             if (ix.gt.lbin) STOP 'Error ddfdfdfdf4433'
               igran=ipair2gran(ibin)
               istart=igran2pair(1,igran)
               istop =igran2pair(2,igran)
               do iab=istart,istop
                 do ii=ifill_level(iab)+1,lbin
                   iindex(1,ii,iab)=0
                 enddo
                 ifill_level(iab)=0
               enddo
               call amp_bin_wrt(nrec,igran,irecord(igran),ndisk_bin,
     *                      isiz,bin(1,istart),iindex(1,1,istart),af,
     *                      igid,2)
             endif
c
             ic_ind=ic_ind+1
             ix=ifill_level(ibin)+1
             ifill_level(ibin)=ix
c
             bin(ix,ibin)=bl(itmp+ic_ind)
             if (ia.ge.ib) then
               iindex(1,ix,ibin)=i
               iindex(2,ix,ibin)=j
             else
               iindex(1,ix,ibin)=j
               iindex(2,ix,ibin)=i
             endif
            enddo
          enddo
        enddo
      enddo
      call dynamic_retmem(1)
      end
c===============================================================================
      subroutine coll_Kext_bins(lbin,bin,iindex,ndisk_s,ndisk_t,irecord,
     *                        idimen,xmat,nval,af,isize,npass,store,
     *                        itable,igranules,igranulesize,igran2pair,
     *                        ipair2gran)
      use memory
      implicit none
      integer lbin,ndisk_s,ndisk_t,idimen,nval,isize,npass,igranules
      integer igranulesize
      integer*4 igran2pair(2,*),ipair2gran(*)
      real*8 bin(lbin,igranulesize),xmat(nval,nval,igranulesize)
      real*8 store(nval,nval,isize)
      integer*2 iindex(2,lbin,igranulesize)
      integer irecord(*)
      integer*4 itable(*)
      logical af
c
      integer ia,nrec,nrec1,ii,i,j,ib,iab,id,igid,info,ind(4)
      integer ib_ind,istart,istop,ipass,isiz,igr,inde,irecij,iadd
      logical equal,reset
c
      reset=.true.
      isiz=igranulesize*lbin
      if (af) then
      do igr=1,igranules
        call zeroit(xmat,nval*nval*igranulesize)
        istart=igran2pair(1,igr)
        istop =igran2pair(2,igr)
        irecij=irecord(igr)
        do
          if (irecij.eq.0) exit
          call bufpointer(iadd,bl(1),isiz*12+8)
          call fafread(ndisk_s,bl(iadd),1,isiz*12+8,1,irecij,info)
          if (info.lt.12*isiz+8) then
            ERROR
            call flush(6)
          endif
          call unpacker(3,bin,isiz*8,iindex,isiz*4,irecij,8)
          if (irecij.lt.0) STOP 'Error fdjh5e6533k2thwwh'
        iab=0
        do ia=1,idimen
          do ib=1,ia
          iab=iab+1
          if (iab.lt.istart) cycle
          if (iab.gt.istop ) exit
          inde=iab-istart+1
          equal=.false.
          if (ia.eq.ib) equal=.true.
          do ii=1,lbin
            i = iindex(1,ii,inde)
            j = iindex(2,ii,inde)
            if (i.le.0) then
              if (i.lt.0) STOP 'Error 33k3nfn5rh'
              exit ! exit loop do ii=1,lbin
            endif
            xmat(i,j,inde)=bin(ii,inde)
            if (equal) xmat(j,i,inde)=xmat(i,j,inde)
          enddo
        enddo
        enddo
        enddo
        iab=0
        do ia=1,idimen
          do ib=1,ia
          iab=iab+1
          if (iab.lt.istart) cycle
          if (iab.gt.istop ) exit
          inde=iab-istart+1
          call acc_and_write_K(idimen,npass,isize,iab,xmat(1,1,inde),
     *                         store,ndisk_t,af,nval*nval,reset,
     *                         itable)
          enddo
        enddo
      enddo
      else
c
      do igr=1,igranules
        call zeroit(xmat,nval*nval*igranulesize)
        istart=igran2pair(1,igr)
        istop =igran2pair(2,igr)
        nrec=irecord(igr)
        do
          read(ndisk_s,rec=nrec) nrec1,bin,iindex
          iab=0
          do ia=1,idimen
            do ib=1,ia
              iab=iab+1
              if (iab.lt.istart) cycle
              if (iab.gt.istop ) exit
              inde=iab-istart+1
              equal=.false.
              if (ia.eq.ib) equal=.true.
              do ii=1,lbin
                i = iindex(1,ii,inde)
                j = iindex(2,ii,inde)
                if (i.le.0) then
                  if (i.lt.0) STOP 'Error 33k3nfn5rh'
                  exit ! exit loop ii=1,lbin
                endif
                xmat(i,j,inde)=bin(ii,inde)
                if (equal) xmat(j,i,inde)=xmat(i,j,inde)
              enddo
            enddo
          enddo
          nrec=nrec1
          if (nrec.eq.0) exit
        enddo
        iab=0
        do ia=1,idimen
          do ib=1,ia
            iab=iab+1
            if (iab.lt.istart) cycle
            if (iab.gt.istop ) exit
            inde=iab-istart+1
            call acc_and_write_K(idimen,npass,isize,iab,xmat(1,1,inde),
     *                           store,ndisk_t,af,nval*nval,reset,
     *                           itable)
          enddo
        enddo
      enddo
      endif
      end
c===============================================================================
c Cannot be parallelized in the current form. Assumes continuous 'iab' flow.
      subroutine acc_and_write_K(idimen,npass,isize,iab,xmat,
     *                           store,ndisk,af,mat_size,reset,
     *                           itable)
      use ccounters
      implicit none
      integer idimen,npass,isize,iab,ndisk,mat_size
      logical af,reset
      real*8 xmat(mat_size),store(mat_size,isize)
      integer*4 itable(2,*)
c
      integer max_dim
      integer curr_size,ia,ipass
      integer istart,istop,irec,iacc,info,i
c
      max_dim=(1+npass)*npass*isize/2
      if (reset) then
        curr_rec=1
        curr_pos=1
        reset=.false.
        irec=0
        iacc=0
        do ia=1,idimen
          do ipass=1,npass
            istart=(ipass-1)*isize+1
            istop =(ipass  )*isize
            if (istart.gt.idimen) STOP '34k5hfwskwir4ut'
            if (istart.gt.ia) exit
            if (istop .gt.ia) istop=ia
            irec=irec+1
            if (irec.gt.max_dim) STOP 'nfhk6wfh54eeddd'
            itable(1,irec)=iacc+istart
            itable(2,irec)=iacc+istop
            max_rec_K=irec
          enddo
          iacc=iacc+ia
        enddo
      endif
      if (iab.lt.itable(1,curr_rec)) then
        do i=1,max_rec_K
          print *,itable(1,i),itable(2,i)
        enddo
          print *,iab,curr_rec,curr_pos,itable(1,curr_rec),
     *            itable(2,curr_rec)
        call flush(6)
        STOP '5ffffffdfvfff'
      endif
      if (iab.gt.itable(2,curr_rec)) STOP 'nbjgtiorld'
      call tfer(xmat,store(1,curr_pos),mat_size)
      curr_pos=curr_pos+1
      curr_size=itable(2,curr_rec)-itable(1,curr_rec)+1
      if (curr_pos.gt.curr_size) then
        if (af) then
          call fafwrite(ndisk,store,8,mat_size*isize,1,curr_rec,info)
        else
          write(ndisk,rec=curr_rec) store
        endif
        curr_pos=1
        curr_rec=curr_rec+1
      endif
      end
c===============================================================================
      subroutine check_kext(nval,idimen,ndisk)
      implicit none
      integer nval,idimen,ndisk
      real*8 xmat(nval,nval)
      real*8 x(nval,idimen,nval,idimen)
      integer iab,i,j,ia,ib
      iab=0
      do ia=1,idimen
        do ib=1,ia
          iab=iab+1
          call general_read(ndisk,iab,xmat,nval*nval*8)
          do i=1,nval
            do j=1,nval
              x(i,ia,j,ib)=xmat(i,j)
              x(j,ib,i,ia)=xmat(i,j)
            enddo
          enddo
        enddo
      enddo
      do i=1,nval
        do ia=1,idimen
          do j=1,nval
            do ib=1,nval
              write(99,'(4I4,F25.10)') i,ia,j,ib,x(i,ia,j,ib)
            enddo
          enddo
        enddo
      enddo
      end
c===============================================================================
c===============================================================================
      subroutine check_kext1(nval,idimen,ndisk,irecadrx,npairs,ndiskx,
     *                       ncf,lbinx,thresh,byt8,nmo,vorb,x)
      implicit none
      integer nval,idimen,ndisk,irecadrx(*),npairs,ndiskx,ncf,lbinx,nmo
      real*8 xmat(idimen,idimen),thresh
      logical byt8,vorb
      real*8 x(nval,idimen,nval,idimen)
      integer iab,i,j,ia,ib
      iab=0
      do i=1,nval
        do j=1,i
          iab=iab+1
          call ExtrOne(i,    j,     irecadrx,  npairs,    ndiskx,
     *                 ncf,   lbinx,thresh,     byt8,      'x',
     *                 'mo',  nmo,   vorb,       xmat)
          do ia=1,idimen
            do ib=1,idimen
              x(i,ia,j,ib)=xmat(ia,ib)
              x(j,ib,i,ia)=xmat(ia,ib)
            enddo
          enddo
        enddo
      enddo
      do i=1,nval
        do ia=1,idimen
          do j=1,nval
            do ib=1,nval
              write(99,'(4I4,F25.10)') i,ia,j,ib,x(i,ia,j,ib)
            enddo
          enddo
        enddo
      enddo
      call flush(99)
      end
c===============================================================================
      subroutine return_addresses(ia,ib,ic,irec_K,iKab,iKac,iKbc)
      implicit none
      integer ia,ib,ic,irec_K(*),iKab,iKac,iKbc
      iKab=irec_K(  ia*(ia-1)/2+ib   )
      iKac=irec_K(  ia*(ia-1)/2+ic   )
      iKbc=irec_K(  ib*(ib-1)/2+ic   )
      if (iKab.eq.1.or.iKac.eq.1.or.iKbc.eq.1) STOP 'Error 54023cjasl4'
      end
c===============================================================================
c Integrals - 3ext J(d,i)(a,b) type -> K(b,d)(a,i) type, direct sort
c===============================================================================
      subroutine simple_3ext(nval,idimen,ndisk,irecadrtt,npairs,ndisktt,
     *                       ncf,lbintt,thresh,byt8,nmo,vorb,af,igid,
     *                       nrep,vtable,m_table,vrevtable,otable)
      use memory
      use kinds
      implicit none
      integer nval,idimen,ndisk,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nrep,vtable(*),m_table(*),vrevtable(*),otable(*)
      integer nmo,igid
      real*8 thresh
      logical byt8,vorb,af
c
      character*256 scrfile,filename
      integer len,len1,ibigmat,ixmat,nslv,istart,istop,ipasssize,npass
      integer isize,ia
      logical fileopen
      real*8 c_ratio
c
c
      nslv=3
      c_ratio=0d0
      call matmark()
      call dynamic_mmark()
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len)
      filename=scrfile(1:len)//'.Kbc'
      len1=len+4
      ndisk=52
      if (af) then
        call fafCreatem(filename(1:len1)//char(0),ndisk,nval*idimen*8)
      else
        inquire(ndisk,opened=fileopen)
        if (fileopen) STOP 'Unit 52 opened'
        OPEN(UNIT=ndisk,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*8)
      endif
c
      do isize=idimen,1,-1
        npass=idimen/isize
        if (mod(idimen,isize).ne.0) npass=npass+1
        if (npass.ge.nslv) then
          ipasssize=isize
          exit
        endif
      enddo
      call optimal_isize(ipasssize,idimen)
      npass=idimen/ipasssize
      if (mod(idimen,ipasssize).ne.0) npass=npass+1
c
      call dynamic_getmem(idimen*idimen*nval,ibigmat)
      call dynamic_getmem(idimen*idimen,ixmat)
      do ia=1,npass
        istart=(ia-1)*ipasssize+1
        istop= (ia  )*ipasssize
        if (istart.gt.idimen) STOP 'Error 392jfallaer'
        if (istop .gt.idimen) istop=idimen
        call direct_ext3_sort(idimen,nval,bl(ixmat),ndisk,irecadrtt,
     *                        npairs,ndisktt,ncf,lbintt,thresh,
     *                        byt8,nmo,vorb,bl(ibigmat),istart,
     *                        istop,af,nrep,vtable,m_table,vrevtable,
     *                        otable,c_ratio)
      enddo
      call matremark()
      call dynamic_retmark()
      end
c===============================================================================
      subroutine direct_ext3_sort(idimen,nval,xmat,ndisk,irecadrtt,
     *                           npairs,ndisktt,ncf,lbintt,thresh,
     *                           byt8,nmo,vorb,bigmat,istart,
     *                           istop,af,nrep,vtable,m_table,vrevtable,
     *                           otable,c_ratio)
      use memory
      implicit none
c
      integer idimen,nval,ndisk,irecadrtt(*),npairs,ndisktt,ncf,lbintt
      integer nmo,istart,istop
      integer nrep,vtable(*),m_table(*),vrevtable(*),otable(*)
      logical byt8,vorb,af
      real*8 thresh,xmat(idimen,idimen),bigmat(idimen,nval,idimen)
      real*8 c_ratio
c
      integer ia,ib,i,irec,id,info,iscr,isize
      real*8 local_ratio
c
c===============================================================================
c Integrals - 3ext J(d,i)(a,b) type -> K(b,d)(a,i) type, direct sort
c===============================================================================
      call dynamic_mmark()
      call dynamic_getmem(idimen*nval,iscr)
      do id=istart,istop
        do i=1,nval
          call ExtrOne(id,    i,     irecadrtt,  npairs,    ndisktt,
     *                 ncf,   lbintt,thresh,     byt8,      'tt',
     *                 'mo',  nmo,   vorb,       xmat)
          do ib=1,idimen
            do ia=1,idimen
              bigmat(ia,i,ib)=xmat(ia,ib)
            enddo
          enddo
        enddo
        do ib=1,idimen
          irec=(ib-1)*idimen+id
          if (nrep.gt.1) then
            call symm_compress(bigmat(1,1,ib),bl(iscr),nrep,idimen,nval,
     *                      vtable,m_table,otable,vrevtable,ib,id,isize)
            local_ratio=dble(isize)/dble(idimen*nval)
            if (local_ratio.gt.c_ratio) c_ratio=local_ratio
            if (af) then
              call fafwrite(ndisk,bl(iscr),8,isize,1,irec,info)
            else
              call general_write(ndisk,irec,bl(iscr),idimen*nval*8)
            endif
          else
            c_ratio=1d0
            if (af) then
           call fafwrite(ndisk,bigmat(1,1,ib),8,idimen*nval,1,irec,info)
            else
             call general_write(ndisk,irec,bigmat(1,1,ib),idimen*nval*8)
            endif
          endif
        enddo
      enddo
      call dynamic_retmark()
      end
c===============================================================================
      subroutine make_3ext_pairs(idimen,nval,af,ndisk_s,ndisk_t)
      use memory
      implicit none
      integer idimen,nval,ndisk_s,ndisk_t
      logical af
c
      character*256 scrfile,filename
      integer len,len1,ibase,iab,ia,ib,irec,iaddr,info
      logical fileopen
c
c
      call getchval('scrf',scrfile)
      call rmblan(scrfile,80,len)
      filename=scrfile(1:len)//'.Kbc1'
      len1=len+5
      ndisk_t=55
      if (af) then
       call fafCreatem(filename(1:len1)//char(0),ndisk_t,nval*idimen*16)
      else
        inquire(ndisk_t,opened=fileopen)
        if (fileopen) STOP 'Unit 55 opened'
        OPEN(UNIT=ndisk_t,FILE=filename(1:len1),FORM='UNFORMATTED',
     *       ACCESS='DIRECT',RECL=nval*idimen*16)
      endif
      call dynamic_getmem(2*idimen*nval,ibase)
      iab=0
      do ia=1,idimen
        do ib=1,ia
          iab=iab+1
          irec=(ia-1)*idimen+ib
          iaddr=ibase
          if (af) then
          call fafread(ndisk_s,bl(iaddr),8,idimen*nval,1,irec,info)
          if (info.lt.8*idimen*nval) then
              ERROR
              call flush(6)
          endif
          else
          call general_read(ndisk_s,irec,bl(iaddr),idimen*nval*8)
          endif
          irec=(ib-1)*idimen+ia
          iaddr=ibase+idimen*nval
          if (af) then
          call fafread(ndisk_s,bl(iaddr),8,idimen*nval,1,irec,info)
          if (info.lt.8*idimen*nval) then
              ERROR
              call flush(6)
          endif
          else
          call general_read(ndisk_s,irec,bl(iaddr),idimen*nval*8)
          endif
          if (af) then
          call fafwrite(ndisk_t,bl(ibase),8,idimen*nval*2,1,iab,info)
          else
          call general_write(ndisk_t,iab,bl(ibase),idimen*nval*16)
          endif
        enddo
      enddo
      if (af) then
        call fafClosem(ndisk_s,0,info)
      else
        close(ndisk_s,status='delete')
      endif
      call dynamic_retmem(1)
      end
c============================================================================
      subroutine symm_compress(source,dest,nrep,idimen,nval,vtable,
     *                         m_table,otable,vrevtable,ia,ib,isize)
      implicit none
      integer nrep,idimen,nval,vtable(idimen+1,*),m_table(8,8)
      integer otable(nval+1,*),vrevtable(*),ia,ib,isize
      real*8 source(idimen,*),dest(*)
c
      integer icount,id_rep,id_size,irep,ik_size,iaddr,ids,iks,ik,id
c
      if (nrep.le.1) then 
        isize=nval*idimen
        return
      endif
      icount=0
      do id_rep=1,nrep  ! loop over virtual representations
        id_size=vtable(1,id_rep)
        irep=m_table(id_rep,vrevtable(ia))
        irep=m_table(irep,vrevtable(ib))
        ik_size =otable(1,irep)
        if (id_size.eq.0.or.ik_size.eq.0) cycle
        iaddr=icount
        do ik=1,ik_size
          iks=otable(ik+1,irep)
          do id=1,id_size
            ids=vtable(id+1,id_rep)
            dest(iaddr+id+(ik-1)*id_size)=source(ids,iks)
            icount=icount+1
          enddo
        enddo
      enddo
      isize=icount
      end
c===============================================================================
      subroutine sym_comp_a(source,dest,nval,idimen,nrep,ia,m_table,
     *                      ijtable,vtable,vrevtable,isize)
      implicit none
      integer nval,idimen,nrep,ia,m_table(8,8),ijtable(nval*nval+1,*)
      integer vtable(idimen+1,*),vrevtable(*),isize
      real*8 source(nval*nval,*),dest(*)
c
      integer icount,id_rep,id_size,ijrep,ijsize,iaddr,id,ij,ids,ijs
c
      if (nrep.le.1) then 
        isize=nval*nval*idimen
        return
      endif
      icount=0
      do id_rep=1,nrep  ! loop over virtual representations
        id_size=vtable(1,id_rep)
        ijrep=m_table(id_rep,vrevtable(ia))
        ijsize=ijtable(1,ijrep)
        iaddr=icount
        if (id_size.eq.0.or.ijsize.eq.0) cycle
        do id=1,id_size
          ids=vtable(id+1,id_rep)
          do ij=1,ijsize
            icount=icount+1
            ijs=ijtable(ij+1,ijrep)
c           call ijresolver(nval,ijs,i,j)
            dest(iaddr+ij+(id-1)*ijsize)=source(ijs,ids)
          enddo
        enddo
      enddo
      isize=icount
      end
c===============================================================================
      subroutine ext3_sizer(nrep,idimen,nval,vtable,m_table,
     *                      otable,vrevtable,ia,ib,isize)
      implicit none
      integer nrep,idimen,nval,vtable(idimen+1,*),m_table(8,8)
      integer otable(nval+1,*),vrevtable(*),ia,ib,isize
c
      integer id_rep,irep,ik_size,id_size
c
      if (nrep.le.1) then 
        isize=idimen*nval
        return
      endif
      isize=0
      do id_rep=1,nrep  ! loop over virtual representations
        id_size = vtable(1,id_rep)
        irep    = m_table(id_rep,vrevtable(ia))
        irep    = m_table(irep,vrevtable(ib))
        ik_size = otable(1,irep)
        if (id_size.eq.0.or.ik_size.eq.0) cycle
        isize=isize+id_size*ik_size
      enddo
      end
c===============================================================================
      subroutine ampl_sizer(nval,idimen,nrep,ia,m_table,
     *                      ijtable,vtable,vrevtable,isize)
      implicit none
      integer nval,idimen,nrep,ia,m_table(8,8),ijtable(nval*nval+1,*)
      integer vtable(idimen+1,*),vrevtable(*),isize
c
      integer icount,id_rep,id_size,ijrep,ijsize
c
      if (nrep.le.1) then 
        isize=nval*nval*idimen
        return
      endif
      icount=0
      do id_rep=1,nrep  ! loop over virtual representations
        id_size=vtable(1,id_rep)
        ijrep=m_table(id_rep,vrevtable(ia))
        ijsize=ijtable(1,ijrep)
        if (id_size.eq.0.or.ijsize.eq.0) cycle
        icount=icount+id_size*ijsize
      enddo
      isize=icount
      end
c===============================================================================
      subroutine granules_calculate1(npairs,idimen,lbin,igran2pair,
     *                              ipair2gran,igranules,igranulesize,
     *                              memavail,maxrec,matsize)
      implicit none
      integer npairs,idimen,lbin,igranules,igranulesize,memavail,maxrec
      integer matsize
      integer*4 igran2pair(2,*),ipair2gran(*)
c
      integer isiz,lrec,mem,istart,istop,i,j
c
c Determine the size of granule by bin record size:
c
      do isiz=npairs,0,-1
        if (isiz.eq.0) STOP 'vdssk6g49eghhv3kstttd'
        lrec =12*lbin*isiz
        mem=(matsize+lbin+lbin/2+1)*isiz
        if (mem.lt.memavail.and.lrec.lt.maxrec) then
          igranulesize=isiz
          exit
        endif
      enddo
c
      call optimal_isize(igranulesize,npairs)
      igranules=npairs/igranulesize
      if (mod(npairs,igranulesize).ne.0) igranules=igranules+1
      lrec =12*lbin*igranulesize  ! record length in bytes
c build maps: granule -> ipairstart,ipairstop
c             ipair   -> granule
      istart=1
      do i=1,igranules
        istop=istart+igranulesize-1
        if (istop.gt.npairs) istop=npairs
        if (istart.gt.npairs) call nerror(3,'MP2 module',
     $   'Error in granule maps determination ',istart,npairs)
        igran2pair(1,i)=istart
        igran2pair(2,i)=istop
        do j=istart,istop
          ipair2gran(j)=i
        enddo
        istart=istop+1
      enddo
c
      end
c===============================================================================
      subroutine calc_chunk_size(idimen,nval,amp_ratio,ext3_ratio,mem,
     *                           nslv,npass,isize)
      implicit none
      integer idimen,nval,mem,nslv,npass,isize
      real*8 amp_ratio,ext3_ratio
c
c
      integer i,ip_array,fragmentation,iampl_store,i3int_store
      integer i3ext_store,iK_store,ialloc_amount,nvirpairs
      integer ngranules,irem,ipass,itrysize
c
      nvirpairs=idimen*(idimen+1)/2
c
      do i=1,100000
        isize=i
c
        ip_array      =    6*isize*isize
        fragmentation =    300000
        iampl_store   =int(nval*nval*idimen*amp_ratio+1)*isize*3
        i3int_store   =    nval*nval*nval               *isize*3
        i3ext_store   =int(nval*idimen*ext3_ratio+1)    *isize*isize*6
        iK_store      =    nval*nval                    *isize*isize*3
c
        ialloc_amount= ip_array           +fragmentation
     *                +idimen             +iampl_store
     *                +idimen             +i3int_store
     *                +idimen*idimen      +i3ext_store
     *                +nvirpairs          +iK_store
c
        if (ialloc_amount.lt.mem) then 
          cycle
        else
          if (isize.gt.idimen) isize=idimen+1
          isize=isize-1
          exit
        endif
      enddo
c Now, try to have as much of chunks to keep every slave busy:
      npass=idimen/isize
      if (mod(idimen,isize).ne.0) npass=npass+1
      ngranules=(npass*(npass+1)*(npass+2))/6
c 
      if (ngranules/nslv.lt.20.and.nslv.gt.1) then
        do itrysize=isize,1,-1
          ipass=idimen/itrysize
          if (mod(idimen,itrysize).ne.0) ipass=ipass+1
          ngranules=(ipass*(ipass+1)*(ipass+2))/6
          irem=mod(ngranules,nslv)
c         if (ipass.gt.idimen/4) exit
          if (ngranules/nslv.gt.20) exit
          if ((ngranules.gt.1).and.irem.eq.0
     *              .or.(nslv-irem).lt.nslv/4) exit
        enddo
        isize=itrysize
        if (isize.eq.0) isize=1
        npass=idimen/isize
        if (mod(idimen,isize).ne.0) npass=npass+1
      endif
c try not to have a very small last chunk: homogenization of chunks,
c but keeping the same npass value
c     isize=7
      call optimal_isize(isize,idimen)
      npass=idimen/isize
      if (mod(idimen,isize).ne.0) npass=npass+1
      end
c===============================================================================
      subroutine array_to_scalar(array,scalar)
      implicit none
      integer*8 scalar
      integer array(8)
c
      integer len_mon(12)
      data len_mon /31,28,31,30,31,30,31,31,30,31,30,31/
      logical leap
      integer iyear,imonth,iday,ihour,imin,isec,len_year,len_year_leap,i
c
      scalar=0
      iyear  = array(1)
      imonth = array(2)
      iday   = array(3)
      ihour  = array(5)
      imin   = array(6)
      isec   = array(7)
      len_year=     365*24*3600
      len_year_leap=366*24*3600
      do i=1900,iyear-1
        if (leap(i)) then
          scalar=scalar+len_year_leap
        else
          scalar=scalar+len_year
        endif
      enddo
      do i=1,imonth-1
        scalar=scalar+len_mon(i)*24*3600
        if (i.eq.2.and.leap(iyear)) scalar=scalar+24*3600
      enddo
      do i=1,iday-1
        scalar=scalar+24*3600
      enddo
      do i=0,ihour-1
        scalar=scalar+3600
      enddo
      do i=0,imin-1
        scalar=scalar+60
      enddo
      do i=0,isec-1
        scalar=scalar+1
      enddo
      end
c===============================================================================
      subroutine scalar_to_array(scalar,array)
      implicit none
      integer*8 scalar
      integer array(8)
c
      integer len_year,len_year_leap,sub,i,iyear
      integer imonth,iday,ihour,imin,isec
      logical leap
      integer len_mon(12)
      data len_mon /31,28,31,30,31,30,31,31,30,31,30,31/
c
      len_year=     365*24*3600
      len_year_leap=366*24*3600
      iyear=1900
      do i=1900,5000
        if (leap(i)) then
          sub=len_year_leap
        else
          sub=len_year
        endif
        if (scalar-sub.lt.0) then
          iyear=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      imonth=1
      do i=1,12
        sub=len_mon(i)*24*3600
        if (i.eq.2.and.leap(iyear)) sub=sub+24*3600
        if (scalar-sub.lt.0) then
          imonth=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      iday=1
      do i=1,31
        sub=24*3600
        if (scalar-sub.lt.0) then
          iday=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      ihour=0
      do i=0,59
        sub=3600
        if (scalar-sub.lt.0) then
          ihour=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      imin=0
      do i=0,59
        sub=60
        if (scalar-sub.lt.0) then
          imin=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      isec=0
      do i=0,59
        sub=1
        if (scalar-sub.lt.0) then
          isec=i
          exit
        else
          scalar=scalar-sub
        endif
      enddo
      array(1) = iyear
      array(2) = imonth
      array(3) = iday
      array(5) = ihour
      array(6) = imin
      array(7) = isec
      array(8) = 0
      end
c===============================================================================
      function leap(year)
      implicit none
      integer year
      logical leap
      if ( mod( year, 400 ) .eq. 0 ) then
        leap = .true.
      else if( mod( year, 100 ) .eq. 0 ) then
        leap = .false.
      else if( mod( year, 4 ) .eq. 0 ) then
        leap = .true.
      else
        leap = .false.
      endif
      end
c===============================================================================
      subroutine sort_and_compress(results,itop,iprnt)
      implicit none
      real*8 results(2,*)
      integer itop,iprnt
c
      logical :: break
      integer i,k,ia,ib,itop1
      real*8 x
c
      itop1=itop
      do ! bubble sort
        break=.true.
        do i=2,itop1
          ia=results(1,i-1)
          ib=results(1,i)
          if (ia>ib) then
            break=.false.
            x=results(1,i)
            results(1,i)=results(1,i-1)
            results(1,i-1)=x
            x=results(2,i)
            results(2,i)=results(2,i-1)
            results(2,i-1)=x
          endif
        enddo
        if (break) exit
        itop1=itop1-1
      enddo
c     do i=1,itop
c       write(6,'(A,F5.0,F25.15)') "---results: ",
c    *                             results(1,i),results(2,i)
c     enddo
      itop1=itop+1
      do i=2,itop  ! squeeze
        ia=results(1,1)
        ib=results(1,i)
        if ((ib-ia).le.1) then
          results(1,1)=results(1,i)
          results(2,1)=results(2,1)+results(2,i)
          if (iprnt.gt.1) then
          write(6,'(A,F8.0,E30.16E2,I8)') 
     *    "Part. result: ",results(1,1), results(2,1),itop
          endif
          itop1=i+1
        else
          itop1=i
          exit
        endif
      enddo
      if (itop1.gt.2) then
        k=1
        do i=itop1,itop  !  move
          k=k+1
          results(1,k)=results(1,i)
          results(2,k)=results(2,i)
        enddo
        itop=k
      endif
      end
c===============================================================================

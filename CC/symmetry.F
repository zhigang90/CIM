      subroutine orb_pairs(nsym,     ncf,      nval,   nvirt,   fpair,
     *                     ivalpair, ivirpair, coeff,xoverlap)
c Search for symmetries in coeff array (one column = one MO orbital)
c It outputs ivalpair & ivirpair arrays, for each sym. operation
c and orbital number its orbital image is given. This is an MP2-like
c subroutine. It assumes that orbitals are already symmetrized, e.g. by
c adv_symmetrize subroutine.
c
c Parameters:
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of occupied (correlated) orbitals
c   nvirt    - number of virtual orbitals
c   fpair    - INTEGER(!) array of function pairs
c   ivalpair - occupied orbital images produced here
c   ivirpair - virtual  orbital images produced here
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c   xoverlap - overlap matrix
c
      implicit none
      integer nsym,ncf,nval,nvirt,fpair(7,ncf),ivalpair(7,nval)
      integer ivirpair(7,nvirt)
      real*8 coeff(ncf,nval+nvirt),xoverlap(ncf,ncf)
c
      integer indexx,ns,io,jo,j,k,ii,jj
      real*8 coeff_prim(ncf)
      real*8 overlap,overprim,xovcoe(ncf,ncf),xmax
      real*8, parameter :: eps = 1d-6
      logical done
c
      xmax=0d0
      if (nsym.eq.0) return
      do jj=1,nval+nvirt
        do ii=1,ncf
          xovcoe(ii,jj)=0d0
          do k=1,ncf
            xovcoe(ii,jj)=xovcoe(ii,jj)+xoverlap(ii,k)*coeff(k,jj)
          enddo
        enddo
      enddo
      do ns=1,nsym
        do io=1,nval
          do j=1,ncf
          indexx=fpair(ns,j)
          coeff_prim(abs(indexx))=coeff(j,io)
          if (indexx.lt.0) 
     *        coeff_prim(abs(indexx))=-coeff_prim(abs(indexx))
          enddo
          done=.false.
          do jo=1,nval
            overprim=0.0d0
            do j=1,ncf
              overprim=overprim+coeff_prim(j)*xovcoe(j,jo)
            enddo
            if (dabs(overprim).lt.1d-1) then
              if (dabs(overprim).gt.xmax) xmax=dabs(overprim)
            endif
            if (dabs(dabs(overprim)-1d0).lt.eps) then
              if (done) STOP 'Error 593h69fhso48f'
              if (overprim.lt.0.0d0) then 
                ivalpair(ns,io)=-jo
              else 
                ivalpair(ns,io)=jo
              endif
              done=.true.
            endif
          enddo
          if (.not.done) stop 'Error 594xx333h3hed8ff'
        enddo
      enddo
c-------------------------------------------------------------------------------
c Virtual part:
c
      do ns=1,nsym
        do io=1,nvirt
          do j=1,ncf
          indexx=fpair(ns,j)
          coeff_prim(abs(indexx))=coeff(j,io+nval)
          if (indexx.lt.0) 
     *        coeff_prim(abs(indexx))=-coeff_prim(abs(indexx))
          enddo
          done=.false.
          do jo=1,nvirt
            overprim=0.0d0
            do j=1,ncf
              overprim=overprim+coeff_prim(j)*xovcoe(j,jo+nval)
            enddo
            if (dabs(overprim).lt.1d-1) then
              if (dabs(overprim).gt.xmax) xmax=dabs(overprim)
            endif
            if (dabs(dabs(overprim)-1d0).lt.eps) then
              if (done) STOP 'Error 593h69fuuuhso48f'
              if (overprim.lt.0.0d0) then 
                ivirpair(ns,io)=-jo
              else 
                ivirpair(ns,io)=jo
              endif
              done=.true.
            endif
          enddo
          if (.not.done) then 
          done=.false.
cc          print *, 'Image of: ',io,',  under symmetry op. no: ',ns
          do jo=1,nvirt
            overprim=0.0d0
            do j=1,ncf
              overprim=overprim+coeff_prim(j)*xovcoe(j,jo+nval)
            enddo
cc            print 33,jo,overprim
 33         FORMAT (I4,3F20.10)
          enddo
            stop 'Error 594333h3hed8ff'
          endif
        enddo
      enddo
      end
c===============================================================================
      subroutine print_symmetries(nvirt, ncf)
      use memory
      implicit none
      integer nvirt,ncf
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer io,ns,int_array
c
      do io=1,nval
        print '(A,I2,A)',"Mirror of orbital i=",io,"  is:"
        do ns=1,nsym
          print '(A,I3)','         ',int_array(bl(ivalpair),ns+(io-1)*7)
        enddo
      enddo
      do io=1,nvirt
        print '(A,I2,A)',"Mirror of orbital i=",io,"  is:"
        do ns=1,nsym
          print '(A,I3)','         ',int_array(bl(ivirpair),ns+(io-1)*7)
        enddo
      enddo
      end
c===============================================================================
      subroutine pair_searcher(i,  j,     iprim,  jprim, ijprim,
     *                         ns, trans, signum)
      use memory
      implicit none
      integer i,j,iprim,jprim,ijprim,ns
      logical trans,signum
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
c
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer ij,iprim1,jprim1,ijprim1,isym,ix1,jx1,int_array
      logical tr1,tr2
c
      if (i.ge.j) then 
        ij=i*(i-1)/2+j
        tr1=.false.
      else
        ij=j*(j-1)/2+i
        tr1=.true.
      endif
c
      ijprim=ij
      if (nsym.eq.0.or.i.gt.nval.or.j.gt.nval) then
        iprim=i
        jprim=j
        ijprim=ij
        ns=0
        trans=tr1
        signum=.false.
        return
      endif
      do isym=1,nsym
        ix1=int_array(bl(ivalpair),isym+(i-1)*7)
        jx1=int_array(bl(ivalpair),isym+(j-1)*7)
        iprim1=abs(ix1)
        jprim1=abs(jx1)
        if (iprim1.ge.jprim1) then 
          ijprim1=iprim1*(iprim1-1)/2+jprim1
          tr2=.false.
        else
          ijprim1=jprim1*(jprim1-1)/2+iprim1
          tr2=.true.
        endif
        if (ijprim1.le.ij) then
          if (ijprim1.le.ijprim) then 
            iprim=iprim1
            jprim=jprim1
            ijprim=ijprim1
            ns=isym
            trans=tr2 ! changed...
            if (ix1*jx1.lt.0) then 
              signum=.true.
            else
              signum=.false.
            endif
          endif
        endif
      enddo
      if (ij.eq.ijprim) then
        iprim=i
        jprim=j
        ijprim=ij
        ns=0
        trans=tr1
        signum=.false.
      endif
c Chech whether the same symmetry operation turns prim pair into source pair
      if (ns.ne.0) then
       if (abs(int_array(bl(ivalpair),ns+(iprim-1)*7)).ne.i
     * .or.abs(int_array(bl(ivalpair),ns+(jprim-1)*7)).ne.j) then
cc      print '(F10.x)', i
        STOP 'I do not expect this might happen....'
       endif
      endif
      end
c===============================================================================
      subroutine orb_searcher(i, iprim, ns, signum)
c For given i orbital search for iprim orbital which is smaller (has lower
c number). If no symmetry or given orbital doesn't have an image with
c lower number, ns=0, signum=.false. and iprim=i
c Parameters:
c   i      - the number of orbital being checked
c   iprim  - the lowest image of orbital i
c   ns     - the symmetry operation transforming i to iprim (0 if no
c            symmetry is switched on or if i=iprim)
c   signum - does sign change when i->iprim, if yes, signum=.true.
      use memory
      implicit none
      integer i,iprim,ns
      logical signum
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
c
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer iprim1,isym,ix1,int_array
c
      iprim=i
      if (nsym.eq.0) then
        iprim=i
        ns=0
        signum=.false.
        return
      endif
      do isym=1,nsym
        ix1=int_array(bl(ivalpair),isym+(i-1)*7)
        iprim1=abs(ix1)
        if (iprim1.le.iprim) then
            iprim=iprim1
            ns=isym
            if (ix1.lt.0) then 
              signum=.true.
            else
              signum=.false.
            endif
        endif
      enddo
      if (i.eq.iprim) then
        iprim=i
        ns=0
        signum=.false.
      endif
c Chech whether the same symmetry operation turns prim pair into source pair
      if (ns.ne.0) then
       if (abs(int_array(bl(ivalpair),ns+(iprim-1)*7)).ne.i) then
cc      print '(F10.x)', i
        STOP 'I do not expect this might happen....'
       endif
      endif
      end
c===============================================================================
      subroutine symm_memory_allocator(ncf)
      implicit none
      integer ncf
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      call dynamic_getmem(ncf*ncf,iscr)
      end
c===============================================================================
      subroutine pair_s_initializer(vorb1,    nval1,    nsym1,ifpa,ifpb,
     *                            ivalpair1,ivirpair1,loca1,iorevtable1,
     *                            ivrevtable1,im_table1)
      implicit none
      integer nval1,nsym1,ifpa,ifpb,ivalpair1,ivirpair1
      integer iorevtable1,ivrevtable1,im_table1
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      logical vorb,vorb1,loca,loca1
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      nval       = nval1
      nsym       = nsym1
      ifp        = ifpa
      ifp1       = ifpb
      ivalpair   = ivalpair1
      ivirpair   = ivirpair1
      vorb       = vorb1
      loca       = loca1
      iorevtable = iorevtable1
      ivrevtable = ivrevtable1
      im_table   = im_table1
      end
c===============================================================================
      subroutine full_pair_searcher(i,  j,     iprim, jprim, ijprim,
     *                              ns, signum)
      use memory
      implicit none
      integer i,j,nval,nsym,ifp,ivalpair,iprim,jprim,ijprim,ns,ivirpair
      integer iscr,ifp1
      logical signum,vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer ij,iprim1,jprim1,ijprim1,isym,ix1,jx1,int_array
      integer indi,indj
c
      ij=(i-1)*nval+j
      ijprim=ij
      if (nsym.eq.0) then
        iprim=i
        jprim=j
        ijprim=ij
        ns=0
        signum=.false.
        return
      endif
      do isym=1,nsym
        indi=isym+7*(i-1)
        indj=isym+7*(j-1)
        ix1=int_array(bl(ivalpair),indi)
        jx1=int_array(bl(ivalpair),indj)
        iprim1=abs(ix1)
        jprim1=abs(jx1)
        ijprim1=nval*(iprim1-1)+jprim1
        if (ijprim1.le.ij) then
          if (ijprim1.le.ijprim) then 
            iprim=iprim1
            jprim=jprim1
            ijprim=ijprim1
            ns=isym
            if (ix1*jx1.lt.0) then 
              signum=.true.
            else
              signum=.false.
            endif
          endif
        endif
      enddo
      if (ij.eq.ijprim) then
        iprim=i
        jprim=j
        ijprim=ij
        ns=0
        signum=.false.
      endif
c Chech whether the same symmetry operation turns prim pair into source pair
      if (ns.ne.0) then
       if (abs(int_array(bl(ivalpair),ns+(iprim-1)*7)).ne.i.or.
     *     abs(int_array(bl(ivalpair),ns+(jprim-1)*7)).ne.j)
     *  STOP 'I do not expect this might happen....'
      endif
      end
c===============================================================================
c This is wrapper for c_symm_transform_ao, used because of scratch
      subroutine symm_transform_ao(idimen,  ns, trans1, signum, xmat)
      use memory
      implicit none
      logical trans1,signum
      integer idimen,ns
      real*8 xmat(idimen,idimen)
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      call c_symm_transform_ao(idimen,ns,trans1,signum,xmat,bl(iscr))
      end
c===============================================================================
      subroutine c_symm_transform_ao(idimen,  ns, trans1, signum, xmat,
     *                             scratch)
      use memory
      implicit none
      logical trans1,signum
      integer idimen,ns
      real*8 xmat(idimen,idimen),scratch(idimen,idimen)
c
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,j,isig,jsig,iprim,jprim,int_array
      logical sig
c
c I use sig.neqv.signum as an equivalent to the 'excusive or' operator
      if (ns.ne.0) then
        do i=1,idimen
          iprim=int_array(bl(ifp),ns+(i-1)*7)
          do j=1,idimen
            jprim=int_array(bl(ifp),ns+(j-1)*7)
            sig=.false.
            if (iprim*jprim.lt.0) sig=.true.
            if (sig.neqv.signum) then 
              scratch(abs(iprim),abs(jprim))=-xmat(i,j)
            else
              scratch(abs(iprim),abs(jprim))= xmat(i,j)
            endif
          enddo
        enddo
      endif
c
      if (ns.ne.0) then
        if (trans1) then
          do i=1,idimen
            do j=1,idimen
              xmat(i,j)=scratch(j,i)
            enddo
          enddo
        else
          call tfer(scratch,xmat,idimen*idimen)
        endif
      else
        if (trans1) call trans(xmat,idimen)
      endif
      end
c===============================================================================
      subroutine symm_trans_vec_ao_cop(idimen, ns, signum, source,dest)
      use memory
      implicit none
      logical signum
      integer idimen,ns
      real*8 source(idimen),dest(idimen)
c
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,isig,iprim,int_array
      logical sig
c
c I use sig.neqv.signum as an equivalent to the 'excusive or' operator
      if (ns.ne.0) then ! only if i!=iprim and symmetry exists
        do i=1,idimen
          iprim=int_array(bl(ifp),ns+(i-1)*7)
          sig=.false.
          if (iprim.lt.0) sig=.true.
          if (sig.neqv.signum) then 
            dest(abs(iprim))=-source(i)
          else
            dest(abs(iprim))= source(i)
          endif
        enddo
      endif
      end
c===============================================================================
      subroutine symm_trans_vec(idimen, ns, signum, xmat)
      use memory
      implicit none
      logical signum
      integer idimen,ns
      real*8 xmat(idimen),scratch(idimen)
c
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,isig,iprim,int_array
      logical sig
c
c I use sig.neqv.signum as an equivalent to the 'excusive or' operator
      if (vorb) then
        if (ns.ne.0) then
          do i=1,idimen
            isig=int_array(bl(ivirpair),ns+(i-1)*7)
            sig=.false.
            if (isig.lt.0) sig=.true.
            if (sig.neqv.signum) xmat(i)=-xmat(i)
          enddo
        endif
      else
        if (ns.ne.0) then
          do i=1,idimen
            iprim=int_array(bl(ifp),ns+(i-1)*7)
            sig=.false.
            if (iprim.lt.0) sig=.true.
            if (sig.neqv.signum) then 
              scratch(abs(iprim))=-xmat(i)
            else
              scratch(abs(iprim))= xmat(i)
            endif
          enddo
        endif
        if (ns.ne.0) then
          call tfer(scratch,xmat,idimen)
        endif
      endif
      end
c===============================================================================
      subroutine symm_transform(idimen,  ns, trans1, signum, xmat)
      use memory
      implicit none
      logical trans1,signum
      integer idimen,ns
      real*8 xmat(idimen,idimen)
c
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      call c_symm_transform(idimen,  ns, trans1, signum, xmat,bl(iscr))
      end
c===============================================================================
      subroutine c_symm_transform(idimen,  ns, trans1, signum, xmat,
     *                          scratch)
      use memory
      implicit none
      logical trans1,signum
      integer idimen,ns
      real*8 xmat(idimen,idimen),scratch(idimen,idimen)
c
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,j,isig,jsig,iprim,jprim,int_array
      logical sig
c
c I use sig.neqv.signum as an equivalent to the 'excusive or' operator
      if (vorb) then
        if (ns.ne.0) then
          do i=1,idimen
            isig=int_array(bl(ivirpair),ns+(i-1)*7)
            do j=1,idimen
              jsig=int_array(bl(ivirpair),ns+(j-1)*7)
              sig=.false.
              if (isig*jsig.lt.0) sig=.true.
              if (sig.neqv.signum) xmat(i,j)=-xmat(i,j)
            enddo
          enddo
        endif
      if (trans1) then 
        call trans(xmat,idimen)
      endif
      else
        if (ns.ne.0) then
          do i=1,idimen
            iprim=int_array(bl(ifp),ns+(i-1)*7)
            do j=1,idimen
              jprim=int_array(bl(ifp),ns+(j-1)*7)
              sig=.false.
              if (iprim*jprim.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                scratch(abs(iprim),abs(jprim))=-xmat(i,j)
              else
                scratch(abs(iprim),abs(jprim))= xmat(i,j)
              endif
            enddo
          enddo
        endif
      if (ns.ne.0) then
        if (trans1) then
          do i=1,idimen
            do j=1,idimen
              xmat(i,j)=scratch(j,i)
            enddo
          enddo
        else
          call tfer(scratch,xmat,idimen*idimen)
        endif
      else
        if (trans1) call trans(xmat,idimen)
      endif
      endif
      end
c===============================================================================
      subroutine full_symm_transform(idimen, ns, signum, xmat)
      use memory
      implicit none
      integer idimen,ns
      logical signum
      real*8 xmat(idimen,idimen)
c
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      call c_full_symm_transform(idimen, ns, signum, xmat, bl(iscr))
      end
c===============================================================================
      subroutine c_full_symm_transform(idimen,ns,signum,xmat,scratch)
      use memory
      implicit none
      logical vorb,signum,loca
      integer idimen,ifp,ivirpair,ns,ivalpair,nsym,nval,iscr,ifp1
      real*8 xmat(idimen,idimen),scratch(idimen,idimen)
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,j,isig,jsig,iprim,jprim,int_array
      logical sig
      integer ind(idimen)
c
      do i=1,idimen
        ind(i)=ns+(i-1)*7
      enddo
      if (vorb) then
        if (ns.ne.0) then
          do i=1,idimen
            isig=int_array(bl(ivirpair),ind(i))
            do j=1,idimen
              jsig=int_array(bl(ivirpair),ind(j))
              sig=.false.
              if (isig*jsig.lt.0) sig=.true.
              if (sig.neqv.signum) xmat(i,j)=-xmat(i,j)
            enddo
          enddo
        endif
      else
        if (ns.ne.0) then
          do i=1,idimen
            iprim=int_array(bl(ifp),ind(i))
            do j=1,idimen
              jprim=int_array(bl(ifp),ind(j))
              sig=.false.
              if (iprim*jprim.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                scratch(abs(iprim),abs(jprim))=-xmat(i,j)
              else
                scratch(abs(iprim),abs(jprim))= xmat(i,j)
              endif
            enddo
          enddo
        endif
      if (ns.ne.0) then
        call tfer(scratch,xmat,idimen*idimen)
      endif
      endif
      end
c===============================================================================
      subroutine full_symm_trans_cop(idimen, ns, signum, xmat, dest)
      use memory
      implicit none
      logical vorb,signum,loca
      integer idimen,ifp,ivirpair,ns,ivalpair,nsym,nval,iscr,ifp1
      real*8 xmat(idimen,idimen),dest(idimen,idimen)
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,j,isig,jsig,iprim,jprim,int_array
      logical sig
      integer ind(idimen)
c
      do i=1,idimen
        ind(i)=ns+(i-1)*7
      enddo
      if (vorb) then
        if (ns.ne.0) then
          do i=1,idimen
            isig=int_array(bl(ivirpair),ind(i))
            do j=1,idimen
              jsig=int_array(bl(ivirpair),ind(j))
              sig=.false.
              if (isig*jsig.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                dest(i,j)=-xmat(i,j)
              else
                dest(i,j)= xmat(i,j)
              endif
            enddo
          enddo
        endif
      else
        if (ns.ne.0) then
          do i=1,idimen
            iprim=int_array(bl(ifp),ind(i))
            do j=1,idimen
              jprim=int_array(bl(ifp),ind(j))
              sig=.false.
              if (iprim*jprim.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                dest(abs(iprim),abs(jprim))=-xmat(i,j)
              else
                dest(abs(iprim),abs(jprim))= xmat(i,j)
              endif
            enddo
          enddo
        endif
      endif
      if (ns.eq.0) then
        call tfer(xmat,dest,idimen*idimen)
      endif
      end
c===============================================================================
      subroutine full_symm_tr_add_tr(idimen, ns, signum, xmat, dest)
c Performs symmetry transformation of xmat, transforms xmat and adds to
c dest
      use memory
      implicit none
      logical vorb,signum,loca
      integer idimen,ifp,ivirpair,ns,ivalpair,nsym,nval,iscr,ifp1
      real*8 xmat(idimen,idimen),dest(idimen,idimen)
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,j,isig,jsig,iprim,jprim,int_array
      logical sig
      integer ind(idimen)
c
      do i=1,idimen
        ind(i)=ns+(i-1)*7
      enddo
      if (vorb) then
        if (ns.ne.0) then
          do i=1,idimen
            isig=int_array(bl(ivirpair),ind(i))
            do j=1,idimen
              jsig=int_array(bl(ivirpair),ind(j))
              sig=.false.
              if (isig*jsig.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                dest(j,i)=dest(j,i)-xmat(i,j)
              else
                dest(j,i)=dest(j,i)+xmat(i,j)
              endif
            enddo
          enddo
        endif
      else
        if (ns.ne.0) then
          do i=1,idimen
            iprim=int_array(bl(ifp),ind(i))
            do j=1,idimen
              jprim=int_array(bl(ifp),ind(j))
              sig=.false.
              if (iprim*jprim.lt.0) sig=.true.
              if (sig.neqv.signum) then 
                dest(abs(jprim),abs(iprim))=dest(abs(jprim),abs(iprim))
     *          -xmat(i,j)
              else
                dest(abs(jprim),abs(iprim))=dest(abs(jprim),abs(iprim))
     *          +xmat(i,j)
              endif
            enddo
          enddo
        endif
      endif
      if (ns.eq.0) then
        call trans(xmat,idimen)
        call add(xmat,dest,idimen*idimen)
      endif
      end
c===============================================================================
      subroutine orb_pairs_check(idimen,coeff,eps)
      use memory
      implicit none
      integer idimen
c
      integer indexx,ns,io,jo,j,i
      real*8 coeff_prim(idimen),eps
      real*8 overlap,overprim
      logical done,vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,int_array,jot,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      real*8 coeff(idimen,nval),diffe,diffe1,dif_overlap,signum
      integer iii
c
      if (nsym.eq.0) return
      if (vorb) then
      do ns=1,nsym
        do io=1,nval
          overlap=0.0d0
          do j=1,idimen
          overlap=overlap+coeff(j,io)*coeff(j,io)
          enddo
          do j=1,idimen
          indexx=int_array(bl(ivirpair),ns+(j-1)*7)
          coeff_prim(abs(indexx))=coeff(j,io)
          if (indexx.lt.0) 
     *        coeff_prim(abs(indexx))=-coeff_prim(abs(indexx))
          enddo
          done=.false.
          jot=int_array(bl(ivalpair),ns+(io-1)*7)
          if (jot.ge.0) then
            signum=+1d0
          else
            signum=-1d0
          endif
          do jo=1,nval
            overprim=0.0d0
            diffe=0.0d0
            do j=1,idimen
            overprim=overprim+coeff(j,jo)*coeff_prim(j)
c           diffe1=dabs(dabs(coeff(j,jo))-dabs(coeff_prim(j)))
            diffe1=dabs((coeff(j,jo))-(coeff_prim(j))*signum)
            if (diffe1.gt.diffe) diffe=diffe1
            enddo
            dif_overlap=dabs(dabs(overprim)-overlap)
            if (dif_overlap.lt.eps) then
              if (abs(jot).ne.jo) STOP 'Error 5934jgmfo439rtjg'
              if (done) STOP 'Error 493skdjrh382w'
              done=.true.
              if (diffe.gt.1d-18)
     *        print '(A,E10.2E3,A,I3,A,I3,A,I3,A,E10.2E3)',
     *                'DIFFE: ',diffe,
     *                '  from: ',io,
     *                '  to:',jo,
     *                '  via: ',ns,'  overlap diff: ',dif_overlap
              if (dif_overlap.gt.1.e-10)
     *        print '(A,E10.2E3,A,I3,A,I3,A,I3)','overlap: ',dif_overlap
     *               ,'  from: ',io,
     *                '  to:',jo,
     *                '  via: ',ns
            endif
          enddo
          if (.not.done) then 
            print *, 'Orbital ',io, ' has no image under symmetry: ',ns
c
          do jo=1,nval
            overprim=0.0d0
            diffe=0.0d0
            do j=1,idimen
              overprim=overprim+coeff(j,jo)*coeff_prim(j)
              diffe1=dabs(dabs(coeff(j,jo))-dabs(coeff_prim(j)))
              if (diffe1.gt.diffe) diffe=diffe1
            enddo
            print '(A,E10.2E3,A,I3,A,I3,A,I3)','DIFFE: ',diffe,
     *              '  from: ',io,
     *              '  to:',jo,
     *              '  via: ',ns
            print *,'           ',jo,dabs(overlap)-dabs(overprim)
          enddo
c
            STOP 'wc-wc-wc'
          endif
        enddo
      enddo
      else
      do ns=1,nsym
        do io=1,nval
          overlap=0.0d0
          do j=1,idimen
          overlap=overlap+coeff(j,io)*coeff(j,io)
          enddo
          do j=1,idimen
          indexx=int_array(bl(ifp),ns+(j-1)*7)
          coeff_prim(abs(indexx))=coeff(j,io)
          if (indexx.lt.0) 
     *        coeff_prim(abs(indexx))=-coeff_prim(abs(indexx))
          enddo
          done=.false.
          jot=int_array(bl(ivalpair),ns+(io-1)*7)
          if (jot.ge.0) then
            signum=+1d0
          else
            signum=-1d0
          endif
          do jo=1,nval
            overprim=0.0d0
            diffe=0.0d0
            do j=1,idimen
            overprim=overprim+coeff(j,jo)*coeff_prim(j)
c           diffe1=dabs(dabs(coeff(j,jo))-dabs(coeff_prim(j)))
            diffe1=dabs((coeff(j,jo))-(coeff_prim(j))*signum)
            if (diffe1.gt.diffe) diffe=diffe1
            enddo
            dif_overlap=dabs(dabs(overprim)-overlap)
            if (dif_overlap.lt.eps) then
              if (abs(jot).ne.jo) STOP 'Error 5934jgmfo439rtjg'
              if (done) STOP 'Error 493skdjrh382w'
              done=.true.
              if (diffe.gt.1d-18)
     *        print '(A,E10.2E3,A,I3,A,I3,A,I3,A,E10.2E3)',
     *                'DIFFE: ',diffe,
     *                '  from: ',io,
     *                '  to:',jo,
     *                '  via: ',ns,'  overlap diff: ',dif_overlap
              if (dif_overlap.gt.1.e-10)
     *        print '(A,E10.2E3,A,I3,A,I3,A,I3)','overlap: ',dif_overlap
     *               ,'  from: ',io,
     *                '  to:',jo,
     *                '  via: ',ns
            endif
          enddo
          if (.not.done) then 
            print *, 'Orbital ',io, ' has no image under symmetry: ',ns
c
          do jo=1,nval
            overprim=0.0d0
            diffe=0.0d0
            do j=1,idimen
              overprim=overprim+coeff(j,jo)*coeff_prim(j)
              diffe1=dabs(dabs(coeff(j,jo))-dabs(coeff_prim(j)))
              if (diffe1.gt.diffe) diffe=diffe1
            enddo
            print '(A,E10.2E3,A,I3,A,I3,A,I3)','DIFFE: ',diffe,
     *              '  from: ',io,
     *              '  to:',jo,
     *              '  via: ',ns
            print *,'           ',jo,dabs(overlap)-dabs(overprim)
          enddo
c
            STOP 'wc-wc-wc'
          endif
        enddo
      enddo
      endif
      end
c===============================================================================
      subroutine symmetrize_the_same_coeff(ncf,coeff)
      use memory
      implicit none
      integer ncf
c
      integer ns,io,jo,ii,iiprim
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,int_array,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      real*8 coeff(ncf,nval),diff,oldcoef
c
      if (nsym.eq.0) return
      do io=1,nval
        do ns=1,nsym
          jo=int_array(bl(ivalpair),ns+(io-1)*7)
          if (abs(jo).eq.io) then
            do ii=1,ncf
              iiprim=int_array(bl(ifp),ns+(ii-1)*7)
              if (abs(iiprim).gt.ii) then
                oldcoef=coeff(abs(iiprim),io)
                coeff(abs(iiprim),io)=coeff(ii,io)
              if (iiprim*jo.lt.0)  coeff(abs(iiprim),io)=
     *                            -coeff(abs(iiprim),io)
              diff=dabs(oldcoef)-dabs(coeff(abs(iiprim),io))
              if (diff.gt.1d-8) print *,'Strange: ',diff,ii,iiprim
              endif
            enddo
          endif
        enddo
      enddo
      end
c===============================================================================
      subroutine symmetrize_the_same_coeff_vir(ncf,nvir,coeff)
      use memory
      implicit none
      integer ncf,nvir
c
      integer ns,io,jo,ii,iiprim
      logical vorb,loca
      integer nval,nsym,ifp,ivalpair,ivirpair,int_array,iscr,ifp1
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
      real*8 coeff(ncf,nval),diff,oldcoef
c
      if (nsym.eq.0) return
      do io=1,nvir
        do ns=1,nsym
          jo=int_array(bl(ivirpair),ns+(io-1)*7)
          if (abs(jo).eq.io) then
            do ii=1,ncf
              iiprim=int_array(bl(ifp),ns+(ii-1)*7)
              if (abs(iiprim).gt.ii) then
                oldcoef=coeff(abs(iiprim),io)
                coeff(abs(iiprim),io)=coeff(ii,io)
              if (iiprim*jo.lt.0)  coeff(abs(iiprim),io)=
     *                            -coeff(abs(iiprim),io)
              diff=dabs(oldcoef)-dabs(coeff(abs(iiprim),io))
              if (diff.gt.1d-8) print *,'Strange: ',diff,ii,iiprim
              endif
            enddo
          endif
        enddo
      enddo
      end
c===============================================================================
      subroutine symmetrizer(nsym,     ncf,      nval,   nvirt,   fpair,
     *                       ivalpair, ivirpair, coeff)
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of correlated orbitals
c   nvirt    - number of virtual orbitals
c   fpair    - INTEGER(!) array of function pairs
c   ivalpair - occupied orbital images produced here
c   ivirpair - virtual  orbital images produced here
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c
      implicit none
      integer nsym,ncf,nval,nvirt,fpair(7,ncf),ivalpair(7,nval)
      integer ivirpair(7,nvirt)
      real*8 coeff(ncf,nval+nvirt)
c
      integer indexx,ns,io,jo,j
      real*8 coeff_prim(ncf),xmax
      real*8 overlap,overprim,signum,diff,s,biggest_diff,diffe,diffe1
      real*8, parameter :: eps = 1d-10, eps1 = 1.11023d-16
      logical done,symm
      integer iii,jprim
c
      do iii=1,100
c
      biggest_diff=0d0
c
      symm=.false.
      if (nsym.eq.0) return
      do ns=1,nsym
        do io=1,nval
          overlap=0.0d0
          do j=1,ncf
          overlap=overlap+coeff(j,io)*coeff(j,io)
          enddo
          jo=ivalpair(ns,io)
c         do jo=1,nval
          overprim=0.0d0
          signum=sign(1d0,dble(jo))
          diffe=0d0
          do j=1,ncf
            jprim=fpair(ns,j)
            overprim=overprim+coeff(j,abs(jo))
     *               *coeff(abs(jprim),io)*signum*sign(1d0,dble(jprim))
            diffe1=dabs(coeff(j,abs(jo))-coeff(abs(jprim),io)*signum*
     *                  sign(1d0,dble(jprim)))
          if (diffe1.gt.diffe) diffe=diffe1
          enddo
          print '(A,E20.10E3)','DDIIFFEE: ',diffe
          diff=dabs(dabs(overprim)-overlap)
          if (diff.gt.eps) then
            print *,diff,eps
            print *,overprim,overlap
            STOP 'Should not happen'
          endif
          if (diffe.gt.eps1) then  ! symmetrize
            print *,io,jo,diff
            symm=.true.
            do j=1,ncf
              jprim=fpair(ns,j)
              s=signum*sign(1d0,dble(jprim))
              coeff(j,io)=(coeff(abs(jprim),abs(jo))*s+coeff(j,io))/2d0
            enddo
            if (io.ne.abs(jo)) then ! backward transformation for jo:
              do j=1,ncf
                jprim=fpair(ns,j)
                coeff(abs(jprim),abs(jo))=coeff(j,io)
     *                                   *sign(1d0,dble(jprim))*signum
              enddo
            endif
          endif
c         enddo do jo=1,nval
        enddo
      enddo
c
      print *,'Symmetrization cycle ',iii
      if (.not.symm) exit
c
      enddo
c
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if (symm) stop 'Cannot establish symmetry exact images'
c-------------------------------------------------------------------------------
c Virtual part:
c
      do iii=1,100
c
      biggest_diff=0d0
c
      symm=.false.
      if (nsym.eq.0) return
      do ns=1,nsym
        do io=1,nvirt
          overlap=0.0d0
          do j=1,ncf
          overlap=overlap+coeff(j,io+nval)*coeff(j,io+nval)
          enddo
          jo=ivirpair(ns,io)
          if (abs(jo).ne.io) STOP 'dssfsdf98dsfsd8f73j4'
c
          overprim=0.0d0
          signum=sign(1d0,dble(jo))
          diffe=0d0
          xmax=0d0
          do j=1,ncf
            jprim=fpair(ns,j)
            overprim=overprim+coeff(j,abs(jo)+nval)
     *               *coeff(abs(jprim),io+nval)*
     *                signum*sign(1d0,dble(jprim))
            diffe1=dabs(coeff(j,abs(jo)+nval)-
     *             coeff(abs(jprim),io+nval)*signum*
     *                  sign(1d0,dble(jprim)))
          if (diffe1.gt.diffe) diffe=diffe1
          if (dabs(coeff(j,abs(jo)+nval)).gt.xmax)
     *                xmax=dabs(coeff(j,abs(jo)+nval))
          enddo
          print '(A,E20.10E3)','VIRT DDIIFFEE: ',diffe
          diff=dabs(dabs(overprim)-overlap)
          if (diff.gt.eps) then
            print *,diff,eps
            print *,overprim,overlap
            STOP 'Should not happen'
          endif
          if (diffe/xmax.gt.eps1*1d0) then  ! symmetrize
            print *,io,jo,diff,xmax
            symm=.true.
            do j=1,ncf
              jprim=fpair(ns,j)
              s=signum*sign(1d0,dble(jprim))
              coeff(j,io+nval)=(coeff(abs(jprim),abs(jo)+nval)*
     *                          s+coeff(j,io+nval))/2d0
            enddo
          endif
c         enddo do jo=1,nval
        enddo
      enddo
c
      print *,'Symmetrization cycle ',iii
      if (.not.symm) exit
c
      enddo
c
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
c     if (symm) stop 'Cannot establish symmetry exact images'
      end
c
c===============================================================================
      subroutine check_operator_symmetry(nsym,fpair,idimen,xmat)
      integer nsym,fpair(7,*),idimen
      real*8 xmat (idimen,idimen)
c
      integer i,j,ipri,jprim,sgn1,sgn2,iii
      real*8 diff,xmax
      if (nsym.eq.0) return
      xmax=0d0
      do iii=1,3
      do ns=1,nsym
        do i=1,idimen
          iprim=fpair(ns,i)
          sgn1=sign(1d0,dble(iprim))
          do j=1,idimen
            jprim=fpair(ns,j)
            sgn2=sign(1d0,dble(jprim))
            diff=(xmat(i,j)-xmat(abs(iprim),abs(jprim))
     *                           *sgn1*sgn2)
            if (abs(diff).gt.xmax) xmax=abs(diff)
            if (dabs(diff).gt.5d-16) then
              xmat(i,j)=(xmat(i,j)+xmat(abs(iprim),abs(jprim))
     *                             *sgn1*sgn2)/2d0
              xmat(abs(iprim),abs(jprim))=xmat(i,j)*sgn1*sgn2
            endif
          enddo
        enddo
      enddo
      enddo
cc      print *,'Max deviation of Fock matrix was: ',xmax
      xmax=0d0
      do ns=1,nsym
        do i=1,idimen
          iprim=fpair(ns,i)
          sgn1=sign(1d0,dble(iprim))
          do j=1,idimen
            jprim=fpair(ns,j)
            sgn2=sign(1d0,dble(jprim))
            diff=(xmat(i,j)-xmat(abs(iprim),abs(jprim))
     *                           *sgn1*sgn2)
            if (abs(diff).gt.xmax) xmax=abs(diff)
          enddo
        enddo
      enddo
cc      print *,'Max deviation of Fock matrix is: ',xmax
      call flush(6)
      end
c==============================================================================
      subroutine ampl_symm(ncf, coeff)
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of correlated orbitals
c   nvirt    - number of virtual orbitals
c   fpair    - INTEGER(!) array of function pairs
c   ivalpair - occupied orbital images produced here
c   ivirpair - virtual  orbital images produced here
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c
      use memory
      implicit none
      integer ncf
      real*8 coeff(ncf,*)
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer indexx,ns,io,jo,j
      real*8 coeff_prim(ncf),xmax
      real*8 overlap,overprim,signum,diff,s,biggest_diff,diffe,diffe1
      real*8, parameter :: eps = 1d-8, eps1 = 1.11023d-16
      logical done,symm
      integer iii,jprim,int_array
c
      if (vorb) stop 'Not ready yet for it'
      do iii=1,100
c
      biggest_diff=0d0
c
      symm=.false.
      if (nsym.eq.0) return
      do ns=1,nsym
        do io=1,nval
          overlap=0.0d0
          do j=1,ncf
          overlap=overlap+coeff(j,io)*coeff(j,io)
          enddo
c         jo=ivalpair(ns,io)
          jo=int_array(bl(ivalpair),ns+(io-1)*7)
c
          overprim=0.0d0
          signum=sign(1d0,dble(jo))
          diffe=0d0
          do j=1,ncf
c           jprim=fpair(ns,j)
            jprim=int_array(bl(ifp),ns+(j-1)*7)
            overprim=overprim+coeff(j,abs(jo))
     *               *coeff(abs(jprim),io)*signum*sign(1d0,dble(jprim))
            diffe1=dabs(coeff(j,abs(jo))-coeff(abs(jprim),io)*signum*
     *                  sign(1d0,dble(jprim)))
          if (diffe1.gt.diffe) diffe=diffe1
          enddo
          print '(A,E20.10E3)','DDIIFFEE: ',diffe
          diff=dabs(dabs(overprim)-overlap)
          if (diff.gt.eps) then
            print *,diff,eps
            print *,overprim,overlap
            STOP 'Should not happen'
          endif
          if (.false..and.diffe.gt.eps1) then  ! symmetrize
            print *,'Symmetrize: ',io,jo,diff
            symm=.true.
            do j=1,ncf
c             jprim=fpair(ns,j)
              jprim=int_array(bl(ifp),ns+(j-1)*7)
              s=signum*sign(1d0,dble(jprim))
              coeff(j,io)=(coeff(abs(jprim),abs(jo))*s+coeff(j,io))/2d0
            enddo
            if (io.ne.abs(jo)) then ! backward transformation for jo:
              do j=1,ncf
c               jprim=fpair(ns,j)
                jprim=int_array(bl(ifp),ns+(j-1)*7)
                coeff(abs(jprim),abs(jo))=coeff(j,io)
     *                                   *sign(1d0,dble(jprim))*signum
              enddo
            endif
          endif
c
        enddo
      enddo
c
      print *,'Symmetrization cycle ',iii
      if (.not.symm) exit
c
      enddo
c
c
c !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      STOP 'tmp spr'
      if (symm) stop 'Cannot establish symmetry exact images'
      end
c==============================================================================
      subroutine check_u_u_tr(nsym,fpair,idimen,xmat)
      integer nsym,fpair(7,*),idimen
      real*8 xmat (idimen,idimen)
      real*8 dev (idimen,idimen)
c
      integer i,j,ipri,jprim,sgn1,sgn2,iii
      real*8 diff,xmax
      if (nsym.eq.0) return
      do i=1,idimen
        do j=1,idimen
          dev(i,j)=0d0
        enddo
      enddo
      do ns=1,nsym
        do i=1,idimen
          iprim=fpair(ns,i)
          sgn1=sign(1d0,dble(iprim))
          do j=1,idimen
            jprim=fpair(ns,j)
            sgn2=sign(1d0,dble(jprim))
            diff=(xmat(i,j)-xmat(abs(iprim),abs(jprim))
     *                           *sgn1*sgn2)
            if (abs(diff).gt.dev(i,j)) dev(i,j)=abs(diff)
          enddo
        enddo
      enddo
      do i=1,idimen
        do j=1,idimen
          write(99,'(2I5,E20.10E3)') j,i,dev(j,i)
        enddo
          write(99,'(A)') '========================================'
      enddo
      end
c===============================================================================
      subroutine symmetrizer2(nsym,  ncf,  nval,  fpair,  ivalpair,
     *                        coeff)
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of correlated orbitals
c   nvirt    - number of virtual orbitals
c   fpair    - INTEGER(!) array of function pairs
c   ivalpair - occupied orbital images produced here
c   ivirpair - virtual  orbital images produced here
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c
      implicit none
      integer nsym,ncf,nval,fpair(7,*),ivalpair(7,*)
      real*8 coeff(ncf,*)
c
      logical done(ncf,ncf)
      real*8 srednia,suma,signum
      integer iii,jjj,io,j,iop,jp,ns
      real*8, parameter :: eps = 1d-8, eps1 = 1.11023d-16
c
      if (nsym.eq.0.or.nval.eq.0) return
      do iii=1,ncf
        do jjj=1,ncf
          done(iii,jjj)=.false.
        enddo
      enddo
      do io=1,nval
        do j=1,ncf
          if (done(j,io)) cycle
          suma=coeff(j,io)
          do ns=1,nsym
            if (ivalpair(ns,io)*fpair(ns,j).lt.0d0) signum=-1d0
            if (ivalpair(ns,io)*fpair(ns,j).ge.0d0) signum=+1d0
            iop=abs(ivalpair(ns,io))
            jp =abs(fpair(ns,j))
            suma=suma+coeff(jp,iop)*signum
            if (done(jp,iop)) STOP 'Weird'
          enddo
          srednia=suma/dble(nsym+1)
          if (dabs(srednia).lt.eps1) srednia=0d0
          coeff(j,io)=srednia
          done(j,io)=.true.
          do ns=1,nsym
            if (ivalpair(ns,io)*fpair(ns,j).lt.0d0) signum=-1d0
            if (ivalpair(ns,io)*fpair(ns,j).ge.0d0) signum=+1d0
            iop=abs(ivalpair(ns,io))
            jp =abs(fpair(ns,j))
            coeff(jp,iop)=srednia*signum
            done(jp,iop)=.true.
          enddo
        enddo
      enddo
c
      end
c
c===============================================================================
      subroutine core_pairs(nsym,   ncf,  nval,   fpair,  ivalpair,
     *                      coeff,xoverlap)
c Search for symmetries in coeff array (column = orbital coeff)
c It outputs ivalpair & ivirpair arrays, for each of nsym operation
c and orbital number is orbital image given
c Parameters:
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of correlated orbitals
c   nvirt    - number of virtual orbitals
c   fpair    - INTEGER(!) array of function pairs
c   ivalpair - occupied orbital images produced here
c   ivirpair - virtual  orbital images produced here
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c
      implicit none
      integer nsym,ncf,nval,fpair(7,ncf),ivalpair(7,nval)
      real*8 coeff(ncf,*),xoverlap(ncf,ncf)
c
      integer indexx,ns,io,jo,j,k,ii,jj
      real*8 coeff_prim(ncf)
      real*8 overlap,overprim,xovcoe(ncf,nval)
      real*8, parameter :: eps = 1d-6
      logical done
c
      if (nsym.eq.0.or.nval.eq.0) return
      do jj=1,nval
        do ii=1,ncf
          xovcoe(ii,jj)=0d0
          do k=1,ncf
            xovcoe(ii,jj)=xovcoe(ii,jj)+xoverlap(ii,k)*coeff(k,jj)
          enddo
        enddo
      enddo
      do ns=1,nsym
        do io=1,nval
          do j=1,ncf
          indexx=fpair(ns,j)
          coeff_prim(abs(indexx))=coeff(j,io)
          if (indexx.lt.0) 
     *        coeff_prim(abs(indexx))=-coeff_prim(abs(indexx))
          enddo
          done=.false.
          do jo=1,nval
            overprim=0.0d0
          do j=1,ncf
            overprim=overprim+coeff_prim(j)*xovcoe(j,jo)
          enddo
            if (dabs(dabs(overprim)-1d0).lt.eps) then
              if (done) STOP 'Error 593rereh69fhso48f'
              if (overprim.lt.0.0d0) then 
                ivalpair(ns,io)=-jo
              else 
                ivalpair(ns,io)=jo
              endif
              done=.true.
            endif
          enddo
          if (.not.done) stop 'Error 59433e3h3hed8ff'
        enddo
      enddo
      end
C===============================================================================
      subroutine symmetrize_energies(ncf,nsym,nval,ivalpair,eps)
      implicit none
      integer ncf,nsym,nval,ivalpair(7,*)
      real*8 eps(*)
      logical done(ncf)
c
      integer ns,ip,i
      real*8 aver
c
      if(nval.eq.0) return
c
      do i=1,ncf
        done(i)=.false.
      enddo
      do i=1,nval
        if (done(i)) cycle
        aver=eps(i)
        do ns=1,nsym
          ip=abs(ivalpair(ns,i))
          if (done(ip)) STOP 'Do not expect this to ever happen'
          aver=aver+eps(ip)
        enddo
        aver=aver/(nsym+1)
        eps(i)=aver
        do ns=1,nsym
          ip=abs(ivalpair(ns,i))
          eps(ip)=aver
          done(ip)=.true.
        enddo
      enddo
      end
c===============================================================================
      subroutine transform_orbital(idimen,sig,ns,xorb)
      use memory
      implicit none
      integer idimen,ns
      real*8 xorb(*)
      logical sig
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer i,ip,int_array
      real*8 scr(idimen),signum,ips
c
      if (ns.eq.0) return
      if (sig) then
        signum=-1d0
      else
        signum=+1d0
      endif
      if (vorb) then
        STOP 'Not ready yet'
      else
        do i=1,idimen
          ip=int_array(bl(ifp),ns+(i-1)*7)
          ips=sign(1d0,dble(ip))
          ip=abs(ip)
          scr(ip)=xorb(i)*ips*signum
        enddo
        do i=1,idimen
          xorb(i)=scr(i)
        enddo
      endif
      end
c===============================================================================
      subroutine symm2_s_ampl(idimen, coeff)
      use memory
      implicit none
c
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer idimen
      real*8 coeff(idimen,*)
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer int_array,io,j,ind1,ind2,iop,jp,iii,jjj,ns
      real*8 suma,srednia,signum
      logical done(idimen,nval)
      real*8, parameter :: eps = 1d-8, eps1 = 1.11023d-16
c
      if (nsym.eq.0) return
      do iii=1,idimen
        do jjj=1,nval
          done(iii,jjj)=.false.
        enddo
      enddo
      if (vorb) then
        do io=1,nval
          do j=1,idimen
            if (done(j,io)) cycle
            suma=coeff(j,io)
            do ns=1,nsym
              ind1=int_array(bl(ivirpair),ns+(j-1)*7)
              ind2=int_array(bl(ivalpair),ns+(io-1)*7)
              if (ind2*ind1.lt.0d0) signum=-1d0
              if (ind2*ind1.ge.0d0) signum=+1d0
              iop=abs(ind2)
              jp =abs(ind1)
              suma=suma+coeff(jp,iop)*signum
              if (done(jp,iop)) STOP 'Weird'
            enddo
            srednia=suma/dble(nsym+1)
            if (dabs(srednia).lt.eps1) srednia=0d0
            coeff(j,io)=srednia
            done(j,io)=.true.
            do ns=1,nsym
              ind1=int_array(bl(ivirpair),ns+(j-1)*7)
              ind2=int_array(bl(ivalpair),ns+(io-1)*7)
              if (ind2*ind1.lt.0d0) signum=-1d0
              if (ind2*ind1.ge.0d0) signum=+1d0
              iop=abs(ind2)
              jp =abs(ind1)
              coeff(jp,iop)=srednia*signum
              done(jp,iop)=.true.
            enddo
          enddo
        enddo
      else
        do io=1,nval
          do j=1,idimen
            if (done(j,io)) cycle
            suma=coeff(j,io)
            do ns=1,nsym
              ind1=int_array(bl(ifp),ns+(j-1)*7)
              ind2=int_array(bl(ivalpair),ns+(io-1)*7)
              if (ind2*ind1.lt.0d0) signum=-1d0
              if (ind2*ind1.ge.0d0) signum=+1d0
              iop=abs(ind2)
              jp =abs(ind1)
              suma=suma+coeff(jp,iop)*signum
              if (done(jp,iop)) STOP 'Weird'
            enddo
            srednia=suma/dble(nsym+1)
            if (dabs(srednia).lt.eps1) srednia=0d0
            coeff(j,io)=srednia
            done(j,io)=.true.
            do ns=1,nsym
              ind1=int_array(bl(ifp),ns+(j-1)*7)
              ind2=int_array(bl(ivalpair),ns+(io-1)*7)
              if (ind2*ind1.lt.0d0) signum=-1d0
              if (ind2*ind1.ge.0d0) signum=+1d0
              iop=abs(ind2)
              jp =abs(ind1)
              coeff(jp,iop)=srednia*signum
              done(jp,iop)=.true.
            enddo
          enddo
        enddo
      endif
c
      end
c
c===============================================================================
      subroutine irr_rep_table(orbpairs,virtpairs,nsym,nval,nvirt,
     *                         icharac,char_count,m_table)
c Make tables of characters of irreducible representations from
c MO orbitals image arrays, created by orb_pairs subroutine.
c This means 'recovering' irreps from orbitals. The subroutine assumes
c existence two arrays for both occupied and virtual orbitals. Virtuals
c are nedded because (I guess) occupied orbitals may not represent all
c possible irrreps in this group. It may happen in small molecules.
c
c Parameters:
c     Input:
c orbpairs    - array with symmetry images of occupied MOs
c virtpairs   - array with symmetry images of virtual MOs
c nsym        - number of symmetry operations
c nval        - number of occupied MOs
c nvirt       - number of virtual MOs
c     Output:
c icharac     - character table
c char_count  - number of irreps (irreducible representations) found.
c m_table     - irreps multiplication table. Needed in CCSD, but not
c               in Fock part.
      implicit none
      integer orbpairs(7,*),virtpairs(7,*),nval,nvirt,icharac(7,8)
      integer nsym,char_count,m_table(8,8)
c
      integer characters(7),is,ic,i,i1,i2,ires,iprnt
      logical equal,characters_eq
c
      char_count=1
      if (nsym.eq.0) return
      call izeroit(icharac,7*8)
      call izeroit(m_table,64)
      do is=1,nsym       ! character=1 is fully symmetric
        icharac(is,1)=1
      enddo
c
      do i=1,nval
       call izeroit(characters,7)
       do is=1,nsym
         if (orbpairs(is,i).eq.i) then
         characters(is)=1
         else if (orbpairs(is,i).eq.-i) then
         characters(is)=-1
         else
         print *,orbpairs(is,i),i
          STOP 'Not a member of irreducible representation'
         endif
       enddo
       equal=.false.
       do ic=1,char_count
         if (characters_eq(icharac(1,ic),characters,nsym)) then
           equal=.true.
           exit
         endif
       enddo
       if (.not.equal) then
       char_count=char_count+1
       if (char_count.gt.8) STOP 'char_count.gt.8'
       do is=1,nsym
       icharac(is,char_count)=characters(is)
       enddo
       endif
      enddo
c
      do i=1,nvirt
       call izeroit(characters,7)
       do is=1,nsym
         if (virtpairs(is,i).eq.i) then
         characters(is)=1
         else if (virtpairs(is,i).eq.-i) then
         characters(is)=-1
         else
         print *,virtpairs(is,i),i
          STOP 'Not a member of irreducible representation'
         endif
       enddo
       equal=.false.
       do ic=1,char_count
         if (characters_eq(icharac(1,ic),characters,nsym)) then
           equal=.true.
           exit
         endif
       enddo
       if (.not.equal) then
       char_count=char_count+1
       if (char_count.gt.8) STOP 'char_count.gt.8'
       do is=1,nsym
       icharac(is,char_count)=characters(is)
       enddo
       endif
      enddo
c
      call getival('ccprint',iprnt)
      if (iprnt.gt.3) then
        write(6,*) 'Character tables:'
        do ic=1,char_count
          write(6,'(I10,7I3)') ic,
     *  icharac(1,ic),
     *  icharac(2,ic),
     *  icharac(3,ic),
     *  icharac(4,ic),
     *  icharac(5,ic),
     *  icharac(6,ic),
     *  icharac(7,ic)
        enddo
        write(6,*) ' '
      endif
      if (char_count.ne.nsym+1) then
        print *,'Number of irr. rep not eq. nsym+1'
        print *,char_count,nsym+1
        STOP 'Number of irr. rep not eq. nsym+1'
      endif
c Create characters multiplication table:
c
      do i1=1,char_count
        do i2=1,char_count
          call izeroit(characters,7)
          do is=1,7
            characters(is)=icharac(is,i1)*icharac(is,i2)
          enddo
          equal=.false.
          do ic=1,char_count
            if (characters_eq(icharac(1,ic),characters,nsym)) then
            ires=ic
            equal=.true.
            exit
            endif
          enddo
          if (.not.equal) STOP '394jfkak48thejk'
          m_table(i1,i2)=ires
        enddo
      enddo
      if (iprnt.gt.3) then
        write(6,*) '==============================================='
        do i1=1,char_count
          do i2=1,char_count
          write(6,'(I4)',ADVANCE='NO') m_table(i1,i2)
          enddo
          write(6,*)
        enddo
        write(6,*) '==============================================='
      endif
      end
c===============================================================================
      function characters_eq(a,b,nsym)
      implicit none
      integer a(*),b(*),nsym
      logical characters_eq
      integer i
      characters_eq=.true.
      do i=1,nsym
        if (a(i).ne.b(i)) then
          characters_eq=.false.
          exit
        endif
      enddo
      end
c===============================================================================
      subroutine sort_orbitals(orbpairs,nsym,norb,icharac,char_count,
     *                         energies,cano,scr,ncf)
      implicit none
      integer orbpairs(7,norb),nsym,norb,icharac(7,8)
      integer char_count,ncf
      real*8 cano(ncf,norb),scr(ncf,norb),energies(norb)
c
      integer characters(7),rep_list(norb),ic,i,last,is
      logical equal,characters_eq
c
      if (nsym.eq.0) return
      do i=1,norb
       call izeroit(characters,7)
       do is=1,nsym
         if (orbpairs(is,i).eq.i) then
         characters(is)=1
         else if (orbpairs(is,i).eq.-i) then
         characters(is)=-1
         else
          STOP 'Not a member of irreducible representation'
         endif
       enddo
       equal=.false.
       do ic=1,char_count
         if (characters_eq(icharac(1,ic),characters,nsym)) then
           equal=.true.
           rep_list(i)=ic
           exit
         endif
       enddo
       if (.not.equal) STOP '450923iuogrfj'
      enddo
      last=0
      do ic=1,char_count
        do i=1,norb
          if (rep_list(i).eq.ic) then
            last=last+1
            if (last.gt.norb) STOP '45939ejfgh34k'
            call tfer(cano(1,i),scr(1,last),ncf)
          endif
        enddo
      enddo
      call tfer(scr,cano,ncf*norb)
      last=0
      do ic=1,char_count
        do i=1,norb
          if (rep_list(i).eq.ic) then
            last=last+1
            if (last.gt.norb) STOP '45939ejfgh34k'
            scr(last,1)=energies(i)
          endif
        enddo
      enddo
      call tfer(scr,energies,norb)
      end
c===============================================================================
      subroutine tables_generator(nsym,norb,char_count,icharac,orbpairs,
     *                            table,revtable)
      implicit none
      integer norb,char_count,icharac(7,8),orbpairs(7,*),table(norb+1,8)
      integer revtable(norb),nsym
c
      integer characters(7),is,ic,i,icount,iprnt
      logical equal,characters_eq
c
c
      call izeroit(table,(norb+1)*8)
      call izeroit(revtable,norb)
      if (char_count.le.1) return
      do i=1,norb
        call zeroit(characters,7)
        do is=1,nsym
          if (orbpairs(is,i).eq.i) then
            characters(is)=1
          else if (orbpairs(is,i).eq.-i) then
            characters(is)=-1
          else
            STOP '93404odkfjghj403'
          endif
        enddo
        equal=.false.
        do ic=1,char_count
          if (characters_eq(icharac(1,ic),characters,nsym)) then
            revtable(i)=ic
            equal=.true.
            exit
          endif
        enddo
        if (.not.equal) STOP '340343utsdkfghdsf'
      enddo
      call getival('ccprint',iprnt)
      if (iprnt.gt.3) then
        do i=1,norb
          write(6,*) 'revtable: ',i,revtable(i)
        enddo
      endif
      call dynamic_lock(revtable,i)
c
c Now when we have revtable ready:
c
      do ic=1,char_count
        icount=1
        do i=1,norb
          if (revtable(i).eq.ic) then
          icount=icount+1
          table(icount,ic)=i
          endif
        enddo
        table(1,ic)=icount-1
        if ((icount-1).gt.norb) STOP '4923u5jgdbsbrbrt5555'
      enddo
      if (iprnt.gt.3) then
        write(*,*) '==================================================='
        do ic=1,char_count
          do i=1,min(norb+1,26)
            write(*,'(I3)',ADVANCE='NO') table(i,ic)
          enddo
          write(*,*)
        enddo
        write(*,*) '==================================================='
      endif
      end
c===============================================================================
      subroutine ijtab_generator(nsym,norb,char_count,icharac,orbpairs,
     *                            m_table,revsrc,table)
      implicit none
      integer norb,char_count,icharac(7,8),orbpairs(7,*)
      integer table(norb*norb+1,8),nsym,revsrc(norb),m_table(8,8)
c
      integer ic,i,j,ij,icount
      integer revtable(norb*norb),iprnt
c
c
      call izeroit(table,(norb*norb+1)*8)
      call izeroit(revtable,norb*norb)
      if (char_count.le.1) return
      ij=0
      do j=1,norb
        do i=1,norb
          ij=ij+1
          revtable(ij)=m_table(revsrc(i),revsrc(j))
        enddo
      enddo
c
      do ic=1,char_count
        icount=1
        ij=0
        do j=1,norb
          do i=1,norb
            ij=ij+1
            if (revtable(ij).eq.ic) then
            icount=icount+1
            table(icount,ic)=ij
            endif
          enddo
        enddo
        if ((icount-1).gt.norb*norb) STOP '4923u5jgdbsbrbrt5555'
        table(1,ic)=icount-1
      enddo
      call getival('ccprint',iprnt)
      if (iprnt.gt.3) then
        write(*,*) '==================================================='
        do ic=1,char_count
          do i=1,min(norb*norb+1,20)
            write(*,'(I4)',ADVANCE='NO') table(i,ic)
          enddo
          write(*,*)
        enddo
        write(*,*) '==================================================='
      endif
      end
c===============================================================================
      subroutine adv_symmetrize(nsym,   ncf,  nval,   fpair,
     *                          coeff,xoverlap,xovcoe,
     *                          coefprim,S,ene)
c Symmetrize orbitals so that they become one dimensional basis vectors for
c irreps of an Abelian group. They can be mixed in the input "coeff",
c because a higher symmetry is present in the system. This subroutine
c should be able to handle up to five orbitals mixed, although testing
c wasn't performed thoroughly for such systems with very high symmetry.
c
c Input parameters:
c   nsym     - number of symmetry operations, 0 means no symetry at all
c   ncf      - number of contracted basis functions
c   nval     - number of orbitals being symmetrized.
c   fpair    - INTEGER(!) array of basis functions symmetry pairs
c   coeff    - the HF coefficient matrix, but without frozen orbitals
c   xoverlap - overlap matrix
c   ene      - array of orbital energies.
c Scratch space:
c   xovcoe   - SC product (overlap*coeff)
c   coefprim - C` = Symmetry operation * C (coeff)
c   S        - the product C`SC will be held here
c
      implicit none
      integer nsym,ncf,nval,fpair(7,ncf)
      real*8 coeff(ncf,nval),xoverlap(ncf,ncf),xovcoe(ncf,nval)
      real*8 coefprim(ncf,nval),S(nval,nval),ene(nval)
c
      integer indexx,ns,io,jo,j,k,ii,jj,i,itry,ll,isize,kk,ns1
      real*8 xsum
      real*8 overprim(2,nval)
      integer, parameter :: max_block = 5
      real*8 values(max_block), SS(max_block*max_block)
      real*8 vectors(max_block*max_block),newco(ncf,max_block),el
      real*8, parameter :: eps = 1d-12 ! square of max nondiagonal element of S
      real*8 :: eps1 = eps
      logical done
      integer iz,nn,maxc,mm,icf,iover,ires,icoepr,ixov
      real*8 xmax,xval
c
      if (nsym.eq.0.or.nval.eq.0) return
c
c calculate the largest non-diagonal element:
      call pointer2bl(coeff,icf)
      call pointer2bl(xoverlap,iover)
      call pointer2bl(S,ires)
      call pointer2bl(xovcoe,ixov)
      call pointer2bl(coefprim,icoepr)
      call matconn('x1cof','r',ncf,nval,icf)
      call matconn('x1ove','r',ncf,ncf,iover)
      call matconn('x1part','r',nval,nval,ires)
      call matconn('xovcoe','r',ncf,nval,ixov)
      call matconn('xcoefpr','r',ncf,nval,icoepr)

c Build S*C matrix
      do ns=1,nsym
        call matmmul2('x1ove','x1cof','xovcoe','n','n','n')
c Build C' matrix, which is a symmetry image of C, under operation "ns"
        do io=1,nval
          do j=1,ncf
          indexx=fpair(ns,j)
          coefprim(abs(indexx),io)=coeff(j,io)
          if (indexx.lt.0) 
     *        coefprim(abs(indexx),io)=-coefprim(abs(indexx),io)
          enddo
        enddo
c Calculate the overlap C'SC
        call matmmul2('xcoefpr','xovcoe','x1part','t','n','n')
c Seek for non-diagonal submatrices to be diagonalized (tricky task)
c They are blocs along diagonal
      ii=0
      ww:do jj=1,nval
        ii=ii+1
        if (ii.ge.nval) exit
 270    continue
        maxc=ii               ! Max diagonal element connected via
c                               nodiagonal with ii, init to ii itself
        do nn=ii,min(ii+20,nval)     ! n goes over columns
          !investigate column nn, below diag. AND maxc
          do kk=max(nn+1,maxc),min(ii+20,nval)
            el=S(kk,nn)*S(kk,nn)
            if (el.gt.eps1) then 
              maxc=kk
            endif
          enddo
          if (maxc.le.nn) exit
          if (maxc.ge.ii+5) then 
            if (eps1.gt.1d-7) then 
              write(6,*)ii,maxc
              do mm=ii,min(maxc,nval)
                do kk=ii,min(maxc,nval)
                  write(6,'(E10.2E2)',ADVANCE='no') S(mm,kk)
                enddo
                write(6,*)
              enddo
              call nerror(1,'adv_symmetrize',
     $                    'Too large representation detected',ii,maxc)
            else
              print '(A,2I5,F15.10)',
     *         'Warninig, increasing eps due to rep size. ',
     *           ii,maxc,dabs(ene(maxc)-ene(ii))
              eps1=eps1*2d0
              goto 270
            endif
          endif
          if (dabs(ene(maxc)-ene(ii)).gt.1d-4) then
            print '(A,2I5,F15.10)',
     *         'Warninig, increasing eps due to energy dif. ',
     *           ii,maxc,dabs(ene(maxc)-ene(ii))
            eps1=eps1*2d0
            goto 270
          endif
        enddo
        if (maxc.eq.ii) cycle
        isize=maxc-ii+1
c
        if (isize.eq.1) STOP '5893kwskdngndsdj47'
        if (isize.gt.5) STOP '589eee3kwskdngndsdj47'
c
c Now we know that the block size is at least: isize
c
        print*, '============================================'
        write(*,'(A)',advance='NO') 'ENES: '
        write(*,'(I5)',advance='NO') isize
        write(*,'(F13.8)',advance='NO') ene(maxc)-ene(ii)
        do kk=ii,maxc
          write(*,'(F12.6)',advance='NO') ene(kk)
        enddo
        write(*,*)
        do kk=1,isize
          do ll=1,isize
            SS(kk+(ll-1)*isize)=S(ii+kk-1,ii+ll-1)
            write(*,'(F15.8)',advance='NO') SS(kk+(ll-1)*isize)
          enddo
          write(*,*)
        enddo
c Diagonalize the submatrix:
        call sdiag2(isize,isize,SS,values,vectors)
        write(*,*)
      print '(5F15.8)',values(1),values(2),values(3),values(4),values(5)
        print *, '==========================================='
c
c create new linear combination of orbitals:
c
        do i=1,isize
          do j=1,ncf
            newco(j,i)=0d0
            do k=1,isize
            newco(j,i)=newco(j,i)+coeff(j,ii+k-1)*vectors(k+(i-1)*isize)
            enddo
          enddo
        enddo
c
c Replace old orbitals by new ones
c
        do k=1,isize
          do j=1,ncf
            coeff(j,ii+k-1)=newco(j,k)
          enddo
        enddo
        ii=ii+isize-1
      enddo ww
      enddo
      call matdisc('xcoefpr')
      call matdisc('xovcoe')
      call matdisc('x1part')
      call matdisc('x1ove')
      call matdisc('x1cof')
cc      print *,'Largest eps1:  ',eps1
      end
C===============================================================================
      function image(i,j) ! if ij pair is an image of a pair before ij, 
c                         return .true., if there is no image with
c                         smaller number, return .false.
      use memory
      implicit none
      logical image
      integer i,j
c
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer iprim,jprim,ijprim,int_array,ns,ij
c
      image=.false.
      if (loca.or.nsym.eq.0) return
      if (i.ge.j) then 
        ij=i*(i-1)/2+j
      else
        ij=j*(j-1)/2+i
      endif
      do ns=1,nsym
        iprim=int_array(bl(ifp1),ns+(i-1)*7)
        jprim=int_array(bl(ifp1),ns+(j-1)*7)
        if (iprim.ge.jprim) then 
          ijprim=iprim*(iprim-1)/2+jprim
        else
          ijprim=jprim*(jprim-1)/2+iprim
        endif
        if (ij.gt.ijprim) then ! Only the first pair is unique
          image=.true.
          return
        endif
      enddo
      end
C===============================================================================
      subroutine images_list(i,j,mu,nu,list,ex3)
      use memory
      implicit none
      integer i,j,mu,nu
      integer list(3,10) ! 1: iprim, 2: jprim, 3: sign= +1 or -1
      logical ex3
c
      integer nval,nsym,ifp,ivalpair,ivirpair,iscr,ifp1
      logical vorb,loca
      integer iorevtable,ivrevtable,im_table
      common /ccsymm/ vorb,nval,nsym,ifp,ivalpair,ivirpair,iscr,loca,
     *                ifp1,iorevtable,ivrevtable,im_table
c
      integer ns,iprim,jprim,int_array,muprim,nuprim
c
      call izeroit(list,30)
      if (loca.or.nsym.eq.0) return
      do ns=1,nsym
        if (ex3) then
          iprim =int_array(bl(ivirpair),ns+(i-1)*7)
        else
          iprim =int_array(bl(ivalpair),ns+(i-1)*7)
        endif
        jprim =int_array(bl(ivalpair),ns+(j-1)*7)
        muprim=int_array(bl(ifp),ns+(mu-1)*7)
        nuprim=int_array(bl(ifp),ns+(nu-1)*7)
        list(1,ns)=abs(muprim)
        list(2,ns)=abs(nuprim)
        if (iprim*jprim*muprim*nuprim.lt.0) then 
          list(3,ns)=-1
        else
          list(3,ns)=+1
        endif
      enddo
      end
C===============================================================================
      subroutine check_S(S,ncf,fpair,nsym)
      implicit none
      integer ncf,fpair(7,ncf),nsym
      real*8 S(ncf,ncf)
c
      integer i,j,iprim,jprim,ns
      real*8 Sprim
c
      do ns=1,nsym
      do i=1,ncf
        iprim=fpair(ns,i)
        do j=1,ncf
          jprim=fpair(ns,j)
          Sprim=S(abs(iprim),abs(jprim))
          if (iprim*jprim<0) Sprim=-Sprim
          print '(F30.17)',Sprim-S(i,j)
        enddo
      enddo
      enddo
      end
C===============================================================================
      subroutine symm_adapted(ncf,  nsym, nrep,icharacter,fpair,
     *                        xmat, xmat1,ibassym)
c
c Create a matrix which transforms AO basis set into symmetry adapted basis.
c This matrix, at least for Abelian groups, must have its columns or
c rows mutually orthogonal (I have proof somewhere in my notes), but
c they are not normalized, so it is not an unitary matrix.
c
c Parameters:
c        Input:
c ncf        - number of basis functions
c nsym       - number of symmetries
c nrep       - number of irreducible representations in the group
c icharacter - character table
c fpair      - AOs' symmetry pairs array.
c        Output:
c xmat       - transforming matrix. It is sorted, so it produces symmetry
c              adapted AOs grouped by type of irrep they generate. This
c              way Fock and overlap matrix produced by similarity
c              transformation are in block-diagonal form.
c ibassym    - an array which indicates first and last symmetry adapted
c              orbital for a given irrep.
c              e.g. first orbital in irrep no 3: ibassym(1,3)
c                   last in the irrep 3:         ibassym(2,3)
c        Scratch:
c xmat1      - huge scratch matrix with size ncf x ncf*8
c
      implicit none
      integer ncf,nsym,nrep,icharacter(7,8),fpair(7,*),ibassym(2,8)
      real*8 xmat(ncf,ncf),xmat1(ncf,ncf*8)
c
      logical done(ncf)
      integer i,j,ik,ikp,iks,ic,ns,irep,istart,istop,iprnt
      real*8 xsum
c
      call zeroit(xmat,ncf*ncf)
      call izeroit(ibassym,16)
      if (nsym.eq.0) then
        do i=1,ncf
          xmat(i,i)=1d0
        enddo
        ibassym(1,1)=1
        ibassym(2,1)=ncf
        return
      endif
      call zeroit(xmat1,ncf*ncf*8)
c
      do i=1,ncf
        done(i)=.false.
      enddo
      ik=0
      do i=1,ncf
        if (done(i)) cycle
        do irep=1,nrep
          ik=ik+1
          do ns=0,nsym
            if (ns.gt.0) then
              ikp=abs(fpair(ns,i))
              iks=sign(1,fpair(ns,i))
              ic=icharacter(ns,irep)
            else
              ikp=i
              iks=1
              ic=1
            endif
            done(ikp)=.true.
            xmat1(ikp,ik)=xmat1(ikp,ik)+dble(iks*ic)
          enddo
        enddo
      enddo
c
      ik=0
      do irep=1,nrep
        do j=irep,ncf*8,nrep
          xsum=0d0
          do i=1,ncf
            xsum=xsum+dabs(xmat1(i,j))
          enddo
          if (xsum.lt.1d-5) cycle
          ik=ik+1
          ibassym(2,irep)=ik
          call tfer(xmat1(1,j),xmat(1,ik),ncf)
        enddo
      enddo
      istart=1
      do irep=1,nrep
        istop=ibassym(2,irep)
        if (istop.eq.0) cycle
        ibassym(1,irep)=istart
        istart=istop+1
      enddo
      call getival('ccprint',iprnt)
      if (iprnt.gt.3) then
        do irep=1,8
          write(6,'(3I5)')irep,ibassym(1,irep),ibassym(2,irep)
        enddo
        call flush(6)
      endif
      end
C===============================================================================
      subroutine canonical_symm_mult(xmat1, xmat2,  xmat3, id1, id2,
     *                               id3,   ir1,    ir2,   nvirt,nrep,
     *                               vtable,mtable, vorb,  op)
      use memory
      implicit none
      integer id1,id2,id3,ir1,ir2,nrep,nvirt,vtable(nvirt+1,*)
      integer mtable(8,8)
      logical vorb
      character*1 op
      real*8 xmat1(id1,id2),xmat2(id2,id3),xmat3(id1,id3)
c
      real*8 beta
      integer imat1,imat2,imat3,irep,irep_1,irep_3,isymdim1,isymdim2
      integer isymdim3,i,ib,k,kb,l,lb,itmp_index
c
      beta=-1d0
      if (op.eq.'a') beta=1d0
      if (op.eq.'n') beta=0d0
      if (beta.lt.-0.5d0) call abort
      if (nrep.le.1.or..not.vorb) then
        call dgemm('N',   ! Do not transpose xmat1 ('T' = transpose)
     *             'N',   ! Do not transpose xmat2
     *             id1,   ! Rows of xmat1 and xmat3
     *             id3,   ! Columns of xmat2 and xmat3
     *             id2,   ! Columns of xmat1 & rows of xmat2
     *             1.0d0, ! alpha of: C = alpha*op( A )*op( B ) + beta*C
     *             xmat1, ! A
     *             id1,   ! Leading dimension of A
     *             xmat2, ! B
     *             id2,   ! Leading dimension of B
     *             beta,  ! beta of: C = alpha*op( A )*op( B ) + beta*C
     *             xmat3, ! C
     *             id1)   ! Leading dimension of C
      else !  if (nrep.le.1.or..not.vorb) then
        if (op.eq.'n') call zeroit(xmat3,id1*id3)
        call dynamic_getmem(id1*id2,imat1)
        call dynamic_getmem(id2*id3,imat2)
        call dynamic_getmem(id1*id3,imat3)
        do irep=1,nrep
          irep_1=mtable(irep,ir1)    ! rep of rows of matrix 1
          irep_3=mtable(irep,ir2)    ! rep of rows of matrix 1
          isymdim1=vtable(1,irep_1)  ! row dim of mat 1
          isymdim2=vtable(1,irep)    ! col & row dim of mat 1 & 2
          isymdim3=vtable(1,irep_3)  ! col dim of mat 2
          if (isymdim1.eq.0.or.isymdim2.eq.0.or.isymdim3.eq.0) cycle
c
c                  Extract submatrix of A
c
         do k=1,isymdim2
           kb=vtable(k+1,irep)
           itmp_index=imat1+(k-1)*isymdim1
           do i=1,isymdim1
             ib=vtable(i+1,irep_1)
             bl(itmp_index+i-1)=xmat1(ib,kb)
           enddo
         enddo
c
c                  Extract submatrix of B
c
         do l=1,isymdim3
           lb=vtable(l+1,irep_3)
           itmp_index=imat2+(l-1)*isymdim2
           do k=1,isymdim2
             kb=vtable(k+1,irep)
             bl(itmp_index+k-1)=xmat2(kb,lb)
           enddo
         enddo
c
        beta=0d0
        call dgemm('N',      ! Do not transpose xmat1 ('T' = transpose)
     *             'N',      ! Do not transpose xmat2
     *             isymdim1, ! Rows of xmat1 and xmat3
     *             isymdim3, ! Columns of xmat2 and xmat3
     *             isymdim2, ! Columns of xmat1 & rows of xmat2
     *             1.0d0,    ! alpha of: C = alpha*op( A )*op( B ) + beta*C
     *             bl(imat1),! A
     *             isymdim1, ! Leading dimension of A
     *             bl(imat2),! B
     *             isymdim2, ! Leading dimension of B
     *             beta,     ! beta of: C = alpha*op( A )*op( B ) + beta*C
     *             bl(imat3),! C
     *             isymdim1) ! Leading dimension of C
c
c                  Distribute result over C:
c
         if (op.eq.'a') then
           do l=1,isymdim3
             lb=vtable(l+1,irep_3)
             itmp_index=imat3+(l-1)*isymdim1
             do i=1,isymdim1
               ib=vtable(i+1,irep_1)
               xmat3(ib,lb)=xmat3(ib,lb)+bl(itmp_index+i-1)
             enddo
           enddo
         else
           do l=1,isymdim3
             lb=vtable(l+1,irep_3)
             itmp_index=imat3+(l-1)*isymdim1
             do i=1,isymdim1
               ib=vtable(i+1,irep_1)
               xmat3(ib,lb)=bl(itmp_index+i-1)
             enddo
           enddo
         endif
c
        enddo
        call dynamic_retmem(3)
      endif
      end
C===============================================================================

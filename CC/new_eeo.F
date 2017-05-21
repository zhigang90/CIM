#include "maxi.h"
C==============================================================================
      subroutine NEW_EEO_INT(ncs,    ncf,     inx,   nval,  nmo,
     *                       nfirst, nlast,   thresh,
     *                       vorb,   ndisktr, npairs, af,iprnt,nbf)
      use memory
      use kinds
      use ccounters
      implicit none
c
      integer ncs,ncf,nval,nmo,nfirst,nlast,iprnt,nbf
      integer ndisktr,npairs,inx(12,*)
      logical vorb,af
      real*8  thresh
c
      integer, parameter :: igsize=15,jlchunkmax=400,ikchunkmax=5000
      integer jlarray(4,ncs*ncs),ikarray(4,ncs*ncs)
c     integer ish_tab_ik(2,ncs*(ncs+1)/2+1),! one more, must be terminated by 0
c    *        ish_tab_jl(2,ncs*(ncs+1)/2+1),! like "C" string ;)
c    *        iatpair(3,ncs*(ncs+1)/2+1)    ! This is because of symmetry,
c                                           ! the contents if this array
c                                           ! may be shorter
      integer ish_tab_ik,ish_tab_jl,iatpair
      real*8 timings(10),xinttime,t0,t1,xmax
      integer iout,igetival,ischwarz,ikchunk,ikchunkno,ndisk1,ndisk2
      integer ikstart,ikstop,iksize,ist,jlchunkno,ijchunk1,ijchunk2
      integer iskip,i,iscr_s,iscr_a,nskip,itot_int,irowskip,na
c
      call dynamic_mmark()
      call matmark()
      do i=1,10
        timings(i)=0d0
      enddo
c allocate memory for integer*2 arrays:
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_ik)
      call dynamic_getmem((ncs*(ncs+1)/2+1)/2    +1,ish_tab_jl)
c allocate memory for integer*4 array:
      call dynamic_getmem((ncs*(ncs+1)/2+1)*3/2  +1,iatpair)
c ish_tab_ik(1,ik) i shell of pair ik
c ish_tab_ik(1,ik) k shell of pair ik.
c pairs are randonly permutted, so this table is absolutely needed
c permutation allows equal load balance in parallel program
      call generate_ikpairs(bl(ish_tab_ik),ncs)
      call getival('na',na)
      call generate_jlpairs(bl(ish_tab_jl),ncs,inx,na,bl(iatpair),iprnt)
      iskip=0
      ic_NEW_EEO_INT=ic_NEW_EEO_INT+1
      ijchunk1=300 ! ij chunk of maps
      ijchunk2=300 ! ij chunk of results, as they are written and then read
      call  matsub('occu','cano',nfirst,nlast)
c
      iout=igetival('iout')
      if (iprnt.gt.1) then
      write(iout,*) '                            The EEOij Module  '
      write(iout,*) ' EEOij integral thresh    = ',thresh
      call flush(iout)
      endif
c
c Perform schwartz integrals calculation:
      call dynamic_getmem(ncs*ncs,ischwarz)
      call schwarz_calc(inx,ncs,bl(ischwarz))
c
c Space for screen matrices:
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_s)
      call dynamic_getmem(ncf*(ncf+1)/2,iscr_a)
      call zeroit(bl(iscr_s),ncf*(ncf+1)/2)
      call zeroit(bl(iscr_a),ncf*(ncf+1)/2)
c
c  Prepare array with jlchunk info
      call int_chunks_jl(inx,ncs,jlchunkmax,jlarray,jlchunkno,
     *                   bl(ish_tab_jl),bl(iatpair),iprnt)
      call int_chunks_ik(inx,ncs,ikchunkmax,ikarray,ikchunkno,
     *                   bl(ish_tab_ik))
c
      xmax=0d0
      call elapsec(t0)
      call prepare_ampl(inx,ijchunk1,ncf,nmo,ndisk1,af,vorb,
     *                  npairs,nval,ncs,jlchunkmax,jlarray,
     *                  jlchunkno,xmax,bl(iscr_s),bl(iscr_a),
     *                  bl(ish_tab_jl),nbf)
      call elapsec(t1)
      timings(1)=timings(1)+t1-t0
      if (iprnt.gt.2) write(6,*) 'The largest amplitude: ',xmax
c
c
      ndisk2=43
      call Open_Scratch_File(ndisk2,ikchunkmax,ijchunk2,af,'EEOb')
c    HERE CALCULATE ics, kcs
      nskip=0
c
      irowskip=0
      do ikchunk=1,ikchunkno
        ikstart=ikarray(1,ikchunk)
        ikstop =ikarray(2,ikchunk)
        iksize =ikarray(3,ikchunk)
c       ifirst =ikarray(4,ikchunk)
        call jl_loop(bl(ish_tab_ik),inx,  ncs,        thresh,  npairs,
     *               ischwarz, jlchunkno, jlchunkmax, jlarray, ijchunk1,
     *               ijchunk2, ikchunkno, ikchunkmax, ikstart, ikstop,
     *               iksize,   ikchunk,   ndisk1 ,    ndisk2 , af,
     *               iskip,    xinttime,  timings,bl(iscr_s),bl(iscr_a),
     *               nskip,    irowskip,  bl(ish_tab_jl))
      enddo
      itot_int=2*(ncf*(ncf+1)/2)**2
      if (iprnt.gt.2) then
      write(6,'(A,I8)')    "Total amount skipped: ",nskip
      write(6,'(A,F5.1)') "Percent skipped:      ",
     *                    dble(nskip)/itot_int*1d2
      write(6,'(A,F5.1)') "Row percent skipped:      ",
     *                    dble(irowskip)/itot_int*1d2
      endif
      if (af) then
        call fafClosem(ndisk1,0,ist)
      else
        close(ndisk1,STATUS='DELETE')
      endif
      if (ic_NEW_EEO_INT.gt.1) then 
        if (af) then
          call fafClosem(ndisktr,0,ist)
        else
          close(ndisktr,STATUS='DELETE')
        endif
      endif
      ndisktr=64
      call Open_Scratch_File(ndisktr,ncf,ncf,af,'EEOc')
c
c RESULTS
      call results_collector(ncf,npairs,af,ijchunk2,ikchunkno,
     *                       ikchunkmax,ndisk2,ndisktr,timings,
     *                       inx,bl(ish_tab_ik),ikarray)
c
      if (af) then
        call fafClosem(ndisk2,0,ist)
      else
        close(ndisk2,STATUS='DELETE')
      endif
      if (iprnt.gt.2) write(6,*) 'Shell pairs skipped: ',iskip
      call dynamic_retmark()
      call matremark()
 33   FORMAT (A,F10.3)
      if (iprnt.gt.1) then
      write(6,*) 'Amplitudes transformation: ', timings(1)
      write(6,*) 'Integral generation:       ', timings(2)
      write(6,*) 'Multiplication:            ', timings(3)
      write(6,*) 'Amplitudes reads:          ', timings(4)
      write(6,*) 'Results write:             ', timings(5)
      write(6,*) 'Results collection IO      ', timings(6)
      write(6,*) 'Screening:                 ', timings(7)
      endif
      end
C==============================================================================
      subroutine prepare_ampl(inx,ijchunkmax,ncf,nmo,ndisk,
     *                        af,vorb,npairs,nval,ncs,
     *                        jlchunkmax,jlarray,jlchunkno,xmax,
     *                        scr_s,scr_a,ish_tab,nbf)
      use memory
      implicit none
      integer inx(12,*),ijchunkmax,ncf,nmo,ndisk,npairs,nval,ncs
      integer jlchunkmax,jlarray(4,*),jlchunkno,nbf
      integer*2 ish_tab(2,*)
      logical af,vorb
      real*8 xmax,scr_s(*),scr_a(*)
c
      integer idimen,ijpasses,ijpass,iistart,iistop,ijsize
c
      call dynamic_mmark()
      call matmark()
      if (vorb) then
        idimen=nbf-nmo
      else
        idimen=ncf
      endif
      ndisk=78
      xmax=0d0
      call Open_Scratch_File(ndisk,jlchunkmax,ijchunkmax,af,'EEOa')
      ijpasses=npairs/ijchunkmax
      if (mod(npairs,ijchunkmax).ne.0) ijpasses=ijpasses+1
      do ijpass=1,ijpasses
        iistart= (ijpass-1)*ijchunkmax+1
        iistop = (ijpass  )*ijchunkmax
        if (iistart.gt.npairs) STOP '4rdfgtgbcs'
        if (iistop .gt.npairs) iistop=npairs
        ijsize=iistop-iistart+1
        call singl_ij(nval,  ncf,   ncs,       idimen,   inx,
     *                af,    vorb,  ndisk,     iistart,  iistop,
     *                ijpass,ijsize,ijchunkmax,jlchunkno,jlchunkmax,
     *                jlarray, xmax,scr_s,     scr_a,    ish_tab)
      enddo !  do ijpass=1,ijpasses
      call dynamic_retmark()
      call matremark()
      end
c===============================================================================
      subroutine ampl_sym_ant(ncf,workAA,sym,ant,xmax,scr_s,scr_a)
      implicit none
      integer ncf
      real*8 workAA(ncf,ncf),sym(*),ant(*),xmax,scr_s(*),scr_a(*)
c
      integer ij,i,j
      real*8 a,b
c
      ij=0
      do i=1,ncf
        do j=1,i
          ij=ij+1
          if (i.ne.j) then
            a=(workAA(i,j)+workAA(j,i))*0.5d0
            b=(workAA(i,j)-workAA(j,i))*0.5d0
            sym(ij)=a
            ant(ij)=b
            if (dabs(a).gt.xmax) xmax=dabs(a)
            if (dabs(b).gt.xmax) xmax=dabs(b)
            if (dabs(a).gt.scr_s(ij)) scr_s(ij)=dabs(a)
            if (dabs(b).gt.scr_a(ij)) scr_a(ij)=dabs(b)
          else
            a=workAA(i,j)
            sym(ij)=a
            ant(ij)=0d0
            if (dabs(a).gt.xmax) xmax=dabs(a)
            if (dabs(a).gt.scr_s(ij)) scr_s(ij)=dabs(a)
          endif
        enddo
      enddo
      end
c===============================================================================
      subroutine int_chunks_ik(inx,ncs,imaxchunk,iarray,ichunk,ipairs)
      implicit none
      integer inx(12,*),ncs,imaxchunk,iarray(4,*),ichunk
      integer*2 ipairs(2,*)
c
      integer ishpairs,ipart_chunk_size,ij,isize,jsize,i,j,ijsize,ijp
c
c iarray(4,chunk_number):
c 1: first ik pair of this chunk
c 2: last  ik pair
c 3: chunk size
c 4: first i index in this chunk (i member of the first ik pair)
c    - probably not used (I MUST CHECK IT)
c
      ishpairs=ncs*(ncs+1)/2
c
      ipart_chunk_size=0
      iarray(1,1)=1
      iarray(4,1)=1
      ichunk=0
      do ij=1,ishpairs+1
        ijp=ij
        i=ipairs(1,ij)   ! pairs are random! resulution array is used here
        j=ipairs(2,ij)
        if (i.eq.0) exit ! the last element, it may be shorter due to symmetry
                         ! the ipairs array ought to be ishpairs+1 long
        isize=inx(3,i)
        jsize=inx(3,j)
        if (i.ne.j) then
          ijsize=isize*jsize
        else
          ijsize=isize*(isize+1)/2
        endif
        if (ijsize.gt.imaxchunk) STOP 'Error xxfssrer'
        if ((ipart_chunk_size+ijsize).gt.imaxchunk) then
          ichunk=ichunk+1
          iarray(2,ichunk)=ij-1
          iarray(3,ichunk)=ipart_chunk_size
          ipart_chunk_size=ijsize
          iarray(1,ichunk+1)=ij
          iarray(4,ichunk+1)=i
        else
          ipart_chunk_size=ipart_chunk_size+ijsize
        endif
      enddo
      ichunk=ichunk+1
      iarray(2,ichunk)=ijp-1
      iarray(3,ichunk)=ipart_chunk_size
      end
c===============================================================================
      subroutine int_chunks_iktst(inx,ncs,imaxchunk,itail,ichunk,ipairs)
      implicit none
      integer inx(12,*),ncs,imaxchunk,itail,ichunk
      integer*2 ipairs(2,*)
c
      integer ishpairs,ipart_chunk_size,ij,isize,jsize,i,j,ijsize,ijp
c
c
      ishpairs=ncs*(ncs+1)/2
c
      ipart_chunk_size=0
      ichunk=0
      do ij=1,ishpairs
        ijp=ij
        i=ipairs(1,ij)
        j=ipairs(2,ij)
        if (i.eq.0) exit
        isize=inx(3,i)
        jsize=inx(3,j)
        if (i.ne.j) then
          ijsize=isize*jsize
        else
          ijsize=isize*(isize+1)/2
        endif
        if (ijsize.gt.imaxchunk) STOP 'Error fsstgreree'
        if ((ipart_chunk_size+ijsize).gt.imaxchunk) then
          ichunk=ichunk+1
          ipart_chunk_size=ijsize
        else
          ipart_chunk_size=ipart_chunk_size+ijsize
        endif
      enddo
      ichunk=ichunk+1
      itail=ipart_chunk_size
      end
c===============================================================================
      subroutine int_chunks_jl(inx,ncs,imaxchunk,iarray,ichunk,ipairs,
     *                         iatpair,iprnt)
      implicit none
      integer inx(12,*),ncs,imaxchunk,iarray(4,*),ichunk,iprnt
      integer*4 iatpair(3,*)
      integer*2 ipairs(2,*)
c
      integer ishpairs,ipart_chunk_size,ij,isize,jsize,i,j,ijsize
c
c iarray(4,chunk_number):
c 1: first jl pair of this chunk
c 2: last  jl pair
c 3: chunk size
c 4: first j index in this chunk (j component of the first jl pair)
c
      ishpairs=ncs*(ncs+1)/2
c
      ipart_chunk_size=0
      ij=0
      iarray(1,1)=1
      iarray(4,1)=1
      ichunk=0
      do ij=1,ishpairs
        i=ipairs(1,ij)
        j=ipairs(2,ij)
        isize=inx(3,i)
        jsize=inx(3,j)
        if (i.ne.j) then
          ijsize=isize*jsize
        else
          ijsize=isize*(isize+1)/2
        endif
        if (ijsize.gt.imaxchunk) STOP 'Error fssrer'
c       if ((ipart_chunk_size+ijsize).gt.imaxchunk .or.
c    *       ij.ne.1.and.iatpair(3,ij).ne.iatpair(3,ij-1).and.
c    *       (ipart_chunk_size.gt.10.or.iatpair(2,ij).gt.100)   ) then
        if ((ipart_chunk_size+ijsize).gt.imaxchunk             .or.
     *       ij.ne.1.and.iatpair(3,ij).ne.iatpair(3,ij-1)
     *       .and.(ipart_chunk_size.gt.60)                     .or.
     *       ij.ne.1.and.iatpair(1,ij).ne.iatpair(1,ij-1)
     *       .and.iatpair(2,ij)  .ge.100                       .or.
     *       ij.ne.1.and.iatpair(1,ij).ne.iatpair(1,ij-1)
     *       .and.iatpair(2,ij-1).ge.100             ) then
c       if ((ipart_chunk_size+ijsize).gt.imaxchunk             .or.
c    *       ij.ne.1.and.iatpair(3,ij).ne.iatpair(3,ij-1).and.
c    *       (ipart_chunk_size.gt.120)                         .or.
c    *       ij.ne.1.and.iatpair(1,ij).ne.iatpair(1,ij-1).and.
c    *       (ipart_chunk_size.gt.30.or.iatpair(2,ij).gt.100)   ) then
          ichunk=ichunk+1
          iarray(2,ichunk)=ij-1
          iarray(3,ichunk)=ipart_chunk_size
          ipart_chunk_size=ijsize
          iarray(1,ichunk+1)=ij
          iarray(4,ichunk+1)=i
        else
          ipart_chunk_size=ipart_chunk_size+ijsize
        endif
        if (ij.eq.ishpairs) then 
          ichunk=ichunk+1
          iarray(2,ichunk)=ij
          iarray(3,ichunk)=ipart_chunk_size
          exit
        endif
      enddo
      if (iprnt.gt.2) then
      do ij=1,ichunk
        write(6,'(5I8)') ij,iarray(1,ij),iarray(2,ij),iarray(3,ij),
     *                   iatpair(2,ij)
      enddo
      endif
      end
c===============================================================================
      subroutine int_chunks_jltst(inx,ncs,imaxchunk,itail,ichunk)
      implicit none
      integer inx(12,*),ncs,imaxchunk,itail,ichunk
c
      integer ishpairs,ipart_chunk_size,ij,isize,jsize,i,j,ijsize
c
c
      ishpairs=ncs*(ncs+1)/2
c
      ipart_chunk_size=0
      ij=0
      ichunk=0
      x:do i=1,ncs
          isize=inx(3,i)
          do j=1,i
            ij=ij+1
            jsize=inx(3,j)
            if (i.ne.j) then
              ijsize=isize*jsize
            else
              ijsize=isize*(isize+1)/2
            endif
            if (ijsize.gt.imaxchunk) STOP 'Error fsstgreree'
            if ((ipart_chunk_size+ijsize).gt.imaxchunk) then
              ichunk=ichunk+1
              ipart_chunk_size=ijsize
            else
              ipart_chunk_size=ipart_chunk_size+ijsize
            endif
            if (ij.eq.ishpairs) then 
              ichunk=ichunk+1
              itail=ipart_chunk_size
              exit x
            endif
          enddo
        enddo x
      end
c===============================================================================
c===============================================================================
      subroutine permute(iarray,n)
      implicit none
      integer n,iarray(n)
c
      integer iseed(10),in1,in2,ii,i
      real*8 x
c
      iseed(1)=34958475
      iseed(2)=34893759
      iseed(3)=6554984574
      iseed(4)=438374545
      call random_seed(SIZE=i)
      if (i.gt.4) STOP 'Increase seed'
      call random_seed(PUT=iseed)
      do i=1,n*50
        call random_number(x)
        in1=int(x*n)+1
        call random_number(x)
        in2=int(x*n)+1
        if (ii.eq.n+1) STOP 'Random number out of range'
c Swap one pair
        ii=iarray(in1)
        iarray(in1)=iarray(in2)
        iarray(in2)=ii
      enddo
      end
c===============================================================================
      subroutine permute_inx(inx1,inx,ncs)
      implicit none
      integer ncs,inx(12,ncs),inx1(12,ncs)
c
      integer list(ncs),i,j,inew
      do i=1,ncs
        list(i)=i
      enddo
      call permute(list,ncs)
      do i=1,ncs
        inew=list(i)
        do j=1,12
          inx(j,inew)=inx1(j,i)
        enddo
      enddo
      end
c===============================================================================
      subroutine itfer(ia,ib,n)
      integer n,ia(n),ib(n),i
      do i=1,n
        ib(i)=ia(i)
      enddo
      end
c===============================================================================
      subroutine generate_ikpairs(ikp,ncs)
      implicit none
      integer ncs
      integer*2 ikp(2,*)
c
      integer list(ncs*ncs),i,j,ij,ijp,ishpairs,ijj
      logical image
c
      ishpairs=ncs*(ncs+1)/2
      call i2zeroit(ikp,ishpairs*2+2)
      do i=1,ishpairs
        list(i)=i
      enddo
      call permute(list,ishpairs)
      ij=0
      do ijj=1,ishpairs
        ijp=list(ijj)
        call get_ij_half(ijp,i,j)
        if (image(i,j)) cycle
        ij=ij+1
        ikp(1,ij)=i
        ikp(2,ij)=j
      enddo
      end
C===============================================================================
      subroutine results_collector(ncf,npairs,af,ijchunk,ikchunkno,
     *                             ikchunkmax,ndisks,ndiskd,timings,
     *                             inx,ish_tab_ik,ikarray)
      use memory
      implicit none
      integer ncf,npairs,ijchunk,ikchunkno,ikchunkmax,ndisks,ndiskd
      integer inx(12,*),ish_tab_ik(2,*),ikarray(4,*)
      logical af
      real*8 timings(*)
c
      integer ijpasses,ijpass,ijstart,ijstop,ijsize
c
c
      ijpasses=npairs/ijchunk
      if (mod(npairs,ijchunk).ne.0) ijpasses=ijpasses+1
c
      do ijpass=1,ijpasses
        ijstart= (ijpass-1)*ijchunk+1
        ijstop = (ijpass  )*ijchunk
        if (ijstart.gt.npairs) STOP '4rdfgtgbcs'
        if (ijstop .gt.npairs) ijstop=npairs
        ijsize=ijstop-ijstart+1
        call ij_res(ncf, af, inx, ish_tab_ik, ikchunkno,
     *              ikchunkmax, ikarray, ndisks, ndiskd, ijsize,
     *              ijpass, ijstart, ijstop, ijchunk,  timings)
      enddo
      end
C===============================================================================
      subroutine jl_loop(ish_tab_ik,inx,    ncs,       thresh, npairs,
     *                   ischwarz,jlchunkno,jlchunkmax,jlarray,ijchunk1,
     *                   ijchunk2,ikchunkno,ikchunkmax,ikstart,ikstop,
     *                   iksize,  ikchunk,  ndisk_a,   ndisk_r,af,
     *                   iskip,   xinttime, timings,scr_s,scr_a,nskip,
     *                   irowskip,ish_tab_jl)
      use memory
      implicit none
      integer inx(12,*),ischwarz,jlchunkno,ikstart
      integer ikstop
      integer*2 ish_tab_jl(2,*),ish_tab_ik(2,*)
      integer iksize,ncs,iskip,ikchunkmax,jlchunkmax,ijchunk1,ijchunk2
      integer npairs,ndisk_a,ndisk_r,ikchunk,ikchunkno,jlarray(4,*)
      integer nskip,irowskip
      real*8 timings(*),thresh,xinttime,scr_s(*),scr_a(*)
      logical af
c
      integer jlchunk,ibuf1,ibuf2,ibuf3,isym,iant,ijpasses1,iisize1
      integer ijpasses2,iisize2,iparts,iparta,ibigress,ibigresa
      integer iress,iresa,irecs,ireca,ist,ijsize,ijstart,ijstop,ijpass
      integer jlstart,jlstop,jlsize,imax_res_memory,izer_ints,izer_inta
      integer*1 isparses(ikchunkmax),isparsea(ikchunkmax)
      integer i,imat_s,imat_a,mataddr,ioutiksizes,ioutiksizea
      logical do_sym_compr,do_ant_compr
      real*8 t0,t1,t2
      integer, parameter :: igsize=15
      logical screen
c
      call dynamic_mmark()
      call matmark()
c
      screen=.true.
      ijpasses1=npairs/ijchunk1
      if (mod(npairs,ijchunk1).ne.0) ijpasses1=ijpasses1+1
      iisize1=ijchunk1*jlchunkmax*8
      ijpasses2=npairs/ijchunk2
      if (mod(npairs,ijchunk2).ne.0) ijpasses2=ijpasses2+1
      iisize2=ijchunk2*ikchunkmax*8
      imax_res_memory=max(ijpasses1*ijchunk1*ikchunkmax,
     *                    ijpasses2*ijchunk2*ikchunkmax)
c
      call dynamic_getmem(igsize**4,ibuf1) ! will be always fixed size
      call dynamic_getmem(igsize**4,ibuf2) ! will be always fixed size
      call dynamic_getmem(igsize**4,ibuf3) ! will be always fixed size
      call dynamic_getmem(jlchunkmax*ikchunkmax,isym) ! integral chunk symm
      call dynamic_getmem(jlchunkmax*ikchunkmax,iant) ! integral chunk asym.
      if (screen) then
        call dynamic_getmem(ijchunk1*ikchunkmax,imat_s) ! compressed result
        call dynamic_getmem(ijchunk1*ikchunkmax,imat_a) ! compressed result
      endif
      call dynamic_getmem(ijchunk1*jlchunkmax,iparts) !for symm. ampl chunk
      call dynamic_getmem(ijchunk1*jlchunkmax,iparta) !for asym. ampl chunk
      call dynamic_getmem(imax_res_memory,ibigress)
      call dynamic_getmem(imax_res_memory,ibigresa)
      call zeroit(bl(ibigress),imax_res_memory)
      call zeroit(bl(ibigresa),imax_res_memory)
c
      do jlchunk=1,jlchunkno
        jlstart=jlarray(1,jlchunk)
        jlstop =jlarray(2,jlchunk)
        jlsize =jlarray(3,jlchunk)
        do i=1,ikchunkmax
          isparses(i)=0
        enddo
        do i=1,ikchunkmax
          isparsea(i)=0
        enddo
        call elapsec(t0)
        call integral_chunk(ish_tab_ik,inx,bl(isym),bl(iant),bl(ibuf1),
     *                      bl(ibuf2),bl(ibuf3),bl(ischwarz),bl,
     *                      ikstart,ikstop,iksize,jlstart,jlstop,
     *                      jlsize,ncs,iskip,thresh,xinttime,
     *                      scr_s,scr_a,nskip,isparses,isparsea,
     *                      ish_tab_jl)
        call elapsec(t1)
        timings(2)=timings(2)+t1-t0
        izer_ints=0
        izer_inta=0
        do i=1,iksize
          if (isparses(i).eq.0) izer_ints=izer_ints+1
        enddo
        do i=1,iksize
          if (isparsea(i).eq.0) izer_inta=izer_inta+1
        enddo
        do_sym_compr=.false.
        do_ant_compr=.false.
        ioutiksizes=iksize
        ioutiksizea=iksize
        if (screen.and.dble(izer_ints)/iksize.gt.0.2d0) then
          irowskip=irowskip+jlsize*izer_ints
          do_sym_compr=.true.
          call compress_mat(bl(isym),isparses,jlsize,iksize,ioutiksizes)
          if (iksize-izer_ints.ne.ioutiksizes) then
            print *,iksize,izer_ints,ioutiksizes
            STOP 'sdfsdf4433s'
          endif
        endif
        if (screen.and.dble(izer_inta)/iksize.gt.0.2d0) then
          irowskip=irowskip+jlsize*izer_inta
          do_ant_compr=.true.
          call compress_mat(bl(iant),isparsea,jlsize,iksize,ioutiksizea)
          if (iksize-izer_inta.ne.ioutiksizea) then
            print *,iksize,izer_inta,ioutiksizea
            STOP 'sdfsdf4433a'
          endif
        endif
c
c       if (ioutiksizes.ne.0)
               call matconn('intsym','r',jlsize,ioutiksizes,isym)
c       if (ioutiksizea.ne.0)
               call matconn('intant','r',jlsize,ioutiksizea,iant)
        do ijpass=1,ijpasses1
          ijstart= (ijpass-1)*ijchunk1+1
          ijstop = (ijpass  )*ijchunk1
          if (ijstart.gt.npairs) STOP '4rdfgtgbcs'
          if (ijstop .gt.npairs) ijstop=npairs
          ijsize=ijstop-ijstart+1
          iress=ibigress+(ijstart-1)*iksize
          iresa=ibigresa+(ijstart-1)*iksize
          if (do_sym_compr) then
            if (ioutiksizes.ne.0)
     *           call matconn('ressym','r',ioutiksizes,ijsize,imat_s)
          else
            call matconn('ressym','r',iksize,ijsize,iress)
          endif
          if (do_ant_compr) then
            if (ioutiksizea.ne.0)
     *           call matconn('resant','r',ioutiksizea,ijsize,imat_a)
          else
            call matconn('resant','r',iksize,ijsize,iresa)
          endif
          irecs=((ijpass-1)*jlchunkno+jlchunk)*2-1
          ireca=((ijpass-1)*jlchunkno+jlchunk)*2
          call elapsec(t0)
          if (af) then
            call fafread(ndisk_a,bl(iparts),8,jlsize*ijsize,1,irecs,ist)
            if (ist.lt.8*jlsize*ijsize) then
              ERROR
              call flush(6)
            endif
            call fafread(ndisk_a,bl(iparta),8,jlsize*ijsize,1,ireca,ist)
            if (ist.lt.8*jlsize*ijsize) then
              ERROR
              call flush(6)
            endif
          else
            call general_read(ndisk_a,irecs,bl(iparts),iisize1)
            call general_read(ndisk_a,ireca,bl(iparta),iisize1)
          endif
          call elapsec(t1)
          timings(4)=timings(4)+t1-t0
          call matconn('amsym','r',jlsize,ijsize,iparts)
          call matconn('amant','r',jlsize,ijsize,iparta)
          if (do_sym_compr) then
            if (ioutiksizes.ne.0) then
              call elapsec(t0)
              call matmmul2('intsym','amsym','ressym','t','n','n')
              call elapsec(t1)
              timings(3)=timings(3)+t1-t0
              call add_result(bl(imat_s),bl(iress),isparses,ioutiksizes,
     *                        iksize,    ijsize)
              call elapsec(t2)
              timings(7)=timings(7)+t2-t1
            endif
          else
            call elapsec(t0)
            call matmmul2('intsym','amsym','ressym','t','n','a')
            call elapsec(t1)
            timings(3)=timings(3)+t1-t0
          endif
          if (do_ant_compr) then
            if (ioutiksizea.ne.0) then
              call elapsec(t0)
              call matmmul2('intant','amant','resant','t','n','n')
              call elapsec(t1)
              timings(3)=timings(3)+t1-t0
              call add_result(bl(imat_a),bl(iresa),isparsea,ioutiksizea,
     *                        iksize,    ijsize)
              call elapsec(t2)
              timings(7)=timings(7)+t2-t1
            endif
          else
            call elapsec(t0)
            call matmmul2('intant','amant','resant','t','n','a')
            call elapsec(t1)
            timings(3)=timings(3)+t1-t0
          endif
          call matdisc('amant')
          call matdisc('amsym')
          if (do_ant_compr) then
            if (ioutiksizea.ne.0) call matdisc('resant')
          else
            call matdisc('resant')
          endif
          if (do_sym_compr) then
            if (ioutiksizes.ne.0) call matdisc('ressym')
          else
            call matdisc('ressym')
          endif
        enddo
        call matdisc('intant')
        call matdisc('intsym')
      enddo
c
      do ijpass=1,ijpasses2
        ijstart= (ijpass-1)*ijchunk2+1
        ijstop = (ijpass  )*ijchunk2
        if (ijstart.gt.npairs) STOP '4rdfgtgbcs'
        if (ijstop .gt.npairs) ijstop=npairs
        ijsize=ijstop-ijstart+1
        iress=ibigress+(ijstart-1)*iksize
        iresa=ibigresa+(ijstart-1)*iksize
        irecs=((ijpass-1)*ikchunkno+ikchunk)*2-1
        ireca=((ijpass-1)*ikchunkno+ikchunk)*2
        call elapsec(t0)
        if (af) then
         call fafwrite(ndisk_r,bl(iress),8,iksize*ijsize,1,irecs,ist)
         call fafwrite(ndisk_r,bl(iresa),8,iksize*ijsize,1,ireca,ist)
        else
          call general_write(ndisk_r,irecs,bl(iress),iisize2)
          call general_write(ndisk_r,ireca,bl(iresa),iisize2)
        endif
        call elapsec(t1)
        timings(5)=timings(5)+t1-t0
      enddo
      call dynamic_retmark()
      call matremark()
      end
C===============================================================================
      subroutine singl_ij(nval,  ncf,   ncs,       idimen,   inx,
     *                    af,    vorb,  ndisk,     iistart,  iistop,
     *                    ijpass,ijsize,ijchunkmax,jlchunkno,jlchunkmax,
     *                    jlarray, xmax,scr_s,     scr_a,    ish_tab)
      use memory
      implicit none
      integer nval,iistart,iistop,ijsize,idimen,ncf,jlchunkno,ijchunkmax
      integer jlchunkmax,jlarray(4,*),ijpass,ndisk,ncs,inx(12,*)
      integer*2 ish_tab(2,*)
      logical af,vorb
      real*8 xmax,scr_s(*),scr_a(*)
c
      integer ij,i,j,ino,icoeff,isym,iant,iparts,iparta,iworkAA
      integer isingleMO,isingleAO,mataddr,is,ia,jlchunk,irecs,ireca
      integer ikd,iss,iaa,istart,istop,kstart,kstop,ii,kk,ist,iks
      integer ics,kcs,ikshno,iksize,isize,ntrincf
c
c
      call dynamic_mmark()
      call matmark()
      isize=ijchunkmax*jlchunkmax*8
      ntrincf =ncf*(ncf+1)/2
      call dynamic_matdef('workAA','q',ncf,ncf)    ! atomic x atomic
      call dynamic_matdef('workAM','r',ncf,idimen) ! atomic x molecular
      iworkAA=mataddr('workAA')
      call dynamic_matdef('singleMO','r',idimen,1)
      isingleMO=mataddr('singleMO')
      call dynamic_matdef('singleAO','r',ncf,1)
      isingleAO=mataddr('singleAO')
      call dynamic_getmem(ijchunkmax*ntrincf,isym)
      call dynamic_getmem(ijchunkmax*ntrincf,iant)
      call dynamic_getmem(ijchunkmax*jlchunkmax,iparts)
      call dynamic_getmem(ijchunkmax*jlchunkmax,iparta)
c
      ij=0
      do i=1,nval
        do j=1,i
          ij=ij+1
          if (ij.lt.iistart) cycle
          if (ij.gt.iistop ) exit
          ino=ij-iistart
          call CoefRead('tc',i,j,icoeff)
          call matconn('coeff','q',idimen,idimen,icoeff)
          if (vorb) then
            call matmmult('genvirt','coeff','workAM')
            call matmmul2('workAM','genvirt','workAA','n','t','n')
          else
            call matcopy('coeff','workAA')
          endif
          call matdisc('coeff')
c Singles part of EEO:
          if (vorb) then
          call readsingles(j,bl(isingleMO))
          call matmmult('genvirt','singleMO','singleAO')
          else
          call readsingles(j,bl(isingleAO))
          endif
          call matsub('singlocc','occu',i,i)
          call matmmul2('singlocc','singleAO','workAA','n','t','a')
          call matdisc('singlocc')
c
          if (vorb) then
          call readsingles(i,bl(isingleMO))
          call matmmult('genvirt','singleMO','singleAO')
          else
          call readsingles(i,bl(isingleAO))
          endif
          call matsub('singlocc','occu',j,j)
          call matmmul2('singleAO','singlocc','workAA','n','t','a')
          call matdisc('singlocc')
c Singles part of EEO end
          is=isym+ntrincf*ino
          ia=iant+ntrincf*ino
          call ampl_sym_ant(ncf,bl(iworkAA),bl(is),bl(ia),xmax,
     *                      scr_s,scr_a)
        enddo
      enddo
c Split into shell chunks
      do jlchunk=1,jlchunkno
        call zeroit(bl(iparts),ijchunkmax*jlchunkmax) ! maybe remove? check it
        call zeroit(bl(iparta),ijchunkmax*jlchunkmax)
        iksize=jlarray(3,jlchunk)
        ikd=-1
        do ij=iistart,iistop
          iss=isym+(ij-iistart)*ntrincf
          iaa=iant+(ij-iistart)*ntrincf
c         ikshno=0
          do ikshno=jlarray(1,jlchunk),jlarray(2,jlchunk)
              ics=ish_tab(1,ikshno)
              kcs=ish_tab(2,ikshno)
              istart=inx(11,ics)+1
              istop =inx(10,ics)
              kstart=inx(11,kcs)+1
              kstop =inx(10,kcs)
c         do ics=1,ncs
c           istart=inx(11,ics)+1
c           istop =inx(10,ics)
c           do kcs=1,ics
c             ikshno=ikshno+1
c             if (ikshno.lt.jlarray(1,jlchunk)) cycle
c             if (ikshno.gt.jlarray(2,jlchunk)) exit
c             kstart=inx(11,kcs)+1
c             kstop =inx(10,kcs)
              do ii=istart,istop
                do kk=kstart,kstop
                  if (ics.eq.kcs.and.kk.gt.ii) exit ! for ics.ne.kcs kk can 
                                                    ! be .lt. ii, if inx
                                                    ! is rearranged
                  if (ii.ge.kk) then          ! see above
                    iks=ii*(ii-1)/2+kk-1
                    ikd=ikd+1
                    bl(iparts+ikd)=bl(iss+iks)
                    bl(iparta+ikd)=bl(iaa+iks)
                  else
                    iks=kk*(kk-1)/2+ii-1
                    ikd=ikd+1
                    bl(iparts+ikd)=bl(iss+iks)
                    bl(iparta+ikd)=-bl(iaa+iks)
                  endif
                enddo
              enddo
c           enddo
          enddo
        enddo
        irecs=((ijpass-1)*jlchunkno+jlchunk)*2-1
        ireca=((ijpass-1)*jlchunkno+jlchunk)*2
        if (af) then
          call fafwrite(ndisk,bl(iparts),8,ijsize*iksize,1,irecs,ist)
          call fafwrite(ndisk,bl(iparta),8,ijsize*iksize,1,ireca,ist)
        else
          call general_write(ndisk,irecs,bl(iparts),isize)
          call general_write(ndisk,ireca,bl(iparta),isize)
        endif
      enddo ! over chunks
      call dynamic_retmark()
      call matremark()
      end
C===============================================================================
      subroutine ij_res(ncf, af, inx, ish_tab_ik, ikchunkno,
     *                  ikchunkmax, ikarray, ndisks, ndiskd, ijsize,
     *                  ijpass, ijstart, ijstop, ijchunk,  timings)
      use memory
      implicit none
      integer ncf,inx(12,*),ikchunkno,ikchunkmax
      integer*2 ish_tab_ik(2,*)
      integer ikarray(4,*),ndisks,ndiskd,ijsize,ijchunk,ijstart,ijpass
      integer ijstop
      logical af
      real*8 timings(8)
c
      integer ikchunk,ikstart,ikstop,iksize,irecs,ireca,ist,iikk,ij
      integer iorb,jorb,ionemat,iparts,iparta,iresult,ics,kcs,istart
      integer istop,kstart,kstop,ii,kk,list1(3,10),iip,kkp,iii,ikcc
      integer iisize
      real*8 t0,t1,el_ik,el_ki
c
      call dynamic_mmark()
      call matmark()
c
      iisize=ijchunk*ikchunkmax*8
      call dynamic_getmem(ijchunk*ikchunkmax,iparts) !for res. chunk
      call dynamic_getmem(ijchunk*ikchunkmax,iparta) !for res. chunk
      call dynamic_getmem(ijchunk*ncf*ncf,iresult)
c
      do ikchunk=1,ikchunkno
        ikstart=ikarray(1,ikchunk)
        ikstop =ikarray(2,ikchunk)
        iksize =ikarray(3,ikchunk)
        irecs=((ijpass-1)*ikchunkno+ikchunk)*2-1
        ireca=((ijpass-1)*ikchunkno+ikchunk)*2
        call elapsec(t0)
        if (af) then
          call fafread(ndisks,bl(iparts),8,iksize*ijsize,1,irecs,ist)
          if (ist.lt.8*iksize*ijsize) then
            ERROR
            call flush(6)
          endif
          call fafread(ndisks,bl(iparta),8,iksize*ijsize,1,ireca,ist)
          if (ist.lt.8*iksize*ijsize) then
            ERROR
            call flush(6)
          endif
        else
          call general_read(ndisks,irecs,bl(iparts),iisize)
          call general_read(ndisks,ireca,bl(iparta),iisize)
        endif
        call elapsec(t1)
        timings(6)=timings(6)+t1-t0
        iikk=0
        do ij=ijstart,ijstop
          call get_ij_half(ij,iorb,jorb)
          ionemat=iresult+(ij-ijstart)*ncf*ncf
          do ikcc=ikstart,ikstop
            ics=ish_tab_ik(1,ikcc)
            kcs=ish_tab_ik(2,ikcc)
            istart=inx(11,ics)+1
            istop =inx(10,ics)
            kstart=inx(11,kcs)+1
            kstop =inx(10,kcs)
            do ii=istart,istop
              do kk=kstart,kstop
                if (kk.gt.ii.and.ics.eq.kcs) exit
                iikk=iikk+1
                el_ik=bl(iparts+iikk-1)+bl(iparta+iikk-1)
                el_ki=bl(iparts+iikk-1)-bl(iparta+iikk-1)
                bl(ionemat+ii-1+(kk-1)*ncf)=el_ik
                bl(ionemat+kk-1+(ii-1)*ncf)=el_ki
                call images_list(iorb,jorb,ii,kk,list1,.false.)
                do iii=1,10
                  if (list1(3,iii).eq.0) exit
                  iip=list1(1,iii)
                  kkp=list1(2,iii)
                  bl(ionemat+(iip-1)+(kkp-1)*ncf)=el_ik*list1(3,iii)
                enddo
                call images_list(iorb,jorb,kk,ii,list1,.false.)
                do iii=1,10
                  if (list1(3,iii).eq.0) exit
                  iip=list1(1,iii)
                  kkp=list1(2,iii)
                  bl(ionemat+(iip-1)+(kkp-1)*ncf)=el_ki*list1(3,iii)
                enddo
              enddo
            enddo
c
          enddo

        enddo
      enddo
      call elapsec(t0)
      do ij=ijstart,ijstop
        ionemat=iresult+(ij-ijstart)*ncf*ncf
        if (af) then
          call fafwrite(ndiskd,bl(ionemat),8,ncf*ncf,1,ij,ist)
        else
          call general_write(ndiskd,ij,bl(ionemat),ncf*ncf*8)
        endif
      enddo
      call elapsec(t1)
      timings(6)=timings(6)+t1-t0
c
      call dynamic_retmark()
      call matremark()
      end
c===============================================================================
      subroutine compress_mat(xmat,icompr_info,irow,icol,ioutcol)
      implicit none
      integer irow,icol,ioutcol
      integer*1 icompr_info(*)
      real*8 xmat(irow,icol)
c
      integer itar,i
c
      itar=0
      do i=1,icol
        if (icompr_info(i).eq.1) then
          itar=itar+1
          if (i.gt.itar) call tfer(xmat(1,i),xmat(1,itar),irow)
        endif
      enddo
      ioutcol=itar
      end
c===============================================================================
      subroutine xdd_result(compr,uncompr,icomp_data,icomprow,irow,icol)
      implicit none
      integer icomprow,irow,icol
      integer*1 icomp_data(*)
      real*8 compr(icomprow,icol),uncompr(irow,icol)
c
      integer idst,isrc,j
c
      isrc=0
      do idst=1,irow
        if (icomp_data(idst).eq.1) then
          isrc=isrc+1
          do j=1,icol
            uncompr(idst,j)=uncompr(idst,j)+compr(isrc,j)
          enddo
        endif
      enddo
      if (isrc.ne.icomprow) STOP '3434zsdkfdj35yfhw'
      end
c===============================================================================
      subroutine add_result(compr,uncompr,icomp_data,icomprow,irow,icol)
      implicit none
      integer icomprow,irow,icol
      integer*1 icomp_data(*)
      real*8 compr(icomprow,icol),uncompr(irow,icol)
c
      integer idst,isrc,j
c
      do j=1,icol
        isrc=0
        do idst=1,irow
          if (icomp_data(idst).eq.1) then
            isrc=isrc+1
            uncompr(idst,j)=uncompr(idst,j)+compr(isrc,j)
          endif
        enddo
      enddo
      if (isrc.ne.icomprow) STOP '3434zsdkfdj35yfhw'
      end
c===============================================================================
      subroutine xenerate_jlpairs(ikp,ncs)
      implicit none
      integer ncs
      integer*2 ikp(2,*)
c
      integer list(ncs*ncs),i,j,ij,ijp,ishpairs
      logical image
c
      ishpairs=ncs*(ncs+1)/2
      call i2zeroit(ikp,ishpairs*2+2)
      do i=1,ishpairs
        list(i)=i
      enddo
      call permute(list,ishpairs)
      do ij=1,ishpairs
        ijp=list(ij)
        call get_ij_half(ijp,i,j)
        ikp(1,ij)=i
        ikp(2,ij)=j
        if (mod(ij,3).eq.0) then
          ikp(1,ij)=j
          ikp(2,ij)=i
        endif
      enddo
      end
c===============================================================================
      subroutine generate_jlpairs(ikp,ncs,inx,na,iato,iprnt)
      use memory
      implicit none
      integer ncs,inx(12,*),na,iprnt
      integer*2 ikp(2,*)
      integer*4 iato(3,*)
c
      integer ij,ishpairs,iat,jat,ish,jsh,iatpair,ibas,idif,jdif,i,ishp
      integer ia,ja,isize,jsize,ijsize,istored,inuc,j
      real*8  l1,l2,r1,r2,expon,expon1,expon2
      integer*1 aa(ncs*(ncs+1)/2)
      integer grouplist(10,na),nogroup,iiat,jjat,igr,jgr,igrpair
c
      ishpairs=ncs*(ncs+1)/2
      iatpair=0
      igrpair=0
      call getival('ibas',ibas)
      call i1zeroit(aa,ishpairs)
      call getival('inuc',inuc)
      call atom_groups(na,bl(inuc),grouplist,nogroup)
      if (iprnt.gt.2) then
      do i=1,na
        if (grouplist(1,i).eq.0) exit
        do j=1,10
          if (grouplist(j,i).eq.0) exit
          write(6,'(I3)',ADVANCE='NO') grouplist(j,i)
        enddo
        write(6,*)
      enddo
      endif
c
c
      ij=0
c     do ish=1,ncs
c       do jsh=1,ish
c         isize=inx(3,ish)
c         jsize=inx(3,jsh)
c         if (ish.ne.jsh) then
c           ijsize=isize*jsize
c         else
c           ijsize=isize*(isize+1)/2
c         endif
c         if (ijsize.eq.1) then
c           ishp=ish*(ish-1)/2+jsh
c           ij=ij+1
c           ikp(1,ij)=ish
c           ikp(2,ij)=jsh
c           iato(1,ij)=iatpair
c           if (aa(ishp).eq.1) STOP 'wawa'
c           aa(ishp)=1
c         endif
c       enddo
c     enddo
c
c
      do idif=1,2
        if (idif.eq.1) then
          l1=0.5d0
          l2=1d20
        else
          l1=0d0
          l2=0.5d0
        endif
        do jdif=1,idif
          if (jdif.eq.1) then
            r1=0.5d0
            r2=1d20
          else
            r1=0d0
            r2=0.5d0
          endif
c         iatpair=0
          do igr=1,nogroup
            do jgr=1,nogroup
              if (idif.eq.jdif.and.jgr.gt.igr) cycle
              igrpair=igrpair+1
              do iiat=1,grouplist(1,igr)
                iat=grouplist(iiat+1,igr)
                if (iat.eq.0) STOP 'wer4tzffg'
                do jjat=1,grouplist(1,jgr)
                  jat=grouplist(jjat+1,jgr)
                  if (jat.eq.0) STOP 'wer54jxgl4tzffg'
                  if (idif.eq.jdif.and.igr.eq.jgr.and.jat.gt.iat) cycle
                  iatpair=iatpair+1
                  do ish=1,ncs
                    call lowest_exponent(inx(1,ish),bl(ibas),expon)
                    if (expon.le.l1.or.expon.gt.l2) cycle
                    if (inx(2,ish).ne.iat) cycle
                    do jsh=1,ncs
                      call lowest_exponent(inx(1,jsh),bl(ibas),expon)
                      if (expon.le.r1.or.expon.gt.r2) cycle
                      if (inx(2,jsh).ne.jat) cycle
                      if (idif.eq.jdif.and.igr.eq.jgr.and.iat.eq.jat
     *                                            .and.jsh.gt.ish) cycle
                      isize=inx(3,ish)
                      jsize=inx(3,jsh)
                      if (ish.ne.jsh) then
                        ijsize=isize*jsize
                      else
                        ijsize=isize*(isize+1)/2
                      endif
c                     if (ijsize.eq.1) cycle
                      if (ish.gt.jsh) then
                        ishp=ish*(ish-1)/2+jsh
                      else
                        ishp=jsh*(jsh-1)/2+ish
                      endif
                      ij=ij+1
                      ikp(1,ij)=ish
                      ikp(2,ij)=jsh
                      iato(1,ij)=iatpair
                      iato(3,ij)=igrpair
                      if (aa(ishp).eq.1) STOP 'wawa'
                      aa(ishp)=1
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
      enddo
      if (ij.ne.ishpairs) STOP 'Llesdf9f88u4yh3terk'
      istored=0
      do ij=1,ishpairs
        ish=ikp(1,ij)
        jsh=ikp(2,ij)
        isize=inx(3,ish)
        jsize=inx(3,jsh)
        if (ish.ne.jsh) then
          ijsize=isize*jsize
        else
          ijsize=isize*(isize+1)/2
        endif
        iato(2,ij)=0
        if (ij.ne.1.and.iato(1,ij).ne.iato(1,ij-1)) then
          iato(2,ij-1)=istored
          istored=0
        endif
        istored=istored+ijsize
      enddo
      iato(2,ishpairs)=istored
      do ij=ishpairs,2,-1
        if (iato(2,ij-1).eq.0)  iato(2,ij-1)=iato(2,ij)
      enddo
      end
C===============================================================================
      subroutine lowest_exponent(inx,basdat,expon)
      implicit none
      integer inx(12)
      real*8 basdat(13,*),expon
c
      integer i
c
      expon=1d20
      do i=inx(1)+1,inx(5)
        if (basdat(1,i).lt.expon) expon=basdat(1,i)
      enddo
      end
C===============================================================================
      subroutine atom_groups(na,natoms,grouplist,ngr)
      implicit none
      integer na,grouplist(10,na),ngr
      real*8 natoms(5,*)
c
      integer inuc,i,nx,ny,nz,iat,ij,nat
      real*8 x1,x2,y1,y2,z1,z2,xa,xb,ya,yb,za,zb,x_start,x_stop,x,y,z
      real*8 y_start,y_stop,z_start,z_stop,tmp,signx,signy,xm,ym,zm
      real*8, parameter :: box_size=5d0
      logical added
c
      call izeroit(grouplist,10*na)
      x1=+0d5
      x2=-0d5
      y1=+0d5
      y2=-0d5
      z1=+0d5
      z2=-0d5
c Calculate lowest and biggest x,y,z for the molecule:
      do i=1,na
        if (natoms(2,i).lt.x1) x1=natoms(2,i)
        if (natoms(2,i).gt.x2) x2=natoms(2,i)
        if (natoms(3,i).lt.y1) y1=natoms(3,i)
        if (natoms(3,i).gt.y2) y2=natoms(3,i)
        if (natoms(4,i).lt.z1) z1=natoms(4,i)
        if (natoms(4,i).gt.z2) z2=natoms(4,i)
      enddo
      x1=x1-1d-3
      x2=x2+1d-3
      y1=y1-1d-3
      y2=y2+1d-3
      z1=z1-1d-3
      z2=z2+1d-3
      nx=ceiling((x2-x1)/box_size)
      ny=ceiling((y2-y1)/box_size)
      nz=ceiling((z2-z1)/box_size)
      xm=(x1+x2)/2d0
      ym=(y1+y2)/2d0
      zm=(z1+z2)/2d0
      if (mod(nx,2).eq.0) then
        x1=xm-nx/2*box_size
        x2=xm+nx/2*box_size
      else
        x1=xm-(nx-1)/2*box_size-0.5d0*box_size
        x2=xm+(nx-1)/2*box_size+0.5d0*box_size
      endif
      if (mod(ny,2).eq.0) then
        y1=ym-ny/2*box_size
        y2=ym+ny/2*box_size
      else
        y1=ym-(ny-1)/2*box_size-0.5d0*box_size
        y2=ym+(ny-1)/2*box_size+0.5d0*box_size
      endif
      if (mod(nz,2).eq.0) then
        z1=zm-nz/2*box_size
        z2=zm+nz/2*box_size
      else
        z1=zm-(nz-1)/2*box_size-0.5d0*box_size
        z2=zm+(nz-1)/2*box_size+0.5d0*box_size
      endif
      zb=z1
      yb=y1
      xb=x1
      signy=1d0
      signx=1d0
      ij=1
      do          ! over z
        za=zb
        zb=za+box_size
        if (za.gt.z2) exit
        z_start=za
        z_stop =zb
        do          ! over y
          ya=yb
          yb=ya+box_size*signy
          if (signy*ya.gt.signy*y2) exit
          if (ya.lt.yb) then
            y_start=ya
            y_stop =yb
          else
            y_start=yb
            y_stop =ya
          endif
          do          ! over x
            xa=xb
            xb=xa+box_size*signx
            if (signx*xa.gt.signx*x2) exit
            if (xa.lt.xb) then
              x_start=xa
              x_stop =xb
            else
              x_start=xb
              x_stop =xa
            endif
            added=.false.
            do iat=1,na
              x=natoms(2,iat)
              y=natoms(3,iat)
              z=natoms(4,iat)
              if (  x_start.le.x .and. x.lt.x_stop  .and.
     *              y_start.le.y .and. y.lt.y_stop  .and.
     *              z_start.le.z .and. z.lt.z_stop         ) then
                nat=grouplist(1,ij)
                grouplist(1,ij)=nat+1
                grouplist(nat+2,ij)=iat
                added=.true.
              endif
            enddo
            if (added) ij=ij+1
          enddo
          signx=-signx
          tmp=x2
          x2=x1
          x1=tmp
          xb=xa
        enddo
        signy=-signy
        tmp=y2
        y2=y1
        y1=tmp
        yb=ya
      enddo
      ngr=ij-1
      end

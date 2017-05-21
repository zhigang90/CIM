c===============================================================================
      subroutine integral_chunk(ish_tab_ik, inx,   xints,   xinta, buf,
     *                          xbuf1,   xbuf2, schwartz,bl,
     *                          ikstart,ikstop,iksize,jlstart,jlstop,
     *                          jlsize,ncs,iskip,thresh,time,
     *                          scr_s,scr_a,nskip,isparses,isparsea,
     *                          ish_tab_jl)
      implicit none
      integer inx(12,*),ncs,iskip,ifirst,jfirst,nskip
      integer*2 ish_tab_jl(2,*),ish_tab_ik(2,*)
      real*8 thresh,xints(*),xinta(*),buf(*),xbuf1(*),xbuf2(*),time
      real*8 bl,scr_s(*),scr_a(*)
      integer ikstart,ikstop,iksize
      integer jlstart,jlstop,jlsize
      integer*1 isparses(*),isparsea(*)
c     common/eeo_integ_pass/ikstart,ikstop,iksize,jlstart,jlstop,jlsize,
c    *                      ncs,iskip,thresh,time
      real*8 schwartz(ncs,ncs)
c
      integer max_shell_size
      parameter (max_shell_size=15)
      integer istart,istop,jstart,jstop,kstart,kstop,lstart,lstop
      integer isize,jsize,ksize,lsize,i,j,k,l,ik,jl
      integer js,ls,kls,kjs,jkls,lkjs,ip,jp,kp,lp,jcs,lcs,ics,kcs
      real*8 eps,perm1,perm2,perm,thresh2,t0,t1,s,a
      integer jlc,jls,jlbegin,ikc,iks,ikbegin,jlcst
      integer kpjs,kpls,jpkls,ik_1_jlsize,ipjkls,jl_abs
c
      logical nonz1,nonz2,i1,i2
c
      call dynamic_mmark()
c
      eps=thresh*1.d-1
      thresh2=thresh*thresh
c
      call zeroit(xints,iksize*jlsize)
      call zeroit(xinta,iksize*jlsize)
      ikbegin=0
      do ikc=ikstart,ikstop
        ics=ish_tab_ik(1,ikc)
        kcs=ish_tab_ik(2,ikc)
        istart=inx(11,ics)+1
        istop =inx(10,ics)
        isize=inx(3,ics)
c
        kstart=inx(11,kcs)+1
        kstop =inx(10,kcs)
        ksize=inx(3,kcs)
c
        if (ics.ne.kcs) then
          iks=isize*ksize
        else
          iks=isize*(isize+1)/2
        endif
c
        jlbegin=0
        jlc=0
        do jlc=jlstart,jlstop
          jcs=ish_tab_jl(1,jlc)
          lcs=ish_tab_jl(2,jlc)
          jstart=inx(11,jcs)+1
          jstop =inx(10,jcs)
          jsize=inx(3,jcs)
          lstart=inx(11,lcs)+1
          lstop =inx(10,lcs)
          lsize=inx(3,lcs)
c       do jcs=1,ncs
c         jstart=inx(11,jcs)+1
c         jstop =inx(10,jcs)
c         jsize=inx(3,jcs)
c       do lcs=1,jcs
c         jlc=jlc+1
c         if (jlc.lt.jlstart) cycle
c         if (jlc.gt.jlstop ) exit
c         lstart=inx(11,lcs)+1
c         lstop =inx(10,lcs)
c         lsize=inx(3,lcs)
          if (jcs.ne.lcs) then
            jls=jsize*lsize
          else
            jls=jsize*(jsize+1)/2
          endif
          i1=.true.
          i2=.true.
          if (sqrt(schwartz(ics,jcs)*schwartz(kcs,lcs)).lt.thresh) 
     *        i1=.false.
          if (sqrt(schwartz(ics,lcs)*schwartz(kcs,jcs)).lt.thresh)
     *        i2=.false.
          if (i1) then
            call secund(t0)
            call twoelmain(ics,jcs,kcs,lcs,eps,nonz1,xbuf1,buf)
            call secund(t1)
            time=time+t1-t0
          else
            nonz1=.false.
            iskip=iskip+1
          endif
          if (i2) then
            if (jcs.eq.lcs) then
              call tfer(xbuf1,xbuf2,isize*jsize*ksize*lsize)
              nonz2=nonz1
            else
              call secund(t0)
              call twoelmain(ics,lcs,kcs,jcs,eps,nonz2,xbuf2,buf)
              call secund(t1)
              time=time+t1-t0
            endif
          else
            nonz2=.false.
            iskip=iskip+1
          endif
          if (.not.(nonz1.or.nonz2)) then 
            jlbegin=jlbegin+jls
            nskip=nskip+iks*jls*2   ! one "shellquad" skipped
            cycle
          endif
          if (.not.nonz1) call zeroit(xbuf1,isize*jsize*ksize*lsize)
          if (.not.nonz2) call zeroit(xbuf2,isize*jsize*ksize*lsize)
          call array_files(23)
c
          perm1=perm(ics,jcs,kcs,lcs)
          perm2=perm(ics,lcs,kcs,jcs)
          js=jsize
          ls=lsize
          kls =ksize*lsize
          kjs =ksize*jsize
          jkls=jsize*ksize*lsize
          lkjs=jkls
c
          ik=ikbegin
          do i=istart,istop
            ip=i-istart
            ipjkls=ip*jkls
            do k=kstart,kstop
              kp=k-kstart
              if (k.gt.i.and.ics.eq.kcs) exit
              ik=ik+1
              ik_1_jlsize=(ik-1)*jlsize
              jl=jlbegin
              kpls=kp*ls
              kpjs=kp*js
              do j=jstart,jstop
                jp=j-jstart
                jpkls=jp*kls
                do l=lstart,lstop
                 lp=l-lstart
                 if (l.gt.j.and.lcs.eq.jcs) exit
                 jl=jl+1
                 if (j.ge.l) then
                   jl_abs=j*(j-1)/2+l
                 else
                   jl_abs=l*(l-1)/2+j
                 endif
                 if (jl.gt.jlsize) STOP 'eeee'
                 if (j.ne.l) then
                   s=(xbuf1(ipjkls+jpkls+kpls+lp+1)*perm1+
     *                xbuf2(ipjkls+lp*kjs+kpjs+jp+1)*perm2)
                   a=(xbuf1(ipjkls+jpkls+kpls+lp+1)*perm1-
     *                xbuf2(ipjkls+lp*kjs+kpjs+jp+1)*perm2)
                 else
                   s=xbuf1(ipjkls+jpkls+kpls+lp+1)*perm1
                   a=0d0
                 endif
                 if (dabs(s*scr_s(jl_abs)).lt.thresh) then
                   xints(jl+ik_1_jlsize)=0d0
                   nskip=nskip+1
                 else
                   xints(jl+ik_1_jlsize)=s
                   isparses(ik)=1
                 endif
                 if (dabs(a*scr_a(jl_abs)).lt.thresh) then
                   xinta(jl+ik_1_jlsize)=0d0
                   nskip=nskip+1
                 else
                   xinta(jl+ik_1_jlsize)=a
                   isparsea(ik)=1
                 endif
                enddo
              enddo
            enddo
          enddo
c
          jlbegin=jlbegin+jls
c     enddo
      enddo
c
          ikbegin=ikbegin+iks
      enddo
      call dynamic_retmark()
      end
c===============================================================================
      function perm(i,j,k,l)
      implicit none
      real*8 perm
      integer i,j,k,l
      perm=1d0
      if (i.eq.j) perm=perm*2d0
      if (k.eq.l) perm=perm*2d0
      if (i.eq.k.and.j.eq.l.or.i.eq.l.and.j.eq.k) perm=perm*2d0
      end
c===============================================================================

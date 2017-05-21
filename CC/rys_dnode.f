c==============================================================================
      subroutine dnode(a,rt,k)
      implicit double precision (a-h,o-z)
c        *****  routine returns in rt(i) the ith root of        *****
c        *****  a polynomial of order k whose mth coefficient   *****
c        *****  is stored in a(m+1).  it is assumed that the    *****
c        *****  initial values in rt bracket the final values.  *****
      dimension a(10),rt(10)
      logical xone
c
c **** this program originally was caught in an infinite loop
c      between statements 45 and 60. a change in the condition
c      at the statement if... go to 60 was made which cured the
c      problem but the accuracy may be not as high as desirable
c      P. Pulay, march 13, 1987.
c      Note 9/20/2005: Now I get another infinite loop, between
c      "if(prod.lt.tol) go to 30" and ? for n=9
c     tol = 2.0d-15
c   Limit the number of cycles
      maxcycle=500
      tol= 1.0d-9
      k1=k+1
      r2=0.0d+00
      p2=a(1)
      icycle=0
      do 100 m=1,k
      r1=r2
      p1=p2
      r2=rt(m)
      p2=a(k1)
      do 10 i=1,k
   10 p2=p2*r2+a(k1-i)
      prod=p1*p2
      if(prod.lt.0.0d+00) go to 20
      call nerror(1,'dnode',' was unable to find root no. ',m,k)
c      write(8,15) m,k
c      write(9,15) m,k
   15 format(/' root number',i3,' was not found for polynomial of order'
     1 ,i3,/)
      call exit
   20 r5=r1
      p5=p1
      r6=r2
      p6=p2
   30 r3=r5
      p3=p5
      r4=r6
      p4=p6
      r =(r3*p4-r4*p3)/(p4-p3)
      dr=r4-r3
      delta=dr
      if(dabs(delta).lt.tol) go to 90
      dr=0.0625d+00*dr
      r5=r-dr
      if(r5.lt.r3) r5=r3
      r6=r+dr
      if(r6.gt.r4) r6=r4
      p5=a(k1)
      p6=p5
      do 40 i=1,k
      p5=p5*r5+a(k1-i)
   40 p6=p6*r6+a(k1-i)
   45 prod=p5*p6
CPP
      icycle=icycle+1
      if(icycle.gt.maxcycle) then
c        print *, 'Stopping at root',m,' cycle=',icycle
      call nerror(2,'dnode','maximum iterations reached, root, cycle=',
     1            m,icycle)
        rt(m)=r
        return
      end if
      if(prod.lt.tol) go to 30
      prod=p3*p5
      if(prod.gt.tol) go to 60
      r5=0.25d+00*r3+0.75d+00*r5
      p5=a(k1)
      do 50 i=1,k
   50 p5=p5*r5+a(k1-i)
      go to 45
   60 r6=0.25d+00*r4+0.75d+00*r6
      p6=a(k1)
      do 70 i=1,k
   70 p6=p6*r6+a(k1-i)
      go to 45
   90 rt(m)=r
  100 continue
c  Polish roots
      do m=1,k
        x=rt(m)
        xone=.false.
        if (dabs(x).ge.1d0) xone=.true.
        icount=0
1000    p=a(k+1)*x+a(k)
        p1=a(k+1)
        do i=k-1,1,-1
          p1=p+p1*x
          p=a(i)+p*x
        end do
        del=p/p1
        icount=icount+1
        if (xone) then
          xdiff=abs(del/x)
        else
          xdiff=abs(del)
        endif
        if(xdiff.gt.1d-13) then
c          print *, k,m,del
          x=x-del
          if (icount.le.4) go to 1000
        end if
      rt(m)=x
      end do
      return
      end

c==============================================================================

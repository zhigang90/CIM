      subroutine grpsmb(sflies,csf,cla,nn,nprt)
      implicit real*8(a-h,o-z)
      character csf,cla,zn(0:9),sflies*4
      data zn /'0','1','2','3','4','5','6','7','8','9'/
c
      nn=0
c     nn is the symmetry number of the z-axis
      csf=sflies(1:1)
      do 5 i=1,9
         if (sflies(2:2).eq.zn(i)) nn=i
    5    continue
      if (nn.eq.0) then
         cla=sflies(2:2)
         if (nprt.gt.0) write (6,'(/,
     1      '' symmetry group of the molecule :   '',2a1)') csf,cla
         if (csf.eq.'c'.and.(cla.eq.'s'.or.(nn.eq.0.and.cla.eq.'h')))
     1   then
            nn=1
            cla='h'
         elseif (csf.eq.'c' .and. cla.eq.'i') then
            nn=2
            cla=' '
            csf='s'
         endif
      else
         np=-1
         do 10 i=0,9
            if (sflies(3:3).eq.zn(i)) np=i
   10       continue
         if (np.ge.0) then
            nn=10*nn+np
            if ( csf.eq.'s'.and. mod(nn,2).eq.1 ) then
               csf = 'c'
               cla = 'h'
            else
               cla=sflies(4:4)
            endif
            if (nprt.gt.0)
     1        write (6,'(/,'' symmetry group of the molecule :   '',
     2           a1,i2,a1)') csf,nn,cla
         else
            if ( csf.eq.'s'.and. mod(nn,2).eq.1 ) then
               csf = 'c'
               cla = 'h'
            else
               cla=sflies(3:3)
            endif
            if (nprt.gt.0)
     1        write (6,'(/,'' symmetry group of the molecule :   '',
     2           a1,i1,a1)') csf,nn,cla
         endif
      endif
      return
      end
c =========================================================================
c
      subroutine getgen (csf,nn,cla,gen,ngen,ierr,nprt)
      implicit real*8(a-h,o-y)
      implicit character (z)
      dimension gen(9,3),unit(9)
      character csf,cla
      data zc,zd,zs,zt,zo,zi,zv,zh/'c','d','s','t','o','i','v','h'/
      a0=0.0d0
      a1=1.0d0
      a2=2.0d0
      a3=3.0d0
      a4=4.0d0
      a5=5.0d0
      ierr=0
      pi=a4*atan(a1)
      twopi=pi+pi
c -- generate unit matrix
      do 10 i=2,8
  10  unit(i) = a0
      unit(1) = a1
      unit(5) = a1
      unit(9) = a1
      ngen=0
c     write(6,*) ' this program sets generators for point groups'
c     write(6,*) ' specified by schoenflies symbols, e.g. d2d'
c     write(6,*) ' please use s1 for cs  and s2 for ci'
c     write(6,*) ' main axis is always the  z - axis'
c     write(6,*) ' secondary axis is always x, e.g. for d3'
c     write(6,*) ' sigma(v) plane is always xz, e.g. for c3v'
c     write(6,*) ' input  schoenflies symbol: c2v  or td etc'
      if (nprt.gt.0) write(6,20)
   20 format(/,' the group has the following generators :')
c     note that in ih the mirror plane is perpendicular to the c2-axis
c     but not to the c5- or c3-axis
      if (csf.eq.zi.and.cla.eq.zh) cla=zv
c
c       set generators for input group
c
c       cn(z)
c
      n=0
      if (nn.gt.0.and.(csf.eq.zc.or.csf.eq.zd).and.cla.ne.zd) n=nn
      if (csf.eq.zt.and.cla.ne.zd) n=2
      if (csf.eq.zo) n=4
      if(csf.eq.zi) n=5
      if (n.le.0) go to 250
      if (n.le.9 .and. nprt.gt.0) write(6,30) n
      if (n.gt.9 .and. nprt.gt.0) write(6,31) n
   30 format(3x,'c',i1,'(z)')
   31 format(3x,'c',i2,'(z)')
      ngen = ngen+1
      if (ngen.gt.3) go to 1000
      do 200 i=1,9
  200 gen(i,ngen) = unit(i)
      an = n
      alpha = twopi/an
      x = sin(alpha)
      y = cos(alpha)
      gen(1,ngen) = y
      gen(5,ngen) = y
      gen(4,ngen) = -x
      gen(2,ngen) = x
  250 continue
c
c       sn(z)
c
      n=0
      if (csf.eq.zs) n=nn
      if (csf.eq.cla) n=2*nn
      if (csf.eq.zt.and.cla.eq.zd) n=4
      if (n.le.0) go to 400
      if (n.le.9 .and. nprt.gt.0) write(6,300) n
      if (n.gt.9 .and. nprt.gt.0) write(6,301) n
  300 format(3x,'s',i1,'(z)')
  301 format(3x,'s',i2,'(z)')
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 350 i=1,9
  350 gen(i,ngen) =-unit(i)
      an = n
      alpha = twopi/an
      x = sin(alpha)
      y = cos(alpha)
      gen(1,ngen) = y
      gen(5,ngen) = y
      gen(4,ngen) = -x
      gen(2,ngen) = x
  400 continue
c
c        c3(1,1,1)  for t,td,th,o,oh
c
      if (csf.ne.zt.and.csf.ne.zo) go to 450
      if (nprt.gt.0) write(6,410)
  410 format(3x,'c3(1,1,1)')
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 420 i=1,9
  420 gen(i,ngen) = a0
      gen(2,ngen) = a1
      gen(6,ngen) = a1
      gen(7,ngen) = a1
  450 continue
c
c       c3 for i or ih
c
      if (csf.ne.zi) go to 500
      if (nprt.gt.0) write(6,480)
  480 format(3x,'c3 for i or ih')
      ngen = ngen + 1
      if (ngen.gt.3) go to 1000
      sqr3 = sqrt(a3)
      cosd = sin(a2*pi/a5)/((a1-cos(a2*pi/a5))*sqr3)
      sind = sqrt(a1-cosd**2)
      gen(1,ngen) = a1 - a3*cosd**2/a2
      gen(2,ngen) = sqr3*cosd/a2
      gen(3,ngen) = a3*sind*cosd/a2
      gen(4,ngen) = -gen(2,ngen)
      gen(5,ngen) = -a1/a2
      gen(6,ngen) = sqr3*sind/a2
      gen(7,ngen) = gen(3,ngen)
      gen(8,ngen) = -gen(6,ngen)
      gen(9,ngen) = a1-a3*sind**2/a2
  500 continue
c
c       c2(x)
c
      if(csf.ne.zd) go to 600
      if (nprt.gt.0) write(6,520)
  520 format(3x,'c2(x)')
      i=1
      ngen = ngen +1
      if (ngen.gt.3) go to 1000
      do 550 j=1,9
  550 gen(j,ngen) = -unit(j)
      gen(4*i-3,ngen) = a1
  600 continue
c
c       mirror plane sigma(xz) or sigma(xy)
c
      if (cla.ne.zv.and.cla.ne.zh) go to 750
      if (cla.eq.zv) then
         if (nprt.gt.0) write(6,610)
  610    format(3x,'mirror plane sigma(xz)')
         i=2
      else
         if (nprt.gt.0) write(6,620)
  620    format(3x,'mirror plane sigma(xy)')
         i=3
      endif
      ngen = ngen + 1
      if (ngen.gt.3) go to 1000
      do 700 j=1,9
  700 gen(j,ngen) = unit(j)
      gen(4*i-3,ngen) = -a1
  750 return
 1000 write(6,*) ' more than ',3,' generators'
      ierr=1
      end
c =========================================================================
c
      subroutine groupL(gen,    ngen,   u,      nt,     max14,
     $                  thrsym, inv,    ipath,  jpath,  iprnt,
     1                  ierr)
      implicit real*8(a-h,o-z)
      dimension u(9,max14),gen(9,ngen),inv(max14),scr(9)
      dimension ipath(max14),jpath(max14)
c ...............................................................
c      generates 3x3 matrices forming group
c      from ngen generators on gen
c      output:
c      u = 3x3 matrices
c      nt = number of symmetry operations
c      max14 = maximum order of group
c      inv  containes index of invers operation
c ...............................................................
      a0=0.0d0
      a1=1.0d0
      ierr=0
c -- set identity element
c -- (every point group has the identity as first operation)
      nt = 1
      do 200 i=2,8
 200  u(i,nt) = a0
      u(1,nt) = a1
      u(5,nt) = a1
      u(9,nt) = a1
      ipath(1)=0
      jpath(1)=0
      if (ngen.eq.0) goto 700
      do 500 igen=1,ngen
c -- loop over generators
         ncoset = 1
         nto = nt
         do 220 i=1,9
 220        scr(i) = gen(i,igen)
c -- check if generator is already contained in (sub)group
c -- of order nto :
         do 240 j=1,nt
            if(dif(u(1,j),scr(1),9).lt.thrsym) go to 500
 240        continue
         ipath(nt+1)=-igen
         jpath(nt+1)=0
c -- for new group element form coset
 245     do 250 it=1,nto
            nt = nt + 1
            if (nt.gt.max14) goto 900
            if (it.eq.1) goto 250
            ipath(nt)=nto*ncoset+1
            jpath(nt)=it
 250        call mult3 (u(1,nto*ncoset+it),scr,u(1,it),3)
         ncoset = ncoset + 1
c -- check : do subgroup plus cosets form a group ?
         do 330 it = 1,nt
            do 320 jcoset = 1,ncoset-1
               call mult3 (scr,u(1,it),u(1,1+nto*jcoset),3)
               do 310 kt = 1,nt
                  if (dif(u(1,kt),scr(1),9).lt.thrsym) goto 320
 310              continue
c -- new symmetry operation found
               ipath(nt+1)=it
               jpath(nt+1)=1+nto*jcoset
               goto 245
 320           continue
 330        continue
c -- cosets now form group - take next generator
 500     continue
c
c -- find inverse operators
      inv(1) = 1
      if (nt.eq.1) goto 700
      do 660 i=2,nt
      inv(i) = 0
      do 640 j=2,nt
      call mult3 (scr,u(1,i),u(1,j),3)
      s = dif (scr,u(1,1),9)
      if (s.gt.thrsym) go to 640
      inv(i) = j
      go to 660
  640 continue
      go to 950
  660 continue
  700 continue
c
      if (iprnt.ge.2) then
         write(6,750)
  750    format(/,10x,'xx',6x,'yx',6x,'zx',6x,'xy',6x,'yy',6x,
     1          'zy',6x,'xz',6x,'yz',6x,'zz')
         do 5 it=1,nt
    5       write(6,'(i4,1x,9f8.4)') it,(u(i,it),i=1,9)
         write(6,*)
         write(6,'(/,'' group element generation record'',
     1          '' (negative numbers correspond to generators) :'')')
         write (6,'(6(2x,i3,''='',i3,''*'',i3))')
     1            (it,ipath(it),jpath(it),it=1,nt)
      end if
c
      return
c
  900 write(6,'(a,i4)') ' size of group larger than parameter ndi14 =',
     1                  max14
      ierr=1
      return
c
  950 write(6,*) ' inversion problem '
      ierr=1
      return
      end
c =========================================================================
c
      subroutine getrep(csf,nn,cla,nsym,ityp,idim,max17,nsymsq,nsumid)
      implicit real*8(a-h,o-y)

      implicit character (z)
      character csf,cla,ityp*4
      dimension ityp(max17),idim(max17),zn(0:9)
      data zn /'0','1','2','3','4','5','6','7','8','9'/

      z1='g'
      z2='u'

      if (csf.eq.'s' .and. mod(nn,2).eq.1) then
         csf='c'
         cla='h'
      endif

      zsf=csf
      nnz=nn
      zla=cla

      if (csf.eq.'c' .or. csf.eq.'d' .or. csf.eq.'s') then
c     non-polyhedral point groups :
         if (cla.eq.'h' .and. mod(nn,2).eq.1) then
            z1=''''
            z2='"'
         endif
         if (cla.eq.'d') then
            if (mod(nn,2).eq.1) then
               zla='h'
            else
               nnz=2*nn
               zla='v'
            endif
         elseif (csf.eq.'s') then
            if (mod(nn,4).eq.2) then
               nnz=nn/2
               zla='h'
            else
               zsf='c'
            endif
         endif
         if (zsf.eq.'d' .and. nnz.eq.2) then
            nsym=4
            if (nsym.gt.max17) then
               write(6,*)
               write(6,*) ' more than ',max17,' representations '
               write(6,*) ' increase parameter ndi17 to ',nsym
               write(6,*)
               call nerror(1,'getrep','error',0,0)
            endif
            ityp(1)='a   '
            ityp(2)='b1  '
            ityp(3)='b2  '
            ityp(4)='b3  '
         elseif (zsf.eq.'d' .or. zla.eq.'v') then
            ityp(1)='a1  '
            ityp(2)='a2  '
            if (mod(nnz,2).eq.1) then
               nsym=(nnz+3)/2
               if (nsym.gt.max17) then
                  write(6,*)
                  write(6,*) ' more than ',max17,' representations '
                  write(6,*) ' increase parameter ndi17 to ',nsym
                  write(6,*)
                  call nerror(1,'getrep','error',0,0)
               endif
               if (nsym.eq.3) then
                  ityp(3)='e   '
               elseif (nsym.gt.3) then
                  do 1200 i=3,nsym
                     if (i-2.le.9) then
                        ityp(i)='e'//zn(i-2)//'  '
                     else
                        ityp(i)='e'//zn((i-2)/10)//zn(mod(i-2,10))//' '
                     endif
 1200                continue
               endif
            else
               nsym=nnz/2+3
               if (nsym.gt.max17) then
                  write(6,*)
                  write(6,*) ' more than ',max17,' representations '
                  write(6,*) ' increase parameter ndi17 to ',nsym
                  write(6,*)
                  call nerror(1,'getrep','error',0,0)
               endif
               ityp(3)='b1  '
               ityp(4)='b2  '
               if (nsym.eq.5) then
                  ityp(5)='e   '
               elseif (nsym.gt.5) then
                  do 1400 i=5,nsym
                     if (i-4.le.9) then
                        ityp(i)='e'//zn(i-4)//'  '
                     else
                        ityp(i)='e'//zn((i-4)/10)//zn(mod(i-4,10))//' '
                     endif
 1400                continue
               endif
            endif
         else
            ityp(1)='a   '
            if (mod(nnz,2).eq.1) then
               nsym=(nnz+1)/2
               if (nsym.gt.max17) then
                  write(6,*)
                  write(6,*) ' more than ',max17,' representations '
                  write(6,*) ' increase parameter ndi17 to ',nsym
                  write(6,*)
                  call nerror(1,'getrep','error',0,0)
               endif
               if (nsym.eq.2) then
                  ityp(2)='e   '
               elseif (nsym.gt.2) then
                  do 1600 i=2,nsym
                     if (i-1.le.9) then
                        ityp(i)='e'//zn(i-1)//'  '
                     else
                        ityp(i)='e'//zn((i-1)/10)//zn(mod(i-1,10))//' '
                     endif
 1600                continue
               endif
            else
               nsym=nnz/2+1
               if (nsym.gt.max17) then
                  write(6,*)
                  write(6,*) ' more than ',max17,' representations '
                  write(6,*) ' increase parameter ndi17 to ',nsym
                  write(6,*)
                  call nerror(1,'getrep','error',0,0)
               endif
               ityp(2)='b   '
               if (nsym.eq.3) then
                  ityp(3)='e   '
               elseif (nsym.gt.3) then
                  do 1800 i=3,nsym
                     if (i-2.le.9) then
                        ityp(i)='e'//zn(i-2)//'  '
                     else
                        ityp(i)='e'//zn((i-2)/10)//zn(mod(i-2,10))//' '
                     endif
 1800                continue
               endif
            endif
         endif
      elseif (csf.eq.'t' .or. csf.eq.'o') then
c     cubic point groups :
         if ((zsf.eq.'t' .and. zla.eq.'d') .or. zsf.eq.'o') then
            nsym=5
            if (nsym.gt.max17) then
               write(6,*)
               write(6,*) ' more than ',max17,' representations '
               write(6,*) ' increase parameter ndi17 to ',nsym
               write(6,*)
               call nerror(1,'getrep','error',0,0)
            endif
            ityp(1)='a1  '
            ityp(2)='a2  '
            ityp(3)='e   '
            ityp(4)='t1  '
            ityp(5)='t2  '
         else
            nsym=3
            if (nsym.gt.max17) then
               write(6,*)
               write(6,*) ' more than ',max17,' representations '
               write(6,*) ' increase parameter ndi17 to ',nsym
               write(6,*)
               call nerror(1,'getrep','error',0,0)
            endif
            ityp(1)='a   '
            ityp(2)='e   '
            ityp(3)='t   '
         endif
      elseif (csf.eq.'i') then
c     icosahedral point groups :
         nsym=5
         if (nsym.gt.max17) then
            write(6,*)
            write(6,*) ' more than ',max17,' representations '
            write(6,*) ' increase parameter ndi17 to ',nsym
            write(6,*)
            call nerror(1,'getrep','error',0,0)
         endif
         ityp(1)='a   '
         ityp(2)='t1  '
         ityp(3)='t2  '
         ityp(4)='g   '
         ityp(5)='h   '
      endif
c     extra representations in presence of sigma-h mirror planes :
      if (zla.eq.'h' .or. csf.eq.'i' .and. cla.eq.'v') then
         if (2*nsym.gt.max17) then
            write(6,*)
            write(6,*) ' more than ',max17,' representations '
            write(6,*) ' increase parameter ndi17 to ',2*nsym
            write(6,*)
            call nerror(1,'getrep','error',0,0)
         endif
         do 3000 i=1,nsym
            ityp(i+nsym)=ityp(i)
            if (ityp(i)(2:2).eq.' ') then
               ityp(i)(2:2)=z1
               ityp(i+nsym)(2:2)=z2
            elseif (ityp(i)(3:3).eq.' ') then
               ityp(i)(3:3)=z1
               ityp(i+nsym)(3:3)=z2
            else
               ityp(i)(4:4)=z1
               ityp(i+nsym)(4:4)=z2
            endif
 3000       continue
         nsym=2*nsym
      endif
      do 4000 i=1,nsym
         if (ityp(i)(1:1).eq.'a' .or. ityp(i)(1:1).eq.'b') idim(i)=1
         if (ityp(i)(1:1).eq.'e') idim(i)=2
         if (ityp(i)(1:1).eq.'t') idim(i)=3
         if (ityp(i)(1:1).eq.'g') idim(i)=4
         if (ityp(i)(1:1).eq.'h') idim(i)=5
 4000    continue

c     output real representations of symmetry group :

      if(nsym.le.8) then
        write(6,5000) nsym,(ityp(i),i=1,nsym)
 5000   format(/,' There are ',i1,' real representations:  ',8(1x,a4))
      else
        write(6,5100) nsym
 5100   format(/,' There are ',i2,' real representations:')
        write(6,5200) (ityp(i),i=1,nsym)
 5200   format(15(1x,a4))
      endif

      nsumid=0
      nsymsq=0
      do 6000 i=1,nsym
         nsumid=nsumid+idim(i)
 6000    nsymsq=nsymsq+idim(i)*idim(i)

      return
      end
c =========================================================================
c
      subroutine grep(nt,ipath,jpath,nsym,nsymsq,idim,charac,
     1                repmat,genrep)
      implicit real*8(a-h,o-z)

c =====================================================================
c
c      calculates all representation matrices from representation
c      matrices of generating symmetry operations (on genrep)
c      repmat = matrix elements of representation matrices
c      charac = characters of real representations
c       index of repmat :
c             i+(j-1)*idim(isym)+(itrans-1)*idim(isym)**2+nt*isymsq
c             where   itrans = index of symmetry operation
c                      isym  = index of representation
c                      i,j   = rep.matrix indices
c                  idim(isym)= dimension of representation
c                     isymsq = sum over jsym<isym (idim(jsym)**2)
c                       nt   = number of symmetry operations
c      nsymsq = sum over all jsym (idim(jsym)**2)
c      ipath,jpath = defines the path of group element generation used
c                    in subroutine group
c
c =====================================================================

      dimension ipath(nt),jpath(nt),idim(nsym),charac(nt*nsym),
     1          repmat(nt*nsymsq),genrep(3*nsymsq)

      a0=0
      a1=1

c     --- initialize repmat
      Call ZeroIT(repmat,nsymsq*nt)

c     --- loop over all representations
      isymsq = 0
      do 8000 isym=1,nsym
      ntidsq = nt*isymsq
      idisym=idim(isym)
c      set neutral element :
c      (every point group has the neutral element as first operation)
      itrans = 1
      do 1000 ij=1,idisym**2,idisym+1
 1000   repmat(ij+ntidsq)=a1
      charac((isym-1)*nt+1)=dble(idim(isym))
c     loop over symmetry operations:
      if (nt.eq.1) return
      do 5000 itrans=2,nt
         itidsq=(itrans-1)*idisym**2+ntidsq
         if (jpath(itrans).eq.0) then
c           (group element is a generator)
            igen=-ipath(itrans)
            do 2000 ij=1,idisym**2
 2000          repmat(ij+itidsq)=genrep(ij+isymsq+(igen-1)*nsymsq)
         else
c           (group element itrans is the product of ipath(itrans)
c            and jpath(itrans) )
            ltidsq=(ipath(itrans)-1)*idisym**2+ntidsq
            ktidsq=(jpath(itrans)-1)*idisym**2+ntidsq
            do 4000 i0=0,idisym*(idisym-1),idisym
               do 3750 k=1,idisym
                  ki=k+i0+ktidsq
                  k0=(k-1)*idisym+ltidsq
                  do 3500 j=1,idisym
                     repmat(j+i0+itidsq)=repmat(j+i0+itidsq)
     1                                  +repmat(j+k0)*repmat(ki)
 3500                continue
 3750             continue
 4000          continue
         endif
c        calculate characters :
         chara=0
         do 4500 ij=1,idisym**2,idisym+1
 4500       chara=chara+repmat(ij+itidsq)
         charac((isym-1)*nt+itrans)=chara
 5000    continue
 8000 isymsq=isymsq+idisym**2

      return
      end
c =========================================================================
c
      subroutine repgen (csf,nn,cla,gen,ngen,nsymsq,nsym,ityp,idim,
     1                   genrep)
      implicit real*8(a-h,o-y)
      implicit character (z)

c =====================================================================
c
c     this subroutine calculates the representation matrices genrep
c     of all real representations (ityp(1..nsym)) for the generating
c     elements of the group (gen(1..ngen))
c
c =====================================================================

      dimension gen(9,3),ityp(nsym),idim(nsym),zn(0:9),rmat(5,5),
     1          genrep(3*nsymsq)
      character csf,cla,ityp*4

      data zd,zs,zt,zo,zi,zv,zh/'d','s','t','o','i','v','h'/
      data zn /'0','1','2','3','4','5','6','7','8','9'/

      a0=0
      a1=1
      a2=2
      a3=3
      a4=4
      a5=5
      pi=a4*atan(a1)
      twopi=pi+pi

      do 180 i=0,ngen-1
  180    genrep(1+i*nsymsq)=a1
      if (nsym.eq.1) return

      igen=0
c     non-polyhedral point groups :
      if (csf.eq.zt .or. csf.eq.zo .or. csf.eq.zi) goto 3000
      n=nn
c
c       cn(z)
c
      if (csf.eq.zs .or. csf.eq.cla) goto 250
      igen = igen +1
      isymsq=1
      do 240 isym=2,nsym
         if (ityp(isym)(1:1).eq.'a') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (ityp(isym)(1:1).eq.'b') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=-a1
            if (csf.eq.zd .and. nn.eq.2 .and. ityp(isym)(2:2).eq.'1')
     1         genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (ityp(isym)(1:1).eq.'e') then
            esym=a1
            do 230 i=2,9
               if (ityp(isym)(2:2).eq.zn(i)) esym=dble(i)
  230          continue
            do 235 i=0,9
               if (ityp(isym)(3:3).eq.zn(i)) then
                do 234 j=1,9
                  if (ityp(isym)(2:2).eq.zn(j)) esym=dble(10*j+i)
  234             continue
               endif
  235          continue
            an = n
            alpha = twopi*esym/an
            x = sin(alpha)
            y = cos(alpha)
            genrep(isymsq+1+(igen-1)*nsymsq)=y
            genrep(isymsq+2+(igen-1)*nsymsq)=x
            genrep(isymsq+3+(igen-1)*nsymsq)=-x
            genrep(isymsq+4+(igen-1)*nsymsq)=y
            isymsq=isymsq+4
         endif
  240    continue
  250 continue
c
c       sn(z)
c
      if (csf.ne.zs .and. cla.ne.csf) goto 400
      if (csf.eq.zs) n=nn/2
      if (n.eq.(n/2)*2) n=2*n
      if (csf.eq.cla) nn=2*nn
      igen = igen +1
      isymsq=1
      do 340 isym=2,nsym
         if (ityp(isym)(1:1).eq.'a') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
            if (ityp(isym)(2:2).eq.'u' .or. ityp(isym)(3:3).eq.'u')
     1         genrep(isymsq+(igen-1)*nsymsq)=-a1
         elseif (ityp(isym)(1:1).eq.'b') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=-a1
         elseif (ityp(isym)(1:1).eq.'e') then
            esym=a1
            do 330 i=2,9
               if (ityp(isym)(2:2).eq.zn(i)) esym=dble(i)
  330          continue
            do 335 i=0,9
               if (ityp(isym)(3:3).eq.zn(i)) then
                do 334 j=1,9
                  if (ityp(isym)(2:2).eq.zn(j)) esym=dble(10*j+i)
  334             continue
               endif
  335          continue
            an = n
            if (ityp(isym)(2:2).eq.'u' .or. ityp(isym)(3:3).eq.'u'
     1          .or. ityp(isym)(4:4).eq.'u')
     1         esym=(a2*esym-a1)/a2
            alpha = twopi*esym/an
            x = sin(alpha)
            y = cos(alpha)
            genrep(isymsq+1+(igen-1)*nsymsq)=y
            genrep(isymsq+2+(igen-1)*nsymsq)=x
            genrep(isymsq+3+(igen-1)*nsymsq)=-x
            genrep(isymsq+4+(igen-1)*nsymsq)=y
            isymsq=isymsq+4
         endif
  340    continue
  400 continue
c
c       c2(x)
c
      if(csf.ne.zd) go to 600
      igen = igen +1
      isymsq=1
      do 540 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
            if (ityp(isym)(2:2).eq.'2' .or. (ityp(isym)(2:2).eq.'1'
     1          .and. nn.eq.2 .and. cla.ne.csf))
     2         genrep(isymsq+(igen-1)*nsymsq)=-a1
         elseif (ityp(isym)(1:1).eq.'e') then
            xy=a1
            if (cla.eq.zh) then
               if (ityp(isym)(2:2).eq.'"' .or. ityp(isym)(3:3).eq.'"'
     1             .or. ityp(isym)(4:4).eq.'"') then
                  xy=-a1
               elseif (nn.eq.(nn/2)*2) then
                  iesym=1
                  do 520 i=2,9
                     if (ityp(isym)(2:2).eq.zn(i)) iesym=i
  520                continue
                  do 525 i=0,9
                     if (ityp(isym)(3:3).eq.zn(i)) then
                      do 524 j=1,9
                        if (ityp(isym)(2:2).eq.zn(j)) iesym=10*j+i
  524                   continue
                     endif
  525                continue
                  if(mod(iesym,2).eq.0 .and. index(ityp(isym),'u').gt.1
     1                   .or.
     2               mod(iesym,2).eq.1 .and. index(ityp(isym),'g').gt.1)
     3                  xy=-a1
               endif
            endif
            genrep(isymsq+1+(igen-1)*nsymsq)=xy
            genrep(isymsq+2+(igen-1)*nsymsq)=a0
            genrep(isymsq+3+(igen-1)*nsymsq)=a0
            genrep(isymsq+4+(igen-1)*nsymsq)=-xy
            isymsq=isymsq+4
         endif
  540    continue
  600 continue
c
c       mirror plane sigma(xz)
c
      if (cla.ne.zv) go to 700
      igen = igen + 1
      isymsq=1
      do 640 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
            if (ityp(isym)(2:2).eq.'2')
     1         genrep(isymsq+(igen-1)*nsymsq)=-a1
         elseif (ityp(isym)(1:1).eq.'e') then
            genrep(isymsq+1+(igen-1)*nsymsq)=a1
            genrep(isymsq+2+(igen-1)*nsymsq)=a0
            genrep(isymsq+3+(igen-1)*nsymsq)=a0
            genrep(isymsq+4+(igen-1)*nsymsq)=-a1
            isymsq=isymsq+4
         endif
  640    continue
  700 continue
c
c       mirror plane sigma(xy)
c
      if (cla.ne.zh) goto 750
      igen = igen + 1
      isymsq=1
      do 740 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
            if (ityp(isym)(2:2).eq.'"' .or. ityp(isym)(3:3).eq.'"' .or.
     1          nn.eq.(nn/4)*4 .and.
     2          (ityp(isym)(2:2).eq.'u'.or.ityp(isym)(3:3).eq.'u') .or.
     3          nn.eq.(nn/4)*4+2 .and.ityp(isym)(1:1).eq.'b' .and.
     4          (ityp(isym)(2:2).eq.'g'.or.ityp(isym)(3:3).eq.'g') .or.
     5          nn.eq.(nn/4)*4+2 .and.ityp(isym)(1:1).eq.'a' .and.
     6          (ityp(isym)(2:2).eq.'u'.or.ityp(isym)(3:3).eq.'u'))
     7         genrep(isymsq+(igen-1)*nsymsq)=-a1
            if (nn.eq.2 .and. csf.eq.zd) then
               if (ityp(isym).eq.'b2g ' .or. ityp(isym).eq.'b3g ' .or.
     1            ityp(isym).eq.'b1u ' .or. ityp(isym).eq.'au  ') then
                  genrep(isymsq+(igen-1)*nsymsq)=-a1
               else
                  genrep(isymsq+(igen-1)*nsymsq)=a1
               endif
            endif
         elseif (ityp(isym)(1:1).eq.'e') then
            xy=a1
            if (index(ityp(isym),'"').gt.1) then
               xy=-a1
            elseif (nn.eq.(nn/2)*2) then
               iesym=1
               do 720 i=2,9
                  if (ityp(isym)(2:2).eq.zn(i)) iesym=i
  720             continue
               do 725 i=0,9
                  if (ityp(isym)(3:3).eq.zn(i)) then
                   do 724 j=1,9
                     if (ityp(isym)(2:2).eq.zn(j)) iesym=10*j+i
  724                continue
                  endif
  725             continue
               if (mod(iesym,2).eq.0 .and. index(ityp(isym),'u').gt.1
     1                  .or.
     2             mod(iesym,2).eq.1 .and. index(ityp(isym),'g').gt.1)
     3                 xy=-a1
            endif
            genrep(isymsq+1+(igen-1)*nsymsq)=xy
            genrep(isymsq+2+(igen-1)*nsymsq)=a0
            genrep(isymsq+3+(igen-1)*nsymsq)=a0
            genrep(isymsq+4+(igen-1)*nsymsq)=xy
            isymsq=isymsq+4
         endif
  740    continue
  750 return
 3000 continue
c     polyhedral point groups :
      if (csf.eq.zi) goto 6000
c     cubic point groups :
c
c       cn(z) or s4(z)
c
      igen = igen +1
      isymsq=1
      do 3500 isym=2,nsym
         if (ityp(isym)(1:2).eq.'a2') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=-a1
         elseif (ityp(isym)(1:1).eq.'a') then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (idim(isym).eq.2) then
            if (csf.eq.zt .and. cla.ne.zd) then
               genrep(isymsq+1+(igen-1)*nsymsq)=a1
               genrep(isymsq+2+(igen-1)*nsymsq)=a0
               genrep(isymsq+3+(igen-1)*nsymsq)=a0
               genrep(isymsq+4+(igen-1)*nsymsq)=a1
            else
               genrep(isymsq+1+(igen-1)*nsymsq)=a1
               genrep(isymsq+2+(igen-1)*nsymsq)=a0
               genrep(isymsq+3+(igen-1)*nsymsq)=a0
               genrep(isymsq+4+(igen-1)*nsymsq)=-a1
            endif
            isymsq=isymsq+4
         elseif ((ityp(isym)(1:2).eq.'t2'.and.cla.ne.zd) .or.
     1           (ityp(isym)(1:2).eq.'t1'.and.cla.eq.zd)) then
            do 3300 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=-gen(i,igen)
 3300          continue
            isymsq=isymsq+9
         else
            do 3400 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=gen(i,igen)
 3400          continue
            isymsq=isymsq+9
         endif
 3500    continue
c
c        c3(1,1,1)  for t,td,o,oh
c
      igen = igen +1
      isymsq=1
      do 3700 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (idim(isym).eq.2) then
            genrep(isymsq+1+(igen-1)*nsymsq)=-a1/a2
            genrep(isymsq+2+(igen-1)*nsymsq)=sqrt(a3)/a2
            genrep(isymsq+3+(igen-1)*nsymsq)=-sqrt(a3)/a2
            genrep(isymsq+4+(igen-1)*nsymsq)=-a1/a2
            isymsq=isymsq+4
         elseif (idim(isym).eq.3) then
            do 3600 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=gen(i,igen)
 3600          continue
            isymsq=isymsq+9
         endif
 3700    continue
c
c       mirror plane sigma(xy)
c
      if (cla.ne.zh) go to 3900
      igen = igen + 1
      isymsq=1
      do 3850 isym=2,nsym
         xy=a1
         if (ityp(isym)(2:2).eq.'u' .or. ityp(isym)(3:3).eq.'u') xy=-a1
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=xy
         elseif (idim(isym).eq.2) then
            genrep(isymsq+1+(igen-1)*nsymsq)=xy
            genrep(isymsq+2+(igen-1)*nsymsq)=a0
            genrep(isymsq+3+(igen-1)*nsymsq)=a0
            genrep(isymsq+4+(igen-1)*nsymsq)=xy
            isymsq=isymsq+4
         elseif (ityp(isym)(1:1).eq.'t') then
            do 3800 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=a0
 3800          continue
            genrep(isymsq+1+(igen-1)*nsymsq)=-xy
            genrep(isymsq+5+(igen-1)*nsymsq)=-xy
            genrep(isymsq+9+(igen-1)*nsymsq)=xy
            isymsq=isymsq+9
         endif
 3850    continue
 3900 return
 6000 continue
c     icosahedral point groups :
c
c       c5(z)
c
      c=sqrt(a2)*sqrt(a3)*sqrt(a5+sqrt(a5))*(a5+sqrt(a5))/(a3*a4*a5)
      igen = igen +1
      isymsq=1
      do 6300 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (ityp(isym)(1:2).eq.'t1') then
            do 6050 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=gen(i,igen)
 6050          continue
            isymsq=isymsq+9
         elseif (ityp(isym)(1:2).eq.'t2') then
            rmat(1,1)=a1
            rmat(1,2)=a0
            rmat(1,3)=a0
            rmat(2,1)=a0
            rmat(2,2)=-(sqrt(a5)+a1)/a4
            rmat(2,3)=sqrt(sqrt(a5)+5)*sqrt(a2)*(sqrt(a5)-a1)/(a2*a4)
            rmat(3,1)=a0
            rmat(3,2)=-rmat(2,3)
            rmat(3,3)=rmat(2,2)
            do 6100 i=1,3
               do 6100 j=1,3
 6100             genrep(isymsq+i+(j-1)*3+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+9
         elseif (idim(isym).eq.4) then
            rmat(1,1)=(sqrt(a5)-a1)/a4
            rmat(1,2)=(-sqrt(sqrt(a5)+a5)*sqrt(a2))/a4
            rmat(1,3)=a0
            rmat(1,4)=a0
            rmat(2,1)=-rmat(1,2)
            rmat(2,2)=rmat(1,1)
            rmat(2,3)=a0
            rmat(2,4)=a0
            rmat(3,1)=a0
            rmat(3,2)=a0
            rmat(3,3)=-(sqrt(a5)+a1)/a4
            rmat(3,4)=(sqrt(-sqrt(a5)+a5)*sqrt(a2))/a4
            rmat(4,1)=a0
            rmat(4,2)=a0
            rmat(4,3)=-rmat(3,4)
            rmat(4,4)=rmat(3,3)
            do 6150 i=1,4
               do 6150 j=1,4
 6150             genrep(isymsq+i+(j-1)*4+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+16
         elseif (idim(isym).eq.5) then
            rmat(1,1)=a1
            rmat(1,2)=a0
            rmat(1,3)=a0
            rmat(1,4)=a0
            rmat(1,5)=a0
            rmat(2,1)=a0
            rmat(2,2)=(sqrt(a5)-a1)/a4
            rmat(2,3)=(-sqrt(a2)*sqrt(sqrt(a5)+a5))/a4
            rmat(2,4)=a0
            rmat(2,5)=a0
            rmat(3,1)=a0
            rmat(3,2)=-rmat(2,3)
            rmat(3,3)=rmat(2,2)
            rmat(3,4)=a0
            rmat(3,5)=a0
            rmat(4,1)=a0
            rmat(4,2)=a0
            rmat(4,3)=a0
            rmat(4,4)=-(sqrt(a5)+a1)/a4
            rmat(4,5)=sqrt(a2)*sqrt(sqrt(a5)+a5)*(sqrt(a5)-a1)/(a2*a4)
            rmat(5,1)=a0
            rmat(5,2)=a0
            rmat(5,3)=a0
            rmat(5,4)=-rmat(4,5)
            rmat(5,5)=rmat(4,4)
            do 6200 i=1,5
               do 6200 j=1,5
 6200             genrep(isymsq+i+(j-1)*5+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+25
         endif
 6300    continue
c
c       c3 for i or ih
c
      igen = igen +1
      isymsq=1
      do 6600 isym=2,nsym
         if (idim(isym).eq.1) then
            isymsq=isymsq+1
            genrep(isymsq+(igen-1)*nsymsq)=a1
         elseif (ityp(isym)(1:2).eq.'t1') then
            do 6350 i=1,9
               genrep(isymsq+i+(igen-1)*nsymsq)=gen(i,igen)
 6350          continue
            isymsq=isymsq+9
         elseif (ityp(isym)(1:2).eq.'t2') then
            rmat(1,1)=(-sqrt(a5))/a5
            rmat(1,2)=(sqrt(-sqrt(a5)+a3)*sqrt(a5)*sqrt(a2))/(a2*a5)
            rmat(1,3)=-c*sqrt(-sqrt(a5)+a3)*sqrt(a3)*sqrt(a2)/a2
            rmat(2,1)=rmat(1,2)
            rmat(2,2)=(a2*sqrt(a5)+a5)/(a2*a5)
            rmat(2,3)=c*sqrt(a3)*(sqrt(a5)-a2)/a2
            rmat(3,1)=-rmat(1,3)
            rmat(3,2)=-rmat(2,3)
            rmat(3,3)=(-a1)/a2
            do 6400 i=1,3
               do 6400 j=1,3
 6400             genrep(isymsq+i+(j-1)*3+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+9
         elseif (idim(isym).eq.4) then
            rmat(1,1)=(a3*sqrt(a5)+a5)/(a4*a5)
            rmat(1,2)=c*sqrt(a3)*(sqrt(a5)-a1)/a4
            rmat(1,3)=(a3*sqrt(a5))/(a2*a5)
            rmat(1,4)=(-sqrt(-a2*sqrt(a5)+a5)*sqrt(a5))/(a2*a5)
            rmat(2,1)=-rmat(1,2)
            rmat(2,2)=(-sqrt(a5)+a1)/a4
            rmat(2,3)=(c*sqrt(a3))/a2
            rmat(2,4)=(c*sqrt(-a2*sqrt(a5)+a5)*sqrt(a3))/a2
            rmat(3,1)=rmat(1,3)
            rmat(3,2)=-rmat(2,3)
            rmat(3,3)=(-a3*sqrt(a5)+a5)/(a4*a5)
            rmat(3,4)=sqrt(-a2*sqrt(a5)+a5)*(sqrt(a5)+a5)/(a4*a5)
            rmat(4,1)=-rmat(1,4)
            rmat(4,2)=rmat(2,4)
            rmat(4,3)=-rmat(3,4)
            rmat(4,4)=(sqrt(a5)+a1)/a4
            do 6450 i=1,4
               do 6450 j=1,4
 6450             genrep(isymsq+i+(j-1)*4+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+16
         elseif (idim(isym).eq.5) then
            rmat(1,1)=(-a1)/a5
            rmat(1,2)=sqrt(a3)*(sqrt(a5)+a1)/(a2*a5)
            rmat(1,3)=a3*c*(a3*sqrt(a5)-a5)/(a2*a5)
            rmat(1,4)=a3*c*(-sqrt(a5)+a5)/(a2*a5)
            rmat(1,5)=sqrt(a3)*(sqrt(a5)-a1)/(a2*a5)
            rmat(2,1)=rmat(1,2)
            rmat(2,2)=(a2*sqrt(a5)+a1)/(a2*a5)
            rmat(2,3)=sqrt(a3)*c*(-a2*sqrt(a5)+a5)/(a2*a5)
            rmat(2,4)=sqrt(a3)*c*(-a3*sqrt(a5)+a5)/a5
            rmat(2,5)=a2/a5
            rmat(3,1)=-rmat(1,3)
            rmat(3,2)=-rmat(2,3)
            rmat(3,3)=(-a1)/a2
            rmat(3,4)=a0
            rmat(3,5)=sqrt(a3)*c*(-sqrt(a5)+a5)/a5
            rmat(4,1)=-rmat(1,4)
            rmat(4,2)=-rmat(2,4)
            rmat(4,3)=a0
            rmat(4,4)=(-a1)/a2
            rmat(4,5)=(-sqrt(a3)*sqrt(a5)*c)/(a2*a5)
            rmat(5,1)=rmat(1,5)
            rmat(5,2)=rmat(2,5)
            rmat(5,3)=-rmat(3,5)
            rmat(5,4)=-rmat(4,5)
            rmat(5,5)=(-a2*sqrt(a5)+a1)/(a2*a5)
            do 6500 i=1,5
               do 6500 j=1,5
 6500             genrep(isymsq+i+(j-1)*5+(igen-1)*nsymsq)=rmat(i,j)
            isymsq=isymsq+25
         endif
 6600    continue
c
c       mirror plane sigma(xz)
c
      if (cla.ne.zv) go to 7000
      igen = igen + 1
      isymsq=1
      do 6900 isym=2,nsym/2
         if (ityp(isym)(1:2).eq.'t1') then
            do 6700 i=1,9
 6700          genrep(isymsq+i+(igen-1)*nsymsq)=a0
            genrep(isymsq+1+(igen-1)*nsymsq)=-a1
            genrep(isymsq+5+(igen-1)*nsymsq)=a1
            genrep(isymsq+9+(igen-1)*nsymsq)=-a1
            isymsq=isymsq+9
         elseif (ityp(isym)(1:2).eq.'t2') then
            do 6750 i=1,9
 6750          genrep(isymsq+i+(igen-1)*nsymsq)=a0
            genrep(isymsq+1+(igen-1)*nsymsq)=-a1
            genrep(isymsq+5+(igen-1)*nsymsq)=-a1
            genrep(isymsq+9+(igen-1)*nsymsq)=a1
            isymsq=isymsq+9
         elseif (idim(isym).eq.4) then
            do 6800 i=1,16
 6800          genrep(isymsq+i+(igen-1)*nsymsq)=a0
            genrep(isymsq+1+(igen-1)*nsymsq)=-a1
            genrep(isymsq+6+(igen-1)*nsymsq)=a1
            genrep(isymsq+11+(igen-1)*nsymsq)=-a1
            genrep(isymsq+16+(igen-1)*nsymsq)=a1
            isymsq=isymsq+16
         elseif (idim(isym).eq.5) then
            do 6850 i=1,25
 6850          genrep(isymsq+i+(igen-1)*nsymsq)=a0
            genrep(isymsq+1+(igen-1)*nsymsq)=a1
            genrep(isymsq+7+(igen-1)*nsymsq)=a1
            genrep(isymsq+13+(igen-1)*nsymsq)=-a1
            genrep(isymsq+19+(igen-1)*nsymsq)=-a1
            genrep(isymsq+25+(igen-1)*nsymsq)=a1
            isymsq=isymsq+25
         endif
 6900    continue
      do 6980 isym=nsym/2+1,nsym
         idimsq=idim(isym)**2
         do 6950 i=1,idimsq
 6950       genrep(isymsq+i+(igen-1)*nsymsq)=
     1                           -genrep(isymsq+i+(2*igen-3)*nsymsq/2)
         isymsq=isymsq+idimsq
 6980    continue
 7000 return
      end

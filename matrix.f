C PP 03/09/96 fixed the argument numbers in calls to matcheck
      subroutine matreset
c  this must be the first call to the object-oriented matrix system
c  it removes all matrix names but does not releases memory
c  it should be used carefully later
c  the combination
c
c      call mmark !mark a stage in memory allocation
c  ... (various calls to the matrix system, matdef calls mainly)
c      call matreset ! remove all matrices
c      call retmark
c...
c  can be used to eliminate all matrix names and reset the memory to
c  a previous stage
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      common /mamark/mark,marks(20)
      icur=0
      mark=0
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matmark
c  this routine sets a marker in the matrix definitions
c  a subsequent call to matremark removes all matrices allocated
c  since the last matmark
c  it does not allocate memory
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      common /mamark/mark,marks(20)
      mark=mark+1
      if(mark.gt.20) call nerror(0,'matmark',
     1  'more than 20 matrix markers',mark,20)
      marks(mark)=icur
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matremark
c  this routine releases a marker in the matrix definitions
c  a mark must have been set before
c  a subsequent call to matremark removes all matrices allocated
c  since the last matmark
c  it does not release memory - that must be done using  mmark and retmark
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      common /mamark/mark,marks(20)
      if(mark.lt.1) call nerror(0,'matremark',
     1  'no preceding matrix marker',mark,20)
      icur=marks(mark)
      mark=mark-1
      end
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matdef(namei, shape,idm1,idm2)
c  thus routine defines a matrix  and reserves memory for it.
c  parameters:
c  input:
c  namei (character*(*)) the name, e.g. 'A MATRIX'. 8 characters are
c  retained only
c  shape may be
c    'v' (vector)
c    'q' (square or quadratic)
c    's' (symmetrical, stored in triangular form
c    'a' (antisymmetric in triangular form, ((a(i,j),j=1,i),i=1,n)
c    'r' rectangular
c    'd' diagonal
c   idm1 and idm2 are the number of rows and columns.
c   idm2 is significant only for rectangular matrices
c   example:
c      m=5
c      n=7
c      mx='5 x 7'
c      call matdef(mx,'r',m,n)
c
c   the matrix itself is located in common /big/ at an address which
c   the system knoews, and you may inquire (see matinfo)
c

      use memory

      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6),shapec(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat),jtype(nmat)
c    icur is the current pointer - the last occupied position
c    in the list. If a matrix is released at the end of the list, the
c    list is made shorter.
      data shapes/'v','q','s','a','r','d'/
      data shapec/'V','Q','S','A','R','D'/
c     convert the dimensions to local variables
      name=namei
      id1=idm1
      id2=idm2
      if(id1.le.0.or.id2.le.0) then
      call nerror(0,'matdef','matrix dimension less than 1',id1,id2)
      end if
      mshape=0
      do 100 i=1,6
      if(shape.eq.shapes(i).or.shape.eq.shapec(i)) then
        mshape=i
       go to 200
       end if
 100  continue
 200  if(mshape.eq.0)
     1  call nerror(1,'matdef','matrix shape not v,q,s,a,r,d',id1,id2)
      call matsear(name,jcur)
      if(jcur.ne.0) then
      call nerror(2,'matdef','name multiply defined '//'name',id1,id2)
          end if
      if(id1.eq.0) call nerror(3,'matdef','null matrix defined',
     1   id1,id2)
c  new pointer
      icur=icur+1
c  name, shape (as integer), dimensions
      names(icur)=name
c  if a rectangular matrix is a vector, define it as such
Css     if (mshape.eq.5.and.id2.eq.1) mshape=1
c  if a vector is a row vector,define it as a rectangular matrix
      if(mshape.eq.1.and.id1.eq.1.and.id2.gt.1) mshape=5
c  if a rectangular matrix is square. define it as such
      if (mshape.eq.5.and.id1.eq.id2) mshape=2
      ishape(icur)=mshape
      call matleng(mshape,id1,id2,leng)
      idim1(icur)=id1
      idim2(icur)=id2
      ileng(icur)=leng
      call getmem(leng,madr)
c  address (first element of the matrix)
      iaddr(icur)=madr
c  no submatrix attached
      isub(icur)=icur
      jtype(icur)=8     ! default type is r8
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matconn(namei,shape,idm1,idm2,maddr)
c  this subroutine is similar to matdef but it does not reserve memory
c  instead, it uses a memory area already reserved in common /big/,
c  beginning with address maddr. the other parameters are identical
c  with matdef.
c  matconn is more dangerous to use than matdef because it may
c  overwrite something.
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6),shapec(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat),jtype(nmat)
c    icur is the current pointer - the last occupied position
c    in the list. If a matrix is released at the end of the list, the
c    list is made shorter.
      data shapes/'v','q','s','a','r','d'/
      data shapec/'V','Q','S','A','R','D'/
      name=namei
      id1=idm1
      id2=idm2
      mshape=0
      do 100 i=1,6
      if(shape.eq.shapes(i).or.shape.eq.shapec(i)) then
        mshape=i
       go to 200
       end if
 100  continue
 200  if(mshape.eq.0)
     1  call nerror(1,'matconn','matrix shape not v,q,s,a,r,d',id1,id2)
      call matsear(name,jcur)
      if(jcur.ne.0) then
      call nerror(2,'matdef','name multiply defined //name',id1,id2)
          end if
      if(id1.eq.0) call nerror(3,'matconn','null matrix defined',
     1   id1,id2)
cc      if(maddr.lt.1) call nerror(4,'matconn','negative start address',
cc     1 maddr,id1)
      icur=icur+1
      names(icur)=name
      ishape(icur)=mshape
      idim1(icur)=id1
      idim2(icur)=id2
      call matleng(mshape,id1,id2,leng)
      ileng(icur)=leng
      iaddr(icur)=maddr
      jtype(icur)=8           ! default type is r8
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matfrom(namei,block,leng)
c  this routine takes the array in "block" and transfers it into the
c  matrix "name". This way, the user does not have to have the common
c  bl present or manipulate it. The length of the matrix must be the same
c  as the rarray length.
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name
          dimension block(leng)
c         common /big/bl(1)
      name=namei
      call matinfo(name,mshape,id1,id2,maddr,length)
      if(leng.ne.length) then
       call nerror(1,'matfrom','matrix lengths should agree',len,length)
      end if
      call tfer(block,bl(maddr),leng)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine mattobl(namei,block,leng)
c  this routine takes the matrix "name" and transfers it into the
c  array "blokc". This way, the user does not have to have the common
c  bl present or manipulate it. The length of the matrix must be the same
c  as the rarray length.
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name
      dimension block(leng)
c     common /big/bl(1)
      name=namei
      call matinfo(name,mshape,id1,id2,maddr,length)
      if(leng.ne.length) then
       call nerror(1,'matfrom','matrix lengths should agree',len,length)
      end if
      call tfer(bl(maddr),block,leng)
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      integer function mataddr(namei)
c  finds the object namei and returns its address. -1 is returned if the
c  object is not found.  The address that is returned is relative to the
c  beginning of the common block in units of the size of the base type
c  of the object.  A real*4 value stored at location 10 in the real*8
c  common block is indicated as being positioned at location 20, i.e. where
c     it would have been located in a real*4 common block.
      character*(*) namei
      character*8 name,names,blank
      integer jcur
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat),jtype(nmat)
      name=namei
c     first, find the object
      call matsear(name,jcur)
      if(jcur.eq.0) then
       mataddr=-1
      else
       mataddr=iaddr(jcur)
      end if
      end
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matredef(oldni,newni,shape,idm1,idm2)
c  this routine changes the definition of the matrix
c  oldni is the old name (character*(*)). it must exist already
c  newni is the new name - it may be t<he same as the old
c  for the other parameters see matdef
c  the new matrix may not occupy more space than the old one
c  the length is n*m for all, except for symmatric
c  and antisymmetric (stored in n*(n+1)/2 locations, and diagonal
c  (stored in n locations), where n=idm1, m=idm2
c
      character*(*) oldni,newni
      character*8 oldn,newn,names,blank
      character*1 shape,shapes(6),shapec(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c    icur is the current pointer - the last occupied position
c    in the list. If a matrix is released at the end of the list, the
c    list is made shorter.
      data shapes/'v','q','s','a','r','d'/
      data shapec/'V','Q','S','A','R','D'/
      oldn=oldni
      newn=newni
      id1=idm1
      id2=idm2
      mshape=0
      do 100 i=1,6
      if(shape.eq.shapes(i).or.shape.eq.shapec(i)) then
        mshape=i
       go to 200
       end if
 100  continue
 200  if(mshape.eq.0)
     1  call nerror(1,'matredef','matrix shape not v,q,s,a,r',id1,id2)
      call matchec(oldn,jcur,2,'matredef')
      if(id1.eq.0) call nerror(3,'matredef','null matrix defined',
     1   id1,id2)
      lenold=ileng(jcur)
      call matleng(mshape,id1,id2,leng)
      if(leng.gt.lenold) then
      call nerror(4,'matredef','matrix cannot be extended',id1,id2)
      end if
      names(jcur)=newn
      ishape(jcur)=mshape
      idim1(jcur)=id1
      idim2(jcur)=id2
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matleng(itype,id1,id2,leng)
c     used only internally
c     itype is 1,2,3,4,5 for Vector, Quadratic,Symmetric,Antis., Rectang
c     id2 is set to the appropriate value and the length is returned in
c     leng
      go to (100,200,300,300,400,450), itype
c vector
 100  leng=id1
      id2=1
      go to 500
c  quadratic (square)
 200  leng=id1**2
      id2=id1
      go to 500
c  symmetric or antisymmetric
 300  leng=id1*(id1+1)/2
      id2=id1
      go to 500
c  rectangular
 400  leng=id1*id2
c  diagonal
      go to 500
 450  leng=id1
      id2=id1
 500  continue
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matrem(namei)
c  this routine remove the matrix namei, and releases the memory
c  it occupies. However, if the matrix is a submatrix of an existing
c  matrix, the memory is not released, and the parent matrix is not
c  removed
c  important: only the last matrix allocated can be removed.
c  therefore, it is best to remove a temporary matrix as soon as it
c  is not needed any more
c
c  parameter (input): namei (character*(*))
c

      use memory

      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      call matchec(name,jcur,1,'matrem')
      id1=idim1(jcur)
      id2=idim2(jcur)
c Do not use matleng - it gives misleading info if the matrix has been redefined
c      call matleng(ishape(jcur),id1,id2,leng)
      call matinfo(name,kshap,idd1,idd2,maddr,leng)
      call memstat(nreq,nmark,lastadr,mtot,mxmem,ioffset)
      iendmat=iaddr(jcur)+leng-1
c     if this is not a submatrix of a matrix defined earlier
c (either it has no submatrix associated with it or it is the main
c  matrix)
c  if the matrix is not submatrix of something, isub(jcur)=jcur,
c  otherwise isub(jcur)<jcur
      if(isub(jcur).ge.jcur) then
      if(lastadr.ne.iendmat.or.icur.gt.jcur) then
        call nerror(2,'matrem','there was a memory allocation later',
     1    iendmat,lastadr)
      end if
      call retmem(1)
      else
      if(icur.gt.jcur) then
        call nerror(3,'matrem','only the last matrix can be removed',
     1     icur,jcur)
      end if
      end if
c  do not return memory if this is a submatrix of a matrix defined
c  earlier
      icur=icur-1
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matdisc(namei)
c  this routine is similar to matrem. It removes the matrix namei
c  from the list but it does not release any memory
c  if the matrix is a submatrix of an existing
c  matrix, the memory is not released, and the parent matrix is not
c  removed
c  important: only the last matrix allocated can be removed.
c  therefore, it is best to remove a temporary matrix as soon as it
c  is not needed any more
c
c  parameter (input): namei (character*(*))
c
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      call matchec(name,jcur,1,'matrem')
      id1=idim1(jcur)
      id2=idim2(jcur)
      call matleng(ishape(jcur),id1,id2,leng)
      iendmat=iaddr(jcur)+leng-1
c     if this is not a submatrix of a matrix defined earlier
c (either it has no submatrix associated with it or it is the main
c  matrix)
      if(icur.gt.jcur) then
      call nerror(3,'matrem','only the last matrix can be removed',
     1   icur,jcur)
      end if
      icur=icur-1
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matsear(namei,jcur)
c  this routine searches the matrix list and locates a name
c   which is identical with "name".  The result is returned in
c   jcur which is zero if an identical name was not found, and
c   contains the number of the matrix with a previously defined
c   identical name if there is one
c   it is mostly used internally
c
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      jcur=0
      do 200 j=1,icur
      if(name.eq.names(j)) then
        jcur=j
      end if
 200  continue
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matinfo(namei,mshape,id1,id2,maddr,length)
c  this routine return information about the matrix namei
c  parameters:
c  input:
c  namei (character*(*)) name of the matrix, 8 characters are kept
c  output:
c  mshape: this is 1,2,3,4,5 or 6;
c     1 - vector
c     2 - quadratic
c     3 - symmetric
c     4 - antisymmetric
c     5 - rectangular
c     6 - diagonal
c  id1,id2: number of rows and columns
c  maddr: address, the matrix is stored from bl(maddr) in common /big/bl
c  length: length of the array, i.e. the number of memory locations
c  occupied
c  this can be used to get hold of the matrix for operations which
c  are not included here
c
      character*(*) namei
      character*8 name,names,blank
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      call matchec(name,jcur,1,'matinfo')
      mshape=ishape(jcur)
      id1=idim1(jcur)
      id2=idim2(jcur)
      length=ileng(jcur)
      maddr=iaddr(jcur)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matsub(namei,nameoli,icl1,icl2)
c  this subroutine defines a matrix as a submatrix of another
c  only the beginning and ending columns can be redefined in
c  Fortran. For a vector, it can redefine the beginning and
c  ending rows. A symmetric matrix can be truncated to smaller
c  size but it is not possible to take out a matrix from the middle
c  no new memory assignment is done, and the two matrices have common
c  memory locations. Therefore, do not use a matrix and its submatrix as
c  two parameters in a subroutine call - it will mix up most optimizers
c  the program checks, for matrix multiplies, if a result is identical
c  with one of the factors or one of its submatrices, and will not let
c  you do such a bad thing
c  parameters:
c  input:
c  namei=new name for the submatrix (character*(*))
c  nameoli: the old matrix, of which this is a part
c    icl1= starting column of the submatrix (starting row for a vector).
c    if this is zero, then the submatrix is identical with the original
c  icl2=last columns of the submatrix
c  the parent matrix' number is kept in isub(icure) where icur is the
c  number of the submatrix
      character*(*) namei,nameoli
      character*8 name,nameold,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      nameold=nameoli
      call matchec(nameold,jcur,1,'matsub')
      icol1=icl1
      icol2=icl2
      call matsear(name,jjj)
      if(jjj.ne.0) then
      call nerror(2,'matsub','new name already in use '//name,
     1  icol1,icol1)
      end if
c  for symmetric or antisymm. matrices, only submatrices 1 to n1 are
c  permitted
      if(icol1.lt.0.or.icol1.gt.idim1(jcur)) then
      call nerror(3,'matsub','dimension icol1 out if range',icol1,
     1   icol2)
      end if
      if(icol1.eq.0) then
      icol1=1
      icol2=idim1(jcur)
      end if
c   column dimension of the parent matrix
      maxcol=idim2(jcur)
      if (ishape(jcur).eq.1) maxcol=max0(idim1(jcur),idim2(jcur))
      if(icol2.lt.icol1.or.icol2.gt.maxcol) then
      call nerror(4,'matsub','dimension icol2 out if range',icol1,
     1   icol2)
      end if
c  vector or diagonal
      if (ishape(jcur).eq.1.or.ishape(jcur).eq.6) then
      icur=icur+1
      names(icur)=name
      iaddr(icur)=iaddr(jcur)+icol1-1
      ishape(icur)=ishape(jcur)
      idim1(icur)=icol2-icol1+1
      if (ishape(jcur).eq.1) then
        idim2(icur)=1
      else
        idim2(icur)=idim1(icur)
      end if
      isub(icur)=jcur
      call matleng(ishape(icur),idim1(icur),idim2(icur),leng)
      ileng(icur)=leng
      return
      else if (ishape(jcur).eq.3.or.ishape(jcur).eq.4) then
      if(icol1.gt.1) then
         call nerror(3,'matsub',
     1   'symmetric submatrix must begin at 1',icol1,icol2)
      end if
      icur=icur+1
      names(icur)=name
      iaddr(icur)=iaddr(jcur)
      ishape(icur)=ishape(jcur)
      idim1(icur)=icol2
      idim2(icur)=icol2
      isub(icur)=jcur
      call matleng(ishape(icur),idim1(icur),idim2(icur),leng)
      ileng(icur)=leng
      return
      else
      icur=icur+1
      names(icur)=name
      iaddr(icur)=iaddr(jcur)+(icol1-1)*idim1(jcur)
c  the submatrix is of the type Rectangular in general
      ishape(icur)=5
c  it is not possible to redefine the rows of a matrix, as the
c  storage is not continuous otherwise
      idim1(icur)=idim1(jcur)
      idim2(icur)=icol2-icol1+1
c  if the two dimensions are equal, it is a symmetric matrix
      if(idim1(icur).eq.idim2(icur)) ishape(icur)=2
      isub(icur)=jcur
      call matleng(ishape(icur),idim1(icur),idim2(icur),leng)
      ileng(icur)=leng
      return
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matzero(namei)
c  zeroes out a matrix or vector
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      name=namei
      call matsear(name,jcur)
      if(jcur.eq.0) then
       call nerror(1,'matzero','no such matrix '//name,0,0)
      end if
      id1=idim1(jcur)
      id2=idim2(jcur)
      msh=ishape(jcur)
      iad=iaddr(jcur)
      len=ileng(jcur)
      call zeroit(bl(iad),len)
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matchec(namei,jcur,ierrno,routine)
c  used mostly internally
c  this routine checks if a matrix is present and calls nerror if it is
c  not. jcur is the position of the matrix. nore that matsear does the
c  same search but it does not crash the job if the name is not present
c  parameters:
c  input:
c  name=matrix name, char*8
c  ierrno: error number (will be printed)
c  routine: calling routine char*8
c  output:
c  jcur: the matrix' position
c
      character*(*) namei,routine
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      call matsear(name,jcur)
      if(jcur.eq.0) then
       call nerror(1,'matchec','no such matrix '//name,0,0)
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matscal(namei,const)
c  this routine multiplies a matrix by a constant
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      name=namei
      call matchec(name,jcur,1,'matscal')
      id1=idim1(jcur)
      id2=idim2(jcur)
      len=ileng(jcur)
      iad=iaddr(jcur)
      call dscal(len,const,bl(iad),1)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matcopy(name1i,name2i)
c  this routine copies matrix name1i to matrix name2i
c   besides copying matrices of the same shape, the following is
c   permitted:
c  quadratic to symmetrical and antisymmetric
c  (the symmetric or antisymm. part of the matrix is copied
c  symmetric or antisymmetric to quadratic (the symmetric matrix is
c  expanded to quadratic)
c  diagonal to quadratic or symmetric (it is added to the diagfonal)
c  quadratic or symmetric to diagonal (the diagonal part is copied)
c  vectors may replace diagonals in these operations, and it will
c  not protest but it is best to declare a matrix diagonal, not a
c  vector
c
      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name1,name2
      data zero,half/0.0d0,0.5d0/
c  const has no significance in copying
      name1=name1i
      name2=name2i
      const=zero
      call matmanip(name1,name2,1,const)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matadd(name1i,name2i)
c  this routine adds matrix name1i to name2i
c  the result is in name2i
c  the addition can take place between matrices of different
c  shapes, as long as it makes some sense.
c  see matcopy for the description of the various ways of combining
c  matrices
c
      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name1,name2
      data zero,half/0.0d0,0.5d0/
c  const has no significance in copying
      name1=name1i
      name2=name2i
      call matmanip(name1,name2,2,const)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matadd1(name1i,const,name2i)
c  this is identical with matadd but the first matrix name1i is
c  multiplied with a const before being added to the second.
c  the result is in the second
c  see matadd and matcopy
c
      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name1,name2
      name1=name1i
      name2=name2i
      call matmanip(name1,name2,3,const)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matmanip(name1i,name2i,ioper,const)
c   used only internally
c  this is a common routine for copying, adding a scale-adding
c  matrices, including matrices of different size
c  name1 will be copied to name2 if ioper=1
c  it will be added to name2 if ioper=2
c  const*name1 will be added to name2 if ioper=3

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name,name1,name2,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      data zero,half/0.0d0,0.5d0/
      name1=name1i
      name2=name2i
      call matchec(name1,jcur1,1,'matcopy')
      call matchec(name2,jcur2,2,'matcopy')
      msh1=ishape(jcur1)
      msh2=ishape(jcur2)
      iad1=iaddr(jcur1)
      iad2=iaddr(jcur2)
      ir1=idim1(jcur1)
      ic1=idim2(jcur1)
      ir2=idim1(jcur2)
      ic2=idim2(jcur2)
      if(ir1.ne.ir2) call nerror(1,'matmanip',
     1 'matrix rows do not match',ir1,ir2)
      if (msh1.eq.msh2) then
      len=ileng(jcur1)
      if(msh1.eq.5.and.ic1.ne.ic2) call nerror(2,'matmanip',
     1    'matrix columns do not match',ic1,ic2)
      if (ioper.eq.1) then
        call dcopy(len,bl(iad1),1,bl(iad2),1)
      else if (ioper.eq.2) then
        call add(bl(iad1),bl(iad2),len)
      else if (ioper.eq.3) then
        call add1(bl(iad1),const,bl(iad2),len)
      else
        call nerror(1,'matmanip','wrong operation',ioper,ioper)
      end if
      return
      end if
      go to (100,200,300,400,500,100),msh1
c
c linear or diagonal can be copied the Q,S,A,D;it will go to the diagonal
 100  if (msh2.eq.4.or.msh2.eq.5)
     1  call nerror(3,'matmanip','cannot copy linear to A or R',
     2 ir2,ic2)
      if(msh2.eq.1.or.msh2.eq.6) then
      if(ioper.eq.1) then
        call dcopy(len,bl(iad1),1,bl(iad2),1)
      else if(ioper.eq.2) then
        call add(bl(iad1),bl(iad2),len)
      else if (ioper.eq.3) then
        call add1(bl(iad1),const,bl(iad2),len)
      end if
      go to 600
      end if
      if (msh2.eq.2) then
       call matvtoq(bl(iad1),const,bl(iad2),ir1,ioper)
       go to 600
      end if
      if(msh2.eq.3) then
      call matvtos(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
c
c quadratic can be copied to Q (already done), V (diagonal will be
c copied, S, A and D. it is best to copy it to D and not V
 200  if(msh2.eq.1.or.msh2.eq.6) then
c  copy the diagonal of a quadratic matrix to a vector
      call matqtov(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
      if (msh2.eq.3) then
c  copy a quadratic matrix to a symmetric one
      call matqtos(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
      if (msh2.eq.4) then
c  copy a quadratic matrix to an antisymm. one
      call matqtoa(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
      if (msh2.eq.5) call nerror(6,'matmanip','cannot copy Q to R',
     1  ir2,ic2)
c
c  copy symmetric matrix - can be copied to V, Q, D
 300  if (msh2.eq.1.or.msh2.eq.6) then
      call matstov(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
      if (msh2.eq.2) then
      call matstoq(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
 400  if (msh2.eq.1.or.msh2.eq.6) then
      call nerror(7,'matmanip','copying A to D does not make sense',
     1    ir1,ic1)
      go to 600
      end if
      if (msh2.eq.2) then
      call matatoq(bl(iad1),const,bl(iad2),ir1,ioper)
      go to 600
      end if
 500  continue
c   let a rectangular matrix be copied to a quadratic one if the dimensions match
c   or vice versa
       write(*,*) 'msh1,msh2,ir1,ir2,ic1,ic2',msh1,msh2,ir1,ir2,ic1,ic2
      if((msh1.eq.3.and.msh2.eq.5.or.msh2.eq.3.and.msh1.eq.5).and.
     1   ir1.eq.ir2.and.ic1.eq.ic2) then
      if (ioper.eq.1) then
        call dcopy(len,bl(iad1),1,bl(iad2),1)
      else if (ioper.eq.2) then
        call add(bl(iad1),bl(iad2),len)
      else if (ioper.eq.3) then
        call add1(bl(iad1),const,bl(iad2),len)
      else
        call nerror(1,'matmanip','wrong operation',ioper,ioper)
      end if
      else
      call nerror(8,'matmanip','Rectangular m. can be only copied to R',
     1  ir2,ic2)
      end if
 600  continue
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matprint(namei,ifil)
c  prints matrix -  and its name as title, on unit ifile
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      name=namei
      call matchec(name,jcur,1,'matprint')
      msh=ishape(jcur)
      ir1=idim1(jcur)
      ic1=idim2(jcur)
      len=ileng(jcur)
      iad=iaddr(jcur)
      write(ifil,*) name
      go to (100,200,300,300,200,100), msh
c  print vector
 100  n5=(ir1-1)/5+1
      imax=5
      if(imax.gt.ir1) imax=ir1
      write(ifil,110) (j,j=1,imax)
      imax=10
      if(imax.gt.ir1) imax=ir1
      if(imax.ge.6) write(ifil,110) (j,j=6,10)
 110  format(9x,4(i4,11x),i4)
      i1=0
      i2=4
      do 120 ii=1,n5
      if(i2.gt.ir1-1) i2=ir1-1
      write(ifil,130) i1+1,(bl(iad+i),i=i1,i2)
      i1=i1+5
      i2=i2+5
 120  continue
 130  format(1x,i5,5f14.7)
      go to 600
c print rectangular or quadratic matrix
 200  n5=(ic1-1)/5+1
      j1=0
      j2=4
      do 260 jj=1,n5
      if(j2.gt.ic1-1) j2=ic1-1
      write(ifil,110) (j+1,j=j1,j2)
      do 250 i=1,ir1
        iad1=iad+i-1
        write(ifil,130) i,(bl(j*ir1+iad1),j=j1,j2)
 250    continue
      j1=j1+5
      j2=j2+5
 260  continue
      go to 600
c  print symmetric or antisymm. matrix
 300  imax=5
      if(imax.gt.ir1) imax=ir1
      write(ifil,110) (j,j=1,imax)
      imax=10
      if(imax.gt.ir1) imax=ir1
      if(imax.ge.6) write(ifil,110) (j,j=6,10)
      do 350 i=1,ir1
      ii=i*(i-1)/2
      iad1=iad+ii
      j5=(i-1)/5+1
      j1=0
      j2=4
      do 360 jj=1,j5
        if(j2.gt.i-1) j2=i-1
        write(ifil,130) i,(bl(iad1+j),j=j1,j2)
        j1=j1+5
        j2=j2+5
 360    continue
 350  continue
 600  return
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matvtoq(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds a diagonal matrix to a square one
c   a is diagonal, b is square, ir is the dimension, ioper is 1,2or 3
c   b(diag part)=a for ioper=1
c   b=b+a for ioper=2
c   b=b+const*a for oper=3
      dimension a(ir),b(ir,ir)
      if(ioper.eq.1) then
      do 100 i=1,ir
        b(i,i)=a(i)
 100    continue
      else if (ioper.eq.2) then
      do 200 i=1,ir
        b(i,i)=b(i,i)+a(i)
 200    continue
      else if(ioper.eq.3) then
      do 300 i=1,ir
        b(i,i)=b(i,i)+const*a(i)
 300    continue
      end if
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matvtos(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds a diagonal matrix to diagonal of a
c  symmetric one
c   a is diagonal, b is triangular, ir is the dimension,
c   ioper is 1,2or 3
c   b(diag part)=a for ioper=1
c   b=b+a for ioper=2
c   b=b+const*a for oper=3
      dimension a(ir),b(*)
      if(ioper.eq.1) then
      ii=0
      do 100 i=1,ir
        ii=ii+i
        b(ii)=a(i)
 100    continue
      else if (ioper.eq.2) then
      ii=0
      do 200 i=1,ir
        ii=ii+i
        b(ii)=b(ii)+a(i)
 200    continue
      else if(ioper.eq.3) then
      ii=0
      do 300 i=1,ir
        ii=ii+i
        b(ii)=b(ii)+const*a(i)
 300    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matqtov(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds the diagonal part of matrix a to
c  vector b. b is diagonal, a is square, ir is the dimension,
c  ioper is 1,2or 3
c   b=a(diag) for ioper=1
c   b=b+a for ioper=2
c   b=b+const*a for oper=3
      dimension b(ir),a(ir,ir)
      if(ioper.eq.1) then
      do 100 i=1,ir
        b(i)=a(i,i)
 100    continue
      else if (ioper.eq.2) then
      do 200 i=1,ir
        b(i)=b(i)+a(i,i)
 200    continue
      else if(ioper.eq.3) then
      do 300 i=1,ir
        b(i)=b(i)+const*a(i,i)
 300    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matqtos(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds a quadratic matrix to a symmetrical one
c   a is square, b is triangular, ir is the dimension, ioper is 1,2or 3
c   the square matrix is symmetrized, and its symmetric part is used
c   b=(a+at)/2 for ioper=1
c   b=b+(+at)/2a for ioper=2
c   b=b+const*(a+at)/2 for oper=3
      dimension a(ir,ir),b(*)
      data half/0.5d0/
      const2=const*half
      ii=0
      if(ioper.eq.1) then
      do 100 j=1,ir
      do 100 i=1,j
        ii=ii+1
        b(ii)=(a(i,j)+a(j,i))*half
 100    continue
      else if (ioper.eq.2) then
      do 200 j=1,ir
      do 200 i=1,j
        ii=ii+1
        b(ii)=b(ii)+(a(i,j)+a(j,i))*half
 200    continue
      else if(ioper.eq.3) then
      do 300 j=1,ir
      do 300 i=1,j
        ii=ii+1
        b(ii)=b(ii)+const2*(a(i,j)+a(j,i))
 300    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matqtoa(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds a quadratic matrix to an antisymmetric
c   a is square, b is triangular, ir is the dimension, ioper is 1,2or 3
c   the square matrix is antisymmetrized, and its antisymmetric part is
c   used. The antisymmetric matrix is stored as the lower triangle
c   row-wise, i.e. element i*(i-1)/2+j is B(i,j).
c   b=(a+at)/2 for ioper=1
c   b=b+(+at)/2a for ioper=2
c   b=b+const*(a+at)/2 for oper=3
      dimension a(ir,ir),b(*)
      data half/0.5d0/
      const2=const*half
      ii=0
      if(ioper.eq.1) then
      do 100 j=1,ir
      do 100 i=1,j
        ii=ii+1
        b(ii)=(a(i,j)-a(j,i))*half
 100    continue
      else if (ioper.eq.2) then
      do 200 j=1,ir
      do 200 i=1,j
        b(ii)=b(ii)+(a(i,j)-a(j,i))*half
 200    continue
      else if(ioper.eq.3) then
      do 300 j=1,ir
      do 300 i=1,j
        b(ii)=b(ii)+const2*(a(i,j)-a(j,i))
 300    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matstov(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds the diagonal of a symmetric matrix
c  to a diagonal matrix (a vector)
c   a is triangular, b is diagonal, ir is the dimension,
c   ioper is 1,2or 3
c   b=a(diag part) for ioper=1
c   b=b+a(diag) for ioper=2
c   b=b+const*a(diag) for oper=3
      dimension a(*),b(ir)
      if(ioper.eq.1) then
      ii=0
      do 100 i=1,ir
        ii=ii+i
        b(i)=a(ii)
 100    continue
      else if (ioper.eq.2) then
      ii=0
      do 200 i=1,ir
        ii=ii+i
        b(i)=b(i)+a(ii)
 200    continue
      else if(ioper.eq.3) then
      ii=0
      do 300 i=1,ir
        ii=ii+i
        b(i)=b(i)+const*a(ii)
 300    continue
      end if
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matstoq(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds a symmetric matrix to a quadratic one
c   b is square, a is triangular, ir is the dimension, ioper is 1,2or 3
c   b=a for ioper=1
c   b=b+a for ioper=2
c   b=b+const*a for oper=3
      dimension b(ir,ir),a(*)
      ii=0
      if(ioper.eq.1) then
      do 100 j=1,ir
      do 100 i=1,j
        ii=ii+1
        b(i,j)=a(ii)
        b(j,i)=a(ii)
 100    continue
      else if (ioper.eq.2) then
      do 250 j=1,ir
        do 200 i=1,j
          ii=ii+1
          b(i,j)=b(i,j)+a(ii)
          b(j,i)=b(j,i)+a(ii)
 200      continue
        b(j,j)=b(j,j)-a(ii)
 250    continue
      else if(ioper.eq.3) then
      do 350 j=1,ir
        do 300 i=1,j
          ii=ii+1
          b(i,j)=const*a(ii)+b(i,j)
          b(j,i)=const*a(ii)+b(j,i)
 300      continue
        b(j,j)=b(j,j)-const*a(ii)
 350    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matatoq(a,const,b,ir,ioper)
c  used internally
      implicit real*8 (a-h,o-z)
c  this routine  copies or adds an antisymmetric matrix to a quadratic
c   b is square, a is triangular, ir is the dimension, ioper is 1,2or 3
c   The antisymmetric matrix is stored as the lower triangle
c   row-wise, i.e. element i*(i-1)/2+j is B(i,j).
c   b=a for ioper=1
c   b=b+a for ioper=2
c   b=b+const*a for oper=3
      dimension b(ir,ir),a(*)
      ii=0
      if(ioper.eq.1) then
      do 100 j=1,ir
      do 100 i=1,j
        ii=ii+1
        b(i,j)=a(ii)
        b(j,i)=-a(ii)
 100    continue
      else if (ioper.eq.2) then
      do 250 j=1,ir
        do 200 i=1,j
          ii=ii+1
          b(i,j)=b(i,j)+a(ii)
          b(j,i)=b(j,i)-a(ii)
 200      continue
 250    continue
      else if(ioper.eq.3) then
      do 350 j=1,ir
        do 300 i=1,j
          ii=ii+1
          b(i,j)=const*a(ii)+b(i,j)
          b(j,i)=-const*a(ii)+b(j,i)
 300      continue
 350    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matpose(namei)
c  transposes a quadratic matrix in place

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      data one/1.0d0/
      name=namei
      call matchec(name,jcur,1,'matpose')
      msh=ishape(jcur)
c   transpose has no effect on vectors, symmetric or diagonal matrices
      if(msh.eq.1.or.msh.eq.3.or.msh.eq.6) return
c   an antisymm. matrix will be its own negative upin transpose
      if(msh.eq.4) call matscal(name,-one)
      if(msh.eq.5) call nerror(2,'matpose',
     1  'rectangular matrix cannot be transposed in situ',
     2  idim1(jcur),idim2(jcur))
      if(msh.eq.2) then
      iad=iaddr(jcur)
      id=idim1(jcur)
      call trans(bl(iad),id)
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matpose2(name1i,name2i,add)
c  this routine copies the transpose of matrix name1i to name2i
c  this is the only transpose usable for rectangular matrices
c  parameters: name1i,name2i are the matrix names
c  if add='a' or 'A', the transpose is added to the content of name2i

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name1,name2,names,blank
      character*1 shape,shapes(6),add,add1,add2
      logical ladd
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      data one/1.0d0/
      data add1,add2/'a','A'/
c  transfers in matrix 2 the transpose of matrix 1; if add=a or add=A,
c  it adds it to the matrix there
      name1=name1i
      name2=name2i
      call  matchec(name1,jcur1,1,'matpose2')
      call matchec(name2,jcur2,2,'matpose2')
      ladd=add.eq.add1.or.add.eq.add2
      msh1=ishape(jcur1)
      if(msh1.eq.1.or.msh1.eq.3.or.msh1.eq.6) then
      if(ladd) then
       call matadd(name1,name2)
            else
        call matcopy(name1,name2)
      end if
      else if (msh1.eq.2) then
      if (ladd) then
        call matpose(name2)
        call matadd(name1,name2)
            else
        call matcopy(name1,name2)
            end if
      call matpose(name2)
      else if (msh1.eq.4) then
      if(ladd) then
        call matadd1(name1,-one,name2)
            else
        call matcopy(name1,name2)
        call matscal(name2,-one)
            end if
      else if (msh1.eq.5) then
      ir1=idim1(jcur1)
      ic1=idim2(jcur1)
      ir2=idim1(jcur2)
      ic2=idim2(jcur2)
      if (ir1.ne.ic2.or.ic1.ne.ir2) then
        call nerror(1,'matpose2','dimensions do not match',ir1,ic2)
      end if
      iad1=iaddr(jcur1)
      iad2=iaddr(jcur2)
      call transpo2(bl(iad1),bl(iad2),ir1,ic1,ladd)
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine transpo2(a,b,m,n,ladd)
      implicit real*8 (a-h,o-z)
      logical ladd
      dimension a(m,n),b(n,m)
      if (ladd) then
      do 100 i=1,m
      do 100 j=1,n
        b(j,i)=b(j,i)+a(i,j)
 100  continue
      else
      do 200 i=1,m
      do 200 j=1,n
        b(j,i)=a(i,j)
 200  continue
      end if
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matmmult(name1i,name2i,name3i)
c  this is a general matrix multiplier. It multiplies matrices
c  name1i and name2i and puts the result in name3i
c  the matrices must be of the correct dimension. However, diagonal,
c  symmetric, antisymmetric matrices and vectors (the latter interpreted
c  as column vectors or row vectors according to the context) may be
c  used. This routine does not allow the transposition of matrices, or
c  the addition of the result to the previous content of name3i.
c  If these services are needed, use matmmul2
c  It uses dgemm internally. It may be replaced by MXMB on Crays.
c  See the comments in matmmul2 also  See the comments in matmmul2 also
c
      character*(*) name1i,name2i,name3i
      character*8 name1,name2,name3
      character*1 no
      data no/'n'/
      name1=name1i
      name2=name2i
      name3=name3i
      call matmmul2(name1,name2,name3,no,no,no)
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matmmul2(name1i,name2i,name3i,tr1,tr2,add)
c  multiplies 1 and 2 to give 3
c  tr1 and tr2 may be 't' (or 'T': then the matrix will be first
c  temporarily transposed, otherwise not.
c  add may be 'a' or 'A': then the result will be added to 3
c  the matrices must be densely packed
c  depending on the dimensions, a  vector is interpreted as a
c  n x 1 matrix, a 1 x n matrix, or a diagonal matrix
c  the special case for symmetric results was removed - dgemm is
c  faster even if it calculates the whole matrix
c
cc if the result is symmetric or antisymmetric, only half of the
cc result is calculated and it is not checked if the matrix is
cc really symmetric. this is the responsibility of the user
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i,name3i
      character*8 name1,name2,name3
      character*8 names,blank,temp1,temp2,temp3
      character*1 shape,shapes(6),tr1,tr2,trp1,trp2,nn,nnn,tip1,tip2
      character*1 add,add1,add2
      logical ltp1,ltp2,ladd
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(20000)
      data zero,one/0.0d0,1.0d0/
      data  trp1,trp2,nn,nnn,add1,add2/'T','t','n','N','a','A'/
      name1=name1i
      name2=name2i
      name3=name3i
      call matchec(name1,jcur1,1,'matmmul2')
      call matchec(name2,jcur2,2,'matmmul2')
      call matchec(name3,jcur3,3,'matmmul2')
      msh1=ishape(jcur1)
      msh2=ishape(jcur2)
      msh3=ishape(jcur3)
      iad1=iaddr(jcur1)
      iad2=iaddr(jcur2)
      iad3=iaddr(jcur3)
      ir1=idim1(jcur1)
      ic1=idim2(jcur1)
      ir2=idim1(jcur2)
      ic2=idim2(jcur2)
      ir3=idim1(jcur3)
      ic3=idim2(jcur3)
      if((name3.eq.name1.or.name3.eq.name2).and.msh1.ne.6.and.msh2.ne.6
     1 .or.(msh3.eq.6.and.msh1.ne.6.and.msh2.ne.6)
     2 .or.msh3.eq.3.or.msh3.eq.4) then
c      define temporary result matrix (for symmetric results as well)
         temp3='MxQ27XZb'
         call matdef(temp3,'r',ir3,ic3)
         call matcopy(name3,temp3)
         call matinfo(temp3,msh3,id1,id2,iadr3,length)
      else
         temp3=name3
      end if
c transposes
      ltp1=.false.
      ltp2=.false.
      if (tr1.eq.trp1.or.tr1.eq.trp2) ltp1=.true.
      if (tr2.eq.trp1.or.tr2.eq.trp2) ltp2=.true.
c  add option: if yes, it adds the product to C, otherwise replaces
      if (add.eq.add1.or.add.eq.add2) then
         xadd=one
         ladd=.true.
      else
         xadd=zero
         ladd=.false.
c  according to the dgemm description, the matrix C does not
c  have to be set if xadd=0.0
      end if
c check dimensions
      if(ltp1) then
         irr1=ic1
         icc1=ir1
         tip1=trp1
      else
         irr1=ir1
         icc1=ic1
         tip1=nn
      end if
      if(ltp2) then
         irr2=ic2
         icc2=ir2
         tip2=trp1
      else
         irr2=ir2
         icc2=ic2
         tip2=nn
      end if
      if (icc1.ne.irr2) then
         call nerror(1,'matmmul2','columns of a .ne. rows of b',
     1               icc1,irr2)
      end if
      if (irr1.ne.ir3) then
         call nerror(2,'matmmul2','rows of a .ne. rows of c',
     1               irr1,ir3)
      end if
      if(icc2.ne.ic3) then
         call nerror(3,'matmmul2','columns of b .ne.columns of c',
     1               icc2,ic3)
      end if
c  expand symmeric matrices temporarily
      if(msh1.eq.3.or.msh1.eq.4) then
         temp1='QZ24tmpX'
         call matdef(temp1,'q',ir1,ir1)
         call matcopy(name1,temp1)
      else
         temp1=name1
      end if
      if(msh2.eq.3.or.msh2.eq.4) then
         temp2='QZ24tmpY'
         call matdef(temp2,'q',ir2,ir2)
         call matcopy(name2,temp2)
      else
         temp2=name2
      end if
      call matinfo(temp1,mxx,ixx1,ixx2,iadr1,lxx1)
      call matinfo(temp2,mxx,ixx1,ixx2,iadr2,lxx1)
      call matinfo(temp3,mxx,ixx1,ixx2,iadr3,lxx1)
c   ordinary matrix product - no diagonal matrices
c
      if (msh1.ne.6.and.msh2.ne.6) then
         call dgemm(tip1,tip2,irr1,icc2,icc1,one,bl(iadr1),ir1,
     1              bl(iadr2),ir2,xadd,bl(iadr3),ir3)
         go to 1000
      end if
c     multiply a matrix from the left by a diagonal
      if (msh1.eq.6.and.msh2.ne.6.and.(msh3.ne.3.and.msh3.ne.4)) then
         if (ltp2) then
            call matmuld(2,bl(iadr1),bl(iadr2),bl(iadr3),ir2,ic2,ic2,
     &                   ir2,ladd)
         else
            call matmuld(1,bl(iadr1),bl(iadr2),bl(iadr3),ir2,ic2,ir2,
     &                   ic2,ladd)
         end if
         go to 1000
      end if
      if(msh2.eq.6.and.msh1.ne.6.and.(msh3.ne.3.and.msh3.ne.4)) then
         if(ltp1) then
            call matmuld(4,bl(iadr2),bl(iadr1),bl(iadr3),ir1,ic1,ic1,
     &                   ir1,ladd)
         else
            call matmuld(3,bl(iadr2),bl(iadr1),bl(iadr3),ir1,ic1,ir1,
     &                   ic1,ladd)
         end if
         go to 1000
      end if
      if(msh1.eq.6.and.msh2.eq.6.and.msh3.eq.6) then
         call matmuld2(bl(iadr1),bl(iadr2),bl(iadr3),ir1,ladd)
         go to 1000
      end if
      if(msh1.eq.6.and.msh2.eq.6.and.msh3.eq.2) then
         call matmuld3(bl(iadr1),bl(iadr2),bl(iadr3),ir1,ladd)
         go to 1000
      end if
c     result symmetrical
c     if (msh1.ne.6.and.msh2.ne.6.and.(msh3.eq.3.or.msh3.eq.4)) then
c     if(irr1.ne.icc2) then
c      call nerror(4,'matmmul2','result dimension is not symmetrical',
c    1   irr1,icc2)
c     else
c       call matmuls(ltp1,ltp2,bl(iadr1),ir1,ic1,bl(iadr2),ir2,ic2,
c    1     bl(iadr3),ladd)
c     end if
c     end if
c     if ((msh1.eq.6.or.msh2.eq.6).and.(msh3.eq.3.or.msh3.eq.4)) then
c     call nerror (5,'matmmul2',
c    1  'general*diagonal=symmetric does not make sense',
c    2  ir1,ic2)
c     end if
 1000 if(msh2.eq.3.or.msh2.eq.4) call matrem(temp2)
      if(msh1.eq.3.or.msh1.eq.4) call matrem(temp1)
      if(name3.ne.temp3) then
         call matcopy(temp3,name3)
         call matrem(temp3)
      end if
      end
c
      subroutine matmuld(mode,d,a,c,ir,ic,ir1,ic1,ladd)
c  multiplies a by the diagonal matrix d from left
      implicit real*8 (a-h,o-z)
      logical ladd
      dimension a(ir,ic),c(ir1,ic1),d(*)
c mode=1 C=D*A (+C)
c mode=2 C=D*A(t) (+C)
c mode=3 C=A*D (+C)
c mode=4 C=A(t)*D (+C)
c if ladd=true then the result is added to C
      if(ladd) then
      if (mode.eq.1) then
        do 100 j=1,ic
          do 100 i=1,ir
            c(i,j)=c(i,j)+d(i)*a(i,j)
 100      continue
      else if (mode.eq.2) then
        do 200 j=1,ic
          dd=d(j)
          do 200 i=1,ir
            c(j,i)=c(j,i)+dd*a(i,j)
 200      continue
      else if (mode.eq.3) then
        do 300 j=1,ic
          dd=d(j)
          do 300 i=1,ir
            c(i,j)=c(i,j)+dd*a(i,j)
 300      continue
      else if (mode.eq.4) then
        do 400 j=1,ic
          do 400 i=1,ir
            c(j,i)=c(j,i)+d(i)*a(i,j)
 400      continue
      end if
      else
      if (mode.eq.1) then
        do 500 j=1,ic
          do 500 i=1,ir
            c(i,j)=d(i)*a(i,j)
 500      continue
      else if (mode.eq.2) then
        do 600 j=1,ic
          dd=d(j)
          do 600 i=1,ir
            c(j,i)=dd*a(i,j)
 600      continue
      else if (mode.eq.3) then
        do 700 j=1,ic
          dd=d(j)
          do 700 i=1,ir
            c(i,j)=dd*a(i,j)
 700      continue
      else if (mode.eq.4) then
        do 800 j=1,ic
          do 800 i=1,ir
            c(j,i)=d(i)*a(i,j)
 800      continue
      end if
      end if
      end
c
      subroutine matmuld2(d1,d2,d3,ir,ladd)
c  multiplies two diagonal matrices
      implicit real*8 (a-h,o-z)
      logical ladd
      dimension d1(ir),d2(ir),d3(ir)
      if (ladd) then
       do 100 i=1,ir
       d3(i)=d3(i)+d1(i)*d2(i)
 100   continue
      else
       do 200 i=1,ir
       d3(i)=d1(i)*d2(i)
 200   continue
      end if
      end
c
      subroutine matmuld3(d1,d2,a,ir,ladd)
c  multiplies two diagonal matrices into a square one
      implicit real*8 (a-h,o-z)
      logical ladd
      dimension d1(ir),d2(ir),a(ir,ir)
      data zero /0.0d0/
      if (ladd) then
       do 100 i=1,ir
       a(i,i)=a(i,i)+d1(i)*d2(i)
 100   continue
      else
      do 250 j=1,ir
      do 200 i=1,ir
         a(i,j)=zero
 200    continue
      a(j,j)=d1(j)*d2(j)
 250    continue
      end if
      end
c
      subroutine matmuls(ltp1,ltp2,a,ir1,ic1,b,ir2,ic2,
     1  c,ladd)
c symmetric or antisymmetric matrix out of the product of a  and b
      implicit real*8 (a-h,o-z)
      logical ltp1,ltp2,ladd
      dimension a(ir1,ic1),b(ir2,ic2),c(*)
      data zero /0.0d0/
      idm=ir1
      if(ltp1) idm=ic1
      if (.not.ladd) then
      call zeroit(c,idm*(idm+1)/2)
      end if
      if(.not.ltp1.and..not.ltp2) then
       do 100 k=1,ic1
         ij=0
         do 100 j=1,ic2
           bb=b(k,j)
           do 100 i=1,j
c for antisymmetric matrices, this gives the wrong sign
             ij=ij+1
             c(ij)=bb*a(i,k)+c(ij)
 100     continue
       else if (ltp1.and..not.ltp2) then
       ij=0
       do 200 i=1,ic1
         do 200 j=1,i
        ij=ij+1
        s=c(ij)
           do 250 k=1,ir2
             s=s+a(k,i)*b(k,i)
 250         continue
           c(ij)=s
 200     continue
       else if (.not.ltp1.and.ltp2) then
       do 300 k=1,ic1
         ij=0
        do 300 i=1,ir1
        do 300 j=1,i
          ij=ij+1
          c(ij)=c(ij)+a(i,k)*b(j,k)
 300     continue
       else if (ltp1.and.ltp2) then
      do 400 k=1,ir1
        ij=0
        do 400 i=1,ic1
       do 400 j=1,i
          ij=ij+1
           c(ij)=c(ij)+a(k,i)*b(j,k)
 400    continue
       end if
       end
c
      subroutine mattrace(namei,trace)
c  this subroutine calculates the trace of a matrix.
c  the matrix must be square, symmetric, or diagonal
c  the trace of an antisymmetric matrix is zero
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(5)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      name=namei
      call matinfo(name,msh,id1,id2,iadr,leng)
      if(id1.ne.id2) then
      call nerror(1,'mattrace','matrix is not symmetrical',
     1 id1,id2)
      end if
      call mattrin(bl(iadr),id1,msh,trace)
      end
c
      subroutine mattrin(a,id,ishape,trace)
c  this routine actually calculates the trace of array a
      implicit real*8 (a-h,o-z)
      dimension a(*)
      data zero/0.0d0/
      trace=zero
      if(ishape.eq.6) then
      do 100 i=1,id
        trace=trace+a(i)
 100    continue
      else if (ishape.eq.3.or.ishape.eq.4) then
      ii=0
      do  200 i=1,id
        ii=ii+i
        trace=trace+a(ii)
 200    continue
      else
      ii=-id
      do 300 i=1,id
        ii=ii+id+1
        trace=trace+a(ii)
 300    continue
      endif
      end
c
      subroutine matinv(namei)
c  this routine inverts a square or symmetric or antisymmetric or
c  diagonal matrix. For non-diagonal matrices, it calls OSINV
c  The latter uses, as far as I can see, Gauss-Jordan elimination
c  with pivoting. OSINV is in service.f . I am unsure of its
c  heritage. I hads it for ages (since about 1968). Apologies to the
c  original author.
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank,temp1,diag,diag1
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      data zero,one,tol/0.0d0,1.0d0,1.0d-12/
      name=namei
      call matchec(name,jcur,1,'matinv')
      msh=ishape(jcur)
      ir1=idim1(jcur)
      ic1=idim2(jcur)
      iadr=iaddr(jcur)
      if(msh.ne.2.and.msh.ne.3.and.msh.ne.4.and.msh.ne.6) then
      call nerror(2,'matinv',
     1  'matrix must be quadratic or (anti)symmetric',msh,ir1)
      end if
      if(msh.eq.6) then
       call matdinv(bl(iadr),ir1)
      return
      end if
      if(msh.eq.3.or.msh.eq.4) then
c  expand symmeric matrices temporarily
       temp1='QZ24tmpX'
       call matdef(temp1,'q',ir1,ir1)
       call matcopy(name,temp1)
      else
       temp1=name
      end if
      diag='QZ37tmpY'
      call matdef(diag,'d',ir1,ir1)
      call matcopy(temp1,diag)
      call matinfo(diag,msh1,id1,id2,iadr,leng)
      call matdinv(bl(iadr),id1)
      call matmmult(temp1,diag,temp1)
      call matinfo(temp1,msh1,id1,id2,iadr,leng)
      call getmem(ir1,iad2)
      call getmem(ir1,iad3)
c   matrix diagonal first normalized to one
      call osinv(bl(iadr),ir1,det,tol,bl(iad2),bl(iad3))
      call matmmult(diag,temp1,temp1)
      call matcopy(temp1,name)
      call retmem(2)
      call matrem(diag)
      if (msh.eq.3.or.msh.eq.4) then
      call matrem(temp1)
      end if
      if(det.eq.zero)  then
      call nerror(2,'matinv','matrix is singular',ir1,ic1)
      end if
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matdinv(a,n)
c  inverts a diagonal matrix
      implicit real*8 (a-h,o-z)
      dimension a(n)
      data one/1.0d0/
      do 100 i=1,n
      a(i)=one/a(i)
 100  continue
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matdiag(namei,diagi,eigeni)
c  this routine diagonalizes a symmetric matrix
c  The input matrix must be defined as symmetric.
c  It uses SDIAG2, (see in service.f). Householder transformation
c  followed by QD iteration
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei,eigeni,diagi
      character*8 name,eigen,diag,temp
      character*8 names
      character*1 shape,shapes(5)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      name=namei
      eigen=eigeni
      diag=diagi
      call matchec(name,jcur1,1,'matdiag')
      call matchec(diag,jcur2,2,'matdiag')
      call matchec(eigen,jcur3,3,'matdiag')
      msh1=ishape(jcur1)
      ir1=idim1(jcur1)
      msh2=ishape(jcur2)
      ir2=idim1(jcur2)
      msh3=ishape(jcur3)
      ir3=idim1(jcur3)
      if (msh1.ne.3.and.msh1.ne.2) then
      call nerror(4,'matdiag',
     1  'this can diagonalize only a symmetric matrix',msh1,ir1)
      end if
      if (msh3.ne.2) then
      call nerror(5,'matdiag',
     1  'eigenvector matrix must be quadratic',msh2,ir2)
      end if
      if (msh2.ne.6) then
      call nerror(6,'matdiag',
     1  'eigenvalue matrix must be diagonal',msh3,ir3)
      end if
      if (ir1.ne.ir2.or.ir1.ne.ir3) then
      call nerror(7,'matdiag','dimension mismatch',ir1,ir2)
      end if
c  copy the symmetrical matrix to a a quadratic
      if (msh1.eq.3) then
      temp='MX97sqWH'
      call matdef(temp,'q',ir1,ir1)
      call matcopy(name,temp)
      else
      temp=name
      end if
      call matinfo(temp,mshape,id1,id2,iad1,length)
      call matinfo(diag,mshape,id1,id2,iad2,length)
      call matinfo(eigen,mshape,id1,id2,iad3,length)
      call sdiag2(ir1,ir1,bl(iad1),bl(iad2),bl(iad3))
c  remove temporary quadratic matrix
      if (temp.ne.name) then
      call matrem(temp)
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matsimtr(ai,ui,bi)
c  this routine performs a similarity transformation,
c  B=U(T)AU. The arguments are the symbolic names of the matrices.
c  they are A=ai, B=bi, U=ui
c  The matrices A and B may be symmetric or square. The dimensions must
c  match properly.
c
      character*(*) ai,ui,bi
      character*8 a,u,b,temp
      a=ai
      u=ui
      b=bi
      call matinfo(a,msha,ia1,ia2,maddr,length)
      if(ia1.ne.ia2) then
      call nerror(1,'matsimtr',
     1   'only a square (or symmetric) matrix may be transformed',
     2   id1,id2)
      end if
      call matinfo(u,mshu,iu1,iu2,mad,le)
      temp='XYZqg975'
      call matdef(temp,'r',iu2,ia1)
      call matmmul2(u,a,temp,'t','n','n')
      call matmmult(temp,u,b)
      call matrem(temp)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matread(name1i,ifile,icod)
c this subroutine is an interface to the old read-rite routines
c name1i is the matrix name(input)
c ifile is the INTERNAL file number (e.g. 1 for the general file)
c icod is the code number

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,icod
      character*8 name1,names,icod1
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      common/founderr/nerrf
      name1=name1i
      icod1=icod
      call matchec(name1,jcur1,1,'matread')
      if(nerrf.lt.0) return
      msh1=ishape(jcur1)
      iad1=iaddr(jcur1)
      length=ileng(jcur1)
      call rea(bl(iad1),length,ifile,icod1)
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine matwrite(name1i,ifile,iwhere,icod)
c this subroutine is an interface to the old read-rite routines
c name1i is the matrix name(input)
c ifile is the INTERNAL file number (e.g. 1 for the general file)
c iwhere defines the position of the write:
c   iwhere=0: append
c   iwhere=k: write k-th record (and destroy all records i>=k)
c   iwhere=-k: overwrite the k-th record from the end
c icod is the code number

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,icod
      character*8 name1,names,icod1
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(1)
      common/founderr/nerrf
      name1=name1i
      icod1=icod
      call matchec(name1,jcur1,1,'matwrite')
      if(nerrf.lt.0) return
      msh1=ishape(jcur1)
      iad1=iaddr(jcur1)
      length=ileng(jcur1)
      call wri(bl(iad1),length,ifile,iwhere,icod1)
c      call wri(bl(iad1),length,ifile,icod)

      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matelem(namei,i,j,x)
c     returns the (i,j) element (or the i-th element for vectors)
c     of the matrix in x

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      name=namei
      call matchec(name,jcur,1,'matelem')
      msh=ishape(jcur)
      ir1=idim1(jcur)
      ic1=idim2(jcur)
      len=ileng(jcur)
      iad=iaddr(jcur)-1
c             V   Q   S   A   R   D
      go to (100,200,300,300,200,100), msh
c  vector or diagonal matrix
 100  x=bl(iad+i)
      go to 600
c rectangular or quadratic matrix
 200  x=bl(iad+i+j*ir1-ir1)
      go to 600
c  symmetric or antisymm. matrix
 300  if (j.le.i) then
      x=bl(iad+i*(i-1)/2+j)
      else
      x=bl(iad+j*(j-1)/2+i)
      if(msh.eq.4) x=-x
      end if
 600  continue
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine mateset(namei,i,j,x)
c     sets the (i,j) element (or the i-th element for vectors)
c     of the matrix to x

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(100)
      name=namei
      call matchec(name,jcur,1,'mateset')
      msh=ishape(jcur)
      ir1=idim1(jcur)
      ic1=idim2(jcur)
      len=ileng(jcur)
      iad=iaddr(jcur)-1
c             V   Q   S   A   R   D
      go to (100,200,300,300,200,100), msh
c  vector or diagonal matrix
 100  bl(iad+i)=x
      go to 600
c rectangular or quadratic matrix
 200  bl(iad+i+j*ir1-ir1)=x
      go to 600
c  symmetric or antisymm. matrix
 300  if (j.le.i) then
      bl(iad+i*(i-1)/2+j)=x
      else
      if(msh.eq.3) then
        bl(iad+j*(j-1)/2+i)=x
      else
        bl(iad+j*(j-1)/2+i)=-x
      end if
      end if
 600  continue
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine matprodtr(name1i,name2i,trace)
c  this subroutine calculates the trace of the product of two matrices.
c  the first one is transposed, i.e. Tr(A(t)B)
c  the product matrix must be square in shape.
c

      use memory

      implicit real*8 (a-h,o-z)
      character*(*) name1i,name2i
      character*8 name1,name2,names,blank,temp
      character*1 shape1,shape2,shapes(5)
      parameter (nmat=200)
      parameter(zero=0.0d0,half=0.5d0)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
c     common /big/bl(10000)
      name1=name1i
      name2=name2i
      call matinfo(name1,msh1,id1,id2,iadr1,leng1)
      call matinfo(name2,msh2,id3,id4,iadr2,leng2)
      if(id1.ne.id3.or.id2.ne.id4) then
      call nerror(1,'matprodtr',
     1 'multiplication is impossible or matrix is not symmetrical',
     2 id1,id3)
      end if
c  Tr(AtB)=Tr(BtA) so order them appropriately
c reorder the matrices so that the type of matrix 2 is higher
c  (R is the same as Q here)
c  msh= 2 (quadr), 3 (symm), 4 (antisymm) 5=(rectangular), 6=(diagonal)
      if(msh2.eq.2.or.msh2.eq.5.or.
     1        msh1.eq.6.or.(msh1.eq.4.and.msh2.eq.3)) then
      temp=name1
      name1=name2
      name2=temp
      msh3=msh1
      msh1=msh2
      msh2=msh3
      iadr3=iadr1
      iadr1=iadr2
      iadr2=iadr3
      end if
      trace=zero
      if(msh1.eq.2.or.msh1.eq.5) then
      if(msh2.eq.3.or.msh2.eq.4) then
        temp='XW7862&z'
        call matdef(temp,'q',id1,id1)
c  expand the triangular matrix to square
        call matcopy(name2,temp)
        call matinfo(temp,isht,id5,id6,iadr2,illl)
      end if
      if(msh2.lt.6) then
        trace=ddot(id1*id2,bl(iadr1),1,bl(iadr2),1)
      else
        ii=0
        do 100 i=1,id1
          trace=trace+bl(iadr1+ii)*bl(iadr2+i-1)
          ii=ii+id1+1
 100      continue
      end if
      if(msh2.eq.3.or.msh2.eq.4) call matrem(temp)
      end if
      if(msh1.eq.3.or.msh1.eq.4) then
      if(msh2.eq.3.or.msh2.eq.4) then
        ii=-1
        ij=-1
        do 300 i=1,id1
          do 200 j=1,i
            ij=ij+1
            trace=trace+bl(iadr1+ij)*bl(iadr2+ij)
 200        continue
          ii=ii+i
          trace=trace-half*bl(iadr1+ii)*bl(iadr2+ii)
 300      continue
        trace=trace+trace
        if(msh1.eq.3.and.msh2.eq.4) trace=zero
c  the antisymmetric case gives zero here
      else if(msh2.eq.6) then
        ii=-1
        do 400 i=1,id1
          ii=ii+i
          trace=trace+bl(iadr1+ii)*bl(iadr2+i-1)
 400      continue
      end if
      end if
      if(msh1.eq.6) then
c in this case, msh2 is also diagonal, otherwise they were interchanged
      do 500 i=1,id1
        trace=trace+bl(iadr1+i-1)*bl(iadr2+i-1)
 500    continue
      end if
      end
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine dynamic_matdef(namei, shape,idm1,idm2)
c  this routine defines a matrix and reserves memory for it.
c  parameters:
c  input:
c  namei (character*(*)) the name, e.g. 'A MATRIX'. 8 characters are
c  retained only
c  shape may be
c    'v' (vector)
c    'q' (square or quadratic)
c    's' (symmetrical, stored in triangular form
c    'a' (antisymmetric in triangular for, ((a(i,j),j=1,i),i=1,n)
c    'r' rectangular
c    'd' diagonal
c   idm1 and idm2 are the number of rows and columns.
c   idm2 is significant only for rectangular matrices
c   example:
c      m=5
c      n=7
c      mx='5 x 7'
c      call matdef(mx,'r',m,n)
c
c   the matrix itself is located in common /big/ at an address which
c   the system knows, and you may inquire (see matinfo)
c
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6),shapec(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat),jtype(nmat)
c    icur is the current pointer - the last occupied position
c    in the list. If a matrix is released at the end of the list, the
c    list is made shorter.
      data shapes/'v','q','s','a','r','d'/
      data shapec/'V','Q','S','A','R','D'/
c     convert the dimensions to local variables
      name=namei
      id1=idm1
      id2=idm2
      if(id1.le.0.or.id2.le.0) then
      call nerror(0,'matdef','matrix dimension less than 1',id1,id2)
      end if
      mshape=0
      do 100 i=1,6
      if(shape.eq.shapes(i).or.shape.eq.shapec(i)) then
        mshape=i
       go to 200
       end if
 100  continue
 200  if(mshape.eq.0) then
c             write(*,'(F10.x)') xx
        call nerror(1,'matdef','matrix shape not v,q,s,a,r,d',id1,id2)
      endif
      call matsear(name,jcur)
      if(jcur.ne.0) then
      call nerror(2,'matdef','name multiply defined '//'name',id1,id2)
          end if
      if(id1.eq.0) then
c             write(*,'(F10.x)') xx
              call nerror(3,'matdef','null matrix defined',
     1   id1,id2)
      endif
c  new pointer
      icur=icur+1
c  name, shape (as integer), dimensions
      names(icur)=name
c  if a rectangular matrix is a vector, define it as such
Css     if (mshape.eq.5.and.id2.eq.1) mshape=1
c  if a vector is a row vector,define it as a rectangular matrix
      if(mshape.eq.1.and.id1.eq.1.and.id2.gt.1) mshape=5
c  if a rectangular matrix is square. define it as such
      if (mshape.eq.5.and.id1.eq.id2) mshape=2
      ishape(icur)=mshape
      call matleng(mshape,id1,id2,leng)
      idim1(icur)=id1
      idim2(icur)=id2
      ileng(icur)=leng
      call dynamic_getmem(leng,madr)
c  address (first element of the matrix)
      iaddr(icur)=madr
c  no submatrix attached
      isub(icur)=icur
      jtype(icur)=8     ! default type is r8
      end
c
c++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
      subroutine dynamic_matrem(namei)
c  this routine remove the matrix namei, and releases the memory
c  it occupies. However, if the matrix is a submatrix of an existing
c  matrix, the memory is not released, and the parent matrix is not
c  removed
c  important: only the last matrix allocated can be removed.
cdynamic_  therefore, it is best to remove a temporary matrix as soon as it
c  is not needed any more
c
c  parameter (input): namei (character*(*))
c
      character*(*) namei
      character*8 name,names,blank
      character*1 shape,shapes(6)
      parameter (nmat=200)
      common/matrix/icur,names(nmat),ishape(nmat),isub(nmat),
     1 idim1(nmat),idim2(nmat),iaddr(nmat),ileng(nmat)
      name=namei
      call matchec(name,jcur,1,'matrem')
      id1=idim1(jcur)
      id2=idim2(jcur)
      call matleng(ishape(jcur),id1,id2,leng)
      iendmat=iaddr(jcur)+leng-1
c     if this is not a submatrix of a matrix defined earlier
c (either it has no submatrix associated with it or it is the main
c  matrix)
c  if the matrix is not submatrix of something, isub(jcur)=jcur,
c  otherwise isub(jcur)<jcur
      if(isub(jcur).ge.jcur) then
      call dynamic_retmem_matrix(1,iremoved)
      if (iremoved.ne.iaddr(jcur)) then
c       write(*,'(F20.x)') x
        call nerror(2,'matrem','This is not the last memory allocation',
     1    iaddr(jcur),iremoved)
      endif
      else
      if(icur.gt.jcur) then
c       write(*,'(F20.x)') x
        call nerror(3,'matrem','only the last matrix can be removed',
     1     icur,jcur)
      end if
      end if
c  do not return memory if this is a submatrix of a matrix defined
c  earlier
      icur=icur-1
      end


C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
      subroutine matmul_mkl(A,B,C,m,n,p)
c  This is a simple version of degmm from MKL for matrix multiplication
c  It is convenient to call this subroutine
c  Zhigang_9/14/2016 @ UARK

      implicit none
      integer m,n,p
      real(kind=8) A(m,n),B(n,p),C(m,p)

      call dgemm('N','N',m,p,n,1.0D0,A,m,B,n,0.0D0,C,m)

      end subroutine matmul_mkl


      subroutine FockSymm_SCF(ngener,ncf,nfock,fock,ifp)
      implicit real*8 (a-h,o-z)
c This routine builds the whole Fock matrix from the
c skeleton one
c PARAMETERS:
c  INPUT:
c    ngener= number of group generators
c    ncf=number of contracted basis functions
c    nfock=number of Fock matrices
c    fock=Fock matrices
c    ifp=function pairs, ifp(isymop,k) gives the symmetry
c       partner of basis function k under symm. op. isymop;
c       if ifp(isymop,k) is negative then the function changes
c       sign upon reflection or rotation
c OUTPUT: the full Fock matrix in fock
C
C ---------------------------------------------------------------
C  WARNING!  Modified to remove use of nfock which is assumed
C            to be unity
C ---------------------------------------------------------------
      dimension fock(*)
      dimension ifp(7,ncf)
      Parameter (one=1.0d0,half=0.5d0)
c
c .....................................................
      if(nfock.ne.1)
     $  call nerror(1,'focksymm_scf','illegal use of multifock',0,0)
c .....................................................
      do 350 ns=1,ngener
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if (i1.gt.0) go to 320
            i1=-i1
            fct=-fct
  320       continue
            j1=ifp(ns,j)
            if (j1.gt.0) go to 330
            j1=-j1
            fct=-fct
  330       continue
            ij1=i1*(i1-1)/2+j1
            if (j1.gt.i1) ij1=j1*(j1-1)/2+i1
c ..............................................
            ff=fock(ij)+fct*fock(ij1)
            if (ij.gt.ij1) ff=ff*half
            fock(ij)=ff
            fock(ij1)=fct*ff
c ..............................................
  340    continue
  350 continue
c
      end
c====================================================================
      subroutine OrbSymm(nsym,ncf    ,iorb   ,coef,  smco,
     1                   ifp ,ptgroup,species,nspec,tmp)
c  this routine determines the symmetry characters of an orbital
c  under the Abelian point group
c
c  Most of the computational work is now done before calling this
c  and only once.
c  The symm labels for the guess orbitals are sometimes different
c  but the change is due to the change in the diagonalizer.
c
c  INPUT:
c  nsym=number of symmetry operations, max. 7
c  ncf =number of contracted functions
c  iorb=orbital to be analyzed, forms acolumn if coef(*,iorb)
c  coef=orbital coefficient matrix, first index basis fnct, second orbital
c  smco  = iorb-th column of S*C matrix (overlap*coeffs),ncf long
c  ifp(7,ncf)= basis function symmetry pairs. abs(ifp(ns,k)) gives
c  the image of basis function k under symmetry operation ns.
c  This number is negative if k changes sign under the symm. operation
c  ptgroup: the name of the point group, D2h, D2, C2v, C2h, C2 or Cs
c  OUTPUT
c  species: (char*3) the symmetry species
c  nspec: a number corresponding to "species"
c    this is :20 for ?
c              1 for A,   2 for B
c              3 for A',  4 for A"
c              5 for Ag,  6 for Au
c              7 for A1,  8 for A2
c              9 for B1, 10 for B2, 11 for B3
c             12 for Bg, 13 for Bu
c             14 for B1g, 15 for B1u, 16 for B2g, 17 for B2u, 18 for B3g, 19 for B3u
c  STORAGE:
c  tmp  = temporary storage, ncf long
c
      implicit real*8 (a-h,o-z)
      character*3 ptgroup
      character*8 species
      dimension coef(ncf,*),tmp(ncf),smco(ncf)
      dimension ifp(7,ncf),ichar(7)
      parameter(zero=0.0d0,one=1.0d0,eps=1.0d-6)
c
      nspec=0
      if(nsym.eq.0) return
      do ns=1,nsym
        do k=1,ncf
c generate the symmetry image of orbital k under symmetry op. ns
           ispr=ifp(ns,k)
           isp1=abs(ispr)
           if(ispr.gt.0) then
             tmp(isp1)=coef(k,iorb)
           else if(ispr.lt.0) then
             tmp(isp1)=-coef(k,iorb)
           else
             tmp(isp1)=zero
           end if
         end do
c  Determine the overlap between the orbital and its symmetry image
c  smco is a column from product of the overlap and coefficients matrix
        ovrl= ddot(ncf,tmp,1,smco,1)
c       ovrl=zero
c       do k=1,ncf
c         sum=zero
c         kl=k*(k-1)/2
c         do l=1,ncf
c           if(k.ge.l) then
c             kl=kl+1
c           else
c             kl=kl+l-1
c           end if
c           sum=sum+sm(kl)*tmp(l)
c         end do
c         ovrl=ovrl+sum*coef(k,iorb)
c       end do
        if(abs(ovrl-one).lt.eps) then
          ichar(ns)=1
        else if(abs(ovrl+one).lt.eps) then
          ichar(ns)=-1
        else
          ichar(ns)=0
        end if
        if(ichar(ns).eq.0) then
          species=' ? '
          nspec=20
        end if
c     print *,ns,iorb,ovrl,ichar(ns)
      end do
      if(nsym.eq.1) then
        if(ptgroup.eq.'Cs '.or.ptgroup.eq.'cs ') then
          if(ichar(1).eq.1) then
            species='A'' '
            nspec=3
          else if(ichar(1).eq.-1) then
            species='A" '
            nspec=4
          end if
        else if(ptgroup.eq.'C2 '.or.ptgroup.eq.'c2 ') then
          if(ichar(1).eq.1) then
            species='A  '
            nspec=1
          else if(ichar(1).eq.-1) then
            species='B  '
            nspec=2
          end if
        else if(ptgroup.eq.'Ci '.or.ptgroup.eq.'ci ') then
          if(ichar(1).eq.1) then
            species='Ag '
            nspec=5
          else if(ichar(1).eq.-1) then
            species='Au '
            nspec=6
          end if
        end if
      end if
      if(nsym.eq.3) then
        if(ptgroup.eq.'C2v'.or.ptgroup.eq.'c2v') then
          if(ichar(1).eq.1.and.ichar(3).eq.1) then
            species='A1 '
            nspec=7
          else if(ichar(3).eq.1.and.ichar(1).eq.-1) then
            species='A2 '
            nspec=8
          else if(ichar(3).eq.-1.and.ichar(1).eq.1) then
            species='B1 '
            nspec=9
          else if(ichar(3).eq.-1.and.ichar(1).eq.-1) then
            species='B2 '
            nspec=10
          end if
        else if(ptgroup.eq.'D2 '.or.ptgroup.eq.'d2 ') then
c          print *, 'ptgroup,ichar  ',ptgroup,(ichar(kk),kk=1,3)
          if(ichar(1).eq.1.and.ichar(3).eq.1) then
            species='A  '
            nspec=1
          else if(ichar(1).eq.1.and.ichar(3).eq.-1) then
            species='B1 '
            nspec=9
          else if(ichar(1).eq.-1.and.ichar(3).eq.-1) then
            species='B2 '
            nspec=10
          else if(ichar(1).eq.-1.and.ichar(3).eq.1) then
            species='B3'
            nspec=11
          end if
           else if(ptgroup.eq.'C2h'.or.ptgroup.eq.'c2h') then
          if(ichar(2).eq.1.and.ichar(1).eq.1) then
            species='Ag '
            nspec=5
          else if(ichar(2).eq.-1.and.ichar(1).eq.-1)then
            species='Bg '
            nspec=12
          else if(ichar(2).eq.1.and.ichar(1).eq.-1) then
            species='Au'
            nspec=6
          else if(ichar(2).eq.-1.and.ichar(1).eq.1) then
            species='Bu '
            nspec=13
          end if
        end if
      end if
c  D2h group
      if(nsym.eq.7) then
        if(ichar(4).eq.1) then
          if(ichar(5).eq.1) then
            if(ichar(7).eq.1) then
              species='Ag '
              nspec=5
            else if(ichar(7).eq.-1) then
              species='Au '
              nspec=6
            end if
          else if(ichar(5).eq.-1) then
            if(ichar(7).eq.1) then
              species='B1g'
              nspec=14
            else if(ichar(7).eq.-1) then
              species='B1u'
              nspec=15
            end if
          end if
        else if(ichar(4).eq.-1) then
          if(ichar(5).eq.1) then
            if(ichar(7).eq.1) then
              species='B2g'
              nspec=16
            else if(ichar(7).eq.-1) then
              species='B2u'
              nspec=17
            end if
          else if(ichar(5).eq.-1) then
            if(ichar(7).eq.1) then
              species='B3g'
              nspec=18
            else if(ichar(7).eq.-1) then
              species='B3u'
              nspec=19
            end if
          end if
        end if
      end if
c      print *, 'nspec, species',nspec,species
      end
c=======================================================
      subroutine DensSymm(nsym,ncf,dens,ifp,devmx)
      implicit real*8 (a-h,o-z)
c This routine symmetrizes the density matrix and detects
c  the maximum deviation from symmetry
c PARAMETERS:
c  INPUT:
c    nsym= number of symmetry operations (max. 7)i
c    It should be just the generators of the group
c    ncf=number of contracted basis functions
c    dens=density matrix
c    ifp=function pairs, ifp(isymop,k) gives the symmetry
c       parter of basis function k under symm. op. isymop;
c       if ifp(isymop,k) is negative then the function changes
c       sign upon reflection or rotation
c   OUTPUT:
c     devmx=max. deviation from symmetry
c      and the symmetrized density in dens
      dimension dens(*)
      dimension ifp(7,ncf)
      data zero,one,half /0.0d0, 1.d0 , 0.5d0/
c
      devmx=zero
      if(nsym.eq.0) then
        RETURN
      end if
      do 350 ns=1,nsym
         ij=0
         do 340 i=1,ncf
         do 340 j=1,i
            fct=one
            ij=ij+1
            i1=ifp(ns,i)
            if (i1.lt.0) then
              i1=-i1
              fct=-fct
            end if
            j1=ifp(ns,j)
            if (j1.lt.0) then
              j1=-j1
              fct=-fct
            end if
            if(i1.ge.j1) then
              ij1=i1*(i1-1)/2+j1
            else
              ij1=j1*(j1-1)/2+i1
            end if
            ff=(dens(ij)+fct*dens(ij1))*half
                delta=abs(dens(ij)-ff)
            dens(ij)=ff
            dens(ij1)=fct*ff
            if(delta.gt.devmx) devmx=delta
  340    continue
  350 continue
c
      end
c====================================================================
      subroutine GradSymm(ngen,nsy,ncf,ifp,icoord,G)
      implicit real*8 (a-h,o-z)
c
c  This routine builds the partial DFT gradient matrix from the
c  skeleton one:
c      DXi.V.Xj        Xi,Xj are basis functions
c                      D is a derivative (x,y or z)
c                      V is the (DFT) potential
c
c  ARGUMENTS
c
c  ngen    -  number of generators
c  nsy     -  list of symmetry operations (generators)
c  ncf     -  number of basis functions
c  ifp     -  list of basis function pairs under symmetry operations
c             ifp(isymop,k) gives the symmetry partner of basis
c             function k under symmetry operation isymop;
c             if ifp(isymop,k) is negative, then the function
c             changes sign upon reflection or rotation
c  icoord  -  derivative type  1 - X; 2 - Y; 3 - Z
c  G       -  input: skeleton partial DFT gradient matrix
c             exit:  symmetrized partial DFT gradient matrix
c
      dimension G(ncf,ncf)
      dimension nsy(ngen),ifp(7,ncf)
      dimension mirror(3,7)
      parameter (Half=0.5d0)
      data mirror(1,1),mirror(2,1),mirror(3,1),mirror(1,2),mirror(2,2),m
     1irror(3,2),mirror(1,3),mirror(2,3),mirror(3,3),mirror(1,4),mirror(
     22,4),mirror(3,4),mirror(1,5),mirror(2,5),mirror(3,5),mirror(1,6),m
     3irror(2,6),mirror(3,6),mirror(1,7),mirror(2,7),mirror(3,7)/-1,1,1,
     41,-1,1,1,1,-1,-1,-1,1,-1,1,-1,1,-1,-1,-1,-1,-1/
c
      do 350 ns=1,ngen
         isyop=nsy(ns)
         do 340 i=1,ncf
         do 340 j=1,ncf
            fct=mirror(icoord,isyop)
            i1=ifp(ns,i)
            if (i1.gt.0) go to 320
            i1=-i1
            fct=-fct
  320       continue
            j1=ifp(ns,j)
            if (j1.gt.0) go to 330
            j1=-j1
            fct=-fct
  330       continue
            ff=G(i,j)+fct*G(i1,j1)
            if(i.gt.i1.or.(i.eq.i1.and.j.gt.j1)) ff = Half*ff
            G(i,j) = ff
            G(i1,j1) = fct*ff
  340    continue
  350 continue
c
      return
      end
c===========================================================================
      subroutine PrintOrbSym(nsym,ncf,nmo,nvir1,coef,sm,ifp,printew,ew)
c  This routine calls OrbSymm and prints the orbital symmetry
c Arguments:
c INPUT:
c  nsym=number of symmetry operations
c  ncf=number of contracted basis functions
c  nmo=number of occupied orbitals to be tested
c  nvirt=number of virtuals to be tested
c  coef: SCF coefficients, coef(ncf,norbs)
c  sm: overlap matrix, symmetrical in triangular form
c  ifp: symmetry images of the basis functions under the 7 symmetry operations
c  ifp(k,iorb) is the orbital iorb is transformed under operation k;
c  if it is negative, iorb is transformed into the negative of orb. ifp(k,iorb)
c  printew= if true, eigenvalues are printed
c  ew: orbital energies
c

      use memory

      implicit real*8 (a-h,o-z)
      logical printew,nonabel
      parameter (line=5,nspec=20)
      dimension coef(ncf,*),sm(*),ew(ncf),ifp(7,ncf)
      character abel*3,sflies*4,symbol*3,spec*3
      character*8 species(line)
      dimension noccup(nspec),symbol(nspec),nsum(nspec),spec(nspec)
c     common /big/bl(10000)
      data symbol
     1  /'a  ','b  ','a'' ','a" ','ag ','au ','a1 ','a2 ','b1 ','b2 ',
     2   'b3 ','bg ','bu ' ,'b1g','b1u','b2g','b2u','b3g','b3u',' ? '/
          call getmem(ncf,iadr)
c storage for S*coef
          call getmem(ncf*ncf,iadr2)
          call zeroit(bl(iadr2),ncf*ncf)
c temp. storage for unpacking S, not too elegant
          call getmem(ncf*ncf,iadr3)
          call packed_to_upper(sm,bl(iadr3),ncf)
c calculate S*coef
          call dsymm('L','U',ncf,ncf,1.0d0,bl(iadr3),ncf,coef,ncf,
     *               0.0d0,bl(iadr2),ncf)
          call retmem(1)
          call getchval('Sflies',sflies)
          call getchval('Abelian',abel)
          iout=igetival('iout')
          nonabel = sflies.NE.'c1  '.AND.sflies.NE.'ci  '.AND.
     $              sflies.NE.'cs  '.AND.sflies.NE.'c2  '.AND.
     $              sflies.NE.'c2v '.AND.sflies.NE.'c2h '.AND.
     $              sflies.NE.'d2  '.AND.sflies.NE.'d2h '
c  store nvirt in a new variable
      nvirt=nvir1
c  threshold for degeneracy
      del=1.0d-5
c zero out the occupancy matrix array
      call izeroit(noccup,nspec)
c Print first the occupied ones
      nlines=nmo/line
      nrest=nmo-nlines*line
c -------------------------------------------------------------
          if(printew) then
           write(iout,*) '  orbital symmetries and energies'
          else
           write(iout,*) '  orbital symmetries'
          endif
          write(iout,*) '  occupied orbitals'
c -------------------------------------------------------------
      iorb=0
c  this is the loop over the full lines
      do k=1,nlines
        do i=1,line
          iorb=iorb+1
          if(nsym.gt.0) then
        call OrbSymm(nsym,ncf ,iorb      ,coef,bl(iadr2+ncf*(iorb-1)),
     1               ifp ,abel,species(i),nsp ,bl(iadr))
c          print *, 'iorb,nsp',iorb,nsp
            if(nsp.gt.0) noccup(nsp)=noccup(nsp)+1
            if(printew) then
              iorbp=iorb+1
              iorbm=iorb-1
              if(iorbm.lt.1) iorbm=iorbp
c -- check for possible degeneracies if non-abelian
                  IF(nonabel) THEN
              if(abs(ew(iorb)-ew(iorbp)).lt.del.or.
     1           abs(ew(iorb)-ew(iorbm)).lt.del) species(i)(4:6)='(e)'
                  ENDIF
            end if
          end if
        end do
        if(nsym.gt.0) write(iout,100) iorb+1-line,species
  100   format(i4,3x,(5(a8,5x)))
        if(printew) write(iout,200) (ew(kk),kk=iorb-line+1,iorb)
  200   format(5f13.5)
      end do
c  now loop over the last line - incompletely filled
      do i=1,nrest
        iorb=iorb+1
        if(nsym.gt.0) then
        call OrbSymm(nsym,ncf ,iorb      ,coef,bl(iadr2+ncf*(iorb-1)),
     1               ifp ,abel,species(i),nsp ,bl(iadr))
            if(nsp.gt.0) noccup(nsp)=noccup(nsp)+1
            if(printew) then
              iorbp=iorb+1
              iorbm=iorb-1
              if(iorbm.lt.1) iorbm=iorbp
c -- check for possible degeneracies if non-abelian
                  IF(nonabel) THEN
              if(abs(ew(iorb)-ew(iorbp)).lt.del.or.
     1           abs(ew(iorb)-ew(iorbm)).lt.del) species(i)(4:6)='(e)'
                  ENDIF
            end if
        end if
      end do
      if(nrest.gt.0.and.nsym.gt.0) then
        write(iout,100) iorb+1-nrest,(species(j),j=1,nrest)
      end if
      if(printew) write(iout,200) (ew(k),k=iorb-nrest+1,iorb)
c now do the virtuals  do not add the occupancy to noccup
      if(nmo+nvirt.gt.ncf) nvirt=ncf-nmo
c  in there are no virtuals, just print summary occupancy
      if(nvirt.le.0) go to 300
c
      nlines=nvirt/line
      nrest=nvirt-nlines*line
      write(iout,*) '  virtual orbitals'
      iorb=nmo
      do k=1,nlines
        do i=1,line
          iorb=iorb+1
          if(nsym.gt.0) then
        call OrbSymm(nsym,ncf ,iorb      ,coef,bl(iadr2+ncf*(iorb-1)),
     1               ifp ,abel,species(i),nsp ,bl(iadr))
            if(printew) then
              iorbp=iorb+1
              iorbm=iorb-1
              if(iorbp.gt.ncf) iorbp=iorbm
c -- check for possible degeneracies if non-abelian
                  IF(nonabel) THEN
              if(abs(ew(iorb)-ew(iorbp)).lt.del.or.
     1           abs(ew(iorb)-ew(iorbm)).lt.del) species(i)(4:6)='(e)'
                  ENDIF
            end if
          end if
        end do
       if(nsym.gt.0) write(iout,100) iorb+1-line,species
        if(printew) write(iout,200) (ew(kk),kk=iorb-line+1,iorb)
      end do
      do i=1,nrest
        iorb=iorb+1
        if(nsym.gt.0) then
        call OrbSymm(nsym,ncf ,iorb      ,coef,bl(iadr2+ncf*(iorb-1)),
     1               ifp ,abel,species(i),nsp ,bl(iadr))
            if(printew) then
              iorbp=iorb+1
              iorbm=iorb-1
              if(iorbp.gt.ncf) iorbp=iorbm
c -- check for possible degeneracies if non-abelian
                  IF(nonabel) THEN
              if(abs(ew(iorb)-ew(iorbp)).lt.del.or.
     1           abs(ew(iorb)-ew(iorbm)).lt.del) species(i)(4:6)='(e)'
                  ENDIF
            end if
        end if
      end do
      if(nrest.gt.0.and.nsym.gt.0) then
        write(iout,100) iorb+1-nrest,(species(j),j=1,nrest)
      end if
      if(printew) write(iout,200) (ew(k),k=iorb-nrest+1,iorb)
  300 call retmem(2)
      if(nsym.gt.0) then
c        print *,'noccup',noccup
        ilast=0
        do ii=1,nspec
          if(noccup(ii).gt.0) then
            ilast=ilast+1
            nsum(ilast)=noccup(ii)
            spec(ilast)=symbol(ii)
          end if
        end do
        write(iout,400) (spec(k),k=1,ilast)
        write(iout,500) (nsum(k),k=1,ilast)
  400   format(' Summary occupancy',9(2x,a3,1x))
  500   format(16x,9i6)
      end if
      end

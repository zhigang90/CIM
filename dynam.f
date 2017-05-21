CCP   Original saved in Replaced_Feb
      subroutine dynamain(istep,finished,na)
c  this is just a wrapper; it only allocates memory
c  istep is the dynamics step number; needed because
c  the first two steps are different

      use memory

      implicit real*8(a-h,o-z)
c     common /big/bl(30 000)
c .............................................
c -- automatic allocation of arrays in F90
      CHARACTER*8 AtSymb(na)
c .............................................
      logical finished,exst
      save icoord,iold,inew,imass,iatchrg,
     1     iaccel,iaccold,iaolder
      character*256 jobname
c
      iout=igetival('iout')
      write(iout,*) 'Dynamics program, step=',istep
c
      call setival('isumscf',1)      ! switch off SCF summary print
      finished=.false.
c
c  read in current geometry
      call getchval('jobname',jobname)
      call rmblan2(jobname,256,lenJ)
c  open the trajectory file
      itraj=35
      if(istep.eq.0) then
        inquire(file=jobname(1:lenJ)//'.trajec',exist=exst)
        if(exst) then
          open(itraj,file=jobname(1:lenJ)//'trajec')
          close(itraj,status='delete')
        end if
        open(unit=itraj,file=jobname(1:lenJ)//'.trajec')
c
c  allocate storage
        call getmem(3*na,icoord)
        call getmem(3*na,iold)
        call getmem(3*na,inew)
        call getmem(na,imass)
        call getmem(na,iatchrg)
        call getmem(3*na,iaccel)
        call getmem(3*na,iaccold)
        call getmem(3*na,iaolder)
c  in the zeroth dynamics step, just allocate memory and return
c  read the DYNA card in step 0 but do not act on it
        inp=igetival('inp')
        read(inp,*)
        return
      end if
c  get current coordinates
      open(1,file=jobname(1:lenJ)//'.coord',
     $     form='formatted',status='old')
      call RdCoordF(1,na,AtSymb,bl(icoord),-1,jnk,
     $              bl(iatchrg),bl(imass))
      close(1)
      write(iout,*) 'dynamain cycle=', istep
      write(iout,*) ' '
c  get current gradient
c  atomic mass unit in atomic units (electron mass units)
      call RdGrad(na,bl(iaccel),'save')
c
CTESTS
c  I put in a strict harmonic oscillator for CO in the z direction
C TAKE IT OUT!!!
c      if(jobname(1:lenJ).eq.'co') then
c        dr=abs(bl(icoord+2)-bl(icoord+5))-2.1423257d0
c        bl(iaccel+2)=1.3534105d0*dr
c        bl(iaccel+5)=-1.3534105d0*dr
c      end if
      call dynam(istep,na,AtSymb,bl(icoord),bl(iatchrg),
     1           bl(imass),bl(iold),bl(inew),bl(iaccel),bl(iaccold),
     2           bl(iaolder),finished)
c
c  write to trajectory file
      call WrTraj(na,istep,itraj,AtSymb,bl(icoord))
      end
c==============================================================
      subroutine dynam(istep, na,   atsymb,   xyz,   q,
     1                 xmass, xold, xnew,     accel, accold,
     2                 aolder,finished)
      implicit real*8(a-h,o-z)
      character*20 cdum
c  main routine for direct molecular dynamics
c  Arguments
c      INPUT:
c  istep is the step number. it is needed only
c  because steps 0 and 1 are different from the rest
c  na = number of atoms
c  atsymb(na) = atomic symbols
c  xyz(3,na) = initial coordinates
c  q(na) = nuclear charges
c  xmass = atomic masses (in amu units)
c      STORAGE:
c   the following quantities are used in the modified Verlet routine
c  xold(3,na)
c  xnew(3,na)
c  accel(3,na): on input, the forces. on output, the cceleration
c  accold(3,na)
c  aolder(3,na)
c     OUTPUT
c  there is no output argument. However, the coordinates in each
c  timestep are written on a file jobname.trajec
c
      character*8 atsymb(na)
      character*256 jobname
      dimension xyz(3,na),q(na),xmass(na),xold(na),xnew(na)
      dimension accel(3,na), accold(3,na),aolder(3,na)
      parameter(nopt=4)
      character*4 opnames(nopt)
      character*256 chopv(nopt)
      logical finished
      dimension ioptyp(nopt),iopv(3,nopt),ropv(3,nopt),ifound(nopt)
      save timestep,temp,maxcyc
      data opnames/'step','temp','maxc','seed'/
      data ioptyp /11,11,1,1/
c  read input options
      inp=igetival('inp')
      iout=igetival('iout')
c  transform the forces to accelerations
c     accel at this point contains the forces. It has to be
      amu=rgetrval('amu ')
      xme=rgetrval('me  ')
      ramu=xme/amu
c     transformed into acceleration
      do i=1,na
        xm=ramu/xmass(i)
        accel(1,i)=accel(1,i)*xm
        accel(2,i)=accel(2,i)*xm
        accel(3,i)=accel(3,i)*xm
      end do
c
      if(istep.eq.1) then
        call readopt(inp,nopt,opnames,ioptyp,iopv,ropv,chopv,ifound)
c default timestep is a thousand atomic units
        timestep=1000.0d0
c  default temperature is 298.15 K (25 C)
        temp=298.15d0
c  default number of cycles is 100
        maxcyc=100
        if(ifound(1).ge.1) then
          timestep=ropv(1,1)
        end if
        if(ifound(2).ge.1) then
          temp=ropv(1,2)
        end if
        if(ifound(3).ge.1) then
          maxcyc=iopv(1,3)
        end if
        iseed=0
        if(ifound(4).ge.1) then
          iseed=iopv(1,4)
        end if
c
        call velocity(na,timestep,xyz,xmass,temp,xnew,accel,iseed)
c        call printcoord(6,na,xnew)
        call WrCoord(na,AtSymb,xnew,1,jnk,q,xmass)
        call tfer(xyz,xold,3*na)
        call tfer(xnew,xyz,3*na)
        call zeroit(aolder,3*na)
        call tfer(accel,accold,3*na)
        return
      else
c  read in the (meaningless) DYNA card
        read(inp,*)
c       call modified (accurate through 4th order) Verlet
        call verlet2(na,  istep, timestep,xyz,   xold,
     1              xmass,accel, accold,  aolder,xnew)
        call WrCoord(na,AtSymb,xnew,1,jnk,q,xmass)
c  calculate the kinetic energy
        call kinetic(na,xnew,xold,xmass,timestep,amu/xme,EKin)
c  read the potential energy
        call getchval('jobname',jobname)
        call rmblan2(jobname,256,lenJ)
        IUnit=1
        open(IUnit,file=jobname(1:lenJ)//'.control',status='old')
        call RdCntrl(IUnit,7,'$energy',2,idum,Eelec,cdum)
        close(1)
        ETot=ekin+Eelec
        write(iout,1000) EKin,EElec,Etot
1000    format('T,V,Etot=',3f16.7)
        call tfer(xyz,xold,3*na)
        call tfer(xnew,xyz,3*na)
        call tfer(accold,aolder,3*na)
        call tfer(accel,accold,3*na)
      end if
      if(istep.ge.maxcyc) finished=.true.
      end
c====================================================
      subroutine verlet2(na,  istep, timestep,xyz,   xold,
     1                  xmass,accel, accold,  aolder,xnew)
      implicit real*8(a-h,o-z)
c  this routine makes a modified Verlet step
c  the algorithm is as follows:
c  x(t+dt)=2x(t)-x(t-dt)+a(t)*dt**2+(1/12)[a(t)+a(t-2dt)-2a(t-dt)]
c
c  arguments:
c  (time at the current step is denoted by t. this is not an argument)
c      INPUT:
c  na = number of atoms
c  istep = step number
c  timestep = dt = time increment
c  xyz(3,na) = current (at time t) x,y,z coordinates of the atoms
c  xold(3,na) = coordinates in the previous step
c  xmass(na) = atomic masses (in amu units)
c  accel(3,na) = on : acceleration at time t
c  accold = acceleration at t-dt
c  aolder = acceleration at t-2dt
c      OUTPUT
c  xnew(3,na) = estimated coordinates at t+dt
c
      parameter(zero=0.0d0,one=1.0d0,two=2.0d0,three=3.0d0,four=4.0d0)
      parameter(rtw=one/12.0d0)
      dimension xyz(3,na),xold(3,na),xmass(na),grad(3,na)
      dimension accel(3,na), accold(3,na),aolder(3,na),aest(3,na)
c
c      print *, 'before Verlet'
c      call printcoord(6,na,xyz)
c      print *,'Old coordinates'
c      call printcoord(6,na,xold)
      call tfer(xyz,xnew,3*na)
      call add(xyz,xnew,3*na)
      call add1(xold,-one,xnew,3*na)
c      print *,'2x-xold'
c      call printcoord(6,na,xnew)
c      print *, 'Acceleration'
c      call printcoord(6,na,accel)
      if(istep.eq.2)then
        call add1(accel,timestep**2,xnew,3*na)
c
c        print*, 'Verlet 1st step'
c        call printcoord(6,na,xnew)
c  X(i+1)=2X(i)-X(i-1)+A(i)*dt**2 where A(i) is the acceleration at i
c  this is the normal Verlet
c  in the first step, the 4th order method cannot be applied
        return
      else
c
c  in steps 3rd or higher, we can do better
c  X(i+1)=2X(i)-X(i-1)+[A(i)+(1/12)(A(i)+A(i-2)-2*A(i-1))]*dt**2
c  this formula is accurate to 4th order (inclusive) vs. 3rd for Verlet
c        print *, 'Data in modified Verlet'
cc        print *, 'accel'
c        call printcoord(6,na,accel)
c        print *, 'accold'
c        call printcoord(6,na,accold)
c        print *, 'aolder'
c        call printcoord(6,na,aolder)
c
c        call add1(accel,timestep**2,xnew,3*na)
c
c straight Verlet
        call add1(accel,timestep**2,xnew,3*na)
c  modified
c        call add1(accel,(one+rtw)*timestep**2,xnew,3*na)
c        call add1(aolder,rtw*timestep**2,xnew,3*na)
c        call add1(accold,-two*rtw*timestep**2,xnew,3*na)
c        print *,'Modified Verlet, 2nd or higher steps'
c        call printcoord(6,na,xnew)
      end if
      end
c======================================================================
      subroutine WrTraj(na,istep,iunit,symb,xyz)
      implicit real*8(a-h,o-z)
c   this routine writes the current coordinates to the file
c   <jobname>.trajec
c   Arguments
c  na = number of atoms
c  iunit = unit number of the trajectory file
      character*8 symb
      dimension xyz(3,na),symb(na)
c
      write(iunit,1000) istep
      do k=1,na
      write(iunit,1100) symb(k),xyz(1,k),xyz(2,k),xyz(3,k)
      end do
c
 1000 Format('Dynamics step: ',I8)
 1100 Format(a8,2x,3f15.7)
c
      end
c======================================================================
      subroutine velocity(na,   timestep, xyz, xmass, temp,
     1                    xnew,accel,iseed)
      implicit real*8 (a-h,o-z)
      real*4 rand
      dimension xyz(3,na),xmass(na),xnew(3,na),accel(3,na),totmom(3)
      dimension x(3),y(3),angm(3),tiner(3,3),center(3),z(3)
      parameter(zero=0.0d0,half=0.5d0,one=1.0d0,two=2.0d0)
c  this routine sets the initial velocities
c  Arguments
c     INTENT(IN):
c  na       = number of atoms
c  timestep = time step for dynamics
c  xyz      = initial geometry
c  xmass    = atomic masses (in amu units)
c  temp     = absolute temperature in K
c  xnew     = new coordinates (input and output)
c  accel    = acceleration at the initial geometry
c  iseed    = random seed, if non-zero, it will be used, otherwise generated
c
c  recall the Boltzmann constant
      bolt=rgetrval('bolt')
      ajoule=rgetrval('ajou')
      hartree=1.0d-18/ajoule
      pi=rgetrval('pi  ')
      iout=igetival('iout')
c  mass unit is me, not amu - need conversion factor
      amu=rgetrval('amu ')
      xme=rgetrval('me  ')
      amu=amu/xme
      call elapsec(telap)
      if(iseed.eq.0) then
        iseed=telap*100.0d0
        write(iout,*) 'Random seed = ',iseed
      end if
      call init_random( iseed )
c     call srand(iseed)
c  units are J/(mol*K)
      timestep2=timestep**2
      totmom(1)=zero
      totmom(2)=zero
      totmom(3)=zero
      totmass=zero
c  generate random velocities with a Gaussian (1-D Maxwell) distribution
c  this is not quite trivial, see G.E.P. Box and M.E. Muller,
c  "A note on the generation of random normal variates"
c  Ann. math. Stat. 29 (1958) 610-11
c  See also M.P. Allen and D.J. Tildesley, "Computer Simulation of
c  Liquids", Clarendon Press, 1987, p. 170 and p. 347 (App. G)
c  the algorithm is:
c  (1)generate two uniform random variates x1 and x2 on (0,1)
c  (2) calculate z1=(-2ln(x1))**(1/2)*cos(2*pi*x2) and
c                z2=(-2ln(x1))**(1/2)*sin(2*pi*x2)
c  z1 and z2 both have normal distibution with zero mean and unit variance
c
c  first assign velocities from a normal (Gaussian) distribution
c  with zero mean and unit variance
c  Calculate the total linear momentum
c
      call zeroit(center,3)
      do i=1,na
        totmass=totmass+xmass(i)
        do l=1,3
          center(l)=center(l)+xyz(l,i)*xmass(i)
          call random_number( x1 )
c         x1=rand(0)
          call random_number( x2 )
c         x2=rand(0)
          vli=sqrt(-two*log(x1))*cos(two*pi*x2)
          ekin=ekin+vli**2
          totmom(l)=totmom(l)+xmass(i)*vli
          xnew(l,i)=xyz(l,i)+vli
        end do
      end do
c
      do l=1,3
        totmom(l)=totmom(l)/totmass
        center(l)=center(l)/totmass
      end do
c  totmom at this point is not momentum but average velocity
c  center is the center of mass
c  Now shift the velocities so that the total momentum is zero
c  We still may have net angular momentum
      do i=1,na
        do l=1,3
          xnew(l,i)=xnew(l,i)-totmom(l)
        end do
CPP
c         print *, 'velocity of atom ',i,(xnew(k,i)-xyz(k,i),k=1,3)
      end do
c  calculate the total angular momentum
c  first zero out it and the inertia tensor
      call zeroit(angm,3)
      call zeroit(tiner,9)
      do i=1,na
c  omit the timestep dt; divide later
        r2m=zero
        do l=1,3
          z(l)=xyz(l,i)-center(l)
          r2m=r2m+xmass(i)*z(l)**2
        end do
        tiner(1,1)=tiner(1,1)+r2m
        tiner(2,2)=tiner(2,2)+r2m
        tiner(3,3)=tiner(3,3)+r2m
        call tfer(xnew(1,i),x,3)
        call add1(xyz(1,i),-one,x,3)
c  x at this point contains the velocities times dt
CPP
c        print *, 'Vel. of atom',i,x
c  the angular momentum contribution is the cross product of
c   the radius vector and the velocity times the mass
        call cross(z,x,y)
CPP
c        print *, 'Ang. mom. contr. y of atom', i,y
        call add1(y,xmass(i),angm,3)
c        angm(1)=angm(1)+xmass(i)*(xyz(2,i)*(xnew(3,i)-xyz(3,i))-
c     1           xyz(3,i)*(xnew(2,i)-xyz(2,i)))
c        angm(2)=angm(2)+xmass(i)*(xyz(3,i)*(xnew(1,i)-xyz(1,i))-
c     1           xyz(1,i)*(xnew(3,i)-xyz(3,i)))
c        angm(3)=angm(3)+xmass(i)*(xyz(1,i)*(xnew(2,i)-xyz(2,i))-
c     1           xyz(2,i)*(xnew(1,i)-xyz(1,i)))
        do l=1,3
          do k=1,3
            tiner(k,l)=tiner(k,l)-xmass(i)*z(k)*z(l)
          end do
        end do
      end do  ! i (atoms) loop
CPP
c      print *, 'Inertia tensor', tiner
c      print *, 'Total angular momentum', angm
c  calculate the generalized inverse of the inertial tensor
      call geninv(tiner,3,det,1.0d-4,x,y,ninv)
CPp
c      print *, 'ninv',ninv
c      print *, 'Gen. inv.', tiner
c  multiply the angular momemtum tensor with the generalized inverse
c  of the inertia tensor
      do l=1,3
        sum=zero
        do k=1,3
          sum=sum+tiner(k,l)*angm(l)
        end do
        y(l)=sum
      end do
CPP
c      print *,'Ang. velocity', y
c  y(1:3) contins the (average) angular velocity of the system
c  now change the velocities so that the ang. mom. is zero
c  the contribution of an angular velocity w to the velocity
c  of a point ar R is Rxw (R cross product w)
c  Also calculate the kinetic energy.
      ekin=zero
      do i=1,na
        do l=1,3
          z(l)=xyz(l,i)-center(l)
        end do
        call cross(z,y,x)
CPP
c        print *, 'Corr. vector for atom',i,x
        call add(x,xnew(1,i),3)
        do l=1,3
          vli=xnew(l,i)-xyz(l,i)
          ekin=ekin+vli**2*xmass(i)
        end do
      end do
c  divide by the timestep squared and include the 1/2 factor
c  also include the amu - me conversion factor
      ekin=amu*half*ekin/timestep2
c  average kinetic energy per degree of freedom
      nfreedom=3+ninv
c  ninv was the number of non-zero roots of the inertia matrix
      eav=ekin/dble(3*na-nfreedom)
c  eav is in atomic units. eav*hartree is in Joules
c  It is difficult to set the temperature at the outset
c  the average kinetic energy should be 1/2kT
c  Therefore, the kinetic energy when started at the equilibrium
c    should be kT for vibrations, as in a classical equipartition
c    half of the energy in a collection of oscillators is kinetic E
      scale=sqrt((bolt*temp)/(eav*hartree))
c      print *, 'Ekin, EAv,kT,Scale', ekin,eav,bolt*temp,scale
c  now re-scale the velocities to give the correct kinetic energy
c  and add the acceleration contribution at the starting point
c  note that the potential energy is assumed to be zero
c      print *, 'scale=',scale
c      print *, 'atom, coord,xyz,vli,xnew'
      do i=1,na
        do l=1,3
          vli=xnew(l,i)-xyz(l,i)
          xnew(l,i)=xyz(l,i)+scale*vli+half*accel(l,i)*timestep2
c          print *, i,l,xyz(l,i),vli,xnew(l,i)
        end do
      end do
c      call printcoord(6,na,xnew)
c check the kinetic energy
      call kinetic(na,xnew,xyz,xmass,timestep*half,amu,ekin)
      write(iout,*) 'Initial kinetic energy=',ekin
      end
c=======================================================================
      subroutine printcoord(iunit,na,xyz)
      implicit real*8 (a-h,o-z)
      dimension xyz(3,na)
      write(iunit,*) 'Coordinates'
      do i=1,na
        write(iunit,*) i, (xyz(l,i),l=1,3)
      end do
      end
c======================================================================
      subroutine kinetic(na,xnew,xold,xmass,timestep,amu,ekin)
      implicit real*8 (a-h,o-z)
c  This routine calculates the kinetic energy in a Verlet step
c  or modified Verlet step
c  Arguments:
c  INPUT
c  na = number of atoms
c  xnew(3,na) = new Cartesian coordinates
c  xold(3,na) = old Cartesian coordinates
c  timestep = time step in atomic units (appr. 0.025 fs)
c  xmass = atomic masses IN AMU (not me) units
c  amu = conversion factor for amu
c  OUTPUT:
c  ekin = kinetic energy in hartrees
      dimension xnew(3,na),xold(3,na),xmass(na)
      parameter (zero=0.0d0)
      ekin=zero
      do i=1,na
        do l=1,3
          ekin=ekin+xmass(i)*(xnew(l,i)-xold(l,i))**2
        end do
      end do
      ekin=0.125d0*ekin*amu/timestep**2
c  the factor of 0.125 comes from 1/2 in the kinetic energy formula,
c  plus the fact that the timestep from xold to xnew is 2*timestep
      end
c======================================================================
      subroutine geninv (c,nq,d,tol,x,y,ninv)
C     GENERALIZED INVERSE ROUTINE
C     DIAGONALIZE, INVERT IF THE DIAGONAL IS NOT TOO SMALL,
C     INVERT BACK
c     c= matrix to get the gen,. inverse
c      nq=dimension of it
c     d it the determinant - no function here
c     tol a small number, if an eigenvalue of c is less than 100*tol
c     then this root is not inverted
c     x and y are two temporary vectors nq long
c     ninv is a result: the number of non-zero (not small) eigenvalues
      IMPLICIT REAL*8 (A-H,O-Z)
      dimension c(nq,nq),x(*),y(*)
      tol1=100.0d0*tol
      call sdiag2(nq,nq,c,x,c)
      det =1.0d0
      ninv=0
      do 100 i=1,nq
         if(abs(x(i)).gt.tol1) then
           x(i)=1.0d0/x(i)
               det=det*x(i)
           ninv=ninv+1
         end if
      d=det
 100  continue
      do 400 i=1,nq
        do 300 j=i,nq
          sum=0.0d0
            do 200 k=1,nq
              sum=sum+c(i,k)*c(j,k)*x(k)
 200        continue
            y(j)=sum
 300     continue
         do 400 j=i,nq
           c(i,j)=y(j)
 400     continue
 500  continue
      do 600 i=1,nq
        do 600 j=i+1,nq
          c(j,i)=c(i,j)
 600  continue
      end

      SUBROUTINE para_cmp2
      IMPLICIT REAL*8(A-H,O-Z)
C
C
C  ............................................................
C    Wrapper for MP2 module
C    Depending on number of alpha/beta electrons and
C    multiplicity calls Closed or Open-shell MP2
C  ............................................................
C
      Character*20 cdum
      Character*256 jobname
      PARAMETER (IUnit=1)          ! unit for checkpoint read
c
      Common /job/jobname,lenJ
C
C
C  read from <control> file
C  multiplicity
C  number of alpha electrons
C  number of beta electrons
C
      OPEN (UNIT=IUnit,FILE=jobname(1:lenJ)//'.control',
     $      FORM='FORMATTED',STATUS='OLD')
      call rdcntrl(IUnit,13,'$multiplicity',1,IMult,rdum,cdum)
      call rdcntrl(IUnit,7,'$nalpha',1,NAlpha,dum,cdum)
      call rdcntrl(IUnit,6,'$nbeta',1,NBeta,dum,cdum)
      CLOSE (UNIT=IUnit,STATUS='KEEP')
c
c -- check there are electrons!
      If(NAlpha.EQ.0) Call nerror(1,'MP2 Module',
     $                     'There Are No Electrons!',0,0)
c
      If(NBeta.EQ.0.AND.IMult.EQ.1) Then
        Call para_rmp2(NAlpha)
      Else
        Call para_ump2(NAlpha,NBeta)
      EndIF
c
      RETURN
      END

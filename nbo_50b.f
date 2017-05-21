C***********************************************************************
      SUBROUTINE INITWT(ERRREF,ERR,RHOSTR,OCC,WGT,Q,IST,IRANK,NRES,
     +                  IRESET,MAXRES,NREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION ERR(NREF),RHOSTR(MAXRES,NREF),OCC(NDIM,MAXRES,NREF),
     +          WGT(MAXRES,NREF),Q(MAXRES,MAXRES),IST(NDIM,NREF),
     +          IRANK(MAXRES),NRES(NREF),ERRREF(NREF),IRESET(NREF)
C
      SAVE ONE,THREE
      DATA ONE,THREE/1.0D0,3.0D0/
C
C  Set RHOW threshold:
C
      RHOW = ONE / THREE
      IF(ISPIN.NE.0) RHOW = RHOW / 2
C
C  Initialize the pointer array, IRANK:
C
      DO 10 IRES = 1,MAXRES
        IRANK(IRES) = IRES
   10 CONTINUE
C
C  Loop over reference resonance structures:
C
      DO 20 IREF = 1,NREF
        IF(IRESET(IREF).GT.0) THEN
          IF(NRES(IREF).GT.0) THEN
C
C  Prepare the Q matrix:
C
            CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,0)
C
C  First compute the reference density error (i.e. the density error
C  when only the parent structure is considered):
C
            CALL EXPWGT(RHOW,RHOSTR,WGT,IRANK,NRES,MAXRES,IREF,NREF,1)
            CALL GETERR(ERRREF(IREF),WGT,Q,NRES,MAXRES,IREF,NREF)
C
C  Then compute the density error with all resonance considered:
C
            CALL EXPWGT(RHOW,RHOSTR,WGT,IRANK,NRES,MAXRES,IREF,NREF,
     +                  NRES(IREF))
            CALL GETERR(ERR(IREF),WGT,Q,NRES,MAXRES,IREF,NREF)
          END IF
        END IF
   20 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE WGTOPT(ERRREF,ERR,RHOSTR,OCC,WGT,T,SCR1,SCR2,SCR3,SCR4,
     +                  SCR5,SCR6,SCR7,SCR8,SCR9,Q,IST,IPTR,KPAR,NRES,
     +                  IRESET,MAXRES,NREF,NRTCTL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION ERR(NREF),RHOSTR(MAXRES,NREF),OCC(NDIM,MAXRES,NREF),
     +          WGT(MAXRES,NREF),T(MAXRES),SCR1(MAXRES,MAXRES),
     +          SCR2(MAXRES,MAXRES),SCR3(MAXRES),SCR4(MAXRES),
     +          SCR5(MAXRES),SCR6(MAXRES),SCR7(MAXRES),SCR8(MAXRES),
     +          SCR9(MAXRES),Q(MAXRES,MAXRES),IST(NDIM,NREF),
     +          IPTR(MAXRES),ERRREF(NREF),NRES(NREF),IRESET(NREF),
     +          KPAR(MAXRES,NREF),NRTCTL(10)
C
      SAVE ZERO,THRE,THRW,EPS,EPSX,EPSY,ONE
      SAVE IB,IP
      DATA ZERO,THRE,THRW,EPS,EPSX,EPSY,ONE
     +    /0.0D0,1.0D-10,1.0D-3,1.0D-2,1.0D-2,0.0D0,1.0D0/
      DATA IB,IP/1HB,1HP/
C
C  Loop over reference resonance structures, optimizing a set of resonance
C  weights for each:
C
      DO 50 IREF = 1,NREF
        IF(IRESET(IREF).GT.0) THEN
          XLAM = ZERO
          YLAM = ZERO
          ITER = 0
          IF(NRES(IREF).GT.1) THEN
C
C  Prepare the Q matrix:
C
            CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,0)
C
C  Optimize with analytic derivatives:
C
            IF(JPRINT(55).EQ.0.OR.JPRINT(55).EQ.IB) THEN
              IF(JPRINT(57).NE.0) WRITE(LFNPR,910)
              IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
              CALL SETOPT(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
              CALL BFGS(ERR(IREF),ERRREF(IREF),OCC,WGT,Q,T,SCR1,SCR2,
     +                  SCR3,SCR4,SCR5,SCR6,SCR7,SCR8,SCR9,IPTR,
     +                  NRES,MAXRES,IREF,NREF,NVAR,ITER)
C
C  Or optimize with numerical derivatives:
C
            ELSE IF(JPRINT(55).EQ.IP) THEN
              IF(JPRINT(57).NE.0) WRITE(LFNPR,920)
              IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
              CALL SETOPT(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
              CALL NBPWLL(ERR(IREF),ERRREF(IREF),WGT,Q,T,SCR1,SCR2,
     +                    SCR3,SCR4,SCR5,SCR6,IPTR,NRES,MAXRES,IREF,
     +                    NREF,NVAR,ITER)
C
C  Or optimize with a simulated annealing algorithm:
C
            ELSE
              ISEED = -ABS(JPRINT(55)) - IREF
              IF(JPRINT(55).GT.0) XLAM = EPSX
              IF(JPRINT(57).NE.0) WRITE(LFNPR,900)
              IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
              CALL ANNEAL(XLAM,YLAM,ERRREF,ERR,RHOSTR,WGT,T,SCR2,Q,
     +                   IPTR,NRES,MAXRES,IREF,NREF,ISEED,NVAR,ITER)
            END IF
            IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
            IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
C
C  Evaluate the Hessian, and if zero eigenvalues are encountered, repeat
C  the optimization of the weights with ANNEAL and the penalty function:
C
            IF(XLAM.EQ.ZERO) THEN
              CALL HESSIN(RHOSTR,OCC,WGT,Q,T,SCR1,SCR2,SCR3,SCR4,SCR5,
     +                    IPTR,NRES,MAXRES,IREF,NREF,IND)
              DO 10 IV = 1,NVAR
                IF(SCR3(IV).LT.THRE) XLAM = EPSX
   10         CONTINUE
              IF(XLAM.NE.ZERO) THEN
C
C  So, reoptimize with the an added penalty function.  But, be careful.
C  If the weight of the reference structure has dropped to zero during
C  the regular optimization, ANNEAL is going to get upset here.  So, if
C  needed, set the weight of the reference structure to 1, and all others
C  to 0:
C
                IF(WGT(1,IREF).LT.EPS) THEN
                  WGT(1,IREF) = ONE
                  DO 20 IRES = 2,NRES(IREF)
                    WGT(IRES,IREF) = ZERO
   20             CONTINUE
                END IF
                IF(JPRINT(57).NE.0) WRITE(LFNPR,950)
                IF(JPRINT(55).EQ.0)  ISEED = -1234567
                IF(JPRINT(55).EQ.IB) ISEED = -1234567
                IF(JPRINT(55).EQ.IP) ISEED = -1234567
                IF(JPRINT(57).NE.0) WRITE(LFNPR,900)
                IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,
     +                                          NREF)
                CALL ANNEAL(XLAM,YLAM,ERRREF,ERR,RHOSTR,WGT,T,SCR2,
     +                 Q,IPTR,NRES,MAXRES,IREF,NREF,ISEED,NVAR,ITER)
                IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,
     +                                          NREF)
                IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
              END IF
            END IF
C
C  Write out results of optimization:
C
            FW = ONE - ERR(IREF) / ERRREF(IREF)
            IF(ITER.LT.0) THEN
              WRITE(LFNPR,930) IREF,RHOSTR(1,IREF),FW,ABS(ITER)
            ELSE
              WRITE(LFNPR,940) IREF,RHOSTR(1,IREF),FW,ITER
            END IF
          ELSE
            WRITE(LFNPR,940) IREF,RHOSTR(1,IREF),ZERO,ITER
          END IF
C
C  Neglect this manifold if the weight of the reference is less than
C  threshold and mark the structure so that it isn't selected as
C  reference in the future:
C
          IF(WGT(1,IREF).LT.THRW.AND.NRTCTL(3).NE.0) THEN
            WRITE(LFNPR,970)
            IRESET(IREF) = -2
            KREF = KPAR(1,IREF) / 10000
            KRES = MOD(KPAR(1,IREF),10000)
            IF(KRES.NE.0.AND.KREF.NE.0) KPAR(KRES,KREF) = 1
C
C  If the weight of the reference structure is less than any of the
C  secondaries, reoptimize with simulated annealing and second penalty
C  function:
C
          ELSE
            N = 0
            DO 30 IRES = 2,NRES(IREF)
              IF(WGT(IRES,IREF).GT.WGT(1,IREF)) N = N + 1
   30       CONTINUE
            IF(N.NE.0) THEN
              WRITE(LFNPR,990) N
              YLAM = EPSY
              IF(YLAM.EQ.ZERO) GOTO 45
              IF(WGT(1,IREF).LT.EPS) THEN
                WGT(1,IREF) = ONE
                DO 40 IRES = 2,NRES(IREF)
                  WGT(IRES,IREF) = ZERO
   40           CONTINUE
              END IF
              IF(JPRINT(57).NE.0) WRITE(LFNPR,1000)
              IF(JPRINT(55).EQ.0)  ISEED = -1234567
              IF(JPRINT(55).EQ.IB) ISEED = -1234567
              IF(JPRINT(55).EQ.IP) ISEED = -1234567
              IF(JPRINT(57).NE.0) WRITE(LFNPR,900)
              IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
              CALL ANNEAL(XLAM,YLAM,ERRREF,ERR,RHOSTR,WGT,T,SCR2,Q,
     +                   IPTR,NRES,MAXRES,IREF,NREF,ISEED,NVAR,ITER)
              IF(JPRINT(57).NE.0) CALL WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
              IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
   45         CONTINUE
            END IF
          END IF
C
C  Other cases:
C
        ELSE IF(IRESET(IREF).EQ.-3) THEN
          IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
          WRITE(LFNPR,980) IREF,RHOSTR(1,IREF)
        ELSE IF(IRESET(IREF).EQ.-4) THEN
          IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
          WRITE(LFNPR,960) IREF,RHOSTR(1,IREF)
        END IF
   50 CONTINUE
      RETURN
C
  900 FORMAT(//1X,'Simulated annealing optimization of weights (SR ',
     + 'ANNEAL):')
  910 FORMAT(//1X,'Variational optimization with analytic gradients',
     + ' (SR BFGS):')
  920 FORMAT(//1X,'Variational optimization with numerical gradients',
     + ' (SR NBPWLL):')
  930 FORMAT(1X,'Reference',I4,':  rho*=',F7.5,', f(w)=',F7.5,
     + ' did not converge after',I4,' iterations')
  940 FORMAT(1X,'Reference',I4,':  rho*=',F7.5,', f(w)=',F7.5,
     + ' converged after',I4,' iterations')
  950 FORMAT(1X,'Linear dependency detected:  Repeating optimization',
     + ' w/penalty function enabled')
  960 FORMAT(1X,'Reference',I4,':  rho*=',F7.5,', ionic; deleted')
  970 FORMAT(1X,'                Low weight, neglected in multi-',
     + 'reference treatment')
  980 FORMAT(1X,'Reference',I4,':  rho*=',F7.5,', ionic; replaced ',
     + 'by covalent form')
  990 FORMAT(1X,'Warning:  Reference structure has lower weight than',
     + I3,' of the secondaries')
 1000 FORMAT(1X,'          Repeating optimization w/penalty function',
     + ' enabled')
      END
C***********************************************************************
      SUBROUTINE HESSIN(RHOSTR,OCC,WGT,Q,T,HESS,EVEC,EVAL,DER,G,
     +                  IPTR,NRES,MAXRES,IREF,NREF,IND)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION RHOSTR(MAXRES,NREF),OCC(NDIM,MAXRES,NREF),
     +          WGT(MAXRES,NREF),Q(MAXRES,MAXRES),T(MAXRES),
     +          HESS(MAXRES,MAXRES),EVEC(MAXRES,MAXRES),EVAL(MAXRES),
     +          DER(MAXRES),G(MAXRES),IPTR(MAXRES),
     +          NRES(NREF)
C
      SAVE ZERO,EPS,ONE,HUNDRD
      DATA ZERO,EPS,ONE,HUNDRD/0.0D0,1.0D-5,1.0D0,1.0D2/
C
C  Setup pointer array:
C
      DO 10 IRES = 1,NRES(IREF)
        IPTR(IRES) = -1
   10 CONTINUE
C
C  The weight of the reference structure will be determined by normalization:
C
      IPTR(1) = 0
      WGT(1,IREF) = -ONE
C
C  Constrain resonance weights with equivalent rho* to have equal weights:
C
      NVAR = 0
   30 JRES = 0
      RHO = HUNDRD
      DO 40 IRES = 1,NRES(IREF)
        IF(RHOSTR(IRES,IREF).GE.ZERO.AND.RHOSTR(IRES,IREF).LT.RHO) THEN
          IF(WGT(IRES,IREF).GE.ZERO) THEN
            JRES = IRES
            RHO = RHOSTR(JRES,IREF)
          END IF
        END IF
   40 CONTINUE
C
C  Compute the optimizable parameter T and set pointers appropriately:
C
      IF(JRES.NE.0) THEN
        NVAR = NVAR + 1
        T(NVAR) = -LOG(WGT(JRES,IREF))
        DO 50 IRES = 1,NRES(IREF)
          IF(ABS(RHOSTR(IRES,IREF)-RHO).LT.EPS) THEN
            IF(WGT(IRES,IREF).GE.ZERO) THEN
              IPTR(IRES) = NVAR
              WGT(IRES,IREF) = -ONE
            END IF
          END IF
   50   CONTINUE
        GOTO 30
      END IF
C
C  Recompute the resonance weights:
C
      NLEAD = 0
      QQ = ZERO
      DO 60 IRES = 1,NRES(IREF)
        WGT(IRES,IREF) = ZERO
        IF(IPTR(IRES).EQ.0) THEN
          NLEAD = NLEAD + 1
        ELSE IF(IPTR(IRES).GT.0) THEN
          WGT(IRES,IREF) = EXP(-T(IPTR(IRES)))
          QQ = QQ + WGT(IRES,IREF)
        END IF
   60 CONTINUE
      QQ = (ONE - QQ) / NLEAD
      DO 70 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) THEN
          WGT(IRES,IREF) = QQ
        END IF
   70 CONTINUE
C
C  Compute the hessian:
C
      CALL GETHES(HESS,DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,IREF,
     +            NREF,NVAR,IFLG)
      IF(IFLG.EQ.1) THEN
        IND = 0
        RETURN
      END IF
C
C  Compute the eigenvalues of the hessian:
C
      CALL NBJACOBI(NVAR,HESS,EVAL,EVEC,MAXRES,MAXRES,0)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SUPPL(FW,WGT,WGTM,MAP,WGTP,IDXRES,KPAR,NRES,IRESET,
     +                 IPTR,MAXRES,ICNT,NREF,MAXREF,NRTCTL,IVALSP)
C***********************************************************************
C 10-Jun-08  EDG  Initialize routine to correct errors in beta analysis
C 25-Oct-01  CMM  Don't simultaneously delete & add a low weight structure.
C                 Don't add a reference structure more than once.
C  2-Jun-98  EDG  Fix several bugs
C 28-Jan-94  JKB  Add IVALSP=2 to restart NRT if no covalent structures and
C                    truncated matrix, else only ionic structures possible
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,KOPT,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
C
      DIMENSION WGT(MAXRES,MAXREF),WGTP(MAXREF),WGTM(MAXRES*MAXREF),
     +          MAP(MAXRES*MAXREF),IDXRES(MAXRES,MAXREF),NRES(MAXREF),
     +          IRESET(MAXREF),IPTR(MAXREF),KPAR(MAXRES,MAXREF),
     +          NRTCTL(10),IDXTMP(9999)
C
      SAVE ZERO,EPS,FWMAX,KREF,IDXTMP,NTMP,NABTMP
      DATA ZERO,EPS,NTMP,IDXTMP,NABTMP/0.0D0,1.0D-5,0,9999*0,0/
C
C  Don't allow for cycles of the reference list if $NRTSTR is found
C  (NRTCTL(1)=-2):
C
      IF(NRTCTL(1).EQ.-2) RETURN
C
C  Initialize routine:
C
      IF(NRTCTL(1).EQ.-1) THEN
        FWMAX = ZERO
        KREF = 0
      END IF
C
C  Abort search for reference structures if oscillating:
C
      IVALSP = 0
      KCNT = 0
      KERR = 0
      IFORCE = 0
      DO 5 IREF = 1,MAXREF
        IF(IRESET(IREF).GE.0) KCNT = KCNT + 1
        IF(IRESET(IREF).LT.-1) KERR = KERR + 1
    5 CONTINUE
      IF(NRTCTL(1).LT.-1) THEN
        FWMAX = FW
        KREF  = KCNT
        IF(KCNT.EQ.1) FWMAX = ZERO
      ELSE
        NRTCTL(1) = 0
        IF(ABS(FW-FWMAX).LT.EPS.AND.KCNT.LE.KREF) THEN
          IF(KERR.EQ.0) RETURN
        END IF
        IF(FW.GT.FWMAX) THEN
          FWMAX = FW
          KREF  = KCNT
        END IF
      END IF
C
C  Reinitialize NTMP for beta spin analysis:
C
      IF(BETA) THEN
        IF(NABTMP.EQ.0) NTMP = 0
        NABTMP = 1
      END IF
C
C  Initialize search for additional reference structures:
C
C  IRESET(IREF):  Reference structure status flag
C        1 : new structure from last pass (accept according to wgt)
C        0 : old structure (accept according to wgt)
C       -1 : unoccupied position in reference list
C       -2 : structure with low weight (delete)
C       -3 : covalent structure from CONDNS (accept for next pass)
C       -4 : an ionic structure from CONDNS (delete)
C
      DO 10 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.0) IRESET(IREF) = 1
   10 CONTINUE
C
C  Find the resonance structure of highest weight:
C
      WGTMAX = ZERO
      DO 15 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          DO 14 IRES = 1,NRES(IREF)
            WGTMAX = MAX(WGTMAX,WGT(IRES,IREF)*WGTP(IREF))
   14     CONTINUE
        END IF
   15 CONTINUE
C
C  Count the number of structures in MAP (WGTM) having weight larger
C  than THRESH:
C
      JCNT = 0
      THRESH = ABS(REFTHR) * WGTMAX
      DO 20 I = 1,ICNT
        IF(WGTM(I).LT.THRESH) GOTO 30
        JCNT = I
   20 CONTINUE
   30 CONTINUE
C
C  Accept all reference structures with weight greater than THRESH:
C
      DO 50 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.1) THEN
          DO 40 I = 1,JCNT
            IF(IDXRES(1,IREF).EQ.ABS(MAP(I))) THEN
              IF(WGTM(I).GT.THRESH) THEN
                IRESET(IREF) = 0
                MAP(I) = -ABS(MAP(I))
                GOTO 50
              END IF
            END IF
   40     CONTINUE
          IFORCE = 1
        END IF
   50 CONTINUE
      N = 0
      DO 60 IREF = 1,MAXREF
        IF(IRESET(IREF).NE.0.AND.IRESET(IREF).NE.-1) THEN
          N = N + 1
          IPTR(N) = IREF
        END IF
   60 CONTINUE
      IF(N.GT.0) THEN
        IF(JPRINT(57).NE.0) THEN
          WRITE(LFNPR,910)
          WRITE(LFNPR,920) (IPTR(I),I=1,N)
        END IF
      END IF
C
C  CMM  If a structure is deleted, don't add it back.
C
      DO 62 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.-2) THEN
          DO 61 I = 1,JCNT
            IF(IDXRES(1,IREF).EQ.MAP(I)) MAP(I) = -ABS(MAP(I))
   61     CONTINUE
        END IF
   62 CONTINUE
C
C  Add covalent partners of the ionic structures deleted in CONDNS:
C  CMM Don't add in a structure for the second time.
C
      DO 65 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.-3) THEN
          IDLFLG = 0
          DO 63 ITMP = 1,NTMP
            IF(IDXRES(1,IREF).EQ.IDXTMP(ITMP)) IDLFLG = 1
   63     CONTINUE
          IF(IDLFLG.EQ.0) THEN
            NTMP = NTMP + 1
            IDXTMP(NTMP) = IDXRES(1,IREF)
            IRESET(IREF) = 2
            DO 64 I = 1,JCNT
              IF(MAP(I).EQ.IDXRES(1,IREF)) MAP(I) = -ABS(MAP(I))
   64       CONTINUE
          ELSE
            IRESET(IREF) = -4
          END IF
        END IF
   65 CONTINUE
C
C  CMM If a structure has been added in the past, don't add it again.
C
          DO 67 ITMP = 1,NTMP
            DO 66 I=1,JCNT
              IF(MAP(I).EQ.IDXTMP(ITMP)) MAP(I) = -ABS(MAP(I))
   66       CONTINUE
   67     CONTINUE
C
C  Add secondary structures of high weight to the reference list
C  (neglect secondary structures from KEKULE/APPEND):
C
      MORE = 0
      DO 80 I = 1,JCNT
        IF(MAP(I).GT.0) THEN
          DO 70 IREF = 1,MAXREF
            IF(IRESET(IREF).NE.-2.AND.
     +         IRESET(IREF).NE. 0.AND.
     +         IRESET(IREF).NE. 1.AND.
     +         IRESET(IREF).NE. 2) THEN
              DO 69 JREF = 1,MAXREF
                IF(IRESET(JREF).EQ.-2.OR.
     +             IRESET(JREF).EQ. 0.OR.
     +             IRESET(JREF).EQ. 1) THEN
                  DO 68 JRES = 1,NRES(JREF)
                    IF(KPAR(JRES,JREF).EQ.0) THEN
                      IF(IDXRES(JRES,JREF).EQ.MAP(I)) THEN
                        IF(WGTM(I).GT.THRESH) THEN
                          NTMP = NTMP + 1
                          IDXTMP(NTMP) = MAP(I)
                          IDXRES(1,IREF) = MAP(I)
                          IRESET(IREF) = 2
                          GOTO 80
                        END IF
                      END IF
                    END IF
   68             CONTINUE
                END IF
   69         CONTINUE
              GOTO 80
            END IF
   70     CONTINUE
          MORE = MORE + 1
        END IF
   80 CONTINUE
      IF(MORE.GT.0) THEN
        WRITE(LFNPR,900) MAXREF,MAXREF+MORE
        NRTCTL(1) = 0
        RETURN
      END IF
C
      DO 85 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.-2) IRESET(IREF) = -1
        IF(IRESET(IREF).EQ. 1) IRESET(IREF) = -1
   85 CONTINUE
C
C  Now search the secondary lists for structures of high weight:
C
      DO 140 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.2) THEN
C
C  First locate the secondary structure of highest weight:
C
          KRES = 0
          KREF = 0
          WGTMAX = ZERO
          DO 100 JREF = 1,MAXREF
            IF(IRESET(JREF).EQ.0) THEN
              DO 90 JRES = 2,NRES(JREF)
                IF(KPAR(JRES,JREF).EQ.0) THEN
                  IF(IDXRES(JRES,JREF).EQ.IDXRES(1,IREF)) THEN
                    WGTRES = WGTP(JREF) * WGT(JRES,JREF)
                    IF(WGTRES.GE.THRESH.AND.WGTRES.GT.WGTMAX) THEN
                      WGTMAX = WGTRES
                      KRES = JRES
                      KREF = JREF
                    END IF
                  END IF
                END IF
   90         CONTINUE
            END IF
  100     CONTINUE
C
C  Insert the KRESth structure of the KREFth reference manifold into
C  the IREFth location of the reference list:
C
          IRESET(IREF) = 1
          NRES(IREF)   = 1
          KPAR(1,IREF) = KREF*10000 + KRES
          IF(KPAR(1,IREF).EQ.0) GOTO 140
C
C  If the KRESth structure of KREF came from KEKULE, mark it so that
C  it won't be used again:
C
          IF(KPAR(KRES,KREF).EQ.0) KPAR(KRES,KREF) = 1
C
C  And search the KREFth manifold for other structures of identical
C  TOPO matrix with high weight.  Also add these to the reference list:
C
          DO 130 JRES = 2,NRES(KREF)
            IF(JRES.NE.KRES) THEN
              IF(KPAR(JRES,KREF).EQ.0) THEN
                IF(IDXRES(JRES,KREF).EQ.IDXRES(1,IREF)) THEN
                  WGTRES = WGTP(KREF) * WGT(JRES,KREF)
                  IF(WGTRES.GT.THRESH) THEN
                    DO 110 JREF = 1,MAXREF
                      IF(IRESET(JREF).LE.-1) THEN
                        IDXRES(1,JREF) = IDXRES(1,IREF)
                        NRES(JREF)   = 1
                        IRESET(JREF) = 1
                        KPAR(1,JREF) = KREF*10000 + JRES
                        GOTO 120
                      END IF
  110               CONTINUE
                    CALL NBHALT('Increase NRTMEM in the $NBO keylist.')
  120               CONTINUE
                  END IF
                END IF
              END IF
            END IF
  130     CONTINUE
        END IF
  140 CONTINUE
C
      N = 0
      DO 150 IREF = 1,MAXREF
        IF(IRESET(IREF).EQ.1) THEN
          NRTCTL(1) = 1
          N = N + 1
          IPTR(N) = IREF
        END IF
  150 CONTINUE
      IF(IFORCE.EQ.1) NRTCTL(1) = 1
      IF(N.GT.0) THEN
        IF(JPRINT(57).NE.0) THEN
          WRITE(LFNPR,930)
          WRITE(LFNPR,920) (IPTR(I),I=1,N)
        END IF
      END IF
C
C  Reset:
C
      DO 160 IREF = 1,MAXREF
        IRESET(IREF) = MAX(IRESET(IREF),-1)
  160 CONTINUE
      NREF = 0
      DO 170 I = 1,ICNT
        MAP(I) = ABS(MAP(I))
  170 CONTINUE
      DO 180 IREF = MAXREF,1,-1
        IF(NREF.EQ.0.AND.IRESET(IREF).GE.0) THEN
          NREF = IREF
          RETURN
        END IF
  180 CONTINUE
C
C If SUPPL got this far without RETURNing, only ionic structures were
C    found and no corresponding covalent structure was generated.  Try
C    restarting with full AO matrix if not done already, else quit.
C
      IF(JPRINT(77).EQ.0) THEN
        WRITE(LFNPR,940)
        IVALSP = 2
      ELSE
        WRITE(LFNPR,950)
        IVALSP = 3
      END IF
      RETURN
C
  900 FORMAT(/1X,'The NRT program is currently configured to handle ',
     + I2,' reference structures',/1X,'which is insufficient for the',
     + ' present calculation.  Increase this number',/1X,'to least ',
     + I2,' using the NRTMEM keyword of the $NBO keylist.')
  910 FORMAT(/1X,'Deleting the following reference structures:')
  920 FORMAT(1X,19I4)
  930 FORMAT(/1X,'Adding the following reference structures:')
  940 FORMAT(/1X,'No covalent structures available within minimal ',
     +   'valence space (SR SUPPL).',/1X,'Reallocate scratch vector',
     +   ' and restart using full AO density matrix.')
  950 FORMAT(/1X,'No covalent reference structures found; ',
     +   'abandon NRT analysis.')
      END
C***********************************************************************
      SUBROUTINE RESWGT(RHOSTR,WGT,WGTP,WGTM,ISCR,MAP,IRANK,NRES,IRESET,
     +                  IDXRES,MAXRES,ICNT,NREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
C
      DIMENSION RHOSTR(MAXRES,NREF),WGT(MAXRES,NREF),WGTP(NREF),
     +          WGTM(MAXRES*NREF),ISCR(MAXRES*NREF),MAP(MAXRES*NREF),
     +          IRANK(MAXRES*NREF),NRES(NREF),IRESET(NREF),
     +          IDXRES(MAXRES,NREF)
C
      SAVE ZERO
      DATA ZERO/0.0D0/
C
C  Form a list of unique resonance structures and sum weights:
C
      ICNT = 0
      DO 30 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          DO 20 IRES = 1,NRES(IREF)
            IF(RHOSTR(IRES,IREF).GE.ZERO) THEN
              IFLG = 0
              DO 10 IC = 1,ICNT
                IF(IDXRES(IRES,IREF).EQ.ISCR(IC)) IFLG = IC
   10         CONTINUE
              IF(IFLG.EQ.0) THEN
                ICNT = ICNT + 1
                ISCR(ICNT) = IDXRES(IRES,IREF)
                WGTM(ICNT) = WGTP(IREF) * WGT(IRES,IREF)
              ELSE
                WGTM(IFLG) = WGTM(IFLG) + WGTP(IREF) * WGT(IRES,IREF)
              END IF
            END IF
   20     CONTINUE
        END IF
   30 CONTINUE
      IF(ICNT.EQ.0) RETURN
C
C  Sort the resonance structures by weight:
C
      CALL RANK(WGTM,ICNT,ICNT,IRANK)
      DO 40 IC = 1,ICNT
        MAP(IC) = ISCR(IRANK(IC))
   40 CONTINUE
      IF(JPRINT(57).NE.0) CALL WGTPR(WGTM,MAP,ICNT)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SYMWGT(RHOSTR,WGTP,WGTM,MAP,NRES,IRESET,IDXRES,MAXRES,
     +                  ICNT,NREF,NRTCTL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
C
      DIMENSION RHOSTR(MAXRES,NREF),WGTP(NREF),WGTM(MAXRES*NREF),
     +          MAP(MAXRES*NREF),NRES(NREF),IRESET(NREF),
     +          IDXRES(MAXRES,NREF),NRTCTL(10)
C
      SAVE ZERO,ONE
      SAVE THR1,THR2,THR3
      DATA ZERO,ONE/0.0D0,1.0D0/
      DATA THR1,THR2,THR3/1.0D-4,1.0D-3,1.0D-3/
C
C  Loop over pairs of reference structures checking for symmetry
C  equivalency:
C
      DO 30 IREF = 1,NREF-1
        IF(IRESET(IREF).EQ.0) THEN
          DO 20 JREF = IREF+1,NREF
            IF(IRESET(JREF).EQ.0) THEN
              IF(ABS(RHOSTR(1,IREF)-RHOSTR(1,JREF)).LT.THR1) THEN
                IF(ABS(WGTP(IREF)-WGTP(JREF)).LT.THR2) THEN
                  ITMP = IREF
                  WGTI = -ONE
                  WGTJ = -ONE
                  DO 10 I = 1,ICNT
                    IF(MAP(I).EQ.IDXRES(1,IREF)) WGTI = WGTM(I)
                    IF(MAP(I).EQ.IDXRES(1,JREF)) WGTJ = WGTM(I)
   10             CONTINUE
                  IF(WGTI.LT.ZERO) CALL NBHALT('Error in SR SYMWGT.')
                  IF(WGTJ.LT.ZERO) CALL NBHALT('Error in SR SYMWGT.')
                  IF(ABS(WGTI-WGTJ).GT.THR3) GOTO 100
                END IF
              END IF
            END IF
   20     CONTINUE
        END IF
   30 CONTINUE
      NRTCTL(2) = 0
      RETURN
C
C  Symmetry equivalent structures found of differing weights:
C
  100 NRTCTL(2) = NRTCTL(2) + 1
      IF(NRTCTL(2).GT.1) THEN
        WRITE(LFNPR,910)
        NRTCTL(2) = 0
        RETURN
      END IF
      RHOTMP = RHOSTR(1,ITMP)
      WGTTMP = WGTP(ITMP)
      DO 110 IREF = 1,NREF
        IF(IRESET(IREF).EQ.0) THEN
          IF(ABS(RHOSTR(1,IREF)-RHOTMP).LT.THR1) THEN
            IF(ABS(WGTP(IREF)-WGTTMP).LT.THR2) THEN
              IRESET(IREF) = 1
              NRES(IREF)   = 1
            END IF
          END IF
        END IF
  110 CONTINUE
      DO 120 IREF = 1,NREF
        IF(IRESET(IREF).NE.1) THEN
          IRESET(IREF) = -1
          NRES(IREF)   = 0
        END IF
  120 CONTINUE
      IF(JPRINT(57).NE.0) WRITE(LFNPR,*)
      WRITE(LFNPR,900)
      RETURN
C
  900 FORMAT(17X,'Symmetry broken solution encountered;',
     + ' NRT restart')
  910 FORMAT(/1X,'Warning:  NRT weights may not reflect the ',
     + 'symmetry of the wavefunction.')
      END
C***********************************************************************
      SUBROUTINE NRTOUT(RHOSTR,OCC,WGT,WGTP,SCR,WGTM,Q,MAP,IST,NRES,
     +                  IRESET,IDXRES,LSTRES,MAXRES,ICNT,NREF,LEN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBTOPO/IORDER(MAXATM),JORDER(MAXATM),NTOPO(MAXATM,MAXATM),
     +            N3CTR,I3CTR(10,3)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION RHOSTR(MAXRES,NREF),OCC(NDIM,MAXRES,NREF),
     +          WGT(MAXRES,NREF),WGTP(NREF),SCR(MAXRES),
     +          WGTM(MAXRES*NREF),Q(MAXRES,MAXRES),MAP(MAXRES*NREF),
     +          IST(NDIM,NREF),NRES(NREF),IDXRES(MAXRES,NREF),
     +          LSTRES(LEN),IRESET(NREF)
      CHARACTER*9 OTHERS,MASK
C
      SAVE ZERO,THRESH,ONE,HUNDRD
      SAVE MASK
      DATA ZERO,THRESH,ONE,HUNDRD/0.0D0,1.0D-3,1.0D0,1.0D2/
      DATA MASK/'    -    '/
C
C  Compute the fractional reduction over all NBOs, valence only, and
C  valence+core only:
C
      WRITE(LFNPR,900)
      DO 20 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          CALL COPY(WGT(1,IREF),SCR,MAXRES,NRES(IREF),1)
          WGT(1,IREF) = ONE
          DO 10 IRES = 2,NRES(IREF)
            WGT(IRES,IREF) = ZERO
   10     CONTINUE
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,0)
          CALL GETERR(DR,WGT,Q,NRES,MAXRES,IREF,NREF)
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,1)
          CALL GETERR(DRV,WGT,Q,NRES,MAXRES,IREF,NREF)
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,2)
          CALL GETERR(DRVC,WGT,Q,NRES,MAXRES,IREF,NREF)
          CALL COPY(SCR,WGT(1,IREF),MAXRES,NRES(IREF),1)
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,0)
          CALL GETERR(D,WGT,Q,NRES,MAXRES,IREF,NREF)
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,1)
          CALL GETERR(DV,WGT,Q,NRES,MAXRES,IREF,NREF)
          CALL FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,2)
          CALL GETERR(DVC,WGT,Q,NRES,MAXRES,IREF,NREF)
          FW = ONE - D / DR
          FWV = ONE - DV / DRV
          FWVC = ONE - DVC / DRVC
          WRITE(LFNPR,910) IREF,WGTP(IREF),RHOSTR(1,IREF),DR,FW,FWVC,FWV
        END IF
   20 CONTINUE
C
C  Fetch the TOPO matrix for the leading structure and print:
C
      IP = MAP(1)
      DO 70 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          DO 60 IRES = 1,NRES(IREF)
            IF(IDXRES(IRES,IREF).EQ.IP) THEN
              JRES = IRES
              JREF = IREF
            END IF
   60     CONTINUE
        END IF
   70 CONTINUE
      CALL TOPGET(JRES,IDXRES,MAXRES,JREF,NREF,LSTRES,LEN)
      WRITE(LFNPR,920)
      CALL TOPOUT
C
C  Store the TOPO matrix in NTOPO for future reference:
C
      DO 90 I = 1,NATOMS
        DO 80 J = 1,NATOMS
          NTOPO(J,I) = ITOPO(J,I)
   80   CONTINUE
   90 CONTINUE
C
C  Print resonance weights:
C
      SUM = ZERO
      SUMO = ZERO
      IOTHER = 0
      EPS = THRESH
      IF(JPRINT(57).NE.0) EPS = ZERO
      WRITE(LFNPR,930)
      CALL PRTWGT(0,0,0,W,0,0,0,0,0)
      DO 170 IC = 1,ICNT
        W = WGTM(IC)
        IF(W.GE.EPS) THEN
          IP = MAP(IC)
          NUM = 0
          DO 95 IREF = 1,NREF
            IF(IRESET(IREF).GE.0) THEN
              IF(IDXRES(1,IREF).EQ.IP) THEN
                WGTRES = WGT(1,IREF)*WGTP(IREF)
                NUM = NUM + 1
              END IF
            END IF
   95     CONTINUE
          DO 99 IREF = 1,NREF
            NTMP = 0
            IF(IRESET(IREF).GE.0) THEN
              DO 98 IRES = 1,NRES(IREF)
                IF(IDXRES(IRES,IREF).EQ.IP) THEN
                  WGTRES = WGT(IRES,IREF)*WGTP(IREF)
                  IF(WGTRES.GE.EPS) NTMP = NTMP + 1
                END IF
   98         CONTINUE
            END IF
            NUM = MAX(NUM,NTMP)
   99     CONTINUE
          IR = 0
          DO 110 IREF = 1,NREF
            IF(IRESET(IREF).GE.0) THEN
              DO 100 IRES = 1,NRES(IREF)
                IF(IDXRES(IRES,IREF).EQ.IP) THEN
                  JRES = IRES
                  JREF = IREF
                  IF(JRES.EQ.1) IR = 1
                END IF
  100         CONTINUE
            END IF
  110     CONTINUE
          SUM = SUM + W
          CALL PRTWGT(IC,IR,NUM,W,0,0,0,0,0)
          CALL TOPGET(JRES,IDXRES,MAXRES,JREF,NREF,LSTRES,LEN)
          DO 140 I = 1,NATOMS
            DO 130 J = I+1,NATOMS
              ND = ITOPO(I,J) - NTOPO(I,J)
              DO 120 II = 1,ABS(ND)
                CALL PRTWGT(IC,IR,NUM,W,I,IATNO(I),J,IATNO(J),
     +                      SIGN(1,ND))
  120         CONTINUE
  130       CONTINUE
  140     CONTINUE
          DO 160 I = 1,NATOMS
            ND = ITOPO(I,I) - NTOPO(I,I)
            DO 150 II = 1,ABS(ND)
              CALL PRTWGT(IC,IR,NUM,W,I,IATNO(I),I,IATNO(I),SIGN(1,ND))
  150       CONTINUE
  160     CONTINUE
        ELSE
          IOTHER = IOTHER + 1
          SUMO = SUMO + W
        END IF
  170 CONTINUE
      CALL PRTWGT(-1,0,0,W,0,0,0,0,0)
      IF(IOTHER.EQ.1) THEN
        WRITE(LFNPR,940) ICNT,SUMO*HUNDRD
      ELSE IF(IOTHER.GT.1) THEN
        OTHERS = MASK
        WRITE(OTHERS(1:4),'(I4)') ICNT-IOTHER+1
        IF(ICNT.LT.10) THEN
          WRITE(OTHERS(6:6),'(I1)') ICNT
        ELSE IF(ICNT.LT.100) THEN
          WRITE(OTHERS(6:7),'(I2)') ICNT
        ELSE IF(ICNT.LT.1000) THEN
          WRITE(OTHERS(6:8),'(I3)') ICNT
        ELSE
          WRITE(OTHERS(6:9),'(I4)') ICNT
        END IF
        WRITE(LFNPR,950) OTHERS,SUMO*HUNDRD
      END IF
      WRITE(LFNPR,960) (SUM + SUMO) * HUNDRD
      RETURN
C
  900 FORMAT(/1X,'                                               fract',
     + 'ional accuracy f(w)',/1X,'                 non-Lewis         ',
     + '    -------------------------------------',/1X,' Ref     Wgt ',
     + '     density      d(0)      all NBOs     val+core     valence',
     + /1X,'---------------------------------------------------------',
     + '-------------------')
  910 FORMAT(1X,I3,F11.5,F11.5,F11.5,F12.5,F13.5,F13.5)
  920 FORMAT(//1X,'TOPO matrix for the leading resonance structure:')
  930 FORMAT(/1X,'        Resonance',/1X,'   RS  ',
     + ' Weight(%)                  Added(Removed)',/1X,'-------------',
     + '--------------------------------------------------------------')
  940 FORMAT(1X,I4,5X,F6.2)
  950 FORMAT(1X,A9,F6.2)
  960 FORMAT(1X,'-----------------------------------------------------',
     + '----------------------',/7X,F9.2,3X,'* Total *                ',
     + '[* = reference structure]')
      END
C***********************************************************************
      SUBROUTINE NBDORD(BORDER,VALENZ,TOPO,WGT,XION,WGTP,NRES,IRESET,
     +                  IDXRES,MAXRES,NREF,LSTRES,LEN,NLOW)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*80 TITLE,STRING,BLANKS
      CHARACTER*7  NULL
      CHARACTER*6  DASHES
      CHARACTER*2  CHARAT
      CHARACTER*1  BLANK
      CHARACTER*3  TCI
C
      PARAMETER(MAXATM = 200)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION BORDER(NATOMS,NATOMS),VALENZ(NATOMS),
     +          TOPO(NATOMS,NATOMS),WGT(MAXRES,NREF),
     +          XION(NLOW,MAXRES,NREF),WGTP(NREF),NRES(NREF),
     +          IDXRES(MAXRES,NREF),LSTRES(LEN),IRESET(NREF)
      DIMENSION BOALFA(MAXATM,MAXATM),VALFA(MAXATM,3)
C
      SAVE BOALFA,VALFA
      SAVE ZERO,CUTOFF,PT5,ONE
      SAVE BLANK,DASHES,NULL
      SAVE TCI
C
      DATA ZERO,CUTOFF,PT5,ONE/0.0D0,1.0D-8,0.5D0,1.0D0/
      DATA BLANK,DASHES,NULL/' ','------','   --- '/
      DATA TCI/'tci'/
C
      DO 5 I = 1,80
        BLANKS(I:I) = BLANK
    5 CONTINUE
C
C  Initialize the BORDER and TOPO arrays:
C
      DO 20 IAT = 1,NATOMS
        DO 10 JAT = 1,NATOMS
          BORDER(JAT,IAT) = ZERO
          TOPO(JAT,IAT)   = ZERO
   10   CONTINUE
   20 CONTINUE
C
C  Compute natural bond orders and atomic valencies:
C
      ETA = ONE
      IF(ISPIN.NE.0) ETA = PT5
      DO 60 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          DO 50 IRES = 1,NRES(IREF)
            CALL TOPGET(IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
            DO 40 JAT = 1,NATOMS
              DO 30 IAT = JAT,NATOMS
                BORDER(IAT,JAT) = BORDER(IAT,JAT) +
     +                  WGTP(IREF) * WGT(IRES,IREF) * ITOPO(IAT,JAT)
   30         CONTINUE
   40       CONTINUE
   50     CONTINUE
        END IF
   60 CONTINUE
      DO 80 JAT = 1,NATOMS
        DO 70 IAT = JAT,NATOMS
          BORDER(IAT,JAT) = ETA * BORDER(IAT,JAT)
          BORDER(JAT,IAT) = BORDER(IAT,JAT)
   70   CONTINUE
   80 CONTINUE
      DO 100 IAT = 1,NATOMS
        VALENZ(IAT) = ZERO
        DO 90 JAT = 1,NATOMS
          IF(IAT.NE.JAT) VALENZ(IAT) = VALENZ(IAT) + BORDER(JAT,IAT)
   90   CONTINUE
  100 CONTINUE
C
C  Compute the ionic contributions to the bond order:
C
      DO 140 IREF = 1,NREF
        IF(IRESET(IREF).GE.0) THEN
          DO 130 IRES = 1,NRES(IREF)
            DO 120 IAT = 2,NATOMS
              IPTR = (IAT - 1) * (IAT - 2) / 2
              DO 110 JAT = 1,IAT-1
                JPTR = IPTR + JAT
                TOPO(IAT,JAT) = TOPO(IAT,JAT) +
     +               WGTP(IREF) * WGT(IRES,IREF) * XION(JPTR,IRES,IREF)
  110         CONTINUE
  120       CONTINUE
  130     CONTINUE
        END IF
  140 CONTINUE
      DO 160 IAT = 2,NATOMS
        DO 150 JAT = 1,IAT-1
          TOPO(IAT,JAT) = ETA * TOPO(IAT,JAT)
          TOPO(JAT,IAT) = TOPO(IAT,JAT)
  150   CONTINUE
  160 CONTINUE
C
C  Compute the fractional ionic character:
C
      DO 180 IAT = 1,NATOMS-1
        DO 170 JAT = IAT+1,NATOMS
          IF(BORDER(JAT,IAT).GT.CUTOFF) THEN
            TOPO(JAT,IAT) = TOPO(JAT,IAT) / BORDER(JAT,IAT)
          ELSE
            TOPO(JAT,IAT) = ZERO
          END IF
          TOPO(IAT,JAT) = TOPO(JAT,IAT)
  170   CONTINUE
  180 CONTINUE
C
C  Print out results:
C
      TITLE = 'Natural Bond Order:  (total/covalent/ionic)'
      WRITE(LFNPR,900) TITLE
      DO 200 IL = 1,NATOMS,9
        IU = MIN(NATOMS,IL+8)
        WRITE(LFNPR,940) (I,I=IL,IU)
        WRITE(LFNPR,950) (DASHES,I=IL,IU)
        DO 190 IAT = 1,NATOMS
          STRING = BLANKS
          IF(IAT.NE.1) WRITE(LFNPR,*)
          WRITE(LFNPR,960) IAT,CHARAT(IATNO(IAT)),TCI(1:1),
     +                     (BORDER(I,IAT),I=IL,IU)
          STRING = BLANKS
          WRITE(STRING,'(7X,2X,A1,9F7.4)') TCI(2:2),
     +                     ((BORDER(I,IAT)*(ONE-TOPO(I,IAT))),I=IL,IU)
          IF(IAT.GE.IL.AND.IAT.LE.IU) THEN
            IC = 7*(IAT-IL+1) + 4
            STRING(IC:IC+6) = NULL
          END IF
          WRITE(LFNPR,970) STRING
          STRING = BLANKS
          WRITE(STRING,'(7X,2X,A1,9F7.4)') TCI(3:3),
     +                     ((BORDER(I,IAT)*TOPO(I,IAT)),I=IL,IU)
          IF(IAT.GE.IL.AND.IAT.LE.IU) THEN
            IC = 7*(IAT-IL+1) + 4
            STRING(IC:IC+6) = NULL
          END IF
          WRITE(LFNPR,970) STRING
  190  CONTINUE
  200 CONTINUE
C
C  Compute and print valencies (save alpha-spin valencies, sum with beta
C  spin valencies):
C
      TITLE = 'Natural Atomic Valencies:'
      WRITE(LFNPR,900) TITLE
      WRITE(LFNPR,910)
      DO 220 IAT = 1,NATOMS
        EVAL = ZERO
        DO 210 JAT = 1,NATOMS
          IF(IAT.NE.JAT) THEN
            EVAL = EVAL + BORDER(JAT,IAT) * TOPO(JAT,IAT)
          END IF
  210   CONTINUE
        CVAL = VALENZ(IAT) - EVAL
        WRITE(LFNPR,920) IAT,NAMEAT(IATNO(IAT)),VALENZ(IAT),CVAL,EVAL
        IF(ISPIN.EQ.2) THEN
          VALFA(IAT,1) = VALENZ(IAT)
          VALFA(IAT,2) = CVAL
          VALFA(IAT,3) = EVAL
        ELSE IF(ISPIN.EQ.-2) THEN
          VALFA(IAT,1) = VALFA(IAT,1) + VALENZ(IAT)
          VALFA(IAT,2) = VALFA(IAT,2) + CVAL
          VALFA(IAT,3) = VALFA(IAT,3) + EVAL
        END IF
  220 CONTINUE
C
C  If alpha-spin, save the bond orders:
C
      IF(ISPIN.EQ.2) THEN
        DO 240 JAT = 1,NATOMS
          DO 230 IAT = 1,NATOMS
            BOALFA(IAT,JAT) = BORDER(IAT,JAT)
  230     CONTINUE
  240   CONTINUE
      END IF
C
C  If beta-spin, compute the total bond orders:
C
      IF(ISPIN.EQ.-2) THEN
        DO 260 JAT = 1,NATOMS
          DO 250 IAT = 1,NATOMS
            BOALFA(IAT,JAT) = BOALFA(IAT,JAT) + BORDER(IAT,JAT)
  250     CONTINUE
  260   CONTINUE
C
C  Print total bond orders and valencies:
C
        TITLE = 'Natural Bond Order (total):'
        CALL AOUT(BOALFA,MAXATM,NATOMS,NATOMS,TITLE,0,NATOMS)
        TITLE = 'Natural Atomic Valencies (total):'
        WRITE(LFNPR,900) TITLE
        WRITE(LFNPR,910)
        DO 270 IAT = 1,NATOMS
          WRITE(LFNPR,920) IAT,NAMEAT(IATNO(IAT)),(VALFA(IAT,I),I=1,3)
  270   CONTINUE
      END IF
      RETURN
C
  900 FORMAT(//1X,A78)
  910 FORMAT(/1X,'                     Co-    Electro-',
     +       /1X,'    Atom  Valency  Valency  Valency',
     +       /1X,'    ----  -------  -------  -------')
  920 FORMAT(1X,I3,'. ',A2,1X,3F9.4)
  940 FORMAT(/1X,4X,'Atom',3X,I5,8(2X,I5))
  950 FORMAT(5X,'----',3X,A6,8(1X,A6))
  960 FORMAT(1X,I3,'. ',A2,2X,A1,9F7.4)
  970 FORMAT(1X,A73)
      END
C***********************************************************************
      SUBROUTINE NRTLST(MAP,WGTM,ICNT,NRES,IRESET,IDXRES,MAXRES,NREF,
     +                  LSTRES,LEN)
C***********************************************************************
C 24-May-07  FAW  Added 'P','H' TYPEs for pentuple, hextuple bonds
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
C
      DIMENSION MAP(ICNT),WGTM(ICNT),IRESET(NREF),IDXRES(MAXRES,NREF),
     +          LSTRES(LEN),NRES(NREF)
      CHARACTER*78 STRING
      CHARACTER*8 LONE,BOND,BLANKS
      CHARACTER*4 TERM
      CHARACTER*2 TYPE(6)
C
      SAVE LONE,BOND,BLANKS,TERM,TYPE
      SAVE ZERO,HUNDRD
      DATA LONE,BOND,BLANKS/'    LONE','    BOND','        '/
      DATA TERM/' END'/
      DATA TYPE/' S',' D',' T',' Q',' P',' H'/
      DATA ZERO,HUNDRD/0.0D0,1.0D2/
C
      IF(ISPIN.EQ.0)  WRITE(LFNPR,900)
      IF(ISPIN.EQ.2)  WRITE(LFNPR,910)
      IF(ISPIN.EQ.-2) WRITE(LFNPR,920)
      DO 90 I = 1,ICNT
        IP = MAP(I)
        IF(WGTM(I).GT.PRTHR) THEN
          JRES = 0
          JREF = 0
          IF(PRTHR.LT.ZERO) THEN
            DO 10 IREF = 1,NREF
              IF(IRESET(IREF).GE.0) THEN
                IF(IP.EQ.IDXRES(1,IREF)) THEN
                  JREF = IREF
                  JRES = 1
                  GOTO 20
                END IF
              END IF
   10       CONTINUE
            GOTO 90
          ELSE
            DO 15 IREF = 1,NREF
              IF(IRESET(IREF).GE.0) THEN
                DO 14 IRES = 1,NRES(IREF)
                  IF(IP.EQ.IDXRES(IRES,IREF)) THEN
                    JREF = IREF
                    JRES = IRES
                    GOTO 20
                  END IF
   14           CONTINUE
              END IF
   15       CONTINUE
          END IF
   20     CONTINUE
          IF(JRES.EQ.0) CALL NBHALT('No structure found in NRTLST.')
C
C  Punch this structure:
C
          CALL TOPGET(JRES,IDXRES,MAXRES,JREF,NREF,LSTRES,LEN)
          WRITE(LFNPR,930) WGTM(I)*HUNDRD
C
C Write lone pair specifications:
C
          ILONE = 0
          DO 30 IAT = 1,NATOMS
            IF(ITOPO(IAT,IAT).NE.0) ILONE = 1
   30     CONTINUE
          IF(ILONE.NE.0) THEN
            N = 8
            STRING(1:N) = LONE
            DO 40 IAT = 1,NATOMS
              IT = ITOPO(IAT,IAT)
              IF(IT.NE.0) THEN
                IF(N+6.GT.78) THEN
                  WRITE(LFNPR,940) STRING(1:N)
                  N = 8
                  STRING(1:N) = BLANKS
                END IF
                IF(IAT.LT.10) THEN
                  N = N + 2
                  WRITE(STRING(N-1:N),'(I2)') IAT
                ELSE IF(IAT.LT.100) THEN
                  N = N + 3
                  WRITE(STRING(N-2:N),'(I3)') IAT
                ELSE
                  N = N + 4
                  WRITE(STRING(N-3:N),'(I4)') IAT
                END IF
                N = N + 2
                WRITE(STRING(N-1:N),'(I2)') IT
              END IF
   40       CONTINUE
            IF(N+4.GT.78) THEN
              WRITE(LFNPR,940) STRING(1:N)
              N = 8
              STRING(1:N) = BLANKS
            END IF
            N = N + 4
            STRING(N-3:N) = TERM
            WRITE(LFNPR,940) STRING(1:N)
          END IF
C
C Write bond specifications:
C
          IBOND = 0
          DO 60 IAT = 1,NATOMS-1
            DO 50 JAT = IAT+1,NATOMS
              IF(ITOPO(JAT,IAT).NE.0) IBOND = 1
   50       CONTINUE
   60     CONTINUE
          IF(IBOND.NE.0) THEN
            N = 8
            STRING(1:N) = BOND
            DO 80 IAT = 1,NATOMS-1
              DO 70 JAT = IAT+1,NATOMS
                IT = ITOPO(JAT,IAT)
                IF(IT.NE.0) THEN
                  IF(N+10.GT.78) THEN
                    WRITE(LFNPR,940) STRING(1:N)
                    N = 8
                    STRING(1:N) = BLANKS
                  END IF
                  N = N + 2
                  WRITE(STRING(N-1:N),'(A2)') TYPE(IT)
                  IF(IAT.LT.10) THEN
                    N = N + 2
                    WRITE(STRING(N-1:N),'(I2)') IAT
                  ELSE IF(IAT.LT.100) THEN
                    N = N + 3
                    WRITE(STRING(N-2:N),'(I3)') IAT
                  ELSE
                    N = N + 4
                    WRITE(STRING(N-3:N),'(I4)') IAT
                  END IF
                  IF(JAT.LT.10) THEN
                    N = N + 2
                    WRITE(STRING(N-1:N),'(I2)') JAT
                  ELSE IF(JAT.LT.100) THEN
                    N = N + 3
                    WRITE(STRING(N-2:N),'(I3)') JAT
                  ELSE
                    N = N + 4
                    WRITE(STRING(N-3:N),'(I4)') JAT
                  END IF
                END IF
   70         CONTINUE
   80       CONTINUE
            IF(N+4.GT.78) THEN
              WRITE(LFNPR,940) STRING(1:N)
              N = 8
              STRING(1:N) = BLANKS
            END IF
            N = N + 4
            STRING(N-3:N) = TERM
            WRITE(LFNPR,940) STRING(1:N)
          END IF
          WRITE(LFNPR,950)
        END IF
   90 CONTINUE
      WRITE(LFNPR,960)
      RETURN
C
  900 FORMAT(//1X,'$NRTSTR')
  910 FORMAT(/1X,'$NRTSTRA')
  920 FORMAT(/1X,'$NRTSTRB')
  930 FORMAT(1X,'  STR        ! Wgt =',F6.2,'%')
  940 FORMAT(1X,A)
  950 FORMAT(1X,'  END')
  960 FORMAT(1X,'$END')
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NRT OPTIMIZATION ROUTINES
C
C      SUBROUTINE SETOPT(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE BFGS(FRET,FREF,OCC,WGT,Q,P,HESS,XT,XIT,G,DG,HDG,DF,W,
C     +                SCR,IPTR,NRES,MAXRES,IREF,NREF,NVAR,ITER)
C      SUBROUTINE NBPWLL(FRET,FREF,WGT,Q,P,XI,PT,XIT,PTT,W,WP,IPTR,
C     +                  NRES,MAXRES,IREF,NREF,NVAR,ITER)
C      SUBROUTINE NLINMN(FRET,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
C     +                  NREF,NVAR)
C      SUBROUTINE DLINMN(FRET,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,
C     +                  MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE NBMNBK(AX,BX,CX,FA,FB,FC,WGT,Q,XT,P,XIT,IPTR,
C     +                  NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE NBBRNT(FRET,AX,BX,CX,FB,TOL,XMIN,WGT,Q,XT,PT,XIT,IPTR,
C     +                  NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE DBRENT(FRET,AX,BX,CX,FB,TOL,XMIN,OCC,WGT,Q,XT,P,XIT,DF,
C     +                  SCR,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE NBFDIM(FX,X,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
C     +                 NREF,NVAR)
C      SUBROUTINE DF1DIM(DX,X,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,
C     +                  MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE NBPHI(FT,T,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE NBDPHI(DER,T,OCC,WGT,Q,SCR,IPTR,NRES,MAXRES,IREF,
C     +                  NREF,NVAR)
C
C***********************************************************************
      SUBROUTINE SETOPT(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION RHOSTR(MAXRES,NREF),WGT(MAXRES,NREF),T(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,DELTA,EPS,ONE,HUNDRD
      DATA ZERO,DELTA,EPS,ONE,HUNDRD/0.0D0,1.0D-4,1.0D-5,1.0D0,1.0D2/
C
C  Setup pointer array:
C
      DO 10 IRES = 1,NRES(IREF)
        IPTR(IRES) = -1
   10 CONTINUE
C
C  Do not optimize the weight of the reference structure that we're
C  working with.  Its weight will be determined by normalization:
C  (T is the array of parameters to be optimized.)
C
      IPTR(1) = 0
      WGT(1,IREF) = -ONE
C
C  Optimize the weight of the rest of the structures, constraining
C  resonance weights with equivalent rho* to have equal weights:
C
      NVAR = 0
   30 JRES = 0
      RHO = HUNDRD
      DO 40 IRES = 1,NRES(IREF)
        IF(RHOSTR(IRES,IREF).GE.ZERO.AND.RHOSTR(IRES,IREF).LT.RHO) THEN
          IF(WGT(IRES,IREF).GE.ZERO) THEN
            JRES = IRES
            RHO = RHOSTR(JRES,IREF)
          END IF
        END IF
   40 CONTINUE
C
C  Compute the optimizable parameter T and set pointers appropriately:
C
      IF(JRES.NE.0) THEN
        NVAR = NVAR + 1
        IF(WGT(JRES,IREF).GT.DELTA) THEN
          T(NVAR) = -LOG(WGT(JRES,IREF))
        ELSE
          T(NVAR) = -LOG(DELTA)
        END IF
        DO 50 IRES = 1,NRES(IREF)
          IF(ABS(RHOSTR(IRES,IREF)-RHO).LT.EPS) THEN
            IF(WGT(IRES,IREF).GE.ZERO) THEN
              IPTR(IRES) = NVAR
              WGT(IRES,IREF) = -ONE
            END IF
          END IF
   50   CONTINUE
        GOTO 30
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE BFGS(FRET,FREF,OCC,WGT,Q,P,HESS,XT,XIT,G,DG,HDG,DF,W,
     +                SCR,IPTR,NRES,MAXRES,IREF,NREF,NVAR,ITER)
C***********************************************************************
C
C  Modified version of DFPMIN Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 307-311.
C
C  Given a starting point P that is a vector of length N, the Broyden-
C  Fletcher-Goldfarb-Shanno variant of Davidon-Fletcher-Powell minimization
C  is performed on a subroutine NBPHI, using its gradient as calculated by
C  routine NBDPHI.  The convergence requirement on the function value is
C  input as TOL.  Returned quantities are P (the location on the minimum),
C  ITER (the number of iterations that were performed), and FRET (the
C  minimum value of the function).
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ITMAX = 200,EPS = 1.0D-10)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION OCC(NDIM,MAXRES,NREF),WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          P(MAXRES),HESS(MAXRES,MAXRES),XT(MAXRES),XIT(MAXRES),
     +          G(MAXRES),DG(MAXRES),HDG(MAXRES),DF(MAXRES),W(MAXRES),
     +          SCR(MAXRES),IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,TOL,ONE,TWO
      DATA ZERO,TOL,ONE,TWO/0.0D0,1.0D-6,1.0D0,2.0D0/
C
C  Calculate starting function value and gradient:
C
      CALL NBPHI(FP,P,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
      CALL NBDPHI(G,P,OCC,WGT,Q,SCR,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  Initialize the inverse Hessian to the unit matrix:
C
      DO 20 I = 1,NVAR
        DO 10 J = 1,NVAR
          HESS(I,J) = ZERO
   10   CONTINUE
        HESS(I,I) = ONE
        XIT(I) = -G(I)
   20 CONTINUE
C
C  Main loop over the iterations
C
      IF(JPRINT(57).NE.0) WRITE(LFNPR,900)
      DO 160 ITER = 1,ITMAX
C
C  Store weights for future reference:
C
        DO 30 IRES = 1,NRES(IREF)
          W(IRES) = WGT(IRES,IREF)
   30   CONTINUE
C
C  Minimize along the gradient vector:
C
        CALL DLINMN(FRET,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,MAXRES,
     +              IREF,NREF,NVAR)
C
C  Find the maximum and rms gradients:
C
        GMAX = ZERO
        GRMS = ZERO
        DO 40 IV = 1,NVAR
          IF(ABS(G(IV)).GT.ABS(GMAX)) GMAX = G(IV)
          GRMS = GRMS + G(IV) * G(IV)
   40   CONTINUE
        GRMS = SQRT(GRMS / NVAR)
C
C  Also, compute the maximum and rms step sizes:
C
        SMAX = ZERO
        SRMS = ZERO
        DO 50 IRES = 1,NRES(IREF)
          STEP = WGT(IRES,IREF) - W(IRES)
          IF(ABS(STEP).GT.ABS(SMAX)) SMAX = STEP
          SRMS = SRMS + STEP * STEP
   50   CONTINUE
        SRMS = SQRT(SRMS / NRES(IREF))
C
C  Print information about convergence:
C
        FW = ONE - FP / FREF
        IF(JPRINT(57).NE.0) WRITE(LFNPR,910) ITER,FW,GMAX,GRMS,SMAX,SRMS
C
C  Check for convergence of FRET, the gradients, and step sizes:
C
        IF(TWO*ABS(FRET-FP).LE.TOL*(ABS(FRET)+ABS(FP)+EPS)) RETURN
C
C  Save the old function value and gradient:
C
        FP = FRET
        DO 60 I = 1,NVAR
          DG(I) = G(I)
   60   CONTINUE
C
C  Get the new function value and gradient:
C
        CALL NBPHI(FRET,P,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
        CALL NBDPHI(G,P,OCC,WGT,Q,SCR,IPTR,NRES,MAXRES,IREF,NREF,
     +              NVAR)
C
C  Compute the difference in gradients:
C
        DO 70 I = 1,NVAR
          DG(I) = G(I) - DG(I)
   70   CONTINUE
C
C  Multiply the current Hessian by the gradient differences:
C
        DO 90 I = 1,NVAR
          HDG(I) = ZERO
          DO 80 J = 1,NVAR
            HDG(I) = HDG(I) + HESS(I,J) * DG(J)
   80     CONTINUE
   90   CONTINUE
C
C  Compute dot products for the denominators:
C
        FAC = ZERO
        FAE = ZERO
        DO 100 I = 1,NVAR
          FAC = FAC + DG(I) * XIT(I)
          FAE = FAE + DG(I) * HDG(I)
  100   CONTINUE
C
C  Make the denominators multiplicative:
C
        FAC = ONE / FAC
        FAD = ONE / FAE
C
C  This vector makes the BFGS different from the DFP optimization scheme:
C
        DO 110 I = 1,NVAR
          DG(I) = FAC * XIT(I) - FAD * HDG(I)
  110   CONTINUE
C
C  The BFGS updating formula:
C
        DO 130 I = 1,NVAR
          DO 120 J = 1,NVAR
            HESS(I,J) = HESS(I,J) + FAC * XIT(I) * XIT(J)
     +        - FAD * HDG(I) * HDG(J) + FAE * DG(I) * DG(J)
  120     CONTINUE
  130   CONTINUE
C
C  Now, determine the direction to go next:
C
        DO 150 I = 1,NVAR
          XIT(I) = ZERO
          DO 140 J = 1,NVAR
            XIT(I) = XIT(I) - HESS(I,J) * G(J)
  140     CONTINUE
  150   CONTINUE
  160 CONTINUE
      ITER = -ITMAX
      RETURN
C
  900 FORMAT(/1X,'  iter      fw     max grad   rms grad   max step',
     + '   rms step',/1X,'-----------------------------------------',
     + '---------------------')
  910 FORMAT(1X,I5,5F11.5)
      END
C***********************************************************************
      SUBROUTINE NBPWLL(FRET,FREF,WGT,Q,P,XI,PT,XIT,PTT,W,WP,IPTR,
     +                  NRES,MAXRES,IREF,NREF,NVAR,ITER)
C***********************************************************************
C 23-Mar-2006  EDG  Name change from POWELL to avoid GAMESS conflict
C-----------------------------------------------------------------------
C
C  Modified version of POWELL, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 294-300.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          P(MAXRES),XI(MAXRES,MAXRES),PT(MAXRES),XIT(MAXRES),
     +          PTT(MAXRES),W(MAXRES),WP(MAXRES),IPTR(MAXRES),
     +          NRES(NREF)
C
      SAVE ZERO,TOL,ONE,TWO
      SAVE ITMAX
      DATA ZERO,TOL,ONE,TWO/0.0D0,1.0D-5,1.0D0,2.0D0/
      DATA ITMAX/100/
C
C  Compute PHI for the initial point P:
C
      CALL NBPHI(FRET,P,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  Initialize XI:
C
      DO 20 J = 1,NVAR
        DO 10 I = 1,NVAR
          XI(I,J) = ZERO
   10   CONTINUE
        XI(J,J) = ONE
   20 CONTINUE
C
C  Save the initial point:
C
      DO 30 J = 1,NVAR
        PT(J) = P(J)
   30 CONTINUE
C
C  Loop until the resonance weights are optimized:
C
      ITER = 0
      IF(JPRINT(57).NE.0) WRITE(LFNPR,900)
   40 ITER = ITER + 1
      FP = FRET
      IBIG = 0
C
C  Save the current resonance weights:
C
      DO 45 IRES = 1,NRES(IREF)
        W(IRES) = WGT(IRES,IREF)
   45 CONTINUE
C
C  DEL will be the largest function decrease:
C
      DEL = ZERO
C
C  In each iteration, loop over all directions in the set, XI:
C
      DO 60 I = 1,NVAR
C
C  Copy the direction:
C
        DO 50 J = 1,NVAR
          XIT(J) = XI(J,I)
   50   CONTINUE
C
C  Minimize along the direction:
C
        FPTT = FRET
        CALL NLINMN(FRET,WGT,Q,PTT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +              NREF,NVAR)
C
C  Record the direction if it is the largest so far:
C
        IF(ABS(FPTT-FRET).GT.DEL) THEN
          DEL = ABS(FPTT - FRET)
          IBIG = I
          DO 55 IRES = 1,NRES(IREF)
            WP(IRES) = WGT(IRES,IREF)
   55     CONTINUE
        END IF
   60 CONTINUE
C
C  Compute the maximum and rms step sizes:
C
        SMAX = ZERO
        SRMS = ZERO
        DO 65 IRES = 1,NRES(IREF)
          STEP = WP(IRES) - W(IRES)
          IF(ABS(STEP).GT.ABS(SMAX)) SMAX = STEP
          SRMS = SRMS + STEP * STEP
   65   CONTINUE
        SRMS = SQRT(SRMS / NRES(IREF))
C
C  Print information about convergence:
C
        FW = ONE - FP / FREF
        IF(JPRINT(57).NE.0) WRITE(LFNPR,910) ITER,FW,SMAX,SRMS
C
C  Should we terminate the optimization:
C
      IF(TWO*ABS(FP-FRET).LE.TOL*(ABS(FP)+ABS(FRET))) RETURN
C
C  Have we taken too many iterations:
C
      IF(ITER.EQ.ITMAX) THEN
        ITER = -ITMAX
        RETURN
      END IF
C
C  Construct the extrapolated point and average direction moved.
C  Also, save the old starting point:
C
      DO 70 J = 1,NVAR
        PTT(J) = TWO * P(J) - PT(J)
        XIT(J) = P(J) - PT(J)
        PT(J) = P(J)
   70 CONTINUE
C
C  Compute the function value at the extrapolated point:
C
      CALL NBPHI(FPTT,PTT,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  Should this new point (or direction) be used:
C
      IF(FPTT.GE.FP) GOTO 40
      T = TWO * (FP - TWO * FRET + FPTT) * (FP - FRET - DEL)**2 -
     +    DEL * (FP - FPTT)**2
      IF(T.GE.ZERO) GOTO 40
C
C  Move to the minimum in this direction:
C
      CALL NLINMN(FRET,WGT,Q,PTT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +            NREF,NVAR)
C
C  And save the new direction:
C
      DO 80 J = 1,NVAR
        XI(J,IBIG) = XIT(J)
   80 CONTINUE
      GO TO 40
C
  900 FORMAT(/1X,'  iter      fw     max step   rms step',
     + /1X,'----------------------------------------')
  910 FORMAT(1X,I5,3F11.5)
      END
C***********************************************************************
      SUBROUTINE NLINMN(FRET,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,
     +                  IREF,NREF,NVAR)
C***********************************************************************
C
C  Modified version of LINMIN, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 300-301.
C
C  Given an N dimensional point P and an N dimensional direction XIT, moves
C  and resets P to where the function FUNC(P) takes on a minimum along the
C  direction XIT from P, and replaces XIT by the actual vector displacement
C  the P moved.  Also, returns as FRET the value of FUNC at the new location
C  P.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          XT(MAXRES),P(MAXRES),XIT(MAXRES),IPTR(MAXRES),
     +          NRES(NREF)
C
      SAVE ZERO,TOL
      DATA ZERO,TOL/0.0D0,1.0D-3/
C
C  Initial guess for brackets:
C
      AX = -ZERO
      XX = -1.0D-2
      CALL NBMNBK(AX,XX,BX,FA,FX,FB,WGT,Q,XT,P,XIT,IPTR,NRES,
     +            MAXRES,IREF,NREF,NVAR)
      CALL NBBRNT(FRET,AX,XX,BX,FX,TOL,XMIN,WGT,Q,XT,P,XIT,IPTR,
     +           NRES,MAXRES,IREF,NREF,NVAR)
C
C  Contruct vector results to return:
C
      DO 10 J = 1,NVAR
        XIT(J) = XMIN * XIT(J)
        P(J) = P(J) + XIT(J)
   10 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE DLINMN(FRET,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,
     +                  MAXRES,IREF,NREF,NVAR)
C***********************************************************************
C
C  Modified version of LINMIN, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 300-301,306.
C
C  Given an N dimensional point P and an N dimensional direction XIT, moves
C  and resets P to where the function FUNC(P) takes on a minimum along the
C  direction XIT from P, and replaces XIT by the actual vector displacement
C  the P moved.  Also, returns as FRET the value of FUNC at the new location
C  P.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION OCC(NDIM,MAXRES,NREF),WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          XT(MAXRES),P(MAXRES),XIT(MAXRES),DF(MAXRES),SCR(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,TOL
      DATA ZERO,TOL/0.0D0,1.0D-3/
C
C  Initial guess for brackets:
C
      AX =  ZERO
      XX =  0.1D0
      CALL NBMNBK(AX,XX,BX,FA,FX,FB,WGT,Q,XT,P,XIT,IPTR,NRES,
     +            MAXRES,IREF,NREF,NVAR)
      CALL DBRENT(FRET,AX,XX,BX,FX,TOL,XMIN,OCC,WGT,Q,XT,P,XIT,DF,SCR,
     +            IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  Contruct vector results to return:
C
      DO 10 J = 1,NVAR
        XIT(J) = XMIN * XIT(J)
        P(J) = P(J) + XIT(J)
   10 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBMNBK(AX,BX,CX,FA,FB,FC,WGT,QQ,XT,P,XIT,IPTR,
     +                  NRES,MAXRES,IREF,NREF,NVAR)
C***********************************************************************
C 23-Mar-2006  EDG  Name change from MNBRAK to avoid GAMESS conflict
C-----------------------------------------------------------------------
C
C  Modified version of MNBRAK, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 277-282.
C
C  Given distinct initial points AX and BX, this routine searches in the
C  downhill direction and returns new points AX, BX, and CX which bracket
C  the minimum of the function NBFDIM.  Also returned to the calling routine
C  are the function values at the three points, FA, FB, and FC.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(GOLD = 1.618034D0,GLIMIT = 10.0D0, TINY=1.0D-10)
C
      DIMENSION WGT(MAXRES,NREF),
     +          QQ(MAXRES,MAXRES),XT(MAXRES),P(MAXRES),XIT(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,TWO
      DATA ZERO,TWO/0.0D0,2.0D0/
C
      CALL NBFDIM(FA,AX,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +            NREF,NVAR)
      CALL NBFDIM(FB,BX,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +            NREF,NVAR)
C
C  Switch roles of A and B (if necessary) so that we go down hill in
C  the direction from A to B:
C
      IF(FB.GT.FA) THEN
        DUM = AX
        AX  = BX
        BX  = DUM
        DUM = FB
        FB  = FA
        FA  = DUM
      END IF
C
C  First guess for C:
C
      CX = BX + GOLD * (BX - AX)
      CALL NBFDIM(FC,CX,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +            NREF,NVAR)
C
C  Keep returning here until we bracket the solution:
C
   10 IF(FB.GT.FC) THEN
C
C  Compute U by parabolic extrapolation from A,B,C.  TINY is used to
C  prevent any possible division by zero:
C
        R = (BX - AX) * (FB - FC)
        Q = (BX - CX) * (FB - FA)
        U = BX - ((BX - CX) * Q - (BX - AX) * R)
        U = U / (TWO * SIGN(MAX(ABS(Q - R),TINY),Q - R))
C
C  Watch the step size:
C
        ULIM = BX + GLIMIT * (CX - BX)
C
C  Try a parabolic fit between B and C:
C
        IF((BX-U)*(U-CX).GT.ZERO) THEN
          CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
          IF(FU.LT.FC) THEN
            AX = BX
            FA = FB
            BX = U
            FB = FU
            GOTO 10
C
C  We have a minimum between A and U:
C
          ELSE IF(FU.GT.FB) THEN
            CX = U
            FC = FU
            GOTO 10
          END IF
C
C  Otherwise the parabolic fit is useless.  Try default magnification:
C
          U = CX + GOLD * (CX - BX)
          CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
C
C  Parabolic fit between C and its allowed limit:
C
        ELSE IF((CX-U)*(U-ULIM).GT.ZERO) THEN
          CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
          IF(FU.LT.FC) THEN
            BX = CX
            CX = U
            U  = CX + GOLD * (CX - BX)
            FB = FC
            FC = FU
            CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,
     +                  IREF,NREF,NVAR)
          END IF
C
C  Limit parabolic U to maximum allowed value:
C
        ELSE IF((U-ULIM)*(ULIM-CX).GE.ZERO) THEN
          U  = ULIM
          CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,
     +                IREF,NREF,NVAR)
C
C  Reject parabolic U, use default magnification:
C
        ELSE
          U  = CX + GOLD * (CX - BX)
          CALL NBFDIM(FU,U,WGT,QQ,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
        END IF
C
C  Eliminate the oldest point an continue:
C
        AX = BX
        BX = CX
        CX = U
        FA = FB
        FB = FC
        FC = FU
        GOTO 10
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBBRNT(FRET,AX,BX,CX,FB,TOL,XMIN,WGT,QQ,XT,PT,XIT,
     +                 IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C***********************************************************************
C 23-Mar-06  EDG  Name change from BRENT to avoid GAMESS conflict
C 22-Jun-95  EDG  BRENT changed from FUNCTION to SUBROUTINE
C-----------------------------------------------------------------------
C
C  Modified version of BRENT, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), pp. 283-286.
C
C  Given a function NBFDIM, and given a bracketing triplet of abscissas
C  (such that BX is between AX and CX, and NBFDIM(BX) is less than both
C  NBFDIM(AX) and NBFDIM(CX)), this routine isolates the minimum to a
C  fractional precision of about TOL using Brent's method.  The abscissa
C  of the minimum is returned as XMIN, and the minimum function value is
C  returned as FRET.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(ITMAX = 100,CGOLD = 0.3819660D0,ZEPS = 1.0D-10)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION WGT(MAXRES,NREF),
     +          QQ(MAXRES,MAXRES),XT(MAXRES),PT(MAXRES),XIT(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,PT5,TWO
      DATA ZERO,PT5,TWO/0.0D0,0.5D0,2.0D0/
C
C  Initializations:
C
      A = MIN(AX,CX)
      B = MAX(AX,CX)
      V = BX
      W = BX
      X = BX
      D = ZERO
      E = ZERO
      FX = FB
      FV = FX
      FW = FX
C
C  Main program loop:
C
      DO 30 ITER = 1,ITMAX
        XM = PT5 * (A + B)
        TOL1 = TOL * ABS(X) + ZEPS
        TOL2 = TWO * TOL1
C
C  Determine whether we are finished:
C
        IF(ABS(X-XM).LE.(TOL2-PT5*(B-A))) GOTO 40
C
C  Construct a trial parabolic fit:
C
        IF(ABS(E).GT.TOL1) THEN
          R = (X - W) * (FX - FV)
          Q = (X - V) * (FX - FW)
          P = (X - V) * Q - (X - W) * R
          Q = TWO * (Q - R)
          IF(Q.GT.ZERO) P = -P
          Q = ABS(Q)
          ETEMP = E
          E = D
C
C  Is the parabolic fit OK:
C
          IF(ABS(P).GE.ABS(PT5*Q*ETEMP).OR.P.LE.Q*(A-X).OR.
     +        P.GE.Q*(B-X)) GOTO 10
C
C  Take a parabolic step:
C
          D = P / Q
          U = X + D
          IF(U-A.LT.TOL2.OR.B-U.LT.TOL2) D = SIGN(TOL1,XM - X)
          GOTO 20
        END IF
C
C  Take a golden section step, rather than parabolic:
C
   10   IF(X.GE.XM) THEN
          E = A - X
        ELSE
          E = B - X
        END IF
        D = CGOLD * E
C
C  Step D has been taken:
C
   20   IF(ABS(D).GE.TOL1) THEN
          U = X + D
        ELSE
          U = X + SIGN(TOL1,D)
        END IF
C
C  Evaluate the function at the point U:
C
        CALL NBFDIM(FU,U,WGT,QQ,XT,PT,XIT,IPTR,NRES,MAXRES,IREF,
     +              NREF,NVAR)
C
C  Do housekeeping:
C
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A = X
          ELSE
            B = X
          END IF
          V  = W
          FV = FW
          W  = X
          FW = FX
          X  = U
          FX = FU
        ELSE
          IF(U.LT.X) THEN
            A = U
          ELSE
            B = U
          END IF
          IF(FU.LE.FW.OR.W.EQ.X) THEN
            V  = W
            FV = FW
            W  = U
            FW = FU
          ELSE IF(FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
            V  = U
            FV = FU
          END IF
        END IF
   30 CONTINUE
      WRITE(LFNPR,900)
C
C  Ready to exit with best value:
C
   40 XMIN = X
      FRET = FX
      RETURN
C
  900 FORMAT(1X,'NBBRNT exceeded maximum iterations.')
      END
C***********************************************************************
      SUBROUTINE DBRENT(FRET,AX,BX,CX,FB,TOL,XMIN,OCC,WGT,Q,XT,P,XIT,DF,
     +                  SCR,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C***********************************************************************
C 22-Jun-95  EDG  DBRENT changed from FUNCTION to SUBROUTINE
C-----------------------------------------------------------------------
C
C  Modified version of DBRENT Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), p. 286-289.
C
C  Given a bracketing triplet of abscissas AX, BX, CX [such that BX is
C  between AX and CX, and NBFDIM(BX) is less than both NBFDIM(AX) and
C  NBFDIM(CX)], this routine isolates the minimum to a fractional precision
C  of about TOL using a modification of Brent's method that uses derivative
C  information.  The abscissa of the minimum is returned as XMIN, and the
C  minimum function value is returned as FRET.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL OK1,OK2
C
      PARAMETER(ITMAX = 100,ZEPS = 1.0D-10)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION OCC(NDIM,MAXRES,NREF),WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          XT(MAXRES),P(MAXRES),XIT(MAXRES),DF(MAXRES),SCR(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,PT5,TWO
      DATA ZERO,PT5,TWO/0.0D0,0.5D0,2.0D0/
C
      A = MIN(AX,CX)
      B = MAX(AX,CX)
      V = BX
      W = V
      X = V
      D = ZERO
      E = ZERO
      FX = FB
      FV = FX
      FW = FX
      CALL DF1DIM(DX,X,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,MAXRES,
     +            IREF,NREF,NVAR)
      DV = DX
      DW = DX
C
C  Main loop over iterations:
C
      DO 30 ITER = 1,ITMAX
        XM = PT5 * (A + B)
        TOL1 = TOL * ABS(X) + ZEPS
        TOL2 = TWO * TOL1
        IF(ABS(X-XM).LE.(TOL2-PT5*(B-A))) GOTO 40
C
C  Initialize the D's to an out of bracket value:
C
        IF(ABS(E).GT.TOL1) THEN
          D1 = TWO * (B - A)
          D2 = D1
C
C  Secant method:
C
          IF(DW.NE.DX) D1 = (W - X) * DX / (DX - DW)
C
C  Secant method with the other stored point:
C
          IF(DV.NE.DX) D2 = (V - X) * DX / (DX - DV)
C
C  Which of these two estimates of D
C
          U1 = X + D1
          U2 = X + D2
          OK1 = ((A-U1)*(U1-B).GT.ZERO).AND.(DX*D1.LE.ZERO)
          OK2 = ((A-U2)*(U2-B).GT.ZERO).AND.(DX*D2.LE.ZERO)
C
C  Movement on the step before last:
C
          OLDE = E
          E = D
C
C  Take only an acceptable D.  If both are acceptable, take the
C  smaller:
C
          IF(.NOT.(OK1.OR.OK2)) THEN
            GO TO 10
          ELSE IF (OK1.AND.OK2) THEN
            IF(ABS(D1).LT.ABS(D2)) THEN
              D = D1
            ELSE
              D = D2
            END IF
          ELSE IF (OK1) THEN
            D = D1
          ELSE
            D = D2
          END IF
          IF(ABS(D).GT.ABS(PT5*OLDE)) GO TO 10
          U = X + D
          IF(U-A.LT.TOL2.OR.B-U.LT.TOL2) D = SIGN(TOL1,XM-X)
          GO TO 20
        END IF
C
C  Decide which segment by the sign of the derivative:
C
   10   IF(DX.GE.0.0D0) THEN
          E = A - X
        ELSE
          E = B - X
        END IF
C
C  Bisect, not by golden section:
C
        D = PT5 * E
   20   IF(ABS(D).GE.TOL1) THEN
          U = X + D
          CALL NBFDIM(FU,U,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
        ELSE
          U = X + SIGN(TOL1,D)
          CALL NBFDIM(FU,U,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,IREF,
     +                NREF,NVAR)
C
C  If the minimum step in the downhill direction takes us uphill, then
C  we are done:
C
          IF(FU.GT.FX) GO TO 40
        END IF
        CALL DF1DIM(DU,U,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,MAXRES,
     +              IREF,NREF,NVAR)
C
C  Now, do all the housecleaning (sigh):
C
        IF(FU.LE.FX) THEN
          IF(U.GE.X) THEN
            A = X
          ELSE
            B = X
          END IF
          V  = W
          FV = FW
          DV = DW
          W  = X
          FW = FX
          DW = DX
          X  = U
          FX = FU
          DX = DU
        ELSE
          IF(U.LT.X) THEN
            A = U
          ELSE
            B = U
          END IF
          IF(FU.LE.FW.OR.W.EQ.X) THEN
            V  = W
            FV = FW
            DV = DW
            W  = U
            FW = FU
            DW = DU
          ELSE IF(FU.LE.FV.OR.V.EQ.X.OR.V.EQ.W) THEN
            V  = U
            FV = FU
            DV = DU
          END IF
        END IF
   30 CONTINUE
      WRITE(LFNPR,900)
C
   40 XMIN = X
      FRET = FX
      RETURN
C
  900 FORMAT(1X,'DBRENT exceeded maximum iterations.')
      END
C***********************************************************************
      SUBROUTINE NBFDIM(FX,X,WGT,Q,XT,P,XIT,IPTR,NRES,MAXRES,
     +                  IREF,NREF,NVAR)
C***********************************************************************
C 23-Mar-06  EDG  Name change from F1DIM to avoid GAMESS conflict
C 22-Jun-95  EDG  Converted from FUNCTION to SUBROUTINE
C-----------------------------------------------------------------------
C
C  Modified version of NBFDIM, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), p. 301.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          XT(MAXRES),P(MAXRES),XIT(MAXRES),IPTR(MAXRES),
     +          NRES(NREF)
C
      DO 10 J = 1,NVAR
        XT(J) = P(J) + X * XIT(J)
   10 CONTINUE
      CALL NBPHI(FX,XT,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
      RETURN
      END
C***********************************************************************
      SUBROUTINE DF1DIM(DX,X,OCC,WGT,Q,XT,P,XIT,DF,SCR,IPTR,NRES,
     +                  MAXRES,IREF,NREF,NVAR)
C***********************************************************************
C 22-Jun-95  EDG  Converted from FUNCTION to SUBROUTINE
C-----------------------------------------------------------------------
C
C  Modified version of DF1DIM, Press, Flannery, Teukolsky, and Vetterling,
C  Numerical Recipes (FORTRAN edition, 1988), p. 306.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION OCC(NDIM,MAXRES,NREF),WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          XT(MAXRES),P(MAXRES),XIT(MAXRES),DF(MAXRES),SCR(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO
      DATA ZERO/0.0D0/
C
      DO 10 J = 1,NVAR
        XT(J) = P(J) + X * XIT(J)
   10 CONTINUE
C
      CALL NBDPHI(DF,XT,OCC,WGT,Q,SCR,IPTR,NRES,MAXRES,IREF,NREF,
     +            NVAR)
C
      DX = ZERO
      DO 20 J = 1,NVAR
        DX = DX + DF(J) * XIT(J)
   20 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBPHI(FT,T,WGT,Q,IPTR,NRES,MAXRES,IREF,NREF,
     +                 NVAR)
C***********************************************************************
C 22-Jun-95  EDG  Converted FUNCTION to SUBROUTINE
C 03-Jan-01  FAW  Added MINJOB=1 option for HYBDIR call
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          T(MAXRES),IPTR(MAXRES),NRES(NREF)
      COMMON/BENDS/CX,CY,CZ,C1,C2,C3,C4,C5,MINJOB
C
      SAVE ZERO,EPS,ONE,TWO,THIRTY,PI
      DATA ZERO,EPS,ONE,TWO,THIRTY/0.0D0,1.0D-6,1.0D0,2.0D0,3.0D1/
      DATA PI/3.14159265357979D0/FOUR/4.0D0/
C
C  Bend job?
C
      IF(MINJOB.EQ.1) THEN
        THETA=T(1)
        PHI=T(2)
        IF(THETA.LT.ZERO) THEN
          FT = 1.D4*(ONE-THETA)
        ELSE IF(THETA.GT.PI) THEN
          FT = 1.D4*(ONE+THETA-PI)
          RETURN
        ENDIF
        IF(PHI.LT.ZERO) THEN
          FT = 1.D4*(ONE-PHI)
        ELSE IF(PHI.GT.TWO*PI) THEN
          FT = 1.D4*(ONE+PHI-TWO*PI)
          RETURN
        ENDIF
        ST=SIN(THETA)
        CT=COS(THETA)
        SP=SIN(PHI)
        CP=COS(PHI)
        X=ST*CP
        Y=ST*SP
        Z=CT
        HYB = CX*X + CY*Y + CZ*Z + C1*(Z*Z) + C2*(X*X - Y*Y)
     +        + C3*X*Y + C4*X*Z + C5*Y*Z
c        FT = -HYB*HYB
        FT = ONE/(ONE + ABS(HYB))
        IF(ABS(HYB).LT.EPS) THEN
          FT = FT*(ONE + ABS(THETA - PI/FOUR) + ABS(PHI - PI/FOUR))
        ENDIF
        IF(ABS(THETA).LT.0.1D0.OR.ABS(PI-THETA).LT.0.1D0) THEN
          FT = FT*(ONE + ABS(PHI))
        ENDIF
        RETURN
      ENDIF
C
C  If any of the T's are negative, reject the choice:
C
      DO 5 IV = 1,NVAR
        IF(T(IV).LT.ZERO) THEN
          FT = THIRTY
          RETURN
        END IF
    5 CONTINUE
C
C  Compute weights for each resonance structure:
C
      QQ = ZERO
      NLEAD = 0
      DO 10 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).GT.0) THEN
          WGT(IRES,IREF) = EXP(-T(IPTR(IRES)))
          QQ = QQ + WGT(IRES,IREF)
        ELSE
          WGT(IRES,IREF) = ZERO
          IF(IPTR(IRES).EQ.0) NLEAD = NLEAD + 1
        END IF
   10 CONTINUE
C
C  Compute the weights for the leading resonance structures, assuming
C  normalization of the weights:
C
      QQ = (ONE - QQ) / NLEAD
      DO 20 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) WGT(IRES,IREF) = QQ
   20 CONTINUE
C
C  Compute the deviation, ERROR:
C
      CALL GETERR(ERROR,WGT,Q,NRES,MAXRES,IREF,NREF)
C
C  Add constraints if WGTs are out-of-bounds.  Note that the constraint
C  is squared to assure continuous derivatives:
C
      DO 30 IRES = 1,NRES(IREF)
        IF(WGT(IRES,IREF).GT.ONE) THEN
          ERROR = ERROR + (WGT(IRES,IREF)-ONE) * (WGT(IRES,IREF)-ONE)
        ELSE IF(WGT(IRES,IREF).LT.ZERO) THEN
          ERROR = ERROR + THIRTY * WGT(IRES,IREF) * WGT(IRES,IREF)
        END IF
   30 CONTINUE
C
C  Finally, if one of the T's becomes large, add a small amount to
C  ERROR.  This prevents floating point overflows in MNBRAK:
C
      DO 40 IV = 1,NVAR
        IF(T(IV).GT.THIRTY) THEN
          ERROR = ERROR + EPS * (T(IV) - THIRTY) * (T(IV) - THIRTY)
        END IF
   40 CONTINUE
      FT = ERROR
      END
C***********************************************************************
      SUBROUTINE NBDPHI(DER,T,OCC,WGT,Q,SCR,IPTR,NRES,MAXRES,IREF,
     +                  NREF,NVAR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION DER(MAXRES),T(MAXRES),OCC(NDIM,MAXRES,NREF),
     +          WGT(MAXRES,NREF),Q(MAXRES,MAXRES),SCR(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,EPS,ONE,TWO,THIRTY
      DATA ZERO,EPS,ONE,TWO,THIRTY/0.0D0,1.0D-6,1.0D0,2.0D0,3.0D1/
C
C  Compute weights for each resonance structure:
C
      QQ = ZERO
      NLEAD = 0
      DO 10 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).GT.0) THEN
          WGT(IRES,IREF) = EXP(-T(IPTR(IRES)))
          QQ = QQ + WGT(IRES,IREF)
        ELSE
          WGT(IRES,IREF) = ZERO
          IF(IPTR(IRES).EQ.0) NLEAD = NLEAD + 1
        END IF
   10 CONTINUE
C
C  Compute the weights for the leading resonance structures, assuming
C  normalization of the weights:
C
      QQ = (ONE - QQ) / NLEAD
      DO 20 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) WGT(IRES,IREF) = QQ
   20 CONTINUE
C
C  Compute the derivative at the `point' WGT:
C
      CALL GETGRD(DER,T,WGT,Q,OCC,SCR,IPTR,NRES,MAXRES,IREF,NREF,
     +            NVAR)
C
C  If a weight is larger than one, adjust its corresponding derivative:
C  (What's done here is not quite correct, but should do the trick.)
C
      DO 30 IRES = 1,NRES(IREF)
        IF(WGT(IRES,IREF).GT.ONE.AND.IPTR(IRES).NE.0) THEN
          DER(IPTR(IRES)) = DER(IPTR(IRES)) - TWO * (WGT(IRES,IREF)-ONE)
        END IF
   30 CONTINUE
C
C  Finally, if any T becomes too large, adjust its derivative to correct
C  the situation:
C
      DO 40 IV = 1,NVAR
        IF(T(IV).GT.THIRTY) THEN
          DER(IV) = DER(IV) + TWO * EPS * (T(IV) - THIRTY)
        END IF
   40 CONTINUE
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  SIMULATED ANNEALING ROUTINES:
C
C      SUBROUTINE ANNEAL(XLAM,YLAM,ERRREF,ERR,RHOSTR,WGT,T,SCR,Q,
C     +                  IPTR,NRES,MAXRES,IREF,NREF,ISEED,NVAR,ITT)
C      SUBROUTINE SETANN(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C      SUBROUTINE ADJUST(T,IT,NVAR,SSIZE,STEP,ISEED)
C      SUBROUTINE METROP(NSUC,FW2,FW1,FT,ISEED)
C      SUBROUTINE RESTOR(T,IT,NVAR,SSIZE)
C
C***********************************************************************
      SUBROUTINE ANNEAL(XLAM,YLAM,ERRREF,ERR,RHOSTR,WGT,T,SCR,Q,
     +                  IPTR,NRES,MAXRES,IREF,NREF,ISEED,NVAR,ITT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION ERRREF(NREF),ERR(NREF),RHOSTR(MAXRES,NREF),
     +          WGT(MAXRES,NREF),T(MAXRES),
     +          SCR(MAXRES),Q(MAXRES,MAXRES),IPTR(MAXRES),
     +          NRES(NREF)
C
      SAVE ZERO,THRFW,THRS,PT5,THREE,TEN
      SAVE MAXIT
      DATA ZERO,THRFW,THRS,PT5,THREE,TEN/0.0D0,1.0D-5,2.0D-5,0.5D0,
     +     3.0D0,1.0D1/
      DATA MAXIT/1000/
C
C  Get a list of optimizable parameters:
C
      CALL SETANN(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  Set a few initial parameters.  These will be changed throughout the
C  optimization:
C
      STHR = 0.03D0
      STEP = STHR
      STEPMX = ZERO
      GLOBAL = ZERO
C
C  Compute the fractional reduction for this set of variational parameters:
C
      CALL CMPFW(FW1,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,NRES,
     +           MAXRES,IREF,NREF)
C
C  Take a few iterations to determine an appropriate ``temperature''
C  FT to start with:
C
      FT = ZERO
      DO 10 I = 1,10*NVAR
C
C  Make an arbitrary adjustment to one of the resonance weights and
C  compute the new fractional reduction:
C
        CALL ADJUST(T,IT,NVAR,SSIZE,STEP,ISEED)
        CALL CMPFW(FW2,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,NRES,
     +             MAXRES,IREF,NREF)
        STEPMX = MAX(STEPMX,ABS(SSIZE))
        CALL RESTOR(T,IT,NVAR,SSIZE)
C
C  Keep track of the largest change in FW, which will determine FT:
C
        FT = MAX(FT,ABS(FW2-FW1))
        FT = MIN(FT,PT5)
   10 CONTINUE
      IF(JPRINT(57).NE.0) WRITE(LFNPR,900) FW1,FT,STEPMX
C
C  Take up to MAXIT grand iterations:
C
      ITOT = 0
      NTST = 10 * NVAR
      DO 50 IG = 1,MAXIT
C
C  Loop over micro iterations.  On each iteration, adjust one of the
C  resonance weights, compute the new fractional reduction, and decide
C  whether the adjustment is acceptable (NSUC.GT.0).  If not, restore
C  the previous weights.  Also, keep track of the largest acceptable
C  adjustment made and the largest and smallest accepted FW:
C
        ITT = IG
        ITER = 0
        ICNT = 0
        NSUC = 0
        FWMAX = FW1
        FWMIN = FW1
        STEPMX = ZERO
        DO 20 IM = 1,100*NVAR
          ITER = ITER + 1
          CALL ADJUST(T,IT,NVAR,SSIZE,STEP,ISEED)
          CALL CMPFW(FW2,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,NRES,
     +               MAXRES,IREF,NREF)
          CALL METROP(NSUC,FW2,FW1,FT,ISEED)
          IF(NSUC.EQ.0) THEN
            CALL RESTOR(T,IT,NVAR,SSIZE)
          ELSE
            ICNT = ICNT + 1
            FW1 = FW2
            FWMAX = MAX(FWMAX,FW1)
            FWMIN = MIN(FWMIN,FW1)
            STEPMX = MAX(STEPMX,ABS(SSIZE))
            IF(GLOBAL.LT.FWMAX) THEN
              ITG = IG
              GLOBAL = FWMAX
              CALL COPY(T,SCR,MAXRES,NVAR,1)
            END IF
            IF(NSUC.EQ.NTST) GOTO 30
          END IF
   20   CONTINUE
C
C  Anneal Scheduling - Adjust ``temperature'' and maximum step size:
C
   30   CONTINUE
        FT = 0.9D0 * FT
        STEP = 1.5D0 * STEPMX
C
C  Print info about optimization:
C
        ITOT = ITOT + ITER
        IF(JPRINT(57).NE.0) WRITE(LFNPR,910) IG,ITER,ICNT,ITOT,FW1,
     +                                       FWMAX,FWMIN,FT,STEPMX
C
C  Has the optimization converged?
C
        IF(FWMAX-FWMIN.LT.THRFW) THEN
C
C  If so, check to see that we have the best FW found so far.  If not,
C  restart the optimization:
C
          IF(ABS(FWMAX-GLOBAL).GT.THRFW) THEN
            CALL COPY(SCR,T,MAXRES,NVAR,1)
            STHR = STHR / TEN
            STEP = MIN(STHR,THREE*STEP)
            FT = ZERO
            CALL CMPFW(FW1,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,
     +                 NRES,MAXRES,IREF,NREF)
            IF(JPRINT(57).NE.0) WRITE(LFNPR,920) ITG,FW1
            DO 40 I = 1,10*NVAR
              CALL ADJUST(T,IT,NVAR,SSIZE,STEP,ISEED)
              CALL CMPFW(FW2,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,
     +                   NRES,MAXRES,IREF,NREF)
              CALL RESTOR(T,IT,NVAR,SSIZE)
              FT = MAX(FT,ABS(FW2-FW1))
              FT = MIN(FT,PT5)
   40       CONTINUE
C
C  Otherwise, abort the iterations.  But first make sure that STEPMX has
C  converged:
C
          ELSE IF(STEPMX.LT.THRS) THEN
            GOTO 60
          END IF
        END IF
   50 CONTINUE
      ITT = -MAXIT
C
C  We're finished (at last)!!  Recompute the fractional reduction
C  and return to the calling routine:
C
   60 CONTINUE
      CALL CMPFW(FW1,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,NRES,
     +           MAXRES,IREF,NREF)
      RETURN
C
  900 FORMAT(/1X,' Grand  Micro          Total    ',
     + '         max      min       ft      max',/1X,'  iter   iter',
     + '   Acc     iter    f(w)     f(w)     f(w)    (temp)    step',
     + /1X,'-------------------------------------------------------',
     + '-------------------',/1X,'    0',23X,F9.5,18X,2F9.5)
  910 FORMAT(1X,I5,I7,I7,I8,1X,5F9.5)
  920 FORMAT(1X,'        Restarting optimization at iteration',
     + I4,', [f(w) =',F8.5,']')
      END
C***********************************************************************
      SUBROUTINE SETANN(RHOSTR,WGT,T,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION RHOSTR(MAXRES,NREF),WGT(MAXRES,NREF),T(MAXRES),
     +          IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,EPS,ONE
      DATA ZERO,EPS,ONE/0.0D0,1.0D-5,1.0D0/
C
C  Setup pointer array:
C
      DO 10 IRES = 1,NRES(IREF)
        IPTR(IRES) = -1
   10 CONTINUE
C
C  Do not optimize the weight of the reference resonance structure;
C  rather determine its weight by normalization:
C
      IPTR(1) = 0
      WGT(1,IREF) = -ONE
C
C  Optimize the weights of the rest of the structures, constraining
C  resonance structures with equivalent rho* to have equal weights:
C
      NVAR = 0
   30 JRES = 0
      WREF = ZERO
      DO 40 IRES = 1,NRES(IREF)
        IF(RHOSTR(IRES,IREF).GE.ZERO.AND.WGT(IRES,IREF).GE.WREF) THEN
          JRES = IRES
          WREF = WGT(IRES,IREF)
          RHO  = RHOSTR(IRES,IREF)
        END IF
   40 CONTINUE
C
C  Compute the optimizable parameter T and set pointers appropriately:
C
      IF(JRES.NE.0) THEN
        NVAR = NVAR + 1
        T(NVAR) = WGT(JRES,IREF)
        DO 50 IRES = 1,NRES(IREF)
          IF(ABS(RHOSTR(IRES,IREF)-RHO).LT.EPS) THEN
            IF(WGT(IRES,IREF).GE.ZERO) THEN
              IPTR(IRES) = NVAR
              WGT(IRES,IREF) = -ONE
            END IF
          END IF
   50   CONTINUE
        GOTO 30
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE ADJUST(T,IT,NVAR,SSIZE,STEP,ISEED)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION T(NVAR)
C
      SAVE ONE
      DATA ONE/1.0D0/
C
C  Determine which T to adjust and the magnitude of the adjustment:
C
      IT = 1
      IF(NVAR.GT.1) IT = INTRN(1,NVAR,ISEED)
      SMAX = MIN(STEP,ONE-T(IT))
      SMIN = -MIN(STEP,T(IT))
      SSIZE = REALRN(SMIN,SMAX,ISEED)
      T(IT) = T(IT) + SSIZE
      RETURN
      END
C***********************************************************************
      SUBROUTINE METROP(NSUC,FW2,FW1,FT,ISEED)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      SAVE ZERO,ONE,BIG
      DATA ZERO,ONE,BIG/0.0D0,1.0D0,-1.0D5/
C
C  Use the Metropolis algorithm to determine whether FW2 is acceptable:
C
C  If FW2 is very large (and negative), then it's unacceptable:
C
      IF(FW2.LT.BIG) THEN
        NSUC = 0
        RETURN
      END IF
C
C  If FW2 is greater than FW1, it's acceptable:
C
      NSUC = NSUC + 1
      IF(FW2.GT.FW1) RETURN
C
C  Otherwise, use the Boltzmann probability to determine whether FW2 is
C  acceptable:
C
      TEST = EXP((FW2-FW1)/FT)
      IF(TEST.GT.REALRN(ZERO,ONE,ISEED)) RETURN
      NSUC = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE RESTOR(T,IT,NVAR,SSIZE)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION T(NVAR)
C
      T(IT) = T(IT) - SSIZE
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NRT BRANCHING ROUTINES:
C
C      SUBROUTINE IDCTRL(IREF)
C      SUBROUTINE IDINIT(ICNT,NC,NH,NUM)
C      SUBROUTINE IDINCR(ICNT,NC,NH,NUM)
C      FUNCTION IDTEST(JDXT,IC,NUM)
C      SUBROUTINE IDCLR
C      SUBROUTINE IDCOPY
C      SUBROUTINE IDREST
C
C***********************************************************************
      SUBROUTINE IDCTRL(IREF)
C***********************************************************************
C
C  Control loops over identical TOPO matrices:
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
C
      JREF = IREF
      IREF = 0
      IF(JCNT.EQ.0) RETURN
C
C  Trim IDX and JDX lists if necessary:
C
      DO 10 I = JCNT,1,-1
        IF(JDX(I).EQ.IDX(I)) THEN
          IDX(I) = 0
          JDX(I) = 0
          JCNT = JCNT - 1
        ELSE
          GOTO 20
        END IF
   10 CONTINUE
   20 CONTINUE
      IF(JCNT.EQ.0) RETURN
      IREF = JREF
      RETURN
      END
C***********************************************************************
      SUBROUTINE IDINIT(ICNT,NC,NH,NUM)
C***********************************************************************
C
C  Branch at this point -- initialize the indices IDX(ICNT) and JDX(ICNT).
C  The NBCHSE routines have found NC acceptable orbitals of which we only
C  want to select NUM.  Use the NH high occupancy orbitals whenever possible.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
C
      IF(ICNT.GT.MAXATM) THEN
        CALL NBHALT('Redimension arrays IDX and JDX in COMMON/NBALT/.')
      ELSE
        JCNT = ICNT
      END IF
C
C  Initialize the JDX index:
C
      MULT = 1
      ITMP = 1
      JDX(JCNT) = 0
      DO 10 I = 1,NUM
        JDX(JCNT) = JDX(JCNT) + MULT * ITMP
        MULT = MULT * 10
        ITMP = ITMP + 1
   10 CONTINUE
C
C  Temporarily store JDX while we search for the limiting index, IDX:
C
      LDX = JDX(ICNT)
C
C  Determine the limiting index:
C
   20 IDX(ICNT) = JDX(ICNT)
      CALL IDINCR(ICNT,NC,NH,NUM)
      IF(JDX(ICNT).GT.0) GOTO 20
C
C  Restore the JDX index:
C
      JDX(ICNT) = LDX
      RETURN
      END
C***********************************************************************
      SUBROUTINE IDINCR(ICNT,NC,NH,NUM)
C***********************************************************************
C
C  Increment the branching index JDX(ICNT), given NC orbitals of which we
C  only want to choose NUM (make sure that the NH high occupancy orbitals
C  appear in the index):
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
C
      DIMENSION ITMP(10)
C
C  Halt if we are illegally trying to increment an index (can only
C  increment the last index in the list):
C
      IF(ICNT.NE.JCNT)
     +  CALL NBHALT('IDINCR: Error in the NBCHSE branching algorithm.')
C
C  Setup ITMP array [containing the digits in JDX(ICNT)]:
C
   10 CONTINUE
      LDX = JDX(ICNT)
      DO 20 I = 1,NUM
        ITMP(I) = MOD(LDX,10)
        LDX = LDX / 10
   20 CONTINUE
C
      DO 50 I = 1,NUM
        IF(I.LT.NUM.AND.ITMP(I)+1.LT.ITMP(I+1)) THEN
          ITMP(I) = ITMP(I) + 1
          DO 30 J = 1,I-1
            ITMP(J) = J
   30     CONTINUE
          GOTO 60
        ELSE IF(I.EQ.NUM.AND.ITMP(I).LT.NC) THEN
          ITMP(I) = ITMP(I) + 1
          DO 40 J = 1,I-1
            ITMP(J) = J
   40     CONTINUE
          GOTO 60
        END IF
   50 CONTINUE
C
C  Trying to increment JDX(ICNT) beyond its limiting value.  So, set
C  the index to zero and return to the calling routine:
C
      JDX(ICNT) = 0
      RETURN
C
C  Put JDX(ICNT) back together:
C
   60 CONTINUE
      MULT = 1
      JDX(ICNT) = 0
      DO 70 I = 1,NUM
        JDX(ICNT) = JDX(ICNT) + MULT * ITMP(I)
        MULT = MULT * 10
   70 CONTINUE
C
C  Test JDX(ICNT), making sure that all high occupancy orbitals appear
C  in it:
C
      IF(NH.NE.0.AND.NH.NE.NC) THEN
        JDXT = JDX(ICNT)
        DO 80 I = 1,NH
          IF(IDTEST(JDXT,I,NUM).EQ.0) GOTO 10
   80   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      FUNCTION IDTEST(JDXT,IC,NUM)
C***********************************************************************
C
C  Test index JDXT to determine whether the ICth orbital should be
C  accepted (IDTEST=1) or rejected (IDTEST=0).  To be accepted, at
C  least one of the following must be true:
C
C  a)  JDXT is negative and IC is less than or equal to NUM;
C
C  b)  the digit IC appears in the integer JDXT.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      IDTEST = 1
      IF(JDXT.LT.0) THEN
        IF(IC.GT.NUM) IDTEST = 0
        RETURN
      END IF
C
C  Search for IC in JDXT:
C
      LDX = JDXT
      DO 10 I = 1,NUM
        IF(MOD(LDX,10).EQ.IC) RETURN
        LDX = LDX / 10
   10 CONTINUE
C
C  IC was not found in JDXT:
C
      IDTEST = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE IDCLR
C***********************************************************************
C
C  Clear the indices used to guide the selection of candidate NBO's
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
C
C  Note that IDX and JDX are dimensioned MAXATM since the number of branches
C  encountered in the NBCHSE routines should scale as the number of atoms.
C  This dimension may need modification (warning messages will be printed
C  if so):
C
      DO 10 I = 1,MAXATM
        IDX(I) = 0
        JDX(I) = 0
   10 CONTINUE
      DO 20 I = 1,MAXBAS
        KDX(I) = 0
   20 CONTINUE
      JCNT = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE IDCOPY
C***********************************************************************
C
C  Make a backup copy of the branching indices:
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
      COMMON/NBALT1/IIDX(MAXATM),JJDX(MAXATM),KKDX(MAXBAS),JJCNT
C
      DO 10 I = 1,MAXATM
        IIDX(I) = IDX(I)
        JJDX(I) = JDX(I)
   10 CONTINUE
      DO 20 I = 1,MAXBAS
        KKDX(I) = KDX(I)
   20 CONTINUE
      JJCNT = JCNT
      RETURN
      END
C***********************************************************************
      SUBROUTINE IDREST
C***********************************************************************
C
C  Restore the contents of COMMON/NBALT/ from COMMON/NBALT1/:
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBALT/IDX(MAXATM),JDX(MAXATM),KDX(MAXBAS),JCNT
      COMMON/NBALT1/IIDX(MAXATM),JJDX(MAXATM),KKDX(MAXBAS),JJCNT
C
      DO 10 I = 1,MAXATM
        IDX(I) = IIDX(I)
        JDX(I) = JJDX(I)
   10 CONTINUE
      DO 20 I = 1,MAXBAS
        KDX(I) = KKDX(I)
   20 CONTINUE
      JCNT = JJCNT
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NRT UTILITY ROUTINES:
C
C      SUBROUTINE ATMORD(WIBERG,NRT)
C      CHARACTER*2 FUNCTION CHARAT(IZ)
C      SUBROUTINE CMPFW(FW,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,NRES,
C     +               MAXRES,IREF,NREF)
C      SUBROUTINE EXPWGT(RHOW,RHOSTR,WGT,IRANK,NRES,MAXRES,IREF,NREF,N)
C      SUBROUTINE FETREF(T,LABEL,IBXM,IREF)
C      SUBROUTINE GETDM(T,DM,TR,SCR1,SCR2,SCR3,IST,IOCC,NVAL,IERR)
C      SUBROUTINE GETDW(DW,GAMMAW,WGT,LVAL,NREF,NVAL)
C      SUBROUTINE GETERR(ERROR,WGT,Q,NRES,MAXRES,IREF,NREF)
C      SUBROUTINE GETFDM(DM,DELOC,TNAO,TBO,SCR)
C      SUBROUTINE GETGRD(DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,IREF,NREF,
C     +                  NVAR)
C      SUBROUTINE GETHES(HESS,DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,
C     +                  IREF,NREF,NVAR,IFLG)
C      SUBROUTINE GETREL(IBO,IBOP,IUNIT1,IUNIT2,IREL,IA,IB,IC,ID)
C      FUNCTION INTRN(ILLIM,IULIM,ISEED)
C      SUBROUTINE MULTI(FW,GAMMAW,WGT,SCR,IRESET,LVAL,NREF,NVAL)
C      SUBROUTINE NRTINP(IESS,INRT)
C      SUBROUTINE PRTWGT(IRES,IR,NUM,WGT,II,IAT,JJ,JAT,IFLG)
C      FUNCTION RANNB(IDUM)
C      FUNCTION REALRN(XLLIM,XULIM,ISEED)
C      SUBROUTINE SVTREF(T,LABEL,IBXM,IREF)
C      SUBROUTINE TOPCMP(IDUP,IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C      SUBROUTINE TOPGET(IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C      SUBROUTINE TOPOUT
C      SUBROUTINE TOPSTR(IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C      SUBROUTINE TOPZER(NATOMS)
C      SUBROUTINE WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
C      SUBROUTINE WGTPR(WGTM,MAP,NWGT)
C      SUBROUTINE WNORM(WGT,IRESET,IW,IP,SKIP)
C
C***********************************************************************
      SUBROUTINE ATMORD(WIBERG,NRT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBTOPO/IORDER(MAXATM),JORDER(MAXATM),NTOPO(MAXATM,MAXATM),
     +            N3CTR,I3CTR(10,3)
C
      DIMENSION WIBERG(NATOMS,NATOMS),SCR(MAXATM),IPT(MAXATM)
C
      SAVE ZERO
      DATA ZERO/0.0D0/
C
C  Determine the order in which atoms will be searched for electron pairs.
C  Special case of a single atoms:
C
      IF(NATOMS.EQ.1) THEN
        IORDER(1) = 1
        GOTO 100
      END IF
C
C  Copy ITOPO to NTOPO if NRT:
C
      IF(NRT.EQ.1) THEN
        DO 5 I = 1,NATOMS
          DO 4 J = 1,NATOMS
            NTOPO(J,I) = ITOPO(J,I)
    4     CONTINUE
    5   CONTINUE
      END IF
C
C  Otherwise, start where we'll have the toughest time finding bond orbitals
C  (i.e. where the formal bond order of the TOPO matrix is greater than the
C  corresponding Wiberg bond index):
C
      KAT = 1
      LAT = 2
      GMAX = ZERO
      DO 20 JAT = 2,NATOMS
        DO 10 IAT = 1,JAT-1
          GTMP = DFLOAT(NTOPO(IAT,JAT)) - WIBERG(IAT,JAT)
          IF(GTMP.GT.GMAX) THEN
            GMAX = GTMP
            KAT  = IAT
            LAT  = JAT
          END IF
   10   CONTINUE
   20 CONTINUE
      IORDER(1) = KAT
      IORDER(2) = LAT
C
C  Now add atoms to IORDER according to atom connectivity:
C
C  Pointers:   IPTR  --  current location in the IORDER array
C              JPTR  --  current location in the IPT array
C              INXT  --  next atom to search in the IORDER list
C              ICNT  --  number of atoms already in the IORDER list
C              JCNT  --  next atom to search in the IORDER list if
C                        we've run into a deadend in the connectivity
C
      ICNT = 2
      INXT = 1
      JCNT = 1
C
C  IORDER(IPTR) is the `current' atom that we're working on:
C
   30 IPTR = INXT
C
C  Rank the differences betweeen the NTOPO and Wiberg matrix elements:
C
        DO 40 IAT = 1,NATOMS
          SCR(IAT) = NTOPO(IAT,IORDER(IPTR)) - WIBERG(IAT,IORDER(IPTR))
   40   CONTINUE
        CALL RANK(SCR,NATOMS,NATOMS,IPT)
C
C  Loop over atoms:
C
        I1ST = 1
        DO 60 JPTR = 1,NATOMS
C
C  Make sure that the atom IPT(JPTR) is bonded to the current atom:
C
          IF(NTOPO(IORDER(IPTR),IPT(JPTR)).NE.0) THEN
C
C  And also make sure that we're not looking at a diagonal element of
C  the TOPO matrix:
C
            IF(IORDER(IPTR).NE.IPT(JPTR)) THEN
C
C  Check to see whether atom IPT(JPTR) is already in the list IORDER:
C
              IFLG = 1
              DO 50 I = 1,ICNT
                IF(IORDER(I).EQ.IPT(JPTR)) IFLG = 0
   50         CONTINUE
C
C  If it's not there, add it to the list:
C
              IF(IFLG.EQ.1) THEN
                ICNT = ICNT + 1
                IORDER(ICNT) = IPT(JPTR)
C
C  Are we finished?
C
                IF(ICNT.EQ.NATOMS) GOTO 100
                IF(I1ST.EQ.1) THEN
                  I1ST = 0
                  INXT = ICNT
                END IF
              END IF
            END IF
          END IF
   60   CONTINUE
C
C  Before looking at the next atom, make sure that we have a new atom
C  to search on:
C
        IF(I1ST.EQ.1) THEN
C
C  If we don't have one, take the next atom in the IORDER list:
C
          JCNT = JCNT + 1
          INXT = JCNT
C
C  Check to see whether we're finished:
C
          IF(INXT.GT.NATOMS) GOTO 100
C
C  If there are multiple molecular units, INXT can surpass ICNT.  If so,
C  arbitrarily choose the next atom to be examined:
C
          IF(INXT.GT.ICNT) THEN
            KPTR = 0
   70       KPTR = KPTR + 1
            JFLG = 1
            DO 80 I = 1,ICNT
              IF(IORDER(I).EQ.KPTR) JFLG = 0
   80       CONTINUE
            IF(JFLG.EQ.0) GOTO 70
            ICNT = ICNT + 1
            IORDER(ICNT) = KPTR
C
C  Again, check to see whether we're finished:
C
            IF(ICNT.EQ.NATOMS) GOTO 100
          END IF
        END IF
      GOTO 30
C
  100 CONTINUE
      RETURN
      END
C***********************************************************************
      CHARACTER*2 FUNCTION CHARAT(IZ)
C***********************************************************************
C 25-May-2005  EDG  Atomic symbols to Z=118
C-----------------------------------------------------------------------
C
C  RETURN ATOMIC SYMBOL FOR NUCLEAR CHARGE IZ (.LE. 118):
C
      CHARACTER*2 NAME(118),BLANK,GHOST
C
      SAVE NAME,BLANK,GHOST
      DATA NAME/' H','He','Li','Be',' B',' C',' N',' O',' F','Ne',
     + 'Na','Mg','Al','Si',' P',' S','Cl','Ar',' K','Ca','Sc','Ti',
     + ' V','Cr','Mn','Fe','Co','Ni','Cu','Zn','Ga','Ge','As',
     + 'Se','Br','Kr','Rb','Sr',' Y','Zr','Nb','Mo','Tc','Ru',
     + 'Rh','Pd','Ag','Cd','In','Sn','Sb','Te',' I','Xe','Cs',
     + 'Ba','La','Ce','Pr','Nd','Pm','Sm','Eu','Gd','Tb','Dy',
     + 'Ho','Er','Tm','Yb','Lu','Hf','Ta',' W','Re','Os','Ir',
     + 'Pt','Au','Hg','Tl','Pb','Bi','Po','At','Rn','Fr','Ra',
     + 'Ac','Th','Pa',' U','Np','Pu','Am','Cm','Bk','Cf','Es',
     + 'Fm','Md','No','Lr','Rf','Db','Sg','Bh','Hs','Mt','10',
     + '11','12','13','14','15','16','17','18'/
      DATA BLANK,GHOST/'  ','gh'/
C
      IF(IZ.LT.0.OR.IZ.GT.118) CHARAT = BLANK
      IF(IZ.GT.0) CHARAT = NAME(IZ)
      IF(IZ.EQ.0) CHARAT = GHOST
      RETURN
      END
C***********************************************************************
      SUBROUTINE CMPFW(FW,XLAM,YLAM,ERRREF,ERR,T,WGT,Q,IPTR,
     +                 NRES,MAXRES,IREF,NREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION ERRREF(NREF),ERR(NREF),T(MAXRES),WGT(MAXRES,NREF),
     +        Q(MAXRES,MAXRES),IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,ONE,BIG
      DATA ZERO,ONE,BIG/0.0D0,1.0D0,-2.0D5/
C
C  Compute weights for each resonance structure:
C
      QQ = ZERO
      NLEAD = 0
      DO 10 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).GT.0) THEN
          WGT(IRES,IREF) = T(IPTR(IRES))
          QQ = QQ + WGT(IRES,IREF)
        ELSE
          WGT(IRES,IREF) = ZERO
          IF(IPTR(IRES).EQ.0) NLEAD = NLEAD + 1
        END IF
   10 CONTINUE
C
C  Compute the weights for the leading resonance structures, assuming
C  normalization of the weights:
C
      QQ = (ONE - QQ) / NLEAD
      IF(QQ.LT.ZERO) THEN
        FW = BIG
        RETURN
      END IF
      DO 20 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) WGT(IRES,IREF) = QQ
   20 CONTINUE
C
C  Compute the deviations, ERR:
C
      CALL GETERR(ERR(IREF),WGT,Q,NRES,MAXRES,IREF,NREF)
C
C  Add a penalty function to maximize the weight of the reference:
C
      ERR(IREF) = ERR(IREF) + YLAM * (ONE - WGT(1,IREF))
C
C  Compute the fractional reduction for this set of resonance weights:
C
      FW = ONE - ERR(IREF) / ERRREF(IREF)
C
C  Finally, add a penalty that tends to maximize a handful of leading
C  structures:
C
      IF(XLAM.GT.ZERO) THEN
        PENLTY = ZERO
        DO 30 IRES = 1,NRES(IREF)
          PENLTY = PENLTY + WGT(IRES,IREF) * WGT(IRES,IREF)
   30   CONTINUE
        FW = FW - XLAM * (ONE - PENLTY)
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE EXPWGT(RHOW,RHOSTR,WGT,IRANK,NRES,MAXRES,IREF,NREF,N)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION RHOSTR(MAXRES,NREF),WGT(MAXRES,NREF),IRANK(MAXRES),
     +          NRES(NREF)
C
      SAVE ZERO,EPS,ONE,HUNDRD
      DATA ZERO,EPS,ONE,HUNDRD/0.0D0,5.0D-5,1.0D0,1.0D2/
C
C  Compute the unnormalized resonance weights:
C
      IF(RHOW.EQ.ZERO) THEN
        RHO = HUNDRD
        DO 5 IRES = 1,NRES(IREF)
          IF(RHOSTR(IRES,IREF).GE.ZERO) THEN
            IF(RHOSTR(IRES,IREF).LT.RHO) THEN
              RHO = RHOSTR(IRES,IREF)
            END IF
          END IF
    5   CONTINUE
        DO 10 IRES = 1,NRES(IREF)
          IF(ABS(RHOSTR(IRES,IREF)-RHO).LE.EPS) THEN
            WGT(IRES,IREF) = ONE
          ELSE
            WGT(IRES,IREF) = ZERO
          END IF
   10   CONTINUE
      ELSE
        FACTOR = ONE / RHOW
        DO 20 IRES = 1,NRES(IREF)
          IF(RHOSTR(IRES,IREF).GE.ZERO) THEN
            WGT(IRES,IREF) = EXP(-FACTOR * RHOSTR(IRES,IREF))
          ELSE
            WGT(IRES,IREF) = ZERO
          END IF
   20   CONTINUE
      END IF
C
C  Eliminate unwanted structures
C
      DO 30 IRES = N+1,NRES(IREF)
        WGT(IRANK(IRES),IREF) = ZERO
   30 CONTINUE
C
C  Compute the normalization Q:
C
      Q = ZERO
      DO 40 IRES = 1,NRES(IREF)
        Q = Q + WGT(IRES,IREF)
   40 CONTINUE
C
C  If Q is zero, halt program execution:
C
      IF(Q.EQ.ZERO) THEN
        CALL NBHALT('All resonance weights vanish in SR EXPWGT.')
C
C  Otherwise, normalize the resonance weights:
C
      ELSE
        Q = ONE / Q
        DO 50 IRES = 1,NRES(IREF)
          WGT(IRES,IREF) = WGT(IRES,IREF) * Q
   50   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE FETREF(T,LABEL,IBXM,IREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXBAS = 2000)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION T(NDIM,NDIM),LABEL(MAXBAS,6),IBXM(MAXBAS),SCR(MAXBAS,7),
     +          LSCR(MAXBAS,7)
      EQUIVALENCE (LSCR(1,1),SCR(1,1))
C
C  FETREF:  FETCH THE NAO-NBO TRANSFORMATION AND LABELS FOR REFERENCE
C           STRUCTURE IREF
C
      NFILE = 2 * IREF + 99
      L3 = NDIM*NDIM
      CALL NBREAD(T,L3,NFILE)
      NFILE = NFILE + 1
      CALL NBREAD(SCR,MAXBAS*7,NFILE)
      DO 20 J = 1,6
        DO 10 I = 1,MAXBAS
          LABEL(I,J) = LSCR(I,J)
   10   CONTINUE
   20 CONTINUE
      DO 30 I = 1,MAXBAS
        IBXM(I) = LSCR(I,7)
   30 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE FORMQ(Q,OCC,IST,NRES,MAXRES,IREF,NREF,ICTRL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION Q(MAXRES,MAXRES),OCC(NDIM,MAXRES,NREF)
      DIMENSION IST(NDIM,NREF),NRES(NREF)
C
      SAVE ZERO,TWO
      DATA ZERO,TWO/0.0D0,2.0D0/
C
C  If ICTRL is zero, then sum over all orbitals:
C
      IF(ICTRL.EQ.0) THEN
        ICNT = NNAO
        DO 30 JRES = 1,NRES(IREF)
          DO 20 IRES = JRES,NRES(IREF)
            Q(IRES,JRES) = ZERO
            DO 10 IBAS = 1,NNAO
              Q(IRES,JRES) = Q(IRES,JRES) +
     +                       OCC(IBAS,IRES,IREF) * OCC(IBAS,JRES,IREF)
   10       CONTINUE
   20     CONTINUE
   30   CONTINUE
C
C  If ICTRL is one, then only sum over valence orbitals:
C
      ELSE IF(ICTRL.EQ.1) THEN
        ICNT = 0
        DO 40 IBAS = 1,NNAO
          IF(IST(IBAS,IREF).EQ.1) ICNT = ICNT + 1
   40   CONTINUE
        DO 70 JRES = 1,NRES(IREF)
          DO 60 IRES = JRES,NRES(IREF)
            Q(IRES,JRES) = ZERO
            DO 50 IBAS = 1,NNAO
              IF(IST(IBAS,IREF).EQ.1) THEN
                Q(IRES,JRES) = Q(IRES,JRES) +
     +                         OCC(IBAS,IRES,IREF) * OCC(IBAS,JRES,IREF)
              END IF
   50       CONTINUE
   60     CONTINUE
   70   CONTINUE
C
C  If ICTRL is two, then only sum over core and valence orbitals:
C
      ELSE IF(ICTRL.EQ.2) THEN
        ICNT = 0
        DO 80 IBAS = 1,NNAO
          IF(IST(IBAS,IREF).EQ.0.OR.IST(IBAS,IREF).EQ.1) ICNT = ICNT + 1
   80   CONTINUE
        DO 110 JRES = 1,NRES(IREF)
          DO 100 IRES = JRES,NRES(IREF)
            Q(IRES,JRES) = ZERO
            DO 90 IBAS = 1,NNAO
              IF(IST(IBAS,IREF).EQ.0.OR.IST(IBAS,IREF).EQ.1) THEN
                Q(IRES,JRES) = Q(IRES,JRES) +
     +                         OCC(IBAS,IRES,IREF) * OCC(IBAS,JRES,IREF)
              END IF
   90       CONTINUE
  100     CONTINUE
  110   CONTINUE
      END IF
C
C  Scale Q and double off-diagonal elements:
C
      DO 130 JRES = 1,NRES(IREF)
        Q(JRES,JRES) = Q(JRES,JRES) / ICNT
        DO 120 IRES = JRES+1,NRES(IREF)
          Q(IRES,JRES) = TWO * Q(IRES,JRES) / ICNT
          Q(JRES,IRES) = Q(IRES,JRES)
  120   CONTINUE
  130 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETDM(T,DM,TR,SCR1,SCR2,SCR3,IST,IOCC,NVAL,IERR)
C***********************************************************************
C  2-Jun-98  EDG  Print orthogonalization error with NRTDTL [JPRINT(57)]
C-----------------------------------------------------------------------
C
C  Calculate the valence NAO density matrix (DM):
C
C  Input data:
C
C        T:  Full transformation
C      IST:  Orbital type index (0:core, 1:valence, 2:extra-valence)
C     IOCC:  Orbital occupancy index (-1:unoccupied, 1:occupied) (destroyed)
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXBAS = 2000)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBNAO/NAOC(MAXBAS),NAOA(MAXBAS),LTYP(MAXBAS),IPRIN(MAXBAS)
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
C
      DIMENSION T(NDIM,NDIM),DM(NVAL*(NVAL+1)/2),TR(NVAL,NVAL),
     + SCR1(NVAL,NVAL),SCR2(NVAL,NVAL),SCR3(NVAL),IST(NDIM),IOCC(NDIM)
C
      SAVE ZERO,ONE,TWO,ISTR
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
      DATA ISTR/3HVal/
C
C  Modify IST index (to be restored later):
C
      DO 10 I = 1,NNAO
        IST(I) = IST(I) * IOCC(I)
   10 CONTINUE
C
C  Setup the truncated valence NAO-NBO transformation:
C
      IVAL = 0
      NOCC = 0
      DO 30 I = 1,NNAO
        IF(ABS(IST(I)).EQ.1) THEN
          IVAL = IVAL + 1
          IF(IST(I).EQ.1) THEN
            NOCC = NOCC + 1
            IOCC(NOCC) = IVAL
          END IF
          JVAL = 0
          DO 20 J = 1,NNAO
            IF(LTYP(J).EQ.ISTR) THEN
              JVAL = JVAL + 1
              TR(JVAL,IVAL) = T(J,I)
            END IF
   20     CONTINUE
        END IF
   30 CONTINUE
C
C  Symmetrically orthogonalize TR:
C
      IERR  = 0
      IPFLG = 0
      IF(JPRINT(57).NE.0) IPFLG = 1
      CALL SYMORT(SCR1,TR,SCR2,NVAL,NVAL,SCR3,IERR,IPFLG)
      IF(IERR.GT.0) RETURN
C
C  Calculate valence NAO density:
C
      ETA = TWO
      IF(ISPIN.NE.0) ETA = ONE
      IJ = 0
      DO 60 J = 1,NVAL
        DO 50 I = 1,J
          IJ = IJ + 1
          DMIJ = ZERO
          DO 40 K = 1,NOCC
            L = IOCC(K)
            DMIJ = DMIJ + ETA * TR(I,L) * TR(J,L)
   40     CONTINUE
          DM(IJ) = DMIJ
   50   CONTINUE
   60 CONTINUE
C
C  Restore IST:
C
      DO 70 I = 1,NNAO
        IST(I) = ABS(IST(I))
   70 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETDW(DW,GAMMAW,WGT,LVAL,NREF,NVAL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C     COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
C    +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
C    +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      DIMENSION GAMMAW(LVAL,NREF),WGT(NREF)
C
      SAVE ZERO,TWO
      DATA ZERO,TWO/0.0D0,2.0D0/
C
C  Calculate error (doubling off-diagonal elements):
C
      IJ = 0
      DW = ZERO
      DO 30 J = 1,NVAL
        DO 20 I = 1,J
          IJ = IJ + 1
          DW2 = ZERO
          DO 10 IREF = 1,NREF
            IF(WGT(IREF).GT.ZERO) THEN
              TMP = WGT(IREF) * GAMMAW(IJ,IREF)
              DW2 = DW2 + TMP
            END IF
   10     CONTINUE
          IF(I.NE.J) THEN
            DW = DW + TWO * DW2 * DW2
          ELSE
            DW = DW + DW2 * DW2
          END IF
   20   CONTINUE
   30 CONTINUE
      DW = SQRT(DW / NVAL / NVAL)
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETERR(ERROR,WGT,Q,NRES,MAXRES,IREF,NREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION WGT(MAXRES,NREF),Q(MAXRES,MAXRES)
      DIMENSION NRES(NREF)
C
      SAVE ZERO
      DATA ZERO/0.0D0/
C
      ERROR = ZERO
      DO 20 JRES = 1,NRES(IREF)
        TMP = ZERO
        DO 10 IRES = 1,JRES
          TMP = TMP + WGT(IRES,IREF) * Q(IRES,JRES)
   10   CONTINUE
        ERROR = ERROR + WGT(JRES,IREF) * TMP
   20 CONTINUE
      ERROR = SQRT(ERROR)
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETFDM(DM,DELOC,TNAO,TBO,SCR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION DM(NDIM,NDIM),DELOC(NDIM,NDIM),TNAO(NDIM,NDIM),
     +          TBO(NDIM,NDIM),SCR(NDIM)
C
C  Fetch the AO Fock matrix and transform it to the BO basis:
C
      CALL FEFAO(DELOC,IWFOCK)
      IF(IWFOCK.EQ.1) THEN
        CALL FETNAO(TNAO)
        CALL SIMTR1(DELOC,TNAO,SCR,NDIM,NBAS,NBAS,NNAO)
        CALL SIMTRS(DELOC,TBO,SCR,NDIM,NNAO)
C
C  If the Fock matrix is not available, store the BO density in DELOC:
C
      ELSE
        CALL COPY(DM,DELOC,NDIM,NBAS,NBAS)
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETGRD(DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,IREF,
     +                  NREF,NVAR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION DER(MAXRES),T(MAXRES),WGT(MAXRES,NREF),Q(MAXRES,MAXRES),
     +          OCC(NDIM,MAXRES,NREF),G(MAXRES),IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO
      DATA ZERO/0.0D0/
C
C  Compute the derivative with respect to the resonance weights:
C
C  Count the number of parent structures (IM):
C
      IM = 0
      DO 10 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) IM = IM + 1
   10 CONTINUE
C
C  Zero derivatives and compute G factors:
C
      DO 30 IV = 1,NVAR
        DER(IV) = ZERO
        IN = 0
        DO 20 IRES = 1,NRES(IREF)
          IF(IPTR(IRES).EQ.IV) IN = IN + 1
   20   CONTINUE
        G(IV) = DFLOAT(IN) / DFLOAT(IM)
   30 CONTINUE
C
C  Loop over orbitals, evaluating contributions to the derivatives:
C
      DO 70 IB = 1,NNAO
C
C  Compute DELTA, the difference between the resonance-weighted
C  occupancy and the actual occupancy of the orbital IB:
C
        DELTA = ZERO
        DO 40 IRES = 1,NRES(IREF)
          DELTA = DELTA + WGT(IRES,IREF) * OCC(IB,IRES,IREF)
   40   CONTINUE
C
C  If this difference is zero, we can skip the rest:
C
        IF(DELTA.NE.ZERO) THEN
C
C  Loop over the derivatives:
C
          DO 60 IV = 1,NVAR
C
C  Evaluate the contribution to the derivative from each resonance
C  structure:
C
            TEMP = ZERO
            DO 50 IRES = 1,NRES(IREF)
              IF(IPTR(IRES).EQ.0) THEN
                TEMP = TEMP + G(IV) * OCC(IB,IRES,IREF)
              ELSE IF(IPTR(IRES).EQ.IV) THEN
                TEMP = TEMP - OCC(IB,IRES,IREF)
              END IF
   50       CONTINUE
            DER(IV) = DER(IV) + DELTA * TEMP
   60     CONTINUE
        END IF
   70 CONTINUE
C
C  Finish derivatives:
C
      CALL GETERR(ERROR,WGT,Q,NRES,MAXRES,IREF,NREF)
      DO 80 IV = 1,NVAR
        IF(DER(IV).EQ.ZERO.OR.ERROR.EQ.ZERO) THEN
          DER(IV) = ZERO
        ELSE
          DER(IV) = DER(IV) * EXP(-T(IV)) / (ERROR * NNAO)
        END IF
   80 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETHES(HESS,DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,
     +                  IREF,NREF,NVAR,IFLG)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION HESS(MAXRES,MAXRES),DER(MAXRES),T(MAXRES),
     +          WGT(MAXRES,NREF),Q(MAXRES,MAXRES),OCC(NDIM,MAXRES,NREF),
     +          G(MAXRES),IPTR(MAXRES),NRES(NREF)
C
      SAVE ZERO,THRESH
      DATA ZERO,THRESH/0.0D0,1.0D-9/
C
C  Compute the hessian (the second derivatives with respect to the
C  resonance weights):
C
C  First, evaluate the function value and gradients at T:
C
      CALL GETERR(ERROR,WGT,Q,NRES,MAXRES,IREF,NREF)
      CALL GETGRD(DER,T,WGT,Q,OCC,G,IPTR,NRES,MAXRES,IREF,NREF,NVAR)
C
C  If error is too small, the second derivatives are ill-defined:
C
      IFLG = 0
      IF(ERROR.LT.THRESH) THEN
        IFLG = 1
        RETURN
      END IF
C
C  Count the number of parent structures (IM):
C
      IM = 0
      DO 10 IRES = 1,NRES(IREF)
        IF(IPTR(IRES).EQ.0) IM = IM + 1
   10 CONTINUE
C
C  Zero matrix elements and compute G factors:
C
      DO 40 IV = 1,NVAR
        DO 20 JV = 1,NVAR
          HESS(JV,IV) = ZERO
   20   CONTINUE
        IN = 0
        DO 30 IRES = 1,NRES(IREF)
          IF(IPTR(IRES).EQ.IV) IN = IN + 1
   30   CONTINUE
        G(IV) = DFLOAT(IN) / DFLOAT(IM)
   40 CONTINUE
C
C  Compute the hessian matrix elements:
C
      DO 90 JV = 1,NVAR
        DO 80 IV = JV,NVAR
          HESS(IV,JV) = ZERO
          DO 70 IB = 1,NNAO
            TEMP = ZERO
            DO 50 IRES = 1,NRES(IREF)
              IF(IPTR(IRES).EQ.0) THEN
                TEMP = TEMP + G(IV) * OCC(IB,IRES,IREF)
              ELSE IF(IPTR(IRES).EQ.IV) THEN
                TEMP = TEMP - OCC(IB,IRES,IREF)
              END IF
   50       CONTINUE
            TEMP1 = ZERO
            DO 60 IRES = 1,NRES(IREF)
              IF(IPTR(IRES).EQ.0) THEN
                TEMP1 = TEMP1 + G(JV) * OCC(IB,IRES,IREF)
              ELSE IF(IPTR(IRES).EQ.JV) THEN
                TEMP1 = TEMP1 - OCC(IB,IRES,IREF)
              END IF
   60       CONTINUE
            HESS(IV,JV) = HESS(IV,JV) + TEMP * TEMP1
   70     CONTINUE
          HESS(IV,JV) = HESS(IV,JV) * EXP(-T(IV)-T(JV))
          HESS(IV,JV) = HESS(IV,JV) / (NNAO * ERROR)
   80   CONTINUE
   90 CONTINUE
C
      DO 110 JV = 1,NVAR
        HESS(JV,JV) = HESS(JV,JV) - DER(IV)
        DO 100 IV = JV,NVAR
          HESS(IV,JV) = HESS(IV,JV) - DER(IV) * DER(JV) / ERROR
          HESS(JV,IV) = HESS(IV,JV)
  100   CONTINUE
  110 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE GETREL(IBO,IBOP,IUNIT1,IUNIT2,IREL,IA,IB,IC,ID)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  PROGRAM TO PASS RELATIONSHIP INDEX (IREL) BETWEEN BOND ORBITALS IBO, IBOP,
C  LISTING ATOMS (IA,IB,IC,ID) IN LINKED ORDER: IBO=(IA,IB), IBOP=(IC,ID).
C   IREL = 0 IF SEPARATED BY MORE THAN ONE LINKING BOND,
C   IREL = 1 IF GEMINAL, IA-(IB=IC)-ID,
C   IREL = 2 IF VICINAL, IA-IB-IC-ID.
C   IREL = 3 IF DOUBLE-VICINAL (A-B-C-D OR B-A-C-D, TRIANGLE A-B-C*)
C   IREL = 4 IF DOUBLE-VICINAL (A-B-C-D OR A-B-D-C, TRIANGLE B*-C-D)
C   IREL = 5 IF DOUBLE-VICINAL (A-B-C-D OR B-A-D-C, SQUARE A*-B-C-D*)
C  IA = 0 IF IBO IS A LONE PAIR; ID = 0 IF IBOP IS A LONE PAIR.
C  [UNIT NUMBERS (IUNIT1,IUNIT2) ARE NOT CURRENTLY IMPLEMENTED.]
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBL(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
      COMMON/NBTOPO/IORDER(MAXATM),JORDER(MAXATM),NTOPO(MAXATM,MAXATM),
     +            N3CTR,I3CTR(10,3)
C
      IREL = 0
      IUNIT1 = NBOUNI(IBO)
      IUNIT2 = NBOUNI(IBOP)
      IA = LABEL(IBXM(IBO),4)
      IB = LABEL(IBXM(IBO),5)
      IC = LABEL(IBXM(IBOP),4)
      ID = LABEL(IBXM(IBOP),5)
C
C  IBO IS A LONE PAIR?
C
      IF(IA.NE.0.AND.IB.EQ.0) THEN
        IB = IA
        IA = 0
C
C  GEMINAL?
C
        IF(IB.EQ.IC) THEN
          IREL = 1
          RETURN
        ELSE IF(IB.EQ.ID) THEN
          IT = IC
          IC = ID
          ID = IT
          IREL = 1
          RETURN
C
C  VICINAL?
C
        ELSE IF(ID.EQ.0) THEN
          IF(NTOPO(IB,IC).NE.0) THEN
            IREL = 2
          ELSE
            IREL = 0
          END IF
          RETURN
        ELSE IF(NTOPO(IB,IC).NE.0) THEN
          IREL = 2
          IF(NTOPO(IB,ID).NE.0) IREL = 4
          RETURN
        ELSE IF(NTOPO(IB,ID).NE.0) THEN
          IT = IC
          IC = ID
          ID = IT
          IREL = 2
          RETURN
        ELSE
          IREL = 4
          RETURN
        END IF
      END IF
C
C  OR IBO IS A BOND
C  GEMINAL?
C
      IF(IA.EQ.IC.OR.IA.EQ.ID) THEN
        IT = IA
        IA = IB
        IB = IT
        IF(IB.EQ.ID) THEN
          IT = IC
          IC = ID
          ID = IT
        END IF
        IREL = 1
      ELSE IF(IB.EQ.IC) THEN
        IREL = 1
      ELSE IF(IB.EQ.ID) THEN
        IT = IC
        IC = ID
        ID = IT
        IREL = 1
      END IF
C
C  SPECIAL GEMINAL CASE: IF IA AND ID ARE BONDED, SWITCH RELATIONSHIP
C  TO VICINAL
C
      IF(IB.EQ.IC.AND.IREL.EQ.1) THEN
        IF(IA.NE.ID.AND.IA.NE.0.AND.ID.NE.0) THEN
          IF(NTOPO(IA,ID).GT.0) THEN
            IT = IA
            IA = IB
            IB = IT
            IT = IC
            IC = ID
            ID = IT
            IREL = 2
          END IF
        END IF
        RETURN
      END IF
C
C  VICINAL?
C
      IF(ID.EQ.0) THEN
        IF(NTOPO(IA,IC).NE.0) THEN
          IT = IA
          IA = IB
          IB = IT
          IREL = 2
          IF(NTOPO(IA,IC).NE.0) IREL = 3
        ELSE IF(NTOPO(IB,IC).NE.0) THEN
          IREL = 2
        ELSE
          IREL = 3
        END IF
        RETURN
      END IF
      IF(NTOPO(IA,IC).NE.0.OR.NTOPO(IA,ID).NE.0) THEN
        IT = IA
        IA = IB
        IB = IT
        IREL = 2
        IF(NTOPO(IB,ID).NE.0) THEN
          IT = IC
          IC = ID
          ID = IT
        END IF
      ELSE IF(NTOPO(IB,IC).NE.0) THEN
        IREL = 2
      ELSE IF(NTOPO(IB,ID).NE.0) THEN
        IT = IC
        IC = ID
        ID = IT
        IREL = 2
      END IF
C
C  HIGHER CONNECTIVITY?
C
      IF(IREL.EQ.2) THEN
        IF(NTOPO(IA,IC).NE.0) IREL = 3
        IF(NTOPO(IB,ID).NE.0) IREL = 4
        IF(NTOPO(IA,ID).NE.0) IREL = 5
      END IF
      RETURN
      END
C***********************************************************************
      FUNCTION INTRN(ILLIM,IULIM,ISEED)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Compute a random integer between ILLIM and IULIM (inclusive):
C
      ITMP = ILLIM + INT((IULIM - ILLIM + 1) * RANNB(ISEED))
      INTRN = MIN(ITMP,IULIM)
      RETURN
      END
C***********************************************************************
      SUBROUTINE MULTI(FW,GAMMAW,WGT,SCR,IRESET,LVAL,NREF,NVAL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL SKIP
C
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION GAMMAW(LVAL,NREF),WGT(NREF),SCR(NREF),IRESET(NREF)
      DIMENSION W(5),IR(5),IC(5)
C
      SAVE ZERO,THRFW,THRS,PT5,ONE,THREE,TEN
      SAVE NPAR,NPARO,MAXIT,IBLNK,IPLUS
      DATA ZERO,THRFW,THRS,PT5,ONE,THREE,TEN/0.0D0,1.0D-5,2.0D-5,0.5D0,
     +     1.0D0,3.0D0,1.0D1/
      DATA NPAR,MAXIT/0,1000/
      DATA IBLNK,IPLUS/1H ,1H+/
C
C  Determine relative weights for each of the reference manifolds:
C
C  Get list of optimizable parameters:
C
      IP    = 0
      NP    = 0
      NPARO = NPAR
      NPAR  = 0
      DO 5 IREF = NREF,1,-1
        IF(IRESET(IREF).GE.0) THEN
          NPAR = NPAR + 1
          IF(IP.EQ.0) THEN
            IP = IREF
          ELSE IF(NP.EQ.0) THEN
            NP = IREF
          END IF
        END IF
    5 CONTINUE
      IF(NPAR.EQ.0) RETURN
C
C  Skip multi-ref optimization?
C
      IF(NPAR.EQ.1) THEN
        WGT(IP) = ONE
        RETURN
      END IF
      ISKIP = 0
      IF(NPAR.EQ.NPARO) ISKIP = 1
      DO 10 IREF = 1,NREF
        IF(IRESET(IREF).GT.0) ISKIP = 0
   10 CONTINUE
      IF(ISKIP.EQ.1) RETURN
      IF(JPRINT(57).NE.0) WRITE(LFNPR,890)
C
C  Determine the reference DWREF:
C
      IPP = 1
      DWREF = TEN
      DO 15 IREF = 1,IP
        WGT(IREF) = ZERO
   15 CONTINUE
      DO 20 IREF = 1,IP
        IF(IRESET(IREF).GE.0) THEN
          WGT(IREF) = ONE
          CALL GETDW(DW1,GAMMAW,WGT,LVAL,NREF,NVAL)
          IF(DW1.LT.DWREF) THEN
            DWREF = DW1
            CALL COPY(WGT,SCR,NREF,NREF,1)
            IPP = IREF
          END IF
          WGT(IREF) = ZERO
        END IF
   20 CONTINUE
      WGT(IPP) = ONE
C
C  Set a few initial parameters.  These will be changed throughout the
C  optimization:
C
      ISEED = -ABS(JPRINT(55)) - NREF
      STHR = 0.1D0
      STEP = STHR
      STEPMX = ZERO
      GLOBAL = DW1
C
C  Take a few iterations to determine an appropriate ``temperature''
C  DT to start with:
C
      DT = ZERO
      DO 40 I = 1,10*NPAR
C
C  Make an arbitrary adjustment to one of the weights and compute DW:
C
   30   CALL ADJUST(WGT,IW,NP,SSIZE,STEP,ISEED)
        CALL WNORM(WGT,IRESET,IW,IP,SKIP)
        IF(SKIP) THEN
          CALL RESTOR(WGT,IW,NP,SSIZE)
          GOTO 30
        END IF
        CALL GETDW(DW2,GAMMAW,WGT,LVAL,NREF,NVAL)
        STEPMX = MAX(STEPMX,ABS(SSIZE))
        CALL RESTOR(WGT,IW,NP,SSIZE)
C
C  Keep track of the largest change in DW, which will determine DT:
C
        DT = MAX(DT,ABS(DW2-DW1))
        DT = MIN(DT,PT5)
   40 CONTINUE
      FW1 = ONE - DW1 / DWREF
      IF(JPRINT(57).NE.0) WRITE(LFNPR,900) FW1,DT,STEPMX
C
C  Take up to MAXIT grand iterations:
C
      ITOT = 0
      NTST = 10 * NPAR
      DO 100 IG = 1,MAXIT
C
C  Loop over micro iterations.  On each iteration, adjust one of the
C  resonance weights, compute the new DW, and decide whether the adjustment
C  is acceptable (NSUC.GT.0).  If not, restore the previous weights.  Also,
C  keep track of the largest acceptable adjustment made and the largest and
C  smallest accepted DW:
C
        ITT = IG
        ITER = 0
        ICNT = 0
        NSUC = 0
        DWMAX = DW1
        DWMIN = DW1
        STEPMX = ZERO
        DO 60 IM = 1,100*NPAR
          ITER = ITER + 1
   50     CALL ADJUST(WGT,IW,NP,SSIZE,STEP,ISEED)
          CALL WNORM(WGT,IRESET,IW,IP,SKIP)
          IF(SKIP) THEN
            CALL RESTOR(WGT,IW,NP,SSIZE)
            GOTO 50
          END IF
          CALL GETDW(DW2,GAMMAW,WGT,LVAL,NREF,NVAL)
          CALL METROP(NSUC,DW1,DW2,DT,ISEED)
          IF(NSUC.EQ.0) THEN
            CALL RESTOR(WGT,IW,NP,SSIZE)
          ELSE
            ICNT = ICNT + 1
            DW1 = DW2
            DWMAX = MAX(DWMAX,DW1)
            DWMIN = MIN(DWMIN,DW1)
            STEPMX = MAX(STEPMX,ABS(SSIZE))
            IF(GLOBAL.GT.DWMIN) THEN
              ITG = IG
              GLOBAL = DWMIN
              CALL COPY(WGT,SCR,NREF,NREF,1)
            END IF
            IF(NSUC.EQ.NTST) GOTO 70
          END IF
   60   CONTINUE
C
C  Anneal Scheduling - Adjust ``temperature'' and maximum step size:
C
   70   CONTINUE
        DT = 0.9D0 * DT
        STEP = 1.5D0 * STEPMX
C
C  Print info about optimization:
C
        ITOT = ITOT + ITER
        FW1   = ONE - DW1   / DWREF
        FWMAX = ONE - DWMIN / DWREF
        FWMIN = ONE - DWMAX / DWREF
        IF(JPRINT(57).NE.0) WRITE(LFNPR,910) IG,ITER,ICNT,ITOT,FW1,
     +                                       FWMAX,FWMIN,DT,STEPMX
C
C  Has the optimization converged?
C
        IF(DWMAX-DWMIN.LT.THRFW) THEN
C
C  If so, check to see that we have the best DW found so far.  If not,
C  restart the optimization:
C
          IF(ABS(DWMIN-GLOBAL).GT.THRFW) THEN
            CALL COPY(SCR,WGT,NREF,NREF,1)
            STHR = STHR / TEN
            STEP = MIN(STHR,THREE*STEP)
            DT = ZERO
            CALL GETDW(DW1,GAMMAW,WGT,LVAL,NREF,NVAL)
            FW1 = ONE - DW1 / DWREF
            IF(JPRINT(57).NE.0) WRITE(LFNPR,920) ITG,FW1
            DO 90 I = 1,10*NPAR
   80         CALL ADJUST(WGT,IW,NP,SSIZE,STEP,ISEED)
              CALL WNORM(WGT,IRESET,IW,IP,SKIP)
              IF(SKIP) THEN
                CALL RESTOR(WGT,IW,NP,SSIZE)
                GOTO 80
              END IF
              CALL GETDW(DW2,GAMMAW,WGT,LVAL,NREF,NVAL)
              CALL RESTOR(WGT,IW,NP,SSIZE)
              DT = MAX(DT,ABS(DW2-DW1))
              DT = MIN(DT,PT5)
   90       CONTINUE
C
C  Otherwise, abort the iterations.  But first make sure that STEPMX has
C  converged:
C
          ELSE IF(STEPMX.LT.THRS) THEN
            GOTO 110
          END IF
        END IF
  100 CONTINUE
      ITT = -MAXIT
C
C  We're finished (at last)!!  Recompute DW:
C
  110 CONTINUE
      CALL WNORM(WGT,IRESET,IW,IP,SKIP)
      CALL GETDW(DW1,GAMMAW,WGT,LVAL,NREF,NVAL)
C
C  Write optimized weights:
C
      IF(JPRINT(57).NE.0) THEN
        WRITE(LFNPR,930)
        IFLG = 0
        IREF = 0
  120   DO 130 I = 1,5
          IC(I) = IBLNK
  130   CONTINUE
        ICNT = 0
  140   IREF = IREF + 1
          IF(IRESET(IREF).GE.0) THEN
            ICNT = ICNT + 1
            IF(IFLG.EQ.1) IC(ICNT) = IPLUS
            IFLG = 1
            W(ICNT) = WGT(IREF)
            IR(ICNT) = IREF
          END IF
        IF(IREF.LT.NREF.AND.ICNT.LT.5) GOTO 140
        WRITE(LFNPR,940) (IC(I),W(I),IR(I),I=1,ICNT)
        IF(IREF.LT.NREF) GOTO 120
        WRITE(LFNPR,*)
      END IF
      FW = ONE - DW1 / DWREF
      IF(ITT.GT.0) THEN
        WRITE(LFNPR,950) NPAR,DW1,FW,ITT
      ELSE
        WRITE(LFNPR,960) NPAR,DW1,FW,ABS(ITT)
      END IF
      RETURN
C
  890 FORMAT(/1X,'Variational optimization of the multi-reference ',
     + 'weights (SR MULTI):')
  900 FORMAT(/1X,' Grand  Micro          Total    ',
     + '         max      min       dt      max',/1X,'  iter   iter',
     + '   Acc     iter    F(W)     F(W)     F(W)    (temp)    step',
     + /1X,'-------------------------------------------------------',
     + '-------------------',/1X,'    0',23X,F9.5,18X,2F9.5)
  910 FORMAT(1X,I5,I7,I7,I8,1X,5F9.5)
  920 FORMAT(1X,'        Restarting optimization at iteration',
     + I4,', (Fw =',F8.5,')')
  930 FORMAT(/1X,'Optimized multi-reference weights:')
  940 FORMAT(1X,5(A1,F8.5,'(',I4,')',:,1X))
  950 FORMAT(1X,'Multi-ref(',I2,'):  D(W)=',F7.5,', F(W)=',F7.5,
     + ' converged after',I4,' iterations')
  960 FORMAT(1X,'Multi-ref(',I2,'):  D(W)=',F7.5,', F(W)=',F7.5,
     + ' did not converge after',I4,' iterations')
      END
C***********************************************************************
      SUBROUTINE NRTINP(IESS,INRT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL END,EQUAL
      DIMENSION KEYWD(8),KNRTA(8),KNRTB(8),KDEL(4),KNBO(4),KCRD(6)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      SAVE KDEL,KNBO,KNRTA,KNRTB,KCRD
      DATA KDEL/1H$,1HD,1HE,1HL/,KNBO/1H$,1HN,1HB,1HO/,
     + KNRTA/1H$,1HN,1HR,1HT,1HS,1HT,1HR,1HA/,
     + KNRTB/1H$,1HN,1HR,1HT,1HS,1HT,1HR,1HB/,
     + KCRD/1H$,1HC,1HO,1HO,1HR,1HD/
C
C  If this is the GAMESS, HONDO, or general version of the NBO program,
C  rewind the input file before searching for $NRTSTR:
C
      IREP = 1
      IF(IESS.EQ.0) IREP = 0
      IF(IESS.EQ.6) IREP = 0
      IF(IESS.EQ.7) IREP = 0
      IF(IREP.EQ.0) REWIND(LFNIN)
C
C  Search for the $NRTSTR keylist:
C
   10 CALL STRTIN(LFNIN)
      LEN = 8
      CALL HFLD(KEYWD,LEN,END)
      IF(EQUAL(KEYWD,KNRTA,4)) GOTO 50
      IF(EQUAL(KEYWD,KNBO,4)) GOTO 60
      IF(EQUAL(KEYWD,KDEL,4)) GOTO 60
      IF(EQUAL(KEYWD,KCRD,6)) GOTO 70
      IF(LEN.EQ.0.AND.END) GOTO 80
      GOTO 10
C
C  $NRTSTR found:
C
   50 CONTINUE
      INRT = 1
      IF(.NOT.OPEN) RETURN
      IF(ALPHA.AND.EQUAL(KEYWD,KNRTA,8)) RETURN
      IF(BETA.AND.EQUAL(KEYWD,KNRTB,8)) RETURN
      INRT = 0
      GOTO 10
C
C  $NBO, $DEL found -- discontinue the search for $NRTSTR (GAUSSIAN, AMPAC)
C                      continue searching for $NRTSTR (GENNBO, GAMESS, HONDO)
C
   60 CONTINUE
      IF(IREP.EQ.0) GOTO 10
      BACKSPACE(LFNIN)
      INRT = 0
      RETURN
C
C  $COORD -- discontinue search for $NRTSTR (GENNBO)
C
   70 CONTINUE
      IF(IESS.EQ.0) THEN
        INRT = 0
        RETURN
      END IF
      GOTO 10
C
C  End of file encountered:
C
   80 CONTINUE
      INRT = 0
      RETURN
      END
C***********************************************************************
      SUBROUTINE PRTWGT(IRES,IR,NUM,WGT,II,IAT,JJ,JAT,IFLG)
C***********************************************************************
C 12-May-02  FAW  3-digit atom numbers in ADDED(REMOVED) strings
C 22-Jun-95  EDG  Initialize ISOLEW, REF, and WRES
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*58 STRING
      CHARACTER*3 ISOLEW,BLANKS,MASK
      CHARACTER*2 CHARAT
      CHARACTER*1 BLANK,COMMA,LPARA,RPARA,DASH,STAR,REF,S1,S2,S3
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      SAVE MRES,IPOS,WRES,REF,LINE,STRING,ISOLEW
      SAVE MAXPOS,BLANK,COMMA,LPARA,RPARA,DASH,STAR,BLANKS,MASK
      SAVE ZERO,HUNDRD
C
      DATA MAXPOS/58/
      DATA BLANK,COMMA,LPARA,RPARA,DASH,STAR/' ',',','(',')','-','*'/
      DATA BLANKS,MASK/'   ','( )'/
      DATA ZERO,HUNDRD/0.0D0,1.0D2/
C
C  If IRES = 0, initialize everything:
C
      IF(IRES.EQ.0) THEN
        MRES = 0
        DO 10 I = 1,MAXPOS
          STRING(I:I) = BLANK
   10   CONTINUE
        IPOS = 0
        LINE = 1
        WRES = ZERO
        REF  = BLANK
        ISOLEW = BLANKS
        RETURN
      END IF
C
C  Force the printing of resonance weights if IRES is negative:
C
      IF(IRES.LT.0) THEN
        IF(MRES.NE.0) THEN
          IF(LINE.EQ.1) THEN
            WRITE(LFNPR,900) MRES,REF,ISOLEW,WRES*HUNDRD,STRING
          ELSE
            WRITE(LFNPR,910) STRING
          END IF
        END IF
        DO 20 I = 1,MAXPOS
          STRING(I:I) = BLANK
   20   CONTINUE
        IPOS = 0
        LINE = 1
        RETURN
      END IF
C
C  Build STRING from info passed from SR NRTOUT.  First, if this is
C  a new IRES, write out the remainder of the info from the previous
C  IRES (stored in MRES):
C
      IF(MRES.EQ.0) THEN
        MRES = IRES
        WRES = WGT
      END IF
      IF(IRES.NE.MRES) THEN
        IF(LINE.EQ.1) THEN
          WRITE(LFNPR,900) MRES,REF,ISOLEW,WRES*HUNDRD,STRING
        ELSE
          WRITE(LFNPR,910) STRING
        END IF
        DO 30 I = 1,MAXPOS
          STRING(I:I) = BLANK
   30   CONTINUE
        IPOS = 0
        LINE = 1
        MRES = IRES
        WRES = WGT
      END IF
C
C  Is this structure a reference in the calculation of the fractional
C  reduction:
C
      REF = BLANK
      IF(IR.NE.0) REF = STAR
      ISOLEW = BLANKS
      IF(NUM.GT.1) THEN
        ISOLEW = MASK
        ISOLEW(2:2) = CHAR(48+NUM)
      END IF
C
C  If IFLG = 0, no new info is to be added to this line:
C
      IF(IFLG.EQ.0) RETURN
C
C  Add comma if this is not the start of a new IRES:
C
      IF(IPOS.GT.0) THEN
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = COMMA
      END IF
C
C  Print the line if it is too long for the new info to be added:
C
      NCTR = 1
      IF(II.NE.JJ) NCTR = 2
      IF(NCTR.EQ.1) THEN
        IF(IFLG.GT.0) JPOS = 6
        IF(IFLG.LT.0) JPOS = 8
      END IF
      IF(NCTR.EQ.2) THEN
        IF(IFLG.GT.0) JPOS = 12
        IF(IFLG.LT.0) JPOS = 14
      END IF
      IF(IPOS+JPOS.GT.MAXPOS-1) THEN
        IF(LINE.EQ.1) THEN
          WRITE(LFNPR,900) MRES,REF,ISOLEW,WRES*HUNDRD,STRING
        ELSE
          WRITE(LFNPR,910) STRING
        END IF
        DO 40 I = 1,MAXPOS
          STRING(I:I) = BLANK
   40   CONTINUE
        IPOS = 0
        LINE = LINE + 1
      ENDIF
C
C  Add the new info to the string:
C
      IPOS = IPOS + 1
      STRING(IPOS:IPOS) = BLANK
      IF(IFLG.LT.0) THEN
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = LPARA
      END IF
      IPOS = IPOS + 2
      STRING(IPOS-1:IPOS) = CHARAT(IAT)
      CALL CONVRT3S(II,S1,S2,S3)
      IPOS = IPOS + 1
      STRING(IPOS:IPOS) = S1
      IPOS = IPOS + 1
      STRING(IPOS:IPOS) = S2
      IPOS = IPOS + 1
      STRING(IPOS:IPOS) = S3
      IF(NCTR.EQ.2) THEN
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = DASH
        IPOS = IPOS + 2
        STRING(IPOS-1:IPOS) = CHARAT(JAT)
        CALL CONVRT3S(JJ,S1,S2,S3)
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = S1
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = S2
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = S3
      END IF
      IF(IFLG.LT.0) THEN
        IPOS = IPOS + 1
        STRING(IPOS:IPOS) = RPARA
      END IF
      RETURN
C
  900 FORMAT(1X,I4,A1,A3,F7.2,2X,A58)
  910 FORMAT(17X,A58)
      END
C***********************************************************************
      FUNCTION RANNB(IDUM)
C***********************************************************************
C
C  Computes a uniform random deviate between 0.0 and 1.0.  Set IDUM
C  to any negative value to initialize or reinitialize the sequence:
C
C  From Numerical Recipes: FUNCTION RAN3(IDUM)
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MBIG=1000000000,MSEED=161803398,MZ=0,FAC=1.0D0/MBIG)
C
      DIMENSION MA(55)
C
      SAVE IFF,INEXT,INEXTP,MA
      DATA IFF/0/
C
      IF(IDUM.LT.0.OR.IFF.EQ.0) THEN
        IFF = 1
        MJ = MSEED - IABS(IDUM)
        MJ = MOD(MJ,MBIG)
        MA(55) = MJ
        MK = 1
        DO I = 1,54
          II = MOD(21*I,55)
          MA(II) = MK
          MK = MJ - MK
          IF(MK.LT.MZ) MK = MK + MBIG
          MJ = MA(II)
        END DO
        DO K = 1,4
          DO I = 1,55
            MA(I) = MA(I) - MA(1+MOD(I+30,55))
            IF(MA(I).LT.MZ) MA(I) = MA(I) + MBIG
          END DO
        END DO
        INEXT = 0
        INEXTP = 31
        IDUM = 1
      END IF
      INEXT = INEXT + 1
      IF(INEXT.EQ.56) INEXT = 1
      INEXTP = INEXTP + 1
      IF(INEXTP.EQ.56) INEXTP = 1
      MJ = MA(INEXT) - MA(INEXTP)
      IF(MJ.LT.MZ) MJ = MJ + MBIG
      MA(INEXT) = MJ
      RANNB = MJ * FAC
      RETURN
      END
C***********************************************************************
      FUNCTION REALRN(XLLIM,XULIM,ISEED)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Compute a real random number between XLLIM and XULIM (inclusive):
C
      REALRN = XLLIM + (XULIM - XLLIM) * RANNB(ISEED)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVTREF(T,LABEL,IBXM,IREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXBAS = 2000)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
      DIMENSION T(NDIM,NDIM),LABEL(MAXBAS,6),IBXM(MAXBAS),SCR(MAXBAS,7),
     +          LSCR(MAXBAS,7)
      EQUIVALENCE (LSCR(1,1),SCR(1,1))
C
C  SVTREF:  SAVE THE NAO-NBO TRANSFORMATION AND LABELS FOR REFERENCE
C           STRUCTURE IREF
C
      NFILE = 2 * IREF + 99
      L3 = NDIM*NDIM
      CALL NBWRIT(T,L3,NFILE)
      DO 20 J = 1,6
        DO 10 I = 1,MAXBAS
          LSCR(I,J) = LABEL(I,J)
   10   CONTINUE
   20 CONTINUE
      DO 30 I = 1,MAXBAS
        LSCR(I,7) = IBXM(I)
   30 CONTINUE
      NFILE = NFILE + 1
      CALL NBWRIT(SCR,MAXBAS*7,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE TOPCMP(IDUP,IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Compare the most recent (IRES) entry with previous entries
C  and delete if duplicate:
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION IDXRES(MAXRES,NREF),LSTRES(LEN)
C
      IF(IRES.LT.1) THEN
        WRITE(LFNPR,900) IRES
        RETURN
      END IF
C
C  Determine the number of elements of LSTRES used to store resonance
C  structure IRES:
C
      NB = IDXRES(IRES,IREF) - 1
   10 NB = NB + 1
      IF(LSTRES(NB).NE.-1) GOTO 10
      NB = NB - IDXRES(IRES,IREF) + 1
C
C  Loop over other resonance structure stored in LSTRES:
C
      IA = 1
      NA = 0
      IDUP = 0
   20 IA = IA + NA
      IF(IA.EQ.IDXRES(IRES,IREF)) RETURN
C
C  How many elements are used to store this resonance structure:
C
      NA = IA - 1
   30 NA = NA + 1
      IF(LSTRES(NA).NE.-1) GOTO 30
      NA = NA - IA + 1
C
C  Compare the two resonance structures:
C
      IF(NA.NE.NB) GOTO 20
      ILOC = IA - 1
      DO 40 I = IDXRES(IRES,IREF),IDXRES(IRES,IREF)+NB-1
        ILOC = ILOC + 1
        IF(LSTRES(ILOC).NE.LSTRES(I)) GOTO 20
   40 CONTINUE
C
C  A duplicate of IRES appears already appears in LSTRES.  Remove
C  the copy:
C
      DO 50 I = IDXRES(IRES,IREF),IDXRES(IRES,IREF)+NB-1
        LSTRES(I) = 0
   50 CONTINUE
      IDXRES(IRES,IREF) = IA
      IDUP = 1
      RETURN
C
  900 FORMAT(/1X,'Variable IRES (=',I4,') is out-of-bounds in SR ',
     + 'TOPCMP.')
      END
C***********************************************************************
      SUBROUTINE TOPGET(IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Fetch the TOPO matrix for resonance structure IRES (reference structure
C  IREF) from LSTRES storage:
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION IDXRES(MAXRES,NREF),LSTRES(LEN)
C
C  Find starting point in LSTRES array:
C
      IF(IDXRES(IRES,IREF).EQ.0) THEN
        WRITE(LFNPR,900) IRES,IREF
        CALL NBHALT('TOPO matrix of structure not in LSTRES array.')
      END IF
      ILOC = IDXRES(IRES,IREF) - 1
C
C  Clear the ITOPO array:
C
      CALL TOPZER(NATOMS)
C
C  Loop over atoms:
C
      DO 20 IA = 1,NATOMS
C
C  Read number of lone pairs (NLP) on atom IA:
C
        ILOC = ILOC + 1
        NLP = LSTRES(ILOC)
        ITOPO(IA,IA) = NLP
C
C  Read number of bond pairs (NBD) on atom IA:
C
        ILOC = ILOC + 1
        NBD = LSTRES(ILOC)
        IF(NBD.EQ.0) GOTO 20
C
C  Insert bond entries into TOPO matrix:
C
        DO 10 IBD = 1,NBD
          ILOC = ILOC + 1
          IB = LSTRES(ILOC)
          ITOPO(IA,IB) = ITOPO(IA,IB) + 1
          ITOPO(IB,IA) = ITOPO(IA,IB)
   10   CONTINUE
   20 CONTINUE
      RETURN
C
  900 FORMAT(/1X,'FATAL ERROR:  The TOPO matrix for resonance ',
     + 'structure',I4,' (reference',I2,')',/1X,'              ',
     + 'does not appear in the LSTRES array.')
      END
C***********************************************************************
      SUBROUTINE TOPOUT
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      SAVE DASHES,MAXCOL
      DATA DASHES/4H ---/
      DATA MAXCOL/17/
C
C  Print the TOPO matrix:
C
      NCU = 0
      NLOOPS = (NATOMS - 1) / MAXCOL + 1
      DO 20 L = 1,NLOOPS
        NCL = NCU + 1
        NCU = MIN(NCU+MAXCOL,NATOMS)
        WRITE(LFNPR,900) (J,J=NCL,NCU)
        WRITE(LFNPR,910) (DASHES,J=NCL,NCU)
        DO 10 I = 1,NATOMS
          WRITE(LFNPR,920) I,NAMEAT(IATNO(I)),
     +                      (ITOPO(I,K),K=NCL,NCU)
   10   CONTINUE
   20 CONTINUE
      RETURN
C
  900 FORMAT(/1X,'    Atom',17(I3,1X))
  910 FORMAT(1X,'    ----',17A4)
  920 FORMAT(1X,I3,'. ',A2,17I4)
      END
C***********************************************************************
      SUBROUTINE TOPSTR(IRES,IDXRES,MAXRES,IREF,NREF,LSTRES,LEN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Store the full TOPO matrix (NATOMS x NATOMS) in condensed
C  form in vector LSTRES, with entries indexed by IDXRES
C
C  IDXRES(I,J) = First element of LSTRES for resonance structure I
C                (reference structure J)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION IDXRES(MAXRES,NREF),LSTRES(LEN)
C
      SAVE ILOC
C
C  Initialize the arrays if IRES<0:
C
      IF(IRES.LT.0) THEN
        ILOC = 0
        DO 20 J = 1,NREF
          DO 10 I = 1,MAXRES
            IDXRES(I,J) = 0
   10     CONTINUE
   20   CONTINUE
        DO 30 I = 1,LEN
          LSTRES(I) = 0
   30   CONTINUE
        RETURN
      END IF
C
C  Find the starting location in LSTRES array:
C
      IF(IRES.EQ.0) THEN
        WRITE(LFNPR,900) IRES
        CALL NBHALT('TOPSTR: Index IRES out of bounds.')
      ELSE IF(IRES.GT.MAXRES) THEN
        WRITE(LFNPR,910) MAXRES
        CALL NBHALT('TOPSTR: Too many resonance structures.')
      ELSE
        IF(ILOC.GT.0.AND.LSTRES(ILOC).NE.-1) THEN
   40     ILOC = ILOC - 1
          IF(ILOC.GT.0.AND.LSTRES(ILOC).NE.-1) GOTO 40
        END IF
        IDXRES(IRES,IREF) = ILOC + 1
      END IF
C
C  Loop over atoms:
C
      DO 70 IA = 1,NATOMS
C
C  Store the number of lone pairs as first entry:
C
        ILOC = ILOC + 1
        LSTRES(ILOC) = ITOPO(IA,IA)
C
C  Save place for bonds:
C
        ILOC  = ILOC + 1
        ILOCB = ILOC
C
C  Store the bonds:
C
        NBONDS = 0
        DO 60 IB = IA+1,NATOMS
          NB = ITOPO(IA,IB)
   50     CONTINUE
          IF(NB.GT.0) THEN
            ILOC = ILOC + 1
            LSTRES(ILOC) = IB
            NBONDS = NBONDS + 1
            NB = NB - 1
            GOTO 50
          END IF
   60   CONTINUE
        LSTRES(ILOCB) = NBONDS
   70 CONTINUE
      ILOC = ILOC + 1
      LSTRES(ILOC) = -1
C
C  Did the LSTRES array overflow?
C
      IF(ILOC.GT.LEN) THEN
        WRITE(LFNPR,920) IRES,LEN
        CALL NBHALT('TOPSTR: Fatal memory overflow.')
      END IF
      RETURN
C
  900 FORMAT(/1X,'FATAL ERROR:  Index IRES (=',I3,') is out-of-bounds.')
  910 FORMAT(/1X,'FATAL ERROR:  Too many resonance structures, MAXRES ',
     + '= ',I5,'.')
  920 FORMAT(/1X,'Fatal memory overflow for resonance structure ',I4,
     +  ' (exceeds LEN =',I7,')',/1X,'See comments about scratch ',
     +  'vector partitioning in SR NBODRV.')
      END
C***********************************************************************
      SUBROUTINE TOPZER(NATOMS)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200)
      COMMON/NBSTR/ITOPO(MAXATM,MAXATM)
C
C  Zero the elements of the ITOPO matrix:
C
      DO 20 J = 1,NATOMS
        DO 10 I = 1,NATOMS
          ITOPO(I,J) = 0
   10   CONTINUE
   20 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE WGTOUT(WGT,NRES,MAXRES,IREF,NREF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION WGT(MAXRES,NREF),NRES(NREF)
      DIMENSION W(5),IR(5),IC(5)
C
      SAVE IBLNK,IPLUS
      DATA IBLNK,IPLUS/1H ,1H+/
C
C  Print the optimized weights for the reference structure IREF:
C
      WRITE(LFNPR,900) IREF
C
      IFLG = 0
      IRES = 0
   10 DO 20 I = 1,5
        IC(I) = IBLNK
   20 CONTINUE
      ICNT = 0
   30 IRES = IRES + 1
        ICNT = ICNT + 1
        IF(IFLG.EQ.1) IC(ICNT) = IPLUS
        IFLG = 1
        W(ICNT) = WGT(IRES,IREF)
        IR(ICNT) = IRES
      IF(IRES.LT.NRES(IREF).AND.ICNT.LT.5) GOTO 30
      WRITE(LFNPR,910) (IC(I),W(I),IR(I),I=1,ICNT)
      IF(IRES.LT.NRES(IREF)) GOTO 10
      RETURN
C
  900 FORMAT(/1X,'Resonance weights for reference structure',I3,':')
  910 FORMAT(1X,5(A1,F8.5,'(',I4,')',:,1X))
      END
C***********************************************************************
      SUBROUTINE WGTPR(WGTM,MAP,NWGT)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION WGTM(NWGT),MAP(NWGT)
      DIMENSION W(4),IR(4),IC(4)
C
      SAVE IBLNK,IPLUS
      DATA IBLNK,IPLUS/1H ,1H+/
C
C  Print the multi-ref optimized resonance weights:
C
      WRITE(LFNPR,900)
C
      IFLG = 0
      IRES = 0
   10 DO 20 I = 1,4
        IC(I) = IBLNK
   20 CONTINUE
      ICNT = 0
   30 IRES = IRES + 1
        ICNT = ICNT + 1
        IF(IFLG.EQ.1) IC(ICNT) = IPLUS
        IFLG = 1
        W(ICNT) = WGTM(IRES)
        IR(ICNT) = MAP(IRES)
      IF(IRES.LT.NWGT.AND.ICNT.LT.4) GOTO 30
      WRITE(LFNPR,910) (IC(I),W(I),IR(I),I=1,ICNT)
      IF(IRES.LT.NWGT) GOTO 10
      RETURN
C
  900 FORMAT(/1X,'Resonance weights by IDXRES index:')
  910 FORMAT(1X,4(A1,F8.5,'[',I6,']',:,1X))
      END
C***********************************************************************
      SUBROUTINE WNORM(WGT,IRESET,IW,IP,SKIP)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL SKIP
C
      DIMENSION WGT(IP),IRESET(IP)
C
      SAVE ZERO,ONE
      DATA ZERO,ONE/0.0D0,1.0D0/
C
      SKIP = .TRUE.
      IF(IRESET(IW).LT.0) RETURN
      WGT(IP) = ONE
      DO 10 IREF = 1,IP-1
        IF(IRESET(IREF).GE.0) THEN
          WGT(IP) = WGT(IP) - WGT(IREF)
        END IF
   10 CONTINUE
      IF(WGT(IP).LT.ZERO) RETURN
      SKIP = .FALSE.
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NBBP MODULE ROUTINES:
C
C      SUBROUTINE NBBP(T,TMO,BBP,SCR,LFN)
C      SUBROUTINE FNDFIL(LFNFIL,ERROR)
C      SUBROUTINE NBPARS(STR,NV,IVAL)
C      FUNCTION NBINTV(STR)
C      FUNCTION LENNB(STRING)
C
C***********************************************************************
      SUBROUTINE NBBP(T,TMO,BBP,SCR,LFN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER LINE*80,CRS*1,CTU*1
C
C  LFNBBP=file of additional NBBP indices (written by FNDFIL from
C         file-input bracket), 4 indices (free format) per line
C
      PARAMETER(MAXBAS = 2000)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBL(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
C
      DIMENSION T(NDIM,NDIM),TMO(NDIM,NDIM),BBP(NDIM,NDIM),
     +          SCR(NDIM*(NDIM+5)),EORB(MAXBAS)
      DIMENSION IHVAL(4)
C
      SAVE ZERO,ONE,TWO,KLEW,ISTAR
      DATA ZERO/0.0D0/ONE/1.0D0/TWO/2.0D0/
      DATA KLEW/3HLEW/
      DATA ISTAR/1H*/IRYD/2HRY/
C
C  Is this an open-shell calculation?  Set OCC = 1 or 2
C
      IF(OPEN) THEN
        OCC=ONE
      ELSE
        OCC=TWO
      END IF
C
C  Fetch the Fock matrix in the AO basis, store it in T:
C
      CALL FEFAO(T,IWFOCK)
C
C  Fetch the AO to MO transformation matrix (using BBP as storage):
C
      CALL FEAOMO(BBP,IT)
      IF(IT.EQ.0) RETURN
C
C  Diagonalize T (using matrix BBP of eigenvectors) to find orbital energies,
C
      CALL SIMTR1(T,BBP,EORB,NDIM,NBAS,NBAS,NNAO)
      DO I=1,NNAO
        EORB(I)=T(I,I)
      ENDDO
C
C  Find the NHO-MO transformation matrix:
C
C  First, find the AO-NHO transformation, store in T:
      CALL FETNAO(T)
      CALL FETNHO(BBP)
      CALL MATML1(T,BBP,SCR,NDIM,NDIM,NDIM,NBAS,NNAO,NNAO)
C  Form the NHO-MO transformation, store in TMO
      CALL FRMTMO(T,TMO,BBP,SCR,3,0)
C
C  Evaluate the bond-bond polarizability [Coulson & Longuet-Higgins,
C  Proc. Roy. Soc. A191, 39 (1947); A192, 16 (1947)] for bonded hybrids.
C  Columns labelled as NBOs to specify the hybrids NBO1=IAB=(IH1,IH2)
C  and NBO2=ICD=(IH3,IH4) of the Pi(IH1,IH2;IH3,IH4) array.  This array
C  includes some "atom-bond" polarizabilities (if one NBO is a lone pair)
C  or "atom-atom" polarizabilities (if both NBOs are lone pairs).
C
C  Find the number of occupied (NOCC=NLEW) and minimal-basis (NVAL) orbitals
      NOCC=0
      NVAL=0
      DO I=1,NNAO
        IF(LABEL(I,2).NE.ISTAR.AND.LABEL(I,3).NE.IRYD) NOCC=NOCC+1
        IF(LABEL(I,3).NE.IRYD) NVAL=NVAL+1
      ENDDO
      NLEW=NOCC
C  Search first the bonded hybrid pairs
      IH1=0
      IH2=0
      DO IBO1=1,NLEW
        IB1=IBXM(IBO1)
        NCTR=1
        IF(LABEL(IB1,5).NE.0) NCTR=2
        IF(LABEL(IB1,6).NE.0) CALL NBHALT('Erroneous NBBP label.')
        IF(NCTR.EQ.1) THEN
          IH1=IH2+1
          IH2=IH1
        ELSE
          IH1=IH2+1
          IH2=IH1+1
        ENDIF
        IH3=0
        IH4=0
        DO IBO2=1,NLEW
          IB2=IBXM(IBO2)
          NCTR=1
          IF(LABEL(IB2,5).NE.0) NCTR=2
          IF(LABEL(IB2,6).NE.0) CALL NBHALT('Erroneous NBBP label.')
          IF(NCTR.EQ.1) THEN
            IH3=IH4+1
            IH4=IH3
          ELSE
            IH3=IH4+1
            IH4=IH3+1
          ENDIF
          PIVAL=ZERO
C
C  Main sum over occupied (J) and virtual (K) MOs
C
          DO J=1,NOCC
            DO K=NOCC+1,NNAO
              V12=TMO(IH1,J)*TMO(IH2,K)+TMO(IH1,K)*TMO(IH2,J)
              V34=TMO(IH3,J)*TMO(IH4,K)+TMO(IH3,K)*TMO(IH4,J)
              PIVAL=PIVAL+V12*V34/(EORB(K)-EORB(J))
            ENDDO
          ENDDO
          PIVAL=OCC*PIVAL
          BBP(IBO1,IBO2)=PIVAL
        ENDDO
      ENDDO
C
C  Print out the matrix of bond-bond (NBO-NBO) bond polarizabilities:
C
      LINE='NBBP: NBO bond-bond polarizability matrix'
      IF(OPEN) THEN
        NL=LENNB(LINE)
        IF(ALPHA) THEN
          LINE=LINE(1:NL)//' (alpha spin)'
        ELSE
          LINE=LINE(1:NL)//' (beta spin)'
        END IF
      END IF
      IFLG=KLEW
      INDEX=4
      CALL AOUT(BBP,NDIM,NLEW,NLEW,LINE,INDEX,IFLG)
C
C  Compute user-requested nonbonded elements of Pi(IH1,IH2;IH3,IH4):
C
      IF(LFN.LE.0) RETURN
      WRITE(LFNPR,1000)
      WRITE(LFNPR,1010)
      WRITE(LFNPR,1020)
      WRITE(LFNPR,1030)
      WRITE(LFNPR,1040)
      IERR=0
C 4-index (bond-bond polarizabilities, etc.) r, s, t, u:
   10 CONTINUE
      READ(LFN,1100,END=2000) LINE
      CALL NBPARS(LINE,NV,IHVAL)
      IF(NV.EQ.4) THEN
        IH1=IHVAL(1)
        IH2=IHVAL(2)
        IH3=IHVAL(3)
        IH4=IHVAL(4)
      ELSE
        IERR=1
        WRITE(LFNPR,1200) NV
        GOTO 10
      END IF
      DO IV=1,NV
        IF(IHVAL(IV).LE.0.OR.IHVAL(IV).GT.NNAO) THEN
          IERR=1
          NL=LENNB(LINE)
          WRITE(LFNPR,1300) IHVAL(IV),(LINE(K:K),K=1,NL)
          GOTO 10
        END IF
      ENDDO
C Evaluate the bond-bond polarizability for requested indices
      PIVAL=ZERO
      DO J=1,NOCC
        DO K=NOCC+1,NNAO
          V12=TMO(IH1,J)*TMO(IH2,K)+TMO(IH1,K)*TMO(IH2,J)
          V34=TMO(IH3,J)*TMO(IH4,K)+TMO(IH3,K)*TMO(IH4,J)
          PIVAL=PIVAL+V12*V34/(EORB(K)-EORB(J))
        ENDDO
      ENDDO
      PIVAL=OCC*PIVAL
C Evaluate the signs of Fock matrix elements F(r,s), F(t,u):
      FRS=ZERO
      FTU=ZERO
      DO K=1,NNAO
        FRS=FRS+TMO(IH1,K)*TMO(IH2,K)*EORB(K)
        FTU=FTU+TMO(IH3,K)*TMO(IH4,K)*EORB(K)
      ENDDO
      IF(FRS.GE.ZERO) THEN
        CRS='+'
      ELSE
        CRS='-'
      END IF
      IF(FTU.GE.ZERO) THEN
        CTU='+'
      ELSE
        CTU='-'
      END IF
      IF((FRS*FTU).GE.ZERO) THEN
        CPIVAL=PIVAL
      ELSE
        CPIVAL=-PIVAL
      END IF
      WRITE(LFNPR,1400) IH1,IH2,IH3,IH4,PIVAL,CRS,CTU,CPIVAL
      GOTO 10
C ERROR CONDITION READING LFN?
 2000 CONTINUE
      IF(IERR.EQ.1) THEN
        WRITE(LFNPR,2010)
      END IF
      CLOSE(LFN)
      RETURN
C
 1000 FORMAT(/,1X,'User-requested Pi(r,s;t,u) NBBP values:',/)
 1010 FORMAT(1X,'    NHO indices    uncorrected      sign',
     + '        "corrected"')
 1020 FORMAT(1X,'   -------------   -----------  -------------',
     + '   -----------')
 1030 FORMAT(1X,'   r   s   t   u   Pi(r,s;t,u)  F(r,s) F(t,u)',
     + '   Pi(r,s;t,u)')
 1040 FORMAT(1X,'   -------------   -----------  ------ ------',
     + '   -----------')
 1100 FORMAT(A80)
 1200 FORMAT(1X,'Wrong number (',I2,') of polarizability indices',
     + ' (should be 4).')
 1300 FORMAT(1X,'Illegal value',i4,' in polarizability indices: ',
     +   50A1)
 1400 FORMAT(1X,4I4,F13.6,6X,A1,6X,A1,2X,F13.6)
 2010 FORMAT(1X,'WARNING:  Erroneous <input> for NBBP indices.')
      END
C***********************************************************************
      SUBROUTINE FNDFIL(LFNFIL,ERROR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ERROR
      DIMENSION LINE(256)
C
      COMMON/NBCRD1/ICD(256),LOOK(256),LENGTH,IPT,LFN,NEXP
      COMMON/NBCRD2/POINT,LEND,NEXT,EXP
      LOGICAL POINT,LEND,NEXT,EXP
C
      SAVE NBLA,NLBR,NSL,NRBR,NEXC,LFNTRY,IFIL
      DATA NBLA,NLBR,NSL,NRBR,NEXC/1H ,1H<,1H/,1H>,1H!/
      DATA IFIL,LFNTRY/0,61/
C
C READ THE LINES OF <line 1/line 2/...> FILE INPUT AND WRITE TO FILE LFNFIL.
C LFNFIL IS INCREMENTED EACH TIME A NEW FILE-FIELD IS FOUND, STARTING FROM
C 61, 62,...  IF NO <...> FILE-FIELD IS FOUND, RETURN WITH LFNFIL=0.
C ERROR IS .TRUE. IF INCORRECT SYNTAX IS ENCOUNTERED AFTER STARTING '<' (IN
C THIS CASE, THE POINTER IS PROBABLY LEFT AT THE END OF THE FILE).
C
      IF(LOOK(1).NE.NLBR) THEN
        LFNFIL=0
        ERROR=.FALSE.
        RETURN
      END IF
C '<' WAS FOUND; INITIALIZE FOR READING ENTRIES
      NEWLINE=1
      LFNTRY=LFNTRY+1
      IFIL=IFIL+1
      IF(IFIL.GT.20) IFIL=1
      LFNFIL=LFNTRY
      OPEN(LFNFIL,STATUS='SCRATCH')
      IPR=0
C ANYTHING STILL IN LOOK?
      IF(LENGTH.GT.1) THEN
        DO ITRY=2,LENGTH
          IPR=IPR+1
          LINE(IPR)=LOOK(ITRY)
        END DO
        IPR=IPR+1
C ADD A BLANK SPACE FOR SAFETY
        LINE(IPR)=NBLA
      END IF
C CONTINUE READING THE LINE
      NCOL=IPT-1
   10 IF(NCOL.GE.256) THEN
        CALL NXTCRD
        IF(LEND) THEN
          LFNFIL=0
          ERROR=.FALSE.
          RETURN
        END IF
        NCOL=0
      END IF
   15 NCOL=NCOL+1
      ICARD=ICD(NCOL)
      IF(ICARD.EQ.NEXC) THEN
        CALL NXTCRD
        IF(LEND) THEN
          LFNFIL=0
          ERROR=.FALSE.
          RETURN
        ELSE
          NCOL=0
          GOTO 15
        END IF
C IS THIS THE END '>'?
      ELSE IF(ICARD.EQ.NRBR) THEN
        IF(IPR.GT.0) THEN
          WRITE(LFNFIL,1000) (LINE(K),K=1,IPR)
        END IF
        REWIND LFNFIL
        IPR=0
        ERROR=.FALSE.
        IPT=NCOL+1
        CALL FNDFLD
        RETURN
C IS THIS A LINE SEPARATOR '/'?
      ELSE IF(ICARD.EQ.NSL) THEN
        IF(IPR.GT.0) THEN
          WRITE(LFNFIL,1000) (LINE(K),K=1,IPR)
        END IF
        IPR=0
        NEWLINE=1
C OR JUST ANOTHER CHARACTER?
      ELSE
        IF(ICARD.NE.NBLA.OR.NEWLINE.EQ.0) THEN
          IPR=IPR+1
          LINE(IPR)=ICARD
          NEWLINE=0
        END IF
      END IF
C ON TO NEXT CHARACTER...
      GOTO 10
 1000 FORMAT(1X,256A1)
      END
C***********************************************************************
      SUBROUTINE NBPARS(STR,NV,IVAL)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER STR*(*)
      DIMENSION IVAL(*)
C Parse string 'STR' into NV integer values IVAL(I)
      NV=0
      NL=LENNB(STR)
      K2=0
      K1=0
   10 CONTINUE
      NV=NV+1
      K1=K2+1
   20 CONTINUE
      IF(STR(K1:K1).EQ.' '.OR.STR(K1:K1).EQ.',') THEN
        K1=K1+1
        IF(K1.LE.NL) THEN
          GOTO 20
        ELSE
          NV=NV-1
          RETURN
        END IF
      END IF
      K2=K1
   30 CONTINUE
      IF(STR(K2:K2).EQ.' '.OR.STR(K2:K2).EQ.',') THEN
        K2=K2-1
        IVAL(NV)=NBINTV(STR(K1:K2))
        GOTO 10
      ELSE
        IF(K2.LT.NL) THEN
          K2=K2+1
          GOTO 30
        ELSE
          IVAL(NV)=NBINTV(STR(K1:K2))
        END IF
      END IF
      RETURN
      END
C***********************************************************************
      FUNCTION NBINTV(STR)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      CHARACTER STR*(*)
      DIMENSION IPOW(7)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      SAVE IPOW
      DATA IPOW/1,10,100,1000,10000,100000,1000000/
C
C CONVERTS A STRING VARIABLE TO INTEGER VALUE
C
      NL=LENNB(STR)
      ITMP=0
      IMAG=0
   10 IF(NL.LE.0) THEN
        NBINTV=ITMP
        RETURN
      ELSE
        IC=ICHAR(STR(NL:NL))
        I=IC-48
        IF(I.LT.0.OR.I.GT.9) THEN
          WRITE(LFNPR,1000) STR(NL:NL),IC,(STR(K:K),K=1,LENNB(STR))
 1000     FORMAT(1X,'*** Illegal text character (',A1,', IC=',I3,
     +     ') in numerical value: ',80A1)
          NBINTV=0
          RETURN
        END IF
        IMAG=IMAG+1
        ITMP=ITMP+I*IPOW(IMAG)
        NL=NL-1
      END IF
      GOTO 10
      END
C***********************************************************************
      FUNCTION LENNB(STRING)
C***********************************************************************
      IMPLICIT REAL (A-H,O-Z)
      CHARACTER*(*) STRING
C
C  FIND LENGTH OF NON-BLANK LEFT SUB-STRING OF STRING:
C
      NL = LEN(STRING)
      IF(NL.LE.0) THEN
        LENNB = 0
        RETURN
      ELSE IF (NL.GT.511) THEN
        NL = 511
      END IF
      DO 10 I = NL,1,-1
        IF(ICHAR(STRING(I:I)).NE.32) THEN
          LENNB = I
          RETURN
        END IF
   10 CONTINUE
      LENNB = 0
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NBSTER (STERIC) MODULE ROUTINES:
C
C      SUBROUTINE NBSTER(T,S,FNLMO,SCR,MOLNBO,LFNSTR)
C      SUBROUTINE WORTUN(FO,S,W,FP,IERR)
C
C***********************************************************************
      SUBROUTINE NBSTER(T,S,FNLMO,SCR,MOLNBO,LFNSTR)
C***********************************************************************
C  4-Jan-94 JKB  Fixed up user-defined indices section
C 15-Oct-93 JKB  New subroutine to calculate steric interaction energies
C 21-Nov-00 FAW  Correct AUKCAL, EVKCAL conversion factors
C 21-Feb-01 FAW  Replace NBOs by NLMOs throughout
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z),INTEGER(I-N)
      EXTERNAL UNPACK
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      LOGICAL DETAIL
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),MOLLST(MAXBAS),IATHY(MAXBAS,3)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBMOL/NMOLEC,MOLAT(MAXATM),MOLEC(MAXATM,MAXATM),
     +              NMOLA,MOLATA(MAXATM),MOLECA(MAXATM,MAXATM)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
      COMMON/ESTER/ESTREA,ESTR2A
      COMMON/IOCWAR/IOCWAR
      DIMENSION T(NDIM,NDIM),S(NDIM,NDIM),
     +         FNLMO(NDIM,NDIM),SCR(NDIM)
      DIMENSION PNBEN(MAXBAS)
      DIMENSION FNLMO2(2,2),FPB2(2,2),SB2(2,2),WN2(2)
      DIMENSION MOLNBO(2,NBAS,NMOLEC)
      DIMENSION INAM(3),JNAM(3),ICH(3,2),JCH(3,2),ISP(3),JSP(3)
C
      SAVE LBD,L3C,LBLNK1,LBLNK2,LHYP,LSTAR
      SAVE AUKCAL,EVKCAL,ZERO,ONE,TWO,TEN
      DATA LBD/2HBD/,L3C/2H3C/,LBLNK1/1H /,LBLNK2/2H  /,LHYP/1H-/
      DATA AUKCAL/627.5093314D0/,EVKCAL/23.061D0/
      DATA ZERO,ONE,TWO,TEN,HALF/0.0D0,1.0D0,2.0D0,1.0D1,0.5D0/
      DATA LSTAR/1H*/
C
C  PERFORM NBO STERIC ANALYSIS:
C
C  STTHR1 IS THE INTRAMOLECULAR THRESHOLD FOR PRINTING STERIC INTERACTION
C  ENERGIES (DEFAULT: 0.5 KCAL/MOL FOR CLOSED SHELL, 0.25 KCAL/MOL FOR OPEN
C  SHELL).  STTHR2=STTHR1/10 IS THE INTERMOLECULAR THRESHOLD.
C
      WRITE(LFNPR,2000)
      IF(IOCWAR.GT.0) WRITE(LFNPR,2005)
      NELE = 0
      STTHR1 = ABS(STTHR)
      IF(ISPIN.NE.0.AND.STTHR.LT.ZERO) STTHR1 = STTHR1/TWO
      STTHR2 = ABS(STTHR)/TEN
      IF(ISPIN.NE.0.AND.STTHR.LT.ZERO) STTHR2 = STTHR2/TWO
      OCC=TWO
      IF(ISPIN.NE.0) OCC=ONE
      DETAIL=.FALSE.
      IF(IWDETL.NE.0) DETAIL=.TRUE.
C
C  Form the AO-PNBO matrix (store in T, using S, SCR scratch)
C
      CALL FEPNAO(T)
      CALL FETNAB(S)
      CALL MATML1(T,S,SCR,NDIM,NDIM,NDIM,NBAS,NNAO,NNAO)
C
C  Fetch the NBO Fock matrix (store in FNLMO)
C
      CALL FEFNBO(FNLMO)
      CALL UNPACK(FNLMO,NDIM,NNAO)
C  Fetch the NBO-NLMO transform in S
      CALL FETLMO(S)
C  Transform Fock matrix to NLMO basis
      CALL SIMTRS(FNLMO,S,SCR,NDIM,NNAO)
C  Form the AO-PNLMO transform
      CALL MATML1(T,S,SCR,NDIM,NDIM,NDIM,NBAS,NNAO,NNAO)
C  T is now the AOPNLMO matrix; check normalization
      CALL FESRAW(S)
      CALL NORMLZ(T,S,NDIM,NBAS,NNAO)
C  Form the PNLMO Fock matrix (store in S; SCR scratch)
      CALL FEFAO(S,IWFOCK)
      CALL SIMTR1(S,T,SCR,NDIM,NBAS,NBAS,NNAO)
C  Store the PNLMO orbital energies
      DO 10 KK=1,NNAO
   10 PNBEN(KK)=S(KK,KK)
C  Form the PNLMO overlap matrix (store in S; SCR scratch)
      CALL FESRAW(S)
      CALL SIMTR1(S,T,SCR,NDIM,NBAS,NBAS,NNAO)
C  Currently,
C     S = PNLMO overlap matrix
C     T = AO-PNLMO transformation
C FNLMO = NLMO Fock matrix
C
C  MAKE UP LIST MOLNBO(1,IBAS,IMOL) OF CORE/LP/BOND NBOS IN MOLEC. UNIT IMOL
C               MOLNBO(2,IBAS,IMOL) OF RYDBERG/ANTIBOND NBOS IN MOLEC. IMOL
C
      DO 40 IMOL = 1,NMOLEC
        DO 20 IBAS = 1,NNAO
          MOLNBO(1,IBAS,IMOL) = 0
          MOLNBO(2,IBAS,IMOL) = 0
   20   CONTINUE
        NOCC = 0
        NSTAR = 0
        DO 30 IBAS = 1,NNAO
          IF(IMOL.EQ.NBOUNI(IBAS)) THEN
            IB = IBXM(IBAS)
            IF(LABEL(IB,2).NE.LSTAR) THEN
              NOCC = NOCC + 1
              MOLNBO(1,NOCC,IMOL) = IBAS
            ELSE
              NSTAR = NSTAR + 1
              MOLNBO(2,NSTAR,IMOL) = IBAS
            ENDIF
          ENDIF
   30   CONTINUE
   40 CONTINUE
C
C  DETERMINE THE CONVERSION FROM INPUT ENERGY UNITS TO KCAL:
C
      IF(MUNIT.EQ.0) THEN
        CONV = AUKCAL
      ELSE IF(MUNIT.EQ.1) THEN
        CONV = EVKCAL
      ELSE
        CONV = ONE
      ENDIF
C
C  If DETAIL, Write orbital energies in a.u.
C
      IF(DETAIL) THEN
        WRITE(LFNPR,2010)
        WRITE(LFNPR,2030) (PNBEN(K),K=1,NNAO)
        WRITE(LFNPR,2020)
        WRITE(LFNPR,2030) (FNLMO(K,K),K=1,NNAO)
      ENDIF
C
C  Write orbital energy differences/shifts by unit and calculate
C    sums by unit and overall
C
      WRITE(LFNPR,2040)
      SUM1A=ZERO
      NOCC=0
      DO 70 IMOL = 1,NMOLEC
        WRITE(LFNPR,2050) IMOL
        SUM1U=ZERO
        DO 60 IOCC = 1,NNAO
          IBAS = MOLNBO(1,IOCC,IMOL)
          IF(IBAS.EQ.0) GO TO 60
          NOCC=NOCC+1
          IB = IBXM(IBAS)
          LBL0 = LABEL(IB,1)
          IF(LBL0.EQ.LBD) THEN
            NCTR = 2
          ELSE IF(LBL0.EQ.L3C) THEN
            NCTR = 3
          ELSE
            NCTR = 1
          ENDIF
          DO 50 I = 1,3
            IA = LABEL(IB,I+3)
            CALL CONVRT2N(IA,ICH(I,1),ICH(I,2))
            INAM(I) = LBLNK2
            IF(IA.GT.0) INAM(I) = NAMEAT(IATNO(IA))
            ISP(I) = LHYP
            IF(I.GE.NCTR) ISP(I) = LBLNK1
   50     CONTINUE
          CENRG=HALF*OCC*CONV*(FNLMO(IBAS,IBAS)-PNBEN(IBAS))
          WRITE(LFNPR,2060) IBAS,(LABEL(IB,K),K=1,3),
     *           (INAM(K),ICH(K,1),ICH(K,2),ISP(K),K=1,2),
     *            INAM(3),ICH(3,1),ICH(3,2),CENRG
          SUM1U=SUM1U+CENRG
   60   CONTINUE
        WRITE(LFNPR,2070) IMOL,SUM1U
        SUM1A=SUM1A+SUM1U
   70 CONTINUE
      WRITE(LFNPR,2080) SUM1A
C
C  Breakdown steric effects pairwise and print table
C
        WRITE(LFNPR,2090) STTHR1
        IF(NMOLEC.GT.1) WRITE(LFNPR,2100) STTHR2
        WRITE(LFNPR,2110)
        IERR3=0
        SUM2A=ZERO
        DO 150 IMOL = 1,NMOLEC
          DO 150 JMOL = IMOL,NMOLEC
            SUM2U=ZERO
            IF(IMOL.EQ.JMOL) THEN
              WRITE(LFNPR,2120) IMOL
              STTHRS = STTHR1
            ELSE
              WRITE(LFNPR,2130) IMOL,JMOL
              STTHRS = STTHR2
            ENDIF
C
C  Loop through basis functions in order to temporarily assign unit 1
C     as the one containing the basis function of lower number
C
            DO 140 MOCC = 1,NNAO
              IMATCH=0
              DO 80 KOCC = 1,NNAO
                KBAS = MOLNBO(1,KOCC,IMOL)
                IF(KBAS.EQ.MOCC) THEN
                  IOCC=KOCC
                  KMOL=IMOL
                  LMOL=JMOL
                  IMATCH=1
                ENDIF
   80         CONTINUE
              DO 90 LOCC = 1,NNAO
                LBAS = MOLNBO(1,LOCC,JMOL)
                IF(LBAS.EQ.MOCC) THEN
                  IOCC=LOCC
                  KMOL=JMOL
                  LMOL=IMOL
                  IMATCH=1
                ENDIF
   90         CONTINUE
              IF(IMATCH.EQ.0) GO TO 140
              DO 130 JOCC = 1,NNAO
                IBAS = MOLNBO(1,IOCC,KMOL)
                JBAS = MOLNBO(1,JOCC,LMOL)
                IF(JBAS.EQ.0) GO TO 130
                IF(IBAS.GE.JBAS) GO TO 130
C
C  Check if orbitals on different centers
C
                IX=IBXM(IBAS)
                JX=IBXM(JBAS)
                KFLG=1
                DO 100 ICAT=4,6
                  DO 100 JCAT=4,6
                    IF(LABEL(IX,ICAT).EQ.LABEL(JX,JCAT).AND.
     *                 LABEL(IX,ICAT).NE.0) KFLG=0
  100           CONTINUE
                IF(KFLG.EQ.0) GO TO 130
C
C  2x2 de-orthogonalization of NBOs
C
                FNLMO2(1,1) = FNLMO(IBAS,IBAS)
                FNLMO2(1,2) = FNLMO(IBAS,JBAS)
                FNLMO2(2,1) = FNLMO(JBAS,IBAS)
                FNLMO2(2,2) = FNLMO(JBAS,JBAS)
                SB2(1,1)   = S(IBAS,IBAS)
                SB2(1,2)   = S(IBAS,JBAS)
                SB2(2,1)   = S(JBAS,IBAS)
                SB2(2,2)   = S(JBAS,JBAS)
                WN2(1)     = OCC
                WN2(2)     = OCC
                CALL WORTUN(FNLMO2,SB2,WN2,FPB2,IERR3)
                IF(IERR3.EQ.1) THEN
                  WRITE(LFNPR,2140) IBAS,JBAS
                  DER=ZERO
                ELSE
                  DER=HALF*OCC*CONV*((FNLMO2(1,1)-FPB2(1,1))+
     *                               (FNLMO2(2,2)-FPB2(2,2)))
                  SUM2U=SUM2U+DER
                ENDIF
C
C  If steric E above threshold, fetch the two NBO labels and print
C
                IF(ABS(DER).LT.STTHRS) GO TO 130
                NELE = NELE + 1
                IB = IBXM(IBAS)
                LBL0 = LABEL(IB,1)
                IF(LBL0.EQ.LBD) THEN
                  NCTR = 2
                ELSE IF(LBL0.EQ.L3C) THEN
                  NCTR = 3
                ELSE
                  NCTR = 1
                ENDIF
                DO 110 I = 1,3
                  IA = LABEL(IB,I+3)
                  CALL CONVRT2N(IA,ICH(I,1),ICH(I,2))
                  INAM(I) = LBLNK2
                  IF(IA.GT.0) INAM(I) = NAMEAT(IATNO(IA))
                  ISP(I) = LHYP
                  IF(I.GE.NCTR) ISP(I) = LBLNK1
  110           CONTINUE
                JB = IBXM(JBAS)
                LBL0 = LABEL(JB,1)
                IF(LBL0.EQ.LBD) THEN
                  NCTR = 2
                ELSE IF(LBL0.EQ.L3C) THEN
                  NCTR = 3
                ELSE
                  NCTR = 1
                ENDIF
                DO 120 J = 1,3
                  JA = LABEL(JB,J+3)
                  CALL CONVRT2N(JA,JCH(J,1),JCH(J,2))
                  JNAM(J) = LBLNK2
                  IF(JA.GT.0) JNAM(J) = NAMEAT(IATNO(JA))
                  JSP(J) = LHYP
                  IF(J.GE.NCTR) JSP(J) = LBLNK1
  120           CONTINUE
                WRITE(LFNPR,2150) IBAS,(LABEL(IB,K),K=1,3),
     *             (INAM(K),ICH(K,1),ICH(K,2),ISP(K),K=1,2),
     *              INAM(3),ICH(3,1),ICH(3,2),
     *                           JBAS,(LABEL(JB,K),K=1,3),
     *             (JNAM(K),JCH(K,1),JCH(K,2),JSP(K),K=1,2),
     *              JNAM(3),JCH(3,1),JCH(3,2),
     *                           S(IBAS,JBAS),DER
  130         CONTINUE
  140       CONTINUE
            IF(NELE.EQ.0) WRITE(LFNPR,2160)
            IF(IMOL.EQ.JMOL) WRITE(LFNPR,2170) IMOL,SUM2U
            IF(IMOL.NE.JMOL) WRITE(LFNPR,2180) IMOL,JMOL,SUM2U
            SUM2A=SUM2A+SUM2U
  150   CONTINUE
        WRITE(LFNPR,2190) SUM2A
C
C  If open shell alpha spin, save steric energies;
C              if beta spin, compute overall steric energy
C
      IF(ISPIN.EQ.2) THEN
        ESTREA = SUM1A
        ESTR2A = SUM2A
      END IF
      IF(ISPIN.EQ.-2) THEN
        ESTRET = ESTREA + SUM1A
        ESTR2T = ESTR2A + SUM2A
        WRITE(LFNPR,2220) ESTRET,ESTR2T
      END IF
C
C  Read in user-requested indices from LFNSTR, calculate steric energy
C
      IUSR = 1
C 16-Oct-03 FAW Fix for Compaq compiler
      IF(LFNSTR.EQ.0) GOTO 190
  160 READ(LFNSTR,*,ERR=190,END=190) IUBAS,JUBAS
      IF(IUSR.EQ.1) WRITE(LFNPR,2200)
      IF(IUSR.EQ.1) WRITE(LFNPR,2110)
      IUSR = 0
      IF(IUBAS.LE.0.OR.JUBAS.LE.0.OR.IUBAS.GT.NOCC.OR.
     *  JUBAS.GT.NOCC) THEN
        WRITE(LFNPR,2210) IUBAS,JUBAS
      ELSE
        FNLMO2(1,1) = FNLMO(IUBAS,IUBAS)
        FNLMO2(1,2) = FNLMO(IUBAS,JUBAS)
        FNLMO2(2,1) = FNLMO(JUBAS,IUBAS)
        FNLMO2(2,2) = FNLMO(JUBAS,JUBAS)
        SB2(1,1)   = S(IUBAS,IUBAS)
        SB2(1,2)   = S(IUBAS,JUBAS)
        SB2(2,1)   = S(JUBAS,IUBAS)
        SB2(2,2)   = S(JUBAS,JUBAS)
        WN2(1)     = OCC
        WN2(2)     = OCC
        CALL WORTUN(FNLMO2,SB2,WN2,FPB2,IERR3)
        IF(IERR3.EQ.1) THEN
          WRITE(LFNPR,2140) IUBAS,JUBAS
          DER=ZERO
        ELSE
          DER=HALF*OCC*CONV*((FNLMO2(1,1)-FPB2(1,1))+
     *                       (FNLMO2(2,2)-FPB2(2,2)))
        ENDIF
        IB = IBXM(IUBAS)
        LBL0 = LABEL(IB,1)
        IF(LBL0.EQ.LBD) THEN
          NCTR = 2
        ELSE IF(LBL0.EQ.L3C) THEN
          NCTR = 3
        ELSE
          NCTR = 1
        ENDIF
        DO 170 I = 1,3
          IA = LABEL(IB,I+3)
          CALL CONVRT2N(IA,ICH(I,1),ICH(I,2))
          INAM(I) = LBLNK2
          IF(IA.GT.0) INAM(I) = NAMEAT(IATNO(IA))
          ISP(I) = LHYP
          IF(I.GE.NCTR) ISP(I) = LBLNK1
  170   CONTINUE
        JB = IBXM(JUBAS)
        LBL0 = LABEL(JB,1)
        IF(LBL0.EQ.LBD) THEN
          NCTR = 2
        ELSE IF(LBL0.EQ.L3C) THEN
          NCTR = 3
        ELSE
          NCTR = 1
        ENDIF
        DO 180 J = 1,3
          JA = LABEL(JB,J+3)
          CALL CONVRT2N(JA,JCH(J,1),JCH(J,2))
          JNAM(J) = LBLNK2
          IF(JA.GT.0) JNAM(J) = NAMEAT(IATNO(JA))
          JSP(J) = LHYP
          IF(J.GE.NCTR) JSP(J) = LBLNK1
  180   CONTINUE
        WRITE(LFNPR,2150) IUBAS,(LABEL(IB,K),K=1,3),
     *      (INAM(K),ICH(K,1),ICH(K,2),ISP(K),K=1,2),
     *       INAM(3),ICH(3,1),ICH(3,2),
     *                   JUBAS,(LABEL(JB,K),K=1,3),
     *      (JNAM(K),JCH(K,1),JCH(K,2),JSP(K),K=1,2),
     *       JNAM(3),JCH(3,1),JCH(3,2),
     *                   S(IUBAS,JUBAS),DER
      ENDIF
      GOTO 160
  190 CONTINUE
      WRITE(LFNPR,2230)
      RETURN
C
 2000 FORMAT(//1X,'NBO/NLMO STERIC ANALYSIS')
 2005 FORMAT(/1X,'WARNING:  Non-physical occupancies were found in ',
     * 'NAO or NBO search.',/11X,'NBSTER will attempt to continue.')
 2010 FORMAT(/2X,'PNLMO Fock matrix, diagonal elements (a.u.)',/)
 2020 FORMAT(/2X,'NLMO Fock matrix, diagonal elements (a.u.)',/)
 2030 FORMAT(4X,6F12.5)
 2040 FORMAT(/2X,'Occupied NLMO contributions dE(i) (kcal/mol) to ',
     *     'total steric exchange energy')
 2050 FORMAT(/2X,'NLMOs (i) in unit',I3,14X,'dE(i)')
 2060 FORMAT(1X,I4,'. ',A2,A1,'(',I2,')',A2,3A1,A2,3A1,A2,2A1,
     *          5X,F8.2)
 2070 FORMAT(/2X,'Steric exchange energy, unit',I2,':',F8.2,
     *     ' kcal/mol')
 2080 FORMAT(2X,48('-'),/2X,'Total steric exchange energy:',
     *     F10.2,' kcal/mol',/2X,48('-'))
 2090 FORMAT(/,2X,'Pairwise steric exchange energies dE(i,j)',
     * ' (kcal/mol) and associated',/,2X,'pre-NLMO overlaps S(i,j)',
     * ' for disjoint (no common atoms) interactions',
     * /,2X,'between NLMOs i,j:',
     * //,10X,'Threshold for printing:  ',F7.2,' kcal/mol')
 2100 FORMAT(9X,'(Intermolecular threshold:',F7.2,' kcal/mol)')
 2110 FORMAT(57X,'PNLMO',4X,'dE(i,j)',/,9X,'NLMO (i)',19X,'NLMO (j)',
     *     12X,'S(i,j)',4X,'kcal/mol',/,2X,72('='))
 2120 FORMAT(/2X,'within unit',I3)
 2130 FORMAT(/2X,'between units',I3,' and',I3)
 2140 FORMAT(5X,'on interaction of NLMOs',2I4,';',/5X,
     *    'Setting steric exchange energy to zero',/)
 2150 FORMAT(1X,I4,'. ',A2,A1,'(',I2,')',A2,3A1,A2,3A1,A2,2A1,
     *   I4,'. ',A2,A1,'(',I2,')',A2,3A1,A2,3A1,A2,2A1,
     *       F8.4,F9.2,F7.2)
 2160 FORMAT(2X,'      None above threshold')
 2170 FORMAT(/38X,'sum within unit',I2,':',7X,F9.2,F7.2)
 2180 FORMAT(/38X,'sum between units',I2,' and',I2,':',F8.2,F7.2)
 2190 FORMAT(/2X,70('-'),/,2X,'Total disjoint NLMO steric exchange ',
     *    'energy from pairwise sum:',F9.2,F7.2,/,2X,70('-'),/)
 2200 FORMAT(/2X,'User-requested steric interactions:')
 2210 FORMAT(/2X,'Requested indices',I4,' and',I4,
     *    ' invalid for occupied orbitals')
 2220 FORMAT(2X,'Overall (spin-averaged) total steric exchange energy:',
     *    F8.2,' kcal/mol',/2X,'Overall disjoint NLMO steric ',
     *    'exchange energy:',5X,F11.2,' kcal/mol',/,2X,70('-'),/)
 2230 FORMAT(2X,72('-'),/)
      END
C***********************************************************************
      SUBROUTINE WORTUN(FO,S,W,FP,IERR)
C***********************************************************************
C 15-Oct-93  JKB  New subroutine, calculates 2x2 "pseudo-pre-orthogonal"
C                 Fock matrix by doing occupancy-weighted
C                 "deorthogonalization" of two natural orbitals
C 31-Oct-00  FAW  Added declaration of SCR(2,2)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION W(2),W2(2,2),S(2,2),WSW(2,2),FP(2,2),FO(2,2)
      DIMENSION OM(2,2),OMT(2,2),TMP(2,2),SCR(2,2)
      DIMENSION EVAL(2),EVALH(2,2),EVEC(2,2),EVECT(2,2)
C
      SAVE ZERO,XLIM
      DATA ZERO/0.0D00/XLIM/1.0D-13/
C
      IERR=0
      DO 10 I=1,2
        DO 10 J=1,2
          EVALH(I,J)=ZERO
          IF(I.EQ.J) THEN
            W2(I,J)=W(I)
          ELSE
            W2(I,J)=ZERO
          ENDIF
   10 CONTINUE
      CALL COPY(W2,TMP,2,2,2)
      CALL MATMLT(TMP,S,SCR,2,2)
      CALL MATMLT(TMP,W2,SCR,2,2)
      CALL COPY(TMP,WSW,2,2,2)
      CALL NBJACOBI(2,TMP,EVAL,EVEC,2,2,0)
      DO 20 J=1,2
        IF(EVAL(J).LT.XLIM) THEN
          IERR=1
          EVALH(J,J)=ZERO
        ELSE
          EVALH(J,J)=SQRT(EVAL(J))
        ENDIF
   20 CONTINUE
      IF(IERR.EQ.1) THEN
        WRITE(LFNPR,110) EVAL(1),EVAL(2)
        RETURN
      ENDIF
      CALL COPY(EVEC,EVECT,2,2,2)
      CALL NBTRSP(EVECT,2,2)
      IF(ABS(W2(1,1)).LT.XLIM.OR.ABS(W2(2,2)).LT.XLIM) THEN
        IERR=1
        WRITE(LFNPR,120) W2(1,1),W2(2,2)
        RETURN
      ELSE
C
C  Multiply ULUt where U = WSW eigvecs, L = WSW eigvals(+1/2)
C
        CALL COPY(EVEC,TMP,2,2,2)
        CALL MATMLT(TMP,EVALH,SCR,2,2)
        CALL MATMLT(TMP,EVECT,SCR,2,2)
C
C  Multiply this by W(-1) on right to get Ow(-1)
C
        W2(1,1)=1.0D0/W2(1,1)
        W2(2,2)=1.0D0/W2(2,2)
        CALL MATMLT(TMP,W2,SCR,2,2)
        CALL COPY(TMP,OMT,2,2,2)
        CALL COPY(OMT,OM,2,2,2)
        CALL NBTRSP(OM,2,2)
C
C  Finally, take Ow(-1)*FNLMO*Ow(-1)t'
C
        CALL COPY(OM,TMP,2,2,2)
        CALL MATMLT(TMP,FO,SCR,2,2)
        CALL MATMLT(TMP,OMT,SCR,2,2)
        CALL COPY(TMP,FP,2,2,2)
      ENDIF
  110 FORMAT(/3X,'Negative eigenvalue found =',2F9.5,
     *      /3X,'Cannot form Ow(-1) matrix; aborting SR WORTUN')
  120 FORMAT(/3X,'Small overlap matrix element found =',2D13.5,
     *      /3X,'Cannot form inverse W matrix; aborting SR WORTUN')
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  NEDA ROUTINES:
C
C      SUBROUTINE NBOEDA(T,DM,SCR1,SCR2)
C      SUBROUTINE EDAOUT(ETOT,E,SNRG,FENRG,FPNRG,FSNRG,EDEF,ECP,EFLD,
C     +                  EDFF,DEF,LISTA)
C      SUBROUTINE SVNEWT(T)
C      SUBROUTINE FENEWT(T)
C      SUBROUTINE SVTEDA(T)
C      SUBROUTINE FETEDA(T)
C      SUBROUTINE SVEIG(U)
C      SUBROUTINE FEEIG(U)
C      SUBROUTINE SVDDEF(SCR,NFRG)
C      SUBROUTINE FEDDEF(SCR,NFRG)
C      SUBROUTINE SVDFLD(SCR,NFRG)
C      SUBROUTINE FEDFLD(SCR,NFRG)
C      SUBROUTINE SVDCP(SCR,NFRG)
C      SUBROUTINE FEDCP(SCR,NFRG)
C      SUBROUTINE SVVNUC(SCR,NFRG)
C      SUBROUTINE FEVNUC(SCR,NFRG)
C
C***********************************************************************
      SUBROUTINE NBOEDA(T,DM,SCR1,SCR2)
C***********************************************************************
C  7-Jul-94  EDG  Write all density matrices to record 86,87 of the DAF
C  3-Mar-94  EDG  Add UHF capabilities
C  2-Mar-94  EDG  Number of occupied orbs determined from NFE rather
C                   than LABEL to be consistent with NBODEL
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      EXTERNAL PACK,UNPACK
C
      PARAMETER(MAXBAS = 2000)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       IATNO(MAXBAS),IBXM(MAXBAS),LOCC(2*MAXBAS),ISCR(2*MAXBAS)
      COMMON/NBEDA/KFRG,NFRG,IFRG(MAXBAS),NFE(MAXBAS),NFEB(MAXBAS)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION T(NDIM,NDIM),DM(*),SCR1(NDIM,NDIM),SCR2(*)
C
      SAVE ZERO,ONE,TWO,BIG
      DATA ZERO,ONE,TWO,BIG/0.0D0,1.0D0,2.0D0,1.0D6/
C
C  Prepare the density matrix corresponding to Psi(def) for current fragment:
C
C  Increment the current fragment:
C
      IF(ISPIN.NE.-2) KFRG = KFRG + 1
C
C  Restore COMMON/NBBAS/ and retrieve the NBO-MO transformation from
C  the NBO DAF:
C
      CALL FENBLB(T)
      CALL FENEWT(T)
C
C  Prepare AO-MO transformation:
C
      CALL FETNBO(SCR1)
      CALL MATML1(SCR1,T,SCR2,NDIM,NDIM,NDIM,NBAS,NNAO,NNAO)
C
C  Count the number of occupied orbitals on this fragment:
C
      IF(ISPIN.EQ.0)  NOCC = NFE(KFRG) / 2
      IF(ISPIN.EQ.2)  NOCC = NFE(KFRG)
      IF(ISPIN.EQ.-2) NOCC = NFEB(KFRG)
C
C  Classify the MOs by fragment:
C
      DO 30 IBAS = 1,NNAO
        TMAX = ZERO
        ISCR(IBAS) = 0
        DO 20 JBAS = 1,NNAO
          IF(ABS(T(JBAS,IBAS)).GT.TMAX) THEN
            TMAX = ABS(T(JBAS,IBAS))
            ISCR(IBAS) = JBAS
          END IF
   20   CONTINUE
        IF(ISCR(IBAS).EQ.0) CALL NBHALT('Error in SR NBOEDA.')
        ISCR(IBAS) = IFRG(NBOUNI(ISCR(IBAS)))
   30 CONTINUE
C
C  Get the eigenvalues of the truncated Fock matrix:
C
      CALL FEEIG(SCR2)
C
C  Prepare the list (LOCC) of occupied orbitals on this fragment:
C
      DO 50 IOCC = 1,NOCC
        TMP = BIG
        DO 40 IBAS = 1,NNAO
          IF(ISCR(IBAS).EQ.KFRG) THEN
            IF(SCR2(IBAS).LT.TMP) THEN
              TMP = SCR2(IBAS)
              LOCC(IOCC) = IBAS
            END IF
          END IF
   40   CONTINUE
        IF(LOCC(IOCC).EQ.0) CALL NBHALT('Error in SR NBOEDA.')
        SCR2(LOCC(IOCC)) = BIG
   50 CONTINUE
C
C  Construct the NBO density matrix:
C
      ETA = TWO
      IF(ISPIN.NE.0) ETA = ONE
      IJ = 0
      DO 80 I = 1,NNAO
        DO 70 J = 1,I
          SUM = ZERO
          DO 60 K = 1,NOCC
            SUM = SUM + T(I,LOCC(K)) * T(J,LOCC(K))
   60     CONTINUE
          IJ = IJ + 1
          DM(IJ) = SUM * ETA
   70   CONTINUE
   80 CONTINUE
C
C  Sort the vectors in SCR1 (AO-MO transform) so that the occupied orbitals
C  are in the first few columns and store on DAF:
C
      DO 100 IOCC = 1,NOCC
        I = LOCC(IOCC)
        DO 90 J = 1,NNAO
          T(J,IOCC) = SCR1(J,I)
   90   CONTINUE
  100 CONTINUE
      I = 0
      DO 140 IUNOCC = NOCC+1,NNAO
  110   I = I + 1
        DO 120 IOCC = 1,NOCC
          IF(I.EQ.LOCC(IOCC)) GOTO 110
  120   CONTINUE
        IF(I.GT.NNAO) CALL NBHALT('Error in NBOEDA.')
        DO 130 J = 1,NNAO
          T(J,IUNOCC) = SCR1(J,I)
  130   CONTINUE
  140 CONTINUE
      CALL SVTEDA(T)
C
C  Backtransform DM to the AO basis and save on the NBO DAF:
C
      CALL FETNBO(T)
      CALL NBTRSP(T,NDIM,NBAS)
      CALL UNPACK(DM,NDIM,NNAO)
      CALL SIMTR1(DM,T,SCR1,NDIM,NNAO,NNAO,NBAS)
      CALL PACK(DM,NDIM,NBAS)
      CALL SVNEWD(DM)
C
C  Also store this density on records 86,87 of the DAF:
C
      L2 = NBAS * (NBAS+1) / 2
      IOFF = (KFRG - 1) * L2 + 1
      IF(KFRG.EQ.1) THEN
        CALL COPY(DM,SCR2(IOFF),L2,L2,1)
        CALL SVDDEF(SCR2,NFRG)
      ELSE
        CALL FEDDEF(SCR2,NFRG)
        CALL COPY(DM,SCR2(IOFF),L2,L2,1)
        CALL SVDDEF(SCR2,NFRG)
      END IF
C
C  Set IATNO elements negative for atomic centers not on current fragment:
C
      DO 150 I = 1,NATOMS
        IF(ISPIN.NE.-2) IATNO(I) = ABS(IATNO(I))
        ISCR(I)  = -1
  150 CONTINUE
      DO 170 IBAS = 1,NNAO
        IF(IFRG(NBOUNI(IBAS)).EQ.KFRG) THEN
          IB = IBXM(IBAS)
          DO 160 J = 4,6
            IAT = LABEL(IB,J)
            IF(IAT.GT.0) ISCR(IAT) = 1
  160     CONTINUE
        END IF
  170 CONTINUE
C
      IF(ISPIN.NE.-2) THEN
        DO 180 IAT = 1,NATOMS
          IATNO(IAT) = IATNO(IAT) * ISCR(IAT)
  180   CONTINUE
      ELSE
        DO 190 IAT = 1,NATOMS
          J = IATNO(IAT) * ISCR(IAT)
          IF(J.LT.0) THEN
            WRITE(LFNPR,900)
            CALL NBHALT('Different fragments for alpha/beta spins.')
          END IF
  190   CONTINUE
      END IF
      RETURN
C
  900 FORMAT(/1X,'NEDA is not applicable.')
      END
C***********************************************************************
      SUBROUTINE EDAOUT(ETOT,E,SNRG,FENRG,FPNRG,FSNRG,EDEF,ECP,EFLD,
     +                  EDFF,DEF,LISTA)
C***********************************************************************
C 12-Jun-05  EDG  Added DFT/NEDA functionality
C 18-May-01  EDG  Warn when not using extended NEDA
C 21-Nov-00  FAW  Correct AUKCAL, EVKCAL conversion factors
C 30-May-98  EDG  Add electric field NEDA treatment
C  2-Jun-95  EDG  Add evaluation of self energies
C 17-Jan-95  EDG  NEDA extended to calculate ES, POL, and EX
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL EXTND,FIELD,DFT
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBATOM/IAT(MAXATM),INO(MAXATM),NORBS(MAXATM),ILL(MAXATM),
     +       ILU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       IATNO(MAXBAS),IBXM(MAXBAS),LOCC(2*MAXBAS),ISCR2(2*MAXBAS)
      COMMON/NBEDA/KFRG,NFRG,IFRG(MAXBAS),NFE(MAXBAS),NFEB(MAXBAS)
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
C
      DIMENSION E(4+7*NFRG),SNRG(NFRG),FENRG(NFRG),FPNRG(NFRG),
     +          FSNRG(NFRG),EDEF(NFRG),ECP(NFRG),EFLD(NFRG),EDFF(NFRG),
     +          DEF(NFRG),LISTA(NATOMS,2),ISTR(80),IK(5)
      DIMENSION DFF(MAXATM)
C
      SAVE ZERO,AUKCAL
      SAVE IBLNK,LL,LR,LP,LM
C
      DATA ZERO,AUKCAL/0.0D0,627.5093314D0/
      DATA IBLNK,LL,LR,LP,LM/1H ,1H(,1H),1H+,1H-/
C
C  Reset IATNO array:
C
      DO 10 I = 1,NATOMS
        IATNO(I) = ABS(IATNO(I))
   10 CONTINUE
C
C  Get energies from NBO DAF:
C
      NFILE = 70
      LEN = 5 + 8 * NFRG
      CALL NBREAD(E,LEN,NFILE)
C
C  Evaluate NEDA components (extended treatment-1996):
C
      ELX  = E(1)
      DE   = ZERO
      SE   = ZERO
      FE   = ZERO
      FP   = ZERO
      FS   = ZERO
      FD   = ZERO
      EINT = ETOT
      FIELD = .FALSE.
      DO 20 I = 1,NFRG
        SNRG(I)  = E(5       +I)
        EDEF(I)  = E(5+  NFRG+I)
        ECP(I)   = E(5+2*NFRG+I)
        EFLD(I)  = E(5+3*NFRG+I)
        EDFF(I)  = E(5+4*NFRG+I)
        FENRG(I) = E(5+5*NFRG+I)
        FPNRG(I) = E(5+6*NFRG+I)
        FSNRG(I) = E(5+7*NFRG+I)
        IF(ECP(I).NE.EFLD(I)) FIELD = .TRUE.
        ELX      = ELX - EDEF(I) + EFLD(I) - ECP(I)
        DEF(I)   = (EDEF(I) - ECP(I)) * AUKCAL
        DFF(I)   = (EDFF(I) - ECP(I)) * AUKCAL
        DE       = DE + DEF(I)
        FD       = FD + DFF(I)
        SNRG(I)  =  SNRG(I) * AUKCAL
        FENRG(I) =  FENRG(I) * AUKCAL
        FPNRG(I) =  FPNRG(I) * AUKCAL
        FSNRG(I) =  FSNRG(I) * AUKCAL
        SE       = SE + SNRG(I)
        FE       = FE + FENRG(I)
        FP       = FP + FPNRG(I)
        FS       = FS + FSNRG(I)
        EINT     = EINT - ECP(I)
   20 CONTINUE
      EINT = EINT * AUKCAL
      ELX  = ELX  * AUKCAL
      ES   = E(2) * AUKCAL
      POL  = E(3) * AUKCAL
      EX   = E(4) * AUKCAL
      DFT  = E(5).NE.ZERO
C     POL  = ELX - ES - EX
      EL   = ES + POL + SE
      CT   = (ETOT - E(1)) * AUKCAL
      CORE = EINT - EL - CT
      EC   = CORE - (EX+DE-SE)
      XC   = EX + EC
C
C  Original NEDA treatment (1994):
C
      EXTND = .TRUE.
      IF(ES.EQ.ZERO.AND.EX.EQ.ZERO) THEN
        EXTND = .FALSE.
        ES = ELX
        EINT = ES + CT + DE
      END IF
C
C  If the last fragment is being processed, write out details of
C  the analysis:
C
      WRITE(LFNPR,900)
C
C  Form chemical formula for complex:
C
      NAT = 0
      DO 40 I = 1,NATOMS
        DO 30 J = 1,NAT
          IF(IATNO(I).EQ.LISTA(J,1)) THEN
            LISTA(J,2) = LISTA(J,2) + 1
            GOTO 40
          END IF
   30   CONTINUE
        NAT = NAT + 1
        LISTA(NAT,1) = IATNO(I)
        LISTA(NAT,2) = 1
   40 CONTINUE
      CALL CHEM(NAT,NATOMS,LISTA,NL,ISTR)
      NL = NL - 2
      DO 50 I = 1,NL
        ISTR(I) = ISTR(I+1)
   50 CONTINUE
      IC = 0
      DO 60 I = 1,NATOMS
        IC = IC + IZNUC(I)
   60 CONTINUE
      DO 70 I = 1,NFRG
        IC = IC - NFE(I)
        IF(UHF) IC = IC - NFEB(I)
   70 CONTINUE
      IF(IC.NE.0) THEN
        NL = NL + 1
        ISTR(NL) = LL
        IF(IC.GT.0) THEN
          NL = NL + 1
          ISTR(NL) = LP
        ELSE
          NL = NL + 1
          ISTR(NL) = LM
          IC = ABS(IC)
        END IF
        IF(IC.GT.1) THEN
          CALL IDIGIT(IC,IK,ND,5)
          DO 80 I = 1,ND
            NL = NL + 1
            ISTR(NL) = IK(I)
   80     CONTINUE
        END IF
        NL = NL + 1
        ISTR(NL) = LR
      END IF
      NL = NL + 1
      ISTR(NL) = IBLNK
      DO 90 I = NL+1,80
        ISTR(I) = IBLNK
   90 CONTINUE
C
C  Write info about complex:
C
      WRITE(LFNPR,910) (ISTR(I),I=1,13),ETOT,E(1),CT
      IF(EXTND) THEN
        WRITE(LFNPR,920) ES
        WRITE(LFNPR,930) POL
        IF(.NOT.DFT) WRITE(LFNPR,940) EX
        IF(     DFT) WRITE(LFNPR,945) XC
      ELSE
        WRITE(LFNPR,920) ES
      END IF
C
C  Loop over fragments, forming chemical formulas for each and writing
C  energies:
C
      DO 170 I = 1,NFRG
        IC = -NFE(I)
        IF(UHF) IC = IC - NFEB(I)
        NAT = 0
        DO 130 J = 1,NATOMS
          DO 120 IBAS = 1,NNAO
            IB = IBXM(IBAS)
            DO 110 K = 4,6
              IF(LABEL(IB,K).NE.0) THEN
                IF(LABEL(IB,K).EQ.J) THEN
                  IF(IFRG(NBOUNI(IBAS)).EQ.I) THEN
                    IC = IC + IZNUC(J)
                    DO 100 L = 1,NAT
                      IF(LISTA(L,1).EQ.IATNO(J)) THEN
                        LISTA(L,2) = LISTA(L,2) + 1
                        GOTO 130
                      END IF
  100               CONTINUE
                    NAT = NAT + 1
                    LISTA(NAT,1) = IATNO(J)
                    LISTA(NAT,2) = 1
                  END IF
                  GOTO 130
                END IF
              END IF
  110       CONTINUE
  120     CONTINUE
  130   CONTINUE
        CALL CHEM(NAT,NATOMS,LISTA,NL,ISTR)
        NL = NL - 2
        DO 140 J = 1,NL
          ISTR(J) = ISTR(J+1)
  140   CONTINUE
        DO 150 J = NL+1,80
          ISTR(J) = IBLNK
  150   CONTINUE
        IF(IC.NE.0) THEN
          NL = NL + 1
          ISTR(NL) = LL
          IF(IC.GT.0) THEN
            NL = NL + 1
            ISTR(NL) = LP
          ELSE
            NL = NL + 1
            ISTR(NL) = LM
            IC = ABS(IC)
          END IF
          IF(IC.GT.1) THEN
            CALL IDIGIT(IC,IK,ND,5)
            DO 160 J = 1,ND
              NL = NL + 1
              ISTR(NL) = IK(J)
  160       CONTINUE
          END IF
          NL = NL + 1
          ISTR(NL) = LR
        END IF
        IF(EXTND) THEN
          IF(.NOT.FIELD) THEN
            WRITE(LFNPR,950) I,(ISTR(J),J=1,10),EDEF(I),ECP(I),
     +                       DEF(I),SNRG(I)
          ELSE
            WRITE(LFNPR,960) I,(ISTR(J),J=1,10),EDEF(I),ECP(I),
     +                       DEF(I),SNRG(I),EFLD(I)
          END IF
        ELSE
          WRITE(LFNPR,970) I,(ISTR(J),J=1,10),EDEF(I),ECP(I),
     +                     DEF(I)
        END IF
  170 CONTINUE
      WRITE(LFNPR,980) EINT
      IF(EXTND.AND..NOT.DFT) WRITE(LFNPR,990) EL,CT,CORE,EINT
      IF(EXTND.AND.    DFT) WRITE(LFNPR,995) EL,CT,CORE,EINT
      IF(.NOT.EXTND) WRITE(LFNPR,1000)
C
C  Check interaction energy:
C
      IF(EXTND.AND..NOT.DFT) THEN
        IF(ABS(EC).GT.1.0D-4) WRITE(LFNPR,1010) EINT-EC,EINT,EC
      END IF
C
C  Write extended field output:
C
      IF(FIELD.AND.EXTND) THEN
        WRITE(LFNPR,1020)
C
C  Form chemical formula for complex:
C
        NAT = 0
        DO 240 I = 1,NATOMS
          DO 230 J = 1,NAT
            IF(IATNO(I).EQ.LISTA(J,1)) THEN
              LISTA(J,2) = LISTA(J,2) + 1
              GOTO 240
            END IF
  230     CONTINUE
          NAT = NAT + 1
          LISTA(NAT,1) = IATNO(I)
          LISTA(NAT,2) = 1
  240   CONTINUE
        CALL CHEM(NAT,NATOMS,LISTA,NL,ISTR)
        NL = NL - 2
        DO 250 I = 1,NL
          ISTR(I) = ISTR(I+1)
  250   CONTINUE
        IC = 0
        DO 260 I = 1,NATOMS
          IC = IC + IZNUC(I)
  260   CONTINUE
        DO 270 I = 1,NFRG
          IC = IC - NFE(I)
          IF(UHF) IC = IC - NFEB(I)
  270   CONTINUE
        IF(IC.NE.0) THEN
          NL = NL + 1
          ISTR(NL) = LL
          IF(IC.GT.0) THEN
            NL = NL + 1
            ISTR(NL) = LP
          ELSE
            NL = NL + 1
            ISTR(NL) = LM
            IC = ABS(IC)
          END IF
          IF(IC.GT.1) THEN
            CALL IDIGIT(IC,IK,ND,5)
            DO 280 I = 1,ND
              NL = NL + 1
              ISTR(NL) = IK(I)
  280       CONTINUE
          END IF
          NL = NL + 1
          ISTR(NL) = LR
        END IF
        NL = NL + 1
        ISTR(NL) = IBLNK
        DO 290 I = NL+1,80
          ISTR(I) = IBLNK
  290   CONTINUE
C
C  Write info about complex:
C
        WRITE(LFNPR,910) (ISTR(I),I=1,13),ETOT,E(1),CT
        WRITE(LFNPR,920) ES-FE
        WRITE(LFNPR,930) POL-FP
        WRITE(LFNPR,940) EX
C
C  Loop over fragments, forming chemical formulas for each and writing
C  energies:
C
        DO 370 I = 1,NFRG
          IC = -NFE(I)
          IF(UHF) IC = IC - NFEB(I)
          NAT = 0
          DO 330 J = 1,NATOMS
            DO 320 IBAS = 1,NNAO
              IB = IBXM(IBAS)
              DO 310 K = 4,6
                IF(LABEL(IB,K).NE.0) THEN
                  IF(LABEL(IB,K).EQ.J) THEN
                    IF(IFRG(NBOUNI(IBAS)).EQ.I) THEN
                      IC = IC + IZNUC(J)
                      DO 300 L = 1,NAT
                        IF(LISTA(L,1).EQ.IATNO(J)) THEN
                          LISTA(L,2) = LISTA(L,2) + 1
                          GOTO 330
                        END IF
  300                 CONTINUE
                      NAT = NAT + 1
                      LISTA(NAT,1) = IATNO(J)
                      LISTA(NAT,2) = 1
                    END IF
                    GOTO 330
                  END IF
                END IF
  310         CONTINUE
  320       CONTINUE
  330     CONTINUE
          CALL CHEM(NAT,NATOMS,LISTA,NL,ISTR)
          NL = NL - 2
          DO 340 J = 1,NL
            ISTR(J) = ISTR(J+1)
  340     CONTINUE
          DO 350 J = NL+1,80
            ISTR(J) = IBLNK
  350     CONTINUE
          IF(IC.NE.0) THEN
            NL = NL + 1
            ISTR(NL) = LL
            IF(IC.GT.0) THEN
              NL = NL + 1
              ISTR(NL) = LP
            ELSE
              NL = NL + 1
              ISTR(NL) = LM
              IC = ABS(IC)
            END IF
            IF(IC.GT.1) THEN
              CALL IDIGIT(IC,IK,ND,5)
              DO 360 J = 1,ND
                NL = NL + 1
                ISTR(NL) = IK(J)
  360         CONTINUE
            END IF
            NL = NL + 1
            ISTR(NL) = LR
          END IF
          WRITE(LFNPR,1030) I,(ISTR(J),J=1,10),EDEF(I),ECP(I),
     +                 DEF(I)-DFF(I),SNRG(I)-FSNRG(I),EDFF(I),EFLD(I)
  370   CONTINUE
        EL   = (ES-FE) + (POL-FP) + (SE-FS)
        CORE = EX + (DE-FD) - (SE-FS)
        EINT = EL + CT + CORE
        WRITE(LFNPR,980) EINT
        WRITE(LFNPR,990) EL,CT,CORE,EINT
C
C  Check interaction energy:
C
        CHK = ETOT
        DO 380 I = 1,NFRG
          CHK = CHK - EFLD(I)
  380   CONTINUE
        CHK = CHK * AUKCAL
        IF(ABS(CHK-EINT).GT.1D-4) WRITE(LFNPR,1010) EINT,CHK,EINT-CHK
      END IF
C
C  and reset KFRG and NFRG
C
      KFRG = 0
      NFRG = 0
      RETURN
C
  900 FORMAT(//1X,'Natural Energy Decomposition Analysis (Summary):',
     + //1X,'                                                        ',
     + '   Component',/1X,'                   Energy(wfn)        Ener',
     + 'gy(wfn)          (kcal/mol)',/1X,'---------------------------',
     + '---------------------------------------------------')
  910 FORMAT(1X,13A1,F14.7,'(scf)',F14.7,'(loc)      CT =',F9.2)
  920 FORMAT(1X,'                                                    ',
     + '     ES =',F9.2)
  930 FORMAT(1X,'                                                    ',
     + '    POL =',F9.2)
  940 FORMAT(1X,'                                                    ',
     + '     EX =',F9.2)
  945 FORMAT(1X,'                                                    ',
     + '     XC =',F9.2)
  950 FORMAT(I2,'. ',10A1,F14.7,'(def)',F14.7,'(cp)  DEF(SE) =',F9.2,
     + '(',F7.2,')')
  960 FORMAT(I2,'. ',10A1,F14.7,'(def)',F14.7,'(cp)  DEF(SE) =',F9.2,
     + '(',F7.2,')',/33X,F14.7,'(fld)')
  970 FORMAT(I2,'. ',10A1,F14.7,'(def)',F14.7,'(cp)      DEF =',F9.2)
  980 FORMAT(1X,'                                                    ',
     + '         ---------',/1X,'                                    ',
     + '                      E =',F9.2)
  990 FORMAT(//1X,'Electrical (ES+POL+SE) :',F11.2,/1X,'  Charge Tran',
     + 'sfer (CT) :',F11.2,/1X,'      Core (EX+DEF-SE) :',F11.2,/1X,
     + '                       ------------',/1X,' Total Interaction ',
     + '(E) :',F11.2)
  995 FORMAT(//1X,'Electrical (ES+POL+SE) :',F11.2,/1X,'  Charge Tran',
     + 'sfer (CT) :',F11.2,/1X,'      Core (XC+DEF-SE) :',F11.2,/1X,
     + '                       ------------',/1X,' Total Interaction ',
     + '(E) :',F11.2)
 1000 FORMAT(/1X,'Extended NEDA requires NOPK=1 in $INTGRL and NOSYM=',
     + '1 in $CONTRL.')
 1010 FORMAT(/1X,'Warning: Error evaluating E; NEDA components sum to',
     + F12.5,' kcal/mol',/40X,'but expected',F12.5,' kcal/mol.',/39X,
     + '(Difference =',F12.5,' kcal/mol)')
 1020 FORMAT(//1X,'Natural Energy Decomposition Analysis (fragment-fi',
     + 'eld contributions deleted):',//1X,'                          ',
     + '                                 Component',/1X,'            ',
     + '       Energy(wfn)        Energy(wfn)          (kcal/mol)',
     + /1X,'---------------------------------------------------------',
     + '---------------------')
 1030 FORMAT(I2,'. ',10A1,F14.7,'(def)',F14.7,'(cp)  DEF(SE) =',F9.2,
     + '(',F7.2,')',/14X,F14.7,'(dff)',F14.7,'(fld)')
      END
C***********************************************************************
      SUBROUTINE SVNEWT(T)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION T(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVNEWT:  SAVE NAO/NHO/NBO/NLMO TO MO TRANSFORMATION FOR NEDA/NBO DELETIONS
C
      NFILE = 71
      IF(BETA) NFILE = 72
      L3 = NDIM * NDIM
      CALL NBWRIT(T,L3,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FENEWT(T)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION T(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FENEWT:  FETCH NAO/NHO/NBO/NLMO TO MO TRANSFORMATION FOR NEDA/NBO DELETIONS
C
      NFILE = 71
      IF(BETA) NFILE = 72
      L3 = NDIM * NDIM
      CALL NBREAD(T,L3,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVTEDA(T)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION T(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVTEDA:  SAVE AO TO MO TRANSFORMATION FOR NEDA
C
      NFILE = 84
      IF(BETA) NFILE = 85
      L3 = NDIM * NDIM
      CALL NBWRIT(T,L3,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FETEDA(T)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION T(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FETEDA:  FETCH AO TO MO TRANSFORMATION FOR NEDA/NBO DELETIONS
C
      NFILE = 84
      IF(BETA) NFILE = 85
      L3 = NDIM * NDIM
      CALL NBREAD(T,L3,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVEIG(U)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION U(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVEIG:  SAVE EIGENVALUES OF THE TRUNCATED FOCK MATRIX
C
      NFILE = 82
      IF(BETA) NFILE = 83
      CALL NBWRIT(U,NDIM,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEEIG(U)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION U(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FEEIG:  FETCH EIGENVALUES OF THE TRUNCATED FOCK MATRIX
C
      NFILE = 82
      IF(BETA) NFILE = 83
      CALL NBREAD(U,NDIM,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVDDEF(SCR,NFRG)
C***********************************************************************
C  7-Jul-94  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVDDEF:  SAVE FRAGMENT DENSITY MATRICES (DEF) ON DAF
C
      NFILE = 86
      IF(BETA) NFILE = 87
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBWRIT(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEDDEF(SCR,NFRG)
C***********************************************************************
C  7-Jul-94  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FEDDEF:  FETCH FRAGMENT DENSITY MATRICES (DEF) FROM DAF
C
      NFILE = 86
      IF(BETA) NFILE = 87
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBREAD(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVDFLD(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVDFLD:  SAVE FRAGMENT DENSITY MATRICES (FLD) ON DAF
C
      NFILE = 88
      IF(BETA) NFILE = 89
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBWRIT(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEDFLD(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FEDFLD:  FETCH FRAGMENT DENSITY MATRICES (FLD) FROM DAF
C
      NFILE = 88
      IF(BETA) NFILE = 89
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBREAD(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVDCP(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  Moved to records 90,91
C 17-Jan-95  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVDCP:  SAVE FRAGMENT DENSITY MATRICES (CP) ON DAF
C
      NFILE = 90
      IF(BETA) NFILE = 91
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBWRIT(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEDCP(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  Moved to records 90,91
C 17-Jan-95  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      DIMENSION SCR(*)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FEDCP:  FETCH FRAGMENT DENSITY MATRICES (CP) FROM DAF
C
      NFILE = 90
      IF(BETA) NFILE = 91
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBREAD(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE SVVNUC(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  Moved to record 92
C 17-Jan-95  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SCR(*)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C SVVNUC:  SAVE FRAGMENT NUCLEAR ATTRACTION INTEGRALS ON DAF
C
      NFILE = 92
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBWRIT(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEVNUC(SCR,NFRG)
C***********************************************************************
C 18-Jun-97  EDG  Moved to record 92
C 17-Jan-95  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SCR(*)
C
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
C
C FEVNUC:  FETCH FRAGMENT NUCLEAR ATTRACTION INTEGRALS FROM DAF
C
      NFILE = 92
      LEN = NFRG * NDIM * (NDIM + 1) / 2
      CALL NBREAD(SCR,LEN,NFILE)
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  CHECKPOINTING ROUTINES:
C
C      SUBROUTINE SRTCHK(DM,T,S,SCR,IRNK)
C      SUBROUTINE RDPERM(T,S,IRNK,IARC,ITMP,IORG,LIST)
C      SUBROUTINE RDPRM2(TA,TB)
C      SUBROUTINE FILFLD(LFN,STR,LENGTH,IAMNEW,END)
C
C***********************************************************************
      SUBROUTINE SRTCHK(DM,T,S,SCR,IRNK)
C***********************************************************************
C 11-Feb-93  EDG  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXBAS = 2000)
C
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
C
      DIMENSION T(NDIM,NDIM),DM(NDIM,NDIM),S(NDIM,NDIM),SCR(NDIM)
      DIMENSION IRNK(NDIM),MAP(9)
C
      SAVE MAP
      DATA MAP/5,2,2,3,3,4,4,4,4/
C
      IF(ABS(JPRINT(63)).EQ.1) WRITE(LFNPR,901)
      IF(ABS(JPRINT(63)).EQ.2) WRITE(LFNPR,902)
      IF(ABS(JPRINT(63)).EQ.3) WRITE(LFNPR,903)
      IF(ABS(JPRINT(63)).EQ.4) WRITE(LFNPR,904)
      IF(ABS(JPRINT(63)).EQ.5) WRITE(LFNPR,905)
      IF(ABS(JPRINT(63)).EQ.6) WRITE(LFNPR,906)
      IF(ABS(JPRINT(63)).EQ.7) WRITE(LFNPR,907)
      IF(ABS(JPRINT(63)).EQ.8) WRITE(LFNPR,908)
      IF(ABS(JPRINT(63)).EQ.9) WRITE(LFNPR,909)
      N = NNAO
      IF(ABS(JPRINT(63)).EQ.1) N = NBAS
C
C  Retrieve the checkpoint matrices from the DAF:
C
      ISPIN = -2
   10 ISPIN = -ISPIN
      CALL FECHK(T,ISPIN,IT)
      IF(IT.NE.0) THEN
C
C  And get the density matrix and compute orbital occupancies:
C
        IF(OPEN) THEN
          ALPHA = ISPIN.GT.0
          BETA  = ISPIN.LT.0
        ELSE
          ALPHA = .FALSE.
          BETA  = .FALSE.
        END IF
        CALL FEDRAW(DM,SCR)
        IF(IWDM.NE.0) THEN
          CALL FESRAW(S)
          CALL SIMTRS(DM,S,SCR,NDIM,NBAS)
        END IF
        CALL SIMTR1(DM,T,SCR,NDIM,NBAS,NBAS,N)
        DO 20 I = 1,N
          SCR(I) = DM(I,I)
   20   CONTINUE
C
C  Sort orbitals by occupancy:
C
        CALL RANK(SCR,N,NDIM,IRNK)
C
C  Stuff sorted orbitals in S and write to DAF:
C
        DO 40 J = 1,N
          IR = IRNK(J)
          DO 30 I = 1,NBAS
            S(I,J) = T(I,IR)
   30     CONTINUE
   40   CONTINUE
        CALL SVCHK(S,ISPIN)
C
C  Write orbital sequence:
C
        IF(ALPHA) WRITE(LFNPR,930)
        IF(BETA)  WRITE(LFNPR,940)
        K = MAP(ABS(JPRINT(63)))
        WRITE(LFNPR,950) (J,(LBL(I,IRNK(J),K),I=1,10),J=1,N)
      END IF
      IF(ISPIN.GT.0) GOTO 10
      RETURN
C
  901 FORMAT(//1X,'Sorting PAOs by occupancy for checkpoint file:',/)
  902 FORMAT(//1X,'Sorting PNAOs by occupancy for checkpoint file:',/)
  903 FORMAT(//1X,'Sorting NAOs by occupancy for checkpoint file:',/)
  904 FORMAT(//1X,'Sorting PNHOs by occupancy for checkpoint file:',/)
  905 FORMAT(//1X,'Sorting NHOs by occupancy for checkpoint file:',/)
  906 FORMAT(//1X,'Sorting PNBOs by occupancy for checkpoint file:',/)
  907 FORMAT(//1X,'Sorting NBOs by occupancy for checkpoint file:',/)
  908 FORMAT(//1X,'Sorting PNLMOs by occupancy for checkpoint file:',/)
  909 FORMAT(//1X,'Sorting NLMOs by occupancy for checkpoint file:',/)
  930 FORMAT(1X,'Alpha orbitals:',/)
  940 FORMAT(1X,'Beta orbitals:',/)
  950 FORMAT(5(I4,'. ',10A1))
      END
C***********************************************************************
      SUBROUTINE RDPERM(T,S,IRNK,IARC,ITMP,IORG,LIST)
C***********************************************************************
C 28-Sep-93  FAW  Read permutation list from LFNPRM (rather than $PERM)
C 15-Feb-93  FAW  New subroutine to permute order of orbitals (uses LFNX = 55)
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C
      LOGICAL END
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      PARAMETER(NLB = 21)
      CHARACTER TLBL*21,STR*21,STR1*21,STR2*21,BL*21,ERRFRAG*12
C
C Read the LFNPRM permutation list and store the permuted order in IRNK
C
      COMMON/NBBAS/LABEL(MAXBAS,6),NAOCTR(MAXBAS),NAOL(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LTYP(MAXBAS),IATHY(MAXBAS,3)
      COMMON/NBLBL/NLEW,NVAL,LBLT(10,MAXBAS,5)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBPRMC/LFNPRM
C
      DIMENSION T(NDIM,NDIM),S(NDIM,NDIM)
      DIMENSION IRNK(NDIM),MAP(9),IVEC(10),ITMP(NDIM),IARC(NDIM)
      DIMENSION LIST(NDIM),TLBL(MAXBAS),IORG(NDIM)
      DIMENSION ISP(3),NAM(3),ICH(3,2)
C
      SAVE LLP,LBD,L3C,LCR,LRY
      SAVE LHYP,LBLNK,L2BLNK
      SAVE MAP
      DATA LLP,LBD,L3C,LCR,LRY/2HLP,2HBD,2H3C,2HCR,2HRY/
      DATA LHYP,LBLNK,L2BLNK/1H-,1H ,2H  /
      DATA MAP/5,2,2,3,3,4,4,4,4/
C
      BL='                     '
C Set initial values
      IFERMI=0
      IHFERMI=0
      IERR=0
      STR=BL
      STR1=BL
      STR2=BL
      IAMNEW=1
C Initialize IRNK, ITMP list
      DO 10 IB=1,NNAO
      IRNK(IB)=0
      IARC(IBXM(IB))=IB
      IORG(IB)=0
   10 ITMP(IB)=IBXM(IB)
C Look for proper spin heading?
      IF(OPEN) THEN
   20   CONTINUE
        CALL FILFLD(LFNPRM,STR1,LENG,IAMNEW,END)
        IF(LENG.EQ.0.OR.END) THEN
          IERR=2
          ERRFRAG='spin header'
          GOTO 130
        END IF
        IF(ALPHA) THEN
C Alpha spin?
          IF(STR1(1:LENG).NE.'ALPHA') GOTO 20
        ELSE IF(BETA) THEN
C Or Beta spin?
          IF(STR1(1:LENG).NE.'BETA') GOTO 20
        END IF
      END IF
C Grand loop over LIST entries
      IFLD=0
      IL=0
   30 CONTINUE
C Read the next entry: LIST(IL)
      NL=0
C Read another field of the entry
   40 CONTINUE
      NL1=NLB
      IFLD=IFLD+1
      CALL FILFLD(LFNPRM,STR1,NL1,IAMNEW,END)
      IF(END) GOTO 60
      IF(IERR.NE.0) THEN
        IERR=3
        ERRFRAG=STR1(1:12)
        GOTO 130
      END IF
      STR1=STR1(1:NL1)//BL(1:15-NL1)
C Is this the Fermi level ("|") marker?
      IF(NL1.EQ.1.AND.STR1(1:1).EQ.'|') THEN
        IHFERMI=1
        IFERMI=IL
        IFLD=IFLD-1
        GOTO 30
      END IF
C Test the first entry: Is it an integer?
      IF(IFLD.EQ.1) THEN
        IWLBL=0
        DO K=1,NL1
          IC=ICHAR(STR1(K:K))
          IF(IC.LT.48.OR.IC.GT.57) IWLBL=1
        END DO
        IF(IWLBL.EQ.1) THEN
C Is this a (P)NBO or (P)NLMO transformation?
          IF(ABS(JPRINT(63)).GE.6) THEN
C Prepare NBO/NLMO label LIST entries.  First, write labels to file LFNX...
            LFNX=55
            OPEN(LFNX,STATUS='SCRATCH')
            DO NBOND = 1,NNAO
              IB = ITMP(NBOND)
              LBL = LABEL(IB,1)
              IF(LBL.EQ.LLP.OR.LBL.EQ.LCR.OR.LBL.EQ.LRY) NCTR = 1
              IF(LBL.EQ.LBD) NCTR = 2
              IF(LBL.EQ.L3C) NCTR = 3
              DO I = 1,3
                IA = LABEL(IB,I+3)
                CALL CONVRT2N(IA,ICH(I,1),ICH(I,2))
                NAM(I) = L2BLNK
                IF(IA.GT.0) NAM(I) = NAMEAT(IATNO(IA))
                ISP(I) = LHYP
                IF(I.GE.NCTR) ISP(I) = LBLNK
              END DO
              IF(NCTR.EQ.1) THEN
                WRITE(LFNX,1210)
     +            (LABEL(IB,K),K=1,3),NAM(1),ICH(1,1),ICH(1,2)
              ELSE
                WRITE(LFNX,1220)  (LABEL(IB,K),K=1,3),
     +            (NAM(K),ICH(K,1),ICH(K,2),ISP(K),K=1,3)
              END IF
            END DO
C Then read them back in (removing blanks) and store in TLBL array
            REWIND LFNX
            DO IBAS=1,NNAO
              READ(LFNX,960) STR2
              K=0
              DO IS2=1,NLB
                IF(ICHAR(STR2(IS2:IS2)).GT.32) THEN
                  K=K+1
                  TLBL(IBAS)(K:K)=STR2(IS2:IS2)
                END IF
              END DO
              TLBL(IBAS)=TLBL(IBAS)(1:K)//BL(1:NLB-K)
            END DO
          ELSE
C This is a (P)AO, (P)NAO, or (P)NHO transformation; use the LBLT label
            KBAS = MAP(ABS(JPRINT(63)))
            DO IBAS=1,NNAO
              DO I=1,10
                IVEC(I)=LBLT(I,IBAS,KBAS)
              END DO
              NL=10
              CALL SCNVRT(IVEC,STR2,NL,IERR)
              K=0
              DO IS2=1,NLB
                IF(ICHAR(STR2(IS2:IS2)).GT.32) THEN
                  K=K+1
                  TLBL(IBAS)(K:K)=STR2(IS2:IS2)
                END IF
              END DO
              TLBL(IBAS)=TLBL(IBAS)(1:K)//BL(1:NLB-K)
            END DO
          END IF
        END IF
      END IF
C End of keylist?
      IF(END.OR.STR1(1:1).EQ.'$'.OR.STR1(1:3).EQ.'END'.OR.
     + (ISPIN.EQ.1.AND.STR1(1:4).EQ.'BETA').OR.
     + (ISPIN.EQ.2.AND.STR1(1:5).EQ.'ALPHA')) THEN
        NL1=1
        STR1(1:1)=' '
        END=.TRUE.
        GOTO 60
      END IF
      IF(IWLBL.EQ.1) THEN
C Convert to upper case and delete blanks; append to STR
        DO 50 IL1=1,NL1
        IC=ICHAR(STR1(IL1:IL1))
        IF(IC.EQ.32) THEN
          GOTO 50
        ELSE IF(IC.GE.97.AND.IC.LE.122) THEN
          IC=IC-32
        END IF
        NL=NL+1
        STR(NL:NL)=CHAR(IC)
   50   CONTINUE
C Can we find it in the label LIST?
        DO IB=1,NNAO
          NLT=LENNB(TLBL(IB))
          IF(NLT.EQ.NL.AND.STR(1:NL).EQ.TLBL(IB)(1:NL)) THEN
            IL=IL+1
            LIST(IL)=ITMP(IB)
            IORG(IL)=IB
C Is LIST(IL) different from any previous entry?
            IF(IL.GT.1) THEN
              DO JL=1,IL-1
                IF(LIST(JL).EQ.LIST(IL)) THEN
                  IERR=10
                  ERRFRAG='duplicate'
                  GOTO 130
                END IF
              END DO
            END IF
C Successful LIST entry, zero the ITMP entry
            ITMP(IB)=0
            GOTO 30
          END IF
        END DO
C No match yet.  Try to complete the label with another field?
        IF(NL.LT.NLB.AND..NOT.END) GOTO 40
        IERR=7
        ERRFRAG=STR(1:12)
        GOTO 130
      ELSE
        ITRY=NBINTV(STR1)
        IF(ITRY.GE.1.AND.ITRY.LE.NNAO) THEN
          IL=IL+1
          LIST(IL)=ITMP(ITRY)
          IORG(IL)=ITRY
C Is LIST(IL) different from any previous entry?
          IF(IL.GT.1) THEN
            DO JL=1,IL-1
              IF(LIST(JL).EQ.LIST(IL)) THEN
                IERR=10
                ERRFRAG='duplicate'
                GOTO 130
              END IF
            END DO
          END IF
C Successful LIST entry, remove from ITMP array
          ITMP(ITRY)=0
          GOTO 30
        END IF
      END IF
      IERR=5
      GOTO 130
C Ready to alter the permutation vector
   60 CONTINUE
C Loop over list elements
      DO 70 I=1,IL
      IF(IHFERMI.EQ.0) THEN
        IP=I
      ELSE
        IP=NLEW+I-IFERMI
        IF(IP.LT.1.OR.IP.GT.NNAO) THEN
          IERR=9
          ERRFRAG='Bad | loc.'
          GOTO 130
        END IF
      END IF
      IRNK(IP)=LIST(I)
      IARC(IP)=IORG(I)
   70 CONTINUE
C Load in remaining elements, in order.
      IP=0
      DO 90 I=1,NNAO
      IF(IRNK(I).EQ.0) THEN
   80   IP=IP+1
        IF(ITMP(IP).NE.0) THEN
          IRNK(I)=ITMP(IP)
          ITMP(IP)=0
          IARC(I)=IP
        ELSE IF(IP.LT.NNAO) THEN
          GOTO 80
        ELSE
          IERR=12
          ERRFRAG='Lost ITMP'
          GOTO 130
        END IF
      END IF
   90 CONTINUE
C
C Successful LFNPRM list: ready to permute the transformation array
C
      IF(ABS(JPRINT(63)).EQ.1) WRITE(LFNPR,901)
      IF(ABS(JPRINT(63)).EQ.2) WRITE(LFNPR,902)
      IF(ABS(JPRINT(63)).EQ.3) WRITE(LFNPR,903)
      IF(ABS(JPRINT(63)).EQ.4) WRITE(LFNPR,904)
      IF(ABS(JPRINT(63)).EQ.5) WRITE(LFNPR,905)
      IF(ABS(JPRINT(63)).EQ.6) WRITE(LFNPR,906)
      IF(ABS(JPRINT(63)).EQ.7) WRITE(LFNPR,907)
      IF(ABS(JPRINT(63)).EQ.8) WRITE(LFNPR,908)
      IF(ABS(JPRINT(63)).EQ.9) WRITE(LFNPR,909)
C
C  Retrieve the checkpoint matrices from the DAF:
C
      ISPIN = -2
  100 ISPIN = -ISPIN
      CALL FECHK(T,ISPIN,IT)
      IF(IT.NE.0) THEN
C
C  Stuff sorted orbitals in S and write to DAF:
C
        DO J = 1,NNAO
          JP=IARC(J)
          DO I = 1,NNAO
            S(I,J) = T(I,JP)
          ENDDO
        ENDDO
        CALL SVCHK(S,ISPIN)
C
C  Write orbital sequence:
C
        IF(ALPHA) WRITE(LFNPR,930)
        IF(BETA)  WRITE(LFNPR,940)
        K = MAP(ABS(JPRINT(63)))
        WRITE(LFNPR,950) (J,(LBLT(I,IARC(J),K),I=1,10),J=1,NNAO)
      END IF
      IF(ISPIN.GT.0) GOTO 100
      RETURN
C
C Fatal error! Ignore <PERM> list and exit
C
  130 CONTINUE
      WRITE(LFNPR,1100) IERR,ERRFRAG
      RETURN
  901 FORMAT(//1X,'Reorder PAOs for checkpoint file ',
     + 'from <PERM> list.',/)
  902 FORMAT(//1X,'Reorder PNAOs for checkpoint file ',
     + 'from <PERM> list.',/)
  903 FORMAT(//1X,'Reorder NAOs for checkpoint file ',
     + 'from <PERM> list.',/)
  904 FORMAT(//1X,'Reorder PNHOs for checkpoint file ',
     + 'from <PERM> list.',/)
  905 FORMAT(//1X,'Reorder NHOs for checkpoint file ',
     + 'from <PERM> list.',/)
  906 FORMAT(//1X,'Reorder PNBOs for checkpoint file ',
     + 'from <PERM> list.',/)
  907 FORMAT(//1X,'Reorder NBOs for checkpoint file ',
     + 'from <PERM> list.',/)
  908 FORMAT(//1X,'Reorder PNLMOs for checkpoint file ',
     + 'from <PERM> list.',/)
  909 FORMAT(//1X,'Reorder NLMOs for checkpoint file ',
     + 'from <PERM> list.',/)
  930 FORMAT(1X,'Alpha orbitals:',/)
  940 FORMAT(1X,'Beta orbitals:',/)
  950 FORMAT(5(I4,'. ',10A1))
  960 FORMAT(A21)
 1100 FORMAT(1X,'*** Fatal error ',I3,' in <PERM> list ',
     + 'at ',A12,'; keylist ignored.')
 1210 FORMAT(1X,A2,A1,'(',I2,')',A2,2A1)
 1220 FORMAT(1X,A2,A1,'(',I2,')',3(A2,3A1))
      END
C***********************************************************************
      SUBROUTINE RDPRM2(TA,TB)
C***********************************************************************
C 24-Feb-01  FAW  New subroutine to perform spin-bracket AOPNBO permutations
      IMPLICIT REAL*8 (A-H,O-Z)
C  Perform the AOPNBO permutations for spin bracket list
      DIMENSION TA(NBAS,NBAS),TB(NBAS,NBAS)
      COMMON/NBELEM/IJELEM(2,100),NELEM
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      COMMON/NBPRMC/LFNPRM
C  Fetch the checkpointed (AOPNBO) transforms
      IF(OPEN.AND.BETA) THEN
        WRITE(LFNPR,1100)
        CALL FECHK(TA,2,ITA)
        CALL FECHK(TB,-2,ITB)
        IF(ITA.EQ.0.OR.ITB.EQ.0) GOTO 100
      ELSE
        IF(OPEN) RETURN
        GOTO 100
      ENDIF
      IF(NELEM.LE.0) RETURN
      DO IELEM = 1,NELEM
        I = IJELEM(1,IELEM)
        J = IJELEM(2,IELEM)
        IC = ABS(I)
        JC = ABS(J)
C  Exchange two alpha orbitals?
        IF(I.GT.0.AND.J.GT.0) THEN
          DO K = 1,NBAS
            T = TA(K,IC)
            TA(K,IC) = TA(K,JC)
            TA(K,JC) = T
          ENDDO
C  Exchange two beta orbitals?
        ELSE IF(I.LT.0.AND.J.LT.0) THEN
          DO K = 1,NBAS
            T = TB(K,IC)
            TB(K,IC) = TB(K,JC)
            TB(K,JC) = T
          ENDDO
C  Exchange alpha(I) for beta(J)?
        ELSE IF(I.GT.0.AND.J.LT.0) THEN
          DO K = 1,NBAS
            T = TA(K,IC)
            TA(K,IC) = TB(K,JC)
            TB(K,JC) = T
          ENDDO
C  Exchange alpha(J) for beta(I)?
        ELSE IF(I.LT.0.AND.J.GT.0) THEN
          DO K = 1,NBAS
            T = TB(K,IC)
            TB(K,IC) = TA(K,JC)
            TA(K,JC) = T
          ENDDO
        ENDIF
      ENDDO
      CALL SVCHK(TA,2)
      CALL SVCHK(TB,-2)
      WRITE(LFNPR,1200)
C  Switch off the alternative LFNPRM permutation
      LFNPRM = 0
      RETURN
  100 WRITE(6,1000)
 1000 FORMAT(1X,'*** Permutation failed.')
 1100 FORMAT(/,1X,'Attempt user-requested permutations of ',
     + 'checkpointed AOPNBO transform.')
 1200 FORMAT(1X,'Permutation completed.')
      RETURN
      END
C***********************************************************************
      SUBROUTINE FILFLD(LFN,STR,LENGTH,IAMNEW,END)
C***********************************************************************
C 31-Oct-96 FW  SAVEd ICD array to correct unix bug (per EDG)
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER ICD*80,STR*(*),C*1
      LOGICAL END
C
      SAVE IPT,ICD
      DATA IPT/1/
C
C  FIND NEXT NON-BLANK STRING OF CHARACTERS IN LFN.  READ IN ANOTHER LINE
C  OF LFN UNTIL NON-BLANK CHARACTERS ARE FOUND AND PLACE THEM IN "STR",
C  OF LENGTH "LENGTH".  SET IAMNEW=1 FOR THE FIRST READ FROM A NEW FILE.
C
C  FIRST ENTRY?
C
      IF(IAMNEW.EQ.1) THEN
        IAMNEW = 0
        IPT = 80
        REWIND LFN
        END = .FALSE.
      END IF
      IF(END) GO TO 35
      IF(IPT.GE.80) THEN
        READ(LFN,1000,END=35) ICD
        IPT = 1
      END IF
 1000 FORMAT(A)
C
C  LOOK FOR START OF FIELD.  SKIP TO NEXT CARD IF "!" IS ENCOUNTERED
C  (COMMENT FIELD):
C
   10 CONTINUE
      DO 20 NCOL = IPT,80
        C = ICD(NCOL:NCOL)
        IF(C.EQ.'!') GOTO 30
        IF(C.NE.' '.AND.C.NE.','.AND.C.NE.'=') GOTO 40
   20 CONTINUE
C
C  NOTHING ADDITIONAL FOUND ON THIS CARD, CONTINUE WITH THE NEXT CARD:
C
   30 READ(LFN,1000,END=35) ICD
      IPT = 1
      IF(.NOT.END) GO TO 10
C
C  END OF FILE FOUND:
C
   35 LENGTH = 0
      END = .TRUE.
      RETURN
C
C  LOOK FOR THE END OF THIS FIELD, COUNTING CHARACTERS AS WE GO AND
C  STORING THESE CHARACTERS IN STR:
C
   40 M = 0
      DO 80 MCOL = NCOL,80
        C = ICD(MCOL:MCOL)
        IF(C.EQ.' '.OR.C.EQ.','.OR.C.EQ.'=') GOTO 100
        M = M + 1
        STR(M:M) = C
   80 CONTINUE
C
C  SET LENGTH TO THE LENGTH OF THE NEW STRING IN STR AND RESET IPT TO
C  THE NEXT SPACE AFTER THIS STRING:
C
  100 LENGTH = M
      IPT = MCOL
      RETURN
      END
C***********************************************************************
C NBO 5.G -- Natural Bond Orbital Analysis Programs
C (c) Copyright 1996-2008 Board of Regents of the University of Wisconsin System
C     on behalf of the Theoretical Chemistry Institute.  All Rights Reserved.
C***********************************************************************
C
C  ANALYSIS OF CANONICAL MOLECULAR ORBITALS:
C
C      SUBROUTINE CMOANL(TMO,TMP,SCR,IFLG)
C
C***********************************************************************
      SUBROUTINE CMOANL(TMO,TMP,SCR,IFLG)
C***********************************************************************
C 24-Jul-1997 FAW  New subroutine
C  5-Feb-2003 FAW  Correct MO numbering error in atom-atom bond 
C                  character table
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
C
C  Perform NBO analysis of canonical molecular orbitals
C
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LTYP(MAXBAS),IATHY(MAXBAS,3)
      COMMON/NBOPT/IWDM,IW3C,IWAPOL,IWHYBS,IWPNAO,IWTNAO,IWTNAB,
     + IWTNBO,IWFOCK,IWCUBF,IPSEUD,IPRINT,IWDETL,IWMULP,ICHOOS,
     + IWNBBP,IWMSP,IWFIXDM,IW3CHB,IWNJC,JCORE,JPRINT(100)
      COMMON/NBIO/LFNIN,LFNPR,LFNAO,LFNPNA,LFNNAO,LFNPNH,LFNNHO,LFNPNB,
     +            LFNNBO,LFNPNL,LFNNLM,LFNMO,LFNDM,LFNNAB,LFNPPA,LFNARC,
     +            LFNDAF,LFNLBL,LFNDEF,LFNBRK(100)
      DIMENSION TMO(NDIM,NDIM),TMP(NDIM,NDIM),
     +          SCR(NDIM,NDIM),EORB(MAXBAS)
      DIMENSION BCBD(MAXATM),BCNB(MAXATM),BCAB(MAXATM),
     + IATBD(3,MAXATM),IATNB(MAXATM),IATAB(3,MAXATM),
     + ILBD(MAXBAS),ILNB(MAXBAS),ILAB(MAXBAS)
      DIMENSION IARCBD(MAXBAS),IARCNB(MAXBAS),IARCAB(MAXBAS),
     + IRNK(MAXBAS),IARC(MAXBAS),LIST(MAXBAS),IDIG(9),ISHELL(4)
      CHARACTER*80 TITLE,TYPE*3
      DATA IFULL,IVAL,ILEW/4HFULL,3HVAL,3HLEW/
      DATA THRESH/0.05D0/ZERO/0.0D0/
      DATA IHS/1H*/IHX/1H /IHLP/1H(/IHRP/1H)/
     + IHO/1Ho/IHV/1Hv/IHE/1H=/
      SAVE IFULL,IVAL,ILEW,THRESH,ZERO,IHS,IHX,
     + IHLP,IHRP,IHO,IHV,IHE
C Initialize arrays:
      DO I = 1,NATOMS
        BCBD(I) = ZERO
        BCNB(I) = ZERO
        BCAB(I) = ZERO
        DO J = 1,3
          IATBD(J,I) = 0
          IATAB(J,I) = 0
        END DO
        IATNB(I) = 0
      END DO
C How many MOs to be analyzed?
      NOCC=NLEW
      IO = IFLG
      IF(IO.EQ.ILEW) THEN
        NMOPR = NOCC
      ELSE IF(IO.EQ.IVAL) THEN
C
C  If IO = IVAL, output only the valence orbitals, determined from the
C  core and valence tables:
C
        IF(NVAL.LT.0) THEN
          IECP = 0
          IO = 0
          DO IAT = 1,NATOMS
            CALL CORTBL(IAT,ISHELL,IECP)
            DO I = 1,4
              MULT = 2*(I-1) + 1
              IO = IO + ISHELL(I)*MULT
            ENDDO
            CALL VALTBL(IAT,ISHELL,0)
            DO I = 1,4
              MULT = 2*(I-1) + 1
              IO = IO + ISHELL(I)*MULT
            ENDDO
          ENDDO
        ELSE
          IO = NVAL
        ENDIF
        NMOPR = IO
      ELSE IF(IO.EQ.IFULL) THEN
        NMOPR = NNAO
      ELSE
        NMOPR = IO
        IF(NMOPR.GT.NNAO) NMOPR = NNAO
      ENDIF
C Print title and headings
      WRITE(LFNPR,8000)
      WRITE(LFNPR,8000)
 8000 FORMAT(1X)
      TITLE='CMO: NBO Analysis of Canonical Molecular Orbitals'
      NT=LENNB(TITLE)
      WRITE(LFNPR,8100) (TITLE(K:K),K=1,NT)
 8100 FORMAT(1X,80A1)
C
C  Prepare orbital energies in EORB
C
C  Fetch the Fock matrix in the AO basis, store it in SCR:
C
      CALL FEFAO(SCR,IWFOCK)
C
C  Fetch the AO to MO transformation matrix (using TMP as storage):
C
      CALL FEAOMO(TMP,IT)
      IF(IT.EQ.0) RETURN
C
C  Diagonalize SCR (using matrix TMP of eigenvectors) to find orbital energies,
C
      CALL SIMTR1(SCR,TMP,EORB,NDIM,NBAS,NBAS,NNAO)
      DO I=1,NNAO
        EORB(I)=SCR(I,I)
      ENDDO
C
C  Find the NBO-MO transformation matrix in TMO:
C
C Print out the leading NBO contributions to CMO:
C
      TITLE='Leading (> 5%) NBO Contributions to Molecular Orbitals'
      NT=LENNB(TITLE)
      WRITE(LFNPR,8000)
      WRITE(LFNPR,8100) (TITLE(K:K),K=1,NT)
      WRITE(LFNPR,8100) (IHE,K=1,NT)
C Grand loop over MOs
      DO IMO = 1,NMOPR
        IF(IMO.LE.NOCC) THEN
          TYPE = 'occ'
        ELSE
          TYPE = 'vir'
        ENDIF
C Prepare list of squared LCNBO-MO coefficients
        NBIG = 0
        DO IBO = 1,NNAO
          T2 = TMO(IBO,IMO)**2
          IF(T2.GE.THRESH) NBIG = NBIG + 1
          LIST(IBO) = INT(10000*T2)
          IF(LIST(IBO).EQ.0) LIST(IBO) = 1
        ENDDO
        IF(NBIG.EQ.0) NBIG = 1
        CALL NBORDR(IRNK,LIST,NNAO,MAXBAS,IARC)
        WRITE(LFNPR,8200) IMO,TYPE,EORB(IMO)
 8200   FORMAT(2X,'MO',I4,' (',A3,'): orbital energy =',F11.6,' a.u.')
        DO I = NNAO,NNAO-NBIG+1,-1
          WRITE(LFNPR,8300) TMO(IARC(I),IMO),IARC(I),
     +     (LABEL(IBXM(IARC(I)),J),J=1,3),(LBL(J,IARC(I),4),J=1,10)
        ENDDO
 8300   FORMAT(15X,F6.3,'*[',I3,']: ',A2,A1,'(',I2,')',10A1)
      ENDDO
C
C  Print table of MO atom-atom bonding character
C
      WRITE(LFNPR,8400)
 8400 FORMAT(/,
     +'           Molecular Orbital Atom-Atom Bonding Character   ',/,
     +'          ================================================',/,
     +'              bonding        nonbonding      antibonding',/,
     +'    MO        (2c, 3c)        (1c, 1c*)       (2c*, 3c*)',/,
     +'  ------  ---------------    ----------    ---------------')
C Grand loop over MOs
      WFBD = ZERO
      WFNB = ZERO
      WFAB = ZERO
      DO IMO = 1,NMOPR
        CALL IDIGIT(IMO,IDIG,ND,4)
        IF(ND.LT.4) THEN
          DO ID = ND,1,-1
            IDIG(ID+4-ND) = IDIG(ID)
          ENDDO
          DO ID = 1,4-ND
            IDIG(ID) = IHX
          ENDDO
        ENDIF
        IDIG(5) = IHLP
        IF(IMO.LE.NOCC) THEN
          IDIG(6) = IHO
        ELSE
          IDIG(6) = IHV
        ENDIF
        IDIG(7) = IHRP
        IDIG(8) = IHX
        IDIG(9) = IHX
C Determine the list of bonding, nonbonding, and antibonding
C indices (BCBD, BCNB, BCAB) for this MO
        TOTBD = ZERO
        IBD=0
        TOTNB = ZERO
        INB=0
        TOTAB = ZERO
        IAB=0
        DO IBO = 1,NNAO
          IBOL = IBXM(IBO)
          T2 = TMO(IBO,IMO)**2
C Nonbonding (single-center) contribution?
          IF(LABEL(IBOL,5).EQ.0.AND.LABEL(IBOL,6).EQ.0) THEN
            TOTNB = TOTNB + T2
            IF(T2.GE.THRESH) THEN
              IDONE = 0
              DO I = 1,INB
                IF(LABEL(IBOL,4).EQ.IATNB(I)) THEN
                  IDONE = 1
                  BCNB(I) = BCNB(I) + T2
                ENDIF
              ENDDO
              IF(IDONE.EQ.0) THEN
                INB = INB + 1
                BCNB(INB) = T2
                IATNB(INB) = LABEL(IBOL,4)
                ILNB(INB) = IBO
              ENDIF
            ENDIF
C Bonding (2c, 3c) contribution?
          ELSE IF(LABEL(IBOL,2).NE.IHS) THEN
            TOTBD = TOTBD + T2
            IF(T2.GE.THRESH) THEN
              IDONE = 0
              DO I = 1,IBD
                IF(IATBD(1,I).EQ.LABEL(IBOL,4).AND.
     +             IATBD(2,I).EQ.LABEL(IBOL,5).AND.
     +             IATBD(3,I).EQ.LABEL(IBOL,6)) THEN
                  IDONE = 1
                  BCBD(I) = BCBD(I) + T2
                ENDIF
              ENDDO
              IF(IDONE.EQ.0) THEN
                IBD = IBD + 1
                BCBD(IBD) = T2
                IATBD(1,IBD) = LABEL(IBOL,4)
                IATBD(2,IBD) = LABEL(IBOL,5)
                IATBD(3,IBD) = LABEL(IBOL,6)
                ILBD(IBD) = IBO
              ENDIF
            ENDIF
C Antibonding (2c*, 3c*) contribution?
          ELSE
            TOTAB = TOTAB + T2
            IF(T2.GE.THRESH) THEN
              IDONE = 0
              DO I = 1,IAB
                IF(IATAB(1,I).EQ.LABEL(IBOL,4).AND.
     +             IATAB(2,I).EQ.LABEL(IBOL,5).AND.
     +             IATAB(3,I).EQ.LABEL(IBOL,6)) THEN
                  IDONE = 1
                  BCAB(I) = BCAB(I) + T2
                ENDIF
              ENDDO
              IF(IDONE.EQ.0) THEN
                IAB = IAB + 1
                BCAB(IAB) = T2
                IATAB(1,IAB) = LABEL(IBOL,4)
                IATAB(2,IAB) = LABEL(IBOL,5)
                IATAB(3,IAB) = LABEL(IBOL,6)
                ILAB(IAB) = IBO
              ENDIF
            ENDIF
          ENDIF
        ENDDO
        MAXNUM = MAX(IBD,INB,IAB)
C Rank the entries for each type [I = IARCXX(N) is the entry of rank N,
C lowest = rank 1,..., highest = rank(IXX)]
C
C Order the bonding characters
        DO I = 1,IBD
          LIST(I) = INT(10000*BCBD(I))
        ENDDO
        CALL NBORDR(IRNK,LIST,IBD,MAXBAS,IARC)
        DO I = 1,IBD
          IARCBD(I) = IARC(IBD-I+1)
        ENDDO
C Order the nonbonding characters
        DO I = 1,INB
          LIST(I) = INT(10000*BCNB(I))
        ENDDO
        CALL NBORDR(IRNK,LIST,INB,MAXBAS,IARC)
        DO I = 1,INB
          IARCNB(I) = IARC(INB-I+1)
        ENDDO
C Order the antibonding characters
        DO I = 1,IAB
          LIST(I) = INT(10000*BCAB(I))
        ENDDO
        CALL NBORDR(IRNK,LIST,IAB,MAXBAS,IARC)
        DO I = 1,IAB
          IARCAB(I) = IARC(IAB-I+1)
        ENDDO
C Print the ordered bonding characters for this MO
        DO INUM = 1,MAXNUM
          IF(INUM.LE.IBD) THEN
            JBD = IARCBD(INUM)
            JBDP = ILBD(JBD)
          ELSE
            JBD = 0
          ENDIF
          IF(INUM.LE.INB) THEN
            JNB = IARCNB(INUM)
            JNBP = ILNB(JNB)
          ELSE
            JNB = 0
          ENDIF
          IF(INUM.LE.IAB) THEN
            JAB = IARCAB(INUM)
            JABP = ILAB(JAB)
          ELSE
            JAB = 0
          ENDIF
C Blank the IDIG (IMO) label after the first line
          IF(INUM.GT.1) THEN
            DO I = 1,9
              IDIG(I) = IHX
            ENDDO
          ELSE IF(IMO.NE.1) THEN
            WRITE(LFNPR,8500)
 8500 FORMAT(2X,56('-'))
          ENDIF
C Print the table entry
          IF(JBD.GT.0.AND.JNB.GT.0.AND.JAB.GT.0) THEN
            WRITE(LFNPR,8600) (IDIG(J),J=1,9),
     +        BCBD(JBD),(LBL(J,JBDP,4),J=1,9),
     +        BCNB(JNB),(LBL(J,JNBP,4),J=1,4),
     +        BCAB(JAB),(LBL(J,JABP,4),J=1,9)
 8600 FORMAT(1X,9A1,F5.3,1X,9A1,4X,F5.3,1X,4A1,4X,F5.3,1X,9A1)
          ELSE IF(JBD.GT.0.AND.JNB.GT.0.AND.JAB.EQ.0) THEN
            WRITE(LFNPR,8700) (IDIG(J),J=1,9),
     +        BCBD(JBD),(LBL(J,JBDP,4),J=1,9),
     +        BCNB(JNB),(LBL(J,JNBP,4),J=1,4)
 8700 FORMAT(1X,9A1,F5.3,1X,9A1,4X,F5.3,1X,4A1,4X,15X)
          ELSE IF(JBD.GT.0.AND.JNB.EQ.0.AND.JAB.GT.0) THEN
            WRITE(LFNPR,8800) (IDIG(J),J=1,9),
     +        BCBD(JBD),(LBL(J,JBDP,4),J=1,9),
     +        BCAB(JAB),(LBL(J,JABP,4),J=1,9)
 8800 FORMAT(1X,9A1,F5.3,1X,9A1,4X,10X,4X,F5.3,1X,9A1)
          ELSE IF(JBD.EQ.0.AND.JNB.GT.0.AND.JAB.GT.0) THEN
            WRITE(LFNPR,8900) (IDIG(J),J=1,9),
     +        BCNB(JNB),(LBL(J,JNBP,4),J=1,4),
     +        BCAB(JAB),(LBL(J,JABP,4),J=1,9)
 8900 FORMAT(1X,9A1,15X,4X,F5.3,1X,4A1,4X,F5.3,1X,9A1)
          ELSE IF(JBD.GT.0.AND.JNB.EQ.0.AND.JAB.EQ.0) THEN
            WRITE(LFNPR,9000) (IDIG(J),J=1,9),
     +        BCBD(JBD),(LBL(J,JBDP,4),J=1,9)
 9000 FORMAT(1X,9A1,F5.3,1X,9A1)
          ELSE IF(JBD.EQ.0.AND.JNB.GT.0.AND.JAB.EQ.0) THEN
            WRITE(LFNPR,9100) (IDIG(J),J=1,9),
     +        BCNB(JNB),(LBL(J,JNBP,4),J=1,4)
 9100 FORMAT(1X,9A1,15X,4X,F5.3,1X,4A1,4X,15X)
          ELSE IF(JBD.EQ.0.AND.JNB.EQ.0.AND.JAB.GT.0) THEN
            WRITE(LFNPR,9200) (IDIG(J),J=1,9),
     +        BCAB(JAB),(LBL(J,JABP,4),J=1,9)
 9200 FORMAT(1X,9A1,15X,4X,10X,4X,F5.3,1X,9A1)
          ENDIF
        ENDDO
        WFBD = WFBD + TOTBD
        WFNB = WFNB + TOTNB
        WFAB = WFAB + TOTAB
        WRITE(LFNPR,9300) TOTBD,TOTNB,TOTAB
 9300 FORMAT(10X,'_____',14X,'_____',9X,'_____',/,10X,
     + F5.3,'(b)',11X,F5.3,'(n)',6X,F5.3,'(a)  total')
      ENDDO
      WRITE(LFNPR,8500)
      WRITE(LFNPR,9400) WFBD,WFNB,WFAB
 9400 FORMAT(8X,F7.3,'(b)',9X,F7.3,'(n)',4X,F7.3,
     + '(a)  Sum total for MOs')
      RETURN
      END
C***********************************************************************
C
C  NATURAL NMR CHEMCIAL SHIELDING TENSOR ANALYSIS:
C     (LFN 26-30,91 ARE USED BY NCS AS SCRATCH FILES)
C
C      SUBROUTINE NMRWR(D,IX,JY,KZ,IFI,NCF)
C      SUBROUTINE FEH01(IAT,IY,H01I,NBAS,IERR)
C      SUBROUTINE FEH11(IAT,IX,JY,H11I,NBAS,IERR)
C      SUBROUTINE NCSATM(IAT,A,B,D,T,C,C10X,C10Y,C10Z,NBAS,NOCC,LFNU,
C     +           LFNI,SCR1,SCR2,SCR3,SCR4,H01J,IWCMO)
C      SUBROUTINE NCSUAJ(IAT,JY,A,B,D,C,NBAS,NOCC,LFNU,H11XJ,H11YJ,
C     +           H11ZJ,SCRA,IWCMO)
C      SUBROUTINE NCSIAJ(IAT,JY,A,B,T,C,C10X,C10Y,C10Z,NBAS,NOCC,LFNI,
C     +           H01JY,SCRX,SCRY,SCRZ,SCR,IWCMO)
C      FUNCTION TRCNCS(A,B,NDIM)
C      SUBROUTINE NCSXYZ(IAT,NBAS,NOCC,SIGXL,SIGYL,SIGZL,LFNU,LFNI,
C     +           LFNPR,LFNISO,IWXYZ,IWCSA,IWD,IWP,IZ,R,IORDER,IWCMO)
C      SUBROUTINE NCSCSA(IAT,SIGTOT,LFNX,LFNY,LFNZ,LFNPR,
C     +           IWCSA,IZ,ITYP,R,IORDER,IWCMO)
C      SUBROUTINE NCSISO(NOCC,NATM,NUC,SIGISO,LFNISO,LFNPR,IWCMO)
C
C***********************************************************************
      SUBROUTINE NMRWR (D,IX,JY,KZ,IFI,NCF)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C Called by NCSMAT
C Writes files for H11 and H01 matrices for NMR
C i,j,k are integer codes associated with  a record.
C  i=atom number ,j = 1st cartesian component, k=second cartesian
C  component
C the data for upper diagonal h11 and h01 matrices
      REAL*8 D(NCF*(NCF+1)/2)
      N=NCF
      JA=1
      WRITE (IFI) IX,JY,KZ
      DO 10 I=1,N
      JE=JA+I-1
      WRITE (IFI) (D(J),J=JA,JE)
      JA=JE+1
   10 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEH01(IAT,IY,H01I,NBAS,IERR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      COMMON/NMRDAT/LFNH11,LFNH01
      DIMENSION H01I(NBAS,NBAS)
C
      IERR = 0
      J = 1
      J = J
      READ(LFNH01) I,J,K
C
      DO 10 II=1,NBAS
        READ (LFNH01) (H01I(II,JJ),JJ=1,II)
   10 CONTINUE
      if (K.ne.IY.OR.I.NE.IAT)then
         ierr=1
         CALL NBHALT('Coordinate mismatch in FEH01.')
      endif
      DO 30 II=1,NBAS-1
        DO 20 JJ=II+1,NBAS
          H01I(II,JJ)=-H01I(JJ,II)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE FEH11(IAT,IX,JY,H11I,NBAS,IERR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NMRDAT/LFNH11,LFNH01
      DIMENSION H11I(NBAS,NBAS)
      IERR = 0
      READ(LFNH11) I,J,K
      DO 10 II=1,NBAS
        READ (LFNH11)  (H11I(II,JJ),JJ=1,II)
   10 CONTINUE
      if (J.NE.IX.OR.K.ne.jy.OR.I.NE.IAT)then
         ierr=1
         CALL NBHALT('Coordinate mismatch in FEH11.')
      endif
C compute upper half
      DO 30 II=1,NBAS-1
        DO 20 JJ=II+1,NBAS
          H11I(II,JJ)=H11I(JJ,II)
C           PRINT*, 'H11I ',II,JJ,H11I(II,JJ)
   20   CONTINUE
   30 CONTINUE
      RETURN
      END
C***********************************************************************
      SUBROUTINE NCSATM(IAT,A,B,D,T,C,C10X,C10Y,C10Z,NBAS,NOCC,LFNU,
     +LFNI,SCR1,SCR2,SCR3,SCR4,H01J,IWCMO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION B(NBAS,NBAS),A(NBAS,NBAS),D(NBAS,NBAS),T(NBAS,NBAS)
      DIMENSION C10X(NBAS,NOCC),C10Y(NBAS,NOCC),C10Z(NBAS,NOCC)
      DIMENSION SCR1(NBAS,NBAS),SCR2(NBAS,NBAS),SCR3(NBAS,NBAS),
     + SCR4(NBAS,NBAS),H01J(NBAS,NBAS),C(NBAS,NBAS)
C Rewind the scratch files for a new atom
      REWIND (LFNU)
      REWIND (LFNI)
      DO JY=1,3
C Compute unperturbed contributions for atom IAT (store in LFNU)
        CALL NCSUAJ (IAT,JY,A,B,D,C,NBAS,NOCC,LFNU,SCR1,SCR2,SCR3,SCR4,
     +  IWCMO)
      END DO
C Compute induced contributions for atom IAT (store in LFNI)
      DO JY=1,3
        CALL FEH01(IAT,JY,H01J,NBAS,IERR)
        CALL NCSIAJ (IAT,JY,A,B,T,C,C10X,C10Y,C10Z,NBAS,NOCC,LFNI,
     +  H01J,SCR1,SCR2,SCR3,SCR4,IWCMO)
      END DO
      RETURN
      END
C***********************************************************************
      SUBROUTINE NCSUAJ(IAT,JY,A,B,D,C,NBAS,NOCC,LFNU,H11XJ,H11YJ,
     +H11ZJ,SCRA,IWCMO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C Compute "unperturbed" components SIGD(IX,JY), IX=1,3, for atom IAT;
C Store results in LFNU
      DIMENSION B(NBAS,NBAS),A(NBAS,NBAS),D(NBAS,NBAS)
      DIMENSION H11XJ(NBAS,NBAS),H11YJ(NBAS,NBAS),H11ZJ(NBAS,NBAS),
     + SCRA(NBAS,NBAS),C(NBAS,NBAS)
C Get the H11 tensor components (IX,JY) in the AO basis for this nucleus
      CALL FEH11(IAT,1,JY,H11XJ,NBAS,IERR)
      CALL FEH11(IAT,2,JY,H11YJ,NBAS,IERR)
      CALL FEH11(IAT,3,JY,H11ZJ,NBAS,IERR)
C Calculate CMO contributions if desired
C This is somewhat laborious, but is saves having to store
C intermediate results.
      IF (IWCMO.EQ.1) THEN
        DO J=1,NOCC
C Evaluate the D00L array for this MO
C (Diagonalizing D00H11 Matrix might be faster)
          DO IP=1,NBAS
            DO IQ=1,NBAS
              SCRA(IP,IQ)=C(IP,J)*C(IQ,J)
            END DO
          END DO
C Contract with H11 matrices to get the localized MO contributions
C Include factor of two for density matrix in printing routine,not here
          SIGDXJN=TRCNCS(SCRA,H11XJ,NBAS)
          SIGDYJN=TRCNCS(SCRA,H11YJ,NBAS)
          SIGDZJN=TRCNCS(SCRA,H11ZJ,NBAS)
C Write out the MO contribution(N=0)
          WRITE (LFNU) IAT,JY,J,0,SIGDXJN,SIGDYJN,SIGDZJN
        END DO
      END IF
C Loop over NBO indices to evaluate localized contributions (J<-L,N<-NL)
C need to calculate elements of the matrix B(transpose)*H11*D
C   calculate H11*D, store result in H11
      CALL MATMLT(H11XJ,D,SCRA,NBAS,NBAS)
      CALL MATMLT(H11YJ,D,SCRA,NBAS,NBAS)
      CALL MATMLT(H11ZJ,D,SCRA,NBAS,NBAS)
C   Transpose the quantity H11*D stored in H11
      CALL NBTRSP(H11XJ,NBAS,NBAS)
      CALL NBTRSP(H11YJ,NBAS,NBAS)
      CALL NBTRSP(H11ZJ,NBAS,NBAS)
C   calculate the transpose of our desired product,
C  ( B(transpose)*H11*D) quantity transposed, store result in H11
      CALL MATMLT(H11XJ,B,SCRA,NBAS,NBAS)
      CALL MATMLT(H11YJ,B,SCRA,NBAS,NBAS)
      CALL MATMLT(H11ZJ,B,SCRA,NBAS,NBAS)
C calculate orbital contributions.
C H11 matrices are now the transposed representation,
C so reverse indices when using H11.
      DO J=1,NOCC
        DO N=1,NBAS
            ANJ=A(N,J)
            SIGDXJN=ANJ*H11XJ(J,N)
            SIGDYJN=ANJ*H11YJ(J,N)
            SIGDZJN=ANJ*H11ZJ(J,N)
C Write out the localized contribution
            WRITE (LFNU) IAT,JY,J,N,SIGDXJN,SIGDYJN,SIGDZJN
        END DO
      END DO
C END OF NBO/NLMO ANALYSIS
      RETURN
      END
C***********************************************************************
      SUBROUTINE NCSIAJ(IAT,JY,A,B,T,C,C10X,C10Y,C10Z,NBAS,NOCC,LFNI,
     +H01JY,SCRX,SCRY,SCRZ,SCR,IWCMO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C Compute "induced" components SIGP(IX,JY), IX=1,3 for atom IAT;
C Store results in LFNI
C As of 3/28/96, this is the routine which seems to take the most time
C when NCSANL runs.
      DIMENSION B(NBAS,NBAS),A(NBAS,NBAS),T(NBAS,NBAS),C(NBAS,NBAS)
      DIMENSION H01JY(NBAS,NBAS),SCR(NBAS,NBAS)
      DIMENSION SCRX(NBAS,NBAS),SCRY(NBAS,NBAS),SCRZ(NBAS,NBAS)
      DIMENSION C10X(NBAS,NOCC),C10Y(NBAS,NOCC),C10Z(NBAS,NOCC)
      DATA ZERO/0.0D0/
C Calculate canonical MO contributions if desired
       IF (IWCMO.EQ.1) THEN
       DO K = 1,NOCC
C Evaluate the D10 matrix for each component, store as SCRX,etc.
        DO IP = 1,NBAS
          SCRX(IP,IP)=ZERO
          SCRY(IP,IP)=ZERO
          SCRZ(IP,IP)=ZERO
          CPK=C(IP,K)
          CPKX=C10X(IP,K)
          CPKY=C10Y(IP,K)
          CPKZ=C10Z(IP,K)
          DO IQ = IP+1,NBAS
            CQK=C(IQ,K)
            TEMPX = CPK*C10X(IQ,K)-CPKX*CQK
            TEMPY = CPK*C10Y(IQ,K)-CPKY*CQK
            TEMPZ = CPK*C10Z(IQ,K)-CPKZ*CQK
            SCRX(IP,IQ) = TEMPX
            SCRX(IQ,IP) = -TEMPX
            SCRY(IP,IQ) = TEMPY
            SCRY(IQ,IP) = -TEMPY
            SCRZ(IP,IQ) = TEMPZ
            SCRZ(IQ,IP) = -TEMPZ
          ENDDO
        ENDDO
C Contract with H01 matrices to get the MO contributions
C The factor of two for the density matrix is done in the
C printing
        SIGPXJN=TRCNCS(SCRX,H01JY,NBAS)
        SIGPYJN=TRCNCS(SCRY,H01JY,NBAS)
        SIGPZJN=TRCNCS(SCRZ,H01JY,NBAS)
C Write out the MO contribution N=0
        WRITE(LFNI) IAT,JY,K,0,SIGPXJN,SIGPYJN,SIGPZJN
      ENDDO
      ENDIF
C START NBO/NLMO analysis
C Sigma_P(J->N)= A(N,J)*(T*C10I*(H01-H01(transpose))*B)(J,N)
C Thanks to P. Pulay for the suggestion
C First, calculate the matrix product (H01-H01(transpose))*B
      DO I=1,NBAS
        DO K=1,NBAS
          HBIK=ZERO
          DO L=1,NBAS
            HBIK=HBIK+(H01JY(I,L)-H01JY(L,I))*B(L,K)
          ENDDO
          SCR(I,K)=HBIK
        ENDDO
      ENDDO
C NEXT CALCULATE THE (NON-SQUARE) PRODUCTS T*C10
C The product has dimension (NOCC,NBAS)
      DO J=1,NOCC
        DO IQ=1,NBAS
          TCXJQ=ZERO
          TCYJQ=ZERO
          TCZJQ=ZERO
          DO K = 1,NOCC
            TJK=T(J,K)
            TCXJQ=TCXJQ+TJK*C10X(IQ,K)
            TCYJQ=TCYJQ+TJK*C10Y(IQ,K)
            TCZJQ=TCZJQ+TJK*C10Z(IQ,K)
          ENDDO
          SCRX(J,IQ)=TCXJQ
          SCRY(J,IQ)=TCYJQ
          SCRZ(J,IQ)=TCZJQ
        ENDDO
      ENDDO
C
C     multiply the two to get TCHB ,store in SCRX,Y,Z
C     this will overwrite H01JY (used now as scratch storage)
      DO J=1,NOCC
         DO N=1,NBAS
           TCHB=ZERO
           DO IQ=1,NBAS
             TCHB=TCHB+SCRX(J,IQ)*SCR(IQ,N)
           ENDDO
           H01JY(J,N)=TCHB
         ENDDO
      ENDDO
      CALL COPY(H01JY,SCRX,NBAS,NBAS,NBAS)
      DO J=1,NOCC
         DO N=1,NBAS
           TCHB=ZERO
           DO IQ=1,NBAS
             TCHB=TCHB+SCRY(J,IQ)*SCR(IQ,N)
           ENDDO
           H01JY(J,N)=TCHB
         ENDDO
      ENDDO
      CALL COPY(H01JY,SCRY,NBAS,NBAS,NBAS)
      DO J=1,NOCC
         DO N=1,NBAS
           TCHB=ZERO
           DO IQ=1,NBAS
             TCHB=TCHB+SCRZ(J,IQ)*SCR(IQ,N)
           ENDDO
           H01JY(J,N)=TCHB
         ENDDO
      ENDDO
      CALL COPY(H01JY,SCRZ,NBAS,NBAS,NBAS)
C calculate final products
C
C Temporary change for G03 : We have changed the signs of SIGPXJN, etc.
C contributions below until 
C repaired in new calculation of C10 Coeffients by NCS independently.
      DO J=1,NOCC
         DO N=1,NBAS
           ANJ=A(N,J)
           SIGPXJN=-ANJ*SCRX(J,N)
           SIGPYJN=-ANJ*SCRY(J,N)
           SIGPZJN=-ANJ*SCRZ(J,N)
C Write out the localized contribution
           WRITE (LFNI) IAT,JY,J,N,SIGPXJN,SIGPYJN,SIGPZJN
        END DO
      END DO
C END OF NBO/NLMO ANALYSIS
      RETURN
      END
C***********************************************************************
      FUNCTION TRCNCS(A,B,NDIM)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z)
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM)
      DATA ZERO/0.0D0/
      TMP = ZERO
      DO I = 1,NDIM
        DO J = 1,NDIM
          TMP = TMP + A(I,J)*B(J,I)
        ENDDO
      ENDDO
      TRCNCS = TMP
      RETURN
      END
C***********************************************************************
      SUBROUTINE NCSXYZ(IAT,NBAS,NOCC,SIGXL,SIGYL,SIGZL,LFNU,LFNI,
     + LFNPR,LFNISO,IWXYZ,IWCSA,IWD,IWP,IZ,R,IORDER,IWCMO)
C***********************************************************************
C 22-Jun-01  FAW  Altered LFNZ to 91 to avoid conflict with plotfile 31
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
      CHARACTER CHARAT*2
C IWCMO=1 use printing style appropriate for canonical MOs instead.
C This routine accumulates total shielding tensor so that it can be
C diagonalized by the routine NCSCSA and printed
C R and IORDER are used by NCSCSA but they need to be kept for later.
      DIMENSION SIGXL(NBAS,NBAS),SIGYL(NBAS,NBAS),SIGZL(NBAS,NBAS)
      DIMENSION SIGT(9),SIGL(9),R(3,3),SIGTL(9),SIGTNL(9),IORDER(3)
      SAVE LFNX,LFNY,LFNZ
      DATA ZERO/0.0D0/THREE/3.0D0/TWO/2.0D0/
      DATA LFNX/29/LFNY/30/LFNZ/91/
C constant to convert A. U. to ppm
C      DATA PPM/26.6256691D+00/
C For G98, the following is needed
       DATA PPM/53.2513382D+00/
C Open the cartesian scratch files
       JLAST=0
      OPEN (UNIT=LFNX,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     + FORM='UNFORMATTED')
      OPEN (UNIT=LFNY,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     + FORM='UNFORMATTED')
      OPEN (UNIT=LFNZ,STATUS='SCRATCH',ACCESS='SEQUENTIAL',
     + FORM='UNFORMATTED')
C Add the "unperturbed" and "induced" contributions; store in
C external scratch files LFNX, LFNY, LFNZ
C
C IWXYZ=1  required to print data out from this subroutine
C ITYP is internal code to this subroutine for type of output
C ITYP=1 add ind+unpert., 2= unperturbed only
C ITYP=3 induced only.
      IF (IWD.EQ.0.AND.IWP.EQ.0) ITYP=1
      IF (IWD.EQ.1.AND.IWP.EQ.0) ITYP=2
      IF (IWD.EQ.0.AND.IWP.EQ.1) ITYP=3
C this last case would be an error-- print out a total and keep going
      IF (IWD.EQ.1.AND.IWP.EQ.1) ITYP=1
C Print header
      IF (IWXYZ.GE.1) THEN
        WRITE (LFNPR,*)
        IF (ITYP.EQ.1) WRITE (LFNPR,80) CHARAT(IZ),IAT
        IF (ITYP.EQ.2) WRITE (LFNPR,90) CHARAT(IZ),IAT
        IF (ITYP.EQ.3) WRITE (LFNPR,100) CHARAT(IZ),IAT
        IF (IWCMO.EQ.1) THEN
          WRITE (LFNPR,210)
        ELSE
          WRITE (LFNPR,200)
        END IF
        WRITE (LFNPR,*)
        IF (IWCMO.EQ.1) THEN
          WRITE (LFNPR,130)
        ELSE
          WRITE (LFNPR,120)
        END IF
        WRITE (LFNPR,140)
      END IF
C Grand loop over cartesian axes IX
      DO 50 IX=1,3
C Rewind the scratch files for reading
        REWIND (LFNU)
        REWIND (LFNI)
C Zero the SIGL arrays
        DO I=1,NBAS
          DO J=1,NBAS
            SIGXL(I,J)=ZERO
            SIGYL(I,J)=ZERO
            SIGZL(I,J)=ZERO
          END DO
        END DO
        IF (ITYP.EQ.3) GO TO 20
C Read another "unperturbed" contribution
   10   READ (LFNU,END=20) IDUM,IXD,J,N,SIGDXJN,SIGDYJN,SIGDZJN
        IF (IDUM.NE.IAT)
     +    CALL NBHALT('NCSXYZ: Bad atom number, unperturbed table.')
C Canonical MO contribution N=0
        IF (N.NE.0.AND.IWCMO.EQ.1) GO TO 10
        IF (N.EQ.0.AND.IWCMO.NE.1) GO TO 10
        IF (N.EQ.0.AND.IWCMO.EQ.1) N=J
        IF (IXD.EQ.IX) THEN
          SIGXL(J,N)=SIGDXJN
          SIGYL(J,N)=SIGDYJN
          SIGZL(J,N)=SIGDZJN
          GO TO 10
        ELSE IF (IXD.LT.IX) THEN
          GO TO 10
        END IF
C Read another "induced" contribution, add to SIGL arrays
   20   IF (ITYP.EQ.2) GO TO 40
   30   READ (LFNI,END=40) IDUM,IXD,J,N,SIGPXJN,SIGPYJN,SIGPZJN
        IF (IDUM.NE.IAT)
     +    CALL NBHALT('NCSXYZ: Bad atom number, induced table.')
C Canonical MO contribution N=0
        IF (N.NE.0.AND.IWCMO.EQ.1) GO TO 30
        IF (N.EQ.0.AND.IWCMO.NE.1) GO TO 30
        IF (N.EQ.0.AND.IWCMO.EQ.1) N=J
        IF (IXD.EQ.IX) THEN
          SIGXL(J,N)=SIGXL(J,N)-SIGPXJN
          SIGYL(J,N)=SIGYL(J,N)-SIGPYJN
          SIGZL(J,N)=SIGZL(J,N)-SIGPZJN
          GO TO 30
        ELSE IF (IXD.LT.IX) THEN
          GO TO 30
        END IF
C  WRITE OUT THE TOTAL--NBO/NLMO contributions
   40   IF (IWCMO.EQ.1) GO TO 220
        DO J=1,NOCC
          DO N=1,NBAS
            IF (IX.EQ.1) THEN
              WRITE (LFNX) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
            ELSE IF (IX.EQ.2) THEN
              WRITE (LFNY) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
            ELSE IF (IX.EQ.3) THEN
              WRITE (LFNZ) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
            END IF
          END DO
        END DO
        GO TO 50
C  WRITE OUT THE TOTAL--MO contributions if needed, instead.
C
  220   DO J=1,NOCC
          N=J
          IF (IX.EQ.1) THEN
            WRITE (LFNX) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
          ELSE IF (IX.EQ.2) THEN
            WRITE (LFNY) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
          ELSE IF (IX.EQ.3) THEN
            WRITE (LFNZ) J,N,SIGXL(J,N),SIGYL(J,N),SIGZL(J,N)
          END IF
        END DO
   50 CONTINUE
C Rewind the X,Y,Z tables
      REWIND (LFNX)
      REWIND (LFNY)
      REWIND (LFNZ)
      DO I=1,9
        SIGT(I)=ZERO
        SIGTL(I)=ZERO
        SIGTNL(I)=ZERO
      END DO
C Read another set of J=>N delocalization terms
C Just accumulate total shielding
   60 READ (LFNX,END=70) JX,NX,SIGL(1),SIGL(4),SIGL(7)
      READ (LFNY,END=70) JY,NY,SIGL(2),SIGL(5),SIGL(8)
      READ (LFNZ,END=70) JZ,NZ,SIGL(3),SIGL(6),SIGL(9)
      IF (JX.NE.JY.OR.JX.NE.JZ) GO TO 230
      IF (NX.NE.NY.OR.NX.NE.NZ) GO TO 230
      J=JX
      N=NX
C Accumulate the total shielding for each tensor component
C convert to PPM
      DO I=1,9
        SIGL(I)=SIGL(I)*TWO*PPM
      END DO
      DO I=1,9
        SIGT(I)=SIGT(I)+SIGL(I)
      END DO
      IF (N.EQ.J) THEN
        DO I=1,9
          SIGTL(I)=SIGTL(I)+SIGL(I)
        END DO
      ELSE
        DO I=1,9
          SIGTNL(I)=SIGTNL(I)+SIGL(I)
        END DO
      END IF
      SIGI=(SIGL(1)+SIGL(5)+SIGL(9))/THREE
C Write out the L and NL contributions to isotropic shielding (LFNISO)
      IF (ITYP.EQ.1.AND.IWCMO.NE.1) WRITE (LFNISO) IAT,J,N,SIGI
      IF (ITYP.EQ.1.AND.IWCMO.EQ.1) WRITE (LFNISO) IAT,J,0,SIGI
C Is the contribution above the printing threshold?
      SMAX=ZERO
      DO I=1,9
        IF (ABS(SIGL(I)).GT.SMAX) SMAX=ABS(SIGL(I))
      END DO
      IF (SMAX.LT.THRNCS) GO TO 60
C Add the J=>N delocalization to the output table
C
C New diagonal "local" term (J = J)?
      IF (J.NE.JLAST.AND.J.EQ.N) THEN
        IF (IWXYZ.GE.1) WRITE (LFNPR,150) J,(SIGL(I),I=1,9)
        JLAST=J
        GO TO 60
C New off-diagonal delocalization term (J=>N)
      ELSE IF (J.EQ.JLAST.AND.N.NE.J) THEN
        IF (IWXYZ.GE.1) WRITE (LFNPR,160) N,(SIGL(I),I=1,9)
        JLAST=J
        GO TO 60
      ELSE
        GO TO 60
      END IF
C Print the total tensorial and isotropic shielding
   70 IF (IWXYZ.GE.1) THEN
        WRITE (LFNPR,110)
        IF (IWCMO.NE.1) THEN
          WRITE (LFNPR,180) (SIGTL(I),I=1,9)
          WRITE (LFNPR,190) (SIGTNL(I),I=1,9)
          WRITE (LFNPR,110)
        END IF
        WRITE (LFNPR,170) (SIGT(I),I=1,9)
        WRITE (LFNPR,110)
      END IF
C Call the routine to diagonalize and print out total tensor
      IF (IWCSA.GT.0) CALL NCSCSA (IAT,SIGT,LFNX,LFNY,LFNZ,
     +LFNPR,IWCSA,IZ,ITYP,R,IORDER,IWCMO)
      CLOSE (LFNX,STATUS='DELETE')
      CLOSE (LFNY,STATUS='DELETE')
      CLOSE (LFNZ,STATUS='DELETE')
      RETURN
C Error condition: Abort
  230 CALL NBHALT('NCSXYZ: Bad match of delocalization terms.')
   80 FORMAT (/,'Full Cartesian NMR shielding tensor (ppm) for atom ',
     + A2,'(',I2,'):')
   90 FORMAT (/,'Unperturbed  Cartesian NMR shielding ',
     + 'tensor (ppm) for atom ',A2,'(',I2,'):')
  100 FORMAT (/,'Induced   Cartesian NMR shielding ',
     + 'tensor (ppm) for atom ',A2,'(',I2,'):')
  110 FORMAT (1x,79('-'))
  120 FORMAT (3x,'NBO',/,'  L  NL      XX',6X,'YX',6X,'ZX',6X,'XY',
     +6X,'YY', 6X,'ZY',6X,'XZ',6X,'YZ',6X,'ZZ')
  130 FORMAT (' MO   ',7X,'XX',6X,'YZ',6X,'ZX',6X,'XY',6X,'YY',6X,'ZY',
     +6X,'XZ',6X,'YZ',6X,'ZZ')
c  120 FORMAT (3x,'NBO',/,'  L  NL      XX',6X,'XY',6X,'XZ',6X,'YX',
c     +6X,'YY', 6X,'YZ',6X,'ZX',6X,'ZY',6X,'ZZ')
c  130 FORMAT (' MO   ',7X,'XX',6X,'XY',6X,'XZ',6X,'YX',6X,'YY',6X,'YZ',
c     +6X,'ZX',6X,'ZY',6X,'ZZ')
  140 FORMAT (1x,79('='))
  150 FORMAT (I3,'.',4X,9(F8.2))
  160 FORMAT (4X,I3,'.',9(F8.2))
  170 FORMAT ('   Total',9(F8.2))
  180 FORMAT ('       L',9(F8.2))
  190 FORMAT ('      NL',9(F8.2))
  200 FORMAT (' Lewis (L) and non-Lewis (NL) contributions')
  210 FORMAT (' Canonical MO contributions')
      END
C***********************************************************************
      SUBROUTINE NCSCSA(IAT,SIGTOT,LFNX,LFNY,LFNZ,LFNPR,
     +IWCSA,IZ,ITYP,R,IORDER,IWCMO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
C Stuff need for making NBO labels
      PARAMETER (MAXBAS = 2000)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBLD(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      CHARACTER CHARAT*2
C Print out the table of symmetric shielding tensor components, with CSA
C for atom IAT (if IWCSA > 0)
C This subroutine also treats antisymmetric components.
C ITYP=1 is Full,  ITYP=2 is unperturbed,  ITYP=3 is induced, from
C SIGXYZ.
C SIGT WILL INITIALLY BE THE TOTAL OF THE FULL TENSOR.
C The contents of SIGT will be destroyed, to be replaced by the
C diagonalized tensor components.
C THE DIAGONAL ELEMENTS will be the diagonal components of the tensor.
      DIMENSION SIGTOT(9), SIGT(3,3), SIG(3,3), SCR(3)
C SIGAN will keep the full  tensor temporarily for printing out
C of the antisymmetric components later.
C R is the direction cosine matrix.
      DIMENSION SIGAN(3,3)
C The sum of the lewis L and non--Lewis NL contributions
      DIMENSION SIGTL(3), SIGTNL(3)
C EIGENVALUE and EIGENVECTOR OF PRINCIPAL AXIS SYSTEM
C EIGENVALUES ARE PRINCIPAL COMPONENTS OF THE TENSOR,
C EIGENVECTORS ARE THE DIRECTION COSINES WHICH TRANSFORM THE
C TENSOR FROM THE MOLECULAR FRAME TO THE PRINCIPAL AXIS FRAME
      DIMENSION SIGVAL(3), R(3,3)
C To make the program "bullet proof" from the ordering of the
C eigenvalues of  a particular
C diagonalizer, array IORDER will keep track of the ordering of the magn
C of the principal components of the shielding tensor.  IORDER(3) is
C index of component of greatest positive magnitude.
      DIMENSION IORDER(3)
      DATA ZERO/0.0D0/THREE/3.0D0/TWO/2.0D0/
      DATA HALF/0.5D0/
C constant to convert A. U. to ppm
C      DATA PPM/26.6256691D+00/
C For G98, the following is needed
      DATA PPM/53.2513382D+00/
      SIGT(1,1)=SIGTOT(1)
      SIGT(1,2)=SIGTOT(2)
      SIGT(1,3)=SIGTOT(3)
      SIGT(2,1)=SIGTOT(4)
      SIGT(2,2)=SIGTOT(5)
      SIGT(2,3)=SIGTOT(6)
      SIGT(3,1)=SIGTOT(7)
      SIGT(3,2)=SIGTOT(8)
      SIGT(3,3)=SIGTOT(9)
C CHECK FOR ANTISYMMETRY
C IWANTI is I Want ANTI-symmetric part of tensor printed
      IWANTI=0
C TURN ON PRINTING OF ANTISYMMETRIC TENSOR IF NEEDED
C IWANTI is I want  anti-symmetric part of tensor printed
      DIFF=HALF*(SIGT(1,2)-SIGT(2,1))
      IF (ABS(DIFF).GT.THRNCS) IWANTI=1
      DIFF=HALF*(SIGT(2,3)-SIGT(3,2))
C TURN ON PRINTING OF ANTISYMMETRIC TENSOR IF NEEDED
      IF (ABS(DIFF).GT.THRNCS) IWANTI=1
      DIFF=HALF*(SIGT(3,1)-SIGT(1,3))
C TURN ON PRINTING OF ANTISYMMETRIC TENSOR IF NEEDED
      IF (ABS(DIFF).GT.THRNCS) IWANTI=1
C Symmetrize the tensor (only if we have the full tensor, not
C if we have only the unperturbed or induced part)
      IF (ITYP.NE.1) GO TO 10
      SIGT(1,2)=HALF*(SIGT(1,2)+SIGT(2,1))
      SIGT(2,1)=SIGT(1,2)
      SIGT(2,3)=HALF*(SIGT(2,3)+SIGT(3,2))
      SIGT(3,2)=SIGT(2,3)
      SIGT(1,3)=HALF*(SIGT(1,3)+SIGT(3,1))
      SIGT(3,1)=SIGT(1,3)
C NOW, DIAGONALIZE  SYMMETRIC PART OF TENSOR
      CALL NBJACOBI (3,SIGT,SIGVAL,R,3,3,0)
C REZERO TOTAL
C Figure out order of eigenvalues
      IF (SIGVAL(1).GE.SIGVAL(2)) THEN
        IORDER(3)=1
        IORDER(2)=2
      ELSE
        IORDER(3)=2
        IORDER(2)=1
      END IF
      IF (SIGVAL(3).GE.SIGVAL(IORDER(3))) THEN
        IORDER(1)=IORDER(2)
        IORDER(2)=IORDER(3)
        IORDER(3)=3
      ELSE IF (SIGVAL(3).GE.SIGVAL(IORDER(2))) THEN
        IORDER(1)=IORDER(2)
        IORDER(2)=3
      ELSE
        IORDER(1)=3
      END IF
C IF NCSCSA was called with ityp=2,3, R was already calculated and
C order of eigenvalues was saved in iorder.
   10 I3=IORDER(3)
      I2=IORDER(2)
      I1=IORDER(1)
C ZERO VARIABLES FOR ACCUMULATION OF TOTALS
      DO I=1,3
        SIGT(I,I)=ZERO
        SIGTL(I)=ZERO
        SIGTNL(I)=ZERO
      END DO
      SIG12L=ZERO
      SIG12NL=ZERO
      SIG13L=ZERO
      SIG13NL=ZERO
      SIG23L=ZERO
      SIG23NL=ZERO
      JLAST=0
      IF (IWCSA.GT.0) THEN
        IF (ITYP.EQ.1) WRITE (LFNPR,50) CHARAT(IZ),IAT
        IF (ITYP.EQ.2) WRITE (LFNPR,60) CHARAT(IZ),IAT
        IF (ITYP.EQ.3) WRITE (LFNPR,70) CHARAT(IZ),IAT
        IF (IWCMO.EQ.1) WRITE (LFNPR,370)
        IF (IWANTI.EQ.1) THEN
          WRITE (LFNPR,360)
          WRITE (LFNPR,*)
          WRITE (LFNPR,80)
          IF (IWCMO.EQ.1) THEN
            WRITE (LFNPR,140)
          ELSE
            WRITE (LFNPR,120)
          END IF
          WRITE (LFNPR,160)
        ELSE
          WRITE (LFNPR,*)
          IF (IWCMO.EQ.1) THEN
            WRITE (LFNPR,130)
          ELSE
            WRITE (LFNPR,110)
          END IF
          WRITE (LFNPR,150)
        END IF
      END IF
      REWIND (LFNX)
      REWIND (LFNY)
      REWIND (LFNZ)
C Read another set of J=>N delocalization terms
C Note that the cartesian components were stored in
C transpose order.
   20 READ (LFNX,END=40) JX,NX,SIG(1,1),SIG(2,1),SIG(3,1)
      READ (LFNY,END=40) JY,NY,SIG(1,2),SIG(2,2),SIG(3,2)
      READ (LFNZ,END=40) JZ,NZ,SIG(1,3),SIG(2,3),SIG(3,3)
      IF (JX.NE.JY.OR.JX.NE.JZ) GO TO 380
      IF (NX.NE.NY.OR.NX.NE.NZ) GO TO 380
      J=JX
      N=NX
C Accumulate the total shielding for each tensor component
C convert to PPM
      IF (IWANTI.EQ.1) THEN
        SIGAN(1,1)=ZERO
        SIGAN(2,2)=ZERO
        SIGAN(3,3)=ZERO
        SIGAN(1,2)=HALF*(SIG(1,2)-SIG(2,1))
        SIGAN(2,1)=-SIGAN(1,2)
        SIGAN(2,3)=HALF*(SIG(2,3)-SIG(3,2))
        SIGAN(3,2)=-SIGAN(2,3)
        SIGAN(1,3)=HALF*(SIG(1,3)-SIG(3,1))
        SIGAN(3,1)=-SIGAN(1,3)
C TRANSFORM TO SAME FRAME AS 11,22,33
        CALL SIMTRS (SIGAN,R,SCR,3,3)
      END IF
C  SYMMETRIZE THE TENSOR
      SIG(1,2)=HALF*(SIG(1,2)+SIG(2,1))
      SIG(2,1)=SIG(1,2)
      SIG(2,3)=HALF*(SIG(2,3)+SIG(3,2))
      SIG(3,2)=SIG(2,3)
      SIG(1,3)=HALF*(SIG(1,3)+SIG(3,1))
      SIG(3,1)=SIG(1,3)
C DO SIMILARITY TRANSFORM
      CALL SIMTRS (SIG,R,SCR,3,3)
C only worry about diagonal components
      DO I=1,3
        SIG(I,I)=SIG(I,I)*TWO*PPM
      END DO
      DO I=1,3
        DO K=1,3
          SIGAN(I,K)=SIGAN(I,K)*TWO*PPM
        END DO
      END DO
      DO I=1,3
        SIGT(I,I)=SIGT(I,I)+SIG(I,I)
      END DO
C Is the contribution above the printing threshold?
      SMAX=ZERO
      DO I=1,3
        IF (ABS(SIG(I,I)).GT.SMAX) SMAX=ABS(SIG(I,I))
      END DO
      IF (J.EQ.N) THEN
C ACCUMULATE ANTI SYMMETRIC COMPONENTS
        SIG12L=SIG12L+SIGAN(I1,I2)
        SIG13L=SIG13L+SIGAN(I1,I3)
        SIG23L=SIG23L+SIGAN(I2,I3)
C ACCUMULATE SYMMETRIC COMPONENTS
        DO I=1,3
          SIGTL(I)=SIGTL(I)+SIG(I,I)
        END DO
      ELSE
C ACCUMULATE ANTI-SYMMETRIC, non-LEWIS components
        SIG12NL=SIG12NL+SIGAN(I1,I2)
        SIG13NL=SIG13NL+SIGAN(I1,I3)
        SIG23NL=SIG23NL+SIGAN(I2,I3)
C Accumulate symmetric non-lewis components
        SIGTNL(1)=SIGTNL(1)+SIG(1,1)
        SIGTNL(2)=SIGTNL(2)+SIG(2,2)
        SIGTNL(3)=SIGTNL(3)+SIG(3,3)
      END IF
      IF (SMAX.LT.THRNCS) GO TO 20
C Add the J=>N delocalization to the output table
      CSA=SIG(I3,I3)-(SIG(I2,I2)+SIG(I1,I1))*HALF
      SISO=(SIG(1,1)+SIG(2,2)+SIG(3,3))/THREE
C jump to canonical MO format ?
      IF (IWCMO.EQ.1) GO TO 30
C New diagonal "local" term (J = J)?
      IF (J.NE.JLAST.AND.J.EQ.N) THEN
        IF (IWANTI.EQ.1) THEN
          IF (IWCSA.GT.0) WRITE (LFNPR,180) J,(LBL(I,N,4),I=1,10),
     +     (SIG(IORDER(I),IORDER(I)),I=1,3),SIGAN(I1,I2),SIGAN(I1,I3),
     +     SIGAN(I2,I3),CSA,SISO
        ELSE
          IF (IWCSA.GT.0) WRITE (LFNPR,170) J,(LBL(I,N,4),I=1,10),
     +     (SIG(IORDER(I),IORDER(I)),I=1,3),CSA,SISO
        END IF
        JLAST=J
        GO TO 20
C New off-diagonal delocalization term (J=>N)
      ELSE IF (J.EQ.JLAST.AND.N.NE.J) THEN
        IF (IWANTI.EQ.1) THEN
          IF (IWCSA.GT.0) WRITE (LFNPR,220) N,(LABEL(N,I),I=1,2),
     +     (SIG(IORDER(I),IORDER(I)),I=1,3),SIGAN(I1,I2),SIGAN(I1,I3),
     +     SIGAN(I2,I3),CSA,SISO
        ELSE
          IF (IWCSA.GT.0) WRITE (LFNPR,210) N,(LBL(I,N,4),I=1,10),
     +     (SIG(IORDER(I),IORDER(I)),I=1,3),CSA,SISO
        END IF
        JLAST=J
        GO TO 20
      ELSE
        GO TO 20
      END IF
C END OF NBO FORMAT STYLE
C CANONICAL MO format style
   30 IF (J.NE.JLAST.AND.J.EQ.N) THEN
        IF (IWANTI.EQ.1) THEN
          IF (IWCSA.GT.0) WRITE (LFNPR,200) J,(SIG(IORDER(I),IORDER(I)),
     +     I=1,3),SIGAN(I1,I2),SIGAN(I1,I3),SIGAN(I2,I3),CSA,SISO
        ELSE
          IF (IWCSA.GT.0) WRITE (LFNPR,190) J,(SIG(IORDER(I),IORDER(I)),
     +     I=1,3),CSA,SISO
        END IF
        JLAST=J
        GO TO 20
      END IF
C Print the total
C Do we want antisymmetric components?
   40 IF (IWANTI.EQ.1) THEN
        IF (IWCSA.GT.0) THEN
          WRITE (LFNPR,100)
          IF (IWCMO.NE.1) THEN
C PRINT  TOTAL LEWIS
            CSA=SIGTL(I3)-(SIGTL(I2)+SIGTL(I1))*HALF
            SISO=(SIGTL(1)+SIGTL(2)+SIGTL(3))/THREE
            WRITE (LFNPR,280) (SIGTL(IORDER(I)),I=1,3),SIG12L,SIG13L,
     +       SIG23L,CSA,SISO
C PRINT TOTAL NON--LEWIS
            CSA=SIGTNL(I3)-(SIGTNL(I2)+SIGTNL(I1))*HALF
            SISO=(SIGTNL(1)+SIGTNL(2)+SIGTNL(3))/THREE
            WRITE (LFNPR,270) (SIGTNL(IORDER(I)),I=1,3),SIG12NL,SIG13NL,
     +       SIG23NL,CSA,SISO
C  PRINT TOTAL
            WRITE (LFNPR,100)
          END IF
          CSA=SIGT(I3,I3)-(SIGT(I2,I2)+SIGT(I1,I1))*HALF
          SISO=(SIGT(1,1)+SIGT(2,2)+SIGT(3,3))/THREE
          WRITE (LFNPR,290) (SIGT(IORDER(I),IORDER(I)),I=1,3),SIG12L+
     +     SIG12NL,SIG13L+SIG13NL,SIG23L+SIG23NL,CSA,SISO
          WRITE (LFNPR,100)
        END IF
C NO ANTISYMMETRIC COMPONENTS
      ELSE
        IF (IWCSA.GT.0) THEN
          WRITE (LFNPR,90)
          IF (IWCMO.NE.1) THEN
C NBO/NLMO TOTAL
C PRINT  TOTAL LEWIS
            CSA=SIGTL(I3)-(SIGTL(I2)+SIGTL(I1))*HALF
            SISO=(SIGTL(1)+SIGTL(2)+SIGTL(3))/THREE
            WRITE (LFNPR,250) (SIGTL(IORDER(I)),I=1,3),CSA,SISO
C PRINT TOTAL NON--LEWIS
            CSA=SIGTNL(I3)-(SIGTNL(I2)+SIGTNL(I1))*HALF
            SISO=(SIGTNL(1)+SIGTNL(2)+SIGTNL(3))/THREE
            WRITE (LFNPR,240) (SIGTNL(IORDER(I)),I=1,3),CSA,SISO
C  PRINT TOTAL
            WRITE (LFNPR,90)
            CSA=SIGT(I3,I3)-(SIGT(I2,I2)+SIGT(I1,I1))*HALF
            SISO=(SIGT(1,1)+SIGT(2,2)+SIGT(3,3))/THREE
            WRITE (LFNPR,230) (SIGT(IORDER(I),IORDER(I)),I=1,3),CSA,
     +       SISO
            WRITE (LFNPR,90)
          ELSE
C Canonical MO total
            CSA=SIGT(I3,I3)-(SIGT(I2,I2)+SIGT(I1,I1))*HALF
            SISO=(SIGT(1,1)+SIGT(2,2)+SIGT(3,3))/THREE
            WRITE (LFNPR,260) (SIGT(IORDER(I),IORDER(I)),I=1,3),CSA,
     +       SISO
            WRITE (LFNPR,90)
          END IF
        END IF
      END IF
C PRINT DIRECTION COSINES
      IF (ITYP.EQ.1.AND.IWCMO.EQ.0) THEN
        WRITE (LFNPR,300) CHARAT(IZ),IAT
        WRITE (LFNPR,350)
        WRITE (LFNPR,340)
        WRITE (LFNPR,310) (R(1,IORDER(I)),I=1,3)
        WRITE (LFNPR,320) (R(2,IORDER(I)),I=1,3)
        WRITE (LFNPR,330) (R(3,IORDER(I)),I=1,3)
      END IF
      RETURN
C Error condition: Abort
  380 CALL NBHALT('NCSCSA: Bad match of delocalization terms.')
   50 FORMAT (/,'Principal components of the tensor (ppm) for atom ',A2,
     +'(',I2,'):')
   60 FORMAT (/,'Principal components of the unperturbed',
     +' tensor (ppm) for atom ',A2,'(',I2,'):')
   70 FORMAT (/,'Principal components of the induced tensor (ppm)',
     +' for atom ',A2,'(',I2,'):')
   80 FORMAT (21(' '),'Principal Components',4(' '),'Antisymmetric Part'
     +,17(' '))
   90 FORMAT (1x,60('-'))
  100 FORMAT (1x,78('-'))
  110 FORMAT (7X,'NBO',16X,'11',6X,'22',6X,'33',5X,'CSA',5X,'ISO')
  120 FORMAT (7X,'NBO',10X,'11',6X,'22',6X,'33',6X,'12',6X,'13',6X,'23',
     +5X,'CSA',5X,'ISO')
  130 FORMAT (7X,'MO',17X,'11',6X,'22',6X,'33',5X,'CSA',5X,'ISO')
  140 FORMAT (7X,'MO',11X,'11',6X,'22',6X,'33',6X,'12',6X,'13',6X,'23',
     +5X,'CSA',5X,'ISO')
  150 FORMAT (1x,60('='))
  160 FORMAT (1x,78('='))
  170 FORMAT (I3,'.',10A1,7X,5(F8.2))
  180 FORMAT (I3,'.',10A1,1X,8(F8.2))
  190 FORMAT (I3,'.',10X,7X,5(F8.2))
  200 FORMAT (I3,'.',10X,1X,8(F8.2))
  210 FORMAT (4X,'NL ',I3,'.',10A1,5(F8.2))
  220 FORMAT (4X,'NL ',I3,'.',' ',A2,A1,8(F8.2))
  230 FORMAT (14X,'  Total',5(F8.2))
  240 FORMAT (10X,'  non-Lewis',5(F8.2))
  250 FORMAT (14X,'  Lewis',5(F8.2))
  260 FORMAT (14X,'  Total',5(F8.2))
  270 FORMAT (5X,' non-Lewis',8(F8.2))
  280 FORMAT (9X,' Lewis',8(F8.2))
  290 FORMAT (9X,' Total',8(F8.2))
  300 FORMAT (/,'Cartesian XYZ to principal shielding axes ',
     + 'for atom ',A2, 1x,'(',I2,'):',/)
  310 FORMAT (' X',1X,3(F9.6,1X))
  320 FORMAT (' Y',1X,3(F9.6,1X))
  330 FORMAT (' Z',1X,3(F9.6,1X))
  340 FORMAT (1x,31('-'))
  350 FORMAT (8X,'1',9X,'2',9X,'3')
  360 FORMAT (' This tensor is non-symmetric. The antisymmetric part',
     + ' will be printed.')
  370 FORMAT (' Canonical MO contributions')
      END
C***********************************************************************
      SUBROUTINE NCSISO(NOCC,NATM,NUC,SIGISO,LFNISO,LFNPR,IWCMO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER DASH*7,CHARAT*2
C Write out the table of localized isotropic chemical shieldings
      DIMENSION SIGISO(NOCC,2,NATM), SIGT(9), SIGL(9), SIGNL(9)
      DIMENSION NUC(99)
C Stuff need for making NBO labels
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      DATA ZERO/0.0D0/
      DATA DASH/'-------'/
      REWIND (LFNISO)
C Count how many atoms we have shielding values for
C put the atoms we have shieldings for in a linear list
C Zero the SIGISO array
      DO I=1,NOCC
        DO IAT=1,NATM
          SIGISO(I,1,IAT)=ZERO
          SIGISO(I,2,IAT)=ZERO
        END DO
      END DO
   10 READ (LFNISO,END=20) IAT,J,N,SIGI
      IF (IWCMO.EQ.1) THEN
        IF (N.EQ.0) SIGISO(J,1,IAT)=SIGI
      ELSE
        IF (J.EQ.N) THEN
          SIGISO(J,1,IAT)=SIGI
        ELSE
          IF (N.NE.0) SIGISO(J,2,IAT)=SIGISO(J,2,IAT)+SIGI
        END IF
      END IF
      GO TO 10
C Print the table of SIGISO values
   20 NCOL=7
      IF (IWCMO.EQ.1) THEN
        WRITE (LFNPR,50)
      ELSE
        WRITE (LFNPR,40)
      END IF
      IAT2=0
   30 IAT1=IAT2+1
      IAT2=IAT1+NCOL-1
      IF (IAT2.GT.NATM) IAT2=NATM
      NPR=IAT2-IAT1+1
      WRITE (LFNPR,*)
      IF (IWCMO.EQ.1) THEN
        WRITE (LFNPR,70) (CHARAT(IATNO(NUC(IAT))),NUC(IAT),IAT=IAT1,
     +   IAT2)
      ELSE
        WRITE (LFNPR,60) (CHARAT(IATNO(NUC(IAT))),NUC(IAT),IAT=IAT1,
     +   IAT2)
      END IF
      WRITE (LFNPR,110) (DASH,ID=1,NPR)
C Zero the accumulation vector
      DO I=1,9
        SIGT(I)=ZERO
        SIGL(I)=ZERO
        SIGNL(I)=ZERO
      END DO
      IF (IWCMO.EQ.1) THEN
        DO J=1,NOCC
C Write out the contributions from MO J
          WRITE (LFNPR,90) J,(SIGISO(J,1,NUC(IAT)),IAT=IAT1,IAT2)
C Add contribution to total isotropic shielding
          IP=0
          DO IAT=IAT1,IAT2
            IP=IP+1
            SIGT(IP)=SIGT(IP)+SIGISO(J,1,NUC(IAT))
          END DO
        END DO
      ELSE
        DO J=1,NOCC
C Write out the contributions from NBO J
          WRITE (LFNPR,80) J,(LBL(I,J,4),I=1,10),(SIGISO(J,1,NUC(IAT)),
     +     IAT=IAT1,IAT2)
          WRITE (LFNPR,100) (SIGISO(J,2,NUC(IAT)),IAT=IAT1,IAT2)
C Add contribution to total isotropic shielding
          IP=0
          DO IAT=IAT1,IAT2
            IP=IP+1
            SIGT(IP)=SIGT(IP)+SIGISO(J,1,NUC(IAT))+SIGISO(J,2,NUC(IAT))
            SIGNL(IP)=SIGNL(IP)+SIGISO(J,2,NUC(IAT))
            SIGL(IP)=SIGL(IP)+SIGISO(J,1,NUC(IAT))
          END DO
        END DO
      END IF
      IF (IWCMO.NE.1) THEN
        WRITE (LFNPR,110) (DASH,ID=1,NPR)
        WRITE (LFNPR,130) (SIGL(IP),IP=1,NPR)
        WRITE (LFNPR,140) (SIGNL(IP),IP=1,NPR)
      END IF
      WRITE (LFNPR,110) (DASH,ID=1,NPR)
      WRITE (LFNPR,120) (SIGT(IP),IP=1,NPR)
      IF (IAT2.LT.NATM) GO TO 30
      WRITE (LFNPR,*)
      RETURN
   40 FORMAT (/,'Summary of isotropic NMR chemical shielding ',/,
     + 'Total Lewis (L) and non-Lewis (NL) contributions: (ppm)')
   50 FORMAT (/,'Summary of isotropic NMR chemical shielding ',/,
     + 'Canonical MO contributions: (ppm)')
   60 FORMAT (5X,'NBO',7X,7(2X,A2,'(',I2,')'))
   70 FORMAT (5X,'MO ',7X,7(2X,A2,'(',I2,')'))
   80 FORMAT (I3,'.',10(A1),'L ',7F8.2)
   90 FORMAT (I3,'.',10X,'  ',7F8.2)
  100 FORMAT (14X,'NL',7F8.2)
  110 FORMAT (1x,12('-'),4X,7(A7,1X))
  120 FORMAT (8X,'Total',3X,7F8.2)
  130 FORMAT (8X,'Lewis',3X,7F8.2)
  140 FORMAT (4X,'non-Lewis',3X,7F8.2)
      END
C***********************************************************************
C
C  FIX DENSITY MATRIX (FIXDM) ROUTINES:
C
C      SUBROUTINE FIXDM(DM,T,V,S,EV,S1,S2,MAP,NDIM,NBAS,NNAO,LFNPR,IFIX)
C      SUBROUTINE AXWBX(NDIM,N,A,B,W,IORTH,IWMATX,X,V1,V2,IERR,LFN)
C      SUBROUTINE REBAKA(NDIM,N,B,DL,Z)
C      SUBROUTINE REDUC1(NDIM,N,A,B,DL,IERR,LFNPR)
C      SUBROUTINE TQLRAT(N,D,EZ,IERR)
C      SUBROUTINE TQL2(NDIM,N,D,E,Z,IERR)
C      SUBROUTINE TRED1(NDIM,N,A,D,E,E2)
C      SUBROUTINE TRED2(NDIM,N,A,D,E,Z)
C
C***********************************************************************
      SUBROUTINE FIXDM(DM,T,V,S,EV,S1,S2,MAP,NDIM,NBAS,NNAO,LFNPR,IFIX)
C***********************************************************************
C  1-Jul-01  EDG  Modified to work with linear dependencies
C 27-Nov-00  FAW  Add FIXDM keyword [IWFIXDM]
C-----------------------------------------------------------------------
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
C
C  Fix input density matrix (FIXDM keyword) to remove unphysical
C  populations (negative or Pauli-violating) and make total number
C  of electrons an integer, "close" to input density matrix.
C
C  Input:
C      DM(NBAS,NBAS) : AO density matrix
C      T(NBAS,NNAO)  : AO-NAO transformation matrix
C  Output:
C      DM(NBAS,NBAS) : "fixed" AO density
C      IFIX          : 1(fixed)/0(not fixed)
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
C
      DIMENSION DM(NDIM,NDIM),T(NDIM,NDIM),V(NDIM,NDIM)
      DIMENSION S(NDIM,NDIM),EV(NDIM),S1(NDIM),S2(NDIM),MAP(NDIM)
C
      DATA ZERO/0.0D0/ONE/1.0D0/TWO/2.0D0/FIVE/5.0D0/THRSH/1.D-6/
      SAVE ZERO,ONE,TWO,FIVE,THRSH
C
C  Calculate the NAO density:
C
      CALL SIMTR1(DM,T,S1,NDIM,NBAS,NBAS,NNAO)
C
C  Calculate eigenvalues of the NAO density:
C
      CALL AXWBX(NDIM,NNAO,DM,S,EV,1,1,V,S1,S2,IERR,LFNPR)
      CALL RANK(EV,NNAO,NDIM,MAP)
      DO I = 1,NNAO
        CALL COPY(V(1,MAP(I)),DM(1,I),NDIM,NNAO,1)
      END DO
      CALL COPY(DM,V,NDIM,NNAO,NNAO)
C
C  Adjust occupancies that are negative (to zero) or greater than
C  occmax (to occmax):
C
      IF(ALPHA.OR.BETA) THEN
        OCCMAX = ONE
      ELSE
        OCCMAX = TWO
      ENDIF
      NFIXM = 0
      NFIXP = 0
      TRACE = ZERO
      FULL = ZERO
      DO I = 1,NNAO
        IF(EV(I).LT.-THRSH) THEN
          EV(I) = ZERO
          NFIXM = NFIXM + 1
        ELSE IF(EV(I).LT.ZERO) THEN
          EV(I) = ZERO
        ELSE IF(EV(I).GT.(OCCMAX+THRSH)) THEN
          EV(I) = OCCMAX
          FULL  = FULL + EV(I)
          NFIXP = NFIXP + 1
        ELSE IF(EV(I).GT.OCCMAX) THEN
          EV(I) = OCCMAX
          FULL  = FULL + EV(I)
        END IF
        TRACE = TRACE + EV(I)
      END DO
C
C  Find nearest integer number of electrons:
C
      INTTR  = NINT(TRACE)
      WANTTR = DFLOAT(INTTR)
      DEV = TRACE - WANTTR
C
C  Print diagnostics:
C
      IFIX = 0
      IF(NFIXM.GT.0.OR.NFIXP.GT.0.OR.ABS(DEV).GT.THRSH*FIVE) THEN
        IFIX = 1
        IF(NNAO.EQ.NBAS) WRITE(LFNPR,*)
        IF(NFIXM.GT.0.OR.NFIXP.GT.0) WRITE(LFNPR,1300) NFIXM,NFIXP
        IF(ABS(DEV).GT.THRSH*FIVE) WRITE(LFNPR,1400) DEV
C
C  Rescale remaining occupancies to integer trace (nearest integer):
C
        TRPART = TRACE - FULL
        WAPART = WANTTR - FULL
        SC = WAPART/TRPART
        DO I = 1,NNAO
          IF(EV(I).LT.OCCMAX.AND.EV(I).GT.ZERO) THEN
            EVTMP = SC*EV(I)
            EV(I) = MIN(EVTMP,OCCMAX)
            WAPART = WAPART - EV(I)
            TRPART = TRPART - EV(I)
            IF(WAPART.LE.ZERO) THEN
              SC = ZERO
            ELSE IF(TRPART.GT.THRSH.AND.EV(I).NE.EVTMP) THEN
              SC = WAPART/TRPART
            END IF
          END IF
        END DO
C
C  Multiply T*V, store transposed result in S:
C
        CALL COPY(T,S,NDIM,NBAS,NNAO)
        CALL MATML1(S,V,S1,NDIM,NDIM,NDIM,NBAS,NNAO,NNAO)
C
C  Re-compose DM = Sum(k) EV(k) |T*v(k)> <T*v(k)| with modified EV(k)'s:
C
        DO I = 1,NBAS
          DO J = 1,I
            TMP = ZERO
            DO K = 1,NNAO
              TMP = TMP + EV(K)*S(I,K)*S(J,K)
            END DO
            DM(I,J) = TMP
            DM(J,I) = TMP
          END DO
        END DO
      END IF
      RETURN
C
 1300 FORMAT(1X,'FIXDM:',I4,' negative populations and',I4,' Pauli-',
     + 'violating populations; fixed')
 1400 FORMAT(1X,'FIXDM: density matrix trace differs from integer ',
     + 'value by',F9.5,'e; fixed')
      END
C***********************************************************************
      SUBROUTINE AXWBX(NDIM,N,A,B,W,IORTH,IWMATX,X,V1,V2,IERR,LFN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),W(NDIM),X(NDIM,NDIM),
     + V1(NDIM),V2(NDIM)
C
C  Find eigenvalues and eigenvectors of real symmetric generalized
C  eigenvalue problem (A - W*B)X = 0
C  [Martin & Wilkinson, Num. Math. 11, 99 (1968)]
C
C   NDIM = declared row dimension of A, B, X
C      N = actual dimension
C      A = real symmetrix matrix
C      B = positive definite matrix (not used for IORTH.EQ.1)
C      X = eigenvector matrix
C      W = ordered eigenvalues (lowest first)
C  IORTH = 0 (if basis functions are nonorthogonal) or 1
C IWMATX = 1 (if eigenvectors wanted) or 0 (eigenvalues only)
C   IERR = error code = 7*N + 1 if B is not positive definite
C                     = N if convergence failure on eigenvalue N
C  V1,V2 = scratch vectors
C
      IF(IORTH.NE.1) THEN
        CALL REDUC1(NDIM,N,A,B,V2,IERR,LFN)
        IF(IERR.NE.0) RETURN
        IF(IWMATX.EQ.0) THEN
          CALL TRED1(NDIM,N,A,W,V1,V2)
          CALL TQLRAT(N,W,V2,IERR)
        ELSE
          CALL TRED2(NDIM,N,A,W,V1,X)
          CALL NBTQL2(NDIM,N,W,V1,X,IERR)
          IF(IERR.NE.0) RETURN
          CALL REBAKA(NDIM,N,B,V2,X)
        ENDIF
      ELSE
        IF(IWMATX.EQ.0) THEN
          CALL TRED1(NDIM,N,A,W,V1,V2)
          CALL TQLRAT(N,W,V2,IERR)
        ELSE
          CALL TRED2(NDIM,N,A,W,V1,X)
          CALL NBTQL2(NDIM,N,W,V1,X,IERR)
        ENDIF
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE REBAKA(NDIM,N,B,DL,Z)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION B(NDIM,NDIM),DL(NDIM),Z(NDIM,NDIM)
C
C  Called by AXWBX
C  [Martin & Wilkinson, Num. Math. 11, 99 (1968)]
C  Form generalized eigenvectors of (A - W*B)Z = 0 by back-transforming
C  reduced ordinary eigenvectors determined by REDUC1.
C
      DO J = 1, N
        DO II = 1, N
          I = N + 1 - II
          I1 = I + 1
          TMP = Z(I,J)
          IF(I.NE.N) THEN
            DO K = I1,N
              TMP = TMP - B(K,I)*Z(K,J)
            ENDDO
          ENDIF
          Z(I,J) = TMP/DL(I)
        ENDDO
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE REDUC1(NDIM,N,A,B,DL,IERR,LFNPR)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION A(NDIM,NDIM),B(NDIM,NDIM),DL(NDIM)
      DATA ZERO/0.0D0/
C
C  Called by AXWBX
C  [Martin & Wilkinson, Num. Math. 11, 99 (1968)]
C
C  Reduce the generalized symmetric eigenvalue problem (A - W*B)X = 0
C  (B = positive definite) to standard form with Cholesky factorization,
C  storing reduced result in lower triangle of A.
C
C  If Cholesky factor L is available, set N negative on input
C  and provide lower triangle of L in lower triangle of B,
C  diagonal elements in DL.
C
C  IERR = 7*N + 1 if B is not positive definite.
C
      IERR = 0
      NN = IABS(N)
      IF (N .LT. 0) GO TO 10
C Form Cholesky factor L in arrays B, DL
      DO I = 1,N
        I1 = I - 1
        DO J = 1,N
          TMP = B(I,J)
          IF(I.NE.1) THEN
            DO K = 1,I1
              TMP = TMP - B(I,K)*B(J,K)
            ENDDO
          ENDIF
          IF(J.EQ.I) THEN
            IF(TMP.LE.ZERO) THEN
C B is not positive definite!
              WRITE(LFNPR,1000)
              IERR = 7*N + 1
              RETURN
            ENDIF
            RTMP = SQRT(TMP)
            DL(I) = RTMP
          ELSE
            B(J,I) = TMP/RTMP
          ENDIF
        ENDDO
      ENDDO
C Form transpose of INV(L)*A (upper triangle) in lower triangle of A
   10 DO I = 1,NN
        I1 = I - 1
        Y = DL(I)
        DO J = I,NN
          TMP = A(I,J)
          IF(I.NE.1) THEN
            DO K = 1,I1
              TMP = TMP - B(I,K)*A(J,K)
            ENDDO
          ENDIF
          A(J,I) = TMP/Y
        ENDDO
      ENDDO
C Pre-multiply by INV(L)
      DO J = 1,NN
        J1 = J - 1
        DO I = J,NN
          TMP = A(I,J)
          IF(I.NE.J) THEN
            I1 = I - 1
            DO K = J,I1
              TMP = TMP - A(K,J)*B(I,K)
            ENDDO
          ENDIF
          IF(J.NE.1) THEN
            DO K = 1,J1
              TMP = TMP - A(J,K)*B(I,K)
            ENDDO
          ENDIF
          A(I,J) = TMP/DL(I)
        ENDDO
      ENDDO
      RETURN
 1000 FORMAT(1X,'REDUC1: Overlap matrix is not positive definite')
      END
C***********************************************************************
      SUBROUTINE TQLRAT(N,D,E2,IERR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION D(N),E2(N)
      SAVE ZERO,ONE,TWO
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C
C  Called by AXWBX
C  [Martin et al., Num. Math. 11, 181 (1968)]
C
C  Rational QL method for eigenvalues and eigenvectors of symmetric
C  tridiagonal matrix (diagonal elements in D, squared subdiagonal
C  elements in E2).
C
C  On return, D contains eigenvalues in ascending order (if IERR = 0).
C
C  IERR = J if J-th eigenvalue failed to converge in 30 iterations.
C
C  RMACHEP is a machine-dependent parameter specifying the relative
C  precision of floating point arithmetic.
C
      RMACHEP = TWO**(-52)
      IERR = 0
      IF (N.EQ.1) RETURN
      DO I = 2,N
        E2(I-1) = E2(I)
      ENDDO
      F = ZERO
      B = ZERO
      E2(N) = ZERO
      DO L = 1,N
        J = 0
        H = RMACHEP*(ABS(D(L)) + SQRT(E2(L)))
        IF(B.LE.H) THEN
          B = H
          C = B*B
        ENDIF
        DO M = L, N
          IF (E2(M) .LE. C) GO TO 10
        ENDDO
   10   IF (M .EQ. L) GO TO 30
   20   IF (J .EQ. 30) GO TO 1000
        J = J + 1
        L1 = L + 1
        S = SQRT(E2(L))
        G = D(L)
        P = (D(L1) - G) / (TWO*S)
        R = SQRT(P*P+ONE)
        D(L) = S / (P + SIGN(R,P))
        H = G - D(L)
        DO I = L1, N
          D(I) = D(I) - H
        ENDDO
        F = F + H
C  Rational QL transformation
        G = D(M)
        IF (G .EQ. ZERO) G = B
        H = G
        S = ZERO
        MML = M - L
        DO II = 1, MML
          I = M - II
          P = G*H
          R = P + E2(I)
          E2(I+1) = S*R
          S = E2(I)/R
          D(I+1) = H + S*(H + D(I))
          G = D(I) - E2(I)/G
          IF (G.EQ.ZERO) G = B
          H = G*P/R
        ENDDO
        E2(L) = S*G
        D(L) = H
C  Underflow guard in convergence test
        IF((H.NE.ZERO).OR.(ABS(E2(L)).GT.ABS(C/H))) THEN
          E2(L) = H*E2(L)
          IF(E2(L).NE.ZERO) GOTO 20
        ENDIF
        P = D(L) + F
C Order the eigenvalues
        IF(L.NE.1) THEN
          DO II = 2,L
            I = L + 2 - II
            IF(P.GE.D(I-1)) GOTO 30
            D(I) = D(I-1)
          ENDDO
        ENDIF
        I = 1
   30   D(I) = P
      ENDDO
      RETURN
C Convergence failure (> 30 interations)
 1000 IERR = L
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBTQL2(NDIM,N,D,E,Z,IERR)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION D(N),E(N),Z(NDIM,N)
      DATA ZERO,ONE,TWO/0.0D0,1.0D0,2.0D0/
C
C  Called by AXWBX
C  [Martin et al., Num. Math. 11, 181 (1968)]
C
C  QL method for eigenvalues and eigenvectors of symmetric tridiagonal
C  matrix (diagonal elements in D, subdiagonal elements in E).
C
C  On return, D contains eigenvalues in ascending order (if IERR = 0),
C  and Z contains the orthonormal eigenvectors.
C
C  IERR = J if J-th eigenvalue failed to converge in 30 iterations.
C
C  RMACHEP is a machine-dependent parameter specifying the relative
C  precision of floating point arithmetic.
C
      RMACHEP = TWO**(-52)
      IERR = 0
      IF (N .EQ. 1) RETURN
      DO I = 2,N
        E(I-1) = E(I)
      ENDDO
      F = ZERO
      B = ZERO
      E(N) = ZERO
      DO L = 1,N
        J = 0
        H = RMACHEP*(ABS(D(L)) + ABS(E(L)))
        IF (B.LT.H) B = H
        DO M = L,N
          IF (ABS(E(M)).LE.B) GO TO 10
        ENDDO
   10   IF (M.EQ.L) GO TO 30
   20   IF (J.EQ.30) GO TO 1000
        J = J + 1
        L1 = L + 1
        G = D(L)
        P = (D(L1) - G)/(TWO*E(L))
        R = SQRT(P*P+ONE)
        D(L) = E(L)/(P + SIGN(R,P))
        H = G - D(L)
        DO I = L1,N
          D(I) = D(I) - H
        ENDDO
        F = F + H
C Perform the QL transformation
        P = D(M)
        C = ONE
        S = ZERO
        MML = M - L
        DO II = 1,MML
          I = M - II
          G = C*E(I)
          H = C*P
          IF(ABS(P).GE.ABS(E(I))) THEN
            C = E(I)/P
            R = SQRT(C*C + ONE)
            E(I+1) = S*P*R
            S = C/R
            C = ONE/R
          ELSE
            C = P/E(I)
            R = SQRT(C*C + ONE)
            E(I+1) = S*E(I)*R
            S = ONE/R
            C = C*S
          ENDIF
          P = C*D(I) - S*G
          D(I+1) = H + S*(C*G + S*D(I))
          DO K = 1,N
            H = Z(K,I+1)
            Z(K,I+1) = S*Z(K,I) + C*H
            Z(K,I) = C*Z(K,I) - S*H
          ENDDO
        ENDDO
        E(L) = S*P
        D(L) = C*P
        IF (ABS(E(L)).GT.B) GO TO 20
   30   D(L) = D(L) + F
      ENDDO
C Order the eigenvalues and eigenvectors
      DO II = 2,N
        I = II - 1
        K = I
        P = D(I)
        DO J = II,N
          IF(D(J).LT.P) THEN
            K = J
            P = D(J)
          ENDIF
        ENDDO
        IF(K.NE.I) THEN
          D(K) = D(I)
          D(I) = P
          DO J = 1, N
            P = Z(J,I)
            Z(J,I) = Z(J,K)
            Z(J,K) = P
          ENDDO
        ENDIF
      ENDDO
      RETURN
C Convergence failure (> 30 iterations)
 1000 IERR = L
      RETURN
      END
C***********************************************************************
      SUBROUTINE TRED1(NDIM,N,A,D,E,E2)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION A(NDIM,N),D(N),E(N),E2(N)
      DATA ZERO/0.0D0/
C
C  Called by AXWBX
C  [Martin et al., Num. Math. 11, 181 (1968)]
C  Reduce real symmetric A to symmetric triadiagonal form (diagonal
C  elements in D, subdiagonal elements in E, squares in E2), using
C  orthogonal similarity transformations.
C
      DO I = 1,N
        D(I) = A(I,I)
      ENDDO
      DO II = 1,N
        I = N + 1 - II
        L = I - 1
        H = ZERO
        SCALE = ZERO
        IF(L.GE.1) THEN
          DO K = 1,L
            SCALE = SCALE + ABS(A(I,K))
          ENDDO
          IF (SCALE.NE.ZERO) GO TO 10
        ENDIF
        E(I) = ZERO
        E2(I) = ZERO
        GO TO 30
   10   DO K = 1,L
          A(I,K) = A(I,K)/SCALE
          H = H + A(I,K)*A(I,K)
        ENDDO
        E2(I) = SCALE*SCALE*H
        F = A(I,L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE*G
        H = H - F*G
        A(I,L) = F - G
        IF (L.EQ.1) GO TO 20
        F = ZERO
        DO J = 1,L
          G = ZERO
          DO K = 1,J
            G = G + A(J,K)*A(I,K)
          ENDDO
          JP1 = J + 1
          IF(L.GE.JP1) THEN
            DO K = JP1,L
              G = G + A(K,J)*A(I,K)
            ENDDO
          ENDIF
          E(J) = G/H
          F = F + E(J)*A(I,J)
        ENDDO
        H = F/(H + H)
        DO J = 1, L
          F = A(I,J)
          G = E(J) - H*F
          E(J) = G
          DO K = 1,J
            A(J,K) = A(J,K) - F*E(K) - G*A(I,K)
          ENDDO
        ENDDO
   20   DO K = 1,L
          A(I,K) = SCALE*A(I,K)
        ENDDO
   30   H = D(I)
        D(I) = A(I,I)
        A(I,I) = H
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE TRED2(NDIM,N,A,D,E,Z)
C***********************************************************************
      IMPLICIT REAL*8(A-H,O-Z),INTEGER(I-N)
      DIMENSION A(NDIM,N),D(N),E(N),Z(NDIM,N)
      DATA ZERO,ONE/0.0D0,1.0D0/
C
C  After Martin et al., Num. Math. 11, 181 (1968)
C
C  Reduce a real symmetric matrix A to symmetric tridiagonal form
C  (with diagonal elements in D, subdiagonal elements in E).
C  Z is the orthogonal transformation matrix produced in reduction.
C  (Used when eigenvectors are wanted.)
C
      DO I = 1,N
        DO J = 1,I
          Z(I,J) = A(I,J)
        ENDDO
      ENDDO
      IF (N .EQ. 1) GO TO 30
      DO II = 2,N
        I = N + 2 - II
        L = I - 1
        H = ZERO
        SCALE = ZERO
        IF(L.GE.2) THEN
          DO K = 1,L
            SCALE = SCALE + ABS(Z(I,K))
          ENDDO
          IF(SCALE.NE.ZERO) GOTO 10
        ENDIF
        E(I) = Z(I,L)
        GOTO 20
   10   DO K = 1,L
          Z(I,K) = Z(I,K)/SCALE
          H = H + Z(I,K)*Z(I,K)
        ENDDO
        F = Z(I,L)
        G = -SIGN(SQRT(H),F)
        E(I) = SCALE*G
        H = H - F*G
        Z(I,L) = F - G
        F = ZERO
        DO J = 1,L
          Z(J,I) = Z(I,J)/H
          G = ZERO
          DO K = 1,J
            G = G + Z(J,K)*Z(I,K)
          ENDDO
          JP1 = J + 1
          IF(L.GE.JP1) THEN
            DO K = JP1,L
              G = G + Z(K,J)*Z(I,K)
            ENDDO
          ENDIF
          E(J) = G/H
          F = F + E(J)*Z(I,J)
        ENDDO
        HH = F/(H + H)
        DO J = 1,L
          F = Z(I,J)
          G = E(J) - HH*F
          E(J) = G
          DO K = 1,J
            Z(J,K) = Z(J,K) - F*E(K) - G*Z(I,K)
          ENDDO
        ENDDO
   20   D(I) = H
      ENDDO
   30 D(1) = ZERO
      E(1) = ZERO
C Accumulate transformation matrices
      DO I = 1,N
        L = I - 1
        IF(D(I).GT.ZERO) THEN
          DO J = 1,L
            G = ZERO
            DO K = 1,L
              G = G + Z(I,K)*Z(K,J)
            ENDDO
            DO K = 1,L
              Z(K,J) = Z(K,J) - G*Z(K,I)
            ENDDO
          ENDDO
        ENDIF
        D(I) = Z(I,I)
         Z(I,I) = ONE
         IF(L.GE.1) THEN
           DO J = 1,L
             Z(I,J) = ZERO
             Z(J,I) = ZERO
           ENDDO
         ENDIF
      ENDDO
      RETURN
      END
C***********************************************************************
C
C  3-CENTER, 4-ELECTRON HYPERBOND (3CHB) SEARCH:
C
C      SUBROUTINE NB3CHB(DM,LFNPR)
C
C***********************************************************************
      SUBROUTINE NB3CHB(DM,LFNPR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C SEARCH FOR 3-CENTER, 4-ELECTRON HYPERBONDS (3CHB KEYWORD)
C
C  REQUIRED INPUT INCLUDES:
C        DM = DENSITY MATRIX IN NHO BASIS (NDIM,NDIM)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORBS(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBL(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
      PARAMETER(MAXHB=MAXATM/3)
      COMMON/NBHB/N3CHB,I3CHB(MAXHB,3)
      DIMENSION DM(NDIM,NDIM)
      DATA LBD/2HBD/LLP/2HLP/LSTAR/1H*/
      DATA HALF/0.5D0/HUNDRD/100.0D0/
      DATA HBTHR/33.33D0/COSTHR/0.4D0/
C
      WRITE(LFNPR,1000) HBTHR
      IF(OPEN.AND.(ALPHA.OR.BETA)) THEN
        OCCTHR = 0.5D0
      ELSE
        OCCTHR = 1.5D0
      ENDIF
      INHO = 0
C FIND A-B :C TRIAD
      DO 50 IBAS = 1,NNAO
        IBDAB = IBXM(IBAS)
        IF(LABEL(IBDAB,2).EQ.LSTAR) GOTO 50
        IF(LABEL(IBDAB,1).EQ.LBD) THEN
          IA = LABEL(IBDAB,4)
          IB = LABEL(IBDAB,5)
C INHOA = NHO HA; INHOB = NHO HB
          INHO = INHO + 1
          INHOA = INHO
          INHO = INHO + 1
          INHOB = INHO
        ELSE
          INHO = INHO + 1
          GOTO 50
        ENDIF
        D11 = DM(INHOA,INHOA)
        D22 = DM(INHOB,INHOB)
        D12 = DM(INHOA,INHOB)
        DIFF = HALF*(D11 - D22)
        AV = HALF*(D11 + D22)
        RAD = SQRT(DIFF*DIFF + D12*D12)
        OCCAB = AV + RAD
        OCCABS = AV - RAD
        JNHO = 0
        DO 40 JBAS = 1,NNAO
          ILPC = IBXM(JBAS)
          IF(LABEL(ILPC,2).EQ.LSTAR) GOTO 40
          IF(LABEL(ILPC,1).EQ.LBD) THEN
             JNHO = JNHO + 2
          ELSE
             JNHO = JNHO + 1
          ENDIF
          IF(LABEL(ILPC,1).EQ.LLP) THEN
            IC = LABEL(ILPC,4)
            IF(IC.EQ.IA.OR.IC.EQ.IB) GOTO 40
            INHOC = JNHO
            OCCLPC = DM(INHOC,INHOC)
          ELSE
            GOTO 40
          ENDIF
C DIAGONALIZE 2X2 NHO DENSITY MATRIX FOR B-C (ITRY=1) OR A-C (ITRY=2)
C TRY A: B-C
          D11 = DM(INHOB,INHOB)
          D22 = DM(INHOC,INHOC)
          D12 = DM(INHOB,INHOC)
          DIFF = HALF*(D11 - D22)
          AV = HALF*(D11 + D22)
          RAD = SQRT(DIFF*DIFF + D12*D12)
          OCCBC = AV + RAD
          OCCBCS = AV - RAD
          COSBC = ABS(D12)/SQRT(D11*D22)
C IS B-C OCCUPANCY ABOVE THRESHOLD?
C TRY B: A-C
          D11 = DM(INHOA,INHOA)
          D22 = DM(INHOC,INHOC)
          D12 = DM(INHOA,INHOC)
          DIFF = HALF*(D11 - D22)
          AV = HALF*(D11 + D22)
          RAD = SQRT(DIFF*DIFF + D12*D12)
          OCCAC = AV + RAD
          OCCACS = AV - RAD
          COSAC = ABS(D12)/SQRT(D11*D22)
C IS A-C OCCUPANCY ABOVE THRESHOLD?
C SELECT BEST ALTERNATIVE BOND
          IF(COSBC.GT.COSAC) THEN
            IF(COSBC.LT.COSTHR) GOTO 40
            IF(OCCBC.LT.OCCTHR) GOTO 40
            IHB1 = IA
            IHB2 = IB
            IHB3 = IC
            PCTAB = HUNDRD*OCCBCS/(OCCABS + OCCBCS)
          ELSE
            IF(COSAC.LT.COSTHR) GOTO 40
            IF(OCCAC.LT.OCCTHR) GOTO 40
            IHB1 = IB
            IHB2 = IA
            IHB3 = IC
            PCTAB = HUNDRD*OCCACS/(OCCABS + OCCACS)
          ENDIF
          PCTALT = HUNDRD - PCTAB
          IF(PCTALT.GT.HBTHR) THEN
            N3CHB = N3CHB + 1
            I3CHB(N3CHB,1) = IHB1
            I3CHB(N3CHB,2) = IHB2
            I3CHB(N3CHB,3) = IHB3
            TOCC = OCCAB + OCCABS + OCCLPC
C Print another line of the hyperbond table
            WRITE(LFNPR,1100) N3CHB,NAMEAT(IATNO(IHB1)),IHB1,
     +         NAMEAT(IATNO(IHB2)),IHB2,NAMEAT(IATNO(IHB3)),
     +         IHB3,PCTAB,PCTALT,TOCC,IBAS,JBAS,INHOA,INHOB,INHOC
          ENDIF
   40   CONTINUE
   50 CONTINUE
      IF(N3CHB.LE.0) WRITE(LFNPR,1200)
C                                                  NBOs       3-center hybrids
C                                              -------------  ----------------
C       Hyperbond A:-B-:C   %A-B/%B-C   OCC.   BD(A-B) LP(C)  h(A)  h(B)  h(C)
C      -------------------  ---------  ------  ------- -----  ----  ----  ----
C   1.  F  1:- F  2-:Cl  3  53.0/47.0  3.9988      1     18      1     2    19
 1000 FORMAT(//,1X,'3-Center, 4-Electron A:-B-:C Hyperbonds ',
     + '(A-B :C <=> A: B-C)',
     + /,14X,'[threshold for detection:',F5.1,'%]',/,
     + /,50X,'NBOs',7X,'3-center hybrids',
     + /,46X,'-------------  ----------------',
     + /,1X,'      Hyperbond A:-B-:C   %A-B/%B-C   OCC. ',
     + '  BD(A-B) LP(C)  h(A)  h(B)  h(C)',
     + /,1X,'     -------------------  ---------  ------',
     + '  ------- -----  ----  ----  ----')
 1100 FORMAT(1X,I3,'.',1X,A2,I3,':-',A2,I3,'-:',A2,I3,2X,
     + F4.1,'/',F4.1,1X,F7.4,I8,3X,I4,3(2X,I4))
 1200 FORMAT(1X,'     (none detected)')
      RETURN
      END
C***********************************************************************
C
C  NATURAL J-COUPLING ROUTINES (LFN 91, 92 USED AS SCRATCH)
C
C      SUBROUTINE NJC(TBL,FCG,FCB,FCL,TBLA,TBLB,FCBA,FCBB,FCLA,FCLB,
C                 DJI_NL,DI_L,DPJI_NL,DPI_L,SCR1,SCR2,V1,V2,V3,V4,LFNPR)
C      SUBROUTINE NJCANL(IATNP,TBLA,TBLB,FCBA,FCBB,FCLA,FCLB,
C     +           DJI_NL,DI_L,DPJI_NL,DPI_L,SCR1,SCR2,V1,V2,V3,V4,LFNPR)
C      SUBROUTINE NJCTAB(IATN,IATNP,DJI_NL,DI_L,DELOC,REPOL,
C     +           SUMBI,SUMBJ,DM,T,LFNPR)
C      SUBROUTINE NJCATM(IAT)
C      SUBROUTINE NJCLNK(IA,IB,LINKS,NLINKS,IERR)
C      SUBROUTINE NJCLAB(IATN,IATNP,LINKS,NLINKS,MOLLOC,TITLE)
C      SUBROUTINE FNDHYBS(INBO,IH1,IH2,IH3)
C      SUBROUTINE CONSTFC(CGYRO,IZ1,IZ2,IPR,LFNPR)
C      SUBROUTINE RMBLNK(STR,NS)
C      FUNCTION IPHTYP(IBO,JBO)
C      SUBROUTINE FCAOEL(FCAO,SCR,ISCR,NBAS,NDIM,NATOMS)
C      SUBROUTINE NBRSTR(FVAL,STR,NS)
C      SUBROUTINE NBISTR(INT0,STR,NS)
C      SUBROUTINE AOVAL2(AOAMP,R,ORIGIN)
C      FUNCTION AOAFAC(R,ITYPE)
C      SUBROUTINE ITYPERR(ITYPE)
C      SUBROUTINE STRUCF(ITYPE,NKAL,CAL,NXAL,NYAL,NZAL)
C      SUBROUTINE BNORM
C      FUNCTION AONORM(ITYPE,ZETA1,ZETA2)
C      FUNCTION ATMFAC(R,N)
C      FUNCTION OLAPCG(N1,ZETA1,C1,N2,ZETA2,C2)
C      FUNCTION OLAP1D(NX1,ZETA1,X1,NX2,ZETA2,X2,ALFA)
C
C***********************************************************************
      SUBROUTINE NJC(TBL,FCG,FCB,FCL,
     +   TBLA,TBLB,FCBA,FCBB,FCLA,FCLB,
     +   DJI_NL,DI_L,DPJI_NL,DPI_L,SCR1,SCR2,V1,V2,V3,V4,LFNPR)
C***********************************************************************
C 17-Jan-01  FAW  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C TBL = NBONLMO transformation
C FCG = Gaussian AO Fermi contact amplitudes
C FCB = NBO Fermi contact amplitudes
C FCL = NLMO Fermi contact amplitudes
C
      PARAMETER(MAXBAS = 2000)
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      DIMENSION TBL(NDIM,NDIM),FCG(NDIM,NATOMS),
     +  FCB(NDIM,NATOMS),FCL(NDIM,NATOMS)
      DIMENSION TBLA(NDIM,NDIM),TBLB(NDIM,NDIM),FCBA(NDIM,NATOMS),
     +  FCBB(NDIM,NATOMS),FCLA(NDIM,NATOMS),FCLB(NDIM,NATOMS),
     +  SCR1(NDIM,NDIM),SCR2(NDIM,NDIM),
     +  V1(NDIM),V2(NDIM),V3(NDIM),V4(NDIM),
     +  DJI_NL(NDIM,NLEW),DI_L(NLEW)
      DATA ZERO/0.0D0/
      SAVE ZERO
C  Find the perturbed atom (from job filename)
      CALL NJCATM(IATNP)
C  Get the Fermi-contact amplitudes (FCG) in the Gaussian AO basis
C  (using FCB a scratch)
      CALL FCAOEL(FCG,FCB,FCB,NBAS,NDIM,NATOMS)
C  Get the AO to NBO transform (store in TBL, for now)
      CALL FETNBO(TBL)
C  Form the NBO Fermi-contact matrix FCB from (TBL) and FCG
      DO IQ = 1,NNAO
        DO IATN = 1,NATOMS
          FCBQN = ZERO
          DO IP = 1,NBAS
            FCBQN = FCBQN + TBL(IP,IQ)*FCG(IP,IATN)
          ENDDO
          FCB(IQ,IATN) = FCBQN
        ENDDO
      ENDDO
C  Get the AO to NLMO transform (store in TBL, for now)
      NFILE = 46
      IF(ISPIN.EQ.-2) NFILE = 47
      L3 = NDIM*NDIM
      CALL NBREAD(TBL,L3,NFILE)
C  Form the NLMO Fermi-contact matrix FCL from (TBL) and FCG
      DO IQ = 1,NNAO
        DO IATN = 1,NATOMS
          FCLQN = ZERO
          DO IP = 1,NBAS
            FCLQN = FCLQN + TBL(IP,IQ)*FCG(IP,IATN)
          ENDDO
          FCL(IQ,IATN) = FCLQN
        ENDDO
      ENDDO
C  Fetch the NBO-NLMO transformation TBL
      CALL FETLMO(TBL)
C  Write the quantities for this spin set to disk (LFN 91-93)
      IF(OPEN.AND.ALPHA) THEN
        OPEN(91,FILE='NJCA.91',STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(91) TBL,FCB,FCL
        CLOSE(91)
      ELSE IF(OPEN.AND.BETA) THEN
        OPEN(92,FILE='NJCB.92',STATUS='UNKNOWN',FORM='UNFORMATTED')
        WRITE(92) TBL,FCB,FCL
        CLOSE(92)
      ENDIF
C  Analyze the results?
      CALL NJCANL(IATNP,TBLA,TBLB,FCBA,FCBB,FCLA,FCLB,
     +   DJI_NL,DI_L,DPJI_NL,DPI_L,SCR1,SCR2,V1,V2,V3,V4,LFNPR)
      RETURN
      END
C***********************************************************************
      SUBROUTINE NJCANL(IATNP,TBLA,TBLB,FCBA,FCBB,FCLA,FCLB,
     +   DJI_NL,DI_L,DPJI_NL,DPI_L,SCR1,SCR2,V1,V2,V3,V4,LFNPR)
C***********************************************************************
C 17-Jan-01  FAW  New subroutine to analyze Fermi contact J-couplings
C                 calculated by Barfield finite-field method
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER STR*3,FNAME*256,FILENM*256
      PARAMETER(MAXBAS = 2000)
      PARAMETER(MAXFIL = 40)
      COMMON/NBNAME/FILENM,NFILE,IFILE(MAXFIL)
C
C  Analyze the J-coupling arrays stored in NCSA.91, NCSB.92
C
      COMMON/NBFLAG/ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      LOGICAL ROHF,UHF,CI,OPEN,COMPLX,ALPHA,BETA,MCSCF,AUHF,ORTHO
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      DIMENSION TBLA(NDIM,NDIM),TBLB(NDIM,NDIM),FCBA(NDIM,NATOMS),
     +  FCBB(NDIM,NATOMS),FCLA(NDIM,NATOMS),FCLB(NDIM,NATOMS),
     +  SCR1(NDIM,NDIM),SCR2(NDIM,NDIM),
     +  V1(NDIM),V2(NDIM),V3(NDIM),V4(NDIM),
     +  DJI_NL(NDIM,NLEW),DI_L(NLEW),DPJI_NL(NDIM,NLEW),DPI_L(NLEW)
      DATA ZERO/0.0D0/TWO/2.0D0/HALF/0.5D0/
      SAVE ZERO,TWO,HALF
      IF(.NOT.OPEN.OR..NOT.BETA) RETURN
C
C  Read transformation and Fermi contact spin matrices from disk (LFN 91, 92)
C
      OPEN(91,FILE='NJCA.91',STATUS='OLD',FORM='UNFORMATTED')
      OPEN(92,FILE='NJCB.92',STATUS='OLD',FORM='UNFORMATTED')
C  Alpha spin-perturbed
      READ(91,ERR=500) TBLA,FCBA,FCLA
      CLOSE(91)
C  Beta spin-perturbed
      READ(92,ERR=500) TBLB,FCBB,FCLB
      CLOSE(92)
C  Open the file (FILENM_n.njc) of Delta matrices for this (current) atom
      NL = LENNB(FILENM)
      FNAME = FILENM(1:NL)//'.njc'
      OPEN(91,FILE=FNAME(1:NL+4),STATUS='UNKNOWN',
     +     FORM='UNFORMATTED')
C  Main loop over atoms
      WRITE(LFNPR,1000)
      NL = LENNB(FNAME)
      WRITE(LFNPR,1100) (FNAME(K:K),K=1,NL)
      DO IATN = 1,NATOMS
C  Evaluate TOTALC = total coupling as sum over occupied NLMOs (I)
        TOTALC = ZERO
        DO I = 1,NLEW
          TOTALC = TOTALC + FCLA(I,IATN)**2 - FCLB(I,IATN)**2
C  Evaluate the diagonal DI_L (Lewis NBO) contributions
          DI_L(I) = FCBA(I,IATN)**2 - FCBB(I,IATN)**2
        ENDDO
C  Evaluate the off-diagonal DJI_NL contributions
        DO I = 1,NLEW
          DO J=1,NNAO
            DJI_NL(J,I) = ZERO
            IF(J.NE.I) THEN
C  Alpha spin
              T1A=(FCBA(J,IATN)**2 - FCBA(I,IATN)**2)*TBLA(J,I)**2
              T2A=TWO*TBLA(I,I)*TBLA(J,I)*FCBA(I,IATN)*FCBA(J,IATN)
              T3A=ZERO
              DO K = 1,NNAO
                IF(K.NE.I.AND.K.NE.J) T3A=T3A+TBLA(K,I)*FCBA(K,IATN)
              ENDDO
              T3A=T3A*TBLA(J,I)*FCBA(J,IATN)
              TA=T1A+T2A+T3A
C  Beta spin
              T1B=(FCBB(J,IATN)**2 - FCBB(I,IATN)**2)*TBLB(J,I)**2
              T2B=TWO*TBLB(I,I)*TBLB(J,I)*FCBB(I,IATN)*FCBB(J,IATN)
              T3B=ZERO
              DO K = 1,NNAO
                IF(K.NE.I.AND.K.NE.J) T3B=T3B+TBLB(K,I)*FCBB(K,IATN)
              ENDDO
              T3B=T3B*TBLB(J,I)*FCBB(J,IATN)
              TB=T1B+T2B+T3B
C  New spin-density contribution
              DJI_NL(J,I) = TA - TB
            ENDIF
          ENDDO
C  Store total NLMO contribution in diagonal term of DJI_NL array
          DJI_NL(I,I) = FCLA(I,IATN)**2 - FCLB(I,IATN)**2
        ENDDO
C  Write the Delta matrices for this atom to disk file FILENM_IATNP.njc
        WRITE(91) DJI_NL,DI_L
C  Was this atom previously perturbed (saved in FILENM_IATN.njc)?
        CALL NBISTR(IATN,STR,NS)
        NL = LENNB(FILENM)
        K = NL
C  Patch for Sun OS (which passes OPEN(92,...,ERR=20) without error
C  even if the same file is opened for two LFNs).  FW 2/6/02
        IF(IATN.EQ.IATNP) GOTO 20
   10   IF(FILENM(K:K).EQ.'_') THEN
          FNAME = FILENM(1:K)//STR(1:NS)//'.njc'
          NF = K+NS+4
          OPEN(92,FILE=FNAME(1:NF),STATUS='OLD',FORM='UNFORMATTED',
     +         ERR=20)
        ELSE IF(K.GT.0) THEN
          K = K - 1
          GOTO 10
        ELSE
          GOTO 20
        ENDIF
C  Read the previous Delta matrices for this atom, store in DPJI_NL, DP_I
        DO I=1,IATNP
          READ(92) DPJI_NL,DPI_L
        ENDDO
C  Symmetrize the Delta matrices
        DO KC = 1,NLEW
          DI_L(KC) = HALF*(DI_L(KC) + DPI_L(KC))
          DO KR = 1,NNAO
            DJI_NL(KR,KC) = HALF*(DJI_NL(KR,KC) + DPJI_NL(KR,KC))
          ENDDO
        ENDDO
C  Print the analysis table for this pair of nuclei
        CALL NJCTAB(IATN,IATNP,DJI_NL,DI_L,V1,V2,V3,V4,SCR1,SCR2,LFNPR)
   20   CONTINUE
        CLOSE(92)
      ENDDO
      CLOSE(91)
      RETURN
C  Error: no spin matrices
  500 CONTINUE
      WRITE(LFNPR,1500)
      RETURN
 1000 FORMAT(//,1X,'NATURAL J-COUPLING ANALYSIS')
 1100 FORMAT(1X,'Writing disk file ',80A1)
 1500 FORMAT(/,1X,'*** Error: No spin-matrix info in LFN = 91, 92')
      END
C***********************************************************************
      SUBROUTINE NJCTAB(IATN,IATNP,DJI_NL,DI_L,DELOC,REPOL,
     + SUMBI,SUMBJ,DM,T,LFNPR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER DASHES*11,TABROW*8,STR*10,BLANK*8,BUFFER*10,TITLE*80
      CHARACTER JLABEL*80,CHARAT*2
C
C  Print the NJC table for atom IATN
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      DIMENSION DJI_NL(NDIM,NLEW),DI_L(NLEW),TABROW(7),
     + IBIG(MAXBAS),DELOC(NLEW),JBIG(MAXBAS),SUMBI(NDIM),SUMBJ(NLEW)
      DIMENSION REPOL(NLEW),DJ(MAXBAS),MOLLOC(MAXATM),V(MAXBAS),
     + DM(NDIM,NDIM),T(NDIM,NDIM),VALROW(7)
      DIMENSION LTOTAL(10),LOTHERS(10),IFLAG(100),LINKS(12),LINKS2(12)
C Intermolecular donor-acceptor interactions
      PARAMETER(NMXDA=100)
      COMMON/DAINT/EDA(NMXDA),NDA,INBOD(NMXDA),INBOA(NMXDA),
     +       IATMD(NMXDA),IATMA(NMXDA),IMOLD(NMXDA),IMOLA(NMXDA)
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBL0(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORB(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      COMMON/NBMOL/NMOLEC,MOLAT(MAXATM),MOLEC(MAXATM,MAXATM),
     +              NMOLA,MOLATA(MAXATM),MOLECA(MAXATM,MAXATM)
      COMMON/NBTHR/THRSET,PRJSET,ACCTHR,CRTSET,E2THR,ATHR,PTHR,ETHR,
     +             DTHR,DLTHR,CHSTHR,REFTHR,STTHR,PRTHR,THRNCS,THRNJC
      DATA DASHES/'-----------'/ZERO/0.0D0/
      DATA LOTHERS/1H ,1H ,1HO,1Ht,1Hh,1He,1Hr,1Hs,1H ,1H|/
      DATA LTOTAL/1H ,1H ,1H ,1HT,1Ho,1Ht,1Ha,1Hl,1H ,1H /
      DATA BLANK/'        '/IBLNK/1H /IINT/1Hi/IBAR/1H|/ISW/0/
      SAVE DASHES,ZERO,BLANK,IBLNK,IINT,IBAR,ISW,MOLLOC
      IPR = 0
C  Prepare the array (MOLLOC) of molecular unit locations
C  for each atom
      IF(ISW.EQ.0) THEN
        ISW = 1
        IF(NMOLEC.GT.0) THEN
          DO IMOL = 1,NMOLEC
            NAT = MOLAT(IMOL)
            DO IAT = 1,NAT
              IA = MOLEC(IMOL,IAT)
              MOLLOC(IA) = IMOL
            ENDDO
          ENDDO
        ELSE
          DO IAT = 1,NATOMS
            MOLLOC(IAT) = 1
          ENDDO
        ENDIF
C  Fill the IATMD, IATMA vectors of coordinated D-A atoms in COMMON/DAINT/
        IF(NDA.GT.0) THEN
C  Fetch the density matrix in the NAO basis
          CALL FEDNAO(DM)
C  Fetch the NAO-NHO transformation
          CALL FETNHO(T)
C  Transform DM to NHO basis
          CALL SIMTRS(DM,T,V,NDIM,NNAO)
C  Search the list of donor-acceptor interactions
          DO IDA = 1,NDA
C  Find the hybrids for each NBO of donor (IBD) or acceptor (IBA) type
            IBD = INBOD(IDA)
            IBA = INBOA(IDA)
            CALL FNDHYBS(IBD,IHD1,IHD2,IHD3)
            CALL FNDHYBS(IBA,IHA1,IHA2,IHA3)
C  Which hybrids have the strongest interaction between units?
            D1A1 = ABS(DM(IHD1,IHA1))
            D1A2 = ZERO
            D2A1 = ZERO
            D2A2 = ZERO
            IF(IHA2.NE.0) D1A2 = ABS(DM(IHD1,IHA1))
            IF(IHD2.NE.0) THEN
              D2A1 = ABS(DM(IHD2,IHA1))
              IF(IHA2.NE.0) D2A2 = ABS(DM(IHD2,IHA2))
            ENDIF
C  Donor hybrid (1 or 2)?
            IF((D1A1+D1A2).GT.(D2A1+D2A2)) THEN
              NHD = 1
            ELSE
              NHD = 2
            ENDIF
C  Acceptor hybrid (1 or 2)?
            IF((D1A1+D2A1).GT.(D1A2+D2A2)) THEN
              NHA = 1
            ELSE
              NHA = 2
            ENDIF
C  What atoms go with these hybrids?
            IATMD(IDA) = LABEL(IBD,NHD+3)
            IATMA(IDA) = LABEL(IBA,NHA+3)
          ENDDO
        ENDIF
      ENDIF
C  Evaluate the multiplicative constant of gyromagnetic ratios
      CALL CONSTFC(CGYRO,IATNO(IATN),IATNO(IATNP),IPR,LFNPR)
C  Multiply all entries by CGYRO
      DO ICOL = 1,NLEW
        DI_L(ICOL) = CGYRO*DI_L(ICOL)
        DELOC(ICOL) = ZERO
        DO JROW = 1,NNAO
          DJI_NL(JROW,ICOL) = CGYRO*DJI_NL(JROW,ICOL)
C  Assemble the total non-Lewis DELOC (including repolarization
C  type, for now) in each column
          IF(JROW.GT.NLEW) THEN
            DELOC(ICOL) = DELOC(ICOL) + DJI_NL(JROW,ICOL)
C  Reassign Lewis-type DJI_NL ("non-Lewis") entries to DI_L
          ELSE IF(JROW.NE.ICOL) THEN
            DI_L(ICOL) = DI_L(ICOL) + DJI_NL(JROW,ICOL)
            DJI_NL(JROW,ICOL) = ZERO
          ENDIF
        ENDDO
      ENDDO
C  Assemble the total spin coupling from DI_L and DJI_NL elements
      TOTALC = ZERO
      TOTL = ZERO
      DO ICOL = 1,NLEW
        TOTL = TOTL + DI_L(ICOL)
        TOTALC = TOTALC + DJI_NL(ICOL,ICOL)
        DJI_NL(ICOL,ICOL) = ZERO
      ENDDO
C  Print the total spin coupling for this atom
      IF(ABS(TOTALC).LT.THRNJC) RETURN
C  Intermolecular spin-coupling?
      IMA = MOLLOC(IATN)
      IMB = MOLLOC(IATNP)
      IF(IMA.NE.IMB) THEN
C  Search the coordination interactions in COMMON/DAINT/ for largest EDA
        IF(NDA.LE.0) THEN
          TITLE = ' '
          NT = 1
          GOTO 5
        ENDIF
        IBIGDA = 0
        EBIGDA = ZERO
        DO IDA = 1,NDA
          IF(IMA.EQ.IMOLD(IDA).AND.IMB.EQ.IMOLA(IDA).OR.
     +       IMB.EQ.IMOLD(IDA).AND.IMA.EQ.IMOLA(IDA)) THEN
            IF(EDA(IDA).GT.EBIGDA) THEN
              EBIGDA = EDA(IDA)
              IBIGDA = IDA
            ENDIF
          ENDIF
        ENDDO
C  Find connecting links from IATNP to the coordination
C  bridge atom IATA
        IDA = IBIGDA
        IATA = IATMA(IDA)
        IATD = IATMD(IDA)
        IF(MOLLOC(IATA).EQ.MOLLOC(IATNP)) THEN
          IA1 = IATNP
          IA2 = IATN
        ELSE
          IA1 = IATN
          IA2 = IATNP
        ENDIF
        CALL NJCLNK(IA1,IATA,LINKS,NLINKS,IERR1)
C  Find connecting links from coordination bridge atom IATD to IATN
        CALL NJCLNK(IATD,IA2,LINKS2,NLINKS2,IERR2)
        IF(IERR1.NE.0.OR.IERR2.NE.0.OR.(NLINKS+NLINKS2).GT.12) THEN
          TITLE = ' '
          NT = 1
          GOTO 5
        ENDIF
C  Combine the two LINKS lists
        IF(IATD.NE.IA2) THEN
          NLINKS = NLINKS + 1
          IF(NLINKS.GT.1) LINKS(NLINKS-1) = IATA
          NLINKS = NLINKS + 1
          LINKS(NLINKS-1) = IATD
          IF(NLINKS2.GT.1) THEN
            DO IL2 = 1,NLINKS2-1
              NLINKS = NLINKS +1
              LINKS(NLINKS-1) = LINKS2(IL2)
            ENDDO
          ENDIF
        ELSE
          NLINKS = NLINKS + 1
          IF(NLINKS.GT.1) LINKS(NLINKS-1) = IATA
        ENDIF
      ELSE
C  or intramolecular spin-coupling?
        CALL NJCLNK(IATN,IATNP,LINKS,NLINKS,IERR)
        IF(IERR.NE.0) THEN
          TITLE = ' '
          NT = 1
          GOTO 5
        ENDIF
        IA1 = IATN
        IA2 = IATNP
      ENDIF
C  Form the n-J[...] label
      CALL NJCLAB(IA1,IA2,LINKS,NLINKS,MOLLOC,TITLE)
      NT = LENNB(TITLE)
    5 CONTINUE
      CALL NBISTR(IATN,STR,NS)
      JLABEL = 'J['//CHARAT(IATNO(IATN))//STR(1:NS)//','
      NJL = NS + 5
      CALL NBISTR(IATNP,STR,NS)
      JLABEL = JLABEL(1:NJL)//CHARAT(IATNO(IATNP))//STR(1:NS)//']'
      NJL = NJL +NS + 3
C  Remove blanks
      CALL RMBLNK(JLABEL,NJL)
C Compose label, value, and TB title
      CALL NBRSTR(TOTALC,STR,NS)
      JLABEL = JLABEL(1:NJL)//' = '//STR(1:NS)//' Hz:  '//TITLE(1:NT)
      NJL = NJL + NT + NS + 9
      WRITE(LFNPR,1000) (JLABEL(K:K),K=1,NJL)
      WRITE(LFNPR,1150) THRNJC
C   What are the significant non-Lewis contributions J?  Store in JBIG
      NJBIG = 0
      DO J = NLEW+1,NNAO
        TJMAX = ZERO
        DJ(J) = ZERO
        DO IP = 1,NLEW
          IF(IPHTYP(IP,J).NE.IINT) THEN
            DJ(J) = DJ(J) + DJI_NL(J,IP)
            TMX = ABS(DJI_NL(J,IP))
            IF(TMX.GT.TJMAX) TJMAX = TMX
          ENDIF
        ENDDO
        IF(TJMAX.GT.THRNJC) THEN
          NJBIG = NJBIG + 1
          JBIG(NJBIG) = J
        ENDIF
      ENDDO
C   What are the significant parent contributions I?  Store in IBIG
      NIBIG = 0
      SUMBL = ZERO
      DO I = 1,NLEW
        TIMAX = ZERO
        DO J = NLEW+1,NNAO
          IF(J.NE.I) THEN
            TMX = ABS(DJI_NL(J,I))
            IF(TMX.GT.TIMAX) TIMAX = TMX
          ENDIF
        ENDDO
        IF(TIMAX.GT.THRNJC) THEN
          NIBIG = NIBIG + 1
          IBIG(NIBIG) = I
          SUMBL = SUMBL + DI_L(I)
        ENDIF
      ENDDO
C  Determine the sum (SUMBJ) of big-J (JBIG) entries in
C  each column ICOL = 1,NLEW
      DO ICOL = 1,NLEW
        SUMBJ(ICOL) = ZERO
        DO JB=1,NJBIG
          J = JBIG(JB)
          IF(IPHTYP(ICOL,J).NE.IINT)
     +      SUMBJ(ICOL) = SUMBJ(ICOL) + DJI_NL(J,ICOL)
        ENDDO
      ENDDO
C  Determine the sum (SUMBI) of big-I (IBIG) entries in
C  each row JROW = 1, NNAO
      DO JROW = NLEW+1,NNAO
        SUMBI(JROW) = ZERO
        DO IB = 1,NIBIG
          I = IBIG(IB)
          IF(IPHTYP(I,JROW).NE.IINT)
     +      SUMBI(JROW) = SUMBI(JROW) + DJI_NL(JROW,I)
        ENDDO
      ENDDO
C  Determine the repolarization correction REPOL(I) for each I = 1,NLEW
C  Also, determine sum (TREPOL) of REPOL(I) and sum (TBREPOL)
C  of big-I REPOL(I)
      DO I = 1,NLEW
        REPOL(I) = ZERO
        DO J = NLEW+1,NNAO
          IF(IPHTYP(I,J).EQ.IINT) THEN
            REPOL(I) = REPOL(I) + DJI_NL(J,I)
          ENDIF
        ENDDO
        DELOC(I) = DELOC(I) - REPOL(I)
      ENDDO
C Gather the sum (TREPOL) of all repolarizations
      TREPOL = ZERO
      DO I = 1,NLEW
        TREPOL = TREPOL + REPOL(I)
      ENDDO
C Gather the sum (TBREPOL) of big repolarization
      TBREPOL = ZERO
      DO IB = 1,NIBIG
        I = IBIG(IB)
        TBREPOL = TBREPOL + REPOL(I)
      ENDDO
C  Determine the "Other Others" and "Total Others" correction
      TOTHER = ZERO
      DO I = 1,NLEW
        TOTHER = TOTHER + DELOC(I) - SUMBJ(I)
      ENDDO
      TBOTHER = ZERO
      DO I = 1,NIBIG
        TBOTHER = TBOTHER + DELOC(IBIG(I)) - SUMBJ(IBIG(I))
      ENDDO
C  Determine the sum TDELOC of all DELOC(I) values
      TDELOC = ZERO
      DO I = 1,NLEW
        TDELOC = TDELOC + DELOC(I)
      ENDDO
C  Determine the sum TBDELOC of big-I DELOC(I) values
      TBDELOC = ZERO
      DO IB = 1,NIBIG
        I = IBIG(IB)
        TBDELOC = TBDELOC + DELOC(I)
      ENDDO
C  Begin main loop over occupied NLMOs (columns ILEW)
      MAXCOL = 5
      IB2 = 0
   60 CONTINUE
      IB1 = IB2 + 1
      IB2 = IB1 + MAXCOL - 1
      IF(IB2.GT.NIBIG) THEN
        ILAST = 1
        IB2 = NIBIG
        NCOL = IB2 - IB1 + 1
      ELSE
        ILAST = 0
        NCOL = MAXCOL
      ENDIF
C  Print the headings
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1200) (IBIG(IB),IB=IB1,IB2)
        WRITE(LFNPR,1250) ((LBL(K,IBIG(IB),4),K=1,10),IB=IB1,IB2)
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL)
      ELSE
        WRITE(LFNPR,1200) (IBIG(IB),IB=IB1,IB2)
        WRITE(LFNPR,1250) ((LBL(K,IBIG(IB),4),K=1,10),IB=IB1,IB2),
     +   (LOTHERS(K),K=1,10),(LTOTAL(K),K=1,8)
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL+2)
      ENDIF
C  Print the diagonal Lewis contribution
      ICOL = 0
      DO IB = IB1,IB2
        ICOL = ICOL + 1
        CALL NBRSTR(DI_L(IBIG(IB)),STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
      ENDDO
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1300) (TABROW(IP),IFLAG(IP),IP=1,NCOL)
      ELSE
        OTHERL = TOTL - SUMBL
        ICOL = ICOL + 1
        CALL NBRSTR(OTHERL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBAR
        ICOL = ICOL + 1
        CALL NBRSTR(TOTL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
        WRITE(LFNPR,1305) (TABROW(IP),IFLAG(IP),IP=1,NCOL+1),
     +        TABROW(NCOL+2)
      ENDIF
C  Print the "Repol." contribution
      ICOL = 0
      DO IB = IB1,IB2
        ICOL = ICOL + 1
        CALL NBRSTR(REPOL(IBIG(IB)),STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
      ENDDO
      OTHREPOL = TREPOL - TBREPOL
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1310) (TABROW(IP),IFLAG(IP),IP=1,NCOL)
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL)
      ELSE
        ICOL = ICOL + 1
        CALL NBRSTR(OTHREPOL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBAR
        ICOL = ICOL + 1
        CALL NBRSTR(TREPOL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
        WRITE(LFNPR,1315) (TABROW(IP),IFLAG(IP),IP=1,NCOL+1),
     +        TABROW(NCOL+2)
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL+2)
      ENDIF
C  Print the significant non-Lewis contributions (rows JBIG(IB))
      IPR = 1
      DO ICOL = 1,7
        VALROW(ICOL) = ZERO
      ENDDO
      TOTHBJ = ZERO
      OTHBJ = ZERO
      DJACC = ZERO
      DO 10 JB = 1,NJBIG
        JP = JBIG(JB)
C  Load (add?) the next DJI_NL entries into VALROW
        ICOL = 0
        DO IB = IB1,IB2
          ICOL = ICOL + 1
          VALROW(ICOL) = VALROW(ICOL) + DJI_NL(JP,IBIG(IB))
        ENDDO
        OTHBJ = OTHBJ + DJ(JP) - SUMBI(JP)
        DJACC = DJACC + DJ(JP)
C  Will this be further accumulated (RY* labels)?
        IF(JB.LT.NJBIG.AND.THRNJC.GT.ZERO) THEN
          JP1 = JBIG(JB+1)
          IACC = 1
          DO K = 1,10
            IF(LBL(K,JP,4).NE.LBL(K,JP1,4)) IACC = 0
          ENDDO
          IF(IACC.EQ.1) THEN
            IPR = 0
            GOTO 10
          ENDIF
        ELSE
          IACC = 0
          IPR = 1
          TOTHBJ = TOTHBJ + OTHBJ
        ENDIF
C  Prepare the VALROW entries for printing
        ICOL = 0
        DO IB = IB1,IB2
          ICOL = ICOL + 1
          CALL NBRSTR(VALROW(ICOL),STR,NS)
          IF(NS.LT.8) THEN
            TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
          ELSE IF(NS.EQ.8) THEN
            TABROW(ICOL) = STR(1:NS)
          ELSE
            TABROW(ICOL) = '********'
          ENDIF
          IFLAG(ICOL) = IPHTYP(IBIG(IB),JP)
          IF(IFLAG(ICOL).EQ.IINT) THEN
            TABROW(ICOL) = '        '
            IFLAG(ICOL) = IBLNK
          ENDIF
        ENDDO
        IF(ILAST.EQ.0) THEN
          IF(IPR.EQ.1) THEN
            WRITE(LFNPR,1700) JP,(LBL(K,JP,4),K=1,10),
     +         (TABROW(IP),IFLAG(IP),IP=1,NCOL)
          ELSE
            IPR = 1
            WRITE(LFNPR,1710) (LBL(K,JP,4),K=1,10),
     +         (TABROW(IP),IFLAG(IP),IP=1,NCOL)
          ENDIF
        ELSE
          CALL NBRSTR(OTHBJ,STR,NS)
          IF(NS.LT.8) THEN
            TABROW(NCOL+1) = BLANK(1:8-NS)//STR(1:NS)
          ELSE IF(NS.EQ.8) THEN
            TABROW(NCOL+1) = STR(1:NS)
          ELSE
            TABROW(NCOL+1) = '********'
          ENDIF
          IFLAG(NCOL+1) = IBAR
          CALL NBRSTR(DJACC,STR,NS)
          IF(NS.LT.8) THEN
            TABROW(NCOL+2) = BLANK(1:8-NS)//STR(1:NS)
          ELSE IF(NS.EQ.8) THEN
            TABROW(NCOL+2) = STR(1:NS)
          ELSE
            TABROW(NCOL+2) = '********'
          ENDIF
          IFLAG(NCOL+2) = IBLNK
          IF(IPR.EQ.1) THEN
            WRITE(LFNPR,1700) JP,(LBL(K,JP,4),K=1,10),
     +         (TABROW(IP),IFLAG(IP),IP=1,NCOL),
     +         TABROW(NCOL+1),IFLAG(NCOL+1),
     +         TABROW(NCOL+2)
          ELSE
            IPR = 1
            WRITE(LFNPR,1710) (LBL(K,JP,4),K=1,10),
     +         (TABROW(IP),IFLAG(IP),IP=1,NCOL),
     +         TABROW(NCOL+1),IFLAG(NCOL+1),
     +         TABROW(NCOL+2)
          ENDIF
        ENDIF
        DO IC = 1,7
          VALROW(IC) = ZERO
        ENDDO
        OTHBJ = ZERO
        DJACC = ZERO
   10 CONTINUE
C  Write the "others" contribution for these columns
      ICOL = 0
      DO IB = IB1,IB2
        ICOL = ICOL + 1
        VAL = DELOC(IBIG(IB))-SUMBJ(IBIG(IB))
        CALL NBRSTR(VAL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
      ENDDO
      OTHERO = TOTHER - TBOTHER
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1750) (TABROW(IP),IFLAG(IP),IP=1,NCOL)
      ELSE
        ICOL = ICOL + 1
        CALL NBRSTR(OTHERO,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBAR
        ICOL = ICOL + 1
        CALL NBRSTR(TOTHER,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
        WRITE(LFNPR,1755) (TABROW(IP),IFLAG(IP),IP=1,NCOL+1),
     +        TABROW(NCOL+2)
      ENDIF
C  Write the "(Total deloc.)" entries for these columns
      ICOL = 0
      DO IB = IB1,IB2
        ICOL = ICOL + 1
        VAL = DELOC(IBIG(IB))
        CALL NBRSTR(VAL,STR,NS)
        BUFFER='('//STR(1:NS)//')'
        NS = NS + 2
        STR=BUFFER(1:NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
      ENDDO
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1770) (TABROW(IP),IFLAG(ICOL),IP=1,NCOL)
      ELSE
        ICOL = ICOL + 1
        VAL = TDELOC - TBDELOC
        CALL NBRSTR(VAL,STR,NS)
        BUFFER='('//STR(1:NS)//')'
        NS = NS + 2
        STR=BUFFER(1:NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBAR
        ICOL = ICOL + 1
        VAL = TDELOC
        CALL NBRSTR(VAL,STR,NS)
        BUFFER='('//STR(1:NS)//')'
        NS = NS + 2
        STR=BUFFER(1:NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        WRITE(LFNPR,1775) (TABROW(IP),IFLAG(IP),IP=1,NCOL+1),
     +        TABROW(NCOL+2)
      ENDIF
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL)
      ELSE
        WRITE(LFNPR,1260) (DASHES,IP=1,NCOL+2)
      ENDIF
C  Write the total NLMO contribution for these columns
      ICOL = 0
      DO IB = IB1,IB2
        ICOL = ICOL + 1
        VAL = DI_L(IBIG(IB))+REPOL(IBIG(IB))+DELOC(IBIG(IB))
        CALL NBRSTR(VAL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBLNK
      ENDDO
      IF(ILAST.EQ.0) THEN
        WRITE(LFNPR,1800) (TABROW(IP),IFLAG(ICOL),IP=1,NCOL)
      ELSE
        ICOL = ICOL + 1
        VAL = TOTALC - (SUMBL + TBREPOL + TBDELOC)
        CALL NBRSTR(VAL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        IFLAG(ICOL) = IBAR
        ICOL = ICOL + 1
        VAL = TOTALC
        CALL NBRSTR(VAL,STR,NS)
        IF(NS.LT.8) THEN
          TABROW(ICOL) = BLANK(1:8-NS)//STR(1:NS)
        ELSE IF(NS.EQ.8) THEN
          TABROW(ICOL) = STR(1:NS)
        ELSE
          TABROW(ICOL) = '********'
        ENDIF
        WRITE(LFNPR,1800) (TABROW(IP),IFLAG(IP),IP=1,NCOL+1),
     +        TABROW(NCOL+2)
      ENDIF
      WRITE(LFNPR,1900)
      IF(ILAST.EQ.0) GOTO 60
      RETURN
 1000 FORMAT(/,1X,50A1)
 1150 FORMAT(5X,'Threshold for printing:',F7.2,' Hz')
 1200 FORMAT(12X,6(I10,'.'))
 1250 FORMAT(15X,'|',5(10A1,' '),10A1)
 1260 FORMAT('----------------',5A11,A8)
 1300 FORMAT(9X,'Lewis |',5(A8,1X,A1,1X))
 1305 FORMAT(9X,'Lewis |',5(A8,1X,A1,1X),A8)
 1310 FORMAT(9X,'Repol.|',5(A8,1X,A1,1X))
 1315 FORMAT(9X,'Repol.|',5(A8,1X,A1,1X),A8)
 1700 FORMAT(I4,'.',10A1,'|',5(A8,1X,A1,1X),A8)
 1710 FORMAT(5X,10A1,'|',5(A8,1X,A1,1X),A8)
 1750 FORMAT(9X,'Others','|',5(A8,1X,A1,1X))
 1755 FORMAT(9X,'Others','|',5(A8,1X,A1,1X),A8)
 1770 FORMAT(1X,'(Total deloc.)|',5(1X,A8,A1,1X))
 1775 FORMAT(1X,'(Total deloc.)|',5(1X,A8,A1,1X),A8)
 1800 FORMAT(9X,' NLMO |',5(A8,1X,A1,1X),A8)
 1900 FORMAT(1X)
      END
C***********************************************************************
      SUBROUTINE NJCATM(IAT)
C***********************************************************************
C 03-Feb-01  FAW  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER*256 FILENM,STR*3
      PARAMETER(MAXFIL = 40)
      COMMON/NBNAME/FILENM,NFILE,IFILE(MAXFIL)
      NL = LENNB(FILENM)
      DO I = NL,1,-1
        IF(FILENM(I:I).EQ.'_') THEN
          STR = FILENM(I+1:NL)
          IAT = INTVAL(STR)
          RETURN
        ENDIF
      ENDDO
      CALL NBHALT('Bad FILENM to find atom number in NJCATM.')
      END
C***********************************************************************
      SUBROUTINE NJCLNK(IA,IB,LINKS,NLINKS,IERR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BDFIND
C
C  Find shortest bond-connected linkage, LINKS(I), I = 1,NLINKS,
C  for atoms IA, IB in the same molecular unit (up to 7 links)
C
      DIMENSION LINKS(12)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      DIMENSION N1A(6),N1B(6),N2A(2,36),N2B(2,36),
     +  N3A(3,216),N3B(3,216)
      IERR = 0
      N1TOA = 0
      N2TOA = 0
      N3TOA = 0
      N1TOB = 0
      N2TOB = 0
      N3TOB = 0
C  Same atom?
      IF(IA.EQ.IB) THEN
        NLINKS = 0
        RETURN
      ENDIF
C  Directly bonded 1-J?
      IF(BDFIND(IA,IB)) THEN
        NLINKS = 1
        RETURN
      ENDIF
C  Find 1st-neighbors to A:  N1A(I), I=1,N1TOA
      N1TOA = 0
      DO IAP = 1,NATOMS
        IF(BDFIND(IA,IAP).AND.IAP.NE.IA) THEN
          N1TOA = N1TOA + 1
          N1A(N1TOA) = IAP
C  Geminal 2-J?  (NLINKS = 2)
          IF(BDFIND(IAP,IB)) THEN
            NLINKS = 2
            LINKS(1) = IAP
            RETURN
          ENDIF
        ENDIF
      ENDDO
      IF(N1TOA.EQ.0) GOTO 100
C  Find 1st-neighbors to B:  N1B(I), I=1,N1TOB
      N1TOB = 0
      DO IBP = 1,NATOMS
        IF(BDFIND(IB,IBP).AND.IBP.NE.IB) THEN
          N1TOB = N1TOB + 1
          N1B(N1TOB) = IBP
        ENDIF
      ENDDO
      IF(N1TOB.EQ.0) GOTO 100
C  Vicinal 3-J?  (NLINKS = 3)
      DO IL = 1,N1TOA
        IAP = N1A(IL)
        DO ILP = 1,N1TOB
          IBP = N1B(ILP)
          IF(BDFIND(IAP,IBP)) THEN
            NLINKS = 3
            LINKS(1) = IAP
            LINKS(2) = IBP
            RETURN
          ENDIF
        ENDDO
      ENDDO
C  Find 2nd-neighbors to A:  N2A(I), I=1,N2TOA
      N2TOA = 0
      DO IL = 1,N1TOA
        ILA = N1A(IL)
        DO ILP = 1,NATOMS
          IF(ILP.NE.ILA.AND.BDFIND(ILP,ILA)) THEN
            N2TOA = N2TOA + 1
            N2A(1,N2TOA) = ILA
            N2A(2,N2TOA) = ILP
          ENDIF
        ENDDO
      ENDDO
      IF(N2TOA.EQ.0) GOTO 100
C  4-bond coupling?
      DO IL = 1,N2TOA
        IAP = N2A(2,IL)
        DO ILP = 1,N1TOB
          IBP = N1B(ILP)
          IF(IAP.NE.IBP.AND.BDFIND(IAP,IBP)) THEN
            NLINKS = 4
            LINKS(1) = N2A(1,IL)
            LINKS(2) = N2A(2,IL)
            LINKS(3) = IBP
            RETURN
          ENDIF
        ENDDO
      ENDDO
C  Find 2nd-neighbors to B:  N2B(I), I=1,N2TOB
      N2TOB = 0
      DO IL = 1,N1TOB
        ILB = N1B(IL)
        DO ILP = 1,NATOMS
          IF(ILP.NE.ILB.AND.BDFIND(ILP,ILB)) THEN
            N2TOB = N2TOB + 1
            N2B(1,N2TOB) = ILB
            N2B(2,N2TOB) = ILP
          ENDIF
        ENDDO
      ENDDO
      IF(N2TOB.EQ.0) GOTO 100
C  5-bond coupling?
      DO IL = 1,N2TOA
        IAP = N2A(2,IL)
        DO ILP = 1,N2TOB
          IBP = N2B(2,ILP)
          IF(IAP.NE.IBP.AND.BDFIND(IAP,IBP)) THEN
            NLINKS = 5
            LINKS(1) = N2A(1,IL)
            LINKS(2) = IAP
            LINKS(3) = IBP
            LINKS(4) = N2B(1,ILP)
            RETURN
          ENDIF
        ENDDO
      ENDDO
C  Find 3rd-neighbors to A:  N3A(I), I=1,N3TOA
      N3TOA = 0
      DO IL = 1,N2TOA
        ILA = N2A(2,IL)
        DO ILP = 1,NATOMS
          IF(ILP.NE.ILA.AND.ILP.NE.N2A(1,IL).AND.BDFIND(ILP,ILA)) THEN
            N3TOA = N3TOA + 1
            N3A(1,N3TOA) = N2A(1,IL)
            N3A(2,N3TOA) = N2A(2,IL)
            N3A(3,N3TOA) = ILP
          ENDIF
        ENDDO
      ENDDO
      IF(N3TOA.EQ.0) GOTO 100
C  6-bond coupling?
      DO IL = 1,N3TOA
        IAP = N3A(3,IL)
        DO ILP = 1,N2TOB
          IBP = N2B(2,ILP)
          IF(IAP.NE.IBP.AND.BDFIND(IAP,IBP)) THEN
            NLINKS = 6
            LINKS(1) = N3A(1,IL)
            LINKS(2) = N3A(2,IL)
            LINKS(3) = IAP
            LINKS(4) = IBP
            LINKS(5) = N2B(1,ILP)
            RETURN
          ENDIF
        ENDDO
      ENDDO
C  Find 3rd-neighbors to B:  N3B(I), I=1,N3TOB
      N3TOB = 0
      DO IL = 1,N2TOB
        ILB = N2B(2,IL)
        DO ILP = 1,NATOMS
          IF(ILP.NE.ILB.AND.ILP.NE.N2B(1,IL).AND.BDFIND(ILP,ILB)) THEN
            N3TOB = N3TOB + 1
            N3B(1,N3TOB) = N2B(1,IL)
            N3B(2,N3TOB) = N2B(2,IL)
            N3B(3,N3TOB) = ILP
          ENDIF
        ENDDO
      ENDDO
      IF(N3TOB.EQ.0) GOTO 100
C  7-bond coupling?
      DO IL = 1,N3TOA
        IAP = N3A(3,IL)
        DO ILP = 1,N3TOB
          IBP = N3B(3,ILP)
          IF(IAP.NE.IBP.AND.BDFIND(IAP,IBP)) THEN
            NLINKS = 7
            LINKS(1) = N3A(1,IL)
            LINKS(2) = N3A(2,IL)
            LINKS(3) = IAP
            LINKS(4) = IBP
            LINKS(5) = N3B(2,ILP)
            LINKS(6) = N3B(1,ILP)
            RETURN
          ENDIF
        ENDDO
      ENDDO
  100 CONTINUE
C  Error condition encountered
      IERR=1
      RETURN
      END
C***********************************************************************
      SUBROUTINE NJCLAB(IATN,IATNP,LINKS,NLINKS,MOLLOC,TITLE)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER TITLE*(*),CHARAT*2,STR*10
      PARAMETER(MAXATM = 200)
      DIMENSION LINKS(12),MOLLOC(MAXATM)
      COMMON/NBATOM/IATNO(MAXATM),INO(MAXATM),NORB(MAXATM),LL(MAXATM),
     +       LU(MAXATM),IZNUC(MAXATM),IATCR(MAXATM)
      CALL NBISTR(NLINKS,STR,NS)
      IF(MOLLOC(IATN).NE.MOLLOC(IATNP)) THEN
        TITLE=STR(1:NS)//'h-J['
        NT = NS + 4
      ELSE
        TITLE=STR(1:NS)//'-J['
        NT = NS + 3
      ENDIF
      TITLE = TITLE(1:NT)//CHARAT(IATNO(IATN))
      NT = NT + 2
      CALL NBISTR(IATN,STR,NS)
      TITLE = TITLE(1:NT)//STR(1:NS)
      NT = NT + NS
      IALM = IATN
      IF(NLINKS.EQ.0) THEN
        TITLE = TITLE(1:NT)//']'
        NT = NT + 1
        GOTO 100
      ELSE IF(NLINKS.EQ.1) THEN
        CALL NBISTR(IATNP,STR,NS)
        IF(MOLLOC(IATN).EQ.MOLLOC(IATNP)) THEN
          TITLE=TITLE(1:NT)//'-'//CHARAT(IATNO(IATNP))//STR(1:NS)//']'
        ELSE
          TITLE=TITLE(1:NT)//':'//CHARAT(IATNO(IATNP))//STR(1:NS)//']'
        ENDIF
        NT = NT + NS + 4
        GOTO 100
      ENDIF
      DO IL = 1,NLINKS - 1
        IAL = LINKS(IL)
        CALL NBISTR(IAL,STR,NS)
        IZL = IATNO(IAL)
        IF(MOLLOC(IAL).EQ.MOLLOC(IALM)) THEN
          TITLE = TITLE(1:NT)//'-'//CHARAT(IZL)//STR(1:NS)
        ELSE
          TITLE = TITLE(1:NT)//':'//CHARAT(IZL)//STR(1:NS)
        ENDIF
        NT = NT + NS + 3
        IALM = IAL
      ENDDO
      CALL NBISTR(IATNP,STR,NS)
      IF(MOLLOC(IATNP).EQ.MOLLOC(IALM)) THEN
        TITLE = TITLE(1:NT)//'-'//CHARAT(IATNO(IATNP))//STR(1:NS)//']'
      ELSE
        TITLE = TITLE(1:NT)//':'//CHARAT(IATNO(IATNP))//STR(1:NS)//']'
      ENDIF
      NT = NT + NS + 4
  100 CONTINUE
C  Remove blanks
      I = 0
      DO IT = 1,NT
        IF(TITLE(IT:IT).NE.' ') THEN
          I = I + 1
          IF(I.NE.IT) TITLE(I:I) = TITLE(IT:IT)
        ENDIF
      ENDDO
      TITLE = TITLE(1:I)
      RETURN
      END
C***********************************************************************
      SUBROUTINE FNDHYBS(INBO,IH1,IH2,IH3)
C***********************************************************************
C 31-Jan-01  FAW  New subroutine
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
C Find the NHOs (IH1, IH2, IH3) for input NBO (INBO)
      PARAMETER(MAXBAS = 2000)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),LBL0(MAXBAS),
     +       LORBC(MAXBAS),LORB(MAXBAS)
      COMMON/NBLBL/NLEW,NVAL,LBL(10,MAXBAS,5)
      COMMON/NBINFO/ISPIN,NATOMS,NDIM,NBAS,NNAO,MXBO,MXAO,MXAOLM,MUNIT
      DIMENSION IHOLST(MAXBAS)
      DATA ISW/0/ICR/2HCR/ILP/2HLP/IBD/2HBD/I3C/2H3C/IRY/2HRY/
      SAVE ISW,ICR,ILP,IBD,I3C,IRY
C  Load the IHOLST vector with the lead NHO for each NBO
      IF(ISW.EQ.0) THEN
        ISW = 1
        IH = 0
        DO IB = 1,NLEW
          IBO = IBXM(IB)
          LIBO = LABEL(IBO,1)
          IH1 = IH + 1
          IF(LIBO.EQ.ICR.OR.LIBO.EQ.ILP) THEN
            IH = IH + 1
          ELSE IF(LIBO.EQ.IBD) THEN
            IH = IH + 2
C  Search for corresponding antibond
            DO 10 IBP = NLEW+1,NNAO
              IBOP = IBXM(IBP)
              DO K = 1,6
                IF(K.NE.2.AND.LABEL(IBOP,K).NE.LABEL(IBO,K)) GOTO 10
              ENDDO
              IHOLST(IBP) = IH1
   10       CONTINUE
C  Antibonds not implemented for 3-center bonds
          ELSE IF(LIBO.EQ.I3C) THEN
            IH = IH + 3
          ENDIF
          IHOLST(IB) = IH1
        ENDDO
        DO IB = NLEW+1,NNAO
          IBO = IBXM(IB)
          IF(IHOLST(IBO).EQ.0.AND.LABEL(IBO,1).EQ.IRY) THEN
            IH1 = IH + 1
            IHOLST(IBO) = IH1
          ENDIF
        ENDDO
      ENDIF
C  Look up the lead NHO for this NBO
      IH1 = IHOLST(INBO)
      IH2 = 0
      IH3 = 0
      IF(LABEL(INBO,1).EQ.IBD) THEN
        IH2 = IH1 + 1
      ELSE IF(LABEL(INBO,1).EQ.I3C) THEN
        IH2 = IH1 + 1
        IH3 = IH2 + 1
      ENDIF
      RETURN
      END
C***********************************************************************
      SUBROUTINE CONSTFC(CGYRO,IZ1,IZ2,IPR,LFNPR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER STR1*3,STR2*3,TITLE1*20,TITLE2*20,CHARAT*2
      DIMENSION GAMMA(103),ISO(103)
C Table of gyromagnetic ratios (elements 1-103).  Value is set
C to "1.0" when no experimental value is available.  (All values
C should be multiplied by 10^7 to get the actual gyromagnetic
C ratio for the formula below.)  Values from www.webelements.com
      DATA GAMMA/
     +  26.7522, -20.3802,  10.3977,  -3.7597,   2.8747,   6.7283,
     +  -2.7126,  -3.6281,  25.1815,  -2.1131,   7.0808,  -1.6389,
     +   6.9763,  -5.3190,  10.8394,   2.0557,   2.6242,   1.0000,
     +   1.2501,  -1.8031,   6.5088,  -1.5105,   7.0455,  -1.5152,
     +   6.6453,   0.8681,   6.3320,  -2.3948,   7.1118,   1.6767,
     +   6.4389,  -0.9360,   4.5962,   5.1254,   6.7256,  -1.0331,
     +   2.5927,  -1.1639,  -1.3163,  -2.4974,   6.5674,  -1.7880,
     +   6.0460,  -1.3770,  -0.8468,  -1.2300,  -1.0889,  -5.6983,
     +   5.8972, -10.0317,   6.4435,  -8.5108,   5.3896,  -7.4521,
     +   3.5333,   2.9930,   3.8083,   0.2189,   8.1907,  -1.4570,
     +   3.6130,  -1.1150,   2.9369,  -1.0769,   6.4310,   1.2890,
     +   5.7100,  -0.7716,  -2.2180,  -1.3025,   3.0552,   1.0860,
     +   3.2438,   1.1282,   6.1057,   2.1071,   0.5227,   5.8385,
     +   0.4731,   4.8458,  15.6922,   5.5805,   4.3750,   7.4000,
     +   1.0000,   1.0000,   1.0000,   1.0000,   3.5000,   0.4000,
     +   3.2100,  -0.5200,   3.1000,   0.9720,   1.4000,   0.2000,
     +   1.0000,   1.0000,   1.0000,   1.0000,   1.0000,   1.0000,
     +   1.0000/
C Table of "standard" isotopes (most abundant odd-spin nuclei).
C This is the mass number of the isotope for which the gyromagnetic
C ratio is given in GAMMA.
      DATA ISO/
     +  1,  3,  7,  9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31,
     + 33, 35, 39, 39, 43, 45, 47, 51, 53, 55, 57, 59, 61, 63, 67,
     + 69, 73, 75, 77, 79, 83, 85, 87, 89, 91, 93, 95, 99,101,103,
     +105,107,111,115,119,121,125,127,129,133,137,139,141,141,143,
     +147,147,153,157,159,163,165,167,169,173,175,177,181,183,187,
     +189,193,195,197,199,205,207,209,209,  0,  0,  0,  0,227,229,
     +231,235,237,239,241,247,247,251,  0,  0,  0,  0,  0/
      DATA ZERO/0.0D0/ONE/1.0D0/FOUR/4.0D0/D7/1.D7/
      SAVE ZERO,ONE,FOUR,D7,COEF,PI,PERT
      DATA COEF/4.17787445134D-13/PI/3.14159265357979D0/PERT/0.02D0/
C here's the equation to translate spin density into Hz:
C
C CONST = CONSTFC * gamma1 * gamma2
C
C where
C
C CONSTFC = 4.17787445134 * 10**-13 * (1/(4*PI*PI)) * (1 / PERT);
C
C and PERT is the perturbation (e.g. PERT = .02 for a perturbation of
C 200 in the gaussian input file).
      IF(IZ1.LT.0.OR.IZ1.GT.103.OR.IZ2.LT.0.OR.IZ2.GT.103) THEN
        WRITE(LFNPR,1000) IZ1,IZ2
 1000   FORMAT(1X,'*** Illegal atomic number(s):',2I5)
        G1G2 = ZERO
      ELSE
        G1 = GAMMA(IZ1)*D7
        G2 = GAMMA(IZ2)*D7
        G1G2 = G1*G2
        IF(IPR.GT.0) THEN
          CALL NBISTR(ISO(IZ1),STR1,NS1)
          TITLE1='gamma('//STR1(1:NS1)//'-'//CHARAT(IZ1)//')'
          NT1 = LENNB(TITLE1)
          CALL RMBLNK(TITLE1,NT1)
          CALL NBISTR(ISO(IZ2),STR2,NS2)
          TITLE2='gamma('//STR2(1:NS2)//'-'//CHARAT(IZ2)//')'
          NT2 = LENNB(TITLE2)
          CALL RMBLNK(TITLE2,NT2)
          IF(IZ1.NE.IZ2) THEN
            WRITE(LFNPR,1100) TITLE1(1:NT1),G1,TITLE2(1:NT2),G2
 1100       FORMAT(1X,A13,' =',D12.4,'; ',A13,' =',D12.4)
          ELSE
            WRITE(LFNPR,1200) TITLE1(1:NT1),G1
 1200       FORMAT(1X,A13,' =',D12.4)
          ENDIF
        ENDIF
      ENDIF
      CGYRO = COEF*(ONE/(FOUR*PI*PI*PERT))*G1G2
      RETURN
      END
C***********************************************************************
      SUBROUTINE RMBLNK(STR,NS)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER STR*(*)
C  Remove blank characters from STR(1:NS)
      IF(NS.LE.0) RETURN
      I = 0
      DO IT = 1,NS
        IF(STR(IT:IT).NE.' ') THEN
          I = I + 1
          IF(I.NE.IT) STR(I:I) = STR(IT:IT)
        ENDIF
      ENDDO
      STR = STR(1:I)
      NS = I
      RETURN
      END
C***********************************************************************
      FUNCTION IPHTYP(IBO,JBO)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      LOGICAL BDFIND
C
      PARAMETER(MAXBAS = 2000)
      COMMON/NBBAS/LABEL(MAXBAS,6),NBOUNI(MAXBAS),NBOTYP(MAXBAS),
     +       LSTOCC(MAXBAS),IBXM(MAXBAS),LARC(MAXBAS),IATHY(MAXBAS,3)
C
      SAVE IV,IG,IR,II
      DATA IV,IG,IR,II/1Hv,1Hg,1Hr,1Hi/
C
C  Determine whether the IBO->JBO delocalization is vicinal (IPHTYP='v'),
C  geminal (IPHTYP='g'), remote (IPHTYP='r'), or internal (IPHTYP='i'):
C  The relationship is "internal" if atoms of JBO are a proper subset
C  of those of IBO
C
      IPHTYP = IR
      IF(NBOUNI(IBO).EQ.NBOUNI(JBO)) THEN
        ICTR = MOD(NBOTYP(IBO),10)
        IB   = IBXM(IBO)
        JCTR = MOD(NBOTYP(JBO),10)
        JB   = IBXM(JBO)
        DO 20 I = 1,ICTR
          IAT = LABEL(IB,I+3)
          DO 10 J = 1,JCTR
            JAT = LABEL(JB,J+3)
            IF(IAT.EQ.JAT) THEN
              IPHTYP = IG
              GOTO 25
            ELSE IF(BDFIND(IAT,JAT)) THEN
              IPHTYP = IV
            END IF
   10     CONTINUE
   20   CONTINUE
C Are atoms of JBO a proper subset of those in IBO? (internal code 'i')
   25   CONTINUE
        DO J = 1,JCTR
          IMATCH = 0
          JAT = LABEL(JB,J+3)
          DO I = 1,ICTR
            IAT = LABEL(IB,I+3)
            IF(JAT.EQ.IAT) IMATCH = 1
          ENDDO
          IF(IMATCH.EQ.0) GOTO 30
        ENDDO
        IPHTYP = II
   30   CONTINUE
      END IF
      RETURN
      END
C***********************************************************************
      SUBROUTINE FCAOEL(FCAO,SCR,ISCR,NBAS,NDIM,NATOMS)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  Evaluate the array of Fermi-contact amplitudes of each basis AO (j)
C  at each nuclear center (N): FCAO(j,N).
C  (Note that scratch arrays SCR, ISCR share the same memory)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      COMMON/NBAO/LCTR(MAXBAS),LANG(MAXBAS)
      DIMENSION FCAO(NDIM,NATOMS),SCR(*),ISCR(*)
      DIMENSION AOAMP(MAXBAS),R(3),ORIGIN(3)
C Set MAXSHL = MAXEXP = 3.6 * MAXATM if re-dimensioned
      PARAMETER (MAXSHL=720,MAXEXP=720,MAXCMP=20)
C  Cf. COMMON/AOINFO/ of ORBPLOT
      COMMON/NBORBP/NATMS,NSHELL,NEXP,NORB,IAN(MAXATM),ATCRD(3,MAXATM),
     +              JCTR(MAXSHL),JNCMP(MAXSHL),JPTR(MAXSHL),JNG(MAXSHL),
     +              JORBS(MAXCMP,MAXSHL),EXX(MAXEXP),CSPDF(MAXEXP,4),
     +              ONORM(MAXBAS)
      DATA ISW/0/ZERO/0.0D0/
      IF(ISW.NE.0) RETURN
      NATMS = NATOMS
C-----------------------------------------------------------------------
C  Load the AO info into COMMON/NBORBP/ (a la WRBAS)
      LFN = 99
      OPEN(LFN,STATUS='SCRATCH')
C-----------------------------------------------------------------------
C  Fetch the number of shells NSHELL, the number of exponents NEXP,
C  the NCOMP, NPRIM, and NPTR arrays, and the orbital exponents and
C  coefficients from the NBO DAF:
C
      CALL FEBAS(NSHELL,NEXP,ISCR)
C
C  Partition the scratch arrays:  (Note that SCR and ISCR occupy the same
C  space in memory)
C
C  ISCR: (integer)
C
C   NSHELL NEXP  NCOMP  NPRIM  NPTR
C  +------+-----+------+------+-----+-----------------------------------
C               I1     I2     I3
C
C  SCR: (real)
C
C                                     EXP   CS   CP   CD   CF   CG
C  ---------------------------------+-----+----+----+----+----+----+
C                                   I4    I5   I6   I7   I8   I9
C
C  ISCR(I1) : NCOMP(1..NSHELL)
C  ISCR(I2) : NPRIM(1..NSHELL)
C  ISCR(I3) : NPTR(1..NSHELL)
C  SCR(I4)  : EXP(1..NEXP)
C  SCR(I5)  : CS(1..NEXP)
C  SCR(I6)  : CP(1..NEXP)
C  SCR(I7)  : CD(1..NEXP)
C  SCR(I8)  : CF(1..NEXP)
C  SCR(I9)  : CG(1..NEXP)
C  SCR(I10) : TITLE(10) or ATCOOR(3*NATOMS)
C
      I1   = 3
      I2   = I1 + NSHELL
      I3   = I2 + NSHELL
      I4   = I3 + NSHELL
      I5   = I4 + NEXP
      I6   = I5 + NEXP
      I7   = I6 + NEXP
      I8   = I7 + NEXP
      I9   = I8 + NEXP
C     IEND = I9 + NEXP
C
C  Fetch the atomic coordinates:
C
      CALL FECOOR(ATCRD)
C
C  Write out information about each shell in the basis set:
C
C     NCTR(I)  --  atomic center of the Ith shell
C
C     NCOMP(I) --  number of components in the Ith shell
C
C     NPTR(I)  --  pointer for the Ith shell into the primitive parameters
C                  of EXP, CS, CP, CD, CF, and CG
C
C     NPRIM(I) --  number of primitive functions in the Ith shell
C
C     LABEL(1..NCOMP(I)) -- symmetry labels for the orbitals of this shell
C
      J1 = 1
      J2 = I1
      J3 = I3
      J4 = I2
      DO 20 I = 1,NSHELL
        NCOMP = ISCR(J2)
        NPRIM = ISCR(J3)
        NPTR  = ISCR(J4)
        WRITE(LFN,950) LCTR(J1),NCOMP,NPRIM,NPTR
        WRITE(LFN,950) (LANG(J1+J),J=0,NCOMP-1)
        J1 = J1 + NCOMP
        J2 = J2 + 1
        J3 = J3 + 1
        J4 = J4 + 1
   20 CONTINUE
C
C  Write out the primitive parameters:
C
      WRITE(LFN,960) (SCR(I4+I),I=0,NEXP-1)
      WRITE(LFN,960) (SCR(I5+I),I=0,NEXP-1)
      WRITE(LFN,960) (SCR(I6+I),I=0,NEXP-1)
      WRITE(LFN,960) (SCR(I7+I),I=0,NEXP-1)
      WRITE(LFN,960) (SCR(I8+I),I=0,NEXP-1)
      WRITE(LFN,960) (SCR(I9+I),I=0,NEXP-1)
  950 FORMAT(1X,10I6)
  960 FORMAT(2X,4E18.9)
C-----------------------------------------------------------------------
C  Code from AOIN (ORBPLOT) to read WRBAS info and write to COMMON/NBORBP/
C-----------------------------------------------------------------------
C
      REWIND(LFN)
      NORB = 0
      DO I = 1,NSHELL
        READ (LFN,*) JCTR(I),JNCMP(I),JPTR(I),JNG(I)
        READ (LFN,*) (JORBS(J,I),J=1,JNCMP(I))
        NORB = NORB + JNCMP(I)
      END DO
C
      READ (LFN,*) (EXX(K),K=1,NEXP)
      READ (LFN,*) (CSPDF(K,1),K=1,NEXP)
      READ (LFN,*) (CSPDF(K,2),K=1,NEXP)
      READ (LFN,*) (CSPDF(K,3),K=1,NEXP)
      READ (LFN,*,ERR=500) (CSPDF(K,4),K=1,NEXP)
  500 CONTINUE
      CALL BNORM
      CLOSE(LFN)
      NBAS=NORB
C-----------------------------------------------------------------------
C  Loop over nuclear positions
      R(1) = ZERO
      R(2) = ZERO
      R(3) = ZERO
      DO IATN = 1,NATOMS
        ORIGIN(1) = ATCRD(1,IATN)
        ORIGIN(2) = ATCRD(2,IATN)
        ORIGIN(3) = ATCRD(3,IATN)
C  Evaluate the amplitude AOAMP(IBAS) at each basis AO
        CALL AOVAL2(AOAMP,R,ORIGIN)
        DO IJ = 1,NBAS
          FCAO(IJ,IATN) = AOAMP(IJ)
        ENDDO
      ENDDO
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBRSTR(FVAL,STR,NS)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER STR*(*),FISTR*10,FDSTR*10
C Convert real FVAL to STR(1:NS) (nearest 0.01)
      DATA ZERO/0.0D0/
      INTF=INT(FVAL)
      CALL NBISTR(INTF,FISTR,NIS)
      DECF=ABS(FVAL-INTF)
      IDEC=INT(1.0D2*DECF+0.49D0)
      IF(IDEC.GT.99) IDEC=99
      CALL NBISTR(IDEC,FDSTR,NDS)
      IF(NDS.GE.2) THEN
        FDSTR=FDSTR(1:2)
      ELSE IF(NDS.EQ.1) THEN
        FDSTR='0'//FDSTR(1:NDS)
      ELSE IF(NDS.EQ.0) THEN
        FDSTR='00'
      ENDIF
      NDS=2
C   10 CONTINUE
C      IF(FDSTR(NDS:NDS).EQ.'0') THEN
C        NDS=NDS-1
C        IF(NDS.GT.3) GOTO 10
C      ENDIF
      IF(FVAL.LT.ZERO.AND.INTF.EQ.0) THEN
        STR='-'//FISTR(1:NIS)//'.'//FDSTR(1:NDS)
        NS=NIS+NDS+2
      ELSE
        STR=FISTR(1:NIS)//'.'//FDSTR(1:NDS)
        NS=NIS+NDS+1
      ENDIF
      RETURN
      END
C***********************************************************************
      SUBROUTINE NBISTR(INT0,STR,NS)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C  Convert interger to STR(1:NS)
      CHARACTER STR*(*),C*1,STRP*50
      IF(INT0.LT.0) THEN
        ISIGN=-1
        INT=-INT0
      ELSE
        ISIGN=1
        INT=INT0
      ENDIF
      NS=0
   10 CONTINUE
      IF(INT.GE.10) THEN
        INTD=INT/10
        IF(INTD.GT.0) THEN
          IREM=INT-10*INTD
          C=CHAR(IREM+48)
          IF(NS.EQ.0) THEN
            STRP=C
          ELSE
            STRP=STRP(1:NS)//C
          ENDIF
          NS=NS+1
          INT=INTD
          GOTO 10
        ENDIF
      ELSE
        C=CHAR(INT+48)
        IF(NS.EQ.0) THEN
          STR=C
          STRP=C
        ELSE
          STRP=STRP(1:NS)//C
        ENDIF
        NS=NS+1
      ENDIF
      IF(ISIGN.LT.0) THEN
        STRP=STRP(1:NS)//'-'
        NS=NS+1
      ENDIF
      IF(NS.GT.1) THEN
        DO I=1,NS
          J=NS+1-I
          STR(J:J)=STRP(I:I)
        ENDDO
      ENDIF
      RETURN
      END
C***********************************************************************
      FUNCTION INTVAL(STR)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      CHARACTER STR*(*)
      DIMENSION IPOW(7)
C CONVERTS A STRING VARIABLE TO POSITIVE INTEGER VALUE (< 10,000,000)
      DATA IPOW/1,10,100,1000,10000,100000,1000000/
      NL=LENNB(STR)
      ITMP=0
      IMAG=0
   10 IF(NL.LE.0) THEN
        INTVAL = ITMP
        RETURN
      ELSE
        IC=ICHAR(STR(NL:NL))
        I=IC-48
        IF(I.LT.0.OR.I.GT.9) THEN
          WRITE(6,1000)
 1000     FORMAT(1X,'*** Illegal INTVAL string variable')
          INTVAL = 0
          RETURN
        ENDIF
        IMAG=IMAG+1
        ITMP=ITMP+I*IPOW(IMAG)
        NL=NL-1
        GOTO 10
      ENDIF
      END
C***********************************************************************
      SUBROUTINE AOVAL2(AOAMP,R,ORIGIN)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
      PARAMETER (TOL=1.0D-8,CVTOAU=0.529177D0)
C Set MAXSHL = MAXEXP = 3.6 * MAXATM if re-dimensioned
      PARAMETER (MAXSHL=720,MAXEXP=720,MAXCMP=20)
C
      DIMENSION R(3),RD(3),IAOS(MAXSHL),ORIGIN(3),AOAMP(MAXBAS)
      COMMON/NBORBP/NATOMS,NSHELL,NEXP,NORB,IAN(MAXATM),ATCRD(3,MAXATM),
     +              JCTR(MAXSHL),JNCMP(MAXSHL),JPTR(MAXSHL),JNG(MAXSHL),
     +              JORBS(MAXCMP,MAXSHL),EXX(MAXEXP),CSPDF(MAXEXP,4),
     +              ONORM(MAXBAS)
C
      IAOS(1) = 1
      DO 50 I = 2,NSHELL
        IAOS(I) = IAOS(I-1) + JNCMP(I-1)
   50 CONTINUE
C
      DO 60 I = 1,NORB
         AOAMP(I) = 0.0D0
   60 CONTINUE
C
      DO 100 I = 1,NSHELL
        IPTR = JPTR(I) - 1
        NG   = JNG(I)
        NCMP = JNCMP(I)
        RD(1) = (R(1) - ATCRD(1,JCTR(I)) + ORIGIN(1))/CVTOAU
        RD(2) = (R(2) - ATCRD(2,JCTR(I)) + ORIGIN(2))/CVTOAU
        RD(3) = (R(3) - ATCRD(3,JCTR(I)) + ORIGIN(3))/CVTOAU
        DO 90 K = 1,NG
          EP   = EXP(-EXX(K+IPTR)*(RD(1)*RD(1)+RD(2)*RD(2)+RD(3)*RD(3)))
          IF (EP.GT.TOL) THEN
            IORB  = IAOS(I)
            DO 80 L = 1,NCMP
              ITYPE = JORBS(L,I)
              IT    = ITYPE/100 + 1
              C1    = CSPDF(K+IPTR,IT)
              AOAMP(IORB) = AOAMP(IORB) + ONORM(IORB)*C1*
     +                                    AOAFAC(RD,ITYPE)*EP
              IORB  = IORB + 1
   80       CONTINUE
          END IF
   90   CONTINUE
  100 CONTINUE
C
      RETURN
      END
C***********************************************************************
      FUNCTION AOAFAC(R,ITYPE)
C***********************************************************************
C
C  Return the value of the atomic orbital structure factor.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(3),CAL(3),N(3)
      DIMENSION NXAL(3),NYAL(3),NZAL(3)
      DATA ZERO/0.0D0/
      SAVE ZERO
C
      IDUM = ITYPE
      CALL STRUCF(IDUM,NKAL,CAL,NXAL,NYAL,NZAL)
C
      TMP = ZERO
      DO I = 1,NKAL
        N(1) = NXAL(I)
        N(2) = NYAL(I)
        N(3) = NZAL(I)
        TMP = TMP+CAL(I)*ATMFAC(R,N)
      END DO
      AOAFAC = TMP
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE ITYPERR(ITYPE)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      WRITE (6,900) ITYPE
  900 FORMAT (2X,'ITYPE = ',I6)
      CALL NBHALT('Illegal value of ITYPE encountered in STRUCF.')
      END
C***********************************************************************
      SUBROUTINE STRUCF(ITYPE,NKAL,CAL,NXAL,NYAL,NZAL)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C   This routine returns the parameters making up the structure factors
C  for the atomic orbital labeled by ITYPE.  Any atomic orbital can be
C  written as
C
C                         NKAL
C         ChiALJ(X,Y,Z) = SUM  dAL(K)(X,Y,Z)*EXP(-ZETAJ*R**2)
C                         K=1
C
C         where dAL(K)(X,Y,Z) is the Kth structure factor of the orbital
C                             labelled by alpha (A) and lambda (L).
C                             (the labels alpha and lambda are similar to
C                             the correct, but restrictive labels l and m)
C
C         dAL(K)(X,Y,Z) = CAL(K)*(X**NXAL(K))*(Y**NYAL(K))*(Z**NZAL(K))
C
C
      PARAMETER (MAXA=6,MAXL=10)  !MAXA and MAXL are limits on alpha and lambda
      PARAMETER (NUMSFS=44)
C
      DIMENSION CAL(1),NXAL(1),NYAL(1),NZAL(1)
C
      DIMENSION LSTNKAL(MAXL,MAXA),LSTPTR(MAXL,MAXA)
      DIMENSION CALLST(NUMSFS),LSTPOW(3,NUMSFS)
C
C  Data structure LSTNKAL contains the number of terms in each structure factor
C                         A (alpha) labels the column and L (lambda) the row.
C
C                    1  2  3  4  5  6  7  8  9 10     A (alpha)     subshell
C                    -  -  -  -  -  -  -  -  - --     ---------     --------
C
      DATA LSTNKAL / 1, 0, 0, 0, 0, 0, 0, 0, 0, 0,   !    1            S
     +               1, 1, 1, 0, 0, 0, 0, 0, 0, 0,   !    2            P
     +               1, 1, 1, 1, 1, 1, 0, 0, 0, 0,   !    3         cart D
     +               1, 1, 1, 2, 3, 0, 0, 0, 0, 0,   !    4         pure D
     +               1, 1, 1, 1, 1, 1, 1, 1, 1, 1,   !    5         cart F
     +               3, 3, 3, 2, 1, 2, 2, 0, 0, 0 /  !    6         pure F
C
C
C  Data structure LSTPTR contains pointers into the LSTCAL and LSTPOW arrays
C                        for each A (columns) and L (rows).
C
C                  1  2  3  4  5  6  7  8  9 10       A (alpha)     subshell
C                  -  -  -  -  -  -  -  -  - --       ---------     --------
C
      DATA LSTPTR /1, 0, 0, 0, 0, 0, 0, 0, 0, 0,     !    1            S
     +             2, 3, 4, 0, 0, 0, 0, 0, 0, 0,     !    2            P
     +             5, 6, 7, 8, 9, 10,0, 0, 0, 0,     !    3         cart D
     +            11,12,13,14,16, 0, 0, 0, 0, 0,     !    4         pure D
     +            19,20,21,22,23,24,25,26,27,28,     !    5         cart F
     +            29,32,35,38,40,41,43, 0, 0, 0 /    !    6         pure F
C
C
C  Data structure CALLST contains the coefficients CAL for each A and L
C
C                    CAL            K    A (alpha)   L (lambda)   orbital
C                    ---           ---   ---------   ----------   -------
C
      DATA CALLST / 1.0D0,       !  1        1           1           S
     +              1.0D0,       !  1        2           1           Px
     +              1.0D0,       !  1        2           2           Py
     +              1.0D0,       !  1        2           3           Pz
     +              1.0D0,       !  1        3           1           Dxx
     +              1.0D0,       !  1        3           2           Dxy
     +              1.0D0,       !  1        3           3           Dxz
     +              1.0D0,       !  1        3           4           Dyy
     +              1.0D0,       !  1        3           5           Dyz
     +              1.0D0,       !  1        3           6           Dzz
     +              1.0D0,       !  1        4           1           Dxy
     +              1.0D0,       !  1        4           2           Dxz
     +              1.0D0,       !  1        4           3           Dyz
     +              1.0D0,       !  1        4           4           Dx2-y2
     +             -1.0D0,       !  2        4           4
     +              2.0D0,       !  1        4           5           Dz2
     +             -1.0D0,       !  2        4           5
     +             -1.0D0,       !  3        4           5
     +              1.0D0,       !  1        5           1           Fxxx
     +              1.0D0,       !  1        5           2           Fxxy
     +              1.0D0,       !  1        5           3           Fxxz
     +              1.0D0,       !  1        5           4           Fxyy
     +              1.0D0,       !  1        5           5           Fxyz
     +              1.0D0,       !  1        5           6           Fxzz
     +              1.0D0,       !  1        5           7           Fyyy
     +              1.0D0,       !  1        5           8           Fyyz
     +              1.0D0,       !  1        5           9           Fyzz
     +              1.0D0,       !  1        5          10           Fzzz
     +              2.0D0,       !  1        6           1           F(0)
     +             -3.0D0,       !  2        6           1
     +             -3.0D0,       !  3        6           1
     +              4.0D0,       !  1        6           2           F(C1)
     +             -1.0D0,       !  2        6           2
     +             -1.0D0,       !  3        6           2
     +              4.0D0,       !  1        6           3           F(S1)
     +             -1.0D0,       !  2        6           3
     +             -1.0D0,       !  3        6           3
     +              1.0D0,       !  1        6           4           F(C2)
     +             -1.0D0,       !  2        6           4
     +              1.0D0,       !  1        6           5           F(S2)
     +              1.0D0,       !  1        6           6           F(C3)
     +             -3.0D0,       !  2        6           6
     +              3.0D0,       !  1        6           7           F(S3)
     +             -1.0D0 /      !  2        6           7
C
C   Data structure LSTPOW contains the powers NXAL(K), NYAL(K), NZAL(K)
C                         for each dAL(K)
C
C                   NX  NY  NZ        K    A (alpha)   L (lambda)   orbital
C                   --  --  --       ---   ---------   ----------   -------
C
      DATA LSTPOW / 0,  0,  0,      ! 1        1            1          S
     +              1,  0,  0,      ! 1        2            1          Px
     +              0,  1,  0,      ! 1        2            2          Py
     +              0,  0,  1,      ! 1        2            3          Pz
     +              2,  0,  0,      ! 1        3            1          Dxx
     +              1,  1,  0,      ! 1        3            2          Dxy
     +              1,  0,  1,      ! 1        3            3          Dxz
     +              0,  2,  0,      ! 1        3            4          Dyy
     +              0,  1,  1,      ! 1        3            5          Dyz
     +              0,  0,  2,      ! 1        3            6          Dzz
     +              1,  1,  0,      ! 1        4            1          Dxy
     +              1,  0,  1,      ! 1        4            2          Dxz
     +              0,  1,  1,      ! 1        4            3          Dyz
     +              2,  0,  0,      ! 1        4            4          Dx2-y2
     +              0,  2,  0,      ! 2        4            4
     +              0,  0,  2,      ! 1        4            5          Dz2
     +              2,  0,  0,      ! 2        4            5
     +              0,  2,  0,      ! 3        4            5
     +              3,  0,  0,      ! 1        5            1          Fxxx
     +              2,  1,  0,      ! 1        5            2          Fxxy
     +              2,  0,  1,      ! 1        5            3          Fxxz
     +              1,  2,  0,      ! 1        5            4          Fxyy
     +              1,  1,  1,      ! 1        5            5          Fxyz
     +              1,  0,  2,      ! 1        5            6          Fxzz
     +              0,  3,  0,      ! 1        5            7          Fyyy
     +              0,  2,  1,      ! 1        5            8          Fyyz
     +              0,  1,  2,      ! 1        5            9          Fyzz
     +              0,  0,  3,      ! 1        5           10          Fzzz
     +              0,  0,  3,      ! 1        6            1          F(0)
     +              2,  0,  1,      ! 2        6            1
     +              0,  2,  1,      ! 3        6            1
     +              1,  0,  2,      ! 1        6            2          F(C1)
     +              3,  0,  0,      ! 2        6            2
     +              1,  2,  0,      ! 3        6            2
     +              0,  1,  2,      ! 1        6            3          F(S1)
     +              2,  1,  0,      ! 2        6            3
     +              0,  3,  0,      ! 3        6            3
     +              2,  0,  1,      ! 1        6            4          F(C2)
     +              0,  2,  1,      ! 2        6            4
     +              1,  1,  1,      ! 1        6            5          F(S2)
     +              3,  0,  0,      ! 1        6            6          F(C3)
     +              1,  2,  0,      ! 2        6            6
     +              2,  1,  0,      ! 1        6            7          F(S3)
     +              0,  3,  0 /     ! 2        6            7
C
C  Begin execution
C
C   First, find the alpha (IALPHA) and lambda (LAMBDA) that correspond to
C  the value of ITYPE.  Check for incorrect values.
C
      IT     = ITYPE/10
      ITA   = IT*10
      LAMBDA = ITYPE - ITA
C
      IF ((ITA.EQ.0).OR.(ITA.EQ.50)) THEN
       IALPHA = 1    ! S orbital
       IF (LAMBDA.NE.1) CALL ITYPERR(ITYPE)
      ELSE IF ((ITA.EQ.100).OR.(ITA.EQ.150)) THEN
       IALPHA = 2    ! P orbital
       IF ((LAMBDA.LT.1).OR.(LAMBDA.GT.3)) CALL ITYPERR(ITYPE)
      ELSE IF (ITA.EQ.200) THEN
       IALPHA = 3    ! cartesian D orbital
       IF ((LAMBDA.LT.1).OR.(LAMBDA.GT.6)) CALL ITYPERR(ITYPE)
      ELSE IF (ITA.EQ.250) THEN
       IALPHA = 4    ! pure D orbital
       IF ((LAMBDA.LT.1).OR.(LAMBDA.GT.5)) CALL ITYPERR(ITYPE)
      ELSE IF (ITA.EQ.300) THEN
       IALPHA = 5    ! cartesian F orbital
       IF ((LAMBDA.LT.1).OR.(LAMBDA.GT.10)) CALL ITYPERR(ITYPE)
      ELSE IF (ITA.EQ.350) THEN
       IALPHA = 6    ! pure F orbital
      ELSE
       CALL ITYPERR(ITYPE)
      END IF
C
C  Now, using the values of alpha and lambda, look up the correct structure
C factor parameters, store them in the call arguments and return.
C
      NKAL = LSTNKAL(LAMBDA,IALPHA)
      IPTR = LSTPTR(LAMBDA,IALPHA) - 1
      DO 10 K=1,NKAL
       CAL(K)  = CALLST(IPTR+K)
       NXAL(K) = LSTPOW(1,IPTR+K)
       NYAL(K) = LSTPOW(2,IPTR+K)
       NZAL(K) = LSTPOW(3,IPTR+K)
   10 CONTINUE
C
      RETURN
      END
C***********************************************************************
      SUBROUTINE BNORM
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
      PARAMETER(MAXATM = 200,MAXBAS = 2000)
C Set MAXSHL = MAXEXP = 3.6 * MAXATM if re-dimensioned
      PARAMETER (MAXSHL=720,MAXEXP=720,MAXCMP=20)
C
      COMMON/NBORBP/NATOMS,NSHELL,NEXP,NORB,IAN(MAXATM),
     +                C(3,MAXATM),JCTR(MAXSHL),JNCMP(MAXSHL),
     +                JPTR(MAXSHL),JNG(MAXSHL),JORBS(MAXCMP,MAXSHL),
     +                EXX(MAXEXP),CSPDF(MAXEXP,4),ONORM(MAXBAS)
      DATA ZERO/0.0D0/ONE/1.0D0/
      SAVE ZERO,ONE
C
C
      IORB = 0
      DO I = 1,NSHELL
        IPTR = JPTR(I) - 1
        NG   = JNG(I)
        NCMP = JNCMP(I)
        DO J = 1,NCMP
          IORB  = IORB + 1
          ITYPE = JORBS(J,I)
          IT    = ITYPE/100 + 1
          SUM   = ZERO
          DO K = 1,NG
            KG    = IPTR + K
            ZETA1 = EXX(KG)
            C1    = CSPDF(KG,IT)
            DO L = 1,NG
              LG    = IPTR + L
              ZETA2 = EXX(LG)
              C2    = CSPDF(LG,IT)
              SUM = SUM + C1*C2*AONORM(ITYPE,ZETA1,ZETA2)
            END DO
          END DO
          ONORM(IORB) = ONE/SQRT(SUM)
        END DO
      END DO
C
      RETURN
      END
C***********************************************************************
      FUNCTION AONORM(ITYPE,ZETA1,ZETA2)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
C
C  This routine computes the overlap between two gaussian atomic orbitals
C located on the same center.  Each has the same angular symmetry specified
C by ITYPE, but different gaussian exponents (ZETA1 and ZETA2).
C
C  See the routine STRUCF for the definition of an atomic orbital in terms
C of geometric structure factors.
C
      DIMENSION CAL(3),NXAL(3),NYAL(3),NZAL(3)
      DIMENSION C1(3),C2(3),N1(3),N2(3)
C
      DATA C1/0.0D0,0.0D0,0.0D0/
      DATA C2/0.0D0,0.0D0,0.0D0/
C
C  Begin execution
C
C    Get structure factor information
C
      IDUM = ITYPE
      CALL STRUCF(IDUM,NKAL,CAL,NXAL,NYAL,NZAL)
C
C  Find value of overlap integral.
C
      TMP = 0.0D0
      DO K1 = 1,NKAL
        DO K2 = 1,NKAL
          N1(1) = NXAL(K1)
          N1(2) = NYAL(K1)
          N1(3) = NZAL(K1)
          N2(1) = NXAL(K2)
          N2(2) = NYAL(K2)
          N2(3) = NZAL(K2)
          TMP   = TMP + CAL(K1)*CAL(K2)
     +            * OLAPCG(N1,ZETA1,C1,N2,ZETA2,C2)
        END DO
      END DO
      AONORM = TMP
C
      RETURN
      END
C***********************************************************************
      FUNCTION ATMFAC(R,N)
C***********************************************************************
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION R(1),N(1)
C
      A = 1.0D0
      DO I = 1,3
        IF (N(I).GT.0) THEN
          IF (N(I).EQ.1) THEN
            A = A*R(I)
          ELSE IF (N(I).EQ.2) THEN
            A = A*R(I)*R(I)
          ELSE IF (N(I).EQ.3) THEN
            A = A*R(I)*R(I)*R(I)
          ELSE IF (N(I).GE.4) THEN
            A = A*R(I)**N(I)
          END IF
        END IF
      END DO
      ATMFAC = A
C
      RETURN
      END
C***********************************************************************
      FUNCTION OLAPCG(N1,ZETA1,C1,N2,ZETA2,C2)
C***********************************************************************
C
C   Cartesian gaussian overlap evaluator : Compute the overlap of
C   two cartesian gaussians, one described by X,Y,Z exponents in N1,
C   gaussian exponent ZETA1, and centered at C1 (in atomic units),
C   and a second gaussian described by N2, ZETA2, and C2.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      INTEGER N1(3),N2(3)
      DIMENSION C1(3),C2(3)
C
C
      ALFA = 1.0D0/(ZETA1+ZETA2)**0.5D0
      OTMP = 1.0D0
      DO I = 1,3
         TEMP = OLAP1D(N1(I),ZETA1,C1(I),N2(I),ZETA2,C2(I),ALFA)
         IF(TEMP.EQ.0.0D0) GOTO 10
         OTMP = OTMP*TEMP
      END DO
      OLAPCG = OTMP
      RETURN
C
   10 CONTINUE
      OLAPCG = 0.0D0
      RETURN
C
      END
C***********************************************************************
      FUNCTION OLAP1D(NX1,ZETA1,X1,NX2,ZETA2,X2,ALFA)
C***********************************************************************
C
C   Evaluate the 1-dimensional overlap of two gaussians.
C
C-----------------------------------------------------------------------
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION C(11),ETA1(6),ETA2(6)
      INTEGER BINOM(6,6)
C
C-----------------------------------------------------------------------
C   Table of binomial coefficients:
C  ---------------------------------
C                            (I-1)!
C            BINOM(I,J) = -------------            (1 <= J <= I)
C                         (I-J)! (J-1)!
C
C                       1       2       3       4       5       6     I/J
C                      ---     ---     ---     ---     ---     ---
      DATA BINOM /      1,      1,      1,      1,      1,      1,    !  1
     +                  0,      1,      2,      3,      4,      5,    !  2
     +                  0,      0,      1,      3,      6,     10,    !  3
     +                  0,      0,      0,      1,      4,     10,    !  4
     +                  0,      0,      0,      0,      1,      5,    !  5
     +                  0,      0,      0,      0,      0,      1/    !  6
C
C-----------------------------------------------------------------------
C   Table of Integrals:
C  ---------------------
C
C    C(I) =  integral of  Z^(I-1) * EXP(Z^2)  over all Z.
C
      DATA C /            1.77245390415192D0,                           !  0
     +                    0.0D0,                                        !  1
     +                    0.886226952075958D0,                          !  2
     +                    0.0D0,                                        !  3
     +                    1.32934042811394D0,                           !  4
     +                    0.0D0,                                        !  5
     +                    3.32335107028484D0,                           !  6
     +                    0.0D0,                                        !  7
     +                   11.6317287459969D0,                            !  8
     +                    0.0D0,                                        !  9
     +                   52.3427793569863D0/                            ! 10
C-----------------------------------------------------------------------
C
      DATA ETA1(1),ETA2(1) / 1.0D0, 1.0D0/
C
C
      NX1P1 = NX1+1
      NX2P1 = NX2+1
C
      BETA = (ZETA1*X1+ZETA2*X2)*ALFA
      GAMA = (ZETA1*X1*X1+ZETA2*X2*X2)
      DEL1 = BETA-X1/ALFA
      DEL2 = BETA-X2/ALFA
      EPSN = EXP(BETA*BETA-GAMA)
C
      DO I = 2,NX1P1
         ETA1(I) = ETA1(I-1)*DEL1
      END DO
      DO I = 2,NX2P1
         ETA2(I) = ETA2(I-1)*DEL2
      END DO
C
      THETA = 1.0D0
      DO I = 1,NX1+NX2P1
         THETA = THETA*ALFA
      END DO
C
      OLAPSUM = 0.0D0
      DO I = 1,NX1P1
         N1RED = NX1P1-I
         DO J = 1,NX2P1
            N2RED = NX2P1-J
            OLAPSUM = OLAPSUM+BINOM(NX1P1,I)*BINOM(NX2P1,J)
     +                *ETA1(I)*ETA2(J)*C(N1RED+N2RED+1)
         END DO
      END DO
C
      OLAP1D = THETA*EPSN*OLAPSUM
C
      RETURN
      END

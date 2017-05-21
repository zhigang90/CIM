c ======================================================================
c
      SUBROUTINE AbunMASS(NATOMS,IAN,AtMASS)
      IMPLICIT REAL*8(A-H,O-Z)
C
C  Assign atomic masses using the most adundant isotope
C
      DIMENSION IAN(NATOMS),AtMASS(NATOMS)
      parameter (nelem=92)
      dimension amass(nelem)
c
c --- define atomic weights for the first <nelem> elements ---
c --- data available on-line from NIST Standard Reference Database 144
c --- Atomic Weights and Isotopic Compositions   July 2010
c
c                         H - Zr
      data (amass(i), i = 1 , 40) /
     1  1.007825d0,  4.002603d0,  7.016005d0,  9.012182d0, 11.009305d0,
     2 12.000000d0, 14.003074d0, 15.994915d0, 18.998403d0, 19.992440d0,
     3 22.989769d0, 23.985042d0, 26.981539d0, 27.976926d0, 30.973762d0,
     4 31.972071d0, 34.968853d0, 39.962383d0, 38.963707d0, 39.962591d0,
     5 44.955912d0, 47.947946d0, 50.943960d0, 51.940508d0, 54.938045d0,
     6 55.934937d0, 58.933195d0, 57.935343d0, 62.929598d0, 63.929142d0,
     7 68.925574d0, 73.921178d0, 74.921597d0, 79.916521d0, 78.918337d0,
     8 83.911507d0, 84.911790d0, 87.905612d0, 88.905848d0, 89.904704d0 /
c                         Nb - Hg
      data (amass(i), i = 41 , 80) /
     1  92.906378d0, 97.905408d0, 97.907216d0,101.904349d0,102.905504d0,
     2 105.903486d0,106.905097d0,113.903359d0,114.903878d0,119.902195d0,
     3 120.903816d0,129.906224d0,126.904473d0,131.904154d0,132.905452d0,
     4 137.905247d0,138.906353d0,139.905439d0,140.907653d0,141.907723d0,
     5 144.912749d0,151.919732d0,152.921230d0,157.924104d0,158.925347d0,
     6 163.929175d0,164.930322d0,165.930293d0,168.934213d0,173.938862d0,
     7 174.940772d0,179.946550d0,180.947996d0,183.950931d0,186.955753d0,
     8 191.961481d0,192.962926d0,194.964791d0,196.966569d0,201.970643d0/
c                         Tl - U
      data (amass(i), i = 81 , 92) /
     1 204.974428d0,207.976652d0,208.980399d0,208.982430d0,208.987148d0,
     2 222.017578d0,223.019736d0,226.025410d0,227.027752d0,232.038055d0,
     3 231.035884d0,238.050788d0 /
c
c
      DO 10 IATM=1,NATOMS
      IATNO = IAN(IATM)
      IF(IATNO.LT.1.OR.IATNO.GT.nelem) THEN
        Call nerror(2,'GEOMETRY/FREQUENCY module',
     $    'Unknown Atom  Cannot assign atomic mass',IATM,0)
      ENDIF
      AtMASS(IATM) = amass(IATNO)
 10   CONTINUE
c
      RETURN
      END
c ======================================================================

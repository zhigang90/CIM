      subroutine opt2(inp,optcyc,convgd)
c  Main routine for the alternative optimizer
      implicit real*8(a-h,o-z)
      integer optcyc
      logical convgd
c  Arguments
c  INTENT(IN)
c  inp       = input file unit number
c  optcyc    = optimization cycle number
c  INTENT(OUT)
c  convgd    = true only if convergence has been achieved
c
c  This routine sets the default parameters, and reads the command line
c  The parameters are taken from the command line
c  The geometries and forces are read from the standard files
      end
c=======================================================================

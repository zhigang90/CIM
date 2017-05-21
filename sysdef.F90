module sysdef

! system dependent values

  implicit none

              ! directory separator character

#ifdef WINDOWS
  character(len=1), parameter ::  DIR_SEP = '\'
#else
  character(len=1), parameter ::  DIR_SEP = '/'
#endif

end module sysdef

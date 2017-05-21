!     ----------------------------------
!     pqsmsg.h
!
!     Messages for PQS program
!     ----------------------------------

module newpara

  use kinds
  implicit none
  private

! It is not a good practice to reuse tags in different contexts.
! This can cause very subtle bugs in the parallel parts as the code
! avoids the use of synchronizing barriers for the sake of speed.

! --------------------------
  integer(kind=i_4), public, parameter :: TxFailure     =  666
  integer(kind=i_4), public, parameter :: TxSlaveError  =  777
  integer(kind=i_4), public, parameter :: TxKissofDeath =  888
  integer(kind=i_4), public, parameter :: TxMemSize     =    1
  integer(kind=i_4), public, parameter :: TxJobType     =    2
  integer(kind=i_4), public, parameter :: TxNslave      =    3

  integer(kind=i_4), public, parameter :: TxNext        =    5
  integer(kind=i_4), public, parameter :: TxTime        =    6

  integer(kind=i_4), public, parameter :: TxStore       =    8
  integer(kind=i_4), public, parameter :: TxDone        =    9
  integer(kind=i_4), public, parameter :: TxContinue    =   10
  integer(kind=i_4), public, parameter :: TxId          =   11
  integer(kind=i_4), public, parameter :: TxMask        =   12
!
! TxJobType values are not proper tags, but integers to be sent as a msg,
! thus they are not declared parameters, as some MPI calls require
! writable variables
! NZG_6/28/2016 Add tags for CIM calculation

  integer, public :: TxDoScf       =   40
  integer, public :: TxDoPost      =   80
  integer, public :: TxDoProp      =  120
  integer, public :: TxDoMP2       =  140
  integer, public :: TxDoHess      =  160
  integer, public :: TxDoCosmoForc =  200
  integer, public :: TxCIMGen      =  350
  integer, public :: TxCIMSubMP2   =  360
  integer, public :: TxCIMSubCC    =  370
  integer, public :: TxFinish      =  999

  integer(kind=i_4), public, parameter :: TxJobInit     =   22
  integer(kind=i_4), public, parameter :: TxScfInit     =   42
  integer(kind=i_4), public, parameter :: TxDftInit     =   62
  integer(kind=i_4), public, parameter :: TxPostInit    =   82
  integer(kind=i_4), public, parameter :: TxPropInit    =  122
  integer(kind=i_4), public, parameter :: TxMP2Init     =  142
  integer(kind=i_4), public, parameter :: TxHessInit    =  162

  integer(kind=i_4), public, parameter :: TxInitVar     =   30
  integer(kind=i_4), public, parameter :: TxInitArr     =   31
  integer(kind=i_4), public, parameter :: TxSymData     =   32

  integer(kind=i_4), public, parameter :: TxScfDens     =   44
  integer(kind=i_4), public, parameter :: TxDftDat      =   64
  integer(kind=i_4), public, parameter :: TxDftNDat     =   74
  integer(kind=i_4), public, parameter :: TxPostDens    =   84
  integer(kind=i_4), public, parameter :: TxPropDat     =  124
  integer(kind=i_4), public, parameter :: TxMP2Dat      =  144
  integer(kind=i_4), public, parameter :: TxHessDat     =  164

  integer(kind=i_4), public, parameter :: TxBlockReq    =   46
  integer(kind=i_4), public, parameter :: TxDftReq      =   66
  integer(kind=i_4), public, parameter :: TxDftNReq     =   76
  integer(kind=i_4), public, parameter :: TxStorReq     =   56
  integer(kind=i_4), public, parameter :: TxShReq       =   96

  integer(kind=i_4), public, parameter :: TxBlockAssign =   48
  integer(kind=i_4), public, parameter :: TxDftAssign   =   68
  integer(kind=i_4), public, parameter :: TxDftNAssign  =   78
  integer(kind=i_4), public, parameter :: TxStorAssign  =   58
  integer(kind=i_4), public, parameter :: TxShAssign    =   98
  integer(kind=i_4), public, parameter :: TxMP2Assign   =  108

  integer(kind=i_4), public, parameter :: TxScfFA       =   50
  integer(kind=i_4), public, parameter :: TxScfFB       =   51
  integer(kind=i_4), public, parameter :: TxDftF1       =   61
  integer(kind=i_4), public, parameter :: TxDftF2       =   63
  integer(kind=i_4), public, parameter :: TxDftF3       =   65
  integer(kind=i_4), public, parameter :: TxDftF4       =   67

  integer(kind=i_4), public, parameter :: TxDftF1b      =  131
  integer(kind=i_4), public, parameter :: TxDftF2b      =  133
  integer(kind=i_4), public, parameter :: TxDftF3b      =  135
  integer(kind=i_4), public, parameter :: TxDftF4b      =  137
  integer(kind=i_4), public, parameter :: TxDftHmp      =  138

  integer(kind=i_4), public, parameter :: TxDftN        =   71
  integer(kind=i_4), public, parameter :: TxDftN2       =   72
  integer(kind=i_4), public, parameter :: TxPostFock    =   90

  integer(kind=i_4), public, parameter :: TxSpinDat     =  125
  integer(kind=i_4), public, parameter :: TxEFGDat      =  126
  integer(kind=i_4), public, parameter :: TxProp1       =  127
  integer(kind=i_4), public, parameter :: TxProp2       =  128
  integer(kind=i_4), public, parameter :: TxProp3       =  129

  integer(kind=i_4), public, parameter :: TxMP2File     =  145
  integer(kind=i_4), public, parameter :: TxMP2S1       =  146
  integer(kind=i_4), public, parameter :: TxMP2S2       =  147
  integer(kind=i_4), public, parameter :: TxMP2S3       =  148
  integer(kind=i_4), public, parameter :: TxMP2S4       =  149
  integer(kind=i_4), public, parameter :: TxMP2S5       =  150
  integer(kind=i_4), public, parameter :: TxMP2S6       =  151

  integer(kind=i_4), public, parameter :: TxBinInit     =  152
  integer(kind=i_4), public, parameter :: TxBinDat      =  153
  integer(kind=i_4), public, parameter :: TxBinDat1     =  154
  integer(kind=i_4), public, parameter :: TxBinDat2     =  155
  integer(kind=i_4), public, parameter :: TxBinS1       =  156
  integer(kind=i_4), public, parameter :: TxBinS2       =  157

  integer(kind=i_4), public, parameter :: TxHess1       =  163
  integer(kind=i_4), public, parameter :: TxHess2       =  165
  integer(kind=i_4), public, parameter :: TxHess3       =  167
  integer(kind=i_4), public, parameter :: TxHess4       =  169

  integer(kind=i_4), public, parameter :: TxWDens1      =  171
  integer(kind=i_4), public, parameter :: TxWDens2      =  172
  integer(kind=i_4), public, parameter :: TxWDens3      =  173
  integer(kind=i_4), public, parameter :: TxWDens4      =  174
  integer(kind=i_4), public, parameter :: TxWDens5      =  175

  integer(kind=i_4), public, parameter :: TxCosmoDat    =  201
  integer(kind=i_4), public, parameter :: TxCosmoReq    =  202
  integer(kind=i_4), public, parameter :: TxCosmoAss    =  203
  integer(kind=i_4), public, parameter :: TxCosmoDen    =  204
  integer(kind=i_4), public, parameter :: TxCosmoPot    =  205
  integer(kind=i_4), public, parameter :: TxCosmoCha    =  206
  integer(kind=i_4), public, parameter :: TxCosmoH0     =  207
  integer(kind=i_4), public, parameter :: TxCosmoForc   =  208
  integer(kind=i_4), public, parameter :: TxCosmoTime   =  209 

  integer(kind=i_4), public, parameter :: TxFTCInit     =  182
  integer(kind=i_4), public, parameter :: TxFTCInit0    =  183
  integer(kind=i_4), public, parameter :: TxFTCInit1    =  184
  integer(kind=i_4), public, parameter :: TxFTCAssign   =  185
  integer(kind=i_4), public, parameter :: TxFTCDen      =  186
  integer(kind=i_4), public, parameter :: TxFTCBTM      =  187
  integer(kind=i_4), public, parameter :: TxFTCFA       =  188
  integer(kind=i_4), public, parameter :: TxFTCFB       =  189
  integer(kind=i_4), public, parameter :: TxFTCFFT      =  190
  integer(kind=i_4), public, parameter :: TxFTCMult     =  191
  integer(kind=i_4), public, parameter :: TxFTCGX       =  192
  integer(kind=i_4), public, parameter :: TxFTCGY       =  193
  integer(kind=i_4), public, parameter :: TxFTCGZ       =  194
!
! CC
!
  integer          , public            :: TxDoCCSD      =  230
  integer(kind=i_4), public, parameter :: TxCCSDInit    =  231
  integer(kind=i_4), public, parameter :: TxCCSDJob     =  232
  integer(kind=i_4), public, parameter :: TxCCSDReq     =  233
  integer(kind=i_4), public, parameter :: TxCCSDRes     =  234
  integer(kind=i_4), public, parameter :: TxCCSDIter    =  235
  integer(kind=i_4), public, parameter :: TxCCSDSing    =  236
  integer(kind=i_4), public, parameter :: TxCoefInit    =  237
  integer(kind=i_4), public, parameter :: TxCCSDMP2Req  =  238
  integer(kind=i_4), public, parameter :: TxCCSDMP2Job  =  239
  integer(kind=i_4), public, parameter :: TxCCSDMP2Ite  =  240
  integer(kind=i_4), public, parameter :: TxCCSDJob1    =  241
  integer(kind=i_4), public, parameter :: TxCCSDReq1    =  242
  integer(kind=i_4), public, parameter :: TxCCSDRes1    =  243
  integer(kind=i_4), public, parameter :: TxCCSDJobQ    =  244
  integer(kind=i_4), public, parameter :: TxCCSDReqQ    =  245
  integer(kind=i_4), public, parameter :: TxCCSDResQ    =  246
  integer(kind=i_4), public, parameter :: TxCCSDAReq    =  247
  integer(kind=i_4), public, parameter :: TxCCSDAJob    =  248
  integer(kind=i_4), public, parameter :: TxCCSDARes    =  249
  integer(kind=i_4), public, parameter :: TxCCSDQuit    =  250
  integer(kind=i_4), public, parameter :: TxCCSDRedu    =  251
  integer(kind=i_4), public, parameter :: TxCCSDRedu1   =  252
  integer(kind=i_4), public, parameter :: TxCCSDRedu2   =  253
  integer(kind=i_4), public, parameter :: TxCCSDRedu3   =  254
  integer(kind=i_4), public, parameter :: MP2file       =  255
  integer(kind=i_4), public, parameter :: MP2_G_GEN_req =  256
  integer(kind=i_4), public, parameter :: MP2_G_GEN_job =  257
  integer(kind=i_4), public, parameter :: TxDiisReq     =  258
  integer(kind=i_4), public, parameter :: TxDiisJob     =  259
  integer(kind=i_4), public, parameter :: TxDiisRes     =  260
  integer(kind=i_4), public, parameter :: TxDIISFlow    =  261
  integer(kind=i_4), public, parameter :: TxTriplesInit =  262
  integer(kind=i_4), public, parameter :: TxTriplesReq  =  263
  integer(kind=i_4), public, parameter :: TxTriplesJob  =  264
  integer(kind=i_4), public, parameter :: TxTriplesRes  =  265
  integer(kind=i_4), public, parameter :: TxCCSDJob2    =  246
  integer(kind=i_4), public, parameter :: TxCCSDReq2    =  247
  integer(kind=i_4), public, parameter :: TxCCSDRes2    =  248
  integer(kind=i_4), public, parameter :: TxCCSDJob3    =  249
  integer(kind=i_4), public, parameter :: TxCCSDReq3    =  250
  integer(kind=i_4), public, parameter :: TxCCSDRes3    =  251
  integer(kind=i_4), public, parameter :: TxCCSDJob4    =  252
  integer(kind=i_4), public, parameter :: TxCCSDReq4    =  253
  integer(kind=i_4), public, parameter :: TxCCSDRes4    =  254
  integer(kind=i_4), public, parameter :: TxCCSDIter1   =  277
  integer(kind=i_4), public, parameter :: CCSD_Rec_Req  =  278
  integer(kind=i_4), public, parameter :: CCSD_Rec_Off  =  279
  integer          , public            :: TxDoROHF      =  300
  integer(kind=i_4), public, parameter :: ROHFiter      =  301
  integer(kind=i_4), public, parameter :: ROHFdata      =  302
  integer(kind=i_4), public, parameter :: ROHFwork      =  303
  integer(kind=i_4), public, parameter :: ROHFred       =  304

! CIM
  integer(kind=i_4), public, parameter :: TxCIMVir      =  351
  integer(kind=i_4), public, parameter :: TxCIMDat      =  352
  integer(kind=i_4), public, parameter :: TxCIMMP2Int   =  353
  integer(kind=i_4), public, parameter :: TxCIMInfo     =  354
  integer(kind=i_4), public, parameter :: TxCIMCCSDInt  =  371

  
     ! number of slaves

  integer, public :: NSLV = 0

     ! array for slave ids, slave hostnames, master hostname

  integer(kind=i_4), public, allocatable :: SLVID(:)
  character(len=256), public, allocatable :: HOSTSL(:)
  character(len=256), public :: HOSTM

  character(len=256), public :: MY_HOSTNM

  integer, public :: MASTER_ID, MY_ID, MY_GID

    ! no. of hosts, hostnames, no. of processes on each host,
    ! no. of cpus on each host, id of processes per host

  integer, public :: NHOST, MAXPROCH, MAXCPUH, MY_HOSTID
  character(len=256), public, allocatable :: HOSTNAMES(:)
  integer, public, allocatable :: NPROCH(:),NCPUH(:),IDH(:,:)

    ! initial and current cpu affinity mask

  integer*4, public, allocatable :: ICPUMASK(:), ICPU(:)
  integer*4, public, allocatable :: MASKSL(:,:)

    ! process ID

  integer*4, public :: MY_PID
  integer*4, public, allocatable :: PIDSL(:)

     ! mpi buffer location, length, position

  integer*4, public ::  MPI_IBUF, MPI_IPOS
  integer, public :: MPI_IBFP

     ! signal number used by the signal handler

  integer(kind=i_4), public, parameter :: SIGINT = 2  ! CTRL+C

     ! type of default integer. This variable will contain either
     ! the PVM of MPI type for the default integer. 
     ! 
     ! It is initialized to 9999, so that it should produce a 
     ! message error if it is used without being set by calling
     ! either set_pvm_def_int or set_mpi_def_int

  integer(kind=i_4), public :: INTEGER_DEF = 9999  
  integer(kind=i_4), public :: NEW_WORLD

end module newpara

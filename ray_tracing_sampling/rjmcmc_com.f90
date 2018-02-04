!==============================================================================
MODULE RJMCMC_COM
   USE MPI
   USE DATA_TYPE
   IMPLICIT NONE

!!
!! General switches
!!
!! Below, SWD is replaced with RT (ray tracing)
!!
   INTEGER(KIND=IB) :: IMAP       !! WRITE REPLICA AND EXIT
   INTEGER(KIND=IB) :: ICOV       !! 0 = Sample implicit over sigma
                                  !! 1 = Sample over sigma

   INTEGER(KIND=IB) :: I_RT       !! Use the RT data ?       
   INTEGER(KIND=IB) :: I_VARPAR       !! Type of Parametrization       

   INTEGER(KIND=IB) :: ENOS       !! 1 = Turn on even numbered order stats
   INTEGER(KIND=IB) :: IPOIPR     !! 1 = Turn on Poisson prior on k
   INTEGER(KIND=IB) :: IAR        !! 1 = Use Autoregressive error model
   INTEGER(KIND=IB) :: IBD_SINGLE !! 1 = include BD for single parameters onto nodes
   INTEGER(KIND=IB) :: iraysum    !! 1 = invert use raysum 0 = use ray3d
   INTEGER(KIND=IB) :: I_VREF     !! 1 = sample VpVs ratio
   INTEGER(KIND=IB) :: I_VPVS     !! 1 = sample VpVs ratio
   INTEGER(KIND=IB) :: ISMPPRIOR  !! Sample from the prior (sets logL = 1._RP)
   INTEGER(KIND=IB) :: ISETSEED   !! Fix the random seed 
   INTEGER(KIND=IB) :: IEXCHANGE  !! 1 = turn on exchange moves (parallel tempering)

!!
!! Model and data dimensions
!!
   INTEGER(KIND=IB)            :: NLMN         ! Min number of layers
   INTEGER(KIND=IB)            :: NLMX         ! Max number of layers
   INTEGER(KIND=IB)            :: NMODE        ! Number of Phases in my case ( or modes for SWD)

   
   INTEGER(KIND=IB)            :: NPL          ! No. parameters per layer
   INTEGER(KIND=IB)            :: NDAT_RT        ! No. of data points

   ! Forward model specifics for Ray Tracer
   INTEGER(KIND=IB)            :: NRAYS        ! Number of rays to compute in Forward Model
   INTEGER(KIND=IB)            :: NSRC        ! Number of rays to compute in Forward Model
   REAL(KIND=RP), ALLOCATABLE, DIMENSION(:):: src_offset      ! The offset (source-receiver distance projected on the surface)
   REAL(KIND=RP), ALLOCATABLE, DIMENSION(:):: src_depth      ! The offset (source-receiver distance projected on the surface)

   CHARACTER(len=64) :: filebasefile      = 'filebase.txt'

!! 
!! Forward RF specific 
!!
  !!                                                      thick      rho      alph     beta    %P    %S    tr    pl    st    di
  REAL(KIND=RP)             :: hmx                !! Max crustal depth in km
  REAL(KIND=RP),DIMENSION(2):: sdmn               !! Min standard deviation
  REAL(KIND=RP),DIMENSION(2):: sdmx               !! Max standard deviation

  !REAL(KIND=SP),ALLOCATABLE,DIMENSION(:)       :: baz2,slow2,sta_dx,sta_dy
  !INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:)    :: nseg
  ! INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:,:):: phaselist
  !INTEGER(KIND=IB),PARAMETER:: mults   = 2  !! Parameter for multiples in layers
  !INTEGER(KIND=IB),PARAMETER:: out_rot = 1  !! 
  !INTEGER(KIND=IB),PARAMETER:: align   = 1  !! 
  !REAL(KIND=SP)             :: shift2
  !REAL(KIND=SP)             :: sampling_dt
  !REAL(KIND=SP),PARAMETER   :: sig   = 0.01_SP
  !REAL(KIND=SP),PARAMETER   :: VPVS  = 1.75_RP
  !REAL(KIND=SP)             :: width2             !! Set < 0 to return impulse response
  !REAL(KIND=SP)             :: wl                 !! water-level for ray3d
  
!!
!!  Prior variables and good seeding model
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlim
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxlim
   INTEGER(KIND=IB)            :: kmin               ! Min number of layers
   INTEGER(KIND=IB)            :: kmax               ! Max number of layers
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pk       ! Poisson prior on k
   REAL(KIND=RP)               :: lambda             ! Lambda parameter for Poisson prior on k
   REAL(KIND=RP)               :: hmin               ! Min allowed layer thickness
   REAL(KIND=RP),PARAMETER     :: fact     = 1.00_RP ! factor for rotated space perturbation
   REAL(KIND=RP),PARAMETER     :: factdelay= 1.50_RP ! shrinking factor for delayed rejection (>1.)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: maxpert
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsc
   REAL(KIND=RP)               :: vs01,vs02,vsh1,vsh2,numin,numax
   REAL(KIND=RP)               :: area_bn,area_kskp,area_cr
!!
!!  Autoregressive model prior variables:
!! Below, SWD is replaced with RT (ray tracing)
!!
   REAL(KIND=RP)              :: armxH      = 0.2_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP)              :: armxV      = 0.2_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP)              :: armxRT    = 0.5_RP           ! Max AR and ARI model range (amplitude units)
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimar, maxlimar, maxpertar
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsd, pertarsdsc
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimarRT, maxlimarRT, maxpertarRT
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertarsdRT, pertarsdscRT

!!
!!  Standard deviation prior variables:
!!
!! Below, SWD is replaced with RT (ray tracing)
!!
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsd, maxlimsd, maxpertsd
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsd, pertsdsdsc
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: minlimsdRT, maxlimsdRT, maxpertsdRT
   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: pertsdsdRT, pertsdsdscRT

!   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE:: fr,fstep   ! Total frequency array
!   REAL(KIND=RP)               :: z_t
!   REAL(KIND=RP)               :: cw
!   REAL(KIND=RP)               :: rw
   CHARACTER(len=64) :: filebase
   INTEGER(KIND=IB)  :: filebaselen
   !CHARACTER(LEN=100) :: infileV
   !CHARACTER(LEN=100) :: infileH
   !CHARACTER(LEN=100) :: infileT
   CHARACTER(LEN=100) :: infileRT
   CHARACTER(LEN=100) :: infileref
   CHARACTER(LEN=100) :: source_data_file

   !CHARACTER(LEN=100) :: repfileR
   !CHARACTER(LEN=100) :: repfileT
   !CHARACTER(LEN=100) :: repfileV
   !CHARACTER(LEN=100) :: repfileS
   !CHARACTER(LEN=100) :: repfileRF
   CHARACTER(LEN=100) :: repfileRT
   CHARACTER(LEN=100) :: repfileRTar
   CHARACTER(LEN=100) :: parfile
   CHARACTER(LEN=64) :: logfile
   CHARACTER(LEN=64) :: seedfile
   CHARACTER(len=64) :: mapfile
   CHARACTER(len=64) :: covfile
  ! CHARACTER(LEN=64)  :: obsfile
   CHARACTER(LEN=64)  :: arfile
   !CHARACTER(LEN=64)  :: predfile
   CHARACTER(LEN=64)  :: obsfileRT
   CHARACTER(LEN=64)  :: arfileRT
   CHARACTER(LEN=64)  :: predfileRT
   CHARACTER(len=64) :: sdfile
   CHARACTER(len=64) :: samplefile
   CHARACTER(len=64) :: stepsizefile
   CHARACTER(len=64) :: lincovfile
   CHARACTER(LEN=64)  :: modname = 'sample.geom'

!!
!! Velocity reference model
!!
  INTEGER(KIND=IB)                        :: NVELREF, NPREM
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: vel_ref, vel_prem

!!
!! Parallel Tempering parameters
!!
  INTEGER(KIND=IB)                         :: NPTCHAINS1            !! # chains T=1
  REAL(KIND=RP)                            :: dTlog                 ! Temperature increment
  INTEGER(KIND=IB)                         :: NT                    !! # tempering levels (temperatures)
  INTEGER(KIND=IB)                         :: NPTCHAINS             !! # parallel tempering chains
  INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: NCHAINT
  INTEGER(KIND=IB)                         :: ncswap     = 0_IB     !! Temp swap accepted counter
  INTEGER(KIND=IB)                         :: ncswapprop = 0_IB     !! Temp swap proposed counter
  REAL(KIND=RP),ALLOCATABLE,DIMENSION(:)   :: beta_pt               !! Temperature array parallel tempering
  INTEGER(KIND=IB)                         :: ibirth  = 0,ideath  = 0
  INTEGER(KIND=IB)                         :: ibirths = 0,ideaths = 0

!!
!!  Sampling specific parameters
!!
   INTEGER(KIND=IB)           :: NFPMX
   INTEGER(KIND=IB)           :: NFPMX2
   INTEGER(KIND=IB)           :: ioutside   = 0
   INTEGER(KIND=IB)           :: ireject    = 0, iaccept = 0, iaccept_delay = 0, ireject_delay = 0
   INTEGER(KIND=IB)           :: i_bd,i_bds     ! Birth-Death track (0=MCMC, 1=birth, 2=death)
   INTEGER(KIND=IB)           :: i_sdpert = 0   ! if sigma is perturbed, don't compute forward model
   INTEGER(KIND=IB)           :: ishearfail = 0 ! if sigma is perturbed, don't compute forward model
   INTEGER(KIND=IB)           :: i_ref_nlay = 0 ! if sigma is perturbed, don't compute forward model

!!
!!  Convergence parameters
!!
   INTEGER(KIND=IB)       :: iconv    = 0       ! Convergence switch slaves
   INTEGER(KIND=IB)       :: iconv2   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iconv3   = 0       ! Convergence switch master
   INTEGER(KIND=IB)       :: iarfail  = 0       ! Tracks failure of AR model when predicted AR series too large

!!
!! RJMCMC parameters
!!
   INTEGER(KIND=IB),PARAMETER    :: NCHAIN     = 1E5_IB  ! # iterations (max # MCMC steps)
   INTEGER(KIND=IB)              :: ICHAINTHIN = 5E0_IB  ! Chain thinning interval
   !INTEGER(KIND=IB)              :: NKEEP      = 1E1_IB  ! Number models to keep before writing
   INTEGER(KIND=IB)              :: NKEEP      = 1E1_IB  ! Number models to keep before writing
   
   INTEGER(KIND=IB),PARAMETER    :: NAP        = 10      ! Misc parameters in sample (for bookeeping)
   INTEGER(KIND=IB),PARAMETER    :: NDM        = 100     ! No. steps in lin rot est
   INTEGER(KIND=IB)          :: TCHCKPT              !! No. seconds (integer value) between checkpoints
   INTEGER(KIND=IB)          :: icheckpoint          !! No. of checkpoints to data (read from checkpoint/status.txt)

   REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:)   :: sdevm  ! Std dev for perturbations
!!
!!  Structures for objects and data  a
!! Below, SWD is replaced with RT (ray tracing)
!!
  INTEGER :: imcmc1 = 1   !! Counter for models at T=1 (needs to survive checkpointing!)
  INTEGER :: imcmc2 = 1   !! Counter for mcmc steps to scale diminishing adaptation (needs to survive checkpointing!)
  INTEGER :: NFIELD = 27  !! The number of fields in objstruc
  INTEGER :: objtype1     !! Name of objtype for MPI sending
  INTEGER :: objtype2     !! Name of objtype for MPI sending
  INTEGER :: objtype3     !! Name of objtype for MPI sending
   TYPE :: objstruc
      SEQUENCE 
      INTEGER(KIND=IB)                        :: k          ! No. nodes
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: voro       ! 1D layer nodes
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: voroidx    ! 1D layer nodes index
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     ! Forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: hiface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: ziface
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparRT      !! Std dev RT data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveRT      !! Std dev RT data
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arpar         !! AR model forward parameters
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparRT      !! AR model forward parameters
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxar      !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarRT      !! AR on/off index (1=on)
      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: gvoroidx   !! Index of live parameters on birth/death node
      INTEGER(KIND=IB)                        :: nunique     !! 
      INTEGER(KIND=IB)                        :: NFP         !! Number forward parameters
      REAL(KIND=RP)                           :: beta
      REAL(KIND=RP)                           :: logL        !! log likelihood
      REAL(KIND=RP)                           :: logPr       !! log Prior probability ratio
      REAL(KIND=RP)                           :: tcmp
      INTEGER(KIND=IB)                        :: ireject_bd = 0
      INTEGER(KIND=IB)                        :: iaccept_bd = 0
      INTEGER(KIND=IB)                        :: ireject_bds = 0
      INTEGER(KIND=IB)                        :: iaccept_bds = 0
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DobsRT     !! Observed data RT
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DpredRT    !! Predicted RT data for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DresRT     !! RT data residuals for trial model
      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: DarRT      !! RT autoregressive model predicted data
   END TYPE objstruc
!!
!! Structure for covariance matrices (only applies for ICOV >= 2)
!!
!  TYPE :: covstruc
!    REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE :: Cdi   ! Inverse covariance matrix
!  END TYPE covstruc
!!
!! File units for IO
!!
  INTEGER(KIND=IB) :: usample, ustep, ulog
  INTEGER(KIND=IB) :: reclen

   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE      :: icount

!!
!!  Buffer for sdave for AR discrimination in CHECKBOUNDS_AR
!!
   INTEGER(KIND=IB), PARAMETER               :: NBUF = 100
   REAL(KIND=RP),DIMENSION(:,:,:),ALLOCATABLE:: sdbuf
!!
!!  Global variables
!!
   REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE   :: sample                    ! Posterior sample
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: tmpmap                    ! temporary for reading map
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buf_save_snd,buf_save_rcv ! Buffers for MPI sending
   INTEGER(KIND=IB),DIMENSION(:),ALLOCATABLE  :: buffer1                   !
   REAL(KIND=RP),DIMENSION(:),ALLOCATABLE     :: buffer2,buffer3           !

!!
!!  MPI global variables
!!
   INTEGER(KIND=IB)            :: rank,NTHREAD,ierr
   INTEGER(KIND=IB)            :: ncount1,ncount2
   INTEGER(KIND=IB), PARAMETER :: src = 0_IB
   INTEGER                     :: to,from,COMM
   INTEGER                     :: status(MPI_STATUS_SIZE)
   INTEGER(KIND=IB)            :: isize1,isize2,isize3
   INTERFACE
      FUNCTION RANDPERM(num)
         USE data_type, ONLY : IB
         IMPLICIT NONE
         INTEGER(KIND=IB), INTENT(IN) :: num
         INTEGER(KIND=IB), DIMENSION(num) :: RANDPERM
      END FUNCTION RANDPERM
   END INTERFACE
   REAL(KIND=RP) :: tsave1, tsave2              ! Overall time 

!!
!!  FFTW stuff
!!
   INTEGER*8 :: planR,planC

  CONTAINS
  !==============================================================================
  integer function newunit(unit)
  !==============================================================================
    integer, intent(out), optional :: unit
    integer, parameter :: LUN_MIN=10, LUN_MAX=1000
    logical :: opened
    integer :: lun
    newunit=-1
    do lun=LUN_MIN,LUN_MAX
      inquire(unit=lun,opened=opened)
      if (.not. opened) then
        newunit=lun
        exit
      end if
    end do
    if (present(unit)) unit=newunit
  end function newunit

END MODULE RJMCMC_COM
!=======================================================================

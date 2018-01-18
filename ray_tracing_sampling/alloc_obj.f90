!=======================================================================

SUBROUTINE MAKE_MPI_STRUC_SP(obj,objtype)
!!=======================================================================
!!
!! 
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj
INTEGER(KIND=IB) :: ifield,objtype
INTEGER(KIND=IB) :: oldtypes(NFIELD), blockcounts(NFIELD)
INTEGER(KIND=MPI_ADDRESS_KIND) :: offsets(NFIELD)
INTEGER(KIND=IB) :: iextent,rextent,dextent

!  Need to first figure offset by getting size of MPI_REAL etc 
call MPI_TYPE_EXTENT(MPI_INTEGER, iextent, ierr)
call MPI_TYPE_EXTENT(MPI_REAL, rextent, ierr)
call MPI_TYPE_EXTENT(MPI_DOUBLE_PRECISION, dextent, ierr)

!!
!! Lengths of arrays
!!
!!                k   voro     voroidx           par          hiface     ziface      
blockcounts = (/ 1,  NLMX*NPL, NLMX*NPL , (NLMX+1)*NPL*NPL,  NPL*NLMX,   NPL*NLMX,&
!!               sdparR    sdparV   sdparT   sdparSWD   sdaveH    sdaveV   sdaveT   sdaveSWD
                 NRF1,      NRF1,    NRF1,    NMODE,     NRF1,     NRF1,    NRF1,    NMODE,&
!!               arpar   arparSWD idxar   idxarSWD  gvoroidx nunique NFP
                 3*NRF1,  NMODE,  3*NRF1,  NMODE,     NPL-1,     1,    1, & 
!!               beta    logL     logPr tcmp  ireject_bd  iaccept_bd ireject_bds  iaccept_bds 
                  1,       1,       1,    1,      1,          1,         1,          1,&
!!               DobsR           DpredR        DobsV       DpredV         DobsT        DpredT
                 NRF1*NTIME2,   NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME, &
!!                S            DresR         DresV        DresT        DarR        DarV        DarT
                 NRF1*NSRC, NRF1*NTIME,   NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME, &
!!                DobsSWD           DpredSWD        DresSWD       DarSWD           periods
                 NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD, NMODE*NDAT_SWD /)

!  INTEGER :: objtype3     !! Name of objtype for MPI sending
!   TYPE :: objstruc
!      INTEGER(KIND=IB)                        :: k          ! No. nodes
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: voro       ! 1D layer nodes
!      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:,:):: voroidx    ! 1D layer nodes index
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: par     ! Forward parameters
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: hiface
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: ziface
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdparSWD      !! Std dev SWD data
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: sdaveSWD      !! Std dev SWD data
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arpar         !! AR model forward parameters
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:):: arparSWD      !! AR model forward parameters
!      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxar      !! AR on/off index (1=on)
!      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: idxarSWD      !! AR on/off index (1=on)
!      INTEGER(KIND=IB),ALLOCATABLE,DIMENSION(:):: gvoroidx   !! Index of live parameters on birth/death node
!      INTEGER(KIND=IB)                        :: nunique     !! 
!      INTEGER(KIND=IB)                        :: NFP         !! Number forward parameters
!      REAL(KIND=RP)                           :: beta
!      REAL(KIND=RP)                           :: logL        !! log likelihood
!      REAL(KIND=RP)                           :: logPr       !! log Prior probability ratio
!      REAL(KIND=RP)                           :: tcmp
!      INTEGER(KIND=IB)                        :: ireject_bd = 0
!      INTEGER(KIND=IB)                        :: iaccept_bd = 0
!      INTEGER(KIND=IB)                        :: ireject_bds = 0
!      INTEGER(KIND=IB)                        :: iaccept_bds = 0
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: Dobs_tt     !! Observed data SWD
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: Dpred_tt    !! Predicted SWD data for trial model
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: Dres_tt     !! SWD data residuals for trial model
!      REAL(KIND=RP),ALLOCATABLE,DIMENSION(:,:):: Dar_tt      !! SWD autoregressive model predicted data
!   END TYPE objstruc

!! Oldtypes
oldtypes(1)     = MPI_INTEGER
oldtypes(2)     = MPI_DOUBLE_PRECISION
oldtypes(3)     = MPI_INTEGER
oldtypes(4:16)   = MPI_DOUBLE_PRECISION
oldtypes(17:21) = MPI_INTEGER
oldtypes(22:25)    = MPI_DOUBLE_PRECISION
oldtypes(26:29) = MPI_INTEGER
oldtypes(30:47) = MPI_DOUBLE_PRECISION

call mpi_get_address(obj%k,offsets(1),ierr)
call mpi_get_address(obj%voro,offsets(2),ierr)
call mpi_get_address(obj%voroidx,offsets(3),ierr)
call mpi_get_address(obj%par,offsets(4),ierr)
call mpi_get_address(obj%hiface,offsets(5),ierr)
call mpi_get_address(obj%ziface,offsets(6),ierr)
call mpi_get_address(obj%sdparR,offsets(7),ierr)
call mpi_get_address(obj%sdparV,offsets(8),ierr)
call mpi_get_address(obj%sdparT,offsets(9),ierr)
call mpi_get_address(obj%sdparSWD,offsets(10),ierr)
call mpi_get_address(obj%sdaveH,offsets(11),ierr)
call mpi_get_address(obj%sdaveV,offsets(12),ierr)
call mpi_get_address(obj%sdaveT,offsets(13),ierr)
call mpi_get_address(obj%sdaveSWD,offsets(14),ierr)
call mpi_get_address(obj%arpar,offsets(15),ierr)
call mpi_get_address(obj%arparSWD,offsets(16),ierr)
call mpi_get_address(obj%idxar,offsets(17),ierr)
call mpi_get_address(obj%idxarSWD,offsets(18),ierr)
call mpi_get_address(obj%gvoroidx,offsets(19),ierr)
call mpi_get_address(obj%nunique,offsets(20),ierr)
call mpi_get_address(obj%NFP,offsets(21),ierr)
call mpi_get_address(obj%beta,offsets(22),ierr)
call mpi_get_address(obj%logL,offsets(23),ierr)
call mpi_get_address(obj%logPr,offsets(24),ierr)
call mpi_get_address(obj%tcmp,offsets(25),ierr)
call mpi_get_address(obj%ireject_bd,offsets(26),ierr)
call mpi_get_address(obj%iaccept_bd,offsets(27),ierr)
call mpi_get_address(obj%ireject_bds,offsets(28),ierr)
call mpi_get_address(obj%iaccept_bds,offsets(29),ierr)
call mpi_get_address(obj%DobsR,offsets(30),ierr)
call mpi_get_address(obj%DpredR,offsets(31),ierr)
call mpi_get_address(obj%DobsV,offsets(32),ierr)
call mpi_get_address(obj%DpredV,offsets(33),ierr)
call mpi_get_address(obj%DobsT,offsets(34),ierr)
call mpi_get_address(obj%DpredT,offsets(35),ierr)
call mpi_get_address(obj%S,offsets(36),ierr)
call mpi_get_address(obj%DresR,offsets(37),ierr)
call mpi_get_address(obj%DresV,offsets(38),ierr)
call mpi_get_address(obj%DresT,offsets(39),ierr)
call mpi_get_address(obj%DarR,offsets(40),ierr)
call mpi_get_address(obj%DarV,offsets(41),ierr)
call mpi_get_address(obj%DarT,offsets(42),ierr)
!! SWD data
call mpi_get_address(obj%DobsSWD,offsets(43),ierr)
call mpi_get_address(obj%DpredSWD,offsets(44),ierr)
call mpi_get_address(obj%DresSWD,offsets(45),ierr)
call mpi_get_address(obj%DarSWD,offsets(46),ierr)
call mpi_get_address(obj%periods,offsets(47),ierr)

DO ifield=2,SIZE(offsets)
  offsets(ifield) = offsets(ifield) - offsets(1)
ENDDO
offsets(1) = 0
!IF(rank == src)PRINT*,'offsets:',offsets

!  Now define structured type and commit it
call MPI_TYPE_CREATE_STRUCT( NFIELD, blockcounts, offsets, oldtypes, objtype,ierr)
call MPI_TYPE_COMMIT(objtype, ierr)

RETURN
END SUBROUTINE MAKE_MPI_STRUC_SP
!!==============================================================================

SUBROUTINE ALLOC_OBJ(obj)
!!==============================================================================
!!
!! Allocates memory.
!!
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc) :: obj

ALLOCATE( obj%voro(NLMX,NPL),obj%voroidx(NLMX,NPL),obj%par((NLMX+1)*NPL*NPL) )
ALLOCATE( obj%hiface(NLMX*NPL),obj%ziface(NLMX*NPL) )
ALLOCATE( obj%sdparR(NRF1),obj%sdparV(NRF1),obj%sdparT(NRF1),obj%arpar(3*NRF1),obj%idxar(3*NRF1) )
ALLOCATE( obj%sdaveH(NRF1),obj%sdaveV(NRF1),obj%sdaveT(NRF1),obj%sdaveSWD(NMODE) )
ALLOCATE( obj%sdparSWD(NMODE),obj%arparSWD(NMODE),obj%idxarSWD(NMODE) )
ALLOCATE( obj%gvoroidx(NPL-1) )
!! R, T, and V seismograms:
ALLOCATE( obj%DobsR(NRF1,NTIME2),obj%DpredR(NRF1,NTIME) )
ALLOCATE( obj%DobsV(NRF1,NTIME2),obj%DpredV(NRF1,NTIME) )
ALLOCATE( obj%DobsT(NRF1,NTIME2),obj%DpredT(NRF1,NTIME) )
ALLOCATE( obj%S(NRF1,NSRC) )
ALLOCATE( obj%DresR(NRF1,NTIME),obj%DresV(NRF1,NTIME),obj%DresT(NRF1,NTIME) )
ALLOCATE( obj%DarR(NRF1,NTIME),obj%DarT(NRF1,NTIME),obj%DarV(NRF1,NTIME) )
!! SWD data:
ALLOCATE( obj%DobsSWD(NMODE,NDAT_SWD),obj%DpredSWD(NMODE,NDAT_SWD),obj%DresSWD(NMODE,NDAT_SWD) )
ALLOCATE( obj%DarSWD(NMODE,NDAT_SWD),obj%periods(NMODE,NDAT_SWD) )

obj%voro     = 0._RP
obj%voroidx  = 0
obj%par      = 0._RP
obj%hiface   = 0._RP
obj%ziface   = 0._RP
obj%sdparR   = 0._RP
obj%sdparV   = 0._RP
obj%sdparT   = 0._RP
obj%sdparSWD = 0._RP
obj%sdaveH   = 0._RP
obj%sdaveV   = 0._RP
obj%sdaveT   = 0._RP
obj%sdaveSWD = 0._RP
obj%arpar    = 0._RP
obj%idxar    = 0
obj%gvoroidx = 0
obj%DobsR    = 0._RP
obj%DpredR   = 0._RP
obj%DobsV    = 0._RP
obj%DpredV   = 0._RP
obj%DobsT    = 0._RP
obj%DpredT   = 0._RP
obj%S        = 0._RP
obj%DresR    = 0._RP
obj%DresV    = 0._RP
obj%DresT    = 0._RP
obj%DarR     = 0._RP
obj%DarV     = 0._RP
obj%DarT     = 0._RP
obj%DobsSWD  = 0._RP
obj%DpredSWD = 0._RP
obj%DresSWD  = 0._RP
obj%DarSWD   = 0._RP
obj%periods  = 0._RP
END SUBROUTINE ALLOC_OBJ
!=======================================================================
!EOF

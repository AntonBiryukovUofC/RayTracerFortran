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
!!               sdparRT   sdaveRT
                 NMODE,    NMODE,&
!!               arpar arparRT idxar idxarRT  gvoroidx nunique NFP
                 NMODE, NMODE,  NMODE,  NMODE,     NPL-1,     1,    1, & 
!!               beta    logL     logPr tcmp  ireject_bd  iaccept_bd ireject_bds  iaccept_bds 
                  1,       1,       1,    1,      1,          1,         1,          1,&
!!               DobsR           DpredR        DobsV       DpredV         DobsT        DpredT
                ! NRF1*NTIME2,   NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME,  NRF1*NTIME2,  NRF1*NTIME, &
!!                S            DresR         DresV        DresT        DarR        DarV        DarT
                ! NRF1*NSRC, NRF1*NTIME,   NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME,  NRF1*NTIME, &
!!                DobsRT           DpredRT        DresRT       DarRT           periods
                 NMODE*NDAT_RT, NMODE*NDAT_RT, NMODE*NDAT_RT, NMODE*NDAT_RT, NMODE*NDAT_RT /)


!! Oldtypes
oldtypes(1)     = MPI_INTEGER
oldtypes(2)     = MPI_DOUBLE_PRECISION
oldtypes(3)     = MPI_INTEGER
oldtypes(4:10)   = MPI_DOUBLE_PRECISION
oldtypes(11:15) = MPI_INTEGER
oldtypes(16:19)    = MPI_DOUBLE_PRECISION
oldtypes(20:23) = MPI_INTEGER
oldtypes(24:27) = MPI_DOUBLE_PRECISION

call mpi_get_address(obj%k,offsets(1),ierr)
call mpi_get_address(obj%voro,offsets(2),ierr)
call mpi_get_address(obj%voroidx,offsets(3),ierr)
call mpi_get_address(obj%par,offsets(4),ierr)
call mpi_get_address(obj%hiface,offsets(5),ierr)
call mpi_get_address(obj%ziface,offsets(6),ierr)
call mpi_get_address(obj%sdparRT,offsets(7),ierr)
call mpi_get_address(obj%sdaveRT,offsets(8),ierr)
call mpi_get_address(obj%arpar,offsets(9),ierr)
call mpi_get_address(obj%arparRT,offsets(10),ierr)
call mpi_get_address(obj%idxar,offsets(11),ierr)
call mpi_get_address(obj%idxarRT,offsets(12),ierr)
call mpi_get_address(obj%gvoroidx,offsets(13),ierr)
call mpi_get_address(obj%nunique,offsets(14),ierr)
call mpi_get_address(obj%NFP,offsets(15),ierr)
call mpi_get_address(obj%beta,offsets(16),ierr)
call mpi_get_address(obj%logL,offsets(17),ierr)
call mpi_get_address(obj%logPr,offsets(18),ierr)
call mpi_get_address(obj%tcmp,offsets(19),ierr)
call mpi_get_address(obj%ireject_bd,offsets(20),ierr)
call mpi_get_address(obj%iaccept_bd,offsets(21),ierr)
call mpi_get_address(obj%ireject_bds,offsets(22),ierr)
call mpi_get_address(obj%iaccept_bds,offsets(23),ierr)
!! RT data
call mpi_get_address(obj%DobsRT,offsets(24),ierr)
call mpi_get_address(obj%DpredRT,offsets(25),ierr)
call mpi_get_address(obj%DresRT,offsets(26),ierr)
call mpi_get_address(obj%DarRT,offsets(27),ierr)

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

ALLOCATE( obj%voro(NLMX,NPL),obj%voroidx(NLMX,NPL))
ALLOCATE(obj%par((NLMX+1)*NPL*NPL) )
ALLOCATE( obj%hiface(NLMX*NPL),obj%ziface(NLMX*NPL) )
print *,' yaii'
ALLOCATE( obj%sdaveRT(NMODE) )

ALLOCATE( obj%sdparRT(NMODE),obj%arparRT(NMODE),obj%arpar(NMODE))
ALLOCATE(obj%idxarRT(NMODE),obj%idxar(NMODE) )
ALLOCATE( obj%gvoroidx(NPL-1) )
!! RT data:

ALLOCATE( obj%DobsRT(NMODE,NDAT_RT),obj%DpredRT(NMODE,NDAT_RT),obj%DresRT(NMODE,NDAT_RT) )
!ALLOCATE( obj%DarRT(NMODE,NDAT_RT),obj%periods(NMODE,NDAT_RT) )
ALLOCATE( obj%DarRT(NMODE,NDAT_RT) )

obj%voro     = 0._RP
obj%voroidx  = 0
obj%par      = 0._RP
obj%hiface   = 0._RP
obj%ziface   = 0._RP
obj%sdparRT = 0._RP
obj%sdaveRT = 0._RP
obj%arpar    = 0._RP
obj%arparRT    = 0._RP
obj%idxar    = 0
obj%idxarRT    = 0

obj%gvoroidx = 0
print *,' yaee'

END SUBROUTINE ALLOC_OBJ
!=======================================================================
!EOF

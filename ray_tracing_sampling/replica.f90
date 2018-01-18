!==============================================================================
!
!  Reversible Jump MCMC Sampling with parallel tempering for 
!  reflection coefficient inversion (including Buckingham porous media model)
!
!------------------------------------------------------------------------------
!
!  Jan Dettmer, University of Victoria, May 9 2014
!  jand@uvic.ca                       (250) 472 4026
!  http://web.uvic.ca~/jand/
!  Last change: May 9 2014
!
!  Based on Green 1995, Malinverno 2002, Bodin Sambridge 2009, 
!  Agostinetti Malinverno 2010, Dettmer etal 2010, 2012, 1013
!
!==============================================================================

PROGRAM REPLICA

!=======================================================================
USE RJMCMC_COM
USE NR
!USE SAC_I_O
USE ieee_arithmetic
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: NBURNIN,NSMP,NSUB,NSAVE,ntmp
INTEGER(KIND=IB)  :: iglob,ismp,NSMP2,ivo,ipar,ipar2,isub,iaz,idat,io

TYPE (objstruc)                         :: obj      ! Object on slave node
!TYPE (covstruc),ALLOCATABLE,DIMENSION(:):: cov      ! Structure for inverse data Covariance Matrices
REAL(KIND=RP),DIMENSION(:),ALLOCATABLE  :: tmp
REAL(KIND=RP)                           :: vref

CHARACTER(len=64):: plotparfile       = 'plotpar.txt'

OPEN(UNIT=20,FILE=filebasefile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)

repfileR     = filebase(1:filebaselen) // '_predR.txt'
repfileT     = filebase(1:filebaselen) // '_predT.txt'
repfileRF    = filebase(1:filebaselen) // '_predRF.txt'
repfileV     = filebase(1:filebaselen) // '_predV.txt'
repfileS     = filebase(1:filebaselen) // '_predS.txt'
repfileSWD   = filebase(1:filebaselen) // '_predSWD.txt'
repfileSWDar = filebase(1:filebaselen) // '_predSWDar.txt'

!!
!! Read parameter file and allocate global prior limits:
!!
CALL READPARFILE()

CALL ALLOC_OBJ(obj)

OPEN(UNIT=32,FILE=plotparfile,FORM='formatted',ACTION='READ')
READ(32,*) NBURNIN
READ(32,*) NSMP
READ(32,*) NSUB
READ(32,*) NSAVE
CLOSE(32)

ncount1 = NFPMX+NAP+3*NRF1+3*NRF1+2*NMODE
ncount2 = NFPMX+1+3*NRF1+3*NRF1+2*NMODE
ALLOCATE( sample(1,ncount1),tmp(ncount1) )
ALLOCATE( tmpmap(ncount2) )
buf_save_snd = 0._RP
buf_save_rcv = 0._RP
sample       = 0._RP
WRITE(6,*) 'SHAPE    = ',SHAPE(sample)
WRITE(6,*) ''

OPEN(UNIT=33,FILE=samplefile,FORM='formatted',ACTION='READ')
iglob = 1
DO ismp = 1,NBURNIN
   READ(33,*) tmp
   iglob = iglob + 1
ENDDO
WRITE(6,*) 'Discarded burn-in.'

IF(I_RV  == 1)THEN
  OPEN(UNIT=50,FILE=repfileR,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=1024)
ELSE
  OPEN(UNIT=50,FILE=repfileRF,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=1024)
ENDIF
OPEN(UNIT=51,FILE=repfileV,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
OPEN(UNIT=52,FILE=repfileS,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
OPEN(UNIT=53,FILE=repfileSWD,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
OPEN(UNIT=54,FILE=repfileSWDar,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)
IF(I_T   == 1)THEN
  OPEN(UNIT=55,FILE=repfileT,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=1024)
ENDIF

NSMP2 = FLOOR(REAL(NSMP-iglob,RP)/REAL(NSUB,RP))
WRITE(6,*)   'NBURNIN    = ', NBURNIN
WRITE(6,*)   'NSMP       = ', NSMP
WRITE(6,*)   'NSMP2      = ', NSMP2
WRITE(6,*)   'NSUB       = ', NSUB
WRITE(6,*)   'NSAVE      = ', NSAVE
!------------------------------------------------------------------------
!  Read in data
!------------------------------------------------------------------------
CALL READDATA(obj)

WRITE(6,*)'sampling_dt = ', sampling_dt
WRITE(6,*)'NTIME  = ', NTIME
WRITE(6,*)'NTIME2 = ', NTIME2
WRITE(6,*)'NSRC   = ', NSRC

!!
!! Read mapfile to start:
!!
OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) tmpmap
CLOSE(20)

obj%k       = INT(tmpmap(1),IB)
obj%NFP     = (obj%k * NPL)
obj%voro    = 0._RP
obj%voroidx = 0
obj%sdparR  = 0._RP
obj%sdparV  = 0._RP
obj%arpar   = 0._RP
obj%sdparSWD= 0._RP
obj%arparSWD= 0._RP

DO ivo = 1,obj%k
  ipar = (ivo-1)*NPL+2
  obj%voro(ivo,:) = tmpmap(ipar:ipar+NPL-1)
  obj%voroidx(ivo,:) = 1
  DO ipar2 = 1,NPL
    IF(obj%voro(ivo,ipar2) < -99._RP)THEN
      obj%voroidx(ivo,ipar2) = 0
    ELSE
      obj%voroidx(ivo,ipar2) = 1
    ENDIF
  ENDDO
ENDDO
IF(ICOV == 1) obj%sdparR   = tmpmap(NFPMX+1+1:NFPMX+1+NRF1)
IF(ICOV == 1) obj%sdparV   = tmpmap(NFPMX+1+1+NRF1:NFPMX+1+2*NRF1)
IF(ICOV == 1) obj%sdparT   = tmpmap(NFPMX+1+1+2*NRF1:NFPMX+1+3*NRF1)
IF(ICOV == 1) obj%sdparSWD = tmpmap(NFPMX+1+1+3*NRF1:NFPMX+1+3*NRF1+NMODE)
IF(IAR == 1)THEN
  obj%arpar = tmpmap(NFPMX+1+3*NRF1+NMODE+1:NFPMX+3*NRF1+NMODE+3*NRF1+1)
  obj%arparSWD = tmpmap(NFPMX+1+3*NRF1+NMODE+3*NRF1+1:NFPMX+3*NRF1+NMODE+3*NRF1+NMODE+1)
  obj%idxar = 1
  DO ipar = 1,NRF1
    IF(obj%arpar(ipar) < minlimar(ipar)) obj%idxar(ipar) = 0
  ENDDO
  obj%idxarSWD = 1
  DO ipar = 1,NMODE
    IF(obj%arparSWD(ipar) < minlimarSWD(ipar)) obj%idxarSWD(ipar) = 0
  ENDDO
ENDIF
!CALL INTERPLAYER_novar(obj)
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(obj)
IF(I_VARPAR == 1)CALL INTERPLAYER(obj)

CALL LOGLHOOD(obj,1)
CALL PRINTPAR(obj)
WRITE(6,*) obj%logL

WRITE(6,*) 'Computing replica',NSMP2,'samples'
DO ismp = 1,NSMP2-1
  READ(33,*) sample(1,:)
  iglob = iglob + 1
  DO isub = 1,NSUB-1
    READ(33,*) tmp
    iglob = iglob + 1
  ENDDO

  obj%k        = INT(sample(1,4))
  obj%NFP      = (obj%k * NPL) + (NPL-1)
  obj%voro     = 0._RP
  obj%voroidx  = 0
  obj%sdparR   = 0._RP
  obj%sdparV   = 0._RP
  obj%arpar    = 0._RP
  obj%sdparSWD = 0._RP
  obj%arparSWD = 0._RP

  DO ivo = 1,obj%k
    ipar = (ivo-1)*NPL+5
    obj%voro(ivo,:) = sample(1,ipar:ipar+NPL-1)
    DO ipar2 = 1,NPL
      IF(obj%voro(ivo,ipar2) < -99._RP)THEN
        obj%voroidx(ivo,ipar2) = 0
      ELSE
        obj%voroidx(ivo,ipar2) = 1
      ENDIF
    ENDDO
  ENDDO
  IF(ICOV == 1) obj%sdparR  = sample(1,NFPMX+1+4:NFPMX+4+NRF1)
  IF(ICOV == 1) obj%sdparV  = sample(1,NFPMX+1+4+NRF1:NFPMX+4+2*NRF1)
  IF(ICOV == 1) obj%sdparT  = sample(1,NFPMX+1+4+2*NRF1:NFPMX+4+3*NRF1)
  IF(ICOV == 1) obj%sdparSWD= sample(1,NFPMX+1+4+3*NRF1:NFPMX+4+3*NRF1+NMODE)
  IF(IAR == 1)THEN
    obj%arpar = sample(1,NFPMX+4+3*NRF1+NMODE+1:NFPMX+4+6*NRF1+NMODE)
    obj%arparSWD = sample(1,NFPMX+4+6*NRF1+NMODE+1:NFPMX+4+6*NRF1+2*NMODE)
    obj%idxar = 1
    DO ipar = 1,NRF1
      IF(obj%arpar(ipar) < minlimar(ipar)) obj%idxar(ipar) = 0
    ENDDO
    obj%idxarSWD = 1
    DO ipar = 1,NMODE
      IF(obj%arparSWD(ipar) < minlimar(ipar)) obj%idxarSWD(ipar) = 0
    ENDDO
  ENDIF

  !CALL INTERPLAYER_novar(obj)
  IF(I_VARPAR == 0)CALL INTERPLAYER_novar(obj)
  IF(I_VARPAR == 1)CALL INTERPLAYER(obj)

  CALL LOGLHOOD(obj,1)
  !CALL PRINTPAR(obj)
  !PRINT*,''
  !PRINT*,''
  IF(MOD(ismp,100)==0)PRINT*,'logL=',sample(1,1),obj%logL
  CALL SAVEREPLICA(obj)
  IF(MOD(ismp,100)==0)WRITE(6,*)ismp

  IF(iglob >= NSMP)EXIT
ENDDO
CLOSE(49)

CALL FLUSH(6)

201 FORMAT(200F12.4)
202 FORMAT(I8,2F16.6,I8,F16.6,2I8,1F13.3,I8,1F13.3)
203 FORMAT(A119)
204 FORMAT(A56)
205 FORMAT(10F10.4)
209 FORMAT(A26,A40)
210 FORMAT(A26,20I4)
211 FORMAT(20000ES12.4)
212 FORMAT(I4,A8,F12.4,A5,I3)
213 FORMAT(I4,I6,A8,F12.4,A5,I3,A14,F12.4,A11,I3)
214 FORMAT(a21,100I4)
215 FORMAT(a21,F8.4)
216 FORMAT(a21,F8.4)
217 FORMAT(a21,F8.4)

END PROGRAM REPLICA
!=======================================================================

SUBROUTINE SAVEREPLICA(obj)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
TYPE (objstruc)  :: obj      ! Best object

IF(I_RV >= 0)THEN
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredR(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(51,208) obj%DpredV(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(52,208) obj%S(i,:)
  ENDDO
  IF(I_T == 1)THEN
    DO i = 1,NRF1
      WRITE(55,208) obj%DpredT(i,:)
    ENDDO
  ENDIF
ELSE
  DO i = 1,NRF1
    WRITE(50,208) obj%DpredR(i,:)
  ENDDO
ENDIF

IF(I_SWD == 1)THEN
  DO i = 1,NMODE
    WRITE(53,208) obj%DpredSWD(i,:)
  ENDDO
  DO i = 1,NRF1
    WRITE(54,208) obj%DarSWD(i,:)
  ENDDO
ENDIF 

208 FORMAT(5000ES20.10)
RETURN
END SUBROUTINE SAVEREPLICA
!==============================================================================
function factorial (n) result (res)
!==============================================================================
 
implicit none
integer, intent (in) :: n
integer :: res
integer :: i
res = product ((/(i, i = 1, n)/))
end function factorial
!==============================================================================
function choose (n, k) result (res)
!==============================================================================
implicit none
integer, intent (in) :: n
integer, intent (in) :: k
integer :: res,factorial
res = factorial (n) / (factorial (k) * factorial (n - k))
end function choose
!=======================================================================
SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
!!=======================================================================
!!
!! 
!!
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: Nx,Ny,Nz
INTEGER(KIND=IB):: irow,icol1,icol2,ix1,ix2
REAL(KIND=RP)   :: x(Nx), y(Ny), z(Nz),xr(Nx)
REAL(KIND=RP)   :: XMAT(Nz,Ny)

XMAT = 0._RP
xr = x(Nx:1:-1)
DO irow=1,Nz
  icol1 = MAX(1,irow-Nx+1)
  icol2 = MIN(Ny,irow)
  ix1   = MAX(1,Nx-irow+1)
  ix2   = MIN(Nz-irow+1,Nx)
  XMAT(irow,icol1:icol2) = xr(ix1:ix2)
ENDDO

z = MATMUL(XMAT,y)
RETURN
END SUBROUTINE CONVT
!=======================================================================

SUBROUTINE DCONVT(z,x,y,Nz,Nx,Ny)
!!=======================================================================
!!
!! Time-domain deconvolution in matrix formulation:
!!    If z=x*y is the  Nz=Nx+Ny-1 length 
!!    convolution, compute y by deconvolving z 
!!    by x. Considering the convolution via
!!    matrix multiplication z=Xy, where X is an
!!    Nz by Ny matrix, the deconvolution is 
!!    carried out using the pseudo-inverse of X.
!!    (S. Dosso, August 2014) 
!!
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: Nx,Ny,Nz
INTEGER(KIND=IB):: irow,icol1,icol2,ix1,ix2
REAL(KIND=RP)   :: x(Nx), y(Ny), z(Nz), b(Nz),xr(Nx)
REAL(KIND=RP)   :: XMAT(Nz,Ny)
!! Lapack variables:
INTEGER(KIND=IB)          :: MRANK, LWORK, INFO,LDA, LDB
INTEGER(KIND=IB),PARAMETER:: LWMAX = 4000, NRHS = 1
REAL(KIND=RP),PARAMETER   :: RCOND = 1.e-12_RP
REAL(KIND=RP)             :: SV(Ny),WORK(LWMAX)

b = z
LDA = Nz
LDB = Nz
!!
!! Build linear system of equations:
XMAT = 0._RP
xr = x(Nx:1:-1)
DO irow=1,Nz
  icol1 = MAX(1,irow-Nx+1)
  icol2 = MIN(Ny,irow)
  ix1   = MAX(1,Nx-irow+1)
  ix2   = MIN(Nz-irow+1,Nx)
  XMAT(irow,icol1:icol2) = xr(ix1:ix2)
ENDDO

!! In Matlab, can use Moore-Penrose pseudo inverse:
!y = pinv(X)*z;
!! Here, use Lapack to solve LLS via SGELSS:
!! SGELSS computes the minimum norm solution to a real linear least
!! squares problem:
!!
!! Minimize 2-norm(| b - A*x |).
!!
!! using the singular value decomposition (SVD) of A. A is an M-by-N
!! matrix which may be rank-deficient.
!!
!! Several right hand side vectors b and solution vectors x can be
!! handled in a single call; they are stored as the columns of the
!! M-by-NRHS right hand side matrix B and the N-by-NRHS solution matrix
!! X.
!!
!! The effective rank of A is determined by treating as zero those
!! singular values which are less than RCOND times the largest singular
!! value.
!! subroutine sgelss(integer M,integer N,integer NRHS,real, dimension( lda, * ) A,
!!                   integer LDA,real, dimension( ldb, * ) B,integer LDB,
!!                   real, dimension( * ) S,real RCOND,integer RANK,
!!                   real, dimension( * ) WORK,integer LWORK,integer INFO )

!! first, figure out optimal work length...
LWORK = -1
CALL DGELSS(Nz, Ny, NRHS, XMAT, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!! Carry out DGELSS with optimal length
CALL DGELSS(Nz, Ny, NRHS, XMAT, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
IF(INFO /= 0)WRITE(*,*)'WARNING DGELSS unstable!'

y = b(1:Ny)
RETURN
END SUBROUTINE DCONVT
!=======================================================================
! This is the end my fiend...
! EOF

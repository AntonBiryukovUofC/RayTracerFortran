!==============================================================================

SUBROUTINE READPARFILE()
!==============================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=RP):: ntmp,iaz,ik
REAL(KIND=RP)   :: vref, dVs, VpVsmin, VpVsmax, dVpVs
INTERFACE
   FUNCTION LOGFACTORIAL(n)
     USE DATA_TYPE
     REAL(KIND=RP) :: LOGFACTORIAL
     REAL(KIND=RP),INTENT(IN):: n
   END FUNCTION LOGFACTORIAL
END INTERFACE

parfile        = filebase(1:filebaselen) // '_parameter.dat'!! Inversion parameter file
!!
!! Read parameter file
!!
OPEN(UNIT=20,FILE=parfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) IMAP
READ(20,*) ICOV
READ(20,*) ENOS       !! Even numbered order statistic on k (avoids thin layers)
READ(20,*) IPOIPR     !! Applies Poisson prior on k
READ(20,*) IAR
READ(20,*) NPL !! Number of parameters per layer
READ(20,*) I_VARPAR
READ(20,*) IBD_SINGLE
READ(20,*) I_RT      !! Invert RT data
READ(20,*) I_VREF
READ(20,*) I_VPVS
READ(20,*) ISMPPRIOR
READ(20,*) ISETSEED
READ(20,*) IEXCHANGE
READ(20,*) NDAT_RT   !! No. RT data
READ(20,*) NMODE      !! No. RT modes (Phases in our case)
READ(20,*) NSRC       !! No. of raypaths computed 
READ(20,*) NLMN       !! Max number of layers
READ(20,*) NLMX       !! Max number of layers
READ(20,*) ICHAINTHIN !! Chain thinning interval
READ(20,*) NKEEP      !! No. models to keep before writing
READ(20,*) NPTCHAINS1 !! Chain thinning interval
READ(20,*) dTlog      !! Temperature increment
READ(20,*) lambda     !! Lambda for Poisson prior on k
READ(20,*) hmx        !! Max. partition depth
READ(20,*) hmin       !! Min. layer thickness (must be small enough to not violate detailed balance)
print *,'1st'
READ(20,*) TCHCKPT    !! Checkpointing interval in s
READ(20,*) dVs        !! Vs one sided prior width (relative to background model)
READ(20,*) dVpVs      !! VpVs ratio one sided prior width
print *,'2ndd'
READ(20,*) sdmn       !! data (residual) error standard deviation prior lower limit
print *,'HI'
READ(20,*) sdmx       !! data (residual) error standard deviation prior upper limit
print *,'3rdd'
READ(20,*) VpVsmin   !! minimum VpVs ratio
print *,'4th'
READ(20,*) VpVsmax   !! maximum VpVs ratio
print *,'Done!'

CLOSE(20)

!! Allocate raysum oarameters
!! (these parameters are also used in ray3d)
!CALL ALLOC_RAYSUM()
!! Read geometry file:
!CALL readgeom(modname,baz2,slow2,sta_dx,sta_dy,ntr)
  
!IF(iraysum /= 1)THEN
!  slow2 = slow2 * 1000._RP
!  baz2 = baz2 * 180._RP / PI2
!ENDIF

!kmin = NLMN
!kmax = NLMX
!! HARDWIRE TO FIXED-D PROBLEM
!! k always refers to No. layers (halfspace is extra) 
kmin = 8
kmax = 8

!! Poisson Prior on k:
ALLOCATE(pk(kmax))
DO ik = kmin,kmax
  pk(ik)  = EXP(-lambda)*lambda**REAL(ik,RP)/EXP(LOGFACTORIAL(REAL(ik,RP)))
ENDDO

!infileV        = filebase(1:filebaselen) // '_V_b.txt'
!IF(I_RV == 1)THEN
!  infileH      = filebase(1:filebaselen) // '_R_b.txt'
!ELSEIF(I_RV == -1)THEN
!  infileH      = filebase(1:filebaselen) // '_RF.txt'
!ENDIF
!infileT        = filebase(1:filebaselen) // '_T_b.txt'
infileRT      = filebase(1:filebaselen) // '_RT.txt'
infileref      = filebase(1:filebaselen) // '_vel_ref.txt'
logfile        = filebase(1:filebaselen) // '_RJMH.log'
seedfile       = filebase(1:filebaselen) // '_seeds.log'
mapfile        = filebase(1:filebaselen) // '_map.dat'
arfile         = filebase(1:filebaselen) // '_ar.dat'
!predfile       = filebase(1:filebaselen) // '_mappred.dat'
!obsfile        = filebase(1:filebaselen) // '_obs.dat'
arfileRT      = filebase(1:filebaselen) // '_maparRT.dat'
predfileRT    = filebase(1:filebaselen) // '_mappredRT.dat'
obsfileRT     = filebase(1:filebaselen) // '_obsRT.dat'
!covfile        = filebase(1:filebaselen) // '_cov.txt'
sdfile         = filebase(1:filebaselen) // '_sigma.txt'
samplefile     = filebase(1:filebaselen) // '_voro_sample.txt'
stepsizefile     = filebase(1:filebaselen) // '_stepsize.txt'
source_data_file  = filebase(1:filebaselen) // '_src_data.txt'



!IF(IDIP == 0) NPL=3      ! No. parameters per layer
!IF(IDIP == 1) NPL=4      ! No. parameters per layer with dip
!IF(IDIP == 2) NPL=5      ! No. parameters per layer with strike and dip

!NTIME2 = NTIME+NSRC-1  ! Number time samples for zero padded observations
!NRRG   = NTIME2 + NTIME - 1 ! No. points for convolution RRG
!NRGRG  = NTIME  + NTIME - 1 ! No. points for convolution RGRG
NFPMX  = NLMX * NPL
NFPMX2 = NLMX * NPL * NPL

ioutside = 0;ireject  = 0; iaccept = 0; iaccept_delay = 0; ireject_delay = 0
i_sdpert = 0;ishearfail = 0 ;i_ref_nlay = 0

!!
!!
!! Read velocity reference file
!!
201 FORMAT(a64)
202 FORMAT(a28,I4,I4)
203 FORMAT(4F12.3)
IF(I_VREF == 1)THEN
  OPEN(UNIT=20,FILE=infileref,FORM='formatted',STATUS='OLD',ACTION='READ')
  READ(20,*) NVELREF, ntmp
  NPREM = ntmp - NVELREF
  ALLOCATE(vel_ref(4,NVELREF),vel_prem(4,NPREM))
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
  IF(rank == src)WRITE(6,202) ' Velocity reference model:  ',NVELREF,NPREM
  IF(rank == src)WRITE(6,201) ' Reference:                                                     '
  IF(rank == src)WRITE(6,201) '   Depth(km)   Vs (km/s)      VpVs        Density               '
  !!  Read reference model to max sampling depth
  DO iaz = 1,NVELREF
    READ(20,*) vel_ref(:,iaz)
    IF(rank == src)WRITE(6,203) vel_ref(:,iaz)
  ENDDO
  !!  Read PREM beyond that
  IF(rank == src)WRITE(6,201) ' Prem (deep reference):                                         '
  DO iaz = 1,NPREM
    READ(20,*) vel_prem(:,iaz)
    IF(rank == src)WRITE(6,203) vel_prem(:,iaz)
  ENDDO
  CLOSE(20)
  IF(rank == src)WRITE(6,201) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ '
ENDIF

!!
!! Set width of Gaussian for raysum. -1. means get impulse response.
!!
!IF(I_RV >= 0)THEN
!  width2 = -1._SP
!! The width for the RF case is read from input parameter file
!ELSE
!  width = 2._SP
!ENDIF

ALLOCATE( sdbuf(3,1,NBUF) )  ! NRF1 replaced with 1
ALLOCATE(sdevm((NLMX*NPL)+NPL-1,NLMX))
sdevm = 0._RP

ALLOCATE(minlim(NPL),maxlim(NPL),maxpert(NPL),pertsd(NPL),pertsdsc(NPL))
minlim   = 0._RP;maxlim   = 0._RP;maxpert  = 0._RP
pertsd   = 0._RP;pertsdsc = 30._RP

ALLOCATE(minlimar(3*1),maxlimar(3*1),maxpertar(3*1)) ! NRF1 replaced with 1
ALLOCATE(pertarsd(3*1),pertarsdsc(3*1))! NRF1 replaced with 1
minlimar  = 0._RP;maxlimar  = 0._RP;maxpertar = 0._RP
pertarsd  = 0._RP;pertarsdsc= 18._RP

ALLOCATE(minlimarRT(NMODE),maxlimarRT(NMODE),maxpertarRT(NMODE))
ALLOCATE(pertarsdRT(NMODE),pertarsdscRT(NMODE))
minlimarRT  = 0._RP;maxlimarRT  = 0._RP;maxpertarRT = 0._RP
pertarsdRT  = 0._RP;pertarsdscRT= 18._RP

ALLOCATE(minlimsd(1),maxlimsd(1),maxpertsd(1))
ALLOCATE(pertsdsd(1),pertsdsdsc(1))
minlimsd  = 0._RP;maxlimsd  = 0._RP;maxpertsd = 0._RP
pertsdsd  = 0._RP;pertsdsdsc= 18._RP

ALLOCATE(minlimsdRT(NMODE),maxlimsdRT(NMODE),maxpertsdRT(NMODE))
ALLOCATE(pertsdsdRT(NMODE),pertsdsdscRT(NMODE))
minlimsdRT  = 0._RP;maxlimsdRT  = 0._RP;maxpertsdRT = 0._RP
pertsdsdRT  = 0._RP;pertsdsdscRT= 18._RP

! Allocate the source-related arrays:
ALLOCATE(src_offset(NSRC))
ALLOCATE(src_depth(NSRC))

!!
!!  Prior bounds
!! (Note: Density is empirical through Birch's Law)
!!
!! Without dip:
!!           h     vs     vp/vs
IF(I_VPVS == 1)THEN
  !! Sample Vs and VpVs ratio
  minlim(1:3) = (/ hmin, -dVs, -dVpVs/)
  maxlim(1:3) = (/ hmx,   dVs,  dVpVs/)
ELSEIF(I_VPVS == 0)THEN
  !! Sample Vs and VpVs ratio
  minlim(1:2) = (/ hmin, 1500._RP/)
  maxlim(1:2) = (/ hmx,   10000._RP/)
ENDIF



maxpert = maxlim-minlim
pertsd = maxpert/pertsdsc

IF(IAR == 1)THEN
  !! Set prior and proposal scaling for AR model:
  minlimar   = -0.5000_RP
  maxlimar   =  0.90_RP
  pertarsdsc =  10._RP
  maxpertar  = maxlimar-minlimar
  pertarsd   = maxpertar/pertarsdsc
  minlimarRT   = -0.5000_RP
  maxlimarRT   =  0.90_RP
  pertarsdscRT =  10._RP
  maxpertarRT  = maxlimarRT-minlimarRT
  pertarsdRT   = maxpertarRT/pertarsdscRT
ENDIF

IF(ICOV == 1)THEN
  !! Set prior and proposal scaling for data error standard deviations:
  minlimsd   = sdmn(1)
  maxlimsd   = sdmx(1)
  pertsdsdsc = 10._RP
  maxpertsd  = maxlimsd-minlimsd
  pertsdsd   = maxpertsd/pertsdsdsc
  minlimsdRT   = sdmn(2)
  maxlimsdRT   = sdmx(2)
  pertsdsdscRT = 10._RP
  maxpertsdRT  = maxlimsdRT-minlimsdRT
  pertsdsdRT   = maxpertsdRT/pertsdsdscRT
ENDIF

!!
!! Write some info:
!!
IF(rank == src)THEN
  WRITE(6,*) 'IMAP      = ', IMAP
  WRITE(6,*) 'ICOV      = ', ICOV
  !WRITE(6,*) 'I_RV      = ', I_RV
  !WRITE(6,*) 'I_T       = ', I_T
  WRITE(6,*) 'I_RT     = ', I_RT
  WRITE(6,*) 'IAR       = ', IAR
  WRITE(6,*) 'I_VPVS    = ', I_VPVS
  WRITE(6,*) 'ISMPPRIOR = ', ISMPPRIOR
  WRITE(6,*) 'ISETSEED  = ', ISETSEED
  WRITE(6,*) 'IEXCHANGE = ', IEXCHANGE
  WRITE(6,*) 'I_VREF = ', I_VREF

  !WRITE(6,*) 'NTIME     = ', NTIME      !! No. time samples
  WRITE(6,*) 'NSRC      = ', NSRC       !! No. of rays (sources)
  WRITE(6,*) 'NLMN      = ', NLMN       !! Max number of layers
  WRITE(6,*) 'NLMX      = ', NLMX       !! Max number of layers
  WRITE(6,*) 'ICHAINTHIN= ', ICHAINTHIN !! Chain thinning interval
  WRITE(6,*) 'NKEEP     = ', NKEEP      !! No. models to keep before writing
  WRITE(6,*) 'NPTCHAINS1= ', NPTCHAINS1 !! Chain thinning interval
  WRITE(6,*) 'dTlog     = ', dTlog      !! Temperature increment
  WRITE(203,*) 'hmx       = ', hmx        !! Max. partition depth
  WRITE(203,*) 'hmin      = ', hmin       !! Min. layer thickness (must be small enough to not violate detailed balance)
  WRITE(6,*) 'armxH     = ', armxH      !! Max. AR prediction size
  WRITE(6,*) 'armxV     = ', armxV      !! Max. AR prediction size
  WRITE(6,*) 'TCHCKPT   = ', TCHCKPT    !! Checkpointing interval in s
  !WRITE(6,*) 'shift2    = ', shift2     !! time series shift for raysum
  !WRITE(6,*) 'width2    = ', width2     !! time series shift for raysum
  !WRITE(6,*) 'wter level= ', wl         !! water level for ray 3D
  !WRITE(6,*) 'sampling_dt= ', sampling_dt!! time series sampling rate raysum
  WRITE(6,*) 'NRF1       = ', 1      !! No. azimuth bins
  WRITE(6,*) 'Sample file: ',samplefile
  WRITE(6,*) ''
  WRITE(6,*) 'minlim:  '
  WRITE(6,203) minlim
  WRITE(6,*) 'maxlim:  '
  WRITE(6,203) maxlim
  WRITE(6,*) 'pertsdsc:'
  WRITE(6,203) pertsdsc
  WRITE(6,*) 'minlim sigma:  '
  WRITE(6,203) minlimsd
  WRITE(6,*) 'maxlim sigma:  '
  WRITE(6,203) maxlimsd
  WRITE(6,*) 'Done reading data.'
  WRITE(6,*) ''
  WRITE(6,*) ' ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
  WRITE(6,*) ''
ENDIF
CALL FLUSH(6)

RETURN
END SUBROUTINE READPARFILE
!==============================================================================

SUBROUTINE READDATA(obj)
!==============================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=RP):: iaz,idat,io
TYPE(objstruc)  :: obj

!IF(I_RV == 1)THEN
  !! Radial and Vertical components:
  !OPEN(UNIT=20,FILE=infileH,FORM='formatted',STATUS='OLD',ACTION='READ')
  !OPEN(UNIT=21,FILE=infileV,FORM='formatted',STATUS='OLD',ACTION='READ')
  !!
  !! Observed are NTIME long but need to zero pad to NTIME2 for convolution
  !!
  !obj%DobsR(1,1:NTIME2) = 0._RP
  !obj%DobsV(1,1:NTIME2) = 0._RP
  !DO iaz = 1,NRF1
  !  READ(20,*) obj%DobsR(iaz,1:NTIME)
  !  READ(21,*) obj%DobsV(iaz,1:NTIME)
  !ENDDO
  !CLOSE(20)
  !CLOSE(21)
 ! IF(I_T == 1)THEN
    !! Transverse components:
    !OPEN(UNIT=20,FILE=infileT,FORM='formatted',STATUS='OLD',ACTION='READ')
    !obj%DobsT(1,1:NTIME2) = 0._RP
    !DO iaz = 1,NRF1
    !  READ(20,*) obj%DobsT(iaz,1:NTIME)
    !ENDDO
    !CLOSE(20)
 ! ENDIF
!ELSEIF(I_RV == -1)THEN
  !! Radial and Vertical components:
 ! OPEN(UNIT=20,FILE=infileH,FORM='formatted',STATUS='OLD',ACTION='READ')
  !!
  !! Observed RF are NTIME long
  !!
  !obj%DobsR(1,:) = 0._RP
  !DO iaz = 1,NRF1
  !  READ(20,*) obj%DobsR(iaz,1:NTIME)
  !ENDDO
  !CLOSE(20)
!ENDIF
IF(I_RT == 1)THEN
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Ray traced travel time data: save in one column !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  OPEN(20,FILE=infileRT,FORM='formatted',STATUS='OLD',ACTION='READ')
  !ndat = 0
  DO idat=1,NDAT_RT
     READ(20,*,IOSTAT=io) obj%DobsRT(1,idat)
     !READ(20,*,IOSTAT=io) src_offset(idat), src_depth(idat)
     
     IF (io > 0) THEN
       STOP "Check input.  Something was wrong"
     ELSEIF (io < 0) THEN
       EXIT
     ELSE
     !  ndatad=ndatad+1
     ENDIF
     !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
  ENDDO
  CLOSE(20)! close the file
ENDIF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! RAY TRACER Source depths & offsets- read from the file: save in two columns, no delimeter please!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
OPEN(20,FILE=source_data_file,FORM='formatted',STATUS='OLD',ACTION='READ')
!ndat = 0
DO idat=1,NSRC
   !READ(20,*,IOSTAT=io) obj%periods(1,idat),obj%DobsRT(1,idat)
   !READ(20,*,IOSTAT=io) sour(1,idat)
   READ(20,*,IOSTAT=io) src_offset(idat), src_depth(idat)

   IF (io > 0) THEN
     STOP "Check input.  Something was wrong"
   ELSEIF (io < 0) THEN
     EXIT
   ELSE
   !  ndatad=ndatad+1
   ENDIF
   !if (i==ndatadmax) stop "number of Dispersion data >= ndatadmax"
ENDDO
CLOSE(20)! close the file













RETURN
END SUBROUTINE READDATA
!==============================================================================

SUBROUTINE PRINTPAR(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER :: ivo
TYPE(objstruc) :: obj

WRITE(6,*) 'Voronoi nodes:',obj%k
WRITE(6,*) 'NFP:',obj%NFP
DO ivo=1,obj%k
   WRITE(6,201) obj%voro(ivo,1:NPL)
   WRITE(6,205) obj%voroidx(ivo,1:NPL)
ENDDO
WRITE(6,*) 'Layer parameter vector:',obj%nunique,'layers'
DO ivo=1,obj%nunique
   WRITE(6,201) obj%par((ivo-1)*NPL+1:ivo*NPL)
ENDDO
WRITE(6,202) '            ',obj%par(obj%nunique*NPL+1:obj%nunique*NPL+(NPL-1))
WRITE(6,*) 'Partition:'
WRITE(6,203) obj%ziface(1:obj%nunique)
WRITE(6,*) 'Layers:'
WRITE(6,203) obj%hiface(1:obj%nunique)
IF(ICOV == 1)THEN
   WRITE(6,*) 'SD parameters:'
   !WRITE(6,206) 'sigma H   = ',obj%sdparR
   !WRITE(6,206) 'sigma V   = ',obj%sdparV
   !WRITE(6,206) 'sigma T   = ',obj%sdparT
   WRITE(6,206) 'sigma RT = ',obj%sdparRT
ENDIF
IF(IAR == 1)THEN
   WRITE(6,*) 'AR parameters:'
   WRITE(6,206) 'R and V: ',obj%arpar
   WRITE(6,206) 'RT:',obj%arparRT
ENDIF

201 FORMAT(6F12.4)
205 FORMAT(6I12)
202 FORMAT(A12,8F12.4)
203 FORMAT(11F12.4)
204 FORMAT(4F16.4)
206 FORMAT(a,128F12.4)
END SUBROUTINE PRINTPAR
!!=======================================================================
RECURSIVE FUNCTION LOGFACTORIAL(n)  RESULT(fact)
!-----Factorial------------------------------------------------------
!!=======================================================================

USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: fact
REAL(KIND=RP), INTENT(IN) :: n

IF (n == 0) THEN
   fact = 0
ELSE
   fact = LOG(n) + LOGFACTORIAL(n-1)
END IF

END FUNCTION LOGFACTORIAL
!==============================================================================
!EOF

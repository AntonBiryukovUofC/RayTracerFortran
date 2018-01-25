!!=======================================================================

SUBROUTINE LOGLHOOD(obj,ipred)
!!=======================================================================
!! This is LOGLHOOD function that combines various data sets for 
!! joint inversion. 
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB):: ipred
TYPE (objstruc) :: obj
REAL(KIND=RP)   :: logL_RT

!! log likelihood RT data
IF(I_RT == 1)THEN
  CALL LOGLHOOD_RT(obj,ipred,logL_RT)
ELSE
  logL_RT = 0._RP
ENDIF

!!
!! Joint likelihood (assumes independent errors on the vasious data sets)
!!
obj%logL = logL_RT
!IF(obj%logL == 0.)THEN
!  CALL PRINTPAR(obj)
!  STOP
!ENDIF

RETURN
END SUBROUTINE LOGLHOOD
!!=======================================================================

SUBROUTINE LOGLHOOD_RT(obj,ipred,logL)
!!=======================================================================
!!
!!  Compute predicted RT data and compute logL_RT
!!
USE RJMCMC_COM
USE MPI
USE RAYMOD
USE ieee_arithmetic
IMPLICIT NONE

INTEGER(KIND=IB)                      :: ipred,ierr_rt,imod,ibadlogL,ilay
TYPE (objstruc)                       :: obj
REAL(KIND=SP),DIMENSION(maxlay,10)    :: curmod
REAL(KIND=SP),DIMENSION(maxlay+NPREM,10):: curmod2
REAL(KIND=RP),DIMENSION(NDAT_RT)     :: DpredRT
REAL(KIND=RP),DIMENSION(NMODE)        :: EtmpRT
REAL(KIND=RP)                         :: logL,factvs,factvpvs
REAL(KIND=RP)                         :: tstart, tend, tcmp   ! Overall time 
REAL(KIND=RP),DIMENSION(obj%k+1)     :: vels
REAL(KIND=RP),DIMENSION(obj%k)     :: thickness

LOGICAL :: ISNAN

curmod = 0.
IF(IMAP == 1)THEN
  PRINT*,'CURMOD IN RT'
  DO ilay=1,obj%nunique+1
    WRITE(*,206)ilay,curmod(ilay,1:10)
  ENDDO
  206   FORMAT(I3,10F12.4)
ENDIF

!!
!!  curmod = thick  rho  alph  beta  %P  %S  tr  pl  st  dip 
!!

!! Find lowest valid entry for perturbations
!factvs   = obj%par(obj%nunique*NPL+1)  !! Vs is in km/s here
!factvpvs = obj%par(obj%nunique*NPL+2)
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,2) == 1)THEN
    factvs   = obj%voro(ilay,2)  !! Vs is in km/s here
    EXIT
  ENDIF
ENDDO
DO ilay=obj%k,1,-1
  IF(obj%voroidx(ilay,3) == 1)THEN
    factvpvs = obj%voro(ilay,3)
    EXIT
  ENDIF
ENDDO

curmod2  = 0.
curmod2(1:obj%nunique+1,:) = curmod(1:obj%nunique+1,:)
!! last layer thickness
!curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%par((obj%nunique-1)*NPL+1)*1000.
curmod2(obj%nunique+1,1)   = (vel_prem(1,1)*1000.)-obj%voro(obj%k,1)*1000.

!PRINT*,'obj voro'
!DO ilay=1,obj%k+1
!  PRINT*,obj%voro(ilay,:)
!ENDDO
!PRINT*,'obj par',obj%par(1:obj%NFP+2)
!PRINT*,'fact',factvs,factvpvs,obj%par((obj%nunique-1)*NPL+1)

!curmod2(obj%nunique+2:obj%nunique+NPREM,1)   = (vel_prem(1,2:NPREM)-vel_prem(1,1:NPREM-1)) * 1000. !! thickness in km
!curmod2(obj%nunique+1+NPREM,1)               = 0.                                                  !! HS thickness is 0 
!curmod2(obj%nunique+2:obj%nunique+1+NPREM,2) = vel_prem(4,1:NPREM) * 1000.            !! Density
!curmod2(obj%nunique+2:obj%nunique+1+NPREM,4) = (vel_prem(2,1:NPREM) + factvs) * 1000. !! Vs
!curmod2(obj%nunique+2:obj%nunique+1+NPREM,3) = (vel_prem(2,1:NPREM)+factvs)*(vel_prem(3,1:NPREM)+factvpvs)*1000. !! Vp

IF(IMAP == 1)THEN
  WRITE(*,*) 'curmod2 (including PREM)'
  DO ilay=1,obj%nunique+1+NPREM+1
    WRITE(*,206)ilay,curmod2(ilay,1:10)
  ENDDO
  !206   FORMAT(I3,10F12.4)
ENDIF

!periods = REAL(obj%periods(1,1:NDAT_RT),SP)
!!
!!  Need to append PREM perturbed by half-space perturbation here to 
!!  ensure that long period RT can be properly modelled. 
!!
!CALL dispersion(obj%nunique+1+NPREM,curmod2(1:obj%nunique+1+NPREM,2)/1000., & 
!     curmod2(1:obj%nunique+1+NPREM,3)/1000.,curmod2(1:obj%nunique+1+NPREM,4)/1000.,&
!     curmod2(1:obj%nunique+1+NPREM,1)/1000.,DpredRT,&
!     periods,NDAT_RT,ierr_rt)

vels = obj%voro(1:(obj%k+1),2) ! Retrieve P-wave velocities here (alphas)
thickness = obj%hiface(1:obj%k) ! Retrieve Layer thicknesses


CALL TraceRays(vels,thickness,obj%k,src_offset,src_depth,NSRC,DpredRT,1)


IF(ierr_rt /= 0)THEN
  logL = -HUGE(1._RP)
  RETURN
ENDIF

obj%DpredRT(1,1:NDAT_RT) = REAL(DpredRT,RP)
obj%DresRT(1,1:NDAT_RT) = obj%DobsRT(1,1:NDAT_RT)-obj%DpredRT(1,1:NDAT_RT)

ibadlogL = 0
IF(IAR == 1)THEN
  obj%DarRT  = 0._RP
  !!
  !!  Compute autoregressive model
  !!
  CALL ARPRED_RT(obj,1,1,NDAT_RT)
  !! Recompute predicted data as ith autoregressive model
  obj%DresRT = obj%DresRT-obj%DarRT

  !! Check if predicted AR model data are outside max allowed bounds
  CALL CHECKBOUNDS_ARMXRT(obj,ibadlogL)
ENDIF

!!
!!  Compute log likelihood
!!
IF(ibadlogL == 0)THEN
  EtmpRT = 0._RP
  IF(ICOV == 1)THEN
    !!
    !! Sample over sigma (one per mode)
    !!
    DO imod = 1,NMODE
      EtmpRT(imod) = LOG(1._RP/(2._RP*PI2)**(REAL(NDAT_RT,RP)/2._RP)) &
                    -(SUM(obj%DresRT(imod,:)**2._RP)/(2._RP*obj%sdparRT(imod)**2._RP)&
                    +REAL(NDAT_RT,RP)*LOG(obj%sdparRT(imod)))
    ENDDO
  ENDIF
  logL = SUM(EtmpRT)
  IF(ieee_is_nan(logL))THEN
    logL = -HUGE(1._RP)
  ENDIF
ELSE
   logL = -HUGE(1._RP)
   ibadlogL = 0
ENDIF

RETURN
207   FORMAT(500ES18.8)
END SUBROUTINE LOGLHOOD_RT
!!=======================================================================

SUBROUTINE INTERPLAYER_novar(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! This does not allow for variable layer complexity.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP)                    :: vref,vpvsref
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!! Sort node according to increasing depth:
obj%nunique = obj%k-1

tmpsort = 0._DRP
tmpsort = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort,RP)

obj%ziface = 0._RP
obj%ziface(1:obj%k-1) = obj%voro(2:obj%k,1)

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%k-1
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%k,1) = obj%ziface(1:obj%k)
!partmp(1:obj%k,2:NPL) = obj%voro(1:obj%k,2:NPL)

!!
!! Apply reference profile:
!!
IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,obj%voro(ivo,1))
    partmp(ivo,2) = vref + obj%voro(ivo,2)
    partmp(ivo,3) = vpvsref + obj%voro(ivo,3)
    !PRINT*,ivo,obj%voro(ivo,1),obj%voro(ivo,2),vref
  ENDDO
ENDIF

obj%par = 0._RP
DO ilay = 1,obj%k-1
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par((obj%k-1)*NPL+1:(obj%k*NPL)-1) = partmp(ilay,2:)

203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER_novar
!!=======================================================================

SUBROUTINE INTERPLAYER(obj)
!!=======================================================================
!!
!! This interpolates 1D layer nodes onto obj%par array for forward model
!! I.e., some layers are duplicates when nodes are not populated:
!!  Parameter 1:  Parameter 2:     Parameter 3:
!!   --o------      ---o---       ------o--------
!!     |               |                |
!!   ------o--      -------       --o------------
!!         |           |            |
!!         |           |            |
!!         |           |            |
!!   ---------      -------       ------------o--
!!         |           |                      |
!!   ---o-----      -------       ---------------
!!      |              |                      |
!!      |              |                      |
!!   ---------      -------       ---------------
!!
!! The layer node always defines the volume partition below the 
!! node position. The node position defines the interface.
!! The first layer is always fixed at 0 and all nodes are populated.
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,ivo2,ipar,ilay,iface,itmp,ntot
INTEGER(KIND=IB) :: NPL_tmp
TYPE(objstruc) :: obj
INTEGER(KIND=IB),DIMENSION(NPL-1)  :: niface
INTEGER(KIND=RP),DIMENSION(NPL*NLMX):: ifaceidx
REAL(KIND=RP),DIMENSION(NPL*NLMX,NPL):: partmp
REAL(KIND=RP),DIMENSION(NLMX,2)  :: vorotmp
REAL(KIND=RP),DIMENSION(NLMX,NPL-1)  :: voroh
REAL(KIND=RP),DIMENSION(NLMX+1,NPL-1,2):: voroz

REAL(KIND=RP),DIMENSION(NLMX,NPL)::  voro(NLMX,NPL),voroidx(NLMX,NPL)
REAL(KIND=RP)                    :: vref,vpvsref
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort
REAL(KIND=DRP),DIMENSION(NLMX*NPL):: tmpsort2

!  obj%k       = 3
!  obj%NFP     = (obj%k * NPL) + (NPL-1)
!  obj%voro    = 0._RP
!  obj%voroidx = 0
!  !! True parameters for simulation:
!  obj%voro(1,:) = (/  0.0_RP, 6.0_RP, 1.8_RP,  0.0_RP/)
!  obj%voro(2,:) = (/ 10.0_RP, 7.0_RP, 0.0_RP,  5.0_RP/)
!  obj%voro(3,:) = (/ 30.0_RP, 8.0_RP, 0.0_RP, 30.0_RP/)
!  obj%voroidx(1,:) = (/ 1, 1, 1, 0/)
!  obj%voroidx(2,:) = (/ 1, 1, 0, 1/)
!  obj%voroidx(3,:) = (/ 1, 1, 0, 1/)

!PRINT*,''
!PRINT*,''
!PRINT*,'Before:',rank
!DO ilay = 1,obj%k
!  WRITE(*,203)ilay,obj%voro(ilay,:)
!ENDDO
!PRINT*,''
!PRINT*,''

!obj%voroidx(1,:) = (/ 1, 1, 1/)
!obj%voroidx(2,:) = (/ 1, 1, 0/)
!obj%voroidx(3,:) = (/ 1, 1, 0/)

IF(IDIP == 1)THEN
  NPL_tmp = NPL-1
ELSE
  NPL_tmp = NPL
ENDIF

tmpsort = 0._DRP
tmpsort = REAL(obj%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),obj%voroidx(1:obj%k,:))
obj%voro(1:obj%k,:) = REAL(tmpsort,RP)

!!
!! Use local variable for voro and voroidx to allow easy change
!! from perturbation value to perturbation+vref
!!
voro    = 0._RP
voroidx = 0
voro    = obj%voro
voroidx = obj%voroidx

!!
!! Apply reference profile:
!!
IF(I_VREF == 1)THEN
  DO ivo = 1,obj%k
    CALL GETREF(vref,vpvsref,voro(ivo,1))
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3)
    IF(voroidx(ivo,2) == 1)voro(ivo,2) = vref + voro(ivo,2)
    IF(voroidx(ivo,3) == 1)voro(ivo,3) = vpvsref + voro(ivo,3)
    !WRITE(*,201)ivo,voro(ivo,1),voro(ivo,2),voro(ivo,3),vref,vpvsref
    !WRITE(*,*)''
  ENDDO
ENDIF
!201 FORMAT(I,5F12.6)
!!
!!  Find interfaces for each parameter
!!
niface = 0._RP
DO ipar = 2,NPL
  niface(ipar-1) = SUM(voroidx(:,ipar))-1
ENDDO
obj%ziface = 0._RP
voroz  = 0._RP
voroh  = 0._RP
iface  = 0
DO ipar = 2,NPL
  vorotmp = 0._RP
  itmp = 0
  DO ivo = 1,obj%k
    IF(voroidx(ivo,ipar) == 1)THEN
      itmp = itmp + 1
      vorotmp(itmp,:) = (/voro(ivo,1),voro(ivo,ipar)/)
    ENDIF
  ENDDO
  IF(niface(ipar-1) > 0)THEN
    DO ivo = 1,niface(ipar-1)
      iface = iface + 1
      obj%ziface(iface) = vorotmp(ivo+1,1)
      voroz(ivo,ipar-1,1) = obj%ziface(iface)
      voroz(ivo,ipar-1,2) = vorotmp(ivo,2)
    ENDDO
    voroz(itmp,ipar-1,2) = vorotmp(itmp,2)
  ELSE
    voroz(1,ipar-1,2) = vorotmp(1,2)
  ENDIF
ENDDO

ntot = SUM(niface)
tmpsort2 = 0._DRP
tmpsort2 = REAL(obj%ziface(1:ntot),DRP)
CALL QSORTC1D(tmpsort2(1:ntot))
obj%ziface(1:ntot) = REAL(tmpsort2,RP)

!!
!!  Find unique interfaces
!!
obj%nunique = obj%k-1
IF(ntot > 1)THEN
  DO ivo = 1,ntot
    itmp = 0
    ifaceidx = 0
    DO ivo2 = ivo+1,ntot
      IF(obj%ziface(ivo) /= obj%ziface(ivo2))THEN
        itmp = itmp + 1
        ifaceidx(itmp) = ivo2
      ENDIF
    ENDDO
    IF(itmp > 0)THEN
      obj%ziface(ivo+1:ivo+itmp) = obj%ziface(ifaceidx(1:itmp))
      obj%ziface(ivo+itmp+1:) = 0._RP
    ENDIF
  ENDDO
ENDIF

obj%hiface = 0._RP
obj%hiface(1) = obj%ziface(1)
DO ivo = 2,obj%nunique
  obj%hiface(ivo) = obj%ziface(ivo)-obj%ziface(ivo-1)
ENDDO

partmp = 0._RP
partmp(1:obj%nunique,1) = obj%ziface(1:obj%nunique)
DO ipar = 2,NPL_tmp
  ivo = 1
  IF(niface(ipar-1) /= 0)THEN
    DO ilay = 1,obj%nunique
      IF(obj%ziface(ilay) <= voroz(ivo,ipar-1,1))THEN
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ELSE
        ivo = ivo + 1
        partmp(ilay,ipar) = voroz(ivo,ipar-1,2)
      ENDIF
      IF(obj%ziface(ilay) >= voroz(niface(ipar-1),ipar-1,1))THEN
        partmp(ilay+1:obj%nunique+1,ipar) = voroz(ivo+1,ipar-1,2)
        EXIT
      ENDIF
    ENDDO
  ELSE
    partmp(1:obj%nunique+1,ipar) = voroz(1,ipar-1,2)
  ENDIF
ENDDO
partmp(1:obj%nunique,1) = obj%hiface(1:obj%nunique)
IF(IDIP == 1)THEN
  ipar = NPL
  ivo = 1
  partmp(:,ipar) = 0._RP
  IF(SUM(voroidx(:,ipar)) > 0)THEN
    DO ilay = 2,obj%nunique+1
      IF(voroidx(ilay,ipar) == 1)THEN
        partmp(ilay,ipar) = voro(ilay,ipar)
      ENDIF
    ENDDO
  ENDIF
ENDIF

obj%par = 0._RP
DO ilay = 1,obj%nunique
  obj%par((ilay-1)*NPL+1:ilay*NPL) = partmp(ilay,:)
ENDDO
obj%par(obj%nunique*NPL+1:((obj%nunique+1)*NPL)-1) = partmp(ilay,2:)

!IF(obj%voro(1,1)<0._RP)THEN
!PRINT*,''
!PRINT*,'After:',rank
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,partmp(ilay,:)
!ENDDO
!CALL PRINTPAR(obj)
!PRINT*,''
!PRINT*,obj%par
!PRINT*,''
!STOP
!PRINT*,'ntot',ntot
!PRINT*,'unique',obj%nunique
!PRINT*,'niface',niface
!PRINT*,'ziface',obj%ziface(1:ntot)
!DO ilay = 1,obj%nunique+1
!  WRITE(*,203)ilay,obj%par((ilay-1)*NPL+1:ilay*NPL)
!ENDDO
!PRINT*,''
!PRINT*,''
!  PRINT*,'INTERP'
!  STOP
!ENDIF
!
!STOP
!200 FORMAT(2I4,20F8.2)
!201 FORMAT(a,20F8.2)
!202 FORMAT(a,20I4)
203 FORMAT(I4,20F8.2)

END SUBROUTINE INTERPLAYER
!=======================================================================

SUBROUTINE GETREF(vref,vpvsref,z)
!!==============================================================================
!!
!! Reads observed data.
!!
USE RJMCMC_COM
IMPLICIT NONE
REAL(KIND=RP)   :: z,vref,vpvsref,grad,dz
INTEGER(KIND=IB):: ipar,iint
IF(z >= vel_ref(1,NVELREF))THEN
  vref = vel_ref(2,NVELREF)
  vpvsref = vel_ref(3,NVELREF)
ELSEIF(z == 0.)THEN
  vref = vel_ref(2,1)
  vpvsref = vel_ref(3,1)
ELSE
  iint = 0
  DO ipar=1,NVELREF
    IF((z-vel_ref(1,ipar)) <= 0.)THEN
      EXIT
    ENDIF
    iint = iint + 1
  ENDDO
  grad = (vel_ref(2,iint+1)-vel_ref(2,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  dz = (z-vel_ref(1,iint))
  vref = vel_ref(2,iint) + dz*grad
  grad = (vel_ref(3,iint+1)-vel_ref(3,iint))/(vel_ref(1,iint+1)-vel_ref(1,iint))
  vpvsref = vel_ref(3,iint) + dz*grad
ENDIF
RETURN
END SUBROUTINE GETREF
!!=======================================================================

SUBROUTINE ARPRED_RF(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxar(ifr) == 1)THEN
   k = 1
   !obj%DarR(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   !dres1 = obj%DresR(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arpar((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arpar((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   !obj%DarR(ifr,idata:idatb) = dar1
   !obj%DarR(ifr,idata) = 0._RP
   !obj%DarR(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_RF
!!=======================================================================

SUBROUTINE ARPRED_RT(obj,ifr,idata,idatb)
!=======================================================================
!!
!! Autoregressive model to model data error correlations.
!! AR process is computed forward and backward and average is used.
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj
INTEGER          :: i,j,k,ifr,idata,idatb
REAL(KIND=RP),DIMENSION(idatb-idata+1)::dres1,dar1
IF(obj%idxarRT(ifr) == 1)THEN
   k = 1
   obj%DarRT(ifr,idata)=0._RP          ! Matlab sets first point to zero...

   !!
   !! Real part:
   !!
   dres1 = 0._RP
   dres1 = obj%DresRT(ifr,idata:idatb)

   dar1(1)=0._RP          ! Matlab sets first point to zero...
   DO i=2,idatb-idata+1
      dar1(i) = 0
      IF(k >= i)THEN
         DO j=1,i-1
            dar1(i) = dar1(i) + obj%arparRT((ifr-1)+j) * dres1(i-j)
         ENDDO
      ELSE
         DO j=1,k
            dar1(i) = dar1(i) + obj%arparRT((ifr-1)+j) * dres1(i-j)
         ENDDO
      ENDIF
   ENDDO
   obj%DarRT(ifr,idata:idatb) = dar1
   obj%DarRT(ifr,idata) = 0._RP
   obj%DarRT(ifr,idatb) = 0._RP
ENDIF
END SUBROUTINE ARPRED_RT
!!=======================================================================

!SUBROUTINE CHECKBOUNDS_ARMXRF(obj,ibadlogL)
!!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: ifr,ibadlogL
!TYPE(objstruc):: obj

!DO ifr = 1,NRF1
!   IF(MAXVAL(obj%DarR(ifr,:)) > armxH)THEN
!      iarfail = iarfail + 1
!      ibadlogL = 1
!   ENDIF
!   IF(MINVAL(obj%DarR(ifr,:)) < -armxH)THEN
!      iarfail = iarfail + 1
!      ibadlogL = 1
!   ENDIF
!ENDDO

!RETURN
!END SUBROUTINE CHECKBOUNDS_ARMXRF
!!=======================================================================

SUBROUTINE CHECKBOUNDS_ARMXRT(obj,ibadlogL)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ifr,ibadlogL
TYPE(objstruc):: obj

DO ifr = 1,NMODE
   IF(MAXVAL(obj%DarRT(ifr,:)) > armxRT)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
   IF(MINVAL(obj%DarRT(ifr,:)) < -armxRT)THEN
      iarfail = iarfail + 1
      ibadlogL = 1
   ENDIF
ENDDO

RETURN
END SUBROUTINE CHECKBOUNDS_ARMXRT
!=======================================================================

SUBROUTINE LOGLHOOD2(obj)
!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc)  :: obj

!!
!!  Compute log likelihood
!!
obj%logL = 1._RP

RETURN
END SUBROUTINE LOGLHOOD2
!!=======================================================================
!!EOF

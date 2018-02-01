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
!  Agostinetti Malinverno 2010, Dettmer etal 2010, 2012, 2013
!
!==============================================================================

PROGRAM  RJMCMC_PLANE

!=======================================================================
USE RJMCMC_COM
USE NR
USE ieee_arithmetic
IMPLICIT NONE
INTRINSIC RANDOM_NUMBER, RANDOM_SEED

INTEGER(KIND=IB)  :: i,j,ipar,ipar2,ifr,ilay,ifreq,imcmc,ithin,isource,ikeep,ik,iaz
INTEGER(KIND=IB)  :: idat,iang,it,ic,it2,ivo,isource1,isource2,io

TYPE (objstruc),DIMENSION(2)            :: objm     ! Objects on master node
TYPE (objstruc)                         :: obj      ! Object on slave node
TYPE (objstruc)                         :: objnew   ! Object on slave node
!TYPE (covstruc),ALLOCATABLE,DIMENSION(:):: cov      ! Structure for inverse data Covariance Matrices
REAL(KIND=RP)                           :: ran_uni  ! Uniform random number
REAL(KIND=RP)                           :: ran_nor  ! Normal random number
REAL(KIND=RP),DIMENSION(:),ALLOCATABLE  :: tmpvoro


!!---------------------------------------------------------------------!
!!     MPI stuff:
!!---------------------------------------------------------------------!
!!
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds

REAL(KIND=RP)               :: tstart, tend              ! Overall time 
REAL(KIND=RP)               :: tstart2, tend2            ! Time for one forward model computation
REAL(KIND=RP)               :: tstartsnd, tendsnd        ! Communication time
REAL(KIND=RP)               :: tstartcmp, tendcmp, tcmp  ! Forward computation time

CALL MPI_INIT( ierr )
CALL MPI_COMM_RANK( MPI_COMM_WORLD, rank, ierr )
CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NTHREAD, ierr )
OPEN(UNIT=20,FILE=filebasefile,STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)

IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~             Reversible Jump MCMC Sampling              ~~~  '
IF(rank == src)WRITE(6,*) '~~~                                                        ~~~  '
IF(rank == src)WRITE(6,*) '~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  '
IF(rank == src)WRITE(6,*) '...running on ',NTHREAD,' cores'

!!
!! Read parameter file and allocate global prior limits:
!!
CALL READPARFILE()
!! Allocate objects for sampling
CALL ALLOC_OBJ(objm(1))
CALL ALLOC_OBJ(objm(2))
CALL ALLOC_OBJ(obj)
CALL ALLOC_OBJ(objnew)
IF(rank == src)THEN
  !!
  !!  File to save posterior samples
  !!
  OPEN(NEWUNIT(usample),FILE=samplefile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE')
  OPEN(NEWUNIT(ustep),FILE=stepsizefile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE')
ENDIF
ncount1 = NFPMX+NAP+3*1+3*1+2*NMODE
ncount2 = NFPMX+1+3*1+3*1+2*NMODE
ALLOCATE( sample(NKEEP,ncount1),tmpvoro(NFPMX) )
ALLOCATE( tmpmap(ncount2) )
sample       = 0._RP

!!print *, "blya"

!!  Make MPI structure to match objstruc
!!
! Specify NCHAIN number here:




CALL MAKE_MPI_STRUC_SP(objm(1),objtype1)
CALL MAKE_MPI_STRUC_SP(objm(2),objtype2)
CALL MAKE_MPI_STRUC_SP(obj,objtype3)

!!------------------------------------------------------------------------
!!  Read in data
!!------------------------------------------------------------------------
CALL READDATA(obj)

!!------------------------------------------------------------------------
!!
!! Population / parallel tempering:
!!
!!------------------------------------------------------------------------
IF(NTHREAD > 1)THEN
  NPTCHAINS  = NTHREAD    ! Number of chains (equal to number of slaves)
  NT = NPTCHAINS-NPTCHAINS1+1
  ALLOCATE( NCHAINT(NT),beta_pt(NPTCHAINS) )
  NCHAINT = 1
  NCHAINT(1) = NPTCHAINS1

!!------------------------------------------------------------------------
!!  Set up tempering schedule
!!------------------------------------------------------------------------
  IF(rank == src)WRITE(6,*)'NT',NT
  IF(rank == src)WRITE(6,*)'NCHAINT',NCHAINT
  it2 = 1_IB
  DO it=1,NT
    DO ic=1,NCHAINT(it)
      beta_pt(it2) = 1._RP/dTlog**REAL(it-1_IB,RP)
      IF(rank == src)WRITE(6,218) 'Chain ',it2,':  beta = ',beta_pt(it2),'  T = ',1._RP/beta_pt(it2)
      it2 = it2 + 1
    ENDDO
  ENDDO
ENDIF
IF(rank == src)WRITE(6,*) ''
218 FORMAT(a6,i4,a10,F8.6,a6,F8.2)

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!! Initialize random seeds on each core (Call RANDOM_SEED only once in the whole code. PARALLEL_SEED calls it)
!!
!CALL RANDOM_SEED


CALL PARALLEL_SEED()

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!! Read mapfile to start:
!!
OPEN(UNIT=20,FILE=mapfile,FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) tmpmap
CLOSE(20)

obj%k       = INT(tmpmap(1),IB)   !! No. interfaces + 1 (e.g. 3 makes 2 layers and a half space)
obj%NFP     = (obj%k * NPL)       !! No. forward parameters; NPL -> max No. parameters per layer
obj%voro    = 0._RP               !! Main array of layer nodes
obj%voroidx = 1
obj%arpar   = 0._RP
obj%sdparRT= 0._RP
obj%arparRT= 0._RP

obj%sdaveRT= 0._RP

DO ivo = 1,obj%k
  ipar = (ivo-1)*NPL+2
  obj%voro(ivo,:) = tmpmap(ipar:ipar+NPL-1)

  IF (I_VARPAR == 0) THEN
  obj%voroidx(ivo,:) = 1
  ELSE
      obj%voroidx(ivo,:) = 1
      DO ipar2 = 1,NPL
        IF(obj%voro(ivo,ipar2) < -99._RP)THEN
          obj%voroidx(ivo,ipar2) = 0
        ELSE
          obj%voroidx(ivo,ipar2) = 1
        ENDIF
      ENDDO
  ENDIF
ENDDO
! Temporary change to voroidx!!!!!!!!
obj%voroidx = 1

!IF(ICOV == 1) obj%sdparR   = tmpmap(NFPMX+1+1:NFPMX+1+NRF1)
!IF(ICOV == 1) obj%sdparV   = tmpmap(NFPMX+1+1+NRF1:NFPMX+1+2*NRF1)
!IF(ICOV == 1) obj%sdparT   = tmpmap(NFPMX+1+1+2*NRF1:NFPMX+1+3*NRF1)
IF(ICOV == 1) obj%sdparRT = tmpmap(NFPMX+1+1:NFPMX+1+NMODE)

IF(IAR == 1)THEN ! These indices are totally wrong!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  obj%arpar = tmpmap(NFPMX+1+3*1+NMODE+1:NFPMX+3*1+NMODE+3*1+1)
  obj%arparRT = tmpmap(NFPMX+1+3*1+NMODE+3*1+1:NFPMX+3*1+NMODE+3*1+NMODE+1)
  obj%idxar = 1
  DO ipar = 1,1
    IF(obj%arpar(ipar) < minlimar(ipar)) obj%idxar(ipar) = 0
  ENDDO
  obj%idxarRT = 1
  DO ipar = 1,NMODE
    IF(obj%arparRT(ipar) < minlimarRT(ipar)) obj%idxarRT(ipar) = 0
  ENDDO
ENDIF


IF(I_VARPAR == 0)CALL INTERPLAYER_novar(obj)
IF(I_VARPAR == 1)CALL INTERPLAYER(obj)
IF(rank /= src)obj%beta = beta_pt(rank)
IF(rank ==src) CALL PRINTPAR(obj)

ALLOCATE( icount(NTHREAD) )
icount = 0
tstart = MPI_WTIME()

CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)

!CALL INTERPLAYER(obj)
!CALL LOGLHOOD(obj,1)
!CALL PRINTPAR(obj)
!WRITE(6,*) 'logL = ',obj%logL
!IF(IMAP == 1)THEN
!  CALL CHECKBOUNDS(obj)
!  tstart2 = MPI_WTIME()
!  DO ic=1,1000
!    !CALL PROPOSAL(obj,objnew,2,4,1._RP)
!    obj%voro(2,4) = 0. + REAL(ic,RP)*30./1000.
!    CALL INTERPLAYER(obj)
!    CALL LOGLHOOD(obj,1)
!!    CALL PRINTPAR(obj)
!    PRINT*,'in loop',ic,obj%logL,obj%voro(2,4)
!    WRITE(77,*) obj%logL,obj%voro(2,4)
!!    PRINT*,''
!!    PRINT*,''
!!    PRINT*,''
!!    PRINT*,''
!  ENDDO
!  tend2 = MPI_WTIME()
!  CALL PRINTPAR(obj)
!  WRITE(6,*) 'time = ',tend2-tstart2
!  WRITE(6,*) 'logL = ',obj%logL
!  !CALL SAVEREPLICA(obj,predfile,obsfile,arfile,predfileRT,obsfileRT,arfileRT)
!  !IF(ioutside == 1)WRITE(*,*)'FAILED in starting model'
!  STOP
!ENDIF


IF(IMAP == 1)THEN
  IF(rank == src)THEN 
   WRITE(6,*) 'IMAP activated, exiting after computing replica for MAP.'
   CALL CHECKBOUNDS(obj)
   tstart2 = MPI_WTIME()
   !DO ic=1,100
   !CALL PRINTPAR(obj)
   !CALL PROPOSAL(obj,objnew,12,1,1.)
   !CALL PRINTPAR(objnew)
   !STOP
   IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
   IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)
   !ENDDO
   tend2 = MPI_WTIME()
   WRITE(6,*) 'time = ',tend2-tstart2
   WRITE(6,*) 'logL = ',obj%logL
   CALL SAVEREPLICA(obj,predfileRT,obsfileRT)
   IF(ioutside == 1)WRITE(*,*)'FAILED in starting model'
  ENDIF
  CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
  CALL MPI_FINALIZE( ierr )

  STOP
ELSE
  IF(rank == 1) THEN ! changed rank == 1 to rank == src
    WRITE(6,*) 'Starting model:'
    IF(rank ==src) CALL PRINTPAR(obj)
    CALL CHECKBOUNDS(obj)
    IF(ioutside == 1)WRITE(*,*)'FAILED in starting model'
    tstart2 = MPI_WTIME()
    IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
    IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)
    tend2 = MPI_WTIME()
    WRITE(6,*) 'logL = ',obj%logL
    WRITE(6,*) 'time = ',tend2-tstart2
    CALL CHECKBOUNDS(obj)
    IF(ioutside == 1)THEN 
      PRINT*,'Outside bounds.'
      STOP
    ENDIF
    tend2 = MPI_WTIME()
    WRITE(6,*) 'time for linrot = ',tend2-tstart2
  ENDIF
ENDIF


CALL FLUSH(6)
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
!!
!!  Need this call to initiate all chains!
!!
IF(ISMPPRIOR == 0)CALL LOGLHOOD(obj,1)
IF(ISMPPRIOR == 1)CALL LOGLHOOD2(obj)

icount(rank+1) = icount(rank+1) + 1

!------------------------------------------------------------------------
!
!          ************ RJMCMC Sampling ************
!
! -----------------------------------------------------------------------
IF(rank == src)WRITE(6,*) 'Starting RJMCMC sampling...'

!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!     MASTER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
IF(rank==src)THEN
ncswap = 0
ncswapprop = 0
ikeep = 1
tsave1 = MPI_WTIME()
DO imcmc = 1,NCHAIN
  print *, 'imcmc1'
  !!
  !!  Receiving samples from any two slave (no particular order -> auto load
  !!  balancing) and propose swap:
  !!
  DO ithin = 1,ICHAINTHIN
    print *, 'imcmc2'

    CALL MPI_RECV(objm(1), 1, objtype1, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
    isource1 = status(MPI_SOURCE)  !! This saves the slave id for the following communication
    tstart = MPI_WTIME()
    print *, 'imcmc3'

    CALL MPI_RECV(objm(2), 1, objtype2, MPI_ANY_SOURCE,MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr )
    isource2 = status(MPI_SOURCE)  !! This saves the slave id for the following communication
    print *, 'imcmc4'
    print *, 'For objm1'
    print *, objm(1)%voroidx

    print *, 'For objm2'
    print *, objm(2)%voroidx
    
    IF(IEXCHANGE == 1) THEN
      CALL TEMPSWP_MH(objm(1),objm(2))
    !IF(IAR == 1)CALL UDATE_SDAVE(objm(1),objm(2))
    print *, 'imcmc5'
    ENDIF
    CALL MPI_SEND(objm(1), 1,objtype1, isource1, rank, MPI_COMM_WORLD, ierr)
    CALL MPI_SEND(objm(2), 1,objtype2, isource2, rank, MPI_COMM_WORLD, ierr)
    tend = MPI_WTIME()
    print *, 'imcmc6'

  ENDDO

!  t_chckpt2 = MPI_WTIME()
!  !! Checkpoint every ... seconds:
!  IF(NINT(t_chckpt2-t_chckpt1)>TCHCKPT)THEN
!    !PRINT*,rank,'dumps object at',t_chckpt2-t_chckpt1,'s'
!    !PRINT*,NINT(t_chckpt2-t_chckpt1),TCHCKPT
!    CALL SAVEOBJECT(obj)
!    t_chckpt1 = MPI_WTIME()
!  ENDIF

  CALL SAVESAMPLE(objm,ikeep,isource1,isource2,REAL(tend-tstart,RP))

!   tstartsnd = MPI_WTIME()
!   CALL SAVESAMPLE(logLG,isource,tcmp)
!   tendsnd = MPI_WTIME()
!   IF(MOD(imcmc,30+1)==0)THEN
!      WRITE(*,*) ''
!      WRITE(6,203) '   imcmc,           logL,       % failed,      k,       iacc/irej, iacc_bd, iacc_bds,  time(send), source,  time(comp)'
!      WRITE(6,203) '----------------------------------------------------------------------------------------------------------------------'
!   ENDIF
!   WRITE(6,202) imcmc,sample(1,(/1,2/)),INT(sample(1,4)),sample(1,5+NFPMX+NRF1+(NARFP*NRF1)),&
!                INT(sample(1,6+NFPMX+NRF1+(NARFP*NRF1))),INT(sample(1,6+NFPMX+NRF1+(NARFP*NRF1)+2)),tendsnd-tstartsnd,isource,tcmp
!!     sample(ikeep,:) =  (/ obj%logL, obj%logPr, REAL(ishearfail,RP)/REAL(i_ref_nlay,RP), REAL(obj%k,RP), & ! 4 parameters
!!                           tmpvoro,obj%sdpar,obj%arpar, & 
!!                           REAL(iaccept,RP)/REAL(iaccept+ireject,RP),REAL(obj%iaccept_bd,RP),&
!!                           REAL(obj%ireject_bd,RP),REAL(i_bd,RP),REAL(ic,RP),REAL(rank,RP) /)
!   CALL FLUSH(6)
!
ENDDO
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!
!!    WORKER PART
!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
!!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
ELSE

tstart = MPI_WTIME()

DO imcmc = 1,NCHAIN

  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!
  !!                      EXPLORE CALLS
 ! print *,'slave1'
  DO ithin = 1,1
    !! BD MH for fixed complexity:
    IF (I_VARPAR == 0)THEN 
      CALL EXPLORE_MH_NOVARPAR(obj,objnew,obj%beta)
    ENDIF
    !! BD MH for variable complexity:
    IF (I_VARPAR == 1)THEN
      print *, ' Got IVARPAR accident'
      STOP
      CALL EXPLORE_MH_VARPAR(obj,objnew,obj%beta)

    ENDIF
!print *,'slave2'
  ENDDO
  !! MH for non partition parameters:
  CALL EXPLORE_MH(obj,objnew,obj%beta)
!print *,'slave3'
  !!                    END EXPLORE CALLS                                    !!
  !!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!!

  tend = MPI_WTIME()
  obj%tcmp = REAL(tend-tstart,RP)
  !print *, 'on slave:before'
  !IF (rank ==1) print *,obj%voroidx
  CALL MPI_SEND(obj, 1,objtype3, src, rank, MPI_COMM_WORLD, ierr)

  CALL MPI_RECV(obj, 1,objtype3, src, MPI_ANY_TAG, MPI_COMM_WORLD, status, ierr)
  tstart = MPI_WTIME()
  !print *, 'on slave:after'
  !IF (rank ==1) print *,obj%voroidx

!  t_chckpt2 = MPI_WTIME()
  !! Checkpoint every once in a while...
!  IF(NINT(t_chckpt2-t_chckpt1)>TCHCKPT)THEN
!    PRINT*,rank,'dumps object at',obj%logL
!    !PRINT*,NINT(t_chckpt2-t_chckpt1),TCHCKPT
!    CALL SAVEOBJECT(obj)
!    t_chckpt1 = MPI_WTIME()
!  ENDIF

  !! Scale proposals (diminishing adaptation after Handbook of MCMC)
!  obj%iacceptvoro(:,1:NCMT)  = SUM(iacceptcmt_buf,3)
!  obj%iproposevoro(:,1:NCMT) = SUM(iproposecmt_buf,3)
!  obj%iacceptvoro(:,NCMT+1:NPL)  = SUM(iacceptloc_buf,3)
!  obj%iproposevoro(:,NCMT+1:NPL) = SUM(iproposeloc_buf,3)
!  !IF(rank == 1)PRINT*, 'propose:',obj%iproposevoro
!  !IF(rank == 1)PRINT*, 'accept:',obj%iacceptvoro
!  IF(IADAPT == 1)THEN
!    !IF(MOD(imcmc,NBUF) == 0)THEN
!    IF(iadaptcmt == 1)THEN
!      iadaptcmt = 0
!      DO ivo = 1,obj%k
!      DO ipar = 1,NCMT
!        IF(REAL(SUM(iacceptcmt_buf(ivo,ipar,:)),RP)/REAL(SUM(iproposecmt_buf(ivo,ipar,:)),RP) < 0.20_RP)THEN
!          obj%pertsd(ivo,ipar) = obj%pertsd(ivo,ipar)*MAX(0.90_RP,1._RP-1._RP/SQRT(REAL(imcmc2,RP)))
!        ENDIF
!        IF(REAL(SUM(iacceptcmt_buf(ivo,ipar,:)),RP)/REAL(SUM(iproposecmt_buf(ivo,ipar,:)),RP) > 0.30_RP)THEN
!          obj%pertsd(ivo,ipar) = obj%pertsd(ivo,ipar)/MAX(0.90_RP,1._RP-1._RP/SQRT(REAL(imcmc2,RP)))
!        ENDIF
!      ENDDO
!      ENDDO
!    ENDIF
!  ENDIF
  !! MCMC counter for diminishing adaptation
!  imcmc2 = imcmc2 + 1_IB
ENDDO
ENDIF !! MPI ENDIF

!CLOSE(6)
CALL MPI_FINALIZE( ierr )

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

END PROGRAM RJMCMC_PLANE
!!==============================================================================

SUBROUTINE EXPLORE_MH(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for all hierarchical and non partition parameters (no birth or 
!! death)
!!
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
INTEGER(KIND=IB)                 :: ncra,iparst
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio,logy,logq1_1,logq1_2
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!!
!! Do Metropolis-Hastings on data-error standard deviations
!!
IF(1 == 1)THEN
IF(ICOV == 1)THEN
  !!
  !! Sample standard deviation of RT
  !!
  IF(I_RT == 1)THEN
  CALL RANDOM_NUMBER(ran_uni_ar)
  IF(ran_uni_ar>=0.10_RP)THEN
    DO ipar = 1,NMODE
      CALL PROPOSAL_SDRT(obj,objnew1,ipar)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
    ENDDO
  ENDIF
  ENDIF
ENDIF ! ICOV if
!!
!! Do Metropolis-Hastings on autoregressive model
!!
!!
!! Do Metropolis-Hastings on autoregressive model RT
!!
IF(IAR == 1)THEN
IF(I_RT == 1)THEN
  !! Perturb AR model with .25 probability
  CALL RANDOM_NUMBER(ran_uni_ar)
  !IF(ran_uni_ar>=0.25_RP)THEN
    DO ipar = 1,NMODE
      IF(obj%idxarRT(ipar) == 0)THEN
        !! Propose birth
        arptype = 1
        logarp = LOG(0.5_RP)
      ELSE
        CALL RANDOM_NUMBER(ran_uni_ar)
        IF(ran_uni_ar>=0.5_RP)THEN
          !! Propose death
          arptype = 2
          logarp = LOG(2._RP)
        ELSE
          !! Propose perturb
          arptype = 3
          logarp = 0._RP
        ENDIF
      ENDIF
      CALL PROPOSAL_ARRT(obj,objnew1,ipar,arptype)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,0)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !!logPLratio = (objnew1%logL - obj%logL)*beta_mh
        !! Input Birth Death AR here:
        logPLratio = logarp + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF
      ELSE
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF
      i_sdpert = 0
    ENDDO
  !ENDIF ! ran_uni_ar if
ENDIF ! I_RT if
ENDIF ! AR if

ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH
!!==============================================================================

SUBROUTINE EXPLORE_MH_NOVARPAR(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for partition parameters with Birth and Death for the case 
!! of fixed complexity in layers
!!
USE DATA_TYPE
USE RJMCMC_COM
USE MPI
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
INTEGER(KIND=IB)                 :: ncra,iparst
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio,logPratio
REAL(KIND=RP)                    :: logy,logq1_1,logq1_2
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: znew,beta_mh,Lr1,Lr2,PROP
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)
i_bd = 0
!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(kmin /= kmax)THEN
  STOP
  !! Perturbing k:
  CALL RANDOM_NUMBER(ran_unik)
  IF(obj%k == kmax)THEN  
    !! If k == kmax, no birth allowed, so 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 2
  ELSEIF(obj%k == kmin)THEN  
    !! If k == kmin, no death allowed, so 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 1
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF((ran_unik <= 0.3333_RP))i_bd = 1
    IF(ran_unik > 0.6666_RP)i_bd = 2
  ENDIF
  !! Proposal for Birth or Death:
  IF(i_bd == 1) CALL BIRTH_FULL(obj,objnew1)
  IF(i_bd == 2) CALL DEATH_FULL(obj,objnew1)
  IF(obj%k /= objnew1%k)THEN
    !!
    !! If k changed, check BD acceptance
    !!
    CALL CHECKBOUNDS(objnew1)
    IF(ioutside == 0)THEN
      IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
      IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
      !logPLratio =  (objnew1%logL - obj%logL)*beta_mh
      logPratio  = objnew1%logPr
      logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh

      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
        objnew1 = obj
        obj%ireject_bd = obj%ireject_bd + 1
      ELSE
        obj = objnew1
        obj%iaccept_bd = obj%iaccept_bd + 1
      ENDIF
    ELSE
      objnew1 = obj
      obj%ireject_bd = obj%ireject_bd + 1
      ioutside = 0
    ENDIF
  ENDIF  ! k-change if
ENDIF
!!
!! Carry out MH sweep
!!
!!
!! Do Metropolis-Hastings update on c, rho, alpha
!!
IF(1 == 1)THEN
DO ivo = 1,obj%k
  idxrand = 0
  idxrand(1:NPL) = RANDPERM(NPL)
  DO ipar = 1,NPL
    !iwhich = idxrand(ipar)
    iwhich = ipar
    IF(ivo == 1 .AND. iwhich == 1) CYCLE
    !print *, 'Before'
    !print *, obj%voroidx(ivo,iwhich)
    IF(obj%voroidx(ivo,iwhich) == 1)THEN


      CALL PROPOSAL(obj,objnew1,ivo,iwhich,1._RP)
      CALL CHECKBOUNDS2(objnew1,ivo,iwhich)
      !print *,obj%voroidx
      !print *, 'After'

      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,1)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        !logPLratio = (objnew1%logL - obj%logL)*beta_mh
        logPratio  = objnew1%logPr
        logPLratio = logPratio + (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF  !! MH if
      ELSE !! outside before delayed rejection
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF !! ioutside if
    ENDIF !!  voroidx if
  ENDDO
ENDDO
ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH_NOVARPAR
!!==============================================================================
SUBROUTINE EXPLORE_MH_VARPAR(obj,objnew1,beta_mh)
!!==============================================================================
!!
!! This is MH for partition parameters with Birth and Death for the case 
!! of variable complexity in layers (most general case with 2 different 
!! birth and death steps.
!!
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iwhich,idel,ipar,ilay,ivo
TYPE(objstruc)                   :: obj,objnew1,objnew2
!TYPE (covstruc),DIMENSION(NRF1)  :: cov
REAL(KIND=RP)                    :: logPLratio
REAL(KIND=RP)                    :: ran_uni,ran_uni_BD, ran_unik,ran_uni_ar
REAL(KIND=RP)                    :: beta_mh
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxrand
INTEGER(KIND=IB)                 :: ipick
INTEGER(KIND=IB),DIMENSION(NLMX) :: idxpick
REAL(KIND=RP)                    :: logarp
INTEGER(KIND=IB)                 :: arptype,choose

objnew1 = obj
!! Draw uniform Birth-Death probability
CALL RANDOM_NUMBER(ran_uni_BD)
i_bd = 0
!! Do BIRTH-DEATH MCMC with 0.5 probability
IF(kmin /= kmax)THEN
  !! Perturbing k:
  CALL RANDOM_NUMBER(ran_unik)
  IF(obj%k == kmax)THEN
    !! If k == kmax, no birth allowed, so 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 2
  ELSEIF(obj%k == kmin)THEN
    !! If k == kmin, no death allowed, so 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP)i_bd = 1
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF(ran_unik <= 0.3333_RP)i_bd = 1
    IF(ran_unik > 0.6666_RP)i_bd = 2
  ENDIF
  IF(i_bd == 1) CALL BIRTH_VARPAR(obj,objnew1)
  IF(i_bd == 2) CALL DEATH_VARPAR(obj,objnew1)
  IF(obj%k /= objnew1%k)THEN
    !!
    !! If k changed, check BD acceptance
    !!
    CALL CHECKBOUNDS(objnew1)
    IF(ioutside == 0)THEN
      IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
      IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
      logPLratio =  (objnew1%logL - obj%logL)*beta_mh

      CALL RANDOM_NUMBER(ran_uni)
      IF(ran_uni >= EXP(logPLratio))THEN
        objnew1 = obj
        obj%ireject_bd = obj%ireject_bd + 1
      ELSE
        obj = objnew1
        obj%iaccept_bd = obj%iaccept_bd + 1
      ENDIF
    ELSE
      objnew1 = obj
      obj%ireject_bd = obj%ireject_bd + 1
      ioutside = 0
    ENDIF
  ENDIF  ! k-change if
ENDIF
IF(IBD_SINGLE == 1)THEN
IF(obj%k > 1)THEN
  objnew1 = obj
  i_bds = 0
  !! Pick random node (except for first one):
  idxpick = 0
  idxpick(1:obj%k-1) = RANDPERM(obj%k-1)
  ipick = idxpick(1)+1

  CALL RANDOM_NUMBER(ran_unik)
  IF(SUM(obj%voroidx(ipick,:)) == 2)THEN
    !! If node has min parameters (2), 1/3 birth, 2/3 stay local
    IF(ran_unik<=0.3333_RP) i_bds = 1
  ELSEIF(SUM(obj%voroidx(ipick,:)) == NPL)THEN
    !! If node has max parameters, 1/3 death, 2/3 stay local
    IF(ran_unik<=0.3333_RP) i_bds = 2
  ELSE
    !! 1/3 birth, 1/3 death, 1/3 stay local
    IF(ran_unik <= 0.3333_RP) i_bds = 1
    IF(ran_unik >  0.6666_RP) i_bds = 2
  ENDIF
  IF(i_bds == 1) CALL BIRTH_SINGLE(obj,objnew1,ipick)
  IF(i_bds == 2) CALL DEATH_SINGLE(obj,objnew1,ipick)
  CALL CHECKBOUNDS(objnew1)
  IF(ioutside == 0)THEN
    IF(ISMPPRIOR == 0) CALL LOGLHOOD(objnew1,1)
    IF(ISMPPRIOR == 1) CALL LOGLHOOD2(objnew1)
    !! Likelihood: ratio
    logPLratio =  (objnew1%logL - obj%logL)*beta_mh

    CALL RANDOM_NUMBER(ran_uni)
    IF(ran_uni >= EXP(logPLratio))THEN
      objnew1 = obj
      obj%ireject_bds = obj%ireject_bds + 1
    ELSE
      obj = objnew1
      obj%iaccept_bds = obj%iaccept_bds + 1
    ENDIF
  ELSE
    objnew1 = obj
    obj%ireject_bds = obj%ireject_bds + 1
    ioutside = 0
  ENDIF
ENDIF ! k > 1 if
ENDIF ! BD single if
!!
!! Carry out MH sweep
!!
!!
!! Do Metropolis-Hastings update on c, rho, alpha
!!
IF(1 == 1)THEN
DO ivo = 1,obj%k
  idxrand = 0
  idxrand(1:NPL) = RANDPERM(NPL)
  DO ipar = 1,NPL
    iwhich = idxrand(ipar)
    !! Skip the topmost node position, since it's fixed at 0.
    IF(ivo == 1 .AND. iwhich == 1) CYCLE
    IF(obj%voroidx(ivo,iwhich) == 1)THEN
      CALL PROPOSAL(obj,objnew1,ivo,iwhich,1._RP)
      CALL CHECKBOUNDS2(objnew1,ivo,iwhich)
      IF(ioutside == 0)THEN
        IF(ISMPPRIOR == 0)CALL LOGLHOOD(objnew1,1)
        IF(ISMPPRIOR == 1)CALL LOGLHOOD2(objnew1)
        logPLratio = (objnew1%logL - obj%logL)*beta_mh
        CALL RANDOM_NUMBER(ran_uni)
        IF(ran_uni >= EXP(logPLratio))THEN
          objnew1 = obj
          ireject = ireject + 1
        ELSE
          obj = objnew1
          iaccept = iaccept + 1
        ENDIF  !! MH if
      ELSE !! outside before delayed rejection
        objnew1 = obj
        ireject = ireject + 1
        ioutside = 0
      ENDIF !! ioutside if
    ENDIF !!  voroidx if
  ENDDO
ENDDO
ENDIF !! 1 == 2 if
END SUBROUTINE EXPLORE_MH_VARPAR
!=======================================================================

SUBROUTINE DEATH_FULL(obj,objnew)
!!=======================================================================
!!
!! This DEATH is the reverse for Birth from prior.
!!
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel,ivo,ipar
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort
REAL(KIND=RP)                    :: zdel,zj,zjp1

objnew = obj
objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!! 1) Pick random node:
idxdeath = 0
idxdeath(1:obj%k-1) = RANDPERM(obj%k-1)
idxdeath = idxdeath + 1
idel = idxdeath(1)
zdel = objnew%voro(idel,1)
zj   = objnew%voro(idel-1,1)
IF(idel == obj%k)THEN
  zjp1 = hmx
ELSE
  zjp1 = objnew%voro(idel+1,1)
ENDIF

!! 2) Save location and parameters of node
!objnew%gvoroidx = 1

!! 3) delete node and re-interpolate Voronoi model
objnew%voro(idel,:) = 0._RP
objnew%voroidx(idel,:) = 0
tmpsort = 0._DRP
tmpsort = REAL(objnew%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),objnew%voroidx(1:obj%k,:))
objnew%voro(1:obj%k,:) = REAL(tmpsort,RP)
!CALL QSORTC2D(objnew%voro(1:obj%k,:),objnew%voroidx(1:obj%k,:))
!! 4) Quicksort stores deleted array (all zeros) first.
!!    Hence, move all one up unless deleted node is last node.
objnew%voro(1:objnew%k,:) = objnew%voro(2:obj%k,:)
objnew%voroidx(1:objnew%k,:) = objnew%voroidx(2:obj%k,:)
objnew%voro(objnew%k+1,:) = 0._RP
objnew%voroidx(objnew%k+1,:) = 0

IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)

!PRINT*,'zj',zj
!PRINT*,'znew',zdel,idel
!PRINT*,'zjp1',zjp1
!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
!! Apply only even-numbered order statistics in prior:
  objnew%logPr = 2._RP*LOG(hmx-hmin)-LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
!! Apply only Poisson prior:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
!! Apply Poisson prior with ENOS:
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+2._RP*LOG(hmx-hmin)- &
                 LOG(2._RP*obj%k*(2._RP*obj%k+1._RP)) + &
                 LOG(zjp1-zj)-LOG(zdel-zj)-LOG(zjp1-zdel)
ENDIF

END SUBROUTINE DEATH_FULL
!=======================================================================

SUBROUTINE BIRTH_FULL(obj,objnew)
!!
!! This birthes a full node with all parameters. Birth is based on the 
!! background value at the position of the new node which are then 
!! perturbed with a Gaussian proposal.
!!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)                 :: iznew, ipar, ivo, nnew
TYPE(objstruc)                   :: obj,objnew
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=RP),DIMENSION(obj%k+1) :: ztmp
REAL(KIND=RP)                    :: znew,zj,zjp1
INTEGER(KIND=IB),DIMENSION(NPL-1):: idxran
INTEGER(KIND=IB),DIMENSION(NPL)  :: voroidx

objnew = obj
objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!!
!! Sample No. new parameters for node from [1,NPL-1]
!! THIS CODE IS FIXED TO ALWAYS BITH ALL PARAMETERS!
!!
!idxran = RANDPERM(NPL-1)
!nnew = idxran(1)
nnew = NPL-1
idxran = 0
idxran = RANDPERM(NPL-1)+1  !! First nnew elements give index of 
                            !! parameters to be perturbed.
!!
!! voroidx identifies live parameters
!! THIS CODE IS FIXED TO ALWAYS BITH ALL PARAMETERS!
!!
voroidx = 1
!voroidx(1) = 1
!voroidx(idxran(1:nnew)) = 1

!!
!! Draw new z
!!
CALL RANDOM_NUMBER(ran_uni)
znew = maxpert(1)*ran_uni
!!
!! Insert new node at bottom of stack...
!!
!! 1) Sample depth from uniform prior:
objnew%voro(objnew%k,1)   = znew
objnew%voroidx(objnew%k,:)= voroidx
!objnew%gvoroidx           = voroidx(2:NPL)

!! 2) Sample new node parameters from prior
DO ipar = 2,NPL
  !!
  !! Gaussian proposal
  !!
  IF(objnew%voroidx(objnew%k,ipar) == 1)THEN
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(objnew%k,ipar) = minlim(ipar) + maxpert(ipar)*ran_uni
  ENDIF
ENDDO
!IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
!IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)

CALL INTERPLAYER_novar(objnew)

ztmp = objnew%ziface(1:objnew%k) - znew
!!
!! Find new interface index
!!
iznew = 0
DO ivo = 1,obj%k
   IF(ztmp(ivo) == 0._RP) iznew = ivo+1
ENDDO
zj = objnew%voro(iznew-1,1)
IF(iznew > obj%k)THEN
  zjp1 = hmx
ELSE
  zjp1 = objnew%voro(iznew+1,1)
ENDIF

!PRINT*,'zj',zj
!PRINT*,'znew',znew,iznew
!PRINT*,'zjp1',zjp1
!iznew = iznew + 1
!!
!! Compute prior for even-numbered order statistics (Green 1995)
!! znew is new interface, zj (j) is interface above and zjp1 is interface below (j+1)
!!
IF(ENOS == 0 .AND. IPOIPR == 0)THEN
  objnew%logPr = 0._RP
ELSEIF(ENOS == 1 .AND. IPOIPR == 0)THEN
  objnew%logPr = LOG(2._RP*obj%k+2._RP)+LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ELSEIF(ENOS == 0 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))
ELSEIF(ENOS == 1 .AND. IPOIPR == 1)THEN
  objnew%logPr = LOG(pk(objnew%k))-LOG(pk(obj%k))+LOG(2._RP*obj%k+2._RP)+ &
                 LOG(2._RP*obj%k+3._RP)-2._RP*LOG(hmx-hmin) + &
                 LOG(znew-zj)+LOG(zjp1-znew)-LOG(zjp1-zj)
ENDIF


RETURN
END SUBROUTINE BIRTH_FULL
!=======================================================================

SUBROUTINE DEATH_VARPAR(obj,objnew)
!!=======================================================================
!!
!! This DEATH is the reverse for Birth from prior.
!!
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: idel,ivo,ipar
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NFPMX):: idxdeath
REAL(KIND=RP)                    :: ran_uni
REAL(KIND=DRP),DIMENSION(NLMX,NPL):: tmpsort

objnew = obj
objnew%k   = obj%k - 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!! 1) Pick random node:
idxdeath = 0
idxdeath(1:obj%k-1) = RANDPERM(obj%k-1)
idxdeath = idxdeath + 1
idel = idxdeath(1)

!! 2) Save location and parameters of node
objnew%gvoroidx = obj%voroidx(idel,2:NPL)

!! 3) delete node and re-interpolate Voronoi model
objnew%voro(idel,:) = 0._RP
objnew%voroidx(idel,:) = 0

tmpsort = 0._DRP
tmpsort = REAL(objnew%voro(1:obj%k,:),DRP)
CALL QSORTC2D(tmpsort(1:obj%k,:),objnew%voroidx(1:obj%k,:))
objnew%voro(1:obj%k,:) = REAL(tmpsort,RP)
!CALL QSORTC2D(objnew%voro(1:obj%k,:),objnew%voroidx(1:obj%k,:))

!! 4) Quicksort stores deleted array (all zeros) first.
!!    Hence, move all one up unless deleted node is last node.
objnew%voro(1:objnew%k,:) = objnew%voro(2:obj%k,:)
objnew%voroidx(1:objnew%k,:) = objnew%voroidx(2:obj%k,:)
objnew%voro(objnew%k+1,:) = 0._RP
objnew%voroidx(objnew%k+1,:) = 0

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
END SUBROUTINE DEATH_VARPAR
!=======================================================================

SUBROUTINE BIRTH_VARPAR(obj,objnew)
!!
!! This birthes a full node with variable No. parameters. Birth is 
!! based on the background value at the position of the new node which 
!! are then perturbed with a Gaussian proposal.
!!
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)              :: iznew, ipar, ivo, nnew
TYPE(objstruc)                :: obj,objnew
REAL(KIND=RP)                 :: znew,ran_uni
REAL(KIND=RP),DIMENSION(obj%nunique):: ztmp
INTEGER(KIND=IB),DIMENSION(NPL-1):: idxran
INTEGER(KIND=IB),DIMENSION(NPL):: voroidx

objnew = obj
objnew%k   = obj%k + 1
objnew%NFP = (objnew%k * NPL) + (NPL-1)

!!
!! Sample No. new parameters for node from [1,NPL-1]
!!
idxran = RANDPERM(NPL-1)
nnew = idxran(1)
idxran = 0
idxran = RANDPERM(NPL-1)+1  !! First nnew elements give index of 
                            !! parameters to be perturbed.
!!
!! voroidx identifies live parameters
!!
voroidx = 0
voroidx(1) = 1
voroidx(idxran(1:nnew)) = 1

!!
!! Draw new z
!!
CALL RANDOM_NUMBER(ran_uni)
znew = maxpert(1)*ran_uni
ztmp = obj%ziface(1:obj%nunique) - znew
!!
!! Insert new node at bottom of stack...
!!
!! 1) Sample depth from uniform prior:
objnew%voro(objnew%k,1)   = znew
objnew%voroidx(objnew%k,:)= voroidx
objnew%gvoroidx           = voroidx(2:NPL)

!! 2) Sample new node parameters from prior
DO ipar = 2,NPL
  !!
  !! Gaussian proposal
  !!
  IF(objnew%voroidx(objnew%k,ipar) == 1)THEN
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(objnew%k,ipar) = minlim(ipar) + maxpert(ipar)*ran_uni
  ENDIF
ENDDO
!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
RETURN
END SUBROUTINE BIRTH_VARPAR
!!=======================================================================

SUBROUTINE DEATH_SINGLE(obj,objnew,ipick)
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
USE qsort_c_module
IMPLICIT NONE
INTEGER(KIND=IB)                 :: ipar,itmp,ipick,idel
INTEGER(KIND=IB)                 :: NPLIVE
TYPE(objstruc)                   :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NPL)  :: idxtmp,idxran

!! 1) Pick random live parameter
idxtmp = 0
itmp = 1
DO ipar = 2,NPL
  IF(obj%voroidx(ipick,ipar) == 1)THEN
    idxtmp(itmp) = ipar
    NPLIVE = itmp
    itmp = itmp + 1
  ENDIF
ENDDO
idxran = 0
idxran(1:NPLIVE) = RANDPERM(NPLIVE)
idel = idxtmp(idxran(1))

!! 2) delete parameter on node
objnew%voro(ipick,idel) = 0._RP
objnew%voroidx(ipick,idel) = 0

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
!PRINT*,'DEATH:',idel,ipick
!PRINT*,obj%voro(ipick,:)
!PRINT*,obj%voroidx(ipick,:)
!PRINT*,objnew%voro(ipick,:)
!PRINT*,objnew%voroidx(ipick,:)
!PRINT*,''

END SUBROUTINE DEATH_SINGLE
!!=======================================================================

SUBROUTINE BIRTH_SINGLE(obj,objnew,ipick)
!!
!! This birthes a single parameter onto node. Birth is based on the 
!! background value at the position of the new node which is then
!! perturbed with a Gaussian proposal.
!!
!!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ipar,itmp,ipick,inew
INTEGER(KIND=IB)               :: NPFREE
TYPE(objstruc)                 :: obj,objnew
INTEGER(KIND=IB),DIMENSION(NPL):: idxtmp,idxran
REAL(KIND=RP)                  :: ran_uni

!! 1) Randomly pick free parameter type
idxtmp = 0
itmp = 1
DO ipar = 2,NPL
  IF(obj%voroidx(ipick,ipar) == 0)THEN
    idxtmp(itmp) = ipar
    NPFREE = itmp
    itmp = itmp + 1
  ENDIF
ENDDO
idxran = 0
idxran(1:NPFREE) = RANDPERM(NPFREE)
inew = idxtmp(idxran(1))

!! 2) Sample new value from prior
CALL RANDOM_NUMBER(ran_uni)
objnew%voro(ipick,inew) = minlim(inew) + maxpert(inew)*ran_uni
objnew%voroidx(ipick,inew) = 1

!IF(IVORO == 0)THEN
!  CALL INTERPLAYER(objnew)
!ELSE
!  CALL INTERPVORO(objnew)
!ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
!PRINT*,'BIRTH:',inew
!PRINT*,obj%voro(ipick,:)
!PRINT*,obj%voroidx(ipick,:)
!PRINT*,objnew%voro(ipick,:)
!PRINT*,objnew%voroidx(ipick,:)
!PRINT*,''
RETURN
END SUBROUTINE BIRTH_SINGLE
!=======================================================================
SUBROUTINE TEMPSWP_MH(obj1,obj2)
!==============================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
TYPE(objstruc)              :: obj1
TYPE(objstruc)              :: obj2
TYPE(objstruc)              :: objtmp1,objtmp2
REAL(KIND=RP)               :: ran_uni1
REAL(KIND=RP)               :: logratio
REAL(KIND=RP)               :: betaratio

betaratio = obj2%beta-obj1%beta
logratio  = betaratio*(obj1%logL-obj2%logL)
CALL RANDOM_NUMBER(ran_uni1)
IF(ran_uni1 <= EXP(logratio))THEN
  !! ACCEPT SWAP
!  IF(obj(ic1idx)%k /= obj(ic2idx)%k)THEN
!    PRINT*,'ACCEPTED',EXP(logratio*betaratio),ran_uni2
!    WRITE(ulog,204)'ic1:',ic1idx,obj(ic1idx)%k,logP1,obj(ic1idx)%logL,logratio,betaratio,EXP(logratio*betaratio),ran_uni2
!    WRITE(ulog,204)'ic2:',ic2idx,obj(ic2idx)%k,logP2,obj(ic2idx)%logL
!  ENDIF
  objtmp1       = obj1
  objtmp2       = obj2
  obj1   = objtmp2
  obj2   = objtmp1
  !! Temperature does not swap
  obj1%beta = objtmp1%beta
  obj2%beta = objtmp2%beta
  !! BD acceptance counters do not swap
!  obj1%iaccept_bd = objtmp1%iaccept_bd
!  obj2%iaccept_bd = objtmp2%iaccept_bd
!  obj1%ipropose_bd = objtmp1%ipropose_bd
!  obj2%ipropose_bd = objtmp2%ipropose_bd
  !! Voro acceptance counters do not swap
!  obj1%iacceptvoro = objtmp1%iacceptvoro
!  obj2%iacceptvoro = objtmp2%iacceptvoro
!  obj1%iproposevoro = objtmp1%iproposevoro
!  obj2%iproposevoro = objtmp2%iproposevoro
  !! Step sizes do not swap
!  obj1%pertsd = objtmp1%pertsd
!  obj2%pertsd = objtmp2%pertsd
  ncswap    = ncswap+1_IB
ENDIF
ncswapprop    = ncswapprop+1_IB

RETURN
END SUBROUTINE TEMPSWP_MH
!=======================================================================

SUBROUTINE PROPOSAL(obj,objnew,ivo,iwhich,factor)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: ivo,iwhich,ipar
TYPE(objstruc) :: obj,objnew,objtmp
REAL(KIND=RP)  :: ran_uni, factor
REAL(KIND=RP)                    :: zp,zj,zjp1,zjm1

objnew = obj
objnew%logPr = 0._RP
!!
!! CAUCHY proposal
!!
CALL RANDOM_NUMBER(ran_uni)
!! Cauchy:
IF(iwhich /= 1)THEN
  objnew%voro(ivo,iwhich) = obj%voro(ivo,iwhich) + fact/factor*pertsd(iwhich)*TAN(PI2*(ran_uni-0.5_RP))
ELSE
  !!
  !! Compute prior for even-numbered order statistics (Green 1995) where 
  !! znew is perturbed interface, zj (j) is original interface, 
  !! zjp1 is interface below (j+1) and zjm1 is interface above.
  !!
  IF(ENOS == 0)THEN
    !! CAUCHY proposal
    CALL RANDOM_NUMBER(ran_uni)
    objnew%voro(ivo,iwhich) = obj%voro(ivo,iwhich) + fact/factor*pertsd(iwhich)*TAN(PI2*(ran_uni-0.5_RP))
    objnew%logPr = 0._RP
  ELSE
    CALL RANDOM_NUMBER(ran_uni)
    zj = obj%voro(ivo,iwhich)
    zjm1 = obj%voro(ivo-1,iwhich)
    IF(ivo == obj%k)THEN
      zjp1 = hmx
    ELSE
      zjp1 = obj%voro(ivo+1,iwhich)
    ENDIF
    !! sample uniform:
    zp = zjm1+ran_uni*(zjp1-zjm1)
    objnew%voro(ivo,iwhich) = zp
    !! Apply even-numbered order statistics in prior:
    objnew%logPr = LOG(zjp1-zp)+LOG(zp-zjm1)-LOG(zjp1-zj)-LOG(zj-zjm1)
!    PRINT*,''
!    PRINT*,'ic',iwhich
!    PRINT*,'zjm1',zjm1
!    PRINT*,'zj',zj
!    PRINT*,'zp',zp
!    PRINT*,'zjp1',zjp1
!    PRINT*,'logPr',objnew%logPr
!    PRINT*,''
  ENDIF
ENDIF
!!
!! If node position comes up as negative, INTERPLAYER breaks...
!!
IF(iwhich == 1)THEN
  objnew%voro(ivo,iwhich) = ABS(objnew%voro(ivo,iwhich))
ENDIF
IF(I_VARPAR == 0)CALL INTERPLAYER_novar(objnew)
IF(I_VARPAR == 1)CALL INTERPLAYER(objnew)
RETURN
END SUBROUTINE PROPOSAL
!=======================================================================

SUBROUTINE PROPOSAL_AR(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arpar(iwhich) = ran_uni*(maxlimar(iwhich)-minlimar(iwhich))+minlimar(iwhich)
  objnew%idxar(iwhich) = 1
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arpar(iwhich) = minlimar(iwhich)-1._RP
  objnew%idxar(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arpar(iwhich) = obj%arpar(iwhich) + pertarsd(iwhich)*ran_nor
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_AR
!=======================================================================

SUBROUTINE PROPOSAL_ARRF(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arpar(iwhich) = ran_uni*(maxlimar(iwhich)-minlimar(iwhich))+minlimar(iwhich)
  objnew%idxar(iwhich) = 1
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arpar(iwhich) = minlimar(iwhich)-1._RP
  objnew%idxar(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arpar(iwhich) = obj%arpar(iwhich) + pertarsd(iwhich)*ran_nor
  IF(((objnew%arpar(iwhich) - minlimar(iwhich)) < 0._RP).OR. &
     ((maxlimar(iwhich) - objnew%arpar(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_ARRF
!=======================================================================

SUBROUTINE PROPOSAL_ARRT(obj,objnew,iwhich,arptype)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich,arptype
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor,ran_uni

!! Birth: sample uniform from prior
IF(arptype == 1)THEN
  CALL RANDOM_NUMBER(ran_uni)
  objnew%arparRT(iwhich) = ran_uni*(maxlimarRT(iwhich)-minlimarRT(iwhich))+minlimarRT(iwhich)
  objnew%idxarRT(iwhich) = 1
  IF(((objnew%arparRT(iwhich) - minlimarRT(iwhich)) < 0._RP).OR. &
     ((maxlimarRT(iwhich) - objnew%arparRT(iwhich)) < 0._RP))ioutside = 1
ENDIF
!! Death
IF(arptype == 2)THEN
  objnew%arparRT(iwhich) = minlimarRT(iwhich)-1._RP
  objnew%idxarRT(iwhich) = 0
ENDIF
!! Perturb
IF(arptype == 3)THEN
  CALL GASDEVJ(ran_nor)
  objnew%arparRT(iwhich) = obj%arparRT(iwhich) + pertarsdRT(iwhich)*ran_nor
  IF(((objnew%arparRT(iwhich) - minlimarRT(iwhich)) < 0._RP).OR. &
     ((maxlimarRT(iwhich) - objnew%arparRT(iwhich)) < 0._RP))ioutside = 1
ENDIF

RETURN
END SUBROUTINE PROPOSAL_ARRT
!=======================================================================

!SUBROUTINE PROPOSAL_SDH(obj,objnew,iwhich)
!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!!INTEGER(KIND=IB) :: iwhich
!TYPE(objstruc) :: obj,objnew
!REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
!CALL GASDEVJ(ran_nor)
!objnew%sdparR(iwhich) = obj%sdparR(iwhich) + pertsdsd(iwhich)*ran_nor
!IF(((objnew%sdparR(iwhich) - minlimsd(iwhich)) < 0._RP).OR.((maxlimsd(iwhich) - objnew%sdparR(iwhich)) < 0._RP))ioutside = 1

!RETURN
!END SUBROUTINE PROPOSAL_SDH
!=======================================================================

!SUBROUTINE PROPOSAL_SDV(obj,objnew,iwhich)
!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: iwhich
!TYPE(objstruc) :: obj,objnew
!REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
!CALL GASDEVJ(ran_nor)
!objnew%sdparV(iwhich) = obj%sdparV(iwhich) + pertsdsd(iwhich)*ran_nor
!IF(((objnew%sdparV(iwhich) - minlimsd(iwhich)) < 0._RP).OR. & 
!   ((maxlimsd(iwhich) - objnew%sdparV(iwhich)) < 0._RP))ioutside = 1

!RETURN
!END SUBROUTINE PROPOSAL_SDV
!=======================================================================

!SUBROUTINE PROPOSAL_SDT(obj,objnew,iwhich)
!!=======================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: iwhich
!TYPE(objstruc) :: obj,objnew
!REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
!CALL GASDEVJ(ran_nor)
!objnew%sdparT(iwhich) = obj%sdparT(iwhich) + pertsdsd(iwhich)*ran_nor
!IF(((objnew%sdparT(iwhich) - minlimsd(iwhich)) < 0._RP).OR.((maxlimsd(iwhich) - objnew%sdparT(iwhich)) < 0._RP))ioutside = 1

!RETURN
!END SUBROUTINE PROPOSAL_SDT
!=======================================================================

SUBROUTINE PROPOSAL_SDRT(obj,objnew,iwhich)
!=======================================================================

USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB) :: iwhich
TYPE(objstruc) :: obj,objnew
REAL(KIND=RP)  :: ran_nor

!!
!! Gaussian proposal
!!
CALL GASDEVJ(ran_nor)
objnew%sdparRT(iwhich) = obj%sdparRT(iwhich) + pertsdsdRT(iwhich)*ran_nor
IF(((objnew%sdparRT(iwhich) - minlimsdRT(iwhich)) < 0._RP).OR. &
   ((maxlimsdRT(iwhich) - objnew%sdparRT(iwhich)) < 0._RP))ioutside = 1

RETURN
END SUBROUTINE PROPOSAL_SDRT

!==============================================================================

SUBROUTINE CHECKBOUNDS(obj)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ivo,ipar,ilay,ncra
INTEGER(KIND=IB)               :: ih,OS,ihamrej
TYPE(objstruc)                 :: obj
REAL(KIND=RP)                  :: vspmin,vspmax,zhere
REAL(KIND=RP),DIMENSION(NFPMX2):: par

DO ilay = 1,obj%nunique
   IF(hmin > obj%hiface(ilay))ioutside = 1
   IF(maxlim(1) < obj%ziface(ilay))ioutside = 1
ENDDO

DO ivo = 1,obj%k
!! Starts at ivo = 2, since 1st voro node is fixed and gives half space if k = 0
  IF(ivo > 1)THEN
    IF((obj%voro(ivo,1) < 0._RP).OR.(obj%voro(ivo,1) > maxlim(1)))THEN
      ioutside = 1
      IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=', &
                     ivo,1,minlim(1),maxlim(1),obj%voro(ivo,1)
    ENDIF
  ENDIF
  DO ipar = 2,NPL
    IF(obj%voroidx(ivo,ipar) == 1)THEN
    IF(((obj%voro(ivo,ipar) - minlim(ipar)) < 0._RP).OR. & 
      ((maxlim(ipar) - obj%voro(ivo,ipar)) < 0._RP))THEN
      ioutside = 1
      IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=',ivo,ipar,&
                     minlim(ipar),maxlim(ipar),obj%voro(ivo,ipar)
    ENDIF
    ENDIF
  ENDDO 
ENDDO 
201 FORMAT(a27,2I4,3F18.4)

RETURN
END SUBROUTINE CHECKBOUNDS
!=======================================================================

SUBROUTINE CHECKBOUNDS2(obj,ivo,ipar)
!=======================================================================
USE DATA_TYPE
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB)               :: ivo,ipar,ilay,ncra
INTEGER(KIND=IB)               :: ih,OS,ihamrej
TYPE(objstruc)                 :: obj
REAL(KIND=RP)                  :: vspmin,vspmax,zhere
REAL(KIND=RP),DIMENSION(NFPMX2):: par

DO ilay = 1,obj%nunique
   IF(hmin > obj%hiface(ilay))ioutside = 1
   IF(maxlim(1) < obj%ziface(ilay))ioutside = 1
ENDDO

!! Starts at ivo = 2, since 1st voro node is fixed and gives half space if k = 0
IF(ivo > 1)THEN
  IF((obj%voro(ivo,1) < 0._RP).OR.(obj%voro(ivo,1) > maxlim(1)))THEN
    ioutside = 1
    IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=', &
                   ivo,1,minlim(1),maxlim(1),obj%voro(ivo,1)
  ENDIF
ENDIF
IF(obj%voroidx(ivo,ipar) == 1)THEN
  IF(((obj%voro(ivo,ipar) - minlim(ipar)) < 0._RP).OR. & 
    ((maxlim(ipar) - obj%voro(ivo,ipar)) < 0._RP))THEN
    ioutside = 1
    IF(rank == src)WRITE(6,201) 'ivo, ipar, min, max, value=',ivo,ipar,&
                   minlim(ipar),maxlim(ipar),obj%voro(ivo,ipar)
  ENDIF
ENDIF
201 FORMAT(a27,2I4,3F18.4)

RETURN
END SUBROUTINE CHECKBOUNDS2
!=======================================================================
!SUBROUTINE UPDATE_SDAVE(obj1,obj2)
!!==============================================================================
!USE DATA_TYPE
!USE RJMCMC_COM
!IMPLICIT NONE
!INTEGER(KIND=IB) :: irf
!TYPE(objstruc)   :: obj1
!TYPE(objstruc)   :: obj2
!
!DO irf = 1,NRF1
!  sdbuf(1,isdbuf)  = obj1%sdparR()
!ENDDO
!DO irf = 1,NRF1
!  obj1%sdaveH(irf)   = SUM(sdbuf(1,:))/REAL(NBUF,RP)
!  obj1%sdaveV(irf)   = SUM(sdbuf(2,:))/REAL(NBUF,RP)
!  obj1%sdaveRT(irf) = SUM(sdbuf(3,:))/REAL(NBUF,RP)
!  obj2%sdaveH(irf)   = SUM(sdbuf(1,:))/REAL(NBUF,RP)
!  obj2%sdaveV(irf)   = SUM(sdbuf(2,:))/REAL(NBUF,RP)
!  obj2%sdaveRT(irf) = SUM(sdbuf(3,:))/REAL(NBUF,RP)
!ENDDO
!
!RETURN
!END SUBROUTINE UPDATE_SDAVE
!=======================================================================

SUBROUTINE PARALLEL_SEED()
!!
!!  Ensure unique random seed for each CPU
!!
!=======================================================================
USE MPI
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed
INTEGER(KIND=IB)                              :: iseedsize
INTEGER(KIND=IB), DIMENSION(:,:), ALLOCATABLE :: iseeds
REAL(KIND=RP),    DIMENSION(:,:), ALLOCATABLE :: rseeds
INTEGER(KIND=IB), DIMENSION(:),   ALLOCATABLE :: iseed1
REAL(KIND=RP) :: ran_uni

CALL RANDOM_SEED
CALL RANDOM_SEED(SIZE=iseedsize)
ALLOCATE( iseed1(iseedsize) )
IF(ISETSEED == 1)THEN
   iseed1 = (/2303055,     2435432,     5604058,     4289794,     3472290, &
      7717070,      141180,     3783525,     3087889,     4812786,     3028075, &
      3712062,     6316731,      436800,     7957708,     2055697,     1944360, &
      1222992,     7537775,     7769874,     5588112,     7590383,     1426393, &
      1753301,     7681841,     2842400,     4411488,     7304010,      497639, &
      4978920,     5345495,      754842,     7360599,     5776102/)
   CALL RANDOM_SEED(PUT=iseed1)
ELSE
   CALL RANDOM_SEED(GET=iseed1)
!   IF(rank==src)WRITE(6,*) 'Master seed:',iseed1
ENDIF

ALLOCATE( iseed(iseedsize), rseeds(iseedsize,NTHREAD), iseeds(iseedsize,NTHREAD) )

iseed = 0
rseeds = 0._RP
iseeds = 0
IF(rank == src)THEN
   CALL RANDOM_NUMBER(rseeds)
   iseeds = -NINT(rseeds*1000000._RP)
ENDIF
DO i = 1,iseedsize
   CALL MPI_BCAST( iseeds(i,:), NTHREAD, MPI_INTEGER, src, MPI_COMM_WORLD, ierr )
ENDDO
iseed = iseeds(:,rank+1)

!!
!! Write seeds to seed logfile:
!!
IF(rank == src)THEN
   OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
   ACTION='WRITE',POSITION='REWIND',RECL=1024)
   WRITE(50,*) 'Rank: ',rank
   WRITE(50,201) iseed
   WRITE(50,*) ''
   CLOSE(50)
ENDIF
CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
DO i = 1,NTHREAD-1
   IF(rank == i)THEN
      OPEN(UNIT=50,FILE=seedfile,FORM='formatted',STATUS='UNKNOWN', &
      ACTION='WRITE',POSITION='APPEND',RECL=1024)
      WRITE(50,*) 'Rank: ',rank
      WRITE(50,201) iseed
      WRITE(50,*) ''
      CLOSE(50)
   ENDIF
   CALL MPI_BARRIER(MPI_COMM_WORLD,ierr)
ENDDO
CALL RANDOM_SEED(PUT=iseed)

DO i = 1,100
   CALL RANDOM_NUMBER(ran_uni)
ENDDO
DO i = 1,CEILING(ran_uni*10000._RP)
   CALL RANDOM_NUMBER(ran_uni)
ENDDO

201   FORMAT(50I10)
END SUBROUTINE PARALLEL_SEED
!=======================================================================

SUBROUTINE SAVESAMPLE(objm,ikeep,isource1,isource2,tmaster)
!!=======================================================================
!!
!! Exchanging and saving posterior samples
!!
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
TYPE (objstruc),DIMENSION(2)     :: objm     ! Objects on master node
INTEGER(KIND=IB)                 :: ipar,ipar2,ivo,ic,imcmc,ikeep,jidx
INTEGER(KIND=IB)                 :: isource,isource1,isource2
REAL(KIND=RP),DIMENSION(NFPMX)   :: tmpvoro  ! Temporary Voronoi cell array for saving
!INTEGER(KIND=RP),DIMENSION(NFPMX):: tmpprop,tmpacc
REAL(KIND=RP),DIMENSION(NPL)     :: tmp  ! Temporary Voronoi cell array for saving
REAL(KIND=RP) :: tmaster
!!
!! Save object into global sample array
!!
DO ic = 1,2
  IF(ic == 1)isource = isource1
  IF(ic == 2)isource = isource2
  IF(objm(ic)%beta > 1._RP/dTlog)THEN
    !!
    !! Write to stdout
    !!
    IF(MOD(imcmc1,20*NKEEP) == 0)THEN
      tsave2 = MPI_WTIME()
                        WRITE(*,213)'          iaccept_bd = ',objm(ic)%iaccept_bd,objm(ic)%beta
      IF(IEXCHANGE == 1)WRITE(*,215)'       T swap accept = ',REAL(ncswap,RP)/REAL(ncswapprop,RP),isource1,isource2
                        WRITE(*,216)'Total time for block = ',tsave2-tsave1,'s'
      WRITE(*,*) ''
      WRITE(*,203) '     imcmc1,         logL,  k,  iacc_bd, iacc_bds, source, time(slave)'
      WRITE(*,203) '----------------------------------------------------------------------'
      tsave1 = MPI_WTIME()
    ENDIF
    tmpvoro = 0._RP
     DO ivo = 1,objm(ic)%k
       ipar = (ivo-1)*NPL+1
       DO ipar2 = 1,NPL
         IF(objm(ic)%voroidx(ivo,ipar2)  == 1)THEN
           tmpvoro(ipar) = objm(ic)%voro(ivo,ipar2)
         ELSE
           tmpvoro(ipar) = -100._RP
         ENDIF
         ipar = ipar + 1
       ENDDO
!      tmpprop(ipar:ipar+NPL-1) = objm(ic)%iproposevoro(ivo,:)
!      tmpacc(ipar:ipar+NPL-1)  = objm(ic)%iacceptvoro(ivo,:)
     ENDDO
!    sample(ikeep,:) =  (/ objm(ic)%logL,objm(ic)%logPr,objm(ic)%lognorm,REAL(objm(ic)%k,RP),tmpvoro,&
!                          objm(ic)%arpar,objm(ic)%sdpar,&
!                          REAL(tmpacc,RP)/REAL(tmpprop,RP),REAL(objm(ic)%iaccept_bd,RP),&
!                          REAL(objm(ic)%ipropose_bd,RP),REAL(i_bd,RP),objm(ic)%tcmp,REAL(isource,RP) /)
    print *, 'Filling sample array'
     sample(ikeep,:) =  (/ objm(ic)%logL, objm(ic)%logPr, objm(ic)%tcmp, REAL(objm(ic)%k,RP), & ! 4 parameters
                           tmpvoro,objm(ic)%sdparRT,objm(ic)%arpar,objm(ic)%arparRT, &
                           REAL(iaccept,RP)/REAL(iaccept+ireject,RP),REAL(objm(ic)%iaccept_bd,RP),&
                           REAL(objm(ic)%ireject_bd,RP),REAL(objm(ic)%iaccept_bds,RP),REAL(ic,RP),REAL(rank,RP) /)

    IF(MOD(imcmc1,NKEEP)==0)THEN
      WRITE(*,202) '       ',imcmc1,objm(ic)%logL,objm(ic)%k,objm(ic)%ireject_bd,objm(ic)%iaccept_bds,isource,objm(ic)%tcmp
    ENDIF
    ikeep = ikeep + 1_IB
    imcmc1 = imcmc1 + 1_IB

    !!
    !! Write to file
    !!
    IF(ikeep > NKEEP)THEN
      DO jidx=1,NKEEP

        WRITE(usample,207) sample(jidx,:)
      ENDDO
      CALL FLUSH(usample)
!      WRITE(ustep,207) (objm(ic)%pertsd(ivo,:),ivo=1,objm(ic)%k)
!      CALL FLUSH(ustep)
      sample  = 0._RP
      ikeep   = 1
    ENDIF
  ENDIF
ENDDO
print *, ' Savesample end'

CALL FLUSH(6)
202 FORMAT(A3,I8,1F14.4,I4,I10,I10,I8,1f12.4)
203 FORMAT(A69)
207 FORMAT(500ES18.8)
213 FORMAT(a23,I4,F8.4) ! Not sure here why this happens to be an error !
214 FORMAT(a23,I6,a9,I4)
215 FORMAT(a23,F8.4,I4,I4)
216 FORMAT(a23,F8.2,a)


RETURN
END SUBROUTINE SAVESAMPLE

!=======================================================================
SUBROUTINE SAVEREPLICA(obj,predfile,obsfile)

!=======================================================================
USE RJMCMC_COM
IMPLICIT NONE

INTEGER(KIND=IB) :: i
TYPE (objstruc)  :: obj      ! Best object
CHARACTER(64) :: predfile, obsfile

WRITE(6,*) 'Global best model:'
IF(rank ==src) CALL PRINTPAR(obj)
WRITE(6,*) 'Global best logL = ',obj%logL

IF(I_RT == 1)THEN
  OPEN(UNIT=50,FILE=predfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE
    WRITE(50,208) obj%DpredRT(i,:)
  ENDDO
  WRITE(6,*)'Done writing predicted Ray Tracer TTimes.'
  CLOSE(50)

  OPEN(UNIT=50,FILE=obsfile,FORM='formatted',STATUS='REPLACE', &
  ACTION='WRITE',RECL=8192)
  DO i = 1,NMODE
    WRITE(50,208) obj%DobsRT(i,:)
  ENDDO
  WRITE(6,*)'Done writing observed Ray Tracer TTimes.'
  CLOSE(50)

  !OPEN(UNIT=50,FILE=arfileRT,FORM='formatted',STATUS='REPLACE', &
  !ACTION='WRITE',RECL=8192)
  !DO i = 1,NMODE
  !  WRITE(50,208) obj%DarRT(i,:)
  !ENDDO
  !CLOSE(50)

ENDIF
208 FORMAT(5000ES20.10)
RETURN
END SUBROUTINE SAVEREPLICA
!=======================================================================


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   MISCELLANEOUS MATH FUNCTIONS BELOW                            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1


FUNCTION RANDPERM(num)
!==============================================================================
USE data_type, ONLY : IB, RP
IMPLICIT NONE
INTEGER(KIND=IB), INTENT(IN) :: num
INTEGER(KIND=IB) :: numb, i, j, k
INTEGER(KIND=IB), DIMENSION(num) :: randperm
REAL(KIND=RP), DIMENSION(num) :: rand2
INTRINSIC RANDOM_NUMBER
CALL RANDOM_NUMBER(rand2)
DO i=1,num
   numb=1
   DO j=1,num
      IF (rand2(i) > rand2(j)) THEN
           numb=numb+1
      END IF
   END DO
   DO k=1,i-1
      IF (rand2(i) <= rand2(k) .AND. rand2(i) >= rand2(k)) THEN
           numb=numb+1
      END IF
   END DO
   randperm(i)=numb
END DO
RETURN
END FUNCTION RANDPERM
!====================================================================

SUBROUTINE GASDEVJ(harvest)
!====================================================================
USE nrtype
USE nr
IMPLICIT NONE
REAL(DP), INTENT(OUT) :: harvest
REAL(DP) :: rsq,v1,v2
REAL(DP), SAVE :: g
LOGICAL, SAVE :: gaus_stored=.FALSE.
IF (gaus_stored) THEN
   harvest=g
   gaus_stored=.FALSE.
ELSE
   DO
      CALL RANDOM_NUMBER(v1)
      CALL RANDOM_NUMBER(v2)
      v1=2.0_DP*v1-1.0_DP
      v2=2.0_DP*v2-1.0_DP
      rsq=v1**2+v2**2
      IF (rsq > 0.0_DP .AND. rsq < 1.0_DP) EXIT
   END DO
   rsq=SQRT(-2.0_DP*LOG(rsq)/rsq)
   harvest=v1*rsq
   g=v2*rsq
   gaus_stored=.TRUE.
END IF
END SUBROUTINE GASDEVJ
!====================================================================

Function ASINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,ii,asinc
ii    = cmplx(0.,1.)
asinc = -ii*LOG(ii*z+SQRT(1.-z**2))
RETURN
END FUNCTION
!=======================================================================

Function COSC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,cosc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
cosc = CMPLX(COS(x)*COSH(y),-SIN(x)*SINH(y),RP)
RETURN
END FUNCTION
!=======================================================================

Function SINC(z)
!=======================================================================
USE DATA_TYPE, ONLY : IB, RP
COMPLEX(KIND=RP) :: z,sinc
REAL(KIND=RP)    :: x,y
x    = REAL(z)
y    = AIMAG(z)
sinc = CMPLX(SIN(x)*COSH(y),COS(x)*SINH(y),RP)
RETURN
END FUNCTION
!==============================================================================

FUNCTION CACOS(z)
!==============================================================================

USE DATA_TYPE
COMPLEX(KIND=RP) :: CACOS
COMPLEX(KIND=RP) :: z
REAL(KIND=RP) :: zrp1,zrm1,zi,zizi,a1,a2,a,b

!CACOS = -CMPLX(0._RP,1._RP,RP)*LOG(z+CMPLX(0._RP,1._RP,RP)*SQRT(1._RP-z*z))
!!
!! This version from IDL; much faster than above
!!
zrp1 = REAL(z,RP)+1._RP
zrm1 = zrp1-2._RP
zi = AIMAG(z)
zizi = zi*zi
a1 = 0.5_RP*SQRT(zrp1*zrp1 + zizi)
a2 = 0.5_RP*SQRT(zrm1*zrm1 + zizi)
a = a1+a2
b = a1- a2
IF(zi >= 0._RP)THEN
   CACOS = ACOS(b) - CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ELSE
   CACOS = ACOS(b) + CMPLX(0._RP,1._RP,RP)*LOG(a + SQRT(a*a - 1))
ENDIF

RETURN
END FUNCTION CACOS
!==============================================================================

FUNCTION ASINH(x)
!==============================================================================
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) :: ASINH
REAL(KIND=RP) :: x

ASINH = LOG(x+SQRT(x**2._RP+1))

RETURN
END FUNCTION ASINH
!==============================================================================

FUNCTION CSIN(z)
!==============================================================================
!! Complex sine (Jan's version)

USE DATA_TYPE
COMPLEX(KIND=RP) :: CSIN
COMPLEX(KIND=RP) :: z

CSIN =  (EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /CMPLX(0._RP,2._RP,RP)
RETURN
END FUNCTION CSIN
!==============================================================================

FUNCTION CTAN(z)
!==============================================================================
!! Complex TAN

USE DATA_TYPE
COMPLEX(KIND=RP) :: CTAN
COMPLEX(KIND=RP) :: z

CTAN =  -CMPLX(0._RP,1._RP,RP)*(EXP( CMPLX(0._RP,1._RP,RP)*z) -EXP(-CMPLX(0._RP,1._RP,RP)*z)) &
                            /(EXP( CMPLX(0._RP,1._RP,RP)*z)+EXP(-CMPLX(0._RP,1._RP,RP)*z))
RETURN
END FUNCTION CTAN

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

LOGICAL FUNCTION ISNAN(a)
!=======================================================================
USE DATA_TYPE
IMPLICIT NONE
REAL(KIND=RP) a
IF (a.NE.a) THEN
ISNAN = .TRUE.
ELSE
ISNAN = .FALSE.
END IF
RETURN
END
!=======================================================================

SUBROUTINE CONVT(x,y,z,Nx,Ny,Nz)
!!=======================================================================
!!
!! Time domain decon
!!
USE MPI
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
USE MPI
USE RJMCMC_COM
IMPLICIT NONE
INTEGER(KIND=IB):: Nx,Ny,Nz
INTEGER(KIND=IB):: irow,icol1,icol2,ix1,ix2
REAL(KIND=RP)   :: x(Nx), y(Ny), z(Nz),xr(Nx)
REAL(KIND=RP)   :: XMAT(Nz,Ny)
!! Lapack variables:
INTEGER(KIND=IB)          :: MRANK, LWORK, INFO,LDA, LDB
INTEGER(KIND=IB),PARAMETER:: LWMAX = 4000, NRHS = 1

!! These are always double precision for DGELLS to work
REAL(KIND=DRP)          :: SV(Ny),WORK(LWMAX)
REAL(KIND=DRP)          :: XMAT2(Nz,Ny),b(NZ)
REAL(KIND=DRP),PARAMETER:: RCOND = 1.e-12_DRP

b = REAL(z,DRP)
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
!! y = pinv(X)*z;
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
XMAT2 = REAL(XMAT,DRP)
LWORK = -1
!CALL DGELSS(Nz, Ny, NRHS, XMAT2, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
!LWORK = MIN( LWMAX, INT( WORK( 1 ) ) )
!! Carry out DGELSS with optimal length
!CALL DGELSS(Nz, Ny, NRHS, XMAT2, LDA, b, LDB, SV, RCOND, MRANK, WORK, LWORK, INFO)
!IF(INFO /= 0)WRITE(*,*)'WARNING DGELSS unstable!'

y = REAL(b(1:Ny),RP)
RETURN
END SUBROUTINE DCONVT
!=======================================================================
! This is the end my fiend...
! EOF

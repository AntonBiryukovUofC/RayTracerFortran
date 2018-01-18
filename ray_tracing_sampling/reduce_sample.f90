!==============================================================================

PROGRAM REDUCE_SAMPLE

!=======================================================================
IMPLICIT NONE
INTEGER(KIND=4), PARAMETER :: RP=KIND(0.0)
INTEGER :: stat,iexit,filebaselen
INTEGER :: iglob,ismp,isub,ibuf,ibuf2
INTEGER :: NBURNIN,NSMP,NSUB,NSAVE,ncount1,NBUF
REAL(KIND=RP),DIMENSION(:,:),ALLOCATABLE  :: sample
CHARACTER(len=64):: plotparfile       = 'plotpar.txt'
CHARACTER(len=64):: filebase
CHARACTER(len=64):: samplefilein,samplefileout

OPEN(UNIT=20,FILE='filebase.txt',FORM='formatted',STATUS='OLD',ACTION='READ')
READ(20,*) filebaselen
READ(20,*) filebase
CLOSE(20)
NBUF = 2000
samplefilein  = filebase(1:filebaselen) // '_voro_sample.txt'
samplefileout = filebase(1:filebaselen) // '_voro_sample2.txt'

OPEN(UNIT=32,FILE=plotparfile,FORM='formatted',ACTION='READ')
READ(32,*) NBURNIN
READ(32,*) NSMP
READ(32,*) NSUB
READ(32,*) NSAVE
READ(32,*) ncount1
CLOSE(32)

WRITE(*,*) 'ncount:',ncount1
ALLOCATE( sample(NBUF,ncount1) )

OPEN(UNIT=33,FILE=samplefilein,FORM='formatted',ACTION='READ')
iglob = 1
DO ismp = 1,NBURNIN
   READ(33,*)
   iglob = iglob + 1
ENDDO
WRITE(6,*) 'Discarded burn-in.'

OPEN(UNIT=51,FILE=samplefileout,FORM='formatted',STATUS='REPLACE', &
ACTION='WRITE',RECL=1024)

ismp = 0
ibuf = 0
iexit = 0
DO 
  ismp = ismp + 1
  ibuf = ibuf + 1
  READ(33,*,iostat=stat) sample(ibuf,:)
  if (stat /= 0) exit
  if(ibuf == NBUF)then
    DO ibuf2 = 1,NBUF
      WRITE(51,201) sample(ibuf2,:)
    ENDDO
    ibuf = 0
  endif
  iglob = iglob + 1
  DO isub = 1,NSUB-1
    READ(33,*,iostat=stat)
    if (stat /= 0)then
      iexit = 1
      exit
    endif
    iglob = iglob + 1
  ENDDO
  if (iexit /= 0)exit
  IF(MOD(ismp,10000)==0)WRITE(6,*)ismp
ENDDO
CALL FLUSH(51)
close(33)
close(51)

201 FORMAT(500ES18.8)
END PROGRAM REDUCE_SAMPLE
!!EOF

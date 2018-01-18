!==============================================================================
MODULE DATA_TYPE
   IMPLICIT NONE
   INTEGER(KIND=4), PARAMETER :: IB=4, RP=KIND(0.D0), DRP=KIND(0.D0), SP=KIND(0.0)
   REAL(KIND=RP),   PARAMETER :: PI2  = 3.141592653589793238462643383279502884197_RP
   COMPLEX(KIND=RP),PARAMETER :: cir  = CMPLX(0._RP,1._RP,KIND=RP)
   COMPLEX(KIND=SP),PARAMETER :: cis  = CMPLX(0._SP,1._SP,KIND=SP)
END MODULE DATA_TYPE
!==============================================================================

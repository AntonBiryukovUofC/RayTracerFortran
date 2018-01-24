!Fortran Test
program divisors
!This program finds the divisors of an integer input by the user.
!The divisors are printed to a file.
integer :: NSRC,idat
real,dimension(6) :: src_offset,src_depth
NSRC=5
!
!! Source depths - read from the file:
!!
OPEN(20,FILE="src.dat",FORM='formatted',STATUS='OLD',ACTION='READ')
!ndat = 0
DO idat=1,NSRC
   !READ(20,*,IOSTAT=io) obj%periods(1,idat),obj%DobsRT(1,idat)
   !READ(20,*,IOSTAT=io) sour(1,idat)
   READ(20,*) src_offset(idat), src_depth(idat)

ENDDO
print *,src_offset
print *,src_depth

CLOSE(20)! close the file
end
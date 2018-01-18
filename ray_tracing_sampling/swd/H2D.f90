subroutine H2D(V0h,V0d,nptref)
implicit none
include 'params.h'

real, intent(in) :: V0h(maxlay,2)
real, intent(out) :: V0d(maxlay,2)
integer, intent(in) :: nptref

real summ
integer i

V0d(1,1)=0
V0d(1,2)=V0h(1,2)

summ =0

do i=2,nptref
      summ=summ+V0h(i-1,1)
      V0d(i,1)=summ
      V0d(i,2)=V0h(i,2)
enddo
return
end

! Recursive Fortran 95 quicksort routine
! sorts real numbers into ascending numerical order
! Author: Juli Rew, SCD Consulting (juliana@ucar.edu), 9/03
! Based on algorithm from Cormen et al., Introduction to Algorithms,
! 1997 printing

! Made F conformant by Walt Brainerd

module qsort_c_module

implicit none
public :: QSORTC1D,QSORTC2D
private :: PARTITION1D,PARTITION2D
INTEGER(KIND=4), PARAMETER :: DP=KIND(0.D0), SP=KIND(0.0)

contains

recursive subroutine QSORTC1D(A)
  real(KIND=DP), intent(in out), dimension(:) :: A
  integer :: iq

  if(size(A) > 1) then
     call PARTITION1D(A, iq)
     call QSORTC1D(A(:iq-1))
     call QSORTC1D(A(iq:))
  endif
end subroutine QSORTC1D

subroutine PARTITION1D(A, marker)
  real(KIND=DP), intent(in out), dimension(:) :: A
  integer, intent(out) :: marker
  integer :: i, j, ilen
  real(KIND=DP):: temp
  real(KIND=DP) :: x      ! pivot point
  temp = 0._DP
  x = A(1)
  i= 0
  j= size(A) + 1

  do
     j = j-1
     do
        if (A(j) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i) and A(j)
        temp = A(i)
        A(i) = A(j)
        A(j) = temp
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine PARTITION1D

recursive subroutine QSORTC2D(A,B)
  real(KIND=DP), intent(in out), dimension(:,:) :: A
  integer(KIND=4), intent(in out), dimension(:,:) :: B
  integer :: iq

  if(size(A,1) > 1) then
     call PARTITION2D(A,B, iq)
     call QSORTC2D(A(:iq-1,:),B(:iq-1,:))
     call QSORTC2D(A(iq:,:),B(iq:,:))
  endif
end subroutine QSORTC2D

subroutine PARTITION2D(A,B, marker)
  real(KIND=DP), intent(in out), dimension(:,:) :: A
  integer(KIND=4), intent(in out), dimension(:,:) :: B
  integer, intent(out) :: marker
  integer :: i, j, ilen
  real(KIND=DP),DIMENSION(:),ALLOCATABLE :: temp
  integer,DIMENSION(:),ALLOCATABLE :: temp2
  real(KIND=DP) :: x      ! pivot point
  ilen = size(A,2)
  ALLOCATE(temp(ilen),temp2(ilen))
  temp = 0._DP
  temp2 = 0
  x = A(1,1)
  i= 0
  j= size(A,1) + 1

  do
     j = j-1
     do
        if (A(j,1) <= x) exit
        j = j-1
     end do
     i = i+1
     do
        if (A(i,1) >= x) exit
        i = i+1
     end do
     if (i < j) then
        ! exchange A(i,:) and A(j,:)
        temp = A(i,:)
        A(i,:) = A(j,:)
        A(j,:) = temp
        temp2 = B(i,:)
        B(i,:) = B(j,:)
        B(j,:) = temp2
     elseif (i == j) then
        marker = i+1
        return
     else
        marker = i
        return
     endif
  end do

end subroutine PARTITION2D

end module qsort_c_module

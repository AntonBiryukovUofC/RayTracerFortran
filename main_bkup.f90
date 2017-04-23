module raymod
   implicit none
   private
   public :: GETSR
   contains
    subroutine GETSR(src,st,src_offset,src_depth)

        use csv_module
        use iso_fortran_env, only: wp => real64

        type(csv_file) :: feq,fst
        character(len=30),dimension(:),allocatable :: header
        real(wp),dimension(:),allocatable :: xEq,yEq,zEq,xSt,ySt,zSt
        !real(wp), DIMENSION(:, :), ALLOCATABLE :: src,st

        real(wp),allocatable, intent(out) :: src(:, :),st(:, :)
        real(wp),allocatable,dimension(:), intent(out) :: src_offset,src_depth

        logical :: status_ok
        integer,dimension(:),allocatable :: itypes
        integer i,j,k,Neq,Nst,AllocateStatus
        ! Read the station file
        ! read the file
        call fst%read('stdf.csv',header_row=1,status_ok=status_ok)

        ! get the header and type info
        call fst%get_header(header,status_ok)
        call fst%variable_types(itypes,status_ok)

        ! get some data
        call fst%get(1,xSt,status_ok)
        call fst%get(2,ySt,status_ok)
        call fst%get(3,zSt,status_ok)

        ! destroy the file
        call fst%destroy()

        ! Read the sources file.

        ! read the file
        call feq%read('eqdf.csv',header_row=1,status_ok=status_ok)

        ! get the header and type info
        call feq%get_header(header,status_ok)
        call feq%variable_types(itypes,status_ok)

        ! get some data
        call feq%get(1,xEq,status_ok)
        call feq%get(2,yEq,status_ok)
        call feq%get(3,zEq,status_ok)
        Neq=size(xEq)
        Nst=size(xSt)
        ! destroy the file

        ALLOCATE ( src(Neq*Nst, 3), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        ALLOCATE ( st(Neq*Nst, 3), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        ALLOCATE ( src_offset(Neq*Nst), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        ALLOCATE ( src_depth(Neq*Nst), STAT = AllocateStatus)
        IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

        call feq%destroy()
        k=1
        do i=1,Neq,1
            do j=1,Nst,1
            src(k,1)=xEq(i)
            src(k,2)=yEq(i)
            src(k,3)=zEq(i)

            st(k,1)=xSt(j)
            st(k,2)=xSt(j)
            st(k,3)=xSt(j)

            src_offset(k) = sqrt( (st(k,1) - src(k,1))**2 + (st(k,2) - src(k,2))**2 )
            src_depth(k) = src(k,3)

            k=k+1
            end do
        end do
        print *,k
        ! By now, all pairs of source-receivers are flattened and sit in SRC, and ST.
    end subroutine


end module raymod




program RAYTRC
    use iso_fortran_env, only: wp => real64
    use raymod, only: GETSR
    real(wp), DIMENSION(:, :), ALLOCATABLE :: src,st

    real(wp), DIMENSION(:), allocatable :: src_offset,src_depth

    ! Here is the model description:
    real(wp), DIMENSION(:), allocatable :: vels,depths
    integer NLayers, k, AllocateStatus
    NLayers=2
    k=0

! Allocate the velocities and depths
    ALLOCATE ( vels(NLayers+1), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE ( depths(NLayers), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    vels=(/3100,4470,6200/)
    depths=(/2000, 4000/)


    ! Get the arrays src (Neq x 3) , and st (Nst x 3) with coordinates of the stations and the receivers:
    call GETSR(src,st,src_offset,src_depth)
    !




    print *,' offsets '
    do k=1,size(src_offset)
            print *, src_offset(k)
    end do

    print *,' depths '

    do k=1,size(src_depth)
            print *, src_depth(k)
    end do


end program



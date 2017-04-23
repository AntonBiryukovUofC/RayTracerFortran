module raymod
   use iso_fortran_env, only: wp => real64

   implicit none
   private
   public :: GETSR,whichLayer,GetPTime, InsertLayer,costFunc,costFunc_Prime,solve
   contains
    subroutine GETSR(src,st,src_offset,src_depth)

        use csv_module
        use iso_fortran_env, only: wp => real64

        type(csv_file) :: feq,fst
        character(len=30),dimension(:),allocatable :: header
        real(wp),dimension(:),allocatable :: xEq,yEq,zEq,xSt,ySt,zSt
        !real(wp), DIMENSION(:, :), ALLOCATABLE :: src,st

        real,allocatable, intent(out) :: src(:, :),st(:, :)
        real,allocatable,dimension(:), intent(out) :: src_offset,src_depth

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
        print *,Nst
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
            src(k,1)=real(xEq(i))
            src(k,2)=real(yEq(i))
            src(k,3)=real(zEq(i))

            st(k,1)=real(xSt(j))
            st(k,2)=real(ySt(j))
            st(k,3)=real(zSt(j))
            print *,'Station'
            print *,st(k,:)
            print *,'End'
            !print *,src(k,:)
            src_offset(k) = sqrt( (st(k,1) - src(k,1))**2 + (st(k,2) - src(k,2))**2 )
            src_depth(k) = src(k,3)

            k=k+1
            end do
        end do
        print *,k
        ! By now, all pairs of source-receivers are flattened and sit in SRC, and ST.
    end subroutine

    integer function whichLayer(depths,dph)
        use iso_fortran_env, only: wp => real64

        real, intent(in) :: depths(:)
        real, intent(in) :: dph
        real :: diff

        integer inN
        integer NLayers
        integer i
        i=0
        inN=0
        NLayers =size(depths)
! By adding a fake interface really deep we should guarantee this condition always gives a layer
        do i=1,NLayers
            inN=i
            diff =depths(i) - dph
            IF (diff > 0.0) exit

        end do

        whichLayer = inN
        if (diff <0) whichLayer = NLayers+1

    end function whichLayer



! This routine inserts a layer into the velocity model at the depth of the receiver

    subroutine InsertLayer(vels,depths,NLayers,vels_new,depths_new,NNew,nl,dph)

    use iso_fortran_env, only: wp => real64

    integer, intent(in) :: NLayers,NNew,nl
    real, intent(in) :: depths(NLayers)
    real, intent(in) :: vels(NLayers+1)
    real, intent(out) :: vels_new(NNew+1),depths_new(NNew)
    real, intent(in) :: dph
    if (nl<NLayers+1) then
            vels_new(1:nl)=vels(1:nl)
            !vels_new(nl+1:NNew+1) = vels(nl:NLayers+1)
            depths_new(1:nl-1)=depths(1:nl-1)
            depths_new(nl)=dph
            !depths_new(nl+1:NNew)=depths(nl:NLayers)
        else
            vels_new(1:NLayers+1)=vels
            depths_new(1:NLayers)=depths
            !vels_new(NNew)=vels(NLayers+1)
            !vels_new()=vels(NLayers+1)

            depths_new(nl)=dph
    endif
    ! Differentiate the depths here:
    depths_new(2:NNew) = depths_new(2:NNew) - depths_new(1:NNew-1)


    end subroutine




! This function actually calculates the travel time given the velocities ,depths


    real(wp) function GetPTime(src_depth,src_offset,nl_src,nl_total,vp,depths,p)

    real, intent(in) :: src_depth,src_offset
    integer, intent(in) :: nl_src,nl_total
    real, intent(in) :: depths(nl_src)
    real, intent(in) :: vp(nl_src)
    real(wp) :: timeP,p0,p_final,p
    real :: cos_t,c_harmonic,weight
    integer :: iters
    logical :: debug,P0_bad    ! set to .true. or .false.

    if (nl_src == 1) then
        print *, 'Source in the top layer'
        timeP = sqrt(src_depth**2+src_offset**2)/vp(nl_src)
    else
        ! Here we actually need to do ray-tracing
        weight=0.9
        c_harmonic = src_depth/sum(depths/vp)
        cos_t = src_depth/sqrt(src_offset**2 + src_depth**2)
        ! Comment / uncomment this line
        p0=p
        p0 = cos_t / c_harmonic*weight
        print *, c_harmonic
        print *, p0

        print *, ' '  ! blank line
		! Check if the p0 estimate is any good :
        P0_bad=.TRUE.

        do while (P0_bad)
            if (ISNAN(sum(sqrt(1-(p0**2)*vp**2)))) then
                p0=p0/2
                print *,'p0 divided by 2...'
            else
                P0_bad = .FALSE.
                print *,'p0 is good !'
            end if
        end do


        call solve( p0, p_final, iters, debug,depths,vp,src_offset)

!        print 11, p_final, iters
!11      format('solver returns x = ', e22.15, ' after', i3, ' iterations')
 !       fx = costFunc(p_final,depths,vp,src_offset)
 !       print 12, fx
!12      format('the value of f(x) is ', e22.15)



    end if



    GetPTime = timeP
    end function GetPTime




 ! Here we define the cost function to find the root for and its derivative
    real(wp) function costFunc(x,H,V,R)
        real, intent(in) :: H(:)
        real, intent(in) :: V(:)
        real, intent(in) :: R
        real(wp), intent(in) :: x
        real(wp) :: sum_term(size(H)),sum_all
        print *, 'Vi is ', V
        print *, 'Hi is ', H
        print *, 'R is ', R

        !a = 1-(x**2)*V**2
        sum_term = H*V*x/sqrt(1-(x**2)*V**2)
        print *,'SumTerm',sum_term

        sum_all =  sum(sum_term)
        print *,'SumAll',sum_all
        costFunc =  R - sum_all
        !costFunc =  0
    end function costFunc



    real(wp) function costFunc_Prime(x,H,V)
        real, intent(in) :: H(:)
        real, intent(in) :: V(:)
        real(wp), intent(in) :: x
        real(wp) :: sum_term(size(H))
        real(wp) :: denom(size(H))


        denom=sqrt((1-(x**2)*V**2)**3 )

        !if np.isnan(denom).any():
        !    raise ValueError('Denominator is nan')

        sum_term = (H*V)/denom
        costFunc_Prime = -sum(sum_term)
    end function costFunc_Prime




subroutine solve(x0, x, iters, debug,H,V,R)
    implicit none
    ! The parameters of the ray-path here
    real, intent(in) :: H(:)
    real, intent(in) :: V(:)
    real, intent(in) :: R
    ! The solver parameters here
    integer, parameter :: maxiter = 8
    real(kind=8), parameter :: tol = 1.d-2
    ! Estimate the zero of f(x) using Newton's method.
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!)
    !   the number of iterations iters
    real(kind=8), intent(in) :: x0
    logical, intent(in) :: debug
    real(kind=8), intent(out) :: x
    integer, intent(out) :: iters

    ! Declare any local variables:
    real(kind=8) :: deltax, fx, fxprime
    integer :: k

    ! Check the HVR:
   ! print *, 'Vi is ', V
   ! print *, 'Hi is ', H
   ! print *, 'R is ', R



    ! initial guess
    x = x0

    if (debug) then
        print 11, x
 11     format('Initial guess: x = ', e22.15)
        endif

    ! Newton iteration to find a zero of f(x)

    do k=1,maxiter

        ! evaluate function and its derivative:
        fx = costFunc(x,H,V,R)
        print *,'FX=',fx
        fxprime = costFunc_Prime(x,H,V)
        print *,'Fprime=',fxprime

        if (abs(fx) < tol) then
            exit  ! jump out of do loop
            endif

        ! compute Newton increment x:
        deltax = fx/fxprime

        ! update x:
        x = x - deltax

        if (debug) then
            print 12, k,x
 12         format('After', i3, ' iterations, x = ', e22.15)
            endif

        enddo


    if (k > maxiter) then
        ! might not have converged

        fx = costFunc(x,H,V,R)
        if (abs(fx) > tol) then
            print *, '*** Warning: has not yet converged'
            endif
        endif

    ! number of iterations taken:
    iters = k-1


end subroutine solve





end module raymod




program RAYTRC
    use iso_fortran_env, only: wp => real64

    use raymod, only: GETSR,whichLayer,GetPTime,InsertLayer
    implicit none

    real, DIMENSION(:, :), ALLOCATABLE :: src,st

    real, DIMENSION(:), allocatable :: src_offset,src_depth
    real :: dph
    real(wp) :: timeP,p
    ! Here is the model description:
    real, DIMENSION(:), allocatable :: vels,depths,vels_new,depths_new
    integer NLayers, k, AllocateStatus,nl
    integer NNew
    k=0
    NLayers=2

! Allocate the velocities and depths
    ALLOCATE ( vels(NLayers+1), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE ( depths(NLayers), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    ! Allocate the velocities and depths
    ALLOCATE ( vels_new(NLayers+2), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"
    ALLOCATE ( depths_new(NLayers+1), STAT = AllocateStatus)
    IF (AllocateStatus /= 0) STOP "*** Not enough memory ***"

    vels=(/3100,3270,5000/)
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

    print *,' depths '

    do k=1,size(src_depth)
            dph = src_depth(k)
            nl = whichLayer(depths,dph)
            if (nl == 1) then
                    !print *,'Source in ',nl, ' layer'
                    p=0.011
                    timeP = GetPTime(dph,src_offset(k),nl,NLayers,vels,depths,p)
                else
                    ! Here we need to pick Hi,Vi by introducing a new layer
                    NNew = NLayers+1
                    p=1e-5
                    call InsertLayer(vels,depths,NLayers,vels_new,depths_new,NNew,nl,dph)
                    print *, 'Source in ', nl, ' layer'

                    print *, 'Vels is ', vels_new
                    print *, 'Depths is ', depths_new

                    timeP = GetPTime(dph,src_offset(k),nl,NLayers,vels_new(1:nl),depths_new(1:nl),p)

            endif
            !timeP=GetPTime(dph,src_offset(k),nl,NLayers,vels,depths)
            print *, 'Time is ', timeP



    end do


end program



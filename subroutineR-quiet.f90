module raymod

   implicit none
   private
   public :: whichLayer,GetPTime, InsertLayer,costFunc,costFunc_Prime,solve, dofullforwardproblem,solvebst,TraceRays
   contains


    integer function whichLayer(depths,dph)

        double precision, intent(in) :: depths(:)
        double precision, intent(in) :: dph
        double precision :: diff

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

    !integer AllocateStatus
    integer, intent(in) :: NLayers,NNew,nl
    double precision, intent(in) :: depths(NLayers)
    double precision, intent(in) :: vels(NLayers+1)
    double precision, intent(in) :: dph
    double precision, DIMENSION(:)::  vels_new,depths_new



    !print *,vels,depths,NLayers,vels_new,depths_new,NNew,nl,dph
  !  print *,vels,depths
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


    double precision function GetPTime(src_depth,src_offset,nl_src,nl_total,vp,depths,p,keep_delta)

    double precision, intent(in) :: src_depth,src_offset
    integer, intent(in) :: nl_src,nl_total
    double precision, intent(in) :: depths(nl_src)
    double precision, intent(in) :: vp(nl_src)
    double precision :: cosV(nl_src),t_int(nl_src), d(nl_src)
    integer, intent(in) :: keep_delta
    double precision :: timeP,p0,p_final,p,cf,cfp,check_x
    double precision :: cos_t,c_harmonic,weight=1
    integer :: iters
    logical :: debug,P0_bad,conv    ! set to .true. or .false.
   ! print *, 'Vsegments'
   ! print *,vp
   ! print *, 'Depth Segments'
   ! print *,depths
   ! print *, 'Calculating T'
    if (nl_src == 1) then
       ! print *, 'Source in the top layer'
        timeP = sqrt(src_depth**2+src_offset**2)/vp(nl_src)
        conv=.TRUE.
    else
        ! Here we actually need to do ray-tracing
        conv=.FALSE.

        c_harmonic = src_depth/sum(depths/vp)
        cos_t = src_depth/sqrt(src_offset**2 + src_depth**2)
        ! Comment / uncomment this line - division by 2 is weird here
        !p0=p
        p0 = cos_t / c_harmonic*weight
        !p0 = 1/maxval(vp)-1d-10
       ! print *, c_harmonic
       ! print *, p0

       ! print *, ' '  ! blank line
        ! Check if the p0 estimate is any good :
        P0_bad=.TRUE.

        do while (P0_bad)
            if (ISNAN(sum(sqrt(1-(p0**2)*(vp+1)**2)))) then
                p0=p0/2
               ! print *,'p0 divided by 2...'
            else
                P0_bad = .FALSE.
               ! print *,'p0 is good !'
            end if
        end do

        ! Do a few safety checks - first, if the function is already negative-valued, then Newton is appropriate:
        cf = costFunc(p0,depths,vp,src_offset)
        cfp =  costFunc_Prime(p0,depths,vp)
        check_x = p0 - cf/cfp
        if (cf<0) then
           ! print *,'Cost function check: negative; Doing Newton..'
            call solve( p0, p_final, iters, debug,depths,vp,src_offset,conv)
        elseif (check_x<1/maxval(vp)-(1d-10)) then
           ! print *,'First x jump is safe; Doing Newton..'
            call solve( p0, p_final, iters, debug,depths,vp,src_offset,conv)
        else
            !print *,'Checks failed to pass...Will do bisect for a little bit, until the jump is safe'
            ! Basically, that means there will be negative terms under sqrt
            call solvebst( 1d-10,1/maxval(vp)-1d-12,p_final, iters, debug,depths,vp,src_offset)
            !print *,'Bisection moved initial point to new point',p_final
           ! print *,'Bisect took iterations:', iters
            p0=p_final
            call solve( p0, p_final, iters, debug,depths,vp,src_offset,conv)
            endif
       ! print *,conv
        ! Here we use p_final and calculate the travel times:
        cosV = sqrt(1-(p_final**2)*vp**2)
        if (keep_delta > 0) then
           d = sqrt((depths/cosV)**2 - depths**2)
           open(unit=1,file='rays.dat',form="FORMATTED",status='OLD',action='READWRITE',position='append')
           write(unit=1,FMT=*) d
           write(unit=1,FMT=*) depths
           !write(unit=1,FMT=*) ""
           close(1)
        end if
        t_int = depths/(vp*cosV)
        timeP = sum(t_int)
        if (.not.(conv)) then
            timeP = -999
        end if


    end if



    GetPTime = timeP

    end function GetPTime




 ! Here we define the cost function to find the root for and its derivative
    double precision function costFunc(x,H,V,R)
        double precision, intent(in) :: H(:)
        double precision, intent(in) :: V(:)
        double precision, intent(in) :: R
        double precision, intent(in) :: x
        double precision :: sum_term(size(H)),sum_all
        !print *, 'Vi is ', V
        !print *, 'Hi is ', H
        !print *, 'R is ', R

        !a = 1-(x**2)*V**2
        sum_term = H*V*x/sqrt(1-(x**2)*V**2)
      !  print *,'SumTerm',sum_term

        sum_all =  sum(sum_term)
       ! print *,'SumAll',sum_all
        costFunc =  R - sum_all
        !costFunc =  0
    end function costFunc



    double precision function costFunc_Prime(x,H,V)
        double precision, intent(in) :: H(:)
        double precision, intent(in) :: V(:)
        double precision, intent(in) :: x
        double precision :: sum_term(size(H))
        double precision :: denom(size(H))


        denom=sqrt((1-(x**2)*V**2))**3

        !if np.isnan(denom).any():
        !    raise ValueError('Denominator is nan')

        sum_term = (H*V)/denom
        costFunc_Prime = -sum(sum_term)
    end function costFunc_Prime




subroutine solve(x0, x, iters, debug,H,V,R,conv)
    implicit none
    ! The parameters of the ray-path here
    double precision, intent(in) :: H(:)
    double precision, intent(in) :: V(:)
    double precision, intent(in) :: R
    ! The solver parameters here
    integer, parameter :: maxiter = 15
    double precision, parameter :: tol = 1.d-1
    ! Estimate the zero of f(x) using Newton's method.
    ! Input:
    !   f:  the function to find a root of
    !   fp: function returning the derivative f'
    !   x0: the initial guess
    !   debug: logical, prints iterations if debug=.true.
    ! Returns:
    !   the estimate x satisfying f(x)=0 (assumes Newton converged!)
    !   the number of iterations iters
    double precision, intent(in) :: x0
    logical, intent(in) :: debug
    double precision, intent(out) :: x
    integer, intent(out) :: iters
    logical,intent(out) :: conv
    ! Declare any local variables:
    double precision :: deltax=0, fx=0, fxprime=0,weight=0,rr=1d-9
    integer :: k

    ! Check the HVR:
   ! print *, 'Vi is ', V
   ! print *, 'Hi is ', H
   ! print *, 'R is ', R
   ! 
    conv=.FALSE.
    ! initial guess
    x = x0

  !  if (debug) then
  !      print 11, x
 !11     format('Initial guess: x = ', e22.15)
 !       endif

    ! Newton iteration to find a zero of f(x)

    do k=1,maxiter
       ! print *, 'Vsegments - Cost'
       ! print *,V
       ! print *, 'Depth Segments - Cost'
       ! print *,H
        ! evaluate function and its derivative:
        fx = costFunc(x,H,V,R)
       ! print *,'FX=',fx
      !  if (ISNAN(fx)) then
      !      weight = 0.5
      !      print *,'Modifying x by division'
      !      x = x * weight
      !      cycle
      !  end if

        fxprime = costFunc_Prime(x,H,V)
        !print *,'Fprime=',fxprime

        if (abs(fx) < tol) then
            conv=.TRUE.
            exit  ! jump out of do loop

            endif

        ! compute Newton increment x:
        deltax = fx/fxprime
        !print *,'DX=',deltax

        ! update x:
        x = x - deltax
        if (x > 1/maxval(V)) then
          !  print *,'Modifying x by shift to boundary'

            x = 1/maxval(V) - rr

        end if

        !if (debug) then
         !   print 12, k,x
 !12      !   format('After', i3, ' iterations, x = ', e22.15)
         !   endif

        enddo


    if (k > maxiter) then
        ! might not have converged

        fx = costFunc(x,H,V,R)


       ! if (abs(fx) > tol) then
          !  print *, '*** Warning: has not yet converged'
           ! endif
        endif

    ! number of iterations taken:
    iters = k-1
    if (abs(fx) > tol) then
        conv=.TRUE.
       ! print *,'Converged'
    end if

end subroutine solve





! Bisection search for the zero; to help with the large Newton Jumps:
subroutine solvebst(x1,x2,x,iters, debug,H,V,R)
    implicit none
    ! The parameters of the ray-path here
    double precision, intent(in) :: H(:)
    double precision, intent(in) :: V(:)
    double precision, intent(in) :: R
    ! The solver parameters here
    integer, parameter :: maxiter = 20
    double precision, intent(in) :: x1,x2
    logical, intent(in) :: debug
    double precision, intent(out) :: x
    integer, intent(out) :: iters
    double precision :: dx=0, fx=0, fxprime=0,fmid=0,xmid=0,f=0,check_x=0,deltax=0,tol=1d-1
    integer :: k
    fmid = costFunc(x2,H,V,R)
    f = costFunc(x1,H,V,R)
   ! if (f*fmid >=0) then
     !   print *,'Function is not changing sign!!!'
   ! end if
    ! orienting the search so that f >0 lies at x + dx
    if (f<0) then
        x=x1
        dx = x2-x1
    else
        x=x2
        dx=x1-x2
    endif

    do k=1,maxiter
        dx=dx*0.5
        xmid = x+dx
        fmid = costFunc(xmid,H,V,R)
        if (fmid <0) then
        x=xmid
        endif

        if (fmid == 0) then
        exit
        endif
        fx = costFunc(xmid,H,V,R)
        !print *,'Bisect FX=',fx
        fxprime = costFunc_Prime(xmid,H,V)
       ! print *,'Bisect Fprime=',fxprime

        ! compute Newton increment x:
        deltax = fx/fxprime
        !print *,'DX=',deltax
        check_x = x - deltax
        ! if the jump is safe:
        if(check_x<1/maxval(V)-(1d-10)) then
          !  print *,'Bisect: X jump is safe; returning X and doing Newton..'
            x=xmid
            exit
        endif

        if (abs(fmid) < tol) then
            exit  ! jump out of do loop
        endif
        enddo



    ! number of iterations taken:
    iters = k-1


end subroutine solvebst


subroutine dofullforwardproblem(vels,depths,NLayers,src_offset,src_depth,NSrc,timeP,keep_delta) bind(C, name="dff_")


    use, intrinsic :: iso_c_binding, only : c_double, c_int

    integer(c_int),intent(in) :: NLayers,NSrc
    integer ::  k,nl
    integer(c_int), intent(in) :: keep_delta
    real(c_double), DIMENSION(NSrc):: src_offset,src_depth
    double precision :: dph
    double precision :: p
    real(c_double), DIMENSION(NSrc),intent(out) :: timeP
    ! Here is the model description:
    real(c_double), DIMENSION(NLayers+1),intent(in):: vels
    real(c_double), DIMENSION(NLayers),intent(in):: depths
    double precision, DIMENSION(NLayers+2):: vels_new
    double precision, DIMENSION(NLayers+1):: depths_new
    integer NNew,AllocateStatus

    depths_new=0
    vels_new=0
    k=0
    timeP=-1
!    print *,' depths '
    open(unit=1,file='rays.dat',form="FORMATTED",status='REPLACE',action='READWRITE')
    close(1)
    do k=1,size(src_depth)
            dph = src_depth(k)
            nl = whichLayer(depths,dph)

            if (nl == 1) then
                   ! print *,'Source in ',nl, ' layer'
                    p=0.011
                    timeP(k) = GetPTime(dph,src_offset(k),nl,NLayers,vels,depths,p,keep_delta)

                else
                    ! Here we need to pick Hi,Vi by introducing a new layer
                    NNew = NLayers+1
                    p=1e-5
                   ! print *, 'Vels is ', vels
                  !  print *, 'Depths is ', depths
                  !  print *,NLayers,NNew
                    call InsertLayer(vels,depths,NLayers,vels_new,depths_new,NNew,nl,dph)
                   ! print *, 'Source in ', nl, ' layer'
                   ! print *, 'Vels is ', vels_new
                    !print *, 'Depths is ', depths_new
                   ! print *, 'Calculating time'
                    timeP(k) = GetPTime(dph,src_offset(k),nl,NLayers,vels_new(1:nl),depths_new(1:nl),p,keep_delta)
            endif
            !timeP=GetPTime(dph,src_offset(k),nl,NLayers,vels,depths)
           ! print *, 'Time is ', timeP(k)
           ! print *, '--------------------------------------'
           ! print *, '--------------------------------------'
           ! print *, '--------------------------------------'
    end do
end subroutine



subroutine TraceRays(vels,depths,NLayers,src_offset,src_depth,NSrc,timeP,keep_delta)


    integer,intent(in) :: NLayers,NSrc
    integer ::  k,nl
    integer, intent(in) :: keep_delta
    double precision, DIMENSION(NSrc):: src_offset,src_depth
    double precision :: dph
    double precision :: p
    double precision, DIMENSION(NSrc),intent(out) :: timeP
    ! Here is the model description:
    double precision, DIMENSION(NLayers+1),intent(in):: vels
    double precision, DIMENSION(NLayers),intent(in):: depths
    double precision, DIMENSION(NLayers+2):: vels_new
    double precision, DIMENSION(NLayers+1):: depths_new
    integer NNew,AllocateStatus

    depths_new=0
    vels_new=0
    k=0
    timeP=-1
!    print *,' depths '
    open(unit=1,file='rays.dat',form="FORMATTED",status='REPLACE',action='READWRITE')
    close(1)
    do k=1,size(src_depth)
            dph = src_depth(k)
            nl = whichLayer(depths,dph)

            if (nl == 1) then
                   ! print *,'Source in ',nl, ' layer'
                    p=0.011
                    timeP(k) = GetPTime(dph,src_offset(k),nl,NLayers,vels,depths,p,keep_delta)

                else
                    ! Here we need to pick Hi,Vi by introducing a new layer
                    NNew = NLayers+1
                    p=1e-5
                   ! print *, 'Vels is ', vels
                  !  print *, 'Depths is ', depths
                  !  print *,NLayers,NNew
                    call InsertLayer(vels,depths,NLayers,vels_new,depths_new,NNew,nl,dph)
                   ! print *, 'Source in ', nl, ' layer'
                   ! print *, 'Vels is ', vels_new
                    !print *, 'Depths is ', depths_new
                   ! print *, 'Calculating time'
                    timeP(k) = GetPTime(dph,src_offset(k),nl,NLayers,vels_new(1:nl),depths_new(1:nl),p,keep_delta)
            endif
            !timeP=GetPTime(dph,src_offset(k),nl,NLayers,vels,depths)
           ! print *, 'Time is ', timeP(k)
           ! print *, '--------------------------------------'
           ! print *, '--------------------------------------'
           ! print *, '--------------------------------------'
    end do
end subroutine


end module raymod


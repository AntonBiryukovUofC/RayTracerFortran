subroutine makevoro2c(voro,npt,V0d,nptref,voro2c,V0d2c,npt2c)
implicit none
include 'params.h'

real, intent(in) :: voro(malay,4)
real, intent(in) :: V0d(maxlay,2)
integer, intent(in) :: nptref,npt

integer, intent(out) :: npt2c
real, intent(out) :: voro2c(maxlay,4), V0d2c(maxlay)

integer num,j,i,order(npt+nptref-2),ind
real pos(npt+nptref-2),maxx,minn,temp(4)



num = nptref+npt-2


pos=[voro(2:npt,1), V0d(2:nptref,1)]

!%%%%%%%%%%%%%%%%%

ind=1
!!!dismin=99999999999999
!print *, npt
do i=1,num
	!write(*,*)'ordre',i
	maxx=99999999
	if (i==1) then
		minn=-9999999
	else
		minn=pos(order(i-1))
	endif
	
	!write(*,*)'min','max',minn,maxx
	do j=1,num
		if ((minn<pos(j)).and.(pos(j)<=maxx)) then
	!		write(*,*)j,minn,voro(j,1),maxx
			ind=j
			maxx=pos(j)
	!		write(*,*)ind,maxx
		endif
		
! ! 		if ((abs(voro(j,1)-voro(i,1))<dismin).and.(i.ne.j)) then
! ! 		dismin=abs(voro(j,1)-voro(i,1))
! ! 		endif
		
	enddo
	order(i)=ind
 	!write(*,*)ind
 	!write(*,*)
!   	write(*,*)

enddo

!5555555555555555555555555555555555555555555555555555555
voro2c(1,1)=0
voro2c(1,2)=V0d(1,2)*(1+voro(1,2))
voro2c(1,3)=voro(1,3)
voro2c(1,4)=voro(1,4) !!CHANGE!!
V0d2c(1)=V0d(1,2)

temp(1)=V0d(1,2)
temp(2)=voro(1,2)
temp(3)=voro(1,3)
temp(4)=voro(1,4)  !!CHANGE!!

      
do i=1,num
      !write(*,*)i,num
      if (order(i)<npt) then ! it is a dv/v
		
		ind=order(i)+1
		!write(*,*)'1',ind
		if ((ind<2).or.(ind>npt)) then
			stop "stop2"
		endif

		voro2c(i+1,1)=voro(ind,1)
		voro2c(i+1,2)=temp(1)*(1+voro(ind,2))
		voro2c(i+1,3)=voro(ind,3)
		voro2c(i+1,4)=voro(ind,4)  !!CHANGE!!
		V0d2c(i+1)=temp(1)
		
		temp(2)=voro(ind,2)
		temp(3)=voro(ind,3)
		temp(4)=voro(ind,4) !!CHANGE!!

      elseif ((order(i)>=npt).and.(order(i)<=num)) then ! it is a V0ref
		
		ind=order(i)-npt+2
		!write(*,*)ind,order(i),npt
		if ((ind<2).or.(ind>nptref)) then
			  stop "stop1"
		endif
	      
		voro2c(i+1,1)=V0d(ind,1)
		voro2c(i+1,2)=V0d(ind,2)*(1+temp(2))
		voro2c(i+1,3)=temp(3)
		voro2c(i+1,4)=temp(4) !!CHANGE!!
		V0d2c(i+1)=V0d(ind,2)
		
		temp(1)=V0d(ind,2)	
      else
      
      write(*,*) 'stop2'
      stop
      endif 
      
enddo

npt2c=num+1
! ! ! ! ! do i=1,npt2c-1
! ! ! ! !       if (voro2c(i,1)==voro2c(i+1,1)) then 
! ! ! ! !       
! ! ! ! !      
! ! ! ! ! 	    write(*,*)'----Y--Y--Y--Y--Y--Y-------------'
! ! ! ! ! 	    write(*,*)nptref
! ! ! ! ! 	    do j=1,nptref
! ! ! ! ! 	    write(*,*)V0d(j,1),V0d(j,2)
! ! ! ! ! 	    enddo
! ! ! ! ! 	    write(*,*)'---------------------------'
! ! ! ! ! 	    write(*,*)'npt',npt
! ! ! ! ! 	    do j=1,npt
! ! ! ! ! 	    write(*,*)voro(j,1),voro(j,2)
! ! ! ! ! 	    enddo
! ! ! ! ! 	    
! ! ! ! ! ! ! 	    write(*,*)
! ! ! ! ! ! ! 	    write(*,*)'---------------------------'
! ! ! ! ! ! ! 	    do j=1,npt2c
! ! ! ! ! ! !             write(*,*)voro2c(j,1),voro2c(j,2),V0d2c(j)
! ! ! ! ! ! !             enddo
! ! ! ! ! ! !             write(*,*)'---------------------------'
! ! ! ! ! ! ! 	    stop "same depth"
! ! ! ! !       endif
! ! ! ! ! enddo

! ! ! if (voro2c(12,2)>50) then
! ! ! write(*,*)'----Y--Y--Y--Y--Y--Y-------------'
! ! ! 	    write(*,*)nptref
! ! ! 	    do j=1,nptref
! ! ! 	    write(*,*)V0d(j,1),V0d(j,2)
! ! ! 	    enddo
! ! ! 	    write(*,*)'---------------------------'
! ! ! 	    write(*,*)'npt',npt
! ! ! 	    do j=1,npt
! ! ! 	    write(*,*)voro(j,1),voro(j,2)
! ! ! 	    enddo
! ! ! 
! ! ! 	    write(*,*)'---------------------------'
! ! ! 	    do j=1,npt2c
! ! !             write(*,*)voro2c(j,1),voro2c(j,2),V0d2c(j)
! ! !             enddo
! ! ! endif

return
end

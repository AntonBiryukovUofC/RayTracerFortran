subroutine dis2qmodel(voro2c,V0d2c,npt2c,d_min,d_max,h,beta,beta0,xi,trend,dismin)
implicit none
include 'params.h'

real, intent(in) :: voro2c(maxlay,4),V0d2c(maxlay),d_min,d_max
real, intent(out) :: beta(maxlay), beta0(maxlay),h(maxlay)
integer, intent(in) :: npt2c

real, intent(out) :: xi(maxlay),trend(maxlay),dismin

real maxx,minn,summ
integer ind,i,j,order(npt2c),k
!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
beta=0
beta0=0
h=0
ind=1
dismin=d_max*1000
do i=1,npt2c

	maxx=d_max*1000
	if (i==1) then
		minn=d_min-1
	else
		minn=voro2c(order(i-1),1)
	endif

	do j=1,npt2c
		if ((minn<voro2c(j,1)).and.(voro2c(j,1)<=maxx)) then
			ind=j
			maxx=voro2c(j,1)
		endif
		
		if ((abs(voro2c(j,1)-voro2c(i,1))<dismin).and.(i.ne.j)) then
		dismin=abs(voro2c(j,1)-voro2c(i,1))
		endif
		
	enddo
	order(i)=ind
enddo
!----------------------------------------------------------------------------

if (order(1)/=1) then
	  write(*,*)'Problem'
	  write(*,*)
	  do i=1,npt2c
		write(*,*)i,voro2c(i,1)
	  enddo
	  write(*,*)
	  do i=1,npt2c
		write(*,*)order(i),voro2c(order(i),1)
	  enddo
	  stop
endif
!----------------------------------------------------------------------

summ=0
do i=2,npt2c
		h(i-1)= voro2c(order(i),1) - summ
		summ = summ + h(i-1)
		beta(i-1) = voro2c(order(i-1),2)
		beta0(i-1) = V0d2c(order(i-1))
		xi(i-1)=voro2c(order(i-1),3)
		trend(i-1)=voro2c(order(i-1),4)!!CHANGE!!
enddo

h(npt2c)=0
beta(npt2c) = voro2c(order(npt2c),2)
beta0(npt2c) = V0d2c(order(npt2c))
xi(npt2c)=  voro2c(order(npt2c),3)
trend(npt2c)=  voro2c(order(npt2c),4)!!CHANGE!!

return
end

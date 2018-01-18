subroutine misfitCC(traces_o1,traces_o2,ntr,nb, dt, npt, beta, h, pct_b, trend, plunge, baz, slow,&
	Ar,data11,data21,data13,data23,&
	LSrt1,like1,LSrt3,like3,Ene1,Ene2,zerotr,dmax)
	
! liker3 is for P2S receiver functions
! liker1 is for S2P receiver functions

implicit none
include 'params.h'
!include 'data_params.h'

!real, parameter ::    pi = 3.14159265      !don't touch
!real, parameter ::    rad = 0.017453292    !don't touch
!integer, parameter :: nb2 = nb/2+1 
integer nb,nb2
integer npt,ntr,j,zerotr
real traces_o1(3,maxsamp,maxtr), traces_o2(3,maxsamp,maxtr),dt,baz(maxtr),slow(maxtr),tracesP_e(3,maxsamp,maxtr),&
tracesSV_e(3,maxsamp,maxtr),NS_e(nb,ntr),EW_e(nb,ntr),NS_o(nb,ntr),EW_o(nb,ntr)
real beta(maxlay),h(maxlay),vs(maxlay),thick(maxlay),trendd(maxlay),&
fs,angle,gauss_a,vt(nb),ht(nb),pp,din, plungee(maxlay),&
alpha(maxlay),pct_a(maxlay),pct_b(maxlay),rho(maxlay),&
plunge(maxlay),strike(maxlay),dip(maxlay),LSrt1,LSrt3

double precision wpic(nb),x(nb),Ar,br

real like1,LSr,sumEW(ntr),sumNS(ntr),lambda(maxtr),like3,dmax

real dat(2*nb),rs(npt),va,dato(2*nb),Ene1(maxtr),Ene2(maxtr) 
real frini,frint,fr,suma, trend(maxlay),thicktot
integer k,fk,i,indx,itr,nptCC

real rin,wf(nb),t,v_est(nb),h_est(nb),w(nb/2+1),wi(nb/2+1),&
ww(nb),wmax,fai(nb),rft(nb),rvmax,wdata1(nb,ntr),wdata2(nb,ntr),&
data11(maxsamp,maxtr),data21(maxsamp,maxtr),&
data13(maxsamp,maxtr),data23(maxsamp,maxtr)

!```````````````
complex cs,yi,rff(nb),fNS_e(nb,ntr),fEW_e(nb,ntr),& 
fwdata1(nb,ntr),fwdata2(nb,ntr),fNS_o(nb,ntr),fEW_o(nb,ntr)


!write(*,*)'HHHHHHHHHHH'
!write(*,*) 

!*************************************************************************
!    COMPUTE THE 
!*************************************************************************
nb2 = nb/2+1 
fs = 1/dt

frini=0.
frint=fs/nb
yi=(0,1)

zerotr=0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!WRITE THE FULL MODEL (density Vp, etc ...) @ metre ds .inc
   !write(*,*)'npt2c',npt
	do i=1,npt
	thick(i)=1000*h(i)  
	vs(i) = 1000*beta(i)
        alpha(i) = vpvs*vs(i)
        pct_a(i) = 0.5 * pct_b(i)
        rho(i)   = 1000*(2.35+0.036*(alpha(i)/1000-3.0)**2)
        trendd(i) = trend(i) * pi /180
        plungee(i) = plunge(i) * pi /180
        strike(i) =  0
        dip(i) =  0
	enddo

! ! do itr=1,ntr
! ! baz(itr)=baz(itr)* pi /180  
! ! enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! VECTOR OF FREQUENCIES

do k=1,nb2
        fk=k-1
        fr=frini+frint*fk
        w(k)=2*pi*fr
enddo

wi=w

!!!!!!!!!!!!!!!!! Compute the SEISMOGRAMS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!print *,  'strike1'
!call seis_spread(npt,ntr,nb,dt,baz,slow,h,beta,alpha,isoflag,pct_a,pct_b,rho,trend,plunge,strike,dip,traces_e)

!!!!!call xi2pctb(vs,rho,xi,pct_a,pct_b)
!print *,  'strike2'

! write(*,*)'--------------'
! do i=1,npt
! write(*,*)thick(i),vs(i),pct_b(i),trendd(i),plungee(i)
! enddo

! write(*,*)thick
nptCC=npt
thicktot=0
do i=1,npt
	thicktot=thick(i)+thicktot
	if(thicktot.ge.(dmax*1000)) then
	      nptCC=i
	      thick(i)=0
	      exit
	endif
enddo
! write(*,*)'---------------'
! write(*,*)thicktot,npt,nptCC

call forward_anirec(nptCC,ntr,nb,dt,baz,slow,thick,vs,alpha,&
pct_a,pct_b,rho,trendd,plungee,strike,dip,tracesP_e,tracesSV_e)

! write(*,*)
! write(*,*)traces_o(1,218:223,12)!,traces_o(1,218:223,12)
! write(*,*)
! write(*,*)traces_e(1,200:305,1)
!print *,  'strike2'

!##################################################
!##################################################
!##################################################
!                               FOR S RECEIVER FUNCTIONS
!##################################################


!!!!!!!!!!!!!!!!! Normalize to Unit Energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do itr = 1,ntr
      Ene1(itr)=0
      do k=1,nb
	      do j=1,3
		      Ene1(itr)=Ene1(itr)+ tracesSV_e(j,k,itr)**2
	      enddo
      enddo
      Ene1(itr)=sqrt(Ene1(itr))
      if (Ene1(itr)==0) then
	    zerotr=1
      endif
enddo

do itr = 1,ntr
do k=1,nb
      NS_e(k,itr) = tracesSV_e(1,k,itr)/Ene1(itr)
      EW_e(k,itr) = tracesSV_e(3,k,itr)/Ene1(itr)
      NS_o(k,itr) = traces_o2(1,k,itr)
      EW_o(k,itr) = traces_o2(3,k,itr)
      enddo
enddo
! write(*,*)NS_e(25:27,1)
! write(*,*)NS_o(25:27,1)

!  print *,  EW_o(662:670,1)
!  print *, 'a'
  !print *,  EW_e(662:670,1)

sumNS=0
do itr = 1,ntr
do k=1,nb ! We sum over nb !!
	sumNS(itr)=sumNS(itr)+(NS_e(k,itr))**2
enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Only take the first ndatar in dat(2*i-1)/nb and pass in the frequency domain

!vshort=v
do itr = 1, ntr

    do k=1,nb
	dato(2*k-1) = NS_e(k,itr)
	dato(2*k)=0
	!if (i>460) dato(2*i-1) = 0 !!! TAPERING !!!!!
    enddo

    call four1(dato,nb,1)
    do k=1,nb
	fNS_e(k,itr)=cmplx(dato(2*k-1),dato(2*k))
    enddo

enddo

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!! Horizontal ***********************************************

sumEW=0
do itr = 1,ntr
    do k=1,nb ! We sum over nb !!
	    sumEW(itr)=sumEW(itr)+(EW_e(k,itr))**2
    enddo
enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Only take the first ndatar in dat(2*i-1)/nb and pass in the frequency domain

do itr = 1, ntr

    do k=1,nb
	dato(2*k-1) = EW_e(k,itr)
	dato(2*k)=0
	!if (i>460) dato(2*i-1) = 0 !!! TAPERING !!!!!
    enddo

    call four1(dato,nb,1)
    do k=1,nb
	fEW_e(k,itr)=cmplx(dato(2*k-1),dato(2*k))
    enddo

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sumh=sqrt(sumh)
lambda=0
do itr=1,ntr
    lambda(itr)=sqrt(sumNS(itr)+sumEW(itr))
enddo
 

 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !*************************************************************
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!!!!!!!!!!!! PASS THE OBSERVATIONS IN THE FREQUENCY DOMAIN !!!!!!!!!
! 
! 
! !!!!!!!! NS !!!!!!!!!!!!!!!!!!!
do itr = 1, ntr
      do k=1,nb
	      dato(2*k-1)=NS_o(k,itr)
	      dato(2*k)=0
      enddo

      call four1(dato,nb,1)
      do k=1,nb
	  fNS_o(k,itr)=cmplx(dato(2*k-1),dato(2*k))
      enddo
enddo

!!!!!!!! EW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr = 1, ntr
      do k=1,nb
	      dato(2*k-1)=EW_o(k,itr)
	      dato(2*k)=0
      enddo

      call four1(dato,nb,1)
      do k=1,nb
	  fEW_o(k,itr)=cmplx(dato(2*k-1),dato(2*k))
      enddo
enddo

! !******************************************************************
! !!!  COMPUTE THE CROSS PRODUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !******************************************************************
! 
! !!!!!!!!!!! NS_o * EW_e  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr=1,ntr
do k=1,nb
    fwdata1(k,itr)=fNS_o(k,itr)*fEW_e(k,itr)
enddo
enddo

!Pass in time
do itr=1,ntr

do i=1,nb
	dato(2*i-1)=real(fwdata1(i,itr))
	dato(2*i)=aimag(fwdata1(i,itr))
enddo

call four1(dato,nb,-1)
do i=1,nb
	wdata1(i,itr)=dato(2*i-1)/nb
enddo

enddo

! !!!!!!!!!!!  EW_o * NS_e !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr=1,ntr
do i=1,nb
    fwdata2(i,itr)=fEW_o(i,itr)*fNS_e(i,itr)
enddo
enddo

!Pass in time
do itr=1,ntr

do i=1,nb
	dato(2*i-1)=real(fwdata2(i,itr))
	dato(2*i)=aimag(fwdata2(i,itr))
enddo

call four1(dato,nb,-1)
do i=1,nb
	wdata2(i,itr)=dato(2*i-1)/nb
enddo

enddo

! !**********************************************************************
! !!!!!!!!!!!!!!!!  Compute the MISFIT AND LIKE !!!!!!!!!!!!!!!!!!!!!!!
! !********************************************************************

!@ LSr est en output est ne represente que le dernier
!-----------------------------
like1=0
LSrt1=0
do itr=1,ntr

      LSr=0
      do i=1,nb
	    br=(wdata1(i,itr)-wdata2(i,itr))
	    LSr=LSr+br**2
	    !write(*,*)br
      enddo

      !br=br/((Ar*lambda(itr))**2)

      !like=like+LSr/(2*(Ar*lambda(itr))**2)
      like1=like1+LSr/(2*(Ar)**2)
      !write(*,*)like,'--------------------------------------'
      LSrt1=LSrt1+LSr
enddo

do  itr=1,ntr
      do i=1,nb
  data11(i,itr)=wdata1(i,itr)
  data21(i,itr)=wdata2(i,itr)
      enddo
enddo

!##################################################
!##################################################
!##################################################
!                               FOR Ps RECEIVER FUNCTIONS
!##################################################


!!!!!!!!!!!!!!!!! Normalize to Unit Energy !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do itr = 1,ntr
      Ene2(itr)=0
      do k=1,nb
	      do j=1,3
		      Ene2(itr)=Ene2(itr)+ tracesP_e(j,k,itr)**2
	      enddo
      enddo
      Ene2(itr)=sqrt(Ene2(itr))
      if (Ene2(itr)==0) then
	    zerotr=1
      endif
enddo

do itr = 1,ntr
do k=1,nb
      NS_e(k,itr) = tracesP_e(1,k,itr)/Ene2(itr)
      EW_e(k,itr) = tracesP_e(3,k,itr)/Ene2(itr)
      NS_o(k,itr) = traces_o1(1,k,itr)
      EW_o(k,itr) = traces_o1(3,k,itr)
      enddo
enddo
! write(*,*)NS_e(25:27,1)
! write(*,*)NS_o(25:27,1)

!  print *,  EW_o(662:670,1)
!  print *, 'a'
  !print *,  EW_e(662:670,1)

sumNS=0
do itr = 1,ntr
do k=1,nb ! We sum over nb !!
	sumNS(itr)=sumNS(itr)+(NS_e(k,itr))**2
enddo
enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Only take the first ndatar in dat(2*i-1)/nb and pass in the frequency domain

!vshort=v
do itr = 1, ntr

    do k=1,nb
	dato(2*k-1) = NS_e(k,itr)
	dato(2*k)=0
	!if (i>460) dato(2*i-1) = 0 !!! TAPERING !!!!!
    enddo

    call four1(dato,nb,1)
    do k=1,nb
	fNS_e(k,itr)=cmplx(dato(2*k-1),dato(2*k))
    enddo

enddo

! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !!! Horizontal ***********************************************

sumEW=0
do itr = 1,ntr
    do k=1,nb ! We sum over nb !!
	    sumEW(itr)=sumEW(itr)+(EW_e(k,itr))**2
    enddo
enddo
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ! Only take the first ndatar in dat(2*i-1)/nb and pass in the frequency domain

do itr = 1, ntr

    do k=1,nb
	dato(2*k-1) = EW_e(k,itr)
	dato(2*k)=0
	!if (i>460) dato(2*i-1) = 0 !!! TAPERING !!!!!
    enddo

    call four1(dato,nb,1)
    do k=1,nb
	fEW_e(k,itr)=cmplx(dato(2*k-1),dato(2*k))
    enddo

enddo

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!sumh=sqrt(sumh)
lambda=0
do itr=1,ntr
    lambda(itr)=sqrt(sumNS(itr)+sumEW(itr))
enddo
 

 
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !*************************************************************
! !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! 
! !!!!!!!!!!!!!!! PASS THE OBSERVATIONS IN THE FREQUENCY DOMAIN !!!!!!!!!
! 
! 
! !!!!!!!! NS !!!!!!!!!!!!!!!!!!!
do itr = 1, ntr
      do k=1,nb
	      dato(2*k-1)=NS_o(k,itr)
	      dato(2*k)=0
      enddo

      call four1(dato,nb,1)
      do k=1,nb
	  fNS_o(k,itr)=cmplx(dato(2*k-1),dato(2*k))
      enddo
enddo

!!!!!!!! EW !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr = 1, ntr
      do k=1,nb
	      dato(2*k-1)=EW_o(k,itr)
	      dato(2*k)=0
      enddo

      call four1(dato,nb,1)
      do k=1,nb
	  fEW_o(k,itr)=cmplx(dato(2*k-1),dato(2*k))
      enddo
enddo

! !******************************************************************
! !!!  COMPUTE THE CROSS PRODUCT !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! !******************************************************************
! 
! !!!!!!!!!!! NS_o * EW_e  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr=1,ntr
do k=1,nb
    fwdata1(k,itr)=fNS_o(k,itr)*fEW_e(k,itr)
enddo
enddo

!Pass in time
do itr=1,ntr

do i=1,nb
	dato(2*i-1)=real(fwdata1(i,itr))
	dato(2*i)=aimag(fwdata1(i,itr))
enddo

call four1(dato,nb,-1)
do i=1,nb
	wdata1(i,itr)=dato(2*i-1)/nb
enddo

enddo

! !!!!!!!!!!!  EW_o * NS_e !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do itr=1,ntr
do i=1,nb
    fwdata2(i,itr)=fEW_o(i,itr)*fNS_e(i,itr)
enddo
enddo

!Pass in time
do itr=1,ntr

do i=1,nb
	dato(2*i-1)=real(fwdata2(i,itr))
	dato(2*i)=aimag(fwdata2(i,itr))
enddo

call four1(dato,nb,-1)
do i=1,nb
	wdata2(i,itr)=dato(2*i-1)/nb
enddo

enddo

! !**********************************************************************
! !!!!!!!!!!!!!!!!  Compute the MISFIT AND LIKE !!!!!!!!!!!!!!!!!!!!!!!
! !********************************************************************

!@ LSr est en output est ne represente que le dernier
!-----------------------------
like3=0
LSrt3=0
do itr=1,ntr

      LSr=0
      do i=1,nb
	    br=(wdata1(i,itr)-wdata2(i,itr))
	    LSr=LSr+br**2
	    !write(*,*)br
      enddo

      !br=br/((Ar*lambda(itr))**2)

      !like=like+LSr/(2*(Ar*lambda(itr))**2)
      like3=like3+LSr/(2*(Ar)**2)
      !write(*,*)like,'--------------------------------------'
      LSrt3=LSrt3+LSr
enddo

do  itr=1,ntr
      do i=1,nb
  data13(i,itr)=wdata1(i,itr)
  data23(i,itr)=wdata2(i,itr)
      enddo
enddo


end

c Calculates a spread of seismograms for a given model.
c Requires readmodel.f,...
c Usage: seis-spread modelfile geometryfile phasefile arrivalfile
C response 
c modified to use Park's reflectivity subroutine instead of ray 
c theoretical. Still uses seis-spread interface (model, geometry, 
c and output file format)
c 

c####&

	subroutine forward_anirec(nlay,ntr,nsamp,dt,baz,slow,thick,
     & beta,alpha,pct_a,pct_b,rho,trend,plunge,strike,dip,
     & tracesP_e,tracesSV_e)

c Anyone who fails to start a Fortran program with this line
c should be severely beaten:
        implicit none
    
c Include constants:
        include 'params.h'
	
	
	real , EXTERNAL    ::    gasdev,ran3
	
c Scratch variables:
        integer i,j,iargc,ii,itr,nsamp,pulse,out_rot
        character modname*(namelen),geomname*(namelen),phname*(namelen)
        character arrname*(namelen),tracename*(namelen)
        character iphname*(namelen)
        real amp_in
c Model parameters:
        integer nlay
        real thick(maxlay),rho(maxlay),alpha(maxlay),beta(maxlay)
        real pct_a(maxlay),pct_b(maxlay),trend(maxlay),plunge(maxlay)
        real strike(maxlay),dip(maxlay),t1,t2,dt,shift
        logical isoflag(maxlay)
c Geometry parameters
        real baz(maxtr),slow(maxtr),sta_dx(maxtr),sta_dy(maxtr)
        integer ntr,ra
c Phase parameters
c        integer phaselist(maxseg,2,maxph),nseg(maxph),numph,
        integer iphase
c Arrivals
c        real travel_time(maxph,maxtr),amplitude(3,maxph,maxtr)
c Traces
        real Tr_cart(3,maxsamp,maxtr),tracesP_e(3,maxsamp,maxtr)
	real Tr_cart3(3,maxsamp,maxtr),tracesSV_e(3,maxsamp,maxtr)
        real Tr_ph(3,maxsamp,maxtr), Tr_cart2(3,maxsamp,maxtr) 
	real width,align


c Determine initial phase

	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PARAMETERS
	!@@ a faire avec parameters
        width=widthh!2*dt !Gaussian pulse width (seconds)
        !! a =sqrt(2)/w  or w =sqrt(2)/a
        out_rot=out_rott !Rotation to output: 0 is EW/Z/NS, 1 is T/Z/R, 2 is SV/SH/P
        !c Determine initial phase
	iphase= 1   !P:1 SV:2 SH:3
	pulse  = 1
	align = 0
     
!         do i=1,nlay
!  	write(*,*)thick(i),rho(i),alpha(i),beta(i),
!      &	pct_a(i),pct_b(i)
!  	enddo
     
!  	write(*,*)'---------------'
!  	write(*,*) ntr, nsamp,dt
!          write(*,*)trend(1:ntr)
!          write(*,*)plunge(1:ntr)
!          write(*,*)nlay,2.*width, shift
!          write(*,*)baz(1:ntr)
!          write(*,*)'---------------'
! !  
! ! 	  write(*,*)out_rot,pulse,slow(1:ntr)
	!do i=1,10
	!write(*,*)'fa1'
	pct_a=0
	pct_b=0
	trend=0
	plunge=0
        call anirec(ntr, nsamp,dt,Tr_cart,Tr_cart2,Tr_cart3,
     c	                thick,rho,alpha, beta, pct_a, pct_b,
     c                  trend,plunge,nlay,2.*width, shifttp,shiftts,
     c			baz,slow, out_rot,pulse)
	
        !enddo
! 	write(*,*)'TRACE'
! 	write(*,*)Tr_cart(3,30:325,1)
! 	write(*,*)'TRACE'
	!if (iphase .eq. 1) then
          if (out_rot .ne. 2) then
!             call writetraces(iounit2,Tr_cart,ntr,nsamp,dt,align,
!      &           	    shift)
	    tracesP_e = Tr_cart
          else
c              Rotate to wavevector coordinates
            call fs_traces_anirec(Tr_cart,slow,alpha(1),beta(1),
     &                       ntr,nsamp,Tr_ph)
            tracesP_e= Tr_ph
            !call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
          end if
        !end if
	
	!if (iphase .eq. 2) then
          if (out_rot .ne. 2) then
!             call writetraces(iounit2,Tr_cart2,ntr,nsamp,dt,align,
!      c           	    shift)
	      tracesSV_e = Tr_cart2
c              Rotate to wavevector coordinates
          else
            call fs_traces_anirec(Tr_cart2,slow,alpha(1),beta(1),
     &                       ntr,nsamp,Tr_ph)
	    tracesSV_e = Tr_ph
!   call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
          end if
        !end if
	
	if (iphase .eq. 3) then
	write(*,*)'iphase==3'
	stop
!           if (out_rot .ne. 2) then
!   !call writetraces(iounit2,Tr_cart3,ntr,nsamp,dt,align,
!   !c             	     shift)
! 	    traces_e = Tr_cart3
!           else
! 
! c              Rotate to wavevector coordinates
!             call fs_traces_anirec(Tr_cart3,slow,alpha(1),beta(1),
!      &                       ntr,nsamp,Tr_ph)
!              traces_e = Tr_ph
!      
!       !      call writetraces(iounit2,Tr_ph,ntr,nsamp,dt,align,shift)
!           end if
        end if
  
                
      end
      
      
      subroutine readparams(filename,mults,nsamp,dt,width,align,shift,
     &                      out_rot,pulse)
      
        implicit none
        include 'params.h'
        
        character filename*(namelen),buffer*(buffsize),firstnonblank
        integer mults,nsamp,align,ios,eof,out_rot,pulse
        real dt,width,shift
        
c          Default values
        mults=0
        nsamp=500
        dt=0.1
        width=2.
        align=0
        shift=0.
        out_rot=0
	pulse=1

        
        open(unit=iounit1,status='old',file=filename,iostat=ios)
        
        if (ios .eq. 0) then
          write(*,*) 'Reading parameters from ',filename
          call getline(iounit1,buffer,eof)
          read (buffer,*) mults
          call getline(iounit1,buffer,eof)
          read (buffer,*) nsamp
          call getline(iounit1,buffer,eof)
          read (buffer,*) dt
          call getline(iounit1,buffer,eof)
          read (buffer,*) width
          call getline(iounit1,buffer,eof)
          read (buffer,*) align
          call getline(iounit1,buffer,eof)
          read (buffer,*) shift
          call getline(iounit1,buffer,eof)
          read (buffer,*) out_rot

creading the source pulse type	
	  read(iounit1,'(A)',iostat=eof) buffer
          if  (((firstnonblank(buffer) .eq. '#') .or.
     c        (firstnonblank(buffer) .eq. ' ')) .and. (eof .eq. 0)) then
            read(iounit1,'(A)',iostat=eof) buffer
            read (buffer,*) pulse
	  else
                write(*,*) 'source type not found in ',filename, 
     c                        '. use default ', pulse
          endif 


        else
          write (*,*) 'Parameter file ',filename,' does not exist.'
          write (*,*) 'Writing default values.'
          close(unit=iounit1)
          open(unit=iounit1,status='unknown',file=filename)
          write(iounit1,*) '# Multiples: 0 for none, 1 for Moho, 2 ',
     &                      'for all first-order, 3 to read file'
          write(iounit1,*) mults
          write(iounit1,*) '# Number of samples per trace'
          write(iounit1,*) nsamp
          write(iounit1,*) '# Sample rate (seconds)'
          write(iounit1,*) dt
          write(iounit1,*) '# Source Pulse width (seconds)'
          write(iounit1,*) width
          write(iounit1,*) '# Alignment: 0 is none, 1 aligns on P'
          write(iounit1,*) align
          write(iounit1,*) '# Shift of traces -- t=0 at this time (sec)'
          write(iounit1,*) shift
          write(iounit1,*) '# Rotation to output: 0 is NS/EW/Z, ',
     &                     '1 is R/T/Z, 2 is P/SV/SH'
          write(iounit1,*) out_rot
  	  write(iounit1,*) '# Source pulse wavelet:,',
     c            ' 1-onesidedpulse 2-oscillation 3-bbpulse 4-hfbb'
	  write(iounit1,*) pulse
  
        end if
        close(unit=iounit1)
                  
      end




c Rotate traces from R-T-Z coordinates to P-SV-SH system.
c Tr_cart and Tr_ph can be the same variable.
      subroutine fs_traces_anirec(Tr_cart,slow,alpha,beta,ntr,
     &                     nsamp,Tr_ph)
      
        implicit none
        include 'params.h'
        
        real Tr_cart(3,maxsamp,maxtr),slow(maxtr),alpha,beta
        real Tr_ph(3,maxsamp,maxtr)
        integer ntr,nsamp
        
        integer itr,isamp,i,j
        real M(3,3),qalpha, qbeta
	
	
	do i = 1,3
	  do j = 1,3
	     M(i,j) = 0.
	  enddo
	enddo
	
	M(3,2) = 0.5
	
        do i = 1, ntr
	   qalpha = sqrt(1. / (alpha*alpha) - slow(i)*slow(i))
	   qbeta =  sqrt(1. / (beta*beta) - slow(i)*slow(i))
	   M(1,1) = slow(i) * beta*beta/alpha
	   M(1,3) = (0.5 -beta*beta*slow(i)*slow(i))/alpha/qalpha
	   M(2,1) = (0.5 -beta*beta*slow(i)*slow(i))/beta/qbeta
	   M(2,3) = -slow(i) * beta
	   do j = 1, nsamp
	       Tr_ph(1,j,i) = Tr_cart(1,j,i) * M(1,1) + 
     c            Tr_cart(3,j,i) * M(1,3)
	       Tr_ph(2,j,i) = Tr_cart(1,j,i) * M(2,1) + 
     c            Tr_cart(3,j,i) * M(2,3)
	       Tr_ph(3,j,i) = 0.5 * Tr_cart(2,j,i)
	   enddo
	enddo
      
      end       
c note the 3 components go as P, SV, and SH




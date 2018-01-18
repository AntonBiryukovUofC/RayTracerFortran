! #########################################
!#########################################

!            PARAMETERS OFR RJ_MCMC

!#########################################
! Declaration
      integer namelen, maxlay,maxtr,out_rott
      integer buffsize,milay,malay,invertype
      integer iounit1,iounit2,maxsamp
      real pi,widthh, shifttp,shiftts,thickmin,vpvs
      integer ndatadmax
      
!#########################################

! IMPORTANT PARAMETERS

! Type of inversion
! 0 : sampling the prior. 1: Ps+SWD. 2: Sp+SWD.
! 3: Ps+Sp+SWD. 4, Only SWD
      parameter (invertype=4)

! min & max number of layers allowed in the \deltaV model
      parameter (milay = 3, malay=45)
  
! Width of the filter when forward modelling traces.
! Should be small (i.e. equal to the sampling rate in observed traces)
! So we don't violate the assumption of uncorrelated noise
      parameter (widthh=0.1)

! Shift in P & S traces (i.e. time before main phase)      
      parameter (shiftts=30,shifttp=30)
      
! Rotation to output: 0 is EW/Z/NS, 1 is T/Z/R, 2 is SV/SH/P   
      parameter (out_rott=1)     
      
! Sacling Parameters      
       parameter  (vpvs = 1.7)
      
!#########################################

! LESS IMPORTANT PARAMETERS

! maxtr is the maximum number of traces allowed
! maxlay is the maximum of layers in \deltaV + Reference model
      parameter (maxlay=100, maxtr=10)

! namelen is the length of filenames (in characters)
      parameter (namelen=40)

! buffsize is the max. line length assumed for reading files.
      parameter (buffsize=120)   

! Units for reading and writing
      parameter (iounit1=1,iounit2=2)   
!PI
      parameter (pi=3.141592653589793)
!Max number of samp per trace
      parameter (maxsamp=5000)

! Thickness minimum for one layer   
! (not activate d at the moment)
      parameter (thickmin=2)  
      
! For Dispersion
       parameter (ndatadmax = 150)
       


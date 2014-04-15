!## D.Pietrobon

program alm_mcmc

  use healpix_modules
  use alm_tools
  use rngmod, only: rand_init, planck_rng 

  implicit none

  type(planck_rng) :: rng_handle 

  integer(i4b), parameter :: nside = 512, &
	&                             pol = 0, &
	&                             nlheader = 60, &
	&                             ncl = 1, &
	&                             nmaps = 1, &
	&                             nfreq = 5!, &
!        &                             nmc = 500

  real(dp), parameter :: logzero = 1.d30

!##  character(len=filenamelen), parameter :: cl_file = '/global/homes/d/dpietrob/myscratch/dx9/maps/ns0128/camb_ctp3_init_cls.fits'

  integer(i4b)            :: nlmax, ifreq, nplm, npix, naccept, nmc
  integer(i4b)            :: i, j, l, m, narg, isample, iseed=-1
!  integer(i8b)            :: n ,ngood,istat 

  real(sp)                :: t1 ,t2 , t3, t4
  real(dp) :: hsqrt2, rms_tt, zeta1_r, zeta1_i, thres, previous_loglike, current_loglike, nullval

  real(dp), allocatable, dimension(:,:) :: cl, tmpbl, plm
  real(dp), dimension(nfreq) :: fwhm
  real(dp), allocatable, dimension(:,:) :: map, mask, rmsmap, chi2, tmpmap
  real(dp), allocatable, dimension(:,:) :: filter, gbeam, pwf
  real(dp), allocatable, dimension(:) :: loglike, tmap, tplm

  complex(dpc), allocatable, dimension(:,:,:) :: current_alms, previous_alms, tmpalms

  character(len=3)        :: sfreq(1:nfreq)
  character(len=4)        :: istring
  character(len=80)       :: header(nlheader)
  character(len=filenamelen)       :: parameter
  character(len=filenamelen), parameter :: root = '/global/scratch/sd/dpietrob/Tools/src/alm_mcmc/'
  character(len=filenamelen)      :: beam_file, tmpfile, maskfile, cl_file


  logical                 :: anynull

  narg = iargc()
  if (narg /= 2) stop 'Calling sequence: ./alm_mcmc nmc sigmal_file.fits'
  call getarg(1,parameter)
  read(parameter,*) nmc

  call getarg(2,parameter)
  read(parameter,*) cl_file

!  call getarg(3,parameter)
!  read(parameter,*) fwhm

   print*, ''
   print*, ' --------------------------------------------------------------'
   print*, '|                                                              |'
   print*, '|                          alms MCMC                           |'
   print*, '|                                                              |'
   print*, ' --------------------------------------------------------------'

  call cpu_time(t3)

  nlmax = 1000
  sfreq = (/'K','K','Q','V','Q'/)
  sfreq(2) = sfreq(2)//'a'
!## print*, ' --> ', sfreq, '<--'

  fwhm = (/0.88,0.66,0.51,0.35,0.22/)
  fwhm = fwhm * 60.d0

  print*, "Computing beams from FWHM..."
!## print*, fwhm
  allocate( gbeam(0:nlmax,nfreq) )
  allocate( pwf(0:nlmax,nfreq) )
  allocate( tmpbl(0:nlmax,1) )
  do i=1,nfreq
     call gaussbeam(fwhm(i), nlmax, tmpbl )
     gbeam(:,i) = tmpbl(:,1)
!## print*, sum(gbeam(:,i))
     call pixel_window( tmpbl, nside )
     pwf(:,i) = tmpbl(:,1)
!## print*, sum(pwf(:,i))
  enddo
  deallocate(tmpbl)
  print*, "Done!"

  print*, "Generating plm..."
  npix = nside2npix(nside) 
  nplm = nside*(nlmax+1) * (2*nlmax-nlmax+2) 
  allocate( plm(0:nplm-1,1) )
  allocate( tplm(0:nplm-1) )
  call plm_gen(nside, nlmax, nlmax, plm)
  tplm = plm(:,1)
  print*, "Done!"

  print*, "Reading Cls in...", trim(cl_file)
  allocate( cl(0:nlmax,ncl) )
  call fits2cl(cl_file, cl, nlmax, ncl, header)
!## print*, sum(cl(:,1))
!## print*, cl(0:10,1)
  print*, "Done!"

  print*, "Read maps in..."
  allocate(map(0:npix-1,nfreq))
  allocate(rmsmap(0:npix-1,nfreq))
  allocate( tmpmap(0:npix-1,1) )
  do i=1,nfreq
     tmpfile = trim(root)//'sims/map_'//trim(sfreq(i))//'_ns0512.fits'
     print*, trim(tmpfile)
     call read_bintab(tmpfile, tmpmap, Npix, Nmaps, nullval, anynull)
     map(:,i) = tmpmap(:,1)
!## print*, sum(map(:,i))/npix

     tmpfile = trim(root)//'sims/rms_'//trim(sfreq(i))//'.fits'
     print*, trim(tmpfile)
     call read_bintab(tmpfile, tmpmap, Npix, Nmaps, nullval, anynull)
     rmsmap(:,i) = tmpmap(:,1)
     call convert_nest2ring(Nside, rmsmap(:,i) )
!## print*, sum(rmsmap(:,i))/npix
  enddo

  maskfile=trim(root)//'sims/wmap_ext_temperature_analysis_mask_r9_7yr_v4.fits'
  print*, "Reading mask in...", trim(maskfile)
  allocate(mask(0:npix-1,1))
  call read_bintab(maskfile, mask, Npix, Nmaps, nullval, anynull)
  print*, "Done!"

  print*, "Generating first alms..."
  allocate( current_alms(1,0:nlmax,0:nlmax) )
  current_alms(1, 0:nlmax, 0:nlmax) = CMPLX(0.0, 0.0, kind=DP)

  allocate( previous_alms(1,0:nlmax,0:nlmax) )
  previous_alms(1, 0:nlmax, 0:nlmax) = CMPLX(0.0, 0.0, kind=DP)

  allocate( tmpalms(1,0:nlmax,0:nlmax) )
  tmpalms(1, 0:nlmax, 0:nlmax) = CMPLX(0.0, 0.0, kind=DP)

!## print*, iseed
  call rand_init(rng_handle, iseed) ! start new sequence
!##  call create_alm(Nside, nlmax, nlmax, pol, cl_file, &
!##	& rng_handle, 0., previous_alms, header) 

    hsqrt2 = SQRT2 / 2.0_dp

    do l = 0, nlmax
       rms_tt = 0.0_dp
       if (cl(l,1) .ne. 0) then
          rms_tt   = sqrt( cl(l,1) )
       endif
!## print*, rms_tt
       !        ------ m = 0 ------
       zeta1_r = rand_gauss(rng_handle)
       zeta1_i = 0.0_dp
       current_alms(1, l, 0)   = previous_alms(1, l, 0) + CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt ! T

       !        ------ m > 0 ------
       do m = 1,l
          zeta1_r = rand_gauss(rng_handle) * hsqrt2
          zeta1_i = rand_gauss(rng_handle) * hsqrt2
!## if (l < 10) print*, zeta1_r, zeta1_i
          current_alms(1, l, m) = previous_alms(1, l, m) + CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt
       enddo
    enddo

!##    call alm2cl( nlmax, nlmax, current_alms, cl )

!## print*, sum(cl(:,1))

  allocate( chi2(0:npix-1,nfreq) )
  chi2 = 0.

  allocate( loglike(nmc) )
  loglike = 0.

  allocate( tmap(0:npix-1) )
  tmap = 0.d0

!## call alm2map( nside, nlmax, nlmax, current_alms, tmap, tplm )
!## print*, sum(tmap)/npix
!## if (sum(tmap) == 0.) stop "Null map!"

  current_loglike = logzero
  previous_loglike = logzero

  print*, "Starting MCMC: ", nmc, " samples"
  naccept = 0

  mcmc: do isample=1,nmc
     call rand_init(rng_handle, iseed-isample) ! start new sequence
  call cpu_time(t1)

     print*, "Sample ",isample,"...", nmc

     print*, "Drawing alms..."
!## --- Taken from alm_map_dd_inc.F90
    !     --- generates randomly the alm according to their power spectrum ---
    !     alm_T = zeta1 * rms_tt
    !     alm_G = zeta1 * rms_g1 + zeta2 * rms_g2
    !     alm_C = zeta3 * rms_cc

    hsqrt2 = SQRT2 / 2.0_dp

    do l = 0, nlmax
       rms_tt = 0.0_dp
!##       rms_g1 = 0.0_dp
       if (cl(l,1) .ne. 0) then
          rms_tt   = sqrt( cl(l,1) )
!##          rms_g1   = cls_tg(l) / rms_tt
       endif

       !        ------ m = 0 ------
       zeta1_r = rand_gauss(rng_handle)
       zeta1_i = 0.0_dp
!       current_alms(1, l, 0)   = previous_alms(1, l, 0) + CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt ! T
       current_alms(1, l, 0)   = CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt ! T
!##       if (polarisation) then
!##          alm_TGC(2, l, 0) = CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_g1 ! G
!##       endif

       !        ------ m > 0 ------
       do m = 1,l
          zeta1_r = rand_gauss(rng_handle) * hsqrt2
          zeta1_i = rand_gauss(rng_handle) * hsqrt2
!          current_alms(1, l, m) = previous_alms(1, l, m) + CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt
          current_alms(1, l, m) = CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_tt
!##          if (polarisation) then
!##             alm_TGC(2, l, m) = CMPLX(zeta1_r, zeta1_i, kind=DP) * rms_g1
!##          endif
       enddo
    enddo
!##
    print*, "Done!"

    print*, "Generating map for each channel..."
    do j=1,nfreq
       tmpalms(1, 0:nlmax, 0:nlmax) = CMPLX(0.0, 0.0, kind=DP)
       do l=0,nlmax
!##          do m=0,l  
             tmpalms(1,l, 0:l) = current_alms(1,l,0:l) * gbeam(l,j) * pwf(l,j) 
!##          enddo
       enddo
!## print*,"Smoothed. Making map..."
       call alm2map( nside, nlmax, nlmax, tmpalms, tmap, tplm )
!## print*, sum(tmap)/npix
!## if (sum(tmap) == 0.) stop "Null map!"
!## print*, "Computing chi2..."
       chi2(0:npix-1,j) = ( (map(0:npix-1,j)-tmap(0:npix-1) ) / rmsmap(0:npix-1,j) )**2
   
       loglike(isample) = loglike(isample) + sum( chi2(0:npix-1,j) * mask(0:npix-1,1) )

    enddo
	
    current_loglike = loglike(isample)
print*, current_loglike, previous_loglike

    if (current_loglike <= previous_loglike) then 
        print*, "Lower chi2: accepted"
        previous_loglike = current_loglike
	previous_alms = current_alms
        print*, "Computing Cls..."
!##        call alm2cl( nlmax, nlmax, previous_alms, cl )
        print*, "Done!"
        do l=0,nlmax
           tmpalms(1,l, 0:l) = previous_alms(1,l,0:l) * pwf(l,j) 
        enddo
        call alm2map( nside, nlmax, nlmax, tmpalms, tmap, tplm)
        print*, "Writing sample map..."
	tmpmap(:,1) = tmap(:)
        call write_minimal_header(header,'MAP',nside=nside,ordering='RING',coordsys='GALACTIC', units='uK')
        write(istring,'(i4.4)') isample
	tmpfile = '!'//trim(root)//'samples/map_no'//istring//'.fits'
	print*, 'Saving into: '//trim(tmpfile)
        call output_map(tmpmap, header,tmpfile)
	naccept = naccept + 1
	print*, 'Saved!'
        
    else
       thres = rand_uni(rng_handle)
       if (exp(-current_loglike/2.)/exp(-previous_loglike/2.) < thres) then
          print*, "Accepted anyway"
        previous_loglike = current_loglike
	previous_alms = current_alms
        print*, "Computing Cls..."
!##        call alm2cl( nlmax, nlmax, previous_alms, cl )
        print*, "Done!"
        do l=0,nlmax
           tmpalms(1,l, 0:l) = previous_alms(1,l,0:l) * pwf(l,j) 
        enddo
        call alm2map( nside, nlmax, nlmax, tmpalms, tmap, tplm)
        print*, "Writing sample map..."
        tmpmap(:,1) = tmap(:)
        call write_minimal_header(header,'MAP',nside=nside,ordering='RING',coordsys='GALACTIC', units='uK')
        write(istring,'(i4.4)') isample
        call output_map(tmpmap, header,'!'//trim(root)//'samples/map_no'//istring//'.fits')
	print*, 'Saved!'
        naccept = naccept + 1
       endif
    endif

    
  call cpu_time(t2)
  print*, 'Sample time:', t2-t1

  enddo mcmc


  call cpu_time(t4)
  print*, 'Job done. See you :)', t4-t3

!##  deallocate( current_alms, gbeam )

end program alm_mcmc

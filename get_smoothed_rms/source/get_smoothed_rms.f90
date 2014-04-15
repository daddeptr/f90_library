!## Kindly provided by L.P.L. Colombo

program get_smoothed_rms

  use healpix_modules
  implicit none

  integer(i4b), parameter :: nside = 2048, &
                             pol = 0, &
                             nlheader = 60, &
                             nside_lfi = 1024, &	
                             ncl = 1, &
                             nfreq = 7, &
                             nsims = 1000

  integer(i4b)            :: nlmax, ifreq

  integer(i8b) ,parameter :: npix = 12*nside**2 ,npix_lfi = 12*nside_lfi**2

  integer(i4b)            :: i,j,l,idum, rnside, rnpix, narg

  integer(i8b)            :: n ,ngood,istat 

  real(dp)                :: fwhm
  real(sp)                :: t1 ,t2 ,dum, t3, t4
  real(dp) ,allocatable   :: w_hfi(:,:),w_lfi(:,:), rmap(:,:)
  real(dp)                :: map(0:npix-1,1) ,tmpmap(0:npix-1,1) &
       ,tmpmap_lfi(0:npix_lfi-1,3), avg(0:npix-1,1), avg2(0:npix-1,1), var(0:npix-1,1) 

!##  real(dp)                :: filter(0:nlmax,1), gbeam(0:nlmax,1)
  real(dp), allocatable, dimension(:,:) :: filter, gbeam
  real(dp)                :: zbounds(2) =[-1,1] ,wrings(1:2*nside,1) = 1.d0
!##  complex(dpc)            :: alms(1,0:nlmax,0:nlmax) 
  complex(dpc), allocatable, dimension(:,:,:) :: alms

  character(len=3)        :: sfreq(1:nfreq)
  character(len=4)        :: istring
  character(len=5)        :: tag
  character(len=80)       :: header(nlheader)
  character(len=256)      :: healpix
  character(len=20)      :: parameter
  character(len=*) ,parameter :: root = '/global/project/projectdirs/planck/data/ffp6/mc_noise/'
  character(len=*) ,parameter :: beam_root = '/global/scratch/sd/paganol/run_FEBeCoP/output_dx9/transFn/Bl/'
  character(len=filenamelen)      :: lfi_tag = beam_root//'Bl_GB_cent_'
  character(len=filenamelen)      :: hfi_tag = beam_root//'Bl_BS_Mars12_cent_'
  character(len=filenamelen)      :: beam_file

  logical                 :: letsgo

  narg = iargc()
  if (narg /= 3) stop 'Calling sequence: ./get_smoothed_rms ifreq nlmax fwhm'
  call getarg(1,parameter)
  read(parameter,*) ifreq 

  call getarg(2,parameter)
  read(parameter,*) nlmax

  call getarg(3,parameter)
  read(parameter,*) fwhm

  sfreq = (/'030','044','070','100','143','217','353'/)

  print*, ''
  print*, ' fwhm  =', fwhm
  print*, ' ifreq =', ifreq
  print*, ' nlmax =', nlmax
  print*, ' sfreq =', sfreq(ifreq)

  allocate( filter(0:nlmax,1) )
  allocate( gbeam(0:nlmax,1) )
  allocate( alms(1,0:nlmax,0:nlmax) )

!##  allocate(w_hfi(0:nlmax,1))
!##  allocate(w_lfi(0:nlmax,1))

!##  call pixel_window(w_lfi,nside_lfi)
!##  call pixel_window(w_hfi,nside)
  
!##  filter = w_hfi/w_lfi
  if (ifreq .lt. 4) then
     beam_file = trim(lfi_tag)//sfreq(ifreq)//'_dx9_CMB.fits' 
     rnside = nside_lfi
     rnpix = 12 * rnside**2
     allocate(rmap(0:rnpix-1,3))
  endif

  if (ifreq .ge. 4) then
     beam_file = trim(hfi_tag)//sfreq(ifreq)//'_dx9_CMB.fits' 
     rnside = nside
     rnpix = 12 * rnside**2
     allocate(rmap(0:rnpix-1,3))
  endif

  print*, 'Reading '//trim(beam_file)

  call gaussbeam(fwhm, nlmax, gbeam) 
  call fits2cl(beam_file, filter, nlmax, ncl, header)
!  deallocate(w_lfi,w_hfi)

  filter = gbeam / filter

  write(tag,'(f5.2)') fwhm 

  call cpu_time(t3)

  ngood = 0
  do i = 1, nsims
!     print*, 'Processing sim ', i
!read rms
     write(istring,'(i4.4)') i-1

     call cpu_time(t1)

!##     call input_map(root // freq //'/ffp6_noise_'//freq //&
!##          '_nominal_map_mc_' // istring // '.fits',tmpmap_lfi,npix_lfi,1)
     call input_map(root // sfreq(ifreq) //'/ffp6_noise_'//sfreq(ifreq) //&
          '_nominal_map_mc_' // istring // '.fits', rmap, rnpix, 1)

     if (ifreq .lt. 4) then
        call udgrade_nest(rmap(:,1),rnside,map(:,1),nside)
     else
         map = rmap
     endif	

     call convert_nest2ring(nside,map(:,1))

     call map2alm(nside,nlmax,nlmax,map(:,1),alms,zbounds,wrings)
     call alter_alm(nside,nlmax,nlmax,0._dp,alms,window=filter)
     call alm2map(nside,nlmax,nlmax,alms,map(:,1))

!##     call write_minimal_header(header,'MAP',nside=nside,ordering='RING',coordsys='GALACTIC')
!##     call output_map(map,header,'!ffp6_noise_'//freq //'_nominal_map_mc_'&
!##          //istring//'.fits')

     avg = avg + map
     avg2 = avg2 + map**2
     
     call cpu_time(t2)
     
     write(*,*) i, t2-t1

     if (mod(i,100) .eq. 0) then
        print*, 'Saving var'
        var = sqrt(avg2/i-(avg/i)**2) * 1.d6
        call write_minimal_header(header,'MAP',nside=nside,ordering='RING',coordsys='GALACTIC', units='uK')
        call output_map(var,header,'!ffp6_rms_'//sfreq(ifreq) //'_'//tag//'arcmin.fits')
	print*, 'Saved!'
     endif
  end do

        print*, 'Saving final var'
        var = sqrt(avg2/nsims-(avg/nsims)**2) * 1.d6
        call write_minimal_header(header,'MAP',nside=nside,ordering='RING',coordsys='GALACTIC', units='uK')
        call output_map(var,header,'!ffp6_rms_'//sfreq(ifreq) //'_'//tag//'arcmin.fits')
	print*, 'Saved!'

  call cpu_time(t4)
  print*, 'Job done. See you :)', t4-t3

  deallocate( alms, gbeam, filter )

end program get_smoothed_rms

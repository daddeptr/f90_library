Module service

   USE healpix_types
   USE fitstools
   USE utilities
   USE pix_tools
   USE paramfile_io

   implicit NONE

   character(len=19), parameter :: code = ' >>> APODIZE_MASK: '
   character(len=4), parameter  :: version='v1.0'

   character(len=filenamelen)   :: in_maskfile, out_maskfile, healpix_dir
   integer(i4b)                 :: feedback, filter
   real(dp)                     :: sigma
   
   contains

   subroutine parsing_hpx( paramfile, code )

      implicit none

      type(paramfile_handle)         :: handle
      character(len=filenamelen), intent(in) :: paramfile
      character(len=*), intent(in)  :: code

      handle = parse_init(paramfile)

      feedback       = parse_int(handle, 'feedback', default=1, descr=code//'Amount of comments printed by the code')
      healpix_dir = get_healpix_main_dir()
      healpix_dir = parse_string(handle, 'healpix_dir', default=TRIM(ADJUSTL(healpix_dir)), descr=code//'Healpix package path')

      sigma       = parse_real(handle, 'sigma', default=60., descr=code//'Apodization angle (arcmin)', vmin=1., vmax=5400.)
      in_maskfile = parse_string(handle, 'in_maskfile', default='', descr='Maskfile')
      out_maskfile = parse_string(handle, 'out_maskfile', default='', descr='Output maskfile')
      filter = parse_int(handle, 'filter_type', default=1, descr='Filter profile: 1-Gaussian, 2-cosine', vmin=1, vmax=2)

      call parse_summarize(handle, code=code)
      call parse_check_unused(handle, code=code)
      call parse_finish(handle)

      return

   end subroutine parsing_hpx 


end module service

program apodize_mask

   USE service
   USE head_fits
   USE mask_tools
   USE extension, only: nArguments, getArgument

   character(len=filenamelen) :: paramfile
   character(len=4) :: io_order
   character(len=8) :: io_sigma
   character(len=80), dimension(80) :: header

   integer(i4b) :: npix, nside, order, nmap
   real(dp), dimension(:,:), allocatable :: mask
   real(dp), dimension(:), allocatable :: distance

   paramfile = ''
   if (nArguments() == 1) call getArgument(1, paramfile)
   call parsing_hpx(paramfile, code)

   if (feedback >= 1) write(*,*) code//'sigma: ', sigma
   npix = getsize_fits(in_maskfile, ordering=order, nside=nside, nmaps=nmap)
   if (feedback >= 1) then
      write(*,*) code//'maskfile: ', trim(adjustl(in_maskfile))
      write(*,*) code//'ordering: ', order
      write(*,*) code//'Nside: ', nside
      write(*,*) code//'Nmaps: ', nmap
   endif

   allocate( mask(0:npix-1,1:nmap) )
   mask = 0.
   if (feedback >= 1) write(*,*) code//'mask allocated.'
   call input_map( in_maskfile, mask, npix, nmap )
   if (feedback >= 1) write(*,*) code//'mask read.'
   if (feedback >= 1) write(*,*) code//'fsky = ',sum(mask)/npix

   if (order==1) call convert_ring2nest( nside, mask )
   if ((feedback >= 1) .and. (order==1) ) write(*,*) code//'mask reordered.'
   allocate( distance(0:npix-1) )
   distance = 0.
   if ((feedback >= 1) .and. (order==1) ) write(*,*) code//'distance allocated.'
   call dist2holes_nest( nside, int(mask(:,1)), distance )
   if (filter == 1)    mask(:,1) = mask(:,1) * (1.-exp(-0.5*(distance/(sigma/60.d0*deg2rad))**2 ) )
   if (feedback >= 1) write(*,*) code//'new fsky = ',sum(mask)/npix

   write(io_sigma,'(f8.3)') sigma
   if (len(trim(adjustl(out_maskfile))) == 0) out_maskfile=trim(adjustl(in_maskfile))//'.apo'//trim(adjustl(io_sigma))

   io_order = 'NEST'
   if (order==1) then
      call convert_nest2ring( nside, mask )
      io_order = 'RING'
   endif
   if (feedback >= 1) then
      write(*,*) code//'ordering: ', io_order
      write(*,*) code//'sigma: ', io_sigma
      write(*,*) code//'saving apodized mask in "'//trim(adjustl(out_maskfile))//'"'
   endif

   call write_minimal_header(header, 'MAP', nside=nside, ordering=io_order, coordsys='G', creator=code, version=version)
   call add_card(header, 'APO_BEAM', sigma)
   call output_map(mask, header, out_maskfile)

   deallocate(distance, mask)

end program


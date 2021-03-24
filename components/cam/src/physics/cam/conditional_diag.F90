module conditional_diag
!-------------------------------------------------
! Conditional diagnostics.
! History:
!  First version by Hui Wan, PNNL, March 2021
!-------------------------------------------------
  use shr_kind_mod, only: r8 => shr_kind_r8
  use spmd_utils,   only: masterproc
  use cam_logfile,  only: iulog

  implicit none

  private

  ! Derived types

  public cnd_diag_info_t
  public cnd_diag_t, cnd_diag

  ! Subroutines

  public conditional_diag_alloc
 !public conditional_diag_dealloc

  ! module parameters

  integer, parameter :: nmetric_max = 10
  integer, parameter :: metric_name_maxlen = 8

  integer, parameter :: nfld_max = 10
  integer, parameter :: fld_name_maxlen = 8 

  integer, parameter :: nphysproc_max  = 100
  integer, parameter :: physproc_name_maxlen  = 20

  integer, parameter :: GEQ = 10
  integer, parameter :: GT  = 11
  integer, parameter :: LEQ = -10
  integer, parameter :: LT  = -11

  !-------------------------------------------------------------------------------
  type cnd_diag_info_t

    ! Do we want to write out the field value after different physical processes?
    logical :: l_output_state = .false.

    ! Do we want to write out tendencies associated with different physicall processes? 
    logical :: l_output_tend  = .false.

    ! Metrics used for conditional sampling.
    ! The current implementation allows the user to define multiple metrics,
    ! each of which will correspond to its own sample and output.
    ! But to keep it simple (at least as a start), we assume that
    ! the physical processes and physical fields to monitor are the same
    ! for different metrics

    integer                       :: nmetric = 0          ! total # of metrics used in this simulation
    character(len=16),allocatable :: metric_name(:)       ! shape = (nmetric); name of the metric
    integer,allocatable           :: metric_nver(:)       ! shape = (nmetric); # of vertical levels
    real(r8),allocatable          :: metric_threshold(:)  ! shape = (nmetric); threshold value for conditional sampling 
    integer,allocatable           :: metric_cmpr_type(:)  ! shape = (nmetric); see module parameters

    ! Physical processes to be monitored
    integer                        :: nphysproc = 0        ! total # of processes
    character(len=16), allocatable :: physproc_name(:)     ! process labels 

    ! Physical fields to be monitored. Each field can have 1, nlev, or nlev+1 vertical levels

    integer                       ::     nfld_1lev = 0
    character(len=16),allocatable :: fld_name_1lev(:)     ! shape = (nfld_1lev)

    integer                       ::     nfld_nlev   = 0
    character(len=16),allocatable :: fld_name_nlev (:)    ! shape = (nfld_nlev)

    integer                       ::     nfld_nlevp  = 0
    character(len=16),allocatable :: fld_name_nlevp(:)    ! shape = (nfld_nlevp)

  end type cnd_diag_info_t

  !-------------------------------------------------------------------------------
  type cnd_diag_t

    real(r8),                   allocatable :: metric_val(:,:)  ! shape = (pcols, info%metric_nver(imetric))
    type(snapshot_and_tendency),allocatable :: fld_1lev (:)     ! shape = (info%nfld_1lev)
    type(snapshot_and_tendency),allocatable :: fld_nlev (:)     ! shape = (info%nfld_nlev)
    type(snapshot_and_tendency),allocatable :: fld_nlevp(:)     ! shape = (info%nfld_nlevp1)

  end type cnd_diag_t

  !---
  type snapshot_and_tendency_t

    real(r8), allocatable :: val(:,:,:) ! shape = (pcols,nver,info%nphysproc) field values after different processes
    real(r8), allocatable :: tnd(:,:,:) ! shape = (pcols,nver,info%nphysproc) tendencies caused by different processes
    real(r8), allocatable :: cur(:,:)   ! shape = (pcols,nver,info%nphysproc) current field values
                                        ! nver is expected to be 1, nlev, or nlev+1
  end type snapshot_and_tendency_t


!===============================================================================
! Module variables
!===============================================================================
  type(cnd_diag_info_t) :: cnd_diag_info

!===============================================================================
! Procedures
!===============================================================================
subroutine conditiona_diag_readnl(nlfile)

!  use namelist_utils,  only: find_group_name
!  use units,           only: getunit, freeunit
!  use mpishorthand
   use infnan, only : nan, assignment(=)

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conditional_diag_readnl'

   ! Local variables for reading namelist

   character(len=metric_name_maxlen)   :: metric_name(nmetric_max) = ' '
   integer  :: metric_nver(nmetric_max)
   integer  :: metric_cmpr_type(nmetric_max)
   real(r8) :: metric_threshold(nmetric_max)

   character(len=physproc_name_maxlen) :: physproc_name(nphysproc_max) = ' '

   character(len=fld_name_maxlen) :: fld_name_1lev (nfld_max) = ' '
   character(len=fld_name_maxlen) :: fld_name_nlev (nfld_max) = ' '
   character(len=fld_name_maxlen) :: fld_name_nlevp(nfld_max) = ' '

   integer :: nmetric = 0
   integer :: nphysproc = 0
   integer :: nfld = 0

   integer :: ii

   namelist /conditional_diag_nl/  metric_name, metric_nver, metric_cmpr_type, metric_threshold, &
                                   physproc_name, fld_name_1lev, fld_name_nlev, fld_name_nlevp

   !----------------------------------------
   !  Default values
   !----------------------------------------
   metric_name      = ' '
   metric_nver      = 0
   metric_cmpr_type = 0
   metric_threshold = nan

   physproc_name  = ' '
   fld_name_1lev  = ' '
   fld_name_nlev  = ' '
   fld_name_nlevp = ' '

   l_output_state = .false.
   l_output_tend  = .false.

   !----------------------------------------
   !  place holder: read namelist !!!!!!!!!!!!!!!!!!!!
   !  BGN: tmp code for testing !!!!!!!!!!!!!!!!!!!!
   !----------------------------------------
   metric_name     (1:3) = (/"CAPE","SSATW","SSATI"/)
   metric_nver     (1:3) = (/1,     72,     72/)
   metric_cmpr_type(1:3) = (/GT,    GT,     GT/)
   metric_metric_threshold(1:3) = (/0._wp, 0._wp, 0._wp/)

   nphysproc = 2
   physproc_name(1) = "zm_conv_tend"
   physproc_name(2) = "cam_radheat" 

   fld_name_1lev(1) = "CAPE"

   ii = 0
   ii = ii + 1 ; fld_name_nlev(ii) = "Q"
   ii = ii + 1 ; fld_name_nlev(ii) = "SSATW"
   ii = ii + 1 ; fld_name_nlev(ii) = "QSATW"
   ii = ii + 1 ; fld_name_nlev(ii) = "SSATI"
   ii = ii + 1 ; fld_name_nlev(ii) = "QSATI"

   !  END: tmp code for testing !!!!!!!!!!!!!!!!!!!!

   !-------------------------------------------
   !  Pack information into to cnd_diag_info
   !-------------------------------------------
   ! metrics for conditional sampling

   ii = 0
   do while ( (ii+1) <= nmetric_max .and. metric_name(ii+1) /= ' ')
      ii = ii + 1
   end do
   nmetric = ii

   cnd_diag_info%nmetric = nmetric

   allocate( cnd_diag_info%metric_name(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%metric_name')
   cnd_diag_info%metric_name(1:nmetric) = trim(adjustl(metric_name(1:nmetric)))

   allocate( cnd_diag_info%metric_nver(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%metric_nver')
   cnd_diag_info%metric_nver(1:nmetric) = metric_nver(1:nmetric)

   allocate( cnd_diag_info%metric_cmpr_type(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%metric_cmpr_type')
   cnd_diag_info%metric_cmpr_type(1:nmetric) = metric_cmpr_type(1:nmetric)

   allocate( cnd_diag_info%metric_threshold(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%metric_threshold')
   cnd_diag_info%metric_threshold(1:nmetric) = metric_threshold(1:nmetric)

   ! physical processes to monitor 

   ii = 0
   do while ( (ii+1) <= nphysproc_max .and. physproc_name(ii+1) /= ' ')
      ii = ii + 1
   end do
   nphysproc = ii

   cnd_diag_info%nphysproc = nphysproc

   allocate( cnd_diag_info%physproc_name(nphysproc), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%physproc_name')
   cnd_diag_info%physproc_name(1:nphysproc) = trim(adjustl(physproc_name(1:nphysproc)))

   ! fields with 1 vertitical level

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_1lev(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_1lev = nfld

   allocate( cnd_diag_info%fld_name_1lev(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%fld_name_1lev')
   cnd_diag_info%fld_name_1lev(1:nfld) = trim(adjustl(fld_name_1lev(1:nfld)))

   ! fields with nlev vertitical levels 

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_nlev(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_nlev = nfld

   allocate( cnd_diag_info%fld_name_nlev(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%fld_name_nlev')
   cnd_diag_info%fld_name_nlev(1:nfld) = trim(adjustl(fld_name_nlev(1:nfld)))

   ! fields with nlev+1 vertitical levels 

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_nlevp(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_nlevp = nfld

   allocate( cnd_diag_info%fld_name_nlevp(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun('conditional_diag_readnl: allocation error for cnd_diag_info%fld_name_nlevp')
   cnd_diag_info%fld_name_nlevp(1:nfld) = trim(adjustl(fld_name_nlevp(1:nfld)))

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

    if (cnd_diag_info%nmetric == 0) then

      write(iulog,*)' -----------------------------------------------------'
      write(iulog,*)'     *** Conditional diagnostics NOT requested ***'
      write(iulog,*)' -----------------------------------------------------'

    else

      write(iulog,*)' --------------------------------------------------'
      write(iulog,*)'     *** Conditional diagnostics requested ***'
      write(iulog,*)' --------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,a12,i12,i12,e20.10)')'metric','nlev','cmpr type','threshold'
      do ii = 1,cnd_diag_info%nmetric
         write(iulog,'(i4.3,a12,i12,i12,e20.10)') ii, cnd_diag_info%metric_name(ii), &
                                                      cnd_diag_info%metric_nver(ii), &
                                                      cnd_diag_info%metric_cmpr_type(ii), &
                                                      cnd_diag_info%metric_threshold(ii)   
      end do

      write(iulog,*)
      write(iulog,'(4x,a20)')'physical process'
      do ii = 1,cnd_diag_info%nphysproc
         write(iulog,'(i4.3,a20)') ii, cnd_diag_info%physproc_name(ii)
      end do

      write(iulog,*)
      write(iulog,'(4x,a40)')'field w/ 1 vertical level'
      do ii = 1,cnd_diag_info%nfld_1lev
         write(iulog,'(i4.3,a40)') ii, cnd_diag_info%fld_name_1lev(ii)
      end do

      write(iulog,*)
      write(iulog,'(4x,a40)')'field w/ nlev vertical levels'
      do ii = 1,cnd_diag_info%nfld_nlev
         write(iulog,'(i4.3,a40)') ii, cnd_diag_info%fld_name_nlev(ii)
      end do

      write(iulog,*)
      write(iulog,'(4x,a40)')'field w/ nlev+1 vertical levels'
      do ii = 1,cnd_diag_info%nfld_nlevp1
         write(iulog,'(i4.3,a40)') ii, cnd_diag_info%fld_name_nlevp(ii)
      end do

      write(iulog,*)  'l_output_state = ',l_output_state
      write(iulog,*)  'l_output_tend  = ',l_output_tend
      write(iulog,*)' --------------------------------------------------'

end subroutine conditiona_diag_readnl

!===============================================================================
subroutine conditional_diag_alloc( psetcols, pver, info, diag )

  use infnan, only : inf, assignment(=)

  integer, intent(in) :: psetcols, pver
  type(cnd_diag_info_t), intent(in) :: info
  type(cnd_diag_t),   intent(inout) :: diag

  integer :: nmetric, im
  integer :: nfld_1lev, nfld_nlev, nfld_nlevp, ifld
  integer :: nphysproc
  integer :: ierr

  nmetric    = info% nmetric
  nfld_1lev  = info% nfld_1lev
  nfld_nlev  = info% nfld_nlev
  nfld_nlevp = info% nfld_nlevp
  nphysproc = info% nphysproc

  ! different groups of diagnostics correspond to different conditions

  allocate( diag(nmetric), stat=ierr)
  if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag')

  ! the metric field, which might have 1, pver, or pver+1 vertical levels

  do im = 1,nmetric
     allocate( diag(im)% metric( psetcols, info%metric_nver(im) ), stat=ierr)
     if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag%metric')
  end do

  ! diagnostics with only vertical level

  if (nfld_1lev > 0) then

     do im = 1,nmetric

        allocate( diag(im)% fld_1lev( nfld_1lev ), stat=ierr)
        if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag%fld_1lev')

        ! field values and tendencies of each diagnostic variables with 1 vertical levels
        do ifld = 1, nfld_1lev

           allocate( diag(im)%fld_1lev(ifld)% cur(npsetcols,1), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag%fld_1lev%cur')

           allocate( diag(im)%fld_1lev(ifld)% val(npsetcols,1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag%fld_1lev%val')

           allocate( diag(im)%fld_1lev(ifld)% tnd(npsetcols,1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_allocate error: allocation error for diag%fld_1lev%tnd')

           diag(im)%fld_1lev(ifld)% cur(:,:) = inf
           diag(im)%fld_1lev(ifld)% val(:,:) = inf
           diag(im)%fld_1lev(ifld)% tnd(:,:) = inf

        end do !ifld

      end do   !im
   end if

   ! diagnostics with pver vertical levels (typically physical variables located at layer midpoints)

   if (nfld_nlev > 0) then

      do im = 1,nmetric

         allocate( diag(im)% fld_nlev( nfld_nlev ), stat=ierr)
         if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlev')

         ! field values and tendencies of each diagnostic variables with 1 vertical levels
         do ifld = 1, nfld_nlev

            allocate( diag(im)%fld_nlev(ifld)% cur(npsetcols,pver), stat=ierr)
            if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlev%cur')

            allocate( diag(im)%fld_nlev(ifld)% val(npsetcols,pver,nphysproc), stat=ierr)
            if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlev%val')

            allocate( diag(im)%fld_nlev(ifld)% tnd(npsetcols,pver,nphysproc), stat=ierr)
            if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlev%tnd')

            diag(im)%fld_nlev(ifld)% cur(:,:) = inf
            diag(im)%fld_nlev(ifld)% val(:,:) = inf
            diag(im)%fld_nlev(ifld)% tnd(:,:) = inf

         end do !ifld

      end do   !im
   end if

   ! diagnostics with pver vertical levels (typically physical variables located at layer interfaces)
   if (nfld_nlevp > 0) then

      do im = 1,nmetric

        allocate( diag(im)% fld_nlevp( nfld_nlevp ), stat=ierr)
        if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlevp')

        ! field values and tendencies of each diagnostic variables with 1 vertical levels
        do ifld = 1, nfld_nlevp

           allocate( diag(im)%fld_nlevp(ifld)% cur(npsetcols,pver+1), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlevp%cur')

           allocate( diag(im)%fld_nlevp(ifld)% val(npsetcols,pver+1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlevp%val')

           allocate( diag(im)%fld_nlevp(ifld)% tnd(npsetcols,pver+1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun('conditional_diag_alloc error: allocation error for diag%fld_nlevp%tnd')

           diag(im)%fld_nlevp(ifld)% cur(:,:) = inf
           diag(im)%fld_nlevp(ifld)% val(:,:) = inf
           diag(im)%fld_nlevp(ifld)% tnd(:,:) = inf

        end do !ifld
      end do   !im

   end if

end subroutine conditional_diag_alloc

!subroutine conditional_diag_dealloc
!end subroutine conditional_diag_alloc



end module conditional_diag

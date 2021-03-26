module conditional_diag
!-------------------------------------------------
! Conditional diagnostics.
! History:
!  First version by Hui Wan, PNNL, March 2021
!-------------------------------------------------
  use shr_kind_mod,   only: r8 => shr_kind_r8
  use spmd_utils,     only: masterproc
  use cam_logfile,    only: iulog
  use cam_abortutils, only: endrun

  implicit none

  private

  ! Derived types

  public cnd_diag_info_t
  public cnd_diag_t, cnd_diag_info

  ! Subroutines

  public conditional_diag_readnl
  public conditional_diag_alloc
 !public conditional_diag_dealloc
  public conditional_diag_output_init

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

    real(r8),                     allocatable :: metric (:,:)     ! shape = (pcols, info%metric_nver(imetric))
    type(snapshot_and_tendency_t),allocatable :: fld_1lev (:)     ! shape = (info%nfld_1lev)
    type(snapshot_and_tendency_t),allocatable :: fld_nlev (:)     ! shape = (info%nfld_nlev)
    type(snapshot_and_tendency_t),allocatable :: fld_nlevp(:)     ! shape = (info%nfld_nlevp1)

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

contains
!===============================================================================
! Procedures
!===============================================================================
subroutine conditional_diag_readnl(nlfile)

   use infnan,          only: nan, assignment(=)
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conditional_diag_readnl'

   ! Local variables for reading namelist

   character(len=metric_name_maxlen)   :: metric_name(nmetric_max)
   integer  :: metric_nver(nmetric_max)
   integer  :: metric_cmpr_type(nmetric_max)
   real(r8) :: metric_threshold(nmetric_max)

   character(len=physproc_name_maxlen) :: physproc_name(nphysproc_max)

   character(len=fld_name_maxlen) :: fld_name_1lev (nfld_max)
   character(len=fld_name_maxlen) :: fld_name_nlev (nfld_max)
   character(len=fld_name_maxlen) :: fld_name_nlevp(nfld_max)

   logical :: l_output_state, l_output_tend

   ! other misc local variables
   integer :: nmetric
   integer :: nphysproc
   integer :: nfld

   integer :: ii

   !-------
   namelist /conditional_diag_nl/  &
            metric_name, metric_nver, metric_cmpr_type, metric_threshold, &
            physproc_name, fld_name_1lev, fld_name_nlev, fld_name_nlevp,  &
            l_output_state, l_output_tend

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
   ! Read namelist
   !----------------------------------------
   if (masterproc) then
      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conditional_diag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conditional_diag_nl, iostat=ierr)
         if (ierr /= 0) then
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)
   end if ! masterproc

#ifdef SPMD
   ! Broadcast namelist variables

   call mpibcast(metric_name,      nmetric_max*len(metric_name(1)), mpichar, 0, mpicom)
   call mpibcast(metric_nver,      nmetric_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_cmpr_type, nmetric_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_threshold, nmetric_max,                     mpir8,   0, mpicom)

   call mpibcast(physproc_name,  nphysproc_max*len(physproc_name(1)), mpichar, 0, mpicom)

   call mpibcast(fld_name_1lev,  nfld_max*len(fld_name_1lev(1)),  mpichar, 0, mpicom)
   call mpibcast(fld_name_nlev,  nfld_max*len(fld_name_nlev(1)),  mpichar, 0, mpicom)
   call mpibcast(fld_name_nlevp, nfld_max*len(fld_name_nlevp(1)), mpichar, 0, mpicom)

   call mpibcast(l_output_state, 1, mpilog, 0, mpicom)
   call mpibcast(l_output_tend,  1, mpilog, 0, mpicom)
#endif

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
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%metric_name')

   do ii = 1,nmetric
      cnd_diag_info%metric_name(ii) = trim(adjustl(metric_name(ii)))
   end do

   allocate( cnd_diag_info%metric_nver(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%metric_nver')
   cnd_diag_info%metric_nver(1:nmetric) = metric_nver(1:nmetric)

   allocate( cnd_diag_info%metric_cmpr_type(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%metric_cmpr_type')
   cnd_diag_info%metric_cmpr_type(1:nmetric) = metric_cmpr_type(1:nmetric)

   allocate( cnd_diag_info%metric_threshold(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%metric_threshold')
   cnd_diag_info%metric_threshold(1:nmetric) = metric_threshold(1:nmetric)

   ! physical processes to monitor 

   ii = 0
   do while ( (ii+1) <= nphysproc_max .and. physproc_name(ii+1) /= ' ')
      ii = ii + 1
   end do
   nphysproc = ii

   cnd_diag_info%nphysproc = nphysproc

   allocate( cnd_diag_info%physproc_name(nphysproc), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%physproc_name')
   do ii = 1,nphysproc
      cnd_diag_info%physproc_name(ii) = trim(adjustl(physproc_name(ii)))
   end do

   ! fields with 1 vertitical level

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_1lev(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_1lev = nfld

   allocate( cnd_diag_info%fld_name_1lev(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%fld_name_1lev')
   do ii = 1,nfld
      cnd_diag_info%fld_name_1lev(ii) = trim(adjustl(fld_name_1lev(ii)))
   end do

   ! fields with nlev vertitical levels 

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_nlev(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_nlev = nfld

   allocate( cnd_diag_info%fld_name_nlev(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%fld_name_nlev')
   do ii = 1,nfld
      cnd_diag_info%fld_name_nlev(ii) = trim(adjustl(fld_name_nlev(ii)))
   end do

   ! fields with nlev+1 vertitical levels 

   ii = 0
   do while ( (ii+1) <= nfld_max .and. fld_name_nlevp(ii+1) /= ' ')
      ii = ii + 1
   end do
   nfld = ii

   cnd_diag_info%nfld_nlevp = nfld

   allocate( cnd_diag_info%fld_name_nlevp(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation error for cnd_diag_info%fld_name_nlevp')
   do ii = 1,nfld
      cnd_diag_info%fld_name_nlevp(ii) = trim(adjustl(fld_name_nlevp(ii)))
   end do

   ! output to history file(s)

   cnd_diag_info%l_output_state = l_output_state
   cnd_diag_info%l_output_tend  = l_output_tend 

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

    if (cnd_diag_info%nmetric == 0) then

      write(iulog,*)'==========================================================='
      write(iulog,*)'       *** Conditional diagnostics NOT requested ***'
      write(iulog,*)'==========================================================='

    else

      write(iulog,*)'==========================================================='
      write(iulog,*)'       *** Conditional diagnostics requested ***'
      write(iulog,*)'-----------------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,2x,a10,a6,a12,a20)')'metric','nlev','cmpr type','threshold'
      do ii = 1,cnd_diag_info%nmetric
         write(iulog,'(i4.3,2x,a10,i6,i12,e20.10)') ii, cnd_diag_info%metric_name(ii), &
                                                        cnd_diag_info%metric_nver(ii), &
                                                        cnd_diag_info%metric_cmpr_type(ii), &
                                                        cnd_diag_info%metric_threshold(ii)
      end do

      write(iulog,*)
      write(iulog,'(4x,a20)')'physical process'
      do ii = 1,cnd_diag_info%nphysproc
         write(iulog,'(i4.3,a20)') ii, cnd_diag_info%physproc_name(ii)
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,a30)')'field w/ 1 vertical level'
      do ii = 1,cnd_diag_info%nfld_1lev
         write(iulog,'(i4.3,a30)') ii, cnd_diag_info%fld_name_1lev(ii)
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,a30)')'field w/ nlev vertical levels'
      do ii = 1,cnd_diag_info%nfld_nlev
         write(iulog,'(i4.3,a30)') ii, cnd_diag_info%fld_name_nlev(ii)
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,a30)')'field w/ nlev+1 vertical levels'
      do ii = 1,cnd_diag_info%nfld_nlevp
         write(iulog,'(i4.3,a30)') ii, cnd_diag_info%fld_name_nlevp(ii)
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,*)' l_output_state = ',l_output_state
      write(iulog,*)' l_output_tend  = ',l_output_tend
      write(iulog,*)
      write(iulog,*)'==========================================================='
      write(iulog,*)

   end if ! cnd_diag_info%nmetric == 0
  end if  ! masterproc

end subroutine conditional_diag_readnl

!===============================================================================
subroutine conditional_diag_alloc( psetcols, pver, metric_nver, nphysproc, &
                                   nfld_1lev, nfld_nlev, nfld_nlevp, diag )

  use infnan, only : inf, assignment(=)

  integer, intent(in) :: psetcols, pver
  integer, intent(in) :: metric_nver, nphysproc
  integer, intent(in) :: nfld_1lev, nfld_nlev, nfld_nlevp

  type(cnd_diag_t), intent(inout) :: diag

  integer :: ifld
  integer :: ierr

  character(len=*), parameter :: subname = 'conditional_diag_alloc'

  ! the metric fields, which might have 1, pver, or pver+1 vertical levels

  allocate( diag% metric(psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%metric')

  ! diagnostical fields with only vertical level

  if (nfld_1lev > 0) then

        allocate( diag% fld_1lev( nfld_1lev ), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_1lev')

        do ifld = 1, nfld_1lev  ! values and tendencies of each field

           allocate( diag%fld_1lev(ifld)% cur(psetcols,1), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_1lev%cur')

           allocate( diag%fld_1lev(ifld)% val(psetcols,1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_1lev%val')

           allocate( diag%fld_1lev(ifld)% tnd(psetcols,1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_1lev%tnd')

           diag%fld_1lev(ifld)% cur(:,:)   = inf
           diag%fld_1lev(ifld)% val(:,:,:) = inf
           diag%fld_1lev(ifld)% tnd(:,:,:) = inf

        end do !ifld

   end if

   ! diagnostical fields with pver vertical levels (typically physical variables located at layer midpoints)

   if (nfld_nlev > 0) then

         allocate( diag% fld_nlev( nfld_nlev ), stat=ierr)
         if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlev')

         do ifld = 1, nfld_nlev  ! values and tendencies of each field

            allocate( diag%fld_nlev(ifld)% cur(psetcols,pver), stat=ierr)
            if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlev%cur')

            allocate( diag%fld_nlev(ifld)% val(psetcols,pver,nphysproc), stat=ierr)
            if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlev%val')

            allocate( diag%fld_nlev(ifld)% tnd(psetcols,pver,nphysproc), stat=ierr)
            if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlev%tnd')

            diag%fld_nlev(ifld)% cur(:,:)   = inf
            diag%fld_nlev(ifld)% val(:,:,:) = inf
            diag%fld_nlev(ifld)% tnd(:,:,:) = inf

         end do !ifld

   end if

   ! diagnostical fields with pver+1 vertical levels (typically physical variables located at layer interfaces)
   if (nfld_nlevp > 0) then

        allocate( diag% fld_nlevp( nfld_nlevp ), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlevp')

        do ifld = 1, nfld_nlevp  ! values and tendencies of each field

           allocate( diag%fld_nlevp(ifld)% cur(psetcols,pver+1), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlevp%cur')

           allocate( diag%fld_nlevp(ifld)% val(psetcols,pver+1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlevp%val')

           allocate( diag%fld_nlevp(ifld)% tnd(psetcols,pver+1,nphysproc), stat=ierr)
           if ( ierr /= 0 ) call endrun(subname//': allocation error for diag%fld_nlevp%tnd')

           diag%fld_nlevp(ifld)% cur(:,:)   = inf
           diag%fld_nlevp(ifld)% val(:,:,:) = inf
           diag%fld_nlevp(ifld)% tnd(:,:,:) = inf

        end do !ifld

   end if

end subroutine conditional_diag_alloc

!subroutine conditional_diag_dealloc
!end subroutine conditional_diag_alloc

subroutine conditional_diag_output_init(pver)

  use cam_history,         only: addfld, horiz_only, add_default
  use cam_history_support, only: max_fieldname_len

  integer,intent(in) :: pver

  integer          :: im, ifld, iphys, ii
  character(len=2) :: imstr
  character(len=4) :: val_tnd_suff(2), suff
  logical          :: l_cycle(2)

  character(len=max_fieldname_len) :: output_fld_name

  character(len=*),parameter :: subname = 'conditional_diag_output_init'

  if (cnd_diag_info%nmetric==0) return

  do im = 1,cnd_diag_info%nmetric

     ! imstr is the metric index as a string; will be appended to output field names
     write(imstr,'(i2.2)') im

     !----------------------------------------
     ! register the metric itself for output
     !----------------------------------------
     output_fld_name = trim(cnd_diag_info%metric_name(im))//'_cnd'//imstr

     if (cnd_diag_info%metric_nver(im)==1) then

       call addfld(trim(output_fld_name), horiz_only, 'A',' ',' ') 

     elseif(cnd_diag_info%metric_nver(im)==pver) then

       call addfld(trim(output_fld_name), (/'lev'/),  'A',' ',' ') 

     elseif(cnd_diag_info%metric_nver(im)==pver+1) then

       call addfld(trim(output_fld_name), (/'ilev'/), 'A',' ',' ') 

     else 
       call endrun(subname//': invalid number of vertical levers')
     end if 

     call add_default(trim(output_fld_name),1,' ')

     !----------------------------------------
     ! register the diagnostics for output
     !----------------------------------------
     val_tnd_suff = (/"_val","_tnd"/)
     l_cycle      = .not.(/cnd_diag_info%l_output_state, cnd_diag_info%l_output_tend/)

     do ii = 1,2 ! field value (state) or tendency

        if (l_cycle(ii)) cycle

        suff = val_tnd_suff(ii)

        ! diagnostic fields with no vertical distribution

        do ifld = 1,cnd_diag_info%nfld_1lev
        do iphys = 1,cnd_diag_info%nphysproc

           output_fld_name = trim(cnd_diag_info%fld_name_1lev(ifld))//'_cnd'//imstr//'_'// &
                             trim(cnd_diag_info%physproc_name(iphys))//suff

           call addfld(trim(output_fld_name), horiz_only,  'A',' ',' ') 
           call add_default(trim(output_fld_name),1,' ')
        end do
        end do

        ! diagnostic fields located at layer midpoints

        do ifld = 1,cnd_diag_info%nfld_nlev
        do iphys = 1,cnd_diag_info%nphysproc

           output_fld_name = trim(cnd_diag_info%fld_name_nlev(ifld))//'_cnd'//imstr//'_'// &
                             trim(cnd_diag_info%physproc_name(iphys))//suff

           call addfld(trim(output_fld_name), (/'lev'/),  'A',' ',' ') 
           call add_default(trim(output_fld_name),1,' ')
        end do
        end do

        ! diagnostic fields located at layer interfaces

        do ifld = 1,cnd_diag_info%nfld_nlevp
        do iphys = 1,cnd_diag_info%nphysproc

           output_fld_name = trim(cnd_diag_info%fld_name_nlevp(ifld))//'_cnd'//imstr//'_'// &
                             trim(cnd_diag_info%physproc_name(iphys))//suff

           call addfld(trim(output_fld_name), (/'ilev'/),  'A',' ',' ') 
           call add_default(trim(output_fld_name),1,' ')
        end do
        end do

     end do ! ii = 1,2, field value (state) or tendency

  end do ! im = 1,nmetric

end subroutine conditional_diag_output_init


end module conditional_diag

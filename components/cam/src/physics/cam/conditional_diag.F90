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
  public cnd_diag_t

  ! Variable(s) of derived type

  public cnd_diag_info

  ! Subroutines

  public conditional_diag_readnl
  public conditional_diag_alloc
 !public conditional_diag_dealloc
  public conditional_diag_output_init

  ! module parameters

  integer, parameter :: nmetric_max = 10
  integer, parameter :: mname_maxlen = 8

  integer, parameter :: nfld_max = 20
  integer, parameter :: fname_maxlen = 8 

  integer, parameter :: nphysproc_max   = 100
  integer, parameter :: ptndname_maxlen = 20
  integer, parameter :: poutname_maxlen = 3

  integer, parameter :: GT  = 1
  integer, parameter :: GEQ = 2
  integer, parameter :: LT  = -1
  integer, parameter :: LEQ = -2

  !-------------------------------------------------------------------------------
  type cnd_diag_info_t

    ! Do we want to write out the field value after different physical processes?
    logical :: l_output_state = .false.

    ! Do we want to write out increments caused by different physicall processes? 
    logical :: l_output_incrm = .false.

    ! Metrics used for conditional sampling.
    ! The current implementation allows the user to define multiple metrics,
    ! each of which will correspond to its own sample and output.
    ! But to keep it simple (at least as a start), we assume that
    ! the physical processes and physical fields to monitor are the same
    ! for different metrics

    integer                      :: nmetric = 0          ! total # of metrics used in this simulation
    character(len=mname_maxlen),&
                     allocatable :: metric_name(:)       ! shape = (nmetric); name of the metric
    integer,allocatable          :: metric_nver(:)       ! shape = (nmetric); # of vertical levels
    real(r8),allocatable         :: metric_fillvalue(:)  ! shape = (nmetric); fill value used in conditional sampling 
    real(r8),allocatable         :: metric_threshold(:)  ! shape = (nmetric); threshold value for conditional sampling 
    integer,allocatable          :: metric_cmpr_type(:)  ! shape = (nmetric); see module parameters
    character(len=ptndname_maxlen),&
                     allocatable :: sample_after(:)      ! shape = (nmetric); after which atmospheric process
                                                         ! will conditional sampling be applied? The process names
                                                         ! need to match ptend%name)

    ! Physical processes to be monitored
    integer                          :: nphysproc = 0    ! total # of processes
    character(len=ptndname_maxlen),&
                         allocatable :: ptend_name(:)    ! process labels (need to match ptend%name)
    character(len=poutname_maxlen),&
                         allocatable :: proc_outname(:)  ! process labels (user specified 3-character 
                                                         ! labels for keeping output variable names short 

    ! Physical fields to be monitored. Each field can have 1, nlev, or nlev+1 vertical levels

    integer                                 :: nfld = 0
    character(len=fname_maxlen),allocatable :: fld_name(:)     ! shape = (nfld)
    integer,allocatable                     :: fld_nver(:)     ! shape = (nfld); # of vertical levels

  end type cnd_diag_info_t

  !-------------------------------------------------------------------------------
  type cnd_diag_t

    real(r8),                      allocatable :: metric (:,:)     ! shape = (pcols, info%metric_nver(imetric))
    real(r8),                      allocatable :: flag   (:,:)     ! shape = (pcols, info%metric_nver(imetric))
    type(snapshot_and_increment_t),allocatable :: fld    (:)       ! shape = (info%nfld)

  end type cnd_diag_t

  !---
  type snapshot_and_increment_t

    real(r8), allocatable :: val(:,:,:) ! shape = (pcols,info%fld_nver(ifld),info%nphysproc) field values after different processes
    real(r8), allocatable :: inc(:,:,:) ! shape = (pcols,info%fld_nver(ifld),info%nphysproc) increments caused by different processes
    real(r8), allocatable :: old(:,:)   ! shape = (pcols,info%fld_nver(ifld),info%nphysproc) old field values

  end type snapshot_and_increment_t


!===============================================================================
! Module variables
!===============================================================================
  type(cnd_diag_info_t) :: cnd_diag_info

contains
!===============================================================================
! Procedures
!===============================================================================
subroutine conditional_diag_readnl(nlfile)

   use infnan,          only: nan, assignment(=), isnan
   use namelist_utils,  only: find_group_name
   use units,           only: getunit, freeunit
   use mpishorthand

   character(len=*), intent(in) :: nlfile  ! filepath for file containing namelist input

   ! Local variables
   integer :: unitn, ierr
   character(len=*), parameter :: subname = 'conditional_diag_readnl'

   ! Local variables for reading namelist

   character(len=mname_maxlen)   :: metric_name     (nmetric_max)
   integer                       :: metric_nver     (nmetric_max)
   integer                       :: metric_cmpr_type(nmetric_max)
   real(r8)                      :: metric_threshold(nmetric_max)
   real(r8)                      :: metric_fillvalue(nmetric_max)
   character(len=ptndname_maxlen):: sample_after    (nmetric_max)

   character(len=ptndname_maxlen) :: ptend_name(nphysproc_max)
   character(len=poutname_maxlen) :: proc_outname(nphysproc_max)

   character(len=fname_maxlen)    :: fld_name (nfld_max)
   integer                        :: fld_nver (nfld_max)

   logical :: l_output_state, l_output_incrm

   ! other misc local variables
   integer :: nmetric
   integer :: nphysproc
   integer :: nfld

   integer :: ii

   !-------
   namelist /conditional_diag_nl/  &
            metric_name, metric_nver, metric_cmpr_type, metric_threshold, metric_fillvalue, sample_after, &
            ptend_name, proc_outname, &
            fld_name, fld_nver,  &
            l_output_state, l_output_incrm

   !----------------------------------------
   !  Default values
   !----------------------------------------
   metric_name      = ' '
   metric_nver      = 0
   metric_cmpr_type = 0
   metric_threshold = nan
   metric_fillvalue = nan
   sample_after     = ' '

   ptend_name     = ' '
   proc_outname   = ' '

   fld_name       = ' '
   fld_nver       = 0 

   l_output_state = .false.
   l_output_incrm = .false.

   !----------------------------------------
   ! Read namelist and check validity
   !----------------------------------------
   if (masterproc) then

      ! Read namelist

      unitn = getunit()
      open( unitn, file=trim(nlfile), status='old' )
      call find_group_name(unitn, 'conditional_diag_nl', status=ierr)
      if (ierr == 0) then
         read(unitn, conditional_diag_nl, iostat=ierr)
         if (ierr /= 0) then
            write(iulog,*) 'read error ',ierr
            call endrun(subname // ':: ERROR reading namelist')
         end if
      end if
      close(unitn)

      ! Check validity of namelist variables for user-specified metrics

      ii = 0
      do while ( (ii+1) <= nmetric_max .and. metric_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      nmetric = ii

      if (any( metric_nver     (1:nmetric)<=0     )) call endrun(subname//' error: need non-zero metric_nver for each metric_name')
      if (any( metric_cmpr_type(1:nmetric)==0     )) call endrun(subname//' error: need valid metric_cmpr_type for each metric_name')
      if (any( isnan(metric_threshold(1:nmetric)) )) call endrun(subname//' error: need valid metric_threshold for each metric_name')
      if (any( isnan(metric_fillvalue(1:nmetric)) )) call endrun(subname//' error: need valid metric_fillvalue for each metric_name')
      if (any( sample_after    (1:nmetric)==' '   )) call endrun(subname//' error: be sure to specify sample_after for each metric_name')

      ! Check validity of namelist variables for atmospheric processes to monitor

      ii = 0
      do while ( (ii+1) <= nphysproc_max .and. ptend_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      nphysproc = ii

      if (any(proc_outname(1:nphysproc)==' ')) call endrun(subname//'error: be sure to specify proc_outname for each ptend_name.')

      ! Check validity of namelist variables for physical fields to monitor

      ii = 0
      do while ( (ii+1) <= nfld_max .and. fld_name(ii+1) /= ' ')
         ii = ii + 1
      end do
      nfld = ii

      if (any(fld_nver(1:nfld)<=0)) call endrun(subname//'error: need positive fld_nver for each fld_name')

   end if ! masterproc
   !--------------------------------------

#ifdef SPMD
   call mpibcast(nmetric,  1, mpiint, 0, mpicom)
   call mpibcast(nphysproc,1, mpiint, 0, mpicom)
   call mpibcast(nfld,     1, mpiint, 0, mpicom)
#endif

   if (nmetric==0) then

      if (masterproc) then
         write(iulog,*)'==========================================================='
         write(iulog,*)'       *** Conditional diagnostics NOT requested ***'
         write(iulog,*)'==========================================================='
      end if

      return

   end if

#ifdef SPMD
   !--------------------------------------
   ! Broadcast namelist variables
   !--------------------------------------
   call mpibcast(metric_name,      nmetric_max*len(metric_name(1)), mpichar, 0, mpicom)
   call mpibcast(metric_nver,      nmetric_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_cmpr_type, nmetric_max,                     mpiint,  0, mpicom)
   call mpibcast(metric_threshold, nmetric_max,                     mpir8,   0, mpicom)
   call mpibcast(metric_fillvalue, nmetric_max,                     mpir8,   0, mpicom)
   call mpibcast(sample_after,     nmetric_max*len(sample_after(1)),mpichar, 0, mpicom)

   call mpibcast(ptend_name,    nphysproc_max*len(ptend_name(1)),   mpichar, 0, mpicom)
   call mpibcast(proc_outname,  nphysproc_max*len(proc_outname(1)), mpichar, 0, mpicom)

   call mpibcast(fld_name,  nfld_max*len(fld_name(1)),  mpichar, 0, mpicom)
   call mpibcast(fld_nver,  nfld_max,                   mpiint,  0, mpicom)

   call mpibcast(l_output_state, 1, mpilog, 0, mpicom)
   call mpibcast(l_output_incrm, 1, mpilog, 0, mpicom)
#endif

   !-------------------------------------------
   !  Pack information into to cnd_diag_info
   !-------------------------------------------
   ! metrics for conditional sampling

   cnd_diag_info%nmetric = nmetric

   allocate( cnd_diag_info%metric_name(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_name')
   do ii = 1,nmetric
      cnd_diag_info%metric_name(ii) = trim(adjustl(metric_name(ii)))
   end do

   allocate( cnd_diag_info%metric_nver(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_nver')
   cnd_diag_info%metric_nver(1:nmetric) = metric_nver(1:nmetric)

   allocate( cnd_diag_info%metric_cmpr_type(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_cmpr_type')
   cnd_diag_info%metric_cmpr_type(1:nmetric) = metric_cmpr_type(1:nmetric)

   allocate( cnd_diag_info%metric_threshold(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_threshold')
   cnd_diag_info%metric_threshold(1:nmetric) = metric_threshold(1:nmetric)

   allocate( cnd_diag_info%metric_fillvalue(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%metric_fillvalue')
   cnd_diag_info%metric_fillvalue(1:nmetric) = metric_fillvalue(1:nmetric)

   allocate( cnd_diag_info%sample_after(nmetric), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%sample_after')
   cnd_diag_info%sample_after(1:nmetric) = sample_after(1:nmetric)

   ! atmospheric processes to monitor 

   cnd_diag_info%nphysproc = nphysproc

   allocate( cnd_diag_info%ptend_name(nphysproc), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%ptend_name')
   do ii = 1,nphysproc
      cnd_diag_info%ptend_name(ii) = trim(adjustl(ptend_name(ii)))
   end do

   allocate( cnd_diag_info%proc_outname(nphysproc), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%proc_outname')
   do ii = 1,nphysproc
      cnd_diag_info%proc_outname(ii) = trim(adjustl(proc_outname(ii)))
   end do

   ! snapshots and increments of physical fields

   cnd_diag_info%nfld = nfld

   allocate( cnd_diag_info%fld_name(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%fld_name')
   do ii = 1,nfld
      cnd_diag_info%fld_name(ii) = trim(adjustl(fld_name(ii)))
   end do

   allocate( cnd_diag_info%fld_nver(nfld), stat=ierr)
   if ( ierr /= 0 ) call endrun(subname//': allocation of cnd_diag_info%fld_nver')
   cnd_diag_info%fld_nver(1:nfld) = fld_nver(1:nfld)

   ! output to history file(s)

   cnd_diag_info%l_output_state = l_output_state
   cnd_diag_info%l_output_incrm = l_output_incrm

   !-----------------------------------------------
   ! Send information to log file
   !-----------------------------------------------
   if (masterproc) then

      write(iulog,*)'==========================================================='
      write(iulog,*)'       *** Conditional diagnostics requested ***'
      write(iulog,*)'-----------------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,2x,a10,a6,a12,a20,a20,a20)')'metric','nlev','cmpr_type','threshold','fill_value', 'sample_after'
      do ii = 1,cnd_diag_info%nmetric
         write(iulog,'(i4.3,2x,a10,i6,i12,e20.10,e20.10,a20)') ii,                          &
                                                adjustr(cnd_diag_info%metric_name(ii)),     &
                                                        cnd_diag_info%metric_nver(ii),      &
                                                        cnd_diag_info%metric_cmpr_type(ii), &
                                                        cnd_diag_info%metric_threshold(ii), &
                                                        cnd_diag_info%metric_fillvalue(ii), &
                                                adjustr(cnd_diag_info%sample_after(ii))
      end do

      write(iulog,*)
      write(iulog,'(4x,a20,a20)')'ptend_name', 'proc_outname'
      do ii = 1,cnd_diag_info%nphysproc
         write(iulog,'(i4.3,a20,a20)') ii, adjustr(cnd_diag_info%ptend_name(ii)), adjustr(cnd_diag_info%proc_outname(ii))
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,'(4x,a30,a6)')'physical fields','nlev'
      do ii = 1,cnd_diag_info%nfld
         write(iulog,'(i4.3,a30,i6)') ii, adjustr(cnd_diag_info%fld_name(ii)), cnd_diag_info%fld_nver(ii)
      end do
      write(iulog,*)'--------------------------------------------------'

      write(iulog,*)
      write(iulog,*)' l_output_state = ',l_output_state
      write(iulog,*)' l_output_incrm = ',l_output_incrm
      write(iulog,*)
      write(iulog,*)'==========================================================='
      write(iulog,*)

  end if  ! masterproc

end subroutine conditional_diag_readnl

!===============================================================================
subroutine conditional_diag_alloc( psetcols, metric_nver, nphysproc, nfld, fld_nver, diag )

  use infnan, only : inf, assignment(=)

  integer, intent(in) :: psetcols
  integer, intent(in) :: metric_nver, nphysproc
  integer, intent(in) :: nfld
  integer, intent(in) :: fld_nver(nfld)

  type(cnd_diag_t), intent(inout) :: diag

  integer :: ifld
  integer :: ierr

  character(len=*), parameter :: subname = 'conditional_diag_alloc'

  ! the metric fields, which might have 1, pver, or pver+1 vertical levels

  allocate( diag% metric(psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of diag%metric')

  allocate( diag% flag  (psetcols,metric_nver), stat=ierr)
  if ( ierr /= 0 ) call endrun(subname//': allocation of diag%flag')

  ! diagnostical fields

  if (nfld > 0) then

     allocate( diag% fld( nfld ), stat=ierr)
     if ( ierr /= 0 ) call endrun(subname//': allocation of diag%fld')

     do ifld = 1, nfld  ! snapshots and increments of each field

        allocate( diag%fld(ifld)% old(psetcols,fld_nver(ifld)), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of diag%fld%old')

        allocate( diag%fld(ifld)% val(psetcols,fld_nver(ifld),nphysproc), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of diag%fld%val')

        allocate( diag%fld(ifld)% inc(psetcols,fld_nver(ifld),nphysproc), stat=ierr)
        if ( ierr /= 0 ) call endrun(subname//': allocation of diag%fld%inc')

        diag%fld(ifld)% old(:,:)   = inf
        diag%fld(ifld)% val(:,:,:) = inf
        diag%fld(ifld)% inc(:,:,:) = inf

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
  character(len=4) :: val_inc_suff(2), suff
  logical          :: l_output(2)

  character(len=max_fieldname_len) :: output_fld_name
  character(len=max_fieldname_len) :: output_fld_name2

  character(len=256) :: fld_long_name

  character(len=*),parameter :: subname = 'conditional_diag_output_init'

  if (cnd_diag_info%nmetric==0) return

  do im = 1,cnd_diag_info%nmetric

     !----------------------------------------
     ! register the metric itself for output
     !----------------------------------------
     output_fld_name  = metric_name_in_output( im, cnd_diag_info )
     output_fld_name2 =   flag_name_in_output( im, cnd_diag_info )

     if (cnd_diag_info%metric_nver(im)==1) then

       call addfld(trim(output_fld_name ), horiz_only, 'A',' ',' ') 
       call addfld(trim(output_fld_name2), horiz_only, 'A',' ',' ') 

     elseif(cnd_diag_info%metric_nver(im)==pver) then

       call addfld(trim(output_fld_name ), (/'lev'/),  'A',' ',' ') 
       call addfld(trim(output_fld_name2), horiz_only, 'A',' ',' ') 

     elseif(cnd_diag_info%metric_nver(im)==pver+1) then

       call addfld(trim(output_fld_name ), (/'ilev'/), 'A',' ',' ') 
       call addfld(trim(output_fld_name2), (/'ilev'/), 'A',' ',' ') 

     else 
       call endrun(subname//': invalid number of vertical levels for metric '//trim(cnd_diag_info%metric_name(im)))
     end if 

     call add_default(trim(output_fld_name ),1,' ')
     call add_default(trim(output_fld_name2),1,' ')

     !----------------------------------------
     ! register the diagnostics for output
     !----------------------------------------
     val_inc_suff = (/"_val","_inc"/)
     l_output     = (/cnd_diag_info%l_output_state, cnd_diag_info%l_output_incrm/)

     do ii = 1,2 ! field value (state) or increment

        if (.not.l_output(ii)) cycle
        suff = val_inc_suff(ii)

        do ifld  = 1,cnd_diag_info%nfld
        do iphys = 1,cnd_diag_info%nphysproc

           output_fld_name = fld_name_in_output( im, ifld, iphys, suff, cnd_diag_info)
           fld_long_name   = fld_long_name_in_output( im, ifld, iphys, suff, cnd_diag_info)

           if (cnd_diag_info%fld_nver(ifld)==1) then
              call addfld(trim(output_fld_name), horiz_only, 'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%fld_nver(ifld)==pver) then
              call addfld(trim(output_fld_name), (/'lev'/),  'A',' ',trim(fld_long_name)) 

           elseif (cnd_diag_info%fld_nver(ifld)==pver+1) then
              call addfld(trim(output_fld_name), (/'ilev'/), 'A',' ',trim(fld_long_name)) 
           else
              call endrun(subname//': invalid number of vertical levels for '//cnd_diag_info%fld_name(ifld))
           end if

           call add_default(trim(output_fld_name),1,' ')

        end do ! iphys
        end do ! ifld

     end do ! ii = 1,2, field value (state) or tendency

  end do ! im = 1,nmetric

end subroutine conditional_diag_output_init

!======================================================
function metric_name_in_output( im, cnd_diag_info )

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: im
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=max_fieldname_len) :: metric_name_in_output

   character(len=2) :: imstr ! metric index as a string

   write(imstr,'(i2.2)') im
   metric_name_in_output = 'cnd'//imstr//'_'//trim(cnd_diag_info%metric_name(im))

end function metric_name_in_output

!======================================================
function flag_name_in_output( im, cnd_diag_info )

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: im
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=max_fieldname_len) :: flag_name_in_output

   character(len=2) :: imstr ! metric index as a string

   write(imstr,'(i2.2)') im
   flag_name_in_output = 'cnd'//imstr//'_'//trim(cnd_diag_info%metric_name(im))//'_flag'

end function flag_name_in_output

!======================================================
function fld_name_in_output( im, ifld, iphys, suff, cnd_diag_info )

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: im, ifld, iphys
   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=max_fieldname_len) :: fld_name_in_output 

   character(len=2) :: imstr ! metric index as a string

   write(imstr,'(i2.2)') im
   fld_name_in_output = 'cnd'//imstr//'_'// &
                        trim(cnd_diag_info%fld_name(ifld))//'_'// &
                        trim(cnd_diag_info%proc_outname(iphys))//suff

end function fld_name_in_output

!======================================================
function fld_long_name_in_output( im, ifld, iphys, suff, cnd_diag_info )

   use cam_history_support, only: max_fieldname_len

   integer,               intent(in)  :: im, ifld, iphys
   character(len=*),      intent(in)  :: suff
   type(cnd_diag_info_t), intent(in)  :: cnd_diag_info

   character(len=256) :: fld_long_name_in_output 

   character(len=2) :: imstr ! metric index as a string

   write(imstr,'(i2.2)') im

   fld_long_name_in_output = trim(cnd_diag_info%fld_name(ifld))//suff// &
                             ' sampled under condition '//imstr// &
                             ' ('//trim(cnd_diag_info%metric_name(im))//')' 

end function fld_long_name_in_output


end module conditional_diag

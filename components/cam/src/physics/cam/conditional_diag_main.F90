module conditional_diag_main

  use shr_kind_mod,   only: r8 => shr_kind_r8
  use cam_abortutils, only: endrun

  use physics_types,    only: physics_state
  use camsrfexch,       only: cam_in_t

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none

  private

  public conditional_diag_cal_and_output

  integer, parameter :: GT  =  1
  integer, parameter :: GE  =  2
  integer, parameter :: LT  = -1
  integer, parameter :: LE  = -2

  real(r8),parameter :: ON  = 1._r8
  real(r8),parameter :: OFF = 0._r8

  real(r8),parameter :: FILLVALUE = 0._r8

contains

!======================================================
subroutine conditional_diag_cal_and_output( state, pname, cam_in )

  use ppgrid,              only: pcols
  use cam_history_support, only: max_fieldname_len
  use cam_history,         only: outfld

  use conditional_diag_hist_util, only: metric_name_in_output, &
                                        flag_name_in_output, &
                                        fld_name_in_output

  type(physics_state), intent(inout), target :: state
  character(len=*),    intent(in)    :: pname

  type(cam_in_t),      intent(in),optional :: cam_in

  integer :: nmetric, nphysproc, nfld
  integer :: im, iphys, ii, ifld
  integer :: ncol, lchnk

  real(r8),pointer :: metric(:,:), flag(:,:), new(:,:), inc(:,:), old(:,:)
  real(r8),pointer :: fld(:,:)

  character(len=max_fieldname_len) :: outfldname

  if (cnd_diag_info%nmetric == 0) return  ! no conditional diagnostics requested 

  nmetric   = cnd_diag_info%nmetric
  nphysproc = cnd_diag_info%nphysproc
  nfld      = cnd_diag_info%nfld

  lchnk    = state%lchnk
  ncol     = state%ncol

  !---------------------------------------------------------------------------
  ! Check if this atmospheric process is being monitored 
  !---------------------------------------------------------------------------
  iphys = 0
  do ii = 1,nphysproc
     if ( trim(cnd_diag_info%ptend_name(ii)) == trim(pname) ) then 
        iphys = ii
        exit
     end if
  end do

  if (iphys>0) then 
  ! This atmospheric process is being monitored. Calculate the diagnostics
  ! (and their increments if needed), then save to cnd_diag. Note that here 
  ! we only calculate and save the diagnostics and/or increments. 
  ! Conditional sampling won't be applied until the metrics and flags 
  ! have been obtained.

     do ifld = 1,nfld

        ! Note: The current implementation is such that the same sets of 
        ! atmospheric processes and diagnostics are monitored for all 
        ! different metrics. Here we calculate the diagnostics and/or
        ! increments just once, for im = 1, and then copy the values
        ! for other metrics (im = 2,nmetric). When conditional sampling
        ! is applied, these different copies will likely be sampled different.

        ! Obtain the most up-to-date values of the diagnostic variable 

        im = 1
        new => state%cnd_diag(im)%fld(ifld)% val(:,:,iphys)
        call get_values( trim(cnd_diag_info%fld_name(ifld)), state, new ) !in, in, out

        do im = 2,nmetric
           state%cnd_diag(im)%fld(ifld)% val(1:ncol,:,iphys) = new(1:ncol,:)
        end do

        ! Calculate increments

        if (cnd_diag_info%l_output_incrm) then

           im = 1
           inc => state%cnd_diag(im)%fld(ifld)% inc(:,:,iphys)
           old => state%cnd_diag(im)%fld(ifld)% old

           inc(1:ncol,:) = new(1:ncol,:) - old(1:ncol,:)
           old(1:ncol,:) = new(1:ncol,:)

           ! Save increments for other metrics; update "old" value

           do im = 2,nmetric
              state%cnd_diag(im)%fld(ifld)% inc(1:ncol,:,iphys) = inc(1:ncol,:)
              state%cnd_diag(im)%fld(ifld)% old(1:ncol,:)       = new(1:ncol,:)
           end do
          
        end if ! l_output_incrm

     end do ! ifld = 1,nfld
  end if ! iphys > 0

  !-------------------------------------------------------------------------------
  ! Conditional sampling 
  !-------------------------------------------------------------------------------
  do im = 1,nmetric
     if (trim(pname).eq.trim(cnd_diag_info%sample_after(im))) then

        !----------------------------------------
        ! Get metric values and send to history

        metric => state%cnd_diag(im)%metric
        call get_values( trim(cnd_diag_info%metric_name(im)), state, metric )

        outfldname = metric_name_in_output( im, cnd_diag_info )
        call outfld( trim(outfldname), metric, pcols, lchnk )

        !----------------------------------------
        ! Get flag values and send to history

        flag => state%cnd_diag(im)%flag
        call get_flags( metric, im, ncol, cnd_diag_info, flag )

        outfldname = flag_name_in_output( im, cnd_diag_info )
        call outfld( trim(outfldname), flag, pcols, lchnk )

        !--------------------------------------------
        ! Apply conditional sampling to diagnostics;
        ! send diagnostics to history

        if (cnd_diag_info%l_output_state) then        
        do ifld = 1,nfld
        do ii   = 1,nphysproc

           fld => state%cnd_diag(im)%fld(ifld)% val(:,:,ii)
           where(flag.eq.OFF)  fld = FILLVALUE 

           outfldname = fld_name_in_output( im, ifld, ii, '_val', cnd_diag_info)
           call outfld( trim(outfldname), fld, pcols, lchnk)
        end do
        end do
        end if

        if (cnd_diag_info%l_output_incrm) then        
        do ifld = 1,nfld
        do ii   = 1,nphysproc

           fld => state%cnd_diag(im)%fld(ifld)% inc(:,:,ii)
           where(flag.eq.OFF)  fld =  FILLVALUE

           outfldname = fld_name_in_output( im, ifld, ii, '_inc', cnd_diag_info)
           call outfld( trim(outfldname), fld, pcols, lchnk)
        end do
        end do
        end if

     end if  !trim(pname).eq.trim(cnd_diag_info%sample_after(im))
  end do ! im = 1,nmetric
     

end subroutine conditional_diag_cal_and_output

!========================================================
subroutine get_values( varname, state, arrayout )

  character(len=*),   intent(in)    :: varname
  type(physics_state),intent(inout) :: state
  real(r8),           intent(inout) :: arrayout(:,:)

  character(len=*),parameter :: subname = 'get_values'

  integer :: ncol

  ncol = state%ncol

  select case (trim(adjustl(varname)))
  case('T')
     arrayout(1:ncol,:) = state%t(1:ncol,:)

  case('U')
     arrayout(1:ncol,:) = state%u(1:ncol,:)

  case('V')
     arrayout(1:ncol,:) = state%v(1:ncol,:)

  case('PMID')
     arrayout(1:ncol,:) = state%pmid(1:ncol,:)

  case('PINT')
     arrayout(1:ncol,:) = state%pmid(1:ncol,:)

  case('PS')
     arrayout(1:ncol,1) = state%ps(1:ncol)

 !elseif (varname.eq.'QSATW') then
 !   call qsatw()

 !elseif (varname.eq.'CAPE') then
 !   call cape()

  case default 
     call endrun(subname//': unknow varname - '//trim(varname))
  end select

end subroutine get_values

subroutine get_flags( metric, im, ncol, cnd_diag_info, flag )

  real(r8),              intent(in) :: metric(:,:)
  integer,               intent(in) :: im
  integer,               intent(in) :: ncol
  type(cnd_diag_info_t), intent(in) :: cnd_diag_info

  real(r8), intent(inout) :: flag(:,:)

  character(len=*),parameter :: subname = 'get_flags'

  flag(:,:) = OFF

  select case (cnd_diag_info%metric_cmpr_type(im))
  case (GT)

    where (metric(1:ncol,:).gt.cnd_diag_info%metric_threshold(im))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (GE)

    where (metric(1:ncol,:).ge.cnd_diag_info%metric_threshold(im))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LT)

    where (metric(1:ncol,:).lt.cnd_diag_info%metric_threshold(im))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case (LE)

    where (metric(1:ncol,:).le.cnd_diag_info%metric_threshold(im))
      flag(1:ncol,:) = ON
    elsewhere
      flag(1:ncol,:) = OFF
    end where

  case default
    call endrun(subname//': unknown cnd_diag_info%metric_cmpr_type')
  end select

end subroutine get_flags

end module conditional_diag_main

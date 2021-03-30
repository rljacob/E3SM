module conditional_diag_hist_util
!-------------------------------------------------
! Conditional diagnostics.
! History:
!  First version by Hui Wan, PNNL, March 2021
!-------------------------------------------------
  use cam_abortutils, only: endrun

  use conditional_diag, only: cnd_diag_info, cnd_diag_info_t

  implicit none
  public

contains

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

       call addfld(trim(output_fld_name ), horiz_only, 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), horiz_only, 'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(im)==pver) then

       call addfld(trim(output_fld_name ), (/'lev'/),  'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'lev'/),  'A',' ','Flags  used in conditional sampling')

     elseif(cnd_diag_info%metric_nver(im)==pver+1) then

       call addfld(trim(output_fld_name ), (/'ilev'/), 'A',' ','Metric used in conditional sampling')
       call addfld(trim(output_fld_name2), (/'ilev'/), 'A',' ','Flags  used in conditional sampling')

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
                             ' after '//trim(cnd_diag_info%ptend_name(iphys))// &
                             ' sampled under condition '//imstr// &
                             ' ('//trim(cnd_diag_info%metric_name(im))//')' 

end function fld_long_name_in_output


end module conditional_diag_hist_util

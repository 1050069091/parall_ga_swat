subroutine parall_pre_calibration
use mymodule
use mpi

integer :: allocate_err1=0,allocate_err2=0,allocate_err3=0,allocate_err4=0,lentrim=1,all_lentrim
integer :: tmp_int,ierr,tmp_num=0

lentrim = 1
do i=1,file_kind_rank
    lentrim = lentrim + len_trim(detail_each_file_include_para_arr(i)%postfile_name)
    all_lentrim = len_trim(swat_input_file_postfix_names)
    !write(*,*) 'wwwww:',lentrim,all_lentrim
    swat_input_file_postfix_names(all_lentrim+1:all_lentrim+lentrim) = detail_each_file_include_para_arr(i)%postfile_name(1:lentrim)
    !write(*,*) 'wwwww:',swat_input_file_postfix_names, trim(detail_each_file_include_para_arr(i)%postfile_name)

end do


files_need_modi=0
do i=1,para_rank
    files_need_modi = para_start_value_3arr(myid+1,i)%size + files_need_modi
end do

!write(*,*) "mod(seed_num,numprocs)",mod(seed_rank,numprocs),myid,interval
if((mod(seed_rank,numprocs)/=0).and.(myid+1 > mod(seed_rank,numprocs))) then
    my_interval = interval - 1
else
    my_interval = interval
end if

allocate(para_values_2arr_in_one_process(interval,files_need_modi),stat=allocate_err1)
allocate(obj_fn_valuses_in_one_process(interval),stat=allocate_err2)
allocate(seed_nums_in_one_process(interval),stat=allocate_err3)

if(allocate_err1 + allocate_err2 + allocate_err3 /= 0) then
    write(*,*) 'in parall_pre_calibration.f90:22,allocate err!'
  stop
end if


end subroutine parall_pre_calibration

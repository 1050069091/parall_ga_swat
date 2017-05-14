subroutine pre_calibration_handle
use mymodule
      
call date_and_time(a(1), a(2), a(3), date_time1)

!获取要率定的参数信息
call read_para_info()
!如果进程数大于种子数，无并行的必要，报错停止
if(seed_rank < numprocs) then
	write(*,"('error: process_num(',I4,' ) > seed_num(',I4,'  )')") numprocs,seed_rank
	stop 
end if

interval = seed_rank / numprocs
if(seed_rank >  numprocs*interval) then
    interval = interval + 1
end if

call format_para_fn(len(myid_str),myid+1,myid_str)

end subroutine pre_calibration_handle

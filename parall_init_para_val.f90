subroutine parall_init_para_value()
use mymodule

integer :: seed_num,allocate_err = 0

allocate(para_start_value_3arr(seed_rank,para_rank),stat=allocate_err)
if(allocate_err/= 0) then
    write(*,*) 'allocate para_start_value_arr3 error'
	stop
end if

!随机初始化种群中个体的位置和速度
call init_random_seed(myid+1)
do i=0,interval-1
	seed_num = myid+1 + numprocs*i
	if(seed_num > seed_rank) then
		exit
	end if
	!初始化需要率定参数的值：在各参数范围内随机生成
	call initialize_para_val(seed_num)

end do

allocate(files_included_by_each_para_arr(para_rank),stat=allocate_err)
if(allocate_err/= 0) then
    write(*,*) 'allocate files_included_by_each_para_arr error'
	stop
end if

do i=1,para_rank

	if(sub_or_hru_arr(i) == 0) then !该参数是每sub就一个
    	files_included_by_each_para_arr(i) = substream_rank
    else 
    	files_included_by_each_para_arr(i) = hru_rank
    end if

end do

end subroutine parall_init_para_value

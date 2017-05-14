subroutine initialize_para_val(seed_num)
use mymodule
integer,intent(in) :: seed_num
integer :: allocate_err = 0,para_file_rank = 0
real,dimension(:),pointer :: tmp_place_p




!随机初始化种群中个体的位置和速度
!call init_random_seed()

! do k=1,seed_rank
	do i=1,para_rank

		if(sub_or_hru_arr(i) == 0) then !该参数是每sub就一个
    		para_file_rank = substream_rank
    	else 
    		para_file_rank = hru_rank
    	end if

    	allocate(tmp_place_p(para_file_rank),stat=allocate_err)
    	if(allocate_err/= 0) then
    		write(*,*) 'allocate tmp_place_p error'
			stop
		end if

    	do j=1,para_file_rank
           call random_number(tmp_place_p(j))   !随机化每个文件中的对应参数值
           tmp_place_p(j) = min_para_arr(i) + tmp_place_p(j)*(max_para_arr(i)-min_para_arr(i))
    	end do
        para_start_value_3arr(seed_num,i)%p => tmp_place_p
        para_start_value_3arr(seed_num,i)%size = para_file_rank

	end do
! end do

end subroutine initialize_para_val

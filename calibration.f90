subroutine calibration()
use mymodule

integer :: seed_num,int_j,int_index
real :: fn_obj_val = 0.0
real,dimension(observed_data_rank) :: simu_datas_arr
real :: best_obj_fn_val_in_one_process = -99999999.9


outter:do i=1,recur_rank
    int_j = 0
    inner:do 

        seed_num = myid+1 + int_j * numprocs
        if(seed_num > seed_rank) exit inner
        if(.not.(i == 1 .and. seed_num == 1 .and. is_use_befor == 1)) then
		    call modify_para_in_inputfile(seed_num)
            !write(*,*) modify,myid
        end if
        
        call system('./swat2012_627_paral '//myid_str//' '//trim(swat_input_file_postfix_names)//' > /dev/null')
		!获得模拟值，并写入文件中
		call extract_swat_output(simu_datas_arr)

		if(fn_type == 1) then
			call solve_r2(observed_data_val_arr,simu_datas_arr,observed_data_rank,fn_obj_val)
            !write(*,"(' recur num:',I4,'    seed num:',I4,A10,A6,F10.4,'    <-------- by node',I3)") &
            !                                i,seed_num,'------->','  R2=',fn_obj_val,myid
		else
			call solve_ns(observed_data_val_arr,simu_datas_arr,observed_data_rank,fn_obj_val)
            !write(*,"(' recur num:',I4,'    seed num:',I4,A10,A6,F10.4,'    <-------- by node',I3)") &
            !                                i,seed_num,'------->','  NS=',fn_obj_val,myid
		end if
		
        if(fn_obj_val*10000 > best_obj_fn_val_in_one_process*10000) then
            best_obj_fn_val_in_one_process = fn_obj_val

            call system('mv -f out/'//trim(obser_val_name)&
                //myid_str//'.txt out/best_sim_data.txt'//trim(myid_str))
            call system("cp -f `ls | grep '[.]"//myid_str//"$'` ./out")
        end if

        int_j = int_j + 1
        obj_fn_valuses_in_one_process(int_j) = fn_obj_val
        seed_nums_in_one_process(int_j) = seed_num

        int_index = 0
        do m=1,para_rank
            do n=1,para_start_value_3arr(seed_num,m)%size
                para_values_2arr_in_one_process(int_j,int_index+n) = para_start_value_3arr(seed_num,m)%p(n)
            end do
            int_index = int_index + para_start_value_3arr(seed_num,m)%size
        end do
    end do inner
    !开始通信
    call communication_to_get_best_values(i)

end do outter

end subroutine calibration

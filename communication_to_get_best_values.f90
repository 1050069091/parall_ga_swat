subroutine communication_to_get_best_values(recur_num)
use mymodule
use mpi
implicit none

integer,intent(in) :: recur_num

integer ierr,recv_buff_size,recv_myid,int_index,seed_num,int_j
real recv_obj_fn_values(interval*numprocs)

real send_para_values(interval*files_need_modi),rate_obj_fn_val_arr(seed_rank)
real recv_para_values(numprocs*interval*files_need_modi)

integer int_obj_fn_val(seed_rank)&
    ,plus_int_obj_fn_val(seed_rank),seed_nums_sorted(seed_rank)

real new_all_para_values_arr(numprocs*interval*files_need_modi)

real random_real_val
integer parent_index1,parent_index2,myid1,myid2,interval1,interval2&
    ,count_ga,intersect_para_index,intersect_index1,intersect_index2

real para_val_in_one_seed_arr1(files_need_modi),para_val_in_one_seed_arr2(files_need_modi),tmp_real

integer seed_count,tmp_interval,tmp_myid,i,j,m,n,inter_time

recv_buff_size=interval*files_need_modi
recv_myid=numprocs-1

do i=1,interval
   send_para_values((i-1)*files_need_modi+1:i*files_need_modi) = &
                        para_values_2arr_in_one_process(i,1:files_need_modi)

end do

!聚集各个体目标函数值到最后进程
call MPI_GATHER(obj_fn_valuses_in_one_process,interval,MPI_REAL&
                ,recv_obj_fn_values,interval,MPI_REAL&
                ,recv_myid,mycomm,ierr)
!聚集各个体染色体上的基因到目标进程
call MPI_GATHER(send_para_values,recv_buff_size&
                ,MPI_REAL,recv_para_values,recv_buff_size&
                ,MPI_REAL,recv_myid,mycomm,ierr)

if(myid==recv_myid) then
    int_index = 1
    int_j = 1
    do i=1,numprocs
        do j=1,interval
            seed_num = (j-1)*numprocs+i 

            if(seed_num <= seed_rank) then 
                seed_nums_sorted(int_j) = seed_num
                int_obj_fn_val(int_j) = 10000*recv_obj_fn_values(int_index) 
                int_j = int_j + 1
            end if
            int_index = int_index + 1
        end do
    end do

    call quick_sort(int_obj_fn_val,seed_nums_sorted,1,seed_rank)

     call date_and_time(b(1), b(2), b(3), date_time2)
     inter_time = date_time2(7)-date_time1(7) + 60*(date_time2(6)-date_time1(6)) &
            + 60*60*(date_time2(5)-date_time1(5)) + 60*60*60*(date_time2(4)-date_time1(4))&
            + 60*60*60*12*(date_time2(3)-date_time1(3))
     write(*,"(' ***************************遗传算法进化代数:',I4,&
        ' 最大目标函数值:',F8.4,'******',' 耗时:',I6,'s*****************************************')") &
        recur_num,(int_obj_fn_val(seed_rank))/10000.0,inter_time

    if(int_obj_fn_val(seed_rank) > g_best_int_obj_fn_val) then

        g_best_int_obj_fn_val = int_obj_fn_val(seed_rank)
        best_sim_myid = mod(seed_nums_sorted(seed_rank)-1,numprocs)
        call format_para_fn(5,best_sim_myid+1,best_sim_myid_str)

        best_sim_recur_id = recur_num
        best_sim_seed_id = seed_nums_sorted(seed_rank)

    end if

    plus_int_obj_fn_val(1) = 0
    do i=1,seed_rank
        if(i==1) then
            plus_int_obj_fn_val(1) = int_obj_fn_val(1)
        else
            plus_int_obj_fn_val(i) = int_obj_fn_val(i) + plus_int_obj_fn_val(i-1)
        end if
    end do

    do i=1,seed_rank
        rate_obj_fn_val_arr(i) = plus_int_obj_fn_val(i) / (plus_int_obj_fn_val(seed_rank)*1.0)
    end do

    !write(*,*)int_obj_fn_val
    !write(*,*)plus_int_obj_fn_val
    !write(*,*)rate_obj_fn_val_arr

    !call init_random_seed(numprocs) 
    seed_count = 1
    outter:do 

        !选择两条染色体
        call random_number(random_real_val)  
        random_real_val = 0.85 + 0.15 * random_real_val
        call find_index_choosed(random_real_val,rate_obj_fn_val_arr,seed_rank,parent_index1)
        !write(*,*) int_obj_fn_val(parent_index1),parent_index1
        !inner12:do 
            call random_number(random_real_val)  
        random_real_val = 0.85 + 0.15 * random_real_val
            call find_index_choosed(random_real_val,rate_obj_fn_val_arr,seed_rank,parent_index2)
        !write(*,*) int_obj_fn_val(parent_index2),parent_index2
        !    if(parent_index2 /= parent_index1) exit inner12

        !end do inner12
        !write(*,*) random_real_val,rate_obj_fn_val_arr,seed_rank,parent_index2

        myid1 = mod(seed_nums_sorted(parent_index1)-1,numprocs)
        interval1 = floor((seed_nums_sorted(parent_index1)-1)/(numprocs*1.0))
        myid2 = mod(seed_nums_sorted(parent_index2)-1,numprocs)
        interval2 = floor((seed_nums_sorted(parent_index2)-1)/(numprocs*1.0))
        
        intersect_index1 = myid1*interval*files_need_modi+interval1*files_need_modi
        intersect_index2 = myid2*interval*files_need_modi+interval2*files_need_modi

        !write(*,*) intersect_index1+1,intersect_index1+files_need_modi
        !write(*,*) intersect_index2+1,intersect_index2+files_need_modi

        para_val_in_one_seed_arr1 =&
            recv_para_values(intersect_index1+1:intersect_index1+files_need_modi)
        para_val_in_one_seed_arr2 =&
            recv_para_values(intersect_index2+1:intersect_index2+files_need_modi)

        !交叉操作
        call random_number(random_real_val)
        if(random_real_val <= intersect_rate) then
            count_ga = files_need_modi * intersect_rate
            do i=1,count_ga
               call random_number(random_real_val)
               intersect_para_index = ceiling(files_need_modi * random_real_val)

               tmp_real = para_val_in_one_seed_arr1(intersect_para_index)
               para_val_in_one_seed_arr1(intersect_para_index) = para_val_in_one_seed_arr2(intersect_para_index)
               para_val_in_one_seed_arr2(intersect_para_index) = tmp_real

            end do
        end if

        !变异操作
        call varition_para_values(para_val_in_one_seed_arr1,files_need_modi,recur_num)
        call varition_para_values(para_val_in_one_seed_arr2,files_need_modi,recur_num)

        !把新个体（染色体）添加到新种群中
        tmp_myid = mod(seed_count-1,numprocs)
        tmp_interval = floor((seed_count-1)/(numprocs*1.0))
        
        !write(*,*) tmp_myid*interval*files_need_modi+tmp_interval*files_need_modi+1&
        !    ,tmp_myid*interval*files_need_modi+(tmp_interval+1)*files_need_modi
        
        new_all_para_values_arr(tmp_myid*interval*files_need_modi+tmp_interval*files_need_modi+1&
            :tmp_myid*interval*files_need_modi+(tmp_interval+1)*files_need_modi)&
            = para_val_in_one_seed_arr1(1:files_need_modi)
        seed_count = seed_count + 1
        if(seed_count > seed_rank) exit outter

        tmp_myid = mod(seed_count-1,numprocs)
        tmp_interval = floor((seed_count-1)/(numprocs*1.0))
        
        !write(*,*) tmp_myid*interval*files_need_modi+tmp_interval*files_need_modi+1&
        !    ,tmp_myid*interval*files_need_modi+(tmp_interval+1)*files_need_modi

        new_all_para_values_arr(tmp_myid*interval*files_need_modi+tmp_interval*files_need_modi+1&
            :tmp_myid*interval*files_need_modi+(tmp_interval+1)*files_need_modi)&
            = para_val_in_one_seed_arr2(1:files_need_modi)
        seed_count = seed_count + 1
        if(seed_count > seed_rank) exit outter
        
    end do outter
    !write(*,*) max_para_arr
    !write(*,*) min_para_arr
    !write(*,*) "-------------------->new:",new_all_para_values_arr
    !write(*,*) "-------------------->old:",recv_para_values
end if

!write(*,*) recur_num,myid,"older----:",recv_para_values

call MPI_SCATTER(new_all_para_values_arr,recv_buff_size&
                ,MPI_REAL,recv_para_values,recv_buff_size&
                ,MPI_REAL,numprocs-1,mycomm,ierr)

!write(*,*) recur_num,myid,"later----:",recv_para_values
!write(*,*) myid,'---->',para_values_2arr_in_one_process
do i=1,interval
    para_values_2arr_in_one_process(i,1:files_need_modi) = &
        recv_para_values((i-1)*files_need_modi+1:i*files_need_modi)
end do
!write(*,*) myid,'---->',para_values_2arr_in_one_process

do i=1,my_interval
    int_index = 0
    seed_num = seed_nums_in_one_process(i)
    do m=1,para_rank
        !write(*,*) recur_num,myid,seed_num,"old---->",m,para_start_value_3arr(seed_num,m)%p
        do n=1,para_start_value_3arr(seed_num,m)%size
!            if(para_values_2arr_in_one_process(i,int_index+n) < max_para_arr(m) &
!                .and. para_values_2arr_in_one_process(i,int_index+n) > min_para_arr(m)) then
                para_start_value_3arr(seed_num,m)%p(n) = para_values_2arr_in_one_process(i,int_index+n) 
!            else
!                 call random_number(random_real_val)
!                para_start_value_3arr(seed_num,m)%p(n) = min_para_arr(m) + random_real_val*(max_para_arr(m)-min_para_arr(m))
!            end if
        end do
        !write(*,*) recur_num,myid,seed_num,"new---->",m,para_start_value_3arr(seed_num,m)%p
        int_index = int_index + para_start_value_3arr(seed_num,m)%size
    end do
end do

end subroutine communication_to_get_best_values

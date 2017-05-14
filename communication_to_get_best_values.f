subroutine communication_to_get_best_values()
use mymodule
use mpi
integer :: buff_position,ierr
real tmp_arr(files_need_modi)
character(1) pack_buff_recv(buff_size)

!allocate(tmp_arr(files_need_modi),stat=ierr)

buff_position = 1
if(myid /= numprocs-1) then

    !write(*,*) myid,":seed_nums_in_one_process:",seed_nums_in_one_process
    call MPI_PACK(seed_nums_in_one_process,interval&
        ,MPI_INTEGER,buff_pack,buff_size,buff_position,mycomm,ierr)

    write(*,*) buff_position,myid,":obj_fn_valuses_in_one_provess:",obj_fn_valuses_in_one_process
    call MPI_PACK(obj_fn_valuses_in_one_process,interval&
        ,MPI_REAL,buff_pack,buff_size,buff_position,mycomm,ierr)

    write(*,*) buff_positio,nmyid,":obj_fn_valuses_in_one_provess:",obj_fn_valuses_in_one_process
    do i=1,interval
        tmp_arr(1:files_need_modi) = para_values_2arr_in_one_process(i,1:files_need_modi)
        write(*,*) buff_size,files_need_modi
        call MPI_PACK(tmp_arr,files_need_modi,&
            MPI_REAL,buff_pack,buff_size,buff_position,mycomm,ierr)
    end do

    call MPI_SEND(pack_buff,buff_position-1,MPI_PACKED,numprocs-1,1,mycomm,ierr)
else

    do i=0,numprocs-2

        write(*,*) myid,":::::",buff_size

        call MPI_RECV(pack_buff,buff_size,MPI_PACKED,i,1,mycomm,ierr)

        call MPI_UNPACK(buff_pack,buff_size,buff_position,seed_nums_in_one_process&
            ,interval,MPI_INTEGER,mycomm,ierr)
        write(*,*) myid,":seed_nums_in_one_process:",seed_nums_in_one_process

        call MPI_UNPACK(buff_pack,buff_size,buff_position,obj_fn_valuses_in_one_process&
            ,interval,MPI_REAL,mycomm,ierr)
        write(*,*) myid,":obj_fn_valuses_in_one_provess:",obj_fn_valuses_in_one_process

        do j=1,interval
            call MPI_UNPACK(buff_pack,buff_size,buff_position,tmp_arr&
                ,files_need_modi,MPI_REAL,mycomm,ierr)
            write(*,*) myid,":para_values_2arr_in_one_process(j):",tmp_arr
        end do

    end do


end if


end subroutine communication_to_get_best_values

program test_fchl_read_job

    use mpi
    use ffchl_wrapper, only: kernel_wrapper_fchl
    use ffchl_reader, only: fread_fchl_sizes, fread_fchl_collection

    implicit none

    ! fchl args
    integer :: nm1
    integer :: nm2
    integer :: max_size
    integer :: max_neighbors

    ! end fchl args

    ! fchl mpi
    double precision, dimension(:,:), allocatable :: collection_x
    integer :: collection_size
    ! end fchl mpi


    ! Read the sizes
    call fread_fchl_sizes("jobname/", &
        max_size, &
        max_neighbors, &
        nm1, &
        nm2)

    write(*,*) "max_neighbors", max_neighbors
    write(*,*) "max_size", max_size
    write(*,*) "nm1", nm1

    !
    collection_size = max_size*5*max_neighbors + 1 + max_size
    !

    ! allocate collections
    allocate(collection_x( nm1, collection_size ))
    call fread_fchl_collection("jobname/_fchl__a", collection_x, nm1, max_size, max_neighbors)

    write(*,*) collection_x

    write(*,*) "collection * size", collection_size*nm1

    deallocate(collection_x)

end program test_fchl_read_job


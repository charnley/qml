program test_fchl_gatherv_kernel

    use mpi
    use ffchl_wrapper, only: kernel_wrapper_fchl
    use ffchl_reader, only: fread_fchl_collection
    use freaders, only: fread_1d_integer, fread_2d_integer
    use fprinters, only: print_matrix

    implicit none

    integer :: i
    integer :: ratio, left
    integer :: fr, to

    integer :: local_size


    ! fchl args
    integer :: nm1
    integer :: nm2
    integer :: max_size
    integer :: max_neighbors

    integer, dimension(:), allocatable :: n1
    integer, dimension(:), allocatable :: n2
    integer, dimension(:,:), allocatable :: nneigh1
    integer, dimension(:,:), allocatable :: nneigh2

    double precision, dimension(:,:,:), allocatable :: kernels
    double precision, dimension(:,:), allocatable :: kernel
    double precision, dimension(:,:), allocatable :: local_kernel

    integer :: nsigmas
    ! end fchl args

    ! fchl mpi
    double precision, dimension(:,:), allocatable :: collection_x
    double precision, dimension(:,:), allocatable :: collection_y
    double precision, dimension(:,:), allocatable :: local_collection_x
    double precision, dimension(:,:), allocatable :: local_collection_y

    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end fchl mpi

    ! mpi idx
    integer :: irank
    integer :: ierror
    integer :: procs

    integer :: groupsize
	integer :: displacement
	integer :: scount
    integer, dimension(:), allocatable :: displacements
    integer, dimension(:), allocatable :: rcounts

    integer :: local_ksize_x
    integer :: local_ksize_y
    ! end mpi idx

    integer :: K


    ! hello mpi
    call MPI_Init(ierror)
    call MPI_Comm_rank(MPI_COMM_WORLD, irank, ierror)
    call MPI_Comm_size(MPI_COMM_WORLD, procs, ierror)

    !
    nm1 = 11
    nm1 = 11
    nm2 = nm1
    max_size = 23
    max_neighbors = 23
    nsigmas = 1
    !

    ! Collection dimensions
    collection_size = max_size*5*max_neighbors + 1 + max_size
    idx_x = 1
    idx_n = idx_x + max_size*5*max_neighbors
    idx_nneigh = idx_n + 1
    !


    ! allocate collections
    allocate(collection_x( nm1, collection_size ))
    allocate(collection_y( nm2, collection_size ))
    !


    if(irank.eq.0) then

        ! read representations
        call fread_fchl_collection("jobname_dev/_fchl__a", collection_x, nm1, max_size, max_neighbors)
        call fread_fchl_collection("jobname_dev/_fchl__b", collection_y, nm2, max_size, max_neighbors)

    end if


    ! Distribute collections
    call MPI_Bcast(collection_x, nm1*collection_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)
    call MPI_Bcast(collection_y, nm2*collection_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierror)


    ! TODO JCK MPI ADMINSTRATION
    ratio = nm1/procs
    left = nm1 - procs*ratio
    fr = 1 + ratio*irank
    to = fr + ratio -1

    if(irank.eq.procs-1) then
        to = to + left
    end if


    ! local collection size
    local_size = to-fr+1


    ! allocate local collections
    allocate(local_collection_y(local_size, collection_size))
    allocate(local_collection_x(nm2, collection_size))
    local_collection_y = collection_y(fr:to,:)
    local_collection_x = collection_x

    local_ksize_x = nm1
    local_ksize_y = local_size


    ! allocate local kernel
    allocate(kernels(nsigmas, local_ksize_x, local_ksize_y))
    allocate(local_kernel(local_ksize_x, local_ksize_y))
    kernels = 0.0d0
	local_kernel = 0.0d0

    ! Fill out kernel
    call kernel_wrapper_fchl(local_collection_x, local_collection_y, kernels)
	local_kernel = kernels(1,:,:)

    ! Collect kernel
    if ( irank .eq. 0 ) then
        allocate(kernel(nm1, nm2))
		kernel = 0.0d0
    end if

    ! Calculate gatherv buffer sizes
	scount = local_size*nm2
    groupsize = procs
    allocate(rcounts(groupsize))
    allocate(displacements(groupsize))

    do i=1,groupsize
        rcounts(i) = ratio*nm2
        displacements(i) = (i-1)*ratio*nm2
    end do

    rcounts(size(rcounts)) = ratio*nm2 + left*nm2

    ! Gather distributed kernel
    call MPI_Gatherv( &
		local_kernel, scount, MPI_DOUBLE_PRECISION, &
		kernel, rcounts, displacements, MPI_DOUBLE_PRECISION, &
		0, MPI_COMM_WORLD, ierror)


    ! Print kernel for testing
    if ( irank .eq. 0) then

        open(unit=7, file="log/test_gatherv_kernel")

        do i=1,nm1
            write(7,"(10000f)") kernel(i, :)
        end do

        close(7)

    end if


    ! Deallocate mem
    deallocate(collection_x)
    deallocate(collection_y)
    deallocate(kernels)
    deallocate(local_kernel)

    if ( irank .eq. 0) then
        deallocate(kernel)
    end if


    ! goodbye mpi
    call MPI_Finalize(ierror)

end program test_fchl_gatherv_kernel


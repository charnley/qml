
program test_fchl_gatherv_kernel

    use mpi

    use freaders, only: fread_2d_double
    use fprinters, only: print_matrix


    use ffchl_wrapper, only: kernel_wrapper_fchl
    use ffchl_reader, only: fread_fchl_sizes, fread_fchl_collection

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
    double precision, dimension(:,:,:), allocatable :: local_kernels
    double precision, dimension(:,:), allocatable :: kernel
    double precision, dimension(:,:), allocatable :: local_kernel

    integer :: nsigmas
    ! end fchl args

    ! prop args
    integer :: n_properties
    double precision, allocatable, dimension(:,:) :: properties
    ! end prop args

    ! arg result
    double precision, allocatable, dimension(:, :) :: alphas


    ! end arg result


    ! fchl mpi
    double precision, dimension(:,:), allocatable :: collection_x
    double precision, dimension(:,:), allocatable :: collection_y
    double precision, dimension(:,:), allocatable :: local_collection_x
    double precision, dimension(:,:), allocatable :: local_collection_y
    integer, dimension(:), allocatable :: local_view_x
    integer, dimension(:), allocatable :: local_view_y

    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end fchl mpi

    ! mpi
    ! pblacs data structures
    integer :: info
    integer :: na, numroc
    integer :: local_i, local_j, global_i, global_j

    integer :: lwork
    double precision, dimension(:), allocatable :: work
    double precision, dimension(:,:), allocatable :: local_K, local_B

    integer :: local_K_rows, local_K_cols
    integer :: local_B_rows, local_B_cols
    integer :: local_id, num_ranks, ranks_rows, ranks_cols, context
    integer :: local_rank_col, local_rank_row, block_size
    integer :: ierr
    integer, dimension(9) :: desca, descb
    integer, dimension(2) :: dims
    double precision, dimension(2) :: work_query
    ! end mpi idx


    ! we have everything we need, now, lets go


    ! hello mpi
    ! init pblacs and set cartesian grid of MPI ranks
    call blacs_pinfo(local_id, num_ranks)

    dims = 0
    call MPI_Dims_create(num_ranks, 2, dims, ierr)
    ranks_rows = dims(1)
    ranks_cols = dims(2)

    block_size = 100

    ! create BLACS context
    call blacs_get(0, 0, context)
    call blacs_gridinit(context, 'R', ranks_rows, ranks_cols)
    call blacs_gridinfo(context, ranks_rows, ranks_cols, &
        local_rank_row, local_rank_col)


    ! write(*,*) local_rank_row, local_rank_col, ranks_rows, ranks_cols, context

    ! Hardcoded test
    nsigmas = 1
    n_properties = 2


    ! read the sizes
    if (local_id .eq. 0) then

        ! Read the sizes
        call fread_fchl_sizes("jobname_dev/", &
            max_size, &
            max_neighbors, &
            nm1, &
            nm2)

    endif

    ! Distribute args
    call MPI_Bcast(nm1, sizeof(nm1), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(nm2, sizeof(nm1), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(n_properties, sizeof(n_properties), MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)

    ! Allocate
    collection_size = max_size*5*max_neighbors + 1 + max_size
    allocate(collection_x( nm1, collection_size ))
    allocate(collection_y( nm2, collection_size ))
    allocate(properties(nm1, n_properties))
    allocate(alphas(nm1, n_properties))

    ! READ (on the first rank)
    if (local_id .eq. 0) then

        ! Collection dimensions and allocate collections, then read collections
        call fread_fchl_collection("jobname_dev/_fchl__a", collection_x, nm1, max_size, max_neighbors)
        call fread_fchl_collection("jobname_dev/_fchl__b", collection_y, nm2, max_size, max_neighbors)

        ! Allocate and read properties
        call fread_2d_double("jobname_dev/_properties", nm1, 1, properties(:,1:1))
        call fread_2d_double("jobname_dev/_properties", nm1, 1, properties(:,2:2))

    endif ! rank 0
    ! END READ


    ! Distribute collections
    call MPI_Bcast(collection_x, nm1*collection_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(collection_y, nm2*collection_size, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
    call MPI_Bcast(properties, nm1*n_properties, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)

    ! Allocate kernel and output
    local_K_rows = numroc(nm1, block_size, local_rank_row, 0, ranks_rows)
    local_K_cols = numroc(nm2, block_size, local_rank_col, 0, ranks_cols)

    allocate(local_K(local_K_rows, local_K_cols))
    allocate(local_kernels(1, local_K_rows, local_K_cols))

    local_B_rows = numroc(nm1, block_size, local_rank_row, 0, ranks_rows)
    local_B_cols = numroc(n_properties, block_size, local_rank_col, 0, ranks_cols)
    allocate(local_B(local_B_rows, local_B_cols))


    ! TODO calculate local collection

    ! write(*,*) "local collection", local_K_rows, collection_size
    allocate(local_collection_x(local_K_rows, collection_size))
    allocate(local_collection_y(local_K_cols, collection_size))

    allocate(local_view_x(local_K_rows*local_K_cols))
    allocate(local_view_y(local_K_rows*local_K_cols))

    ! cycklic block distribution pblas

    do local_i = 1, local_K_rows
        call l2g(local_i, local_rank_row, ranks_rows, &
            block_size, global_i)
        local_collection_x(local_i,:) = collection_x(global_i,:)
    end do

    do local_j = 1, local_K_cols
        call l2g(local_j, local_rank_col, ranks_cols, &
            block_size, global_j)
        local_collection_y(local_j,:) = collection_x(global_j,:)
    end do

    ! compute local kernel
    call kernel_wrapper_fchl(local_collection_x, local_collection_y, local_kernels)
	local_K = local_kernels(1,:,:)

    ! Solve KRR
    if (local_id .eq. 0) then

        ! call print_matrix(local_K, 11, 11)

    end if


    ! Setup variables for LAPACK
    call descinit(desca, nm1, nm2, block_size, block_size, 0, 0, context, MAX(1, local_K_rows), info)
    call descinit(descb, nm1, n_properties, block_size, block_size, 0, 0, context, MAX(1, local_B_rows), info)

    ! Allocate local work arrays
    call pdgels("N", nm1, nm2, n_properties, local_K, 1, 1, desca, local_B, 1, 1, DESCB, work_query, -1, info)
    lwork = INT(work_query(1))
    allocate(work(lwork))

    ! copy data
    local_B = 0.0d0
    do local_j = 1, local_B_cols
        do local_i = 1, local_B_rows
            call l2g(local_i, local_rank_row, ranks_rows, block_size, global_i)
            call l2g(local_j, local_rank_col, ranks_cols, block_size, global_j)
            local_B(local_i, local_j) = properties(global_i, global_j)
        enddo
    enddo


    ! Solver
    call pdgels("N", nm1, nm2, n_properties, local_K, 1, 1, desca, local_B, 1, 1, DESCB, work, lwork, info)
    deallocate(work)


    ! Copy LAPACK output
    alphas = 0.0d0
    do local_j = 1, local_B_cols
        do local_i = 1, local_B_rows
            call l2g(local_i, local_rank_row, ranks_rows, block_size, global_i)
            call l2g(local_j, local_rank_col, ranks_cols, block_size, global_j)
            alphas(global_i, global_j) = local_B(local_i, local_j)
        enddo
    enddo


    ! TODO Return alphas
    call DGSUM2D(context, "All", "1-tree", nm1, n_properties, alphas, 1, -1, -1)

    ! Save alphas to file
    if (local_id .eq. 0) then
        open(unit = 9, file = "jobname_dev/test_qm7_nn_2p_alpha", form="formatted")
        call print_matrix(alphas, n_properties, nm1, 9)
        close(9)
    end if


    ! Shutdown

    ! Deallocate reader
    deallocate(collection_x)
    deallocate(collection_y)
    deallocate(properties)
    deallocate(local_collection_x)
    deallocate(local_collection_y)
    deallocate(local_view_x)
    deallocate(local_view_y)

    ! Deallocate kernels
    deallocate(local_B)
    deallocate(local_K)
    deallocate(local_kernels)

    ! Tear down MPI
    call blacs_exit(0)

end program test_fchl_gatherv_kernel

subroutine l2g(il,p,np,nb,i)

   implicit none
   integer :: il   ! local array index, input
   integer :: p    ! processor array index, input
   integer :: np   ! processor array dimension, input
   integer :: nb   ! block size, input
   integer :: i    ! global array index, output
   integer :: ilm1

   ilm1 = il-1
   i    = (((ilm1/nb) * np) + p)*nb + mod(ilm1,nb) + 1

   return
end subroutine l2g


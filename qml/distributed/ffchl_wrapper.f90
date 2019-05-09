module ffchl_wrapper

    implicit none

contains

subroutine kernel_wrapper_fchl(collection_x, collection_y, kernels)

    use ffchl_reader, only: fread_fchl_args, collection2representation
    use ffchl_scalar, only: fget_kernels

    implicit none

    ! fchl collection of representations
    double precision, dimension(:,:), intent(in) :: collection_x
    double precision, dimension(:,:), intent(in) :: collection_y
    double precision, dimension(:,:,:), intent(out) :: kernels
    !

    ! fchl args
    double precision, allocatable, dimension(:,:,:,:) :: x1
    double precision, allocatable, dimension(:,:,:,:) :: x2

    integer :: nm1
    integer :: nm2
    integer :: max_size
    integer :: max_neighbors

    integer, dimension(:), allocatable :: n1
    integer, dimension(:), allocatable :: n2
    integer, dimension(:,:), allocatable :: nneigh1
    integer, dimension(:,:), allocatable :: nneigh2
    double precision, dimension(:), allocatable:: sigmas
    integer :: n_kernels
    double precision :: two_body_power
    double precision :: three_body_power
    double precision :: t_width
    double precision :: d_width
    double precision :: cut_start
    double precision :: cut_distance
    integer :: order
    double precision :: distance_scale
    double precision :: angular_scale
    double precision, dimension(:,:), allocatable :: pd
    logical :: alchemy

    ! Kernel ID and corresponding parameters
    integer :: kernel_idx
    double precision, dimension(:,:), allocatable :: parameters

    ! end fchl args

    ! collection sizes
    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end collection sizes

    integer, dimension(:), allocatable :: jk

    ! size of collections
    nm1 = size(collection_x, 1)
    nm2 = size(collection_y, 1)

    ! hardcode fchl sizes
    n_kernels = 1
    ! end hardcode fchl sizes

    ! read fchl args
    ! TODO Broadcast this stuff
    call fread_fchl_args("jobname_dev/_fchl", max_size, max_neighbors, &
       & n_kernels, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
       & kernel_idx, parameters)


    ! collection index
    collection_size = max_size*5*max_neighbors + 1 + max_size
    idx_x = 1
    idx_n = idx_x + max_size*5*max_neighbors
    idx_nneigh = idx_n + 1
    ! end collection index

    ! allocate fchl
    allocate(x1(nm1, max_size, 5, max_neighbors))
    allocate(x2(nm1, max_size, 5, max_neighbors))

    allocate(n1(nm1))
    allocate(n2(nm2))

    allocate(nneigh1(nm1, max_size))
    allocate(nneigh2(nm2, max_size))
    ! end allocate fchl


    ! Reshape collection to representations
    call collection2representation(collection_x, x1, n1, nneigh1, max_size, max_neighbors)
    call collection2representation(collection_y, x2, n2, nneigh2, max_size, max_neighbors)

    ! call the ffchl kernel function
    call fget_kernels(x1, x2, .False., n1, n2, nneigh1, nneigh2, nm1, nm2, n_kernels, &
         & t_width, d_width, cut_start, cut_distance, order, pd, &
         & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, &
         & kernel_idx, parameters, kernels)

    deallocate(x1)
    deallocate(x2)

    deallocate(n1)
    deallocate(n2)
    deallocate(nneigh1)
    deallocate(nneigh2)
    deallocate(pd)

end subroutine

end module ffchl_wrapper

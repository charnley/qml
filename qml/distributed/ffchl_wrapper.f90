module ffchl_wrapper

    implicit none

contains

subroutine collection2representation()

    implicit none

    ! allocate(representations(nm1, max_size, 5, max_size))
    ! call fread_fchl_representations("data/qm7_fchl_representations", nm1, max_size, 5, max_size, representations)
    !
    ! allocate(n1(nm1))
    ! allocate(n2(nm2))
    ! call fread_1d_integer("data/qm7_fchl_n1", nm1, n1)
    ! call fread_1d_integer("data/qm7_fchl_n1", nm2, n2)
    !
    ! allocate(nneigh1(nm1, max_size))
    ! allocate(nneigh2(nm2, max_size))
    ! call fread_2d_integer("data/qm7_fchl_neighbors1", nm1, max_size, nneigh1)
    ! call fread_2d_integer("data/qm7_fchl_neighbors2", nm2, max_size, nneigh2)
    !
    ! ! Translate representations to global 2d collection
    ! !
    ! !       x1                       n1        nneigh1
    ! ! max_size*5*max_neighbors  +    1    +    max_size
    ! !
    !
    ! collection_size = max_size*5*max_neighbors + 1 + max_size
    ! idx_x = 1
    ! idx_n = idx_x + max_size*5*max_neighbors
    ! idx_nneigh = idx_n + 1
    !
    ! allocate(collection_x( nm1, collection_size ))
    ! collection_x = 0.0d0
    ! collection_x(:,idx_x:idx_n-1) = reshape(source=representations, shape=[nm1, max_size*5*max_neighbors])
    ! collection_x(:,idx_n) = n1
    ! collection_x(:,idx_nneigh:collection_size) = nneigh1
    !
    ! allocate(collection_y( nm2, collection_size ))
    ! collection_y = 0.0d0
    ! collection_y(:,idx_x:idx_n-1) = reshape(source=representations, shape=[nm2, max_size*5*max_neighbors])
    ! collection_y(:,idx_n) = n2
    ! collection_y(:,idx_nneigh:collection_size) = nneigh2



end subroutine

subroutine representation2collection()

    implicit none

end subroutine

subroutine kernel_wrapper_fchl(collection_x, collection_y, kernel)

    use ffchl_reader, only: fread_fchl_args
    use ffchl_scalar_kernels, only: fget_kernels_fchl

    implicit none

    ! fchl collection of representations
    double precision, dimension(:,:), intent(in) :: collection_x
    double precision, dimension(:,:), intent(in) :: collection_y
    double precision, dimension(:,:,:), intent(out) :: kernel
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
    integer :: nsigmas
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
    ! end fchl args

    ! collection sizes
    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end collection sizes

    ! hardcode fchl sizes
    nm1 = 11 ! TODO read from dimension of collection
    nm2 = nm1
    max_size = 23 ! TODO read from broadcast / files
    max_neighbors = 23
    nsigmas = 1
    ! end hardcode fchl sizes

    ! collection index
    collection_size = max_size*5*max_neighbors + 1 + max_size
    idx_x = 1
    idx_n = idx_x + max_size*5*max_neighbors
    idx_nneigh = idx_n + 1
    ! end collection index

    ! allocate fchl
    allocate(x1(nm1, max_size, 5, max_neighbors))
    allocate(x2(nm1, max_size, 5, max_neighbors))

    allocate(sigmas(nsigmas))
    allocate(n1(nm1))
    allocate(n2(nm2))
    allocate(nneigh1(nm1, max_size))
    allocate(nneigh2(nm2, max_size))

    allocate(pd(100,100))
    ! end allocate fchl

    ! read fchl args
    ! TODO Broadcast this stuff
    call fread_fchl_args(max_size, max_neighbors, &
       & sigmas, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power)


    ! Reshape collection to representations
    x1 = reshape(source=collection_x(:,idx_x:idx_n-1), shape=[nm1, max_size, 5, max_neighbors])
    x2 = reshape(source=collection_y(:,idx_x:idx_n-1), shape=[nm2, max_size, 5, max_neighbors])

    nneigh1 = reshape(source=collection_x(:,idx_nneigh:collection_size), shape=[nm1, max_size])
    nneigh2 = reshape(source=collection_y(:,idx_nneigh:collection_size), shape=[nm2, max_size])

    n1 = collection_x(:,idx_n)
    n2 = collection_y(:,idx_n)

    ! call the kernel function
    ! call fjckget_kernels_fchl(x1, x2, n1, n2, nneigh1, nneigh2, &
    call fget_kernels_fchl(x1, x2, n1, n2, nneigh1, nneigh2, &
       & sigmas, nm1, nm2, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power, kernel)

    deallocate(x1)
    deallocate(x2)

    deallocate(sigmas)
    deallocate(n1)
    deallocate(n2)
    deallocate(nneigh1)
    deallocate(nneigh2)
    deallocate(pd)

end subroutine

end module ffchl_wrapper

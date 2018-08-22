module ffchl_reader

    implicit none

contains


subroutine collection2representation( &
    collection, &
    xx, nx, nneighx, &
    max_size, max_neighbors)

    !
    ! Translate representations to global 2d collection
    !
    !       x1                       n1        nneigh1
    ! max_size*5*max_neighbors  +    1    +    max_size
    !

    implicit none

    ! fchl collection of representations
    double precision, dimension(:,:), intent(in) :: collection
    !

    ! fchl args
    double precision, dimension(:,:,:,:), intent(out) :: xx
    integer, dimension(:), intent(out) :: nx
    integer, dimension(:,:), intent(out) :: nneighx

    integer :: nmx
    integer, intent(in) :: max_size
    integer, intent(in) :: max_neighbors
    ! end fchl args

    ! collection sizes
    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end collection sizes


    ! Size of first dimension
    nmx = size(collection, 1)

    ! collection idx
    collection_size = max_size*5*max_neighbors + 1 + max_size
    idx_x = 1
    idx_n = idx_x + max_size*5*max_neighbors
    idx_nneigh = idx_n + 1
    !

    ! Reshape
    xx = reshape(source=collection(:,idx_x:idx_n-1), shape=[nmx, max_size, 5, max_neighbors])
    nneighx = reshape(source=collection(:,idx_nneigh:collection_size), shape=[nmx, max_size])
    nx = collection(:,idx_n)

end subroutine

subroutine representation2collection( &
    collection, &
    xx, nx, nneighx, &
    max_size, max_neighbors)

    !
    ! Translate representations to global 2d collection
    !
    !       x1                       n1        nneigh1
    ! max_size*5*max_neighbors  +    1    +    max_size
    !

    implicit none

    ! fchl collection of representations
    double precision, dimension(:,:), intent(out) :: collection
    !

    ! fchl args
    double precision, dimension(:,:,:,:), intent(in) :: xx
    integer, dimension(:), intent(in) :: nx
    integer, dimension(:,:), intent(in) :: nneighx

    integer :: nmx
    integer, intent(in) :: max_size
    integer, intent(in) :: max_neighbors
    ! end fchl args

    ! collection sizes
    integer :: collection_size
    integer :: idx_x
    integer :: idx_n
    integer :: idx_nneigh
    ! end collection sizes

    ! Size of first dimension
    nmx = size(collection, 1)

    ! collection idx
    collection_size = max_size*5*max_neighbors + 1 + max_size
    idx_x = 1
    idx_n = idx_x + max_size*5*max_neighbors
    idx_nneigh = idx_n + 1
    !

    collection = 0.0d0
    collection(:,idx_x:idx_n-1) = reshape(source=xx, shape=[nmx, max_size*5*max_neighbors])
    collection(:,idx_n) = nx
    collection(:,idx_nneigh:collection_size) = nneighx

end subroutine


subroutine fread_fchl_collection(fname, collection, nmx, max_size, max_neighbors)

    use freaders, only: fread_1d_integer, fread_2d_integer

    implicit none

    character(len=*), intent(in) :: fname
    double precision, dimension(:,:), intent(out):: collection

    integer, intent(in) :: nmx
    integer, intent(in) :: max_size
    integer, intent(in) :: max_neighbors

    double precision, dimension(:,:,:,:), allocatable :: representations

    integer, dimension(:), allocatable :: nx
    integer, dimension(:,:), allocatable :: nneighx

    integer :: ni, nj, nk, nl


    ! Read all fchl representations
    allocate(representations(nmx, max_size, 5, max_neighbors))
    call fread_fchl_representations(trim(fname) // trim("_x"), representations)

    allocate(nx(nmx))
    call fread_1d_integer(trim(fname) // trim("_n"), nmx, nx)

    allocate(nneighx(nmx, max_size))
    call fread_2d_integer(trim(fname) // trim("_neighbors"), nmx, max_size, nneighx)

    ! Translate
    call representation2collection(collection, representations, nx, nneighx, max_size, max_neighbors)

end subroutine


subroutine fread_fchl_representations(fname, representations)

    implicit none

    character(len=*), intent(in) :: fname
    double precision, dimension(:,:,:,:), intent(out):: representations

    ! write(*,*) shape(representations)

    open(21, file=fname, form='unformatted')
    read(21) representations
    close(21)

end subroutine


subroutine fread_fchl_args(fname, max_size, max_neighbors, &
       & sigmas, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power)

    use freaders, only: fread_logical, fread_integer, fread_double, fread_1d_integer, fread_2d_integer, fread_1d_double, fread_2d_double

    implicit none

    !jck
    character(len=*), intent(in) :: fname
    integer, intent(out) :: max_size
    integer, intent(out) :: max_neighbors

    ! Sigma in the Gaussian kernel
    double precision, dimension(:), allocatable, intent(out) :: sigmas

    ! Number of sigmas
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

    double precision, allocatable, dimension(:,:), intent(out) :: pd

    logical :: alchemy

    max_size = 23
    max_neighbors = 23

    allocate(sigmas(nsigmas))

    allocate(pd(100,100))

    call fread_logical(trim(fname) // trim("_doalchemy"), alchemy)

    call fread_integer(trim(fname) // trim("_fourier_order"), order)

    call fread_double(trim(fname) // trim("_two_body_power"), two_body_power)
    call fread_double(trim(fname) // trim("_two_body_width"), d_width)
    call fread_double(trim(fname) // trim("_three_body_power"), three_body_power)
    call fread_double(trim(fname) // trim("_three_body_width"), t_width)

    call fread_double(trim(fname) // trim("_cut_start"), cut_start)
    call fread_double(trim(fname) // trim("_cut_distance"), cut_distance)

    call fread_double(trim(fname) // trim("_two_body_scaling"), distance_scale)
    call fread_double(trim(fname) // trim("_three_body_scaling"), angular_scale)

    call fread_1d_double(trim(fname) // trim("_sigmas"), nsigmas, sigmas)

    call fread_2d_double(trim(fname) // trim("_pd"), 100, 100, pd)


end subroutine


end module ffchl_reader

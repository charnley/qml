module ffchl_reader

    implicit none

contains

subroutine fread_fchl_collection(fname, collection, nmx, max_size, max_neighbors)

    use ffchl_wrapper, only: representation2collection
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


subroutine fread_fchl_args(max_size, max_neighbors, &
       & sigmas, nsigmas, &
       & t_width, d_width, cut_start, cut_distance, order, pd, &
       & distance_scale, angular_scale, alchemy, two_body_power, three_body_power)

    use freaders, only: fread_logical, fread_integer, fread_double, fread_1d_integer, fread_2d_integer, fread_1d_double, fread_2d_double

    implicit none

    !jck
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

    call fread_logical("data/qm7_fchl_doalchemy", alchemy)

    call fread_integer("data/qm7_fchl_fourier_order", order)

    call fread_double("data/qm7_fchl_two_body_power", two_body_power)
    call fread_double("data/qm7_fchl_two_body_width", d_width)
    call fread_double("data/qm7_fchl_three_body_power", three_body_power)
    call fread_double("data/qm7_fchl_three_body_width", t_width)

    call fread_double("data/qm7_fchl_cut_start", cut_start)
    call fread_double("data/qm7_fchl_cut_distance", cut_distance)

    call fread_double("data/qm7_fchl_two_body_scaling", distance_scale)
    call fread_double("data/qm7_fchl_three_body_scaling", angular_scale)

    call fread_1d_double("data/qm7_fchl_sigmas", nsigmas, sigmas)

    call fread_2d_double("data/qm7_fchl_pd", 100, 100, pd)


end subroutine


end module ffchl_reader

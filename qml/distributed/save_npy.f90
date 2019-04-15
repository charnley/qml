
subroutine fsave_fchl_representations(fname, representations)

    implicit none

    character(len=*), intent(in) :: fname
    double precision, dimension(:,:,:,:), intent(in) :: representations

    integer :: ni, nj, nk, nl
    integer :: i,j,k,l

    ni = size(representations, dim=1)
    nj = size(representations, dim=2)
    nk = size(representations, dim=3)
    nl = size(representations, dim=4)

    ! do i=1,ni
    !     do j=1,nj
    !         do k=1,nk
    !             do l=1,nl
    !                 write(*,*) representations(i,j,k,l)
    !             end do
    !         end do
    !     end do
    ! end do

    open(21, file=fname, form='unformatted')
    write(21) representations
    close(21)

end subroutine


subroutine fread_fchl_representations(fname, ni, nj, nk, nl)

    implicit none

    character(len=*), intent(in) :: fname
    double precision, allocatable, dimension(:,:,:,:) :: representations

    integer :: ni, nj, nk, nl
    integer :: i,j,k,l


    allocate(representations(ni, nj, nk, nl))

    open(21, file=fname, form='unformatted')
    read(21) representations
    close(21)

    ! ni = size(representations, dim=1)
    ! nj = size(representations, dim=2)
    ! nk = size(representations, dim=3)
    ! nl = size(representations, dim=4)

    write(*,*) "jck ni", ni


    do i=1,ni
        do j=1,nj
            do k=1,nk
                do l=1,nl
                    write(*,*) representations(i,j,k,l)
                end do
            end do
        end do
    end do

end subroutine

subroutine fsave_fchl_kernel(fname, kernel)

    implicit none

    character(len=*), intent(in) :: fname
    double precision, dimension(:,:), intent(in) :: kernel

    integer :: ni, nj
    integer :: i,j

    ni = size(kernel, dim=1)
    nj = size(kernel, dim=2)

    do i=1,ni
        do j=1,nj
            write(*,*) kernel(i,j)
        end do
    end do

    open(21, file=fname, form='unformatted')
    write(21) kernel
    close(21)

end subroutine


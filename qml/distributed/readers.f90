module freaders

    implicit none

contains

subroutine fread_logical(fname, arg)
    implicit none
    character(len=*), intent(in) :: fname
    integer :: i
    logical, intent(out) :: arg
    open(2, file=fname)
    read(2, *) i
    close(2)
    if (i .eq. 0) then
      arg = .T.
    else
      arg = .F.
    endif
end subroutine

subroutine fread_integer(fname, arg)
    implicit none
    character(len=*), intent(in) :: fname
    integer, intent(out) :: arg
    open(2, file=fname)
    read(2, *) arg
    close(2)
end subroutine

subroutine fread_double(fname, arg)
    implicit none
    character(len=*), intent(in) :: fname
    double precision, intent(out) :: arg
    open(2, file=fname, form='formatted')
    read(2, '(f)') arg
    close(2)
end subroutine

subroutine fread_1d_integer(fname, n, arg)
    implicit none
    integer :: n, i
    character(len=*), intent(in) :: fname
    integer, dimension(:), intent(out) :: arg
    open(2, file=fname)
    do i=1,n
        read(2, *) arg(i)
    end do
    close(2)
end subroutine

subroutine fread_2d_integer(fname, n, m, arg)
    implicit none
    integer :: n, m, i
    character(len=*), intent(in) :: fname
    integer, dimension(:,:), intent(out) :: arg
    open(2, file=fname)
    do i=1,n
        read(2, *) arg(i,:)
    end do
    close(2)
end subroutine

subroutine fread_1d_double(fname, n, arg)
    implicit none
    integer :: n, i
    character(len=*), intent(in) :: fname
    double precision, dimension(:), intent(out) :: arg
    open(2, file=fname)
    do i=1,n
        read(2, *) arg(i)
    end do
    close(2)
end subroutine

subroutine fread_2d_double(fname, n, m, arg)
    implicit none
    integer :: n, m, i
    character(len=*), intent(in) :: fname
    double precision, dimension(:,:), intent(out) :: arg
    open(2, file=fname)
    do i=1,n
        read(2, *) arg(i,:)
    end do
    close(2)
end subroutine

end module freaders

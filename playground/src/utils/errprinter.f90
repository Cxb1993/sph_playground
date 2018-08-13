module errprinter
  implicit none

  public :: error

  private

  interface error
    module procedure errstr, errint
  end interface

contains
  subroutine errstr(str, itt, file, line)
    integer, intent(in) :: line
    character(len=*), intent(in) :: file, str, itt
    character(len=30) :: linestr

    write(linestr, *) line

    print*, '# <!> ' // str // " | " // itt
    print*, '# <!> ' // file //":"// adjustl(linestr)
    stop
  end subroutine errstr

  subroutine errint(str, itt, file, line)
    integer, intent(in) :: line, itt
    character(len=*), intent(in) :: file, str
    character(len=30) :: linestr, ittstr

    write(linestr, *) line
    write(ittstr, *) itt

    print*, '# <!> ' // str // " | " // ittstr
    print*, '# <!> ' // file //":"// adjustl(linestr)
    stop
  end subroutine errint

end module errprinter

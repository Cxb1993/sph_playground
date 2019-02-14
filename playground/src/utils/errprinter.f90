module errprinter
  implicit none

  public :: error, warning

  private

  interface error
    module procedure errstr, errint
  end interface
  interface warning
    module procedure warstr, warint, warreal
  end interface

contains
  subroutine errstr(str, itt, file, line)
    integer, intent(in) :: line
    character(len=*), intent(in) :: file, str, itt
    character(len=30) :: linestr

    write(linestr, *) line

    print*, '# <!> ' // str // " | " // itt // file //":"// adjustl(linestr)
    stop
  end subroutine errstr

  subroutine errint(str, itt, file, line)
    integer, intent(in) :: line, itt
    character(len=*), intent(in) :: file, str
    character(len=30) :: linestr, ittstr

    write(linestr, *) line
    write(ittstr, *) itt

    print*, '# <!> ' // str // " | " // trim(adjustl(trim(ittstr))) // " | " // file //":"// adjustl(linestr)
    stop
  end subroutine errint

  subroutine warstr(str, itt, file, line)
    integer, intent(in) :: line
    character(len=*), intent(in) :: file, str, itt
    character(len=30) :: linestr

    write(linestr, *) line

    if (itt == "") then
      print*, '# <?> '//str//" | "//file//":"//adjustl(linestr)
    else
      print*, '# <?> '//str//" | "//adjustl(trim(itt))//" | "//file//":"//adjustl(linestr)
    end if
  end subroutine warstr

  subroutine warint(str, itt, file, line)
    integer, intent(in) :: line, itt
    character(len=*), intent(in) :: file, str
    character(len=30) :: linestr, ittstr

    write(linestr, *) line
    write(ittstr, *) itt

    print*, '# <?> ' // str // " | " // ittstr // file //":"// adjustl(linestr)
  end subroutine warint

  subroutine warreal(str, itt, file, line)
    integer, intent(in) :: line
    real, intent(in) :: itt
    character(len=*), intent(in) :: file, str
    character(len=30) :: linestr, ittstr

    write(linestr, *) line
    write(ittstr, *) itt

    print*, '# <?> ' // str // " | " // ittstr // file //":"// adjustl(linestr)
  end subroutine warreal

end module errprinter

module n2ext
  use const
  implicit none

  public :: n2f, n2df, n2ddf, n2C, n2R, n2Name

  private

    real :: n2C(3,2) = reshape((/ -0.261008192263701 + 1.0/sqrt(pi), 1.0/pi, 1.0*pi**(-1.5),&
                                  -0.261008192263701 + 1.0/sqrt(pi), 1.0/pi, 1.0*pi**(-1.5) /),(/3,2/))
    character (len=10) :: n2Name=' genesis '
    real :: n2R = 2.0
 contains

  subroutine n2f(r, h, f)
    real, intent(in)  :: r, h
    real, intent(out) :: f
    real              :: q

    q = r / h
    if (q < 0.0) then
      print *, 'f: something went wrong, q =', q
      stop
    elseif (q < 2.0) then
      f  = 0.0183127964683465*q**3 - 0.174477294027598*q**2 - 1.7053968728626432*q - (-2.2035977856942233*q + 1.8169193558843038)*erf(1.5690821350931794*q - 1.2937459461663592) + 0.14859319297412515*exp(q*(-2.4620187466685706*q + 4.059987302957713)) + 1.369093749870065
    else
      f = .0
    end if
  end subroutine n2f

  subroutine n2df(r, h, df)
    real, intent(in)  :: r, h
    real, intent(out) :: df
    real              :: q

    q = r / h
    if (q < 0.0) then
      print *, 'df: something went wrong, q =', q
      stop
    elseif (q < 2.0) then
      df  = 0.0549383894050395*q**2 - 0.348954588055196*q + 2.2035977856942233*erf(1.5690821350931794*q - 1.2937459461663592) - 1.7053968728626432
    else
      df = .0
    end if
  end subroutine n2df

  subroutine n2ddf(r, h, ddf)
    real, intent(in)  :: r, h
    real, intent(out) :: ddf
    real              :: q

    q = r / h
    if (q < 0.0) then
      print *, 'ddf: something went wrong, q =', q
      stop
    elseif (q < 2.0) then
      ddf  = 0.109876778810079*q - 0.348954588055196 + 3.90151305400392*exp(-2.46201874666857*(q - 0.824524043216852)**2)
    else
      ddf = .0
    end if
  end subroutine n2ddf
end module n2ext

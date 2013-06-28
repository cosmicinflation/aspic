module infprec
  implicit none

  public

!quad precision


!double precision
#ifdef QUADPREC
  integer, parameter :: kp = kind(1.0_16)
#else
  integer, parameter :: kp = kind(1.0_8)
#endif

!home made precision: p number of digit
! integer, parameter :: kp = selected_real_kind(p=32)

!default integration accuracy
  real(kp), parameter :: tolkp = 10000._kp * epsilon(1._kp)

!increased integration accuracy
!  real(kp), parameter :: tolkp = 1.d-36

  real(kp), parameter :: pi = 3.1415926535897932384626433832795_kp
  real(kp), parameter :: euler = 0.5772156649015328606065120900824024310422_kp
  real(kp), parameter :: CConst = euler + log(2._kp) - 2._kp

!workaround for passing argument to old f77 functions. Only pointer
!can be deferred shape in derived data type.
!Allows to stop integration from conditions coming from called
!functions (find the end of inflation)
  type transfert
     logical :: yesno1,yesno2, yesno3, yesno4
     integer :: int1, int2, int3
     real(kp) :: real1, real2, real3, real4, real5, real6
     real(kp), dimension(:), pointer :: ptrvector1 => null()
     real(kp), dimension(:), pointer :: ptrvector2 => null()
!reserved
     logical :: check,update
     real(kp) :: xend
  end type transfert


end module infprec

module infprec
  implicit none

  public

!quad precision
!  integer, parameter :: kp = kind(1.0_16)

!double precision
  integer, parameter :: kp = kind(1.0_8)

!home made precision: p number of digit
! integer, parameter :: kp = selected_real_kind(p=32)

!default integration accuracy
  real(kp), parameter :: tolkp = 1.d-12

!workaround for passing argument to old f77 functions. Only pointer
!can be deferred shape in derived data type.
!Allows to stop integration from conditions coming from called
!functions (find the end of inflation)
  type transfert
     logical :: yesno1,yesno2, yesno3, yesno4
     integer :: int1, int2, int3
     real(kp) :: real1, real2, real3, real4, real5
     real(kp), dimension(:), pointer :: ptrvector1 => null()
     real(kp), dimension(:), pointer :: ptrvector2 => null()
!reserved
     logical :: check,update
     real(kp) :: xend
  end type transfert


end module infprec

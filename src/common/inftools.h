  interface
     subroutine fcn(n, x, y, yprime, otherdata)
       use infprec, only : kp,transfert         
       implicit none          
       integer :: n
       real(kp) :: x
       real(kp), dimension(n) :: y, yprime
       type(transfert), optional, intent(inout) :: otherdata
     end subroutine fcn
  end interface


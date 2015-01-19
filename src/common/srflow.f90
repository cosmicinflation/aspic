module srflow
!this is in k* defined by k* eta* = -1 and eps() stands for the Hubble
!flow parameters evaluated at eta*
  use infprec, only : kp, pi, CConst

  implicit none

  private

  real(kp), parameter :: pi2 = pi*pi

  logical, parameter :: display = .false.

  public slowroll_violated
  public tensor_to_scalar_ratio
  public scalar_spectral_index, tensor_spectral_index
  public scalar_running, tensor_running
  public scalar_running_running, tensor_running_running
  public slowroll_corrections, ln_slowroll_corrections

! because 1/slowroll_corrections >< inverse_slowroll_corrections
! it is better to numerically stick to one choice only
!  public inverse_slowroll_corrections, ln_inverse_slowroll_corrections
contains


  function slowroll_violated(eps)
    implicit none
    logical :: slowroll_violated
    real(kp), dimension(:), intent(in) :: eps

    slowroll_violated =  any(abs(eps).gt.1._kp)       

  end function slowroll_violated



  function scalar_spectral_index(eps)
    implicit none
    real(kp) :: scalar_spectral_index
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: nsm1
    integer :: neps

    neps = size(eps,1)

    select case (neps)
              
    case (2)
       nsm1 = -2._kp*eps(1) - eps(2)

    case (3)
       nsm1 = -2._kp*eps(1) - eps(2)
       nsm1 = nsm1 - 2._kp*eps(1)**2 - (3._kp+2._kp*CConst)*eps(1)*eps(2) &
            - CConst*eps(2)*eps(3)
      
    case default
       stop 'scalar_spectral_index: neps not implemented!'
    end select


    scalar_spectral_index = nsm1 + 1._kp

  end function scalar_spectral_index



  function tensor_spectral_index(eps)
    implicit none
    real(kp) :: tensor_spectral_index
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: nt
    integer :: neps

    neps = size(eps,1)

    select case (neps)
              
    case (1)
       nt = -2._kp*eps(1)

    case (2)
       nt = -2._kp*eps(1) &
            - 2._kp*eps(1)**2 - 2._kp*(1._kp+CConst)*eps(1)*eps(2)

    case(3)
       nt = -2._kp*eps(1) &
            - 2._kp*eps(1)**2 - 2._kp*(1._kp+CConst)*eps(1)*eps(2) &
            - 2._kp*eps(1)**3 - (14._kp+6._kp*CConst-pi2)*eps(1)**2*eps(2) &
            - (2._kp + 2._kp*CConst*(1._kp+CConst) -pi2/12._kp)*eps(1)*eps(2)*(eps(2)+eps(3))
            
    case default
       stop 'tensor_spectral_index: neps not implemented!'
    end select


    tensor_spectral_index = nt

  end function tensor_spectral_index



  function tensor_to_scalar_ratio(eps)
    implicit none
    real(kp) :: tensor_to_scalar_ratio
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: r
    integer :: neps

    neps = size(eps,1)

    select case (neps)
       
    case (1)
       r = 1._kp

    case (2)
       r = 1._kp + CConst*eps(2)

    case (3)
       r = 1._kp + CConst*eps(2) &
            + (CConst - pi2/2._kp + 4._kp)*eps(1)*eps(2) &
            + (0.5_kp*CConst**2 - pi2/8._kp + 1._kp)*eps(2)**2 &
            + (0.5_kp*CConst**2 - pi2/24._kp)*eps(2)*eps(3)

    case default
       stop 'tensor_to_scalar_ratio: neps not implemented!'
    end select

!that may happen if eps2 > 1, i.e. when slow-roll is violated. If
!this is the case, we use the zero order formula. If this is not the
!case, this is nasty and we abort.
    if (r.lt.0._kp) then
       if (slowroll_violated(eps)) then
          if (display) write(*,*),'tensor_to_scalar_ratio: eps(:)= ',eps(:)
          tensor_to_scalar_ratio = 16._kp*eps(1)
       else
          stop 'tensor_to_scalar_ratio: r < 0!'
       endif
    else
       tensor_to_scalar_ratio = 16._kp*eps(1)*r
    end if


    
  end function tensor_to_scalar_ratio


  function scalar_running(eps)
    implicit none
    real(kp) :: scalar_running
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: alpha
    integer :: neps

    neps = size(eps,1)

    select case (neps)
       
    case (1,2)
       alpha = 0._kp

    case (3)
       alpha = -2._kp*eps(1)*eps(2) - eps(2)*eps(3)
      
    case default
       stop 'scalar_running: neps not implemented!'
    end select
    
    scalar_running = alpha

  end function scalar_running




  function tensor_running(eps)
    implicit none
    real(kp) :: tensor_running
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: alpha
    integer :: neps

    neps = size(eps,1)

    select case (neps)
       
    case (1)
       alpha = 0._kp

    case (2)
       alpha = -2._kp*eps(1)*eps(2)

    case (3)
       alpha = -2._kp*eps(1)*eps(2) &
            - 6._kp*eps(1)**2*eps(2) - 2._kp*(1._kp+CConst)*eps(1)*eps(2)**2 &
            - 2._kp*(1._kp+CConst)*eps(1)*eps(2)*eps(3)
      
    case default
       stop 'scalar_running: neps not implemented!'
    end select
    
    tensor_running = alpha

  end function tensor_running




  function scalar_running_running(eps)
    implicit none
    real(kp) :: scalar_running_running
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: beta
    integer :: neps

    neps = size(eps,1)

    select case (neps)
       
    case (1,2,3)
       beta = 0._kp

    case (4)
       beta = -2._kp*eps(1)*eps(2)**2 - 2._kp*eps(1)*eps(2)*eps(3) &
            - eps(2)*eps(3)**2 - eps(2)*eps(3)*eps(4)
      
    case default
       stop 'tensor_running_running: neps not implemented!'
    end select
    
    scalar_running_running = beta


  end function scalar_running_running




  function tensor_running_running(eps)
    implicit none
    real(kp) :: tensor_running_running
    real(kp), intent(in), dimension(:) :: eps
    real(kp) :: beta
    integer :: neps

    neps = size(eps,1)

    select case (neps)
       
    case (1,2)
       beta = 0._kp

    case (3)
       beta = -2._kp*eps(1)*eps(2)*(eps(2)+eps(3))
      
    case default
       stop 'tensor_running_running: neps not implemented!'
    end select
    
    tensor_running_running = beta


  end function tensor_running_running




!this gives the normal slow-roll expansion of P(k*)/[H*^2/(8pi^2 eps1*)]
  function slowroll_corrections(eps)
    implicit none
    real(kp) :: slowroll_corrections
    real(kp), dimension(:), intent(in) :: eps
    integer :: neps

    neps = size(eps,1)

    if (slowroll_violated(eps)) then
       write(*,*)'eps= ',eps
       stop 'slowroll_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       slowroll_corrections = 1._kp
       
    case (2)
       slowroll_corrections = 1._kp - 2._kp*(1._kp+Cconst)*eps(1) &
            - CConst*eps(2)

    case (3)
       slowroll_corrections = 1._kp - 2._kp*(1._kp+Cconst)*eps(1) &
            - CConst*eps(2) &
            + (-3._kp + 2._kp*CConst + 2._kp*CConst**2 + pi2/2._kp) &
            *eps(1)**2 &
            + (-6._kp - Cconst + CConst**2 + 7._kp*pi2/12._kp) &
            *eps(1)*eps(2) &
            + (-1._kp + CConst**2/2._kp + pi2/8._kp)*eps(2)**2 &
            + (-CConst**2/2._kp + pi2/24._kp)*eps(2)*eps(3)

    case default
       stop 'inverse_slow_corrections: order not implemented!'

    end select

  end function slowroll_corrections



!this gives [H*^2/(8pi^2 eps1*)]/P(k*) consitently expanded in slow-roll
!hence the name "inverse" as opposed to normal slow-roll expansion
  function inverse_slowroll_corrections(eps)
    implicit none
    real(kp) :: inverse_slowroll_corrections
    real(kp), dimension(:), intent(in) :: eps
    integer :: neps

    neps = size(eps,1)

    if (slowroll_violated(eps)) then
       write(*,*)'eps= ',eps
       stop 'inverse_slowroll_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       inverse_slowroll_corrections = 1._kp
       
    case (2)
       inverse_slowroll_corrections = 1._kp + 2._kp*(1._kp+Cconst)*eps(1) &
            + CConst*eps(2)

    case (3)
       inverse_slowroll_corrections = 1._kp + 2._kp*(1._kp+Cconst)*eps(1) &
            + CConst*eps(2) &
            + (7._kp + 6._kp*CConst + 2._kp*CConst**2 - pi2/2._kp) &
            *eps(1)**2 &
            + (6._kp + 5._kp*Cconst + 3._kp*CConst**2 - 7._kp*pi2/12._kp) &
            *eps(1)*eps(2) &
            + (1._kp + CConst**2/2._kp - pi2/8._kp)*eps(2)**2 &
            + (CConst**2/2._kp - pi2/24._kp)*eps(2)*eps(3)

    case default
       stop 'inverse_slow_corrections: order not implemented!'

    end select

  end function inverse_slowroll_corrections



!this gives ln[P(k*)] - ln[H*^2/(8pi^2 eps1*)] consistently expanded
  function ln_slowroll_corrections(eps)
    implicit none
    real(kp) :: ln_slowroll_corrections
    real(kp), dimension(:), intent(in) :: eps
    integer :: neps
    
    neps = size(eps,1)

    if (slowroll_violated(eps)) then
       write(*,*)'eps= ',eps
       stop 'inverse_slowroll_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       ln_slowroll_corrections = 0._kp
       
    case (2)
       ln_slowroll_corrections = -2._kp*(1._kp+Cconst)*eps(1) &
            - CConst*eps(2)

    case (3)
       ln_slowroll_corrections = -2._kp*(1._kp+Cconst)*eps(1) &
            - CConst*eps(2) &
            - (5._kp + 2._kp*CConst - pi2/2._kp)*eps(1)**2 &
            - (6._kp + 3._kp*Cconst + CConst**2 - 7._kp*pi2/12._kp) &
            *eps(1)*eps(2) &
            - (1._kp - pi2/8._kp)*eps(2)**2 &
            - (CConst**2/2._kp - pi2/24._kp)*eps(2)*eps(3)

    case default
       stop 'ln_slow_corrections: order not implemented!'

    end select


  end function ln_slowroll_corrections


!this gives ln[H*^2/(8pi^2 eps1*)] - ln[P(k*)] consistently expanded
  function ln_inverse_slowroll_corrections(eps)
    implicit none
    real(kp) :: ln_inverse_slowroll_corrections
    real(kp), dimension(:), intent(in) :: eps
    integer :: neps

    ln_inverse_slowroll_corrections = -ln_slowroll_corrections(eps)

  end function ln_inverse_slowroll_corrections




end module srflow

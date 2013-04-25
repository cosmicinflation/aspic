module srflow
!this is in k* defined by k* eta* = -1
  use infprec, only : kp, pi, CConst

  implicit none

  private

  public slowroll_violated
  public scalar_spectral_index, tensor_to_scalar_ratio, scalar_running
  public inverse_slowroll_corrections

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
            + (CConst - pi**2/2._kp + 4._kp)*eps(1)*eps(2) &
            + (0.5_kp*CConst**2 - pi**2/8._kp + 1._kp)*eps(2)**2 &
            + (0.5_kp*CConst**2 - pi**2/24._kp)*eps(2)*eps(3)

    case default
       stop 'tensor_to_scalar_ratio: neps not implemented!'
    end select

    tensor_to_scalar_ratio = 16._kp*eps(1)*r

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


!this gives [H*^2/(8pi^2 eps1*)]/P(k*)
!hence the name "inverse" as opposed to normal slow-roll expansion
!giving P(k*)/[H*^2/(8pi^2 eps1*)]
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

    case default
       stop 'inverse_slow_corrections: order not implemented!'

    end select

  end function inverse_slowroll_corrections




end module srflow

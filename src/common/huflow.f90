module huflow
!this is in k* defined by k* eta* = -1 and epsH() stands for the Hubble
!flow parameters evaluated at eta*
  use infprec, only : kp, pi, CConst

  implicit none

  private

  real(kp), parameter :: pi2 = pi*pi

  logical, parameter :: display = .false.

  public hf_slowroll_violated
  public hf_tensor_to_scalar_ratio
  public hf_scalar_spectral_index, hf_tensor_spectral_index
  public hf_scalar_running, hf_tensor_running
  public hf_scalar_running_running, hf_tensor_running_running
  public hubbleflow_corrections, ln_hubbleflow_corrections
  public inverse_hubbleflow_corrections, ln_inverse_hubbleflow_corrections

contains


  function hf_slowroll_violated(epsH)
    implicit none
    logical :: hf_slowroll_violated
    real(kp), dimension(:), intent(in) :: epsH

    hf_slowroll_violated =  any(abs(epsH).gt.1._kp)       

  end function hf_slowroll_violated



  function hf_scalar_spectral_index(epsH)
    implicit none
    real(kp) :: hf_scalar_spectral_index
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: nsm1
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
              
    case (2)
       nsm1 = -2._kp*epsH(1) - epsH(2)

    case (3)
       nsm1 = -2._kp*epsH(1) - epsH(2)
       nsm1 = nsm1 - 2._kp*epsH(1)**2 - (3._kp+2._kp*CConst)*epsH(1)*epsH(2) &
            - CConst*epsH(2)*epsH(3)
      
    case default
       stop 'hf_scalar_spectral_index: neps not implemented!'
    end select


    hf_scalar_spectral_index = nsm1 + 1._kp

  end function hf_scalar_spectral_index



  function hf_tensor_spectral_index(epsH)
    implicit none
    real(kp) :: hf_tensor_spectral_index
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: nt
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
              
    case (1)
       nt = -2._kp*epsH(1)

    case (2)
       nt = -2._kp*epsH(1) &
            - 2._kp*epsH(1)**2 - 2._kp*(1._kp+CConst)*epsH(1)*epsH(2)

    case(3)
       nt = -2._kp*epsH(1) &
            - 2._kp*epsH(1)**2 - 2._kp*(1._kp+CConst)*epsH(1)*epsH(2) &
            - 2._kp*epsH(1)**3 - (14._kp+6._kp*CConst-pi2)*epsH(1)**2*epsH(2) &
            - (2._kp + 2._kp*CConst*(1._kp+CConst) -pi2/12._kp)*epsH(1)*epsH(2)*(epsH(2)+epsH(3))
            
    case default
       stop 'hf_tensor_spectral_index: neps not implemented!'
    end select


    hf_tensor_spectral_index = nt

  end function hf_tensor_spectral_index



  function hf_tensor_to_scalar_ratio(epsH)
    implicit none
    real(kp) :: hf_tensor_to_scalar_ratio
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: r
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
       
    case (1)
       r = 1._kp

    case (2)
       r = 1._kp + CConst*epsH(2)

    case (3)
       r = 1._kp + CConst*epsH(2) &
            + (CConst - pi2/2._kp + 4._kp)*epsH(1)*epsH(2) &
            + (0.5_kp*CConst**2 - pi2/8._kp + 1._kp)*epsH(2)**2 &
            + (0.5_kp*CConst**2 - pi2/24._kp)*epsH(2)*epsH(3)

    case default
       stop 'hf_tensor_to_scalar_ratio: neps not implemented!'
    end select

!that may happen if epsH2 > 1, i.e. when slow-roll is violated. If
!this is the case, we use the zero order formula. If this is not the
!case, this is nasty and we abort.
    if (r.lt.0._kp) then
       if (hf_slowroll_violated(epsH)) then
          if (display) write(*,*)'hf_tensor_to_scalar_ratio: epsH(:)= ',epsH(:)
          hf_tensor_to_scalar_ratio = 16._kp*epsH(1)
       else
          stop 'hf_tensor_to_scalar_ratio: r < 0!'
       endif
    else
       hf_tensor_to_scalar_ratio = 16._kp*epsH(1)*r
    end if


    
  end function hf_tensor_to_scalar_ratio


  function hf_scalar_running(epsH)
    implicit none
    real(kp) :: hf_scalar_running
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: alpha
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
       
    case (1,2)
       alpha = 0._kp

    case (3)
       alpha = -2._kp*epsH(1)*epsH(2) - epsH(2)*epsH(3)
      
    case default
       stop 'hf_scalar_running: neps not implemented!'
    end select
    
    hf_scalar_running = alpha

  end function hf_scalar_running




  function hf_tensor_running(epsH)
    implicit none
    real(kp) :: hf_tensor_running
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: alpha
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
       
    case (1)
       alpha = 0._kp

    case (2)
       alpha = -2._kp*epsH(1)*epsH(2)

    case (3)
       alpha = -2._kp*epsH(1)*epsH(2) &
            - 6._kp*epsH(1)**2*epsH(2) - 2._kp*(1._kp+CConst)*epsH(1)*epsH(2)**2 &
            - 2._kp*(1._kp+CConst)*epsH(1)*epsH(2)*epsH(3)
      
    case default
       stop 'hf_scalar_running: neps not implemented!'
    end select
    
    hf_tensor_running = alpha

  end function hf_tensor_running




  function hf_scalar_running_running(epsH)
    implicit none
    real(kp) :: hf_scalar_running_running
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: beta
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
       
    case (1,2,3)
       beta = 0._kp

    case (4)
       beta = -2._kp*epsH(1)*epsH(2)**2 - 2._kp*epsH(1)*epsH(2)*epsH(3) &
            - epsH(2)*epsH(3)**2 - epsH(2)*epsH(3)*epsH(4)
      
    case default
       stop 'hf_tensor_running_running: neps not implemented!'
    end select
    
    hf_scalar_running_running = beta


  end function hf_scalar_running_running




  function hf_tensor_running_running(epsH)
    implicit none
    real(kp) :: hf_tensor_running_running
    real(kp), intent(in), dimension(:) :: epsH
    real(kp) :: beta
    integer :: neps

    neps = size(epsH,1)

    select case (neps)
       
    case (1,2)
       beta = 0._kp

    case (3)
       beta = -2._kp*epsH(1)*epsH(2)*(epsH(2)+epsH(3))
      
    case default
       stop 'hf_tensor_running_running: neps not implemented!'
    end select
    
    hf_tensor_running_running = beta


  end function hf_tensor_running_running




!this gives the normal slow-roll expansion of P(k*)/[H*^2/(8pi^2 eps1*)]
  function hubbleflow_corrections(epsH)
    implicit none
    real(kp) :: hubbleflow_corrections
    real(kp), dimension(:), intent(in) :: epsH
    integer :: neps

    neps = size(epsH,1)

    if (hf_slowroll_violated(epsH)) then
       write(*,*)'epsH= ',epsH
       stop 'hubbleflow_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       hubbleflow_corrections = 1._kp
       
    case (2)
       hubbleflow_corrections = 1._kp - 2._kp*(1._kp+Cconst)*epsH(1) &
            - CConst*epsH(2)

    case (3)
       hubbleflow_corrections = 1._kp - 2._kp*(1._kp+Cconst)*epsH(1) &
            - CConst*epsH(2) &
            + (-3._kp + 2._kp*CConst + 2._kp*CConst**2 + pi2/2._kp) &
            *epsH(1)**2 &
            + (-6._kp - Cconst + CConst**2 + 7._kp*pi2/12._kp) &
            *epsH(1)*epsH(2) &
            + (-1._kp + CConst**2/2._kp + pi2/8._kp)*epsH(2)**2 &
            + (-CConst**2/2._kp + pi2/24._kp)*epsH(2)*epsH(3)

    case default
       stop 'inverse_slow_corrections: order not implemented!'

    end select

  end function hubbleflow_corrections



!this gives [H*^2/(8pi^2 eps1*)]/P(k*) consitently expanded in slow-roll
!hence the name "inverse" as opposed to normal slow-roll expansion
  function inverse_hubbleflow_corrections(epsH)
    implicit none
    real(kp) :: inverse_hubbleflow_corrections
    real(kp), dimension(:), intent(in) :: epsH
    integer :: neps

    neps = size(epsH,1)

    if (hf_slowroll_violated(epsH)) then
       write(*,*)'epsH= ',epsH
       stop 'inverse_hubbleflow_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       inverse_hubbleflow_corrections = 1._kp
       
    case (2)
       inverse_hubbleflow_corrections = 1._kp + 2._kp*(1._kp+Cconst)*epsH(1) &
            + CConst*epsH(2)

    case (3)
       inverse_hubbleflow_corrections = 1._kp + 2._kp*(1._kp+Cconst)*epsH(1) &
            + CConst*epsH(2) &
            + (7._kp + 6._kp*CConst + 2._kp*CConst**2 - pi2/2._kp) &
            *epsH(1)**2 &
            + (6._kp + 5._kp*Cconst + 3._kp*CConst**2 - 7._kp*pi2/12._kp) &
            *epsH(1)*epsH(2) &
            + (1._kp + CConst**2/2._kp - pi2/8._kp)*epsH(2)**2 &
            + (CConst**2/2._kp - pi2/24._kp)*epsH(2)*epsH(3)

    case default
       stop 'inverse_slow_corrections: order not implemented!'

    end select

  end function inverse_hubbleflow_corrections



!this gives ln[P(k*)] - ln[H*^2/(8pi^2 eps1*)] consistently expanded
  function ln_hubbleflow_corrections(epsH)
    implicit none
    real(kp) :: ln_hubbleflow_corrections
    real(kp), dimension(:), intent(in) :: epsH
    integer :: neps
    
    neps = size(epsH,1)

    if (hf_slowroll_violated(epsH)) then
       write(*,*)'epsH= ',epsH
       stop 'inverse_hubbleflow_corrections not reliable!'
    endif
    
    select case (neps)

    case (1)
       ln_hubbleflow_corrections = 0._kp
       
    case (2)
       ln_hubbleflow_corrections = -2._kp*(1._kp+Cconst)*epsH(1) &
            - CConst*epsH(2)

    case (3)
       ln_hubbleflow_corrections = -2._kp*(1._kp+Cconst)*epsH(1) &
            - CConst*epsH(2) &
            - (5._kp + 2._kp*CConst - pi2/2._kp)*epsH(1)**2 &
            - (6._kp + 3._kp*Cconst + CConst**2 - 7._kp*pi2/12._kp) &
            *epsH(1)*epsH(2) &
            - (1._kp - pi2/8._kp)*epsH(2)**2 &
            - (CConst**2/2._kp - pi2/24._kp)*epsH(2)*epsH(3)

    case default
       stop 'ln_slow_corrections: order not implemented!'

    end select


  end function ln_hubbleflow_corrections


!this gives ln[H*^2/(8pi^2 eps1*)] - ln[P(k*)] consistently expanded
  function ln_inverse_hubbleflow_corrections(epsH)
    implicit none
    real(kp) :: ln_inverse_hubbleflow_corrections
    real(kp), dimension(:), intent(in) :: epsH
    integer :: neps

    ln_inverse_hubbleflow_corrections = -ln_hubbleflow_corrections(epsH)

  end function ln_inverse_hubbleflow_corrections




end module huflow

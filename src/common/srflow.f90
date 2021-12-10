module srflow
!epsV() stands for the slow-roll parameters derived from the potential
!at *leading order*, i.e. the ones given by the aspic functions
  use infprec, only : kp

  implicit none

  private


  logical, parameter :: display = .false.

  public slowroll_to_hubble
  public slowroll_violated
  public tensor_to_scalar_ratio
  public scalar_spectral_index, tensor_spectral_index
  public scalar_running, tensor_running
  public scalar_running_running, tensor_running_running
  public slowroll_corrections, ln_slowroll_corrections

! because 1/slowroll_corrections >< inverse_slowroll_corrections
! it is better to numerically stick to one choice only
! public inverse_slowroll_corrections, ln_inverse_slowroll_corrections

contains



  function slowroll_to_hubble(epsV)
    implicit none
    real(kp), dimension(:), intent(in) :: epsV
    real(kp), dimension(size(epsV,1)) :: slowroll_to_hubble
    real(kp), dimension(size(epsV,1)) :: epsH

!epsH(n) depends on epsV(1:n+1); but epsV(n+1) appears in the power
!spectra expanded at order n, at order n+1, so it can be neglected
!(for the power spectra!!!!!!)
    real(kp), parameter :: epsVnp1 = 0._kp

    integer :: neps

    neps = size(epsV,1)

!do not add corrections which are valid only if epsV are small    
    if (slowroll_violated(epsV)) then
       if (display) write(*,*)'slowroll_to_hubble: no correction, epsV violates slow-roll!'
       slowroll_to_hubble = epsV
       return
    endif
    
    select case (neps)
    
       case (1,2)
          epsH = epsV

       case (3)
          epsH(1) = epsV(1) * (1._kp - epsV(2)/3._kp)
          epsH(2) = epsV(2) * (1._kp - epsV(2)/6._kp - epsV(3)/3._kp)
          epsH(3) = epsV(3) * (1._kp - epsV(2)/3._kp - epsVnp1/3._kp)

       case (4)
          epsH(1) = epsV(1) * (1._kp - epsV(2)/3._kp) &
               - epsV(1)**2*epsV(2)/9._kp + (5._kp/36._kp)*epsV(1)*epsV(2)**2 &
               + epsV(1)*epsV(2)*epsV(3)/9._kp
          epsH(2) = epsV(2) * (1._kp - epsV(2)/6._kp - epsV(3)/3._kp) &
               - epsV(1)*epsV(2)**2/6._kp + epsV(2)**3/18._kp - epsV(1)*epsV(2)*epsV(3)/9._kp &
               + (5._kp/18._kp)*epsV(2)**2*epsV(3) + epsV(2)*epsV(3)**2/9._kp &
               + epsV(2)*epsV(3)*epsV(4)/9._kp
          epsH(3) = epsV(3) * (1._kp - epsV(2)/3._kp - epsV(4)/3._kp) &
               - epsV(1)*epsV(2)**2/6 - epsV(1)*epsV(2)*epsV(3)/3._kp + epsV(2)**2*epsV(3)/6._kp &
               + (5._kp/18._kp)*epsV(2)*epsV(3)**2 - epsV(1)*epsV(3)*epsV(4)/9._kp &
               + (5._kp/18._kp)*epsV(2)*epsV(3)*epsV(4) + epsV(3)**2*epsV(4)/9._kp &
               + epsV(3)*epsV(4)**2/9._kp + epsV(3)*epsV(4)*epsVnp1/9._kp
       case default
          stop 'slowroll_to_hubble: order not implemented!'
       end select

       slowroll_to_hubble =  epsH


  end function slowroll_to_hubble




  function slowroll_violated(epsV)
    implicit none
    logical :: slowroll_violated
    real(kp), dimension(:), intent(in) :: epsV

    slowroll_violated = any(abs(epsV).gt.1._kp)

  end function slowroll_violated


  function slowroll_deeply_violated(epsV)
    implicit none
    logical :: slowroll_deeply_violated
    real(kp), dimension(:), intent(in) :: epsV

    slowroll_deeply_violated = any(abs(epsV).gt.3._kp)

  end function slowroll_deeply_violated



  function scalar_spectral_index(epsV)
    use huflow, only : hf_scalar_spectral_index
    implicit none
    real(kp) :: scalar_spectral_index
    real(kp), intent(in), dimension(:) :: epsV
    
    scalar_spectral_index =  hf_scalar_spectral_index(slowroll_to_hubble(epsV))

  end function scalar_spectral_index




  function tensor_spectral_index(epsV)
    use huflow, only : hf_tensor_spectral_index
    implicit none
    real(kp) :: tensor_spectral_index
    real(kp), intent(in), dimension(:) :: epsV
   
    tensor_spectral_index = hf_tensor_spectral_index(slowroll_to_hubble(epsV))

  end function tensor_spectral_index




  function tensor_to_scalar_ratio(epsV)
    use huflow, only : hf_tensor_to_scalar_ratio
    implicit none
    real(kp) :: tensor_to_scalar_ratio
    real(kp), intent(in), dimension(:) :: epsV
    
    tensor_to_scalar_ratio = hf_tensor_to_scalar_ratio(slowroll_to_hubble(epsV))
    
  end function tensor_to_scalar_ratio




  function scalar_running(epsV)
    use huflow, only : hf_scalar_running
    implicit none
    real(kp) :: scalar_running
    real(kp), intent(in), dimension(:) :: epsV

    scalar_running = hf_scalar_running(slowroll_to_hubble(epsV))

  end function scalar_running




  function tensor_running(epsV)
    use huflow, only : hf_tensor_running
    implicit none
    real(kp) :: tensor_running
    real(kp), intent(in), dimension(:) :: epsV
    
    tensor_running = hf_tensor_running(slowroll_to_hubble(epsV))

  end function tensor_running




  function scalar_running_running(epsV)
    use huflow, only : hf_scalar_running_running
    implicit none
    real(kp) :: scalar_running_running
    real(kp), intent(in), dimension(:) :: epsV
  
    scalar_running_running = hf_scalar_running_running(slowroll_to_hubble(epsV))

  end function scalar_running_running




  function tensor_running_running(epsV)
    use huflow, only : hf_tensor_running_running
    implicit none
    real(kp) :: tensor_running_running
    real(kp), intent(in), dimension(:) :: epsV
    
    tensor_running_running = hf_tensor_running_running(slowroll_to_hubble(epsV))

  end function tensor_running_running



!this gives the normal slow-roll expansion of P(k*)/[H*^2/(8pi^2 eps1*)]
  function slowroll_corrections(epsV)
    use huflow, only : hubbleflow_corrections
    implicit none
    real(kp) :: slowroll_corrections
    real(kp), dimension(:), intent(in) :: epsV
   
    slowroll_corrections = hubbleflow_corrections(slowroll_to_hubble(epsV))

  end function slowroll_corrections



!this gives [H*^2/(8pi^2 eps1*)]/P(k*) consitently expanded in slow-roll
!hence the name "inverse" as opposed to normal slow-roll expansion
  function inverse_slowroll_corrections(epsV)
    use huflow, only : inverse_hubbleflow_corrections
    implicit none
    real(kp) :: inverse_slowroll_corrections
    real(kp), dimension(:), intent(in) :: epsV
    
    inverse_slowroll_corrections = inverse_hubbleflow_corrections(slowroll_to_hubble(epsV))

  end function inverse_slowroll_corrections



!this gives ln[P(k*)] - ln[H*^2/(8pi^2 eps1*)] consistently expanded
  function ln_slowroll_corrections(epsV)
    use huflow, only : ln_hubbleflow_corrections
    implicit none
    real(kp) :: ln_slowroll_corrections
    real(kp), dimension(:), intent(in) :: epsV
    
    ln_slowroll_corrections = ln_hubbleflow_corrections(slowroll_to_hubble(epsV))

  end function ln_slowroll_corrections


!this gives ln[H*^2/(8pi^2 eps1*)] - ln[P(k*)] consistently expanded
  function ln_inverse_slowroll_corrections(epsV)
    use huflow, only : ln_inverse_hubbleflow_corrections
    implicit none
    real(kp) :: ln_inverse_slowroll_corrections
    real(kp), dimension(:), intent(in) :: epsV

    ln_inverse_slowroll_corrections = ln_inverse_hubbleflow_corrections(slowroll_to_hubble(epsV))

  end function ln_inverse_slowroll_corrections



end module srflow

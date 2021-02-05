!Common functions for the hybrid natural inflation potential
!
!V(phi) = M**4 [ 1 + alpha cos(x) )
!
!x = phi/f

module hni1sr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent

  use hnicommon, only : hni_check_params, hni_norm_potential
  use hnicommon, only : hni_norm_deriv_potential, hni_norm_deriv_second_potential
  use hnicommon, only : hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  use hnicommon, only : hni_x_epsoneunity, hni_alpha, hni_f
  use hnicommon, only : find_hni_x_trajectory, hni_efold_primitive
  use hnicommon, only : hni_numacc_xinimin
  
  implicit none

  private

  public hni1_norm_potential, hni1_norm_deriv_potential, hni1_norm_deriv_second_potential
  public hni1_epsilon_one, hni1_epsilon_two, hni1_epsilon_three
  public hni1_x_trajectory, hni1_x_endinf, hni1_check_params
  public hni1_alphamin, hni1_fmax, hni1_numacc_efoldmax



contains
      

  function hni1_alphamin(f)
    implicit none
    real(kp) :: hni1_alphamin
    real(kp), intent(in) :: f
    
    hni1_alphamin = hni_alpha(f)
    
  end function hni1_alphamin

  function hni1_fmax(alpha)
    implicit none
    real(kp) :: hni1_fmax
    real(kp), intent(in) :: alpha
    
    hni1_fmax = hni_f(alpha)
    
  end function hni1_fmax
  
  
  function hni1_check_params(alpha,f)
    implicit none
    logical :: hni1_check_params
    real(kp), intent(in) :: alpha,f

    hni1_check_params = hni_check_params(alpha,f) .and. (alpha.ge.hni1_alphamin(f))

  end function hni1_check_params


!returns V/M**4
  function hni1_norm_potential(x,alpha,f)
    implicit none
    real(kp) :: hni1_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni1_norm_potential = hni_norm_potential(x,alpha,f)

  end function hni1_norm_potential


!returns the first derivative of the potential with respect to x=phi/f, divided by M**4
  function hni1_norm_deriv_potential(x,alpha,f)
    implicit none
    real(kp) :: hni1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni1_norm_deriv_potential = hni_norm_deriv_potential(x,alpha,f)

  end function hni1_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/f, divided by M**4
  function hni1_norm_deriv_second_potential(x,alpha,f)
    implicit none
    real(kp) :: hni1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni1_norm_deriv_second_potential = hni_norm_deriv_second_potential(x,alpha,f)

  end function hni1_norm_deriv_second_potential


!epsilon1(x)
  function hni1_epsilon_one(x,alpha,f)    
    implicit none
    real(kp) :: hni1_epsilon_one
    real(kp), intent(in) :: x,alpha,f

    hni1_epsilon_one = hni_epsilon_one(x,alpha,f)    

  end function hni1_epsilon_one


!epsilon2(x)
  function hni1_epsilon_two(x,alpha,f)    
    implicit none
    real(kp) :: hni1_epsilon_two
    real(kp), intent(in) :: x,alpha,f

    hni1_epsilon_two = hni_epsilon_two(x,alpha,f) 

  end function hni1_epsilon_two

!epsilon3(x)
  function hni1_epsilon_three(x,alpha,f)    
    implicit none
    real(kp) :: hni1_epsilon_three
    real(kp), intent(in) :: x,alpha,f

    hni1_epsilon_three = hni_epsilon_three(x,alpha,f)

  end function hni1_epsilon_three


  !returns x at the end of inflation defined as epsilon1=1
  function hni1_x_endinf(alpha,f)
    implicit none
    real(kp) :: hni1_x_endinf
    real(kp), intent(in) :: alpha,f
    real(kp), dimension(2) :: xepsone

    xepsone = hni_x_epsoneunity(alpha,f)
    
    hni1_x_endinf = xepsone(1)

  end function hni1_x_endinf



  function hni1_numacc_efoldmax(alpha,f)
    implicit none
    real(kp) :: hni1_numacc_efoldmax
    real(kp), intent(in) :: alpha,f

    real(kp) :: xinimin, xend

    xend = hni1_x_endinf(alpha,f)
    xinimin = hni_numacc_xinimin(alpha,f)


    hni1_numacc_efoldmax = -hni_efold_primitive(xend,alpha,f) &
         + hni_efold_primitive(xinimin,alpha,f)


  end function hni1_numacc_efoldmax
  

!returns x at bfold=-efolds before the end of inflation
  function hni1_x_trajectory(bfold,xend,alpha,f)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,f
    real(kp) :: hni1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: hni1Data

    mini = epsilon(1._kp) !potential maximum
    maxi = xend

    hni1Data%real1 = alpha
    hni1Data%real2 = f
    hni1Data%real3 = -bfold + hni_efold_primitive(xend,alpha,f)

    hni1_x_trajectory = zbrent(find_hni_x_trajectory,mini,maxi,tolFind,hni1Data)

  end function hni1_x_trajectory



end module hni1sr

!Common functions for the hybrid natural inflation potential
!
!V(phi) = M**4 [ 1 + alpha cos(x) )
!
!x = phi/f

module hnicommon
  use infprec, only : kp, tolkp,transfert, pi
  use inftools, only : zbrent
  implicit none

  private

  public hni_norm_potential, hni_norm_deriv_potential, hni_norm_deriv_second_potential
  public hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  public hni_efold_primitive, find_hni_x_trajectory, hni_x_epsoneunity
  public hni_alpha, hni_f, hni_check_params, hni_x_potmin, hni_numacc_x_epsonenull
  public hni_numacc_xinimin
  
  real(kp), parameter :: hniSmall = epsilon(1._kp)

  public hniSmall
  
contains


  function hni_check_params(alpha,f)
    implicit none
    logical :: hni_check_params
    real(kp), intent(in) :: alpha,f

    hni_check_params = (alpha .gt. 0._kp).and.(alpha.lt.1._kp)

  end function hni_check_params


!Minimum value of alpha (depending on f) so that inflation ends naturally
  function hni_alpha(f)
    implicit none
    real(kp) :: hni_alpha
    real(kp), intent(in) :: f

    hni_alpha = 1._kp/sqrt(1._kp+0.5_kp/f**2)

  end function hni_alpha

!Maximum value of f (depending on alpha) so that inflation ends naturally
  function hni_f(alpha)
    implicit none
    real(kp) :: hni_f
    real(kp), intent(in) :: alpha

    hni_f = alpha/(sqrt(2._kp*(1._kp - alpha**2)))

  end function hni_f

!Returns x at the bottom of the potential (first minimum, Vmin >=0 for alpha <=1)
  function hni_x_potmin(alpha,f)
    implicit none
    real(kp) :: hni_x_potmin
    real(kp), intent(in) :: alpha,f

    hni_x_potmin = pi

  end function hni_x_potmin
  

!returns V/M**4
  function hni_norm_potential(x,alpha,f)
    implicit none
    real(kp) :: hni_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni_norm_potential = 1._kp+alpha*cos(x)

  end function hni_norm_potential


!returns the first derivative of the potential with respect to x=phi/f, divided by M**4
  function hni_norm_deriv_potential(x,alpha,f)
    implicit none
    real(kp) :: hni_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni_norm_deriv_potential = -alpha*sin(x)

  end function hni_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/f, divided by M**4
  function hni_norm_deriv_second_potential(x,alpha,f)
    implicit none
    real(kp) :: hni_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: f

    hni_norm_deriv_second_potential = -alpha*cos(x)

  end function hni_norm_deriv_second_potential


!epsilon1(x)
  function hni_epsilon_one(x,alpha,f)    
    implicit none
    real(kp) :: hni_epsilon_one
    real(kp), intent(in) :: x,alpha,f

    hni_epsilon_one = ((alpha**2*sin(x)**2)/(2._kp*(1._kp+alpha*cos(x))**2))/f**2

  end function hni_epsilon_one


!epsilon2(x)
  function hni_epsilon_two(x,alpha,f)    
    implicit none
    real(kp) :: hni_epsilon_two
    real(kp), intent(in) :: x,alpha,f

    hni_epsilon_two = 2._kp/f**2*(alpha*(alpha+cos(x)))/(1._kp+alpha*cos(x))**2

  end function hni_epsilon_two

!epsilon3(x)
  function hni_epsilon_three(x,alpha,f)    
    implicit none
    real(kp) :: hni_epsilon_three
    real(kp), intent(in) :: x,alpha,f

    hni_epsilon_three = ((alpha*(-1._kp+2._kp*alpha**2+alpha*cos(x))* &
                        sin(x)**2)/((alpha+cos(x))*(1._kp+alpha*cos(x))**2))/f**2

  end function hni_epsilon_three

!returns the two solutions of epsilon1=1 in [0,pi]
  function hni_x_epsoneunity(alpha,f)
    implicit none
    real(kp), dimension(2) :: hni_x_epsoneunity
    real(kp), intent(in) :: alpha,f


    !check the solution exist
    if (alpha.lt.hni_alpha(f)) then
       write(*,*)'alpha= f= ',alpha,f
       stop 'hni_x_epsoneunity: no solution to epsilon1=1!'
    endif
    
    hni_x_epsoneunity(1) = acos( ( sqrt(alpha**2/(4._kp*f**2)+alpha**2/2._kp-0.5_kp)/f- &
         1._kp)/(alpha*(1._kp+1._kp/(2._kp*f**2))))

    hni_x_epsoneunity(2) = acos( ( -sqrt(alpha**2/(4._kp*f**2)+alpha**2/2._kp-0.5_kp)/f- &
                  1._kp)/(alpha*(1._kp+1._kp/(2._kp*f**2))))

  end function hni_x_epsoneunity


!returns the two solutions of epsilon1=epsMini in [0,pi]
  function hni_numacc_x_epsonenull(alpha,f)
    implicit none
    real(kp), dimension(2) ::  hni_numacc_x_epsonenull
    real(kp), intent(in) :: alpha,f

    real(kp) :: sqrteps

    sqrteps = sqrt(hniSmall)
    
!let's be clever
    hni_numacc_x_epsonenull = hni_x_epsoneunity(alpha,f*sqrteps)

  end function hni_numacc_x_epsonenull



  !return the minimal positive value of xini for ensuring eps1 >
!numerical accuracy, this is 0 + smallterm
  function hni_numacc_xinimin(alpha,f)
    implicit none
    real(kp) :: hni_numacc_xinimin
    real(kp), intent(in) :: alpha,f

    real(kp), dimension(2) :: xnumacc
      
    xnumacc = hni_numacc_x_epsonenull(alpha,f)

!keep the smallest one (close to the top of the potential)
    hni_numacc_xinimin = max(epsilon(1._kp), xnumacc(1))

    
    
  end function hni_numacc_xinimin

  

!this is integral[V(phi)/V'(phi) dphi]
  function hni_efold_primitive(x,alpha,f)
    implicit none
    real(kp), intent(in) :: x,alpha,f
    real(kp) :: hni_efold_primitive, test

!-f**2/alpha*(log(tan(x/2._kp))+alpha*log(sin(x)))

    hni_efold_primitive = f**2/alpha*((1._kp-alpha)*log(cos(x/2._kp))-(1._kp+alpha)*log(sin(x/2))) &
         - f**2*log(2._kp)
    
  end function hni_efold_primitive




  function find_hni_x_trajectory(x,hniData)
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: hniData
    real(kp) :: find_hni_x_trajectory
    real(kp) :: alpha,f,NplusPrimEnd

    alpha = hniData%real1
    f = hniData%real2
    NplusPrimEnd = hniData%real3

    find_hni_x_trajectory = hni_efold_primitive(x,alpha,f) - NplusPrimEnd

  end function find_hni_x_trajectory





end module hnicommon

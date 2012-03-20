!slow-roll functions for the pseudo natural inflation potential
!
!V(phi) = M**4 ( 1 - alpha ln(cos(phi/mu)) )
!
!x = phi/mu

module psnisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : polylog
  implicit none

  private

  public  psni_norm_potential, psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  public  psni_x_endinf, psni_efold_primitive, psni_x_trajectory
  public  psni_norm_deriv_potential, psni_norm_deriv_second_potential

 
contains
!returns V/M**4 as function of x=phi/mu
  function psni_norm_potential(x,alpha)
    implicit none
    real(kp) :: psni_norm_potential
    real(kp), intent(in) :: x,alpha

    psni_norm_potential = 1._kp-alpha*log(cos(x))

  end function psni_norm_potential



!returns the first derivative of the potential with respect to phi/Mp, divided by M**4
  function psni_norm_deriv_potential(x,alpha,mu)
    implicit none
    real(kp) :: psni_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,mu

   psni_norm_deriv_potential = -alpha*tan(x)/mu

  end function psni_norm_deriv_potential



!returns the second derivative of the potential with respect to phi/Mp, divided by M**4
  function psni_norm_deriv_second_potential(x,alpha,mu)
    implicit none
    real(kp) :: psni_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,mu

    psni_norm_deriv_second_potential = -alpha/(cos(x)*mu)**2

  end function psni_norm_deriv_second_potential



!epsilon_one(x)
  function psni_epsilon_one(x,alpha,mu)    
    implicit none
    real(kp) :: psni_epsilon_one
    real(kp), intent(in) :: x,alpha,mu
    
    psni_epsilon_one =0.5_kp*(alpha/mu*tan(x)/(1._kp+alpha*log(cos(x))))**2
    
  end function psni_epsilon_one


!epsilon_two(x)
  function psni_epsilon_two(x,alpha,mu)    
    implicit none
    real(kp) :: psni_epsilon_two
    real(kp), intent(in) :: x,alpha,mu
    
    psni_epsilon_two =(2._kp*alpha*((1._kp+alpha*log(cos(x)))/(cos(x)**2)+ &
                      alpha*tan(x)**2))/(1._kp+alpha*log(cos(x)))**2
    
  end function psni_epsilon_two


!epsilon_three(x)
  function psni_epsilon_three(x,alpha,mu)    
    implicit none
    real(kp) :: psni_epsilon_three
    real(kp), intent(in) :: x,alpha,mu
    
    psni_epsilon_three = (alpha*(2._kp+3._kp*alpha+alpha**2-alpha**2*cos(2._kp*x)+ &
                         alpha*(4._kp+3._kp*alpha)*log(cos(x))+2._kp*alpha**2*log(cos(x))**2)/ &
                         (cos(x)**2)*tan(x)**2)/((1._kp+alpha*log(cos(x)))**2*((1._kp+alpha*log(cos(x)))/ &
                         (cos(x)**2)+alpha*tan(x)**2))
    
  end function psni_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function psni_x_endinf(alpha,mu)
    implicit none
    real(kp), intent(in) :: alpha,mu
    real(kp) :: psni_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: psniData

    mini = epsilon(1._kp)
    maxi = acos(-1._kp)/2._kp*(1._kp-epsilon(1._kp))
  

    psniData%real1 = alpha
    psniData%real2 = mu
    
    psni_x_endinf = zbrent(find_psni_x_endinf,mini,maxi,tolFind,psniData)
   
   
  end function psni_x_endinf

  function find_psni_x_endinf(x,psniData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: psniData
    real(kp) :: find_psni_x_endinf
    real(kp) :: alpha,mu

    alpha = psniData%real1
    mu = psniData%real2

    find_psni_x_endinf = psni_epsilon_one(x,alpha,mu)-1._kp
   
  end function find_psni_x_endinf


!this is integral(V(phi)/V'(phi) dphi)
  function psni_efold_primitive(x,alpha,mu)
    implicit none
    real(kp), intent(in) :: x,alpha,mu
    real(kp) :: psni_efold_primitive

    if (alpha.eq.0._kp) stop 'psni_efold_primitive: alpha=0 !'

    psni_efold_primitive = -mu**2/alpha*((1._kp+alpha*log(cos(x)))*log(sin(x))+ &
                            0.25_kp*alpha*real(polylog(complex(cos(x)**2,0._kp),complex(2._kp,0._kp)),kp))


  end function psni_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function psni_x_trajectory(bfold,xend,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold,alpha,mu,xend
    real(kp) :: psni_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: psniData

    mini=epsilon(1._kp)
    maxi = xend

    psniData%real1 = alpha
    psniData%real2 = mu
    psniData%real3 = -bfold + psni_efold_primitive(xend,alpha,mu)
    
    psni_x_trajectory = zbrent(find_psnitraj,mini,maxi,tolFind,psniData)
       
  end function psni_x_trajectory

  function find_psnitraj(x,psniData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: psniData
    real(kp) :: find_psnitraj
    real(kp) :: alpha,mu,NplusNuend

    alpha= psniData%real1
    mu = psniData%real2
    NplusNuend = psniData%real3

    find_psnitraj = psni_efold_primitive(x,alpha,mu) - NplusNuend
   
  end function find_psnitraj


  
end module psnisr

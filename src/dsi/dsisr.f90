!slow-roll functions for the dynamical supersymmetrical inflation potential
!
!V(phi) = M^4 [ 1 + (phi/mu)^(-p) ]
!
!x = phi/mu
!mu=mu/Mp

module dsisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : hypergeom_2F1
  implicit none

  private

  public  dsi_norm_potential, dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three
  public  dsi_x_min, dsi_efold_primitive, dsi_x_trajectory
  public  dsi_norm_deriv_potential, dsi_norm_deriv_second_potential

 
contains
!returns V/M**4 as function of x=phi/mu
  function dsi_norm_potential(x,p)
    implicit none
    real(kp) :: dsi_norm_potential
    real(kp), intent(in) :: x,p

    dsi_norm_potential = 1._kp+x**(-p)

  end function dsi_norm_potential



!returns the first derivative of the potential with respect to x=phi/mu, divided by M**4
  function dsi_norm_deriv_potential(x,p)
    implicit none
    real(kp) :: dsi_norm_deriv_potential
    real(kp), intent(in) :: x,p

   dsi_norm_deriv_potential = -p/(x**(p+1._kp))

  end function dsi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/mu, divided by M**4
  function dsi_norm_deriv_second_potential(x,p)
    implicit none
    real(kp) :: dsi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p

    dsi_norm_deriv_second_potential = p*(p+1._kp)/(x**(p+2._kp))

  end function dsi_norm_deriv_second_potential



!epsilon_one(x)
  function dsi_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: dsi_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    dsi_epsilon_one = p**2/(2._kp*mu**2)*x**(-2._kp*p-2._kp)/((1+x**(-p))**2)
    
  end function dsi_epsilon_one


!epsilon_two(x)
  function dsi_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: dsi_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    dsi_epsilon_two =-2._kp*p/(mu**2)*x**(-p-2._kp)*(x**(-p)+p+1._kp)/((1+x**(-p))**2)
    
  end function dsi_epsilon_two


!epsilon_three(x)
  function dsi_epsilon_three(x,p,mu)    
    implicit none
    real(kp) :: dsi_epsilon_three
    real(kp), intent(in) :: x,p,mu
    
    dsi_epsilon_three = -p/(mu**2)*x**(-p-2._kp)/ &
                        ((1._kp+x**(-p))**2*(x**(-p)+p+1._kp))* &
                        (2._kp*x**(-2._kp)*(p+1._kp)*(p-4._kp)*x**(-p)+(p+1._kp)*(p+2._kp))
    
  end function dsi_epsilon_three


!returns the minimum value for x, defined as epsilon1=1, such that inflation takes place
  function dsi_x_min(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: dsi_x_min
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: dsiData

    mini = epsilon(1._kp)
    maxi = (p/(sqrt(2._kp)*mu))**(1._kp/(p+1._kp))*100._kp !Asymptotic Solution when x_min >> 1, times some safety number
  

    dsiData%real1 = p
    dsiData%real2 = mu
    
    dsi_x_min = zbrent(find_dsi_x_min,mini,maxi,tolFind,dsiData)
   
   
  end function dsi_x_min

  function find_dsi_x_min(x,dsiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: dsiData
    real(kp) :: find_dsi_x_min
    real(kp) :: p,mu

    p = dsiData%real1
    mu = dsiData%real2

    find_dsi_x_min = dsi_epsilon_one(x,p,mu)-1._kp
   
  end function find_dsi_x_min


!this is integral(V(phi)/V'(phi) dphi)
  function dsi_efold_primitive(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: dsi_efold_primitive

    if (p.eq.0._kp) stop 'dsi_efold_primitive: p=0 !'

    dsi_efold_primitive = -mu**2/(2._kp*p)*(x**2+2._kp/(p+2._kp)*x**(p+2._kp))

  end function dsi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function dsi_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold,p,mu,xend
    real(kp) :: dsi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: dsiData

    mini=dsi_x_min(p,mu)
    maxi = xend

    dsiData%real1 = p
    dsiData%real2 = mu
    dsiData%real3 = -bfold + dsi_efold_primitive(xend,p,mu)
    
    dsi_x_trajectory = zbrent(find_dsitraj,mini,maxi,tolFind,dsiData)
       
  end function dsi_x_trajectory

  function find_dsitraj(x,dsiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: dsiData
    real(kp) :: find_dsitraj
    real(kp) :: p,mu,NplusNuend

    p= dsiData%real1
    mu = dsiData%real2
    NplusNuend = dsiData%real3

    find_dsitraj = dsi_efold_primitive(x,p,mu) - NplusNuend
   
  end function find_dsitraj


  
end module dsisr

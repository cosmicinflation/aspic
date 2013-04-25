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
  use cosmopar, only : powerAmpScalar
  implicit none

  private

  public dsi_norm_potential, dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three
  public dsi_efold_primitive, dsi_x_trajectory
  public dsi_norm_deriv_potential, dsi_norm_deriv_second_potential
  public dsi_xinimin, dsi_xendmin, dsi_x_epsoneunity, dsi_xendmax, dsi_mumax

    
 
contains
!returns V/M**4 as function of x=phi/mu
  function dsi_norm_potential(x,p,mu)
    implicit none
    real(kp) :: dsi_norm_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    dsi_norm_potential = 1._kp+x**(-p)

  end function dsi_norm_potential



!returns the first derivative of the potential with respect to x=phi/mu, divided by M**4
  function dsi_norm_deriv_potential(x,p,mu)
    implicit none
    real(kp) :: dsi_norm_deriv_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    dsi_norm_deriv_potential = -p/(x**(p+1._kp))

  end function dsi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/mu, divided by M**4
  function dsi_norm_deriv_second_potential(x,p,mu)
    implicit none
    real(kp) :: dsi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

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


!returns the field value at which eps1=1
  function dsi_x_epsoneunity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: dsi_x_epsoneunity
    real(kp) :: xeps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: dsiData


    dsiData%real1 = p
    dsiData%real2 = mu

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    dsiData%real1 = p
    dsiData%real2 = mu
    
    dsi_x_epsoneunity = zbrent(find_dsi_x_epsoneunity,mini,maxi,tolFind,dsiData)
    
  end function dsi_x_epsoneunity
  
  function find_dsi_x_epsoneunity(x,dsiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: dsiData
    real(kp) :: find_dsi_x_epsoneunity
    real(kp) :: p,mu

    p = dsiData%real1
    mu = dsiData%real2
    
    find_dsi_x_epsoneunity = dsi_epsilon_one(x,p,mu) - 1._kp
  
  end function find_dsi_x_epsoneunity



!the min value for xini
  function dsi_xinimin(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: dsi_xinimin
    
    dsi_xinimin = dsi_x_epsoneunity(p,mu)
    
  end function dsi_xinimin



!the min value for xend if we want efolds inflation
  function dsi_xendmin(efold,p,mu)
    implicit none
    real(kp), intent(in) :: efold, p, mu
    real(kp) :: dsi_xendmin, xini

    xini = dsi_xinimin(p,mu)    

    dsi_xendmin = dsi_x_trajectory(efold,xini,p,mu)

  end function dsi_xendmin

!the max value for xend in order for an extra non renormalizable term x^q in the potential to have no effect. The result depends on q.
  function dsi_xendmax(p,mu,q)
    implicit none
    real(kp), intent(in) :: p,mu,q
    real(kp) :: dsi_xendmax, xini

    xini = dsi_xinimin(p,mu)    

    dsi_xendmax = (720._kp*acos(-1._kp)**2*p**3/(q+4._kp)* &
         mu**(-q-6._kp)*60._kp*powerAmpScalar)**(1._kp/(3._kp*p+q+6._kp))

  end function dsi_xendmax

!the max value for mu in order to have dsi_xendmin<dsi_xendmax
  function dsi_mumax(p,q,efolds)
    implicit none
    real(kp), intent(in) :: p,q,efolds
    real(kp) :: dsi_mumax  

    dsi_mumax = (720._kp*acos(-1._kp)**2*p**3/(q+4._kp)* &
         60._kp*powerAmpScalar)**((p+2._kp)/(p*q))/ &
         ((p*(p+2._kp)*efolds)**((3._kp*p+q+6._kp)/(p*q)))

  end function dsi_mumax

 


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
    real(kp) :: dsi_x_trajectory,xiniMin
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: dsiData

    xiniMin = dsi_xinimin(p,mu)

    if (bfold.gt.0._kp) then
       mini=xiniMin
       maxi = 1._kp/epsilon(1._kp)
    else
       mini= xiniMin
       maxi = xend
    endif

    dsiData%real1 = p
    dsiData%real2 = mu
    dsiData%real3 = -bfold + dsi_efold_primitive(xend,p,mu)
    
    dsi_x_trajectory = zbrent(find_dsi_x_trajectory,mini,maxi,tolFind,dsiData)
       
  end function dsi_x_trajectory

  function find_dsi_x_trajectory(x,dsiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: dsiData
    real(kp) :: find_dsi_x_trajectory
    real(kp) :: p,mu,NplusNuend

    p= dsiData%real1
    mu = dsiData%real2
    NplusNuend = dsiData%real3

    find_dsi_x_trajectory = dsi_efold_primitive(x,p,mu) - NplusNuend
   
  end function find_dsi_x_trajectory


  
end module dsisr

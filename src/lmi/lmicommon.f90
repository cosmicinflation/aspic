!slow-roll functions for the Logamediate inflation 1 and 2 potential
!("1" means that inflation proceeds from the right to the left)
!
!V(phi) = M^4 x^[4(1-gam)] exp[-beta x^gam]
!
!x = (phi-phi0)/Mp


module lmicommon
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  implicit none

  private

  public lmi_alpha, lmi_norm_potential
  public lmi_norm_deriv_potential, lmi_norm_deriv_second_potential
  public lmi_epsilon_one, lmi_epsilon_two, lmi_epsilon_three
  public lmi_epsilon_one_max, lmi_x_epsonemax, lmi_x_potmax  
  public lmi_efold_primitive, find_lmi_x_trajectory
  public lmi_epstwo_potmax, lmi_x_epstwounity
  public lmi_numacc_x_big, lmi_numacc_dx_potmax

contains



  function lmi_alpha(gam)
    implicit none
    real(kp), intent(in) :: gam
    real(kp) :: lmi_alpha

    lmi_alpha = 4._kp*(1._kp-gam)

  end function lmi_alpha


!returns V/M^4
  function lmi_norm_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi_norm_potential
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)

    lmi_norm_potential = x**alpha*exp(-beta*x**gam)

  end function lmi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function lmi_norm_deriv_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi_norm_deriv_potential
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)

    lmi_norm_deriv_potential = (alpha*x**(alpha-1._kp)- &
         beta*gam*x**(alpha+gam-1._kp))*exp(-beta*x**gam)

  end function lmi_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function lmi_norm_deriv_second_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)
   
    lmi_norm_deriv_second_potential = (x**(-2._kp + alpha)*((-1._kp + alpha)*alpha &
         - beta*gam*(-1._kp+2._kp*alpha + gam)*x**gam &
         + beta**2._kp*gam**2._kp*x**(2._kp*gam)))*exp(-beta*x**gam)
    

  end function lmi_norm_deriv_second_potential



!epsilon_one(x)
  function lmi_epsilon_one(x,gam,beta)    
    implicit none
    real(kp) :: lmi_epsilon_one
    real(kp), intent(in) :: x,gam,beta

    real(kp) ::alpha
    alpha = lmi_alpha(gam)

    lmi_epsilon_one = (alpha-beta*gam*x**gam)**2._kp/(2._kp*x**2._kp)
    
  end function lmi_epsilon_one


!epsilon_two(x)
  function lmi_epsilon_two(x,gam,beta)    
    implicit none
    real(kp) :: lmi_epsilon_two
    real(kp), intent(in) :: x,gam,beta

    real(kp) ::alpha
    alpha = lmi_alpha(gam)
    
    lmi_epsilon_two = 2._kp*(alpha+beta*(-1._kp+gam) &
         *gam*x**gam)/(x**2._kp)
    
  end function lmi_epsilon_two


!epsilon_three(x)
  function lmi_epsilon_three(x,gam,beta)    
    implicit none
    real(kp) :: lmi_epsilon_three
    real(kp), intent(in) :: x,gam,beta

    real(kp) ::alpha

    alpha = lmi_alpha(gam)
    
    lmi_epsilon_three = ((alpha-beta*gam*x**gam)*(2._kp*alpha- & 
         beta*(-2._kp+gam)*(-1._kp+gam)*gam* &
         x**gam))/(x**2._kp*(alpha+beta*(-1._kp+gam)* & 
         gam*x**gam))
    
  end function lmi_epsilon_three

!the maximal value of eps1 in the domain x>xvmax
  function lmi_epsilon_one_max(gam,beta)
    implicit none
    real(kp) :: lmi_epsilon_one_max
    real(kp), intent(in) :: gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)

    lmi_epsilon_one_max= 0.5_kp*(alpha*gam/(1._kp-gam))**2 &
         *(beta*gam*(1._kp-gam)/alpha)**(2._kp/gam)

  end function lmi_epsilon_one_max


!the x value at which eps1 is maximal in the domain x > xvmax
  function lmi_x_epsonemax(gam,beta)
    implicit none
    real(kp) :: lmi_x_epsonemax
    real(kp), intent(in) :: gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)

    lmi_x_epsonemax = (alpha/(beta*gam)/(1._kp-gam))**(1._kp/gam)

  end function lmi_x_epsonemax


!xvmax, the x value at which the potential is maximal
  function lmi_x_potmax(gam,beta)
    implicit none
    real(kp) :: lmi_x_potmax
    real(kp), intent(in) :: gam,beta
    real(kp) :: alpha

    alpha = lmi_alpha(gam)

    lmi_x_potmax = (alpha/(beta*gam))**(1._kp/gam)

  end function lmi_x_potmax



!eps2 at the top of the potential (xVmax)
  function lmi_epstwo_potmax(gam,beta)
    implicit none
    real(kp) :: lmi_epstwo_potmax
    real(kp), intent(in) :: gam,beta
    real(kp) :: xVmax

    xVmax = lmi_x_potmax(gam,beta)

    lmi_epstwo_potmax = lmi_epsilon_two(xVmax,gam,beta)

  end function lmi_epstwo_potmax




!returns the value xeps2 at which eps2=1 (it always exists)
  function lmi_x_epstwounity(gam,beta)
    implicit none
    real(kp), intent(in) :: gam, beta
    real(kp) :: lmi_x_epstwounity
    real(kp) :: alpha
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmiData

    if (gam.gt.1._kp) then
       write(*,*)'gamma= ',gam
       stop 'lmi_x_epstwounity: gam>1'
    endif

    alpha = lmi_alpha(gam)

    mini = epsilon(1._kp)
    maxi = lmi_x_epsonemax(gam,beta)
       
    lmiData%real1 = gam
    lmiData%real2 = beta
    
    lmi_x_epstwounity &
         = zbrent(find_lmi_x_epstwounity,mini,maxi,tolFind,lmiData)
        
  end function lmi_x_epstwounity



  function find_lmi_x_epstwounity(x,lmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmiData
    real(kp) :: find_lmi_x_epstwounity
    real(kp) :: gam,beta

    gam = lmiData%real1
    beta = lmiData%real2
    
    find_lmi_x_epstwounity = lmi_epsilon_two(x,gam,beta) - 1._kp
   
  end function find_lmi_x_epstwounity




!this is integral[V(phi)/V'(phi) dphi]
  function lmi_efold_primitive(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi_efold_primitive

    real(kp) ::alpha
    alpha = lmi_alpha(gam)

    if (alpha.eq.0._kp) stop 'lmi_efold_primitive: gam=1!  (PLI)'
    if (gam.eq.0._kp) stop 'lmi_efold_primitive: gam=0!'

    lmi_efold_primitive = x**2._kp/(2._kp*alpha)*hypergeom_2F1(1._kp,2._kp/gam, &
         2._kp/gam+1._kp,beta*gam/alpha*x**gam)

  end function lmi_efold_primitive



  function find_lmi_x_trajectory(x,lmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmiData
    real(kp) :: find_lmi_x_trajectory
    real(kp) :: gam,beta,NplusNuend

    gam = lmiData%real1
    beta = lmiData%real2
    NplusNuend = lmiData%real3

    find_lmi_x_trajectory = lmi_efold_primitive(x,gam,beta) - NplusNuend
   
  end function find_lmi_x_trajectory


!we use exp(-beta x^gam), the argument cannot be too large
  function lmi_numacc_x_big(gam,beta)
    implicit none
    real(kp) :: lmi_numacc_x_big
    real(kp), intent(in) :: gam,beta
    real(kp), parameter :: big = epsilon(1._kp)*huge(1._kp)

    if (gam.eq.0._kp) stop 'lmi_numacc_x_big: gam=0!'

    lmi_numacc_x_big = (log(big)/beta)**(1._kp/gam)

  end function lmi_numacc_x_big

!at x=xVmax, eps1=0 => cannot be integrated. This function returns dx
!such that at x=xVmax+-dx, eps1 = min(numprec,lmi_epsonemax)
  function lmi_numacc_dx_potmax(gam,beta)
    implicit none
    real(kp) :: lmi_numacc_dx_potmax
    real(kp), intent(in) :: gam,beta
    real(kp) :: alpha,xVmax,delta
    real(kp) :: epsmin
    
    alpha = lmi_alpha(gam)

    xVmax = lmi_x_potmax(gam,beta)
    epsmin = min(epsilon(1._kp), lmi_epsilon_one_max(gam,beta))

    delta = xVmax/(alpha*gam)*sqrt(2*epsmin)
    
    if (delta.ge.1._kp) then
       print *,'test',gam,beta
       write(*,*)'xVmax= dxVmax= ',xVmax, delta*xVmax,epsmin
       stop 'lmi_numacc_dx_potmax: xVmax too large to compute dxVmax!'
    endif
        
    lmi_numacc_dx_potmax = delta*xVmax

  end function lmi_numacc_dx_potmax

end module lmicommon

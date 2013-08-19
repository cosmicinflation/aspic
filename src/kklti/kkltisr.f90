!slow-roll functions for the exact KKLT brane inflation potential
!
!V(phi) = M**4 / [1 + x**(-p)]
!
!x = phi/mu
!mu = mu/Mp

module kkltisr
  use infprec, only : pi, kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public kklti_norm_potential, kklti_norm_deriv_potential, kklti_norm_deriv_second_potential
  public kklti_epsilon_one, kklti_epsilon_two,kklti_epsilon_three, kklti_x_epstwounity
  public kklti_x_epsoneunity, kklti_efold_primitive, kklti_x_trajectory
  public kklti_ln_xstg, kklti_ln_xuv
 
contains
!returns V/M**4
  function kklti_norm_potential(x,p,mu)
    implicit none
    real(kp) :: kklti_norm_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in) :: mu

    kklti_norm_potential = 1._kp/(1._kp + x**(-p))
  end function kklti_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function kklti_norm_deriv_potential(x,p,mu)
    implicit none
    real(kp) :: kklti_norm_deriv_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

   kklti_norm_deriv_potential = p*x**(-1._kp-p)/(1+x**(-p))**2

  end function kklti_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function kklti_norm_deriv_second_potential(x,p,mu)
    implicit none
    real(kp) :: kklti_norm_deriv_second_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    kklti_norm_deriv_second_potential = 2*p**2*x**(-2._kp-2*p)/(1._kp+x**(-p))**3 &
         - p*(1._kp+p)*x**(-2._kp-p)/(1._kp+x**(-p))**2

  end function kklti_norm_deriv_second_potential

!epsilon1(x)
  function kklti_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: kklti_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    kklti_epsilon_one = 0.5_kp*(p/(x+x**(p+1._kp)))**2/mu**2
    
  end function kklti_epsilon_one


!epsilon2(x)
  function kklti_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: kklti_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    kklti_epsilon_two = 2._kp*p*(1._kp + x**p*(p+1._kp)) &
         / (x+x**(1._kp+p))**2 /mu**2

  end function kklti_epsilon_two

!epsilon3(x)
  function kklti_epsilon_three(x,p,mu)    
    implicit none
    real(kp) :: kklti_epsilon_three
    real(kp), intent(in) :: x,p,mu
    
    kklti_epsilon_three = (p*(2._kp + (1._kp + p)*x**p*(4._kp - p + (2._kp + p)*x**p))) &
         /(x**2*(1._kp + x**p)**2*(1._kp + (1._kp + p)*x**p)) / mu**2
    
  end function kklti_epsilon_three


!this is integral[V(phi)/V'(phi) dphi]
  function kklti_efold_primitive(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: kklti_efold_primitive

    if ((p.eq.0._kp).or.(p.eq.-2._kp)) stop 'kklti_efold_primitivec: p=0/-2 is singular'

    kklti_efold_primitive = 0.5_kp * mu**2 * x**2 * (2._kp*x**p+p+2._kp)/p/(2._kp+p)

  end function kklti_efold_primitive


!returns x at epsilon1=1
  function kklti_x_epsoneunity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: kklti_x_epsoneunity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltiData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltiData%real1 = p
    kkltiData%real2 = mu

    kklti_x_epsoneunity = zbrent(find_kklti_x_epsoneunity,mini,maxi,tolFind,kkltiData)
   
  end function kklti_x_epsoneunity
  
  function find_kklti_x_epsoneunity(x,kkltiData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kkltiData
    real(kp) :: find_kklti_x_epsoneunity
    real(kp) :: p,mu
    
    p=kkltiData%real1
    mu=kkltiData%real2
!this is epsilon1(x)=1
    find_kklti_x_epsoneunity = p-sqrt(2._kp)*mu*(x+x**(1._kp+p))
    
  end function find_kklti_x_epsoneunity
 

!returns x at epsilon2=1
  function kklti_x_epstwounity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: kklti_x_epstwounity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltiData

    mini = epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    kkltiData%real1 = p
    kkltiData%real2 = mu

    kklti_x_epstwounity = zbrent(find_kklti_x_epstwounity,mini,maxi,tolFind,kkltiData)
   
  end function kklti_x_epstwounity
  
  function find_kklti_x_epstwounity(x,kkltiData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: kkltiData
    real(kp) :: find_kklti_x_epstwounity
    real(kp) :: p,mu
    
    p=kkltiData%real1
    mu=kkltiData%real2
!this is epsilon2(x)=1
    find_kklti_x_epstwounity = 2._kp*p*(1._kp+(p+1._kp)*x**p) &
         - (mu*(x+x**(p+1._kp)))**2
    
  end function find_kklti_x_epstwounity
 



!returns x at bfold=-efolds before the end of inflation
  function kklti_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: kklti_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: kkltiData

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    kkltiData%real1 = p
    kkltiData%real2 = mu
    kkltiData%real3 = -bfold + kklti_efold_primitive(xend,p,mu)
    
    kklti_x_trajectory = zbrent(find_kklti_x_trajectory,mini,maxi,tolFind,kkltiData)
    
  end function kklti_x_trajectory

  function find_kklti_x_trajectory(x,kkltiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: kkltiData
    real(kp) :: find_kklti_x_trajectory
    real(kp) :: p,mu,NplusPrimEnd

    p=kkltiData%real1
    mu = kkltiData%real2
    NplusPrimEnd = kkltiData%real3

    find_kklti_x_trajectory = kklti_efold_primitive(x,p,mu) - NplusPrimEnd
   
  end function find_kklti_x_trajectory

  

!return xstg in terms of the universal string variable y,vbar, flux
!number calN et string parama alpha'
  function kklti_ln_xstg(p,mu,y,vbar,calN,alpha)
    implicit none
    real(kp) :: kklti_ln_xstg
    real(kp), intent(in) :: p,mu,y,vbar,calN,alpha

    if (p.ne.4._kp) stop 'kklti_ln_xstg: p >< 4!'

    kklti_ln_xstg = 0.25_kp*calN + y**(-0.25_kp)
    
  end function kklti_ln_xstg

!return xUV in terms of the universal string variable y,vbar, flux
!number calN et string parama alpha'
  function kklti_ln_xuv(p,mu,y,vbar,calN,alpha)
    implicit none
    real(kp) :: kklti_ln_xuv
    real(kp), intent(in) :: p,mu,y,vbar,calN,alpha

    if (p.ne.4._kp) stop 'kklti_ln_xuv: p >< 4!'

    kklti_ln_xuv = -log(mu) + 0.75_kp*log(y) + 0.5_kp*log(vbar) &
         -log(pi) - 0.5_kp*log(2._kp*alpha*calN)

    if (kklti_ln_xuv.gt.log(2._kp)-log(mu)) then
       write(*,*)'kklti_ln_xuv: throat volume larger than bulk volume!'
       write(*,*)'phiUV= phiUVmax= ',exp(kklti_ln_xuv)*mu,2._kp
       stop
    end if

  end function kklti_ln_xuv

end module kkltisr

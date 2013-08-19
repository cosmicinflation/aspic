!slow-roll functions for the brane inflation potential
!
!V(phi) = M**4 [1 - x**(-p)]
!
!x = phi/mu
!mu = mu/Mp

module bisr
  use infprec, only : pi, kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public bi_norm_potential, bi_norm_deriv_potential, bi_norm_deriv_second_potential
  public bi_epsilon_one, bi_epsilon_two,bi_epsilon_three, bi_x_epstwounity
  public bi_x_epsoneunity, bi_efold_primitive, bi_x_trajectory
  public bi_ln_xstg, bi_ln_xuv
 
contains
!returns V/M**4
  function bi_norm_potential(x,p,mu)
    implicit none
    real(kp) :: bi_norm_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in) :: mu

    bi_norm_potential = 1._kp - x**(-p)
  end function bi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function bi_norm_deriv_potential(x,p,mu)
    implicit none
    real(kp) :: bi_norm_deriv_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

   bi_norm_deriv_potential = p*x**(-1._kp-p)

  end function bi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function bi_norm_deriv_second_potential(x,p,mu)
    implicit none
    real(kp) :: bi_norm_deriv_second_potential
    real(kp), intent(in) :: x,p
    real(kp), intent(in), optional :: mu

    bi_norm_deriv_second_potential = -p*(1._kp+p)*x**(-2._kp-p)

  end function bi_norm_deriv_second_potential

!epsilon1(x)
  function bi_epsilon_one(x,p,mu)    
    implicit none
    real(kp) :: bi_epsilon_one
    real(kp), intent(in) :: x,p,mu
    
    bi_epsilon_one = p**2/(2._kp*mu**2*x**2*(-1._kp+x**p)**2)
    
  end function bi_epsilon_one


!epsilon2(x)
  function bi_epsilon_two(x,p,mu)    
    implicit none
    real(kp) :: bi_epsilon_two
    real(kp), intent(in) :: x,p,mu
    
    bi_epsilon_two = (2._kp*p*(-1._kp+(1._kp+p)*x**p))/ &
                     (mu**2*x**2*(-1._kp+x**p)**2)
    
  end function bi_epsilon_two

!epsilon3(x)
  function bi_epsilon_three(x,p,mu)    
    implicit none
    real(kp) :: bi_epsilon_three
    real(kp), intent(in) :: x,p,mu
    
    bi_epsilon_three = (p*(2._kp+(1._kp+p)*x**p*(-4._kp+p+(2._kp+p)*x**p)))/ &
                       (mu**2*x**2*(-1._kp+x**p)**2*(-1._kp+(1._kp+p)*x**p))
    
  end function bi_epsilon_three



!this is integral[V(phi)/V'(phi) dphi]
  function bi_efold_primitive(x,p,mu)
    implicit none
    real(kp), intent(in) :: x,p,mu
    real(kp) :: bi_efold_primitive

    if (p.eq.0._kp) stop 'sf_nufunc: p=0 is singular'

    bi_efold_primitive = mu**2/p*(x**(p+2._kp)/(p+2._kp)-x**2/2._kp)

  end function bi_efold_primitive



!returns x for which epsilon1=1
  function bi_x_epsoneunity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: bi_x_epsoneunity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: biData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    biData%real1 = p
    biData%real2 = mu

    bi_x_epsoneunity = zbrent(find_bi_x_epsoneunity,mini,maxi,tolFind,biData)
   
  end function bi_x_epsoneunity
  
  function find_bi_x_epsoneunity(x,biData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: biData
    real(kp) :: find_bi_x_epsoneunity
    real(kp) :: p,mu
    
    p=biData%real1
    mu=biData%real2
!this is epsilon1(x)=1
    find_bi_x_epsoneunity = p-sqrt(2._kp)*mu*x*(x**p-1._kp)
    
  end function find_bi_x_epsoneunity
 


!returns x for which epsilon2=1
  function bi_x_epstwounity(p,mu)
    implicit none
    real(kp), intent(in) :: p,mu
    real(kp) :: bi_x_epstwounity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: biData

    mini = 1._kp + epsilon(1._kp)
    maxi = 1._kp/epsilon(1._kp)

    biData%real1 = p
    biData%real2 = mu

    bi_x_epstwounity = zbrent(find_bi_x_epstwounity,mini,maxi,tolFind,biData)
   
  end function bi_x_epstwounity
  
  function find_bi_x_epstwounity(x,biData)
    use infprec, only : transfert    
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: biData
    real(kp) :: find_bi_x_epstwounity
    real(kp) :: p,mu
    
    p=biData%real1
    mu=biData%real2
!this is epsilon2(x)=1
    find_bi_x_epstwounity = 2._kp*p*((1+p)*x**p-1._kp)- (mu*x*(x**p-1))**2
    
  end function find_bi_x_epstwounity
 



!returns x at bfold=-efolds before the end of inflation
  function bi_x_trajectory(bfold,xend,p,mu)
    implicit none
    real(kp), intent(in) :: bfold, p, mu, xend
    real(kp) :: bi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: biData

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    biData%real1 = p
    biData%real2 = mu
    biData%real3 = -bfold + bi_efold_primitive(xend,p,mu)
    
    bi_x_trajectory = zbrent(find_bi_x_trajectory,mini,maxi,tolFind,biData)
    
  end function bi_x_trajectory

  function find_bi_x_trajectory(x,biData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: biData
    real(kp) :: find_bi_x_trajectory
    real(kp) :: p,mu,NplusPrimEnd

    p=biData%real1
    mu = biData%real2
    NplusPrimEnd = biData%real3

    find_bi_x_trajectory = bi_efold_primitive(x,p,mu) - NplusPrimEnd
   
  end function find_bi_x_trajectory

!return xstg in terms of the universal string variable y,vbar, flux
!number calN et string parama alpha'
  function bi_ln_xstg(p,mu,y,vbar,calN,alpha)
    implicit none
    real(kp) :: bi_ln_xstg
    real(kp), intent(in) :: p,mu,y,vbar,calN,alpha

    if (p.ne.4._kp) stop 'bi_ln_xstg: p >< 4!'

    bi_ln_xstg = 0.25_kp*calN + y**(-0.25_kp)
    
  end function bi_ln_xstg

!return xUV in terms of the universal string variable y,vbar, flux
!number calN et string parama alpha'
  function bi_ln_xuv(p,mu,y,vbar,calN,alpha)
    implicit none
    real(kp) :: bi_ln_xuv
    real(kp), intent(in) :: p,mu,y,vbar,calN,alpha

    if (p.ne.4._kp) stop 'bi_ln_xuv: p >< 4!'

    bi_ln_xuv = -log(mu) + 0.75_kp*log(y) + 0.5_kp*log(vbar) &
         -log(pi) - 0.5_kp*log(2._kp*alpha*calN)

    if (bi_ln_xuv.gt.log(2._kp)-log(mu)) then
       write(*,*)'bi_ln_xuv: throat volume larger than bulk volume!'
       write(*,*)'phiUV= phiUVmax= ',exp(bi_ln_xuv)*mu,2._kp
       stop
    end if

  end function bi_ln_xuv


end module bisr

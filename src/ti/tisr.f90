!slow-roll functions for the Tip inflation potential
!
!V(phi) = M^4 [ 1 + cos(x) + alpha sin^2(x) ]
!
!x = phi/mu

module tisr
  use infprec, only : kp, tolkp,transfert
  use inftools, only : zbrent
  implicit none

  private

  public ti_norm_potential, ti_norm_deriv_potential, ti_norm_deriv_second_potential
  public ti_epsilon_one, ti_epsilon_two,ti_epsilon_three
  public ti_efold_primitive, ti_x_trajectory,ti_x_endinf

 
contains
!returns V/M**4
  function ti_norm_potential(x,alpha)
    implicit none
    real(kp) :: ti_norm_potential
    real(kp), intent(in) :: x,alpha

    ti_norm_potential = 1._kp+cos(x)+alpha*sin(x)**2
  end function ti_norm_potential


!returns the first derivative of the potential with respect to x=phi/mu, divided by M**4
  function ti_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ti_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   ti_norm_deriv_potential = sin(x)*(2._kp*alpha*cos(x)-1._kp)

  end function ti_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function ti_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ti_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    ti_norm_deriv_second_potential = 2._kp*alpha*cos(2._kp*x)-cos(x)

  end function ti_norm_deriv_second_potential

!epsilon1(x)
  function ti_epsilon_one(x,alpha,mu)    
    implicit none
    real(kp) :: ti_epsilon_one
    real(kp), intent(in) :: x,alpha,mu
    
    ti_epsilon_one = ((1._kp-2._kp*alpha*cos(x))**2*sin(x)**2)/ &
                     (2._kp*mu**2*(1._kp+cos(x)+alpha*sin(x)**2)**2)
    
  end function ti_epsilon_one


!epsilon2(x)
  function ti_epsilon_two(x,alpha,mu)    
    implicit none
    real(kp) :: ti_epsilon_two
    real(kp), intent(in) :: x,alpha,mu
    
    ti_epsilon_two = -((2._kp*cos(x/2._kp)**2*(-2._kp-alpha*(3._kp+4._kp*alpha) + & 
                     alpha*((6._kp+4._kp*alpha)*cos(x)+cos(2._kp*x))))/ &
                     (mu*(1._kp+cos(x)+alpha*sin(x)**2))**2)
    
  end function ti_epsilon_two

!epsilon3(x)
  function ti_epsilon_three(x,alpha,mu)    
    implicit none
    real(kp) :: ti_epsilon_three
    real(kp), intent(in) :: x,alpha,mu
    
    ti_epsilon_three = (-2._kp+(-2._kp-4._kp*alpha)/(1._kp+alpha-alpha*cos(x))**2 + &
                       (5._kp+3._kp*alpha)/(1._kp+alpha-alpha*cos(x))+(4._kp*(1._kp+alpha+3._kp*alpha**2) &
                       -2._kp*alpha*(7._kp+4._kp*alpha)*cos(x))/(-2._kp-alpha*(3._kp+4._kp*alpha)+ &
                       alpha*((6._kp+4._kp*alpha)*cos(x)+cos(2._kp*x)))+1._kp/(cos(x/2._kp)**2))/mu**4
    
  end function ti_epsilon_three

!returns x at the end of inflation defined as epsilon1=1
  function ti_x_endinf(alpha,mu)
    implicit none
    real(kp) :: ti_x_endinf
    real(kp), intent(in) :: alpha,mu
    complex(kp) :: Delta,smalldelta,Sigma,smallsigma,smallsigmaprime

    Delta=-864._kp*alpha**6*(1._kp+2._kp*alpha)**3*mu**2*(2._kp+mu**2)**2*((-1._kp+2._kp*alpha)**3 + & 
         2._kp*(1._kp+2._kp*alpha)*(-2._kp+(-10._kp+alpha)*alpha)*mu**2-4._kp*(1._kp+2._kp*alpha)**2*mu**4)
    smalldelta=8._kp*alpha**3*(2._kp*(-1._kp+2._kp*alpha)**3-3._kp*(1._kp+2._kp*alpha)*(5._kp+2._kp*alpha)* &
               (1._kp+4._kp*alpha)*mu**2-15._kp*(1._kp+alpha)*(1._kp+2._kp*alpha)**2*mu**4-2._kp* &
               (1._kp+2._kp*alpha)**3*mu**6)
    smallsigma=3._kp+4._kp*alpha-4._kp*alpha**2-2._kp*(mu+2._kp*alpha*mu)**2-8._kp/(2._kp+mu**2)
    Sigma=(2._kp+2._kp*alpha+2._kp*mu**2+alpha*mu**2)/(6._kp*alpha+3._kp*alpha*mu**2)
    smallsigmaprime=1/(2._kp*alpha**2*(2._kp+mu**2))

    ti_x_endinf =acos(real(Sigma+(1._kp+(0._kp,1._kp)*sqrt(3._kp))*smallsigma/(3._kp*2._kp**(2._kp/3._kp)* &
                 (smalldelta+sqrt(Delta))**(1._kp/3._kp))-(1._kp-(0._kp,1._kp)*sqrt(3._kp))*smallsigmaprime* &
                 (smalldelta+sqrt(Delta))**(1._kp/3._kp)/(6._kp*2._kp**(1._kp/3._kp)),kp))


    ti_x_endinf =min(ti_x_endinf,acos(-1._kp)*(1._kp-epsilon(1._kp))) !to avoid numerical NaN when ti_x_endinf is very close to pi

  end function ti_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function ti_efold_primitive(x,alpha,mu)
    implicit none
    real(kp), intent(in) :: x,alpha,mu
    real(kp) :: ti_efold_primitive
    
     if (alpha.eq.0.5_kp) stop 'ti_efold_primitive: alpha=1/2 is singular'

    ti_efold_primitive = mu**2*((2._kp*alpha+1._kp)/(2._kp*(1._kp-2._kp*alpha))* &
                         log(-2._kp*alpha*cos(x)+1._kp)+log(1-cos(x))/(2._kp*alpha-1._kp))

  end function ti_efold_primitive
 


!returns x at bfold=-efolds before the end of inflation
  function ti_x_trajectory(bfold,xend,alpha,mu)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,mu
    real(kp) :: ti_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: tiData

    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))

    tiData%real1 = alpha
    tiData%real2 = mu
    tiData%real3 = -bfold + ti_efold_primitive(xend,alpha,mu)
    
    ti_x_trajectory = zbrent(find_ti_x_trajectory,mini,maxi,tolFind,tiData)
    
  end function ti_x_trajectory

  function find_ti_x_trajectory(x,tiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: tiData
    real(kp) :: find_ti_x_trajectory
    real(kp) :: alpha,mu,NplusPrimEnd

    alpha=tiData%real1
    mu=tiData%real2
    NplusPrimEnd = tiData%real3

    find_ti_x_trajectory = ti_efold_primitive(x,alpha,mu) - NplusPrimEnd
   
  end function find_ti_x_trajectory



end module tisr

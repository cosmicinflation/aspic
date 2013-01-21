!slow-roll functions for the orientifold inflation potential
!
!V(phi) = M**4 x**4 [ alpha ln**2(x) - 1 ]
!
!x = phi/phi0

module oisr
  use infprec, only : kp, tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public oi_norm_potential, oi_norm_deriv_potential, oi_norm_deriv_second_potential
  public oi_epsilon_one, oi_epsilon_two,oi_epsilon_three
  public oi_efold_primitive, oi_x_trajectory,oi_x_endinf

 
contains
!returns V/M**4
  function oi_norm_potential(x,alpha,phi0)
    implicit none
    real(kp) :: oi_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: phi0

    oi_norm_potential = x**4*(alpha*log(x)**2-1._kp)

  end function oi_norm_potential


!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function oi_norm_deriv_potential(x,alpha,phi0)
    implicit none
    real(kp) :: oi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: phi0

   oi_norm_deriv_potential = 2._kp*x**3*(-2._kp+alpha*log(x)*(1._kp+2._kp*log(x)))

  end function oi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function oi_norm_deriv_second_potential(x,alpha,phi0)
    implicit none
    real(kp) :: oi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in), optional :: phi0

    oi_norm_deriv_second_potential = 2._kp*x**2*(-6._kp+alpha+alpha* &
         log(x)*(7._kp+6._kp*log(x)))

  end function oi_norm_deriv_second_potential

!epsilon1(x)
  function oi_epsilon_one(x,alpha,phi0)    
    implicit none
    real(kp) :: oi_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0
    
    oi_epsilon_one = ((2._kp*(-2._kp+alpha*log(x)*(1._kp+2._kp*log(x)))**2)/ &
         (phi0**2*x**2*(-1._kp+alpha*log(x)**2)**2))
    
  end function oi_epsilon_one


!epsilon2(x)
  function oi_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: oi_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0
    
    oi_epsilon_two = ((4._kp*(2._kp+alpha+alpha*log(x)*(-1._kp+ &
         log(x)*(-4._kp+alpha+alpha*log(x)*(1._kp+2._kp* &
         log(x))))))/(phi0**2*x**2*(-1._kp+alpha*log(x)**2)**2))
    
  end function oi_epsilon_two

!epsilon3(x)
  function oi_epsilon_three(x,alpha,phi0)    
    implicit none
    real(kp) :: oi_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0
    
    oi_epsilon_three = ((2._kp*(8._kp+6._kp*alpha-alpha*(8._kp+15._kp*alpha)* &
         log(x)+2._kp*alpha*(-16._kp-2._kp*alpha+3._kp*alpha**2)* &
         log(x)**2+8._kp*alpha**2*(3._kp+alpha)*log(x)**3+ & 
         2._kp*alpha**2*(24._kp-5._kp*alpha+alpha**2)*log(x)**4+ & 
         alpha**3*(-24._kp+7._kp*alpha)*log(x)**5+8._kp*(-4._kp+alpha)* &
         alpha**3*log(x)**6+8._kp*alpha**4*log(x)**7+8._kp*alpha**4* &
         log(x)**8))/(phi0**2*x**2*(-1._kp+alpha*log(x)**2)**2*(2._kp+ & 
         alpha-alpha*log(x)+(-4._kp+alpha)*alpha*log(x)**2+alpha**2* &
         log(x)**3+2._kp*alpha**2*log(x)**4)))
    
  end function oi_epsilon_three

!returns x at the end of inflation defined as epsilon1=1
  function oi_x_endinf(alpha,phi0)
    implicit none
    real(kp), intent(in) :: alpha,phi0
    real(kp) :: oi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: oiData
!Position where the potential vanishes and epsilon1 diverges    
    mini = exp(1._kp/sqrt(alpha))*(1._kp+epsilon(1._kp)) 
    maxi = mini/epsilon(1._kp)

    oiData%real1 = alpha
    oiData%real2 = phi0
    
    oi_x_endinf = zbrent(find_oi_x_endinf,mini,maxi,tolFind,oiData)


  end function oi_x_endinf

  function find_oi_x_endinf(x,oiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: oiData
    real(kp) :: find_oi_x_endinf
    real(kp) :: alpha,phi0

    alpha = oiData%real1
    phi0 = oiData%real2
    
    find_oi_x_endinf = oi_epsilon_one(x,alpha,phi0) - 1._kp
  
  end function find_oi_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function oi_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: oi_efold_primitive
    real(kp) :: xminus,xplus
    
     if (alpha.eq.0._kp) stop 'oi_efold_primitive: alpha=0'

!Lower position where the derivative of the potential vanishes
     xminus=exp(-0.25_kp-sqrt(1._kp/16._kp+1._kp/alpha))
!Upper position where the derivative of the potential vanishes
     xplus=exp(-0.25_kp+sqrt(1._kp/16._kp+1._kp/alpha))

    oi_efold_primitive = phi0**2*(x**2/8._kp &
         +(alpha*log(xplus)**2-1._kp)/(2._kp*sqrt(alpha**2+16._kp*alpha)) &
         *xplus**2* ei(2._kp*log(x/xplus)) -(alpha*log(xminus)**2-1._kp) &
         /(2._kp*sqrt(alpha**2+16._kp*alpha))*xminus**2* ei(2._kp*log(x/xminus)))

  end function oi_efold_primitive
 


!returns x at bfold=-efolds before the end of inflation
  function oi_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold,xend,alpha,phi0
    real(kp) :: oi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: oiData

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    oiData%real1 = alpha
    oiData%real2 = phi0
    oiData%real3 = -bfold + oi_efold_primitive(xend,alpha,phi0)
    
    oi_x_trajectory = zbrent(find_oi_x_trajectory,mini,maxi,tolFind,oiData)
    
  end function oi_x_trajectory

  function find_oi_x_trajectory(x,oiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: oiData
    real(kp) :: find_oi_x_trajectory
    real(kp) :: alpha,phi0,NplusPrimEnd

    alpha=oiData%real1
    phi0=oiData%real2
    NplusPrimEnd = oiData%real3

    find_oi_x_trajectory = oi_efold_primitive(x,alpha,phi0) - NplusPrimEnd
   
  end function find_oi_x_trajectory



end module oisr

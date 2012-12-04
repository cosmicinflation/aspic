!slow-roll functions for the constant ns B inflation potential
!
!V(phi) = M**4 ((3-alpha**2) tan**2(alpha x / sqrt(2) ) -2 )
!
!x = phi/Mp

module cnbisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  use cosmopar, only : pi
  implicit none

  private

  public  cnbi_norm_potential, cnbi_epsilon_one, cnbi_epsilon_two, cnbi_epsilon_three
  public  cnbi_x_endinf, cnbi_efold_primitive, cnbi_x_trajectory, cnbi_x_max
  public  cnbi_norm_deriv_potential, cnbi_norm_deriv_second_potential
 
contains
!returns V/M**4
  function cnbi_norm_potential(x,alpha)
    implicit none
    real(kp) :: cnbi_norm_potential
    real(kp), intent(in) :: x,alpha

    cnbi_norm_potential = (3._kp-alpha**2)*tan(alpha/sqrt(2._kp)*x)**2-3._kp

  end function cnbi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function cnbi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: cnbi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   cnbi_norm_deriv_potential = -sqrt(2._kp)*alpha*(-3._kp+alpha**2)*1._kp/ &
                               cos((alpha*x)/sqrt(2._kp))**2*tan((alpha*x)/sqrt(2._kp))

  end function cnbi_norm_deriv_potential



!returns the *1._kp/cosond derivative of the potential with respect to x, divided by M**4
  function cnbi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: cnbi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    cnbi_norm_deriv_second_potential = alpha**2*(-3._kp+alpha**2)*(-2._kp+cos(sqrt(2._kp)* &
                                       alpha*x))*1._kp/cos((alpha*x)/sqrt(2._kp))**4

  end function cnbi_norm_deriv_second_potential



!epsilon_one(x)
  function cnbi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_one = (4._kp*alpha**2*(-3._kp+alpha**2)**2*tan((alpha*x)/sqrt(2._kp))**2)/ &
                       (alpha**2-(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x))**2
    
  end function cnbi_epsilon_one


!epsilon_two(x)
  function cnbi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_two =-(alpha**2*(-3._kp+alpha**2)*(6._kp+alpha**2-2._kp*(-6._kp+alpha**2)* &
                      cos(sqrt(2._kp)*alpha*x)+(-6._kp+alpha**2)*cos(2._kp*sqrt(2._kp)*alpha*x))* &
                      1._kp/cos((alpha*x)/sqrt(2._kp))**6)/(2._kp*(3._kp+(-3._kp+alpha**2)* &
                      tan((alpha*x)/sqrt(2._kp))**2)**2)

  end function cnbi_epsilon_two


!epsilon_three(x)
  function cnbi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_three = (2._kp*alpha**2*(-3._kp+alpha**2)*(-6._kp*(72._kp-14._kp*alpha**2+alpha**4)+ &
                         (-6._kp+alpha**2)*(78._kp+7._kp*alpha**2)*cos(sqrt(2._kp)*alpha*x)- &
                         2._kp*(72._kp-18._kp*alpha**2+alpha**4)*cos(2._kp*sqrt(2._kp)*alpha*x)+ &
                         (-6._kp+alpha**2)**2*cos(3._kp*sqrt(2._kp)*alpha*x))*tan((alpha*x)/ &
                         sqrt(2._kp))**2)/((alpha**2-(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x))**2* &
                         (6._kp+alpha**2-2._kp*(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x)+ &
                         (-6._kp+alpha**2)*cos(2._kp*sqrt(2._kp)*alpha*x)))

  end function cnbi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function cnbi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnbiData

    mini = sqrt(2._kp)/alpha*atan(sqrt(3._kp/(3._kp-alpha**2)))*(1._kp+epsilon(1._kp)) !Position where the potential vanishes and eps1 diverges
    maxi = acos((alpha**2-6._kp+sqrt(alpha**4-36._kp*alpha**2+180._kp))/(2._kp*(alpha**2-6._kp)))/(alpha*sqrt(2._kp)) !Position where eps2 vanishes and eps1 is minimum

    cnbiData%real1 = alpha

    cnbi_x_endinf = zbrent(find_cnbiendinf,mini,maxi,tolFind,cnbiData)
   
  end function cnbi_x_endinf


!returns the maximum position x for the slow roll to be valid defined as epsilon1=1
  function cnbi_x_max(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnbiData

    mini = acos((alpha**2-6._kp+sqrt(alpha**4-36._kp*alpha**2+180._kp))/(2._kp*(alpha**2-6._kp)))/(alpha*sqrt(2._kp)) !Position where eps2 vanishes and eps1 is minimum
    maxi = pi/(sqrt(2._kp)*alpha) !Position where the potential diverges

    cnbiData%real1 = alpha

    cnbi_x_max = zbrent(find_cnbiendinf,mini,maxi,tolFind,cnbiData)
   
  end function cnbi_x_max


  function find_cnbiendinf(x,cnbiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: cnbiData
    real(kp) :: find_cnbiendinf
    real(kp) :: alpha
    
    alpha = cnbiData%real1
    
    find_cnbiendinf = cnbi_epsilon_one(x,alpha) - 1._kp
    
  end function find_cnbiendinf


!this is integral(V(phi)/V'(phi) dphi)
  function cnbi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: cnbi_efold_primitive

    if (alpha.eq.0._kp) stop 'cnbi_efold_primitive: alpha=0!'

    cnbi_efold_primitive = -1._kp/(alpha**2*(3._kp-alpha**2))* &
                           (3._kp*log(sin(alpha*x/sqrt(2._kp)))- &
                            (6._kp-alpha**2)/2._kp*sin(alpha*x/sqrt(2._kp))**2)

  end function cnbi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function cnbi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: cnbi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnbiData

  
    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = cnbi_x_max(alpha)*(1._kp-epsilon(1._kp)) !Position where slow roll stops being valid
  


    cnbiData%real1 = alpha
    cnbiData%real2 = -bfold + cnbi_efold_primitive(xend,alpha)
    
    cnbi_x_trajectory = zbrent(find_cnbitraj,mini,maxi,tolFind,cnbiData)
       
  end function cnbi_x_trajectory

  function find_cnbitraj(x,cnbiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cnbiData
    real(kp) :: find_cnbitraj
    real(kp) :: alpha,NplusNuend

    alpha = cnbiData%real1
    NplusNuend = cnbiData%real2

    find_cnbitraj = cnbi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_cnbitraj


  
end module cnbisr

!slow-roll functions for the constant ns A inflation potential
!
!V(phi) = M**4 [ 3 - (3+alpha**2) tanh**2( alpha/sqrt(2) x) )
!
!x = phi/Mp

module cnaisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert, atanh
  use inftools, only : zbrent
  implicit none

  private

  public  cnai_norm_potential, cnai_epsilon_one, cnai_epsilon_two, cnai_epsilon_three
  public  cnai_x_endinf, cnai_efold_primitive, cnai_x_trajectory
  public  cnai_norm_deriv_potential, cnai_norm_deriv_second_potential
 
contains
!returns V/M**4
  function cnai_norm_potential(x,alpha)
    implicit none
    real(kp) :: cnai_norm_potential
    real(kp), intent(in) :: x,alpha

    cnai_norm_potential = 3._kp-(3._kp+alpha**2)*tanh(alpha/sqrt(2._kp)*x)**2

  end function cnai_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function cnai_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: cnai_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   cnai_norm_deriv_potential = -sqrt(2._kp)*alpha*(3._kp+alpha**2)* &
                               1._kp/cosh((alpha*x)/sqrt(2._kp))**2* &
                               tanh((alpha*x)/sqrt(2._kp))

  end function cnai_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function cnai_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: cnai_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    cnai_norm_deriv_second_potential = alpha**2*(3._kp+alpha**2)* &
                                       (-2._kp+cosh(sqrt(2._kp)*alpha*x))* & 
                                       1._kp/cosh((alpha*x)/sqrt(2._kp))**4

  end function cnai_norm_deriv_second_potential



!epsilon_one(x)
  function cnai_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: cnai_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    cnai_epsilon_one = (4._kp*alpha**2*(3._kp+alpha**2)**2*tanh((alpha*x)/ &
                       sqrt(2._kp))**2)/(6._kp+alpha**2-alpha**2* &
                       cosh(sqrt(2._kp)*alpha*x))**2
    
  end function cnai_epsilon_one


!epsilon_two(x)
  function cnai_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: cnai_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    cnai_epsilon_two =(2._kp*alpha**2*(3._kp+alpha**2)*(12._kp+alpha**2+ &
                      alpha**2*(-2._kp*cosh(sqrt(2._kp)*alpha*x)+ &
                      cosh(2._kp*sqrt(2._kp)*alpha*x)))* &
                      1._kp/cosh((alpha*x)/sqrt(2._kp))**2)/ &
                      (6._kp+alpha**2-alpha**2*cosh(sqrt(2._kp)*alpha*x))**2

  end function cnai_epsilon_two


!epsilon_three(x)
  function cnai_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: cnai_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    cnai_epsilon_three = (2._kp*alpha**2*(3._kp+alpha**2)* &
                         (-6._kp*(24._kp-2._kp*alpha**2+alpha**4)+ &
                         alpha**2*((120._kp+7._kp*alpha**2)*cosh(sqrt(2._kp)*alpha*x)- &
                         2._kp*(-6._kp+alpha**2)*cosh(2._kp*sqrt(2._kp)*alpha*x)+ &
                         alpha**2*cosh(3._kp*sqrt(2._kp)*alpha*x)))* &
                         tanh((alpha*x)/sqrt(2._kp))**2)/((6._kp+alpha**2-alpha**2* &
                         cosh(sqrt(2._kp)*alpha*x))**2*(12._kp+alpha**2+alpha**2* &
                         (-2._kp*cosh(sqrt(2._kp)*alpha*x)+cosh(2._kp*sqrt(2._kp)*alpha*x))))


  end function cnai_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function cnai_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnai_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnaiData

    mini = epsilon(1._kp)
    maxi = sqrt(2._kp)/alpha*atanh(sqrt(3._kp/(3._kp+alpha**2)))*(1._kp-epsilon(1._kp)) !Position where the potential vanishes

    cnaiData%real1 = alpha

    cnai_x_endinf = zbrent(find_cnaiendinf,mini,maxi,tolFind,cnaiData)
   
  end function cnai_x_endinf


  function find_cnaiendinf(x,cnaiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: cnaiData
    real(kp) :: find_cnaiendinf
    real(kp) :: alpha
    
    alpha = cnaiData%real1
    
    find_cnaiendinf = cnai_epsilon_one(x,alpha) - 1._kp
    
  end function find_cnaiendinf


!this is integral[V(phi)/V'(phi) dphi)
  function cnai_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: cnai_efold_primitive

    if (alpha.eq.0._kp) stop 'cnai_efold_primitive: alpha=0!'

    cnai_efold_primitive = -1._kp/(alpha**2*(3._kp+alpha**2))* &
                           (3._kp*log(sinh(alpha*x/sqrt(2._kp)))- &
                            alpha**2/2._kp*sinh(alpha*x/sqrt(2._kp))**2)

  end function cnai_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function cnai_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: cnai_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnaiData

  
    mini = epsilon(1._kp)
    maxi = xEnd*(1._kp-epsilon(1._kp))
  


    cnaiData%real1 = alpha
    cnaiData%real2 = -bfold + cnai_efold_primitive(xend,alpha)
    
    cnai_x_trajectory = zbrent(find_cnaitraj,mini,maxi,tolFind,cnaiData)
       
  end function cnai_x_trajectory

  function find_cnaitraj(x,cnaiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cnaiData
    real(kp) :: find_cnaitraj
    real(kp) :: alpha,NplusNuend

    alpha = cnaiData%real1
    NplusNuend = cnaiData%real2

    find_cnaitraj = cnai_efold_primitive(x,alpha) - NplusNuend
   
  end function find_cnaitraj


  
end module cnaisr

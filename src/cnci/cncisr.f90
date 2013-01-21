!slow-roll functions for the constant ns C inflation potential
!
!V(phi) = M**4 ( (3+alpha**2) 1._kp/tanh( alpha/sqrt(2) x )**2 - 3 )
!
!x = phi/Mp

module cncisr
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  implicit none

  private

  public cnci_norm_potential, cnci_norm_deriv_potential, cnci_norm_deriv_second_potential
  public cnci_epsilon_one, cnci_epsilon_two, cnci_epsilon_three
  public cnci_efold_primitive, cnci_x_trajectory
  public cnci_xendmin, cnci_x_epsoneunity

 
contains
!returns V/M**4
  function cnci_norm_potential(x,alpha)
    implicit none
    real(kp) :: cnci_norm_potential
    real(kp), intent(in) :: x,alpha

    cnci_norm_potential = (3._kp+alpha**2)*1._kp/tanh(alpha/sqrt(2._kp)*x)**2-3._kp
  end function cnci_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function cnci_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: cnci_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   cnci_norm_deriv_potential = -sqrt(2._kp)*alpha*(3._kp+alpha**2)*1._kp/tanh((alpha*x) &
        /sqrt(2._kp))*1._kp/sinh((alpha*x)/sqrt(2._kp))**2

  end function cnci_norm_deriv_potential


!returns the second derivative of the potential with respect to x, divided by M**4
  function cnci_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: cnci_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    cnci_norm_deriv_second_potential = alpha**2*(3._kp+alpha**2) &
         *(2._kp+cosh(sqrt(2._kp)*alpha*x)) * 1._kp/sinh((alpha*x)/sqrt(2._kp))**4

  end function cnci_norm_deriv_second_potential

!epsilon1(x)
  function cnci_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: cnci_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    cnci_epsilon_one = (4._kp*alpha**2*(3._kp+alpha**2)**2*1._kp/tanh((alpha*x) &
         /sqrt(2._kp))**2)/(6._kp+alpha**2+alpha**2*cosh(sqrt(2._kp)*alpha*x))**2
    
  end function cnci_epsilon_one


!epsilon2(x)
  function cnci_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: cnci_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    cnci_epsilon_two = -(2._kp*alpha**2*(3._kp+alpha**2)*(12._kp+alpha**2+alpha**2* &
         (2._kp*cosh(sqrt(2._kp)*alpha*x)+cosh(2.*sqrt(2._kp)*alpha*x)))* &
         1._kp/sinh((alpha*x)/sqrt(2._kp))**2)/(6._kp+alpha**2+alpha**2* &
         cosh(sqrt(2._kp)*alpha*x))**2
    
  end function cnci_epsilon_two

!epsilon3(x)
  function cnci_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: cnci_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    cnci_epsilon_three = -(2._kp*alpha**2*(3._kp+alpha**2)*(6._kp*(24._kp-2._kp* &
         alpha**2+alpha**4)+alpha**2*((120._kp+7._kp*alpha**2)* &
         cosh(sqrt(2._kp)*alpha*x)+2._kp*(-6._kp+alpha**2)* &
         cosh(2._kp*sqrt(2._kp)*alpha*x)+alpha**2*cosh(3.*sqrt(2._kp)* &
         alpha*x)))*1._kp/tanh((alpha*x)/sqrt(2._kp))**2)/((6._kp+alpha**2+ &
         alpha**2*cosh(sqrt(2._kp)*alpha*x))**2*(12._kp+alpha**2+alpha**2* &
         (2._kp*cosh(sqrt(2._kp)*alpha*x)+cosh(2._kp*sqrt(2._kp)*alpha*x))))
    
  end function cnci_epsilon_three



!returns the minimal value for xend such that there are efold number
!of inflation from xiniMin=xepsone
  function cnci_xendmin(efold,alpha)
    implicit none
    real(kp), intent(in) :: efold,alpha
    real(kp) :: cnci_xendmin, xiniMin
    
    xiniMin = cnci_x_epsoneunity(alpha)       

    cnci_xendmin = cnci_x_trajectory(efold,xiniMin,alpha)
       
  end function cnci_xendmin


 
! Returns the value of x such that epsilon_one=1
  function cnci_x_epsoneunity(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnci_x_epsoneunity
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnciData

    mini = epsilon(1._kp)
    maxi = (log(12._kp/alpha)+log(4._kp*alpha))/(sqrt(2._kp*alpha))*10._kp**(3.)

    cnciData%real1 = alpha

    cnci_x_epsoneunity=zbrent(find_cnci_x_epsoneunity,mini,maxi,tolFind,cnciData)


  end function cnci_x_epsoneunity

  function find_cnci_x_epsoneunity(x,cnciData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cnciData
    real(kp) :: find_cnci_x_epsoneunity
    real(kp) :: alpha

    alpha= cnciData%real1

    find_cnci_x_epsoneunity = cnci_epsilon_one(x,alpha) - 1._kp
   
  end function find_cnci_x_epsoneunity


!this is integral(V(phi)/V'(phi) dphi)
  function cnci_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: cnci_efold_primitive
    
     if (alpha.eq.0._kp) stop 'cnci_efold_primitive: alpha=0 is singular'

    cnci_efold_primitive = -1._kp/(alpha**2*(3._kp+alpha**2)) &
         *(3._kp*log(cosh(alpha*x/sqrt(2._kp))) &
         + alpha**2/2._kp*cosh(alpha*x/sqrt(2._kp))**2)

  end function cnci_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function cnci_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold,alpha,xend
    real(kp) :: cnci_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnciData

    if (bfold .lt. 0._kp) then
    maxi = xend*(1._kp-epsilon(1._kp))
    mini = cnci_x_epsoneunity(alpha)
    else !Used to obtain xend_min
    mini=xend*(1._kp+epsilon(1._kp))
    maxi=10._kp**(6._kp)*mini
    endif

    cnciData%real1 = alpha
    cnciData%real2 = -bfold + cnci_efold_primitive(xend,alpha)
    
    cnci_x_trajectory = zbrent(find_cncitraj,mini,maxi,tolFind,cnciData)
       
  end function cnci_x_trajectory

  function find_cncitraj(x,cnciData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cnciData
    real(kp) :: find_cncitraj
    real(kp) :: alpha,NplusNuend

    alpha= cnciData%real1
    NplusNuend = cnciData%real2

    find_cncitraj = cnci_efold_primitive(x,alpha) - NplusNuend
   
  end function find_cncitraj

 

end module cncisr

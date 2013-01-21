!slow-roll functions for the constant ns B inflation potential
!
!V(phi) = M**4 ((3-alpha**2) tan**2(alpha x / sqrt(2) ) -3 )
!
!x = phi/Mp

module cnbisr
  use infprec, only : kp,tolkp,transfert,pi
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public cnbi_norm_potential, cnbi_epsilon_one, cnbi_epsilon_two, cnbi_epsilon_three
  public cnbi_x_endinf, cnbi_efold_primitive, cnbi_x_trajectory
  public cnbi_epsilon_one_min, cnbi_x_epsonemin
  public cnbi_x_potzero, cnbi_x_potinfty, cnbi_x_epsoneunity
  public cnbi_norm_deriv_potential, cnbi_norm_deriv_second_potential
 
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

    cnbi_norm_deriv_second_potential = alpha**2*(-3._kp+alpha**2) &
         *(-2._kp+cos(sqrt(2._kp)* &
         alpha*x))*1._kp/cos((alpha*x)/sqrt(2._kp))**4

  end function cnbi_norm_deriv_second_potential



!epsilon_one(x)
  function cnbi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_one = (4._kp*alpha**2*(-3._kp+alpha**2)**2*tan((alpha*x) &
         /sqrt(2._kp))**2)/ &
         (alpha**2-(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x))**2
    
  end function cnbi_epsilon_one


!epsilon_two(x)
  function cnbi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_two =-(alpha**2*(-3._kp+alpha**2) &
         *(6._kp+alpha**2-2._kp*(-6._kp+alpha**2)* &
         cos(sqrt(2._kp)*alpha*x)+(-6._kp+alpha**2) &
         *cos(2._kp*sqrt(2._kp)*alpha*x))* &
         1._kp/cos((alpha*x)/sqrt(2._kp))**6)/(2._kp*(3._kp+(-3._kp+alpha**2)* &
         tan((alpha*x)/sqrt(2._kp))**2)**2)

  end function cnbi_epsilon_two


!epsilon_three(x)
  function cnbi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: cnbi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    cnbi_epsilon_three = (2._kp*alpha**2*(-3._kp+alpha**2) &
         *(-6._kp*(72._kp-14._kp*alpha**2+alpha**4)+ &
         (-6._kp+alpha**2)*(78._kp+7._kp*alpha**2)*cos(sqrt(2._kp)*alpha*x)- &
         2._kp*(72._kp-18._kp*alpha**2+alpha**4)*cos(2._kp*sqrt(2._kp)*alpha*x)+ &
         (-6._kp+alpha**2)**2*cos(3._kp*sqrt(2._kp)*alpha*x))*tan((alpha*x)/ &
         sqrt(2._kp))**2)/((alpha**2-(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x))**2* &
         (6._kp+alpha**2-2._kp*(-6._kp+alpha**2)*cos(sqrt(2._kp)*alpha*x)+ &
         (-6._kp+alpha**2)*cos(2._kp*sqrt(2._kp)*alpha*x)))

  end function cnbi_epsilon_three



!returns the first field value at which the potential vanishes
  function cnbi_x_potzero(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_potzero
    
    cnbi_x_potzero = sqrt(2._kp)/alpha &
         *atan(sqrt(3._kp/(3._kp-alpha**2)))

  end function cnbi_x_potzero



!returns the first field value at which the potential diverges
  function cnbi_x_potinfty(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_potinfty

    cnbi_x_potinfty = pi/(sqrt(2._kp)*alpha)

  end function cnbi_x_potinfty



!returns the field value at which eps2=0 et eps1 is minimal
  function cnbi_x_epsonemin(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_epsonemin

    cnbi_x_epsonemin = acos((alpha**2-6._kp+sqrt(alpha**4-36._kp*alpha**2+180._kp)) &
         /(2._kp*(alpha**2-6._kp)))/(alpha*sqrt(2._kp)) 

  end function cnbi_x_epsonemin



!returns the actual value of eps1 when it is minimal
  function cnbi_epsilon_one_min(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_epsilon_one_min
    real(kp) :: xeps1min

    xeps1min = cnbi_x_epsonemin(alpha)

    cnbi_epsilon_one_min = cnbi_epsilon_one(xeps1min,alpha)

  end function cnbi_epsilon_one_min




!returns a vector containing the two roots of eps1=1
  function cnbi_x_epsoneunity(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp), dimension(2) :: cnbi_x_epsoneunity

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: cnbiData

    real(kp) :: xZero, xInfty, xeps1Min

    if (cnbi_epsilon_one_min(alpha).gt.1._kp) then
       stop 'cnbi_x_epsoneunity: alpha>0.2975; eps1=1 has no solution!'
    endif

    xeps1Min = cnbi_x_epsonemin(alpha)

    if (cnbi_epsilon_one_min(alpha).eq.1._kp) then
       cnbi_x_epsoneunity(1) = xeps1Min
       cnbi_x_epsoneunity(2) = xeps1Min
       return
    endif

   
    xZero = cnbi_x_potzero(alpha)
    xInfty = cnbi_x_potinfty(alpha)
    

!first root
    mini = xZero + epsilon(1._kp)
    maxi = xeps1Min

    cnbiData%real1 = alpha

    cnbi_x_epsoneunity(1) = zbrent(find_cnbi_x_epsoneunity,mini,maxi,tolFind,cnbiData)

!second root
    mini = xeps1Min + epsilon(1._kp)
    maxi = xInfty

    cnbi_x_epsoneunity(2) = zbrent(find_cnbi_x_epsoneunity,mini,maxi,tolFind,cnbiData)


  end function cnbi_x_epsoneunity


  function find_cnbi_x_epsoneunity(x,cnbiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: cnbiData
    real(kp) :: find_cnbi_x_epsoneunity
    real(kp) :: alpha
    
    alpha = cnbiData%real1
    
    find_cnbi_x_epsoneunity = cnbi_epsilon_one(x,alpha) - 1._kp
    
  end function find_cnbi_x_epsoneunity



  function cnbi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: cnbi_x_endinf

    real(kp), dimension(2) :: xeps1
    
    xeps1 = cnbi_x_epsoneunity(alpha)
    
    cnbi_x_endinf = xeps1(1)

  end function cnbi_x_endinf




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

    real(kp) :: efoldMax
    real(kp), dimension(2) :: xeps1

    xeps1 = cnbi_x_epsoneunity(alpha)

    efoldMax = -cnbi_efold_primitive(xEnd,alpha) &
         + cnbi_efold_primitive(xeps1(2),alpha)

    if (-bfold.gt.efoldMax) then
       write(*,*)'cnbi_x_trajectory: not enough efolds!'
       write(*,*)'efold requested=   efold maxi= ',-bfold,efoldMax
       stop
    endif


    mini = xEnd+epsilon(1._kp)
    maxi = xeps1(2) - epsilon(1._kp)

    cnbiData%real1 = alpha
    cnbiData%real2 = -bfold + cnbi_efold_primitive(xend,alpha)
    
    cnbi_x_trajectory = zbrent(find_cnbi_x_trajectory,mini,maxi,tolFind,cnbiData)
       
  end function cnbi_x_trajectory

  function find_cnbi_x_trajectory(x,cnbiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: cnbiData
    real(kp) :: find_cnbi_x_trajectory
    real(kp) :: alpha,NplusNuend

    alpha = cnbiData%real1
    NplusNuend = cnbiData%real2

    find_cnbi_x_trajectory = cnbi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_cnbi_x_trajectory


  
end module cnbisr

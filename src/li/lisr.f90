!slow-roll functions for the loop inflation potential
!
!V(phi) = M^4 [1 + alpha ln(x) ]
!
!x = phi/Mp

module lisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public  li_norm_potential, li_epsilon_one, li_epsilon_two, li_epsilon_three
  public  li_x_endinf, li_efold_primitive, li_x_trajectory, li_alphamin
  public  li_norm_deriv_potential, li_norm_deriv_second_potential, li_x_epsoneunity
 
!under which eps1=1 has no solution
  real(kp), parameter :: alphaMin = 2._kp/(log(2._kp)-2)

!for 0 < alpha < alphaNumAccPos, exp(1/alpha) > Infinity while for
! alpha < alphaNumAccNeg < 0 xend^2 ln(xend) is too large to extract
! bfoldstar due to machine accuracy. The solution is a lambert
! function of the Bound number that we bound by one step recurrence
  real(kp), parameter :: alphaNumAccPos = 1._kp/log(huge(1._kp))
  real(kp), parameter :: bigBound = 1._kp/epsilon(1._kp) 
  real(kp), parameter :: alphaNumAccNeg = -2._kp/log(2._kp/log(bigBound)*bigBound)

contains

  subroutine li_check_accuracy(alpha)
    implicit none
    real(kp), intent(in) :: alpha

    if ((alpha.lt.alphaNumAccPos).and.(alpha.gt.alphaNumAccNeg)) then
       write(*,*)'li_check_accuracy:'
       write(*,*)'alpha= ',alpha
       write(*,*)'out of precision interval= ',alphaNumAccNeg,alphaNumAccPos
       stop
    endif

  end subroutine li_check_accuracy

!returns V/M^4
  function li_norm_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_potential
    real(kp), intent(in) :: x,alpha

    li_norm_potential = 1._kp+alpha*log(x)

  end function li_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function li_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   li_norm_deriv_potential = alpha/x

  end function li_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function li_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: li_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    li_norm_deriv_second_potential = -alpha/(x**2)

  end function li_norm_deriv_second_potential



!epsilon_one(x)
  function li_epsilon_one(x,alpha)
    implicit none
    real(kp) :: li_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_one = alpha**2/(2._kp*x**2) &
         /(1._kp+alpha*log(x))**2
    
  end function li_epsilon_one


!epsilon_two(x)
  function li_epsilon_two(x,alpha)
    implicit none
    real(kp) :: li_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_two = 2._kp*alpha/(x**2)*(1._kp+alpha+alpha*log(x)) &
         /(1._kp+alpha*log(x))**2
    
  end function li_epsilon_two


!epsilon_three(x)
  function li_epsilon_three(x,alpha)
    implicit none
    real(kp) :: li_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    li_epsilon_three = 2._kp*alpha/(x**2) &
         /(1._kp+alpha+alpha*log(x)) &
         /(1._kp+alpha*log(x))**2 &
         *(1._kp+3._kp*alpha/2._kp+alpha**2+ &
         (2._kp*alpha+3._kp*alpha**2/2._kp)*log(x) &
         +alpha**2*(log(x))**2)
    
  end function li_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function li_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: li_x_endinf
    real(kp), dimension(2) :: xepsones
    
    xepsones = li_x_epsoneunity(alpha)

    if (alpha .gt. 0._kp) then
       li_x_endinf = xepsones(1)
    else 
       li_x_endinf = xepsones(2)
    endif

  end function li_x_endinf


!returns the the roots of epsilon_1=1
  function li_x_epsoneunity(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp), dimension(2) :: li_x_epsoneunity

    call li_check_accuracy(alpha)
    
    if (alpha .gt. 0._kp) then
       li_x_epsoneunity = 1._kp/sqrt(2._kp) &
            /lambert(exp(1._kp/alpha)/(sqrt(2._kp)),0)
    elseif (alpha.lt.0._kp) then
       
       li_x_epsoneunity(1) = -1._kp/sqrt(2._kp) &
            /lambert(-exp(1._kp/alpha)/(sqrt(2._kp)),-1)

       li_x_epsoneunity(2) = -1._kp/sqrt(2._kp) &
            /lambert(-exp(1._kp/alpha)/(sqrt(2._kp)),0)

    else
       stop 'li_x_epsoneunity: ill defined alpha'
    endif
   
  end function li_x_epsoneunity


 
!this is integral[V(phi)/V'(phi) dphi]
  function li_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: li_efold_primitive

    call li_check_accuracy(alpha)

    li_efold_primitive = (-1._kp/4._kp+1._kp/(2._kp*alpha))*x**2 &
         +1._kp/2._kp*x**2*log(x)

  end function li_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function li_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: li_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: liData

    real(kp), dimension(2) :: xepsones

    xepsones = li_x_epsoneunity(alpha)
    
    if (alpha .gt. 0._kp) then

       mini = xEnd
       maxi = huge(1._kp)

    elseif (alpha.lt.0._kp) then
       
       mini = xepsones(1) + epsilon(1._kp)
       maxi = xEnd

    else
       stop 'li_x_trajectory: alpha ill defined!'
    endif

    liData%real1 = alpha
    liData%real2 = -bfold + li_efold_primitive(xend,alpha)
    
    li_x_trajectory = zbrent(find_li_x_trajectory,mini,maxi,tolFind,liData)
       
  end function li_x_trajectory

  function find_li_x_trajectory(x,liData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: liData
    real(kp) :: find_li_x_trajectory
    real(kp) :: alpha,mu,NplusNuend

    alpha = liData%real1
    NplusNuend = liData%real2

    find_li_x_trajectory = li_efold_primitive(x,alpha) - NplusNuend
   
  end function find_li_x_trajectory



!return the minimum value of alpha to get efold number of inflation
  function li_alphamin(efold)
    implicit none
    real(kp) :: li_alphamin
    real(kp), intent(in) :: efold
    real(kp), save :: efoldSave = 0._kp
    real(kp), save :: alphaSave = alphaMin
    type(transfert) :: liData
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi
    

    if (efold.eq.efoldSave) then
       li_alphamin = alphaSave
       return
    endif

    mini = 2._kp/(log(2._kp)-2._kp)
    maxi = alphaNumAccNeg
    
    liData%real1 = efold

    li_alphamin = zbrent(find_li_alphamin,mini,maxi,tolFind,liData)

    efoldSave = efold
    alphaSave = li_alphamin

  end function li_alphamin


  function find_li_alphamin(x,liData)
    implicit none
    real(kp), intent(in) :: x
    real(kp) :: find_li_alphamin
    type(transfert), optional, intent(inout) :: liData
    real(kp) :: efold
    real(kp), dimension(2) :: xepsones

    efold = liData%real1

    xepsones = li_x_epsoneunity(x)
   
    find_li_alphamin = efold - li_efold_primitive(xepsones(1),x) &
         +  li_efold_primitive(xepsones(2),x)

  end function find_li_alphamin

  
end module lisr

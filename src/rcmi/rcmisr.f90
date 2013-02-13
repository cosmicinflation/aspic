!slow-roll functions for the radiatively corrected massive inflation potential
!
!V(phi) = M^4 x^2 [ 1 - 2 alpha x^2 ln(x) ]
!
!x = phi/Mp

module rcmisr
  use infprec, only : kp, tolkp, transfert
  use specialinf, only : lambert
  use inftools, only : zbrent, easydverk
  implicit none

  private

  public rcmi_norm_potential, rcmi_epsilon_one, rcmi_epsilon_two
  public rcmi_epsilon_three
  public rcmi_x_endinf, rcmi_efold_primitive, rcmi_x_trajectory
  public rcmi_norm_deriv_potential, rcmi_norm_deriv_second_potential
  public rcmi_x_potmax
 
contains
!returns V/M^4
  function rcmi_norm_potential(x,alpha)
    implicit none
    real(kp) :: rcmi_norm_potential
    real(kp), intent(in) :: x,alpha

    rcmi_norm_potential = x**2*(1._kp-2._kp*alpha*x**2*log(x))

  end function rcmi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M^4
  function rcmi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: rcmi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha


   rcmi_norm_deriv_potential = 2._kp*x-8._kp*alpha*x**3*log(x) &
         -2._kp*alpha*x**3

  end function rcmi_norm_deriv_potential


!returns the second derivative of the potential with respect to x, divided by M^4
  function rcmi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: rcmi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    rcmi_norm_deriv_second_potential = 2._kp-14._kp*alpha*x**2 &
         -24._kp*alpha*x**2*log(x)

  end function rcmi_norm_deriv_second_potential


!epsilon_one(x)
  function rcmi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_one = 2/(x**2) &
         *(1._kp-alpha*x**2-4._kp*alpha*x**2*log(x))**2 &
         /(1._kp-2._kp*alpha*x**2*log(x))**2
    
  end function rcmi_epsilon_one


!epsilon_two(x)
  function rcmi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_two = 4._kp/(x**2) &
         /(1._kp-2._kp*alpha*x**2*log(x))**2 &
         *(1._kp+3._kp*alpha*x**2 &
         -2._kp*alpha*x**2*log(x)+2._kp*alpha**2*x**4+2._kp*alpha**2*x**4*log(x) &
         +8._kp*alpha**2*x**4*(log(x))**2)
    
  end function rcmi_epsilon_two 


!epsilon_three(x)
  function rcmi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: rcmi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    rcmi_epsilon_three = 4._kp/(x**2) &
         /(1._kp-2._kp*alpha*x**2*log(x))**2 &
         *(1._kp-alpha*x**2-4._kp*alpha*x**2*log(x)) &
         *(1._kp-alpha*x**2-9._kp*alpha**2*x**4-4._kp*alpha**3*x**6 &
         -6._kp*alpha*x**2*log(x)-20._kp*alpha**2*x**4*log(x) &
         -6._kp*alpha**3*x**6*log(x) &
         -4._kp*alpha**3*x**6*(log(x))**2-16._kp*alpha**3*x**6*(log(x))**3) &
         *(1._kp+3._kp*alpha*x**2-2._kp*alpha*x**2*log(x)+2._kp*alpha**2*x**4 &
         +2._kp*alpha**2*x**4*log(x)+8._kp*alpha**2*x**4*(log(x))**2)**(-1)
    
  end function rcmi_epsilon_three


!return the field value at which V is maximal
  function rcmi_x_potmax(alpha)
    use specialinf, only : lambert
    implicit none
    real(kp) :: rcmi_x_potmax
    real(kp), intent(in) :: alpha
    
    rcmi_x_potmax = 1._kp/sqrt(2._kp*alpha*lambert(0.5_kp*exp(0.5_kp)/alpha,0))

  end function rcmi_x_potmax


!returns x at the end of inflation defined as epsilon1=1
  function rcmi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: rcmi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rcmiData

    mini = epsilon(1._kp)
    maxi = rcmi_x_potmax(alpha) - epsilon(1._kp)

    rcmiData%real1 = alpha

    rcmi_x_endinf = zbrent(find_rcmi_x_endinf,mini,maxi,tolFind,rcmiData)
   
  end function rcmi_x_endinf


  function find_rcmi_x_endinf(x,rcmiData)
    implicit none
    real(kp), intent(in) :: x    
    type(transfert), optional, intent(inout) :: rcmiData
    real(kp) :: find_rcmi_x_endinf
    real(kp) :: alpha
    
    alpha = rcmiData%real1
    
    find_rcmi_x_endinf = rcmi_epsilon_one(x,alpha) - 1._kp
    
  end function find_rcmi_x_endinf

 


!this is integral[V(phi)/V'(phi) dphi]
  function rcmi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: rcmi_efold_primitive
    type(transfert) :: rcmiData

    real(kp), parameter :: tolInt = tolkp
    integer, parameter :: neq = 1

    real(kp) :: xvar
    real(kp), dimension(neq) :: yvar
    real(kp) :: primapprox

!first order approximation in alpha
!    primapprox = x**2/4._kp &
!         +alpha**2/16._kp*x**4 &
!         +alpha/4._kp*x**4*log(x)


!initial values that matches the first order expansion in alpha
    xvar = epsilon(1._kp)
    yvar(1) = 0._kp

    rcmiData%real1 = alpha

    call easydverk(neq,find_rcmi_efold_primitive,xvar,yvar,x,tolInt,rcmiData)

    rcmi_efold_primitive = yvar(1)

!    print *,'test',rcmi_efold_primitive, primapprox
    
  end function rcmi_efold_primitive


  subroutine find_rcmi_efold_primitive(n,x,y,yprime,rcmiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: rcmiData
    real(kp) :: alpha

    alpha = rcmiData%real1
!regularized to avoid eps1=0
    yprime(1) = 1._kp/sqrt(epsilon(1._kp) + 2._kp*rcmi_epsilon_one(x,alpha))

  end subroutine find_rcmi_efold_primitive



  
  function rcmi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: rcmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: rcmiData

    real(kp) :: efoldMax, xpotmax

    xpotmax = rcmi_x_potmax(alpha)

    efoldMax = -rcmi_efold_primitive(xEnd,alpha) &
         + rcmi_efold_primitive(xpotmax,alpha)

    if (-bfold.gt.efoldMax) then
       write(*,*)'rcmi_x_trajectory: not enough efolds!'
       write(*,*)'efold requested=   efold maxi= ',-bfold,efoldMax
       stop
    endif

    mini = xEnd 
    maxi = rcmi_x_potmax(alpha)

    rcmiData%real1 = alpha
    rcmiData%real2 = -bfold + rcmi_efold_primitive(xend,alpha)
    
    rcmi_x_trajectory = zbrent(find_rcmi_x_trajectory,mini,maxi,tolFind,rcmiData)
       
  end function rcmi_x_trajectory


  function find_rcmi_x_trajectory(x,rcmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: rcmiData
    real(kp) :: find_rcmi_x_trajectory
    real(kp) :: alpha,NplusNuend

    alpha = rcmiData%real1
    NplusNuend = rcmiData%real2

    find_rcmi_x_trajectory = rcmi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_rcmi_x_trajectory


  
end module rcmisr

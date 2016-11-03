!slow-roll functions for the smeared Higgs inflation potential
!
!V(phi) = M**4 * [(1-x**2)**2 + alpha x**4 (log(x)-1/4) +alpha/4 ]
!
!x = phi/phi0
!alpha is a free positive parameter

module shisr
  use infprec, only : kp,tolkp,toldp,transfert
  use inftools, only : zbrent, easydverk
  implicit none

  private

  public  shi_norm_potential, shi_epsilon_one, shi_epsilon_two, shi_epsilon_three
  public  shi_x_endinf, shi_efold_primitive, shi_x_trajectory
  public  shi_norm_deriv_potential, shi_norm_deriv_second_potential
 
contains
!returns V/M**4
  function shi_norm_potential(x,alpha,phi0)
    implicit none
    real(kp) :: shi_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

    shi_norm_potential = (1._kp-x**2)**2+alpha*x**4*(log(x)-0.25_kp)+alpha/4._kp

  end function shi_norm_potential


!returns the first derivative of the potential with respect to x=phi/phi0, divided by M**4
  function shi_norm_deriv_potential(x,alpha,phi0)
    implicit none
    real(kp) :: shi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

   shi_norm_deriv_potential = 4._kp*x*(-1._kp+x**2+x**2*alpha*log(x))

  end function shi_norm_deriv_potential



!returns the second derivative of the potential with respect to x=phi/phi0, divided by M**4
  function shi_norm_deriv_second_potential(x,alpha,phi0)
    implicit none
    real(kp) :: shi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp), intent(in) :: phi0

    shi_norm_deriv_second_potential = 4._kp*(-1._kp+x**2*(3._kp+alpha)+3._kp*x**2*alpha*log(x))

  end function shi_norm_deriv_second_potential



!epsilon_one(x)
  function shi_epsilon_one(x,alpha,phi0)    
    implicit none
    real(kp) :: shi_epsilon_one
    real(kp), intent(in) :: x,alpha,phi0

  
    shi_epsilon_one = ((128._kp*(-x+x**3+x**3*alpha*log(x))**2)/(4._kp-8._kp*x**2-x**4* &
                        (-4._kp+alpha)+alpha+4._kp*x**4*alpha*log(x))**2)/phi0**2

    
  end function shi_epsilon_one


!epsilon_two(x)
  function shi_epsilon_two(x,alpha,phi0)    
    implicit none
    real(kp) :: shi_epsilon_two
    real(kp), intent(in) :: x,alpha,phi0
    
    shi_epsilon_two = ((32._kp*((-1._kp+x**2)*(-4._kp-alpha+x**2*alpha*(6._kp+alpha)+ &
                        x**4*(4._kp-alpha+alpha**2))-x**2*alpha*(4._kp*x**2+x**4*(-8._kp+ &
                        alpha)+3._kp*(4+alpha))*log(x)+4._kp*x**6*alpha**2*log(x)**2))/ &
                        (4._kp-8._kp*x**2-x**4*(-4._kp+alpha)+alpha+4._kp*x**4*alpha*log(x))**2) &
                        /phi0**2
   
  end function shi_epsilon_two


!epsilon_three(x)
  function shi_epsilon_three(x,alpha,phi0)    
    implicit none
    real(kp) :: shi_epsilon_three
    real(kp), intent(in) :: x,alpha,phi0
    
    shi_epsilon_three = ((16._kp*x**2*(-1._kp+x**2+x**2._kp*alpha*log(x))*((-1._kp+x**2)**2* &
                        (-96._kp+80._kp*alpha+46._kp*alpha**2+5._kp*alpha**3+x**4*(32._kp- &
                        16._kp*alpha+14._kp*alpha**2+5._kp*alpha**3)+2._kp*x**2*(32._kp+64._kp* &
                        alpha+30._kp*alpha**2+5._kp*alpha**3))+2._kp*alpha*(80._kp*x**6*alpha &
                        +48._kp*x**2*(4._kp+alpha)+3._kp*(4._kp+alpha)**2-2._kp*x**4*(144._kp+ &
                        68._kp*alpha+5._kp*alpha**2)+x**8*(48._kp-16._kp*alpha+7._kp*alpha**2))* &
                        log(x)-16._kp*x**4*alpha**2*(x**4*(-6._kp+alpha)+6._kp*(4._kp+alpha))* &
                        log(x)**2+32._kp*x**8*alpha**3*log(x)**3))/((4._kp-8._kp*x**2-x**4*(- &
                        4._kp+alpha)+alpha+4._kp*x**4*alpha*log(x))**2*((-1._kp+x**2)*(-4._kp- &
                        alpha+x**2*alpha*(6._kp+alpha)+x**4*(4._kp-alpha+alpha**2))-x**2*alpha* &
                        (4._kp*x**2+x**4*(-8._kp+alpha)+3._kp*(4._kp+alpha))*log(x)+4._kp*x &
                        **6*alpha**2*log(x)**2)))/phi0**2
    
  end function shi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function shi_x_endinf(alpha,phi0)
    implicit none
    real(kp), intent(in) :: alpha,phi0
    real(kp) :: shi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: shiData  

    mini = epsilon(1._kp)
    maxi = 1._kp-epsilon(1._kp)

    shiData%real1 = alpha
    shiData%real2 = phi0	
    
    shi_x_endinf = zbrent(find_shi_x_endinf,mini,maxi,tolFind,shiData)

  end function shi_x_endinf



  function find_shi_x_endinf(x,shiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: shiData
    real(kp) :: find_shi_x_endinf
    real(kp) :: alpha,phi0

    alpha = shiData%real1
    phi0 = shiData%real2
    
    find_shi_x_endinf = shi_epsilon_one(x,alpha,phi0) - 1._kp
   
  end function find_shi_x_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function shi_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: shi_efold_primitive

    type(transfert) :: shiData
!avoids prohibitive integration time in QUADPREC
    real(kp), parameter :: tolInt = max(toldp,tolkp)
    integer, parameter :: neq = 1
    real(kp) :: xvar, xinf
    real(kp), dimension(neq) :: yvar

    !let us start at xend vanishes
    xvar = shi_x_endinf(alpha,phi0)
    yvar(1) = 0._kp

    shiData%real1 = alpha
    shiData%real2 = phi0

    call easydverk(neq,find_shi_efold_primitive,xvar,yvar,x,tolInt,shiData)

    shi_efold_primitive = yvar(1)

  end function shi_efold_primitive

  subroutine find_shi_efold_primitive(n,x,y,yprime,shiData)
    implicit none          
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: shiData
    real(kp) :: alpha, phi0, x4p2n, phi4p2n

    alpha = shiData%real1
    phi0 = shiData%real2

    yprime(1) = phi0**2*((1.-x**2)**2+alpha*x**4*(log(x)-0.25_kp)+ &
                alpha/4._kp)/(4._kp*x*(-1._kp+x**2+x**2*alpha*log(x)))

  end subroutine find_shi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function shi_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold, alpha,phi0, xend
    real(kp) :: shi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: shiData

  
    mini = tolkp
    maxi = shi_x_endinf(alpha,phi0)*(1._kp-epsilon(1._kp))
  

    shiData%real1 = alpha
    shiData%real2 = phi0	
    shiData%real3 = -bfold + shi_efold_primitive(xend,alpha,phi0)
    
    shi_x_trajectory = zbrent(find_shi_x_trajectory,mini,maxi,tolFind,shiData)
       
  end function shi_x_trajectory

  function find_shi_x_trajectory(x,shiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: shiData
    real(kp) :: find_shi_x_trajectory
    real(kp) :: alpha,phi0,NplusNuend

    alpha = shiData%real1
    phi0 = shiData%real2
    NplusNuend = shiData%real3

    find_shi_x_trajectory = shi_efold_primitive(x,alpha,phi0) - NplusNuend
   
  end function find_shi_x_trajectory



end module shisr

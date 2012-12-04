!slow-roll functions for the renormalizable inflection point inflation potential
!
!V(phi) = M**4 (x**2 - alpha x**3 + 9 alpha**2/32 x**4)
!
!x = phi/Mp

module ripisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use cosmopar, only : pi
  implicit none

  private

  public  ripi_norm_potential, ripi_epsilon_one, ripi_epsilon_two, ripi_epsilon_three
  public  ripi_x_endinf, ripi_efold_primitive, ripi_x_trajectory
  public  ripi_norm_deriv_potential, ripi_norm_deriv_second_potential
 
contains
!returns V/M**4
  function ripi_norm_potential(x,alpha)
    implicit none
    real(kp) :: ripi_norm_potential
    real(kp), intent(in) :: x,alpha

    ripi_norm_potential = x**2-alpha*x**3+9._kp/32._kp*alpha**2*x**4

  end function ripi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ripi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha

   ripi_norm_deriv_potential = 1._kp/8._kp*x*(4._kp-3._kp*alpha*x)**2

  end function ripi_norm_deriv_potential



!returns the *1._kp/cosond derivative of the potential with respect to x, divided by M**4
  function ripi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ripi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha

    ripi_norm_deriv_second_potential = 1._kp/8._kp*(-4._kp+3._kp*alpha*x)* &
                                       (-4._kp+9._kp*alpha*x)

  end function ripi_norm_deriv_second_potential



!epsilon_one(x)
  function ripi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: ripi_epsilon_one
    real(kp), intent(in) :: x,alpha
    
    ripi_epsilon_one = (8._kp*(4._kp-3._kp*alpha*x)**4)/(x**2*(32._kp+ &
                       alpha*x*(-32._kp+9._kp*alpha*x))**2)
    
  end function ripi_epsilon_one


!epsilon_two(x)
  function ripi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: ripi_epsilon_two
    real(kp), intent(in) :: x,alpha
    
    ripi_epsilon_two =(8._kp*(-4._kp+3._kp*alpha*x)*(-128._kp+alpha*x* &
                      (160._kp+27._kp*alpha*x*(-4._kp+alpha*x))))/(x**2*(32._kp+ &
                      alpha*x*(-32._kp+9._kp*alpha*x))**2)

  end function ripi_epsilon_two


!epsilon_three(x)
  function ripi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: ripi_epsilon_three
    real(kp), intent(in) :: x,alpha
    
    ripi_epsilon_three = (8._kp*(-4._kp+3._kp*alpha*x)*(16384._kp+3._kp* &
                         alpha*x*(-16384._kp+alpha*x*(20992._kp+alpha*x* &
                         (-15104._kp+27._kp*alpha*x*(256._kp+9._kp*alpha*x* &
                         (-8._kp+alpha*x)))))))/(x**2*(32._kp+alpha*x* &
                         (-32._kp+9._kp*alpha*x))**2*(-128._kp+alpha*x* &
                         (160._kp+27._kp*alpha*x*(-4._kp+alpha*x))))

  end function ripi_epsilon_three


!returns x at the end of inflation defined as epsilon1=1
  function ripi_x_endinf(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: ripi_x_endinf

    ripi_x_endinf = (1._kp/(27._kp*alpha))*2._kp*(16._kp+9._kp*sqrt(2._kp)* &
                    alpha+(2._kp**(2._kp/3._kp)*(-20._kp+9._kp*(2._kp*sqrt(2._kp)- &
                    9._kp*alpha)*alpha))/(544._kp-27._kp*alpha*(11._kp*sqrt(2._kp)+ &
                    9._kp*alpha*(-2._kp+3._kp*sqrt(2._kp)*alpha))+27._kp*sqrt(6._kp)* &
                    sqrt(64._kp+alpha*(-64._kp*sqrt(2._kp)+9._kp*alpha*(11._kp-18._kp* &
                    sqrt(2._kp)*alpha))))**(1._kp/3._kp)-2._kp**(1._kp/3._kp)*(544._kp- &
                    27._kp*alpha*(11._kp*sqrt(2._kp)+9._kp*alpha*(-2._kp+3._kp*sqrt(2._kp)* &
                    alpha))+27._kp*sqrt(6._kp)*sqrt(64._kp+alpha*(-64._kp*sqrt(2._kp)+9._kp* &
                    alpha*(11._kp-18._kp*sqrt(2._kp)*alpha))))**(1._kp/3._kp))
   
  end function ripi_x_endinf



!this is integral(V(phi)/V'(phi) dphi)
  function ripi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ripi_efold_primitive

    if (alpha.eq.0._kp) stop 'ripi_efold_primitive: alpha=0!'

    ripi_efold_primitive = -(2._kp*x/(9._kp*alpha)-x**2/8._kp+16._kp/(27._kp*alpha**2* &
                           (3._kp*alpha*x-4._kp))+4._kp*log(4._kp-3._kp*alpha*x)/(27._kp*alpha**2))

  end function ripi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ripi_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: ripi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ripiData

  
    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = 4._kp/(3._kp*alpha)*(1._kp-epsilon(1._kp)) !Position of the flat inflection point
  


    ripiData%real1 = alpha
    ripiData%real2 = -bfold + ripi_efold_primitive(xend,alpha)
    
    ripi_x_trajectory = zbrent(find_ripitraj,mini,maxi,tolFind,ripiData)
       
  end function ripi_x_trajectory

  function find_ripitraj(x,ripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ripiData
    real(kp) :: find_ripitraj
    real(kp) :: alpha,NplusNuend

    alpha = ripiData%real1
    NplusNuend = ripiData%real2

    find_ripitraj = ripi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_ripitraj


  
end module ripisr

!slow-roll functions for the intermediate inflation potential
!
!V(phi) = M^4 [ x^(-beta) - beta^2/6 x^(-beta-2) ]
!
!x = (phi-phi0)/Mp

module iisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public  ii_norm_potential, ii_epsilon_one, ii_epsilon_two, ii_epsilon_three
  public  ii_efold_primitive, ii_x_trajectory
  public  ii_norm_deriv_potential, ii_norm_deriv_second_potential
  public  ii_prior_xendmin
 
contains
!returns V/M^4
  function ii_norm_potential(x,beta)
    implicit none
    real(kp) :: ii_norm_potential
    real(kp), intent(in) :: x,beta

    ii_norm_potential = x**(-beta)-beta**2/6._kp*x**(-beta-2._kp)

  end function ii_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function ii_norm_deriv_potential(x,beta)
    implicit none
    real(kp) :: ii_norm_deriv_potential
    real(kp), intent(in) :: x,beta

   ii_norm_deriv_potential = -beta*x**(-beta-1._kp) &
         +beta**2*(beta+2._kp)/6._kp *x**(-beta-3._kp)

  end function ii_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function ii_norm_deriv_second_potential(x,beta)
    implicit none
    real(kp) :: ii_norm_deriv_second_potential
    real(kp), intent(in) :: x,beta

    ii_norm_deriv_second_potential = beta*(beta+1._kp)*x**(-beta-2._kp) &
         -beta**2*(beta+2._kp)*(beta+3._kp)/6._kp*x**(-beta-4._kp)

  end function ii_norm_deriv_second_potential



!epsilon_one(x)
  function ii_epsilon_one(x,beta)    
    implicit none
    real(kp) :: ii_epsilon_one
    real(kp), intent(in) :: x,beta
    
    ii_epsilon_one = 0.5_kp*(beta**2*(beta+2._kp)/6._kp-beta*x**2)**2 &
         /(-beta**2*x/6._kp+x**3)**2
    
  end function ii_epsilon_one


!epsilon_two(x)
  function ii_epsilon_two(x,beta)    
    implicit none
    real(kp) :: ii_epsilon_two
    real(kp), intent(in) :: x,beta
    
    ii_epsilon_two =(-2._kp*beta*x**4+beta**2/3._kp &
         *(2._kp*beta+6._kp)*x**2-beta**4/18._kp*(beta+2._kp)) &
         /(x**3-beta**2*x/6._kp)**2
    
  end function ii_epsilon_two


!epsilon_three(x)
  function ii_epsilon_three(x,beta)    
    implicit none
    real(kp) :: ii_epsilon_three
    real(kp), intent(in) :: x,beta
    
    ii_epsilon_three = (beta**5/18._kp*(2._kp+beta) &
         -beta**3*(2._kp+beta)*x**2 +6._kp*beta*(4._kp+beta)*x**4 &
         -12._kp*x**6)*beta*(6._kp*x**2-beta*(2._kp+beta)) &
         /((x**3-beta**2/6._kp*x)**2 &
         *(beta**3*(beta+2._kp)-12._kp*beta*(beta+3._kp)*x**2+36._kp*x**4))
    
  end function ii_epsilon_three

!returns x at the end of inflation
  function ii_x_endinf(beta,xend)
    implicit none
    real(kp), intent(in) :: beta,xend
    real(kp) :: ii_x_endinf
   
    ii_x_endinf = xend
   
  end function ii_x_endinf


!this is integral[V(phi)/V'(phi) dphi]
  function ii_efold_primitive(x,beta)
    implicit none
    real(kp), intent(in) :: x,beta
    real(kp) :: ii_efold_primitive

    if (beta.eq.0._kp) stop 'ii_efold_primitive: beta=0!'

    ii_efold_primitive = -1._kp/(2._kp*beta)*x**2 &
         -1._kp/6._kp*log(x**2-beta*(beta+2._kp)/6._kp)

  end function ii_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ii_x_trajectory(bfold,xend,beta)
    implicit none
    real(kp), intent(in) :: bfold, beta, xend
    real(kp) :: ii_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: iiData

  
    mini = sqrt(beta/2._kp*(1._kp+beta/3._kp+sqrt(1._kp+4._kp*beta/9._kp)))
    maxi = xend
  


    iiData%real1 = beta
    iiData%real2 = -bfold + ii_efold_primitive(xend,beta)
    
    ii_x_trajectory = zbrent(find_ii_x_trajectory,mini,maxi,tolFind,iiData)
       
  end function ii_x_trajectory

  function find_ii_x_trajectory(x,iiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: iiData
    real(kp) :: find_ii_x_trajectory
    real(kp) :: beta,NplusNuend

    beta = iiData%real1
    NplusNuend = iiData%real2

    find_ii_x_trajectory = ii_efold_primitive(x,beta) - NplusNuend
   
  end function find_ii_x_trajectory

 
!Returns the minimum value of xend in order to have at least -bfolds possible
  function ii_prior_xendmin(beta,bfold)
    implicit none
    real(kp), intent(in) :: beta,bfold
    real(kp) :: ii_prior_xendmin
    complex(kp) ::xin


    if (beta .gt. 9._kp/2._kp*(1._kp+sqrt(2._kp))) then
      xin=beta/(3._kp*sqrt(2._kp)) &
         +cmplx(1._kp,sqrt(3._kp),kp) &
         /(3._kp*sqrt(2._kp)) &
         *beta**(4._kp/3._kp) &
         /((cmplx(9._kp+2._kp*beta,sqrt(-81._kp-36._kp*beta+4._kp*beta**2)))**(1._kp/3._kp),kp) &
         +cmplx(1._kp,-sqrt(3._kp),kp) &
         *beta**(2._kp/3._kp)/(6._kp*sqrt(2._kp))* &
         (cmplx(9._kp+2._kp*beta,sqrt(-81._kp-36._kp*beta+4._kp*beta**2),kp))**(1._kp/3._kp)


    else
      xin=sqrt(beta/2._kp*(1._kp+beta/3._kp+sqrt(1._kp+4._kp*beta/9._kp)))

    endif

    ii_prior_xendmin=sqrt((real(xin,kp))**2-2._kp*beta*bfold)

   end function ii_prior_xendmin

  
end module iisr

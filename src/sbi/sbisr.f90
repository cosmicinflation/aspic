!slow-roll functions for the supergravity brane inflation potential
!
!V(phi) = M**4 { 1 + [-alpha + beta ln(x) ] x**4 }
!
!x = phi/Mp

module sbisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : ei,lambert
  implicit none

  private

  public sbi_norm_potential, sbi_epsilon_one, sbi_epsilon_two, sbi_epsilon_three
  public sbi_x_endinf, sbi_efold_primitive, sbi_x_trajectory
  public sbi_norm_deriv_potential, sbi_norm_deriv_second_potential
  public sbi_alphamin, sbi_x_potmin, sbi_x_potzero

 
contains
!returns V/M**4 as function of x
  function sbi_norm_potential(x,alpha,beta)
    implicit none
    real(kp) :: sbi_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    sbi_norm_potential = 1._kp+(-alpha+beta*log(x))*x**4

  end function sbi_norm_potential



!returns the first derivative of the potential with respect to x
  function sbi_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: sbi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta

   sbi_norm_deriv_potential = x**3*(-4._kp*alpha+beta+4._kp*beta*log(x))

  end function sbi_norm_deriv_potential



!returns the second derivative of the potential with respect to x
  function sbi_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: sbi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta

    sbi_norm_deriv_second_potential = x**2*(-12._kp*alpha+7._kp*beta+12._kp*beta*log(x))

  end function sbi_norm_deriv_second_potential



!epsilon_one(x)
  function sbi_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: sbi_epsilon_one
    real(kp), intent(in) :: x,alpha,beta
    
    sbi_epsilon_one = (x**6*(-4._kp*alpha+beta+4._kp*beta*log(x))**2)/ &
         (2._kp*(1._kp-alpha*x**4+beta*x**4*log(x))**2)
    
  end function sbi_epsilon_one


!epsilon_two(x)
  function sbi_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: sbi_epsilon_two
    real(kp), intent(in) :: x,alpha,beta
    
    sbi_epsilon_two = (2._kp*(x**6*(-4._kp*alpha+beta+4._kp*beta*log(x))**2- &
         x**2*(-12._kp*alpha+7._kp*beta+12._kp*beta*log(x))* &
         (1._kp+x**4*(-alpha+beta*log(x)))))/(1._kp+x**4*(-alpha+beta*log(x)))**2
    
  end function sbi_epsilon_two


!epsilon_three(x)
  function sbi_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: sbi_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
    
    sbi_epsilon_three = (1._kp/(x**2))*(8._kp+(2._kp*(-4._kp+beta*x**4)**2)/ &
         (1._kp-alpha*x**4+beta*x**4*log(x))**2+(-52._kp+ &
         9._kp*beta*x**4)/(1._kp-alpha*x**4+beta*x**4*log(x))+ &
         (144._kp*alpha-84._kp*beta+(28._kp*alpha-11._kp*beta)* &
         beta*x**4-4._kp*beta*(36._kp+7._kp*beta*x**4)*log(x))/ &
         (12._kp*alpha-7._kp*beta+(4._kp*alpha**2-alpha*beta+beta**2)* & 
         x**4+beta*log(x)*(-12._kp+(-8._kp*alpha+beta)*x**4+4._kp*beta*x**4*log(x))))
    
  end function sbi_epsilon_three


!field value at which the potential is minimal
  function sbi_x_potmin(alpha,beta)
    implicit none
    real(kp) , intent(in) :: alpha,beta
    real(kp) :: sbi_x_potmin

    sbi_x_potmin = exp(-0.25_kp+alpha/beta)

  end function sbi_x_potmin


!Returns the minimum value for alpha such that the minimum value of the potential is negative
  function sbi_alphamin(beta)
    implicit none
    real(kp) , intent(in) :: beta
    real(kp) :: sbi_alphamin

    sbi_alphamin = beta/4._kp*(1._kp-log(beta/4._kp))

  end function sbi_alphamin


!field value at which the potential vanishes, for 0<x<x_{V'=0}
  function sbi_x_potzero(alpha,beta)
    implicit none
    real(kp) , intent(in) :: alpha,beta
    real(kp) :: sbi_x_potzero,L1,L2,W

    if (alpha .lt. sbi_alphamin(beta)) stop 'sbi_x_potzero: alpha < alphamin !'

    if (alpha .gt. sbi_alphamin(beta)) then

      sbi_x_potzero = (-4._kp/(beta*lambert(-4._kp/beta*exp(-4._kp*alpha/beta),-1)))**(0.25_kp)
   
    else

      sbi_x_potzero = (0.25_kp*beta)**(-0.25_kp)

    endif


  end function sbi_x_potzero



!returns x at the end of inflation defined as epsilon1=1
  function sbi_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: sbi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sbiData

    mini = epsilon(1._kp)

    if (alpha.eq.sbi_alphamin(beta)) then
       maxi = sbi_x_potzero(alpha,beta)*(1-sqrt(epsilon(1._kp)))
    else
       maxi = sbi_x_potzero(alpha,beta)*(1._kp-epsilon(1._kp))
    endif

    sbiData%real1 = alpha
    sbiData%real2 = beta
    
    sbi_x_endinf = zbrent(find_sbi_x_endinf,mini,maxi,tolFind,sbiData)
   
  end function sbi_x_endinf

  function find_sbi_x_endinf(x,sbiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: sbiData
    real(kp) :: find_sbi_x_endinf
    real(kp) :: alpha,beta

    alpha = sbiData%real1
    beta = sbiData%real2

    find_sbi_x_endinf = sbi_epsilon_one(x,alpha,beta)-1._kp
   
  end function find_sbi_x_endinf


!this is integral(V(phi)/V'(phi) dphi)
  function sbi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: sbi_efold_primitive

    if (beta.eq.0._kp) stop 'sbi_efold_primitive: beta=0 !'

    sbi_efold_primitive = exp(0.5_kp-2._kp*alpha/beta)/(4._kp*beta) &
         * ei(-0.5_kp+2._kp*alpha/beta-2._kp*log(x)) &
         - exp(2._kp*alpha/beta-0.5_kp)/16._kp*ei(0.5_kp-2._kp*alpha/beta+2._kp*log(x)) &
         +x**2/8._kp

  end function sbi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function sbi_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold,alpha,beta,xend
    real(kp) :: sbi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: sbiData


    maxi = xend*(1._kp-epsilon(1._kp))
    mini = epsilon(1._kp)

    sbiData%real1 = alpha
    sbiData%real2 = beta
    sbiData%real3 = -bfold + sbi_efold_primitive(xend,alpha,beta)
    
    sbi_x_trajectory = zbrent(find_sbi_x_trajectory,mini,maxi,tolFind,sbiData)
       
  end function sbi_x_trajectory

  function find_sbi_x_trajectory(x,sbiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: sbiData
    real(kp) :: find_sbi_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha= sbiData%real1
    beta = sbiData%real2
    NplusNuend = sbiData%real3

    find_sbi_x_trajectory = sbi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_sbi_x_trajectory


  
end module sbisr

!slow-roll functions for the arctan potential
!
!V(phi) = M**4 * [1 - 2/pi arctan(x)]
!
!x = phi/mu
!mu=mu/Mp

module aisr
  use infprec, only : kp,tolkp,transfert, pi
  use inftools, only : zbrent
  implicit none

  private

  public ai_norm_potential, ai_epsilon_one, ai_epsilon_two, ai_epsilon_three
  public ai_x_endinf, ai_efold_primitive, ai_x_trajectory
  public ai_norm_deriv_potential, ai_norm_deriv_second_potential
  public ai_mumax, ai_x_epsonemax, ai_numacc_xinimin
 
contains
!returns V/M**4
  function ai_norm_potential(x,mu)
    implicit none
    real(kp) :: ai_norm_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: mu

    ai_norm_potential = 1._kp-2._kp/pi*atan(x)

  end function ai_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function ai_norm_deriv_potential(x,mu)
    implicit none
    real(kp) :: ai_norm_deriv_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: mu

   ai_norm_deriv_potential = -(2._kp/(pi+pi*x**2))

  end function ai_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function ai_norm_deriv_second_potential(x,mu)
    implicit none
    real(kp) :: ai_norm_deriv_second_potential
    real(kp), intent(in) :: x
    real(kp), intent(in) :: mu

    ai_norm_deriv_second_potential = (4._kp*x)/(pi*(1._kp+x**2)**2)

  end function ai_norm_deriv_second_potential



!epsilon_one(x=phi/mu)
  function ai_epsilon_one(x,mu)    
    implicit none
    real(kp) :: ai_epsilon_one
    real(kp), intent(in) :: x,mu

  
    ai_epsilon_one = 2._kp/(mu**2*(1._kp+x**2)**2*(pi-2._kp*atan(x))**2)

    
  end function ai_epsilon_one


!epsilon_two(x=phi/mu)
  function ai_epsilon_two(x,mu)    
    implicit none
    real(kp) :: ai_epsilon_two
    real(kp), intent(in) :: x,mu
    
    ai_epsilon_two = (8._kp/mu**2*(1._kp-pi*x+2._kp*x*atan(x)))/ &
         ((1._kp+x**2)**2*(pi-2._kp*atan(x))**2)
    
  end function ai_epsilon_two


!epsilon_three(x=phi/mu)
  function ai_epsilon_three(x,mu)    
    implicit none
    real(kp) :: ai_epsilon_three
    real(kp), intent(in) :: x,mu
    
    ai_epsilon_three = (2._kp/mu**2*(-4._kp+pi*(pi+6._kp*x-3._kp*pi*x**2)- & 
         4._kp*atan(x)*(pi+3._kp*x-3._kp*pi*x**2+(-1._kp+3._kp*x**2)* &
         atan(x))))/((1._kp+x**2)**2*(pi-2._kp*atan(x))**2* &
         (-1._kp+pi*x-2._kp*x*atan(x)))
    
  end function ai_epsilon_three



  function ai_x_epsonemax(mu)
    implicit none
    real(kp) :: ai_x_epsonemax
    real(kp), intent(in) :: mu

    ai_x_epsonemax= 0.428977908964179282343424050497337_kp

  end function ai_x_epsonemax



  function ai_mumax()
    implicit none
    real(kp) :: ai_mumax

    ai_mumax = 0.512377601980274977255003526679718_kp

  end function ai_mumax



  function ai_numacc_xinimin(mu)
    implicit none
    real(kp) :: ai_numacc_xinimin
    real(kp), intent(in) :: mu

!at first order in 0
!sqrt(eps1) = sqrt(2)/pi (1+2/pi x) + O(x^2)
    
    ai_numacc_xinimin = -sqrt( 1._kp/ (mu*pi*sqrt(2._kp*epsilon(1._kp))))

  end function ai_numacc_xinimin



!returns x=phi/mu at the end of inflation defined as epsilon1=1
  function ai_x_endinf(mu)
    implicit none
    real(kp), intent(in) ::mu
    real(kp) :: ai_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: aiData
    real(kp) :: xepsonemax

    xepsonemax = ai_x_epsonemax(mu)

    if (ai_epsilon_one(xepsonemax,mu).lt.1._kp) then
       stop 'ai_x_endinf: epsonemax<1, inflation never stops!'
    endif
   
! Minimal numerical bound ensuring eps1>numacc
    mini = ai_numacc_xinimin(mu)    
!Position where epsilon1 is maximum
    maxi= xepsonemax

    aiData%real1 = mu

    ai_x_endinf=zbrent(find_ai_x_endinf,mini,maxi,tolFind,aiData)

   
  end function ai_x_endinf

  function find_ai_x_endinf(x,aiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: aiData
    real(kp) :: find_ai_x_endinf
    real(kp) :: mu

    mu = aiData%real1
    
    find_ai_x_endinf = ai_epsilon_one(x,mu) - 1._kp
  
  end function find_ai_x_endinf



!tais is integral[V(phi)/V'(phi) dphi]
  function ai_efold_primitive(x,mu)
    implicit none
    real(kp), intent(in) :: x,mu
    real(kp) :: ai_efold_primitive

    ai_efold_primitive = -mu**2*(pi*x/2._kp+x**2/6._kp+ &
                         pi*x**3/6._kp-(1._kp+x**2/3._kp)*x*atan(x)+ & 
                         log(1._kp+x**2)/3._kp)

  end function ai_efold_primitive


!returns x=phi/mi at bfold=-efolds before the end of inflation, ie N-Nend
  function ai_x_trajectory(bfold,xend,mu)
    implicit none
    real(kp), intent(in) :: bfold, xend, mu
    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: ai_x_trajectory,mini,maxi
    real(kp) :: xinimin, efoldMax
    type(transfert) :: aiData

    aiData%real1 = bfold
    aiData%real2 = xend
    aiData%real3 = mu


    xinimin = ai_numacc_xinimin(mu)

    efoldMax = - ai_efold_primitive(xend,mu) &
         + ai_efold_primitive(xinimin,mu)

    if (-bfold.gt.efoldMax) then
       write(*,*)'ai_x_trajectory: not enough efolds!'
       write(*,*)'efold requested= efold maxi= ',-bfold,efoldMax
       write(*,*) 'xinimin (numacc)= ',xinimin
       stop 
    endif

!    mini = -mu**(-2._kp/3._kp)*10._kp**(10._kp)
    mini = xinimin
    maxi = xend
    
    ai_x_trajectory =zbrent(find_ai_x_trajectory,mini,maxi,tolzbrent,aiData)
       
  end function ai_x_trajectory

  function find_ai_x_trajectory(x,aiData)   
    implicit none
    real(kp) :: find_ai_x_trajectory
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: aiData

    real(kp) :: bfold,xend,mu

    bfold = aiData%real1
    xend = aiData%real2
    mu = aiData%real3

    find_ai_x_trajectory = ai_efold_primitive(x,mu) - ai_efold_primitive(xend,mu) + bfold
  
  end function find_ai_x_trajectory


end module aisr

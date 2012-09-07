!slow-roll functions for the Twisted inflation potential
!
!V(phi) = M^4 [ 1 - A (x/phi0)**2 exp(-x/phi0) ]
!
!x = phi/Mp
!y = phi/phi0
!A=32/(92 ζ(5))=0.33183220315845589991447280462623125031833134154035

module twisr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : ei
  use inftools, only : zbrent
  implicit none

  private

  public  twi_norm_potential, twi_epsilon_one, twi_epsilon_two, twi_epsilon_three
  public  twi_x_endinf, twi_efold_primitive, twi_x_trajectory
  public  twi_norm_deriv_potential, twi_norm_deriv_second_potential

  real(kp), parameter ::A=0.33183220315845589991447280462623125031833134154035_kp

!if bigger numerical errors
  real(kp), parameter :: ymax = 100._kp
 
contains



!returns V/M**4
  function twi_norm_potential(x,phi0)
    implicit none
    real(kp) :: twi_norm_potential
    real(kp), intent(in) :: x,phi0

    twi_norm_potential = 1._kp-A*(x/phi0)**2*exp(-x/phi0)

  end function twi_norm_potential


!returns the first derivative of the potential with respect to x, divided by M**4
  function twi_norm_deriv_potential(x,phi0)
    implicit none
    real(kp) :: twi_norm_deriv_potential
    real(kp), intent(in) :: x,phi0

   twi_norm_deriv_potential = (A*exp(-x/phi0)*x*(-2._kp*phi0+x))/phi0**3

  end function twi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M**4
  function twi_norm_deriv_second_potential(x,phi0)
    implicit none
    real(kp) :: twi_norm_deriv_second_potential
    real(kp), intent(in) :: x,phi0

    twi_norm_deriv_second_potential = -((A*exp(-x/phi0)*(2._kp*phi0**2-4._kp*phi0*x &
                                      +x**2))/phi0**4)

  end function twi_norm_deriv_second_potential



!epsilon_one(x)
  function twi_epsilon_one(x,phi0)    
    implicit none
    real(kp) :: twi_epsilon_one
    real(kp), intent(in) :: x,phi0

  
    twi_epsilon_one = (A**2*x**2*(-2._kp*phi0+x)**2)/ &
                      (2._kp*(exp(x/phi0)*phi0**3-A*phi0*x**2)**2)

    
  end function twi_epsilon_one


!epsilon_two(x)
  function twi_epsilon_two(x,phi0)    
    implicit none
    real(kp) :: twi_epsilon_two
    real(kp), intent(in) :: x,phi0
    
    twi_epsilon_two = (2._kp*A*(2._kp*A*x**2+exp(x/phi0)*(2._kp*phi0**2-4._kp*phi0*x &
                      +x**2)))/(exp(x/phi0)*phi0**2-A*x**2)**2
    
  end function twi_epsilon_two


!epsilon_three(x)
  function twi_epsilon_three(x,phi0)    
    implicit none
    real(kp) :: twi_epsilon_three
    real(kp), intent(in) :: x,phi0
    
    twi_epsilon_three = (A*(2._kp*phi0-x)*x*(4._kp*A**2*phi0*x**3 - & 
                        exp((2._kp*x)/phi0)*phi0**2*(6._kp*phi0**2-6*phi0*x+x**2)- & 
                        A*exp(x/phi0)*x*(-12._kp*phi0**3+18._kp*phi0**2*x-6._kp*phi0 &
                        *x**2+x**3)))/((exp(x/phi0)*phi0**3-A*phi0*x**2)**2 &
                        *(2._kp*A*x**2+exp(x/phi0)*(2._kp*phi0**2-4._kp*phi0*x+x**2)))
    
  end function twi_epsilon_three


!returns x at the end of slow-roll epsilon2=2
  function twi_x_endinf(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: twi_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: twiData
    
    if (phi0.gt.0.04228_kp) stop 'twi_x_endinf: no slow roll solution for this value of phi0!'

    maxi=ymax*phi0 !if bigger, <numaccuracy errors
    mini = twi_x_eps1max(phi0) !second maximum of eps1

    twiData%real1 = phi0
    
    twi_x_endinf = zbrent(find_twiendinf,mini,maxi,tolFind,twiData)


  end function twi_x_endinf



  function find_twiendinf(x,twiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: twiData
    real(kp) :: find_twiendinf
    real(kp) :: phi0

    phi0 = twiData%real1
    
    find_twiendinf = twi_epsilon_one(x,phi0) - 1._kp
  
  end function find_twiendinf


!returns the position x of the second maximum of epsilon1
  function twi_x_eps1max(phi0)
    implicit none
    real(kp), intent(in) :: phi0
    real(kp) :: twi_x_eps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: twiData

    maxi=ymax*phi0 !if bigger, <numaccuracy errors
    mini = 2._kp*phi0*(1._kp+epsilon(1._kp)) !minimum of the potential (where eps1=0)

    twiData%real1 = phi0
    
    twi_x_eps1max = zbrent(find_twieps1max,mini,maxi,tolFind,twiData)


  end function twi_x_eps1max



  function find_twieps1max(x,twiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: twiData
    real(kp) :: find_twieps1max
    real(kp) :: phi0

    phi0 = twiData%real1
    
    find_twieps1max = 2._kp*A*(x/phi0)**2+exp(x/phi0)*((x/phi0)**2-4._kp*x/phi0+2._kp)
   
  end function find_twieps1max




!this is integral[V(phi)/V'(phi) dphi]
  function twi_efold_primitive(x,phi0)
    implicit none
    real(kp), intent(in) :: x,phi0
    real(kp) :: twi_efold_primitive


    twi_efold_primitive = phi0**2*(-ei(x/phi0)/(2._kp*A)+exp(2._kp)*ei(x/phi0-2._kp) &
                          /(2._kp*A)-x/phi0-2._kp*log(x/phi0-2._kp))

  end function twi_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function twi_x_trajectory(bfold,xend,phi0)
    implicit none
    real(kp), intent(in) :: bfold, phi0, xend
    real(kp) :: twi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: twiData

  
    mini = xend
    maxi=ymax*phi0 !if bigger, <numaccuracy errors

    twiData%real1 = phi0
    twiData%real2 = -bfold + twi_efold_primitive(xend,phi0)
    
    twi_x_trajectory = zbrent(find_twitraj,mini,maxi,tolFind,twiData)
       
  end function twi_x_trajectory

  function find_twitraj(x,twiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: twiData
    real(kp) :: find_twitraj
    real(kp) :: phi0,NplusNuend

    phi0 = twiData%real1
    NplusNuend = twiData%real2

    find_twitraj = twi_efold_primitive(x,phi0) - NplusNuend
   
  end function find_twitraj



end module twisr

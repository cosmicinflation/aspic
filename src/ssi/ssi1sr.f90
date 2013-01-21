!slow-roll functions for the sneutrino supersymmetric 1 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!1: alpha>0, beta>0
!
!x = phi/Mp

module ssi1sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential
  use ssicommon, only : ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssi_x_trajectory, ssi136_x_epstwozero


  implicit none

  private

  public ssi1_norm_potential
  public ssi1_epsilon_one, ssi1_epsilon_two, ssi1_epsilon_three
  public ssi1_x_endinf, ssi1_efold_primitive, ssi1_x_trajectory
  public ssi1_norm_deriv_potential, ssi1_norm_deriv_second_potential
  public ssi1_alphamin, ssi1_x_epsonemax
  
contains

!returns V/M**4
  function ssi1_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi1_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi1_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi1_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi1_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi1_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi1_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi1_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi1_norm_deriv_second_potential



!epsilon_one(x)
  function ssi1_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi1_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi1_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi1_epsilon_one


!epsilon_two(x)
  function ssi1_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi1_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi1_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi1_epsilon_two


!epsilon_three(x)
  function ssi1_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi1_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi1_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi1_epsilon_three


!returns the position x where epsilon_one is maximum
  function ssi1_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssi1_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssi1_x_epsonemax = ssi136_x_epstwozero(alpha,beta)
    
  end function ssi1_x_epsonemax


!returns the minimum value of alpha (given beta) in order for
!inflation to end by slow roll violation (eps1max>1)
  function ssi1_alphamin(beta)    
    implicit none
    real(kp) :: ssi1_alphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi1Data


    mini=1._kp
    maxi=10._kp**(3._kp)

    ssi1Data%real1 = beta


!In that case inflation ends by slow roll violation for any value of alpha
    if(beta .gt. 0.251_kp) then

       ssi1_alphamin = 0._kp
 
    else

       ssi1_alphamin = zbrent(find_ssi1alphamin,mini,maxi,tolFind,ssi1Data)

    endif
    
  end function ssi1_alphamin

  function find_ssi1alphamin(alpha,ssi1Data)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: ssi1Data
    real(kp) :: find_ssi1alphamin
    real(kp) :: beta

    beta = ssi1Data%real1
    
    find_ssi1alphamin = ssi1_epsilon_one(ssi1_x_epsonemax(alpha,beta),alpha,beta)-1._kp
   
  end function find_ssi1alphamin


!returns x at the end of inflation defined as epsilon1=1
  function ssi1_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi1Data

    if (alpha .lt. ssi1_alphamin(beta)) then
       stop 'ssi1_x_endinf: epsilon1max < 1, inflation cannot stop by slow roll violation!'
    endif

    mini = ssi1_x_epsonemax(alpha,beta)
    maxi = mini/epsilon(1._kp)

!    print*,'ss1_xend:  mini=',mini,'   maxi=',maxi,'   epsOne(mini)=',ssi1_epsilon_one(mini,alpha,beta), &
!                 '   epsOne(maxi)=',ssi1_epsilon_one(maxi,alpha,beta)
!    pause

    ssi1Data%real1 = alpha
    ssi1Data%real2 = beta
    
    ssi1_x_endinf = zbrent(find_ssi1endinf,mini,maxi,tolFind,ssi1Data)

  end function ssi1_x_endinf



  function find_ssi1endinf(x,ssi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi1Data
    real(kp) :: find_ssi1endinf
    real(kp) :: alpha,beta

    alpha = ssi1Data%real1
    beta = ssi1Data%real2
    
    find_ssi1endinf = ssi1_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi1endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi1_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi1_efold_primitive

    ssi1_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi1_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi1Data


    mini = ssi1_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)
  
    ssi1Data%real1 = alpha
    ssi1Data%real2 = beta
    ssi1Data%real3 = -bfold + ssi1_efold_primitive(xend,alpha,beta)
    
    ssi1_x_trajectory = zbrent(find_ssi_x_trajectory,mini,maxi,tolFind,ssi1Data)
       
  end function ssi1_x_trajectory
 

end module ssi1sr

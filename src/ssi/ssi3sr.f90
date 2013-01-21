!slow-roll functions for the sneutrino supersymmetric 3 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!3: alpha>0, beta<0
!
!x = phi/Mp

module ssi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential
  use ssicommon, only : ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssi_x_trajectory, ssi136_x_epstwozero
  use ssicommon, only : ssi3456_x_derivpotzero


  implicit none

  private

  public ssi3_norm_potential
  public ssi3_epsilon_one, ssi3_epsilon_two, ssi3_epsilon_three
  public ssi3_x_endinf, ssi3_efold_primitive, ssi3_x_trajectory
  public ssi3_norm_deriv_potential, ssi3_norm_deriv_second_potential
  public ssi3_alphamin, ssi3_x_epsonemax, ssi3_x_potmax
  
contains

!returns V/M**4
  function ssi3_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi3_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi3_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi3_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi3_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi3_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi3_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi3_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi3_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi3_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi3_norm_deriv_second_potential



!epsilon_one(x)
  function ssi3_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi3_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi3_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi3_epsilon_one


!epsilon_two(x)
  function ssi3_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi3_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi3_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi3_epsilon_two


!epsilon_three(x)
  function ssi3_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi3_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi3_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi3_epsilon_three

!returns the position x where the potential is maximal
  
  function ssi3_x_potmax(alpha,beta)    
    implicit none
    real(kp) :: ssi3_x_potmax
    real(kp), intent(in) :: alpha,beta

    ssi3_x_potmax = ssi3456_x_derivpotzero(alpha,beta)
    
  end function ssi3_x_potmax


!returns the position x where epsilon_one is maximum
  function ssi3_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssi3_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssi3_x_epsonemax = ssi136_x_epstwozero(alpha,beta)
    
  end function ssi3_x_epsonemax



!returns the minimum value of alpha (given beta) in order for inflation to end by slow roll violation (eps1max>1)
  function ssi3_alphamin(beta)    
    implicit none
    real(kp) :: ssi3_alphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi3Data

    mini=epsilon(1._kp)
    maxi=10._kp**(4._kp)

    ssi3Data%real1 = beta

!    print*,'ssi3_alphamin:  beta=',beta,'mini=',mini,'  xeps2mini=',ssi3_x_epsonemax(mini,beta), &
!           'maxi=',maxi,'  xeps2maxi=',ssi3_x_epsonemax(maxi,beta), &
!           '  eps1(xeps2mini)=',ssi3_epsilon_one(ssi3_x_epsonemax(mini,beta),mini,beta), &
!           '  eps1(xeps2maxi)=',ssi3_epsilon_one(ssi3_x_epsonemax(maxi,beta),maxi,beta)
!    pause

       ssi3_alphamin = zbrent(find_ssi3_alphamin,mini,maxi,tolFind,ssi3Data)
    
  end function ssi3_alphamin

  function find_ssi3_alphamin(alpha,ssi3Data)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: ssi3Data
    real(kp) :: find_ssi3_alphamin
    real(kp) :: beta

    beta = ssi3Data%real1
    
    find_ssi3_alphamin = ssi3_epsilon_one(ssi3_x_epsonemax(alpha,beta),alpha,beta)-1._kp
   
  end function find_ssi3_alphamin


!returns x at the end of inflation defined as epsilon1=1
  function ssi3_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi3_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi3Data

    if (alpha .lt. ssi3_alphamin(beta)) then
       stop 'ssi3_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'
    endif

    mini = ssi3_x_epsonemax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi3_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))


    ssi3Data%real1 = alpha
    ssi3Data%real2 = beta
    
    ssi3_x_endinf = zbrent(find_ssi3_x_endinf,mini,maxi,tolFind,ssi3Data)

  end function ssi3_x_endinf



  function find_ssi3_x_endinf(x,ssi3Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi3Data
    real(kp) :: find_ssi3_x_endinf
    real(kp) :: alpha,beta

    alpha = ssi3Data%real1
    beta = ssi3Data%real2
    
    find_ssi3_x_endinf = ssi3_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi3_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi3_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi3_efold_primitive

    ssi3_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi3_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi3_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi3Data


    mini = ssi3_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssi3_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssi3Data%real1 = alpha
    ssi3Data%real2 = beta
    ssi3Data%real3 = -bfold + ssi3_efold_primitive(xend,alpha,beta)
    
    ssi3_x_trajectory = zbrent(find_ssi_x_trajectory,mini,maxi,tolFind,ssi3Data)
       
  end function ssi3_x_trajectory
 

end module ssi3sr

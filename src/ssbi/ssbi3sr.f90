!slow-roll functions for the spontaneous symmetry breaking 3 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!3: alpha>0, beta<0
!
!x = phi/Mp

module ssbi3sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory, ssbi136_x_epstwozero
  use ssbicommon, only : ssbi3456_x_derivpotzero


  implicit none

  private

  public ssbi3_norm_potential
  public ssbi3_epsilon_one, ssbi3_epsilon_two, ssbi3_epsilon_three
  public ssbi3_x_endinf, ssbi3_efold_primitive, ssbi3_x_trajectory
  public ssbi3_norm_deriv_potential, ssbi3_norm_deriv_second_potential
  public ssbi3_alphamin, ssbi3_x_epsonemax, ssbi3_x_potmax, ssbi3_eps2_x_potmax
  
contains

!returns V/M**4
  function ssbi3_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi3_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi3_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi3_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi3_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi3_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi3_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi3_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi3_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi3_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi3_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi3_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi3_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi3_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi3_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi3_epsilon_one


!epsilon_two(x)
  function ssbi3_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi3_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi3_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi3_epsilon_two


!epsilon_three(x)
  function ssbi3_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi3_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi3_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi3_epsilon_three

!returns the position x where the potential is maximal
  
  function ssbi3_x_potmax(alpha,beta)    
    implicit none
    real(kp) :: ssbi3_x_potmax
    real(kp), intent(in) :: alpha,beta

    ssbi3_x_potmax = ssbi3456_x_derivpotzero(alpha,beta)
    
  end function ssbi3_x_potmax

! returns the value of eps2 at the maximum of the potential
  function ssbi3_eps2_x_potmax(alpha,beta)    
    implicit none
    real(kp) :: ssbi3_eps2_x_potmax
    real(kp), intent(in) :: alpha,beta

    ssbi3_eps2_x_potmax = ssbi3_epsilon_two(ssbi3_x_potmax(alpha,beta),alpha,beta)
    
  end function ssbi3_eps2_x_potmax


!returns the position x where epsilon_one is maximum
  function ssbi3_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssbi3_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssbi3_x_epsonemax = ssbi136_x_epstwozero(alpha,beta)
    
  end function ssbi3_x_epsonemax



!returns the minimum value of alpha (given beta) in order for inflation to end by slow roll violation (eps1max>1)
  function ssbi3_alphamin(beta)    
    implicit none
    real(kp) :: ssbi3_alphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi3Data

    mini=epsilon(1._kp)
    maxi=10._kp**(4._kp)

    ssbi3Data%real1 = beta


       ssbi3_alphamin = zbrent(find_ssbi3_alphamin,mini,maxi,tolFind,ssbi3Data)
    
  end function ssbi3_alphamin

  function find_ssbi3_alphamin(alpha,ssbi3Data)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: ssbi3Data
    real(kp) :: find_ssbi3_alphamin
    real(kp) :: beta

    beta = ssbi3Data%real1
    
    find_ssbi3_alphamin = ssbi3_epsilon_one(ssbi3_x_epsonemax(alpha,beta),alpha,beta)-1._kp
   
  end function find_ssbi3_alphamin


!returns x at the end of inflation defined as epsilon1=1
  function ssbi3_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi3_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi3Data

    if (alpha .lt. ssbi3_alphamin(beta)) then
       stop 'ssbi3_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'
    endif

    mini = ssbi3_x_epsonemax(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi3_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))


    ssbi3Data%real1 = alpha
    ssbi3Data%real2 = beta
    
    ssbi3_x_endinf = zbrent(find_ssbi3_x_endinf,mini,maxi,tolFind,ssbi3Data)

  end function ssbi3_x_endinf



  function find_ssbi3_x_endinf(x,ssbi3Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi3Data
    real(kp) :: find_ssbi3_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi3Data%real1
    beta = ssbi3Data%real2
    
    find_ssbi3_x_endinf = ssbi3_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi3_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi3_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi3_efold_primitive

    ssbi3_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi3_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi3_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi3_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi3Data


    mini = ssbi3_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = ssbi3_x_potmax(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssbi3Data%real1 = alpha
    ssbi3Data%real2 = beta
    ssbi3Data%real3 = -bfold + ssbi3_efold_primitive(xend,alpha,beta)
    
    ssbi3_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi3Data)
       
  end function ssbi3_x_trajectory
 

end module ssbi3sr

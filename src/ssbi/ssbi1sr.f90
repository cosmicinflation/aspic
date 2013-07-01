!slow-roll functions for the spontaneous symmetry breaking 1 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!1: alpha>0, beta>0
!
!x = phi/Mp

module ssbi1sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory, ssbi136_x_epstwozero


  implicit none

  private

  public ssbi1_norm_potential
  public ssbi1_epsilon_one, ssbi1_epsilon_two, ssbi1_epsilon_three
  public ssbi1_x_endinf, ssbi1_efold_primitive, ssbi1_x_trajectory
  public ssbi1_norm_deriv_potential, ssbi1_norm_deriv_second_potential
  public ssbi1_alphamin, ssbi1_x_epsonemax
  
contains

!returns V/M**4
  function ssbi1_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi1_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi1_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi1_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi1_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi1_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi1_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi1_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi1_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi1_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi1_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi1_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi1_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi1_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi1_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi1_epsilon_one


!epsilon_two(x)
  function ssbi1_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi1_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi1_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi1_epsilon_two


!epsilon_three(x)
  function ssbi1_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi1_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi1_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi1_epsilon_three


!returns the position x where epsilon_one is maximum
  function ssbi1_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssbi1_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssbi1_x_epsonemax = ssbi136_x_epstwozero(alpha,beta)
    
  end function ssbi1_x_epsonemax


!returns the minimum value of alpha (given beta) in order for
!inflation to end by slow roll violation (eps1max>1)
  function ssbi1_alphamin(beta)    
    implicit none
    real(kp) :: ssbi1_alphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi1Data


    mini=1._kp
    maxi=10._kp**(3._kp)

    ssbi1Data%real1 = beta


!In that case inflation ends by slow roll violation for any value of alpha
    if(beta .gt. 0.251_kp) then

       ssbi1_alphamin = 0._kp
 
    else

       ssbi1_alphamin = zbrent(find_ssbi1_alphamin,mini,maxi,tolFind,ssbi1Data)

    endif
    
  end function ssbi1_alphamin

  function find_ssbi1_alphamin(alpha,ssbi1Data)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: ssbi1Data
    real(kp) :: find_ssbi1_alphamin
    real(kp) :: beta

    beta = ssbi1Data%real1
    
    find_ssbi1_alphamin = ssbi1_epsilon_one(ssbi1_x_epsonemax(alpha,beta),alpha,beta)-1._kp
   
  end function find_ssbi1_alphamin


!returns x at the end of inflation defined as epsilon1=1
  function ssbi1_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi1_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, eps1max
    type(transfert) :: ssbi1Data

    eps1max = ssbi1_epsilon_one(ssbi1_x_epsonemax(alpha,beta),alpha,beta)

    if (eps1max.lt.1._kp) then
       write(*,*) 'eps1max=',eps1max
       stop 'ssbi1_x_endinf: epsilon1max< 1, inflation cannot stop by slow roll violation!'
    endif

    mini = ssbi1_x_epsonemax(alpha,beta)
    maxi = mini/epsilon(1._kp)

!    print*,'ss1_xend:  mini=',mini,'   maxi=',maxi,'   epsOne(mini)=',ssbi1_epsilon_one(mini,alpha,beta), &
!                 '   epsOne(maxi)=',ssbi1_epsilon_one(maxi,alpha,beta)
!    pause

    ssbi1Data%real1 = alpha
    ssbi1Data%real2 = beta
    
    ssbi1_x_endinf = zbrent(find_ssbi1_x_endinf,mini,maxi,tolFind,ssbi1Data)

  end function ssbi1_x_endinf



  function find_ssbi1_x_endinf(x,ssbi1Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi1Data
    real(kp) :: find_ssbi1_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi1Data%real1
    beta = ssbi1Data%real2
    
    find_ssbi1_x_endinf = ssbi1_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi1_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi1_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi1_efold_primitive

    ssbi1_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi1_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi1_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi1_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi1Data


    mini = ssbi1_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)
  
    ssbi1Data%real1 = alpha
    ssbi1Data%real2 = beta
    ssbi1Data%real3 = -bfold + ssbi1_efold_primitive(xend,alpha,beta)
    
    ssbi1_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi1Data)
       
  end function ssbi1_x_trajectory
 

end module ssbi1sr

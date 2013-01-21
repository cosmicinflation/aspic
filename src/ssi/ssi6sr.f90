!slow-roll functions for the sneutrino supersymmetric 6 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!5: alpha<0, beta>0, x**2 > -alpha / ( 2 beta )
!
!x = phi/Mp

module ssi6sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential
  use ssicommon, only : ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssi_x_trajectory 
  use ssicommon, only : ssi3456_x_derivpotzero,ssi136_x_epstwozero


  implicit none

  private

  public ssi6_norm_potential
  public ssi6_epsilon_one, ssi6_epsilon_two, ssi6_epsilon_three
  public ssi6_x_endinf, ssi6_efold_primitive, ssi6_x_trajectory
  public ssi6_norm_deriv_potential, ssi6_norm_deriv_second_potential
  public ssi6_alphamax, ssi6_x_epsonemax, ssi6_x_potzero
  
contains

!returns V/M**4
  function ssi6_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi6_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi6_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi6_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi6_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi6_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi6_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi6_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi6_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi6_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi6_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi6_norm_deriv_second_potential



!epsilon_one(x)
  function ssi6_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi6_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi6_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi6_epsilon_one


!epsilon_two(x)
  function ssi6_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi6_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi6_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi6_epsilon_two


!epsilon_three(x)
  function ssi6_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi6_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi6_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi6_epsilon_three


!returns the position x where epsilon_one is maximum when alpha**2 < 4 * beta
  function ssi6_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssi6_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssi6_x_epsonemax = ssi136_x_epstwozero(alpha,beta)
    
  end function ssi6_x_epsonemax



 function ssi6_alphamax(beta)
    implicit none
    real(kp) :: ssi6_alphamax
    real(kp), intent(in) :: beta
    
    ssi6_alphamax = -ssi6_absalphamin(beta)

  end function ssi6_alphamax


!returns the minimum value of the absolute value of alpha (given beta) in order for inflation to end by slow roll violation (eps1max>1)
  function ssi6_absalphamin(beta)    
    implicit none
    real(kp) :: ssi6_absalphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi6Data

    mini=epsilon(1._kp)
    maxi=2._kp*sqrt(beta)*(1._kp-1000._kp*epsilon(1._kp))

    ssi6Data%real1 = beta


     if (ssi6_epsilon_one(ssi6_x_epsonemax(-mini,beta),-mini,beta) .gt. 1._kp ) then !In that case inflation ends by slow roll violation for any value of alpha
        ssi6_absalphamin = 0._kp
     else    
        ssi6_absalphamin = zbrent(find_ssi6_absalphamin,mini,maxi,tolFind,ssi6Data)
     endif

!    print*,'ssi6_alphamin: ssi6_absalphamin=',ssi6_absalphamin
!    pause
    
    
  end function ssi6_absalphamin

  function find_ssi6_absalphamin(abs_alpha,ssi6Data)    
    implicit none
    real(kp), intent(in) :: abs_alpha   
    type(transfert), optional, intent(inout) :: ssi6Data
    real(kp) :: find_ssi6_absalphamin
    real(kp) :: beta

    beta = ssi6Data%real1
    
    find_ssi6_absalphamin = ssi6_epsilon_one(ssi6_x_epsonemax(-abs_alpha,beta) &
         ,-abs_alpha,beta)-1._kp
   
  end function find_ssi6_absalphamin


! Return the position x at which the potential vanishes
  function ssi6_x_potzero(alpha,beta)
  implicit none
    real(kp) :: ssi6_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssi6_x_potzero = sqrt(-(alpha-sqrt(alpha**2-4._kp*beta))/(2._kp*beta))

   end function ssi6_x_potzero

!returns x at the end of inflation defined as epsilon1=1
  function ssi6_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi6_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi6Data

    if (abs(alpha) .lt. ssi6_absalphamin(beta)) then
       stop 'ssi6_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'
    endif


    if (alpha**2 .lt. 4._kp*beta) then
       	mini=ssi6_x_epsonemax(alpha,beta)*(1._kp+epsilon(1._kp))
!        print*,'ssi6_xend: maxi in the first case=',maxi
    else
	mini=ssi6_x_potzero(alpha,beta)*(1._kp+epsilon(1._kp))
!        print*,'ssi6_xend: maxi in the second case=',maxi
    endif
    maxi=mini/epsilon(1._kp)

    ssi6Data%real1 = alpha
    ssi6Data%real2 = beta
    
    ssi6_x_endinf = zbrent(find_ssi6_x_endinf,mini,maxi,tolFind,ssi6Data)

  end function ssi6_x_endinf



  function find_ssi6_x_endinf(x,ssi6Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi6Data
    real(kp) :: find_ssi6_x_endinf
    real(kp) :: alpha,beta

    alpha = ssi6Data%real1
    beta = ssi6Data%real2
    
    find_ssi6_x_endinf = ssi6_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi6_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi6_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi6_efold_primitive

    ssi6_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi6_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi6_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi6_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi6Data


    mini = ssi6_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = 10._kp**(6._kp)*mini
  
    ssi6Data%real1 = alpha
    ssi6Data%real2 = beta
    ssi6Data%real3 = -bfold + ssi6_efold_primitive(xend,alpha,beta)
    
    ssi6_x_trajectory = zbrent(find_ssi_x_trajectory,mini,maxi,tolFind,ssi6Data)
       
  end function ssi6_x_trajectory
 

end module ssi6sr

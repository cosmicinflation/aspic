!slow-roll functions for the sspontaneous symmetry breaking 6 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!5: alpha<0, beta>0, x**2 > -alpha / ( 2 beta )
!
!x = phi/Mp

module ssbi6sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory 
  use ssbicommon, only : ssbi3456_x_derivpotzero,ssbi136_x_epstwozero


  implicit none

  private

  public ssbi6_norm_potential
  public ssbi6_epsilon_one, ssbi6_epsilon_two, ssbi6_epsilon_three
  public ssbi6_x_endinf, ssbi6_efold_primitive, ssbi6_x_trajectory
  public ssbi6_norm_deriv_potential, ssbi6_norm_deriv_second_potential
  public ssbi6_alphamax, ssbi6_x_epsonemax, ssbi6_x_potzero
  
contains

!returns V/M**4
  function ssbi6_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi6_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi6_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi6_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi6_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi6_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi6_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi6_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi6_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi6_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi6_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi6_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi6_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi6_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi6_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi6_epsilon_one


!epsilon_two(x)
  function ssbi6_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi6_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi6_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi6_epsilon_two


!epsilon_three(x)
  function ssbi6_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi6_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi6_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi6_epsilon_three


!returns the position x where epsilon_one is maximum when alpha**2 < 4 * beta
  function ssbi6_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssbi6_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    ssbi6_x_epsonemax = ssbi136_x_epstwozero(alpha,beta)
    
  end function ssbi6_x_epsonemax



 function ssbi6_alphamax(beta)
    implicit none
    real(kp) :: ssbi6_alphamax
    real(kp), intent(in) :: beta
    
    ssbi6_alphamax = -ssbi6_absalphamin(beta)

  end function ssbi6_alphamax


!returns the minimum value of the absolute value of alpha (given beta) in order for inflation to end by slow roll violation (eps1max>1)
  function ssbi6_absalphamin(beta)    
    implicit none
    real(kp) :: ssbi6_absalphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi6Data

    mini=epsilon(1._kp)
    maxi=2._kp*sqrt(beta)*(1._kp-1000._kp*epsilon(1._kp))

    ssbi6Data%real1 = beta


     if (ssbi6_epsilon_one(ssbi6_x_epsonemax(-mini,beta),-mini,beta) .gt. 1._kp ) then !In that case inflation ends by slow roll violation for any value of alpha
        ssbi6_absalphamin = 0._kp
     else    
        ssbi6_absalphamin = zbrent(find_ssbi6_absalphamin,mini,maxi,tolFind,ssbi6Data)
     endif

!    print*,'ssbi6_alphamin: ssbi6_absalphamin=',ssbi6_absalphamin
!    pause
    
    
  end function ssbi6_absalphamin

  function find_ssbi6_absalphamin(abs_alpha,ssbi6Data)    
    implicit none
    real(kp), intent(in) :: abs_alpha   
    type(transfert), optional, intent(inout) :: ssbi6Data
    real(kp) :: find_ssbi6_absalphamin
    real(kp) :: beta

    beta = ssbi6Data%real1
    
    find_ssbi6_absalphamin = ssbi6_epsilon_one(ssbi6_x_epsonemax(-abs_alpha,beta) &
         ,-abs_alpha,beta)-1._kp
   
  end function find_ssbi6_absalphamin


! Return the position x at which the potential vanishes
  function ssbi6_x_potzero(alpha,beta)
  implicit none
    real(kp) :: ssbi6_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssbi6_x_potzero = sqrt(-(alpha-sqrt(alpha**2-4._kp*beta))/(2._kp*beta))

   end function ssbi6_x_potzero

!returns x at the end of inflation defined as epsilon1=1
  function ssbi6_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi6_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi6Data

    if (abs(alpha) .lt. ssbi6_absalphamin(beta)) then
       stop 'ssbi6_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'
    endif


    if (alpha**2 .lt. 4._kp*beta) then
       	mini=ssbi6_x_epsonemax(alpha,beta)*(1._kp+epsilon(1._kp))
!        print*,'ssbi6_xend: maxi in the first case=',maxi
    else
	mini=ssbi6_x_potzero(alpha,beta)*(1._kp+epsilon(1._kp))
!        print*,'ssbi6_xend: maxi in the second case=',maxi
    endif
    maxi=mini/epsilon(1._kp)

    ssbi6Data%real1 = alpha
    ssbi6Data%real2 = beta
    
    ssbi6_x_endinf = zbrent(find_ssbi6_x_endinf,mini,maxi,tolFind,ssbi6Data)

  end function ssbi6_x_endinf



  function find_ssbi6_x_endinf(x,ssbi6Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi6Data
    real(kp) :: find_ssbi6_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi6Data%real1
    beta = ssbi6Data%real2
    
    find_ssbi6_x_endinf = ssbi6_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi6_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi6_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi6_efold_primitive

    ssbi6_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi6_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi6_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi6_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi6Data


    mini = ssbi6_x_endinf(alpha,beta)*(1._kp+epsilon(1._kp))
    maxi = 10._kp**(6._kp)*mini
  
    ssbi6Data%real1 = alpha
    ssbi6Data%real2 = beta
    ssbi6Data%real3 = -bfold + ssbi6_efold_primitive(xend,alpha,beta)
    
    ssbi6_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi6Data)
       
  end function ssbi6_x_trajectory
 

end module ssbi6sr

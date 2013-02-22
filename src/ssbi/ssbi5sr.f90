!slow-roll functions for the spontaneous symmetry breaking 5 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!5: alpha<0, beta>0, x**2 < -alpha / ( 2 beta )
!
!x = phi/Mp

module ssbi5sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssbicommon, only : ssbi_norm_potential, ssbi_norm_deriv_potential
  use ssbicommon, only : ssbi_norm_deriv_second_potential
  use ssbicommon, only : ssbi_epsilon_one, ssbi_epsilon_two, ssbi_epsilon_three
  use ssbicommon, only : ssbi_efold_primitive, find_ssbi_x_trajectory 
  use ssbicommon, only : ssbi245_x_potzero


  implicit none

  private

  public ssbi5_norm_potential
  public ssbi5_epsilon_one, ssbi5_epsilon_two, ssbi5_epsilon_three
  public ssbi5_x_endinf, ssbi5_efold_primitive, ssbi5_x_trajectory
  public ssbi5_norm_deriv_potential, ssbi5_norm_deriv_second_potential
  public ssbi5_x_potzero, ssbi5_alphamax
  
contains

!returns V/M**4
  function ssbi5_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi5_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssbi5_norm_potential = ssbi_norm_potential(x,alpha,beta)

  end function ssbi5_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssbi5_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi5_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssbi5_norm_deriv_potential = ssbi_norm_deriv_potential(x,alpha,beta)

  end function ssbi5_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssbi5_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssbi5_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssbi5_norm_deriv_second_potential &
         = ssbi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssbi5_norm_deriv_second_potential



!epsilon_one(x)
  function ssbi5_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi5_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssbi5_epsilon_one = ssbi_epsilon_one(x,alpha,beta)
    
  end function ssbi5_epsilon_one


!epsilon_two(x)
  function ssbi5_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi5_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssbi5_epsilon_two = ssbi_epsilon_two(x,alpha,beta)
    
  end function ssbi5_epsilon_two


!epsilon_three(x)
  function ssbi5_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssbi5_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssbi5_epsilon_three = ssbi_epsilon_three(x,alpha,beta)
    
  end function ssbi5_epsilon_three


!returns the position x where epsilon_one is maximum when alpha**2 < 4 * beta
  function ssbi5_x_epsonemax(alpha,beta)    
    implicit none
    real(kp) :: ssbi5_x_epsonemax
    real(kp), intent(in) :: alpha,beta

    if (alpha**2 .gt. 4._kp*beta) then
       write(*,*)'(alpha=', alpha,'  beta=',beta,')'
       stop 'ssbi5_x_epsonemax: error: alpha^2 > 4 beta, epsilon_2 never vanishes'
    else

    ssbi5_x_epsonemax = sqrt(-(alpha/(6._kp*beta))+((1._kp-complex(0._kp,1._kp)*sqrt(3._kp))* &
         (5._kp*alpha**2*beta**2-36._kp*beta**3))/(6._kp*2._kp**(2._kp/3._kp)* &
         beta**2*(16._kp*alpha**3*beta**3+sqrt(complex(256._kp*alpha**6*beta**6+ &
         4._kp*(5._kp*alpha**2*beta**2-36._kp*beta**3)**3,0._kp)))**(1._kp/3._kp))- &
         ((1._kp+complex(0._kp,1._kp)*sqrt(3._kp))*(16._kp*alpha**3*beta**3+ &
         sqrt(complex(256._kp*alpha**6*beta**6+4._kp*(5._kp*alpha**2*beta**2- &
         36._kp*beta**3)**3,0._kp)))**(1._kp/3._kp))/(12._kp*2._kp**(1._kp/3._kp)*beta**2))

    endif
    
  end function ssbi5_x_epsonemax


  function ssbi5_x_potzero(alpha,beta)    
    implicit none
    real(kp) :: ssbi5_x_potzero
    real(kp), intent(in) :: alpha,beta

    ssbi5_x_potzero = ssbi245_x_potzero(alpha,beta)
    
  end function ssbi5_x_potzero


  function ssbi5_alphamax(beta)
    implicit none
    real(kp) :: ssbi5_alphamax
    real(kp), intent(in) :: beta
    
    ssbi5_alphamax = -ssbi5_absalphamin(beta)

  end function ssbi5_alphamax



!returns the minimum value of the absolute value of alpha (given beta)
!in order for inflation to end by slow roll violation (eps1max>1)
  function ssbi5_absalphamin(beta)    
    implicit none
    real(kp) :: ssbi5_absalphamin
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi5Data

    mini=epsilon(1._kp)
    maxi=2._kp*sqrt(beta)*(1._kp-100._kp*epsilon(1._kp))

    ssbi5Data%real1 = beta

    ssbi5_absalphamin = zbrent(find_ssbi5_absalphamin,mini,maxi,tolFind,ssbi5Data)

    
  end function ssbi5_absalphamin

  function find_ssbi5_absalphamin(abs_alpha,ssbi5Data)    
    implicit none
    real(kp), intent(in) :: abs_alpha   
    type(transfert), optional, intent(inout) :: ssbi5Data
    real(kp) :: find_ssbi5_absalphamin
    real(kp) :: beta

    beta = ssbi5Data%real1
    
    find_ssbi5_absalphamin = ssbi5_epsilon_one(ssbi5_x_epsonemax(-abs_alpha,beta) &
         ,- abs_alpha,beta) - 1._kp
   
  end function find_ssbi5_absalphamin


!returns x at the end of inflation defined as epsilon1=1
  function ssbi5_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssbi5_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi5Data

    if (alpha .gt. ssbi5_alphamax(beta)) then
       stop 'ssbi5_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'
    endif

    mini=epsilon(1._kp)
    if (alpha**2 .lt. 4._kp*beta) then
       	maxi=ssbi5_x_epsonemax(alpha,beta)*(1._kp-epsilon(1._kp))
    else
	maxi=ssbi5_x_potzero(alpha,beta)*(1._kp-epsilon(1._kp))
    endif

    ssbi5Data%real1 = alpha
    ssbi5Data%real2 = beta
    
    ssbi5_x_endinf = zbrent(find_ssbi5_x_endinf,mini,maxi,tolFind,ssbi5Data)

  end function ssbi5_x_endinf



  function find_ssbi5_x_endinf(x,ssbi5Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssbi5Data
    real(kp) :: find_ssbi5_x_endinf
    real(kp) :: alpha,beta

    alpha = ssbi5Data%real1
    beta = ssbi5Data%real2
    
    find_ssbi5_x_endinf = ssbi5_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssbi5_x_endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssbi5_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssbi5_efold_primitive

    ssbi5_efold_primitive = ssbi_efold_primitive(x,alpha,beta)
  
  end function ssbi5_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssbi5_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssbi5_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssbi5Data


    mini = epsilon(1._kp)
    maxi = ssbi5_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssbi5Data%real1 = alpha
    ssbi5Data%real2 = beta
    ssbi5Data%real3 = -bfold + ssbi5_efold_primitive(xend,alpha,beta)
    
    ssbi5_x_trajectory = zbrent(find_ssbi_x_trajectory,mini,maxi,tolFind,ssbi5Data)
       
  end function ssbi5_x_trajectory
 

end module ssbi5sr

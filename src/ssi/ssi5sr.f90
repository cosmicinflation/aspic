!slow-roll functions for the sneutrino supersymmetric 5 potential
!
!
!V(phi) = M**4 [ 1 + alpha x**2 + beta x**4 ]
!
!5: alpha<0, beta>0, x**2 < -alpha / ( 2 beta )
!
!x = phi/Mp

module ssi5sr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use ssicommon, only : ssi_norm_potential, ssi_norm_deriv_potential
  use ssicommon, only : ssi_norm_deriv_second_potential
  use ssicommon, only : ssi_epsilon_one, ssi_epsilon_two, ssi_epsilon_three
  use ssicommon, only : ssi_efold_primitive, find_ssitraj 
  use ssicommon, only : ssi3456_x_Vprime_Equals_0, ssi245_x_V_Equals_0


  implicit none

  private

  public ssi5_norm_potential
  public ssi5_epsilon_one, ssi5_epsilon_two, ssi5_epsilon_three
  public ssi5_x_endinf, ssi5_efold_primitive, ssi5_x_trajectory
  public ssi5_norm_deriv_potential, ssi5_norm_deriv_second_potential
  public ssi5_abs_alpha_min
  
contains

!returns V/M**4
  function ssi5_norm_potential(x,alpha,beta)    
    implicit none
    real(kp) :: ssi5_norm_potential
    real(kp), intent(in) :: x,alpha,beta

    ssi5_norm_potential = ssi_norm_potential(x,alpha,beta)

  end function ssi5_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M**4
  function ssi5_norm_deriv_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi5_norm_deriv_potential
    real(kp), intent(in) :: x,alpha,beta
    
    ssi5_norm_deriv_potential = ssi_norm_deriv_potential(x,alpha,beta)

  end function ssi5_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M**4
  function ssi5_norm_deriv_second_potential(x,alpha,beta)
    implicit none
    real(kp) :: ssi5_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha,beta    

    ssi5_norm_deriv_second_potential &
         = ssi_norm_deriv_second_potential(x,alpha,beta)
    

  end function ssi5_norm_deriv_second_potential



!epsilon_one(x)
  function ssi5_epsilon_one(x,alpha,beta)    
    implicit none
    real(kp) :: ssi5_epsilon_one
    real(kp), intent(in) :: x,alpha,beta

    ssi5_epsilon_one = ssi_epsilon_one(x,alpha,beta)
    
  end function ssi5_epsilon_one


!epsilon_two(x)
  function ssi5_epsilon_two(x,alpha,beta)    
    implicit none
    real(kp) :: ssi5_epsilon_two
    real(kp), intent(in) :: x,alpha,beta

    ssi5_epsilon_two = ssi_epsilon_two(x,alpha,beta)
    
  end function ssi5_epsilon_two


!epsilon_three(x)
  function ssi5_epsilon_three(x,alpha,beta)    
    implicit none
    real(kp) :: ssi5_epsilon_three
    real(kp), intent(in) :: x,alpha,beta
   
    ssi5_epsilon_three = ssi_epsilon_three(x,alpha,beta)
    
  end function ssi5_epsilon_three


!returns the position x where epsilon_one is maximum when alpha**2 < 4 * beta
  function ssi5_x_epsilon2_Equals_0(alpha,beta)    
    implicit none
    real(kp) :: ssi5_x_epsilon2_Equals_0
    real(kp), intent(in) :: alpha,beta

    if (alpha**2 .gt. 4._kp*beta) then
               print*,'(alpha=', alpha,'  beta=',beta,')'
               stop 'ssi5_x_epsilon2_Equals_0: error: alpha^2 > 4 beta, epsilon_2 never vanishes'
    else

    ssi5_x_epsilon2_Equals_0 = sqrt(-(alpha/(6._kp*beta))+((1._kp-complex(0._kp,1._kp)*sqrt(3._kp))* &
                               (5._kp*alpha**2*beta**2-36._kp*beta**3))/(6._kp*2._kp**(2._kp/3._kp)* &
                               beta**2*(16._kp*alpha**3*beta**3+sqrt(complex(256._kp*alpha**6*beta**6+ &
                               4._kp*(5._kp*alpha**2*beta**2-36._kp*beta**3)**3,0._kp)))**(1._kp/3._kp))- &
                               ((1._kp+complex(0._kp,1._kp)*sqrt(3._kp))*(16._kp*alpha**3*beta**3+ &
                               sqrt(complex(256._kp*alpha**6*beta**6+4._kp*(5._kp*alpha**2*beta**2- &
                               36._kp*beta**3)**3,0._kp)))**(1._kp/3._kp))/(12._kp*2._kp**(1._kp/3._kp)*beta**2))

    endif
    
  end function ssi5_x_epsilon2_Equals_0


!returns the minimum value of the absolute value of alpha (given beta) in order for inflation to end by slow roll violation (eps1max>1)
  function ssi5_abs_alpha_min(beta)    
    implicit none
    real(kp) :: ssi5_abs_alpha_min
    real(kp), intent(in) :: beta
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi5Data

    mini=epsilon(1._kp)
    maxi=2._kp*sqrt(beta)*(1._kp-100._kp*epsilon(1._kp))

    ssi5Data%real1 = beta

!    print*,'ssi5_alphamin:  beta=',beta,'mini=',mini,'  xeps2mini=',ssi5_x_epsilon2_Equals_0(-mini,beta), &
!           'maxi=',maxi,'  xeps2maxi=',ssi5_x_epsilon2_Equals_0(-maxi,beta), &
!           '  eps1(xeps2mini)=',ssi5_epsilon_one(ssi5_x_epsilon2_Equals_0(-mini,beta),-mini,beta), &
!           '  eps1(xeps2maxi)=',ssi5_epsilon_one(ssi5_x_epsilon2_Equals_0(-maxi,beta),-maxi,beta)
!    pause

       ssi5_abs_alpha_min = zbrent(find_ssi5_abs_alpha_min,mini,maxi,tolFind,ssi5Data)

!    print*,'ssi5_alphamin: ssi5_abs_alpha_min=',ssi5_abs_alpha_min
!    pause
    
    
  end function ssi5_abs_alpha_min

  function find_ssi5_abs_alpha_min(abs_alpha,ssi5Data)    
    implicit none
    real(kp), intent(in) :: abs_alpha   
    type(transfert), optional, intent(inout) :: ssi5Data
    real(kp) :: find_ssi5_abs_alpha_min
    real(kp) :: beta

    beta = ssi5Data%real1
    
    find_ssi5_abs_alpha_min = ssi5_epsilon_one(ssi5_x_epsilon2_Equals_0(-abs_alpha,beta),-abs_alpha,beta)-1._kp
   
  end function find_ssi5_abs_alpha_min


!returns x at the end of inflation defined as epsilon1=1
  function ssi5_x_endinf(alpha,beta)
    implicit none
    real(kp), intent(in) :: alpha,beta
    real(kp) :: ssi5_x_endinf
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi5Data

    if (abs(alpha) .lt. ssi5_abs_alpha_min(beta)) stop 'ssi5_x_endinf: epsilon1max<1, inflation cannot stop by slow roll violation!'

    mini=epsilon(1._kp)
    if (alpha**2 .lt. 4._kp*beta) then
       	maxi=ssi5_x_epsilon2_Equals_0(alpha,beta)*(1._kp-epsilon(1._kp))
    else
	maxi=ssi245_x_V_Equals_0(alpha,beta)*(1._kp-epsilon(1._kp))
    endif


!    print*,'ssi5_xend:  mini=',mini,'   maxi=',maxi,'   epsOne(mini)=',ssi5_epsilon_one(mini,alpha,beta), &
!                 '   epsOne(maxi)=',ssi5_epsilon_one(maxi,alpha,beta)
!    pause

    ssi5Data%real1 = alpha
    ssi5Data%real2 = beta
    
    ssi5_x_endinf = zbrent(find_ssi5endinf,mini,maxi,tolFind,ssi5Data)

  end function ssi5_x_endinf



  function find_ssi5endinf(x,ssi5Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: ssi5Data
    real(kp) :: find_ssi5endinf
    real(kp) :: alpha,beta

    alpha = ssi5Data%real1
    beta = ssi5Data%real2
    
    find_ssi5endinf = ssi5_epsilon_one(x,alpha,beta) - 1._kp
   
  end function find_ssi5endinf

!this is integral[V(phi)/V'(phi) dphi]
  function ssi5_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: ssi5_efold_primitive

    ssi5_efold_primitive = ssi_efold_primitive(x,alpha,beta)
  
  end function ssi5_efold_primitive


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function ssi5_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend,beta
    real(kp) :: ssi5_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: ssi5Data


    mini = epsilon(1._kp)
    maxi = ssi5_x_endinf(alpha,beta)*(1._kp-epsilon(1._kp))
  
    ssi5Data%real1 = alpha
    ssi5Data%real2 = beta
    ssi5Data%real3 = -bfold + ssi5_efold_primitive(xend,alpha,beta)
    
    ssi5_x_trajectory = zbrent(find_ssitraj,mini,maxi,tolFind,ssi5Data)
       
  end function ssi5_x_trajectory
 

end module ssi5sr

!slow-roll functions for the GMSSM inflation potential
!
!V(phi) = M^4 [ x^2 - alpha x^6 + beta x^10 ]
!
!x = phi/Mp

module gmssmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : atan_ito_log
  use gmssmicommon, only : gmssmi_norm_potential, gmssmi_norm_deriv_potential, gmssmi_norm_deriv_second_potential
  use gmssmicommon, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmicommon, only : gmssmi_x_endinf, gmssmi_x_epsonemin
  implicit none

  private

  public  gmssmi_norm_potential, gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  public  gmssmi_x_endinf, gmssmi_efold_primitive, gmssmi_x_trajectory
  public  gmssmi_norm_deriv_potential, gmssmi_norm_deriv_second_potential
  public  gmssmi_x_epsonemin,gmssmi_alphamin

 
contains


!this is integral[V(phi)/V'(phi) dphi]
  function gmssmi_efold_primitive(x,alpha,beta)
    implicit none
    real(kp), intent(in) :: x,alpha,beta
    real(kp) :: gmssmi_efold_primitive
    complex(kp) ::aplus,aminus,bplus,bminus

    if (beta.eq.0._kp.and.alpha.eq.0._kp) stop 'gmssmi_efold_primitive: alpha=0 and beta=0!'

    if (alpha/beta.lt.10._kp**(-5.).and.beta.lt.10._kp**(-10.)) then
    ! This is an approximated formula, valid when alpha<<1 , beta<<1 , and which leads to better numerical convergence
       gmssmi_efold_primitive = real(x**2/4._kp)
    else

    aplus=-1.5_kp*alpha+0.5_kp*sqrt(complex(9._kp*alpha**2-20._kp*beta,0._kp))
    aminus=-1.5_kp*alpha-0.5_kp*sqrt(complex(9._kp*alpha**2-20._kp*beta,0._kp))
    bplus=(2._kp*aplus+alpha)/(aplus-aminus)
    bminus=(2._kp*aminus+alpha)/(aminus-aplus)

    gmssmi_efold_primitive = real(x**2/20._kp+bplus/(10._kp*sqrt(aplus))*atan_ito_log(sqrt(aplus)*x**2) &
         +bminus/(10._kp*sqrt(aminus))*atan_ito_log(sqrt(aminus)*x**2))

    endif

  end function gmssmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gmssmi_x_trajectory(bfold,xend,alpha,beta)
    implicit none
    real(kp), intent(in) :: bfold, alpha, beta, xend
    real(kp) :: gmssmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

  
    mini = xend

    if (alpha**2/beta>20._kp/9._kp) then
	maxi = gmssmi_x_epsonemin(alpha,beta) !local maximum of the potential
    else
	maxi=100._kp*gmssmi_x_epsonemin(alpha,beta)
    endif

    gmssmiData%real1 = alpha
    gmssmiData%real2 = beta
    gmssmiData%real3 = -bfold + gmssmi_efold_primitive(xend,alpha,beta)
    
    gmssmi_x_trajectory = zbrent(find_gmssmi_x_trajectory,mini,maxi,tolFind,gmssmiData)
       
  end function gmssmi_x_trajectory

  function find_gmssmi_x_trajectory(x,gmssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_x_trajectory
    real(kp) :: alpha,beta,NplusNuend

    alpha= gmssmiData%real1
    beta = gmssmiData%real2
    NplusNuend = gmssmiData%real3

    find_gmssmi_x_trajectory = gmssmi_efold_primitive(x,alpha,beta) - NplusNuend
   
  end function find_gmssmi_x_trajectory


!Returns the prior alphamin(beta) such that the first local minimum of epsilon1 is less than one
  function gmssmi_alphamin(beta)    
    implicit none
    real(kp), intent(in) :: beta   
    real(kp) :: gmssmi_alphamin
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

    if (beta<0.000796037) then 
	gmssmi_alphamin = 0._kp
    else

    mini = 0._kp
    maxi = 10._kp
  

    gmssmiData%real1 = beta

    gmssmi_alphamin = zbrent(find_gmssmi_alphamin,mini,maxi,tolFind,gmssmiData)

    endif
   
  end function gmssmi_alphamin

  function find_gmssmi_alphamin(alpha,gmssmiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_alphamin
    real(kp) :: beta

    beta= gmssmiData%real1

    find_gmssmi_alphamin = gmssmi_epsilon_one(gmssmi_x_epsonemin(alpha,beta),alpha,beta)-1._kp
   
  end function find_gmssmi_alphamin

  
end module gmssmisr

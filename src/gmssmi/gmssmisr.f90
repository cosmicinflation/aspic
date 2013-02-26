!slow-roll functions for the GMSSM inflation potential
!
!V(phi) = M^4 [ x^2 - 2/3 alpha x^6 + alpha/5 x^10 ]
!
!x = phi/phi0

module gmssmisr
  use infprec, only : kp,tolkp,transfert
  use inftools, only : zbrent
  use specialinf, only : atan_ito_log
  use gmssmicommon, only : gmssmi_norm_potential, gmssmi_norm_deriv_potential
  use gmssmicommon, only : gmssmi_norm_deriv_second_potential
  use gmssmicommon, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmicommon, only : gmssmi_x_endinf, gmssmi_x_epsonemin
  implicit none

  private

  public  gmssmi_norm_potential, gmssmi_epsilon_one, gmssmi_epsilon_two
  public  gmssmi_epsilon_three,gmssmi_x_endinf
  public  gmssmi_efold_primitive, gmssmi_x_trajectory
  public  gmssmi_norm_deriv_potential, gmssmi_norm_deriv_second_potential
  public  gmssmi_alphamin, gmssmi_alphamax, gmssmi_x_epsonemin

 
contains


!this is integral[V(phi)/V'(phi) dphi]
  function gmssmi_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: gmssmi_efold_primitive
    real(kp) :: xVprime,a,b
    complex(kp) ::aplus,aminus,bplus,bminus,test1,test2

    if (alpha.eq.0._kp) then 
       gmssmi_efold_primitive = x**2/4._kp

    elseif (alpha.eq.1._kp) then

       stop 'gmssmi_efold_primitive: you must use mssmi for alpha=1!'

    else
       
       aplus=-alpha+sqrt((alpha**2-alpha)*(1._kp,0._kp))
       aminus=-alpha-sqrt((alpha**2-alpha)*(1._kp,0._kp))
       bplus=2._kp*(aplus+alpha/3._kp)/(aplus-aminus)
       bminus=2._kp*(aminus+alpha/3._kp)/(aminus-aplus)

       gmssmi_efold_primitive = phi0**2*(real(x**2/20._kp &
            +bplus/(10._kp*sqrt(aplus))*atan_ito_log(sqrt(aplus)*x**2) &
            +bminus/(10._kp*sqrt(aminus))*atan_ito_log(sqrt(aminus)*x**2)))

    endif

  end function gmssmi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gmssmi_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold, alpha, phi0, xend
    real(kp) :: gmssmi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

  
    mini = xend

    if (alpha .lt. 1._kp) then

	maxi = gmssmi_x_epsonemin(alpha)*100._kp 

    elseif (alpha.eq.1._kp) then

       stop 'gmssmi_x_trajectory: you must use mssmi for alpha=1!'

    else

	maxi = gmssmi_x_epsonemin(alpha) !local maximum of the potential

    endif


    gmssmiData%real1 = alpha
    gmssmiData%real2 = phi0
    gmssmiData%real3 = -bfold + gmssmi_efold_primitive(xend,alpha,phi0)
    
    gmssmi_x_trajectory = zbrent(find_gmssmi_x_trajectory,mini,maxi,tolFind,gmssmiData)
       
  end function gmssmi_x_trajectory

  function find_gmssmi_x_trajectory(x,gmssmiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_x_trajectory
    real(kp) :: alpha,phi0,NplusNuend

    alpha= gmssmiData%real1
    phi0 = gmssmiData%real2
    NplusNuend = gmssmiData%real3

    find_gmssmi_x_trajectory = gmssmi_efold_primitive(x,alpha,phi0) - NplusNuend
   
  end function find_gmssmi_x_trajectory


!Returns the prior alphamin(phi0) such that, when alpha<1, the minimum of epsilon1 is less than one
  function gmssmi_alphamin(phi0)    
    implicit none
    real(kp), intent(in) :: phi0   
    real(kp) :: gmssmi_alphamin
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

    mini = epsilon(1._kp)
    mini = 0.1_kp
    maxi = 1._kp*(1._kp-epsilon(1._kp))
  
    gmssmiData%real1 = phi0

    gmssmi_alphamin = zbrent(find_gmssmi_alphamin,mini,maxi,tolFind,gmssmiData)

   
  end function gmssmi_alphamin

  function find_gmssmi_alphamin(alpha,gmssmiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_alphamin
    real(kp) :: phi0

    phi0= gmssmiData%real1

    find_gmssmi_alphamin = gmssmi_epsilon_one(gmssmi_x_epsonemin(alpha),alpha,phi0)-1._kp
   
  end function find_gmssmi_alphamin

!Returns the prior alphamax(phi0,DeltaNstar) such that, when alpha>1,
!at least efolds e-folds can be realized
  function gmssmi_alphamax(phi0,efold)    
    implicit none
    real(kp), intent(in) :: phi0,efold 
    real(kp) :: gmssmi_alphamax
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gmssmiData

    mini = 1._kp*(1._kp+epsilon(1._kp))
    mini = 1._kp*(1._kp+5._kp*epsilon(1._kp))
    maxi = 2._kp
  
    gmssmiData%real1 = phi0
    gmssmiData%real2 = efold

    gmssmi_alphamax = zbrent(find_gmssmi_alphamax,mini,maxi,tolFind,gmssmiData)

!    print*,'phi0=',phi0,'mini=',mini,'f(mini)=', &
!         gmssmi_efold_primitive(gmssmi_x_epsonemin(mini),mini,phi0) &
!                          -gmssmi_efold_primitive(gmssmi_x_endinf(mini,phi0),mini,phi0) &
!                          -efold, 'maxi=',maxi,'f(maxi)=',&
!                gmssmi_efold_primitive(gmssmi_x_epsonemin(maxi),maxi,phi0) &
!                          -gmssmi_efold_primitive(gmssmi_x_endinf(maxi,phi0),maxi,phi0) &
!                          -efold

   
  end function gmssmi_alphamax

  function find_gmssmi_alphamax(alpha,gmssmiData)    
    implicit none
    real(kp), intent(in) :: alpha   
    type(transfert), optional, intent(inout) :: gmssmiData
    real(kp) :: find_gmssmi_alphamax
    real(kp) :: phi0,efold

    phi0= gmssmiData%real1
    efold= gmssmiData%real2

    find_gmssmi_alphamax = gmssmi_efold_primitive(gmssmi_x_epsonemin(alpha),alpha,phi0) &
         -gmssmi_efold_primitive(gmssmi_x_endinf(alpha,phi0),alpha,phi0) &
         -efold
   
  end function find_gmssmi_alphamax

  
end module gmssmisr

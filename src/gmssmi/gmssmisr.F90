!slow-roll functions for the GMSSM inflation potential
!
!V(phi) = M^4 [ x^2 - 2/3 alpha x^6 + alpha/5 x^10 ]
!
!x = phi/phi0

module gmssmisr
  use infprec, only : kp,tolkp,pi,transfert
  use inftools, only : zbrent
#ifdef NOF08
  use specialinf, only : atan
#endif
  use gmssmicommon, only : gmssmi_norm_potential, gmssmi_norm_deriv_potential
  use gmssmicommon, only : gmssmi_norm_deriv_second_potential
  use gmssmicommon, only : gmssmi_epsilon_one, gmssmi_epsilon_two, gmssmi_epsilon_three
  use gmssmicommon, only : gmssmi_x_endinf, gmssmi_x_epsonemin
  use gmssmicommon, only : gmssmi_x_epstwomin, gmssmi_epstwomin
  implicit none

  private

  public gmssmi_norm_potential, gmssmi_epsilon_one, gmssmi_epsilon_two
  public gmssmi_epsilon_three,gmssmi_x_endinf
  public gmssmi_efold_primitive, gmssmi_x_trajectory
  public gmssmi_norm_deriv_potential, gmssmi_norm_deriv_second_potential
  public gmssmi_alphamin, gmssmi_x_epsonemin
  public gmssmi_epstwomin, gmssmi_x_epstwomin

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
            +bplus/(10._kp*sqrt(aplus))*atan(sqrt(aplus)*x**2) &
            +bminus/(10._kp*sqrt(aminus))*atan(sqrt(aminus)*x**2),kp))

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

       maxi = 1._kp/epsilon(1._kp)

    elseif (alpha.eq.1._kp) then

       stop 'gmssmi_x_trajectory: you must use mssmi for alpha=1!'

    else

	maxi = gmssmi_x_epsonemin(alpha,phi0) !local maximum of the potential

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


!for alpha < 1, returns alphamin to get at least efold number of inflation
  function gmssmi_alphamin(efold,phi0)
    implicit none
    real(kp) :: gmssmi_alphamin
    real(kp), intent(in) :: efold, phi0

    integer, save :: counter = 0
    integer, parameter :: repeat = 10
    
    real(kp) :: crazySmall

    crazySmall = phi0**4*pi**2/(900._kp*efold**2)

    if ((counter.le.repeat).and.(crazySmall.le.epsilon(1._kp))) then
       write(*,*)'gmssmi_alphamin: 1-alphamin < machine precision!'
       write(*,*)'1-alphamin= ',crazySmall
       counter = counter + 1
    endif

    gmssmi_alphamin = 1._kp - crazySmall 
    
  end function gmssmi_alphamin

  
end module gmssmisr

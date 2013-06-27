!slow-roll functions for the GRIP inflation potential
!
!V(phi) = M**4 [ x**2 - 4/3 alpha x**3 + alpha/2 x**4 ]
!
!x = phi/phi0

module gripisr
  use infprec, only : kp,pi,tolkp,transfert
  use inftools, only : zbrent
#ifdef NOF08
  use specialinf, only : atan
#endif
  use gripicommon, only : gripi_norm_potential, gripi_norm_deriv_potential
  use gripicommon, only : gripi_norm_deriv_second_potential
  use gripicommon, only : gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three
  use gripicommon, only : gripi_x_epsonemin, gripi_x_endinf
  use gripicommon, only : gripi_x_epstwozero,gripi_x_epsonezero
  use gripicommon, only : gripi_x_epstwomin, gripi_epstwomin
  implicit none

  private

  public gripi_norm_potential, gripi_norm_deriv_potential
  public gripi_norm_deriv_second_potential, gripi_alphamin
  public gripi_epsilon_one, gripi_epsilon_two, gripi_epsilon_three  
  public gripi_x_epsonemin, gripi_x_endinf
  public gripi_x_epstwozero,gripi_x_epsonezero
  public gripi_x_epstwomin, gripi_epstwomin
  public gripi_x_trajectory, gripi_efold_primitive


contains

 
!this is integral[V(phi)/V'(phi) dphi]
  function gripi_efold_primitive(x,alpha,phi0)
    implicit none
    real(kp), intent(in) :: x,alpha,phi0
    real(kp) :: gripi_efold_primitive
    complex(kp) :: carg

    if (alpha.eq.0._kp) then 
       gripi_efold_primitive = x**2/4._kp

    elseif (alpha.eq.1._kp) then

       stop 'gripi_efold_primitive: you must use ripi for alpha=1!'

    else

!       carg = (x-1._kp)/(sqrt(cmplx(1._kp/alpha-1._kp,0._kp,kp)))
       carg = (x-1._kp)/sqrt((1._kp/alpha-1._kp)*(1._kp,0._kp))

       gripi_efold_primitive = phi0**2*real((5._kp-4._kp)/ &
            (12._kp*sqrt((alpha*(1._kp-alpha))*(1._kp,0._kp)))* &
            atan(carg)+0.5_kp*x*(0.25_kp*x-1._kp/3._kp)+ &
            (1._kp/(8._kp*alpha)-1._kp/6._kp) &
            *log((1._kp+alpha*x*(x-2._kp))*(1._kp,0._kp) ),kp)

    endif

  end function gripi_efold_primitive



!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function gripi_x_trajectory(bfold,xend,alpha,phi0)
    implicit none
    real(kp), intent(in) :: bfold, alpha, phi0, xend
    real(kp) :: gripi_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: gripiData

    mini = xend

    if (alpha .lt. 1._kp) then

       maxi = 1._kp/epsilon(1._kp)

    elseif (alpha.eq.1._kp) then

       stop 'gripi_x_trajectory: you must use ripi for alpha=1!'

    else

       maxi = gripi_x_epsonemin(alpha,phi0) !local maximum of the potential

    endif

    gripiData%real1 = alpha
    gripiData%real2 = phi0
    gripiData%real3 = -bfold + gripi_efold_primitive(xend,alpha,phi0)

    gripi_x_trajectory = zbrent(find_gripi_x_trajectory,mini,maxi,tolFind,gripiData)

  end function gripi_x_trajectory

  function find_gripi_x_trajectory(x,gripiData)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: gripiData
    real(kp) :: find_gripi_x_trajectory
    real(kp) :: alpha,phi0,NplusNuend

    alpha= gripiData%real1
    phi0 = gripiData%real2
    NplusNuend = gripiData%real3

    find_gripi_x_trajectory = gripi_efold_primitive(x,alpha,phi0) - NplusNuend

  end function find_gripi_x_trajectory


!for alpha < 1, returns alphamin to get at least efold number of inflation
  function gripi_alphamin(efold,phi0)
    implicit none
    real(kp) :: gripi_alphamin
    real(kp), intent(in) :: efold, phi0

    integer, save :: counter = 0
    integer, parameter :: repeat = 10

    real(kp) :: crazySmall

    crazySmall = phi0**4*pi**2/(576._kp*efold**2)

    if ((counter.le.repeat).and.(crazySmall.le.epsilon(1._kp))) then
       write(*,*)'gripi_alphamin: 1-alphamin < machine precision!'
       write(*,*)'1-alphamin= ',crazySmall
       counter = counter + 1
    endif

    gripi_alphamin = 1._kp - crazySmall 
    
  end function gripi_alphamin

  

end module gripisr

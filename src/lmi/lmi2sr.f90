!slow-roll functions for the Logamediate inflation 2 potential ("2"
!means that inflation proceeds from the left to the right)
!
!V(phi) = M^4 x^[4(1-gam)] exp[-beta * x^gam]
!
!x = (phi-phi0)/Mp


module lmi2sr
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : hypergeom_2F1
  use inftools, only : zbrent
  use lmicommon, only : lmi_alpha
  use lmicommon, only : lmi_norm_potential, lmi_norm_deriv_potential
  use lmicommon, only : lmi_norm_deriv_second_potential
  use lmicommon, only : lmi_epsilon_one_max, lmi_x_potmax, lmi_x_epsonemax
  use lmicommon, only : lmi_epsilon_one, lmi_epsilon_two, lmi_epsilon_three
  use lmicommon, only : lmi_efold_primitive, find_lmi_x_trajectory
  use lmicommon, only : lmi_x_epstwounity, lmi_numacc_x_big, lmi_numacc_dx_potmax
  implicit none

  private

  public lmi2_norm_potential
  public lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  public lmi2_efold_primitive, lmi2_x_trajectory, lmi2_numacc_betamin
  public lmi2_norm_deriv_potential, lmi2_norm_deriv_second_potential
  public lmi2_x_epsoneunity, lmi2_xinimin, lmi2_xendmin
 
contains

!returns V/M^4
  function lmi2_norm_potential(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_norm_potential
    real(kp), intent(in) :: x,gam,beta

    lmi2_norm_potential = lmi_norm_potential(x,gam,beta)

  end function lmi2_norm_potential



!returns the first derivative of the potential with respect to x,
!divided by M^4
  function lmi2_norm_deriv_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_potential
    real(kp), intent(in) :: x,gam,beta
    
    lmi2_norm_deriv_potential = lmi_norm_deriv_potential(x,gam,beta)

  end function lmi2_norm_deriv_potential



!returns the second derivative of the potential with respect to x,
!divided by M^4
  function lmi2_norm_deriv_second_potential(x,gam,beta)
    implicit none
    real(kp) :: lmi2_norm_deriv_second_potential
    real(kp), intent(in) :: x,gam,beta    

    lmi2_norm_deriv_second_potential &
         = lmi_norm_deriv_second_potential(x,gam,beta)
    

  end function lmi2_norm_deriv_second_potential



!epsilon_one(x)
  function lmi2_epsilon_one(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_one
    real(kp), intent(in) :: x,gam,beta

    lmi2_epsilon_one = lmi_epsilon_one(x,gam,beta)
    
  end function lmi2_epsilon_one


!epsilon_two(x)
  function lmi2_epsilon_two(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_two
    real(kp), intent(in) :: x,gam,beta

    lmi2_epsilon_two = lmi_epsilon_two(x,gam,beta)
    
  end function lmi2_epsilon_two


!epsilon_three(x)
  function lmi2_epsilon_three(x,gam,beta)    
    implicit none
    real(kp) :: lmi2_epsilon_three
    real(kp), intent(in) :: x,gam,beta
   
    lmi2_epsilon_three = lmi_epsilon_three(x,gam,beta)
    
  end function lmi2_epsilon_three



!returns the minimum value for xini: if eps1<1, in the whole x>xVmax
!interval, returns max(xVmax,xeps2), otherwise, returns the highest
!solution for eps1=1
  function lmi2_xinimin(gam,beta)
    implicit none
    real(kp), intent(in) :: gam,beta
    real(kp) :: lmi2_xinimin,xeps2one
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    real(kp) ::alpha,xVmax,xeps1max,dxVmax

    alpha = lmi_alpha(gam)

    xVmax = lmi_x_potmax(gam,beta)
    xeps1max = lmi_x_epsonemax(gam,beta)

!at x=xVmax, eps1=0 => cannot be integrated. At x=xVmax+dxVmax, eps1 = numprec
    dxVmax = lmi_numacc_dx_potmax(gam,beta)

    if (lmi_epsilon_one_max(gam,beta).lt.1._kp) then

       xeps2one = lmi_x_epstwounity(gam,beta)

       lmi2_xinimin = max(xeps2one, xVmax+dxVmax)

    else

       lmi2_xinimin = lmi2_x_epsoneunity(gam,beta) &
            * (1._kp+epsilon(1._kp))

    endif

  end function lmi2_xinimin

!returns the minimum value for xend such that lmi2_xinimin and
!xend be searated by efolds e-folds
  function lmi2_xendmin(efolds,gam,beta)
    implicit none
    real(kp), intent(in) :: gam,beta,efolds
    real(kp) :: lmi2_xendmin
    real(kp) :: xendMin, xiniMin

    xiniMin = lmi2_xinimin(gam,beta)
        xendMin = lmi2_x_trajectory(efolds,xiniMin,gam,beta)

    lmi2_xendmin = xendMin

  end function lmi2_xendmin



!returns the value xeps1 at which eps1=1 in the domain x>xvmax (if it
!exists, otherwise stop)
  function lmi2_x_epsoneunity(gam,beta)
    implicit none
    real(kp), intent(in) :: gam, beta
    real(kp) :: lmi2_x_epsoneunity
    real(kp) :: alpha, xeps1max
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: lmi2Data

    
    if (lmi_epsilon_one_max(gam,beta).lt.1._kp) then
       stop 'lmi2_x_epsoneunity: no solution for eps1=1'
    endif

    xeps1max = lmi_x_epsonemax(gam,beta)
    alpha = lmi_alpha(gam)

    mini = xeps1max * (1._kp+epsilon(1._kp))
!    maxi = lmi_numacc_x_big(gam,beta) * max(alpha,(beta*gam)**(1._kp/(1._kp-gam)) &
!         ,(alpha*beta*gam)**(1._kp/(2._kp-gam)))
    maxi = epsilon(1._kp)*huge(1._kp)
      
    lmi2Data%real1 = gam
    lmi2Data%real2 = beta
    
    lmi2_x_epsoneunity &
         = zbrent(find_lmi2_x_epsoneunity,mini,maxi,tolFind,lmi2Data)
        
  end function lmi2_x_epsoneunity



  function find_lmi2_x_epsoneunity(x,lmi2Data)    
    implicit none
    real(kp), intent(in) :: x   
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2_x_epsoneunity
    real(kp) :: gam,beta

    gam = lmi2Data%real1
    beta = lmi2Data%real2
    
    find_lmi2_x_epsoneunity = lmi2_epsilon_one(x,gam,beta) - 1._kp
   
  end function find_lmi2_x_epsoneunity




!this is integral[V(phi)/V'(phi) dphi]
  function lmi2_efold_primitive(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi2_efold_primitive

    lmi2_efold_primitive = lmi_efold_primitive(x,gam,beta)
    
  end function lmi2_efold_primitive


!this is integral[V(phi)/V'(phi) dphi], approximated in the limit x/x0>>1
  function lmi2_efold_primitive_approximated(x,gam,beta)
    implicit none
    real(kp), intent(in) :: x,gam,beta
    real(kp) :: lmi2_efold_primitive_approximated

    real(kp) ::alpha
    alpha = lmi_alpha(gam)

    if (gam.eq.0._kp) stop 'lmi2_efold_primitive: gam=0!'

    lmi2_efold_primitive_approximated = x**(2._kp-gam) &
         /(beta*gam*(gam-2._kp))

  end function lmi2_efold_primitive_approximated


!returns x at bfold=-efolds before the end of inflation, ie N-Nend
  function lmi2_x_trajectory(bfold,xend,gam,beta)
    implicit none
    real(kp), intent(in) :: bfold, gam, xend,beta
    real(kp) :: lmi2_x_trajectory
    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi, xiniMin, deltaN
    type(transfert) :: lmiData

    real(kp) ::alpha

    alpha = lmi_alpha(gam)
    xiniMin = lmi2_xinimin(gam,beta)

    if (xend.le.xiniMin .and. bfold.lt.0._kp) then
       write(*,*)'xiniMin= xend= ',xiniMin, xend
       stop 'lmi2_x_trajectory: xend < xiniMin'
    endif

    if (bfold.lt.0._kp) then
      maxi = xend * (1._kp+epsilon(1._kp))
      mini = xiniMin * (1._kp-epsilon(1._kp))
    else
      mini=xend
      deltaN = 10._kp*bfold
!Uses the approximated formula of the trajectory when x/x0>>1 note:
!this is only to get a reasonable value of maxi, then the real
!trajectory is solved
      maxi = mini*( &
           1._kp + beta*gam*(2._kp-gam)*deltaN/mini**(2._kp-gam) &
           )**(1._kp/(2._kp-gam))
    endif

    lmiData%real1 = gam
    lmiData%real2 = beta
    lmiData%real3 = -bfold + lmi_efold_primitive(xend,gam,beta)
    
    lmi2_x_trajectory = zbrent(find_lmi_x_trajectory,mini,maxi,tolFind,lmiData)

       
  end function lmi2_x_trajectory



 
! Returns the minimum value for beta in order to end inflation with
! slow roll violation ( beta>betamin(gam) <=> epsOneMax>1 )
  function lmi2_betamin(gam)
    implicit none
    real(kp), intent(in) :: gam
    real(kp) :: lmi2_betamin
    real(kp) ::alpha

    alpha = lmi_alpha(gam)

    lmi2_betamin=(sqrt(2._kp)*(1._kp-gam)/(alpha*gam))**gam &
         *alpha/(gam*(1._kp-gam))

  end function lmi2_betamin

! Returns the minimum value for gam in order to end inflation with
! slow roll violation ( gam>gammin(beta) <=> epsOneMax>1 )
  function lmi2_gammin(beta)
    implicit none
    real(kp), intent(in) :: beta
    real(kp) :: lmi2_gammin
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: lmi2Data

    if (beta.lt.sqrt(2._kp)) then
       stop 'lmi2_gammin: beta<sqrt(2): inflation cannot end by slow-roll violation!'
    endif

    lmi2Data%real1 = beta

    lmi2_gammin = zbrent(find_lmi2_gammin,epsilon(1._kp) &
         ,1._kp-epsilon(1._kp),tolFind,lmi2Data)

  end function lmi2_gammin

  function find_lmi2_gammin(x,lmi2Data)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data
    real(kp) :: find_lmi2_gammin
    real(kp)  ::  beta

    beta = lmi2Data%real1
    
    find_lmi2_gammin = lmi2_betamin(x)-beta

  end function find_lmi2_gammin


!return betamin such that epsonemax > machineprecision
  function lmi2_numacc_betamin(gam)
    implicit none
    real(kp) :: lmi2_numacc_betamin
    real(kp), intent(in) :: gam
    real(kp) :: epsmin

    epsmin = tolkp

    lmi2_numacc_betamin = 4._kp/gam * (epsmin/8._kp/gam/gam)**(gam/2._kp)

  end function lmi2_numacc_betamin


end module lmi2sr

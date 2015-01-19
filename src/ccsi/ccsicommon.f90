!slow-roll functions for the R + R^2/m^2 + alpha R^3/m^4 inflation potential
!
!V(phi) = M^4 * exp(-2x) * [ exp(x) - 1 ]^2/{1 + sqrt[1+3 alpha(exp(x)-1)]}^3 
!  * {1 + sqrt[1 + 3 alpha (exp(x)-1)] + 2 alpha (exp(x)-1) }
!
!x = phi/Mp * sqrt(2/3)


module ccsicommon
  use infprec, only : kp,tolkp,transfert
  use specialinf, only : lambert
  use inftools, only : zbrent
  implicit none

  private

  public ccsi_norm_potential, ccsi_epsilon_one, ccsi_epsilon_two, ccsi_epsilon_three
  public ccsi_norm_deriv_potential, ccsi_norm_deriv_second_potential
  public ccsih_x_trajectory, ccsih_efold_primitive, ccsi_xmax, ccsi_x_potmax
  public ccsi_x_epsoneunity, ccsi_x_epstwonull, ccsi_alphamin
  public ccsi_efold_primitive, find_ccsi_x_trajectory

  public ccsi_numacc_x_epsonenull

!because we use exp(x) everywhere
  real(kp), parameter :: ccsiBig = log(epsilon(1._kp)*huge(1._kp))
  real(kp), parameter :: ccsiSmall = epsilon(1._kp)

  public ccsiBig, ccsiSmall

contains

!returns V/M^4
  function ccsi_norm_potential(x,alpha)
    implicit none
    real(kp) :: ccsi_norm_potential
    real(kp), intent(in) :: x,alpha
    real(kp) :: aexm1,s1p3aexm1

    aexm1 = alpha*(exp(x)-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_norm_potential = (1._kp - exp(-x))**2/(1._kp+s1p3aexm1)**3 &
         * (1._kp + 2._kp*aexm1 + s1p3aexm1)

  end function ccsi_norm_potential



!returns the first derivative of the potential with respect to x, divided by M^4
  function ccsi_norm_deriv_potential(x,alpha)
    implicit none
    real(kp) :: ccsi_norm_deriv_potential
    real(kp), intent(in) :: x,alpha
    real(kp) :: ex,aexm1,s1p3aexm1

    ex = exp(x)
    aexm1 = alpha*(ex-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_norm_deriv_potential = exp(-2._kp*x)*(4._kp - 4._kp*s1p3aexm1 &
         -  3._kp*alpha*(6._kp - 4._kp*s1p3aexm1 &
         +  ex*(-3._kp + s1p3aexm1)))/(27._kp*alpha**2)

  end function ccsi_norm_deriv_potential



!returns the second derivative of the potential with respect to x, divided by M^4
  function ccsi_norm_deriv_second_potential(x,alpha)
    implicit none
    real(kp) :: ccsi_norm_deriv_second_potential
    real(kp), intent(in) :: x,alpha
    real(kp) :: ex, aexm1, s1p3aexm1

    ex = exp(x)
    aexm1 = alpha*(ex-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_norm_deriv_second_potential = exp(-2._kp*x) * (-16._kp*(-1._kp + s1p3aexm1) &
         + 3._kp*alpha*(ex*(14._kp + 3._kp*(-14._kp + ex)*alpha &
         - 6._kp*s1p3aexm1) + 8._kp*(-4._kp + 6._kp*alpha &
         + 3._kp*s1p3aexm1))) &
         / (54._kp*alpha**2*s1p3aexm1)
     
  end function ccsi_norm_deriv_second_potential



!epsilon_one(x)
  function ccsi_epsilon_one(x,alpha)    
    implicit none
    real(kp) :: ccsi_epsilon_one
    real(kp), intent(in) :: x,alpha
    real(kp) :: emx,aexm1,s1p3aexm1

    emx = exp(-x)
    aexm1 = alpha*(1._kp/emx-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_epsilon_one = (emx*(2._kp - 8._kp*alpha) - 1._kp - 2._kp*aexm1 &
         + 8._kp*alpha + s1p3aexm1)**2/( (1._kp - emx)*(1._kp + 4._kp*aexm1))**2/3._kp
    
    
  end function ccsi_epsilon_one


!epsilon_two(x)
  function ccsi_epsilon_two(x,alpha)    
    implicit none
    real(kp) :: ccsi_epsilon_two
    real(kp), intent(in) :: x,alpha
    real(kp) :: ex,aexm1,s1p3aexm1

    ex = exp(x)
    aexm1 = alpha*(ex-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_epsilon_two = (2._kp*ex*(2._kp*(1._kp + s1p3aexm1) &
         + 12._kp*aexm1**2*(2._kp + ex + 4._kp*s1p3aexm1) &
         - aexm1*(ex*(-5._kp + 4._kp*s1p3aexm1) - 2._kp*(7._kp + 10._kp*s1p3aexm1)))) &
         /(3._kp*sqrt(1._kp + 3._kp*aexm1)*(1._kp + 4._kp*aexm1)**2*(-1._kp + ex)**2)
   
  end function ccsi_epsilon_two


!epsilon_three(x)
  function ccsi_epsilon_three(x,alpha)    
    implicit none
    real(kp) :: ccsi_epsilon_three
    real(kp), intent(in) :: x,alpha
    real(kp) :: ex,aexm1,s1p3aexm1

    ex = exp(x)
    aexm1 = alpha*(ex-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)

    ccsi_epsilon_three = -alpha**2 * ((-2._kp + 8._kp*alpha + ex*(1._kp - s1p3aexm1 &
         + 2._kp*(-5._kp + ex)*alpha))*(144._kp*ex**6*alpha**4 + 24._kp*ex**5*alpha**3 &
         *(3._kp - 4._kp*s1p3aexm1 + 12._kp*(3._kp + 4._kp*s1p3aexm1)*alpha) &
         + 4._kp*(1._kp - 4._kp*alpha)**2*(-1._kp + 3._kp*alpha) &
         *(-1._kp - s1p3aexm1 + (3._kp + 6._kp*s1p3aexm1)*alpha) - ex**4*alpha**2 &
         *(23._kp + 8._kp*s1p3aexm1 + 24._kp*alpha*(-55._kp - 48._kp*s1p3aexm1 &
         + 36._kp*(5._kp + 4._kp*s1p3aexm1)*alpha)) + 3._kp*ex**2*alpha &
         *(-1._kp + 4._kp*alpha)*(-4._kp*(8._kp + 3._kp*s1p3aexm1) + alpha &
         *(165._kp - 16._kp*s1p3aexm1 + 12._kp*(-15._kp + 16._kp*s1p3aexm1)*alpha)) &
         + 2._kp*ex**3*alpha*(-5._kp + 4._kp*s1p3aexm1 + 2._kp*alpha*(155._kp + 68._kp*s1p3aexm1 &
         + 6._kp*alpha*(-155._kp - 72._kp*s1p3aexm1 + 48._kp*(5._kp + 2._kp*s1p3aexm1)*alpha))) &
         - 2._kp*ex*(-1._kp + 4._kp*alpha)*(2._kp*(1._kp + s1p3aexm1) + alpha &
         *(-7._kp + 20._kp*s1p3aexm1 + 3._kp*alpha*(-11._kp - 72._kp*s1p3aexm1 &
         + 36._kp*(1._kp + 4._kp*s1p3aexm1)*alpha))))) &
         / &
        (3._kp*aexm1**2*(1._kp + 3._kp*aexm1)*(1._kp + 4._kp*aexm1)**2 &
        *(2._kp*(1._kp + s1p3aexm1) + 12._kp*aexm1**2*(2._kp + ex + 4._kp*s1p3aexm1) &
        - aexm1*(ex*(-5._kp + 4._kp*s1p3aexm1) - 2._kp*(7._kp + 10._kp*s1p3aexm1))))
 

  end function ccsi_epsilon_three

!returns alpha such that the min value of epsone is unity. For value
!of alpha less than this, epsilon1 is always greater than unity and
!there is no inflation
  function ccsi_alphamin()
    implicit none
    real(kp) :: ccsi_alphamin

    ccsi_alphamin = 3._kp/44._kp * (1._kp - 2._kp*sqrt(3._kp))

  end function ccsi_alphamin

!above this value the potential becomes complex for alpha < 0
  function ccsi_xmax(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: ccsi_xmax

    if (alpha.ge.0._kp) stop 'ccsi_xmax: x not bounded for alpha>=0'

    ccsi_xmax = log(1._kp - 1._kp/3._kp/alpha)

  end function ccsi_xmax



!returns the solution of eps1=1
  function ccsi_x_epsoneunity(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp), dimension(2) :: ccsi_x_epsoneunity
    real(kp) :: xmax

    if (alpha.le.ccsi_alphamin()) then
       stop 'ccsi_x_epsoneunity: alpha < alphamin, eps1>1 everywhere!'
    endif

!smaller root
    ccsi_x_epsoneunity(1) = log((-15._kp - 14._kp*sqrt(3._kp) + 176._kp*alpha &
         + 132._kp*sqrt(3._kp)*alpha + sqrt(813._kp + 420._kp*sqrt(3._kp) &
         + 4444._kp*alpha + 2728*sqrt(3._kp)*alpha))/(242._kp*alpha))

    if (alpha.ge.0._kp) then
       ccsi_x_epsoneunity(2) = ccsi_x_epsoneunity(1)
       return
    endif
    
!the other root exists only for alphamin < alpha < 0 provided eps1(xmax)>1
    xmax = ccsi_xmax(alpha)

    if (ccsi_epsilon_one(xmax,alpha).le.1._kp) then
       ccsi_x_epsoneunity(2) = ccsi_x_epsoneunity(1)
    else
       ccsi_x_epsoneunity(2) = log(-(15._kp + 14._kp*sqrt(3._kp) - (176._kp &
            + 132._kp*sqrt(3._kp))*alpha + sqrt(813._kp + 420._kp*sqrt(3._kp) &
            + 4444._kp*alpha + 2728._kp*sqrt(3._kp)*alpha))/(242._kp*alpha))
!safeguard
       if (ccsi_x_epsoneunity(2).gt.xmax) stop 'ccsi_x_epsoneunity: internal error!'

    endif
    
  end function ccsi_x_epsoneunity


!field value at which epsilon two vanishes (epsilon one extremum)
  function ccsi_x_epstwonull(alpha)
    implicit none
    real(kp), intent(in) :: alpha
    real(kp) :: ccsi_x_epstwonull

    if (alpha.ge.0) then
       stop 'ccsi_x_epstwonull: no vanishing eps2 for alpha > 0!'
    endif

    ccsi_x_epstwonull = log(-6._kp + 24._kp*alpha - (2._kp*sqrt(alpha &
         *(-1._kp + 4._kp*alpha)*(-1._kp + 6._kp*alpha)**2))/alpha)

  end function ccsi_x_epstwonull



!field value at which the potential is maximal (only for alpha>0)
  function ccsi_x_potmax(alpha)
    implicit none
    real(kp) , intent(in) :: alpha
    real(kp) :: ccsi_x_potmax

    if (alpha.le.0._kp) stop 'ccsi_x_potmax: no maximum, alpha<=0!'

    ccsi_x_potmax = log((2._kp*(1._kp + 2._kp*sqrt(alpha)))/sqrt(alpha))

  end function ccsi_x_potmax


!return x such that epsilon_1 = ccsiSmall for checking numerical accuracy  
  function ccsi_numacc_x_epsonenull(alpha)
    implicit none
    real(kp), dimension(2) :: ccsi_numacc_x_epsonenull
    real(kp), intent(in) :: alpha
    real(kp) :: epsMini


    epsMini = ccsiSmall
    
    ccsi_numacc_x_epsonenull(1) = log((3._kp*(8._kp*alpha + 3._kp*EpsMini*(-1._kp + 8._kp*alpha)&
         + 2._kp*sqrt(3._kp*EpsMini)*(-1._kp + 10._kp*alpha) &
         + sqrt(9._kp*EpsMini**2 + 16._kp*alpha + 48._kp*sqrt(3._kp*EpsMini)*alpha &
         + 12._kp*sqrt(3._kp*EpsMini**3)*(1._kp + 2._kp*alpha) + 12._kp*EpsMini &
         * (1._kp + 9._kp*alpha))))/(2._kp*(sqrt(3._kp) + 6._kp*sqrt(EpsMini))**2*alpha))

    ccsi_numacc_x_epsonenull(2) = log((3._kp*(2._kp*sqrt(3._kp*EpsMini)*(1._kp - 10._kp*alpha) &
         + 8._kp*alpha +  3._kp*EpsMini*(-1._kp + 8._kp*alpha) &
         +  sqrt(9._kp*EpsMini**2 + 16._kp*alpha - 48._kp*sqrt(3._kp*EpsMini)*alpha &
         - 12._kp*sqrt(3._kp*EpsMini**3)*(1._kp + 2._kp*alpha) &
         + 12._kp*EpsMini*(1._kp + 9._kp*alpha)))) &
         / (2._kp*(sqrt(3._kp) - 6._kp*sqrt(EpsMini))**2*alpha))


  end function ccsi_numacc_x_epsonenull
  

!this is integral[V(phi)/V'(phi) dphi]
  function ccsi_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ccsi_efold_primitive

    real(kp) :: ex, aexm1, s1p3aexm1
    complex(kp) :: sqralpha, cprim

    ex = exp(x)
    aexm1 = alpha*(ex-1._kp)
    s1p3aexm1 = sqrt(1._kp + 3._kp*aexm1)
    sqralpha = sqrt(cmplx(alpha,0._kp,kp))

!this is a complex function with non-vanishing imaginary
!part. However, provided the trajectory do not cross singular point,
!the imaginary part is constant for the same branch cut and cancels
!out in the trajectory. Therefore, we only return the real part of
!this complex function

    cprim = -2._kp*x + (2*atanh((sqralpha*(-4._kp + ex))/2._kp))/sqralpha &
         + 6._kp*atanh(s1p3aexm1/(1._kp - 3._kp*sqralpha)) &
         + (2._kp*atanh(s1p3aexm1/(-1._kp + 3._kp*sqralpha)))/sqralpha &
         + 6._kp*atanh(s1p3aexm1/(1._kp + 3._kp*sqralpha)) &
         + (2._kp*atanh(s1p3aexm1/(1._kp + 3._kp*sqralpha)))/sqralpha &
         - 3._kp*log(abs(4._kp - alpha*(-4._kp + ex)**2))
  
    ccsi_efold_primitive = real(cprim,kp) * 3._kp/8._kp

  end function ccsi_efold_primitive

  function find_ccsi_x_trajectory(x,ccsiData)    
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsiData
    real(kp) :: find_ccsi_x_trajectory
    real(kp) :: alpha,NplusNuend

    alpha = ccsiData%real1
    NplusNuend = ccsiData%real2

    find_ccsi_x_trajectory = ccsi_efold_primitive(x,alpha) - NplusNuend
   
  end function find_ccsi_x_trajectory


!Higgs Inflation Model (HI)
  function ccsih_efold_primitive(x,alpha)
    implicit none
    real(kp), intent(in) :: x,alpha
    real(kp) :: ccsih_efold_primitive

    if (alpha.ne.0._kp) stop 'ccsih_efold_primitive: alpha is not vanishing!'
   
    ccsih_efold_primitive = 3._kp/4._kp*(exp(x)-x) 

  end function ccsih_efold_primitive


  
  function ccsih_x_trajectory(bfold,xend,alpha)
    implicit none
    real(kp), intent(in) :: bfold, alpha, xend
    real(kp) :: ccsih_x_trajectory

    if (alpha.ne.0._kp) stop 'ccsih_x_trajectory: alpha is not vanishing!'

    ccsih_x_trajectory = (4._kp/3._kp*bfold+xend-exp(xend)- &
         lambert(-exp(4._kp/3._kp*bfold+xend-exp(xend)),-1))
      
  end function ccsih_x_trajectory



  


end module ccsicommon

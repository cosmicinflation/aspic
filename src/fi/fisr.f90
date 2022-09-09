!slow-roll functions for the fibre inflation potential
!
!V(phi) =M^4 [ (1+2 delta/3) Exp(-2 x/sqrt(3)) - 4(1+delta/6) Exp(-x/sqrt(3)) + delta/(1+n)Exp(2(1+n)x/sqrt(3)) + 3 - delta/(1+n) )
!
!x = phi/Mp

module fisr
  use infprec, only : kp, toldp, tolkp, transfert
  use inftools, only : zbrent, easydverk
  implicit none

  private

  public fi_norm_potential, fi_norm_deriv_potential, fi_norm_deriv_second_potential
  public fi_epsilon_one, fi_epsilon_two, fi_epsilon_three
  public fi_efold_primitive, fi_x_trajectory, fi_x_endinf, fi_x_epstwozero
  public fi_check_params, fi_x_epsoneunity
  public fi_efoldmax, fi_epsilon_one_min

  real(kp), parameter :: sqr3 = sqrt(3._kp)

contains


  function fi_check_params(efoldNum,epsilonMax,delta,n)
    implicit none
    logical :: fi_check_params
    real(kp), intent(in) :: efoldNum, epsilonMax
    real(kp), intent(in) :: delta,n

    real(kp) :: Nmax,epsilon1min

    Nmax = fi_efoldmax(delta,n)
    epsilon1min = fi_epsilon_one_min(delta,n)

    fi_check_params = ((epsilon1min .le. epsilonMax) .and. (Nmax .ge. efoldNum) )

  end function fi_check_params



!returns V/M**4
  function fi_norm_potential(x,delta,n)
    implicit none
    real(kp) :: fi_norm_potential
    real(kp), intent(in) :: x,delta,n

    fi_norm_potential = (1._kp+2._kp*delta/3._kp)*exp(-4._kp*x/sqr3)- &
                        4._kp*(1._kp+delta/6._kp)*exp(-x/sqr3)+delta/(1._kp+n)* &
                        exp(2._kp*(1._kp+n)*x/sqr3)+3._kp-delta/(1._kp+n)

  end function fi_norm_potential


!returns the first derivative of the potential with respect to x=phi/Mp, divided by M**4
  function fi_norm_deriv_potential(x,delta,n)
    implicit none
    real(kp) :: fi_norm_deriv_potential
    real(kp), intent(in) :: x,delta,n

    fi_norm_deriv_potential = (2._kp*exp(-((4._kp*x)/sqr3))* &
                                (-6._kp-4._kp*delta+3._kp*exp((2._kp*(3._kp+n)*x)/ &
                                sqr3)*delta+exp(x*sqr3)*(6._kp+delta)))/ &
                                (3._kp*sqr3)

  end function fi_norm_deriv_potential


!returns the second derivative of the potential with respect to x=phi/Mp, divided by M**4
  function fi_norm_deriv_second_potential(x,delta,n)
    implicit none
    real(kp) :: fi_norm_deriv_second_potential
    real(kp), intent(in) :: x,delta,n

    fi_norm_deriv_second_potential = 2._kp/9._kp*exp(-((4._kp*x)/sqr3))* &
                                    (24._kp+16._kp*delta+6._kp*exp((2._kp*(3._kp+n)*x)/sqr3)* &
                                    (1._kp+n)*delta-exp(x*sqr3)*(6._kp+delta))

  end function fi_norm_deriv_second_potential



!returns the third derivative of the potential with respect to x=phi/Mp, divided by M**4
  function fi_norm_deriv_third_potential(x,delta,n)
    implicit none
    real(kp) :: fi_norm_deriv_third_potential
    real(kp), intent(in) :: x,delta,n

    fi_norm_deriv_third_potential = (2._kp*exp(-((4._kp*x)/sqr3))* &
                                    (12._kp*exp((2._kp*(3._kp+n)*x)/sqr3)* &
                                    (1._kp+n)**2*delta+exp(x*sqr3)*(6._kp+delta) &
                                    -32._kp*(3._kp+2._kp*delta)))/(9._kp*sqr3)

  end function fi_norm_deriv_third_potential




  function fi_epsilon_one(x,delta,n)
    implicit none
    real(kp) :: fi_epsilon_one
    real(kp), intent(in) :: x,delta,n

    fi_epsilon_one = (2*(-6 - 4*delta + (6 + delta)*exp(sqr3*x) &
         + 3*delta*exp((2*(3 + n)*x)/sqr3))**2*(1 + n)**2) &
         / (3*(3*delta*exp((2*(3 + n)*x)/sqr3) + (3 + 2*delta)*(1 + n) &
         - 2*(6 + delta)*exp(sqr3*x)*(1 + n) + exp((4*x)/sqr3)*(9 - 3*delta + 9*n))**2)

!   print *,'test',fi_epsilon_one, 0.5_kp*(fi_norm_deriv_potential(x,delta,n)/ &
!                    fi_norm_potential(x,delta,n))**2

  end function fi_epsilon_one



  function fi_epsilon_two(x,delta,n)
    implicit none
    real(kp) :: fi_epsilon_two
    real(kp), intent(in) :: x,delta,n


    fi_epsilon_two = (4*(1 + n)*(3*(18 + 15*delta + 2*delta**2) &
         * exp(sqr3*x)*(1 + n) + 6*delta*exp((2*(5 + n)*x)/sqr3) &
         * (-3 + delta - 3*n)*(1 + n) - 2*delta*(3 + 2*delta) &
         * exp((2*(3 + n)*x)/sqr3)*(3 + n)**2 + delta*(6 + delta) &
         * exp(((9 + 2*n)*x)/sqr3)*(3 + 2*n)**2 - 8*(3 + 2*delta)*exp((4*x)/sqr3) &
         * (3 - delta + 3*n) + (6 + delta)*exp((7*x)/sqr3)*(3 - delta + 3*n))) &
         / (3*delta*exp((2*(3 + n)*x)/sqr3) + (3 + 2*delta)*(1 + n) - 2*(6 + delta) &
         * exp(sqr3*x)*(1 + n) + exp((4*x)/sqr3)*(9 - 3*delta + 9*n))**2

!    print *,'test2',fi_epsilon_two,2._kp*((fi_norm_deriv_potential(x,delta,n)/ &
!                     fi_norm_potential(x,delta,n))**2- &
!                     fi_norm_deriv_second_potential(x,delta,n)/ &
!                     fi_norm_potential(x,delta,n))

  end function fi_epsilon_two


  function fi_epsilon_three(x,delta,n)
    implicit none
    real(kp) :: fi_epsilon_three
    real(kp), intent(in) :: x,delta,n
!    real(kp) :: V,Vp,Vpp,Vppp

    fi_epsilon_three = (-2*(-6 - 4*delta + (6 + delta)*exp(sqr3*x) &
         + 3*delta*exp((2*(3 + n)*x)/sqr3))*(1 + n)*(-4*(-3*(6 + delta)*exp(sqr3*x)*(1 + n)&
         + 3*delta*exp((2*(3 + n)*x)/sqr3)*(3 + n) + 6*exp((4*x)/sqr3)*(3 - delta + 3*n)) &
         * (3*(18 + 15*delta + 2*delta**2)* exp(sqr3*x)*(1 + n) + 6*delta*exp((2*(5 + n)*x)/sqr3) &
         * (-3 + delta - 3*n)*(1 + n) - 2*delta*(3 + 2*delta)*exp((2*(3 + n)*x)/sqr3)*(3 + n)**2 &
         + delta*(6 + delta)*exp(((9 + 2*n)*x)/sqr3)*(3 + 2*n)**2 - 8*(3 + 2*delta)*exp((4*x)/sqr3) &
         * (3 - delta + 3*n) + (6 + delta)*exp((7*x)/sqr3)*(3 - delta + 3*n)) &
         + (9*(18 + 15*delta + 2*delta**2)*exp(sqr3*x)*(1 + n) - 4*delta*(3 + 2*delta) &
         * exp((2*(3 + n)*x)/sqr3)*(3 + n)**3 + 12*delta*exp((2*(5 + n)*x)/sqr3) &
         * (-3 + delta - 3*n)*(1 + n)*(5 + n) + delta*(6 + delta)*exp(((9 + 2*n)*x)/sqr3) &
         *(3 + 2*n)**2*(9 + 2*n) - 32*(3 + 2*delta)*exp((4*x)/sqr3)*(3 - delta + 3*n) &
         + 7*(6 + delta)*exp((7*x)/sqr3)*(3 - delta + 3*n))*(3*delta*exp((2*(3 + n)*x)/sqr3)&
         + (3 + 2*delta)*(1 + n) - 2*(6 + delta)*exp(sqr3*x)*(1 + n) + exp((4*x)/sqr3) &
         *(9 - 3*delta + 9*n))))/(3*(3*(18 + 15*delta + 2*delta**2)*exp(sqr3*x)*(1 + n) &
         + 6*delta*exp((2*(5 + n)*x)/sqr3)*(-3 + delta - 3*n)*(1 + n) - 2*delta*(3 + 2*delta) &
         * exp((2*(3 + n)*x)/sqr3)*(3 + n)**2 + delta*(6 + delta)*exp(((9 + 2*n)*x)/sqr3) &
         * (3 + 2*n)**2 - 8*(3 + 2*delta)*exp((4*x)/sqr3)*(3 - delta + 3*n) + (6 + delta)*exp((7*x)/sqr3) &
         * (3 - delta + 3*n))*(3*delta*exp((2*(3 + n)*x)/sqr3) + (3 + 2*delta)*(1 + n) - 2*(6 + delta) &
         * exp(sqr3*x)*(1 + n) + exp((4*x)/sqr3)*(9 - 3*delta + 9*n))**2)


!    V = fi_norm_potential(x,delta,n)
!    Vp = fi_norm_deriv_potential(x,delta,n)
!    Vpp = fi_norm_deriv_second_potential(x,delta,n)
!    Vppp = fi_norm_deriv_third_potential(x,delta,n)
!    print *,'test3',fi_epsilon_three, 2*(Vppp*Vp/V**2-3._kp*Vpp*Vp**2/V**3+2._kp*(Vp/V)**4)/ &
!                        fi_epsilon_two(x,delta,n)

  end function fi_epsilon_three



!this is x at which the third branch of the potential becomes dominant
  function fi_x_potthirdbranch(delta,n)
    implicit none
    real(kp) :: fi_x_potthirdbranch
    real(kp), intent(in) :: delta,n

    fi_x_potthirdbranch = ((sqr3*log(((1._kp+n)*(1._kp+(2._kp*delta)/3._kp))/delta)) &
         /(2._kp*(3._kp+n)))/10._kp

  end function fi_x_potthirdbranch

!maximal value of x above which the potential overflows
  function fi_numacc_xmax(delta,n)
    implicit none
    real(kp) :: fi_numacc_xmax
    real(kp), intent(in) :: delta,n

    fi_numacc_xmax = (sqr3*(log(1._kp/epsilon(1._kp))))/(2._kp*(1._kp+n))

  end function fi_numacc_xmax




!returns x where epsilon2=0 and epsilon1 is minimal
  function fi_x_epstwozero(delta,n)
    implicit none
    real(kp) :: fi_x_epstwozero, mini, maxi
    real(kp), intent(in) :: delta,n
    real(kp), parameter :: tolFind=tolkp
    type(transfert) :: fiData

    mini = fi_x_potthirdbranch(delta,n)
    maxi = fi_numacc_xmax(delta,n)

    fiData%real1 = delta
    fiData%real2 = n

    fi_x_epstwozero = zbrent(find_fi_x_epstwozero,mini,maxi,tolFind,fiData)

  end function fi_x_epstwozero

  function find_fi_x_epstwozero(x,fiData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData
    real(kp) :: find_fi_x_epstwozero
    real(kp) :: delta,n

    delta = fiData%real1
    n = fiData%real2

    find_fi_x_epstwozero = fi_epsilon_two(x,delta,n)

  end function find_fi_x_epstwozero



!returns the minimum value of epsilon1
  function fi_epsilon_one_min(delta,n)
    implicit none
    real(kp), intent(in) :: delta,n
    real(kp) :: fi_epsilon_one_min

    fi_epsilon_one_min = fi_epsilon_one(fi_x_epstwozero(delta,n),delta,n)

  end function fi_epsilon_one_min




  function fi_x_epsoneunity(delta,n)
    implicit none
    real(kp), dimension(2) :: fi_x_epsoneunity
    real(kp), intent(in) :: delta,n

    real(kp) :: mini,maxi
    real(kp) :: xepstwonull, xnumaccMax, xthirdbranch

    real(kp), parameter :: tolFind = tolkp
    type(transfert) :: fiData

    xepstwonull = fi_x_epstwozero(delta,n)
    xthirdbranch = fi_x_potthirdbranch(delta,n)
    xnumaccMax = fi_numacc_xmax(delta,n)

    fiData%real1 = delta
    fiData%real2 = n

    mini = xthirdbranch
    maxi = xepstwonull

    fi_x_epsoneunity(1) = zbrent(find_fi_x_epsoneunity,mini,maxi,tolFind,fiData)


    mini = xepstwonull
    maxi = xnumaccMax

! in this case epsilon1 reaches 2/3<1 at infinity, and the function
! returns the maximal value of x above which the potential becomes
! larger than numerical accuracy
    if (n .eq. 0._kp) then

       fi_x_epsoneunity(2) = xnumaccMax

    else

       fi_x_epsoneunity(2) = zbrent(find_fi_x_epsoneunity,mini,maxi,tolFind,fiData)

    endif

  end function fi_x_epsoneunity



  function find_fi_x_epsoneunity(x,fiData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData
    real(kp) :: find_fi_x_epsoneunity
    real(kp) :: delta,n

    delta = fiData%real1
    n = fiData%real2

    find_fi_x_epsoneunity = fi_epsilon_one(x,delta,n) - 1._kp

  end function find_fi_x_epsoneunity




!returns x at the end of inflation defined as epsilon1=1
  function fi_x_endinf(delta,n)
    implicit none
    real(kp) :: fi_x_endinf
    real(kp), intent(in) :: delta,n
    real(kp), dimension(2) :: xepsone

    xepsone = fi_x_epsoneunity(delta,n)

    fi_x_endinf = xepsone(1)

  end function fi_x_endinf




!this is integral[V(phi)/V'(phi) dphi]
  function fi_efold_primitive(x,delta,n)
    implicit none
    real(kp), intent(in) :: x,delta,n
    real(kp) :: fi_efold_primitive
    type(transfert) :: fiData
!avoids prohibitive integration time in QUADPREC
    real(kp), parameter :: tolInt = max(toldp,tolkp)
    integer, parameter :: neq = 1
    real(kp) :: xvar, xinf
    real(kp), dimension(neq) :: yvar

    !let us start where epsilon2 vanishes
    xvar = fi_x_epstwozero(delta,n)
    yvar(1) = 0._kp

    fiData%real1 = delta
    fiData%real2 = n

    call easydverk(neq,find_fi_efold_primitive,xvar,yvar,x,tolInt,fiData)

    fi_efold_primitive = yvar(1)

  end function fi_efold_primitive

  subroutine find_fi_efold_primitive(n,x,y,yprime,fiData)
    implicit none
    integer :: n
    real(kp) :: x
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: fiData
    real(kp) :: delta, nfi, x4p2n, phi4p2n

    delta = fiData%real1
    nfi = fiData%real2

    yprime(1) = ((1._kp+2._kp*delta/3._kp)*exp(-4._kp*x/sqr3) &
         - 4._kp*(1._kp+delta/6._kp)*exp(-x/sqr3)+delta/(1._kp+nfi) &
         * exp(2._kp*(1._kp+nfi)*x/sqr3)+3._kp-delta/(1._kp+nfi)) &
         / ((2._kp*exp(-((4._kp*x)/sqr3)) &
         * (-6._kp-4._kp*delta+3._kp*exp((2._kp*(3._kp+nfi)*x) &
         /sqr3)*delta+exp(x*sqr3)*(6._kp+delta)))/(3._kp*sqr3) )

  end subroutine find_fi_efold_primitive



!returns x at bfold=-efolds before the end of inflation
  function fi_x_trajectory(bfold,xend,delta,n)
    implicit none
    real(kp), intent(in) :: bfold,xend,delta,n
    real(kp) :: fi_x_trajectory

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi
    type(transfert) :: fiData

    real(kp), dimension(2) :: xepsone

    xepsone = fi_x_epsoneunity(delta,n)

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xepsone(2)*(1._kp-epsilon(1._kp))

    fiData%real1 = delta
    fiData%real2 = n
    fiData%real3 = -bfold + fi_efold_primitive(xend,delta,n)

    fi_x_trajectory = zbrent(find_fi_x_trajectory,mini,maxi,tolFind,fiData)

  end function fi_x_trajectory

  function find_fi_x_trajectory(x,fiData)
    implicit none
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData
    real(kp) :: find_fi_x_trajectory
    real(kp) :: delta,n,NplusPrimEnd

    delta = fiData%real1
    n = fiData%real2
    NplusPrimEnd = fiData%real3

    find_fi_x_trajectory = fi_efold_primitive(x,delta,n) - NplusPrimEnd

  end function find_fi_x_trajectory


!Maximal number of slow-roll efolds
  function fi_efoldmax(delta,n)
    implicit none
    real(kp) :: fi_efoldmax
    real(kp), intent(in) :: delta,n

    real(kp), dimension(2) :: xepsone

    xepsone = fi_x_epsoneunity(delta,n)

    fi_efoldmax = -fi_efold_primitive(xepsone(1),delta,n) + fi_efold_primitive(xepsone(2),delta,n)

  end function fi_efoldmax



end module fisr

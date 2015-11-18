module dicommon
  use infprec, only : kp, pi, tolkp, transfert
  use specialinf, only : ellipticK, ellipticE
  use inftools, only : zbrent
  implicit none

  private

  public di_x
  public di_deriv_x, di_deriv_second_x, di_deriv_third_x
  public di_norm_parametric_potential, di_norm_deriv_parametric_potential
  public di_norm_deriv_second_parametric_potential
  public di_norm_deriv_third_parametric_potential

contains

!returns the field value x from k2 by direct integration
  function di_x(k2)
    implicit none
    real(kp) :: di_x
    real(kp), intent(in) :: k2
    real(kp), parameter :: tolint = tolkp

    integer, parameter :: neq = 1
    real(kp), dimension(neq) :: yvar
    real(kp) :: xvar

    xvar = k2
    yvar(1) = 0._kp

    call easydverk(neq,find_di_x,xvar,yvar,1._kp,tolint)

    di_x = yvar(1)

  end function di_x

  subroutine find_di_x(n,k2,y,yprime)
    implicit none          
    integer :: n
    real(kp) :: k2
    real(kp), dimension(n) :: y, yprime

    yprime(1) = -di_deriv_x(k2)

  end subroutine find_di_x



!returns the derivative of the field x with respect to k2; dx/dk2
  function di_deriv_x(k2)
    implicit none
    real(kp) :: di_deriv_x
    real(kp), intent(in) :: k2

    di_deriv_x = -2._kp*sqrt(2._kp)/pi * sqrt(ellipticK(k2) &
         *ellipticK(1._kp-k2))/ k2**1.5_kp

  end function di_deriv_x


!d^2x/dk2^2
  function di_deriv_second_x(k2)
    implicit none
    real(kp) :: di_deriv_second
    real(kp), intent(in) :: k2
    real(kp) :: elKo, elEo, elEp, elKp

    elKo = ellipticK(k2)
    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)

    di_deriv_second_x = (-(elEp*elKo) + elKp*(elEo + elKo*(-7._kp+8._kp*k2))) &
         /(sqrt(2._kp)*sqrt(elKo*elKp)*(-1._kp + k2)*k2**2.5*pi)

  end function di_deriv_second_x

!d^3x/dk2^3
  function di_deriv_third_x(k2)
    implicit none
    real(kp) :: di_deriv_third_x
    real(kp), intent(in) :: k2
    real(kp) :: elKo, elEo, elEp, elKp

    elKo = ellipticK(k2)
    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)

    di_deriv_third_x = (elEo**2*elKp**2 + 2*elEo*elKo*elKp*(elEp + elKp*(7._kp-10._kp*k2)) &
         + elKo**2*(elEp**2 + 2*elEp*elKp*(-9._kp+10._kp*k2) - 3._kp*elKp**2 &
         *(25._kp+8._kp*k2*(-7._kp+4._kp*k2)))) &
         /(4._kp*sqrt(2._kp)*(elKo*elKp)**1.5_kp*(-1._kp+k2)**2*k2**3.5_kp*pi)

  end function di_deriv_third_x



!returns the potential in terms of the parameter k2
  function di_norm_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_parametric_potential
    real(kp), intent(in) :: k2,f

    di_norm_parametric_potential = di_norm_parametric_adspot(k2,f) &
         + di_norm_uplifting(f)

  end function di_norm_parametric_potential



  function di_norm_parametric_adspot(k2,f)
    implicit none
    real(kp) :: di_norm_parametric_adspot
    real(kp), intent(in) :: k2,f
    
    real(kp) :: elKo, elEo, elKp
    real(kp) :: nu

    elKo = ellipticK(k2)
    elEo = ellipticE(k2)

    di_norm_parametric_adspot = 1._kp - 2._kp*(elKo-elEo)/(k2*elKo)

    nu = di_norm_parametric_nu(k2,f)

    if (nu.gt.0._kp) then
           
       elKp = ellipticK(1._kp-k2)

       di_norm_parametric_adspot = di_norm_parametric_adspot &
            - pi/(k2*elKo*elKp) * nu*nu

    endif
 

  end function di_norm_parametric_adspot




  function di_norm_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_parametric_nu
    real(kp), intent(in) :: k2,f

    real(kp) :: elEp, elKp, elKo

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elKo = ellipticK(k2)

    di_norm_parametric_nu = 1._kp - 8._kp*f*sqrt(2._kp)/pi/pi &
         * elKo /sqrt(k2) * (elEp - elKp)**2
    
  end function di_norm_parametric_nu


!returns k2 at which the nu term changes signs, i.e. the monopole
!contribution appears (k2 of the kink)
  function di_k2_nunull(f)
    implicit none
    real(kp) :: di_k2_nunull
    real(kp), intent(in) :: f
    type(transfert) :: diData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: min, max
    
    mini = epsilon(1._kp)
    maxi = 1._kp

    diData%real1 = f
    di_k2_nunull = zbrent(find_di_k2_nunull,mini,maxi,tolFind,diData)

  end function di_k2_nunull



  function find_di_k2_nunull(k2,diData)
    implicit none
    real(kp) :: find_di_k2_nunull
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: f

    f = diData%real1
    find_di_k2_nunull = di_norm_parametric_nu(k2,f)

  end function find_di_k2_nunull



!this is dnu/dk2
  function di_norm_deriv_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_parametric_nu
    real(kp), intent(in) :: k2, f

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_parametric_numterm = (4._kp*Sqrt(2._kp)*(elEp - elKp) &
         *(elEo*elEp - (elEo + 2*elKo*(-1._kp + k2))*elKp)) &
         /(f*(-1._kp + k2)*k2**1.5_kp*pi**2)

  end function di_norm_deriv_parametric_nu


!this is d^2nu/dk2^2
  function di_norm_deriv_second_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_second_parametric_nu
    real(kp), intent(in) :: k2, f

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_second_parametric_numterm =  (2._kp*sqrt(2._kp)*(elEp**2*(-3*elKo*(-1 + k2) &
         - 2*elEo*k2) + 2*elEp*elKp*(elKo + elKo*k2*(-5 + 4*k2) + elEo*(-2 + 4*k2)) &
         + elKp**2*(elEo*(4 - 6*k2) - elKo*(-1 + k2)*(-7 + 10*k2)))) &
         /(f*(-1 + k2)**2*k2**2.5*pi**2)

  end function di_norm_deriv_second_parametric_nu


!this is d^3nu/dk2^3
  function di_norm_deriv_third_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_third_parametric_nu
    real(kp), intent(in) :: k2, f

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_third_parametric_numterm = (sqrt(2._kp)*(elKp**2*(2*elKo*(-1 + k2) &
         * (19 - 51*k2 + 36*k2**2) + elEo*(23 - 65*k2 + 50*k2**2)) &
         + elEp**2*(4*elKo*(-1 + k2)*(-5 + 7*k2) + elEo*(-7 + k2*(7 + 8*k2))) &
         - 2*elEp*elKp*(elKo*(-1 + k2)*(-3 + k2*(-13 + 24*k2)) &
         + elEo*(5 + k2*(-23 + 26*k2)))))/(f*(-1 + k2)**3*k2**3.5*pi**2)

  end function di_norm_deriv_third_parametric_nu




!returns the constant term to be added to the pure ads potential such
!that the minimum is Minkowski
  function di_norm_uplifting(f)
    implicit none
    real(kp) :: di_norm_uplifting
    real(kp), intent(in) :: f

    real(kp), save :: upliftcte = 0._kp
    real(kp), save :: fdone = 0._kp
    
    real(kp) :: k2atmin

    if (f.eq.fdone) then
       di_norm_uplifting = upliftcte
       return
    endif

    k2atmin = di_k2_potmin(f)
    upliftcte = - di_norm_parametric_adspot(k2atmin,f)

    fdone = f
    di_norm_uplifting = upliftcte

  end function di_norm_uplifting


!k2 at which the potential is minimum
  function di_k2_potmin(f)
    implicit none
    real(kp), intent(in) :: f
    type(transfert) :: diData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: mini, maxi

    mini = 0.5_kp
    maxi = 1._kp

    diData%real1 = f
    di_k2_potmin = zbrent(find_di_k2_potmin,mini,maxi,tolFind,diData)
    
  end function di_k2_potmin


  function find_di_k2_potmin(k2,diData)
    implicit none
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: find_di_k2_potmin
    real(kp) :: f

    f = diData%real1

    find_di_k2_potmin = di_norm_deriv_parametric_potential(k2,f)

  end function find_di_k2_potmin




!returns the first derivative of the potential with respect to k2
  function di_norm_deriv_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu
    real(kp) :: elEo, elEp, elKo, elKp

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)
    dnu = di_norm_deriv_parametric_nu(k2,f)

    if (nu.gt.0._kp) then

       di_norm_deriv_parametric_potential = (2*(elEo**2 + (-1._kp + k2)*elKo**2)*elKp**2 & 
            + nu*(-((-(elEp*elKo) + (elEo + elKo)*elKp)*nu) &
            - 4*(-1 + k2)*k2*elKo*elKp*dnu)*pi)/(2*(-1 + k2)*k2**2*elKo**2*elKp**2)
    else

       di_norm_deriv_parametric_potential = (elEo**2 + (-1 + k2)*elKo**2) &
            /((-1 + k2)*k2**2*elKo**2)

    endif

  end function di_norm_deriv_parametric_potential


!returns the 2nd derivative of the potential with respect to k2
  function di_norm_deriv_second_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_second_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu, d2nu
    real(kp) :: elEo, elEp, elKo, elKp

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)
    dnu = di_norm_deriv_parametric_nu(k2,f)
    d2nu = di_norm_deriv_second_parametric(k2,f)

    if (nu.gt.0._kp) then

       di_norm_deriv_second_parametric_potential = (2*elKp**3*(elEo**3 &
            - elEo*elKo**2*(-1 + k2) - 2*elKo**3*(-1 + k2)**2 & 
            - elEo**2*elKo*k2) - (4*dnu**2*elKo**2*elKp**2*(-1 + k2)**2*k2**2 &
            + 4*elKo*elKp*(-1 + k2)*k2*(dnu*elEo*elKp + elKo*(-(dnu*elEp) &
            + elKp*(dnu + d2nu*(-1 + k2)*k2)))*nu + (elEo**2*elKp**2 &
            - elEo*elKo*elKp*(elEp + elKp*(-1 + k2)) + elKo**2*(elEp**2 &
            + elEp*elKp*(-2 + k2) - 2*elKp**2*(-1 + k2)))*nu**2)*Pi) &
            / (2*elKo**3*elKp**3*(-1 + k2)**2*k2**3)

    else

       di_norm_deriv_second_parametric_potential = (elEo**3 - elEo*elKo**2*(-1 + k2)&
            - 2*elKo**3*(-1 + k2)**2 - elEo**2*elKo*k2)/(elKo**3*(-1 + k2)**2*k2**3)

    endif

  end function di_norm_deriv_second_parametric_potential


!returns the 3rd derivative of the potential with respect to k2
  function di_norm_deriv_third_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_third_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu, d2nu, d3nu
    real(kp) :: elEo, elEp, elKo, elKp

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)
    dnu = di_norm_deriv_parametric_nu(k2,f)
    d2nu = di_norm_deriv_second_parametric(k2,f)
    d3nu = di_norm_deriv_third_parametric(k2,f)

    if (nu.gt.0._kp) then

       di_norm_deriv_third_parametric_potential = ((3*elEp*elKo + elKp*(-3*elEo &
            + elKo*(-3 + 4*k2)))*(2*elKp**3*(-elEo**3 + elEo*elKo**2*(-1 + k2) &
            + 2*elKo**3*(-1 + k2)**2 + elEo**2*elKo*k2) &
            + (4*dnu**2*elKo**2*elKp**2*(-1 + k2)**2*k2**2 + 4*elKo*elKp*(-1 + k2)*k2 &
            * (dnu*elEo*elKp + elKo*(-(dnu*elEp) + elKp*(dnu + d2nu*(-1 + k2)*k2)))*nu &
            + (elEo**2*elKp**2 - elEo*elKo*elKp*(elEp + elKp*(-1 + k2)) &
            + elKo**2*(elEp**2 + elEp*elKp*(-2 + k2) - 2*elKp**2*(-1 + k2)))*nu**2)*pi) &
            + elKo*elKp*(2*elKp**2*(3*elEp*(elEo**3 - elEo*elKo**2*(-1 + k2) &
            - 2*elKo**3*(-1 + k2)**2 - elEo**2*elKo*k2) + elKp*(elEo**3*(-3 + k2) &
            + elEo**2*elKo*(1 + 2*k2) + elKo**3*(-1 + k2)**2*(-5 + 4*k2) &
            + elEo*elKo**2*(-1 + k2)*(-7 + 10*k2))) + (-8*dnu*elKo**2*elKp**2*(2*dnu &
            + 3*d2nu*(-1 + k2))*(-1 + k2)**2*k2**3 - 4*elKo*elKp*(-1 + k2)*k2 &
            * (2*elKo*elKp*(2*d2nu + d3nu*(-1 + k2))*(-1 + k2)*k2**2 &
            + dnu*(elKp*(elKo*(3 - 2*k2) + elEo*(-3 + k2)) + elEp*(3*elEo - elKo*k2)))*nu &
            + (elEp**2*elKo*(3*elEo - elKo*k2) + elEp*elKp*(-3*elEo**2 + elKo**2*(1 + k2) &
            + 2*elEo*elKo*(-3 + 2*k2)) - elKp**2*(elEo**2*(-3 + k2) + elEo*elKo*(-2 + 3*k2) &
            + elKo**2*(5 + 4*(-2 + k2)*k2)))*nu**2)*pi))/(4.*elKo**4*elKp**4*(-1 + k2)**3*k2**4)
    else

       di_norm_deriv_third_parametric_potential = (3*elEo**4 - 6*elEo**3*elKo*k2 &
            + elKo**4*(-1 + k2)**2*(-11 + 12*k2) + 4*elEo**2*elKo**2*(1 - k2 + k2**2) &
            + 4*elEo*elKo**3*(1 - 3*k2 + 2*k2**2))/(2.*elKo**4*(-1 + k2)**3*k2**4)
    endif

  end function di_norm_deriv_third_parametric_potential

!WARNING: this is 1/2 (dlnV/dk2)^2, not 1/2 (dlnV/dx)^2
  function di_parametric_epsilon_one(k2,f)
    implicit none
    real(kp) :: di_parametric_epsilon_one
    real(kp), intent(in) :: k2, f
    real(kp) :: pot

    pot = di_norm_parametric_potential(k2,f)

    if (pot.le.0._kp) stop 'di_parametric_epsilon_one: V<=0!'

    di_parametric_epsilon_one = 0.5_kp * (di_norm_deriv_parametric_potential(k2,f) &
         / pot)**2

  end function di_parametric_epsilon_one

    



end module dicommon

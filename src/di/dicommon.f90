module dicommon
  use infprec, only : kp, pi, toldp, tolkp, transfert
  use specialinf, only : ellipticK, ellipticE
  use inftools, only : zbrent, easydverk
  implicit none

  real(kp), parameter :: sqr2 = sqrt(2._kp)
  real(kp), parameter :: ln2 = log(2._kp)
  real(kp), parameter :: ln16 = 4._kp*ln2
  real(kp), parameter :: ln256 = log(256._kp)

  private


  interface di_direct_x
     module procedure di_x
  end interface di_direct_x

  interface di_direct_k2
     module procedure di_k2
  end interface di_direct_k2


  public di_direct_x, di_direct_k2, di_k2_nunull, di_k2_potmin
  public di_deriv_x, di_deriv_second_x, di_deriv_third_x, di_deriv_x_ln
  public di_norm_parametric_potential, di_norm_uplifting
  public di_norm_deriv_second_parametric_potential
  public di_norm_deriv_third_parametric_potential
  public di_norm_deriv_parametric_potential, di_norm_deriv_ln_parametric_potential
  public di_norm_deriv_parametric_potential_regularized
  public di_parametric_epsilon_one, di_parametric_epsilon_two
  public di_parametric_epsilon_three, di_parametric_efold_primitive


contains


!returns the field value x from k2 by direct integration / serie expansion
  function di_x(k2)
    implicit none
    real(kp) :: di_x
    real(kp), intent(in) :: k2
    real(kp), parameter :: tolint = tolkp

    integer, parameter :: neq = 1
    real(kp), dimension(neq) :: yvar
    real(kp) :: xvar
    real(kp), parameter :: k2cutmax = 1._kp-tolkp
    real(kp), parameter :: k2cutmin = tolkp
    real(kp), parameter :: erfisqrtln4 = 2.3111740399905110619206103738822_kp

    if (k2.eq.1._kp) then
       di_x=0._kp
       return
    endif

!integration of first order expansion of di_deriv_x in k2=1
    if (k2.gt.k2cutmax) then
       di_x = (8._kp*sqrt(2._kp*pi)*erfc(sqrt(log(16._kp/(1._kp - k2)))) &
            -(-1._kp + k2)*sqrt(log(256._kp/(-1._kp + k2)**2)))/sqrt(pi)       
       return
    endif

!integration of first order expansion of di_deriv_x in k2=0
    if (k2.lt.k2cutmin) then
       di_x = (erfisqrtln4 - 8._kp*sqrt(log(4._kp)/pi))/2._kp &
            + (log(65536._kp/k2**4)*log(256._kp/k2**2)**2 &
            - 16._kp*(3._kp + 16._kp*ln2**2 + ln16 &
            + log(k2)*(-1._kp - 8._kp*ln2 + log(k2)))) &
            / (sqrt(pi*k2)*log(256._kp/k2**2)**2.5_kp)
       return
    end if


    xvar = log(k2)

    yvar(1) = (8._kp*sqrt(2._kp*pi)*erfc(sqrt(log(16._kp/(1._kp - k2cutmax)))) &
            -(-1._kp + k2cutmax)*sqrt(log(256._kp/(-1._kp + k2cutmax)**2)))/sqrt(pi)

    call easydverk(neq,find_di_x,xvar,yvar,log(k2cutmax),tolint)

    di_x = yvar(1)

  end function di_x

  subroutine find_di_x(n,lnk2,y,yprime,unused)
    implicit none          
    integer :: n
    real(kp) :: lnk2
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: unused

    yprime(1) = -di_deriv_x_ln(lnk2)

  end subroutine find_di_x




!returns k2 from the field value x by inverting direct
!integration. Holy cows, this is going to be slow!
  function di_k2(x)
    implicit none
    real(kp) :: di_k2
    real(kp), intent(in) :: x

    type(transfert) :: diData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: mini, maxi

    if (x.eq.0._kp) then
       di_k2 = 1._kp
       return
    endif

    if (x.lt.0._kp) stop 'di_k2: x < 0!'

    mini = tiny(1._kp)
    maxi = 1._kp

    diData%real1 = x
    di_k2 = zbrent(find_di_k2,mini,maxi,tolFind,diData)

  end function di_k2

  function find_di_k2(k2,diData)
    real(kp) :: find_di_k2
    real(kp), intent(in) :: k2
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: x

    x = diData%real1
    find_di_k2 = di_x(k2) - x

  end function find_di_k2




!returns the derivative of the field x with respect to lnk2; dx/dlnk2,
!more robust against integration than dx/dk2 in k2->0
  function di_deriv_x_ln(lnk2)
    implicit none
    real(kp) :: di_deriv_x_ln
    real(kp), intent(in) :: lnk2
    real(kp) :: k2

    k2 = exp(lnk2)

    di_deriv_x_ln = -2._kp/pi * sqrt(2._kp*ellipticK(k2) &
         *ellipticK(1._kp-k2)/k2)

  end function di_deriv_x_ln



!returns the derivative of the field x with respect to k2; dx/dk2
  function di_deriv_x(k2)
    implicit none
    real(kp) :: di_deriv_x
    real(kp), intent(in) :: k2

!for k2-->1 di_deriv_x --> -sqrt(log(256._kp/(1._kp-k2)**2)/pi)

    if (k2.eq.1._kp) stop 'di_deriv_x: Infinity in k2=1!'

    di_deriv_x = -2._kp/pi * sqrt(2._kp*ellipticK(k2) &
         *ellipticK(1._kp-k2))/ k2**1.5_kp

  end function di_deriv_x


!d^2x/dk2^2
  function di_deriv_second_x(k2)
    implicit none
    real(kp) :: di_deriv_second_x
    real(kp), intent(in) :: k2
    real(kp) :: elKo, elEo, elEp, elKp

    if (k2.eq.1._kp) stop 'di_deriv_second_x: Infinity in k2=1!'

    elKo = ellipticK(k2)
    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)

    di_deriv_second_x = (-(elEp*elKo) + elKp*(elEo + elKo*(-7._kp+8._kp*k2))) &
         /(sqrt(2._kp*elKo*elKp)*(-1._kp + k2)*k2**2.5*pi)

  end function di_deriv_second_x

!d^3x/dk2^3
  function di_deriv_third_x(k2)
    implicit none
    real(kp) :: di_deriv_third_x
    real(kp), intent(in) :: k2
    real(kp) :: elKo, elEo, elEp, elKp

    if (k2.eq.1._kp) stop 'di_deriv_third_x: Infinity in k2=1!'

    elKo = ellipticK(k2)
    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)

    di_deriv_third_x = (elEo**2*elKp**2 + 2*elEo*elKo*elKp*(elEp + elKp*(7._kp-10._kp*k2)) &
         + elKo**2*(elEp**2 + 2*elEp*elKp*(-9._kp+10._kp*k2) - 3._kp*elKp**2 &
         *(25._kp+8._kp*k2*(-7._kp+4._kp*k2)))) &
         /(4._kp*sqr2*(elKo*elKp)**1.5_kp*(-1._kp+k2)**2*k2**3.5_kp*pi)

  end function di_deriv_third_x



!returns the potential in terms of the parameter k2
  function di_norm_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_parametric_potential
    real(kp), intent(in) :: k2,f
    real(kp) :: ppot

     ppot = di_norm_parametric_adspot(k2,f) &
         + di_norm_uplifting(f)

     if (ppot.lt.-10._kp*epsilon(1._kp)) then
        stop 'di_norm_parametric_potential < 0!'
     elseif (ppot.lt.0._kp) then
        ppot = 0._kp
     endif

    di_norm_parametric_potential = ppot


  end function di_norm_parametric_potential



  function di_norm_parametric_adspot(k2,f)
    implicit none
    real(kp) :: di_norm_parametric_adspot
    real(kp), intent(in) :: k2,f
    
    real(kp) :: elKo, elEo, elKp
    real(kp) :: nu

    if (k2.eq.1._kp) then
       di_norm_parametric_adspot = -1._kp
       return
    elseif (k2.eq.0._kp) then
       di_norm_parametric_adspot = 0._kp
       return
    endif

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

    di_norm_parametric_nu = 1._kp - 8._kp*sqr2/pi/pi &
         * elKo /sqrt(k2)/f * (elEp - elKp)**2
    
  end function di_norm_parametric_nu


!returns k2 at which the nu term changes signs, i.e. the monopole
!contribution appears (k2 of the kink)
  function di_k2_nunull(f)
    implicit none
    real(kp) :: di_k2_nunull
    real(kp), intent(in) :: f
    type(transfert) :: diData
    real(kp), parameter :: tolFind = tolkp
    real(kp) :: mini, maxi

    real(kp), save :: fdone = -1._kp
    real(kp), save :: k2null = 1._kp
!$omp threadprivate(fdone,k2null)
    
    if (f.eq.fdone) then
       di_k2_nunull = k2null
       return
    endif
       

    mini = epsilon(1._kp)
    maxi = 1._kp

    diData%real1 = f
    di_k2_nunull = zbrent(find_di_k2_nunull,mini,maxi,tolFind,diData)

    fdone = f
    k2null = di_k2_nunull

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
    real(kp) :: elEp, elKp, elEo, elKo

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_parametric_nu = (4._kp*sqr2*(elEp - elKp) &
         *(elEo*elEp - (elEo + 2*elKo*(-1._kp + k2))*elKp)) &
         /(f*(-1._kp + k2)*k2**1.5_kp*pi**2)

  end function di_norm_deriv_parametric_nu


!this is d^2nu/dk2^2
  function di_norm_deriv_second_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_second_parametric_nu
    real(kp), intent(in) :: k2, f
    real(kp) :: elEp, elKp, elEo, elKo

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_second_parametric_nu =  (2._kp*sqr2*(elEp**2*(-3*elKo*(-1 + k2) &
         - 2*elEo*k2) + 2*elEp*elKp*(elKo + elKo*k2*(-5 + 4*k2) + elEo*(-2 + 4*k2)) &
         + elKp**2*(elEo*(4 - 6*k2) - elKo*(-1 + k2)*(-7 + 10*k2)))) &
         /(f*(-1 + k2)**2*k2**2.5*pi**2)

  end function di_norm_deriv_second_parametric_nu


!this is d^3nu/dk2^3
  function di_norm_deriv_third_parametric_nu(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_third_parametric_nu
    real(kp), intent(in) :: k2, f
    real(kp) :: elEp, elKp, elEo, elKo

    elEp = ellipticE(1._kp-k2)
    elKp = ellipticK(1._kp-k2)
    elEo = ellipticE(k2)
    elKo = ellipticK(k2)
    
    di_norm_deriv_third_parametric_nu = (sqr2*(elKp**2*(2*elKo*(-1 + k2) &
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
    real(kp), save :: fdone = -1._kp
!$omp threadprivate(upliftcte,fdone)

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
    real(kp) :: di_k2_potmin
    real(kp), intent(in) :: f
    type(transfert) :: diData
    real(kp), parameter :: tolFind = epsilon(1._kp)
    real(kp) :: mini, maxi

    mini = (1._kp + epsilon(1._kp))*di_k2_nunull(f)
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


!return the logarithmic derivative of the parametric potential w.r.t k2; i.e. dV/V/dk2
  function di_norm_deriv_ln_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_ln_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp), parameter :: epsCut = tolkp**(1._kp/3._kp)

    real(kp) :: nu, dnu, ucte, mk2m1
    real(kp) :: elEo, elEp, elKo, elKp

    ucte = di_norm_uplifting(f)

!analytic extension
    if (k2.eq.1._kp) then

       di_norm_deriv_ln_parametric_potential = 1._kp/(-1._kp + ucte)
       return

    elseif (1._kp-k2.lt.epsCut) then

       mk2m1 = 1._kp-k2

       di_norm_deriv_ln_parametric_potential = 1._kp/(-1._kp + ucte) &
            + (mk2m1*(-3._kp - (16._kp*sqr2*(-1._kp + ucte))/f + 7._kp*ucte)) &
            / (4._kp*(-1._kp + ucte)**2) + (mk2m1**2*(-96._kp*sqr2*(-1._kp + ucte) &
            *(-1._kp + 2._kp*ucte) + f*(13._kp - 36._kp*ucte + 39._kp*ucte**2))) &
            /(16._kp*f*(-1._kp + ucte)**3) + (mk2m1**3*(-64._kp*sqr2*f*(-1._kp + ucte) &
            *(35._kp + ucte*(-98._kp + 95._kp*ucte)) + f**2*(-207._kp + ucte*(821._kp &
            + ucte*(-1149._kp + 791._kp*ucte))) + 128._kp*(-1._kp + ucte)**2 &
            *(33._kp - ucte + 16._kp*(-1 + ucte)*ln2 &
            - 4._kp*(-1._kp + ucte)*log(mk2m1))))/(256._kp*f**2*(-1._kp + ucte)**4)

       return

    endif
    
    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)


    if (nu.gt.0._kp) then

       dnu = di_norm_deriv_parametric_nu(k2,f)

       di_norm_deriv_ln_parametric_potential = (2*elKp**2*(elEo**2 + elKo**2*(-1._kp + k2)) &
            -  nu*(elEo*elKp*nu + elKo*(-(elEp*nu) + elKp*(4*dnu*(-1._kp + k2)*k2 + nu)))*pi) &
            / (2._kp*elKo*elKp*(-1._kp + k2)*k2*(-(nu**2*Pi) &
            + elKp*(2*elEo + elKo*(-2._kp + k2 + k2*ucte))))
    else

       di_norm_deriv_ln_parametric_potential = (elEo**2 + elKo**2*(-1._kp + k2)) &
            /(elKo*(-1 + k2)*k2*(2*elEo + elKo*(-2._kp + k2 + k2*ucte)))

    endif
    

  end function di_norm_deriv_ln_parametric_potential



!returns the first derivative of the potential with respect to k2
  function di_norm_deriv_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp), parameter :: epsCut = tolkp**(1._kp/3._kp)
    
    real(kp) :: nu, dnu, mk2m1
    real(kp) :: elEo, elEp, elKo, elKp

!analytic extension
    if (k2.eq.1._kp) then
       di_norm_deriv_parametric_potential = 1._kp
       return

    elseif (1._kp-k2.lt.epsCut) then

       mk2m1 = 1._kp-k2

       di_norm_deriv_parametric_potential =  1._kp + (1.75_kp - (4._kp*sqr2)/f)*mk2m1 + (2.4375_kp &
            - (12._kp*sqr2)/f)*mk2m1**2 + (mk2m1**4*(-79960._kp*sqrt(2._kp)*f + 7615._kp*f**2 &
            + 512._kp*(-12._kp + 140._kp*ln2 - 35._kp*log(mk2m1)))) &
            / (2048._kp*f**2) + mk2m1**3*(3.08984375_kp + (-2._kp - 95._kp*sqr2*f &
            + 32._kp*ln2)/(4._kp*f**2) - (2._kp*log(mk2m1))/f**2)
       return
    end if

    di_norm_deriv_parametric_potential = di_norm_deriv_parametric_potential_regularized(k2,f) &
         /k2/k2/(1._kp-k2)

  end function di_norm_deriv_parametric_potential



!returns the first derivative of the potential with respect to k2 and
!multiplied k2^2 x (1-k2) [and explicit singularity at k2=0 and k2=1]
  function di_norm_deriv_parametric_potential_regularized(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_parametric_potential_regularized
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu
    real(kp) :: elEo, elEp, elKo, elKp

!analytic extension
    if (k2.eq.1._kp) then
       stop 'di_norm_deriv_parametric_potential_regularized: do not use in k2=1!'
    endif

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)

    if (nu.gt.0._kp) then

       dnu = di_norm_deriv_parametric_nu(k2,f)

       di_norm_deriv_parametric_potential_regularized = (2*(elEo**2 + (-1._kp + k2)*elKo**2)*elKp**2 & 
            + nu*(-((-(elEp*elKo) + (elEo + elKo)*elKp)*nu) &
            - 4*(-1 + k2)*k2*elKo*elKp*dnu)*pi)/(-2*elKo**2*elKp**2)
    else

       di_norm_deriv_parametric_potential_regularized = (elEo**2 + (-1 + k2)*elKo**2) &
            /(-elKo**2)

    endif

  end function di_norm_deriv_parametric_potential_regularized



!returns the 2nd derivative of the potential with respect to k2
  function di_norm_deriv_second_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_second_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu, d2nu
    real(kp) :: elEo, elEp, elKo, elKp

!analytic extension
    if (k2.eq.1._kp) then
       di_norm_deriv_second_parametric_potential = 4._kp*sqr2/f - 7._kp/4._kp
       return
    endif

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)    

    if (nu.gt.0._kp) then

       dnu = di_norm_deriv_parametric_nu(k2,f)
       d2nu = di_norm_deriv_second_parametric_nu(k2,f)

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
  recursive function di_norm_deriv_third_parametric_potential(k2,f)
    implicit none
    real(kp) :: di_norm_deriv_third_parametric_potential
    real(kp), intent(in) :: k2,f

    real(kp) :: nu, dnu, d2nu, d3nu
    real(kp) :: elEo, elEp, elKo, elKp

!analytic extension
    if (k2.eq.1._kp) then
       di_norm_deriv_third_parametric_potential = -24._kp*sqr2/f + 39._kp/8._kp
       return
    endif

    elEo = ellipticE(k2)
    elEp = ellipticE(1._kp - k2)
    elKo = ellipticK(k2)
    elKp = ellipticK(1._kp - k2)

    nu = di_norm_parametric_nu(k2,f)    

    if (nu.gt.0._kp) then
       
       dnu = di_norm_deriv_parametric_nu(k2,f)
       d2nu = di_norm_deriv_second_parametric_nu(k2,f)
       d3nu = di_norm_deriv_third_parametric_nu(k2,f)

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
            + elKo**2*(5 + 4*(-2 + k2)*k2)))*nu**2)*pi))/(4._kp*elKo**4*elKp**4*(-1 + k2)**3*k2**4)

    else
       di_norm_deriv_third_parametric_potential = (3*elEo**4 - 6*elEo**3*elKo*k2 &
            + elKo**4*(-1 + k2)**2*(-11 + 12*k2) + 4*elEo**2*elKo**2*(1 - k2 + k2**2) &
            + 4*elEo*elKo**3*(1 - 3*k2 + 2*k2**2))/(2._kp*elKo**4*(-1 + k2)**3*k2**4)
    endif

  end function di_norm_deriv_third_parametric_potential


!eps1 in terms of k2 and multiplied by lambda^2
  function di_parametric_epsilon_one(k2,f)
    real(kp) :: di_parametric_epsilon_one
    real(kp), intent(in) :: k2, f
    real(kp) :: dx, dlnV, dxk2, dV, V

!    dlnV = di_norm_deriv_ln_parametric_potential(k2,f)
!    dx = di_deriv_x(k2)

    dxk2 = di_deriv_x_ln(log(k2))
    dV = di_norm_deriv_parametric_potential(k2,f)
    V = di_norm_parametric_potential(k2,f) + epsilon(1._kp)
    
!    di_parametric_epsilon_one = 0.5_kp* k2 * k2 * (dlnV/dxk2)**2

    di_parametric_epsilon_one = 0.5_kp* k2 * k2 * (dV/V/dxk2)**2

  end function di_parametric_epsilon_one


!eps2 in terms of k2 and multiplied by lambda^2
  function di_parametric_epsilon_two(k2,f)
    implicit none
    real(kp) :: di_parametric_epsilon_two
    real(kp), intent(in) :: k2,f

    real(kp) :: V, dV, d2V, dx, d2x

    dx = di_deriv_x(k2)
    d2x = di_deriv_second_x(k2)
    V = di_norm_parametric_potential(k2,f)
    dV = di_norm_deriv_parametric_potential(k2,f)
    d2V = di_norm_deriv_second_parametric_potential(k2,f)

    di_parametric_epsilon_two = (2*d2x*dV*V + 2*dx*(dV**2 - d2V*V)) &
         /(dx**3*V**2)

  end function di_parametric_epsilon_two


!eps2 in terms of k2 and multiplied by lambda^2
  function di_parametric_epsilon_three(k2,f)
    implicit none
    real(kp) :: di_parametric_epsilon_three
    real(kp), intent(in) :: k2,f

    real(kp) :: V, dV, d2V, d3V, dx, d2x, d3x

    dx = di_deriv_x(k2)
    d2x = di_deriv_second_x(k2)
    d3x = di_deriv_third_x(k2)
    V = di_norm_parametric_potential(k2,f)
    dV = di_norm_deriv_parametric_potential(k2,f)
    d2V = di_norm_deriv_second_parametric_potential(k2,f)
    d3V = di_norm_deriv_third_parametric_potential(k2,f)

    di_parametric_epsilon_three = (dV*(-2*dV**3*dx**2 - 3*d2x*dV**2*dx*V &
         + dx*(3*d2V*d2x - d3V*dx)*V**2 &
         + dV*V*(3*d2V*dx**2 - 3*d2x**2*V + d3x*dx*V))) &
         / (dx**3*V**2*(-(d2x*dV*V) + dx*(-dV**2 + d2V*V)))

  end function di_parametric_epsilon_three




!efold primitive in terms of the parameter k2 and divided by
!lambda^2. That beast becomes wild at k2=0 and k2=1 such that I had to
!pull out some numacc sticks to keep it quiet in its hole.
  function di_parametric_efold_primitive(k2,f)
    implicit none
    real(kp) :: di_parametric_efold_primitive
    real(kp), intent(in) :: k2, f
    type(transfert) :: diData

!avoid prohibitive integration time in QUADPREC
    real(kp), parameter :: tolint=max(tolkp,toldp)
!the expansion are done at O(x^3), not need to go further than    
    real(kp), parameter :: epsCut = tolkp**(1._kp/3._kp)
    real(kp), parameter :: lnepsCut = log(epsCut)
    real(kp), parameter :: ln1mepsCut = log(1._kp-epsCut)

    real(kp), save :: k2null = -1._kp
    real(kp), save :: fdone = -1._kp
    real(kp), save :: ucte = -1._kp
!$omp threadprivate(k2null,fdone,ucte)
    real(kp), save :: ycutAtZero = 0._kp
    real(kp), save :: ycutAtOne = 0._kp
!$omp threadprivate(ycutAtZero,ycutAtOne)
    integer, parameter :: neq = 1
    real(kp) :: xvar ,lnk2
    real(kp), dimension(neq) :: yvar


    if (f.ne.fdone) then
       fdone = f      
       ucte = di_norm_uplifting(f)
       ycutAtZero = di_numacc_parametric_efold_primitive_atzero(epsCut,ucte)
       ycutAtOne = di_numacc_parametric_efold_primitive_atone(1._kp-epsCut,ucte,f)
    endif

!integration w.r.t lnk2, faster and more robust
    lnk2 = log(k2)

!this set the overall constant such that the expansion in k2=1 and the
!exact integration match to the best
    xvar = ln1mepsCut
    yvar(1) = ycutAtOne

    diData%real1 = f
   
    if ((lnk2.gt.lnepsCut).and.(lnk2.lt.ln1mepsCut)) then

       call easydverk(neq,find_di_parametric_efold_primitive,xvar,yvar,lnk2,tolint,diData)

    elseif (lnk2.le.lnepsCut) then

       call easydverk(neq,find_di_parametric_efold_primitive,xvar,yvar,lnepsCut,tolint,diData)

       yvar = yvar + di_numacc_parametric_efold_primitive_atzero(k2,ucte) &
            - ycutAtZero

    elseif (lnk2.ge.ln1mepsCut) then
       
       yvar = di_numacc_parametric_efold_primitive_atone(k2,ucte,f)

    else
       stop 'di_parametric_efold_primitive: internal error!'

    endif

    di_parametric_efold_primitive = yvar(1)
     
  end function di_parametric_efold_primitive


  subroutine find_di_parametric_efold_primitive(n,lnk2,y,yprime,diData)
    implicit none
    integer :: n
    real(kp) :: lnk2
    real(kp), dimension(n) :: y, yprime
    type(transfert), optional, intent(inout) :: diData
    real(kp) :: f, V, dV, dxk2, k2

    f = diData%real1

    k2 = exp(lnk2)

    V = di_norm_parametric_potential(k2,f)
!this is dV/dx * k2^2 * (1-k2)
    dV = di_norm_deriv_parametric_potential_regularized(k2,f)
!this is dx/dlnk2 = dx/dk2 * k2
    dxk2 = di_deriv_x_ln(lnk2)
    
    yprime = (dxk2/(dV/k2)) *V * dxk2  * (1._kp - k2)

  end subroutine find_di_parametric_efold_primitive



!analytic integration using an Laurent expansion in k2 towards k2=0
  function di_numacc_parametric_efold_primitive_atzero(k2,ucte)
    implicit none
    real(kp) :: di_numacc_parametric_efold_primitive_atzero
    real(kp), intent(in) :: k2, ucte
    real(kp) :: lnk2
    real(kp), parameter :: ln512 = ln256 + ln2

    lnk2 = log(k2)

    di_numacc_parametric_efold_primitive_atzero = -((4._kp*(1._kp - 8._kp*ln2 + 2._kp*lnk2)*ucte) &
         /(k2**2*pi) + (-2._kp + ln256 + 32*ln2*ucte - 2._kp*lnk2*(1._kp + 4*ucte))/(k2*pi) &
         + (k2**2*(736._kp + 36._kp*ln2*(32._kp - 203._kp*ucte) + 589._kp*ucte &
         + 9._kp*lnk2*(-32._kp + 203._kp*ucte)))/(6144._kp*pi) &
         + (lnk2*(15._kp*lnk2*ucte + 4._kp*(4._kp + (3._kp - 30._kp*ln2)*ucte)))/(16._kp*pi) &
         + (k2*(lnk2*(-9._kp + 90._kp*ucte) &
         + 4._kp*(12._kp + ln512 - (7._kp + 90._kp*ln2)*ucte)))/(96._kp*pi))
    
  end function di_numacc_parametric_efold_primitive_atzero



!analytic integration using an Laurent expansion in 1-k2 towards k2=1
!(implicitely assume that the monopole terms are on: never call this
!function with k2 < k2null.
  function di_numacc_parametric_efold_primitive_atone(k2,ucte,f)
    implicit none
    real(kp) :: di_numacc_parametric_efold_primitive_atone
    real(kp), intent(in) :: k2, ucte, f
    real(kp) :: mk2m1, lnmk2m1
    real(kp), parameter :: ln4096 = 12._kp*ln2

    real(kp), save :: fdone = -1._kp
    real(kp), save :: k2null = -1._kp
!$omp threadprivate(fdone,k2null)
    if (f.ne.fdone) then
       k2null = di_k2_nunull(f)
       fdone = f
    endif

    if (k2.lt.k2null) then
       write(*,*)'di_numacc_parametric_efold_primitive_atone: k2 < k2null!'
       write(*,*)'k2= k2null= ',k2,k2null
    endif

!minus(k2 minus 1) = 1-k2
    mk2m1 = 1._kp - k2
    lnmk2m1 = log(mk2m1)    

    di_numacc_parametric_efold_primitive_atone = (-2*(1._kp + ln16 - lnmk2m1)*mk2m1*(-1._kp + ucte))/pi &
         - (mk2m1**2*(16._kp*(1._kp + ln256 - 2._kp*lnmk2m1)*sqr2*(-1._kp + ucte) &
         + f*(-7._kp + lnmk2m1*(22._kp - 14._kp*ucte) + 3._kp*ucte &
         + 8._kp*ln2*(-11._kp + 7._kp*ucte))))/(8.*f*pi) - (mk2m1**3*(2048._kp*(1 + ln4096 &
         - 3._kp*lnmk2m1)*(-1._kp + ucte) + 128*f*sqr2*(-4._kp + 3*lnmk2m1*(7._kp - 6*ucte) &
         + 3._kp*ucte + 12._kp*ln2*(-7._kp + 6._kp*ucte)) + 3._kp*f**2*(-5._kp &
         + 6._kp*lnmk2m1*(53._kp - 25._kp*ucte) - 19._kp*ucte &
         + 24._kp*ln2*(-53._kp + 25._kp*ucte))))/(288._kp*f**2*pi)


  end function di_numacc_parametric_efold_primitive_atone

end module dicommon

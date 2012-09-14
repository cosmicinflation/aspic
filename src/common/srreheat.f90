!generic reheating functions assuming slow-roll evolution and a
!mean values for wreh
module srreheat
  use infprec, only : kp,pi
  use cosmopar, only : lnMpcToKappa, HubbleSquareRootOf2OmegaRad
  use cosmopar, only : QrmsOverT, powerAmpScalar, lnMpinGeV
  implicit none

  private
  

  real(kp), parameter :: Nzero = log(kstar) - lnMpcToKappa &
       - 0.5_kp*log(sqrt(1.5_kp)*HubbleSquareRootOf2OmegaRad)
   
  
  logical, parameter :: display = .false.
    
  public display, pi, Nzero
  public slowroll_validity, get_calfconst
  public quadrupole_to_primscalar, find_reheat
  public ln_rho_endinf, ln_rho_reheat,log_energy_reheat_ingev

contains


  function slowroll_validity(eps1,eps2)
    implicit none
    real(kp), intent(in) :: eps1
    real(kp), intent(in), optional :: eps2
    logical :: slowroll_validity

    slowroll_validity = .true.
    
    if ((eps1.lt.0._kp).or.(eps1.gt.3._kp)) then
       write(*,*)'epsilon_1= ',eps1
       write(*,*) 'slowroll_validity: slow-roll is inconsistent'
       stop
    elseif (eps1.gt.1_kp) then       
       slowroll_validity = .false.
    elseif (eps1.lt.epsilon(1._kp)) then
       write(*,*)
       write(*,*)'epsilon_1 < numaccuracy!',eps1,epsilon(1._kp)
       write(*,*)
    endif

    if (present(eps2)) then
       if (abs(eps2).gt.1._kp) then          
          slowroll_validity = .false.
       endif
    endif

  end function slowroll_validity


  function get_calfconst(lnRhoReh,Pstar,w,epsEnd,potEnd)
    implicit none
    real(kp) :: get_calfconst
    real(kp), intent(in) :: w,Pstar, epsEnd,potEnd,lnRhoReh

    real(kp) :: lnHoverSqrtEps,calF

    lnHoverSqrtEps = 0.5_kp*log(Pstar*8*pi*pi)

    calF = -Nzero + (1._kp+3._kp*w)/(3._kp+3._kp*w)*lnHoverSqrtEps &
         - 1._kp/(3._kp+3._kp*w)*log(9._kp*2._kp**(1.5*w+0.5)*potEnd/(3._kp-epsEnd)) &
         +(1._kp-3._kp*w)/(12._kp+12._kp*w)*lnRhoReh
        
    get_calfconst = calF


  end function get_calfconst




  function quadrupole_to_primscalar(QoverT)
    implicit none
    real(kp) :: quadrupole_to_primscalar
    real(kp), intent(in) :: QoverT
    real(kp) :: H2OverEpsOneOverPi2,Pstar

    H2OverEpsOneOverPi2 = 480._kp*(QoverT)**2

    Pstar = H2OverEpsoneOverPi2/8._kp

    quadrupole_to_primscalar = Pstar

  end function quadrupole_to_primscalar



  function find_reheat(nuStar,calFplusNuEnd,w,epsStar,Vstar)
    implicit none
    real(kp) :: find_reheat
    real(kp), intent(in) :: nuStar, calFplusNuEnd
    real(kp), intent(in) :: w, epsStar, Vstar

    find_reheat = nuStar - calFplusNuEnd + 1._kp/(3._kp + 3._kp*w) &
         * log( (9._kp-3._kp*epsStar)/( 9._kp*(2._kp*epsStar)**(0.5_kp+1.5_kp*w) &
         *Vstar) )

  end function find_reheat



  function ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_endinf
    real(kp), intent(in) :: Pstar, epsOneStar, epsOneEnd
    real(kp), intent(in) :: VendOverVstar

    real(kp) :: H2OverEps

    H2OverEps = Pstar*8*pi*pi

    ln_rho_endinf = log(H2OverEps) + log(3._kp*epsOneStar*(3._kp - epsOneStar) &
         /(3._kp - epsOneEnd)) + log(VendOverVstar)

  end function ln_rho_endinf

 
  function ln_rho_reheat(w,Pstar,epsOneStar,epsOneEnd,deltaNstar,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_reheat
    real(kp), intent(in) :: w, Pstar, epsOneStar, epsOneEnd, deltaNstar
    real(kp), intent(in) :: VendOverVstar

    real(kp) :: lnH2OverEps, oneMinusThreeWlnEreh

     if (w.eq.1._kp/3._kp) then
        write(*,*)'ln_rho_reheat: w = 1/3!'
        stop
     end if

     lnH2OverEps = log(Pstar*8*pi*pi)


     oneMinusThreeWlnEreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsOneStar) &
         + log(3._kp*(3._kp - epsOneStar)/(epsOneStar*(3._kp - epsOneEnd))) &
         + log(VendOverVstar)

     ln_rho_reheat = oneMinusThreeWlnEreh*4._kp/(1._kp - 3._kp*w)

  end function ln_rho_reheat

  function log_energy_reheat_ingev(lnRhoReh)

    implicit none
    real(kp) :: log_energy_reheat_ingev
    real(kp), intent(in) :: lnRhoReh

    log_energy_reheat_ingev=0.25_kp*(lnRhoReh + 4._kp*lnMpinGeV)/log(10._kp)

  end function log_energy_reheat_ingev











end module srreheat


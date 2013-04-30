!generic reheating functions assuming slow-roll evolution and a mean
!values for wreh or input lnRrad
module srreheat
  use infprec, only : kp,pi
  use srflow, only : slowroll_violated, inverse_slowroll_corrections
  use cosmopar, only : lnMpcToKappa, HubbleSquareRootOf3OmegaRad
  use cosmopar, only : QrmsOverT, kstar, powerAmpScalar, lnMpinGeV
  use cosmopar, only : RelatDofRatio
  implicit none

  private

!you certainly want to use all _leadorder functions only. _anyorder
!functions are meant for accuracy test by consistently adding first
!order slow-roll corrections that change Reheating related quantities
!by typically Delta N* ~ 0.01

  interface find_reheat
     module procedure find_reheat_rhow_leadorder
     module procedure find_reheat_rhow_anyorder
  end interface find_reheat

  interface find_reheat_rrad
     module procedure find_reheat_rrad_leadorder
     module procedure find_reheat_rrad_anyorder
  end interface find_reheat_rrad

  interface find_reheat_rreh
     module procedure find_reheat_rreh_leadorder
     module procedure find_reheat_rreh_anyorder
  end interface find_reheat_rreh

  interface get_calfconst
     module procedure get_calfconst_rhow
  end interface get_calfconst

  interface ln_rho_endinf
     module procedure ln_rho_endinf_leadorder
     module procedure ln_rho_endinf_anyorder
  end interface ln_rho_endinf

  interface ln_rho_reheat
     module procedure ln_rho_reheat_leadorder
     module procedure ln_rho_reheat_anyorder
  end interface ln_rho_reheat
  
  interface ln_potential_normalization
     module procedure ln_potential_normalization_leadorder
     module procedure ln_potential_normalization_anyorder
  end interface ln_potential_normalization


!  real(kp), parameter :: Nzero = log(kstar) - lnMpcToKappa &
!       - 0.5_kp*log(sqrt(1.5_kp)*HubbleSquareRootOf2OmegaRad)
   
  real(kp), parameter :: Nzero = log(kstar) - lnMpcToKappa &
       - 0.5_kp*log(HubbleSquareRootOf3OmegaRad) &
       -0.25_kp*log(RelatDofRatio)


!calculations are at this slow-roll order accuracy max
  integer, parameter :: norder = 1
  integer, parameter :: neps = norder + 1

  
  logical, parameter :: display = .false.
    
  public display, pi, Nzero
  public slowroll_validity, get_calfconst
  public quadrupole_to_primscalar, primscalar_to_quadrupole, find_reheat
  public ln_rho_endinf, ln_rho_reheat,log_energy_reheat_ingev
  public ln_potential_normalization

  public get_calfconst_rrad, find_reheat_rrad
  public get_calfconst_rreh, find_reheat_rreh
  public get_lnrrad_rreh, get_lnrreh_rrad
  public get_lnrrad_rhow, get_lnrreh_rhow

contains


  function slowroll_validity(eps1,eps2)
    implicit none
    real(kp), intent(in) :: eps1
    real(kp), intent(in), optional :: eps2
    logical :: slowroll_validity

    if (eps1.lt.epsilon(1._kp)) then
       write(*,*)
       write(*,*)'epsilon_1 < numaccuracy!',eps1,epsilon(1._kp)
       write(*,*)
    endif

    if (eps1.lt.0._kp) then
       stop 'Congratulations: eps1 < 0!'
    endif

    if (present(eps2)) then
       slowroll_validity = .not.slowroll_violated((/eps1,eps2/))
    else
       slowroll_validity = .not.slowroll_violated((/eps1/))
    endif

  end function slowroll_validity

  

  function find_reheat_rhow_leadorder(nuStar,calFplusNuEnd,w,epsStar,Vstar)
    implicit none
    real(kp) :: find_reheat_rhow_leadorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd
    real(kp), intent(in) :: w, epsStar, Vstar

    find_reheat_rhow_leadorder = nuStar - calFplusNuEnd + 1._kp/(3._kp + 3._kp*w) &
         * log( (9._kp)/( 9._kp*(2._kp*epsStar)**(0.5_kp+1.5_kp*w) &
         *Vstar) )

  end function find_reheat_rhow_leadorder

  
!for testing accuracy, add slow-roll corrections. Not robust if called with eps~1
  function find_reheat_rhow_anyorder(nuStar,calFplusNuEnd,w,epsStarVec,Vstar)
    implicit none
    real(kp) :: find_reheat_rhow_anyorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd
    real(kp), intent(in) :: w, Vstar
    real(kp), dimension(neps), intent(in) :: epsStarVec

    find_reheat_rhow_anyorder = nuStar - calFplusNuEnd + 1._kp/(3._kp + 3._kp*w) &
         * log( (9._kp - 3._kp*epsStarVec(1)) &
         /( 9._kp*(2._kp*epsStarVec(1))**(0.5_kp+1.5_kp*w) &
         *Vstar) ) - 0.25_kp*log(inverse_slowroll_corrections(epsStarVec))

  end function find_reheat_rhow_anyorder



  function get_calfconst_rhow(lnRhoReh,Pstar,w,epsEnd,potEnd)
    implicit none
    real(kp) :: get_calfconst_rhow
    real(kp), intent(in) :: lnRhoReh,Pstar,w,epsEnd,potEnd

    real(kp) :: cmbMeasuredHalf

    cmbMeasuredHalf = 0.5_kp*log(Pstar*8._kp*pi**2)

    get_calfconst_rhow = -Nzero + (1._kp+3._kp*w)/(3._kp+3._kp*w)*cmbMeasuredHalf &
         - 1._kp/(3._kp+3._kp*w)*log(9._kp*2._kp**(1.5_kp*w+0.5_kp) &
         *potEnd/(3._kp-epsEnd)) + (1._kp-3._kp*w)/(12._kp+12._kp*w)*lnRhoReh
        
  end function get_calfconst_rhow



  function find_reheat_rrad_leadorder(nuStar,calFplusNuEnd,epsOneStar,Vstar)
    implicit none
    real(kp) :: find_reheat_rrad_leadorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd, Vstar
    real(kp), intent(in) :: epsOneStar

    find_reheat_rrad_leadorder = nuStar - calFplusNuEnd &
         + 0.25_kp*log(9._kp/(epsOneStar*Vstar))

  end function find_reheat_rrad_leadorder



  function find_reheat_rrad_anyorder(nuStar,calFplusNuEnd,epsStarVec,Vstar)
    implicit none
    real(kp) :: find_reheat_rrad_anyorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd, Vstar
    real(kp), intent(in), dimension(neps) :: epsStarVec

    find_reheat_rrad_anyorder = nuStar - calFplusNuEnd &
         + 0.25_kp*log((9._kp-3._kp*epsStarVec(1))/(epsStarVec(1)*Vstar)) &
         - 0.25_kp*log(inverse_slowroll_corrections(epsStarVec))

  end function find_reheat_rrad_anyorder


  
  function get_calfconst_rrad(lnRrad,Pstar,epsEnd,potEnd)
    implicit none
    real(kp) :: get_calfconst_rrad
    real(kp), intent(in) :: lnRrad,Pstar,epsEnd,potEnd

    real(kp) :: cmbMeasuredQuart

    cmbMeasuredQuart = 0.25_kp*log(Pstar*8._kp*pi**2)

    get_calfconst_rrad = -Nzero + cmbMeasuredQuart &
         - 0.25_kp*log(potEnd/(3._kp-epsEnd)) + lnRrad
            
  end function get_calfconst_rrad
  

  function find_reheat_rreh_leadorder(nuStar,calFplusNuEnd,Vstar)
    implicit none
    real(kp) :: find_reheat_rreh_leadorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd, Vstar    

    find_reheat_rreh_leadorder = nuStar - calFplusNuEnd &
         + 0.5_kp*log((9._kp)/Vstar)

  end function find_reheat_rreh_leadorder


  function find_reheat_rreh_anyorder(nuStar,calFplusNuEnd,epsOneStar,Vstar)
    implicit none
    real(kp) :: find_reheat_rreh_anyorder
    real(kp), intent(in) :: nuStar, calFplusNuEnd, Vstar
    real(kp), intent(in) :: epsOneStar

    find_reheat_rreh_anyorder = nuStar - calFplusNuEnd &
         + 0.5_kp*log((9._kp-3._kp*epsOneStar)/Vstar)

  end function find_reheat_rreh_anyorder



  function get_calfconst_rreh(lnR,epsEnd,potEnd)
    implicit none
    real(kp) :: get_calfconst_rreh
    real(kp), intent(in) :: lnR,epsEnd,potEnd

    get_calfconst_rreh = -Nzero - 0.5_kp*log(potEnd/(3._kp-epsEnd)) + lnR
            
  end function get_calfconst_rreh
  
  
  


  function ln_rho_endinf_leadorder(Pstar,epsOneStar,epsOneEnd,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_endinf_leadorder
    real(kp), intent(in) :: Pstar, epsOneStar, epsOneEnd
    real(kp), intent(in) :: VendOverVstar

    real(kp) :: lnH2OverEps

    lnH2OverEps = ln_hsquare_by_epsone(Pstar)
    
    ln_rho_endinf_leadorder = lnH2OverEps + log(9._kp*epsOneStar &
         /(3._kp - epsOneEnd)) + log(VendOverVstar)

  end function ln_rho_endinf_leadorder


 
  function ln_rho_endinf_anyorder(Pstar,epsStarVec,epsOneEnd,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_endinf_anyorder
    real(kp), intent(in) :: Pstar, epsOneEnd, VendOverVstar
    real(kp), dimension(neps), intent(in) :: epsStarVec

    real(kp) :: lnH2OverEps

    lnH2OverEps = ln_hsquare_by_epsone(Pstar,epsStarVec)

    ln_rho_endinf_anyorder = lnH2OverEps &
         + log(3._kp*epsStarVec(1)*(3._kp - epsStarVec(1))/(3._kp - epsOneEnd)) &
         + log(VendOverVstar)

  end function ln_rho_endinf_anyorder



  function ln_rho_reheat_leadorder(w,Pstar,epsOneStar,epsOneEnd,deltaNstar,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_reheat_leadorder
    real(kp), intent(in) :: w, Pstar, epsOneStar, epsOneEnd, deltaNstar
    real(kp), intent(in) :: VendOverVstar

    real(kp) :: lnH2OverEps, oneMinusThreeWlnEreh

     if (w.eq.1._kp/3._kp) then
        write(*,*)'ln_rho_reheat: w = 1/3!'
        stop
     end if

     lnH2OverEps = ln_hsquare_by_epsone(Pstar)


     oneMinusThreeWlnEreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsOneStar) &
         + log(9._kp/(epsOneStar*(3._kp - epsOneEnd))) &
         + log(VendOverVstar)

     ln_rho_reheat_leadorder = oneMinusThreeWlnEreh*4._kp/(1._kp - 3._kp*w)

   end function ln_rho_reheat_leadorder


   
   function ln_rho_reheat_anyorder(w,Pstar,epsStarVec,epsOneEnd &
        ,deltaNstar,VendOverVstar)
    implicit none
    real(kp) :: ln_rho_reheat_anyorder
    real(kp), intent(in) :: w, Pstar, epsOneEnd, deltaNstar
    real(kp), intent(in) :: VendOverVstar
    real(kp), dimension(neps), intent(in) :: epsStarVec

    real(kp) :: lnH2OverEps, oneMinusThreeWlnEreh

     if (w.eq.1._kp/3._kp) then
        write(*,*)'ln_rho_reheat: w = 1/3!'
        stop
     end if

     lnH2OverEps = ln_hsquare_by_epsone(Pstar,epsStarVec)

     oneMinusThreeWlnEreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsStarVec(1)) &
         + log(3._kp*(3._kp - epsStarVec(1))/(epsStarVec(1)*(3._kp - epsOneEnd))) &
         + log(VendOverVstar)

     ln_rho_reheat_anyorder = oneMinusThreeWlnEreh*4._kp/(1._kp - 3._kp*w)

   end function ln_rho_reheat_anyorder


   function get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd)
     implicit none
     real(kp) :: get_lnrrad_rhow
     real(kp), intent(in) :: lnRhoReh,w,lnRhoEnd

     get_lnrrad_rhow = (1._kp-3._kp*w)/(12._kp+12._kp*w)*(lnRhoReh-lnRhoEnd)

   end function get_lnrrad_rhow


   function get_lnrreh_rhow(lnRhoReh,w,lnRhoEnd)
     implicit none
     real(kp) :: get_lnrreh_rhow
     real(kp), intent(in) :: lnRhoReh,w,lnRhoEnd
     real(kp) :: lnRrad

     lnRrad = get_lnrrad_rhow(lnRhoReh,w,lnRhoEnd)

     get_lnrreh_rhow = lnRrad + 0.25_kp*lnRhoEnd
     
   end function get_lnrreh_rhow



   function get_lnrrad_rreh(lnR,lnRhoEnd)
     implicit none
     real(kp) :: get_lnrrad_rreh
     real(kp), intent(in) :: lnR, lnRhoEnd

     get_lnrrad_rreh = lnR - 0.25_kp*lnRhoEnd

   end function get_lnrrad_rreh


   function get_lnrreh_rrad(lnRrad,lnRhoEnd)
     implicit none
     real(kp) :: get_lnrreh_rrad
     real(kp), intent(in) :: lnRrad, lnRhoEnd

     get_lnrreh_rrad = lnRrad + 0.25_kp*lnRhoEnd

   end function get_lnrreh_rrad

  


  function log_energy_reheat_ingev(lnRhoReh)
    implicit none
    real(kp) :: log_energy_reheat_ingev
    real(kp), intent(in) :: lnRhoReh

    log_energy_reheat_ingev=0.25_kp*(lnRhoReh + 4._kp*lnMpinGeV)/log(10._kp)

  end function log_energy_reheat_ingev


function ln_hsquare_by_epsone(Pstar,epsStarVec)
    implicit none
    real(kp) :: ln_hsquare_by_epsone
    real(kp), intent(in) :: Pstar
    real(kp), dimension(neps), intent(in), optional :: epsStarVec

    ln_hsquare_by_epsone = log(Pstar*8*pi*pi)

    if (present(epsStarVec)) then
       ln_hsquare_by_epsone = ln_hsquare_by_epsone &
            + log(inverse_slowroll_corrections(epsStarVec))
    endif
         

  end function ln_hsquare_by_epsone

 

!this is ln(M)
  function ln_potential_normalization_leadorder(Pstar,epsOneStar,Vstar)
    implicit none
    real(kp) :: ln_potential_normalization_leadorder
    real(kp), intent(in) :: Pstar,epsOneStar,Vstar
    
    real(kp) :: lnM4, lnH2overEps

    lnH2OverEps = ln_hsquare_by_epsone(Pstar)

    lnM4 = log(3._kp*epsOneStar/Vstar) + lnH2OverEps

    ln_potential_normalization_leadorder = 0.25_kp*lnM4

  end function ln_potential_normalization_leadorder


  function ln_potential_normalization_anyorder(Pstar,epsStarVec,Vstar)
    implicit none
    real(kp) :: ln_potential_normalization_anyorder
    real(kp), dimension(neps), intent(in) :: epsStarVec
    real(kp), intent(in) :: Pstar,Vstar
    
    real(kp) :: lnM4, lnH2overEps

    lnH2OverEps = ln_hsquare_by_epsone(Pstar,epsStarVec)

    lnM4 = log((3._kp-epsStarVec(1))*epsStarVec(1)/Vstar) + lnH2OverEps

    ln_potential_normalization_anyorder = 0.25_kp*lnM4

  end function ln_potential_normalization_anyorder



  function quadrupole_to_primscalar(QoverT)
    implicit none
    real(kp) :: quadrupole_to_primscalar
    real(kp), intent(in) :: QoverT
    real(kp) :: H2OverEpsOneOverPi2,Pstar

    H2OverEpsOneOverPi2 = 480._kp*(QoverT)**2

    Pstar = H2OverEpsoneOverPi2/8._kp

    quadrupole_to_primscalar = Pstar

  end function quadrupole_to_primscalar
  

  function primscalar_to_quadrupole(Pstar)
    implicit none
    real(kp) :: primscalar_to_quadrupole
    real(kp), intent(in) :: Pstar
    real(kp) :: QoverT

    QoverT = sqrt(60._kp*Pstar)

    primscalar_to_quadrupole = QoverT

  end function primscalar_to_quadrupole





end module srreheat


!generic reheating functions assuming slow-roll evolution and a mean
!values for wreh or input lnRrad
module srreheat
  use infprec, only : kp,pi
  use srflow, only : slowroll_violated
  use srflow, only : slowroll_corrections, ln_slowroll_corrections
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
  
  interface potential_normalization
     module procedure potential_normalization_leadorder
     module procedure potential_normalization_anyorder
  end interface potential_normalization

  interface primscalar
     module procedure primscalar_leadorder
     module procedure primscalar_anyorder
  end interface primscalar



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
  public potential_normalization, primscalar

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
         * log( 9._kp/( epsStar**(0.5_kp+1.5_kp*w)*Vstar) )

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
         /( epsStarVec(1)**(0.5_kp+1.5_kp*w)*Vstar) ) &
         + 0.25_kp*log(slowroll_corrections(epsStarVec))

  end function find_reheat_rhow_anyorder



  function get_calfconst_rhow(lnRhoReh,Pstar,w,epsEnd,potEnd,lnOmega4End)
    implicit none
    real(kp) :: get_calfconst_rhow
    real(kp), intent(in) :: lnRhoReh,Pstar,w,epsEnd,potEnd
    real(kp), intent(in), optional :: lnOmega4End

    real(kp) :: cmbMeasuredHalf

    cmbMeasuredHalf = 0.5_kp*log(Pstar*8._kp*pi**2)

    get_calfconst_rhow = -Nzero + (1._kp+3._kp*w)/(3._kp+3._kp*w)*cmbMeasuredHalf &
         - 1._kp/(3._kp+3._kp*w)*log(potEnd/(3._kp-epsEnd)) &
         + (1._kp-3._kp*w)/(12._kp+12._kp*w)*lnRhoReh

    if (present(lnOmega4End)) then
       get_calfconst_rhow = get_calfconst_rhow - (1._kp-3._kp*w)/(12._kp+12._kp*w)*lnOmega4End
    endif
    
        
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
         + 0.25_kp*log(slowroll_corrections(epsStarVec))

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



  function get_calfconst_rreh(lnR,epsEnd,potEnd,lnOmega4End)
    implicit none
    real(kp) :: get_calfconst_rreh
    real(kp), intent(in) :: lnR,epsEnd,potEnd
!for scalar-tensor theories, conformal factor ln(Omega^4)
    real(kp), intent(in), optional :: lnOmega4End
    
    get_calfconst_rreh = -Nzero - 0.5_kp*log(potEnd/(3._kp-epsEnd)) + lnR

    if (present(lnOmega4End)) then
       get_calfconst_rreh = get_calfconst_rreh - 0.25_kp*lnOmega4End
    endif
            
  end function get_calfconst_rreh
  
  
  

!if lnOmega4End is provided, this function assumes that slow-roll
!calculations are done in the Einstein Frame and return rhoEndInfBar,
!in the Jordan Frame!
  function ln_rho_endinf_leadorder(Pstar,epsOneStar,epsOneEnd,VendOverVstar,lnOmega4End)
    implicit none
    real(kp) :: ln_rho_endinf_leadorder
    real(kp), intent(in) :: Pstar, epsOneStar, epsOneEnd
    real(kp), intent(in) :: VendOverVstar
    real(kp), intent(in), optional :: lnOmega4End

    real(kp) :: lnH2OverEps

    lnH2OverEps = log(primscalar_to_hsquare(Pstar))
    
    ln_rho_endinf_leadorder = lnH2OverEps + log(9._kp*epsOneStar &
         /(3._kp - epsOneEnd)) + log(VendOverVstar)

    if (present(lnOmega4End)) then
       ln_rho_endinf_leadorder = ln_rho_endinf_leadorder + lnOmega4End
    endif
    
  end function ln_rho_endinf_leadorder


!if lnOmega4End is provided, this function assumes that slow-roll
!calculations are done in the Einstein Frame and return rhoEndInfBar,
!in the Jordan Frame!
  function ln_rho_endinf_anyorder(Pstar,epsStarVec,epsOneEnd,VendOverVstar,lnOmega4End)
    implicit none
    real(kp) :: ln_rho_endinf_anyorder
    real(kp), intent(in) :: Pstar, epsOneEnd, VendOverVstar
    real(kp), dimension(neps), intent(in) :: epsStarVec
    real(kp), intent(in), optional :: lnOmega4End
    
    real(kp) :: lnH2OverEps

    lnH2OverEps = log(primscalar_to_hsquare(Pstar,epsStarVec))

    ln_rho_endinf_anyorder = lnH2OverEps &
         + log(3._kp*epsStarVec(1)*(3._kp - epsStarVec(1))/(3._kp - epsOneEnd)) &
         + log(VendOverVstar)

    if (present(lnOmega4End)) then
       ln_rho_endinf_anyorder = ln_rho_endinf_anyorder + lnOmega4End
    endif
    
  end function ln_rho_endinf_anyorder



!if lnOmega4End is provided, this function assumes that slow-roll
!calculations are done in the Einstein Frame and return rhoEndInfBar,
!in the Jordan Frame!
  function ln_rho_reheat_leadorder(w,Pstar,epsOneStar,epsOneEnd,deltaNstar,VendOverVstar &
       ,lnOmega4End)
    implicit none
    real(kp) :: ln_rho_reheat_leadorder
    real(kp), intent(in) :: w, Pstar, epsOneStar, epsOneEnd, deltaNstar
    real(kp), intent(in) :: VendOverVstar
    real(kp), intent(in), optional :: lnOmega4End
    
    real(kp) :: lnH2OverEps, oneMinusThreeWlnEreh

     if (w.eq.1._kp/3._kp) then
        write(*,*)'ln_rho_reheat: w = 1/3!'
        stop
     end if

     lnH2OverEps = log(primscalar_to_hsquare(Pstar))


     oneMinusThreeWlnEreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsOneStar) &
         + log(9._kp/(epsOneStar*(3._kp - epsOneEnd))) &
         + log(VendOverVstar)

     ln_rho_reheat_leadorder = oneMinusThreeWlnEreh*4._kp/(1._kp - 3._kp*w)

     if (present(lnOmega4End)) then
        ln_rho_reheat_leadorder = ln_rho_reheat_leadorder + lnOmega4End
     endif
     
   end function ln_rho_reheat_leadorder


!if lnOmega4End is provided, this function assumes that slow-roll
!calculations are done in the Einstein Frame and return rhoEndInfBar,
!in the Jordan Frame!
   function ln_rho_reheat_anyorder(w,Pstar,epsStarVec,epsOneEnd,deltaNstar,VendOverVstar &
        , lnOmega4End)
    implicit none
    real(kp) :: ln_rho_reheat_anyorder
    real(kp), intent(in) :: w, Pstar, epsOneEnd, deltaNstar
    real(kp), intent(in) :: VendOverVstar
    real(kp), dimension(neps), intent(in) :: epsStarVec
    real(kp), intent(in), optional :: lnOmega4End
    
    real(kp) :: lnH2OverEps, oneMinusThreeWlnEreh

     if (w.eq.1._kp/3._kp) then
        write(*,*)'ln_rho_reheat: w = 1/3!'
        stop
     end if

     lnH2OverEps = log(primscalar_to_hsquare(Pstar,epsStarVec))

     oneMinusThreeWlnEreh &
         = (3._kp + 3._kp*w)*(Nzero + deltaNstar) - 0.5_kp*(1._kp+3._kp*w)*lnH2OverEps &
         + 0.5_kp*(1._kp - 3._kp*w)*log(epsStarVec(1)) &
         + log(3._kp*(3._kp - epsStarVec(1))/(epsStarVec(1)*(3._kp - epsOneEnd))) &
         + log(VendOverVstar)

     ln_rho_reheat_anyorder = oneMinusThreeWlnEreh*4._kp/(1._kp - 3._kp*w)

     if (present(lnOmega4End)) then
        ln_rho_reheat_anyorder = ln_rho_reheat_anyorder + lnOmega4End
     endif
     
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


!returns H^2/epsilon_1 fropm Pstar
  function primscalar_to_hsquare(Pstar,epsStarVec)
    implicit none
    real(kp) :: primscalar_to_hsquare
    real(kp), intent(in) :: Pstar
    real(kp), dimension(neps), intent(in), optional :: epsStarVec

    real(kp) :: H2overEps

    H2overEps = Pstar*8*pi*pi

    if (present(epsStarVec)) then
        H2overEps = H2overEps / slowroll_corrections(epsStarVec)
    endif

    primscalar_to_hsquare = H2overEps

  end function primscalar_to_hsquare

!for test, that is not strictly ln(primscalar_to_hsquare) at second
!order, that's why this function exists
  function primscalar_to_lnhsquare(Pstar,epsStarVec)
    implicit none
    real(kp) :: primscalar_to_lnhsquare
    real(kp), intent(in) :: Pstar
    real(kp), dimension(neps), intent(in), optional :: epsStarVec

    real(kp) :: lnH2overEps

    stop 'primscalar_to_lnhsquare: do not use this function!'

    lnH2overEps = log(Pstar*8*pi*pi)

    if (present(epsStarVec)) then
       lnH2overEps = lnH2overEps - ln_slowroll_corrections(epsStarVec)
    endif

    primscalar_to_lnhsquare = lnH2overEps

  end function primscalar_to_lnhsquare



!returns Pstar from H^2/epsilon_1
  function hsquare_to_primscalar(H2overEps,epsStarVec)
    implicit none
    real(kp) :: hsquare_to_primscalar
    real(kp), intent(in) :: H2overEps
    real(kp), dimension(neps), intent(in), optional :: epsStarVec

    real(kp) :: Pstar

    Pstar = H2overEps/(8*pi*pi)

    if (present(epsStarVec)) then
       Pstar = Pstar * slowroll_corrections(epsStarVec)
    endif

    hsquare_to_primscalar = Pstar

  end function hsquare_to_primscalar

!for test, this is slightly different than ln(hsquare_to_primscalar)
  function hsquare_to_lnprimscalar(H2overEps,epsStarVec)
    implicit none
    real(kp) :: hsquare_to_lnprimscalar
    real(kp), intent(in) :: H2overEps
    real(kp), dimension(neps), intent(in), optional :: epsStarVec

    real(kp) :: lnPstar

    stop 'primscalar_to_lnhsquare: do not use this function!'

    lnPstar = log(H2overEps/(8*pi*pi))

    if (present(epsStarVec)) then
       lnPstar = lnPstar + ln_slowroll_corrections(epsStarVec)
    endif

    hsquare_to_lnprimscalar = lnPstar

  end function hsquare_to_lnprimscalar


  
!this is M from Pstar and eps
  function potential_normalization_leadorder(Pstar,epsOneStar,Vstar)
    implicit none
    real(kp) :: potential_normalization_leadorder
    real(kp), intent(in) :: Pstar,epsOneStar,Vstar
    
    real(kp) :: M4, H2overEps

    H2OverEps = primscalar_to_hsquare(Pstar)

    M4 = 3._kp*epsOneStar/Vstar * H2overEps

    potential_normalization_leadorder = M4**0.25_kp

  end function potential_normalization_leadorder


  function potential_normalization_anyorder(Pstar,epsStarVec,Vstar)
    implicit none
    real(kp) :: potential_normalization_anyorder
    real(kp), dimension(neps), intent(in) :: epsStarVec
    real(kp), intent(in) :: Pstar,Vstar
    
    real(kp) :: M4, H2overEps

    H2overEps = primscalar_to_hsquare(Pstar,epsStarVec)

    M4 = (3._kp-epsStarVec(1))*epsStarVec(1)/Vstar * H2overEps

    potential_normalization_anyorder = M4**0.25_kp

  end function potential_normalization_anyorder



!this is Pstar from M and eps
  function primscalar_leadorder(M,epsOneStar,Vstar)
    implicit none
    real(kp) :: primscalar_leadorder
    real(kp), intent(in) :: M,epsOneStar,Vstar
    
    real(kp) :: H2overEps

    H2overEps = M**4*Vstar/3._kp/epsOneStar

    primscalar_leadorder = hsquare_to_primscalar(H2overEps)
    
  end function primscalar_leadorder
  

  function primscalar_anyorder(M,epsStarVec,Vstar)
    implicit none
    real(kp) :: primscalar_anyorder
    real(kp), dimension(neps), intent(in) :: epsStarVec
    real(kp), intent(in) :: M,Vstar
    
    real(kp) :: H2overEps

    H2overEps = M**4*Vstar/(3._kp-epsStarVec(1))/epsStarVec(1)

    primscalar_anyorder = hsquare_to_primscalar(H2overEps,epsStarVec)
    
  end function primscalar_anyorder


!Pstar from Q/T effective
  function quadrupole_to_primscalar(QoverT)
    implicit none
    real(kp) :: quadrupole_to_primscalar
    real(kp), intent(in) :: QoverT
    real(kp) :: H2OverEpsOneOverPi2,Pstar

    H2OverEpsOneOverPi2 = 480._kp*(QoverT)**2

    Pstar = H2OverEpsoneOverPi2/8._kp

    quadrupole_to_primscalar = Pstar

  end function quadrupole_to_primscalar
  

!Q/T effective from Pstar
  function primscalar_to_quadrupole(Pstar)
    implicit none
    real(kp) :: primscalar_to_quadrupole
    real(kp), intent(in) :: Pstar
    real(kp) :: QoverT

    QoverT = sqrt(Pstar/60._kp)

    primscalar_to_quadrupole = QoverT

  end function primscalar_to_quadrupole





end module srreheat


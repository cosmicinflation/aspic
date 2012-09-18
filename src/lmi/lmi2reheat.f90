!Logamediate inflation 2 reheating functions in the slow-roll approximations

module lmi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  use lmi2sr, only : lmi2_norm_potential
  use lmi2sr, only : lmi2_xini_min, lmi2_efold_primitive
  use cosmopar, only : QrmsOverT
  implicit none

  private

  public lmi2_x_star, lmi2_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi2_x_star(gamma,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi2_x_star
    real(kp), intent(in) :: gamma,beta,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: lmi2Data

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = lmi2_epsilon_one(xEnd,gamma,beta)
    potEnd = lmi2_norm_potential(xEnd,gamma,beta)

    primEnd = lmi2_efold_primitive(xEnd,gamma,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi2Data%real1 = gamma
    lmi2Data%real2 = beta
    lmi2Data%real3 = w
    lmi2Data%real4 = calF + primEnd

    mini = lmi2_xini_min(gamma,beta)
    maxi = xEnd*(1._kp-epsilon(1._kp))


    x = zbrent(find_lmi2_x_star,mini,maxi,tolzbrent,lmi2Data)
    lmi2_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi2_efold_primitive(x,gamma,beta) - primEnd)
    endif

  end function lmi2_x_star

  function find_lmi2_x_star(x,lmi2Data)   
    implicit none
    real(kp) :: find_lmi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data

    real(kp) :: primStar,gamma,beta,w,CalFplusprimEnd,potStar,epsOneStar

    gamma=lmi2Data%real1
    beta=lmi2Data%real2
    w = lmi2Data%real3
    CalFplusprimEnd = lmi2Data%real4

    primStar = lmi2_efold_primitive(x,gamma,beta)
    epsOneStar = lmi2_epsilon_one(x,gamma,beta)
    potStar = lmi2_norm_potential(x,gamma,beta)

    find_lmi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi2_x_star



  function lmi2_lnrhoend(gamma,beta,xEnd,Pstar) 
    implicit none
    real(kp) :: lmi2_lnrhoend
    real(kp), intent(in) :: gamma,beta,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd

    potEnd  = lmi2_norm_potential(xEnd,gamma,beta)
    epsOneEnd = lmi2_epsilon_one(xEnd,gamma,beta)


!   Trick to return x such that rho_reh=rho_end

    x = lmi2_x_star(gamma,beta,xEnd,wrad,junk,Pstar)  

    potStar = lmi2_norm_potential(x,gamma,beta)
    epsOneStar = lmi2_epsilon_one(x,gamma,beta)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi2_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi2_lnrhoend = lnRhoEnd

  end function lmi2_lnrhoend

  
end module lmi2reheat

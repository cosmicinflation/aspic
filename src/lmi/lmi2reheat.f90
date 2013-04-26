!Logamediate inflation 2 reheating functions in the slow-roll approximations

module lmi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use lmi2sr, only : lmi2_epsilon_one, lmi2_epsilon_two, lmi2_epsilon_three
  use lmi2sr, only : lmi2_norm_potential
  use lmi2sr, only : lmi2_xini_min, lmi2_efold_primitive

  implicit none

  private

  public lmi2_x_star, lmi2_lnrhoreh_max
  public lmi2_x_rrad, lmi2_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi2_x_star(gam,beta,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi2_x_star
    real(kp), intent(in) :: gam,beta,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: lmi2Data

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = lmi2_epsilon_one(xEnd,gam,beta)
    potEnd = lmi2_norm_potential(xEnd,gam,beta)

    primEnd = lmi2_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi2Data%real1 = gam
    lmi2Data%real2 = beta
    lmi2Data%real3 = w
    lmi2Data%real4 = calF + primEnd

    mini = lmi2_xini_min(gam,beta)
    maxi = xEnd*(1._kp-epsilon(1._kp))


    x = zbrent(find_lmi2_x_star,mini,maxi,tolzbrent,lmi2Data)
    lmi2_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi2_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi2_x_star

  function find_lmi2_x_star(x,lmi2Data)   
    implicit none
    real(kp) :: find_lmi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data

    real(kp) :: primStar,gam,beta,w,CalFplusprimEnd,potStar,epsOneStar

    gam=lmi2Data%real1
    beta=lmi2Data%real2
    w = lmi2Data%real3
    CalFplusprimEnd = lmi2Data%real4

    primStar = lmi2_efold_primitive(x,gam,beta)
    epsOneStar = lmi2_epsilon_one(x,gam,beta)
    potStar = lmi2_norm_potential(x,gam,beta)

    find_lmi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi2_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function lmi2_x_rrad(gam,beta,xEnd,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi2_x_rrad
    real(kp), intent(in) :: gam,beta,xEnd,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: lmi2Data

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = lmi2_epsilon_one(xEnd,gam,beta)
    potEnd = lmi2_norm_potential(xEnd,gam,beta)

    primEnd = lmi2_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    lmi2Data%real1 = gam
    lmi2Data%real2 = beta
    lmi2Data%real3 = calF + primEnd

    mini = lmi2_xini_min(gam,beta)
    maxi = xEnd*(1._kp-epsilon(1._kp))


    x = zbrent(find_lmi2_x_rrad,mini,maxi,tolzbrent,lmi2Data)
    lmi2_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi2_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi2_x_rrad

  function find_lmi2_x_rrad(x,lmi2Data)   
    implicit none
    real(kp) :: find_lmi2_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data

    real(kp) :: primStar,gam,beta,CalFplusprimEnd,potStar,epsOneStar

    gam=lmi2Data%real1
    beta=lmi2Data%real2
    CalFplusprimEnd = lmi2Data%real3

    primStar = lmi2_efold_primitive(x,gam,beta)
    epsOneStar = lmi2_epsilon_one(x,gam,beta)
    potStar = lmi2_norm_potential(x,gam,beta)

    find_lmi2_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_lmi2_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function lmi2_x_rreh(gam,beta,xEnd,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: lmi2_x_rreh
    real(kp), intent(in) :: gam,beta,xEnd,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: lmi2Data

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = lmi2_epsilon_one(xEnd,gam,beta)
    potEnd = lmi2_norm_potential(xEnd,gam,beta)

    primEnd = lmi2_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    lmi2Data%real1 = gam
    lmi2Data%real2 = beta
    lmi2Data%real3 = calF + primEnd

    mini = lmi2_xini_min(gam,beta)
    maxi = xEnd*(1._kp-epsilon(1._kp))


    x = zbrent(find_lmi2_x_rreh,mini,maxi,tolzbrent,lmi2Data)
    lmi2_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi2_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi2_x_rreh

  function find_lmi2_x_rreh(x,lmi2Data)   
    implicit none
    real(kp) :: find_lmi2_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi2Data

    real(kp) :: primStar,gam,beta,CalFplusprimEnd,potStar

    gam=lmi2Data%real1
    beta=lmi2Data%real2
    CalFplusprimEnd = lmi2Data%real3

    primStar = lmi2_efold_primitive(x,gam,beta)
    potStar = lmi2_norm_potential(x,gam,beta)

    find_lmi2_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_lmi2_x_rreh


  function lmi2_lnrhoreh_max(gam,beta,xEnd,Pstar) 
    implicit none
    real(kp) :: lmi2_lnrhoreh_max
    real(kp), intent(in) :: gam,beta,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd

    potEnd  = lmi2_norm_potential(xEnd,gam,beta)
    epsOneEnd = lmi2_epsilon_one(xEnd,gam,beta)


!   Trick to return x such that rho_reh=rho_end

    x = lmi2_x_star(gam,beta,xEnd,wrad,junk,Pstar)  

    potStar = lmi2_norm_potential(x,gam,beta)
    epsOneStar = lmi2_epsilon_one(x,gam,beta)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi2_lnrhoreh_max = lnRhoEnd

  end function lmi2_lnrhoreh_max

  
end module lmi2reheat

!radiatively corrected quartic inflation reheating functions in the slow-roll approximations

module rchireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rchisr, only : rchi_epsilon_one, rchi_epsilon_two, rchi_epsilon_three
  use rchisr, only : rchi_norm_potential, rchi_potmax_exists, rchi_x_potmax
  use rchisr, only : rchi_x_endinf, rchi_efold_primitive, RchiDeltaXendMax
  implicit none

  private

  public rchi_x_star, rchi_lnrhoreh_max
  public rchi_x_rrad, rchi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rchi_x_star(AI,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rchi_x_star
    real(kp), intent(in) :: AI,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rchiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rchi_epsilon_one(xEnd,AI)
    potEnd = rchi_norm_potential(xEnd,AI)
    primEnd = rchi_efold_primitive(xEnd,AI)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rchiData%real1 = AI    
    rchiData%real2 = w
    rchiData%real3 = calF + primEnd

    mini=xend*(1._kp+epsilon(1._kp))
    if (rchi_potmax_exists(AI)) then !maximum of the potential
      maxi=min(rchi_x_potmax(AI),xend + RchiDeltaXendMax)
    else
      maxi = xend + RchiDeltaXendMax
    endif

    x = zbrent(find_rchi_x_star,mini,maxi,tolzbrent,rchiData)
    rchi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rchi_efold_primitive(x,AI) - primEnd)
    endif

  end function rchi_x_star

  function find_rchi_x_star(x,rchiData)   
    implicit none
    real(kp) :: find_rchi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rchiData

    real(kp) :: primStar,AI,w,CalFplusprimEnd,potStar,epsOneStar

    AI=rchiData%real1
    w = rchiData%real2
    CalFplusprimEnd = rchiData%real3

    primStar = rchi_efold_primitive(x,AI)
    epsOneStar = rchi_epsilon_one(x,AI)
    potStar = rchi_norm_potential(x,AI)

    find_rchi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rchi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rchi_x_rrad(AI,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rchi_x_rrad
    real(kp), intent(in) :: AI,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rchiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rchi_epsilon_one(xEnd,AI)
    potEnd = rchi_norm_potential(xEnd,AI)
    primEnd = rchi_efold_primitive(xEnd,AI)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rchiData%real1 = AI    
    rchiData%real2 = calF + primEnd

    mini=xend*(1._kp+epsilon(1._kp))
    if (rchi_potmax_exists(AI)) then !maximum of the potential
      maxi=min(rchi_x_potmax(AI),xend + RchiDeltaXendMax)
    else
      maxi = xend + RchiDeltaXendMax
    endif

    x = zbrent(find_rchi_x_rrad,mini,maxi,tolzbrent,rchiData)
    rchi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (rchi_efold_primitive(x,AI) - primEnd)
    endif

  end function rchi_x_rrad

  function find_rchi_x_rrad(x,rchiData)   
    implicit none
    real(kp) :: find_rchi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rchiData

    real(kp) :: primStar,AI,CalFplusprimEnd,potStar,epsOneStar

    AI=rchiData%real1
    CalFplusprimEnd = rchiData%real2

    primStar = rchi_efold_primitive(x,AI)
    epsOneStar = rchi_epsilon_one(x,AI)
    potStar = rchi_norm_potential(x,AI)

    find_rchi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rchi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rchi_x_rreh(AI,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rchi_x_rreh
    real(kp), intent(in) :: AI,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rchiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rchi_epsilon_one(xEnd,AI)
    potEnd = rchi_norm_potential(xEnd,AI)
    primEnd = rchi_efold_primitive(xEnd,AI)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rchiData%real1 = AI    
    rchiData%real2 = calF + primEnd

    mini=xend*(1._kp+epsilon(1._kp))
    if (rchi_potmax_exists(AI)) then !maximum of the potential
      maxi=min(rchi_x_potmax(AI),xend + RchiDeltaXendMax)
    else
      maxi = xend + RchiDeltaXendMax
    endif

    x = zbrent(find_rchi_x_rreh,mini,maxi,tolzbrent,rchiData)
    rchi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (rchi_efold_primitive(x,AI) - primEnd)
    endif

  end function rchi_x_rreh

  function find_rchi_x_rreh(x,rchiData)   
    implicit none
    real(kp) :: find_rchi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rchiData

    real(kp) :: primStar,AI,CalFplusprimEnd,potStar

    AI=rchiData%real1
    CalFplusprimEnd = rchiData%real2

    primStar = rchi_efold_primitive(x,AI)
    potStar = rchi_norm_potential(x,AI)

    find_rchi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rchi_x_rreh



  function rchi_lnrhoreh_max(AI,xend,Pstar) 
    implicit none
    real(kp) :: rchi_lnrhoreh_max
    real(kp), intent(in) :: AI,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = rchi_norm_potential(xEnd,AI)
    epsOneEnd = rchi_epsilon_one(xEnd,AI)

!   Trick to return x such that rho_reh=rho_end

    x = rchi_x_star(AI,xend,wrad,junk,Pstar)    
    potStar = rchi_norm_potential(x,AI)
    epsOneStar = rchi_epsilon_one(x,AI)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rchi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rchi_lnrhoreh_max = lnRhoEnd

  end function rchi_lnrhoreh_max

  
end module rchireheat

!double well inflation reheating functions in the slow-roll approximations

module wrireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use wrisr, only : wri_epsilon_one, wri_epsilon_two, wri_epsilon_three
  use wrisr, only : wri_norm_potential, wri_x_trajectory
  use wrisr, only : wri_x_endinf, wri_efold_primitive
  implicit none

  private

  public wri_x_star, wri_lnrhoreh_max 
  public wri_x_rrad, wri_x_rreh

contains

!returns x =phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function wri_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: wri_x_star
    real(kp), intent(in) :: phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: wriData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = wri_x_endinf(phi0)

    epsOneEnd = wri_epsilon_one(xEnd,phi0)
    potEnd = wri_norm_potential(xEnd,phi0)

    primEnd = wri_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    wriData%real1 = phi0
    wriData%real2 = w
    wriData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi =  wri_x_trajectory(-200._kp,xend,phi0) !Position 200 efolds before the end of inflation 

    x = zbrent(find_wri_x_star,mini,maxi,tolzbrent,wriData)
    wri_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (wri_efold_primitive(x,phi0) - primEnd)
    endif

  end function wri_x_star


  function find_wri_x_star(x,wriData)   
    implicit none
    real(kp) :: find_wri_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: wriData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=wriData%real1
    w = wriData%real2
    CalFplusprimEnd = wriData%real3

    primStar = wri_efold_primitive(x,phi0)
    epsOneStar = wri_epsilon_one(x,phi0)
    potStar = wri_norm_potential(x,phi0)

    find_wri_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_wri_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function wri_x_rrad(phi0,lnRRad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: wri_x_rrad
    real(kp), intent(in) :: phi0,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: wriData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = wri_x_endinf(phi0)

    epsOneEnd = wri_epsilon_one(xEnd,phi0)
    potEnd = wri_norm_potential(xEnd,phi0)

    primEnd = wri_efold_primitive(xEnd,phi0)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)


    wriData%real1 = phi0
    wriData%real2 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi =  wri_x_trajectory(-200._kp,xend,phi0) !Position 200 efolds before the end of inflation 

    x = zbrent(find_wri_x_rrad,mini,maxi,tolzbrent,wriData)
    wri_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (wri_efold_primitive(x,phi0) - primEnd)
    endif

  end function wri_x_rrad

  function find_wri_x_rrad(x,wriData)   
    implicit none
    real(kp) :: find_wri_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: wriData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar,epsOneStar

    phi0 = wriData%real1
    CalFplusprimEnd = wriData%real2

    primStar = wri_efold_primitive(x,phi0)
    epsOneStar = wri_epsilon_one(x,phi0)
    potStar = wri_norm_potential(x,phi0)

    find_wri_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_wri_x_rrad


  !returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function wri_x_rreh(phi0,lnRReh,bfoldstar)    
    implicit none
    real(kp) :: wri_x_rreh
    real(kp), intent(in) :: phi0,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: wriData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = wri_x_endinf(phi0)

    epsOneEnd = wri_epsilon_one(xEnd,phi0)
    potEnd = wri_norm_potential(xEnd,phi0)

    primEnd = wri_efold_primitive(xEnd,phi0)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)


    wriData%real1 = phi0
    wriData%real2 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi =  wri_x_trajectory(-200._kp,xend,phi0) !Position 200 efolds before the end of inflation 

    x = zbrent(find_wri_x_rreh,mini,maxi,tolzbrent,wriData)
    wri_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (wri_efold_primitive(x,phi0) - primEnd)
    endif

  end function wri_x_rreh

  function find_wri_x_rreh(x,wriData)   
    implicit none
    real(kp) :: find_wri_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: wriData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar

    phi0 = wriData%real1
    CalFplusprimEnd = wriData%real2

    primStar = wri_efold_primitive(x,phi0)    
    potStar = wri_norm_potential(x,phi0)

    find_wri_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_wri_x_rreh


  function wri_lnrhoreh_max(phi0,Pstar) 
    implicit none
    real(kp) :: wri_lnrhoreh_max
    real(kp), intent(in) :: phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = wri_x_endinf(phi0)
    potEnd  = wri_norm_potential(xEnd,phi0)
    epsOneEnd = wri_epsilon_one(xEnd,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = wri_x_star(phi0,wrad,junk,Pstar)  
 
    potStar = wri_norm_potential(x,phi0)
    epsOneStar = wri_epsilon_one(x,phi0)
   
    if (.not.slowroll_validity(epsOneStar)) stop 'wri_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    wri_lnrhoreh_max = lnRhoEnd

  end function wri_lnrhoreh_max
  
end module wrireheat

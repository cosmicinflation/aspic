!Twisted inflation reheating functions in the slow-roll approximations

module twireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use twisr, only : twi_epsilon_one, twi_epsilon_two, twi_epsilon_three
  use twisr, only : twi_norm_potential
  use twisr, only : twi_efold_primitive
  implicit none

  private

!if bigger, <numaccuracy errors
  real(kp), parameter :: TwiMaxFactor = 100._kp

  real(kp), parameter :: epsVclamp = 2._kp
  
  public twi_x_star, twi_lnrhoreh_max
  public twi_x_rrad, twi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function twi_x_star(phi0,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: twi_x_star
    real(kp), intent(in) :: phi0,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: twiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

!there are possibly slow-roll violations at the end of inflation, do
!not feed the reheating finder with crap
    epsOneEnd = min(twi_epsilon_one(xEnd,phi0),epsVclamp)
    potEnd = twi_norm_potential(xEnd,phi0)

    primEnd = twi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    twiData%real1 = phi0
    twiData%real2 = w
    twiData%real3 = calF + primEnd

    mini = xEnd
    maxi=TwiMaxFactor*phi0 

    x = zbrent(find_twi_x_star,mini,maxi,tolzbrent,twiData)
    twi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (twi_efold_primitive(x,phi0) - primEnd)
    endif

  end function twi_x_star

  function find_twi_x_star(x,twiData)   
    implicit none
    real(kp) :: find_twi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: twiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=twiData%real1
    w = twiData%real2
    CalFplusprimEnd = twiData%real3

    primStar = twi_efold_primitive(x,phi0)
    epsOneStar = twi_epsilon_one(x,phi0)
    potStar = twi_norm_potential(x,phi0)

    find_twi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_twi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function twi_x_rrad(phi0,xEnd,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: twi_x_rrad
    real(kp), intent(in) :: phi0,xEnd,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: twiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

!there are possibly slow-roll violations at the end of inflation, do
!not feed the reheating finder with crap
    epsOneEnd = min(twi_epsilon_one(xEnd,phi0),epsVclamp)
    potEnd = twi_norm_potential(xEnd,phi0)

    primEnd = twi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    twiData%real1 = phi0
    twiData%real2 = calF + primEnd

    mini = xEnd
    maxi=TwiMaxFactor*phi0 

    x = zbrent(find_twi_x_rrad,mini,maxi,tolzbrent,twiData)
    twi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (twi_efold_primitive(x,phi0) - primEnd)
    endif

  end function twi_x_rrad

  function find_twi_x_rrad(x,twiData)   
    implicit none
    real(kp) :: find_twi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: twiData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar,epsOneStar

    phi0=twiData%real1
    CalFplusprimEnd = twiData%real2

    primStar = twi_efold_primitive(x,phi0)
    epsOneStar = twi_epsilon_one(x,phi0)
    potStar = twi_norm_potential(x,phi0)

    find_twi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_twi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function twi_x_rreh(phi0,xEnd,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: twi_x_rreh
    real(kp), intent(in) :: phi0,xEnd,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: twiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

!there are possibly slow-roll violations at the end of inflation, do
!not feed the reheating finder with crap
    epsOneEnd = min(twi_epsilon_one(xEnd,phi0),epsVclamp)
    potEnd = twi_norm_potential(xEnd,phi0)

    primEnd = twi_efold_primitive(xEnd,phi0)
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    twiData%real1 = phi0
    twiData%real2 = calF + primEnd

    mini = xEnd
    maxi=TwiMaxFactor*phi0 

    x = zbrent(find_twi_x_rreh,mini,maxi,tolzbrent,twiData)
    twi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (twi_efold_primitive(x,phi0) - primEnd)
    endif

  end function twi_x_rreh

  function find_twi_x_rreh(x,twiData)   
    implicit none
    real(kp) :: find_twi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: twiData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar

    phi0=twiData%real1
    CalFplusprimEnd = twiData%real2

    primStar = twi_efold_primitive(x,phi0)
    potStar = twi_norm_potential(x,phi0)

    find_twi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_twi_x_rreh



  function twi_lnrhoreh_max(phi0,xEnd,Pstar) 
    implicit none
    real(kp) :: twi_lnrhoreh_max
    real(kp), intent(in) :: phi0,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd


    potEnd  = twi_norm_potential(xEnd,phi0)

!there are possibly slow-roll violations at the end of inflation, do
!not feed the rho reheating equation with crap
    epsOneEnd = min(twi_epsilon_one(xEnd,phi0),epsVclamp)


!   Trick to return x such that rho_reh=rho_end

    x = twi_x_star(phi0,xEnd,wrad,junk,Pstar)  

 
    potStar = twi_norm_potential(x,phi0)
    epsOneStar = twi_epsilon_one(x,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'twi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    twi_lnrhoreh_max = lnRhoEnd

  end function twi_lnrhoreh_max

  
end module twireheat

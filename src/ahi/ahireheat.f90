!Axion Hilltop Inflation reheating functions

module ahireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use ahisr, only : ahi_epsilon_one, ahi_epsilon_two, ahi_epsilon_three
  use ahisr, only : ahi_norm_potential
  use ahisr, only : ahi_x_endinf, ahi_efold_primitive
  implicit none

  private

  public ahi_x_star, ahi_lnrhoreh_max 
  public ahi_x_rrad, ahi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ahi_x_star(f,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ahi_x_star
    real(kp), intent(in) :: f,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ahiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ahi_epsilon_one(xEnd,f)
    potEnd = ahi_norm_potential(xEnd,f)
    primEnd = ahi_efold_primitive(xEnd,f)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ahiData%real1 = f    
    ahiData%real2 = w
    ahiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+tolkp)
    maxi = pi*(1._kp-tolkp)

    x = zbrent(find_ahi_x_star,mini,maxi,tolzbrent,ahiData)
    ahi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ahi_efold_primitive(x,f) - primEnd)
    endif

  end function ahi_x_star

  function find_ahi_x_star(x,ahiData)   
    implicit none
    real(kp) :: find_ahi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ahiData

    real(kp) :: primStar,f,w,CalFplusprimEnd,potStar,epsOneStar

    f=ahiData%real1
    w = ahiData%real2
    CalFplusprimEnd = ahiData%real3

    primStar = ahi_efold_primitive(x,f)
    epsOneStar = ahi_epsilon_one(x,f)
    potStar = ahi_norm_potential(x,f)

    find_ahi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ahi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ahi_x_rrad(f,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ahi_x_rrad
    real(kp), intent(in) :: f,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ahiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ahi_epsilon_one(xEnd,f)
    potEnd = ahi_norm_potential(xEnd,f)
    primEnd = ahi_efold_primitive(xEnd,f)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    ahiData%real1 = f
    ahiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+tolkp)
    maxi = pi*(1._kp-tolkp)

    x = zbrent(find_ahi_x_rrad,mini,maxi,tolzbrent,ahiData)
    ahi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ahi_efold_primitive(x,f) - primEnd)
    endif

  end function ahi_x_rrad

  function find_ahi_x_rrad(x,ahiData)   
    implicit none
    real(kp) :: find_ahi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ahiData

    real(kp) :: primStar,f,CalFplusprimEnd,potStar,epsOneStar

    f=ahiData%real1
    CalFplusprimEnd = ahiData%real2

    primStar = ahi_efold_primitive(x,f)
    epsOneStar = ahi_epsilon_one(x,f)
    potStar = ahi_norm_potential(x,f)

    find_ahi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_ahi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ahi_x_rreh(f,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ahi_x_rreh
    real(kp), intent(in) :: f,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ahiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ahi_epsilon_one(xEnd,f)
    potEnd = ahi_norm_potential(xEnd,f)
    primEnd = ahi_efold_primitive(xEnd,f)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    ahiData%real1 = f
    ahiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+tolkp)
    maxi = pi*(1._kp-tolkp)

    x = zbrent(find_ahi_x_rreh,mini,maxi,tolzbrent,ahiData)
    ahi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ahi_efold_primitive(x,f) - primEnd)
    endif

  end function ahi_x_rreh

  function find_ahi_x_rreh(x,ahiData)   
    implicit none
    real(kp) :: find_ahi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ahiData

    real(kp) :: primStar,f,CalFplusprimEnd,potStar

    f=ahiData%real1
    CalFplusprimEnd = ahiData%real2

    primStar = ahi_efold_primitive(x,f)
    potStar = ahi_norm_potential(x,f)

    find_ahi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_ahi_x_rreh




  function ahi_lnrhoreh_max(f,xend,Pstar) 
    implicit none
    real(kp) :: ahi_lnrhoreh_max
    real(kp), intent(in) :: f,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ahi_norm_potential(xEnd,f)
    epsOneEnd = ahi_epsilon_one(xEnd,f)

    x = ahi_x_star(f,xend,wrad,junk,Pstar)    
    potStar = ahi_norm_potential(x,f)
    epsOneStar = ahi_epsilon_one(x,f)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ahi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ahi_lnrhoreh_max = lnRhoEnd

  end function ahi_lnrhoreh_max

  
end module ahireheat

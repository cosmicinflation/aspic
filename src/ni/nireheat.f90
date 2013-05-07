!Natural Inflation reheating functions in the slow-roll

module nireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use nisr, only : ni_epsilon_one, ni_epsilon_two
  use nisr, only : ni_norm_potential
  use nisr, only : ni_x_endinf, ni_efold_primitive
  implicit none

  private

  public ni_x_star, ni_lnrhoreh_max
  public ni_x_rrad, ni_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function ni_x_star(f,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ni_x_star
    real(kp), intent(in) :: f,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: niData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ni_x_endinf(f)
    epsOneEnd = ni_epsilon_one(xEnd,f)
    potEnd = ni_norm_potential(xEnd,f)
    primEnd = ni_efold_primitive(xEnd,f)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    niData%real1 = f
    niData%real2 = w
    niData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd

    x = zbrent(find_ni_x_star,mini,maxi,tolzbrent,niData)
    ni_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ni_efold_primitive(x,f) - primEnd)
    endif
    
  end function ni_x_star

  function find_ni_x_star(x,niData)   
    implicit none
    real(kp) :: find_ni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: niData

    real(kp) :: primStar,f,w,CalFplusprimEnd,potStar,epsOneStar

    f = niData%real1
    w = niData%real2
    CalFplusprimEnd = niData%real3

    primStar = ni_efold_primitive(x,f)
    epsOneStar = ni_epsilon_one(x,f)
    potStar = ni_norm_potential(x,f)

    find_ni_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_ni_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ni_x_rrad(f,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ni_x_rrad
    real(kp), intent(in) :: f,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: niData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ni_x_endinf(f)
    epsOneEnd = ni_epsilon_one(xEnd,f)
    potEnd = ni_norm_potential(xEnd,f)
    primEnd = ni_efold_primitive(xEnd,f)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    niData%real1 = f
    niData%real2 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd

    x = zbrent(find_ni_x_rrad,mini,maxi,tolzbrent,niData)
    ni_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ni_efold_primitive(x,f) - primEnd)
    endif
    
  end function ni_x_rrad

  function find_ni_x_rrad(x,niData)   
    implicit none
    real(kp) :: find_ni_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: niData

    real(kp) :: primStar,f,CalFplusprimEnd,potStar,epsOneStar

    f = niData%real1
    CalFplusprimEnd = niData%real2

    primStar = ni_efold_primitive(x,f)
    epsOneStar = ni_epsilon_one(x,f)
    potStar = ni_norm_potential(x,f)

    find_ni_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)

  end function find_ni_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ni_x_rreh(f,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ni_x_rreh
    real(kp), intent(in) :: f,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: niData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ni_x_endinf(f)
    epsOneEnd = ni_epsilon_one(xEnd,f)
    potEnd = ni_norm_potential(xEnd,f)
    primEnd = ni_efold_primitive(xEnd,f)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    niData%real1 = f
    niData%real2 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd

    x = zbrent(find_ni_x_rreh,mini,maxi,tolzbrent,niData)
    ni_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ni_efold_primitive(x,f) - primEnd)
    endif
    
  end function ni_x_rreh

  function find_ni_x_rreh(x,niData)   
    implicit none
    real(kp) :: find_ni_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: niData

    real(kp) :: primStar,f,CalFplusprimEnd,potStar

    f = niData%real1
    CalFplusprimEnd = niData%real2

    primStar = ni_efold_primitive(x,f)
    potStar = ni_norm_potential(x,f)

    find_ni_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)

  end function find_ni_x_rreh




  function ni_lnrhoreh_max(f,Pstar) 
    implicit none
    real(kp) :: ni_lnrhoreh_max
    real(kp), intent(in) :: f,Pstar

    real(kp) :: xEnd, potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = ni_x_endinf(f)      
    potEnd  = ni_norm_potential(xEnd,f)
    epsEnd = ni_epsilon_one(xEnd,f)

!   Trick to return x such that rho_reh=rho_end
       
    x = ni_x_star(f,wrad,junk,Pstar)    
    potStar = ni_norm_potential(x,f)
    epsOneStar = ni_epsilon_one(x,f)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ni_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    ni_lnrhoreh_max = lnRhoEnd

  end function ni_lnrhoreh_max

  

    
end module nireheat

! S-dual inflation reheating functions in the slow-roll approximations

module sdireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sdisr, only : sdi_epsilon_one, sdi_epsilon_two, sdi_epsilon_three
  use sdisr, only : sdi_norm_potential
  use sdisr, only : sdi_efold_primitive
  implicit none

  private

  public sdi_x_star, sdi_lnrhoreh_max
  public sdi_x_rrad, sdi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function sdi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sdi_x_star
    real(kp), intent(in) :: phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sdiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sdi_epsilon_one(xEnd,phi0)
    potEnd = sdi_norm_potential(xEnd)
    primEnd = sdi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    sdiData%real1 = phi0 
    sdiData%real2 = xEnd
    sdiData%real3 = w
    sdiData%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_sdi_x_star,mini,maxi,tolzbrent,sdiData)
    sdi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (sdi_efold_primitive(x,phi0) - primEnd)
    endif

  end function sdi_x_star

  function find_sdi_x_star(x,sdiData)   
    implicit none
    real(kp) :: find_sdi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sdiData

    real(kp) :: primStar,phi0,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=sdiData%real1
    xEnd=sdiData%real2
    w = sdiData%real3
    CalFplusprimEnd = sdiData%real4

    primStar = sdi_efold_primitive(x,phi0)
    epsOneStar = sdi_epsilon_one(x,phi0)
    potStar = sdi_norm_potential(x)

    find_sdi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_sdi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sdi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sdi_x_rrad
    real(kp), intent(in) :: phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sdiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sdi_epsilon_one(xEnd,phi0)
    potEnd = sdi_norm_potential(xEnd)
    primEnd = sdi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    sdiData%real1 = phi0 
    sdiData%real2 = xEnd
    sdiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_sdi_x_rrad,mini,maxi,tolzbrent,sdiData)
    sdi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (sdi_efold_primitive(x,phi0) - primEnd)
    endif

  end function sdi_x_rrad

  function find_sdi_x_rrad(x,sdiData)   
    implicit none
    real(kp) :: find_sdi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sdiData

    real(kp) :: primStar,phi0,xEnd,CalFplusprimEnd,potStar,epsOneStar

    phi0=sdiData%real1
    xEnd=sdiData%real2
    CalFplusprimEnd = sdiData%real3

    primStar = sdi_efold_primitive(x,phi0)
    epsOneStar = sdi_epsilon_one(x,phi0)
    potStar = sdi_norm_potential(x)

    find_sdi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_sdi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sdi_x_rreh(phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: sdi_x_rreh
    real(kp), intent(in) :: phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sdiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = sdi_epsilon_one(xEnd,phi0)
    potEnd = sdi_norm_potential(xEnd)
    primEnd = sdi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    sdiData%real1 = phi0 
    sdiData%real2 = xEnd
    sdiData%real3 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_sdi_x_rreh,mini,maxi,tolzbrent,sdiData)
    sdi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (sdi_efold_primitive(x,phi0) - primEnd)
    endif

  end function sdi_x_rreh

  function find_sdi_x_rreh(x,sdiData)   
    implicit none
    real(kp) :: find_sdi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sdiData

    real(kp) :: primStar,phi0,xEnd,CalFplusprimEnd,potStar

    phi0=sdiData%real1
    xEnd=sdiData%real2
    CalFplusprimEnd = sdiData%real3

    primStar = sdi_efold_primitive(x,phi0)    
    potStar = sdi_norm_potential(x)

    find_sdi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_sdi_x_rreh




  function sdi_lnrhoreh_max(phi0,xend,Pstar) 
    implicit none
    real(kp) :: sdi_lnrhoreh_max
    real(kp), intent(in) :: phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
       
    potEnd  = sdi_norm_potential(xEnd)
    epsOneEnd = sdi_epsilon_one(xEnd,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = sdi_x_star(phi0,xEnd,wrad,junk,Pstar)    
    potStar = sdi_norm_potential(x)
    epsOneStar = sdi_epsilon_one(x,phi0)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'sdi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sdi_lnrhoreh_max = lnRhoEnd

  end function sdi_lnrhoreh_max

  
end module sdireheat

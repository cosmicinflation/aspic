!string axion I inflation common reheating functions in the slow-roll
!approximations

module saiicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use saiicommon, only : saii_epsilon_one
  use saiicommon, only : saii_norm_potential, saii_efold_primitive
  
  implicit none

  private
  
  public saii_x_star, saii_x_rrad, saii_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saii_x_star(alpha,mu,w,lnRhoReh,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saii_x_star
    real(kp), intent(in) :: alpha,mu,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saii_epsilon_one(xEnd,alpha,mu)
    potEnd = saii_norm_potential(xEnd,alpha,mu)

    primEnd = saii_efold_primitive(xEnd,alpha,mu)
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    saiiData%real1 = alpha
    saiiData%real2 = mu
    saiiData%real3 = w
    saiiData%real4 = calF + primEnd
    
    x = zbrent(find_saii_x_star,xmin,xmax,tolzbrent,saiiData)
    saii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (saii_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function saii_x_star

  function find_saii_x_star(x,saiiData)   
    implicit none
    real(kp) :: find_saii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData

    real(kp) :: primStar,alpha,mu,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=saiiData%real1
    mu = saiiData%real2
    w = saiiData%real3
    CalFplusprimEnd = saiiData%real4

    primStar = saii_efold_primitive(x,alpha,mu)
    epsOneStar = saii_epsilon_one(x,alpha,mu)
    potStar = saii_norm_potential(x,alpha,mu)

    find_saii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
    
  end function find_saii_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saii_x_rrad(alpha,mu,lnRrad,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saii_x_rrad
    real(kp), intent(in) :: alpha,mu,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saii_epsilon_one(xEnd,alpha,mu)
    potEnd = saii_norm_potential(xEnd,alpha,mu)

    primEnd = saii_efold_primitive(xEnd,alpha,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    saiiData%real1 = alpha
    saiiData%real2 = mu
    saiiData%real3 = calF + primEnd

    x = zbrent(find_saii_x_rrad,xmin,xmax,tolzbrent,saiiData)
    saii_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (saii_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function saii_x_rrad

  function find_saii_x_rrad(x,saiiData)   
    implicit none
    real(kp) :: find_saii_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData

    real(kp) :: primStar,alpha,mu,CalFplusprimEnd,potStar,epsOneStar

    alpha=saiiData%real1
    mu=saiiData%real2
    CalFplusprimEnd = saiiData%real3

    primStar = saii_efold_primitive(x,alpha,mu)
    epsOneStar = saii_epsilon_one(x,alpha,mu)
    potStar = saii_norm_potential(x,alpha,mu)

    find_saii_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_saii_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function saii_x_rreh(alpha,mu,lnRreh,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saii_x_rreh
    real(kp), intent(in) :: alpha,mu,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saii_epsilon_one(xEnd,alpha,mu)
    potEnd = saii_norm_potential(xEnd,alpha,mu)

    primEnd = saii_efold_primitive(xEnd,alpha,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    saiiData%real1 = alpha
    saiiData%real2 = mu
    saiiData%real3 = calF + primEnd

    x = zbrent(find_saii_x_rreh,xmin,xmax,tolzbrent,saiiData)
    saii_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (saii_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function saii_x_rreh

  function find_saii_x_rreh(x,saiiData)   
    implicit none
    real(kp) :: find_saii_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiData

    real(kp) :: primStar,alpha,mu,CalFplusprimEnd,potStar

    alpha=saiiData%real1
    mu=saiiData%real2
    CalFplusprimEnd = saiiData%real3

    primStar = saii_efold_primitive(x,alpha,mu)
    potStar = saii_norm_potential(x,alpha,mu)

    find_saii_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_saii_x_rreh


  
end module saiicomreh

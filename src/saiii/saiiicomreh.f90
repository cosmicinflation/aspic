!string axion II inflation common reheating functions in the slow-roll
!approximations

module saiiicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use saiiicommon, only : saiii_epsilon_one, saiiiXBig
  use saiiicommon, only : saiii_norm_potential, saiii_efold_primitive
  
  implicit none

  private

  public saiiiXBig
  public saiii_x_star, saiii_x_rrad, saiii_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saiii_x_star(alpha,beta,mu,w,lnRhoReh,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saiii_x_star
    real(kp), intent(in) :: alpha,beta,mu,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saiii_epsilon_one(xEnd,alpha,beta,mu)
    potEnd = saiii_norm_potential(xEnd,alpha,beta,mu)

    primEnd = saiii_efold_primitive(xEnd,alpha,beta,mu)
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = w
    saiiiData%real5 = calF + primEnd
    
    x = zbrent(find_saiii_x_star,xmin,xmax,tolzbrent,saiiiData)
    saiii_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (saiii_efold_primitive(x,alpha,beta,mu) - primEnd)
    endif

  end function saiii_x_star

  function find_saiii_x_star(x,saiiiData)   
    implicit none
    real(kp) :: find_saiii_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: primStar,alpha,beta,mu,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=saiiiData%real1
    beta=saiiiData%real2
    mu = saiiiData%real3
    w = saiiiData%real4
    CalFplusprimEnd = saiiiData%real5

    primStar = saiii_efold_primitive(x,alpha,beta,mu)
    epsOneStar = saiii_epsilon_one(x,alpha,beta,mu)
    potStar = saiii_norm_potential(x,alpha,beta,mu)

    find_saiii_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
    
  end function find_saiii_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saiii_x_rrad(alpha,beta,mu,lnRrad,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saiii_x_rrad
    real(kp), intent(in) :: alpha,beta,mu,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saiii_epsilon_one(xEnd,alpha,beta,mu)
    potEnd = saiii_norm_potential(xEnd,alpha,beta,mu)

    primEnd = saiii_efold_primitive(xEnd,alpha,beta,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = calF + primEnd

    x = zbrent(find_saiii_x_rrad,xmin,xmax,tolzbrent,saiiiData)
    saiii_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (saiii_efold_primitive(x,alpha,beta,mu) - primEnd)
    endif

  end function saiii_x_rrad

  function find_saiii_x_rrad(x,saiiiData)   
    implicit none
    real(kp) :: find_saiii_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: primStar,alpha,beta,mu,CalFplusprimEnd,potStar,epsOneStar

    alpha=saiiiData%real1
    beta=saiiiData%real2
    mu=saiiiData%real3
    CalFplusprimEnd = saiiiData%real4

    primStar = saiii_efold_primitive(x,alpha,beta,mu)
    epsOneStar = saiii_epsilon_one(x,alpha,beta,mu)
    potStar = saiii_norm_potential(x,alpha,beta,mu)

    find_saiii_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_saiii_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function saiii_x_rreh(alpha,beta,mu,lnRreh,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: saiii_x_rreh
    real(kp), intent(in) :: alpha,beta,mu,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: saiiiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = saiii_epsilon_one(xEnd,alpha,beta,mu)
    potEnd = saiii_norm_potential(xEnd,alpha,beta,mu)

    primEnd = saiii_efold_primitive(xEnd,alpha,beta,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    saiiiData%real1 = alpha
    saiiiData%real2 = beta
    saiiiData%real3 = mu
    saiiiData%real4 = calF + primEnd

    x = zbrent(find_saiii_x_rreh,xmin,xmax,tolzbrent,saiiiData)
    saiii_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (saiii_efold_primitive(x,alpha,beta,mu) - primEnd)
    endif

  end function saiii_x_rreh

  function find_saiii_x_rreh(x,saiiiData)   
    implicit none
    real(kp) :: find_saiii_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saiiiData

    real(kp) :: primStar,alpha,beta,mu,CalFplusprimEnd,potStar

    alpha=saiiiData%real1
    beta=saiiiData%real2
    mu=saiiiData%real3
    CalFplusprimEnd = saiiiData%real4

    primStar = saiii_efold_primitive(x,alpha,beta,mu)
    potStar = saiii_norm_potential(x,alpha,beta,mu)

    find_saiii_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_saiii_x_rreh


  
end module saiiicomreh

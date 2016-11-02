!Generalized Mixed large field inflation reheating functions in the
!slow-roll approximations

module rpqdireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rpqdisr, only : rpqdi_epsilon_one, rpqdi_epsilon_two, rpqdi_epsilon_three
  use rpqdisr, only : rpqdi_norm_potential, rpqdi_x_endinf
  use rpqdisr, only :  rpqdi_efold_primitive, rpqdi_xstar_brackets
  implicit none

  private

  public rpqdi_x_star, rpqdi_lnrhoreh_max
  public rpqdi_x_rrad, rpqdi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpqdi_x_star(phi0,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpqdi_x_star
    real(kp), intent(in) :: phi0,alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqdiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqdi_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rpqdi_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqdi_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqdi_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqdi_efold_primitive(xEnd,phi0,alpha,beta) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpqdiData%real1 = phi0
    rpqdiData%real2 = alpha
    rpqdiData%real3 = beta
    rpqdiData%real4 = xEnd
    rpqdiData%real5 = w
    rpqdiData%real6 = calF + primEnd

    x = zbrent(find_rpqdi_x_star,mini,maxi,tolzbrent,rpqdiData)
    rpqdi_x_star = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqdi_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqdi_x_star

  function find_rpqdi_x_star(x,rpqdiData)   
    implicit none
    real(kp) :: find_rpqdi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqdiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=rpqdiData%real1
    alpha=rpqdiData%real2
    beta=rpqdiData%real3
    xEnd=rpqdiData%real4
    w = rpqdiData%real5
    CalFplusprimEnd = rpqdiData%real6

    primStar = rpqdi_efold_primitive(x,phi0,alpha,beta)
    epsOneStar = rpqdi_epsilon_one(x,phi0,alpha,beta)
    potStar = rpqdi_norm_potential(x,phi0,alpha,beta)

    find_rpqdi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpqdi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rpqdi_x_rrad(phi0,alpha,beta,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpqdi_x_rrad
    real(kp), intent(in) :: phi0,alpha,beta,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqdiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqdi_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    xEnd = rpqdi_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqdi_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqdi_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqdi_efold_primitive(xEnd,phi0,alpha,beta) 
    
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rpqdiData%real1 = phi0
    rpqdiData%real2 = alpha
    rpqdiData%real3 = beta
    rpqdiData%real4 = xEnd
    rpqdiData%real5 = calF + primEnd

    x = zbrent(find_rpqdi_x_rrad,mini,maxi,tolzbrent,rpqdiData)
    rpqdi_x_rrad = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqdi_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqdi_x_rrad

  function find_rpqdi_x_rrad(x,rpqdiData)   
    implicit none
    real(kp) :: find_rpqdi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqdiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,CalFplusprimEnd,potStar,epsOneStar

    phi0=rpqdiData%real1
    alpha=rpqdiData%real2
    beta=rpqdiData%real3
    xEnd=rpqdiData%real4
    CalFplusprimEnd = rpqdiData%real5

    primStar = rpqdi_efold_primitive(x,phi0,alpha,beta)
    epsOneStar = rpqdi_epsilon_one(x,phi0,alpha,beta)
    potStar = rpqdi_norm_potential(x,phi0,alpha,beta)

    find_rpqdi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rpqdi_x_rrad

  
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rpqdi_x_rreh(phi0,alpha,beta,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rpqdi_x_rreh
    real(kp), intent(in) :: phi0,alpha,beta,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqdiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqdi_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    xEnd = rpqdi_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqdi_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqdi_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqdi_efold_primitive(xEnd,phi0,alpha,beta) 
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rpqdiData%real1 = phi0
    rpqdiData%real2 = alpha
    rpqdiData%real3 = beta
    rpqdiData%real4 = xEnd
    rpqdiData%real5 = calF + primEnd

    x = zbrent(find_rpqdi_x_rreh,mini,maxi,tolzbrent,rpqdiData)
    rpqdi_x_rreh = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqdi_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqdi_x_rreh

  function find_rpqdi_x_rreh(x,rpqdiData)   
    implicit none
    real(kp) :: find_rpqdi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqdiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,CalFplusprimEnd,potStar

    phi0=rpqdiData%real1
    alpha=rpqdiData%real2
    beta=rpqdiData%real3
    xEnd=rpqdiData%real4
    CalFplusprimEnd = rpqdiData%real5

    primStar = rpqdi_efold_primitive(x,phi0,alpha,beta)
    potStar = rpqdi_norm_potential(x,phi0,alpha,beta)

    find_rpqdi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rpqdi_x_rreh


  function rpqdi_lnrhoreh_max(phi0,alpha,beta,Pstar) 
    implicit none
    real(kp) :: rpqdi_lnrhoreh_max
    real(kp), intent(in) :: phi0,alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = rpqdi_x_endinf(phi0,alpha,beta)
    potEnd  = rpqdi_norm_potential(xEnd,phi0,alpha,beta)
    epsOneEnd = rpqdi_epsilon_one(xEnd,phi0,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = rpqdi_x_star(phi0,alpha,beta,wrad,junk,Pstar)    
    potStar = rpqdi_norm_potential(x,phi0,alpha,beta)
    epsOneStar = rpqdi_epsilon_one(x,phi0,alpha,beta)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rpqdi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpqdi_lnrhoreh_max = lnRhoEnd

  end function rpqdi_lnrhoreh_max

  
end module rpqdireheat

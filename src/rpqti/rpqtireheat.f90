!Generalized Mixed large field inflation reheating functions in the
!slow-roll approximations

module rpqtireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rpqtisr, only : rpqti_epsilon_one, rpqti_epsilon_two, rpqti_epsilon_three
  use rpqtisr, only : rpqti_norm_potential, rpqti_x_endinf
  use rpqtisr, only :  rpqti_efold_primitive, rpqti_xstar_brackets
  implicit none

  private

  public rpqti_x_star, rpqti_lnrhoreh_max
  public rpqti_x_rrad, rpqti_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rpqti_x_star(phi0,alpha,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpqti_x_star
    real(kp), intent(in) :: phi0,alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqtiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqti_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rpqti_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqti_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqti_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqti_efold_primitive(xEnd,phi0,alpha,beta) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rpqtiData%real1 = phi0
    rpqtiData%real2 = alpha
    rpqtiData%real3 = beta
    rpqtiData%real4 = xEnd
    rpqtiData%real5 = w
    rpqtiData%real6 = calF + primEnd

    x = zbrent(find_rpqti_x_star,mini,maxi,tolzbrent,rpqtiData)
    rpqti_x_star = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqti_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqti_x_star

  function find_rpqti_x_star(x,rpqtiData)   
    implicit none
    real(kp) :: find_rpqti_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqtiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=rpqtiData%real1
    alpha=rpqtiData%real2
    beta=rpqtiData%real3
    xEnd=rpqtiData%real4
    w = rpqtiData%real5
    CalFplusprimEnd = rpqtiData%real6

    primStar = rpqti_efold_primitive(x,phi0,alpha,beta)
    epsOneStar = rpqti_epsilon_one(x,phi0,alpha,beta)
    potStar = rpqti_norm_potential(x,phi0,alpha,beta)

    find_rpqti_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rpqti_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rpqti_x_rrad(phi0,alpha,beta,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rpqti_x_rrad
    real(kp), intent(in) :: phi0,alpha,beta,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqtiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqti_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    xEnd = rpqti_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqti_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqti_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqti_efold_primitive(xEnd,phi0,alpha,beta) 
    
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rpqtiData%real1 = phi0
    rpqtiData%real2 = alpha
    rpqtiData%real3 = beta
    rpqtiData%real4 = xEnd
    rpqtiData%real5 = calF + primEnd

    x = zbrent(find_rpqti_x_rrad,mini,maxi,tolzbrent,rpqtiData)
    rpqti_x_rrad = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqti_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqti_x_rrad

  function find_rpqti_x_rrad(x,rpqtiData)   
    implicit none
    real(kp) :: find_rpqti_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqtiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,CalFplusprimEnd,potStar,epsOneStar

    phi0=rpqtiData%real1
    alpha=rpqtiData%real2
    beta=rpqtiData%real3
    xEnd=rpqtiData%real4
    CalFplusprimEnd = rpqtiData%real5

    primStar = rpqti_efold_primitive(x,phi0,alpha,beta)
    epsOneStar = rpqti_epsilon_one(x,phi0,alpha,beta)
    potStar = rpqti_norm_potential(x,phi0,alpha,beta)

    find_rpqti_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rpqti_x_rrad

  
!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rpqti_x_rreh(phi0,alpha,beta,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rpqti_x_rreh
    real(kp), intent(in) :: phi0,alpha,beta,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rpqtiData

    real(kp), dimension(1:2) :: xstar_brackets

    xstar_brackets = rpqti_xstar_brackets(phi0,alpha,beta)
    mini = xstar_brackets(1)
    maxi = xstar_brackets(2)
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    xEnd = rpqti_x_endinf(phi0,alpha,beta)
    epsOneEnd = rpqti_epsilon_one(xEnd,phi0,alpha,beta)
    potEnd = rpqti_norm_potential(xEnd,phi0,alpha,beta)
    primEnd = rpqti_efold_primitive(xEnd,phi0,alpha,beta) 
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rpqtiData%real1 = phi0
    rpqtiData%real2 = alpha
    rpqtiData%real3 = beta
    rpqtiData%real4 = xEnd
    rpqtiData%real5 = calF + primEnd

    x = zbrent(find_rpqti_x_rreh,mini,maxi,tolzbrent,rpqtiData)
    rpqti_x_rreh = x  

    if (present(bfoldstar)) then
       bfoldstar = - (rpqti_efold_primitive(x,phi0,alpha,beta) - primEnd)
    endif

  end function rpqti_x_rreh

  function find_rpqti_x_rreh(x,rpqtiData)   
    implicit none
    real(kp) :: find_rpqti_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rpqtiData

    real(kp) :: primStar,phi0,alpha,beta,xEnd,CalFplusprimEnd,potStar

    phi0=rpqtiData%real1
    alpha=rpqtiData%real2
    beta=rpqtiData%real3
    xEnd=rpqtiData%real4
    CalFplusprimEnd = rpqtiData%real5

    primStar = rpqti_efold_primitive(x,phi0,alpha,beta)
    potStar = rpqti_norm_potential(x,phi0,alpha,beta)

    find_rpqti_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rpqti_x_rreh


  function rpqti_lnrhoreh_max(phi0,alpha,beta,Pstar) 
    implicit none
    real(kp) :: rpqti_lnrhoreh_max
    real(kp), intent(in) :: phi0,alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = rpqti_x_endinf(phi0,alpha,beta)
    potEnd  = rpqti_norm_potential(xEnd,phi0,alpha,beta)
    epsOneEnd = rpqti_epsilon_one(xEnd,phi0,alpha,beta)

!   Trick to return x such that rho_reh=rho_end

    x = rpqti_x_star(phi0,alpha,beta,wrad,junk,Pstar)    
    potStar = rpqti_norm_potential(x,phi0,alpha,beta)
    epsOneStar = rpqti_epsilon_one(x,phi0,alpha,beta)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rpqti_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rpqti_lnrhoreh_max = lnRhoEnd

  end function rpqti_lnrhoreh_max

  
end module rpqtireheat

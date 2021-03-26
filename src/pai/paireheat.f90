! pure arctan inflation reheating functions in the slow-roll approximations

module paireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use paisr, only : pai_epsilon_one, pai_epsilon_two, pai_epsilon_three
  use paisr, only : pai_norm_potential, pai_efold_primitive, pai_numacc_xinimax
  implicit none

  private

  public pai_x_star, pai_lnrhoreh_max
  public pai_x_rrad, pai_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function pai_x_star(mu,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pai_x_star
    real(kp), intent(in) :: mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: paiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = pai_epsilon_one(xEnd,mu)
    potEnd = pai_norm_potential(xEnd,mu)
    primEnd = pai_efold_primitive(xEnd,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    paiData%real1 = mu 
    paiData%real2 = xEnd
    paiData%real3 = w
    paiData%real4 = calF + primEnd

    mini = xend
    maxi = pai_numacc_xinimax(mu)


    x = zbrent(find_pai_x_star,mini,maxi,tolzbrent,paiData)
    pai_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (pai_efold_primitive(x,mu) - primEnd)
    endif

  end function pai_x_star

  function find_pai_x_star(x,paiData)   
    implicit none
    real(kp) :: find_pai_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: paiData

    real(kp) :: primStar,mu,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    mu=paiData%real1
    xEnd=paiData%real2
    w = paiData%real3
    CalFplusprimEnd = paiData%real4

    primStar = pai_efold_primitive(x,mu)
    epsOneStar = pai_epsilon_one(x,mu)
    potStar = pai_norm_potential(x,mu)

    find_pai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_pai_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function pai_x_rrad(mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pai_x_rrad
    real(kp), intent(in) :: mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: paiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = pai_epsilon_one(xEnd,mu)
    potEnd = pai_norm_potential(xEnd,mu)
    primEnd = pai_efold_primitive(xEnd,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    paiData%real1 = mu 
    paiData%real2 = xEnd
    paiData%real3 = calF + primEnd

    mini = xEnd
    maxi = pai_numacc_xinimax(mu)


    x = zbrent(find_pai_x_rrad,mini,maxi,tolzbrent,paiData)
    pai_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (pai_efold_primitive(x,mu) - primEnd)
    endif

  end function pai_x_rrad

  function find_pai_x_rrad(x,paiData)   
    implicit none
    real(kp) :: find_pai_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: paiData

    real(kp) :: primStar,mu,xEnd,CalFplusprimEnd,potStar,epsOneStar

    mu=paiData%real1
    xEnd=paiData%real2
    CalFplusprimEnd = paiData%real3

    primStar = pai_efold_primitive(x,mu)
    epsOneStar = pai_epsilon_one(x,mu)
    potStar = pai_norm_potential(x,mu)

    find_pai_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_pai_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function pai_x_rreh(mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: pai_x_rreh
    real(kp), intent(in) :: mu,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: paiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = pai_epsilon_one(xEnd,mu)
    potEnd = pai_norm_potential(xEnd,mu)
    primEnd = pai_efold_primitive(xEnd,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    paiData%real1 = mu 
    paiData%real2 = xEnd
    paiData%real3 = calF + primEnd

    mini = xEnd
    maxi = pai_numacc_xinimax(mu)


    x = zbrent(find_pai_x_rreh,mini,maxi,tolzbrent,paiData)
    pai_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (pai_efold_primitive(x,mu) - primEnd)
    endif

  end function pai_x_rreh

  function find_pai_x_rreh(x,paiData)   
    implicit none
    real(kp) :: find_pai_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: paiData

    real(kp) :: primStar,mu,xEnd,CalFplusprimEnd,potStar

    mu=paiData%real1
    xEnd=paiData%real2
    CalFplusprimEnd = paiData%real3

    primStar = pai_efold_primitive(x,mu)    
    potStar = pai_norm_potential(x,mu)

    find_pai_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_pai_x_rreh




  function pai_lnrhoreh_max(mu,xend,Pstar) 
    implicit none
    real(kp) :: pai_lnrhoreh_max
    real(kp), intent(in) :: mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
       
    potEnd  = pai_norm_potential(xEnd,mu)
    epsOneEnd = pai_epsilon_one(xEnd,mu)


!   Trick to return x such that rho_reh=rho_end

    x = pai_x_star(mu,xEnd,wrad,junk,Pstar)    
    potStar = pai_norm_potential(x,mu)
    epsOneStar = pai_epsilon_one(x,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'pai_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    pai_lnrhoreh_max = lnRhoEnd

  end function pai_lnrhoreh_max

  
end module paireheat

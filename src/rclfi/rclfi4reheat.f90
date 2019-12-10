!Radiatively Corrected Large Field Inflation 4 reheating functions in the
!slow-roll approximations

module rclfi4reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rclfi4sr, only : rclfi4_epsilon_one, rclfi4_epsilon_two, rclfi4_epsilon_three
  use rclfi4sr, only : rclfi4_norm_potential, rclfi4_x_endinf, rclfi4_efold_primitive
  use rclfi4sr, only : rclfixBig
  
  implicit none

  private

  public rclfi4_x_star, rclfi4_lnrhoreh_max
  public rclfi4_x_rrad, rclfi4_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rclfi4_x_star(alpha,p,mu,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rclfi4_x_star
    real(kp), intent(in) :: alpha,p,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rclfi4_epsilon_one(xEnd,alpha,p,mu)
    potEnd = rclfi4_norm_potential(xEnd,alpha,p,mu)
    primEnd = rclfi4_efold_primitive(xEnd,alpha,p,mu) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = w
    rclfiData%real5 = calF + primEnd

    mini = xend
    maxi = rclfixBig
    
    rclfi4_x_star = zbrent(find_rclfi4_x_star,mini,maxi,tolzbrent,rclfiData)
    
    if (present(bfoldstar)) then
       x = rclfi4_x_star
       bfoldstar = - (rclfi4_efold_primitive(x,alpha,p,mu) - primEnd)
    endif

  end function rclfi4_x_star

  function find_rclfi4_x_star(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi4_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,alpha,p,mu,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    w = rclfiData%real4
    CalFplusprimEnd = rclfiData%real5

    primStar = rclfi4_efold_primitive(x,alpha,p,mu)
    epsOneStar = rclfi4_epsilon_one(x,alpha,p,mu)
    potStar = rclfi4_norm_potential(x,alpha,p,mu)

    find_rclfi4_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
    
  end function find_rclfi4_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rclfi4_x_rrad(alpha,p,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rclfi4_x_rrad
    real(kp), intent(in) :: alpha,p,mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData


    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rclfi4_epsilon_one(xEnd,alpha,p,mu)
    potEnd = rclfi4_norm_potential(xEnd,alpha,p,mu)
    primEnd = rclfi4_efold_primitive(xEnd,alpha,p,mu) 
    

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = calF + primEnd

    mini = xend
    maxi = rclfixBig

    rclfi4_x_rrad = zbrent(find_rclfi4_x_rrad,mini,maxi,tolzbrent,rclfiData)


    if (present(bfoldstar)) then
       x = rclfi4_x_rrad
       bfoldstar = - (rclfi4_efold_primitive(x,alpha,p,mu) - primEnd)
    endif

  end function rclfi4_x_rrad

  function find_rclfi4_x_rrad(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi4_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,alpha,p,mu,CalFplusprimEnd,potStar,epsOneStar


    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    CalFplusprimEnd = rclfiData%real4

    primStar = rclfi4_efold_primitive(x,alpha,p,mu)
    epsOneStar = rclfi4_epsilon_one(x,alpha,p,mu)
    potStar = rclfi4_norm_potential(x,alpha,p,mu)

    find_rclfi4_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rclfi4_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rclfi4_x_rreh(alpha,p,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rclfi4_x_rreh
    real(kp), intent(in) :: alpha,p,mu,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData


    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rclfi4_epsilon_one(xEnd,alpha,p,mu)
    potEnd = rclfi4_norm_potential(xEnd,alpha,p,mu)
    primEnd = rclfi4_efold_primitive(xEnd,alpha,p,mu) 
    

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = calF + primEnd

    mini = xend
    maxi = rclfixBig

    rclfi4_x_rreh = zbrent(find_rclfi4_x_rreh,mini,maxi,tolzbrent,rclfiData)

    if (present(bfoldstar)) then
       x = rclfi4_x_rreh
       bfoldstar = - (rclfi4_efold_primitive(x,alpha,p,mu) - primEnd)
    endif

  end function rclfi4_x_rreh

  function find_rclfi4_x_rreh(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi4_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,alpha,p,mu,CalFplusprimEnd,potStar


    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    CalFplusprimEnd = rclfiData%real4

    primStar = rclfi4_efold_primitive(x,alpha,p,mu)
    potStar = rclfi4_norm_potential(x,alpha,p,mu)

    find_rclfi4_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rclfi4_x_rreh



  function rclfi4_lnrhoreh_max(alpha,p,mu,xend,Pstar) 
    implicit none
    real(kp) :: rclfi4_lnrhoreh_max
    real(kp), intent(in) :: alpha,p,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = rclfi4_norm_potential(xEnd,alpha,p,mu)
    epsOneEnd = rclfi4_epsilon_one(xEnd,alpha,p,mu)

!   Trick to return x such that rho_reh=rho_end

    x = rclfi4_x_star(alpha,p,mu,xend,wrad,junk,Pstar)    

    potStar = rclfi4_norm_potential(x,alpha,p,mu)
    epsOneStar = rclfi4_epsilon_one(x,alpha,p,mu)

    
    if (.not.slowroll_validity(epsOneStar)) then
       write(*,*)'WARNING:'
       write(*,*)'rclfi4_lnrhoreh_max: slow-roll violated!'
    endif
       
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rclfi4_lnrhoreh_max = lnRhoEnd

  end function rclfi4_lnrhoreh_max

  
end module rclfi4reheat

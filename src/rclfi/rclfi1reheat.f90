!Radiatively Corrected Large Field Inflation 1 reheating functions in the
!slow-roll approximations

module rclfi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rclfi1sr, only : rclfi1_epsilon_one, rclfi1_epsilon_two, rclfi1_epsilon_three
  use rclfi1sr, only : rclfi1_norm_potential, rclfi1_x_endinf, rclfi1_efold_primitive
  use rclfi1sr, only : rclfi1_numacc_xinimax
  implicit none

  private

  public rclfi1_x_star, rclfi1_lnrhoreh_max
  public rclfi1_x_rrad, rclfi1_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rclfi1_x_star(p,alpha,mu,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rclfi1_x_star
    real(kp), intent(in) :: p,alpha,mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = rclfi1_epsilon_one(xEnd,p,alpha,mu)
    potEnd = rclfi1_norm_potential(xEnd,p,alpha,mu)
    primEnd = rclfi1_efold_primitive(xEnd,p,alpha,mu) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = w
    rclfiData%real5 = calF + primEnd

    mini = xend
    maxi = rclfi1_numacc_xinimax(p,alpha,mu)

    rclfi1_x_star = zbrent(find_rclfi1_x_star,mini,maxi,tolzbrent,rclfiData)
    
    if (present(bfoldstar)) then
       x = rclfi1_x_star
       bfoldstar = - (rclfi1_efold_primitive(x,p,alpha,mu) - primEnd)
    endif

  end function rclfi1_x_star

  function find_rclfi1_x_star(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,p,alpha,mu,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    w = rclfiData%real4
    CalFplusprimEnd = rclfiData%real5

    primStar = rclfi1_efold_primitive(x,p,alpha,mu)
    epsOneStar = rclfi1_epsilon_one(x,p,alpha,mu)
    potStar = rclfi1_norm_potential(x,p,alpha,mu)

    find_rclfi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rclfi1_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rclfi1_x_rrad(p,alpha,mu,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rclfi1_x_rrad
    real(kp), intent(in) :: p,alpha,mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData


    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rclfi1_epsilon_one(xEnd,p,alpha,mu)
    potEnd = rclfi1_norm_potential(xEnd,p,alpha,mu)
    primEnd = rclfi1_efold_primitive(xEnd,p,alpha,mu) 
    

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = calF + primEnd

    mini = xend
    maxi = rclfi1_numacc_xinimax(p,alpha,mu)

    rclfi1_x_rrad = zbrent(find_rclfi1_x_rrad,mini,maxi,tolzbrent,rclfiData)


    if (present(bfoldstar)) then
       x = rclfi1_x_rrad
       bfoldstar = - (rclfi1_efold_primitive(x,p,alpha,mu) - primEnd)
    endif

  end function rclfi1_x_rrad

  function find_rclfi1_x_rrad(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi1_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,p,alpha,mu,CalFplusprimEnd,potStar,epsOneStar


    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    CalFplusprimEnd = rclfiData%real4

    primStar = rclfi1_efold_primitive(x,p,alpha,mu)
    epsOneStar = rclfi1_epsilon_one(x,p,alpha,mu)
    potStar = rclfi1_norm_potential(x,p,alpha,mu)

    find_rclfi1_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rclfi1_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rclfi1_x_rreh(p,alpha,mu,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rclfi1_x_rreh
    real(kp), intent(in) :: p,alpha,mu,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rclfiData


    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
        
    epsOneEnd = rclfi1_epsilon_one(xEnd,p,alpha,mu)
    potEnd = rclfi1_norm_potential(xEnd,p,alpha,mu)
    primEnd = rclfi1_efold_primitive(xEnd,p,alpha,mu) 
    

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)


    rclfiData%real1 = alpha
    rclfiData%real2 = p
    rclfiData%real3 = mu
    rclfiData%real4 = calF + primEnd

    mini = xend
    maxi = rclfi1_numacc_xinimax(p,alpha,mu)

    rclfi1_x_rreh = zbrent(find_rclfi1_x_rreh,mini,maxi,tolzbrent,rclfiData)

    if (present(bfoldstar)) then
       x = rclfi1_x_rreh
       bfoldstar = - (rclfi1_efold_primitive(x,p,alpha,mu) - primEnd)
    endif

  end function rclfi1_x_rreh

  function find_rclfi1_x_rreh(x,rclfiData)   
    implicit none
    real(kp) :: find_rclfi1_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rclfiData

    real(kp) :: primStar,p,alpha,mu,CalFplusprimEnd,potStar


    alpha = rclfiData%real1
    p = rclfiData%real2
    mu = rclfiData%real3
    CalFplusprimEnd = rclfiData%real4

    primStar = rclfi1_efold_primitive(x,p,alpha,mu)
    potStar = rclfi1_norm_potential(x,p,alpha,mu)

    find_rclfi1_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rclfi1_x_rreh



  function rclfi1_lnrhoreh_max(p,alpha,mu,xend,Pstar) 
    implicit none
    real(kp) :: rclfi1_lnrhoreh_max
    real(kp), intent(in) :: p,alpha,mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = rclfi1_norm_potential(xEnd,p,alpha,mu)
    epsOneEnd = rclfi1_epsilon_one(xEnd,p,alpha,mu)

!   Trick to return x such that rho_reh=rho_end

    x = rclfi1_x_star(p,alpha,mu,xend,wrad,junk,Pstar)    
    potStar = rclfi1_norm_potential(x,p,alpha,mu)
    epsOneStar = rclfi1_epsilon_one(x,p,alpha,mu)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rclfi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rclfi1_lnrhoreh_max = lnRhoEnd

  end function rclfi1_lnrhoreh_max

  
end module rclfi1reheat

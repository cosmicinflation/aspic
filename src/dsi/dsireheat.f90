!dynamical supersymmetric inflation reheating functions in the slow-roll approximations

module dsireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use dsisr, only : dsi_epsilon_one, dsi_epsilon_two, dsi_epsilon_three
  use dsisr, only : dsi_norm_potential
  use dsisr, only : dsi_xinimin, dsi_efold_primitive
  implicit none

  private

  public dsi_x_star, dsi_lnrhoreh_max 
  public dsi_x_rrad, dsi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function dsi_x_star(p,mu,xEnd,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: dsi_x_star
    real(kp), intent(in) :: p,mu,xEnd,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: dsiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = dsi_epsilon_one(xEnd,p,mu)
    potEnd = dsi_norm_potential(xEnd,p,mu)
    primEnd = dsi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    dsiData%real1 = p 
    dsiData%real2 = mu
    dsiData%real3 = xEnd
    dsiData%real4 = w
    dsiData%real5 = calF + primEnd

    mini=dsi_xinimin(p,mu)
    maxi = xend

    x = zbrent(find_dsi_x_star,mini,maxi,tolzbrent,dsiData)
    dsi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (dsi_efold_primitive(x,p,mu) - primEnd)
    endif

  end function dsi_x_star

  function find_dsi_x_star(x,dsiData)   
    implicit none
    real(kp) :: find_dsi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: dsiData

    real(kp) :: primStar,p,mu,xEnd,w,CalFplusprimEnd,potStar,epsOneStar

    p=dsiData%real1
    mu=dsiData%real2
    xEnd=dsiData%real3
    w = dsiData%real4
    CalFplusprimEnd = dsiData%real5

    primStar = dsi_efold_primitive(x,p,mu)
    epsOneStar = dsi_epsilon_one(x,p,mu)
    potStar = dsi_norm_potential(x,p,mu)

    find_dsi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_dsi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function dsi_x_rrad(p,mu,xEnd,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: dsi_x_rrad
    real(kp), intent(in) :: p,mu,xEnd,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: dsiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = dsi_epsilon_one(xEnd,p,mu)
    potEnd = dsi_norm_potential(xEnd,p,mu)
    primEnd = dsi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    dsiData%real1 = p 
    dsiData%real2 = mu
    dsiData%real3 = xEnd
    dsiData%real4 = calF + primEnd

    mini=dsi_xinimin(p,mu)
    maxi = xend

    x = zbrent(find_dsi_x_rrad,mini,maxi,tolzbrent,dsiData)
    dsi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (dsi_efold_primitive(x,p,mu) - primEnd)
    endif

  end function dsi_x_rrad

  function find_dsi_x_rrad(x,dsiData)   
    implicit none
    real(kp) :: find_dsi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: dsiData

    real(kp) :: primStar,p,mu,xEnd,CalFplusprimEnd,potStar,epsOneStar

    p=dsiData%real1
    mu=dsiData%real2
    xEnd=dsiData%real3
    CalFplusprimEnd = dsiData%real4

    primStar = dsi_efold_primitive(x,p,mu)
    epsOneStar = dsi_epsilon_one(x,p,mu)
    potStar = dsi_norm_potential(x,p,mu)

    find_dsi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_dsi_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function dsi_x_rreh(p,mu,xEnd,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: dsi_x_rreh
    real(kp), intent(in) :: p,mu,xEnd,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: dsiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = dsi_epsilon_one(xEnd,p,mu)
    potEnd = dsi_norm_potential(xEnd,p,mu)
    primEnd = dsi_efold_primitive(xEnd,p,mu)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    dsiData%real1 = p 
    dsiData%real2 = mu
    dsiData%real3 = xEnd
    dsiData%real4 = calF + primEnd

    mini=dsi_xinimin(p,mu)
    maxi = xend

    x = zbrent(find_dsi_x_rreh,mini,maxi,tolzbrent,dsiData)
    dsi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (dsi_efold_primitive(x,p,mu) - primEnd)
    endif

  end function dsi_x_rreh

  function find_dsi_x_rreh(x,dsiData)   
    implicit none
    real(kp) :: find_dsi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: dsiData

    real(kp) :: primStar,p,mu,xEnd,CalFplusprimEnd,potStar

    p=dsiData%real1
    mu=dsiData%real2
    xEnd=dsiData%real3
    CalFplusprimEnd = dsiData%real4

    primStar = dsi_efold_primitive(x,p,mu)    
    potStar = dsi_norm_potential(x,p,mu)

    find_dsi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_dsi_x_rreh




  function dsi_lnrhoreh_max(p,mu,xEnd,Pstar) 
    implicit none
    real(kp) :: dsi_lnrhoreh_max
    real(kp), intent(in) :: p,mu,xEnd,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    potEnd  = dsi_norm_potential(xEnd,p,mu)
    epsOneEnd = dsi_epsilon_one(xEnd,p,mu)

!   Trick to return x such that rho_reh=rho_end

    x = dsi_x_star(p,mu,xEnd,wrad,junk,Pstar)    
    potStar = dsi_norm_potential(x,p,mu)
    epsOneStar = dsi_epsilon_one(x,p,mu)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'dsi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    dsi_lnrhoreh_max = lnRhoEnd

  end function dsi_lnrhoreh_max

  
end module dsireheat

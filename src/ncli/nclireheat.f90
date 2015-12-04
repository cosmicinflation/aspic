!Non-Renormalizable Corrected Loop Inflation reheating functions in the slow-roll approximations

module nclireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use nclisr, only : ncli_epsilon_one, ncli_epsilon_two, ncli_epsilon_three
  use nclisr, only : ncli_norm_potential, ncli_efold_primitive, ncli_x_endinf
  use nclisr, only : xNumAccMax
  implicit none

  private

  public ncli_x_star, ncli_lnrhoreh_max 
  public ncli_x_rrad, ncli_x_rreh

contains

!returns x such potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ncli_x_star(alpha,phi0,n,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ncli_x_star
    real(kp), intent(in) :: alpha,phi0,n,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xEnd
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ncliData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ncli_x_endinf(alpha,phi0,n)
    epsOneEnd = ncli_epsilon_one(xEnd,alpha,phi0,n)
    potEnd = ncli_norm_potential(xEnd,alpha,phi0,n)
    primEnd = ncli_efold_primitive(xEnd,alpha,phi0,n)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    ncliData%real1 = alpha 
    ncliData%real2 = phi0
    ncliData%real3 = n
    ncliData%real4 = w
    ncliData%real5 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xNumAccMax

    x = zbrent(find_ncli_x_star,mini,maxi,tolzbrent,ncliData)
    ncli_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ncli_efold_primitive(x,alpha,phi0,n) - primEnd)
    endif

  end function ncli_x_star

  function find_ncli_x_star(x,ncliData)   
    implicit none
    real(kp) :: find_ncli_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData

    real(kp) :: primStar,alpha,phi0,n,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ncliData%real1
    phi0=ncliData%real2
    n=ncliData%real3
    w = ncliData%real4
    CalFplusprimEnd = ncliData%real5

    primStar = ncli_efold_primitive(x,alpha,phi0,n)
    epsOneStar = ncli_epsilon_one(x,alpha,phi0,n)
    potStar = ncli_norm_potential(x,alpha,phi0,n)

    find_ncli_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ncli_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ncli_x_rrad(alpha,phi0,n,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ncli_x_rrad
    real(kp), intent(in) :: alpha,phi0,n,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xEnd
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ncliData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ncli_x_endinf(alpha,phi0,n)
    epsOneEnd = ncli_epsilon_one(xEnd,alpha,phi0,n)


    potEnd = ncli_norm_potential(xEnd,alpha,phi0,n)
    primEnd = ncli_efold_primitive(xEnd,alpha,phi0,n)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    ncliData%real1 = alpha 
    ncliData%real2 = phi0
    ncliData%real3 = n
    ncliData%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = xNumAccMax

    x = zbrent(find_ncli_x_rrad,mini,maxi,tolzbrent,ncliData)
    ncli_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ncli_efold_primitive(x,alpha,phi0,n) - primEnd)
    endif

  end function ncli_x_rrad

  function find_ncli_x_rrad(x,ncliData)   
    implicit none
    real(kp) :: find_ncli_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData

    real(kp) :: primStar,alpha,phi0,n,CalFplusprimEnd,potStar,epsOneStar

    alpha=ncliData%real1
    phi0=ncliData%real2
    n=ncliData%real3
    CalFplusprimEnd = ncliData%real4

    primStar = ncli_efold_primitive(x,alpha,phi0,n)
    epsOneStar = ncli_epsilon_one(x,alpha,phi0,n)
    potStar = ncli_norm_potential(x,alpha,phi0,n)

    find_ncli_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_ncli_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ncli_x_rreh(alpha,phi0,n,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: ncli_x_rreh
    real(kp), intent(in) :: alpha,phi0,n,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xEnd
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: ncliData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ncli_x_endinf(alpha,phi0,n)
    epsOneEnd = ncli_epsilon_one(xEnd,alpha,phi0,n)
    potEnd = ncli_norm_potential(xEnd,alpha,phi0,n)
    primEnd = ncli_efold_primitive(xEnd,alpha,phi0,n)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    ncliData%real1 = alpha 
    ncliData%real2 = phi0
    ncliData%real3 = n
    ncliData%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = xNumAccMax

    x = zbrent(find_ncli_x_rreh,mini,maxi,tolzbrent,ncliData)
    ncli_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ncli_efold_primitive(x,alpha,phi0,n) - primEnd)
    endif

  end function ncli_x_rreh

  function find_ncli_x_rreh(x,ncliData)   
    implicit none
    real(kp) :: find_ncli_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ncliData

    real(kp) :: primStar,alpha,phi0,n,CalFplusprimEnd,potStar

    alpha=ncliData%real1
    phi0=ncliData%real2
    n=ncliData%real3
    CalFplusprimEnd = ncliData%real4

    primStar = ncli_efold_primitive(x,alpha,phi0,n)
    potStar = ncli_norm_potential(x,alpha,phi0,n)

    find_ncli_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_ncli_x_rreh



  function ncli_lnrhoreh_max(alpha,phi0,n,Pstar) 
    implicit none
    real(kp) :: ncli_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,n,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar, xEnd

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd

    xEnd = ncli_x_endinf(alpha,phi0,n)
    potEnd  = ncli_norm_potential(xEnd,alpha,phi0,n)
    epsOneEnd = ncli_epsilon_one(xEnd,alpha,phi0,n)

!   Trick to return x such that rho_reh=rho_end

    x = ncli_x_star(alpha,phi0,n,wrad,junk,Pstar)    
    potStar = ncli_norm_potential(x,alpha,phi0,n)
    epsOneStar = ncli_epsilon_one(x,alpha,phi0,n)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'ncli_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ncli_lnrhoreh_max = lnRhoEnd

  end function ncli_lnrhoreh_max

  
end module nclireheat

!tip inflation reheating functions in the slow-roll approximations

module tireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use tisr, only : ti_epsilon_one, ti_epsilon_two, ti_epsilon_three
  use tisr, only : ti_norm_potential,ti_efold_primitive,ti_x_endinf,ti_x_potmax
  implicit none

  private

  public ti_x_star, ti_lnrhoreh_max
  public ti_x_rrad, ti_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function ti_x_star(alpha,mu,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: ti_x_star
    real(kp), intent(in) :: alpha,mu,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: tiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ti_x_endinf(alpha,mu)
    
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)
    potEnd = ti_norm_potential(xEnd,alpha)
    primEnd = ti_efold_primitive(xEnd,alpha,mu)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    tiData%real1 = alpha
    tiData%real2 = mu
    tiData%real3 = w
    tiData%real4 = calF + primEnd

    mini=ti_x_potmax(alpha,mu) *(1._kp+epsilon(1._kp)) !potential maximum
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_ti_x_star,mini,maxi,tolFind,tiData)
    ti_x_star = x

    if (present(bfold)) then
       bfold = -(ti_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function ti_x_star

  function find_ti_x_star(x,tiData)   
    implicit none
    real(kp) :: find_ti_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: tiData

    real(kp) :: primStar,alpha,mu,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=tiData%real1
    mu=tiData%real2
    w = tiData%real3
    CalFplusPrimEnd = tiData%real4

    primStar = ti_efold_primitive(x,alpha,mu)
    epsOneStar = ti_epsilon_one(x,alpha,mu)
    potStar = ti_norm_potential(x,alpha)

    find_ti_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_ti_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ti_x_rrad(alpha,mu,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: ti_x_rrad
    real(kp), intent(in) :: alpha,mu,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: tiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ti_x_endinf(alpha,mu)
    
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)
    potEnd = ti_norm_potential(xEnd,alpha)
    primEnd = ti_efold_primitive(xEnd,alpha,mu)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    tiData%real1 = alpha
    tiData%real2 = mu
    tiData%real3 = calF + primEnd

    mini=ti_x_potmax(alpha,mu) *(1._kp+epsilon(1._kp)) !potential maximum
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_ti_x_rrad,mini,maxi,tolFind,tiData)
    ti_x_rrad = x

    if (present(bfold)) then
       bfold = -(ti_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function ti_x_rrad

  function find_ti_x_rrad(x,tiData)   
    implicit none
    real(kp) :: find_ti_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: tiData

    real(kp) :: primStar,alpha,mu,CalFplusPrimEnd,potStar,epsOneStar

    alpha=tiData%real1
    mu=tiData%real2
    CalFplusPrimEnd = tiData%real3

    primStar = ti_efold_primitive(x,alpha,mu)
    epsOneStar = ti_epsilon_one(x,alpha,mu)
    potStar = ti_norm_potential(x,alpha)

    find_ti_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_ti_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ti_x_rreh(alpha,mu,lnRreh,bfold)    
    implicit none
    real(kp) :: ti_x_rreh
    real(kp), intent(in) :: alpha,mu,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: tiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=ti_x_endinf(alpha,mu)
    
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)
    potEnd = ti_norm_potential(xEnd,alpha)
    primEnd = ti_efold_primitive(xEnd,alpha,mu)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    tiData%real1 = alpha
    tiData%real2 = mu
    tiData%real3 = calF + primEnd

    mini=ti_x_potmax(alpha,mu) *(1._kp+epsilon(1._kp)) !potential maximum
    maxi = xEnd*(1._kp-epsilon(1._kp))

    x = zbrent(find_ti_x_rreh,mini,maxi,tolFind,tiData)
    ti_x_rreh = x

    if (present(bfold)) then
       bfold = -(ti_efold_primitive(x,alpha,mu) - primEnd)
    endif

  end function ti_x_rreh

  function find_ti_x_rreh(x,tiData)   
    implicit none
    real(kp) :: find_ti_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: tiData

    real(kp) :: primStar,alpha,mu,CalFplusPrimEnd,potStar

    alpha=tiData%real1
    mu=tiData%real2
    CalFplusPrimEnd = tiData%real3

    primStar = ti_efold_primitive(x,alpha,mu)
    potStar = ti_norm_potential(x,alpha)

    find_ti_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_ti_x_rreh




  function ti_lnrhoreh_max(alpha,mu,Pstar) 
    implicit none
    real(kp) :: ti_lnrhoreh_max
    real(kp), intent(in) :: alpha,mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    xEnd=ti_x_endinf(alpha,mu) 
    potEnd  = ti_norm_potential(xEnd,alpha)
    epsOneEnd = ti_epsilon_one(xEnd,alpha,mu)

       
    x = ti_x_star(alpha,mu,wrad,junk,Pstar)    
    potStar = ti_norm_potential(x,alpha)
    epsOneStar = ti_epsilon_one(x,alpha,mu)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'ti_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ti_lnrhoreh_max = lnRhoEnd

  end function ti_lnrhoreh_max

  
  
end module tireheat

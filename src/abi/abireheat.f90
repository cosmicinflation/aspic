module abireheat
  use infprec, only : kp
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use abisr, only : abi_norm_potential, abi_epsilon_one
  use abisr, only : abi_x_trajectory, abi_efold_primitive

  implicit none

  private

  public abi_x_star, abi_lnrhoreh_max
  public abi_x_rrad, abi_x_rreh
  
contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function abi_x_star(alpha,beta,w,lnRhoReh,Pstar,bfold)
    implicit none
    real(kp) :: abi_x_star
    real(kp), intent(in) :: alpha,beta,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini, maxi, calF, x
    real(kp) :: primEnd, epsOneEnd, xEnd, potEnd

    type(transfert) :: abiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    xEnd = abi_x_endinf(alpha,beta)

!should be one    
    epsOneEnd = abi_epsilon_one(xEnd,alpha,beta)
    potEnd = abi_norm_potential(xEnd,alpha,beta)
!should be zero
    primEnd = abi_efold_primitive(xEnd,alpha,beta)

    print *,'zero one= ',primEnd,epsOneEnd

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    abiData%real1 = alpha
    abiData%real2 = beta
    abiData%real3 = w
    abiData%real4 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = abi_x_trajectory((-1._kp/epsilon(1._kp)),xend,alpha,beta)
       
    else
       mini = epsilon(1._kp)
       maxi = xEnd*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_abi_x_star,mini,maxi,tolFind,abiData)
    abi_x_star = x

    if (present(bfold)) then
       bfold = -(abi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function abi_x_star

  function find_abi_x_star(x,abiData)   
    implicit none
    real(kp) :: find_abi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: abiData

    real(kp) :: primStar,alpha,beta,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=abiData%real1
    beta=abiData%real2
    w = abiData%real3
    CalFplusPrimEnd = abiData%real4

    primStar = abi_efold_primitive(x,alpha,beta)
    epsOneStar = abi_epsilon_one(x,alpha,beta)
    potStar = abi_norm_potential(x,alpha,beta)

    find_abi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_abi_x_star



!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function abi_x_rrad(alpha,beta,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: abi_x_rrad
    real(kp), intent(in) :: alpha,beta,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: abiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=abi_x_endinf(alpha,beta)
    
    epsOneEnd = abi_epsilon_one(xEnd,alpha,beta)
    potEnd = abi_norm_potential(xEnd,alpha,beta)
    primEnd = abi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    abiData%real1 = alpha
    abiData%real2 = beta
    abiData%real3 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = abi_x_trajectory((-1._kp/epsilon(1._kp)),xend,alpha,beta)
       
    else
       mini = epsilon(1._kp)
       maxi = xEnd*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_abi_x_rrad,mini,maxi,tolFind,abiData)
    abi_x_rrad = x

    if (present(bfold)) then
       bfold = -(abi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function abi_x_rrad

  function find_abi_x_rrad(x,abiData)   
    implicit none
    real(kp) :: find_abi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: abiData

    real(kp) :: primStar,alpha,beta,CalFplusPrimEnd,potStar,epsOneStar

    alpha=abiData%real1
    beta=abiData%real2
    CalFplusPrimEnd = abiData%real3

    primStar = abi_efold_primitive(x,alpha,beta)
    epsOneStar = abi_epsilon_one(x,alpha,beta)
    potStar = abi_norm_potential(x,alpha,beta)

    find_abi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_abi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function abi_x_rreh(alpha,beta,lnRreh,bfold)    
    implicit none
    real(kp) :: abi_x_rreh
    real(kp), intent(in) :: alpha,beta,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xEnd,potEnd

    type(transfert) :: abiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd=abi_x_endinf(alpha,beta)
    
    epsOneEnd = abi_epsilon_one(xEnd,alpha,beta)
    potEnd = abi_norm_potential(xEnd,alpha,beta)
    primEnd = abi_efold_primitive(xEnd,alpha,beta)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    abiData%real1 = alpha
    abiData%real2 = beta
    abiData%real3 = calF + primEnd

    if (alpha.le.2._kp) then

       mini = xEnd*(1._kp+epsilon(1._kp))
       maxi = abi_x_trajectory((-1._kp/epsilon(1._kp)),xend,alpha,beta)
       
    else
       mini = epsilon(1._kp)
       maxi = xEnd*(1._kp-epsilon(1._kp))
    endif

    x = zbrent(find_abi_x_rreh,mini,maxi,tolFind,abiData)
    abi_x_rreh = x

    if (present(bfold)) then
       bfold = -(abi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function abi_x_rreh

  function find_abi_x_rreh(x,abiData)   
    implicit none
    real(kp) :: find_abi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: abiData

    real(kp) :: primStar,alpha,beta,CalFplusPrimEnd,potStar

    alpha=abiData%real1
    beta=abiData%real2
    CalFplusPrimEnd = abiData%real3

    primStar = abi_efold_primitive(x,alpha,beta)
    potStar = abi_norm_potential(x,alpha,beta)

    find_abi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_abi_x_rreh


  function abi_lnrhoreh_max(alpha,beta,Pstar) 
    implicit none
    real(kp) :: abi_lnrhoreh_max
    real(kp), intent(in) :: alpha,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    xEnd = abi_x_endinf(alpha,beta)
    potEnd  = abi_norm_potential(xEnd,alpha,beta)
    epsOneEnd = abi_epsilon_one(xEnd,alpha,beta)
       
    x = abi_x_star(alpha,beta,wrad,junk,Pstar)

    potStar = abi_norm_potential(x,alpha,beta)
    epsOneStar = abi_epsilon_one(x,alpha,beta)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'abi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    abi_lnrhoreh_max = lnRhoEnd

  end function abi_lnrhoreh_max



end module abireheat

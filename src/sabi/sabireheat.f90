!SuperConformal alpha-attractor B reheating functions in the slow-roll approximation

module sabireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sabisr, only : sabi_epsilon_one, sabi_epsilon_two, sabi_epsilon_three
  use sabisr, only : sabi_norm_potential, sabi_efold_primitive, sabi_x_endinf
  implicit none

  private

  public sabi_x_star, sabi_lnrhoreh_max
  public sabi_x_rrad, sabi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function sabi_x_star(alpha,n,xend,w,lnRhoReh,Pstar,bfold)
    implicit none
    real(kp) :: sabi_x_star
    real(kp), intent(in) :: alpha,n,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sabiData

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = sabi_epsilon_one(xEnd,alpha,n)
    potEnd = sabi_norm_potential(xEnd,alpha,n)
    primEnd = sabi_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    sabiData%real1 = alpha
    sabiData%real2 = n
    sabiData%real3 = w
    sabiData%real4 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(6._kp*alpha)/2._kp*log(4._kp*n*200._kp/(3._kp*alpha))

    x = zbrent(find_sabi_x_star,mini,maxi,tolFind,sabiData)
    sabi_x_star = x

    if (present(bfold)) then
       bfold = -(sabi_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sabi_x_star

  function find_sabi_x_star(x,sabiData)
    implicit none
    real(kp) :: find_sabi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sabiData

    real(kp) :: primStar,alpha,n,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=sabiData%real1
    n=sabiData%real2
    w = sabiData%real3
    CalFplusPrimEnd = sabiData%real4

    primStar = sabi_efold_primitive(x,alpha,n)
    epsOneStar = sabi_epsilon_one(x,alpha,n)
    potStar = sabi_norm_potential(x,alpha,n)

    find_sabi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_sabi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sabi_x_rrad(alpha,n,xend,lnRrad,Pstar,bfold)
    implicit none
    real(kp) :: sabi_x_rrad
    real(kp), intent(in) :: alpha,n,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sabiData

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = sabi_epsilon_one(xEnd,alpha,n)
    potEnd = sabi_norm_potential(xEnd,alpha,n)
    primEnd = sabi_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    sabiData%real1 = alpha
    sabiData%real2 = n
    sabiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(6._kp*alpha)/2._kp*log(4._kp*n*200._kp/(3._kp*alpha))

    x = zbrent(find_sabi_x_rrad,mini,maxi,tolFind,sabiData)
    sabi_x_rrad = x

    if (present(bfold)) then
       bfold = -(sabi_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sabi_x_rrad

  function find_sabi_x_rrad(x,sabiData)
    implicit none
    real(kp) :: find_sabi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sabiData

    real(kp) :: primStar,alpha,n,CalFplusPrimEnd,potStar,epsOneStar

    alpha=sabiData%real1
    n=sabiData%real2
    CalFplusPrimEnd = sabiData%real3

    primStar = sabi_efold_primitive(x,alpha,n)
    epsOneStar = sabi_epsilon_one(x,alpha,n)
    potStar = sabi_norm_potential(x,alpha,n)

    find_sabi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_sabi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sabi_x_rreh(alpha,n,xend,lnRreh,bfold)
    implicit none
    real(kp) :: sabi_x_rreh
    real(kp), intent(in) :: alpha,n,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: sabiData

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = sabi_epsilon_one(xEnd,alpha,n)
    potEnd = sabi_norm_potential(xEnd,alpha,n)
    primEnd = sabi_efold_primitive(xEnd,alpha,n)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    sabiData%real1 = alpha
    sabiData%real2 = n
    sabiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    !maxi = sqrt(6._kp*alpha)/2._kp*log(4._kp*n*200._kp/(3._kp*alpha))
    maxi = 100._kp*sqrt(alpha)

    x = zbrent(find_sabi_x_rreh,mini,maxi,tolFind,sabiData)
    sabi_x_rreh = x

    if (present(bfold)) then
       bfold = -(sabi_efold_primitive(x,alpha,n) - primEnd)
    endif

  end function sabi_x_rreh

  function find_sabi_x_rreh(x,sabiData)
    implicit none
    real(kp) :: find_sabi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sabiData

    real(kp) :: primStar,alpha,n,CalFplusPrimEnd,potStar

    alpha=sabiData%real1
    n=sabiData%real2
    CalFplusPrimEnd = sabiData%real3

    primStar = sabi_efold_primitive(x,alpha,n)
    potStar = sabi_norm_potential(x,alpha,n)

    find_sabi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_sabi_x_rreh


  function sabi_lnrhoreh_max(alpha,n,xend,Pstar)
    implicit none
    real(kp) :: sabi_lnrhoreh_max
    real(kp), intent(in) :: alpha,n,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd

    potEnd  = sabi_norm_potential(xEnd,alpha,n)
    epsOneEnd = sabi_epsilon_one(xEnd,alpha,n)

    x = sabi_x_star(alpha,n,xend,wrad,junk,Pstar)

    potStar = sabi_norm_potential(x,alpha,n)
    epsOneStar = sabi_epsilon_one(x,alpha,n)

    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar
        stop 'sabi_lnrhoreh_max: slow-roll violated!'
    endif

    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sabi_lnrhoreh_max = lnRhoEnd

  end function sabi_lnrhoreh_max



end module sabireheat

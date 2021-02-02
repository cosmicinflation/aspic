!hybrid natural reheating common functions in the slow-roll approximation

module hnicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hnicommon, only : hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  use hnicommon, only : hni_norm_potential, hni_efold_primitive
  implicit none

  private

  public hni_x_star
  public hni_x_rrad, hni_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hni_x_star(alpha,f,w,lnRhoReh,Pstar,xend,xmin,xmax,bfold)    
    implicit none
    real(kp) :: hni_x_star
    real(kp), intent(in) :: alpha,f,w,lnRhoReh,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,f)
    potEnd = hni_norm_potential(xEnd,alpha,f)
    primEnd = hni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = f
    hniData%real3 = w
    hniData%real4 = calF + primEnd

    x = zbrent(find_hni_x_star,xmin,xmax,tolFind,hniData)
    hni_x_star = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function hni_x_star

  function find_hni_x_star(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,f,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=hniData%real1
    f=hniData%real2
    w = hniData%real3
    CalFplusPrimEnd = hniData%real4

    primStar = hni_efold_primitive(x,alpha,f)
    epsOneStar = hni_epsilon_one(x,alpha,f)
    potStar = hni_norm_potential(x,alpha,f)

    find_hni_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_hni_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hni_x_rrad(alpha,f,lnRrad,Pstar,xend,xmin,xmax,bfold)    
    implicit none
    real(kp) :: hni_x_rrad
    real(kp), intent(in) :: alpha,f,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,f)
    potEnd = hni_norm_potential(xEnd,alpha,f)
    primEnd = hni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = f
    hniData%real3 = calF + primEnd

    x = zbrent(find_hni_x_rrad,xmin,xmax,tolFind,hniData)
    hni_x_rrad = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function hni_x_rrad

  function find_hni_x_rrad(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,f,CalFplusPrimEnd,potStar,epsOneStar

    alpha=hniData%real1
    f=hniData%real2
    CalFplusPrimEnd = hniData%real3

    primStar = hni_efold_primitive(x,alpha,f)
    epsOneStar = hni_epsilon_one(x,alpha,f)
    potStar = hni_norm_potential(x,alpha,f)

    find_hni_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_hni_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hni_x_rreh(alpha,f,lnRreh,xend,xmin,xmax,bfold)    
    implicit none
    real(kp) :: hni_x_rreh
    real(kp), intent(in) :: alpha,f,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,f)
    potEnd = hni_norm_potential(xEnd,alpha,f)
    primEnd = hni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = f
    hniData%real3 = calF + primEnd

    x = zbrent(find_hni_x_rreh,xmin,xmax,tolFind,hniData)
    hni_x_rreh = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function hni_x_rreh

  function find_hni_x_rreh(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,f,CalFplusPrimEnd,potStar

    alpha=hniData%real1
    f=hniData%real2
    CalFplusPrimEnd = hniData%real3

    primStar = hni_efold_primitive(x,alpha,f)
    potStar = hni_norm_potential(x,alpha,f)

    find_hni_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_hni_x_rreh


  
  
end module hnicomreh

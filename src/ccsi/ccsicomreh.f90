!Cubicly corrected Starobinski inflation common reheating functions in
!the slow-roll approximations

module ccsicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use ccsicommon, only : ccsi_epsilon_one
  use ccsicommon, only : ccsi_norm_potential, ccsi_efold_primitive
  
  implicit none

  private
  
  public ccsi_x_star, ccsi_x_rrad, ccsi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ccsi_x_star(alpha,w,lnRhoReh,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ccsi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd, lnOmega4End
    type(transfert) :: ccsiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ccsi_epsilon_one(xEnd,alpha)
    potEnd = ccsi_norm_potential(xEnd,alpha)

    primEnd = ccsi_efold_primitive(xEnd,alpha)

    lnOmega4End = 2._kp*xend
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd,lnOmega4End)

    ccsiData%real1 = alpha
    ccsiData%real2 = w
    ccsiData%real3 = calF + primEnd

    x = zbrent(find_ccsi_x_star,xmin,xmax,tolzbrent,ccsiData)
    ccsi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi_x_star

  function find_ccsi_x_star(x,ccsiData)   
    implicit none
    real(kp) :: find_ccsi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ccsiData%real1
    w = ccsiData%real2
    CalFplusprimEnd = ccsiData%real3

    primStar = ccsi_efold_primitive(x,alpha)
    epsOneStar = ccsi_epsilon_one(x,alpha)
    potStar = ccsi_norm_potential(x,alpha)

    find_ccsi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ccsi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ccsi_x_rrad(alpha,lnRrad,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ccsi_x_rrad
    real(kp), intent(in) :: alpha,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: ccsiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ccsi_epsilon_one(xEnd,alpha)
    potEnd = ccsi_norm_potential(xEnd,alpha)

    primEnd = ccsi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    ccsiData%real1 = alpha
    ccsiData%real2 = calF + primEnd

    x = zbrent(find_ccsi_x_rrad,xmin,xmax,tolzbrent,ccsiData)
    ccsi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi_x_rrad

  function find_ccsi_x_rrad(x,ccsiData)   
    implicit none
    real(kp) :: find_ccsi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=ccsiData%real1
    CalFplusprimEnd = ccsiData%real2

    primStar = ccsi_efold_primitive(x,alpha)
    epsOneStar = ccsi_epsilon_one(x,alpha)
    potStar = ccsi_norm_potential(x,alpha)

    find_ccsi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_ccsi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ccsi_x_rreh(alpha,lnRreh,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ccsi_x_rreh
    real(kp), intent(in) :: alpha,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd,lnOmega4End
    type(transfert) :: ccsiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ccsi_epsilon_one(xEnd,alpha)
    potEnd = ccsi_norm_potential(xEnd,alpha)

    primEnd = ccsi_efold_primitive(xEnd,alpha)

    lnOmega4End = 2._kp*xend
    
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd,lnOmega4End)

    ccsiData%real1 = alpha
    ccsiData%real2 = calF + primEnd

    x = zbrent(find_ccsi_x_rreh,xmin,xmax,tolzbrent,ccsiData)
    ccsi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ccsi_efold_primitive(x,alpha) - primEnd)
    endif

  end function ccsi_x_rreh

  function find_ccsi_x_rreh(x,ccsiData)   
    implicit none
    real(kp) :: find_ccsi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ccsiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=ccsiData%real1
    CalFplusprimEnd = ccsiData%real2

    primStar = ccsi_efold_primitive(x,alpha)
    potStar = ccsi_norm_potential(x,alpha)

    find_ccsi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_ccsi_x_rreh


  
end module ccsicomreh

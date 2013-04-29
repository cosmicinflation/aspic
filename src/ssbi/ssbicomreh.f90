!spontaneous symmetry breaking common reheating functions in the
!slow-roll approximations

module ssbicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use ssbicommon, only : ssbi_epsilon_one
  use ssbicommon, only : ssbi_norm_potential, ssbi_efold_primitive

  
  implicit none

  private

  public ssbi_x_star, ssbi_x_rrad, ssbi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ssbi_x_star(alpha,beta,w,lnRhoReh,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ssbi_x_star
    real(kp), intent(in) :: alpha,beta,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: ssbiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ssbi_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ssbiData%real1 = alpha
    ssbiData%real2 = beta
    ssbiData%real3 = w
    ssbiData%real4 = calF + primEnd

    x = zbrent(find_ssbi_x_star,xmin,xmax,tolzbrent,ssbiData)
    ssbi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi_x_star

  function find_ssbi_x_star(x,ssbiData)   
    implicit none
    real(kp) :: find_ssbi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbiData

    real(kp) :: primStar,alpha,beta,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbiData%real1
    beta=ssbiData%real2
    w = ssbiData%real3
    CalFplusprimEnd = ssbiData%real4

    primStar = ssbi_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi_epsilon_one(x,alpha,beta)
    potStar = ssbi_norm_potential(x,alpha,beta)

    find_ssbi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ssbi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ssbi_x_rrad(alpha,beta,lnRrad,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ssbi_x_rrad
    real(kp), intent(in) :: alpha,beta,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: ssbiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ssbi_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    ssbiData%real1 = alpha
    ssbiData%real2 = beta
    ssbiData%real3 = calF + primEnd

    x = zbrent(find_ssbi_x_rrad,xmin,xmax,tolzbrent,ssbiData)
    ssbi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi_x_rrad

  function find_ssbi_x_rrad(x,ssbiData)   
    implicit none
    real(kp) :: find_ssbi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbiData

    real(kp) :: primStar,alpha,beta,CalFplusprimEnd,potStar,epsOneStar

    alpha=ssbiData%real1
    beta=ssbiData%real2
    CalFplusprimEnd = ssbiData%real3

    primStar = ssbi_efold_primitive(x,alpha,beta)
    epsOneStar = ssbi_epsilon_one(x,alpha,beta)
    potStar = ssbi_norm_potential(x,alpha,beta)

    find_ssbi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_ssbi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function ssbi_x_rreh(alpha,beta,lnRreh,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: ssbi_x_rreh
    real(kp), intent(in) :: alpha,beta,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: ssbiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = ssbi_epsilon_one(xEnd,alpha,beta)
    potEnd = ssbi_norm_potential(xEnd,alpha,beta)

    primEnd = ssbi_efold_primitive(xEnd,alpha,beta)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    ssbiData%real1 = alpha
    ssbiData%real2 = beta
    ssbiData%real3 = calF + primEnd

    x = zbrent(find_ssbi_x_rreh,xmin,xmax,tolzbrent,ssbiData)
    ssbi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ssbi_efold_primitive(x,alpha,beta) - primEnd)
    endif

  end function ssbi_x_rreh

  function find_ssbi_x_rreh(x,ssbiData)   
    implicit none
    real(kp) :: find_ssbi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ssbiData

    real(kp) :: primStar,alpha,beta,CalFplusprimEnd,potStar

    alpha=ssbiData%real1
    beta=ssbiData%real2
    CalFplusprimEnd = ssbiData%real3

    primStar = ssbi_efold_primitive(x,alpha,beta)
    potStar = ssbi_norm_potential(x,alpha,beta)

    find_ssbi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_ssbi_x_rreh


  
end module ssbicomreh

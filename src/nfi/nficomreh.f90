!N-formalism inflation common reheating functions in the slow-roll
!approximations

module nficomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use nficommon, only : nfi_epsilon_one
  use nficommon, only : nfi_norm_potential, nfi_efold_primitive

  
  implicit none

  private
  
  public nfi_x_star, nfi_x_rrad, nfi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function nfi_x_star(a,b,w,lnRhoReh,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: nfi_x_star
    real(kp), intent(in) :: a,b,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: nfiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nfi_epsilon_one(xEnd,a,b)
    potEnd = nfi_norm_potential(xEnd,a,b)

    primEnd = nfi_efold_primitive(xEnd,a,b)
    
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    nfiData%real1 = a
    nfiData%real2 = b
    nfiData%real3 = w
    nfiData%real4 = calF + primEnd

    x = zbrent(find_nfi_x_star,xmin,xmax,tolzbrent,nfiData)
    nfi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (nfi_efold_primitive(x,a,b) - primEnd)
    endif

  end function nfi_x_star

  function find_nfi_x_star(x,nfiData)   
    implicit none
    real(kp) :: find_nfi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: nfiData

    real(kp) :: primStar,a,b,w,CalFplusprimEnd,potStar,epsOneStar

    a=nfiData%real1
    b=nfiData%real2
    w = nfiData%real3
    CalFplusprimEnd = nfiData%real4

    primStar = nfi_efold_primitive(x,a,b)
    epsOneStar = nfi_epsilon_one(x,a,b)
    potStar = nfi_norm_potential(x,a,b)

    find_nfi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_nfi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function nfi_x_rrad(a,b,lnRrad,Pstar,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: nfi_x_rrad
    real(kp), intent(in) :: a,b,lnRrad,Pstar
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: nfiData

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nfi_epsilon_one(xEnd,a,b)
    potEnd = nfi_norm_potential(xEnd,a,b)

    primEnd = nfi_efold_primitive(xEnd,a,b)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    nfiData%real1 = a
    nfiData%real2 = b
    nfiData%real3 = calF + primEnd

    x = zbrent(find_nfi_x_rrad,xmin,xmax,tolzbrent,nfiData)
    nfi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (nfi_efold_primitive(x,a,b) - primEnd)
    endif

  end function nfi_x_rrad

  function find_nfi_x_rrad(x,nfiData)   
    implicit none
    real(kp) :: find_nfi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: nfiData

    real(kp) :: primStar,a,b,CalFplusprimEnd,potStar,epsOneStar

    a=nfiData%real1
    b=nfiData%real2
    CalFplusprimEnd = nfiData%real3

    primStar = nfi_efold_primitive(x,a,b)
    epsOneStar = nfi_epsilon_one(x,a,b)
    potStar = nfi_norm_potential(x,a,b)

    find_nfi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_nfi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function nfi_x_rreh(a,b,lnRreh,xend,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: nfi_x_rreh
    real(kp), intent(in) :: a,b,lnRreh
    real(kp), intent(in) :: xend, xmin, xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: nfiData

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = nfi_epsilon_one(xEnd,a,b)
    potEnd = nfi_norm_potential(xEnd,a,b)

    primEnd = nfi_efold_primitive(xEnd,a,b)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    nfiData%real1 = a
    nfiData%real2 = b
    nfiData%real3 = calF + primEnd

    x = zbrent(find_nfi_x_rreh,xmin,xmax,tolzbrent,nfiData)
    nfi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (nfi_efold_primitive(x,a,b) - primEnd)
    endif

  end function nfi_x_rreh

  function find_nfi_x_rreh(x,nfiData)   
    implicit none
    real(kp) :: find_nfi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: nfiData

    real(kp) :: primStar,a,b,CalFplusprimEnd,potStar

    a=nfiData%real1
    b=nfiData%real2
    CalFplusprimEnd = nfiData%real3

    primStar = nfi_efold_primitive(x,a,b)
    potStar = nfi_norm_potential(x,a,b)

    find_nfi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_nfi_x_rreh


  
end module nficomreh

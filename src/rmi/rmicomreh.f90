!common function for rmi1234 reheat modules

module rmicomreh
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rmicommon, only : rmi_epsilon_one, rmi_epsilon_two, rmi_epsilon_three
  use rmicommon, only : rmi_norm_potential, rmi_efold_primitive

  implicit none

  private

  public rmi_x_star, rmi_x_rrad, rmi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi_x_star(c,phi0,xend,w,lnRhoReh,Pstar,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: rmi_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(in) :: xmin,xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmiData

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi_epsilon_one(xEnd,c,phi0)
    potEnd = rmi_norm_potential(xEnd,c,phi0)

    primEnd = rmi_efold_primitive(xEnd,c,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rmiData%real1 = c
    rmiData%real2 = phi0
    rmiData%real3 = w
    rmiData%real4 = calF + primEnd

    x = zbrent(find_rmi_x_star,xmin,xmax,tolzbrent,rmiData)
    rmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi_x_star

  function find_rmi_x_star(x,rmiData)   
    implicit none
    real(kp) :: find_rmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmiData

    real(kp) :: primStar,c,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    c=rmiData%real1
    phi0=rmiData%real2
    w = rmiData%real3
    CalFplusprimEnd = rmiData%real4

    primStar = rmi_efold_primitive(x,c,phi0)
    epsOneStar = rmi_epsilon_one(x,c,phi0)
    potStar = rmi_norm_potential(x,c,phi0)

    find_rmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rmi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rmi_x_rrad(c,phi0,xend,lnRrad,Pstar,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: rmi_x_rrad
    real(kp), intent(in) :: c,phi0,xend,lnRrad,Pstar
    real(kp), intent(in) :: xmin,xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmiData


    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi_epsilon_one(xEnd,c,phi0)
    potEnd = rmi_norm_potential(xEnd,c,phi0)

    primEnd = rmi_efold_primitive(xEnd,c,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rmiData%real1 = c
    rmiData%real2 = phi0
    rmiData%real3 = calF + primEnd

    x = zbrent(find_rmi_x_rrad,xmin,xmax,tolzbrent,rmiData)
    rmi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi_x_rrad

  function find_rmi_x_rrad(x,rmiData)   
    implicit none
    real(kp) :: find_rmi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmiData

    real(kp) :: primStar,c,phi0,CalFplusprimEnd,potStar,epsOneStar

    c=rmiData%real1
    phi0=rmiData%real2
    CalFplusprimEnd = rmiData%real3

    primStar = rmi_efold_primitive(x,c,phi0)
    epsOneStar = rmi_epsilon_one(x,c,phi0)
    potStar = rmi_norm_potential(x,c,phi0)

    find_rmi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_rmi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rmi_x_rreh(c,phi0,xend,lnRreh,xmin,xmax,bfoldstar)    
    implicit none
    real(kp) :: rmi_x_rreh
    real(kp), intent(in) :: c,phi0,xend,lnRreh
    real(kp), intent(in) :: xmin,xmax
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmiData


    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi_epsilon_one(xEnd,c,phi0)
    potEnd = rmi_norm_potential(xEnd,c,phi0)

    primEnd = rmi_efold_primitive(xEnd,c,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rmiData%real1 = c
    rmiData%real2 = phi0
    rmiData%real3 = calF + primEnd

    x = zbrent(find_rmi_x_rreh,xmin,xmax,tolzbrent,rmiData)
    rmi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi_x_rreh

  function find_rmi_x_rreh(x,rmiData)   
    implicit none
    real(kp) :: find_rmi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmiData

    real(kp) :: primStar,c,phi0,CalFplusprimEnd,potStar

    c=rmiData%real1
    phi0=rmiData%real2
    CalFplusprimEnd = rmiData%real3

    primStar = rmi_efold_primitive(x,c,phi0)   
    potStar = rmi_norm_potential(x,c,phi0)

    find_rmi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_rmi_x_rreh


end module rmicomreh

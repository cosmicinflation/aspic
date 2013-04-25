!running mass 3 reheating functions in the slow-roll approximations

module rmi3reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rmi3sr, only : rmi3_epsilon_one, rmi3_epsilon_two, rmi3_epsilon_three
  use rmi3sr, only : rmi3_norm_potential, rmi3_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public rmi3_x_star, rmi3_lnrhoreh_max

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi3_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi3_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmi3Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi3_epsilon_one(xEnd,c,phi0)
    potEnd = rmi3_norm_potential(xEnd,c,phi0)

    primEnd = rmi3_efold_primitive(xEnd,c,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rmi3Data%real1 = c
    rmi3Data%real2 = phi0
    rmi3Data%real3 = w
    rmi3Data%real4 = calF + primEnd

    mini = epsilon(1._kp)
    maxi = xend*(1._kp-epsilon(1._kp))

    x = zbrent(find_rmi3_x_star,mini,maxi,tolzbrent,rmi3Data)
    rmi3_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi3_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi3_x_star

  function find_rmi3_x_star(x,rmi3Data)   
    implicit none
    real(kp) :: find_rmi3_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmi3Data

    real(kp) :: primStar,c,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    c=rmi3Data%real1
    phi0=rmi3Data%real2
    w = rmi3Data%real3
    CalFplusprimEnd = rmi3Data%real4

    primStar = rmi3_efold_primitive(x,c,phi0)
    epsOneStar = rmi3_epsilon_one(x,c,phi0)
    potStar = rmi3_norm_potential(x,c,phi0)

    find_rmi3_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rmi3_x_star



  function rmi3_lnrhoreh_max(c,phi0,xend,Pstar) 
    implicit none
    real(kp) :: rmi3_lnrhoreh_max
    real(kp), intent(in) :: c, phi0, xend, Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd

    potEnd  = rmi3_norm_potential(xEnd,c,phi0)

    epsOneEnd = rmi3_epsilon_one(xEnd,c,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = rmi3_x_star(c,phi0,xend,wrad,junk,Pstar)  


    potStar = rmi3_norm_potential(x,c,phi0)
    epsOneStar = rmi3_epsilon_one(x,c,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rmi3_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rmi3_lnrhoreh_max = lnRhoEnd

  end function rmi3_lnrhoreh_max

  
end module rmi3reheat

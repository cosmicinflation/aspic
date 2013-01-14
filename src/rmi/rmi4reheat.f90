!running mass 4 reheating functions in the slow-roll approximations

module rmi4reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rmi4sr, only : rmi4_epsilon_one, rmi4_epsilon_two, rmi4_epsilon_three
  use rmi4sr, only : rmi4_norm_potential, rmi4_efold_primitive
  use cosmopar, only : QrmsOverT

  implicit none

  private

  public rmi4_x_star, rmi4_lnrhoend

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rmi4_x_star(c,phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rmi4_x_star
    real(kp), intent(in) :: c,phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    type(transfert) :: rmi4Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = rmi4_epsilon_one(xEnd,c,phi0)
    potEnd = rmi4_norm_potential(xEnd,c,phi0)

    primEnd = rmi4_efold_primitive(xEnd,c,phi0)


   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rmi4Data%real1 = c
    rmi4Data%real2 = phi0
    rmi4Data%real3 = w
    rmi4Data%real4 = calF + primEnd


    mini = xend*(1._kp+epsilon(1._kp))
    maxi = 10._kp

    x = zbrent(find_rmi4_x_star,mini,maxi,tolzbrent,rmi4Data)
    rmi4_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rmi4_efold_primitive(x,c,phi0) - primEnd)
    endif

  end function rmi4_x_star

  function find_rmi4_x_star(x,rmi4Data)   
    implicit none
    real(kp) :: find_rmi4_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rmi4Data

    real(kp) :: primStar,c,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    c=rmi4Data%real1
    phi0=rmi4Data%real2
    w = rmi4Data%real3
    CalFplusprimEnd = rmi4Data%real4

    primStar = rmi4_efold_primitive(x,c,phi0)
    epsOneStar = rmi4_epsilon_one(x,c,phi0)
    potStar = rmi4_norm_potential(x,c,phi0)

    find_rmi4_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rmi4_x_star



  function rmi4_lnrhoend(c,phi0,xend,Pstar) 
    implicit none
    real(kp) :: rmi4_lnrhoend
    real(kp), intent(in) :: c, phi0, xend, Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    

    potEnd  = rmi4_norm_potential(xEnd,c,phi0)

    epsOneEnd = rmi4_epsilon_one(xEnd,c,phi0)


!   Trick to return x such that rho_reh=rho_end

    x = rmi4_x_star(c,phi0,xend,wrad,junk,Pstar)  


    potStar = rmi4_norm_potential(x,c,phi0)
    epsOneStar = rmi4_epsilon_one(x,c,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rmi4_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rmi4_lnrhoend = lnRhoEnd

  end function rmi4_lnrhoend

  
end module rmi4reheat

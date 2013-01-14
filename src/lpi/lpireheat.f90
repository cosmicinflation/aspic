!Logarithmic Potential inflation reheating functions in the
!slow-roll approximations

module lpireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lpisr, only : lpi_epsilon_one, lpi_epsilon_two, lpi_epsilon_three
  use lpisr, only : lpi_norm_potential, lpi_x_endinf, lpi_efold_primitive
  implicit none

  private

  public lpi_x_star, lpi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lpi_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lpi_x_star
    real(kp), intent(in) :: phi0,p,q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: lpiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lpi_x_endinf(p,q,phi0)
    epsOneEnd = lpi_epsilon_one(xEnd,p,q,phi0)
    potEnd = lpi_norm_potential(xEnd,p,q,phi0)
    primEnd = lpi_efold_primitive(xEnd,p,q,phi0) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    lpiData%real1 = p
    lpiData%real2 = q
    lpiData%real3 = phi0
    lpiData%real4 = xEnd
    lpiData%real5 = w
    lpiData%real6 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = sqrt(xend**2+2._kp*p/(phi0**2)*10._kp**(4.)) !Uses an approximate solution for the slow roll trajectory when x>>1 and sets maxi at ~ 10^4 efolds away from xend

    x = zbrent(find_lpi_x_star,mini,maxi,tolzbrent,lpiData)
    lpi_x_star = x  


    if (present(bfoldstar)) then
       bfoldstar = - (lpi_efold_primitive(x,p,q,phi0) - primEnd)
    endif

  end function lpi_x_star

  function find_lpi_x_star(x,lpiData)   
    implicit none
    real(kp) :: find_lpi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lpiData

    real(kp) :: primStar,phi0,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar


    p=lpiData%real1
    q=lpiData%real2
    phi0=lpiData%real3
    xEnd=lpiData%real4
    w = lpiData%real5
    CalFplusprimEnd = lpiData%real6

    primStar = lpi_efold_primitive(x,p,q,phi0)
    epsOneStar = lpi_epsilon_one(x,p,q,phi0)
    potStar = lpi_norm_potential(x,p,q,phi0)

    find_lpi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lpi_x_star



  function lpi_lnrhoend(p,q,phi0,Pstar) 
    implicit none
    real(kp) :: lpi_lnrhoend
    real(kp), intent(in) :: phi0,p,q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lpi_x_endinf(p,q,phi0)
    potEnd  = lpi_norm_potential(xEnd,p,q,phi0)
    epsOneEnd = lpi_epsilon_one(xEnd,p,q,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = lpi_x_star(p,q,phi0,wrad,junk,Pstar)    
    potStar = lpi_norm_potential(x,p,q,phi0)
    epsOneStar = lpi_epsilon_one(x,p,q,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lpi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lpi_lnrhoend = lnRhoEnd

  end function lpi_lnrhoend

  
end module lpireheat

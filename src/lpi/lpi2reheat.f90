!Logarithmic Potential inflation reheating functions in the
!slow-roll approximations

module lpi2reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lpicommon, only : lpi_x_potmax
  use lpi2sr, only : lpi2_epsilon_one, lpi2_epsilon_two, lpi2_epsilon_three
  use lpi2sr, only : lpi2_norm_potential, lpi2_x_endinf, lpi2_efold_primitive
  implicit none

  private

  public lpi2_x_star, lpi2_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lpi2_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lpi2_x_star
    real(kp), intent(in) :: phi0,p,q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: lpi2Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lpi2_x_endinf(p,q,phi0)
    epsOneEnd = lpi2_epsilon_one(xEnd,p,q,phi0)
    potEnd = lpi2_norm_potential(xEnd,p,q,phi0)
    primEnd = lpi2_efold_primitive(xEnd,p,q,phi0) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    lpi2Data%real1 = p
    lpi2Data%real2 = q
    lpi2Data%real3 = phi0
    lpi2Data%real4 = xEnd
    lpi2Data%real5 = w
    lpi2Data%real6 = calF + primEnd

    mini = lpi_x_potmax(p,q)
    maxi = xend

    x = zbrent(find_lpi2_x_star,mini,maxi,tolzbrent,lpi2Data)
    lpi2_x_star = x  


    if (present(bfoldstar)) then
       bfoldstar = - (lpi2_efold_primitive(x,p,q,phi0) - primEnd)
    endif

  end function lpi2_x_star

  function find_lpi2_x_star(x,lpi2Data)   
    implicit none
    real(kp) :: find_lpi2_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lpi2Data

    real(kp) :: primStar,phi0,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar


    p=lpi2Data%real1
    q=lpi2Data%real2
    phi0=lpi2Data%real3
    xEnd=lpi2Data%real4
    w = lpi2Data%real5
    CalFplusprimEnd = lpi2Data%real6

    primStar = lpi2_efold_primitive(x,p,q,phi0)
    epsOneStar = lpi2_epsilon_one(x,p,q,phi0)
    potStar = lpi2_norm_potential(x,p,q,phi0)

    find_lpi2_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lpi2_x_star



  function lpi2_lnrhoreh_max(p,q,phi0,Pstar) 
    implicit none
    real(kp) :: lpi2_lnrhoreh_max
    real(kp), intent(in) :: phi0,p,q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lpi2_x_endinf(p,q,phi0)
    potEnd  = lpi2_norm_potential(xEnd,p,q,phi0)
    epsOneEnd = lpi2_epsilon_one(xEnd,p,q,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = lpi2_x_star(p,q,phi0,wrad,junk,Pstar)    
    potStar = lpi2_norm_potential(x,p,q,phi0)
    epsOneStar = lpi2_epsilon_one(x,p,q,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lpi2_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lpi2_lnrhoreh_max = lnRhoEnd

  end function lpi2_lnrhoreh_max

  
end module lpi2reheat

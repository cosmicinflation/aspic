!Logarithmic Potential inflation reheating functions in the
!slow-roll approximations

module lpi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lpi1sr, only : lpi1_epsilon_one, lpi1_epsilon_two, lpi1_epsilon_three
  use lpi1sr, only : lpi1_norm_potential, lpi1_x_endinf, lpi1_efold_primitive
  implicit none

  private

  public lpi1_x_star, lpi1_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lpi1_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lpi1_x_star
    real(kp), intent(in) :: phi0,p,q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: lpi1Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lpi1_x_endinf(p,q,phi0)
    epsOneEnd = lpi1_epsilon_one(xEnd,p,q,phi0)
    potEnd = lpi1_norm_potential(xEnd,p,q,phi0)
    primEnd = lpi1_efold_primitive(xEnd,p,q,phi0) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    lpi1Data%real1 = p
    lpi1Data%real2 = q
    lpi1Data%real3 = phi0
    lpi1Data%real4 = xEnd
    lpi1Data%real5 = w
    lpi1Data%real6 = calF + primEnd

    mini = xend
    maxi = 1._kp/epsilon(1._kp)

    x = zbrent(find_lpi1_x_star,mini,maxi,tolzbrent,lpi1Data)
    lpi1_x_star = x  


    if (present(bfoldstar)) then
       bfoldstar = - (lpi1_efold_primitive(x,p,q,phi0) - primEnd)
    endif

  end function lpi1_x_star

  function find_lpi1_x_star(x,lpi1Data)   
    implicit none
    real(kp) :: find_lpi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lpi1Data

    real(kp) :: primStar,phi0,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar


    p=lpi1Data%real1
    q=lpi1Data%real2
    phi0=lpi1Data%real3
    xEnd=lpi1Data%real4
    w = lpi1Data%real5
    CalFplusprimEnd = lpi1Data%real6

    primStar = lpi1_efold_primitive(x,p,q,phi0)
    epsOneStar = lpi1_epsilon_one(x,p,q,phi0)
    potStar = lpi1_norm_potential(x,p,q,phi0)

    find_lpi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lpi1_x_star



  function lpi1_lnrhoreh_max(p,q,phi0,Pstar) 
    implicit none
    real(kp) :: lpi1_lnrhoreh_max
    real(kp), intent(in) :: phi0,p,q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lpi1_x_endinf(p,q,phi0)
    potEnd  = lpi1_norm_potential(xEnd,p,q,phi0)
    epsOneEnd = lpi1_epsilon_one(xEnd,p,q,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = lpi1_x_star(p,q,phi0,wrad,junk,Pstar)    
    potStar = lpi1_norm_potential(x,p,q,phi0)
    epsOneStar = lpi1_epsilon_one(x,p,q,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lpi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lpi1_lnrhoreh_max = lnRhoEnd

  end function lpi1_lnrhoreh_max

  
end module lpi1reheat

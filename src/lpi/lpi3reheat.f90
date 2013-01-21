!Logarithmic Potential inflation reheating functions in the
!slow-roll approximations

module lpi3reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use lpicommon, only : lpi_x_potmax
  use lpi3sr, only : lpi3_epsilon_one, lpi3_epsilon_two, lpi3_epsilon_three
  use lpi3sr, only : lpi3_norm_potential, lpi3_x_endinf, lpi3_efold_primitive
  implicit none

  private

  public lpi3_x_star, lpi3_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lpi3_x_star(p,q,phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lpi3_x_star
    real(kp), intent(in) :: phi0,p,q,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: lpi3Data
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lpi3_x_endinf(p,q,phi0)
    epsOneEnd = lpi3_epsilon_one(xEnd,p,q,phi0)
    potEnd = lpi3_norm_potential(xEnd,p,q,phi0)
    primEnd = lpi3_efold_primitive(xEnd,p,q,phi0) 
    

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    lpi3Data%real1 = p
    lpi3Data%real2 = q
    lpi3Data%real3 = phi0
    lpi3Data%real4 = xEnd
    lpi3Data%real5 = w
    lpi3Data%real6 = calF + primEnd

    mini = xend
    maxi = lpi_x_potmax(p,q)

    x = zbrent(find_lpi3_x_star,mini,maxi,tolzbrent,lpi3Data)
    lpi3_x_star = x  


    if (present(bfoldstar)) then
       bfoldstar = - (lpi3_efold_primitive(x,p,q,phi0) - primEnd)
    endif

  end function lpi3_x_star

  function find_lpi3_x_star(x,lpi3Data)   
    implicit none
    real(kp) :: find_lpi3_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lpi3Data

    real(kp) :: primStar,phi0,p,q,xEnd,w,CalFplusprimEnd,potStar,epsOneStar


    p=lpi3Data%real1
    q=lpi3Data%real2
    phi0=lpi3Data%real3
    xEnd=lpi3Data%real4
    w = lpi3Data%real5
    CalFplusprimEnd = lpi3Data%real6

    primStar = lpi3_efold_primitive(x,p,q,phi0)
    epsOneStar = lpi3_epsilon_one(x,p,q,phi0)
    potStar = lpi3_norm_potential(x,p,q,phi0)

    find_lpi3_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lpi3_x_star



  function lpi3_lnrhoend(p,q,phi0,Pstar) 
    implicit none
    real(kp) :: lpi3_lnrhoend
    real(kp), intent(in) :: phi0,p,q,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = lpi3_x_endinf(p,q,phi0)
    potEnd  = lpi3_norm_potential(xEnd,p,q,phi0)
    epsOneEnd = lpi3_epsilon_one(xEnd,p,q,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = lpi3_x_star(p,q,phi0,wrad,junk,Pstar)    
    potStar = lpi3_norm_potential(x,p,q,phi0)
    epsOneStar = lpi3_epsilon_one(x,p,q,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'lpi3_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lpi3_lnrhoend = lnRhoEnd

  end function lpi3_lnrhoend

  
end module lpi3reheat

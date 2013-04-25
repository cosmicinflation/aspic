!double well inflation reheating functions in the slow-roll approximations

module wrhireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use wrhisr, only : wrhi_epsilon_one, wrhi_epsilon_two, wrhi_epsilon_three
  use wrhisr, only : wrhi_norm_potential, wrhi_x_trajectory
  use wrhisr, only : wrhi_x_endinf, wrhi_efold_primitive
  implicit none

  private

  public wrhi_x_star, wrhi_lnrhoreh_max 

contains

!returns x =phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function wrhi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: wrhi_x_star
    real(kp), intent(in) :: phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: wrhiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = wrhi_x_endinf(phi0)

    epsOneEnd = wrhi_epsilon_one(xEnd,phi0)
    potEnd = wrhi_norm_potential(xEnd)

    primEnd = wrhi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    wrhiData%real1 = phi0
    wrhiData%real2 = w
    wrhiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi =  wrhi_x_trajectory(-200._kp,xend,phi0) !Position 200 efolds before the end of inflation 

    x = zbrent(find_wrhi_x_star,mini,maxi,tolzbrent,wrhiData)
    wrhi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (wrhi_efold_primitive(x,phi0) - primEnd)
    endif

  end function wrhi_x_star


  function find_wrhi_x_star(x,wrhiData)   
    implicit none
    real(kp) :: find_wrhi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: wrhiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=wrhiData%real1
    w = wrhiData%real2
    CalFplusprimEnd = wrhiData%real3

    primStar = wrhi_efold_primitive(x,phi0)
    epsOneStar = wrhi_epsilon_one(x,phi0)
    potStar = wrhi_norm_potential(x)

    find_wrhi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_wrhi_x_star


  function wrhi_lnrhoreh_max(phi0,Pstar) 
    implicit none
    real(kp) :: wrhi_lnrhoreh_max
    real(kp), intent(in) :: phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = wrhi_x_endinf(phi0)
    potEnd  = wrhi_norm_potential(xEnd)
    epsOneEnd = wrhi_epsilon_one(xEnd,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = wrhi_x_star(phi0,wrad,junk,Pstar)  
 
    potStar = wrhi_norm_potential(x)
    epsOneStar = wrhi_epsilon_one(x,phi0)
   
    if (.not.slowroll_validity(epsOneStar)) stop 'wrhi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    wrhi_lnrhoreh_max = lnRhoEnd

  end function wrhi_lnrhoreh_max
  
end module wrhireheat

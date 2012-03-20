!Colemann Weinberg inflation reheating functions in the slow-roll approximations

module twireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use twisr, only : twi_epsilon_one, twi_epsilon_two, twi_epsilon_three
  use twisr, only : twi_norm_potential
  use twisr, only : twi_x_endinf, twi_efold_primitive
  implicit none

  private

  public twi_x_star, twi_lnrhoend 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function twi_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: twi_x_star
    real(kp), intent(in) :: phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: twiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = twi_x_endinf(phi0)

    epsOneEnd = twi_epsilon_one(xEnd,phi0)
    potEnd = twi_norm_potential(xEnd,phi0)

    primEnd = twi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    twiData%real1 = phi0
    twiData%real2 = w
    twiData%real3 = calF + primEnd

    mini = twi_x_endinf(phi0)
    maxi=100._kp*phi0 !if bigger, <numaccuracy errors

    x = zbrent(find_twi_x_star,mini,maxi,tolzbrent,twiData)
    twi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (twi_efold_primitive(x,phi0) - primEnd)
    endif

  end function twi_x_star

  function find_twi_x_star(x,twiData)   
    implicit none
    real(kp) :: find_twi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: twiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=twiData%real1
    w = twiData%real2
    CalFplusprimEnd = twiData%real3

    primStar = twi_efold_primitive(x,phi0)
    epsOneStar = twi_epsilon_one(x,phi0)
    potStar = twi_norm_potential(x,phi0)

    find_twi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_twi_x_star



  function twi_lnrhoend(phi0,Pstar) 
    implicit none
    real(kp) :: twi_lnrhoend
    real(kp), intent(in) :: phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = twi_x_endinf(phi0)


    potEnd  = twi_norm_potential(xEnd,phi0)

    epsOneEnd = twi_epsilon_one(xEnd,phi0)



!   Trick to return x such that rho_reh=rho_end

    x = twi_x_star(phi0,wrad,junk,Pstar)  

 
    potStar = twi_norm_potential(x,phi0)
    epsOneStar = twi_epsilon_one(x,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'twi_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    twi_lnrhoend = lnRhoEnd

  end function twi_lnrhoend

  
end module twireheat

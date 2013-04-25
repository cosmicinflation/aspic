!double well inflation reheating functions in the slow-roll approximations

module ostireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use ostisr, only : osti_epsilon_one, osti_epsilon_two, osti_epsilon_three
  use ostisr, only : osti_norm_potential
  use ostisr, only : osti_x_endinf, osti_efold_primitive
  implicit none

  private

  public osti_x_star, osti_lnrhoreh_max 

contains

!returns x =phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function osti_x_star(phi0,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: osti_x_star
    real(kp), intent(in) :: phi0,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: ostiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = osti_x_endinf(phi0)

    epsOneEnd = osti_epsilon_one(xEnd,phi0)
    potEnd = osti_norm_potential(xEnd)

    primEnd = osti_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    ostiData%real1 = phi0
    ostiData%real2 = w
    ostiData%real3 = calF + primEnd

    mini = xend*(1._kp+epsilon(1._kp))
    maxi = exp(-0.5_kp)*(1._kp-epsilon(1._kp))

    x = zbrent(find_osti_x_star,mini,maxi,tolzbrent,ostiData)
    osti_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (osti_efold_primitive(x,phi0) - primEnd)
    endif

  end function osti_x_star


  function find_osti_x_star(x,ostiData)   
    implicit none
    real(kp) :: find_osti_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: ostiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0=ostiData%real1
    w = ostiData%real2
    CalFplusprimEnd = ostiData%real3

    primStar = osti_efold_primitive(x,phi0)
    epsOneStar = osti_epsilon_one(x,phi0)
    potStar = osti_norm_potential(x)

    find_osti_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_osti_x_star


  function osti_lnrhoreh_max(phi0,Pstar) 
    implicit none
    real(kp) :: osti_lnrhoreh_max
    real(kp), intent(in) :: phi0,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = osti_x_endinf(phi0)
    potEnd  = osti_norm_potential(xEnd)
    epsOneEnd = osti_epsilon_one(xEnd,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = osti_x_star(phi0,wrad,junk,Pstar)  
 
    potStar = osti_norm_potential(x)
    epsOneStar = osti_epsilon_one(x,phi0)
   
    if (.not.slowroll_validity(epsOneStar)) stop 'osti_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    osti_lnrhoreh_max = lnRhoEnd

  end function osti_lnrhoreh_max
  
end module ostireheat

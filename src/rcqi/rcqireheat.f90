!radiatively corrected quartic inflation reheating functions in the slow-roll approximations

module rcqireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use rcqisr, only : rcqi_epsilon_one, rcqi_epsilon_two, rcqi_epsilon_three
  use rcqisr, only : rcqi_norm_potential
  use rcqisr, only : rcqi_x_endinf, rcqi_efold_primitive
  implicit none

  private

  public rcqi_x_star, rcqi_lnrhoreh_max 

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function rcqi_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rcqi_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: rcqiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = rcqi_x_endinf(alpha)
    epsOneEnd = rcqi_epsilon_one(xEnd,alpha)
    potEnd = rcqi_norm_potential(xEnd,alpha)
    primEnd = rcqi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rcqiData%real1 = alpha    
    rcqiData%real2 = w
    rcqiData%real3 = calF + primEnd

    mini = xEnd
    maxi = min(1._kp/epsilon(1._kp),exp(-0.25_kp+1._kp/alpha))

    x = zbrent(find_rcqi_x_star,mini,maxi,tolzbrent,rcqiData)
    rcqi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rcqi_efold_primitive(x,alpha) - primEnd)
    endif

  end function rcqi_x_star

  function find_rcqi_x_star(x,rcqiData)   
    implicit none
    real(kp) :: find_rcqi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcqiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=rcqiData%real1
    w = rcqiData%real2
    CalFplusprimEnd = rcqiData%real3

    primStar = rcqi_efold_primitive(x,alpha)
    epsOneStar = rcqi_epsilon_one(x,alpha)
    potStar = rcqi_norm_potential(x,alpha)

    find_rcqi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_rcqi_x_star



  function rcqi_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: rcqi_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = rcqi_x_endinf(alpha)
    potEnd  = rcqi_norm_potential(xEnd,alpha)
    epsOneEnd = rcqi_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = rcqi_x_star(alpha,wrad,junk,Pstar)    
    potStar = rcqi_norm_potential(x,alpha)
    epsOneStar = rcqi_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'rcqi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    rcqi_lnrhoreh_max = lnRhoEnd

  end function rcqi_lnrhoreh_max

  
end module rcqireheat

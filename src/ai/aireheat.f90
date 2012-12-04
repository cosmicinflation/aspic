!ArcTan inflation reheating functions in the slow-roll approximations

module aireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use aisr, only : ai_epsilon_one, ai_epsilon_two, ai_epsilon_three
  use aisr, only : ai_norm_potential
  use aisr, only : ai_x_endinf, ai_efold_primitive, ai_x_trajectory
  implicit none

  private

  public ai_x_star, ai_lnrhoend 

contains

!returns x=phi/mu such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ai_x_star(mu,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ai_x_star
    real(kp), intent(in) :: mu,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: aiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = ai_x_endinf(mu)

    epsOneEnd = ai_epsilon_one(xEnd,mu)
    potEnd = ai_norm_potential(xEnd)

    primEnd = ai_efold_primitive(xEnd,mu)

    
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)


    aiData%real1 = mu
    aiData%real2 = w
    aiData%real3 = calF + primEnd

    mini = -mu**(-2._kp/3._kp)*10._kp**(10._kp)
    maxi = xend*(1._kp-epsilon(1._kp))


    x = zbrent(find_ai_x_star,mini,maxi,tolzbrent,aiData)
    ai_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (ai_efold_primitive(x,mu) - primEnd)
    endif

  end function ai_x_star

  function find_ai_x_star(x,aiData)   
    implicit none
    real(kp) :: find_ai_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: aiData

    real(kp) :: primStar,mu,w,CalFplusprimEnd,potStar,epsOneStar

    mu = aiData%real1
    w = aiData%real2
    CalFplusprimEnd = aiData%real3

    primStar = ai_efold_primitive(x,mu)
    epsOneStar = ai_epsilon_one(x,mu)
    potStar = ai_norm_potential(x)

    find_ai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ai_x_star



  function ai_lnrhoend(mu,Pstar) 
    implicit none
    real(kp) :: ai_lnrhoend
    real(kp), intent(in) :: mu,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1_kp/3_kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = ai_x_endinf(mu)

    potEnd  = ai_norm_potential(xEnd)

    epsOneEnd = ai_epsilon_one(xEnd,mu)



!   Trick to return x such that rho_reh=rho_end

    x = ai_x_star(mu,wrad,junk,Pstar)  

    potStar = ai_norm_potential(x)
    epsOneStar = ai_epsilon_one(x,mu)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'ai_lnrhoend: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ai_lnrhoend = lnRhoEnd

  end function ai_lnrhoend

  
end module aireheat

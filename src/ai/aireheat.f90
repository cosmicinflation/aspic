!ArcTan inflation reheating functions in the slow-roll approximations

module aireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use aisr, only : ai_epsilon_one, ai_epsilon_two, ai_epsilon_three
  use aisr, only : ai_norm_potential, ai_numacc_xinimin
  use aisr, only : ai_x_endinf, ai_efold_primitive
  implicit none

  private

  public ai_x_star, ai_lnrhoreh_max
  public ai_x_rrad, ai_x_rreh

contains

!returns x=phi/mu such potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function ai_x_star(mu,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ai_x_star
    real(kp), intent(in) :: mu,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: aiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ai_epsilon_one(xEnd,mu)
    potEnd = ai_norm_potential(xEnd,mu)
    primEnd = ai_efold_primitive(xEnd,mu)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    aiData%real1 = mu
    aiData%real2 = w
    aiData%real3 = calF + primEnd

!    mini = -mu**(-2._kp/3._kp)*10._kp**(10._kp)
!    maxi = xend*(1._kp-epsilon(1._kp))
    mini = ai_numacc_xinimin(mu)
    maxi = xend


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
    potStar = ai_norm_potential(x,mu)

    find_ai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_ai_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ai_x_rrad(mu,xend,lnRRad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: ai_x_rrad
    real(kp), intent(in) :: mu,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: aiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ai_epsilon_one(xEnd,mu)
    potEnd = ai_norm_potential(xEnd,mu)

    primEnd = ai_efold_primitive(xEnd,mu)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)


    aiData%real1 = mu
    aiData%real2 = calF + primEnd

!    mini = -mu**(-2._kp/3._kp)*10._kp**(10._kp)
!    maxi = xend*(1._kp-epsilon(1._kp))
    mini = ai_numacc_xinimin(mu)
    maxi = xend

    x = zbrent(find_ai_x_rrad,mini,maxi,tolzbrent,aiData)
    ai_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (ai_efold_primitive(x,mu) - primEnd)
    endif

  end function ai_x_rrad

  function find_ai_x_rrad(x,aiData)   
    implicit none
    real(kp) :: find_ai_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: aiData

    real(kp) :: primStar,mu,CalFplusprimEnd,potStar,epsOneStar

    mu = aiData%real1
    CalFplusprimEnd = aiData%real2

    primStar = ai_efold_primitive(x,mu)
    epsOneStar = ai_epsilon_one(x,mu)
    potStar = ai_norm_potential(x,mu)

    find_ai_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd &
         ,epsOneStar,potStar)
  
  end function find_ai_x_rrad


  !returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function ai_x_rreh(mu,xend,lnRReh,bfoldstar)    
    implicit none
    real(kp) :: ai_x_rreh
    real(kp), intent(in) :: mu,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: aiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = ai_epsilon_one(xEnd,mu)
    potEnd = ai_norm_potential(xEnd,mu)

    primEnd = ai_efold_primitive(xEnd,mu)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)


    aiData%real1 = mu
    aiData%real2 = calF + primEnd

!    mini = -mu**(-2._kp/3._kp)*10._kp**(10._kp)
!    maxi = xend*(1._kp-epsilon(1._kp))
    mini = ai_numacc_xinimin(mu)
    maxi = xend

    x = zbrent(find_ai_x_rreh,mini,maxi,tolzbrent,aiData)
    ai_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (ai_efold_primitive(x,mu) - primEnd)
    endif

  end function ai_x_rreh

  function find_ai_x_rreh(x,aiData)   
    implicit none
    real(kp) :: find_ai_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: aiData

    real(kp) :: primStar,mu,CalFplusprimEnd,potStar

    mu = aiData%real1
    CalFplusprimEnd = aiData%real2

    primStar = ai_efold_primitive(x,mu)    
    potStar = ai_norm_potential(x,mu)

    find_ai_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd &
         ,potStar)
  
  end function find_ai_x_rreh


  function ai_lnrhoreh_max(mu,xend,Pstar) 
    implicit none
    real(kp) :: ai_lnrhoreh_max
    real(kp), intent(in) :: mu,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = ai_norm_potential(xEnd,mu)
    epsOneEnd = ai_epsilon_one(xEnd,mu)

!   Trick to return x such that rho_reh=rho_end

    x = ai_x_star(mu,xend,wrad,junk,Pstar)  

    potStar = ai_norm_potential(x,mu)
    epsOneStar = ai_epsilon_one(x,mu)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'ai_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    ai_lnrhoreh_max = lnRhoEnd

  end function ai_lnrhoreh_max

  
end module aireheat

!Higgs inflation reheating functions in the slow-roll approximations

module hireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hisr, only : hi_epsilon_one, hi_epsilon_two, hi_epsilon_three
  use hisr, only : hi_norm_potential
  use hisr, only : hi_x_endinf, hi_efold_primitive
  implicit none

  private

  public hi_x_star, hi_lnrhoreh_max
  public hi_x_rrad, hi_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function hi_x_star(w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: hi_x_star
    real(kp), intent(in) :: lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: hiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = hi_x_endinf()

    epsOneEnd = hi_epsilon_one(xEnd)
    potEnd = hi_norm_potential(xEnd)

    primEnd = hi_efold_primitive(xEnd)

    
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    hiData%real1 = w
    hiData%real2 = calF + primEnd

    mini = hi_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_hi_x_star,mini,maxi,tolzbrent,hiData)
    hi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (hi_efold_primitive(x) - primEnd)
    endif

  end function hi_x_star

  function find_hi_x_star(x,hiData)   
    implicit none
    real(kp) :: find_hi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: primStar,w,CalFplusprimEnd,potStar,epsOneStar

    w = hiData%real1
    CalFplusprimEnd = hiData%real2

    primStar = hi_efold_primitive(x)
    epsOneStar = hi_epsilon_one(x)
    potStar = hi_norm_potential(x)

    find_hi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_hi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hi_x_rrad(lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: hi_x_rrad
    real(kp), intent(in) :: lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: hiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    
    xEnd = hi_x_endinf()

    epsOneEnd = hi_epsilon_one(xEnd)
    potEnd = hi_norm_potential(xEnd)

    primEnd = hi_efold_primitive(xEnd)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    hiData%real1 = calF + primEnd

    mini = hi_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_hi_x_rrad,mini,maxi,tolzbrent,hiData)
    hi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (hi_efold_primitive(x) - primEnd)
    endif

  end function hi_x_rrad

  function find_hi_x_rrad(x,hiData)   
    implicit none
    real(kp) :: find_hi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: primStar,CalFplusprimEnd,potStar,epsOneStar

    CalFplusprimEnd = hiData%real1

    primStar = hi_efold_primitive(x)
    epsOneStar = hi_epsilon_one(x)
    potStar = hi_norm_potential(x)

    find_hi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_hi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hi_x_rreh(lnRreh,bfoldstar)    
    implicit none
    real(kp) :: hi_x_rreh
    real(kp), intent(in) :: lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: hiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
  
    xEnd = hi_x_endinf()

    epsOneEnd = hi_epsilon_one(xEnd)
    potEnd = hi_norm_potential(xEnd)

    primEnd = hi_efold_primitive(xEnd)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    hiData%real1 = calF + primEnd

    mini = hi_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_hi_x_rreh,mini,maxi,tolzbrent,hiData)
    hi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (hi_efold_primitive(x) - primEnd)
    endif

  end function hi_x_rreh

  function find_hi_x_rreh(x,hiData)   
    implicit none
    real(kp) :: find_hi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hiData

    real(kp) :: primStar,CalFplusprimEnd,potStar

    CalFplusprimEnd = hiData%real1

    primStar = hi_efold_primitive(x)    
    potStar = hi_norm_potential(x)

    find_hi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_hi_x_rreh



  function hi_lnrhoreh_max(Pstar) 
    implicit none
    real(kp) :: hi_lnrhoreh_max
    real(kp), intent(in) :: Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = hi_x_endinf()

    potEnd  = hi_norm_potential(xEnd)

    epsOneEnd = hi_epsilon_one(xEnd)

!   Trick to return x such that rho_reh=rho_end

    x = hi_x_star(wrad,junk,Pstar)  

    potStar = hi_norm_potential(x)
    epsOneStar = hi_epsilon_one(x)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'hi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hi_lnrhoreh_max = lnRhoEnd

  end function hi_lnrhoreh_max

  
end module hireheat

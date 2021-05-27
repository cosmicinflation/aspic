!Starobinsky inflation reheating functions in the slow-roll approximations

module sireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sisr, only : si_epsilon_one, si_epsilon_two, si_epsilon_three
  use sisr, only : si_norm_potential
  use sisr, only : si_x_endinf, si_efold_primitive
  implicit none

  private

  public si_x_star, si_lnrhoreh_max
  public si_x_rrad, si_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function si_x_star(xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: si_x_star
    real(kp), intent(in) :: xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd, lnOmega4End

    type(transfert) :: siData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = si_epsilon_one(xEnd)
    potEnd = si_norm_potential(xEnd)

    primEnd = si_efold_primitive(xEnd)
!ln(Omega^4) = ln(F^2)
    lnOmega4End = 2._kp*sqrt(2._kp/3._kp)*xend
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd,lnOmega4End)

    siData%real1 = w
    siData%real2 = calF + primEnd

    mini = si_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_si_x_star,mini,maxi,tolzbrent,siData)
    si_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (si_efold_primitive(x) - primEnd)
    endif

  end function si_x_star

  function find_si_x_star(x,siData)   
    implicit none
    real(kp) :: find_si_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: siData

    real(kp) :: primStar,w,CalFplusprimEnd,potStar,epsOneStar

    w = siData%real1
    CalFplusprimEnd = siData%real2

    primStar = si_efold_primitive(x)
    epsOneStar = si_epsilon_one(x)
    potStar = si_norm_potential(x)

    find_si_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_si_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function si_x_rrad(xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: si_x_rrad
    real(kp), intent(in) :: xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: siData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
        

    epsOneEnd = si_epsilon_one(xEnd)
    potEnd = si_norm_potential(xEnd)

    primEnd = si_efold_primitive(xEnd)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    siData%real1 = calF + primEnd

    mini = si_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_si_x_rrad,mini,maxi,tolzbrent,siData)
    si_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (si_efold_primitive(x) - primEnd)
    endif

  end function si_x_rrad

  function find_si_x_rrad(x,siData)   
    implicit none
    real(kp) :: find_si_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: siData

    real(kp) :: primStar,CalFplusprimEnd,potStar,epsOneStar

    CalFplusprimEnd = siData%real1

    primStar = si_efold_primitive(x)
    epsOneStar = si_epsilon_one(x)
    potStar = si_norm_potential(x)

    find_si_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_si_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function si_x_rreh(xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: si_x_rreh
    real(kp), intent(in) :: xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd, lnOmega4End

    type(transfert) :: siData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif  

    epsOneEnd = si_epsilon_one(xEnd)
    potEnd = si_norm_potential(xEnd)

    primEnd = si_efold_primitive(xEnd)

!ln(Omega^4) = ln(F^2)
    lnOmega4End = 2._kp*sqrt(2._kp/3._kp)*xend

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd,lnOmega4End)

    siData%real1 = calF + primEnd

    mini = si_x_endinf()
    maxi = 1./epsilon(1._kp)

    x = zbrent(find_si_x_rreh,mini,maxi,tolzbrent,siData)
    si_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (si_efold_primitive(x) - primEnd)
    endif

  end function si_x_rreh

  function find_si_x_rreh(x,siData)   
    implicit none
    real(kp) :: find_si_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: siData

    real(kp) :: primStar,CalFplusprimEnd,potStar

    CalFplusprimEnd = siData%real1

    primStar = si_efold_primitive(x)    
    potStar = si_norm_potential(x)

    find_si_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_si_x_rreh



  function si_lnrhoreh_max(xend,Pstar) 
    implicit none
    real(kp) :: si_lnrhoreh_max
    real(kp), intent(in) :: xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnOmega4End, lnRhoEndJF

    potEnd  = si_norm_potential(xEnd)

    epsOneEnd = si_epsilon_one(xEnd)

!   Trick to return x such that rho_reh=rho_end

    x = si_x_star(xend,wrad,junk,Pstar)  

    potStar = si_norm_potential(x)
    epsOneStar = si_epsilon_one(x)

    !ln(Omega^4) = ln(F^2)
    lnOmega4End = 2._kp*sqrt(2._kp/3._kp)*xend
    
    if (.not.slowroll_validity(epsOneStar)) stop 'si_lnrhoreh_max: slow-roll violated!'
       
    lnRhoEndJF = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar,lnOmega4End)

    si_lnrhoreh_max = lnRhoEndJF

  end function si_lnrhoreh_max

  
end module sireheat

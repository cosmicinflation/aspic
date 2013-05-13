!power law inflation reheating functions in the slow-roll approximations

module plireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use plisr, only : pli_epsilon_one, pli_epsilon_two, pli_epsilon_three
  use plisr, only : pli_norm_potential, pli_x_endinf
  use plisr, only : pli_efold_primitive
  implicit none

  private

  public pli_x_star, pli_lnrhoreh_max
  public pli_x_rrad, pli_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function pli_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pli_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: pliData    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = pli_x_endinf(alpha)

    epsOneEnd = pli_epsilon_one(xEnd,alpha)
    potEnd = pli_norm_potential(xEnd,alpha)
    primEnd = pli_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    pliData%real1 = alpha    
    pliData%real2 = w
    pliData%real3 = calF + primEnd

    mini = 1.
    maxi = xend

    x = zbrent(find_pli_x_star,mini,maxi,tolzbrent,pliData)
    pli_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (pli_efold_primitive(x,alpha) - primEnd)
    endif

  end function pli_x_star

  function find_pli_x_star(x,pliData)   
    implicit none
    real(kp) :: find_pli_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: pliData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=pliData%real1
    w = pliData%real2
    CalFplusprimEnd = pliData%real3

    primStar = pli_efold_primitive(x,alpha)
    epsOneStar = pli_epsilon_one(x,alpha)
    potStar = pli_norm_potential(x,alpha)

    find_pli_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_pli_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function pli_x_rrad(alpha,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: pli_x_rrad
    real(kp), intent(in) :: alpha,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: pliData    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = pli_x_endinf(alpha)

    epsOneEnd = pli_epsilon_one(xEnd,alpha)
    potEnd = pli_norm_potential(xEnd,alpha)
    primEnd = pli_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    pliData%real1 = alpha    
    pliData%real2 = calF + primEnd

    mini = 1.
    maxi = xend

    x = zbrent(find_pli_x_rrad,mini,maxi,tolzbrent,pliData)
    pli_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (pli_efold_primitive(x,alpha) - primEnd)
    endif

  end function pli_x_rrad

  function find_pli_x_rrad(x,pliData)   
    implicit none
    real(kp) :: find_pli_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: pliData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=pliData%real1
    CalFplusprimEnd = pliData%real2

    primStar = pli_efold_primitive(x,alpha)
    epsOneStar = pli_epsilon_one(x,alpha)
    potStar = pli_norm_potential(x,alpha)

    find_pli_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_pli_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function pli_x_rreh(alpha,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: pli_x_rreh
    real(kp), intent(in) :: alpha,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: pliData    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = pli_x_endinf(alpha)

    epsOneEnd = pli_epsilon_one(xEnd,alpha)
    potEnd = pli_norm_potential(xEnd,alpha)
    primEnd = pli_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    pliData%real1 = alpha    
    pliData%real2 = calF + primEnd

    mini = 1.
    maxi = xend

    x = zbrent(find_pli_x_rreh,mini,maxi,tolzbrent,pliData)
    pli_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (pli_efold_primitive(x,alpha) - primEnd)
    endif

  end function pli_x_rreh

  function find_pli_x_rreh(x,pliData)   
    implicit none
    real(kp) :: find_pli_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: pliData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=pliData%real1
    CalFplusprimEnd = pliData%real2

    primStar = pli_efold_primitive(x,alpha)
    potStar = pli_norm_potential(x,alpha)

    find_pli_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_pli_x_rreh



  function pli_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: pli_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp


    real(kp) :: lnRhoEnd
    

    xEnd=pli_x_endinf(alpha)
    potEnd  = pli_norm_potential(xEnd,alpha)
    epsOneEnd = pli_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = pli_x_star(alpha,wrad,junk,Pstar)    
    potStar = pli_norm_potential(x,alpha)
    epsOneStar = pli_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'pli_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    pli_lnrhoreh_max = lnRhoEnd

  end function pli_lnrhoreh_max

  
end module plireheat

!Logamediate inflation 1 reheating functions in the slow-roll approximations

module lmi1reheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use lmicommon, only : lmi_x_potmax, lmi_numacc_dx_potmax
  use lmi1sr, only : lmi1_epsilon_one, lmi1_epsilon_two, lmi1_epsilon_three
  use lmi1sr, only : lmi1_norm_potential
  use lmi1sr, only : lmi1_x_endinf, lmi1_efold_primitive

  implicit none

  private

  public lmi1_x_star, lmi1_lnrhoreh_max
  public lmi1_x_rrad, lmi1_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function lmi1_x_star(gam,beta,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi1_x_star
    real(kp), intent(in) :: gam,beta,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xVmax
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: lmi1Data

  
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lmi1_x_endinf(gam,beta)
    xVmax = lmi_x_potmax(gam,beta)

    epsOneEnd = lmi1_epsilon_one(xEnd,gam,beta)
    potEnd = lmi1_norm_potential(xEnd,gam,beta)

    primEnd = lmi1_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    lmi1Data%real1 = gam
    lmi1Data%real2 = beta
    lmi1Data%real3 = w
    lmi1Data%real4 = calF + primEnd

    maxi = xvMax
    mini = xEnd

    x = zbrent(find_lmi1_x_star,mini,maxi,tolzbrent,lmi1Data)
    lmi1_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi1_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi1_x_star

  function find_lmi1_x_star(x,lmi1Data)   
    implicit none
    real(kp) :: find_lmi1_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi1Data

    real(kp) :: primStar,gam,beta,w,CalFplusprimEnd,potStar,epsOneStar

    gam=lmi1Data%real1
    beta=lmi1Data%real2
    w = lmi1Data%real3
    CalFplusprimEnd = lmi1Data%real4

    primStar = lmi1_efold_primitive(x,gam,beta)
    epsOneStar = lmi1_epsilon_one(x,gam,beta)
    potStar = lmi1_norm_potential(x,gam,beta)

    find_lmi1_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_lmi1_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function lmi1_x_rrad(gam,beta,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: lmi1_x_rrad
    real(kp), intent(in) :: gam,beta,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xVmax
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: lmi1Data

  
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lmi1_x_endinf(gam,beta)
    xVmax = lmi_x_potmax(gam,beta)

    epsOneEnd = lmi1_epsilon_one(xEnd,gam,beta)
    potEnd = lmi1_norm_potential(xEnd,gam,beta)

    primEnd = lmi1_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    lmi1Data%real1 = gam
    lmi1Data%real2 = beta
    lmi1Data%real3 = calF + primEnd

    maxi = xvMax
    mini = xEnd

    x = zbrent(find_lmi1_x_rrad,mini,maxi,tolzbrent,lmi1Data)
    lmi1_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi1_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi1_x_rrad

  function find_lmi1_x_rrad(x,lmi1Data)   
    implicit none
    real(kp) :: find_lmi1_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi1Data

    real(kp) :: primStar,gam,beta,CalFplusprimEnd,potStar,epsOneStar

    gam=lmi1Data%real1
    beta=lmi1Data%real2
    CalFplusprimEnd = lmi1Data%real3

    primStar = lmi1_efold_primitive(x,gam,beta)
    epsOneStar = lmi1_epsilon_one(x,gam,beta)
    potStar = lmi1_norm_potential(x,gam,beta)

    find_lmi1_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_lmi1_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function lmi1_x_rreh(gam,beta,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: lmi1_x_rreh
    real(kp), intent(in) :: gam,beta,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xVmax
    real(kp) :: primEnd,epsOneEnd,xend,potEnd
    type(transfert) :: lmi1Data

  
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = lmi1_x_endinf(gam,beta)
    xVmax = lmi_x_potmax(gam,beta)

    epsOneEnd = lmi1_epsilon_one(xEnd,gam,beta)
    potEnd = lmi1_norm_potential(xEnd,gam,beta)

    primEnd = lmi1_efold_primitive(xEnd,gam,beta)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    lmi1Data%real1 = gam
    lmi1Data%real2 = beta
    lmi1Data%real3 = calF + primEnd

    maxi = xvMax
    mini = xEnd

    x = zbrent(find_lmi1_x_rreh,mini,maxi,tolzbrent,lmi1Data)
    lmi1_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (lmi1_efold_primitive(x,gam,beta) - primEnd)
    endif

  end function lmi1_x_rreh

  function find_lmi1_x_rreh(x,lmi1Data)   
    implicit none
    real(kp) :: find_lmi1_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: lmi1Data

    real(kp) :: primStar,gam,beta,CalFplusprimEnd,potStar

    gam=lmi1Data%real1
    beta=lmi1Data%real2
    CalFplusprimEnd = lmi1Data%real3

    primStar = lmi1_efold_primitive(x,gam,beta)
    potStar = lmi1_norm_potential(x,gam,beta)

    find_lmi1_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_lmi1_x_rreh



  function lmi1_lnrhoreh_max(gam,beta,Pstar) 
    implicit none
    real(kp) :: lmi1_lnrhoreh_max
    real(kp), intent(in) :: gam,beta,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = lmi1_x_endinf(gam,beta)

    potEnd  = lmi1_norm_potential(xEnd,gam,beta)

    epsOneEnd = lmi1_epsilon_one(xEnd,gam,beta)


!   Trick to return x such that rho_reh=rho_end

    x = lmi1_x_star(gam,beta,wrad,junk,Pstar)  


    potStar = lmi1_norm_potential(x,gam,beta)
    epsOneStar = lmi1_epsilon_one(x,gam,beta)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'lmi1_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    lmi1_lnrhoreh_max = lnRhoEnd

  end function lmi1_lnrhoreh_max

  
end module lmi1reheat

!renormalizable inflection point reheating functions in the slow-roll approximations

module oripireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use oripisr, only : oripi_epsilon_one, oripi_epsilon_two, oripi_epsilon_three
  use oripisr, only : oripi_norm_potential, oripi_x_derivpotzero
  use oripisr, only : oripi_x_endinf, oripi_efold_primitive, oripi_x_trajectory
  implicit none

  private

  public oripi_x_star, oripi_lnrhoreh_max
  public oripi_x_rrad, oripi_x_rreh

contains

!returns x=phi/phi0 such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function oripi_x_star(phi0,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: oripi_x_star
    real(kp), intent(in) :: phi0,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oripiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = oripi_epsilon_one(xEnd,phi0)
    potEnd = oripi_norm_potential(xEnd,phi0)

    primEnd = oripi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    oripiData%real1 = phi0
    oripiData%real2 = w
    oripiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = oripi_x_derivpotzero(phi0)*(1._kp-epsilon(1._kp)) !Position of the flat inflection point


    x = zbrent(find_oripi_x_star,mini,maxi,tolzbrent,oripiData)
    oripi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (oripi_efold_primitive(x,phi0) - primEnd)
    endif

  end function oripi_x_star

  function find_oripi_x_star(x,oripiData)   
    implicit none
    real(kp) :: find_oripi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oripiData

    real(kp) :: primStar,phi0,w,CalFplusprimEnd,potStar,epsOneStar

    phi0 = oripiData%real1
    w = oripiData%real2
    CalFplusprimEnd = oripiData%real3

    primStar = oripi_efold_primitive(x,phi0)
    epsOneStar = oripi_epsilon_one(x,phi0)
    potStar = oripi_norm_potential(x,phi0)

    find_oripi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_oripi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function oripi_x_rrad(phi0,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: oripi_x_rrad
    real(kp), intent(in) :: phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oripiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = oripi_epsilon_one(xEnd,phi0)
    potEnd = oripi_norm_potential(xEnd,phi0)

    primEnd = oripi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    oripiData%real1 = phi0
    oripiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = oripi_x_derivpotzero(phi0)*(1._kp-epsilon(1._kp)) !Position of the flat inflection point


    x = zbrent(find_oripi_x_rrad,mini,maxi,tolzbrent,oripiData)
    oripi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (oripi_efold_primitive(x,phi0) - primEnd)
    endif

  end function oripi_x_rrad

  function find_oripi_x_rrad(x,oripiData)   
    implicit none
    real(kp) :: find_oripi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oripiData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar,epsOneStar

    phi0 = oripiData%real1
    CalFplusprimEnd = oripiData%real2

    primStar = oripi_efold_primitive(x,phi0)
    epsOneStar = oripi_epsilon_one(x,phi0)
    potStar = oripi_norm_potential(x,phi0)

    find_oripi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_oripi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function oripi_x_rreh(phi0,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: oripi_x_rreh
    real(kp), intent(in) :: phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oripiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    

    epsOneEnd = oripi_epsilon_one(xEnd,phi0)
    potEnd = oripi_norm_potential(xEnd,phi0)

    primEnd = oripi_efold_primitive(xEnd,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    oripiData%real1 = phi0
    oripiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = oripi_x_derivpotzero(phi0)*(1._kp-epsilon(1._kp)) !Position of the flat inflection point


    x = zbrent(find_oripi_x_rreh,mini,maxi,tolzbrent,oripiData)
    oripi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (oripi_efold_primitive(x,phi0) - primEnd)
    endif

  end function oripi_x_rreh

  function find_oripi_x_rreh(x,oripiData)   
    implicit none
    real(kp) :: find_oripi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oripiData

    real(kp) :: primStar,phi0,CalFplusprimEnd,potStar

    phi0 = oripiData%real1
    CalFplusprimEnd = oripiData%real2

    primStar = oripi_efold_primitive(x,phi0)
    potStar = oripi_norm_potential(x,phi0)

    find_oripi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_oripi_x_rreh



  function oripi_lnrhoreh_max(phi0,xend,Pstar) 
    implicit none
    real(kp) :: oripi_lnrhoreh_max
    real(kp), intent(in) :: phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = oripi_norm_potential(xEnd,phi0)

    epsOneEnd = oripi_epsilon_one(xEnd,phi0)

!   Trick to return x such that rho_reh=rho_end

    x = oripi_x_star(phi0,xend,wrad,junk,Pstar)  

    potStar = oripi_norm_potential(x,phi0)
    epsOneStar = oripi_epsilon_one(x,phi0)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'oripi_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    oripi_lnrhoreh_max = lnRhoEnd

  end function oripi_lnrhoreh_max

  
end module oripireheat

!radiatively corrected massive inflation reheating functions in the slow-roll

module rcmireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use rcmisr, only : rcmi_epsilon_one, rcmi_epsilon_two
  use rcmisr, only : rcmi_norm_potential, rcmi_x_potmax
  use rcmisr, only : rcmi_x_endinf, rcmi_efold_primitive
  use specialinf, only : lambert
  implicit none

  private

  public rcmi_x_star, rcmi_lnrhoreh_max
  public rcmi_x_rrad, rcmi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function rcmi_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rcmi_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xPotMax
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rcmiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xPotMax = rcmi_x_potmax(alpha)
    epsOneEnd = rcmi_epsilon_one(xEnd,alpha)
    potEnd = rcmi_norm_potential(xEnd,alpha)
    primEnd = rcmi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    rcmiData%real1 = alpha
    rcmiData%real2 = w
    rcmiData%real3 = calF + primEnd

    mini = xEnd
    maxi = min(xPotMax,1._kp/epsilon(1._kp))

    x = zbrent(find_rcmi_x_star,mini,maxi,tolzbrent,rcmiData)
    rcmi_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (rcmi_efold_primitive(x,alpha) - primEnd)
    endif
    
  end function rcmi_x_star

  function find_rcmi_x_star(x,rcmiData)   
    implicit none
    real(kp) :: find_rcmi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcmiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha = rcmiData%real1
    w = rcmiData%real2
    CalFplusprimEnd = rcmiData%real3

    primStar = rcmi_efold_primitive(x,alpha)
    epsOneStar = rcmi_epsilon_one(x,alpha)
    potStar = rcmi_norm_potential(x,alpha)

    find_rcmi_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)

  end function find_rcmi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function rcmi_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: rcmi_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xPotMax
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rcmiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xPotMax = rcmi_x_potmax(alpha)
    epsOneEnd = rcmi_epsilon_one(xEnd,alpha)
    potEnd = rcmi_norm_potential(xEnd,alpha)
    primEnd = rcmi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    rcmiData%real1 = alpha
    rcmiData%real2 = calF + primEnd

    mini = xEnd
    maxi = min(xPotMax,1._kp/epsilon(1._kp))

    x = zbrent(find_rcmi_x_rrad,mini,maxi,tolzbrent,rcmiData)
    rcmi_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (rcmi_efold_primitive(x,alpha) - primEnd)
    endif
    
  end function rcmi_x_rrad

  function find_rcmi_x_rrad(x,rcmiData)   
    implicit none
    real(kp) :: find_rcmi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcmiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha = rcmiData%real1
    CalFplusprimEnd = rcmiData%real2

    primStar = rcmi_efold_primitive(x,alpha)
    epsOneStar = rcmi_epsilon_one(x,alpha)
    potStar = rcmi_norm_potential(x,alpha)

    find_rcmi_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)

  end function find_rcmi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function rcmi_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: rcmi_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xPotMax
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: rcmiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xPotMax = rcmi_x_potmax(alpha)
    epsOneEnd = rcmi_epsilon_one(xEnd,alpha)
    potEnd = rcmi_norm_potential(xEnd,alpha)
    primEnd = rcmi_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    rcmiData%real1 = alpha
    rcmiData%real2 = calF + primEnd

    mini = xEnd
    maxi = min(xPotMax,1._kp/epsilon(1._kp))

    x = zbrent(find_rcmi_x_rreh,mini,maxi,tolzbrent,rcmiData)
    rcmi_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (rcmi_efold_primitive(x,alpha) - primEnd)
    endif
    
  end function rcmi_x_rreh

  function find_rcmi_x_rreh(x,rcmiData)   
    implicit none
    real(kp) :: find_rcmi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: rcmiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha = rcmiData%real1
    CalFplusprimEnd = rcmiData%real2

    primStar = rcmi_efold_primitive(x,alpha)
    potStar = rcmi_norm_potential(x,alpha)

    find_rcmi_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)

  end function find_rcmi_x_rreh



  function rcmi_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: rcmi_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk = 0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = rcmi_norm_potential(xEnd,alpha)
    epsEnd = rcmi_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end
       
    x = rcmi_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = rcmi_norm_potential(x,alpha)
    epsOneStar = rcmi_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'rcmi_lnrhoreh_max: slow-roll violated!'

    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsEnd,potEnd/potStar)

    rcmi_lnrhoreh_max = lnRhoEnd

  end function rcmi_lnrhoreh_max

  

    
end module rcmireheat

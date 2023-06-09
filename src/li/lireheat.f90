!loop inflation reheating functions in the slow-roll approximations

module lireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use specialinf, only : lambert
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use lisr, only : li_epsilon_one, li_epsilon_two, li_epsilon_three
  use lisr, only : li_norm_potential
  use lisr, only : li_x_endinf, li_efold_primitive, li_x_epsoneunity
  implicit none

  private

  public li_x_star, li_lnrhoreh_max
  public li_x_rrad, li_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function li_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: li_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsones
    type(transfert) :: liData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    potEnd = li_norm_potential(xEnd,alpha)
    primEnd = li_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    liData%real1 = alpha
    liData%real2 = w
    liData%real3 = calF + primEnd

    xepsones = li_x_epsoneunity(alpha)

    if (alpha .gt. 0._kp) then
       mini = xEnd
       maxi = huge(1._kp)*epsilon(1._kp)
    else
       mini = xepsones(1)
       maxi = xEnd
    endif

    x = zbrent(find_li_x_star,mini,maxi,tolzbrent,liData)
    li_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (li_efold_primitive(x,alpha) - primEnd)
    endif

  end function li_x_star

  function find_li_x_star(x,liData)   
    implicit none
    real(kp) :: find_li_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: liData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=liData%real1
    w = liData%real2
    CalFplusprimEnd = liData%real3

    primStar = li_efold_primitive(x,alpha)
    epsOneStar = li_epsilon_one(x,alpha)
    potStar = li_norm_potential(x,alpha)

    find_li_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_li_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function li_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: li_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsones
    type(transfert) :: liData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif    
    
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    potEnd = li_norm_potential(xEnd,alpha)
    primEnd = li_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    liData%real1 = alpha
    liData%real2 = calF + primEnd

    xepsones = li_x_epsoneunity(alpha)

    if (alpha .gt. 0._kp) then
       mini = xEnd
       maxi = huge(1._kp)*epsilon(1._kp)
    else
       mini = xepsones(1)
       maxi = xEnd
    endif   

    x = zbrent(find_li_x_rrad,mini,maxi,tolzbrent,liData)
    li_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (li_efold_primitive(x,alpha) - primEnd)
    endif

  end function li_x_rrad

  function find_li_x_rrad(x,liData)   
    implicit none
    real(kp) :: find_li_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: liData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=liData%real1
    CalFplusprimEnd = liData%real2

    primStar = li_efold_primitive(x,alpha)
    epsOneStar = li_epsilon_one(x,alpha)
    potStar = li_norm_potential(x,alpha)

    find_li_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_li_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function li_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: li_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsones
    type(transfert) :: liData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif    
    
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    potEnd = li_norm_potential(xEnd,alpha)
    primEnd = li_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    liData%real1 = alpha
    liData%real2 = calF + primEnd

    xepsones = li_x_epsoneunity(alpha)

    if (alpha .gt. 0._kp) then
       mini = xEnd
       maxi = huge(1._kp)*epsilon(1._kp)
    else
       mini = xepsones(1)
       maxi = xEnd
    endif

    x = zbrent(find_li_x_rreh,mini,maxi,tolzbrent,liData)
    li_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (li_efold_primitive(x,alpha) - primEnd)
    endif

  end function li_x_rreh

  function find_li_x_rreh(x,liData)   
    implicit none
    real(kp) :: find_li_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: liData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=liData%real1
    CalFplusprimEnd = liData%real2

    primStar = li_efold_primitive(x,alpha)
    potStar = li_norm_potential(x,alpha)

    find_li_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_li_x_rreh



  function li_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: li_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = li_norm_potential(xEnd,alpha)
    epsOneEnd = li_epsilon_one(xEnd,alpha)
    

!   Trick to return x such that rho_reh=rho_end

    x = li_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = li_norm_potential(x,alpha)
    epsOneStar = li_epsilon_one(x,alpha)

    
    if (.not.slowroll_validity(epsOneStar)) stop 'li_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    li_lnrhoreh_max = lnRhoEnd

  end function li_lnrhoreh_max

  
end module lireheat

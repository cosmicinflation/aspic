!SuperConformal alpha Attractor A reheating functions in the slow-roll approximations

module saaireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use saaisr, only : saai_epsilon_one, saai_epsilon_two, saai_epsilon_three
  use saaisr, only : saai_norm_potential
  use saaisr, only : saai_x_endinf, saai_efold_primitive
  implicit none

  private

  public saai_x_star, saai_lnrhoreh_max 
  public saai_x_rrad, saai_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function saai_x_star(alpha,xend,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saai_x_star
    real(kp), intent(in) :: alpha,xend,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: saaiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = saai_epsilon_one(xEnd,alpha)
    potEnd = saai_norm_potential(xEnd,alpha)
    primEnd = saai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    saaiData%real1 = alpha    
    saaiData%real2 = w
    saaiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = 2._kp*sqrt(2._kp/(3._kp*alpha))*200._kp

    x = zbrent(find_saai_x_star,mini,maxi,tolzbrent,saaiData)
    saai_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (saai_efold_primitive(x,alpha) - primEnd)
    endif

  end function saai_x_star

  function find_saai_x_star(x,saaiData)   
    implicit none
    real(kp) :: find_saai_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saaiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=saaiData%real1
    w = saaiData%real2
    CalFplusprimEnd = saaiData%real3

    primStar = saai_efold_primitive(x,alpha)
    epsOneStar = saai_epsilon_one(x,alpha)
    potStar = saai_norm_potential(x,alpha)

    find_saai_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_saai_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function saai_x_rrad(alpha,xend,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: saai_x_rrad
    real(kp), intent(in) :: alpha,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: saaiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = saai_epsilon_one(xEnd,alpha)
    potEnd = saai_norm_potential(xEnd,alpha)
    primEnd = saai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    saaiData%real1 = alpha
    saaiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = 2._kp*sqrt(2._kp/(3._kp*alpha))*200._kp

    x = zbrent(find_saai_x_rrad,mini,maxi,tolzbrent,saaiData)
    saai_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (saai_efold_primitive(x,alpha) - primEnd)
    endif

  end function saai_x_rrad

  function find_saai_x_rrad(x,saaiData)   
    implicit none
    real(kp) :: find_saai_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saaiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=saaiData%real1
    CalFplusprimEnd = saaiData%real2

    primStar = saai_efold_primitive(x,alpha)
    epsOneStar = saai_epsilon_one(x,alpha)
    potStar = saai_norm_potential(x,alpha)

    find_saai_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_saai_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function saai_x_rreh(alpha,xend,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: saai_x_rreh
    real(kp), intent(in) :: alpha,xend,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: saaiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = saai_epsilon_one(xEnd,alpha)
    potEnd = saai_norm_potential(xEnd,alpha)
    primEnd = saai_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    saaiData%real1 = alpha
    saaiData%real2 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = 2._kp*sqrt(2._kp/(3._kp*alpha))*200._kp

    x = zbrent(find_saai_x_rreh,mini,maxi,tolzbrent,saaiData)
    saai_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (saai_efold_primitive(x,alpha) - primEnd)
    endif

  end function saai_x_rreh

  function find_saai_x_rreh(x,saaiData)   
    implicit none
    real(kp) :: find_saai_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: saaiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=saaiData%real1
    CalFplusprimEnd = saaiData%real2

    primStar = saai_efold_primitive(x,alpha)
    potStar = saai_norm_potential(x,alpha)

    find_saai_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_saai_x_rreh




  function saai_lnrhoreh_max(alpha,xend,Pstar) 
    implicit none
    real(kp) :: saai_lnrhoreh_max
    real(kp), intent(in) :: alpha,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    potEnd  = saai_norm_potential(xEnd,alpha)
    epsOneEnd = saai_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = saai_x_star(alpha,xend,wrad,junk,Pstar)    
    potStar = saai_norm_potential(x,alpha)
    epsOneStar = saai_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'saai_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    saai_lnrhoreh_max = lnRhoEnd

  end function saai_lnrhoreh_max

  
end module saaireheat

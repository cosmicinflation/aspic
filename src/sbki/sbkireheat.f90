!exponential SUSY inflation reheating functions in the slow-roll approximations

module sbkireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use sbkisr, only : sbki_epsilon_one, sbki_epsilon_two, sbki_epsilon_three
  use sbkisr, only : sbki_norm_potential, sbki_x_max
  use sbkisr, only : sbki_x_endinf, sbki_efold_primitive
  implicit none

  private

  public sbki_x_star, sbki_lnrhoreh_max 
  public sbki_x_rrad, sbki_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function sbki_x_star(alpha,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sbki_x_star
    real(kp), intent(in) :: alpha,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: sbkiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = sbki_x_endinf(alpha)
    epsOneEnd = sbki_epsilon_one(xEnd,alpha)
    potEnd = sbki_norm_potential(xEnd,alpha)
    primEnd = sbki_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    sbkiData%real1 = alpha    
    sbkiData%real2 = w
    sbkiData%real3 = calF + primEnd

    xplus = sbki_x_max(alpha)

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = xplus*(1._kp-epsilon(1._kp))

    x = zbrent(find_sbki_x_star,mini,maxi,tolzbrent,sbkiData)
    sbki_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbki_efold_primitive(x,alpha) - primEnd)
    endif

  end function sbki_x_star

  function find_sbki_x_star(x,sbkiData)   
    implicit none
    real(kp) :: find_sbki_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbkiData

    real(kp) :: primStar,alpha,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=sbkiData%real1
    w = sbkiData%real2
    CalFplusprimEnd = sbkiData%real3

    primStar = sbki_efold_primitive(x,alpha)
    epsOneStar = sbki_epsilon_one(x,alpha)
    potStar = sbki_norm_potential(x,alpha)

    find_sbki_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_sbki_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function sbki_x_rrad(alpha,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: sbki_x_rrad
    real(kp), intent(in) :: alpha,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: sbkiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = sbki_x_endinf(alpha)
    epsOneEnd = sbki_epsilon_one(xEnd,alpha)
    potEnd = sbki_norm_potential(xEnd,alpha)
    primEnd = sbki_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    sbkiData%real1 = alpha
    sbkiData%real2 = calF + primEnd

    xplus = sbki_x_max(alpha)

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = xplus*(1._kp-epsilon(1._kp))

    x = zbrent(find_sbki_x_rrad,mini,maxi,tolzbrent,sbkiData)
    sbki_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbki_efold_primitive(x,alpha) - primEnd)
    endif

  end function sbki_x_rrad

  function find_sbki_x_rrad(x,sbkiData)   
    implicit none
    real(kp) :: find_sbki_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbkiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar,epsOneStar

    alpha=sbkiData%real1
    CalFplusprimEnd = sbkiData%real2

    primStar = sbki_efold_primitive(x,alpha)
    epsOneStar = sbki_epsilon_one(x,alpha)
    potStar = sbki_norm_potential(x,alpha)

    find_sbki_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_sbki_x_rrad



!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function sbki_x_rreh(alpha,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: sbki_x_rreh
    real(kp), intent(in) :: alpha,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x,xplus
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: sbkiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = sbki_x_endinf(alpha)
    epsOneEnd = sbki_epsilon_one(xEnd,alpha)
    potEnd = sbki_norm_potential(xEnd,alpha)
    primEnd = sbki_efold_primitive(xEnd,alpha)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    sbkiData%real1 = alpha
    sbkiData%real2 = calF + primEnd

    xplus = sbki_x_max(alpha)

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = xplus*(1._kp-epsilon(1._kp))

    x = zbrent(find_sbki_x_rreh,mini,maxi,tolzbrent,sbkiData)
    sbki_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (sbki_efold_primitive(x,alpha) - primEnd)
    endif

  end function sbki_x_rreh

  function find_sbki_x_rreh(x,sbkiData)   
    implicit none
    real(kp) :: find_sbki_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: sbkiData

    real(kp) :: primStar,alpha,CalFplusprimEnd,potStar

    alpha=sbkiData%real1
    CalFplusprimEnd = sbkiData%real2

    primStar = sbki_efold_primitive(x,alpha)
    potStar = sbki_norm_potential(x,alpha)

    find_sbki_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_sbki_x_rreh




  function sbki_lnrhoreh_max(alpha,Pstar) 
    implicit none
    real(kp) :: sbki_lnrhoreh_max
    real(kp), intent(in) :: alpha,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0._kp

    real(kp) :: lnRhoEnd
    
    xEnd = sbki_x_endinf(alpha)
    potEnd  = sbki_norm_potential(xEnd,alpha)
    epsOneEnd = sbki_epsilon_one(xEnd,alpha)

!   Trick to return x such that rho_reh=rho_end

    x = sbki_x_star(alpha,wrad,junk,Pstar)    
    potStar = sbki_norm_potential(x,alpha)
    epsOneStar = sbki_epsilon_one(x,alpha)
    
    if (.not.slowroll_validity(epsOneStar)) stop 'sbki_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    sbki_lnrhoreh_max = lnRhoEnd

  end function sbki_lnrhoreh_max

  
end module sbkireheat

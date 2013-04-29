!pseudo natural inflation reheating functions in the slow-roll approximations

module psnireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf
  use srreheat, only : ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use psnisr, only : psni_epsilon_one, psni_epsilon_two, psni_epsilon_three
  use psnisr, only : psni_norm_potential
  use psnisr, only : psni_x_endinf, psni_efold_primitive
  implicit none

  private

  public psni_x_star, psni_lnrhoreh_max
  public psni_x_rrad, psni_x_rreh

contains

!returns x such given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the corresponding bfoldstar
  function psni_x_star(alpha,f,w,lnRhoReh,Pstar,bfoldstar)    
    implicit none
    real(kp) :: psni_x_star
    real(kp), intent(in) :: alpha,f,lnRhoReh,w,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: psniData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    xEnd = psni_x_endinf(alpha,f)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,f)
    potEnd = psni_norm_potential(xEnd,alpha)
    primEnd = psni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    psniData%real1 = alpha 
    psniData%real2 = f
    psniData%real3 = w
    psniData%real4 = calF + primEnd

    mini=epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_psni_x_star,mini,maxi,tolzbrent,psniData)
    psni_x_star = x

    if (present(bfoldstar)) then
       bfoldstar = - (psni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function psni_x_star

  function find_psni_x_star(x,psniData)   
    implicit none
    real(kp) :: find_psni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: psniData

    real(kp) :: primStar,alpha,f,w,CalFplusprimEnd,potStar,epsOneStar

    alpha=psniData%real1
    f=psniData%real2
    w = psniData%real3
    CalFplusprimEnd = psniData%real4

    primStar = psni_efold_primitive(x,alpha,f)
    epsOneStar = psni_epsilon_one(x,alpha,f)
    potStar = psni_norm_potential(x,alpha)

    find_psni_x_star = find_reheat(primStar,calFplusprimEnd,w,epsOneStar,potStar)
  
  end function find_psni_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function psni_x_rrad(alpha,f,lnRrad,Pstar,bfoldstar)    
    implicit none
    real(kp) :: psni_x_rrad
    real(kp), intent(in) :: alpha,f,lnRrad,Pstar
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: psniData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = psni_x_endinf(alpha,f)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,f)
    potEnd = psni_norm_potential(xEnd,alpha)
    primEnd = psni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    psniData%real1 = alpha 
    psniData%real2 = f
    psniData%real3 = calF + primEnd

    mini=epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_psni_x_rrad,mini,maxi,tolzbrent,psniData)
    psni_x_rrad = x

    if (present(bfoldstar)) then
       bfoldstar = - (psni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function psni_x_rrad

  function find_psni_x_rrad(x,psniData)   
    implicit none
    real(kp) :: find_psni_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: psniData

    real(kp) :: primStar,alpha,f,CalFplusprimEnd,potStar,epsOneStar

    alpha=psniData%real1
    f=psniData%real2
    CalFplusprimEnd = psniData%real3

    primStar = psni_efold_primitive(x,alpha,f)
    epsOneStar = psni_epsilon_one(x,alpha,f)
    potStar = psni_norm_potential(x,alpha)

    find_psni_x_rrad = find_reheat_rrad(primStar,calFplusprimEnd,epsOneStar,potStar)
  
  end function find_psni_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function psni_x_rreh(alpha,f,lnRreh,bfoldstar)    
    implicit none
    real(kp) :: psni_x_rreh
    real(kp), intent(in) :: alpha,f,lnRreh
    real(kp), intent(out), optional :: bfoldstar

    real(kp), parameter :: tolzbrent=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,xend,potEnd

    type(transfert) :: psniData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    xEnd = psni_x_endinf(alpha,f)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,f)
    potEnd = psni_norm_potential(xEnd,alpha)
    primEnd = psni_efold_primitive(xEnd,alpha,f)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    psniData%real1 = alpha 
    psniData%real2 = f
    psniData%real3 = calF + primEnd

    mini=epsilon(1._kp)
    maxi = xEnd


    x = zbrent(find_psni_x_rreh,mini,maxi,tolzbrent,psniData)
    psni_x_rreh = x

    if (present(bfoldstar)) then
       bfoldstar = - (psni_efold_primitive(x,alpha,f) - primEnd)
    endif

  end function psni_x_rreh

  function find_psni_x_rreh(x,psniData)   
    implicit none
    real(kp) :: find_psni_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: psniData

    real(kp) :: primStar,alpha,f,CalFplusprimEnd,potStar

    alpha=psniData%real1
    f=psniData%real2
    CalFplusprimEnd = psniData%real3

    primStar = psni_efold_primitive(x,alpha,f)    
    potStar = psni_norm_potential(x,alpha)

    find_psni_x_rreh = find_reheat_rreh(primStar,calFplusprimEnd,potStar)
  
  end function find_psni_x_rreh



  function psni_lnrhoreh_max(alpha,f,Pstar) 
    implicit none
    real(kp) :: psni_lnrhoreh_max
    real(kp), intent(in) :: alpha,f,Pstar

    real(kp) :: xEnd, potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp),parameter :: wrad=1._kp/3._kp
    real(kp),parameter :: junk=0_kp

    real(kp) :: lnRhoEnd
    
    xEnd = psni_x_endinf(alpha,f)
    potEnd  = psni_norm_potential(xEnd,alpha)
    epsOneEnd = psni_epsilon_one(xEnd,alpha,f)

!   Trick to return x such that rho_reh=rho_end

    x = psni_x_star(alpha,f,wrad,junk,Pstar)    
    potStar = psni_norm_potential(x,alpha)
    epsOneStar = psni_epsilon_one(x,alpha,f)


    
    if (.not.slowroll_validity(epsOneStar)) stop 'psni_lnrhoreh_max: slow-roll violated!'
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    psni_lnrhoreh_max = lnRhoEnd

  end function psni_lnrhoreh_max

  
end module psnireheat

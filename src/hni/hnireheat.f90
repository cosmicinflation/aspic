!hybrid natural reheating functions in the slow-roll approximation

module hnireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use hnisr, only : hni_epsilon_one, hni_epsilon_two, hni_epsilon_three
  use hnisr, only : hni_norm_potential,hni_efold_primitive,hni_x_endinf
  implicit none

  private

  public hni_x_star, hni_lnrhoreh_max
  public hni_x_rrad, hni_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function hni_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: hni_x_star
    real(kp), intent(in) :: alpha,phi0,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,phi0)
    potEnd = hni_norm_potential(xEnd,alpha,phi0)
    primEnd = hni_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = phi0
    hniData%real3 = w
    hniData%real4 = calF + primEnd

    mini = epsilon(1._kp) !potential maximum
    maxi = xend*(1._kp-epsilon(1._kp))

    x = zbrent(find_hni_x_star,mini,maxi,tolFind,hniData)
    hni_x_star = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function hni_x_star

  function find_hni_x_star(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,phi0,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=hniData%real1
    phi0=hniData%real2
    w = hniData%real3
    CalFplusPrimEnd = hniData%real4

    primStar = hni_efold_primitive(x,alpha,phi0)
    epsOneStar = hni_epsilon_one(x,alpha,phi0)
    potStar = hni_norm_potential(x,alpha,phi0)

    find_hni_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_hni_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function hni_x_rrad(alpha,phi0,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: hni_x_rrad
    real(kp), intent(in) :: alpha,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,phi0)
    potEnd = hni_norm_potential(xEnd,alpha,phi0)
    primEnd = hni_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = phi0
    hniData%real3 = calF + primEnd

    mini = epsilon(1._kp) !potential maximum
    maxi = xend*(1._kp-epsilon(1._kp))

    x = zbrent(find_hni_x_rrad,mini,maxi,tolFind,hniData)
    hni_x_rrad = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function hni_x_rrad

  function find_hni_x_rrad(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,phi0,CalFplusPrimEnd,potStar,epsOneStar

    alpha=hniData%real1
    phi0=hniData%real2
    CalFplusPrimEnd = hniData%real3

    primStar = hni_efold_primitive(x,alpha,phi0)
    epsOneStar = hni_epsilon_one(x,alpha,phi0)
    potStar = hni_norm_potential(x,alpha,phi0)

    find_hni_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_hni_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function hni_x_rreh(alpha,phi0,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: hni_x_rreh
    real(kp), intent(in) :: alpha,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: hniData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif
    
    epsOneEnd = hni_epsilon_one(xEnd,alpha,phi0)
    potEnd = hni_norm_potential(xEnd,alpha,phi0)
    primEnd = hni_efold_primitive(xEnd,alpha,phi0)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    hniData%real1 = alpha
    hniData%real2 = phi0
    hniData%real3 = calF + primEnd

    mini = epsilon(1._kp) !potential maximum
    maxi = xend*(1._kp-epsilon(1._kp))

    x = zbrent(find_hni_x_rreh,mini,maxi,tolFind,hniData)
    hni_x_rreh = x

    if (present(bfold)) then
       bfold = -(hni_efold_primitive(x,alpha,phi0) - primEnd)
    endif

  end function hni_x_rreh

  function find_hni_x_rreh(x,hniData)   
    implicit none
    real(kp) :: find_hni_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: hniData

    real(kp) :: primStar,alpha,phi0,CalFplusPrimEnd,potStar

    alpha=hniData%real1
    phi0=hniData%real2
    CalFplusPrimEnd = hniData%real3

    primStar = hni_efold_primitive(x,alpha,phi0)
    potStar = hni_norm_potential(x,alpha,phi0)

    find_hni_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_hni_x_rreh


  function hni_lnrhoreh_max(alpha,phi0,xend,Pstar) 
    implicit none
    real(kp) :: hni_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = hni_norm_potential(xEnd,alpha,phi0)
    epsOneEnd = hni_epsilon_one(xEnd,alpha,phi0)
       
    x = hni_x_star(alpha,phi0,xend,wrad,junk,Pstar)

    potStar = hni_norm_potential(x,alpha,phi0)
    epsOneStar = hni_epsilon_one(x,alpha,phi0)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'hni_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    hni_lnrhoreh_max = lnRhoEnd

  end function hni_lnrhoreh_max

  
  
end module hnireheat

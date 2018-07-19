!oip inflation reheaoing functions in the slow-roll approximaoions

module oireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf,ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use oisr, only : oi_epsilon_one, oi_epsilon_two, oi_epsilon_three
  use oisr, only : oi_norm_potential,oi_efold_primitive,oi_x_endinf
  implicit none

  private

  public oi_x_star, oi_lnrhoreh_max
  public oi_x_rrad, oi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function oi_x_star(alpha,phi0,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: oi_x_star
    real(kp), intent(in) :: alpha,phi0,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oiData
    

    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = oi_epsilon_one(xEnd,alpha,phi0)
    potEnd = oi_norm_potential(xEnd,alpha,phi0)
    primEnd = oi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    oiData%real1 = alpha
    oiData%real2 = phi0
    oiData%real3 = w
    oiData%real4 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    x = zbrent(find_oi_x_star,mini,maxi,tolFind,oiData)
    oi_x_star = x

    if (present(bfold)) then
       bfold = -(oi_efold_primitive(x,alpha,phi0) - primEnd)
    endif


  end function oi_x_star

  function find_oi_x_star(x,oiData)   
    implicit none
    real(kp) :: find_oi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oiData

    real(kp) :: primStar,alpha,phi0,w,CalFplusPrimEnd,potStar,epsOneStar

    alpha=oiData%real1
    phi0=oiData%real2
    w = oiData%real3
    CalFplusPrimEnd = oiData%real4

    primStar = oi_efold_primitive(x,alpha,phi0)
    epsOneStar = oi_epsilon_one(x,alpha,phi0)
    potStar = oi_norm_potential(x,alpha,phi0)

    find_oi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_oi_x_star


!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function oi_x_rrad(alpha,phi0,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: oi_x_rrad
    real(kp), intent(in) :: alpha,phi0,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oiData
    

    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = oi_epsilon_one(xEnd,alpha,phi0)
    potEnd = oi_norm_potential(xEnd,alpha,phi0)
    primEnd = oi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    oiData%real1 = alpha
    oiData%real2 = phi0
    oiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    x = zbrent(find_oi_x_rrad,mini,maxi,tolFind,oiData)
    oi_x_rrad = x

    if (present(bfold)) then
       bfold = -(oi_efold_primitive(x,alpha,phi0) - primEnd)
    endif


  end function oi_x_rrad

  function find_oi_x_rrad(x,oiData)   
    implicit none
    real(kp) :: find_oi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oiData

    real(kp) :: primStar,alpha,phi0,CalFplusPrimEnd,potStar,epsOneStar

    alpha=oiData%real1
    phi0=oiData%real2
    CalFplusPrimEnd = oiData%real3

    primStar = oi_efold_primitive(x,alpha,phi0)
    epsOneStar = oi_epsilon_one(x,alpha,phi0)
    potStar = oi_norm_potential(x,alpha,phi0)

    find_oi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_oi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function oi_x_rreh(alpha,phi0,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: oi_x_rreh
    real(kp), intent(in) :: alpha,phi0,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd

    type(transfert) :: oiData
    

    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = oi_epsilon_one(xEnd,alpha,phi0)
    potEnd = oi_norm_potential(xEnd,alpha,phi0)
    primEnd = oi_efold_primitive(xEnd,alpha,phi0)
   
    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    oiData%real1 = alpha
    oiData%real2 = phi0
    oiData%real3 = calF + primEnd

    mini = xEnd*(1._kp+epsilon(1._kp))
    maxi = mini/epsilon(1._kp)

    x = zbrent(find_oi_x_rreh,mini,maxi,tolFind,oiData)
    oi_x_rreh = x

    if (present(bfold)) then
       bfold = -(oi_efold_primitive(x,alpha,phi0) - primEnd)
    endif


  end function oi_x_rreh

  function find_oi_x_rreh(x,oiData)   
    implicit none
    real(kp) :: find_oi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: oiData

    real(kp) :: primStar,alpha,phi0,CalFplusPrimEnd,potStar

    alpha=oiData%real1
    phi0=oiData%real2
    CalFplusPrimEnd = oiData%real3

    primStar = oi_efold_primitive(x,alpha,phi0)
    potStar = oi_norm_potential(x,alpha,phi0)

    find_oi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_oi_x_rreh



  function oi_lnrhoreh_max(alpha,phi0,xend,Pstar) 
    implicit none
    real(kp) :: oi_lnrhoreh_max
    real(kp), intent(in) :: alpha,phi0,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
        
    potEnd  = oi_norm_potential(xEnd,alpha,phi0)
    epsOneEnd = oi_epsilon_one(xEnd,alpha,phi0)

    x = oi_x_star(alpha,phi0,xend,wrad,junk,Pstar)    
    potStar = oi_norm_potential(x,alpha,phi0)
    epsOneStar = oi_epsilon_one(x,alpha,phi0)


    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'oi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    oi_lnrhoreh_max = lnRhoEnd

  end function oi_lnrhoreh_max

  
  
end module oireheat

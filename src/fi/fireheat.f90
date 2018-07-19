!hybrid natural reheating functions in the slow-roll approximation

module fireheat
  use infprec, only : kp, tolkp, transfert
  use inftools, only : zbrent
  use srreheat, only : get_calfconst, find_reheat, slowroll_validity
  use srreheat, only : display, pi, Nzero, ln_rho_endinf, ln_rho_reheat
  use srreheat, only : find_reheat_rrad, find_reheat_rreh
  use srreheat, only : get_calfconst_rrad, get_calfconst_rreh
  use fisr, only : fi_epsilon_one, fi_epsilon_two, fi_epsilon_three
  use fisr, only : fi_norm_potential, fi_efold_primitive, fi_x_endinf
  use fisr, only : fi_x_epsoneunity
  implicit none

  private

  public fi_x_star, fi_lnrhoreh_max
  public fi_x_rrad, fi_x_rreh

contains

!returns x given potential parameters, scalar power, wreh and
!lnrhoreh. If present, returns the correspoding bfoldstar
  function fi_x_star(delta,n,xend,w,lnRhoReh,Pstar,bfold)    
    implicit none
    real(kp) :: fi_x_star
    real(kp), intent(in) :: delta,n,xend,w,lnRhoReh,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsone
    
    type(transfert) :: fiData
    
    if (w.eq.1._kp/3._kp) then
       if (display) write(*,*)'w = 1/3 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = fi_epsilon_one(xEnd,delta,n)
    potEnd = fi_norm_potential(xEnd,delta,n)
    primEnd = fi_efold_primitive(xEnd,delta,n)

    calF = get_calfconst(lnRhoReh,Pstar,w,epsOneEnd,potEnd)

    fiData%real1 = delta
    fiData%real2 = n
    fiData%real3 = w
    fiData%real4 = calF + primEnd

    xepsone = fi_x_epsoneunity(delta,n)
    
    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xEpsOne(2)*(1._kp-epsilon(1._kp))

    x = zbrent(find_fi_x_star,mini,maxi,tolFind,fiData)
    fi_x_star = x

    if (present(bfold)) then
       bfold = -(fi_efold_primitive(x,delta,n) - primEnd)
    endif

  end function fi_x_star

  function find_fi_x_star(x,fiData)   
    implicit none
    real(kp) :: find_fi_x_star
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData

    real(kp) :: primStar,delta,n,w,CalFplusPrimEnd,potStar,epsOneStar

    delta=fiData%real1
    n=fiData%real2
    w = fiData%real3
    CalFplusPrimEnd = fiData%real4

    primStar = fi_efold_primitive(x,delta,n)
    epsOneStar = fi_epsilon_one(x,delta,n)
    potStar = fi_norm_potential(x,delta,n)

    find_fi_x_star = find_reheat(PrimStar,calFplusPrimEnd,w,epsOneStar,potStar)

  end function find_fi_x_star

!returns x given potential parameters, scalar power, and lnRrad.
!If present, returns the corresponding bfoldstar
  function fi_x_rrad(delta,n,xend,lnRrad,Pstar,bfold)    
    implicit none
    real(kp) :: fi_x_rrad
    real(kp), intent(in) :: delta,n,xend,lnRrad,Pstar
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsone
    
    type(transfert) :: fiData
    
    if (lnRrad.eq.0._kp) then
       if (display) write(*,*)'Rrad=1 : solving for rhoReh = rhoEnd'
    endif
   
    epsOneEnd = fi_epsilon_one(xEnd,delta,n)
    potEnd = fi_norm_potential(xEnd,delta,n)
    primEnd = fi_efold_primitive(xEnd,delta,n)

    calF = get_calfconst_rrad(lnRrad,Pstar,epsOneEnd,potEnd)

    fiData%real1 = delta
    fiData%real2 = n
    fiData%real3 = calF + primEnd

    xepsone = fi_x_epsoneunity(delta,n)
    
    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xEpsOne(2)*(1._kp-epsilon(1._kp))

    x = zbrent(find_fi_x_rrad,mini,maxi,tolFind,fiData)
    fi_x_rrad = x

    if (present(bfold)) then
       bfold = -(fi_efold_primitive(x,delta,n) - primEnd)
    endif

  end function fi_x_rrad

  function find_fi_x_rrad(x,fiData)   
    implicit none
    real(kp) :: find_fi_x_rrad
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData

    real(kp) :: primStar,delta,n,CalFplusPrimEnd,potStar,epsOneStar

    delta=fiData%real1
    n=fiData%real2
    CalFplusPrimEnd = fiData%real3

    primStar = fi_efold_primitive(x,delta,n)
    epsOneStar = fi_epsilon_one(x,delta,n)
    potStar = fi_norm_potential(x,delta,n)

    find_fi_x_rrad = find_reheat_rrad(PrimStar,calFplusPrimEnd,epsOneStar,potStar)

  end function find_fi_x_rrad


!returns x given potential parameters, scalar power, and lnRreh.
!If present, returns the corresponding bfoldstar
  function fi_x_rreh(delta,n,xend,lnRreh,bfold)    
    implicit none
    real(kp) :: fi_x_rreh
    real(kp), intent(in) :: delta,n,xend,lnRreh
    real(kp), intent(out), optional :: bfold

    real(kp), parameter :: tolFind=tolkp
    real(kp) :: mini,maxi,calF,x
    real(kp) :: primEnd,epsOneEnd,potEnd
    real(kp), dimension(2) :: xepsone
    
    type(transfert) :: fiData
    
    if (lnRreh.eq.0._kp) then
       if (display) write(*,*)'Rreh=1 : solving for rhoReh = rhoEnd'
    endif

    epsOneEnd = fi_epsilon_one(xEnd,delta,n)
    potEnd = fi_norm_potential(xEnd,delta,n)
    primEnd = fi_efold_primitive(xEnd,delta,n)

    calF = get_calfconst_rreh(lnRreh,epsOneEnd,potEnd)

    fiData%real1 = delta
    fiData%real2 = n
    fiData%real3 = calF + primEnd

    xepsone = fi_x_epsoneunity(delta,n)
        
    mini = xend*(1._kp+epsilon(1._kp))
    maxi = xEpsOne(2)*(1._kp-epsilon(1._kp))

    x = zbrent(find_fi_x_rreh,mini,maxi,tolFind,fiData)
    fi_x_rreh = x

    if (present(bfold)) then
       bfold = -(fi_efold_primitive(x,delta,n) - primEnd)
    endif

  end function fi_x_rreh

  function find_fi_x_rreh(x,fiData)   
    implicit none
    real(kp) :: find_fi_x_rreh
    real(kp), intent(in) :: x
    type(transfert), optional, intent(inout) :: fiData

    real(kp) :: primStar,delta,n,CalFplusPrimEnd,potStar

    delta=fiData%real1
    n=fiData%real2
    CalFplusPrimEnd = fiData%real3

    primStar = fi_efold_primitive(x,delta,n)
    potStar = fi_norm_potential(x,delta,n)

    find_fi_x_rreh = find_reheat_rreh(PrimStar,calFplusPrimEnd,potStar)

  end function find_fi_x_rreh


  function fi_lnrhoreh_max(delta,n,xend,Pstar) 
    implicit none
    real(kp) :: fi_lnrhoreh_max
    real(kp), intent(in) :: delta,n,xend,Pstar

    real(kp) :: potEnd, epsOneEnd
    real(kp) :: x, potStar, epsOneStar

    real(kp), parameter :: wrad = 1._kp/3._kp
    real(kp), parameter :: junk= 0._kp
    real(kp) :: lnRhoEnd
    
    potEnd  = fi_norm_potential(xEnd,delta,n)
    epsOneEnd = fi_epsilon_one(xEnd,delta,n)
       
    x = fi_x_star(delta,n,xend,wrad,junk,Pstar)

    potStar = fi_norm_potential(x,delta,n)
    epsOneStar = fi_epsilon_one(x,delta,n)
    
    if (.not.slowroll_validity(epsOneStar)) then
        print*,'xstar=',x,'  epsOneStar=',epsOneStar 
        stop 'fi_lnrhoreh_max: slow-roll violated!'
    endif
    
    lnRhoEnd = ln_rho_endinf(Pstar,epsOneStar,epsOneEnd,potEnd/potStar)

    fi_lnrhoreh_max = lnRhoEnd

  end function fi_lnrhoreh_max

  
  
end module fireheat
